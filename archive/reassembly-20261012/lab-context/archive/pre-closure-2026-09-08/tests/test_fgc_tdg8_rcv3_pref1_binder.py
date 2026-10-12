from __future__ import annotations

import ast
from copy import deepcopy
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock
import warnings
import zipfile

import numpy as np

from recursive_horizons.fgc.evolution import tdg8_rcv3_pref1_binder as binder


ROOT = Path(__file__).resolve().parents[1]


class TDG8RCV3PREF1BinderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = (ROOT / binder.CONFIG_PATH).read_bytes()
        cls.result = (ROOT / binder.RESULT_PATH).read_bytes()

    def test_compact_result_and_exact_contract_pass(self) -> None:
        result = binder.validate_compact_result(self.config, self.result)
        self.assertEqual(result["artifact_id"], binder.ARTIFACT_ID)
        claims = result["artifact_payload"]["claims"]
        self.assertIs(claims["installed_projection_authenticated"], True)
        self.assertIs(
            result["artifact_payload"]["scope"]["bounded_event_execution_authorized"],
            False,
        )
        self.assertIs(claims["common_event_completed"], False)
        self.assertIs(claims["candidate_execution_authorized"], False)
        self.assertIs(claims["physical_result_earned"], False)

    def test_live_projection_passes_independent_full_binding(self) -> None:
        observed = binder.bind_installed_projection(self.config, ROOT)
        self.assertEqual(observed, binder.expected_anchor(self.config))
        self.assertEqual(observed["receipt"]["receipt_sha256"], binder.RECEIPT_SHA256)
        self.assertEqual(observed["receipt"]["raw_sha256"], binder.RECEIPT_RAW_SHA256)
        self.assertEqual(observed["members"][0]["retry_count"], 2)
        self.assertEqual(observed["members"][0]["pending_cap_hex"], "0x1.aaa9612df8000p-11")

    def test_binder_has_no_bootstrap_writer_or_campaign_runtime_import(self) -> None:
        source = (ROOT / "src/recursive_horizons/fgc/evolution/tdg8_rcv3_pref1_binder.py").read_text()
        tree = ast.parse(source)
        imports: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imports.append(node.module or "")
        self.assertFalse(
            any("tdg8_rcv3_fork_runtime" in module for module in imports),
            imports,
        )
        self.assertNotIn("HLT16CampaignStore(", source)
        self.assertNotIn("advance_member(", source)
        self.assertNotIn("solve_accelerations(", source)

    def test_duplicate_json_key_is_rejected(self) -> None:
        with self.assertRaises(binder.TDG8RCV3PREF1Error) as caught:
            binder._json(b'{"a":1,"a":1}', "duplicate")
        self.assertEqual(caught.exception.stop_id, "PREF1_JSON_DRIFT")

    def test_claim_promotion_and_config_mutation_are_rejected(self) -> None:
        promoted = self.config.replace(
            b"candidate_execution_authorized = false",
            b"candidate_execution_authorized = true",
        )
        with self.assertRaises(binder.TDG8RCV3PREF1Error) as caught:
            binder.validate_compact_result(promoted, self.result)
        self.assertEqual(caught.exception.stop_id, "PREF1_CONFIG_DRIFT")

        value = json.loads(self.result)
        value["artifact_payload"]["claims"]["physical_result_earned"] = True
        with self.assertRaises(binder.TDG8RCV3PREF1Error) as caught:
            binder.validate_compact_result(self.config, binder.canonical_result(value))
        self.assertEqual(caught.exception.stop_id, "PREF1_COMPACT_DRIFT")

    def _regular_leaf(self, origin: Path, relative: str, pattern=None) -> Path:
        name = Path(relative).name
        if "/" not in relative or (pattern is not None and pattern.fullmatch(name) is None):
            raise AssertionError(f"generation-9 leaf is malformed: {relative}")
        path = origin / relative
        if path.parent.is_symlink() or path.is_symlink() or not path.is_file():
            raise AssertionError(f"generation-9 leaf is missing or a symlink: {relative}")
        return path

    def _prefix_indexed(self, origin: Path, directory: str, pattern, limit: int):
        selected = []
        numbers = []
        root = origin / directory
        if root.is_symlink() or not root.is_dir():
            raise AssertionError(f"generation-9 directory is missing or a symlink: {directory}")
        for entry in sorted(root.iterdir(), key=lambda item: item.name):
            match = pattern.fullmatch(entry.name)
            if match is None or int(match.group("number")) > limit:
                continue
            relative = f"{directory}/{entry.name}"
            selected.append((relative, self._regular_leaf(origin, relative, pattern)))
            numbers.append(int(match.group("number")))
        if numbers != list(range(limit + 1)):
            raise AssertionError(f"{directory} 0..{limit} is incomplete")
        return selected

    def _generation9_prefix_files(self, origin: Path) -> list[tuple[str, Path]]:
        if origin.is_symlink() or not origin.is_dir():
            raise AssertionError("generation-9 origin is missing or a symlink")
        files = [
            *self._prefix_indexed(
                origin, "checkpoints", binder._CHECKPOINT_NAME, binder.ANCHOR_GENERATION,
            ),
            *self._prefix_indexed(
                origin, "journal", binder._JOURNAL_NAME, binder.ANCHOR_JOURNAL_SEQUENCE,
            ),
        ]
        states: set[str] = set()
        for relative, path in files:
            if not relative.startswith("checkpoints/"):
                continue
            checkpoint = binder._json(path.read_bytes(), relative)
            states.update(
                binder._sha(digest, relative) for digest in checkpoint.get("states", {}).values()
            )
            states.update(
                binder._sha(member["descriptor_sha256"], relative)
                for member in checkpoint.get("members", {}).values()
            )
        payloads: set[str] = set()
        for digest in sorted(states):
            relative = f"states/{digest}.json"
            path = self._regular_leaf(origin, relative, binder._STATE_NAME)
            files.append((relative, path))
            semantic = binder._json(path.read_bytes(), relative).get("semantic_sha256")
            if semantic is not None:
                payloads.add(binder._sha(semantic, relative))
        for digest in sorted(payloads):
            relative = f"payloads/{digest}.npz"
            files.append((relative, self._regular_leaf(origin, relative, binder._PAYLOAD_NAME)))
        receipt = "receipts/generation-zero.json"
        files.append((receipt, self._regular_leaf(origin, receipt)))
        return sorted(files, key=lambda item: item[0])

    def _generation9_projection_leaves(self, origin: Path | None = None) -> dict[str, binder._Leaf]:
        leaves: dict[str, binder._Leaf] = {}
        root = ROOT / binder.DESTINATION_STORE if origin is None else origin
        for relative, entry in self._generation9_prefix_files(root):
            raw = entry.read_bytes()
            leaves[relative] = binder._Leaf(
                relative, len(raw), sha256(raw).hexdigest(), raw
            )
        return dict(sorted(leaves.items()))

    def _copy_generation9_prefix(self, origin: Path, destination: Path) -> None:
        for directory in binder.STORE_DIRECTORIES:
            (destination / directory).mkdir(parents=True, exist_ok=True)
        for relative, entry in self._generation9_prefix_files(origin):
            shutil.copy2(entry, destination / relative)

    def _copy_repository_evidence(self, temporary: Path, origin: Path | None = None) -> None:
        source = temporary / binder.SOURCE_STORE
        wrapper = temporary / binder.DESTINATION_WRAPPER
        destination = temporary / binder.DESTINATION_STORE
        source.parent.mkdir(parents=True)
        wrapper.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(ROOT / binder.SOURCE_STORE, source)
        self._copy_generation9_prefix(origin or ROOT / binder.DESTINATION_STORE, destination)
        shutil.copy2(ROOT / binder.RECEIPT_PATH, wrapper / "bootstrap-receipt.json")

    def _bind_without_git(self, temporary: Path) -> None:
        with mock.patch.object(binder, "_bind_authority", return_value=None):
            binder.bind_installed_projection(self.config, temporary)

    def test_destination_byte_mutation_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="fgc-rcv3-pref1-") as directory:
            root = Path(directory)
            self._copy_repository_evidence(root)
            payload = next((root / binder.DESTINATION_STORE / "payloads").glob("*.npz"))
            raw = payload.read_bytes()
            payload.write_bytes(raw[:-1] + bytes([raw[-1] ^ 1]))
            with self.assertRaises(binder.TDG8RCV3PREF1Error) as caught:
                self._bind_without_git(root)
            self.assertIn(caught.exception.stop_id, {"PREF1_TREE_DRIFT", "PREF1_PROJECTION_DRIFT"})

    def test_source_mutation_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="fgc-rcv3-pref1-") as directory:
            root = Path(directory)
            self._copy_repository_evidence(root)
            receipt = root / binder.SOURCE_STORE / "receipts/generation-zero.json"
            raw = receipt.read_bytes()
            receipt.write_bytes(raw[:-1] + b" ")
            with self.assertRaises(binder.TDG8RCV3PREF1Error) as caught:
                self._bind_without_git(root)
            self.assertEqual(caught.exception.stop_id, "PREF1_TREE_DRIFT")

    def test_terminal_suffix_or_lock_in_projection_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="fgc-rcv3-pref1-") as directory:
            root = Path(directory)
            self._copy_repository_evidence(root)
            source = root / binder.SOURCE_STORE / binder.EXCLUDED_PATHS[0]
            destination = root / binder.DESTINATION_STORE / binder.EXCLUDED_PATHS[0]
            shutil.copy2(source, destination)
            with self.assertRaises(binder.TDG8RCV3PREF1Error) as caught:
                self._bind_without_git(root)
            self.assertEqual(caught.exception.stop_id, "PREF1_TREE_DRIFT")

    def test_nested_symlink_is_rejected_without_following(self) -> None:
        with tempfile.TemporaryDirectory(prefix="fgc-rcv3-pref1-") as directory:
            root = Path(directory)
            self._copy_repository_evidence(root)
            payload = next((root / binder.DESTINATION_STORE / "payloads").glob("*.npz"))
            outside = root / "outside"
            outside.write_bytes(payload.read_bytes())
            payload.unlink()
            payload.symlink_to(outside)
            with self.assertRaises(binder.TDG8RCV3PREF1Error) as caught:
                self._bind_without_git(root)
            self.assertEqual(caught.exception.stop_id, "PREF1_PATH_UNSAFE")

    def test_symlinked_store_root_is_rejected_without_following(self) -> None:
        with tempfile.TemporaryDirectory(prefix="fgc-rcv3-pref1-") as directory:
            root = Path(directory)
            self._copy_repository_evidence(root)
            store = root / binder.DESTINATION_STORE
            outside = root / "outside-calibration"
            store.rename(outside)
            store.symlink_to(outside, target_is_directory=True)
            with self.assertRaises(binder.TDG8RCV3PREF1Error) as caught:
                self._bind_without_git(root)
            self.assertEqual(caught.exception.stop_id, "PREF1_PATH_UNSAFE")

    def test_symlinked_inner_directory_is_rejected_without_following(self) -> None:
        with tempfile.TemporaryDirectory(prefix="fgc-rcv3-pref1-") as directory:
            root = Path(directory)
            self._copy_repository_evidence(root)
            payloads = root / binder.DESTINATION_STORE / "payloads"
            outside = root / "outside-payloads"
            payloads.rename(outside)
            payloads.symlink_to(outside, target_is_directory=True)
            with self.assertRaises(binder.TDG8RCV3PREF1Error) as caught:
                self._bind_without_git(root)
            self.assertEqual(caught.exception.stop_id, "PREF1_PATH_UNSAFE")

    def test_receipt_authority_mutation_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="fgc-rcv3-pref1-") as directory:
            root = Path(directory)
            self._copy_repository_evidence(root)
            path = root / binder.RECEIPT_PATH
            value = json.loads(path.read_bytes())
            value["external_install_authority"]["authority_commit_sha"] = "0" * 40
            body = dict(value)
            body.pop("receipt_sha256")
            value["receipt_sha256"] = binder._digest(body)
            path.write_bytes(binder._canonical(value))
            with self.assertRaises(binder.TDG8RCV3PREF1Error) as caught:
                self._bind_without_git(root)
            self.assertIn(caught.exception.stop_id, {"PREF1_RECEIPT_DRIFT", "PREF1_PATH_UNSAFE"})

    def test_depth_two_retry_mutation_is_rejected(self) -> None:
        leaves = self._generation9_projection_leaves()
        checkpoints, _ = binder._validate_lineage(leaves)
        descriptors = binder._descriptor_objects(leaves, checkpoints)
        payloads = binder._load_payloads(leaves, descriptors)
        altered = deepcopy(checkpoints[binder.ANCHOR_GENERATION])
        retry = altered["members"]["RK4-2049"]["cursor"]["retry_successor_payload_or_none"]
        retry["retry_count"] = 1
        with self.assertRaises(binder.TDG8RCV3PREF1Error) as caught:
            binder._member_evidence(altered, descriptors, payloads)
        self.assertIn(caught.exception.stop_id, {"PREF1_MEMBER_DRIFT", "PREF1_RETRY_DRIFT"})

    def test_duplicate_zip_member_and_object_array_are_rejected(self) -> None:
        leaves = self._generation9_projection_leaves()
        checkpoints, _ = binder._validate_lineage(leaves)
        descriptors = binder._descriptor_objects(leaves, checkpoints)
        descriptor = next(iter(descriptors.values()))

        duplicate = BytesIO()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            with zipfile.ZipFile(duplicate, "w") as archive:
                archive.writestr("u.npy", b"x")
                archive.writestr("u.npy", b"x")
        with self.assertRaises(binder.TDG8RCV3PREF1Error) as caught:
            binder._decode_npz(duplicate.getvalue(), descriptor, "duplicate")
        self.assertEqual(caught.exception.stop_id, "PREF1_PAYLOAD_DRIFT")

        malicious = BytesIO()
        with zipfile.ZipFile(malicious, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for name in binder.ARRAY_NAMES:
                encoded = BytesIO()
                value = np.asarray([object()], dtype=object) if name == "u" else np.asarray([0.0], dtype="<f8")
                np.save(encoded, value, allow_pickle=True)
                archive.writestr(f"{name}.npy", encoded.getvalue())
        with self.assertRaises(binder.TDG8RCV3PREF1Error) as caught:
            binder._decode_npz(malicious.getvalue(), descriptor, "object")
        self.assertEqual(caught.exception.stop_id, "PREF1_PAYLOAD_DRIFT")

    def test_unreferenced_successor_blobs_are_excluded_from_generation9_prefix(self) -> None:
        extra_state = f"{'f' * 64}.json"
        extra_payload = f"{'e' * 64}.npz"
        extra_receipt = "future-successor.json"
        extras = {
            f"states/{extra_state}",
            f"payloads/{extra_payload}",
            f"receipts/{extra_receipt}",
        }
        with tempfile.TemporaryDirectory(prefix="fgc-rcv3-pref1-") as directory:
            root = Path(directory)
            origin = root / "contaminated-origin"
            self._copy_generation9_prefix(ROOT / binder.DESTINATION_STORE, origin)
            (origin / "states" / extra_state).write_bytes(b"{}")
            (origin / "payloads" / extra_payload).write_bytes(b"npz")
            (origin / "receipts" / extra_receipt).write_bytes(b"{}")
            leaves = self._generation9_projection_leaves(origin)
            snapshot = binder._Snapshot(leaves, binder.STORE_DIRECTORIES)
            self.assertEqual(list(leaves), sorted(leaves))
            self.assertEqual(len(leaves), binder.PROJECTED_LEAF_COUNT)
            self.assertEqual(snapshot.byte_count(), binder.PROJECTED_BYTE_COUNT)
            self.assertEqual(snapshot.manifest_sha256(), binder.PROJECTED_TREE_SHA256)
            self.assertTrue(extras.isdisjoint(leaves))
            self._copy_repository_evidence(root, origin)
            destination = root / binder.DESTINATION_STORE
            self.assertFalse((destination / "states" / extra_state).exists())
            self.assertFalse((destination / "payloads" / extra_payload).exists())
            self.assertFalse((destination / "receipts" / extra_receipt).exists())
            with mock.patch.object(binder, "_bind_authority", return_value=None):
                observed = binder.bind_installed_projection(self.config, root)
            self.assertEqual(observed, binder.expected_anchor(self.config))


if __name__ == "__main__":
    unittest.main()
