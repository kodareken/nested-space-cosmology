from __future__ import annotations

from copy import deepcopy
import ast
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
from types import SimpleNamespace
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]

from scripts import bootstrap_fgc_tdg8_rcv3 as bootstrap
from recursive_horizons.fgc.evolution import tdg8_rcv3_freeze as freeze
from recursive_horizons.fgc.evolution import tdg8_rcv3_fork_runtime as runtime


class TDG8RCV3FreezeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = freeze.read_nofollow(
            ROOT, freeze.CONFIG_PATH, "test RCV3 FRZ1 config"
        )
        cls.result = freeze.read_nofollow(
            ROOT, freeze.RESULT_PATH, "test RCV3 FRZ1 result"
        )
        cls.parsed = freeze.validate_compact(cls.config, cls.result)

    def mutated_result(self) -> dict[str, object]:
        return deepcopy(self.parsed)

    def assertCompactRejected(self, value: object) -> None:
        with self.assertRaises(freeze.TDG8RCV3FreezeError):
            freeze.validate_compact(self.config, freeze.canonical_result(value))

    def test_tracked_compact_result_passes_without_live_store(self) -> None:
        self.assertEqual(self.parsed["gate_status"], "pass")
        payload = self.parsed["artifact_payload"]
        self.assertTrue(payload["claims"]["bootstrap_authorized"])
        self.assertFalse(payload["claims"]["bootstrap_completed"])
        self.assertFalse(payload["claims"]["bounded_execution_authorized"])

    def test_live_terminal_and_prefix_match_the_frozen_result(self) -> None:
        observed = freeze.derive_source_evidence(ROOT)
        frozen = self.parsed["artifact_payload"]["source"]
        self.assertEqual(observed, frozen)
        self.assertEqual(observed["full_manifest"]["leaf_count"], 59)
        self.assertEqual(observed["prefix_manifest"]["leaf_count"], 50)

    def test_prefix_is_exactly_full_tree_minus_nine_frozen_leaves(self) -> None:
        source = self.parsed["artifact_payload"]["source"]
        full = {row["path"] for row in source["full_manifest"]["leaves"]}
        prefix = {row["path"] for row in source["prefix_manifest"]["leaves"]}
        self.assertEqual(full - prefix, set(freeze.EXCLUDED_PATHS))
        self.assertEqual(len(full - prefix), 9)
        self.assertFalse(any(path.startswith("locks/") for path in prefix))

    def test_coordination_and_retired_locks_are_bound_but_not_projected(self) -> None:
        source = self.parsed["artifact_payload"]["source"]
        full = {row["path"] for row in source["full_manifest"]["leaves"]}
        excluded_locks = {path for path in freeze.EXCLUDED_PATHS if path.startswith("locks/")}
        self.assertEqual(len(excluded_locks), 7)
        self.assertTrue(excluded_locks.issubset(full))
        self.assertIn("locks/bootstrap.guard", excluded_locks)
        self.assertIn("locks/writer.guard", excluded_locks)
        self.assertIn("locks/terminal.lock", excluded_locks)

    def test_generation_nine_retry_boundary_is_exact(self) -> None:
        anchor = self.parsed["artifact_payload"]["source"]["anchor"]
        self.assertEqual(anchor["generation9_checkpoint_sha256"], freeze.GENERATION9_CHECKPOINT_SHA256)
        self.assertEqual(anchor["sequence10_journal_sha256"], freeze.SEQUENCE10_JOURNAL_SHA256)
        self.assertEqual(anchor["member_mode"], "RETRY_PENDING")
        self.assertEqual(anchor["retry_depth"], 2)
        self.assertEqual(anchor["pending_cap_binary64_hex"], "0x1.aaa9612df8000p-11")
        self.assertEqual(anchor["generation9_disposition"], "nonterminal")

    def test_repair_sources_are_still_exact_commit_bytes(self) -> None:
        observed = freeze.derive_repair_evidence(ROOT)
        self.assertEqual(observed["commit"], freeze.REPAIR_COMMIT)
        self.assertEqual(observed["head_at_observation"], freeze.REPAIR_COMMIT)
        self.assertTrue(observed["repair_commit_is_ancestor"])
        self.assertTrue(observed["live_sources_equal_commit"])

    def test_observation_head_mutation_is_rejected(self) -> None:
        value = self.mutated_result()
        value["artifact_payload"]["repair"]["head_at_observation"] = "f" * 40
        self.assertCompactRejected(value)

    def test_bootstrap_implementation_hash_mutations_are_rejected(self) -> None:
        for field in ("runtime_sha256", "script_sha256"):
            with self.subTest(field=field):
                value = self.mutated_result()
                value["artifact_payload"]["bootstrap"][field] = "0" * 64
                self.assertCompactRejected(value)
        contract = freeze.expected_bootstrap_contract()
        contract["runtime_sha256"] = "0" * 64
        raw_by_path = {
            path: bootstrap._read_nofollow(ROOT, path)
            for path in (
                bootstrap.CONFIG_PATH,
                bootstrap.RESULT_PATH,
                bootstrap.FREEZE_MODULE_PATH,
                bootstrap.RUNTIME_PATH,
                bootstrap.SCRIPT_PATH,
            )
        }
        raw_by_path[bootstrap.CONFIG_PATH] = raw_by_path[bootstrap.CONFIG_PATH].replace(
            freeze.BOOTSTRAP_RUNTIME_SHA256.encode("ascii"), b"0" * 64
        )
        with self.assertRaises(bootstrap.TDG8RCV3BootstrapError):
            bootstrap._prevalidate_authority_bytes(raw_by_path)

    def test_one_bit_manifest_hash_mutation_is_rejected(self) -> None:
        value = self.mutated_result()
        rows = value["artifact_payload"]["source"]["prefix_manifest"]["leaves"]
        rows[0]["sha256"] = ("0" if rows[0]["sha256"][0] != "0" else "1") + rows[0]["sha256"][1:]
        self.assertCompactRejected(value)

    def test_leaf_path_and_count_mutations_are_rejected(self) -> None:
        value = self.mutated_result()
        value["artifact_payload"]["source"]["full_manifest"]["leaves"][0]["path"] = "states/unsafe.json"
        self.assertCompactRejected(value)
        value = self.mutated_result()
        value["artifact_payload"]["source"]["prefix_manifest"]["leaf_count"] = 51
        self.assertCompactRejected(value)

    def test_claim_promotion_is_rejected(self) -> None:
        value = self.mutated_result()
        value["artifact_payload"]["claims"]["bounded_execution_authorized"] = True
        self.assertCompactRejected(value)
        value = self.mutated_result()
        value["artifact_payload"]["claims"]["candidate_execution_authorized"] = True
        self.assertCompactRejected(value)

    def test_noncanonical_or_mutated_config_is_rejected(self) -> None:
        parsed = json.loads(self.result)
        compact = json.dumps(parsed, sort_keys=True, separators=(",", ":")).encode("ascii")
        with self.assertRaises(freeze.TDG8RCV3FreezeError):
            freeze.validate_compact(self.config, compact)
        mutated = self.config.replace(b"full_manifest_leaf_count = 59", b"full_manifest_leaf_count = 58")
        with self.assertRaises(freeze.TDG8RCV3FreezeError):
            freeze.parse_config(mutated)

    def _empty_store(self, root: Path) -> Path:
        store = root / "store"
        store.mkdir()
        for name in freeze.EXPECTED_DIRECTORIES:
            (store / name).mkdir()
        return store

    def test_tree_walk_rejects_symlink_leaf_and_directory(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            store = self._empty_store(root)
            target = root / "target"
            target.write_bytes(b"x")
            os.symlink(target, store / "states" / "linked.json")
            with self.assertRaises(freeze.TDG8RCV3FreezeError):
                freeze._walk_tree(store)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            store = root / "store"
            store.mkdir()
            target = root / "target-directory"
            target.mkdir()
            for name in freeze.EXPECTED_DIRECTORIES:
                if name == "states":
                    os.symlink(target, store / name)
                else:
                    (store / name).mkdir()
            with self.assertRaises(freeze.TDG8RCV3FreezeError):
                freeze._walk_tree(store)

    def test_tree_walk_rejects_hardlinked_and_foreign_leaves(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            store = self._empty_store(root)
            first = store / "states" / "first.json"
            first.write_bytes(b"{}")
            os.link(first, store / "states" / "second.json")
            with self.assertRaises(freeze.TDG8RCV3FreezeError):
                freeze._walk_tree(store)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            store = self._empty_store(root)
            (store / "foreign").mkdir()
            with self.assertRaises(freeze.TDG8RCV3FreezeError):
                freeze._walk_tree(store)

    def test_destination_absence_is_nofollow_and_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "runs" / "fgc-2-sf1").mkdir(parents=True)
            self.assertTrue(
                freeze._path_absent_nofollow(root, freeze.DESTINATION_WRAPPER_PATH)
            )
            os.symlink(root / "elsewhere", root / freeze.DESTINATION_WRAPPER_PATH)
            with self.assertRaises(freeze.TDG8RCV3FreezeError):
                freeze._path_absent_nofollow(root, freeze.DESTINATION_STORE_PATH)

    def _committed_execution_image(self, repository: Path) -> None:
        paths = (
            freeze.CONFIG_PATH,
            freeze.RESULT_PATH,
            bootstrap.FREEZE_MODULE_PATH,
            freeze.BOOTSTRAP_RUNTIME_PATH,
            freeze.BOOTSTRAP_SCRIPT_PATH,
        )
        for relative in paths:
            destination = repository / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / relative, destination)
        (repository / "tracked-sentinel.txt").write_text("clean\n", "ascii")
        subprocess.run(["git", "init", "-q"], cwd=repository, check=True)
        subprocess.run(["git", "add", "."], cwd=repository, check=True)
        subprocess.run(
            [
                "git", "-c", "user.name=RCV3 Test", "-c",
                "user.email=rcv3@example.invalid", "commit", "-q", "-m", "fixture",
            ],
            cwd=repository,
            check=True,
        )

    def test_bootstrap_execution_image_requires_committed_clean_tracked_state(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            repository = Path(temporary)
            self._committed_execution_image(repository)
            observed = bootstrap.preauthenticate_execution_image(repository)
            self.assertEqual(observed["runtime_sha256"], freeze.BOOTSTRAP_RUNTIME_SHA256)
            self.assertEqual(observed["script_sha256"], freeze.BOOTSTRAP_SCRIPT_SHA256)
            (repository / "tracked-sentinel.txt").write_text("dirty\n", "ascii")
            with self.assertRaisesRegex(
                bootstrap.TDG8RCV3BootstrapError, "tracked Git state is dirty"
            ):
                bootstrap.preauthenticate_execution_image(repository)

    def test_mutated_project_module_is_rejected_before_it_can_execute(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            repository = Path(temporary)
            self._committed_execution_image(repository)
            module = repository / bootstrap.FREEZE_MODULE_PATH
            module.write_text(
                module.read_text("utf-8") + "\nraise RuntimeError('must not execute')\n",
                "utf-8",
            )
            with self.assertRaisesRegex(
                bootstrap.TDG8RCV3BootstrapError, "tracked Git state is dirty"
            ):
                bootstrap.preauthenticate_execution_image(repository)

    def test_bootstrap_calls_only_fixed_production_installer_with_external_authority(self) -> None:
        result_sha = sha256(self.result).hexdigest()
        image = {
            "head": "a" * 40,
            "config": freeze.parse_config(self.config),
            "config_raw": self.config,
            "result_raw": self.result,
            "config_sha256": sha256(self.config).hexdigest(),
            "result_sha256": result_sha,
            "freeze_module_sha256": "c" * 64,
            "runtime_sha256": freeze.BOOTSTRAP_RUNTIME_SHA256,
            "script_sha256": freeze.BOOTSTRAP_SCRIPT_SHA256,
        }
        installed = SimpleNamespace(
            projection_id=freeze.PROJECTION_ID,
            receipt_relative="runs/fgc-2-sf1/tdg8-rcv3/bootstrap-receipt.json",
            receipt_sha256="b" * 64,
            anchor_checkpoint_sha256=freeze.GENERATION9_CHECKPOINT_SHA256,
            installed_now=True,
        )
        with (
            patch.object(
                bootstrap,
                "preauthenticate_execution_image",
                side_effect=(image, image),
            ),
            patch.object(freeze, "destination_absence_status", return_value=(True, True)),
            patch.object(runtime, "install_recovery_fork", return_value=installed) as call,
        ):
            result = bootstrap.bootstrap(ROOT)
        self.assertTrue(result["installed_now"])
        self.assertFalse(result["bounded_execution_authorized"])
        self.assertFalse(result["candidate_branch_opened"])
        authority = call.call_args.kwargs["authority"]
        self.assertEqual(authority.authority_commit_sha, "a" * 40)
        self.assertEqual(authority.frz1_result_sha256, result_sha)
        self.assertEqual(authority.runtime_sha256, freeze.BOOTSTRAP_RUNTIME_SHA256)
        self.assertEqual(authority.bootstrap_script_sha256, freeze.BOOTSTRAP_SCRIPT_SHA256)
        self.assertIs(call.call_args.args[1], runtime.PRODUCTION_SPEC)

    def test_bootstrap_imports_no_candidate_or_evolution_runtime(self) -> None:
        tree = ast.parse((ROOT / freeze.BOOTSTRAP_SCRIPT_PATH).read_text("utf-8"))
        imports: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imports.append(node.module or "")
        joined = "\n".join(imports).lower()
        self.assertNotIn("recursive_horizons", joined)
        for forbidden in (
            "fgcqr", "fgc_qr", "sgb", "progression_attempt", "numerical_engine",
        ):
            self.assertNotIn(forbidden, joined)

    def test_bootstrap_main_requires_isolated_no_bytecode_launch(self) -> None:
        process = subprocess.run(
            [sys.executable, str(ROOT / freeze.BOOTSTRAP_SCRIPT_PATH)],
            cwd=ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertNotEqual(process.returncode, 0)
        self.assertIn(b"requires isolated no-bytecode launch", process.stderr)


if __name__ == "__main__":
    unittest.main()
