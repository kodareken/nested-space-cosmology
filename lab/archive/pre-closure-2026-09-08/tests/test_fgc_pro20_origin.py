"""Synthetic, store-blind tests for the read-only PRO20 origin capture."""

from __future__ import annotations

import ast
from dataclasses import FrozenInstanceError, replace
from hashlib import sha256
import inspect
import json
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.evidence_io import (  # noqa: E402
    canonical_json_bytes,
    load_canonical_json,
)
from recursive_horizons.fgc.evolution.hlt16_member_codec import (  # noqa: E402
    ARRAY_NAMES,
    GEN0_PROTOCOL,
)
from recursive_horizons.fgc.evolution.hlt16_state_store import (  # noqa: E402
    STATE_SCHEMA,
    canonical as hlt16_canonical,
    digest as hlt16_digest,
    encode_payload,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.protocol_v17 import (  # noqa: E402
    MEMBER_KEYS,
    STATE_ARRAY_NAMES,
)
from recursive_horizons.fgc.evolution import pro20_origin as origin  # noqa: E402
from recursive_horizons.fgc.evolution.pro20_origin import (  # noqa: E402
    Pro20OriginError,
    capture_pro20_historical_origin,
    cleanup_baseline_tree_digest,
    compare_bridge_member_to_original,
    full_cursor_sha256,
    inventory_source_tree,
    physical_array_sha256,
    revalidate_pro20_origin_capture,
)


def _digest(label: str) -> str:
    return sha256(label.encode("ascii")).hexdigest()


def _write(path: Path, payload: bytes) -> bytes:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return payload


def _canonical(value: object) -> bytes:
    return canonical_json_bytes(value)


def _time() -> dict[str, str]:
    return {
        "rational": origin.ORIGIN_TIME_RATIONAL,
        "binary64_hex": origin.ORIGIN_TIME_HEX,
    }


def _counters(key: str) -> dict[str, int]:
    return dict(origin._OWNER_PINS.member_counters[key])


def _arrays(key: str) -> dict[str, np.ndarray]:
    seed = float(MEMBER_KEYS.index(key) + 1)
    labels = 2
    events = 2
    points = 2
    return {
        "u": np.full((points, 6), seed, dtype="<f8"),
        "p": np.full((points, 6), seed + 0.25, dtype="<f8"),
        "q": np.full((points, 6), seed + 0.5, dtype="<f8"),
        "grid_coordinates": np.linspace(0.0, 1.0, points, dtype="<f8"),
        "tracer_labels": np.arange(labels, dtype="<f8"),
        "tracer_positions": np.full((labels,), seed + 0.75, dtype="<f8"),
        "tracer_proper_times": np.full((labels,), seed + 1.0, dtype="<f8"),
        "event_proper_times": np.full((events, labels), seed + 1.25, dtype="<f8"),
        "event_fields": np.full((events, labels, 6), seed + 1.5, dtype="<f8"),
    }


def _state_object(key: str, arrays: dict[str, np.ndarray]) -> dict[str, object]:
    counters = _counters(key)
    physical = []
    for name in STATE_ARRAY_NAMES:
        physical.append(
            {
                "logical_name": name,
                "source_bundle_key": f"bundle-{key}",
                "npz_storage_key": name,
                "dtype": "<f8",
                "layout": "C",
                "shape": list(arrays[name].shape),
                "little_endian_c_bytes_sha256": physical_array_sha256(arrays[name]),
            }
        )
    return {
        "schema_id": "FGC-2-SF1-PROTO17-state-object-v1",
        "member_key": key,
        "method": key.split("-", 1)[0],
        "point_count": 2,
        "coordinate_time": _time(),
        "accepted_boundary_time": _time(),
        "step_index": counters["step_index"],
        "transaction_serial": counters["transaction_serial"],
        "input_hash": _digest(f"input-{key}"),
        "source_retry_count": counters["source_retry_count"],
        "CFL_retry_count": counters["CFL_retry_count"],
        "runtime_monitor_state": {
            "accepted_stage_count": counters["step_index"],
            "first_failed_premise": None,
            "first_failed_time": None,
            "first_failed_transaction_serial": None,
            "last_accepted_time": 1.4375,
            "last_transaction_serial": counters["transaction_serial"] - 1,
            "member_key": key,
        },
        "causal_state": {
            "accepted_time": 0,
            "accumulated_characteristic_distance": 0,
            "previous_speed_upper": 0,
            "member_key": key,
        },
        "tracer_state": {"label_count": 2},
        "event_history_state": {"sample_count": 2},
        "physical_arrays": physical,
        "physical_state_sha256": array_content_sha256(
            arrays["u"], arrays["p"], arrays["q"]
        ),
        "restart_payload_sha256": _digest(f"restart-{key}"),
    }


def _cursor(key: str, state_content_id: str) -> dict[str, object]:
    counters = _counters(key)
    return {
        "member_key": key,
        "accepted_boundary_time": _time(),
        "accepted_state_sha256": state_content_id,
        "committed_common_event_index": origin.ORIGIN_EVENT,
        "cursor_chain_sha256": _digest(f"chain-{key}"),
        "previous_step_index": counters["step_index"],
        "previous_transaction_serial": counters["transaction_serial"],
    }


def _metadata(key: str, cursor: dict[str, object], state: dict[str, object]) -> dict[str, object]:
    return {
        "runtime_identity": {
            "member_key": key,
            "generation": 0,
            "accepted_time": _time(),
        },
        "provenance_cursor": {
            "source_protocol": GEN0_PROTOCOL,
            "bridge": "authenticated_GEN0_predecessor",
            "payload": cursor,
        },
        "counters": {
            "step_index": state["step_index"],
            "transaction_serial": state["transaction_serial"],
            "source_retry_count": state["source_retry_count"],
            "CFL_retry_count": state["CFL_retry_count"],
        },
        "template": {"input_hash": state["input_hash"]},
        "monitor": state["runtime_monitor_state"],
        "causal": state["causal_state"],
    }


def _hlt16_descriptor(arrays: dict[str, np.ndarray], metadata: dict[str, object]):
    payload, semantic, manifest = encode_payload({name: arrays[name] for name in ARRAY_NAMES})
    raw_archive = sha256(payload).hexdigest()
    body = {
        "schema": STATE_SCHEMA,
        "semantic_sha256": semantic,
        "raw_archive_sha256": raw_archive,
        "arrays": list(manifest),
        "metadata": metadata,
    }
    address = hlt16_digest(body)
    descriptor = {**body, "descriptor_sha256": address}
    return address, hlt16_canonical(descriptor), payload, semantic, raw_archive


class LoadedMember:
    def __init__(
        self,
        *,
        key: str,
        arrays: dict[str, np.ndarray],
        metadata: dict[str, object],
        descriptor: str,
        payload_raw: str,
        semantic: str,
        descriptor_obj: dict[str, object],
    ) -> None:
        self.arrays = arrays
        self.metadata = metadata
        self.descriptor = {**descriptor_obj, "metadata": metadata}
        self.descriptor_sha256 = descriptor
        self.raw_archive_sha256 = payload_raw
        self.semantic_sha256 = semantic


class FakeStore:
    def __init__(self, path: Path, checkpoint, loaded: dict[str, LoadedMember]) -> None:
        self.path = path
        self.checkpoint = checkpoint
        self.loaded = loaded
        self.calls: list[tuple[object, ...]] = []

    def authenticated_checkpoint_at_generation(self, generation: int, **kwargs):
        self.calls.append(("authenticated_checkpoint_at_generation", generation, kwargs))
        return self.checkpoint

    def load_state(self, descriptor: str):
        self.calls.append(("load_state", descriptor))
        return self.loaded[descriptor]

    def recover(self, *args, **kwargs):
        raise AssertionError("recover")

    def persist_snapshot(self, *args, **kwargs):
        raise AssertionError("persist_snapshot")

    def publish_journal(self, *args, **kwargs):
        raise AssertionError("publish_journal")

    def publish_checkpoint(self, *args, **kwargs):
        raise AssertionError("publish_checkpoint")

    def adopt_generation_zero(self, *args, **kwargs):
        raise AssertionError("adopt_generation_zero")

    def inspect_recovery(self, *args, **kwargs):
        raise AssertionError("inspect_recovery")


class FakeCheckpoint:
    def __init__(self, *, generation: int, sha256_value: str, members) -> None:
        self.generation = generation
        self.sha256 = sha256_value
        self.event = origin.ORIGIN_EVENT
        self.members = members


class FakeMember:
    def __init__(self, descriptor_sha256: str) -> None:
        self.descriptor_sha256 = descriptor_sha256


BRIDGE_ID = _digest("generation-1-bridge")


def _install_world(root: Path, *, extras: dict[str, bytes] | None = None):
    members = {}
    loaded_by_descriptor: dict[str, LoadedMember] = {}
    checkpoint_members = {}
    for key in MEMBER_KEYS:
        arrays = _arrays(key)
        state = _state_object(key, arrays)
        raw = _canonical(state)
        content_id = sha256(raw).hexdigest()
        relative = f"{origin.SOURCE_STORE_RELATIVE}/states/{content_id}.json"
        _write(root / relative, raw)
        cursor = _cursor(key, content_id)
        metadata = _metadata(key, cursor, state)
        descriptor, descriptor_bytes, payload, semantic, payload_raw = _hlt16_descriptor(
            arrays, metadata
        )
        _write(
            root / f"{origin.SOURCE_STORE_RELATIVE}/states/{descriptor}.json",
            descriptor_bytes,
        )
        _write(
            root / f"{origin.SOURCE_STORE_RELATIVE}/payloads/{semantic}.npz",
            payload,
        )
        descriptor_obj = json.loads(descriptor_bytes.decode("ascii"))
        loaded = LoadedMember(
            key=key,
            arrays=arrays,
            metadata=metadata,
            descriptor=descriptor,
            payload_raw=payload_raw,
            semantic=semantic,
            descriptor_obj=descriptor_obj,
        )
        members[key] = {
            "arrays": arrays,
            "state": state,
            "cursor": cursor,
            "content_id": content_id,
            "path": relative,
            "raw": raw,
            "descriptor": descriptor,
            "descriptor_bytes": descriptor_bytes,
            "payload_bytes": payload,
            "payload_raw": payload_raw,
            "semantic": semantic,
            "loaded": loaded,
        }
        loaded_by_descriptor[descriptor] = loaded
        checkpoint_members[key] = FakeMember(descriptor)

    checkpoint_body = {
        "checkpoint_sha256": _digest("hlt15-checkpoint-content"),
        "committed_common_event_index": origin.ORIGIN_EVENT,
        "states": {key: members[key]["content_id"] for key in MEMBER_KEYS},
        "cursors": {key: members[key]["cursor"] for key in MEMBER_KEYS},
    }
    checkpoint_raw = _canonical(checkpoint_body)
    checkpoint_path = (
        f"{origin.SOURCE_STORE_RELATIVE}/checkpoints/"
        f"00000000000000000000-{checkpoint_body['checkpoint_sha256']}.json"
    )
    _write(root / checkpoint_path, checkpoint_raw)
    receipt_body = {
        "receipt_sha256": _digest("hlt15-receipt-content"),
        "checkpoint_sha256": checkpoint_body["checkpoint_sha256"],
    }
    receipt_raw = _canonical(receipt_body)
    receipt_path = f"{origin.SOURCE_STORE_RELATIVE}/receipts/generation-zero.json"
    _write(root / receipt_path, receipt_raw)

    gen1_body = {
        "authorization_commit": "a" * 40,
        "campaign_id": "SYNTHETIC",
        "checkpoint_sha256": BRIDGE_ID,
        "disposition": "nonterminal",
        "event": origin.ORIGIN_EVENT,
        "generation": origin.SELECTED_GENERATION,
        "journal_sequence": 0,
        "journal_tip_sha256": "b" * 64,
        "members": {
            key: {"descriptor_sha256": members[key]["descriptor"]}
            for key in MEMBER_KEYS
        },
        "parent_sha256": "c" * 64,
        "plan_sha256": "d" * 64,
        "protocol": "FGC-2-SF1-PROTO18",
        "target": {"binary64_hex": (3 / 2).hex(), "rational": "3/2"},
        "terminal": None,
    }
    gen1_raw = _canonical(gen1_body)
    gen1_path = (
        f"{origin.SOURCE_STORE_RELATIVE}/checkpoints/"
        f"{origin.SELECTED_GENERATION:020d}-{BRIDGE_ID}.json"
    )
    extra_payloads = extras if extras is not None else {
        gen1_path: gen1_raw,
        f"{origin.SOURCE_STORE_RELATIVE}/journal/00000000000000000000-{_digest('journal')}.journal": b"journal",
    }
    if extras is None:
        extra_payloads = {
            gen1_path: gen1_raw,
            f"{origin.SOURCE_STORE_RELATIVE}/journal/00000000000000000000-{_digest('journal')}.journal": b"journal",
        }
    else:
        extra_payloads = dict(extras)
        extra_payloads.setdefault(gen1_path, gen1_raw)
    for relative, payload in extra_payloads.items():
        _write(root / relative, payload)

    leaf_rows = [
        {
            "content_sha256": checkpoint_body["checkpoint_sha256"],
            "path": checkpoint_path,
            "raw_sha256": sha256(checkpoint_raw).hexdigest(),
        },
        {
            "content_sha256": receipt_body["receipt_sha256"],
            "path": receipt_path,
            "raw_sha256": sha256(receipt_raw).hexdigest(),
        },
    ]
    for key in MEMBER_KEYS:
        leaf_rows.append(
            {
                "content_sha256": members[key]["content_id"],
                "path": members[key]["path"],
                "raw_sha256": sha256(members[key]["raw"]).hexdigest(),
            }
        )
    ordered = sorted(leaf_rows, key=lambda item: item["path"])
    tree = sha256(_canonical(ordered)).hexdigest()
    compact = {
        "artifact_id": "FGC-1-PRO18-PREF27",
        "artifact_payload": {
            "store_leaves": ordered,
            "store_tree_sha256": tree,
        },
    }
    compact_raw = (json.dumps(compact, indent=2, sort_keys=True) + "\n").encode("utf-8")
    _write(root / origin.PREF27_COMPACT_RELATIVE, compact_raw)

    store = FakeStore(
        root / origin.SOURCE_STORE_RELATIVE,
        FakeCheckpoint(
            generation=origin.SELECTED_GENERATION,
            sha256_value=BRIDGE_ID,
            members=checkpoint_members,
        ),
        loaded_by_descriptor,
    )
    inventory = inventory_source_tree(root, origin.SOURCE_TREE_RELATIVE)
    pins = origin._OwnerPins(
        selected_method=origin.SELECTED_METHOD,
        selected_generation=origin.SELECTED_GENERATION,
        origin_event=origin.ORIGIN_EVENT,
        origin_time_rational=origin.ORIGIN_TIME_RATIONAL,
        origin_time_hex=origin.ORIGIN_TIME_HEX,
        generation1_bridge_content_id=BRIDGE_ID,
        expanded_leaf_count=inventory.leaf_count,
        expanded_byte_count=inventory.byte_count,
        expanded_cleanup_digest=inventory.cleanup_digest,
        pref27_compact_raw_sha256=sha256(compact_raw).hexdigest(),
        pref27_tree_sha256=tree,
        hlt15_checkpoint_content_id=checkpoint_body["checkpoint_sha256"],
        hlt15_checkpoint_raw_sha256=sha256(checkpoint_raw).hexdigest(),
        hlt15_receipt_content_id=receipt_body["receipt_sha256"],
        hlt15_receipt_raw_sha256=sha256(receipt_raw).hexdigest(),
        member_counters=origin._OWNER_PINS.member_counters,
    )
    return members, store, pins, inventory, receipt_body, checkpoint_body


class Pro20OriginHelperTests(unittest.TestCase):
    def test_cleanup_digest_changes_with_leaf_hash_or_path(self) -> None:
        first = cleanup_baseline_tree_digest(
            (("runs/a", 4, _digest("a")), ("runs/b", 5, _digest("b")))
        )
        changed_hash = cleanup_baseline_tree_digest(
            (("runs/a", 4, _digest("A")), ("runs/b", 5, _digest("b")))
        )
        changed_path = cleanup_baseline_tree_digest(
            (("runs/c", 4, _digest("a")), ("runs/b", 5, _digest("b")))
        )
        self.assertNotEqual(first, changed_hash)
        self.assertNotEqual(first, changed_path)
        self.assertTrue(first.startswith(tuple("0123456789abcdef")))

    def test_physical_array_hash_is_not_content_id_or_shape_prefixed_hash(self) -> None:
        array = np.arange(12, dtype="<f8").reshape(2, 6)
        physical = physical_array_sha256(array)
        shaped = array_content_sha256(array)
        self.assertNotEqual(physical, shaped)
        self.assertEqual(physical, sha256(array.tobytes(order="C")).hexdigest())

    def test_full_cursor_hash_is_not_cursor_chain_id(self) -> None:
        cursor = _cursor("RK4-2049", _digest("state"))
        self.assertNotEqual(full_cursor_sha256(cursor), cursor["cursor_chain_sha256"])


class Pro20OriginCaptureTests(unittest.TestCase):
    def test_public_capture_has_no_synthetic_success_bypass(self) -> None:
        signature = inspect.signature(capture_pro20_historical_origin)
        self.assertEqual(list(signature.parameters), ["repository_root"])
        tree = ast.parse(
            (ROOT / "src/recursive_horizons/fgc/evolution/pro20_origin.py").read_text(
                "utf-8"
            )
        )
        names: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                names.append(node.module or "")
            elif isinstance(node, ast.FunctionDef):
                names.append(node.name)
        joined = " ".join(names)
        self.assertNotIn("verify_generation_zero_store", joined)
        self.assertNotIn("run_fgc_gr0_calibration", joined)
        self.assertNotIn("proto19_gr0_static_factory", joined)
        self.assertNotIn("reconstruct_gr0_members", joined)
        source = inspect.getsource(origin)
        self.assertNotIn("repair_stages=True", source)
        capture_source = inspect.getsource(capture_pro20_historical_origin)
        self.assertNotIn(".recover(", capture_source)
        self.assertNotIn(".publish_", capture_source)
        self.assertNotIn(".adopt_generation_zero(", capture_source)
        self.assertNotIn("success", inspect.signature(capture_pro20_historical_origin).parameters)

    def _capture(self, root: Path, store: FakeStore, pins: origin._OwnerPins):
        with (
            patch.object(origin, "_OWNER_PINS", pins),
            patch.object(origin, "_open_historical_store", return_value=store),
            patch(
                "recursive_horizons.fgc.evolution.proto18_pref27_binder.verify_generation_zero_store",
                side_effect=AssertionError("exact-eight checker"),
            ) as exact_eight,
            patch(
                "recursive_horizons.fgc.evolution.proto19_progression_inputs.reconstruct_gr0_members",
                side_effect=AssertionError("reconstruct"),
            ),
            patch(
                "recursive_horizons.fgc.evolution.proto19_gr0_static_factory.build_static_gr0_shells",
                side_effect=AssertionError("gr0 shells"),
            ),
        ):
            captured = capture_pro20_historical_origin(root)
        exact_eight.assert_not_called()
        return captured

    def test_successful_synthetic_capture_pins_required_fields(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            members, store, pins, inventory, receipt_body, checkpoint_body = _install_world(root)
            captured = self._capture(root, store, pins)
            self.assertEqual(captured.selected_method, "TDG11-IMP1")
            self.assertEqual(captured.origin_event, 23)
            self.assertEqual(captured.origin_time["rational"], "23/16")
            self.assertEqual(captured.origin_time["binary64_hex"], origin.ORIGIN_TIME_HEX)
            self.assertEqual(
                captured.generation1_bridge_content_id, pins.generation1_bridge_content_id
            )
            self.assertEqual(
                tuple(item.member_key for item in captured.members), MEMBER_KEYS
            )
            self.assertEqual(captured.source_inventory.leaf_count, inventory.leaf_count)
            self.assertGreater(captured.source_inventory.leaf_count, 8)
            self.assertEqual(captured.pref27_tree_sha256, pins.pref27_tree_sha256)
            self.assertEqual(
                captured.pref27_compact.raw_sha256, pins.pref27_compact_raw_sha256
            )
            receipt = next(
                leaf for leaf in captured.original_leaves if leaf.path.endswith(
                    "/receipts/generation-zero.json"
                )
            )
            self.assertEqual(receipt.content_id, receipt_body["receipt_sha256"])
            self.assertNotEqual(receipt.content_id, receipt_body["checkpoint_sha256"])
            self.assertEqual(receipt_body["checkpoint_sha256"], checkpoint_body["checkpoint_sha256"])
            first = captured.members[0]
            self.assertNotEqual(first.descriptor_content_id, first.payload_raw_sha256)
            self.assertNotEqual(first.full_cursor_sha256, first.cursor_chain_sha256)
            self.assertNotEqual(
                first.physical_array_sha256["u"],
                array_content_sha256(members["RK4-2049"]["arrays"]["u"]),
            )
            self.assertIn("grid_coordinates", first.physical_array_sha256)
            self.assertIn("tracer_labels", first.physical_array_sha256)
            descriptor_bytes, payload_bytes = captured.retained_bridge_payload("RK4-2049")
            self.assertEqual(descriptor_bytes, members["RK4-2049"]["descriptor_bytes"])
            self.assertEqual(payload_bytes, members["RK4-2049"]["payload_bytes"])
            self.assertEqual(
                store.calls[0][0], "authenticated_checkpoint_at_generation"
            )
            self.assertEqual(store.calls[0][1], 1)
            self.assertEqual(store.calls[0][2], {})
            self.assertTrue(all(call[0] != "recover" for call in store.calls))
            projected = captured.fresh_physical_arrays("RK4-2049")
            np.testing.assert_array_equal(
                projected["u"], members["RK4-2049"]["arrays"]["u"]
            )
            bridge = captured.fresh_bridge_arrays("RK4-2049")
            np.testing.assert_array_equal(
                bridge["grid_coordinates"],
                members["RK4-2049"]["arrays"]["grid_coordinates"],
            )
            with patch.object(origin, "_OWNER_PINS", pins):
                revalidate_pro20_origin_capture(captured)

    def test_changed_original_leaf_hash_is_rejected(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            _members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            receipt = root / f"{origin.SOURCE_STORE_RELATIVE}/receipts/generation-zero.json"
            receipt.write_bytes(receipt.read_bytes() + b" ")
            with self.assertRaisesRegex(Pro20OriginError, "raw hash differs"):
                self._capture(root, store, pins)

    def test_selected_later_generation_is_rejected(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            _members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            store.checkpoint.generation = 2
            with self.assertRaisesRegex(Pro20OriginError, "generation 1"):
                self._capture(root, store, pins)

    def test_swapped_members_are_rejected(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            first = members["RK4-2049"]["descriptor"]
            second = members["RK4-4097"]["descriptor"]
            store.loaded[first], store.loaded[second] = (
                store.loaded[second],
                store.loaded[first],
            )
            with self.assertRaisesRegex(Pro20OriginError, "differs from retained payload"):
                self._capture(root, store, pins)

    def test_altered_inherited_counters_are_rejected(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            path = root / (
                f"{origin.SOURCE_STORE_RELATIVE}/states/"
                f"{members['RK4-2049']['descriptor']}.json"
            )
            descriptor = json.loads(path.read_bytes().decode("ascii"))
            descriptor["metadata"]["counters"]["step_index"] += 1
            path.write_bytes(_canonical(descriptor))
            with self.assertRaisesRegex(Pro20OriginError, "descriptor address differs"):
                self._capture(root, store, pins)

    def test_altered_original_cursor_is_rejected(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            first = members["RK4-2049"]
            loaded = LoadedMember(
                key="RK4-2049",
                arrays=first["arrays"],
                metadata=first["loaded"].metadata,
                descriptor=first["descriptor"],
                payload_raw=first["payload_raw"],
                semantic=first["semantic"],
                descriptor_obj=json.loads(first["descriptor_bytes"].decode("ascii")),
            )
            payload = dict(loaded.metadata["provenance_cursor"]["payload"])
            payload["cursor_chain_sha256"] = _digest("mutated-cursor")
            loaded.metadata = dict(loaded.metadata)
            loaded.metadata["provenance_cursor"] = dict(
                loaded.metadata["provenance_cursor"]
            )
            loaded.metadata["provenance_cursor"]["payload"] = payload
            loaded.descriptor = {**loaded.descriptor, "metadata": loaded.metadata}
            with self.assertRaisesRegex(Pro20OriginError, "original cursor differs"):
                compare_bridge_member_to_original(
                    "RK4-2049",
                    loaded=loaded,
                    original_state=first["state"],
                    original_cursor=first["cursor"],
                    original_state_content_id=first["content_id"],
                    descriptor_bytes=first["descriptor_bytes"],
                    payload_bytes=first["payload_bytes"],
                )

    def test_altered_event_arrays_are_rejected(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            loaded = members["RK4-2049"]["loaded"]
            loaded.arrays = dict(loaded.arrays)
            mutated = loaded.arrays["event_fields"].copy()
            mutated[0, 0, 0] += 1.0
            loaded.arrays["event_fields"] = mutated
            with self.assertRaisesRegex(
                Pro20OriginError, "event_fields|retained payload"
            ):
                self._capture(root, store, pins)

    def test_absent_source_tree_is_rejected(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            with self.assertRaisesRegex(Pro20OriginError, "absent"):
                capture_pro20_historical_origin(root)

    def test_partial_original_leaves_are_rejected(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            _members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            receipt = root / f"{origin.SOURCE_STORE_RELATIVE}/receipts/generation-zero.json"
            receipt.unlink()
            with self.assertRaisesRegex(Pro20OriginError, "missing original PREF27 leaf"):
                self._capture(root, store, pins)

    def test_forked_extra_source_path_is_rejected(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            _members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            _write(root / f"{origin.SOURCE_TREE_RELATIVE}/fork/extra.json", b"fork")
            with self.assertRaisesRegex(Pro20OriginError, "forked or has extra paths"):
                self._capture(root, store, pins)

    def test_unsafe_extra_symlink_is_rejected(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            _members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            target = root / f"{origin.SOURCE_STORE_RELATIVE}/receipts/generation-zero.json"
            link = root / f"{origin.SOURCE_STORE_RELATIVE}/states/link.json"
            os.symlink(target, link)
            with self.assertRaisesRegex(Pro20OriginError, "unsafe"):
                self._capture(root, store, pins)

    def test_no_repair_or_write_calls_on_store(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            _members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            captured = self._capture(root, store, pins)
            self.assertEqual(
                [call[0] for call in store.calls],
                ["authenticated_checkpoint_at_generation"]
                + ["load_state"] * len(MEMBER_KEYS),
            )
            self.assertEqual(captured.origin_event, 23)

    def test_returned_identity_is_immutable(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            captured = self._capture(root, store, pins)
            before = captured.captured_bytes
            with self.assertRaises(FrozenInstanceError):
                captured.generation1_bridge_content_id = "0" * 64
            with self.assertRaises(TypeError):
                captured.members[0].physical_array_sha256["u"] = "0" * 64
            projected = captured.fresh_physical_arrays("RK4-2049")
            projected["u"][0, 0] = 99.0
            self.assertEqual(captured.captured_bytes, before)
            self.assertEqual(
                captured.members[0].physical_array_sha256["u"],
                physical_array_sha256(members["RK4-2049"]["arrays"]["u"]),
            )
            self.assertEqual(
                load_canonical_json(captured.captured_bytes),
                captured.identity_record(),
            )

    def test_compare_helper_rejects_counter_mismatch(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            members, _store, _pins, _inventory, _receipt, _checkpoint = _install_world(root)
            first = members["RK4-2049"]
            loaded = LoadedMember(
                key="RK4-2049",
                arrays=first["arrays"],
                metadata=first["loaded"].metadata,
                descriptor=first["descriptor"],
                payload_raw=first["payload_raw"],
                semantic=first["semantic"],
                descriptor_obj=json.loads(first["descriptor_bytes"].decode("ascii")),
            )
            loaded.metadata = dict(loaded.metadata)
            loaded.metadata["counters"] = dict(loaded.metadata["counters"])
            loaded.metadata["counters"]["step_index"] += 1
            loaded.descriptor = {**loaded.descriptor, "metadata": loaded.metadata}
            with self.assertRaisesRegex(Pro20OriginError, "step_index"):
                compare_bridge_member_to_original(
                    "RK4-2049",
                    loaded=loaded,
                    original_state=first["state"],
                    original_cursor=first["cursor"],
                    original_state_content_id=first["content_id"],
                    descriptor_bytes=first["descriptor_bytes"],
                    payload_bytes=first["payload_bytes"],
                )

    def test_compare_helper_rejects_swapped_event_history(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            members, _store, _pins, _inventory, _receipt, _checkpoint = _install_world(root)
            first = members["RK4-2049"]
            second = members["RK4-4097"]
            arrays = dict(first["arrays"])
            arrays["event_fields"] = second["arrays"]["event_fields"]
            with self.assertRaisesRegex(Pro20OriginError, "event_fields"):
                compare_bridge_member_to_original(
                    "RK4-2049",
                    loaded=LoadedMember(
                        key="RK4-2049",
                        arrays=arrays,
                        metadata=first["loaded"].metadata,
                        descriptor=first["descriptor"],
                        payload_raw=first["payload_raw"],
                        semantic=first["semantic"],
                        descriptor_obj=json.loads(first["descriptor_bytes"].decode("ascii")),
                    ),
                    original_state=first["state"],
                    original_cursor=first["cursor"],
                    original_state_content_id=first["content_id"],
                    descriptor_bytes=first["descriptor_bytes"],
                    payload_bytes=first["payload_bytes"],
                )

    def test_receipt_content_id_is_not_embedded_checkpoint_id(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            _members, store, pins, _inventory, receipt_body, checkpoint_body = _install_world(root)
            captured = self._capture(root, store, pins)
            receipt = next(
                leaf
                for leaf in captured.original_leaves
                if origin._is_receipt_path(leaf.path)
            )
            self.assertEqual(receipt.content_id, pins.hlt15_receipt_content_id)
            self.assertEqual(receipt.content_id, receipt_body["receipt_sha256"])
            self.assertEqual(checkpoint_body["checkpoint_sha256"], receipt_body["checkpoint_sha256"])
            self.assertNotEqual(receipt.content_id, checkpoint_body["checkpoint_sha256"])
            self.assertIn("checkpoint_sha256", load_canonical_json(receipt.raw_bytes))
            self.assertIn("receipt_sha256", load_canonical_json(receipt.raw_bytes))

    def test_compact_modification_is_rejected(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            _members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            compact = root / origin.PREF27_COMPACT_RELATIVE
            compact.write_bytes(compact.read_bytes() + b"\n")
            with self.assertRaisesRegex(Pro20OriginError, "PREF27 compact"):
                self._capture(root, store, pins)

    def test_coordinated_rehash_replace_is_not_a_verified_origin(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            _members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            captured = self._capture(root, store, pins)
            forged_members = []
            for member in captured.members:
                if member.member_key == "RK4-2049":
                    member = replace(
                        member,
                        step_index=0,
                        physical_array_shapes={
                            **member.physical_array_shapes,
                            "u": (12,),
                        },
                    )
                forged_members.append(member)
            identity = origin._identity_record(
                selected_method=captured.selected_method,
                origin_event=99,
                origin_time={"rational": "3/2", "binary64_hex": (3 / 2).hex()},
                generation1_bridge_content_id="0" * 64,
                generation1_checkpoint=captured.generation1_checkpoint,
                pref27_compact=captured.pref27_compact,
                pref27_tree_sha256=captured.pref27_tree_sha256,
                original_leaves=captured.original_leaves,
                source_inventory=captured.source_inventory,
                members=tuple(forged_members),
            )
            with patch.object(origin, "_OWNER_PINS", pins):
                with self.assertRaisesRegex(Pro20OriginError, "retained origin bytes"):
                    replace(
                        captured,
                        origin_event=99,
                        origin_time={"rational": "3/2", "binary64_hex": (3 / 2).hex()},
                        generation1_bridge_content_id="0" * 64,
                        members=tuple(forged_members),
                        captured_bytes=canonical_json_bytes(identity),
                    )

    def test_descriptor_and_payload_mismatch_are_rejected(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            _members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            captured = self._capture(root, store, pins)
            first, second = captured.members[0], captured.members[1]
            swapped = replace(
                first,
                descriptor_bytes=second.descriptor_bytes,
                descriptor_content_id=second.descriptor_content_id,
                payload_bytes=second.payload_bytes,
                payload_raw_sha256=second.payload_raw_sha256,
            )
            members = (swapped,) + captured.members[1:]
            identity = origin._identity_record(
                selected_method=captured.selected_method,
                origin_event=captured.origin_event,
                origin_time=captured.origin_time,
                generation1_bridge_content_id=captured.generation1_bridge_content_id,
                generation1_checkpoint=captured.generation1_checkpoint,
                pref27_compact=captured.pref27_compact,
                pref27_tree_sha256=captured.pref27_tree_sha256,
                original_leaves=captured.original_leaves,
                source_inventory=captured.source_inventory,
                members=members,
            )
            with patch.object(origin, "_OWNER_PINS", pins):
                with self.assertRaisesRegex(
                    Pro20OriginError,
                    "generation-1 member|member identity|physical-array",
                ):
                    replace(
                        captured,
                        members=members,
                        captured_bytes=canonical_json_bytes(identity),
                    )

    def test_fresh_hlt16_snapshot_is_rebuilt_from_retained_bytes(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            _members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            captured = self._capture(root, store, pins)
            first = captured.fresh_hlt16_snapshot("RK4-2049")
            expected = captured.fresh_bridge_arrays("RK4-2049")
            self.assertEqual(tuple(first.arrays), origin.ARRAY_NAMES)
            for name in origin.ARRAY_NAMES:
                np.testing.assert_array_equal(first.arrays[name], expected[name])
            first.arrays["u"].flat[0] += 1
            second = captured.fresh_hlt16_snapshot("RK4-2049")
            self.assertNotEqual(first.arrays["u"].flat[0], second.arrays["u"].flat[0])

    def test_generation_one_checkpoint_path_is_bound_to_source_inventory(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            _members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            captured = self._capture(root, store, pins)
            moved = replace(
                captured.generation1_checkpoint,
                path=(
                    origin.SOURCE_STORE_RELATIVE
                    + "/checkpoints/foreign-"
                    + captured.generation1_bridge_content_id
                    + ".json"
                ),
            )
            identity = origin._identity_record(
                selected_method=captured.selected_method,
                origin_event=captured.origin_event,
                origin_time=captured.origin_time,
                generation1_bridge_content_id=captured.generation1_bridge_content_id,
                generation1_checkpoint=moved,
                pref27_compact=captured.pref27_compact,
                pref27_tree_sha256=captured.pref27_tree_sha256,
                original_leaves=captured.original_leaves,
                source_inventory=captured.source_inventory,
                members=captured.members,
            )
            with patch.object(origin, "_OWNER_PINS", pins), self.assertRaisesRegex(
                Pro20OriginError, "source-inventory leaf"
            ):
                replace(
                    captured,
                    generation1_checkpoint=moved,
                    captured_bytes=canonical_json_bytes(identity),
                )

    def test_source_drift_across_capture_is_rejected(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            root = Path(tmp).resolve()
            _members, store, pins, _inventory, _receipt, _checkpoint = _install_world(root)
            real_inventory = origin.inventory_source_tree
            calls = {"n": 0}

            def drifting_inventory(path, relative=origin.SOURCE_TREE_RELATIVE):
                calls["n"] += 1
                inventory = real_inventory(path, relative)
                if calls["n"] == 1:
                    journal = root / (
                        f"{origin.SOURCE_STORE_RELATIVE}/journal/"
                        f"00000000000000000000-{_digest('journal')}.journal"
                    )
                    journal.write_bytes(b"changed-during-capture")
                return inventory

            with (
                patch.object(origin, "_OWNER_PINS", pins),
                patch.object(origin, "_open_historical_store", return_value=store),
                patch.object(origin, "inventory_source_tree", drifting_inventory),
            ):
                with self.assertRaisesRegex(Pro20OriginError, "drifted during capture"):
                    capture_pro20_historical_origin(root)

    def test_symlink_root_is_refused(self) -> None:
        with TemporaryDirectory(prefix="fgc-pro20-origin-") as tmp:
            base = Path(tmp).resolve()
            real = base / "repo"
            real.mkdir()
            members, store, pins, _inventory, _receipt, _checkpoint = _install_world(real)
            link = base / "link"
            os.symlink(real, link)
            with self.assertRaisesRegex(Pro20OriginError, "symlink"):
                with (
                    patch.object(origin, "_OWNER_PINS", pins),
                    patch.object(origin, "_open_historical_store", return_value=store),
                ):
                    capture_pro20_historical_origin(link)
            captured = self._capture(real, store, pins)
            self.assertEqual(tuple(item.member_key for item in captured.members), MEMBER_KEYS)
            del members


if __name__ == "__main__":
    unittest.main()
