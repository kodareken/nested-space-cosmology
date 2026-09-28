from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
from io import BytesIO
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.boundary_domain import CausalBudgetState
from recursive_horizons.fgc.evolution.hlt16_member_codec import (
    ARRAY_NAMES,
    HLT16MemberSnapshot,
    validate_codec_metadata,
)
from recursive_horizons.fgc.evolution.hlt16_state_store import (
    HLT16StateStore,
    HLT16StateStoreError,
    authorize_gen0_import,
    decode_payload,
    digest,
    encode_payload,
)
from recursive_horizons.fgc.evolution.proto15_runtime import _hash, _ledger_mapping
from recursive_horizons.fgc.evolution.proto17_pure_construction import cursor
from recursive_horizons.fgc.evolution.proto5_runtime import GR0RuntimeMonitorState
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_runtime import (
    TDG6TemporalLedger,
    tdg6_checkpoint_extension,
)


def arrays(points: int = 2049, labels: int = 3, events: int = 2) -> dict[str, np.ndarray]:
    return {
        "u": np.zeros((points, 6), dtype="<f8"),
        "p": np.ones((points, 6), dtype="<f8"),
        "q": np.full((points, 6), 2.0, dtype="<f8"),
        "grid_coordinates": np.linspace(0.0, 128.0, points, dtype="<f8"),
        "tracer_labels": np.linspace(1.0, 3.0, labels, dtype="<f8"),
        "tracer_positions": np.linspace(1.0, 3.0, labels, dtype="<f8"),
        "tracer_proper_times": np.zeros(labels, dtype="<f8"),
        "event_proper_times": np.zeros((events, labels), dtype="<f8"),
        "event_fields": np.zeros((events, labels, 6), dtype="<f8"),
    }


def snapshot() -> HLT16MemberSnapshot:
    ledger = TDG6TemporalLedger.zero(initial_time=23 / 16)
    extension, debit = tdg6_checkpoint_extension(ledger)
    source_cursor = cursor(
        {
            "protocol_artifact_id": "FGC-2-SF1-PROTO17",
            "campaign_id": "HLT16-TEST",
            "member_key": "RK4-2049",
            "method": "RK4",
            "point_count": 2049,
            "committed_common_event_index": 23,
            "accepted_boundary_time": {
                "rational": "23/16",
                "binary64_hex": (23 / 16).hex(),
            },
            "accepted_state_sha256": "a" * 64,
            "previous_step_index": 0,
            "previous_transaction_serial": 0,
            "TDG6_ledger_sha256": _hash(_ledger_mapping(ledger)),
            "cursor_generation": 0,
            "event_target_time": {
                "rational": "3/2",
                "binary64_hex": (3 / 2).hex(),
            },
            "attempt_serial": 0,
            "attempt_id": "HLT16-TEST:RK4-2049:23:0",
            "mode": "FRESH_READY",
            "retry_successor_payload_or_none": None,
            "journal_tip_sha256": "0" * 64,
            "cursor_chain_parent_sha256": "0" * 64,
        }
    )
    template = {"amplitude": "3", "method": "RK4", "point_count": 2049}
    metadata = validate_codec_metadata(
        {
            "runtime_identity": {
                "protocol_artifact_id": "FGC-2-SF1-PROTO18",
                "campaign_id": "HLT16-TEST",
                "branch": "GR-0",
                "amplitude": "3",
                "member_key": "RK4-2049",
                "method": "RK4",
                "point_count": 2049,
                "accepted_time": {
                    "rational": "23/16",
                    "binary64_hex": (23 / 16).hex(),
                },
                "generation": 0,
            },
            "template": template,
            "template_sha256": digest(template),
            "accepted_plan": None,
            "counters": {
                "step_index": 0,
                "transaction_serial": 0,
                "accepted_macro_steps": 0,
                "source_retry_count": 0,
                "CFL_retry_count": 0,
            },
            "tdg6": {
                "extension": extension,
                "accumulated_debit_hex": [float(item).hex() for item in debit],
            },
            "monitor": asdict(GR0RuntimeMonitorState()),
            "causal": asdict(CausalBudgetState(accepted_time=23 / 16)),
            "tracer": {
                "label_count": 3,
                "event_count": 2,
                "cutoff_hex": float(16).hex(),
                "outer_radius_hex": float(128).hex(),
            },
            "provenance_cursor": {
                "source_protocol": "FGC-2-SF1-PROTO17",
                "bridge": "authenticated_GEN0_predecessor",
                "payload": source_cursor,
            },
        }
    )
    return HLT16MemberSnapshot(arrays(), metadata)


class HLT16StateStoreTests(unittest.TestCase):
    def test_payload_round_trip_is_semantic_and_bounded(self) -> None:
        values = arrays(points=17)
        raw, semantic, manifest = encode_payload(values)
        decoded, observed, rebuilt = decode_payload(
            raw, expected_semantic_sha256=semantic, expected_manifest=manifest
        )
        self.assertEqual(observed, semantic)
        self.assertEqual(rebuilt, manifest)
        self.assertEqual(
            {name: value.tobytes() for name, value in decoded.items()},
            {name: value.tobytes() for name, value in values.items()},
        )
        with patch(
            "recursive_horizons.fgc.evolution.hlt16_state_store.MAX_ARRAY_MEMBER_BYTES",
            8,
        ):
            with self.assertRaises(HLT16StateStoreError):
                decode_payload(raw)

    def test_wrong_inventory_dtype_and_nonfinite_are_rejected(self) -> None:
        value = arrays(points=17)
        value["u"] = value["u"].astype(">f8")
        with self.assertRaises(HLT16StateStoreError):
            encode_payload(value)
        value = arrays(points=17)
        value["u"][0, 0] = np.nan
        with self.assertRaises(HLT16StateStoreError):
            encode_payload(value)
        value = arrays(points=17)
        value.pop("q")
        with self.assertRaises(HLT16StateStoreError):
            encode_payload(value)

    def test_gen0_authority_persist_load_and_idempotent_reuse(self) -> None:
        item = snapshot()
        authority = authorize_gen0_import(item)
        with TemporaryDirectory() as directory:
            store = HLT16StateStore(Path(directory))
            with self.assertRaises(HLT16StateStoreError):
                store.persist_snapshot(item)
            first = store.persist_snapshot(item, gen0_authority=authority)
            second = store.persist_snapshot(item, gen0_authority=authority)
            self.assertEqual(first.descriptor_sha256, second.descriptor_sha256)
            loaded = store.load_state(first.descriptor_sha256)
            self.assertEqual(loaded.raw_archive_sha256, first.raw_archive_sha256)
            self.assertEqual(
                {name: value.tobytes() for name, value in loaded.arrays.items()},
                {name: value.tobytes() for name, value in item.arrays.items()},
            )

    def test_payload_descriptor_and_authority_mutations_fail_closed(self) -> None:
        item = snapshot()
        authority = authorize_gen0_import(item)
        with TemporaryDirectory() as directory:
            root = Path(directory)
            store = HLT16StateStore(root)
            persisted = store.persist_snapshot(item, gen0_authority=authority)
            payload = root / "payloads" / f"{persisted.semantic_sha256}.npz"
            payload.write_bytes(payload.read_bytes() + b"x")
            with self.assertRaises(HLT16StateStoreError):
                store.load_state(persisted.descriptor_sha256)
        changed = type(authority)(
            **{**authority.body(), "source_state_object_sha256": "b" * 64},
            authority_sha256=authority.authority_sha256,
        )
        with TemporaryDirectory() as directory:
            with self.assertRaises(HLT16StateStoreError):
                HLT16StateStore(Path(directory)).persist_snapshot(
                    item, gen0_authority=changed
                )

    def test_symlink_parent_and_publication_fault_are_not_adopted(self) -> None:
        item = snapshot()
        authority = authorize_gen0_import(item)
        with TemporaryDirectory() as directory:
            root = Path(directory)
            outside = root / "outside"
            outside.mkdir()
            os.symlink(outside, root / "payloads")
            with self.assertRaises(HLT16StateStoreError):
                HLT16StateStore(root).persist_snapshot(item, gen0_authority=authority)
        for phase in (
            "before_write",
            "after_write",
            "after_file_fsync",
            "after_link",
            "after_parent_fsync",
            "after_stage_cleanup",
        ):
            with self.subTest(phase=phase), TemporaryDirectory() as directory:
                def fault(observed: str) -> None:
                    if observed == phase:
                        raise RuntimeError(phase)

                root = Path(directory)
                with self.assertRaises(RuntimeError):
                    HLT16StateStore(root, fault_hook=fault).persist_snapshot(
                        item, gen0_authority=authority
                    )
                # The only permitted residue is the deterministic task-owned
                # publication stage, possibly hard-linked to its final inode.
                self.assertEqual(list(root.rglob(".*.tmp-*")), [])
                for residue in root.rglob(".*"):
                    if residue.is_file():
                        self.assertIn(".hlt16-stage-", residue.name)
                recovered = HLT16StateStore(root).persist_snapshot(
                    item, gen0_authority=authority
                )
                loaded = HLT16StateStore(root).load_state(
                    recovered.descriptor_sha256
                )
                self.assertEqual(
                    {name: value.tobytes() for name, value in loaded.arrays.items()},
                    {name: value.tobytes() for name, value in item.arrays.items()},
                )
                self.assertEqual(
                    [
                        path
                        for path in root.rglob(".*")
                        if path.is_file() and ".hlt16-stage-" in path.name
                    ],
                    [],
                )


if __name__ == "__main__":
    unittest.main()
