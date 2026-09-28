"""Focused tests for the exact RCV3 recovery-successor authority."""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution import hlt16_campaign_recovery as recovery
from recursive_horizons.fgc.evolution import tdg8_rcv3_rec1_authority as authority
from recursive_horizons.fgc.evolution.hlt16_campaign_store import HLT16CampaignStore


ROOT = Path(__file__).resolve().parents[1]
LIVE_STORE = ROOT / authority.DESTINATION_PATH


def _receipt() -> authority.TDG8RCV3REC1Authority:
    return authority.TDG8RCV3REC1Authority(
        authority_commit="a" * 40,
        prior_auth1_commit=authority.PRIOR_AUTH1_COMMIT,
        normalization_fix_commit=authority.NORMALIZATION_FIX_COMMIT,
        projection_id=authority.PROJECTION_ID,
        original_execution_commit=authority.ORIGINAL_EXECUTION_COMMIT,
        original_plan_sha256=authority.ORIGINAL_PLAN_SHA256,
        campaign_id=authority.CAMPAIGN_ID,
        destination_path=authority.DESTINATION_PATH,
        entry_checkpoint_sha256=authority.ENTRY_CHECKPOINT_SHA256,
        predicted_checkpoint_sha256=authority.PREDICTED_CHECKPOINT_SHA256,
        config_sha256="b" * 64,
        result_sha256="c" * 64,
        implementation_inventory=(),
        environment={},
    )


def _compact(config_raw: bytes) -> bytes:
    config = authority._parse_config(config_raw)
    inventory = [
        {"role": role, "path": path, "sha256": f"{index + 1:064x}"}
        for index, (role, path) in enumerate(authority.IMPLEMENTATION_INVENTORY)
    ]
    result = {
        "schema_version": authority.SCHEMA_VERSION,
        "artifact_id": authority.ARTIFACT_ID,
        "project_version": authority.PROJECT_VERSION,
        "target_protocol": authority.TARGET_PROTOCOL,
        "classification": authority.CLASSIFICATION,
        "gate_status": "pass",
        "source_config_sha256": authority._sha(config_raw),
        "artifact_payload": authority._payload(
            config,
            inventory,
            {"implementation": "CPython", "python_version": "3.13.0",
             "numpy_version": "2.3.0", "system": "Darwin", "machine": "arm64",
             "platform": "darwin", "byteorder": "little"},
        ),
    }
    return authority.canonical(result)


class TDG8RCV3REC1AuthorityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config_raw = (ROOT / authority.CONFIG_PATH).read_bytes()
        self.result_raw = _compact(self.config_raw)

    def _copy_store(self, directory: str) -> Path:
        destination = Path(directory) / "calibration"
        shutil.copytree(LIVE_STORE, destination)
        return destination

    def test_compact_contract_is_store_blind_and_premise_only(self) -> None:
        with patch.object(
            authority, "HLT16CampaignStore",
            side_effect=AssertionError("compact verifier opened a store"),
        ):
            result = authority.validate_compact(self.config_raw, self.result_raw)
        payload = result["artifact_payload"]
        self.assertTrue(payload["scope"]["exact_predicted_generation10_required_before_PDE"])
        self.assertTrue(payload["claims"]["exact_recovery_reconciliation_authorized"])
        self.assertFalse(payload["claims"]["recovery_completed"])
        self.assertFalse(payload["claims"]["event24_completed"])
        self.assertFalse(payload["claims"]["candidate_execution_authorized"])

    def test_compact_one_bit_or_promotion_mutations_fail(self) -> None:
        original = json.loads(self.result_raw)
        mutations = []
        checkpoint = deepcopy(original)
        checkpoint["artifact_payload"]["predicted_recovery"]["checkpoint_sha256"] = "0" * 64
        mutations.append(checkpoint)
        suffix = deepcopy(original)
        suffix["artifact_payload"]["recovery_entry"]["seq12_sha256"] = "0" * 64
        mutations.append(suffix)
        retry = deepcopy(original)
        retry["artifact_payload"]["bounded_policy"]["temporal_retry_cap_per_macro_step"] = 33
        mutations.append(retry)
        promotion = deepcopy(original)
        promotion["artifact_payload"]["claims"]["event24_completed"] = True
        mutations.append(promotion)
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                with self.assertRaises(authority.TDG8RCV3REC1AuthorityError):
                    authority.validate_compact(self.config_raw, authority.canonical(mutation))

    def test_exact_entry_and_two_phase_recovery_use_only_a_copy(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            destination = self._copy_store(temporary)
            before = authority.inspect_store(destination, _receipt())
            self.assertEqual(before.phase, authority.REC1Phase.EXACT_PRE_RECOVERY)
            store = HLT16CampaignStore(destination)
            checkpoint = store.authenticated_snapshot().checkpoint
            token = store.take_over_stale_writer(
                checkpoint, host=authority.os.uname().nodename,
                owner_token="1" * 64,
            )
            pending = authority.inspect_store(
                destination, _receipt(), owned_writer_pid=authority.os.getpid(),
            )
            self.assertEqual(pending.phase, authority.REC1Phase.GEN9_RECONCILE_PENDING)
            decision = recovery.reconcile(store)
            self.assertEqual(decision.checkpoint_sha256, authority.PREDICTED_CHECKPOINT_SHA256)
            handoff = authority.inspect_store(
                destination, _receipt(), owned_writer_pid=authority.os.getpid(),
            )
            self.assertEqual(handoff.phase, authority.REC1Phase.GEN10_HANDOFF_PENDING)
            self.assertTrue(handoff.physical_state_unchanged_at_recovery)
            store.release_active_writer(
                owner_token=token, host=authority.os.uname().nodename,
                pid=authority.os.getpid(),
            )
            clean = authority.inspect_store(destination, _receipt())
            self.assertEqual(clean.phase, authority.REC1Phase.EXACT_CLEAN_GEN10)
            self.assertFalse(clean.active_write)

    def test_repeated_takeover_crash_cuts_remain_bounded(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            destination = self._copy_store(temporary)
            store = HLT16CampaignStore(destination)
            checkpoint = store.authenticated_snapshot().checkpoint
            store.take_over_stale_writer(
                checkpoint, host=authority.os.uname().nodename,
                pid=91_000_001, owner_token="2" * 64,
            )
            store.take_over_stale_writer(
                checkpoint, host=authority.os.uname().nodename,
                pid=91_000_002, owner_token="3" * 64,
            )
            boundary = authority.inspect_store(destination, _receipt())
            self.assertEqual(boundary.phase, authority.REC1Phase.GEN9_RECONCILE_PENDING)
            self.assertGreaterEqual(boundary.retired_writer_count, 2)

    def test_foreign_leaf_or_one_bit_suffix_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            destination = self._copy_store(temporary)
            (destination / "states" / ("f" * 64 + ".json")).write_bytes(b"{}")
            with self.assertRaises(authority.TDG8RCV3REC1AuthorityError):
                authority.inspect_store(destination, _receipt())
        with tempfile.TemporaryDirectory() as temporary:
            destination = self._copy_store(temporary)
            path = destination / authority.SEQ12_PATH
            raw = bytearray(path.read_bytes())
            raw[-2] ^= 1
            path.write_bytes(raw)
            with self.assertRaises(authority.TDG8RCV3REC1AuthorityError):
                authority.inspect_store(destination, _receipt())

    def test_owner_token_is_absent_from_tracked_contract(self) -> None:
        config = authority._parse_config(self.config_raw)
        self.assertFalse(config["bounded_policy"]["owner_token_tracked"])
        tracked = self.config_raw.decode("utf-8")
        self.assertNotIn("owner_token =", tracked)
        live = (LIVE_STORE / authority.ENTRY_WRITER_PATH).read_text("utf-8")
        token = json.loads(live)["owner_token"]
        self.assertNotIn(token, tracked)


if __name__ == "__main__":
    unittest.main()
