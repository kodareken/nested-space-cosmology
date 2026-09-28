"""Focused tests for the prospective RCV3 bounded-event authority."""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from recursive_horizons.fgc.evolution import tdg8_rcv3_execution_authority as authority


ROOT = Path(__file__).resolve().parents[1]


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
            {
                "implementation": "CPython", "python_version": "3.13.0",
                "numpy_version": "2.3.0", "system": "Darwin",
                "machine": "arm64", "platform": "darwin", "byteorder": "little",
            },
        ),
    }
    return authority.canonical(result)


def _receipt() -> authority.TDG8RCV3ExecutionAuthority:
    return authority.TDG8RCV3ExecutionAuthority(
        authority_commit="a" * 40,
        pref1_authority_commit=authority.PREF1_AUTHORITY_COMMIT,
        pref1_result_sha256=authority.PREF1_RESULT_SHA256,
        projection_id=authority.PROJECTION_ID,
        receipt_sha256=authority.RECEIPT_SHA256,
        original_execution_commit=authority.ORIGINAL_EXECUTION_COMMIT,
        original_plan_sha256=authority.ORIGINAL_PLAN_SHA256,
        config_sha256="b" * 64,
        result_sha256="c" * 64,
        campaign_id=authority.CAMPAIGN_ID,
        destination_path=authority.DESTINATION_PATH,
        entry_checkpoint_sha256=authority.ENTRY_CHECKPOINT_SHA256,
        entry_journal_sha256=authority.ENTRY_JOURNAL_SHA256,
        temporal_retry_cap=authority.TEMPORAL_RETRY_CAP,
        minimum_macro_step_hex=authority.MINIMUM_MACRO_STEP_HEX,
        implementation_inventory=(),
        environment={},
    )


class TDG8RCV3ExecutionAuthorityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config_raw = (ROOT / authority.CONFIG_PATH).read_bytes()
        self.result_raw = _compact(self.config_raw)

    def test_compact_contract_is_store_blind_and_premise_only(self) -> None:
        with patch.object(
            authority, "HLT16CampaignStore",
            side_effect=AssertionError("compact verification opened a store"),
        ):
            result = authority.validate_compact(self.config_raw, self.result_raw)
        payload = result["artifact_payload"]
        self.assertTrue(payload["scope"]["compact_verifier_store_blind"])
        self.assertTrue(payload["claims"]["bounded_retry_chain_authorized"])
        self.assertFalse(payload["claims"]["execution_started"])
        self.assertFalse(payload["claims"]["common_event_completed"])
        self.assertFalse(payload["claims"]["candidate_execution_authorized"])
        self.assertFalse(payload["claims"]["physical_result_earned"])

    def test_compact_mutations_cannot_change_boundary_policy_or_claims(self) -> None:
        original = json.loads(self.result_raw)
        mutations = []
        retry = deepcopy(original)
        retry["artifact_payload"]["bounded_policy"]["temporal_retry_cap_per_macro_step"] = 33
        mutations.append(retry)
        anchor = deepcopy(original)
        anchor["artifact_payload"]["internal_boundary"]["checkpoint_generation"] = 10
        mutations.append(anchor)
        receipt = deepcopy(original)
        receipt["artifact_payload"]["projection"]["receipt_sha256"] = "0" * 64
        mutations.append(receipt)
        promoted = deepcopy(original)
        promoted["artifact_payload"]["claims"]["physical_result_earned"] = True
        mutations.append(promoted)
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                with self.assertRaises(authority.TDG8RCV3ExecutionAuthorityError):
                    authority.validate_compact(self.config_raw, authority.canonical(mutation))

    def test_frozen_scope_preserves_original_identity_and_excludes_candidates(self) -> None:
        config = authority._parse_config(self.config_raw)
        self.assertEqual(config["predecessor"]["pref1_authority_commit"], authority.PREF1_AUTHORITY_COMMIT)
        self.assertEqual(config["internal_boundary"]["authorization_commit"], authority.ORIGINAL_EXECUTION_COMMIT)
        self.assertEqual(config["internal_boundary"]["plan_sha256"], authority.ORIGINAL_PLAN_SHA256)
        self.assertEqual(config["internal_boundary"]["checkpoint_generation"], 9)
        self.assertEqual(config["internal_boundary"]["retry_depth"], 2)
        self.assertEqual(config["bounded_policy"]["temporal_retry_cap_per_macro_step"], 32)
        self.assertEqual(config["bounded_policy"]["event_successor"], 24)
        self.assertTrue(config["scope"]["candidate_branches_forbidden"])
        self.assertFalse(config["scope"]["generation_zero_reimport"])

    def test_live_projection_is_a_safe_generation_nine_retry_boundary(self) -> None:
        boundary = authority.inspect_descendant(ROOT, _receipt())
        self.assertEqual(boundary.checkpoint_generation, 9)
        self.assertEqual(boundary.checkpoint_sha256, authority.ENTRY_CHECKPOINT_SHA256)
        self.assertEqual(boundary.journal_sequence, 10)
        self.assertEqual(boundary.event, 23)
        self.assertEqual(boundary.target, authority.TARGET)
        self.assertTrue(boundary.external_receipt_authenticated)
        self.assertTrue(boundary.entry_ancestor_present)
        self.assertEqual(boundary.temporal_pending_members, ("RK4-2049",))
        self.assertTrue(boundary.safe_to_restart)
        self.assertFalse(boundary.terminal)

    def test_receipt_raw_or_semantic_mutation_is_rejected(self) -> None:
        real = (ROOT / authority.RECEIPT_PATH).read_bytes()
        with patch.object(authority.rcv2.rcv1, "_read_leaf", return_value=real[:-1] + b" "):
            with self.assertRaisesRegex(
                authority.TDG8RCV3ExecutionAuthorityError, "receipt raw identity",
            ):
                authority._authenticate_external_receipt(ROOT)

    def test_read_only_inspection_rejects_a_moving_descendant(self) -> None:
        first = SimpleNamespace(
            checkpoint=SimpleNamespace(sha256="1" * 64),
            suffix_kinds=(), suffix_records=(), staging_paths=(),
            orphan_descriptor_sha256=None,
            orphan_payload_semantic_sha256=None,
            terminal_lock_present=False,
        )
        second = SimpleNamespace(
            checkpoint=SimpleNamespace(sha256="2" * 64),
            suffix_kinds=(), suffix_records=(), staging_paths=(),
            orphan_descriptor_sha256=None,
            orphan_payload_semantic_sha256=None,
            terminal_lock_present=False,
        )
        store = SimpleNamespace(
            authenticated_snapshot=Mock(side_effect=(first, second)),
            inspect_recovery=lambda: SimpleNamespace(),
        )
        with (
            patch.object(authority, "_real_directory", return_value=ROOT),
            patch.object(authority, "_authenticate_external_receipt", return_value=None),
            patch.object(authority, "HLT16CampaignStore", return_value=store),
        ):
            with self.assertRaisesRegex(
                authority.TDG8RCV3ExecutionAuthorityError,
                "changed during read-only inspection",
            ):
                authority.inspect_descendant(ROOT, _receipt())

    def test_read_only_inspection_rejects_writer_lock_movement(self) -> None:
        snapshot = SimpleNamespace(
            checkpoint=SimpleNamespace(sha256="1" * 64),
            suffix_kinds=(), suffix_records=(), staging_paths=(),
            orphan_descriptor_sha256=None,
            orphan_payload_semantic_sha256=None,
            terminal_lock_present=False,
        )
        clean = SimpleNamespace(
            active_write=False, writer_state=None, state="clean_checkpoint",
            safe_to_restart=True,
        )
        live = SimpleNamespace(
            active_write=True, writer_state="live_local_writer",
            state="live_local_writer", safe_to_restart=False,
        )
        store = SimpleNamespace(
            authenticated_snapshot=Mock(side_effect=(snapshot, snapshot)),
            inspect_recovery=Mock(side_effect=(clean, live)),
        )
        with (
            patch.object(authority, "_real_directory", return_value=ROOT),
            patch.object(authority, "_authenticate_external_receipt", return_value=None),
            patch.object(authority, "HLT16CampaignStore", return_value=store),
        ):
            with self.assertRaisesRegex(
                authority.TDG8RCV3ExecutionAuthorityError,
                "changed during read-only inspection",
            ):
                authority.inspect_descendant(ROOT, _receipt())

    def test_typed_receipt_cannot_switch_projection_or_plan(self) -> None:
        projection = deepcopy(_receipt())
        object.__setattr__(projection, "projection_id", "other")
        with self.assertRaises(authority.TDG8RCV3ExecutionAuthorityError):
            authority._require_typed_receipt(projection)
        plan = deepcopy(_receipt())
        object.__setattr__(plan, "original_plan_sha256", "0" * 64)
        with self.assertRaises(authority.TDG8RCV3ExecutionAuthorityError):
            authority._require_typed_receipt(plan)

    def test_authority_delta_is_exact_and_contains_no_runtime_engine_edits(self) -> None:
        expected_core = {
            authority.CONFIG_PATH, authority.RESULT_PATH, authority.OWNER_DOCUMENT,
            "scripts/reproduce_fgc_tdg8_rcv3_auth1.py",
            "scripts/run_fgc_tdg8_rcv3_event.py",
            "src/recursive_horizons/fgc/evolution/tdg8_rcv3_execution_authority.py",
            "tests/test_fgc_tdg8_rcv3_execution_authority.py",
            "tests/test_fgc_tdg8_rcv3_event_runner.py",
        }
        self.assertTrue(expected_core <= set(authority.AUTHORITY_DELTA_PATHS))
        self.assertFalse(any(
            path.endswith(("hlt16_campaign_runtime.py", "tdg8_successor_runtime.py"))
            for path in authority.AUTHORITY_DELTA_PATHS
        ))

    def test_private_predecessor_authority_helpers_are_hash_bound(self) -> None:
        inventory = set(authority.IMPLEMENTATION_INVENTORY)
        self.assertIn(
            (
                "predecessor_authority",
                "src/recursive_horizons/fgc/evolution/tdg8_bounded_retry_authority.py",
            ),
            inventory,
        )
        self.assertIn(
            (
                "predecessor_authority",
                "src/recursive_horizons/fgc/evolution/tdg8_retry_recovery_authority.py",
            ),
            inventory,
        )


if __name__ == "__main__":
    unittest.main()
