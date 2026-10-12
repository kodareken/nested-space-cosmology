"""Focused tests for the prospective bounded TDG8 retry authority."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution import tdg8_bounded_retry_authority as authority


ROOT = Path(__file__).resolve().parents[1]


def _compact(config_raw: bytes) -> bytes:
    config = authority._parse_config(config_raw)
    inventory = [
        {"role": role, "path": path, "sha256": f"{index + 1:064x}"}
        for index, (role, path) in enumerate(authority.IMPLEMENTATION_INVENTORY)
    ]
    live = {
        "entry": authority._expected_entry(),
        "first_reconciliation": authority._expected_first_reconciliation(),
        "entry_suffix_kinds": ["tdg6_rejection", "cursor_transition"],
        "entry_physical_state_advanced": False,
    }
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
            live,
            {
                "implementation": "CPython",
                "python_version": "3.13.0",
                "numpy_version": "2.3.0",
                "system": "Darwin",
                "machine": "arm64",
                "platform": "darwin",
                "byteorder": "little",
            },
        ),
    }
    return authority.canonical(result)


def _receipt() -> authority.TDG8BoundedRetryExecutionAuthority:
    return authority.TDG8BoundedRetryExecutionAuthority(
        authority_commit="a" * 40,
        predecessor_authority_commit=authority.RCV1_AUTHORITY_COMMIT,
        original_execution_commit=authority.ORIGINAL_EXECUTION_COMMIT,
        original_plan_sha256=authority.ORIGINAL_PLAN_SHA256,
        config_sha256="b" * 64,
        result_sha256="c" * 64,
        campaign_id=authority.CAMPAIGN_ID,
        destination_path=authority.DESTINATION_PATH,
        entry_checkpoint_sha256=authority.ENTRY_CHECKPOINT_SHA256,
        entry_rejection_sha256=authority.ENTRY_REJECTION_SHA256,
        entry_transition_sha256=authority.ENTRY_TRANSITION_SHA256,
        first_recovery_checkpoint_sha256=authority.FIRST_RECOVERY_CHECKPOINT_SHA256,
        temporal_retry_cap=authority.TEMPORAL_RETRY_CAP,
        minimum_macro_step_hex=authority.MINIMUM_MACRO_STEP_HEX,
        implementation_inventory=(),
        environment={},
    )


class TDG8BoundedRetryAuthorityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config_raw = (ROOT / authority.CONFIG_PATH).read_bytes()
        self.result_raw = _compact(self.config_raw)

    def test_compact_contract_is_store_blind_and_premise_only(self) -> None:
        with (
            patch.object(
                authority, "inspect_live_entry",
                side_effect=AssertionError("compact verification opened the live store"),
            ),
            patch.object(
                authority, "HLT16CampaignStore",
                side_effect=AssertionError("compact verification constructed a live store"),
            ),
        ):
            result = authority.validate_compact(self.config_raw, self.result_raw)
        payload = result["artifact_payload"]
        self.assertTrue(payload["scope"]["compact_verifier_store_blind"])
        self.assertTrue(payload["claims"]["bounded_retry_chain_authorized"])
        self.assertFalse(payload["claims"]["physical_state_advanced_by_recovery"])
        self.assertFalse(payload["claims"]["candidate_execution_authorized"])
        self.assertFalse(payload["claims"]["physical_result_earned"])

    def test_compact_mutations_cannot_change_bounds_entry_or_claims(self) -> None:
        original = json.loads(self.result_raw.decode("ascii"))
        mutations = []
        cap = deepcopy(original)
        cap["artifact_payload"]["bounded_policy"][
            "temporal_retry_cap_per_macro_step"
        ] = 33
        mutations.append(cap)
        entry = deepcopy(original)
        entry["artifact_payload"]["entry"]["checkpoint_generation"] = 9
        mutations.append(entry)
        claim = deepcopy(original)
        claim["artifact_payload"]["claims"]["physical_result_earned"] = True
        mutations.append(claim)
        first = deepcopy(original)
        first["artifact_payload"]["first_reconciliation"][
            "physical_state_advanced"
        ] = True
        mutations.append(first)
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                with self.assertRaises(authority.TDG8BoundedRetryAuthorityError):
                    authority.validate_compact(
                        self.config_raw, authority.canonical(mutation),
                    )

    def test_policy_is_generic_but_strictly_bounded_to_one_event(self) -> None:
        config = authority._parse_config(self.config_raw)
        policy = config["bounded_policy"]
        self.assertTrue(policy["authenticated_descendants_allowed"])
        self.assertTrue(policy["recognized_publication_cuts_allowed"])
        self.assertEqual(policy["temporal_retry_cap_per_macro_step"], 32)
        self.assertEqual(
            policy["minimum_macro_step_binary64_hex"],
            authority.MINIMUM_MACRO_STEP_HEX,
        )
        self.assertTrue(policy["same_event_only"])
        self.assertEqual(policy["event_successor"], 24)
        self.assertFalse(config["scope"]["generation_zero_reimport"])
        self.assertTrue(config["scope"]["candidate_branches_forbidden"])

    def test_descendant_rejects_escape_before_loading_physical_state(self) -> None:
        checkpoint = SimpleNamespace(
            generation=10,
            authorization_commit=authority.ORIGINAL_EXECUTION_COMMIT,
            plan_sha256=authority.ORIGINAL_PLAN_SHA256,
            campaign_id=authority.CAMPAIGN_ID,
            protocol=authority.TARGET_PROTOCOL,
            members={key: object() for key in authority.MEMBER_KEYS},
            event=25,
            target={"rational": "13/8", "binary64_hex": (13 / 8).hex()},
            disposition="nonterminal",
        )
        snapshot = SimpleNamespace(checkpoint=checkpoint)
        store = SimpleNamespace(
            authenticated_snapshot=lambda: snapshot,
            inspect_recovery=lambda: SimpleNamespace(),
            load_state=unittest.mock.Mock(),
        )
        with patch.object(authority, "HLT16CampaignStore", return_value=store):
            with self.assertRaisesRegex(
                authority.TDG8BoundedRetryAuthorityError,
                "escaped the one-event scope",
            ):
                authority.inspect_descendant(ROOT, _receipt())
        store.load_state.assert_not_called()

    def test_event_24_is_only_an_exact_event_complete_boundary(self) -> None:
        checkpoint = SimpleNamespace(
            event=24,
            disposition="nonterminal",
            target=dict(authority.SUCCESSOR_TARGET),
        )
        with self.assertRaisesRegex(
            authority.TDG8BoundedRetryAuthorityError,
            "event-24 boundary differs",
        ):
            authority._require_same_event(checkpoint)
        checkpoint.disposition = "event_complete"
        authority._require_same_event(checkpoint)

    def test_any_authenticated_bounded_descendant_is_not_hash_enumerated(self) -> None:
        checkpoint = SimpleNamespace(
            generation=17,
            sha256="d" * 64,
            journal_sequence=31,
            journal_tip_sha256="e" * 64,
            authorization_commit=authority.ORIGINAL_EXECUTION_COMMIT,
            plan_sha256=authority.ORIGINAL_PLAN_SHA256,
            campaign_id=authority.CAMPAIGN_ID,
            protocol=authority.TARGET_PROTOCOL,
            members={key: object() for key in authority.MEMBER_KEYS},
            event=23,
            target=dict(authority.TARGET),
            disposition="nonterminal",
        )
        snapshot = SimpleNamespace(
            checkpoint=checkpoint,
            suffix_kinds=("tdg6_rejection",),
            suffix_records=(SimpleNamespace(),),
            suffix_classification="tdg6_rejection",
            terminal_lock_present=False,
        )
        status = SimpleNamespace(
            active_write=True,
            writer_state="stale_verified_lock",
            terminal=False,
            safe_to_restart=False,
        )
        store = SimpleNamespace(
            authenticated_snapshot=lambda: snapshot,
            inspect_recovery=lambda: status,
            load_state=lambda _descriptor: SimpleNamespace(
                arrays={"u": object(), "p": object(), "q": object()},
            ),
        )
        with (
            patch.object(authority, "HLT16CampaignStore", return_value=store),
            patch.object(authority, "_require_entry_ancestors", return_value=True),
            patch.object(
                authority, "_require_temporal_pending",
                return_value=("SSPRK3-16385",),
            ),
            patch.object(
                authority, "array_content_sha256",
                return_value=authority.ENTRY_PHYSICAL_STATE_SHA256,
            ),
        ):
            boundary = authority.inspect_descendant(ROOT, _receipt())
        self.assertEqual(boundary.checkpoint_generation, 17)
        self.assertTrue(boundary.first_recovery_ancestor_present)
        self.assertTrue(boundary.recoverable_under_rcv2)
        self.assertEqual(boundary.temporal_pending_members, ("SSPRK3-16385",))

    def test_authority_delta_includes_only_the_recovery_fix_and_rcv2_surface(self) -> None:
        expected = {
            "Makefile",
            authority.CONFIG_PATH,
            authority.RESULT_PATH,
            authority.OWNER_DOCUMENT,
            "scripts/check_repo.py",
            "scripts/reproduce_fgc_tdg8_rcv2_auth1.py",
            "scripts/run_fgc_tdg8_bounded_retry.py",
            "src/recursive_horizons/fgc/evolution/hlt16_campaign_recovery.py",
            "src/recursive_horizons/fgc/evolution/tdg8_bounded_retry_authority.py",
            "tests/test_fgc_hlt16_campaign_recovery.py",
            "tests/test_fgc_tdg8_bounded_retry_authority.py",
            "tests/test_fgc_tdg8_bounded_retry_runner.py",
        }
        self.assertEqual(set(authority.AUTHORITY_DELTA_PATHS), expected)


if __name__ == "__main__":
    unittest.main()
