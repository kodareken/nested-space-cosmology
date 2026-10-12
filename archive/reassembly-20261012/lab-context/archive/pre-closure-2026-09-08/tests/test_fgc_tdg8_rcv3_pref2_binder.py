"""Adversarial tests for the outcome-neutral RCV3 terminal binder."""

from __future__ import annotations

import ast
from copy import deepcopy
import json
from pathlib import Path
import unittest
from unittest import mock

from recursive_horizons.fgc.evolution import tdg8_rcv3_pref1_binder as pref1
from recursive_horizons.fgc.evolution import tdg8_rcv3_pref2_binder as binder


ROOT = Path(__file__).resolve().parents[1]


class TDG8RCV3PREF2BinderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = (ROOT / binder.CONFIG_PATH).read_bytes()
        cls.result = (ROOT / binder.RESULT_PATH).read_bytes()
        cls.snapshot = pref1._snapshot(ROOT, binder.STORE_PATH)

    def test_compact_result_and_nonpromotion_contract_pass(self) -> None:
        result = binder.validate_compact_result(self.config, self.result)
        payload = result["artifact_payload"]
        self.assertEqual(result["artifact_id"], binder.ARTIFACT_ID)
        self.assertEqual(
            payload["terminal_evidence"]["terminal"]["reason"],
            "temporal_retry_exhausted",
        )
        self.assertIs(
            payload["claims"]["accepted_physical_state_preserved_through_terminal"],
            True,
        )
        self.assertIs(payload["claims"]["common_event_completed"], False)
        self.assertIs(payload["claims"]["GR0_calibration_completed"], False)
        self.assertIs(payload["claims"]["candidate_execution_authorized"], False)
        self.assertIs(payload["claims"]["mechanism_result_earned"], False)
        self.assertIs(payload["scope"]["successor_remedy_selected"], False)

    def test_live_terminal_passes_complete_independent_binding(self) -> None:
        observed = binder.bind_terminal_store(self.config, ROOT)
        self.assertEqual(observed, binder.expected_evidence(self.config))
        self.assertEqual(
            observed["recovery"]["checkpoint_sha256"],
            binder.RECOVERY["checkpoint_sha256"],
        )
        self.assertEqual(observed["terminal"]["retry_count"], 22)

    def test_compact_validation_is_store_blind(self) -> None:
        with mock.patch.object(
            pref1, "_snapshot", side_effect=AssertionError("store was read")
        ):
            binder.validate_compact_result(self.config, self.result)

    def test_binder_import_graph_has_no_recovery_runner_or_evolution_runtime(self) -> None:
        source = (
            ROOT
            / "src/recursive_horizons/fgc/evolution/tdg8_rcv3_pref2_binder.py"
        ).read_text(encoding="utf-8")
        tree = ast.parse(source)
        imports: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imports.append(node.module or "")
        forbidden = (
            "tdg8_rcv3_rec1_authority", "hlt16_campaign_recovery",
            "hlt16_campaign_runtime", "hlt16_campaign_store",
            "hlt16_progression_attempt", "tdg8_successor_runtime",
        )
        self.assertFalse(any(any(name in module for name in forbidden) for module in imports))
        self.assertNotIn("advance_member", source)
        self.assertNotIn("solve_accelerations", source)

    def _journal(self, sequence: int) -> dict[str, object]:
        prefix = f"journal/{sequence:020d}-"
        path = next(path for path in self.snapshot.leaves if path.startswith(prefix))
        return json.loads(self.snapshot.leaves[path].raw)

    def _member(self, generation: int) -> dict[str, object]:
        prefix = f"checkpoints/{generation:020d}-"
        path = next(path for path in self.snapshot.leaves if path.startswith(prefix))
        return json.loads(self.snapshot.leaves[path].raw)["members"][binder.MEMBER_KEY]

    def test_accepted_state_accounting_mutations_are_rejected_directly(self) -> None:
        anchor = self._member(9)
        mutations = (
            ("debit", lambda member: member["ledger"]["accumulated_debit_vector_hex"].__setitem__(0, "0x0.0p+0")),
            ("accepted_count", lambda member: member["ledger"].__setitem__("accepted_macro_step_count", 6)),
            ("accepted_time", lambda member: member["ledger"].__setitem__("last_accepted_time_hex", "0x1.8p+0")),
            ("last_accepted_retry", lambda member: member["ledger"].__setitem__("last_accepted_macro_step_temporal_retry_count", 1)),
            ("boundary_time", lambda member: member["cursor"]["accepted_boundary_time"].__setitem__("binary64_hex", "0x1.8p+0")),
            ("step_index", lambda member: member["cursor"].__setitem__("previous_step_index", 363)),
        )
        for label, mutate in mutations:
            with self.subTest(label=label):
                member = deepcopy(self._member(10))
                mutate(member)
                with self.assertRaises(binder.TDG8RCV3PREF2Error) as caught:
                    binder._validate_accepted_state_preservation(member, anchor, 10)
                self.assertEqual(caught.exception.stop_id, "PREF2_STATE_DRIFT")

    def test_terminal_evidence_mutations_are_rejected(self) -> None:
        for field, replacement in (
            ("retry_count_for_current_macro_step", 21),
            ("classification_is_numerical_not_physical", False),
            ("accepted_state_bitwise_preserved", False),
            ("initial_state_sha256", "0" * 64),
        ):
            with self.subTest(field=field):
                record = deepcopy(self._journal(49))
                record["payload"]["evidence"][field] = replacement
                with self.assertRaises(binder.TDG8RCV3PREF2Error) as caught:
                    binder._rejection_evidence(record, 22)
                self.assertEqual(caught.exception.stop_id, "PREF2_LINEAGE_DRIFT")

    def test_channel_inventory_and_failed_set_are_bound(self) -> None:
        record = deepcopy(self._journal(49))
        record["payload"]["evidence"]["channel_admissions"].reverse()
        with self.assertRaises(binder.TDG8RCV3PREF2Error) as caught:
            binder._rejection_evidence(record, 22)
        self.assertEqual(caught.exception.stop_id, "PREF2_LINEAGE_DRIFT")

        record = deepcopy(self._journal(49))
        record["payload"]["evidence"]["failed_channels"] = []
        with self.assertRaises(binder.TDG8RCV3PREF2Error) as caught:
            binder._rejection_evidence(record, 22)
        self.assertEqual(caught.exception.stop_id, "PREF2_LINEAGE_DRIFT")

    def test_terminal_lock_inventory_rejects_foreign_or_active_writer(self) -> None:
        leaves = dict(self.snapshot.leaves)
        leaves["locks/active-write.lock"] = next(
            leaf for path, leaf in leaves.items()
            if path.startswith("locks/.active-write.lock.hlt16-quarantine-")
        )
        altered = pref1._Snapshot(leaves, self.snapshot.directories)
        with self.assertRaises(binder.TDG8RCV3PREF2Error) as caught:
            binder._validate_tree_grammar(altered)
        self.assertEqual(caught.exception.stop_id, "PREF2_TREE_DRIFT")

    def test_claim_promotion_and_config_threshold_mutation_are_rejected(self) -> None:
        promoted = self.config.replace(
            b"mechanism_result_earned = false",
            b"mechanism_result_earned = true",
        )
        with self.assertRaises(binder.TDG8RCV3PREF2Error) as caught:
            binder.validate_compact_result(promoted, self.result)
        self.assertEqual(caught.exception.stop_id, "PREF2_CONFIG_DRIFT")

        changed = self.config.replace(
            b"retry_count = 22", b"retry_count = 21", 1,
        )
        with self.assertRaises(binder.TDG8RCV3PREF2Error) as caught:
            binder.validate_compact_result(changed, self.result)
        self.assertEqual(caught.exception.stop_id, "PREF2_CONFIG_DRIFT")

    def test_result_terminal_or_successor_promotion_mutation_is_rejected(self) -> None:
        value = json.loads(self.result)
        value["artifact_payload"]["terminal_evidence"]["terminal"]["reason"] = "physics"
        with self.assertRaises(binder.TDG8RCV3PREF2Error) as caught:
            binder.validate_compact_result(self.config, binder.canonical_result(value))
        self.assertEqual(caught.exception.stop_id, "PREF2_COMPACT_DRIFT")

        value = json.loads(self.result)
        value["artifact_payload"]["scope"]["successor_remedy_selected"] = True
        with self.assertRaises(binder.TDG8RCV3PREF2Error) as caught:
            binder.validate_compact_result(self.config, binder.canonical_result(value))
        self.assertEqual(caught.exception.stop_id, "PREF2_COMPACT_DRIFT")

    def test_snapshot_failure_is_fail_closed_and_never_reclassified_as_physics(self) -> None:
        error = pref1.TDG8RCV3PREF1Error("PREF1_PATH_UNSAFE", "synthetic")
        with (
            mock.patch.object(binder, "_bind_predecessor", return_value={}),
            mock.patch.object(pref1, "_snapshot", side_effect=error),
            self.assertRaises(binder.TDG8RCV3PREF2Error) as caught,
        ):
            binder.bind_terminal_store(self.config, ROOT)
        self.assertEqual(caught.exception.stop_id, "PREF2_PATH_UNSAFE")
        self.assertNotIn("physical", caught.exception.detail.lower())


if __name__ == "__main__":
    unittest.main()
