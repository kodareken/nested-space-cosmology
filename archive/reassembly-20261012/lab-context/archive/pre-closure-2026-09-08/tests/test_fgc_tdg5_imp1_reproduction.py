from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_tdg5_imp1 as imp1  # noqa: E402


class FGCTDG5IMP1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = imp1.verify_canonical()

    def test_canonical_implementation_record_reproduces(self) -> None:
        self.assertEqual(self.record["artifact_id"], "FGC-1-TDG5-IMP1")
        self.assertEqual(
            self.record["classification"],
            "implemented_and_synthetically_validated_TDG5_atomic_runtime_pair_threshold_design_only_authorized",
        )
        status = self.record["gate_status"]
        self.assertTrue(status["runtime_refinement_pair_implemented"])
        self.assertTrue(status["runtime_refinement_pair_synthetic_validation_passed"])
        self.assertTrue(status["runtime_threshold_design_authorized"])
        self.assertFalse(status["runtime_thresholds_frozen"])
        self.assertFalse(status["replacement_temporal_admission_defined"])
        self.assertFalse(status["PROTO14_frozen"])

    def test_both_method_transactions_and_continuous_controls_are_bound(self) -> None:
        synthetic = self.record["artifact_payload"]["synthetic_validation"]
        methods = synthetic["method_owned_pair_controls"]
        self.assertEqual(methods["RK4"]["committed_stage_count"], 10)
        self.assertEqual(methods["RK4"]["committed_member_transaction_serial"], 10)
        self.assertEqual(methods["SSPRK3"]["committed_stage_count"], 8)
        self.assertEqual(
            methods["SSPRK3"]["committed_member_transaction_serial"], 8
        )
        for evidence in methods.values():
            self.assertTrue(evidence["coarse_endpoint_differs_from_fine_endpoint"])
            self.assertFalse(evidence["coarse_path_committed"])
            self.assertTrue(evidence["fine_half_steps_committed_atomically"])
            continuous = evidence["continuous_difference"]
            self.assertFalse(continuous["declared_cubic_is_exact_PDE_history"])
            self.assertFalse(continuous["trajectory_asymptotic_regime_proved"])
            self.assertFalse(continuous["rigorous_global_PDE_error_bound_proved"])

        controls = synthetic["continuous_envelope_controls"]
        bump = controls["endpoint_only_bump"]
        cubic = controls["genuine_cubic_two_root_control"]
        self.assertEqual(bump["raw_candidate_maximum"], 1.0)
        self.assertEqual(bump["real_interior_root_count"], 1)
        self.assertEqual(cubic["real_interior_root_count"], 2)
        self.assertTrue(
            controls["dense_audit"]["certified_bound_dominates_dense_audit"]
        )
        self.assertTrue(controls["one_bit_mutation"]["evidence_changed"])

    def test_late_fine_failure_and_commit_drift_are_fail_closed(self) -> None:
        atomic = self.record["artifact_payload"]["synthetic_validation"][
            "atomic_failure_controls"
        ]
        self.assertEqual(atomic["fine_right_injected_stop_path"], "fine_right")
        self.assertEqual(
            atomic["fine_right_injected_stop_reason"], "nonpositive_lapse"
        )
        self.assertTrue(atomic["real_state_monitor_and_causal_ledgers_unchanged"])
        self.assertTrue(atomic["post_prepare_causal_ledger_drift_rejected"])
        self.assertTrue(atomic["no_coarse_or_partial_fine_commit_on_failure"])

    def test_claim_promotion_and_runtime_mutations_fail_closed(self) -> None:
        controls = self.record["artifact_payload"]["mutation_controls"]
        self.assertTrue(controls)
        self.assertEqual(set(controls.values()), {True})
        config = imp1.load_config()
        for mutation in (
            lambda value: value["claims"].__setitem__(
                "runtime_thresholds_frozen", True
            ),
            lambda value: value["claims"].__setitem__("PROTO14_frozen", True),
            lambda value: value["claims"].__setitem__(
                "GR0_case_eligible", True
            ),
            lambda value: value["claims"].__setitem__(
                "FGCQR_holdout_execution_authorized", True
            ),
            lambda value: value["runtime_contract"].__setitem__(
                "coarse_path_is_evidence_only", False
            ),
        ):
            attacked = deepcopy(config)
            mutation(attacked)
            with self.assertRaises(ValueError):
                imp1.validate_config_data(attacked)

    def test_noncanonical_or_promoted_result_is_rejected(self) -> None:
        promoted = deepcopy(self.record)
        promoted["gate_status"]["PROTO14_frozen"] = True
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "promoted.json"
            path.write_text(
                json.dumps(promoted, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            with self.assertRaises(ValueError):
                imp1.verify_canonical(output_path=path)

        noncanonical = deepcopy(self.record)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "noncanonical.json"
            path.write_text(json.dumps(noncanonical), encoding="utf-8")
            with self.assertRaises(ValueError):
                imp1.verify_canonical(output_path=path)


if __name__ == "__main__":
    unittest.main()
