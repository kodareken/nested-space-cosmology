from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

from scripts import reproduce_fgc_tdg2_pref17 as pref17


class FGCTDG2PREF17ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = pref17.verify_canonical()

    def test_canonical_certificate_reproduces_every_history_ladder(self) -> None:
        self.assertEqual(self.record["artifact_id"], pref17.ARTIFACT_ID)
        self.assertEqual(
            self.record["classification"],
            "completed_TDG2_mixed_absolute_tail_diagnosis_TDG3_required",
        )
        diagnosis = self.record["artifact_payload"]["actual_history_diagnosis"]
        self.assertEqual(
            sum(len(records) for records in diagnosis["signal_records"].values()),
            576,
        )
        self.assertEqual(diagnosis["frozen_class_names"], list(pref17.TDG2_CLASS_NAMES))
        self.assertEqual(
            diagnosis["frozen_measure_names"], list(pref17.TDG2_MEASURE_NAMES)
        )

    def test_compact_mixed_outcome_and_claim_boundary_are_exact(self) -> None:
        payload = self.record["artifact_payload"]
        compact = payload["compact_findings"]
        self.assertEqual(compact, pref17._expected_observed())
        self.assertEqual(compact["legacy_finest_individual_failures"], 152)
        self.assertEqual(compact["all_three_below_DFT_floor_ladders"], 96)
        self.assertEqual(
            compact["enclosure_dominated_finest_failures_field"], 101
        )
        self.assertEqual(
            compact["enclosure_dominated_finest_failures_derivative"], 102
        )
        self.assertEqual(compact["nonzero_convergent_classifications"], 0)
        self.assertEqual(compact["unresolved_classifications"], 0)
        self.assertTrue(payload["conclusions"]["absolute_tail_outcome_is_mixed"])
        self.assertFalse(
            payload["claim_boundary"]["replacement_temporal_admission_defined"]
        )
        self.assertFalse(payload["claim_boundary"]["new_trajectory_authorized"])
        self.assertFalse(payload["claim_boundary"]["candidate_branch_opened"])
        self.assertFalse(payload["claim_boundary"]["physical_question_answered"])

    def test_arithmetic_does_not_own_a_finest_legacy_failure(self) -> None:
        diagnosis = self.record["artifact_payload"]["actual_history_diagnosis"]
        for method in ("RK4", "SSPRK3"):
            for measure in pref17.TDG2_MEASURE_NAMES:
                causes = diagnosis["methods"][method][
                    "finest_failure_enclosure_cause_counts"
                ][measure]
                self.assertEqual(
                    causes.get("finest_arithmetic_interval_reaches_zero", 0), 0
                )
        self.assertTrue(
            self.record["artifact_payload"]["conclusions"][
                "interpolation_debit_is_the_dominant_unresolved_owner"
            ]
        )

    def test_tracked_result_is_canonical(self) -> None:
        raw = pref17.DEFAULT_OUTPUT.read_bytes()
        stored = json.loads(raw)
        self.assertEqual(raw, pref17._canonical_bytes(stored))

    def test_complete_absence_is_allowed_but_partial_raw_bundle_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = tuple(root / name for name in ("a", "b", "c", "d"))
            self.assertFalse(pref17._raw_bundle_is_complete_or_absent(paths))
            paths[0].write_bytes(b"present")
            with self.assertRaisesRegex(ValueError, "bundle is partial"):
                pref17._raw_bundle_is_complete_or_absent(paths)
            for path in paths[1:]:
                path.write_bytes(b"present")
            self.assertTrue(pref17._raw_bundle_is_complete_or_absent(paths))

    def test_config_mutations_fail_closed(self) -> None:
        config = pref17.load_config()
        attacks = []
        value = deepcopy(config)
        value["immutable_lineage"]["freeze_result_sha256"] = "0" * 64
        attacks.append(value)
        value = deepcopy(config)
        value["observed_diagnosis"]["legacy_finest_individual_failures"] = 151
        attacks.append(value)
        value = deepcopy(config)
        value["conclusions"][
            "interpolation_debit_is_the_dominant_unresolved_owner"
        ] = False
        attacks.append(value)
        value = deepcopy(config)
        value["claims"]["replacement_temporal_admission_defined"] = True
        attacks.append(value)
        value = deepcopy(config)
        value["claims"]["GR0_case_eligible"] = True
        attacks.append(value)
        for attacked in attacks:
            with self.assertRaises(ValueError):
                pref17.validate_config_data(attacked)

    def test_all_embedded_mutation_controls_pass(self) -> None:
        controls = self.record["artifact_payload"]["mutation_controls"]
        self.assertTrue(controls)
        self.assertEqual(set(controls.values()), {True})


if __name__ == "__main__":
    unittest.main()
