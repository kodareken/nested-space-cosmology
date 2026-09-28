from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

from scripts import reproduce_fgc_tdg3_pref18 as pref18


class FGCTDG3PREF18ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = pref18.verify_canonical()

    def test_canonical_certificate_reproduces_every_native_history_ladder(self) -> None:
        self.assertEqual(self.record["artifact_id"], pref18.ARTIFACT_ID)
        self.assertEqual(
            self.record["classification"],
            "completed_TDG3_mixed_native_grid_history_diagnosis_"
            "sampling_identifiability_theorem_required",
        )
        diagnosis = self.record["artifact_payload"]["actual_history_diagnosis"]
        self.assertEqual(
            sum(len(records) for records in diagnosis["signal_records"].values()),
            576,
        )
        self.assertEqual(
            diagnosis["frozen_estimator_names"], list(pref18.TDG3_ESTIMATOR_NAMES)
        )
        self.assertEqual(
            diagnosis["frozen_power_class_names"],
            list(pref18.TDG3_POWER_CLASS_NAMES),
        )
        self.assertEqual(
            diagnosis["frozen_combined_class_names"],
            list(pref18.TDG3_COMBINED_CLASS_NAMES),
        )

    def test_compact_mixed_outcome_and_claim_boundary_are_exact(self) -> None:
        payload = self.record["artifact_payload"]
        compact = payload["compact_findings"]
        self.assertEqual(compact, pref18._expected_observed())
        self.assertEqual(compact["native_estimator_agreement_count"], 1119)
        self.assertEqual(compact["native_estimator_disagreement_count"], 33)
        self.assertEqual(compact["native_estimators_agree_unresolved"], 554)
        self.assertEqual(compact["TDG2_finest_zero_survive_native_zero"], 92)
        self.assertEqual(
            compact["TDG2_finest_nonzero_survive_native_nonzero"], 0
        )
        self.assertEqual(
            compact["TDG2_finest_enclosure_remain_native_unresolved"], 101
        )
        self.assertEqual(
            compact["TDG2_finest_enclosure_native_disagreement"], 3
        )
        self.assertTrue(
            payload["conclusions"][
                "common_grid_interpolation_was_a_real_but_not_sole_owner"
            ]
        )
        self.assertFalse(
            payload["claim_boundary"]["replacement_temporal_admission_defined"]
        )
        self.assertFalse(payload["claim_boundary"]["new_trajectory_authorized"])
        self.assertFalse(payload["claim_boundary"]["candidate_branch_opened"])
        self.assertFalse(payload["claim_boundary"]["physical_question_answered"])

    def test_all_TDG2_zero_contracting_finest_failures_survive(self) -> None:
        compact = self.record["artifact_payload"]["compact_findings"]
        self.assertEqual(compact["TDG2_finest_zero_classifications"], 92)
        self.assertEqual(
            compact["TDG2_finest_zero_survive_native_zero"],
            compact["TDG2_finest_zero_classifications"],
        )
        self.assertTrue(
            self.record["artifact_payload"]["conclusions"][
                "all_TDG2_zero_contracting_finest_failures_survive_both_native_estimators"
            ]
        )

    def test_no_signal_claims_resampling_or_continuum_enclosure(self) -> None:
        diagnosis = self.record["artifact_payload"]["actual_history_diagnosis"]
        self.assertFalse(diagnosis["common_grid_resampling_performed"])
        self.assertFalse(diagnosis["continuum_surrogate_enclosure_claimed"])
        for records in diagnosis["signal_records"].values():
            for item in records:
                self.assertFalse(item["common_grid_resampling_performed"])
                self.assertFalse(item["continuum_surrogate_enclosure_claimed"])
                for estimator in pref18.TDG3_ESTIMATOR_NAMES:
                    for observation in item["estimators"][estimator]["observations"]:
                        self.assertFalse(
                            observation["common_grid_resampling_performed"]
                        )
                        self.assertFalse(
                            observation["continuum_surrogate_enclosure_claimed"]
                        )

    def test_tracked_result_is_canonical(self) -> None:
        raw = pref18.DEFAULT_OUTPUT.read_bytes()
        stored = json.loads(raw)
        self.assertEqual(raw, pref18._canonical_bytes(stored))

    def test_complete_absence_is_allowed_but_partial_raw_bundle_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = tuple(root / name for name in ("a", "b", "c", "d"))
            self.assertFalse(pref18._raw_bundle_is_complete_or_absent(paths))
            paths[0].write_bytes(b"present")
            with self.assertRaisesRegex(ValueError, "bundle is partial"):
                pref18._raw_bundle_is_complete_or_absent(paths)
            for path in paths[1:]:
                path.write_bytes(b"present")
            self.assertTrue(pref18._raw_bundle_is_complete_or_absent(paths))

    def test_config_mutations_fail_closed(self) -> None:
        config = pref18.load_config()
        attacks = []
        value = deepcopy(config)
        value["immutable_lineage"]["freeze_result_sha256"] = "0" * 64
        attacks.append(value)
        value = deepcopy(config)
        value["observed_diagnosis"]["native_estimator_disagreement_count"] = 32
        attacks.append(value)
        value = deepcopy(config)
        value["scope"]["common_grid_resampling_performed"] = True
        attacks.append(value)
        value = deepcopy(config)
        value["claims"][
            "piecewise_linear_or_quadrature_surrogate_is_a_continuum_enclosure"
        ] = True
        attacks.append(value)
        value = deepcopy(config)
        value["claims"]["replacement_temporal_admission_defined"] = True
        attacks.append(value)
        value = deepcopy(config)
        value["claims"]["GR0_case_eligible"] = True
        attacks.append(value)
        for attacked in attacks:
            with self.assertRaises(ValueError):
                pref18.validate_config_data(attacked)

    def test_all_embedded_mutation_controls_pass(self) -> None:
        controls = self.record["artifact_payload"]["mutation_controls"]
        self.assertTrue(controls)
        self.assertEqual(set(controls.values()), {True})


if __name__ == "__main__":
    unittest.main()
