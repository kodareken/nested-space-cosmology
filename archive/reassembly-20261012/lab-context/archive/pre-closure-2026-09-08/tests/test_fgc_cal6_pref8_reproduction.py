from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_cal6_pref8 as pref8  # noqa: E402


class FGCCAL6PREF8ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.canonical = json.loads(pref8.DEFAULT_OUTPUT.read_text(encoding="utf-8"))

    def test_canonical_result_reproduces_and_stays_pretrajectory(self) -> None:
        self.assertEqual(pref8.record(), self.canonical)
        payload = self.canonical["artifact_payload"]
        aggregate = payload["aggregate"]
        self.assertEqual(aggregate["raw_failed_finest_pair_metric_count"], 16)
        self.assertEqual(
            aggregate["prospective_diagnostically_saturated_metric_count"], 16
        )
        self.assertEqual(
            aggregate["prospective_directly_resolved_finest_pair_metric_count"],
            32,
        )
        self.assertTrue(aggregate["all_raw_PROTO9_admissions_failed"])
        self.assertTrue(
            aggregate["all_prospective_resolved_or_saturated_admissions_passed"]
        )
        self.assertLess(aggregate["maximum_top_band_erasure_over_round_trip"], 0.006)
        self.assertFalse(payload["decision"]["current_campaign_authorized"])
        self.assertFalse(payload["epistemic_boundary"]["trajectory_or_candidate_outcome_read"])
        self.assertTrue(all(value is False for value in self.canonical["nonclaims"].values()))

    def test_every_case_preserves_the_raw_veto_and_requires_guarded_saturation(
        self,
    ) -> None:
        cases = self.canonical["artifact_payload"]["case_diagnostics"]
        self.assertEqual(
            [(item["amplitude"], item["method"]) for item in cases],
            [("5/2", "RK4"), ("5/2", "SSPRK3"), ("3", "RK4"), ("3", "SSPRK3")],
        )
        for item in cases:
            self.assertFalse(item["raw_PROTO9_admission"]["admission_passed"])
            prospective = item["prospective_resolved_or_saturated_admission"]
            self.assertTrue(prospective["admission_passed"])
            self.assertTrue(prospective["saturation_used"])
            self.assertTrue(prospective["every_profile_contracted"])
            self.assertTrue(
                all(
                    prospective["coarse_to_medium_direct_passed_by_field"].values()
                )
            )
            self.assertTrue(
                all(
                    prospective[
                        "coarse_to_medium_direct_passed_by_derivative"
                    ].values()
                )
            )

    def test_contract_broadening_and_result_promotion_fail_closed(self) -> None:
        source = pref8.DEFAULT_CONFIG.read_text(encoding="utf-8")
        broadened = source.replace(
            "coarse_to_medium_field_and_derivative_ratios_must_pass_directly = true",
            "coarse_to_medium_field_and_derivative_ratios_must_pass_directly = false",
        )
        self.assertNotEqual(source, broadened)
        with tempfile.TemporaryDirectory() as directory:
            config_path = Path(directory) / "mutated.toml"
            config_path.write_text(broadened, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "diagnostic definition differs"):
                pref8.load_config(config_path)

            promoted = deepcopy(self.canonical)
            promoted["gate_status"][
                "PROTO10_fresh_GR0_dynamic_calibration_authorized"
            ] = True
            result_path = Path(directory) / "promoted.json"
            result_path.write_text(
                json.dumps(promoted, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "differs from a fresh reproduction"):
                pref8.verify_canonical(result_path)


if __name__ == "__main__":
    unittest.main()
