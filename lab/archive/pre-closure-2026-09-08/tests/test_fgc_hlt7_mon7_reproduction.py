from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_hlt4_mon4 import _canonical, _serial  # noqa: E402
from scripts.reproduce_fgc_hlt7_mon7 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    load_canonical_result,
    load_config,
    record,
    validate_run_plan,
    verify_canonical,
)


class FGCHLT7MON7ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.observed = load_canonical_result(DEFAULT_OUTPUT)
        cls.reproduced = _serial(record(DEFAULT_CONFIG))

    def test_canonical_result_reproduces_and_stays_pretrajectory(self) -> None:
        verify_canonical(DEFAULT_OUTPUT, DEFAULT_CONFIG)
        self.assertEqual(self.observed, self.reproduced)
        gates = self.observed["gate_status"]
        self.assertTrue(gates["PROTO9_successor_runtime_implemented"])
        self.assertFalse(gates["PROTO9_fresh_GR0_dynamic_calibration_authorized"])
        self.assertFalse(gates["FGCQR_holdout_execution_authorized"])
        self.assertFalse(gates["retained_EFT_evolution_authorized"])
        self.assertTrue(all(value is False for value in self.observed["nonclaims"].values()))

    def test_direct_inputs_pass_source_constraints_and_absolute_budgets_only(self) -> None:
        payload = self.observed["artifact_payload"]
        inputs = payload["frozen_run_inputs"]
        self.assertEqual(len(inputs), 12)
        finest = [item for item in inputs if item["point_count"] == 8193]
        self.assertEqual(len(finest), 4)
        self.assertTrue(all(item["new_8193_state_constructed"] for item in finest))
        self.assertTrue(all(not item["conditional_projection_used"] for item in inputs))
        overlap = [item for item in inputs if item["point_count"] in {2049, 4097}]
        self.assertTrue(
            all(item["overlap_projected_state_matches_HLT6"] for item in overlap)
        )
        self.assertTrue(
            all(
                item["accepted_state_source_gate_passed"]
                for item in payload["accepted_state_source_prechecks"]
            )
        )
        events = payload["initial_common_events"]
        self.assertEqual(len(events), 4)
        self.assertTrue(all(item["constraint_admission_passed"] for item in events))
        self.assertTrue(
            all(item["finest_pair_spatial_individual_budgets_passed"] for item in events)
        )
        self.assertTrue(
            all(not item["all_adjacent_spatial_tail_ratios_passed"] for item in events)
        )
        self.assertTrue(all(not item["admission_passed"] for item in events))

    def test_config_and_run_plan_broadening_fail_closed(self) -> None:
        config = DEFAULT_CONFIG.read_text(encoding="utf-8")
        plan = (REPOSITORY / "configs/fgc/fgc-1-cal6-run1.toml").read_text(
            encoding="utf-8"
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            promoted = root / "promoted.toml"
            promoted.write_text(
                config.replace(
                    "PROTO9_fresh_GR0_dynamic_calibration_authorized = false",
                    "PROTO9_fresh_GR0_dynamic_calibration_authorized = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "claims"):
                load_config(promoted)
            weakened = root / "weakened.toml"
            weakened.write_text(
                plan.replace(
                    'source_residual_infinity_max = "1/1000000000000"',
                    'source_residual_infinity_max = "1/100000000000"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "changes more"):
                validate_run_plan(weakened)
            relabelled = root / "relabelled.toml"
            relabelled.write_text(
                plan.replace(
                    "resolutions = [2049, 4097, 8193]",
                    "resolutions = [1025, 2049, 4097]",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "resolution ladder"):
                validate_run_plan(relabelled)

    def test_noncanonical_and_duplicate_results_fail(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            noncanonical = root / "noncanonical.json"
            noncanonical.write_text(
                _canonical(self.reproduced).rstrip(), encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(noncanonical)
            duplicate = root / "duplicate.json"
            duplicate.write_text('{"a":1,"a":2}\n', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate JSON key"):
                load_canonical_result(duplicate)


if __name__ == "__main__":
    unittest.main()
