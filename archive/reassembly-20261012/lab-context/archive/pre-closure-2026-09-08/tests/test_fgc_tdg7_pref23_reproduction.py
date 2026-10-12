from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from scripts import reproduce_fgc_tdg7_pref23 as pref23  # noqa: E402


class TDG7PREF23ReproductionTests(unittest.TestCase):
    def test_canonical_binder_authorizes_only_runtime_repair_implementation(self) -> None:
        record = pref23.build(pref23.load_config())
        tracked = json.loads(pref23.OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(record, tracked)
        self.assertEqual(record["artifact_id"], "FGC-1-TDG7-PREF23")
        self.assertEqual(record["gate_status"], "PASS_TDG7_RUNTIME_REPAIR_IMPLEMENTATION_AUTHORIZED")
        binder = record["artifact_payload"]["independent_binder"]
        self.assertTrue(binder["all_independent_binder_controls_pass"])
        self.assertFalse(binder["theorem_module_imports_TDG7_design_helper"])
        self.assertTrue(binder["TDG7_independent_binder_completed"])
        self.assertTrue(binder["TDG7_runtime_repair_implementation_authorized"])
        self.assertFalse(binder["TDG7_runtime_repair_implemented"])
        boundary = record["artifact_payload"]["claim_boundary"]
        for name in ("TDG7_runtime_repair_implemented", "PROTO14_runtime_mutation_authorized",
                     "fresh_GR0_calibration_authorized", "GR0_case_eligible",
                     "SGBL_execution_authorized", "FGCQR_holdout_execution_authorized",
                     "DEF1_execution_authorized"):
            self.assertFalse(boundary[name], name)

    def test_independent_witness_and_source_stage_contract_are_exact(self) -> None:
        binder = pref23.build(pref23.load_config())["artifact_payload"]["independent_binder"]
        witness = binder["independent_CAL11_event24_witness"]
        self.assertEqual(witness["macro_quantum_hex"], "0x1.0000000000000p-49")
        self.assertEqual(witness["fine_width_over_Q"], 7329968570070)
        self.assertTrue(
            witness["endpoint_only_four_q_counterexample_derived_from_historical_one_step"]
        )
        self.assertFalse(witness["endpoint_only_four_q_midpoint_exact"])
        self.assertEqual(witness["stage_triplet_count"], 7)
        engine = binder["independent_engine_stage_abscissa_contract"]
        self.assertTrue(engine["candidate_endpoint_uses_final_time"])
        self.assertIn(["rk4_k2", "midpoint"], engine["stage_time_kinds"])
        self.assertIn(["ssprk3_s2", "midpoint"], engine["stage_time_kinds"])
        historical = binder["independent_historical_TDG6_source_contract"]
        self.assertTrue(historical["step_is_width_over_count"])
        self.assertTrue(historical["boundaries_are_start_plus_index_times_step"])
        self.assertTrue(historical["expected_final_is_start_plus_width"])
        self.assertTrue(historical["macro_endpoint_guard_present"])
        self.assertTrue(historical["adjacent_uniformity_guard_present"])
        controls = pref23.build(pref23.load_config())["artifact_payload"]["mutation_controls"]
        self.assertTrue(all(controls.values()))
        self.assertTrue(
            controls["historical_source_boundary_expression_mutation_rejected"]
        )

    def test_config_mutations_fail_closed(self) -> None:
        config = pref23.load_config()
        attacks = []
        value = deepcopy(config); value["immutable_lineage"]["freeze_result_sha256"] = "0" * 64; attacks.append(value)
        value = deepcopy(config); value["independent_stage_safe_theorem"]["stage_safe_macro_quantum"] = "4*Q"; attacks.append(value)
        value = deepcopy(config); value["historical_event24_witness"]["stage_safe_eight_Q_tick_count"] = 1; attacks.append(value)
        value = deepcopy(config); value["runtime_repair_authorization"]["fresh_GR0_calibration_authorized"] = True; attacks.append(value)
        value = deepcopy(config); value["claims"]["TDG7_runtime_repair_implemented"] = True; attacks.append(value)
        value = deepcopy(config); value["claims"]["FGCQR_holdout_execution_authorized"] = True; attacks.append(value)
        for attacked in attacks:
            with self.subTest(attacked=attacked):
                with self.assertRaises(ValueError):
                    pref23.validate_config_data(attacked)

    def test_cli_verify_matches_tracked_result(self) -> None:
        completed = subprocess.run(
            [sys.executable, "scripts/reproduce_fgc_tdg7_pref23.py", "--verify"],
            cwd=ROOT, check=True, capture_output=True, text=True,
        )
        self.assertIn('"artifact_id": "FGC-1-TDG7-PREF23"', completed.stdout)


if __name__ == "__main__":
    unittest.main()
