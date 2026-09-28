from __future__ import annotations

import math
import unittest

from scripts import reproduce_fgc_rsp2_frz1 as rsp2


class RSP2FreezeReproductionTests(unittest.TestCase):
    def test_canonical_record_reproduces(self) -> None:
        record = rsp2.verify_canonical()
        self.assertEqual(record["artifact_id"], rsp2.ARTIFACT_ID)
        self.assertTrue(record["gate_status"]["RSP2_execution_authorized"])
        self.assertFalse(record["gate_status"]["RSP2_outcome_read"])

    def test_claim_mutation_fails_closed(self) -> None:
        plan = rsp2.validate_run_plan()
        plan["claims"]["FGCQR_holdout_execution_authorized"] = True
        with self.assertRaises(ValueError):
            rsp2.validate_run_plan_data(plan)

    def test_subthreshold_reduction_last_bits_do_not_mutate_study_identity(self) -> None:
        first = rsp2._source_precheck_certificate(5.306866057708248e-14, 1e-12)
        second = rsp2._source_precheck_certificate(5.3235194030776256e-14, 1e-12)
        self.assertEqual(first, second)
        self.assertFalse(first["raw_binary64_reduction_used_in_study_identity"])
        self.assertTrue(first["raw_value_preserved_by_runtime_diagnostics"])
        for invalid in (-1.0, math.inf, math.nan, 2e-12):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    rsp2._source_precheck_certificate(invalid, 1e-12)

    def test_method_time_component_and_grid_mutations_fail_closed(self) -> None:
        mutations = {
            "method": ("method", "method_label", "RK4"),
            "time": ("numerics", "final_coordinate_time", "3/2"),
            "component": ("question", "target_component", "gauge_r"),
            "grid": ("numerics", "new_resolution", 8193),
        }
        for name, (table, key, value) in mutations.items():
            with self.subTest(name=name):
                plan = rsp2.validate_run_plan()
                plan[table][key] = value
                with self.assertRaises(ValueError):
                    rsp2.validate_run_plan_data(plan)


if __name__ == "__main__":
    unittest.main()
