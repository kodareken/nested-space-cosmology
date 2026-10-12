"""Compact and adversarial controls for LOC1 ERR1."""

from __future__ import annotations

from copy import deepcopy
import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("loc1_err1", ROOT / "scripts/reproduce_fgc_tdg9_loc1_err1.py")
assert SPEC is not None and SPEC.loader is not None
err1 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(err1)


class LOC1ERR1Tests(unittest.TestCase):
    def test_compact_result_is_exact_and_store_blind(self) -> None:
        raw = (ROOT / err1.RESULT_PATH).read_bytes()
        value = err1.validate(raw)
        self.assertEqual(value["exception"]["detail"], [3, "u:alpha", "value_V", "D01"])
        self.assertEqual(value["diagnosis"]["primary_candidate_count"], 8176)
        self.assertEqual(value["diagnosis"]["independent_candidate_count"], 8544)
        self.assertFalse(value["scope"]["physical_result_earned"])

    def test_mutations_fail_closed(self) -> None:
        value = err1._expected()
        for mutation in (
            lambda item: item["diagnosis"].__setitem__("affected_cubic_count", 183),
            lambda item: item["diagnosis"]["first_witness"].__setitem__("exact_derivative_roots", ["0"]),
            lambda item: item["scope"].__setitem__("mechanism_result_earned", True),
            lambda item: item["construction_evidence"].__setitem__("LOC1_terminal_absent", False),
        ):
            changed = deepcopy(value)
            mutation(changed)
            with self.assertRaises(err1.ERR1Error):
                err1.validate(err1._canonical(changed))

    def test_sealed_tuple_is_reconstructed_without_loc2(self) -> None:
        source = (ROOT / "scripts/reproduce_fgc_tdg9_loc1_err1.py").read_text()
        self.assertNotIn("tdg9_local_extrema_independent_v2", source)
        reproduced = err1._reconstruct_diagnosis()
        expected = err1._expected()["diagnosis"]
        self.assertEqual(reproduced["primary_candidate_count"], 8176)
        self.assertEqual(reproduced["independent_candidate_count"], 8544)
        self.assertEqual(reproduced["affected_cubic_count"], 184)
        self.assertEqual(reproduced["false_extra_candidate_count"], 368)
        self.assertEqual(reproduced["co_maximizers"], expected["co_maximizers"])
        self.assertEqual(reproduced["first_witness"], expected["first_witness"])


if __name__ == "__main__":
    unittest.main()
