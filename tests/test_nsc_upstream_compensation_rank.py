"""The rank gate checks algebra and exact margins without field solves."""
from fractions import Fraction as Q
import importlib.util
import json
from pathlib import Path
import unittest

PATH = Path(__file__).resolve().parents[1] / "scripts/derive_nsc_upstream_compensation_rank.py"
SPEC = importlib.util.spec_from_file_location("rank_gate", PATH)
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)


class RankGateTests(unittest.TestCase):
    def test_symbolic_operator_inverse_and_derivative(self):
        self.assertTrue(all(gate.exact_identities().values()))

    def test_exact_margin_and_authenticated_receipt(self):
        record = gate.build_record()
        self.assertEqual(json.loads(gate.OUTPUT.read_text()), record)
        self.assertEqual(Q(record["bounds"]["normalized_response_sigma_min_strict_lower"]), Q(51, 200))
        self.assertEqual(record["scope"]["full_C0_matching"], "OPEN")
        self.assertFalse(record["scope"]["amplitudes_selected"])

    def test_wrong_label_domain_cannot_reuse_gate(self):
        for labels in ((Q(0), Q(9, 4), Q(11, 20)),
                       (Q(31, 20), Q(0), Q(11, 20)),
                       (Q(31, 20), Q(11, 5), Q(0))):
            with self.assertRaises(ValueError):
                gate.guard_labels(*labels)


if __name__ == "__main__":
    unittest.main()
