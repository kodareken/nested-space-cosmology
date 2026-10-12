from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_cal1_pref3 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    load_canonical_result,
    load_config,
    record,
)


class FGCCAL1PREF3ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_scope_remains_pretrajectory_and_fail_closed(self) -> None:
        gate = self.payload["gate_status"]
        self.assertTrue(gate["PROTO4_semidiscrete_common_event_contract_obstructed"])
        self.assertTrue(gate["PROTO5_premise_revision_required"])
        self.assertFalse(gate["PROTO4_fresh_GR0_dynamic_calibration_authorized"])
        boundary = self.payload["artifact_payload"]["epistemic_boundary"]
        self.assertFalse(boundary["trajectory_read"])
        self.assertFalse(boundary["mechanism_question_answered"])

    def test_config_mutation_is_rejected(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        mutated = source.replace(
            'comparator_finest_guard_max = "1/200"',
            'comparator_finest_guard_max = "1/100"',
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mutated.toml"
            path.write_text(mutated, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "prospective repair"):
                load_config(path)


if __name__ == "__main__":
    unittest.main()
