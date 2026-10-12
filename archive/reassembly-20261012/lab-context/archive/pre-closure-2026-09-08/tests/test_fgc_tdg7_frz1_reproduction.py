from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from scripts import reproduce_fgc_tdg7_frz1 as reproduce  # noqa: E402


class TDG7FRZ1ReproductionTests(unittest.TestCase):
    def test_canonical_tracked_record_and_closed_boundary(self) -> None:
        record = reproduce.build(reproduce.load_config())
        tracked = json.loads(reproduce.OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(record, tracked)
        self.assertEqual(record["artifact_id"], "FGC-1-TDG7-FRZ1")
        self.assertEqual(record["gate_status"], "PASS_EXACT_BINARY64_SUBDIVISION_LATTICE_DESIGN_FROZEN")
        boundary = record["artifact_payload"]["claim_boundary"]
        self.assertTrue(boundary["TDG7_subdivision_lattice_design_frozen"])
        for name in ("TDG7_runtime_repair_implemented", "fresh_GR0_calibration_authorized",
                     "GR0_case_eligible", "FGCQR_holdout_execution_authorized"):
            self.assertFalse(boundary[name])
        plan = record["artifact_payload"]["selected_lattice_plan"]
        self.assertEqual(plan["event_quantum_hex"], "0x1.0000000000000p-52")
        self.assertEqual(plan["macro_quantum_hex"], "0x1.0000000000000p-49")
        self.assertEqual(plan["macro_eight_Q_tick_count"], 3664984285035)
        self.assertEqual(plan["fine_width_over_Q"], 7329968570070)
        self.assertTrue(plan["stage_times_exact_and_strictly_interior"])
        self.assertEqual(
            record["artifact_payload"]["historical_vs_selected_first_macro_endpoint"],
            {
                "historical_expected_endpoint_hex": "0x1.71aaa90b0fb5cp+0",
                "selected_aligned_endpoint_hex": "0x1.71aaa90b0fb58p+0",
                "equal": False,
                "event_target_unchanged": True,
            },
        )
        attack = record["artifact_payload"]["four_q_endpoint_only_midpoint_attack"]
        self.assertEqual(attack["endpoint_only_fine_width_over_Q"], 7329968570071)
        self.assertFalse(attack["midpoint_exact_binary64"])

    def test_cli_verify_matches_result(self) -> None:
        completed = subprocess.run([sys.executable, "scripts/reproduce_fgc_tdg7_frz1.py", "--verify"],
                                   cwd=ROOT, check=True, capture_output=True, text=True)
        self.assertIn('"artifact_id": "FGC-1-TDG7-FRZ1"', completed.stdout)

    def test_mutations_fail_closed(self) -> None:
        config = reproduce.load_config()
        variants = []
        changed = deepcopy(config); changed["immutable_lineage"]["source_result_sha256"] = "0" * 64; variants.append(changed)
        changed = deepcopy(config); changed["selected_route"]["rounding_up_is_forbidden"] = False; variants.append(changed)
        changed = deepcopy(config); changed["historical_event24_witness"]["macro_quantum_hex"] = "0x1.0000000000000p-50"; variants.append(changed)
        changed = deepcopy(config); changed["historical_event24_witness"]["requested_macro_width_hex"] = "0x1.0000000000000p-8"; variants.append(changed)
        changed = deepcopy(config); changed["selected_route"]["macro_quantum"] = "4*event_quantum"; variants.append(changed)
        changed = deepcopy(config); changed["claims"]["FGCQR_holdout_execution_authorized"] = True; variants.append(changed)
        for value in variants:
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    reproduce.validate_config_data(value)

    def test_optional_checkpoint_presence_does_not_change_certificate_contract(self) -> None:
        config = reproduce.load_config()
        expected = {
            "whole_file_hash_required_before_metadata_read": True,
            "metadata_utf8_only": True,
            "physical_state_arrays_consumed": False,
            "clean_clone_absence_permitted": True,
        }
        self.assertEqual(reproduce._optional_checkpoint(config), expected)
        with mock.patch("pathlib.Path.is_file", return_value=False):
            self.assertEqual(reproduce._optional_checkpoint(config), expected)


if __name__ == "__main__":
    unittest.main()
