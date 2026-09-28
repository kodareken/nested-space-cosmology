from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_cal2_pref4 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    load_canonical_result,
    load_config,
    record,
)


class FGCCAL2PREF4ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_failure_is_a_test_contract_obstruction_not_a_mechanism_result(self) -> None:
        gates = self.payload["gate_status"]
        self.assertTrue(gates["PROTO5_GR0_source_stop_ownership_contract_obstructed"])
        self.assertTrue(gates["PROTO6_premise_revision_required"])
        self.assertTrue(gates["last_accepted_fine_state_passes_original_source_gate"])
        self.assertTrue(gates["unaccepted_full_proposal_misses_original_source_gate"])
        self.assertTrue(gates["one_half_proposal_passes_unchanged_source_gate"])
        self.assertFalse(gates["fresh_GR0_dynamic_calibration_completed"])
        self.assertFalse(gates["FGCQR_holdout_execution_authorized"])
        self.assertFalse(gates["FGCQR_mechanism_rejected"])

    def test_replay_keeps_the_raw_threshold_and_localizes_the_failure(self) -> None:
        replay = self.payload["artifact_payload"]["diagnostic_replay"]
        self.assertEqual(replay["first_retry_factor_that_passes_every_stage"], "1/2")
        self.assertEqual(replay["failed_unit_proposal_stage_names"], ["rk4_k4"])
        self.assertTrue(replay["failure_localizes_to_first_annular_node"])
        self.assertFalse(replay["proposals"][0]["all_stage_raw_source_gates_passed"])
        self.assertTrue(replay["proposals"][1]["all_stage_raw_source_gates_passed"])

    def test_campaign_hash_mutation_is_rejected(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        mutated = source.replace(
            'manifest_sha256 = "f47c9c54c9ddccb9c1b97807512c3ed1ae93f6692030e0d583ac8b71c6ba55d0"',
            'manifest_sha256 = "047c9c54c9ddccb9c1b97807512c3ed1ae93f6692030e0d583ac8b71c6ba55d0"',
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mutated.toml"
            path.write_text(mutated, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "immutable campaign"):
                load_config(path)


if __name__ == "__main__":
    unittest.main()
