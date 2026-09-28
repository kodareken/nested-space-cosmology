from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import tomllib
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.protocol_v5 import (  # noqa: E402
    validate_sf1_protocol_v5,
)
from scripts.reproduce_fgc_pro5_frz1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)


class FGCProtocolV5Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.protocol_path = REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v5.toml"
        cls.protocol = tomllib.loads(cls.protocol_path.read_text(encoding="utf-8"))

    def test_overlay_validates_and_remains_fail_closed(self) -> None:
        certificate = validate_sf1_protocol_v5(self.protocol)
        self.assertEqual(certificate["artifact_id"], "FGC-2-SF1-PROTO5")
        self.assertEqual(certificate["protocol_version"], 5)
        self.assertTrue(certificate["outcome_neutral_contract_validated"])
        self.assertTrue(certificate["native_semidiscrete_q_initialization_validated"])
        self.assertTrue(certificate["physical_u_and_p_unchanged"])
        self.assertTrue(certificate["method_owned_raw_constraint_guards_validated"])
        self.assertTrue(certificate["roundoff_zero_classification_only_validated"])
        self.assertTrue(certificate["monotone_three_grid_finest_pair_order_validated"])
        self.assertEqual(certificate["eligible_calibration_amplitudes"], ["5/2", "3"])
        self.assertFalse(certificate["resolved_holdout_manifest_present"])
        self.assertTrue(all(value is False for value in certificate["claims"].values()))

    def test_overlay_mutations_fail_closed(self) -> None:
        mutations = (
            (
                ("replacement", "initial_state", "physical_u_and_p_must_be_bitwise_preserved"),
                False,
                "initial-state replacement",
            ),
            (
                ("replacement", "initial_state", "auxiliary_q_initialization"),
                "analytic_q",
                "initial-state replacement",
            ),
            (
                ("replacement", "constraint_admission", "comparator_finest_normalized_guard_max"),
                "1/100",
                "constraint-admission replacement",
            ),
            (
                ("replacement", "constraint_admission", "roundoff_operation_budget"),
                8192,
                "constraint-admission replacement",
            ),
            (
                ("replacement", "provenance", "holdout_output_root"),
                "runs/fgc-2-sf1/proto4/holdout",
                "provenance replacement",
            ),
            (
                ("claims", "FGCQR_holdout_execution_authorized"),
                True,
                "claims",
            ),
        )
        for path, value, message in mutations:
            with self.subTest(path=path):
                mutated = deepcopy(self.protocol)
                target = mutated
                for key in path[:-1]:
                    target = target[key]
                target[path[-1]] = value
                with self.assertRaisesRegex(ValueError, message):
                    validate_sf1_protocol_v5(mutated)

    def test_freeze_config_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            weakened = root / "weakened.toml"
            weakened.write_text(
                source.replace(
                    "only_q_may_be_projected_and_u_p_must_be_bitwise_preserved = true",
                    "only_q_may_be_projected_and_u_p_must_be_bitwise_preserved = false",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "proof contract"):
                load_config(weakened)

            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "PROTO5_fresh_GR0_dynamic_calibration_authorized = false",
                    "PROTO5_fresh_GR0_dynamic_calibration_authorized = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "claims"):
                load_config(promoted)

    def test_canonical_freeze_reproduction(self) -> None:
        payload = load_canonical_result(DEFAULT_OUTPUT)
        self.assertEqual(payload, record(DEFAULT_CONFIG))
        gates = payload["gate_status"]
        self.assertTrue(gates["PROTO5_outcome_neutral_protocol_frozen"])
        self.assertTrue(gates["PROTO5_immutable_lineage_verified"])
        self.assertTrue(gates["PROTO5_native_semidiscrete_initialization_contract_frozen"])
        self.assertFalse(gates["PROTO5_fresh_GR0_dynamic_calibration_authorized"])
        self.assertFalse(gates["FGCQR_holdout_execution_authorized"])
        boundary = payload["artifact_payload"]["epistemic_boundary"]
        self.assertTrue(boundary["initial_common_event_diagnostics_were_seen"])
        self.assertFalse(boundary["fresh_calibration_trajectory_read"])
        self.assertFalse(boundary["FGCQR_outcome_read"])
        self.assertTrue(all(value is False for value in payload["nonclaims"].values()))

    def test_noncanonical_and_duplicate_results_fail(self) -> None:
        payload = record(DEFAULT_CONFIG)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            noncanonical = root / "noncanonical.json"
            noncanonical.write_text(_canonical(payload).rstrip(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(noncanonical)

            duplicate = root / "duplicate.json"
            duplicate.write_text('{"a":1,"a":2}\n', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate JSON key"):
                load_canonical_result(duplicate)


if __name__ == "__main__":
    unittest.main()
