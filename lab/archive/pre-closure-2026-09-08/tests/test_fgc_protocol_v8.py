from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import tomllib
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.protocol_v8 import (  # noqa: E402
    validate_sf1_protocol_v8,
)
from scripts.reproduce_fgc_pro8_frz1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)


class FGCProtocolV8Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.protocol_path = (
            REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v8.toml"
        )
        cls.protocol = tomllib.loads(cls.protocol_path.read_text(encoding="utf-8"))

    def test_overlay_validates_and_remains_fail_closed(self) -> None:
        certificate = validate_sf1_protocol_v8(self.protocol)
        self.assertEqual(certificate["artifact_id"], "FGC-2-SF1-PROTO8")
        self.assertEqual(certificate["protocol_version"], 8)
        self.assertTrue(certificate["outcome_neutral_contract_validated"])
        self.assertEqual(certificate["projector_fixed_outer_rows"], 4)
        self.assertEqual(certificate["RK4_excluded_outer_rows"], 7)
        self.assertEqual(certificate["SSPRK3_excluded_outer_rows"], 5)
        self.assertTrue(certificate["accumulated_roundoff_changes_order_only"])
        self.assertTrue(certificate["raw_constraint_magnitude_guards_unchanged"])
        self.assertTrue(certificate["coarsest_spectrum_is_public_convergence_witness"])
        self.assertTrue(certificate["finest_pair_absolute_spectral_budgets_required"])
        self.assertTrue(certificate["all_adjacent_nested_tail_ratios_required"])
        self.assertTrue(certificate["all_physical_and_numerical_thresholds_unchanged"])
        self.assertFalse(certificate["resolved_holdout_manifest_present"])
        self.assertTrue(all(value is False for value in certificate["claims"].values()))

    def test_overlay_mutations_fail_closed(self) -> None:
        mutations = (
            (
                (
                    "replacement",
                    "common_event_constraint_ownership",
                    "projector_fixed_outer_rows",
                ),
                3,
                "constraint-ownership replacement",
            ),
            (
                (
                    "replacement",
                    "common_event_constraint_ownership",
                    "accumulated_roundoff_may_change_only_convergence_order_classification",
                ),
                False,
                "constraint-ownership replacement",
            ),
            (
                (
                    "replacement",
                    "common_event_spectral_ownership",
                    "medium_and_finest_grids_each_must_pass_every_unchanged_absolute_budget",
                ),
                False,
                "spectral-ownership replacement",
            ),
            (
                (
                    "replacement",
                    "common_event_spectral_ownership",
                    "maximum_nested_tail_ratio",
                ),
                "1/2",
                "spectral-ownership replacement",
            ),
            (
                ("replacement", "provenance", "holdout_output_root"),
                "runs/fgc-2-sf1/proto7/holdout",
                "provenance replacement",
            ),
            (("claims", "FGCQR_holdout_execution_authorized"), True, "claims"),
        )
        for path, value, message in mutations:
            with self.subTest(path=path):
                mutated = deepcopy(self.protocol)
                target = mutated
                for key in path[:-1]:
                    target = target[key]
                target[path[-1]] = value
                with self.assertRaisesRegex(ValueError, message):
                    validate_sf1_protocol_v8(mutated)

    def test_freeze_config_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            weakened = root / "weakened.toml"
            weakened.write_text(
                source.replace(
                    "every_adjacent_nested_tail_must_still_contract = true",
                    "every_adjacent_nested_tail_must_still_contract = false",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "proof contract"):
                load_config(weakened)

            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "PROTO8_fresh_GR0_dynamic_calibration_authorized = false",
                    "PROTO8_fresh_GR0_dynamic_calibration_authorized = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "claims"):
                load_config(promoted)

    def test_canonical_freeze_reproduction(self) -> None:
        payload = load_canonical_result(DEFAULT_OUTPUT)
        self.assertEqual(payload, record(DEFAULT_CONFIG))
        gates = payload["gate_status"]
        self.assertTrue(gates["PROTO8_outcome_neutral_protocol_frozen"])
        self.assertTrue(gates["PROTO8_immutable_lineage_verified"])
        self.assertTrue(gates["PROTO8_evolution_owned_constraint_contract_frozen"])
        self.assertTrue(gates["PROTO8_accumulated_roundoff_order_contract_frozen"])
        self.assertTrue(gates["PROTO8_finest_pair_nested_spectral_contract_frozen"])
        self.assertFalse(gates["PROTO8_fresh_GR0_dynamic_calibration_authorized"])
        self.assertFalse(gates["FGCQR_holdout_execution_authorized"])
        boundary = payload["artifact_payload"]["epistemic_boundary"]
        self.assertTrue(boundary["PROTO7_GR0_calibration_trajectory_was_seen"])
        self.assertTrue(boundary["first_nonzero_common_event_completed"])
        self.assertFalse(boundary["PROTO8_trajectory_read"])
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
