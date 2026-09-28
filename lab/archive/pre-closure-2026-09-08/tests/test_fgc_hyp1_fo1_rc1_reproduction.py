from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from scripts.reproduce_fgc_hyp1_fo1_rc1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    load_config,
    record,
)


class FGCHYP1FO1RC1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = record()

    def test_frozen_config_and_canonical_result(self) -> None:
        config = load_config()
        self.assertEqual(config.artifact_id, "FGC-1-HYP1-FO1-RC1")
        self.assertEqual(
            [config.flat_fixture["fixture_id"], config.activated_fixture["fixture_id"]],
            ["FGCQR_flat_reference_vacuum", "FGCQR_activated_generic"],
        )
        committed = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(self.payload, committed)
        certificate = self.payload["first_order_certificate"]
        self.assertTrue(certificate["all_declared_exact_checks_pass"])
        self.assertTrue(certificate["exact_local_first_order_dae_lift_derived"])
        self.assertTrue(
            certificate["complete_flat_root_linearized_first_order_map_derived"]
        )
        self.assertTrue(
            certificate["kinematic_radial_reduction_constraint_identity_derived"]
        )
        self.assertEqual(
            certificate["controls"]["activated_off_shell"]["solution_status"],
            "off_shell_control_not_a_solution",
        )
        self.assertFalse(self.payload["gate_status"]["evolution_authorized"])
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_unknown_key_traversal_order_and_promoted_gates_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            unknown = root / "unknown.toml"
            unknown.write_text(source + "\nunknown_key = true\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "keys differ"):
                load_config(unknown)

            traversal = root / "traversal.toml"
            traversal.write_text(
                source.replace(
                    'propagation_config = "configs/fgc/fgc-1-hyp1-mhg-propagation.toml"',
                    'propagation_config = "../fgc-1-hyp1-mhg-propagation.toml"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "stay inside|does not exist"):
                load_config(traversal)

            noncanonical = root / "noncanonical.toml"
            noncanonical.write_text(
                source.replace(
                    'radial_domain_minimum = "1/2"',
                    'radial_domain_minimum = "2/4"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "canonical rational"):
                load_config(noncanonical)

            reordered = root / "reordered.toml"
            reordered.write_text(
                source.replace(
                    'state_order = ["u.h_tt", "u.h_tr"',
                    'state_order = ["u.h_tr", "u.h_tt"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "frozen ordered string list"):
                load_config(reordered)

            promoted_solver = root / "promoted-solver.toml"
            promoted_solver.write_text(
                source.replace(
                    "nonlinear_acceleration_solver_implemented = false",
                    "nonlinear_acceleration_solver_implemented = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "cannot promote"):
                load_config(promoted_solver)

            promoted_solution = root / "promoted-solution.toml"
            promoted_solution.write_text(
                source.replace(
                    'activated_background_solution_status = "off_shell_control_not_a_solution"',
                    'activated_background_solution_status = "solution"',
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "off shell"):
                load_config(promoted_solution)

            promoted_gate = root / "promoted-gate.toml"
            promoted_gate.write_text(
                source.replace(
                    "evolution_authorized = false", "evolution_authorized = true"
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "open gate"):
                load_config(promoted_gate)


if __name__ == "__main__":
    unittest.main()
