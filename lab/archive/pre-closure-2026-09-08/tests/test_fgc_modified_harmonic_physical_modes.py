from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


sys.path[:0] = [
    str(Path(__file__).resolve().parents[1]),
    str(Path(__file__).resolve().parents[1] / "src"),
]

from scripts.reproduce_fgc_hyp1_fo1_rc1 import load_config  # noqa: E402
from recursive_horizons.fgc.modified_harmonic_constraints import (  # noqa: E402
    activated_compatible_state,
)
from recursive_horizons.fgc.modified_harmonic_modes import (  # noqa: E402
    exact_comp1_acceleration_root,
)
from recursive_horizons.fgc.modified_harmonic_physical_modes import (  # noqa: E402
    PHYSICAL_REQUIRED_CHART,
    physical_mode_certificate,
    quadratic_mode_chart,
    quotient_determinant,
)


class PhysicalModeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        config = load_config()
        datum = activated_compatible_state(config.flat_fixture)
        root = exact_comp1_acceleration_root(
            datum["state"],
            reference=datum["reference"],
            coordinate_radius=Q(4),
            tilde_normal_factor=Q(4),
            hat_normal_factor=Q(9),
        )
        cls.certificate = physical_mode_certificate(root["solved_state"])

    def test_comp1_has_complete_exact_physical_mode_charts(self) -> None:
        certificate = self.certificate
        self.assertTrue(certificate["all_exact_checks_pass"])
        self.assertTrue(certificate["scalar_factorization_exact"])
        self.assertFalse(certificate["uniform_domain_eigenframe_proven"])
        self.assertFalse(certificate["strong_hyperbolicity_proven"])
        physical = certificate["charts"]["physical_metric"]
        regulator = certificate["charts"]["regulator"]
        self.assertEqual(
            (physical["omitted_row"], physical["free_column"]),
            PHYSICAL_REQUIRED_CHART,
        )
        self.assertEqual(physical["second_order_vector"][-1], (Q(1), Q(0)))
        for chart in (physical, regulator):
            self.assertGreater(chart["quadratic_discriminant"], 0)
            self.assertNotEqual(chart["pivot_unit_norm"], 0)
            self.assertTrue(chart["second_order_residue_zero"])
            self.assertTrue(chart["first_order_residue_zero"])
            self.assertTrue(chart["free_coordinate_normalized"])
            self.assertTrue(
                chart["kernel_dimension_one_at_both_roots_by_unit_minor"]
            )

    def test_regulator_chart_selection_is_deterministic(self) -> None:
        regulator = self.certificate["charts"]["regulator"]
        self.assertEqual(
            (regulator["omitted_row"], regulator["free_column"]), (0, 0)
        )

    def test_malformed_factors_symbols_and_nonunits_fail_closed(self) -> None:
        zero = (((Q(0),),) * 6,) * 6
        with self.assertRaisesRegex(ValueError, "unit required"):
            quadratic_mode_chart(
                zero,
                (Q(-1), Q(0), Q(1)),
                required_chart=PHYSICAL_REQUIRED_CHART,
            )
        with self.assertRaisesRegex(ValueError, "degree exactly two"):
            quadratic_mode_chart(zero, (Q(-1), Q(1)))
        with self.assertRaisesRegex(ValueError, "six-by-six"):
            quadratic_mode_chart(zero[:5], (Q(-1), Q(0), Q(1)))
        with self.assertRaisesRegex(ValueError, "square"):
            quotient_determinant((), b=Q(0), d=Q(-1))


if __name__ == "__main__":
    unittest.main()
