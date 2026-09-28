from __future__ import annotations

from copy import deepcopy
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc.regular_center import (  # noqa: E402
    LaurentSeries,
    SeriesJet2,
    validate_regular_profile,
)
from recursive_horizons.fgc.sgb1_ctl1_center import (  # noqa: E402
    EVEN_JET_MINIMUM_ERROR_RATIO,
    FIRST_GRID_RADII,
    FORMAL_SERIES_IDENTITY_KEYS,
    REGULAR_EQUATION_ORDER,
    sgbl_elementary_flatness_defect,
    sgbl_empty_minkowski_center_profile,
    sgbl_initial_center_first_grid,
    sgbl_initial_center_series,
    sgbl_nontrivial_regular_center_profile,
    sgbl_pointwise_center_state,
    sgbl_pointwise_regular_equations,
    validate_sgbl_initial_center_profile,
)
from recursive_horizons.fgc.sgb1_ctl1_source import (  # noqa: E402
    REFERENCE_RADIAL_MINIMUM,
    SGBLSourceInputs,
    sgbl_source_residual,
)
from recursive_horizons.fgc.spherical_reduction import residuals  # noqa: E402


class SGBLInitialCenterTests(unittest.TestCase):
    def test_empty_minkowski_buffer_is_an_exact_initial_centre(self) -> None:
        profile = sgbl_empty_minkowski_center_profile()
        parsed = validate_sgbl_initial_center_profile(profile)
        self.assertTrue(parsed["elementary_flatness_value_dt_dtt_valid"])
        self.assertFalse(parsed["evolving_center_after_matter_arrives"])
        series = sgbl_initial_center_series(profile)
        self.assertEqual(series["center_limits"], (Q(0),) * 6)
        self.assertEqual(series["curvature_center"]["Ricci_scalar"], Q(0))
        self.assertEqual(series["curvature_center"]["Gauss_Bonnet"], Q(0))
        self.assertTrue(series["certified_negative_powers_absent"])
        self.assertTrue(series["initial_empty_center"])
        self.assertFalse(series["evolving_center_after_matter_arrives"])
        self.assertFalse(series["nonclaims"]["evolving_center_after_matter_arrives_qualified"])
        self.assertFalse(series["nonclaims"]["SGBL_branch_owned_and_healthy"])
        grid = sgbl_initial_center_first_grid(profile, formal_series=series)
        self.assertTrue(grid["all_regular_equation_controls_passed"])
        self.assertFalse(grid["uses_source_annulus_minimum_one_half"])
        self.assertTrue(grid["not_a_PDE_or_time_method_convergence_claim"])
        for record in grid["records"]:
            self.assertTrue(all(value == 0 for value in record["values"]))

    def test_qr_profile_validator_refuses_the_sgbl_action(self) -> None:
        profile = sgbl_empty_minkowski_center_profile()
        with self.assertRaisesRegex(ValueError, "FGC-QR"):
            validate_regular_profile(profile)
        source = Path(
            __import__("recursive_horizons.fgc.sgb1_ctl1_center", fromlist=["dummy"]).__file__
        ).read_text()
        self.assertNotIn("regular_center_series_certificate", source)
        self.assertNotIn("validate_regular_profile", source)
        self.assertIn("LaurentSeries", source)
        self.assertIn("SeriesJet2", source)

    def test_nontrivial_regular_first_grid_converges(self) -> None:
        profile = sgbl_nontrivial_regular_center_profile()
        series = sgbl_initial_center_series(profile)
        self.assertTrue(any(limit != 0 for limit in series["center_limits"]))
        self.assertFalse(series["initial_empty_center"])
        grid = sgbl_initial_center_first_grid(profile, formal_series=series)
        self.assertEqual(grid["radii"], FIRST_GRID_RADII)
        self.assertEqual(grid["minimum_required_error_ratio"], EVEN_JET_MINIMUM_ERROR_RATIO)
        self.assertTrue(grid["all_regular_equation_controls_passed"])
        self.assertLess(FIRST_GRID_RADII[-1], REFERENCE_RADIAL_MINIMUM)
        self.assertFalse(grid["evolving_center_after_matter_arrives"])
        for record in grid["records"]:
            self.assertTrue(record["at_least_second_order_control_passed"])
            for ratio in record["coarse_to_fine_error_ratios"]:
                if ratio is not None:
                    self.assertGreaterEqual(ratio, EVEN_JET_MINIMUM_ERROR_RATIO)

    def test_elementary_flatness_and_odd_scalar_fail_closed(self) -> None:
        conical = deepcopy(sgbl_nontrivial_regular_center_profile())
        conical["fields"]["A"]["value"][0] = Q(11, 10)
        self.assertEqual(
            sgbl_elementary_flatness_defect(conical)["Ricci_scalar_r_minus_2_from_conical_defect"],
            -Q(42, 121),
        )
        with self.assertRaisesRegex(ValueError, "elementary flatness"):
            validate_sgbl_initial_center_profile(conical)

        odd = deepcopy(sgbl_nontrivial_regular_center_profile())
        odd["fields"]["phi"]["value"][1] = Q(1, 17)
        self.assertEqual(
            sgbl_elementary_flatness_defect(odd)["box_phi_r_minus_1_from_odd_scalar_term"],
            Q(2, 17),
        )
        with self.assertRaisesRegex(ValueError, "even centre parity"):
            validate_sgbl_initial_center_profile(odd)

        time_defect = deepcopy(sgbl_nontrivial_regular_center_profile())
        time_defect["fields"]["A"]["dt"][0] = Q(1, 39)
        with self.assertRaisesRegex(ValueError, "elementary flatness"):
            validate_sgbl_initial_center_profile(time_defect)

    def test_pointwise_evaluation_refuses_the_origin_and_source_annulus_is_not_a_center(self) -> None:
        profile = sgbl_empty_minkowski_center_profile()
        with self.assertRaisesRegex(ValueError, "r>0"):
            sgbl_pointwise_center_state(profile, coordinate_radius=0)
        with self.assertRaisesRegex(ValueError, "r>0"):
            sgbl_pointwise_regular_equations(profile, coordinate_radius=0)
        self.assertEqual(REFERENCE_RADIAL_MINIMUM, Q(1, 2))
        self.assertLess(FIRST_GRID_RADII[0], REFERENCE_RADIAL_MINIMUM)
        empty_at_first_grid = sgbl_pointwise_regular_equations(
            profile, coordinate_radius=FIRST_GRID_RADII[0]
        )
        self.assertEqual(empty_at_first_grid, (Q(0),) * 6)
        point = SGBLSourceInputs(
            coordinate_radius=FIRST_GRID_RADII[0],
            alpha=(1, 0, 0, 0, 0),
            shift=(0, 0, 0, 0, 0),
            radial_metric=(1, 0, 0, 0, 0),
            areal_radius=(FIRST_GRID_RADII[0], 0, 1, 0, 0),
            phi=(0, 0, 0, 0, 0),
            chi=(0, 0, 0, 0, 0),
            planck_mass=2,
            scalar_mass=3,
            quartic_coupling=Q(1, 2),
            alpha_gb=-Q(1, 4),
        )
        with self.assertRaisesRegex(ValueError, "strictly inside the annulus"):
            sgbl_source_residual(point)
        state = sgbl_pointwise_center_state(profile, coordinate_radius=Q(1, 16))
        self.assertEqual(residuals(state)["R"], Q(0))
        self.assertEqual(type(LaurentSeries.constant(1)), LaurentSeries)
        self.assertEqual(type(SeriesJet2.constant(1)), SeriesJet2)

    def test_initial_empty_centre_is_not_an_evolving_centre_claim(self) -> None:
        empty = sgbl_initial_center_series(sgbl_empty_minkowski_center_profile())
        nontrivial = sgbl_initial_center_series(sgbl_nontrivial_regular_center_profile())
        self.assertNotEqual(empty["profile_id"], nontrivial["profile_id"])
        self.assertTrue(empty["initial_empty_center"])
        self.assertFalse(nontrivial["initial_empty_center"])
        self.assertFalse(empty["evolving_center_after_matter_arrives"])
        self.assertFalse(nontrivial["evolving_center_after_matter_arrives"])

    def test_emptiness_is_derived_from_coefficients_not_the_profile_name(self) -> None:
        renamed_nonempty = deepcopy(sgbl_nontrivial_regular_center_profile())
        renamed_nonempty["profile_id"] = "sgbl_empty_not_actually_empty"
        nonempty = sgbl_initial_center_series(renamed_nonempty)
        self.assertEqual(nonempty["profile_id"], "sgbl_empty_not_actually_empty")
        self.assertFalse(nonempty["initial_empty_center"])
        self.assertTrue(any(limit != 0 for limit in nonempty["center_limits"]))

        renamed_empty = deepcopy(sgbl_empty_minkowski_center_profile())
        renamed_empty["profile_id"] = "definitely_not_empty_by_name"
        empty = sgbl_initial_center_series(renamed_empty)
        self.assertTrue(empty["initial_empty_center"])
        self.assertEqual(empty["center_limits"], (Q(0),) * 6)

    def test_first_grid_refuses_incomplete_foreign_and_altered_formal_caches(self) -> None:
        profile = sgbl_nontrivial_regular_center_profile()
        computed = sgbl_initial_center_series(profile)
        self.assertEqual(len(computed["center_limits"]), len(REGULAR_EQUATION_ORDER))
        self.assertEqual(
            set(FORMAL_SERIES_IDENTITY_KEYS),
            {
                "profile_id",
                "regular_equation_order",
                "center_limits",
                "regular_equation_series",
                "curvature_center",
            },
        )

        with self.assertRaisesRegex(ValueError, "incomplete"):
            sgbl_initial_center_first_grid(
                profile,
                formal_series={
                    "profile_id": profile["profile_id"],
                    "center_limits": (),
                },
            )
        with self.assertRaisesRegex(ValueError, "complete 6-equation exact vector"):
            sgbl_initial_center_first_grid(
                profile,
                formal_series={
                    **computed,
                    "center_limits": computed["center_limits"][:5],
                },
            )
        with self.assertRaisesRegex(ValueError, "different profile"):
            sgbl_initial_center_first_grid(
                profile,
                formal_series={**computed, "profile_id": "foreign-cache"},
            )
        with self.assertRaisesRegex(ValueError, "center_limits differ"):
            sgbl_initial_center_first_grid(
                profile,
                formal_series={**computed, "center_limits": (Q(0),) * 6},
            )
        accepted = sgbl_initial_center_first_grid(profile, formal_series=computed)
        self.assertEqual(len(accepted["records"]), 6)
        self.assertTrue(accepted["all_regular_equation_controls_passed"])


if __name__ == "__main__":
    unittest.main()
