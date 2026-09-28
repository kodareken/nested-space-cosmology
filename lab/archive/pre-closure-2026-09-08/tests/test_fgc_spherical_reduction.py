from __future__ import annotations

from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import recursive_horizons.fgc.spherical_reduction as reduction_module  # noqa: E402
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    Jet2,
    SphericalState,
    adm_pg_kinematics,
    adm_pg_principal_matrix,
    direct_4d_curvature,
    flat_flrw_curvature_control,
    principal_matrix,
    red1_certificate,
    residuals,
    spherical_reduction_certificate,
    state_from_generalized_adm_pg_fixture,
    warped_2plus2_curvature,
)


def jet(value, dt=0, dr=0, dtt=0, dtr=0, drr=0) -> Jet2:
    return Jet2(Q(value), Q(dt), Q(dr), Q(dtt), Q(dtr), Q(drr))


def generic_state(branch: str = "FGC-QR") -> SphericalState:
    return SphericalState(
        h_tt=jet(-2, 1, 1, 2, 1, -1),
        h_tr=jet(1, 0, 1, 1, 0, 1),
        h_rr=jet(3, 1, 0, 0, 1, 2),
        areal_radius=jet(2, 1, 1, 1, 1, 2),
        phi=jet(1, 1, 2, 3, 1, 2),
        chi=jet(2, 2, 1, 1, 2, 1),
        planck_mass=Q(2),
        beta=Q(-1, 2) if branch == "FGC-QR" else Q(0),
        mu=Q(3),
        g4=Q(1, 2),
        alpha=Q(-1, 4) if branch == "SGB-L" else Q(0),
        eta=Q(1, 2) if branch == "FGC-QR" else Q(0),
        branch=branch,
    )


def schwarzschild_fixture() -> dict[str, object]:
    zero = "0"
    return {
        "fixture_id": "schwarzschild",
        "model_id": "FGC-QR",
        "purpose": "zero_regulator_ricci_flat_schwarzschild_exterior_control",
        "action_parameters": {
            "planck_mass": "2",
            "scalar_mass": "3",
            "quartic_coupling": "1/2",
            "ricci_coupling": "-1/4",
            "linear_gb_coupling": "0",
            "quadratic_gb_coupling": "1/2",
        },
        "state": {
            "alpha": "1",
            "alpha_t": zero,
            "alpha_r": zero,
            "alpha_tt": zero,
            "alpha_tr": zero,
            "alpha_rr": zero,
            "shift": "1/2",
            "shift_t": zero,
            "shift_r": "-1/32",
            "shift_tt": zero,
            "shift_tr": zero,
            "shift_rr": "3/512",
            "lambda": "1",
            "lambda_t": zero,
            "lambda_r": zero,
            "lambda_tt": zero,
            "lambda_tr": zero,
            "lambda_rr": zero,
            "areal_radius": "8",
            "areal_radius_t": zero,
            "areal_radius_r": "1",
            "areal_radius_tt": zero,
            "areal_radius_tr": zero,
            "areal_radius_rr": zero,
            "phi": zero,
            "phi_t": zero,
            "phi_r": zero,
            "phi_tt": zero,
            "phi_tr": zero,
            "phi_rr": zero,
            "chi": zero,
            "chi_t": zero,
            "chi_r": zero,
            "chi_tt": zero,
            "chi_tr": zero,
            "chi_rr": zero,
        },
    }


class FGCRED1Tests(unittest.TestCase):
    def test_jet_arithmetic_and_chain_rule(self) -> None:
        value = jet(2, 3, 5, 7, 11, 13)
        self.assertEqual((value * value).dtr, 2 * 3 * 5 + 2 * 2 * 11)
        self.assertEqual(value.reciprocal().value, Q(1, 2))
        self.assertEqual(value.compose(4, 6, 8).dtt, 6 * 7 + 8 * 9)
        self.assertEqual((value / value).value, 1)

    def test_direct_four_dimensional_and_warped_routes_agree(self) -> None:
        for state in (
            generic_state("GR-0"),
            generic_state("SGB-L"),
            generic_state("FGC-QR"),
        ):
            with self.subTest(branch=state.branch):
                direct_riemann, _, _ = direct_4d_curvature(state)
                self.assertEqual(direct_riemann, warped_2plus2_curvature(state))
                direct = residuals(state)
                warped = residuals(state, use_warped=True)
                for key in (
                    "ricci",
                    "R",
                    "GB",
                    "hessian_phi",
                    "hessian_chi",
                    "metric",
                    "phi",
                    "chi",
                ):
                    self.assertEqual(direct[key], warped[key])

    def test_branch_sources_and_uneliminated_principal_matrix(self) -> None:
        gr = residuals(generic_state("GR-0"))
        self.assertTrue(
            all(value == 0 for row in gr["nonminimal_term"] for value in row)
        )
        self.assertTrue(
            all(value == 0 for row in gr["gb_residual_term"] for value in row)
        )

        constant_phi = replace(generic_state("SGB-L"), phi=jet(Q(3, 2)))
        constant_f = residuals(constant_phi)
        self.assertTrue(
            all(value == 0 for row in constant_f["gb_residual_term"] for value in row)
        )

        activated = generic_state("FGC-QR")
        matrix = principal_matrix(activated)
        self.assertEqual(len(matrix["matrix"]), 6)
        self.assertTrue(all(len(row) == 18 for row in matrix["matrix"]))
        phi_columns = [
            index
            for index, name in enumerate(matrix["column_order"])
            if name.startswith("phi.")
        ]
        metric_columns = [
            index
            for index, name in enumerate(matrix["column_order"])
            if name.split(".")[0] in {"h_tt", "h_tr", "h_rr", "areal_radius"}
        ]
        self.assertTrue(
            any(
                matrix["matrix"][row][column]
                for row in range(4)
                for column in phi_columns
            )
        )
        self.assertTrue(any(matrix["matrix"][4][column] for column in metric_columns))

    def test_exact_schwarzschild_pg_control(self) -> None:
        state = state_from_generalized_adm_pg_fixture(schwarzschild_fixture())
        result = residuals(state)
        self.assertEqual(result["R"], 0)
        self.assertEqual(result["GB"], Q(3, 16384))
        self.assertTrue(all(value == 0 for row in result["metric"] for value in row))
        self.assertEqual(result["phi"], 0)
        self.assertEqual(result["chi"], 0)

        certificate = spherical_reduction_certificate(
            {"fixtures": [schwarzschild_fixture()]}
        )
        control = certificate["fixtures"][0]["schwarzschild_vacuum_control"]
        self.assertTrue(control["metric_residual_zero"])
        self.assertTrue(control["scalar_residuals_zero"])
        self.assertTrue(control["ricci_scalar_zero"])
        self.assertEqual(control["maximum_abs_metric_residual"], "0")
        self.assertEqual(control["absolute_phi_residual"], "0")
        self.assertEqual(control["absolute_chi_residual"], "0")

    def test_flat_flrw_curvature_control(self) -> None:
        # At the frozen point: a=2, H=3/2, H_dot=-2/3, and radial x=5.
        # Therefore a_dot=3, a_ddot=19/6,
        # R_4=6(H_dot+2H^2)=23, and G=24H^2(H_dot+H^2)=171/2.
        state = SphericalState(
            h_tt=jet(-1),
            h_tr=jet(0),
            h_rr=jet(4, 12, 0, Q(92, 3), 0, 0),
            areal_radius=jet(10, 15, 2, Q(95, 6), 3, 0),
            phi=jet(0),
            chi=jet(0),
            planck_mass=Q(2),
            mu=Q(3),
            g4=Q(1, 2),
            branch="GR-0",
        )
        result = residuals(state)
        self.assertEqual(result["R"], 23)
        self.assertEqual(result["GB"], Q(171, 2))
        self.assertTrue(flat_flrw_curvature_control()["exact_match"])

    def test_adm_pg_mapping_metric_speeds_and_chi_principal_factor(self) -> None:
        fixture = schwarzschild_fixture()
        state = state_from_generalized_adm_pg_fixture(fixture)
        self.assertEqual(state.h_tt.value, Q(-3, 4))
        self.assertEqual(state.h_tr.value, Q(1, 2))
        self.assertEqual(state.h_rr.value, Q(1))

        kinematics = adm_pg_kinematics(Q(1), Q(1), Q(1, 2))
        self.assertEqual(kinematics["v_out"], Q(1, 2))
        self.assertEqual(kinematics["v_in"], Q(-3, 2))
        self.assertEqual(kinematics["base_metric_determinant"], -1)

        matrix = adm_pg_principal_matrix(fixture)
        columns = matrix["column_order"]
        chi_columns = [columns.index(f"chi.{slot}") for slot in ("dtt", "dtr", "drr")]
        self.assertEqual(
            tuple(matrix["matrix"][5][index] for index in chi_columns),
            (Q(-1), Q(1), Q(3, 4)),
        )
        self.assertTrue(
            all(
                matrix["matrix"][row][index] == 0
                for row in range(5)
                for index in chi_columns
            )
        )

    def test_fail_closed_validation_and_curvature_mutation(self) -> None:
        with self.assertRaisesRegex(ValueError, "alpha and lambda must be positive"):
            bad = schwarzschild_fixture()
            bad_state = bad["state"]
            assert isinstance(bad_state, dict)
            bad_state["alpha"] = "0"
            state_from_generalized_adm_pg_fixture(bad)
        with self.assertRaisesRegex(ValueError, "Lorentzian"):
            SphericalState(
                h_tt=jet(1),
                h_tr=jet(0),
                h_rr=jet(1),
                areal_radius=jet(1),
                phi=jet(0),
                chi=jet(0),
                branch="GR-0",
            )
        with self.assertRaisesRegex(ValueError, "beta!=0"):
            replace(generic_state("FGC-QR"), beta=Q(0))
        with self.assertRaisesRegex(ValueError, "F must be positive"):
            replace(generic_state("FGC-QR"), planck_mass=Q(1), beta=Q(-2))

        original = warped_2plus2_curvature

        def mutated(state):
            tensor = original(state)
            tensor[0][2][0][2] += 1
            return tensor

        with patch(
            "recursive_horizons.fgc.spherical_reduction.warped_2plus2_curvature",
            side_effect=mutated,
        ):
            certificate = red1_certificate(generic_state("FGC-QR"))
        self.assertNotEqual(
            certificate["direct_warped_max_exact_residuals"]["riemann"], "0"
        )

        original_connection = reduction_module._warped_2plus2_metric_connection

        def mutated_connection(state):
            metric, inverse, connection = original_connection(state)
            connection[0][2][2] += 1
            return metric, inverse, connection

        with patch(
            "recursive_horizons.fgc.spherical_reduction._warped_2plus2_metric_connection",
            side_effect=mutated_connection,
        ):
            certificate = red1_certificate(generic_state("FGC-QR"))
        self.assertNotEqual(
            certificate["direct_warped_max_exact_residuals"]["connection"], "0"
        )
        self.assertNotEqual(
            certificate["direct_warped_max_exact_residuals"]["hessian_phi"], "0"
        )


if __name__ == "__main__":
    unittest.main()
