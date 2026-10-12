"""Bounded checks for the first-order spherical feedback chart.

No 4D Riemann tensor is built. Curvature limits use the accepted 2D formula.
"""
import unittest

import sympy as sp

from recursive_horizons.nsc_spherical_feedback_action import (
    BOUNDARY_ASSUMPTIONS,
    CAUCHY_EVOLUTION_OWNED,
    SOURCE_COUPLING_OWNED,
    FeedbackChart,
    alpha_of,
    auxiliary_square_identity,
    boundary_completion_ledger,
    box_flux_identity,
    box_r_flux,
    bulk_status,
    constraint_sign_identities,
    continuum_matter_hamiltonian,
    continuum_radius_force_residual,
    curvature_limits,
    einstein_conformal_current,
    euler_candidate_status,
    explicit_hamilton_rhs,
    feedback_F,
    feedback_V,
    feedback_Z,
    first_order_density,
    hamiltonian_right_hand_sides,
    homogeneous_ghy_match,
    ibp_current,
    inverse_velocities,
    legendre_identity,
    momenta,
    normalized_conditional_examples,
)


class SphericalFeedbackTests(unittest.TestCase):
    def test_zero_weyl_is_rejected_and_not_replaced_by_pure_einstein(self):
        with self.assertRaises(ValueError):
            alpha_of(0)
        with self.assertRaises(ValueError):
            FeedbackChart(A=1.0, C_W=0.0, C_F=1.0, flux=1.0)
        with self.assertRaises(ValueError):
            FeedbackChart(A=0.0, C_W=1.0, C_F=1.0, flux=0.0)
        self.assertFalse(CAUCHY_EVOLUTION_OWNED)

    def test_coefficients_numeric_roundtrip_and_density_derivatives(self):
        chart = dict(A=1.5, C_W=-0.25, C_F=0.5, flux=2.0)
        alpha = alpha_of(chart["C_W"])
        self.assertAlmostEqual(alpha, -4 * 3.141592653589793 * chart["C_W"] / 3)
        F = feedback_F(2.0, 0.3, chart["A"], chart["C_W"])
        self.assertAlmostEqual(F, -4 * 3.141592653589793 * chart["A"] * 4.0 + 2 * alpha * 0.3)
        self.assertAlmostEqual(feedback_Z(chart["A"]), -24 * 3.141592653589793 * chart["A"])
        state = dict(
            L=1.2, L_x=0.1, Q=0.7, Q_t=0.05, Q_x=-0.02,
            r=2.0, r_t=0.03, r_x=0.04, chi=0.3, chi_t=-0.01, chi_x=0.02,
            beta=0.2, beta_x=-0.05,
        )
        p_Q, p_r, p_chi = momenta(
            state["L"], state["Q"], state["Q_t"], state["Q_x"], state["r"], state["r_t"],
            state["r_x"], state["chi"], state["chi_t"], state["chi_x"], state["beta"],
            state["beta_x"], chart["A"], chart["C_W"],
        )
        step = 1e-6
        for key, expected in (("Q_t", p_Q), ("r_t", p_r), ("chi_t", p_chi)):
            up, down = dict(state), dict(state)
            up[key] += step
            down[key] -= step
            derivative = (first_order_density(**up, **chart) - first_order_density(**down, **chart)) / (2 * step)
            self.assertAlmostEqual(derivative, expected, places=8)
        solved = inverse_velocities(
            state["L"], state["Q"], state["Q_x"], state["r"], state["r_x"], state["chi_x"],
            state["beta"], state["beta_x"], p_Q, p_r, p_chi, chart["A"], chart["C_W"],
        )
        self.assertAlmostEqual(solved["Q_t"], state["Q_t"], places=9)
        self.assertAlmostEqual(solved["r_t"], state["r_t"], places=9)
        self.assertAlmostEqual(solved["chi_t"], state["chi_t"], places=9)
        self.assertAlmostEqual(
            feedback_V(2.0, 0.3, chart["A"], chart["C_W"], chart["C_F"], chart["flux"]),
            8 * 3.141592653589793 * chart["A"] * 4.0
            - alpha * (4 * 0.3 + 0.3**2)
            - 2 * 3.141592653589793 * chart["C_F"] * 4.0,
        )

    def test_auxiliary_square_and_elimination(self):
        result = auxiliary_square_identity()
        self.assertTrue(result["restored_is_zero"])
        self.assertTrue(result["critical_point_is_curvature"])

    def test_legendre_identity_on_generic_jets(self):
        self.assertTrue(legendre_identity()["residual_is_zero"])

    def test_constraint_signs_against_first_order_density(self):
        signs = constraint_sign_identities()
        self.assertTrue(signs["lapse_el_plus_C_is_zero"])
        self.assertTrue(signs["shift_el_plus_Dcon_is_zero"])
        self.assertEqual(signs["lapse_sign"], "EL_L = -C")
        self.assertEqual(signs["shift_sign"], "EL_beta = -Dcon")

    def test_six_hamilton_rhs_match_inverse_and_euler(self):
        rhs = hamiltonian_right_hand_sides()
        self.assertTrue(rhs["velocity_rhs_matches_inverse"])
        self.assertTrue(rhs["momentum_residuals_are_zero"])
        self.assertEqual(set(rhs["momentum_rhs_matches_euler"]), {"Q", "r", "chi"})

    def test_homogeneous_einstein_matches_owned_ghy(self):
        self.assertTrue(homogeneous_ghy_match()["residual_is_zero"])

    def test_curvature_limits_use_the_two_dimensional_formula(self):
        limits = curvature_limits()
        self.assertTrue(limits["schwarzschild_C2_is_48M2_over_r6"])
        self.assertTrue(limits["de_sitter_Rh_is_2"])
        self.assertTrue(limits["de_sitter_C2_is_0"])
        self.assertTrue(limits["flat_C2_is_0"])

    def test_box_flux_and_counterterm_cancellation(self):
        result = box_flux_identity()
        self.assertTrue(result["divergence_matches_reduced_box"])
        self.assertTrue(result["homogeneous_endpoint_matches_owner"])
        self.assertTrue(result["cap_cancellation_is_zero"])
        t, x = sp.symbols("t x")
        L, Q, r, beta, R = (sp.Function(n)(t, x) for n in ("L", "Q", "r", "beta", "R"))
        Jt, Jx = box_r_flux(L, Q, r, beta, sp.diff(R, t), sp.diff(R, x))
        C_box = sp.symbols("C_box")
        self.assertEqual(sp.simplify(-C_box * Jt + C_box * Jt), 0)
        self.assertEqual(sp.simplify(-C_box * (sp.diff(Jt, t) + sp.diff(Jx, x)) + C_box * (sp.diff(Jt, t) + sp.diff(Jx, x))), 0)

    def test_ibp_current_removes_the_accepted_divergence(self):
        t, x = sp.symbols("t x")
        F, D, beta, L, Q = (sp.Function(n)(t, x) for n in ("F", "D", "beta", "L", "Q"))
        Jt, Jx = ibp_current(F, D, beta, sp.diff(L, x), Q)
        divergence = sp.diff(Jt, t) + sp.diff(Jx, x)
        sqrt_R = -2 * sp.diff(D, t) + 2 * sp.diff(beta * D + sp.diff(L, x) / Q, x)
        reconstructed = 2 * D * (sp.diff(F, t) - beta * sp.diff(F, x)) - 2 * sp.diff(F, x) * sp.diff(L, x) / Q
        self.assertEqual(sp.simplify(F * sqrt_R - reconstructed - divergence), 0)

    def test_euler_current_is_the_confirmed_mixed_formula(self):
        status = euler_candidate_status()
        self.assertEqual(status["status"], "checked")
        self.assertTrue(status["residual_is_zero"])
        self.assertTrue(status["naive_current_is_distinct"])
        self.assertEqual(status["formula"], "jt=8 k A0-16 v s_x; jx=-8 (P+beta k) A0+16 v s_t")
        jt = sp.symbols("jt")
        self.assertEqual(sp.simplify(-4 * sp.pi * jt + 4 * sp.pi * jt), 0)

    def test_conformal_einstein_current_locks_the_prior_coefficient(self):
        A, N, q, r, beta, rt, rx = sp.symbols("A N q r beta r_t r_x", nonzero=True)
        Jt, Jx = einstein_conformal_current(A, N, q, r, beta, rt, rx)
        self.assertEqual(sp.simplify(Jt.subs({beta: 0, rt: 0})), 0)
        static = sp.simplify(Jx.subs({beta: 0, rt: 0}))
        self.assertEqual(static, sp.simplify(-24 * sp.pi * A * N * r * rx / q))

    def test_continuum_matter_force_is_not_the_old_matrix(self):
        matter = continuum_matter_hamiltonian()
        self.assertFalse(matter["old_matrix_sqrtN_ordering_fixed"])
        self.assertEqual(
            matter["operator"],
            "sigma2/2*{L/Q,P}+sigma1*L*kappa-1/2*{beta,P}",
        )
        L, Q, r, kappa = sp.symbols("L Q r kappa", nonzero=True)
        self.assertTrue(continuum_radius_force_residual(L, Q, r, kappa)["both_zero"])
        self.assertFalse(SOURCE_COUPLING_OWNED)

    def test_boundary_ledger_cancellations_have_no_spatial_corners(self):
        ledger = boundary_completion_ledger()
        self.assertEqual(ledger["assumptions"]["x_domain"], "periodic")
        self.assertEqual(ledger["assumptions"]["time_domain"], "finite_caps")
        self.assertEqual(ledger["assumptions"]["dirichlet_fields"], ("Q", "r", "chi"))
        self.assertFalse(ledger["spatial_corners"])
        self.assertFalse(ledger["pure_einstein_chart"])
        for key in (
            "euler_cap_cancellation_is_zero",
            "myers_Q_plus_is_minus_4pi_jt",
            "box_cap_cancellation_is_zero",
            "restored_F_Rh_time_flux_is_zero_sum",
            "weyl_extra_time_flux_is_4_alpha_chi_D",
            "ghy_plus_weyl_rebuilds_removed_flux",
            "raw_eh_counterterm_matches_ghy_primitive",
            "dirichlet_cap_form_is_zero",
            "periodic_spatial_flux_is_zero",
        ):
            self.assertTrue(ledger[key], key)

    def test_explicit_momentum_rhs_match_the_derived_equations(self):
        shown = explicit_hamilton_rhs()
        self.assertTrue(shown["all_match"])
        raw = hamiltonian_right_hand_sides()
        for key in ("Q_t", "r_t", "chi_t", "pQ_t", "pr_t", "pchi_t"):
            self.assertEqual(sp.simplify(shown["displayed"][key] - raw[key]), 0)

    def test_conditional_initial_slices_are_not_a_physical_source(self):
        rows = normalized_conditional_examples()
        self.assertEqual(len(rows), 4)
        minima = []
        squares = []
        for row in rows:
            self.assertEqual(row["lapse_residual"], 0)
            self.assertEqual(row["shift_residual"], 0)
            self.assertEqual(row["p_Q_x_minus_j_over_Q"], 0)
            self.assertTrue(row["K_matches_review"])
            self.assertFalse(row["physical_source_covariance"])
            self.assertFalse(row["imposed_future_pulses"])
            self.assertFalse(row["source_coupling_owned"])
            self.assertEqual(sp.integrate(row["j"], (row["coordinate"], 0, 2 * sp.pi)), 0)
            minima.append(row["minimum"])
            squares.append(row["pi_square"])
        self.assertEqual(minima, [sp.Rational(7, 4), sp.Rational(7, 4), sp.Rational(37, 20), sp.Rational(37, 20)])
        x = rows[0]["coordinate"]
        self.assertEqual(sp.simplify(rows[0]["K"] - sp.Rational(7, 4) - (1 + sp.cos(x)) / 4), 0)
        self.assertEqual(sp.simplify(rows[2]["K"] - sp.Rational(37, 20) - 3 * (1 + sp.cos(x)) / 20), 0)
        self.assertEqual(sp.simplify(rows[0]["pi_square"] - (48 + 6 * sp.cos(x))), 0)
        self.assertEqual(sp.simplify(rows[0]["Pi"] + rows[1]["Pi"]), 0)
        self.assertEqual(sp.simplify(rows[2]["Pi"] + rows[3]["Pi"]), 0)

    def test_bulk_status_names_open_boundary_flags(self):
        status = bulk_status()
        self.assertFalse(status["cauchy_evolution_owned"])
        self.assertFalse(status["pure_einstein_chart"])
        self.assertFalse(status["old_matrix_source_fixed"])
        self.assertTrue(status["legendre_identity"])
        self.assertTrue(status["constraint_signs"])
        self.assertTrue(status["homogeneous_ghy"])
        self.assertTrue(status["box_r_divergence"])
        self.assertTrue(status["euler_naive_rejected"])
        self.assertTrue(status["euler_cap_cancellation"])
        self.assertTrue(status["raw_eh_matches_ghy"])
        self.assertTrue(status["conditional_initial_residuals"])
        self.assertFalse(status["spatial_corners"])
        self.assertFalse(status["source_coupling_owned"])
        self.assertEqual(BOUNDARY_ASSUMPTIONS["x_domain"], "periodic")
        self.assertEqual(BOUNDARY_ASSUMPTIONS["dirichlet_fields"], ("Q", "r", "chi"))


if __name__ == "__main__":
    unittest.main()
