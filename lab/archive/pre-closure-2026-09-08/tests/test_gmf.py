from __future__ import annotations

from math import inf, nan, pi
import unittest

from recursive_horizons.gmf import (
    ClosedDeSitterChild,
    KottlerRegion,
    TimelikeShell,
    child_region_sample,
    classify_parent_round_sphere,
    darmois_dynamic_gate,
    darmois_mass_jump,
    os_boundary_geodesic_residual,
    os_boundary_radius,
    os_dust_density,
    os_dust_ricci_scalar,
    os_endpoint_status,
    os_friedmann_residual,
    os_mass,
    parent_expansion_product_residual,
    shell_conservation_residual,
    shell_effective_potential,
    shell_junction_residual,
    shell_orientation_gate,
    shell_static_gate,
)


class GMFTests(unittest.TestCase):
    def setUp(self) -> None:
        self.parent = KottlerRegion(mass=1.0)
        self.child = ClosedDeSitterChild(radius=10.0)

    def test_same_lambda_nonzero_mass_has_darmois_obstruction(self) -> None:
        same = ClosedDeSitterChild(radius=10.0)
        parent = KottlerRegion(mass=2.0, cosmological_constant=same.cosmological_constant)
        self.assertAlmostEqual(darmois_mass_jump(parent, same, 4.0), 2.0)
        gate = darmois_dynamic_gate(parent, same, 3.0, 4.0)
        self.assertFalse(gate["mass_gate_passes_at_both_radii"])
        self.assertFalse(gate["input_is_trivial_same_geometry_control"])
        self.assertFalse(gate["full_darmois_match_proven"])
        trivial = KottlerRegion(mass=0.0, cosmological_constant=same.cosmological_constant)
        trivial_gate = darmois_dynamic_gate(trivial, same, 3.0, 4.0)
        self.assertTrue(trivial_gate["input_is_trivial_same_geometry_control"])
        self.assertTrue(trivial_gate["mass_gate_passes_at_both_radii"])
        self.assertFalse(trivial_gate["full_darmois_match_proven"])

    def test_unequal_lambda_can_only_balance_one_fixed_radius(self) -> None:
        parent = KottlerRegion(mass=1.0, cosmological_constant=0.0)
        radius = (6.0 / self.child.cosmological_constant) ** (1.0 / 3.0)
        self.assertAlmostEqual(darmois_mass_jump(parent, self.child, radius), 0.0)
        self.assertNotAlmostEqual(darmois_mass_jump(parent, self.child, 1.1 * radius), 0.0)

    def test_oppenheimer_snyder_exact_matching_and_singular_endpoint(self) -> None:
        scale = 10.0
        density = 3.0 / (8.0 * pi * scale**2)
        chi = pi / 3.0
        mass = os_mass(density, scale, chi)
        self.assertGreater(mass, 0.0)
        self.assertAlmostEqual(os_friedmann_residual(density, scale, 0.0), 0.0)
        self.assertAlmostEqual(os_boundary_geodesic_residual(density, scale, 0.0, chi), 0.0)
        self.assertAlmostEqual(
            os_dust_density(density, scale, scale / 100.0) / density,
            1.0e6,
        )
        self.assertAlmostEqual(
            os_dust_ricci_scalar(density, scale, scale / 100.0)
            / os_dust_ricci_scalar(density, scale, scale),
            1.0e6,
        )
        self.assertTrue(os_endpoint_status()["singularity_present_at_a_zero"])

    def test_parent_ef_expansion_regions_and_product_identity(self) -> None:
        horizon = 2.0
        self.assertEqual(classify_parent_round_sphere(self.parent, 3.0), "normal")
        self.assertEqual(classify_parent_round_sphere(self.parent, horizon), "marginal")
        self.assertEqual(classify_parent_round_sphere(self.parent, 1.0), "trapped")
        self.assertAlmostEqual(parent_expansion_product_residual(self.parent, 1.5), 0.0)
        self.assertEqual(child_region_sample(2.0, pi / 2.0, self.child), "anti_trapped")

    def test_timelike_shell_requires_unsquared_orientation_and_is_not_static(self) -> None:
        radius = 4.0
        kappa = ((sqrt_child := (1.0 - self.child.cosmological_constant * radius**2 / 3.0) ** 0.5) - (1.0 - 2.0 / radius) ** 0.5) / radius
        shell = TimelikeShell(surface_density=kappa / (4.0 * pi))
        self.assertAlmostEqual(shell_junction_residual(self.parent, self.child, shell, radius, 0.0), 0.0)
        self.assertAlmostEqual(shell_effective_potential(self.parent, self.child, shell, radius), 0.0)
        accepted = shell_orientation_gate(self.parent, self.child, shell, radius, 0.0)
        self.assertTrue(accepted["accepted"])
        wrong = TimelikeShell(surface_density=shell.surface_density, epsilon_child=-1)
        self.assertFalse(shell_orientation_gate(self.parent, self.child, wrong, radius, 0.0)["accepted"])
        self.assertFalse(
            shell_static_gate(self.parent, self.child, wrong, radius)[
                "orientation_accepted_at_rest"
            ]
        )
        static = shell_static_gate(self.parent, self.child, shell, radius)
        self.assertTrue(static["orientation_accepted_at_rest"])
        self.assertFalse(static["actually_static"])
        self.assertTrue(static["fixed_kappa_during_derivative"])
        self.assertFalse(static["general_shell_stability_claimed"])
        self.assertAlmostEqual(shell_conservation_residual(shell.surface_density, -shell.surface_density, radius, 0.3), 0.0)

    def test_invalid_inputs_are_rejected(self) -> None:
        for invalid in (True, "1.0", None, inf, nan):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    KottlerRegion(mass=invalid)  # type: ignore[arg-type]
        for invalid in (0, 1.0, True, "1"):
            with self.subTest(orientation=invalid):
                with self.assertRaises(ValueError):
                    TimelikeShell(
                        surface_density=1.0,
                        epsilon_parent=invalid,  # type: ignore[arg-type]
                    )
        for chi in (0.0, pi / 2.0, 0.75 * pi):
            with self.subTest(chi=chi):
                with self.assertRaises(ValueError):
                    os_boundary_radius(1.0, chi)
        with self.assertRaises(ValueError):
            darmois_dynamic_gate(self.parent, self.child, 1.0, 1.0)


if __name__ == "__main__":
    unittest.main()
