#!/usr/bin/env python3
"""Regenerate the bounded GMF-1 spherical-matching artifact."""

from __future__ import annotations

import argparse
import json
from math import pi
from pathlib import Path
import sys

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.gmf import (
    ClosedDeSitterChild,
    KottlerRegion,
    TimelikeShell,
    child_region_sample,
    classify_parent_round_sphere,
    darmois_dynamic_gate,
    darmois_mass_jump,
    os_boundary_geodesic_residual,
    os_dust_density,
    os_dust_ricci_scalar,
    os_endpoint_status,
    os_friedmann_residual,
    os_mass,
    parent_ef_null_expansions,
    parent_expansion_product_residual,
    shell_conservation_residual,
    shell_orientation_gate,
    shell_static_gate,
)  # noqa: E402


def record() -> dict[str, object]:
    parent = KottlerRegion(mass=1.0)
    child = ClosedDeSitterChild(radius=10.0)
    same_lambda_parent = KottlerRegion(mass=1.0, cosmological_constant=child.cosmological_constant)
    radius = 4.0
    f_child = 1.0 - child.cosmological_constant * radius**2 / 3.0
    f_parent = 1.0 - 2.0 * parent.mass / radius
    shell = TimelikeShell(surface_density=(f_child**0.5 - f_parent**0.5) / (4.0 * pi * radius))
    scale = 10.0
    density = 3.0 / (8.0 * pi * scale**2)
    chi = pi / 3.0
    unequal_radius = (6.0 * parent.mass / child.cosmological_constant) ** (1.0 / 3.0)
    plus, minus = parent_ef_null_expansions(parent, 1.5)
    return {
        "schema_version": 1,
        "project_version": "0.11.0",
        "artifact": "GMF-1_spherical_matching_and_flux_closure_audit",
        "model_id": "GMF-1",
        "classification": "spherical_matching_obstruction_and_local_shell_diagnostic",
        "global_solution_constructed": False,
        "child_topology_proven": False,
        "derived_exchange_law": False,
        "thin_shell_distributional": True,
        "scope": "Exact spherical GR matching identities and a prescribed local Israel shell; not a global parent-to-child spacetime, smooth transition action, or dark-energy derivation.",
        "darmois_mass_gate": {
            "same_lambda_parent_mass": same_lambda_parent.mass,
            "same_lambda_mass_jump_at_R_4": darmois_mass_jump(same_lambda_parent, child, radius),
            "same_lambda_nonzero_parent_mass_direct_match_allowed": False,
            "unequal_lambda_balance_radius": unequal_radius,
            "unequal_lambda_jump_at_balance_radius": darmois_mass_jump(parent, child, unequal_radius),
            "unequal_lambda_jump_at_1p1_balance_radius": darmois_mass_jump(parent, child, 1.1 * unequal_radius),
            "two_radius_gate": darmois_dynamic_gate(parent, child, unequal_radius, 1.1 * unequal_radius),
        },
        "oppenheimer_snyder": {
            "scale_factor": scale,
            "density": density,
            "chi_boundary": chi,
            "finite_parent_mass": os_mass(density, scale, chi),
            "friedmann_residual": os_friedmann_residual(density, scale, 0.0),
            "boundary_geodesic_residual": os_boundary_geodesic_residual(density, scale, 0.0, chi),
            "density_at_a_over_100": os_dust_density(
                density,
                scale,
                scale / 100.0,
            ),
            "density_growth_a_to_a_over_100": os_dust_density(
                density,
                scale,
                scale / 100.0,
            )
            / density,
            "ricci_scalar_at_a": os_dust_ricci_scalar(density, scale, scale),
            "ricci_scalar_at_a_over_100": os_dust_ricci_scalar(
                density,
                scale,
                scale / 100.0,
            ),
            "endpoint_status": os_endpoint_status(),
        },
        "parent_null_expansions": {
            "normal_R_3": classify_parent_round_sphere(parent, 3.0),
            "marginal_R_2GM": classify_parent_round_sphere(parent, 2.0),
            "trapped_R_1": classify_parent_round_sphere(parent, 1.0),
            "theta_plus_R_1p5": plus,
            "theta_minus_R_1p5": minus,
            "product_identity_residual_R_1p5": parent_expansion_product_residual(parent, 1.5),
            "child_anti_trapped_sample_not_joined": child_region_sample(2.0, pi / 2.0, child),
        },
        "prescribed_timelike_shell": {
            "radius": radius,
            "radial_velocity": 0.0,
            "surface_density": shell.surface_density,
            "orientation_gate": shell_orientation_gate(parent, child, shell, radius, 0.0),
            "constant_tension_conservation_residual": shell_conservation_residual(shell.surface_density, -shell.surface_density, radius, 0.3),
            "static_gate": shell_static_gate(parent, child, shell, radius),
            "q_derived": False,
            "warning": "A thin Israel shell has distributional curvature and is only a local compatibility diagnostic.",
        },
        "closure_status": {
            "surface_action_supplied": False,
            "smooth_transition_action_supplied": False,
            "global_trajectory_supplied": False,
            "global_parent_child_solution": False,
            "bulk_Q_derived": False,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPOSITORY / "results" / "gmf-1-spherical-matching.json",
    )
    args = parser.parse_args()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(
        record(),
        indent=2,
        sort_keys=True,
        allow_nan=False,
    ) + "\n"
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(output)
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
