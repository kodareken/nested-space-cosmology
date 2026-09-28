#!/usr/bin/env python3
"""Complete the charged neck gate with positive compact-level CTP response."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from multiprocessing import get_all_start_methods, get_context
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.nsc_compact_ctp_neck import (
    CompactCTPConfig,
    _adiabatic_bloch_function,
    compact_ctp_source,
)


OUTPUT = ROOT / "results/development/charged-compact-ctp-completion.json"
INPUTS = (
    "results/development/charged-ctp-neck-source.json",
    "results/development/scale-binding.json",
    "results/development/horizon-source.json",
    "results/development/canonical-spectral-bridge.json",
    "results/development/compact-interaction.json",
)
SOURCES = (
    "src/recursive_horizons/nsc_compact_ctp_neck.py",
    "scripts/derive_nsc_compact_ctp_completion.py",
    "docs/nsc-compact-ctp-completion.md",
)


def hashes(paths):
    return {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in paths}


def compare(expected, actual, path="$", atol=2e-6, rtol=2e-5):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys():
            raise AssertionError(f"keys differ at {path}")
        for key in expected:
            compare(expected[key], actual[key], f"{path}/{key}", atol, rtol)
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise AssertionError(f"list differs at {path}")
        for index, (left, right) in enumerate(zip(expected, actual)):
            compare(left, right, f"{path}/{index}", atol, rtol)
    elif isinstance(expected, float):
        if isinstance(actual, bool) or not isinstance(actual, (int, float)):
            raise AssertionError(f"numeric type differs at {path}")
        if abs(expected-actual) > atol+rtol*abs(expected):
            raise AssertionError(f"number differs at {path}")
    elif type(expected) is not type(actual) or expected != actual:
        raise AssertionError(f"value differs at {path}")


def _run(spec, scale, horizon):
    config = CompactCTPConfig(
        magnetic_flux=scale["magnetic_flux"],
        omega=scale["omega"],
        cutoff=scale["omega"],
        horizon_rho=horizon["reused_horizon_rho"],
        surface_gravity=horizon["reused_surface_gravity"],
        compact_levels=tuple(spec["compact_levels"]),
        angular_levels=spec["angular_levels"],
        frequency_points=spec["frequency_points"],
        frequency_max=spec["frequency_max"],
        evolution_steps=spec["evolution_steps"],
        scattering_tolerance=spec["scattering_tolerance"],
    )
    return compact_ctp_source(config)


def calculate():
    retained, scale_record, horizon_record, allocation, interaction = (
        json.loads((ROOT/path).read_text()) for path in INPUTS
    )
    scale = scale_record["development_branch"]
    horizon = horizon_record["benchmark_Unruh_per_abs_q"]
    if scale["magnetic_flux"] != 4 or scale["coefficients"]["V_full_relational"] != 0:
        raise ValueError("locked q=4 scale and relational vacuum required")
    masses = allocation["allocation"]["canonical_fields"]["mass_formula"]
    if masses != "m_n=n*pi/ell":
        raise ValueError("unexpected compact mass convention")
    specifications = [
        {
            "name": "base",
            "compact_levels": [1, 2],
            "angular_levels": 9,
            "frequency_points": 8,
            "frequency_max": 40.0,
            "evolution_steps": 6000,
            "scattering_tolerance": 1e-9,
        },
        {
            "name": "angular_coarse",
            "compact_levels": [1, 2],
            "angular_levels": 7,
            "frequency_points": 8,
            "frequency_max": 40.0,
            "evolution_steps": 6000,
            "scattering_tolerance": 1e-9,
        },
        {
            "name": "frequency_extent_coarse",
            "compact_levels": [1, 2],
            "angular_levels": 9,
            "frequency_points": 8,
            "frequency_max": 32.0,
            "evolution_steps": 6000,
            "scattering_tolerance": 1e-9,
        },
        {
            "name": "frequency_time_refinement",
            "compact_levels": [1, 2],
            "angular_levels": 7,
            "frequency_points": 10,
            "frequency_max": 40.0,
            "evolution_steps": 8000,
            "scattering_tolerance": 5e-10,
        },
    ]
    # Build the symbolic superadiabatic projector once; forked workers share it.
    _adiabatic_bloch_function()
    context = get_context("fork") if "fork" in get_all_start_methods() else get_context()
    with ProcessPoolExecutor(max_workers=4, mp_context=context) as pool:
        futures = [pool.submit(_run, spec, scale, horizon) for spec in specifications]
        runs = [future.result() for future in futures]
    named = {spec["name"]: run for spec, run in zip(specifications, runs)}
    for name, run in named.items():
        checks = run["checks"]
        if checks["minimum_initial_covariance_eigenvalue"] < -3e-12:
            raise AssertionError(f"{name} covariance lost positivity")
        if checks["maximum_initial_covariance_eigenvalue"] > 1+3e-12:
            raise AssertionError(f"{name} covariance exceeded one")
        if checks["maximum_scattering_current_defect"] > 2e-9:
            raise AssertionError(f"{name} current normalization failed")
    base = named["base"]
    old = retained["runs"]["base"]["totals"]
    fields = ("rho", "T01", "p_parallel", "p_sphere",
              "radial_null_plus", "radial_null_minus", "parent_Killing_power")
    total = {field: old[field]+base["totals"][field] for field in fields}
    total["both_radial_null_components_negative"] = (
        total["radial_null_plus"] < 0 and total["radial_null_minus"] < 0
    )
    if not total["both_radial_null_components_negative"]:
        raise AssertionError("completed compact channel changed the null-sign gate")
    compact_changes = {
        "angular_7_to_9": {
            field: base["totals"][field]-named["angular_coarse"]["totals"][field]
            for field in fields
        },
        "frequency_extent_32_to_40": {
            field: base["totals"][field]-named["frequency_extent_coarse"]["totals"][field]
            for field in fields
        },
        "quadrature_time_at_angular_7": {
            field: (
                named["frequency_time_refinement"]["totals"][field]
                - named["angular_coarse"]["totals"][field]
            )
            for field in fields
        },
    }
    pressure = (total["p_parallel"]+2*total["p_sphere"])/3
    background = {
        "rho_neck": total["rho"],
        "isotropic_pressure_projection": pressure,
        "w_isotropic_at_neck": pressure/total["rho"],
        "anisotropic_stress_p_parallel_minus_p_sphere": (
            total["p_parallel"]-total["p_sphere"]
        ),
        "momentum_density": total["T01"],
        "homogeneous_DE_pressure_candidate_available": True,
        "cosmological_density_rate_or_H_of_z_derived": False,
        "remaining_background_inputs": [
            "covariant conservation on the evolved geometry",
            "proper-volume and clock projection",
            "deposition/component assignment",
        ],
    }
    return {
        "schema": "NSC-CHARGED-COMPACT-CTP-COMPLETION-v1",
        "status": (
            "PASS: positive compact levels complete the declared free Gaussian "
            "charged angular+compact CTP neck source and preserve both negative null signs"
        ),
        "source_hashes": hashes(SOURCES),
        "input_hashes": hashes(INPUTS),
        "locked_inputs": {
            "q": scale["magnetic_flux"],
            "Omega": scale["omega"],
            "zeta": scale["zeta"],
            "V_full": scale["coefficients"]["V_full_relational"],
            "compact_mass_formula": "m_j=j*pi/(2*L_star)",
            "positive_levels_below_cutoff": [1, 2],
            "scale_or_flux_refitted": False,
        },
        "renormalized_CTP_allocation": {
            "local_owner": "locked Wilsonian compact response already present once in charged-ctp-neck-source",
            "nonlocal_owner": "fourth-order superadiabatically subtracted canonical compact covariance",
            "adiabatic_projector": "Bloch recursion through derivative order four",
            "instantaneous_vacuum_subtraction_used": False,
            "massless_limit_matches_existing_E2_E4_and_P2_P4": True,
            "positive_level_copy_count": 2,
            "interacting_compact_vertex_included": False,
            "interaction_reason": interaction["normalization"]["s_T_meaning"],
        },
        "runs": named,
        "refinement_changes": compact_changes,
        "completed_tensor": total,
        "background_projection": background,
        "gate": {
            "both_neck_null_components_negative": total["both_radial_null_components_negative"],
            "charged_angular_plus_positive_compact_CTP_complete_in_declared_free_Gaussian_realization": True,
            "accepted_for_homogeneous_background_projection": True,
            "H_of_z_or_likelihood_authorized": False,
            "next_owner": (
                "covariant conservation/proper-volume background projection of this same completed tensor"
            ),
        },
        "scope": {
            "old_charged_angular_record_recomputed": False,
            "old_MMP_or_historical_generator_rerun": False,
            "new_metric_equation_or_ansatz": False,
            "new_compensator": False,
            "q_Omega_or_zeta_changed": False,
            "dark_fraction_target_used": False,
            "angular_and_compact_infinite_remainder_rigorous_bound": False,
            "measured_refinements_recorded": True,
            "compact_interactions_or_torsion_enabled": False,
        },
        "comparison": {
            "fields": "all",
            "exact": "keys, types, strings and source/input hashes",
            "float_atol": 2e-6,
            "float_rtol": 2e-5,
            "exceptions": [],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = calculate()
    if args.check:
        compare(json.loads(OUTPUT.read_text()), result)
        print("positive compact CTP completion reproduced from authenticated inputs")
    elif args.output:
        if args.output.exists():
            raise FileExistsError("refusing to overwrite recorded evidence")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
