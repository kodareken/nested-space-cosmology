#!/usr/bin/env python3
"""Evaluate the charged angular plus compact-local CTP source at the neck."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.nsc_charged_ctp_neck import (
    ChargedCTPNeckConfig,
    charged_ctp_neck_source,
)


OUTPUT = ROOT / "results/development/charged-ctp-neck-source.json"
INPUTS = (
    "results/development/scale-binding.json",
    "results/development/horizon-source.json",
    "results/development/adm-neck-source-map.json",
    "results/development/charged-sector.json",
    "results/development/relational-vacuum-normalization.json",
)
SOURCES = (
    "src/recursive_horizons/nsc_charged_ctp_neck.py",
    "scripts/derive_nsc_charged_ctp_neck.py",
    "docs/nsc-charged-ctp-neck-source.md",
)


def hashes(paths):
    return {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in paths}


def compare(expected, actual, path="$", atol=8e-8, rtol=8e-7):
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
        if abs(expected - actual) > atol + rtol * abs(expected):
            raise AssertionError(f"number differs at {path}")
    elif type(expected) is not type(actual) or expected != actual:
        raise AssertionError(f"value differs at {path}")


def _run(spec, scale, horizon):
    config = ChargedCTPNeckConfig(
        magnetic_flux=scale["magnetic_flux"],
        omega=scale["omega"],
        horizon_rho=horizon["reused_horizon_rho"],
        surface_gravity=horizon["reused_surface_gravity"],
        angular_levels=spec["angular_levels"],
        points_per_frequency_interval=spec["points_per_frequency_interval"],
        evolution_steps=spec["evolution_steps"],
        reflection_tolerance=spec["reflection_tolerance"],
    )
    return charged_ctp_neck_source(
        config,
        compact_weyl_coefficient=scale["coefficients"]["C_Weyl"],
    ).to_dict()


def calculate():
    scale_record, horizon_record, neck, charged, vacuum = (
        json.loads((ROOT / path).read_text()) for path in INPUTS
    )
    scale = scale_record["development_branch"]
    horizon = horizon_record["benchmark_Unruh_per_abs_q"]
    if scale["magnetic_flux"] != 4 or vacuum["gate"]["V_full"] != 0:
        raise ValueError("the locked q=4 relational-vacuum branch is required")
    if not any(row["anomaly_free_zero_sector"] for row in charged["parity_and_anomaly_rows"]):
        raise ValueError("the anomaly-free opposite-parity sector is required")
    specifications = [
        {
            "name": "base",
            "angular_levels": 12,
            "points_per_frequency_interval": 12,
            "evolution_steps": 5000,
            "reflection_tolerance": 2e-10,
        },
        {
            "name": "angular_refinement",
            "angular_levels": 16,
            "points_per_frequency_interval": 12,
            "evolution_steps": 5000,
            "reflection_tolerance": 2e-10,
        },
        {
            "name": "frequency_time_refinement",
            "angular_levels": 12,
            "points_per_frequency_interval": 16,
            "evolution_steps": 8000,
            "reflection_tolerance": 1e-10,
        },
    ]
    with ProcessPoolExecutor(max_workers=3) as pool:
        futures = [pool.submit(_run, spec, scale, horizon) for spec in specifications]
        runs = [future.result() for future in futures]
    named = {spec["name"]: run for spec, run in zip(specifications, runs)}
    base = named["base"]
    angular = named["angular_refinement"]
    temporal = named["frequency_time_refinement"]
    for name, run in named.items():
        massive = run["massive_angular_source"]
        if massive["minimum_initial_covariance_eigenvalue"] < -3e-12:
            raise AssertionError(f"{name} covariance lost positivity")
        if massive["maximum_initial_covariance_eigenvalue"] > 1 + 3e-12:
            raise AssertionError(f"{name} covariance exceeded the CAR bound")
        if massive["mode_norm_defect"] > 2e-11:
            raise AssertionError(f"{name} evolution lost norm")
        if massive["maximum_scattering_current_defect"] > 3e-9:
            raise AssertionError(f"{name} scattering current did not close")
        if not run["totals"]["both_radial_null_components_negative"]:
            raise AssertionError(f"{name} retained null-sign result changed")
    fields = ("rho", "T01", "p_parallel", "p_sphere",
              "radial_null_plus", "radial_null_minus", "parent_Killing_power")
    changes = {
        "angular": {
            key: angular["totals"][key] - base["totals"][key] for key in fields
        },
        "frequency_time": {
            key: temporal["totals"][key] - base["totals"][key] for key in fields
        },
    }
    lll = base["lll_source"]["child_frame"]
    return {
        "schema": "NSC-CHARGED-CTP-NECK-SOURCE-v1",
        "status": (
            "charged angular CTP covariance plus locked compact local response "
            "passes the retained null-sign gate; full compact nonlocal CTP gate remains open"
        ),
        "source_hashes": hashes(SOURCES),
        "input_hashes": hashes(INPUTS),
        "locked_inputs": {
            "magnetic_flux": scale["magnetic_flux"],
            "Omega": scale["omega"],
            "zeta": scale["zeta"],
            "V_full": scale["coefficients"]["V_full_relational"],
            "scale_refitted": False,
            "black_universe_neck": neck["geometry"],
        },
        "operator_and_state": {
            "angular_spectrum": "lambda_n=sqrt(n*(n+abs(q))), degeneracy=2*(abs(q)+2*n), n>=1; LLL handled exactly",
            "compact_allocation": "one charged canonical zero field plus the locked Wilsonian compact local response",
            "state": "parent affine occupation plus scale-inherited child occupation at kappa_c=Omega*kappa_p",
            "domain": "actual exterior reflection and horizon frame continued to the stored Bronnikov neck",
            "real_time_owner": "canonical Lorentzian covariance; no Euclidean heat continuation through the horizon",
            "opposite_compact_parities": True,
        },
        "runs": named,
        "refinement_changes": changes,
        "channel_signs": {
            "scale_inherited_LLL_null_plus": lll["radial_null_plus"],
            "scale_inherited_LLL_null_minus": lll["radial_null_minus"],
            "LLL_channel_sign": "positive in both radial null directions",
            "massive_angular_mean_null": base["massive_angular_source"]["source"]["radial_null"],
            "compact_local_mean_null": base["compact_local_source"]["radial_null"],
            "retained_total_null_plus": base["totals"]["radial_null_plus"],
            "retained_total_null_minus": base["totals"]["radial_null_minus"],
        },
        "gate": {
            "retained_both_neck_null_components_negative": base["totals"]["both_radial_null_components_negative"],
            "retained_four_metric_projections_available": base["diagnostics"]["retained_four_metric_projections_returned"],
            "full_charged_angular_plus_compact_CTP_complete": False,
            "accepted_for_background_projection": False,
            "failure_owner": (
                "nonlocal CTP covariance and four ADM variations of the first positive compact Dirac levels; "
                "the locked local compact coefficient is not their full horizon-domain response"
            ),
        },
        "scope": {
            "new_metric_ansatz_or_equation": False,
            "old_MMP_angular_or_historical_generator_rerun": False,
            "new_compensator_added": False,
            "q_Omega_or_zeta_refitted": False,
            "dark_fraction_target_used": False,
            "cosmological_projection_or_likelihood_run": False,
        },
        "comparison": {
            "fields": "all",
            "exact": "keys, types, strings and source/input hashes",
            "float_atol": 8e-8,
            "float_rtol": 8e-7,
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
        print("charged CTP neck source reproduced from authenticated inputs")
    elif args.output:
        if args.output.exists():
            raise FileExistsError("refusing to overwrite recorded evidence")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

