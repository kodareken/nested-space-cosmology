#!/usr/bin/env python3
"""Bind canonical CTP evolution to induced source accounting in one finite owner."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.nsc_causal_common import (
    BoundaryDomain,
    CausalCommonFunctional,
    GaugeHistory,
    LinkHistory,
    MetricHistory,
    SpectralInducedSource,
)
from recursive_horizons.nsc_influence import ground_covariance


OUTPUT = ROOT / "results/development/causal-common-functional.json"
INPUTS = (
    "results/development/adm-source-constraints.json",
    "results/development/warp-local-neck-source.json",
    "results/nsc-11-response-matching.json",
    "results/development/mmp-embedding.json",
)
SOURCES = (
    "src/recursive_horizons/nsc_causal_common.py",
    "scripts/derive_nsc_causal_common.py",
    "docs/nsc-causal-common-functional.md",
)


def hashes(paths):
    return {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in paths}


def compare(expected, actual, path="$"):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys():
            raise AssertionError(f"keys differ at {path}")
        for key in expected:
            compare(expected[key], actual[key], f"{path}/{key}")
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise AssertionError(f"list differs at {path}")
        for index, (left, right) in enumerate(zip(expected, actual)):
            compare(left, right, f"{path}/{index}")
    elif isinstance(expected, float):
        if isinstance(actual, bool) or not isinstance(actual, (int, float)):
            raise AssertionError(f"numeric type differs at {path}")
        if abs(expected - actual) > 3e-12 + 3e-11 * abs(expected):
            raise AssertionError(f"numeric value differs at {path}")
    elif type(expected) is not type(actual) or expected != actual:
        raise AssertionError(f"exact value differs at {path}")


def finite_control(induced_forces):
    times = np.array([0.0, 0.13, 0.31])
    parent = np.repeat(np.array([[[-0.8]]], complex), len(times), axis=0)
    child = np.repeat(np.array([[[0.9]]], complex), len(times), axis=0)
    links = np.repeat(np.array([[[0.23 + 0.04j]]], complex), len(times), axis=0)
    full = np.array([[-0.8, 0.23 + 0.04j], [0.23 - 0.04j, 0.9]])
    _, _, covariance = ground_covariance(full)
    vertices = {
        "N": np.array([[0.4, 0.1], [0.1, -0.2]]),
        "beta": np.array([[0.0, 0.2j], [-0.2j, 0.0]]),
        "q": np.diag([0.3, -0.1]),
        "r": np.array([[0.0, 0.07], [0.07, 0.0]]),
    }
    metric = MetricHistory(times, parent, child, vertices)
    gauge = GaugeHistory(np.zeros_like(parent), np.zeros_like(child))
    link = LinkHistory(links)
    spectral = SpectralInducedSource(
        forces=induced_forces,
        coefficient_owner="stored local warp source; magnitude used only for one-counting control",
        unresolved_terms=(
            "full common-action vacuum coefficient V_full",
            "physical charged parent-child covariance and return history",
        ),
    )
    domain = BoundaryDomain(
        "finite two-mode CTP accounting control", 1, 1,
        "common canonical frame", "hbar=1; control matrices dimensionless",
        stress_factors={
            "rho": ("N", 2.0), "T_01": ("beta", 3.0),
            "p_parallel": ("q", -4.0), "p_perp": ("r", -5.0),
        },
        power_weights={"N": -0.5, "beta": 0.25, "q": 0.75},
    )
    return CausalCommonFunctional().evaluate(
        metric, gauge, link, covariance, spectral, domain
    )


def calculate():
    adm, warp, matching, embedding = (
        json.loads((ROOT / path).read_text()) for path in INPUTS
    )
    local = warp["counted_source_total"]["force_densities"]
    induced = {
        "N": local["F_N"], "beta": local["F_beta"],
        "q": local["F_q"], "r": local["F_r"],
    }
    source = finite_control(induced)
    identity = matching["finite_grid_frequency_identity"]
    max_binding = max(abs(row["difference"]) for row in identity)
    return {
        "schema": "NSC-CAUSAL-COMMON-FUNCTIONAL-v1",
        "status": (
            "canonical CTP owner implemented with one-counted induced forces; "
            "actual charged parent-child history remains the first physical input"
        ),
        "source_hashes": hashes(SOURCES),
        "input_hashes": hashes(INPUTS),
        "selected_real_time_owner": {
            "functional": "S_induced[+]-S_induced[-]-i log det(I-C0+C0 Uminus_dagger Uplus)",
            "state_dependent_nonlocal_response": "canonical Lorentzian Dirac evolution",
            "state_independent_local_coefficients": "same spectral coefficient account, added once",
            "continued_bare_Euclidean_kernel_used_as_retarded_response": False,
        },
        "interface": {
            "inputs": [
                "MetricHistory", "GaugeHistory", "LinkHistory", "initial_covariance",
                "SpectralInducedSource", "BoundaryDomain",
            ],
            "outputs": [
                "Gamma_CTP", "F_N", "F_beta", "F_q", "F_r", "stress",
                "boundary_power", "retarded_response", "noise_kernel",
                "ward_residuals", "unresolved_terms",
            ],
        },
        "finite_control": source.to_dict(),
        "equilibrium_binding_reused": {
            "frequencies": [row["nu"] for row in identity],
            "maximum_direct_loop_minus_retarded_at_i_nu": max_binding,
            "finite_frequency_loop_matches_retarded": matching["gate"]["finite_frequency_loop_matches_retarded"],
            "covariant_static_owner_difference": matching["covariant_static_owner_match"]["difference"],
            "finite_proper_time_equals_canonical_Kubo": matching["finite_cutoff_distinction"]["equal_at_finite_cutoff"],
            "meaning": matching["finite_cutoff_distinction"]["meaning"],
        },
        "ADM_binding_reused": {
            "exact_Ward_and_stress_identities": adm["exact_identities"],
            "full_absolute_state_source_evaluated": adm["scope"]["full_absolute_state_source_evaluated"],
        },
        "gate": {
            "finite_interface_executable": True,
            "equal_history_normalized": source.ward_residuals["equal_history_action"] < 2e-12,
            "finite_energy_work_balance": abs(source.ward_residuals["energy_work_balance"]) < 2e-12,
            "induced_source_counted_once": all(
                abs(source.forces[name] - source.matter_forces[name] - source.induced_forces[name]) < 2e-14
                for name in ("N", "beta", "q", "r")
            ),
            "actual_physical_CausalSource_complete": False,
            "MMP_embedding_unblocked": embedding["decision"]["dependent_recursive_solve_unblocked"],
        },
        "first_missing_dependency": (
            "full common-action V_full, followed by the charged parent-child Hamiltonian, "
            "covariance and closed magnetic return history"
        ),
        "old_scientific_generators_rerun": False,
        "comparison": {
            "fields": "all", "exact": "keys, types, strings and source/input hashes",
            "float_atol": 3e-12, "float_rtol": 3e-11, "exceptions": [],
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
        print("causal common functional reproduced from authenticated inputs")
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
