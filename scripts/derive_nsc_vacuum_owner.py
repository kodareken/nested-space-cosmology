#!/usr/bin/env python3
"""Identify which existing common-action term can own the MMP vacuum mismatch."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/development/mmp-vacuum-owner.json"
INPUTS = (
    "results/development/mmp-embedding.json",
    "results/development/causal-common-functional.json",
    "results/development/compact-matching.json",
    "results/development/compact-boundary-action.json",
    "results/development/compact-casimir.json",
    "results/development/canonical-spectral-bridge.json",
)
SOURCES = (
    "scripts/derive_nsc_vacuum_owner.py",
    "docs/nsc-mmp-vacuum-owner.md",
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
        if abs(expected - actual) > 3e-13 + 3e-12 * abs(expected):
            raise AssertionError(f"numeric value differs at {path}")
    elif type(expected) is not type(actual) or expected != actual:
        raise AssertionError(f"exact value differs at {path}")


def calculate():
    embedding, causal, matching, boundary, casimir, bridge = (
        json.loads((ROOT / path).read_text()) for path in INPUTS
    )
    row = next(item for item in matching["matched_coefficients"]
               if item["matching_cutoff"] == 1.0)
    vacuum, einstein, gauge = row["V_Dirac"], row["A_Dirac"], row["C_gauge_Dirac"]
    maximum = 2 * einstein * einstein / gauge
    ratio = maximum / vacuum
    assert abs(embedding["coefficient_map"]["Xi_per_q_squared"]
               - vacuum * gauge / (8 * einstein * einstein)) < 2e-13
    assert bridge["scope"]["new_independent_coefficients"] is False
    assert boundary["scope"]["full_renormalized_vacuum_fixed"] is False
    assert casimir["holonomy"]["complete_state_or_recursive_return_path_derived"] is False
    assert causal["selected_real_time_owner"]["continued_bare_Euclidean_kernel_used_as_retarded_response"] is False

    return {
        "schema": "NSC-MMP-VACUUM-OWNER-v1",
        "status": (
            "the canonical-CTP revision leaves the state-independent bulk vacuum coefficient "
            "unfixed; no existing same-action contribution owns the required cancellation"
        ),
        "source_hashes": hashes(SOURCES),
        "input_hashes": hashes(INPUTS),
        "retained_partial_coefficients": {
            "V": vacuum, "A": einstein, "C": gauge,
            "Xi_per_q_squared": vacuum * gauge / (8 * einstein * einstein),
        },
        "charged_seed_requirement": {
            "V_full_max_for_abs_q_1": maximum,
            "V_full_over_V_partial_max_for_abs_q_1": ratio,
            "minimum_fraction_of_partial_V_that_must_be_removed": 1 - ratio,
            "general": "V_full <= 2*A_full^2/(q^2*C_full)",
            "asymptotically_flat_MMP": "V_full=0",
        },
        "existing_owner_classification": [
            {
                "term": "retained compact/light spectral volume coefficient",
                "can_change_asymptotic_bulk_V": True,
                "value_known": True,
                "value": vacuum,
            },
            {
                "term": "compact boundary heat terms",
                "can_change_asymptotic_bulk_V": False,
                "reason": "boundary-supported variation is not a homogeneous four-volume term",
            },
            {
                "term": "Weyl, Euler, box-R and higher-curvature terms",
                "can_change_asymptotic_bulk_V": False,
                "reason": "their bulk densities vanish or remain curvature-dependent in the flat MMP asymptotic region",
            },
            {
                "term": "finite Casimir interaction",
                "can_change_asymptotic_bulk_V": False,
                "reason": "topology and return-length dependent interaction, not a homogeneous vacuum coefficient",
            },
            {
                "term": "canonical CTP state determinant",
                "can_change_asymptotic_bulk_V": False,
                "reason": "equal-history normalization fixes causal state response but not the state-independent local vacuum counterterm",
            },
            {
                "term": "recursive tail",
                "can_change_asymptotic_bulk_V": None,
                "reason": "physical return domain and stationary recursive solution are downstream and not yet defined",
            },
        ],
        "one_revision_attempt": {
            "realization": "canonical Lorentzian CTP with spectral induced coefficients counted once",
            "causal_owner_resolved": causal["gate"]["finite_interface_executable"],
            "bulk_V_owner_resolved": False,
            "reason": "the causal/state split does not supply a normalization condition for V_full",
        },
        "decision": {
            "current_parameter_free_MMP_embedding": False,
            "dependent_recursive_scale_solve_started": False,
            "coefficient_fitted_or_set_to_zero": False,
            "additional_sector_added": False,
            "blocking_definition": (
                "a same-action ultraviolet matching or normalization condition that determines V_full"
            ),
            "next_allowed_input": (
                "an explicit microscopic spectrum/measure completing the supertrace, or a derived recursive "
                "normalization condition; an imposed cosmological counterterm is not sufficient"
            ),
        },
        "comparison": {
            "fields": "all", "exact": "keys, types, strings and source/input hashes",
            "float_atol": 3e-13, "float_rtol": 3e-12, "exceptions": [],
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
        print("MMP vacuum owner gate reproduced from authenticated inputs")
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
