#!/usr/bin/env python3
"""Evaluate the NSC-to-MMP embedding gate without rerunning the throat solve."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.nsc_mmp_embedding import evaluate_mmp_embedding


OUTPUT = ROOT / "results/development/mmp-embedding.json"
INPUTS = (
    "results/development/charged-sector.json",
    "results/development/compact-matching.json",
    "results/development/compact-casimir.json",
    "results/development/vacuum-charge-matching.json",
)
SOURCES = (
    "src/recursive_horizons/nsc_mmp_embedding.py",
    "scripts/derive_nsc_mmp_embedding.py",
    "docs/nsc-mmp-embedding.md",
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


def calculate():
    records = [json.loads((ROOT / path).read_text()) for path in INPUTS]
    result = evaluate_mmp_embedding(*records).to_dict()
    result.update({
        "schema": "NSC-MMP-EMBEDDING-v1",
        "status": (
            "field, parity, state and coefficient conventions bind to the imported MMP sector; "
            "the retained partial coefficients fail the seed bound and V_full remains unresolved"
        ),
        "source_hashes": hashes(SOURCES),
        "input_hashes": hashes(INPUTS),
        "imported_relations": {
            "action": "integral sqrt(-g)*(A*R-C*F^2-V)",
            "G_N": "1/(16*pi*A)",
            "g_squared": "1/(4*C)",
            "r_e_squared": "q^2*C/(4*A)",
            "ell_MMP": "16*r_e^3/(abs(q)*G_N)",
            "Casimir_energy": "-abs(q)/(8*ell_MMP) in the MMP short-exterior regime",
        },
        "decision": {
            "MMP_gravity_ODE_rerun": False,
            "generic_fermion_throat_support_rederived": False,
            "current_partial_realization_embeds_MMP_seed": result["retained_partial_seed_passes"],
            "next_owner": result["first_unresolved_owner"],
            "dependent_recursive_solve_unblocked": result["full_embedding_decidable"]
                and result["retained_partial_seed_passes"],
        },
        "comparison": {
            "fields": "all",
            "exact": "keys, types, strings and source/input hashes",
            "float_atol": 3e-12,
            "float_rtol": 3e-11,
            "exceptions": [],
        },
    })
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = calculate()
    if args.check:
        compare(json.loads(OUTPUT.read_text()), result)
        print("MMP embedding gate reproduced from authenticated stored inputs")
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
