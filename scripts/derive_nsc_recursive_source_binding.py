#!/usr/bin/env python3
"""Apply the scale-inherited transparent LLL state to the actual NSC neck."""
from __future__ import annotations

import argparse
import hashlib
import json
from math import pi, sqrt
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.nsc_recursive_source_binding import (
    bind_recursive_lll_to_neck,
)


OUTPUT = ROOT / "results/development/recursive-source-binding.json"
INPUTS = (
    "results/development/scale-binding.json",
    "results/development/horizon-source.json",
    "results/development/adm-neck-source-map.json",
)
SOURCES = (
    "src/recursive_horizons/nsc_recursive_source_binding.py",
    "scripts/derive_nsc_recursive_source_binding.py",
    "docs/nsc-recursive-source-binding.md",
)


def hashes(paths):
    return {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in paths}


def compare(expected, actual, path="$", atol=4e-11, rtol=4e-10):
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


def calculate():
    scale, horizon, neck = (json.loads((ROOT / path).read_text()) for path in INPUTS)
    branch = scale["development_branch"]
    q = branch["magnetic_flux"]
    omega = branch["omega"]
    kappa = horizon["benchmark_Unruh_per_abs_q"]["reused_surface_gravity"]
    geometry = neck["geometry"]
    source = bind_recursive_lll_to_neck(
        magnetic_flux=q,
        omega=omega,
        surface_gravity=kappa,
        metric_A=1.0 - 3.0 * pi / 2.0,
        metric_A_prime=6.0,
        metric_A_second=-3.0 * pi,
        beta=geometry["beta"],
        radius=geometry["r"],
    )
    expected_power = source.outgoing_state_constant - source.incoming_state_constant
    if abs(source.parent_killing_power - expected_power) > 5e-14:
        raise AssertionError("state constants and projected Killing power differ")
    if source.source_null_signs_match:
        raise AssertionError("stored source-sign obstruction unexpectedly disappeared")
    return {
        "schema": "NSC-RECURSIVE-SOURCE-BINDING-v1",
        "status": (
            "scale-inherited transparent LLL state evaluated on the actual neck; "
            "it does not supply the required black-universe null source"
        ),
        "source_hashes": hashes(SOURCES),
        "input_hashes": hashes(INPUTS),
        "state_law": {
            "parent": "t_u=abs(q)*kappa_p^2/(48*pi)",
            "child_in_parent_energy_frame": "kappa_c=Omega*kappa_p",
            "recursive_incoming": "t_v=abs(q)*kappa_c^2/(48*pi)=Omega^2*t_u",
            "transmission": "one for the imported free massless magnetic LLL",
            "external_incoming_reservoir_inserted": False,
        },
        "scale_branch": {
            "q": q,
            "Omega": omega,
            "zeta": branch["zeta"],
        },
        "source": source.to_dict(),
        "identities": {
            "t_v_over_t_u_minus_zeta": (
                source.incoming_state_constant / source.outgoing_state_constant
                - branch["zeta"]
            ),
            "power_minus_tu_plus_tv": (
                source.parent_killing_power
                - source.outgoing_state_constant
                + source.incoming_state_constant
            ),
            "beta_squared_minus_three_pi_over_two": geometry["beta"]**2 - 3*pi/2,
        },
        "gate": {
            "recursive_state_is_directed": source.parent_killing_power != 0,
            "static_flux_balance": source.parent_killing_power == 0,
            "both_required_null_signs_supplied": source.source_null_signs_match,
            "MMP_LLL_can_be_substituted_for_full_black_universe_CTP_source": False,
            "next_source_owner": (
                "full charged angular/compact CTP covariance and metric variation "
                "on the horizon-penetrating domain"
            ),
        },
        "scope": {
            "new_metric_equation_derived": False,
            "MMP_geometry_recomputed": False,
            "old_angular_or_scattering_generator_rerun": False,
            "interroom_clock_law_from_global_geometry_proved": False,
            "cosmological_Q_identified_with_Killing_power": False,
            "physical_solution_claimed": False,
        },
        "comparison": {
            "fields": "all",
            "exact": "keys, types, strings and source/input hashes",
            "float_atol": 4e-11,
            "float_rtol": 4e-10,
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
        print("recursive source binding reproduced from authenticated inputs")
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
