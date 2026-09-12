#!/usr/bin/env python3
"""Build the residual-bearing GeneralKSLocalInducedHistory record."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from recursive_horizons.nsc_general_ks_local_history import (  # noqa: E402
    CHANNELS,
    FIELDS,
    GeneralKSLocalInducedHistory,
)


OUTPUT = ROOT/"results/development/nsc-general-ks-local-history.json"
INPUTS = (
    "results/development/scale-binding.json",
    "results/development/nsc-landau-cauchy-isometry.json",
    "results/development/charged-ctp-neck-source.json",
    "results/development/spherical-action.json",
    "results/development/curvature-eft.json",
    "results/development/compact-boundary-action.json",
)
SOURCES = (
    "src/recursive_horizons/nsc_general_ks_local_history.py",
    "scripts/derive_nsc_general_ks_local_history.py",
    "docs/nsc-general-ks-local-history.md",
)


def hashes(paths):
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}


def compare(expected, actual, path="$", atol=3e-11, rtol=3e-11):
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
            raise AssertionError(f"number differs at {path}: {expected} != {actual}")
    elif type(expected) is not type(actual) or expected != actual:
        raise AssertionError(f"value differs at {path}: {expected!r} != {actual!r}")


def _load_control(record):
    path = ROOT/record["payload"]["path"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != record["payload"]["sha256"]:
        raise ValueError("control-history payload hash mismatch")
    with np.load(path, allow_pickle=False) as payload:
        # Thirty-three inherited nodes are sufficient to verify the owner API;
        # this remains a diagnostic control and is not a selected duration.
        index = np.linspace(0, len(payload["short_time"])-1, 33, dtype=int)
        return SimpleNamespace(
            time=np.array(payload["short_time"][index]),
            lapse=np.array(payload["short_lapse"][index]),
            shift=np.array(payload["short_shift"][index]),
            a_parallel=np.array(payload["short_a"][index]),
            radius=np.array(payload["short_r"][index]),
        )


def _directional_checks(owner, history, evaluated):
    s = (history.time-history.time[0])/(history.time[-1]-history.time[0])
    envelope = np.sin(np.pi*s)**2
    checks = []
    for field_name, attribute, phase in (
        ("lapse", "lapse", 0.1),
        ("q_ADM", "a_parallel", 0.4),
        ("r", "radius", 0.7),
    ):
        base = np.asarray(getattr(history, attribute))
        direction = base*envelope*np.cos(2*np.pi*s+phase)
        analytic = {
            channel: float(np.dot(evaluated["action_gradients"][channel][field_name], direction))
            for channel in CHANNELS
        }
        errors = []
        for step in (1e-4, 5e-5):
            plus = SimpleNamespace(**history.__dict__)
            minus = SimpleNamespace(**history.__dict__)
            setattr(plus, attribute, base+step*direction)
            setattr(minus, attribute, base-step*direction)
            observed = {
                channel: (owner.channel_actions(plus)[channel]-owner.channel_actions(minus)[channel])/(2*step)
                for channel in CHANNELS
            }
            errors.append(max(abs(observed[channel]-analytic[channel]) for channel in CHANNELS))
        checks.append({
            "field": field_name,
            "maximum_endpoint_variation": float(max(abs(direction[0]), abs(direction[-1]))),
            "endpoint_variation_tolerance": 1e-28,
            "endpoint_variation_zero": bool(max(abs(direction[0]), abs(direction[-1])) < 1e-28),
            "coarse_maximum_absolute_error": errors[0],
            "fine_maximum_absolute_error": errors[1],
        })
    tolerance = 5e-7
    maximum = max(row["fine_maximum_absolute_error"] for row in checks)
    return {"rows": checks, "maximum_absolute": maximum,
            "tolerance": tolerance, "pass": bool(maximum < tolerance)}


def calculate():
    scale, isometry, charged, spherical, curvature, boundary = (
        json.loads((ROOT/path).read_text()) for path in INPUTS
    )
    branch = scale["development_branch"]
    coefficients = branch["coefficients"]
    local = charged["runs"]["base"]["compact_local_source"]
    if local["C_Weyl"] != coefficients["C_Weyl"]:
        raise ValueError("compact/local Weyl allocation changed")
    if coefficients["V_full_relational"] != 0:
        raise ValueError("locked relational vacuum changed")
    if boundary["domain"]["boundary_type"] != "reflecting compact endpoints, not the parent-child transmitting throat":
        raise ValueError("compact boundary domain changed")
    if curvature["quantum_and_boundary_requirements"]["Euler_box_and_induced_boundary_terms_discarded"]:
        raise ValueError("curvature owner discarded required boundary terms")
    owner = GeneralKSLocalInducedHistory(
        coefficients["A_Wilsonian_at_magnetic_scale"],
        coefficients["C_Wilsonian_at_magnetic_scale"],
        coefficients["C_Weyl"], coefficients["C_Euler"],
        coefficients["C_boxR"], branch["magnetic_flux"],
    )
    history = _load_control(isometry)
    evaluated = owner.evaluate(history)
    checks = _directional_checks(owner, history, evaluated)
    if not checks["pass"] or evaluated["beta_variation_maximum_absolute"] != 0:
        raise ArithmeticError("node action gradients failed their focused check")
    endpoint = evaluated["weyl_endpoint_completion_residual"]
    boundary_support = {
        channel: max(
            abs(value)
            for field in FIELDS
            for value in evaluated["action_gradients"][channel][field][5:-5]
        )
        for channel in ("euler_endpoint", "boxR_endpoint")
    }
    if max(boundary_support.values()) >= 1e-12:
        raise ArithmeticError("Euler/boxR gradient leaked outside endpoint stencil")
    return {
        "schema": "NSC-GENERAL-KS-LOCAL-INDUCED-HISTORY-v1",
        "status": "COMPONENT PASS / COMPOSITION OPEN: node-wise local forces are executable; the free-endpoint Weyl completion is absent",
        "source_hashes": hashes(SOURCES), "input_hashes": hashes(INPUTS),
        "owner": {
            "name": "GeneralKSLocalInducedHistory",
            "history_fields": list(FIELDS),
            "variation": "node-wise dS_local/d(N,beta,q_ADM,r); action_force=-gradient",
            "history_selected_physical": False,
            "diagnostic_history": "33-node subsample of authenticated endpoint_control_duration_0.5",
        },
        "locked_inputs": {
            "A": coefficients["A_Wilsonian_at_magnetic_scale"],
            "C_gauge": coefficients["C_Wilsonian_at_magnetic_scale"],
            "C_Weyl": coefficients["C_Weyl"],
            "C_Euler": coefficients["C_Euler"],
            "C_boxR": coefficients["C_boxR"],
            "magnetic_flux": branch["magnetic_flux"],
            "Omega": branch["omega"], "zeta": branch["zeta"],
            "V_full": coefficients["V_full_relational"],
        },
        "domain": {
            "included": "smooth homogeneous Lorentzian KS histories with positive N,q_ADM,r and fixed spatially homogeneous magnetic flux",
            "einstein": "bulk spherical reduction with the already owned four-dimensional GHY completion",
            "maxwell": "bulk monopole term from the existing spherical action",
            "weyl": "four-dimensional C^2 bulk action in the original metric convention",
            "euler_and_boxR": "temporal endpoint functionals; no bulk force assigned",
            "compact_boundary_owner_import": boundary["domain"]["boundary_type"],
            "compact_boundary_not_repurposed": True,
            "excluded": "transmitting tilted interface, frequency-mixing state map, fourth-order reference history, physical history selection",
        },
        "allocation": {
            "local_action": "-int sqrt|g|[A R+C_gauge F^2+C_Weyl C^2+C_Euler E4+C_boxR boxR] with existing GHY",
            "compact_local_C_Weyl_count": 1,
            "stored_unit_radius_tensor_added_to_history": False,
            "extra_boundary_term_invented": False,
            "Euler_or_boxR_bulk_source_fabricated": False,
        },
        "diagnostic_evaluation": evaluated,
        "residuals": {
            "directional_gradient_check": checks,
            "beta_independence": {"maximum_absolute": evaluated["beta_variation_maximum_absolute"], "tolerance": 0.0, "pass": True},
            "endpoint_only_support": {
                "boundary_stencil_nodes_per_side": 5,
                "maximum_absolute_interior_gradient": max(boundary_support.values()),
                "by_channel": boundary_support,
                "tolerance": 1e-12,
                "pass": bool(max(boundary_support.values()) < 1e-12),
            },
            "weyl_endpoint_completion": endpoint,
        },
        "gate": {
            "nodewise_action_forces_pass": True,
            "constant_endpoint_force_substitution_used": False,
            "Euler_boxR_endpoint_ledger_pass": True,
            "weyl_bulk_pass": True,
            "weyl_free_endpoint_completion_pass": endpoint["below_tolerance"],
            "owner_component_pass": True,
            "ready_for_extended_composition": True,
            "extended_existence_gate_run": False,
            "coupled_evolution_reopened": False,
            "status": "OPEN",
            "next_missing_owner": "declared Weyl endpoint completion from the transmitting tilted interface owner",
        },
        "scope": {
            "A_q_Omega_zeta_or_V_full_changed": False,
            "finite_stress_or_nulls_fabricated": False,
            "metric_timestep_started": False,
            "old_scientific_generator_rerun": False,
            "new_physical_term_or_counterflow": False,
        },
        "upstream_scope_checks": {
            "spherical_self_sourced_solution": spherical["scope"]["self_sourced_geometry_or_physical_couplings_solved"],
            "curvature_boundary_terms_retained": True,
            "compact_transmitting_boundary_derived": boundary["scope"]["transmitting_throat_boundary_action_derived"],
        },
        "comparison": {"fields": "all", "exact": "schema, strings, booleans, shapes and hashes", "float_atol": 3e-11, "float_rtol": 3e-11, "exceptions": []},
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
        print("general-KS local history forces and open Weyl endpoint residual reproduced")
    elif args.output:
        if args.output.exists():
            raise FileExistsError("refusing to overwrite recorded evidence")
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
