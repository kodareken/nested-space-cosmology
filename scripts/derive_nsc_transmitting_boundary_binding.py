#!/usr/bin/env python3
"""Run steps A--D of the full extended gate using authenticated existing inputs.

Only the covariance freedom probe and endpoint binding are new computations.
Historical generators, source integrals and time evolution are never called.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_extended_tilted_history_gate import compose_extended_gate
from recursive_horizons.nsc_mode_resolved_cauchy_state import array_digest
from recursive_horizons.nsc_transmitting_boundary_history import (
    ModeResolvedTransmittingBoundaryHistoryAction, compose_full_extended_gate,
)

INPUTS = {
    "reference": "results/development/nsc-general-ks-reference.json",
    "local": "results/development/nsc-general-ks-local-history.json",
    "interface": "results/development/nsc-tilted-landau-interface.json",
    "homogeneous": "results/development/nsc-general-ks-same-action-history.json",
    "boundary": "results/development/compact-boundary-action.json",
    "state": "results/development/nsc-mode-resolved-cauchy-state.json",
    "scale": "results/development/scale-binding.json",
}
SOURCES = (
    "src/recursive_horizons/nsc_transmitting_boundary_history.py",
    "scripts/derive_nsc_transmitting_boundary_binding.py",
    "tests/test_nsc_transmitting_boundary_history.py",
    "docs/nsc-transmitting-boundary-selection.md",
    "docs/nsc-weyl-endpoint-match.md",
    "docs/nsc-full-extended-history-gate.md",
)
EVIDENCE = (
    ("src/recursive_horizons/nsc_causal_common.py", "LinkHistory", "values is supplied; no transmitting kernel law is encoded by this data container"),
    ("src/recursive_horizons/nsc_causal_common.py", "CausalCommonFunctional._assemble", "consumes metric, gauge and link histories; does not determine the missing link history"),
    ("src/recursive_horizons/nsc_causal_common.py", "CausalCommonFunctional.evaluate", "requires histories, initial covariance, induced coefficients and boundary domain"),
    ("src/recursive_horizons/nsc_boundary_state.py", "GaussianBoundaryState.kernels", "two-time kernels require link_t, child_t, link_s, child_s"),
    ("src/recursive_horizons/nsc_landau_cauchy_isometry.py", "propagate_blockwise", "propagates a supplied history in diagonal frequency blocks"),
    ("src/recursive_horizons/nsc_boundary.py", "joined_operator", "fixed spatial interval operator, not a selected tilted spacetime embedding"),
    ("src/recursive_horizons/nsc_tilted_landau_interface.py", "TransmittingTiltedLandauInterface.coordinate_kernel", "requires a canonical_unitary argument; cannot select it from the norm equation"),
    ("src/recursive_horizons/nsc_general_ks_local_history.py", "GeneralKSLocalInducedHistory.evaluate", "provides nodal local action gradients and a diagnostic endpoint projection"),
)
DOCUMENTS = (
    "docs/nsc-compact-boundary-action.md", "docs/nsc-spherical-action.md",
    "docs/nsc-curvature-eft.md", "docs/nsc-tilted-landau-interface.md",
    "src/recursive_horizons/nsc_extended_tilted_history_gate.py",
    "src/recursive_horizons/nsc_mode_resolved_cauchy_state.py",
)
OUTPUT_NAMES = (
    "nsc-transmitting-boundary-selection.json", "nsc-weyl-endpoint-match.json",
    "nsc-full-extended-history-gate.json",
)


def digest(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def inventory():
    rows = []
    for path, qualified, scope in EVIDENCE:
        node = ast.parse((ROOT/path).read_text())
        for part in qualified.split("."):
            node = next(n for n in node.body if isinstance(n, (ast.ClassDef, ast.FunctionDef)) and n.name == part)
        row = {"path": path, "symbol": qualified, "line": node.lineno,
               "end_line": node.end_lineno, "binding_scope": scope}
        if isinstance(node, ast.FunctionDef):
            row["arguments"] = [a.arg for a in node.args.args + node.args.kwonlyargs]
        else:
            row["declared_fields"] = [n.target.id for n in node.body if isinstance(n, ast.AnnAssign) and isinstance(n.target, ast.Name)]
        rows.append(row)
    if rows[0]["declared_fields"] != ["values"]:
        raise ValueError("LinkHistory changed; inspect its new selection capabilities")
    if "canonical_unitary" not in rows[6]["arguments"]:
        raise ValueError("tilted kernel provider changed")
    return rows


def compare(expected, actual, path="$"):
    if type(expected) is not type(actual):
        raise AssertionError(f"type mismatch at {path}")
    if isinstance(expected, dict):
        if expected.keys() != actual.keys():
            raise AssertionError(f"key mismatch at {path}")
        for k in expected:
            compare(expected[k], actual[k], path+"/"+k)
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise AssertionError(f"length mismatch at {path}")
        for i, (a, b) in enumerate(zip(expected, actual)):
            compare(a, b, path+"/"+str(i))
    elif isinstance(expected, float):
        if not math.isfinite(expected) or not math.isfinite(actual) or not math.isclose(expected, actual, abs_tol=3e-13, rel_tol=3e-13):
            raise AssertionError(f"numeric mismatch at {path}")
    elif expected != actual:
        raise AssertionError(f"value mismatch at {path}")


def calculate():
    records = {k: json.loads((ROOT/p).read_text()) for k, p in INPUTS.items()}
    state, scale, local = (records[k] for k in ("state", "scale", "local"))
    locked = local["locked_inputs"]
    branch = scale["development_branch"]
    for name in ("A", "magnetic_flux", "Omega", "zeta", "V_full"):
        if locked[name] != state["locked_inputs"][name]:
            raise ValueError(f"locked input differs: {name}")
    if (locked["magnetic_flux"], locked["Omega"], locked["zeta"], locked["V_full"]) != (4, branch["omega"], branch["zeta"], 0):
        raise ValueError("locked scale branch differs")
    payload = state["payload"]
    if digest(payload["path"]) != payload["sha256"]:
        raise ValueError("mode-state artifact authentication failed")
    with np.load(ROOT/payload["path"], allow_pickle=False) as f:
        arrays = {key: np.array(f[key]) for key in ("sample_offsets", "covariance_seed", "frequency")}
    for key, value in arrays.items():
        spec = payload["arrays"][key]
        if list(value.shape) != spec["shape"] or value.dtype.str != spec["dtype"] or array_digest(key, value) != spec["sha256"]:
            raise ValueError(f"mode-state array differs: {key}")
    owner = ModeResolvedTransmittingBoundaryHistoryAction()
    selection = owner.selection(records["interface"], records["boundary"], arrays)
    endpoint = owner.endpoint(local, local_record_sha256=digest(INPUTS["local"]))
    probe = selection["covariance_freedom_probe"]
    if probe["maximum_isometry_tangent_residual"] > probe["tolerance"]:
        raise ArithmeticError("covariance probe leaves the admissible interface tangent space")
    if probe["channels_with_state_tangent_above_tolerance"] == 0:
        raise ArithmeticError("the proposed probe does not establish state transport ambiguity")
    base = compose_extended_gate(records["reference"], local, records["interface"], records["homogeneous"])
    gate = compose_full_extended_gate(base, selection, endpoint)
    if not gate["old_homogeneous_nonexistence"]:
        raise ValueError("imported homogeneous regression certificate changed")
    evidence_paths = sorted(set(p for p, _, _ in EVIDENCE) | set(DOCUMENTS))
    common = {
        "source_hashes": {p: digest(p) for p in SOURCES},
        "input_hashes": {p: digest(p) for p in INPUTS.values()},
        "action_evidence_hashes": {p: digest(p) for p in evidence_paths},
        "mode_state_artifact": {"path": payload["path"], "sha256": payload["sha256"]},
        "locked_inputs": locked,
        "scope": {"parameters_refitted": False, "new_physical_term": False,
                  "physical_duration_chosen": False, "finite_stress_fabricated": False,
                  "metric_timestep_started": False, "old_generator_called": False},
        "comparison": {"fields": "all", "float_atol": 3e-13, "float_rtol": 3e-13,
                       "exact": "types, keys, strings, integers, booleans, nulls and hashes", "exceptions": []},
    }
    a = dict(common, schema="NSC-TRANSMITTING-BOUNDARY-SELECTION-v1", status="OPEN", owner="ModeResolvedTransmittingBoundaryHistoryAction", selection=selection, action_inventory=inventory())
    b = dict(common, schema="NSC-WEYL-ENDPOINT-MATCH-v1", status="OPEN", endpoint=endpoint,
             action_declaration=records["boundary"]["metric_boundary_matching"],
             interpretation="93.54264532195464 is a diagnostic end-node local gradient; the two-sided physical mismatch is unevaluated until the boundary action and pullback are supplied")
    c = dict(common, schema="NSC-FULL-EXTENDED-HISTORY-GATE-v1", status=gate["decision"],
             steps={"A": selection["status"], "B": endpoint["match"]["status"], "C": "residual-bearing OPEN composition", "D": gate["decision"]},
             composition=gate, component_gates=base,
             upstream_numerical_controls={
                 "reference_bloch": records["reference"]["residuals"]["maximum_bloch_recursion_residual"],
                 "reference_parity_beta": records["reference"]["residuals"]["maximum_parity_completed_beta_projection"],
                 "local_node_gradient": local["residuals"]["directional_gradient_check"]["maximum_absolute"],
                 "local_node_gradient_tolerance": local["residuals"]["directional_gradient_check"]["tolerance"],
                 "diagnostic_local_endpoint_gradient": endpoint["local_maximum_absolute"],
                 "interface_unitarity": records["interface"]["residuals"]["maximum_canonical_unitarity_max_abs"],
                 "interface_tolerance": records["interface"]["residuals"]["declared_tolerance"],
                 "source": "unchanged authenticated component records; full controls remain there",
             })
    return dict(zip(OUTPUT_NAMES, (a, b, c)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    group = p.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--write", action="store_true")
    args = p.parse_args()
    results = calculate()
    folder = ROOT/"results/development"
    if args.check:
        for name, value in results.items():
            compare(json.loads((folder/name).read_text()), value)
        print("A/B binding records and full extended OPEN gate: all fields verified")
    elif args.write:
        if any((folder/name).exists() for name in results):
            raise FileExistsError("refusing to overwrite existing records")
        for name, value in results.items():
            (folder/name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+"\n")
            print(name+": "+value["status"])
    else:
        print(json.dumps(results, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
