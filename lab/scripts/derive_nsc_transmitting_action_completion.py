#!/usr/bin/env python3
"""Apply the existing Gaussian CTP action differential to locked NSC channels."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from recursive_horizons.nsc_influence import influence
from recursive_horizons.nsc_mode_resolved_cauchy_state import array_digest
from recursive_horizons.nsc_tilted_landau_interface import full_covariance
from recursive_horizons.nsc_transmitting_boundary_history import EndpointVariation, FIELDS
from recursive_horizons.nsc_transmitting_ctp_variation import (
    ctp_first_variation, endpoint_ctp_pullback, completion_gate,
)
from derive_nsc_transmitting_boundary_binding import compare

INPUTS = {
    "state": "results/development/nsc-mode-resolved-cauchy-state.json",
    "selection": "results/development/nsc-transmitting-boundary-selection.json",
    "endpoint": "results/development/nsc-weyl-endpoint-match.json",
    "extended": "results/development/nsc-full-extended-history-gate.json",
}
SOURCES = (
    "src/recursive_horizons/nsc_transmitting_ctp_variation.py",
    "scripts/derive_nsc_transmitting_action_completion.py",
    "tests/test_nsc_transmitting_ctp_variation.py",
    "docs/nsc-transmitting-ctp-variation.md",
    "docs/nsc-transmitting-metric-pullback.md",
    "docs/nsc-extended-action-completion-gate.md",
)
DEPENDENCIES = (
    "src/recursive_horizons/nsc_influence.py",
    "src/recursive_horizons/nsc_boundary_state.py",
    "src/recursive_horizons/nsc_causal_common.py",
    "src/recursive_horizons/nsc_mode_resolved_cauchy_state.py",
    "src/recursive_horizons/nsc_tilted_landau_interface.py",
    "src/recursive_horizons/nsc_transmitting_boundary_history.py",
    "scripts/derive_nsc_transmitting_boundary_binding.py",
)
OUTPUTS = ("nsc-transmitting-ctp-variation.json", "nsc-transmitting-metric-pullback.json", "nsc-extended-action-completion-gate.json")


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def calculate():
    inputs = {key: json.loads((ROOT/path).read_text()) for key, path in INPUTS.items()}
    state = inputs["state"]
    locked = inputs["selection"]["locked_inputs"]
    for key in ("A", "magnetic_flux", "Omega", "zeta", "V_full"):
        if locked[key] != state["locked_inputs"][key]:
            raise ValueError("locked parameter drift: "+key)
    payload = state["payload"]
    if sha(payload["path"]) != payload["sha256"]:
        raise ValueError("Cauchy-state artifact mismatch")
    with np.load(ROOT/payload["path"], allow_pickle=False) as f:
        covariance = np.array(f["covariance_seed"])
        offsets = np.array(f["sample_offsets"])
    for name, array in (("covariance_seed", covariance), ("sample_offsets", offsets)):
        if array_digest(name, array) != payload["arrays"][name]["sha256"]:
            raise ValueError("Cauchy-state array mismatch: "+name)
    rows = []
    eps = 2e-5  # Dimensionless differentiation step, never a history duration.
    for channel, (start, stop) in enumerate(zip(offsets[:-1], offsets[1:])):
        c = full_covariance(covariance[start:stop])
        size = len(c)
        u = np.eye(size, dtype=complex)
        x = u/np.sqrt(size)
        common = ctp_first_variation(c, u, u, 1j*x, 1j*x)
        relative = ctp_first_variation(c, u, u, 0.5j*x, -0.5j*x)
        def action(s):
            return influence(c, np.exp(0.5j*s/np.sqrt(size))*u,
                             np.exp(-0.5j*s/np.sqrt(size))*u)["principal_action"]
        numeric = (action(eps)-action(-eps))/(2*eps)
        exact = float(np.trace(c).real/np.sqrt(size))
        rows.append({
            "channel": channel, "canonical_dimension": size,
            "common_branch_derivative": [common["derivative"].real, common["derivative"].imag],
            "relative_branch_derivative": [relative["derivative"].real, relative["derivative"].imag],
            "normalized_trace_identity": exact,
            "existing_action_finite_difference": [numeric.real, numeric.imag],
            "first_differential_error": float(abs(relative["derivative"]-numeric)),
            "trace_identity_residual": float(abs(relative["derivative"]-exact)),
            "unitarity_and_tangent_residual": max(common["unitarity_residual"], common["tangent_residual"], relative["tangent_residual"]),
            "linear_solve_residual": relative["linear_solve_residual"],
        })
    residuals = {
        "maximum_first_differential_error": max(r["first_differential_error"] for r in rows),
        "first_differential_tolerance": 3e-8,
        "maximum_trace_identity_residual": max(r["trace_identity_residual"] for r in rows),
        "maximum_unitarity_and_tangent_residual": max(r["unitarity_and_tangent_residual"] for r in rows),
        "algebra_tolerance": 3e-11,
        "maximum_common_branch_derivative": max(abs(complex(*r["common_branch_derivative"])) for r in rows),
        "minimum_relative_identity_derivative": min(r["normalized_trace_identity"] for r in rows),
        "maximum_relative_identity_derivative": max(r["normalized_trace_identity"] for r in rows),
    }
    if residuals["maximum_first_differential_error"] > residuals["first_differential_tolerance"] or residuals["maximum_trace_identity_residual"] > residuals["algebra_tolerance"]:
        raise ArithmeticError("first differential failed the existing-action comparison")
    endpoint = inputs["endpoint"]["endpoint"]
    basis = EndpointVariation(endpoint["basis"]["id"], tuple(tuple(endpoint["local_components"][field]) for field in FIELDS))
    # A missing physical Jacobian remains missing. The conditional pullback is
    # tested with explicitly test-only jets, never filled using diagnostic U.
    pullback = endpoint_ctp_pullback(None, None, None, basis, None)
    common = {
        "source_hashes": {path: sha(path) for path in SOURCES},
        "input_hashes": {path: sha(path) for path in (*INPUTS.values(), *DEPENDENCIES)},
        "locked_inputs": locked,
        "mode_state_payload": {"path": payload["path"], "sha256": payload["sha256"]},
        "comparison": inputs["selection"]["comparison"],
        "scope": {"parameters_refitted": False, "new_physical_action_term": False,
                  "source_integration_run": False, "duration_selected": False,
                  "metric_timestep_started": False, "old_generator_rerun": False,
                  "finite_stress_fabricated": False},
    }
    selection = {
        "status": "OPEN",
        "conditional_differential_status": "PASS",
        "functional": "Gamma_G=-i log det Q; Q=I-C0+C0 Vminus^dagger Vplus",
        "differential": "dGamma_G=-i Tr[Q^-1 C0(dVminus^dagger Vplus+Vminus^dagger dVplus)]",
        "generator": "X=I/sqrt(d_c) in each retained channel; ||X||_F=1",
        "difference_convention": "Vplus(s)=exp(i s X/2)V; Vminus(s)=exp(-i s X/2)V; V=I is an algebra control only",
        "dimensionless_derivative_step": eps,
        "rows": rows, "residuals": residuals,
        "physical_selection_residual": None, "selected_Vc": None,
        "reason": "common-branch variation is identically flat; the relative derivative is a source differential, not an independent equation extremizing unrestricted V",
        "central_derivative_invariant_under_V": True,
        "normalization_note": "per-channel finite canonical determinant; not a continuum angular sum, renormalized stress, or total same-action derivative",
        "previous_covariance_freedom_probe_reused": inputs["selection"]["selection"]["covariance_freedom_probe"]["maximum_state_tangent"],
        "missing_action_ingredient": "the physical map U[g,B,embedding] and its admissible metric/link variations, including the same-action non-Gaussian boundary completion",
        "nonexistence_of_full_extended_class_claimed": False,
    }
    metric = {
        "status": "OPEN", "conditional_gaussian_pullback": pullback,
        "endpoint_basis": endpoint["basis"],
        "locked_local_components": endpoint["local_components"],
        "locked_extraction_residual": endpoint["binding_residual"],
        "preserved_diagnostic_local_maximum": endpoint["local_maximum_absolute"],
        "coordinate_aliases": {"N": "N", "beta": "beta", "q": "q_ADM", "r": "r"},
        "coordinate_aliases_are_not_geometric_embedding": True,
        "missing_derivative_slots": {
            "dVplus_dgDelta": None, "dVminus_dgDelta": None,
            "boundary_metric_and_normal_jet_pullback": None,
            "remaining_same_action_boundary_derivative": None,
        },
        "metric_mismatch": {"components": None, "maximum_absolute": None, "tolerance": 3e-11, "status": "OPEN"},
        "independent_first_jet_derivative": endpoint["independent_first_jet_derivative"],
        "source_counting": "Gaussian differential only; existing reference/local pieces counted by their owners, never added as an independent cancelling term",
    }
    gate = completion_gate(inputs["extended"], selection, metric)
    records = (
        dict(common, schema="NSC-TRANSMITTING-CTP-VARIATION-v1", status="OPEN", selection=selection),
        dict(common, schema="NSC-TRANSMITTING-METRIC-PULLBACK-v1", status="OPEN", metric_pullback=metric),
        dict(common, schema="NSC-EXTENDED-ACTION-COMPLETION-GATE-v1", status=gate["decision"], composition=gate,
             prerequisites=inputs["extended"]["component_gates"], new_residuals=residuals),
    )
    return dict(zip(OUTPUTS, records))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--write", action="store_true")
    args = parser.parse_args()
    records = calculate()
    folder = ROOT/"results/development"
    if args.check:
        for name, record in records.items(): compare(json.loads((folder/name).read_text()), record)
        print("F/G action differential and extended OPEN gate: all fields verified")
    elif args.write:
        if any((folder/name).exists() for name in records): raise FileExistsError("refusing to overwrite locked records")
        for name, record in records.items():
            (folder/name).write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+"\n")
            print(name+": "+record["status"])
    else:
        print(json.dumps(records,indent=2,sort_keys=True,allow_nan=False))


if __name__ == "__main__": main()
