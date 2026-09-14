#!/usr/bin/env python3
"""Test whether T's complete trace pair determines the requested Hamiltonian B.

This is a representation certificate, not a physical-link constructor. It
uses the authenticated T artifact and source h, with no source/jet replay.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from recursive_horizons.nsc_dirac_endpoint_jets import bulk_hamiltonian_and_vertices
from recursive_horizons.nsc_influence import _covariance
from derive_nsc_transmitting_boundary_binding import compare

OUTPUT = "results/development/nsc-hamiltonian-trace-representation.json"
INPUTS = (
    "results/development/nsc-transmitting-dirac-domain.json",
    "results/development/nsc-mode-resolved-cauchy-state.json",
)
SOURCES = (
    "scripts/check_nsc_hamiltonian_trace_representation.py",
    "docs/nsc-hamiltonian-trace-representation.md",
)
DEPENDENCIES = (
    "src/recursive_horizons/nsc_causal_common.py",
    "src/recursive_horizons/nsc_influence.py",
    "src/recursive_horizons/nsc_dirac_endpoint_jets.py",
    "scripts/derive_nsc_transmitting_boundary_binding.py",
)


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def payload(record):
    item = record["payload"]
    if sha(item["path"]) != item["sha256"]:raise ValueError("authenticated artifact differs")
    with np.load(ROOT/item["path"], allow_pickle=False) as f:
        return {key: np.array(f[key]) for key in f.files}


def calculate():
    domain, state = (json.loads((ROOT/path).read_text()) for path in INPUTS)
    traces, modes = payload(domain), payload(state)
    seed = state["cauchy_surfaces"]["seed"]
    for key in ("A", "magnetic_flux", "Omega", "zeta", "V_full"):
        if domain["locked_inputs"][key] != state["locked_inputs"][key]:raise ValueError("locked ledger differs")
    tol = 3e-11
    fields = ("N", "beta", "q_ADM", "r")
    residuals = {name: 0. for name in (
        "artifact_canonicalization", "normalized_graph_isometry",
        "graph_complement_orthogonality", "extension_hermiticity",
        "extension_graph_invariance", "physical_h_compression",
        "physical_vertex_compression", "CAR_distance_from_exact_one",
    )}
    car_errors, ranks, b_differences = [], [], []
    jet_differences = {field: 0. for field in fields}
    covariance_checks_pass = True
    samples = len(modes["frequency"])
    for index in range(samples):
        canonical = traces["PG_trace_to_seed"][index]@traces["seed_to_PG_trace"][index]
        e = np.vstack((canonical, canonical))/np.sqrt(2)
        f = np.vstack((canonical, -canonical))/np.sqrt(2)
        residuals["artifact_canonicalization"] = max(residuals["artifact_canonicalization"], float(np.max(abs(canonical-np.eye(2)))))
        residuals["normalized_graph_isometry"] = max(residuals["normalized_graph_isometry"], float(np.max(abs(e.conj().T@e-np.eye(2)))))
        residuals["graph_complement_orthogonality"] = max(residuals["graph_complement_orthogonality"], float(np.max(abs(e.conj().T@f))))
        car = e@e.conj().T  # Operator anticommutator of b=E a; NOT <b^dagger b>.
        error = float(np.linalg.norm(car-np.eye(4), 2))
        car_errors.append(error);ranks.append(int(np.linalg.matrix_rank(car, tol)))
        residuals["CAR_distance_from_exact_one"] = max(residuals["CAR_distance_from_exact_one"], abs(error-1))
        # Positivity of a supplied occupation covariance does not establish CAR.
        _covariance(e@modes["covariance_seed"][index]@e.conj().T)
        h, vertices = bulk_hamiltonian_and_vertices(
            modes["hx_seed"][index], modes["hy_seed"][index],
            modes["hz_seed"][index]*seed["a_parallel"],
            1., 0., seed["a_parallel"], seed["r"], seed["r"])
        # Algebraic witnesses only. The complementary block is not a new
        # physical sector, source, Hamiltonian, or term adopted by the model.
        lift_zero = e@h@e.conj().T                  # unphysical complement F=0
        lift_equal = lift_zero+f@h@f.conj().T       # complement chosen equal to h
        b_differences.append(float(np.max(abs(lift_zero[:2, 2:]-lift_equal[:2, 2:]))))
        for lift in (lift_zero, lift_equal):
            residuals["extension_hermiticity"] = max(residuals["extension_hermiticity"], float(np.max(abs(lift-lift.conj().T))))
            residuals["extension_graph_invariance"] = max(residuals["extension_graph_invariance"], float(np.max(abs(lift@e-e@h))))
            residuals["physical_h_compression"] = max(residuals["physical_h_compression"], float(np.max(abs(e.conj().T@lift@e-h))))
        for field, vertex in zip(fields, vertices):
            v0 = e@vertex@e.conj().T
            v1 = v0+f@vertex@f.conj().T
            jet_differences[field] = max(jet_differences[field], float(np.max(abs(v0[:2, 2:]-v1[:2, 2:]))))
            for v in (v0, v1):
                residuals["physical_vertex_compression"] = max(residuals["physical_vertex_compression"], float(np.max(abs(e.conj().T@v@e-vertex))))
    if any(x > tol for x in residuals.values()) or any(rank != 2 for rank in ranks):raise ArithmeticError("representation certificate failed")
    dimensions = domain["finite_domain_dimensions"]
    if sum(ranks) != dimensions["independent_coefficient_dimension"]:raise ValueError("domain rank differs")
    return {
        "schema": "NSC-HAMILTONIAN-TRACE-REPRESENTATION-v1",
        "status": "U OPEN: literal independent-trace embedding rejected; constrained lifts do not determine physical B",
        "source_hashes": {path: sha(path) for path in SOURCES},
        "input_hashes": {path: sha(path) for path in (*INPUTS, *DEPENDENCIES)},
        "authenticated_T_payload": domain["payload"],
        "locked_inputs": domain["locked_inputs"],
        "domain": "two complete T traces of one canonical seed coefficient vector; no bulk spatial partition is supplied",
        "operator_CAR_certificate": {
            "embedding": "E=(I,I)^T/sqrt(2) after the authenticated T inverse",
            "physical_relation": "{a,a^dagger}=I_d; b=E a implies {b,b^dagger}=E E^dagger",
            "required_by_independent_room_assembly": "I_(2d)",
            "independent_coefficients": sum(ranks),
            "proposed_independent_trace_slots": 4*samples,
            "rank_defect": 2*samples,
            "exact_spectral_norm_lower_bound": 1,
            "measured_minimum_spectral_norm_mismatch": min(car_errors),
            "measured_maximum_spectral_norm_mismatch": max(car_errors),
            "target_tolerance": tol,
            "literal_trace_pair_as_independent_rooms": "FAIL",
            "occupation_covariance_positivity_checks_pass": covariance_checks_pass,
            "positivity_does_not_change_operator_CAR": True,
            "not_a_no_go_for_true_orthogonal_bulk_partitions": True,
        },
        "constrained_extension_certificate": {
            "definition": "H_J=E h E^dagger+F J F^dagger; F=(I,-I)^T/sqrt(2); J arbitrary Hermitian on the redundant complement",
            "blocks": "H_pp=H_cc=(h+J)/2; B=(h-J)/2 in ideal canonical coordinates",
            "witnesses": ["J=0 gives B=h/2", "J=h gives B=0"],
            "fixed_data": "the physical single-field h and its vertices; T does not supply separately fixed independent H_parent and H_child",
            "redundant_diagonal_blocks_also_vary_between_witnesses": True,
            "physical_h_and_vertices_unchanged": True,
            "maximum_witness_B_difference": max(b_differences),
            "maximum_witness_dB_difference_by_field": jet_differences,
            "witnesses_adopted_as_physical_links": False,
            "complement_included_in_physical_spectral_trace": False,
            "meaning": "B and its derivatives depend on an unphysical extension unless a genuine room partition is specified",
        },
        "residuals": residuals, "certificate_tolerance": tol,
        "missing_physical_representation": {
            "required": "equal-time independent bulk subspaces and embeddings P_parent,P_child of the physical Dirac Hilbert space, compatible with T's trace condition",
            "conditions": "P_parent P_child=0; P_parent+P_child=I on the physical space; domains and common clock specified",
            "only_then": "B=P_parent H_D[g] P_child with derivatives of H_D and the actual projectors",
            "these_projectors_are_not_Ts_two_complete_trace_maps": True,
        },
        "physical_outputs": {"B": None, "dB_dgDelta": None, "full_EndpointBranchJets": None,
                             "Gamma_rest_boundary_derivative": None, "two_sided_Weyl_mismatch": None,
                             "physical_Vc": None, "stationary_history": None},
        "downstream": {"U": "OPEN", "V": "moving-surface jets not evaluated", "W": "OPEN", "X": "OPEN"},
        "scope": {"T_preserved": True, "new_physical_terms_or_copies": False,
                  "given_B_used": False, "parameters_refitted": False,
                  "physical_history_selected": False, "old_generators_rerun": False,
                  "metric_timestep_started": False, "optimizer_started": False,
                  "finite_stress_fabricated": False, "extended_nonexistence_claimed": False},
        "publication_gate": "physical U minimum not met; retain v0.22.0 PDF and curated release",
        "comparison": domain["comparison"],
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    group = p.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--write", action="store_true")
    args = p.parse_args();record = calculate();out = ROOT/OUTPUT
    if args.check:
        compare(json.loads(out.read_text()), record)
        print("Hamiltonian trace-representation certificate verified; physical U remains OPEN")
    elif args.write:
        if out.exists():raise FileExistsError("refusing overwrite")
        out.write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+"\n")
        print(record["status"])
    else:print(json.dumps(record,indent=2,sort_keys=True,allow_nan=False))


if __name__ == "__main__":main()
