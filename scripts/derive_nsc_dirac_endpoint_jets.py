#!/usr/bin/env python3
"""Verify known KS Hamiltonian-to-CTP jets; leave physical transmission OPEN."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from recursive_horizons.nsc_dirac_endpoint_jets import bulk_endpoint_jets
from recursive_horizons.nsc_landau_cauchy_isometry import KSCauchyHistory, propagate_blockwise
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_transmitting_boundary_history import EndpointVariation, FIELDS
from recursive_horizons.nsc_transmitting_ctp_variation import endpoint_ctp_pullback
from recursive_horizons.nsc_influence import influence
from derive_nsc_transmitting_boundary_binding import compare

INPUTS = {
    "state": "results/development/nsc-mode-resolved-cauchy-state.json",
    "isometry": "results/development/nsc-landau-cauchy-isometry.json",
    "endpoint": "results/development/nsc-weyl-endpoint-match.json",
    "differential": "results/development/nsc-transmitting-ctp-variation.json",
    "extended": "results/development/nsc-extended-action-completion-gate.json",
}
SOURCES = (
    "src/recursive_horizons/nsc_dirac_endpoint_jets.py",
    "scripts/derive_nsc_dirac_endpoint_jets.py",
    "docs/nsc-dirac-endpoint-jets.md",
    "docs/nsc-transmitting-boundary-remainder.md",
    "docs/nsc-physical-jet-extended-gate.md",
)
DEPENDENCIES = (
    "src/recursive_horizons/nsc_landau_cauchy_isometry.py",
    "src/recursive_horizons/nsc_transmitting_ctp_variation.py",
    "src/recursive_horizons/nsc_transmitting_boundary_history.py",
    "src/recursive_horizons/nsc_mode_resolved_cauchy_state.py",
    "src/recursive_horizons/nsc_influence.py",
    "scripts/derive_nsc_transmitting_boundary_binding.py",
)
OUTPUTS = ("nsc-dirac-endpoint-jets.json", "nsc-transmitting-boundary-remainder.json", "nsc-physical-jet-extended-gate.json")


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def load_payload(record):
    p = record["payload"]
    if sha(p["path"]) != p["sha256"]:
        raise ValueError("authenticated payload changed")
    with np.load(ROOT/p["path"], allow_pickle=False) as archive:
        return {key: np.array(archive[key]) for key in archive.files}


def calculate():
    records = {key: json.loads((ROOT/path).read_text()) for key, path in INPUTS.items()}
    state = records["state"]
    locked = records["differential"]["locked_inputs"]
    if any(locked[k] != state["locked_inputs"][k] for k in ("A", "magnetic_flux", "Omega", "zeta", "V_full")):
        raise ValueError("locked scale or coefficient differs")
    all_arrays, control = load_payload(state), load_payload(records["isometry"])
    index = np.linspace(0, len(control["short_time"])-1, 33, dtype=int)
    history = KSCauchyHistory(*(control[key][index] for key in
        ("short_time", "short_a", "short_r", "short_lapse", "short_shift")),
        label="authenticated 33-node diagnostic; not a physical history")
    # One fixed middle quadrature node per retained symmetry channel. No scan.
    representatives = np.array([int(c["sample_offset"])+int(c["sample_count"])//2 for c in state["channels"]])
    names = ("hx_seed", "hy_seed", "hz_seed", "covariance_seed")
    arrays = {key: all_arrays[key][representatives] for key in names}
    seed = state["cauchy_surfaces"]["seed"]
    endpoint = records["endpoint"]["endpoint"]
    basis = EndpointVariation(endpoint["basis"]["id"], tuple(tuple(endpoint["local_components"][f]) for f in FIELDS))
    if endpoint["basis"]["history_nodes"] != len(history.time):
        raise ValueError("endpoint node basis changed")
    computed = bulk_endpoint_jets(arrays, history, seed_a=seed["a_parallel"], seed_r=seed["r"], basis_id=basis.basis_id)
    base = propagate_blockwise(arrays, history, seed_a=seed["a_parallel"], seed_r=seed["r"])
    u = computed["unitary"]
    contractions = [endpoint_ctp_pullback(c, v, v, basis, j) for c, v, j in
                    zip(arrays["covariance_seed"], u, computed["branch_jets"])]
    numeric_jets = np.empty_like(computed["dU_dendpoint"])
    numeric_action = np.empty((len(representatives), 4, 2), dtype=complex)
    epsilon = 1e-5
    for field, attribute in enumerate(("lapse", "shift", "a_parallel", "radius")):
        for side, node in enumerate((0, len(history.time)-1)):
            variants = []
            for sign in (1, -1):
                varied = KSCauchyHistory(*(np.array(getattr(history, k), copy=True) for k in
                                          ("time", "a_parallel", "radius", "lapse", "shift")),
                                          label="new-derivative perturbation of frozen control")
                getattr(varied, attribute)[node] += sign*epsilon/2
                variants.append(propagate_blockwise(arrays, varied, seed_a=seed["a_parallel"], seed_r=seed["r"])["unitary"])
            numeric_jets[:, field, side] = (variants[0]-variants[1])/epsilon
            for mode, c in enumerate(arrays["covariance_seed"]):
                forward = influence(c, variants[0][mode], variants[1][mode])["principal_action"]
                backward = influence(c, variants[1][mode], variants[0][mode])["principal_action"]
                numeric_action[mode, field, side] = (forward-backward)/(2*epsilon)
    action = np.array([[row["gaussian_endpoint_covector"][field] for field in FIELDS] for row in contractions])
    residuals = {
        "unitarity": computed["unitarity_residual"],
        "step_tangent": computed["step_tangent_residual"],
        "endpoint_tangent": computed["final_tangent_residual"],
        "U_against_owned_propagator": float(np.max(np.abs(u-base["unitary"]))),
        "dU_against_perturbed_owned_propagator": float(np.max(np.abs(computed["dU_dendpoint"]-numeric_jets))),
        "endpoint_contraction_against_owned_action": float(np.max(np.abs(action-numeric_action))),
        "algebra_tolerance": 3e-11, "finite_difference_tolerance": 3e-8,
    }
    for k in ("unitarity", "step_tangent", "endpoint_tangent", "U_against_owned_propagator"):
        if residuals[k] > residuals["algebra_tolerance"]: raise ArithmeticError(k)
    for k in ("dU_against_perturbed_owned_propagator", "endpoint_contraction_against_owned_action"):
        if residuals[k] > residuals["finite_difference_tolerance"]: raise ArithmeticError(k)
    payload_arrays = {"representative_sample_indices": representatives, "time": history.time,
                      "lapse": history.lapse, "shift": history.shift, "a_parallel": history.a_parallel,
                      "radius": history.radius, "bulk_unitary": u,
                      "bulk_dU_dendpoint": computed["dU_dendpoint"],
                      "bulk_gaussian_endpoint_covector": action}
    data = deterministic_npz_bytes(payload_arrays)
    digest = hashlib.sha256(data).hexdigest()
    payload = {"path": f"results/development/artifacts/nsc-dirac-endpoint-jets.{digest}.npz", "sha256": digest,
               "bytes": len(data), "physical": False, "description": "partial diagonal bulk jets on authenticated diagnostic control; no transmitting jets"}
    common = {
        "source_hashes": {p: sha(p) for p in SOURCES},
        "input_hashes": {p: sha(p) for p in (*INPUTS.values(), *DEPENDENCIES)},
        "locked_inputs": locked, "comparison": records["differential"]["comparison"],
        "scope": {"parameters_refitted": False, "new_physical_term": False, "old_generator_rerun": False,
                  "metric_timestep_started": False, "vacuum_stress_integration_run": False,
                  "physical_duration_selected": False, "finite_stress_fabricated": False},
    }
    j = dict(common, schema="NSC-DIRAC-ENDPOINT-JETS-v1", status="OPEN", payload=payload,
        bulk_jet_component={"status": "PASS", "residuals": residuals,
            "hamiltonian": "H=N(hx_seed sigma1+hy_seed*r_seed/r sigma2+k/a sigma3)-beta*k I; k=hz_seed*a_seed",
            "propagation": "same midpoint ordered product as nsc_landau_cauchy_isometry",
            "derivative": "Frechet derivative of each matrix exponential, transported by the product rule",
            "branch_convention": "dUplus/dgDelta=+dU/dg/2; dUminus/dgDelta=-dU/dg/2",
            "controls": {"channels": len(representatives), "blocks_per_channel": 1, "representative_indices": representatives.tolist(),
                         "history_nodes": len(history.time), "source": records["isometry"]["payload"],
                         "purpose": "verify new derivative transport only; old homogeneous class remains excluded",
                         "raw_metric_finite_difference_step": epsilon}},
        physical_jets={"status": "OPEN", "jets": None, "residual": None,
                       "dB_dgDelta": None, "transmitting_embedding_and_normal_derivative": None,
                       "reason": "known diagonal bulk vertices do not supply the transmitting link, frequency-mixing, or hypersurface derivative"})
    k = dict(common, schema="NSC-TRANSMITTING-BOUNDARY-REMAINDER-v1", status="OPEN",
        endpoint_basis=endpoint["basis"], locked_local_components=endpoint["local_components"],
        locked_extraction_residual=endpoint["binding_residual"], preserved_Weyl_value=endpoint["local_maximum_absolute"],
        existing_EndpointBranchJets_used=True,
        partial_bulk_contraction={"record": "nsc-dirac-endpoint-jets.json", "status": "diagnostic derivative only; not an allocated transmitting boundary term"},
        physical_boundary_remainder=None, full_physical_jets=None,
        two_sided_mismatch={"status": "OPEN", "value": None, "tolerance": 3e-11},
        normal_jet_and_boundary_allocation=endpoint["independent_first_jet_derivative"],
        reason="neither a transmitting boundary action density nor its endpoint/normal-jet variation is supplied; the bulk control cannot cancel the Weyl vector")
    previous = records["extended"]
    if not previous["composition"]["homogeneous_nonexistence_preserved"]: raise ValueError("old regression changed")
    gate = dict(common, schema="NSC-PHYSICAL-JET-EXTENDED-GATE-v1", status="OPEN",
        steps={"J": "bulk derivative PASS / physical jets OPEN", "K": "OPEN", "L": "OPEN", "M": "OPEN"},
        new_component_residuals=residuals,
        missing_residuals=["transmitting dB/dgDelta and embedding/normal derivatives in dU/dgDelta",
                           "same-action boundary remainder on the existing endpoint and normal-jet basis",
                           "physical kernel selection and full history stationarity using those derivatives"],
        inherited_prerequisites=previous["prerequisites"],
        locked_Gaussian_differential_checks=records["differential"]["selection"]["residuals"],
        homogeneous_nonexistence_preserved=True, optimizer_started=False,
        existence_claimed=False, extended_nonexistence_claimed=False,
        physical_kernel=None, finite_stress=None, nulls=None, updated_constraints=None,
        metric_timestep_started=False, coupled_evolution_reopened=False)
    return dict(zip(OUTPUTS, (j, k, gate))), data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--write", action="store_true")
    args = parser.parse_args()
    records, payload_data = calculate()
    payload = records[OUTPUTS[0]]["payload"]
    folder = ROOT/"results/development"
    if args.check:
        for name, record in records.items():compare(json.loads((folder/name).read_text()), record)
        if (ROOT/payload["path"]).read_bytes() != payload_data:raise AssertionError("partial jet artifact differs")
        print("J/K records and extended OPEN gate verified; partial bulk jet artifact identical")
    elif args.write:
        targets = [folder/name for name in records]+[ROOT/payload["path"]]
        if any(p.exists() for p in targets):raise FileExistsError("refusing to overwrite records or artifact")
        (ROOT/payload["path"]).write_bytes(payload_data)
        for name, record in records.items():(folder/name).write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+"\n")
        print("Wrote three OPEN records with verified partial bulk jets")
    else:print(json.dumps(records,indent=2,sort_keys=True,allow_nan=False))


if __name__ == "__main__":main()
