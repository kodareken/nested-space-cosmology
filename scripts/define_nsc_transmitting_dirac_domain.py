#!/usr/bin/env python3
"""Define and verify the retained NSC Dirac trace domain at the actual seed seam."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from recursive_horizons.nsc_transmitting_dirac_domain import (
    TransmittingDiracSeamDomain, DomainRepresentationError,
    MODE_TO_CURRENT, I2, S1, S2, S3, require_transmitting_domain,
)
from recursive_horizons.nsc_causal_common import LinkHistory
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from derive_nsc_transmitting_boundary_binding import compare

OUTPUT = "results/development/nsc-transmitting-dirac-domain.json"
INPUTS = {
    "state": "results/development/nsc-mode-resolved-cauchy-state.json",
    "frame": "results/development/adm-neck-source-map.json",
    "previous": "results/development/nsc-physical-jet-extended-gate.json",
    "interface": "results/development/nsc-tilted-landau-interface.json",
}
SOURCES = (
    "src/recursive_horizons/nsc_transmitting_dirac_domain.py",
    "scripts/define_nsc_transmitting_dirac_domain.py",
    "tests/test_nsc_transmitting_dirac_domain.py",
    "docs/nsc-transmitting-dirac-domain.md",
)
DEPENDENCIES = (
    "src/recursive_horizons/nsc_causal_common.py",
    "src/recursive_horizons/nsc_chiral_boundary.py",
    "src/recursive_horizons/nsc_lorentzian.py",
    "src/recursive_horizons/nsc_mode_resolved_cauchy_state.py",
    "scripts/map_nsc_adm_neck_source.py",
    "scripts/derive_nsc_transmitting_boundary_binding.py",
)


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def calculate():
    records = {key: json.loads((ROOT/path).read_text()) for key, path in INPUTS.items()}
    state, frame, previous = (records[k] for k in ("state", "frame", "previous"))
    for key in ("A", "magnetic_flux", "Omega", "zeta", "V_full"):
        if previous["locked_inputs"][key] != state["locked_inputs"][key]:
            raise ValueError("locked ledger differs: "+key)
    p = state["payload"]
    if sha(p["path"]) != p["sha256"]: raise ValueError("stored state payload differs")
    with np.load(ROOT/p["path"], allow_pickle=False) as f:
        arrays = {key: np.array(f[key]) for key in f.files}
    geometry = frame["geometry"]
    if geometry["rho_coordinate"] != 0.: raise ValueError("source seam is not rho=0")
    domain = TransmittingDiracSeamDomain(geometry["N"], geometry["q"], geometry["beta"], geometry["r"])
    require_transmitting_domain(domain)
    domain.bind_seed(state)
    normal = domain.normal_geometry()
    seed = state["cauchy_surfaces"]["seed"]
    rotation = MODE_TO_CURRENT
    # The representation rotation is on the KS modes. The Lorentz coframe
    # conversion is a separate same-surface trace map, not a PG time generator.
    rotated_h = rotation@(arrays["hx_seed"][:, None, None]*S1
                          +arrays["hy_seed"][:, None, None]*S2
                          +arrays["hz_seed"][:, None, None]*S3)@rotation.conj().T
    expected_h = arrays["hx_seed"][:, None, None]*S1-arrays["hy_seed"][:, None, None]*S3+arrays["hz_seed"][:, None, None]*S2
    channels = []
    maps, inverses, grams, covariances = [], [], [], []
    for channel in state["channels"]:
        start = channel["sample_offset"];stop = start+channel["sample_count"]
        result = domain.evaluate_channel(arrays["quadrature_weight"][start:stop], arrays["covariance_seed"][start:stop])
        channels.append({"channel": channel, "residuals": result["residuals"],
                         "independent_coefficients": result["independent_coefficients"],
                         "doubled_trace_dimension": result["doubled_trace_dimension"],
                         "matching_constraint_rank": result["matching_constraint_rank"],
                         "allowed_graph_dimension": result["allowed_graph_dimension"]})
        for collection, key in ((maps, "trace_map"), (inverses, "trace_inverse"),
                                (grams, "trace_metric"), (covariances, "coordinate_trace_covariance")):
            collection.append(result[key])
    residuals = {key: max(row["residuals"][key] for row in channels) for key in channels[0]["residuals"]}
    residuals.update({
        "spin_rotation_unitarity": float(np.max(np.abs(rotation.conj().T@rotation-I2))),
        "current_intertwining": float(np.max(np.abs(rotation@S3@rotation.conj().T-S2))),
        "rotated_mode_H": float(np.max(np.abs(rotated_h-expected_h))),
        "stored_axial_momentum_plus_frequency": float(np.max(np.abs(seed["a_parallel"]*arrays["hz_seed"]+arrays["frequency"]))),
        "normal_unit": normal["normal_unit_residual"],
        "normal_tangent": normal["normal_tangent_residual"],
        "induced_metric": normal["induced_metric_residual"],
        "seed_a_match": abs(domain.induced_axial_scale-seed["a_parallel"]),
    })
    tolerance = 3e-11
    if any(value > tolerance for value in residuals.values()):raise ArithmeticError("domain definition residual failed")
    rejected = {}
    for label, call in (
        ("given_B_only", lambda: require_transmitting_domain(LinkHistory(np.zeros((1, 1, 1))))),
        ("instantaneous_B_from_temporal_sewing", domain.instantaneous_link_history),
        ("copy_seed_to_rstar", lambda: TransmittingDiracSeamDomain(domain.lapse, domain.radial_scale,
            domain.shift, state["cauchy_surfaces"]["selected"]["r"]).bind_seed(state)),
    ):
        try: call()
        except DomainRepresentationError as error: rejected[label] = str(error)
        else: raise AssertionError("invalid domain path accepted: "+label)
    payload_arrays = {
        "mode_to_current": rotation, "normal_spin_metric": normal["normal_spin_metric"],
        "frequency": arrays["frequency"], "spatial_momentum_on_seam": seed["a_parallel"]*arrays["hz_seed"],
        "quadrature_weight": arrays["quadrature_weight"], "channel_index": arrays["channel_index"],
        "seed_to_PG_trace": np.concatenate(maps), "PG_trace_to_seed": np.concatenate(inverses),
        "trace_metric": np.concatenate(grams), "coordinate_trace_covariance": np.concatenate(covariances),
    }
    data = deterministic_npz_bytes(payload_arrays)
    digest = hashlib.sha256(data).hexdigest()
    record = {
        "schema": "NSC-TRANSMITTING-DIRAC-DOMAIN-v1",
        "status": "DEFINITION PASS: finite retained same-surface seed-seam domain; dynamical link and full closure OPEN",
        "source_hashes": {path: sha(path) for path in SOURCES},
        "input_hashes": {path: sha(path) for path in (*INPUTS.values(), *DEPENDENCIES)},
        "locked_inputs": previous["locked_inputs"],
        "payload": {"path": f"results/development/artifacts/nsc-transmitting-dirac-domain.{digest}.npz",
                    "sha256": digest, "bytes": len(data), "description": "same-surface spectral trace maps and coordinate covariances; not U_L0 or a physical stress"},
        "domain": {
            "surface": "Sigma0: rho=0 with coordinates (tau,theta,phi); tau=z on the seam",
            "surface_type": "spacelike three-surface; not the two-sphere rho=0 intersect tau=constant",
            "in": "retained reduced charged angular/compact channels; one Dirac field restricted from rho>0 and rho<0; common gauge and spin frame; smooth shared induced metric",
            "out": "continuum spectral completeness, four-spinor reconstruction of unresolved degeneracies, null surfaces, a moving neck or tilted SigmaL, physical kernel selection and independent simultaneous room Hamiltonians",
            "transmission_condition": "psi_parent|Sigma0=psi_child|Sigma0; same physical trace from both sides",
            "parent_outward": "future normal toward decreasing rho",
            "child_outward": "opposite past normal toward increasing rho",
            "no_independent_particle_copies_added": True,
            "compact_reflecting_endpoints_unchanged": True,
            "finite_projection_definition_not_continuum_wellposedness_proof": True,
        },
        "geometry": {
            "N": domain.lapse, "q_PG": domain.radial_scale, "beta": domain.shift, "r": domain.radius,
            "induced_a": domain.induced_axial_scale, "F_gtautau": float(normal["orbit_metric"][0, 0]),
            "coframe_rapidity": domain.coframe_rapidity,
            "coframe_rapidity_is_not_Landau_eta": True,
            "future_normal_tau_rho": normal["future_normal"].tolist(),
            "embedding_rule": "fixed coordinate cut X_Sigma0(tau,theta,phi)=(tau,0,theta,phi) in the inherited chart",
            "normal_rule": "aSigma=sqrt(q^2 beta^2-N^2); n=(q beta/(N aSigma),-aSigma/(N q)); future normal for q beta>N",
            "moving_embedding_or_history_selected": False,
        },
        "channel_map": {
            "stored_current": "sigma3 in the reduced KS mode representation",
            "common_current": "sigma2 after R=(I+i sigma1)/sqrt2",
            "operator_rotation": "hx sigma1+hy sigma2+hz sigma3 -> hx sigma1-hy sigma3+hz sigma2",
            "normal_current_metric": "J=(q beta I-N sigma2)/aSigma=exp(-xi sigma2)",
            "trace_map": "T_j=(r sqrt(aSigma*w_j))^-1 exp(+xi sigma2/2) R",
            "trace_metric": "G_j=aSigma*r^2*w_j*J",
            "canonical_identity": "T_j^dagger G_j T_j=I",
            "sewing_graph": "E_j a_j=(T_j a_j,T_j a_j); Green form diag(+G_j,-G_j)",
            "canonical_CAR_space": "one independent canonical a_j per retained node; the two traces are redundant, not a doubled CAR Hilbert space",
            "frequency_convention": "stored parent-neck omega; k=-omega on Sigma0; weights and occupations unchanged; no second Omega rescaling",
            "multiplicity_convention": "existing channel labels, degeneracies, copy counts and source factors preserved; no new angular-state basis or sector selected",
            "channels": channels,
        },
        "residuals": residuals, "tolerance": tolerance,
        "rejected_paths": rejected,
        "finite_domain_dimensions": {
            "channels": len(channels), "nodes": len(arrays["frequency"]),
            "independent_coefficient_dimension": sum(c["independent_coefficients"] for c in channels),
            "doubled_trace_dimension": sum(c["doubled_trace_dimension"] for c in channels),
            "matching_rank": sum(c["matching_constraint_rank"] for c in channels),
            "maximal_isotropic_graph_dimension": sum(c["allowed_graph_dimension"] for c in channels),
        },
        "substrate": {
            "channel_to_seed_boundary": "DEFINED and evaluated in the finite retained domain",
            "fixed_seam_embedding_and_normal": "DEFINED from the inherited geometry; no dynamical selection",
            "Hamiltonian_B_g_XSigma": "OPEN: temporal trace sewing does not specify an independent simultaneous parent/child Hamiltonian partition",
            "Gamma_rest_metric_and_normal_jet_dependence": "OPEN: not defined by fermionic trace matching",
        },
        "steps": {"T": "PASS in finite retained seed-seam domain", "U": "OPEN", "V": "fixed-seam normal defined; moving/selected-surface jets not evaluated", "W": "OPEN", "X": "OPEN"},
        "preserved": {"homogeneous_nonexistence": previous["homogeneous_nonexistence_preserved"],
                      "Gaussian_differential_checks": previous["locked_Gaussian_differential_checks"],
                      "bulk_jet_checks": previous["new_component_residuals"],
                      "Weyl_eight_coefficients_and_93_54264532195464": "unchanged; physical two-sided mismatch not evaluated"},
        "scope": {"new_physical_interaction": False, "parameters_refitted": False,
                  "physical_Vc_or_U_L0_selected": False, "coordinate_boost_used_as_propagation": False,
                  "source_tensor_copied_to_rstar": False, "finite_stress": None, "nulls": None,
                  "physical_two_sided_mismatch": None, "stationary_history": None,
                  "optimizer_started": False, "metric_timestep_started": False, "old_generators_rerun": False},
        "comparison": previous["comparison"],
    }
    return record, data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--write", action="store_true")
    args = parser.parse_args()
    record, data = calculate();output = ROOT/OUTPUT;artifact = ROOT/record["payload"]["path"]
    if args.check:
        compare(json.loads(output.read_text()), record)
        if artifact.read_bytes() != data:raise AssertionError("domain artifact differs")
        print("Transmitting seed-seam domain verified: all fields and trace-map artifact agree")
    elif args.write:
        if output.exists() or artifact.exists():raise FileExistsError("refusing overwrite")
        artifact.write_bytes(data);output.write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+"\n")
        print(record["status"])
    else:print(json.dumps(record,indent=2,sort_keys=True,allow_nan=False))


if __name__ == "__main__":main()
