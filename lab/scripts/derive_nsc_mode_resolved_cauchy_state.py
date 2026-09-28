#!/usr/bin/env python3
"""Persist the actual retained NSC Gaussian mode state and its seed reducer."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_all_start_methods, get_context
import hashlib
import json
from math import cosh, exp, sinh
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

from recursive_horizons.nsc_charged_ctp_neck import ChargedCTPNeckConfig
from recursive_horizons.nsc_compact_ctp_neck import CompactCTPConfig
from recursive_horizons.nsc_mode_resolved_cauchy_state import (
    ModeResolvedCauchyState,
    array_digest,
    deterministic_npz_bytes,
    extract_charged_angular,
    extract_lll,
    extract_positive_compact,
    pack_channels,
    reconstruct_seed_tensor,
)


OUTPUT = ROOT/"results/development/nsc-mode-resolved-cauchy-state.json"
ARTIFACT_DIRECTORY = ROOT/"results/development/artifacts"
INPUTS = (
    "results/development/nsc-coupled-ctp-metric-evolution.json",
    "results/development/nsc-constraint-complete-neck.json",
    "results/development/charged-compact-ctp-completion.json",
    "results/development/charged-ctp-neck-source.json",
    "results/development/scale-binding.json",
)
SOURCES = (
    "src/recursive_horizons/nsc_mode_resolved_cauchy_state.py",
    "scripts/derive_nsc_mode_resolved_cauchy_state.py",
    "docs/nsc-mode-resolved-cauchy-state.md",
    "src/recursive_horizons/nsc_charged_ctp_neck.py",
    "src/recursive_horizons/nsc_compact_ctp_neck.py",
    "src/recursive_horizons/nsc_angular_stress.py",
    "src/recursive_horizons/nsc_horizon_source.py",
    "src/recursive_horizons/nsc_unruh_state.py",
)


def hashes(paths):
    return {
        path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
        for path in paths
    }


def compare(expected, actual, path="$", atol=2e-6, rtol=2e-5):
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


def _configs(completion, angular):
    ac = angular["runs"]["base"]["config"]
    cc = completion["runs"]["base"]["config"]
    charged = ChargedCTPNeckConfig(
        magnetic_flux=ac["magnetic_flux"], omega=ac["omega"],
        horizon_rho=ac["horizon_rho"], surface_gravity=ac["surface_gravity"],
        angular_levels=ac["angular_levels"],
        points_per_frequency_interval=ac["points_per_frequency_interval"],
        frequency_edges=tuple(ac["frequency_edges"]),
        evolution_steps=ac["evolution_steps"], phase_cutoff=ac["phase_cutoff"],
        horizon_offset=ac["horizon_offset"],
        reflection_tolerance=ac["reflection_tolerance"],
    )
    compact = CompactCTPConfig(
        magnetic_flux=cc["magnetic_flux"], omega=cc["omega"],
        cutoff=cc["cutoff"], horizon_rho=cc["horizon_rho"],
        surface_gravity=cc["surface_gravity"],
        compact_levels=tuple(cc["compact_levels"]),
        angular_levels=cc["angular_levels"],
        frequency_points=cc["frequency_points"],
        frequency_max=cc["frequency_max"], evolution_steps=cc["evolution_steps"],
        horizon_offset=cc["horizon_offset"],
        scattering_tolerance=cc["scattering_tolerance"],
        outer_floor=cc["outer_floor"],
    )
    return charged, compact


def _extract(charged, compact):
    context = get_context("fork") if "fork" in get_all_start_methods() else get_context()
    with ProcessPoolExecutor(max_workers=2, mp_context=context) as pool:
        angular_future = pool.submit(extract_charged_angular, charged)
        compact_future = pool.submit(extract_positive_compact, compact)
        angular_result, angular_channels, _jets = angular_future.result()
        compact_result, compact_channels = compact_future.result()
    lll = extract_lll(charged)
    arrays, channels = pack_channels([lll, *angular_channels, *compact_channels])
    return arrays, channels, angular_result, compact_result


def _array_manifest(arrays):
    return {
        name: {
            "dtype": value.dtype.str,
            "shape": list(value.shape),
            "sha256": array_digest(name, value),
        }
        for name, value in sorted(arrays.items())
    }


def _local_spin_current_sqrt(eta):
    # B_L^(1/2)=exp(-eta sigma2/2), in [upper,lower] order.
    c, s = cosh(eta/2), sinh(eta/2)
    return [
        [[c, 0.0], [0.0, s]],
        [[0.0, -s], [c, 0.0]],
    ]


def _record_from_payload(path, arrays, channels, records, generated=None):
    coupled, selected, completion, angular, scale = records
    state = ModeResolvedCauchyState(arrays, tuple(channels))
    validation = state.validate()
    compact_local = angular["runs"]["base"]["compact_local_source"]["child_frame"]
    reconstructed, details = reconstruct_seed_tensor(arrays, channels, compact_local)
    expected = completion["completed_tensor"]
    fields = ("rho", "T01", "p_parallel", "p_sphere", "parent_Killing_power")
    residuals = {name: reconstructed[name]-expected[name] for name in fields}
    maximum = max(abs(value) for value in residuals.values())
    reconstruction_tolerance = 3e-6
    if maximum >= reconstruction_tolerance:
        raise ArithmeticError(f"seed moment reconstruction failed: {residuals}")
    eta = selected["source_selected_frame"]["rapidity"]
    local_sqrt = _local_spin_current_sqrt(eta)
    metric = coupled["prestep_metric_vertex"]
    payload_bytes = path.read_bytes()
    digest = hashlib.sha256(payload_bytes).hexdigest()
    locked = selected["locked_inputs"]
    channel_counts = {}
    for row in channels:
        key = str(row["family"])
        channel_counts[key] = channel_counts.get(key, 0)+1
    return {
        "schema": "NSC-MODE-RESOLVED-CAUCHY-STATE-v1",
        "artifact_id": "NSC-MODE-RESOLVED-CAUCHY-STATE",
        "status": (
            "PARTIAL PASS: all 1904 seed covariance blocks are serialized and "
            "reconstruct the unit-radius tensor; the physical Landau Cauchy "
            "isometry and new-history renormalized stress remain undefined"
        ),
        "source_hashes": hashes(SOURCES),
        "input_hashes": hashes(INPUTS),
        "locked_inputs": {
            "magnetic_flux": locked["q"], "Omega": locked["Omega"],
            "zeta": locked["zeta"], "A": locked["A"],
            "V_full": locked["V_full"],
            "occupations_refitted": False,
        },
        "basis": {
            "canonical_half_density": "chi=r*sqrt(a_parallel)*psi at N=1",
            "index_order": "channel, frequency, spinor; upper before lower",
            "Pauli_axes": "H=-m_j*sigma1+(lambda_n/r)*sigma2+(k/a_parallel)*sigma3",
            "covariance_basis": "canonical transported seed basis used by the immutable source kernels",
        },
        "cauchy_surfaces": {
            "seed": {"r": 1.0, "a_parallel": selected["source_selected_initial_geometry"]["a_parallel"]},
            "selected": {
                "r": selected["source_selected_initial_geometry"]["r"],
                "a_parallel": selected["source_selected_initial_geometry"]["a_parallel"],
                "H_parallel": selected["source_selected_initial_geometry"]["H_parallel"],
                "H_sphere": selected["source_selected_initial_geometry"]["H_sphere"],
                "a_parallel_and_H_parallel": "seed retention",
            },
        },
        "channels": channels,
        "channel_summary": {
            "count": len(channels), "sample_count": validation["samples"],
            "family_codes": {"0": "LLL", "1": "massive angular", "2": "positive compact"},
            "counts_by_family": channel_counts,
        },
        "payload": {
            "path": str(path.relative_to(ROOT)), "sha256": digest,
            "bytes": len(payload_bytes), "format": "deterministic uncompressed NPZ",
            "allow_pickle": False, "arrays": _array_manifest(arrays),
        },
        "seed_validation": {
            **validation,
            "reconstructed_tensor": reconstructed,
            "expected_tensor": {name: expected[name] for name in fields},
            "moment_residuals": residuals,
            "maximum_absolute_moment_residual": maximum,
            "declared_reconstruction_tolerance": reconstruction_tolerance,
            "allocation_details": details,
        },
        "slice_map": {
            "Landau_v": selected["source_selected_frame"]["velocity_old_child_frame"],
            "Landau_eta": eta,
            "local_spin_current_sqrt_complex_pairs": local_sqrt,
            "local_matrix_is_full_Cauchy_map": False,
            "required_full_map": "U_L0=J_L Res_L E_g Res_0^-1 J_0^-1",
            "required_isometry": "U_L0^dagger U_L0=I in quadrature-normalized half-density basis",
            "frequency_mixing_allowed": True,
            "serialized": False,
            "reason": (
                "endpoint rapidity does not specify the Landau hypersurface "
                "embedding or intervening Dirac metric/gauge history"
            ),
        },
        "renormalization": {
            "LLL": "analytic conformal local term plus serialized transparent covariance",
            "massive_angular": "serialized physical C with E2+E4/P2+P4 seed subtraction",
            "positive_compact": "serialized physical C and fourth-order superadiabatic seed reference",
            "compact_Wilsonian_local": "existing completion allocation included exactly once",
            "instantaneous_vacuum_subtraction": False,
            "general_KS_order4_reference_history_implemented": False,
            "general_KS_local_induced_history_implemented": False,
        },
        "selected_geometry_response": {
            "angular_Hamiltonian_fractional_change": metric["angular_Hamiltonian_fractional_change"],
            "fixed_covariance_density_vertex": metric["fixed_covariance_density_log_radius_vertex"],
            "direction_check_nonzero": metric["fixed_covariance_density_log_radius_vertex"] < 0,
            "finite_renormalized_stress_evaluated": False,
            "linearized_diagnostic_used_as_stress": False,
        },
        "api_binding": {
            "loader": "ModeResolvedCauchyState.load(path, channels)",
            "consumer_shape": "blockwise weighted 2x2 channel iterator",
            "dense_3808_by_3808_embedding_used": False,
            "CausalCommonFunctional_direct_hook_complete": False,
            "ADM_metric_step_authorized": False,
        },
        "gate": {
            "physical_C_eigenvalues_in_unit_interval": (
                validation["minimum_eigenvalue"] >= -2e-10
                and validation["maximum_eigenvalue"] <= 1+2e-10
            ),
            "seed_serialization_pass": True,
            "old_surface_moment_reconstruction_pass": maximum < reconstruction_tolerance,
            "selected_surface_direction_check_pass": metric["fixed_covariance_density_log_radius_vertex"] < 0,
            "Landau_slice_CAR_isometry_pass": False,
            "selected_surface_renormalized_stress_pass": False,
            "full_contract_pass": False,
            "coupled_evolution_reopened": False,
            "next_owner": "physical Dirac Cauchy isometry plus general-KS order-4/local history providers",
        },
        "extraction": {
            "one_authorized_base_run": generated is not None,
            "immutable_source_files_modified": False,
            "adapter": "profiled return-frame capture of exact final kernel vectors",
            "old_generator_invoked_as_standalone": False,
        },
        "scope": {
            "new_source_channel_or_occupation_fit": False,
            "A_q_Omega_zeta_or_V_full_changed": False,
            "four_moment_surrogate_or_w_of_a": False,
            "metric_time_step_started": False,
        },
        "comparison": {
            "fields": "all manifest fields and payload arrays",
            "exact": "schema, hashes, dtypes, shapes, offsets, labels and payload bytes",
            "float_atol": 2e-6, "float_rtol": 2e-5,
            "exceptions": ["LLL finite quadrature replaces its analytic integral within 2e-12"],
        },
    }


def _load_records():
    return tuple(json.loads((ROOT/path).read_text()) for path in INPUTS)


def _load_payload(record):
    relative = Path(record["payload"]["path"])
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("unsafe payload path")
    path = ROOT/relative
    if hashlib.sha256(path.read_bytes()).hexdigest() != record["payload"]["sha256"]:
        raise ValueError("payload digest mismatch")
    state = ModeResolvedCauchyState.load(path, record["channels"])
    return path, state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = _load_records()
    if args.check:
        expected = json.loads(OUTPUT.read_text())
        if expected["source_hashes"] != hashes(SOURCES) or expected["input_hashes"] != hashes(INPUTS):
            raise AssertionError("source or input hash changed")
        path, state = _load_payload(expected)
        actual = _record_from_payload(path, state.arrays, list(state.channels), records)
        # Extraction metadata records that creation performed the authorized run.
        actual["extraction"]["one_authorized_base_run"] = expected["extraction"]["one_authorized_base_run"]
        compare(expected, actual)
        print("mode-resolved seed state, payload and old-surface moments verified without source rerun")
        return
    output = args.output or OUTPUT
    if output.exists():
        raise FileExistsError("refusing to overwrite recorded evidence")
    completion, angular = records[2], records[3]
    charged_config, compact_config = _configs(completion, angular)
    arrays, channels, angular_result, compact_result = _extract(charged_config, compact_config)
    # The legacy reducers ran on the same in-memory arrays; require their stored totals.
    compare(angular["runs"]["base"]["massive_angular_source"], angular_result, "$/angular")
    compare(completion["runs"]["base"], compact_result, "$/compact")
    payload_bytes = deterministic_npz_bytes(arrays)
    digest = hashlib.sha256(payload_bytes).hexdigest()
    ARTIFACT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    payload_path = ARTIFACT_DIRECTORY/f"nsc-mode-resolved-cauchy-state.{digest}.npz"
    if payload_path.exists() and payload_path.read_bytes() != payload_bytes:
        raise FileExistsError("content-addressed payload collision")
    payload_path.write_bytes(payload_bytes)
    record = _record_from_payload(
        payload_path, arrays, channels, records,
        generated={"angular": angular_result, "compact": compact_result},
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(record, indent=2, sort_keys=True)+"\n")
    print(output)
    print(payload_path)


if __name__ == "__main__":
    main()
