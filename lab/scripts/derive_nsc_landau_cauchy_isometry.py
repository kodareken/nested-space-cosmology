#!/usr/bin/env python3
"""Evaluate the Landau Cauchy-isometry history-selection gate."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))

from recursive_horizons.nsc_landau_cauchy_isometry import (
    endpoint_control_history,
    propagate_blockwise,
)
from recursive_horizons.nsc_mode_resolved_cauchy_state import (
    ModeResolvedCauchyState,
    array_digest,
    deterministic_npz_bytes,
)


OUTPUT = ROOT/"results/development/nsc-landau-cauchy-isometry.json"
ARTIFACT_DIRECTORY = ROOT/"results/development/artifacts"
INPUTS = (
    "results/development/nsc-mode-resolved-cauchy-state.json",
    "results/development/nsc-constraint-complete-neck.json",
    "results/development/nsc-coupled-ctp-metric-evolution.json",
    "results/development/adm-source-constraints.json",
)
SOURCES = (
    "src/recursive_horizons/nsc_landau_cauchy_isometry.py",
    "scripts/derive_nsc_landau_cauchy_isometry.py",
    "docs/nsc-landau-cauchy-isometry.md",
)


def hashes(paths):
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}


def compare(expected, actual, path="$", atol=3e-12, rtol=3e-12):
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


def _load_inputs():
    return tuple(json.loads((ROOT/path).read_text()) for path in INPUTS)


def _load_seed(manifest):
    path = ROOT/manifest["payload"]["path"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != manifest["payload"]["sha256"]:
        raise ValueError("seed payload hash mismatch")
    return ModeResolvedCauchyState.load(path, manifest["channels"])


def _history_dict(history):
    return {
        "label": history.label, "duration": float(history.time[-1]),
        "points": len(history.time),
        "endpoints": {
            "a_parallel": [float(history.a_parallel[0]), float(history.a_parallel[-1])],
            "r": [float(history.radius[0]), float(history.radius[-1])],
            "lapse": [float(history.lapse[0]), float(history.lapse[-1])],
            "shift": [float(history.shift[0]), float(history.shift[-1])],
        },
    }


def _payload_manifest(path, arrays):
    raw = path.read_bytes()
    return {
        "path": str(path.relative_to(ROOT)), "bytes": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "format": "deterministic uncompressed NPZ; diagnostic control histories",
        "arrays": {
            name: {"dtype": value.dtype.str, "shape": list(value.shape),
                   "sha256": array_digest(name, value)}
            for name, value in sorted(arrays.items())
        },
    }


def calculate(path, arrays, records):
    state_manifest, selected, coupled, adm = records
    state = _load_seed(state_manifest)
    geometry = selected["source_selected_initial_geometry"]
    frame = selected["source_selected_frame"]
    short = endpoint_control_history(
        duration=0.5, points=257, seed_radius=geometry["seed_r"],
        selected_radius=geometry["r"], a0=geometry["a_parallel"],
        H_parallel=geometry["H_parallel"],
        landau_velocity=frame["velocity_old_child_frame"],
        landau_gamma=frame["gamma"], label="endpoint_control_duration_0.5",
    )
    long = endpoint_control_history(
        duration=1.0, points=513, seed_radius=geometry["seed_r"],
        selected_radius=geometry["r"], a0=geometry["a_parallel"],
        H_parallel=geometry["H_parallel"],
        landau_velocity=frame["velocity_old_child_frame"],
        landau_gamma=frame["gamma"], label="endpoint_control_duration_1.0",
    )
    # Recomputed here and required to match the authenticated artifact.
    short_result = propagate_blockwise(
        state.arrays, short, seed_a=geometry["a_parallel"], seed_r=geometry["seed_r"]
    )
    long_result = propagate_blockwise(
        state.arrays, long, seed_a=geometry["a_parallel"], seed_r=geometry["seed_r"]
    )
    calculated = {
        "U_short": short_result["unitary"], "C_short": short_result["covariance"],
        "U_long": long_result["unitary"], "C_long": long_result["covariance"],
        "short_time": short.time, "short_a": short.a_parallel,
        "short_r": short.radius, "short_lapse": short.lapse, "short_shift": short.shift,
        "long_time": long.time, "long_a": long.a_parallel,
        "long_r": long.radius, "long_lapse": long.lapse, "long_shift": long.shift,
    }
    for name in calculated:
        if name not in arrays or arrays[name].shape != calculated[name].shape:
            raise ValueError(f"control artifact shape differs: {name}")
        if np.max(abs(arrays[name]-calculated[name])) > 3e-13:
            raise ValueError(f"control artifact value differs: {name}")
    difference = float(np.max(abs(short_result["unitary"]-long_result["unitary"])))
    if difference <= 1e-4:
        raise ArithmeticError("endpoint-identical histories did not expose path dependence")
    tolerance = 3e-11
    conditional_pass = (
        max(short_result["maximum_unitarity_residual"], long_result["maximum_unitarity_residual"]) < tolerance
        and min(short_result["minimum_covariance_eigenvalue"], long_result["minimum_covariance_eigenvalue"]) >= -tolerance
        and max(short_result["maximum_covariance_eigenvalue"], long_result["maximum_covariance_eigenvalue"]) <= 1+tolerance
    )
    if not conditional_pass:
        raise ArithmeticError("conditional history evolution lost unitarity or CAR")
    locked = selected["locked_inputs"]
    return {
        "schema": "NSC-LANDAU-CAUCHY-ISOMETRY-v1",
        "status": (
            "CONDITIONAL API PASS / PHYSICAL SELECTION FAIL: each declared KS "
            "history gives a unitary CAR-preserving map, but endpoint-identical "
            "histories give different U and no same-action history selects one"
        ),
        "source_hashes": hashes(SOURCES), "input_hashes": hashes(INPUTS),
        "locked_inputs": {
            "magnetic_flux": locked["q"], "Omega": locked["Omega"],
            "zeta": locked["zeta"], "A": locked["A"], "V_full": locked["V_full"],
            "occupations_refitted": False,
        },
        "conditional_isometry_owner": {
            "formula": "U[history]=T exp[-i integral H_ADM[history] dt] in canonical half-density blocks",
            "Hamiltonian": "N[-m sigma1+(lambda/r)sigma2+(k/a)sigma3]-shift*k*I",
            "frequency_mixing_interface_allowed": True,
            "controls_are_frequency_diagonal_because_they_are_homogeneous": True,
            "full_dense_embedding_used": False,
        },
        "control_histories": [_history_dict(short), _history_dict(long)],
        "control_results": {
            "short": {key: value for key, value in short_result.items() if key not in {"unitary", "covariance"}},
            "long": {key: value for key, value in long_result.items() if key not in {"unitary", "covariance"}},
            "maximum_U_short_minus_U_long": difference,
            "same_endpoint_data": True,
        },
        "payload": _payload_manifest(path, arrays),
        "physical_history_gate": {
            "selected": False,
            "reason": (
                "the constraint record fixes endpoint data but no lapse/shift/metric "
                "history or Landau hypersurface embedding between non-isometric slices"
            ),
            "circular_dependency": (
                "the requested history would be produced by coupled metric evolution, "
                "while that evolution was forbidden until this map exists"
            ),
            "local_B_sqrt_used_as_U": False,
            "conditional_control_U_claimed_physical": False,
        },
        "finite_selected_stress": {
            "evaluated": False,
            "reason": (
                "different conditional U already give different target covariance; "
                "the general-KS fourth-order reference and local induced history are also absent"
            ),
            "old_tensor_copied": False,
            "linearized_density_used": False,
            "null_signs": None,
        },
        "constraints": {
            "fixed_tensor_source_selected_residuals": selected["constraints"]["residuals"],
            "updated_with_finite_target_stress": False,
            "ADM_owner_rederived": False,
        },
        "gate": {
            "conditional_U_unitarity_pass": conditional_pass,
            "conditional_CAR_pass": conditional_pass,
            "endpoint_data_select_unique_U": False,
            "physical_U_L0_pass": False,
            "finite_rstar_stress_pass": False,
            "constraints_with_updated_stress_pass": False,
            "full_contract_pass": False,
            "coupled_evolution_reopened": False,
            "next_owner": "joint history/state boundary-value solve or an explicit same-action history-selection equation",
        },
        "scope": {
            "new_physical_metric_ansatz_selected": False,
            "metric_time_step_started": False,
            "A_q_Omega_zeta_or_V_full_changed": False,
            "pointwise_spin_boost_used_as_isometry": False,
            "four_moment_surrogate_or_dark_fit": False,
        },
        "comparison": {
            "fields": "all JSON fields and control artifact arrays",
            "exact": "hashes, schema, labels, dtypes and shapes",
            "float_atol": 3e-12, "float_rtol": 3e-12, "exceptions": [],
        },
    }


def _load_artifact(record):
    path = ROOT/record["payload"]["path"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != record["payload"]["sha256"]:
        raise ValueError("isometry control payload hash mismatch")
    with np.load(path, allow_pickle=False) as data:
        arrays = {name: np.array(data[name], copy=True) for name in data.files}
    return path, arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args();records = _load_inputs()
    if args.check:
        expected = json.loads(OUTPUT.read_text());path, arrays = _load_artifact(expected)
        actual = calculate(path, arrays, records);compare(expected, actual)
        print("conditional Cauchy isometries verified; physical history remains unselected")
        return
    output = args.output or OUTPUT
    if output.exists(): raise FileExistsError("refusing to overwrite recorded evidence")
    state = _load_seed(records[0]);selected=records[1]
    geometry=selected["source_selected_initial_geometry"];frame=selected["source_selected_frame"]
    histories=[
        endpoint_control_history(duration=d,points=p,seed_radius=geometry["seed_r"],selected_radius=geometry["r"],a0=geometry["a_parallel"],H_parallel=geometry["H_parallel"],landau_velocity=frame["velocity_old_child_frame"],landau_gamma=frame["gamma"],label=f"endpoint_control_duration_{d}")
        for d,p in ((0.5,257),(1.0,513))
    ]
    results=[propagate_blockwise(state.arrays,h,seed_a=geometry["a_parallel"],seed_r=geometry["seed_r"]) for h in histories]
    arrays={"U_short":results[0]["unitary"],"C_short":results[0]["covariance"],"U_long":results[1]["unitary"],"C_long":results[1]["covariance"]}
    for prefix,h in zip(("short","long"),histories):
        arrays.update({f"{prefix}_time":h.time,f"{prefix}_a":h.a_parallel,f"{prefix}_r":h.radius,f"{prefix}_lapse":h.lapse,f"{prefix}_shift":h.shift})
    raw=deterministic_npz_bytes(arrays);digest=hashlib.sha256(raw).hexdigest();ARTIFACT_DIRECTORY.mkdir(parents=True,exist_ok=True)
    path=ARTIFACT_DIRECTORY/f"nsc-landau-cauchy-control-isometries.{digest}.npz";path.write_bytes(raw)
    record=calculate(path,arrays,records);output.write_text(json.dumps(record,indent=2,sort_keys=True)+"\n")
    print(output);print(path)


if __name__ == "__main__": main()
