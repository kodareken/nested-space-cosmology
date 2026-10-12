#!/usr/bin/env python3
"""Discovery atlas of sealed frames. No new trajectory and no initial solve.

``--check`` hashes the sealed preparations and recomputes the nf256 T=0.3
reference-window numbers. ``--write`` evaluates the saved baseline frames
and creates one successor JSON+NPZ. It refuses to replace either file.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from dataclasses import replace
from pathlib import Path
import sys

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np

from recursive_horizons import nsc_nested_parent_child as nested
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons.nsc_discovery_observables import (
    API_VERSION,
    frozen_geometry_policy,
    observe,
    reference_window_integral,
)
from recursive_horizons.nsc_spherical_coupling import OCCUPATIONS, CauchyState


LAB = Path(__file__).resolve().parents[1]
OUT = LAB / "results/development/nsc-discovery-atlas-v1.json"
SCHEMA = "NSC-DISCOVERY-ATLAS-v1"
CHECK_SCHEMA = "NSC-DISCOVERY-ATLAS-CHECK-v1"
CAUCHY_NAMES = ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1")
BASIS_JSON = "results/development/nsc-nested-parent-child-replay-basis-v1.json"
BASIS_NPZ = "results/development/nsc-nested-parent-child-replay-basis-v1.npz"
BASELINE_CASES = (
    {"preparation": "original_conformal_v2", "case": "nf256_dt_0.0005", "nf": 256, "store": "conformal",
     "npz": "results/development/nsc-spherical-conformal-episode-v2.npz"},
    {"preparation": "original_conformal_v2", "case": "nf512_dt_0.0005", "nf": 512, "store": "conformal",
     "npz": "results/development/nsc-spherical-conformal-episode-v2.npz"},
    {"preparation": "separated_pair_v1", "case": "nf256_baseline_dt0.0005", "nf": 256, "store": "nested",
     "npz": "results/development/nsc-nested-parent-child-v1.npz"},
    {"preparation": "separated_pair_v2", "case": "nf512_baseline_dt0.0005", "nf": 512, "store": "nested",
     "npz": "results/development/nsc-nested-parent-child-confirmation-v2.npz"},
)
FINGERPRINT_KEYS = ("schema", "gauge", "source_layout", "status", "parent_interval", "child_interval")
DEPENDENCIES = {
    "producer": ("scripts/derive_nsc_discovery_atlas.py",),
    "observer": ("src/recursive_horizons/nsc_discovery_observables.py",),
    "backend": (
        "src/recursive_horizons/nsc_spherical_galerkin_coupling.py",
        "src/recursive_horizons/nsc_regional_energy_exchange.py",
        "src/recursive_horizons/nsc_spherical_episode_assessment.py",
        "scripts/derive_nsc_spherical_conformal_curvature.py",
    ),
    "core": (
        "src/recursive_horizons/nsc_nested_parent_child.py",
        "src/recursive_horizons/nsc_spherical_coupling.py",
        "src/recursive_horizons/nsc_spherical_feedback_action.py",
        "src/recursive_horizons/nsc_conformal_adm_source.py",
    ),
    "consumer": ("scripts/derive_nsc_spherical_conformal_frames.py",),
}
PREPARATIONS = (
    {
        "name": "original_conformal_v2",
        "schema": "NSC-SPHERICAL-CONFORMAL-EPISODE-v2",
        "json": "results/development/nsc-spherical-conformal-episode-v2.json",
        "npz": "results/development/nsc-spherical-conformal-episode-v2.npz",
    },
    {
        "name": "separated_pair_v1",
        "schema": "NSC-NESTED-PARENT-CHILD-v1",
        "json": "results/development/nsc-nested-parent-child-v1.json",
        "npz": "results/development/nsc-nested-parent-child-v1.npz",
    },
    {
        "name": "separated_pair_v2",
        "schema": "NSC-NESTED-PARENT-CHILD-CONFIRMATION-v2",
        "json": "results/development/nsc-nested-parent-child-confirmation-v2.json",
        "npz": "results/development/nsc-nested-parent-child-confirmation-v2.npz",
    },
    {
        "name": "frozen_basis",
        "schema": "NSC-NESTED-PARENT-CHILD-REPLAY-BASIS-v1",
        "json": "results/development/nsc-nested-parent-child-replay-basis-v1.json",
        "npz": "results/development/nsc-nested-parent-child-replay-basis-v1.npz",
    },
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _header(path):
    record = json.loads(Path(path).read_text())
    header = {key: record.get(key) for key in ("schema", "gauge", "source_layout", "status", "parent_interval", "child_interval", "payload_sha256")}
    saved = {}
    validation = record.get("validation")
    if isinstance(validation, dict):
        saved["validation"] = {
            "SVD_reconstructed_for_replay": validation.get("SVD_reconstructed_for_replay"),
            "maximum_metric_gap": validation.get("maximum_metric_gap"),
        }
    bindings = record.get("source_bindings")
    if isinstance(bindings, dict):
        saved["source_bindings"] = bindings
    return header, saved


def _fingerprint(header):
    return tuple(json.dumps(header.get(key), sort_keys=True, allow_nan=False) for key in FINGERPRINT_KEYS)


def _geometry_comparison(spec, json_hash, basis_saved):
    policy = frozen_geometry_policy(spec["name"])
    comparison = {
        "status": policy["status"],
        "control": policy["control"],
        "reason": policy["reason"],
        "saved_control_applied": policy["status"] == "saved_control",
        "svd_reconstructed": False,
        "saved_validation_svd_reconstructed": None,
        "saved_maximum_metric_gap": None,
        "hash_matches_saved_control": None,
    }
    if policy["status"] != "saved_control":
        return comparison
    validation = basis_saved.get("validation", {})
    comparison["saved_validation_svd_reconstructed"] = validation.get("SVD_reconstructed_for_replay")
    comparison["saved_maximum_metric_gap"] = validation.get("maximum_metric_gap")
    if spec["name"] == "frozen_basis":
        comparison["hash_matches_saved_control"] = True
        return comparison
    expected = basis_saved.get("source_bindings", {}).get(spec["json"])
    comparison["hash_matches_saved_control"] = expected is not None and expected == json_hash
    return comparison


def check_inputs(lab=LAB):
    """Hash and classify sealed inputs. Does not open trajectory arrays."""
    lab = Path(lab)
    loaded = {}
    rows = []
    for spec in PREPARATIONS:
        json_path = lab / spec["json"]
        npz_path = lab / spec["npz"]
        if not json_path.is_file() or not npz_path.is_file():
            raise FileNotFoundError(f"missing sealed preparation {spec['name']}")
        header, saved = _header(json_path)
        if header["schema"] != spec["schema"]:
            raise ValueError(f"{spec['name']} schema is {header['schema']}, expected {spec['schema']}")
        payload_hash = sha256(npz_path)
        json_hash = sha256(json_path)
        if header.get("payload_sha256") != payload_hash:
            raise ValueError(f"{spec['name']} payload hash does not match its sealed npz")
        loaded[spec["name"]] = {"header": header, "saved": saved, "json_sha256": json_hash, "npz_sha256": payload_hash}
    basis_saved = loaded["frozen_basis"]["saved"]
    if "source_bindings" not in basis_saved or "validation" not in basis_saved:
        raise ValueError("frozen replay basis has no saved control metadata")
    for spec in PREPARATIONS:
        item = loaded[spec["name"]]
        header = item["header"]
        rows.append({
            "name": spec["name"],
            "schema": header["schema"],
            "gauge": header.get("gauge"),
            "source_layout": header.get("source_layout"),
            "status": header.get("status"),
            "parent_interval": header.get("parent_interval"),
            "child_interval": header.get("child_interval"),
            "json": spec["json"],
            "npz": spec["npz"],
            "json_sha256": item["json_sha256"],
            "npz_sha256": item["npz_sha256"],
            "payload_sha256_matches": True,
            "fingerprint": list(_fingerprint(header)),
            "frozen_geometry_comparison": _geometry_comparison(spec, item["json_sha256"], basis_saved),
            "trajectory_arrays_copied": False,
        })
    new_pair = frozen_geometry_policy("new_pair")
    rows.append({
        "name": "new_pair",
        "schema": None,
        "input_record": False,
        "frozen_geometry_comparison": {
            "status": new_pair["status"],
            "control": new_pair["control"],
            "reason": new_pair["reason"],
            "saved_control_applied": False,
            "svd_reconstructed": False,
        },
        "trajectory_arrays_copied": False,
    })
    fingerprints = [tuple(row["fingerprint"]) for row in rows if "fingerprint" in row]
    if len(set(fingerprints)) != len(PREPARATIONS):
        raise ValueError("preparation metadata fingerprints are not distinct")
    if len({row["schema"] for row in rows if row.get("schema")}) != len(PREPARATIONS):
        raise ValueError("preparation schemas are not distinct")
    for row in rows:
        comparison = row["frozen_geometry_comparison"]
        if comparison["svd_reconstructed"] is not False:
            raise ValueError("atlas check must not reconstruct a geometry frame")
        if row["name"] in ("separated_pair_v1", "separated_pair_v2", "frozen_basis") and comparison["hash_matches_saved_control"] is not True:
            raise ValueError(row["name"] + " does not match the saved frozen-geometry control")
        if row["name"] in ("original_conformal_v2", "new_pair") and comparison["status"] != "missing":
            raise ValueError(row["name"] + " must keep a missing frozen-geometry control")
    return {
        "schema": CHECK_SCHEMA,
        "api_version": API_VERSION,
        "writer_executed": False,
        "evolution_performed": False,
        "trajectories_loaded": False,
        "npz_arrays_opened": False,
        "observe_called": False,
        "metadata_distinct": True,
        "output": OUT.relative_to(lab).as_posix() if lab == LAB else str(OUT),
        "output_exists": OUT.exists(),
        "preparations": rows,
    }


def atlas_record(report):
    preparations = []
    for row in report["preparations"]:
        item = {
            "name": row["name"],
            "schema": row.get("schema"),
            "gauge": row.get("gauge"),
            "source_layout": row.get("source_layout"),
            "status": row.get("status"),
            "parent_interval": row.get("parent_interval"),
            "child_interval": row.get("child_interval"),
            "frozen_geometry_comparison": row["frozen_geometry_comparison"],
            "trajectory_arrays_copied": False,
        }
        if "json_sha256" in row:
            item["json"] = row["json"]
            item["json_sha256"] = row["json_sha256"]
            item["npz_sha256"] = row["npz_sha256"]
        preparations.append(item)
    return {
        "schema": SCHEMA,
        "api_version": API_VERSION,
        "creation_only": True,
        "evolution_performed": False,
        "trajectories_loaded": False,
        "observe_called": False,
        "metadata_distinct": True,
        "successor_of": [spec["json"] for spec in PREPARATIONS],
        "preparations": preparations,
    }


def _readonly(value, dtype=None):
    out = np.array(value, dtype=dtype, copy=True)
    out.setflags(write=False)
    return out


def load_frozen_pair(nf, lab=LAB, weights=None):
    """Saved replay-basis frame. The geometry map is not rebuilt."""
    entry = json.loads((Path(lab) / BASIS_JSON).read_text())["bases"][str(int(nf))]
    weights = np.array(OCCUPATIONS if weights is None else weights, dtype=float)
    with np.load(Path(lab) / BASIS_NPZ) as stored:
        matrix = stored[f"nf{int(nf)}_W"]
        reference = stored[f"nf{int(nf)}_reference_columns"]
        source = stored[f"nf{int(nf)}_source_columns"]
        grid = galerkin.build_grid(int(nf), gauge="conformal")
        grid.fine = replace(grid.fine, occupations=weights.copy())
        return nested.NestedPair(
            grid, _readonly(matrix), _readonly(entry["coarse_indices"], int),
            _readonly(entry["child_indices"], int), _readonly(entry["parent_indices"], int),
            _readonly(reference[:int(nf)]), _readonly(reference[int(nf):]), _readonly(reference),
            _readonly(source[:int(nf)]), _readonly(source[int(nf):]), _readonly(source),
            _readonly(weights), entry["source"], entry["geometry"],
        )


def _frame_state(spec, stored, pair, index):
    if spec["store"] == "conformal":
        cauchy = CauchyState(**{name: np.array(stored[f"{spec['case']}_frames_{name}"][index], copy=True) for name in CAUCHY_NAMES})
        state = nested.encode_state(pair, cauchy)
        gap = float(np.max(np.abs(state.p_Q - pair.grid.dx_g * pair.geometry_map.T @ cauchy.p_Q)))
        return state, gap
    state = nested.NestedState(*(np.array(stored[f"{spec['case']}_{name}"][index], copy=True) for name in nested.STATE_NAMES))
    saved_reference = stored[f"{spec['case']}_reference_columns"]
    if not np.array_equal(saved_reference, pair.reference_columns):
        raise ValueError("frozen basis observer differs from the saved pair reference")
    return state, 0.0


def evaluate_case(spec, lab=LAB, indices=None):
    """Observe saved frames. The Hamiltonian is evaluated; nothing is integrated forward."""
    pair = load_frozen_pair(spec["nf"], lab)
    with np.load(Path(lab) / spec["npz"]) as stored:
        time_key = spec["case"] + ("_frame_times" if spec["store"] == "conformal" else "_times")
        times = np.array(stored[time_key], dtype=float, copy=True)
        chosen = list(range(len(times))) if indices is None else [int(index) % len(times) for index in indices]
        frames = []
        for index in chosen:
            state, momentum_gap = _frame_state(spec, stored, pair, index)
            nodal = nested.reconstruct_state(pair, state)
            rate, bundle = galerkin.compose_fine_hamiltonian(pair.grid, nodal)
            coupled = observe(
                pair, state, float(times[index]), bundle=bundle, nodal_rate=rate,
                preparation=spec["preparation"], control_mode="coupled",
                capture_profiles=(index == chosen[-1]),
            )
            item = {"time": float(times[index]), "row": coupled, "momentum_conversion_gap": momentum_gap}
            if index == chosen[-1]:
                item["frozen"] = observe(
                    pair, state, float(times[index]), bundle=bundle, nodal_rate=rate,
                    preparation=spec["preparation"], control_mode="frozen_geometry",
                )
            frames.append(item)
    return {"spec": spec, "frames": frames}


def _summands_for_saved_frame(spec, lab=LAB):
    """One saved frame's regional summands, using the same evaluation as the atlas row."""
    import derive_nsc_spherical_conformal_frames as frames
    from recursive_horizons import nsc_regional_energy_exchange as regional
    pair = load_frozen_pair(spec["nf"], lab)
    with np.load(Path(lab) / spec["npz"]) as stored:
        time_key = spec["case"] + ("_frame_times" if spec["store"] == "conformal" else "_times")
        times = np.array(stored[time_key], dtype=float, copy=True)
        state, momentum_gap = _frame_state(spec, stored, pair, -1)
    nodal = nested.reconstruct_state(pair, state)
    rate, bundle = galerkin.compose_fine_hamiltonian(pair.grid, nodal)
    coupled = observe(
        pair, state, float(times[-1]), bundle=bundle, nodal_rate=rate,
        preparation=spec["preparation"], control_mode="coupled", capture_profiles=True,
    )
    frozen = observe(
        pair, state, float(times[-1]), bundle=bundle, nodal_rate=rate,
        preparation=spec["preparation"], control_mode="frozen_geometry",
    )
    fine = bundle["fine_state"]
    lifted_q = bundle["lifted_Q"]
    coordinate_work = bundle["source"]["force_Q"] * lifted_q + bundle["source"]["force_L"] * bundle["lapse_dot"]
    ledger = regional.matter_ledger(bundle["fine_system"], fine)
    samples = {
        "normal_energy": np.asarray(ledger["normal_energy_nodal"], dtype=float),
        "coordinate_metric_work": np.asarray(coordinate_work, dtype=float),
        "flux_shift": np.asarray(ledger["flux_shift"], dtype=float),
        "flux_proper": np.asarray(ledger["flux_proper"], dtype=float),
        "flux_shift_pressure": np.asarray(ledger["flux_shift_pressure"], dtype=float),
    }
    reported = {
        "normal_energy": coupled["normal_energy_total"],
        "coordinate_metric_work": coupled["coordinate_metric_work_total"],
        "flux_shift": coupled["flux_projection_shift_total"],
        "flux_proper": coupled["flux_projection_proper_total"],
        "flux_shift_pressure": coupled["flux_projection_shift_pressure_total"],
    }
    channels = {}
    for name, values in samples.items():
        raw = float(np.sum(values))
        window = float(np.dot(frames.interval_weights(pair.grid.nq, pair.grid.length, 0.0, pair.grid.length), values))
        again = reference_window_integral(pair.grid, values, (0.0, float(pair.grid.length)))
        if max(abs(reported[name] - raw), abs(reported[name] - window), abs(reported[name] - again)) > 1e-8 * max(1.0, abs(raw)):
            raise ValueError(name + " does not match the reference sum on the saved frame")
        channels[name] = {"raw_sum": raw, "reference_window": window, "reported": reported[name]}
    return {
        "time": float(times[-1]),
        "count": int(pair.grid.nq),
        "channels": channels,
        "normal_energy_summand": samples["normal_energy"],
        "momentum_conversion_gap": momentum_gap,
        "frozen_Q_dot_max": frozen["curvature"]["Q_dot_max"],
        "frozen_coordinate_metric_work_total": frozen["coordinate_metric_work_total"],
        "frozen_child_lapse_work": frozen["child_window"]["momentum_lapse_work"],
        "coupled_child_lapse_work": coupled["child_window"]["momentum_lapse_work"],
        "coupled_child_r_mean": coupled["proper"]["child_r_proper_mean"],
        "frozen_child_r_mean": frozen["proper"]["child_r_proper_mean"],
        "coupled_Q_dot_max": coupled["curvature"]["Q_dot_max"],
    }


def saved_frame_reference(lab=LAB):
    spec = next(row for row in BASELINE_CASES if row["case"] == "nf256_dt_0.0005")
    return _summands_for_saved_frame(spec, lab)


def baseline_array_shapes(lab=LAB):
    shapes = {}
    for spec in BASELINE_CASES:
        with np.load(Path(lab) / spec["npz"]) as stored:
            time_key = spec["case"] + ("_frame_times" if spec["store"] == "conformal" else "_times")
            phi_key = spec["case"] + ("_frames_phi0" if spec["store"] == "conformal" else "_phi0")
            times = stored[time_key]
            phi = stored[phi_key]
            shapes[spec["preparation"] + ":" + spec["case"]] = {
                "frames": int(times.shape[0]), "final_time": float(times[-1]), "phi": list(phi.shape),
            }
    return shapes


def bound_input_paths():
    return tuple(spec[kind] for spec in PREPARATIONS for kind in ("json", "npz"))


def _repository():
    return Path(LAB).parent


def evidence_binding(lab=LAB):
    """Hashes of the producing code and the sealed inputs. Nothing is written."""
    lab = Path(lab)
    dependencies = {}
    watched = []
    for role, relatives in DEPENDENCIES.items():
        dependencies[role] = {}
        for relative in relatives:
            path = lab / relative
            if not path.is_file():
                raise FileNotFoundError(relative)
            dependencies[role][relative] = sha256(path)
            watched.append("lab/" + relative)
    inputs = {}
    for relative in bound_input_paths():
        path = lab / relative
        if not path.is_file():
            raise FileNotFoundError(relative)
        inputs[relative] = sha256(path)
        watched.append("lab/" + relative)
    commit = subprocess.check_output(
        ["git", "-C", str(_repository()), "rev-parse", "HEAD"], text=True,
    ).strip()
    status = subprocess.check_output(
        ["git", "-C", str(_repository()), "status", "--porcelain", "--untracked-files=all", "--", *watched],
        text=True,
    )
    dirty = {}
    for line in status.splitlines():
        if len(line) < 4:
            continue
        name = line[3:]
        if " -> " in name:
            name = name.split(" -> ", 1)[1]
        dirty[name] = line[:2]
    return {
        "producing_commit": commit,
        "dirty_files": dirty,
        "dependencies": dependencies,
        "inputs": inputs,
        "sealed_inputs_regenerated": False,
    }


def _series_row(row):
    curvature = row["curvature"]
    localization = row["localization"]
    return {
        "time": row["time"],
        "R_h_min": curvature["R_h_min"],
        "R_h_max": curvature["R_h_max"],
        "chi_shell_gap_max": curvature["chi_shell_gap_max"],
        "Q_dot_max": curvature["Q_dot_max"],
        "Q_ddot_max": curvature["Q_ddot_max"],
        "localization_peak_x": localization["peak_x"],
        "localization_proper_width": localization["proper_width"],
        "localization_resultant": localization["resultant"],
        "normal_energy_total": row["normal_energy_total"],
        "coordinate_metric_work_total": row["coordinate_metric_work_total"],
        "flux_projection_proper_total": row["flux_projection_proper_total"],
        "child_normal_energy": row["child_window"]["normal_energy"],
        "parent_normal_energy": row["parent_window"]["normal_energy"],
        "child_boundary_flux": row["child_window"]["normal_boundary_flux"],
        "child_pressure_work": row["child_window"]["proper_pressure_work"],
        "child_lapse_work": row["child_window"]["momentum_lapse_work"],
        "child_metric_work": row["child_window"]["coordinate_metric_work"],
    }


def write_atlas(path=OUT, lab=LAB, *, cases=None, indices=None):
    """Create one successor JSON+NPZ of saved-frame observations. No evolution."""
    destination = Path(path)
    payload = destination.with_suffix(".npz")
    sealed = {(Path(lab) / spec[kind]).resolve() for spec in PREPARATIONS for kind in ("json", "npz")}
    sealed.add((Path(lab) / BASIS_NPZ).resolve())
    if destination.resolve() in sealed or payload.resolve() in sealed:
        raise FileExistsError("refusing to overwrite a sealed preparation")
    if destination.exists() or payload.exists():
        raise FileExistsError("refusing to overwrite a discovery atlas successor")
    selected = BASELINE_CASES if cases is None else tuple(cases)
    series, arrays = [], {}
    for spec in selected:
        evaluated = evaluate_case(spec, lab, indices=indices)
        rows = [_series_row(frame["row"]) for frame in evaluated["frames"]]
        final = evaluated["frames"][-1]
        frozen = final["frozen"]
        prefix = spec["preparation"] + "_" + spec["case"]
        arrays[prefix + "_times"] = np.array([frame["time"] for frame in evaluated["frames"]])
        arrays[prefix + "_normal_energy_total"] = np.array([row["normal_energy_total"] for row in rows])
        arrays[prefix + "_final_R_h"] = np.array(final["row"]["profiles"]["R_h"], copy=True)
        arrays[prefix + "_final_eta_N"] = np.array(final["row"]["profiles"]["normal_energy_density"], copy=True)
        arrays[prefix + "_final_Q_dot"] = np.array(final["row"]["profiles"]["Q_dot"], copy=True)
        series.append({
            "preparation": spec["preparation"], "case": spec["case"], "nf": spec["nf"], "store": spec["store"],
            "control_mode": "coupled", "rows": rows, "momentum_conversion_gap": final["momentum_conversion_gap"],
            "final_frozen_geometry": {
                "Q_dot_max": frozen["curvature"]["Q_dot_max"],
                "Q_ddot_max": frozen["curvature"]["Q_ddot_max"],
                "coordinate_metric_work_total": frozen["coordinate_metric_work_total"],
                "child_momentum_lapse_work": frozen["child_window"]["momentum_lapse_work"],
                "child_proper_pressure_work": frozen["child_window"]["proper_pressure_work"],
                "child_r_proper_mean": frozen["proper"]["child_r_proper_mean"],
            },
        })
    record = {
        "schema": SCHEMA, "api_version": API_VERSION, "creation_only": True,
        "evolution_performed": False, "initial_solve_performed": False, "trajectories_created": False,
        "evidence": evidence_binding(lab),
        "series": series,
    }
    with payload.open("xb") as handle:
        np.savez_compressed(handle, **arrays)
    record["payload_sha256"] = sha256(payload)
    with destination.open("x", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2, allow_nan=False)
        handle.write("\n")
    return record


def authenticate_successor(path, lab=LAB):
    """Read a successor and its payload. Compare hashes and stored rows. Do not write."""
    destination = Path(path)
    payload = destination.with_suffix(".npz")
    before = (destination.stat().st_mtime_ns, payload.stat().st_mtime_ns, destination.stat().st_size, payload.stat().st_size)
    record = json.loads(destination.read_text())
    if record.get("schema") != SCHEMA or record.get("evolution_performed") is not False:
        raise ValueError("successor is not an unevolved discovery atlas")
    if sha256(payload) != record.get("payload_sha256"):
        raise ValueError("successor payload hash does not match")
    evidence = record.get("evidence") or {}
    dependencies = evidence.get("dependencies") or {}
    if tuple(dependencies) != tuple(DEPENDENCIES):
        raise ValueError("successor dependency roles differ")
    for role, files in DEPENDENCIES.items():
        bound = dependencies[role]
        if tuple(bound) != files:
            raise ValueError("successor " + role + " dependency list differs")
        for relative, expected in bound.items():
            if sha256(Path(lab) / relative) != expected:
                raise ValueError("dependency bytes changed: " + relative)
    inputs = evidence.get("inputs") or {}
    if tuple(inputs) != bound_input_paths():
        raise ValueError("successor input list differs")
    for relative, expected in inputs.items():
        if sha256(Path(lab) / relative) != expected:
            raise ValueError("sealed input changed: " + relative)
    if not evidence.get("producing_commit") or evidence.get("sealed_inputs_regenerated") is not False:
        raise ValueError("successor evidence is missing its commit or regenerated an input")
    with np.load(payload) as stored:
        for series in record["series"]:
            prefix = series["preparation"] + "_" + series["case"]
            times = stored[prefix + "_times"]
            energy = stored[prefix + "_normal_energy_total"]
            curvature = stored[prefix + "_final_R_h"]
            if len(series["rows"]) != len(times):
                raise ValueError("stored row count differs from the payload")
            for index, row in enumerate(series["rows"]):
                if abs(float(times[index]) - float(row["time"])) > 1e-12:
                    raise ValueError("stored time differs from the payload")
                if abs(float(energy[index]) - float(row["normal_energy_total"])) > 1e-9:
                    raise ValueError("stored normal energy differs from the payload")
                if not np.isfinite(row["R_h_max"]) or not np.isfinite(row["localization_peak_x"]):
                    raise ValueError("stored curvature or localisation is not finite")
            if abs(float(np.max(curvature)) - float(series["rows"][-1]["R_h_max"])) > 1e-8:
                raise ValueError("final curvature profile does not match the last row")
    after = (destination.stat().st_mtime_ns, payload.stat().st_mtime_ns, destination.stat().st_size, payload.stat().st_size)
    if before != after:
        raise ValueError("authentication rewrote the successor")
    return {
        "authenticated": True,
        "series": len(record["series"]),
        "producing_commit": evidence["producing_commit"],
        "dirty_files": evidence.get("dirty_files", {}),
        "rewritten": False,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="validate sealed inputs and saved-frame numbers; write nothing")
    parser.add_argument("--write", action="store_true", help="create the successor JSON+NPZ; refuse to replace it")
    parser.add_argument("--output", type=Path, default=None, help="successor path used only with --write")
    parser.add_argument("--successor", type=Path, default=None, help="existing successor to authenticate with --check")
    args = parser.parse_args(argv)
    if args.check == args.write:
        parser.error("choose exactly one of --check or --write")
    if args.check:
        if args.output is not None:
            parser.error("--output is only valid with --write")
        report = check_inputs()
        report["baseline_shapes"] = baseline_array_shapes()
        report["saved_frame"] = {
            key: value for key, value in saved_frame_reference().items() if key != "normal_energy_summand"
        }
        successor = args.successor if args.successor is not None else (OUT if OUT.exists() else None)
        if successor is not None:
            report["successor"] = authenticate_successor(successor)
        report["writer_executed"] = False
        json.dump(report, sys.stdout, indent=2, allow_nan=False)
        sys.stdout.write("\n")
        return 0
    destination = OUT if args.output is None else args.output
    write_atlas(destination)
    print("created", destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
