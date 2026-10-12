#!/usr/bin/env python3
"""Read immutable parent initial slices: constraints, momenta, and projected metric jets.

Default JSON goes to stdout. --output PATH creates one immutable summary once.
Raw same-point Hamiltonian jets are diagnostic references, never acceptance or
a replacement of the actual projected metric. No preparation or evolution runs.
"""
from __future__ import annotations

import argparse
import ast
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from types import SimpleNamespace

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

import numpy as np
from threadpoolctl import threadpool_limits
from recursive_horizons import nsc_discovery_backend as backend
from recursive_horizons import nsc_discovery_leading_einstein as leading
from recursive_horizons import nsc_discovery_tidal as tidal
from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons import provenance

SCHEMA = "NSC-DISCOVERY-PARENT-INITIAL-ASSESSMENT-v1"
DEFAULT = leading.LAB/"results/development/nsc-discovery-parent-v1/feasibility-nf64-128"


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_hash(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def authenticate(directory):
    """Authenticate stored data and historical preparation without calling its current constructor."""
    directory = Path(directory).expanduser().resolve()
    paths = (directory/"parent.json", directory/"parent.npz")
    if any(path.stat().st_mode & 0o222 for path in paths):
        raise ValueError("parent input must be immutable JSON/NPZ")
    record = json.loads(paths[0].read_text())
    if record["schema"] != "NSC-DISCOVERY-PARENT-v1" or file_hash(paths[1]) != record["payload_sha256"]:
        raise ValueError("parent record/payload identity mismatch")
    with np.load(paths[1], allow_pickle=False) as saved:
        arrays = {key: np.array(saved[key], copy=True) for key in saved.files}
    if set(arrays) != set(record["array_sha256"]):
        raise ValueError("parent array inventory mismatch")
    if any(array_hash(arrays[key]) != digest for key, digest in record["array_sha256"].items()):
        raise ValueError("stored parent array changed")
    for path, digest in record["input_hashes"].items():
        if file_hash(leading.ROOT/path) != digest:
            raise ValueError("historical parent input changed: "+path)
    historical_parent = None
    for path, digest in record["producers"].items():
        raw = provenance.resolve_pinned_source_bytes(leading.ROOT, path, digest, commit=record["producing_commit"])
        if path == "lab/src/recursive_horizons/nsc_discovery_parent.py":
            historical_parent = raw
    profile, profile_binding = declared_profile(directory, record, historical_parent)
    binding = {"directory": str(directory), "record_sha256": file_hash(paths[0]),
               "payload_sha256": record["payload_sha256"], "producing_commit": record["producing_commit"],
               "historical_producers_authenticated": True, "current_preparation_required": False,
               "declared_profile": profile, "profile_binding": profile_binding}
    return record, arrays, binding


def declared_profile(directory, record, historical_parent):
    """Use bound record/pilot metadata, or the proven no-override historical v1 pilot defaults."""
    if isinstance(record.get("profile"), dict):
        return record["profile"], {"source": "parent_record.profile"}
    pilot_path = directory.parent/"pilot.json"
    if not pilot_path.is_file() or pilot_path.stat().st_mode & 0o222:
        raise ValueError("parent profile needs its immutable bound pilot metadata")
    pilot = json.loads(pilot_path.read_text())
    member = next((item for item in pilot["records"] if item["directory"] == directory.name), None)
    if member is None or member["parent_json_sha256"] != file_hash(directory/"parent.json") or \
            pilot["producing_commit"] != record["producing_commit"]:
        raise ValueError("pilot profile is not bound to this parent record")
    binding = {"path": str(pilot_path), "sha256": file_hash(pilot_path)}
    if isinstance(pilot.get("profile"), dict):
        return pilot["profile"], dict(binding, source="pilot.profile")
    if historical_parent is None:
        raise ValueError("historical parent defaults are unavailable")
    tree = ast.parse(historical_parent)
    pilot_function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "pilot")
    if any(arg.arg == "profile" for arg in pilot_function.args.args+pilot_function.args.kwonlyargs):
        raise ValueError("override-capable pilot omitted its profile; refusing a silent region default")
    constants = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in ("TRANSITION_END", "ENVELOPE_ZERO", "EXTERIOR_Q"):
                    constants[target.id] = ast.literal_eval(node.value)
    if set(constants) != {"TRANSITION_END", "ENVELOPE_ZERO", "EXTERIOR_Q"}:
        raise ValueError("historical pilot region constants are incomplete")
    return {"transition_end": constants["TRANSITION_END"], "envelope_zero": constants["ENVELOPE_ZERO"],
            "exterior_Q": constants["EXTERIOR_Q"]}, dict(binding, source="authenticated_no_override_v1_pilot_constants")


def summary(values, mask=None):
    array = np.asarray(values, dtype=float)
    if mask is not None:
        array = array[mask]
    if not array.size or not np.isfinite(array).all():
        raise ValueError("nonfinite or empty regional diagnostic")
    scale = float(np.max(abs(array)))
    return {"min": float(np.min(array)), "max": float(np.max(array)), "max_abs": scale,
            "rms": 0. if scale == 0 else scale*float(np.sqrt(np.mean((array/scale)**2)))}


def highband(grid, values):
    modes = np.fft.fftfreq(grid.nq)*grid.nq
    power = abs(np.fft.fft(values)/grid.nq)**2
    total = float(np.sum(power))
    fraction = lambda mask: 0. if total == 0 else float(np.sum(power[mask])/total)
    return {"power_including_zero_mode": total, "geometry_max_mode": grid.ng//2,
            "top_quarter_geometry_band_fraction": fraction(abs(modes) >= .75*(grid.ng//2)),
            "outside_geometry_band_fraction": fraction(abs(modes) > grid.ng//2),
            "above_common_mode31_fraction": fraction(abs(modes) > 31)}


def reconstruct(record, arrays):
    """Reconstruct only the stored v1 carrier and canonical state; no physical initial solve."""
    nf = int(record["nf"])
    grid = backend.make_fft_grid(galerkin.build_grid(nf, length=8., gauge="conformal"))
    grid.fine = replace(grid.fine, occupations=np.array(arrays["weights"], copy=True))
    if arrays["W"].shape != (grid.ng, grid.ng) or np.max(abs(arrays["W"].T@arrays["W"]-np.eye(grid.ng))) > 1e-9:
        raise ValueError("stored canonical frame does not span this carrier")
    if arrays["phi0"].shape != (nf, 2) or arrays["phi1"].shape != (nf, 2):
        raise ValueError("parent-v1 initial assessment requires its stored rank-two field")
    if not np.array_equal(np.vstack((arrays["phi0"], arrays["phi1"])), arrays["source_columns"]):
        raise ValueError("saved initial field is not its bound source")
    pair = SimpleNamespace(grid=grid, geometry_map=arrays["W"], weights=arrays["weights"],
                           reference_columns=arrays["reference_columns"], source_columns=arrays["source_columns"],
                           child_interval=tuple(record["intervals"]["child"]),
                           parent_interval=tuple(record["intervals"]["parent"]),
                           clock_locations=tuple(record["clock_locations"]))
    state = leading.State(*(arrays[key] for key in ("Q", "r", "pi_Q", "pi_r", "phi0", "phi1")))
    return pair, state


def measure(directory):
    record, arrays, binding = authenticate(directory)
    pair, state = reconstruct(record, arrays)
    grid = pair.grid
    tokens = {name: array_hash(value) for name, value in arrays.items()}
    rate, bundle = leading.rates(pair, state, return_bundle=True)
    fine, system, source = bundle["fine_state"], bundle["fine_system"], bundle["source"]
    curvature, velocity, acceleration = leading.metric_jets(pair, state, rate, bundle)
    C, D = leading.constraint_arrays(grid, fine, system, source)
    dx = lambda values: tidal.metric.spectral_dx(values, grid.length)
    def metric(v, acc):
        return tidal.curvature_from_local_jets(
            fine.Q, fine.r, v.Q, dx(fine.Q), acc.Q, dx(v.Q), dx(dx(fine.Q)),
            v.r, dx(fine.r), acc.r, dx(v.r), dx(dx(fine.r)))
    raw_velocity = leading.State(*bundle["unprojected_rates"])
    raw_acceleration = leading.fine_jvp(system, fine, raw_velocity)
    raw_curvature = metric(raw_velocity, raw_acceleration)
    no_final_projection = metric(velocity, leading.fine_jvp(system, fine, velocity))
    normal_radial = (velocity.r/fine.r+velocity.Q/fine.Q)/(fine.r*fine.Q)
    normal_angular = velocity.r/(fine.r**2*fine.Q)
    h, source_h = fine.Q*C, fine.Q*source["force_L"]/grid.dx_q
    sigma = abs((grid.xi_q-4.+grid.length/2)%grid.length-grid.length/2)
    transition = float(binding["declared_profile"]["transition_end"])
    envelope_zero = float(binding["declared_profile"]["envelope_zero"])
    masks = {"child": sigma <= .5, "collar_only": (sigma > .5)&(sigma <= 1.),
             "splice": (sigma > 1.)&(sigma <= transition),
             "source_annulus": (sigma > transition)&(sigma <= envelope_zero),
             "outer_parent": (sigma > envelope_zero)&(sigma <= 3.), "outside_parent": sigma > 3.}
    fields = {"r": fine.r, "Q": fine.Q, "p_Q_nodal": fine.p_Q, "p_r_nodal": fine.p_r,
              "normal_radial_H": normal_radial, "normal_angular_H": normal_angular,
              "h_c": h, "D": D, "source_Q_rho": source_h,
              "R_0101": curvature["tides"]["R_0101"], "R_0202": curvature["tides"]["R_0202"],
              "R4": curvature["tides"]["R4"]}
    projected = {key: galerkin.prolong_geometry(grid, galerkin.pull_geometry(grid, value)) for key, value in
                 (("h_c", h), ("D", D))}
    peaks = {key: {"x": float(grid.xi_q[np.argmax(abs(value))]),
                   "distance_from_centre": float(sigma[np.argmax(abs(value))]), "max_abs": float(np.max(abs(value)))}
             for key, value in fields.items()}
    factor = abs(system.C_W)/system.A
    riemann = np.maximum.reduce([abs(curvature["base"][key]) for key in ("A", "B00", "B11", "B01", "S")])
    raw_riemann = np.maximum.reduce([abs(raw_curvature["base"][key]) for key in ("A", "B00", "B11", "B01", "S")])
    trace = {"actual_projected_R4": summary(curvature["tides"]["R4"]),
             "projected_velocity_without_final_acceleration_projection_R4": summary(no_final_projection["tides"]["R4"]),
             "raw_same_point_ODE_R4": summary(raw_curvature["tides"]["R4"]),
             "raw_same_point_ODE_R_0101": summary(raw_curvature["tides"]["R_0101"]),
             "raw_same_point_ODE_R_0202": summary(raw_curvature["tides"]["R_0202"]),
             "final_acceleration_projection_R4_difference_max": float(np.max(abs(
                 curvature["tides"]["R4"]-no_final_projection["tides"]["R4"]))),
             "velocity_path_R4_difference_max": float(np.max(abs(
                 no_final_projection["tides"]["R4"]-raw_curvature["tides"]["R4"]))),
             "raw_reference_is_acceptance_or_replacement": False,
             "raw_reference_uses_same_SBP_hamiltonian": True,
             "continuum_tracefree_reference": "massless Dirac plus Maxwell gives R4=0 on shell",
             "metric_trace_variation_residual": summary(-8*np.pi*system.A*fine.r**3*fine.Q**2*curvature["tides"]["R4"])}
    clocks = leading.clock_rates(pair, state).tolist()
    saved_clocks = [row["N"] for row in record["diagnosis"]["proper_clocks"]]
    kinetic = float(np.mean(fine.p_Q**2/fine.r**3))
    gram = state.phi0.conj().T@state.phi0+state.phi1.conj().T@state.phi1
    roots = np.sqrt(pair.weights)
    result = {"nf": grid.nf, "binding": binding, "record_status": record["status"],
              "record_converged_flag": record["converged"], "k_common_declared": record["k"],
              "actual_kinetic_anchor": kinetic, "declared_kinetic_anchor": record["diagnosis"]["anchor"],
              "actual_anchor_difference": kinetic-record["diagnosis"]["anchor"],
              "regions": {name: {key: summary(value, mask) for key, value in fields.items()} for name, mask in masks.items()},
              "global": {key: summary(value) for key, value in fields.items()}, "peaks": peaks,
              "centre_values": {key: float(value[grid.nq//2]) for key, value in fields.items()},
              "analysis_region_distances": {"child": [0., .5], "collar_only": [.5, 1.],
                                            "splice": [1., transition], "source_annulus": [transition, envelope_zero],
                                            "outer_parent": [envelope_zero, 3.], "outside_parent": [3., 4.]},
              "highbands": {key: highband(grid, value) for key, value in fields.items()},
              "constraints": {"projected_metric": {key: summary(value) for key, value in projected.items()},
                              "source_h_max": float(np.max(abs(source_h))),
                              "full_h_over_source_h_max": float(np.max(abs(h))/np.max(abs(source_h))),
                              "current": summary(source["force_beta"]/grid.dx_q),
                              "observable_error_bound": None, "physical_admissibility_inferred": False},
              "trace_projection_diagnostic": trace,
              "curvature_benchmark": {"abs_CW_over_A": factor, "actual_physical_Riemann_ratio_max": factor*float(np.max(riemann)),
                                      "regional_actual_Riemann_ratios": {name: factor*float(np.max(riemann[mask])) for name, mask in masks.items()},
                                      "raw_reference_Riemann_ratio_max": factor*float(np.max(raw_riemann)),
                                      "CW_is_in_the_leading_action": False, "EFT_range_certified": False,
                                      "units": "abs(CW)/A length^2 times physical tetrad curvature length^-2; never raw Rh"},
              "proper_clock_locations": list(pair.clock_locations), "proper_clock_rates": clocks,
              "saved_clock_rate_replay_gap": float(np.max(abs(np.asarray(clocks)-saved_clocks))),
              "energies": leading.energy(pair, state), "CAR_nonzero_eigenvalues": np.linalg.eigvalsh(roots[:,None]*gram*roots[None,:]).tolist(),
              "source_and_clock_array_hashes_unchanged": tokens == {name: array_hash(value) for name, value in arrays.items()}}
    return result, fields


def report(directories):
    started = time.process_time()
    sources = {str(Path(module.__file__).resolve()): file_hash(module.__file__) for module in
               (backend, leading, tidal, coupling, galerkin, provenance)}
    sources[str(Path(__file__).resolve())] = file_hash(__file__)
    rows, profiles, input_hashes = [], [], {}
    with threadpool_limits(limits=1), backend.fft_thread_limit(1):
        for directory in directories:
            directory = Path(directory).expanduser().resolve()
            for filename in ("parent.json", "parent.npz"):
                input_hashes[str(directory/filename)] = file_hash(directory/filename)
            measured, profile = measure(directory)
            profile_binding = measured["binding"]["profile_binding"]
            if "path" in profile_binding:
                input_hashes[profile_binding["path"]] = profile_binding["sha256"]
            rows.append(measured); profiles.append(profile)
    comparisons = []
    for left, right, p, q in zip(rows[:-1], rows[1:], profiles[:-1], profiles[1:]):
        nl, nr = len(p["r"]), len(q["r"])
        if nr % nl:
            comparisons.append({"nf": [left["nf"], right["nf"]], "common_nodes_available": False})
            continue
        stride = nr//nl
        comparisons.append({"nf": [left["nf"], right["nf"]], "common_nodes_available": True,
                            "gaps": {key: summary(q[key][::stride]-p[key]) for key in p},
                            "refinement_error_certificate": False})
    if any(file_hash(path) != value for path, value in {**sources, **input_hashes}.items()):
        raise ValueError("input or evaluator changed during read-only assessment")
    return {"schema": SCHEMA, "rows": rows, "common_node_comparisons": comparisons,
            "input_hashes": input_hashes, "measurement_source_hashes": sources,
            "input_source_clock_bytes_unchanged": True, "evolved": False, "reprepared": False, "healed": False,
            "scope": "initial finite leading-Hamiltonian slices; projected curvature remains the actual readout; raw jets only diagnose Euler/projection defects",
            "threshold_gate": None, "physical_instability_claimed": False,
            "suggested_countercheck_only": "at the same resolution broaden the existing collar-to-exterior transition from 1.2 to 1.4 (<pi/2), preserving the sealed collar and actual Dirac extension; compare full h,D, high-pr power and trace defect with the separately repaired kinetic anchor before any evolution",
            "cpu_seconds": time.process_time()-started}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", nargs="+", type=Path, default=[DEFAULT/"nf64", DEFAULT/"nf128"])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.output and args.output.exists():
            raise FileExistsError("refusing to overwrite exact assessment output")
        if args.output and any(args.output.resolve().is_relative_to(path.resolve()) for path in args.records):
            raise ValueError("assessment output must stay outside its immutable input record directories")
        result = report(args.records)
        result["record_written"] = bool(args.output)
        payload = json.dumps(result, indent=2, allow_nan=False)+"\n"
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open("x") as stream: stream.write(payload)
            args.output.chmod(0o444)
        print(payload, end="")
        return 0
    except (ValueError, RuntimeError, FileNotFoundError, FileExistsError, KeyError) as error:
        print(type(error).__name__+": "+str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
