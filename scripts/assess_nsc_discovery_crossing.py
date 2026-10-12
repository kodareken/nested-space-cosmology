#!/usr/bin/env python3
"""Read-only T8 metric-jet assessment of the authenticated crossing campaign."""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import resource
import time

for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

import numpy as np
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_discovery_episode as episode
from recursive_horizons import nsc_discovery_response as response
from recursive_horizons import nsc_discovery_tidal as tidal
from recursive_horizons import nsc_nested_parent_child as nested
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin

REPO = episode.REPO
DIRECTORY = episode.LAB / "results/development/nsc-discovery-crossing-v1"
PREDECESSOR = episode.LAB / "results/development/nsc-discovery-episode-v1"
OUTPUT = episode.LAB / "results/development/nsc-discovery-crossing-assessment-v1.json"
SCHEMA = "NSC-DISCOVERY-CROSSING-ASSESSMENT-v1"
CROSSING_COMMIT = "b7f0dab8de62d9a3472ac4e7b33351ab65d25847"
CPU_BUDGET = 30.0
CHANNELS = ("R_0101", "R_0202", "Ricci2", "K", "R_h", "R4", "owned_W")


def relative_path(path):
    return str(Path(path).resolve().relative_to(REPO.resolve()))


def bound_path(name):
    """All newly recorded source/input references are portable repository paths."""
    name = Path(name)
    if name.is_absolute():
        raise ValueError("bound paths must be repository relative")
    result = (REPO / name).resolve()
    result.relative_to(REPO.resolve())
    return result


def hashes(paths):
    return {relative_path(path): tidal.sha256_file(path) for path in paths}


def authenticate(directory=DIRECTORY, predecessor=PREDECESSOR):
    """Verify the crossing's post-run observation and original physics envelope."""
    path = Path(directory) / "observed-run-binding.json"
    crossing = json.loads(path.read_text())
    if crossing.get("producing_commit") != CROSSING_COMMIT:
        raise ValueError("unexpected crossing producing commit")
    artifacts = crossing.get("artifact_paths") or {}
    producers = crossing.get("producer_paths") or {}
    if not artifacts or not producers:
        raise ValueError("crossing binding lacks artifacts or producer closure")
    for group in (artifacts, producers):
        for name, declaration in group.items():
            if tidal.sha256_file(bound_path(name)) != declaration["sha256"]:
                raise ValueError("crossing binding hash changed: " + name)
    original = tidal.frozen_binding(predecessor, "episode")
    original = {
        **original, "envelope": relative_path(original["envelope"]),
        "physics_hashes": {relative_path(name): digest for name, digest in original["physics_hashes"].items()},
    }
    return {
        "crossing": {"envelope": relative_path(path), "sha256": tidal.sha256_file(path),
                     "producing_commit": crossing["producing_commit"],
                     "binding_timing": crossing["binding_timing"],
                     "artifact_hashes": {name: value["sha256"] for name, value in artifacts.items()},
                     "producer_hashes": {name: value["sha256"] for name, value in producers.items()}},
        "original_episode": original,
    }


def validate_handoff(handoff, original):
    """Reject any reset, mismatched predecessor pointer, or partial hash dictionary."""
    pointer = handoff.get("predecessor_chunk") or {}
    left, right = handoff.get("array_sha256"), original.get("array_sha256")
    if not left or not right or left != right:
        raise ValueError("T3 handoff array hash dictionaries differ")
    if (handoff.get("case_id") != original.get("case_id")
            or pointer.get("case_id") != original.get("case_id")
            or pointer.get("arrays_sha256") != original.get("arrays_sha256")
            or int(pointer.get("ordinal", -1)) != int(original.get("ordinal", -2))
            or abs(float(handoff.get("coordinate_time", -1)) - 3.0) > 1e-8
            or abs(float(original.get("coordinate_time", -1)) - 3.0) > 1e-8):
        raise ValueError("T3 handoff predecessor or time mismatch")
    return {"all_array_hashes_equal": True, "coordinate_time": 3.0,
            "predecessor_arrays_sha256": pointer["arrays_sha256"],
            "no_state_momentum_clock_basis_reset": True}


def match_observations(left, right):
    """Join physical times within 1e-7; never pair rows by their positions."""
    def keyed(rows):
        result = {}
        for row in rows:
            key = round(float(row["time"]), 7)
            if key in result:
                raise ValueError("duplicate observation timestamp")
            result[key] = row
        return result
    a, b = keyed(left), keyed(right)
    common = sorted(a.keys() & b.keys())
    pairs = [(a[key], b[key]) for key in common]
    return pairs, {"matched_count": len(common), "left_count": len(left), "right_count": len(right),
                   "left_only_times": sorted(a.keys() - b.keys()), "right_only_times": sorted(b.keys() - a.keys()),
                   "join": "coordinate time rounded to 1e-7; no positional zip"}


def cadence_gaps(rows):
    times = sorted(float(row["time"]) for row in rows)
    return [{"left": a, "right": b, "duration": b-a}
            for a, b in zip(times, times[1:]) if b-a > 1.5 * episode.OBSERVATION_CADENCE]


def comparison(first, second):
    result = {}
    for channel in CHANNELS:
        a, b = first[channel], second[channel]
        shared = min(a.size, b.size)
        if a.size % shared or b.size % shared:
            raise ValueError("profiles do not share periodic nodes")
        a, b = a[::a.size // shared], b[::b.size // shared]
        delta = b-a
        world = int(round(2 * shared / 8))
        result[channel] = {
            "shape_max_relative_to_second": float(np.max(abs(delta)) / max(np.max(abs(b)), 1e-300)),
            "shape_l2_relative_to_second": float(np.linalg.norm(delta) / max(np.linalg.norm(b), 1e-300)),
            "worldline_relative_to_second": float(abs(delta[world]) / max(abs(b[world]), 1e-300)),
        }
    return result


def measure_endpoint(directory, case_id):
    record, arrays = episode.load_checkpoint(directory, case_id)
    if abs(record["coordinate_time"] - 8.0) > 1e-8:
        raise ValueError("crossing endpoint must be the saved T8 station")
    pair = episode.pair_from_arrays(arrays, record)
    state = episode.state_from_arrays(arrays, record["momentum_representation"])
    jets = tidal.analytic_accelerations(pair, state, record["control_mode"])
    profile = tidal.profiles_on_grid(pair, jets)
    observed = tidal.station_observables(pair, jets, arrays["normal_clocks"])
    tides, base = profile["tides"], profile["base"]
    q, r = tides["Q"], tides["r"]
    N = q*r
    dx = lambda values: tidal.metric.spectral_dx(values, pair.grid.length)
    qt, rt, qtt, rtt = (jets[name] for name in ("Q_t", "r_t", "Q_tt", "r_tt"))
    qx, rx = dx(q), dx(r)
    radial_terms = np.array([rtt/r, -rt**2/r**2, qtt/q, -qt**2/q**2,
                             -dx(rx)/r, rx**2/r**2, -dx(qx)/q, qx**2/q**2]) / N**2
    world = int(round(2/pair.grid.dx_q))
    columns = np.vstack((state.phi0, state.phi1))
    gram = columns.conj().T @ columns
    eigen = np.linalg.eigvalsh(np.sqrt(pair.weights[:, None]*pair.weights[None, :])*gram)
    logs = episode.read_observations(directory, case_id)
    matches = [row for row in logs if abs(row["time"] - record["coordinate_time"]) <= 1e-7]
    if len(matches) != 1:
        raise ValueError("endpoint requires exactly one matching observation")
    last = matches[0]
    row = {
        "case_id": case_id, "T": record["coordinate_time"], "arrays_sha256": record["arrays_sha256"],
        "Q_min": float(q.min()), "Q_max": float(q.max()), "r_min": float(r.min()), "r_max": float(r.max()),
        "N_min": float(N.min()), "N_max": float(N.max()), "positive_metric": bool(np.min(N) > 0),
        "clocks": arrays["normal_clocks"].tolist(), "clock_rates": arrays["clock_rates"].tolist(),
        "gram_gap": float(np.max(abs(gram-np.eye(6)))),
        "occupied_covariance_eigenvalues": eigen.tolist(), "remaining_covariance_eigenvalues": "zero",
        "child_probability": last["child_probability"], "child_proper_length": last["child_proper_length"],
        "localization_width": last["localization_width"],
        "regional_stocks": {name: last[name] for name in
                            ("child_probability", "parent_probability", "child_normal_energy", "parent_normal_energy")},
        "worldline_x2": observed["worldline_x2"], "conditioning": observed["conditioning"],
        "radial_condition_x2": float(np.sum(abs(radial_terms[:, world])) / max(abs(tides["R_0101"][world]), 1e-300)),
        "component_maxima": {name: observed["summaries"][name]["max_abs"] for name in CHANNELS},
        "stored_scalar_diagnostics": {name: last[name] for name in ("R_h_max", "weyl_C2_max", "aux_discrepancy_max")},
        "chi_substituted_for_curvature": False, "cadence_gaps": cadence_gaps(logs),
    }
    if record["control_mode"] == "coupled":
        rate = nested.rates(pair, state)
        tangent = response.StateTangent(*(getattr(rate, name) for name in nested.STATE_NAMES), np.zeros(6))
        jv = response.nested_rate_jacobian_vector(pair, state, tangent)
        row["full_Jv_geometry_relative_gaps"] = {
            name: float(np.max(abs(galerkin.prolong_geometry(pair.grid, pair.geometry_map @ getattr(jv, name))
                                   - jets[name+"_tt"])) / max(np.max(abs(jets[name+"_tt"])), 1e-300))
            for name in ("Q", "r")
        }
    channels = {name: tides[name] if name in tides else base[name] for name in CHANNELS}
    return row, channels


def assess(directory=DIRECTORY, predecessor=PREDECESSOR):
    started = time.process_time()
    bindings = authenticate(directory, predecessor)
    watched_inputs = list(Path(directory).glob("*.json")) + list(Path(directory).glob("*.npz"))
    watched_inputs += list(Path(directory).glob("*.jsonl"))
    watched_inputs += list(Path(predecessor).glob("*.json")) + list(Path(predecessor).glob("*.npz"))
    watched_inputs += list(Path(predecessor).glob("*-observations.jsonl"))
    sources = list(bindings["crossing"]["producer_hashes"])
    sources += list(bindings["original_episode"]["physics_hashes"])
    sources += [relative_path(path) for path in tidal.producer_hashes()]
    sources += [relative_path(response.__file__), "lab/src/recursive_horizons/nsc_causal_common.py",
                relative_path(__file__), "lab/tests/test_nsc_discovery_crossing.py"]
    source_before = hashes(bound_path(path) for path in set(sources))
    input_before = hashes(watched_inputs)
    rows, profiles, handoffs = [], {}, {}
    with threadpool_limits(limits=1):
        for spec in episode.first_batch_specs():
            case_id = spec["case_id"]
            handoff, _ = episode.load_checkpoint(directory, case_id, 0)
            pointer = handoff.get("predecessor_chunk") or {}
            original, _ = episode.load_checkpoint(predecessor, case_id, int(pointer.get("ordinal", -1)))
            handoffs[case_id] = validate_handoff(handoff, original)
            row, channels = measure_endpoint(directory, case_id)
            rows.append(row)
            profiles[case_id] = channels
            if time.process_time()-started > CPU_BUDGET:
                raise RuntimeError("crossing assessment exceeded thirty CPU seconds")
    comparisons = {}
    for cap in ("0.001", "0.0005"):
        comparisons["nf128_to_nf256_dt"+cap] = comparison(
            profiles["nf128_coupled_dt"+cap], profiles["nf256_coupled_dt"+cap])
    for nf in (128, 256):
        comparisons["dt0.001_to_dt0.0005_nf"+str(nf)] = comparison(
            profiles[f"nf{nf}_coupled_dt0.001"], profiles[f"nf{nf}_coupled_dt0.0005"])
    history, onset, matches = {}, {}, {}
    for nf in (128, 256):
        left = episode.read_observations(directory, f"nf{nf}_coupled_dt0.001")
        right = episode.read_observations(directory, f"nf{nf}_coupled_dt0.0005")
        pairs, matches[str(nf)] = match_observations(left, right)
        history[str(nf)] = [{
            "T": b["time"],
            "Rh_max_step_relative": abs(a["R_h_max"]-b["R_h_max"])/max(abs(b["R_h_max"]), 1e-300),
            "W_max_step_relative": abs(a["weyl_C2_max"]-b["weyl_C2_max"])/max(abs(b["weyl_C2_max"]), 1e-300),
            "aux_gap": b["aux_discrepancy_max"],
            "aux_gap_over_Rh_max": b["aux_discrepancy_max"]/max(abs(b["R_h_max"]), 1e-300),
        } for a, b in pairs]
        combined = episode.read_observations(predecessor, f"nf{nf}_coupled_dt0.0005") + right
        combined.sort(key=lambda row: row["time"])
        onset[str(nf)] = {}
        for ratio in (0.001, 0.01, 0.1, 1.0):
            hit = next((row for row in combined if row["aux_discrepancy_max"] /
                        max(abs(row["R_h_max"]), 1e-300) > ratio), None)
            onset[str(nf)][str(ratio)] = None if hit is None else {
                "T": hit["time"], "aux_gap": hit["aux_discrepancy_max"], "Rh_max": hit["R_h_max"]}
    source_after = hashes(bound_path(path) for path in source_before)
    input_after = hashes(bound_path(path) for path in input_before)
    if source_after != source_before or input_after != input_before:
        raise ValueError("assessment source or input changed during consumption")
    cpu = time.process_time()-started
    if cpu > CPU_BUDGET:
        raise RuntimeError("crossing assessment exceeded thirty CPU seconds")
    return {
        "schema": SCHEMA, "bindings": bindings, "T3_handoffs": handoffs, "rows": rows,
        "comparisons": comparisons, "scalar_history": history, "timestamp_matches": matches,
        "auxiliary_gap_onset": onset,
        "onset_protocol": "Illustrative first-sample aux-gap/Rh-max ratios, not acceptance thresholds or physical curvature breakpoints",
        "scope": {
            "actual_tidal_agreement": "supported on the tested nf128/256 grids and both saved timestep caps; no continuum certificate",
            "R4_Rh_Weyl_resolved": False,
            "return": "carrier-period return candidate; frozen geometry also returns",
            "maintained_physical_scale": False, "ambient_echo_proven": False,
            "regional_stock_change_is_transport_closure": False, "singularity_declared": False,
            "geometry_or_memory_reset": False, "new_ode_or_trajectory": False,
        },
        "profile_comparison_protocol": "Fine arrays restricted to shared coarse periodic nodes; norms relative to the second case",
        "source_hashes": source_after, "input_hashes": input_after, "input_sources_unchanged": True,
        "cpu_seconds": cpu, "cpu_budget_seconds": CPU_BUDGET, "numerical_threads": 1,
        "code_API": "episode.load_checkpoint -> pair_from_arrays/state_from_arrays -> tidal.analytic_accelerations -> profiles_on_grid/station_observables; independent response.nested_rate_jacobian_vector along nested.rates",
    }


def write_record(report, path=OUTPUT):
    destination = Path(path).expanduser().resolve()
    if destination != OUTPUT.resolve():
        raise ValueError("creation is restricted to the crossing-assessment-v1 owner")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(report, indent=2, allow_nan=False)+"\n")


def _stable_values(report):
    return np.array([value for row in report["rows"] for value in (
        row["worldline_x2"]["R_0101"], row["worldline_x2"]["R_0202"],
        row["Q_min"], row["Q_max"], row["r_min"], row["r_max"],
        row["N_min"], row["N_max"], *row["clocks"], *row["clock_rates"],
        row["child_probability"], row["child_proper_length"], *row["occupied_covariance_eigenvalues"],
        *(row["component_maxima"][name] for name in ("R_0101", "R_0202", "Ricci2", "K")),
    )], dtype=float)


def check_record(path=OUTPUT):
    record = json.loads(Path(path).read_text())
    if record.get("schema") != SCHEMA:
        raise ValueError("unexpected crossing assessment schema")
    for group in ("source_hashes", "input_hashes"):
        if not record.get(group):
            raise ValueError("record lacks source/input hashes")
        for name, digest in record[group].items():
            if tidal.sha256_file(bound_path(name)) != digest:
                raise ValueError("bound hash changed: "+name)
    current = assess()
    before, after = _stable_values(record), _stable_values(current)
    # Numerical replay tolerance only; no new scientific acceptance threshold.
    if before.shape != after.shape or not np.allclose(before, after, rtol=1e-10, atol=1e-20):
        raise ValueError("stable endpoint quantities do not reproduce")
    if record.get("scope") != current["scope"] or record.get("T3_handoffs") != current["T3_handoffs"]:
        raise ValueError("assessment scope or handoff authentication changed")
    return {"ok": True, "recomputed_metric_jets": True, "evolved": False, "bytes_written": 0,
            "cpu_seconds": current["cpu_seconds"], "scalar_replay_is_not_a_physical_certificate": True}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--write", action="store_true", help="create the fixed assessment owner once")
    modes.add_argument("--check", action="store_true", help="authenticate and recompute metric jets without evolution")
    args = parser.parse_args(argv)
    if args.write and OUTPUT.exists():
        raise FileExistsError("refusing to overwrite "+str(OUTPUT))
    # CLI process CPU cap, preserving a stricter inherited limit.
    soft, hard = resource.getrlimit(resource.RLIMIT_CPU)
    budget_end = math.ceil(time.process_time()+CPU_BUDGET)
    if hard != resource.RLIM_INFINITY:
        budget_end = min(budget_end, hard)
    if soft == resource.RLIM_INFINITY or soft > budget_end:
        resource.setrlimit(resource.RLIMIT_CPU, (budget_end, hard))
    if args.check:
        print(json.dumps(check_record(), indent=2))
        return 0
    report = assess()
    if args.write:
        write_record(report)
    print(json.dumps({"schema": report["schema"], "cpu_seconds": report["cpu_seconds"],
                      "scope": report["scope"], "rows": report["rows"],
                      "comparisons": report["comparisons"], "auxiliary_gap_onset": report["auxiliary_gap_onset"]},
                     indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
