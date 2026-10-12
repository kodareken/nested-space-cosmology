#!/usr/bin/env python3
"""Physical-window ledgers from finalized conformal episode Cauchy frames.

No trajectory evolves here. Exact integrals of the periodic trigonometric
interpolants over [0,2],[2,4],[4,6],[6,8] define width-2 windows. Boundary
fluxes, normal shell pressure/lapse work and directional balance defects
use the active L=Q system and actual projected lifted rate. The legacy
FQ-only coordinate work is deliberately not called for this dynamic lapse.
"""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import sys
import time

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
from dataclasses import replace
import derive_nsc_spherical_conformal_episode as episode
from recursive_horizons import nsc_regional_energy_exchange as regional
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons.nsc_spherical_coupling import CauchyState, _combine, source_from_columns

LAB = Path(__file__).resolve().parents[1]
OUT = LAB / "results/development/nsc-spherical-conformal-frames-v1.json"
NPZ = LAB / "results/development/nsc-spherical-conformal-frames-v1.npz"
DRIVER = Path(__file__).resolve()
BOUNDS = (0., 2., 4., 6., 8.)


def interval_weights(count, length, start, end):
    """Integral of a periodic nodal trigonometric interpolant, divided by dx."""
    if not 0 <= start < end <= length:
        raise ValueError("ordered interval inside one period required")
    wave = 2 * np.pi * np.fft.fftfreq(count, d=length / count)
    integrals = np.empty(count, dtype=complex)
    nonzero = wave != 0
    integrals[nonzero] = (np.exp(1j * wave[nonzero] * end) - np.exp(1j * wave[nonzero] * start)) / (1j * wave[nonzero])
    integrals[~nonzero] = end - start
    return (np.fft.fft(integrals) / length).real


def stage_normal_slope(system, fine, rate, epsilon=1e-6):
    """Independent derivative of FL/r with BOTH state and lapse varied."""
    forward, backward = _combine(fine, rate, epsilon), _combine(fine, rate, -epsilon)
    forward_system = replace(system, length_density=forward.Q, shift=np.zeros_like(forward.Q))
    backward_system = replace(system, length_density=backward.Q, shift=np.zeros_like(backward.Q))
    plus = source_from_columns(forward_system, forward)["force_L"] / forward.r
    minus = source_from_columns(backward_system, backward)["force_L"] / backward.r
    return (plus - minus) / (2 * epsilon)


def frame_ledger(grid, state, time_value):
    coarse, bundle = galerkin.compose_fine_hamiltonian(grid, state)
    fine, system, source = bundle["fine_state"], bundle["fine_system"], bundle["source"]
    lifted = episode.legacy.lifted_rate(grid, coarse, bundle)
    ledger = regional.matter_ledger(system, fine)
    terms = regional.proper_balance_terms(system, fine, lifted, ledger)
    geometry = regional.geometric_ledger(system, fine, lifted)
    dynamic_work = source["force_Q"] * lifted.Q + source["force_L"] * bundle["lapse_dot"]
    normal_slope = stage_normal_slope(system, fine, lifted)
    continuity_residual = normal_slope + system.derivative @ terms["flux_nodal"] - terms["source_nodal"]
    weights = [interval_weights(grid.nq, grid.length, start, end) for start, end in zip(BOUNDS[:-1], BOUNDS[1:])]
    probability = np.sum((abs(fine.phi0) ** 2 + abs(fine.phi1) ** 2) * system.occupations[None, :], axis=1)
    probability_flux = 2 * np.sum(np.imag(np.conjugate(fine.phi0) * fine.phi1) * system.occupations[None, :], axis=1)
    phi0_t, phi1_t = lifted.phi0, lifted.phi1
    probability_slope = 2 * np.sum(np.real(np.conjugate(fine.phi0) * phi0_t + np.conjugate(fine.phi1) * phi1_t) * system.occupations[None, :], axis=1)
    probability_residual = probability_slope + system.derivative @ probability_flux
    total_probability = float(np.sum(probability))
    circular = np.sum(probability * np.exp(2j * np.pi * grid.xi_q / grid.length)) / total_probability
    rows = []
    for index, window in enumerate(weights):
        start, end = BOUNDS[index:index + 2]
        a = int(round(start / grid.dx_q)) % grid.nq
        b = int(round(end / grid.dx_q)) % grid.nq
        boundary_flux = float((terms["flux_nodal"][a] - terms["flux_nodal"][b]) / grid.dx_q)
        discrete_flux = float(np.dot(system.derivative @ window, terms["flux_nodal"]))
        probability_boundary = float((probability_flux[a] - probability_flux[b]) / grid.dx_q)
        inside = (grid.xi_q >= start) & (grid.xi_q < end)
        rows.append({
            "interval": [start, end], "width": end - start,
            "normal_shell": float(np.dot(window, ledger["normal_energy_nodal"])),
            "coordinate_matter_energy": float(np.dot(window, ledger["hamiltonian_integrand"])),
            "coordinate_gravity_energy": float(np.dot(window, geometry["nodal_energy"])),
            "normal_boundary_flux": boundary_flux,
            "normal_discrete_boundary_flux": discrete_flux,
            "boundary_flux_identity_gap": abs(boundary_flux - discrete_flux),
            "normal_left_outward_flux": float(-terms["flux_nodal"][a] / grid.dx_q),
            "normal_right_outward_flux": float(terms["flux_nodal"][b] / grid.dx_q),
            "proper_pressure_work": float(np.dot(window, terms["proper_pressure_work"])),
            "momentum_lapse_work": float(np.dot(window, terms["momentum_lapse_work"])),
            "coordinate_metric_work": float(np.dot(window, dynamic_work)),
            "normal_directional_slope": float(np.dot(window, normal_slope)),
            "normal_balance_defect": float(np.dot(window, continuity_residual)),
            "probability": float(np.dot(window, probability)),
            "probability_share": float(np.dot(window, probability) / total_probability),
            "probability_boundary_flux": probability_boundary,
            "probability_directional_slope": float(np.dot(window, probability_slope)),
            "probability_balance_defect": float(np.dot(window, probability_residual)),
            "normal_internal_flux_l2": float(np.sqrt(max(0., np.dot(window, terms["flux_nodal"] ** 2) / grid.dx_q))),
            "normal_internal_flux_max_abs": float(np.max(abs(terms["flux_nodal"][inside])) / grid.dx_q),
            "probability_internal_flux_l2": float(np.sqrt(max(0., np.dot(window, probability_flux ** 2) / grid.dx_q))),
            "probability_internal_flux_max_abs": float(np.max(abs(probability_flux[inside])) / grid.dx_q),
        })
    shell = np.array([row["normal_shell"] for row in rows])
    shell_leader = int(np.argmax(shell))
    probability_leader = int(np.argmax([row["probability"] for row in rows]))
    fields = {"normal_energy_nodal": ledger["normal_energy_nodal"], "normal_flux_nodal": terms["flux_nodal"],
              "proper_pressure_work": terms["proper_pressure_work"], "momentum_lapse_work": terms["momentum_lapse_work"],
              "normal_balance_defect": continuity_residual, "probability_nodal": probability,
              "probability_flux_nodal": probability_flux, "probability_balance_defect": probability_residual}
    return {
        "time": float(time_value), "windows": rows, "shell_leader": shell_leader, "probability_leader": probability_leader,
        "shell_leader_share": float(shell[shell_leader] / np.sum(shell)),
        "probability_leader_share": rows[probability_leader]["probability_share"],
        "canonical_probability": total_probability,
        "circular_probability_location": float(np.mod(np.angle(circular), 2 * np.pi) * grid.length / (2 * np.pi)),
        "circular_concentration": float(abs(circular)),
        "dynamic_work_gap": float(abs(np.sum(dynamic_work) - coarse.fieldwork_power)),
        "normal_balance_defect_l2_density": float(np.sqrt(np.sum(continuity_residual ** 2) / grid.dx_q)),
        "probability_balance_defect_l2_density": float(np.sqrt(np.sum(probability_residual ** 2) / grid.dx_q)),
        "partition_gap": float(np.max(abs(np.sum(weights, axis=0) - 1.))),
        "boundary_flux_identity_gap_max": max(row["boundary_flux_identity_gap"] for row in rows),
        "full_normal_energy": float(np.sum(ledger["normal_energy_nodal"])),
        "source_kernel_gap_max": max(ledger["kernel_gap"].values()),
        "shell_identity_gap_max": float(np.max(abs(ledger["stresses"]["shell_gap"]))),
    }, fields


def summarize(rows):
    first, last = rows[0], rows[-1]
    dt = float(rows[1]["time"] - rows[0]["time"])
    windows = []
    for index in range(4):
        interval = [row["windows"][index] for row in rows]
        channels = {name: episode.integrate([row[name] for row in interval], dt) for name in
                    ("normal_boundary_flux", "proper_pressure_work", "momentum_lapse_work", "coordinate_metric_work", "normal_balance_defect", "probability_boundary_flux", "probability_balance_defect")}
        shell_change = interval[-1]["normal_shell"] - interval[0]["normal_shell"]
        probability_change = interval[-1]["probability"] - interval[0]["probability"]
        predicted_shell = sum(channels[name]["trapezoid"] for name in ("normal_boundary_flux", "proper_pressure_work", "momentum_lapse_work", "normal_balance_defect"))
        predicted_probability = channels["probability_boundary_flux"]["trapezoid"] + channels["probability_balance_defect"]["trapezoid"]
        windows.append({"interval": interval[0]["interval"], "normal_shell_change": shell_change, "probability_change": probability_change,
                        "integrals": channels, "normal_shell_integrated_gap": shell_change - predicted_shell,
                        "probability_integrated_gap": probability_change - predicted_probability,
                        "min_normal_shell": min(row["normal_shell"] for row in interval),
                        "max_abs_normal_boundary_flux": max(abs(row["normal_boundary_flux"]) for row in interval),
                        "max_internal_normal_flux": max(row["normal_internal_flux_max_abs"] for row in interval),
                        "max_internal_probability_flux": max(row["probability_internal_flux_max_abs"] for row in interval),
                        "max_abs_normal_throughflow": max(max(abs(row["normal_left_outward_flux"]), abs(row["normal_right_outward_flux"])) for row in interval)})
    return {"initial": first, "final": last, "windows": windows,
            "shell_leader_unchanged": all(row["shell_leader"] == first["shell_leader"] for row in rows),
            "probability_leader_unchanged": all(row["probability_leader"] == first["probability_leader"] for row in rows),
            "normal_shell_leader_share_min": min(row["shell_leader_share"] for row in rows),
            "probability_leader_share_min": min(row["probability_leader_share"] for row in rows),
            "probability_total_drift": max(abs(row["canonical_probability"] - first["canonical_probability"]) for row in rows),
            "max_dynamic_work_gap": max(row["dynamic_work_gap"] for row in rows),
            "max_boundary_flux_identity_gap": max(row["boundary_flux_identity_gap_max"] for row in rows),
            "normal_continuity_defect_l2_max": max(row["normal_balance_defect_l2_density"] for row in rows),
            "probability_continuity_defect_l2_max": max(row["probability_balance_defect_l2_density"] for row in rows),
            "renewal": False, "continuum_balance_certified": False}


def run():
    episode.refuse_existing_outputs(OUT, NPZ)
    start_cpu = time.process_time()
    core = json.loads(episode.OUT.read_text())
    original_hashes = {"episode_json": episode.sha256(episode.OUT), "episode_npz": episode.sha256(episode.NPZ)}
    if not core["full_target_reached"] or core["gauge"] != "conformal":
        raise ValueError("finalized conformal full-target episode required")
    record = {"schema": "NSC-SPHERICAL-CONFORMAL-FRAMES-v1", "gauge": "conformal", "evolution_performed": False,
              "source_trajectory_hashes": original_hashes, "bound_sources": {"driver": episode.sha256(DRIVER), "regional": episode.sha256(regional._MODULE_PATH), "galerkin": episode.sha256(galerkin._MODULE_PATH)},
              "normal_shell": "FL/r nodal integral; pressure/lapse ledger evaluated at stage-local L=Q",
              "windows": [list(pair) for pair in zip(BOUNDS[:-1], BOUNDS[1:])],
              "window_quadrature": "exact integrals of periodic nodal trigonometric interpolants, width 2 each; positive continuum density not certified",
              "balance_directional_test": "central derivative of FL/r with L=Q varied and actual projected fine state/rate; defects retained",
              "time_quadrature": ".005 saved frames, Simpson/trapezoid indicators; no certified integral", "results": {}, "series": {}}
    arrays = {}
    with np.load(episode.NPZ, allow_pickle=False) as payload:
        for nf in (256, 512):
            grid = galerkin.build_grid(nf, quadrature=4 * nf, gauge="conformal")
            for spec in (row for row in episode.RUN_PLAN if row["nf"] == nf):
                times = payload[spec["name"] + "_frame_times"]
                frames = {name: payload[spec["name"] + "_frames_" + name] for name in episode.STATE_NAMES}
                rows, fields = [], {}
                for index, t in enumerate(times):
                    if core["cpu_seconds"] + time.process_time() - start_cpu > episode.CPU_BUDGET - 5:
                        raise RuntimeError("shared diagnostic CPU budget would be exceeded")
                    state = CauchyState(**{name: frames[name][index] for name in episode.STATE_NAMES})
                    row, profiles = frame_ledger(grid, state, t)
                    rows.append(row)
                    for name, value in profiles.items():
                        fields.setdefault(name, []).append(value)
                record["results"][spec["name"]] = summarize(rows)
                record["series"][spec["name"]] = rows
                arrays[spec["name"] + "_times"] = times
                for name, value in fields.items():
                    arrays[spec["name"] + "_" + name] = np.array(value)
                print("frames", spec["name"], "shell_leader", rows[-1]["shell_leader"], "shell_change", record["results"][spec["name"]]["windows"][0]["normal_shell_change"], "flux", rows[-1]["windows"][0]["normal_boundary_flux"], flush=True)
    record["cpu_seconds"] = time.process_time() - start_cpu
    record["combined_episode_and_frames_cpu"] = record["cpu_seconds"] + core["cpu_seconds"]
    record["source_trajectory_unchanged"] = original_hashes == {"episode_json": episode.sha256(episode.OUT), "episode_npz": episode.sha256(episode.NPZ)}
    total_before = episode.OUT.stat().st_size + episode.NPZ.stat().st_size
    with NPZ.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    record["payload_sha256"] = episode.sha256(NPZ)
    OUT.write_text(json.dumps(episode.jsonable(record), indent=2, allow_nan=False) + "\n")
    record["combined_payload_bytes"] = total_before + OUT.stat().st_size + NPZ.stat().st_size
    if record["combined_payload_bytes"] > episode.PAYLOAD_LIMIT:
        raise RuntimeError("combined payload exceeds 64MiB")
    OUT.write_text(json.dumps(episode.jsonable(record), indent=2, allow_nan=False) + "\n")
    print("saved", OUT, "combined_cpu", record["combined_episode_and_frames_cpu"], "combined_payload", record["combined_payload_bytes"], flush=True)
    return record


def verify_saved(*, source_ref=None):
    before = {path: episode.sha256(path) for path in (OUT, NPZ)}
    record = json.loads(OUT.read_text())
    if record.get("schema") != "NSC-SPHERICAL-CONFORMAL-FRAMES-v1" or record.get("payload_sha256") != episode.sha256(NPZ):
        raise ValueError("frame schema or payload differs")
    episode.check_recorded_sources(record.get("source_trajectory_hashes"),
        {"episode_json": episode.OUT, "episode_npz": episode.NPZ})
    episode.check_recorded_sources(record.get("bound_sources"),
        {"driver": DRIVER, "regional": regional._MODULE_PATH, "galerkin": galerkin._MODULE_PATH}, source_ref=source_ref)
    with np.load(NPZ, allow_pickle=False) as payload:
        for spec in episode.RUN_PLAN:
            name = spec["name"]
            rows = record["series"][name]
            episode.check_saved_values(summarize(rows), record["results"][name], name)
            normal = payload[name + "_normal_energy_nodal"]
            flux = payload[name + "_normal_flux_nodal"]
            for index, (start, end) in enumerate(zip(BOUNDS[:-1], BOUNDS[1:])):
                weight = interval_weights(normal.shape[1], 8., start, end)
                episode.check_saved_values(normal @ weight,
                    [row["windows"][index]["normal_shell"] for row in rows], name + ".shell")
                left, right = index * normal.shape[1] // 4, ((index + 1) * normal.shape[1] // 4) % normal.shape[1]
                net = (flux[:, left] - flux[:, right]) / (8. / normal.shape[1])
                episode.check_saved_values(net,
                    [row["windows"][index]["normal_boundary_flux"] for row in rows], name + ".boundary")
    if before != {path: episode.sha256(path) for path in before}:
        raise RuntimeError("read-only check changed frame bytes")
    return {"status": "VERIFIED_STORED_CONFORMAL_FRAMES", "wrote": False, "source_ref": source_ref}


if __name__ == "__main__":
    if sys.argv[1:] == ["--check"]:
        print(json.dumps(verify_saved(source_ref=episode.SEALED_SOURCE_REF)))
    elif len(sys.argv) == 1:
        run()
    else:
        raise SystemExit("Use --check for sealed evidence, or no arguments to create new outputs")
