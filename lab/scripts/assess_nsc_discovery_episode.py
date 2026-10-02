#!/usr/bin/env python3
"""Independent consumer of the sealed discovery-episode v1 batch.

Reads the finished six-case chunks, their observation rows, the run-binding
envelope, and atlas v1. It does not import episode producers, evolve a state,
or rewrite sealed bytes. ``--write`` creates one new record outside the
episode directory and only after the sealed hashes match.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
LAB = ROOT / "lab"
DEFAULT_EPISODE = LAB / "results" / "development" / "nsc-discovery-episode-v1"
DEFAULT_ATLAS = LAB / "results" / "development" / "nsc-discovery-atlas-v1.json"

SCHEMA = "NSC-DISCOVERY-EPISODE-ASSESSMENT-v1"
PRODUCING_COMMIT = "1c9e7705c53e123497ffd40f94fb365e250e8e79"
PERIOD = 8.0
STATIONS = (1.0, 3.0)
HANDOFF = 0.3
TIME_TOL = 1e-8
LEDGER_KEYS = (
    "coordinate_fieldwork",
    "pressure_work",
    "lapse_work",
    "boundary_child",
    "boundary_parent",
)
DEFECT_KEYS = (
    "normal_energy_projection_defect",
    "projection_defect",
    "normal_energy_defect",
)
MISSING_ENERGY_ACCOUNTS = (
    "energy_Q",
    "energy_r",
    "energy_chi",
    "energy_by_source_column",
)
CHILD_WINDOW = (1.0, 3.0)
PARENT_WINDOW = (0.0, 4.0)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_hashed(path: Path) -> tuple[bytes, str]:
    data = path.read_bytes()
    return data, sha256_bytes(data)


def relative_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return str(path.resolve())


def relative_gap(left: float, right: float) -> float:
    scale = max(abs(left), abs(right))
    if scale == 0.0:
        return 0.0
    return abs(left - right) / scale


def first_gap_exceeding_step(times, gaps, steps):
    """First time at which a spatial gap exceeds the finer case's own step change.

    This is a scale comparison on the stored cadence. It is not an acceptance gate.
    """
    found = None
    for time, gap, step in zip(times, gaps, steps):
        if step is None or step == 0.0:
            continue
        if gap > step and found is None:
            found = float(time)
    return found


def recorded_defect(row: dict):
    for key in DEFECT_KEYS:
        if key in row:
            return row[key]
    return None


def trapezoid(times, values) -> float:
    total = 0.0
    for index in range(1, len(times)):
        dt = float(times[index] - times[index - 1])
        total += 0.5 * dt * (float(values[index]) + float(values[index - 1]))
    return total


def time_key(time: float) -> str:
    return f"{float(time):.1f}"


def row_at(rows: list[dict], time: float) -> dict:
    hits = [row for row in rows if abs(float(row["time"]) - float(time)) <= TIME_TOL]
    if len(hits) != 1:
        raise RuntimeError(f"expected one observation at {time}, found {len(hits)}")
    return hits[0]


def match_rows(coarse_rows: list[dict], fine_rows: list[dict]) -> list[tuple[dict, dict]]:
    """Match every recorded coordinate frame once, irrespective of file order."""
    if len(coarse_rows) != len(fine_rows):
        raise RuntimeError("refinement observation counts differ")
    pairs = [(row, row_at(fine_rows, row["time"])) for row in sorted(coarse_rows, key=lambda row: row["time"])]
    if len({float(left["time"]) for left, _ in pairs}) != len(pairs):
        raise RuntimeError("duplicate coarse observation times")
    return pairs


def json_equal(left, right) -> bool:
    if isinstance(left, float) or isinstance(right, float):
        return isinstance(left, (int, float)) and isinstance(right, (int, float)) and float(left) == float(right)
    return left == right


def occupation_trace(phi0, phi1, weights) -> float:
    column = np.sum(np.abs(phi0) ** 2 + np.abs(phi1) ** 2, axis=0)
    return float(np.sum(np.asarray(weights, dtype=float) * column))


def nodal_window(phi0, phi1, weights, start: float, end: float) -> float:
    """Box sum on x_i = i * period / n_f. Not the observer's trigonometric integral."""
    count = int(phi0.shape[0])
    coordinate = np.arange(count) * (PERIOD / count)
    density = np.sum(np.asarray(weights, dtype=float) * (np.abs(phi0) ** 2 + np.abs(phi1) ** 2), axis=1)
    return float(np.sum(density[(coordinate >= start) & (coordinate < end)]))


def load_json_bytes(data: bytes) -> dict:
    return json.loads(data.decode("utf-8"))


def load_jsonl_bytes(data: bytes) -> list[dict]:
    rows = [json.loads(line) for line in data.decode("utf-8").splitlines() if line.strip()]
    if not rows:
        raise RuntimeError("observation file is empty")
    return rows


def arrays_from_npz(data: bytes) -> dict:
    with np.load(io.BytesIO(data)) as payload:
        return {name: payload[name] for name in payload.files}


def periodic_integral(values, start: float, end: float) -> float:
    """Integral of the finite periodic interpolant, independently of producers."""
    values = np.asarray(values, dtype=float)
    coefficients = np.fft.fft(values) / len(values)
    waves = 2.0 * np.pi * np.fft.fftfreq(len(values), d=PERIOD / len(values))
    nonzero = waves != 0.0
    factors = (np.exp(1j * waves[nonzero] * end) - np.exp(1j * waves[nonzero] * start)) / (1j * waves[nonzero])
    return float((coefficients[0] * (end - start) + np.dot(coefficients[nonzero], factors)).real)


def independent_metric_content(arrays: dict) -> dict:
    nf = arrays["phi0"].shape[0]
    count = 4 * nf  # Declared quadrature carrier of this six-case batch.
    ng = len(arrays["Q"])
    if ng % 2 != 1 or count < ng:
        raise RuntimeError("expected odd geometry band inside the field grid")
    modes = np.rint(np.fft.fftfreq(ng) * ng).astype(int)

    def prolong(name):
        coarse = arrays["W"] @ arrays[name]
        coefficients = np.zeros(count, dtype=complex)
        coefficients[modes % count] = np.fft.fft(coarse) / ng
        return (np.fft.ifft(coefficients) * count).real

    q, r = prolong("Q"), prolong("r")
    radial = q * r
    clock_coefficients = np.fft.fft(radial) / count
    clock_modes = np.fft.fftfreq(count) * count
    clock_phase = np.exp(2j * np.pi * np.array([1., 2., 3.])[:, None] * clock_modes[None, :] / PERIOD)
    def prolong_columns(name):
        # Remove the AP half mode, interpolate periodic integer modes, restore it.
        demodulated = arrays[name] * np.exp(-1j * np.pi * np.arange(nf) / nf)[:, None]
        coefficients = np.zeros((count, arrays[name].shape[1]), dtype=complex)
        field_modes = np.rint(np.fft.fftfreq(nf) * nf).astype(int)
        coefficients[field_modes % count] = np.fft.fft(demodulated, axis=0) / nf
        return np.fft.ifft(coefficients, axis=0) * count * np.exp(1j * np.pi * np.arange(count) / count)[:, None] * np.sqrt(nf / count)

    phi0, phi1 = prolong_columns("phi0"), prolong_columns("phi1")
    mass = np.sum(arrays["source_weights"] * (np.abs(phi0) ** 2 + np.abs(phi1) ** 2), axis=1)
    return {
        "child_proper_length": periodic_integral(radial, *CHILD_WINDOW),
        "parent_proper_length": periodic_integral(radial, *PARENT_WINDOW),
        "child_probability": periodic_integral(mass, *CHILD_WINDOW) * count / PERIOD,
        "parent_probability": periodic_integral(mass, *PARENT_WINDOW) * count / PERIOD,
        "clock_rates": (clock_phase @ clock_coefficients).real.tolist(),
        "Q_min": float(np.min(q)),
        "r_min": float(np.min(r)),
    }


def chunk_paths(episode: Path, case_id: str) -> list[Path]:
    paths = sorted(episode.glob(f"{case_id}-[0-9][0-9][0-9][0-9][0-9][0-9].json"))
    if len(paths) != 3:
        raise RuntimeError(f"{case_id} has {len(paths)} chunks, expected 3")
    return paths


def assert_creation_path(path: Path, episode: Path, atlas_json: Path) -> Path:
    destination = path.resolve()
    sealed = episode.resolve()
    atlas = atlas_json.resolve()
    if destination == sealed or sealed in destination.parents:
        raise PermissionError(f"refusing to write inside the sealed episode: {destination}")
    if destination in {atlas, atlas.with_suffix(".npz")}:
        raise PermissionError(f"refusing to rewrite atlas: {destination}")
    if destination.exists():
        raise FileExistsError(destination)
    return destination


def write_exclusive(path: Path, report: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    blob = json.dumps(report, indent=2, sort_keys=True, allow_nan=False).encode("utf-8") + b"\n"
    descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    with os.fdopen(descriptor, "wb") as stream:
        stream.write(blob)


def _station_view(row: dict, arrays: dict | None, stability: dict | None) -> dict:
    weights = None if arrays is None else np.asarray(arrays["source_weights"], dtype=float)
    phi0 = None if arrays is None else arrays["phi0"]
    phi1 = None if arrays is None else arrays["phi1"]
    trace = None if arrays is None else occupation_trace(phi0, phi1, weights)
    chart = {} if stability is None else stability.get("chart") or {}
    curvature = {} if stability is None else stability.get("actual_metric_curvature") or {}
    chi = {} if stability is None else stability.get("chi_shell") or {}
    energy = {} if stability is None else stability.get("energy") or {}
    work = {} if stability is None else stability.get("work") or {}
    algebra = {} if stability is None else stability.get("constraint_algebra_residual") or {}
    defect = recorded_defect(row)
    independent = None if arrays is None else independent_metric_content(arrays)
    accounts = {
        "energy_field": row.get("energy_field"),
        "energy_gravity": row.get("energy_gravity"),
        "energy_total": row.get("energy_total"),
    }
    for name in MISSING_ENERGY_ACCOUNTS:
        accounts[name] = row[name] if name in row else None
    return {
        "time": row["time"],
        "sample_kind": row["sample_kind"],
        "child_probability": row["child_probability"],
        "parent_probability": row["parent_probability"],
        "occupation_trace": trace,
        "child_probability_over_trace": None if trace in (None, 0.0) else row["child_probability"] / trace,
        "parent_probability_over_trace": None if trace in (None, 0.0) else row["parent_probability"] / trace,
        "exterior_probability_over_trace": None if trace in (None, 0.0) else 1.0 - row["parent_probability"] / trace,
        "child_probability_per_proper_length": row["child_probability"] / row["child_proper_length"],
        "nodal_child_probability": None if arrays is None else nodal_window(phi0, phi1, weights, *CHILD_WINDOW),
        "nodal_parent_probability": None if arrays is None else nodal_window(phi0, phi1, weights, *PARENT_WINDOW),
        "nodal_minus_observer_child": None
        if arrays is None
        else nodal_window(phi0, phi1, weights, *CHILD_WINDOW) - row["child_probability"],
        "nodal_exterior_probability": None
        if trace is None
        else trace - nodal_window(phi0, phi1, weights, *PARENT_WINDOW),
        "child_proper_length": row["child_proper_length"],
        "parent_proper_length": row["parent_proper_length"],
        "parent_annulus_proper_length": row["parent_annulus_proper_length"],
        "length_sum_matches_parent": abs(
            row["child_proper_length"] + row["parent_annulus_proper_length"] - row["parent_proper_length"]
        )
        <= 1e-9,
        "localization_peak_x": row["localization_peak_x"],
        "child_normal_energy": row["child_normal_energy"],
        "parent_normal_energy": row["parent_normal_energy"],
        "R_h_max": row["R_h_max"],
        "weyl_C2_max": row["weyl_C2_max"],
        "coordinate_fieldwork_rate": row["coordinate_fieldwork"],
        "pressure_work_rate": row["pressure_work"],
        "lapse_work_rate": row["lapse_work"],
        "boundary_child_net": row["boundary_child"],
        "boundary_parent_net": row["boundary_parent"],
        "boundary_incoming": row["boundary_incoming"] if "boundary_incoming" in row else None,
        "boundary_outgoing": row["boundary_outgoing"] if "boundary_outgoing" in row else None,
        "energy_accounts": accounts,
        "normal_energy_projection_defect": defect,
        "independent_metric_content": independent,
        "chart_Q_min": chart.get("Q_min", None if independent is None else independent["Q_min"]),
        "chart_r_min": chart.get("r_min", None if independent is None else independent["r_min"]),
        "chart_admissible": chart.get("admissible"),
        "chart_clamped": chart.get("clamped"),
        "Q_dot_max": curvature.get("Q_dot_max"),
        "Q_ddot_max": curvature.get("Q_ddot_max") if "Q_ddot_max" in curvature else None,
        "L_dot_max": curvature.get("L_dot_max"),
        "R_h_minmax_l2": curvature.get("R_h"),
        "weyl_minmax_l2": curvature.get("weyl_C2"),
        "chi_proxy": chi.get("proxy"),
        "R_h_minus_chi_minus_2": chi.get("R_h_minus_chi_minus_2"),
        "actual_minus_proxy": chi.get("actual_minus_proxy"),
        "chi_substituted": curvature.get("chi_substituted"),
        "continuum_curvature_certified": curvature.get("continuum_curvature_certified"),
        "field_closure_error": energy.get("field_closure_error"),
        "normal_shell_max": energy.get("normal_shell_max"),
        "metric_coordinate_work": work.get("metric_coordinate_work"),
        "projection_forcing_norm": algebra.get("projection_forcing_norm"),
        "quadrature_forcing_norm": algebra.get("quadrature_forcing_norm"),
        "constraint_norm": algebra.get("constraint_norm"),
        "pi_Q_abs_max": None if arrays is None else float(np.max(np.abs(arrays["pi_Q"]))),
        "chi_coefficient_abs_max": None if arrays is None else float(np.max(np.abs(arrays["chi"]))),
        "clock_rates": None if arrays is None else [float(value) for value in arrays["clock_rates"]],
        "normal_clocks": None if arrays is None else [float(value) for value in arrays["normal_clocks"]],
        "stability_certificate": None if stability is None else stability.get("stability_certificate"),
        "coupled_stability_classification": None
        if stability is None
        else stability.get("coupled_stability_classification"),
    }


def _pair_station(coarse: dict, fine: dict, names: tuple[str, ...]) -> dict:
    compared = {}
    for name in names:
        left, right = coarse.get(name), fine.get(name)
        if isinstance(left, (int, float)) and isinstance(right, (int, float)):
            compared[name] = {
                "coarse": left,
                "fine": right,
                "absolute_gap": abs(float(left) - float(right)),
                "relative_gap": relative_gap(float(left), float(right)),
            }
        else:
            compared[name] = {"coarse": left, "fine": right, "absolute_gap": None, "relative_gap": None}
    return compared


def _cadence(coarse_rows: list[dict], fine_rows: list[dict]) -> dict:
    series = []
    previous = None
    for row, other in match_rows(coarse_rows, fine_rows):
        step = None if previous is None else abs(other["R_h_max"] - previous["R_h_max"])
        weyl_step = None if previous is None else abs(other["weyl_C2_max"] - previous["weyl_C2_max"])
        gap = abs(row["R_h_max"] - other["R_h_max"])
        weyl_gap = abs(row["weyl_C2_max"] - other["weyl_C2_max"])
        series.append(
            {
                "time": row["time"],
                "R_h_gap": gap,
                "R_h_relative_gap": relative_gap(row["R_h_max"], other["R_h_max"]),
                "fine_R_h_step": step,
                "R_h_gap_over_fine_step": None if not step else gap / step,
                "weyl_gap": weyl_gap,
                "weyl_relative_gap": relative_gap(row["weyl_C2_max"], other["weyl_C2_max"]),
                "fine_weyl_step": weyl_step,
                "weyl_gap_over_fine_step": None if not weyl_step else weyl_gap / weyl_step,
                "child_length_relative_gap": relative_gap(row["child_proper_length"], other["child_proper_length"]),
                "child_probability_relative_gap": relative_gap(row["child_probability"], other["child_probability"]),
                "fine_child_probability": other["child_probability"],
                "fine_child_proper_length": other["child_proper_length"],
                "fine_R_h_max": other["R_h_max"],
                "fine_weyl_C2_max": other["weyl_C2_max"],
            }
        )
        previous = other
    return {
        "samples": series,
        "first_time_R_h_relative_gap_increases": next(
            (right["time"] for left, right in zip(series, series[1:]) if right["R_h_relative_gap"] > left["R_h_relative_gap"]), None
        ),
        "first_time_weyl_relative_gap_increases": next(
            (right["time"] for left, right in zip(series, series[1:]) if right["weyl_relative_gap"] > left["weyl_relative_gap"]), None
        ),
        "first_time_R_h_gap_exceeds_fine_step": first_gap_exceeding_step(
            [item["time"] for item in series],
            [item["R_h_gap"] for item in series],
            [item["fine_R_h_step"] for item in series],
        ),
        "first_time_weyl_gap_exceeds_fine_step": first_gap_exceeding_step(
            [item["time"] for item in series],
            [item["weyl_gap"] for item in series],
            [item["fine_weyl_step"] for item in series],
        ),
    }


def _clock_interval(start: dict, end: dict) -> dict:
    start_clock = np.asarray(start["normal_clocks"], dtype=float)
    end_clock = np.asarray(end["normal_clocks"], dtype=float)
    return {
        "coordinate_dt": float(end["time"] - start["time"]),
        "proper_dt": [float(value) for value in (end_clock - start_clock)],
        "coordinate_time_end": float(end["time"]),
        "clock_rates_end": end["clock_rates"],
        "normal_clocks_end": end["normal_clocks"],
        "observer_coordinates": [1.0, 2.0, 3.0],
        "scope": "normal clocks at fixed coordinate observers; coordinate stations are not matched proper times",
    }


def _control_effect(coupled_initial: dict, coupled: dict, frozen_initial: dict, frozen: dict, names) -> dict:
    """Difference of movements, retaining each control's own handoff value."""
    result = {}
    for name in names:
        coupled_change = float(coupled[name] - coupled_initial[name])
        frozen_change = float(frozen[name] - frozen_initial[name])
        effect = coupled_change - frozen_change
        result[name] = {
            "coupled_initial": coupled_initial[name], "frozen_initial": frozen_initial[name],
            "coupled_final": coupled[name], "frozen_final": frozen[name],
            "coupled_change": coupled_change, "frozen_change": frozen_change,
            "initial_subtracted_effect": effect,
            "effect_over_larger_movement": None if max(abs(coupled_change), abs(frozen_change)) == 0 else effect / max(abs(coupled_change), abs(frozen_change)),
        }
    return result


def _atlas_handoff(atlas: dict) -> dict:
    series = next(
        item
        for item in atlas["series"]
        if item["preparation"] == "separated_pair_v1" and item["case"] == "nf256_baseline_dt0.0005"
    )
    rows = [row for row in series["rows"] if row.get("time") is not None]
    end = rows[-1]
    previous = rows[-2]
    dt = float(end["time"] - previous["time"])
    energy_rate = (end["child_normal_energy"] - previous["child_normal_energy"]) / dt
    flux = 0.5 * (end["child_boundary_flux"] + previous["child_boundary_flux"])
    pressure = 0.5 * (end["child_pressure_work"] + previous["child_pressure_work"])
    lapse = 0.5 * (end["child_lapse_work"] + previous["child_lapse_work"])
    defect = recorded_defect(end)
    residual = energy_rate - flux - pressure - lapse
    return {
        "time": end["time"],
        "R_h_max": end["R_h_max"],
        "Q_dot_max": end["Q_dot_max"],
        "Q_ddot_max": end["Q_ddot_max"],
        "chi_shell_gap_max": end["chi_shell_gap_max"],
        "child_normal_energy": end["child_normal_energy"],
        "child_boundary_flux": end["child_boundary_flux"],
        "frozen_Q_ddot_max": series["final_frozen_geometry"]["Q_ddot_max"],
        "normal_energy_projection_defect": defect,
        "child_balance_residual_without_recorded_defect": residual,
        "shell_closure_claimed": False,
    }


def _scale_indicator(station: dict) -> dict:
    q_min = station.get("chart_Q_min")
    r_h = station.get("R_h_max")
    if not isinstance(q_min, (int, float)) or not isinstance(r_h, (int, float)) or q_min == 0.0:
        return {"available": False}
    return {
        "available": True,
        "colocated": False,
        "meaning": (
            "2/Q_min^2 and the parenthesis that would yield R_h_max at that same Q_min. "
            "The chart minimum and the curvature maximum are not a stored common node."
        ),
        "prefactor_2_over_Q_min_squared": 2.0 / float(q_min) ** 2,
        "bracket_if_R_h_max_sat_at_Q_min": float(r_h) * float(q_min) ** 2 / 2.0,
        "pi_Q_abs_max": station.get("pi_Q_abs_max"),
        "chi_coefficient_abs_max": station.get("chi_coefficient_abs_max"),
    }


def assess(episode: Path | None = None, atlas_json: Path | None = None) -> dict:
    episode = Path(episode or DEFAULT_EPISODE)
    atlas_json = Path(atlas_json or DEFAULT_ATLAS)
    atlas_npz = atlas_json.with_suffix(".npz")
    read_hashes: dict[str, str] = {}
    problems: list[str] = []
    case_rows: dict[str, list[dict]] = {}

    def take(path: Path) -> bytes:
        data, digest = read_hashed(path)
        read_hashes[relative_path(path)] = digest
        return data

    manifest = load_json_bytes(take(episode / "manifest.json"))
    ledger = load_json_bytes(take(episode / "cpu-ledger.json"))
    binding = load_json_bytes(take(episode / "run-binding-envelope.json"))
    closure = load_json_bytes(take(episode / "authenticated-closure.json"))
    atlas = load_json_bytes(take(atlas_json))
    atlas_payload, atlas_payload_hash = read_hashed(atlas_npz)
    read_hashes[relative_path(atlas_npz)] = atlas_payload_hash
    del atlas_payload
    if atlas.get("payload_sha256") != atlas_payload_hash:
        problems.append("atlas payload hash does not match nsc-discovery-atlas-v1.npz")
    if binding.get("producing_commit") != PRODUCING_COMMIT:
        problems.append("producing commit is not 1c9e7705c53e123497ffd40f94fb365e250e8e79")
    if closure.get("producing_commit") != PRODUCING_COMMIT or closure.get("completed") != "STATION_REACHED":
        problems.append("batch closure does not bind a completed producing commit")

    case_ids = [item["case_id"] for item in manifest["results"]]
    if len(case_ids) != 6 or len(set(case_ids)) != 6:
        problems.append("manifest does not record six distinct cases")
    cpu_sum = float(sum(item["child_cpu_seconds"] for item in manifest["results"]))
    if abs(cpu_sum - float(manifest["child_cpu_seconds"])) > 1e-9 or abs(cpu_sum - float(ledger["spent"])) > 1e-9:
        problems.append("child CPU sum does not match the manifest and ledger")

    cases: dict[str, dict] = {}
    chunk_count = 0
    for case_id in case_ids:
        observation_path = episode / f"{case_id}-observations.jsonl"
        rows = load_jsonl_bytes(take(observation_path))
        case_rows[case_id] = rows
        chronological = all(rows[index]["time"] <= rows[index + 1]["time"] for index in range(len(rows) - 1))
        ordered = sorted(rows, key=lambda row: float(row["time"]))
        chunks = []
        for json_path in chunk_paths(episode, case_id):
            record = load_json_bytes(take(json_path))
            npz_path = episode / record["npz"] if "npz" in record else json_path.with_suffix(".npz")
            if not npz_path.is_file():
                npz_path = json_path.with_suffix(".npz")
            payload, payload_hash = read_hashed(npz_path)
            read_hashes[relative_path(npz_path)] = payload_hash
            chunk_count += 1
            if record.get("immutable") is not True:
                problems.append(f"{json_path.name} is not marked immutable")
            if record.get("arrays_sha256") != payload_hash:
                problems.append(f"{npz_path.name} hash does not match its json")
            before = record.get("source_pins_before") or {}
            after = record.get("source_pins_after") or {}
            for pin in ("W_sha256", "observer_columns_sha256", "source_columns_sha256", "weights_sha256"):
                if before.get(pin) != after.get(pin) or not before.get(pin):
                    problems.append(f"{json_path.name} {pin} is not stable")
            if after.get("phi_sha256") == after.get("source_columns_sha256"):
                problems.append(f"{json_path.name} stores Phi as the source columns")
            if "clock_rates_sha256" in after and after.get("clock_rates_sha256") == after.get("weights_sha256"):
                problems.append(f"{json_path.name} stores clock rates as the source weights")
            if "channel_sample" in record:
                sample = record["channel_sample"]
                matched = row_at(ordered, record["coordinate_time"])
                for key, value in sample.items():
                    if key not in matched or not json_equal(matched[key], value):
                        problems.append(f"{case_id} station {record['coordinate_time']} row disagrees on {key}")
                        break
                if matched.get("sample_kind") != "station":
                    problems.append(f"{case_id} coordinate time {record['coordinate_time']} is not a station row")
            arrays = arrays_from_npz(payload)
            chunks.append({"record": record, "arrays": arrays})
        by_station = {}
        for time in (HANDOFF, *STATIONS):
            record = next(item["record"] for item in chunks if abs(float(item["record"]["coordinate_time"]) - time) <= TIME_TOL)
            arrays = next(item["arrays"] for item in chunks if item["record"] is record)
            row = row_at(ordered, time)
            by_station[time_key(time)] = _station_view(row, arrays, record.get("stability"))
            independent = by_station[time_key(time)]["independent_metric_content"]
            for name in ("child_proper_length", "parent_proper_length", "child_probability", "parent_probability"):
                if abs(independent[name] - row[name]) > 1e-9 * max(1.0, abs(row[name])):
                    problems.append(f"{case_id} at {time} independent {name} disagrees with row")
            if not np.allclose(independent["clock_rates"], arrays["clock_rates"], atol=1e-9, rtol=1e-9):
                problems.append(f"{case_id} at {time} independent clock rates disagree with payload")
        integrals = {key: trapezoid([row["time"] for row in ordered], [row[key] for row in ordered]) for key in LEDGER_KEYS}
        final = next(item["record"] for item in chunks if abs(float(item["record"]["coordinate_time"]) - STATIONS[-1]) <= TIME_TOL)
        discrepancies = {
            key: abs(integrals[key] - float(final["work_ledger"][key])) for key in LEDGER_KEYS
        }
        if max(discrepancies.values()) > 1e-6:
            problems.append(f"{case_id} sample trapezoid does not reproduce the work ledger")
        trace_values = [view["occupation_trace"] for view in by_station.values()]
        if any(abs(value - 3.0) > 1e-8 for value in trace_values):
            problems.append(f"{case_id} occupation trace is not 3")
        cases[case_id] = {
            "control_mode": final["control_mode"],
            "nf": final["nf"],
            "step_cap": final["step_cap"],
            "stations_reached": final.get("stations_reached"),
            "stability_certificate": final.get("stability_certificate"),
            "observation_count": len(ordered),
            "observation_file_order_is_chronological": chronological,
            "stations": by_station,
            "coordinate_work_integrals": integrals,
            "work_ledger": {key: final["work_ledger"][key] for key in LEDGER_KEYS},
            "work_integral_abs_discrepancy": discrepancies,
            "normal_energy_projection_defect": recorded_defect(ordered[-1]),
        }

    if chunk_count != 18:
        problems.append(f"sealed chunk count is {chunk_count}, expected 18")

    space_coarse = "nf128_coupled_dt0.0005"
    space_fine = "nf256_coupled_dt0.0005"
    time_coarse = "nf256_coupled_dt0.001"
    compared_names = (
        "child_probability",
        "parent_probability",
        "child_probability_over_trace",
        "parent_probability_over_trace",
        "exterior_probability_over_trace",
        "child_probability_per_proper_length",
        "child_proper_length",
        "parent_proper_length",
        "R_h_max",
        "weyl_C2_max",
        "chart_Q_min",
        "chart_r_min",
        "child_normal_energy",
        "parent_normal_energy",
        "Q_dot_max",
        "Q_ddot_max",
    )
    space_stations = {}
    time_stations = {}
    for time in (HANDOFF, *STATIONS):
        key = time_key(time)
        space_stations[key] = _pair_station(cases[space_coarse]["stations"][key], cases[space_fine]["stations"][key], compared_names)
        time_stations[key] = _pair_station(cases[time_coarse]["stations"][key], cases[space_fine]["stations"][key], compared_names)
        for label, station in (("space", space_stations[key]), ("time", time_stations[key])):
            for moment in ("R_h_minmax_l2", "weyl_minmax_l2", "chi_proxy"):
                coarse_moment = cases[space_coarse if label == "space" else time_coarse]["stations"][key][moment]
                fine_moment = cases[space_fine]["stations"][key][moment]
                if isinstance(coarse_moment, dict) and isinstance(fine_moment, dict):
                    station[moment] = {
                        part: {
                            "coarse": coarse_moment.get(part),
                            "fine": fine_moment.get(part),
                            "relative_gap": None
                            if coarse_moment.get(part) is None or fine_moment.get(part) is None
                            else relative_gap(coarse_moment[part], fine_moment[part]),
                        }
                        for part in ("min", "max", "l2")
                    }

    cadence = _cadence(case_rows[space_coarse], case_rows[space_fine])
    atlas_tie = _atlas_handoff(atlas)
    handoff = cases[space_fine]["stations"]["0.3"]
    atlas_tie["episode_R_h_max_absolute_difference"] = abs(handoff["R_h_max"] - atlas_tie["R_h_max"])
    atlas_tie["episode_child_normal_energy_absolute_difference"] = abs(
        handoff["child_normal_energy"] - atlas_tie["child_normal_energy"]
    )

    clocks = {}
    for case_id in (space_fine, "nf256_frozen_geometry_dt0.0005", space_coarse, "nf128_frozen_geometry_dt0.0005"):
        ordered_views = [cases[case_id]["stations"][time_key(time)] for time in (HANDOFF, *STATIONS)]
        clocks[case_id] = {
            f"{time_key(ordered_views[index]['time'])}->{time_key(ordered_views[index + 1]['time'])}": _clock_interval(
                ordered_views[index], ordered_views[index + 1]
            )
            for index in range(len(ordered_views) - 1)
        }

    fine_late = cases[space_fine]["stations"]["3.0"]
    frozen_late = cases["nf256_frozen_geometry_dt0.0005"]["stations"]["3.0"]
    fine_handoff = cases[space_fine]["stations"]["0.3"]
    control_effects = {}
    control_names = (
        "child_probability", "parent_probability", "child_probability_over_trace",
        "parent_probability_over_trace", "exterior_probability_over_trace",
        "child_probability_per_proper_length", "child_proper_length", "parent_proper_length",
        "chart_Q_min", "chart_r_min", "R_h_max", "weyl_C2_max",
        "child_normal_energy", "parent_normal_energy",
    )
    for nf in (128, 256):
        coupled = cases[f"nf{nf}_coupled_dt0.0005"]["stations"]
        frozen = cases[f"nf{nf}_frozen_geometry_dt0.0005"]["stations"]
        control_effects[str(nf)] = {
            key: _control_effect(coupled["0.3"], coupled[key], frozen["0.3"], frozen[key], control_names)
            for key in ("1.0", "3.0")
        }
    control_refinement = {}
    for key in ("1.0", "3.0"):
        control_refinement[key] = {
            name: {
                "coarse_effect": control_effects["128"][key][name]["initial_subtracted_effect"],
                "fine_effect": control_effects["256"][key][name]["initial_subtracted_effect"],
                "relative_gap": relative_gap(control_effects["128"][key][name]["initial_subtracted_effect"], control_effects["256"][key][name]["initial_subtracted_effect"]),
            } for name in control_names
        }
    for path_text, digest in read_hashes.items():
        path = Path(ROOT, path_text)
        if path.parent.resolve() == episode.resolve() and path.name != "authenticated-closure.json":
            if closure.get("artifacts", {}).get(path.name) != digest:
                problems.append(f"{path.name} does not match the authenticated batch closure")
    changed = []
    for path_text, digest in read_hashes.items():
        again = sha256_bytes(Path(ROOT, path_text).read_bytes())
        if again != digest:
            changed.append(path_text)
    if changed:
        problems.append("sealed input hashes changed while they were being read")

    defects = [cases[case_id]["normal_energy_projection_defect"] for case_id in case_ids]
    incoming_present = any(
        view["boundary_incoming"] is not None or view["boundary_outgoing"] is not None
        for case in cases.values()
        for view in case["stations"].values()
    )
    report = {
        "schema": SCHEMA,
        "ok": not problems,
        "problems": problems,
        "producing_commit": binding.get("producing_commit"),
        "producing_commit_abbreviation": "1c9e770",
        "consumer_sha256": sha256_bytes(Path(__file__).read_bytes()),
        "read_hashes": read_hashes,
        "hashes_unchanged": not changed,
        "sealed_bytes_rewritten": False,
        "n_cases": len(case_ids),
        "n_chunks": chunk_count,
        "cpu_seconds": cpu_sum,
        "cpu_accounting": manifest.get("cpu_accounting"),
        "claims": {
            "stability_theorem": False,
            "regeneration_claimed": False,
            "shell_closure_claimed": False,
            "normal_energy_closure_claimed": False,
            "continuum_curvature_certified": False,
            "one_percent_gate_installed": False,
            "normal_energy_projection_defect": None if all(item is None for item in defects) else defects,
            "defect_filled_with_zero": False,
            "incoming_outgoing_separated": incoming_present,
            "incoming_outgoing_split": "recorded" if incoming_present else "needed",
            "analytic_qddot_closes_T3_spatial_gap": False,
            "coordinate_station_is_common_proper_time": False,
            "continuous_evolution_error_bound": None,
            "tidal_response_computed_by_this_consumer": False,
        },
        "definitions": {
            "relative_gap": "|a-b|/max(|a|,|b|)",
            "gap_over_fine_step": "spatial |coarse-fine| divided by the fine case change over the previous 0.05 cadence step",
            "occupation_trace": "sum of stored source weights times nodal |phi0|^2+|phi1|^2",
            "nodal_window": "period-8 box sum on x_i=i*8/n_f; the stored child and parent probabilities remain the observer integrals",
            "normal_energy_projection_defect": "null when the sample row has no defect field",
            "initial_subtracted_control_effect": "(coupled(T)-coupled(0.3))-(frozen(T)-frozen(0.3)), at matched coordinate time",
            "pressure_lapse_scope": "stored rate sums child [1,3] and parent [0,4], overlapping windows; not per-child or a global disjoint energy ledger",
            "curvature_scale": "prefactor and bracket are not a co-located jet evaluation",
        },
        "physical_outcome": {
            "interpretation": (
                "Coupled evolution through T=3 contracts the coordinate-window proper lengths, "
                "raises the areal-radius floor, and moves probability out of the child window. "
                "Late R_h and Weyl maxima have substantial spatial sensitivity; these scalars alone do not classify tidal response. This is not a regeneration "
                "or stability theorem."
            ),
            "nf256_coupled_child_proper_length": {
                key: cases[space_fine]["stations"][key]["child_proper_length"] for key in ("0.3", "1.0", "3.0")
            },
            "nf256_coupled_parent_proper_length": {
                key: cases[space_fine]["stations"][key]["parent_proper_length"] for key in ("0.3", "1.0", "3.0")
            },
            "nf256_coupled_child_probability": {
                key: cases[space_fine]["stations"][key]["child_probability"] for key in ("0.3", "1.0", "3.0")
            },
            "nf256_coupled_child_probability_over_trace": {
                key: cases[space_fine]["stations"][key]["child_probability_over_trace"] for key in ("0.3", "1.0", "3.0")
            },
            "nf256_coupled_localization_peak_x": {
                key: cases[space_fine]["stations"][key]["localization_peak_x"] for key in ("0.3", "1.0", "3.0")
            },
            "occupation_trace_nf256_coupled_T3": fine_late["occupation_trace"],
            "areal_r_min_frozen_T3": frozen_late["chart_r_min"],
            "areal_r_min_coupled_T3": fine_late["chart_r_min"],
            "areal_expansion_factor_r_min": fine_late["chart_r_min"] / frozen_late["chart_r_min"],
            "proper_length_contraction_factor": fine_late["child_proper_length"] / fine_handoff["child_proper_length"],
            "chart_Q_min_coupled_T3": fine_late["chart_Q_min"],
            "chart_Q_min_frozen_T3": frozen_late["chart_Q_min"],
            "T1_R_h_max_space_relative_gap": space_stations["1.0"]["R_h_max"]["relative_gap"],
            "T3_R_h_max_space_relative_gap": space_stations["3.0"]["R_h_max"]["relative_gap"],
            "T3_R_h_l2_space_relative_gap": space_stations["3.0"]["R_h_minmax_l2"]["l2"]["relative_gap"],
            "T3_weyl_max_space_relative_gap": space_stations["3.0"]["weyl_C2_max"]["relative_gap"],
            "T3_Q_min_space_relative_gap": space_stations["3.0"]["chart_Q_min"]["relative_gap"],
            "T3_r_min_space_relative_gap": space_stations["3.0"]["chart_r_min"]["relative_gap"],
            "T3_child_length_space_relative_gap": space_stations["3.0"]["child_proper_length"]["relative_gap"],
            "T3_R_h_max_time_relative_gap": time_stations["3.0"]["R_h_max"]["relative_gap"],
            "late_R_h_gap_exceeds_fine_step_at": cadence["first_time_R_h_gap_exceeds_fine_step"],
            "late_weyl_gap_exceeds_fine_step_at": cadence["first_time_weyl_gap_exceeds_fine_step"],
            "stability_certificate": False,
        },
        "cases": cases,
        "space_refinement": {
            "coarse": space_coarse,
            "fine": space_fine,
            "stations": space_stations,
            "cadence": cadence,
        },
        "time_refinement": {
            "larger_step": time_coarse,
            "smaller_step": space_fine,
            "stations": time_stations,
        },
        "clocks": clocks,
        "coupled_vs_frozen": {
            "comparison_scope": "matched coordinate stations; differences are not evaluated at matched proper clocks",
            "stations": control_effects,
            "space_refinement_of_initial_subtracted_effect": control_refinement,
        },
        "atlas_handoff": atlas_tie,
        "curvature_scales_T3": {
            "nf256_coupled": _scale_indicator(fine_late),
            "nf128_coupled": _scale_indicator(cases[space_coarse]["stations"]["3.0"]),
            "coupled_Q_ddot_max": fine_late["Q_ddot_max"],
            "frozen_Q_ddot_max": frozen_late["Q_ddot_max"],
            "frozen_R_h_max": frozen_late["R_h_max"],
            "atlas_handoff_Q_ddot_max": atlas_tie["Q_ddot_max"],
        },
        "balance": {
            "episode_boundary_samples": "net",
            "incoming_outgoing_split": "needed",
            "normal_energy_projection_defect": None,
            "sample_trapezoid_matches_coordinate_work_ledger": all(
                max(case["work_integral_abs_discrepancy"].values()) <= 1e-6 for case in cases.values()
            ),
            "max_work_integral_abs_discrepancy": max(
                value for case in cases.values() for value in case["work_integral_abs_discrepancy"].values()
            ),
            "shell_closure_claimed": False,
            "constraint_projection_forcing_is_a_different_ledger": True,
            "normal_energy_closure_claimed": False,
            "atlas_child_balance_residual_without_recorded_defect": atlas_tie[
                "child_balance_residual_without_recorded_defect"
            ],
        },
        "next_experiment": "Combine this coordinate-station assessment with the separately computed native-Jv radial/angular tidal response, then choose matched proper-clock comparisons and a finite-projection-defect normal-energy ledger before extending numerical trajectories.",
    }
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Assess the sealed discovery-episode v1 batch")
    parser.add_argument("--check", action="store_true", help="authenticate inputs and print the assessment")
    parser.add_argument("--write", type=Path, default=None, help="create one new JSON record outside the episode")
    parser.add_argument("--episode", type=Path, default=DEFAULT_EPISODE)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    args = parser.parse_args(argv)
    if not args.check and args.write is None:
        parser.error("choose --check, --write, or both")
    destination = None
    if args.write is not None:
        destination = assert_creation_path(args.write, args.episode, args.atlas)
    report = assess(args.episode, args.atlas)
    if args.check:
        json.dump(report, sys.stdout, indent=2, sort_keys=True, allow_nan=False)
        sys.stdout.write("\n")
    if destination is not None:
        if not report["ok"] or not report["hashes_unchanged"]:
            raise RuntimeError("refusing to write before the sealed inputs authenticate")
        write_exclusive(destination, report)
        print(f"wrote {destination}", file=sys.stderr)
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, PermissionError, FileExistsError) as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1)
