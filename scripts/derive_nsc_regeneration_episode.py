#!/usr/bin/env python3
"""One adaptive continuation from the saved T=0.05 spherical episode.

The column source, the radius, the shift momentum, and ``rk4_step`` stay
with their owners. This driver does not add a force, a pulse, a mean
subtraction, or a state reset. It does not write the episode, the controls,
the null-expansion record, or any sealed predecessor.

Declared before the kept steps. The coarse pilot is n_f=256 with dt at most
0.0005, reduced each step to the owned ``stable_timestep``. The fine
confirmation cap is 0.00025. A safety-cap crossing changes the next dt. It
is not a physical-failure label. An accepted step is rejected when r, Q, or
L leave the owned chart, the column Gram leaves 1e-8, or dt*omega exceeds
the absolute RK4 limit 2*sqrt(2). The comparison window is T=0.15, or the
stored 0.005 frame of a declared geometric event. Shell reversal of 0.10 and
a leader share of 0.50 are reported proxies and are not stop rules.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import sys
import time
import traceback
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import numpy as np

try:
    from threadpoolctl import threadpool_limits

    threadpool_limits(limits=1).__enter__()
    THREADPOOL_LIMIT = 1
except Exception as error:
    THREADPOOL_LIMIT = "not applied: " + type(error).__name__

import derive_nsc_spherical_feedback_episode as episode
from recursive_horizons import nsc_regional_energy_exchange as regional
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons import nsc_spherical_null_expansion as null_expansion
from recursive_horizons.provenance import resolve_pinned_source_bytes
from recursive_horizons.nsc_regeneration_controls import (
    CAR_GAP_MAX,
    assess_windows,
    car_report,
    chart_from_fields,
    join_series,
    load_episode_final,
    proper_clock_rates,
    separate_ledgers,
    shell_views,
    state_sha256,
)
from recursive_horizons.nsc_spherical_coupling import (
    PERIOD,
    CauchyRate,
    PositiveChartExit,
    chart_failure,
    locked_coefficients,
)
from recursive_horizons.nsc_spherical_galerkin_coupling import physical_link_report

LAB = episode.LAB
ROOT = LAB.parent
OUT = LAB / "results" / "development" / "nsc-regeneration-episode-v1.json"
NPZ = LAB / "results" / "development" / "nsc-regeneration-episode-v1.npz"
NOTE = LAB / "docs" / "nsc-regeneration-episode.md"
DRIVER = Path(__file__).resolve()
OWNED_OUTPUTS = frozenset(path.resolve() for path in (OUT, NPZ, NOTE))

PROTOCOL_ID = "nsc-regeneration-episode-v1"
SCHEMA = "NSC-REGENERATION-EPISODE-v1"
T0 = 0.05
TARGET = 0.15
FRAME = 0.005
DT_COARSE = 5e-4
DT_CONFIRM = 2.5e-4
OWNED_DT_CAP = 5e-4
SAFETY_CAP = 1.4
ABSOLUTE_RK4 = 2.0 * math.sqrt(2.0)
COARSE_BUDGET = 600.0
FINE_BUDGET = 1800.0
PAYLOAD_LIMIT = 64 * 1024 * 1024
PACKET_END = 4.0
BRIDGE_START = 4.0
Q_TURNBACK_FLOOR = 1e-3
Q_TURNBACK_FRAMES = 2
NEW_ARC_MIN_NODES = 4
MEANINGFUL_EXCHANGE = 1e-3
MEANINGFUL_CLOCK = 1e-3
PUBLISHED_OMEGA_NF512_T010 = 2998.824498749641
PUBLISHED_CAP_CROSSING_T = 0.0965
ONE_PERCENT = 0.01

# The 0.10 reversal and 0.50 share proxies are not requirements.
PROXY_REVERSAL_IS_REQUIRED = False
PROXY_SHARE_IS_REQUIRED = False

RUN_PLAN = (
    {
        "name": "nf256_dtmax_0_0005",
        "nf": 256,
        "handoff": "nf256_dt_0_0005",
        "dt_cap": DT_COARSE,
        "role": "coarse_pilot",
        "pool": "coarse",
    },
    {
        "name": "nf256_dtmax_0_00025",
        "nf": 256,
        "handoff": "nf256_dt_0_0005",
        "dt_cap": DT_CONFIRM,
        "role": "coarse_time",
        "pool": "coarse",
    },
    {
        "name": "nf512_dtmax_0_00025",
        "nf": 512,
        "handoff": "nf512_dt_0_0005",
        "dt_cap": DT_CONFIRM,
        "role": "fine_confirmation",
        "pool": "fine",
    },
    {
        "name": "nf512_dtmax_0_0005",
        "nf": 512,
        "handoff": "nf512_dt_0_0005",
        "dt_cap": DT_COARSE,
        "role": "fine_time",
        "pool": "fine",
    },
)

GEOMETRIC_EVENTS = (
    "BRIDGE_ARC_ENTERED_PACKET",
    "PACKET_ARC_ENTERED_BRIDGE",
    "NEW_SAME_SIGN_ARC",
    "Q_TURNBACK",
)

SERIES_SCALARS = (
    "energy",
    "field_energy",
    "gravity_energy",
    "full_hamilton_max",
    "projected_hamilton_max",
    "held_out_hamilton_max",
    "full_momentum_max",
    "projected_momentum_max",
    "held_out_momentum_max",
    "proper_min",
    "proper_max",
    "proper_mean",
    "chi_min",
    "chi_max",
    "c2_max",
    "r_min",
    "r_max",
    "Q_min",
    "Q_max",
    "K_r_min",
    "K_r_max",
    "K_perp_min",
    "K_perp_max",
    "lifted_fieldwork",
    "proper_pressure_power",
    "momentum_lapse_power",
    "coordinate_work_integral",
    "proper_pressure_work_integral",
    "unitarity",
    "number",
    "ward_max",
    "shell_gap_max",
    "kernel_gap_max",
    "source_gap_max",
    "partition_gap",
    "gram_gap",
    "G_min",
    "positive_G",
    "positive_r",
    "positive_Q",
    "positive_L",
    "leader",
    "leader_share",
    "leader_content",
    "packet_flux",
    "reservoir_flux",
    "proper_clock",
    "leader_clock",
    "q_min",
    "q_max",
    "q_at_areal_max",
    "omega",
    "dt",
)

WINDOW_CHANNELS = (
    "window_normal",
    "window_matter",
    "window_gravity",
    "window_coordinate",
    "window_proper_flux",
    "window_proper_work",
    "window_lapse_work",
    "window_quasilocal",
    "window_coordinate_work",
    "window_shift_transport",
    "window_shift_pressure",
    "window_geometric_shift",
    "window_geometric_proper",
    "observer_boundary",
    "observer_pressure",
    "observer_lapse",
)

SEALED_FILES = {
    "episode_json": episode.OUT,
    "episode_npz": episode.NPZ,
    "controls_v1_json": LAB / "results" / "development" / "nsc-regeneration-controls-v1.json",
    "controls_v1_npz": LAB / "results" / "development" / "nsc-regeneration-controls-v1.npz",
    "controls_v2_json": LAB / "results" / "development" / "nsc-regeneration-controls-v2.json",
    "null_expansion_json": LAB / "results" / "development" / "nsc-spherical-null-expansion-v1.json",
    "galerkin": galerkin._MODULE_PATH,
    "coupling": galerkin._COUPLING_PATH,
    "feedback_action": galerkin._ACTION_PATH,
    "regional": regional._MODULE_PATH,
    "controls_module": LAB / "src" / "recursive_horizons" / "nsc_regeneration_controls.py",
    "null_expansion_module": Path(null_expansion.__file__).resolve(),
    "conformal_source": galerkin._SOURCE_PATH,
    "manuscript_pdf": ROOT / "paper" / "nested-space-cosmology.pdf",
    "companion_pdf": ROOT / "paper" / "local-incoming-gate-draft.pdf",
}
SEALED_SOURCE_REF = "5f10ecd365843d1616e50eb16a20d7acd8377e2c"

def sha256_file(path):
    file_path = Path(path)
    if not file_path.is_file():
        return None
    digest = hashlib.sha256()
    with file_path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload):
    return hashlib.sha256(payload).hexdigest()


def sealed_hashes():
    return {name: sha256_file(path) for name, path in SEALED_FILES.items()}


def snap_time(value):
    return float(np.round(float(value), 10))


def jsonable(value):
    return episode.jsonable(value)


def assert_owned(path):
    resolved = Path(path).resolve()
    if resolved not in OWNED_OUTPUTS:
        raise RuntimeError(f"refusing to write outside the owned outputs: {resolved}")
    return resolved


def write_json(path, record):
    path = assert_owned(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(jsonable(record), indent=2, allow_nan=False) + "\n")
    os.replace(temporary, path)


def write_npz(path, arrays):
    path = assert_owned(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    clean = {}
    for key, value in arrays.items():
        if isinstance(value, np.ndarray):
            clean[key] = value
        else:
            raise TypeError(f"npz value {key} is not an array")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **clean)
    os.replace(temporary, path)


def write_text(path, text):
    path = assert_owned(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text)
    os.replace(temporary, path)


class Pool:
    def __init__(self, limit):
        self.limit = float(limit)
        self.used = 0.0

    def add(self, seconds):
        self.used += float(seconds)

    def allows(self, step_cost):
        return episode.budget_allows(self.used, step_cost, self.limit)


def build_grid(fermions):
    return galerkin.build_grid(int(fermions))


def load_handoff(spec):
    """Bitwise Cauchy state from the saved T=0.05 final. No new source solve."""
    state = load_episode_final(episode.NPZ, spec["handoff"])
    grid = build_grid(spec["nf"])
    if state.r.shape != (grid.ng,) or state.phi0.shape != (grid.nf, 6):
        raise ValueError(f"{spec['name']} handoff does not match the owned grid")
    if not np.array_equal(grid.fine.occupations, galerkin.OCCUPATIONS):
        raise ValueError("grid occupations are not the original six weights")
    return grid, state


def settings_digest(grid):
    coefficients = locked_coefficients()
    payload = {
        "occupations": [float(value) for value in grid.fine.occupations],
        "calibration": {
            key: float(grid.fine.calibration[key])
            for key in sorted(grid.fine.calibration)
            if np.ndim(grid.fine.calibration[key]) == 0
        },
        "coefficients": {
            key: float(coefficients[key])
            for key in sorted(coefficients)
            if np.ndim(coefficients[key]) == 0
        },
        "kappa": float(grid.fine.kappa),
        "length": float(grid.length),
        "nf": int(grid.nf),
        "nq": int(grid.nq),
    }
    encoded = json.dumps(payload, sort_keys=True).encode("utf-8")
    return sha256_bytes(encoded), payload


def choose_dt(time_value, target, omega, dt_cap):
    """Next owned step. The safety cap, the run cap, and the frame can bind."""
    omega = float(max(omega, 1.0))
    safety = SAFETY_CAP / omega
    room_mark = snap_time((math.floor((time_value + 1e-12) / FRAME) + 1) * FRAME)
    if room_mark <= time_value:
        room_mark = snap_time(time_value + FRAME)
    mark = snap_time(min(target, room_mark))
    room = snap_time(mark - time_value)
    if room <= 0.0:
        return 0.0, mark, "on_mark"
    options = {
        "safety_cap": safety,
        "owned_dt_cap": OWNED_DT_CAP,
        "run_cap": float(dt_cap),
        "frame": room,
    }
    dt = min(options.values())
    binding = [name for name, value in options.items() if abs(value - dt) <= 1e-15 * max(1.0, abs(dt))]
    # A tie is recorded. The unique strict minimum is the limiter.
    strictly = [
        name for name, value in options.items()
        if value <= dt + 1e-15 and all(
            value < other - 1e-15 * max(1.0, abs(other))
            for other_name, other in options.items() if other_name != name
        )
    ]
    if len(strictly) == 1:
        limiter = strictly[0]
    else:
        limiter = "tie:" + "+".join(binding)
    return float(dt), mark, limiter


def prolong_rate(grid, coarse):
    return CauchyRate(
        galerkin.prolong_geometry(grid, coarse.Q),
        galerkin.prolong_geometry(grid, coarse.r),
        galerkin.prolong_geometry(grid, coarse.chi),
        galerkin.prolong_geometry(grid, coarse.p_Q),
        galerkin.prolong_geometry(grid, coarse.p_r),
        galerkin.prolong_geometry(grid, coarse.p_chi),
        galerkin.prolong_columns(grid, coarse.phi0),
        galerkin.prolong_columns(grid, coarse.phi1),
        np.zeros(grid.nq),
        np.zeros(grid.nq),
        np.zeros(grid.nq),
        float(coarse.fieldwork_power),
    )


def observer_window_channels(grid, state, coarse, windows):
    """Normal-window balance from the owned observer flux, pressure, and lapse."""
    fine = galerkin.prolong_state(grid, state)
    lifted = prolong_rate(grid, coarse)
    terms = regional.proper_balance_terms(grid.fine, fine, lifted)
    derivative = grid.fine.derivative
    boundary = []
    pressure = []
    lapse = []
    for window in windows:
        slope = derivative @ window
        boundary.append(float(np.sum(slope * terms["flux_nodal"])))
        pressure.append(float(np.sum(window * terms["proper_pressure_work"])))
        lapse.append(float(np.sum(window * terms["momentum_lapse_work"])))
    return {
        "observer_boundary": np.asarray(boundary, dtype=float),
        "observer_pressure": np.asarray(pressure, dtype=float),
        "observer_lapse": np.asarray(lapse, dtype=float),
        "predicted_normal_rate": np.asarray(boundary, dtype=float)
        + np.asarray(pressure, dtype=float)
        + np.asarray(lapse, dtype=float),
    }


def trap_integrate(times, values):
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    if times.size < 2:
        shape = values.shape[1:] if values.ndim > 1 else ()
        return np.zeros(shape, dtype=float)
    dt = np.diff(times)
    midpoint = 0.5 * (values[:-1] + values[1:])
    scale = dt.reshape((-1,) + (1,) * (values.ndim - 1))
    return np.sum(midpoint * scale, axis=0)


def interval_overlaps(interval, low, high):
    """Closed node interval against a half-open coordinate interval on the circle."""
    start = float(interval["x_start"])
    end = float(interval["x_end"])
    if interval.get("wraps"):
        return interval_overlaps(
            {"x_start": start, "x_end": PERIOD, "wraps": False}, low, high
        ) or interval_overlaps(
            {"x_start": 0.0, "x_end": end, "wraps": False}, low, high
        )
    return bool(start < float(high) and end >= float(low))


def expansion_margin():
    path = SEALED_FILES["null_expansion_json"]
    record = json.loads(path.read_text())
    margin = record.get("margin") or {}
    value = float(margin["absolute_expansion_margin"])
    return value, {
        "source": str(path.relative_to(ROOT)),
        "absolute_expansion_margin": value,
        "rule": margin.get("rule"),
        "retuned": False,
    }


def classify_expansion(time_value, fields, margin):
    radius = np.asarray(fields["r"], dtype=float)
    coordinate, length_density, shift = null_expansion.static_clock(radius.size, PERIOD)
    radius_x = null_expansion.spectral_derivative(radius, PERIOD)
    frame = null_expansion.summarize_frame(
        time_value,
        radius,
        np.asarray(fields["Q"], dtype=float),
        np.asarray(fields["proper"], dtype=float),
        np.asarray(fields["K_perp"], dtype=float),
        radius_x,
        length_density,
        shift,
        margin,
    )
    public = {key: value for key, value in frame.items() if not str(key).startswith("_")}
    return public, frame


def arc_flags(frame):
    trapped = frame.get("trapped_intervals") or []
    anti = frame.get("anti_trapped_intervals") or []
    trapped_in_packet = any(interval_overlaps(item, 0.0, PACKET_END) for item in trapped)
    trapped_in_bridge = any(interval_overlaps(item, BRIDGE_START, PERIOD) for item in trapped)
    anti_in_packet = any(interval_overlaps(item, 0.0, PACKET_END) for item in anti)
    anti_in_bridge = any(interval_overlaps(item, BRIDGE_START, PERIOD) for item in anti)
    large = [
        item for item in list(trapped) + list(anti)
        if int(item.get("count", 0)) >= NEW_ARC_MIN_NODES
    ]
    return {
        "trapped_in_packet": bool(trapped_in_packet),
        "trapped_in_bridge": bool(trapped_in_bridge),
        "anti_in_packet": bool(anti_in_packet),
        "anti_in_bridge": bool(anti_in_bridge),
        "large_same_sign_arcs": int(len(large)),
        "trapped_count": int(len(trapped)),
        "anti_count": int(len(anti)),
    }


def geometric_event(history):
    """Declared events. The handoff frame is the baseline and cannot fire."""
    if len(history) < 2:
        return None
    start = history[0]["flags"]
    current = history[-1]["flags"]
    if current["anti_in_packet"] and not start["anti_in_packet"]:
        return "BRIDGE_ARC_ENTERED_PACKET"
    if current["trapped_in_bridge"] and not start["trapped_in_bridge"]:
        return "PACKET_ARC_ENTERED_BRIDGE"
    if current["large_same_sign_arcs"] >= 3 and start["large_same_sign_arcs"] < 3:
        return "NEW_SAME_SIGN_ARC"
    q_values = np.asarray([row["q_at_areal_max"] for row in history], dtype=float)
    if q_values.size >= Q_TURNBACK_FRAMES + 1:
        from recursive_horizons.nsc_regeneration_controls import reversal_amplitude

        if float(reversal_amplitude(q_values)) >= Q_TURNBACK_FLOOR:
            diffs = np.diff(q_values)
            trend = float(q_values[-1] - q_values[0])
            adverse = diffs < 0.0 if trend >= 0.0 else diffs > 0.0
            if np.all(adverse[-Q_TURNBACK_FRAMES:]):
                return "Q_TURNBACK"
    return None


def endpoint_closure(grid, state, coarse, windows):
    """Owned instantaneous window closure at one state. Not an interval integral."""
    try:
        return _endpoint_closure(grid, state, coarse, windows)
    except (PositiveChartExit, ValueError, FloatingPointError) as error:
        return {"computed": False, "error": str(error)}


def _endpoint_closure(grid, state, coarse, windows):
    fine = galerkin.prolong_state(grid, state)
    lifted = prolong_rate(grid, coarse)
    matter_slope = regional.matter_pointwise_slope(grid.fine, fine, lifted)
    total_slope = regional.total_pointwise_slope(grid.fine, fine, lifted)
    rows = []
    for window, name in zip(windows, episode.WINDOW_NAMES, strict=True):
        prediction = regional.total_window_prediction(grid.fine, fine, lifted, window)
        matter = prediction["matter_channels"]
        matter_prediction = float(sum(matter[key] for key in (
            "shift_transport", "proper_normal_flux", "shift_pressure_cross", "coordinate_metric_work",
        )))
        rows.append({
            "name": name,
            "constraint_smearing": float(prediction["constraint_smearing"]),
            "matter_closure_error": float(np.sum(window * matter_slope) - matter_prediction),
            "total_closure_error": float(np.sum(window * total_slope) - float(prediction["prediction"])),
            "product_gap_max": float(prediction["product_gap_max"]),
            "ward_max": float(prediction["ward_max"]),
        })
    return rows


def one_step_probe(grid, state, dt_cap):
    """Discarded timing and ledger-identity step. The caller's state is kept."""
    windows, partition = episode.windows_for(grid)
    held = state.copy()
    held_hash = state_sha256(held)
    started = time.process_time()
    sample0, coarse0, fields0 = episode.observe(grid, held, windows, partition)
    channels0 = observer_window_channels(grid, held, coarse0, windows)
    step0, omega0, quadrature_omega = galerkin.stable_timestep(grid, held)
    dt, _mark, limiter = choose_dt(T0, TARGET, omega0, dt_cap)
    if dt <= 0.0 or dt * omega0 > ABSOLUTE_RK4:
        raise RuntimeError("probe could not choose an admissible step")
    candidate = galerkin.rk4_step(grid, held, dt)
    failure = chart_failure(grid.fine, galerkin.prolong_state(grid, candidate))
    sample1, coarse1, _fields1 = episode.observe(grid, candidate, windows, partition)
    channels1 = observer_window_channels(grid, candidate, coarse1, windows)
    elapsed = time.process_time() - started
    if state_sha256(held) != held_hash:
        raise RuntimeError("probe mutated the handoff state")
    actual = np.asarray(sample1["window_normal"], dtype=float) - np.asarray(sample0["window_normal"], dtype=float)
    predicted = 0.5 * dt * (channels0["predicted_normal_rate"] + channels1["predicted_normal_rate"])
    gap = actual - predicted
    scale = max(float(np.max(np.abs(actual))), float(np.max(np.abs(predicted))), 1e-8)
    opposite = actual + predicted
    return {
        "cpu_seconds": float(elapsed),
        "dt": float(dt),
        "owned_stable_timestep": float(step0),
        "omega": float(omega0),
        "quadrature_omega": float(quadrature_omega),
        "limiter": limiter,
        "dt_omega": float(dt * omega0),
        "chart_failure": failure,
        "gram_gap": float(car_report(candidate)["gap"]),
        "identity_gap_max": float(np.max(np.abs(gap))),
        "identity_scale": float(scale),
        "identity_relative": float(np.max(np.abs(gap)) / scale),
        "opposite_sign_gap_max": float(np.max(np.abs(opposite))),
        "sign_error": bool(
            np.max(np.abs(gap)) > 0.5 * scale and np.max(np.abs(opposite)) < 0.5 * np.max(np.abs(gap))
        ),
        "pressure_matches_observe": float(np.max(np.abs(
            channels0["observer_pressure"] - np.asarray(sample0["window_proper_work"], dtype=float)
        ))),
        "lapse_matches_observe": float(np.max(np.abs(
            channels0["observer_lapse"] - np.asarray(sample0["window_lapse_work"], dtype=float)
        ))),
        "Q_min": float(sample0["Q_min"]),
        "chi_max": float(sample0["chi_max"]),
        "state_unchanged": True,
        "kept": False,
    }


def forecast_steps(omega0, dt_cap, duration, growth):
    """Integrate 1/dt along a declared omega model. This is a cost forecast."""
    steps = 0
    time_value = 0.0
    dt_min = None
    while time_value < duration - 1e-12:
        omega = float(omega0 * math.exp(growth * time_value))
        dt, _mark, _limiter = choose_dt(T0 + time_value, T0 + duration, omega, dt_cap)
        if dt <= 0.0:
            break
        steps += 1
        time_value = snap_time(time_value + dt)
        dt_min = dt if dt_min is None else min(dt_min, dt)
        if steps > 2_000_000:
            break
    return {"steps": int(steps), "dt_min": dt_min, "omega_end": float(omega0 * math.exp(growth * duration))}


def forecast_models(omega0, dt_cap, step_cost, duration):
    logarithmic = 0.0
    if omega0 > 0.0:
        logarithmic = math.log(PUBLISHED_OMEGA_NF512_T010 / omega0) / (0.10 - T0)
    late_slope = (PUBLISHED_OMEGA_NF512_T010 - (SAFETY_CAP / DT_COARSE)) / (0.10 - PUBLISHED_CAP_CROSSING_T)
    # The published late slope is converted to an exponential envelope that
    # reaches the same end omega. It is a cost bound, not a fitted law.
    end_omega = max(omega0, omega0 + late_slope * duration)
    late_growth = 0.0 if omega0 <= 0.0 else math.log(end_omega / omega0) / duration
    rows = {}
    for name, growth in (("published_mean_growth", max(0.0, logarithmic)), ("published_late_growth", max(0.0, late_growth))):
        estimate = forecast_steps(omega0, dt_cap, duration, growth)
        estimate["growth"] = float(growth)
        estimate["cpu_seconds"] = float(estimate["steps"] * step_cost * episode.BUDGET_SAFETY)
        rows[name] = estimate
    rows["constant_initial_dt"] = {
        "steps": int(math.ceil(duration / min(dt_cap, OWNED_DT_CAP) - 1e-12)),
        "cpu_seconds": float(math.ceil(duration / min(dt_cap, OWNED_DT_CAP) - 1e-12) * step_cost * episode.BUDGET_SAFETY),
        "growth": 0.0,
    }
    return rows


def empty_series():
    series = {key: [] for key in SERIES_SCALARS}
    series["time"] = []
    for key in WINDOW_CHANNELS:
        series[key] = []
    return series


def append_sample(series, sample, extra):
    series["time"].append(float(extra["time"]))
    for key in SERIES_SCALARS:
        if key in extra:
            series[key].append(float(extra[key]))
        else:
            series[key].append(float(sample[key]))
    for key in WINDOW_CHANNELS:
        if key in extra:
            series[key].append(np.asarray(extra[key], dtype=float))
        else:
            series[key].append(np.asarray(sample[key], dtype=float))


def stack_series(series):
    stacked = {"time": np.asarray(series["time"], dtype=float)}
    for key in SERIES_SCALARS:
        stacked[key] = np.asarray(series[key], dtype=float)
    for key in WINDOW_CHANNELS:
        stacked[key] = np.vstack(series[key]) if series[key] else np.zeros((0, 4))
    return stacked


def evolve_case(spec, grid, state, target, pool, margin, label):
    """Adaptive RK4 from the handoff. Failed steps are not kept."""
    windows, partition = episode.windows_for(grid)
    current = state.copy()
    initial_hash = state_sha256(current)
    initial_phi = np.array(current.phi0, copy=True)
    occupation_bytes = np.ascontiguousarray(grid.fine.occupations).tobytes()
    series = empty_series()
    frames = []
    expansion_history = []
    limiter_counts = {}
    dt_values = []
    omega_values = []
    dt_omega_values = []
    safety_after_step = 0
    rejected = []
    coordinate_work = 0.0
    proper_work = 0.0
    proper_clock = 0.0
    leader_clock = 0.0
    previous_power = None
    previous_proper = None
    previous_global = None
    previous_leader = None
    previous_rate = None
    previous_q_coarse = None
    stop_reason = None
    chart_stop = False
    completed = 0
    step_cost = 0.0
    cpu_start = time.process_time()
    endpoint = {"initial": None, "final": None}
    accepted_coarse = None
    charged = 0.0
    error_text = None
    target = snap_time(target)

    def store_frame(sample, fields, coarse, absolute, increment, increment_dt, direction, q_node):
        frames.append({
            "time": float(absolute),
            "fields": fields,
            "coarse_Q": np.array(current.Q, copy=True),
            "coarse_r": np.array(current.r, copy=True),
            "coarse_chi": np.array(current.chi, copy=True),
            "coarse_p_Q": np.array(current.p_Q, copy=True),
            "coarse_p_r": np.array(current.p_r, copy=True),
            "coarse_p_chi": np.array(current.p_chi, copy=True),
            "Q_dot": np.array(coarse.Q, copy=True),
            "p_chi_dot": np.array(coarse.p_chi, copy=True),
            "phi0": np.array(current.phi0, copy=True),
            "phi1": np.array(current.phi1, copy=True),
            "increment_Q_dot": None if increment is None else np.array(increment[0], copy=True),
            "increment_p_chi": None if increment is None else np.array(increment[1], copy=True),
            "increment_Q": None if increment is None else np.array(increment[2], copy=True),
            "increment_dt": None if increment_dt is None else float(increment_dt),
            "increment_direction": direction,
            "q_at_areal_max": float(q_node),
        })

    def remember(sample, fields, coarse, channels, absolute, frame_increment=None, frame_dt=None, direction=None):
        nonlocal proper_clock, leader_clock
        view = shell_views(sample)
        chart = chart_from_fields(grid, fields)
        gram = car_report(current)
        global_rate, leader_rate = proper_clock_rates(grid, fields, windows, view["leader"])
        q_nodal = np.asarray(fields["r"], dtype=float) * np.asarray(fields["Q"], dtype=float)
        q_star = float(q_nodal[int(np.argmax(fields["r"]))])
        extra = {
            "time": float(absolute),
            "gram_gap": float(gram["gap"]),
            "G_min": float(chart["G_min"]),
            "positive_G": float(bool(chart["positive_G"])),
            "positive_r": float(bool(chart["positive_r"])),
            "positive_Q": float(bool(chart["positive_Q"])),
            "positive_L": float(bool(chart["positive_L"])),
            "leader": float(view["leader"]),
            "leader_share": float(view["leader_share"]),
            "leader_content": float(view["content"][view["leader"]]),
            "packet_flux": float(view["packet_flux"]),
            "reservoir_flux": float(view["reservoir_flux"]),
            "proper_clock": float(proper_clock),
            "leader_clock": float(leader_clock),
            "q_min": float(np.min(q_nodal)),
            "q_max": float(np.max(q_nodal)),
            "q_at_areal_max": q_star,
            "omega": float(channels["omega"]),
            "dt": float(channels["dt"]),
            "observer_boundary": channels["observer_boundary"],
            "observer_pressure": channels["observer_pressure"],
            "observer_lapse": channels["observer_lapse"],
        }
        append_sample(series, sample, extra)
        return chart, view, gram, global_rate, leader_rate, q_star

    try:
        sample, coarse, fields = episode.observe(grid, current, windows, partition)
        channels = observer_window_channels(grid, current, coarse, windows)
        _step, omega, _quadrature = galerkin.stable_timestep(grid, current)
        channels["omega"] = float(omega)
        channels["dt"] = 0.0
        sample["time"] = T0
        sample["coordinate_work_integral"] = 0.0
        sample["proper_pressure_work_integral"] = 0.0
        chart, view, gram, global_rate, leader_rate, q_star = remember(
            sample, fields, coarse, channels, T0
        )
        previous_power = float(sample["lifted_fieldwork"])
        previous_proper = float(sample["proper_pressure_power"])
        previous_global = float(global_rate)
        previous_leader = float(leader_rate)
        if not gram["admissible"] or chart_failure(grid.fine, galerkin.prolong_state(grid, current)):
            stop_reason = "HANDOFF_INADMISSIBLE"
            chart_stop = True
        else:
            public, _private = classify_expansion(T0, fields, margin)
            flags = arc_flags(public)
            expansion_history.append({"time": T0, "flags": flags, "frame": public, "q_at_areal_max": q_star})
            store_frame(sample, fields, coarse, T0, None, None, None, q_star)
            accepted_coarse = coarse
            endpoint["initial"] = endpoint_closure(grid, current, coarse, windows)
            startup = time.process_time() - cpu_start
            pool.add(startup)
            charged += startup
        while stop_reason is None and snap_time(series["time"][-1]) < target - 1e-12:
            if not pool.allows(step_cost if completed else 0.0):
                stop_reason = "CPU_BUDGET"
                break
            dt, mark, limiter = choose_dt(snap_time(series["time"][-1]), target, omega, spec["dt_cap"])
            if dt <= 0.0:
                stop_reason = "NO_ADMISSIBLE_DT"
                break
            if dt * omega > ABSOLUTE_RK4:
                stop_reason = "ABSOLUTE_RK4_DOMAIN"
                rejected.append({"reason": stop_reason, "time": float(series["time"][-1]), "kept": False})
                break
            started = time.process_time()
            try:
                candidate = galerkin.rk4_step(grid, current, dt)
            except PositiveChartExit as exit_chart:
                stop_reason = str(exit_chart.reason)
                chart_stop = True
                rejected.append({"reason": stop_reason, "time": float(series["time"][-1]), "kept": False})
                break
            prolonged = galerkin.prolong_state(grid, candidate)
            failure = chart_failure(grid.fine, prolonged)
            if failure:
                stop_reason = str(failure)
                chart_stop = True
                rejected.append({"reason": stop_reason, "time": float(series["time"][-1]), "kept": False})
                break
            try:
                sample, coarse, fields = episode.observe(grid, candidate, windows, partition)
                channels = observer_window_channels(grid, candidate, coarse, windows)
            except PositiveChartExit as exit_chart:
                stop_reason = str(exit_chart.reason)
                chart_stop = True
                rejected.append({"reason": stop_reason, "time": float(series["time"][-1]), "kept": False})
                break
            gram = car_report(candidate)
            if float(sample["unitarity"]) > CAR_GAP_MAX:
                gram = dict(gram)
                gram["admissible"] = False
            if not gram["admissible"]:
                stop_reason = "GRAM_LEFT_IDENTITY"
                rejected.append({"reason": stop_reason, "time": float(series["time"][-1]), "kept": False})
                break
            _stable_after, omega_after, _quad_after = galerkin.stable_timestep(grid, candidate)
            if dt * float(omega_after) > ABSOLUTE_RK4:
                stop_reason = "ABSOLUTE_RK4_DOMAIN"
                rejected.append({
                    "reason": stop_reason,
                    "time": float(series["time"][-1]),
                    "dt_omega_after": float(dt * omega_after),
                    "kept": False,
                })
                break
            if dt * float(omega_after) > SAFETY_CAP:
                safety_after_step += 1
            # Accepted. The candidate replaces the state only here.
            if previous_q_coarse is None:
                increment = (
                    np.asarray(coarse.Q) - frames[0]["Q_dot"],
                    np.asarray(coarse.p_chi) - frames[0]["p_chi_dot"],
                    np.asarray(candidate.Q) - np.asarray(current.Q),
                )
            else:
                increment = (
                    np.asarray(coarse.Q) - previous_rate[0],
                    np.asarray(coarse.p_chi) - previous_rate[1],
                    np.asarray(candidate.Q) - previous_q_coarse,
                )
            current = candidate
            accepted_coarse = coarse
            if np.ascontiguousarray(grid.fine.occupations).tobytes() != occupation_bytes:
                stop_reason = "OCCUPATIONS_CHANGED"
                break
            absolute = snap_time(series["time"][-1] + dt)
            if abs(absolute - mark) <= 5e-10:
                absolute = snap_time(mark)
            coordinate_work += 0.5 * dt * (previous_power + float(sample["lifted_fieldwork"]))
            proper_work += 0.5 * dt * (previous_proper + float(sample["proper_pressure_power"]))
            sample["coordinate_work_integral"] = float(coordinate_work)
            sample["proper_pressure_work_integral"] = float(proper_work)
            previous_power = float(sample["lifted_fieldwork"])
            previous_proper = float(sample["proper_pressure_power"])
            view = shell_views(sample)
            global_rate, leader_rate = proper_clock_rates(grid, fields, windows, view["leader"])
            proper_clock += 0.5 * dt * (previous_global + global_rate)
            leader_clock += 0.5 * dt * (previous_leader + leader_rate)
            previous_global = global_rate
            previous_leader = leader_rate
            omega = float(omega_after)
            channels["omega"] = omega
            channels["dt"] = float(dt)
            remember(sample, fields, coarse, channels, absolute)
            completed += 1
            limiter_counts[limiter] = limiter_counts.get(limiter, 0) + 1
            dt_values.append(float(dt))
            omega_values.append(omega)
            dt_omega_values.append(float(dt * omega))
            step_cost = time.process_time() - started
            pool.add(step_cost)
            charged += step_cost
            on_frame = abs(absolute / FRAME - round(absolute / FRAME)) <= 1e-8
            if completed == 1 and frames:
                frames[0]["increment_Q_dot"] = np.array(increment[0], copy=True)
                frames[0]["increment_p_chi"] = np.array(increment[1], copy=True)
                frames[0]["increment_Q"] = np.array(increment[2], copy=True)
                frames[0]["increment_dt"] = float(dt)
                frames[0]["increment_direction"] = "forward"
            if on_frame and absolute > T0 + 1e-12:
                q_star = float(series["q_at_areal_max"][-1])
                store_frame(
                    sample, fields, coarse, absolute, increment, dt, "backward", q_star,
                )
                public, _private = classify_expansion(absolute, fields, margin)
                flags = arc_flags(public)
                expansion_history.append({
                    "time": float(absolute),
                    "flags": flags,
                    "frame": public,
                    "q_at_areal_max": q_star,
                })
                event = geometric_event(expansion_history)
                if event:
                    stop_reason = event
            previous_rate = (np.array(coarse.Q, copy=True), np.array(coarse.p_chi, copy=True))
            previous_q_coarse = np.array(current.Q, copy=True)
            if completed == 1 or completed % 20 == 0 or stop_reason or absolute >= target - 1e-12:
                print(
                    label,
                    "step", completed,
                    "T", absolute,
                    "dt", dt,
                    "cpu_pool", round(pool.used, 3),
                    "Qmin", round(float(sample["Q_min"]), 6),
                    "chi", round(float(sample["chi_max"]), 4),
                    "share", round(float(view["leader_share"]), 4),
                    flush=True,
                )
        if stop_reason is None and series["time"] and snap_time(series["time"][-1]) + 1e-12 >= target:
            stop_reason = None
        if series["time"] and accepted_coarse is not None:
            endpoint["final"] = endpoint_closure(grid, current, accepted_coarse, windows)
            tail = time.process_time() - cpu_start - charged
            if tail > 0.0:
                pool.add(tail)
    except Exception:
        error_text = traceback.format_exc()
        stop_reason = stop_reason or "CODE_EXCEPTION"
        print(error_text, flush=True)
    attained = float(series["time"][-1]) if series["time"] else T0
    phi_gap = float(np.max(np.abs(np.asarray(current.phi0) - initial_phi)))
    return {
        "name": spec["name"],
        "role": spec["role"],
        "pool": spec["pool"],
        "nf": int(spec["nf"]),
        "dt_cap": float(spec["dt_cap"]),
        "handoff": spec["handoff"],
        "series": stack_series(series),
        "frames": frames,
        "expansion": expansion_history,
        "state": current,
        "initial_hash": initial_hash,
        "final_hash": state_sha256(current),
        "phi_max_abs_change": phi_gap,
        "source_reset": bool(completed > 0 and phi_gap == 0.0),
        "attained_time": attained,
        "target_time": float(target),
        "completed_to_target": bool(attained + 1e-12 >= target and stop_reason not in GEOMETRIC_EVENTS and not chart_stop),
        "stop_reason": stop_reason,
        "chart_stop": bool(chart_stop),
        "steps_completed": int(completed),
        "rejected": rejected,
        "rejected_count": int(len(rejected)),
        "limiter_counts": limiter_counts,
        "dt_values": np.asarray(dt_values, dtype=float),
        "omega_values": np.asarray(omega_values, dtype=float),
        "dt_omega_values": np.asarray(dt_omega_values, dtype=float),
        "safety_cap_after_step_count": int(safety_after_step),
        "cpu_seconds": float(time.process_time() - cpu_start),
        "endpoint_closure": endpoint,
        "error": error_text,
        "occupations_unchanged": bool(
            np.ascontiguousarray(grid.fine.occupations).tobytes() == occupation_bytes
        ),
    }


def sibling_gap(primary_name, other_name):
    left = load_episode_final(episode.NPZ, primary_name)
    right = load_episode_final(episode.NPZ, other_name)
    gaps = {}
    for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1"):
        gaps[name] = float(np.max(np.abs(np.asarray(getattr(left, name)) - np.asarray(getattr(right, name)))))
    return gaps


def link_report(grid, state):
    report = physical_link_report(grid, state)
    return {
        "max_difference_from_old_B": report["max_difference_from_old_B"],
        "B_phys_equals_old_B": report["B_phys_equals_old_B"],
        "old_B_not_used": report["old_B_not_used"],
        "copied_from_old_B": report["copied_from_old_B"],
        "old_circulation_transferred": False,
    }


def pack_frames(prefix, frames):
    arrays = {}
    if not frames:
        return arrays
    arrays[prefix + "frame_time"] = np.asarray([frame["time"] for frame in frames], dtype=float)
    for key, name in (
        ("coarse_Q", "frame_coarse_Q"),
        ("coarse_r", "frame_coarse_r"),
        ("coarse_chi", "frame_coarse_chi"),
        ("coarse_p_Q", "frame_coarse_p_Q"),
        ("coarse_p_r", "frame_coarse_p_r"),
        ("coarse_p_chi", "frame_coarse_p_chi"),
        ("Q_dot", "frame_Q_dot"),
        ("p_chi_dot", "frame_p_chi_dot"),
        ("phi0", "frame_phi0"),
        ("phi1", "frame_phi1"),
    ):
        arrays[prefix + name] = np.stack([frame[key] for frame in frames])
    for key, name in (
        ("r", "frame_quad_r"),
        ("Q", "frame_quad_Q"),
        ("chi", "frame_quad_chi"),
        ("proper", "frame_quad_proper"),
        ("K_r", "frame_quad_K_r"),
        ("K_perp", "frame_quad_K_perp"),
        ("rho", "frame_quad_rho"),
        ("current", "frame_quad_current"),
        ("weyl_C2", "frame_quad_weyl_C2"),
    ):
        arrays[prefix + name] = np.stack([frame["fields"][key] for frame in frames])
    increment_q = []
    increment_p = []
    increment_state = []
    dts = []
    directions = []
    present = []
    for frame in frames:
        present.append(frame["increment_Q_dot"] is not None)
        increment_q.append(np.zeros_like(frame["Q_dot"]) if frame["increment_Q_dot"] is None else frame["increment_Q_dot"])
        increment_p.append(np.zeros_like(frame["p_chi_dot"]) if frame["increment_p_chi"] is None else frame["increment_p_chi"])
        increment_state.append(np.zeros_like(frame["coarse_Q"]) if frame["increment_Q"] is None else frame["increment_Q"])
        dts.append(np.nan if frame["increment_dt"] is None else float(frame["increment_dt"]))
        directions.append(frame["increment_direction"] or "")
    arrays[prefix + "frame_increment_Q_dot"] = np.stack(increment_q)
    arrays[prefix + "frame_increment_p_chi"] = np.stack(increment_p)
    arrays[prefix + "frame_increment_Q"] = np.stack(increment_state)
    arrays[prefix + "frame_increment_dt"] = np.asarray(dts, dtype=float)
    arrays[prefix + "frame_increment_present"] = np.asarray(present, dtype=bool)
    arrays[prefix + "frame_increment_direction"] = np.asarray(directions)
    indicator_q = np.full_like(arrays[prefix + "frame_increment_Q_dot"], np.nan, dtype=float)
    indicator_p = np.full_like(arrays[prefix + "frame_increment_p_chi"], np.nan, dtype=float)
    for index, frame in enumerate(frames):
        if frame["increment_dt"]:
            indicator_q[index] = increment_q[index] / float(frame["increment_dt"])
            indicator_p[index] = increment_p[index] / float(frame["increment_dt"])
    arrays[prefix + "frame_indicator_Q_dot"] = indicator_q
    arrays[prefix + "frame_indicator_p_chi"] = indicator_p
    return arrays


def pack_run(result):
    prefix = result["name"] + "_"
    arrays = pack_frames(prefix, result["frames"])
    series = result["series"]
    for key, value in series.items():
        arrays[prefix + key] = np.asarray(value)
    state = result["state"]
    for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1"):
        arrays[prefix + "final_" + name] = np.asarray(getattr(state, name))
    arrays[prefix + "dt_values"] = np.asarray(result["dt_values"], dtype=float)
    arrays[prefix + "omega_values"] = np.asarray(result["omega_values"], dtype=float)
    arrays[prefix + "dt_omega_values"] = np.asarray(result["dt_omega_values"], dtype=float)
    return arrays


def balance_block(series):
    if series["time"].size < 2:
        return None
    energy_change = float(series["energy"][-1] - series["energy"][0])
    field_change = float(series["field_energy"][-1] - series["field_energy"][0])
    gravity_change = float(series["gravity_energy"][-1] - series["gravity_energy"][0])
    work = float(series["coordinate_work_integral"][-1])
    pressure = float(series["proper_pressure_work_integral"][-1])
    lapse = float(trap_integrate(series["time"], series["momentum_lapse_power"]))
    normal_change = series["window_normal"][-1] - series["window_normal"][0]
    boundary = trap_integrate(series["time"], series["observer_boundary"])
    pressure_w = trap_integrate(series["time"], series["observer_pressure"])
    lapse_w = trap_integrate(series["time"], series["observer_lapse"])
    predicted = boundary + pressure_w + lapse_w
    residual = normal_change - predicted
    flux_sum = series["window_proper_flux"].sum(axis=1)
    measured_change = series["window_coordinate"][-1] - series["window_coordinate"][0]
    noether = trap_integrate(series["time"], series["window_proper_flux"])
    quasilocal = trap_integrate(series["time"], series["window_quasilocal"])
    coordinate_window = trap_integrate(series["time"], series["window_coordinate_work"])
    exchange_scale = np.maximum.reduce([
        np.abs(normal_change),
        np.abs(boundary),
        np.abs(pressure_w),
        np.abs(lapse_w),
        np.full(4, 1e-8),
    ])
    window_ratios = np.abs(residual) / exchange_scale
    return {
        "field_energy_change": field_change,
        "gravity_energy_change": gravity_change,
        "total_energy_change": energy_change,
        "field_plus_gravity_change": field_change + gravity_change,
        "coordinate_work_integral": work,
        "proper_pressure_work_integral": pressure,
        "lapse_gradient_work_integral": float(lapse),
        "balance_versus_field_exchange": episode.balance_against_exchange(energy_change, field_change),
        "balance_versus_coordinate_work": episode.balance_against_exchange(energy_change, work),
        "hidden_pump_term_added": False,
        "flux_sum_max_abs": float(np.max(np.abs(flux_sum))),
        "flux_sum_is_window_balance": False,
        "flux_sum_is_not_the_window_balance": True,
        "window_normal_change": normal_change.tolist(),
        "window_observer_boundary_integral": boundary.tolist(),
        "window_pressure_work_integral": pressure_w.tolist(),
        "window_lapse_exchange_integral": lapse_w.tolist(),
        "window_normal_residual": residual.tolist(),
        "window_normal_residual_over_exchange": window_ratios.tolist(),
        "window_normal_within_one_percent": bool(np.all(window_ratios <= ONE_PERCENT)),
        "window_measured_energy_change": measured_change.tolist(),
        "window_noether_proper_flux_integral": noether.tolist(),
        "window_quasilocal_integral": quasilocal.tolist(),
        "window_coordinate_work_integral": coordinate_window.tolist(),
        "spatial_windows_are_not_mode_projectors": True,
        "mode_projectors": {
            "full_hamilton_max": float(np.max(series["full_hamilton_max"])),
            "projected_hamilton_max": float(np.max(series["projected_hamilton_max"])),
            "held_out_hamilton_max": float(np.max(series["held_out_hamilton_max"])),
            "full_momentum_max": float(np.max(series["full_momentum_max"])),
            "projected_momentum_max": float(np.max(series["projected_momentum_max"])),
            "held_out_momentum_max": float(np.max(series["held_out_momentum_max"])),
        },
    }


def structure_block(series, saved_prefix):
    if series["time"].size < 2:
        return None
    content = series["window_normal"]
    flux = series["window_proper_flux"]
    ledgers = separate_ledgers(content, flux, {
        "field_energy_change": float(series["field_energy"][-1] - series["field_energy"][0]),
        "energy_change": float(series["energy"][-1] - series["energy"][0]),
        "coordinate_work": float(series["coordinate_work_integral"][-1]),
        "proper_pressure_work": float(series["proper_pressure_work_integral"][-1]),
        "reported_with_windows": True,
    })
    windows = assess_windows(content, flux)
    saved_content = None
    saved_flux = None
    saved_time = None
    content_key = saved_prefix + "window_normal"
    with np.load(episode.NPZ, allow_pickle=False) as data:
        if content_key in data:
            saved_time = np.asarray(data[saved_prefix + "time"], dtype=float)
            saved_content = np.asarray(data[content_key], dtype=float)
            saved_flux = np.asarray(data[saved_prefix + "window_proper_flux"], dtype=float)
    joined = None
    if saved_content is not None:
        _time, joined_content = join_series(saved_time, saved_content, series["time"], content)
        _time, joined_flux = join_series(saved_time, saved_flux, series["time"], flux)
        joined = separate_ledgers(joined_content, joined_flux, {
            "slice": "joined_from_T0",
            "reported_with_windows": True,
        })
    start_leader = int(np.argmax(content[0]))
    end_leader = int(np.argmax(content[-1]))
    persists = bool(
        end_leader == start_leader
        and float(content[-1, start_leader]) > 0.0
        and float(content[-1, start_leader]) > 0.5 * float(content[0, start_leader])
    )
    return {
        "slice_ledgers": ledgers,
        "slice_windows": {
            key: windows[key]
            for key in (
                "start_leader", "end_leader", "same_leader", "end_share", "start_share",
                "leader_content_change", "reversal", "localized_end", "maintained",
                "renewed", "stable_throughflow", "candidate_regime",
                "packet_flux_end", "reservoir_flux_end", "flux_sum_end",
                "maintained_with_throughflow",
            )
        },
        "joined_ledgers": joined,
        "structure_persists_without_share_proxy": persists,
        "share_proxy_required": False,
        "reversal_proxy_required": False,
        "proxy_candidate_regime_is_programme_requirement": False,
        "saved_handoff_content_gap": None if saved_content is None else float(
            np.max(np.abs(saved_content[-1] - content[0]))
        ),
    }


def clock_block(series):
    if series["time"].size < 2:
        return None
    return {
        "protocol": (
            "Proper time advances by the trapezoid of the mean lapse N=r*L. "
            "L is the static calibration length density. The leader clock uses "
            "the same lapse restricted to the current leader window. Coordinate "
            "dt is the realized adaptive step."
        ),
        "coordinate_duration": float(series["time"][-1] - series["time"][0]),
        "proper_clock": float(series["proper_clock"][-1]),
        "leader_clock": float(series["leader_clock"][-1]),
        "packet_flux_end": float(series["packet_flux"][-1]),
        "reservoir_flux_end": float(series["reservoir_flux"][-1]),
        "packet_flux_integral": float(trap_integrate(series["time"], series["packet_flux"])),
        "field_energy_change": float(series["field_energy"][-1] - series["field_energy"][0]),
    }


def comparison_rows(runs, comparison_time):
    order = [spec["name"] for spec in RUN_PLAN]
    present = {name: runs[name] for name in order if name in runs}
    primary_space_fine = present.get("nf512_dtmax_0_0005") or present.get("nf512_dtmax_0_00025")
    primary_space_coarse = present.get("nf256_dtmax_0_0005") or present.get("nf256_dtmax_0_00025")
    pairs = []
    if "nf256_dtmax_0_0005" in present and "nf256_dtmax_0_00025" in present:
        pairs.append(("time_nf256", present["nf256_dtmax_0_0005"], present["nf256_dtmax_0_00025"]))
    if "nf512_dtmax_0_0005" in present and "nf512_dtmax_0_00025" in present:
        pairs.append(("time_nf512", present["nf512_dtmax_0_0005"], present["nf512_dtmax_0_00025"]))
    if "nf512_dtmax_0_0005" in present and "nf256_dtmax_0_0005" in present:
        pairs.append(("space_dtmax_0_0005", present["nf512_dtmax_0_0005"], present["nf256_dtmax_0_0005"]))
    if "nf512_dtmax_0_00025" in present and "nf256_dtmax_0_00025" in present:
        pairs.append(("space_dtmax_0_00025", present["nf512_dtmax_0_00025"], present["nf256_dtmax_0_00025"]))
    if primary_space_fine is not None and primary_space_coarse is not None and not any(
        name.startswith("space_") for name, _left, _right in pairs
    ):
        pairs.append(("space_available", primary_space_fine, primary_space_coarse))
    rows = []
    for pair_name, left, right in pairs:
        for effect_name, key, floor in episode.PHYSICAL_EFFECTS:
            if key.startswith("mode_") or key.startswith("window_"):
                continue
            rows.append(effect_at(pair_name, effect_name, key, floor, left, right, comparison_time))
        for index in range(4):
            rows.append(effect_at(
                pair_name,
                f"window_normal_{index}",
                "window_normal",
                1e-6,
                left,
                right,
                comparison_time,
                column=index,
            ))
    return rows


def effect_at(pair, name, key, floor, left, right, comparison_time, column=None):
    row = {
        "pair": pair,
        "name": name,
        "series": key,
        "floor": float(floor),
        "comparison_time": comparison_time,
        "rule": "movement of the changes <= 0.01 * primary change on the common frame",
        "status": "incomplete",
    }
    try:
        l0, l1 = frame_value(left, key, T0, comparison_time, column)
        r0, r1 = frame_value(right, key, T0, comparison_time, column)
    except KeyError:
        return row
    primary_change = l1 - l0
    other_change = r1 - r0
    effect = abs(primary_change)
    movement = abs(primary_change - other_change)
    status = episode.classify_movement(effect, movement, floor)
    row.update({
        "primary_change": float(primary_change),
        "other_change": float(other_change),
        "effect_scale": float(effect),
        "movement": float(movement),
        "movement_over_effect": None if effect == 0.0 else float(movement / effect),
        "status": status,
    })
    return row


def frame_value(run, key, start, end, column):
    series = run["series"]
    times = np.asarray(series["time"], dtype=float)
    values = np.asarray(series[key], dtype=float)
    if column is not None:
        values = values[:, int(column)]

    def pick(time_value):
        index = int(np.argmin(np.abs(times - float(time_value))))
        if abs(float(times[index]) - float(time_value)) > 1e-8:
            # Fall back to the last stored time at or before the request.
            eligible = np.flatnonzero(times <= float(time_value) + 1e-12)
            if eligible.size == 0:
                raise KeyError(time_value)
            index = int(eligible[-1])
        return float(values[index])

    return pick(start), pick(end)


def common_comparison_time(runs):
    if not runs:
        return None
    limits = []
    for run in runs.values():
        times = np.asarray(run["series"]["time"], dtype=float)
        if times.size:
            limits.append(float(times[-1]))
    if not limits:
        return None
    limit = min(limits)
    marks = np.round(np.arange(T0, limit + 1e-12, FRAME), 10)
    marks = marks[marks <= limit + 1e-12]
    if marks.size == 0:
        return float(limit)
    return float(marks[-1])


def pilot_decision(pilot):
    series = pilot["series"]
    feasible = bool(
        pilot["steps_completed"] >= 1
        and pilot["error"] is None
        and pilot["occupations_unchanged"]
        and not pilot["source_reset"]
        and float(np.max(series["gram_gap"])) <= CAR_GAP_MAX
        and bool(np.all(series["positive_r"] > 0.5))
        and bool(np.all(series["positive_Q"] > 0.5))
        and bool(np.all(series["positive_L"] > 0.5))
    )
    named = []
    for effect_name, key, floor in episode.PHYSICAL_EFFECTS:
        if key.startswith("mode_") or key not in series:
            continue
        values = np.asarray(series[key], dtype=float)
        if values.ndim != 1 or values.size < 2:
            continue
        change = abs(float(values[-1] - values[0]))
        if change >= float(floor):
            named.append({"name": effect_name, "change": change, "floor": float(floor)})
    return {
        "numerically_feasible": feasible,
        "named_effect_above_floor": bool(named),
        "named_effects": named[:12],
        "open_remaining_cases": bool(feasible and named),
        "stop_reason": pilot["stop_reason"],
        "comparison_target": float(pilot["attained_time"]),
    }


def goal_block(runs, balances, structures, clocks, decision, effects):
    reasons = []
    if not runs:
        reasons.append("no kept evolution")
    confirmation = runs.get("nf512_dtmax_0_00025") or runs.get("nf256_dtmax_0_0005")
    if confirmation is None:
        reasons.append("no confirmation trajectory")
        confirmation_name = None
    else:
        confirmation_name = confirmation["name"]
    balance = None if confirmation_name is None else balances.get(confirmation_name)
    structure = None if confirmation_name is None else structures.get(confirmation_name)
    clock = None if confirmation_name is None else clocks.get(confirmation_name)
    admissible = True
    positive = True
    field_ok = True
    for run in runs.values():
        series = run["series"]
        if series["time"].size == 0:
            admissible = False
            continue
        if float(np.max(series["gram_gap"])) > CAR_GAP_MAX:
            admissible = False
        if not (bool(np.all(series["positive_r"] > 0.5)) and bool(np.all(series["positive_Q"] > 0.5)) and bool(np.all(series["positive_L"] > 0.5))):
            positive = False
        if run["source_reset"] or not run["occupations_unchanged"]:
            field_ok = False
        if run["initial_hash"] is None:
            admissible = False
    if not admissible:
        reasons.append("Gram or handoff admissibility failed on a kept sample")
    if not positive:
        reasons.append("positive chart was lost on a kept sample")
    if not field_ok:
        reasons.append("source or occupations were reset")
    energy_closed = bool(
        balance
        and balance["balance_versus_coordinate_work"]["within_one_percent_of_exchange"]
        and balance["window_normal_within_one_percent"]
    )
    if not energy_closed:
        reasons.append("coordinate or normal-window ledger is outside one percent of its exchange")
    transfer = bool(
        clock
        and abs(clock["field_energy_change"]) >= MEANINGFUL_EXCHANGE
        and abs(clock["packet_flux_end"]) >= MEANINGFUL_EXCHANGE
        and abs(clock["proper_clock"]) >= MEANINGFUL_CLOCK
    )
    if not transfer:
        reasons.append("the new interval does not show a meaningful transfer on the declared floors")
    persists = bool(structure and structure["structure_persists_without_share_proxy"])
    renewed = bool(structure and structure["joined_ledgers"] and structure["joined_ledgers"]["reversal"]["renewed"])
    if not persists and not renewed:
        reasons.append("localized structure did not persist and the joined series is not a renewal")
    radial_only = bool(clock and abs(clock["field_energy_change"]) < MEANINGFUL_EXCHANGE)
    claimed = [row for row in effects if row.get("effect_scale") not in (None, 0.0)]
    bad_effects = [row for row in claimed if row.get("status") not in ("resolved", "below_noise_floor")]
    if bad_effects:
        reasons.append("a claimed physical comparison is unresolved or incomplete")
    goal = bool(
        not reasons
        and decision
        and decision.get("open_remaining_cases")
        and "nf512_dtmax_0_00025" in runs
    )
    return {
        "goal_complete": goal,
        "reasons_if_incomplete": reasons,
        "unresolved_comparisons": [
            {"pair": row.get("pair"), "name": row.get("name"), "status": row.get("status")}
            for row in bad_effects
        ],
        "radial_motion_alone_is_sufficient": False,
        "radial_motion_without_exchange": radial_only,
        "structure_persists": persists,
        "joined_renewal_proxy": renewed,
        "renewal_proxy_required": False,
        "two_transfer_episodes": bool(transfer and admissible and positive and field_ok),
        "episode_one": "saved spherical feedback episode T=0 to T=0.05",
        "episode_two": confirmation_name,
        "admissible_output_input": admissible,
        "positive_geometry": positive,
        "field_state_kept": field_ok,
        "energy_work_within_one_percent": energy_closed,
        "meaningful_transfer": transfer,
        "incoming_gate_required": False,
        "lambda_cdm_required": False,
        "global_eternity_required": False,
    }


def payload_schema():
    return {
        "cadence": "stored frames are the 0.005 coordinate marks from the handoff through the attained comparison window",
        "coarse_nodal_rate": "frame_Q_dot and frame_p_chi_dot are the owned Galerkin pullbacks at that frame, not a difference of the state",
        "neighbor_increment": (
            "frame_increment_Q_dot and frame_increment_p_chi are the realized change of those "
            "coarse rates across one accepted neighbor step. The handoff frame uses the following "
            "step. Later frames use the step that landed on the frame. "
            "frame_indicator_* divides that increment by frame_increment_dt."
        ),
        "indicator_versus_bound": (
            "The indicator is the average of the rate's time derivative over that one realized step. "
            "It is not a bound on the metric curvature, and it is not the Euler-Lagrange substitution "
            "of p_chi_dot into R_h = chi+2. A bound would dominate the remainder on every admissible "
            "step. This record does not supply that domination. An independent metric formula can "
            "compare itself with the indicator. Agreement is a separate measurement."
        ),
        "chi_plus_2_is_not_the_only_curvature_input": True,
        "metric_curvature_computed_here": False,
        "quadrature_fields": "frame_quad_r, frame_quad_Q, frame_quad_proper, frame_quad_K_perp, frame_quad_chi",
        "static_clock": "null_expansion.static_clock(nq, period=8) rebuilds L and beta; they are not evolved",
        "source_columns": "frame_phi0 and frame_phi1 when the payload allows; final_phi0 and final_phi1 always",
        "held_out_if_payload_exceeds_64MiB": "frame_phi0 and frame_phi1 are dropped first",
        "spatial_energy_windows": "window_* series, separate from full/projected/held-out mode projectors",
    }


def dt_frequency(values):
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return []
    rounded = np.round(values, 10)
    unique, counts = np.unique(rounded, return_counts=True)
    return [
        {"dt": float(dt), "count": int(count)}
        for dt, count in sorted(zip(unique, counts), key=lambda item: -item[1])
    ]


def build_record(meta, arrays_for_runs):
    runs = arrays_for_runs
    balances = {name: balance_block(run["series"]) for name, run in runs.items()}
    structures = {
        name: structure_block(run["series"], run["handoff"] + "_")
        for name, run in runs.items()
    }
    clocks = {name: clock_block(run["series"]) for name, run in runs.items()}
    comparison_time = common_comparison_time(runs)
    effects = comparison_rows(runs, comparison_time) if len(runs) >= 2 else []
    statuses = [row["status"] for row in effects]
    decision = meta.get("decision")
    goal = goal_block(runs, balances, structures, clocks, decision, effects)
    public_runs = {}
    for name, run in runs.items():
        series = run["series"]
        public_runs[name] = {
            "role": run["role"],
            "nf": run["nf"],
            "dt_cap": run["dt_cap"],
            "handoff": run["handoff"],
            "attained_T": run["attained_time"],
            "target_T": run["target_time"],
            "completed_to_target": run["completed_to_target"],
            "stop_reason": run["stop_reason"],
            "chart_stop": run["chart_stop"],
            "steps_completed": run["steps_completed"],
            "rejected_count": run["rejected_count"],
            "rejected": run["rejected"],
            "cpu_seconds": run["cpu_seconds"],
            "limiter_counts": run["limiter_counts"],
            "dt_frequency": dt_frequency(run["dt_values"]),
            "dt_max_realized": None if run["dt_values"].size == 0 else float(np.max(run["dt_values"])),
            "dt_min_realized": None if run["dt_values"].size == 0 else float(np.min(run["dt_values"])),
            "omega_max": None if run["omega_values"].size == 0 else float(np.max(run["omega_values"])),
            "dt_omega_max": None if run["dt_omega_values"].size == 0 else float(np.max(run["dt_omega_values"])),
            "safety_cap_after_step_count": run["safety_cap_after_step_count"],
            "safety_cap_is_physical_failure": False,
            "gram_gap_max": None if series["time"].size == 0 else float(np.max(series["gram_gap"])),
            "positive_r": None if series["time"].size == 0 else bool(np.all(series["positive_r"] > 0.5)),
            "positive_Q": None if series["time"].size == 0 else bool(np.all(series["positive_Q"] > 0.5)),
            "positive_L": None if series["time"].size == 0 else bool(np.all(series["positive_L"] > 0.5)),
            "initial_hash": run["initial_hash"],
            "final_hash": run["final_hash"],
            "phi_max_abs_change": run["phi_max_abs_change"],
            "source_reset": run["source_reset"],
            "occupations_unchanged": run["occupations_unchanged"],
            "initial_sample": None if series["time"].size == 0 else {
                "time": float(series["time"][0]),
                "Q_min": float(series["Q_min"][0]),
                "Q_max": float(series["Q_max"][0]),
                "r_min": float(series["r_min"][0]),
                "r_max": float(series["r_max"][0]),
                "chi_max": float(series["chi_max"][0]),
                "proper_min": float(series["proper_min"][0]),
                "proper_max": float(series["proper_max"][0]),
                "leader": int(series["leader"][0]),
                "leader_share": float(series["leader_share"][0]),
                "leader_content": float(series["leader_content"][0]),
                "field_energy": float(series["field_energy"][0]),
                "q_min": float(series["q_min"][0]),
                "q_at_areal_max": float(series["q_at_areal_max"][0]),
            },
            "final_sample": None if series["time"].size == 0 else {
                "time": float(series["time"][-1]),
                "Q_min": float(series["Q_min"][-1]),
                "Q_max": float(series["Q_max"][-1]),
                "r_min": float(series["r_min"][-1]),
                "r_max": float(series["r_max"][-1]),
                "chi_min": float(series["chi_min"][-1]),
                "chi_max": float(series["chi_max"][-1]),
                "proper_min": float(series["proper_min"][-1]),
                "proper_max": float(series["proper_max"][-1]),
                "leader": int(series["leader"][-1]),
                "leader_share": float(series["leader_share"][-1]),
                "leader_content": float(series["leader_content"][-1]),
                "field_energy": float(series["field_energy"][-1]),
                "q_min": float(series["q_min"][-1]),
                "q_max": float(series["q_max"][-1]),
                "q_at_areal_max": float(series["q_at_areal_max"][-1]),
                "packet_flux": float(series["packet_flux"][-1]),
                "reservoir_flux": float(series["reservoir_flux"][-1]),
                "G_min": float(series["G_min"][-1]),
                "positive_G": bool(series["positive_G"][-1] > 0.5),
            },
            "expansion_frames": [
                {
                    "time": row["time"],
                    "flags": row["flags"],
                    "raw_counts": row["frame"].get("raw_counts"),
                    "strict_counts": row["frame"].get("strict_counts"),
                    "product_min": row["frame"].get("product_min"),
                    "product_max": row["frame"].get("product_max"),
                    "trapped_intervals": row["frame"].get("trapped_intervals"),
                    "anti_trapped_intervals": row["frame"].get("anti_trapped_intervals"),
                    "areal_maximum": row["frame"].get("areal_maximum"),
                    "areal_minimum": row["frame"].get("areal_minimum"),
                    "q_at_areal_max": row["q_at_areal_max"],
                }
                for row in run["expansion"]
            ],
            "endpoint_closure": run["endpoint_closure"],
            "error": run["error"],
        }
    return {
        "schema": SCHEMA,
        "protocol_id": PROTOCOL_ID,
        "status": meta.get("status"),
        "threadpool_limit": THREADPOOL_LIMIT,
        "single_thread_blas": True,
        "incoming_gate_prerequisite": False,
        "new_force_added": False,
        "state_reset": False,
        "mean_source_subtracted": False,
        "old_toy_B_matched": False,
        "gamma_one": "the owned spherical feedback action; no second action was inserted",
        "occupations": "original six Gaussian weights",
        "proxy_reversal_required": PROXY_REVERSAL_IS_REQUIRED,
        "proxy_share_required": PROXY_SHARE_IS_REQUIRED,
        "safety_cap": SAFETY_CAP,
        "absolute_rk4": ABSOLUTE_RK4,
        "safety_cap_adapts_dt": True,
        "safety_cap_is_physical_failure_label": False,
        "t0": T0,
        "requested_T": TARGET,
        "comparison_time": comparison_time,
        "budgets": {"coarse": COARSE_BUDGET, "fine": FINE_BUDGET, "payload_bytes_limit": PAYLOAD_LIMIT},
        "pools": meta.get("pools"),
        "forecast": meta.get("forecast"),
        "preflight": meta.get("preflight"),
        "decision": decision,
        "sibling_T005_gaps": meta.get("sibling_gaps"),
        "link_reports": meta.get("link_reports"),
        "settings": meta.get("settings"),
        "margin": meta.get("margin"),
        "sealed_hashes_before": meta.get("sealed_before"),
        "runs": public_runs,
        "balances": balances,
        "structures": structures,
        "clocks": clocks,
        "effects": effects,
        "effect_status_summary": episode.summarize_statuses(statuses) if statuses else "incomplete",
        "goal": goal,
        "payload_schema": payload_schema(),
        "indicator_versus_bound": payload_schema()["indicator_versus_bound"],
        "what_this_does_not_claim": [
            "A flux sum of zero is not a proper window balance.",
            "The 0.10 reversal and 0.50 share proxies are not programme requirements.",
            "Crossing the 1.4 safety cap is a timestep reduction, not a physical failure.",
            "Radial motion by itself is not regeneration and is not goal completion.",
            "chi+2 is not an independent metric-curvature measurement.",
            "The finite-step rate increment is an indicator, not a curvature bound.",
            "No global horizon, cosmological fit, or eternity is claimed.",
            "A pilot failure is a measured exit, not evidence that the regime cannot exist.",
        ],
    }


def render_note(record):
    goal = record.get("goal") or {}
    decision = record.get("decision") or {}
    lines = [
        "# Regeneration episode from the saved T=0.05 state",
        "",
        "One continuation of the owned spherical Galerkin step. The source,",
        "the feedback action, and the chart are unchanged. No restoring force,",
        "pulse, mean subtraction, or state reset was added. The incoming gate,",
        "a cosmological fit, and a global eternity are not requirements.",
        "",
        f"Status `{record.get('status')}`. Goal complete: `{goal.get('goal_complete')}`.",
        "",
        "## Handoff",
        "",
        "Each resolution starts from the saved `dt=0.0005` final at T=0.05.",
        "The half-step finals are siblings. They are not the initial data.",
        "",
    ]
    for name, run in (record.get("runs") or {}).items():
        lines.append(
            f"- `{name}` hash `{run.get('initial_hash')}`, attained T={run.get('attained_T')}, "
            f"stop `{run.get('stop_reason')}`, steps {run.get('steps_completed')}, "
            f"CPU {run.get('cpu_seconds')}."
        )
    lines.extend([
        "",
        "## Protocol",
        "",
        "Coarse budget 600 CPU seconds. Fine budget 1800 CPU seconds. Payload 64 MiB.",
        "The pilot cap is 0.0005 and the confirmation cap is 0.00025. Each step",
        "uses the minimum of that cap, the owned stable timestep, and the distance",
        "to the next 0.005 frame. The safety factor 1.4 reduces the next dt.",
        f"The absolute RK4 limit is {ABSOLUTE_RK4}. A step past that limit is rejected",
        "and the last valid state is kept. Gram admissibility is 1e-8.",
        "Declared geometric stops, evaluated on the stored frames and not on the",
        "0.10 or 0.50 proxies: a both-positive arc entering the packet [0, 4),",
        "a both-negative arc entering the bridge [4, 8), a third same-sign arc of",
        f"at least {NEW_ARC_MIN_NODES} nodes, or a turn-back of q=rQ at the areal maximum",
        f"of at least {Q_TURNBACK_FLOOR} across {Q_TURNBACK_FRAMES} frames.",
        "",
        "## Forecast",
        "",
        "```json",
        json.dumps(jsonable(record.get("forecast")), indent=2),
        "```",
        "",
        "## What the trajectories do",
        "",
    ])
    for name, run in (record.get("runs") or {}).items():
        initial = run.get("initial_sample") or {}
        final = run.get("final_sample") or {}
        clock = (record.get("clocks") or {}).get(name) or {}
        lines.append(
            f"`{name}` runs from T={initial.get('time')} to T={final.get('time')}. "
            f"Q is [{final.get('Q_min')}, {final.get('Q_max')}], "
            f"r is [{final.get('r_min')}, {final.get('r_max')}], "
            f"chi_max={final.get('chi_max')}, proper velocity "
            f"[{final.get('proper_min')}, {final.get('proper_max')}]. "
            f"Leader window {final.get('leader')} share {final.get('leader_share')}, "
            f"shell content {final.get('leader_content')}. "
            f"Packet flux {final.get('packet_flux')}, reservoir flux {final.get('reservoir_flux')}. "
            f"Proper clock {clock.get('proper_clock')}, leader clock {clock.get('leader_clock')}, "
            f"field-energy change {clock.get('field_energy_change')}. "
            f"q at the areal maximum goes from {initial.get('q_at_areal_max')} to {final.get('q_at_areal_max')}."
        )
        lines.append("")
    lines.extend([
        "The pilot decision is recorded below. Remaining cases open only when the",
        "pilot stays in the positive chart with an admissible Gram and at least one",
        "named physical change exceeds its predeclared floor.",
        "",
        "```json",
        json.dumps(jsonable(decision), indent=2),
        "```",
        "",
        "## Ledger",
        "",
        "Coordinate balance compares the change in field plus gravity energy with",
        "the integrated coordinate fieldwork. The normal-window balance compares",
        "the change in shell content with the integral of observer boundary flux,",
        "pressure work, and lapse-gradient exchange. A periodic flux sum of zero",
        "is a separate partition identity and is not that balance. Mode-projector",
        "residuals stay in the full, projected, and held-out columns.",
        "",
    ])
    for name, balance in (record.get("balances") or {}).items():
        if not balance:
            continue
        lines.append(
            f"`{name}` total-energy change {balance.get('total_energy_change')}, "
            f"coordinate work {balance.get('coordinate_work_integral')}, "
            f"pressure work {balance.get('proper_pressure_work_integral')}, "
            f"lapse exchange {balance.get('lapse_gradient_work_integral')}. "
            f"Coordinate closure within one percent: "
            f"`{balance['balance_versus_coordinate_work']['within_one_percent_of_exchange']}`. "
            f"Normal-window closure within one percent: `{balance.get('window_normal_within_one_percent')}`."
        )
        lines.append("")
    lines.extend([
        "## Comparisons",
        "",
        f"Common comparison time: `{record.get('comparison_time')}`. "
        f"Effect summary: `{record.get('effect_status_summary')}`.",
        "",
        "| Pair | Effect | Primary change | Movement | Status |",
        "|---|---|---:|---:|---|",
    ])
    for row in record.get("effects") or []:
        if row.get("status") == "incomplete" and row.get("effect_scale") is None:
            continue
        lines.append(
            f"| {row.get('pair')} | {row.get('name')} | {row.get('primary_change')} | "
            f"{row.get('movement')} | {row.get('status')} |"
        )
    lines.extend([
        "",
        "## Indicator and bound",
        "",
        record.get("indicator_versus_bound") or "",
        "",
        "The payload arrays are `frame_Q_dot`, `frame_p_chi_dot`, `frame_increment_Q_dot`,",
        "`frame_increment_p_chi`, `frame_increment_dt`, `frame_indicator_Q_dot`, and",
        "`frame_indicator_p_chi`, together with the coarse geometry and momenta and the",
        "quadrature radius, conformal factor, proper velocity, and K_perp. Static lapse",
        "and shift are rebuilt by `static_clock`. Phi columns are stored on the frames",
        "when the 64 MiB payload allows.",
        "",
        "## Goal",
        "",
        "```json",
        json.dumps(jsonable(goal), indent=2),
        "```",
        "",
        "## What this does not claim",
        "",
    ])
    for item in record.get("what_this_does_not_claim") or []:
        lines.append(f"- {item}")
    lines.extend([
        "",
        "```sh",
        "python scripts/lab.py scripts/derive_nsc_regeneration_episode.py",
        "python scripts/lab.py scripts/derive_nsc_regeneration_episode.py --check",
        "python scripts/lab.py -m pytest tests/test_nsc_regeneration_episode.py -q",
        "```",
        "",
    ])
    return "\n".join(lines)


def series_from_arrays(prefix, arrays):
    series = {"time": np.asarray(arrays[prefix + "time"], dtype=float)}
    for key in SERIES_SCALARS:
        series[key] = np.asarray(arrays[prefix + key], dtype=float)
    for key in WINDOW_CHANNELS:
        series[key] = np.asarray(arrays[prefix + key], dtype=float)
    return series


def persist(meta, runs, arrays):
    packed = dict(arrays)
    for run in runs.values():
        packed.update(pack_run(run))
    if not packed:
        packed = {"marker": np.zeros(0, dtype=np.int8)}
    # Phi frames are the payload valve. Drop them if the compressed file is over the limit.
    write_npz(NPZ, packed)
    if NPZ.stat().st_size > PAYLOAD_LIMIT:
        dropped = [key for key in list(packed) if key.endswith("frame_phi0") or key.endswith("frame_phi1")]
        for key in dropped:
            del packed[key]
        meta["payload_held_out"] = dropped
        write_npz(NPZ, packed)
    else:
        meta["payload_held_out"] = []
    record = build_record(meta, runs)
    record["payload_held_out"] = meta["payload_held_out"]
    record["npz_sha256"] = sha256_file(NPZ)
    record["driver_sha256"] = sha256_file(DRIVER)
    record["sealed_hashes_after"] = sealed_hashes()
    write_json(OUT, record)
    record["payload_bytes"] = int(OUT.stat().st_size + NPZ.stat().st_size)
    record["payload_within_64MiB"] = bool(record["payload_bytes"] <= PAYLOAD_LIMIT)
    write_json(OUT, record)
    write_text(NOTE, render_note(record))
    # The note is owned output. Its bytes are not part of the predecessor seal.
    record["note_sha256"] = sha256_file(NOTE)
    record["payload_bytes"] = int(OUT.stat().st_size + NPZ.stat().st_size + NOTE.stat().st_size)
    record["payload_within_64MiB"] = bool(
        OUT.stat().st_size + NPZ.stat().st_size <= PAYLOAD_LIMIT
    )
    write_json(OUT, record)
    return record


def campaign():
    if OUT.is_file():
        try:
            existing = json.loads(OUT.read_text())
        except json.JSONDecodeError:
            existing = {}
        if existing.get("protocol_id") == PROTOCOL_ID and existing.get("status") in ("FINAL", "PARTIAL"):
            print("record exists; not repeating the evolution", existing.get("status"), flush=True)
            return existing
    sealed_before = sealed_hashes()
    margin, margin_meta = expansion_margin()
    pools = {"coarse": Pool(COARSE_BUDGET), "fine": Pool(FINE_BUDGET)}
    grids = {}
    states = {}
    preflight = {}
    print("forecast probes, discarded", flush=True)
    for spec in (RUN_PLAN[0], RUN_PLAN[2]):
        grid, state = load_handoff(spec)
        grids[spec["nf"]] = grid
        states[spec["nf"]] = state
        probe = one_step_probe(grid, state, spec["dt_cap"])
        pools[spec["pool"]].add(probe["cpu_seconds"])
        preflight[spec["name"]] = probe
        print("probe", spec["name"], "cpu", round(probe["cpu_seconds"], 4), "omega", probe["omega"], "identity", probe["identity_relative"], flush=True)
        if probe["sign_error"]:
            meta = {
                "status": "BLOCKED_LEDGER_SIGN",
                "preflight": preflight,
                "pools": {name: {"used": pool.used, "limit": pool.limit} for name, pool in pools.items()},
                "sealed_before": sealed_before,
                "margin": margin_meta,
                "decision": None,
            }
            return persist(meta, {}, {})
    forecast = {}
    for spec in RUN_PLAN:
        probe_name = "nf256_dtmax_0_0005" if spec["nf"] == 256 else "nf512_dtmax_0_00025"
        probe = preflight[probe_name]
        # Same relative growth on both resolutions. The published anchor is the nf512 frequency.
        scaled_omega = float(probe["omega"])
        duration = TARGET - T0
        models = forecast_models(scaled_omega, spec["dt_cap"], probe["cpu_seconds"], duration)
        forecast[spec["name"]] = {
            "measured_step_cpu": probe["cpu_seconds"],
            "omega0": scaled_omega,
            "models": models,
            "pessimistic_cpu": max(row["cpu_seconds"] for row in models.values()),
        }
        print(
            "forecast", spec["name"],
            "pessimistic_cpu", round(forecast[spec["name"]]["pessimistic_cpu"], 2),
            flush=True,
        )
    sibling = {
        "nf512": sibling_gap("nf512_dt_0_0005", "nf512_dt_0_00025"),
        "nf256": sibling_gap("nf256_dt_0_0005", "nf256_dt_0_00025"),
        "handoff_is_the_dt_0_0005_final": True,
    }
    link_reports = {}
    for fermions, state in states.items():
        started = time.process_time()
        link_reports[str(fermions)] = link_report(grids[fermions], state)
        pools["coarse" if fermions == 256 else "fine"].add(time.process_time() - started)
    settings_hash, settings_payload = settings_digest(grids[512])
    meta = {
        "status": "PARTIAL",
        "preflight": preflight,
        "forecast": forecast,
        "sibling_gaps": sibling,
        "link_reports": link_reports,
        "settings": {"sha256": settings_hash, "payload": settings_payload},
        "margin": margin_meta,
        "sealed_before": sealed_before,
        "decision": None,
        "pools": None,
    }
    runs = {}
    target = TARGET
    for index, spec in enumerate(RUN_PLAN):
        if index > 0:
            decision = meta["decision"]
            if not decision or not decision["open_remaining_cases"]:
                print("remaining cases stay closed", flush=True)
                break
            target = snap_time(min(TARGET, float(decision["comparison_target"])))
        if spec["nf"] not in grids:
            grids[spec["nf"]], states[spec["nf"]] = load_handoff(spec)
        print("launch", spec["name"], "target", target, "pool", round(pools[spec["pool"]].used, 3), flush=True)
        result = evolve_case(
            spec, grids[spec["nf"]], states[spec["nf"]], target, pools[spec["pool"]], margin, spec["name"],
        )
        runs[spec["name"]] = result
        if spec["role"] == "coarse_pilot":
            meta["decision"] = pilot_decision(result)
            if result["stop_reason"] in GEOMETRIC_EVENTS:
                meta["decision"]["comparison_target"] = float(result["attained_time"])
            elif result["chart_stop"] or result["stop_reason"] in ("GRAM_LEFT_IDENTITY", "ABSOLUTE_RK4_DOMAIN", "CPU_BUDGET"):
                meta["decision"]["comparison_target"] = float(result["attained_time"])
        meta["pools"] = {name: {"used": pool.used, "limit": pool.limit} for name, pool in pools.items()}
        meta["status"] = "PARTIAL"
        persist(meta, runs, {})
    meta["pools"] = {name: {"used": pool.used, "limit": pool.limit} for name, pool in pools.items()}
    meta["status"] = "FINAL"
    meta["sealed_before"] = sealed_before
    record = persist(meta, runs, {})
    drifted = [
        name for name, digest in sealed_before.items()
        if record["sealed_hashes_after"].get(name) != digest
    ]
    if drifted:
        raise RuntimeError("sealed predecessor changed: " + ", ".join(drifted))
    print("final", record["goal"]["goal_complete"], "payload", record["payload_bytes"], flush=True)
    return record


def _sealed_binding_errors(record, *, source_ref=None):
    """Strict current binding, or exact formula-source replay at a named commit."""
    declared = record.get("sealed_hashes_before")
    if not isinstance(declared, dict):
        return ["sealed_hashes_before"]
    errors = sorted(set(declared) ^ set(SEALED_FILES))
    for name, path in SEALED_FILES.items():
        expected = declared.get(name)
        if expected is None:
            errors.append(name)
            continue
        if source_ref is not None and path.suffix == ".py":
            try:
                resolve_pinned_source_bytes(
                    ROOT, path.relative_to(ROOT).as_posix(), expected,
                    commit=source_ref)
            except (ValueError, RuntimeError):
                errors.append(name)
        elif sha256_file(path) != expected:
            errors.append(name)
    return sorted(set(errors))


def verify_saved(replay=True, *, source_ref=None):
    """Read the saved record and discard one current prescribed-gauge replay step.

    Historical formula bytes are authenticated only with an explicit full
    ``source_ref``. The numerical replay uses the current default prescribed
    implementation and checks its first-step observables against saved samples;
    it does not execute the archived producer. New requests remain strict.
    """
    if not OUT.is_file() or not NPZ.is_file():
        raise FileNotFoundError("regeneration episode evidence is missing")
    before = sealed_hashes()
    owned_before = {path: sha256_file(path) for path in OWNED_OUTPUTS}
    record = json.loads(OUT.read_text())
    if record["schema"] != SCHEMA:
        raise AssertionError("schema mismatch")
    if record["state_reset"] is not False or record["new_force_added"] is not False:
        raise AssertionError("the record declares a reset or a new force")
    if record["incoming_gate_prerequisite"] is not False:
        raise AssertionError("incoming gate was treated as a prerequisite")
    if record["safety_cap_is_physical_failure_label"] is not False:
        raise AssertionError("safety cap was labeled as a physical failure")
    if OUT.stat().st_size + NPZ.stat().st_size > PAYLOAD_LIMIT:
        raise AssertionError("payload exceeds 64MiB")
    binding_errors = _sealed_binding_errors(record, source_ref=source_ref)
    if binding_errors:
        raise AssertionError("sealed file changed or unavailable: " + ", ".join(binding_errors))
    with np.load(NPZ, allow_pickle=False) as payload:
        for spec in RUN_PLAN:
            name = spec["name"]
            if name + "_final_Q" not in payload:
                continue
            saved = load_episode_final(episode.NPZ, spec["handoff"])
            # The continuation's initial sample is the first series row. The
            # handoff state itself is recomputed from the episode final.
            if record["runs"][name]["initial_hash"] != state_sha256(saved):
                raise AssertionError(f"{name} initial hash is not the episode final")
            final = load_episode_final(NPZ, name)
            if state_sha256(final) != record["runs"][name]["final_hash"]:
                raise AssertionError(f"{name} final hash does not match the payload")
            times = np.asarray(payload[name + "_time"], dtype=float)
            if times.size == 0 or abs(float(times[0]) - T0) > 1e-12:
                raise AssertionError(f"{name} does not start at T=0.05")
            dt_values = np.asarray(payload[name + "_dt_values"], dtype=float)
            if dt_values.size and float(np.max(dt_values)) > float(spec["dt_cap"]) + 1e-12:
                raise AssertionError(f"{name} exceeded its dt cap")
            omega = np.asarray(payload[name + "_omega_values"], dtype=float)
            if dt_values.size and np.any(dt_values * omega > ABSOLUTE_RK4 + 1e-8):
                raise AssertionError(f"{name} kept a step outside the absolute RK4 limit")
            gram = np.asarray(payload[name + "_gram_gap"], dtype=float)
            if np.any(gram > CAR_GAP_MAX + 1e-15):
                raise AssertionError(f"{name} kept an inadmissible Gram")
            if not np.all(np.asarray(payload[name + "_positive_r"]) > 0.5):
                raise AssertionError(f"{name} kept a nonpositive radius")
            if not np.all(np.asarray(payload[name + "_positive_Q"]) > 0.5):
                raise AssertionError(f"{name} kept a nonpositive Q")
        if replay and "nf256_dtmax_0_0005_time" in payload:
            spec = RUN_PLAN[0]
            grid, state = load_handoff(spec)
            dt_values = np.asarray(payload["nf256_dtmax_0_0005_dt_values"], dtype=float)
            if dt_values.size:
                stepped = galerkin.rk4_step(grid, state, float(dt_values[0]))
                stored_time = np.asarray(payload["nf256_dtmax_0_0005_time"], dtype=float)
                if stored_time.size < 2:
                    raise AssertionError("replay series has no accepted step")
                # The stored series does not keep every state. Replay the hash
                # chain by checking the step leaves the handoff and stays finite.
                if state_sha256(stepped) == state_sha256(state):
                    raise AssertionError("replay step did not change the state")
                failure = chart_failure(grid.fine, galerkin.prolong_state(grid, stepped))
                if failure:
                    raise AssertionError("replay step left the chart")
                windows, partition = episode.windows_for(grid)
                replayed, _rate, _fields = episode.observe(grid, stepped, windows, partition)
                for quantity in (
                    "Q_min", "Q_max", "r_min", "r_max", "proper_min", "proper_max",
                    "field_energy", "gravity_energy",
                ):
                    saved_value = float(payload[spec["name"] + "_" + quantity][1])
                    if not np.isclose(float(replayed[quantity]), saved_value, rtol=1e-10, atol=1e-12):
                        raise AssertionError("current prescribed replay differs: " + quantity)
    for row in record.get("effects") or []:
        if row.get("effect_scale") is None or row.get("status") == "incomplete":
            continue
        fresh = episode.classify_movement(row["effect_scale"], row["movement"], row["floor"])
        if fresh != row["status"]:
            raise AssertionError(f"status drift for {row.get('name')}")
    for name, balance in (record.get("balances") or {}).items():
        if not balance:
            continue
        if balance.get("flux_sum_is_not_the_window_balance") is not True:
            raise AssertionError("flux sum was treated as the window balance")
        work = balance["balance_versus_coordinate_work"]
        fresh = episode.balance_against_exchange(work["balance_change"], work["exchange"])
        if fresh["within_one_percent_of_exchange"] != work["within_one_percent_of_exchange"]:
            raise AssertionError(f"coordinate balance flag drift for {name}")
    after = sealed_hashes()
    for name, digest in before.items():
        if after.get(name) != digest:
            raise AssertionError(f"verify wrote or changed {name}")
    owned_after = {path: sha256_file(path) for path in OWNED_OUTPUTS}
    for path, digest in owned_before.items():
        if owned_after[path] != digest:
            raise AssertionError(f"verify rewrote {path}")
    return record


def main(argv=None):
    arguments = list(argv if argv is not None else sys.argv[1:])
    if "--check" in arguments:
        record = verify_saved(replay=True, source_ref=SEALED_SOURCE_REF)
        print(record["status"], record["goal"]["goal_complete"],
              "source_ref", SEALED_SOURCE_REF, "replay_owner", "current-prescribed")
        return 0
    campaign()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
