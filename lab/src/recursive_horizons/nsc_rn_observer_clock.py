"""Transported PG-normal clock beside the frozen RN producers.

The chart is the owners' areal line element

    ds² = N² dt² − (dr + β dt)² − r² dΩ².

``nsc_rn_pg.stage_rates`` already stores two passive scalars. This module
does not write those scalars, does not advance φ, and does not replace the
geometry. It reports what those scalars are, and it integrates a different
object: the normal observer

    dr/dt = −β(t, r_obs),    dτ/dt = N(t, r_obs),

with the first-order tangent

    d(δr)/dt = −δβ − β_r δr,    d(δτ)/dt = δN + N_r δr.

A fixed-areal static clock is accepted only where g_tt = N² − β² > 0.
Leaving the closed radial node interval, or losing that timelike condition,
returns a named exit and the last valid state. Rates are not clamped to zero
and the metric is not extrapolated past excision.
"""
from __future__ import annotations

import math
import os

for _thread_name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_name, "1")

import numpy as np

from . import nsc_rn_source as source


API_VERSION = "nsc-rn-observer-clock-v1"
INTERPOLATION = "piecewise_linear_half_open_right_secant"
TIME_INTERPOLATION = "piecewise_linear_nodal_between_supplied_samples"
STRIDE_RATIO_TOLERANCE = 1.0e-6

EQUATIONS = {
    "chart": "ds^2 = N^2 dt^2 - (dr + beta dt)^2 - r^2 dOmega^2",
    "g_tt": "N^2 - beta^2",
    "g_tr": "-beta",
    "g_rr": "-1",
    "normal_radius": "dr/dt = -beta(t, r_obs)",
    "normal_proper_time": "dtau/dt = N(t, r_obs)",
    "normal_tangent_radius": "d(delta_r)/dt = -delta_beta - beta_r * delta_r",
    "normal_tangent_proper_time": "d(delta_tau)/dt = delta_N + N_r * delta_r",
    "static_proper_time": "dtau/dt = sqrt(N^2 - beta^2) only when N^2 - beta^2 > 0",
    "normal_four_velocity": "u = N^{-1} (partial_t - beta partial_r)",
    "interpolation": INTERPOLATION,
    "tie": (
        "at an interior node the value is the nodal value and the slope used "
        "is the right secant; the left secant is reported and is not averaged in"
    ),
    "edge": (
        "the left endpoint uses the first interval's secant and the right "
        "endpoint uses the last interval's secant; no exterior node is invented"
    ),
}

HISTORICAL_CLOCKS = {
    "killing_rate": "sqrt(N_out^2 - beta_out^2) at the fixed outer node if that g_tt > 0, else 0",
    "killing_meaning": (
        "static proper-time rate of the fixed outer areal node when g_tt > 0; "
        "the producer clamps the stored rate to 0 when g_tt <= 0"
    ),
    "normal_rate": "N[0] at the fixed excision node",
    "normal_meaning": (
        "lapse integral along the fixed-excision curve r = r_excision; "
        "on the vacuum RN chart that curve is spacelike because g_tt < 0 there; "
        "not proper time along the PG normal, and not a child clock"
    ),
    "normal_radius_used": False,
    "normal_is_pg_normal_proper_time": False,
    "normal_is_child_clock": False,
    "affects_phi_or_metric": False,
    "producer": "nsc_rn_pg.stage_rates",
    "producer_lines": "1002-1005, 1052, 1072-1075, 1119-1128",
    "values_are_replaced_by_this_adapter": False,
}

SAVED_HISTORY_SCOPE = {
    "station_spacing_over_rm": 0.5,
    "archive_writer": "lab/results/development/nsc-rn-neutral-family-v1/case-runner.py",
    "snapshot_constant": "SNAPSHOT_OVER_RM = 0.5",
    "stored_metric_arrays": (
        "frozen_lapse and frozen_shift of the initial sourced metric, "
        "copied at each station; coupled_mass; radius"
    ),
    "evolved_lapse_and_shift_stored": False,
    "rk4_stage_history_stored": False,
    "stored_clocks": (
        "coupled_normal, coupled_killing, fixed_normal, fixed_killing, "
        "and the matching coordinate times"
    ),
    "exact_old_run_proper_clock_recoverable": False,
    "reason": (
        "the four coupled RK4 metric stages were not saved, and the station "
        "files do not store the evolved lapse and shift"
    ),
    "recoverable_without_new_evolution": (
        "the historical normal and killing scalars already stored at each station, "
        "with the labels in HISTORICAL_CLOCKS; the frozen initial metric; "
        "an approximate transported clock only if lapse and shift samples are supplied"
    ),
    "stride_comparison": (
        "0.5 against 1 is a sampling indicator for a supplied metric history, "
        "not an error bar on the unsaved stages of an old run"
    ),
}


class ObserverClockError(ValueError):
    """The observer input is not a well-posed sample. This is not a physical exit."""

    def __init__(self, blocker):
        super().__init__(blocker)
        self.blocker = str(blocker)


def formulas():
    """Continuum relations owned by this adapter. No new field equation."""
    return {
        "api_version": API_VERSION,
        "measurement_correction": True,
        "physical_exclusion": False,
        "publication_barrier": False,
        "new_field_dynamics": False,
        "writes_producer_clocks": False,
        "affects_phi_or_metric": False,
        "equations": dict(EQUATIONS),
        "historical_clocks": dict(HISTORICAL_CLOCKS),
        "saved_history": dict(SAVED_HISTORY_SCOPE),
        "interpolation": INTERPOLATION,
        "time_interpolation": TIME_INTERPOLATION,
        "primary_initial_radius": "packet centre 12 r_m at the declared initial t0",
        "optional_second_initial_radius": "child outer edge; inward normal motion exits at excision and stops",
        "static_acceptance": "g_tt > 0 only",
        "domain": "closed interval of the supplied radial nodes; no continuation beyond excision",
    }


def _finite(value, name):
    number = float(value)
    if not math.isfinite(number):
        raise ObserverClockError(f"{name} must be finite")
    return number


def _radius_nodes(radius):
    nodes = np.array(np.asarray(radius, dtype=float), dtype=float, copy=True)
    if nodes.ndim != 1 or nodes.size < 2 or not np.all(np.isfinite(nodes)):
        raise ObserverClockError("radius nodes must be a finite array of at least two points")
    if np.any(np.diff(nodes) <= 0.0):
        raise ObserverClockError("radius nodes must be strictly increasing")
    return nodes


def _field_on_nodes(values, nodes, name):
    array = np.array(np.asarray(values, dtype=float), dtype=float, copy=True)
    if array.shape != nodes.shape or not np.all(np.isfinite(array)):
        raise ObserverClockError(f"{name} must be finite and match the radius nodes")
    return array


def _coordinates(observer):
    return {
        "t": float(observer["t"]),
        "r": float(observer["r"]),
        "tau": float(observer["tau"]),
        "delta_r": float(observer["delta_r"]),
        "delta_tau": float(observer["delta_tau"]),
    }


def _copy_observer(observer):
    copied = dict(observer)
    copied["domain"] = (float(observer["domain"][0]), float(observer["domain"][1]))
    copied["last_valid"] = None if observer["last_valid"] is None else dict(observer["last_valid"])
    copied["exit"] = None if observer["exit"] is None else dict(observer["exit"])
    return copied


def historical_producer_rates(lapse, shift):
    """Scalars stored by ``stage_rates``. The numbers are not worldline clocks.

    The killing branch matches the producer, including the clamp to zero.
    This function does not advance a state.
    """
    nodes_lapse = np.asarray(lapse, dtype=float)
    nodes_shift = np.asarray(shift, dtype=float)
    if nodes_lapse.ndim != 1 or nodes_lapse.shape != nodes_shift.shape or nodes_lapse.size < 2:
        raise ObserverClockError("historical rates need matching lapse and shift profiles")
    if not np.all(np.isfinite(nodes_lapse)) or not np.all(np.isfinite(nodes_shift)):
        raise ObserverClockError("historical lapse and shift must be finite")
    metric_tt = float(nodes_lapse[-1] ** 2 - nodes_shift[-1] ** 2)
    if metric_tt > 0.0:
        killing_rate = math.sqrt(metric_tt)
        killing_branch = "sqrt_g_tt_at_fixed_outer_node"
    else:
        killing_rate = 0.0
        killing_branch = "producer_clamped_to_zero"
    return {
        "killing_rate": float(killing_rate),
        "killing_branch": killing_branch,
        "killing_g_tt_outer": metric_tt,
        "normal_rate": float(nodes_lapse[0]),
        "normal_location": "fixed_excision_node",
        "normal_radius_used": False,
        "normal_is_pg_normal_proper_time": False,
        "normal_is_child_clock": False,
        "zero_clamp_is_not_a_timelike_refusal": True,
        "affects_phi_or_metric": False,
        "formula_killing": HISTORICAL_CLOCKS["killing_rate"],
        "formula_normal": HISTORICAL_CLOCKS["normal_rate"],
        "meaning_killing": HISTORICAL_CLOCKS["killing_meaning"],
        "meaning_normal": HISTORICAL_CLOCKS["normal_meaning"],
    }


def label_saved_clocks(clocks):
    """Attach the producer meaning without changing the stored numbers."""
    if not isinstance(clocks, dict):
        raise ObserverClockError("saved clocks must be a mapping")
    labeled = {
        "t": _finite(clocks["t"], "t"),
        "killing": _finite(clocks["killing"], "killing"),
        "normal": _finite(clocks["normal"], "normal"),
    }
    labeled.update({
        "values_unchanged": True,
        "normal_meaning": HISTORICAL_CLOCKS["normal_meaning"],
        "killing_meaning": HISTORICAL_CLOCKS["killing_meaning"],
        "normal_is_pg_normal_proper_time": False,
        "normal_is_child_clock": False,
        "killing_zero_is_a_clamp_when_present": True,
        "affects_phi_or_metric": False,
        "replaced": False,
    })
    return labeled


def _bracket(nodes, radius):
    """Interval containing ``radius``. None when the point is off the closed domain."""
    if radius < float(nodes[0]) or radius > float(nodes[-1]):
        return None
    index = int(np.searchsorted(nodes, radius, side="right"))
    if index >= nodes.size:
        left = int(nodes.size - 2)
        right = int(nodes.size - 1)
        at_node = False
        edge = "right"
    else:
        right = index
        left = index - 1
        at_node = left > 0 and radius == float(nodes[left])
        if radius == float(nodes[0]):
            edge = "left"
        else:
            edge = None
    return left, right, at_node, edge


def sample_metric(radius, lapse, shift, r, *, delta_lapse=None, delta_shift=None):
    """Piecewise-linear sample. The slope is the secant of that same piece.

    An interior node uses the right-hand piece. Both one-sided secants are
    returned on a tie. Endpoints use the only adjacent interval. Points outside
    the closed node interval are a domain exit with no extrapolated fields.
    """
    nodes = _radius_nodes(radius)
    lapse_nodes = _field_on_nodes(lapse, nodes, "lapse")
    shift_nodes = _field_on_nodes(shift, nodes, "shift")
    delta_n = None if delta_lapse is None else _field_on_nodes(delta_lapse, nodes, "delta_lapse")
    delta_b = None if delta_shift is None else _field_on_nodes(delta_shift, nodes, "delta_shift")
    point = _finite(r, "r")
    found = _bracket(nodes, point)
    if found is None:
        return {
            "status": "domain_exit",
            "reason": "radius is outside the closed nodal interval",
            "r": point,
            "domain": (float(nodes[0]), float(nodes[-1])),
            "extrapolated": False,
            "clamped": False,
            "interpolation": INTERPOLATION,
            "N": None,
            "beta": None,
            "N_r": None,
            "beta_r": None,
        }
    left, right, tie, edge = found
    span = float(nodes[right] - nodes[left])
    weight = (point - float(nodes[left])) / span

    def _piece(values):
        slope = float(values[right] - values[left]) / span
        value = float(values[left] + slope * (point - float(nodes[left])))
        return value, slope

    lapse_value, lapse_slope = _piece(lapse_nodes)
    shift_value, shift_slope = _piece(shift_nodes)
    if delta_n is None:
        delta_n_value, delta_n_slope = 0.0, 0.0
    else:
        delta_n_value, delta_n_slope = _piece(delta_n)
    if delta_b is None:
        delta_b_value, delta_b_slope = 0.0, 0.0
    else:
        delta_b_value, delta_b_slope = _piece(delta_b)
    left_lapse_slope = None
    left_shift_slope = None
    if tie:
        previous = float(nodes[left] - nodes[left - 1])
        left_lapse_slope = float(lapse_nodes[left] - lapse_nodes[left - 1]) / previous
        left_shift_slope = float(shift_nodes[left] - shift_nodes[left - 1]) / previous
    slope_side = "left_secant_of_last_interval" if edge == "right" else "right_secant"
    return {
        "status": "active",
        "reason": None,
        "r": point,
        "domain": (float(nodes[0]), float(nodes[-1])),
        "extrapolated": False,
        "clamped": False,
        "interpolation": INTERPOLATION,
        "interval": (left, right),
        "tie": bool(tie),
        "edge": edge,
        "slope_side": slope_side,
        "N": lapse_value,
        "beta": shift_value,
        "N_r": lapse_slope,
        "beta_r": shift_slope,
        "N_r_left": left_lapse_slope,
        "beta_r_left": left_shift_slope,
        "N_r_right": lapse_slope,
        "beta_r_right": shift_slope,
        "delta_N": delta_n_value,
        "delta_beta": delta_b_value,
        "delta_N_r": delta_n_slope,
        "delta_beta_r": delta_b_slope,
        "g_tt": lapse_value ** 2 - shift_value ** 2,
        "g_tr": -shift_value,
        "g_rr": -1.0,
        "piece_weight": weight,
    }


def _require_active_sample(sample):
    if sample.get("status") != "active":
        return {
            "status": "domain_exit",
            "reason": sample.get("reason", "metric sample is outside the radial domain"),
            "clamped_to_zero": False,
            "extrapolated": False,
        }
    return None


def normal_observer_rate(sample, delta_r=0.0):
    """Coordinate rates of the PG normal and its tangent at one metric sample."""
    refused = _require_active_sample(sample)
    if refused is not None:
        return refused
    lapse = float(sample["N"])
    shift = float(sample["beta"])
    if not math.isfinite(lapse) or not math.isfinite(shift):
        raise ObserverClockError("normal sample fields must be finite")
    if lapse <= 0.0:
        return {
            "status": "lapse_not_positive",
            "reason": "N <= 0, so dτ/dt = N is not an oriented proper-time rate",
            "clamped_to_zero": False,
            "extrapolated": False,
            "N": lapse,
            "beta": shift,
        }
    displacement = _finite(delta_r, "delta_r")
    lapse_r = _finite(sample["N_r"], "N_r")
    shift_r = _finite(sample["beta_r"], "beta_r")
    delta_lapse = _finite(sample.get("delta_N", 0.0), "delta_N")
    delta_shift = _finite(sample.get("delta_beta", 0.0), "delta_beta")
    radius_rate = -shift
    proper_rate = lapse
    tangent_radius = -delta_shift - shift_r * displacement
    tangent_proper = delta_lapse + lapse_r * displacement
    along = lapse ** 2 - (radius_rate + shift) ** 2
    return {
        "status": "active",
        "dr_dt": radius_rate,
        "dtau_dt": proper_rate,
        "d_delta_r_dt": tangent_radius,
        "d_delta_tau_dt": tangent_proper,
        "N": lapse,
        "beta": shift,
        "N_r": lapse_r,
        "beta_r": shift_r,
        "delta_N": delta_lapse,
        "delta_beta": delta_shift,
        "g_tt": lapse ** 2 - shift ** 2,
        "coordinate_interval": along,
        "normalization": along / (proper_rate ** 2),
        "equations": {
            "dr_dt": EQUATIONS["normal_radius"],
            "dtau_dt": EQUATIONS["normal_proper_time"],
            "d_delta_r_dt": EQUATIONS["normal_tangent_radius"],
            "d_delta_tau_dt": EQUATIONS["normal_tangent_proper_time"],
        },
        "clamped_to_zero": False,
        "extrapolated": False,
    }


def static_areal_rate(sample):
    """Proper-time rate of a fixed areal radius. Accepted only for g_tt > 0."""
    refused = _require_active_sample(sample)
    if refused is not None:
        refused["status"] = "domain_exit"
        return refused
    lapse = float(sample["N"])
    shift = float(sample["beta"])
    metric_tt = lapse ** 2 - shift ** 2
    if lapse <= 0.0 or metric_tt <= 0.0:
        return {
            "status": "static_not_timelike",
            "reason": "fixed-areal curve has g_tt <= 0 or N <= 0",
            "rate": None,
            "g_tt": metric_tt,
            "N": lapse,
            "beta": shift,
            "r": sample.get("r"),
            "clamped_to_zero": False,
            "extrapolated": False,
        }
    return {
        "status": "active",
        "reason": None,
        "rate": math.sqrt(metric_tt),
        "g_tt": metric_tt,
        "N": lapse,
        "beta": shift,
        "r": sample.get("r"),
        "clamped_to_zero": False,
        "extrapolated": False,
        "formula": EQUATIONS["static_proper_time"],
    }


def _observer_shell(*, kind, radius, domain, t, tau, delta_r, delta_tau, label):
    point = _finite(radius, "r")
    start = _finite(t, "t")
    proper = _finite(tau, "tau")
    dr = _finite(delta_r, "delta_r")
    dtau = _finite(delta_tau, "delta_tau")
    lower = _finite(domain[0], "domain_min")
    upper = _finite(domain[1], "domain_max")
    if upper <= lower:
        raise ObserverClockError("observer domain must be a positive radial interval")
    observer = {
        "kind": kind,
        "label": str(label),
        "status": "active",
        "t": start,
        "r": point,
        "tau": proper,
        "delta_r": dr,
        "delta_tau": dtau,
        "domain": (lower, upper),
        "initial_t": start,
        "initial_r": point,
        "initial_tau": proper,
        "exit": None,
        "last_valid": None,
        "interpolation": INTERPOLATION,
        "writes_producer_clocks": False,
        "affects_phi_or_metric": False,
        "clock_advanced": False,
    }
    if point < lower or point > upper:
        observer["status"] = "domain_exit"
        observer["exit"] = {
            "kind": "domain_exit",
            "reason": "declared initial radius is outside the closed radial domain",
            "attempted_r": point,
            "attempted_t": start,
        }
        observer["last_valid"] = None
        return observer
    observer["last_valid"] = _coordinates(observer)
    return observer


def make_normal_observer(radius, domain, *, t=0.0, tau=0.0, delta_r=0.0, delta_tau=0.0,
                         label="declared"):
    """Normal observer at a declared areal radius and coordinate time."""
    observer = _observer_shell(
        kind="pg_normal", radius=radius, domain=domain, t=t, tau=tau,
        delta_r=delta_r, delta_tau=delta_tau, label=label,
    )
    observer["worldline"] = EQUATIONS["normal_radius"]
    observer["proper_time"] = EQUATIONS["normal_proper_time"]
    return observer


def make_static_observer(radius, domain, *, t=0.0, tau=0.0, label="fixed_areal"):
    """Fixed-areal clock. The radius never moves; timelike checks happen on update."""
    observer = _observer_shell(
        kind="static_areal", radius=radius, domain=domain, t=t, tau=tau,
        delta_r=0.0, delta_tau=0.0, label=label,
    )
    observer["worldline"] = "dr/dt = 0"
    observer["proper_time"] = EQUATIONS["static_proper_time"]
    return observer


def initial_observers(grid, *, t0=0.0, tau0=0.0):
    """Parent packet centre, and the optional child-outer-edge normal observer.

    The centre is ``packet_request(r_m)``, which is 12 r_m on this family.
    The second observer starts at the child interval's outer edge. Both use
    the grid's closed nodal domain. An initial radius outside that domain is
    returned already stopped, without being moved onto a node.
    """
    if "r_m" not in grid or "radius" not in grid or "child_interval" not in grid:
        raise ObserverClockError("initial observers need r_m, radius, and child_interval")
    magnetic = _finite(grid["r_m"], "r_m")
    if magnetic <= 0.0:
        raise ObserverClockError("r_m must be positive")
    centre = float(source.packet_request(magnetic)["center"])
    child_outer = _finite(grid["child_interval"][1], "child_outer_edge")
    nodes = _radius_nodes(grid["radius"])
    domain = (float(nodes[0]), float(nodes[-1]))
    start = _finite(t0, "t0")
    proper = _finite(tau0, "tau0")
    parent = make_normal_observer(
        centre, domain, t=start, tau=proper, label="parent_packet_centre",
    )
    child = make_normal_observer(
        child_outer, domain, t=start, tau=proper, label="child_outer_edge",
    )
    return {
        "t0": start,
        "tau0": proper,
        "r_m": magnetic,
        "packet_centre_over_rm": float(source.PACKET_CENTER_OVER_RM),
        "parent_packet_centre": parent,
        "child_outer_edge": child,
        "domain": domain,
        "child_outer_edge_is_optional": True,
        "child_outer_edge_expected_fate_on_inward_normal": "domain_exit_at_excision_and_stop",
        "writes_producer_clocks": False,
    }


def _exit_observer(observer, kind, reason, *, attempted_r, attempted_t, stage=None):
    stopped = _copy_observer(observer)
    stopped["status"] = kind
    stopped["exit"] = {
        "kind": kind,
        "reason": reason,
        "attempted_r": None if attempted_r is None else float(attempted_r),
        "attempted_t": None if attempted_t is None else float(attempted_t),
        "stage": stage,
    }
    stopped["clock_advanced"] = False
    return stopped


def _sample_from_frame(frame, radius, stage_time):
    if callable(frame):
        raw = frame(float(stage_time), float(radius))
        if not isinstance(raw, dict):
            raise ObserverClockError("metric evaluator must return a mapping")
        if raw.get("status") not in (None, "active"):
            return raw
        required = ("N", "beta", "N_r", "beta_r")
        if any(key not in raw for key in required):
            raise ObserverClockError("metric evaluator must return N, beta, N_r, and beta_r")
        return {
            "status": "active",
            "r": float(radius),
            "N": _finite(raw["N"], "N"),
            "beta": _finite(raw["beta"], "beta"),
            "N_r": _finite(raw["N_r"], "N_r"),
            "beta_r": _finite(raw["beta_r"], "beta_r"),
            "delta_N": 0.0 if raw.get("delta_N") is None else _finite(raw["delta_N"], "delta_N"),
            "delta_beta": 0.0 if raw.get("delta_beta") is None else _finite(raw["delta_beta"], "delta_beta"),
            "extrapolated": False,
            "interpolation": raw.get("interpolation", "caller_evaluator"),
        }
    if not isinstance(frame, dict) or "lapse" not in frame or "shift" not in frame:
        raise ObserverClockError("a nodal stage needs lapse and shift")
    nodes = frame.get("radius")
    if nodes is None:
        raise ObserverClockError("a nodal stage needs radius")
    return sample_metric(
        nodes, frame["lapse"], frame["shift"], radius,
        delta_lapse=frame.get("delta_lapse"), delta_shift=frame.get("delta_shift"),
    )


def _with_radius(frame, radius):
    """Couple a stage_rates mapping to the grid radius without copying its arrays."""
    if radius is None or callable(frame) or not isinstance(frame, dict) or frame.get("radius") is not None:
        return frame
    paired = dict(frame)
    paired["radius"] = radius
    return paired


def _prepare_stages(observer, dt, stages, evaluator, radius):
    step = _finite(dt, "dt")
    if step <= 0.0:
        raise ObserverClockError("observer step must be positive")
    if (stages is None) == (evaluator is None):
        raise ObserverClockError("pass either four metric stages or one metric evaluator")
    if evaluator is not None:
        if not callable(evaluator):
            raise ObserverClockError("metric evaluator must be callable")
        frames = (evaluator, evaluator, evaluator, evaluator)
    else:
        frames = tuple(_with_radius(frame, radius) for frame in stages)
        if len(frames) != 4:
            raise ObserverClockError("RK4 needs the four coupled metric stages in order")
    times = (
        float(observer["t"]),
        float(observer["t"]) + 0.5 * step,
        float(observer["t"]) + 0.5 * step,
        float(observer["t"]) + step,
    )
    return step, frames, times


def _domain_holds(observer, radius):
    lower, upper = observer["domain"]
    return lower <= float(radius) <= upper


def rk4_normal_step(observer, dt, stages=None, *, evaluator=None, radius=None):
    """One RK4 step on the normal worldline.

    Nodal ``stages`` are the four metric frames already produced by the
    coupled step, in the same order as ``nsc_rn_pg.rk4_step``: the base
    stage, the two midpoint stages, and the full-step stage. Pass ``radius``
    when those frames are ``stage_rates`` mappings that carry lapse and shift
    but not the grid. A callable ``evaluator(t, r)`` is the alternative when
    the frames are not stored. The observer state is not fed back into φ or
    the metric.
    """
    if observer.get("kind") != "pg_normal":
        raise ObserverClockError("normal RK4 requires a pg_normal observer")
    if observer.get("status") != "active":
        return _copy_observer(observer)
    step, frames, times = _prepare_stages(observer, dt, stages, evaluator, radius)
    state = np.array([
        observer["r"], observer["tau"], observer["delta_r"], observer["delta_tau"],
    ], dtype=float)
    weights = (0.0, 0.5, 0.5, 1.0)
    slopes = []
    for index, (frame, stage_time, weight) in enumerate(zip(frames, times, weights)):
        probe = state if index == 0 else state + weight * step * slopes[-1]
        if not _domain_holds(observer, probe[0]):
            return _exit_observer(
                observer, "domain_exit",
                "RK4 stage radius left the closed nodal domain",
                attempted_r=float(probe[0]), attempted_t=stage_time, stage=index,
            )
        sample = _sample_from_frame(frame, float(probe[0]), stage_time)
        if sample.get("status") != "active":
            return _exit_observer(
                observer, "domain_exit",
                sample.get("reason", "metric sample left the radial domain"),
                attempted_r=float(probe[0]), attempted_t=stage_time, stage=index,
            )
        rate = normal_observer_rate(sample, delta_r=float(probe[2]))
        if rate["status"] != "active":
            return _exit_observer(
                observer, rate["status"], rate["reason"],
                attempted_r=float(probe[0]), attempted_t=stage_time, stage=index,
            )
        slopes.append(np.array([
            rate["dr_dt"], rate["dtau_dt"], rate["d_delta_r_dt"], rate["d_delta_tau_dt"],
        ], dtype=float))
    increment = (step / 6.0) * (slopes[0] + 2.0 * slopes[1] + 2.0 * slopes[2] + slopes[3])
    updated_state = state + increment
    if not _domain_holds(observer, updated_state[0]) or not np.all(np.isfinite(updated_state)):
        return _exit_observer(
            observer, "domain_exit",
            "RK4 update left the closed nodal domain",
            attempted_r=float(updated_state[0]), attempted_t=float(observer["t"] + step),
            stage=4,
        )
    updated = _copy_observer(observer)
    updated["t"] = float(observer["t"] + step)
    updated["r"] = float(updated_state[0])
    updated["tau"] = float(updated_state[1])
    updated["delta_r"] = float(updated_state[2])
    updated["delta_tau"] = float(updated_state[3])
    updated["status"] = "active"
    updated["exit"] = None
    updated["clock_advanced"] = True
    updated["last_valid"] = _coordinates(updated)
    updated["stage_order"] = "base, midpoint, midpoint, full_step"
    return updated


def rk4_static_step(observer, dt, stages=None, *, evaluator=None, radius=None):
    """RK4 quadrature of a fixed-areal clock. Any non-timelike stage refuses the step."""
    if observer.get("kind") != "static_areal":
        raise ObserverClockError("static RK4 requires a static_areal observer")
    if observer.get("status") != "active":
        return _copy_observer(observer)
    step, frames, times = _prepare_stages(observer, dt, stages, evaluator, radius)
    if not _domain_holds(observer, observer["r"]):
        return _exit_observer(
            observer, "domain_exit",
            "static radius is outside the closed nodal domain",
            attempted_r=observer["r"], attempted_t=observer["t"], stage=0,
        )
    rates = []
    for index, (frame, stage_time) in enumerate(zip(frames, times)):
        sample = _sample_from_frame(frame, float(observer["r"]), stage_time)
        if sample.get("status") != "active":
            return _exit_observer(
                observer, "domain_exit",
                sample.get("reason", "static sample left the radial domain"),
                attempted_r=observer["r"], attempted_t=stage_time, stage=index,
            )
        rate = static_areal_rate(sample)
        if rate["status"] != "active":
            return _exit_observer(
                observer, rate["status"], rate["reason"],
                attempted_r=observer["r"], attempted_t=stage_time, stage=index,
            )
        rates.append(float(rate["rate"]))
    proper = (step / 6.0) * (rates[0] + 2.0 * rates[1] + 2.0 * rates[2] + rates[3])
    updated = _copy_observer(observer)
    updated["t"] = float(observer["t"] + step)
    updated["tau"] = float(observer["tau"] + proper)
    updated["status"] = "active"
    updated["exit"] = None
    updated["clock_advanced"] = True
    updated["last_valid"] = _coordinates(updated)
    updated["static_stage_rates"] = rates
    return updated


def integrate_normal(observer, dt, steps, *, evaluator=None, stage_sequence=None, radius=None):
    """Repeated normal RK4. The first named exit stops the loop and is kept."""
    count = int(steps)
    if isinstance(steps, bool) or count < 1:
        raise ObserverClockError("integration needs a positive step count")
    if (evaluator is None) == (stage_sequence is None):
        raise ObserverClockError("pass either an evaluator or one stage bundle per step")
    if stage_sequence is not None and len(tuple(stage_sequence)) != count:
        raise ObserverClockError("stage_sequence length must equal the step count")
    current = _copy_observer(observer)
    accepted = 0
    for index in range(count):
        if current["status"] != "active":
            break
        if evaluator is not None:
            current = rk4_normal_step(current, dt, evaluator=evaluator)
        else:
            current = rk4_normal_step(current, dt, stages=stage_sequence[index], radius=radius)
        if current["status"] != "active":
            break
        accepted += 1
    return {
        "observer": current,
        "accepted_steps": accepted,
        "stopped": current["status"],
        "continues_after_exit": False,
    }


def _history_samples(samples):
    frames = tuple(samples)
    if len(frames) < 2:
        raise ObserverClockError("saved-history reconstruction needs at least two metric samples")
    prepared = []
    for frame in frames:
        if not isinstance(frame, dict):
            raise ObserverClockError("each history sample must be a mapping")
        nodes = _radius_nodes(frame["radius"])
        prepared.append({
            "t": _finite(frame["t"], "sample_t"),
            "radius": nodes,
            "lapse": _field_on_nodes(frame["lapse"], nodes, "lapse"),
            "shift": _field_on_nodes(frame["shift"], nodes, "shift"),
        })
    times = [item["t"] for item in prepared]
    if any(times[index] >= times[index + 1] for index in range(len(times) - 1)):
        raise ObserverClockError("history samples must be strictly increasing in time")
    reference_radius = prepared[0]["radius"]
    for item in prepared[1:]:
        if item["radius"].shape != reference_radius.shape or not np.array_equal(item["radius"], reference_radius):
            raise ObserverClockError("history samples must share one immutable radial grid")
    return prepared


def _on_stride(time_value, origin, stride):
    ratio = (float(time_value) - float(origin)) / float(stride)
    return abs(ratio - round(ratio)) <= STRIDE_RATIO_TOLERANCE


def _subset_stride(samples, stride):
    step = _finite(stride, "stride")
    if step <= 0.0:
        raise ObserverClockError("history stride must be positive")
    origin = samples[0]["t"]
    kept = [item for item in samples if _on_stride(item["t"], origin, step)]
    if len(kept) < 2:
        raise ObserverClockError("stride does not land on two supplied samples")
    return kept


def _reconstruct_one(samples, *, r0, tau0, delta_r, delta_tau, substeps):
    domain = (float(samples[0]["radius"][0]), float(samples[0]["radius"][-1]))
    observer = make_normal_observer(
        r0, domain, t=samples[0]["t"], tau=tau0, delta_r=delta_r, delta_tau=delta_tau,
        label="saved_history_approximation",
    )
    if observer["status"] != "active":
        return observer, 0
    pieces = int(substeps)
    if isinstance(substeps, bool) or pieces < 1:
        raise ObserverClockError("history substeps must be positive")
    accepted = 0
    for left, right in zip(samples, samples[1:]):
        span = float(right["t"] - left["t"])
        piece = span / pieces
        for index in range(pieces):
            start = float(left["t"] + index * piece)

            def _evaluator(time_value, radius, *, _left=left, _right=right, _span=span):
                weight = (float(time_value) - float(_left["t"])) / _span
                if weight < -1.0e-10 or weight > 1.0 + 1.0e-10:
                    raise ObserverClockError("history evaluator was called outside its sample interval")
                weight = min(1.0, max(0.0, weight))
                lapse = (1.0 - weight) * _left["lapse"] + weight * _right["lapse"]
                shift = (1.0 - weight) * _left["shift"] + weight * _right["shift"]
                return sample_metric(_left["radius"], lapse, shift, radius)

            observer = rk4_normal_step(observer, piece, evaluator=_evaluator)
            if observer["status"] != "active":
                return observer, accepted
            if abs(observer["t"] - (start + piece)) > 1.0e-8 * max(1.0, abs(start + piece)):
                raise ObserverClockError("history step did not land on the interpolated time")
            accepted += 1
    return observer, accepted


def reconstruct_saved_history(samples, *, r0, tau0=0.0, strides=(0.5, 1.0), substeps=4,
                              delta_r=0.0, delta_tau=0.0):
    """Approximate PG-normal clock from supplied metric samples.

    Time interpolation is piecewise linear between the samples kept by each
    stride. Radial interpolation is the declared piecewise-linear rule.
    The result is a sampling indicator. It is not the proper clock of an old
    run whose RK4 stages were not saved.
    """
    prepared = _history_samples(samples)
    spacing = [prepared[index + 1]["t"] - prepared[index]["t"] for index in range(len(prepared) - 1)]
    uniform = max(spacing) - min(spacing) <= 1.0e-8 * max(1.0, max(spacing))
    by_stride = {}
    for stride in strides:
        selected = _subset_stride(prepared, stride)
        observer, accepted = _reconstruct_one(
            selected, r0=r0, tau0=tau0, delta_r=delta_r, delta_tau=delta_tau, substeps=substeps,
        )
        by_stride[float(stride)] = {
            "observer": observer,
            "samples_used": len(selected),
            "sample_times": [item["t"] for item in selected],
            "accepted_substeps": accepted,
            "status": observer["status"],
        }
    stride_values = [float(item) for item in strides]
    compared = None
    if len(stride_values) >= 2:
        first = by_stride[stride_values[0]]["observer"]
        second = by_stride[stride_values[1]]["observer"]
        if first["status"] == "active" and second["status"] == "active":
            compared = {
                "stride_a": stride_values[0],
                "stride_b": stride_values[1],
                "delta_r": float(first["r"] - second["r"]),
                "delta_tau": float(first["tau"] - second["tau"]),
                "both_active": True,
            }
        else:
            compared = {
                "stride_a": stride_values[0],
                "stride_b": stride_values[1],
                "delta_r": None,
                "delta_tau": None,
                "both_active": False,
                "status_a": first["status"],
                "status_b": second["status"],
            }
    return {
        "exact_old_run_proper_clock": False,
        "claim": "sampling_indicator_only",
        "reason": SAVED_HISTORY_SCOPE["reason"],
        "time_interpolation": TIME_INTERPOLATION,
        "radial_interpolation": INTERPOLATION,
        "uniform_input_spacing": bool(uniform),
        "input_spacing": spacing,
        "substeps_per_interval": int(substeps),
        "by_stride": by_stride,
        "sampling_difference": compared,
        "archive": dict(SAVED_HISTORY_SCOPE),
        "affects_phi_or_metric": False,
    }


def saved_history_scope():
    """What a station archive can and cannot say about the transported clock."""
    return dict(SAVED_HISTORY_SCOPE)
