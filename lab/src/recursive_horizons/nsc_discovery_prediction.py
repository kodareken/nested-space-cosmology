"""Matched proper-time prediction from the existing full-state Jacobian-vector product.

The state and its tangent advance together by RK4. The tangent map is
``nested_rate_jacobian_vector``. Proper time and its tangent are sums of the
same stage samples. They are not the frozen-slice product ``δτ̇ Δt``.

The baseline Cauchy state is a saved nf=128 or nf=256 initial frame, with
that frame's frozen ``W``. The occupation direction is the trace-free
contrast ``(+,+,0,0,-,-)``. Its radius tangent is one call to
``prepared_nested_radius_tangent``. The saved ``T=0.3`` source is not
changed and its constraint is not solved again.

An independent control rebuilds changed sources with the dense
``initial_state`` constructor and evolves those states on the FFT carrier.
Derivative steps are ``0.001`` and ``0.0005``. The held-out amplitude
``0.01`` is not one of those steps.

The linear transport ``δO|t − Ȯ δτ/τ̇`` stays a Taylor indicator. The
measurement compared with the prediction stops each finite trajectory at
the baseline proper time by a bracketed, decreasing step. ``--run`` writes
that comparison. The campaign executor calls it; the tests stay short.
This module does not form the future effective stress.
"""
from __future__ import annotations

from dataclasses import replace
import io
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np

from . import nsc_discovery_backend as backend
from . import nsc_discovery_episode as episode
from . import nsc_discovery_response as discovery
from . import nsc_nested_parent_child as model
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin


SCHEMA = "NSC-DISCOVERY-PREDICTION-v1"
METHOD = "rk4_analytic_jacobian_vector"
CLOCK_QUADRATURE = "rk4_stages"
FORMULA = "deltaO_t - Odot * delta_tau / tau_dot"
PRIMARY_READOUT = discovery.PRIMARY_READOUT
SECONDARY_READOUT = discovery.SECONDARY_READOUT
CHILD_INTERVAL = tuple(model.CHILD_INTERVAL)
CHILD_CLOCK = float(discovery.CHILD_CLOCK)
HANDOFF_TIME = float(episode.HANDOFF_TIME)
NONPRODUCTION_DURATION = 0.01
DERIVATIVE_H = (0.001, 0.0005)
HELD_OUT_ALPHA = 0.01
OCCUPATION_CONTRAST = np.array([1.0, 1.0, 0.0, 0.0, -1.0, -1.0], dtype=float)
CONTRAST_PACKETS = (
    "left_parent_plus",
    "left_parent_minus",
    "child_plus",
    "child_minus",
    "right_parent_plus",
    "right_parent_minus",
)
EVOLUTION_BACKEND = "fft"
PREPARATION_BACKEND = "dense"
RADIUS_TANGENT = "prepared_nested_radius_tangent"
CONSTRUCTOR = "nsc_nested_parent_child.initial_state"
MANIFEST_NAME = "prediction-manifest.json"
RUN_JSON = "prediction-run.json"
RUN_NPZ = "prediction-run.npz"
LATER_OUTPUTS = (RUN_JSON, RUN_NPZ)
FIRST_PREDICTION_TARGET = 0.3
OPTIONAL_PREDICTION_TARGET = 1.0
EQUAL_TAU_METHOD = "bracketed_decreasing_dt"
TAYLOR_ROLE = "linear_clock_correction_not_the_equal_tau_measurement"
ROOT_BISECTIONS = 20
FINITE_DIFFERENCE_USED_AS_ANALYTIC = False
NEWTON_ON_ANALYTIC_TANGENT = False
T03_SOURCE_CHANGED = False
T03_CONSTRAINT_RESOLVED = False
EFFECTIVE_STRESS_EVALUATED = False
HARD_ACCEPTANCE_GATE = False
MISSING = (
    "ctp_future_effective_stress",
    "late_T3_curvature_unresolved",
)
_STATE_FIELDS = model.STATE_NAMES
_MODULE = Path(__file__).resolve()


def specification():
    """Contract for the prediction wave. Production evolution stays outside it."""
    return {
        "schema": SCHEMA,
        "method": METHOD,
        "clock_quadrature": CLOCK_QUADRATURE,
        "formula": FORMULA,
        "primary_readout": PRIMARY_READOUT,
        "secondary_readout": SECONDARY_READOUT,
        "child_interval": list(CHILD_INTERVAL),
        "child_clock": CHILD_CLOCK,
        "occupation_contrast": OCCUPATION_CONTRAST.tolist(),
        "contrast_packets": list(CONTRAST_PACKETS),
        "contrast_trace": float(np.sum(OCCUPATION_CONTRAST)),
        "derivative_h": list(DERIVATIVE_H),
        "held_out_alpha": HELD_OUT_ALPHA,
        "held_out_distinct_from_derivative_h": True,
        "evolution_backend": EVOLUTION_BACKEND,
        "preparation_backend": PREPARATION_BACKEND,
        "radius_tangent": RADIUS_TANGENT,
        "changed_source_constructor": CONSTRUCTOR,
        "finite_difference_used_as_analytic": FINITE_DIFFERENCE_USED_AS_ANALYTIC,
        "newton_on_analytic_tangent": NEWTON_ON_ANALYTIC_TANGENT,
        "t03_source_changed": T03_SOURCE_CHANGED,
        "t03_constraint_resolved": T03_CONSTRAINT_RESOLVED,
        "effective_stress_evaluated": EFFECTIVE_STRESS_EVALUATED,
        "effective_stress": None,
        "nonproduction_duration": NONPRODUCTION_DURATION,
        "comparison_time": HANDOFF_TIME,
        "first_prediction_target": FIRST_PREDICTION_TARGET,
        "optional_prediction_target": OPTIONAL_PREDICTION_TARGET,
        "equal_proper_time_method": EQUAL_TAU_METHOD,
        "taylor_indicator_role": TAYLOR_ROLE,
        "hard_acceptance_gate": HARD_ACCEPTANCE_GATE,
        "single_process": True,
        "missing": list(MISSING),
        "later_outputs": list(LATER_OUTPUTS),
    }


def _validate_contrast(contrast):
    array = np.asarray(contrast, dtype=float)
    if array.shape != (6,) or not np.isfinite(array).all():
        raise ValueError("occupation contrast must be six finite weights")
    if not np.array_equal(array, OCCUPATION_CONTRAST):
        raise ValueError("occupation contrast is the locked trace-free (+, +, 0, 0, -, -)")
    if float(np.sum(array)) != 0.0:
        raise ValueError("occupation contrast must preserve the occupation trace")
    return array


def assert_admissible_source(phi0, phi1, weights):
    """Six weights in ``(0, 1]`` and covariance eigenvalues in ``[0, 1]``."""
    values = np.asarray(weights, dtype=float)
    if values.shape != (6,) or not np.isfinite(values).all():
        raise ValueError("source weights must be six finite numbers")
    if np.any(values <= 0.0) or np.any(values > 1.0):
        raise ValueError("source weights must lie in (0, 1]")
    eigenvalues = np.asarray(galerkin.occupation_eigenvalues(phi0, phi1, values), dtype=float)
    if not np.isfinite(eigenvalues).all() or np.any(eigenvalues < -1e-8) or np.any(eigenvalues > 1.0 + 1e-8):
        raise ValueError("source eigenvalues left the admissible interval [0, 1]")
    return eigenvalues


def pair_with_occupations(pair, weights):
    """Copy one pair onto a new occupation vector. The caller's pair is unchanged."""
    values = np.asarray(weights, dtype=float).copy()
    if values.shape != (6,) or not np.isfinite(values).all() or np.any(values <= 0.0) or np.any(values > 1.0):
        raise ValueError("replacement occupations must lie in (0, 1]")
    grid = replace(pair.grid, fine=replace(pair.grid.fine, occupations=values.copy()))
    return replace(pair, grid=grid, weights=values.copy())


def _saved_documents():
    v1 = json.loads(episode.V1_JSON.read_text())
    basis = json.loads(episode.BASIS_JSON.read_text())
    payload = episode.file_sha256(episode.V1_NPZ)
    basis_payload = episode.file_sha256(episode.BASIS_NPZ)
    if payload != v1["payload_sha256"]:
        raise ValueError("saved separated-pair payload does not match its record")
    if basis_payload != basis["payload_sha256"]:
        raise ValueError("frozen replay basis payload does not match its record")
    return v1, basis


def _frame_index(times, time):
    time = float(time)
    if not np.isfinite(time) or time < -1e-15 or time > HANDOFF_TIME + 1e-15:
        raise ValueError("saved baseline frames lie on [0, 0.3]")
    matches = np.flatnonzero(np.abs(times - time) <= 1e-12)
    if matches.size != 1:
        raise ValueError("requested time is not a saved baseline frame")
    return int(matches[0])


def load_saved_baseline(case_name, time=0.0):
    """Saved Cauchy frame and frozen ``W``. No initial solve and no new source.

    ``time=0`` is the prediction seed. ``time=0.3`` is the comparison frame.
    Loading the comparison frame does not change its occupations or radius.
    """
    if case_name not in episode.BASELINE_CASES:
        raise ValueError("baseline cases are the saved nf128 and nf256 records")
    v1, basis = _saved_documents()
    case = v1["results"][case_name]
    if abs(float(case["summary"]["final"]["time"]) - HANDOFF_TIME) > 1e-12:
        raise ValueError("saved baseline does not end at T=0.3")
    nf = int(case["nf"])
    entry = basis["bases"][str(nf)]
    with np.load(episode.V1_NPZ, allow_pickle=False) as episode_payload, np.load(
        episode.BASIS_NPZ, allow_pickle=False,
    ) as frozen:
        times = np.array(episode_payload[case_name + "_times"], dtype=float, copy=True)
        index = _frame_index(times, time)
        fields = {
            name: np.array(episode_payload[case_name + "_" + name][index], copy=True)
            for name in _STATE_FIELDS
        }
        initial_fields = {
            name: np.array(episode_payload[case_name + "_" + name][0], copy=True)
            for name in _STATE_FIELDS
        }
        final_fields = {
            name: np.array(episode_payload[case_name + "_" + name][-1], copy=True)
            for name in _STATE_FIELDS
        }
        observer = np.array(episode_payload[case_name + "_reference_columns"], copy=True)
        basis_w = np.array(frozen[f"nf{nf}_W"], copy=True)
        basis_observer = np.array(frozen[f"nf{nf}_reference_columns"], copy=True)
        basis_source = np.array(frozen[f"nf{nf}_source_columns"], copy=True)
    if abs(float(times[0])) > 1e-12 or abs(float(times[-1]) - HANDOFF_TIME) > 1e-12:
        raise ValueError("saved baseline time array does not run from 0 to 0.3")
    if not np.array_equal(observer, basis_observer):
        raise ValueError("saved observer and frozen basis observer differ")
    initial = model.NestedState(*(initial_fields[name] for name in _STATE_FIELDS))
    final = model.NestedState(*(final_fields[name] for name in _STATE_FIELDS))
    state = model.NestedState(*(fields[name] for name in _STATE_FIELDS))
    summary = case["summary"]
    if episode.state_sha256(initial) != summary["initial_state_sha256"]:
        raise ValueError("saved baseline initial state does not match its pin")
    if episode.state_sha256(final) != summary["final_state_sha256"]:
        raise ValueError("saved baseline final state does not match its pin")
    weights = np.array(case["source"]["occupations"], dtype=float)
    if weights.shape != (6,) or not np.array_equal(weights, np.asarray(coupling.OCCUPATIONS, dtype=float)):
        raise ValueError("saved baseline occupations are not the locked six weights")
    assert_admissible_source(state.phi0, state.phi1, weights)
    record = {
        "nf": nf,
        "coarse_indices": entry["coarse_indices"],
        "child_indices": entry["child_indices"],
        "parent_indices": entry["parent_indices"],
        "source_metadata": entry["source"],
        "geometry_metadata": entry["geometry"],
        "clock_locations": list(summary["final"]["metrics"]["clock_locations"]),
    }
    arrays = {
        "W": basis_w,
        "source_phi0": basis_source[:nf],
        "source_phi1": basis_source[nf:],
        "observer_columns": basis_observer,
        "source_weights": weights,
    }
    pair = episode.pair_from_arrays(arrays, record)
    if not np.array_equal(np.asarray(pair.geometry_map), basis_w):
        raise ValueError("pair reconstruction changed frozen W")
    if not np.array_equal(np.asarray(pair.grid.fine.occupations), weights):
        raise ValueError("pair reconstruction changed the saved occupations")
    resolving = abs(float(time) - HANDOFF_TIME) <= 1e-12
    return {
        "case_name": case_name,
        "time": float(times[index]),
        "frame_index": index,
        "nf": nf,
        "step_cap": float(case["step_cap"]),
        "state": state,
        "pair": pair,
        "weights": weights.copy(),
        "times": times,
        "clock_locations": record["clock_locations"],
        "state_sha256": episode.state_sha256(state),
        "initial_state_sha256": summary["initial_state_sha256"],
        "final_state_sha256": summary["final_state_sha256"],
        "pins": {
            "v1_json_sha256": episode.file_sha256(episode.V1_JSON),
            "v1_npz_sha256": v1["payload_sha256"],
            "basis_json_sha256": episode.file_sha256(episode.BASIS_JSON),
            "basis_npz_sha256": basis["payload_sha256"],
            "state_sha256": episode.state_sha256(state),
            "initial_state_sha256": summary["initial_state_sha256"],
            "final_state_sha256": summary["final_state_sha256"],
            "W_sha256": episode.array_sha256(basis_w),
            "weights_sha256": episode.array_sha256(weights),
            "source_columns_sha256": episode.array_sha256(basis_source),
            "observer_columns_sha256": episode.array_sha256(basis_observer),
        },
        "frozen_W": True,
        "initial_state_called": False,
        "constraint_resolved": False,
        "source_changed": False,
        "comparison_frame": resolving,
        "domain": _domain(nf, weights),
    }


def compare_saved_baseline(state, case_name, time):
    """Difference from a saved frame. The frame's source is left as stored."""
    frame = load_saved_baseline(case_name, time)
    gaps = {}
    for name in _STATE_FIELDS:
        difference = np.asarray(getattr(state, name)) - np.asarray(getattr(frame["state"], name))
        gaps[name] = float(np.max(np.abs(difference)))
    return {
        "case_name": case_name,
        "time": frame["time"],
        "gaps": gaps,
        "max_gap": max(gaps.values()),
        "saved_state_sha256": frame["state_sha256"],
        "weights": frame["weights"].tolist(),
        "weights_sha256": frame["pins"]["weights_sha256"],
        "source_changed": False,
        "constraint_resolved": False,
        "initial_state_called": False,
        "pins": frame["pins"],
    }


def _domain(nf, weights):
    return {
        "gauge": "conformal",
        "nf": int(nf),
        "child_interval": list(CHILD_INTERVAL),
        "child_clock": CHILD_CLOCK,
        "parent_interval": list(model.PARENT_INTERVAL),
        "occupation_layout": "separated rank-6",
        "occupations": [float(value) for value in np.asarray(weights, dtype=float)],
        "occupation_contrast": OCCUPATION_CONTRAST.tolist(),
        "contrast_packets": list(CONTRAST_PACKETS),
    }


def _identity(pair, state, weights):
    eigenvalues = assert_admissible_source(state.phi0, state.phi1, weights)
    return {
        "state_sha256": episode.state_sha256(state),
        "W_sha256": episode.array_sha256(np.asarray(pair.geometry_map)),
        "weights_sha256": episode.array_sha256(np.asarray(weights, dtype=float)),
        "eigenvalue_min": float(np.min(eigenvalues)),
        "eigenvalue_max": float(np.max(eigenvalues)),
        "eigenvalues_admissible": True,
        "nf": int(pair.grid.nf),
        "backend": backend.operator_backend(pair.grid),
    }


def contrast_tangent(pair, state, contrast=None):
    """Trace-free occupation contrast plus one linearized radius solve.

    ``Q``, ``chi``, the momenta and the column tangent stay zero. No Newton
    loop is run and the base radius is not replaced.
    """
    if backend.operator_backend(pair.grid) != PREPARATION_BACKEND:
        raise ValueError("the radius tangent is formed on the dense preparation grid")
    direction = _validate_contrast(OCCUPATION_CONTRAST if contrast is None else contrast)
    tangent = discovery.zero_tangent(pair.grid)
    tangent.occupations[:] = direction
    prepared = discovery.prepared_nested_radius_tangent(pair, state, tangent)
    if not prepared.available or prepared.delta_r is None:
        missing = prepared.missing_primitive or "unavailable_radius_tangent"
        raise ValueError("prepared radius tangent is missing: " + missing)
    if prepared.newton_used:
        raise ValueError("analytic radius tangent does not use Newton")
    tangent.r = np.array(prepared.delta_r, dtype=float, copy=True)
    return tangent, prepared


def _allowed_source_amplitude(amplitude):
    allowed = tuple(DERIVATIVE_H) + tuple(-step for step in DERIVATIVE_H) + (HELD_OUT_ALPHA,)
    if not any(abs(amplitude - item) <= 1e-15 for item in allowed):
        raise ValueError("source amplitudes are ± the derivative steps or the held-out alpha")


def prepare_changed_source(pair, alpha):
    """Dense source-consistent initial state at ``c + alpha * contrast``.

    The baseline itself is not sent back through the constructor. A zero
    amplitude is refused so a saved or already solved state is not replaced.
    """
    if backend.operator_backend(pair.grid) != PREPARATION_BACKEND:
        raise ValueError("changed sources are prepared with the dense original constructor")
    amplitude = float(alpha)
    if not np.isfinite(amplitude) or amplitude == 0.0:
        raise ValueError("a changed source needs a nonzero finite amplitude")
    _allowed_source_amplitude(amplitude)
    base = np.array(pair.grid.fine.occupations, dtype=float, copy=True)
    weights = base + amplitude * _validate_contrast(OCCUPATION_CONTRAST)
    if abs(float(np.sum(weights) - np.sum(base))) > 1e-12:
        raise ValueError("changed source did not preserve the occupation trace")
    assert_admissible_source(pair.source_phi0, pair.source_phi1, weights)
    changed = pair_with_occupations(pair, weights)
    state, report = model.initial_state(changed, solve_constraints=True)
    if not report.get("converged"):
        raise ValueError("dense source-consistent constructor did not converge")
    eigenvalues = assert_admissible_source(state.phi0, state.phi1, weights)
    if not np.array_equal(np.asarray(pair.grid.fine.occupations), base):
        raise ValueError("changed-source preparation mutated the caller's occupations")
    return {
        "pair": changed,
        "state": state,
        "alpha": amplitude,
        "weights": weights,
        "constructor": CONSTRUCTOR,
        "constructor_grid": PREPARATION_BACKEND,
        "method": report.get("method"),
        "converged": True,
        "constraint_resolved": True,
        "newton_homotopy": True,
        "baseline_reused": False,
        "eigenvalues_admissible": True,
        "eigenvalue_min": float(np.min(eigenvalues)),
        "eigenvalue_max": float(np.max(eigenvalues)),
    }


def _resolve(pair, state, backend_name):
    resolved, returned, info = episode.resolve_pair(pair, state, backend=backend_name)
    if info.get("W_changed") or info.get("state_changed") or info.get("initial_state_called"):
        raise ValueError("carrier resolution changed W, the state, or called the initializer")
    if backend_name == "fft" and info.get("backend") != "fft":
        raise ValueError("FFT evolution was requested and not selected")
    if not all(np.array_equal(getattr(returned, name), getattr(state, name)) for name in _STATE_FIELDS):
        raise ValueError("carrier resolution changed the Cauchy state")
    return resolved, returned, info


def _shift_tangent(tangent, rate, factor):
    factor = float(factor)
    return discovery.StateTangent(
        tangent.Q + factor * rate.Q,
        tangent.r + factor * rate.r,
        tangent.chi + factor * rate.chi,
        tangent.p_Q + factor * rate.p_Q,
        tangent.p_r + factor * rate.p_r,
        tangent.p_chi + factor * rate.p_chi,
        tangent.phi0 + factor * rate.phi0,
        tangent.phi1 + factor * rate.phi1,
        np.array(tangent.occupations, dtype=float, copy=True),
    )


def _nodal_tangent(pair, tangent):
    nodal = model.reconstruct_state(pair, model.NestedState(
        tangent.Q, tangent.r, tangent.chi, tangent.p_Q, tangent.p_r, tangent.p_chi,
        tangent.phi0, tangent.phi1,
    ))
    return discovery.StateTangent(
        nodal.Q, nodal.r, nodal.chi, nodal.p_Q, nodal.p_r, nodal.p_chi,
        nodal.phi0, nodal.phi1, np.array(tangent.occupations, dtype=float, copy=True),
    )


def clock_rates(pair, state, tangent):
    """Instantaneous ``τ̇`` and ``δτ̇`` at one slice. The increment is not formed here."""
    nodal_state = model.reconstruct_state(pair, state)
    nodal_tangent = _nodal_tangent(pair, tangent)
    tau_dot, delta_tau_dot, frozen = discovery.child_clock_rate(
        pair.grid, nodal_state, nodal_tangent, coordinate_duration=0.0,
    )
    if frozen != 0.0:
        raise ValueError("a zero coordinate duration must not emit a frozen clock increment")
    return float(tau_dot), float(delta_tau_dot)


def _rk_quadrature(samples):
    first, second, third, fourth = (float(value) for value in samples)
    return (first + 2.0 * second + 2.0 * third + fourth) / 6.0


def _chart(pair, state):
    fine = galerkin.prolong_state(pair.grid, model.reconstruct_state(pair, state))
    active = galerkin.active_fine_system(pair.grid, fine)
    reason = coupling.chart_failure(active, fine)
    if reason is not None:
        raise coupling.PositiveChartExit(reason, np.nan, state.copy())


def coupled_rk4_step(pair, state, tangent, dt):
    """One RK4 step of the state and of ``DJ · tangent``, with stage clock sums."""
    dt = float(dt)
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("a positive finite timestep is required")
    with backend.fft_thread_limit(1):
        rate1 = model.rates(pair, state)
        jvp1 = discovery.nested_rate_jacobian_vector(pair, state, tangent)
        state2 = model._combine(state, rate1, dt / 2.0)
        tangent2 = _shift_tangent(tangent, jvp1, dt / 2.0)
        rate2 = model.rates(pair, state2)
        jvp2 = discovery.nested_rate_jacobian_vector(pair, state2, tangent2)
        state3 = model._combine(state, rate2, dt / 2.0)
        tangent3 = _shift_tangent(tangent, jvp2, dt / 2.0)
        rate3 = model.rates(pair, state3)
        jvp3 = discovery.nested_rate_jacobian_vector(pair, state3, tangent3)
        state4 = model._combine(state, rate3, dt)
        tangent4 = _shift_tangent(tangent, jvp3, dt)
        rate4 = model.rates(pair, state4)
        jvp4 = discovery.nested_rate_jacobian_vector(pair, state4, tangent4)
        updated = model.NestedState(*(
            getattr(state, name) + dt / 6.0 * (
                getattr(rate1, name) + 2.0 * getattr(rate2, name)
                + 2.0 * getattr(rate3, name) + getattr(rate4, name)
            )
            for name in _STATE_FIELDS
        ))
        updated_tangent = discovery.StateTangent(*(
            getattr(tangent, name) + dt / 6.0 * (
                getattr(jvp1, name) + 2.0 * getattr(jvp2, name)
                + 2.0 * getattr(jvp3, name) + getattr(jvp4, name)
            )
            for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1")
        ), np.array(tangent.occupations, dtype=float, copy=True))
        tau_samples = []
        delta_samples = []
        for stage_state, stage_tangent in (
            (state, tangent), (state2, tangent2), (state3, tangent3), (state4, tangent4),
        ):
            tau_dot, delta_tau_dot = clock_rates(pair, stage_state, stage_tangent)
            tau_samples.append(tau_dot)
            delta_samples.append(delta_tau_dot)
    _chart(pair, updated)
    return {
        "state": updated,
        "tangent": updated_tangent,
        "tau": dt * _rk_quadrature(tau_samples),
        "delta_tau": dt * _rk_quadrature(delta_samples),
        "stage_tau_dot": tuple(tau_samples),
        "stage_delta_tau_dot": tuple(delta_samples),
        "dt": dt,
        "frozen_rate_interval": False,
        "evolved_clock_integral": True,
    }


def nonlinear_rk4_step(pair, state, dt):
    """State RK4 with ``τ`` summed on the stages. No tangent is formed."""
    dt = float(dt)
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("a positive finite timestep is required")
    zero = discovery.zero_tangent(pair.grid)
    with backend.fft_thread_limit(1):
        rate1 = model.rates(pair, state)
        state2 = model._combine(state, rate1, dt / 2.0)
        rate2 = model.rates(pair, state2)
        state3 = model._combine(state, rate2, dt / 2.0)
        rate3 = model.rates(pair, state3)
        state4 = model._combine(state, rate3, dt)
        rate4 = model.rates(pair, state4)
        updated = model.NestedState(*(
            getattr(state, name) + dt / 6.0 * (
                getattr(rate1, name) + 2.0 * getattr(rate2, name)
                + 2.0 * getattr(rate3, name) + getattr(rate4, name)
            )
            for name in _STATE_FIELDS
        ))
        tau_samples = [
            clock_rates(pair, stage, zero)[0]
            for stage in (state, state2, state3, state4)
        ]
    _chart(pair, updated)
    return {
        "state": updated,
        "tau": dt * _rk_quadrature(tau_samples),
        "delta_tau": 0.0,
        "stage_tau_dot": tuple(tau_samples),
        "dt": dt,
        "frozen_rate_interval": False,
        "evolved_clock_integral": True,
    }


def _assert_duration(duration, production):
    duration = float(duration)
    if not np.isfinite(duration) or duration <= 0.0:
        raise ValueError("coordinate duration must be positive and finite")
    if not production and duration > NONPRODUCTION_DURATION + 1e-15:
        raise PermissionError(
            "coordinate duration above 0.01 is a production trajectory; "
            "the campaign executor launches it after freeze"
        )
    return duration


def _choose_step(limit, remaining, fixed_step):
    if fixed_step is None:
        step = min(float(limit), remaining)
    else:
        step = min(float(fixed_step), remaining)
        if step > float(limit) + 1e-15:
            raise ValueError("fixed step exceeds the RK4 stability restriction")
    if step <= 0.0:
        raise ValueError("timestep vanished before the duration was reached")
    return step


def _advance(pair, state, tangent, *, duration, step_cap, fixed_step, production, backend_name, coupled):
    duration = _assert_duration(duration, production)
    step_cap = float(step_cap)
    if not np.isfinite(step_cap) or step_cap <= 0.0:
        raise ValueError("a positive finite step cap is required")
    resolved, state, info = _resolve(pair, state, backend_name)
    weights = np.array(resolved.grid.fine.occupations, dtype=float, copy=True)
    if tangent is None:
        tangent = discovery.zero_tangent(resolved.grid)
    tau_dot0, delta_tau_dot0 = clock_rates(resolved, state, tangent)
    tau = 0.0
    delta_tau = 0.0
    mark = 0.0
    taken = []
    current = state
    current_tangent = tangent
    while mark < duration - 1e-13:
        limit, _omega, _quadrature = model.stable_timestep(resolved, current, step_cap)
        step = _choose_step(limit, duration - mark, fixed_step)
        if coupled:
            stepped = coupled_rk4_step(resolved, current, current_tangent, step)
            current_tangent = stepped["tangent"]
        else:
            stepped = nonlinear_rk4_step(resolved, current, step)
        current = stepped["state"]
        tau += float(stepped["tau"])
        delta_tau += float(stepped["delta_tau"])
        mark += step
        taken.append(step)
    if not taken:
        raise ValueError("propagation took no steps")
    eigenvalues = assert_admissible_source(current.phi0, current.phi1, weights)
    readout = matched_readout(
        resolved, current, current_tangent, tau=tau, delta_tau=delta_tau, coordinate_time=mark,
    )
    return {
        "state": current,
        "tangent": current_tangent if coupled else None,
        "pair": resolved,
        "weights": weights,
        "coordinate_time": mark,
        "steps": taken,
        "step_count": len(taken),
        "step_cap": step_cap,
        "fixed_step": None if fixed_step is None else float(fixed_step),
        "tau": tau,
        "delta_tau": delta_tau,
        "tau_dot_initial": tau_dot0,
        "delta_tau_dot_initial": delta_tau_dot0,
        "frozen_slice_tau": tau_dot0 * duration,
        "frozen_slice_delta_tau": delta_tau_dot0 * duration,
        "frozen_rate_interval": False,
        "evolved_clock_integral": True,
        "clock_quadrature": CLOCK_QUADRATURE,
        "backend": info.get("backend"),
        "initial_state_called": False,
        "constraint_resolved": False,
        "readout": readout,
        "eigenvalues_admissible": True,
        "eigenvalue_min": float(np.min(eigenvalues)),
        "eigenvalue_max": float(np.max(eigenvalues)),
        "identity": _identity(resolved, current, weights),
        "domain": _domain(resolved.grid.nf, weights),
    }


def matched_readout(pair, state, tangent, *, tau, delta_tau, coordinate_time):
    """``δO|τ = δO|t − Ȯ δτ / τ̇`` with the integrated clock, at the final slice."""
    raw = discovery.nested_readout(pair, state, tangent, coordinate_duration=0.0)
    if raw["evolved_clock_integral"] is not False or float(raw["delta_tau"]) != 0.0:
        raise ValueError("the slice primitive must not supply the propagated clock")
    primary = discovery.matched_tau_correction(
        raw["delta_child_regional_content_t"],
        raw["child_regional_content_dot"],
        delta_tau,
        raw["tau_dot"],
    )
    secondary = discovery.matched_tau_correction(
        raw["delta_child_proper_mean_r_t"],
        raw["child_proper_mean_r_dot"],
        delta_tau,
        raw["tau_dot"],
    )
    return {
        "primary_readout": PRIMARY_READOUT,
        "secondary_readout": SECONDARY_READOUT,
        "formula": FORMULA,
        "child_regional_content": raw["child_regional_content"],
        "delta_child_regional_content_t": raw["delta_child_regional_content_t"],
        "child_regional_content_dot": raw["child_regional_content_dot"],
        "delta_child_regional_content_tau": primary,
        "child_proper_mean_r": raw["child_proper_mean_r"],
        "delta_child_proper_mean_r_t": raw["delta_child_proper_mean_r_t"],
        "child_proper_mean_r_dot": raw["child_proper_mean_r_dot"],
        "delta_child_proper_mean_r_tau": secondary,
        "tau": float(tau),
        "tau_dot": raw["tau_dot"],
        "delta_tau": float(delta_tau),
        "delta_tau_dot": raw["delta_tau_dot"],
        "coordinate_time": float(coordinate_time),
        "slice_coordinate_duration": 0.0,
        "frozen_rate_interval": False,
        "evolved_clock_integral": True,
        "clock_quadrature": CLOCK_QUADRATURE,
    }


def propagate_coupled(pair, state, tangent, *, duration, step_cap, production=False,
                      fixed_step=None, backend_name=EVOLUTION_BACKEND):
    """FFT (or an explicit dense oracle) advance of one full state and its tangent."""
    return _advance(
        pair, state, tangent, duration=duration, step_cap=step_cap, fixed_step=fixed_step,
        production=production, backend_name=backend_name, coupled=True,
    )


def propagate_baseline(pair, state, *, duration, step_cap, production=False,
                       fixed_step=None, backend_name=EVOLUTION_BACKEND):
    """Nonlinear advance used by the saved-frame comparison and by finite differences."""
    return _advance(
        pair, state, None, duration=duration, step_cap=step_cap, fixed_step=fixed_step,
        production=production, backend_name=backend_name, coupled=False,
    )


def predict_occupation_response(pair, state, *, duration, step_cap, production=False,
                                fixed_step=None, backend_name=EVOLUTION_BACKEND):
    """Analytic contrast response. The initializer is not called."""
    tangent, prepared = contrast_tangent(pair, state)
    propagated = propagate_coupled(
        pair, state, tangent, duration=duration, step_cap=step_cap, production=production,
        fixed_step=fixed_step, backend_name=backend_name,
    )
    readout = propagated["readout"]
    return {
        "schema": SCHEMA,
        "method": METHOD,
        "primary_readout": PRIMARY_READOUT,
        "secondary_readout": SECONDARY_READOUT,
        "primary_tau": readout["delta_child_regional_content_tau"],
        "secondary_tau": readout["delta_child_proper_mean_r_tau"],
        "child_regional_content": readout["child_regional_content"],
        "child_proper_mean_r": readout["child_proper_mean_r"],
        "formula": FORMULA,
        "radius_tangent": RADIUS_TANGENT,
        "jacobian_owner": prepared.jacobian_owner,
        "linear_residual_max": prepared.linear_residual_max,
        "newton_used": False,
        "initial_state_called": False,
        "constraint_resolved_on_baseline": False,
        "t03_source_changed": False,
        "t03_constraint_resolved": False,
        "effective_stress": None,
        "finite_difference_used_as_analytic": False,
        "propagated": propagated,
        "readout": readout,
        "eigenvalues_admissible": propagated["eigenvalues_admissible"],
        "evolved_clock_integral": True,
        "frozen_rate_interval": False,
        "domain": propagated["domain"],
        "identity": _identity(pair, state, pair.grid.fine.occupations),
    }


def _independent_matched(baseline, evolved):
    """Taylor indicator at one shared coordinate time.

    This is the linear clock correction. It is not a measurement at equal
    proper time, and it uses the same formula as the analytic prediction.
    """
    content_dot = baseline["readout"]["child_regional_content_dot"]
    mean_dot = baseline["readout"]["child_proper_mean_r_dot"]
    tau_dot = baseline["readout"]["tau_dot"]
    delta_tau = evolved["tau"] - baseline["tau"]
    delta_content = evolved["readout"]["child_regional_content"] - baseline["readout"]["child_regional_content"]
    delta_mean = evolved["readout"]["child_proper_mean_r"] - baseline["readout"]["child_proper_mean_r"]
    return {
        "primary": discovery.matched_tau_correction(delta_content, content_dot, delta_tau, tau_dot),
        "secondary": discovery.matched_tau_correction(delta_mean, mean_dot, delta_tau, tau_dot),
        "delta_content_t": delta_content,
        "delta_mean_t": delta_mean,
        "delta_tau": delta_tau,
        "matched_content": baseline["readout"]["child_regional_content"] + discovery.matched_tau_correction(
            delta_content, content_dot, delta_tau, tau_dot,
        ),
    }


def finite_difference_control(pair, state, *, duration, step_cap, production=False, fixed_step=None):
    """Centred differences at both derivative steps, plus the separate held-out run.

    Changed sources are built by the dense constructor and evolved with FFT.
    The analytic number is not replaced by either difference.
    """
    if abs(HELD_OUT_ALPHA) in {abs(value) for value in DERIVATIVE_H}:
        raise ValueError("held-out alpha collided with a derivative step")
    analytic = predict_occupation_response(
        pair, state, duration=duration, step_cap=step_cap, production=production,
        fixed_step=fixed_step, backend_name=EVOLUTION_BACKEND,
    )
    baseline = propagate_baseline(
        pair, state, duration=duration, step_cap=step_cap, production=production,
        fixed_step=fixed_step, backend_name=EVOLUTION_BACKEND,
    )
    state_gap = max(
        float(np.max(np.abs(getattr(analytic["propagated"]["state"], name) - getattr(baseline["state"], name))))
        for name in _STATE_FIELDS
    )
    amplitudes = tuple(sign * step for step in DERIVATIVE_H for sign in (1.0, -1.0)) + (HELD_OUT_ALPHA,)
    arms = []
    for amplitude in amplitudes:
        prepared = prepare_changed_source(pair, amplitude)
        evolved = propagate_baseline(
            prepared["pair"], prepared["state"], duration=duration, step_cap=step_cap,
            production=production, fixed_step=fixed_step, backend_name=EVOLUTION_BACKEND,
        )
        arms.append({
            "alpha": float(amplitude),
            "evolved": evolved,
            "constructor": prepared["constructor"],
            "constructor_grid": prepared["constructor_grid"],
            "method": prepared["method"],
            "converged": prepared["converged"],
            "constraint_resolved": True,
            "eigenvalues_admissible": prepared["eigenvalues_admissible"],
            "eigenvalue_min": prepared["eigenvalue_min"],
            "eigenvalue_max": prepared["eigenvalue_max"],
            "weights": prepared["weights"].tolist(),
            "final_eigenvalues_admissible": evolved["eigenvalues_admissible"],
        })

    def arm(amplitude):
        matches = [item for item in arms if abs(item["alpha"] - amplitude) <= 1e-15]
        if len(matches) != 1:
            raise ValueError("missing finite-difference arm")
        return matches[0]["evolved"]

    centred = []
    for step in DERIVATIVE_H:
        plus = _independent_matched(baseline, arm(step))
        minus = _independent_matched(baseline, arm(-step))
        centred.append({
            "h": float(step),
            "primary": (plus["primary"] - minus["primary"]) / (2.0 * step),
            "secondary": (plus["secondary"] - minus["secondary"]) / (2.0 * step),
            "delta_tau": (plus["delta_tau"] - minus["delta_tau"]) / (2.0 * step),
            "used_as_analytic": False,
            "comparison_kind": "linear_clock_correction",
            "equal_proper_time": False,
            "same_formula_as_prediction": True,
        })
    held = _independent_matched(baseline, arm(HELD_OUT_ALPHA))
    scale = HELD_OUT_ALPHA
    return {
        "schema": SCHEMA,
        "finite_difference_used_as_analytic": False,
        "evolution_backend": EVOLUTION_BACKEND,
        "preparation_backend": PREPARATION_BACKEND,
        "derivative_h": list(DERIVATIVE_H),
        "held_out_alpha": HELD_OUT_ALPHA,
        "held_out_distinct_from_derivative_h": True,
        "analytic": {
            "primary_tau": analytic["primary_tau"],
            "secondary_tau": analytic["secondary_tau"],
            "delta_tau": analytic["readout"]["delta_tau"],
            "frozen_slice_delta_tau": analytic["propagated"]["frozen_slice_delta_tau"],
            "child_regional_content": baseline["readout"]["child_regional_content"],
            "child_proper_mean_r": baseline["readout"]["child_proper_mean_r"],
            "initial_state_called": False,
            "newton_used": False,
            "eigenvalues_admissible": analytic["eigenvalues_admissible"],
            "step_count": baseline["step_count"],
            "coordinate_time": baseline["coordinate_time"],
            "nonlinear_state_gap": state_gap,
        },
        "centred": centred,
        "held_out": {
            "alpha": HELD_OUT_ALPHA,
            "used_as_derivative_h": False,
            "comparison_kind": "linear_clock_correction",
            "equal_proper_time": False,
            "same_formula_as_prediction": True,
            "primary_scaled": held["primary"] / scale,
            "secondary_scaled": held["secondary"] / scale,
            "predicted_content": baseline["readout"]["child_regional_content"] + scale * analytic["primary_tau"],
            "measured_matched_content": held["matched_content"],
            "eigenvalues_admissible": arm(HELD_OUT_ALPHA)["eigenvalues_admissible"],
        },
        "preparations": [
            {key: value for key, value in item.items() if key != "evolved"}
            for item in arms
        ],
        "effective_stress": None,
        "domain": analytic["domain"],
    }


def accepted_duration(duration, *, production):
    """Short tests, the first target 0.3, or optional 1. Not a residual gate."""
    duration = float(duration)
    if not np.isfinite(duration) or duration <= 0.0:
        raise ValueError("duration must be positive and finite")
    if duration <= NONPRODUCTION_DURATION + 1e-15:
        return {"duration": duration, "production": False, "target": "short_test"}
    first = abs(duration - FIRST_PREDICTION_TARGET) <= 1e-15
    optional = abs(duration - OPTIONAL_PREDICTION_TARGET) <= 1e-15
    if first or optional:
        if not production:
            raise PermissionError(
                "duration above 0.01 is a production trajectory; "
                "the campaign executor launches it"
            )
        return {
            "duration": duration,
            "production": True,
            "target": "T0.3" if first else "optional_T1",
        }
    raise ValueError(
        "duration must be a short test at most 0.01, the first prediction target 0.3, "
        "or optional 1; late T=3 curvature is unresolved and is not a target"
    )


def _observables(pair, state):
    nodal = model.reconstruct_state(pair, state)
    return {
        "child_regional_content": float(discovery.child_regional_content(pair.grid, nodal)),
        "child_proper_mean_r": float(discovery.child_proper_mean_r(pair.grid, nodal)),
        "tau_dot": float(discovery.child_clock_rate(pair.grid, nodal)),
    }


def _append_sample(history, pair, state, coordinate_time, tau):
    observed = _observables(pair, state)
    history["time"].append(float(coordinate_time))
    history["tau"].append(float(tau))
    history["child_regional_content"].append(observed["child_regional_content"])
    history["child_proper_mean_r"].append(observed["child_proper_mean_r"])
    history["tau_dot"].append(observed["tau_dot"])
    return observed


def _march_clock(pair, state, *, duration, step_cap, fixed_step, production, tau_target=None):
    """March to a coordinate duration and, when asked, through a proper-time crossing."""
    duration = _assert_duration(duration, production)
    step_cap = float(step_cap)
    resolved, state, info = _resolve(pair, state, EVOLUTION_BACKEND)
    if info.get("backend") != EVOLUTION_BACKEND:
        raise ValueError("equal-time comparison evolves on the FFT carrier")
    weights = np.array(resolved.grid.fine.occupations, dtype=float, copy=True)
    current = state.copy()
    tau = 0.0
    mark = 0.0
    history = {"time": [], "tau": [], "child_regional_content": [], "child_proper_mean_r": [], "tau_dot": []}
    _append_sample(history, resolved, current, mark, tau)
    bracket = None
    coordinate = None
    steps = []
    ceiling = duration * 1.25
    if not production:
        ceiling = min(ceiling, NONPRODUCTION_DURATION)
    while True:
        reached_time = mark >= duration - 1e-13
        reached_tau = tau_target is None or tau >= float(tau_target)
        if coordinate is None and reached_time:
            coordinate = {
                "state": current.copy(),
                "tau": tau,
                "time": mark,
                "observed": _observables(resolved, current),
            }
        if reached_time and reached_tau:
            break
        if mark >= ceiling - 1e-15:
            break
        limit, _omega, _quadrature = model.stable_timestep(resolved, current, step_cap)
        horizon = duration if mark < duration - 1e-13 else ceiling
        step = _choose_step(limit, horizon - mark, fixed_step)
        previous = current.copy()
        previous_tau = tau
        previous_mark = mark
        stepped = nonlinear_rk4_step(resolved, current, step)
        new_tau = tau + float(stepped["tau"])
        if tau_target is not None and bracket is None and previous_tau < float(tau_target) <= new_tau:
            bracket = {
                "state": previous,
                "tau": previous_tau,
                "time": previous_mark,
                "dt_hi": step,
                "tau_hi": new_tau,
            }
        current = stepped["state"]
        tau = new_tau
        mark += step
        steps.append(step)
        _append_sample(history, resolved, current, mark, tau)
    if not steps:
        raise ValueError("clock march took no steps")
    return {
        "pair": resolved,
        "weights": weights,
        "state": current.copy(),
        "tau": tau,
        "coordinate_time": mark,
        "steps": steps,
        "history": history,
        "bracket": bracket,
        "coordinate": coordinate,
        "backend": info.get("backend"),
        "ceiling": ceiling,
        "tau_target_reached": tau_target is None or tau >= float(tau_target),
    }


def _root_step_to_tau(pair, bracket, tau_target):
    """Decrease the crossing step until the proper-time bracket collapses."""
    target = float(tau_target)
    lo = 0.0
    hi = float(bracket["dt_hi"])
    tau_lo = float(bracket["tau"])
    best = None
    widths = []
    for _index in range(ROOT_BISECTIONS):
        widths.append(hi - lo)
        if hi - lo <= max(hi, bracket["dt_hi"]) * 1e-12:
            break
        mid = 0.5 * (lo + hi)
        if mid <= 0.0:
            break
        trial = nonlinear_rk4_step(pair, bracket["state"], mid)
        tau_mid = float(bracket["tau"]) + float(trial["tau"])
        if tau_mid + 1e-15 < tau_lo:
            raise ValueError("proper time decreased inside the equal-tau bracket")
        if tau_mid < target:
            lo = mid
            tau_lo = tau_mid
        else:
            hi = mid
            best = {
                "state": trial["state"],
                "tau": tau_mid,
                "dt": mid,
                "time": float(bracket["time"]) + mid,
            }
    if best is None:
        raise ValueError("equal-tau bracket did not produce a crossing step")
    observed = _observables(pair, best["state"])
    return {
        "state": best["state"],
        "tau": best["tau"],
        "time": best["time"],
        "dt": best["dt"],
        "tau_residual": best["tau"] - target,
        "bracket_width": float(widths[-1]) if widths else hi - lo,
        "bisections": len(widths),
        "method": EQUAL_TAU_METHOD,
        "observed": observed,
        "uses_linear_clock_correction": False,
    }


def _interpolate_at_tau(history, tau_target, key):
    tau = np.asarray(history["tau"], dtype=float)
    values = np.asarray(history[key], dtype=float)
    target = float(tau_target)
    if tau.size < 2 or np.any(np.diff(tau) <= 0.0):
        return {"available": False, "reason": "proper time samples are not strictly increasing"}
    if target < tau[0] or target > tau[-1]:
        return {"available": False, "reason": "proper time target is outside the stored samples"}
    upper = int(np.searchsorted(tau, target, side="left"))
    upper = min(max(upper, 1), tau.size - 1)
    lower = upper - 1
    span = float(tau[upper] - tau[lower])
    if span <= 0.0:
        return {"available": False, "reason": "proper time bracket has no width"}
    weight = (target - float(tau[lower])) / span
    return {
        "available": True,
        "value": float((1.0 - weight) * values[lower] + weight * values[upper]),
        "lower_time": float(history["time"][lower]),
        "upper_time": float(history["time"][upper]),
        "weight": weight,
    }


def _equal_tau_event(marched, tau_target):
    target = float(tau_target)
    if marched["bracket"] is not None:
        rooted = _root_step_to_tau(marched["pair"], marched["bracket"], target)
        interpolated = _interpolate_at_tau(marched["history"], target, "child_regional_content")
        mean_interp = _interpolate_at_tau(marched["history"], target, "child_proper_mean_r")
        content_gap = None
        if interpolated["available"]:
            content_gap = abs(interpolated["value"] - rooted["observed"]["child_regional_content"])
        return {
            "reached": True,
            "state": rooted["state"],
            "tau": rooted["tau"],
            "time": rooted["time"],
            "tau_residual": rooted["tau_residual"],
            "bracket_width": rooted["bracket_width"],
            "bisections": rooted["bisections"],
            "method": EQUAL_TAU_METHOD,
            "uses_linear_clock_correction": False,
            "child_regional_content": rooted["observed"]["child_regional_content"],
            "child_proper_mean_r": rooted["observed"]["child_proper_mean_r"],
            "tau_dot": rooted["observed"]["tau_dot"],
            "dense_interpolation": interpolated,
            "dense_mean_r_interpolation": mean_interp,
            "dense_interpolation_content_gap": content_gap,
        }
    if abs(marched["tau"] - target) <= 0.0 and marched["coordinate"] is not None:
        observed = marched["coordinate"]["observed"]
        return {
            "reached": True,
            "state": marched["coordinate"]["state"],
            "tau": marched["coordinate"]["tau"],
            "time": marched["coordinate"]["time"],
            "tau_residual": marched["coordinate"]["tau"] - target,
            "bracket_width": 0.0,
            "bisections": 0,
            "method": "coordinate_endpoint_already_at_target",
            "uses_linear_clock_correction": False,
            "child_regional_content": observed["child_regional_content"],
            "child_proper_mean_r": observed["child_proper_mean_r"],
            "tau_dot": observed["tau_dot"],
            "dense_interpolation": {"available": True, "value": observed["child_regional_content"]},
            "dense_mean_r_interpolation": {"available": True, "value": observed["child_proper_mean_r"]},
            "dense_interpolation_content_gap": 0.0,
        }
    return {
        "reached": False,
        "state": None,
        "tau": marched["tau"],
        "time": marched["coordinate_time"],
        "tau_residual": marched["tau"] - target,
        "bracket_width": None,
        "bisections": 0,
        "method": EQUAL_TAU_METHOD,
        "uses_linear_clock_correction": False,
        "child_regional_content": None,
        "child_proper_mean_r": None,
        "tau_dot": None,
        "dense_interpolation": {"available": False, "reason": "proper time target was not bracketed"},
        "dense_mean_r_interpolation": {"available": False, "reason": "proper time target was not bracketed"},
        "dense_interpolation_content_gap": None,
    }


def _coordinate_view(marched):
    """Slice at the requested coordinate duration, for the Taylor indicator only."""
    if marched["coordinate"] is None:
        raise ValueError("coordinate slice was not reached")
    coordinate = marched["coordinate"]
    return {
        "state": coordinate["state"],
        "tau": coordinate["tau"],
        "coordinate_time": coordinate["time"],
        "readout": {
            "child_regional_content": coordinate["observed"]["child_regional_content"],
            "child_proper_mean_r": coordinate["observed"]["child_proper_mean_r"],
            "child_regional_content_dot": None,
            "child_proper_mean_r_dot": None,
            "tau_dot": coordinate["observed"]["tau_dot"],
        },
    }


def _taylor_from_slices(baseline_slice, arm_slice, baseline_dots):
    """Linear correction using baseline rates. Not the equal-tau measurement."""
    baseline = {
        "tau": baseline_slice["tau"],
        "readout": {
            "child_regional_content": baseline_slice["readout"]["child_regional_content"],
            "child_proper_mean_r": baseline_slice["readout"]["child_proper_mean_r"],
            "child_regional_content_dot": baseline_dots["child_regional_content_dot"],
            "child_proper_mean_r_dot": baseline_dots["child_proper_mean_r_dot"],
            "tau_dot": baseline_dots["tau_dot"],
        },
    }
    evolved = {
        "tau": arm_slice["tau"],
        "readout": {
            "child_regional_content": arm_slice["readout"]["child_regional_content"],
            "child_proper_mean_r": arm_slice["readout"]["child_proper_mean_r"],
        },
    }
    matched = _independent_matched(baseline, evolved)
    return {
        "role": TAYLOR_ROLE,
        "comparison_kind": "linear_clock_correction",
        "equal_proper_time": False,
        "same_formula_as_prediction": True,
        "uses_linear_clock_correction": True,
        "primary": matched["primary"],
        "secondary": matched["secondary"],
        "delta_tau": matched["delta_tau"],
        "matched_content": matched["matched_content"],
    }


def _time_error(event, baseline_time):
    return {
        "tau_residual": None if event["tau_residual"] is None else float(event["tau_residual"]),
        "coordinate_offset": None if event["time"] is None else float(event["time"]) - float(baseline_time),
        "bracket_width": event["bracket_width"],
        "bisections": event["bisections"],
    }


def _space_error(measured, predicted):
    if measured is None or predicted is None:
        return None
    return float(measured) - float(predicted)


def equal_proper_time_comparison(pair, state, *, duration, step_cap, production=False, fixed_step=None):
    """Baseline, analytic tangent, both centred steps, and the held-out arm.

    Finite arms are measured where their integrated proper time equals the
    baseline proper time. The linear clock correction is reported separately.
    """
    accepted_duration(duration, production=production)
    analytic = predict_occupation_response(
        pair, state, duration=duration, step_cap=step_cap, production=production,
        fixed_step=fixed_step, backend_name=EVOLUTION_BACKEND,
    )
    baseline_march = _march_clock(
        pair, state, duration=duration, step_cap=step_cap, fixed_step=fixed_step,
        production=production, tau_target=None,
    )
    baseline_slice = _coordinate_view(baseline_march)
    tau_target = float(baseline_slice["tau"])
    baseline_event = _equal_tau_event(baseline_march, tau_target)
    zero = discovery.zero_tangent(baseline_march["pair"].grid)
    dots = discovery.nested_readout(
        baseline_march["pair"], baseline_slice["state"], zero, coordinate_duration=0.0,
    )
    amplitudes = tuple(sign * step for step in DERIVATIVE_H for sign in (1.0, -1.0)) + (HELD_OUT_ALPHA,)
    arms = []
    for amplitude in amplitudes:
        prepared = prepare_changed_source(pair, amplitude)
        marched = _march_clock(
            prepared["pair"], prepared["state"], duration=duration, step_cap=step_cap,
            fixed_step=fixed_step, production=production, tau_target=tau_target,
        )
        event = _equal_tau_event(marched, tau_target)
        coordinate = _coordinate_view(marched)
        taylor = _taylor_from_slices(baseline_slice, coordinate, dots)
        arms.append({
            "alpha": float(amplitude),
            "weights": prepared["weights"].copy(),
            "constructor": prepared["constructor"],
            "constructor_grid": prepared["constructor_grid"],
            "method": prepared["method"],
            "converged": True,
            "eigenvalues_admissible": prepared["eigenvalues_admissible"],
            "event": event,
            "coordinate": coordinate,
            "taylor": taylor,
            "history": marched["history"],
            "pair": marched["pair"],
        })

    def find(amplitude):
        matches = [item for item in arms if abs(item["alpha"] - amplitude) <= 1e-15]
        if len(matches) != 1:
            raise ValueError("missing equal-tau arm")
        return matches[0]

    def centred_row(step):
        plus = find(step)
        minus = find(-step)
        if not plus["event"]["reached"] or not minus["event"]["reached"]:
            return {
                "h": float(step),
                "reached": False,
                "measured_primary": None,
                "measured_secondary": None,
                "predicted_primary": analytic["primary_tau"],
                "predicted_secondary": analytic["secondary_tau"],
                "space_error_primary": None,
                "space_error_secondary": None,
                "time_error_plus": _time_error(plus["event"], baseline_slice["coordinate_time"]),
                "time_error_minus": _time_error(minus["event"], baseline_slice["coordinate_time"]),
                "taylor_indicator": {
                    "primary": (plus["taylor"]["primary"] - minus["taylor"]["primary"]) / (2.0 * step),
                    "secondary": (plus["taylor"]["secondary"] - minus["taylor"]["secondary"]) / (2.0 * step),
                    "role": TAYLOR_ROLE,
                    "equal_proper_time": False,
                    "same_formula_as_prediction": True,
                },
            }
        measured_primary = (
            plus["event"]["child_regional_content"] - minus["event"]["child_regional_content"]
        ) / (2.0 * step)
        measured_secondary = (
            plus["event"]["child_proper_mean_r"] - minus["event"]["child_proper_mean_r"]
        ) / (2.0 * step)
        return {
            "h": float(step),
            "reached": True,
            "method": EQUAL_TAU_METHOD,
            "uses_linear_clock_correction": False,
            "measured_primary": measured_primary,
            "measured_secondary": measured_secondary,
            "predicted_primary": analytic["primary_tau"],
            "predicted_secondary": analytic["secondary_tau"],
            "space_error_primary": _space_error(measured_primary, analytic["primary_tau"]),
            "space_error_secondary": _space_error(measured_secondary, analytic["secondary_tau"]),
            "time_error_plus": _time_error(plus["event"], baseline_slice["coordinate_time"]),
            "time_error_minus": _time_error(minus["event"], baseline_slice["coordinate_time"]),
            "dense_interpolation_content_gap_plus": plus["event"]["dense_interpolation_content_gap"],
            "dense_interpolation_content_gap_minus": minus["event"]["dense_interpolation_content_gap"],
            "taylor_indicator": {
                "primary": (plus["taylor"]["primary"] - minus["taylor"]["primary"]) / (2.0 * step),
                "secondary": (plus["taylor"]["secondary"] - minus["taylor"]["secondary"]) / (2.0 * step),
                "role": TAYLOR_ROLE,
                "equal_proper_time": False,
                "same_formula_as_prediction": True,
            },
        }

    held = find(HELD_OUT_ALPHA)
    predicted_change = HELD_OUT_ALPHA * analytic["primary_tau"]
    predicted_secondary = HELD_OUT_ALPHA * analytic["secondary_tau"]
    if held["event"]["reached"]:
        measured_change = held["event"]["child_regional_content"] - baseline_event["child_regional_content"]
        measured_secondary = held["event"]["child_proper_mean_r"] - baseline_event["child_proper_mean_r"]
        held_report = {
            "alpha": HELD_OUT_ALPHA,
            "reached": True,
            "used_as_derivative_h": False,
            "method": EQUAL_TAU_METHOD,
            "uses_linear_clock_correction": False,
            "measured_change": measured_change,
            "measured_secondary_change": measured_secondary,
            "predicted_change": predicted_change,
            "predicted_secondary_change": predicted_secondary,
            "space_error": _space_error(measured_change, predicted_change),
            "space_error_secondary": _space_error(measured_secondary, predicted_secondary),
            "scaled_primary": measured_change / HELD_OUT_ALPHA,
            "scaled_space_error": _space_error(measured_change / HELD_OUT_ALPHA, analytic["primary_tau"]),
            "time_error": _time_error(held["event"], baseline_slice["coordinate_time"]),
            "dense_interpolation_content_gap": held["event"]["dense_interpolation_content_gap"],
            "child_regional_content": held["event"]["child_regional_content"],
            "predicted_content": baseline_event["child_regional_content"] + predicted_change,
        }
    else:
        held_report = {
            "alpha": HELD_OUT_ALPHA,
            "reached": False,
            "used_as_derivative_h": False,
            "method": EQUAL_TAU_METHOD,
            "uses_linear_clock_correction": False,
            "measured_change": None,
            "predicted_change": predicted_change,
            "space_error": None,
            "time_error": _time_error(held["event"], baseline_slice["coordinate_time"]),
            "taylor_indicator": held["taylor"],
        }
    held_report["taylor_indicator"] = held["taylor"]
    return {
        "schema": SCHEMA,
        "hard_acceptance_gate": False,
        "tau_target": tau_target,
        "baseline_coordinate_time": baseline_slice["coordinate_time"],
        "baseline_content": baseline_event["child_regional_content"],
        "baseline_mean_r": baseline_event["child_proper_mean_r"],
        "analytic_primary": analytic["primary_tau"],
        "analytic_secondary": analytic["secondary_tau"],
        "analytic_state": analytic["propagated"]["state"],
        "analytic_tangent": analytic["propagated"]["tangent"],
        "analytic_readout": analytic["readout"],
        "baseline_state": baseline_event["state"],
        "baseline_history": baseline_march["history"],
        "baseline_pair": baseline_march["pair"],
        "centred": [centred_row(step) for step in DERIVATIVE_H],
        "held_out": held_report,
        "arms": arms,
        "effective_stress": None,
        "newton_used_on_analytic_tangent": False,
        "initial_state_called_on_baseline": False,
        "domain": analytic["domain"],
        "eigenvalues_admissible": analytic["eigenvalues_admissible"],
    }


def expected_cost(duration, step_cap):
    """Operation count when the saved step cap binds. Not a wall-clock certificate."""
    duration = float(duration)
    step_cap = float(step_cap)
    if duration <= 0.0 or step_cap <= 0.0 or not np.isfinite(duration) or not np.isfinite(step_cap):
        raise ValueError("cost needs a positive duration and step cap")
    steps = int(np.ceil(duration / step_cap - 1e-12))
    rates_per_step = 4
    return {
        "assumes_step_cap_binds": True,
        "duration": duration,
        "step_cap": step_cap,
        "steps_at_cap": steps,
        "analytic_rate_evaluations": rates_per_step * steps,
        "analytic_jacobian_vector_evaluations": rates_per_step * steps,
        "analytic_clock_samples": rates_per_step * steps,
        "nonlinear_baseline_rate_evaluations": rates_per_step * steps,
        "saved_frame_comparison_reuses_nonlinear_baseline": True,
        "finite_difference_trajectories": 4,
        "finite_difference_rate_evaluations": 4 * rates_per_step * steps,
        "finite_difference_dense_initial_solves": 4,
        "held_out_trajectories": 1,
        "held_out_rate_evaluations": rates_per_step * steps,
        "held_out_dense_initial_solves": 1,
        "effective_stress_evaluations": 0,
        "t03_resolves": 0,
        "equal_tau_root_bisections_per_arm": ROOT_BISECTIONS,
        "equal_tau_root_rate_evaluations": 5 * ROOT_BISECTIONS * 4,
    }


def cpu_forecast(nf, duration, step_cap, *, probe_step_seconds=None, n_cases=1):
    """Single-process CPU estimate. It does not stop a run and it is not a gate."""
    cost = expected_cost(duration, step_cap)
    rate_evaluations = (
        cost["analytic_rate_evaluations"]
        + cost["analytic_jacobian_vector_evaluations"]
        + cost["nonlinear_baseline_rate_evaluations"]
        + cost["finite_difference_rate_evaluations"]
        + cost["held_out_rate_evaluations"]
        + cost["equal_tau_root_rate_evaluations"]
    ) * int(n_cases)
    # Recorded one-thread sample: three nf256 RK4 steps, twelve rate evaluations, 0.037 s.
    seconds_per_rate_nf256 = 0.037 / 12.0
    scale = (float(nf) / 256.0) * (np.log2(max(float(nf), 2.0)) / np.log2(256.0))
    scaled_seconds = seconds_per_rate_nf256 * scale * rate_evaluations
    probe_seconds = None
    if probe_step_seconds is not None:
        probe_seconds = (float(probe_step_seconds) / 4.0) * rate_evaluations
    return {
        "nf": int(nf),
        "n_cases": int(n_cases),
        "duration": float(duration),
        "step_cap": float(step_cap),
        "rate_evaluations": int(rate_evaluations),
        "forecast_cpu_seconds": float(probe_seconds if probe_seconds is not None else scaled_seconds),
        "scaled_sample_cpu_seconds": float(scaled_seconds),
        "probe_cpu_seconds": None if probe_seconds is None else float(probe_seconds),
        "forecast_source": "probe_step" if probe_step_seconds is not None else "backend_nf256_sample_scaled",
        "jacobian_vector_counted_as_one_rate": True,
        "dense_initial_solves_not_in_rate_sample": cost["finite_difference_dense_initial_solves"] + cost["held_out_dense_initial_solves"],
        "single_process": True,
        "process_pool": False,
        "hard_gate": False,
        "stability_certificate": False,
    }


def _probe_step_seconds(pair, state, step_cap):
    resolved, state, info = _resolve(pair, state, EVOLUTION_BACKEND)
    if info.get("backend") != "fft":
        raise ValueError("CPU probe expected the FFT carrier")
    limit, _omega, _quadrature = model.stable_timestep(resolved, state, step_cap)
    step = min(float(step_cap), float(limit))
    started = time.perf_counter()
    nonlinear_rk4_step(resolved, state, step)
    return time.perf_counter() - started, step


def production_cost():
    """Counts for a later T=0.3 run of the four saved baselines. Not executed here."""
    cases = []
    for case_name in episode.BASELINE_CASES:
        cap = 0.001 if case_name.endswith("dt0.001") else 0.0005
        item = expected_cost(HANDOFF_TIME, cap)
        item["case_name"] = case_name
        cases.append(item)
    totals = {
        name: int(sum(item[name] for item in cases))
        for name in (
            "steps_at_cap",
            "analytic_rate_evaluations",
            "analytic_jacobian_vector_evaluations",
            "nonlinear_baseline_rate_evaluations",
            "finite_difference_rate_evaluations",
            "finite_difference_dense_initial_solves",
            "held_out_rate_evaluations",
            "held_out_dense_initial_solves",
        )
    }
    # The two step caps of one nf share the five dense solves.
    totals["finite_difference_dense_initial_solves_if_shared_across_caps"] = 8
    totals["held_out_dense_initial_solves_if_shared_across_caps"] = 2
    return {
        "executed": False,
        "comparison_time": HANDOFF_TIME,
        "cases": cases,
        "totals_if_caps_are_prepared_separately": totals,
        "rate_step_scale": (
            "backend records 0.037 s for three nf256 RK4 steps at one thread; "
            "Jacobian-vector products and dense initial solves are additional"
        ),
    }


def _producer_pins():
    paths = {
        "prediction": _MODULE,
        "response": Path(discovery.__file__).resolve(),
        "backend": Path(backend.__file__).resolve(),
        "nested": Path(model.__file__).resolve(),
        "cli": episode.LAB / "scripts" / "derive_nsc_discovery_prediction.py",
    }
    pins = {name: {"path": str(path), "sha256": episode.file_sha256(path)} for name, path in paths.items()}
    try:
        revision = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=str(episode.REPO), stderr=subprocess.DEVNULL, timeout=5,
        )
        git_head = revision.decode().strip()
    except (OSError, subprocess.SubprocessError):
        git_head = None
    return {"files": pins, "git_head": git_head}


def _arm_token(alpha):
    if abs(alpha - HELD_OUT_ALPHA) <= 1e-15:
        return "alpha_0p01"
    sign = "plus" if alpha > 0.0 else "minus"
    return f"h_{sign}_{abs(alpha):.4f}".replace(".", "p")


def _state_from_prefix(arrays, prefix):
    return model.NestedState(*(arrays[f"{prefix}_{name}"] for name in _STATE_FIELDS))


def _pack_state(arrays, prefix, state):
    for name in _STATE_FIELDS:
        arrays[f"{prefix}_{name}"] = np.array(getattr(state, name), copy=True)


def _pack_tangent(arrays, prefix, tangent):
    for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1", "occupations"):
        arrays[f"{prefix}_{name}"] = np.array(getattr(tangent, name), copy=True)


def _pack_history(arrays, prefix, history):
    for key, values in history.items():
        arrays[f"{prefix}_history_{key}"] = np.asarray(values, dtype=float)


def _json_event(event):
    return {
        key: event[key]
        for key in (
            "reached", "tau", "time", "tau_residual", "bracket_width", "bisections", "method",
            "uses_linear_clock_correction", "child_regional_content", "child_proper_mean_r",
            "tau_dot", "dense_interpolation", "dense_interpolation_content_gap",
        )
        if key in event
    }


def _exclusive_bytes(path, payload):
    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY
    descriptor = os.open(path, flags, 0o644)
    try:
        view = memoryview(payload)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError("exclusive write made no progress")
            view = view[written:]
    finally:
        os.close(descriptor)


def write_run_record(directory, record, arrays):
    """Exclusive JSON and NPZ. Existing successor bytes are left untouched."""
    directory = episode.assert_campaign_output(directory)
    directory.mkdir(parents=True, exist_ok=True)
    json_path = directory / RUN_JSON
    npz_path = directory / RUN_NPZ
    if json_path.exists() or npz_path.exists():
        raise FileExistsError("refusing to overwrite a prediction run record")
    array_hashes = {name: episode.array_sha256(value) for name, value in arrays.items()}
    buffer = io.BytesIO()
    np.savez(buffer, **arrays)
    blob = buffer.getvalue()
    _exclusive_bytes(npz_path, blob)
    payload_hash = episode.file_sha256(npz_path)
    body = dict(record)
    body["array_sha256"] = array_hashes
    body["payload_sha256"] = payload_hash
    body["payload_bytes"] = int(npz_path.stat().st_size)
    _exclusive_bytes(json_path, (json.dumps(episode.jsonable(body), indent=2) + "\n").encode())
    body["json_path"] = str(json_path)
    body["npz_path"] = str(npz_path)
    return body


def run_saved_comparison(nfs, duration, output, *, step_cap=0.0005, production=None):
    """FFT comparison for selected saved resolutions. Writes one JSON and one NPZ."""
    selected = tuple(int(value) for value in nfs)
    if not selected or any(value not in (128, 256) for value in selected):
        raise ValueError("resolutions are selected from nf 128 and nf 256")
    if production is None:
        production = float(duration) > NONPRODUCTION_DURATION + 1e-15
    decision = accepted_duration(duration, production=production)
    if output is not None:
        destination = episode.assert_campaign_output(output)
        if (destination / RUN_JSON).exists() or (destination / RUN_NPZ).exists():
            raise FileExistsError("refusing to overwrite a prediction run record")
    cap = float(step_cap)
    if abs(cap - 0.001) <= 1e-15:
        cap_token = "0.001"
    elif abs(cap - 0.0005) <= 1e-15:
        cap_token = "0.0005"
    else:
        raise ValueError("step cap must be a saved baseline cap, 0.001 or 0.0005")
    cases = []
    arrays = {}
    probe_seconds = None
    for nf in selected:
        case_name = f"nf{nf}_baseline_dt{cap_token}"
        seed = load_saved_baseline(case_name, 0.0)
        if probe_seconds is None:
            probe_seconds, _probe_dt = _probe_step_seconds(seed["pair"], seed["state"], cap)
        comparison = equal_proper_time_comparison(
            seed["pair"], seed["state"], duration=decision["duration"], step_cap=cap, production=decision["production"],
        )
        saved_gap = None
        if abs(decision["duration"] - HANDOFF_TIME) <= 1e-15:
            saved_gap = compare_saved_baseline(comparison["baseline_state"], case_name, HANDOFF_TIME)
        prefix = f"nf{nf}"
        _pack_state(arrays, f"{prefix}_baseline_event", comparison["baseline_state"])
        _pack_state(arrays, f"{prefix}_analytic_state", comparison["analytic_state"])
        _pack_tangent(arrays, f"{prefix}_analytic_tangent", comparison["analytic_tangent"])
        _pack_history(arrays, f"{prefix}_baseline", comparison["baseline_history"])
        arrays[f"{prefix}_W"] = np.array(seed["pair"].geometry_map, copy=True)
        arm_summaries = []
        for item in comparison["arms"]:
            token = _arm_token(item["alpha"])
            arm_prefix = f"{prefix}_{token}"
            if item["event"]["reached"]:
                _pack_state(arrays, f"{arm_prefix}_event", item["event"]["state"])
            _pack_state(arrays, f"{arm_prefix}_coordinate", item["coordinate"]["state"])
            _pack_history(arrays, arm_prefix, item["history"])
            arm_summaries.append({
                "alpha": item["alpha"],
                "token": token,
                "weights": [float(value) for value in item["weights"]],
                "constructor": item["constructor"],
                "constructor_grid": item["constructor_grid"],
                "method": item["method"],
                "converged": item["converged"],
                "eigenvalues_admissible": item["eigenvalues_admissible"],
                "event_prefix": f"{arm_prefix}_event" if item["event"]["reached"] else None,
                "coordinate_prefix": f"{arm_prefix}_coordinate",
                "equal_tau": _json_event(item["event"]),
                "taylor_indicator": item["taylor"],
                "coordinate_time": item["coordinate"]["coordinate_time"],
                "coordinate_tau": item["coordinate"]["tau"],
                "coordinate_content": item["coordinate"]["readout"]["child_regional_content"],
            })
        cases.append({
            "case_name": case_name,
            "nf": nf,
            "step_cap": cap,
            "duration": decision["duration"],
            "target": decision["target"],
            "tau_target": comparison["tau_target"],
            "baseline_coordinate_time": comparison["baseline_coordinate_time"],
            "baseline_content": comparison["baseline_content"],
            "baseline_mean_r": comparison["baseline_mean_r"],
            "baseline_prefix": f"{prefix}_baseline_event",
            "analytic_primary": comparison["analytic_primary"],
            "analytic_secondary": comparison["analytic_secondary"],
            "analytic_state_prefix": f"{prefix}_analytic_state",
            "analytic_tangent_prefix": f"{prefix}_analytic_tangent",
            "W_key": f"{prefix}_W",
            "W_sha256": seed["pins"]["W_sha256"],
            "seed_state_sha256": seed["state_sha256"],
            "weights_sha256": seed["pins"]["weights_sha256"],
            "source_pins": seed["pins"],
            "centred_equal_tau": comparison["centred"],
            "held_out_equal_tau": comparison["held_out"],
            "arms": arm_summaries,
            "saved_T03_gap": None if saved_gap is None else {
                "max_gap": saved_gap["max_gap"],
                "gaps": saved_gap["gaps"],
                "source_changed": False,
                "constraint_resolved": False,
                "hard_gate": False,
            },
            "t03_source_changed": False,
            "t03_constraint_resolved": False,
            "hard_acceptance_gate": False,
            "domain": comparison["domain"],
        })
    forecast = cpu_forecast(
        selected[-1], decision["duration"], cap, probe_step_seconds=probe_seconds, n_cases=len(selected),
    )
    record = {
        "schema": SCHEMA,
        "status": "RUN",
        "creation_only": True,
        "evolved": True,
        "production_trajectory": bool(decision["production"]),
        "production_launched_by": "campaign_executor" if decision["production"] else "short_test",
        "replay_recomputes_observables_from_saved_states": True,
        "evolution_on_replay": False,
        "hard_acceptance_gate": False,
        "equal_proper_time_method": EQUAL_TAU_METHOD,
        "taylor_indicator_role": TAYLOR_ROLE,
        "duration": decision["duration"],
        "target": decision["target"],
        "nfs": list(selected),
        "step_cap": cap,
        "backend": EVOLUTION_BACKEND,
        "preparation_backend": PREPARATION_BACKEND,
        "first_prediction_target": FIRST_PREDICTION_TARGET,
        "optional_prediction_target": OPTIONAL_PREDICTION_TARGET,
        "late_T3_curvature_unresolved": True,
        "single_process": True,
        "process_pool": False,
        "cpu_forecast": forecast,
        "inputs": saved_input_hashes(),
        "producers": _producer_pins(),
        "cases": cases,
        "effective_stress": None,
        "settings": specification(),
    }
    if output is None:
        record["arrays"] = arrays
        return record
    return write_run_record(output, record, arrays)


def replay_run_record(directory):
    """Recompute saved-state observables. Do not step and do not solve."""
    directory = Path(directory).expanduser().resolve()
    json_path = directory / RUN_JSON
    npz_path = directory / RUN_NPZ
    if not json_path.is_file() or not npz_path.is_file():
        raise FileNotFoundError("prediction run record is incomplete")
    record = json.loads(json_path.read_text())
    if record.get("schema") != SCHEMA or record.get("status") != "RUN":
        raise ValueError("prediction run record schema is wrong")
    if episode.file_sha256(npz_path) != record["payload_sha256"]:
        raise ValueError("prediction payload does not match its record")
    with np.load(npz_path, allow_pickle=False) as stored:
        arrays = {name: np.array(stored[name], copy=True) for name in stored.files}
    for name, digest in record["array_sha256"].items():
        if episode.array_sha256(arrays[name]) != digest:
            raise ValueError("prediction array changed: " + name)
    current_inputs = saved_input_hashes()
    for key in ("v1_json_sha256", "v1_npz_sha256", "basis_json_sha256", "basis_npz_sha256"):
        if record["inputs"][key] != current_inputs[key]:
            raise ValueError("saved input changed after the run: " + key)
    producers = _producer_pins()
    producer_match = all(
        producers["files"][name]["sha256"] == record["producers"]["files"][name]["sha256"]
        for name in record["producers"]["files"]
    )
    if not producer_match:
        raise ValueError("producer bytes changed after the run record was written")
    gaps = []
    residuals = []
    for case in record["cases"]:
        seed = load_saved_baseline(case["case_name"], 0.0)
        if episode.array_sha256(arrays[case["W_key"]]) != seed["pins"]["W_sha256"]:
            raise ValueError("frozen W in the run record does not match the saved pin")
        if case["t03_constraint_resolved"] or case["t03_source_changed"]:
            raise ValueError("run record claims a T=0.3 resolve or a source edit")
        baseline = _state_from_prefix(arrays, case["baseline_prefix"])
        base_pair, _baseline_state, _info = episode.resolve_pair(
            pair_with_occupations(seed["pair"], seed["weights"]), baseline, backend="fft",
        )
        base_observed = _observables(base_pair, baseline)
        gaps.append(abs(base_observed["child_regional_content"] - case["baseline_content"]))
        tangent = arrays[f"{case['analytic_tangent_prefix']}_occupations"]
        if not np.array_equal(np.asarray(tangent, dtype=float), OCCUPATION_CONTRAST):
            raise ValueError("analytic clock tangent lost the occupation contrast")
        for arm in case["arms"]:
            if not arm["equal_tau"]["reached"]:
                continue
            state = _state_from_prefix(arrays, arm["event_prefix"])
            pair, _arm_state, _info = episode.resolve_pair(
                pair_with_occupations(seed["pair"], arm["weights"]), state, backend="fft",
            )
            observed = _observables(pair, state)
            gaps.append(abs(observed["child_regional_content"] - arm["equal_tau"]["child_regional_content"]))
            gaps.append(abs(observed["child_proper_mean_r"] - arm["equal_tau"]["child_proper_mean_r"]))
            residuals.append(abs(float(arm["equal_tau"]["tau_residual"])))
            if arm["equal_tau"]["uses_linear_clock_correction"]:
                raise ValueError("equal-tau event was marked as the linear clock correction")
        for row in case["centred_equal_tau"]:
            if row["reached"] and row["uses_linear_clock_correction"]:
                raise ValueError("centred equal-tau row uses the linear clock correction")
            if row["taylor_indicator"]["equal_proper_time"]:
                raise ValueError("Taylor indicator was stored as the equal-tau measurement")
    return {
        "schema": SCHEMA,
        "ok": True,
        "evolved": True,
        "replayed": True,
        "evolution_called": False,
        "readonly": True,
        "cases": len(record["cases"]),
        "max_observable_gap": float(max(gaps) if gaps else 0.0),
        "max_tau_residual": float(max(residuals) if residuals else 0.0),
        "producer_match": True,
        "hard_acceptance_gate": False,
        "cpu_forecast": record["cpu_forecast"],
    }


def saved_input_hashes():
    v1, basis = _saved_documents()
    return {
        "v1_json": str(episode.V1_JSON),
        "v1_json_sha256": episode.file_sha256(episode.V1_JSON),
        "v1_npz": str(episode.V1_NPZ),
        "v1_npz_sha256": v1["payload_sha256"],
        "basis_json": str(episode.BASIS_JSON),
        "basis_json_sha256": episode.file_sha256(episode.BASIS_JSON),
        "basis_npz": str(episode.BASIS_NPZ),
        "basis_npz_sha256": basis["payload_sha256"],
        "prediction_py_sha256": episode.file_sha256(_MODULE),
        "response_py_sha256": episode.file_sha256(Path(discovery.__file__).resolve()),
    }


def create_successor(output):
    """Write the successor declaration. Do not evolve and do not emit trajectories."""
    directory = episode.assert_campaign_output(output)
    manifest_path = directory / MANIFEST_NAME
    if manifest_path.exists() or any((directory / name).exists() for name in LATER_OUTPUTS):
        raise FileExistsError("refusing to overwrite prediction successor files")
    frames = []
    for case_name in episode.BASELINE_CASES:
        seed = load_saved_baseline(case_name, 0.0)
        final = load_saved_baseline(case_name, HANDOFF_TIME)
        if seed["pins"]["weights_sha256"] != final["pins"]["weights_sha256"]:
            raise ValueError("saved T=0.3 occupations differ from the T=0 source")
        if seed["initial_state_called"] or final["constraint_resolved"] or final["source_changed"]:
            raise ValueError("loading a saved frame must not solve or change the source")
        frames.append({
            "case_name": case_name,
            "nf": seed["nf"],
            "step_cap": seed["step_cap"],
            "seed_time": seed["time"],
            "comparison_time": final["time"],
            "initial_state_sha256": seed["initial_state_sha256"],
            "final_state_sha256": seed["final_state_sha256"],
            "seed_state_sha256": seed["state_sha256"],
            "comparison_state_sha256": final["state_sha256"],
            "W_sha256": seed["pins"]["W_sha256"],
            "weights_sha256": seed["pins"]["weights_sha256"],
            "source_columns_sha256": seed["pins"]["source_columns_sha256"],
            "clock_locations": seed["clock_locations"],
            "domain": seed["domain"],
            "t03_source_changed": False,
            "t03_constraint_resolved": False,
        })
    directory.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema": SCHEMA,
        "status": "CREATED",
        "creation_only": True,
        "evolved": False,
        "production_trajectory": False,
        "production_launched": False,
        "campaign_executor": "sole campaign executor after freeze",
        "settings": specification(),
        "inputs": saved_input_hashes(),
        "frames": frames,
        "later_outputs": list(LATER_OUTPUTS),
        "later_outputs_written": False,
        "expected_cost": production_cost(),
        "effective_stress": None,
    }
    manifest_path.write_text(json.dumps(episode.jsonable(manifest), indent=2) + "\n")
    return manifest


def check_successor(output):
    """Reload a creation manifest, or replay a run record without evolving."""
    directory = Path(output).expanduser().resolve()
    if (directory / RUN_JSON).is_file():
        return replay_run_record(directory)
    manifest_path = directory / MANIFEST_NAME
    if not manifest_path.is_file():
        raise FileNotFoundError("prediction manifest is missing")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema") != SCHEMA or manifest.get("creation_only") is not True:
        raise ValueError("prediction manifest schema or creation flag is wrong")
    if manifest.get("evolved") is not False or manifest.get("production_trajectory") is not False:
        raise ValueError("creation manifest records an evolution")
    present = [name for name in LATER_OUTPUTS if (directory / name).exists()]
    if present:
        raise ValueError("later successor outputs were written: " + ", ".join(present))
    current = saved_input_hashes()
    recorded = manifest["inputs"]
    for key in ("v1_json_sha256", "v1_npz_sha256", "basis_json_sha256", "basis_npz_sha256"):
        if recorded[key] != current[key]:
            raise ValueError("saved input changed after the manifest: " + key)
    settings = manifest["settings"]
    if settings["derivative_h"] != list(DERIVATIVE_H) or settings["held_out_alpha"] != HELD_OUT_ALPHA:
        raise ValueError("manifest derivative steps or held-out alpha changed")
    if HELD_OUT_ALPHA in settings["derivative_h"]:
        raise ValueError("held-out alpha was stored as a derivative step")
    if settings["effective_stress_evaluated"] is not False or settings["t03_constraint_resolved"] is not False:
        raise ValueError("manifest records a stress evaluation or a T=0.3 resolve")
    return {
        "schema": SCHEMA,
        "ok": True,
        "evolved": False,
        "frames": len(manifest["frames"]),
        "later_outputs_written": False,
    }
