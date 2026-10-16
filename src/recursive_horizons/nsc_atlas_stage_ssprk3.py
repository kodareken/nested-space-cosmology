"""Explicit callback SSPRK(3,3) for the existing coupled departure state.

Only time scheduling is owned here. A supplied stage owns all physical
reconstruction and rates, including final diagnostics. Existing state,
configuration, midpoint and combination helpers retain their conventions.
No operator is installed through a module-global replacement.
"""
import copy
import math
from numbers import Real
import time

import numpy as np

from .nsc_atlas_coupled_evolution import _finite, _combine, CoupledEvolutionError
from .nsc_atlas_compatible_evolution import _state
from .nsc_atlas_background_ssprk3 import _times, _configuration, TIME_SCHEME
from .nsc_atlas_characteristic_initialization import Budget


def _shared_budget(budget):
    """Validate the existing clock/count without replacing either one."""
    if (not isinstance(budget, Budget) or type(budget.maximum_evaluations) is not int
            or budget.maximum_evaluations <= 0 or type(budget.evaluations) is not int
            or budget.evaluations < 0):
        raise ValueError("one finite shared caller Budget is required")
    for name in ("started", "cpu_seconds"):
        value = getattr(budget, name)
        if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
                or not math.isfinite(value) or value < 0
                or (name == "cpu_seconds" and value == 0)):
            raise ValueError("shared budget clock and positive CPU limit must be finite")
    if budget.started > time.process_time():
        raise ValueError("shared budget cannot start in future process time")
    budget.check()


def advance_stage_ssprk3(initial, x_left, x_mid, x_right, *, stage, eta, y_faces,
                         source, gauge, alpha, tube, budget, profile,
                         maximum_courant=.4, exterior_quadrature_order=2,
                         cap_quadrature_order=4):
    """Three callback stages at exact supplied left, right, midpoint keys.

    The callback has the existing background-stage keyword contract. Bind
    additional operator choices with an explicit closure or functools.partial.
    All attempted callback calls, including a failed call, are counted. Input
    lookup failure is distinguished from operator failure. No endpoint stage,
    retry, field reset or source replacement is hidden in this function.
    """
    phase, calls, max_courant = "validation", 0, 0.
    state, attempt = copy.deepcopy(initial), copy.deepcopy(initial)
    stages, stage_states, times = [], [], None
    try:
        _shared_budget(budget)
        if not callable(stage):
            raise ValueError("stage must be callable")
        left, mid, right = _times(x_left, x_mid, x_right)
        times = [left, right, mid]
        eta, y_faces, alpha, maximum_courant = _configuration(
            eta, y_faces, source, gauge, alpha, budget, maximum_courant,
            exterior_quadrature_order, cap_quadrature_order)
        state = _state(initial, len(eta)-1, len(y_faces)-1)
        attempt = state
        h = right-left

        def evaluate(value, i):
            nonlocal phase, calls, max_courant, attempt
            attempt = value
            phase = f"stage_{i}_inputs"
            src, ref = source(times[i]), gauge(times[i])
            phase = f"stage_{i}_operator"
            calls += 1
            row = stage(value, x=times[i], eta=eta, y_faces=y_faces,
                source=src, gauge=ref, alpha=alpha, tube=tube, budget=budget, profile=profile,
                exterior_quadrature_order=exterior_quadrature_order,
                cap_quadrature_order=cap_quadrature_order)
            stages.append(row)
            stage_states.append(copy.deepcopy(value))
            courant_rate = float(_finite(row["courant_rate"], "courant_rate", ()))
            if courant_rate < 0:
                raise ValueError("negative characteristic Courant rate")
            max_courant = max(max_courant, h*courant_rate)
            if max_courant > maximum_courant:
                raise ValueError("declared characteristic Courant limit exceeded")
            budget.check()
            return row

        first = evaluate(state, 0)
        phase = "first predictor"
        predictor = _state(_combine(state, first["rates"], 1., h), len(eta)-1, len(y_faces)-1)
        second = evaluate(predictor, 1)
        phase = "midpoint convex combination"
        midpoint = _state(_combine(state, _combine(predictor, second["rates"], 1., h), .75, .25),
                          len(eta)-1, len(y_faces)-1)
        third = evaluate(midpoint, 2)
        phase = "accepted convex combination"
        accepted = _state(_combine(state, _combine(midpoint, third["rates"], 1., h), 1/3, 2/3),
                          len(eta)-1, len(y_faces)-1)
        budget.check()
    except (ValueError, RuntimeError, FloatingPointError, OverflowError, KeyError, TypeError) as error:
        raise CoupledEvolutionError(dict(phase=phase, stage_calls=calls, completed_stages=len(stages),
            stage_times=times, charged_evaluations=getattr(budget, "evaluations", None),
            max_courant=max_courant, reason=str(error), accepted_state=copy.deepcopy(state),
            attempt_state=copy.deepcopy(attempt))) from error
    return dict(accepted_state=accepted, stage_states=stage_states, stages=stages,
                stage_times=times, stage_calls=calls, max_courant=max_courant,
                time_scheme=TIME_SCHEME, retries=0)


def evolve_stage_ssprk3(initial, x_nodes, *, stage, mid_nodes, eta, y_faces,
                        source, gauge, alpha, tube, budget, profile, save_every=1,
                        maximum_courant=.4, exterior_quadrature_order=2,
                        cap_quadrature_order=4):
    """Uniform explicit-callback evolution, with 3*N+1 successful calls.

    The final callback is the same one used at intermediate stages. Saved
    accepted states and physical cap totals preserve the existing schema.
    Courant enforcement does not assert stability for arbitrary callbacks.
    """
    phase, step, calls, max_courant = "validation", 0, 0, 0.
    state, attempt = copy.deepcopy(initial), copy.deepcopy(initial)
    snapshots, local_failure = [], None
    try:
        _shared_budget(budget)
        if not callable(stage):
            raise ValueError("stage must be callable")
        eta, y_faces, alpha, maximum_courant = _configuration(
            eta, y_faces, source, gauge, alpha, budget, maximum_courant,
            exterior_quadrature_order, cap_quadrature_order)
        x, mids = _finite(x_nodes, "x_nodes"), _finite(mid_nodes, "mid_nodes")
        if (x.ndim != 1 or not 2 <= len(x) <= 2049 or mids.shape != (len(x)-1,)
                or np.any(np.diff(x) <= 0)
                or not np.allclose(np.diff(x), np.diff(x)[0], rtol=1e-10, atol=1e-14)):
            raise ValueError("require 1..2048 uniform positive steps and one midpoint per step")
        if type(save_every) is not int or save_every < 1 or (len(x)-1) % save_every:
            raise ValueError("save_every must divide the number of steps")
        for i in range(len(mids)):
            _times(x[i], mids[i], x[i+1])
        state = _state(initial, len(eta)-1, len(y_faces)-1)

        def save(at, value, result):
            snapshots.append(dict(x=float(at), state=copy.deepcopy(value),
                physical=dict(increment=copy.deepcopy(result["physical_increment"]),
                              wall=copy.deepcopy(result["physical_wall"]),
                              cap=copy.deepcopy(result["physical_cap"]))))

        for step in range(len(x)-1):
            phase, attempt = "step", state
            try:
                result = advance_stage_ssprk3(state, x[step], mids[step], x[step+1],
                    stage=stage, eta=eta, y_faces=y_faces, source=source, gauge=gauge, alpha=alpha,
                    tube=tube, budget=budget, profile=profile, maximum_courant=maximum_courant,
                    exterior_quadrature_order=exterior_quadrature_order,
                    cap_quadrature_order=cap_quadrature_order)
            except CoupledEvolutionError as error:
                local_failure = error.diagnostic
                calls += local_failure["stage_calls"]
                max_courant = max(max_courant, local_failure["max_courant"])
                attempt = local_failure["attempt_state"]
                raise
            calls += result["stage_calls"]
            max_courant = max(max_courant, result["max_courant"])
            if step % save_every == 0:
                save(x[step], state, result["stages"][0])
            state = result["accepted_state"]
        phase, attempt = "final_inputs", state
        src, ref = source(float(x[-1])), gauge(float(x[-1]))
        phase = "final_operator"
        calls += 1
        final = stage(state, x=float(x[-1]), eta=eta, y_faces=y_faces,
            source=src, gauge=ref, alpha=alpha, tube=tube, budget=budget, profile=profile,
            exterior_quadrature_order=exterior_quadrature_order,
            cap_quadrature_order=cap_quadrature_order)
        # A final callback still has to return a finite nonnegative speed,
        # although it is not followed by another step.
        if float(_finite(final["courant_rate"], "courant_rate", ())) < 0:
            raise ValueError("negative characteristic Courant rate")
        save(x[-1], state, final)
        budget.check()
    except (ValueError, RuntimeError, FloatingPointError, OverflowError, KeyError, TypeError) as error:
        raise CoupledEvolutionError(dict(phase=phase, step=step, stage_calls=calls,
            charged_evaluations=getattr(budget, "evaluations", None), max_courant=max_courant,
            local_failure=local_failure, reason=str(error), accepted_state=copy.deepcopy(state),
            attempt_state=copy.deepcopy(attempt))) from error
    return dict(snapshots=snapshots, final_stage=final, steps=len(x)-1, stages=calls,
                max_courant=max_courant, time_scheme=TIME_SCHEME, retries=0)
