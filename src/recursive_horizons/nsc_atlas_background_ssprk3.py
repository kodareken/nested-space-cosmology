"""SSPRK(3,3) on the existing coupled characteristic-background equations.

The Shu--Osher coefficients are standard (Gottlieb, Shu and Tadmor, 2001,
equation 4.2). This changes only the time integrator. Raw departures retain
their existing meaning, and each stage restores its source and reference
geometry at the actual supplied time. No new PDE stability or convergence
theorem is asserted; knots and reconstruction switches need measured controls.
"""
import copy
import math
import time

import numpy as np

from .nsc_atlas_characteristic_initialization import Budget
from .nsc_atlas_coupled_evolution import _finite, _grid, _combine, CoupledEvolutionError
from .nsc_atlas_compatible_evolution import _state
from .nsc_atlas_background_evolution import background_stage


TIME_SCHEME = "SSPRK(3,3); nonautonomous stage times 0,1,1/2"


def _times(left, mid, right):
    left, mid, right = (float(_finite(v, name, ())) for v, name in
                        ((left, "x_left"), (mid, "x_mid"), (right, "x_right")))
    h = right-left
    if not math.isfinite(h) or not left < mid < right:
        raise ValueError("require a finite positive interval and interior midpoint")
    expected = left+h/2
    tolerance = 4*max(math.ulp(v) for v in (left, mid, right, expected))
    if abs(mid-expected) > tolerance:
        raise ValueError("supplied midpoint is not the arithmetic midpoint within four ulps")
    # Preserve the caller's floating-point cache key, including an admitted
    # last-bit difference from the recomputed arithmetic midpoint.
    return left, mid, right


def _configuration(eta, y_faces, source, gauge, alpha, budget, maximum_courant,
                   exterior_quadrature_order, cap_quadrature_order):
    if (not isinstance(budget, Budget) or type(budget.maximum_evaluations) is not int
            or budget.maximum_evaluations <= 0 or type(budget.evaluations) is not int
            or budget.evaluations < 0 or isinstance(budget.cpu_seconds, bool)
            or not math.isfinite(budget.cpu_seconds) or budget.cpu_seconds <= 0
            or not math.isfinite(budget.started) or budget.started > time.process_time()):
        raise ValueError("one finite shared caller Budget is required")
    budget.check()
    eta, y_faces = _grid(eta, "eta", 2, 512), _grid(y_faces, "y_faces", 3, 128)
    alpha = float(_finite(alpha, "alpha", ()))
    maximum_courant = float(_finite(maximum_courant, "maximum_courant", ()))
    if alpha < 0 or not 0 < maximum_courant <= .5:
        raise ValueError("require nonnegative alpha and 0<maximum_courant<=.5")
    if not callable(source) or not callable(gauge):
        raise ValueError("source and gauge must be exact-coordinate callables")
    if (type(exterior_quadrature_order) is not int or exterior_quadrature_order not in (2, 4)
            or type(cap_quadrature_order) is not int or cap_quadrature_order not in (4, 8, 16)):
        raise ValueError("unsupported frozen-operator quadrature order")
    return eta, y_faces, alpha, maximum_courant


def advance_background_ssprk3(initial, x_left, x_mid, x_right, *, eta, y_faces,
                              source, gauge, alpha, tube, budget, profile,
                              maximum_courant=.4, exterior_quadrature_order=2,
                              cap_quadrature_order=4):
    """One three-stage step using supplied exact left/mid/right cache keys.

    Returned ``stage_states`` and ``stages`` are ordered left, right, midpoint.
    ``accepted_state`` is the raw departure state; no endpoint field is reset
    or independently prescribed. There is no extra accepted-endpoint stage.
    Failures expose actual attempted ``stage_calls`` (including zero during
    validation) in CoupledEvolutionError.diagnostic, with no retry.
    """
    phase, calls, max_courant = "validation", 0, 0.
    state, attempt = copy.deepcopy(initial), copy.deepcopy(initial)
    stages, stage_states, times = [], [], None
    try:
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
            row = background_stage(value, x=times[i], eta=eta, y_faces=y_faces,
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


def evolve_background_ssprk3(initial, x_nodes, *, mid_nodes, eta, y_faces, source,
                             gauge, alpha, tube, budget, profile, save_every=1,
                             maximum_courant=.4, exterior_quadrature_order=2,
                             cap_quadrature_order=4):
    """Uniform SSPRK3 evolution with exact supplied midpoint keys and no retries.

    Snapshots retain the old background solver's structure. Saved starts use
    the already computed left stage; a single final stage supplies the final
    diagnostics. Thus N successful steps use exactly 3*N+1 stage calls.
    """
    phase, step, calls, max_courant = "validation", 0, 0, 0.
    state, attempt = copy.deepcopy(initial), copy.deepcopy(initial)
    snapshots = []
    local_failure = None
    try:
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
            phase = "step"
            attempt = state
            try:
                result = advance_background_ssprk3(state, x[step], mids[step], x[step+1],
                    eta=eta, y_faces=y_faces, source=source, gauge=gauge, alpha=alpha,
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
        phase = "final_inputs"
        attempt = state
        src, ref = source(float(x[-1])), gauge(float(x[-1]))
        phase = "final_operator"
        calls += 1
        final = background_stage(state, x=float(x[-1]), eta=eta, y_faces=y_faces,
            source=src, gauge=ref, alpha=alpha, tube=tube, budget=budget, profile=profile,
            exterior_quadrature_order=exterior_quadrature_order,
            cap_quadrature_order=cap_quadrature_order)
        save(x[-1], state, final)
        budget.check()
    except (ValueError, RuntimeError, FloatingPointError, OverflowError, KeyError, TypeError) as error:
        raise CoupledEvolutionError(dict(phase=phase, step=step, stage_calls=calls,
            charged_evaluations=getattr(budget, "evaluations", None), max_courant=max_courant,
            local_failure=local_failure, reason=str(error), accepted_state=copy.deepcopy(state),
            attempt_state=copy.deepcopy(attempt))) from error
    return dict(snapshots=snapshots, final_stage=final, steps=len(x)-1, stages=calls,
                max_courant=max_courant, time_scheme=TIME_SCHEME, retries=0)
