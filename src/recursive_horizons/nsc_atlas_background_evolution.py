"""Coupled exterior/wall evolution with characteristic cap departures.

Only the cap representation changes from compatible_evolution. Exterior
vacuum departures, first-order transverse jets, source data, wall geometry
and clocks use the same fixed-step equations. The cap stores departures
in Ip/Iq about one CharacteristicProfile; the total wall primitive remains
its single canonical scalar value. This adapter supplies no physical source
preparation, no profile fitting and no independent wall scalar trajectory.
"""
import copy
import math

import numpy as np

from .nsc_atlas_characteristic_initialization import Budget
from .nsc_atlas_coupled_evolution import _finite, _grid, _combine, CoupledEvolutionError
from .nsc_atlas_compatible_evolution import FIELDS, WALL_FIELDS, _state, quiet_departures
from .nsc_atlas_nonlinear_boundary import hybrid_wall_boundary
from .nsc_atlas_vacuum_departure import vacuum_slice, departure_slice
from .nsc_atlas_background_cap import background_cap_metric, background_cap_rhs


def background_stage(state, *, x, eta, y_faces, source, gauge, alpha, tube,
                     budget, profile, exterior_quadrature_order=2, cap_quadrature_order=4):
    """One coupled stage, retaining independent boundary compatibility defects.

    Only W_ref(x) is a prescribed wall coordinate. Physical radius, lapse and
    both clocks are reconstructed from their evolved departures. The original
    source values and first jets supply the source endpoint without replacement.
    Cap state Ip/Iq are departures; physical_cap contains the actual integrals.
    """
    if not isinstance(budget, Budget):
        raise ValueError("a shared caller Budget is required")
    budget.check()
    eta = _grid(eta, "eta", 2, 512)
    y_faces = _grid(y_faces, "y_faces", 3, 128)
    state = _state(state, len(eta)-1, len(y_faces)-1)
    x, alpha = float(_finite(x, "x", ())), float(_finite(alpha, "alpha", ()))
    if alpha < 0:
        raise ValueError("alpha must be nonnegative")
    W, Wx, boost, boost_x = (float(_finite(gauge[k], "gauge."+k, ()))
                            for k in ("W", "W_x", "eta_plus", "eta_plus_x"))
    if W <= 0 or Wx <= 0 or not math.isclose(boost, -.5*math.log(Wx), rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("require positive W,W_x and their compatible boost")
    reference_wall = {key: float(_finite(gauge[key], "gauge."+key, ()))
                      for key in WALL_FIELDS}
    reference_rates = {key: float(_finite(gauge[name], "gauge."+name, ()))
                       for key, name in (("r", "r_x_total"), ("L", "L_x_total"),
                                         ("T", "T_x"), ("s", "s_x"))}
    wall = {key: reference_wall[key]+float(state["wall"][key]) for key in WALL_FIELDS}
    metric = background_cap_metric(y_faces, wall["r"], state["cap"],
        T=wall["T"], profile=profile, alpha=alpha,
        quadrature_order=cap_quadrature_order, budget=budget)
    wall["c"] = float(metric["wall"]["c"])
    vacuum = vacuum_slice(eta*W, x, budget=budget)
    src = {key: float(_finite(source["log_C" if key == "L" else key], "source."+key, ()))
           for key in FIELDS}
    base = {key: vacuum[key] if key != "c" else np.zeros_like(eta) for key in FIELDS}
    offsets = {key: src[key]-base[key][0] for key in FIELDS}
    departures = {key: np.r_[0., state["departure"][key],
                             wall[key]-base[key][-1]-offsets[key]] for key in FIELDS}
    X = {key: base[key]+offsets[key]+departures[key] for key in FIELDS}
    for key in FIELDS:
        X[key][0], X[key][-1] = src[key], wall[key]
    # Differentiate only departures; negative grid speed enters from wall.
    gradients = {key: np.r_[np.diff(departures[key]), departures[key][-1]-departures[key][-2]]
                       /np.diff(eta)[0] for key in FIELDS}
    delta_jets = {"r_Ubar": gradients["r"]/W,
                  "log_C_Ubar": gradients["L"]/W, "c_Ubar": gradients["c"]/W}
    ext = departure_slice(eta*W, X, delta_jets, source, x=x, alpha=alpha,
                          tube=tube, budget=budget, quadrature_order=exterior_quadrature_order)
    outgoing = {"r_V": float(ext["r_x"][-1]), "c_V": float(ext["c_x"][-1]),
                "log_C_V": float(ext["log_C_x"][-1])}
    common = {"r": wall["r"], "c": wall["c"], "eta_plus": boost}
    cap_wall = {key: metric["wall"][key] for key in ("A", "delta", "w_minus")}
    C = math.exp(wall["L"])
    first = hybrid_wall_boundary(common, outgoing, cap_wall,
                                 {"C_wall": C, "log_C_s": 0.}, alpha=alpha)
    xs = float(first["Vdot_plus"])
    Lx = 2*boost_x-2*float(first["k_plus"])/xs+2*outgoing["log_C_V"]
    hybrid = hybrid_wall_boundary(common, outgoing, cap_wall,
                                  {"C_wall": C, "log_C_s": Lx*xs}, alpha=alpha)
    wall_rates = {"r": float(hybrid["r_s"])/xs, "L": Lx,
                  "T": float(hybrid["T_s"])/xs, "s": 1/xs,
                  "c": float(hybrid["c_s"])/xs}
    cap = background_cap_rhs(metric, T_x=wall_rates["T"], radius_x=wall_rates["r"],
                           incoming_w_plus=float(hybrid["w_plus"]), budget=budget)
    # The total wall primitive includes the moving analytic background.
    # Sum(departure Iq rates) alone omits its physical wall contribution.
    chi_x = float(cap["wall_chi_x"])
    f = float(metric["wall"]["f"])
    f_x = (chi_x-f*wall_rates["r"])/wall["r"]
    u = wall["c"]/wall["r"]
    from_flux_cx = u*wall_rates["r"]+wall["r"]*math.sqrt(1+alpha*u*u/6)*f_x
    delta_q = {"r": ext["delta_r_x"], "L": ext["delta_log_C_x"], "c": ext["delta_c_x"]}
    source_delta_q = {"r": float(source["r_x"])-1., "L": float(source["log_C_x"]),
                      "c": float(source["c_x"])}
    rates = {"departure": {key: (delta_q[key]-source_delta_q[key]
                                   +eta*(Wx/W)*gradients[key])[1:-1] for key in FIELDS},
             "wall": {key: wall_rates[key]-reference_rates[key] for key in WALL_FIELDS},
             "cap": cap["rates"]}
    jets = {key: ext[key] for key in ("r_Ubar", "log_C_Ubar", "c_Ubar")}
    q = {"r": ext["r_x"], "L": ext["log_C_x"], "c": ext["c_x"]}
    p = {"r": jets["r_Ubar"], "L": jets["log_C_Ubar"], "c": jets["c_Ubar"]}
    residual = {key: wall_rates[key]-q[key][-1]-Wx*p[key][-1] for key in FIELDS}
    if not all(np.all(np.isfinite(v)) for group in rates.values() for v in group.values()):
        raise ValueError("nonfinite background coupled rate")
    budget.check(2)
    return {"rates": rates, "physical_wall": wall, "wall_rates": wall_rates,
            "physical_increment": {key: X[key][1:-1]-src[key] for key in FIELDS},
            "physical_cap": {key: value.copy() for key, value in metric["physical_integrals"].items()},
            "exterior": X, "exterior_first_jets": jets, "exterior_longitudinal": q,
            "exterior_prefix": ext, "cap_metric": metric, "cap_rhs": cap, "hybrid": hybrid,
            "wall_kinematic_defect": residual,
            "canonical_wall_flux_defect": chi_x-wall_rates["r"]*f-wall["r"]*float(hybrid["f_s"])/xs,
            "conformal_wall_flux_defect": from_flux_cx-wall_rates["c"],
            "vacuum_wall_anchor_defect": {key: base[key][-1]-reference_wall[key] for key in ("r", "L")},
            "courant_rate": max(float(Wx/W/np.diff(eta)[0]), cap["characteristic_courant_rate"])}


def evolve_background(initial, x_nodes, *, eta, y_faces, source, gauge, alpha,
                       tube, budget, profile, save_every=1, maximum_courant=.4,
                       exterior_quadrature_order=2, cap_quadrature_order=4):
    """Fixed SSPRK2 with one caller-supplied profile and no retries.

    Raw cap state stores departures in physical radial integrals. Snapshots
    also retain the total cap integrals at their actual clock and radius.
    The profile is a subtraction choice, not another incoming boundary datum.
    """
    if not isinstance(budget, Budget):
        raise ValueError("a shared caller Budget is required")
    budget.check()
    x = _finite(x_nodes, "x_nodes")
    if (x.ndim != 1 or not 2 <= len(x) <= 2049 or np.any(np.diff(x) <= 0)
            or not np.allclose(np.diff(x), np.diff(x)[0], rtol=1e-10, atol=1e-14)):
        raise ValueError("use 1..2048 fixed uniform positive time steps")
    eta, y_faces = _grid(eta, "eta", 2, 512), _grid(y_faces, "y_faces", 3, 128)
    state = _state(initial, len(eta)-1, len(y_faces)-1)
    if type(save_every) is not int or save_every < 1 or (len(x)-1) % save_every:
        raise ValueError("save_every must divide the number of steps")
    if not math.isfinite(maximum_courant) or not 0 < maximum_courant <= .5:
        raise ValueError("require 0<maximum_courant<=.5")
    stages, max_courant, saved = 0, 0., []
    phase, step, attempt_state = "initial", 0, state

    def evaluate(value, at, dt):
        nonlocal stages, max_courant, attempt_state
        attempt_state = value
        result = background_stage(value, x=float(at), eta=eta, y_faces=y_faces,
            source=source(float(at)), gauge=gauge(float(at)), alpha=alpha, tube=tube,
            budget=budget, profile=profile, exterior_quadrature_order=exterior_quadrature_order,
            cap_quadrature_order=cap_quadrature_order)
        stages += 1
        max_courant = max(max_courant, dt*result["courant_rate"])
        if max_courant > maximum_courant:
            raise ValueError("declared characteristic Courant limit exceeded")
        return result

    def save(at, result):
        saved.append({"x": float(at), "state": copy.deepcopy(state),
                      "physical": {"increment": copy.deepcopy(result["physical_increment"]),
                                   "wall": copy.deepcopy(result["physical_wall"]),
                                   "cap": copy.deepcopy(result["physical_cap"])}})

    try:
        for step in range(len(x)-1):
            dt = float(x[step+1]-x[step])
            phase = "first_stage"
            first = evaluate(state, x[step], dt)
            if step % save_every == 0:
                save(x[step], first)
            predictor = _combine(state, first["rates"], 1., dt)
            phase = "predictor_stage"
            second = evaluate(predictor, x[step+1], dt)
            corrected = _combine(predictor, second["rates"], 1., dt)
            state = _combine(state, corrected, .5, .5)
        phase = "final_diagnostics"
        final = evaluate(state, x[-1], 0.)
        save(x[-1], final)
        budget.check()
    except (ValueError, RuntimeError, FloatingPointError, OverflowError) as error:
        raise CoupledEvolutionError({"phase": phase, "step": step,
            "x_left": float(x[step]), "x_right": float(x[min(step+1, len(x)-1)]),
            "reason": str(error), "completed_stages": stages,
            "charged_evaluations": budget.evaluations, "max_courant": max_courant,
            "accepted_state": copy.deepcopy(state), "attempt_state": copy.deepcopy(attempt_state)}) from error
    return {"snapshots": saved, "final_stage": final, "stages": stages,
            "steps": len(x)-1, "max_courant": max_courant,
            "time_scheme": "fixed-step SSPRK2", "retries": 0}
