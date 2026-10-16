"""Coupled stage with one explicit departure-gradient operator.

The physical cap, characteristic prefixes and hybrid wall equations remain
owned by their existing modules. Only coupling orchestration is retained
here so frozen historical adapters need not change. One returned gradient
feeds both the transverse prefix jets and the moving-grid rate.
"""
import math

import numpy as np

from .nsc_atlas_coupled_evolution import _finite, _grid
from .nsc_atlas_compatible_evolution import FIELDS, WALL_FIELDS, _state
from .nsc_atlas_nonlinear_boundary import hybrid_wall_boundary
from .nsc_atlas_vacuum_departure import vacuum_slice, departure_slice
from .nsc_atlas_background_cap import background_cap_metric, background_cap_rhs
from .nsc_atlas_departure_gradient import quadratic_departure_gradient
from .nsc_atlas_stage_ssprk3 import _shared_budget


def coupled_stage(state, *, x, eta, y_faces, source, gauge, alpha, tube,
                   budget, profile, departure_gradient=quadratic_departure_gradient,
                   exterior_quadrature_order=2, cap_quadrature_order=4):
    """Restore one physical state and evaluate all coupled rates.

    ``departure_gradient(eta, values, budget=budget)`` receives three rows
    ordered r, L, c, including both endpoint departures. It must return one
    finite array with the same shape and charge its work to the shared budget.
    Read-only callback inputs prevent accidental alteration of source or wall
    data. The returned gradient is not replaced by a required wall derivative.
    The Courant rate is a geometric speed/grid diagnostic, not a stability
    theorem for the selected spatial operator.
    """
    _shared_budget(budget)
    initial_evaluations = budget.evaluations
    if not callable(departure_gradient):
        raise ValueError("departure_gradient must be callable")
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

    gradient_eta, values = eta.copy(), np.stack([departures[key] for key in FIELDS])
    gradient_eta.flags.writeable = values.flags.writeable = False
    G = _finite(departure_gradient(gradient_eta, values, budget=budget),
                "departure gradient", values.shape)
    gradients = {key: G[i] for i, key in enumerate(FIELDS)}
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
        raise ValueError("nonfinite gradient coupled rate")
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
            "courant_rate": max(float(Wx/W/np.diff(eta)[0]), cap["characteristic_courant_rate"]),
            "departure_gradients": gradients,
            "charged_evaluations": budget.evaluations-initial_evaluations}
