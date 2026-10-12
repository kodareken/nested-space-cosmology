"""Regular null embedding of the atlas's fixed isolated thin-wall reference.

Published signed dynamics and Kruskal coordinates are integration inputs.
This adapter does not evolve the assembled many-body initial data.
"""
import math

import numpy as np

from .nsc_atlas_connection_slice import wall_slice_data, vacuum_transition, connection_slice_record

SCHEMA = "nsc-atlas-causal-reference-v1"


def signed_state(x):
    x = np.asarray(x, dtype=float)
    if not np.all(np.isfinite(x)) or np.any(x <= 0):
        raise ValueError("reference radius must be positive and finite")
    bm, bp = 1/(5*x*x)+3*x/4, 1/(5*x*x)-5*x/4
    q = 25*x*x/16-1+3/(10*x)+1/(25*x**4)
    s = np.sqrt(q)
    f = 1-.8/x
    # Rationalize only where subtraction would lose the small difference.
    difference = np.array(s-bp, copy=True)
    np.divide(-f, s+bp, out=difference, where=bp >= 0)
    vr = -1/(s*difference)
    return {"beta_minus": bm, "beta_plus": bp, "q": q, "Rdot": s,
            "f_plus": f, "advanced_derivative": vr}


def _integral(function, lower, upper, order=64):
    if order not in (32, 64):
        raise ValueError("the admitted quadrature orders are 32 and 64")
    if lower == upper:
        return 0.
    sign = 1 if upper > lower else -1
    lo, hi = sorted((float(lower), float(upper)))
    if lo < .02 or hi > 6:
        raise ValueError("numeric reference is restricted to 0.02 <= R/ell <= 6")
    cuts = [lo]+[v for v in (.2, .5, .8, 1., 1.5, 3.) if lo < v < hi]+[hi]
    nodes, weights = np.polynomial.legendre.leggauss(order)
    result = 0.
    for a, b in zip(cuts[:-1], cuts[1:]):
        grid = (a+b)/2+(b-a)/2*nodes
        result += (b-a)/2*float(np.dot(weights, function(grid)))
    return sign*result


def schwarzschild_slice_time(x, order=64):
    """T/ell anchored to the v4 K=0 bridge at T=0, for its left transition."""
    if not math.isfinite(x) or not 1 <= x <= 1.5:
        raise ValueError("left slice transition has 1 <= R/ell <= 1.5")
    wall = wall_slice_data()
    def gradient(radius):
        fields = vacuum_transition(radius, wall)
        return fields["V"]/(1-.8/radius)/fields["Rprime"]
    return _integral(gradient, 1., x, order)


def reference_wall_events(radii_over_ell, order=64):
    x = np.asarray(radii_over_ell, dtype=float)
    if x.ndim != 1 or len(x) == 0 or len(x) > 300 or not np.all(np.isfinite(x)) or np.any(x < .02) or np.any(x > 6):
        raise ValueError("admitted finite wall grid: 1..300 radii in [0.02,6]")
    M = .4
    anchor_T = schwarzschild_slice_time(1.5, order)
    anchor_v = anchor_T+1.5+.8*math.log(1.5/.8-1)
    v = np.array([anchor_v+_integral(lambda y: signed_state(y)["advanced_derivative"], 1.5, y, order) for y in x])
    tau = np.array([_integral(lambda y: 1/signed_state(y)["Rdot"], 1.5, y, order) for y in x])
    state = signed_state(x)
    V = -np.exp(v/(4*M))
    W = (1-x/(2*M))*np.exp(x/(2*M))
    U = W/V
    V_R = V/(4*M)*state["advanced_derivative"]
    W_R = -x/(4*M*M)*np.exp(x/(2*M))
    U_R = W_R/V-W*V_R/(V*V)
    dU, dV = U_R*state["Rdot"], V_R*state["Rdot"]
    norm = -32*M**3/x*np.exp(-x/(2*M))*dU*dV
    return {"R_over_ell": x, "advanced_over_ell": v, "tau_wall_over_ell": tau,
            "U": U, "V": V, "X": (np.arctan(V)-np.arctan(U))/2,
            "T": (np.arctan(V)+np.arctan(U))/2,
            "ell_dU_dtau": dU, "ell_dV_dtau": dV, "metric_norm": norm,
            "beta_minus": state["beta_minus"], "beta_plus": state["beta_plus"]}


def parent_event_from_U(U, radius_over_ell=1.6):
    """Intersection of an outgoing null ray with a right static observer."""
    if not math.isfinite(U) or U >= 0 or not math.isfinite(radius_over_ell) or radius_over_ell <= .8:
        raise ValueError("right-parent null access requires U < 0 and R > 2M")
    R = radius_over_ell
    V = (1-R/.8)*math.exp(R/.8)/U
    t = .8*math.log(V/(-U))
    return {"U": U, "V": V, "X": (math.atan(V)-math.atan(U))/2,
            "T": (math.atan(V)+math.atan(U))/2,
            "schwarzschild_T_over_ell": t, "parent_proper_time_over_ell": math.sqrt(1-.8/R)*t}


def initial_slice_curve():
    c = connection_slice_record()
    rows = []
    for R in c["transition"]["samples"]["R_over_ell"]:
        T = schwarzschild_slice_time(R)
        factor = math.sqrt((R/.8-1)*math.exp(R/.8))
        rows.append((factor*math.exp(-T/1.6), -factor*math.exp(T/1.6)))
    for R, u in zip(c["bridge"]["samples"]["R_over_ell"], c["bridge"]["samples"]["u"]):
        factor = math.sqrt(max(0., (R/.8-1)*math.exp(R/.8)))
        rows.append((factor, -factor) if u <= 0 else (-factor, factor))
    U, V = np.asarray(rows).T
    return {"U": U.tolist(), "V": V.tolist(), "X": ((np.arctan(V)-np.arctan(U))/2).tolist(),
            "T": ((np.arctan(V)+np.arctan(U))/2).tolist()}


def causal_reference_record():
    x = np.unique(np.concatenate((np.linspace(.02, .8, 81), np.linspace(.8, 1.5, 65), np.linspace(1.5, 6., 81), [.5, 1.])))
    wall = reference_wall_events(x)
    independent = reference_wall_events(np.array([.02, .5, .8, 1., 1.5, 3., 6.]), 32)
    control = reference_wall_events(independent["R_over_ell"], 64)
    early, now = (reference_wall_events(np.array([r])) for r in (.5, 1.5))
    arrival = parent_event_from_U(float(early["U"][0]))
    # The unique positive stationary point is analytic, not found by a sweep.
    y = (.3+math.sqrt(.09+2.))/(25/4)
    xmin = y**(1/3)
    return {"schema": SCHEMA, "scope": "isolated spherical reference completion only",
        "status": "MONOTONIC_EXPANDING_REFERENCE_WITH_SINGULAR_PAST",
        "inputs": {"H_ell": 1., "k_ell": 2., "M_over_ell": .4, "R_wall_cut_over_ell": 1.5,
                   "parent_observer_R_over_ell": 1.6, "Schwarzschild_T_origin": "v4 time-symmetric bridge"},
        "admission": {"owner": "docs/nsc-cosmological-atlas-v5.md",
            "borrowed": ["gr-qc/0508093v2 SectionII signed classification", "gr-qc/9403003v1 Section2.1 Eq9-13"],
            "missing": "our fixed initial slice and observer events within the known reference history",
            "decides": "regular null paths and reference causal exclusions",
            "unlocks": "observer access table and choice of response problem",
            "finite_domain": {"R_over_ell": [.02, 6.], "samples": len(x)},
            "quadrature_orders": [32, 64], "cpu_seconds_max": 60,
            "stop": "wrong signed branch, nonfuture null coordinates or failed normalization",
            "priority_claim": False},
        "classification": {"q_has_unique_positive_minimum": True, "q_minimum_at_R_over_ell": xmin,
            "published_branch": "gr-qc/0508093v2 Figure10 solution13, Lambda_plus->0 limit",
            "published_parameters": {"A": 0., "B": .75, "Y": .6},
            "limit_changes": "outer Schwarzschild-deSitter cosmological horizon moves to infinity; exterior sequence becomes IV->III",
            "q_minimum": float(signed_state(xmin)["q"]),
            "all_radii_positivity_proof": "x>=0.8:25*x^2/16>=1; x<=0.8:AM-GM gives25*x^2/16+1/(25*x^4)>3/4 and3/(10*x)>=3/8, so q>1/8",
            "beta_plus_zero_over_ell": (.16)**(1/3), "beta_minus_positive_all_radii": True,
            "Schwarzschild_horizon_over_ell": .8, "de_Sitter_horizon_over_ell": 1.,
            "quadrants": "IV (U<0,V<0) for R<0.8ell; III (U>0,V<0) for R>0.8ell",
            "past": "R->0 in finite wall proper time; retained Schwarzschild curvature48*M^2/R^6 diverges",
            "future": "R->infinity at infinite wall proper time; interior expanding de Sitter reference",
            "future_V_sign": "negative for every finite R; late v derivative is O(R^-2), giving a finite negative limiting V",
            "regular_parent_formation_established": False},
        "wall_samples": {key: value.tolist() for key, value in wall.items()},
        "initial_slice": initial_slice_curve(),
        "events": {"early_wall": {key: float(value[0]) for key, value in early.items()},
            "cut_wall": {key: float(value[0]) for key, value in now.items()}, "early_wall_parent_arrival": arrival},
        "access": [
            {"from": "W_early", "to": "P", "verdict": "allowed_reference_null_path",
             "emission": "wallR/ell=0.5 before the cut", "basis": "constantU<0 from IV to I; V increases",
             "field_transmission_amplitude": "not specified"},
            {"from": "P", "to": "W", "verdict": "excluded_in_reference",
             "emission": "any event on right static parent worldline", "basis": "futureV cannot decrease from V>0 to wallV<0"},
            {"from": "W_cut_or_later", "to": "P", "verdict": "excluded_in_reference",
             "emission": "wallR>=1.5ell", "basis": "futureU cannot decrease from U>0 to parentU<0"}],
        "checks": {"wall_metric_norm_max_error": float(np.max(abs(wall["metric_norm"]+1))),
            "minimum_ell_dU_dtau": float(np.min(wall["ell_dU_dtau"])),
            "minimum_ell_dV_dtau": float(np.min(wall["ell_dV_dtau"])),
            "quadrature_32_64_coordinate_max_difference": max(float(np.max(abs(independent[key]-control[key]))) for key in ("X", "T", "tau_wall_over_ell")),
            "parent_null_U_error": abs(arrival["U"]-float(early["U"][0])),
            "comparison_is_rigorous_bound": False},
        "limits": ["known isolated reference, not the glued many-body future", "finite plots do not prove infinite-time limits",
            "no regular formation or perturbative stability proved", "allowed causal paths do not determine signal amplitudes",
            "clock origins belong to named worldlines, not a common cosmic age"]}
