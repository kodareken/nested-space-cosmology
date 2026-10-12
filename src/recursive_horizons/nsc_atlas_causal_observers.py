"""Event-specific de Sitter access for the fixed isolated vacuum-wall reference.

The complete-future exclusions use analytic inequalities, not the plot cutoff.
The right-parent Schwarzschild decisions belong to the companion reference.
"""

import math

import numpy as np


SCHEMA = "nsc-atlas-causal-observers-v1"
MIN_RADIUS = .05
MAX_RADIUS = 3.
ANCHOR_RADIUS = 1.5


def _scale(length_scale):
    ell = float(length_scale)
    if not math.isfinite(ell) or ell <= 0:
        raise ValueError("length scale must be finite and positive")
    return ell


def de_sitter_wall_state(radii, length_scale=1.):
    """Regular fixed-source kinematics, including both static horizons."""
    ell = _scale(length_scale)
    R = np.asarray(radii, dtype=float)
    x = R/ell
    if not np.all(np.isfinite(x)) or np.any(x < MIN_RADIUS) or np.any(x > MAX_RADIUS):
        raise ValueError("radius outside retained finite reference interval")
    beta_minus = 1/(5*x*x)+3*x/4
    beta_plus = 1/(5*x*x)-5*x/4
    q = 25*x*x/16-1+3/(10*x)+1/(25*x**4)
    if np.any(q <= 0):
        raise ValueError("non-real or turning wall motion")
    speed = np.sqrt(q)
    z = (beta_minus+speed)/(1+x)
    gamma, gamma_v = (z+1/z)/2, (z-1/z)/2
    gamma_other = (1+speed*speed)/(beta_minus+x*speed)
    return {"R": R, "Rdot": speed, "beta_minus": beta_minus, "beta_plus": beta_plus,
            "gamma": gamma, "gamma_v": gamma_v, "v": gamma_v/gamma,
            "gamma_rationalized": gamma_other,
            "f_minus": 1-x*x, "f_plus": 1-.8/x}


def de_sitter_wall_events(radii, length_scale=1., order=32):
    """Two anchored clocks and conformal events on the recorded radius interval.

    The returned arrays follow the input order. Quadrature knots include the
    clock anchor; no singular Schwarzschild/de Sitter static time is used.
    """
    ell = _scale(length_scale)
    R = np.asarray(radii, dtype=float)
    if R.ndim != 1 or len(R) == 0:
        raise ValueError("radii must be a nonempty vector")
    if not isinstance(order, int) or not 4 <= order <= 64:
        raise ValueError("quadrature order must be an integer from 4 to 64")
    state = de_sitter_wall_state(R, ell)
    # Fixed intermediate knots bound individual integration spans even when
    # a caller requests only one event far from the anchor.
    base = ell*np.linspace(MIN_RADIUS, MAX_RADIUS, 117)
    knots = np.unique(np.concatenate((base, R, [ANCHOR_RADIUS*ell])))
    nodes, weights = np.polynomial.legendre.leggauss(order)
    clocks = np.zeros((len(knots), 2))
    for index, (lower, upper) in enumerate(zip(knots[:-1], knots[1:])):
        probe = (upper+lower)/2+(upper-lower)/2*nodes
        data = de_sitter_wall_state(probe, ell)
        factor = (upper-lower)/2
        clocks[index+1] = clocks[index]+factor*np.array([
            np.dot(weights, data["gamma"]/data["Rdot"]),
            np.dot(weights, 1/data["Rdot"]),
        ])
    anchor = int(np.searchsorted(knots, ANCHOR_RADIUS*ell))
    clocks -= clocks[anchor]
    selected = np.searchsorted(knots, R)
    t, tau = clocks[selected, 0], clocks[selected, 1]
    conformal_factor = np.exp(-t/ell)
    return dict(state, t_flat=t, tau_wall=tau, eta=-ell*conformal_factor,
                r_wall=R*conformal_factor,
                quadrature_evaluations=(len(knots)-1)*order)


def wall_to_centre_intercept(radius, flat_time, length_scale=1.):
    """Classify an inward null ray in expanding de Sitter, not a material field."""
    ell = _scale(length_scale)
    x = float(radius)/ell
    if not math.isfinite(x) or x <= 0 or not math.isfinite(flat_time):
        raise ValueError("wall event must have finite positive radius and finite flat time")
    eta = -ell*math.exp(-flat_time/ell)
    r = -x*eta
    if x < 1:
        return {"status": "FINITE_CENTRE_RECEPTION", "centre_flat_time": flat_time-ell*math.log1p(-x),
                "eta_intercept": eta+r}
    if x == 1:
        return {"status": "ONLY_AT_INFINITE_CENTRE_TIME", "centre_flat_time": None,
                "eta_intercept": 0.}
    return {"status": "BEYOND_FUTURE_CONFORMAL_BOUNDARY", "centre_flat_time": None,
            "eta_intercept": eta+r}


def centre_emission_to_wall(radius, flat_time, length_scale=1.):
    """Outgoing null intercept; matched-history coverage is checked separately."""
    ell = _scale(length_scale)
    x = float(radius)/ell
    if not math.isfinite(x) or x <= 0 or not math.isfinite(flat_time):
        raise ValueError("wall event must have finite positive radius and finite flat time")
    centre_time = flat_time-ell*math.log1p(x)
    return {"centre_flat_time": centre_time,
            "eta_intercept": -ell*math.exp(-centre_time/ell)}


def require_retained_time(time, earliest, latest):
    if not all(math.isfinite(x) for x in (time, earliest, latest)) or not earliest <= time <= latest:
        raise ValueError("event lies outside retained matched-history time interval")


def causal_observers_record(length_scale=1.):
    ell = _scale(length_scale)
    x = np.unique(np.concatenate((np.linspace(MIN_RADIUS, MAX_RADIUS, 117), [.5, .8, 1., 1.5, 1.8])))
    data = de_sitter_wall_events(ell*x, ell, 32)
    check = de_sitter_wall_events(ell*x, ell, 64)
    tmin, tmax = float(data["t_flat"][0]), float(data["t_flat"][-1])

    def event(radius):
        index = int(np.where(x == radius)[0][0])
        return {"R_over_ell": radius, "t_flat_over_ell": float(data["t_flat"][index]/ell),
                "tau_wall_over_ell": float(data["tau_wall"][index]/ell),
                "eta_over_ell": float(data["eta"][index]/ell),
                "r_wall_over_ell": float(data["r_wall"][index]/ell),
                "v": float(data["v"][index])}

    early, now, horizon = event(.5), event(1.5), event(1.)
    schwarzschild_horizon = event(.8)
    early_reception = wall_to_centre_intercept(.5*ell, early["t_flat_over_ell"]*ell, ell)
    historical_emission = centre_emission_to_wall(1.5*ell, 0., ell)
    early_emission = centre_emission_to_wall(.5*ell, early["t_flat_over_ell"]*ell, ell)
    escape_emission = centre_emission_to_wall(.8*ell, schwarzschild_horizon["t_flat_over_ell"]*ell, ell)
    # Both ends and the full monotone time intervals lie within the retained
    # matched history; timelike characteristic-intercept monotonicity then
    # keeps the recorded radial rays on the interior side until reception.
    require_retained_time(early_reception["centre_flat_time"], tmin, tmax)
    require_retained_time(historical_emission["centre_flat_time"], tmin, tmax)
    require_retained_time(early_emission["centre_flat_time"], tmin, tmax)
    require_retained_time(escape_emission["centre_flat_time"], tmin, tmax)
    eta_c_from_early = -ell*math.exp(-early_reception["centre_flat_time"]/ell)
    eta_c_to_now = historical_emission["eta_intercept"]
    samples = {"R_over_ell": x.tolist()}
    for source, target in (("t_flat", "t_flat_over_ell"), ("tau_wall", "tau_wall_over_ell"),
                           ("eta", "eta_over_ell"), ("r_wall", "r_wall_over_ell")):
        samples[target] = (data[source]/ell).tolist()
    for source, target in (("v", "velocity_relative_flat"), ("Rdot", "Rdot"),
                           ("gamma", "gamma"), ("beta_minus", "beta_minus"), ("beta_plus", "beta_plus")):
        samples[target] = data[source].tolist()
    outgoing_intercept = data["eta"]-data["r_wall"]
    ingoing_intercept = data["eta"]+data["r_wall"]
    return {
        "schema": SCHEMA, "status": "EVENT_SPECIFIC_REFERENCE_ACCESS",
        "admission": {"owner": "docs/nsc-atlas-causal-observers.md",
            "borrowed": ["0711.2198v1 Eq2", "gr-qc/0508093v2 SectionII", "hep-th/0110007v2 Section2.1"],
            "missing": "centre/wall emission and reception map for the fixed isolated reference",
            "decision": "historical light paths coexist with future observer disconnection",
            "unlocks": "event-specific causal labels and later field/state questions",
            "cpu_seconds_max": 60,
            "quadrature_evaluations": data["quadrature_evaluations"]+check["quadrature_evaluations"],
            "stop": "failed regular motion, null equality or retained-history coverage",
            "priority_claim": False},
        "theory": {"H_ell": 1., "k_ell": 2., "M_over_ell": .4, "fixed_within_test": True,
                   "reference": "isolated expanding de Sitter/Schwarzschild pure-wall branch"},
        "observers": [
            {"id": "C", "worldline": "de Sitter flat-chart centre r=0", "clock": "proper time t_flat"},
            {"id": "W", "worldline": "fixed signed expanding shell", "clock": "proper shell time tau_wall"},
            {"id": "P", "worldline": "right Schwarzschild static R=4M", "clock": "separate parent proper time",
             "access_owner": "companion Schwarzschild causal-reference record"},
        ],
        "conventions": {"time_origin": "W0 at R/ell=1.5 has t_flat=tau_wall=0 and scale factor1",
                        "conformal_time": "eta=-ell exp(-t_flat/ell), eta<0",
                        "flat_chart": "expanding planar chart; time reverse of hep-th/0110007v2 Eq21",
                        "null_paths": "radial eta+r=constant inward; eta-r=constant outward",
                        "clock_origins_do_not_establish_universe_ages": True},
        "trajectory": {"R_over_ell_domain": [MIN_RADIUS, MAX_RADIUS],
                       "t_flat_over_ell_domain": [tmin/ell, tmax/ell],
                       "finite_plot_is_complete_future": False, "samples": samples},
        "named_events": {"W_early": early, "W0": now, "W_dS_horizon": horizon,
                         "W_Schwarzschild_horizon": schwarzschild_horizon,
                         "C0": {"t_flat_over_ell": 0., "eta_over_ell": -1., "r_over_ell": 0.},
                         "C_emits_to_W0": {"t_flat_over_ell": historical_emission["centre_flat_time"]/ell,
                                           "eta_over_ell": eta_c_to_now/ell, "r_over_ell": 0.},
                         "C_emits_to_early_wall": {"t_flat_over_ell": early_emission["centre_flat_time"]/ell,
                                                   "eta_over_ell": early_emission["eta_intercept"]/ell,
                                                   "r_over_ell": 0.},
                         "C_escape_threshold": {"t_flat_over_ell": escape_emission["centre_flat_time"]/ell,
                                                  "eta_over_ell": escape_emission["eta_intercept"]/ell,
                                                  "r_over_ell": 0.},
                         "C_receives_early_wall": {"t_flat_over_ell": early_reception["centre_flat_time"]/ell,
                                                   "eta_over_ell": eta_c_from_early/ell, "r_over_ell": 0.}},
        "access": [
            {"from": "C_emits_to_early_wall", "to": "W_early", "status": "CERTIFIED_RADIAL_NULL_PATH",
             "scope": "finite isolated matched-reference history", "direction": "outward",
             "proof": "emission is after retained start; increasing wall eta-r gives first reception"},
            {"from": "W_early", "to": "C_receives_early_wall", "status": "CERTIFIED_RADIAL_NULL_PATH",
             "scope": "finite isolated matched-reference history", "direction": "inward",
             "proof": "eta+r negative and constant; timelike wall has strictly increasing eta+r"},
            {"from": "C_emits_to_W0", "to": "W0", "status": "CERTIFIED_RADIAL_NULL_PATH",
             "scope": "finite isolated matched-reference history", "direction": "outward",
             "proof": "eta-r constant; strictly increasing wall eta-r makes this its first reception"},
            {"from": "W0", "to": "future_C", "status": "NO_FUTURE_CAUSAL_ACCESS",
             "scope": "isolated expanding de Sitter completion", "proof": "eta_W0+r_W0=.5ell>0"},
            {"from": "C0", "to": "future_W", "status": "NO_FUTURE_CAUSAL_ACCESS",
             "scope": "entire isolated fixed expanding future branch",
             "proof": "future wall comoving radius>=1.5ell; centre future light cone radius<ell"},
            {"from": "future_W_after_W0", "to": "future_C", "status": "NO_FUTURE_CAUSAL_ACCESS",
             "scope": "entire isolated fixed expanding future branch",
             "proof": "eta_W+r_W=(R/ell-1)*ell*exp(-t/ell)>0 for R>=1.5ell"},
        ],
        "future_proof": {"analytic_not_sample_inference": True,
            "domain": "R/ell>=1.5",
            "beta_minus_lower_bound": 1/(5*1.5**2)+3*1.5/4,
            "beta_minus_derivative_positive": "3/4-2/(5x^3)>0 for x>=1.5",
            "velocity_sign_identity": "gamma v=(beta_minus^2-1)/(Rdot+x beta_minus)>0",
            "inward_intercept": "eta+r=ell exp(-t/ell)(x-1)>0",
            "wall_comoving_radius_lower_bound_over_ell": 1.5,
            "C0_remaining_null_radius_over_ell": 1.,
            "complete_future_wall_time": "q>0 and Rdot<2x for x>=1.5 imply integral dR/Rdot diverges",
            "no_past_completeness_claim": True},
        "parent_escape_threshold": {
            "status": "CENTRAL_NULL_RAY_TO_SCHWARZSCHILD_HORIZON",
            "centre_emission_t_over_ell": escape_emission["centre_flat_time"]/ell,
            "eta_minus_r_over_ell": escape_emission["eta_intercept"]/ell,
            "wall_crossing_R_over_ell": .8,
            "matched_path_retained": True,
            "parent_interpretation_requires_companion_reference": True,
            "companion_condition": "wall U_K<0 before R=.8ell, U_K=0 at crossing, U_K>0 afterwards",
            "strict_early_emission_interval_t_over_ell": [tmin/ell, escape_emission["centre_flat_time"]/ell],
            "upper_endpoint_admits_finite_parent_reception": False,
            "scope": "isolated Schwarzschild completion; not the glued many-body future",
            "radial_null_segment": "r=eta-eta_escape from centre to wall horizon",
        },
        "independent_checks": {
            "clock_32_64_max_error_over_ell": max(float(np.max(abs(data[key]-check[key])))/ell for key in ("t_flat", "tau_wall")),
            "regular_boost_formula_max_error": float(np.max(abs(data["gamma"]-data["gamma_rationalized"]))),
            "boost_unit_max_error": float(np.max(abs(data["gamma"]**2-data["gamma_v"]**2-1))),
            "radial_tangent_max_error": float(np.max(abs(data["Rdot"]-data["gamma"]*x-data["gamma_v"]))),
            "radial_normal_max_error": float(np.max(abs(data["beta_minus"]-data["gamma"]-data["gamma_v"]*x))),
            "early_inward_null_error_over_ell": abs(eta_c_from_early/ell-early["eta_over_ell"]-early["r_wall_over_ell"]),
            "historical_outward_null_error_over_ell": abs(eta_c_to_now/ell-now["eta_over_ell"]+now["r_wall_over_ell"]),
            "early_outward_null_error_over_ell": abs(early_emission["eta_intercept"]/ell-early["eta_over_ell"]+early["r_wall_over_ell"]),
            "escape_outward_null_error_over_ell": abs(escape_emission["eta_intercept"]/ell-schwarzschild_horizon["eta_over_ell"]+schwarzschild_horizon["r_wall_over_ell"]),
            "historical_emission_after_retained_start_over_ell": (historical_emission["centre_flat_time"]-tmin)/ell,
            "early_emission_after_retained_start_over_ell": (early_emission["centre_flat_time"]-tmin)/ell,
            "escape_emission_after_retained_start_over_ell": (escape_emission["centre_flat_time"]-tmin)/ell,
            "early_reception_before_retained_end_over_ell": (tmax-early_reception["centre_flat_time"])/ell,
            "minimum_outgoing_intercept_increment_over_ell": float(np.min(np.diff(outgoing_intercept)))/ell,
            "minimum_ingoing_intercept_increment_over_ell": float(np.min(np.diff(ingoing_intercept)))/ell,
            "quadrature_comparison_is_rigorous_error_bound": False,
        },
        "source_reading": [
            {"id": "hep-th/0110007v2", "url": "https://arxiv.org/html/hep-th/0110007v2",
             "depth": "Section2.1 Eq15-21 and observer causal-diamond discussion",
             "use": "planar chart with explicit time reversal; conformal null propagation"},
            {"id": "0711.2198v1", "url": "https://arxiv.org/html/0711.2198v1",
             "depth": "Section2 Eq2 and Section3, reused v3 reading", "use": "fixed signed shell dynamics"},
            {"id": "gr-qc/0508093v2", "url": "https://arxiv.org/html/gr-qc/0508093v2",
             "depth": "SectionII patch and orientation discussion, reused v3 reading",
             "use": "isolated-reference causal regions; no stable formation claim"},
        ],
        "limits": ["isolated reference; not a future solution for the glued many-body parent",
                   "past white-hole completion is not regular formation", "test-light causal paths do not specify field transmission through the wall",
                   "no new clock-based age, inherited constant or sibling-signal prediction",
                   "no complete-future conclusion inferred from plotting samples"],
    }
