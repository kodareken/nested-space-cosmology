"""Exact conditional wall-to-centre transfer for one conformal test scalar.

The exterior-side incoming characteristic is prescribed. No exterior source,
scattering solution, energy flux or backreaction is inferred from this map.
"""

from fractions import Fraction
import math
import time

import numpy as np

from .nsc_atlas_causal_observers import de_sitter_wall_events

SCHEMA = "nsc-atlas-field-transfer-v1"
PULSE_HALF_WIDTH = .03
CENTRE_WALL_RADIUS = .5
INITIAL_WALL_RADIUS = .05
END_WALL_RADIUS = 1.5
MAX_WALL_SAMPLES = 513


def _scale(length_scale):
    ell = float(length_scale)
    if not math.isfinite(ell) or ell <= 0:
        raise ValueError("length scale must be finite and positive")
    return ell


def compact_bump(z):
    """C3 compact primitive B and C2 derivative, including exact edge zeros."""
    z = np.asarray(z, dtype=float)
    if not np.all(np.isfinite(z)):
        raise ValueError("pulse coordinate must be finite")
    active = abs(z) < 1
    safe = np.where(active, z, 0.)
    b = np.where(active, (1-safe*safe)**4, 0.)
    derivative = np.where(active, -8*safe*(1-safe*safe)**3, 0.)
    return b, derivative


def conformal_junction_coefficients(k, xi=1/6):
    """Coefficients multiplying phi or chi, with n from cap to exterior."""
    k, xi = float(k), float(xi)
    if not math.isfinite(k) or k <= 0 or not math.isfinite(xi):
        raise ValueError("finite positive wall k and finite coupling required")
    return {"ricci_delta": 6*k, "normal_phi_jump_per_phi": 6*xi*k,
            "normal_R_jump_per_R": -k, "normal_chi_jump_per_chi": (6*xi-1)*k}


def arrival_from_wall_event(radius, flat_time, gamma, velocity, length_scale=1.):
    """Finite proper-time arrival and its positive derivative on R<ell."""
    ell = _scale(length_scale)
    x, t, gamma, velocity = (np.asarray(value, dtype=float) for value in
                              (np.asarray(radius)/ell, flat_time, gamma, velocity))
    if any(not np.all(np.isfinite(value)) for value in (x, t, gamma, velocity)):
        raise ValueError("wall event must be finite")
    if np.any(x <= 0) or np.any(x >= 1):
        raise ValueError("finite central arrival requires 0<R/ell<1")
    if np.any(gamma < 1) or np.any(abs(velocity) >= 1):
        raise ValueError("wall event must have a future timelike boost")
    return t-ell*np.log1p(-x), gamma*(1+velocity)/(1-x)


def regular_radial_field(eta, radius, primitive, derivative, length_scale=1.):
    """Exact interior solution at comoving radius r; supplied F is dimensionless."""
    ell = _scale(length_scale)
    eta, radius = float(eta), float(radius)
    if not math.isfinite(eta) or eta >= 0 or not math.isfinite(radius) or radius < 0:
        raise ValueError("expanding conformal chart needs eta<0 and radius>=0")
    scale_factor = -ell/eta
    if radius == 0:
        return 2*float(derivative(eta))/scale_factor
    return (float(primitive(eta+radius))-float(primitive(eta-radius)))/(scale_factor*radius)


def _input_squared_norm():
    # Integral[-1,1] 64*z^2*(1-z^2)^6 dz, evaluated exactly as a rational.
    return 128*sum((Fraction((-1)**j*math.comb(6, j), 2*j+3) for j in range(7)), Fraction())


def field_transfer_record(length_scale=1.):
    """Normalized pulse, exact response map and bounded numerical diagnostics."""
    ell = _scale(length_scale)
    started = time.process_time()
    sampled_radii = set()
    quadrature_evaluations = 0

    def clock(x, order=32):
        nonlocal quadrature_evaluations
        x = np.atleast_1d(np.asarray(x, dtype=float))
        sampled_radii.update(float(value) for value in x)
        if len(sampled_radii) > MAX_WALL_SAMPLES:
            raise RuntimeError("fixed wall sample budget exhausted")
        if time.process_time()-started > 60:
            raise RuntimeError("fixed scalar transfer CPU budget exhausted")
        result = de_sitter_wall_events(ell*x, ell, order)
        quadrature_evaluations += result["quadrature_evaluations"]
        return result

    anchors = clock([INITIAL_WALL_RADIUS, CENTRE_WALL_RADIUS, .8, END_WALL_RADIUS])
    tau_c = float(anchors["tau_wall"][1]/ell)
    delta = PULSE_HALF_WIDTH

    def profile(x, order=32):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        event = clock(x, order)
        tau = event["tau_wall"]/ell
        z = (tau-tau_c)/delta
        bump, bump_prime = compact_bump(z)
        arrival, stretch = arrival_from_wall_event(ell*x, event["t_flat"],
                                                   event["gamma"], event["v"], ell)
        incident = bump_prime/delta
        return {"R_over_ell": x, "tau_wall_over_ell": tau,
                "t_wall_flat_over_ell": event["t_flat"]/ell,
                "arrival_t_C_over_ell": arrival/ell,
                "incoming_null_label_over_ell": (event["eta"]+event["r_wall"])/ell,
                "z": z, "primitive_F_over_epsilon": bump/2,
                "ell_I_over_epsilon": incident,
                "ell_phi_C_over_epsilon": incident/stretch,
                "dT_dtau": stretch, "Rdot": event["Rdot"],
                "gamma": event["gamma"], "v_peculiar": event["v"],
                "a_wall": np.exp(event["t_flat"]/ell)}

    # Two finite endpoint solves at fixed source; no parameter search.
    targets = tau_c+np.array([-delta, delta])
    low, high = np.array([.4, .5]), np.array([.5, .6])
    if np.any(clock(low)["tau_wall"]/ell > targets) or np.any(clock(high)["tau_wall"]/ell < targets):
        raise ValueError("fixed pulse support is outside the admitted radius brackets")
    for _ in range(40):
        middle = (low+high)/2
        before = clock(middle)["tau_wall"]/ell < targets
        low, high = np.where(before, middle, low), np.where(before, high, middle)
    support_radii = (low+high)/2
    support = profile(support_radii)
    sample_radii = np.unique(np.concatenate((np.linspace(*support_radii, 129), [.5])))
    samples = profile(sample_radii)
    higher_order = profile(sample_radii, 64)

    def extremum(lower, upper):
        ratio = (math.sqrt(5)-1)/2
        left, right = upper-ratio*(upper-lower), lower+ratio*(upper-lower)
        fleft = abs(float(profile([left])["ell_phi_C_over_epsilon"][0]))
        fright = abs(float(profile([right])["ell_phi_C_over_epsilon"][0]))
        for _ in range(32):
            if fleft < fright:
                lower, left, fleft = left, right, fright
                right = lower+ratio*(upper-lower)
                fright = abs(float(profile([right])["ell_phi_C_over_epsilon"][0]))
            else:
                upper, right, fright = right, left, fleft
                left = upper-ratio*(upper-lower)
                fleft = abs(float(profile([left])["ell_phi_C_over_epsilon"][0]))
        data = profile([(lower+upper)/2])
        return {key: float(value[0]) for key, value in data.items()}

    positive = extremum(support_radii[0], .5)
    negative = extremum(.5, support_radii[1])

    integrals = []
    for order in (32, 64):
        node, weight = np.polynomial.legendre.leggauss(order)
        width = float(support_radii[1]-support_radii[0])/2
        probe = float(np.mean(support_radii))+width*node
        data = profile(probe, order)
        incident, stretch, speed = data["ell_I_over_epsilon"], data["dT_dtau"], data["Rdot"]
        integrals.append({"input_squared": width*float(np.dot(weight, incident**2/speed)),
                          "output_squared": width*float(np.dot(weight, incident**2/(stretch*speed))),
                          "signed_area": width*float(np.dot(weight, incident/speed))})
    exact_input_norm = float(_input_squared_norm())/delta

    # An independent symmetric radial limit uses the actual composed F(v),
    # found by two bracketed null-label inversions, not the transfer formula.
    centre_eta = positive["incoming_null_label_over_ell"]
    radial_probe = 1e-5
    desired_labels = centre_eta+np.array([-radial_probe, radial_probe])
    low = np.full(2, support_radii[0])
    high = np.full(2, support_radii[1])
    for _ in range(38):
        middle = (low+high)/2
        data = clock(middle)
        labels = (data["eta"]+data["r_wall"])/ell
        before = labels < desired_labels
        low, high = np.where(before, middle, low), np.where(before, high, middle)
    primitive = profile((low+high)/2)["primitive_F_over_epsilon"]
    radial_limit = (primitive[1]-primitive[0])/((-1/centre_eta)*radial_probe)

    initial_eta = float(anchors["eta"][0]/ell)
    initial_vmax = float((anchors["eta"][0]+anchors["r_wall"][0])/ell)
    support_v = support["incoming_null_label_over_ell"]
    arrival_times = support["arrival_t_C_over_ell"]
    horizon_wall_eta = float(anchors["eta"][2]/ell)
    reflected_radius_bound = -1-float(support_v[0])
    initial_gap = float(support_v[0])-initial_vmax
    after_horizon_gap = float(support_v[0])-horizon_wall_eta
    no_return_gap = .8-reflected_radius_bound
    if initial_gap <= 0 or after_horizon_gap <= 0 or no_return_gap <= 0 or arrival_times[1] >= 0:
        raise ValueError("fixed pulse leaves the certified initial or no-return window")

    junction = conformal_junction_coefficients(2/ell)
    sigma = (2/ell)/(4*math.pi)
    trace_ricci = -8*math.pi*(-3*sigma)
    positive_peak = positive["ell_phi_C_over_epsilon"]
    negative_peak = negative["ell_phi_C_over_epsilon"]
    input_peak = (8/math.sqrt(7))*(6/7)**3/delta
    zero_checks = {"initial_field": 0., "initial_time_derivative": 0.,
                   "central_field_after_last_arrival_through_t0": 0.,
                   "outgoing_wall_characteristic_before_t0": 0.}
    clock_error = max(float(np.max(abs(samples[key]-higher_order[key])))
                      for key in ("tau_wall_over_ell", "t_wall_flat_over_ell", "arrival_t_C_over_ell"))
    return {
        "schema": SCHEMA, "status": "CONDITIONAL_EARLY_CLASSICAL_SCALAR_RESPONSE",
        "admission": {"owner": "docs/nsc-atlas-field-transfer.md", "cpu_seconds_max": 60,
            "distinct_wall_radius_samples": len(sampled_radii), "wall_samples_max": MAX_WALL_SAMPLES,
            "clock_quadrature_evaluations": quadrature_evaluations,
            "borrowed": ["v5 fixed signed wall and clocks", "1512.02120v2 Eq8,26-34", "gr-qc/0603031v1 conformal map"],
            "missing_closed": "field, incoming state and exact conditional centre response for an early allowed path",
            "unlocks": "compare physically supplied incoming states or alternative field theories",
            "stop": "failed support, junction, time orientation or no-return certificate", "priority_claim": False},
        "theory": {"dimension": 4, "field": "classical real massless scalar", "angular_l": 0,
            "curvature_coupling_xi": 1/6, "additional_scalar_surface_coupling": 0.,
            "H_ell": 1., "k_ell": 2., "M_over_ell": .4, "fixed_within_test": True,
            "amplitude": "formal epsilon; field O(epsilon), backreaction omitted at O(epsilon^2)"},
        "incoming_state": {"side": "exterior limit at the moving wall",
            "characteristic": "I=(d_tau_wall+n)(R phi), n from cap to exterior",
            "wall_clock_centre_over_ell": tau_c, "wall_clock_half_width_over_ell": delta,
            "primitive": "B(z)=(1-z^2)^4 for abs(z)<1; zero otherwise; z=(tau-tau_c)/delta",
            "input": "I=epsilon B'(z)/delta", "regularity": "B is C3; I is C2",
            "prescribed_characteristics": 1, "outgoing_characteristic_independently_prescribed": False,
            "initial_interior": "phi=0 and d_t phi=0 on flat cap at wall R/ell=.05",
            "outside_support": "incoming trace stipulated zero throughout retained response window",
            "remote_source_or_global_exterior_state_solved": False},
        "support": {"wall_R_over_ell": support_radii.tolist(),
            "wall_tau_over_ell": support["tau_wall_over_ell"].tolist(),
            "wall_flat_t_over_ell": support["t_wall_flat_over_ell"].tolist(),
            "incoming_null_label_over_ell": support_v.tolist(),
            "centre_arrival_t_over_ell": arrival_times.tolist(),
            "centre_response_window_end_over_ell": 0.,
            "initial_wall_R_over_ell": INITIAL_WALL_RADIUS,
            "initial_flat_eta_over_ell": initial_eta},
        "transfer": {"arrival": "T=t_wall-ell log(1-R/ell)",
            "jacobian": "T'=gamma(1+v_peculiar)/(1-R/ell)",
            "response": "phi_C(T(tau))=I(tau)/T'", "centre": "phi_C=2F'(eta_C)/a_C",
            "interior": "chi=F(eta+r)-F(eta-r), with F(v_wall(tau))=epsilon B(z)/2",
            "scope": "finite interior initial-boundary problem with supplied incoming trace"},
        "samples": {key: value.tolist() for key, value in samples.items()},
        "peaks": {"positive_centre": positive, "negative_centre": negative,
            "input_abs_peak_ell_I_over_epsilon": input_peak,
            "input_peak_tau_over_ell": [tau_c-delta/math.sqrt(7), tau_c+delta/math.sqrt(7)]},
        "metrics": {"support_stretch_ratio": float(np.diff(arrival_times)[0])/(2*delta),
            "negative_to_positive_peak_magnitude": abs(negative_peak)/positive_peak,
            "input_squared_amplitude_integral_ell_over_epsilon2": integrals[-1]["input_squared"],
            "output_squared_amplitude_integral_ell_over_epsilon2": integrals[-1]["output_squared"],
            "squared_amplitude_ratio": integrals[-1]["output_squared"]/integrals[-1]["input_squared"],
            "signed_output_area_over_epsilon": integrals[-1]["signed_area"],
            "squared_amplitude_is_energy": False, "late_central_field_through_t0": "exactly zero for the stipulated state"},
        "no_return_certificate": {"initial_null_label_gap_over_ell": initial_gap,
            "earliest_arrival_after_wall_R08_eta_gap_over_ell": after_horizon_gap,
            "wall_radius_lower_bound_after_first_arrival_over_ell": .8,
            "reflected_comoving_radius_upper_bound_through_t0_over_ell": reflected_radius_bound,
            "strict_radius_gap_over_ell": no_return_gap,
            "basis": "first centre arrival after R=.8ell; until t=0, a<=1 and r_wall>=R>=.8ell; reflected r<=-ell-v_min",
            "does_not_exclude_uncomputed_exterior_tails": True},
        "checks": {"clock_32_64_max_difference_over_ell": clock_error,
            "support_tau_endpoint_max_error_over_ell": float(np.max(abs(support["tau_wall_over_ell"]-targets))),
            "integral_32_64_max_difference": max(abs(integrals[0][key]-integrals[1][key]) for key in integrals[0]),
            "input_squared_amplitude_analytic": exact_input_norm,
            "input_squared_amplitude_analytic_error": abs(integrals[-1]["input_squared"]-exact_input_norm),
            "surface_trace_vs_extrinsic_ricci_error_times_ell": ell*abs(trace_ricci-junction["ricci_delta"]),
            "conformal_normal_chi_jump_times_ell": ell*junction["normal_chi_jump_per_chi"],
            "symmetric_radial_centre_limit_relative_error": abs(radial_limit-positive_peak)/abs(positive_peak),
            "radial_limit_probe_r_over_ell": radial_probe,
            "minimum_arrival_jacobian": float(np.min(samples["dT_dtau"])),
            "area_identity": "integral phi_C dT = integral I d_tau = epsilon integral B'(z) dz = 0",
            "radial_wave_identity": "d_eta^2 chi=d_r^2 chi=F''(eta+r)-F''(eta-r)",
            "quadrature_comparison_is_rigorous_error_bound": False, "analytic_zero_checks": zero_checks},
        "limits": ["incoming state is supplied at the wall, not generated by a solved exterior source",
            "conformal cancellation removes the shell delta potential, not Schwarzschild bulk scattering",
            "no quantum vacuum, detector transition probability, energy flux or backreaction calculated",
            "absence of a postpulse centre field is conditional; no generic absence of memory claimed",
            "no late-parent driver, inter-sibling signal, dark component or inherited constant inferred",
            "the isolated reference and its prescribed state are not a solved many-body future"],
    }
