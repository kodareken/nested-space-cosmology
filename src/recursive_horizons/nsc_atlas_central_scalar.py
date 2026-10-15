"""Central scalar coefficient from a supplied exterior wall response.

The cap equation is integrated on an incoming characteristic and reflected
regularly at its centre. All quantities are formal perturbation coefficients;
the module selects no physical amplitude and performs no exterior evolution.
"""
import math

import numpy as np

from .nsc_atlas_cap_response import (
    MAX_QUERY, QUADRATURE_ORDERS, _budget, _choice, _jet, _number, _radial_cuts,
    _time, radial_metric,
)
from .nsc_atlas_cap_wall_response import reference_coefficients


SCHEMA = "nsc-atlas-central-scalar-v1"
SCOPE = "formal central scalar coefficient on the retained cap; no finite-amplitude or uniform error bound"


def optical_time_derivative(cap_field, T, r, *, radial_order=16, budget=None):
    """Evaluate D_T using first field jets and the central proper-clock anchor."""
    order = _choice(radial_order, QUADRATURE_ORDERS, "cap radial order")
    budget = _budget(cap_field, budget)
    if any(np.asarray(value).dtype.kind not in "iuf" for value in (T, r)):
        raise ValueError("optical derivative endpoints must be real numbers, not booleans")
    T, r = np.broadcast_arrays(np.asarray(T, dtype=float), np.asarray(r, dtype=float))
    if not 1 <= T.size <= MAX_QUERY:
        raise ValueError("request at most 513 paired optical derivative endpoints")
    endpoint = _jet(cap_field, T, r, budget)
    shape = T.shape
    times, radii = T.ravel(), r.ravel()
    profile = cap_field.profile(times, second=True)
    pc, pcT = 2*np.asarray(profile["f_prime"]), 2*np.asarray(profile["f_second"])
    if pc.shape != times.shape or pcT.shape != times.shape or not np.all(np.isfinite(pc+pcT)):
        raise ValueError("invalid central profile derivatives")
    anchor = -pc*pcT/6
    points, point_times, point_weights, owners = [], [], [], []
    nodes, weights = np.polynomial.legendre.leggauss(order)
    for index, (time, upper) in enumerate(zip(times, radii)):
        cuts = _radial_cuts(float(time), float(upper), cap_field.edges)
        for lower, last in zip(cuts[:-1], cuts[1:]):
            half = (last-lower)/2
            points.extend(lower+half*(1+nodes))
            point_times.extend(np.full(order, time))
            point_weights.extend(half*weights)
            owners.extend(np.full(order, index))
    radial = np.zeros(T.size)
    if points:
        rs, ts, ws = map(np.asarray, (points, point_times, point_weights))
        jets = _jet(cap_field, ts, rs, budget)
        p, pT, pr = (jets[key] for key in ("p", "p_T", "p_r"))
        A = 1-rs*rs
        density = ((1-3*rs*rs)*pT*pr-2*rs*p*pT)/A
        radial = np.bincount(np.asarray(owners, dtype=int), ws*density, minlength=T.size)
    result = {"D_T": anchor+radial, "central_anchor": anchor, "radial_integral": radial,
              "central_profile_defect": pc-endpoint["p_center"].ravel()}
    if not all(np.all(np.isfinite(value)) for value in result.values()):
        raise ValueError("nonfinite optical derivative")
    budget.check()
    return {**{key: value.reshape(shape) for key, value in result.items()}, "radial_nodes": len(points)}


def incoming_response_integrals(cap_field, q, x, *, radial_order=16, ray_order=16, budget=None):
    """Integrate W1*c0 and D_T/A on the full reflected incoming ray.

    The integral of W1*c0 enters the central coefficient with a minus sign.
    Both reflected-profile and central-anchor edges split the ray. The quiet
    pre-arrival mass shortcut is intentionally inapplicable to these points.
    """
    q, x = _number(q, "central arrival label"), _number(x, "wall radius")
    radial_order = _choice(radial_order, QUADRATURE_ORDERS, "cap radial order")
    ray_order = _choice(ray_order, QUADRATURE_ORDERS, "cap ray order")
    budget = _budget(cap_field, budget)
    if not cap_field.domain["qa"] <= q <= cap_field.domain["qb"] or not .45 <= x <= .55:
        raise ValueError("incoming ray leaves the supplied arrival or wall interval")
    u = q-2*math.atanh(x)
    if u < cap_field.domain["qmin"]-8*np.finfo(float).eps or u >= cap_field.domain["qa"]:
        raise ValueError("incoming ray requires the retained no-return wall labels")
    # Existing field guards verify both endpoint events against the actual cap.
    _jet(cap_field, np.array([q, q-math.atanh(x)]), np.array([0., x]), budget)
    cuts = [0., x]
    for edge in cap_field.edges:
        if q > edge:
            for factor in (1., .5):
                point = math.tanh(factor*(q-float(edge)))
                if 0 < point < x:
                    cuts.append(point)
    cuts = np.unique(cuts)
    nodes, weights = np.polynomial.legendre.leggauss(ray_order)
    points, point_weights = [], []
    for lower, last in zip(cuts[:-1], cuts[1:]):
        half = (last-lower)/2
        points.extend(lower+half*(1+nodes))
        point_weights.extend(half*weights)
    radii, weights = np.asarray(points), np.asarray(point_weights)
    integrands = {name: np.empty(len(radii)) for name in ("forcing", "optical_time", "lambda_plus")}
    reflected_defect = central_defect = 0.
    optical_nodes = 0
    center_samples = []
    for first in range(0, len(radii), MAX_QUERY):
        rs = radii[first:first+MAX_QUERY]
        ts = q-np.arctanh(rs)
        metric = radial_metric(cap_field, ts, rs, radial_order=radial_order, budget=budget)
        optical = optical_time_derivative(cap_field, ts, rs, radial_order=radial_order, budget=budget)
        p, pT, pr, mass, lapse = (metric[key] for key in ("p", "p_T", "p_r", "nu_E", "delta_E"))
        A = 1-rs*rs
        W1 = 2*mass/rs**3-p*p/3-pT*pT/(6*A)+A*pr*pr/6
        c0 = rs*p  # Stable regular field, including its reflected term.
        for index in range(min(2-len(center_samples), len(rs))):
            K1 = float(A[index]*W1[index]/4)
            center_samples.append({"r": float(rs[index]), "T": float(ts[index]),
                "W1": float(W1[index]), "K1": K1, "K1_over_r_squared": K1/float(rs[index]**2)})
        D = lapse-2*mass/(rs*A)
        integrands["forcing"][first:first+len(rs)] = W1*c0
        integrands["optical_time"][first:first+len(rs)] = optical["D_T"]/A
        integrands["lambda_plus"][first:first+len(rs)] = -D/A
        # Public profile batches remain bounded even for a full metric batch.
        incoming_f = float(cap_field.profile(np.array([q]))["f"][0])
        outgoing_f = cap_field.profile(q-2*np.arctanh(rs))["f"]
        reflected_defect = max(reflected_defect, float(np.max(abs(c0-incoming_f+outgoing_f))))
        central_defect = max(central_defect, float(np.max(abs(optical["central_profile_defect"]))))
        optical_nodes += optical["radial_nodes"]
    result = {}
    for name, values in integrands.items():
        if not np.all(np.isfinite(values)):
            raise ValueError("nonfinite incoming scalar integral")
        positive = float(weights @ np.maximum(values, 0.))
        negative = float(weights @ np.minimum(values, 0.))
        result.update({name+"_integral": positive+negative, name+"_positive": positive,
                       name+"_negative": negative, name+"_sample_min": float(min(values)),
                       name+"_sample_max": float(max(values))})
    budget.check()
    return {"q": q, "x": x, "u_wall": u, **result, "ray_nodes": len(radii),
            "optical_radial_nodes": optical_nodes, "reflected_profile_defect": reflected_defect,
            "central_profile_defect": central_defect, "analytic_regular_center_K1": 0.,
            "center_samples": center_samples,
            "first_arrival_margin": q-cap_field.domain["qa"],
            "last_arrival_margin": cap_field.domain["qb"]-q,
            "no_return_margin": cap_field.domain["qa"]-u}


def reconstruct_central_scalar(wall, cap_row, profile, integrals, *, y_wall,
                               y_wall_s, J_minus, arrival_spline_derivative):
    """Algebraic reconstruction from supplied compatible cap/wall rows.

    y_wall_s is the material derivative of the exterior areal coefficient,
    i.e. V0_s times the exterior solver's total wall derivative with respect
    to V. Cap T1/e/ray data may be reused only after their preparation and
    numerical settings have been independently authenticated by the caller.
    C_fixedV is recomputed from the supplied current exterior V1.
    """
    y, ys, J, qV = [_number(value, label) for value, label in (
        (y_wall, "wall scalar coefficient"), (y_wall_s, "material wall scalar derivative"),
        (J_minus, "outgoing scalar derivative"), (arrival_spline_derivative, "arrival spline derivative"))]
    names = ("x", "v", "beta_minus", "beta_plus", "V", "Vdot", "V1", "V1_s")
    x, v, bm, bp, V, Vs, V1, V1s = [_number(wall[key], key) for key in names]
    if not .45 <= x <= .55 or Vs <= 0 or V >= 0 or v <= bp or qV <= 0:
        raise ValueError("central reconstruction leaves the signed reference branch")
    q, e, es, T1, T1s, mass, lapse = [_number(cap_row[key], key) for key in
        ("q0", "e", "e_s", "T1", "T1_s", "nu_E_wall", "delta_E_wall")]
    f, fp, fpp = [_number(profile[key], key) for key in ("f", "f_prime", "f_second")]
    I, IT, lam = [_number(integrals[key], key) for key in
        ("forcing_integral", "optical_time_integral", "lambda_plus_integral")]
    for key, expected in (("q", q), ("x", x)):
        if not math.isclose(_number(integrals[key], key), expected, rel_tol=0., abs_tol=1e-13):
            raise ValueError("incoming integral and wall coordinates differ")
    A, As = 1-x*x, -2*x*v
    Q, Qminus = (bm+v)/A, (bm-v)/A
    if Q <= 0 or Qminus <= 0:
        raise ValueError("central reconstruction requires monotone cap null labels")
    for key, expected in (("x", x), ("q0_s", Q), ("V0_s", Vs)):
        if key in cap_row and not math.isclose(_number(cap_row[key], key), expected,
                                               rel_tol=2e-11, abs_tol=2e-12):
            raise ValueError("saved cap and current wall reference data differ: "+key)
    reference = reference_coefficients(x)
    acceleration, bm_x, bp_x = [float(reference[key]) for key in
        ("radial_acceleration", "beta_minus_x", "beta_plus_x")]
    Qs = (bm_x*v+acceleration)/A-Q*As/A
    d, ds = v-bp, acceleration-bp_x*v
    Vss = Vs*(Vs/V-ds/d)
    ratio, ratio_s = Q/Vs, Qs/Vs-Q*Vss/Vs**2
    Dwall = lapse-2*mass/(x*A)
    Delta = T1+e/A+lam
    Delta_s = T1s+es/A-e*As/A**2-v*Dwall/A-Q*IT
    C, Cs = Delta-ratio*V1, Delta_s-ratio_s*V1-ratio*V1s
    H, Hq = y-C*fp, ys/Q-Cs*fp/Q-C*fpp
    wall_incoming = Hq-Qminus*J/Q
    p0, p1 = 2*fp, 2*wall_incoming-I
    result = dict(q=q, x=x, y_wall=y, y_wall_s=ys, J_minus=J, f=f,
        f_prime=fp, f_second=fpp, q_plus_s=Q, q_minus_s=Qminus,
        q_plus_ss=Qs, V0_s=Vs, V0_ss=Vss, D_wall=Dwall,
        Delta_C_fixed_s=Delta, Delta_C_fixed_s_s=Delta_s,
        Delta_C_fixed_V=C, Delta_C_fixed_V_s=Cs, H=H, H_q=Hq,
        incoming_wall_derivative=wall_incoming, p0_center=p0, p1_center=p1,
        exterior_trace_term=2*ys/Q, arrival_rate_term=-2*Cs*fp/Q,
        arrival_shift_term=-2*C*fpp, outgoing_wall_term=-2*Qminus*J/Q,
        incoming_forcing_term=-I, arrival_continuum_derivative=ratio,
        arrival_spline_derivative=qV, arrival_rate_defect=qV*Vs-Q,
        static_center_derivative=p1+p0**3/12,
        jordan_frame_subtraction=-p0**3/12,
        central_frame_algebra_defect=(p1+p0**3/12)-p0**3/12-p1)
    if "ray_integral" in cap_row:
        result["reused_ray_integral_defect"] = lam-_number(cap_row["ray_integral"], "saved ray integral")
    if not all(math.isfinite(value) for value in result.values()):
        raise ValueError("nonfinite central scalar coefficient")
    return result


class CentralScalarResponse:
    """Live-cap adapter; lower-level reconstruction also accepts saved rows."""

    def __init__(self, cap_solution):
        self.cap_solution = cap_solution
        self.cap_field, self.wall_solution = cap_solution.cap_field, cap_solution.wall_solution
        self.budget = _budget(self.cap_field, cap_solution.budget)
        if self.wall_solution.budget is not self.budget:
            raise ValueError("central response requires the same cumulative wall budget")
        if self.cap_field.prepared_field is not self.wall_solution.wedge.prepared_field:
            raise ValueError("central response requires the same prepared leading field")
        self.domain, self.settings = dict(cap_solution.domain), dict(cap_solution.settings)
        self.cache = {}

    def sample(self, t, *, y_wall, y_wall_s, J_minus):
        t = _time(t)
        if t not in self.cache:
            if len(self.cache) >= 129:
                raise ValueError("retain at most 129 central scalar geometric samples")
            cap = self.cap_solution.sample(t)
            wall = self.wall_solution.sample(t)
            q = cap["q0"]
            profile = {key: float(value[0]) for key, value in
                       self.cap_field.profile(np.array([q]), second=True).items()}
            integrals = incoming_response_integrals(self.cap_field, q, cap["x"],
                radial_order=self.settings["radial_order"], ray_order=self.settings["ray_order"], budget=self.budget)
            qV = float(self.cap_field.geometry["arrival"](wall["V"], 1))
            self.cache[t] = (wall, cap, profile, integrals, qV)
        wall, cap, profile, integrals, qV = self.cache[t]
        result = reconstruct_central_scalar(wall, cap, profile, integrals,
            y_wall=y_wall, y_wall_s=y_wall_s, J_minus=J_minus, arrival_spline_derivative=qV)
        self.budget.check()
        return {"t": t, **result, "integrals": dict(integrals),
                "endpoint_is_one_sided_formal_limit": t in (0., 1.),
                "uniform_error_enclosure": False}
