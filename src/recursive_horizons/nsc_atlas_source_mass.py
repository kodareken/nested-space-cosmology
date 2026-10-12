"""Quadratic geometric mass on the existing conformal source characteristic.

The mixed radius equation and the already determined transverse scalar close
this initial-leg question. Neither the field nor the coupled bulk is evolved.
"""
import math

import numpy as np

from .nsc_atlas_characteristic_initialization import (
    Budget, amplitude_parameter, constraint_samples, source_geometry, source_jet,
    validate_geometry, validate_model,
)
from .nsc_atlas_field_preparation import exterior_coefficient, radius_from_null, source_profile

SCHEMA = "nsc-atlas-source-mass-v1"
MAX_SAMPLES = 129
MAX_EVALUATIONS = 100000
CPU_SECONDS = 60
TOLERANCE_ABSOLUTE = 1e-8
TOLERANCE_RELATIVE = 1e-8


def validate_mass_inputs(*, corner_mass_correction=0., corner_transverse_correction=0., **model):
    validate_model(**model)
    if not all(math.isfinite(float(v)) for v in (corner_mass_correction, corner_transverse_correction)):
        raise ValueError("corner mass and transverse data must be finite")
    if corner_mass_correction or corner_transverse_correction:
        raise ValueError("the original corner mass and transverse radius derivative are fixed")


def null_coordinate_boost(geometry, factor):
    """Unew=bU, Vnew=V/b; the geometric cut and dimensionless source stay fixed."""
    factor = float(factor)
    if not math.isfinite(factor) or factor <= 0:
        raise ValueError("null-coordinate boost must be finite and positive")
    validate_geometry(geometry)
    transformed = dict(geometry)
    for key in ("U0", "Ua", "Ub"):
        transformed[key] *= factor
    for key in ("Va", "Vb"):
        transformed[key] /= factor
    transformed["support_V"] = [v/factor for v in geometry["support_V"]]
    validate_geometry(transformed)
    return transformed


def mass_samples(geometry, *, samples=MAX_SAMPLES, order=64, amplitude=1., budget=None):
    """Independent scalar kernels plus frame/mixed-geometry mass expressions."""
    budget = Budget(maximum_evaluations=MAX_EVALUATIONS) if budget is None else budget
    start_count = budget.evaluations
    initialized = constraint_samples(geometry, samples=samples, order=order,
                                     amplitude=amplitude, budget=budget)
    constraint_count = budget.evaluations-start_count
    x, V = initialized["x"], initialized["V"]
    p, dp, ddp = (initialized[key] for key in ("p", "p_prime", "p_second"))
    h, hp, I0 = (initialized[key] for key in ("h", "h_prime", "integral_p_squared"))
    d = (.8/x-1)/4
    C = -3.2/geometry["U0"]
    p_over_t2, mass_integral, kernel_integral = (np.zeros(samples) for _ in range(3))
    nodes, weights = np.polynomial.legendre.leggauss(order)
    support_x0, support_x1 = geometry["support_x"]
    support_V0, support_V1 = geometry["support_V"]
    x_count = v_count = 0
    for i, endpoint in enumerate(x):
        upper = min(endpoint, support_x1)
        if upper <= support_x0:
            continue
        t = (support_x0+upper)/2+(upper-support_x0)*nodes/2
        w = (upper-support_x0)*weights/2
        budget.check(len(t))
        x_count += len(t)
        q, dq, _ = source_jet(t, geometry, amplitude)
        p_over_t2[i] = float(w @ (q/t**2))
        mass_integral[i] = float(w @ (t*(.8-t)*dq*dq))
        # Independent original V-kernel: no substitution of the x primitive.
        upper_V = min(V[i], support_V1)
        v = (support_V0+upper_V)/2+(upper_V-support_V0)*nodes/2
        wv = (upper_V-support_V0)*weights/2
        budget.check(len(v))
        v_count += len(v)
        radius = radius_from_null(geometry["U0"], v)
        b = amplitude*source_profile(v, geometry["Va"], geometry["Vb"])[0]
        kernel_integral[i] = float(wv @ (exterior_coefficient(radius)*b))
    s = -d*p/x-.8*p_over_t2/(4*x)
    s_kernel = (-kernel_integral/C-d*p)/x
    s_prime = -(s+d*dp)/x
    mu_E = -mass_integral/4
    frame_term = .4*p*p/12+x*x*p*(s+d*dp)/3
    mu_J = mu_E+frame_term
    k = -d*h/x+I0/(24*x)+d*p*p/6+x*p*s/6
    mu_J_mixed = .4*h/x+2*x*(k+d*hp)
    T_xx = (2*dp*dp-p*ddp)/3
    T_Ux_over_C = (dp*s-p*s_prime)/3
    mu_J_prime_stress = x*x*(T_Ux_over_C-d*T_xx)
    mu_J_prime_frame = (-x*x*d*dp*dp+.4*p*dp/6
                        +(2*x*p*(s+d*dp)+x*x*dp*(s+d*dp)
                          +x*x*p*(s_prime-.8*dp/(4*x*x)+d*ddp))/3)
    budget.check()
    data = {"x": x, "V": V, "p": p, "p_prime": dp, "p_second": ddp,
            "h": h, "h_prime": hp, "d": d, "s": s, "s_kernel": s_kernel,
            "s_prime": s_prime, "k": k, "mu_einstein": mu_E, "mu_jordan": mu_J,
            "mu_geometric": mu_J_mixed, "frame_mass_difference": frame_term,
            "mu_jordan_prime_stress": mu_J_prime_stress, "mu_jordan_prime_frame": mu_J_prime_frame,
            "T_xx_normalized": T_xx, "T_Ux_over_C_normalized": T_Ux_over_C,
            "integral_p_over_x_squared": p_over_t2, "integral_p_squared": I0,
            "integral_x_a_minus_x_dp_squared": mass_integral,
            "integral_V_kernel_b": kernel_integral}
    return data, {"constraint": constraint_count, "mass_and_transverse_x": x_count,
                  "original_V_kernel": v_count}


def physical_mass_length_coefficient(mu_J, epsilon, kappa_4, length_scale):
    """Formal quadratic mass-length shift, not a finite-amplitude validity bound."""
    if not math.isfinite(float(mu_J)):
        raise ValueError("finite mass coefficient required")
    return float(length_scale)*amplitude_parameter(epsilon, kappa_4, length_scale)*mu_J


def exact_algebra_checks():
    import sympy as sp
    if sp.__version__ != "1.14.0":
        raise RuntimeError("exact mass checks require pinned SymPy 1.14.0")
    x, a, U, C, boost = sp.symbols("x a U C boost", positive=True)
    p, s, h = (sp.Function(name)(x) for name in ("p", "s", "h"))
    I, Q, S, W = (sp.Function(name)(x) for name in ("I", "Q", "S", "W"))
    d = (a/x-1)/4
    sp_rule = -(s+d*sp.diff(p, x))/x
    s_formula = -d*p/x-a*Q/(4*x)
    muEprime = -x*x*d*sp.diff(p, x)**2
    frame = a*p*p/24+x*x*p*(s+d*sp.diff(p, x))/3
    stress = x*x*((sp.diff(p, x)*s-p*sp_rule)/3
                  -d*(2*sp.diff(p, x)**2-p*sp.diff(p, x, 2))/3)
    H = h-x*p*p/12
    k = -d*h/x+I/(24*x)+d*p*p/6+x*p*s/6
    kE = k-d*p*p/12-x*p*s/6
    h_moment = -W/2+x*p*p/12-I/6
    hp_moment = -S/2+x*p*sp.diff(p, x)/6-p*p/12
    mu_mixed = a*h/(2*x)+2*x*(k+d*sp.diff(h, x))
    muE_moment = -((a-x)*S+W)/4
    K = a**4*sp.exp(-x/a)/x**4
    Fprime = -x*sp.exp(x/a)/a**2
    Cexpr = -4*a/U
    residuals = {
        "background_transverse_constraint": sp.diff(x*d, x)+sp.Rational(1, 4),
        "transverse_scalar_transport": (sp.diff(s_formula, x)+(s_formula+d*sp.diff(p, x))/x)
                                       .subs(sp.diff(Q, x), p/x**2),
        "original_kernel_change_of_variable": K*(x*p)*Fprime/U/Cexpr-a*p/(4*x**2),
        "mixed_radius_integral": x*kE+d*H-I/24,
        "mixed_radius_equation": sp.diff(x*kE+d*H, x).subs(sp.diff(I, x), p*p)-p*p/24,
        "jordan_mixed_improved_stress": sp.diff(x*k+d*h, x).subs(
            {sp.diff(I, x): p*p, sp.diff(s, x): sp_rule})
            -x*x*(sp.diff(p, x)*s-p*sp_rule)/6,
        "mixed_and_frame_mass": (mu_mixed-muE_moment-frame).subs(
            {sp.diff(h, x): hp_moment, h: h_moment}),
        "jordan_improved_stress_mass_transport": muEprime+sp.diff(frame, x).subs(sp.diff(s, x), sp_rule)-stress,
        "post_support_mass_from_radius": (h_moment+(a-x)*hp_moment)/2+I/12-muE_moment,
        "null_boost_kernel_normalization": (Q/boost)/(C/boost)-Q/C,
    }
    residuals["post_support_mass_from_radius"] = residuals["post_support_mass_from_radius"].subs(
        {sp.diff(p, x): 0, p: 0})
    zeros = {key: str(sp.simplify(value)) for key, value in residuals.items()}
    if any(value != "0" for value in zeros.values()):
        raise ValueError("exact source mass identity failed: "+str(zeros))
    controls = {
        "omitting_transverse_frame_term_changes_mass": sp.simplify(x*x*p*s/3) != 0,
        "dropping_improved_mixed_stress_changes_transport": sp.simplify(x*x*(sp.diff(p, x)*s-p*sp_rule)/3) != 0,
    }
    if not all(controls.values()):
        raise ValueError("source mass negative control unexpectedly vanished")
    return {"zero_residuals": zeros, "nonzero_controls": controls}


def mass_calculation(budget=None):
    budget = Budget(maximum_evaluations=MAX_EVALUATIONS) if budget is None else budget
    validate_mass_inputs()
    g32, g64 = source_geometry(32, budget), source_geometry(64, budget)
    fine, counts64 = mass_samples(g64, order=64, budget=budget)
    coarse_quad, counts32 = mass_samples(g64, order=32, budget=budget)
    coarse_geo, countsgeo = mass_samples(g32, order=64, budget=budget)
    refinements = {}
    for category, coarse in (("source_quadrature", coarse_quad), ("geometry_quadrature", coarse_geo)):
        refinements[category] = {}
        for key in fine:
            delta = abs(fine[key]-coarse[key])
            refinements[category][key] = float(np.max(delta))
            if np.any(delta > TOLERANCE_ABSOLUTE+TOLERANCE_RELATIVE*abs(fine[key])):
                raise RuntimeError("source mass refinement failed: "+category+"/"+key)
    comparisons = {}
    for left, right in (("s", "s_kernel"), ("mu_jordan", "mu_geometric"),
                        ("mu_jordan_prime_frame", "mu_jordan_prime_stress")):
        delta = abs(fine[left]-fine[right])
        comparisons[left+"_vs_"+right] = float(np.max(delta))
        if np.any(delta > TOLERANCE_ABSOLUTE+TOLERANCE_RELATIVE*abs(fine[left])):
            raise RuntimeError("independent source mass comparison failed: "+left)
    if not fine["mu_jordan"][-1] < 0 or fine["mu_jordan"][-1] != fine["mu_einstein"][-1]:
        raise RuntimeError("analytic post-source mass sign/frame agreement failed")
    count_sums = {key: sum(row[key] for row in (counts64, counts32, countsgeo)) for key in counts64}
    if budget.evaluations != sum(count_sums.values())+1440:
        raise RuntimeError("source mass integrand accounting mismatch")
    budget.check()
    return {"geometry": {**g64, "C": -3.2/g64["U0"]},
            "samples": {key: value.tolist() for key, value in fine.items()},
            "last": {key: float(value[-1]) for key, value in fine.items()},
            "checks": {"quadrature_32_64_max_differences": refinements,
                       "independent_comparison_max_differences": comparisons,
                       "during_source_max_abs_frame_difference": float(np.max(abs(fine["frame_mass_difference"]))),
                       "post_source_mass_max_minus_min": float(np.ptp(fine["mu_jordan"][96:])),
                       "tolerance_absolute": TOLERANCE_ABSOLUTE, "tolerance_relative": TOLERANCE_RELATIVE,
                       "refinement_is_rigorous_enclosure": False},
            "execution": {"displayed_source_samples": MAX_SAMPLES, "quadrature_orders": [32, 64],
                          "geometry_integrand_evaluations": 1440,
                          "reused_constraint_integrand_evaluations": count_sums["constraint"],
                          "mass_and_transverse_x_integrand_evaluations": count_sums["mass_and_transverse_x"],
                          "original_V_kernel_integrand_evaluations": count_sums["original_V_kernel"],
                          "total_integrand_evaluations": budget.evaluations,
                          "maximum_integrand_evaluations": MAX_EVALUATIONS, "cpu_seconds_max": CPU_SECONDS,
                          "field_solves": 0, "bulk_einstein_evolutions": 0}}
