"""Formal cap mass and moving-wall coefficients for the fixed atlas reference.

These functions evaluate algebraic ingredients only. They neither integrate a
mass functional nor solve a bulk spacetime. The cap coefficient is relative to
its vacuum mass at the actual displaced radius. The owning admission limits
the physical application to the original conformal scalar and no-return strip.
"""
import time

import numpy as np

from .nsc_atlas_causal_reference import signed_state


SCHEMA = "nsc-atlas-cap-wall-response-v1"
SYMPY_VERSION = "1.14.0"
CPU_SECONDS = 60
H_ELL = 1.0
K_ELL = 2.0
M_OVER_ELL = 0.4
WALL_LOWER = 0.45
WALL_UPPER = 0.55

FORMULAS = {
    "kinetic_density": "x**2*(beta*A*B-v*(A**2+B**2)/2)/2",
    "cap_potential_density": "h**2*x**2*v*p**2/2",
    "cap_endpoint": "x**2*p*(v*A-beta*B)/6-h**2*x**3*p**2/12",
    "cap_endpoint_trace": "(v-beta)*w*j/6+(1/(6*x)-h**2*x/4)*w**2",
    "exterior_endpoint": "m*p**2/12+x**2*p*(v*A-beta*B)/6",
    "improved_t_un": "2*A*B/3-p*Hun/3",
    "cap_improved_t_nn": "A**2/6+B**2/2-p*Hnn/3+h**2*p**2/6",
    "cap_wave_hessian": "Huu+2*(v*A-beta*B)/x+2*h**2*p",
    "cap_total_mass": "nu+3*h**2*x**2*rho/2",
}

IDENTITIES = (
    "cap_full_improved_normal_projection",
    "transformed_fixed_vacuum_potential",
    "cap_frame_endpoint_subtracts_actual_radius_vacuum",
    "cap_transport_equals_direct_jordan_stress",
    "cap_endpoint_acceleration_cancels",
    "cap_endpoint_trace_form",
    "cap_actual_radius_background_displacement",
    "signed_background_angular_jump",
    "signed_background_mass_difference",
    "background_acceleration_agrees_on_both_sides",
    "linear_angular_jump",
    "linear_exterior_mass_constraint",
    "linear_cap_mass_constraint",
    "two_displacement_formulas_agree",
    "mass_difference_transport_from_equal_scalar_stresses",
    "mass_difference_from_both_first_derivative_functionals",
    "temporal_jump_from_angular_and_mixed_flux",
    "normal_force_from_two_radial_einstein_equations",
    "linear_acceleration_from_transported_mass_difference",
    "complete_exterior_acceleration_numerator_is_beta_times_delta_a",
    "linear_mean_curvature_variation_vanishes",
    "conformal_material_trace_retains_tension_jump",
    "zero_source_mass_ingredients",
    "quadratic_field_scaling",
    "normal_reversal_preserves_cap_mass_ingredients",
    "physical_mass_and_displacement_unit_scaling",
)
CONTROLS = (
    "omitted_transformed_potential_changes_transport",
    "omitted_cap_einstein_tensor_changes_normal_stress",
    "omitted_cap_displacement_changes_total_mass",
    "omitted_endpoint_changes_jordan_transport",
    "wrong_mixed_flux_sign_changes_temporal_jump",
    "changed_cap_mass_without_exterior_changes_wall_forcing",
)


def _finite_arrays(**values):
    arrays = []
    for name, value in values.items():
        raw = np.asarray(value)
        if raw.dtype.kind not in "iuf" or raw.size == 0:
            raise ValueError(name + " must contain finite real numbers, not booleans")
        array = np.asarray(raw, dtype=float)
        if not np.all(np.isfinite(array)):
            raise ValueError(name + " must contain finite real numbers")
        arrays.append(array)
    try:
        return np.broadcast_arrays(*arrays)
    except ValueError as error:
        raise ValueError("input shapes are not broadcast compatible") from error


def reference_coefficients(x):
    """Signed reference data on the admitted expanding strip, including beta=0."""
    (x,) = _finite_arrays(x=x)
    if np.any((x < WALL_LOWER) | (x > WALL_UPPER)):
        raise ValueError("radius exits the admitted 0.45 <= x <= 0.55 wall strip")
    state = signed_state(x)
    v = state["Rdot"]
    if not np.all(np.isfinite(v)) or np.any(v <= 0):
        raise ValueError("the admitted wall branch must have positive proper radial speed")
    m, h, k = M_OVER_ELL, H_ELL, K_ELL
    bp_x = -2*m/(k*x**3)-(k*k+h*h)/(2*k)
    bm_x = -2*m/(k*x**3)+(k*k-h*h)/(2*k)
    return {
        "x": x, "v": v,
        "beta_minus": state["beta_minus"], "beta_plus": state["beta_plus"],
        "beta_minus_x": bm_x, "beta_plus_x": bp_x,
        "beta_second": 6*m/(k*x**4),
        "radial_acceleration": state["beta_plus"]*bp_x-m/x**2,
        "cap_background_mass": h*h*x**3/2,
    }


def field_mass_ingredients(x, w, j):
    """Return first-derivative integrands/endpoints; no mass integral is taken.

    w=chi_W/epsilon and j=dw/ds are continuum leading traces. With quiet
    anchors, mu_plus=integral(L_plus ds)+Q_plus and
    nu_minus=integral(L_minus ds)+Q_minus. nu_minus excludes h^2 r^3/2.
    """
    x, w, j = _finite_arrays(x=x, w=w, j=j)
    geometry = reference_coefficients(x)
    v, bm, bp = (geometry[key] for key in ("v", "beta_minus", "beta_plus"))
    p = w/x
    A, Bm, Bp = (j-v*p)/x, (j-bm*p)/x, (j-bp*p)/x

    def kinetic(beta, B):
        return x*x*(beta*A*B-v*(A*A+B*B)/2)/2

    potential = H_ELL**2*x*x*v*p*p/2
    Qm = x*x*p*(v*A-bm*Bm)/6-H_ELL**2*x**3*p*p/12
    Qp = M_OVER_ELL*p*p/12+x*x*p*(v*A-bp*Bp)/6
    return {
        "p": p, "A": A, "B_minus": Bm, "B_plus": Bp,
        "L_minus": kinetic(bm, Bm)+potential, "Q_minus": Qm,
        "L_plus": kinetic(bp, Bp), "Q_plus": Qp,
        "cap_potential_density": potential,
        "Q_minus_trace": (v-bm)*w*j/6+(1/(6*x)-H_ELL**2*x/4)*w*w,
        "Q_plus_trace": (v-bp)*w*j/6+(1/(6*x)-M_OVER_ELL/(4*x*x))*w*w,
    }


def cap_total_mass_coefficient(x, displacement, excess_mass):
    """Coefficient of cap geometric mass at its actual displaced wall radius."""
    x, displacement, excess_mass = _finite_arrays(
        x=x, displacement=displacement, excess_mass=excess_mass)
    reference_coefficients(x)
    return excess_mass+1.5*H_ELL**2*x*x*displacement


def wall_response_coefficients(x, displacement, exterior_mass, cap_excess_mass):
    """Formal proper-time wall coefficients from compatible one-sided masses.

    The two masses must be transported from the fixed quiet anchors by their
    own mass laws. Algebraic inputs alone do not certify that history. The
    normal-acceleration expression additionally uses equal scalar mixed and
    normal stresses, so D_s=-k*x^3*t_uN/2. No trajectory is integrated here.
    """
    x, rho, mu, nu = _finite_arrays(x=x, displacement=displacement,
                                   exterior_mass=exterior_mass,
                                   cap_excess_mass=cap_excess_mass)
    g = reference_coefficients(x)
    v, bp, bm = (g[key] for key in ("v", "beta_plus", "beta_minus"))
    difference = mu-nu
    beta_input = difference/(K_ELL*x*x)
    forcing = bp*beta_input+mu/x
    rho_s = (g["radial_acceleration"]*rho+forcing)/v
    return {
        "mass_difference": difference,
        "displacement_rate": rho_s,
        "displacement_forcing": forcing/v,
        "delta_beta_plus": beta_input+g["beta_plus_x"]*rho,
        "delta_beta_minus": beta_input+g["beta_minus_x"]*rho,
        "cap_total_mass": nu+1.5*H_ELL**2*x*x*rho,
        "cap_side_displacement_rate": (g["radial_acceleration"]*rho+bm*beta_input+nu/x)/v,
        "delta_normal_acceleration": -2*difference/(K_ELL*x**3)+g["beta_second"]*rho,
    }


def expected_certificate():
    return {
        "sympy_version": SYMPY_VERSION,
        "zero_residuals": dict.fromkeys(IDENTITIES, "0"),
        "nonzero_controls": dict.fromkeys(CONTROLS, True),
        "arithmetic": "exact symbolic expressions; no mass quadrature, trajectory or bulk evolution",
    }


def exact_algebra_checks(started=None):
    """Check independent stress/frame/junction identities within the CPU budget."""
    started = time.process_time() if started is None else started

    def budget():
        if time.process_time()-started > CPU_SECONDS:
            raise RuntimeError("cap-wall exact checks exceeded process CPU budget")

    budget()
    import sympy as sp
    if sp.__version__ != SYMPY_VERSION:
        raise RuntimeError("cap-wall exact checks require pinned SymPy " + SYMPY_VERSION)
    names = "x h b m p A B v beta j w acc Huu Hun Hnn rho mu nu alpha scale q P rs ass am ap bm bp Tm Tp phi kappa sigma ell epsilon"
    symbols = dict(zip(names.split(), sp.symbols(names)))
    (x, h, b, m, p, A, B, v, beta, j, w, acc, Huu, Hun, Hnn,
     rho, mu, nu, alpha, scale, q, P, rs, ass, am, ap, bm, bp,
     Tm, Tp, phi, kappa, sigma, ell, epsilon) = symbols.values()
    e = {name: sp.sympify(value, locals=symbols) for name, value in FORMULAS.items()}
    L, Q = e["kinetic_density"]+e["cap_potential_density"], e["cap_endpoint"]
    # Direct Hilbert projection retains G_NN=-3h^2 and Box p=2h^2 p.
    box_square = 2*(-A*A+B*B)+4*h*h*p*p
    direct_nn = (A*A+B*B)/2+(-3*h*h*p*p+box_square-2*B*B-2*p*Hnn)/6
    rates = {x: v, p: A, A: Huu+acc*B, B: Hun+acc*A,
             v: acc*beta+h*h*x, beta: acc*v}
    Qs = sum(sp.diff(Q, variable)*rate for variable, rate in rates.items())
    jordan = x*x*(beta*e["improved_t_un"]-v*e["cap_improved_t_nn"])/2
    jordan = jordan.subs(Hnn, e["cap_wave_hessian"])
    omega = 1-alpha*p*p/12
    omega_t, omega_n = -alpha*p*A/6, -alpha*p*B/6
    radius = x+alpha*rho
    original_mass = h*h*radius**3/2+alpha*nu
    frame_mass = (omega*original_mass-radius*radius*(-omega_t*v+omega_n*beta)
                  -radius**3*(-omega_t**2+omega_n**2)/(2*omega))
    frame_excess = frame_mass-h*h*(omega*radius)**3/2
    bm0 = m/(b*x*x)+(b*b-h*h)*x/(2*b)
    bp0 = m/(b*x*x)-(b*b+h*h)*x/(2*b)
    bmx, bpx = sp.diff(bm0, x), sp.diff(bp0, x)
    acceleration = bp0*bpx-m/x**2
    D = mu-nu
    dbp, dbm = D/(b*x*x)+bpx*rho, D/(b*x*x)+bmx*rho
    rho_s = (acceleration*rho+bp0*D/(b*x*x)+mu/x)/v
    cap_rho_s = (acceleration*rho+bm0*D/(b*x*x)+nu/x)/v
    # Both scalar mass laws are written after the cap vacuum subtraction.
    mu_s = x*x*(bp0*q-v*P)/2
    nu_s = x*x*(bm0*q-v*P)/2
    D_s = -b*x**3*q/2
    da = (D_s/(b*x*x)-2*v*D/(b*x**3)+sp.diff(bpx, x)*v*rho+x*q/2)/v
    da_closed = -2*D/(b*x**3)+6*m*rho/(b*x**4)
    # Differentiate the displacement equation before imposing v^2=q0.
    # This checks the whole radial-Einstein numerator, avoiding division by
    # beta_plus even at its regular zero crossing.
    rho_ss = (sp.diff(rho_s, x)*v+sp.diff(rho_s, v)*acceleration
              +sp.diff(rho_s, rho)*rho_s+sp.diff(rho_s, mu)*mu_s
              +sp.diff(rho_s, nu)*nu_s)
    acceleration_numerator = rho_ss+mu/x**2-2*m*rho/x**3+x*P/2-bpx*dbp
    q0 = bm0*bm0-1+h*h*x*x

    def on_background(expression):
        numerator, _ = sp.fraction(sp.cancel(expression))
        return sp.rem(numerator, v*v-q0, v)

    # A stronger mass-difference check uses the two actual first-derivative
    # functionals, not merely a subtraction of their stated transport laws.
    Bp = B+b*p
    Qp = e["exterior_endpoint"].subs({beta: bp0, B: Bp}, simultaneous=True)
    Qm = Q.subs(beta, bm0)
    Lp = e["kinetic_density"].subs({beta: bp0, B: Bp}, simultaneous=True)
    Lm = L.subs(beta, bm0)
    trace_rates = {x: v, p: A, A: Huu+bmx*B,
                   B: Hun+bmx*A, v: acceleration}
    difference_endpoint_rate = sum(sp.diff(Qp-Qm, variable)*rate
                                   for variable, rate in trace_rates.items())
    direct_difference_rate = Lp-Lm+difference_endpoint_rate
    mean_beta = (bm0+bp0)/2
    mean_delta_beta = (dbm+dbp)/2
    # Exact temporal identity from d(beta_minus-beta_plus)=b*r_s.
    temporal = ((am*rs-x*q/2)-(ap*rs-x*q/2)-b*rs)
    # Difference of two radial Einstein equations, after angular+temporal law.
    radial_difference = am*bm-ap*bp+(b*x*x*(bm+bp)/2)/x**2+x*(Tp-Tm)/2
    normal_identity = b*x*((am+ap)/2+(bm+bp)/x+(Tp-Tm)/(2*b))
    f = 1/(2*kappa)-phi*phi/12
    fp = sp.diff(f, phi)
    total_mass_coefficient = sp.diff(h*h*radius**3/2+alpha*nu, alpha).subs(alpha, 0)
    residuals = {
        "cap_full_improved_normal_projection": direct_nn-e["cap_improved_t_nn"],
        "transformed_fixed_vacuum_potential": sp.diff(3*h*h/(1-alpha*p*p/6)**2, alpha).subs(alpha, 0)-h*h*p*p,
        "cap_frame_endpoint_subtracts_actual_radius_vacuum": sp.diff(frame_excess, alpha).subs(alpha, 0)-(nu-Q),
        "cap_transport_equals_direct_jordan_stress": L+Qs-jordan,
        "cap_endpoint_acceleration_cancels": sp.diff(Qs, acc),
        "cap_endpoint_trace_form": (Q.subs({A: (j-v*p)/x, B: (j-beta*p)/x}).subs(p, w/x)-e["cap_endpoint_trace"]).expand().subs(v*v, beta*beta-1+h*h*x*x),
        "cap_actual_radius_background_displacement": total_mass_coefficient-e["cap_total_mass"],
        "signed_background_angular_jump": bm0-bp0-b*x,
        "signed_background_mass_difference": bm0*bm0-bp0*bp0-2*m/x+h*h*x*x,
        "background_acceleration_agrees_on_both_sides": acceleration-bm0*bmx-h*h*x,
        "linear_angular_jump": dbm-dbp-b*rho,
        "linear_exterior_mass_constraint": bp0*dbp-v*rho_s-m*rho/x**2+mu/x,
        "linear_cap_mass_constraint": bm0*dbm-v*rho_s+h*h*x*rho+nu/x,
        "two_displacement_formulas_agree": rho_s-cap_rho_s,
        "mass_difference_transport_from_equal_scalar_stresses": mu_s-nu_s-D_s,
        "mass_difference_from_both_first_derivative_functionals": on_background(direct_difference_rate+b*x**3*e["improved_t_un"]/2),
        "temporal_jump_from_angular_and_mixed_flux": temporal-(am-ap-b)*rs,
        "normal_force_from_two_radial_einstein_equations": (radial_difference-normal_identity).subs({am: ap+b, bm: bp+b*x}),
        "linear_acceleration_from_transported_mass_difference": da-da_closed,
        "complete_exterior_acceleration_numerator_is_beta_times_delta_a": on_background(acceleration_numerator-bp0*da_closed),
        "linear_mean_curvature_variation_vanishes": da_closed+2*mean_delta_beta/x-2*mean_beta*rho/x**2,
        "conformal_material_trace_retains_tension_jump": 4*f+12*fp*fp-2/kappa,
        "zero_source_mass_ingredients": (L+Q).subs({p: 0, A: 0, B: 0}),
        "quadratic_field_scaling": (L+Q).subs({p: scale*p, A: scale*A, B: scale*B}, simultaneous=True)-scale*scale*(L+Q),
        "normal_reversal_preserves_cap_mass_ingredients": (L+Q).subs({beta: -beta, B: -B}, simultaneous=True)-(L+Q),
        "physical_mass_and_displacement_unit_scaling": (kappa*epsilon**2*(nu+rho)/ell).subs({ell: scale*ell, kappa: scale**2*kappa}, simultaneous=True)-scale*kappa*epsilon**2*(nu+rho)/ell,
    }
    if set(residuals) != set(IDENTITIES):
        raise ValueError("cap-wall exact identity inventory changed")
    for name, expression in residuals.items():
        budget()
        if sp.factor(expression) != 0:
            raise ValueError("cap-wall exact identity failed: " + name)
    controls = {
        "omitted_transformed_potential_changes_transport": L-e["cap_potential_density"]+Qs-jordan,
        "omitted_cap_einstein_tensor_changes_normal_stress": direct_nn+h*h*p*p/2-e["cap_improved_t_nn"],
        "omitted_cap_displacement_changes_total_mass": nu-total_mass_coefficient,
        "omitted_endpoint_changes_jordan_transport": L-jordan,
        "wrong_mixed_flux_sign_changes_temporal_jump": ((am*rs-x*q/2)-(ap*rs+x*q/2)-b*rs).subs(am, ap+b),
        "changed_cap_mass_without_exterior_changes_wall_forcing": sp.diff(rho_s, nu),
    }
    if set(controls) != set(CONTROLS):
        raise ValueError("cap-wall omission control inventory changed")
    for name, expression in controls.items():
        budget()
        if sp.factor(expression) == 0:
            raise ValueError("cap-wall omission control unexpectedly vanished: " + name)
    budget()
    return expected_certificate()
