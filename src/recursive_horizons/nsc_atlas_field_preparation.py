"""Finite characteristic preparation outside the fixed atlas wall.

This certifies a retained Schwarzschild-IV domain and its conformal-scalar
transfer operator. It neither solves a new waveform nor prepares a signal
from the causally excluded right parent.
"""
import math
import time

import numpy as np
from scipy.special import lambertw

from .nsc_atlas_causal_reference import reference_wall_events
from .nsc_atlas_causal_observers import de_sitter_wall_events
from .nsc_atlas_field_transfer import (
    arrival_from_wall_event, compact_bump, conformal_junction_coefficients,
)

SCHEMA = "nsc-atlas-field-preparation-v1"
M_OVER_ELL = .4
WALL_LOWER = .45
WALL_UPPER = .55
RADIAL_FLOOR = .4
WALL_SAMPLES = 65
MAX_SAMPLES = 513
CPU_SECONDS = 60


def _scale(length_scale):
    ell = float(length_scale)
    if not math.isfinite(ell) or ell <= 0:
        raise ValueError("length scale must be finite and positive")
    return ell


def radius_product(radius_over_ell):
    """UV on the white-hole branch; x is dimensionless areal radius."""
    x = np.asarray(radius_over_ell, dtype=float)
    if not np.all(np.isfinite(x)) or np.any(x <= 0) or np.any(x >= .8):
        raise ValueError("the preparation uses 0<R/ell<0.8 in sector IV")
    return (1-x/.8)*np.exp(x/.8)


def radius_from_null(U, V, length_scale=1.):
    """Principal Lambert branch in IV, returning physical areal radius.

    Signs alone do not establish retained-side membership. The wedge also
    requires U<=U_W(V), as certified by its boundary construction.
    """
    ell = _scale(length_scale)
    U, V = np.broadcast_arrays(np.asarray(U, dtype=float), np.asarray(V, dtype=float))
    product = U*V
    if (not np.all(np.isfinite(U)) or not np.all(np.isfinite(V))
            or np.any(U >= 0) or np.any(V >= 0)
            or np.any(product <= 0) or np.any(product >= 1)):
        raise ValueError("regular sector IV requires U<0, V<0 and 0<UV<1")
    return ell*.8*(1+lambertw(-product/math.e, 0).real)


def exterior_coefficient(radius_over_ell):
    """Dimensionless K in chi_UV+K chi=0, for the frozen scalar model."""
    x = np.asarray(radius_over_ell, dtype=float)
    radius_product(x)
    return 16*M_OVER_ELL**4*np.exp(-x/(2*M_OVER_ELL))/x**4


def coefficient_from_radius_derivatives(radius_over_ell):
    """Independent K=-R_UV/R from implicit differentiation of UV=F(R)."""
    x = np.asarray(radius_over_ell, dtype=float)
    F = radius_product(x)
    a = 2*M_OVER_ELL
    first = -x*np.exp(x/a)/a**2
    second = -np.exp(x/a)*(1+x/a)/a**2
    mixed = 1/first-F*second/first**3
    return -mixed/x


def source_profile(V, Va, Vb):
    """Unit-amplitude independent C3 compact polynomial, and d/dV.

    Its support occupies the middle half of the preparation characteristic.
    This polynomial is C3, not C-infinity; its first derivative is C2.
    """
    Va, Vb = float(Va), float(Vb)
    V = np.asarray(V, dtype=float)
    if (not math.isfinite(Va) or not math.isfinite(Vb) or not Va < Vb < 0
            or not np.all(np.isfinite(V)) or np.any(V < Va) or np.any(V > Vb)):
        raise ValueError("source values must lie on the admitted negative-V interval")
    z = 4*(V-(Va+Vb)/2)/(Vb-Va)
    value, derivative = compact_bump(z)
    return value, derivative*4/(Vb-Va)


def validate_preparation_inputs(*, source_corner=0., bottom_corner=0.,
                                initial_cap_field=0., initial_cap_time_derivative=0.,
                                xi=1/6, surface_scalar_coupling=0.):
    """Reject data or theories incompatible with this zero-cap certificate."""
    values = (source_corner, bottom_corner, initial_cap_field,
              initial_cap_time_derivative, xi, surface_scalar_coupling)
    if not all(math.isfinite(float(value)) for value in values):
        raise ValueError("preparation data must be finite")
    if source_corner != bottom_corner:
        raise ValueError("the two characteristic traces disagree at their corner")
    if any(value != 0 for value in values[:4]):
        raise ValueError("the certified bottom and initial cap data are zero")
    if xi != 1/6 or surface_scalar_coupling != 0:
        raise ValueError("changed coupling requires a separate transmission candidate")


def validate_preparation_domain(U0, Ua, Ub, Va, Vb):
    """Check endpoints used in the analytic retained-wedge proof."""
    values = tuple(float(value) for value in (U0, Ua, Ub, Va, Vb))
    U0, Ua, Ub, Va, Vb = values
    if not all(math.isfinite(value) for value in values):
        raise ValueError("preparation boundaries must be finite")
    if not U0 < Ua < Ub < 0 or not Va < Vb < 0:
        raise ValueError("source must precede the future timelike wall in retained IV")
    for product, radius in ((U0*Va, RADIAL_FLOOR),
                            (Ua*Va, WALL_LOWER), (Ub*Vb, WALL_UPPER)):
        if not math.isclose(product, float(radius_product(radius)), rel_tol=2e-12, abs_tol=2e-12):
            raise ValueError("boundary product does not match its declared areal radius")


def require_observer_time(time, arrival_window, length_scale=1.):
    """Reject extension of this finite preparation to uncomputed clock times."""
    ell = _scale(length_scale)
    value = float(time)
    window = np.asarray(arrival_window, dtype=float)
    if (window.shape != (2,) or not np.all(np.isfinite(window))
            or not window[0] < window[1] <= 0 or not math.isfinite(value)
            or not ell*window[0] <= value <= ell*window[1]):
        raise ValueError("observer time is outside the finite prepared response window")


def preparation_record(length_scale=1.):
    """Return a reproducible certificate and at most 513 boundary samples.

    Numerical quadrature checks locate the boundaries. Existence and no-return
    inequalities also have analytic bounds independent of quadrature accuracy.
    No field PDE, operator iteration, energy flux or central waveform is run.
    """
    ell = _scale(length_scale)
    started = time.process_time()
    validate_preparation_inputs()
    x = np.linspace(WALL_LOWER, WALL_UPPER, WALL_SAMPLES)
    wall = reference_wall_events(x, 64)
    cap = de_sitter_wall_events(ell*x, ell, 64)
    U, V = wall["U"], wall["V"]
    Ua, Ub, Va, Vb = (float(value) for value in (U[0], U[-1], V[0], V[-1]))
    U0 = float(radius_product(RADIAL_FLOOR))/Va
    validate_preparation_domain(U0, Ua, Ub, Va, Vb)
    if np.any(np.diff(U) <= 0) or np.any(np.diff(V) <= 0):
        raise ValueError("wall must be future timelike with monotone null coordinates")

    source_V = np.linspace(Va, Vb, WALL_SAMPLES)
    source_U = np.full(WALL_SAMPLES, U0)
    bottom_U = np.linspace(U0, Ua, 33)
    bottom_V = np.full(len(bottom_U), Va)
    source_value, source_derivative = source_profile(source_V, Va, Vb)
    source_radius = radius_from_null(source_U, source_V, ell)/ell
    bottom_radius = radius_from_null(bottom_U, bottom_V, ell)/ell
    arrivals, stretch = arrival_from_wall_event(ell*x, cap["t_flat"], cap["gamma"], cap["v"], ell)
    arrival_window = (arrivals[[0, -1]]/ell).tolist()

    u = (cap["eta"]-cap["r_wall"])/ell
    v = (cap["eta"]+cap["r_wall"])/ell
    no_return_margin = float(v[0]-u[-1])
    # On this interval beta_- decreases and beta_->1. Since q=Rdot^2>1/8
    # and Rdot^2=beta_-^2-(1-x^2)<beta_-^2, gamma=(1+q)/(beta_-+x Rdot)
    # is below 1+beta_-(a)^2. This bounds the entire flat-time increment.
    beta_max = 1/(5*WALL_LOWER**2)+3*WALL_LOWER/4
    flat_increment_bound = (WALL_UPPER-WALL_LOWER)*math.sqrt(8)*(1+beta_max**2)
    no_return_log_margin = math.log((1+WALL_UPPER)/(1-WALL_LOWER))-flat_increment_bound
    if no_return_margin <= 0 or no_return_log_margin <= 0:
        raise ValueError("zero outgoing wall characteristic lacks a no-return certificate")

    Fc, Fb = (float(radius_product(value)) for value in (RADIAL_FLOOR, WALL_UPPER))
    Kmax = float(exterior_coefficient(RADIAL_FLOOR))
    rectangle_bound = Kmax*(Ub-U0)*(Vb-Va)
    analytic_bound = Kmax*(math.sqrt(Fc)-math.sqrt(Fb))**2
    if not 0 < rectangle_bound <= analytic_bound < 1:
        raise ValueError("finite characteristic operator is outside the contraction certificate")
    deviation = analytic_bound/(1-analytic_bound)

    independent = reference_wall_events(np.array([WALL_LOWER, .5, WALL_UPPER]), 32)
    selection = [0, WALL_SAMPLES//2, WALL_SAMPLES-1]
    quadrature_difference = max(float(np.max(abs(wall[key][selection]-independent[key])))
                                for key in ("U", "V", "tau_wall_over_ell"))
    # Check the orientation against beta_+, independently of a picture of IV.
    Fprime = -x*np.exp(x/.8)/.8**2
    normal_R = (-wall["ell_dU_dtau"]*V+wall["ell_dV_dtau"]*U)/Fprime
    coefficient_grid = np.linspace(RADIAL_FLOOR, WALL_UPPER, 33)
    coefficient_error = float(np.max(abs(exterior_coefficient(coefficient_grid)
                                        -coefficient_from_radius_derivatives(coefficient_grid))))
    counts = {"distinct_wall_radii": len(x),
              "boundary_samples": len(x)+len(source_V)+len(bottom_U),
              "maximum_samples": MAX_SAMPLES, "cpu_seconds_max": CPU_SECONDS}
    if counts["boundary_samples"] > MAX_SAMPLES or time.process_time()-started > CPU_SECONDS:
        raise RuntimeError("finite preparation budget exhausted")

    return {
        "schema": SCHEMA,
        "status": "FINITE_RETAINED_IV_PREPARATION_AND_OPERATOR_CERTIFIED",
        "scope": "fixed isolated reference; finite conditional classical preparation, no solved waveform",
        "inputs": {"length_scale": ell, "H_ell": 1., "k_ell": 2., "M_over_ell": .4,
                   "xi": 1/6, "angular_mode": 0, "linear_probe": True,
                   "additional_surface_scalar_coupling": 0., "wall_branch": "expanding"},
        "geometry": {"wall_radius_window_over_ell": [WALL_LOWER, WALL_UPPER],
                     "radius_lower_bound_over_ell": RADIAL_FLOOR,
                     "radius_upper_bound_over_ell": WALL_UPPER,
                     "U0": U0, "Ua": Ua, "Ub": Ub, "Va": Va, "Vb": Vb,
                     "domain": "Va<=V<=Vb and U0<=U<=U_W(V)",
                     "retained_side": "U<=U_W(V), normal n=(-Udot,Vdot)",
                     "radius_proof": "F(.55)<=UV<=F(.4), F strictly decreasing in IV"},
        "preparation": {"side": "retained Schwarzschild white-hole sector IV, outside the wall",
                        "surface": "U=U0, Va<=V<=Vb",
                        "profile": "g(V)=epsilon B(4(V-Vmid)/(Vb-Va)); B(z)=(1-z^2)^4 on |z|<1, zero otherwise",
                        "regularity": "C3 compact polynomial; derivative C2; not C-infinity",
                        "support_V": [Va+(Vb-Va)/4, Va+3*(Vb-Va)/4],
                        "cap_initial_t_over_ell": float(cap["t_flat"][0]/ell),
                        "cap_data": "phi=0 and partial_t phi=0 on the flat cap slice t=t_a",
                        "bottom_data": "chi=0 on V=Va, U0<=U<=Ua",
                        "corner_compatible": True, "independent_of_target_waveform": True,
                        "remote_generation_solved": False},
        "boundary_samples": {
            "wall": {"U": U.tolist(), "V": V.tolist(), "R_over_ell": x.tolist(),
                     "t_flat_over_ell": (cap["t_flat"]/ell).tolist(),
                     "tau_wall_over_ell": (cap["tau_wall"]/ell).tolist(),
                     "arrival_t_C_over_ell": (arrivals/ell).tolist(), "dT_dtau": stretch.tolist()},
            "source": {"U": source_U.tolist(), "V": source_V.tolist(),
                       "R_over_ell": source_radius.tolist(), "chi_over_epsilon": source_value.tolist(),
                       "dchi_dV_over_epsilon": source_derivative.tolist()},
            "bottom": {"U": bottom_U.tolist(), "V": bottom_V.tolist(),
                       "R_over_ell": bottom_radius.tolist()}},
        "operator": {
            "equation": "chi_UV+K chi=0; K=16 M^4 exp(-R/(2M))/R^4=-R_UV/R",
            "wall_condition": "chi_U=0 from conformal matching and zero outgoing cap characteristic",
            "lower_limit": "L(u)=Va for u<=Ua; L(u)=inverse wall V(u) otherwise",
            "integral_equation": "chi(U,V)=g(V)-integral[U0,U]du integral[L(u),V]dv K(u,v)chi(u,v)",
            "incoming_trace": "I=2 Vdot (g'(V)-integral[U0,U_W(V)]du K(u,V)chi(u,V))",
            "central_map": "phi_C[T(tau)]=I(tau)/T'(tau), T=t_flat-ell log(1-R/ell)",
            "coefficient_max": Kmax, "rectangle_norm_bound": rectangle_bound,
            "analytic_norm_bound": analytic_bound,
            "bound_proof": "DeltaU DeltaV=Fc+Fb-Fc r-Fb/r <= (sqrt(Fc)-sqrt(Fb))^2, r=Vb/Va>0",
            "reduced_field_relative_deviation_bound": deviation,
            "midpoint_wall_reduced_field_min_over_epsilon": 1-deviation,
            "bound_scope": "uniform reduced field chi; not a central phi peak or energy flux",
            "solution_iteration_performed": False, "new_waveform_evaluated": False},
        "observer": {"arrival_window_over_ell": arrival_window,
                     "initial_cap_t_over_ell": float(cap["t_flat"][0]/ell),
                     "first_possible_incoming_label_over_ell": float(v[0]),
                     "last_outgoing_wall_label_over_ell": float(u[-1]),
                     "no_return_margin_over_ell": no_return_margin,
                     "flat_time_increment_upper_bound_over_ell": flat_increment_bound,
                     "analytic_no_return_log_margin_lower_bound": no_return_log_margin,
                     "centre_regularity": "chi=F(eta+r)-F(eta-r); phi_C=2F'(eta_C)/a_C",
                     "after_window": "uncomputed; no zero tail or extension through t_C=0 inferred"},
        "access": [
            {"from": "P", "to": "W or C", "verdict": "excluded_in_reference",
             "basis": "P has V>0, wall has V<0; future causal V cannot decrease",
             "meaning": "no driven transfer from the designated right parent, not an exclusion of initial correlations"},
            {"from": "Gamma_IV", "to": "W then C", "verdict": "finite_conditional_preparation_admitted",
             "basis": "constant-V future rays remain in the retained wedge to the wall; R/ell<1 admits centre reception",
             "field_amplitude": "bounded operator established; new waveform not evaluated"}],
        "checks": {"kernel_geometric_identity_max_error": coefficient_error,
                   "normal_orientation_max_error": float(np.max(abs(normal_R-wall["beta_plus"]))),
                   "wall_metric_norm_max_error": float(np.max(abs(wall["metric_norm"]+1))),
                   "minimum_ell_dU_dtau": float(np.min(wall["ell_dU_dtau"])),
                   "minimum_ell_dV_dtau": float(np.min(wall["ell_dV_dtau"])),
                   "quadrature_32_64_difference": quadrature_difference,
                   "quadrature_is_rigorous_enclosure": False,
                   "analytic_contraction_independent_of_quadrature": True,
                   "analytic_no_return_independent_of_quadrature": True,
                   "reduced_field_jump_coefficient": conformal_junction_coefficients(2/ell)["normal_chi_jump_per_chi"]},
        "budget": counts,
        "source_reading": [
            "1106.4318v1 sections II-III, scalar radial operator, especially equation 11",
            "1408.4995v2 section4.2 Theorem23: characteristic context and hypotheses, not a thin-wall theorem",
            "existing v5 signed Kruskal wall and v6 conformal transmission / regular-centre construction"],
        "limits": ["conditional IV initial data, not a source on the right-parent worldline",
                   "finite classical probe; no global vacuum, source generation, energy flux or backreaction",
                   "no new waveform, numerical PDE solution, formation or perturbative stability result",
                   "not the causal development of the glued many-body parent",
                   "later scattering and centre returns require a new domain and boundary account"]}
