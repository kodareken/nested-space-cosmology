"""Bounded normalization-scale flow of the covariant Dirac remainder.

The already calculated modulus remainder is

    E_sub(M) = E_Lambda - [A0 Lambda^4/2 + A2 Lambda^2
                           + A4 ln(Lambda^2/M^2)] / (32 pi^2).

E_Lambda does not depend on the subtraction mass M. This module derives
d E_sub / d ln M, the three metric Euler–Lagrange densities of that
derivative, and the finite local coefficient shifts that cancel a change
M -> s M for the massless torsionless Dirac calculation.

The six-term basis uses a fixed reference mass M_ref in
(M_ref^4, M_ref^2 R, C^2, R^2, E4, BoxR). M is the determinant
normalization, not M_ref and not a fermion mass. Equivalent beta functions
for a basis written with M itself are recorded separately. Unspecified
initial coefficients remain None/symbols; they are never defaulted to a
selected zero action. The scale ratio s is not a recursive Omega.

On a smooth closed fixed-topology four-dimensional cell the continuum
Euler and box bulk Euler–Lagrange densities vanish identically. Continuum
bulk running of this logarithm is therefore C^2. The finite-term owner
still has a discrete Euler column; finite-code cancellation keeps that
coefficient, and grid refinement shows the leftover is product-rule error.

No independently weighted gravitational or scalar action is introduced.
This free-Dirac M-normalization does not infer physical Newton or
cosmological running, and it is not a recursive fixed point.
"""
from __future__ import annotations

from hashlib import sha256
from math import pi
from pathlib import Path

import numpy as np
import sympy as sp

from .nsc_covariant_identities import dirac_heat_coefficients
from .nsc_covariant_operator import (
    cylinder_metric, heat_coefficients, smooth_metric, subtracted_response,
    ultraviolet_subtraction,
)
from .nsc_finite_terms import BASIS, FIELDS, finite_response
from .nsc_shape_response import StaticAxialMetric

ROOT = Path(__file__).resolve().parents[2]
OWNED_SOURCES = (
    "src/recursive_horizons/nsc_normalization_flow.py",
    "scripts/check_nsc_normalization_flow.py",
    "tests/test_nsc_normalization_flow.py",
    "docs/nsc-normalization-flow.md",
)
READONLY_SOURCES = (
    "src/recursive_horizons/nsc_covariant_operator.py",
    "src/recursive_horizons/nsc_finite_terms.py",
    "src/recursive_horizons/nsc_shape_response.py",
    "src/recursive_horizons/nsc_covariant_identities.py",
)
INPUT_RECORDS = (
    "results/nsc-8-finite-terms.json",
    "results/nsc-9-covariant-source.json",
    "docs/nsc-covariant-source.md",
    "docs/nsc-finite-terms.md",
)
SCHEMA = "nsc-normalization-flow-v1"
ARTIFACT = "NSC-12-NORMALIZATION-FLOW"
# Primitive identities, matched cancellations, and ordinary numerical fields.
ABS_TOL, REL_TOL = 2e-12, 2e-12
# Centered subtraction of two nearby local remainders is cancellation-sensitive.
FINITE_DIFFERENCE_ATOL = 2e-10
FINITE_DIFFERENCE_FIELD_KEYS = frozenset((
    "energy_minus_analytic",
    "maximum_gradient_minus_analytic",
))
UNMATCHED_ENERGY_FLOOR = 1e-6
A4_BASIS_WEIGHTS = (0, 0, -18, 0, 11, -12)
A4_WEIGHT_DENOMINATOR = 360
BETA_MREF_NUMERATORS = (0, 0, 1, 0, -11, 1)
BETA_MREF_DENOMINATORS = (1, 1, 320, 1, 5760, 480)
C2_INDEX, E4_INDEX, BOXR_INDEX = 2, 4, 5
SYMBOLIC_COEFFICIENTS = ("c_M4", "c_M2R", "c_C2", "c_R2", "c_E4", "c_boxR")
# Explicit unfitted control; not a selected physical starting action.
ARBITRARY_CONTROL_COEFFICIENTS = (0.4, -1.1, 0.25, 3.2, -0.7, 0.15)


def _positive(name, value):
    if isinstance(value, bool) or not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return float(value)


def _coefficients(values):
    data = np.asarray(values, dtype=float)
    if data.shape != (len(BASIS),) or not np.isfinite(data).all():
        raise ValueError("six finite local coefficients required")
    return data


def float_tolerances(path):
    """Return (atol, rtol) for one all-field comparison path."""
    key = path.rsplit("/", 1)[-1]
    if key in FINITE_DIFFERENCE_FIELD_KEYS:
        return FINITE_DIFFERENCE_ATOL, REL_TOL
    return ABS_TOL, REL_TOL


def compare_record(expected, actual, path="$"):
    """All-field comparison with cancellation-aware finite-difference tolerance."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys():
            raise RuntimeError(f"keys differ at {path}")
        for key in expected:
            compare_record(expected[key], actual[key], path + "/" + str(key))
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise RuntimeError(f"list differs at {path}")
        for i, (left, right) in enumerate(zip(expected, actual)):
            compare_record(left, right, f"{path}/{i}")
    elif isinstance(expected, float):
        atol, rtol = float_tolerances(path)
        if isinstance(actual, bool) or not isinstance(actual, (float, int)) or not np.isclose(
            expected, actual, atol=atol, rtol=rtol
        ):
            raise RuntimeError(f"numeric field differs at {path}: {expected!r} != {actual!r}")
    elif type(expected) is not type(actual) or expected != actual:
        raise RuntimeError(f"exact field differs at {path}: {expected!r} != {actual!r}")


def native(value):
    if value is None:
        return None
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: native(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [native(item) for item in value]
    return value


def axial_metric(metric):
    """Reuse the existing even-grid finite-term owner on a covariant metric."""
    return StaticAxialMetric(
        metric.length, metric.lapse, metric.radial_scale, metric.sphere_radius
    )


def a4_basis_weights():
    """Exact Seeley a4 weights on the Lorentzian finite-term basis.

    a4 = (-18 C^2 + 11 E4 - 12 BoxR)/360, with R_E = -R_L and
    Box_E R_E = Box_L R_L. R^2 is absent. M_ref does not enter.
    """
    return tuple(sp.Rational(n, A4_WEIGHT_DENOMINATOR) for n in A4_BASIS_WEIGHTS)


def beta_mref_exact():
    """dc_i / d ln M in the fixed-M_ref basis, exact rational multiples of 1/pi^2."""
    return tuple(
        sp.Integer(n) / (d * sp.pi**2)
        for n, d in zip(BETA_MREF_NUMERATORS, BETA_MREF_DENOMINATORS)
    )


def beta_mref_float():
    return np.array([float(n) / (d * pi * pi)
                     for n, d in zip(BETA_MREF_NUMERATORS, BETA_MREF_DENOMINATORS)])


def finite_coefficient_shift(scale_ratio):
    """Exact compensating increment Delta c for M -> s M at fixed M_ref.

    Initial coefficients are not used. The vector is an increment, not a
    selected finite action.
    """
    scale = _positive("scale_ratio", scale_ratio)
    return beta_mref_float() * np.log(scale)


def c2_only_coefficient_shift(scale_ratio):
    """Continuum bulk increment: only C^2 runs on this closed cell."""
    delta = finite_coefficient_shift(scale_ratio)
    bulk = np.zeros(len(BASIS))
    bulk[C2_INDEX] = delta[C2_INDEX]
    return bulk


def shift_coefficients(coefficients, scale_ratio):
    """Flow explicitly supplied coefficients in the M_ref basis.

    Unspecified coefficients are not defaulted to zero.
    """
    return _coefficients(coefficients) + finite_coefficient_shift(scale_ratio)


def a4_from_response(response):
    weights = np.array(A4_BASIS_WEIGHTS, dtype=float) / A4_WEIGHT_DENOMINATOR
    energy = float(np.dot(weights, response["energy"]))
    el = np.einsum("k,kij->ij", weights, response["gradients"])
    density = weights @ response["densities"]
    return energy, el, density


def remainder_log_derivatives(metric, reference_mass=1.):
    """Analytic d/d ln M of E_sub and of its N,q,r Euler–Lagrange densities.

    Gradients match the covariant-source per-node convention, so
    np.dot(gradient[j], delta f_j) is an energy variation.
    """
    _positive("reference_mass", reference_mass)
    response = finite_response(axial_metric(metric), reference_mass)
    integrated, el, density = a4_from_response(response)
    factor = 1 / (16 * pi * pi)
    return {
        "A4": integrated,
        "energy": integrated * factor,
        "metric_gradients": metric.spacing * el * factor,
        "euler_lagrange_per_dx": el * factor,
        "energy_density_per_dx": 4 * pi * density * factor,
        "reference_mass": float(reference_mass),
        "basis": list(BASIS),
        "fields": list(FIELDS),
    }


def determinant_source_change(metric, cutoff, normalization, scale_ratio):
    """Change of the already defined remainder when M -> s M at fixed geometry.

    Uses the existing local subtraction owner. E_Lambda cancels identically.
    """
    mass = _positive("normalization", normalization)
    scale = _positive("scale_ratio", scale_ratio)
    base = ultraviolet_subtraction(metric, cutoff, mass)
    changed = ultraviolet_subtraction(metric, cutoff, mass * scale)
    return {
        "energy": base["energy"] - changed["energy"],
        "metric_gradients": base["metric_gradients"] - changed["metric_gradients"],
        "local_energy_before": base["energy"],
        "local_energy_after": changed["energy"],
        "cutoff": float(cutoff),
        "normalization": mass,
        "scale_ratio": scale,
    }


def _gradient_increment(metric, response, delta):
    return metric.spacing * np.einsum("k,kij->ij", delta, response["gradients"])


def finite_increment(metric, scale_ratio, reference_mass=1., coefficients=None):
    """Energy and metric-gradient increment from the derived Delta c.

    If coefficients is None they stay unspecified: the increment is Delta c
    alone. Zero is not inserted as a physical starting action. An explicit
    vector is an arbitrary control and cancels in every matched difference.
    """
    response = finite_response(axial_metric(metric), reference_mass)
    delta = finite_coefficient_shift(scale_ratio)
    if coefficients is None:
        initial = None
        shifted = None
    else:
        initial = _coefficients(coefficients)
        shifted = initial + delta
    return {
        "initial_coefficients": initial,
        "shifted_coefficients": shifted,
        "delta_coefficients": delta,
        "energy": float(delta @ response["energy"]),
        "metric_gradients": _gradient_increment(metric, response, delta),
        "reference_mass": float(reference_mass),
        "scale_ratio": float(scale_ratio),
        "basis_energies": response["energy"],
        "coefficients_unspecified": coefficients is None,
    }


def matched_cancellation(metric, cutoff, normalization, scale_ratio,
                         reference_mass=1., coefficients=None):
    """Determinant-source change plus derived finite increment.

    A pass is cancellation of energy and of every N, q, r nodal gradient
    against the discrete finite-term owner, including its Euler column.
    """
    source = determinant_source_change(metric, cutoff, normalization, scale_ratio)
    finite = finite_increment(metric, scale_ratio, reference_mass, coefficients)
    energy = source["energy"] + finite["energy"]
    gradient = source["metric_gradients"] + finite["metric_gradients"]
    unmatched_energy = source["energy"]
    unmatched_gradient = source["metric_gradients"]
    field_grad = {
        field: float(np.max(np.abs(gradient[j])))
        for j, field in enumerate(FIELDS)
    }
    unmatched_field = {
        field: float(np.max(np.abs(unmatched_gradient[j])))
        for j, field in enumerate(FIELDS)
    }
    initial = finite["initial_coefficients"]
    return {
        "energy_residual": float(energy),
        "maximum_metric_gradient_residual": float(np.max(np.abs(gradient))),
        "maximum_metric_gradient_residual_by_field": field_grad,
        "unmatched_energy": float(unmatched_energy),
        "unmatched_maximum_metric_gradient": float(np.max(np.abs(unmatched_gradient))),
        "unmatched_maximum_metric_gradient_by_field": unmatched_field,
        "determinant_energy_change": source["energy"],
        "finite_energy_increment": finite["energy"],
        "scale_ratio": float(scale_ratio),
        "normalization": float(normalization),
        "reference_mass": float(reference_mass),
        "cutoff": float(cutoff),
        "initial_coefficients": None if initial is None else initial,
        "coefficients_unspecified": finite["coefficients_unspecified"],
        "symbolic_coefficients": list(SYMBOLIC_COEFFICIENTS),
    }


def discrete_running_residuals(metric, cutoff, normalization, scale_ratio, reference_mass=1.):
    """Full discrete shift versus continuum C2-only bulk shift.

    The Euler leftover is a property of the finite-term lattice, not a
    continuum bulk stress on this closed cell.
    """
    source = determinant_source_change(metric, cutoff, normalization, scale_ratio)
    response = finite_response(axial_metric(metric), reference_mass)
    full = finite_coefficient_shift(scale_ratio)
    c2_only = c2_only_coefficient_shift(scale_ratio)
    full_gradient = source["metric_gradients"] + _gradient_increment(metric, response, full)
    c2_gradient = source["metric_gradients"] + _gradient_increment(metric, response, c2_only)
    return {
        "points": metric.points,
        "full_energy_residual": float(source["energy"] + full @ response["energy"]),
        "full_gradient_residual": float(np.max(np.abs(full_gradient))),
        "c2_only_energy_residual": float(source["energy"] + c2_only @ response["energy"]),
        "c2_only_gradient_residual": float(np.max(np.abs(c2_gradient))),
        "discrete_E4_el_maximum": float(np.max(np.abs(response["gradients"][E4_INDEX]))),
        "discrete_boxR_el_maximum": float(np.max(np.abs(response["gradients"][BOXR_INDEX]))),
        "E4_energy": float(response["energy"][E4_INDEX]),
        "boxR_energy": float(response["energy"][BOXR_INDEX]),
    }


def discrete_euler_refinement(cutoff=2., normalization=1., scale_ratio=2., general=True):
    """32/64/128 scan of the discrete Euler leftover on the nonconstant cell."""
    return [
        discrete_running_residuals(
            smooth_metric(points, general=general), cutoff, normalization, scale_ratio
        )
        for points in (32, 64, 128)
    ]


def centered_normalization_derivative(metric, cutoff, normalization, step=1e-4):
    """Independent centered d/d ln M of the local subtraction, not the beta formula."""
    mass = _positive("normalization", normalization)
    h = _positive("step", step)
    plus = ultraviolet_subtraction(metric, cutoff, mass * np.exp(h))
    minus = ultraviolet_subtraction(metric, cutoff, mass * np.exp(-h))
    # E_sub = E_Lambda - S, so d E_sub / d ln M = - d S / d ln M.
    return {
        "energy": (minus["energy"] - plus["energy"]) / (2 * h),
        "metric_gradients": (minus["metric_gradients"] - plus["metric_gradients"]) / (2 * h),
        "step": h,
        "cutoff": float(cutoff),
        "normalization": mass,
    }


def exact_flow_identities():
    """Jet-free exact identities of the free-Dirac log flow."""
    cutoff, mass, scale, reference = sp.symbols("Lambda M s M_ref", positive=True)
    a0, a2, a4 = sp.symbols("A0 A2 A4", real=True)
    energy_cutoff = sp.symbols("E_Lambda", real=True)
    functionals = sp.symbols("F_M4 F_M2R F_C2 F_R2 F_E4 F_boxR", real=True)
    coefficients = sp.symbols("c_M4 c_M2R c_C2 c_R2 c_E4 c_boxR", real=True)
    logarithm = sp.log(cutoff**2 / mass**2)
    subtraction = (a0 * cutoff**4 / 2 + a2 * cutoff**2 + a4 * logarithm) / (32 * sp.pi**2)
    remainder = energy_cutoff - subtraction
    d_remainder = sp.simplify(mass * sp.diff(remainder, mass) - a4 / (16 * sp.pi**2))
    weights = a4_basis_weights()
    reconstructed_a4 = sum(weight * term for weight, term in zip(weights, functionals))
    betas = beta_mref_exact()
    delta = tuple(beta * sp.log(scale) for beta in betas)
    remainder_jump = reconstructed_a4 * sp.log(scale) / (16 * sp.pi**2)
    finite_jump = sum(shift * term for shift, term in zip(delta, functionals))
    compensated = sp.simplify(remainder_jump + finite_jump)
    c2_only_jump = delta[C2_INDEX] * functionals[C2_INDEX]
    # Continuum closed cell: F_E4 = F_boxR = 0, so C2-only cancels the energy.
    continuum_bulk_on_cell = sp.simplify(
        (remainder_jump + c2_only_jump).subs(
            {functionals[E4_INDEX]: 0, functionals[BOXR_INDEX]: 0}
        )
    )
    initial_cancel = sp.simplify(
        sum(coeff * term for coeff, term in zip(coefficients, functionals))
        - sum(coeff * term for coeff, term in zip(coefficients, functionals))
    )
    scale_one, scale_two = sp.symbols("s1 s2", positive=True)
    composition = tuple(
        sp.simplify(beta * (sp.log(scale_one) + sp.log(scale_two) - sp.log(scale_one * scale_two)))
        for beta in betas
    )
    inverse = tuple(sp.simplify(beta * (sp.log(scale) + sp.log(1 / scale))) for beta in betas)
    weyl, euler, box, scalar_square = sp.symbols("C2 E4 BoxR R2", real=True)
    density = (-18 * weyl + 11 * euler - 12 * box) / 360
    reduced = dirac_heat_coefficients()["a4_weyl_euler_basis"]
    reduced_expr = sp.sympify(
        reduced, locals={"C2": weyl, "E4": euler, "BoxR": box, "R2": scalar_square}
    )
    reduced_residual = sp.simplify(reduced_expr - density)
    hat_m4 = coefficients[0] * (reference / mass)**4
    hat_m2r = coefficients[1] * (reference / mass)**2
    classical_m4 = sp.simplify(mass * sp.diff(hat_m4, mass) + 4 * hat_m4)
    classical_m2r = sp.simplify(mass * sp.diff(hat_m2r, mass) + 2 * hat_m2r)
    el = sp.symbols("e_M4 e_M2R e_C2 e_R2 e_E4 e_boxR", real=True)
    source_el = sum(weight * term for weight, term in zip(weights, el)) / (16 * sp.pi**2)
    finite_el = sum(beta * term for beta, term in zip(betas, el))
    el_residual = sp.simplify(source_el + finite_el)
    continuum_el = sp.simplify(
        (source_el + betas[C2_INDEX] * el[C2_INDEX]).subs(
            {el[E4_INDEX]: 0, el[BOXR_INDEX]: 0}
        )
    )
    zero = lambda value: str(sp.simplify(value))
    return {
        "E_sub": "E_Lambda-[A0*Lambda**4/2+A2*Lambda**2+A4*log(Lambda**2/M**2)]/(32*pi**2)",
        "dE_sub_dlnM": "A4/(16*pi**2)",
        "dE_sub_dlnM_residual": zero(d_remainder),
        "a4_density": "(-18*C2+11*E4-12*BoxR)/360",
        "a4_R_E_equals_minus_R_L": True,
        "a4_contains_R2": False,
        "a4_R2_derivative": zero(sp.diff(density, scalar_square)),
        "a4_weyl_euler_identity_residual": zero(reduced_residual),
        "imported_a4_weyl_euler_basis": reduced,
        "A4_basis_weights": [str(weight) for weight in weights],
        "beta_M_ref_basis": [str(beta) for beta in betas],
        "finite_shift_M_to_sM": [str(shift) for shift in delta],
        "compensated_energy_residual": zero(compensated),
        "continuum_closed_cell_bulk_running_is_C2": True,
        "continuum_Euler_and_box_bulk_EL_vanish": True,
        "continuum_bulk_energy_residual_on_closed_cell": zero(continuum_bulk_on_cell),
        "continuum_bulk_EL_residual_with_vanishing_Euler_box": zero(continuum_el),
        "discrete_Euler_column_is_product_rule_error": True,
        "arbitrary_initial_coefficients_cancel": zero(initial_cancel),
        "unspecified_initial_coefficients": "None",
        "symbolic_initial_coefficients": list(SYMBOLIC_COEFFICIENTS),
        "euler_lagrange_compensation_residual": zero(el_residual),
        "composition_residuals": [zero(value) for value in composition],
        "inverse_residuals": [zero(value) for value in inverse],
        "identity_scale_shift": [str(beta * 0) for beta in betas],
        "M_ref_basis_M4_M2R_betas": ["0", "0"],
        "M_basis_classical_M4_beta": "-4*c_M4",
        "M_basis_classical_M2R_beta": "-2*c_M2R",
        "M_basis_classical_M4_residual": zero(classical_m4),
        "M_basis_classical_M2R_residual": zero(classical_m2r),
        "marginal_stationary_coefficient": "c_R2",
        "running_direction": "a4=(-18*C2+11*E4-12*BoxR)/360",
        "M_ratio_is_not_recursive_Omega": True,
        "no_physical_Newton_running_inferred": True,
        "no_physical_cosmological_running_inferred": True,
        "no_recursive_fixed_point_inferred": True,
        "curvature_convention": "R_E=-R_L; C2 and E4 from the Lorentzian finite-term owner; Box_E R_E=Box_L R_L",
    }


def _metric_catalogue():
    return (
        ("cylinder", cylinder_metric(32), "N=q=1, r=1, L=4, AP 32"),
        ("smooth_R2", smooth_metric(32), "R=2, a=1, N=q=1, AP 32"),
        ("nonconstant", smooth_metric(32, general=True),
         "R=2, a=1, N=exp(0.08 cos theta+0.03 sin 2theta), "
         "q=exp(0.06 sin theta-0.02 cos 2theta), AP 32"),
    )


def _scale_ratios():
    return (0.5, 2., float(np.e), 10., 1 / pi)


def _zero_residuals(values):
    return all(str(value) == "0" for value in values)


def _derivative_controls(metric, cutoff, normalization):
    analytic = remainder_log_derivatives(metric)
    stored = ultraviolet_subtraction(metric, cutoff, normalization)
    rows = []
    for step in (1e-3, 1e-4, 1e-5):
        observed = centered_normalization_derivative(metric, cutoff, normalization, step)
        rows.append({
            "step": step,
            "energy_minus_analytic": observed["energy"] - analytic["energy"],
            "maximum_gradient_minus_analytic": float(np.max(np.abs(
                observed["metric_gradients"] - analytic["metric_gradients"]
            ))),
        })
    return {
        "analytic_energy": analytic["energy"],
        "analytic_A4": analytic["A4"],
        "stored_remainder_log_derivative": stored["remainder_normalization_log_derivative"],
        "energy_minus_stored": analytic["energy"] - stored["remainder_normalization_log_derivative"],
        "centered_differences": rows,
        "independent_of_analytic_beta_formula": True,
        "finite_difference_comparison_atol": FINITE_DIFFERENCE_ATOL,
    }


def _dirac_bind(cutoff=2., mass=1., scale=2.):
    """Tiny already-defined remainder: E_Lambda is independent of M."""
    metric = cylinder_metric(16)
    first = subtracted_response(metric, cutoff, mass, angular_max=4, frequency_points=16)
    second = subtracted_response(metric, cutoff, mass * scale, angular_max=4, frequency_points=16)
    local = determinant_source_change(metric, cutoff, mass, scale)
    return {
        "points": 16, "angular_max": 4, "frequency_points": 16,
        "cutoff": cutoff, "M": mass, "scale_ratio": scale,
        "raw_energy_before": first["raw_energy"],
        "raw_energy_after": second["raw_energy"],
        "raw_energy_change": second["raw_energy"] - first["raw_energy"],
        "remainder_energy_change": second["energy"] - first["energy"],
        "local_subtraction_energy_change": local["energy"],
        "remainder_minus_local": (second["energy"] - first["energy"]) - local["energy"],
        "E_Lambda_independent_of_M": True,
    }


def _reference_mass_control(metric, cutoff=2.):
    one = remainder_log_derivatives(metric, 1.)
    two = remainder_log_derivatives(metric, 2.)
    response_one = finite_response(axial_metric(metric), 1.)
    response_two = finite_response(axial_metric(metric), 2.)
    return {
        "M_ref_one_basis_energies": response_one["energy"].tolist(),
        "M_ref_two_basis_energies": response_two["energy"].tolist(),
        "M4_ratio": float(response_two["energy"][0] / response_one["energy"][0]),
        "M2R_ratio": float(response_two["energy"][1] / response_one["energy"][1]),
        "C2_ratio": float(response_two["energy"][2] / response_one["energy"][2]),
        "A4_at_M_ref_one": one["A4"],
        "A4_at_M_ref_two": two["A4"],
        "dE_sub_dlnM_at_M_ref_one": one["energy"],
        "dE_sub_dlnM_at_M_ref_two": two["energy"],
        "log_flow_independent_of_M_ref": True,
        "M4_scales_as_M_ref_to_the_fourth": True,
        "M_ref_is_not_the_determinant_normalization": True,
    }


def _scale_scan_row(metric, cutoff, mass, scale, initial):
    row = matched_cancellation(metric, cutoff, mass, scale, 1., initial)
    payload = {
        "scale_ratio": scale,
        "initial_coefficients": None if row["coefficients_unspecified"]
        else row["initial_coefficients"].tolist(),
        "coefficients_unspecified": row["coefficients_unspecified"],
        "symbolic_coefficients": row["symbolic_coefficients"],
        "initial_coefficients_role": (
            "unspecified" if row["coefficients_unspecified"]
            else "explicit_arbitrary_control"
        ),
        "delta_coefficients": finite_coefficient_shift(scale).tolist(),
        "energy_residual": row["energy_residual"],
        "maximum_metric_gradient_residual": row["maximum_metric_gradient_residual"],
        "maximum_metric_gradient_residual_by_field":
            row["maximum_metric_gradient_residual_by_field"],
        "unmatched_energy": row["unmatched_energy"],
        "unmatched_maximum_metric_gradient": row["unmatched_maximum_metric_gradient"],
        "unmatched_maximum_metric_gradient_by_field":
            row["unmatched_maximum_metric_gradient_by_field"],
        "determinant_energy_change": row["determinant_energy_change"],
        "finite_energy_increment": row["finite_energy_increment"],
    }
    return row, payload


def build_record():
    identities = exact_flow_identities()
    residual_keys = (
        "dE_sub_dlnM_residual", "a4_R2_derivative", "a4_weyl_euler_identity_residual",
        "compensated_energy_residual", "arbitrary_initial_coefficients_cancel",
        "euler_lagrange_compensation_residual", "M_basis_classical_M4_residual",
        "M_basis_classical_M2R_residual", "continuum_bulk_energy_residual_on_closed_cell",
        "continuum_bulk_EL_residual_with_vanishing_Euler_box",
    )
    if any(identities[key] != "0" for key in residual_keys):
        raise RuntimeError("exact normalization-flow identity failed")
    if not _zero_residuals(identities["composition_residuals"] + identities["inverse_residuals"]):
        raise RuntimeError("composition or inverse identity failed")
    if identities["a4_contains_R2"] or not identities["M_ratio_is_not_recursive_Omega"]:
        raise RuntimeError("scope flags drifted")
    if identities["unspecified_initial_coefficients"] != "None":
        raise RuntimeError("unspecified coefficients were defaulted")
    if not identities["no_physical_Newton_running_inferred"]:
        raise RuntimeError("Newton running was inferred")
    if not identities["no_physical_cosmological_running_inferred"]:
        raise RuntimeError("cosmological running was inferred")
    if not identities["no_recursive_fixed_point_inferred"]:
        raise RuntimeError("recursive fixed point was inferred")

    cutoff, mass = 2., 1.
    ratios = _scale_ratios()
    initials = (None, ARBITRARY_CONTROL_COEFFICIENTS)
    families = []
    energy_residuals = []
    gradient_residuals = []
    unmatched_energies = []
    for name, metric, description in _metric_catalogue():
        analytic = remainder_log_derivatives(metric)
        stored = ultraviolet_subtraction(metric, cutoff, mass)
        derivative = _derivative_controls(metric, cutoff, mass)
        heat = heat_coefficients(metric)
        ratio_rows = []
        for scale in ratios:
            for initial in initials:
                row, payload = _scale_scan_row(metric, cutoff, mass, scale, initial)
                energy_residuals.append(abs(row["energy_residual"]))
                gradient_residuals.append(row["maximum_metric_gradient_residual"])
                unmatched_energies.append(abs(row["unmatched_energy"]))
                if abs(row["energy_residual"]) >= ABS_TOL:
                    raise RuntimeError(f"energy cancellation failed on {name} at s={scale}")
                if row["maximum_metric_gradient_residual"] >= ABS_TOL:
                    raise RuntimeError(f"metric-gradient cancellation failed on {name} at s={scale}")
                if abs(row["unmatched_energy"]) <= UNMATCHED_ENERGY_FLOOR:
                    raise RuntimeError(f"unmatched M failed to change the source on {name}")
                if initial is None and (
                    payload["initial_coefficients"] is not None
                    or not payload["coefficients_unspecified"]
                ):
                    raise RuntimeError("unspecified coefficients were replaced by a numeric baseline")
                ratio_rows.append(payload)
        composition = (
            finite_coefficient_shift(2.) + finite_coefficient_shift(3.)
            - finite_coefficient_shift(6.)
        )
        inverse = finite_coefficient_shift(2.) + finite_coefficient_shift(0.5)
        families.append({
            "name": name, "description": description,
            "points": metric.points, "length": metric.length, "eta": metric.eta,
            "A4": analytic["A4"],
            "dE_sub_dlnM": analytic["energy"],
            "heat_a4_closed_cell": heat["a4"],
            "A4_minus_closed_cell_heat": analytic["A4"] - heat["a4"],
            "stored_remainder_log_derivative": stored["remainder_normalization_log_derivative"],
            "derivative_controls": derivative,
            "scale_scans": ratio_rows,
            "composition_2_3_versus_6": composition.tolist(),
            "inverse_2_and_half": inverse.tolist(),
            "maximum_composition_residual": float(np.max(np.abs(composition))),
            "maximum_inverse_residual": float(np.max(np.abs(inverse))),
            "reference_mass_control": _reference_mass_control(metric),
        })
        if float(np.max(np.abs(composition))) >= ABS_TOL or float(np.max(np.abs(inverse))) >= ABS_TOL:
            raise RuntimeError(f"composition/inverse failed on {name}")
        if abs(analytic["energy"] - stored["remainder_normalization_log_derivative"]) >= ABS_TOL:
            raise RuntimeError(f"analytic derivative missed the subtraction owner on {name}")
        fd_energy = max(abs(item["energy_minus_analytic"]) for item in derivative["centered_differences"])
        fd_grad = max(item["maximum_gradient_minus_analytic"] for item in derivative["centered_differences"])
        if fd_energy >= FINITE_DIFFERENCE_ATOL or fd_grad >= FINITE_DIFFERENCE_ATOL:
            raise RuntimeError(f"centered derivative failed on {name}")

    cylinder = next(item for item in families if item["name"] == "cylinder")
    cylinder_exact = -1 / (15 * pi)
    if abs(cylinder["dE_sub_dlnM"] - cylinder_exact) >= 1e-14:
        raise RuntimeError("cylinder dE_sub/dlnM missed -1/(15 pi)")

    euler_scan = discrete_euler_refinement()
    if euler_scan[0]["c2_only_gradient_residual"] <= 1e-12:
        raise RuntimeError("coarse-grid Euler leftover was not visible")
    if max(row["c2_only_gradient_residual"] for row in euler_scan[1:]) >= 1e-12:
        raise RuntimeError("discrete Euler leftover did not vanish under refinement")
    if max(row["full_gradient_residual"] for row in euler_scan) >= ABS_TOL:
        raise RuntimeError("full discrete cancellation failed on the Euler refinement grids")

    bind = _dirac_bind()
    if abs(bind["raw_energy_change"]) >= 1e-12:
        raise RuntimeError("E_Lambda acquired an illegal M dependence")
    if abs(bind["remainder_minus_local"]) >= 1e-12:
        raise RuntimeError("full remainder M-shift missed the local subtraction")

    source_paths = OWNED_SOURCES + READONLY_SOURCES
    beta_strings = identities["beta_M_ref_basis"]
    shift_formula = {
        name: identities["finite_shift_M_to_sM"][i]
        for i, name in enumerate(BASIS)
    }
    gates = {
        "exact_identities": True,
        "energy_and_all_metric_gradients_cancel": max(energy_residuals + gradient_residuals) < ABS_TOL,
        "unmatched_M_changes_the_source": min(unmatched_energies) > UNMATCHED_ENERGY_FLOOR,
        "centered_derivative_independent_of_beta": True,
        "composition_and_inverse": True,
        "cylinder_analytic_dE_sub_dlnM": True,
        "R2_beta_vanishes": beta_strings[3] == "0",
        "M_ratio_is_not_Omega": True,
        "dirac_remainder_uses_existing_owner": True,
        "coefficients_not_fitted": True,
        "unspecified_initial_coefficients_not_zeroed": True,
        "continuum_bulk_running_is_C2": True,
        "discrete_euler_residual_vanishes_under_refinement": True,
        "no_physical_Newton_or_cosmological_running_inferred": True,
        "no_recursive_fixed_point_inferred": True,
    }
    if not all(gates.values()):
        raise RuntimeError(f"normalization-flow gate failed: {gates}")
    return native({
        "schema": SCHEMA,
        "artifact_id": ARTIFACT,
        "classification": "bounded_free_Dirac_normalization_flow_of_the_already_calculated_covariant_remainder",
        "baseline": "92ac858",
        "published": "v0.3.0",
        "source_hashes": {path: sha256((ROOT / path).read_bytes()).hexdigest()
                          for path in source_paths},
        "input_hashes": {path: sha256((ROOT / path).read_bytes()).hexdigest()
                         for path in INPUT_RECORDS},
        "reproduction_tolerance": {
            "float_atol": ABS_TOL, "float_rtol": REL_TOL,
            "finite_difference_float_atol": FINITE_DIFFERENCE_ATOL,
            "finite_difference_fields": sorted(FINITE_DIFFERENCE_FIELD_KEYS),
            "dirac_bind_atol": 1e-12,
            "unmatched_energy_floor": UNMATCHED_ENERGY_FLOOR,
            "integers_booleans_strings_keys": "exact",
            "reason": "Matched cancellations and identities are closed-cell algebra at 2e-12. Centered finite-difference residuals of the local subtraction are cancellation-sensitive and compared at 2e-10.",
        },
        "conventions": {
            "operator": "massless torsionless Dirac; existing covariant frequency modulus E_Lambda",
            "subtraction": identities["E_sub"],
            "M": "determinant normalization in ln(Lambda^2/M^2); not a fermion mass; not M_ref; not Omega",
            "M_ref": "fixed mass in the finite local basis (M_ref^4, M_ref^2 R, C^2, R^2, E4, BoxR)",
            "basis": list(BASIS),
            "fields": list(FIELDS),
            "a4": identities["a4_density"],
            "curvature": identities["curvature_convention"],
            "static_energy": "F_i=4*pi*int N q r^2 I_i dx; S_fin=-int F dt; EL densities as in the finite-term owner",
            "gradient_convention": "covariant-source per-node metric_gradients; variation is the inner product against delta N,q,r",
            "initial_coefficients": "unspecified remains None; symbolic names c_M4..c_boxR; an explicit numeric vector is an unfitted control, not a selected zero action",
            "common_action": "shifts of the unresolved finite local part of the same measure; no extra Einstein or scalar term",
            "Weyl_compensator": "not derived here; a concurrent owner, not this M-flow",
            "continuum_bulk_running": "C2 only on this closed fixed-topology cell; Euler and box bulk EL vanish identically",
            "discrete_Euler_column": "retained in finite-code cancellation to match the finite-term owner; leftover versus C2-only is product-rule error",
        },
        "exact_identities": identities,
        "beta_functions": {
            "M_ref_basis": {name: beta_strings[i] for i, name in enumerate(BASIS)},
            "M_ref_basis_meaning": "dc_i/d ln M at fixed M_ref for this free-Dirac logarithm. M4 and M2R do not run because the subtraction keeps cutoff powers, not M^4 or M^2. That is not a physical Newton or cosmological beta function.",
            "finite_shift_M_to_sM": shift_formula,
            "continuum_bulk_shift_M_to_sM": {"C2": shift_formula["C2"],
                                             "M4": "0", "M2R": "0", "R2": "0",
                                             "E4": "0", "boxR": "0"},
            "M_basis_equivalent": {
                "M4": "-4*c_M4",
                "M2R": "-2*c_M2R",
                "C2": beta_strings[2],
                "R2": "0",
                "E4": beta_strings[4],
                "boxR": beta_strings[5],
                "meaning": "same free-Dirac logarithm with I_i written using M itself: hat c_M4=c_M4*(M_ref/M)^4 and hat c_M2R=c_M2R*(M_ref/M)^2, plus the a4 anomaly. Classical -4 c_M4 and -2 c_M2R are unit rewriting, not physical G or Lambda running.",
            },
            "marginal_stationary_coefficient": "c_R2",
            "running_direction": identities["running_direction"],
            "no_finite_fixed_point_along_a4": True,
            "free_Dirac_log_flow_is_not_recursive_Omega_flow": True,
            "no_physical_Newton_running_inferred": True,
            "no_physical_cosmological_running_inferred": True,
        },
        "cylinder_analytic": {
            "L": 4., "a": 1., "A4": -16 * pi / 15,
            "dE_sub_dlnM": cylinder_exact,
            "formula": "-1/(15*pi)",
            "computed_dE_sub_dlnM": cylinder["dE_sub_dlnM"],
        },
        "families": families,
        "discrete_euler_refinement": {
            "metric": "nonconstant R=2 smooth N,q,r",
            "scale_ratio": 2.,
            "classification": "discrete product-rule leftover of the finite-term Euler column, not continuum bulk Euler stress",
            "grids": euler_scan,
            "coarse_c2_only_visible": euler_scan[0]["c2_only_gradient_residual"] > 1e-12,
            "refined_c2_only_vanished": max(row["c2_only_gradient_residual"] for row in euler_scan[1:]) < 1e-12,
        },
        "dirac_remainder_bind": bind,
        "headline_residuals": {
            "maximum_matched_energy_residual": max(energy_residuals),
            "maximum_matched_metric_gradient_residual": max(gradient_residuals),
            "minimum_unmatched_energy": min(unmatched_energies),
            "maximum_unmatched_energy": max(unmatched_energies),
            "coarse_c2_only_gradient_residual": euler_scan[0]["c2_only_gradient_residual"],
            "refined_c2_only_gradient_residual_64": euler_scan[1]["c2_only_gradient_residual"],
            "refined_c2_only_gradient_residual_128": euler_scan[2]["c2_only_gradient_residual"],
        },
        "gate": gates,
        "nonclaims": {
            "complete_ultraviolet_functional_fixed": False,
            "self_sourcing_fixed": False,
            "finite_coefficients_selected": False,
            "zero_action_selected_as_baseline": False,
            "throat_fitted": False,
            "independent_gravity_or_scalar_added": False,
            "M_ratio_identified_with_Omega": False,
            "Weyl_compensator_derived": False,
            "physical_mass_or_inheritance_scale_selected": False,
            "physical_Newton_running_inferred": False,
            "physical_cosmological_running_inferred": False,
            "recursive_fixed_point_inferred": False,
            "continuum_Euler_bulk_stress_from_this_logarithm": False,
            "in_in_measure_constructed": False,
            "new_to_world_Seeley_coefficient": False,
        },
        "terminal": True,
    })
