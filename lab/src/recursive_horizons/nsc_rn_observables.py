"""Observers and the boundary ledger for the RN radial parent.

Vacuum geometry, horizons, SBP weights, the weighted Gram and curvature
invariants come from ``nsc_rn_reference``. This module does not rebuild
those formulas, does not form a second derivative by applying a first
derivative twice, and does not multiply the κ=1 shell into the Gram or
the angular CAR a second time.

Canonical φ has shape ``(2, N, rank)``. Occupations are the supplied
Gaussian weights. The Hamiltonian sample is the SBP quadrature of a
supplied density, or the positive exterior Killing weight of the same
quadrature. Multiplicity 4 is metadata and is applied once, outside
those products.
"""
from __future__ import annotations

import math

import numpy as np

from . import nsc_rn_reference as reference

API_VERSION = "nsc-rn-observables-v2"
CALIBRATION_TOLERANCE = 1e-4
SHELL_KAPPA = reference.KAPPA
SHELL_MULTIPLICITY = reference.MULTIPLICITY
OBSERVERS = (
    "child_near_horizon_shell",
    "parent_exterior",
    "outer_absorbing_boundary",
    "excision_outflow_boundary",
    "trapping_horizon",
    "outer_static_killing_clock",
    "interior_pg_normal_clock",
)
LEDGER_KEYS = (
    "probability_remaining",
    "energy_remaining",
    "outer_outflow",
    "excision_outflow",
    "sat_debit",
    "horizon_reynolds",
    "horizon_fixed_flux",
)
FORMULAS = {
    "chart": reference.ACTION["line_element"],
    "magnetic_radius": "r_m = sqrt(P^2), P^2 = magnetic_r2",
    "horizons": "roots of r^2 - 2 mass r + magnetic_r2 = 0, owned by RNReference.horizons",
    "excision": "(r_minus + r_plus) / 2",
    "misner_sharp": "nsc_rn_reference.misner_sharp",
    "charged_mass": "nsc_rn_reference.charged_mass_from_gradient",
    "curvature": "nsc_rn_reference.RNReference.invariants",
    "gram": "nsc_rn_reference.weighted_gram",
    "angular_coupling": "nsc_rn_reference.angular_coupling, once in the CAR",
    "spacing": "min(lambda_min/24, (r_plus-r_minus)/16, sigma/24)",
    "full_cfl": "dt = spacing / max(|N| + |beta|)",
    "endpoint_jet": "a D(2,1) product D@D is not a metric jet; endpoint f'' is half the consistent jet",
}


class ObservableError(ValueError):
    """The supplied sample cannot be observed under the declared contracts."""


def _finite_scalar(value, name):
    number = float(value)
    if not math.isfinite(number):
        raise ObservableError(f"{name} must be finite")
    return number


def _as_float_array(values, name):
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size == 0 or not np.all(np.isfinite(array)):
        raise ObservableError(f"{name} must be a nonempty finite radial array")
    return array


def require_phi(phi, points=None):
    """Canonical field with shape (2, N, rank)."""
    array = np.asarray(phi)
    if array.ndim != 3 or array.shape[0] != 2 or array.shape[1] == 0 or array.shape[2] == 0:
        raise ObservableError("canonical phi shape is (2, N, rank)")
    if points is not None and array.shape[1] != int(points):
        raise ObservableError("canonical phi radial length does not match the grid")
    if not np.all(np.isfinite(array.real)) or not np.all(np.isfinite(array.imag)):
        raise ObservableError("canonical phi must be finite")
    return array


def require_weights(weights, points):
    array = _as_float_array(weights, "SBP weights")
    if array.shape != (int(points),):
        raise ObservableError("SBP weights do not match radial points")
    if np.any(array <= 0.0):
        raise ObservableError("SBP weights must be positive")
    return array


def chart_speeds(lapse, shift):
    """Radial null speeds for a supplied lapse. Vacuum N=1 agrees with the reference."""
    lapse = _as_float_array(np.atleast_1d(lapse), "lapse")
    shift = _as_float_array(np.atleast_1d(shift), "shift")
    if lapse.shape != shift.shape or np.any(lapse <= 0.0):
        raise ObservableError("speeds need a positive lapse and a matching shift")
    return -shift + lapse, -shift - lapse


characteristic_speeds = chart_speeds


def reference_model(mass, charge):
    """Reference RN with magnetic radius square charge². charge is Q, not r₋."""
    charge = _finite_scalar(charge, "charge")
    mass = _finite_scalar(mass, "mass")
    if charge == 0.0 or mass <= 0.0:
        raise ObservableError("magnetic RN needs a positive mass and a nonzero charge")
    target = charge * charge
    audited = reference.magnetic_radius_square(
        reference.AUDITED_A, reference.AUDITED_C_F, reference.AUDITED_FLUX,
    )
    flux = 1.0
    coupling = reference.AUDITED_C_F * (target / audited) * (reference.AUDITED_FLUX / flux) ** 2
    model = reference.RNReference.from_action(reference.AUDITED_A, coupling, flux, mass)
    if abs(model.magnetic_r2 - target) > 1e-8 * max(1.0, target):
        raise ObservableError("reference magnetic radius did not realize P^2 = Q^2")
    return model


def rn_horizons(mass, charge):
    """Reference roots. The charge argument is Q, and r_m is not r₋."""
    model = reference_model(mass, charge)
    horizons = model.horizons
    if not horizons["black_hole"] or horizons["r_minus"] is None:
        raise ObservableError("reference RN has no horizon pair")
    return {
        "mass": float(model.mass),
        "charge": float(charge),
        "magnetic_r2": float(model.magnetic_r2),
        "r_minus": float(horizons["r_minus"]),
        "r_plus": float(horizons["r_plus"]),
        "horizon_gap": float(horizons["r_plus"] - horizons["r_minus"]),
        "r_m_is_inner_horizon": False,
    }


def declare_cuts(*, mass, charge, r_out=None):
    """Excision at the horizon midpoint. r_m = |Q| is not copied into r₋."""
    horizons = rn_horizons(mass, charge)
    model = reference_model(mass, charge)
    r_minus = horizons["r_minus"]
    r_plus = horizons["r_plus"]
    gap = horizons["horizon_gap"]
    magnetic_radius = abs(float(charge))
    outer = float(32.0 * magnetic_radius if r_out is None else r_out)
    excision = 0.5 * (r_minus + r_plus)
    child = (r_plus + gap / 16.0, r_plus + 0.5 * gap)
    parent = (child[1], outer)
    placed = reference.placement(
        model, parent=parent, child=child, excision=excision, outer="absorbing",
    )
    if not placed["supported"]:
        raise ObservableError("reference placement refused the fixed cuts: " + "; ".join(placed["reasons"]))
    return {
        "observers": OBSERVERS,
        "observers_declared_once": True,
        "cuts_declared_once": True,
        "same_areal_radius": False,
        "r_m_identified_with_r_minus": False,
        "r_minus": r_minus,
        "r_plus": r_plus,
        "excision": excision,
        "child": child,
        "parent": parent,
        "outer": outer,
        "child_interval": child,
        "parent_interval": parent,
        "excision_between_horizons": True,
        "pure_outflow_excision": True,
        "outer_boundary": "nonperiodic_absorbing",
        "periodic": False,
        "placement_supported": True,
    }


def weighted_gram(phi, weights):
    """Reference Gram. Occupations are its diagonal and are not renormalized."""
    array = require_phi(phi)
    quadrature = require_weights(weights, array.shape[1])
    gram = reference.weighted_gram(array, quadrature)
    occupations = np.real(np.diag(gram))
    return {
        "gram": gram,
        "occupations": occupations,
        "trace": float(np.real(np.trace(gram))),
        "normalized": False,
        "normalized_survivors": False,
        "multiplicity_in_gram": False,
    }


def shell_factor(rank, multiplicity):
    """One shell factor. Rank 4 together with multiplicity 4 is a second insertion."""
    rank = int(rank)
    multiplicity = int(multiplicity)
    if rank == 1 and multiplicity == SHELL_MULTIPLICITY:
        return float(SHELL_MULTIPLICITY), 1
    if rank == SHELL_MULTIPLICITY and multiplicity == 1:
        return 1.0, 1
    if rank == SHELL_MULTIPLICITY and multiplicity == SHELL_MULTIPLICITY:
        raise ObservableError("shell multiplicity would be counted twice")
    raise ObservableError("rank and shell multiplicity do not form one closed shell")


def misner_sharp(radius, lapse, shift):
    """Bare Misner–Sharp mass through the reference gradient map."""
    return reference.misner_sharp(radius, grad_r_squared(lapse, shift))


def vacuum_pg_rn(radius, mass, charge):
    """Reference vacuum jets. This does not replace a coupled lapse."""
    model = reference_model(mass, charge)
    jets = model.pg_metric_jets(np.asarray(radius, dtype=float))
    return {
        "radius": np.asarray(radius, dtype=float),
        "lapse": np.asarray(jets["N"], dtype=float),
        "shift": np.asarray(jets["beta"], dtype=float),
        "mass": float(model.mass),
        "charge": float(charge),
        "magnetic_r2": float(model.magnetic_r2),
        "coupled_lapse_replaced": False,
        "per_radius_refit": False,
        "vacuum_pg_lapse_identity": True,
    }


def static_killing_redshift(metric_tt):
    """Static clock only where ∂ₜ is timelike."""
    metric_tt = _finite_scalar(metric_tt, "g_tt")
    if metric_tt <= 0.0:
        raise ObservableError("static Killing clock is refused outside the exterior timelike region")
    return math.sqrt(metric_tt)


def pg_normal_redshift(lapse):
    """Normal-observer proper time per coordinate time, including infall."""
    lapse = _finite_scalar(lapse, "lapse")
    if lapse <= 0.0:
        raise ObservableError("PG normal clock needs a positive lapse")
    return lapse


def grad_r_squared(lapse, shift):
    """g^{rr} on this chart, the argument of the reference Misner–Sharp mass."""
    lapse = _as_float_array(np.atleast_1d(lapse), "lapse")
    shift = _as_float_array(np.atleast_1d(shift), "shift")
    if np.any(lapse <= 0.0):
        raise ObservableError("lapse must be positive")
    return (shift / lapse) ** 2 - 1.0


def charged_mass(radius, lapse, shift, magnetic_r2):
    """One enclosed charge through the reference mass map. Not a radial refit."""
    if np.ndim(magnetic_r2) != 0:
        raise ObservableError("per-radius charge refit refused")
    return reference.charged_mass_from_gradient(
        radius, grad_r_squared(lapse, shift), magnetic_r2,
    )


def curvature_from_reference(model, radius):
    """R₄ and the tidal row from the reference reduction, on every supplied radius."""
    radius = _as_float_array(np.atleast_1d(radius), "radius")
    invariants = model.invariants(radius)
    tides = model.tetrad_riemann(radius)
    scalar = np.asarray(invariants["R4"], dtype=float)
    return {
        "R4": scalar,
        "scalar_curvature": scalar,
        "Ricci2": np.asarray(invariants["Ricci2"], dtype=float),
        "K": np.asarray(invariants["K"], dtype=float),
        "R_0101": np.asarray(tides["R_0101"], dtype=float),
        "R_0202": np.asarray(tides["R_0202"], dtype=float),
        "rows": int(radius.size),
        "all_rows": True,
        "jet_provenance": "nsc_rn_reference.invariants",
        "jets_differentiated_here": False,
        "naive_D_at_D": False,
    }


def _pg_christoffel(radius, theta, lapse, shift, lapse_r, shift_r, lapse_t=None, shift_t=None):
    """Γ^λ_μν of the ingoing PG chart at one polar angle.

    Time derivatives of the metric are 2 N N_t - 2 β β_t and -β_t. Second
    time derivatives are not supplied and are not invented.
    """
    radius = np.asarray(radius, dtype=float)
    lapse = np.asarray(lapse, dtype=float)
    shift = np.asarray(shift, dtype=float)
    lapse_r = np.asarray(lapse_r, dtype=float)
    shift_r = np.asarray(shift_r, dtype=float)
    count = radius.size
    if lapse_t is None:
        lapse_t = np.zeros(count)
    if shift_t is None:
        shift_t = np.zeros(count)
    lapse_t = np.asarray(lapse_t, dtype=float)
    shift_t = np.asarray(shift_t, dtype=float)
    sine = math.sin(theta)
    cosine = math.cos(theta)
    inverse = np.zeros((count, 4, 4))
    lapse2 = lapse * lapse
    inverse[:, 0, 0] = 1.0 / lapse2
    inverse[:, 0, 1] = inverse[:, 1, 0] = -shift / lapse2
    inverse[:, 1, 1] = shift * shift / lapse2 - 1.0
    inverse[:, 2, 2] = -1.0 / radius ** 2
    inverse[:, 3, 3] = -1.0 / (radius ** 2 * sine ** 2)
    derivative = np.zeros((count, 4, 4, 4))
    derivative[:, 0, 0, 0] = 2.0 * lapse * lapse_t - 2.0 * shift * shift_t
    derivative[:, 0, 0, 1] = derivative[:, 0, 1, 0] = -shift_t
    derivative[:, 1, 0, 0] = 2.0 * lapse * lapse_r - 2.0 * shift * shift_r
    derivative[:, 1, 0, 1] = derivative[:, 1, 1, 0] = -shift_r
    derivative[:, 1, 2, 2] = -2.0 * radius
    derivative[:, 1, 3, 3] = -2.0 * radius * sine ** 2
    derivative[:, 2, 3, 3] = -2.0 * radius ** 2 * sine * cosine
    gamma = np.zeros((count, 4, 4, 4))
    for lam in range(4):
        for mu in range(4):
            for nu in range(mu, 4):
                accumulator = np.zeros(count)
                for sig in range(4):
                    accumulator = accumulator + inverse[:, lam, sig] * (
                        derivative[:, mu, nu, sig] + derivative[:, nu, mu, sig]
                        - derivative[:, sig, mu, nu]
                    )
                gamma[:, lam, mu, nu] = 0.5 * accumulator
                gamma[:, lam, nu, mu] = gamma[:, lam, mu, nu]
    return gamma


def curvature_from_pg_jets(radius, lapse, shift, *, lapse_r, shift_r, lapse_rr, shift_rr,
                           lapse_t=None, shift_t=None, shift_tr=None):
    """R₄ and Ricci² from the PG Riemann contraction.

    Stationary slices use the radial jets. A moving slice also uses the
    constraint rates N_t and β_t and the radial derivative β_tr of β_t.
    N_tr and the second time derivatives were checked against this same
    contraction and do not enter R₄ or Ricci², so they are not invented.
    The scalar is not replaced by a trace-free condition.
    """
    radius = _as_float_array(np.atleast_1d(radius), "radius")
    lapse = _as_float_array(np.atleast_1d(lapse), "lapse")
    shift = _as_float_array(np.atleast_1d(shift), "shift")
    lapse_r = _as_float_array(np.atleast_1d(lapse_r), "lapse_r")
    shift_r = _as_float_array(np.atleast_1d(shift_r), "shift_r")
    lapse_rr = _as_float_array(np.atleast_1d(lapse_rr), "lapse_rr")
    shift_rr = _as_float_array(np.atleast_1d(shift_rr), "shift_rr")
    fields = (radius, lapse, shift, lapse_r, shift_r, lapse_rr, shift_rr)
    if len({field.shape for field in fields}) != 1:
        raise ObservableError("curvature jets do not share the radial grid")
    if np.any(lapse <= 0.0) or np.any(radius <= 0.0):
        raise ObservableError("curvature needs a positive lapse and radius")
    time_scale = 0.0
    if lapse_t is not None or shift_t is not None:
        if lapse_t is None or shift_t is None:
            raise ObservableError("lapse and shift time jets must be supplied together")
        lapse_t = _as_float_array(np.atleast_1d(lapse_t), "lapse_t")
        shift_t = _as_float_array(np.atleast_1d(shift_t), "shift_t")
        if lapse_t.shape != radius.shape or shift_t.shape != radius.shape:
            raise ObservableError("time jets do not share the radial grid")
        time_scale = float(max(np.max(np.abs(lapse_t)), np.max(np.abs(shift_t))))
        if shift_tr is None:
            if np.max(np.abs(shift_t)) <= 1e-14:
                shift_tr = np.zeros_like(shift_t)
            else:
                raise ObservableError("shift_tr is the radial derivative of shift_t")
        else:
            shift_tr = _as_float_array(np.atleast_1d(shift_tr), "shift_tr")
            if shift_tr.shape != radius.shape:
                raise ObservableError("shift_tr does not share the radial grid")
    else:
        lapse_t = np.zeros_like(radius)
        shift_t = np.zeros_like(radius)
        shift_tr = np.zeros_like(radius)
    moving = time_scale > 1e-8 or float(np.max(np.abs(shift_tr))) > 1e-8
    theta = 0.5 * math.pi
    step = 1.0e-6

    def christoffel_at(radius_value, lapse_value, lapse_slope, shift_value, shift_slope,
                       lapse_time, shift_time):
        return _pg_christoffel(
            radius_value, theta, lapse_value, shift_value, lapse_slope, shift_slope,
            lapse_time, shift_time,
        )

    base = christoffel_at(radius, lapse, lapse_r, shift, shift_r, lapse_t, shift_t)
    radial = (
        christoffel_at(
            radius + step,
            lapse + step * lapse_r + 0.5 * step ** 2 * lapse_rr, lapse_r + step * lapse_rr,
            shift + step * shift_r + 0.5 * step ** 2 * shift_rr, shift_r + step * shift_rr,
            lapse_t, shift_t + step * shift_tr,
        )
        - christoffel_at(
            radius - step,
            lapse - step * lapse_r + 0.5 * step ** 2 * lapse_rr, lapse_r - step * lapse_rr,
            shift - step * shift_r + 0.5 * step ** 2 * shift_rr, shift_r - step * shift_rr,
            lapse_t, shift_t - step * shift_tr,
        )
    ) / (2.0 * step)
    temporal = (
        christoffel_at(
            radius, lapse + step * lapse_t, lapse_r,
            shift + step * shift_t, shift_r + step * shift_tr, lapse_t, shift_t,
        )
        - christoffel_at(
            radius, lapse - step * lapse_t, lapse_r,
            shift - step * shift_t, shift_r - step * shift_tr, lapse_t, shift_t,
        )
    ) / (2.0 * step)
    polar = (
        _pg_christoffel(
            radius, theta + step, lapse, shift, lapse_r, shift_r, lapse_t, shift_t,
        )
        - _pg_christoffel(
            radius, theta - step, lapse, shift, lapse_r, shift_r, lapse_t, shift_t,
        )
    ) / (2.0 * step)
    partial = (temporal, radial, polar, None)
    count = radius.size
    riemann = np.zeros((count, 4, 4, 4, 4))
    for first in range(4):
        for second in range(4):
            for third in range(4):
                for fourth in range(4):
                    quadratic = np.zeros(count)
                    for mid in range(4):
                        quadratic += (
                            base[:, first, third, mid] * base[:, mid, fourth, second]
                            - base[:, first, fourth, mid] * base[:, mid, third, second]
                        )
                    left = 0.0 if partial[third] is None else partial[third][:, first, fourth, second]
                    right = 0.0 if partial[fourth] is None else partial[fourth][:, first, third, second]
                    riemann[:, first, second, third, fourth] = left - right + quadratic
    inverse = np.zeros((count, 4, 4))
    inverse[:, 0, 0] = 1.0 / lapse ** 2
    inverse[:, 0, 1] = inverse[:, 1, 0] = -shift / lapse ** 2
    inverse[:, 1, 1] = shift ** 2 / lapse ** 2 - 1.0
    inverse[:, 2, 2] = -1.0 / radius ** 2
    inverse[:, 3, 3] = -1.0 / radius ** 2
    ricci = np.zeros((count, 4, 4))
    for second in range(4):
        for fourth in range(4):
            for slot in range(4):
                ricci[:, second, fourth] += riemann[:, slot, second, slot, fourth]
    scalar = np.zeros(count)
    ricci2 = np.zeros(count)
    for second in range(4):
        for fourth in range(4):
            scalar += inverse[:, second, fourth] * ricci[:, second, fourth]
    raised = np.zeros((count, 4, 4))
    for second in range(4):
        for fourth in range(4):
            for left in range(4):
                for right in range(4):
                    raised[:, second, fourth] += (
                        inverse[:, second, left] * inverse[:, fourth, right] * ricci[:, left, right]
                    )
    for second in range(4):
        for fourth in range(4):
            ricci2 += ricci[:, second, fourth] * raised[:, second, fourth]
    return {
        "status": "sourced_time_dependent_jets" if moving else "stationary_radial_jets",
        "R4": scalar,
        "Ricci2": ricci2,
        "time_rate_max": time_scale,
        "shift_tr_max": float(np.max(np.abs(shift_tr))),
        "outer_lapse_rate": float(lapse_t[-1]),
        "static_formula_applied": False,
        "trace_free_condition_imposed": False,
        "time_rates_negligible": not moving,
        "jets_differentiated_here": True,
        "naive_D_at_D": False,
        "time_jet_provenance": (
            "constraint_N_t_beta_t_and_nodal_beta_tr"
            if moving else "time_rates_below_1e-8"
        ),
        "lapse_tr_used": False,
        "second_time_derivatives_used": False,
        "suggested_box_formula_imposed": False,
    }


def endpoint_second_derivative_bias(supplied, consistent):
    """Detect a naive D@D endpoint jet. Do not replace it with twice that value."""
    supplied = _as_float_array(supplied, "supplied second jet")
    consistent = _as_float_array(consistent, "consistent second jet")
    if supplied.shape != consistent.shape:
        raise ObservableError("second-jet comparison lengths differ")
    ratios = []
    half = False
    for index in (0, -1):
        actual = float(consistent[index])
        sample = float(supplied[index])
        if abs(actual) <= 1e-8:
            ratios.append(None)
            continue
        ratio = sample / actual
        ratios.append(ratio)
        if abs(ratio - 0.5) <= 0.05 and abs(ratio - 1.0) > 0.2:
            half = True
    return {
        "endpoint_half_actual": half,
        "endpoint_ratios": ratios,
        "artifact": "D_at_D_half_actual" if half else None,
        "repaired": False,
        "used_as_metric_jet": False,
    }


def trapping_horizon(radius, lapse, shift):
    """Outermost simple zero of the outgoing speed on supplied samples."""
    radius = _as_float_array(radius, "radius")
    outgoing, ingoing = chart_speeds(lapse, shift)
    if np.any(np.diff(radius) <= 0.0):
        raise ObservableError("radius must increase")
    crossing = None
    for index in range(1, radius.size):
        left, right = float(outgoing[index - 1]), float(outgoing[index])
        if left == 0.0:
            crossing = float(radius[index - 1])
        elif left < 0.0 <= right:
            span = right - left
            fraction = 0.0 if span == 0.0 else -left / span
            crossing = float(radius[index - 1] + fraction * (radius[index] - radius[index - 1]))
    return {
        "trapping_horizon": crossing,
        "outgoing_speed": outgoing,
        "ingoing_speed": ingoing,
        "located": crossing is not None,
        "event_horizon_copied_from_trapping": False,
    }


def horizon_report(model, radius, lapse, shift, *, stationary):
    """Reference r₊, the located trapping radius, and the event horizon stay distinct."""
    parameter = model.horizons
    located = trapping_horizon(radius, lapse, shift)
    trapping = located["trapping_horizon"]
    if stationary and trapping is not None:
        event = trapping
        status = "stationary_slice_event_horizon_coincides_with_trapping"
        coincides = True
    else:
        event = None
        status = "unresolved_on_dynamical_slice"
        coincides = False
    return {
        "trapping_horizon": trapping,
        "rn_r_plus": parameter["r_plus"],
        "rn_r_minus": parameter["r_minus"],
        "event_horizon": event,
        "event_horizon_status": status,
        "event_horizon_equals_trapping": coincides,
        "stationary_slice": bool(stationary),
        "r_m_identified_with_r_minus": False,
    }


def shell_accounts(phi, weights, radius, *, killing_frequency, hamiltonian_density=None,
                   metric_tt=None):
    """H quadrature and one CAR insertion. Multiplicity is not applied again."""
    array = require_phi(phi)
    points = array.shape[1]
    weights = require_weights(weights, points)
    radius = _as_float_array(radius, "radius")
    if radius.shape != (points,):
        raise ObservableError("radius does not match phi")
    frequency = _finite_scalar(killing_frequency, "Killing frequency")
    if frequency <= 0.0:
        raise ObservableError("exterior Killing frequency must be positive")
    gram = reference.weighted_gram(array, weights)
    density = np.sum(np.abs(array) ** 2, axis=(0, 2))
    probability = float(np.real(np.trace(gram)))
    coupling = np.asarray(reference.angular_coupling(radius), dtype=float)
    car = float(np.dot(weights, coupling * density))
    if metric_tt is None:
        timelike = np.ones(points, dtype=bool)
    else:
        metric_tt = _as_float_array(metric_tt, "g_tt")
        if metric_tt.shape != (points,):
            raise ObservableError("g_tt does not match phi")
        timelike = metric_tt > 0.0
    if hamiltonian_density is None:
        energy = frequency * float(np.dot(weights, density * timelike))
        energy_source = "exterior Killing weight on the SBP quadrature"
    else:
        supplied = _as_float_array(hamiltonian_density, "Hamiltonian density")
        if supplied.shape != (points,):
            raise ObservableError("Hamiltonian density does not match phi")
        energy = float(np.dot(weights, supplied))
        energy_source = "supplied H quadrature"
    return {
        "gram": gram,
        "occupations_normalized_to_survivors": False,
        "probability_remaining": probability,
        "energy_remaining": energy,
        "energy_source": energy_source,
        "car_angular": car,
        "multiplicity": SHELL_MULTIPLICITY,
        "multiplicity_applications_in_car": 1,
        "multiplicity_in_gram": False,
        "multiplicity_in_H": False,
        "kappa": SHELL_KAPPA,
        "angular_is_4d_mass": False,
        "filled_sea": False,
        "killing_frequency": frequency,
        "V4": 0.0,
    }


def _interp(radius, values, at):
    radius = _as_float_array(radius, "radius")
    values = _as_float_array(values, "sample")
    at = _finite_scalar(at, "sample radius")
    if at < float(radius[0]) or at > float(radius[-1]):
        raise ObservableError("horizon sample lies outside the radial nodes")
    index = int(np.searchsorted(radius, at))
    if index <= 0:
        return float(values[0])
    if index >= radius.size or float(radius[index - 1]) == at:
        return float(values[min(index, radius.size - 1)])
    left, right = float(radius[index - 1]), float(radius[index])
    fraction = (at - left) / (right - left)
    return float((1.0 - fraction) * values[index - 1] + fraction * values[index])


def coordinate_current(phi, lapse, shift):
    array = require_phi(phi)
    outgoing_speed, ingoing_speed = chart_speeds(lapse, shift)
    if outgoing_speed.shape[0] != array.shape[1]:
        raise ObservableError("characteristic speeds do not match phi")
    rho_out = np.sum(np.abs(array[0]) ** 2, axis=1)
    rho_in = np.sum(np.abs(array[1]) ** 2, axis=1)
    return outgoing_speed * rho_out + ingoing_speed * rho_in, rho_out + rho_in


def horizon_transport(density, flux, radius_rate):
    """Reynolds transport beside the fixed characteristic flux."""
    density = _finite_scalar(density, "horizon density")
    flux = _finite_scalar(flux, "horizon flux")
    radius_rate = _finite_scalar(radius_rate, "horizon radius rate")
    reynolds = density * radius_rate
    return {
        "horizon_reynolds": reynolds,
        "horizon_fixed_flux": flux,
        "horizon_flux": reynolds + flux,
        "piston_work_included": False,
    }


def boundary_ledger(*, phi, weights, radius, lapse, shift, killing_frequency,
                    metric_tt=None, hamiltonian_density=None, sat_debit=None,
                    sat_enabled=False, horizon_radius=None, horizon_radius_rate=0.0,
                    excision_index=0, outer_index=-1, electric_current=0.0,
                    magnetic_flux_fixed=True):
    """Distinct remaining, outflow and penalty accounts."""
    accounts = shell_accounts(
        phi, weights, radius, killing_frequency=killing_frequency,
        hamiltonian_density=hamiltonian_density, metric_tt=metric_tt,
    )
    array = require_phi(phi)
    points = array.shape[1]
    lapse = _as_float_array(lapse, "lapse")
    shift = _as_float_array(shift, "shift")
    radius = _as_float_array(radius, "radius")
    if not (radius.shape == lapse.shape == shift.shape == (points,)):
        raise ObservableError("ledger geometry does not match phi")
    if float(electric_current) != 0.0:
        raise ObservableError("electric current is not part of the closed neutral shell")
    if magnetic_flux_fixed is not True:
        raise ObservableError("magnetic flux is fixed")
    current, density = coordinate_current(array, lapse, shift)
    outer = int(outer_index)
    excision = int(excision_index)
    if outer < 0:
        outer += points
    if not (0 <= excision < points and 0 <= outer < points and excision != outer):
        raise ObservableError("excision and outer ledger nodes must be distinct")
    outgoing, ingoing = chart_speeds(lapse, shift)
    if sat_enabled:
        if sat_debit is None:
            raise ObservableError("SAT debit must be supplied when SAT is enabled")
        penalty = _finite_scalar(sat_debit, "SAT debit")
    else:
        if sat_debit not in (None, 0, 0.0):
            raise ObservableError("SAT debit was supplied while SAT is disabled")
        penalty = 0.0
    if horizon_radius is None:
        transport = horizon_transport(0.0, 0.0, 0.0)
    else:
        transport = horizon_transport(
            _interp(radius, density, horizon_radius),
            _interp(radius, current, horizon_radius),
            horizon_radius_rate,
        )
    ledger = {
        "probability_remaining": accounts["probability_remaining"],
        "energy_remaining": accounts["energy_remaining"],
        "outer_outflow": float(current[outer]),
        "excision_outflow": float(-current[excision]),
        "sat_debit": penalty,
        "horizon_reynolds": transport["horizon_reynolds"],
        "horizon_fixed_flux": transport["horizon_fixed_flux"],
        "horizon_flux": transport["horizon_flux"],
        "car_angular": accounts["car_angular"],
        "energy_source": accounts["energy_source"],
        "multiplicity": accounts["multiplicity"],
        "multiplicity_applications_in_car": 1,
        "multiplicity_in_gram": False,
        "multiplicity_in_H": False,
        "occupations_are_normalized_survivors": False,
        "pure_outflow_excision": bool(outgoing[excision] < 0.0 and ingoing[excision] < 0.0),
        "sat_enabled": bool(sat_enabled),
        "sat_added_into_outflow": False,
        "electric_current": 0.0,
        "magnetic_flux_fixed": True,
        "V4": 0.0,
        "filled_sea": False,
        "accounts_distinct": True,
    }
    return ledger


def compare_station(sample, initial, enclosed_mass):
    """One initial RN and one enclosed-mass RN. The sample lapse is not replaced."""
    if np.ndim(enclosed_mass) != 0 or np.ndim(initial.mass) != 0:
        raise ObservableError("per-radius refit refused")
    radius = _finite_scalar(sample["radius"], "station radius")
    lapse = _finite_scalar(sample["lapse"], "station lapse")
    shift = _finite_scalar(sample["shift"], "station shift")
    if lapse <= 0.0 or radius <= 0.0:
        raise ObservableError("station sample needs positive radius and lapse")
    matched = reference.RNReference.from_action(
        initial.A, initial.C_F, initial.flux, float(enclosed_mass),
    )
    initial_jets = initial.pg_metric_jets(radius)
    matched_jets = matched.pg_metric_jets(radius)
    initial_lapse = float(np.asarray(initial_jets["N"]).reshape(-1)[0])
    initial_shift = float(np.asarray(initial_jets["beta"]).reshape(-1)[0])
    matched_lapse = float(np.asarray(matched_jets["N"]).reshape(-1)[0])
    matched_shift = float(np.asarray(matched_jets["beta"]).reshape(-1)[0])
    return {
        "radius": radius,
        "sample_lapse": lapse,
        "sample_shift": shift,
        "coupled_lapse_replaced": False,
        "per_radius_refit": False,
        "initial_rn_lapse": initial_lapse,
        "initial_rn_shift": initial_shift,
        "matched_rn_mass": float(enclosed_mass),
        "matched_rn_lapse": matched_lapse,
        "matched_rn_shift": matched_shift,
        "magnetic_r2": float(initial.magnetic_r2),
        "lapse_residual_against_initial": lapse - initial_lapse,
        "shift_residual_against_initial": shift - initial_shift,
        "lapse_residual_against_matched": lapse - matched_lapse,
        "shift_residual_against_matched": shift - matched_shift,
    }


def target_spacing(lambda_min, horizon_gap, sigma):
    """Physical scales. The throat piece is (r₊ − r₋)/16."""
    parts = (
        _finite_scalar(lambda_min, "lambda_min") / 24.0,
        _finite_scalar(horizon_gap, "horizon_gap") / 16.0,
        _finite_scalar(sigma, "sigma") / 24.0,
    )
    if min(parts) <= 0.0:
        raise ObservableError("grid scales must be positive")
    return min(parts)


def throat_spacing(horizon_gap, lambda_min=None, sigma=None):
    """Delta-throat rule, tightened by λ_min and σ only when those scales exist."""
    gap = _finite_scalar(horizon_gap, "horizon_gap")
    if gap <= 0.0:
        raise ObservableError("horizon gap must be positive")
    parts = [gap / 16.0]
    pending = []
    if lambda_min is None:
        pending.append("lambda_min")
    else:
        parts.append(_finite_scalar(lambda_min, "lambda_min") / 24.0)
    if sigma is None:
        pending.append("sigma")
    else:
        parts.append(_finite_scalar(sigma, "sigma") / 24.0)
    if min(parts) <= 0.0:
        raise ObservableError("grid scales must be positive")
    return {"spacing": min(parts), "pending_scales": pending, "derivative_ladder": False}


def route_points(length, spacing, route):
    if int(route) not in (32, 64):
        raise ObservableError("route floor is 32 or the 64 control")
    length = _finite_scalar(length, "length")
    spacing = _finite_scalar(spacing, "spacing")
    if length <= 0.0 or spacing <= 0.0:
        raise ObservableError("length and spacing must be positive")
    return max(int(route), int(math.ceil(length / spacing)) + 1)


def full_cfl_dt(lapse, shift, spacing):
    """Model CFL from the full characteristic speed |N| + |β|."""
    lapse = np.abs(_as_float_array(np.atleast_1d(lapse), "lapse"))
    shift = np.abs(_as_float_array(np.atleast_1d(shift), "shift"))
    spacing = _finite_scalar(spacing, "spacing")
    speed = float(np.max(lapse + shift))
    if spacing <= 0.0 or speed <= 0.0:
        raise ObservableError("full CFL needs a positive spacing and speed")
    return spacing / speed
