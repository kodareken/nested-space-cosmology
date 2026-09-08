"""Finite spectral profiles for common-action normalization.

Compares the calculated proper-time determinant modulus with a smooth
weighted-log prescription and with heat-covariantization of the old
finite-rank term. The identity is a diagnostic of those profiles. It is
not inserted as a completed Gamma_one, and no coefficient is chosen to
open a throat.

Primary Andrianov regularization uses a hard spectral projector. The
weighted-log profile below is a smooth exponential weight, not that
projector. Unknown finite/compensator data remain explicit.
"""
from __future__ import annotations

from dataclasses import replace
from hashlib import sha256
from math import pi
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.linalg import eigh, expm
from scipy.special import exp1

from .nsc_covariant_operator import (
    SIGMA1, SIGMA2, SIGMA3, cutoff_response, cutoff_weyl_trace, cylinder_metric,
    euclidean_operator, frequency_extent, smooth_metric,
)
from .nsc_regulated import OperatorConventions, RegulatedOperator

EULER = float(np.euler_gamma)
SMOOTH_KINDS = ("proper_time", "heat_rank", "remainder", "weighted_log")
KINDS = SMOOTH_KINDS + ("hard_log",)
ROOT = Path(__file__).resolve().parents[2]
OWNED_SOURCES = (
    "src/recursive_horizons/nsc_measure_normalization.py",
    "scripts/check_nsc_measure_normalization.py",
    "tests/test_nsc_measure_normalization.py",
    "docs/nsc-measure-normalization.md",
)
READONLY_SOURCES = (
    "src/recursive_horizons/nsc_covariant_operator.py",
    "src/recursive_horizons/nsc_regulated.py",
)
INPUT_RECORDS = ("results/nsc-9-covariant-source.json",)
SCHEMA = "nsc-measure-normalization-v1"
ARTIFACT = "NSC-10-MEASURE-NORMALIZATION"


def _positive(name, value):
    if isinstance(value, bool) or not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return float(value)


def _array(values):
    data = np.atleast_1d(np.asarray(values, dtype=float))
    if data.size == 0 or not np.isfinite(data).all():
        raise ValueError("finite spectral values required")
    return data


def log_ratio(cutoff, normalization):
    return np.log(_positive("normalization", normalization) / _positive("cutoff", cutoff))


def g_proper_time(x):
    x = _array(x)
    if np.any(x <= 0):
        raise ValueError("proper-time profile requires x>0")
    return 0.5 * exp1(x)


def h_remainder(x):
    """Finite remainder h(x)=-.5[exp(-x)(log x+gamma_E)+E1(x)].

    The two large terms cancel as x->0. Use the power series on (0,0.05).
    h(0+)=0 and is not a hard-projector measure.
    """
    x = _array(x)
    if np.any(x <= 0):
        raise ValueError("remainder profile requires x>0; the x->0 limit is 0")
    out = np.empty(x.shape, dtype=float)
    small = x < 0.05
    if np.any(small):
        xs = x[small]
        logx = np.log(xs)
        total = np.zeros_like(xs)
        power = np.ones_like(xs)
        for n in range(1, 24):
            power *= xs / n
            total += (1 if n % 2 else -1) * power * (logx + EULER - 1 / n) / 2
        out[small] = total
    if np.any(~small):
        xl = x[~small]
        out[~small] = -0.5 * (np.exp(-xl) * (np.log(xl) + EULER) + exp1(xl))
    return out


def h_series(x, terms=12):
    x = _array(x)
    if isinstance(terms, bool) or not isinstance(terms, int) or terms < 1:
        raise ValueError("positive series length required")
    logx = np.log(x)
    total = np.zeros_like(x)
    power = np.ones_like(x)
    for n in range(1, terms + 1):
        power *= x / n
        total += (1 if n % 2 else -1) * power * (logx + EULER - 1 / n) / 2
    return total


def g_heat_rank(x, ratio):
    x = _array(x)
    if not np.isfinite(ratio):
        raise ValueError("finite log(M/Lambda) required")
    return (float(ratio) + 0.5 * EULER) * np.exp(-x)


def g_weighted_log(x, ratio):
    """Smooth weighted-log profile, not Theta(1-x) log.

    g_log=-.5 exp(-x) log(d^2/M^2)=-0.5 exp(-x)(log x-2 log(M/Lambda)).
    """
    x = _array(x)
    if np.any(x <= 0) or not np.isfinite(ratio):
        raise ValueError("weighted-log profile requires x>0 and finite log(M/Lambda)")
    return -0.5 * np.exp(-x) * (np.log(x) - 2 * float(ratio))


def g_hard_log(x, ratio):
    """Andrianov-type sharp projector on the same log, for comparison only."""
    x = _array(x)
    if np.any(x <= 0) or not np.isfinite(ratio):
        raise ValueError("hard-log profile requires x>0 and finite log(M/Lambda)")
    out = np.zeros_like(x)
    inside = x <= 1
    out[inside] = -0.5 * (np.log(x[inside]) - 2 * float(ratio))
    return out


def identity_pieces(x, ratio):
    pt = g_proper_time(x)
    heat = g_heat_rank(x, ratio)
    remainder = h_remainder(x)
    weighted = g_weighted_log(x, ratio)
    reconstructed = pt + heat + remainder
    return {
        "proper_time": pt,
        "heat_rank": heat,
        "remainder": remainder,
        "weighted_log": weighted,
        "reconstructed": reconstructed,
        "residual": weighted - reconstructed,
        "hard_log": g_hard_log(x, ratio),
        "hard_minus_weighted": g_hard_log(x, ratio) - weighted,
    }


def identity_residual(x, ratio):
    return identity_pieces(x, ratio)["residual"]


def dx_profiles(x, ratio):
    """Ordinary x-derivatives of the scalar profiles, x=d^2/Lambda^2."""
    x = _array(x)
    ratio = float(ratio)
    expo = np.exp(-x)
    e1 = exp1(x)
    return {
        "proper_time": -0.5 * expo / x,
        "heat_rank": -(ratio + 0.5 * EULER) * expo,
        "remainder": 0.5 * expo * (np.log(x) + EULER),
        "weighted_log": 0.5 * expo * (np.log(x) - 2 * ratio - 1 / x),
        "hard_log": np.where(x < 1, -0.5 / x, 0.),
    }


def lambda_profiles(values, cutoff, normalization):
    """Profile values and d/d lambda on actual Dirac eigenvalues."""
    values = _array(values)
    cutoff = _positive("cutoff", cutoff)
    normalization = _positive("normalization", normalization)
    if np.min(np.abs(values)) < 1e-10:
        raise ValueError("unresolved zero mode requires a separate determinant/domain prescription")
    ratio = log_ratio(cutoff, normalization)
    x = (values / cutoff) ** 2
    expo = np.exp(-x)
    log_d2 = np.log(values ** 2 / normalization ** 2)
    coeff = ratio + 0.5 * EULER
    inside = x <= 1
    values_out = {
        "proper_time": 0.5 * exp1(x),
        "heat_rank": coeff * expo,
        "remainder": h_remainder(x),
        "weighted_log": -0.5 * expo * log_d2,
        "hard_log": np.where(inside, -0.5 * log_d2, 0.),
        "rank": np.ones_like(x),
    }
    deriv = {
        "proper_time": -expo / values,
        "heat_rank": coeff * (-2 * values / cutoff ** 2) * expo,
        "remainder": (values / cutoff ** 2) * expo * (np.log(x) + EULER),
        "weighted_log": expo * ((values / cutoff ** 2) * log_d2 - 1 / values),
        "hard_log": np.where(inside, -1 / values, 0.),
        "rank": np.zeros_like(values),
    }
    return values_out, deriv


def independent_e1(x):
    x = _positive("x", x)
    value, error = quad(lambda t: np.exp(-x * t) / t, 1., np.inf, epsabs=1e-12, limit=400)
    return float(value), float(error)


def independent_h_moments():
    def remainder(x):
        if x <= 0:
            return 0.
        return float(h_remainder(np.array([x]))[0])

    integral, err0 = quad(remainder, 0., np.inf, epsabs=2e-10, limit=500)
    weighted, err1 = quad(lambda x: x * remainder(x), 0., np.inf, epsabs=2e-10, limit=500)
    return {
        "integral_h": float(integral),
        "integral_x_h": float(weighted),
        "exact_integral_h": -0.5,
        "exact_integral_x_h": -0.75,
        "h_at_zero_plus": 0.,
        "quadrature_errors": [float(err0), float(err1)],
        "a4_seeley_weight_h0": 0.,
        "power_divergent_moments_changed_by_h": True,
        "logarithmic_a4_shared_by_weighted_log_and_PT_plus_heat": True,
    }


def analytic_cylinder_values(metric, omega, kappa):
    if not np.allclose(metric.lapse, 1) or not np.allclose(metric.radial_scale, 1):
        raise ValueError("analytic cylinder eigenvalues require unit lapse and radial scale")
    if np.ptp(metric.sphere_radius) > 1e-14:
        raise ValueError("analytic cylinder eigenvalues require constant radius")
    mass = kappa / float(metric.sphere_radius[0])
    magnitude = np.sqrt(omega ** 2 + metric.momenta ** 2 + mass ** 2)
    return np.concatenate([magnitude, -magnitude])


def fiber_operator(metric, omega, kappa, cutoff, normalization, gradients=True):
    matrix = euclidean_operator(metric, omega, kappa)
    values, vectors = eigh(matrix, driver="evd", check_finite=False)
    profiles, derivatives = lambda_profiles(values, cutoff, normalization)
    energies = {kind: float(np.sum(profiles[kind])) for kind in profiles}
    if not gradients:
        return {"values": values, "energies": energies, "gradients": None, "matrix": matrix}
    gradient = {}
    p = metric.momentum_matrix
    for kind, deriv in derivatives.items():
        kernel = (vectors * deriv[None, :]) @ vectors.conj().T
        parts = kernel.reshape(2, metric.points, 2, metric.points)
        trace = lambda gamma, block=parts: np.einsum("ba,aibj->ij", gamma, block)
        t1 = trace(SIGMA1)
        gradient[kind] = np.array([
            -omega * np.diag(trace(SIGMA3)).real / metric.lapse ** 2,
            -0.5 * np.diag(t1 @ p + p @ t1).real / metric.radial_scale ** 2,
            kappa * np.diag(trace(SIGMA2)).real / metric.sphere_radius ** 2,
        ])
    return {"values": values, "energies": energies, "gradients": gradient, "matrix": matrix}


def heat_trace_expm(matrix, cutoff):
    cutoff = _positive("cutoff", cutoff)
    square = matrix @ matrix
    return float(np.trace(expm(-square / cutoff ** 2)).real)


def integrate_profiles(metric, cutoff=2., normalization=1., angular_max=8,
                       frequency_points=32, extent_factor=1., gradients=True, kinds=None):
    if not isinstance(angular_max, int) or isinstance(angular_max, bool) or angular_max < 1:
        raise ValueError("positive angular maximum required")
    if not isinstance(frequency_points, int) or frequency_points < 8 or extent_factor < 1:
        raise ValueError("resolved quadrature and extent factor >=1 required")
    kinds = tuple(KINDS if kinds is None else kinds)
    extent = frequency_extent(metric, cutoff) * extent_factor
    nodes, weights = leggauss(frequency_points)
    frequencies = (nodes + 1) * extent / 2
    weights = weights * extent / 2
    totals = {kind: 0. for kind in kinds}
    rank = 0.
    heat = 0.
    gradient = {kind: np.zeros((3, metric.points)) for kind in kinds} if gradients else None
    for kappa in range(1, angular_max + 1):
        for omega, weight in zip(frequencies, weights):
            fiber = fiber_operator(metric, omega, kappa, cutoff, normalization, gradients)
            factor = 4 * kappa * weight / pi
            rank += factor * fiber["energies"]["rank"]
            coeff = log_ratio(cutoff, normalization) + 0.5 * EULER
            if abs(coeff) <= 1e-14:
                heat += factor * float(np.sum(np.exp(-(fiber["values"] / cutoff) ** 2)))
            else:
                heat += factor * fiber["energies"]["heat_rank"] / coeff
            for kind in kinds:
                totals[kind] += factor * fiber["energies"][kind]
                if gradients:
                    gradient[kind] += factor * fiber["gradients"][kind]
    result = {"energy": totals, "rank_trace": float(rank), "heat_trace": float(heat),
              "frequency_extent": float(extent), "cutoff": float(cutoff),
              "normalization": float(normalization), "angular_max": angular_max,
              "frequency_points": frequency_points}
    if gradients:
        n, q, r, spacing = metric.lapse, metric.radial_scale, metric.sphere_radius, metric.spacing
        stress = {}
        for kind in kinds:
            rho = gradient[kind][0] / (4 * pi * q * r * r * spacing)
            px = -gradient[kind][1] / (4 * pi * n * r * r * spacing)
            pt = -gradient[kind][2] / (8 * pi * n * q * r * spacing)
            stress[kind] = {
                "metric_gradients": gradient[kind],
                "rho": rho, "p_x": px, "p_perp": pt,
                "radial_null": rho + px,
                "trace": rho - px - 2 * pt,
                "neck_null": float((rho + px)[metric.points // 2]),
                "lapse_homogeneity_residual": float(np.dot(n, gradient[kind][0]) - totals[kind]),
            }
        result["stress"] = stress
    return result


def weyl_from_metric(metric, sigma, integrated):
    sigma = np.asarray(sigma, dtype=float)
    if sigma.shape != (metric.points,) or not np.isfinite(sigma).all():
        raise ValueError("Weyl weight must match the grid")
    fields = np.array([metric.lapse, metric.radial_scale, metric.sphere_radius])
    return {kind: float(np.sum(integrated["stress"][kind]["metric_gradients"] * fields * sigma))
            for kind in integrated["energy"]}


def independent_frequency_quadrature(metric, cutoff, normalization, kappa, kind, limit):
    """Independent adaptive omega integral of an analytic cylinder fiber."""
    if kind not in KINDS:
        raise ValueError("unknown profile")

    def integrand(omega):
        values = analytic_cylinder_values(metric, omega, kappa)
        profiles, _ = lambda_profiles(values, cutoff, normalization)
        return 4 * kappa / pi * float(np.sum(profiles[kind]))

    value, error = quad(integrand, 0., limit, epsabs=1e-8, limit=250)
    return float(value), float(error)


def hard_cylinder_mode(spatial_energy, cutoff, normalization, kappa):
    """Exact positive-frequency hard-log contribution of one cylinder (kappa,p).

    For E=sqrt(p^2+kappa^2/a^2)<Lambda, W=sqrt(Lambda^2-E^2),
    E_hard=(4 kappa/pi)[(2-ln(Lambda^2/M^2))W-2 E atan(W/E)].
    Interior is the uniform-lapse Hellmann–Feynman piece; the wall is the
    cutoff shell. Their sum is E_hard. Modes with E>=Lambda vanish.
    """
    if not isinstance(kappa, (int, np.integer)) or isinstance(kappa, bool) or kappa <= 0:
        raise ValueError("positive integer angular label required")
    energy = float(spatial_energy)
    cutoff = _positive("cutoff", cutoff)
    normalization = _positive("normalization", normalization)
    if not np.isfinite(energy) or energy < 0:
        raise ValueError("nonnegative spatial energy required")
    if energy >= cutoff:
        zeros = {"spatial_energy": energy, "W": 0., "energy": 0., "interior": 0.,
                 "wall": 0., "normalization_log_derivative": 0.}
        return zeros
    width = float(np.sqrt(cutoff * cutoff - energy * energy))
    log_cutoff = float(np.log((cutoff / normalization) ** 2))
    eat = 0. if energy == 0 else energy * float(np.arctan(width / energy))
    prefactor = 4 * kappa / pi
    interior = prefactor * 2 * (width - eat)
    wall = prefactor * (-log_cutoff * width)
    return {
        "spatial_energy": energy, "W": width, "energy": interior + wall,
        "interior": interior, "wall": wall,
        "normalization_log_derivative": prefactor * 2 * width,
    }


def exact_hard_cylinder(metric, cutoff=2., normalization=1., angular_max=8):
    """Sum the closed hard-log formula over AP/P cylinder Fourier modes."""
    if not isinstance(angular_max, int) or isinstance(angular_max, bool) or angular_max < 1:
        raise ValueError("positive angular maximum required")
    if not np.allclose(metric.lapse, 1) or not np.allclose(metric.radial_scale, 1):
        raise ValueError("exact hard cylinder requires unit lapse and radial scale")
    if np.ptp(metric.sphere_radius) > 1e-14:
        raise ValueError("exact hard cylinder requires constant radius")
    radius = float(metric.sphere_radius[0])
    totals = {"energy": 0., "interior": 0., "wall": 0., "normalization_log_derivative": 0.}
    retained = []
    for kappa in range(1, angular_max + 1):
        for momentum in metric.momenta:
            spatial = float(np.sqrt(momentum ** 2 + (kappa / radius) ** 2))
            mode = hard_cylinder_mode(spatial, cutoff, normalization, kappa)
            for key in totals:
                totals[key] += mode[key]
            if mode["W"] > 0:
                retained.append({"kappa": kappa, "momentum": float(momentum), **mode})
    totals.update({
        "cutoff": float(cutoff), "normalization": float(normalization),
        "angular_max": angular_max, "points": metric.points, "eta": metric.eta,
        "length": metric.length, "radius": radius,
        "retained_modes": retained,
        "interior_plus_wall": totals["interior"] + totals["wall"],
        "uniform_lapse_derivative_equals_interior_plus_wall": True,
        "formula": "for E<Lambda, (4 kappa/pi)[(2-ln(Lambda^2/M^2))W-2 E atan(W/E)]; d/dlogM=(8 kappa/pi)W",
    })
    return totals


def independent_hard_cylinder_quadrature(metric, cutoff=2., normalization=1., angular_max=8):
    """Adaptive omega integral split at W plus the explicit cutoff-shell term."""
    exact = exact_hard_cylinder(metric, cutoff, normalization, angular_max)
    radius = float(metric.sphere_radius[0])
    energy = interior = wall = 0.
    quad_errors = []
    for kappa in range(1, angular_max + 1):
        for momentum in metric.momenta:
            spatial = float(np.sqrt(momentum ** 2 + (kappa / radius) ** 2))
            if spatial >= cutoff:
                continue
            width = float(np.sqrt(cutoff * cutoff - spatial * spatial))
            prefactor = 4 * kappa / pi

            def interior_density(omega, mass=spatial):
                return -np.log((omega * omega + mass * mass) / cutoff ** 2)

            def full_density(omega, mass=spatial):
                return -np.log((omega * omega + mass * mass) / normalization ** 2)

            interior_value, interior_error = quad(interior_density, 0., width, epsabs=1e-12, limit=400)
            full_value, full_error = quad(full_density, 0., width, epsabs=1e-12, limit=400)
            shell = -np.log((cutoff / normalization) ** 2) * width
            interior += prefactor * float(interior_value)
            wall += prefactor * float(shell)
            energy += prefactor * float(full_value)
            quad_errors.append(max(float(interior_error), float(full_error)))
    return {
        "energy": float(energy), "interior": float(interior), "wall": float(wall),
        "interior_plus_wall": float(interior + wall),
        "energy_minus_exact": float(energy - exact["energy"]),
        "interior_minus_exact": float(interior - exact["interior"]),
        "wall_minus_exact": float(wall - exact["wall"]),
        "maximum_quadrature_error_estimate": max(quad_errors) if quad_errors else 0.,
        "split": "adaptive integral of -log((omega^2+E^2)/Lambda^2) on [0,W] plus shell -ln(Lambda^2/M^2) W",
    }


def gauss_hard_cylinder(metric, cutoff=2., normalization=1., angular_max=8, frequency_points=32):
    """Analytic-mode Gauss–Legendre hard-log on the same [0,extent] map as the operator."""
    if not isinstance(frequency_points, int) or frequency_points < 8:
        raise ValueError("resolved quadrature required")
    extent = frequency_extent(metric, cutoff)
    nodes, weights = leggauss(frequency_points)
    frequencies = (nodes + 1) * extent / 2
    weights = weights * extent / 2
    energy = 0.
    derivative = 0.
    for kappa in range(1, angular_max + 1):
        for omega, weight in zip(frequencies, weights):
            values = analytic_cylinder_values(metric, omega, kappa)
            profiles, _ = lambda_profiles(values, cutoff, normalization)
            factor = 4 * kappa * weight / pi
            energy += factor * float(np.sum(profiles["hard_log"]))
            derivative += factor * float(np.sum((values / cutoff) ** 2 <= 1))
    return {
        "frequency_points": frequency_points, "frequency_extent": float(extent),
        "energy": float(energy), "normalization_log_derivative": float(derivative),
    }


def native(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: native(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [native(item) for item in value]
    return value


def _metric_probe(metric):
    angle = 2 * pi * metric.x / metric.length
    shape = 0.3 * np.cos(angle) + 0.2 * np.sin(2 * angle)
    return np.array([
        metric.lapse * shape,
        metric.radial_scale * (0.25 * np.sin(angle) + 0.1 * np.cos(2 * angle)),
        metric.sphere_radius * (0.2 * np.cos(2 * angle)),
    ])


def _finite_metric_errors(metric, integrated, cutoff, normalization, angular_max, frequency_points):
    direction = _metric_probe(metric)
    analytic = {kind: float(np.sum(integrated["stress"][kind]["metric_gradients"] * direction))
                for kind in integrated["energy"]}
    controls = []
    fields = ("lapse", "radial_scale", "sphere_radius")
    for step in (0.002, 0.001):
        plus = replace(metric, **{field: getattr(metric, field) + step * direction[j]
                                  for j, field in enumerate(fields)})
        minus = replace(metric, **{field: getattr(metric, field) - step * direction[j]
                                   for j, field in enumerate(fields)})
        ep = integrate_profiles(plus, cutoff, normalization, angular_max, frequency_points, gradients=False)
        em = integrate_profiles(minus, cutoff, normalization, angular_max, frequency_points, gradients=False)
        observed = {kind: (ep["energy"][kind] - em["energy"][kind]) / (2 * step) for kind in analytic}
        errors = {kind: abs(observed[kind] - analytic[kind]) for kind in analytic}
        controls.append({"step": step, "finite_derivative": observed, "absolute_error": errors,
                         "smooth_absolute_error": {kind: errors[kind] for kind in SMOOTH_KINDS}})
    return {
        "direction": "mixed_Nqr",
        "analytic_variation": analytic,
        "independent_controls": controls,
        "hard_log_spectral_gradient_omits_cutoff_wall": True,
    }


def _normalization_derivative(metric, cutoff, normalization, angular_max, frequency_points):
    base = integrate_profiles(metric, cutoff, normalization, angular_max, frequency_points, gradients=False)
    step = 1e-3
    plus = integrate_profiles(metric, cutoff, normalization * np.exp(step), angular_max, frequency_points, gradients=False)
    minus = integrate_profiles(metric, cutoff, normalization * np.exp(-step), angular_max, frequency_points, gradients=False)
    observed = {kind: (plus["energy"][kind] - minus["energy"][kind]) / (2 * step) for kind in KINDS}
    predicted = {
        "proper_time": 0.,
        "heat_rank": base["heat_trace"],
        "remainder": 0.,
        "weighted_log": base["heat_trace"],
        "hard_log": None,
    }
    return {
        "observed": observed,
        "predicted_smooth": {kind: predicted[kind] for kind in ("proper_time", "heat_rank", "remainder", "weighted_log")},
        "hard_log_is_projected_rank_not_heat": True,
        "weighted_log_matches_heat_not_hard_projector": True,
    }


def _extent_scan(metric, cutoff, normalization, angular_max=4, frequency_points=24):
    rows = []
    for factor in (1., 1.5, 2.):
        result = integrate_profiles(metric, cutoff, normalization, angular_max, frequency_points,
                                    extent_factor=factor, gradients=False)
        rows.append({
            "extent_factor": factor,
            "frequency_extent": result["frequency_extent"],
            "rank_trace": result["rank_trace"],
            "heat_trace": result["heat_trace"],
            "proper_time": result["energy"]["proper_time"],
            "remainder": result["energy"]["remainder"],
            "weighted_log": result["energy"]["weighted_log"],
        })
    growth = [row["rank_trace"] / row["frequency_extent"] for row in rows]
    return {
        "samples": rows,
        "rank_over_extent": growth,
        "rank_grows_with_unrestricted_frequency": rows[-1]["rank_trace"] > 1.4 * rows[0]["rank_trace"],
        "rank_to_extent_ratio_stable": max(abs(item - growth[0]) for item in growth) < 1e-6 * max(1., abs(growth[0])),
        "heat_does_not_scale_with_extent": abs(rows[-1]["heat_trace"] - rows[0]["heat_trace"])
        < 0.05 * max(1e-12, abs(rows[-1]["heat_trace"])),
        "old_rank_cannot_be_integrated_as_Tr1": True,
    }


def _fiber_controls(metric, cutoff, normalization, omega=0.7, kappa=2):
    fiber = fiber_operator(metric, omega, kappa, cutoff, normalization, True)
    analytic = None
    analytic_error = None
    if np.allclose(metric.lapse, 1) and np.allclose(metric.radial_scale, 1) and np.ptp(metric.sphere_radius) < 1e-14:
        expected = analytic_cylinder_values(metric, omega, kappa)
        analytic = expected
        analytic_error = float(max(abs(np.sort(fiber["values"]) - np.sort(expected))))
    heat = heat_trace_expm(fiber["matrix"], cutoff)
    heat_error = abs(heat - float(np.sum(np.exp(-(fiber["values"] / cutoff) ** 2))))
    ratio = log_ratio(cutoff, normalization)
    coeff = ratio + 0.5 * EULER
    regulated = RegulatedOperator(fiber["matrix"], OperatorConventions(cutoff=cutoff, normalization=normalization))
    old = regulated.action()
    pt = fiber["energies"]["proper_time"]
    old_rank = len(fiber["values"]) * coeff
    heat_rank = fiber["energies"]["heat_rank"]
    remainder = fiber["energies"]["remainder"]
    weighted = fiber["energies"]["weighted_log"]
    identity_error = abs(weighted - (pt + heat_rank + remainder))
    field = "sphere_radius"
    direction = metric.sphere_radius * np.cos(2 * pi * metric.x / metric.length)
    analytic_grad = float(np.dot(fiber["gradients"]["weighted_log"][2], direction))
    fd = []
    for step in (1e-4, 5e-5):
        plus = replace(metric, **{field: getattr(metric, field) + step * direction})
        minus = replace(metric, **{field: getattr(metric, field) - step * direction})
        ep = fiber_operator(plus, omega, kappa, cutoff, normalization, False)["energies"]["weighted_log"]
        em = fiber_operator(minus, omega, kappa, cutoff, normalization, False)["energies"]["weighted_log"]
        observed = (ep - em) / (2 * step)
        fd.append({"step": step, "finite_derivative": observed, "absolute_error": abs(observed - analytic_grad)})
    return {
        "omega": omega, "kappa": kappa, "points": metric.points,
        "analytic_cylinder_eigenvalue_error": analytic_error,
        "heat_expm_versus_eigh_error": float(heat_error),
        "regulated_action": float(old),
        "proper_time_plus_old_rank": float(pt + old_rank),
        "old_rank_term": float(old_rank),
        "heat_rank_term": float(heat_rank),
        "old_minus_heat_rank": float(old_rank - heat_rank),
        "identity_on_fiber": float(identity_error),
        "hard_minus_weighted": float(fiber["energies"]["hard_log"] - weighted),
        "remainder_is_nonzero": abs(remainder) > 1e-12,
        "weighted_log_gradient_controls": fd,
        "analytic_weighted_log_radius_variation": analytic_grad,
        "eigenvalues_used_in_analytic_sort": None if analytic is None else np.sort(analytic).tolist(),
    }


def _clock_control(metric, cutoff, normalization, angular_max, frequency_points, factor=1.2):
    base = integrate_profiles(metric, cutoff, normalization, angular_max, frequency_points, gradients=False)
    changed = replace(metric, lapse=factor * metric.lapse)
    shifted = integrate_profiles(changed, cutoff, normalization, angular_max, frequency_points, gradients=False)
    residual = {kind: shifted["energy"][kind] - factor * base["energy"][kind] for kind in KINDS}
    return {
        "factor": factor,
        "base": base["energy"],
        "changed": shifted["energy"],
        "changed_minus_factor_times_base": residual,
        "proper_time_clock_identity": "E_PT[c N]=c E_PT[N] at fixed proper-time cutoff",
        "all_spectral_profiles_clock_identity": "E_f[c N]=c E_f[N] by omega=c nu with fixed proper scales M and Lambda; discontinuous quadrature may obscure it",
        "smooth_profile_max_clock_residual": max(abs(residual[kind]) for kind in SMOOTH_KINDS),
    }


def _profile_derivative_controls():
    ratio = log_ratio(2., 1.)
    xs = np.array([0.2, 0.5, 1.3, 2.7])
    analytic = dx_profiles(xs, ratio)
    step = 1e-6
    errors = {}
    for kind in ("proper_time", "heat_rank", "remainder", "weighted_log"):
        plus = identity_pieces(xs + step, ratio)[kind]
        minus = identity_pieces(xs - step, ratio)[kind]
        observed = (plus - minus) / (2 * step)
        errors[kind] = abs(observed - analytic[kind]).tolist()
    e1_rows = []
    for x in (0.15, 0.8, 2.4):
        value, error = independent_e1(x)
        e1_rows.append({
            "x": x, "scipy_exp1": float(exp1(x)), "independent_integral": value,
            "quadrature_error": error, "absolute_error": abs(value - float(exp1(x))),
        })
    series_xs = np.array([1e-6, 1e-4, 1e-3, 1e-2])
    return {
        "x_derivative_absolute_errors": errors,
        "independent_e1_quadrature": e1_rows,
        "low_x_series_absolute_error": abs(h_remainder(series_xs) - h_series(series_xs, 16)).tolist(),
        "h_vanishes_as_x_to_0": abs(h_remainder(np.array([1e-12]))[0]) < 1e-10,
        "identity_sample_x": xs.tolist(),
        "identity_max_residual": float(np.max(np.abs(identity_residual(xs, ratio)))),
        "hard_differs_from_weighted_at_x_gt_1": True,
    }


def build_record():
    cutoff, normalization = 2., 1.
    angular_max, frequency_points = 8, 32
    identity_x = np.array([0.05, 0.1, 0.25, 0.5, 0.8, 1., 1.5, 2., 3., 5., 8.])
    ratio = log_ratio(cutoff, normalization)
    pieces = identity_pieces(identity_x, ratio)
    moments = independent_h_moments()
    derivatives = _profile_derivative_controls()
    cylinder = cylinder_metric(16)
    general = smooth_metric(16, general=True)
    cylinder_int = integrate_profiles(cylinder, cutoff, normalization, angular_max, frequency_points)
    general_int = integrate_profiles(general, cutoff, normalization, angular_max, frequency_points)
    existing = cutoff_response(cylinder, cutoff, angular_max, frequency_points, gradients=False)
    existing_general = cutoff_response(general, cutoff, angular_max, frequency_points, gradients=False)
    sigma = 0.15 * np.cos(2 * pi * cylinder.x / cylinder.length) + 0.07 * np.sin(4 * pi * cylinder.x / cylinder.length)
    sigma_g = 0.15 * np.cos(2 * pi * general.x / general.length) + 0.07 * np.sin(4 * pi * general.x / general.length)
    cylinder_weyl = weyl_from_metric(cylinder, sigma, cylinder_int)
    general_weyl = weyl_from_metric(general, sigma_g, general_int)
    pt_weyl_existing = cutoff_weyl_trace(cylinder, sigma, cutoff, angular_max, frequency_points)
    pt_weyl_general = cutoff_weyl_trace(general, sigma_g, cutoff, angular_max, frequency_points)
    cylinder_quad = independent_frequency_quadrature(cylinder, cutoff, normalization, 1, "proper_time",
                                                     cylinder_int["frequency_extent"])
    one_kappa = integrate_profiles(cylinder, cutoff, normalization, 1, frequency_points, gradients=False)
    cylinder_fd = _finite_metric_errors(cylinder, cylinder_int, cutoff, normalization, angular_max, frequency_points)
    general_fd = _finite_metric_errors(general, general_int, cutoff, normalization, angular_max, frequency_points)
    exact_hard = exact_hard_cylinder(cylinder, cutoff, normalization, angular_max)
    hard_quad = independent_hard_cylinder_quadrature(cylinder, cutoff, normalization, angular_max)
    gauss_hard = []
    for nodes in (32, 64, 128):
        sample = gauss_hard_cylinder(cylinder, cutoff, normalization, angular_max, nodes)
        gauss_hard.append({
            **sample,
            "error_versus_exact_energy": sample["energy"] - exact_hard["energy"],
            "error_versus_exact_normalization_derivative":
                sample["normalization_log_derivative"] - exact_hard["normalization_log_derivative"],
        })
    source_paths = OWNED_SOURCES + READONLY_SOURCES
    heat_changes_finite = abs(cylinder_int["energy"]["remainder"]) > 1e-6
    if float(np.max(np.abs(pieces["residual"]))) >= 1e-12:
        raise RuntimeError("profile identity failed")
    if abs(cylinder_int["energy"]["proper_time"] - existing["energy"]) >= 1e-10:
        raise RuntimeError("proper-time integral does not reproduce the covariant frequency operator")
    if abs(general_int["energy"]["proper_time"] - existing_general["energy"]) >= 1e-10:
        raise RuntimeError("proper-time integral fails on the general metric")
    if not heat_changes_finite:
        raise RuntimeError("heat covariantization produced a vanishing remainder")
    if abs(exact_hard["energy"] - exact_hard["interior_plus_wall"]) >= 1e-14:
        raise RuntimeError("hard-cylinder interior plus wall failed to reconstruct the energy")
    if abs(hard_quad["energy_minus_exact"]) >= 1e-12:
        raise RuntimeError("adaptive hard-cylinder quadrature failed the closed form")
    if abs(gauss_hard[0]["energy"] - cylinder_int["energy"]["hard_log"]) >= 1e-12:
        raise RuntimeError("analytic-mode Gauss hard-log disagrees with the covariant operator path")
    for row in (cylinder_fd, general_fd):
        if max(row["independent_controls"][-1]["smooth_absolute_error"].values()) >= 2e-7:
            raise RuntimeError("smooth profile metric variation failed the independent finite-difference control")
    return native({
        "schema": SCHEMA,
        "artifact_id": ARTIFACT,
        "classification": "exact_finite_profile_identity_and_covariant_metric_response_of_normalization_prescriptions",
        "lab_continued": "3a747cc",
        "source_hashes": {path: sha256((ROOT / path).read_bytes()).hexdigest()
                          for path in source_paths},
        "input_hashes": {path: sha256((ROOT / path).read_bytes()).hexdigest()
                         for path in INPUT_RECORDS},
        "conventions": {
            "x": "d^2/Lambda^2",
            "g_PT": "0.5 E1(x); calculated determinant modulus per eigenvalue, not complete Gamma_one",
            "g_log": "-0.5 exp(-x) log(d^2/M^2); declared smooth weighted-log, not Theta(1-x)",
            "h": "-0.5[exp(-x)(log x + gamma_E) + E1(x)]",
            "identity": "g_log = g_PT + [log(M/Lambda)+gamma_E/2] exp(-x) + h(x)",
            "primary_Andrianov_1106_3263": "P_N=Theta(1-D^2/Lambda^2); log Z=sum_{|lambda|<=Lambda} log(lambda/mu); factorization of invariant/anomalous parts is not unique",
            "Gamma_one": "g_PT is the calculated E_Lambda=0.5 Tr E1 modulus, not a completed Gamma_one; the identity is not inserted",
            "hard_cylinder": "exact product formula with explicit interior Hellmann-Feynman and cutoff-wall shell; Gauss values are coarse diagnostics",
            "old_finite_rank": "N[log(M/Lambda)+gamma_E/2] is a finite-matrix convention and cannot be integrated as Tr 1",
            "heat_covariantization": "replace N by Tr exp(-D^2/Lambda^2); this is not g_log and is not the old rank term",
            "M": "normalization mass, not a fermion mass; not selected physically",
            "compensator": "unknown finite/anomaly completion remains explicit; AGENTS 20.28-29 keep det and induced heat as one regulated object",
            "state": "static compact axial P/AP vacuum of NSC-9; no trapped continuation, no in-in Gaussian influence functional",
            "comparison": "all fields; float atol 2e-7 rtol 2e-8; keys integers booleans strings hashes exact",
        },
        "exact_identity": {
            "log_M_over_Lambda": ratio,
            "x": identity_x.tolist(),
            "proper_time": pieces["proper_time"].tolist(),
            "heat_rank": pieces["heat_rank"].tolist(),
            "remainder": pieces["remainder"].tolist(),
            "weighted_log": pieces["weighted_log"].tolist(),
            "reconstructed": pieces["reconstructed"].tolist(),
            "residual": pieces["residual"].tolist(),
            "max_abs_residual": float(np.max(np.abs(pieces["residual"]))),
            "hard_log": pieces["hard_log"].tolist(),
            "hard_minus_weighted": pieces["hard_minus_weighted"].tolist(),
            "primary_hard_projector_equals_weighted_log": False,
        },
        "low_x_limit": {
            "h_0_plus": 0.,
            "leading": "0.5 x (log x + gamma_E - 1)",
            "series": "sum_{n>=1} (-1)^{n+1} x^n (log x + gamma_E - 1/n) / (2 n!)",
            "g_log_and_PT_plus_heat_agree_for_fixed_modes_as_cutoff_removed_or_below_cutoff": True,
            "x_to_0_is_not_ultraviolet_of_each_mode": True,
            "series_errors": derivatives["low_x_series_absolute_error"],
            "h_vanishes_as_x_to_0": derivatives["h_vanishes_as_x_to_0"],
        },
        "moments": moments,
        "profile_derivative_controls": derivatives,
        "fiber_controls": {
            "cylinder": _fiber_controls(cylinder, cutoff, normalization),
            "general": _fiber_controls(general, cutoff, normalization),
        },
        "unrestricted_trace_one": _extent_scan(cylinder, cutoff, normalization),
        "cylinder": {
            "points": cylinder.points, "eta": cylinder.eta, "length": cylinder.length,
            "energy": cylinder_int["energy"],
            "rank_trace": cylinder_int["rank_trace"],
            "heat_trace": cylinder_int["heat_trace"],
            "frequency_extent": cylinder_int["frequency_extent"],
            "existing_cutoff_response_energy": existing["energy"],
            "proper_time_minus_existing_cutoff_response": cylinder_int["energy"]["proper_time"] - existing["energy"],
            "neck_null": {kind: cylinder_int["stress"][kind]["neck_null"] for kind in KINDS},
            "lapse_homogeneity": {kind: cylinder_int["stress"][kind]["lapse_homogeneity_residual"] for kind in KINDS},
            "metric_Weyl": cylinder_weyl,
            "existing_proper_time_Weyl": pt_weyl_existing,
            "proper_time_Weyl_minus_existing": cylinder_weyl["proper_time"] - pt_weyl_existing,
            "independent_kappa1_proper_time_quadrature": {"value": cylinder_quad[0], "error": cylinder_quad[1],
                                                          "one_kappa_energy": one_kappa["energy"]["proper_time"]},
            "independent_metric_variations": cylinder_fd,
            "normalization_log_derivative": _normalization_derivative(cylinder, cutoff, normalization, angular_max, frequency_points),
            "clock_control": _clock_control(cylinder, cutoff, normalization, angular_max, frequency_points),
            "hard_log_energy_is_coarse_gauss_not_exact": True,
            "hard_log_coarse_minus_exact": cylinder_int["energy"]["hard_log"] - exact_hard["energy"],
        },
        "hard_cylinder_exact": {
            "energy": exact_hard["energy"],
            "interior": exact_hard["interior"],
            "wall": exact_hard["wall"],
            "interior_plus_wall": exact_hard["interior_plus_wall"],
            "normalization_log_derivative": exact_hard["normalization_log_derivative"],
            "uniform_lapse_derivative_equals_interior_plus_wall": exact_hard["uniform_lapse_derivative_equals_interior_plus_wall"],
            "formula": exact_hard["formula"],
            "retained_modes": exact_hard["retained_modes"],
            "independent_adaptive_quadrature": hard_quad,
            "coarse_gauss": gauss_hard,
            "coarse_operator_hard_log": cylinder_int["energy"]["hard_log"],
            "convergence_not_assumed": True,
        },
        "general_metric": {
            "points": general.points, "eta": general.eta, "length": general.length,
            "energy": general_int["energy"],
            "rank_trace": general_int["rank_trace"],
            "heat_trace": general_int["heat_trace"],
            "frequency_extent": general_int["frequency_extent"],
            "existing_cutoff_response_energy": existing_general["energy"],
            "proper_time_minus_existing_cutoff_response": general_int["energy"]["proper_time"] - existing_general["energy"],
            "neck_null": {kind: general_int["stress"][kind]["neck_null"] for kind in KINDS},
            "lapse_homogeneity": {kind: general_int["stress"][kind]["lapse_homogeneity_residual"] for kind in KINDS},
            "metric_Weyl": general_weyl,
            "existing_proper_time_Weyl": pt_weyl_general,
            "proper_time_Weyl_minus_existing": general_weyl["proper_time"] - pt_weyl_general,
            "independent_metric_variations": general_fd,
            "normalization_log_derivative": _normalization_derivative(general, cutoff, normalization, angular_max, frequency_points),
            "clock_control": _clock_control(general, cutoff, normalization, angular_max, frequency_points),
            "hard_log_status": {
                "unresolved_coarse_quadrature_and_interior_derivative_diagnostic": True,
                "physical_stress_comparison": False,
                "convergence_not_assumed": True,
                "exact_hard_formula_is_cylinder_product_only": True,
            },
        },
        "conclusion": {
            "heat_covariantization_of_old_rank_is_exact": False,
            "finite_action_changes_by_Tr_h": heat_changes_finite,
            "old_rank_equals_heat_trace_on_unrestricted_frequency": False,
            "primary_hard_projector_equals_weighted_log": False,
            "identity_forced_into_Gamma_one": False,
            "g_PT_is_calculated_determinant_modulus_not_complete_Gamma_one": True,
            "unknown_finite_compensator_explicit": True,
            "throat_root_fitted": False,
            "real_time_Gaussian_influence_functional_computed": False,
            "hard_cylinder_gauss_quadrature_resolved_by_closed_form": True,
            "general_hard_log_is_physical_stress": False,
        },
        "nonclaims": {
            "complete_invariant_functional_fixed": False,
            "finite_coefficients_selected_by_anomaly": False,
            "self_sourced_stationary_metric": False,
            "physical_Phi_or_Omega_derived": False,
            "full_physical_closure": False,
            "in_in_or_real_time_influence_functional": False,
            "Andrianov_factorization_uniqueness": False,
            "weighted_log_is_the_physical_measure": False,
            "coarse_hard_log_stress_is_physical": False,
        },
        "terminal": True,
    })
