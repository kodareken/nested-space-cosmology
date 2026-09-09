"""Finite compact-boundary interaction of the existing covariant Dirac operator.

The imported interval determinant gives f(D4)=-1/2 log(1-s exp(-2 ell |D4|)).
This module evaluates that nonlocal interaction and its metric variations.
Local determinant terms, the conformal boundary cocycle, phase and link/gauge
backgrounds are not included. The transverse interval ell is distinct from
the axial-circle length and its P/AP spin structure.
"""
from math import pi, sqrt

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh
from scipy.special import erfc, zeta

from .nsc_covariant_operator import (
    SIGMA1, SIGMA2, SIGMA3, euclidean_operator, frequency_extent, radial_momentum,
)


def spectral_profile(values, interval, twist=1):
    if not np.isfinite(interval) or interval <= 0 or twist not in (-1, 1):
        raise ValueError("positive conformal interval and twist +1/-1 required")
    absolute = np.abs(values)
    if np.min(absolute) < 1e-11:
        raise ValueError("zero modes need a separate determinant and variation prescription")
    argument = 2 * interval * absolute
    exponential = np.exp(-argument)
    denominator = -np.expm1(-argument) if twist == 1 else 1 + exponential
    occupancy = twist * exponential / denominator
    if twist == 1:
        profile = -.5 * np.where(argument < .69, np.log(denominator), np.log1p(-exponential))
    else:
        profile = -.5 * np.log1p(exponential)
    derivative = -interval * np.sign(values) * occupancy
    interval_derivative = -absolute * occupancy
    return profile, derivative, interval_derivative


def frequency_response(metric, omega, kappa, interval, twist=1, gradients=True):
    values, vectors = eigh(euclidean_operator(metric, omega, kappa),
                           driver="evd", check_finite=False)
    profile, derivative, interval_derivative = spectral_profile(values, interval, twist)
    if not gradients:
        return float(np.sum(profile)), float(np.sum(interval_derivative)), None
    response = (vectors * derivative[None, :]) @ vectors.conj().T
    blocks = response.reshape(2, metric.points, 2, metric.points)
    trace = lambda gamma: np.einsum("ba,aibj->ij", gamma, blocks)
    gn = -omega * np.diag(trace(SIGMA3)).real / metric.lapse ** 2
    t, p = trace(SIGMA1), metric.momentum_matrix
    gq = -.5 * np.diag(t @ p + p @ t).real / metric.radial_scale ** 2
    gr = kappa * np.diag(trace(SIGMA2)).real / metric.sphere_radius ** 2
    return float(np.sum(profile)), float(np.sum(interval_derivative)), np.array([gn, gq, gr])


def energy_tail_bound(metric, interval, angular_max, extent, copies):
    """Bound the frequency tail at fixed radial/angular truncation.

    |f(lambda)| <= exp(-2 ell |lambda|)/(2(1-exp(-2 ell |lambda|))).
    The existing bound |D_omega| >= sqrt(omega²/Nmax²-c*omega)
    increases at least as fast as omega/Nmax beyond the chosen endpoint.
    The bound is evaluated in binary64, not interval arithmetic. It does
    not bound derivative, radial, angular or domain-size errors.
    """
    pq, inverse_n = radial_momentum(metric), 1 / metric.lapse
    c = float(np.linalg.norm(pq * inverse_n[None, :] - inverse_n[:, None] * pq, 2))
    nmax = float(max(metric.lapse))
    lower = extent ** 2 / nmax ** 2 - c * extent
    if lower <= 0:
        raise ValueError("frequency endpoint has no positive spectral lower bound")
    x = np.exp(-2 * interval * sqrt(lower))
    sectors = angular_max * (angular_max + 1) / 2
    return float(copies * 2 * metric.points * sectors * nmax / (pi * interval) * x / (1 - x))


def compact_response(metric, interval=2., angular_max=8, frequency_points=48,
                     extent_factor=1., copies=2, twist=1, gradients=True):
    """Interaction energy per g4 coordinate time and effective 4D metric source.

    The source is integrated over the transverse interval. It is not a
    pointwise 5D stress or a boundary-localized surface tension.
    """
    if not np.isfinite(interval) or interval <= 0 or twist not in (-1, 1):
        raise ValueError("positive conformal interval and twist +1/-1 required")
    if not isinstance(copies, int) or isinstance(copies, bool) or copies < 1:
        raise ValueError("positive integer number of 5D Dirac copies required")
    if not isinstance(angular_max, int) or angular_max < 1:
        raise ValueError("positive angular maximum required")
    if not isinstance(frequency_points, int) or frequency_points < 8 or extent_factor < 1:
        raise ValueError("resolved frequency quadrature and extent factor >=1 required")
    # Reuse the covariant spectral endpoint estimate: |D| >= 16/ell here.
    extent = frequency_extent(metric, 2 / interval) * extent_factor
    nodes, weights = leggauss(frequency_points)
    frequencies, weights = (nodes + 1) * extent / 2, weights * extent / 2
    energy = derivative = 0.
    gradient = np.zeros((3, metric.points))
    sectors = []
    for kappa in range(1, angular_max + 1):
        sector = 0.
        for omega, weight in zip(frequencies, weights):
            value, dl, g = frequency_response(metric, omega, kappa, interval, twist, gradients)
            factor = copies * 4 * kappa * weight / pi
            sector += factor * value
            derivative += factor * dl
            if gradients:
                gradient += factor * g
        energy += sector
        sectors.append({"kappa": kappa, "energy": float(sector)})
    result = {"energy": float(energy), "interval_derivative": float(derivative),
              "angular_sectors": sectors, "frequency_extent": float(extent),
              "finite_angular_frequency_energy_tail_bound":
                  energy_tail_bound(metric, interval, angular_max, extent, copies)}
    if gradients:
        n, q, r, dx = metric.lapse, metric.radial_scale, metric.sphere_radius, metric.spacing
        rho = gradient[0] / (4 * pi * q * r ** 2 * dx)
        px = -gradient[1] / (4 * pi * n * r ** 2 * dx)
        pt = -gradient[2] / (8 * pi * n * q * r * dx)
        result.update({"metric_gradients": gradient, "rho_4D": rho, "p_radial_4D": px,
                       "p_sphere_4D": pt, "radial_null_4D": rho + px,
                       "lapse_homogeneity_residual": float(np.dot(n, gradient[0]) - energy)})
    return result


def matched_local_terms(metric, interval=2., copies=2):
    """Match the first two curvature terms of this SAME finite interaction.

    E_R uses the inherited (+---) scalar-curvature energy convention.
    Keeping full_response-local_terms is exact accounting, not a truncation
    justified at arbitrary ell²*curvature. The IR-singular higher moments
    are not treated as a local derivative series.
    """
    if not np.isfinite(interval) or interval <= 0 or not isinstance(copies, int) or isinstance(copies, bool) or copies < 1:
        raise ValueError("positive interval and integer copies required")
    from .nsc_finite_terms import finite_response
    from .nsc_shape_response import StaticAxialMetric
    local = finite_response(StaticAxialMetric(
        metric.length, metric.lapse, metric.radial_scale, metric.sphere_radius))
    rank = 4 * copies
    c0 = float(3 * rank * zeta(5) / (128 * pi ** 2 * interval ** 4))
    cr = float(rank * zeta(3) / (768 * pi ** 2 * interval ** 2))
    fields = {
        "energy": c0 * local["energy"][0] + cr * local["energy"][1],
        "rho_4D": c0 * local["rho"][0] + cr * local["rho"][1],
        "p_radial_4D": c0 * local["p_x"][0] + cr * local["p_x"][1],
        "p_sphere_4D": c0 * local["p_perp"][0] + cr * local["p_perp"][1],
        "radial_null_4D": cr * local["null"][1],
    }
    return {"volume_coefficient": c0, "curvature_coefficient": cr,
            "source": fields, "curvature_radial_null": cr * local["null"][1]}


def holonomy_cutoff_energy(metric, phase, interval=2., cutoff=3.,
                          angular_max=32, compact_max=16, copies=2, derivatives=False):
    """Full free fermion determinant's holonomy-dependent energy control.

    This is not just the compact interaction. Reuse the already derived
    ultrastatic frequency integral for every KK level, including the
    half-rank zero level. Local counterterms are independent of the flat
    axial U(1) connection; only energy differences have a removal-of-cutoff
    interpretation here. The phase combines spin structure and Wilson line.
    """
    if not np.all(metric.lapse == 1):
        raise ValueError("this holonomy control requires unit ultrastatic lapse")
    if not all(np.isfinite(x) for x in (phase, interval, cutoff)) or interval <= 0 or cutoff <= 0:
        raise ValueError("finite phase and positive interval/cutoff required")
    if any(not isinstance(x, int) or isinstance(x, bool) or x < 1
           for x in (angular_max, compact_max, copies)):
        raise ValueError("positive integer truncations and copies required")
    masses = np.arange(compact_max + 1) * pi / interval
    weights = np.ones(compact_max + 1)
    weights[0] = .5
    gauge_vertex = np.kron(SIGMA1, np.diag(2 * pi / (metric.length * metric.radial_scale)))
    energy = first = second = 0.
    for kappa in range(1, angular_max + 1):
        matrix = euclidean_operator(metric, 0., kappa) + (phase - metric.eta) * gauge_vertex
        if derivatives:
            values, vectors = eigh(matrix, driver="evd", check_finite=False)
        else:
            values = np.linalg.eigvalsh(matrix)
        absolute = np.sqrt(values[:, None] ** 2 + masses[None, :] ** 2)
        if np.min(absolute) < 1e-11:
            raise ValueError("unresolved zero mode in holonomy control")
        exponential = np.exp(-(absolute / cutoff) ** 2)
        complementary = erfc(absolute / cutoff)
        profile = (cutoff / sqrt(pi) * exponential - absolute * complementary) @ weights
        factor = copies * 2 * kappa  # all signed radial eigenvalues
        energy += factor * np.sum(profile)
        if derivatives:
            fp = (-values[:, None] / absolute * complementary) @ weights
            fpp = (-masses[None, :] ** 2 / absolute ** 3 * complementary
                   + 2 * values[:, None] ** 2 / (cutoff * sqrt(pi) * absolute ** 2) * exponential) @ weights
            vertex = vectors.conj().T @ gauge_vertex @ vectors
            first += factor * np.dot(fp, np.diag(vertex).real)
            separation = values[:, None] - values[None, :]
            near = np.abs(separation) <= 1e-9 * np.maximum(1., np.maximum(np.abs(values[:, None]), np.abs(values[None, :])))
            divided = np.empty_like(separation)
            np.divide(fp[:, None] - fp[None, :], separation, out=divided, where=~near)
            average = .5 * (fpp[:, None] + fpp[None, :])
            divided[near] = average[near]
            second += factor * np.sum(divided * np.abs(vertex) ** 2)
    result = {"energy": float(energy)}
    if derivatives:
        result.update({"phase_derivative": float(first), "phase_second_derivative": float(second)})
    return result
