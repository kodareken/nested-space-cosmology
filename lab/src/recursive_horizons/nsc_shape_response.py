"""UV-finite spin-structure response on a smooth compact axial geometry.

Static signature (+---): ds²=N²dt²-q²dx²-r²dOmega², x~x+L.
The radial normalization u=r*sqrt(q)*psi gives
H=-i sigma2[c partial_x+c'/2]+v sigma1, c=N/q, v=N*kappa/r.
The lattice below discretizes this spatial Hamiltonian and is removed by
resolution checks. Its spectral cutoff is not a covariant physical Lambda
when N varies. No absolute AP stress or stationary metric is computed.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import pi

import numpy as np
from scipy.linalg import svd
import sympy as sp


@lru_cache(maxsize=1)
def radial_reduction_identities():
    """Check the spin connection and spherical/coordinate measure removal."""
    x = sp.symbols("x", real=True)
    n, q, r, u = (sp.Function(name)(x) for name in ("N", "q", "r", "u"))
    c = n / q
    psi = u / (r * sp.sqrt(q))
    # The diagonal tetrad gives c(∂+r'/r+N'/(2N)) before u=r sqrt(q) psi.
    transformed = sp.simplify(r * sp.sqrt(q) * c * (
        sp.diff(psi, x) + (sp.diff(r, x) / r + sp.diff(n, x) / (2 * n)) * psi))
    expected = c * sp.diff(u, x) + sp.diff(c, x) * u / 2
    return {"transformed_derivative": str(transformed),
            "symmetric_derivative_residual": str(sp.simplify(transformed - expected)),
            "measure_residual": str(sp.simplify(q * r**2 * psi**2 - u**2)),
            "Hilbert_measure_before": "q*r^2 dx dOmega for psi",
            "Hilbert_measure_after_angular_separation": "dx for u=r*sqrt(q)*psi"}


@dataclass(frozen=True)
class StaticAxialMetric:
    length: float
    lapse: np.ndarray
    radial_scale: np.ndarray
    sphere_radius: np.ndarray

    def __post_init__(self):
        if isinstance(self.length, bool) or not np.isfinite(self.length) or self.length <= 0:
            raise ValueError("positive finite coordinate period required")
        arrays = []
        for name in ("lapse", "radial_scale", "sphere_radius"):
            array = np.array(getattr(self, name), dtype=float, copy=True)
            if array.ndim != 1 or len(array) < 8 or len(array) % 2:
                raise ValueError("metric arrays must have an equal even grid size >=8")
            if not np.isfinite(array).all() or np.min(array) <= 0:
                raise ValueError("all metric functions must be finite and positive")
            array.setflags(write=False)
            object.__setattr__(self, name, array)
            arrays.append(array)
        if len({len(a) for a in arrays}) != 1:
            raise ValueError("metric grid sizes differ")

    @property
    def points(self):
        return len(self.lapse)

    @property
    def spacing(self):
        return self.length / self.points

    @property
    def x(self):
        return (np.arange(self.points) - self.points // 2) * self.spacing


def smooth_metric(radius_parameter, points, throat_radius=1.):
    """One compact motif, not an unwrapped Bloch array or transverse y circle."""
    if not np.isfinite(radius_parameter) or radius_parameter <= 0:
        raise ValueError("R must be finite and positive")
    if not np.isfinite(throat_radius) or throat_radius <= 0:
        raise ValueError("throat radius must be finite and positive")
    if isinstance(points, bool) or not isinstance(points, int) or points < 8 or points % 2:
        raise ValueError("use an even grid of at least eight points")
    length = 2 * radius_parameter
    x = (np.arange(points) - points // 2) * length / points
    wave = pi / (2 * radius_parameter)
    radius = np.sqrt(throat_radius**2 + (np.sin(wave * x) / wave)**2)
    return StaticAxialMetric(length, np.ones(points), np.ones(points), radius)


def staggered_operator(metric, kappa, eta):
    """Node-to-edge first-order block and coefficients used by its variation."""
    if isinstance(kappa, bool) or not isinstance(kappa, int) or kappa < 1:
        raise ValueError("kappa must be a positive spinor-sphere integer")
    if eta not in (0., .5):
        raise ValueError("only periodic/AP spin structures are declared")
    c = metric.lapse / metric.radial_scale
    v = metric.lapse * kappa / metric.sphere_radius
    # This guards the unresolved regime where the centered mass overwhelms
    # the lattice derivative. It is not a continuum physical UV cutoff.
    if np.max(metric.spacing * v / c) >= 2:
        raise ValueError("resolve angular potential: require h*kappa*q/r < 2")
    cnext = np.roll(c, -1)
    vedge = (v + np.roll(v, -1)) / 4
    cedge = (c + cnext) / 2
    left_derivative = np.sqrt(cedge * c) / metric.spacing
    right_derivative = np.sqrt(cedge * cnext) / metric.spacing
    phase = np.ones(metric.points)
    phase[-1] = 1. if eta == 0 else -1.
    i = np.arange(metric.points)
    j = (i + 1) % metric.points
    matrix = np.zeros((metric.points, metric.points))
    matrix[i, i] = -left_derivative + vedge
    matrix[i, j] = (right_derivative + vedge) * phase
    return matrix, (c, v, cedge, left_derivative, right_derivative, phase)


def _singular_response(matrix, phase):
    """Singular values and the only two polar-factor diagonals needed by HF."""
    u, values, vh = svd(matrix, full_matrices=False, check_finite=False, lapack_driver="gesdd")
    i = np.arange(len(matrix))
    j = (i + 1) % len(matrix)
    # d Tr sqrt(A†A)=Re Tr[(U V†)† dA]. Real P/AP phases are exact.
    left = np.einsum("ij,ji->i", u, vh)
    right = np.einsum("ij,ji->i", u, vh[:, j]) * phase
    return values, left, right


def _coefficient_gradients(metric, coefficients, left, right):
    c, _, cedge, tl, tr, _ = coefficients
    cn = np.roll(c, -1)
    own = left * (-tl / (4 * cedge) - tl / (2 * c)) + right * tr / (4 * cedge)
    following = -left * tl / (4 * cedge) + right * (tr / (4 * cedge) + tr / (2 * cn))
    gradient_c = own + np.roll(following, 1)
    contribution = (left + right) / 4
    gradient_v = contribution + np.roll(contribution, 1)
    return gradient_c, gradient_v


def spin_shape_difference(metric, angular_max=8, include_stress=True):
    """E_P-E_AP and every metric derivative, before setting N=q=1.

    Full four-component counting is E=-4 sum_kappa kappa sum singular(A).
    P/AP values and polar responses are subtracted within each angular sector
    before accumulating, reducing cancellation of their extensive UV terms.
    """
    if isinstance(angular_max, bool) or not isinstance(angular_max, int) or angular_max < 1:
        raise ValueError("angular_max must be a positive integer")
    g_n = np.zeros(metric.points)
    g_q = np.zeros(metric.points)
    g_r = np.zeros(metric.points)
    sectors = []
    for kappa in range(1, angular_max + 1):
        responses = []
        for eta in (0., .5):
            matrix, coeff = staggered_operator(metric, kappa, eta)
            if include_stress:
                values, left, right = _singular_response(matrix, coeff[-1])
                responses.append((values, left, right))
            else:
                values = svd(matrix, compute_uv=False, check_finite=False, lapack_driver="gesdd")
                responses.append((values,))
        weight = -4 * kappa
        energy = weight * float(np.sum(responses[0][0] - responses[1][0]))
        sector = {"kappa": kappa, "energy_difference": energy,
                  "minimum_singular_value": float(min(r[0][-1] for r in responses))}
        if include_stress:
            gc, gv = _coefficient_gradients(metric, coeff,
                                            responses[0][1] - responses[1][1],
                                            responses[0][2] - responses[1][2])
            gn = weight * (gc / metric.radial_scale + gv * kappa / metric.sphere_radius)
            gq = -weight * gc * metric.lapse / metric.radial_scale**2
            gr = -weight * gv * metric.lapse * kappa / metric.sphere_radius**2
            g_n += gn
            g_q += gq
            g_r += gr
            sector["maximum_metric_gradient_magnitude"] = float(max(np.max(np.abs(g)) for g in (gn, gq, gr)))
        sectors.append(sector)
    result = {"energy_difference": float(sum(s["energy_difference"] for s in sectors)),
              "angular_sectors": sectors,
              "maximum_lattice_mass_ratio": float(np.max(metric.spacing * angular_max * metric.radial_scale / metric.sphere_radius))}
    if include_stress:
        n, q, r, h = metric.lapse, metric.radial_scale, metric.sphere_radius, metric.spacing
        density = g_n / (4 * pi * q * r**2 * h)
        radial = -g_q / (4 * pi * n * r**2 * h)
        transverse = -g_r / (8 * pi * n * q * r * h)
        result.update({"gradient_lapse": g_n, "gradient_radial_scale": g_q,
                       "gradient_sphere_radius": g_r,
                       "density_difference": density, "axial_pressure_difference": radial,
                       "transverse_pressure_difference": transverse,
                       "axial_null_difference": density + radial,
                       "local_trace_difference": density - radial - 2 * transverse,
                       "local_weyl_identity": n * g_n + q * g_q + r * g_r,
                       "global_lapse_homogeneity_residual": float(np.dot(n, g_n) - result["energy_difference"])})
    return result


def stress_summary(metric, response):
    """Diagnostics use the (+---) trace and physical stress conservation."""
    h = metric.spacing
    derivative = lambda f: (np.roll(f, -1) - np.roll(f, 1)) / (2 * h)
    rho, px, pt = (response[key] for key in ("density_difference", "axial_pressure_difference", "transverse_pressure_difference"))
    conservation = (derivative(px) + derivative(metric.lapse) / metric.lapse * (rho + px)
                    + 2 * derivative(metric.sphere_radius) / metric.sphere_radius * (px - pt))
    center = metric.points // 2
    null = rho + px
    return {"points": metric.points, "spacing": h, "energy_difference": response["energy_difference"],
            "throat_density_difference": float(rho[center]),
            "throat_axial_pressure_difference": float(px[center]),
            "throat_transverse_pressure_difference": float(pt[center]),
            "throat_axial_null_difference": float(null[center]),
            "proper_axial_integral_null_difference": float(h * np.dot(metric.radial_scale, null)),
            "sphere_weighted_proper_integral_null_difference": float(4 * pi * h * np.dot(metric.radial_scale * metric.sphere_radius**2, null)),
            "maximum_absolute_local_trace": float(np.max(np.abs(response["local_trace_difference"]))),
            "maximum_absolute_local_weyl_identity": float(np.max(np.abs(response["local_weyl_identity"]))),
            "maximum_absolute_conservation_residual": float(np.max(np.abs(conservation))),
            "rms_conservation_residual": float(np.sqrt(np.mean(conservation**2))),
            "global_lapse_homogeneity_residual": response["global_lapse_homogeneity_residual"],
            "maximum_lattice_mass_ratio": response["maximum_lattice_mass_ratio"]}
