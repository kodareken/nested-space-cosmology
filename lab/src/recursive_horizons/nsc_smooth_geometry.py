"""Seam-free periodic areal radius, spatial Dirac gap, and required Einstein stress.

The controlled profile is

    r(x) = sqrt(a^2 + [sin(k x)/k]^2),   k = pi/(2R),   x in [-R, R],

with a and R explicit inputs. For a = 1 the throat has r_min = 1 and r''(0) = 1,
while r'(±R) = 0, so the 2R-periodic extension is C^∞. The previous repeated
sqrt(1+rho^2) motif is continuous but not C^1.

This module solves the kappa = 1 spatial Dirac operator with w = kappa/r on that
profile, and derives the 4D ultrastatic Einstein stress required by
ds^2 = dt^2 - dx^2 - r(x)^2 dOmega^2. The geometry is not a stationary GR
solution. The stress is the Einstein residual of the imposed metric, not the
derived matter stress of S_one or of the full spectral/5D theory. No phantom
field or fitted coupling is introduced.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.special import ellipk
import sympy as sp

from recursive_horizons.nsc_geometric_chain import continuum_first_band_edge as seamed_continuum_first_band_edge


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULT = ROOT / "results/nsc-4-smooth-geometry.json"
SCHEMA = "nsc-smooth-geometry-v1"
ARTIFACT = "NSC-4-SMOOTH-GEOMETRY"
OWNED_SOURCE_PATHS = (
    "src/recursive_horizons/nsc_smooth_geometry.py",
    "scripts/check_nsc_smooth_geometry.py",
    "tests/test_nsc_smooth_geometry.py",
    "docs/nsc-smooth-geometry.md",
)
PUBLISHED_COMMITS = (
    "76975387b9b03e323c1832237e0270dc9d49a382",
    "8e745d2744b1c8251b439d8ddca0db2a9740c653",
)
FROZEN_RESULT_HASH_FIELDS = (
    ("results/nsc-3-geometric-chain.json", "source_hashes"),
    ("results/nsc-3-regulated-recursion.json", "source_hashes"),
)
FROZEN_SINGLE_HASH = (
    ("results/nsc-3-radial-spectrum.json", "source_sha256", "scripts/check_nsc_radial_spectrum.py"),
)
THROAT = 1.0
KAPPA = 1.0
RADII = (2.0, 4.0, 8.0)
INTERVALS = (32, 64, 128)
PHASES = 33
ODE_RTOL = 2.0e-11
ODE_ATOL = 2.0e-12
ODE_MAX_STEP = 0.15


def wave_number(radius: float) -> float:
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("positive finite motif radius required")
    return float(np.pi / (2.0 * radius))


def areal_radius(x, throat: float = THROAT, radius: float = 2.0):
    if not np.isfinite(throat) or throat <= 0:
        raise ValueError("positive finite throat radius a required")
    k = wave_number(radius)
    return np.sqrt(throat * throat + np.square(np.sin(k * x) / k))


def areal_radius_derivative(x, throat: float = THROAT, radius: float = 2.0):
    k = wave_number(radius)
    r = areal_radius(x, throat, radius)
    return np.sin(2.0 * k * x) / (2.0 * k * r)


def areal_radius_second(x, throat: float = THROAT, radius: float = 2.0):
    k = wave_number(radius)
    r = areal_radius(x, throat, radius)
    rp = areal_radius_derivative(x, throat, radius)
    return (np.cos(2.0 * k * x) - rp * rp) / r


def areal_radius_third(x, throat: float = THROAT, radius: float = 2.0):
    k = wave_number(radius)
    r = areal_radius(x, throat, radius)
    rp = areal_radius_derivative(x, throat, radius)
    rpp = areal_radius_second(x, throat, radius)
    return (-2.0 * k * np.sin(2.0 * k * x) - 3.0 * rp * rpp) / r


def superpotential(x, throat: float = THROAT, radius: float = 2.0, kappa: float = KAPPA):
    if not np.isfinite(kappa) or kappa < 0:
        raise ValueError("nonnegative finite kappa required")
    return kappa / areal_radius(x, throat, radius)


def seamed_radius_derivative_jump(radius: float) -> float:
    """Magnitude of the old r' jump; the right-minus-left signed jump is negative."""
    return float(2.0 * radius / np.sqrt(1.0 + radius * radius))


def smoothness_invariants(throat: float = THROAT, radius: float = 2.0) -> dict[str, Any]:
    k = wave_number(radius)
    r0 = float(areal_radius(0.0, throat, radius))
    r_left = float(areal_radius(-radius, throat, radius))
    r_right = float(areal_radius(radius, throat, radius))
    rmax = float(np.sqrt(throat * throat + (2.0 * radius / np.pi) ** 2))
    return {
        "a": float(throat),
        "R": float(radius),
        "k": k,
        "k_times_R": k * float(radius),
        "r_min": r0,
        "r_at_origin": r0,
        "r_second_at_origin": float(areal_radius_second(0.0, throat, radius)),
        "expected_r_second_at_origin": 1.0 / float(throat),
        "r_plus_R": r_right,
        "r_minus_R": r_left,
        "r_max": rmax,
        "r_prime_plus_R": float(areal_radius_derivative(radius, throat, radius)),
        "r_prime_minus_R": float(areal_radius_derivative(-radius, throat, radius)),
        "r_second_plus_R": float(areal_radius_second(radius, throat, radius)),
        "r_second_minus_R": float(areal_radius_second(-radius, throat, radius)),
        "r_third_plus_R": float(areal_radius_third(radius, throat, radius)),
        "r_third_minus_R": float(areal_radius_third(-radius, throat, radius)),
        "periodic_C_infinity": True,
        "seamed_sqrt_one_plus_rho_squared_r_prime_jump": seamed_radius_derivative_jump(radius),
        "smooth_r_prime_jump": 0.0,
    }


def eight_pi_G_stress(x, throat: float = THROAT, radius: float = 2.0) -> dict[str, Any]:
    r = areal_radius(x, throat, radius)
    rp = areal_radius_derivative(x, throat, radius)
    rpp = areal_radius_second(x, throat, radius)
    rho = (1.0 - rp * rp - 2.0 * r * rpp) / (r * r)
    p_radial = (rp * rp - 1.0) / (r * r)
    p_tangential = rpp / r
    return {
        "r": r,
        "r_prime": rp,
        "r_second": rpp,
        "eight_pi_G_rho": rho,
        "eight_pi_G_p_r": p_radial,
        "eight_pi_G_p_t": p_tangential,
        "eight_pi_G_rho_plus_p_r": rho + p_radial,
    }


def riemann_convention() -> dict[str, str]:
    return {
        "signature": "(+,-,-,-)",
        "metric": "ds^2 = dt^2 - dx^2 - r(x)^2 dOmega_2^2",
        "connection": "Gamma^l_{ij} = (1/2) g^{ls}(d_i g_{js} + d_j g_{is} - d_s g_{ij})",
        "riemann": "R^r_{sij} = d_i Gamma^r_{js} - d_j Gamma^r_{is} + Gamma^r_{il} Gamma^l_{js} - Gamma^r_{jl} Gamma^l_{is}",
        "ricci": "R_{sj} = R^i_{sij} (contract first and third index)",
        "einstein": "G_{ij} = R_{ij} - (1/2) R g_{ij}",
        "static_observer": "u^mu = (1,0,0,0),  8 pi G rho = G_tt = G^t_t",
        "radial_pressure": "8 pi G p_r = G_xx = -G^x_x",
        "tangential_pressure": "8 pi G p_t = G_theta_theta / r^2 = -G^theta_theta",
    }


def _christoffel(metric, inverse, coords):
    dim = len(coords)
    connection = [[[0] * dim for _ in range(dim)] for _ in range(dim)]
    for lam in range(dim):
        for i in range(dim):
            for j in range(dim):
                total = 0
                for sigma in range(dim):
                    total += inverse[lam, sigma] * (
                        sp.diff(metric[j, sigma], coords[i])
                        + sp.diff(metric[i, sigma], coords[j])
                        - sp.diff(metric[i, j], coords[sigma])
                    )
                connection[lam][i][j] = sp.simplify(total / 2)
    return connection


def _ricci(connection, coords):
    dim = len(coords)
    ricci = sp.zeros(dim)
    for i in range(dim):
        for j in range(dim):
            value = 0
            for k in range(dim):
                value += sp.diff(connection[k][j][i], coords[k])
                value -= sp.diff(connection[k][k][i], coords[j])
                for m in range(dim):
                    value += connection[k][k][m] * connection[m][j][i]
                    value -= connection[k][j][m] * connection[m][k][i]
            ricci[i, j] = sp.simplify(value)
    return ricci


def derive_einstein_identities() -> dict[str, Any]:
    """Exact curvature, conservation, and lapse identities for generic r(x)."""
    t, x, theta, phi = sp.symbols("t x theta phi", real=True)
    radius = sp.Function("r", positive=True)
    lapse = sp.Function("N", positive=True)
    rp = sp.diff(radius(x), x)
    rpp = sp.diff(radius(x), x, 2)
    expected_rho = (1 - rp**2 - 2 * radius(x) * rpp) / radius(x) ** 2
    expected_p_r = (rp**2 - 1) / radius(x) ** 2
    expected_p_t = rpp / radius(x)
    expected_spatial_scalar = 2 * expected_rho

    spatial_coords = [x, theta, phi]
    spatial_metric = sp.diag(1, radius(x) ** 2, radius(x) ** 2 * sp.sin(theta) ** 2)
    spatial_inverse = spatial_metric.inv()
    spatial_gamma = _christoffel(spatial_metric, spatial_inverse, spatial_coords)
    spatial_ricci = _ricci(spatial_gamma, spatial_coords)
    spatial_scalar = sp.simplify(
        sum(
            spatial_inverse[i, j] * spatial_ricci[i, j]
            for i in range(3)
            for j in range(3)
        )
    )

    spacetime_coords = [t, x, theta, phi]
    metric = sp.diag(1, -1, -radius(x) ** 2, -radius(x) ** 2 * sp.sin(theta) ** 2)
    inverse = metric.inv()
    gamma = _christoffel(metric, inverse, spacetime_coords)
    ricci = _ricci(gamma, spacetime_coords)
    scalar = sp.simplify(
        sum(inverse[i, j] * ricci[i, j] for i in range(4) for j in range(4))
    )
    einstein = sp.zeros(4)
    for i in range(4):
        for j in range(4):
            einstein[i, j] = sp.simplify(ricci[i, j] - sp.Rational(1, 2) * scalar * metric[i, j])
    mixed = sp.zeros(4)
    for i in range(4):
        for j in range(4):
            mixed[i, j] = sp.simplify(sum(inverse[i, k] * einstein[k, j] for k in range(4)))

    conservation = sp.simplify(sp.diff(expected_p_r, x) + 2 * rp / radius(x) * (expected_p_r - expected_p_t))
    null_local = sp.simplify(expected_rho + expected_p_r + 2 * rpp / radius(x))

    lapse_metric = sp.diag(
        lapse(x) ** 2, -1, -radius(x) ** 2, -radius(x) ** 2 * sp.sin(theta) ** 2
    )
    lapse_inverse = lapse_metric.inv()
    lapse_gamma = _christoffel(lapse_metric, lapse_inverse, spacetime_coords)
    lapse_ricci = _ricci(lapse_gamma, spacetime_coords)
    lapse_scalar = sp.simplify(
        sum(lapse_inverse[i, j] * lapse_ricci[i, j] for i in range(4) for j in range(4))
    )
    g_tt = lapse_ricci[0, 0] - sp.Rational(1, 2) * lapse_scalar * lapse_metric[0, 0]
    normal_projection = sp.simplify(g_tt / lapse(x) ** 2)
    lagrangian = sp.simplify(lapse(x) * radius(x) ** 2 * lapse_scalar)
    n_sym, np_sym, npp_sym = sp.symbols("N Np Npp", real=True)
    r_sym, rp_sym, rpp_sym = sp.symbols("r rp rpp", real=True)
    substituted = lagrangian.subs(
        {
            lapse(x): n_sym,
            sp.diff(lapse(x), x): np_sym,
            sp.diff(lapse(x), x, 2): npp_sym,
            radius(x): r_sym,
            rp: rp_sym,
            rpp: rpp_sym,
        }
    )
    d_l_d_n = sp.diff(substituted, n_sym)
    d_l_d_np = sp.diff(substituted, np_sym)
    d_l_d_npp = sp.diff(substituted, npp_sym)
    as_fields = {
        n_sym: lapse(x),
        np_sym: sp.diff(lapse(x), x),
        npp_sym: sp.diff(lapse(x), x, 2),
        r_sym: radius(x),
        rp_sym: rp,
        rpp_sym: rpp,
    }
    euler = sp.simplify(
        d_l_d_n.subs(as_fields)
        - sp.diff(d_l_d_np.subs(as_fields), x)
        + sp.diff(d_l_d_npp.subs(as_fields), x, 2)
    )
    # EL[N r^2 ^4R] = - r^2 ^3R = -2 r^2 G_nn, independent of the lapse.
    euler_identity = sp.simplify(euler + radius(x) ** 2 * spatial_scalar)

    zero = lambda expr: str(sp.simplify(expr))
    return {
        "spatial_Gamma_x_theta_theta": str(spatial_gamma[0][1][1]),
        "spatial_scalar_curvature": str(spatial_scalar),
        "spatial_scalar_minus_twice_eight_pi_G_rho": zero(spatial_scalar - expected_spatial_scalar),
        "four_scalar_curvature": str(scalar),
        "G_tt": str(einstein[0, 0]),
        "G_xx": str(einstein[1, 1]),
        "G_theta_theta": str(einstein[2, 2]),
        "G_tt_minus_eight_pi_G_rho": zero(einstein[0, 0] - expected_rho),
        "G_xx_minus_eight_pi_G_p_r": zero(einstein[1, 1] - expected_p_r),
        "G_theta_theta_minus_r_squared_times_eight_pi_G_p_t": zero(
            einstein[2, 2] - radius(x) ** 2 * expected_p_t
        ),
        "G_up_t_t_minus_eight_pi_G_rho": zero(mixed[0, 0] - expected_rho),
        "G_up_x_x_plus_eight_pi_G_p_r": zero(mixed[1, 1] + expected_p_r),
        "G_up_theta_theta_plus_eight_pi_G_p_t": zero(mixed[2, 2] + expected_p_t),
        "conservation_p_r_prime_plus_two_r_prime_over_r_times_p_r_minus_p_t": zero(conservation),
        "local_null_identity_rho_plus_p_r_plus_two_r_second_over_r": zero(null_local),
        "lapse_normal_projection_minus_eight_pi_G_rho": zero(normal_projection - expected_rho),
        "lapse_normal_projection_minus_half_spatial_scalar": zero(
            normal_projection - spatial_scalar / 2
        ),
        "lapse_euler_plus_r_squared_spatial_scalar": zero(euler_identity),
        "ADM_hamiltonian_eight_pi_G_rho_equals_half_spatial_scalar": True,
        "independent_checks": (
            "3D Christoffel/Ricci scalar",
            "4D Christoffel/Ricci/Einstein tensor",
            "lapse-normal projection G_mu nu n^mu n^nu",
            "Euler-Lagrange of int N r^2 ^4R",
        ),
    }


def finite_difference_curvature_check(
    throat: float = THROAT, radius: float = 2.0, location: float = 0.35, step: float = 1e-4
) -> dict[str, Any]:
    """Metric-derivative Ricci check that does not reuse the closed r', r'' formulas."""
    nodes = location + step * np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    g_theta = np.square(areal_radius(nodes, throat, radius))
    dg = (g_theta[0] - 8.0 * g_theta[1] + 8.0 * g_theta[3] - g_theta[4]) / (12.0 * step)
    d2g = (
        -g_theta[0] + 16.0 * g_theta[1] - 30.0 * g_theta[2] + 16.0 * g_theta[3] - g_theta[4]
    ) / (12.0 * step * step)
    g = float(g_theta[2])
    r = float(np.sqrt(g))
    gamma_x_thth = -0.5 * dg
    ric_xx = -d2g / g + (dg * dg) / (2.0 * g * g)
    # R_theta_theta = 1 - r'^2 - r r''. With g=r^2 those first-derivative pieces cancel,
    # leaving 1 - (1/2) g''.
    ric_thth = 1.0 - 0.5 * d2g
    spatial_scalar = ric_xx + 2.0 * ric_thth / g
    eight_rho = spatial_scalar / 2.0
    eight_p_r = (dg * dg) / (4.0 * g * g) - 1.0 / g
    eight_p_t = d2g / (2.0 * g) - (dg * dg) / (4.0 * g * g)
    analytic = eight_pi_G_stress(location, throat, radius)
    return {
        "location": float(location),
        "step": float(step),
        "r_from_metric": r,
        "Ric_xx": float(ric_xx),
        "Ric_theta_theta": float(ric_thth),
        "spatial_scalar": float(spatial_scalar),
        "eight_pi_G_rho": float(eight_rho),
        "eight_pi_G_p_r": float(eight_p_r),
        "eight_pi_G_p_t": float(eight_p_t),
        "rho_residual": float(eight_rho - analytic["eight_pi_G_rho"]),
        "p_r_residual": float(eight_p_r - analytic["eight_pi_G_p_r"]),
        "p_t_residual": float(eight_p_t - analytic["eight_pi_G_p_t"]),
        "Gamma_x_theta_theta_plus_r_r_prime": float(
            gamma_x_thth + r * areal_radius_derivative(location, throat, radius)
        ),
    }


def conservation_residual_grid(
    throat: float = THROAT, radius: float = 2.0, samples: int = 4096
) -> dict[str, Any]:
    x = np.linspace(-radius, radius, samples, endpoint=False)
    stress = eight_pi_G_stress(x, throat, radius)
    spacing = float(x[1] - x[0])
    p_r = np.asarray(stress["eight_pi_G_p_r"], dtype=float)
    p_t = np.asarray(stress["eight_pi_G_p_t"], dtype=float)
    rp = np.asarray(stress["r_prime"], dtype=float)
    rpp = np.asarray(stress["r_second"], dtype=float)
    r = np.asarray(stress["r"], dtype=float)
    exact_p_r_prime = 2.0 * rp * (r * rpp - rp * rp + 1.0) / r**3
    exact_residual = exact_p_r_prime + 2.0 * rp / r * (p_r - p_t)
    rolled = (np.roll(p_r, -1) - np.roll(p_r, 1)) / (2.0 * spacing)
    discrete_residual = rolled + 2.0 * rp / r * (p_r - p_t)
    return {
        "samples": samples,
        "spacing": spacing,
        "max_abs_exact_residual": float(np.max(np.abs(exact_residual))),
        "max_abs_periodic_central_residual": float(np.max(np.abs(discrete_residual))),
        "rms_periodic_central_residual": float(np.sqrt(np.mean(discrete_residual * discrete_residual))),
        "periodic_gradient": True,
    }


def period_null_integrals(throat: float = THROAT, radius: float = 2.0) -> dict[str, Any]:
    def rho_plus_p(x):
        return float(eight_pi_G_stress(x, throat, radius)["eight_pi_G_rho_plus_p_r"])

    def slope(x):
        return float(
            (areal_radius_derivative(x, throat, radius) / areal_radius(x, throat, radius)) ** 2
        )

    null, null_err = quad(rho_plus_p, -radius, radius, epsabs=1e-12, epsrel=1e-12)
    slope_integral, slope_err = quad(slope, -radius, radius, epsabs=1e-12, epsrel=1e-12)
    throat_stress = eight_pi_G_stress(0.0, throat, radius)
    return {
        "eight_pi_G_integral_rho_plus_p_r": float(null),
        "integral_abserr": float(null_err),
        "minus_two_integral_r_prime_over_r_squared": float(-2.0 * slope_integral),
        "slope_abserr": float(slope_err),
        "identity_residual": float(null + 2.0 * slope_integral),
        "strictly_negative": bool(null < 0),
        "throat_eight_pi_G_rho": float(throat_stress["eight_pi_G_rho"]),
        "throat_eight_pi_G_p_r": float(throat_stress["eight_pi_G_p_r"]),
        "throat_eight_pi_G_p_t": float(throat_stress["eight_pi_G_p_t"]),
        "throat_eight_pi_G_rho_plus_p_r": float(throat_stress["eight_pi_G_rho_plus_p_r"]),
        "expected_throat_eight_pi_G_rho_plus_p_r": -2.0 / float(throat * throat),
    }


def cosmological_constant_null_projection() -> dict[str, Any]:
    return {
        "equation_of_state": "p = -rho",
        "eight_pi_G_rho_plus_p": 0.0,
        "can_supply_required_radial_null_stress": False,
        "reason": (
            "A cosmological constant contributes equally to 8 pi G rho and -8 pi G p, "
            "so its radial-null projection vanishes. The smooth profile has "
            "8 pi G int (rho + p_r) dx = -2 int (r'/r)^2 dx < 0 whenever r' is not identically zero."
        ),
    }


def zero_energy_integral(throat: float = THROAT, radius: float = 2.0, kappa: float = KAPPA) -> float:
    k = wave_number(radius)
    return float(2.0 * kappa / (throat * k) * ellipk(-1.0 / (throat * throat * k * k)))


def continuum_discriminant(
    energy: float, throat: float = THROAT, radius: float = 2.0, kappa: float = KAPPA
) -> tuple[float, float]:
    def rhs(_x, state):
        weight = superpotential(_x, throat, radius, kappa)
        generator = np.array([[-weight, energy], [-energy, weight]], dtype=float)
        return (generator @ state.reshape(2, 2)).ravel()

    solution = solve_ivp(
        rhs,
        (-radius, radius),
        np.eye(2).ravel(),
        method="DOP853",
        rtol=ODE_RTOL,
        atol=ODE_ATOL,
        max_step=ODE_MAX_STEP,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    transfer = solution.y[:, -1].reshape(2, 2)
    return float(np.trace(transfer)), float(np.linalg.det(transfer))


def continuum_first_band_edge(throat: float = THROAT, radius: float = 2.0) -> dict[str, Any]:
    previous = 0.0
    for energy in np.linspace(0.0, 1.0, 129)[1:]:
        trace, _ = continuum_discriminant(energy, throat, radius)
        if trace <= 2.0:
            root = brentq(
                lambda value: continuum_discriminant(value, throat, radius)[0] - 2.0,
                previous,
                energy,
                xtol=1e-12,
            )
            final, determinant = continuum_discriminant(root, throat, radius)
            return {
                "band_edge": float(root),
                "discriminant_residual": float(final - 2.0),
                "monodromy_determinant": float(determinant),
                "bracket": [float(previous), float(energy)],
            }
        previous = float(energy)
    raise RuntimeError("first band edge was not bracketed; no gap value inferred")


def smooth_motif(
    throat: float, radius: float, intervals: int, kappa: float = KAPPA
) -> dict[str, Any]:
    """Node/edge staggered Dirac motif, same stencil family as nsc_geometric_chain."""
    if intervals < 4:
        raise ValueError("at least four intervals required")
    spacing = 2.0 * radius / intervals
    edges = -radius + (np.arange(intervals) + 0.5) * spacing
    weight = superpotential(edges, throat, radius, kappa)
    if np.max(spacing * weight) >= 2.0:
        raise ValueError("refine grid so h*w<2; retain the first-order transfer branch")
    left = -1.0 / spacing + weight / 2.0
    right = 1.0 / spacing + weight / 2.0
    forward = np.diag(left) + np.diag(right[:-1], 1)
    local = np.block(
        [[np.zeros_like(forward), forward.T], [forward, np.zeros_like(forward)]]
    )
    link = np.zeros_like(local)
    link[-1, 0] = right[-1]
    return {
        "H": local,
        "B": link,
        "h": float(spacing),
        "zero_energy_transfer": float(np.prod(-left / right)),
        "continuum_zero_transfer": float(np.exp(-zero_energy_integral(throat, radius, kappa))),
    }


def bloch_gap(
    throat: float, radius: float, intervals: int, phases: int = PHASES
) -> dict[str, Any]:
    motif = smooth_motif(throat, radius, intervals)
    minimum = np.inf
    location = None
    for phase in np.linspace(0.0, np.pi, phases):
        matrix = motif["H"] + np.exp(1j * phase) * motif["B"] + np.exp(-1j * phase) * motif["B"].T
        gap = float(np.min(np.abs(eigh(matrix, eigvals_only=True))))
        if gap < minimum:
            minimum, location = gap, float(phase)
    return {
        "intervals": int(intervals),
        "spacing": motif["h"],
        "phase_samples": int(phases),
        "minimum_sampled_bloch_gap": float(minimum),
        "minimum_phase": location,
        "zero_transfer": motif["zero_energy_transfer"],
        "continuum_zero_transfer": motif["continuum_zero_transfer"],
    }


def source_hashes(paths=OWNED_SOURCE_PATHS) -> dict[str, str]:
    return {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in paths}


def preserved_published_hashes() -> dict[str, Any]:
    inventory: dict[str, Any] = {
        "commits": list(PUBLISHED_COMMITS),
        "result_source_hashes": {},
        "verified_file_hashes": {},
    }
    for relative, field in FROZEN_RESULT_HASH_FIELDS:
        record = json.loads((ROOT / relative).read_text())
        expected = record[field]
        inventory["result_source_hashes"][relative] = expected
        for path, digest in expected.items():
            observed = hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
            if observed != digest:
                raise RuntimeError(f"published source bytes changed: {path}")
            inventory["verified_file_hashes"][path] = observed
    for relative, field, path in FROZEN_SINGLE_HASH:
        record = json.loads((ROOT / relative).read_text())
        digest = record[field]
        observed = hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
        if observed != digest:
            raise RuntimeError(f"published source bytes changed: {path}")
        inventory["result_source_hashes"][relative] = {path: digest}
        inventory["verified_file_hashes"][path] = observed
    chain = hashlib.sha256((ROOT / "src/recursive_horizons/nsc_geometric_chain.py").read_bytes()).hexdigest()
    inventory["read_only_imported_module"] = {
        "path": "src/recursive_horizons/nsc_geometric_chain.py",
        "sha256": chain,
        "use": "seamed continuum band-edge comparison only",
    }
    return inventory


def jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, (np.floating, float)):
        number = float(value)
        if not np.isfinite(number):
            raise RuntimeError("nonfinite float in record")
        return number
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, str):
        return value
    if isinstance(value, bool) or value is None:
        return value
    raise TypeError(f"unserializable {type(value)!r}")


def gap_family(throat: float, radius: float) -> dict[str, Any]:
    continuum = continuum_first_band_edge(throat, radius)
    lattice = [bloch_gap(throat, radius, count) for count in INTERVALS]
    errors = [
        abs(row["minimum_sampled_bloch_gap"] - continuum["band_edge"]) for row in lattice
    ]
    orders = [float(np.log2(errors[index] / errors[index + 1])) for index in range(2)]
    zero_trace, zero_det = continuum_discriminant(0.0, throat, radius)
    integral = zero_energy_integral(throat, radius)
    analytic_trace = float(2.0 * np.cosh(integral))
    if abs(zero_trace - analytic_trace) >= 1e-8 or abs(zero_det - 1.0) >= 1e-8:
        raise RuntimeError("zero-energy monodromy failed the elliptic-integral identity")
    if errors[-1] >= 2e-5 or errors[0] <= errors[1] or errors[1] <= errors[2]:
        raise RuntimeError("Bloch gap is not converging onto the continuum edge")
    if min(orders) <= 1.5 or max(orders) >= 2.4:
        raise RuntimeError("observed lattice order left the second-order window")
    seamed = seamed_continuum_first_band_edge(radius)
    null = period_null_integrals(throat, radius)
    smooth = smoothness_invariants(throat, radius)
    conservation = conservation_residual_grid(throat, radius)
    if abs(null["identity_residual"]) >= 1e-10 or not null["strictly_negative"]:
        raise RuntimeError("period-integrated radial-null identity failed")
    if conservation["max_abs_exact_residual"] >= 1e-10:
        raise RuntimeError("pointwise conservation identity failed on the profile")
    if conservation["max_abs_periodic_central_residual"] >= 5e-3:
        raise RuntimeError("periodic central conservation residual is too large")
    if abs(smooth["r_prime_plus_R"]) >= 1e-12 or abs(smooth["r_third_plus_R"]) >= 1e-10:
        raise RuntimeError("periodic endpoint derivatives are not a smooth seam")
    return {
        "radius": float(radius),
        "throat": float(throat),
        "smoothness": smooth,
        "continuum": continuum,
        "lattice": lattice,
        "absolute_errors": errors,
        "observed_orders": orders,
        "zero_discriminant": zero_trace,
        "analytic_zero_discriminant": analytic_trace,
        "zero_energy_integral": integral,
        "zero_monodromy_determinant": zero_det,
        "seamed_continuum_band_edge": seamed["band_edge"],
        "smooth_minus_seamed_band_edge": float(continuum["band_edge"] - seamed["band_edge"]),
        "null_stress": null,
        "conservation": conservation,
        "gap_versus_null_stress": {
            "continuum_band_edge": continuum["band_edge"],
            "eight_pi_G_period_integral_rho_plus_p_r": null["eight_pi_G_integral_rho_plus_p_r"],
            "throat_eight_pi_G_rho_plus_p_r": null["throat_eight_pi_G_rho_plus_p_r"],
            "same_profile": True,
            "fitted_coupling": False,
            "added_phantom_field": False,
        },
    }


def build_record() -> dict[str, Any]:
    identities = derive_einstein_identities()
    residual_keys = (
        "spatial_scalar_minus_twice_eight_pi_G_rho",
        "G_tt_minus_eight_pi_G_rho",
        "G_xx_minus_eight_pi_G_p_r",
        "G_theta_theta_minus_r_squared_times_eight_pi_G_p_t",
        "G_up_t_t_minus_eight_pi_G_rho",
        "G_up_x_x_plus_eight_pi_G_p_r",
        "G_up_theta_theta_plus_eight_pi_G_p_t",
        "conservation_p_r_prime_plus_two_r_prime_over_r_times_p_r_minus_p_t",
        "local_null_identity_rho_plus_p_r_plus_two_r_second_over_r",
        "lapse_normal_projection_minus_eight_pi_G_rho",
        "lapse_normal_projection_minus_half_spatial_scalar",
        "lapse_euler_plus_r_squared_spatial_scalar",
    )
    for key in residual_keys:
        if identities[key] != "0":
            raise RuntimeError(f"curvature identity failed: {key} = {identities[key]}")
    families = [gap_family(THROAT, radius) for radius in RADII]
    metric_check = finite_difference_curvature_check()
    if max(abs(metric_check[name]) for name in ("rho_residual", "p_r_residual", "p_t_residual")) >= 1e-8:
        raise RuntimeError("finite-difference Einstein residual is too large")
    free = smooth_motif(THROAT, 2.0, 32, kappa=0.0)
    free_gap = float(np.min(np.abs(np.linalg.eigvalsh(free["H"] + free["B"] + free["B"].T))))
    if free_gap >= 1e-12:
        raise RuntimeError("removing the angular potential did not close the zero-phase gap")
    return {
        "schema": SCHEMA,
        "artifact_id": ARTIFACT,
        "classification": (
            "seam_free_periodic_areal_radius_opens_a_spatial_band_gap_and_requires_"
            "NEC_violating_effective_Einstein_stress_that_a_cosmological_constant_cannot_supply"
        ),
        "source_hashes": source_hashes(),
        "preserved_published_bytes": preserved_published_hashes(),
        "declared_geometry": {
            "formula": "r(x)=sqrt(a^2+[sin(k x)/k]^2), k=pi/(2R), x in [-R,R]",
            "a": THROAT,
            "R_values": list(RADII),
            "a_is_input": True,
            "R_is_input": True,
            "independent_mass_inserted": False,
            "added_phantom_field": False,
            "fitted_coupling": False,
            "derivative_seams": "removed; r is C^infinity 2R-periodic",
            "stationary_GR": False,
            "stress_status": (
                "required_effective_Einstein_stress_of_the_imposed_4D_ultrastatic_metric; "
                "not_derived_matter_stress_of_S_one_or_full_spectral_5D_theory"
            ),
        },
        "riemann_convention": riemann_convention(),
        "curvature_identities": identities,
        "finite_difference_metric_check": metric_check,
        "cosmological_constant": cosmological_constant_null_projection(),
        "zero_gap_exclusion": {
            "continuum_transfer_eigenvalues": "exp(plus_or_minus I), I = kappa * int dx/r = 2 kappa/(a k) K(-1/(a^2 k^2))",
            "argument": (
                "I>0 so neither multiplier has modulus 1; zero is outside every real Bloch fiber; "
                "compact phase and continuity give an open spectral gap"
            ),
            "discrete_transfer": "product[(1-h*w/2)/(1+h*w/2)] in (0,1)",
            "scope": "equal-scale periodic spatial family on the smooth profile only",
        },
        "gap_families": families,
        "removed_angular_potential_control": {
            "kappa": 0.0,
            "zero_phase_gap": free_gap,
            "interpretation": "formal control with w removed; kappa=0 is not a spinor-sphere eigenvalue",
        },
        "proof_scope": {
            "spatial_dirac_gap": True,
            "required_4D_ultrastatic_Einstein_stress": True,
            "stationary_Einstein_solution": False,
            "S_one_or_spectral_matter_stress": False,
            "five_dimensional_or_full_heat_kernel_stress": False,
            "particle_species_identification": False,
            "derived_a_or_R": False,
            "new_to_world_mechanism": False,
        },
        "next_equation": {
            "statement": (
                "The geometric residual G_mu nu[r] must be identified with the variational "
                "stress of the same one-action / spectral-boundary functional at frozen Theta, "
                "with no extra phantom field and no fitted coupling."
            ),
            "einstein_matching": "G_mu nu[r] = 8 pi G T^{rest}_mu nu, with T_rest from Gamma_one minus its already included Einstein term",
            "radial_null_constraint": (
                "int_{-R}^{R} ell^mu ell^nu (G_mu nu[r] - 8 pi G T^{rest}_mu nu) dx = 0, "
                "equivalently 8 pi G int (rho + p_r)_one dx = -2 int (r'/r)^2 dx"
            ),
            "owner": (
                "quantum/boundary stress of S_one on this, or a dynamically selected, confining profile; "
                "the present record only computes the geometric right-hand side"
            ),
        },
        "nonclaims": {
            "stationary_geometry_or_a_R_derived": False,
            "periodic_spatial_array_is_a_Lorentzian_nested_cosmology": False,
            "derived_S_one_or_spectral_matter_stress": False,
            "cosmological_constant_accounts_for_the_null_projection": False,
            "particle_species_or_observed_Lambda_predicted": False,
            "physical_metric_stability_proved": False,
            "new_to_world_mechanism_established": False,
        },
        "terminal": True,
    }
