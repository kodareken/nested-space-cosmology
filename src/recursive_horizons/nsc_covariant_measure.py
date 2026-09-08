"""Covariant-measure obstructions and a full-angular vacuum-response control.

Natural units hbar=c=1. The cylinder is the four-dimensional ultrastatic
product R_t x S1_L x S2_a, with one complex four-component massless Dirac
field in its static ground state. It is not the trapped five-dimensional
carrier. The two allowed eta values are spatial spin structures, not fitted
continuous boundary phases. No gravitational counterterm is added here.
"""
from __future__ import annotations

from functools import lru_cache
from math import pi

import numpy as np
from scipy.integrate import quad, quad_vec
from scipy.special import expit, kv
import sympy as sp

from .nsc_regulated import OperatorConventions, RegulatedOperator


@lru_cache(maxsize=1)
def product_geometry_identities():
    """Contract the coordinate Riemann tensor, then vary the C^2 functional."""
    coords = sp.symbols("tau x theta phi", real=True)
    a, length, time, alpha, volume = sp.symbols("a L T alpha V", positive=True)
    metric = sp.diag(1, (length / (2 * sp.pi))**2,
                     a**2, a**2 * sp.sin(coords[2])**2)
    inverse = metric.inv()
    dim = 4
    christoffel = [[[sp.simplify(sum(
        inverse[i, m] * (sp.diff(metric[m, j], coords[k])
                         + sp.diff(metric[m, k], coords[j])
                         - sp.diff(metric[j, k], coords[m]))
        for m in range(dim)) / 2)
        for k in range(dim)] for j in range(dim)] for i in range(dim)]
    riemann = [[[[sp.simplify(
        sp.diff(christoffel[i][j][l], coords[k])
        - sp.diff(christoffel[i][j][k], coords[l])
        + sum(christoffel[i][m][k] * christoffel[m][j][l]
              - christoffel[i][m][l] * christoffel[m][j][k]
              for m in range(dim)))
        for l in range(dim)] for k in range(dim)]
        for j in range(dim)] for i in range(dim)]
    ricci = sp.Matrix(dim, dim, lambda j, l: sp.simplify(
        sum(riemann[i][j][i][l] for i in range(dim))))
    scalar = sp.simplify(sum(inverse[i, j] * ricci[i, j]
                             for i in range(dim) for j in range(dim)))
    ricci2 = sp.simplify(sum(inverse[i, i] * inverse[j, j] * ricci[i, j]**2
                            for i in range(dim) for j in range(dim)))
    riemann2 = sp.simplify(sum(
        inverse[i, i] * inverse[j, j] * inverse[k, k] * inverse[l, l]
        * (metric[i, i] * riemann[i][j][k][l])**2
        for i in range(dim) for j in range(dim)
        for k in range(dim) for l in range(dim)))
    weyl2 = sp.simplify(riemann2 - 2 * ricci2 + scalar**2 / 3)
    euler = sp.simplify(riemann2 - 4 * ricci2 + scalar**2)
    energy = sp.simplify(alpha * 4 * sp.pi * a**2 * length * weyl2)
    action = time * energy
    scale_residual = sp.simplify(sum(v * sp.diff(action, v)
                                     for v in (time, length, a)))
    fixed_volume = sp.simplify(energy.subs(a**2, volume / (4 * sp.pi * length)))
    density = sp.simplify(energy / (4 * sp.pi * a**2 * length))
    radial_pressure = sp.simplify(-sp.diff(energy, length) / (4 * sp.pi * a**2))
    transverse_pressure = sp.simplify(-sp.diff(energy, a) / (8 * sp.pi * a * length))
    # Pointwise Weyl weights in d=4: sqrt(g) -> exp(4 sigma),
    # C_abcd C^abcd -> exp(-4 sigma). This is not just rigid scaling.
    sigma = sp.symbols("sigma", real=True)
    local_weyl_residual = sp.diff(sp.exp(4 * sigma) * sp.exp(-4 * sigma), sigma)
    return {key: str(value) for key, value in {
        "scalar_curvature": scalar, "ricci_squared": ricci2,
        "riemann_squared": riemann2, "weyl_squared": weyl2,
        "euler_density": euler, "delta_energy": energy,
        "spacing_derivative": sp.diff(energy, length),
        "sphere_radius_derivative": sp.diff(energy, a),
        "rigid_four_dimensional_weyl_residual": scale_residual,
        "local_weyl_weight_residual": local_weyl_residual,
        "fixed_volume_delta_energy": fixed_volume,
        "fixed_volume_spacing_derivative": sp.diff(fixed_volume, length),
        "local_C_squared_energy_density": density,
        "local_C_squared_radial_pressure": radial_pressure,
        "local_C_squared_transverse_pressure": transverse_pressure,
        "local_C_squared_radial_null_projection": sp.simplify(density + radial_pressure),
        "local_C_squared_stress_trace": sp.simplify(-density + radial_pressure + 2 * transverse_pressure),
    }.items()}


def finite_weyl_average(extents=(5., 10., 20.)):
    """Integrate the unchanged finite prescription on negative Weyl orbits.

    For any nonzero finite D, E1(e^(-2 phi) D^2/Lambda^2) -> 0 as
    phi -> -infinity. The flat group-average integrand tends to a positive
    constant; the infinite integral diverges. Quadrature checks this example.
    """
    op = RegulatedOperator(np.diag([1., -1.]),
                           OperatorConventions(cutoff=2., normalization=2.,
                                               spatial_dimension=3))
    identity = np.eye(2)

    def integrand(phi):
        return float(np.exp(-op.scale_orbit(identity, phi).action()))

    rows = []
    for extent in extents:
        if not np.isfinite(extent) or extent <= 0 or extent > 100:
            raise ValueError("use a finite orbit extent in (0,100]")
        value, error = quad(integrand, -extent, 0., epsabs=1e-12, epsrel=1e-12)
        rows.append({"extent": float(extent), "negative_orbit_integral": value,
                     "quadrature_error_estimate": error,
                     "left_endpoint_integrand": integrand(-extent)})
    return {"dimension": 2, "dirac_eigenvalues": [1., -1.],
            "cutoff": 2., "normalization": 2.,
            "asymptotic_integrand": float(np.exp(-np.euler_gamma)), "integrals": rows}


def _parameters(length, radius, angular_max, winding_max):
    for value in (length, radius):
        if isinstance(value, bool) or not np.isfinite(value) or value <= 0:
            raise ValueError("length and radius must be finite and positive")
    for value in (angular_max, winding_max):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError("angular and winding truncations must be positive integers")
    return (np.arange(1, angular_max + 1, dtype=float)[:, None],
            np.arange(1, winding_max + 1, dtype=float)[None, :])


def cylinder_energy(length, radius=1., eta=0., angular_max=64, winding_max=64):
    """Return (E, dE/dL) relative to the decompactified cylinder.

    A sphere eigenspace has eigenvalue +/-k/a with degeneracy 2k per sign.
    Tensoring with the two-component (t,x) Dirac factor yields H eigenvalues
    +/-sqrt(p^2+k^2/a^2), each with degeneracy 4k. Hence -Tr|H|/2 is
    -4k sqrt(...) per spatial momentum, not half or twice this expression.
    The zero-winding term is subtracted in this relative observable.
    """
    k, n = _parameters(length, radius, angular_max, winding_max)
    if eta not in (0., .5):
        raise ValueError("eta must be the periodic or antiperiodic spatial spin structure")
    weights = np.ones_like(n) if eta == 0 else (-1.)**n
    z = n * k * length / radius
    energy = 8 / (pi * radius) * np.sum(k**2 * weights * kv(1, z) / n)
    derivative = -4 / (pi * radius**2) * np.sum(k**3 * weights * (kv(0, z) + kv(2, z)))
    return float(energy), float(derivative)


def spin_difference(length, radius=1., angular_max=64, winding_max=64):
    """Counterterm-independent E_periodic-E_antiperiodic on the same metric."""
    k, n = _parameters(length, radius, angular_max, winding_max)
    weights = 1 - (-1.)**n
    z = n * k * length / radius
    energy = 8 / (pi * radius) * np.sum(k**2 * weights * kv(1, z) / n)
    derivative = -4 / (pi * radius**2) * np.sum(k**3 * weights * (kv(0, z) + kv(2, z)))
    return float(energy), float(derivative)


def cylinder_stress(length, radius=1., eta=0., angular_max=64, winding_max=64):
    """Finite compactification stress, relative to the infinite cylinder.

    For the massless field on this constant-radius product, the decompactified
    ground state and all local covariant counterterms are Lorentz invariant in
    the (t,x) factor and hence have rho+p_x=0. This radial-null projection is
    consequently also the unambiguous absolute projection in the compact
    ground state. Individual rho and pressures here are relative quantities.
    No varying-radius or trapped-domain stress is inferred from this control.
    """
    energy, length_derivative = cylinder_energy(length, radius, eta, angular_max, winding_max)
    k, n = _parameters(length, radius, angular_max, winding_max)
    weights = np.ones_like(n) if eta == 0 else (-1.)**n
    z = n * k * length / radius
    radius_derivative = float(8 / (pi * radius**2) * np.sum(
        k**2 * weights / n * (-kv(1, z) + z * (kv(0, z) + kv(2, z)) / 2)))
    area = 4 * pi * radius**2
    density = energy / (area * length)
    radial = -length_derivative / area
    transverse = -radius_derivative / (8 * pi * radius * length)
    return {"energy": energy, "length_derivative": length_derivative,
            "radius_derivative": radius_derivative,
            "relative_energy_density": density, "relative_radial_pressure": radial,
            "relative_transverse_pressure": transverse,
            "radial_null_stress": density + radial,
            "relative_stress_trace": -density + radial + 2 * transverse,
            "homogeneity_residual": energy + length * length_derivative + radius * radius_derivative}


def antiperiodic_fermi_integral(length, radius=1., angular_max=64):
    """Independent AP energy and two metric derivatives from a Fermi log.

    E_AP=-(8/pi) sum_k k int_0^infinity dp log(1+exp(-L omega)),
    omega=sqrt(p^2+k^2/a^2). The integrand gives E<0 and E_L>0
    without an alternating winding sum. E_a is differentiated under this
    integral independently of Bessel differentiation or a scaling identity.
    """
    k, _ = _parameters(length, radius, angular_max, 1)
    k = k[:, 0]

    def integrand(momentum):
        omega = np.sqrt(momentum**2 + (k / radius)**2)
        occupation = expit(-length * omega)
        return np.array([
            -8 / pi * np.sum(k * np.logaddexp(0., -length * omega)),
            8 / pi * np.sum(k * omega * occupation),
            -8 * length / (pi * radius**3) * np.sum(k**3 / omega * occupation),
        ])

    values, error = quad_vec(integrand, 0., np.inf, epsabs=2e-12, epsrel=2e-12)
    return {"energy": float(values[0]), "length_derivative": float(values[1]),
            "radius_derivative": float(values[2]), "quadrature_error_estimate_only": float(error)}


def winding_theta_sums(b, eta=0.):
    """Stable (sum cos(2pi n eta)e^(-b n²), sum n² cos(...)e^(-b n²)).

    Direct winding evaluation is used only for b>=1. For b<1, differentiated
    Poisson identities avoid catastrophic AP cancellation. Twelve direct
    terms or eight dual terms leave exponentially tiny tails in their chosen
    branches (at most exp(-169) or exp(-pi²*8²) before polynomial factors).
    Very small AP second sums can underflow to zero; no numerical sign is
    inferred in that regime. The exact product identity fixes their sign.
    """
    if not np.isfinite(b) or b <= 0 or eta not in (0., .5):
        raise ValueError("use finite b>0 and periodic/antiperiodic eta")
    if b >= 1:
        n = np.arange(1, 13, dtype=float)
        weights = np.ones_like(n) if eta == 0 else (-1.)**n
        terms = weights * np.exp(-b * n**2)
        return float(terms.sum()), float(np.sum(n**2 * terms))
    if eta == .5:
        m = np.arange(8, dtype=float) + .5
        c = pi**2 * m**2
        terms = np.exp(-c / b)
        theta = 2 * np.sqrt(pi / b) * terms.sum()
        second = -np.sqrt(pi) / b**1.5 * np.sum(terms * (c / b - .5))
    else:
        m = np.arange(1, 9, dtype=float)
        c = pi**2 * m**2
        terms = np.exp(-c / b)
        theta = np.sqrt(pi / b) * (1 + 2 * terms.sum())
        second = (np.sqrt(pi) / (4 * b**1.5) * (1 + 2 * terms.sum())
                  - np.sqrt(pi) / b**2.5 * np.sum(c * terms))
    return float((theta - 1) / 2), float(second)


def finite_cutoff_cylinder_stress(length, radius=1., cutoff=1., eta=0., angular_max=64):
    """Proper-time compactification response at fixed dimensional Lambda.

    This uses the same lower proper-time limit Lambda^-2 as the existing
    regulated modulus, after frequency integration and zero-winding
    subtraction. It does not assemble the unresolved full compensator.
    The null direction is the COMPACT axial x direction. It is not the radial
    source of an unwrapped periodic throat array or a transverse-y compact
    carrier. Individual stresses below are compactification contributions.
    """
    k, _ = _parameters(length, radius, angular_max, 1)
    if isinstance(cutoff, bool) or not np.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("cutoff must be finite and positive")
    if eta not in (0., .5):
        raise ValueError("use periodic or antiperiodic spatial spin structure")
    k = k[:, 0]
    area = 4 * pi * radius**2

    def integrand(t):
        angular_terms = k * np.exp(-t * (k / radius)**2)
        angular = angular_terms.sum()
        angular3 = np.sum(k**2 * angular_terms)
        winding, second_winding = winding_theta_sums(length**2 / (4 * t), eta)
        energy = 2 * length / pi * angular * winding / t**2
        derivative = (2 / pi * angular * winding / t**2
                      - length**2 / pi * angular * second_winding / t**3)
        null = length**2 / (4 * pi**2 * radius**2) * angular * second_winding / t**3
        radius_derivative = 4 * length / (pi * radius**3) * angular3 * winding / t
        return np.array([energy, derivative, null, radius_derivative])

    lower = cutoff**-2
    values, error = quad_vec(integrand, lower, np.inf, epsabs=2e-12, epsrel=2e-12)
    energy, derivative, null, radius_derivative = map(float, values)
    density, radial = energy / (area * length), -derivative / area
    transverse = -radius_derivative / (8 * pi * radius * length)
    cutoff_derivative = float(2 / cutoff**3 * integrand(lower)[0])
    return {"energy": energy, "length_derivative_fixed_cutoff": derivative,
            "radius_derivative_fixed_cutoff": radius_derivative,
            "relative_energy_density": density, "relative_axial_pressure": radial,
            "relative_transverse_pressure": transverse,
            "axial_null_stress_integral": null,
            "axial_null_stress_from_energy_variation": density + radial,
            "cutoff_derivative": cutoff_derivative,
            "relative_stress_trace": -density + radial + 2 * transverse,
            "cutoff_homogeneity_residual": energy + length * derivative + radius * radius_derivative
                                             - cutoff * cutoff_derivative,
            "quadrature_error_estimate_only": float(error)}


def spin_difference_heat_integral(length, radius=1., angular_max=64, winding_max=64):
    """Independent frequency-integrated winding heat representation.

    Delta E=(4L/pi) int_0^infty dt/t^2 [sum_k k exp(-t k^2/a^2)]
        [sum_{positive odd n} exp(-n^2 L^2/(4t))].
    Frequency and circle Poisson factors are both present. This quadrature
    uses Gaussian sums, not the Bessel formula. Its reported error estimates
    quadrature only; angular/winding truncations are checked separately.
    """
    k, n = _parameters(length, radius, angular_max, winding_max)
    k = k[:, 0]
    n = n[0, ::2]

    def integrand(t):
        if t <= 0:
            return 0.
        angular = np.sum(k * np.exp(-t * (k / radius)**2))
        winding = np.exp(-(n * length)**2 / (4 * t)).sum()
        return float(4 * length / pi * angular * winding / t**2)

    return quad(integrand, 0., np.inf, epsabs=2e-12, epsrel=2e-12, limit=200)


def spatial_heat_weyl_ratio(proper_time, radius=1., length=1., angular_max=1024):
    """Check full four-component rank against 4 Volume/(4 pi t)^(3/2)."""
    _parameters(length, radius, angular_max, 1)
    if not np.isfinite(proper_time) or proper_time <= 0:
        raise ValueError("proper time must be finite and positive")
    k = np.arange(1, angular_max + 1, dtype=float)
    actual = 8 * np.sum(k * np.exp(-proper_time * (k / radius)**2))
    actual *= length / np.sqrt(4 * pi * proper_time)
    volume = 4 * pi * radius**2 * length
    leading = 4 * volume / (4 * pi * proper_time)**1.5
    return float(actual / leading)
