"""Exact static four-dimensional Dirac identities for the cutoff modulus.

The metric is ds_E^2=N(x)^2 d_tau^2+q(x)^2 dx^2+r(x)^2 d_Omega^2,
with smooth positive N, q, r, torsionless connection and one complex Dirac
field. Angular reduction retains both signs of kappa. These are continuum
differential identities, not proofs that a finite derivative matrix obeys
Leibniz, or that a finite regulator fixes a renormalized physical action.
The fixed-cutoff modulus does not determine its physical finite coefficients,
normalization, determinant phase, compensator or recursive completion. In
particular the old fixed-matrix-rank normalization term cannot be extended
to unrestricted frequency and angular traces as an unregulated Tr(1).

Conventions: Hermitian gamma matrices with {gamma_a,gamma_b}=2 delta_ab;
D=-i gamma^mu nabla_mu; positive scalar curvature on the round sphere.
Lichnerowicz: D^2=nabla* nabla+R/4. Heat coefficients exclude (4*pi*t)^-2.
References: https://arxiv.org/html/1106.3263v1, equations 6--8, and
https://doi.org/10.1140/epjc/s10052-020-7805-1, Appendix B and section 2.
"""
from __future__ import annotations

from functools import lru_cache

import sympy as sp


def _matrices():
    rho1 = sp.Matrix([[0, 1], [1, 0]])
    rho2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    rho3 = sp.diag(1, -1)
    identity = sp.eye(2)
    gamma = (sp.kronecker_product(identity, rho3),
             sp.kronecker_product(identity, rho1),
             sp.kronecker_product(rho2, rho2),
             sp.kronecker_product(rho3, rho2))
    alpha = sp.kronecker_product(identity, rho2)
    angular = sp.kronecker_product(rho3, rho1)
    chirality = -sp.kronecker_product(rho1, rho2)
    return gamma, alpha, angular, chirality


def _strings(matrix):
    return [[str(sp.simplify(entry)) for entry in row]
            for row in matrix.tolist()]


def _zero(matrix):
    """Evaluate every residual entry and fail if an identity does not close."""
    reduced = matrix.applyfunc(sp.simplify)
    if reduced != sp.zeros(*reduced.shape):
        raise ArithmeticError(f"Nonzero exact Dirac residual: {reduced}")
    return _strings(reduced)


def _scalar_zero(expression):
    reduced = sp.simplify(expression)
    if reduced != 0:
        raise ArithmeticError(f"Nonzero exact scalar residual: {reduced}")
    return str(reduced)


@lru_cache(maxsize=1)
def flat_measure_identities():
    """Apply the unrescaled radial connection to a general scalar amplitude."""
    x = sp.symbols("x", real=True)
    lapse = sp.Function("N", positive=True)(x)
    radial = sp.Function("q", positive=True)(x)
    radius = sp.Function("r", positive=True)(x)
    amplitude = sp.Function("f")(x)
    rescaling = radius * sp.sqrt(lapse * radial)
    connection = sp.diff(lapse, x) / (2 * lapse) + sp.diff(radius, x) / radius

    def radial_action(connection_term):
        return -sp.I * rescaling / radial * (
            sp.diff(amplitude / rescaling, x)
            + connection_term * amplitude / rescaling)

    momentum = -sp.I * (sp.diff(amplitude, x) / radial
                       + sp.diff(1 / radial, x) * amplitude / 2)
    omitted_lapse = connection - sp.diff(lapse, x) / (2 * lapse)
    return {
        "hilbert_weight_residual": _scalar_zero(
            rescaling**2 - lapse * radial * radius**2),
        "flat_radial_operator_residual": _scalar_zero(
            radial_action(connection) - momentum),
        "omitted_lapse_connection_defect": str(sp.simplify(
            radial_action(omitted_lapse) - momentum)),
        "rescaling": str(rescaling),
        "unrescaled_radial_connection": str(connection),
    }


@lru_cache(maxsize=1)
def paired_operator_identities():
    """Differentiate arbitrary spinor fields; no finite-grid product rule is used."""
    x = sp.symbols("x", real=True)
    lapse = sp.Function("N", positive=True)(x)
    radial = sp.Function("q", positive=True)(x)
    radius = sp.Function("r", positive=True)(x)
    omega, kappa = sp.symbols("omega kappa", real=True)
    sigma = sp.Function("sigma", real=True)(x)
    gamma, alpha, angular, chirality = _matrices()
    beta = gamma[0]
    eta3_rho3 = sp.I * alpha * angular
    amplitude = sp.Matrix([sp.Function(f"u{j}")(x) for j in range(4)])

    def momentum(vector, metric_q=radial):
        return -sp.I * (sp.diff(vector, x) / metric_q
                       + sp.diff(1 / metric_q, x) * vector / 2)

    def hamiltonian(vector):
        return alpha * momentum(vector) + angular * kappa / radius * vector

    def dirac(vector, metric_n=lapse, metric_q=radial, metric_r=radius):
        return beta * omega / metric_n * vector + sp.I * beta * (
            alpha * momentum(vector, metric_q) + angular * kappa / metric_r * vector)

    frequency, sphere = omega / lapse, kappa / radius
    square = (momentum(momentum(amplitude))
              + (frequency**2 + sphere**2) * amplitude
              - eta3_rho3 * sp.diff(sphere, x) / radial * amplitude
              - alpha * sp.diff(frequency, x) / radial * amplitude)
    # eta3_rho3=i alpha A, so its sign is independently fixed by Pauli products.
    square_residual = _zero(dirac(dirac(amplitude)) - square)
    scale = sp.exp(sigma)
    changed = dirac(amplitude, scale * lapse, scale * radial, scale * radius)
    sandwich = sp.exp(-sigma / 2) * dirac(sp.exp(-sigma / 2) * amplitude)
    ultrastatic = dirac(dirac(amplitude)).subs(lapse, sp.Integer(1)).doit()
    old_square = hamiltonian(hamiltonian(amplitude)) + frequency**2 * amplitude
    missing_lapse = sp.simplify(dirac(dirac(amplitude)) - old_square)
    eps = sp.symbols("epsilon", real=True)
    nq, qq, rq = (sp.Function(name, real=True)(x) for name in ("n", "a", "b"))
    direct = dirac(amplitude, sp.exp(eps * nq) * lapse,
                   sp.exp(eps * qq) * radial, sp.exp(eps * rq) * radius)
    delta_momentum = sp.I * (qq / radial * sp.diff(amplitude, x)
                            + sp.diff(qq / radial, x) * amplitude / 2)
    tangent = (-beta * frequency * nq * amplitude
               + sp.I * beta * (alpha * delta_momentum
                                 - angular * sphere * rq * amplitude))
    constant = sp.symbols("c", positive=True)
    return {
        "clifford_residuals": [_zero(gamma[a] * gamma[b] + gamma[b] * gamma[a]
                                     - 2 * int(a == b) * sp.eye(4))
                                for a in range(4) for b in range(4)],
        "gamma_product_chirality_residual": _zero(
            gamma[0] * gamma[1] * gamma[2] * gamma[3] - chirality),
        "physical_chirality_residual": _zero(
            dirac(chirality * amplitude) + chirality * dirac(amplitude)),
        "beta_hamiltonian_residual": _zero(
            hamiltonian(beta * amplitude) + beta * hamiltonian(amplitude)),
        "frequency_full_square_residual": square_residual,
        "ultrastatic_square_residual": _zero(
            ultrastatic - omega**2 * amplitude - hamiltonian(hamiltonian(amplitude))),
        "local_weyl_sandwich_residual": _zero(changed - sandwich),
        "metric_tangent_residual": _zero(sp.diff(direct, eps).subs(eps, 0) - tangent),
        "uniform_lapse_frequency_residual": _zero(
            dirac(amplitude, constant * lapse) - dirac(amplitude).subs(omega, omega / constant)),
        "omitted_lapse_square_defect": _strings(missing_lapse),
    }


@lru_cache(maxsize=1)
def lichnerowicz_identities():
    """Build the coordinate connection and contract D^2 on all spinor columns."""
    tau, x, theta, phi = sp.symbols("tau x theta phi", real=True)
    coordinates = (tau, x, theta, phi)
    lapse = sp.Function("N", positive=True)(x)
    radial = sp.Function("q", positive=True)(x)
    radius = sp.Function("r", positive=True)(x)
    metric = sp.diag(lapse**2, radial**2, radius**2, radius**2 * sp.sin(theta)**2)
    inverse = metric.inv()
    christoffel = [[[sp.simplify(sum(
        inverse[k, m] * (sp.diff(metric[m, j], coordinates[i])
                         + sp.diff(metric[m, i], coordinates[j])
                         - sp.diff(metric[i, j], coordinates[m]))
        for m in range(4)) / 2) for j in range(4)] for i in range(4)] for k in range(4)]
    ricci = sp.Matrix(4, 4, lambda i, j: sp.simplify(sum(
        sp.diff(christoffel[k][i][j], coordinates[k])
        - sp.diff(christoffel[k][i][k], coordinates[j])
        + sum(christoffel[k][k][m] * christoffel[m][i][j]
              - christoffel[k][j][m] * christoffel[m][i][k] for m in range(4))
        for k in range(4))))
    scalar = sp.simplify(sum(inverse[i, j] * ricci[i, j]
                             for i in range(4) for j in range(4)))
    gamma, _, _, _ = _matrices()
    spin = (
        sp.diff(lapse, x) / (2 * radial) * gamma[0] * gamma[1],
        sp.zeros(4),
        -sp.diff(radius, x) / (2 * radial) * gamma[1] * gamma[2],
        -sp.diff(radius, x) * sp.sin(theta) / (2 * radial) * gamma[1] * gamma[3]
        - sp.cos(theta) / 2 * gamma[2] * gamma[3],
    )
    tetrad_inverse = (1 / lapse, 1 / radial, 1 / radius, 1 / (radius * sp.sin(theta)))
    amplitude = sp.Function("f")(*coordinates) * sp.eye(4)

    def covariant(vector, i):
        return sp.diff(vector, coordinates[i]) + spin[i] * vector

    def dirac(vector):
        return sum((-sp.I * gamma[i] * tetrad_inverse[i] * covariant(vector, i)
                    for i in range(4)), sp.zeros(4))

    bochner = sp.zeros(4)
    for i in range(4):
        second = covariant(covariant(amplitude, i), i) - sum(
            (christoffel[k][i][i] * covariant(amplitude, k) for k in range(4)), sp.zeros(4))
        bochner -= inverse[i, i] * second

    def dot(value):
        return sp.diff(value, x) / radial

    scalar_expected = (2 * (1 - dot(radius)**2) / radius**2 - 4 * dot(dot(radius)) / radius
                       - 2 * dot(dot(lapse)) / lapse - 4 * dot(lapse) * dot(radius) / (lapse * radius))
    contracted = sum((gamma[i] * tetrad_inverse[i] * spin[i] for i in range(4)), sp.zeros(4))
    expected_connection = (gamma[1] / radial * (sp.diff(lapse, x) / (2 * lapse)
                                                + sp.diff(radius, x) / radius)
                           + gamma[2] * sp.cot(theta) / (2 * radius))
    return {
        "scalar_curvature": str(scalar),
        "scalar_curvature_warped_product_residual": _scalar_zero(scalar - scalar_expected),
        "contracted_spin_connection_residual": _zero(contracted - expected_connection),
        "full_unseparated_lichnerowicz_residual": _zero(dirac(dirac(amplitude))
                                                       - bochner - scalar * amplitude / 4),
        "spin_connection_coordinate_matrices": [_strings(value) for value in spin],
    }


@lru_cache(maxsize=1)
def dirac_heat_coefficients():
    """Evaluate Clifford curvature traces and the universal Laplace heat formula.

Sectional curvatures A,B,C,D parameterize the static spherical family. The
universal heat formula is imported, while its rank, endomorphism traces,
spin-curvature contraction and reduction to R,C^2,E4 are evaluated here.
"""
    gamma, _, _, _ = _matrices()
    aa, bb, cc, dd = sp.symbols("A B C D", real=True)
    sections = {(0, 1): aa, (0, 2): bb, (0, 3): bb,
                (1, 2): cc, (1, 3): cc, (2, 3): dd}
    curvature_trace = sp.Integer(0)
    for (a, b), sectional in sections.items():
        spin_curvature = sectional * gamma[a] * gamma[b] / 2
        # Both ordered pairs (a,b),(b,a) occur in Omega_ab Omega^ab.
        curvature_trace += 2 * sp.trace(spin_curvature * spin_curvature)
    riemann_squared_sections = 4 * sum(value**2 for value in sections.values())
    scalar, ricci2, riemann2, box_scalar = sp.symbols("R Ricci2 Riemann2 BoxR", real=True)
    weyl2, euler = sp.symbols("C2 E4", real=True)
    rank = sp.trace(sp.eye(4))
    endomorphism = -scalar * sp.eye(4) / 4
    curvature_coefficient = sp.cancel(curvature_trace / riemann_squared_sections)
    a0 = rank
    a2 = sp.trace(endomorphism + scalar * sp.eye(4) / 6)
    a4 = sp.expand((60 * scalar * sp.trace(endomorphism)
                    + 180 * sp.trace(endomorphism**2)
                    + 30 * curvature_coefficient * riemann2
                    + rank * (5 * scalar**2 - 2 * ricci2 + 2 * riemann2)
                    + (60 * (-rank / 4) + 12 * rank) * box_scalar) / 360)
    basis = sp.solve((sp.Eq(weyl2, riemann2 - 2 * ricci2 + scalar**2 / 3),
                      sp.Eq(euler, riemann2 - 4 * ricci2 + scalar**2)), (ricci2, riemann2))
    a4_weyl = sp.simplify(a4.subs(basis, simultaneous=True))
    radius = sp.symbols("r", positive=True)
    cylinder = {scalar: 2 / radius**2, ricci2: 2 / radius**4,
                riemann2: 4 / radius**4, box_scalar: 0}
    return {
        "spin_curvature_squared_trace": str(sp.expand(curvature_trace)),
        "spin_curvature_riemann_ratio": str(curvature_coefficient),
        "spin_curvature_contraction_residual": _scalar_zero(
            curvature_trace + riemann_squared_sections / 2),
        "a0": str(a0), "a2": str(a2), "a4": str(a4),
        "a4_weyl_euler_basis": str(a4_weyl),
        "cylinder_a0": str(a0),
        "cylinder_a2": str(sp.simplify(a2.subs(cylinder))),
        "cylinder_a4": str(sp.simplify(a4.subs(cylinder))),
    }


@lru_cache(maxsize=1)
def exact_checks():
    """Return deterministic computed expressions and zero residual strings."""
    return {"flat_measure": flat_measure_identities(),
            "paired_operator": paired_operator_identities(),
            "geometry": lichnerowicz_identities(),
            "dirac_heat": dirac_heat_coefficients()}
