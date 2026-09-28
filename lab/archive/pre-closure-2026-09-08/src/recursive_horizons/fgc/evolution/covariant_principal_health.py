"""Covariant all-direction principal symbol for the ACT1 classical system.

The established spherical machinery evaluates the exact unredefined ACT1/
VAR1 equations and their reference-modified-harmonic extension.  It is the
authority for the symmetry-reduced evolution equations.  HYP2 has a different
obligation: before a spherical mechanism experiment is opened, test the full
four-dimensional principal system in directions that the spherical evolution
does not contain.

This module implements that local principal symbol in a physical orthonormal
frame.  It retains all ten symmetric metric polarizations, the regulator
``phi``, and the independent canonical matter scalar ``chi``.  The physical
block is differentiated from the convention-locked VAR1 equation; the gauge
block is the Kovacs--Reall modified-harmonic projector with the same auxiliary
normal factors used by REF1.

The code is a numerical analytic adapter.  It does not, by itself, prove a
uniform run envelope or strong hyperbolicity.  Those claims require the DOM4
and HYP2 certificate layers to add quantified, all-covector bounds and to fail
closed when any bound is unavailable.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, sqrt
from numbers import Real
from typing import Sequence

import numpy as np

from ..spherical_reduction import BASE_FIELD_ORDER, SphericalState, direct_4d_curvature
from .floating_jet import scalar_primal


DIMENSION = 4
METRIC_COMPONENTS = tuple(
    (first, second)
    for first in range(DIMENSION)
    for second in range(first, DIMENSION)
)
FIELD_ORDER = tuple(
    [f"g_{first}{second}" for first, second in METRIC_COMPONENTS]
    + ["phi", "chi"]
)
FIELD_COUNT = len(FIELD_ORDER)
FIRST_ORDER_FIELD_COUNT = 2 * FIELD_COUNT
MINKOWSKI = np.diag((-1.0, 1.0, 1.0, 1.0))
FIELD_FROBENIUS_SCALES = np.asarray(
    [
        sqrt(2.0) if first != second else 1.0
        for first, second in METRIC_COMPONENTS
    ]
    + [1.0, 1.0],
    dtype=np.float64,
)


def _finite(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    return answer


def _array(name: str, value: object, shape: tuple[int, ...]) -> np.ndarray:
    try:
        answer = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite array of shape {shape}") from exc
    if answer.shape != shape or not np.all(np.isfinite(answer)):
        raise ValueError(f"{name} must be a finite array of shape {shape}")
    return answer.copy()


def _symmetric_basis(component: tuple[int, int]) -> np.ndarray:
    """Return the coordinate perturbation with one independent component.

    Off-diagonal components occur in both symmetric tensor slots.  Equation
    rows use the corresponding variational contraction, so their off-diagonal
    response is counted twice.  With that convention the ungauged principal
    matrix is the Hessian of the action and must be symmetric.
    """

    first, second = component
    answer = np.zeros((DIMENSION, DIMENSION), dtype=np.float64)
    answer[first, second] = 1.0
    answer[second, first] = 1.0
    return answer


def _variational_component(value: np.ndarray, component: tuple[int, int]) -> float:
    first, second = component
    return float(value[first, second] if first == second else 2.0 * value[first, second])


def _metric_inverse(metric: np.ndarray) -> np.ndarray:
    try:
        inverse = np.linalg.inv(metric)
    except np.linalg.LinAlgError as exc:
        raise ValueError("metric is singular") from exc
    if not np.all(np.isfinite(inverse)):
        raise ValueError("metric inverse is nonfinite")
    return inverse


def _principal_riemann(metric_perturbation: np.ndarray, covector: np.ndarray) -> np.ndarray:
    """Return the principal variation of ``R_abcd`` in the locked convention."""

    answer = np.zeros((DIMENSION,) * 4, dtype=np.float64)
    for a in range(DIMENSION):
        for b in range(DIMENSION):
            for c in range(DIMENSION):
                for d in range(DIMENSION):
                    answer[a, b, c, d] = 0.5 * (
                        covector[b] * covector[c] * metric_perturbation[a, d]
                        + covector[a] * covector[d] * metric_perturbation[b, c]
                        - covector[b] * covector[d] * metric_perturbation[a, c]
                        - covector[a] * covector[c] * metric_perturbation[b, d]
                    )
    return answer


def _ricci_and_scalar(
    riemann_lower: np.ndarray,
    inverse_metric: np.ndarray,
) -> tuple[np.ndarray, float]:
    ricci = np.einsum("ac,abcd->bd", inverse_metric, riemann_lower)
    scalar = float(np.einsum("bd,bd", inverse_metric, ricci))
    return ricci, scalar


def _double_dual_lower(
    riemann_lower: np.ndarray,
    metric: np.ndarray,
    inverse_metric: np.ndarray,
) -> np.ndarray:
    ricci, scalar = _ricci_and_scalar(riemann_lower, inverse_metric)
    answer = np.zeros((DIMENSION,) * 4, dtype=np.float64)
    for a in range(DIMENSION):
        for b in range(DIMENSION):
            for c in range(DIMENSION):
                for d in range(DIMENSION):
                    answer[a, b, c, d] = (
                        riemann_lower[a, b, c, d]
                        - metric[a, c] * ricci[d, b]
                        + metric[a, d] * ricci[c, b]
                        + metric[b, c] * ricci[d, a]
                        - metric[b, d] * ricci[c, a]
                        + 0.5
                        * scalar
                        * (
                            metric[a, c] * metric[d, b]
                            - metric[a, d] * metric[c, b]
                        )
                    )
    return answer


def _gauss_bonnet_principal_variation(
    background_riemann_lower: np.ndarray,
    principal_riemann_lower: np.ndarray,
    inverse_metric: np.ndarray,
) -> float:
    background_ricci, background_scalar = _ricci_and_scalar(
        background_riemann_lower, inverse_metric
    )
    principal_ricci, principal_scalar = _ricci_and_scalar(
        principal_riemann_lower, inverse_metric
    )
    background_ricci_up = inverse_metric @ background_ricci @ inverse_metric
    background_riemann_up = np.einsum(
        "ae,bf,cg,dh,efgh->abcd",
        inverse_metric,
        inverse_metric,
        inverse_metric,
        inverse_metric,
        background_riemann_lower,
    )
    return float(
        2.0 * background_scalar * principal_scalar
        - 8.0 * np.einsum("ab,ab", background_ricci_up, principal_ricci)
        + 2.0
        * np.einsum(
            "abcd,abcd", background_riemann_up, principal_riemann_lower
        )
    )


def _trace_reversal_projector(
    inverse_metric: np.ndarray,
    alpha: int,
    beta: int,
    mu: int,
    nu: int,
) -> float:
    return 0.5 * (
        float(alpha == mu) * inverse_metric[nu, beta]
        + float(alpha == nu) * inverse_metric[mu, beta]
        - float(alpha == beta) * inverse_metric[mu, nu]
    )


def _auxiliary_inverse(normal_factor: float) -> np.ndarray:
    if normal_factor <= 1.0:
        raise ValueError("auxiliary normal factor must exceed one")
    return np.diag((-normal_factor, 1.0, 1.0, 1.0))


@dataclass(frozen=True, slots=True)
class CovariantPrincipalBackground:
    """Local ACT1 principal data in a physical orthonormal frame."""

    effective_planck_squared: float
    effective_planck_prime: float
    gb_coupling_prime: float
    riemann_lower: np.ndarray
    hessian_gb_lower: np.ndarray

    def __post_init__(self) -> None:
        effective_planck = _finite(
            "effective_planck_squared", self.effective_planck_squared
        )
        if effective_planck <= 0.0:
            raise ValueError("effective Planck coefficient must be positive")
        object.__setattr__(self, "effective_planck_squared", effective_planck)
        object.__setattr__(
            self,
            "effective_planck_prime",
            _finite("effective_planck_prime", self.effective_planck_prime),
        )
        object.__setattr__(
            self,
            "gb_coupling_prime",
            _finite("gb_coupling_prime", self.gb_coupling_prime),
        )
        riemann = _array("riemann_lower", self.riemann_lower, (DIMENSION,) * 4)
        hessian = _array(
            "hessian_gb_lower", self.hessian_gb_lower, (DIMENSION, DIMENSION)
        )
        scale = max(1.0, float(np.max(np.abs(riemann))))
        tolerance = 2048.0 * np.finfo(np.float64).eps * scale
        defects = (
            np.max(np.abs(riemann + np.swapaxes(riemann, 0, 1))),
            np.max(np.abs(riemann + np.swapaxes(riemann, 2, 3))),
            np.max(np.abs(riemann - np.transpose(riemann, (2, 3, 0, 1)))),
            np.max(
                np.abs(
                    riemann
                    + np.transpose(riemann, (0, 2, 3, 1))
                    + np.transpose(riemann, (0, 3, 1, 2))
                )
            ),
        )
        if max(defects) > tolerance:
            raise ValueError("riemann_lower violates its algebraic curvature identities")
        if not np.allclose(hessian, hessian.T, rtol=0.0, atol=tolerance):
            raise ValueError("hessian_gb_lower must be symmetric")
        object.__setattr__(self, "riemann_lower", riemann)
        object.__setattr__(self, "hessian_gb_lower", hessian)


@dataclass(frozen=True, slots=True)
class PrincipalCoefficientTensors:
    """Quadratic covector coefficients in the Frobenius-normalized basis.

    ``P(xi)=A*xi0^2 + sum_i B_i*xi0*n_i + sum_ij C_ij*n_i*n_j``.
    The two spatial tensor indices on ``C`` are explicitly symmetric.
    """

    time_time: np.ndarray
    time_space: np.ndarray
    space_space: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "time_time",
            _array("time_time", self.time_time, (FIELD_COUNT, FIELD_COUNT)),
        )
        object.__setattr__(
            self,
            "time_space",
            _array(
                "time_space",
                self.time_space,
                (3, FIELD_COUNT, FIELD_COUNT),
            ),
        )
        space = _array(
            "space_space",
            self.space_space,
            (3, 3, FIELD_COUNT, FIELD_COUNT),
        )
        tolerance = 32768.0 * np.finfo(np.float64).eps * max(
            1.0, float(np.max(np.abs(space)))
        )
        if not np.allclose(
            space,
            np.swapaxes(space, 0, 1),
            rtol=0.0,
            atol=tolerance,
        ):
            raise ValueError("space_space covector indices must be symmetric")
        object.__setattr__(self, "space_space", space)


def esf_background(planck_mass: Real = 2.0) -> CovariantPrincipalBackground:
    """Return the exact flat Einstein--two-scalar reference background."""

    mass = _finite("planck_mass", planck_mass)
    if mass <= 0.0:
        raise ValueError("planck_mass must be positive")
    return CovariantPrincipalBackground(
        effective_planck_squared=mass * mass,
        effective_planck_prime=0.0,
        gb_coupling_prime=0.0,
        riemann_lower=np.zeros((DIMENSION,) * 4),
        hessian_gb_lower=np.zeros((DIMENSION, DIMENSION)),
    )


def _orthonormal_frame(state: SphericalState) -> np.ndarray:
    h_tt = scalar_primal(state.h_tt.value)
    h_tr = scalar_primal(state.h_tr.value)
    h_rr = scalar_primal(state.h_rr.value)
    radius = scalar_primal(state.areal_radius.value)
    if h_rr <= 0.0 or radius <= 0.0:
        raise ValueError("spherical orthonormal frame requires h_rr>0 and R>0")
    orbit = np.asarray(((h_tt, h_tr), (h_tr, h_rr)), dtype=np.float64)
    inverse = _metric_inverse(orbit)
    if inverse[0, 0] >= 0.0:
        raise ValueError("coordinate-time covector must be timelike")
    lapse = 1.0 / sqrt(-inverse[0, 0])
    normal = -lapse * inverse[:, 0]
    radial = np.asarray((0.0, 1.0 / sqrt(h_rr)))
    frame = np.zeros((DIMENSION, DIMENSION), dtype=np.float64)
    frame[:2, 0] = normal
    frame[:2, 1] = radial
    frame[2, 2] = 1.0 / radius
    frame[3, 3] = 1.0 / radius
    metric = np.zeros((DIMENSION, DIMENSION), dtype=np.float64)
    metric[:2, :2] = orbit
    metric[2, 2] = radius * radius
    metric[3, 3] = radius * radius
    pulled = frame.T @ metric @ frame
    tolerance = 4096.0 * np.finfo(np.float64).eps * max(
        1.0, float(np.max(np.abs(metric)))
    )
    if not np.allclose(pulled, MINKOWSKI, rtol=0.0, atol=tolerance):
        raise ValueError("constructed spherical frame is not orthonormal")
    return frame


def _spherical_scalar_hessian(state: SphericalState, connection: object) -> np.ndarray:
    gamma = np.asarray(connection, dtype=object)
    gradient = np.asarray(
        (
            scalar_primal(state.phi.dt),
            scalar_primal(state.phi.dr),
            0.0,
            0.0,
        ),
        dtype=np.float64,
    )
    second = np.zeros((DIMENSION, DIMENSION), dtype=np.float64)
    second[0, 0] = scalar_primal(state.phi.dtt)
    second[0, 1] = second[1, 0] = scalar_primal(state.phi.dtr)
    second[1, 1] = scalar_primal(state.phi.drr)
    answer = second.copy()
    for a in range(DIMENSION):
        for b in range(DIMENSION):
            answer[a, b] -= sum(
                scalar_primal(gamma[c, a, b]) * gradient[c]
                for c in range(DIMENSION)
            )
    return answer


def background_from_spherical_state(state: SphericalState) -> CovariantPrincipalBackground:
    """Map one complete spherical second jet into local covariant HYP2 data."""

    if not isinstance(state, SphericalState):
        raise TypeError("state must be a SphericalState")
    if state.branch != "FGC-QR":
        raise ValueError("HYP2 principal background currently targets FGC-QR")
    if tuple(BASE_FIELD_ORDER) != (
        "h_tt",
        "h_tr",
        "h_rr",
        "areal_radius",
        "phi",
        "chi",
    ):
        raise RuntimeError("spherical field order changed")
    riemann_coordinate, _inverse, connection = direct_4d_curvature(state)
    riemann_coordinate_array = np.vectorize(scalar_primal)(
        np.asarray(riemann_coordinate, dtype=object)
    ).astype(np.float64)
    hessian_coordinate = _spherical_scalar_hessian(state, connection)
    frame = _orthonormal_frame(state)
    riemann_orthonormal = np.einsum(
        "ma,nb,pc,qd,mnpq->abcd",
        frame,
        frame,
        frame,
        frame,
        riemann_coordinate_array,
    )
    hessian_phi = frame.T @ hessian_coordinate @ frame
    gradient_coordinate = np.asarray(
        (
            scalar_primal(state.phi.dt),
            scalar_primal(state.phi.dr),
            0.0,
            0.0,
        )
    )
    gradient_phi = frame.T @ gradient_coordinate
    phi = scalar_primal(state.phi.value)
    effective_planck = float(state.planck_mass) ** 2 + float(state.beta) * phi * phi
    effective_planck_prime = 2.0 * float(state.beta) * phi
    gb_prime = float(state.eta) * phi / 4.0
    gb_second = float(state.eta) / 4.0
    hessian_gb = (
        gb_prime * hessian_phi
        + gb_second * np.outer(gradient_phi, gradient_phi)
    )
    return CovariantPrincipalBackground(
        effective_planck_squared=effective_planck,
        effective_planck_prime=effective_planck_prime,
        gb_coupling_prime=gb_prime,
        riemann_lower=riemann_orthonormal,
        hessian_gb_lower=hessian_gb,
    )


def ungauged_action_principal_symbol(
    background: CovariantPrincipalBackground,
    covector: Sequence[Real],
) -> np.ndarray:
    """Return the symmetric twelve-field action-Hessian principal matrix.

    VAR1 writes the metric equation with lower indices, while the independent
    metric variables here are the covariant components ``g_ab``.  The metric
    Euler--Lagrange row conjugate to those variables is therefore
    ``-E^ab/2``.  Raising both equation indices and retaining that minus sign
    is essential: using ``E_ab/2`` is not an action Hessian and fails even on
    the flat Einstein reference.
    """

    if not isinstance(background, CovariantPrincipalBackground):
        raise TypeError("background must be CovariantPrincipalBackground")
    xi = _array("covector", covector, (DIMENSION,))
    inverse = MINKOWSKI
    metric = MINKOWSKI
    xi_up = inverse @ xi
    xi_squared = float(xi @ xi_up)
    riemann = background.riemann_lower
    ricci, scalar = _ricci_and_scalar(riemann, inverse)
    double_dual = _double_dual_lower(riemann, metric, inverse)
    hessian_gb_up = inverse @ background.hessian_gb_lower @ inverse
    answer = np.zeros((FIELD_COUNT, FIELD_COUNT), dtype=np.float64)

    metric_scalar_cross = np.zeros((DIMENSION, DIMENSION), dtype=np.float64)
    for a in range(DIMENSION):
        for b in range(DIMENSION):
            metric_scalar_cross[a, b] = 0.5 * background.effective_planck_prime * (
                metric[a, b] * xi_squared - xi[a] * xi[b]
            ) + 4.0 * background.gb_coupling_prime * sum(
                double_dual[a, c, b, d] * xi_up[c] * xi_up[d]
                for c in range(DIMENSION)
                for d in range(DIMENSION)
            )

    for column, component in enumerate(METRIC_COMPONENTS):
        perturbation = _symmetric_basis(component)
        delta_riemann = _principal_riemann(perturbation, xi)
        delta_ricci, delta_scalar = _ricci_and_scalar(delta_riemann, inverse)
        delta_double_dual = _double_dual_lower(delta_riemann, metric, inverse)
        delta_einstein = delta_ricci - 0.5 * metric * delta_scalar
        response_lower = (
            0.5 * background.effective_planck_squared * delta_einstein
            + 4.0
            * np.einsum("acbd,cd->ab", delta_double_dual, hessian_gb_up)
        )
        response = -inverse @ response_lower @ inverse
        for row, equation_component in enumerate(METRIC_COMPONENTS):
            answer[row, column] = _variational_component(
                response, equation_component
            )
        scalar_metric = (
            0.5 * background.effective_planck_prime * delta_scalar
            + background.gb_coupling_prime
            * _gauss_bonnet_principal_variation(riemann, delta_riemann, inverse)
        )
        answer[10, column] = scalar_metric

    metric_scalar_cross = -inverse @ metric_scalar_cross @ inverse
    for row, equation_component in enumerate(METRIC_COMPONENTS):
        answer[row, 10] = _variational_component(
            metric_scalar_cross, equation_component
        )
    answer[10, 10] = xi_squared
    answer[11, 11] = xi_squared

    symmetry_defect = float(np.max(np.abs(answer - answer.T)))
    tolerance = 32768.0 * np.finfo(np.float64).eps * max(
        1.0, float(np.max(np.abs(answer)))
    )
    if symmetry_defect > tolerance:
        raise ValueError(
            "ungauged ACT1 principal matrix violates action-Hessian symmetry: "
            f"defect={symmetry_defect:.17g}, tolerance={tolerance:.17g}"
        )
    return 0.5 * (answer + answer.T)


def modified_harmonic_gauge_principal_symbol(
    background: CovariantPrincipalBackground,
    covector: Sequence[Real],
    *,
    tilde_normal_factor: Real = 4.0,
    hat_normal_factor: Real = 9.0,
) -> np.ndarray:
    """Return the full ten-metric Kovacs--Reall gauge principal block.

    The block is expressed in the same ``-E^ab/2`` metric-row normalization
    as :func:`ungauged_action_principal_symbol`.  Kovacs--Reall's gauge term
    contributes ``-F * total`` to ``E^ab`` in the locked convention, hence
    ``+F * total / 2`` here.
    """

    if not isinstance(background, CovariantPrincipalBackground):
        raise TypeError("background must be CovariantPrincipalBackground")
    xi = _array("covector", covector, (DIMENSION,))
    tilde_factor = _finite("tilde_normal_factor", tilde_normal_factor)
    hat_factor = _finite("hat_normal_factor", hat_normal_factor)
    if not 1.0 < tilde_factor < hat_factor:
        raise ValueError("auxiliary normal factors must satisfy 1<tilde<hat")
    tilde = _auxiliary_inverse(tilde_factor)
    hat = _auxiliary_inverse(hat_factor)
    answer = np.zeros((FIELD_COUNT, FIELD_COUNT), dtype=np.float64)
    scale = 0.5 * background.effective_planck_squared

    for column, component in enumerate(METRIC_COMPONENTS):
        perturbation = _symmetric_basis(component)
        response_up = np.zeros((DIMENSION, DIMENSION), dtype=np.float64)
        for mu in range(DIMENSION):
            for nu in range(DIMENSION):
                total = 0.0
                for alpha in range(DIMENSION):
                    for beta in range(DIMENSION):
                        if MINKOWSKI[alpha, beta] == 0.0:
                            continue
                        hat_xi = sum(
                            _trace_reversal_projector(
                                hat, alpha, gamma, mu, nu
                            )
                            * xi[gamma]
                            for gamma in range(DIMENSION)
                        )
                        tilde_xi_h = sum(
                            _trace_reversal_projector(
                                tilde, beta, delta, rho, sigma
                            )
                            * xi[delta]
                            * perturbation[rho, sigma]
                            for delta in range(DIMENSION)
                            for rho in range(DIMENSION)
                            for sigma in range(DIMENSION)
                        )
                        total += (
                            hat_xi
                            * MINKOWSKI[alpha, beta]
                            * tilde_xi_h
                        )
                response_up[mu, nu] = scale * total
        for row, equation_component in enumerate(METRIC_COMPONENTS):
            answer[row, column] = _variational_component(
                response_up, equation_component
            )
    return answer


def complete_modified_harmonic_principal_symbol(
    background: CovariantPrincipalBackground,
    covector: Sequence[Real],
    *,
    tilde_normal_factor: Real = 4.0,
    hat_normal_factor: Real = 9.0,
) -> np.ndarray:
    """Return the complete 12-by-12 ACT1 modified-harmonic symbol."""

    return ungauged_action_principal_symbol(
        background, covector
    ) + modified_harmonic_gauge_principal_symbol(
        background,
        covector,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )


def frobenius_normalized_principal_symbol(symbol: object) -> np.ndarray:
    """Express a field symbol in an orthonormal symmetric-tensor basis.

    The raw independent coordinates store an off-diagonal metric component
    once even though it occupies two tensor slots.  Multiplication by
    ``sqrt(2)`` makes the field norm equal to the spacetime Frobenius norm.
    In that basis spatial rotations act orthogonally, which is essential for
    direction-independent matrix-norm estimates.
    """

    matrix = _array("symbol", symbol, (FIELD_COUNT, FIELD_COUNT))
    inverse_scale = np.diag(1.0 / FIELD_FROBENIUS_SCALES)
    return inverse_scale @ matrix @ inverse_scale


def frobenius_normalized_first_order_matrix(matrix: object) -> np.ndarray:
    """Apply the corresponding similarity transform to a companion matrix."""

    value = _array(
        "first_order_matrix",
        matrix,
        (FIRST_ORDER_FIELD_COUNT, FIRST_ORDER_FIELD_COUNT),
    )
    scale = np.diag(
        np.concatenate((FIELD_FROBENIUS_SCALES, FIELD_FROBENIUS_SCALES))
    )
    inverse_scale = np.diag(
        1.0 / np.concatenate(
            (FIELD_FROBENIUS_SCALES, FIELD_FROBENIUS_SCALES)
        )
    )
    return scale @ value @ inverse_scale


def principal_coefficient_tensors(
    background: CovariantPrincipalBackground,
    *,
    tilde_normal_factor: Real = 4.0,
    hat_normal_factor: Real = 9.0,
    include_gauge_extension: bool = True,
) -> PrincipalCoefficientTensors:
    """Extract every covector coefficient without directional sampling."""

    if not isinstance(include_gauge_extension, bool):
        raise TypeError("include_gauge_extension must be boolean")

    def evaluate(covector: Sequence[Real]) -> np.ndarray:
        raw = (
            complete_modified_harmonic_principal_symbol(
                background,
                covector,
                tilde_normal_factor=tilde_normal_factor,
                hat_normal_factor=hat_normal_factor,
            )
            if include_gauge_extension
            else ungauged_action_principal_symbol(background, covector)
        )
        return frobenius_normalized_principal_symbol(raw)

    origin = np.zeros(4, dtype=np.float64)
    time = origin.copy()
    time[0] = 1.0
    time_time = evaluate(time)
    time_space = np.zeros((3, FIELD_COUNT, FIELD_COUNT), dtype=np.float64)
    space_space = np.zeros(
        (3, 3, FIELD_COUNT, FIELD_COUNT), dtype=np.float64
    )
    spatial_basis = []
    for index in range(3):
        vector = origin.copy()
        vector[index + 1] = 1.0
        spatial_basis.append(vector)
        space_space[index, index] = evaluate(vector)

        plus = vector.copy()
        plus[0] = 1.0
        minus = vector.copy()
        minus[0] = -1.0
        time_space[index] = 0.5 * (evaluate(plus) - evaluate(minus))

    for first in range(3):
        for second in range(first + 1, 3):
            plus = spatial_basis[first] + spatial_basis[second]
            minus = spatial_basis[first] - spatial_basis[second]
            mixed = 0.25 * (evaluate(plus) - evaluate(minus))
            space_space[first, second] = mixed
            space_space[second, first] = mixed
    return PrincipalCoefficientTensors(
        time_time=time_time,
        time_space=time_space,
        space_space=space_space,
    )


def coefficient_tensor_symbol(
    coefficients: PrincipalCoefficientTensors,
    covector: Sequence[Real],
) -> np.ndarray:
    """Reassemble a normalized symbol from its extracted coefficients."""

    if not isinstance(coefficients, PrincipalCoefficientTensors):
        raise TypeError("coefficients must be PrincipalCoefficientTensors")
    xi = _array("covector", covector, (DIMENSION,))
    answer = coefficients.time_time * xi[0] ** 2
    for index in range(3):
        answer = answer + coefficients.time_space[index] * xi[0] * xi[index + 1]
    for first in range(3):
        for second in range(3):
            answer = (
                answer
                + coefficients.space_space[first, second]
                * xi[first + 1]
                * xi[second + 1]
            )
    return answer


def normalized_first_order_matrix_from_coefficients(
    coefficients: PrincipalCoefficientTensors,
    spatial_covector: Sequence[Real],
) -> tuple[np.ndarray, dict[str, float]]:
    """Build the normalized 24-field companion matrix for one unit direction."""

    if not isinstance(coefficients, PrincipalCoefficientTensors):
        raise TypeError("coefficients must be PrincipalCoefficientTensors")
    spatial = _array("spatial_covector", spatial_covector, (3,))
    norm = float(np.linalg.norm(spatial))
    if abs(norm - 1.0) > 4096.0 * np.finfo(np.float64).eps:
        raise ValueError("spatial_covector must have Euclidean unit norm")
    a_matrix = coefficients.time_time
    b_matrix = sum(
        (
            coefficients.time_space[index] * spatial[index]
            for index in range(3)
        ),
        np.zeros_like(a_matrix),
    )
    c_matrix = np.einsum(
        "i,ijab,j->ab",
        spatial,
        coefficients.space_space,
        spatial,
    )
    try:
        solved_c = np.linalg.solve(a_matrix, c_matrix)
        solved_b = np.linalg.solve(a_matrix, b_matrix)
        condition = float(np.linalg.cond(a_matrix, p=np.inf))
    except np.linalg.LinAlgError as exc:
        raise ValueError("all-direction kinetic block is singular") from exc
    zero = np.zeros((FIELD_COUNT, FIELD_COUNT), dtype=np.float64)
    identity = np.eye(FIELD_COUNT, dtype=np.float64)
    matrix = np.block([[zero, identity], [-solved_c, -solved_b]])
    return matrix, {
        "kinetic_condition_infinity": condition,
        "kinetic_determinant": float(np.linalg.det(a_matrix)),
        "kinetic_smallest_singular_value": float(
            np.linalg.svd(a_matrix, compute_uv=False)[-1]
        ),
    }


def second_order_coefficient_matrices(
    background: CovariantPrincipalBackground,
    spatial_covector: Sequence[Real],
    *,
    tilde_normal_factor: Real = 4.0,
    hat_normal_factor: Real = 9.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``A,B,C`` for ``P((xi0,xi_i))=A xi0^2+B xi0+C``."""

    spatial = _array("spatial_covector", spatial_covector, (3,))
    norm = float(np.linalg.norm(spatial))
    if abs(norm - 1.0) > 4096.0 * np.finfo(np.float64).eps:
        raise ValueError("spatial_covector must have Euclidean unit norm")

    def evaluate(time_component: float) -> np.ndarray:
        return complete_modified_harmonic_principal_symbol(
            background,
            np.concatenate(([time_component], spatial)),
            tilde_normal_factor=tilde_normal_factor,
            hat_normal_factor=hat_normal_factor,
        )

    zero = evaluate(0.0)
    plus = evaluate(1.0)
    minus = evaluate(-1.0)
    time_time = 0.5 * (plus + minus) - zero
    time_space = 0.5 * (plus - minus)
    return time_time, time_space, zero


def first_order_principal_matrix(
    background: CovariantPrincipalBackground,
    spatial_covector: Sequence[Real],
    *,
    tilde_normal_factor: Real = 4.0,
    hat_normal_factor: Real = 9.0,
) -> tuple[np.ndarray, dict[str, float]]:
    """Return the 24-by-24 standard companion matrix and kinetic diagnostics."""

    a_matrix, b_matrix, c_matrix = second_order_coefficient_matrices(
        background,
        spatial_covector,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    try:
        solved_c = np.linalg.solve(a_matrix, c_matrix)
        solved_b = np.linalg.solve(a_matrix, b_matrix)
        condition = float(np.linalg.cond(a_matrix, p=np.inf))
    except np.linalg.LinAlgError as exc:
        raise ValueError("all-direction kinetic block is singular") from exc
    if not isfinite(condition):
        raise ValueError("all-direction kinetic condition is nonfinite")
    zero = np.zeros((FIELD_COUNT, FIELD_COUNT), dtype=np.float64)
    identity = np.eye(FIELD_COUNT, dtype=np.float64)
    matrix = np.block([[zero, identity], [-solved_c, -solved_b]])
    return matrix, {
        "kinetic_condition_infinity": condition,
        "kinetic_determinant": float(np.linalg.det(a_matrix)),
        "kinetic_smallest_singular_value": float(
            np.linalg.svd(a_matrix, compute_uv=False)[-1]
        ),
    }


def point_spectrum_diagnostics(
    background: CovariantPrincipalBackground,
    spatial_covector: Sequence[Real],
    *,
    tilde_normal_factor: Real = 4.0,
    hat_normal_factor: Real = 9.0,
) -> dict[str, object]:
    """Return pointwise exploratory spectrum data without promoting HYP2."""

    matrix, kinetic = first_order_principal_matrix(
        background,
        spatial_covector,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    return {
        "eigenvalues": tuple(complex(value) for value in eigenvalues),
        "maximum_abs_imaginary_part": float(np.max(np.abs(eigenvalues.imag))),
        "eigenvector_condition_2": float(np.linalg.cond(eigenvectors)),
        "kinetic": kinetic,
        "classification": "pointwise_binary64_all_polarization_spectrum_not_a_uniform_HYP2_certificate",
    }
