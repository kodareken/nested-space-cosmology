"""Fast direct-array GR-0 evaluation of the frozen REF1 equations.

The general FGC tensor evaluator deliberately favors exactness and reusable
scalar algebra.  That is the analytic authority, but its object-level tensor
loops are too expensive for the thousands of stages in the outcome-blind
GR-0 calibration.  This module evaluates the *same* GR-0 specialization with
plain NumPy tensors:

* the physical equations remain the unredefined ACT1/VAR1 Einstein--two-scalar
  equations;
* the reference modified-harmonic term is reconstructed from the same
  ``tilde=4`` and ``hat=9`` auxiliary inverse metrics;
* the flat spherical reference connection and its angular derivatives are
  retained exactly at the equator; and
* all six ADM coordinate-time accelerations are obtained from the complete
  six-row residual, not from a substituted evolution system.

NUM1's generic adapter remains the source of truth.  A premise-corrected
successor protocol must compare this direct-array route against that adapter
on nontrivial jets before any dynamic GR-0 calibration result can be admitted.
This module does not inspect or evolve an FGC-QR holdout.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from numbers import Real

import numpy as np

from .health_monitor import FIELD_ORDER


FIELD_COUNT = len(FIELD_ORDER)
FIELD_INDEX = {name: index for index, name in enumerate(FIELD_ORDER)}
SPACETIME_DIMENSION = 4


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _array(name: str, value: object, *, ndim: int) -> np.ndarray:
    answer = np.asarray(value, dtype=np.float64)
    if answer.ndim != ndim or not np.all(np.isfinite(answer)):
        raise ValueError(f"{name} must be a finite rank-{ndim} binary64 array")
    return np.ascontiguousarray(answer)


@dataclass(frozen=True, slots=True)
class GR0ResidualBatch:
    """Complete residual and diagnostic tensors for one acceleration batch."""

    full_residual: np.ndarray
    unredefined_metric_residual: np.ndarray
    scalar_residual: np.ndarray
    gauge_constraint: np.ndarray
    hamiltonian_constraint: np.ndarray
    momentum_constraint: np.ndarray
    ricci_scalar: np.ndarray
    ricci_squared: np.ndarray

    def __post_init__(self) -> None:
        arrays = {
            "full_residual": (self.full_residual, 3),
            "unredefined_metric_residual": (self.unredefined_metric_residual, 4),
            "scalar_residual": (self.scalar_residual, 3),
            "gauge_constraint": (self.gauge_constraint, 3),
            "hamiltonian_constraint": (self.hamiltonian_constraint, 2),
            "momentum_constraint": (self.momentum_constraint, 2),
            "ricci_scalar": (self.ricci_scalar, 2),
            "ricci_squared": (self.ricci_squared, 2),
        }
        for name, (value, ndim) in arrays.items():
            frozen = _array(name, value, ndim=ndim).copy()
            frozen.setflags(write=False)
            object.__setattr__(self, name, frozen)


@dataclass(frozen=True, slots=True)
class GR0AccelerationResult:
    """Unique GR-0 ADM acceleration root and its numerical margins."""

    accelerations: np.ndarray
    residuals: np.ndarray
    residual_infinity: float
    kinetic_condition_infinity_maximum: float
    refinement_iterations: int
    residual_decreased_monotonically: bool
    hamiltonian_constraint: np.ndarray
    momentum_constraint: np.ndarray
    gauge_constraint: np.ndarray
    ricci_scalar: np.ndarray
    ricci_squared: np.ndarray

    def __post_init__(self) -> None:
        for name, ndim in (
            ("accelerations", 2),
            ("residuals", 2),
            ("hamiltonian_constraint", 1),
            ("momentum_constraint", 1),
            ("gauge_constraint", 2),
            ("ricci_scalar", 1),
            ("ricci_squared", 1),
        ):
            frozen = _array(name, getattr(self, name), ndim=ndim).copy()
            frozen.setflags(write=False)
            object.__setattr__(self, name, frozen)
        for name in ("residual_infinity", "kinetic_condition_infinity_maximum"):
            value = _finite(name, getattr(self, name))
            if value < 0.0:
                raise ValueError(f"{name} must be nonnegative")
            object.__setattr__(self, name, value)
        if (
            isinstance(self.refinement_iterations, bool)
            or not isinstance(self.refinement_iterations, int)
            or not 0 <= self.refinement_iterations <= 16
        ):
            raise ValueError("refinement_iterations must lie in [0,16]")
        if not isinstance(self.residual_decreased_monotonically, bool):
            raise TypeError("residual_decreased_monotonically must be bool")


def _primitive_second_derivatives(
    dtt: np.ndarray,
    p_r: np.ndarray,
    q_r: np.ndarray,
) -> np.ndarray:
    batch, points, fields = dtt.shape
    if fields != FIELD_COUNT or p_r.shape != (points, fields) or q_r.shape != (points, fields):
        raise ValueError("primitive second-derivative shapes differ")
    answer = np.zeros((batch, points, 2, 2, fields), dtype=np.float64)
    answer[:, :, 0, 0, :] = dtt
    answer[:, :, 0, 1, :] = p_r[None, :, :]
    answer[:, :, 1, 0, :] = p_r[None, :, :]
    answer[:, :, 1, 1, :] = q_r[None, :, :]
    return answer


def _metric_jet_from_adm(
    u: np.ndarray,
    p: np.ndarray,
    q: np.ndarray,
    dtt: np.ndarray,
    p_r: np.ndarray,
    q_r: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``g,dg,ddg`` at ``theta=pi/2`` for an acceleration batch."""

    batch, points, fields = dtt.shape
    if fields != FIELD_COUNT or u.shape != (points, fields) or p.shape != u.shape or q.shape != u.shape:
        raise ValueError("ADM lower-jet shapes differ")
    primitive_first = np.stack((p, q), axis=1)  # point,direction,field
    primitive_second = _primitive_second_derivatives(dtt, p_r, q_r)

    alpha = u[:, FIELD_INDEX["alpha"]]
    shift = u[:, FIELD_INDEX["shift"]]
    radial = u[:, FIELD_INDEX["lambda"]]
    areal = u[:, FIELD_INDEX["R"]]
    if np.any(alpha <= 0.0) or np.any(radial <= 0.0) or np.any(areal <= 0.0):
        raise ValueError("direct GR-0 annular source requires alpha, lambda, and R positive")

    metric = np.zeros((points, 4, 4), dtype=np.float64)
    metric[:, 0, 0] = -alpha**2 + radial**2 * shift**2
    metric[:, 0, 1] = metric[:, 1, 0] = radial**2 * shift
    metric[:, 1, 1] = radial**2
    metric[:, 2, 2] = metric[:, 3, 3] = areal**2

    derivative = np.zeros((points, 4, 4, 4), dtype=np.float64)
    second = np.zeros((batch, points, 4, 4, 4, 4), dtype=np.float64)
    ia, iv, il, iR = (FIELD_INDEX[name] for name in ("alpha", "shift", "lambda", "R"))
    for x in range(2):
        ax = primitive_first[:, x, ia]
        vx = primitive_first[:, x, iv]
        lx = primitive_first[:, x, il]
        rx = primitive_first[:, x, iR]
        derivative[:, x, 0, 0] = -2.0 * alpha * ax + 2.0 * radial * lx * shift**2 + 2.0 * radial**2 * shift * vx
        derivative[:, x, 0, 1] = derivative[:, x, 1, 0] = 2.0 * radial * lx * shift + radial**2 * vx
        derivative[:, x, 1, 1] = 2.0 * radial * lx
        derivative[:, x, 2, 2] = derivative[:, x, 3, 3] = 2.0 * areal * rx
        for y in range(2):
            ay = primitive_first[:, y, ia]
            vy = primitive_first[:, y, iv]
            ly = primitive_first[:, y, il]
            ry = primitive_first[:, y, iR]
            axy = primitive_second[:, :, x, y, ia]
            vxy = primitive_second[:, :, x, y, iv]
            lxy = primitive_second[:, :, x, y, il]
            rxy = primitive_second[:, :, x, y, iR]
            htt = (
                -2.0 * (ax[None, :] * ay[None, :] + alpha[None, :] * axy)
                + 2.0 * (lx[None, :] * ly[None, :] + radial[None, :] * lxy) * shift[None, :] ** 2
                + 4.0 * radial[None, :] * lx[None, :] * shift[None, :] * vy[None, :]
                + 4.0 * radial[None, :] * ly[None, :] * shift[None, :] * vx[None, :]
                + 2.0 * radial[None, :] ** 2 * (vx[None, :] * vy[None, :] + shift[None, :] * vxy)
            )
            htr = (
                2.0 * (lx[None, :] * ly[None, :] + radial[None, :] * lxy) * shift[None, :]
                + 2.0 * radial[None, :] * lx[None, :] * vy[None, :]
                + 2.0 * radial[None, :] * ly[None, :] * vx[None, :]
                + radial[None, :] ** 2 * vxy
            )
            hrr = 2.0 * (lx[None, :] * ly[None, :] + radial[None, :] * lxy)
            hang = 2.0 * (rx[None, :] * ry[None, :] + areal[None, :] * rxy)
            second[:, :, x, y, 0, 0] = htt
            second[:, :, x, y, 0, 1] = second[:, :, x, y, 1, 0] = htr
            second[:, :, x, y, 1, 1] = hrr
            second[:, :, x, y, 2, 2] = second[:, :, x, y, 3, 3] = hang
    # The only nonzero angular second derivative at the equator.
    second[:, :, 2, 2, 3, 3] = -2.0 * areal[None, :] ** 2
    return metric, derivative, second


def _connection(
    metric: np.ndarray,
    derivative: np.ndarray,
    second: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    inverse = np.linalg.inv(metric)
    inverse_derivative = -np.einsum(
        "nai,neij,njb->neab", inverse, derivative, inverse, optimize=True
    )
    points = metric.shape[0]
    batch = second.shape[0]
    gamma = np.zeros((points, 4, 4, 4), dtype=np.float64)
    gamma_derivative = np.zeros((batch, points, 4, 4, 4, 4), dtype=np.float64)
    for a in range(4):
        for b in range(4):
            for c in range(4):
                first_sum = np.zeros(points, dtype=np.float64)
                for d in range(4):
                    first_sum += inverse[:, a, d] * (
                        derivative[:, b, d, c]
                        + derivative[:, c, d, b]
                        - derivative[:, d, b, c]
                    ) / 2.0
                gamma[:, a, b, c] = first_sum
                for e in range(4):
                    value = np.zeros((batch, points), dtype=np.float64)
                    for d in range(4):
                        first_bracket = (
                            derivative[:, b, d, c]
                            + derivative[:, c, d, b]
                            - derivative[:, d, b, c]
                        )
                        second_bracket = (
                            second[:, :, e, b, d, c]
                            + second[:, :, e, c, d, b]
                            - second[:, :, e, d, b, c]
                        )
                        value += (
                            inverse_derivative[:, e, a, d][None, :] * first_bracket[None, :]
                            + inverse[:, a, d][None, :] * second_bracket
                        ) / 2.0
                    gamma_derivative[:, :, e, a, b, c] = value
    return inverse, inverse_derivative, gamma, gamma_derivative


def _ricci(
    gamma: np.ndarray,
    gamma_derivative: np.ndarray,
) -> np.ndarray:
    batch, points = gamma_derivative.shape[:2]
    answer = np.zeros((batch, points, 4, 4), dtype=np.float64)
    for b in range(4):
        for d in range(4):
            value = np.zeros((batch, points), dtype=np.float64)
            for a in range(4):
                value += gamma_derivative[:, :, a, a, d, b]
                value -= gamma_derivative[:, :, d, a, a, b]
                for e in range(4):
                    value += gamma[:, a, a, e][None, :] * gamma[:, e, d, b][None, :]
                    value -= gamma[:, a, d, e][None, :] * gamma[:, e, a, b][None, :]
            answer[:, :, b, d] = value
    return answer


def _flat_reference(
    radii: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    points = radii.size
    gamma = np.zeros((points, 4, 4, 4), dtype=np.float64)
    derivative = np.zeros((points, 4, 4, 4, 4), dtype=np.float64)
    r, theta, phi = 1, 2, 3
    gamma[:, r, theta, theta] = -radii
    gamma[:, r, phi, phi] = -radii
    gamma[:, theta, r, theta] = gamma[:, theta, theta, r] = 1.0 / radii
    gamma[:, phi, r, phi] = gamma[:, phi, phi, r] = 1.0 / radii
    derivative[:, r, r, theta, theta] = -1.0
    derivative[:, r, r, phi, phi] = -1.0
    derivative[:, r, theta, r, theta] = derivative[:, r, theta, theta, r] = -1.0 / radii**2
    derivative[:, r, phi, r, phi] = derivative[:, r, phi, phi, r] = -1.0 / radii**2
    derivative[:, theta, theta, phi, phi] = 1.0
    derivative[:, theta, phi, theta, phi] = -1.0
    derivative[:, theta, phi, phi, theta] = -1.0
    return gamma, derivative


def _auxiliary_inverse_and_derivative(
    inverse: np.ndarray,
    inverse_derivative: np.ndarray,
    factor: float,
) -> tuple[np.ndarray, np.ndarray]:
    denominator = inverse[:, 0, 0]
    normal = -(
        inverse[:, :, 0, None] * inverse[:, None, :, 0]
    ) / denominator[:, None, None]
    derivative = np.empty((inverse.shape[0], 4, 4, 4), dtype=np.float64)
    for e in range(4):
        denominator_e = inverse_derivative[:, e, 0, 0]
        for a in range(4):
            for b in range(4):
                numerator = inverse[:, a, 0] * inverse[:, b, 0]
                numerator_e = (
                    inverse_derivative[:, e, a, 0] * inverse[:, b, 0]
                    + inverse[:, a, 0] * inverse_derivative[:, e, b, 0]
                )
                derivative[:, e, a, b] = (
                    -numerator_e / denominator
                    + numerator * denominator_e / denominator**2
                )
    auxiliary = inverse - (factor - 1.0) * normal
    auxiliary_derivative = inverse_derivative - (factor - 1.0) * derivative
    return auxiliary, auxiliary_derivative


def gr0_ref1_residual_batch(
    u: object,
    p: object,
    q: object,
    dtt: object,
    p_r: object,
    q_r: object,
    radii: object,
    *,
    planck_mass: Real = 2.0,
    scalar_mass: Real = 3.0,
    quartic_coupling: Real = 0.5,
    tilde_normal_factor: Real = 4.0,
    hat_normal_factor: Real = 9.0,
) -> GR0ResidualBatch:
    """Evaluate the six complete GR-0 REF1 rows for all acceleration seeds."""

    base = _array("u", u, ndim=2)
    velocity = _array("p", p, ndim=2)
    gradient = _array("q", q, ndim=2)
    acceleration = _array("dtt", dtt, ndim=3)
    mixed = _array("p_r", p_r, ndim=2)
    radial_second = _array("q_r", q_r, ndim=2)
    radius = _array("radii", radii, ndim=1)
    points = radius.size
    if base.shape != (points, FIELD_COUNT) or velocity.shape != base.shape or gradient.shape != base.shape:
        raise ValueError("GR-0 state shapes differ")
    if acceleration.shape[1:] != base.shape or mixed.shape != base.shape or radial_second.shape != base.shape:
        raise ValueError("GR-0 second-jet shapes differ")
    if np.any(radius <= 0.0):
        raise ValueError("GR-0 REF1 evaluator requires positive coordinate radii")
    mplanck = _finite("planck_mass", planck_mass, positive=True)
    mass = _finite("scalar_mass", scalar_mass, positive=True)
    quartic = _finite("quartic_coupling", quartic_coupling, positive=True)
    tilde_factor = _finite("tilde_normal_factor", tilde_normal_factor, positive=True)
    hat_factor = _finite("hat_normal_factor", hat_normal_factor, positive=True)
    if tilde_factor <= 1.0 or hat_factor <= tilde_factor:
        raise ValueError("auxiliary normal factors require 1<tilde<hat")

    metric, derivative, second = _metric_jet_from_adm(
        base, velocity, gradient, acceleration, mixed, radial_second
    )
    inverse, inverse_derivative, gamma, gamma_derivative = _connection(
        metric, derivative, second
    )
    ricci = _ricci(gamma, gamma_derivative)
    scalar = np.einsum("nij,bnij->bn", inverse, ricci, optimize=True)
    einstein = ricci - metric[None, :, :, :] * scalar[:, :, None, None] / 2.0
    ricci_squared = np.einsum(
        "nac,nbd,snab,sncd->sn", inverse, inverse, ricci, ricci, optimize=True
    )

    scalar_fields = (FIELD_INDEX["phi"], FIELD_INDEX["chi"])
    hessians: list[np.ndarray] = []
    scalar_rows: list[np.ndarray] = []
    stresses: list[np.ndarray] = []
    primitive_second = _primitive_second_derivatives(acceleration, mixed, radial_second)
    for field_index in scalar_fields:
        first = np.zeros((points, 4), dtype=np.float64)
        first[:, 0] = velocity[:, field_index]
        first[:, 1] = gradient[:, field_index]
        second_field = np.zeros((acceleration.shape[0], points, 4, 4), dtype=np.float64)
        second_field[:, :, :2, :2] = primitive_second[:, :, :, :, field_index]
        hessian = second_field.copy()
        for a in range(4):
            for b in range(4):
                for c in range(4):
                    hessian[:, :, a, b] -= gamma[:, c, a, b][None, :] * first[:, c][None, :]
        box = np.einsum("nij,bnij->bn", inverse, hessian, optimize=True)
        value = base[:, field_index]
        potential = (
            0.5 * mass**2 * value**2 + 0.25 * quartic * value**4
            if field_index == FIELD_INDEX["phi"]
            else np.zeros(points, dtype=np.float64)
        )
        gradient_square = np.einsum("nij,ni,nj->n", inverse, first, first, optimize=True)
        stress = (
            first[None, :, :, None] * first[None, :, None, :]
            - metric[None, :, :, :] * (gradient_square[None, :, None, None] / 2.0 + potential[None, :, None, None])
        )
        row = box - (
            mass**2 * value[None, :] + quartic * value[None, :] ** 3
            if field_index == FIELD_INDEX["phi"]
            else 0.0
        )
        hessians.append(hessian)
        scalar_rows.append(row)
        stresses.append(stress)

    original_metric = mplanck**2 * einstein - stresses[0] - stresses[1]

    reference_gamma, reference_derivative = _flat_reference(radius)
    difference = gamma - reference_gamma
    difference_derivative = gamma_derivative - reference_derivative[None, :, :, :, :, :]
    tilde, tilde_derivative = _auxiliary_inverse_and_derivative(
        inverse, inverse_derivative, tilde_factor
    )
    constraint = -np.einsum("nrs,nars->na", tilde, difference, optimize=True)
    partial_constraint = np.empty((acceleration.shape[0], points, 4, 4), dtype=np.float64)
    for direction in range(4):
        for a in range(4):
            value = np.zeros((acceleration.shape[0], points), dtype=np.float64)
            for rho in range(4):
                for sigma in range(4):
                    value -= (
                        tilde_derivative[:, direction, rho, sigma][None, :]
                        * difference[:, a, rho, sigma][None, :]
                        + tilde[:, rho, sigma][None, :]
                        * difference_derivative[:, :, direction, a, rho, sigma]
                    )
            partial_constraint[:, :, direction, a] = value
    covariant_constraint_derivative = partial_constraint.copy()
    for direction in range(4):
        for a in range(4):
            for c in range(4):
                covariant_constraint_derivative[:, :, direction, a] += (
                    gamma[:, a, direction, c][None, :] * constraint[:, c][None, :]
                )

    hat, _hat_derivative = _auxiliary_inverse_and_derivative(
        inverse, inverse_derivative, hat_factor
    )
    extension_up = np.zeros((acceleration.shape[0], points, 4, 4), dtype=np.float64)
    for mu in range(4):
        for nu in range(4):
            value = np.zeros((acceleration.shape[0], points), dtype=np.float64)
            for alpha in range(4):
                for beta in range(4):
                    projector = (
                        (1.0 if alpha == mu else 0.0) * hat[:, nu, beta]
                        + (1.0 if alpha == nu else 0.0) * hat[:, mu, beta]
                        - (1.0 if alpha == beta else 0.0) * hat[:, mu, nu]
                    ) / 2.0
                    value += projector[None, :] * covariant_constraint_derivative[:, :, beta, alpha]
            extension_up[:, :, mu, nu] = mplanck**2 * value
    extension_down = np.einsum(
        "nam,snmp,npb->snab", metric, extension_up, metric, optimize=True
    )
    full_metric = original_metric + extension_down
    full = np.stack(
        (
            full_metric[:, :, 0, 0],
            full_metric[:, :, 0, 1],
            full_metric[:, :, 1, 1],
            full_metric[:, :, 2, 2],
            scalar_rows[0],
            scalar_rows[1],
        ),
        axis=2,
    )
    shift = metric[:, 0, 1] / metric[:, 1, 1]
    hamiltonian = (
        original_metric[:, :, 0, 0]
        - 2.0 * shift[None, :] * original_metric[:, :, 0, 1]
        + shift[None, :] ** 2 * original_metric[:, :, 1, 1]
    )
    momentum = original_metric[:, :, 0, 1] - shift[None, :] * original_metric[:, :, 1, 1]
    return GR0ResidualBatch(
        full_residual=full,
        unredefined_metric_residual=original_metric,
        scalar_residual=np.stack(scalar_rows, axis=2),
        gauge_constraint=np.broadcast_to(constraint[None, :, :], (acceleration.shape[0], points, 4)),
        hamiltonian_constraint=hamiltonian,
        momentum_constraint=momentum,
        ricci_scalar=scalar,
        ricci_squared=ricci_squared,
    )


def solve_gr0_grid_accelerations(
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radii: object,
    *,
    residual_tolerance: Real = 1.0e-12,
    condition_number_maximum: Real = 1.0e10,
    maximum_refinement_iterations: int = 16,
) -> GR0AccelerationResult:
    """Solve the unique affine GR-0 REF1 acceleration root on an annulus."""

    base = _array("u", u, ndim=2)
    points, fields = base.shape
    if fields != FIELD_COUNT:
        raise ValueError("u must contain the six canonical ADM fields")
    seeds = np.zeros((FIELD_COUNT + 1, points, FIELD_COUNT), dtype=np.float64)
    for field in range(FIELD_COUNT):
        seeds[field + 1, :, field] = 1.0
    evaluated = gr0_ref1_residual_batch(base, p, q, seeds, p_r, q_r, radii)
    constant = evaluated.full_residual[0]
    jacobian = np.empty((points, FIELD_COUNT, FIELD_COUNT), dtype=np.float64)
    for field in range(FIELD_COUNT):
        jacobian[:, :, field] = evaluated.full_residual[field + 1] - constant
    try:
        inverse = np.linalg.inv(jacobian)
        acceleration = np.linalg.solve(jacobian, -constant[..., None])[..., 0]
    except np.linalg.LinAlgError as exc:
        raise ValueError("direct GR-0 kinetic block is singular") from exc
    condition = np.max(np.sum(np.abs(jacobian), axis=2), axis=1) * np.max(
        np.sum(np.abs(inverse), axis=2), axis=1
    )
    condition_limit = _finite(
        "condition_number_maximum", condition_number_maximum, positive=True
    )
    if not np.all(np.isfinite(condition)) or np.any(condition >= condition_limit):
        raise ValueError("direct GR-0 kinetic condition limit reached")
    tolerance = _finite("residual_tolerance", residual_tolerance, positive=True)
    if (
        isinstance(maximum_refinement_iterations, bool)
        or not isinstance(maximum_refinement_iterations, int)
        or not 0 <= maximum_refinement_iterations <= 16
    ):
        raise ValueError("maximum_refinement_iterations must lie in [0,16]")
    if not np.all(np.isfinite(acceleration)):
        raise ValueError("direct GR-0 acceleration root became nonfinite")
    verified = gr0_ref1_residual_batch(
        base,
        p,
        q,
        acceleration[None, :, :],
        p_r,
        q_r,
        radii,
    )
    residual = verified.full_residual[0]
    residual_norm = float(np.max(np.abs(residual), initial=0.0))
    history = [float(np.max(np.abs(constant), initial=0.0)), residual_norm]
    refinements = 0
    while residual_norm > tolerance and refinements < maximum_refinement_iterations:
        correction = np.linalg.solve(jacobian, -residual[..., None])[..., 0]
        acceleration = acceleration + correction
        if not np.all(np.isfinite(acceleration)):
            raise ValueError("direct GR-0 refinement root became nonfinite")
        candidate = gr0_ref1_residual_batch(
            base,
            p,
            q,
            acceleration[None, :, :],
            p_r,
            q_r,
            radii,
        )
        candidate_residual = candidate.full_residual[0]
        candidate_norm = float(
            np.max(np.abs(candidate_residual), initial=0.0)
        )
        refinements += 1
        history.append(candidate_norm)
        # The same affine Jacobian is already exact to binary64 precision.
        # Continuing after a nondecrease would turn roundoff into a fake
        # Newton iteration and can oscillate rather than improve the root.
        if candidate_norm >= residual_norm:
            break
        verified = candidate
        residual = candidate_residual
        residual_norm = candidate_norm
    monotonic = all(
        right < left
        for left, right in zip(history, history[1:], strict=False)
        if left > tolerance
    )
    if residual_norm > tolerance or not monotonic:
        raise ValueError(
            "direct GR-0 acceleration root missed its residual tolerance "
            f"(observed={residual_norm:.17g}, refinements={refinements})"
        )
    return GR0AccelerationResult(
        accelerations=acceleration,
        residuals=residual,
        residual_infinity=residual_norm,
        kinetic_condition_infinity_maximum=float(np.max(condition, initial=0.0)),
        refinement_iterations=refinements,
        residual_decreased_monotonically=monotonic,
        hamiltonian_constraint=verified.hamiltonian_constraint[0],
        momentum_constraint=verified.momentum_constraint[0],
        gauge_constraint=verified.gauge_constraint[0],
        ricci_scalar=verified.ricci_scalar[0],
        ricci_squared=verified.ricci_squared[0],
    )


__all__ = [
    "GR0AccelerationResult",
    "GR0ResidualBatch",
    "gr0_ref1_residual_batch",
    "solve_gr0_grid_accelerations",
]
