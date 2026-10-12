"""Reference-covariant binary64 evaluation of the GR-0 REF1 equations.

The legacy direct evaluator first constructs the physical spherical
Christoffel symbols and their coordinate derivatives and subsequently
subtracts the flat spherical reference.  Near the reference solution, both
terms contain ``1/r`` and ``1/r**2`` contributions even though their tensorial
difference is small.  SRC2/PREF11 proved that this evaluation order loses
enough binary64 information to reject an exact affine root at a captured
GR-0 point.

This module evaluates the same equations in a reference-covariant order.  If
``h = g - g_bar`` and ``bar_nabla`` is the flat spherical reference
connection, the connection difference is formed directly as

``C^a_bc = 1/2 g^ad (bar_nabla_b h_dc + bar_nabla_c h_db - bar_nabla_d h_bc)``.

Its coordinate derivative is obtained by differentiating that identity.
Ricci and the modified-harmonic gauge terms are then assembled from ``C`` and
``partial C`` without ever subtracting two independently rounded
``1/r**2`` connection derivatives.  No equation, reference geometry, source
threshold, field, coupling, or acceleration branch is changed.

This is a numerical instrument, not an evolution authorization.  A separate
prospective SRC3 certificate must bound its domain and controls before any
fresh GR-0 campaign can consume it.
"""

from __future__ import annotations

from numbers import Real

import numpy as np

from .gr0_direct_source import (
    FIELD_COUNT,
    FIELD_INDEX,
    GR0ResidualBatch,
    _array,
    _auxiliary_inverse_and_derivative,
    _finite,
    _flat_reference,
    _metric_jet_from_adm,
    _primitive_second_derivatives,
)
from .cal2_source_diagnosis import AffineRootDiagnosis, ResidualIteration


SPACETIME_DIMENSION = 4


def _flat_reference_metric_jet(
    radii: np.ndarray,
    *,
    acceleration_batch_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the exact spherical Minkowski metric two-jet in binary64."""

    points = radii.size
    reference_u = np.zeros((points, FIELD_COUNT), dtype=np.float64)
    reference_p = np.zeros_like(reference_u)
    reference_q = np.zeros_like(reference_u)
    reference_u[:, FIELD_INDEX["alpha"]] = 1.0
    reference_u[:, FIELD_INDEX["lambda"]] = 1.0
    reference_u[:, FIELD_INDEX["R"]] = radii
    reference_q[:, FIELD_INDEX["R"]] = 1.0
    reference_dtt = np.zeros(
        (acceleration_batch_size, points, FIELD_COUNT), dtype=np.float64
    )
    return _metric_jet_from_adm(
        reference_u,
        reference_p,
        reference_q,
        reference_dtt,
        reference_p,
        reference_p,
    )


def _inverse_and_derivative(
    metric: np.ndarray,
    metric_derivative: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``g^-1`` and its first coordinate derivative."""

    inverse = np.linalg.inv(metric)
    inverse_derivative = -np.einsum(
        "nai,neij,njb->neab",
        inverse,
        metric_derivative,
        inverse,
        optimize=True,
    )
    return inverse, inverse_derivative


def _reference_covariant_connection_difference(
    metric: np.ndarray,
    metric_derivative: np.ndarray,
    metric_second_derivative: np.ndarray,
    inverse: np.ndarray,
    inverse_derivative: np.ndarray,
    radii: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return ``Gamma_bar``, ``d Gamma_bar``, ``C``, and ``partial C``.

    ``C = Gamma(g) - Gamma_bar`` is computed from ``bar_nabla(g-g_bar)``.
    The returned derivative is the coordinate derivative of that tensorial
    connection difference, not a subtraction of two connection derivatives.
    """

    batch, points = metric_second_derivative.shape[:2]
    reference_metric, reference_first, reference_second = (
        _flat_reference_metric_jet(
            radii,
            acceleration_batch_size=batch,
        )
    )
    reference_gamma, reference_gamma_derivative = _flat_reference(radii)
    perturbation = metric - reference_metric
    perturbation_first = metric_derivative - reference_first
    perturbation_second = metric_second_derivative - reference_second

    # bar_nabla_b h_dc and its ordinary coordinate derivative.
    covariant_first = np.zeros(
        (points, SPACETIME_DIMENSION, SPACETIME_DIMENSION, SPACETIME_DIMENSION),
        dtype=np.float64,
    )
    covariant_first_derivative = np.zeros(
        (
            batch,
            points,
            SPACETIME_DIMENSION,
            SPACETIME_DIMENSION,
            SPACETIME_DIMENSION,
            SPACETIME_DIMENSION,
        ),
        dtype=np.float64,
    )
    for b in range(SPACETIME_DIMENSION):
        for d in range(SPACETIME_DIMENSION):
            for c in range(SPACETIME_DIMENSION):
                value = perturbation_first[:, b, d, c].copy()
                for f in range(SPACETIME_DIMENSION):
                    value -= (
                        reference_gamma[:, f, b, d] * perturbation[:, f, c]
                        + reference_gamma[:, f, b, c] * perturbation[:, d, f]
                    )
                covariant_first[:, b, d, c] = value
                for e in range(SPACETIME_DIMENSION):
                    derivative = perturbation_second[:, :, e, b, d, c].copy()
                    for f in range(SPACETIME_DIMENSION):
                        derivative -= (
                            reference_gamma_derivative[:, e, f, b, d][None, :]
                            * perturbation[:, f, c][None, :]
                            + reference_gamma[:, f, b, d][None, :]
                            * perturbation_first[:, e, f, c][None, :]
                            + reference_gamma_derivative[:, e, f, b, c][None, :]
                            * perturbation[:, d, f][None, :]
                            + reference_gamma[:, f, b, c][None, :]
                            * perturbation_first[:, e, d, f][None, :]
                        )
                    covariant_first_derivative[:, :, e, b, d, c] = derivative

    christoffel_bracket = np.zeros(
        (points, SPACETIME_DIMENSION, SPACETIME_DIMENSION, SPACETIME_DIMENSION),
        dtype=np.float64,
    )
    christoffel_bracket_derivative = np.zeros(
        (
            batch,
            points,
            SPACETIME_DIMENSION,
            SPACETIME_DIMENSION,
            SPACETIME_DIMENSION,
            SPACETIME_DIMENSION,
        ),
        dtype=np.float64,
    )
    for d in range(SPACETIME_DIMENSION):
        for b in range(SPACETIME_DIMENSION):
            for c in range(SPACETIME_DIMENSION):
                christoffel_bracket[:, d, b, c] = (
                    covariant_first[:, b, d, c]
                    + covariant_first[:, c, d, b]
                    - covariant_first[:, d, b, c]
                )
                for e in range(SPACETIME_DIMENSION):
                    christoffel_bracket_derivative[:, :, e, d, b, c] = (
                        covariant_first_derivative[:, :, e, b, d, c]
                        + covariant_first_derivative[:, :, e, c, d, b]
                        - covariant_first_derivative[:, :, e, d, b, c]
                    )

    difference = np.zeros(
        (points, SPACETIME_DIMENSION, SPACETIME_DIMENSION, SPACETIME_DIMENSION),
        dtype=np.float64,
    )
    difference_derivative = np.zeros(
        (
            batch,
            points,
            SPACETIME_DIMENSION,
            SPACETIME_DIMENSION,
            SPACETIME_DIMENSION,
            SPACETIME_DIMENSION,
        ),
        dtype=np.float64,
    )
    for a in range(SPACETIME_DIMENSION):
        for b in range(SPACETIME_DIMENSION):
            for c in range(SPACETIME_DIMENSION):
                for d in range(SPACETIME_DIMENSION):
                    difference[:, a, b, c] += (
                        inverse[:, a, d]
                        * christoffel_bracket[:, d, b, c]
                        / 2.0
                    )
                    for e in range(SPACETIME_DIMENSION):
                        difference_derivative[:, :, e, a, b, c] += (
                            inverse_derivative[:, e, a, d][None, :]
                            * christoffel_bracket[:, d, b, c][None, :]
                            + inverse[:, a, d][None, :]
                            * christoffel_bracket_derivative[:, :, e, d, b, c]
                        ) / 2.0
    return (
        reference_gamma,
        reference_gamma_derivative,
        difference,
        difference_derivative,
    )


def _reference_covariant_ricci(
    reference_gamma: np.ndarray,
    difference: np.ndarray,
    difference_derivative: np.ndarray,
) -> np.ndarray:
    """Assemble Ricci from a flat reference and its connection difference."""

    batch, points = difference_derivative.shape[:2]
    answer = np.zeros(
        (batch, points, SPACETIME_DIMENSION, SPACETIME_DIMENSION),
        dtype=np.float64,
    )
    for b in range(SPACETIME_DIMENSION):
        for d in range(SPACETIME_DIMENSION):
            value = np.zeros((batch, points), dtype=np.float64)
            for a in range(SPACETIME_DIMENSION):
                value += difference_derivative[:, :, a, a, d, b]
                value -= difference_derivative[:, :, d, a, a, b]
                for e in range(SPACETIME_DIMENSION):
                    value += (
                        reference_gamma[:, a, a, e][None, :]
                        * difference[:, e, d, b][None, :]
                        + difference[:, a, a, e][None, :]
                        * reference_gamma[:, e, d, b][None, :]
                        + difference[:, a, a, e][None, :]
                        * difference[:, e, d, b][None, :]
                    )
                    value -= (
                        reference_gamma[:, a, d, e][None, :]
                        * difference[:, e, a, b][None, :]
                        + difference[:, a, d, e][None, :]
                        * reference_gamma[:, e, a, b][None, :]
                        + difference[:, a, d, e][None, :]
                        * difference[:, e, a, b][None, :]
                    )
            answer[:, :, b, d] = value
    return answer


def gr0_reference_balanced_ref1_residual_batch(
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
    """Evaluate the complete GR-0 REF1 rows in reference-covariant order."""

    base = _array("u", u, ndim=2)
    velocity = _array("p", p, ndim=2)
    gradient = _array("q", q, ndim=2)
    acceleration = _array("dtt", dtt, ndim=3)
    mixed = _array("p_r", p_r, ndim=2)
    radial_second = _array("q_r", q_r, ndim=2)
    radius = _array("radii", radii, ndim=1)
    points = radius.size
    if (
        base.shape != (points, FIELD_COUNT)
        or velocity.shape != base.shape
        or gradient.shape != base.shape
    ):
        raise ValueError("SRC3 GR-0 state shapes differ")
    if (
        acceleration.shape[1:] != base.shape
        or mixed.shape != base.shape
        or radial_second.shape != base.shape
    ):
        raise ValueError("SRC3 GR-0 second-jet shapes differ")
    if np.any(radius <= 0.0):
        raise ValueError("SRC3 GR-0 REF1 evaluator requires positive radii")
    mplanck = _finite("planck_mass", planck_mass, positive=True)
    mass = _finite("scalar_mass", scalar_mass, positive=True)
    quartic = _finite("quartic_coupling", quartic_coupling, positive=True)
    tilde_factor = _finite(
        "tilde_normal_factor", tilde_normal_factor, positive=True
    )
    hat_factor = _finite("hat_normal_factor", hat_normal_factor, positive=True)
    if tilde_factor <= 1.0 or hat_factor <= tilde_factor:
        raise ValueError("SRC3 auxiliary factors require 1<tilde<hat")

    metric, metric_derivative, metric_second = _metric_jet_from_adm(
        base,
        velocity,
        gradient,
        acceleration,
        mixed,
        radial_second,
    )
    inverse, inverse_derivative = _inverse_and_derivative(
        metric, metric_derivative
    )
    (
        reference_gamma,
        _reference_gamma_derivative,
        difference,
        difference_derivative,
    ) = _reference_covariant_connection_difference(
        metric,
        metric_derivative,
        metric_second,
        inverse,
        inverse_derivative,
        radius,
    )
    gamma = reference_gamma + difference
    ricci = _reference_covariant_ricci(
        reference_gamma,
        difference,
        difference_derivative,
    )
    scalar = np.einsum("nij,bnij->bn", inverse, ricci, optimize=True)
    einstein = ricci - metric[None, :, :, :] * scalar[:, :, None, None] / 2.0
    ricci_squared = np.einsum(
        "nac,nbd,snab,sncd->sn",
        inverse,
        inverse,
        ricci,
        ricci,
        optimize=True,
    )

    scalar_fields = (FIELD_INDEX["phi"], FIELD_INDEX["chi"])
    scalar_rows: list[np.ndarray] = []
    stresses: list[np.ndarray] = []
    primitive_second = _primitive_second_derivatives(
        acceleration, mixed, radial_second
    )
    for field_index in scalar_fields:
        first = np.zeros((points, SPACETIME_DIMENSION), dtype=np.float64)
        first[:, 0] = velocity[:, field_index]
        first[:, 1] = gradient[:, field_index]
        second_field = np.zeros(
            (acceleration.shape[0], points, SPACETIME_DIMENSION, SPACETIME_DIMENSION),
            dtype=np.float64,
        )
        second_field[:, :, :2, :2] = primitive_second[:, :, :, :, field_index]
        hessian = second_field.copy()
        for a in range(SPACETIME_DIMENSION):
            for b in range(SPACETIME_DIMENSION):
                for c in range(SPACETIME_DIMENSION):
                    hessian[:, :, a, b] -= (
                        gamma[:, c, a, b][None, :] * first[:, c][None, :]
                    )
        box = np.einsum("nij,bnij->bn", inverse, hessian, optimize=True)
        field_value = base[:, field_index]
        potential = (
            0.5 * mass**2 * field_value**2
            + 0.25 * quartic * field_value**4
            if field_index == FIELD_INDEX["phi"]
            else np.zeros(points, dtype=np.float64)
        )
        gradient_square = np.einsum(
            "nij,ni,nj->n", inverse, first, first, optimize=True
        )
        stress = (
            first[None, :, :, None] * first[None, :, None, :]
            - metric[None, :, :, :]
            * (
                gradient_square[None, :, None, None] / 2.0
                + potential[None, :, None, None]
            )
        )
        row = box - (
            mass**2 * field_value[None, :]
            + quartic * field_value[None, :] ** 3
            if field_index == FIELD_INDEX["phi"]
            else 0.0
        )
        scalar_rows.append(row)
        stresses.append(stress)

    original_metric = mplanck**2 * einstein - stresses[0] - stresses[1]
    tilde, tilde_derivative = _auxiliary_inverse_and_derivative(
        inverse, inverse_derivative, tilde_factor
    )
    constraint = -np.einsum(
        "nrs,nars->na", tilde, difference, optimize=True
    )
    partial_constraint = np.empty(
        (acceleration.shape[0], points, SPACETIME_DIMENSION, SPACETIME_DIMENSION),
        dtype=np.float64,
    )
    for direction in range(SPACETIME_DIMENSION):
        for a in range(SPACETIME_DIMENSION):
            value = np.zeros((acceleration.shape[0], points), dtype=np.float64)
            for rho in range(SPACETIME_DIMENSION):
                for sigma in range(SPACETIME_DIMENSION):
                    value -= (
                        tilde_derivative[:, direction, rho, sigma][None, :]
                        * difference[:, a, rho, sigma][None, :]
                        + tilde[:, rho, sigma][None, :]
                        * difference_derivative[
                            :, :, direction, a, rho, sigma
                        ]
                    )
            partial_constraint[:, :, direction, a] = value
    covariant_constraint_derivative = partial_constraint.copy()
    for direction in range(SPACETIME_DIMENSION):
        for a in range(SPACETIME_DIMENSION):
            for c in range(SPACETIME_DIMENSION):
                covariant_constraint_derivative[:, :, direction, a] += (
                    gamma[:, a, direction, c][None, :]
                    * constraint[:, c][None, :]
                )

    hat, _hat_derivative = _auxiliary_inverse_and_derivative(
        inverse, inverse_derivative, hat_factor
    )
    extension_up = np.zeros(
        (acceleration.shape[0], points, SPACETIME_DIMENSION, SPACETIME_DIMENSION),
        dtype=np.float64,
    )
    for mu in range(SPACETIME_DIMENSION):
        for nu in range(SPACETIME_DIMENSION):
            value = np.zeros((acceleration.shape[0], points), dtype=np.float64)
            for alpha in range(SPACETIME_DIMENSION):
                for beta in range(SPACETIME_DIMENSION):
                    projector = (
                        (1.0 if alpha == mu else 0.0) * hat[:, nu, beta]
                        + (1.0 if alpha == nu else 0.0) * hat[:, mu, beta]
                        - (1.0 if alpha == beta else 0.0) * hat[:, mu, nu]
                    ) / 2.0
                    value += (
                        projector[None, :]
                        * covariant_constraint_derivative[:, :, beta, alpha]
                    )
            extension_up[:, :, mu, nu] = mplanck**2 * value
    extension_down = np.einsum(
        "nam,snmp,npb->snab",
        metric,
        extension_up,
        metric,
        optimize=True,
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
    momentum = (
        original_metric[:, :, 0, 1]
        - shift[None, :] * original_metric[:, :, 1, 1]
    )
    return GR0ResidualBatch(
        full_residual=full,
        unredefined_metric_residual=original_metric,
        scalar_residual=np.stack(scalar_rows, axis=2),
        gauge_constraint=np.broadcast_to(
            constraint[None, :, :],
            (acceleration.shape[0], points, SPACETIME_DIMENSION),
        ),
        hamiltonian_constraint=hamiltonian,
        momentum_constraint=momentum,
        ricci_scalar=scalar,
        ricci_squared=ricci_squared,
    )


def diagnose_gr0_reference_balanced_accelerations(
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radii: object,
    *,
    raw_tolerance: Real = 1.0e-12,
    condition_number_maximum: Real = 1.0e10,
    maximum_refinement_iterations: int = 16,
) -> AffineRootDiagnosis:
    """Solve the SRC3 affine source and return its unchanged-gate evidence."""

    base = _array("u", u, ndim=2)
    velocity = _array("p", p, ndim=2)
    gradient = _array("q", q, ndim=2)
    mixed = _array("p_r", p_r, ndim=2)
    radial_second = _array("q_r", q_r, ndim=2)
    radius = _array("radii", radii, ndim=1)
    points, fields = base.shape
    if fields != FIELD_COUNT:
        raise ValueError("SRC3 u must contain the six canonical ADM fields")
    if (
        velocity.shape != base.shape
        or gradient.shape != base.shape
        or mixed.shape != base.shape
        or radial_second.shape != base.shape
        or radius.shape != (points,)
    ):
        raise ValueError("SRC3 GR-0 source shapes differ")
    if np.any(radius <= 0.0):
        raise ValueError("SRC3 radii must be positive")
    tolerance = _finite("raw_tolerance", raw_tolerance, positive=True)
    condition_limit = _finite(
        "condition_number_maximum", condition_number_maximum, positive=True
    )
    if (
        isinstance(maximum_refinement_iterations, bool)
        or not isinstance(maximum_refinement_iterations, int)
        or not 0 <= maximum_refinement_iterations <= 16
    ):
        raise ValueError("maximum_refinement_iterations must lie in [0,16]")

    seeds = np.zeros((FIELD_COUNT + 1, points, FIELD_COUNT), dtype=np.float64)
    for field in range(FIELD_COUNT):
        seeds[field + 1, :, field] = 1.0
    evaluated = gr0_reference_balanced_ref1_residual_batch(
        base,
        velocity,
        gradient,
        seeds,
        mixed,
        radial_second,
        radius,
    )
    constant = evaluated.full_residual[0]
    jacobian = np.moveaxis(
        evaluated.full_residual[1:] - constant[None, :, :],
        0,
        2,
    )
    try:
        inverse = np.linalg.inv(jacobian)
        acceleration = np.linalg.solve(
            jacobian, -constant[..., None]
        )[..., 0]
    except np.linalg.LinAlgError as error:
        raise ValueError("SRC3 GR-0 kinetic block is singular") from error
    condition = np.max(np.sum(np.abs(jacobian), axis=2), axis=1) * np.max(
        np.sum(np.abs(inverse), axis=2), axis=1
    )
    condition_maximum = float(np.max(condition, initial=0.0))
    if not np.all(np.isfinite(condition)) or condition_maximum >= condition_limit:
        raise ValueError("SRC3 GR-0 kinetic condition limit reached")
    if not np.all(np.isfinite(acceleration)):
        raise ValueError("SRC3 GR-0 acceleration root became nonfinite")

    iterations: list[ResidualIteration] = []
    residual: np.ndarray | None = None
    residual_norm: float | None = None
    previous_norm: float | None = None
    decreased = True
    floor = False
    for iteration in range(maximum_refinement_iterations + 1):
        evaluated_root = gr0_reference_balanced_ref1_residual_batch(
            base,
            velocity,
            gradient,
            acceleration[None, :, :],
            mixed,
            radial_second,
            radius,
        )
        residual = evaluated_root.full_residual[0]
        residual_norm = float(np.max(np.abs(residual), initial=0.0))
        maximum_index = np.unravel_index(
            int(np.argmax(np.abs(residual))), residual.shape
        )
        correction = np.linalg.solve(
            jacobian, -residual[..., None]
        )[..., 0]
        correction_norm = float(np.max(np.abs(correction), initial=0.0))
        acceleration_scale = max(
            1.0,
            float(np.max(np.abs(acceleration), initial=0.0)),
        )
        iterations.append(
            ResidualIteration(
                iteration=iteration,
                residual_infinity=residual_norm,
                acceleration_correction_infinity=correction_norm,
                relative_acceleration_correction=(
                    correction_norm / acceleration_scale
                ),
                maximum_residual_point_index=int(maximum_index[0]),
                maximum_residual_row_index=int(maximum_index[1]),
                maximum_residual_radius=float(radius[maximum_index[0]]),
            )
        )
        if residual_norm < tolerance or iteration == maximum_refinement_iterations:
            break
        candidate = acceleration + correction
        if not np.all(np.isfinite(candidate)):
            raise ValueError("SRC3 GR-0 refinement root became nonfinite")
        candidate_residual = gr0_reference_balanced_ref1_residual_batch(
            base,
            velocity,
            gradient,
            candidate[None, :, :],
            mixed,
            radial_second,
            radius,
        ).full_residual[0]
        candidate_norm = float(
            np.max(np.abs(candidate_residual), initial=0.0)
        )
        if previous_norm is not None and residual_norm >= previous_norm:
            decreased = False
        if candidate_norm >= residual_norm:
            floor = True
            break
        previous_norm = residual_norm
        acceleration = candidate

    if residual is None or residual_norm is None:
        raise RuntimeError("SRC3 root produced no residual")
    return AffineRootDiagnosis(
        accelerations=acceleration,
        residuals=residual,
        residual_infinity=residual_norm,
        raw_tolerance=tolerance,
        raw_gate_passed=residual_norm < tolerance,
        kinetic_condition_infinity_maximum=condition_maximum,
        residual_decreased_until_floor=decreased,
        roundoff_floor_reached=floor,
        iterations=tuple(iterations),
    )


__all__ = [
    "diagnose_gr0_reference_balanced_accelerations",
    "gr0_reference_balanced_ref1_residual_batch",
]
