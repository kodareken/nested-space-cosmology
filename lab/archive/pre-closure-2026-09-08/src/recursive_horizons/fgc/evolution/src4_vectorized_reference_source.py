"""Tensor-contracted binary64 evaluation of the unchanged SRC3 equations.

SRC3 removed a false source rejection by assembling the spherical connection
difference directly from ``h = g - g_bar``.  Its reference implementation
keeps the tensor indices as explicit Python loops so the derivation is easy to
audit.  Those loops dominate the cost of every Runge--Kutta stage.

This module evaluates the same reference-covariant identities with explicit
NumPy contractions.  It neither changes the continuum equations nor adds a
residual floor, tolerance, fallback precision, fitted coefficient, or root
branch.  The independent SRC3 implementation remains the differential oracle.
"""

from __future__ import annotations

from numbers import Real

import numpy as np

from .cal2_source_diagnosis import AffineRootDiagnosis, ResidualIteration
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
from .src3_reference_balanced_source import (
    SPACETIME_DIMENSION,
    _flat_reference_metric_jet,
    _inverse_and_derivative,
)


def _vectorized_reference_covariant_connection_difference(
    metric: np.ndarray,
    metric_derivative: np.ndarray,
    metric_second_derivative: np.ndarray,
    inverse: np.ndarray,
    inverse_derivative: np.ndarray,
    radii: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return ``Gamma_bar``, ``d Gamma_bar``, ``C``, and ``partial C``.

    Axis labels in the contractions follow the derivation directly:

    ``n`` point, ``s`` acceleration seed, ``e`` derivative direction, and
    ``a,b,c,d,f`` spacetime indices.  No contraction spans the point or seed
    axes, so the operation remains the same independent affine system at each
    radius.
    """

    batch, _points = metric_second_derivative.shape[:2]
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

    # bar_nabla_b h_dc
    covariant_first = (
        perturbation_first
        - np.einsum(
            "nfbd,nfc->nbdc",
            reference_gamma,
            perturbation,
            optimize=True,
        )
        - np.einsum(
            "nfbc,ndf->nbdc",
            reference_gamma,
            perturbation,
            optimize=True,
        )
    )

    # partial_e(bar_nabla_b h_dc).  The connection terms are independent of
    # the acceleration seed and are broadcast over its leading axis.
    covariant_first_derivative_correction = (
        np.einsum(
            "nefbd,nfc->nebdc",
            reference_gamma_derivative,
            perturbation,
            optimize=True,
        )
        + np.einsum(
            "nfbd,nefc->nebdc",
            reference_gamma,
            perturbation_first,
            optimize=True,
        )
        + np.einsum(
            "nefbc,ndf->nebdc",
            reference_gamma_derivative,
            perturbation,
            optimize=True,
        )
        + np.einsum(
            "nfbc,nedf->nebdc",
            reference_gamma,
            perturbation_first,
            optimize=True,
        )
    )
    covariant_first_derivative = (
        perturbation_second
        - covariant_first_derivative_correction[None, ...]
    )

    # B_dbc = bar_nabla_b h_dc + bar_nabla_c h_db - bar_nabla_d h_bc.
    bracket = (
        covariant_first.transpose(0, 2, 1, 3)
        + covariant_first.transpose(0, 2, 3, 1)
        - covariant_first
    )
    bracket_derivative = (
        covariant_first_derivative.transpose(0, 1, 2, 4, 3, 5)
        + covariant_first_derivative.transpose(0, 1, 2, 4, 5, 3)
        - covariant_first_derivative
    )

    difference = np.einsum(
        "nad,ndbc->nabc",
        inverse,
        bracket,
        optimize=True,
    ) / 2.0
    difference_derivative = (
        np.einsum(
            "nead,ndbc->neabc",
            inverse_derivative,
            bracket,
            optimize=True,
        )[None, ...]
        + np.einsum(
            "nad,snedbc->sneabc",
            inverse,
            bracket_derivative,
            optimize=True,
        )
    ) / 2.0
    return (
        reference_gamma,
        reference_gamma_derivative,
        difference,
        difference_derivative,
    )


def _vectorized_reference_covariant_ricci(
    reference_gamma: np.ndarray,
    difference: np.ndarray,
    difference_derivative: np.ndarray,
) -> np.ndarray:
    """Assemble Ricci from the flat reference and ``C`` by contraction."""

    derivative_part = np.einsum(
        "snaadb->snbd",
        difference_derivative,
        optimize=True,
    ) - np.einsum(
        "sndaab->snbd",
        difference_derivative,
        optimize=True,
    )
    reference_trace = np.einsum(
        "naae->ne",
        reference_gamma,
        optimize=True,
    )
    difference_trace = np.einsum(
        "naae->ne",
        difference,
        optimize=True,
    )
    algebraic_part = (
        np.einsum(
            "ne,nedb->nbd",
            reference_trace,
            difference,
            optimize=True,
        )
        + np.einsum(
            "ne,nedb->nbd",
            difference_trace,
            reference_gamma,
            optimize=True,
        )
        + np.einsum(
            "ne,nedb->nbd",
            difference_trace,
            difference,
            optimize=True,
        )
        - np.einsum(
            "nade,neab->nbd",
            reference_gamma,
            difference,
            optimize=True,
        )
        - np.einsum(
            "nade,neab->nbd",
            difference,
            reference_gamma,
            optimize=True,
        )
        - np.einsum(
            "nade,neab->nbd",
            difference,
            difference,
            optimize=True,
        )
    )
    return derivative_part + algebraic_part[None, ...]


def gr0_vectorized_reference_ref1_residual_batch(
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
    """Evaluate the complete unchanged GR-0 REF1 rows via tensor contractions."""

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
        raise ValueError("SRC4 GR-0 state shapes differ")
    if (
        acceleration.shape[1:] != base.shape
        or mixed.shape != base.shape
        or radial_second.shape != base.shape
    ):
        raise ValueError("SRC4 GR-0 second-jet shapes differ")
    if np.any(radius <= 0.0):
        raise ValueError("SRC4 GR-0 REF1 evaluator requires positive radii")
    mplanck = _finite("planck_mass", planck_mass, positive=True)
    mass = _finite("scalar_mass", scalar_mass, positive=True)
    quartic = _finite("quartic_coupling", quartic_coupling, positive=True)
    tilde_factor = _finite(
        "tilde_normal_factor", tilde_normal_factor, positive=True
    )
    hat_factor = _finite("hat_normal_factor", hat_normal_factor, positive=True)
    if tilde_factor <= 1.0 or hat_factor <= tilde_factor:
        raise ValueError("SRC4 auxiliary factors require 1<tilde<hat")

    metric, metric_derivative, metric_second = _metric_jet_from_adm(
        base,
        velocity,
        gradient,
        acceleration,
        mixed,
        radial_second,
    )
    inverse, inverse_derivative = _inverse_and_derivative(
        metric,
        metric_derivative,
    )
    (
        reference_gamma,
        _reference_gamma_derivative,
        difference,
        difference_derivative,
    ) = _vectorized_reference_covariant_connection_difference(
        metric,
        metric_derivative,
        metric_second,
        inverse,
        inverse_derivative,
        radius,
    )
    gamma = reference_gamma + difference
    ricci = _vectorized_reference_covariant_ricci(
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
        acceleration,
        mixed,
        radial_second,
    )
    for field_index in scalar_fields:
        first = np.zeros((points, SPACETIME_DIMENSION), dtype=np.float64)
        first[:, 0] = velocity[:, field_index]
        first[:, 1] = gradient[:, field_index]
        second_field = np.zeros(
            (
                acceleration.shape[0],
                points,
                SPACETIME_DIMENSION,
                SPACETIME_DIMENSION,
            ),
            dtype=np.float64,
        )
        second_field[:, :, :2, :2] = primitive_second[:, :, :, :, field_index]
        hessian = second_field - np.einsum(
            "ncab,nc->nab",
            gamma,
            first,
            optimize=True,
        )[None, ...]
        box = np.einsum("nij,bnij->bn", inverse, hessian, optimize=True)
        field_value = base[:, field_index]
        potential = (
            0.5 * mass**2 * field_value**2
            + 0.25 * quartic * field_value**4
            if field_index == FIELD_INDEX["phi"]
            else np.zeros(points, dtype=np.float64)
        )
        gradient_square = np.einsum(
            "nij,ni,nj->n",
            inverse,
            first,
            first,
            optimize=True,
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
        inverse,
        inverse_derivative,
        tilde_factor,
    )
    constraint = -np.einsum(
        "nrs,nars->na",
        tilde,
        difference,
        optimize=True,
    )
    partial_constraint = -(
        np.einsum(
            "nerq,narq->nea",
            tilde_derivative,
            difference,
            optimize=True,
        )[None, ...]
        + np.einsum(
            "nrq,xnearq->xnea",
            tilde,
            difference_derivative,
            optimize=True,
        )
    )
    covariant_constraint_derivative = partial_constraint + np.einsum(
        "nadc,nc->nda",
        gamma,
        constraint,
        optimize=True,
    )[None, ...]

    hat, _hat_derivative = _auxiliary_inverse_and_derivative(
        inverse,
        inverse_derivative,
        hat_factor,
    )
    identity = np.eye(SPACETIME_DIMENSION, dtype=np.float64)
    projector = (
        np.einsum("am,nvb->nmvab", identity, hat, optimize=True)
        + np.einsum("av,nmb->nmvab", identity, hat, optimize=True)
        - np.einsum("ab,nmv->nmvab", identity, hat, optimize=True)
    ) / 2.0
    extension_up = mplanck**2 * np.einsum(
        "nmvab,snba->snmv",
        projector,
        covariant_constraint_derivative,
        optimize=True,
    )
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


def diagnose_gr0_vectorized_reference_accelerations(
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
    """Solve the tensor-contracted affine source under SRC3's strict gate."""

    base = _array("u", u, ndim=2)
    velocity = _array("p", p, ndim=2)
    gradient = _array("q", q, ndim=2)
    mixed = _array("p_r", p_r, ndim=2)
    radial_second = _array("q_r", q_r, ndim=2)
    radius = _array("radii", radii, ndim=1)
    points, fields = base.shape
    if fields != FIELD_COUNT:
        raise ValueError("SRC4 u must contain the six canonical ADM fields")
    if (
        velocity.shape != base.shape
        or gradient.shape != base.shape
        or mixed.shape != base.shape
        or radial_second.shape != base.shape
        or radius.shape != (points,)
    ):
        raise ValueError("SRC4 GR-0 source shapes differ")
    if np.any(radius <= 0.0):
        raise ValueError("SRC4 radii must be positive")
    tolerance = _finite("raw_tolerance", raw_tolerance, positive=True)
    condition_limit = _finite(
        "condition_number_maximum",
        condition_number_maximum,
        positive=True,
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
    evaluated = gr0_vectorized_reference_ref1_residual_batch(
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
        acceleration = np.linalg.solve(jacobian, -constant[..., None])[..., 0]
    except np.linalg.LinAlgError as error:
        raise ValueError("SRC4 GR-0 kinetic block is singular") from error
    condition = np.max(np.sum(np.abs(jacobian), axis=2), axis=1) * np.max(
        np.sum(np.abs(inverse), axis=2),
        axis=1,
    )
    condition_maximum = float(np.max(condition, initial=0.0))
    if not np.all(np.isfinite(condition)) or condition_maximum >= condition_limit:
        raise ValueError("SRC4 GR-0 kinetic condition limit reached")
    if not np.all(np.isfinite(acceleration)):
        raise ValueError("SRC4 GR-0 acceleration root became nonfinite")

    iterations: list[ResidualIteration] = []
    residual: np.ndarray | None = None
    residual_norm: float | None = None
    previous_norm: float | None = None
    decreased = True
    floor = False
    for iteration in range(maximum_refinement_iterations + 1):
        evaluated_root = gr0_vectorized_reference_ref1_residual_batch(
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
            int(np.argmax(np.abs(residual))),
            residual.shape,
        )
        correction = np.linalg.solve(jacobian, -residual[..., None])[..., 0]
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
            raise ValueError("SRC4 GR-0 refinement root became nonfinite")
        candidate_residual = gr0_vectorized_reference_ref1_residual_batch(
            base,
            velocity,
            gradient,
            candidate[None, :, :],
            mixed,
            radial_second,
            radius,
        ).full_residual[0]
        candidate_norm = float(np.max(np.abs(candidate_residual), initial=0.0))
        if previous_norm is not None and residual_norm >= previous_norm:
            decreased = False
        if candidate_norm >= residual_norm:
            floor = True
            break
        previous_norm = residual_norm
        acceleration = candidate

    if residual is None or residual_norm is None:
        raise RuntimeError("SRC4 root produced no residual")
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
    "diagnose_gr0_vectorized_reference_accelerations",
    "gr0_vectorized_reference_ref1_residual_batch",
]
