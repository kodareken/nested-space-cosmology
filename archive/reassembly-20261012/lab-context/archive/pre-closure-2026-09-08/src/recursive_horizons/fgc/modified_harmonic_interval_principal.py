"""Exact interval REF1 first-order principal-block enclosures.

This module encloses only the complete REF1 implicit acceleration block and
its first-order radial derivative blocks on a compact jet box.  It supplies
no characteristic factorization, eigenvector, symmetrizer, constraint, or
evolution assertion.
"""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
from numbers import Integral
from typing import Any, Sequence

from .exact_interval import (
    Interval,
    interval,
    interval_matrix_infinity_row_sum_bound,
    interval_matrix_multiply,
    interval_matrix_subtract,
    point_matrix_interval_multiply,
)
from .exact_linear_algebra import matrix_det, matrix_inverse
from .interval_tangent import IntervalFirstTangent, primal_and_tangent
from .modified_harmonic_first_order import exact_full_residual_first_order_argument_jacobian
from .modified_harmonic_reference import ReferenceConnection, modified_harmonic_full_residuals
from .modified_harmonic_quantified_domain import (
    _regular_domain_ledger,
    interval_full_residual_acceleration_jacobian,
    krawczyk_acceleration_box,
)
from .spherical_reduction import BASE_FIELD_ORDER, Jet2, SphericalState


Q = Fraction
FIELD_COUNT = len(BASE_FIELD_ORDER)
INTERVAL_PRINCIPAL_STATE_ORDER = tuple(f"dt_{field}" for field in BASE_FIELD_ORDER) + tuple(
    f"dr_{field}" for field in BASE_FIELD_ORDER
)
INTERVAL_PRINCIPAL_ARGUMENT_ORDER = (
    tuple(f"p_t.{field}" for field in BASE_FIELD_ORDER)
    + tuple(f"p_r.{field}" for field in BASE_FIELD_ORDER)
    + tuple(f"q_r.{field}" for field in BASE_FIELD_ORDER)
)
_JET_SLOTS = ("value", "dt", "dr", "dtt", "dtr", "drr")


def _fraction(name: str, value: object, *, nonnegative: bool = False) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational")
    result = value if isinstance(value, Fraction) else Q(int(value))
    if nonnegative and result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _interval_state(
    center: SphericalState,
    *,
    parameter_half_width: Fraction,
    acceleration_half_width: Fraction,
    seed_slot: str | None = None,
    seed_field: str | None = None,
) -> SphericalState:
    if seed_slot is not None and seed_slot not in ("dtt", "dtr", "drr"):
        raise ValueError("interval principal seed slot must be dtt, dtr, or drr")
    if seed_slot is None and seed_field is not None:
        raise ValueError("interval principal seed field requires a seed slot")
    if seed_field is not None and seed_field not in BASE_FIELD_ORDER:
        raise ValueError("unknown interval principal seed field")
    fields: dict[str, Jet2] = {}
    for field in BASE_FIELD_ORDER:
        source = getattr(center, field)
        entries: dict[str, object] = {}
        for slot in _JET_SLOTS:
            width = acceleration_half_width if slot == "dtt" else parameter_half_width
            value = interval(getattr(source, slot) - width, getattr(source, slot) + width)
            if field == seed_field and slot == seed_slot:
                value = IntervalFirstTangent.seed(value)
            entries[slot] = value
        fields[field] = Jet2(**entries)
    return SphericalState(
        **fields,
        planck_mass=center.planck_mass,
        beta=center.beta,
        mu=center.mu,
        g4=center.g4,
        alpha=center.alpha,
        eta=center.eta,
        branch=center.branch,
    )


def _full_residual(
    state: SphericalState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: Fraction,
    tilde_normal_factor: Fraction,
    hat_normal_factor: Fraction,
) -> tuple[object, ...]:
    return tuple(
        modified_harmonic_full_residuals(
            state,
            reference=reference,
            coordinate_radius=coordinate_radius,
            tilde_normal_factor=tilde_normal_factor,
            hat_normal_factor=hat_normal_factor,
        )["full_residual_vector"]
    )


def interval_ref1_first_order_blocks(
    center: SphericalState,
    *,
    parameter_half_width: Fraction | int,
    acceleration_half_width: Fraction | int,
    reference: ReferenceConnection,
    coordinate_radius: Fraction | int,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
) -> dict[str, Any]:
    """Enclose complete REF1 ``J_a,J_p_r,J_q_r`` by interval forward AD.

    Every one of the 30 non-acceleration jet coordinates receives the same
    parameter half-width, while every one of the six dtt coordinates receives
    the acceleration half-width.  Each derivative column is obtained from a
    distinct ``IntervalFirstTangent`` seed of the complete residual evaluator.
    """
    if not isinstance(center, SphericalState):
        raise TypeError("interval principal center must be a SphericalState")
    if not isinstance(reference, ReferenceConnection):
        raise TypeError("interval principal reference must be a ReferenceConnection")
    dz = _fraction("parameter_half_width", parameter_half_width, nonnegative=True)
    da = _fraction("acceleration_half_width", acceleration_half_width, nonnegative=True)
    radius = _fraction("coordinate_radius", coordinate_radius)
    tilde = _fraction("tilde_normal_factor", tilde_normal_factor)
    hat = _fraction("hat_normal_factor", hat_normal_factor)
    if radius <= reference.radial_domain_minimum:
        raise ValueError("coordinate radius lies outside the reference annulus")
    if not 1 < tilde < hat:
        raise ValueError("interval principal requires 1 < tilde factor < hat factor")

    columns: dict[str, list[tuple[Interval, ...]]] = {"p_t": [], "p_r": [], "q_r": []}
    primal_columns: dict[str, list[tuple[Interval, ...]]] = {"p_t": [], "p_r": [], "q_r": []}
    for group, slot in (("p_t", "dtt"), ("p_r", "dtr"), ("q_r", "drr")):
        for field in BASE_FIELD_ORDER:
            seeded = _interval_state(
                center,
                parameter_half_width=dz,
                acceleration_half_width=da,
                seed_slot=slot,
                seed_field=field,
            )
            parts = tuple(
                primal_and_tangent(value)
                for value in _full_residual(
                    seeded,
                    reference=reference,
                    coordinate_radius=radius,
                    tilde_normal_factor=tilde,
                    hat_normal_factor=hat,
                )
            )
            primal_columns[group].append(tuple(item[0] for item in parts))
            columns[group].append(tuple(item[1] for item in parts))
    blocks = {
        group: tuple(
            tuple(columns[group][column][row] for column in range(FIELD_COUNT))
            for row in range(FIELD_COUNT)
        )
        for group in ("p_t", "p_r", "q_r")
    }
    primals_match = all(
        primal_columns[group][column] == primal_columns["p_t"][0]
        for group in ("p_t", "p_r", "q_r")
        for column in range(FIELD_COUNT)
    )
    if not primals_match:
        raise ValueError("interval tangent primal residual changed across derivative seeds")
    return {
        "equation_order": BASE_FIELD_ORDER,
        "argument_order": INTERVAL_PRINCIPAL_ARGUMENT_ORDER,
        "state_order": INTERVAL_PRINCIPAL_STATE_ORDER,
        "parameter_half_width": dz,
        "acceleration_half_width": da,
        "all_30_parameter_axes_have_nonzero_width": dz > 0,
        "all_6_acceleration_axes_have_nonzero_width": da > 0,
        "blocks": blocks,
        "seeded_primal_residual_boxes": primal_columns,
        "seeded_primals_reproduce_same_complete_residual_box": primals_match,
        "differentiation_method": "complete_REF1_exact_rational_IntervalFirstTangent_forward_AD",
    }


def _identity_interval(size: int) -> tuple[tuple[Interval, ...], ...]:
    return tuple(tuple(interval(int(row == column)) for column in range(size)) for row in range(size))


def _point_matrix_infinity_norm(value: Sequence[Sequence[Fraction]]) -> Fraction:
    return max(sum(abs(entry) for entry in row) for row in value)


def kinetic_neumann_inverse_enclosure(
    exact_center_kinetic: Sequence[Sequence[Fraction]],
    kinetic_box: Sequence[Sequence[Interval]],
) -> dict[str, Any]:
    """Enclose a kinetic inverse via a fixed exact-center Neumann bound."""
    if matrix_det(exact_center_kinetic) == 0:
        raise ValueError("exact center kinetic block is singular")
    center_inverse = matrix_inverse(exact_center_kinetic)
    preconditioned = point_matrix_interval_multiply(center_inverse, kinetic_box)
    defect = interval_matrix_subtract(_identity_interval(FIELD_COUNT), preconditioned)
    rho = interval_matrix_infinity_row_sum_bound(defect)
    if rho >= 1:
        raise ValueError("kinetic Neumann enclosure does not contract")
    inverse_norm = _point_matrix_infinity_norm(center_inverse)
    inverse_delta_bound = rho * inverse_norm / (1 - rho)
    inverse_box = tuple(
        tuple(interval(entry - inverse_delta_bound, entry + inverse_delta_bound) for entry in row)
        for row in center_inverse
    )
    return {
        "exact_center_inverse": center_inverse,
        "preconditioned_kinetic_box": preconditioned,
        "identity_minus_preconditioned_kinetic_box": defect,
        "neumann_rho_infinity_upper_bound": rho,
        "exact_center_inverse_infinity_norm": inverse_norm,
        "inverse_difference_infinity_upper_bound": inverse_delta_bound,
        "inverse_box": inverse_box,
        "kinetic_inverse_regular_over_entire_box": True,
    }


def compact_ref1_interval_principal_certificate(
    center: SphericalState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: Fraction | int,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
    parameter_half_width: Fraction = Q(1, 2**60),
    acceleration_half_width: Fraction = Q(1, 2**50),
) -> dict[str, Any]:
    """Return the compact nonflat REF1 branch-principal enclosure only."""
    dz = _fraction("parameter_half_width", parameter_half_width, nonnegative=True)
    da = _fraction("acceleration_half_width", acceleration_half_width, nonnegative=True)
    if dz == 0 or da == 0:
        raise ValueError("compact certificate requires nonzero parameter and acceleration widths")
    blocks = interval_ref1_first_order_blocks(
        center,
        parameter_half_width=dz,
        acceleration_half_width=da,
        reference=reference,
        coordinate_radius=coordinate_radius,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    exact = exact_full_residual_first_order_argument_jacobian(
        center,
        reference=reference,
        coordinate_radius=coordinate_radius,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )["blocks"]
    neumann = kinetic_neumann_inverse_enclosure(exact["p_t"], blocks["blocks"]["p_t"])
    krawczyk_input = interval_full_residual_acceleration_jacobian(
        center,
        parameter_half_width=dz,
        acceleration_half_width=da,
        reference=reference,
        coordinate_radius=coordinate_radius,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    krawczyk = krawczyk_acceleration_box(
        inverse_center_jacobian=neumann["exact_center_inverse"],
        residual_at_center_acceleration=krawczyk_input["residual_at_center_acceleration_parameter_box"],
        acceleration_jacobian_box=krawczyk_input["acceleration_jacobian_box"],
        acceleration_half_width=da,
    )
    radial_box = tuple(
        tuple(blocks["blocks"]["p_r"][row] + blocks["blocks"]["q_r"][row])
        for row in range(FIELD_COUNT)
    )
    top = interval_matrix_multiply(neumann["inverse_box"], radial_box)
    zero = interval(0)
    lower = tuple(
        tuple(interval(-int(row == column)) for column in range(FIELD_COUNT)) + tuple(zero for _ in range(FIELD_COUNT))
        for row in range(FIELD_COUNT)
    )
    branch = tuple(top) + lower
    regularity = _regular_domain_ledger(
        center,
        parameter_half_width=dz,
        reference=reference,
        coordinate_radius=Q(coordinate_radius),
        tilde_normal_factor=Q(tilde_normal_factor),
        hat_normal_factor=Q(hat_normal_factor),
    )
    return {
        "classification": "exact_rational_compact_REF1_implicit_branch_principal_enclosure_not_a_characteristic_or_hyperbolicity_certificate",
        "orders": {
            "equation_order": BASE_FIELD_ORDER,
            "acceleration_block": tuple(f"p_t.{field}" for field in BASE_FIELD_ORDER),
            "radial_blocks": tuple(f"p_r.{field}" for field in BASE_FIELD_ORDER) + tuple(f"q_r.{field}" for field in BASE_FIELD_ORDER),
            "first_order_state": INTERVAL_PRINCIPAL_STATE_ORDER,
        },
        "box": blocks,
        "regularity": regularity,
        "kinetic_neumann_inverse": neumann,
        "complete_REF1_acceleration_Krawczyk": krawczyk,
        "solved_branch_principal_box": branch,
        "nonclaims": {
            "characteristic_factorization_proven": False,
            "real_characteristic_roots_proven": False,
            "uniform_complete_eigenframe_proven": False,
            "uniform_symmetrizer_proven": False,
            "uniform_strong_hyperbolicity_proven": False,
            "constraint_propagation_proven": False,
            "evolution_authorized": False,
            "retained_eft_domain_proven": False,
        },
    }


def exact_comp1_interval_principal_certificate(
    unresolved_comp1_state: SphericalState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: Fraction | int,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
    parameter_half_width: Fraction = Q(1, 2**60),
    acceleration_half_width: Fraction = Q(1, 2**50),
) -> dict[str, Any]:
    """Derive the exact COMP1 root, then enclose its compact REF1 branch.

    The root derivation is delegated to MODE1's exact complete-residual check;
    this wrapper deliberately returns its serializable root facts separately
    from the local interval principal enclosure.
    """
    from .modified_harmonic_modes import exact_comp1_acceleration_root

    root = exact_comp1_acceleration_root(
        unresolved_comp1_state,
        reference=reference,
        coordinate_radius=coordinate_radius,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    interval_certificate = compact_ref1_interval_principal_certificate(
        root["solved_state"],
        reference=reference,
        coordinate_radius=coordinate_radius,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
        parameter_half_width=parameter_half_width,
        acceleration_half_width=acceleration_half_width,
    )
    return {
        "classification": interval_certificate["classification"],
        "exact_comp1_root": {
            key: value for key, value in root.items() if key != "solved_state"
        },
        "interval_principal": interval_certificate,
    }
