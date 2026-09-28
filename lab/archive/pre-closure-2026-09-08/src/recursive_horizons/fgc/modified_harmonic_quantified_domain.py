"""Quantified rational-interval branch box for the complete REF1 residual.

QIFT1 lifts the already frozen REF1 evaluator onto exact rational interval
arithmetic.  Around the IMP1 flat root it applies one fixed exact inverse
Jacobian as a preconditioner and verifies a uniform contraction/Krawczyk
inclusion for every point of a declared full-dimensional local first-jet box.

The result is finite-dimensional and algebraic.  It does not close the
metric-derived constraints, construct a full quasilinear evolution system,
prove strong hyperbolicity, define an EFT cutoff, or authorize evolution.
"""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
from numbers import Integral
from typing import Any, Mapping, Sequence

from .exact_interval import (
    Interval,
    IntervalMatrix,
    IntervalVector,
    interval,
    interval_matrix_determinant,
    interval_matrix_infinity_row_sum_bound,
    interval_matrix_subtract,
    interval_matrix_vector,
    point_matrix_interval_multiply,
)
from .exact_linear_algebra import Matrix, matrix_det, matrix_inverse, matrix_rank
from .interval_tangent import IntervalFirstTangent, primal_and_tangent
from .modified_harmonic import auxiliary_inverse_metric
from .modified_harmonic_implicit import (
    EXPECTED_FGCQR_PARAMETERS,
    exact_full_residual_acceleration_jacobian,
)
from .modified_harmonic_reference import modified_harmonic_full_residuals
from .reference_connection import ReferenceConnection
from .spherical_reduction import (
    BASE_FIELD_ORDER,
    Jet2,
    SphericalState,
    state_from_generalized_adm_pg_fixture,
)


Q = Fraction
FIELD_COUNT = len(BASE_FIELD_ORDER)
QIFT1_PARAMETER_GROUPS = ("u", "p", "q", "p_r", "q_r")
QIFT1_PARAMETER_ORDER = tuple(
    f"{group}.{field}"
    for group in QIFT1_PARAMETER_GROUPS
    for field in BASE_FIELD_ORDER
)
QIFT1_ACCELERATION_ORDER = tuple(f"p_t.{field}" for field in BASE_FIELD_ORDER)

_GROUP_TO_SLOT = {
    "u": "value",
    "p": "dt",
    "q": "dr",
    "p_r": "dtr",
    "q_r": "drr",
}


def _exact_fraction(name: str, value: object, *, positive: bool = False) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational value")
    result = value if isinstance(value, Fraction) else Q(int(value))
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _zero_vector(size: int) -> tuple[Fraction, ...]:
    return tuple(Q(0) for _ in range(size))


def _identity_interval(size: int) -> IntervalMatrix:
    return tuple(
        tuple(interval(int(row == column)) for column in range(size))
        for row in range(size)
    )


def _interval_centered(value: Fraction, half_width: Fraction) -> Interval:
    return interval(value - half_width, value + half_width)


def _interval_state(
    center: SphericalState,
    *,
    parameter_half_width: Fraction,
    acceleration_half_width: Fraction,
    seeded_acceleration_field: str | None = None,
) -> SphericalState:
    """Return the declared 30-parameter box and six-acceleration enclosure."""

    if seeded_acceleration_field is not None and seeded_acceleration_field not in BASE_FIELD_ORDER:
        raise ValueError("unknown acceleration seed field")
    fields: dict[str, Jet2] = {}
    for field in BASE_FIELD_ORDER:
        source = getattr(center, field)
        slots = {
            slot: _interval_centered(getattr(source, slot), parameter_half_width)
            for slot in _GROUP_TO_SLOT.values()
        }
        acceleration: Interval | IntervalFirstTangent = interval(
            source.dtt - acceleration_half_width,
            source.dtt + acceleration_half_width,
        )
        if field == seeded_acceleration_field:
            acceleration = IntervalFirstTangent.seed(acceleration)
        fields[field] = Jet2(dtt=acceleration, **slots)
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


def interval_full_residual_acceleration_jacobian(
    center: SphericalState,
    *,
    parameter_half_width: int | Fraction,
    acceleration_half_width: int | Fraction,
    reference: ReferenceConnection,
    coordinate_radius: int | Fraction,
    tilde_normal_factor: int | Fraction,
    hat_normal_factor: int | Fraction,
) -> dict[str, Any]:
    """Enclose ``R(a0,Z)`` and ``D_a R(A,Z)`` by interval forward AD."""

    if not isinstance(center, SphericalState):
        raise TypeError("center must be a SphericalState")
    if not isinstance(reference, ReferenceConnection):
        raise TypeError("reference must be a ReferenceConnection")
    parameter_radius = _exact_fraction(
        "parameter_half_width", parameter_half_width, positive=True
    )
    acceleration_radius = _exact_fraction(
        "acceleration_half_width", acceleration_half_width, positive=True
    )
    radius = _exact_fraction("coordinate_radius", coordinate_radius, positive=True)
    tilde = _exact_fraction("tilde_normal_factor", tilde_normal_factor, positive=True)
    hat = _exact_fraction("hat_normal_factor", hat_normal_factor, positive=True)
    if not 1 < tilde < hat:
        raise ValueError("QIFT1 requires 1 < tilde factor < hat factor")
    if radius <= reference.radial_domain_minimum:
        raise ValueError("QIFT1 coordinate radius lies outside the reference annulus")

    parameter_state = _interval_state(
        center,
        parameter_half_width=parameter_radius,
        acceleration_half_width=Q(0),
    )
    center_acceleration_residual = tuple(
        primal_and_tangent(value)[0]
        for value in _full_residual(
            parameter_state,
            reference=reference,
            coordinate_radius=radius,
            tilde_normal_factor=tilde,
            hat_normal_factor=hat,
        )
    )

    columns: list[IntervalVector] = []
    seeded_primals: list[IntervalVector] = []
    for field in BASE_FIELD_ORDER:
        state = _interval_state(
            center,
            parameter_half_width=parameter_radius,
            acceleration_half_width=acceleration_radius,
            seeded_acceleration_field=field,
        )
        parts = tuple(
            primal_and_tangent(value)
            for value in _full_residual(
                state,
                reference=reference,
                coordinate_radius=radius,
                tilde_normal_factor=tilde,
                hat_normal_factor=hat,
            )
        )
        seeded_primals.append(tuple(value[0] for value in parts))
        columns.append(tuple(value[1] for value in parts))

    jacobian = tuple(
        tuple(columns[column][row] for column in range(FIELD_COUNT))
        for row in range(FIELD_COUNT)
    )
    return {
        "parameter_order": QIFT1_PARAMETER_ORDER,
        "acceleration_order": QIFT1_ACCELERATION_ORDER,
        "parameter_dimension": len(QIFT1_PARAMETER_ORDER),
        "acceleration_dimension": FIELD_COUNT,
        "parameter_half_width": parameter_radius,
        "acceleration_half_width": acceleration_radius,
        "residual_at_center_acceleration_parameter_box": center_acceleration_residual,
        "seeded_primal_residual_boxes": tuple(seeded_primals),
        "acceleration_jacobian_box": jacobian,
        "differentiation_method": "exact_rational_interval_first_tangent_forward_ad",
        "all_parameter_axes_have_nonzero_width": parameter_radius > 0,
    }


def _regular_domain_ledger(
    center: SphericalState,
    *,
    parameter_half_width: Fraction,
    reference: ReferenceConnection,
    coordinate_radius: Fraction,
    tilde_normal_factor: Fraction,
    hat_normal_factor: Fraction,
) -> dict[str, Any]:
    state = _interval_state(
        center,
        parameter_half_width=parameter_half_width,
        acceleration_half_width=Q(0),
    )
    h_tt = state.h_tt.value
    h_tr = state.h_tr.value
    h_rr = state.h_rr.value
    areal_radius = state.areal_radius.value
    phi = state.phi.value
    base_determinant = h_tt * h_rr - h_tr**2
    inverse_tt = h_rr / base_determinant
    radius_squared = areal_radius**2
    physical_inverse = (
        (h_rr / base_determinant, -h_tr / base_determinant, interval(0), interval(0)),
        (-h_tr / base_determinant, h_tt / base_determinant, interval(0), interval(0)),
        (interval(0), interval(0), 1 / radius_squared, interval(0)),
        (interval(0), interval(0), interval(0), 1 / radius_squared),
    )
    tilde_inverse = auxiliary_inverse_metric(physical_inverse, tilde_normal_factor)
    hat_inverse = auxiliary_inverse_metric(physical_inverse, hat_normal_factor)
    tilde_determinant = interval_matrix_determinant(tilde_inverse)
    hat_determinant = interval_matrix_determinant(hat_inverse)
    effective_planck = center.planck_mass**2 + center.beta * phi**2
    reference_margin = coordinate_radius - reference.radial_domain_minimum
    margins = {
        "areal_radius_positive_lower_margin": areal_radius.lower,
        "base_metric_determinant_negative_lower_margin": -base_determinant.upper,
        "coordinate_time_inverse_metric_negative_lower_margin": -inverse_tt.upper,
        "effective_planck_coefficient_positive_lower_margin": effective_planck.lower,
        "tilde_auxiliary_determinant_negative_lower_margin": -tilde_determinant.upper,
        "hat_auxiliary_determinant_negative_lower_margin": -hat_determinant.upper,
        "reference_annulus_coordinate_margin": reference_margin,
    }
    if any(value <= 0 for value in margins.values()):
        raise ValueError("QIFT1 regular-domain ledger lacks a strict positive margin")
    return {
        "base_metric_determinant": base_determinant,
        "coordinate_time_inverse_metric": inverse_tt,
        "areal_radius": areal_radius,
        "effective_planck_coefficient": effective_planck,
        "tilde_auxiliary_inverse_metric_determinant": tilde_determinant,
        "hat_auxiliary_inverse_metric_determinant": hat_determinant,
        "strict_positive_margins": margins,
        "all_declared_denominators_and_sign_branches_regular": True,
    }


def krawczyk_acceleration_box(
    *,
    residual_at_center_acceleration: Sequence[Interval],
    acceleration_jacobian_box: Sequence[Sequence[Interval]],
    inverse_center_jacobian: Matrix,
    acceleration_half_width: int | Fraction,
) -> dict[str, Any]:
    """Apply the fixed-center contraction/Krawczyk inclusion exactly."""

    acceleration_radius = _exact_fraction(
        "acceleration_half_width", acceleration_half_width, positive=True
    )
    if len(residual_at_center_acceleration) != FIELD_COUNT:
        raise ValueError("QIFT1 residual box must have six rows")
    if len(acceleration_jacobian_box) != FIELD_COUNT or any(
        len(row) != FIELD_COUNT for row in acceleration_jacobian_box
    ):
        raise ValueError("QIFT1 acceleration Jacobian box must be six by six")
    if len(inverse_center_jacobian) != FIELD_COUNT or any(
        len(row) != FIELD_COUNT for row in inverse_center_jacobian
    ):
        raise ValueError("QIFT1 center inverse Jacobian must be six by six")

    preconditioned_jacobian = point_matrix_interval_multiply(
        inverse_center_jacobian, acceleration_jacobian_box
    )
    contraction_matrix = interval_matrix_subtract(
        _identity_interval(FIELD_COUNT), preconditioned_jacobian
    )
    contraction_bound = interval_matrix_infinity_row_sum_bound(contraction_matrix)
    center_correction = tuple(
        -value
        for value in interval_matrix_vector(
            inverse_center_jacobian, residual_at_center_acceleration
        )
    )
    acceleration_box = tuple(
        interval(-acceleration_radius, acceleration_radius)
        for _ in range(FIELD_COUNT)
    )
    linear_image = interval_matrix_vector(contraction_matrix, acceleration_box)
    krawczyk_box = tuple(
        center_correction[index] + linear_image[index]
        for index in range(FIELD_COUNT)
    )
    inclusion_margins = tuple(
        min(
            krawczyk_box[index].lower - acceleration_box[index].lower,
            acceleration_box[index].upper - krawczyk_box[index].upper,
        )
        for index in range(FIELD_COUNT)
    )
    strict_inclusion = all(
        krawczyk_box[index].strictly_inside(acceleration_box[index])
        for index in range(FIELD_COUNT)
    )
    contraction = contraction_bound < 1
    if not contraction:
        raise ValueError("QIFT1 contraction infinity-norm bound is not below one")
    if not strict_inclusion or min(inclusion_margins) <= 0:
        raise ValueError("QIFT1 Krawczyk image is not strictly inside the acceleration box")
    inverse_norm = max(
        sum(abs(entry) for entry in row) for row in inverse_center_jacobian
    )
    uniform_inverse_bound = inverse_norm / (1 - contraction_bound)
    return {
        "norm": "unscaled_componentwise_infinity_norm_in_frozen_field_order",
        "theorem_route": "uniform_banach_contraction_with_krawczyk_box_inclusion",
        "inverse_center_jacobian": inverse_center_jacobian,
        "inverse_center_jacobian_infinity_norm": inverse_norm,
        "preconditioned_acceleration_jacobian_box": preconditioned_jacobian,
        "identity_minus_preconditioned_jacobian_box": contraction_matrix,
        "contraction_infinity_norm_upper_bound": contraction_bound,
        "contraction_strictly_less_than_one": contraction,
        "center_newton_correction_box": center_correction,
        "acceleration_box": acceleration_box,
        "linear_contraction_image_box": linear_image,
        "krawczyk_image_box": krawczyk_box,
        "strict_interior_inclusion_margins": inclusion_margins,
        "minimum_strict_interior_inclusion_margin": min(inclusion_margins),
        "krawczyk_image_strictly_inside_acceleration_box": strict_inclusion,
        "uniform_acceleration_jacobian_inverse_infinity_norm_upper_bound": uniform_inverse_bound,
        "unique_acceleration_root_for_every_declared_parameter_point": True,
        "fixed_point_iterations_converge_from_every_point_in_declared_acceleration_box": True,
    }


def required_qift1_nonclaims() -> dict[str, bool]:
    return {
        name: False
        for name in (
            "exact_closed_form_nonlinear_acceleration_map_derived",
            "full_quasilinear_first_order_evolution_system_derived",
            "metric_derived_gauge_constraint_propagation_proven",
            "complete_reduction_constraint_system_propagation_proven",
            "physical_hamiltonian_momentum_constraints_derived",
            "physical_initial_constraints_solved",
            "constraint_preserving_initial_boundary_value_problem_proven",
            "uniform_open_domain_symmetrizer_proven",
            "uniform_radial_strong_hyperbolicity_box_proven",
            "retained_eft_cutoff_and_omitted_operator_envelope_defined",
            "open_retained_eft_domain_proven",
            "boundary_or_regular_center_system_derived",
            "evolution_authorized",
            "collapse_solution_derived",
            "metric_null_affine_defocusing_derived",
            "singularity_resolution_derived",
        )
    }


def quantified_implicit_branch_certificate(
    configuration: Mapping[str, Any],
) -> dict[str, Any]:
    """Build the exact QIFT1 full-dimensional local branch certificate."""

    expected_keys = {
        "flat_fixture",
        "reference",
        "coordinate_radius",
        "tilde_normal_factor",
        "hat_normal_factor",
        "parameter_half_width",
        "acceleration_half_width",
    }
    if not isinstance(configuration, Mapping) or set(configuration) != expected_keys:
        raise ValueError("QIFT1 configuration adapter has unexpected keys")
    fixture = configuration["flat_fixture"]
    reference = configuration["reference"]
    if not isinstance(fixture, Mapping):
        raise ValueError("QIFT1 flat fixture must be a mapping")
    if not isinstance(reference, ReferenceConnection):
        raise ValueError("QIFT1 reference must be a ReferenceConnection")
    radius = _exact_fraction(
        "coordinate_radius", configuration["coordinate_radius"], positive=True
    )
    tilde = _exact_fraction(
        "tilde_normal_factor", configuration["tilde_normal_factor"], positive=True
    )
    hat = _exact_fraction(
        "hat_normal_factor", configuration["hat_normal_factor"], positive=True
    )
    parameter_half_width = _exact_fraction(
        "parameter_half_width", configuration["parameter_half_width"], positive=True
    )
    acceleration_half_width = _exact_fraction(
        "acceleration_half_width",
        configuration["acceleration_half_width"],
        positive=True,
    )
    if not 1 < tilde < hat:
        raise ValueError("QIFT1 requires 1 < tilde factor < hat factor")
    if radius <= reference.radial_domain_minimum:
        raise ValueError("QIFT1 root must lie inside the reference annulus")

    center = state_from_generalized_adm_pg_fixture(fixture)
    parameters = {
        "planck_mass": center.planck_mass,
        "beta": center.beta,
        "mu": center.mu,
        "g4": center.g4,
        "alpha": center.alpha,
        "eta": center.eta,
    }
    if parameters != EXPECTED_FGCQR_PARAMETERS:
        raise ValueError("QIFT1 requires the frozen IMP1 FGC-QR parameters")
    exact_root = exact_full_residual_acceleration_jacobian(
        center,
        reference=reference,
        coordinate_radius=radius,
        tilde_normal_factor=tilde,
        hat_normal_factor=hat,
    )
    if tuple(exact_root["base_full_residual_vector"]) != _zero_vector(FIELD_COUNT):
        raise ValueError("QIFT1 center is not the exact IMP1 residual root")
    if tuple(exact_root["base_acceleration_vector"]) != _zero_vector(FIELD_COUNT):
        raise ValueError("QIFT1 center acceleration is not the zero IMP1 root")
    center_determinant = matrix_det(exact_root["jacobian"])
    center_rank = matrix_rank(exact_root["jacobian"])
    if center_determinant != Q(-165888) or center_rank != FIELD_COUNT:
        raise ValueError("QIFT1 center kinetic block differs from IMP1")
    inverse = matrix_inverse(exact_root["jacobian"])

    interval_data = interval_full_residual_acceleration_jacobian(
        center,
        parameter_half_width=parameter_half_width,
        acceleration_half_width=acceleration_half_width,
        reference=reference,
        coordinate_radius=radius,
        tilde_normal_factor=tilde,
        hat_normal_factor=hat,
    )
    jacobian_box = interval_data["acceleration_jacobian_box"]
    center_jacobian_contained = all(
        jacobian_box[row][column].lower
        <= exact_root["jacobian"][row][column]
        <= jacobian_box[row][column].upper
        for row in range(FIELD_COUNT)
        for column in range(FIELD_COUNT)
    )
    if not center_jacobian_contained:
        raise ValueError("QIFT1 interval Jacobian does not enclose the IMP1 center")
    regular = _regular_domain_ledger(
        center,
        parameter_half_width=parameter_half_width,
        reference=reference,
        coordinate_radius=radius,
        tilde_normal_factor=tilde,
        hat_normal_factor=hat,
    )
    krawczyk = krawczyk_acceleration_box(
        residual_at_center_acceleration=interval_data[
            "residual_at_center_acceleration_parameter_box"
        ],
        acceleration_jacobian_box=jacobian_box,
        inverse_center_jacobian=inverse,
        acceleration_half_width=acceleration_half_width,
    )
    verified = {
        "frozen_parameter_and_acceleration_orders_exact": len(QIFT1_PARAMETER_ORDER)
        == 30
        and len(QIFT1_ACCELERATION_ORDER) == FIELD_COUNT,
        "all_thirty_parameter_axes_have_nonzero_width": interval_data[
            "all_parameter_axes_have_nonzero_width"
        ],
        "flat_center_complete_ref1_residual_zero_exact": tuple(
            exact_root["base_full_residual_vector"]
        )
        == _zero_vector(FIELD_COUNT),
        "flat_center_acceleration_jacobian_matches_imp1_exact": exact_root[
            "seeded_primals_reproduce_base_residual"
        ]
        and center_determinant == Q(-165888)
        and center_rank == FIELD_COUNT,
        "interval_acceleration_jacobian_contains_exact_center": center_jacobian_contained,
        "all_regular_chart_margins_strictly_positive": regular[
            "all_declared_denominators_and_sign_branches_regular"
        ],
        "contraction_infinity_norm_strictly_below_one": krawczyk[
            "contraction_strictly_less_than_one"
        ],
        "krawczyk_image_strictly_inside_acceleration_box": krawczyk[
            "krawczyk_image_strictly_inside_acceleration_box"
        ],
        "uniform_unique_acceleration_root_for_every_parameter_point": krawczyk[
            "unique_acceleration_root_for_every_declared_parameter_point"
        ],
    }
    if not all(verified.values()):
        failed = sorted(name for name, value in verified.items() if not value)
        raise ValueError(f"QIFT1 exact gates failed: {failed}")
    nonclaims = required_qift1_nonclaims()
    return {
        "artifact_id": "FGC-1-HYP1-DOM1-QIFT1",
        "classification": "exact_rational_interval_uniform_contraction_certificate_for_a_full_dimensional_local_ref1_implicit_acceleration_branch_not_a_pde_health_or_eft_theorem",
        "formulation": {
            "residual": "complete_REF1_residual_R(a;z)=0",
            "parameter_definition": "z=(u,p,q,p_r,q_r)_in_frozen_base_field_order",
            "acceleration_definition": "a=p_t_in_frozen_base_field_order",
            "parameter_order": QIFT1_PARAMETER_ORDER,
            "acceleration_order": QIFT1_ACCELERATION_ORDER,
            "parameter_dimension": len(QIFT1_PARAMETER_ORDER),
            "acceleration_dimension": FIELD_COUNT,
            "coordinate_radius_is_fixed_not_a_box_axis": True,
            "reference_connection": reference.reference_id,
            "coordinate_radius": radius,
            "tilde_normal_factor": tilde,
            "hat_normal_factor": hat,
        },
        "parameter_box": {
            "center": tuple(
                getattr(getattr(center, field), _GROUP_TO_SLOT[group])
                for group in QIFT1_PARAMETER_GROUPS
                for field in BASE_FIELD_ORDER
            ),
            "uniform_half_width": parameter_half_width,
            "dimension": len(QIFT1_PARAMETER_ORDER),
            "all_axes_nonzero_width": True,
            "closed_box_has_nonempty_open_interior_in_local_jet_space": True,
        },
        "acceleration_box": {
            "center": tuple(getattr(center, field).dtt for field in BASE_FIELD_ORDER),
            "uniform_half_width": acceleration_half_width,
            "dimension": FIELD_COUNT,
        },
        "regular_domain_ledger": regular,
        "interval_residual_and_jacobian": interval_data,
        "krawczyk_contraction_certificate": krawczyk,
        "verified_exact_checks": verified,
        "all_declared_exact_checks_pass": all(verified.values()),
        "quantified_full_dimensional_local_implicit_branch_box_derived": True,
        "nonclaims": nonclaims,
    }
