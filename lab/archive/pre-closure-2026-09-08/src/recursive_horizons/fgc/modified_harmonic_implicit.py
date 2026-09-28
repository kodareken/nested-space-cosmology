"""Exact local implicit-acceleration gate for the full REF1 residual.

IMP1 differentiates the complete six-equation reference-gauge residual with
respect to the six base-field coordinate-time accelerations using exact
first-tangent arithmetic.  At one declared flat FGC-QR vacuum root it checks
the hypotheses of the finite-dimensional implicit function theorem.  The
result is a local algebraic branch, not a first-order PDE, propagation theorem,
quantified state domain, or evolution authorization.
"""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
from typing import Any, Mapping, Sequence

from .exact_linear_algebra import (
    Matrix,
    matrix_det,
    matrix_inverse,
    matrix_mul,
    matrix_rank,
)
from .exact_tangent import FirstTangent, primal_and_tangent
from .modified_harmonic import modified_harmonic_symbol
from .modified_harmonic_reference import (
    MHG2_FULL_EQUATION_ORDER,
    modified_harmonic_full_residuals,
)
from .reference_connection import ReferenceConnection
from .spherical_reduction import (
    BASE_FIELD_ORDER,
    Jet2,
    SphericalState,
    state_from_generalized_adm_pg_fixture,
)


Q = Fraction

IMPLICIT_ACCELERATION_ORDER = tuple(
    f"{field}.dtt" for field in BASE_FIELD_ORDER
)
IMPLICIT_EQUATION_ORDER = MHG2_FULL_EQUATION_ORDER
EXPECTED_FGCQR_PARAMETERS = {
    "planck_mass": Q(2),
    "beta": Q(-1, 4),
    "mu": Q(3),
    "g4": Q(1, 2),
    "alpha": Q(0),
    "eta": Q(1, 2),
}


def _tensor_zero(value: Any) -> bool:
    if isinstance(value, Fraction):
        return value == 0
    if isinstance(value, FirstTangent):
        return value.primal == 0 and value.tangent == 0
    if isinstance(value, (tuple, list)):
        return all(_tensor_zero(item) for item in value)
    raise ValueError("exact tensor contains an unsupported scalar")


def _identity(size: int) -> Matrix:
    return tuple(
        tuple(Q(int(row == column)) for column in range(size))
        for row in range(size)
    )


def _jet_is(value: Jet2, *, field_value: Fraction, dr: Fraction = Q(0)) -> bool:
    return value == Jet2(field_value, dr=dr)


def _flat_fgcqr_root_checks(
    state: SphericalState,
    *,
    coordinate_radius: Fraction,
) -> dict[str, bool]:
    parameter_values = {
        "planck_mass": state.planck_mass,
        "beta": state.beta,
        "mu": state.mu,
        "g4": state.g4,
        "alpha": state.alpha,
        "eta": state.eta,
    }
    return {
        "fgcqr_branch_retained": state.branch == "FGC-QR",
        "fgcqr_action_parameters_exact": parameter_values
        == EXPECTED_FGCQR_PARAMETERS,
        "flat_base_metric_two_jet_exact": _jet_is(state.h_tt, field_value=Q(-1))
        and _jet_is(state.h_tr, field_value=Q(0))
        and _jet_is(state.h_rr, field_value=Q(1)),
        "areal_radius_matches_reference_chart_exact": _jet_is(
            state.areal_radius,
            field_value=coordinate_radius,
            dr=Q(1),
        ),
        "vacuum_scalar_two_jets_exact": _jet_is(state.phi, field_value=Q(0))
        and _jet_is(state.chi, field_value=Q(0)),
        "base_coordinate_time_accelerations_zero": all(
            getattr(state, field).dtt == 0 for field in BASE_FIELD_ORDER
        ),
    }


def exact_full_residual_acceleration_jacobian(
    state: SphericalState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: int | Fraction,
    tilde_normal_factor: int | Fraction,
    hat_normal_factor: int | Fraction,
) -> dict[str, Any]:
    """Differentiate the complete REF1 residual by exact forward tangents.

    Each column seeds one base-field ``dtt`` entry with derivative one.  The
    primal computation is required to reproduce the unseeded residual in every
    column, so this path cannot silently substitute the MHG1 principal matrix.
    """

    radius = Q(coordinate_radius)
    base = modified_harmonic_full_residuals(
        state,
        reference=reference,
        coordinate_radius=radius,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    base_vector = tuple(base["full_residual_vector"])
    columns: list[tuple[Fraction, ...]] = []
    primal_columns: list[tuple[Fraction, ...]] = []
    for field in BASE_FIELD_ORDER:
        jet = getattr(state, field)
        if isinstance(jet.dtt, FirstTangent):
            raise ValueError("IMP1 input state must not already contain tangents")
        seeded_state = replace(
            state,
            **{
                field: replace(
                    jet,
                    dtt=FirstTangent.seed(jet.dtt),
                )
            },
        )
        seeded = modified_harmonic_full_residuals(
            seeded_state,
            reference=reference,
            coordinate_radius=radius,
            tilde_normal_factor=tilde_normal_factor,
            hat_normal_factor=hat_normal_factor,
        )["full_residual_vector"]
        parts = tuple(primal_and_tangent(value) for value in seeded)
        primal_columns.append(tuple(value[0] for value in parts))
        columns.append(tuple(value[1] for value in parts))

    jacobian = tuple(
        tuple(columns[column][row] for column in range(len(BASE_FIELD_ORDER)))
        for row in range(len(IMPLICIT_EQUATION_ORDER))
    )
    return {
        "equation_order": IMPLICIT_EQUATION_ORDER,
        "acceleration_order": IMPLICIT_ACCELERATION_ORDER,
        "base_acceleration_vector": tuple(
            getattr(state, field).dtt for field in BASE_FIELD_ORDER
        ),
        "base_full_residual_vector": base_vector,
        "seeded_primal_residual_vectors": tuple(primal_columns),
        "seeded_primals_reproduce_base_residual": all(
            column == base_vector for column in primal_columns
        ),
        "jacobian": jacobian,
        "differentiation_method": "exact_first_tangent_forward_ad",
    }


def _coordinate_time_kinetic_block(
    state: SphericalState,
    *,
    tilde_normal_factor: int | Fraction,
    hat_normal_factor: int | Fraction,
) -> Matrix:
    symbol = modified_harmonic_symbol(
        state,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    return tuple(
        tuple(entry[2] if len(entry) > 2 else Q(0) for entry in row)
        for row in symbol
    )


def required_imp1_nonclaims() -> dict[str, bool]:
    return {
        name: False
        for name in (
            "quantified_implicit_branch_neighborhood_proven",
            "nonzero_activated_solution_on_branch_derived",
            "complete_lower_order_first_order_sources_derived",
            "gauge_constraint_lower_order_propagation_proven",
            "reduction_constraint_propagation_proven",
            "physical_initial_constraints_solved",
            "uniform_open_domain_symmetrizer_proven",
            "open_retained_eft_domain_proven",
            "boundary_or_regular_center_system_derived",
            "evolution_authorized",
            "collapse_solution_derived",
            "metric_null_affine_defocusing_derived",
            "singularity_resolution_derived",
        )
    }


def modified_harmonic_implicit_certificate(
    configuration: Mapping[str, Any],
) -> dict[str, Any]:
    """Build the exact IMP1 local implicit-function certificate."""

    expected_keys = {
        "fixture",
        "reference",
        "coordinate_radius",
        "tilde_normal_factor",
        "hat_normal_factor",
    }
    if not isinstance(configuration, Mapping) or set(configuration) != expected_keys:
        raise ValueError("IMP1 configuration adapter has unexpected keys")
    fixture = configuration["fixture"]
    reference = configuration["reference"]
    if not isinstance(fixture, Mapping):
        raise ValueError("IMP1 fixture must be a mapping")
    if not isinstance(reference, ReferenceConnection):
        raise ValueError("IMP1 reference must be a ReferenceConnection")
    radius = Q(configuration["coordinate_radius"])
    tilde_factor = Q(configuration["tilde_normal_factor"])
    hat_factor = Q(configuration["hat_normal_factor"])
    if not 1 < tilde_factor < hat_factor:
        raise ValueError("IMP1 requires 1 < tilde factor < hat factor")
    if radius <= reference.radial_domain_minimum:
        raise ValueError("IMP1 root must lie strictly inside the reference annulus")

    state = state_from_generalized_adm_pg_fixture(fixture)
    root_checks = _flat_fgcqr_root_checks(state, coordinate_radius=radius)
    if not all(root_checks.values()):
        failed = sorted(name for name, passed in root_checks.items() if not passed)
        raise ValueError(f"IMP1 declared flat FGC-QR root failed: {failed}")

    full = modified_harmonic_full_residuals(
        state,
        reference=reference,
        coordinate_radius=radius,
        tilde_normal_factor=tilde_factor,
        hat_normal_factor=hat_factor,
    )
    tangent = exact_full_residual_acceleration_jacobian(
        state,
        reference=reference,
        coordinate_radius=radius,
        tilde_normal_factor=tilde_factor,
        hat_normal_factor=hat_factor,
    )
    jacobian = tangent["jacobian"]
    mhg1_kinetic = _coordinate_time_kinetic_block(
        state,
        tilde_normal_factor=tilde_factor,
        hat_normal_factor=hat_factor,
    )
    determinant = matrix_det(jacobian)
    rank = matrix_rank(jacobian)
    inverse = matrix_inverse(jacobian) if determinant != 0 else None
    identity = _identity(len(BASE_FIELD_ORDER))
    left_identity = inverse is not None and matrix_mul(inverse, jacobian) == identity
    right_identity = inverse is not None and matrix_mul(jacobian, inverse) == identity

    physical = full["gauge"].physical
    base_determinant = (
        state.h_tt.value * state.h_rr.value - state.h_tr.value**2
    )
    g00_inverse = physical.inverse_metric[0][0]
    effective_planck = full["extension"]["effective_planck_coefficient"]
    tilde_determinant = matrix_det(full["gauge"].tilde_inverse_metric)
    hat_determinant = matrix_det(full["extension"]["hat_inverse_metric"])
    smooth_domain = {
        "base_metric_determinant_negative": base_determinant < 0,
        "coordinate_time_covector_physical_timelike": g00_inverse < 0,
        "areal_radius_positive": state.areal_radius.value > 0,
        "effective_planck_coefficient_positive": effective_planck > 0,
        "reference_radius_strictly_inside_annulus": radius
        > reference.radial_domain_minimum,
        "auxiliary_inverse_metrics_lorentzian": tilde_determinant < 0
        and hat_determinant < 0,
    }
    root_residual_zero = _tensor_zero(full["full_residual_vector"])
    gauge_surface_exact = _tensor_zero(full["gauge"].constraint_up) and _tensor_zero(
        full["gauge"].covariant_constraint_derivative
    )
    extension_zero = _tensor_zero(full["extension_residual_vector"])
    verified = {
        **root_checks,
        "complete_ref1_residual_zero_at_root": root_residual_zero,
        "reference_gauge_surface_exact_at_root": gauge_surface_exact,
        "reference_gauge_extension_zero_at_root": extension_zero,
        "all_rational_smooth_domain_denominators_regular_at_root": all(
            smooth_domain.values()
        ),
        "exact_forward_tangent_primals_reproduce_root": tangent[
            "seeded_primals_reproduce_base_residual"
        ],
        "full_residual_acceleration_jacobian_matches_mhg1_kinetic_block": jacobian
        == mhg1_kinetic,
        "acceleration_jacobian_full_rank": rank == len(BASE_FIELD_ORDER),
        "acceleration_jacobian_determinant_nonzero": determinant != 0,
        "exact_two_sided_inverse_identities": left_identity and right_identity,
    }
    local_branch = (
        all(verified.values())
        and root_residual_zero
        and determinant != 0
        and all(smooth_domain.values())
    )
    return {
        "artifact_id": "FGC-1-HYP1-MHG3-IMP1",
        "classification": "exact_local_implicit_coordinate_time_acceleration_branch_at_flat_fgcqr_reference_gauge_root_not_open_domain_or_evolution",
        "formulation": {
            "physical_equations": "unredefined_ACT1_VAR1",
            "gauge_equations": "REF1_full_reference_connection_modified_harmonic",
            "acceleration_variables": "base_metric_and_scalar_coordinate_time_second_derivatives",
            "differentiation_method": tangent["differentiation_method"],
            "implicit_theorem": "finite_dimensional_real_implicit_function_theorem",
            "reference_connection": reference.reference_id,
            "reference_radial_domain_minimum": reference.radial_domain_minimum,
            "coordinate_radius": radius,
            "tilde_normal_factor": tilde_factor,
            "hat_normal_factor": hat_factor,
        },
        "reference_solution": {
            "fixture_id": str(fixture.get("fixture_id", "unspecified")),
            "model_id": str(fixture.get("model_id", "unspecified")),
            "action_parameters": EXPECTED_FGCQR_PARAMETERS,
            "base_metric_determinant": base_determinant,
            "inverse_metric_g00": g00_inverse,
            "effective_planck_coefficient": effective_planck,
            "constraint_up": full["gauge"].constraint_up,
            "covariant_constraint_derivative": full[
                "gauge"
            ].covariant_constraint_derivative,
            "gauge_extension_residual_vector": full[
                "extension_residual_vector"
            ],
            "full_residual_vector": full["full_residual_vector"],
        },
        "smoothness_domain_at_root": smooth_domain,
        "acceleration_map": {
            **tangent,
            "mhg1_coordinate_time_kinetic_block": mhg1_kinetic,
            "jacobian_determinant": determinant,
            "jacobian_rank": rank,
            "jacobian_inverse": inverse,
            "left_inverse_identity_exact": left_identity,
            "right_inverse_identity_exact": right_identity,
        },
        "verified_exact_checks": verified,
        "all_declared_exact_checks_pass": all(verified.values()),
        "local_smooth_implicit_coordinate_time_acceleration_branch_proven": local_branch,
        "implicit_branch_neighborhood": "exists_but_not_quantified",
        "nonclaims": required_imp1_nonclaims(),
    }
