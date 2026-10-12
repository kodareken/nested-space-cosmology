"""Exact local first-order lift of the complete REF1 spherical residual.

FO1-RC1 rewrites one local two-jet of the six RED1 base fields as a first
jet of the eighteen variables ``U=(u,p,q)``, where ``p=partial_t u`` and
``q=partial_r u``.  The complete REF1 equations remain an implicit,
generally fully nonlinear relation in ``partial_t p``.  At the IMP1 flat
FGC-QR root this module differentiates that relation exactly and publishes the
complete linearized first-order map.  It also derives the kinematic radial
reduction-constraint identity associated with the declared coordinate
derivative reduction.

This module does not implement a nonlinear acceleration solver, quantify the
implicit-function neighborhood, derive the metric-defined gauge or physical
constraint systems, prove open-domain hyperbolicity, or authorize evolution.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction
from numbers import Integral
from typing import Any, Mapping, Sequence

from .exact_linear_algebra import (
    Matrix,
    matrix_det,
    matrix_inverse,
    matrix_mul,
    matrix_rank,
)
from .exact_tangent import FirstTangent, primal_and_tangent
from .modified_harmonic import (
    modified_harmonic_symbol,
    standard_first_order_principal_system,
)
from .modified_harmonic_implicit import (
    EXPECTED_FGCQR_PARAMETERS,
    exact_full_residual_acceleration_jacobian,
)
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
FIELD_COUNT = len(BASE_FIELD_ORDER)

FO1_FIELD_ORDER = BASE_FIELD_ORDER
FO1_U_ORDER = tuple(f"u.{field}" for field in FO1_FIELD_ORDER)
FO1_P_ORDER = tuple(f"p.{field}" for field in FO1_FIELD_ORDER)
FO1_Q_ORDER = tuple(f"q.{field}" for field in FO1_FIELD_ORDER)
FO1_STATE_ORDER = FO1_U_ORDER + FO1_P_ORDER + FO1_Q_ORDER
FO1_TIME_DERIVATIVE_ORDER = tuple(f"dt.{name}" for name in FO1_STATE_ORDER)
FO1_RADIAL_DERIVATIVE_ORDER = tuple(f"dr.{name}" for name in FO1_STATE_ORDER)
FO1_U_DEFINITION_EQUATION_ORDER = tuple(
    f"dt_u_minus_p.{field}" for field in FO1_FIELD_ORDER
)
FO1_PHYSICAL_EQUATION_ORDER = MHG2_FULL_EQUATION_ORDER
FO1_Q_COMPATIBILITY_EQUATION_ORDER = tuple(
    f"dt_q_minus_dr_p.{field}" for field in FO1_FIELD_ORDER
)
FO1_EQUATION_ORDER = (
    FO1_U_DEFINITION_EQUATION_ORDER
    + FO1_PHYSICAL_EQUATION_ORDER
    + FO1_Q_COMPATIBILITY_EQUATION_ORDER
)
FO1_REDUCTION_CONSTRAINT_ORDER = tuple(
    f"C_r.{field}=q.{field}-dr_u.{field}" for field in FO1_FIELD_ORDER
)

# The complete REF1 residual depends on these six local argument families.
# q_t is used only by the kinematic compatibility row and u_t/u_r only by the
# definition/reduction rows.
FO1_PHYSICAL_ARGUMENT_GROUPS = ("u", "p", "q", "p_t", "p_r", "q_r")
FO1_PHYSICAL_ARGUMENT_ORDER = tuple(
    f"{group}.{field}"
    for group in FO1_PHYSICAL_ARGUMENT_GROUPS
    for field in FO1_FIELD_ORDER
)
FO1_IMPLICIT_NONACCELERATION_GROUPS = ("u", "p", "q", "p_r", "q_r")
FO1_IMPLICIT_NONACCELERATION_ORDER = tuple(
    f"{group}.{field}"
    for group in FO1_IMPLICIT_NONACCELERATION_GROUPS
    for field in FO1_FIELD_ORDER
)

_GROUP_TO_JET_SLOT = {
    "u": "value",
    "p": "dt",
    "q": "dr",
    "p_t": "dtt",
    "p_r": "dtr",
    "q_r": "drr",
}


def _exact_fraction(name: str, value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational value")
    return value if isinstance(value, Fraction) else Q(int(value))


def _identity(size: int) -> Matrix:
    return tuple(
        tuple(Q(int(row == column)) for column in range(size))
        for row in range(size)
    )


def _zero_matrix(rows: int, columns: int) -> Matrix:
    return tuple(tuple(Q(0) for _ in range(columns)) for _ in range(rows))


def _matrix_add(left: Matrix, right: Matrix) -> Matrix:
    if len(left) != len(right) or any(
        len(left_row) != len(right_row)
        for left_row, right_row in zip(left, right, strict=True)
    ):
        raise ValueError("matrix shapes differ")
    return tuple(
        tuple(a + b for a, b in zip(left_row, right_row, strict=True))
        for left_row, right_row in zip(left, right, strict=True)
    )


def _matrix_negate(value: Matrix) -> Matrix:
    return tuple(tuple(-entry for entry in row) for row in value)


def _matrix_hstack(*blocks: Matrix) -> Matrix:
    if not blocks or any(len(block) != len(blocks[0]) for block in blocks):
        raise ValueError("horizontal blocks must have a common row count")
    return tuple(
        tuple(entry for block in blocks for entry in block[row])
        for row in range(len(blocks[0]))
    )


def _matrix_slice_columns(value: Matrix, start: int, stop: int) -> Matrix:
    return tuple(tuple(row[start:stop]) for row in value)


def _vector_is_zero(value: Sequence[Fraction]) -> bool:
    return all(entry == 0 for entry in value)


def _state_has_tangent(state: SphericalState) -> bool:
    for field in FO1_FIELD_ORDER:
        jet = getattr(state, field)
        if any(
            isinstance(getattr(jet, slot), FirstTangent)
            for slot in ("value", "dt", "dr", "dtt", "dtr", "drr")
        ):
            return True
    return any(
        isinstance(getattr(state, name), FirstTangent)
        for name in ("planck_mass", "beta", "mu", "g4", "alpha", "eta")
    )


@dataclass(frozen=True, slots=True)
class FirstOrderFieldJet:
    """One local first jet of the first-order state triple ``(u,p,q)``.

    The nine entries are the values ``(u,p,q)`` followed by their coordinate
    derivatives needed by the declared 1+1 relation.  Keeping ``u_t`` and
    ``u_r`` independent from ``p`` and ``q`` lets the definition and reduction
    constraints be evaluated instead of silently assuming them.
    """

    u: Fraction
    p: Fraction
    q: Fraction
    u_t: Fraction
    p_t: Fraction
    q_t: Fraction
    u_r: Fraction
    p_r: Fraction
    q_r: Fraction

    def __post_init__(self) -> None:
        for name in (
            "u",
            "p",
            "q",
            "u_t",
            "p_t",
            "q_t",
            "u_r",
            "p_r",
            "q_r",
        ):
            object.__setattr__(self, name, _exact_fraction(name, getattr(self, name)))


@dataclass(frozen=True, slots=True)
class FirstOrderSphericalState:
    """One exact local first jet of ``U=(u,p,q)`` in frozen base-field order."""

    h_tt: FirstOrderFieldJet
    h_tr: FirstOrderFieldJet
    h_rr: FirstOrderFieldJet
    areal_radius: FirstOrderFieldJet
    phi: FirstOrderFieldJet
    chi: FirstOrderFieldJet
    planck_mass: Fraction = Q(1)
    beta: Fraction = Q(0)
    mu: Fraction = Q(1)
    g4: Fraction = Q(1)
    alpha: Fraction = Q(0)
    eta: Fraction = Q(0)
    branch: str = "FGC-QR"

    def __post_init__(self) -> None:
        for field in FO1_FIELD_ORDER:
            if not isinstance(getattr(self, field), FirstOrderFieldJet):
                raise TypeError(f"{field} must be a FirstOrderFieldJet")
        for name in ("planck_mass", "beta", "mu", "g4", "alpha", "eta"):
            object.__setattr__(
                self, name, _exact_fraction(name, getattr(self, name))
            )
        # Reuse the established exact domain and branch validation.  The
        # physical residual uses p_r for the symmetric mixed second derivative;
        # q_t remains an independently checkable compatibility slot.
        spherical_state_from_first_order_state(
            self, require_kinematic_equations=False
        )


def first_order_state_from_spherical_state(
    state: SphericalState,
) -> FirstOrderSphericalState:
    """Losslessly repack one ordinary-rational two-jet as a first-order jet."""

    if not isinstance(state, SphericalState):
        raise TypeError("state must be a SphericalState")
    if _state_has_tangent(state):
        raise ValueError("FO1 public state conversion does not accept tangent inputs")

    def field_jet(field: str) -> FirstOrderFieldJet:
        jet = getattr(state, field)
        return FirstOrderFieldJet(
            u=jet.value,
            p=jet.dt,
            q=jet.dr,
            u_t=jet.dt,
            p_t=jet.dtt,
            q_t=jet.dtr,
            u_r=jet.dr,
            p_r=jet.dtr,
            q_r=jet.drr,
        )

    return FirstOrderSphericalState(
        **{field: field_jet(field) for field in FO1_FIELD_ORDER},
        planck_mass=state.planck_mass,
        beta=state.beta,
        mu=state.mu,
        g4=state.g4,
        alpha=state.alpha,
        eta=state.eta,
        branch=state.branch,
    )


def spherical_state_from_first_order_state(
    state: FirstOrderSphericalState,
    *,
    require_kinematic_equations: bool = True,
) -> SphericalState:
    """Reconstruct the REF1 two-jet used by the implicit physical rows.

    A genuine lift additionally satisfies ``u_t=p``, ``u_r=q``, and
    ``q_t=p_r``.  The last equality is mixed-partial compatibility.  With
    ``require_kinematic_equations=False`` the physical residual can be
    evaluated off that kinematic shell while the violations remain visible in
    :func:`first_order_kinematic_residuals`.
    """

    if not isinstance(state, FirstOrderSphericalState):
        raise TypeError("state must be a FirstOrderSphericalState")
    kinematic = first_order_kinematic_residuals(state)
    if require_kinematic_equations and not all(
        _vector_is_zero(kinematic[name])
        for name in (
            "u_time_definition_residual",
            "radial_reduction_constraint",
            "mixed_partial_compatibility_residual",
        )
    ):
        raise ValueError("first-order jet violates the declared kinematic equations")

    return SphericalState(
        **{
            field: Jet2(
                value=getattr(state, field).u,
                dt=getattr(state, field).p,
                dr=getattr(state, field).q,
                dtt=getattr(state, field).p_t,
                dtr=getattr(state, field).p_r,
                drr=getattr(state, field).q_r,
            )
            for field in FO1_FIELD_ORDER
        },
        planck_mass=state.planck_mass,
        beta=state.beta,
        mu=state.mu,
        g4=state.g4,
        alpha=state.alpha,
        eta=state.eta,
        branch=state.branch,
    )


def first_order_kinematic_residuals(
    state: FirstOrderSphericalState,
) -> dict[str, tuple[Fraction, ...]]:
    """Evaluate the definition, reduction, and mixed-partial rows exactly."""

    if not isinstance(state, FirstOrderSphericalState):
        raise TypeError("state must be a FirstOrderSphericalState")
    u_time = tuple(
        getattr(state, field).u_t - getattr(state, field).p
        for field in FO1_FIELD_ORDER
    )
    radial_constraint = tuple(
        getattr(state, field).q - getattr(state, field).u_r
        for field in FO1_FIELD_ORDER
    )
    compatibility = tuple(
        getattr(state, field).q_t - getattr(state, field).p_r
        for field in FO1_FIELD_ORDER
    )
    return {
        "field_order": FO1_FIELD_ORDER,
        "u_time_definition_residual": u_time,
        "radial_reduction_constraint": radial_constraint,
        "mixed_partial_compatibility_residual": compatibility,
        # If u_t=p holds as a differentiable field equation, then
        # d_t(q-d_r u)=q_t-d_r p.  This is the on-definition-shell value.
        "radial_reduction_constraint_time_derivative_on_definition_shell": compatibility,
    }


def modified_harmonic_first_order_dae_residuals(
    state: FirstOrderSphericalState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: int | Fraction,
    tilde_normal_factor: int | Fraction,
    hat_normal_factor: int | Fraction,
) -> dict[str, Any]:
    """Evaluate the exact 18-row local differential-algebraic relation."""

    spherical = spherical_state_from_first_order_state(
        state, require_kinematic_equations=False
    )
    physical = modified_harmonic_full_residuals(
        spherical,
        reference=reference,
        coordinate_radius=coordinate_radius,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    kinematic = first_order_kinematic_residuals(state)
    full_vector = (
        kinematic["u_time_definition_residual"]
        + tuple(physical["full_residual_vector"])
        + kinematic["mixed_partial_compatibility_residual"]
    )
    return {
        "state_order": FO1_STATE_ORDER,
        "equation_order": FO1_EQUATION_ORDER,
        "physical_equation_order": FO1_PHYSICAL_EQUATION_ORDER,
        "physical_residual_vector": tuple(physical["full_residual_vector"]),
        "kinematic": kinematic,
        "dae_residual_vector": full_vector,
        "reconstructed_spherical_state": spherical,
        "physical_residual_is_complete_ref1": True,
        "physical_relation_is_implicit_in_p_t": True,
        "nonlinear_p_t_solve_performed": False,
    }


def exact_full_residual_first_order_argument_jacobian(
    state: SphericalState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: int | Fraction,
    tilde_normal_factor: int | Fraction,
    hat_normal_factor: int | Fraction,
) -> dict[str, Any]:
    """Differentiate REF1 in all 36 physical first-order argument slots."""

    if not isinstance(state, SphericalState):
        raise TypeError("state must be a SphericalState")
    if _state_has_tangent(state):
        raise ValueError("FO1 Jacobian input must not already contain tangents")
    base = tuple(
        modified_harmonic_full_residuals(
            state,
            reference=reference,
            coordinate_radius=coordinate_radius,
            tilde_normal_factor=tilde_normal_factor,
            hat_normal_factor=hat_normal_factor,
        )["full_residual_vector"]
    )
    columns: list[tuple[Fraction, ...]] = []
    primals: list[tuple[Fraction, ...]] = []
    for group in FO1_PHYSICAL_ARGUMENT_GROUPS:
        slot = _GROUP_TO_JET_SLOT[group]
        for field in FO1_FIELD_ORDER:
            jet = getattr(state, field)
            seeded_state = replace(
                state,
                **{
                    field: replace(
                        jet,
                        **{slot: FirstTangent.seed(getattr(jet, slot))},
                    )
                },
            )
            residual = modified_harmonic_full_residuals(
                seeded_state,
                reference=reference,
                coordinate_radius=coordinate_radius,
                tilde_normal_factor=tilde_normal_factor,
                hat_normal_factor=hat_normal_factor,
            )["full_residual_vector"]
            parts = tuple(primal_and_tangent(value) for value in residual)
            primals.append(tuple(value[0] for value in parts))
            columns.append(tuple(value[1] for value in parts))

    jacobian = tuple(
        tuple(columns[column][row] for column in range(len(columns)))
        for row in range(len(FO1_PHYSICAL_EQUATION_ORDER))
    )
    blocks = {
        group: _matrix_slice_columns(
            jacobian,
            index * FIELD_COUNT,
            (index + 1) * FIELD_COUNT,
        )
        for index, group in enumerate(FO1_PHYSICAL_ARGUMENT_GROUPS)
    }

    # A simultaneous weighted tangent is an independent assembly check for the
    # stored single-column matrix (the evaluator path is shared, the seeding
    # and column assembly are not).
    weighted = state
    weights = tuple(Q(index + 1) for index in range(len(columns)))
    for group_index, group in enumerate(FO1_PHYSICAL_ARGUMENT_GROUPS):
        slot = _GROUP_TO_JET_SLOT[group]
        for field_index, field in enumerate(FO1_FIELD_ORDER):
            jet = getattr(weighted, field)
            weight = weights[group_index * FIELD_COUNT + field_index]
            current = getattr(jet, slot)
            primal = current.primal if isinstance(current, FirstTangent) else current
            tangent = current.tangent if isinstance(current, FirstTangent) else Q(0)
            weighted = replace(
                weighted,
                **{
                    field: replace(
                        jet,
                        **{slot: FirstTangent(primal, tangent + weight)},
                    )
                },
            )
    weighted_residual = modified_harmonic_full_residuals(
        weighted,
        reference=reference,
        coordinate_radius=coordinate_radius,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )["full_residual_vector"]
    weighted_parts = tuple(primal_and_tangent(value) for value in weighted_residual)
    predicted_weighted_tangent = tuple(
        sum(jacobian[row][column] * weights[column] for column in range(len(weights)))
        for row in range(len(jacobian))
    )
    return {
        "equation_order": FO1_PHYSICAL_EQUATION_ORDER,
        "argument_order": FO1_PHYSICAL_ARGUMENT_ORDER,
        "base_full_residual_vector": base,
        "seeded_primal_residual_vectors": tuple(primals),
        "seeded_primals_reproduce_base_residual": all(
            primal == base for primal in primals
        ),
        "jacobian": jacobian,
        "blocks": blocks,
        "weighted_direction": weights,
        "weighted_direction_primal_reproduces_base": tuple(
            value[0] for value in weighted_parts
        )
        == base,
        "weighted_direction_tangent": tuple(value[1] for value in weighted_parts),
        "weighted_direction_predicted_tangent": predicted_weighted_tangent,
        "weighted_direction_superposition_exact": tuple(
            value[1] for value in weighted_parts
        )
        == predicted_weighted_tangent,
        "differentiation_method": "exact_first_tangent_forward_ad",
    }


def flat_root_linearized_first_order_system(
    state: SphericalState,
    argument_jacobian: Mapping[str, Any],
    *,
    tilde_normal_factor: int | Fraction,
    hat_normal_factor: int | Fraction,
) -> dict[str, Any]:
    """Solve only the exact derivative of the IMP1 branch at its flat root."""

    blocks = argument_jacobian.get("blocks")
    if not isinstance(blocks, Mapping) or set(blocks) != set(
        FO1_PHYSICAL_ARGUMENT_GROUPS
    ):
        raise ValueError("FO1 argument Jacobian blocks are incomplete")
    acceleration = blocks["p_t"]
    determinant = matrix_det(acceleration)
    if determinant == 0:
        raise ValueError("FO1 flat-root acceleration block is singular")
    inverse = matrix_inverse(acceleration)
    nonacceleration = _matrix_hstack(
        *(blocks[group] for group in FO1_IMPLICIT_NONACCELERATION_GROUPS)
    )
    branch_derivative = _matrix_negate(matrix_mul(inverse, nonacceleration))
    derivative_identity = _matrix_add(
        matrix_mul(acceleration, branch_derivative), nonacceleration
    ) == _zero_matrix(FIELD_COUNT, len(FO1_IMPLICIT_NONACCELERATION_ORDER))

    lower_state_derivative = _matrix_slice_columns(
        branch_derivative, 0, 3 * FIELD_COUNT
    )
    radial_derivative = _matrix_slice_columns(
        branch_derivative, 3 * FIELD_COUNT, 5 * FIELD_COUNT
    )
    p_principal = _matrix_negate(radial_derivative)
    zero6 = _zero_matrix(FIELD_COUNT, FIELD_COUNT)
    identity6 = _identity(FIELD_COUNT)
    minus_identity6 = _matrix_negate(identity6)
    principal_12 = tuple(
        tuple(p_principal[row]) for row in range(FIELD_COUNT)
    ) + tuple(
        tuple(minus_identity6[row] + zero6[row]) for row in range(FIELD_COUNT)
    )
    principal_18 = (
        tuple(tuple(Q(0) for _ in range(3 * FIELD_COUNT)) for _ in range(FIELD_COUNT))
        + tuple(
            tuple(
                (Q(0),) * FIELD_COUNT
                + p_principal[row]
            )
            for row in range(FIELD_COUNT)
        )
        + tuple(
            tuple(
                (Q(0),) * FIELD_COUNT
                + minus_identity6[row]
                + zero6[row]
            )
            for row in range(FIELD_COUNT)
        )
    )

    source_18 = (
        tuple(
            tuple(
                Q(int(column == FIELD_COUNT + row))
                for column in range(3 * FIELD_COUNT)
            )
            for row in range(FIELD_COUNT)
        )
        + tuple(tuple(lower_state_derivative[row]) for row in range(FIELD_COUNT))
        + tuple(
            tuple(Q(0) for _ in range(3 * FIELD_COUNT))
            for _ in range(FIELD_COUNT)
        )
    )

    mhg1 = standard_first_order_principal_system(
        modified_harmonic_symbol(
            state,
            tilde_normal_factor=tilde_normal_factor,
            hat_normal_factor=hat_normal_factor,
        )
    )
    return {
        "linearization_point": "IMP1_flat_FGCQR_reference_gauge_root",
        "state_order": FO1_STATE_ORDER,
        "nonacceleration_argument_order": FO1_IMPLICIT_NONACCELERATION_ORDER,
        "acceleration_jacobian": acceleration,
        "acceleration_jacobian_determinant": determinant,
        "acceleration_jacobian_rank": matrix_rank(acceleration),
        "acceleration_jacobian_inverse": inverse,
        "implicit_branch_derivative": branch_derivative,
        "implicit_branch_derivative_identity_exact": derivative_identity,
        "linearized_radial_principal_matrix": principal_18,
        "linearized_lower_order_source_matrix": source_18,
        "p_q_radial_principal_matrix": principal_12,
        "mhg1_standard_first_order_principal_matrix": mhg1[
            "normal_principal_matrix"
        ],
        "p_q_principal_matrix_matches_mhg1_exact": principal_12
        == mhg1["normal_principal_matrix"],
        "classification": "complete_exact_flat_root_linearization_of_local_implicit_branch_not_nonlinear_source_evaluator",
    }


def required_fo1_rc1_nonclaims() -> dict[str, bool]:
    return {
        name: False
        for name in (
            "explicit_nonlinear_acceleration_map_implemented",
            "solved_acceleration_map_away_from_flat_root_derived",
            "quantified_implicit_branch_neighborhood_proven",
            "complete_nonlinear_lower_order_first_order_sources_derived",
            "complete_quasilinear_first_order_evolution_system_derived",
            "metric_derived_gauge_constraint_propagation_proven",
            "complete_reduction_constraint_system_propagation_proven",
            "physical_hamiltonian_momentum_constraints_derived",
            "physical_initial_constraints_solved",
            "constraint_preserving_initial_boundary_value_problem_proven",
            "uniform_open_domain_symmetrizer_proven",
            "open_retained_eft_domain_proven",
            "boundary_or_regular_center_system_derived",
            "evolution_authorized",
            "collapse_solution_derived",
            "metric_null_affine_defocusing_derived",
            "singularity_resolution_derived",
        )
    }


def modified_harmonic_first_order_certificate(
    configuration: Mapping[str, Any],
) -> dict[str, Any]:
    """Build the exact FO1-RC1 two-control first-order-lift certificate."""

    expected_keys = {
        "flat_fixture",
        "activated_fixture",
        "reference",
        "flat_coordinate_radius",
        "activated_coordinate_radius",
        "tilde_normal_factor",
        "hat_normal_factor",
    }
    if not isinstance(configuration, Mapping) or set(configuration) != expected_keys:
        raise ValueError("FO1-RC1 configuration adapter has unexpected keys")
    reference = configuration["reference"]
    if not isinstance(reference, ReferenceConnection):
        raise ValueError("FO1-RC1 reference must be a ReferenceConnection")
    flat_radius = Q(configuration["flat_coordinate_radius"])
    activated_radius = Q(configuration["activated_coordinate_radius"])
    tilde_factor = Q(configuration["tilde_normal_factor"])
    hat_factor = Q(configuration["hat_normal_factor"])
    if not 1 < tilde_factor < hat_factor:
        raise ValueError("FO1-RC1 requires 1 < tilde factor < hat factor")
    if min(flat_radius, activated_radius) <= reference.radial_domain_minimum:
        raise ValueError("FO1-RC1 controls must lie inside the reference annulus")
    flat_fixture = configuration["flat_fixture"]
    activated_fixture = configuration["activated_fixture"]
    if not isinstance(flat_fixture, Mapping) or not isinstance(
        activated_fixture, Mapping
    ):
        raise ValueError("FO1-RC1 fixtures must be mappings")

    flat_state = state_from_generalized_adm_pg_fixture(flat_fixture)
    activated_state = state_from_generalized_adm_pg_fixture(activated_fixture)
    flat_fo = first_order_state_from_spherical_state(flat_state)
    activated_fo = first_order_state_from_spherical_state(activated_state)
    flat_roundtrip = spherical_state_from_first_order_state(flat_fo)
    activated_roundtrip = spherical_state_from_first_order_state(activated_fo)

    def direct(state: SphericalState, radius: Fraction) -> tuple[Fraction, ...]:
        return tuple(
            modified_harmonic_full_residuals(
                state,
                reference=reference,
                coordinate_radius=radius,
                tilde_normal_factor=tilde_factor,
                hat_normal_factor=hat_factor,
            )["full_residual_vector"]
        )

    flat_direct = direct(flat_state, flat_radius)
    activated_direct = direct(activated_state, activated_radius)
    flat_dae = modified_harmonic_first_order_dae_residuals(
        flat_fo,
        reference=reference,
        coordinate_radius=flat_radius,
        tilde_normal_factor=tilde_factor,
        hat_normal_factor=hat_factor,
    )
    activated_dae = modified_harmonic_first_order_dae_residuals(
        activated_fo,
        reference=reference,
        coordinate_radius=activated_radius,
        tilde_normal_factor=tilde_factor,
        hat_normal_factor=hat_factor,
    )
    flat_jacobian = exact_full_residual_first_order_argument_jacobian(
        flat_state,
        reference=reference,
        coordinate_radius=flat_radius,
        tilde_normal_factor=tilde_factor,
        hat_normal_factor=hat_factor,
    )
    activated_jacobian = exact_full_residual_first_order_argument_jacobian(
        activated_state,
        reference=reference,
        coordinate_radius=activated_radius,
        tilde_normal_factor=tilde_factor,
        hat_normal_factor=hat_factor,
    )
    imp1 = exact_full_residual_acceleration_jacobian(
        flat_state,
        reference=reference,
        coordinate_radius=flat_radius,
        tilde_normal_factor=tilde_factor,
        hat_normal_factor=hat_factor,
    )
    linearized = flat_root_linearized_first_order_system(
        flat_state,
        flat_jacobian,
        tilde_normal_factor=tilde_factor,
        hat_normal_factor=hat_factor,
    )
    flat_kinematic = flat_dae["kinematic"]
    activated_kinematic = activated_dae["kinematic"]
    flat_parameters = {
        "planck_mass": flat_state.planck_mass,
        "beta": flat_state.beta,
        "mu": flat_state.mu,
        "g4": flat_state.g4,
        "alpha": flat_state.alpha,
        "eta": flat_state.eta,
    }
    verified = {
        "frozen_eighteen_variable_and_equation_orders_exact": len(FO1_STATE_ORDER)
        == 18
        and len(FO1_EQUATION_ORDER) == 18,
        "flat_fgcqr_action_parameters_match_imp1_exact": flat_parameters
        == EXPECTED_FGCQR_PARAMETERS,
        "flat_first_order_two_jet_roundtrip_exact": flat_roundtrip == flat_state,
        "activated_first_order_two_jet_roundtrip_exact": activated_roundtrip
        == activated_state,
        "flat_complete_ref1_residual_reproduced_exact": flat_dae[
            "physical_residual_vector"
        ]
        == flat_direct,
        "activated_complete_ref1_residual_reproduced_exact": activated_dae[
            "physical_residual_vector"
        ]
        == activated_direct,
        "flat_all_kinematic_rows_and_reduction_constraints_zero_exact": all(
            _vector_is_zero(flat_kinematic[name])
            for name in (
                "u_time_definition_residual",
                "radial_reduction_constraint",
                "mixed_partial_compatibility_residual",
            )
        ),
        "activated_all_kinematic_rows_and_reduction_constraints_zero_exact": all(
            _vector_is_zero(activated_kinematic[name])
            for name in (
                "u_time_definition_residual",
                "radial_reduction_constraint",
                "mixed_partial_compatibility_residual",
            )
        ),
        "kinematic_radial_reduction_constraint_time_derivative_zero_exact": _vector_is_zero(
            flat_kinematic[
                "radial_reduction_constraint_time_derivative_on_definition_shell"
            ]
        )
        and _vector_is_zero(
            activated_kinematic[
                "radial_reduction_constraint_time_derivative_on_definition_shell"
            ]
        ),
        "flat_complete_ref1_residual_zero_exact": _vector_is_zero(flat_direct),
        "flat_all_argument_tangent_primals_reproduce_residual": flat_jacobian[
            "seeded_primals_reproduce_base_residual"
        ],
        "activated_all_argument_tangent_primals_reproduce_residual": activated_jacobian[
            "seeded_primals_reproduce_base_residual"
        ],
        "flat_weighted_tangent_superposition_exact": flat_jacobian[
            "weighted_direction_superposition_exact"
        ],
        "activated_weighted_tangent_superposition_exact": activated_jacobian[
            "weighted_direction_superposition_exact"
        ],
        "flat_acceleration_block_matches_imp1_exact": flat_jacobian["blocks"][
            "p_t"
        ]
        == imp1["jacobian"],
        "flat_acceleration_block_has_imp1_determinant_and_rank": linearized[
            "acceleration_jacobian_determinant"
        ]
        == Q(-165888)
        and linearized["acceleration_jacobian_rank"] == FIELD_COUNT,
        "flat_implicit_branch_derivative_identity_exact": linearized[
            "implicit_branch_derivative_identity_exact"
        ],
        "flat_linearized_p_q_principal_matrix_matches_mhg1_exact": linearized[
            "p_q_principal_matrix_matches_mhg1_exact"
        ],
        "activated_control_is_off_shell_and_nontrivial": not _vector_is_zero(
            activated_direct
        ),
        "no_nonlinear_acceleration_solve_performed": not flat_dae[
            "nonlinear_p_t_solve_performed"
        ]
        and not activated_dae["nonlinear_p_t_solve_performed"],
    }
    return {
        "artifact_id": "FGC-1-HYP1-FO1-RC1",
        "classification": "exact_local_first_order_differential_algebraic_lift_with_flat_root_complete_linearization_and_kinematic_radial_reduction_constraint_identity_not_nonlinear_evolution",
        "formulation": {
            "physical_equations": "unredefined_ACT1_VAR1",
            "gauge_equations": "REF1_full_reference_connection_modified_harmonic",
            "first_order_kind": "exact_local_first_order_differential_algebraic_lift",
            "state_definition": "U=(u,p=partial_t_u,q=partial_r_u)",
            "physical_relation": "R_REF1(u,p,q,partial_t_p,partial_r_p,partial_r_q;r)=0",
            "kinematic_rows": (
                "partial_t_u-p=0",
                "partial_t_q-partial_r_p=0",
            ),
            "radial_reduction_constraint": "C_r=q-partial_r_u",
            "radial_reduction_constraint_identity": "partial_t_C_r=partial_t_q-partial_r_p=0_on_declared_kinematic_rows",
            "state_order": FO1_STATE_ORDER,
            "equation_order": FO1_EQUATION_ORDER,
            "physical_argument_order": FO1_PHYSICAL_ARGUMENT_ORDER,
            "reference_connection": reference.reference_id,
            "reference_radial_domain_minimum": reference.radial_domain_minimum,
            "tilde_normal_factor": tilde_factor,
            "hat_normal_factor": hat_factor,
        },
        "controls": {
            "flat": {
                "fixture_id": str(flat_fixture.get("fixture_id", "unspecified")),
                "model_id": str(flat_fixture.get("model_id", "unspecified")),
                "coordinate_radius": flat_radius,
                "complete_ref1_residual": flat_direct,
                "kinematic": flat_kinematic,
                "argument_jacobian": flat_jacobian,
            },
            "activated_off_shell": {
                "fixture_id": str(
                    activated_fixture.get("fixture_id", "unspecified")
                ),
                "model_id": str(
                    activated_fixture.get("model_id", "unspecified")
                ),
                "coordinate_radius": activated_radius,
                "solution_status": "off_shell_control_not_a_solution",
                "complete_ref1_residual": activated_direct,
                "kinematic": activated_kinematic,
                "argument_jacobian": activated_jacobian,
                "acceleration_block_determinant": matrix_det(
                    activated_jacobian["blocks"]["p_t"]
                ),
                "acceleration_block_rank": matrix_rank(
                    activated_jacobian["blocks"]["p_t"]
                ),
            },
        },
        "flat_root_linearized_first_order_system": linearized,
        "verified_exact_checks": verified,
        "all_declared_exact_checks_pass": all(verified.values()),
        "exact_local_first_order_dae_lift_derived": all(verified.values()),
        "complete_flat_root_linearized_first_order_map_derived": all(
            verified[name]
            for name in (
                "flat_acceleration_block_matches_imp1_exact",
                "flat_acceleration_block_has_imp1_determinant_and_rank",
                "flat_implicit_branch_derivative_identity_exact",
                "flat_linearized_p_q_principal_matrix_matches_mhg1_exact",
            )
        ),
        "kinematic_radial_reduction_constraint_identity_derived": verified[
            "kinematic_radial_reduction_constraint_time_derivative_zero_exact"
        ],
        "nonclaims": required_fo1_rc1_nonclaims(),
    }
