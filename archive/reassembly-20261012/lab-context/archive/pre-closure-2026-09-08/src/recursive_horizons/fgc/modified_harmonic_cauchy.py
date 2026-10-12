"""Conditional boundary-free Cauchy preservation of the REF1 gauge vector.

CON2 supplies the exact metric-derived homogeneous subsidiary operator on the
complete REF1 metric and scalar equations.  This module checks the remaining
finite-dimensional hypotheses at the activated COMP1 witness and records the
standard normally-hyperbolic uniqueness consequence with its existence,
hypersurface, and boundary hypotheses left explicit.

It does not construct a spacetime solution, a nonzero-width initial slice, or
a constraint-preserving boundary map.
"""

from __future__ import annotations

from fractions import Fraction
from numbers import Integral
from typing import Any, Mapping, Sequence

from .modified_harmonic import auxiliary_inverse_metric
from .modified_harmonic_constraints import (
    normal_gauge_extension_map,
    physical_constraint_projections,
)
from .modified_harmonic_reference import (
    modified_harmonic_gauge_constraint,
    physical_connection_data,
)
from .reference_connection import ReferenceConnection
from .spherical_reduction import SphericalState


Q = Fraction
N = 4


def _fraction(name: str, value: Any) -> Fraction:
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an exact rational")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Integral):
        return Q(int(value))
    if isinstance(value, str) and value:
        try:
            result = Q(value)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError(f"{name} must be an exact rational") from exc
        canonical = (
            str(result.numerator)
            if result.denominator == 1
            else f"{result.numerator}/{result.denominator}"
        )
        if value != canonical:
            raise ValueError(f"{name} must use canonical rational encoding")
        return result
    raise TypeError(f"{name} must be an exact rational")


def _all_zero(value: Any) -> bool:
    if isinstance(value, (Fraction, Integral)) and not isinstance(value, bool):
        return Q(value) == 0
    if isinstance(value, str):
        return _fraction("zero tensor entry", value) == 0
    if isinstance(value, (tuple, list)):
        return all(_all_zero(item) for item in value)
    raise TypeError(f"unsupported zero-tensor entry {type(value).__name__}")


def _exact_matrix(name: str, value: Any, size: int) -> tuple[tuple[Fraction, ...], ...]:
    if not isinstance(value, (tuple, list)) or len(value) != size:
        raise ValueError(f"{name} must be {size} by {size}")
    rows = tuple(
        tuple(_fraction(f"{name}[{row}][{column}]", entry) for column, entry in enumerate(source))
        for row, source in enumerate(value)
        if isinstance(source, (tuple, list)) and len(source) == size
    )
    if len(rows) != size:
        raise ValueError(f"{name} must be {size} by {size}")
    return rows


def _validate_metric_subsidiary(record: Mapping[str, Any]) -> None:
    if not isinstance(record, Mapping):
        raise TypeError("metric subsidiary record must be a mapping")
    if record.get("classification") != "exact_local_metric_derived_REF1_gauge_subsidiary_identity_not_Cauchy_or_IBVP":
        raise ValueError("CON2 metric-derived subsidiary classification differs")
    required_true = (
        "both_scalar_equations_required",
        "full_metric_equations_required",
        "closed_linear_second_order_operator_in_metric_derived_C",
        "metric_derived_C_evaluates_second_reference_connection_derivatives",
        "local_metric_derived_gauge_subsidiary_identity_derived",
    )
    if any(record.get(name) is not True for name in required_true):
        raise ValueError("CON2 metric-derived subsidiary premise is absent")
    if record.get("identity") != "L_hat(C[g])=div(E_REF1)+E_phi*grad(phi)+E_chi*grad(chi)":
        raise ValueError("CON2 metric-derived subsidiary identity differs")
    if not _all_zero(record.get("subsidiary_identity_residual")):
        raise ValueError("CON2 metric-derived subsidiary residual is nonzero")
    nonclaims = record.get("nonclaims")
    if not isinstance(nonclaims, Mapping) or any(value is not False for value in nonclaims.values()):
        raise ValueError("CON2 metric-derived subsidiary nonclaim was promoted")


def _validate_comp1_composition(
    record: Mapping[str, Any],
    computed_normal_map: Mapping[str, Any],
) -> None:
    if not isinstance(record, Mapping):
        raise TypeError("COMP1 compatible-constraint record must be a mapping")
    certificate = record.get("compatible_constraint_certificate", record)
    if not isinstance(certificate, Mapping):
        raise ValueError("COMP1 compatible-constraint certificate is absent")
    checks = certificate.get("local_exact_checks")
    if not isinstance(checks, Mapping):
        raise ValueError("COMP1 local exact checks are absent")
    required = (
        "metric_defined_C_zero",
        "metric_defined_all_non_normal_nabla_C_components_zero",
        "nabla_r_Ct_and_Cr_zero",
        "normal_gauge_extension_map_nonsingular",
        "unredefined_H_zero",
        "unredefined_M_zero",
        "unredefined_Ett_and_Etr_zero_for_every_acceleration",
        "conditional_qift_root_implies_normal_gauge_compatibility",
    )
    if any(checks.get(name) is not True for name in required):
        raise ValueError("COMP1 Cauchy-compatibility premise is absent")
    external = certificate.get("normal_gauge_extension_map")
    if not isinstance(external, Mapping):
        raise ValueError("COMP1 normal gauge-extension map is absent")
    external_matrix = _exact_matrix("COMP1 normal gauge-extension map", external.get("matrix"), 2)
    if external_matrix != computed_normal_map.get("matrix"):
        raise ValueError("COMP1 stored and recomputed normal gauge-extension maps differ")
    if _fraction("COMP1 normal map determinant", external.get("determinant")) != computed_normal_map.get("determinant"):
        raise ValueError("COMP1 stored and recomputed normal map determinants differ")
    composition = record.get("theorem_composition", certificate.get("theorem_composition"))
    if not isinstance(composition, Mapping):
        raise ValueError("COMP1 theorem composition is absent")
    if composition.get("full_mhg_root_implies_normal_nabla_C_zero") is not True:
        raise ValueError("COMP1 QIFT-root normal-gauge implication is absent")
    if composition.get("gauge_extension_zero_at_full_mhg_root") is not True:
        raise ValueError("COMP1 QIFT-root extension implication is absent")


def conditional_gauge_cauchy_certificate(
    state: SphericalState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: Fraction | int,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
    metric_subsidiary_record: Mapping[str, Any],
    comp1_constraint_record: Mapping[str, Any],
) -> dict[str, Any]:
    """Check the exact COMP1 witness and state the conditional Cauchy theorem."""

    if not isinstance(state, SphericalState):
        raise TypeError("state must be a SphericalState")
    if not isinstance(reference, ReferenceConnection):
        raise TypeError("reference must be a ReferenceConnection")
    radius = _fraction("coordinate_radius", coordinate_radius)
    tilde_factor = _fraction("tilde_normal_factor", tilde_normal_factor)
    hat_factor = _fraction("hat_normal_factor", hat_normal_factor)
    if not Q(1) < tilde_factor < hat_factor:
        raise ValueError("Cauchy control requires 1 < tilde factor < hat factor")
    _validate_metric_subsidiary(metric_subsidiary_record)

    effective_planck_square = state.planck_mass**2 + state.beta * state.phi.value**2
    if effective_planck_square <= 0:
        raise ValueError("F must be strictly positive for the subsidiary operator")
    physical = physical_connection_data(state)
    hat_inverse = auxiliary_inverse_metric(physical.inverse_metric, hat_factor)
    hat_base_determinant = (
        hat_inverse[0][0] * hat_inverse[1][1]
        - hat_inverse[0][1] * hat_inverse[1][0]
    )
    if not (
        hat_inverse[0][0] < 0
        and hat_base_determinant < 0
        and hat_inverse[2][2] > 0
        and hat_inverse[3][3] > 0
    ):
        raise ValueError("hat auxiliary inverse is not Lorentzian with t slices spacelike")

    gauge = modified_harmonic_gauge_constraint(
        state,
        reference=reference,
        coordinate_radius=radius,
        tilde_normal_factor=tilde_factor,
    )
    if not _all_zero(gauge.constraint_up):
        raise ValueError("COMP1 witness does not have C=0")
    if not _all_zero(gauge.covariant_constraint_derivative):
        raise ValueError("COMP1 witness does not have the complete local nabla C=0 data")
    projections = physical_constraint_projections(state)
    if projections["H"] != 0 or projections["M"] != 0:
        raise ValueError("COMP1 witness does not satisfy both physical normal constraints")
    normal_map = normal_gauge_extension_map(
        state, gauge, hat_normal_factor=hat_factor
    )
    if normal_map["determinant"] == 0:
        raise ValueError("COMP1 normal gauge-extension map is singular")
    _validate_comp1_composition(comp1_constraint_record, normal_map)

    return {
        "classification": "exact_COMP1_point_witness_plus_conditional_boundary_free_normally_hyperbolic_gauge_Cauchy_uniqueness_not_solution_or_IBVP",
        "subsidiary_operator": {
            "equation_on_full_REF1_and_both_scalar_shells": "L_hat(C[g])=0",
            "principal_part": "(F/2)*hat_g^(mu nu)*nabla_mu*nabla_nu*C^alpha",
            "vector_component_principal_factor_is_scalar_hat_cone_identity": True,
            "lower_order_terms_are_linear_homogeneous_in_C_and_nabla_C": True,
            "effective_planck_square_F": effective_planck_square,
            "F_strictly_positive_at_COMP1": True,
            "hat_inverse_metric_at_COMP1": hat_inverse,
            "hat_base_inverse_determinant": hat_base_determinant,
            "hat_metric_Lorentzian_at_COMP1": True,
            "coordinate_t_slice_hat_spacelike_at_COMP1": True,
            "normally_hyperbolic_at_COMP1": True,
        },
        "compatible_point_witness": {
            "coordinate_radius": radius,
            "metric_defined_C_zero": True,
            "tangential_and_normal_nabla_C_zero_at_local_jet": True,
            "physical_H_zero": True,
            "physical_M_zero": True,
            "normal_gauge_extension_map": normal_map,
            "QIFT1_full_REF1_root_conditionally_supplies_normal_Cauchy_data": True,
            "one_local_activated_two_jet_only": True,
        },
        "conditional_theorem": {
            "statement": "On any sufficiently smooth full REF1 plus two-scalar solution, if F>0, hat_g is smooth Lorentzian, Sigma is hat-spacelike and C[g]|Sigma=0 with hat_n.nabla(C[g])|Sigma=0, standard uniqueness for the homogeneous normally hyperbolic system gives C[g]=0 in the boundary-free domain of dependence.",
            "regularity_required": "sufficient_for_metric_derived_C_and_linear_normally_hyperbolic_uniqueness_at_least_C3_metric_scalar_jet_control",
            "complete_REF1_metric_equations_required": True,
            "both_unmodified_scalar_equations_required": True,
            "boundary_free_domain_of_dependence_or_pre_boundary_region_required": True,
            "zero_C_and_zero_hat_normal_derivative_on_hypersurface_required": True,
            "conditional_boundary_free_Cauchy_uniqueness_statement_derived": True,
            "REF1_extension_then_vanishes_and_unredefined_ACT1_metric_equation_is_recovered": True,
        },
        "nonclaims": {
            "smooth_nonlinear_REF1_solution_constructed": False,
            "nonzero_width_compatible_initial_hypersurface_constructed": False,
            "physical_initial_constraint_manifold_solved": False,
            "uniform_normal_gauge_map_on_UHYP1_graph_proven": False,
            "metric_derived_incoming_gauge_boundary_map_to_main_fields_derived": False,
            "physical_constraint_propagation_proven": False,
            "constraint_preserving_ACT1_IBVP_proven": False,
            "quasilinear_local_existence_proven": False,
            "evolution_authorized": False,
            "collapse_solution_derived": False,
            "finite_invariant_transition_surface_derived": False,
            "singularity_resolution_derived": False,
        },
    }
