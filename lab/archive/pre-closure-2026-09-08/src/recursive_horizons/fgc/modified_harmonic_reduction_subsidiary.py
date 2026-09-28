"""Exact kinematic reduction-subsidiary system for the FO1 radial lift.

FO1 uses, for each field, ``U=(u,p,q)`` with the two kinematic equations

``D = partial_t u - p = 0`` and ``K = partial_t q - partial_r p = 0``

and the radial definition constraint ``C = q - partial_r u``.  Before any
physical equation is used, commutation of the exact ``t,r`` partials gives

``partial_t C + partial_r D - K = 0``.

Thus the complete reduction constraint of this one-spatial-dimensional lift
obeys ``partial_t C=0`` on the two declared kinematic equations.  There is no
independent spatial curl constraint in one dimension and the subsidiary
principal speed is zero, so this kinematic constraint supplies no incoming
boundary field.  Gauge and physical Einstein constraints are separate and
are intentionally not promoted by this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from typing import Any, Mapping

from .spherical_reduction import BASE_FIELD_ORDER


Q = Fraction


def _fraction(name: str, value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational")
    return value if isinstance(value, Fraction) else Q(int(value))


@dataclass(frozen=True, slots=True)
class ReductionDifferentialFieldJet:
    """The exact slots needed by the off-shell differential identity."""

    u_t: Fraction
    p: Fraction
    q: Fraction
    u_r: Fraction
    q_t: Fraction
    p_r: Fraction
    u_tr: Fraction

    def __post_init__(self) -> None:
        for name in ("u_t", "p", "q", "u_r", "q_t", "p_r", "u_tr"):
            object.__setattr__(self, name, _fraction(name, getattr(self, name)))

    @classmethod
    def from_mapping(cls, value: Mapping[str, object]) -> "ReductionDifferentialFieldJet":
        order = ("u_t", "p", "q", "u_r", "q_t", "p_r", "u_tr")
        if not isinstance(value, Mapping) or set(value) != set(order):
            raise ValueError("reduction differential jet requires exactly seven named slots")
        return cls(*(_fraction(name, value[name]) for name in order))


@dataclass(frozen=True, slots=True)
class ReductionDifferentialState:
    """Six reduction differential jets in the frozen FO1 field order."""

    h_tt: ReductionDifferentialFieldJet
    h_tr: ReductionDifferentialFieldJet
    h_rr: ReductionDifferentialFieldJet
    areal_radius: ReductionDifferentialFieldJet
    phi: ReductionDifferentialFieldJet
    chi: ReductionDifferentialFieldJet

    def __post_init__(self) -> None:
        for field in BASE_FIELD_ORDER:
            if not isinstance(getattr(self, field), ReductionDifferentialFieldJet):
                raise TypeError(f"{field} must be a ReductionDifferentialFieldJet")

    @classmethod
    def from_mapping(cls, value: Mapping[str, Mapping[str, object]]) -> "ReductionDifferentialState":
        if not isinstance(value, Mapping) or set(value) != set(BASE_FIELD_ORDER):
            raise ValueError("reduction differential state must cover exactly six fields")
        return cls(
            **{
                field: ReductionDifferentialFieldJet.from_mapping(value[field])
                for field in BASE_FIELD_ORDER
            }
        )


def reduction_subsidiary_identity(state: ReductionDifferentialState) -> dict[str, Any]:
    """Return the exact off-shell FO1 reduction differential identity."""

    if not isinstance(state, ReductionDifferentialState):
        raise TypeError("state must be a ReductionDifferentialState")
    definition = []
    radial_constraint = []
    compatibility = []
    radial_definition_derivative = []
    constraint_time_derivative = []
    identity = []
    for field in BASE_FIELD_ORDER:
        jet = getattr(state, field)
        D = jet.u_t - jet.p
        C = jet.q - jet.u_r
        K = jet.q_t - jet.p_r
        d_r_D = jet.u_tr - jet.p_r
        d_t_C = jet.q_t - jet.u_tr
        definition.append(D)
        radial_constraint.append(C)
        compatibility.append(K)
        radial_definition_derivative.append(d_r_D)
        constraint_time_derivative.append(d_t_C)
        identity.append(d_t_C + d_r_D - K)
    exact = all(value == 0 for value in identity)
    if not exact:
        raise ValueError("FO1 reduction subsidiary differential identity failed")
    on_kinematic_shell = all(
        definition[index] == 0
        and compatibility[index] == 0
        and radial_definition_derivative[index] == 0
        for index in range(len(BASE_FIELD_ORDER))
    )
    return {
        "classification": "exact_complete_kinematic_FO1_radial_reduction_subsidiary_identity_not_full_ACT1_constraints",
        "field_order": BASE_FIELD_ORDER,
        "u_time_definition_residual_D": tuple(definition),
        "radial_reduction_constraint_C": tuple(radial_constraint),
        "mixed_partial_compatibility_residual_K": tuple(compatibility),
        "radial_derivative_of_definition_residual": tuple(radial_definition_derivative),
        "time_derivative_of_radial_constraint": tuple(constraint_time_derivative),
        "off_shell_identity_residual_d_t_C_plus_d_r_D_minus_K": tuple(identity),
        "off_shell_identity_exact": True,
        "on_differentiable_kinematic_shell": on_kinematic_shell,
        "on_shell_propagation_residual_d_t_C": (
            tuple(constraint_time_derivative) if on_kinematic_shell else None
        ),
        "on_shell_radial_constraint_is_time_constant": on_kinematic_shell
        and all(value == 0 for value in constraint_time_derivative),
        "spatial_dimension": 1,
        "independent_spatial_curl_reduction_constraints": 0,
        "subsidiary_radial_principal_matrix": tuple(
            tuple(Q(0) for _ in BASE_FIELD_ORDER) for _ in BASE_FIELD_ORDER
        ),
        "subsidiary_characteristic_speeds": (Q(0),) * len(BASE_FIELD_ORDER),
        "incoming_reduction_constraint_fields_at_inner_boundary": (),
        "incoming_reduction_constraint_fields_at_outer_boundary": (),
        "complete_kinematic_reduction_subsidiary_system_derived": True,
        "conditional_zero_constraint_statement": (
            "If C=0 initially and D=K=0 hold as differentiable kinematic equations, then C remains zero."
        ),
        "nonclaims": {
            "metric_derived_gauge_Cauchy_propagation_proven": False,
            "physical_Hamiltonian_momentum_constraint_propagation_proven": False,
            "constraint_preserving_boundary_map_to_main_fields_derived": False,
            "boundary_stability_or_Kreiss_estimate_proven": False,
            "full_ACT1_initial_boundary_value_problem_proven": False,
            "evolution_authorized": False,
            "collapse_solution_derived": False,
            "finite_invariant_transition_surface_derived": False,
        },
    }
