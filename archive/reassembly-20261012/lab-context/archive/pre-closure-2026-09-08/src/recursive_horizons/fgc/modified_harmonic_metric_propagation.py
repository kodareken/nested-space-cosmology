"""Exact local metric-derived gauge-subsidiary identity for REF1.

PROP1 evaluates the complete lower-order divergence of the REF1 gauge
extension for an independent two-jet ``C^mu``.  This module supplies the
missing differential handoff: a physical spherical *third* jet is used to
derive ``C[g]`` through second coordinate derivatives, including the second
derivative of the fixed spherical reference connection.  The resulting jet
is then consumed by PROP1 and compared with a direct exact evaluation of the
unredefined ACT1 Noether identity.

The result is a local differential identity,

``L_hat(C[g]) = div(E_REF1) + E_phi grad(phi) + E_chi grad(chi)``.

Consequently the metric-derived gauge vector obeys the homogeneous linear
subsidiary equation on the full metric and both scalar equations.  This is
not, by itself, a Cauchy uniqueness theorem, a reduction-constraint closure,
or an initial-boundary value problem.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from typing import Any, Mapping

from .exact_tangent import FirstTangent, primal_and_tangent
from .modified_harmonic_propagation import (
    GaugeConstraintJet2,
    VectorJet2,
    reference_gauge_propagation_operator,
)
from .modified_harmonic_reference import (
    ModifiedHarmonicGaugeData,
    PhysicalConnectionData,
    _auxiliary_inverse_with_derivative,
    modified_harmonic_gauge_constraint,
    physical_connection_data,
)
from .reference_connection import (
    ReferenceConnection,
    ReferenceConnectionData,
    flat_spherical_annulus_connection,
)
from .reference_connection_second import (
    flat_spherical_annulus_connection_second_derivative,
)
from .spherical_reduction import (
    BASE_FIELD_ORDER,
    Jet2,
    SphericalState,
    residuals,
)
from .third_jet import THIRD_DERIVATIVE_ORDER, Jet3


Q = Fraction
N = 4
BASE_DIRECTIONS = (0, 1)
JET2_COMPONENT_ORDER = ("value", "dt", "dr", "dtt", "dtr", "drr")


def _fraction(name: str, value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational")
    return value if isinstance(value, Fraction) else Q(int(value))


def _freeze(value: Any) -> Any:
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value


def _zero(rank: int, *, dimension: int = N) -> Any:
    if rank == 1:
        return [Q(0) for _ in range(dimension)]
    return [_zero(rank - 1, dimension=dimension) for _ in range(dimension)]


def _all_zero(value: Any) -> bool:
    if isinstance(value, Fraction):
        return value == 0
    if isinstance(value, (tuple, list)):
        return all(_all_zero(item) for item in value)
    raise TypeError(f"expected an exact rational tensor, got {type(value).__name__}")


def _tensor_parts(value: Any) -> tuple[Any, Any]:
    """Split a nested exact tangent tensor into primal and tangent tensors."""

    if isinstance(value, (FirstTangent, Fraction, Integral)) and not isinstance(value, bool):
        return primal_and_tangent(value)
    if isinstance(value, (tuple, list)):
        parts = tuple(_tensor_parts(item) for item in value)
        return tuple(item[0] for item in parts), tuple(item[1] for item in parts)
    raise TypeError(f"unsupported tangent tensor scalar {type(value).__name__}")


@dataclass(frozen=True, slots=True)
class SphericalThirdJetState:
    """Six exact scalar third jets plus the frozen ACT1 parameters."""

    h_tt: Jet3
    h_tr: Jet3
    h_rr: Jet3
    areal_radius: Jet3
    phi: Jet3
    chi: Jet3
    planck_mass: Fraction = Q(1)
    beta: Fraction = Q(0)
    mu: Fraction = Q(1)
    g4: Fraction = Q(1)
    alpha: Fraction = Q(0)
    eta: Fraction = Q(0)
    branch: str = "FGC-QR"

    def __post_init__(self) -> None:
        for field in BASE_FIELD_ORDER:
            if not isinstance(getattr(self, field), Jet3):
                raise TypeError(f"{field} must be a Jet3")
        for name in ("planck_mass", "beta", "mu", "g4", "alpha", "eta"):
            object.__setattr__(self, name, _fraction(name, getattr(self, name)))
        # Reuse RED1's exact branch, signature, radius, and coupling checks.
        self.as_spherical_state()

    @classmethod
    def from_spherical_state(
        cls,
        state: SphericalState,
        third_derivatives: Mapping[str, Mapping[str, object]],
    ) -> "SphericalThirdJetState":
        """Attach a complete exact third-derivative table to one two-jet."""

        if not isinstance(state, SphericalState):
            raise TypeError("state must be a SphericalState")
        if not isinstance(third_derivatives, Mapping) or set(third_derivatives) != set(BASE_FIELD_ORDER):
            raise ValueError("third derivatives must cover exactly the six base fields")
        fields: dict[str, Jet3] = {}
        for field in BASE_FIELD_ORDER:
            supplied = third_derivatives[field]
            if not isinstance(supplied, Mapping) or set(supplied) != set(THIRD_DERIVATIVE_ORDER):
                raise ValueError(f"{field} must supply exactly the four symmetric third derivatives")
            base = getattr(state, field)
            fields[field] = Jet3(
                *(getattr(base, name) for name in JET2_COMPONENT_ORDER),
                *(_fraction(f"{field}.{name}", supplied[name]) for name in THIRD_DERIVATIVE_ORDER),
            )
        return cls(
            **fields,
            planck_mass=state.planck_mass,
            beta=state.beta,
            mu=state.mu,
            g4=state.g4,
            alpha=state.alpha,
            eta=state.eta,
            branch=state.branch,
        )

    def as_spherical_state(self) -> SphericalState:
        """Return the exact RED1 two-jet truncation."""

        return SphericalState(
            **{field: getattr(self, field).as_jet2() for field in BASE_FIELD_ORDER},
            planck_mass=self.planck_mass,
            beta=self.beta,
            mu=self.mu,
            g4=self.g4,
            alpha=self.alpha,
            eta=self.eta,
            branch=self.branch,
        )

    def directional_tangent_state(self, direction: str) -> SphericalState:
        """Seed the exact coordinate derivative of every RED1 two-jet slot."""

        if direction not in {"t", "r"}:
            raise ValueError("third-jet direction must be t or r")
        fields: dict[str, Jet2] = {}
        for field in BASE_FIELD_ORDER:
            third = getattr(self, field)
            primal = third.as_jet2()
            tangent = third.directional_jet2(direction)
            fields[field] = Jet2(
                *(
                    FirstTangent(getattr(primal, name), getattr(tangent, name))
                    for name in JET2_COMPONENT_ORDER
                )
            )
        return SphericalState(
            **fields,
            planck_mass=self.planck_mass,
            beta=self.beta,
            mu=self.mu,
            g4=self.g4,
            alpha=self.alpha,
            eta=self.eta,
            branch=self.branch,
        )


@dataclass(frozen=True, slots=True)
class PhysicalConnectionSecondData:
    """Base physical connection plus its exact ``t,r`` second partials."""

    physical: PhysicalConnectionData
    inverse_metric_second_partial: tuple[Any, ...]
    christoffel_second_partial: tuple[Any, ...]
    mixed_partials_agree_exactly: bool


@dataclass(frozen=True, slots=True)
class MetricDerivedGaugeJetData:
    """The exact metric-derived two-jet of spherical ``C^mu``."""

    gauge: ModifiedHarmonicGaugeData
    physical_second: PhysicalConnectionSecondData
    reference: ReferenceConnectionData
    tilde_inverse_metric_second_partial: tuple[Any, ...]
    gamma_difference_second_partial: tuple[Any, ...]
    reference_second_connection_constraint_contribution: tuple[Any, ...]
    constraint_second_partial_up: tuple[Any, ...]
    gauge_jet: GaugeConstraintJet2
    mixed_partials_agree_exactly: bool
    reference_second_connection_contribution_zero_by_spherical_contraction: bool


def physical_connection_second_data(state: SphericalThirdJetState) -> PhysicalConnectionSecondData:
    """Differentiate RED1's exact physical connection in both base directions."""

    if not isinstance(state, SphericalThirdJetState):
        raise TypeError("state must be a SphericalThirdJetState")
    base = physical_connection_data(state.as_spherical_state())
    inverse_second = [[[[Q(0) for _ in range(N)] for _ in range(N)] for _ in BASE_DIRECTIONS] for _ in BASE_DIRECTIONS]
    gamma_second = [[[[[Q(0) for _ in range(N)] for _ in range(N)] for _ in range(N)] for _ in BASE_DIRECTIONS] for _ in BASE_DIRECTIONS]
    for direction, name in zip(BASE_DIRECTIONS, ("t", "r"), strict=True):
        tangent_physical = physical_connection_data(state.directional_tangent_state(name))
        inverse_primal, _ = _tensor_parts(tangent_physical.inverse_metric)
        gamma_primal, _ = _tensor_parts(tangent_physical.christoffel)
        if inverse_primal != base.inverse_metric or gamma_primal != base.christoffel:
            raise ValueError("directional physical connection changed the declared primal state")
        derivative_primal, derivative_tangent = _tensor_parts(tangent_physical.inverse_metric_derivative)
        gamma_derivative_primal, gamma_derivative_tangent = _tensor_parts(tangent_physical.christoffel_derivative)
        if derivative_primal != base.inverse_metric_derivative or gamma_derivative_primal != base.christoffel_derivative:
            raise ValueError("directional physical connection derivative changed its primal")
        for second_direction in BASE_DIRECTIONS:
            inverse_second[direction][second_direction] = derivative_tangent[second_direction]
            gamma_second[direction][second_direction] = gamma_derivative_tangent[second_direction]
    mixed = (
        inverse_second[0][1] == inverse_second[1][0]
        and gamma_second[0][1] == gamma_second[1][0]
    )
    if not mixed:
        raise ValueError("physical connection mixed coordinate partials disagree")
    return PhysicalConnectionSecondData(
        physical=base,
        inverse_metric_second_partial=_freeze(inverse_second),
        christoffel_second_partial=_freeze(gamma_second),
        mixed_partials_agree_exactly=True,
    )


def metric_derived_gauge_jet(
    state: SphericalThirdJetState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: Fraction | int,
    tilde_normal_factor: Fraction | int,
) -> MetricDerivedGaugeJetData:
    """Derive ``C[g]`` through its complete spherical ``t,r`` two-jet."""

    if not isinstance(state, SphericalThirdJetState):
        raise TypeError("state must be a SphericalThirdJetState")
    if not isinstance(reference, ReferenceConnection):
        raise TypeError("reference must be a ReferenceConnection")
    radius = _fraction("coordinate_radius", coordinate_radius)
    factor = _fraction("tilde_normal_factor", tilde_normal_factor)
    if factor <= 1:
        raise ValueError("tilde normal factor must exceed one")
    base_state = state.as_spherical_state()
    gauge = modified_harmonic_gauge_constraint(
        base_state,
        reference=reference,
        coordinate_radius=radius,
        tilde_normal_factor=factor,
    )
    physical_second = physical_connection_second_data(state)
    reference_data = flat_spherical_annulus_connection(
        reference, coordinate_radius=radius
    )
    reference_second_derivative = (
        flat_spherical_annulus_connection_second_derivative(
            reference, coordinate_radius=radius
        )
    )

    tilde_second = [[[[Q(0) for _ in range(N)] for _ in range(N)] for _ in BASE_DIRECTIONS] for _ in BASE_DIRECTIONS]
    for direction, name in zip(BASE_DIRECTIONS, ("t", "r"), strict=True):
        tangent_physical = physical_connection_data(state.directional_tangent_state(name))
        tangent_tilde, tangent_tilde_derivative = _auxiliary_inverse_with_derivative(
            tangent_physical, factor
        )
        tilde_primal, _ = _tensor_parts(tangent_tilde)
        derivative_primal, derivative_tangent = _tensor_parts(tangent_tilde_derivative)
        if tilde_primal != gauge.tilde_inverse_metric or derivative_primal != gauge.tilde_inverse_metric_derivative:
            raise ValueError("directional auxiliary inverse changed the declared primal")
        for second_direction in BASE_DIRECTIONS:
            tilde_second[direction][second_direction] = derivative_tangent[second_direction]

    gamma_second = [[[[[Q(0) for _ in range(N)] for _ in range(N)] for _ in range(N)] for _ in BASE_DIRECTIONS] for _ in BASE_DIRECTIONS]
    for first in BASE_DIRECTIONS:
        for second in BASE_DIRECTIONS:
            for upper in range(N):
                for lower_one in range(N):
                    for lower_two in range(N):
                        gamma_second[first][second][upper][lower_one][lower_two] = (
                            physical_second.christoffel_second_partial[first][second][upper][lower_one][lower_two]
                            - reference_second_derivative[first][second][upper][lower_one][lower_two]
                        )

    reference_second_contribution = [[[Q(0) for _ in range(N)] for _ in BASE_DIRECTIONS] for _ in BASE_DIRECTIONS]
    constraint_second = [[[Q(0) for _ in range(N)] for _ in BASE_DIRECTIONS] for _ in BASE_DIRECTIONS]
    for first in BASE_DIRECTIONS:
        for second in BASE_DIRECTIONS:
            for upper in range(N):
                # This is the exact +tilde^rho_sigma*d_first*d_second
                # bar_Gamma^upper_(rho sigma) term inside d^2 C.  For the
                # declared spherical C^t,C^r sector the real flat-reference
                # table contracts to zero, but it is evaluated explicitly so
                # that the cancellation is proved rather than assumed.
                reference_second_contribution[first][second][upper] = sum(
                    gauge.tilde_inverse_metric[rho][sigma]
                    * reference_second_derivative[first][second][upper][rho][sigma]
                    for rho in range(N)
                    for sigma in range(N)
                )
                constraint_second[first][second][upper] = -sum(
                    tilde_second[first][second][rho][sigma]
                    * gauge.gamma_difference[upper][rho][sigma]
                    + gauge.tilde_inverse_metric_derivative[second][rho][sigma]
                    * gauge.gamma_difference_derivative[first][upper][rho][sigma]
                    + gauge.tilde_inverse_metric_derivative[first][rho][sigma]
                    * gauge.gamma_difference_derivative[second][upper][rho][sigma]
                    + gauge.tilde_inverse_metric[rho][sigma]
                    * gamma_second[first][second][upper][rho][sigma]
                    for rho in range(N)
                    for sigma in range(N)
                )

    mixed = (
        tilde_second[0][1] == tilde_second[1][0]
        and gamma_second[0][1] == gamma_second[1][0]
        and constraint_second[0][1] == constraint_second[1][0]
    )
    if not mixed:
        raise ValueError("metric-derived gauge mixed coordinate partials disagree")
    if not _all_zero(gauge.constraint_up[2:]) or not _all_zero(
        tuple(row[2:] for row in gauge.partial_constraint_up)
    ) or not _all_zero(tuple(tuple(row[2:] for row in plane) for plane in constraint_second)):
        raise ValueError("metric-derived gauge jet leaves the spherical C^t,C^r sector")

    components = []
    for upper in BASE_DIRECTIONS:
        components.append(
            VectorJet2(
                gauge.constraint_up[upper],
                gauge.partial_constraint_up[0][upper],
                gauge.partial_constraint_up[1][upper],
                constraint_second[0][0][upper],
                constraint_second[0][1][upper],
                constraint_second[1][1][upper],
            )
        )
    return MetricDerivedGaugeJetData(
        gauge=gauge,
        physical_second=physical_second,
        reference=reference_data,
        tilde_inverse_metric_second_partial=_freeze(tilde_second),
        gamma_difference_second_partial=_freeze(gamma_second),
        reference_second_connection_constraint_contribution=_freeze(reference_second_contribution),
        constraint_second_partial_up=_freeze(constraint_second),
        gauge_jet=GaugeConstraintJet2(*components),
        mixed_partials_agree_exactly=True,
        reference_second_connection_contribution_zero_by_spherical_contraction=_all_zero(
            reference_second_contribution
        ),
    )


def unredefined_noether_identity(state: SphericalThirdJetState) -> dict[str, Any]:
    """Evaluate the exact ACT1 Noether identity on one spherical third jet."""

    if not isinstance(state, SphericalThirdJetState):
        raise TypeError("state must be a SphericalThirdJetState")
    base_state = state.as_spherical_state()
    source = residuals(base_state)
    physical = physical_connection_data(base_state)
    metric_partial = [_zero(2) for _ in BASE_DIRECTIONS]
    scalar_partial: dict[str, list[Fraction]] = {"phi": [], "chi": []}
    for direction, name in zip(BASE_DIRECTIONS, ("t", "r"), strict=True):
        tangent_source = residuals(state.directional_tangent_state(name))
        metric_primal, metric_tangent = _tensor_parts(tangent_source["metric"])
        if metric_primal != source["metric"]:
            raise ValueError("directional unredefined residual changed its primal")
        metric_partial[direction] = metric_tangent
        for scalar in ("phi", "chi"):
            primal, tangent = primal_and_tangent(tangent_source[scalar])
            if primal != source[scalar]:
                raise ValueError("directional scalar residual changed its primal")
            scalar_partial[scalar].append(tangent)

    metric = source["metric"]
    if any(metric[a][b] != 0 for a in range(N) for b in range(N) if (a, b) not in {(0, 0), (0, 1), (1, 0), (1, 1), (2, 2), (3, 3)}):
        raise ValueError("unredefined residual leaves the spherical tensor sector")
    if metric[2][2] != metric[3][3]:
        raise ValueError("equatorial angular residual components disagree")

    divergence_down = [Q(0) for _ in range(N)]
    inverse = physical.inverse_metric
    gamma = physical.christoffel
    for covector in range(N):
        divergence_down[covector] = sum(
            inverse[metric_index][direction]
            * (
                (metric_partial[direction][metric_index][covector] if direction in BASE_DIRECTIONS else Q(0))
                - sum(
                    gamma[sigma][direction][metric_index] * metric[sigma][covector]
                    + gamma[sigma][direction][covector] * metric[metric_index][sigma]
                    for sigma in range(N)
                )
            )
            for metric_index in range(N)
            for direction in range(N)
        )
    scalar_source_down = [Q(0) for _ in range(N)]
    for direction in BASE_DIRECTIONS:
        field_derivative = (
            getattr(base_state.phi, "dt" if direction == 0 else "dr"),
            getattr(base_state.chi, "dt" if direction == 0 else "dr"),
        )
        scalar_source_down[direction] = (
            source["phi"] * field_derivative[0]
            + source["chi"] * field_derivative[1]
        )
    divergence_up = tuple(
        sum(inverse[upper][lower] * divergence_down[lower] for lower in range(N))
        for upper in range(N)
    )
    scalar_source_up = tuple(
        sum(inverse[upper][lower] * scalar_source_down[lower] for lower in range(N))
        for upper in range(N)
    )
    exact_down = tuple(
        divergence_down[index] + scalar_source_down[index] for index in range(N)
    )
    exact_up = tuple(
        divergence_up[index] + scalar_source_up[index] for index in range(N)
    )
    if not _all_zero(exact_down) or not _all_zero(exact_up):
        raise ValueError("ACT1 unredefined Noether identity failed on the supplied third jet")
    return {
        "metric_residual_covariant": metric,
        "metric_residual_t_r_partial": _freeze(metric_partial),
        "scalar_residuals": (source["phi"], source["chi"]),
        "scalar_residual_t_r_partial": _freeze(scalar_partial),
        "divergence_down": tuple(divergence_down),
        "divergence_up": divergence_up,
        "scalar_noether_source_down": tuple(scalar_source_down),
        "scalar_noether_source_up": scalar_source_up,
        "divergence_plus_scalar_source_down": exact_down,
        "divergence_plus_scalar_source_up": exact_up,
        "angular_partial_rule": "spherical tensor coordinate partials vanish at theta=pi/2; connection terms retained",
        "unredefined_ACT1_noether_identity_exact": True,
    }


def metric_derived_gauge_propagation_identity(
    state: SphericalThirdJetState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: Fraction | int,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
) -> dict[str, Any]:
    """Close the exact local REF1 gauge-subsidiary differential identity."""

    tilde = _fraction("tilde_normal_factor", tilde_normal_factor)
    hat = _fraction("hat_normal_factor", hat_normal_factor)
    if not Q(1) < tilde < hat:
        raise ValueError("metric propagation requires 1 < tilde factor < hat factor")
    gauge = metric_derived_gauge_jet(
        state,
        reference=reference,
        coordinate_radius=coordinate_radius,
        tilde_normal_factor=tilde,
    )
    noether = unredefined_noether_identity(state)
    operator = reference_gauge_propagation_operator(
        state.as_spherical_state(),
        gauge.gauge_jet,
        reference=reference,
        coordinate_radius=coordinate_radius,
        hat_normal_factor=hat,
    )
    if not operator["two_routes_agree_exactly"]:
        raise ValueError("metric-derived gauge operator routes disagree")
    full_divergence = tuple(
        noether["divergence_up"][index] + operator["coordinate_divergence"][index]
        for index in range(N)
    )
    subsidiary_identity = tuple(
        operator["coordinate_divergence"][index]
        - full_divergence[index]
        - noether["scalar_noether_source_up"][index]
        for index in range(N)
    )
    if not _all_zero(subsidiary_identity):
        raise ValueError("metric-derived gauge subsidiary identity failed")
    return {
        "classification": "exact_local_metric_derived_REF1_gauge_subsidiary_identity_not_Cauchy_or_IBVP",
        "third_jet_state": state,
        "metric_derived_gauge": gauge,
        "unredefined_noether": noether,
        "metric_derived_extension_operator": operator,
        "full_modified_metric_divergence_up": full_divergence,
        "subsidiary_identity_residual": subsidiary_identity,
        "identity": "L_hat(C[g])=div(E_REF1)+E_phi*grad(phi)+E_chi*grad(chi)",
        "both_scalar_equations_required": True,
        "full_metric_equations_required": True,
        "closed_linear_second_order_operator_in_metric_derived_C": True,
        "metric_derived_C_evaluates_second_reference_connection_derivatives": True,
        "second_reference_connection_contribution_zero_in_declared_spherical_sector": gauge.reference_second_connection_contribution_zero_by_spherical_contraction,
        "local_metric_derived_gauge_subsidiary_identity_derived": True,
        "conditional_homogeneous_subsidiary_equation": (
            "On E_REF1=0, E_phi=0, and E_chi=0, the exact identity gives L_hat(C[g])=0."
        ),
        "nonclaims": {
            "Cauchy_uniqueness_and_zero_constraint_propagation_proven": False,
            "complete_first_order_reduction_constraint_system_proven": False,
            "physical_initial_constraint_hypersurface_solved": False,
            "constraint_preserving_boundary_conditions_derived": False,
            "boundary_stability_or_Kreiss_estimate_proven": False,
            "initial_boundary_value_problem_proven": False,
            "evolution_authorized": False,
            "collapse_solution_derived": False,
            "finite_invariant_transition_surface_derived": False,
            "singularity_resolution_derived": False,
        },
    }
