"""Exact lower-order gauge-propagation operator for the REF1 formulation.

This module evaluates the divergence of the REF1 gauge extension for an
*independent* spherically symmetric contravariant gauge-vector two-jet.  It
does not differentiate the metric-defined connection difference through a
third metric jet and therefore does not claim a direct full-residual-divergence
or a reduction-constraint propagation theorem.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from typing import Any, Mapping

from .exact_linear_algebra import poly
from .modified_harmonic import (
    _trace_reversal_projector,
    auxiliary_inverse_metric,
    gauge_constraint_propagation_principal_symbol,
)
from .modified_harmonic_implicit import (
    modified_harmonic_implicit_certificate,
)
from .modified_harmonic_reference import PhysicalConnectionData, physical_connection_data
from .reference_connection import (
    SPHERICAL_FLAT_REFERENCE_ID,
    ReferenceConnection,
)
from .spherical_reduction import (
    SphericalState,
    _ricci,
    direct_4d_curvature,
    state_from_generalized_adm_pg_fixture,
)


Q = Fraction
N = 4
GAUGE_VECTOR_ORDER = ("C^t", "C^r", "C^theta", "C^phi")
SPHERICAL_GAUGE_VECTOR_ORDER = GAUGE_VECTOR_ORDER[:2]
SECOND_JET_ORDER = ("dtt", "dtr", "drr")


def _q(name: str, value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational value")
    return value if isinstance(value, Fraction) else Q(int(value))


def _zero(rank: int):
    if rank == 1:
        return [Q(0) for _ in range(N)]
    return [_zero(rank - 1) for _ in range(N)]


def _freeze(value: Any) -> Any:
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value


@dataclass(frozen=True, slots=True)
class VectorJet2:
    """One exact scalar component through the declared ``t,r`` two-jet."""

    value: Fraction
    dt: Fraction = Q(0)
    dr: Fraction = Q(0)
    dtt: Fraction = Q(0)
    dtr: Fraction = Q(0)
    drr: Fraction = Q(0)

    def __post_init__(self) -> None:
        for name in ("value", "dt", "dr", "dtt", "dtr", "drr"):
            object.__setattr__(self, name, _q(name, getattr(self, name)))


@dataclass(frozen=True, slots=True)
class GaugeConstraintJet2:
    """Independent spherical contravariant ``C^mu`` two-jet.

    Angular components and their coordinate derivatives are fixed to zero.
    Their *covariant* derivatives are nevertheless retained through the full
    physical spherical connection.
    """

    t: VectorJet2
    r: VectorJet2

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "GaugeConstraintJet2":
        if not isinstance(value, Mapping) or set(value) != {"C^t", "C^r"}:
            raise ValueError("gauge constraint requires exactly C^t and C^r tables")

        def component(name: str) -> VectorJet2:
            raw = value[name]
            if not isinstance(raw, Mapping):
                raise ValueError(f"{name} must be a mapping")
            allowed = {"value", "dt", "dr", "dtt", "dtr", "drr"}
            if set(raw) - allowed or "value" not in raw:
                raise ValueError(f"{name} has invalid two-jet keys")
            return VectorJet2(**{key: _q(f"{name}.{key}", raw.get(key, 0)) for key in allowed})

        return cls(t=component("C^t"), r=component("C^r"))

    def components(self) -> tuple[VectorJet2, VectorJet2, VectorJet2, VectorJet2]:
        return (self.t, self.r, VectorJet2(0), VectorJet2(0))


def covariant_vector_derivatives(
    gauge: GaugeConstraintJet2, physical: PhysicalConnectionData
) -> dict[str, Any]:
    """Return exact first and second physical covariant derivatives of ``C``."""

    if not isinstance(gauge, GaugeConstraintJet2):
        raise TypeError("gauge must be a GaugeConstraintJet2")
    if not isinstance(physical, PhysicalConnectionData):
        raise TypeError("physical must be PhysicalConnectionData")
    jets = gauge.components()
    values = [jet.value for jet in jets]
    partial = [[Q(0) for _ in range(N)] for _ in range(N)]
    second_partial = _zero(3)
    for component, jet in enumerate(jets):
        partial[0][component], partial[1][component] = jet.dt, jet.dr
        second_partial[0][0][component] = jet.dtt
        second_partial[0][1][component] = second_partial[1][0][component] = jet.dtr
        second_partial[1][1][component] = jet.drr

    first = _zero(2)
    for direction in range(N):
        for upper in range(N):
            first[direction][upper] = partial[direction][upper] + sum(
                physical.christoffel[upper][direction][sigma] * values[sigma]
                for sigma in range(N)
            )

    second = _zero(3)
    for mu in range(N):
        for beta in range(N):
            for upper in range(N):
                partial_first = second_partial[mu][beta][upper] + sum(
                    physical.christoffel_derivative[mu][upper][beta][sigma]
                    * values[sigma]
                    + physical.christoffel[upper][beta][sigma]
                    * partial[mu][sigma]
                    for sigma in range(N)
                )
                second[mu][beta][upper] = partial_first + sum(
                    physical.christoffel[upper][mu][sigma] * first[beta][sigma]
                    - physical.christoffel[sigma][mu][beta] * first[sigma][upper]
                    for sigma in range(N)
                )
    return {
        "values": tuple(values),
        "partial": _freeze(partial),
        "second_partial": _freeze(second_partial),
        "covariant_first": _freeze(first),
        "covariant_second": _freeze(second),
    }


def _auxiliary_with_derivative(
    physical: PhysicalConnectionData, normal_factor: Fraction
) -> tuple[Any, Any]:
    inverse, derivative = physical.inverse_metric, physical.inverse_metric_derivative
    hat = auxiliary_inverse_metric(inverse, normal_factor)
    normal = [[-inverse[a][0] * inverse[b][0] / inverse[0][0] for b in range(N)] for a in range(N)]
    normal_d = _zero(3)
    for direction in range(N):
        denom, ddenom = inverse[0][0], derivative[direction][0][0]
        for a in range(N):
            for b in range(N):
                numerator = inverse[a][0] * inverse[b][0]
                dnumerator = derivative[direction][a][0] * inverse[b][0] + inverse[a][0] * derivative[direction][b][0]
                normal_d[direction][a][b] = -dnumerator / denom + numerator * ddenom / denom**2
    hat_d = _zero(3)
    for direction in range(N):
        for a in range(N):
            for b in range(N):
                hat_d[direction][a][b] = derivative[direction][a][b] - (normal_factor - 1) * normal_d[direction][a][b]
    return hat, _freeze(hat_d)


def _projector_and_derivative(hat: Any, hat_d: Any) -> tuple[Any, Any]:
    projector = _zero(4)
    derivative = _zero(5)
    for alpha in range(N):
        for beta in range(N):
            for mu in range(N):
                for nu in range(N):
                    projector[alpha][beta][mu][nu] = _trace_reversal_projector(hat, alpha, beta, mu, nu)
                    for direction in range(N):
                        derivative[direction][alpha][beta][mu][nu] = (
                            Q(int(alpha == mu)) * hat_d[direction][nu][beta]
                            + Q(int(alpha == nu)) * hat_d[direction][mu][beta]
                            - Q(int(alpha == beta)) * hat_d[direction][mu][nu]
                        ) / 2
    return _freeze(projector), _freeze(derivative)


def _covariant_projector_derivative(projector: Any, partial: Any, physical: PhysicalConnectionData) -> Any:
    output = _zero(5)
    gamma = physical.christoffel
    for direction in range(N):
        for alpha in range(N):
            for beta in range(N):
                for mu in range(N):
                    for nu in range(N):
                        output[direction][alpha][beta][mu][nu] = partial[direction][alpha][beta][mu][nu] + sum(
                            -gamma[rho][direction][alpha] * projector[rho][beta][mu][nu]
                            + gamma[beta][direction][rho] * projector[alpha][rho][mu][nu]
                            + gamma[mu][direction][rho] * projector[alpha][beta][rho][nu]
                            + gamma[nu][direction][rho] * projector[alpha][beta][mu][rho]
                            for rho in range(N)
                        )
    return _freeze(output)


def reference_gauge_propagation_operator(
    state: SphericalState,
    gauge: GaugeConstraintJet2,
    *,
    reference: ReferenceConnection,
    coordinate_radius: Fraction | int,
    hat_normal_factor: Fraction | int,
) -> dict[str, Any]:
    """Assemble ``nabla_mu[F P_hat nabla C]`` by two exact routes.

    ``reference`` and ``coordinate_radius`` are validated to lock this operator
    to REF1's annular formulation; the operator input remains independent C.
    """

    if not isinstance(state, SphericalState):
        raise TypeError("state must be a SphericalState")
    if not isinstance(reference, ReferenceConnection):
        raise TypeError("reference must be a ReferenceConnection")
    if reference.reference_id != SPHERICAL_FLAT_REFERENCE_ID:
        raise ValueError("PROP1 requires REF1's flat spherical reference")
    radius, factor = _q("coordinate_radius", coordinate_radius), _q("hat_normal_factor", hat_normal_factor)
    if radius <= reference.radial_domain_minimum or factor <= 1:
        raise ValueError("PROP1 requires a positive annular radius and hat factor > 1")
    physical = physical_connection_data(state)
    derivatives = covariant_vector_derivatives(gauge, physical)
    values = derivatives["values"]
    partial = derivatives["partial"]
    second_partial = derivatives["second_partial"]
    first = derivatives["covariant_first"]
    second = derivatives["covariant_second"]
    F = state.planck_mass**2 + state.beta * state.phi.value**2
    if F <= 0:
        raise ValueError("PROP1 requires positive effective Planck coefficient")
    dF = (2 * state.beta * state.phi.value * state.phi.dt, 2 * state.beta * state.phi.value * state.phi.dr, Q(0), Q(0))
    hat, hat_d = _auxiliary_with_derivative(physical, factor)
    projector, projector_partial = _projector_and_derivative(hat, hat_d)
    cov_projector = _covariant_projector_derivative(projector, projector_partial, physical)

    # Route 1: construct S^{mu nu}, take its ordinary coordinate derivative,
    # then add the two contravariant divergence connection terms.
    extension = _zero(2)
    extension_partial = _zero(3)
    for mu in range(N):
        for nu in range(N):
            extension[mu][nu] = F * sum(projector[a][b][mu][nu] * first[b][a] for a in range(N) for b in range(N))
            for direction in range(N):
                extension_partial[direction][mu][nu] = sum(
                    dF[direction] * projector[a][b][mu][nu] * first[b][a]
                    + F * projector_partial[direction][a][b][mu][nu] * first[b][a]
                    + F
                    * projector[a][b][mu][nu]
                    * (
                        second_partial[direction][b][a]
                        + sum(
                            physical.christoffel_derivative[direction][a][b][rho]
                            * values[rho]
                            + physical.christoffel[a][b][rho]
                            * partial[direction][rho]
                            for rho in range(N)
                        )
                    )
                    for a in range(N) for b in range(N)
                )
    coordinate = _zero(1)
    for nu in range(N):
        coordinate[nu] = sum(
            extension_partial[mu][mu][nu]
            + sum(
                physical.christoffel[mu][mu][rho] * extension[rho][nu]
                + physical.christoffel[nu][mu][rho] * extension[mu][rho]
                for rho in range(N)
            )
            for mu in range(N)
        )

    # Route 2: product rule plus the exact vector commutator identity.
    riemann, inverse, _ = direct_4d_curvature(state)
    ricci = _ricci(riemann, inverse)
    expanded = _zero(1)
    wave = _zero(1)
    connection_wave = _zero(1)
    curvature = _zero(1)
    derivative_F = _zero(1)
    derivative_projector = _zero(1)
    for nu in range(N):
        wave[nu] = F * sum(hat[mu][beta] * second[mu][beta][nu] for mu in range(N) for beta in range(N)) / 2
        raw_second_wave = F * sum(
            hat[mu][beta] * derivatives["second_partial"][mu][beta][nu]
            for mu in range(N) for beta in range(N)
        ) / 2
        # Only t/r partial two-jets exist; this isolates the connection-bearing
        # portion of the covariant hat-wave term on a spherical probe.
        connection_wave[nu] = wave[nu] - raw_second_wave
        curvature[nu] = F * sum(hat[nu][beta] * ricci[sigma][beta] * values[sigma] for beta in range(N) for sigma in range(N)) / 2
        derivative_F[nu] = sum(dF[mu] * projector[a][b][mu][nu] * first[b][a] for mu in range(N) for a in range(N) for b in range(N))
        derivative_projector[nu] = F * sum(cov_projector[mu][a][b][mu][nu] * first[b][a] for mu in range(N) for a in range(N) for b in range(N))
        expanded[nu] = wave[nu] + curvature[nu] + derivative_F[nu] + derivative_projector[nu]
    return {
        "gauge_input_is_independent": True,
        "direct_metric_derived_residual_divergence_evaluated": False,
        "coordinate_route_uses_raw_partial_jet_not_covariant_second": True,
        "expanded_route_uses_covariant_second": True,
        "gauge_vector_order": GAUGE_VECTOR_ORDER,
        "physical": physical,
        "gauge_derivatives": derivatives,
        "effective_planck_coefficient": F,
        "gradient_effective_planck": tuple(dF),
        "hat_inverse_metric": hat,
        "extension_contravariant": _freeze(extension),
        "coordinate_divergence": tuple(coordinate),
        "expanded_wave_term": tuple(wave),
        "expanded_connection_wave_contribution": tuple(connection_wave),
        "expanded_curvature_commutator_term": tuple(curvature),
        "expanded_gradient_F_term": tuple(derivative_F),
        "expanded_projector_derivative_term": tuple(derivative_projector),
        "expanded_divergence": tuple(expanded),
        "two_routes_agree_exactly": tuple(coordinate) == tuple(expanded),
    }


def gauge_propagation_principal_blocks(
    state: SphericalState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: Fraction | int,
    hat_normal_factor: Fraction | int,
) -> dict[str, Any]:
    """Extract all ``t/r`` second-partial-C blocks and regress to MHG1.

    The extraction uses independent gauge-vector inputs with all values and
    first derivatives zero.  It therefore reads the complete lower-order
    operator's principal coefficients, rather than inserting MHG1's block.
    """

    coefficient_names = ("value", "dt", "dr") + SECOND_JET_ORDER
    all_blocks = {
        name: [
            [Q(0) for _ in SPHERICAL_GAUGE_VECTOR_ORDER]
            for _ in range(N)
        ]
        for name in coefficient_names
    }
    blocks = {name: all_blocks[name] for name in SECOND_JET_ORDER}
    for alpha in range(len(SPHERICAL_GAUGE_VECTOR_ORDER)):
        for derivative, kwargs in (
            ("value", {"value": Q(1)}),
            ("dt", {"dt": Q(1)}),
            ("dr", {"dr": Q(1)}),
            ("dtt", {"dtt": Q(1)}),
            ("dtr", {"dtr": Q(1)}),
            ("drr", {"drr": Q(1)}),
        ):
            components = [VectorJet2(0), VectorJet2(0)]
            components[alpha] = (
                VectorJet2(**kwargs)
                if derivative == "value"
                else VectorJet2(0, **kwargs)
            )
            data = reference_gauge_propagation_operator(
                state, GaugeConstraintJet2(*components), reference=reference,
                coordinate_radius=coordinate_radius, hat_normal_factor=hat_normal_factor,
            )
            for nu, value in enumerate(data["coordinate_divergence"]):
                all_blocks[derivative][nu][alpha] = value
    radial_polynomial = tuple(
        tuple(
            poly(
                (
                    blocks["drr"][nu][alpha],
                    -blocks["dtr"][nu][alpha],
                    blocks["dtt"][nu][alpha],
                )
            )
            for alpha in range(len(SPHERICAL_GAUGE_VECTOR_ORDER))
        )
        for nu in range(N)
    )
    expected = gauge_constraint_propagation_principal_symbol(
        state, hat_normal_factor=hat_normal_factor
    )["full_projector_contraction"]
    # The spherical independent-input restriction controls t/r columns; the
    # MHG1 four-index result remains serialized as the authority for all four.
    spherical_exact = all(
        radial_polynomial[nu][alpha] == expected[nu][alpha]
        for nu in range(N) for alpha in range(2)
    )
    return {
        "represented_C_input_order": SPHERICAL_GAUGE_VECTOR_ORDER,
        "coefficient_block_shape": (N, len(SPHERICAL_GAUGE_VECTOR_ORDER)),
        "second_partial_C_order": SECOND_JET_ORDER,
        "C_value_first_second_partial_blocks": {name: _freeze(value) for name, value in all_blocks.items()},
        "radial_principal_polynomial": radial_polynomial,
        "mhg1_full_four_index_principal": expected,
        "spherical_t_r_columns_match_mhg1_exactly": spherical_exact,
        "angular_independent_C_columns_not_represented": True,
    }


def required_prop1_nonclaims() -> dict[str, bool]:
    return {
        name: False
        for name in (
            "direct_metric_derived_full_residual_divergence_evaluated",
            "metric_derived_gauge_constraint_propagation_proven",
            "constraint_preserving_initial_boundary_value_problem_proven",
            "complete_lower_order_first_order_sources_derived",
            "reduction_constraint_propagation_proven",
            "physical_initial_constraints_solved",
            "quantified_implicit_branch_neighborhood_proven",
            "uniform_open_domain_symmetrizer_proven",
            "open_retained_eft_domain_proven",
            "boundary_or_regular_center_system_derived",
            "evolution_authorized",
            "collapse_solution_derived",
            "metric_null_affine_defocusing_derived",
            "singularity_resolution_derived",
        )
    }


def _component_values(component: VectorJet2) -> tuple[Fraction, ...]:
    return tuple(
        getattr(component, name)
        for name in ("value", "dt", "dr", "dtt", "dtr", "drr")
    )


def _combine_gauge_jets(
    left: GaugeConstraintJet2,
    right: GaugeConstraintJet2,
) -> GaugeConstraintJet2:
    def add(a: VectorJet2, b: VectorJet2) -> VectorJet2:
        return VectorJet2(
            *(x + y for x, y in zip(_component_values(a), _component_values(b)))
        )

    return GaugeConstraintJet2(add(left.t, right.t), add(left.r, right.r))


def _scale_gauge_jet(
    gauge: GaugeConstraintJet2,
    factor: Fraction,
) -> GaugeConstraintJet2:
    def scale(component: VectorJet2) -> VectorJet2:
        return VectorJet2(*(factor * value for value in _component_values(component)))

    return GaugeConstraintJet2(scale(gauge.t), scale(gauge.r))


def _independent_companion_probe(
    gauge: GaugeConstraintJet2,
) -> GaugeConstraintJet2:
    """Return a deterministic exact probe not proportional to ``gauge``."""

    t = _component_values(gauge.t)
    r = _component_values(gauge.r)
    return GaugeConstraintJet2(
        VectorJet2(r[0], -t[1], r[2], -t[3], r[4], -t[5]),
        VectorJet2(-t[0], r[1], -t[2], r[3], -t[4], r[5]),
    )


def _operator_summary(record: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "gauge_input_is_independent": record["gauge_input_is_independent"],
        "direct_metric_derived_residual_divergence_evaluated": record[
            "direct_metric_derived_residual_divergence_evaluated"
        ],
        "coordinate_route_uses_raw_partial_jet_not_covariant_second": record[
            "coordinate_route_uses_raw_partial_jet_not_covariant_second"
        ],
        "expanded_route_uses_covariant_second": record[
            "expanded_route_uses_covariant_second"
        ],
        "effective_planck_coefficient": record["effective_planck_coefficient"],
        "gradient_effective_planck": record["gradient_effective_planck"],
        "hat_inverse_metric": record["hat_inverse_metric"],
        "extension_contravariant": record["extension_contravariant"],
        "coordinate_divergence": record["coordinate_divergence"],
        "expanded_wave_term": record["expanded_wave_term"],
        "expanded_connection_wave_contribution": record[
            "expanded_connection_wave_contribution"
        ],
        "expanded_curvature_commutator_term": record[
            "expanded_curvature_commutator_term"
        ],
        "expanded_gradient_F_term": record["expanded_gradient_F_term"],
        "expanded_projector_derivative_term": record[
            "expanded_projector_derivative_term"
        ],
        "expanded_divergence": record["expanded_divergence"],
        "two_routes_agree_exactly": record["two_routes_agree_exactly"],
    }


def _tensor_nonzero(value: Any) -> bool:
    if isinstance(value, Fraction):
        return value != 0
    if isinstance(value, (tuple, list)):
        return any(_tensor_nonzero(item) for item in value)
    raise ValueError("PROP1 tensor contains an unsupported scalar")


def modified_harmonic_propagation_certificate(
    configuration: Mapping[str, Any],
) -> dict[str, Any]:
    """Build the exact conditional lower-order reference-gauge certificate."""

    expected_keys = {
        "flat_fixture",
        "activated_fixture",
        "flat_probe",
        "activated_probe",
        "reference",
        "flat_coordinate_radius",
        "activated_coordinate_radius",
        "tilde_normal_factor",
        "hat_normal_factor",
    }
    if not isinstance(configuration, Mapping) or set(configuration) != expected_keys:
        raise ValueError("PROP1 configuration adapter has unexpected keys")
    flat_probe = configuration["flat_probe"]
    activated_probe = configuration["activated_probe"]
    reference = configuration["reference"]
    if not isinstance(flat_probe, GaugeConstraintJet2) or not isinstance(
        activated_probe, GaugeConstraintJet2
    ):
        raise ValueError("PROP1 probes must be GaugeConstraintJet2 values")
    if not isinstance(reference, ReferenceConnection):
        raise ValueError("PROP1 reference must be a ReferenceConnection")
    tilde_factor = _q("tilde_normal_factor", configuration["tilde_normal_factor"])
    hat_factor = _q("hat_normal_factor", configuration["hat_normal_factor"])
    flat_radius = _q("flat_coordinate_radius", configuration["flat_coordinate_radius"])
    activated_radius = _q(
        "activated_coordinate_radius", configuration["activated_coordinate_radius"]
    )
    if not 1 < tilde_factor < hat_factor:
        raise ValueError("PROP1 requires 1 < tilde factor < hat factor")
    if min(flat_radius, activated_radius) <= reference.radial_domain_minimum:
        raise ValueError("PROP1 fixtures must lie inside the reference annulus")

    flat_fixture = configuration["flat_fixture"]
    activated_fixture = configuration["activated_fixture"]
    if not isinstance(flat_fixture, Mapping) or not isinstance(
        activated_fixture, Mapping
    ):
        raise ValueError("PROP1 fixtures must be mappings")
    flat_state = state_from_generalized_adm_pg_fixture(flat_fixture)
    activated_state = state_from_generalized_adm_pg_fixture(activated_fixture)
    implicit = modified_harmonic_implicit_certificate(
        {
            "fixture": flat_fixture,
            "reference": reference,
            "coordinate_radius": flat_radius,
            "tilde_normal_factor": tilde_factor,
            "hat_normal_factor": hat_factor,
        }
    )
    zero_probe = GaugeConstraintJet2(VectorJet2(0), VectorJet2(0))
    flat_zero = reference_gauge_propagation_operator(
        flat_state,
        zero_probe,
        reference=reference,
        coordinate_radius=flat_radius,
        hat_normal_factor=hat_factor,
    )
    flat_nonzero = reference_gauge_propagation_operator(
        flat_state,
        flat_probe,
        reference=reference,
        coordinate_radius=flat_radius,
        hat_normal_factor=hat_factor,
    )
    activated = reference_gauge_propagation_operator(
        activated_state,
        activated_probe,
        reference=reference,
        coordinate_radius=activated_radius,
        hat_normal_factor=hat_factor,
    )
    activated_zero = reference_gauge_propagation_operator(
        activated_state,
        zero_probe,
        reference=reference,
        coordinate_radius=activated_radius,
        hat_normal_factor=hat_factor,
    )
    principal = gauge_propagation_principal_blocks(
        activated_state,
        reference=reference,
        coordinate_radius=activated_radius,
        hat_normal_factor=hat_factor,
    )
    companion_probe = _independent_companion_probe(activated_probe)
    companion = reference_gauge_propagation_operator(
        activated_state,
        companion_probe,
        reference=reference,
        coordinate_radius=activated_radius,
        hat_normal_factor=hat_factor,
    )
    combined = reference_gauge_propagation_operator(
        activated_state,
        _combine_gauge_jets(activated_probe, companion_probe),
        reference=reference,
        coordinate_radius=activated_radius,
        hat_normal_factor=hat_factor,
    )
    superposition_exact = combined["coordinate_divergence"] == tuple(
        a + b
        for a, b in zip(
            activated["coordinate_divergence"],
            companion["coordinate_divergence"],
        )
    )
    homogeneity_factor = Q(-7, 5)
    scaled = reference_gauge_propagation_operator(
        activated_state,
        _scale_gauge_jet(activated_probe, homogeneity_factor),
        reference=reference,
        coordinate_radius=activated_radius,
        hat_normal_factor=hat_factor,
    )
    homogeneity_exact = scaled["coordinate_divergence"] == tuple(
        homogeneity_factor * value
        for value in activated["coordinate_divergence"]
    )
    coefficients = principal["C_value_first_second_partial_blocks"]
    verified = {
        "implicit_flat_fgcqr_root_predecessor_exact": implicit[
            "all_declared_exact_checks_pass"
        ]
        and implicit[
            "local_smooth_implicit_coordinate_time_acceleration_branch_proven"
        ],
        "flat_zero_gauge_operator_exact": flat_zero["coordinate_divergence"]
        == (Q(0),) * N,
        "flat_zero_two_routes_agree_exactly": flat_zero[
            "two_routes_agree_exactly"
        ],
        "flat_nonzero_probe_exercises_operator": _tensor_nonzero(
            flat_nonzero["coordinate_divergence"]
        ),
        "flat_nonzero_two_routes_agree_exactly": flat_nonzero[
            "two_routes_agree_exactly"
        ],
        "activated_two_routes_agree_exactly": activated[
            "two_routes_agree_exactly"
        ],
        "activated_gradient_F_contribution_nonzero": _tensor_nonzero(
            activated["expanded_gradient_F_term"]
        ),
        "activated_projector_derivative_contribution_nonzero": _tensor_nonzero(
            activated["expanded_projector_derivative_term"]
        ),
        "activated_connection_contribution_nonzero": _tensor_nonzero(
            activated["expanded_connection_wave_contribution"]
        ),
        "activated_curvature_contribution_nonzero": _tensor_nonzero(
            activated["expanded_curvature_commutator_term"]
        ),
        "activated_first_gauge_jet_coefficients_nonzero": _tensor_nonzero(
            (coefficients["dt"], coefficients["dr"])
        ),
        "activated_second_gauge_jet_coefficients_nonzero": _tensor_nonzero(
            (coefficients["dtt"], coefficients["dtr"], coefficients["drr"])
        ),
        "exact_linear_homogeneity": homogeneity_exact,
        "exact_linear_superposition": superposition_exact,
        "spherical_principal_columns_match_mhg1_exactly": principal[
            "spherical_t_r_columns_match_mhg1_exactly"
        ],
        "angular_covariant_effects_retained_without_independent_angular_C": principal[
            "angular_independent_C_columns_not_represented"
        ],
        "operator_input_explicitly_independent_not_metric_derived": activated[
            "gauge_input_is_independent"
        ]
        and not activated["direct_metric_derived_residual_divergence_evaluated"],
        "reference_connection_is_not_an_independent_operator_source": activated_zero[
            "coordinate_divergence"
        ]
        == (Q(0),) * N
        and activated_zero["two_routes_agree_exactly"],
    }
    return {
        "artifact_id": "FGC-1-HYP1-MHG4-PROP1",
        "classification": "exact_lower_order_reference_gauge_operator_conditional_on_noether_and_scalar_equations_not_direct_metric_residual_divergence_or_constraint_propagation",
        "formulation": {
            "physical_equations": "unredefined_ACT1_VAR1",
            "gauge_equations": "REF1_full_reference_connection_modified_harmonic",
            "operator": "nabla_mu[F*hat_P_alpha^(beta_mu_nu)*nabla_beta_C^alpha]",
            "operator_input": "independent_spherically_symmetric_contravariant_gauge_vector_two_jet",
            "noether_identity": "nabla_mu_E^mu_nu+E_phi*nabla^nu_phi+E_chi*nabla^nu_chi=0",
            "off_scalar_shell_source": "E_phi*nabla^nu_phi+E_chi*nabla^nu_chi",
            "homogeneous_only_on_both_scalar_equations": True,
            "direct_full_metric_residual_divergence_evaluated": False,
            "reference_connection_enters_through_metric_derived_C_not_as_independent_forcing": True,
            "gauge_vector_order": SPHERICAL_GAUGE_VECTOR_ORDER,
            "gauge_jet_order": ("value", "dt", "dr", "dtt", "dtr", "drr"),
            "output_order": ("nu=t", "nu=r", "nu=theta", "nu=phi"),
            "flat_coordinate_radius": flat_radius,
            "activated_coordinate_radius": activated_radius,
            "tilde_normal_factor": tilde_factor,
            "hat_normal_factor": hat_factor,
        },
        "flat_root_control": {
            "fixture_id": str(flat_fixture.get("fixture_id", "unspecified")),
            "zero_probe": _operator_summary(flat_zero),
            "nonzero_independent_probe": _operator_summary(flat_nonzero),
        },
        "activated_operator_probe": {
            "fixture_id": str(activated_fixture.get("fixture_id", "unspecified")),
            "probe_is_metric_derived_C": False,
            "operator": _operator_summary(activated),
            "coefficient_blocks": coefficients,
            "zero_independent_probe_control": _operator_summary(activated_zero),
            "homogeneity_factor": homogeneity_factor,
        },
        "principal_regression": principal,
        "verified_exact_checks": verified,
        "all_declared_exact_checks_pass": all(verified.values()),
        "conditional_lower_order_reference_gauge_operator_derived": all(
            verified.values()
        ),
        "nonclaims": required_prop1_nonclaims(),
    }
