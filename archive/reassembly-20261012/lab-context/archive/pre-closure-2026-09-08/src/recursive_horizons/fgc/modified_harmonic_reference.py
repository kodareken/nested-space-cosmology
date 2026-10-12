"""Exact spherical reference-gauge modified-harmonic residuals for ACT1.

MHG1 certifies the frozen principal symbol.  This module supplies the full
reference-gauge completion of that principal residual on a positive-radius
spherical annulus.  It uses
the tensorial gauge vector

``C^a = -tilde_g^(bc) (Gamma^a_bc - bar_Gamma^a_bc)``

and adds ``F hat_P_c^(dab) nabla_d C^c`` to the unredefined VAR1 metric
equation.  All calculations are exact rational evaluations of local two-jets.
No implicit acceleration solve, constraint-propagation theorem, open-domain
bound, boundary system, or evolution is supplied here.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction
from typing import Any, Mapping, Sequence

from .exact_interval import Interval
from .exact_linear_algebra import Polynomial, ZERO, poly
from .interval_tangent import IntervalFirstTangent
from .modified_harmonic import (
    _trace_reversal_projector,
    auxiliary_inverse_metric,
    modified_harmonic_gauge_symbol,
)
from .reference_connection import (
    ReferenceConnection,
    flat_spherical_annulus_connection,
)
from .spherical_reduction import (
    BASE_FIELD_ORDER,
    INDEPENDENT_EQUATION_ORDER,
    SECOND_DERIVATIVE_ORDER,
    SphericalState,
    _direct_metric,
    _inv,
    residuals,
    state_from_generalized_adm_pg_fixture,
)


Q = Fraction
N = 4

MHG2_METRIC_EQUATION_ORDER = (
    "metric_tt_mhg",
    "metric_tr_mhg",
    "metric_rr_mhg",
    "metric_theta_theta_mhg",
)
MHG2_FULL_EQUATION_ORDER = MHG2_METRIC_EQUATION_ORDER + (
    "scalar_phi",
    "scalar_chi",
)


def _freeze(value):
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value


def _zero2() -> list[list[Fraction]]:
    return [[Q(0) for _ in range(N)] for _ in range(N)]


def _zero3() -> list[list[list[Fraction]]]:
    return [
        [[Q(0) for _ in range(N)] for _ in range(N)]
        for _ in range(N)
    ]


def _zero4() -> list[list[list[list[Fraction]]]]:
    return [
        [
            [[Q(0) for _ in range(N)] for _ in range(N)]
            for _ in range(N)
        ]
        for _ in range(N)
    ]


@dataclass(frozen=True, slots=True)
class PhysicalConnectionData:
    metric: tuple[tuple[Fraction, ...], ...]
    inverse_metric: tuple[tuple[Fraction, ...], ...]
    metric_derivative: tuple[tuple[tuple[Fraction, ...], ...], ...]
    inverse_metric_derivative: tuple[tuple[tuple[Fraction, ...], ...], ...]
    christoffel: tuple[tuple[tuple[Fraction, ...], ...], ...]
    christoffel_derivative: tuple[
        tuple[tuple[tuple[Fraction, ...], ...], ...], ...
    ]


@dataclass(frozen=True, slots=True)
class ModifiedHarmonicGaugeData:
    reference_id: str
    coordinate_radius: Fraction
    physical: PhysicalConnectionData
    gamma_difference: tuple[tuple[tuple[Fraction, ...], ...], ...]
    gamma_difference_derivative: tuple[
        tuple[tuple[tuple[Fraction, ...], ...], ...], ...
    ]
    tilde_inverse_metric: tuple[tuple[Fraction, ...], ...]
    tilde_inverse_metric_derivative: tuple[
        tuple[tuple[Fraction, ...], ...], ...
    ]
    constraint_up: tuple[Fraction, ...]
    constraint_down: tuple[Fraction, ...]
    partial_constraint_up: tuple[tuple[Fraction, ...], ...]
    covariant_constraint_derivative: tuple[tuple[Fraction, ...], ...]


def physical_connection_data(state: SphericalState) -> PhysicalConnectionData:
    """Return ``g,dg,g^-1,dg^-1,Gamma,dGamma`` from RED1's exact metric jet."""

    metric, metric_derivative, metric_second_derivative = _direct_metric(state)
    inverse = _inv(metric)
    inverse_derivative = [_zero2() for _ in range(N)]
    for direction in range(N):
        for a in range(N):
            for b in range(N):
                inverse_derivative[direction][a][b] = -sum(
                    inverse[a][m]
                    * metric_derivative[direction][m][n]
                    * inverse[n][b]
                    for m in range(N)
                    for n in range(N)
                )

    gamma = _zero3()
    for a in range(N):
        for b in range(N):
            for c in range(N):
                gamma[a][b][c] = sum(
                    inverse[a][d]
                    * (
                        metric_derivative[b][d][c]
                        + metric_derivative[c][d][b]
                        - metric_derivative[d][b][c]
                    )
                    / 2
                    for d in range(N)
                )

    gamma_derivative = _zero4()
    for direction in range(N):
        for a in range(N):
            for b in range(N):
                for c in range(N):
                    gamma_derivative[direction][a][b][c] = sum(
                        inverse_derivative[direction][a][d]
                        * (
                            metric_derivative[b][d][c]
                            + metric_derivative[c][d][b]
                            - metric_derivative[d][b][c]
                        )
                        / 2
                        + inverse[a][d]
                        * (
                            metric_second_derivative[direction][b][d][c]
                            + metric_second_derivative[direction][c][d][b]
                            - metric_second_derivative[direction][d][b][c]
                        )
                        / 2
                        for d in range(N)
                    )

    return PhysicalConnectionData(
        metric=_freeze(metric),
        inverse_metric=_freeze(inverse),
        metric_derivative=_freeze(metric_derivative),
        inverse_metric_derivative=_freeze(inverse_derivative),
        christoffel=_freeze(gamma),
        christoffel_derivative=_freeze(gamma_derivative),
    )


def _auxiliary_inverse_with_derivative(
    physical: PhysicalConnectionData,
    normal_factor: int | Fraction,
) -> tuple[
    tuple[tuple[Fraction, ...], ...],
    tuple[tuple[tuple[Fraction, ...], ...], ...],
]:
    factor = Q(normal_factor)
    inverse = physical.inverse_metric
    inverse_derivative = physical.inverse_metric_derivative
    auxiliary = auxiliary_inverse_metric(inverse, factor)
    normal_outer = [
        [
            -inverse[a][0] * inverse[b][0] / inverse[0][0]
            for b in range(N)
        ]
        for a in range(N)
    ]
    normal_outer_derivative = [_zero2() for _ in range(N)]
    for direction in range(N):
        denominator = inverse[0][0]
        denominator_derivative = inverse_derivative[direction][0][0]
        for a in range(N):
            for b in range(N):
                numerator = inverse[a][0] * inverse[b][0]
                numerator_derivative = (
                    inverse_derivative[direction][a][0] * inverse[b][0]
                    + inverse[a][0] * inverse_derivative[direction][b][0]
                )
                normal_outer_derivative[direction][a][b] = (
                    -numerator_derivative / denominator
                    + numerator * denominator_derivative / denominator**2
                )
    auxiliary_derivative = [
        [
            [
                inverse_derivative[direction][a][b]
                - (factor - 1) * normal_outer_derivative[direction][a][b]
                for b in range(N)
            ]
            for a in range(N)
        ]
        for direction in range(N)
    ]
    expected = tuple(
        tuple(
            inverse[a][b] - (factor - 1) * normal_outer[a][b]
            for b in range(N)
        )
        for a in range(N)
    )
    if auxiliary != expected:
        raise ValueError("differentiated auxiliary metric disagrees at zeroth order")
    return auxiliary, _freeze(auxiliary_derivative)


def modified_harmonic_gauge_constraint(
    state: SphericalState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: int | Fraction,
    tilde_normal_factor: int | Fraction,
) -> ModifiedHarmonicGaugeData:
    """Evaluate the tensorial contravariant MHG gauge vector and ``nabla C``."""

    physical = physical_connection_data(state)
    reference_data = flat_spherical_annulus_connection(
        reference, coordinate_radius=coordinate_radius
    )
    gamma_difference = _zero3()
    gamma_difference_derivative = _zero4()
    for a in range(N):
        for b in range(N):
            for c in range(N):
                gamma_difference[a][b][c] = (
                    physical.christoffel[a][b][c]
                    - reference_data.christoffel[a][b][c]
                )
                for direction in range(N):
                    gamma_difference_derivative[direction][a][b][c] = (
                        physical.christoffel_derivative[direction][a][b][c]
                        - reference_data.christoffel_derivative[direction][a][b][c]
                    )

    tilde, tilde_derivative = _auxiliary_inverse_with_derivative(
        physical, tilde_normal_factor
    )
    constraint = tuple(
        -sum(
            tilde[rho][sigma] * gamma_difference[a][rho][sigma]
            for rho in range(N)
            for sigma in range(N)
        )
        for a in range(N)
    )
    partial_constraint = tuple(
        tuple(
            -sum(
                tilde_derivative[direction][rho][sigma]
                * gamma_difference[a][rho][sigma]
                + tilde[rho][sigma]
                * gamma_difference_derivative[direction][a][rho][sigma]
                for rho in range(N)
                for sigma in range(N)
            )
            for a in range(N)
        )
        for direction in range(N)
    )
    covariant_derivative = tuple(
        tuple(
            partial_constraint[direction][a]
            + sum(
                physical.christoffel[a][direction][c] * constraint[c]
                for c in range(N)
            )
            for a in range(N)
        )
        for direction in range(N)
    )
    constraint_down = tuple(
        sum(physical.metric[a][b] * constraint[b] for b in range(N))
        for a in range(N)
    )
    return ModifiedHarmonicGaugeData(
        reference_id=reference.reference_id,
        coordinate_radius=Q(coordinate_radius),
        physical=physical,
        gamma_difference=_freeze(gamma_difference),
        gamma_difference_derivative=_freeze(gamma_difference_derivative),
        tilde_inverse_metric=tilde,
        tilde_inverse_metric_derivative=tilde_derivative,
        constraint_up=constraint,
        constraint_down=constraint_down,
        partial_constraint_up=partial_constraint,
        covariant_constraint_derivative=covariant_derivative,
    )


def modified_harmonic_extension_residual(
    state: SphericalState,
    *,
    gauge: ModifiedHarmonicGaugeData,
    hat_normal_factor: int | Fraction,
) -> dict[str, Any]:
    """Return the full MHG gauge extension with both index forms."""

    if gauge.physical != physical_connection_data(state):
        raise ValueError("gauge data and state geometry differ")
    hat = auxiliary_inverse_metric(
        gauge.physical.inverse_metric, Q(hat_normal_factor)
    )
    effective_planck = state.planck_mass**2 + state.beta * state.phi.value**2
    effective_planck_primal = (
        effective_planck.primal
        if isinstance(effective_planck, IntervalFirstTangent)
        else effective_planck
    )
    if isinstance(effective_planck_primal, Interval):
        if not effective_planck_primal.strictly_positive():
            raise ValueError("modified-harmonic gauge extension requires F>0")
    elif effective_planck_primal <= 0:
        raise ValueError("modified-harmonic gauge extension requires F>0")

    extension_up = _zero2()
    for mu in range(N):
        for nu in range(N):
            extension_up[mu][nu] = effective_planck * sum(
                _trace_reversal_projector(hat, alpha, beta, mu, nu)
                * gauge.covariant_constraint_derivative[beta][alpha]
                for alpha in range(N)
                for beta in range(N)
            )
    extension_down = _zero2()
    metric = gauge.physical.metric
    for a in range(N):
        for b in range(N):
            extension_down[a][b] = sum(
                metric[a][mu] * metric[b][nu] * extension_up[mu][nu]
                for mu in range(N)
                for nu in range(N)
            )
    return {
        "hat_inverse_metric": hat,
        "effective_planck_coefficient": effective_planck,
        "contravariant": _freeze(extension_up),
        "covariant": _freeze(extension_down),
        "contravariant_symmetric": all(
            extension_up[a][b] == extension_up[b][a]
            for a in range(N)
            for b in range(N)
        ),
        "covariant_symmetric": all(
            extension_down[a][b] == extension_down[b][a]
            for a in range(N)
            for b in range(N)
        ),
        "equatorial_angular_isotropy": extension_down[2][2]
        == extension_down[3][3],
    }


def _independent_metric_vector(metric: Sequence[Sequence[Fraction]]) -> tuple[Fraction, ...]:
    return (metric[0][0], metric[0][1], metric[1][1], metric[2][2])


def modified_harmonic_full_residuals(
    state: SphericalState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: int | Fraction,
    tilde_normal_factor: int | Fraction,
    hat_normal_factor: int | Fraction,
) -> dict[str, Any]:
    """Evaluate unredefined VAR1, the gauge extension, and their exact sum."""

    gauge = modified_harmonic_gauge_constraint(
        state,
        reference=reference,
        coordinate_radius=coordinate_radius,
        tilde_normal_factor=tilde_normal_factor,
    )
    extension = modified_harmonic_extension_residual(
        state, gauge=gauge, hat_normal_factor=hat_normal_factor
    )
    original = residuals(state)
    full_metric = tuple(
        tuple(
            original["metric"][a][b] + extension["covariant"][a][b]
            for b in range(N)
        )
        for a in range(N)
    )
    original_vector = _independent_metric_vector(original["metric"]) + (
        original["phi"],
        original["chi"],
    )
    extension_vector = _independent_metric_vector(extension["covariant"]) + (
        Q(0),
        Q(0),
    )
    full_vector = _independent_metric_vector(full_metric) + (
        original["phi"],
        original["chi"],
    )
    return {
        "equation_order": MHG2_FULL_EQUATION_ORDER,
        "gauge": gauge,
        "extension": extension,
        "unredefined_metric": original["metric"],
        "full_metric": full_metric,
        "unredefined_residual_vector": original_vector,
        "extension_residual_vector": extension_vector,
        "full_residual_vector": full_vector,
        "full_equals_unredefined_plus_extension": full_vector
        == tuple(a + b for a, b in zip(original_vector, extension_vector)),
        "scalar_equations_unmodified": full_vector[4:] == original_vector[4:],
    }


def reference_extension_principal_symbol(
    state: SphericalState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: int | Fraction,
    tilde_normal_factor: int | Fraction,
    hat_normal_factor: int | Fraction,
) -> tuple[tuple[Polynomial, ...], ...]:
    """Differentiate the full extension in every base-field second-jet slot."""

    columns: dict[tuple[str, str], tuple[Fraction, ...]] = {}
    for field in BASE_FIELD_ORDER:
        for derivative in SECOND_DERIVATIVE_ORDER:
            jet = getattr(state, field)
            plus = replace(
                state,
                **{field: replace(jet, **{derivative: getattr(jet, derivative) + 1})},
            )
            minus = replace(
                state,
                **{field: replace(jet, **{derivative: getattr(jet, derivative) - 1})},
            )
            plus_vector = modified_harmonic_full_residuals(
                plus,
                reference=reference,
                coordinate_radius=coordinate_radius,
                tilde_normal_factor=tilde_normal_factor,
                hat_normal_factor=hat_normal_factor,
            )["extension_residual_vector"]
            minus_vector = modified_harmonic_full_residuals(
                minus,
                reference=reference,
                coordinate_radius=coordinate_radius,
                tilde_normal_factor=tilde_normal_factor,
                hat_normal_factor=hat_normal_factor,
            )["extension_residual_vector"]
            columns[(field, derivative)] = tuple(
                (a - b) / 2 for a, b in zip(plus_vector, minus_vector)
            )

    output = [[ZERO for _ in BASE_FIELD_ORDER] for _ in INDEPENDENT_EQUATION_ORDER]
    for equation in range(len(INDEPENDENT_EQUATION_ORDER)):
        for field_index, field in enumerate(BASE_FIELD_ORDER):
            dtt = columns[(field, "dtt")][equation]
            dtr = columns[(field, "dtr")][equation]
            drr = columns[(field, "drr")][equation]
            output[equation][field_index] = poly((drr, -dtr, dtt))
    return tuple(tuple(row) for row in output)


def _tensor_zero(value: Any) -> bool:
    if isinstance(value, Fraction):
        return value == 0
    if isinstance(value, (tuple, list)):
        return all(_tensor_zero(item) for item in value)
    raise ValueError("exact tensor contains a non-rational value")


def modified_harmonic_reference_fixture_certificate(
    fixture: Mapping[str, Any],
    *,
    reference: ReferenceConnection,
    coordinate_radius: int | Fraction,
    tilde_normal_factor: int | Fraction,
    hat_normal_factor: int | Fraction,
) -> dict[str, Any]:
    """Build one exact REF1 reference-gauge fixture certificate."""

    state = state_from_generalized_adm_pg_fixture(fixture)
    full = modified_harmonic_full_residuals(
        state,
        reference=reference,
        coordinate_radius=coordinate_radius,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    derived_symbol = reference_extension_principal_symbol(
        state,
        reference=reference,
        coordinate_radius=coordinate_radius,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    expected_symbol = modified_harmonic_gauge_symbol(
        state,
        tilde_normal_factor=tilde_normal_factor,
        hat_normal_factor=hat_normal_factor,
    )
    gauge = full["gauge"]
    extension = full["extension"]
    gauge_surface = _tensor_zero(gauge.constraint_up) and _tensor_zero(
        gauge.covariant_constraint_derivative
    )
    return {
        "fixture_id": str(fixture.get("fixture_id", "unspecified")),
        "model_id": str(fixture.get("model_id", "unspecified")),
        "purpose": str(fixture.get("purpose", "unspecified")),
        "coordinate_radius": Q(coordinate_radius),
        "effective_planck_coefficient": extension["effective_planck_coefficient"],
        "constraint_up": gauge.constraint_up,
        "constraint_down": gauge.constraint_down,
        "partial_constraint_up": gauge.partial_constraint_up,
        "covariant_constraint_derivative": gauge.covariant_constraint_derivative,
        "extension_covariant": extension["covariant"],
        "unredefined_residual_vector": full["unredefined_residual_vector"],
        "extension_residual_vector": full["extension_residual_vector"],
        "full_residual_vector": full["full_residual_vector"],
        "gauge_surface_at_evaluation_jet": gauge_surface,
        "gauge_surface_extension_zero": (not gauge_surface)
        or _tensor_zero(extension["covariant"]),
        "gauge_surface_full_equals_unredefined": (not gauge_surface)
        or full["full_residual_vector"] == full["unredefined_residual_vector"],
        "full_equals_unredefined_plus_extension": full[
            "full_equals_unredefined_plus_extension"
        ],
        "scalar_equations_unmodified": full["scalar_equations_unmodified"],
        "extension_tensor_symmetric": extension["contravariant_symmetric"]
        and extension["covariant_symmetric"],
        "equatorial_angular_isotropy": extension[
            "equatorial_angular_isotropy"
        ],
        "reference_extension_principal_symbol": derived_symbol,
        "mhg1_principal_gauge_symbol": expected_symbol,
        "reference_extension_principal_regression_exact": derived_symbol
        == expected_symbol,
        "extension_nonzero": not _tensor_zero(extension["covariant"]),
    }


def required_ref1_nonclaims() -> dict[str, bool]:
    return {
        name: False
        for name in (
            "implicit_second_time_derivative_branch_proven",
            "complete_lower_order_first_order_sources_derived",
            "gauge_constraint_lower_order_propagation_proven",
            "reduction_constraint_propagation_proven",
            "uniform_open_domain_symmetrizer_proven",
            "open_retained_eft_domain_proven",
            "boundary_or_regular_center_system_derived",
            "evolution_authorized",
            "collapse_solution_derived",
            "metric_null_affine_defocusing_derived",
            "singularity_resolution_derived",
        )
    }


def modified_harmonic_reference_certificate(
    configuration: Mapping[str, Any],
) -> dict[str, Any]:
    """Build the aggregate exact REF1 reference-gauge certificate."""

    expected_keys = {
        "fixtures",
        "coordinate_radii",
        "reference",
        "tilde_normal_factor",
        "hat_normal_factor",
        "gauge_surface_fixture_ids",
        "off_gauge_fixture_ids",
    }
    if not isinstance(configuration, Mapping) or set(configuration) != expected_keys:
        raise ValueError("REF1 configuration adapter has unexpected keys")
    fixtures = configuration["fixtures"]
    coordinate_radii = configuration["coordinate_radii"]
    reference = configuration["reference"]
    if not isinstance(fixtures, (list, tuple)) or not fixtures:
        raise ValueError("REF1 requires a nonempty fixture list")
    if not isinstance(coordinate_radii, Mapping):
        raise ValueError("REF1 coordinate radii must be a mapping")
    if not isinstance(reference, ReferenceConnection):
        raise ValueError("REF1 reference must be a ReferenceConnection")
    fixture_ids = tuple(str(fixture.get("fixture_id", "")) for fixture in fixtures)
    if any(not fixture_id for fixture_id in fixture_ids):
        raise ValueError("REF1 fixtures require identifiers")
    if set(coordinate_radii) != set(fixture_ids):
        raise ValueError("REF1 coordinate-radius keys must match fixture identifiers")
    tilde_factor = Q(configuration["tilde_normal_factor"])
    hat_factor = Q(configuration["hat_normal_factor"])
    if not 1 < tilde_factor < hat_factor:
        raise ValueError("REF1 requires 1 < tilde factor < hat factor")
    gauge_surface_ids = tuple(configuration["gauge_surface_fixture_ids"])
    off_gauge_ids = tuple(configuration["off_gauge_fixture_ids"])
    if set(gauge_surface_ids) & set(off_gauge_ids):
        raise ValueError("REF1 gauge-surface and off-gauge fixtures must be disjoint")
    if set(gauge_surface_ids) | set(off_gauge_ids) != set(fixture_ids):
        raise ValueError("REF1 gauge classification must cover every fixture")

    records = tuple(
        modified_harmonic_reference_fixture_certificate(
            fixture,
            reference=reference,
            coordinate_radius=Q(coordinate_radii[fixture_id]),
            tilde_normal_factor=tilde_factor,
            hat_normal_factor=hat_factor,
        )
        for fixture, fixture_id in zip(fixtures, fixture_ids)
    )
    by_id = {record["fixture_id"]: record for record in records}
    schwarzschild = by_id.get("FGCQR_schwarzschild_exterior")
    verified = {
        "reference_annulus_strictly_positive": reference.radial_domain_minimum
        > 0
        and all(
            record["coordinate_radius"] > reference.radial_domain_minimum
            for record in records
        ),
        "positive_effective_planck_coefficients": all(
            record["effective_planck_coefficient"] > 0 for record in records
        ),
        "declared_gauge_surface_jets_exact": all(
            by_id[fixture_id]["gauge_surface_at_evaluation_jet"]
            and by_id[fixture_id]["gauge_surface_extension_zero"]
            and by_id[fixture_id]["gauge_surface_full_equals_unredefined"]
            for fixture_id in gauge_surface_ids
        ),
        "declared_off_gauge_jets_exercise_extension": all(
            not by_id[fixture_id]["gauge_surface_at_evaluation_jet"]
            and by_id[fixture_id]["extension_nonzero"]
            for fixture_id in off_gauge_ids
        ),
        "full_residual_decomposition_exact": all(
            record["full_equals_unredefined_plus_extension"]
            for record in records
        ),
        "scalar_equations_unmodified": all(
            record["scalar_equations_unmodified"] for record in records
        ),
        "extension_tensors_symmetric_and_spherical": all(
            record["extension_tensor_symmetric"]
            and record["equatorial_angular_isotropy"]
            for record in records
        ),
        "reference_extension_directional_symbol_matches_mhg1_gauge_block": all(
            record["reference_extension_principal_regression_exact"]
            for record in records
        ),
        "schwarzschild_unredefined_solution_control_retained": schwarzschild
        is not None
        and _tensor_zero(schwarzschild["unredefined_residual_vector"]),
        "contravariant_spherical_gauge_components_only": all(
            record["constraint_up"][2:] == (Q(0), Q(0))
            for record in records
        ),
    }
    return {
        "artifact_id": "FGC-1-HYP1-MHG2-REF1",
        "classification": "exact_reference_connection_modified_harmonic_residual_not_implicit_solve_propagation_domain_or_evolution",
        "formulation": {
            "physical_equations": "unredefined_ACT1_VAR1",
            "reference_connection": reference.reference_id,
            "reference_domain": "positive_radius_spherical_annulus_without_center",
            "reference_radial_domain_minimum": reference.radial_domain_minimum,
            "gauge_constraint": "C^mu=-tilde_g^rho_sigma*(Gamma^mu_rho_sigma-bar_Gamma^mu_rho_sigma)",
            "metric_extension": "E_mhg^mu_nu=E^mu_nu+F*hat_P_alpha^beta_mu_nu*nabla_beta_C^alpha",
            "tilde_normal_factor": tilde_factor,
            "hat_normal_factor": hat_factor,
            "scalar_equations": "unmodified",
        },
        "fixture_records": records,
        "verified_exact_checks": verified,
        "all_declared_exact_checks_pass": all(verified.values()),
        "nonclaims": required_ref1_nonclaims(),
    }
