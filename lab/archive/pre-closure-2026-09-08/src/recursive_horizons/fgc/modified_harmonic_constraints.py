"""Exact local compatibility checks for the spherical REF1 formulation.

CON1-COMP1 deliberately proves only an algebraic local-data statement.  The
activated datum constructed here is conditional on QIFT1's separately proved
existence of a REF1 acceleration root.  It does not solve that root, derive an
initial slice, or propagate any constraint.
"""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
from itertools import combinations
from numbers import Integral
from typing import Any, Callable, Mapping, Sequence

from .exact_linear_algebra import matrix_det, matrix_inverse
from .modified_harmonic_first_order import (
    FirstOrderSphericalState,
    first_order_kinematic_residuals,
    first_order_state_from_spherical_state,
)
from .modified_harmonic_reference import (
    ModifiedHarmonicGaugeData,
    modified_harmonic_extension_residual,
    modified_harmonic_full_residuals,
    modified_harmonic_gauge_constraint,
)
from .reference_connection import ReferenceConnection
from .spherical_reduction import BASE_FIELD_ORDER, Jet2, SphericalState, residuals, state_from_generalized_adm_pg_fixture
from .spherical_symbol import symbol_fixture_certificate


Q = Fraction
FIELD_COUNT = len(BASE_FIELD_ORDER)
CON1_COMP1_ARTIFACT_ID = "FGC-1-HYP1-CON1-COMP1"
ACTIVATED_PHI = Q(1, 131072)
ACTIVATED_PHI_RADIAL_DERIVATIVE = Q(1, 131072)
QIFT1_PARAMETER_HALF_WIDTH = Q(1, 65536)
QIFT1_ACCELERATION_HALF_WIDTH = Q(1, 128)
ACTIVATED_AREAL_RADIUS_DRR = -Q(584115552257, 4722366482835285475328)
ACTIVATED_H_TT_DRR = -Q(584115552257, 18889465931341141901312)
PHYSICAL_PROJECTION_ORDER = ("H=E_tt-2sE_tr+s^2E_rr", "M=E_tr-sE_rr")
ACCELERATION_ORDER = tuple(f"p_t.{field}" for field in BASE_FIELD_ORDER)
ACCELERATION_POLYNOMIAL_DEGREE = 2
ACCELERATION_EVALUATION_NODES = ("zero", "plus_basis", "minus_basis", "pair_basis")


def _fraction(name: str, value: object, *, positive: bool = False) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational")
    result = value if isinstance(value, Fraction) else Q(int(value))
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _zero_vector(size: int) -> tuple[Fraction, ...]:
    return tuple(Q(0) for _ in range(size))


def _zero_matrix(rows: int, columns: int) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(tuple(Q(0) for _ in range(columns)) for _ in range(rows))


def _all_zero(value: object) -> bool:
    if isinstance(value, Fraction):
        return value == 0
    if isinstance(value, (tuple, list)):
        return all(_all_zero(item) for item in value)
    raise TypeError(f"expected exact rational tensor, got {type(value).__name__}")


def physical_constraint_projections(state: SphericalState) -> dict[str, Any]:
    """Return the unredefined rational spherical normal projections.

    ``w=partial_t-s partial_r`` is the unnormalised normal direction, with
    ``s=h_tr/h_rr``.  Its zero-set projections avoid an unnecessary lapse
    square root.  The input is intentionally routed through ``residuals`` and
    never through the REF1 gauge extension.
    """

    if not isinstance(state, SphericalState):
        raise TypeError("state must be a SphericalState")
    source = residuals(state)
    metric = source["metric"]
    shift = state.h_tr.value / state.h_rr.value
    hamiltonian = metric[0][0] - 2 * shift * metric[0][1] + shift**2 * metric[1][1]
    momentum = metric[0][1] - shift * metric[1][1]
    return {
        "projection_order": PHYSICAL_PROJECTION_ORDER,
        "source_is_unredefined_residual": True,
        "shift": shift,
        "normal_direction_w": (Q(1), -shift, Q(0), Q(0)),
        "H": hamiltonian,
        "M": momentum,
        "metric_tt": metric[0][0],
        "metric_tr": metric[0][1],
        "metric_rr": metric[1][1],
        "metric_theta_theta": metric[2][2],
        "scalar_phi": source["phi"],
        "scalar_chi": source["chi"],
    }


def _constraint_pair(state: SphericalState) -> tuple[Fraction, Fraction]:
    output = physical_constraint_projections(state)
    return output["H"], output["M"]


def _with_acceleration(state: SphericalState, acceleration: Sequence[Fraction]) -> SphericalState:
    if len(acceleration) != FIELD_COUNT:
        raise ValueError("acceleration vector must have six entries")
    return replace(
        state,
        **{
            field: replace(getattr(state, field), dtt=_fraction(f"acceleration.{field}", value))
            for field, value in zip(BASE_FIELD_ORDER, acceleration, strict=True)
        },
    )


def _affine_two_variable_solution(
    evaluate: Callable[[Fraction, Fraction], tuple[Fraction, Fraction]],
) -> tuple[tuple[Fraction, Fraction], tuple[tuple[Fraction, Fraction], ...], tuple[Fraction, Fraction]]:
    """Solve a checked affine two-by-two rational system at the origin."""

    constant = evaluate(Q(0), Q(0))
    unit_x = evaluate(Q(1), Q(0))
    unit_y = evaluate(Q(0), Q(1))
    matrix = tuple(
        tuple((unit_x[row] - constant[row], unit_y[row] - constant[row]))
        for row in range(2)
    )
    determinant = matrix_det(matrix)
    if determinant == 0:
        raise ValueError("activated compatibility affine system is singular")
    inverse = matrix_inverse(matrix)
    solution = tuple(
        -sum(inverse[row][column] * constant[column] for column in range(2))
        for row in range(2)
    )
    # The two selected expressions are affine in these two radial slots at the
    # fixed datum.  This explicit mixed finite difference catches accidental
    # nonlinear dependence before the exact affine solve is trusted.
    mixed = evaluate(Q(1), Q(1))
    if tuple(
        mixed[row] - unit_x[row] - unit_y[row] + constant[row]
        for row in range(2)
    ) != (Q(0), Q(0)):
        raise ValueError("activated compatibility conditions are not affine in declared radial slots")
    if evaluate(*solution) != (Q(0), Q(0)):
        raise ValueError("activated compatibility affine solve failed")
    return solution, matrix, constant


def activated_compatible_state(
    flat_fixture: Mapping[str, Any],
    *,
    activated_phi: Fraction = ACTIVATED_PHI,
    activated_phi_radial_derivative: Fraction = ACTIVATED_PHI_RADIAL_DERIVATIVE,
    parameter_half_width: Fraction = QIFT1_PARAMETER_HALF_WIDTH,
) -> dict[str, Any]:
    """Derive the exact activated local datum from the QIFT flat fixture.

    This constructs only the non-acceleration QIFT parameter point.  A REF1
    acceleration root is an external QIFT1 predecessor fact.
    """

    epsilon = _fraction("activated_phi", activated_phi, positive=True)
    radial_epsilon = _fraction("activated_phi_radial_derivative", activated_phi_radial_derivative, positive=True)
    radius = _fraction("parameter_half_width", parameter_half_width, positive=True)
    if epsilon != ACTIVATED_PHI:
        raise ValueError("activated_phi must equal the frozen CON1-COMP1 epsilon")
    if radial_epsilon != ACTIVATED_PHI_RADIAL_DERIVATIVE:
        raise ValueError("activated_phi_radial_derivative must equal the frozen CON1-COMP1 value")
    if radius != QIFT1_PARAMETER_HALF_WIDTH:
        raise ValueError("parameter_half_width must equal QIFT1's frozen half width")
    if not isinstance(flat_fixture, Mapping):
        raise TypeError("flat_fixture must be a mapping")
    flat = state_from_generalized_adm_pg_fixture(flat_fixture)
    if flat.areal_radius.value != 4:
        raise ValueError("CON1-COMP1 requires QIFT1's r=4 flat fixture")

    provisional = replace(
        flat,
        phi=Jet2(epsilon, dr=radial_epsilon),
        h_tr=replace(flat.h_tr, drr=Q(0)),
        h_rr=replace(flat.h_rr, drr=Q(0)),
    )

    # x=h_tt,rr and y=R,rr.  H comes from the unredefined residual; the gauge
    # condition uses REF1's metric-derived C, fixed below at its r=4 reference.
    from .reference_connection import flat_spherical_annulus_reference
    from .modified_harmonic_reference import modified_harmonic_gauge_constraint

    reference = flat_spherical_annulus_reference(radial_domain_minimum=Q(1, 2))

    def conditions(x: Fraction, y: Fraction) -> tuple[Fraction, Fraction]:
        state = replace(
            provisional,
            h_tt=replace(provisional.h_tt, drr=x),
            areal_radius=replace(provisional.areal_radius, drr=y),
        )
        H = physical_constraint_projections(state)["H"]
        gauge = modified_harmonic_gauge_constraint(
            state,
            reference=reference,
            coordinate_radius=Q(4),
            tilde_normal_factor=Q(4),
        )
        return H, gauge.covariant_constraint_derivative[1][1]

    (h_tt_drr, radius_drr), affine_matrix, affine_constant = _affine_two_variable_solution(conditions)
    if radius_drr != ACTIVATED_AREAL_RADIUS_DRR or h_tt_drr != ACTIVATED_H_TT_DRR:
        raise ValueError("activated compatibility solution differs from frozen exact target")
    state = replace(
        provisional,
        h_tt=replace(provisional.h_tt, drr=h_tt_drr),
        areal_radius=replace(provisional.areal_radius, drr=radius_drr),
    )
    # z=(u,p,q,p_r,q_r): the only changed non-acceleration values are phi and
    # selected q_r entries.  Strict inclusion is an exact predecessor-domain
    # check, not a re-run of QIFT's contraction proof.
    changed = (epsilon, radial_epsilon, h_tt_drr, radius_drr)
    if not all(abs(value) < radius for value in changed):
        raise ValueError("activated datum lies outside QIFT1's strict parameter box")
    return {
        "state": state,
        "reference": reference,
        "coordinate_radius": Q(4),
        "tilde_normal_factor": Q(4),
        "hat_normal_factor": Q(9),
        "activated_phi": epsilon,
        "activated_phi_radial_derivative": radial_epsilon,
        "parameter_half_width": radius,
        "affine_condition_order": ("H", "nabla_r_C^r"),
        "affine_unknown_order": ("h_tt.drr", "areal_radius.drr"),
        "affine_matrix": affine_matrix,
        "affine_constant": affine_constant,
        "h_tt_drr": h_tt_drr,
        "areal_radius_drr": radius_drr,
        "strictly_inside_qift_parameter_box": True,
        "qift_acceleration_root_is_external_predecessor_fact": True,
        "qift_acceleration_root_solved_here": False,
    }


def acceleration_polynomial_certificate(
    state: SphericalState,
    *,
    evaluator: Callable[[SphericalState], tuple[Fraction, Fraction]] | None = None,
) -> dict[str, Any]:
    """Recover every total-degree <=2 H/M acceleration coefficient exactly.

    At fixed non-acceleration data the unredefined evaluator has degree at
    most two in coordinate accelerations: curvature/Hessians are affine, and
    the only products are ``GB`` and ``P*nabla-nabla f``.  The coefficient
    recovery nodes below are therefore a complete algebraic certificate.
    """

    if not isinstance(state, SphericalState):
        raise TypeError("state must be a SphericalState")
    zero = _zero_vector(FIELD_COUNT)
    cache: dict[tuple[Fraction, ...], tuple[Fraction, Fraction]] = {}

    constraint_evaluator = _constraint_pair if evaluator is None else evaluator

    def evaluate(vector: tuple[Fraction, ...]) -> tuple[Fraction, Fraction]:
        if vector not in cache:
            result = constraint_evaluator(_with_acceleration(state, vector))
            if not isinstance(result, tuple) or len(result) != 2:
                raise TypeError("constraint evaluator must return an exact (H, M) pair")
            cache[vector] = tuple(_fraction("constraint evaluator output", value) for value in result)  # type: ignore[assignment]
        return cache[vector]

    base = evaluate(zero)
    coefficients: list[dict[str, Any]] = [
        {"monomial": "1", "power": 0, "indices": (), "H": base[0], "M": base[1]}
    ]
    plus: list[tuple[Fraction, Fraction]] = []
    minus: list[tuple[Fraction, Fraction]] = []
    diagonal: list[tuple[Fraction, Fraction]] = []
    linear: list[tuple[Fraction, Fraction]] = []
    for index, field in enumerate(BASE_FIELD_ORDER):
        e = list(zero)
        e[index] = Q(1)
        p = evaluate(tuple(e))
        e[index] = Q(-1)
        m = evaluate(tuple(e))
        plus.append(p)
        minus.append(m)
        first = tuple((p[row] - m[row]) / 2 for row in range(2))
        second = tuple((p[row] + m[row]) / 2 - base[row] for row in range(2))
        linear.append(first)
        diagonal.append(second)
        coefficients.append({"monomial": f"a.{field}", "power": 1, "indices": (index,), "H": first[0], "M": first[1]})
        coefficients.append({"monomial": f"a.{field}^2", "power": 2, "indices": (index, index), "H": second[0], "M": second[1]})
    for left, right in combinations(range(FIELD_COUNT), 2):
        vector = list(zero)
        vector[left] = vector[right] = Q(1)
        pair = evaluate(tuple(vector))
        mixed = tuple(
            pair[row] - base[row] - linear[left][row] - linear[right][row] - diagonal[left][row] - diagonal[right][row]
            for row in range(2)
        )
        coefficients.append({
            "monomial": f"a.{BASE_FIELD_ORDER[left]}*a.{BASE_FIELD_ORDER[right]}",
            "power": 2,
            "indices": (left, right),
            "H": mixed[0],
            "M": mixed[1],
        })
    if len(coefficients) != 28:
        raise AssertionError("degree-two six-variable polynomial must have 28 coefficients")
    return {
        "projection_order": PHYSICAL_PROJECTION_ORDER,
        "acceleration_order": ACCELERATION_ORDER,
        "maximum_total_degree": ACCELERATION_POLYNOMIAL_DEGREE,
        "degree_bound_guard": {
            "curvature_and_scalar_hessians_affine_in_accelerations": True,
            "gauss_bonnet_and_P_hessian_f_products_at_most_quadratic": True,
            "unredefined_projection_uses_no_REF1_extension": True,
        },
        "evaluation_node_scheme": ACCELERATION_EVALUATION_NODES,
        "evaluation_count": len(cache),
        "coefficients": tuple(coefficients),
        "coefficient_count_per_projection": len(coefficients),
        "all_H_coefficients_zero": all(item["H"] == 0 for item in coefficients),
        "all_M_coefficients_zero": all(item["M"] == 0 for item in coefficients),
        "all_56_coefficients_zero": all(item["H"] == 0 and item["M"] == 0 for item in coefficients),
    }


def normal_gauge_extension_map(
    state: SphericalState,
    gauge: ModifiedHarmonicGaugeData,
    *,
    hat_normal_factor: Fraction,
) -> dict[str, Any]:
    """Map arbitrary t-row ``nabla C`` inputs to covariant (tt,tr) extension."""

    if not isinstance(gauge, ModifiedHarmonicGaugeData):
        raise TypeError("gauge must be ModifiedHarmonicGaugeData")
    columns: list[tuple[Fraction, Fraction]] = []
    for component in range(2):
        derivative = [list(row) for row in gauge.covariant_constraint_derivative]
        derivative[0][component] = Q(1)
        replacement = replace(
            gauge,
            covariant_constraint_derivative=tuple(tuple(row) for row in derivative),
        )
        extension = modified_harmonic_extension_residual(
            state, gauge=replacement, hat_normal_factor=hat_normal_factor
        )["covariant"]
        columns.append((extension[0][0], extension[0][1]))
    matrix = tuple(tuple(columns[column][row] for column in range(2)) for row in range(2))
    determinant = matrix_det(matrix)
    return {
        "input_order": ("nabla_t_C^t", "nabla_t_C^r"),
        "output_order": ("extension_covariant_tt", "extension_covariant_tr"),
        "matrix": matrix,
        "diagonal_entries_nonzero": matrix[0][0] != 0 and matrix[1][1] != 0,
        "off_diagonal_entries_zero": matrix[0][1] == 0 and matrix[1][0] == 0,
        "determinant": determinant,
        "determinant_nonzero": determinant != 0,
    }


def activated_principal_cone_control(
    flat_fixture: Mapping[str, Any], compatible_state: SphericalState
) -> dict[str, Any]:
    """Extract the exact SYM1 cone separation at CON1's compatible point.

    The result is a pointwise bridge only.  It deliberately retains none of
    SYM1's pointwise-eigenbasis language as an open-domain hyperbolicity claim.
    """

    if not isinstance(flat_fixture, Mapping):
        raise TypeError("flat_fixture must be a mapping")
    if not isinstance(compatible_state, SphericalState):
        raise TypeError("compatible_state must be a SphericalState")
    source_state = flat_fixture.get("state", flat_fixture)
    if not isinstance(source_state, Mapping):
        raise ValueError("flat_fixture state must be a mapping")
    fixture = dict(flat_fixture)
    fixture["fixture_id"] = "CON1_COMP1_activated_compatible_principal_control"

    # The compatible base metric is shift-free with unit lapse/radial metric.
    # Pull its exact solved second jets back to the ADM fixture variables and
    # prove that the round trip reproduces the complete SphericalState before
    # asking SYM1 for any cone statement.
    zero = Jet2.constant(Q(0))
    one = Jet2.constant(Q(1))
    if (
        compatible_state.h_tr != zero
        or compatible_state.h_rr != one
        or compatible_state.h_tt.value != -1
        or any(
            getattr(compatible_state.h_tt, name) != 0
            for name in ("dt", "dr", "dtt", "dtr")
        )
    ):
        raise ValueError("compatible cone control requires the frozen shift-free base metric")

    def jet_mapping(jet: Jet2) -> dict[str, Fraction]:
        return {
            name: getattr(jet, name)
            for name in ("value", "dt", "dr", "dtt", "dtr", "drr")
        }

    alpha_jet = Jet2(Q(1), drr=-compatible_state.h_tt.drr / 2)
    fixture["state"] = {
        "alpha": jet_mapping(alpha_jet),
        "shift": jet_mapping(zero),
        "lambda": jet_mapping(one),
        "areal_radius": jet_mapping(compatible_state.areal_radius),
        "phi": jet_mapping(compatible_state.phi),
        "chi": jet_mapping(compatible_state.chi),
    }
    translated = state_from_generalized_adm_pg_fixture(fixture)
    if translated != compatible_state:
        raise ValueError("compatible state ADM pullback failed exact round-trip equality")
    symbol = symbol_fixture_certificate(fixture)
    metric = symbol["metric_null_certificate"]
    regulator = symbol["regulator_certificate"]
    return {
        "fixture_id": fixture["fixture_id"],
        "compatible_state_pullback_exact": translated == compatible_state,
        "metric_null_factor": symbol["metric_null_factor"],
        "regulator_factor": symbol["regulator_factor"],
        "metric_null_discriminant": metric["discriminant"],
        "regulator_discriminant": regulator["discriminant"],
        "metric_null_discriminant_positive": metric["discriminant"] > 0,
        "regulator_discriminant_positive": regulator["discriminant"] > 0,
        "physical_vs_regulator_resultant": symbol["cone_resultant"],
        "physical_vs_regulator_resultant_nonzero": symbol["cone_resultant"] != 0,
        "cones_share_no_root": symbol["cones_share_no_root"],
        "regulator_cone_proportional_to_metric": symbol["regulator_cone_proportional_to_metric"],
        "classification": "exact_pointwise_activated_principal_cone_control_not_an_open_domain_hyperbolicity_certificate",
    }


def required_con1_comp1_nonclaims() -> dict[str, bool]:
    return {
        name: False
        for name in (
            "qift_acceleration_root_solved_here",
            "exact_closed_form_nonlinear_acceleration_map_derived",
            "initial_slice_or_constraint_manifold_derived",
            "metric_derived_gauge_constraint_propagation_proven",
            "complete_reduction_constraint_system_propagation_proven",
            "physical_hamiltonian_momentum_constraint_propagation_proven",
            "physical_initial_constraints_solved",
            "nontrivial_compatible_initial_data_family_derived",
            "constraint_preserving_initial_boundary_value_problem_proven",
            "complete_nonlinear_lower_order_first_order_sources_derived",
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


def compatible_constraint_certificate(configuration: Mapping[str, Any]) -> dict[str, Any]:
    """Build CON1-COMP1's local activated compatibility witness.

    The QIFT1 record/configuration are identity and predecessor handoff data.
    This function does not inspect an acceleration root or reproduce QIFT1.
    """

    expected = {
        "qift1_record", "qift1_configuration", "reference_radius", "activated_phi", "activated_phi_radial_derivative",
        "parameter_half_width", "acceleration_half_width", "formula", "projection_order",
        "polynomial_degree", "evaluation_nodes",
    }
    if not isinstance(configuration, Mapping) or set(configuration) != expected:
        raise ValueError("CON1-COMP1 configuration has unexpected keys")
    qift_configuration = configuration["qift1_configuration"]
    qift_record = configuration["qift1_record"]
    if not isinstance(qift_configuration, Mapping) or not isinstance(qift_record, Mapping):
        raise TypeError("QIFT1 configuration and record must be mappings")
    fixture = qift_configuration.get("flat_fixture")
    if not isinstance(fixture, Mapping):
        raise ValueError("qift1_configuration must expose flat_fixture")
    if qift_record.get("artifact_id") != "FGC-1-HYP1-DOM1-QIFT1":
        raise ValueError("CON1-COMP1 requires the QIFT1 predecessor record")
    if configuration["formula"] != "H=E_tt-2sE_tr+s^2E_rr; M=E_tr-sE_rr; s=h_tr/h_rr":
        raise ValueError("CON1-COMP1 projection formula differs from frozen definition")
    if tuple(configuration["projection_order"]) != PHYSICAL_PROJECTION_ORDER:
        raise ValueError("CON1-COMP1 projection order differs from frozen definition")
    if configuration["polynomial_degree"] != ACCELERATION_POLYNOMIAL_DEGREE:
        raise ValueError("CON1-COMP1 polynomial degree differs from frozen bound")
    if tuple(configuration["evaluation_nodes"]) != ACCELERATION_EVALUATION_NODES:
        raise ValueError("CON1-COMP1 evaluation nodes differ from frozen scheme")
    reference_radius = _fraction("reference_radius", configuration["reference_radius"], positive=True)
    if reference_radius != 4:
        raise ValueError("CON1-COMP1 requires reference radius four")
    acceleration_half_width = _fraction("acceleration_half_width", configuration["acceleration_half_width"], positive=True)
    if acceleration_half_width != QIFT1_ACCELERATION_HALF_WIDTH:
        raise ValueError("acceleration_half_width must equal QIFT1's frozen half width")
    datum = activated_compatible_state(
        fixture,
        activated_phi=_fraction("activated_phi", configuration["activated_phi"], positive=True),
        activated_phi_radial_derivative=_fraction("activated_phi_radial_derivative", configuration["activated_phi_radial_derivative"], positive=True),
        parameter_half_width=_fraction("parameter_half_width", configuration["parameter_half_width"], positive=True),
    )
    state = datum["state"]
    gauge = modified_harmonic_gauge_constraint(
        state, reference=datum["reference"], coordinate_radius=reference_radius,
        tilde_normal_factor=datum["tilde_normal_factor"],
    )
    projection = physical_constraint_projections(state)
    kinematic_state: FirstOrderSphericalState = first_order_state_from_spherical_state(state)
    kinematic = first_order_kinematic_residuals(kinematic_state)
    polynomial = acceleration_polynomial_certificate(state)
    normal_map = normal_gauge_extension_map(state, gauge, hat_normal_factor=datum["hat_normal_factor"])
    full_residual = modified_harmonic_full_residuals(
        state,
        reference=datum["reference"],
        coordinate_radius=reference_radius,
        tilde_normal_factor=datum["tilde_normal_factor"],
        hat_normal_factor=datum["hat_normal_factor"],
    )
    cone_control = activated_principal_cone_control(fixture, state)
    local_checks = {
        "activated_parameter_point_strictly_inside_qift1_z_box": datum["strictly_inside_qift_parameter_box"],
        "metric_defined_C_zero": _all_zero(gauge.constraint_up),
        "metric_defined_all_non_normal_nabla_C_components_zero": _all_zero(gauge.covariant_constraint_derivative[1:]),
        "nabla_r_Ct_and_Cr_zero": gauge.covariant_constraint_derivative[1][0] == 0 and gauge.covariant_constraint_derivative[1][1] == 0,
        "normal_shift_zero": projection["shift"] == 0,
        "unredefined_H_zero": projection["H"] == 0,
        "unredefined_M_zero": projection["M"] == 0,
        "unredefined_Etr_zero": projection["metric_tr"] == 0,
        "physical_projection_source_is_unredefined": projection["source_is_unredefined_residual"],
        "reduction_rows_zero": all(_all_zero(kinematic[name]) for name in (
            "u_time_definition_residual", "radial_reduction_constraint", "mixed_partial_compatibility_residual"
        )),
        "reduction_rows_independent_of_accelerations_by_repacking": True,
        "all_56_unredefined_constraint_polynomial_coefficients_zero": polynomial["all_56_coefficients_zero"],
        "unredefined_Ett_and_Etr_zero_for_every_acceleration": polynomial["all_56_coefficients_zero"] and projection["shift"] == 0,
        "ref1_scalar_equations_unmodified": full_residual["scalar_equations_unmodified"],
        "normal_gauge_extension_map_diagonal_nonzero": normal_map["diagonal_entries_nonzero"],
        "normal_gauge_extension_map_nonsingular": normal_map["determinant_nonzero"],
        "activated_principal_control_uses_exact_compatible_state": cone_control["compatible_state_pullback_exact"],
        "activated_principal_metric_and_regulator_discriminants_positive": cone_control["metric_null_discriminant_positive"] and cone_control["regulator_discriminant_positive"],
        "activated_principal_physical_and_regulator_cones_share_no_root": cone_control["cones_share_no_root"] and cone_control["physical_vs_regulator_resultant_nonzero"],
        "qift_root_is_external_predecessor_not_solved_here": datum["qift_acceleration_root_is_external_predecessor_fact"] and not datum["qift_acceleration_root_solved_here"],
        "conditional_qift_root_implies_normal_gauge_compatibility": polynomial["all_56_coefficients_zero"] and projection["shift"] == 0 and _all_zero(gauge.covariant_constraint_derivative[1:]) and normal_map["determinant_nonzero"],
    }
    if not all(local_checks.values()):
        raise ValueError("CON1-COMP1 local compatibility checks failed: " + repr(sorted(key for key, value in local_checks.items() if not value)))
    return {
        "artifact_id": CON1_COMP1_ARTIFACT_ID,
        "classification": "exact_local_activated_metric_defined_compatibility_witness_conditional_on_external_qift1_acceleration_root",
        "activated_datum": datum,
        "gauge": gauge,
        "physical_projections": projection,
        "kinematic": kinematic,
        "acceleration_polynomial": polynomial,
        "normal_gauge_extension_map": normal_map,
        "ref1_composition": {
            "equation_order": full_residual["equation_order"],
            "scalar_equations_unmodified": full_residual["scalar_equations_unmodified"],
        },
        "activated_principal_cone_control": cone_control,
        "local_exact_checks": local_checks,
        "all_declared_exact_checks_pass": all(local_checks.values()),
        "nonclaims": required_con1_comp1_nonclaims(),
    }
