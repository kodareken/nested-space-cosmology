"""Closed spherical physical, gauge, and reduction constraint contract.

The modified-harmonic REF1 equations evolve six spherical fields, but their
metric rows are not themselves the physical ACT1 equations away from the
metric-defined gauge shell.  This module keeps the distinction explicit.

It supplies four ingredients used by CON4-PHY1:

* the Hamiltonian and radial-momentum projections of the *unredefined*
  ACT1/VAR1 metric equation, decomposed into their signed physical terms;
* the structural (Gauss--Codazzi/double-dual) reason those projections contain
  no coordinate-time accelerations;
* an exact, general spherical map from those two physical constraints to the
  normal derivative of the REF1 gauge vector on the REF1 metric shell; and
* executable scale-free cancellation monitors for physical, gauge, and FO1
  kinematic-reduction constraints.

Together with CON2's homogeneous metric-derived gauge subsidiary identity and
CON3's conditional boundary-free uniqueness theorem, the map closes physical
constraint propagation for any sufficiently smooth compatible solution before
a boundary enters.  It does not construct that solution, a compatible initial
slice, a regular centre, or a constraint-preserving boundary condition.
"""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
from itertools import product
from numbers import Integral
from typing import Any, Mapping, Sequence

from .exact_linear_algebra import matrix_det, matrix_inverse
from .exact_tangent import FirstTangent, primal_and_tangent
from .modified_harmonic import (
    auxiliary_inverse_metric,
    inverse_metric_null_polynomial,
    quadratic_companion_certificate,
)
from .modified_harmonic_constraints import physical_constraint_projections
from .modified_harmonic_reduction_subsidiary import ReductionDifferentialState
from .modified_harmonic_reference import (
    ModifiedHarmonicGaugeData,
    modified_harmonic_extension_residual,
    physical_connection_data,
)
from .spherical_reduction import BASE_FIELD_ORDER, SphericalState, residuals


Q = Fraction
N = 4
ACTIVE_SPHERICAL_GAUGE_COMPONENTS = ("C^t", "C^r")
PHYSICAL_CONSTRAINT_TERM_ORDER = (
    "F_times_Einstein",
    "nonminimal_F",
    "minus_phi_stress",
    "minus_chi_stress",
    "Gauss_Bonnet_response",
)


def _fraction(name: str, value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational")
    return value if isinstance(value, Fraction) else Q(int(value))


def _project_tensor(
    tensor: Sequence[Sequence[Fraction]], shift: Fraction
) -> tuple[Fraction, Fraction]:
    """Project a covariant spherical tensor on ``w=(1,-shift)`` and ``dr``."""

    if len(tensor) != N or any(len(row) != N for row in tensor):
        raise ValueError("physical constraint term must be a four-by-four tensor")
    hamiltonian = tensor[0][0] - 2 * shift * tensor[0][1] + shift**2 * tensor[1][1]
    momentum = tensor[0][1] - shift * tensor[1][1]
    return hamiltonian, momentum


def _negate_tensor(
    tensor: Sequence[Sequence[Fraction]],
) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(tuple(-entry for entry in row) for row in tensor)


def relative_cancellation_monitor(
    terms: Sequence[Fraction | int],
    *,
    expected_value: Fraction | int | None = None,
) -> dict[str, Any]:
    """Return a bounded, dimensionless cancellation monitor.

    The denominator is the sum of the absolute signed contributions.  If all
    contributions vanish exactly, the relative value is defined to be zero
    and the zero-scale convention is recorded.  This monitor detects loss of
    cancellation but deliberately supplies no absolute truncation-error scale;
    HLT1/NUM1 must add that independent scale before using a threshold.
    """

    values = tuple(_fraction(f"terms[{index}]", value) for index, value in enumerate(terms))
    if not values:
        raise ValueError("relative cancellation monitor requires at least one term")
    value = sum(values, Q(0))
    if expected_value is not None and value != _fraction("expected_value", expected_value):
        raise ValueError("signed monitor terms do not reconstruct the expected value")
    denominator = sum((abs(item) for item in values), Q(0))
    zero_scale = denominator == 0
    ratio = Q(0) if zero_scale else abs(value) / denominator
    if not Q(0) <= ratio <= Q(1):
        raise ValueError("relative cancellation monitor escaped its exact unit interval")
    return {
        "raw_value": value,
        "signed_terms": values,
        "absolute_term_sum": denominator,
        "relative_cancellation": ratio,
        "zero_scale_convention_applied": zero_scale,
        "dimensionless_and_bounded_by_one": True,
        "absolute_error_scale_supplied_here": False,
    }


def physical_constraint_term_certificate(state: SphericalState) -> dict[str, Any]:
    """Evaluate and exactly recompose the unredefined physical constraints."""

    if not isinstance(state, SphericalState):
        raise TypeError("state must be a SphericalState")
    source = residuals(state)
    projection = physical_constraint_projections(state)
    shift = projection["shift"]
    tensors = (
        source["einstein_term"],
        source["nonminimal_term"],
        _negate_tensor(source["phi_stress"]),
        _negate_tensor(source["chi_stress"]),
        source["gb_residual_term"],
    )
    projected = tuple(_project_tensor(tensor, shift) for tensor in tensors)
    h_terms = tuple(item[0] for item in projected)
    m_terms = tuple(item[1] for item in projected)
    h_monitor = relative_cancellation_monitor(h_terms, expected_value=projection["H"])
    m_monitor = relative_cancellation_monitor(m_terms, expected_value=projection["M"])

    h_tt = state.h_tt.value
    h_tr = state.h_tr.value
    h_rr = state.h_rr.value
    lapse_squared = -h_tt + h_tr**2 / h_rr
    if h_rr <= 0 or lapse_squared <= 0:
        raise ValueError("constraint projections require a spacelike coordinate slice")
    return {
        "source_equations": "unredefined_ACT1_VAR1",
        "normalization": {
            "unnormalized_normal_w": (Q(1), -shift, Q(0), Q(0)),
            "shift_s": shift,
            "lapse_squared": lapse_squared,
            "radial_spatial_metric_h_rr": h_rr,
            "unit_H_equals_H_over_lapse_squared": True,
            "unit_M_equals_M_over_lapse_times_sqrt_h_rr": True,
            "zero_sets_identical_under_positive_factors": True,
        },
        "term_order": PHYSICAL_CONSTRAINT_TERM_ORDER,
        "Hamiltonian_terms": h_terms,
        "momentum_terms": m_terms,
        "Hamiltonian": projection["H"],
        "momentum": projection["M"],
        "Hamiltonian_exact_recomposition": sum(h_terms, Q(0)) == projection["H"],
        "momentum_exact_recomposition": sum(m_terms, Q(0)) == projection["M"],
        "Hamiltonian_monitor": h_monitor,
        "momentum_monitor": m_monitor,
    }


def _double_dual_support(first_normal: int, second_index: int) -> tuple[tuple[int, ...], ...]:
    """Enumerate nonzero abstract support of ``eps[a,c,e,f] eps[b,d,g,h]``."""

    support: list[tuple[int, ...]] = []
    for c, d, e, f, g, h in product(range(N), repeat=6):
        if len({first_normal, c, e, f}) == N and len({second_index, d, g, h}) == N:
            support.append((c, d, e, f, g, h))
    return tuple(support)


def constraint_acceleration_structure_certificate() -> dict[str, Any]:
    """Record the covariant reason H/M contain no double-normal derivatives.

    For the Gauss--Bonnet term this performs the finite index-support part of
    the proof directly.  The double-dual representation

    ``P_acbd proportional eps_acEF eps_bdGH R^EFGH``

    shows which curvature and scalar-Hessian indices can survive the normal
    projections.  The remaining Einstein and nonminimal statements are the
    standard Gauss--Codazzi and normal/spatial Hessian decompositions written
    explicitly in the returned certificate.
    """

    h_support = _double_dual_support(0, 0)
    m_support = _double_dual_support(0, 1)
    if not h_support or not m_support:
        raise AssertionError("double-dual support enumeration is unexpectedly empty")

    h_hessian_spatial = all(c != 0 and d != 0 for c, d, *_ in h_support)
    h_curvature_spatial = all(
        all(index != 0 for index in (e, f, g, h))
        for _, _, e, f, g, h in h_support
    )
    m_hessian_at_most_one_normal = all(
        sum(index == 0 for index in (c, d)) <= 1
        for c, d, *_ in m_support
    )
    m_curvature_at_most_one_normal = all(
        sum(index == 0 for index in (e, f, g, h)) <= 1
        for _, _, e, f, g, h in m_support
    )
    if not all(
        (
            h_hessian_spatial,
            h_curvature_spatial,
            m_hessian_at_most_one_normal,
            m_curvature_at_most_one_normal,
        )
    ):
        raise ValueError("double-dual support does not establish the declared split")

    return {
        "classification": "covariant_normal_projection_acceleration_structure_theorem",
        "normal_projection_basis": "orthonormal_n_plus_three_spatial_directions",
        "Einstein_sector": {
            "Hamiltonian_identity": "G_nn=(3R+K^2-K_ij*K^ij)/2",
            "momentum_identity": "G_ni=D_j(K^j_i)-D_i(K)",
            "contains_second_normal_metric_derivatives": False,
        },
        "nonminimal_F_sector": {
            "Hamiltonian_identity": "n^a*n^b*(g_ab*boxF-nabla_a_nabla_bF)=-h^ij*nabla_i_nabla_jF",
            "momentum_identity": "n^a*h_i^b*(g_ab*boxF-nabla_a_nabla_bF)=-n^a*h_i^b*nabla_a_nabla_bF",
            "contains_second_normal_scalar_derivatives": False,
        },
        "matter_stress_sector": {
            "depends_on_fields_and_first_derivatives_only": True,
            "contains_coordinate_time_accelerations": False,
        },
        "Gauss_Bonnet_sector": {
            "double_dual_identity": "P_acbd proportional eps_acEF*eps_bdGH*R^EFGH",
            "Hamiltonian_support_count": len(h_support),
            "momentum_support_count_for_one_spatial_direction": len(m_support),
            "Hamiltonian_hessian_indices_all_spatial": h_hessian_spatial,
            "Hamiltonian_curvature_indices_all_spatial": h_curvature_spatial,
            "momentum_hessian_has_at_most_one_normal_index": m_hessian_at_most_one_normal,
            "momentum_curvature_has_at_most_one_normal_index": m_curvature_at_most_one_normal,
            "contains_double_normal_metric_or_scalar_derivatives": False,
        },
        "Hamiltonian_and_momentum_are_initial_data_constraints": True,
        "coordinate_time_accelerations_absent_structurally": True,
        "independent_exact_evaluator_controls_still_required": True,
    }


def acceleration_independence_control(states: Sequence[SphericalState]) -> dict[str, Any]:
    """Differentiate H/M exactly in all six acceleration directions at controls."""

    if not states:
        raise ValueError("at least one exact state control is required")
    records: list[dict[str, Any]] = []
    for state_index, state in enumerate(states):
        if not isinstance(state, SphericalState):
            raise TypeError("every acceleration control must be a SphericalState")
        base = physical_constraint_projections(state)
        directions: list[dict[str, Any]] = []
        for field in BASE_FIELD_ORDER:
            jet = getattr(state, field)
            seeded = replace(
                state,
                **{field: replace(jet, dtt=FirstTangent(jet.dtt, Q(1)))},
            )
            projection = physical_constraint_projections(seeded)
            h_primal, h_tangent = primal_and_tangent(projection["H"])
            m_primal, m_tangent = primal_and_tangent(projection["M"])
            if h_primal != base["H"] or m_primal != base["M"]:
                raise ValueError("acceleration tangent changed its exact primal constraint")
            directions.append(
                {
                    "acceleration": f"partial_t^2.{field}",
                    "dH_dacceleration": h_tangent,
                    "dM_dacceleration": m_tangent,
                    "both_exactly_zero": h_tangent == 0 and m_tangent == 0,
                }
            )
        if not all(item["both_exactly_zero"] for item in directions):
            raise ValueError(f"exact acceleration-independence control {state_index} failed")
        records.append(
            {
                "state_index": state_index,
                "base_H": base["H"],
                "base_M": base["M"],
                "directions": tuple(directions),
                "all_twelve_derivatives_zero": True,
            }
        )
    return {
        "state_count": len(records),
        "direction_count_per_state": len(BASE_FIELD_ORDER),
        "records": tuple(records),
        "all_exact_controls_passed": True,
        "controls_are_regressions_not_the_global_proof": True,
    }


def physical_to_normal_gauge_map(
    state: SphericalState,
    gauge: ModifiedHarmonicGaugeData,
    *,
    hat_normal_factor: Fraction | int,
) -> dict[str, Any]:
    """Derive the exact spherical H/M to normal-``nabla C`` map.

    Assume ``C=0`` on a hypersurface.  Its spatial covariant derivatives then
    vanish there.  With ``x=nabla_t C^t`` and ``y=nabla_t C^r``, evaluate the
    REF1 extension on two exact basis inputs and project ``E=-X``.  The result
    is compared entry-by-entry with the analytic ADM formula valid at every
    spherical spacelike point with ``F>0`` and ``hat_normal_factor>1``.
    """

    if not isinstance(state, SphericalState):
        raise TypeError("state must be a SphericalState")
    if not isinstance(gauge, ModifiedHarmonicGaugeData):
        raise TypeError("gauge must be ModifiedHarmonicGaugeData")
    factor = _fraction("hat_normal_factor", hat_normal_factor)
    if factor <= 1:
        raise ValueError("hat normal factor must exceed one")
    physical = physical_connection_data(state)
    if gauge.physical != physical:
        raise ValueError("gauge and physical state differ")

    h_tt = state.h_tt.value
    h_tr = state.h_tr.value
    h_rr = state.h_rr.value
    if h_rr <= 0:
        raise ValueError("coordinate slice must have positive radial metric")
    shift = h_tr / h_rr
    lapse_squared = -h_tt + h_tr**2 / h_rr
    if lapse_squared <= 0:
        raise ValueError("coordinate slice must be spacelike")
    effective_planck = state.planck_mass**2 + state.beta * state.phi.value**2
    if effective_planck <= 0:
        raise ValueError("effective Planck coefficient must be positive")

    zero = tuple(tuple(Q(0) for _ in range(N)) for _ in range(N))
    columns: list[tuple[Fraction, Fraction]] = []
    extension_columns: list[tuple[tuple[Fraction, ...], ...]] = []
    for component in range(2):
        derivative = [list(row) for row in zero]
        derivative[0][component] = Q(1)
        basis_gauge = replace(
            gauge,
            covariant_constraint_derivative=tuple(tuple(row) for row in derivative),
        )
        extension = modified_harmonic_extension_residual(
            state,
            gauge=basis_gauge,
            hat_normal_factor=factor,
        )["covariant"]
        # REF1 metric shell: E_unredefined + X_extension = 0.
        physical_residual = _negate_tensor(extension)
        columns.append(_project_tensor(physical_residual, shift))
        extension_columns.append(extension)
    computed = tuple(
        tuple(columns[column][row] for column in range(2)) for row in range(2)
    )
    expected = (
        (effective_planck * factor * lapse_squared / 2, Q(0)),
        (-effective_planck * factor * h_tr / 2, -effective_planck * factor * h_rr / 2),
    )
    if computed != expected:
        raise ValueError("direct REF1 extension and analytic physical-constraint map differ")
    determinant = matrix_det(computed)
    expected_determinant = -effective_planck**2 * factor**2 * lapse_squared * h_rr / 4
    if determinant != expected_determinant or determinant >= 0:
        raise ValueError("normal gauge/physical constraint map lost its exact determinant")
    inverse = matrix_inverse(computed)
    return {
        "input_order": ("nabla_t_C^t", "nabla_t_C^r"),
        "output_order": ("unredefined_H_on_REF1_shell", "unredefined_M_on_REF1_shell"),
        "computed_matrix": computed,
        "analytic_matrix": expected,
        "direct_and_analytic_matrices_equal": True,
        "determinant": determinant,
        "analytic_determinant": expected_determinant,
        "determinant_strictly_negative": True,
        "inverse_matrix": inverse,
        "domain": {
            "F_strictly_positive": True,
            "hat_normal_factor_greater_than_one": True,
            "lapse_squared_strictly_positive": True,
            "h_rr_strictly_positive": True,
            "C_zero_on_slice_required": True,
            "spatial_covariant_derivative_of_C_zero_on_slice_follows": True,
            "complete_REF1_metric_equations_required": True,
        },
        "normal_derivative_equivalence": {
            "unit_normal": "n=(partial_t-s*partial_r)/sqrt(lapse_squared)",
            "spatial_nabla_C_zero_makes_nabla_n_C_equal_nabla_t_C_over_sqrt_lapse_squared": True,
            "H_and_M_zero_iff_normal_nabla_C_zero_on_REF1_shell": True,
        },
        "extension_basis_columns": tuple(extension_columns),
    }


def gauge_constraint_characteristics(
    state: SphericalState,
    *,
    hat_normal_factor: Fraction | int,
) -> dict[str, Any]:
    """Return the two radial hat-cone branches of the gauge subsidiary system."""

    if not isinstance(state, SphericalState):
        raise TypeError("state must be a SphericalState")
    factor = _fraction("hat_normal_factor", hat_normal_factor)
    physical = physical_connection_data(state)
    hat = auxiliary_inverse_metric(physical.inverse_metric, factor)
    polynomial = inverse_metric_null_polynomial(hat)
    roots = quadratic_companion_certificate(polynomial)["roots"]

    def sign(root: Mapping[str, Fraction]) -> int:
        if "exact" in root:
            value = root["exact"]
            return -1 if value < 0 else 1 if value > 0 else 0
        lower, upper = root["lower"], root["upper"]
        return -1 if upper < 0 else 1 if lower > 0 else 0

    signs = tuple(sign(root) for root in roots)
    if any(value == 0 for value in signs):
        raise ValueError("gauge characteristic control reaches a zero radial speed")
    negative = sum(value < 0 for value in signs)
    positive = sum(value > 0 for value in signs)
    return {
        "principal_equation": "(F/2)*hat_g^ab*nabla_a*nabla_b*C^mu plus lower order equals zero",
        "radial_covector_convention": "xi_a=(-c,1)",
        "hat_null_polynomial": polynomial,
        "root_isolations": roots,
        "root_signs": signs,
        "active_spherical_component_order": ACTIVE_SPHERICAL_GAUGE_COMPONENTS,
        "each_active_component_has_both_hat_cone_branches": True,
        "outer_incoming_component_count_at_control": len(ACTIVE_SPHERICAL_GAUGE_COMPONENTS) * negative,
        "inner_incoming_component_count_at_control": len(ACTIVE_SPHERICAL_GAUGE_COMPONENTS) * positive,
        "angular_gauge_components_identically_absent_by_spherical_symmetry": True,
        "physical_H_and_M_are_algebraic_normal_gauge_data_on_REF1_shell": True,
        "reduction_constraints_have_zero_radial_principal_speed": True,
        "constraint_preserving_boundary_map_to_main_fields_derived": False,
    }


def gauge_constraint_monitor_bundle(gauge: ModifiedHarmonicGaugeData) -> dict[str, Any]:
    """Reassemble ``C`` and ``nabla C`` into exact cancellation monitors."""

    if not isinstance(gauge, ModifiedHarmonicGaugeData):
        raise TypeError("gauge must be ModifiedHarmonicGaugeData")
    c_records: list[dict[str, Any]] = []
    derivative_records: list[dict[str, Any]] = []
    for component in range(N):
        terms = tuple(
            -gauge.tilde_inverse_metric[rho][sigma]
            * gauge.gamma_difference[component][rho][sigma]
            for rho in range(N)
            for sigma in range(N)
        )
        c_records.append(
            {
                "component": component,
                "monitor": relative_cancellation_monitor(
                    terms, expected_value=gauge.constraint_up[component]
                ),
            }
        )
        for direction in range(N):
            partial_terms = tuple(
                -(
                    gauge.tilde_inverse_metric_derivative[direction][rho][sigma]
                    * gauge.gamma_difference[component][rho][sigma]
                    + gauge.tilde_inverse_metric[rho][sigma]
                    * gauge.gamma_difference_derivative[direction][component][rho][sigma]
                )
                for rho in range(N)
                for sigma in range(N)
            )
            connection_terms = tuple(
                gauge.physical.christoffel[component][direction][source]
                * gauge.constraint_up[source]
                for source in range(N)
            )
            terms = partial_terms + connection_terms
            derivative_records.append(
                {
                    "direction": direction,
                    "component": component,
                    "partial_term_count": len(partial_terms),
                    "connection_term_count": len(connection_terms),
                    "monitor": relative_cancellation_monitor(
                        terms,
                        expected_value=gauge.covariant_constraint_derivative[direction][component],
                    ),
                }
            )
    return {
        "constraint_C": tuple(c_records),
        "covariant_derivative_nabla_C": tuple(derivative_records),
        "all_components_recomposed_exactly": True,
        "absolute_threshold_deferred_to_HLT1_NUM1": True,
    }


def reduction_constraint_monitor_bundle(
    state: ReductionDifferentialState,
) -> dict[str, Any]:
    """Return exact monitors for every FO1 kinematic and subsidiary residual."""

    if not isinstance(state, ReductionDifferentialState):
        raise TypeError("state must be a ReductionDifferentialState")
    records: list[dict[str, Any]] = []
    for field in BASE_FIELD_ORDER:
        jet = getattr(state, field)
        d_terms = (jet.u_t, -jet.p)
        c_terms = (jet.q, -jet.u_r)
        k_terms = (jet.q_t, -jet.p_r)
        d_r_d_terms = (jet.u_tr, -jet.p_r)
        d_t_c_terms = (jet.q_t, -jet.u_tr)
        identity_terms = (
            sum(d_t_c_terms, Q(0)),
            sum(d_r_d_terms, Q(0)),
            -sum(k_terms, Q(0)),
        )
        records.append(
            {
                "field": field,
                "D_equals_partial_t_u_minus_p": relative_cancellation_monitor(d_terms),
                "C_equals_q_minus_partial_r_u": relative_cancellation_monitor(c_terms),
                "K_equals_partial_t_q_minus_partial_r_p": relative_cancellation_monitor(k_terms),
                "subsidiary_identity_d_t_C_plus_d_r_D_minus_K": relative_cancellation_monitor(
                    identity_terms, expected_value=Q(0)
                ),
            }
        )
    return {
        "field_order": BASE_FIELD_ORDER,
        "records": tuple(records),
        "off_shell_subsidiary_identity_recomposed_exactly": True,
        "subsidiary_radial_principal_speed": Q(0),
        "incoming_reduction_constraint_fields_at_either_boundary": (),
        "absolute_threshold_deferred_to_HLT1_NUM1": True,
    }


def constraint_monitor_contract() -> dict[str, Any]:
    """Declare the common numerical interpretation of all CON4 monitors."""

    return {
        "definition": "abs(sum signed terms)/sum(abs(signed terms))",
        "range": "closed_unit_interval",
        "zero_denominator_policy": "if_every_term_is_exactly_zero_return_zero_and_record_zero_scale",
        "raw_residual_and_denominator_always_recorded": True,
        "common_rescaling_invariant": True,
        "not_an_absolute_truncation_or_physical_error_norm": True,
        "HLT1_and_NUM1_must_add_absolute_scale_and_thresholds": True,
        "smallness_on_one_run_is_not_a_propagation_proof": True,
    }


def conditional_constraint_closure_statement() -> dict[str, Any]:
    """Return the exact logical composition used by CON4-PHY1."""

    return {
        "equations_solved_by_future_evolution": (
            "six_complete_REF1_metric_plus_unmodified_scalar_equations",
            "six_u_time_definition_equations",
            "six_q_time_minus_p_radial_compatibility_equations",
        ),
        "constraints_imposed_on_initial_data": (
            "metric_defined_C^t=C^r=0",
            "unredefined_H=0",
            "unredefined_M=0",
            "six_radial_reduction_constraints_q-partial_r_u=0",
        ),
        "constraints_propagated_analytically": (
            "C_by_CON2_homogeneous_hat_wave_plus_CON3_boundary_free_uniqueness",
            "H_and_M_by_invertible_REF1_shell_map_and_vanishing_extension",
            "six_reduction_constraints_by_d_t_C+d_r_D-K=0",
        ),
        "constraints_monitored_numerically": (
            "unredefined_H_and_M",
            "metric_defined_C_and_nabla_C",
            "six_D_C_K_reduction_residuals",
        ),
        "logical_chain": (
            "C_zero_on_initial_slice_implies_spatial_nabla_C_zero",
            "REF1_shell_plus_H_equals_M_equals_zero_implies_normal_nabla_C_zero_by_invertible_map",
            "CON2_plus_both_scalar_equations_gives_homogeneous_hat_cone_subsidiary_system",
            "CON3_uniqueness_gives_C_zero_in_the_boundary_free_domain_of_dependence",
            "extension_then_vanishes_so_REF1_recovers_unredefined_ACT1_VAR1_and_H_equals_M_equals_zero",
            "FO1_kinematic_identity_independently_preserves_six_radial_reduction_constraints",
        ),
        "complete_boundary_free_spherical_constraint_system_closed_conditionally": True,
        "smooth_solution_existence_supplied": False,
        "regular_center_supplied": False,
        "compatible_initial_hypersurface_supplied": False,
        "constraint_preserving_boundary_map_supplied": False,
        "initial_boundary_value_problem_supplied": False,
        "evolution_authorized": False,
    }
