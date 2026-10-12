from __future__ import annotations

from dataclasses import replace
from fractions import Fraction as Q
from math import inf, nextafter
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recursive_horizons.fgc.def1_stab1 import (  # noqa: E402
    AFFINE_NORMALIZATION_RESIDUAL_MAX,
    ActivationAssessment,
    AffineIntervalAssessment,
    AssembledQErrorBudget,
    DEF1_BOOLEAN_NAMES,
    DIRECT_RAYCHAUDHURI_AGREEMENT_MAX,
    Def1BooleanRecord,
    Def1Stab1Error,
    ERROR_BUDGET_COMPONENTS,
    GR0_PHI_POLICY_CANONICAL_ZERO,
    GR0_PHI_POLICY_HISTORICAL_PLANTED_SEED,
    MIN_ACTIVATION_LOWER_BOUND,
    MIN_AFFINE_INTERVAL_OVER_L0,
    MIN_AFFINE_SAMPLES,
    MassFluxAssessment,
    MatchedControlSample,
    PREMISE_STATUS_CONDITIONAL,
    PREMISE_STATUS_PROVEN,
    PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN,
    QComponentPremise,
    RaychaudhuriErrorBudget,
    SOURCE_IMP1_ADMISSION_DEBIT,
    SOURCE_RICHARDSON,
    SOURCE_SUPPLIED_CERTIFIED,
    STOPPED_MASS_FLUX_INCONSISTENCY,
    activation_threshold_passed,
    affine_interval_and_sample_gate,
    affine_normalization_residual_passed,
    assemble_q_error_budget,
    assess_activation,
    assess_mass_flux_ledger,
    complete_q_margin_passed,
    control_dominance_passed,
    declared_method_count_meets_minimum,
    declared_resolution_count_meets_minimum,
    direct_raychaudhuri_routes_agree,
    inject_component_input_error,
    inverse_base_metric,
    inverse_base_metric_from_adm,
    map_mass_flux_inconsistency_token,
    misner_sharp_equation_defect,
    misner_sharp_identity_holds,
    misner_sharp_mass,
    misner_sharp_mass_gradient,
    mix_base_tensor,
    observed_order_meets_minimum,
    record_def1_booleans,
    record_def1_booleans_from_mapping,
    trace_adjusted_radial_projector,
    trappedness_margin_passed,
    zero_error_premises,
)
from recursive_horizons.fgc.exact_interval import Interval  # noqa: E402
from recursive_horizons.fgc.metric_null_observable import (  # noqa: E402
    SphericalADMGeometryJet,
    metric_null_raychaudhuri_point_certificate,
)
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    Jet2,
    SphericalState,
    direct_4d_curvature,
)


N = 4


def _minkowski_geometry() -> SphericalADMGeometryJet:
    return SphericalADMGeometryJet(
        lapse=Jet2.constant(1),
        shift=Jet2.constant(0),
        radial_scale=Jet2.constant(1),
        areal_radius=Jet2(4, dr=1),
    )


def _de_sitter_geometry() -> SphericalADMGeometryJet:
    return SphericalADMGeometryJet(
        lapse=Jet2.constant(1),
        shift=Jet2.constant(0),
        radial_scale=Jet2(1, dt=-1, dtt=1),
        areal_radius=Jet2(2, dt=-2, dr=1, dtt=2, dtr=-1),
    )


def _synthetic_positive_q_geometry() -> SphericalADMGeometryJet:
    return SphericalADMGeometryJet(
        lapse=Jet2.constant(1),
        shift=Jet2.constant(0),
        radial_scale=Jet2.constant(1),
        areal_radius=Jet2(4, dt=-1, dr=1, dtt=1),
    )


def _base_einstein(state: SphericalState):
    riemann, inverse, _ = direct_4d_curvature(state)
    metric = [
        [state.h_tt.value, state.h_tr.value, Q(0), Q(0)],
        [state.h_tr.value, state.h_rr.value, Q(0), Q(0)],
        [Q(0), Q(0), state.areal_radius.value**2, Q(0)],
        [Q(0), Q(0), Q(0), state.areal_radius.value**2],
    ]
    ricci = [
        [
            sum(
                inverse[a][c] * riemann[a][b][c][d]
                for a in range(N)
                for c in range(N)
            )
            for d in range(N)
        ]
        for b in range(N)
    ]
    scalar = sum(inverse[a][b] * ricci[a][b] for a in range(N) for b in range(N))
    covariant = tuple(
        tuple(ricci[a][b] - metric[a][b] * scalar / 2 for b in range(2))
        for a in range(2)
    )
    inverse_base = tuple(tuple(inverse[a][b] for b in range(2)) for a in range(2))
    return covariant, inverse_base, scalar


def _matched_controls(*, gr0_phi: Q = Q(0), sgbl_phi: Q = Q(0), error: Q = Q(0)):
    return (
        MatchedControlSample("GR-0", gr0_phi, error),
        MatchedControlSample("SGB-L", sgbl_phi, error),
    )


class ErrorBudgetContractTests(unittest.TestCase):
    def test_frozen_field_order_and_nonnegative_bounds(self) -> None:
        assembled = assemble_q_error_budget(zero_error_premises())
        self.assertEqual(
            tuple(assembled.budget.__dataclass_fields__),
            ERROR_BUDGET_COMPONENTS + ("total_upper_bound",),
        )
        self.assertEqual(
            tuple(ActivationAssessment.__dataclass_fields__),
            (
                "raw_factor",
                "lower_bound",
                "matched_control_upper_bound",
                "threshold_passed",
                "control_dominance_passed",
            ),
        )
        for name in ERROR_BUDGET_COMPONENTS + ("total_upper_bound",):
            value = getattr(assembled.budget, name)
            self.assertIsInstance(value, float)
            self.assertGreaterEqual(value, 0.0)
            self.assertTrue(value == value and abs(value) != inf)
        self.assertEqual(assembled.exact_total, 0)
        self.assertFalse(assembled.used_measured_q)
        self.assertFalse(assembled.global_pde_error_certified)
        self.assertFalse(assembled.richardson_treated_as_global_pde_error)
        self.assertFalse(assembled.imp1_admission_debit_treated_as_global_pde_error)
        self.assertTrue(assembled.all_components_declared_proven_enclosures)

    def test_missing_and_overlapping_premises_refuse(self) -> None:
        premises = list(zero_error_premises())
        with self.assertRaisesRegex(Def1Stab1Error, "exactly the frozen"):
            assemble_q_error_budget(premises[:-1])
        with self.assertRaisesRegex(Def1Stab1Error, "duplicate Q-error component"):
            assemble_q_error_budget(tuple(premises) + (premises[0],))
        with self.assertRaisesRegex(Def1Stab1Error, "unknown Q-error component"):
            QComponentPremise(
                component="pde_global",
                sensitivity=Interval.singleton(1),
                input_error=Interval.singleton(0),
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
            )

    def test_negative_nonfinite_malformed_and_bool_premises_refuse(self) -> None:
        with self.assertRaisesRegex(Def1Stab1Error, "nonnegative"):
            QComponentPremise(
                component="affine",
                sensitivity=Interval.singleton(1),
                input_error=Interval(-1, 0),
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
            )
        with self.assertRaisesRegex(ValueError, "must not exceed"):
            Interval(1, 0)
        with self.assertRaises(TypeError):
            QComponentPremise(
                component="affine",
                sensitivity=True,
                input_error=Interval.singleton(0),
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
            )
        with self.assertRaises(TypeError):
            assemble_q_error_budget("missing")
        with self.assertRaisesRegex(Def1Stab1Error, "finite and nonnegative"):
            RaychaudhuriErrorBudget(
                **{name: 0.0 for name in ERROR_BUDGET_COMPONENTS},
                total_upper_bound=float("nan"),
            )
        with self.assertRaisesRegex(Def1Stab1Error, "finite and nonnegative"):
            RaychaudhuriErrorBudget(
                **{name: 0.0 for name in ERROR_BUDGET_COMPONENTS},
                total_upper_bound=-0.0 if False else -1.0,
            )

    def test_missing_component_is_not_filled_with_zero(self) -> None:
        premises = [
            premise
            for premise in zero_error_premises()
            if premise.component != "conservation"
        ]
        with self.assertRaisesRegex(Def1Stab1Error, "missing=\\['conservation'\\]"):
            assemble_q_error_budget(premises)
        with self.assertRaisesRegex(Def1Stab1Error, "missing component conservation"):
            inject_component_input_error(premises, "conservation", Q(1, 8))

    def test_bounds_do_not_shrink_with_positive_q_or_measured_q_argument(self) -> None:
        baseline = assemble_q_error_budget(zero_error_premises())
        inflated = assemble_q_error_budget(
            inject_component_input_error(
                zero_error_premises(), "physical_constraint", Q(3, 7)
            )
        )
        self.assertGreater(inflated.exact_total, baseline.exact_total)
        self.assertGreater(
            inflated.budget.total_upper_bound, baseline.budget.total_upper_bound
        )
        for (name, baseline_value), (_, inflated_value) in zip(
            baseline.exact_components, inflated.exact_components, strict=True
        ):
            if name == "physical_constraint":
                self.assertGreater(inflated_value, baseline_value)
            else:
                self.assertEqual(inflated_value, baseline_value)
        with self.assertRaisesRegex(Def1Stab1Error, "measured Q"):
            assemble_q_error_budget(zero_error_premises(), measured_q=Q(1000))
        with self.assertRaisesRegex(Def1Stab1Error, "measured Q"):
            assemble_q_error_budget(zero_error_premises(), measured_q=Q(1, 1000))
        self.assertEqual(
            assemble_q_error_budget(zero_error_premises()).exact_total,
            baseline.exact_total,
        )

    def test_richardson_and_imp1_remain_conditional_not_pde_error(self) -> None:
        with self.assertRaisesRegex(Def1Stab1Error, "conditional premise"):
            QComponentPremise(
                component="spatial_temporal",
                sensitivity=Interval.singleton(1),
                input_error=Interval.singleton(Q(1, 8)),
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_RICHARDSON,
            )
        with self.assertRaisesRegex(Def1Stab1Error, "conditional premise"):
            QComponentPremise(
                component="spatial_temporal",
                sensitivity=Interval.singleton(1),
                input_error=Interval.singleton(Q(1, 8)),
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_IMP1_ADMISSION_DEBIT,
            )
        with self.assertRaisesRegex(Def1Stab1Error, "only be declared on spatial_temporal"):
            QComponentPremise(
                component="conservation",
                sensitivity=Interval.singleton(1),
                input_error=Interval.singleton(Q(1, 8)),
                status=PREMISE_STATUS_CONDITIONAL,
                source=SOURCE_RICHARDSON,
            )
        premises = []
        for premise in zero_error_premises():
            if premise.component == "spatial_temporal":
                premises.append(
                    QComponentPremise(
                        component="spatial_temporal",
                        sensitivity=Interval.singleton(2),
                        input_error=Interval.singleton(Q(1, 5)),
                        status=PREMISE_STATUS_CONDITIONAL,
                        source=SOURCE_IMP1_ADMISSION_DEBIT,
                    )
                )
            else:
                premises.append(premise)
        assembled = assemble_q_error_budget(premises)
        self.assertEqual(assembled.exact_components[0][1], Q(2, 5))
        self.assertFalse(assembled.all_components_declared_proven_enclosures)
        self.assertFalse(assembled.global_pde_error_certified)
        self.assertFalse(assembled.imp1_admission_debit_treated_as_global_pde_error)
        with self.assertRaisesRegex(Def1Stab1Error, "not global PDE error"):
            AssembledQErrorBudget(
                budget=assembled.budget,
                exact_components=assembled.exact_components,
                exact_total=assembled.exact_total,
                component_statuses=assembled.component_statuses,
                component_sources=assembled.component_sources,
                all_components_declared_proven_enclosures=False,
                global_pde_error_certified=True,
                used_measured_q=False,
                richardson_treated_as_global_pde_error=False,
                imp1_admission_debit_treated_as_global_pde_error=False,
            )

    def test_upper_bound_rounds_outward_from_exact_rationals(self) -> None:
        premises = []
        for premise in zero_error_premises():
            if premise.component == "interpolation":
                premises.append(
                    QComponentPremise(
                        component="interpolation",
                        sensitivity=Interval(-Q(2), Q(5, 2)),
                        input_error=Interval(Q(0), Q(1, 3)),
                        status=PREMISE_STATUS_PROVEN,
                        source=SOURCE_SUPPLIED_CERTIFIED,
                    )
                )
            else:
                premises.append(premise)
        assembled = assemble_q_error_budget(premises)
        exact = Q(5, 2) * Q(1, 3)
        self.assertEqual(assembled.exact_total, exact)
        self.assertGreaterEqual(Q(assembled.budget.interpolation), exact)
        self.assertGreaterEqual(Q(assembled.budget.total_upper_bound), exact)
        self.assertGreaterEqual(
            Q(assembled.budget.total_upper_bound),
            Q(assembled.budget.interpolation),
        )
        nearest = float(exact)
        if Q(nearest) < exact:
            self.assertGreater(assembled.budget.interpolation, nearest)


class ActivationPolicyTests(unittest.TestCase):
    def test_canonical_gr0_uses_s_ref_and_forbids_zero_over_zero(self) -> None:
        assessment = assess_activation(
            max_abs_phi=Q(5, 2),
            activation_error=Q(1, 2),
            s_ref=Q(1),
            gr0_phi_policy=GR0_PHI_POLICY_CANONICAL_ZERO,
            matched_controls=_matched_controls(sgbl_phi=Q(1, 4)),
        )
        self.assertEqual(assessment.raw_factor, 2.5)
        self.assertTrue(assessment.threshold_passed)
        self.assertTrue(assessment.control_dominance_passed)
        with self.assertRaisesRegex(Def1Stab1Error, "strictly positive"):
            assess_activation(
                max_abs_phi=Q(0),
                activation_error=Q(0),
                s_ref=Q(0),
                gr0_phi_policy=GR0_PHI_POLICY_CANONICAL_ZERO,
                matched_controls=_matched_controls(),
            )
        with self.assertRaisesRegex(Def1Stab1Error, "nonzero GR-0 numerator"):
            assess_activation(
                max_abs_phi=Q(2),
                activation_error=Q(0),
                s_ref=Q(1),
                gr0_phi_policy=GR0_PHI_POLICY_CANONICAL_ZERO,
                matched_controls=_matched_controls(gr0_phi=Q(1, 131072)),
            )
        with self.assertRaisesRegex(Def1Stab1Error, "historical planted seed"):
            assess_activation(
                max_abs_phi=Q(2),
                activation_error=Q(0),
                s_ref=Q(1),
                gr0_phi_policy=GR0_PHI_POLICY_CANONICAL_ZERO,
                matched_controls=_matched_controls(),
                historical_planted_seed=Q(1, 131072),
            )

    def test_historical_planted_seed_is_a_declared_policy_not_a_silent_default(self) -> None:
        planted = Q(1, 131072)
        assessment = assess_activation(
            max_abs_phi=planted,
            activation_error=Q(0),
            s_ref=2 * planted,
            gr0_phi_policy=GR0_PHI_POLICY_HISTORICAL_PLANTED_SEED,
            matched_controls=_matched_controls(gr0_phi=planted),
            historical_planted_seed=planted,
        )
        self.assertEqual(assessment.raw_factor, 0.5)
        self.assertFalse(assessment.threshold_passed)
        with self.assertRaisesRegex(Def1Stab1Error, "missing"):
            assess_activation(
                max_abs_phi=Q(1),
                activation_error=Q(0),
                s_ref=Q(1),
                gr0_phi_policy=GR0_PHI_POLICY_HISTORICAL_PLANTED_SEED,
                matched_controls=_matched_controls(),
            )
        with self.assertRaisesRegex(Def1Stab1Error, "strictly positive"):
            assess_activation(
                max_abs_phi=Q(1),
                activation_error=Q(0),
                s_ref=Q(1),
                gr0_phi_policy=GR0_PHI_POLICY_HISTORICAL_PLANTED_SEED,
                matched_controls=_matched_controls(),
                historical_planted_seed=Q(0),
            )
        with self.assertRaisesRegex(Def1Stab1Error, "recognized"):
            assess_activation(
                max_abs_phi=Q(1),
                activation_error=Q(0),
                s_ref=Q(1),
                gr0_phi_policy="choose_convenient_gr0_policy",
                matched_controls=_matched_controls(),
            )

    def test_threshold_equality_and_strict_control_dominance(self) -> None:
        self.assertTrue(activation_threshold_passed(MIN_ACTIVATION_LOWER_BOUND))
        self.assertFalse(
            activation_threshold_passed(MIN_ACTIVATION_LOWER_BOUND - Q(1, 2**40))
        )
        self.assertTrue(control_dominance_passed(Q(3), Q(2)))
        self.assertFalse(control_dominance_passed(Q(2), Q(2)))
        equal = assess_activation(
            max_abs_phi=Q(5, 2),
            activation_error=Q(1, 2),
            s_ref=Q(1),
            gr0_phi_policy=GR0_PHI_POLICY_CANONICAL_ZERO,
            matched_controls=_matched_controls(sgbl_phi=Q(2)),
        )
        self.assertTrue(equal.threshold_passed)
        self.assertFalse(equal.control_dominance_passed)
        dominated = assess_activation(
            max_abs_phi=Q(5, 2),
            activation_error=Q(1, 2),
            s_ref=Q(1),
            gr0_phi_policy=GR0_PHI_POLICY_CANONICAL_ZERO,
            matched_controls=_matched_controls(sgbl_phi=Q(1)),
        )
        self.assertTrue(dominated.control_dominance_passed)


class MisnerSharpAlgebraTests(unittest.TestCase):
    def test_minkowski_mass_and_vacuum_defect_vanish(self) -> None:
        geometry = _minkowski_geometry()
        state = geometry.physical_state()
        inverse = inverse_base_metric_from_adm(1, 0, 1)
        gradient = (state.areal_radius.dt, state.areal_radius.dr)
        mass = misner_sharp_mass(state.areal_radius.value, inverse, gradient)
        covariant, inverse_base, scalar = _base_einstein(state)
        mixed = mix_base_tensor(covariant, inverse_base)
        self.assertEqual(inverse, inverse_base)
        self.assertEqual(mass, 0)
        self.assertEqual(scalar, 0)
        self.assertEqual(mixed, ((0, 0), (0, 0)))
        self.assertEqual(
            misner_sharp_mass_gradient(state.areal_radius.value, mixed, gradient),
            (0, 0),
        )
        self.assertEqual(
            misner_sharp_equation_defect(
                state.areal_radius.value, 1, mixed, gradient
            ),
            (0, 0),
        )
        self.assertTrue(
            misner_sharp_identity_holds(
                radius=state.areal_radius.value,
                coupling_F=1,
                mixed_einstein=mixed,
                mixed_effective_source=((0, 0), (0, 0)),
                mixed_equation_residual=mixed,
                radius_derivatives=gradient,
            )
        )

    def test_de_sitter_mass_gradient_is_trace_sensitive(self) -> None:
        geometry = _de_sitter_geometry()
        state = geometry.physical_state()
        inverse = inverse_base_metric_from_adm(1, 0, 1)
        gradient = (state.areal_radius.dt, state.areal_radius.dr)
        mass = misner_sharp_mass(state.areal_radius.value, inverse, gradient)
        covariant, inverse_base, scalar = _base_einstein(state)
        mixed = mix_base_tensor(covariant, inverse_base)
        dm = misner_sharp_mass_gradient(state.areal_radius.value, mixed, gradient)
        cosmological = scalar / 4
        radius = state.areal_radius.value
        self.assertEqual(mass, 4)
        self.assertEqual(cosmological, 3)
        self.assertEqual(mass, cosmological * radius**3 / 6)
        self.assertEqual(dm, (-12, 6))
        self.assertEqual(
            dm,
            (
                (cosmological / 2) * radius**2 * gradient[0],
                (cosmological / 2) * radius**2 * gradient[1],
            ),
        )
        dropped_trace = (
            mixed[0][0] * gradient[0] + mixed[0][1] * gradient[1],
            mixed[1][0] * gradient[0] + mixed[1][1] * gradient[1],
        )
        self.assertNotEqual(
            (radius**2 / 2 * dropped_trace[0], radius**2 / 2 * dropped_trace[1]),
            dm,
        )
        self.assertEqual(
            trace_adjusted_radial_projector(mixed, gradient),
            (dm[0] * 2 / radius**2, dm[1] * 2 / radius**2),
        )
        self.assertTrue(
            misner_sharp_identity_holds(
                radius=radius,
                coupling_F=1,
                mixed_einstein=mixed,
                mixed_effective_source=((0, 0), (0, 0)),
                mixed_equation_residual=mixed,
                radius_derivatives=gradient,
            )
        )

    def test_injected_defect_and_identity_with_compensating_source(self) -> None:
        gradient = (Q(0), Q(1))
        injected = ((Q(0), Q(1)), (Q(0), Q(0)))
        defect = misner_sharp_equation_defect(4, 2, injected, gradient)
        projector = trace_adjusted_radial_projector(injected, gradient)
        self.assertEqual(projector, (Q(1), Q(0)))
        self.assertEqual(defect, (Q(4), Q(0)))
        compensating = ((Q(0), Q(-1)), (Q(0), Q(0)))
        self.assertTrue(
            misner_sharp_identity_holds(
                radius=4,
                coupling_F=2,
                mixed_einstein=((0, 0), (0, 0)),
                mixed_effective_source=compensating,
                mixed_equation_residual=injected,
                radius_derivatives=gradient,
            )
        )
        with self.assertRaisesRegex(Def1Stab1Error, "does not equal F G - T_eff"):
            misner_sharp_identity_holds(
                radius=4,
                coupling_F=2,
                mixed_einstein=((0, 0), (0, 0)),
                mixed_effective_source=((0, 0), (0, 0)),
                mixed_equation_residual=injected,
                radius_derivatives=gradient,
            )

    def test_invalid_f_r_reference_and_hamiltonian_momentum_only(self) -> None:
        inverse = inverse_base_metric_from_adm(1, 0, 1)
        gradient = (Q(0), Q(1))
        mixed = ((0, 0), (0, 0))
        with self.assertRaisesRegex(Def1Stab1Error, "strictly positive"):
            misner_sharp_mass(0, inverse, gradient)
        with self.assertRaisesRegex(Def1Stab1Error, "strictly positive"):
            misner_sharp_equation_defect(4, 0, mixed, gradient)
        with self.assertRaisesRegex(Def1Stab1Error, "strictly positive"):
            misner_sharp_equation_defect(4, -1, mixed, gradient)
        with self.assertRaisesRegex(Def1Stab1Error, "lapse"):
            inverse_base_metric_from_adm(0, 0, 1)
        with self.assertRaisesRegex(Def1Stab1Error, "Hamiltonian/momentum"):
            assess_mass_flux_ledger(
                radius=4,
                inverse_metric=inverse,
                radius_derivatives=gradient,
                coupling_F=1,
                mixed_equation_residual=mixed,
                delta_mass=0,
                integrated_flux=0,
                residual_enclosure=0,
                hamiltonian_momentum_only=True,
            )
        with self.assertRaisesRegex(Def1Stab1Error, "2-by-2"):
            assess_mass_flux_ledger(
                radius=4,
                inverse_metric=inverse,
                radius_derivatives=gradient,
                coupling_F=1,
                mixed_equation_residual=((0, 0),),
                delta_mass=0,
                integrated_flux=0,
                residual_enclosure=0,
            )
        with self.assertRaises(TypeError):
            misner_sharp_mass(4, inverse, (True, 1))

    def test_mass_flux_veto_maps_to_existing_protocol_token(self) -> None:
        inverse = inverse_base_metric_from_adm(1, 0, 1)
        gradient = (Q(0), Q(1))
        mixed = ((0, 0), (0, 0))
        valid = assess_mass_flux_ledger(
            radius=4,
            inverse_metric=inverse,
            radius_derivatives=gradient,
            coupling_F=1,
            mixed_equation_residual=mixed,
            mixed_einstein=mixed,
            delta_mass=Q(1, 2),
            integrated_flux=Q(1, 2),
            residual_enclosure=0,
        )
        self.assertTrue(valid.ledger_valid)
        self.assertIsNone(valid.protocol_token)
        self.assertEqual(valid.mass, 0)
        equal_edge = assess_mass_flux_ledger(
            radius=4,
            inverse_metric=inverse,
            radius_derivatives=gradient,
            coupling_F=1,
            mixed_equation_residual=mixed,
            delta_mass=Q(1, 2),
            integrated_flux=0,
            residual_enclosure=Q(1, 2),
        )
        self.assertTrue(equal_edge.ledger_valid)
        vetoed = assess_mass_flux_ledger(
            radius=4,
            inverse_metric=inverse,
            radius_derivatives=gradient,
            coupling_F=1,
            mixed_equation_residual=mixed,
            delta_mass=Q(1, 2),
            integrated_flux=0,
            residual_enclosure=Q(1, 4),
        )
        self.assertFalse(vetoed.ledger_valid)
        self.assertEqual(vetoed.protocol_token, PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN)
        self.assertEqual(vetoed.veto_class, STOPPED_MASS_FLUX_INCONSISTENCY)
        self.assertEqual(
            map_mass_flux_inconsistency_token("hamiltonian_momentum_only"),
            PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN,
        )
        self.assertEqual(
            map_mass_flux_inconsistency_token("missing_trace_adjusted_projector"),
            PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN,
        )
        with self.assertRaisesRegex(Def1Stab1Error, "outside the typed"):
            map_mass_flux_inconsistency_token("new_mass_token")
        with self.assertRaisesRegex(Def1Stab1Error, "existing protocol token"):
            MassFluxAssessment(
                mass=0,
                mass_gradient=None,
                equation_defect=None,
                residual=1,
                enclosure=0,
                ledger_valid=False,
                protocol_token="stopped_mass_flux_inconsistency",
                veto_class=STOPPED_MASS_FLUX_INCONSISTENCY,
                typed_reason="conservation_residual_exceeds_enclosure",
                used_hamiltonian_momentum_only=False,
            )


class MarginHelperTests(unittest.TestCase):
    def test_trappedness_strict_and_one_bit_edges(self) -> None:
        error = Q(1, 8)
        threshold = -4 * error
        self.assertFalse(
            trappedness_margin_passed(threshold, threshold, error, error)
        )
        below = threshold - Q(1, 2**40)
        above = threshold + Q(1, 2**40)
        self.assertTrue(trappedness_margin_passed(below, below, error, error))
        self.assertFalse(trappedness_margin_passed(above, below, error, error))
        float_threshold = float(threshold)
        ulp_below = Q(nextafter(float_threshold, -inf))
        ulp_above = Q(nextafter(float_threshold, inf))
        self.assertTrue(
            trappedness_margin_passed(ulp_below, ulp_below, error, error)
        )
        self.assertFalse(
            trappedness_margin_passed(ulp_above, ulp_below, error, error)
        )

    def test_complete_q_requires_strict_positive_error_margin(self) -> None:
        error = Q(1, 16)
        self.assertFalse(complete_q_margin_passed(4 * error, error))
        self.assertTrue(complete_q_margin_passed(4 * error + Q(1, 2**40), error))
        self.assertFalse(complete_q_margin_passed(Q(1), 0))
        with self.assertRaisesRegex(Def1Stab1Error, "nonnegative"):
            complete_q_margin_passed(Q(1), Q(-1, 8))

    def test_affine_interval_samples_and_overlapping_refuse(self) -> None:
        scale = Q(1)
        passing = [
            Q(index, MIN_AFFINE_SAMPLES - 1) * MIN_AFFINE_INTERVAL_OVER_L0
            for index in range(MIN_AFFINE_SAMPLES)
        ]
        assessment = affine_interval_and_sample_gate(passing, scale)
        self.assertTrue(assessment.passed)
        self.assertEqual(assessment.interval_over_L0, MIN_AFFINE_INTERVAL_OVER_L0)
        seven = affine_interval_and_sample_gate(passing[:-1], scale)
        self.assertFalse(seven.passed)
        self.assertFalse(seven.samples_passed)
        short = [
            Q(index, MIN_AFFINE_SAMPLES - 1) * Q(1, 128)
            for index in range(MIN_AFFINE_SAMPLES)
        ]
        short_assessment = affine_interval_and_sample_gate(short, scale)
        self.assertTrue(short_assessment.samples_passed)
        self.assertFalse(short_assessment.interval_passed)
        with self.assertRaisesRegex(Def1Stab1Error, "overlapping"):
            affine_interval_and_sample_gate([0, 0, 1, 2, 3, 4, 5, 6], 1)
        with self.assertRaisesRegex(Def1Stab1Error, "strictly positive"):
            affine_interval_and_sample_gate(passing, 0)

    def test_route_affine_resolution_and_method_edges(self) -> None:
        self.assertTrue(direct_raychaudhuri_routes_agree(0, DIRECT_RAYCHAUDHURI_AGREEMENT_MAX))
        self.assertFalse(
            direct_raychaudhuri_routes_agree(
                0, DIRECT_RAYCHAUDHURI_AGREEMENT_MAX + Q(1, 10**12)
            )
        )
        self.assertTrue(
            affine_normalization_residual_passed(AFFINE_NORMALIZATION_RESIDUAL_MAX)
        )
        self.assertFalse(
            affine_normalization_residual_passed(
                AFFINE_NORMALIZATION_RESIDUAL_MAX + Q(1, 10**12)
            )
        )
        self.assertTrue(declared_resolution_count_meets_minimum(3))
        self.assertFalse(declared_resolution_count_meets_minimum(2))
        self.assertTrue(declared_method_count_meets_minimum(2))
        self.assertFalse(declared_method_count_meets_minimum(1))
        self.assertTrue(observed_order_meets_minimum(Q(3, 2)))
        self.assertFalse(observed_order_meets_minimum(Q(3, 2) - Q(1, 100)))


class IndependentBooleanTests(unittest.TestCase):
    def test_nine_booleans_remain_independent_and_complete(self) -> None:
        self.assertEqual(len(DEF1_BOOLEAN_NAMES), 9)
        record = record_def1_booleans(
            resolved_activation=True,
            control_dominance=True,
            resolved_trapped_interval=True,
            resolved_complete_defocusing=True,
            direct_Raychaudhuri_agreement=True,
            all_health_constraints_scales_valid=False,
            mass_flux_ledger_valid=True,
            three_resolution_convergence=False,
            two_method_agreement=False,
        )
        mapping = record.as_mapping()
        self.assertTrue(mapping["resolved_complete_defocusing"])
        self.assertFalse(mapping["all_health_constraints_scales_valid"])
        self.assertFalse(mapping["three_resolution_convergence"])
        self.assertFalse(mapping["two_method_agreement"])
        with self.assertRaisesRegex(Def1Stab1Error, "cannot be inferred"):
            record_def1_booleans_from_mapping(
                {
                    "resolved_activation": True,
                    "control_dominance": True,
                    "resolved_trapped_interval": True,
                    "resolved_complete_defocusing": True,
                    "direct_Raychaudhuri_agreement": True,
                    "mass_flux_ledger_valid": True,
                }
            )
        with self.assertRaises(TypeError):
            record_def1_booleans(
                resolved_activation=True,
                control_dominance=True,
                resolved_trapped_interval=True,
                resolved_complete_defocusing=True,
                direct_Raychaudhuri_agreement=True,
                mass_flux_ledger_valid=True,
                three_resolution_convergence=True,
                two_method_agreement=True,
            )
        self.assertTrue(declared_resolution_count_meets_minimum(3))
        self.assertTrue(declared_method_count_meets_minimum(2))
        self.assertIsInstance(record, Def1BooleanRecord)


class IndependentGeometricControlTests(unittest.TestCase):
    def test_minkowski_de_sitter_and_synthetic_q_remain_independent(self) -> None:
        minkowski = metric_null_raychaudhuri_point_certificate(
            _minkowski_geometry(),
            branch="outgoing",
            affine_scale=Jet2.constant(1),
        )
        de_sitter = metric_null_raychaudhuri_point_certificate(
            _de_sitter_geometry(),
            branch="outgoing",
            affine_scale=Jet2(1, dt=1, dtt=1),
        )
        synthetic = metric_null_raychaudhuri_point_certificate(
            _synthetic_positive_q_geometry(),
            branch="outgoing",
            affine_scale=Jet2.constant(1),
        )
        self.assertEqual(minkowski["null_frame"]["round_sphere_classification"], "normal")
        self.assertFalse(
            trappedness_margin_passed(
                minkowski["null_frame"]["theta_plus"],
                minkowski["null_frame"]["theta_minus"],
                0,
                0,
            )
        )
        self.assertEqual(minkowski["raychaudhuri"]["complete_rhs"], Q(-1, 8))
        self.assertTrue(
            trappedness_margin_passed(
                de_sitter["null_frame"]["theta_plus"],
                de_sitter["null_frame"]["theta_minus"],
                0,
                0,
            )
        )
        self.assertFalse(
            complete_q_margin_passed(de_sitter["raychaudhuri"]["complete_rhs"], Q(1, 16))
        )
        self.assertEqual(de_sitter["raychaudhuri"]["minus_R_ab_k_a_k_b"], 0)
        self.assertEqual(synthetic["null_frame"]["round_sphere_classification"], "marginal")
        self.assertFalse(
            trappedness_margin_passed(
                synthetic["null_frame"]["theta_plus"],
                synthetic["null_frame"]["theta_minus"],
                0,
                0,
            )
        )
        self.assertTrue(
            complete_q_margin_passed(synthetic["raychaudhuri"]["complete_rhs"], Q(1, 16))
        )
        booleans = record_def1_booleans(
            resolved_activation=False,
            control_dominance=False,
            resolved_trapped_interval=trappedness_margin_passed(
                de_sitter["null_frame"]["theta_plus"],
                de_sitter["null_frame"]["theta_minus"],
                0,
                0,
            ),
            resolved_complete_defocusing=complete_q_margin_passed(
                synthetic["raychaudhuri"]["complete_rhs"], Q(1, 16)
            ),
            direct_Raychaudhuri_agreement=True,
            all_health_constraints_scales_valid=False,
            mass_flux_ledger_valid=True,
            three_resolution_convergence=False,
            two_method_agreement=False,
        )
        self.assertTrue(booleans.resolved_trapped_interval)
        self.assertTrue(booleans.resolved_complete_defocusing)
        self.assertFalse(booleans.all_health_constraints_scales_valid)
        self.assertFalse(booleans.three_resolution_convergence)


def _uniform_budget(component: float, total: float) -> RaychaudhuriErrorBudget:
    return RaychaudhuriErrorBudget(
        **{name: component for name in ERROR_BUDGET_COMPONENTS},
        total_upper_bound=total,
    )


class ConstructorContractHoleTests(unittest.TestCase):
    def test_float_dtos_normalize_exact_integers_and_refuse_silent_rounding(self) -> None:
        budget = _uniform_budget(1, len(ERROR_BUDGET_COMPONENTS))
        self.assertTrue(all(type(getattr(budget, name)) is float
                            for name in (*ERROR_BUDGET_COMPONENTS, "total_upper_bound")))
        components = {name: 0 for name in ERROR_BUDGET_COMPONENTS}
        components["arithmetic"] = 2**53 + 1
        with self.assertRaisesRegex(Def1Stab1Error, "losslessly"):
            RaychaudhuriErrorBudget(**components, total_upper_bound=float(2**53))
        with self.assertRaisesRegex(Def1Stab1Error, "lower bound cannot exceed"):
            ActivationAssessment(raw_factor=0, lower_bound=2,
                                 matched_control_upper_bound=0,
                                 threshold_passed=True, control_dominance_passed=True)
        with self.assertRaisesRegex(Def1Stab1Error, "losslessly"):
            ActivationAssessment(raw_factor=float(2**54), lower_bound=2**53 + 1,
                                 matched_control_upper_bound=0,
                                 threshold_passed=True, control_dominance_passed=True)

    def test_budget_total_must_dominate_stored_binary64_sum(self) -> None:
        with self.assertRaisesRegex(Def1Stab1Error, "exact sum of the stored"):
            _uniform_budget(1.0, 1.0)
        tiny = nextafter(0.0, inf)
        stored_sum = Q(0)
        for _ in ERROR_BUDGET_COMPONENTS:
            stored_sum += Q(tiny)
        total = float(stored_sum)
        if Q(total) < stored_sum:
            total = nextafter(total, inf)
        budget = _uniform_budget(tiny, total)
        self.assertGreaterEqual(Q(budget.total_upper_bound), stored_sum)
        with self.assertRaisesRegex(Def1Stab1Error, "exact sum of the stored"):
            _uniform_budget(tiny, tiny)
        huge = 1e308
        with self.assertRaisesRegex(Def1Stab1Error, "exact sum of the stored"):
            _uniform_budget(huge, huge)
        premises = []
        for premise in zero_error_premises():
            if premise.component == "arithmetic":
                premises.append(
                    replace(
                        premise,
                        sensitivity=Interval.singleton(10**400),
                        input_error=Interval.singleton(1),
                    )
                )
            else:
                premises.append(premise)
        with self.assertRaisesRegex(Def1Stab1Error, "finite range"):
            assemble_q_error_budget(premises)

    def test_positive_definite_singular_and_nonsymmetric_inverses_refuse(self) -> None:
        euclidean = ((1, 0), (0, 1))
        singular = ((-1, 0), (0, 0))
        nonsymmetric = ((-1, 1), (0, 1))
        lorentz = ((-1, 0), (0, 1))
        covariant = ((0, 1), (0, 0))
        gradient = (1, 0)
        for matrix in (euclidean, singular, nonsymmetric):
            with self.assertRaisesRegex(Def1Stab1Error, "Lorentzian|symmetric"):
                misner_sharp_mass(2, matrix, gradient)
            with self.assertRaisesRegex(Def1Stab1Error, "Lorentzian|symmetric"):
                mix_base_tensor(((0, 0), (0, 0)), matrix)
            with self.assertRaisesRegex(Def1Stab1Error, "Lorentzian|symmetric"):
                inverse_base_metric(matrix)
            with self.assertRaisesRegex(Def1Stab1Error, "Lorentzian|symmetric"):
                assess_mass_flux_ledger(
                    radius=2,
                    inverse_metric=matrix,
                    radius_derivatives=gradient,
                    coupling_F=1,
                    mixed_equation_residual=((0, 0), (0, 0)),
                    delta_mass=0,
                    integrated_flux=0,
                    residual_enclosure=0,
                )
        self.assertEqual(misner_sharp_mass(2, lorentz, gradient), 2)
        mixed = mix_base_tensor(covariant, lorentz)
        self.assertEqual(mixed[0][1], Q(1))

    def test_mass_flux_replace_residual_attack_is_rejected(self) -> None:
        inverse = inverse_base_metric_from_adm(1, 0, 1)
        valid = assess_mass_flux_ledger(
            radius=4,
            inverse_metric=inverse,
            radius_derivatives=(0, 1),
            coupling_F=1,
            mixed_equation_residual=((0, 0), (0, 0)),
            delta_mass=Q(1),
            integrated_flux=0,
            residual_enclosure=Q(1),
        )
        self.assertTrue(valid.ledger_valid)
        with self.assertRaisesRegex(Def1Stab1Error, "residual <= enclosure"):
            replace(valid, residual=2)
        vetoed = replace(
            valid,
            residual=2,
            ledger_valid=False,
            protocol_token=PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN,
            veto_class=STOPPED_MASS_FLUX_INCONSISTENCY,
            typed_reason="conservation_residual_exceeds_enclosure",
        )
        self.assertFalse(vetoed.ledger_valid)

    def test_affine_interval_direct_constructor_rejects_inconsistent_flags(self) -> None:
        with self.assertRaisesRegex(Def1Stab1Error, "samples_passed"):
            AffineIntervalAssessment(
                sample_count=1,
                interval_over_L0=0,
                samples_passed=True,
                interval_passed=True,
                passed=True,
            )
        failing = AffineIntervalAssessment(
            sample_count=1,
            interval_over_L0=0,
            samples_passed=False,
            interval_passed=False,
            passed=False,
        )
        self.assertFalse(failing.passed)
        passing = affine_interval_and_sample_gate(
            [Q(index, 7) / 64 for index in range(8)],
            1,
        )
        with self.assertRaisesRegex(Def1Stab1Error, "samples_passed"):
            replace(passing, sample_count=1)

    def test_activation_saved_bounds_reject_unsupported_true_flags(self) -> None:
        with self.assertRaisesRegex(Def1Stab1Error, "threshold_passed"):
            ActivationAssessment(
                raw_factor=0,
                lower_bound=0,
                matched_control_upper_bound=1,
                threshold_passed=True,
                control_dominance_passed=True,
            )
        with self.assertRaises(TypeError):
            ActivationAssessment(
                raw_factor=0,
                lower_bound=0,
                matched_control_upper_bound=1,
                threshold_passed=1,
                control_dominance_passed=False,
            )
        consistent = ActivationAssessment(
            raw_factor=0,
            lower_bound=0,
            matched_control_upper_bound=1,
            threshold_passed=False,
            control_dominance_passed=False,
        )
        self.assertFalse(consistent.threshold_passed)
        self.assertFalse(consistent.control_dominance_passed)
        near = assess_activation(
            max_abs_phi=Q(5, 2),
            activation_error=Q(1, 2),
            s_ref=Q(1),
            gr0_phi_policy=GR0_PHI_POLICY_CANONICAL_ZERO,
            matched_controls=_matched_controls(sgbl_phi=Q(2) - Q(1, 2**60)),
        )
        self.assertTrue(near.threshold_passed)
        self.assertFalse(near.control_dominance_passed)
        self.assertFalse(
            Q(near.lower_bound) > Q(near.matched_control_upper_bound)
        )

    def test_assembled_record_rejects_fake_flags_drift_and_malformed_inventories(
        self,
    ) -> None:
        assembled = assemble_q_error_budget(zero_error_premises())
        with self.assertRaisesRegex(Def1Stab1Error, "status declarations"):
            replace(assembled, all_components_declared_proven_enclosures=False)
        with self.assertRaisesRegex(Def1Stab1Error, "not global PDE error"):
            replace(assembled, global_pde_error_certified=True)
        with self.assertRaises(TypeError):
            replace(assembled, used_measured_q=1)
        with self.assertRaisesRegex(Def1Stab1Error, "exact_total"):
            replace(assembled, exact_total=Q(1))
        inflated = assemble_q_error_budget(
            inject_component_input_error(
                zero_error_premises(), "interpolation", Q(1, 3)
            )
        )
        with self.assertRaisesRegex(Def1Stab1Error, "drifted below"):
            replace(inflated, budget=assembled.budget)
        reversed_status = tuple(reversed(assembled.component_statuses))
        with self.assertRaisesRegex(Def1Stab1Error, "frozen component order"):
            replace(assembled, component_statuses=reversed_status)
        bad_label = (
            ("spatial_temporal", "proven"),
            *assembled.component_statuses[1:],
        )
        with self.assertRaisesRegex(Def1Stab1Error, "malformed declaration"):
            replace(assembled, component_statuses=bad_label)


if __name__ == "__main__":
    unittest.main()
