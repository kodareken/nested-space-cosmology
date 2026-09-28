from __future__ import annotations

from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.fgc.def1_geometry_error import (  # noqa: E402
    INPUT_NAMES,
    evaluate_q_box,
    q_input_error_bound,
    q_sensitivity_enclosures,
)
from recursive_horizons.fgc.def1_stab1 import (  # noqa: E402
    COMPLETE_Q_ERROR_FACTOR,
    DEF1_BOOLEAN_NAMES,
    DIRECT_RAYCHAUDHURI_AGREEMENT_MAX,
    ERROR_BUDGET_COMPONENTS,
    MIN_ACTIVATION_LOWER_BOUND,
    MIN_AFFINE_INTERVAL_OVER_L0,
    MIN_AFFINE_SAMPLES,
    MIN_INDEPENDENT_METHODS,
    MIN_NESTED_RESOLUTIONS,
    MIN_OBSERVED_ORDER,
    PREMISE_STATUS_CONDITIONAL,
    PREMISE_STATUS_PROVEN,
    PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN,
    SOURCE_EXACT_ALGEBRA,
    SOURCE_IMP1_ADMISSION_DEBIT,
    SOURCE_RICHARDSON,
    SOURCE_SUPPLIED_CERTIFIED,
    TRAPPEDNESS_ERROR_FACTOR,
    QComponentPremise,
    assemble_q_error_budget,
    assess_mass_flux_ledger,
    inverse_base_metric_from_adm,
    record_def1_booleans,
)
from recursive_horizons.fgc.def1_stab1_providers import (  # noqa: E402
    CONVERSION_GAUGE_EXTENSION,
    CONVERSION_IMP1,
    CONVERSION_PHYSICAL_RESIDUAL,
    CONVERSION_INTERPOLATION_SECOND_ORDER,
    CONVERSION_RICHARDSON,
    CONVERSION_SUPPLIED_RADII,
    FULL_RESIDUAL_LENGTH,
    METRIC_ACCELERATION_INPUTS,
    METRIC_FIRST_DERIVATIVE_INPUTS,
    REMAINING_GATE_WORK,
    AssembledProviderMap,
    ComponentContribution,
    Def1Stab1ProviderError,
    ProviderRecord,
    affine_provider,
    affine_residual_times_step_refused,
    arithmetic_provider,
    assemble_provider_error_map,
    boundary_provider,
    complete_zero_radii,
    conservation_provider,
    conservative_residual_q_debit,
    extraction_provider,
    gauge_constraint_provider,
    initial_data_provider,
    interpolation_provider,
    joint_input_radii,
    nonlinear_source_provider,
    nullness_provider,
    physical_constraint_provider,
    reduction_constraint_provider,
    residual_kk,
    residual_q_debit,
    spatial_temporal_provider,
    trajectory_alignment_provider,
)
from recursive_horizons.fgc.exact_interval import Interval  # noqa: E402


def _minkowski() -> dict[str, Q]:
    inputs = {name: Q(0) for name in INPUT_NAMES}
    inputs.update(
        {
            "alpha.value": Q(1),
            "lambda.value": Q(1),
            "R.value": Q(2),
            "R.dr": Q(1),
            "k.t": Q(1),
            "k.r": Q(1),
        }
    )
    return inputs


def _zero_tensor() -> tuple[tuple[Q, ...], ...]:
    zero = (Q(0), Q(0), Q(0), Q(0))
    return (zero, zero, zero, zero)


SHARED_EVALUATION_CONTEXT = (
    "synthetic joint-box evaluation at declared nominal ADM two-jet"
)
SHARED_UNIT = "geometry-input units and complete-Q units"


def _decl(component: str, extra: str = "synthetic injected control") -> dict[str, str]:
    return {
        "context": SHARED_EVALUATION_CONTEXT,
        "unit": SHARED_UNIT,
        "provenance": f"{component}: {extra}",
    }


def _radii_with(**slots: object) -> dict[str, Q]:
    values = complete_zero_radii()
    for name, raw in slots.items():
        values[name] = Q(raw) if not isinstance(raw, Q) else raw
    return values


def _mass_flux(*, valid: bool = True):
    inverse = inverse_base_metric_from_adm(1, 0, 1)
    if valid:
        return assess_mass_flux_ledger(
            radius=4,
            inverse_metric=inverse,
            radius_derivatives=(Q(0), Q(1)),
            coupling_F=1,
            mixed_equation_residual=((0, 0), (0, 0)),
            delta_mass=0,
            integrated_flux=0,
            residual_enclosure=0,
        )
    return assess_mass_flux_ledger(
        radius=4,
        inverse_metric=inverse,
        radius_derivatives=(Q(0), Q(1)),
        coupling_F=1,
        mixed_equation_residual=((0, 0), (0, 0)),
        delta_mass=Q(1, 2),
        integrated_flux=0,
        residual_enclosure=Q(1, 4),
    )


def _zero_reduction() -> dict[str, Q]:
    return {name: Q(0) for name in METRIC_FIRST_DERIVATIVE_INPUTS}


def _zero_interpolation() -> dict[str, Q]:
    return {name: Q(0) for name in INPUT_NAMES}


def _exact_zero_spatial() -> ProviderRecord:
    return spatial_temporal_provider(
        source=SOURCE_EXACT_ALGEBRA,
        status=PREMISE_STATUS_PROVEN,
        radii=complete_zero_radii(),
        additive_q=0,
        **_decl("spatial_temporal", "explicit exact-algebra zero"),
    )


def _physical(residual=None, **kwargs) -> ProviderRecord:
    return physical_constraint_provider(
        residual_E=residual if residual is not None else _zero_tensor(),
        tangent=(Q(1), Q(1)),
        coupling_F=1,
        status=PREMISE_STATUS_PROVEN,
        source=SOURCE_EXACT_ALGEBRA,
        **_decl("physical_constraint"),
        **kwargs,
    )


def _assembly_records(
    *,
    initial: dict[str, Q] | None = None,
    alignment: dict[str, Q] | None = None,
) -> tuple[ProviderRecord, ...]:
    flux = _mass_flux()
    return (
        _exact_zero_spatial(),
        _physical(),
        gauge_constraint_provider(
            extension_E=_zero_tensor(),
            tangent=(Q(1), Q(1)),
            coupling_F=1,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl("gauge_constraint"),
        ),
        reduction_constraint_provider(
            first_derivative_discrepancies=_zero_reduction(),
            derivative_operator_bound=1,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl("reduction_constraint"),
        ),
        initial_data_provider(
            radii=initial or complete_zero_radii(),
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl("initial_data"),
        ),
        nonlinear_source_provider(
            full_residual=(0,) * FULL_RESIDUAL_LENGTH,
            inverse_jacobian_bound=1,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl("nonlinear_source"),
        ),
        affine_provider(
            transport_defect=0,
            affine_interval=Q(1, 8),
            gronwall_factor=1,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl("affine"),
        ),
        nullness_provider(
            null_residual=0,
            null_to_tangent_bound=1,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl("nullness"),
        ),
        trajectory_alignment_provider(
            radii=alignment or complete_zero_radii(),
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl("trajectory_alignment"),
        ),
        interpolation_provider(
            sample_spacing=Q(1, 8),
            first_derivative_enclosures=_zero_interpolation(),
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl("interpolation"),
        ),
        extraction_provider(
            affine_samples=(0, Q(1, 8), Q(1, 4)),
            derivative_bound=0,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl("extraction"),
        ),
        boundary_provider(
            additional_debit=0,
            no_influence_premise=True,
            coverage_premise=True,
            independent_guard=True,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl("boundary"),
        ),
        conservation_provider(
            additional_debit=0,
            mass_flux=flux,
            coverage_premise=True,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl("conservation"),
        ),
        arithmetic_provider(
            additive_q=0,
            source=SOURCE_EXACT_ALGEBRA,
            status=PREMISE_STATUS_PROVEN,
            **_decl("arithmetic"),
        ),
    )


def _nonzero_records() -> dict[str, ProviderRecord]:
    injected = [list(row) for row in _zero_tensor()]
    injected[0][0] = Q(1)
    residual_E = tuple(tuple(row) for row in injected)
    return {
        "spatial_temporal": spatial_temporal_provider(
            source=SOURCE_RICHARDSON,
            status=PREMISE_STATUS_CONDITIONAL,
            additive_q=Q(1, 16),
            **_decl("spatial_temporal", "conditional Richardson debit"),
        ),
        "physical_constraint": physical_constraint_provider(
            residual_E=residual_E,
            tangent=(Q(1), Q(1)),
            coupling_F=2,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl("physical_constraint"),
        ),
        "gauge_constraint": gauge_constraint_provider(
            extension_E=residual_E,
            tangent=(Q(1), Q(1)),
            coupling_F=2,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl("gauge_constraint"),
        ),
        "reduction_constraint": reduction_constraint_provider(
            first_derivative_discrepancies=_zero_reduction() | {"alpha.dt": Q(1, 5)},
            derivative_operator_bound=3,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl("reduction_constraint"),
        ),
        "initial_data": initial_data_provider(
            radii=_radii_with(**{"R.dtt": Q(1, 32)}),
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl("initial_data"),
        ),
        "nonlinear_source": nonlinear_source_provider(
            full_residual=(Q(1, 8), 0, 0, 0, 0, 0),
            inverse_jacobian_bound=2,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl("nonlinear_source"),
        ),
        "affine": affine_provider(
            transport_defect=Q(1, 10),
            affine_interval=Q(1, 8),
            gronwall_factor=3,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl("affine"),
        ),
        "nullness": nullness_provider(
            null_residual=Q(1, 6),
            null_to_tangent_bound=2,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl("nullness"),
        ),
        "trajectory_alignment": trajectory_alignment_provider(
            radii=_radii_with(**{"k.r": Q(1, 20)}),
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl("trajectory_alignment"),
        ),
        "interpolation": interpolation_provider(
            sample_spacing=Q(1, 8),
            first_derivative_enclosures=_zero_interpolation() | {"R.value": Q(2)},
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl("interpolation"),
        ),
        "extraction": extraction_provider(
            affine_samples=(0, Q(1, 2), 1),
            derivative_bound=Q(1, 3),
            sample_q_values=(1, 1, 1),
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl("extraction"),
        ),
        "boundary": boundary_provider(
            additional_debit=Q(1, 7),
            physical_causality_passed=True,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl("boundary"),
        ),
        "conservation": conservation_provider(
            additional_debit=Q(1, 9),
            mass_flux=_mass_flux(),
            coverage_premise=True,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl("conservation"),
        ),
        "arithmetic": arithmetic_provider(
            additive_q=Q(1, 64),
            source=SOURCE_IMP1_ADMISSION_DEBIT,
            status=PREMISE_STATUS_CONDITIONAL,
            **_decl("arithmetic", "conditional IMP1 local debit"),
        ),
    }


class ProviderContractTests(unittest.TestCase):
    def test_fourteen_components_and_declaration_fields(self) -> None:
        records = _assembly_records()
        self.assertEqual(tuple(record.component for record in records), ERROR_BUDGET_COMPONENTS)
        self.assertEqual(len(ERROR_BUDGET_COMPONENTS), 14)
        assembled = assemble_provider_error_map(
            _minkowski(), records, mass_flux=_mass_flux()
        )
        self.assertEqual(
            tuple(item.component for item in assembled.contributions),
            ERROR_BUDGET_COMPONENTS,
        )
        for record in assembled.records:
            self.assertEqual(tuple(name for name, _ in record.radii), INPUT_NAMES)
            self.assertEqual(record.context, SHARED_EVALUATION_CONTEXT)
            self.assertEqual(record.unit, SHARED_UNIT)
            self.assertTrue(record.provenance)
            self.assertTrue(record.conversion)
        self.assertEqual(assembled.evaluation_context, SHARED_EVALUATION_CONTEXT)
        self.assertEqual(assembled.unit, SHARED_UNIT)
        self.assertEqual(
            tuple(name for name, _ in assembled.nominal), INPUT_NAMES
        )
        self.assertFalse(assembled.def1_error_map_passed)
        self.assertFalse(assembled.used_measured_q)
        self.assertFalse(assembled.global_pde_error_certified)
        self.assertFalse(assembled.richardson_treated_as_global_pde_error)
        self.assertFalse(assembled.imp1_admission_debit_treated_as_global_pde_error)
        self.assertEqual(assembled.remaining_gate_work, REMAINING_GATE_WORK)
        self.assertIsNone(assembled.def1_booleans)

    def test_joint_box_covers_nonlinear_truth_that_independent_boxes_underbound(
        self,
    ) -> None:
        nominal = _minkowski()
        initial = _radii_with(**{"R.dr": Q(1, 8)})
        alignment = _radii_with(**{"k.r": Q(1, 4)})
        assembled = assemble_provider_error_map(
            nominal,
            _assembly_records(initial=initial, alignment=alignment),
            mass_flux=_mass_flux(),
        )
        joint_radii = {name: value for name, value in assembled.joint_radii}
        self.assertEqual(joint_radii["R.dr"], Q(1, 8))
        self.assertEqual(joint_radii["k.r"], Q(1, 4))
        self.assertEqual(joint_radii, joint_input_radii(assembled.records))
        independent = q_input_error_bound(nominal, initial) + q_input_error_bound(
            nominal, alignment
        )
        truth = nominal | {"R.dr": Q(9, 8), "k.r": Q(5, 4)}
        observed = abs(evaluate_q_box(truth)["q"].lower - evaluate_q_box(nominal)["q"].lower)
        joint_total = assembled.assembled.exact_total
        self.assertEqual(independent, Q(29, 64))
        self.assertEqual(observed, Q(1001, 2048))
        self.assertEqual(joint_total, Q(315, 512))
        self.assertGreater(observed, independent)
        self.assertLessEqual(observed, joint_total)
        own_box = {
            name: Interval(nominal[name] - initial[name], nominal[name] + initial[name])
            for name in INPUT_NAMES
        }
        joint_box = {name: interval for name, interval in assembled.joint_box}
        self.assertGreater(
            q_sensitivity_enclosures(joint_box)["R.dr"].abs_upper(),
            q_sensitivity_enclosures(own_box)["R.dr"].abs_upper(),
        )
        contributions = {
            item.component: item.total for item in assembled.contributions
        }
        self.assertEqual(
            contributions["initial_data"],
            assembled.lipschitz[INPUT_NAMES.index("R.dr")][1].abs_upper() * Q(1, 8),
        )
        self.assertEqual(
            contributions["trajectory_alignment"],
            assembled.lipschitz[INPUT_NAMES.index("k.r")][1].abs_upper() * Q(1, 4),
        )

    def test_each_provider_has_a_supported_nonzero_control(self) -> None:
        records = _nonzero_records()
        self.assertEqual(tuple(records), ERROR_BUDGET_COMPONENTS)
        self.assertEqual(records["physical_constraint"].additive_q, Q(1, 2))
        self.assertEqual(records["gauge_constraint"].additive_q, Q(1, 2))
        self.assertEqual(records["spatial_temporal"].additive_q, Q(1, 16))
        self.assertEqual(records["spatial_temporal"].conversion, CONVERSION_RICHARDSON)
        self.assertEqual(records["spatial_temporal"].status, PREMISE_STATUS_CONDITIONAL)
        self.assertEqual(
            records["reduction_constraint"].radii_mapping["alpha.dt"], Q(3, 5)
        )
        self.assertEqual(records["initial_data"].radii_mapping["R.dtt"], Q(1, 32))
        for name in METRIC_ACCELERATION_INPUTS:
            self.assertEqual(records["nonlinear_source"].radii_mapping[name], Q(1, 4))
        self.assertEqual(records["affine"].radii_mapping["k.t"], Q(3, 10))
        self.assertEqual(records["nullness"].radii_mapping["k.r"], Q(1, 3))
        self.assertEqual(records["trajectory_alignment"].radii_mapping["k.r"], Q(1, 20))
        self.assertEqual(records["interpolation"].radii_mapping["R.value"], Q(1, 4))
        self.assertEqual(records["extraction"].additive_q, Q(1, 6))
        self.assertEqual(records["boundary"].additive_q, Q(1, 7))
        self.assertEqual(records["conservation"].additive_q, Q(1, 9))
        self.assertEqual(records["arithmetic"].additive_q, Q(1, 64))
        self.assertEqual(records["arithmetic"].conversion, CONVERSION_IMP1)
        for name, record in records.items():
            with self.subTest(name=name):
                self.assertFalse(record.zero_contribution)
                self.assertEqual(record.component, name)


class ProviderFailureTests(unittest.TestCase):
    def test_complete_residual_inputs_cannot_be_omitted_or_relabelled(self) -> None:
        common = {
            "radii": tuple((name, Q(0)) for name in INPUT_NAMES),
            "additive_q": Q(0),
            "status": PREMISE_STATUS_PROVEN,
            "source": SOURCE_EXACT_ALGEBRA,
        }
        for component, conversion in (
            ("physical_constraint", CONVERSION_PHYSICAL_RESIDUAL),
            ("gauge_constraint", CONVERSION_GAUGE_EXTENSION),
        ):
            with self.subTest(component=component), self.assertRaisesRegex(
                Def1Stab1ProviderError, "omitted its bound inputs"
            ):
                ProviderRecord(
                    component=component,
                    conversion=conversion,
                    **common,
                    **_decl(component),
                )
        initial = _nonzero_records()["initial_data"]
        with self.assertRaisesRegex(Def1Stab1ProviderError, "cannot relabel"):
            replace(
                initial,
                error_residual=_zero_tensor(),
                coupling_F_lower=Q(1),
                declared_tangent=(Q(1), Q(1)),
            )

    def test_missing_context_unit_provenance_and_aliases_refuse(self) -> None:
        factories = _nonzero_records()
        aliases = (True, 1.0, float("inf"), float("nan"), "1")
        for name, record in factories.items():
            for field in ("context", "unit", "provenance"):
                with self.subTest(name=name, field=field):
                    with self.assertRaises(Def1Stab1ProviderError):
                        replace(record, **{field: ""})
            with self.subTest(name=name, field="negative"):
                with self.assertRaises(Def1Stab1ProviderError):
                    replace(record, additive_q=Q(-1, 8))
        with self.assertRaises(TypeError):
            physical_constraint_provider(
                residual_E=_zero_tensor(),
                tangent=(Q(1), Q(1)),
                coupling_F=True,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl("physical_constraint"),
            )
        for bad in aliases:
            with self.subTest(bad=bad):
                radii = complete_zero_radii()
                radii["R.dtt"] = bad
                with self.assertRaises(TypeError):
                    initial_data_provider(
                        radii=radii,
                        status=PREMISE_STATUS_PROVEN,
                        source=SOURCE_SUPPLIED_CERTIFIED,
                        **_decl("initial_data"),
                    )

    def test_richardson_and_imp1_stay_conditional_not_pde(self) -> None:
        with self.assertRaisesRegex(Def1Stab1ProviderError, "conditional premise"):
            spatial_temporal_provider(
                source=SOURCE_RICHARDSON,
                status=PREMISE_STATUS_PROVEN,
                additive_q=Q(1, 8),
                **_decl("spatial_temporal"),
            )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "zero-only stubs"):
            spatial_temporal_provider(
                source=SOURCE_RICHARDSON,
                status=PREMISE_STATUS_CONDITIONAL,
                additive_q=0,
                **_decl("spatial_temporal"),
            )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "only be declared on spatial_temporal"):
            ProviderRecord(
                component="conservation",
                radii=tuple((name, Q(0)) for name in INPUT_NAMES),
                additive_q=Q(1, 8),
                status=PREMISE_STATUS_CONDITIONAL,
                source=SOURCE_RICHARDSON,
                conversion=CONVERSION_RICHARDSON,
                **_decl("conservation"),
            )
        assembled = assemble_provider_error_map(
            _minkowski(),
            (
                spatial_temporal_provider(
                    source=SOURCE_IMP1_ADMISSION_DEBIT,
                    status=PREMISE_STATUS_CONDITIONAL,
                    additive_q=Q(1, 5),
                    **_decl("spatial_temporal"),
                ),
            )
            + _assembly_records()[1:],
            mass_flux=_mass_flux(),
        )
        self.assertEqual(assembled.contributions[0].total, Q(1, 5))
        self.assertFalse(assembled.assembled.all_components_declared_proven_enclosures)
        self.assertFalse(assembled.global_pde_error_certified)
        self.assertFalse(assembled.imp1_admission_debit_treated_as_global_pde_error)
        self.assertEqual(assembled.records[0].conversion, CONVERSION_IMP1)

    def test_physical_and_gauge_require_complete_residual_not_c_or_hm(self) -> None:
        with self.assertRaisesRegex(Def1Stab1ProviderError, "H/M-only"):
            physical_constraint_provider(
                residual_E=((0, 0), (0, 0)),
                tangent=(1, 1),
                coupling_F=1,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl("physical_constraint"),
            )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "unexpected arguments"):
            physical_constraint_provider(
                residual_E=_zero_tensor(),
                tangent=(1, 1),
                coupling_F=1,
                hamiltonian=Q(1, 8),
                momentum=Q(1, 8),
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl("physical_constraint"),
            )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "C-only"):
            gauge_constraint_provider(
                constraint_C=(Q(1, 4), Q(0)),
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl("gauge_constraint"),
            )
        converted = gauge_constraint_provider(
            constraint_C=(Q(1, 4), 0, 0, 0),
            c_to_jet_inverse_bound=2,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl("gauge_constraint"),
        )
        self.assertEqual(converted.radii_mapping["alpha.value"], Q(1, 2))
        self.assertTrue(all(value == Q(1, 2) for value in converted.radii_mapping.values()))

    def test_source_inversion_requires_full_residual_and_inverse_j(self) -> None:
        with self.assertRaisesRegex(Def1Stab1ProviderError, "length 6"):
            nonlinear_source_provider(
                full_residual=(Q(1, 8), 0),
                inverse_jacobian_bound=2,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl("nonlinear_source"),
            )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "cannot default to zero"):
            nonlinear_source_provider(
                full_residual=(Q(1, 8), 0, 0, 0, 0, 0),
                inverse_jacobian_bound=0,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl("nonlinear_source"),
            )
        with self.assertRaises(TypeError):
            nonlinear_source_provider(
                full_residual=(Q(1, 8), 0, 0, 0, 0, 0),
                inverse_jacobian_bound=0.5,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl("nonlinear_source"),
            )

    def test_reduction_requires_explicit_discrepancies_and_operator_bound(self) -> None:
        partial = dict(_zero_reduction())
        del partial["R.dr"]
        with self.assertRaisesRegex(Def1Stab1ProviderError, "missing="):
            reduction_constraint_provider(
                first_derivative_discrepancies=partial,
                derivative_operator_bound=1,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl("reduction_constraint"),
            )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "cannot default to zero"):
            reduction_constraint_provider(
                first_derivative_discrepancies=_zero_reduction() | {"lambda.dt": 1},
                derivative_operator_bound=0,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl("reduction_constraint"),
            )

    def test_affine_requires_gronwall_not_residual_times_step(self) -> None:
        with self.assertRaisesRegex(Def1Stab1ProviderError, "Gronwall"):
            affine_residual_times_step_refused(
                transport_defect=Q(1, 10), affine_interval=Q(1, 8)
            )
        with self.assertRaises(TypeError):
            affine_provider(
                transport_defect=Q(1, 10),
                affine_interval=Q(1, 8),
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl("affine"),
            )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "residual times step"):
            affine_provider(
                transport_defect=Q(1, 10),
                affine_interval=Q(1, 8),
                gronwall_factor=0,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl("affine"),
            )

    def test_interpolation_lipschitz_is_not_a_second_order_remainder(self) -> None:
        with self.assertRaisesRegex(Def1Stab1ProviderError, "derivative bounds are missing"):
            interpolation_provider(
                sample_spacing=Q(1, 8),
                first_derivative_enclosures={},
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl("interpolation"),
            )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "second-order"):
            interpolation_provider(
                sample_spacing=Q(1, 8),
                first_derivative_enclosures=_zero_interpolation() | {"R.value": 1},
                claim_second_order_remainder=True,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl("interpolation"),
            )
        remainder_zeros = complete_zero_radii()
        with self.assertRaisesRegex(Def1Stab1ProviderError, "empty derivative"):
            interpolation_provider(
                sample_spacing=Q(1, 8),
                first_derivative_enclosures=_zero_interpolation(),
                second_derivative_enclosures={},
                remainder_radii=remainder_zeros,
                claim_second_order_remainder=True,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl("interpolation"),
            )
        genuine = interpolation_provider(
            sample_spacing=Q(1, 4),
            first_derivative_enclosures=_zero_interpolation(),
            second_derivative_enclosures={"R.value": Q(2)},
            remainder_radii={
                name: Q(0) for name in INPUT_NAMES if name != "R.value"
            },
            claim_second_order_remainder=True,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_SUPPLIED_CERTIFIED,
            **_decl("interpolation"),
        )
        self.assertEqual(genuine.conversion, CONVERSION_INTERPOLATION_SECOND_ORDER)
        self.assertEqual(genuine.radii_mapping["R.value"], Q(1, 16))
        with self.assertRaisesRegex(Def1Stab1ProviderError, "missing="):
            interpolation_provider(
                sample_spacing=Q(1, 8),
                first_derivative_enclosures={"R.value": 1},
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl("interpolation"),
            )

    def test_extraction_unsampled_gaps_are_not_positive_samples(self) -> None:
        with self.assertRaisesRegex(Def1Stab1ProviderError, "continuous-Q"):
            extraction_provider(
                affine_samples=(0, Q(1, 2), 1),
                derivative_bound=None,
                sample_q_values=(1, 1, 1),
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl("extraction"),
            )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "at least two samples"):
            extraction_provider(
                affine_samples=(0,),
                derivative_bound=1,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl("extraction"),
            )

    def test_boundary_causality_is_not_zero_numerical_error(self) -> None:
        with self.assertRaisesRegex(Def1Stab1ProviderError, "physical causality alone"):
            boundary_provider(
                additional_debit=0,
                physical_causality_passed=True,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl("boundary"),
            )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "numerical-boundary guard"):
            boundary_provider(
                additional_debit=0,
                no_influence_premise=True,
                coverage_premise=True,
                independent_guard=False,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl("boundary"),
            )

    def test_conservation_mass_flux_veto_and_coverage(self) -> None:
        vetoed = _mass_flux(valid=False)
        with self.assertRaisesRegex(
            Def1Stab1ProviderError, PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN
        ):
            conservation_provider(
                additional_debit=0,
                mass_flux=vetoed,
                coverage_premise=True,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl("conservation"),
            )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "coverage premise"):
            conservation_provider(
                additional_debit=0,
                mass_flux=_mass_flux(),
                coverage_premise=False,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl("conservation"),
            )
        assembled_records = _assembly_records()
        with self.assertRaisesRegex(
            Def1Stab1ProviderError, PROTOCOL_MASS_FLUX_INCONSISTENCY_TOKEN
        ):
            assemble_provider_error_map(
                _minkowski(), assembled_records, mass_flux=vetoed
            )

    def test_initial_data_and_alignment_require_complete_radii(self) -> None:
        incomplete = complete_zero_radii()
        del incomplete["R.drr"]
        with self.assertRaisesRegex(Def1Stab1ProviderError, "26 named inputs"):
            initial_data_provider(
                radii=incomplete,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl("initial_data"),
            )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "26 named inputs"):
            trajectory_alignment_provider(
                radii=incomplete,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl("trajectory_alignment"),
            )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "nonnegative"):
            initial_data_provider(
                radii=_radii_with(**{"R.dtt": Q(-1, 8)}),
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                **_decl("initial_data"),
            )


class AssemblyGuardTests(unittest.TestCase):
    def test_contributions_are_outward_and_ignore_positive_q(self) -> None:
        initial = _radii_with(**{"R.dtt": Q(1, 32)})
        minkowski = assemble_provider_error_map(
            _minkowski(),
            _assembly_records(initial=initial),
            mass_flux=_mass_flux(),
        )
        positive = _minkowski() | {"R.dtt": Q(1)}
        self.assertGreater(evaluate_q_box(positive)["q"].lower, 0)
        defocusing = assemble_provider_error_map(
            positive,
            _assembly_records(initial=initial),
            mass_flux=_mass_flux(),
        )
        self.assertEqual(minkowski.assembled.exact_total, defocusing.assembled.exact_total)
        self.assertEqual(minkowski.assembled.exact_total, Q(1, 32))
        for item in minkowski.contributions:
            stored = Q(getattr(minkowski.assembled.budget, item.component))
            self.assertGreaterEqual(item.total, 0)
            self.assertGreaterEqual(stored, item.total)
        self.assertGreaterEqual(
            Q(minkowski.assembled.budget.total_upper_bound),
            minkowski.assembled.exact_total,
        )
        self.assertEqual(
            minkowski.assembled.exact_total,
            sum((item.total for item in minkowski.contributions), Q(0)),
        )
        nonzero = _nonzero_records()
        complete = assemble_provider_error_map(
            _minkowski(),
            tuple(nonzero[name] for name in ERROR_BUDGET_COMPONENTS),
            mass_flux=_mass_flux(),
        )
        self.assertEqual(len(complete.contributions), 14)
        self.assertEqual(
            complete.assembled.exact_total,
            sum((item.total for item in complete.contributions), Q(0)),
        )
        for item in complete.contributions:
            stored = Q(getattr(complete.assembled.budget, item.component))
            self.assertGreater(item.total, 0)
            self.assertGreaterEqual(stored, item.total)
        self.assertGreaterEqual(
            Q(complete.assembled.budget.total_upper_bound),
            complete.assembled.exact_total,
        )
        self.assertFalse(complete.assembled.used_measured_q)
        with self.assertRaisesRegex(Def1Stab1ProviderError, "measured Q"):
            assemble_provider_error_map(
                _minkowski(),
                _assembly_records(),
                mass_flux=_mass_flux(),
                measured_q=Q(1000),
            )

    def test_forged_records_and_false_gate_refuse(self) -> None:
        assembled = assemble_provider_error_map(
            _minkowski(), _assembly_records(), mass_flux=_mass_flux()
        )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "rederivation"):
            replace(assembled, def1_error_map_passed=True)
        with self.assertRaisesRegex(Def1Stab1ProviderError, "rederivation"):
            replace(assembled, used_measured_q=True)
        with self.assertRaisesRegex(Def1Stab1ProviderError, "rederivation"):
            replace(assembled, global_pde_error_certified=True)
        with self.assertRaisesRegex(Def1Stab1ProviderError, "rederivation"):
            replace(assembled, remaining_gate_work=())
        drifted = list(assembled.contributions)
        drifted[0] = replace(drifted[0], total=Q(1), radii_term=Q(1))
        with self.assertRaisesRegex(Def1Stab1ProviderError, "rederivation"):
            replace(assembled, contributions=tuple(drifted))
        with self.assertRaisesRegex(Def1Stab1ProviderError, "exactly the frozen"):
            assemble_provider_error_map(
                _minkowski(),
                _assembly_records()[:-1],
                mass_flux=_mass_flux(),
            )

    def test_frozen_thresholds_and_nine_booleans_remain_independent(self) -> None:
        self.assertEqual(MIN_ACTIVATION_LOWER_BOUND, 2)
        self.assertEqual(TRAPPEDNESS_ERROR_FACTOR, 4)
        self.assertEqual(COMPLETE_Q_ERROR_FACTOR, 4)
        self.assertEqual(MIN_AFFINE_SAMPLES, 8)
        self.assertEqual(MIN_AFFINE_INTERVAL_OVER_L0, Q(1, 64))
        self.assertEqual(DIRECT_RAYCHAUDHURI_AGREEMENT_MAX, Q(1, 100_000_000))
        self.assertEqual(MIN_NESTED_RESOLUTIONS, 3)
        self.assertEqual(MIN_INDEPENDENT_METHODS, 2)
        self.assertEqual(MIN_OBSERVED_ORDER, Q(3, 2))
        booleans = record_def1_booleans(
            resolved_activation=True,
            control_dominance=False,
            resolved_trapped_interval=True,
            resolved_complete_defocusing=False,
            direct_Raychaudhuri_agreement=True,
            all_health_constraints_scales_valid=False,
            mass_flux_ledger_valid=True,
            three_resolution_convergence=False,
            two_method_agreement=True,
        )
        assembled = assemble_provider_error_map(
            _minkowski(),
            _assembly_records(),
            mass_flux=_mass_flux(),
            def1_booleans=booleans,
        )
        self.assertEqual(assembled.def1_booleans, booleans)
        self.assertEqual(tuple(booleans.as_mapping()), DEF1_BOOLEAN_NAMES)
        self.assertFalse(assembled.def1_error_map_passed)
        self.assertFalse(booleans.as_mapping()["three_resolution_convergence"])
        self.assertTrue(assembled.mass_flux.ledger_valid)

    def test_remaining_gate_work_keeps_the_map_unpassed(self) -> None:
        assembled = assemble_provider_error_map(
            _minkowski(), _assembly_records(), mass_flux=_mass_flux()
        )
        self.assertIn("DEF1_error_map_passed remains false", " ".join(assembled.remaining_gate_work))
        self.assertIn("COL1", " ".join(assembled.remaining_gate_work))
        self.assertIn("not a prerequisite", " ".join(assembled.remaining_gate_work))
        self.assertFalse(assembled.def1_error_map_passed)
        self.assertIsInstance(assembled, AssembledProviderMap)
        self.assertGreaterEqual(len(assembled.remaining_gate_work), 8)


def _cancelling_residual() -> tuple[tuple[Q, ...], ...]:
    rows = [list(row) for row in _zero_tensor()]
    rows[0][0] = Q(7)
    rows[1][1] = Q(-7)
    return tuple(tuple(row) for row in rows)


class ReviewedGuardTests(unittest.TestCase):
    def test_replace_joint_box_and_zero_jacobian_cannot_erase_r_dr_bound(self) -> None:
        assembled = assemble_provider_error_map(
            _minkowski(),
            _assembly_records(initial=_radii_with(**{"R.dr": Q(1, 8)})),
            mass_flux=_mass_flux(),
        )
        self.assertEqual(
            assembled.contributions[
                ERROR_BUDGET_COMPONENTS.index("initial_data")
            ].total,
            Q(9, 64),
        )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "rederivation"):
            replace(assembled, joint_box=())
        zero_lipschitz = tuple(
            (name, Interval.singleton(0)) for name in INPUT_NAMES
        )
        zero_contributions = tuple(
            ComponentContribution(item.component, Q(0), Q(0), Q(0))
            for item in assembled.contributions
        )
        zero_budget = assemble_q_error_budget(
            tuple(
                QComponentPremise(
                    component=record.component,
                    sensitivity=Interval.singleton(1),
                    input_error=Interval.singleton(0),
                    status=record.status,
                    source=record.source,
                )
                for record in assembled.records
            )
        )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "rederivation"):
            replace(
                assembled,
                lipschitz=zero_lipschitz,
                contributions=zero_contributions,
                assembled=zero_budget,
            )

    def test_duplicate_radius_names_and_aliases_refuse(self) -> None:
        zeros = tuple((name, Q(0)) for name in INPUT_NAMES)
        duplicated = (("R.dr", Q(1, 8)), ("R.dr", Q(1, 4))) + zeros[2:]
        with self.assertRaisesRegex(Def1Stab1ProviderError, "duplicate radius names"):
            ProviderRecord(
                component="initial_data",
                radii=duplicated,
                additive_q=0,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                conversion=CONVERSION_SUPPLIED_RADII,
                **_decl("initial_data"),
            )
        aliased = (("radius.dr", Q(1, 8)),) + zeros[1:]
        with self.assertRaisesRegex(Def1Stab1ProviderError, "aliases"):
            ProviderRecord(
                component="initial_data",
                radii=aliased,
                additive_q=0,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_SUPPLIED_CERTIFIED,
                conversion=CONVERSION_SUPPLIED_RADII,
                **_decl("initial_data"),
            )

    def test_zero_debit_guards_require_membership_not_unrelated_names(self) -> None:
        boundary = boundary_provider(
            additional_debit=0,
            no_influence_premise=True,
            coverage_premise=True,
            independent_guard=True,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl("boundary"),
        )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "unrelated"):
            replace(boundary, guards=("unrelated",))
        with self.assertRaisesRegex(Def1Stab1ProviderError, "unrelated"):
            replace(
                boundary,
                guards=("coverage_premise", "unrelated"),
            )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "duplicate guards"):
            replace(
                boundary,
                guards=(
                    "no_influence_premise",
                    "coverage_premise",
                    "independent_guard",
                    "independent_guard",
                ),
            )
        conservation = conservation_provider(
            additional_debit=0,
            mass_flux=_mass_flux(),
            coverage_premise=True,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl("conservation"),
        )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "unrelated"):
            replace(conservation, guards=("unrelated",))
        with self.assertRaisesRegex(Def1Stab1ProviderError, "missing"):
            replace(conservation, guards=("mass_flux_ledger_valid",))

    def test_signed_ekk_cancellation_is_not_the_q_error_debit(self) -> None:
        residual = _cancelling_residual()
        self.assertEqual(residual_kk(residual, (0, 0)), 0)
        self.assertEqual(residual_kk(residual, (1, 1)), 0)
        with self.assertRaisesRegex(Def1Stab1ProviderError, "future-directed"):
            residual_q_debit(residual, (0, 0), 1)
        with self.assertRaisesRegex(Def1Stab1ProviderError, "future-directed"):
            physical_constraint_provider(
                residual_E=residual,
                tangent=(0, 0),
                coupling_F=1,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl("physical_constraint"),
            )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "angular"):
            physical_constraint_provider(
                residual_E=residual,
                tangent=(1, 0, 1, 0),
                coupling_F=1,
                status=PREMISE_STATUS_PROVEN,
                source=SOURCE_EXACT_ALGEBRA,
                **_decl("physical_constraint"),
            )
        self.assertEqual(conservative_residual_q_debit(residual, (1, 1), 1), 14)
        record = physical_constraint_provider(
            residual_E=residual,
            tangent=(Q(1), Q(1)),
            coupling_F=1,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl("physical_constraint"),
        )
        self.assertEqual(record.additive_q, 14)
        assembled = assemble_provider_error_map(
            _minkowski(),
            (_assembly_records()[0], record) + _assembly_records()[2:],
            mass_flux=_mass_flux(),
        )
        self.assertEqual(
            assembled.contributions[
                ERROR_BUDGET_COMPONENTS.index("physical_constraint")
            ].additive_q,
            14,
        )
        foreign = physical_constraint_provider(
            residual_E=residual,
            tangent=(Q(2), Q(1)),
            coupling_F=1,
            status=PREMISE_STATUS_PROVEN,
            source=SOURCE_EXACT_ALGEBRA,
            **_decl("physical_constraint"),
        )
        with self.assertRaisesRegex(Def1Stab1ProviderError, "foreign"):
            assemble_provider_error_map(
                _minkowski(),
                (_assembly_records()[0], foreign) + _assembly_records()[2:],
                mass_flux=_mass_flux(),
            )

    def test_mixed_evaluation_context_or_units_refuse(self) -> None:
        records = list(_assembly_records())
        records[0] = replace(records[0], unit="other units")
        with self.assertRaisesRegex(Def1Stab1ProviderError, "shared evaluation-context"):
            assemble_provider_error_map(
                _minkowski(), tuple(records), mass_flux=_mass_flux()
            )
        records = list(_assembly_records())
        records[3] = replace(records[3], context="a different evaluation chart")
        with self.assertRaisesRegex(Def1Stab1ProviderError, "shared evaluation-context"):
            assemble_provider_error_map(
                _minkowski(), tuple(records), mass_flux=_mass_flux()
            )


if __name__ == "__main__":
    unittest.main()
