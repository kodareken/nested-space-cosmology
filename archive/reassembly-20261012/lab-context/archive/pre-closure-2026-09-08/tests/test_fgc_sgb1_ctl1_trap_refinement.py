from __future__ import annotations

import ast
from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.exact_interval import interval  # noqa: E402
from recursive_horizons.fgc.sgb1_ctl1_initial_health import (  # noqa: E402
    CERTIFICATE_CHART_KIND,
    DECLARED_BASE_CELLS,
    DECLARED_C_DOMAIN,
    DECLARED_K_DOMAIN,
    DECLARED_LAMBDA_DOMAIN,
    DECLARED_MAX_PICARD_ITERATIONS,
    DECLARED_MAX_RATIONAL_BITS,
    SGBLExactInitialSlice,
    SGBLODECellEnclosure,
    sgbl_validate_ode_cell_inventory,
    sgbl_validated_lambda_k_constraint_ode,
)
from recursive_horizons.fgc.sgb1_ctl1_trap_refinement import (  # noqa: E402
    DECLARED_AFFINE_EVALUATIONS_PER_RHS_CALL,
    DECLARED_COMPACTNESS_STRICT_UPPER,
    DECLARED_LADDER_DEPTHS,
    DECLARED_PRODUCT_BOX_RESOURCE_LADDER,
    DECLARED_RHS_CALLS_PER_PARTITION,
    INSTRUMENT_ID,
    NOMINAL_CHI_AMPLITUDE,
    PRODUCT_BOX_CHART_KIND,
    SMALL_AMPLITUDE_CONTROL,
    SGBLProductBoxLadderRecord,
    SGBLTrapRefinementStop,
    sgbl_declared_product_box_bit_cap,
    sgbl_declared_product_box_cell_cap,
    sgbl_declared_product_box_rhs_cap,
    sgbl_evaluate_product_box_ladder_level,
    sgbl_evaluate_product_box_resource_ladder,
    sgbl_product_box_cells_cover_support,
    sgbl_product_box_inventory_span,
    sgbl_product_box_shared_coverage_monotone,
    sgbl_trap_refinement_contract_sha256,
    sgbl_trap_refinement_health_gate,
)


DEPTH8_PRODUCT_BOX_C_UPPER = Q(18453842205463973737, 2**64)
DEPTH8_PRODUCT_BOX_C_MARGIN = Q(-7098131754422121, 2**64)
DEPTH10_PRODUCT_BOX_C_UPPER = Q(9223921512318141519, 2**63)
DEPTH10_PRODUCT_BOX_C_MARGIN = Q(-549475463365711, 2**63)
DEPTH12_PRODUCT_BOX_C_UPPER = Q(4611781872439090751, 2**62)
DEPTH12_PRODUCT_BOX_C_MARGIN = Q(-95854011702847, 2**62)
CORRECTED_CELL_COUNTS = (18, 21, 25)
CORRECTED_COVERAGE_RIGHTS = (Q(727, 64), Q(23265, 2048), Q(186123, 16384))
CORRECTED_RHS_EVALUATIONS = (384, 464, 536)
CORRECTED_COMPACTNESS_RHS = (96, 116, 134)
NOMINAL_CONTRACT_SHA256 = (
    "a82363a29b641c538a0e15f14ec2f4c99a49182a9a01c84ae7b1d14a953834c8"
)
SMALL_AMPLITUDE_CONTRACT_SHA256 = (
    "c3cca0ef56f4f1e257bd0d2b3e0ec088505aaae2c4e782950d4c632102afd299"
)


def _dummy_cell(*, radius, compactness, lam=None, k=None) -> SGBLODECellEnclosure:
    zero = interval(0)
    one = interval(1)
    lam = one if lam is None else lam
    k = zero if k is None else k
    return SGBLODECellEnclosure(
        radius=radius,
        lambda_box=lam,
        k_box=k,
        compactness_box=compactness,
        j_box=zero,
        lambda_left=lam,
        k_left=k,
        compactness_left=compactness,
        j_left=zero,
        lambda_right=lam,
        k_right=k,
        compactness_right=compactness,
        j_right=zero,
        compactness=compactness,
        direct_compactness=compactness,
        centered_compactness=compactness,
        propagated_compactness=compactness,
        invariant_residual=zero,
        denominator=one,
        denominator_margin=Q(1),
        chart_coordinates=False,
        chart_kind=PRODUCT_BOX_CHART_KIND,
        depth=0,
        picard_strict_self_map=True,
        residual_contains_origin=True,
        compactness_invariant_contains_zero=True,
        correlated_strictly_tighter=False,
    )


class SGBLProductBoxLadderFreezeTests(unittest.TestCase):
    def test_ladder_is_frozen_before_any_ode_evaluation(self) -> None:
        source = Path(
            ROOT / "src/recursive_horizons/fgc/sgb1_ctl1_trap_refinement.py"
        ).read_text()
        tree = ast.parse(source)
        names: list[str] = []
        for node in tree.body:
            if isinstance(node, ast.Assign):
                names.extend(
                    target.id for target in node.targets if isinstance(target, ast.Name)
                )
            elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
                names.append(node.target.id)
            elif isinstance(node, ast.FunctionDef):
                names.append(node.name)
        self.assertLess(
            names.index("DECLARED_PRODUCT_BOX_RESOURCE_LADDER"),
            names.index("sgbl_evaluate_product_box_resource_ladder"),
        )
        self.assertLess(
            names.index("DECLARED_PRODUCT_BOX_RESOURCE_LADDER"),
            names.index("sgbl_evaluate_product_box_ladder_level"),
        )
        self.assertIn("sgbl_validated_lambda_k_constraint_ode", source)
        self.assertIn("sgbl_validate_ode_cell_inventory", source)
        self.assertNotIn("sgbl_validated_constraint_ode(", source)

    def test_declared_caps_are_combinatorial_and_do_not_change_c_or_domains(self) -> None:
        self.assertEqual(DECLARED_LADDER_DEPTHS, (8, 10, 12))
        self.assertEqual(DECLARED_BASE_CELLS, 16)
        self.assertEqual(DECLARED_COMPACTNESS_STRICT_UPPER, 1)
        self.assertEqual(NOMINAL_CHI_AMPLITUDE, 3)
        self.assertEqual(SGBLExactInitialSlice().chi_amplitude, NOMINAL_CHI_AMPLITUDE)
        self.assertEqual(len(DECLARED_PRODUCT_BOX_RESOURCE_LADDER), 3)
        for depth, level in zip(DECLARED_LADDER_DEPTHS, DECLARED_PRODUCT_BOX_RESOURCE_LADDER):
            with self.subTest(depth=depth):
                self.assertEqual(level.max_bisection_depth, depth)
                self.assertEqual(level.base_cells, DECLARED_BASE_CELLS)
                self.assertEqual(level.max_picard_iterations, DECLARED_MAX_PICARD_ITERATIONS)
                self.assertEqual(level.lambda_domain, DECLARED_LAMBDA_DOMAIN)
                self.assertEqual(level.k_domain, DECLARED_K_DOMAIN)
                self.assertEqual(level.compactness_domain, DECLARED_C_DOMAIN)
                self.assertEqual(
                    level.max_cells, sgbl_declared_product_box_cell_cap(depth)
                )
                self.assertEqual(
                    level.max_rhs_evaluations, sgbl_declared_product_box_rhs_cap(depth)
                )
                self.assertEqual(
                    level.max_rational_bit_length,
                    sgbl_declared_product_box_bit_cap(depth),
                )
                self.assertGreater(level.compactness_domain.upper, 1)
                self.assertLess(level.compactness_domain.lower, 0)
                policy = level.policy()
                self.assertEqual(policy.max_bisection_depth, depth)
                self.assertEqual(policy.base_cells, 16)
        self.assertEqual(sgbl_declared_product_box_cell_cap(8), 16 * (2**9 - 1))
        self.assertEqual(sgbl_declared_product_box_cell_cap(10), 16 * (2**11 - 1))
        self.assertEqual(sgbl_declared_product_box_cell_cap(12), 16 * (2**13 - 1))
        self.assertEqual(
            sgbl_declared_product_box_rhs_cap(8),
            sgbl_declared_product_box_cell_cap(8)
            * DECLARED_RHS_CALLS_PER_PARTITION
            * DECLARED_AFFINE_EVALUATIONS_PER_RHS_CALL,
        )
        self.assertEqual(sgbl_declared_product_box_bit_cap(8), DECLARED_MAX_RATIONAL_BITS)
        self.assertEqual(sgbl_declared_product_box_bit_cap(10), 2 * DECLARED_MAX_RATIONAL_BITS)
        self.assertEqual(sgbl_declared_product_box_bit_cap(12), 4 * DECLARED_MAX_RATIONAL_BITS)
        self.assertGreater(DECLARED_PRODUCT_BOX_RESOURCE_LADDER[1].max_cells, 1024)
        self.assertGreater(DECLARED_PRODUCT_BOX_RESOURCE_LADDER[1].max_rhs_evaluations, 8192)

    def test_changed_domains_or_base_cells_are_refused_by_the_level_contract(self) -> None:
        with self.assertRaises(ValueError):
            replace(DECLARED_PRODUCT_BOX_RESOURCE_LADDER[0], base_cells=8)
        with self.assertRaises(ValueError):
            replace(
                DECLARED_PRODUCT_BOX_RESOURCE_LADDER[0],
                compactness_domain=interval(Q(-1), Q(1)),
            )
        with self.assertRaises(ValueError):
            replace(
                DECLARED_PRODUCT_BOX_RESOURCE_LADDER[0],
                lambda_domain=interval(Q(1, 4), Q(8)),
            )


class SGBLProductBoxLadderNominalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec = SGBLExactInitialSlice()
        cls.record = sgbl_evaluate_product_box_resource_ladder(cls.spec)

    def test_every_declared_level_runs_and_does_not_stop_at_the_first_outcome(self) -> None:
        self.assertEqual(len(self.record.levels), 3)
        self.assertTrue(self.record.evaluated_every_declared_level)
        self.assertEqual(
            tuple(level.level.max_bisection_depth for level in self.record.levels),
            DECLARED_LADDER_DEPTHS,
        )
        for record in self.record.levels:
            self.assertIsNotNone(record.ode)
            self.assertEqual(record.chart_kind, PRODUCT_BOX_CHART_KIND)
            self.assertNotEqual(record.chart_kind, CERTIFICATE_CHART_KIND)
            self.assertTrue(record.product_box_is_not_the_certificate_chart)
            self.assertTrue(record.sampled_nodes_are_not_the_certificate)
            self.assertGreaterEqual(record.rhs_evaluations, 1)
            self.assertGreaterEqual(record.ode_cell_count, 1)
            self.assertIsInstance(record.compactness_upper, Q)
            self.assertIsInstance(record.compactness_margin, Q)
            self.assertEqual(
                record.compactness_margin,
                1 - max(record.compactness_upper, record.exterior_compactness_upper),
            )

    def test_nominal_levels_bind_the_gap_free_prefix_inventory(self) -> None:
        expected_uppers = (
            DEPTH8_PRODUCT_BOX_C_UPPER,
            DEPTH10_PRODUCT_BOX_C_UPPER,
            DEPTH12_PRODUCT_BOX_C_UPPER,
        )
        expected_margins = (
            DEPTH8_PRODUCT_BOX_C_MARGIN,
            DEPTH10_PRODUCT_BOX_C_MARGIN,
            DEPTH12_PRODUCT_BOX_C_MARGIN,
        )
        previous_cells = 0
        for index, record in enumerate(self.record.levels):
            with self.subTest(depth=record.level.max_bisection_depth):
                self.assertEqual(record.classification, "interval_inconclusive")
                self.assertEqual(record.obstruction, "compactness_enclosure_not_below_one")
                self.assertFalse(record.complete_support_coverage)
                self.assertFalse(record.continuous_no_initial_trapped_sphere)
                self.assertFalse(record.ode.tiles_compact_support)
                self.assertFalse(record.ode.covers_complete_support)
                self.assertEqual(record.ode_cell_count, CORRECTED_CELL_COUNTS[index])
                self.assertGreater(record.ode_cell_count, previous_cells)
                self.assertEqual(record.rhs_evaluations, CORRECTED_RHS_EVALUATIONS[index])
                self.assertEqual(
                    record.compactness_rhs_evaluations, CORRECTED_COMPACTNESS_RHS[index]
                )
                self.assertEqual(record.compactness_upper, expected_uppers[index])
                self.assertEqual(record.compactness_margin, expected_margins[index])
                self.assertGreaterEqual(record.compactness_upper, 1)
                self.assertLess(record.compactness_margin, 0)
                self.assertEqual(record.coverage_left, self.spec.support_minimum)
                self.assertEqual(record.coverage_right, CORRECTED_COVERAGE_RIGHTS[index])
                self.assertEqual(record.coverage_right, record.ode.cells[-1].radius.upper)
                self.assertLess(record.coverage_right, self.spec.support_maximum)
                left, right = sgbl_validate_ode_cell_inventory(
                    record.ode.cells, origin=self.spec.support_minimum
                )
                self.assertEqual((left, right), (record.coverage_left, record.coverage_right))
                self.assertEqual(
                    sgbl_product_box_inventory_span(record.ode.cells, self.spec),
                    (record.coverage_left, record.coverage_right),
                )
                self.assertFalse(
                    sgbl_product_box_cells_cover_support(record.ode.cells, self.spec)
                )
                self.assertTrue(
                    all(
                        record.ode.cells[i].radius.upper
                        == record.ode.cells[i + 1].radius.lower
                        for i in range(len(record.ode.cells) - 1)
                    )
                )
                previous_cells = record.ode_cell_count
        independent = sgbl_validated_lambda_k_constraint_ode(
            self.spec, policy=self.record.levels[0].level.policy()
        )
        self.assertEqual(len(independent.cells), 18)
        self.assertEqual(independent.support_compactness.upper, DEPTH8_PRODUCT_BOX_C_UPPER)
        self.assertEqual(independent.obstruction, "compactness_enclosure_not_below_one")
        self.assertTrue(self.record.compactness_upper_nonincreasing)
        self.assertTrue(self.record.shared_coverage_monotone)
        self.assertEqual(self.record.proved_level_names, ())
        self.assertFalse(self.record.any_declared_level_proves_continuous_no_initial_trap)
        self.assertEqual(self.record.contract_sha256, NOMINAL_CONTRACT_SHA256)

    def test_wall_and_resources_are_recorded_for_every_level(self) -> None:
        for record in self.record.levels:
            with self.subTest(depth=record.level.max_bisection_depth):
                self.assertGreaterEqual(record.diagnostic_wall_seconds, 0)
                self.assertLessEqual(
                    record.ode_cell_count, record.level.max_cells
                )
                self.assertLessEqual(
                    record.rhs_evaluations, record.level.max_rhs_evaluations
                )
                self.assertLessEqual(
                    record.max_rational_bit_length_observed,
                    record.level.max_rational_bit_length,
                )
                print(
                    f"product-box depth {record.level.max_bisection_depth}: "
                    f"wall={record.diagnostic_wall_seconds:.4f}s "
                    f"cells={record.ode_cell_count} "
                    f"rhs={record.rhs_evaluations} "
                    f"c_rhs={record.compactness_rhs_evaluations} "
                    f"bits={record.max_rational_bit_length_observed} "
                    f"C={record.compactness_upper} "
                    f"margin={record.compactness_margin} "
                    f"coverage={record.complete_support_coverage} "
                    f"span=[{record.coverage_left},{record.coverage_right}] "
                    f"obstruction={record.obstruction}",
                    flush=True,
                )

    def test_aggregate_health_and_nonpromoting_contract_stay_false(self) -> None:
        gate = sgbl_trap_refinement_health_gate(self.record)
        self.assertIs(self.record.SGBL_branch_owned_and_healthy, False)
        self.assertIs(self.record.FRZ1, False)
        self.assertIs(self.record.PREF1, False)
        self.assertIs(self.record.holdout_authorized, False)
        self.assertIs(self.record.execution_authorized, False)
        self.assertIs(gate["SGBL_branch_owned_and_healthy"], False)
        self.assertIs(gate["FRZ1"], False)
        self.assertIs(gate["PREF1"], False)
        self.assertIs(gate["holdout_authorized"], False)
        self.assertIs(gate["execution_authorized"], False)
        self.assertEqual(self.record.contract_payload["INSTRUMENT_ID"], INSTRUMENT_ID)
        self.assertIs(
            self.record.contract_payload["SGBL_branch_owned_and_healthy"], False
        )
        self.assertEqual(
            self.record.contract_sha256,
            sgbl_trap_refinement_contract_sha256(self.record.contract_payload),
        )
        self.assertEqual(self.record.contract_sha256, NOMINAL_CONTRACT_SHA256)
        self.assertNotIn("wall", str(dict(self.record.contract_payload)))
        self.assertNotIn("diagnostic_wall_seconds", dict(self.record.contract_payload))
        self.assertIn("coverage_left", self.record.contract_payload["levels"][0])
        self.assertIn("coverage_right", self.record.contract_payload["levels"][0])


class SGBLProductBoxLadderPositiveControlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec = SGBLExactInitialSlice(chi_amplitude=SMALL_AMPLITUDE_CONTROL)
        cls.record = sgbl_evaluate_product_box_resource_ladder(cls.spec)

    def test_small_amplitude_passes_every_declared_level_without_stopping_early(self) -> None:
        self.assertEqual(self.spec.chi_amplitude, Q(1, 8))
        self.assertEqual(len(self.record.levels), 3)
        for record in self.record.levels:
            self.assertEqual(record.classification, "continuous_no_initial_trapped_sphere")
            self.assertTrue(record.continuous_no_initial_trapped_sphere)
            self.assertTrue(record.complete_support_coverage)
            self.assertIsNone(record.obstruction)
            self.assertLess(record.compactness_upper, 1)
            self.assertLess(record.exterior_compactness_upper, 1)
            self.assertGreater(record.compactness_margin, 0)
            self.assertEqual(record.ode_cell_count, DECLARED_BASE_CELLS)
            self.assertEqual(record.coverage_left, self.spec.support_minimum)
            self.assertEqual(record.coverage_right, self.spec.support_maximum)
            self.assertTrue(record.ode.tiles_compact_support)
            self.assertTrue(record.ode.covers_complete_support)
            self.assertEqual(record.chart_kind, PRODUCT_BOX_CHART_KIND)
            print(
                f"small-amplitude depth {record.level.max_bisection_depth}: "
                f"wall={record.diagnostic_wall_seconds:.4f}s "
                f"cells={record.ode_cell_count} "
                f"rhs={record.rhs_evaluations} "
                f"bits={record.max_rational_bit_length_observed} "
                f"C={record.compactness_upper} "
                f"margin={record.compactness_margin}",
                flush=True,
            )
        self.assertEqual(
            self.record.proved_level_names, ("depth_8", "depth_10", "depth_12")
        )
        self.assertTrue(self.record.any_declared_level_proves_continuous_no_initial_trap)
        self.assertTrue(self.record.shared_coverage_monotone)
        self.assertEqual(self.record.contract_sha256, SMALL_AMPLITUDE_CONTRACT_SHA256)
        self.assertIs(self.record.SGBL_branch_owned_and_healthy, False)
        self.assertIs(self.record.FRZ1, False)
        self.assertIs(self.record.PREF1, False)
        gate = sgbl_trap_refinement_health_gate(self.record)
        self.assertIs(gate["SGBL_branch_owned_and_healthy"], False)
        self.assertIs(gate["FRZ1"], False)
        self.assertIs(gate["PREF1"], False)
        self.assertIs(gate["holdout_authorized"], False)


class SGBLProductBoxLadderAttackTests(unittest.TestCase):
    def test_insufficient_resource_is_a_typed_resource_result_not_a_compactness_nonpass(
        self,
    ) -> None:
        tiny = replace(
            DECLARED_PRODUCT_BOX_RESOURCE_LADDER[0], max_rhs_evaluations=1
        )
        record = sgbl_evaluate_product_box_ladder_level(SGBLExactInitialSlice(), tiny)
        self.assertEqual(record.classification, "resource_limit")
        self.assertEqual(record.resource_reason, "max_rhs_evaluations")
        self.assertFalse(record.continuous_no_initial_trapped_sphere)
        self.assertFalse(record.complete_support_coverage)
        self.assertNotEqual(record.classification, "interval_inconclusive")
        self.assertIsNone(record.obstruction)
        self.assertEqual(record.missing_theorem, "declared_product_box_resource_budget")

    def test_changed_ladder_is_refused(self) -> None:
        changed = (
            replace(DECLARED_PRODUCT_BOX_RESOURCE_LADDER[0], max_cells=1),
            DECLARED_PRODUCT_BOX_RESOURCE_LADDER[1],
            DECLARED_PRODUCT_BOX_RESOURCE_LADDER[2],
        )
        with self.assertRaises(SGBLTrapRefinementStop) as stopped:
            sgbl_evaluate_product_box_resource_ladder(ladder=changed)
        self.assertEqual(stopped.exception.reason, "changed_ladder")

    def test_reordered_levels_are_refused(self) -> None:
        reordered = tuple(reversed(DECLARED_PRODUCT_BOX_RESOURCE_LADDER))
        with self.assertRaises(SGBLTrapRefinementStop) as stopped:
            sgbl_evaluate_product_box_resource_ladder(ladder=reordered)
        self.assertEqual(stopped.exception.reason, "reordered_levels")

    def test_skipped_level_is_refused(self) -> None:
        skipped = (
            DECLARED_PRODUCT_BOX_RESOURCE_LADDER[0],
            DECLARED_PRODUCT_BOX_RESOURCE_LADDER[2],
        )
        with self.assertRaises(SGBLTrapRefinementStop) as stopped:
            sgbl_evaluate_product_box_resource_ladder(ladder=skipped)
        self.assertEqual(stopped.exception.reason, "skipped_level")

    def test_dropped_or_gapped_inventory_is_refused(self) -> None:
        spec = SGBLExactInitialSlice()
        cells = (
            _dummy_cell(radius=interval(10, 11), compactness=interval(0, Q(1, 8))),
            _dummy_cell(radius=interval(11, 12), compactness=interval(0, Q(1, 4))),
            _dummy_cell(radius=interval(12, 13), compactness=interval(0, Q(1, 3))),
        )
        self.assertEqual(
            sgbl_product_box_inventory_span(cells, spec), (Q(10), Q(13))
        )
        with self.assertRaises(SGBLTrapRefinementStop) as dropped:
            sgbl_product_box_cells_cover_support(cells[1:], spec)
        self.assertEqual(dropped.exception.reason, "gapped_inventory")
        with self.assertRaises(SGBLTrapRefinementStop) as gapped:
            sgbl_product_box_inventory_span((cells[0], cells[2]), spec)
        self.assertEqual(gapped.exception.reason, "gapped_inventory")

    def test_forged_completeness_on_a_prefix_is_refused(self) -> None:
        record = sgbl_evaluate_product_box_resource_ladder(SGBLExactInitialSlice())
        prefix = record.levels[0]
        self.assertFalse(prefix.complete_support_coverage)
        self.assertLess(prefix.coverage_right, record.slice.support_maximum)
        with self.assertRaises(ValueError):
            replace(prefix, complete_support_coverage=True)
        with self.assertRaises(ValueError):
            replace(
                prefix,
                complete_support_coverage=True,
                coverage_right=record.slice.support_maximum,
            )

    def test_nonmonotone_shared_coverage_is_detected(self) -> None:
        coarser = (
            _dummy_cell(radius=interval(10, 11), compactness=interval(Q(1, 8), Q(1, 4))),
        )
        finer = (
            _dummy_cell(radius=interval(10, 11), compactness=interval(Q(1, 16), Q(1, 2))),
        )
        self.assertFalse(sgbl_product_box_shared_coverage_monotone(coarser, finer))
        nested = (
            _dummy_cell(radius=interval(10, Q(21, 2)), compactness=interval(Q(1, 8), Q(1, 5))),
        )
        self.assertTrue(sgbl_product_box_shared_coverage_monotone(coarser, nested))

    def test_tampered_pass_and_promoting_flags_are_refused(self) -> None:
        record = sgbl_evaluate_product_box_resource_ladder(SGBLExactInitialSlice())
        with self.assertRaises(TypeError):
            replace(record, SGBL_branch_owned_and_healthy=True)
        with self.assertRaises(TypeError):
            replace(record, FRZ1=True)
        with self.assertRaises(SGBLTrapRefinementStop) as stopped:
            replace(record, any_declared_level_proves_continuous_no_initial_trap=True)
        self.assertEqual(stopped.exception.reason, "tampered_contract")
        with self.assertRaises(SGBLTrapRefinementStop):
            replace(record, contract_sha256="0" * 64)
        with self.assertRaises(SGBLTrapRefinementStop):
            replace(record, shared_coverage_monotone=not record.shared_coverage_monotone)
        with self.assertRaises(ValueError):
            replace(
                record.levels[0],
                classification="continuous_no_initial_trapped_sphere",
                continuous_no_initial_trapped_sphere=True,
                obstruction=None,
                complete_support_coverage=True,
                missing_theorem=None,
            )
        with self.assertRaises(SGBLTrapRefinementStop):
            SGBLProductBoxLadderRecord(
                slice=record.slice,
                ladder=record.ladder,
                levels=record.levels,
                shared_coverage_monotone=record.shared_coverage_monotone,
                compactness_upper_nonincreasing=record.compactness_upper_nonincreasing,
                proved_level_names=record.proved_level_names,
                any_declared_level_proves_continuous_no_initial_trap=record.any_declared_level_proves_continuous_no_initial_trap,
                contract_payload=dict(record.contract_payload)
                | {"SGBL_branch_owned_and_healthy": True},
                contract_sha256=record.contract_sha256,
            )

    def test_undeclared_amplitude_is_not_a_fitted_retry(self) -> None:
        with self.assertRaises(ValueError):
            sgbl_evaluate_product_box_resource_ladder(
                SGBLExactInitialSlice(chi_amplitude=Q(1, 4))
            )


class SGBLProductBoxLadderImportTests(unittest.TestCase):
    def test_module_does_not_import_qr_health_or_shared_routing(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_trap_refinement as owner

        tree = ast.parse(Path(owner.__file__).read_text())
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imported.update(alias.name for alias in node.names)
                if node.module is not None:
                    imported.add(node.module)
            elif isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
        forbidden = {
            "proto1_compactness_upper_bound",
            "FGCQRActionParameters",
            "solve_initial_data",
            "run_fgc_gr0_calibration_v1",
            "sgbl_continuous_initial_compactness",
            "sgbl_validated_constraint_ode",
            "sgbl_validated_cj_constraint_ode",
            "build_artifact_catalog",
            "artifact-catalog",
        }
        self.assertTrue(forbidden.isdisjoint(imported))
        self.assertIn("sgbl_validate_ode_cell_inventory", imported)
        source = Path(owner.__file__).read_text()
        self.assertNotIn("GB_EFFECTIVE_DENSITY_BOUND", source)
        self.assertNotIn("1e-12", source)
        self.assertNotIn("HYP2", source)
        self.assertIn("sgbl_validated_lambda_k_constraint_ode", source)


if __name__ == "__main__":
    unittest.main()
