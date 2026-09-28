from __future__ import annotations

import ast
from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.sgb1_ctl1_initial_health import (  # noqa: E402
    CERTIFICATE_CHART_KIND,
    DECLARED_BASE_CELLS,
    DECLARED_C_DOMAIN,
    DECLARED_K_DOMAIN,
    DECLARED_LAMBDA_DOMAIN,
    DECLARED_MAX_PICARD_ITERATIONS,
    SGBLExactInitialSlice,
    sgbl_validate_ode_cell_inventory,
)
from recursive_horizons.fgc.sgb1_ctl1_trap_refinement import (  # noqa: E402
    DECLARED_LADDER_DEPTHS,
    DECLARED_PRODUCT_BOX_RESOURCE_LADDER,
    INSTRUMENT_ID as PREDECESSOR_ID,
    SMALL_AMPLITUDE_CONTROL,
    sgbl_evaluate_product_box_ladder_level,
    sgbl_evaluate_product_box_resource_ladder,
)
from recursive_horizons.fgc.sgb1_ctl1_trap_refinement2 import (  # noqa: E402
    DECLARED_LADDER2_BIT_CAP,
    DECLARED_LADDER2_CELL_CAPS,
    DECLARED_LADDER2_DEPTHS,
    DECLARED_LADDER2_RHS_CAPS,
    DECLARED_PRODUCT_BOX_RESOURCE_LADDER2,
    INSTRUMENT_ID,
    PREDECESSOR_INSTRUMENT_ID,
    PREDECESSOR_NOMINAL_CONTRACT_SHA256,
    SGBLProductBoxLadder2Record,
    SGBLTrapRefinement2Stop,
    sgbl_evaluate_product_box_resource_ladder2,
    sgbl_trap_refinement2_health_gate,
)


DEPTH14_C_UPPER = Q(2305843524457714205, 2**61)
DEPTH14_C_MARGIN = Q(-515244020253, 2**61)
DEPTH16_C_UPPER = Q(18446753041991488113, 2**64)
DEPTH16_C_MARGIN = Q(-8968281936497, 2**64)
DEPTH18_C_UPPER = Q(18446754249224814419, 2**64)
DEPTH18_C_MARGIN = Q(-10175515262803, 2**64)
CELL_COUNTS = (28, 32, 37)
COVERAGE_RIGHTS = (Q(372247, 32768), Q(744495, 65536), Q(2977981, 262144))
RHS_EVALUATIONS = (584, 648, 736)
COMPACTNESS_RHS = (146, 162, 184)
NOMINAL_CONTRACT_SHA256 = (
    "2cc8ba39a756afc433f1b04ce6fa3e26f45d795f2ba4c617ec080e7342cc396d"
)
SMALL_AMPLITUDE_CONTRACT_SHA256 = (
    "9a1f2c91b2b3682b88cdad9aedf121ee5bbcd35a13d1d8bfdd0443bd8cb4df1a"
)


class SGBLProductBoxLadder2FreezeTests(unittest.TestCase):
    def test_successor_ladder_is_frozen_before_evaluation(self) -> None:
        source = Path(
            ROOT / "src/recursive_horizons/fgc/sgb1_ctl1_trap_refinement2.py"
        ).read_text()
        tree = ast.parse(source)
        names: list[str] = []
        for node in tree.body:
            if isinstance(node, ast.Assign):
                names.extend(
                    target.id for target in node.targets if isinstance(target, ast.Name)
                )
            elif isinstance(node, ast.FunctionDef):
                names.append(node.name)
        self.assertLess(
            names.index("DECLARED_PRODUCT_BOX_RESOURCE_LADDER2"),
            names.index("sgbl_evaluate_product_box_resource_ladder2"),
        )
        self.assertIn(PREDECESSOR_NOMINAL_CONTRACT_SHA256, source)
        self.assertNotIn("sgbl_evaluate_product_box_resource_ladder(", source)

    def test_practical_caps_are_predeclared_and_do_not_change_family_or_c(self) -> None:
        self.assertEqual(DECLARED_LADDER2_DEPTHS, (14, 16, 18))
        self.assertEqual(DECLARED_LADDER2_CELL_CAPS, (4096, 8192, 16384))
        self.assertEqual(DECLARED_LADDER2_RHS_CAPS, (131072, 262144, 524288))
        self.assertEqual(DECLARED_LADDER2_BIT_CAP, 65536)
        self.assertEqual(DECLARED_LADDER_DEPTHS, (8, 10, 12))
        self.assertEqual(PREDECESSOR_INSTRUMENT_ID, PREDECESSOR_ID)
        self.assertNotEqual(INSTRUMENT_ID, PREDECESSOR_ID)
        self.assertEqual(len(DECLARED_PRODUCT_BOX_RESOURCE_LADDER2), 3)
        self.assertEqual(SGBLExactInitialSlice().chi_amplitude, 3)
        for depth, cells, rhs, level in zip(
            DECLARED_LADDER2_DEPTHS,
            DECLARED_LADDER2_CELL_CAPS,
            DECLARED_LADDER2_RHS_CAPS,
            DECLARED_PRODUCT_BOX_RESOURCE_LADDER2,
        ):
            with self.subTest(depth=depth):
                self.assertEqual(level.max_bisection_depth, depth)
                self.assertEqual(level.base_cells, DECLARED_BASE_CELLS)
                self.assertEqual(level.max_picard_iterations, DECLARED_MAX_PICARD_ITERATIONS)
                self.assertEqual(level.lambda_domain, DECLARED_LAMBDA_DOMAIN)
                self.assertEqual(level.k_domain, DECLARED_K_DOMAIN)
                self.assertEqual(level.compactness_domain, DECLARED_C_DOMAIN)
                self.assertEqual(level.max_cells, cells)
                self.assertEqual(level.max_rhs_evaluations, rhs)
                self.assertEqual(level.max_rational_bit_length, DECLARED_LADDER2_BIT_CAP)
                self.assertGreater(level.compactness_domain.upper, 1)
                self.assertLess(level.compactness_domain.lower, 0)
        self.assertNotEqual(
            DECLARED_PRODUCT_BOX_RESOURCE_LADDER2,
            DECLARED_PRODUCT_BOX_RESOURCE_LADDER,
        )

    def test_predecessor_nominal_hash_is_the_preserved_depth_eight_to_twelve_contract(
        self,
    ) -> None:
        predecessor = sgbl_evaluate_product_box_resource_ladder(SGBLExactInitialSlice())
        self.assertEqual(predecessor.contract_sha256, PREDECESSOR_NOMINAL_CONTRACT_SHA256)
        self.assertFalse(predecessor.any_declared_level_proves_continuous_no_initial_trap)
        self.assertEqual(tuple(level.level.max_bisection_depth for level in predecessor.levels), (8, 10, 12))


class SGBLProductBoxLadder2NominalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec = SGBLExactInitialSlice()
        cls.record = sgbl_evaluate_product_box_resource_ladder2(cls.spec)

    def test_every_declared_successor_level_runs(self) -> None:
        self.assertEqual(len(self.record.levels), 3)
        self.assertTrue(self.record.evaluated_every_declared_level)
        self.assertEqual(
            tuple(level.level.max_bisection_depth for level in self.record.levels),
            DECLARED_LADDER2_DEPTHS,
        )
        self.assertEqual(
            self.record.predecessor_nominal_contract_sha256,
            PREDECESSOR_NOMINAL_CONTRACT_SHA256,
        )

    def test_nominal_levels_remain_inconclusive_on_gap_free_prefixes(self) -> None:
        expected_uppers = (DEPTH14_C_UPPER, DEPTH16_C_UPPER, DEPTH18_C_UPPER)
        expected_margins = (DEPTH14_C_MARGIN, DEPTH16_C_MARGIN, DEPTH18_C_MARGIN)
        for index, record in enumerate(self.record.levels):
            with self.subTest(depth=record.level.max_bisection_depth):
                self.assertEqual(record.classification, "interval_inconclusive")
                self.assertEqual(record.obstruction, "compactness_enclosure_not_below_one")
                self.assertFalse(record.complete_support_coverage)
                self.assertFalse(record.continuous_no_initial_trapped_sphere)
                self.assertEqual(record.ode_cell_count, CELL_COUNTS[index])
                self.assertEqual(record.rhs_evaluations, RHS_EVALUATIONS[index])
                self.assertEqual(record.compactness_rhs_evaluations, COMPACTNESS_RHS[index])
                self.assertEqual(record.compactness_upper, expected_uppers[index])
                self.assertEqual(record.compactness_margin, expected_margins[index])
                self.assertGreaterEqual(record.compactness_upper, 1)
                self.assertLess(record.compactness_margin, 0)
                self.assertEqual(record.coverage_left, self.spec.support_minimum)
                self.assertEqual(record.coverage_right, COVERAGE_RIGHTS[index])
                self.assertLess(record.coverage_right, self.spec.support_maximum)
                left, right = sgbl_validate_ode_cell_inventory(
                    record.ode.cells, origin=self.spec.support_minimum
                )
                self.assertEqual((left, right), (record.coverage_left, record.coverage_right))
                self.assertNotEqual(record.chart_kind, CERTIFICATE_CHART_KIND)
        self.assertEqual(self.record.proved_level_names, ())
        self.assertEqual(self.record.tiled_support_level_names, ())
        self.assertEqual(self.record.tiled_support_with_c_lt_1_level_names, ())
        self.assertFalse(self.record.any_declared_level_proves_continuous_no_initial_trap)
        self.assertFalse(self.record.any_declared_level_tiles_compact_support)
        self.assertFalse(self.record.any_declared_level_tiles_and_c_strictly_below_one)
        self.assertTrue(self.record.shared_coverage_monotone)
        self.assertFalse(self.record.compactness_upper_nonincreasing)
        self.assertEqual(self.record.contract_sha256, NOMINAL_CONTRACT_SHA256)

    def test_wall_and_resources_are_recorded_for_every_successor_level(self) -> None:
        for record in self.record.levels:
            with self.subTest(depth=record.level.max_bisection_depth):
                self.assertGreaterEqual(record.diagnostic_wall_seconds, 0)
                self.assertLessEqual(record.ode_cell_count, record.level.max_cells)
                self.assertLessEqual(record.rhs_evaluations, record.level.max_rhs_evaluations)
                self.assertLessEqual(
                    record.max_rational_bit_length_observed,
                    record.level.max_rational_bit_length,
                )
                print(
                    f"successor depth {record.level.max_bisection_depth}: "
                    f"wall={record.diagnostic_wall_seconds:.4f}s "
                    f"cells={record.ode_cell_count} "
                    f"rhs={record.rhs_evaluations} "
                    f"c_rhs={record.compactness_rhs_evaluations} "
                    f"bits={record.max_rational_bit_length_observed} "
                    f"C={record.compactness_upper} "
                    f"margin={record.compactness_margin} "
                    f"tile={record.complete_support_coverage} "
                    f"span=[{record.coverage_left},{record.coverage_right}] "
                    f"obstruction={record.obstruction}",
                    flush=True,
                )

    def test_aggregate_health_stays_false(self) -> None:
        gate = sgbl_trap_refinement2_health_gate(self.record)
        self.assertIs(self.record.SGBL_branch_owned_and_healthy, False)
        self.assertIs(self.record.FRZ1, False)
        self.assertIs(self.record.PREF1, False)
        self.assertIs(self.record.holdout_authorized, False)
        self.assertIs(gate["SGBL_branch_owned_and_healthy"], False)
        self.assertIs(gate["FRZ1"], False)
        self.assertIs(gate["PREF1"], False)
        self.assertEqual(
            self.record.contract_payload["predecessor_nominal_contract_sha256"],
            PREDECESSOR_NOMINAL_CONTRACT_SHA256,
        )
        self.assertNotIn("diagnostic_wall_seconds", dict(self.record.contract_payload))


class SGBLProductBoxLadder2PositiveControlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec = SGBLExactInitialSlice(chi_amplitude=SMALL_AMPLITUDE_CONTROL)
        cls.record = sgbl_evaluate_product_box_resource_ladder2(cls.spec)

    def test_small_amplitude_passes_every_successor_level_without_stopping_early(self) -> None:
        self.assertEqual(len(self.record.levels), 3)
        for record in self.record.levels:
            self.assertEqual(record.classification, "continuous_no_initial_trapped_sphere")
            self.assertTrue(record.complete_support_coverage)
            self.assertEqual(record.ode_cell_count, DECLARED_BASE_CELLS)
            self.assertEqual(record.coverage_left, self.spec.support_minimum)
            self.assertEqual(record.coverage_right, self.spec.support_maximum)
            self.assertLess(record.compactness_upper, 1)
            print(
                f"successor small-amplitude depth {record.level.max_bisection_depth}: "
                f"wall={record.diagnostic_wall_seconds:.4f}s "
                f"cells={record.ode_cell_count} "
                f"rhs={record.rhs_evaluations} "
                f"C={record.compactness_upper} "
                f"margin={record.compactness_margin}",
                flush=True,
            )
        self.assertEqual(
            self.record.proved_level_names, ("depth_14", "depth_16", "depth_18")
        )
        self.assertTrue(self.record.any_declared_level_tiles_and_c_strictly_below_one)
        self.assertIs(self.record.SGBL_branch_owned_and_healthy, False)
        self.assertIs(self.record.FRZ1, False)
        self.assertEqual(self.record.contract_sha256, SMALL_AMPLITUDE_CONTRACT_SHA256)
        gate = sgbl_trap_refinement2_health_gate(self.record)
        self.assertIs(gate["PREF1"], False)
        self.assertIs(gate["holdout_authorized"], False)


class SGBLProductBoxLadder2AttackTests(unittest.TestCase):
    def test_wrong_predecessor_hash_is_refused(self) -> None:
        with self.assertRaises(SGBLTrapRefinement2Stop) as stopped:
            sgbl_evaluate_product_box_resource_ladder2(
                predecessor_nominal_contract_sha256="0" * 64
            )
        self.assertEqual(stopped.exception.reason, "wrong_predecessor")

    def test_changed_resource_caps_are_refused(self) -> None:
        changed = (
            replace(DECLARED_PRODUCT_BOX_RESOURCE_LADDER2[0], max_cells=1),
            DECLARED_PRODUCT_BOX_RESOURCE_LADDER2[1],
            DECLARED_PRODUCT_BOX_RESOURCE_LADDER2[2],
        )
        with self.assertRaises(SGBLTrapRefinement2Stop) as stopped:
            sgbl_evaluate_product_box_resource_ladder2(ladder=changed)
        self.assertEqual(stopped.exception.reason, "changed_ladder")

    def test_insufficient_resource_is_typed_resource_limit(self) -> None:
        tiny = replace(DECLARED_PRODUCT_BOX_RESOURCE_LADDER2[0], max_rhs_evaluations=1)
        record = sgbl_evaluate_product_box_ladder_level(SGBLExactInitialSlice(), tiny)
        self.assertEqual(record.classification, "resource_limit")
        self.assertFalse(record.continuous_no_initial_trapped_sphere)
        self.assertFalse(record.complete_support_coverage)

    def test_reordered_and_skipped_levels_are_refused(self) -> None:
        with self.assertRaises(SGBLTrapRefinement2Stop) as reordered:
            sgbl_evaluate_product_box_resource_ladder2(
                ladder=tuple(reversed(DECLARED_PRODUCT_BOX_RESOURCE_LADDER2))
            )
        self.assertEqual(reordered.exception.reason, "reordered_levels")
        with self.assertRaises(SGBLTrapRefinement2Stop) as skipped:
            sgbl_evaluate_product_box_resource_ladder2(
                ladder=(
                    DECLARED_PRODUCT_BOX_RESOURCE_LADDER2[0],
                    DECLARED_PRODUCT_BOX_RESOURCE_LADDER2[2],
                )
            )
        self.assertEqual(skipped.exception.reason, "skipped_level")

    def test_pass_promotion_and_forged_tiling_are_refused(self) -> None:
        record = sgbl_evaluate_product_box_resource_ladder2(SGBLExactInitialSlice())
        with self.assertRaises(TypeError):
            replace(record, SGBL_branch_owned_and_healthy=True)
        with self.assertRaises(SGBLTrapRefinement2Stop) as stopped:
            replace(record, any_declared_level_proves_continuous_no_initial_trap=True)
        self.assertEqual(stopped.exception.reason, "tampered_contract")
        with self.assertRaises(SGBLTrapRefinement2Stop):
            replace(record, any_declared_level_tiles_compact_support=True)
        with self.assertRaises(SGBLTrapRefinement2Stop):
            replace(record, contract_sha256="0" * 64)
        with self.assertRaises(ValueError):
            replace(
                record.levels[0],
                classification="continuous_no_initial_trapped_sphere",
                continuous_no_initial_trapped_sphere=True,
                obstruction=None,
                complete_support_coverage=True,
                missing_theorem=None,
            )
        with self.assertRaises(SGBLTrapRefinement2Stop):
            SGBLProductBoxLadder2Record(
                slice=record.slice,
                ladder=record.ladder,
                levels=record.levels,
                predecessor_nominal_contract_sha256=record.predecessor_nominal_contract_sha256,
                shared_coverage_monotone=record.shared_coverage_monotone,
                compactness_upper_nonincreasing=record.compactness_upper_nonincreasing,
                proved_level_names=record.proved_level_names,
                tiled_support_level_names=record.tiled_support_level_names,
                tiled_support_with_c_lt_1_level_names=record.tiled_support_with_c_lt_1_level_names,
                any_declared_level_proves_continuous_no_initial_trap=record.any_declared_level_proves_continuous_no_initial_trap,
                any_declared_level_tiles_compact_support=record.any_declared_level_tiles_compact_support,
                any_declared_level_tiles_and_c_strictly_below_one=record.any_declared_level_tiles_and_c_strictly_below_one,
                contract_payload=dict(record.contract_payload) | {"FRZ1": True},
                contract_sha256=record.contract_sha256,
            )

    def test_undeclared_amplitude_is_not_a_fitted_retry(self) -> None:
        with self.assertRaises(ValueError):
            sgbl_evaluate_product_box_resource_ladder2(
                SGBLExactInitialSlice(chi_amplitude=Q(1, 4))
            )


class SGBLProductBoxLadder2ImportTests(unittest.TestCase):
    def test_module_does_not_import_qr_health_or_shared_routing(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_trap_refinement2 as owner

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
            "build_artifact_catalog",
            "FGCQRActionParameters",
            "sgbl_continuous_initial_compactness",
            "sgbl_validated_constraint_ode",
            "run_fgc_gr0_calibration_v1",
        }
        self.assertTrue(forbidden.isdisjoint(imported))
        source = Path(owner.__file__).read_text()
        self.assertNotIn("GB_EFFECTIVE_DENSITY_BOUND", source)
        self.assertNotIn("1e-12", source)
        self.assertIn("sgbl_evaluate_product_box_ladder_level", source)


if __name__ == "__main__":
    unittest.main()
