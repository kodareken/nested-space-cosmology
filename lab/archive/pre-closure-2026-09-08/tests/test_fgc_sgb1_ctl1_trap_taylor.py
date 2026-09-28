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
    SGBLExactInitialSlice,
    sgbl_interval_affine_coefficients,
)
from recursive_horizons.fgc.sgb1_ctl1_trap_refinement2 import (  # noqa: E402
    PREDECESSOR_NOMINAL_CONTRACT_SHA256 as LADDER1_HASH,
)
from recursive_horizons.fgc.sgb1_ctl1_trap_taylor import (  # noqa: E402
    BUMP_SCALED_DERIVATIVE_BOUNDS,
    DECLARED_COMPACTNESS_STRICT_UPPER,
    DECLARED_JET_CALLS_PER_PARTITION,
    DECLARED_TAYLOR_LEVEL_DEPTHS,
    DECLARED_TAYLOR_ORDER,
    DECLARED_TAYLOR_RESOURCE_LADDER,
    INSTRUMENT_ID,
    IntervalJet,
    NOMINAL_CHI_AMPLITUDE,
    PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256,
    PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256,
    PRODUCT_BOX_CHART_KIND,
    SGBLTaylorCellEnclosure,
    SGBLTaylorLadderRecord,
    SGBLTaylorTrapStop,
    SMALL_AMPLITUDE_CONTROL,
    sgbl_affine_constraint_jets,
    sgbl_compact_bump_interval_jet,
    sgbl_declared_taylor_cell_cap,
    sgbl_declared_taylor_jet_cap,
    sgbl_evaluate_taylor_ladder_level,
    sgbl_evaluate_taylor_resource_ladder,
    sgbl_interval_jet_shift,
    sgbl_lagrange_remainder,
    sgbl_require_complete_jet,
    sgbl_require_honest_remainder,
    sgbl_taylor_cells_cover_support,
    sgbl_taylor_contract_sha256,
    sgbl_taylor_health_gate,
    sgbl_taylor_inventory_span,
    sgbl_taylor_known_linear_ode_control,
    sgbl_taylor_model_from_jet,
    sgbl_taylor_polynomial_ode_control,
    sgbl_taylor_shared_coverage_monotone,
    sgbl_validate_taylor_cell_inventory,
)


MODULE_PATH = ROOT / "src/recursive_horizons/fgc/sgb1_ctl1_trap_taylor.py"
ZERO = interval(0)
ONE = interval(1)
DEPTH2_C_UPPER = Q(9979474004043034349, 2**63)
DEPTH2_C_MARGIN = Q(-756101967188258541, 2**63)
DEPTH4_C_UPPER = Q(9307090083172832251, 2**63)
DEPTH4_C_MARGIN = Q(-83718046318056443, 2**63)
CORRECTED_CELL_COUNTS = (11, 13)
CORRECTED_COVERAGE_RIGHTS = (Q(183, 16), Q(731, 64))
CORRECTED_JET_EVALUATIONS = (190, 250)
CORRECTED_COMPACTNESS_JETS = (80, 90)
SMALL_AMPLITUDE_C_UPPER = Q(28607003284473089, 2**64)
NOMINAL_CONTRACT_SHA256 = (
    "7de1421d200f9d276b95a7e541c0dd8b4b75eddddbe966db6611316b969f820b"
)
SMALL_AMPLITUDE_CONTRACT_SHA256 = (
    "d3cb9b3dd90c267e9851516025376eccce41453bdb6b853fdfcf49e630852d5a"
)


def _dummy_cell(*, radius, compactness, lam=None, k=None) -> SGBLTaylorCellEnclosure:
    lam = ONE if lam is None else lam
    k = ZERO if k is None else k
    return SGBLTaylorCellEnclosure(
        radius=radius,
        lambda_box=lam,
        k_box=k,
        compactness_box=compactness,
        lambda_left=lam,
        k_left=k,
        compactness_left=compactness,
        lambda_right=lam,
        k_right=k,
        compactness_right=compactness,
        compactness=compactness,
        direct_compactness=compactness,
        propagated_compactness=compactness,
        invariant_residual=ZERO,
        taylor_order=DECLARED_TAYLOR_ORDER,
        depth=0,
        taylor_strict_self_map=True,
        residual_contains_origin=True,
        compactness_invariant_contains_zero=True,
        correlated_strictly_tighter=False,
        chart_kind=PRODUCT_BOX_CHART_KIND,
        chart_coordinates=False,
    )


class SGBLTaylorPolicyFreezeTests(unittest.TestCase):
    def test_policy_is_frozen_before_any_evaluation(self) -> None:
        source = MODULE_PATH.read_text()
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
            names.index("DECLARED_TAYLOR_RESOURCE_LADDER"),
            names.index("sgbl_evaluate_taylor_resource_ladder"),
        )
        self.assertLess(
            names.index("DECLARED_TAYLOR_RESOURCE_LADDER"),
            names.index("sgbl_evaluate_taylor_ladder_level"),
        )
        self.assertLess(
            names.index("DECLARED_TAYLOR_ORDER"),
            names.index("DECLARED_TAYLOR_RESOURCE_LADDER"),
        )
        self.assertIn("sgbl_lagrange_remainder", source)
        self.assertIn("Taylor theorem", source)
        self.assertNotIn("sgbl_validated_lambda_k_constraint_ode", source)
        self.assertNotIn("sgbl_validated_constraint_ode(", source)
        self.assertNotIn("finite_difference(", source)
        self.assertNotIn("(f(x+h)", source)
        self.assertIn("sampled finite differences", source)
        self.assertIn("are not a proof", source)

    def test_declared_caps_do_not_change_family_domains_or_c_gate(self) -> None:
        self.assertEqual(DECLARED_TAYLOR_ORDER, 4)
        self.assertEqual(DECLARED_TAYLOR_LEVEL_DEPTHS, (2, 4))
        self.assertEqual(DECLARED_BASE_CELLS, 16)
        self.assertEqual(DECLARED_COMPACTNESS_STRICT_UPPER, 1)
        self.assertEqual(NOMINAL_CHI_AMPLITUDE, 3)
        self.assertEqual(SGBLExactInitialSlice().chi_amplitude, NOMINAL_CHI_AMPLITUDE)
        self.assertEqual(SGBLExactInitialSlice().support_minimum, 10)
        self.assertEqual(SGBLExactInitialSlice().support_maximum, 14)
        self.assertEqual(len(DECLARED_TAYLOR_RESOURCE_LADDER), 2)
        for depth, level in zip(DECLARED_TAYLOR_LEVEL_DEPTHS, DECLARED_TAYLOR_RESOURCE_LADDER):
            with self.subTest(depth=depth):
                self.assertEqual(level.taylor_order, DECLARED_TAYLOR_ORDER)
                self.assertEqual(level.max_bisection_depth, depth)
                self.assertEqual(level.base_cells, DECLARED_BASE_CELLS)
                self.assertEqual(level.lambda_domain, DECLARED_LAMBDA_DOMAIN)
                self.assertEqual(level.k_domain, DECLARED_K_DOMAIN)
                self.assertEqual(level.compactness_domain, DECLARED_C_DOMAIN)
                self.assertEqual(level.max_cells, sgbl_declared_taylor_cell_cap(depth))
                self.assertEqual(
                    level.max_jet_evaluations, sgbl_declared_taylor_jet_cap(depth)
                )
                self.assertGreater(level.compactness_domain.upper, 1)
                self.assertLess(level.compactness_domain.lower, 0)
        self.assertEqual(sgbl_declared_taylor_cell_cap(2), 16 * (2 ** 3 - 1))
        self.assertEqual(sgbl_declared_taylor_cell_cap(4), 16 * (2 ** 5 - 1))
        self.assertEqual(
            sgbl_declared_taylor_jet_cap(2),
            sgbl_declared_taylor_cell_cap(2)
            * DECLARED_JET_CALLS_PER_PARTITION
            * (DECLARED_TAYLOR_ORDER + 2),
        )

    def test_changed_domains_order_or_base_cells_are_refused(self) -> None:
        with self.assertRaises(ValueError):
            replace(DECLARED_TAYLOR_RESOURCE_LADDER[0], base_cells=8)
        with self.assertRaises(ValueError):
            replace(DECLARED_TAYLOR_RESOURCE_LADDER[0], taylor_order=3)
        with self.assertRaises(ValueError):
            replace(
                DECLARED_TAYLOR_RESOURCE_LADDER[0],
                compactness_domain=interval(Q(-1), Q(1)),
            )
        with self.assertRaises(ValueError):
            replace(
                DECLARED_TAYLOR_RESOURCE_LADDER[0],
                lambda_domain=interval(Q(1, 4), Q(8)),
            )

    def test_predecessor_hashes_are_regression_only_and_exact(self) -> None:
        self.assertEqual(PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256, LADDER1_HASH)
        self.assertEqual(
            PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256,
            "a82363a29b641c538a0e15f14ec2f4c99a49182a9a01c84ae7b1d14a953834c8",
        )
        self.assertEqual(
            PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256,
            "2cc8ba39a756afc433f1b04ce6fa3e26f45d795f2ba4c617ec080e7342cc396d",
        )
        source = MODULE_PATH.read_text()
        self.assertIn("predecessors_are_regression_only", source)
        self.assertNotIn("sgbl_evaluate_product_box_resource_ladder(", source)
        self.assertNotIn("sgbl_evaluate_product_box_resource_ladder2(", source)


class SGBLTaylorTheoremAndControlTests(unittest.TestCase):
    def test_polynomial_ivps_have_exact_zero_remainder(self) -> None:
        records = sgbl_taylor_polynomial_ode_control()
        self.assertEqual(len(records), DECLARED_TAYLOR_ORDER + 1)
        for record in records:
            with self.subTest(degree=record["polynomial_degree"]):
                self.assertEqual(record["taylor_order"], DECLARED_TAYLOR_ORDER)
                self.assertTrue(record["remainder_derivative_is_zero"])
                self.assertTrue(record["honest_remainder_is_zero"])
                self.assertTrue(record["exact_at_right_endpoint"])
                self.assertTrue(record["contains_exact_right_endpoint"])

    def test_known_quadratic_ode_is_reconstructed(self) -> None:
        record = sgbl_taylor_known_linear_ode_control()
        self.assertEqual(record["ode"], "y_prime_equals_two_r")
        self.assertTrue(record["remainder_derivative_is_zero"])
        self.assertTrue(record["exact_at_right_endpoint"])
        self.assertEqual(record["reconstructed"], record["exact"])

    def test_lagrange_remainder_of_a_cubic_is_exact(self) -> None:
        # f(x)=x^3, f'''=6, expansion at 0, order 2 leaves remainder 6/6 x^3=x^3.
        remainder = sgbl_lagrange_remainder(interval(6), interval(0, Q(1, 2)), 2)
        self.assertEqual(remainder, interval(0, Q(1, 8)))
        jet = IntervalJet((ZERO, ZERO, ZERO))
        model = sgbl_taylor_model_from_jet(jet, interval(6), Q(0), interval(0, Q(1, 2)))
        self.assertEqual(model.at(Q(1, 2)), interval(Q(1, 8)))

    def test_omitted_derivative_is_refused(self) -> None:
        jet = IntervalJet((ONE, interval(2)))
        with self.assertRaises(SGBLTaylorTrapStop) as stopped:
            sgbl_require_complete_jet(jet, 2)
        self.assertEqual(stopped.exception.reason, "omitted_derivative")
        with self.assertRaises(SGBLTaylorTrapStop) as shifted:
            sgbl_interval_jet_shift(jet, 1, 2)
        self.assertEqual(shifted.exception.reason, "omitted_derivative")
        square = IntervalJet.identity(interval(Q(1, 2)), 2)
        square = square * square
        self.assertEqual(square.derivatives[2], interval(2))
        with self.assertRaises(SGBLTaylorTrapStop):
            sgbl_interval_jet_shift(square.truncate(1), 2, 0)

    def test_wrong_remainder_is_refused(self) -> None:
        independent = sgbl_lagrange_remainder(interval(2), interval(0, 1), 1)
        self.assertNotEqual(independent, ZERO)
        with self.assertRaises(SGBLTaylorTrapStop) as stopped:
            sgbl_require_honest_remainder(ZERO, interval(2), interval(0, 1), 1)
        self.assertEqual(stopped.exception.reason, "wrong_remainder")
        honest = sgbl_require_honest_remainder(
            independent, interval(2), interval(0, 1), 1
        )
        self.assertEqual(honest, independent)

    def test_bump_jet_at_the_peak_and_outside_support(self) -> None:
        peak = sgbl_compact_bump_interval_jet(
            interval(12),
            center=Q(12),
            half_width=Q(2),
            order=2,
        )
        self.assertTrue(ONE.subset_of(peak.primal))
        self.assertTrue(ZERO.subset_of(peak.derivatives[1]))
        self.assertLessEqual(peak.derivatives[2].upper, 0)
        outside = sgbl_compact_bump_interval_jet(
            interval(8),
            center=Q(12),
            half_width=Q(2),
            order=3,
        )
        self.assertEqual(outside, IntervalJet.constant(0, 3))
        edge = sgbl_compact_bump_interval_jet(
            interval(10),
            center=Q(12),
            half_width=Q(2),
            order=3,
        )
        self.assertEqual(edge, IntervalJet.constant(0, 3))
        self.assertGreaterEqual(BUMP_SCALED_DERIVATIVE_BOUNDS[1], 2)
        self.assertGreaterEqual(BUMP_SCALED_DERIVATIVE_BOUNDS[2], 2)

    def test_order_zero_affine_jet_matches_the_interval_kernel(self) -> None:
        spec = SGBLExactInitialSlice()
        radius = interval(12)
        lam = interval(1)
        k = interval(0)
        from recursive_horizons.fgc.sgb1_ctl1_trap_taylor import _TaylorBudget

        jets = sgbl_affine_constraint_jets(
            IntervalJet.identity(radius, 0),
            IntervalJet.constant(lam, 0),
            IntervalJet.constant(k, 0),
            spec,
            budget=_TaylorBudget(),
            jet_limit=16,
            bit_limit=16384,
        )
        independent = sgbl_interval_affine_coefficients(
            radius=radius,
            radial_metric=lam,
            angular_extrinsic_curvature=k,
            spec=spec,
        )
        for name in (
            "hamiltonian_constant",
            "hamiltonian_lambda_r",
            "momentum_constant",
            "momentum_k_r",
        ):
            with self.subTest(name=name):
                self.assertTrue(independent[name].subset_of(jets[name].primal))


class SGBLTaylorNominalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec = SGBLExactInitialSlice()
        cls.record = sgbl_evaluate_taylor_resource_ladder(cls.spec)

    def test_every_declared_level_runs_and_does_not_stop_early(self) -> None:
        self.assertEqual(len(self.record.levels), 2)
        self.assertTrue(self.record.evaluated_every_declared_level)
        self.assertEqual(
            tuple(level.level.max_bisection_depth for level in self.record.levels),
            DECLARED_TAYLOR_LEVEL_DEPTHS,
        )
        self.assertEqual(
            tuple(level.level.taylor_order for level in self.record.levels),
            (DECLARED_TAYLOR_ORDER, DECLARED_TAYLOR_ORDER),
        )
        self.assertEqual(
            self.record.predecessor_trap_refinement_contract_sha256,
            PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256,
        )
        self.assertEqual(
            self.record.predecessor_trap_refinement2_contract_sha256,
            PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256,
        )
        for record in self.record.levels:
            self.assertEqual(record.chart_kind, PRODUCT_BOX_CHART_KIND)
            self.assertNotEqual(record.chart_kind, CERTIFICATE_CHART_KIND)
            self.assertTrue(record.product_box_is_not_the_certificate_chart)
            self.assertTrue(record.sampled_nodes_are_not_the_certificate)
            self.assertIn(
                record.classification,
                {
                    "interval_inconclusive",
                    "resource_limit",
                    "domain_error",
                    "continuous_no_initial_trapped_sphere",
                },
            )
            if record.classification == "resource_limit":
                self.assertIsNotNone(record.resource_reason)
                self.assertIsNone(record.obstruction)
            elif record.classification != "domain_error":
                self.assertGreaterEqual(record.jet_evaluations, 1)
            self.assertIsInstance(record.compactness_upper, Q)
            self.assertEqual(
                record.compactness_margin,
                1 - max(record.compactness_upper, record.exterior_compactness_upper),
            )

    def test_nominal_inventory_is_gap_free_and_health_stays_false(self) -> None:
        for record in self.record.levels:
            with self.subTest(depth=record.level.max_bisection_depth):
                if record.ode is not None and record.ode.cells:
                    left, right = sgbl_validate_taylor_cell_inventory(
                        record.ode.cells, origin=self.spec.support_minimum
                    )
                    self.assertEqual(left, record.coverage_left)
                    self.assertEqual(right, record.coverage_right)
                    self.assertEqual(record.coverage_left, self.spec.support_minimum)
                    self.assertLessEqual(record.coverage_right, self.spec.support_maximum)
                    self.assertTrue(
                        all(
                            record.ode.cells[index].radius.upper
                            == record.ode.cells[index + 1].radius.lower
                            for index in range(len(record.ode.cells) - 1)
                        )
                    )
                    self.assertTrue(
                        all(cell.residual_contains_origin for cell in record.ode.cells)
                    )
                    self.assertTrue(
                        all(
                            cell.compactness_invariant_contains_zero
                            for cell in record.ode.cells
                        )
                    )
                    self.assertEqual(
                        sgbl_taylor_inventory_span(record.ode.cells, self.spec),
                        (record.coverage_left, record.coverage_right),
                    )
                    if record.complete_support_coverage:
                        self.assertEqual(record.coverage_right, self.spec.support_maximum)
                        self.assertTrue(
                            sgbl_taylor_cells_cover_support(record.ode.cells, self.spec)
                        )
                    else:
                        self.assertFalse(
                            sgbl_taylor_cells_cover_support(record.ode.cells, self.spec)
                        )
                self.assertEqual(record.classification, "interval_inconclusive")
                self.assertEqual(record.obstruction, "compactness_enclosure_not_below_one")
                self.assertFalse(record.complete_support_coverage)
                self.assertFalse(record.continuous_no_initial_trapped_sphere)
                self.assertFalse(record.ode.tiles_compact_support)
                self.assertGreaterEqual(record.compactness_upper, 1)
                self.assertLess(record.compactness_margin, 0)
        expected_uppers = (DEPTH2_C_UPPER, DEPTH4_C_UPPER)
        expected_margins = (DEPTH2_C_MARGIN, DEPTH4_C_MARGIN)
        previous_cells = 0
        for index, record in enumerate(self.record.levels):
            with self.subTest(frozen=record.level.max_bisection_depth):
                self.assertEqual(record.ode_cell_count, CORRECTED_CELL_COUNTS[index])
                self.assertGreater(record.ode_cell_count, previous_cells)
                self.assertEqual(record.jet_evaluations, CORRECTED_JET_EVALUATIONS[index])
                self.assertEqual(
                    record.compactness_jet_evaluations, CORRECTED_COMPACTNESS_JETS[index]
                )
                self.assertEqual(record.compactness_upper, expected_uppers[index])
                self.assertEqual(record.compactness_margin, expected_margins[index])
                self.assertEqual(record.coverage_right, CORRECTED_COVERAGE_RIGHTS[index])
                previous_cells = record.ode_cell_count
        self.assertTrue(self.record.compactness_upper_nonincreasing)
        self.assertTrue(self.record.shared_coverage_monotone)
        self.assertEqual(self.record.proved_level_names, ())
        self.assertFalse(self.record.any_declared_level_proves_continuous_no_initial_trap)
        self.assertEqual(self.record.contract_sha256, NOMINAL_CONTRACT_SHA256)
        self.assertIs(self.record.SGBL_branch_owned_and_healthy, False)
        self.assertIs(self.record.FRZ1, False)
        self.assertIs(self.record.PREF1, False)
        self.assertIs(self.record.holdout_authorized, False)
        self.assertIs(self.record.execution_authorized, False)
        gate = sgbl_taylor_health_gate(self.record)
        self.assertIs(gate["SGBL_branch_owned_and_healthy"], False)
        self.assertIs(gate["FRZ1"], False)
        self.assertIs(gate["PREF1"], False)
        self.assertIs(gate["holdout_authorized"], False)
        self.assertIs(gate["copied_gr0_or_fgcqr_health_evidence"], False)
        self.assertTrue(gate["sampled_finite_differences_are_not_the_proof"])
        self.assertEqual(
            self.record.contract_sha256,
            sgbl_taylor_contract_sha256(self.record.contract_payload),
        )
        self.assertNotIn("wall", str(dict(self.record.contract_payload)))
        self.assertNotIn("diagnostic_wall_seconds", dict(self.record.contract_payload))

    def test_wall_and_resources_are_recorded_for_every_level(self) -> None:
        for record in self.record.levels:
            with self.subTest(depth=record.level.max_bisection_depth):
                self.assertGreaterEqual(record.diagnostic_wall_seconds, 0)
                self.assertLessEqual(record.ode_cell_count, record.level.max_cells)
                if record.classification != "resource_limit":
                    self.assertLessEqual(
                        record.jet_evaluations, record.level.max_jet_evaluations
                    )
                    self.assertLessEqual(
                        record.max_rational_bit_length_observed,
                        record.level.max_rational_bit_length,
                    )
                print(
                    f"taylor order {record.level.taylor_order} "
                    f"depth {record.level.max_bisection_depth}: "
                    f"wall={record.diagnostic_wall_seconds:.4f}s "
                    f"cells={record.ode_cell_count} "
                    f"jets={record.jet_evaluations} "
                    f"c_jets={record.compactness_jet_evaluations} "
                    f"bits={record.max_rational_bit_length_observed} "
                    f"C={record.compactness_upper} "
                    f"margin={record.compactness_margin} "
                    f"coverage={record.complete_support_coverage} "
                    f"span=[{record.coverage_left},{record.coverage_right}] "
                    f"class={record.classification} "
                    f"obstruction={record.obstruction} "
                    f"resource={record.resource_reason}",
                    flush=True,
                )


class SGBLTaylorPositiveControlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec = SGBLExactInitialSlice(chi_amplitude=SMALL_AMPLITUDE_CONTROL)
        cls.record = sgbl_evaluate_taylor_resource_ladder(cls.spec)

    def test_small_amplitude_keeps_aggregate_flags_false(self) -> None:
        self.assertEqual(self.spec.chi_amplitude, Q(1, 8))
        self.assertEqual(len(self.record.levels), 2)
        for record in self.record.levels:
            print(
                f"small-amplitude taylor depth {record.level.max_bisection_depth}: "
                f"wall={record.diagnostic_wall_seconds:.4f}s "
                f"cells={record.ode_cell_count} "
                f"jets={record.jet_evaluations} "
                f"C={record.compactness_upper} "
                f"margin={record.compactness_margin} "
                f"class={record.classification} "
                f"obstruction={record.obstruction}",
                flush=True,
            )
            self.assertEqual(record.classification, "continuous_no_initial_trapped_sphere")
            self.assertTrue(record.continuous_no_initial_trapped_sphere)
            self.assertTrue(record.complete_support_coverage)
            self.assertIsNone(record.obstruction)
            self.assertLess(record.compactness_upper, 1)
            self.assertGreater(record.compactness_margin, 0)
            self.assertEqual(record.ode_cell_count, DECLARED_BASE_CELLS)
            self.assertEqual(record.coverage_left, self.spec.support_minimum)
            self.assertEqual(record.coverage_right, self.spec.support_maximum)
            self.assertEqual(record.compactness_upper, SMALL_AMPLITUDE_C_UPPER)
            self.assertTrue(record.ode.tiles_compact_support)
            self.assertIs(self.record.SGBL_branch_owned_and_healthy, False)
        self.assertEqual(
            self.record.proved_level_names,
            ("order_4_depth_2", "order_4_depth_4"),
        )
        self.assertTrue(self.record.any_declared_level_proves_continuous_no_initial_trap)
        self.assertEqual(self.record.contract_sha256, SMALL_AMPLITUDE_CONTRACT_SHA256)
        gate = sgbl_taylor_health_gate(self.record)
        self.assertIs(gate["FRZ1"], False)
        self.assertIs(gate["PREF1"], False)
        self.assertIs(gate["holdout_authorized"], False)


class SGBLTaylorAttackTests(unittest.TestCase):
    def test_insufficient_resource_is_a_typed_resource_result(self) -> None:
        tiny = replace(DECLARED_TAYLOR_RESOURCE_LADDER[0], max_jet_evaluations=1)
        record = sgbl_evaluate_taylor_ladder_level(SGBLExactInitialSlice(), tiny)
        self.assertEqual(record.classification, "resource_limit")
        self.assertEqual(record.resource_reason, "max_jet_evaluations")
        self.assertFalse(record.continuous_no_initial_trapped_sphere)
        self.assertFalse(record.complete_support_coverage)
        self.assertNotEqual(record.classification, "interval_inconclusive")
        self.assertIsNone(record.obstruction)
        self.assertEqual(record.missing_theorem, "declared_taylor_model_resource_budget")

    def test_changed_policy_is_refused(self) -> None:
        changed = (
            replace(DECLARED_TAYLOR_RESOURCE_LADDER[0], max_cells=1),
            DECLARED_TAYLOR_RESOURCE_LADDER[1],
        )
        with self.assertRaises(SGBLTaylorTrapStop) as stopped:
            sgbl_evaluate_taylor_resource_ladder(ladder=changed)
        self.assertEqual(stopped.exception.reason, "changed_policy")

    def test_reordered_levels_are_refused(self) -> None:
        reordered = tuple(reversed(DECLARED_TAYLOR_RESOURCE_LADDER))
        with self.assertRaises(SGBLTaylorTrapStop) as stopped:
            sgbl_evaluate_taylor_resource_ladder(ladder=reordered)
        self.assertEqual(stopped.exception.reason, "reordered_levels")

    def test_skipped_level_is_refused(self) -> None:
        skipped = (DECLARED_TAYLOR_RESOURCE_LADDER[0],)
        with self.assertRaises(SGBLTaylorTrapStop) as stopped:
            sgbl_evaluate_taylor_resource_ladder(ladder=skipped)
        self.assertEqual(stopped.exception.reason, "skipped_level")

    def test_wrong_predecessor_hash_is_refused(self) -> None:
        with self.assertRaises(SGBLTaylorTrapStop) as stopped:
            sgbl_evaluate_taylor_resource_ladder(
                predecessor_trap_refinement_contract_sha256="0" * 64
            )
        self.assertEqual(stopped.exception.reason, "wrong_predecessor")
        with self.assertRaises(SGBLTaylorTrapStop) as stopped2:
            sgbl_evaluate_taylor_resource_ladder(
                predecessor_trap_refinement2_contract_sha256="1" * 64
            )
        self.assertEqual(stopped2.exception.reason, "wrong_predecessor")

    def test_dropped_or_gapped_inventory_is_refused(self) -> None:
        spec = SGBLExactInitialSlice()
        cells = (
            _dummy_cell(radius=interval(10, 11), compactness=interval(0, Q(1, 8))),
            _dummy_cell(radius=interval(11, 12), compactness=interval(0, Q(1, 4))),
            _dummy_cell(radius=interval(12, 13), compactness=interval(0, Q(1, 3))),
        )
        self.assertEqual(sgbl_taylor_inventory_span(cells, spec), (Q(10), Q(13)))
        with self.assertRaises(SGBLTaylorTrapStop) as dropped:
            sgbl_taylor_cells_cover_support(cells[1:], spec)
        self.assertEqual(dropped.exception.reason, "gapped_inventory")
        with self.assertRaises(SGBLTaylorTrapStop) as gapped:
            sgbl_taylor_inventory_span((cells[0], cells[2]), spec)
        self.assertEqual(gapped.exception.reason, "gapped_inventory")

    def test_nonmonotone_shared_coverage_is_detected(self) -> None:
        coarser = (
            _dummy_cell(
                radius=interval(10, 11), compactness=interval(Q(1, 8), Q(1, 4))
            ),
        )
        finer = (
            _dummy_cell(
                radius=interval(10, 11), compactness=interval(Q(1, 16), Q(1, 2))
            ),
        )
        self.assertFalse(sgbl_taylor_shared_coverage_monotone(coarser, finer))
        nested = (
            _dummy_cell(
                radius=interval(10, Q(21, 2)), compactness=interval(Q(1, 8), Q(1, 5))
            ),
        )
        self.assertTrue(sgbl_taylor_shared_coverage_monotone(coarser, nested))

    def test_tampered_pass_and_promoting_flags_are_refused(self) -> None:
        record = sgbl_evaluate_taylor_resource_ladder(SGBLExactInitialSlice())
        with self.assertRaises(TypeError):
            replace(record, SGBL_branch_owned_and_healthy=True)
        with self.assertRaises(TypeError):
            replace(record, FRZ1=True)
        with self.assertRaises(SGBLTaylorTrapStop) as stopped:
            replace(record, any_declared_level_proves_continuous_no_initial_trap=True)
        self.assertEqual(stopped.exception.reason, "tampered_contract")
        with self.assertRaises(SGBLTaylorTrapStop):
            replace(record, contract_sha256="0" * 64)
        with self.assertRaises(SGBLTaylorTrapStop):
            SGBLTaylorLadderRecord(
                slice=record.slice,
                ladder=record.ladder,
                levels=record.levels,
                predecessor_trap_refinement_contract_sha256=(
                    record.predecessor_trap_refinement_contract_sha256
                ),
                predecessor_trap_refinement2_contract_sha256=(
                    record.predecessor_trap_refinement2_contract_sha256
                ),
                shared_coverage_monotone=record.shared_coverage_monotone,
                compactness_upper_nonincreasing=record.compactness_upper_nonincreasing,
                proved_level_names=record.proved_level_names,
                any_declared_level_proves_continuous_no_initial_trap=(
                    record.any_declared_level_proves_continuous_no_initial_trap
                ),
                contract_payload=dict(record.contract_payload)
                | {"SGBL_branch_owned_and_healthy": True},
                contract_sha256=record.contract_sha256,
            )

    def test_undeclared_amplitude_is_not_a_fitted_retry(self) -> None:
        with self.assertRaises(ValueError):
            sgbl_evaluate_taylor_resource_ladder(
                SGBLExactInitialSlice(chi_amplitude=Q(1, 4))
            )


class SGBLTaylorImportTests(unittest.TestCase):
    def test_module_does_not_import_qr_health_or_picard_ode(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_trap_taylor as owner

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
            "sgbl_validated_lambda_k_constraint_ode",
            "sgbl_validated_cj_constraint_ode",
            "build_artifact_catalog",
            "artifact-catalog",
            "numpy",
        }
        self.assertTrue(forbidden.isdisjoint(imported))
        source = Path(owner.__file__).read_text()
        self.assertNotIn("GB_EFFECTIVE_DENSITY_BOUND", source)
        self.assertNotIn("1e-12", source)
        self.assertIn("PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256", source)
        self.assertEqual(INSTRUMENT_ID, "FGC-1-SGB1-CTL1-TRAP-TAYLOR")


if __name__ == "__main__":
    unittest.main()
