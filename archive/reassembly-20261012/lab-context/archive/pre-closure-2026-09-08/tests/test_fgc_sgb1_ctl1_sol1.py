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
from recursive_horizons.fgc.sgb1_ctl1_cell_admission import (  # noqa: E402
    sgbl_cell_admission_health_gate,
    sgbl_whole_cell_admission_from_health,
)
from recursive_horizons.fgc.sgb1_ctl1_family import (  # noqa: E402
    SGBLFamilyParameters,
    sgbl_branch_health_gate,
)
from recursive_horizons.fgc.sgb1_ctl1_family_principal import (  # noqa: E402
    sgbl_family_principal_from_health,
    sgbl_family_principal_health_gate,
)
from recursive_horizons.fgc.sgb1_ctl1_initial_health import (  # noqa: E402
    CERTIFICATE_CHART_KIND,
    SGBLExactInitialSlice,
    SGBLInitialHealthStop,
    sgbl_continuous_initial_compactness,
    sgbl_constraint_rhs_state_jacobian,
    sgbl_interval_affine_coefficients,
    sgbl_validated_constraint_ode,
)
from recursive_horizons.fgc.sgb1_ctl1_sol1 import (  # noqa: E402
    DECLARED_CONTINUATION_STEPS,
    DECLARED_ORBIT_TUBE_POLICY,
    DECLARED_STEP_WIDTH,
    DECLARED_TUBE_RADIUS,
    FALSE_CLAIMS,
    FORBIDDEN_HEALTH_IMPORTS,
    INSTRUMENT_ID,
    NOMINAL_CHI_AMPLITUDE,
    PREDECESSOR_TRAP_BARRIER_CONTRACT_SHA256,
    PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256,
    PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256,
    PREDECESSOR_TRAP_TAYLOR_CONTRACT_SHA256,
    PREFIX_C_MARGIN,
    PREFIX_C_UPPER,
    PREFIX_CELL_COUNT,
    PREFIX_D_MARGIN,
    PREFIX_LAST_CELL_LEFT,
    PREFIX_OBSTRUCTION,
    PREFIX_RADIUS,
    SGBLOrbitLocalRecord,
    SGBLOrbitLocalStop,
    SGBLOrbitTubeEnclosure,
    SMALL_AMPLITUDE_CONTROL,
    sgbl_authenticate_ck_prefix,
    sgbl_bind_orbit_local_picard_lindelof,
    sgbl_classify_on_orbit_compactness,
    sgbl_declared_orbit_step_grid,
    sgbl_interval_matrix_infinity_norm_2x2,
    sgbl_orbit_inventory_span,
    sgbl_orbit_local_contract_sha256,
    sgbl_orbit_local_health_gate,
    sgbl_orbit_local_minkowski_buffer_control,
    sgbl_orbit_local_small_amplitude_control,
    sgbl_orbit_tube_around,
    sgbl_orbit_tube_picard_lindelof,
    sgbl_picard_lindelof_constant_ode_control,
    sgbl_require_affine_jacobian_diagonals,
    sgbl_validate_orbit_tube_inventory,
)


MODULE_PATH = ROOT / "src/recursive_horizons/fgc/sgb1_ctl1_sol1.py"
NOMINAL_CONTRACT_SHA256 = (
    "5fc44f7023eb82632309c6c15120f2e62e5e6483911e18589579576b36654f4b"
)
SMALL_AMPLITUDE_CONTRACT_SHA256 = (
    "661f06f9e6d6fa91d64f1767cf8b8b810c6ce17a4e9bd7ccf8a293b72ed946b9"
)
ZERO = interval(0)
ONE = interval(1)


def _dummy_tube(*, radius, lam=None, k=None) -> SGBLOrbitTubeEnclosure:
    lam = ONE if lam is None else lam
    k = ZERO if k is None else k
    tube = sgbl_orbit_tube_around(lam, DECLARED_TUBE_RADIUS)
    k_tube = sgbl_orbit_tube_around(k, DECLARED_TUBE_RADIUS)
    compactness = ZERO
    return SGBLOrbitTubeEnclosure(
        radius=radius,
        lambda_left=lam,
        k_left=k,
        lambda_tube=tube,
        k_tube=k_tube,
        lambda_image=lam,
        k_image=k,
        lambda_right=lam,
        k_right=k,
        compactness=compactness,
        direct_compactness=compactness,
        invariant_residual=ZERO,
        hamiltonian_lambda_r=ONE,
        momentum_k_r=ONE,
        lipschitz=Q(0),
        max_speed=Q(0),
        contraction_constant=Q(0),
        initial_inside_tube=True,
        picard_self_map=True,
        jacobian_diagonals_exclude_zero=True,
        residual_contains_origin=True,
        compactness_invariant_contains_zero=True,
        unique_local_affine_orbit=True,
        classification="unique_local_affine_orbit",
        obstruction=None,
        missing_theorem=None,
        rhs_evaluations=0,
        compactness_rhs_evaluations=0,
        max_rational_bit_length_observed=1,
        depth=0,
    )


class SGBLOrbitLocalPolicyFreezeTests(unittest.TestCase):
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
        self.assertLess(names.index("PREFIX_RADIUS"), names.index("sgbl_bind_orbit_local_picard_lindelof"))
        self.assertLess(
            names.index("DECLARED_TUBE_RADIUS"),
            names.index("sgbl_orbit_tube_picard_lindelof"),
        )
        self.assertLess(
            names.index("PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256"),
            names.index("sgbl_bind_orbit_local_picard_lindelof"),
        )
        self.assertEqual(INSTRUMENT_ID, "FGC-1-SGB1-CTL1-SOL1")
        self.assertEqual(NOMINAL_CHI_AMPLITUDE, 3)
        self.assertEqual(SMALL_AMPLITUDE_CONTROL, Q(1, 8))
        self.assertEqual(PREFIX_RADIUS, Q(10749, 1024))
        self.assertEqual(PREFIX_LAST_CELL_LEFT, Q(2687, 256))
        self.assertEqual(PREFIX_CELL_COUNT, 9)
        self.assertEqual(PREFIX_OBSTRUCTION, "picard_strict_self_map_failed")
        self.assertEqual(PREFIX_C_UPPER, Q(17053332284082708107, 2**64))
        self.assertEqual(PREFIX_C_MARGIN, 1 - PREFIX_C_UPPER)
        self.assertEqual(PREFIX_D_MARGIN, Q(326806051576253549, 2**62))
        self.assertEqual(DECLARED_CONTINUATION_STEPS, 16)
        self.assertEqual(DECLARED_TUBE_RADIUS, Q(1, 8))
        self.assertEqual(DECLARED_STEP_WIDTH, Q(3587, 16384))
        self.assertEqual(DECLARED_ORBIT_TUBE_POLICY.max_bisection_depth, 0)
        self.assertEqual(DECLARED_ORBIT_TUBE_POLICY.max_rational_bit_length, 16384)
        grid = sgbl_declared_orbit_step_grid()
        self.assertEqual(len(grid), 17)
        self.assertEqual(grid[0], PREFIX_RADIUS)
        self.assertEqual(grid[-1], 14)
        self.assertEqual(grid[1] - grid[0], DECLARED_STEP_WIDTH)

    def test_predecessor_hashes_are_the_frozen_no_trap_identities(self) -> None:
        self.assertEqual(
            PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256,
            "a82363a29b641c538a0e15f14ec2f4c99a49182a9a01c84ae7b1d14a953834c8",
        )
        self.assertEqual(
            PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256,
            "2cc8ba39a756afc433f1b04ce6fa3e26f45d795f2ba4c617ec080e7342cc396d",
        )
        self.assertEqual(
            PREDECESSOR_TRAP_TAYLOR_CONTRACT_SHA256,
            "7de1421d200f9d276b95a7e541c0dd8b4b75eddddbe966db6611316b969f820b",
        )
        self.assertEqual(
            PREDECESSOR_TRAP_BARRIER_CONTRACT_SHA256,
            "a75d564c0ea0658802c138e0575edf386375ebba35728661bd8a11c579e7d689",
        )


class SGBLOrbitLocalControlTests(unittest.TestCase):
    def test_constant_ode_is_a_unique_untrapped_orbit(self) -> None:
        tube = sgbl_picard_lindelof_constant_ode_control()
        self.assertTrue(tube.unique_local_affine_orbit)
        self.assertEqual(tube.classification, "unique_local_affine_orbit")
        self.assertTrue(tube.picard_self_map)
        self.assertEqual(tube.lipschitz, 0)
        self.assertEqual(tube.compactness, ZERO)
        self.assertLess(tube.compactness.upper, 1)
        self.assertEqual(tube.depth, 0)

    def test_minkowski_buffer_reduction_remains_a_unique_orbit_control(self) -> None:
        tube = sgbl_orbit_local_minkowski_buffer_control()
        self.assertTrue(tube.unique_local_affine_orbit)
        self.assertEqual(tube.radius.lower, 10)
        self.assertEqual(tube.lambda_left, ONE)
        self.assertEqual(tube.k_left, ZERO)
        self.assertEqual(tube.compactness, ZERO)

    def test_jacobian_zero_injection_is_an_on_orbit_jacobian_nonpass(self) -> None:
        with self.assertRaises(SGBLInitialHealthStop) as stopped:
            sgbl_require_affine_jacobian_diagonals(interval(-1, 1), interval(1, 2))
        self.assertEqual(stopped.exception.reason, "interval_inconclusive")
        self.assertEqual(
            stopped.exception.payload["obstruction"], "jacobian_diagonal_contains_zero"
        )
        self.assertEqual(
            sgbl_classify_on_orbit_compactness(interval(0, Q(3, 2))),
            "on_orbit_compactness_not_below_one",
        )
        self.assertEqual(
            sgbl_classify_on_orbit_compactness(interval(0, Q(1, 2))),
            "unique_local_affine_orbit",
        )

    def test_lipschitz_infinity_norm_is_the_row_sum_bound(self) -> None:
        norm = sgbl_interval_matrix_infinity_norm_2x2(
            interval(1, 2), interval(-3, -1), interval(0, 1), interval(4, 5)
        )
        self.assertEqual(norm, 6)

    def test_tube_around_a_singleton_has_the_frozen_radius(self) -> None:
        tube = sgbl_orbit_tube_around(ONE, DECLARED_TUBE_RADIUS)
        self.assertEqual(tube, interval(Q(7, 8), Q(9, 8)))
        self.assertTrue(ONE.strictly_inside(tube))


class SGBLOrbitLocalNominalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec = SGBLExactInitialSlice()
        cls.ode = sgbl_validated_constraint_ode(cls.spec)
        cls.record = sgbl_bind_orbit_local_picard_lindelof(
            cls.spec, prefix_ode=cls.ode
        )

    def test_authenticated_prefix_matches_the_frozen_identities(self) -> None:
        last = sgbl_authenticate_ck_prefix(self.ode, spec=self.spec)
        self.assertEqual(self.ode.chart_kind, CERTIFICATE_CHART_KIND)
        self.assertEqual(self.ode.obstruction, PREFIX_OBSTRUCTION)
        self.assertEqual(len(self.ode.cells), PREFIX_CELL_COUNT)
        self.assertEqual(self.ode.coverage_right, PREFIX_RADIUS)
        self.assertEqual(self.ode.support_compactness.upper, PREFIX_C_UPPER)
        self.assertEqual(last.radius.lower, PREFIX_LAST_CELL_LEFT)
        self.assertEqual(last.radius.upper, PREFIX_RADIUS)
        self.assertEqual(last.denominator_margin, PREFIX_D_MARGIN)
        self.assertGreater(last.lambda_right.width(), 2 * DECLARED_TUBE_RADIUS)

    def test_nominal_family_is_wrapping_not_a_compactness_or_jacobian_nonpass(self) -> None:
        self.assertEqual(self.record.classification, "interval_inconclusive")
        self.assertEqual(
            self.record.obstruction, "prefix_endpoint_not_inside_declared_tube"
        )
        self.assertEqual(
            self.record.missing_theorem,
            "authenticated_prefix_endpoint_contained_in_the_declared_orbit_tube",
        )
        self.assertFalse(self.record.unique_affine_orbit_on_remaining_support)
        self.assertFalse(self.record.tiles_remaining_support)
        self.assertEqual(len(self.record.tubes), 1)
        tube = self.record.tubes[0]
        self.assertFalse(tube.initial_inside_tube)
        self.assertFalse(tube.unique_local_affine_orbit)
        self.assertNotEqual(tube.classification, "on_orbit_compactness_not_below_one")
        self.assertNotEqual(tube.classification, "jacobian_diagonal_contains_zero")
        self.assertEqual(tube.depth, 0)
        self.assertEqual(self.record.slice.chi_amplitude, 3)
        self.assertEqual(self.record.contract_sha256, NOMINAL_CONTRACT_SHA256)
        self.assertEqual(
            sgbl_orbit_local_contract_sha256(dict(self.record.contract_payload)),
            NOMINAL_CONTRACT_SHA256,
        )

    def test_aggregate_flags_remain_false(self) -> None:
        gate = sgbl_orbit_local_health_gate(self.record)
        for name in FALSE_CLAIMS:
            self.assertIs(gate[name], False)
            self.assertIs(getattr(self.record, name, False) if hasattr(self.record, name) else False, False)
        self.assertIs(self.record.SGBL_branch_owned_and_healthy, False)
        self.assertIs(self.record.FRZ1, False)
        self.assertIs(self.record.PREF1, False)
        self.assertIs(self.record.holdout_authorized, False)
        self.assertIs(self.record.execution_authorized, False)
        self.assertIs(self.record.actual_orbit_traps, False)
        self.assertIs(self.record.sgbl_model_rejected, False)
        self.assertIs(self.record.whole_domain_nagumo_invariant, False)
        self.assertTrue(self.record.sampled_nodes_are_not_the_certificate)
        self.assertTrue(gate["sampled_nodes_are_not_the_certificate"])

    def test_inherited_affine_apis_remain_the_ivp(self) -> None:
        last = self.ode.cells[-1]
        coefficients = sgbl_interval_affine_coefficients(
            radius=last.radius,
            radial_metric=last.lambda_box,
            angular_extrinsic_curvature=last.k_box,
            spec=self.spec,
        )
        self.assertIn("hamiltonian_lambda_r", coefficients)
        jacobian = sgbl_constraint_rhs_state_jacobian(
            last.radius, last.lambda_box, last.k_box, self.spec
        )
        self.assertIn("d_lambda_r_d_lambda", jacobian)


class SGBLOrbitLocalSmallAmplitudeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.spec = SGBLExactInitialSlice(chi_amplitude=SMALL_AMPLITUDE_CONTROL)
        cls.ode = sgbl_validated_constraint_ode(cls.spec)
        cls.record = sgbl_orbit_local_small_amplitude_control(
            cls.spec, prefix_ode=cls.ode
        )

    def test_small_amplitude_is_a_control_and_never_replaces_the_nominal_family(self) -> None:
        self.assertEqual(self.spec.chi_amplitude, Q(1, 8))
        self.assertNotEqual(self.spec.chi_amplitude, NOMINAL_CHI_AMPLITUDE)
        self.assertEqual(self.record.slice.chi_amplitude, Q(1, 8))
        self.assertIs(
            self.record.contract_payload["replaced_nominal_family_by_small_amplitude"],
            False,
        )
        self.assertGreaterEqual(len(self.record.tubes), 1)
        self.assertTrue(self.record.tubes[0].unique_local_affine_orbit)
        self.assertEqual(self.record.tubes[0].classification, "unique_local_affine_orbit")
        self.assertLess(self.record.tubes[0].compactness.upper, 1)
        self.assertFalse(self.record.unique_affine_orbit_on_remaining_support)
        self.assertEqual(self.record.classification, "interval_inconclusive")
        self.assertEqual(self.record.obstruction, "picard_strict_self_map_failed")
        self.assertEqual(self.record.contract_sha256, SMALL_AMPLITUDE_CONTRACT_SHA256)
        gate = sgbl_orbit_local_health_gate(self.record)
        self.assertIs(gate["SGBL_branch_owned_and_healthy"], False)
        self.assertIs(gate["FRZ1"], False)


class SGBLOrbitLocalAttackTests(unittest.TestCase):
    def test_changed_policy_is_refused(self) -> None:
        changed = replace(DECLARED_ORBIT_TUBE_POLICY, tube_radius=Q(1, 16))
        with self.assertRaises(SGBLOrbitLocalStop) as stopped:
            sgbl_bind_orbit_local_picard_lindelof(policy=changed)
        self.assertEqual(stopped.exception.reason, "changed_policy")
        with self.assertRaises(SGBLOrbitLocalStop) as steps:
            sgbl_bind_orbit_local_picard_lindelof(
                policy=replace(DECLARED_ORBIT_TUBE_POLICY, continuation_steps=8)
            )
        self.assertEqual(steps.exception.reason, "changed_policy")
        with self.assertRaises(SGBLOrbitLocalStop) as depth:
            sgbl_bind_orbit_local_picard_lindelof(
                policy=replace(DECLARED_ORBIT_TUBE_POLICY, max_bisection_depth=1)
            )
        self.assertEqual(depth.exception.reason, "changed_policy")

    def test_wrong_predecessor_hash_is_refused(self) -> None:
        with self.assertRaises(SGBLOrbitLocalStop) as stopped:
            sgbl_bind_orbit_local_picard_lindelof(
                predecessor_trap_refinement_contract_sha256="0" * 64
            )
        self.assertEqual(stopped.exception.reason, "wrong_predecessor")
        with self.assertRaises(SGBLOrbitLocalStop) as barrier:
            sgbl_bind_orbit_local_picard_lindelof(
                predecessor_trap_barrier_contract_sha256="1" * 64
            )
        self.assertEqual(barrier.exception.reason, "wrong_predecessor")

    def test_family_mutation_is_refused(self) -> None:
        with self.assertRaises(SGBLOrbitLocalStop) as stopped:
            sgbl_bind_orbit_local_picard_lindelof(
                SGBLExactInitialSlice(chi_amplitude=Q(1, 4))
            )
        self.assertEqual(stopped.exception.reason, "family_mutation")
        with self.assertRaises(SGBLOrbitLocalStop) as small:
            sgbl_orbit_local_small_amplitude_control(
                SGBLExactInitialSlice(chi_amplitude=Q(1, 4))
            )
        self.assertEqual(small.exception.reason, "family_mutation")

    def test_resource_limit_is_not_a_compactness_nonpass(self) -> None:
        tiny = replace(DECLARED_ORBIT_TUBE_POLICY, max_rational_bit_length=1)
        tube = sgbl_orbit_tube_picard_lindelof(
            PREFIX_RADIUS,
            PREFIX_RADIUS + DECLARED_STEP_WIDTH,
            ONE,
            ZERO,
            SGBLExactInitialSlice(),
            policy=tiny,
        )
        self.assertEqual(tube.classification, "resource_limit")
        self.assertEqual(tube.obstruction, "max_rational_bit_length")
        self.assertFalse(tube.unique_local_affine_orbit)
        self.assertNotEqual(tube.classification, "on_orbit_compactness_not_below_one")

    def test_wrapping_endpoint_is_interval_inconclusive(self) -> None:
        wide = interval(Q(1, 2), Q(4))
        tube = sgbl_orbit_tube_picard_lindelof(
            PREFIX_RADIUS,
            PREFIX_RADIUS + DECLARED_STEP_WIDTH,
            wide,
            ZERO,
            SGBLExactInitialSlice(),
        )
        self.assertEqual(tube.classification, "interval_inconclusive")
        self.assertEqual(tube.obstruction, "prefix_endpoint_not_inside_declared_tube")
        self.assertFalse(tube.initial_inside_tube)
        self.assertNotEqual(tube.classification, "jacobian_diagonal_contains_zero")

    def test_dropped_or_gapped_inventory_is_refused(self) -> None:
        origin = PREFIX_RADIUS
        h = DECLARED_STEP_WIDTH
        tubes = (
            _dummy_tube(radius=interval(origin, origin + h)),
            _dummy_tube(radius=interval(origin + h, origin + 2 * h)),
            _dummy_tube(radius=interval(origin + 2 * h, origin + 3 * h)),
        )
        self.assertEqual(
            sgbl_validate_orbit_tube_inventory(tubes, origin=origin)[0], origin
        )
        with self.assertRaises(SGBLOrbitLocalStop) as dropped:
            sgbl_orbit_inventory_span(tubes[1:], origin=origin)
        self.assertEqual(dropped.exception.reason, "gapped_inventory")
        with self.assertRaises(SGBLOrbitLocalStop) as gapped:
            sgbl_orbit_inventory_span((tubes[0], tubes[2]), origin=origin)
        self.assertEqual(gapped.exception.reason, "gapped_inventory")

    def test_compactness_crossing_cannot_be_recorded_as_a_unique_pass(self) -> None:
        base = sgbl_picard_lindelof_constant_ode_control()
        with self.assertRaises(SGBLOrbitLocalStop) as stopped:
            replace(
                base,
                compactness=interval(0, Q(3, 2)),
                classification="unique_local_affine_orbit",
            )
        self.assertEqual(stopped.exception.reason, "tampered_contract")

    def test_claim_promotion_and_tampered_hash_are_refused(self) -> None:
        with self.assertRaises(SGBLOrbitLocalStop) as claims:
            sgbl_bind_orbit_local_picard_lindelof(claims={"FRZ1": True})
        self.assertEqual(claims.exception.reason, "claim_mutation")
        record = sgbl_picard_lindelof_constant_ode_control()
        with self.assertRaises(TypeError):
            replace(record, FRZ1=True)
        bound = sgbl_bind_orbit_local_picard_lindelof()
        with self.assertRaises(TypeError):
            replace(bound, SGBL_branch_owned_and_healthy=True)
        with self.assertRaises(SGBLOrbitLocalStop):
            replace(bound, contract_sha256="0" * 64)
        with self.assertRaises(SGBLOrbitLocalStop):
            SGBLOrbitLocalRecord(
                slice=bound.slice,
                policy=bound.policy,
                tubes=bound.tubes,
                classification=bound.classification,
                obstruction=bound.obstruction,
                missing_theorem=bound.missing_theorem,
                coverage_left=bound.coverage_left,
                coverage_right=bound.coverage_right,
                compactness_upper=bound.compactness_upper,
                compactness_margin=bound.compactness_margin,
                predecessor_trap_refinement_contract_sha256=(
                    bound.predecessor_trap_refinement_contract_sha256
                ),
                predecessor_trap_refinement2_contract_sha256=(
                    bound.predecessor_trap_refinement2_contract_sha256
                ),
                predecessor_trap_taylor_contract_sha256=(
                    bound.predecessor_trap_taylor_contract_sha256
                ),
                predecessor_trap_barrier_contract_sha256=(
                    bound.predecessor_trap_barrier_contract_sha256
                ),
                contract_payload=dict(bound.contract_payload)
                | {"SGBL_branch_owned_and_healthy": True},
                contract_sha256=bound.contract_sha256,
            )


class SGBLOrbitLocalCompatibilityTests(unittest.TestCase):
    def test_family_and_incomplete_feeders_stay_closed(self) -> None:
        parameters = SGBLFamilyParameters()
        spec = SGBLExactInitialSlice.from_family_parameters(parameters)
        self.assertEqual(spec.chi_amplitude, 3)
        gate = sgbl_branch_health_gate()
        self.assertIs(gate["SGBL_branch_owned_and_healthy"], False)
        health = sgbl_continuous_initial_compactness(spec)
        self.assertEqual(health.classification, "interval_inconclusive")
        principal = sgbl_family_principal_from_health(health, spec=spec)
        self.assertEqual(principal.classification, "family_geometry_incomplete")
        self.assertIs(sgbl_family_principal_health_gate(principal)["FRZ1"], False)
        admission = sgbl_whole_cell_admission_from_health(health, spec=spec)
        self.assertEqual(admission.classification, "family_geometry_incomplete")
        self.assertIs(sgbl_cell_admission_health_gate(admission)["PREF1"], False)


class SGBLOrbitLocalImportTests(unittest.TestCase):
    def test_module_does_not_import_qr_campaign_pro20_or_pass_bits(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_sol1 as owner

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
            "WeakCouplingThresholds",
            "weak_coupling_health_certificate",
            "FGCQRActionParameters",
            "solve_initial_data",
            "run_fgc_gr0_calibration_v1",
            "pro20_ev1_runtime",
            "pro20_ev1_protocol",
            "build_artifact_catalog",
            "artifact-catalog",
            "numpy",
            "sgbl_validated_lambda_k_constraint_ode",
            "sgbl_evaluate_product_box_ladder_level",
            "sgbl_evaluate_taylor_resource_ladder",
            "sgbl_bind_nagumo_barrier",
        }
        self.assertTrue(forbidden.isdisjoint(imported))
        self.assertTrue(set(FORBIDDEN_HEALTH_IMPORTS).isdisjoint(imported))
        source = Path(owner.__file__).read_text()
        self.assertNotIn("GB_EFFECTIVE_DENSITY_BOUND", source)
        self.assertNotIn("1e-12", source)
        self.assertIn("PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256", source)
        self.assertEqual(INSTRUMENT_ID, "FGC-1-SGB1-CTL1-SOL1")


if __name__ == "__main__":
    unittest.main()
