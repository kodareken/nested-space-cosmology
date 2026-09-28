from __future__ import annotations

import ast
from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.modified_harmonic import (  # noqa: E402
    _physical_metric,
    auxiliary_inverse_metric,
    inverse_metric_null_polynomial,
)
from recursive_horizons.fgc.sgb1_ctl1_continuity import (  # noqa: E402
    CONDITIONAL_LOGICAL_CHAIN,
    CONDITIONAL_PROPAGATION_ASSUMPTIONS,
    CONDITIONAL_STATEMENT,
    FORBIDDEN_NEIGHBOR_IMPORTS,
    HOLDOUT_ACCELERATION,
    NAMED_HAT_WAVE_UNIQUENESS_PREMISE,
    SGBLConstraintPropagationContract,
    SGBLContinuityStop,
    UNCONDITIONAL_NONCLAIMS,
    sgbl_constraint_continuity,
    sgbl_constraint_propagation_contract,
    sgbl_fo1_reduction_identity,
    sgbl_hat_cone_characteristics,
    sgbl_mhg_hat_wave_identity,
    sgbl_physical_constraints_acceleration_independent,
    sgbl_physical_to_normal_gauge_map,
    sgbl_scalar_equations_unmodified,
    sgbl_source_jacobian_identity,
)
from recursive_horizons.fgc.sgb1_ctl1_source import (  # noqa: E402
    HAT_NORMAL_FACTOR,
    SGBLSourceInputs,
    sgbl_source_coefficients,
    sgbl_source_residual,
    sgbl_source_state,
)


def _fixture_a() -> SGBLSourceInputs:
    return SGBLSourceInputs(
        coordinate_radius=Q(5, 2),
        alpha=(2, Q(1, 3), Q(-1, 2), Q(1, 5), Q(-1, 7)),
        shift=(Q(1, 3), Q(-2, 5), Q(1, 4), Q(-1, 6), Q(2, 9)),
        radial_metric=(Q(3, 2), Q(1, 8), Q(-1, 3), Q(1, 7), Q(-1, 4)),
        areal_radius=(4, Q(-1, 5), Q(3, 2), Q(1, 9), Q(-2, 5)),
        phi=(Q(1, 2), Q(2, 3), Q(-4, 5), Q(1, 4), Q(-1, 8)),
        chi=(Q(-3, 4), Q(1, 2), Q(3, 5), Q(-2, 7), Q(1, 6)),
        planck_mass=2,
        scalar_mass=3,
        quartic_coupling=Q(1, 2),
        alpha_gb=Q(-1, 4),
    )


def _flat() -> SGBLSourceInputs:
    return SGBLSourceInputs(
        coordinate_radius=Q(5, 2),
        alpha=(1, 0, 0, 0, 0),
        shift=(0, 0, 0, 0, 0),
        radial_metric=(1, 0, 0, 0, 0),
        areal_radius=(Q(5, 2), 0, 1, 0, 0),
        phi=(0, 0, 0, 0, 0),
        chi=(0, 0, 0, 0, 0),
        planck_mass=2,
        scalar_mass=3,
        quartic_coupling=Q(1, 2),
        alpha_gb=Q(-1, 4),
    )


class SGBLContinuityIdentityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.point = _fixture_a()
        cls.record = sgbl_constraint_continuity(cls.point)

    def test_scalar_mhg_and_source_identities_hold_without_opening_health(self) -> None:
        record = self.record
        self.assertEqual(record.classification, "constraint_continuity_compatible")
        self.assertTrue(record.constraint_continuity_compatible)
        self.assertTrue(record.scalar_equations_unmodified)
        self.assertTrue(record.source_affine_identity)
        self.assertTrue(record.physical_constraints_acceleration_independent)
        self.assertTrue(record.hat_wave_identity)
        self.assertTrue(record.mhg_extension_metric_only)
        self.assertEqual(record.action["beta"], 0)
        self.assertEqual(record.action["eta"], 0)
        self.assertEqual(record.scalar_payload["F_prime"], 0)
        self.assertEqual(record.scalar_payload["F"], 4)
        self.assertEqual(record.hat_wave_payload["F"], 4)
        self.assertTrue(record.hat_wave_payload["independent_of_con4_certificate"])
        self.assertIs(record.constraint_propagation_on_a_domain_proven, False)
        self.assertIs(record.SGBL_branch_owned_and_healthy, False)
        self.assertIs(record.strongly_hyperbolic, False)
        self.assertEqual(record.cone_certificate, "unqualified")
        with self.assertRaises(TypeError):
            replace(record, SGBL_branch_owned_and_healthy=True)
        with self.assertRaises(TypeError):
            replace(record, constraint_propagation_on_a_domain_proven=True)

    def test_source_holdout_matches_an_independent_residual(self) -> None:
        payload = sgbl_source_jacobian_identity(self.point)
        coefficients = sgbl_source_coefficients(self.point)
        actual = sgbl_source_residual(self.point, HOLDOUT_ACCELERATION)
        self.assertEqual(coefficients.evaluate(HOLDOUT_ACCELERATION), actual)
        self.assertTrue(payload["affine_identity_holds"])
        flipped = tuple(-entry for entry in actual)
        self.assertNotEqual(flipped, actual)

    def test_physical_constraints_do_not_see_accelerations(self) -> None:
        payload = sgbl_physical_constraints_acceleration_independent(self.point)
        self.assertTrue(payload["physical_constraints_independent_of_accelerations"])
        self.assertTrue(payload["source_is_unredefined_residual"])

    def test_flat_vacuum_scalars_remain_unmodified(self) -> None:
        payload = sgbl_scalar_equations_unmodified(_flat())
        self.assertTrue(payload["scalar_equations_unmodified"])
        self.assertEqual(payload["phi_residual"], 0)
        self.assertEqual(payload["chi_residual"], 0)
        self.assertEqual(payload["extension_scalar_rows"], (Q(0), Q(0)))


class SGBLContinuityFailureTests(unittest.TestCase):
    def test_trace_sign_mutation_is_broken_gauge(self) -> None:
        with self.assertRaises(SGBLContinuityStop) as stopped:
            sgbl_mhg_hat_wave_identity(_fixture_a(), mutate_trace_sign=1)
        self.assertEqual(stopped.exception.reason, "broken_gauge")

    def test_zero_coupling_is_refused_as_a_production_continuity_point(self) -> None:
        with self.assertRaises(ValueError):
            sgbl_constraint_continuity(
                replace(_fixture_a(), alpha_gb=0)  # type: ignore[arg-type]
            )

    def test_module_does_not_import_neighboring_certificates(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_continuity as owner

        tree = ast.parse(Path(owner.__file__).read_text())
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.Import):
                imported.update(alias.name.split(".")[0] for alias in node.names)
        self.assertTrue(set(FORBIDDEN_NEIGHBOR_IMPORTS).isdisjoint(imported))
        self.assertNotIn("weak_coupling_health_certificate", imported)
        self.assertNotIn("activated_compatible_state", imported)
        self.assertNotIn("constraint_system_certificate", imported)
        self.assertNotIn("physical_to_normal_gauge_map", imported)
        self.assertNotIn("gauge_constraint_characteristics", imported)
        self.assertNotIn("modified_harmonic_constraint_propagation_certificate", imported)
        self.assertNotIn("constraint_system_certificate", imported)


class SGBLConstraintPropagationContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.point = _fixture_a()
        cls.contract = sgbl_constraint_propagation_contract(cls.point)

    def test_conditional_theorem_is_proved_without_existence_or_health(self) -> None:
        contract = self.contract
        self.assertEqual(
            contract.classification,
            "conditional_boundary_free_constraint_propagation_proved",
        )
        self.assertIs(
            contract.conditional_boundary_free_constraint_propagation_proved, True
        )
        self.assertIsInstance(contract, SGBLConstraintPropagationContract)
        self.assertTrue(contract.compatibility.constraint_continuity_compatible)
        self.assertEqual(contract.map_determinant, Q(-2916))
        self.assertTrue(contract.map_determinant_strictly_negative)
        self.assertLess(contract.map_determinant, 0)
        self.assertEqual(set(contract.hat_cone_root_signs), {-1, 1})
        self.assertEqual(contract.hat_cone_root_signs[0], -1)
        self.assertEqual(contract.hat_cone_root_signs[1], 1)
        self.assertEqual(contract.hat_cone["exact_roots"], (Q(-7, 9), Q(1, 9)))
        self.assertEqual(
            contract.hat_cone["hat_null_polynomial"],
            (Q(7, 36), Q(-3, 2), Q(-9, 4)),
        )
        self.assertTrue(contract.fo1_reduction["off_shell_identity_exact"])
        self.assertTrue(contract.fo1_reduction["field_by_field_zero"])
        self.assertEqual(
            contract.fo1_reduction["subsidiary_characteristic_speeds"],
            (Q(0),) * 6,
        )
        self.assertTrue(contract.fo1_reduction["off_shell_constraints_need_not_vanish"])
        self.assertNotEqual(contract.fo1_reduction["off_shell_D"], (Q(0),) * 6)
        self.assertEqual(contract.physical_to_normal_map["F"], 4)
        self.assertEqual(contract.physical_to_normal_map["F_prime"], 0)
        self.assertEqual(contract.physical_to_normal_map["hat_normal_factor"], 9)
        self.assertTrue(
            contract.physical_to_normal_map["direct_and_analytic_matrices_equal"]
        )
        self.assertTrue(contract.hat_cone["has_both_radial_branches"])
        self.assertTrue(contract.hat_cone["discriminant_positive"])
        self.assertEqual(contract.logical_chain, CONDITIONAL_LOGICAL_CHAIN)
        self.assertEqual(
            contract.named_uniqueness_premise, NAMED_HAT_WAVE_UNIQUENESS_PREMISE
        )
        self.assertEqual(contract.conditional_statement, CONDITIONAL_STATEMENT)
        self.assertIn("named linear-PDE premise", contract.named_uniqueness_premise)
        self.assertIn(
            "not an imported CON3 certificate", contract.named_uniqueness_premise
        )
        for name in CONDITIONAL_PROPAGATION_ASSUMPTIONS:
            self.assertIs(contract.assumptions[name], True)
        for name in UNCONDITIONAL_NONCLAIMS:
            self.assertIs(getattr(contract, name), False)
        self.assertIs(contract.constraint_propagation_on_a_domain_proven, False)
        self.assertIs(contract.FRZ1, False)
        self.assertIs(contract.PREF1, False)
        self.assertIs(contract.holdout, False)
        self.assertFalse(any(contract.remaining_nonclaims.values()))
        self.assertIn("boundary-free domain of dependence", contract.conditional_statement)

    def test_map_determinant_sign_root_and_inverse_are_exact(self) -> None:
        payload = sgbl_physical_to_normal_gauge_map(self.point)
        state = sgbl_source_state(self.point)
        effective_planck = self.point.planck_mass ** 2
        h_tt = state.h_tt.value
        h_tr = state.h_tr.value
        h_rr = state.h_rr.value
        lapse_squared = -h_tt + h_tr**2 / h_rr
        independent = (
            (effective_planck * HAT_NORMAL_FACTOR * lapse_squared / 2, Q(0)),
            (
                -effective_planck * HAT_NORMAL_FACTOR * h_tr / 2,
                -effective_planck * HAT_NORMAL_FACTOR * h_rr / 2,
            ),
        )
        independent_det = (
            -(effective_planck**2)
            * HAT_NORMAL_FACTOR**2
            * lapse_squared
            * h_rr
            / 4
        )
        self.assertEqual(payload["computed_matrix"], independent)
        self.assertEqual(payload["analytic_matrix"], independent)
        self.assertEqual(independent, ((Q(72), Q(0)), (Q(-27, 2), Q(-81, 2))))
        self.assertEqual(payload["determinant"], Q(-2916))
        self.assertEqual(independent_det, Q(-2916))
        self.assertLess(payload["determinant"], 0)
        inverse = payload["physical_to_normal_derivative_matrix"]
        self.assertEqual(inverse, payload["map_from_physical_HM_to_normal_nabla_C"])
        self.assertEqual(inverse, ((Q(1, 72), Q(0)), (Q(-1, 216), Q(-2, 81))))
        identity = (
            (
                payload["computed_matrix"][0][0] * inverse[0][0]
                + payload["computed_matrix"][0][1] * inverse[1][0],
                payload["computed_matrix"][0][0] * inverse[0][1]
                + payload["computed_matrix"][0][1] * inverse[1][1],
            ),
            (
                payload["computed_matrix"][1][0] * inverse[0][0]
                + payload["computed_matrix"][1][1] * inverse[1][0],
                payload["computed_matrix"][1][0] * inverse[0][1]
                + payload["computed_matrix"][1][1] * inverse[1][1],
            ),
        )
        self.assertEqual(identity, ((Q(1), Q(0)), (Q(0), Q(1))))
        self.assertEqual(
            self.contract.physical_to_normal_derivative_matrix, inverse
        )

        _, physical_inverse = _physical_metric(state)
        hat = auxiliary_inverse_metric(physical_inverse, HAT_NORMAL_FACTOR)
        polynomial = inverse_metric_null_polynomial(hat)
        characteristics = sgbl_hat_cone_characteristics(self.point)
        self.assertEqual(characteristics["hat_null_polynomial"], polynomial)
        self.assertEqual(polynomial, (Q(7, 36), Q(-3, 2), Q(-9, 4)))
        self.assertEqual(characteristics["exact_roots"], (Q(-7, 9), Q(1, 9)))
        self.assertEqual(characteristics["root_signs"], (-1, 1))
        self.assertLess(characteristics["exact_roots"][0], 0)
        self.assertGreater(characteristics["exact_roots"][1], 0)

    def test_broken_map_sign_hat_root_action_assumption_and_replace_attacks(
        self,
    ) -> None:
        with self.assertRaises(SGBLContinuityStop) as broken_map:
            sgbl_physical_to_normal_gauge_map(
                self.point, mutate_computed_scale=0
            )
        self.assertEqual(broken_map.exception.reason, "broken_map")
        with self.assertRaises(SGBLContinuityStop) as flipped_map:
            sgbl_physical_to_normal_gauge_map(
                self.point, mutate_computed_scale=-1
            )
        self.assertEqual(flipped_map.exception.reason, "broken_map")
        with self.assertRaises(SGBLContinuityStop) as sign:
            sgbl_physical_to_normal_gauge_map(
                self.point, mutate_analytic_sign=-1
            )
        self.assertEqual(sign.exception.reason, "broken_sign")
        with self.assertRaises(SGBLContinuityStop) as hat_locked:
            sgbl_physical_to_normal_gauge_map(self.point, mutate_hat_factor=1)
        self.assertEqual(hat_locked.exception.reason, "broken_hat")
        with self.assertRaises(SGBLContinuityStop) as hat_unlocked:
            sgbl_physical_to_normal_gauge_map(self.point, mutate_hat_factor=4)
        self.assertEqual(hat_unlocked.exception.reason, "broken_hat")
        with self.assertRaises(SGBLContinuityStop) as hat_char:
            sgbl_hat_cone_characteristics(self.point, mutate_hat_factor=1)
        self.assertEqual(hat_char.exception.reason, "broken_hat")
        with self.assertRaises(SGBLContinuityStop) as hat_poly:
            sgbl_hat_cone_characteristics(self.point, mutate_hat_factor=4)
        self.assertEqual(hat_poly.exception.reason, "broken_hat")
        with self.assertRaises(SGBLContinuityStop) as root:
            sgbl_hat_cone_characteristics(self.point, mutate_root_signs=True)
        self.assertEqual(root.exception.reason, "broken_root")
        with self.assertRaises(SGBLContinuityStop) as fo1:
            sgbl_fo1_reduction_identity(self.point, mutate_identity=True)
        self.assertEqual(fo1.exception.reason, "broken_assumption")
        with self.assertRaises(SGBLContinuityStop) as action:
            sgbl_physical_to_normal_gauge_map(replace(_fixture_a(), alpha_gb=0))
        self.assertEqual(action.exception.reason, "action_identity_error")

        contract = self.contract
        with self.assertRaises(TypeError):
            replace(contract, SGBL_branch_owned_and_healthy=True)
        with self.assertRaises(TypeError):
            replace(contract, IBVP_proved=True)
        with self.assertRaises(TypeError):
            replace(contract, smooth_solution_exists=True)
        with self.assertRaises(TypeError):
            replace(contract, FRZ1=True)
        with self.assertRaises(ValueError):
            replace(
                contract,
                conditional_boundary_free_constraint_propagation_proved=False,
            )
        broken_assumptions = {
            name: True for name in CONDITIONAL_PROPAGATION_ASSUMPTIONS
        }
        broken_assumptions["boundary_free_domain_of_dependence"] = False
        with self.assertRaises(ValueError):
            replace(contract, assumptions=broken_assumptions)
        with self.assertRaises(ValueError):
            replace(contract, named_uniqueness_premise="imported CON3 certificate")
        with self.assertRaises(ValueError):
            replace(contract, logical_chain=("not the SGB-L chain",))


if __name__ == "__main__":
    unittest.main()
