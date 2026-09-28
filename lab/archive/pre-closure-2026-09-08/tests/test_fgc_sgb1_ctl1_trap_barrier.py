from __future__ import annotations

import ast
from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.exact_interval import Interval, interval  # noqa: E402
from recursive_horizons.fgc.sgb1_ctl1_constraints import (  # noqa: E402
    SGBLAnnularInputs,
    sgbl_constraint_coefficients,
)
from recursive_horizons.fgc.sgb1_ctl1_initial_health import (  # noqa: E402
    SGBLExactInitialSlice,
    sgbl_interval_affine_coefficients,
    sgbl_misner_sharp_compactness_interval,
    sgbl_misner_sharp_compactness_radial_derivative,
)
from recursive_horizons.fgc.sgb1_ctl1_trap_barrier import (  # noqa: E402
    AFFINE_DENOMINATOR,
    BARRIER_FORMULA,
    BARRIER_K_SIGNS,
    CLASSIFICATION,
    DECLARED_SCALAR_MASS_MAJORANT,
    INSTRUMENT_ID,
    MISSING_THEOREM,
    NEIGHBORHOOD_D_R_STRICT_UPPER,
    NEIGHBORHOOD_RADIUS_OPEN_LEFT,
    NEIGHBORHOOD_RADIUS_OPEN_RIGHT,
    NOMINAL_CHI_AMPLITUDE,
    OUTWARD_K_SIGN,
    SMALL_AMPLITUDE_CONTROL,
    SUBSTITUTED_D_R_FORMULA,
    VACUUM_D_R_IDENTITY,
    WITNESS_D_R_STRICT_UPPER,
    WITNESS_K,
    WITNESS_LAMBDA,
    WITNESS_RADIUS,
    SGBLTrapBarrierRecord,
    SGBLTrapBarrierStop,
    sgbl_affine_substituted_d_r,
    sgbl_bind_nagumo_barrier,
    sgbl_exact_barrier_k,
    sgbl_exact_barrier_value,
    sgbl_misner_sharp_barrier_value,
    sgbl_null_expansion_barrier_identity,
    sgbl_parameterize_barrier_k,
    sgbl_polar_areal_null_expansions,
    sgbl_trap_barrier_aggregate_flags,
    sgbl_trap_barrier_contract_sha256,
    sgbl_vacuum_barrier_d_r_identity,
    sgbl_vacuum_substituted_d_r,
)


MODULE_PATH = ROOT / "src/recursive_horizons/fgc/sgb1_ctl1_trap_barrier.py"
NOMINAL_CONTRACT_SHA256 = (
    "a75d564c0ea0658802c138e0575edf386375ebba35728661bd8a11c579e7d689"
)
ZERO = interval(0)
ONE = interval(1)
TWO = interval(2)


def _independent_substituted_d_r(radius, lam, k, spec):
    radius = interval(radius) if type(radius) is not Interval else radius
    lam = interval(lam) if type(lam) is not Interval else lam
    k = interval(k) if type(k) is not Interval else k
    coefficients = sgbl_interval_affine_coefficients(
        radius=radius,
        radial_metric=lam,
        angular_extrinsic_curvature=k,
        spec=spec,
    )
    hamiltonian_lambda_r = coefficients["hamiltonian_lambda_r"]
    momentum_k_r = coefficients["momentum_k_r"]
    lambda_r = -coefficients["hamiltonian_constant"] / hamiltonian_lambda_r
    k_r = -coefficients["momentum_constant"] / momentum_k_r
    compactness_r = sgbl_misner_sharp_compactness_radial_derivative(
        radius, lam, k, lambda_r, k_r
    )
    substituted = (
        -TWO * radius * (k**2)
        + TWO * (radius**2) * k * (coefficients["momentum_constant"] / momentum_k_r)
        + TWO * coefficients["hamiltonian_constant"] / (hamiltonian_lambda_r * (lam**3))
    )
    return {
        "compactness": sgbl_misner_sharp_compactness_interval(radius, lam, k),
        "barrier": ONE - sgbl_misner_sharp_compactness_interval(radius, lam, k),
        "hamiltonian_lambda_r": hamiltonian_lambda_r,
        "momentum_k_r": momentum_k_r,
        "compactness_r": compactness_r,
        "d_r": substituted,
        "chain_rule": -compactness_r,
    }


class SGBLTrapBarrierFormulaTests(unittest.TestCase):
    def test_constants_are_frozen_before_the_bind(self) -> None:
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
            names.index("WITNESS_RADIUS"), names.index("sgbl_bind_nagumo_barrier")
        )
        self.assertLess(
            names.index("SUBSTITUTED_D_R_FORMULA"),
            names.index("sgbl_affine_substituted_d_r"),
        )
        self.assertEqual(INSTRUMENT_ID, "FGC-1-SGB1-CTL1-TRAP-BARRIER")
        self.assertEqual(CLASSIFICATION, "barrier_nonpass_not_trajectory_trapping")
        self.assertEqual(WITNESS_RADIUS, Q(209, 20))
        self.assertEqual(WITNESS_LAMBDA, 1)
        self.assertEqual(WITNESS_K, Q(20, 209))
        self.assertEqual(OUTWARD_K_SIGN, 1)
        self.assertEqual(BARRIER_K_SIGNS, (1, -1))
        self.assertEqual(WITNESS_D_R_STRICT_UPPER, -Q(3, 8))
        self.assertEqual(NEIGHBORHOOD_D_R_STRICT_UPPER, -Q(1, 20))
        self.assertEqual(NEIGHBORHOOD_RADIUS_OPEN_LEFT, Q(52, 5))
        self.assertEqual(NEIGHBORHOOD_RADIUS_OPEN_RIGHT, Q(21, 2))
        self.assertEqual(NOMINAL_CHI_AMPLITUDE, 3)
        self.assertEqual(SMALL_AMPLITUDE_CONTROL, Q(1, 8))
        self.assertEqual(DECLARED_SCALAR_MASS_MAJORANT, 6)
        self.assertEqual(AFFINE_DENOMINATOR, "H_L*M_k")
        self.assertIn("D_r=-2*r*k^2", SUBSTITUTED_D_R_FORMULA)
        self.assertEqual(BARRIER_FORMULA, "D=1-C=lambda^{-2}-r^2*k^2")
        self.assertEqual(VACUUM_D_R_IDENTITY, "D_r=1/r")
        self.assertIn("Nagumo", MISSING_THEOREM)

    def test_both_k_signs_parameterize_exact_D_zero(self) -> None:
        for sign in BARRIER_K_SIGNS:
            k = sgbl_exact_barrier_k(WITNESS_RADIUS, WITNESS_LAMBDA, sign)
            self.assertEqual(k, sign * WITNESS_K)
            self.assertEqual(
                sgbl_exact_barrier_value(WITNESS_RADIUS, WITNESS_LAMBDA, k),
                0,
            )
            self.assertEqual(
                sgbl_parameterize_barrier_k(WITNESS_RADIUS, WITNESS_LAMBDA, sign),
                Interval.singleton(k),
            )
            self.assertEqual(
                sgbl_misner_sharp_barrier_value(WITNESS_RADIUS, WITNESS_LAMBDA, k),
                ZERO,
            )

    def test_substituted_d_r_matches_independent_minus_chain_rule(self) -> None:
        spec = SGBLExactInitialSlice()
        owned = sgbl_affine_substituted_d_r(
            WITNESS_RADIUS, WITNESS_LAMBDA, WITNESS_K, spec
        )
        independent = _independent_substituted_d_r(
            WITNESS_RADIUS, WITNESS_LAMBDA, WITNESS_K, spec
        )
        self.assertEqual(owned["d_r"], independent["d_r"])
        self.assertEqual(owned["d_r"], independent["chain_rule"])
        self.assertEqual(owned["compactness"], ONE)
        self.assertEqual(owned["barrier"], ZERO)
        self.assertFalse(owned["hamiltonian_lambda_r"].contains_zero())
        self.assertFalse(owned["momentum_k_r"].contains_zero())
        self.assertTrue(owned["d_r"].upper < WITNESS_D_R_STRICT_UPPER)

    def test_vacuum_d_r_is_exactly_one_over_r_for_both_signs(self) -> None:
        for sign in BARRIER_K_SIGNS:
            vacuum = sgbl_vacuum_barrier_d_r_identity(
                WITNESS_RADIUS, WITNESS_LAMBDA, sign
            )
            self.assertEqual(vacuum["barrier"], 0)
            self.assertEqual(vacuum["d_r"], 1 / WITNESS_RADIUS)
            self.assertEqual(vacuum["d_r"], Q(20, 209))
            self.assertEqual(vacuum["identity"], vacuum["d_r"])
            k = sgbl_exact_barrier_k(WITNESS_RADIUS, WITNESS_LAMBDA, sign)
            self.assertEqual(
                sgbl_vacuum_substituted_d_r(WITNESS_RADIUS, WITNESS_LAMBDA, k),
                1 / WITNESS_RADIUS,
            )


class SGBLTrapBarrierWitnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = sgbl_bind_nagumo_barrier()
        cls.spec = SGBLExactInitialSlice()

    def test_outward_witness_is_exact_and_independently_reproduced(self) -> None:
        record = self.record
        self.assertEqual(record.witness_radius, Q(209, 20))
        self.assertEqual(record.witness_lambda, 1)
        self.assertEqual(record.witness_k, Q(20, 209))
        self.assertEqual(record.witness_k_sign, 1)
        self.assertEqual(record.barrier_value, ZERO)
        self.assertEqual(record.compactness, ONE)
        self.assertTrue(record.jacobian_diagonals_exclude_zero)
        self.assertTrue(record.d_r.upper < -Q(3, 8))
        self.assertTrue(record.outward_witness)
        independent = _independent_substituted_d_r(
            WITNESS_RADIUS, WITNESS_LAMBDA, WITNESS_K, self.spec
        )
        self.assertEqual(record.d_r, independent["d_r"])
        self.assertEqual(
            record.hamiltonian_lambda_r, independent["hamiltonian_lambda_r"]
        )
        self.assertEqual(record.momentum_k_r, independent["momentum_k_r"])
        self.assertEqual(record.classification, CLASSIFICATION)
        self.assertFalse(record.whole_domain_nagumo_invariant)
        self.assertFalse(record.actual_orbit_traps)
        self.assertFalse(record.trajectory_trapping_claimed)
        self.assertFalse(record.sgbl_model_rejected)

    def test_open_neighborhood_has_d_r_below_negative_one_twentieth(self) -> None:
        record = self.record
        self.assertEqual(record.neighborhood_radius.lower, Q(52, 5))
        self.assertEqual(record.neighborhood_radius.upper, Q(21, 2))
        self.assertTrue(
            record.neighborhood_radius.lower
            < WITNESS_RADIUS
            < record.neighborhood_radius.upper
        )
        self.assertTrue(record.neighborhood_d_r.upper < -Q(1, 20))
        self.assertTrue(record.neighborhood_outward)
        neighborhood_k = sgbl_parameterize_barrier_k(
            record.neighborhood_radius, WITNESS_LAMBDA, OUTWARD_K_SIGN
        )
        independent = _independent_substituted_d_r(
            record.neighborhood_radius,
            WITNESS_LAMBDA,
            neighborhood_k,
            self.spec,
        )
        self.assertEqual(record.neighborhood_d_r, independent["d_r"])
        self.assertTrue(independent["d_r"].upper < NEIGHBORHOOD_D_R_STRICT_UPPER)
        self.assertFalse(independent["hamiltonian_lambda_r"].contains_zero())
        self.assertFalse(independent["momentum_k_r"].contains_zero())

    def test_null_expansion_barrier_equivalence(self) -> None:
        plus, minus = sgbl_polar_areal_null_expansions(
            WITNESS_RADIUS, WITNESS_LAMBDA, WITNESS_K
        )
        self.assertEqual(plus, ZERO)
        self.assertEqual(minus, interval(-Q(80, 209)))
        self.assertEqual(self.record.theta_plus, plus)
        self.assertEqual(self.record.theta_minus, minus)
        identity = sgbl_null_expansion_barrier_identity(
            WITNESS_RADIUS, WITNESS_LAMBDA, WITNESS_K
        )
        self.assertEqual(identity["residual"], ZERO)
        self.assertEqual(identity["expansion_product"], identity["formula"])
        opposite = sgbl_polar_areal_null_expansions(
            WITNESS_RADIUS, WITNESS_LAMBDA, -WITNESS_K
        )
        self.assertEqual(opposite[1], ZERO)
        self.assertNotEqual(opposite[0], ZERO)
        self.assertTrue(self.record.null_expansion_barrier_equivalent)

    def test_opposite_sign_is_parameterized_and_not_the_outward_witness(self) -> None:
        self.assertEqual(self.record.opposite_sign_k, -WITNESS_K)
        self.assertEqual(
            sgbl_exact_barrier_value(
                WITNESS_RADIUS, WITNESS_LAMBDA, self.record.opposite_sign_k
            ),
            0,
        )
        self.assertTrue(self.record.opposite_sign_d_r.strictly_positive())
        self.assertFalse(self.record.opposite_sign_d_r.upper < WITNESS_D_R_STRICT_UPPER)


class SGBLTrapBarrierControlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = sgbl_bind_nagumo_barrier()

    def test_small_amplitude_is_an_inward_positive_control(self) -> None:
        small = sgbl_affine_substituted_d_r(
            WITNESS_RADIUS,
            WITNESS_LAMBDA,
            WITNESS_K,
            SGBLExactInitialSlice(chi_amplitude=SMALL_AMPLITUDE_CONTROL),
        )
        self.assertEqual(self.record.small_amplitude_d_r, small["d_r"])
        self.assertTrue(small["d_r"].strictly_positive())
        self.assertTrue(self.record.small_amplitude_inward)
        self.assertEqual(small["barrier"], ZERO)
        self.assertFalse(small["hamiltonian_lambda_r"].contains_zero())
        self.assertFalse(small["momentum_k_r"].contains_zero())
        self.assertIs(self.record.SGBL_branch_owned_and_healthy, False)
        self.assertIs(self.record.FRZ1, False)

    def test_scalar_mass_majorant_does_not_close_the_barrier(self) -> None:
        dropped = sgbl_affine_substituted_d_r(
            WITNESS_RADIUS,
            WITNESS_LAMBDA,
            WITNESS_K,
            SGBLExactInitialSlice(phi_amplitude=0),
        )
        majorant = sgbl_affine_substituted_d_r(
            WITNESS_RADIUS,
            WITNESS_LAMBDA,
            WITNESS_K,
            SGBLExactInitialSlice(scalar_mass=DECLARED_SCALAR_MASS_MAJORANT),
        )
        self.assertEqual(self.record.phi_dropped_d_r, dropped["d_r"])
        self.assertEqual(self.record.scalar_mass_majorant_d_r, majorant["d_r"])
        self.assertTrue(dropped["d_r"].upper < WITNESS_D_R_STRICT_UPPER)
        self.assertTrue(majorant["d_r"].upper < WITNESS_D_R_STRICT_UPPER)
        self.assertIs(self.record.scalar_mass_majorant_closes, False)
        self.assertGreater(DECLARED_SCALAR_MASS_MAJORANT, 3)

    def test_vacuum_control_is_exactly_the_witness_reciprocal(self) -> None:
        self.assertEqual(self.record.vacuum_d_r, Q(20, 209))
        self.assertEqual(self.record.vacuum_d_r, 1 / WITNESS_RADIUS)
        coefficients = sgbl_constraint_coefficients(
            SGBLAnnularInputs(
                radius=WITNESS_RADIUS,
                radial_metric=WITNESS_LAMBDA,
                angular_extrinsic_curvature=WITNESS_K,
                phi=0,
                phi_r=0,
                phi_rr=0,
                phi_pi=0,
                phi_pi_r=0,
                chi_r=0,
                chi_pi=0,
                planck_mass=2,
                scalar_mass=3,
                quartic_coupling=Q(1, 2),
                alpha_gb=-Q(1, 4),
            )
        )
        independent = (
            -2 * WITNESS_RADIUS * WITNESS_K**2
            + 2
            * WITNESS_RADIUS**2
            * WITNESS_K
            * (coefficients.momentum_constant / coefficients.momentum_k_r)
            + 2
            * coefficients.hamiltonian_constant
            / (coefficients.hamiltonian_lambda_r * WITNESS_LAMBDA**3)
        )
        self.assertEqual(independent, 1 / WITNESS_RADIUS)
        self.assertNotEqual(coefficients.hamiltonian_lambda_r, 0)
        self.assertNotEqual(coefficients.momentum_k_r, 0)


class SGBLTrapBarrierContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = sgbl_bind_nagumo_barrier()

    def test_aggregate_flags_and_nonpromoting_payload_stay_false(self) -> None:
        record = self.record
        gate = sgbl_trap_barrier_aggregate_flags(record)
        for name in (
            "SGBL_branch_owned_and_healthy",
            "FRZ1",
            "PREF1",
            "holdout_authorized",
            "execution_authorized",
            "whole_domain_nagumo_invariant",
            "trajectory_trapping_claimed",
            "actual_orbit_traps",
            "sgbl_model_rejected",
            "copied_gr0_or_fgcqr_health_evidence",
        ):
            with self.subTest(name=name):
                self.assertIs(getattr(record, name), False)
                self.assertIs(gate[name], False)
                self.assertIs(record.contract_payload[name], False)
        self.assertEqual(record.contract_payload["INSTRUMENT_ID"], INSTRUMENT_ID)
        self.assertEqual(record.contract_payload["classification"], CLASSIFICATION)
        self.assertEqual(
            record.contract_sha256,
            sgbl_trap_barrier_contract_sha256(record.contract_payload),
        )
        self.assertEqual(record.contract_sha256, NOMINAL_CONTRACT_SHA256)
        self.assertNotIn("wall", str(dict(record.contract_payload)))
        self.assertNotIn("diagnostic_wall_seconds", dict(record.contract_payload))
        self.assertIs(record.sampled_nodes_are_not_the_certificate, True)

    def test_bind_is_deterministic(self) -> None:
        again = sgbl_bind_nagumo_barrier()
        self.assertEqual(again.contract_sha256, self.record.contract_sha256)
        self.assertEqual(
            dict(again.contract_payload), dict(self.record.contract_payload)
        )
        self.assertEqual(again.d_r, self.record.d_r)


class SGBLTrapBarrierAttackTests(unittest.TestCase):
    def test_formula_mutation_is_refused(self) -> None:
        with self.assertRaises(SGBLTrapBarrierStop) as stopped:
            sgbl_bind_nagumo_barrier(formula="C_r")
        self.assertEqual(stopped.exception.reason, "formula_mutation")
        with self.assertRaises(SGBLTrapBarrierStop) as chain:
            sgbl_affine_substituted_d_r(
                WITNESS_RADIUS,
                WITNESS_LAMBDA,
                WITNESS_K,
                SGBLExactInitialSlice(),
                formula="D_r=-C_r",
            )
        self.assertEqual(chain.exception.reason, "formula_mutation")

    def test_sign_mutation_is_refused(self) -> None:
        with self.assertRaises(SGBLTrapBarrierStop) as stopped:
            sgbl_bind_nagumo_barrier(k_sign=-1)
        self.assertEqual(stopped.exception.reason, "sign_mutation")
        with self.assertRaises(SGBLTrapBarrierStop) as flipped:
            sgbl_bind_nagumo_barrier(angular_extrinsic_curvature=-WITNESS_K)
        self.assertEqual(flipped.exception.reason, "sign_mutation")
        with self.assertRaises(SGBLTrapBarrierStop) as boolean:
            sgbl_bind_nagumo_barrier(k_sign=True)
        self.assertEqual(boolean.exception.reason, "sign_mutation")

    def test_denominator_mutation_is_refused(self) -> None:
        with self.assertRaises(SGBLTrapBarrierStop) as stopped:
            sgbl_bind_nagumo_barrier(denominator="H_L")
        self.assertEqual(stopped.exception.reason, "denominator_mutation")
        with self.assertRaises(SGBLTrapBarrierStop) as dropped:
            sgbl_affine_substituted_d_r(
                WITNESS_RADIUS,
                WITNESS_LAMBDA,
                WITNESS_K,
                SGBLExactInitialSlice(),
                denominator="1",
            )
        self.assertEqual(dropped.exception.reason, "denominator_mutation")

    def test_witness_mutation_is_refused(self) -> None:
        with self.assertRaises(SGBLTrapBarrierStop) as radius:
            sgbl_bind_nagumo_barrier(radius=Q(10))
        self.assertEqual(radius.exception.reason, "witness_mutation")
        with self.assertRaises(SGBLTrapBarrierStop) as lam:
            sgbl_bind_nagumo_barrier(radial_metric=Q(2))
        self.assertEqual(lam.exception.reason, "witness_mutation")

    def test_family_mutation_is_refused(self) -> None:
        with self.assertRaises(SGBLTrapBarrierStop) as small:
            sgbl_bind_nagumo_barrier(
                SGBLExactInitialSlice(chi_amplitude=SMALL_AMPLITUDE_CONTROL)
            )
        self.assertEqual(small.exception.reason, "family_mutation")
        with self.assertRaises(SGBLTrapBarrierStop) as undeclared:
            sgbl_bind_nagumo_barrier(SGBLExactInitialSlice(chi_amplitude=Q(2)))
        self.assertEqual(undeclared.exception.reason, "family_mutation")

    def test_claim_mutation_is_refused(self) -> None:
        record = sgbl_bind_nagumo_barrier()
        for claim in (
            "FRZ1",
            "PREF1",
            "SGBL_branch_owned_and_healthy",
            "holdout_authorized",
            "actual_orbit_traps",
            "sgbl_model_rejected",
            "trajectory_trapping_claimed",
        ):
            with self.subTest(claim=claim):
                with self.assertRaises(SGBLTrapBarrierStop) as stopped:
                    sgbl_bind_nagumo_barrier(claims={claim: True})
                self.assertEqual(stopped.exception.reason, "claim_mutation")
        with self.assertRaises(TypeError):
            replace(record, FRZ1=True)
        with self.assertRaises(SGBLTrapBarrierStop) as classified:
            replace(record, classification="trajectory_trapping")
        self.assertEqual(classified.exception.reason, "claim_mutation")
        with self.assertRaises(SGBLTrapBarrierStop) as majorant:
            replace(record, scalar_mass_majorant_closes=True)
        self.assertEqual(majorant.exception.reason, "claim_mutation")
        with self.assertRaises(SGBLTrapBarrierStop) as digest:
            replace(record, contract_sha256="0" * 64)
        self.assertEqual(digest.exception.reason, "tampered_contract")
        with self.assertRaises(SGBLTrapBarrierStop) as payload:
            SGBLTrapBarrierRecord(
                slice=record.slice,
                classification=record.classification,
                witness_radius=record.witness_radius,
                witness_lambda=record.witness_lambda,
                witness_k=record.witness_k,
                witness_k_sign=record.witness_k_sign,
                compactness=record.compactness,
                barrier_value=record.barrier_value,
                hamiltonian_lambda_r=record.hamiltonian_lambda_r,
                momentum_k_r=record.momentum_k_r,
                d_r=record.d_r,
                opposite_sign_k=record.opposite_sign_k,
                opposite_sign_d_r=record.opposite_sign_d_r,
                neighborhood_radius=record.neighborhood_radius,
                neighborhood_d_r=record.neighborhood_d_r,
                vacuum_d_r=record.vacuum_d_r,
                small_amplitude_d_r=record.small_amplitude_d_r,
                scalar_mass_majorant_d_r=record.scalar_mass_majorant_d_r,
                phi_dropped_d_r=record.phi_dropped_d_r,
                theta_plus=record.theta_plus,
                theta_minus=record.theta_minus,
                scalar_mass_majorant_closes=False,
                contract_payload=dict(record.contract_payload) | {"FRZ1": True},
                contract_sha256=record.contract_sha256,
            )
        self.assertEqual(payload.exception.reason, "tampered_contract")


class SGBLTrapBarrierImportTests(unittest.TestCase):
    def test_module_does_not_import_health_ode_or_qr(self) -> None:
        from recursive_horizons.fgc import sgb1_ctl1_trap_barrier as owner

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
            "integrate_sgbl_radial_constraints",
            "build_artifact_catalog",
            "artifact-catalog",
            "numpy",
        }
        self.assertTrue(forbidden.isdisjoint(imported))
        source = Path(owner.__file__).read_text()
        self.assertNotIn("GB_EFFECTIVE_DENSITY_BOUND", source)
        self.assertNotIn("1e-12", source)
        self.assertNotIn("FRZ1 = True", source)
        self.assertNotIn("PREF1 = True", source)
        self.assertIn("barrier_nonpass_not_trajectory_trapping", source)
        self.assertEqual(INSTRUMENT_ID, "FGC-1-SGB1-CTL1-TRAP-BARRIER")


if __name__ == "__main__":
    unittest.main()
