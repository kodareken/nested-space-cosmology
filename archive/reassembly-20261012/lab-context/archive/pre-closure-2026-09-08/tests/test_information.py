from __future__ import annotations

import math
import sys
import unittest
from decimal import Decimal
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.information import (  # noqa: E402
    TransitionCodeSpec,
    analytic_child_entropy,
    analytic_radiation_entropy,
    apply_isometry,
    binary_entropy,
    coherent_erasure_isometry,
    expected_child_marginal,
    expected_radiation_marginal,
    max_abs_difference,
    no_cloning_overlap_contrast,
    partial_trace_child,
    partial_trace_radiation,
    recover_code_state,
    state_density,
    trace,
    validate_transition,
    von_neumann_entropy,
)


class InformationMapTests(unittest.TestCase):
    def setUp(self) -> None:
        self.pure_qubit = state_density([1 / math.sqrt(2), 1j / math.sqrt(2)])
        self.mixed_qubit = [[0.75 + 0.0j, 0.0j], [0.0j, 0.25 + 0.0j]]

    def test_isometry_has_required_shape_and_preserves_inner_products(self) -> None:
        for spec in (
            TransitionCodeSpec(2, 0.0),
            TransitionCodeSpec(2, 0.37),
            TransitionCodeSpec(3, 1.0),
        ):
            isometry = coherent_erasure_isometry(spec)
            self.assertEqual((len(isometry), len(isometry[0])), (spec.output_dimension, spec.d))
            audit = validate_transition(spec, state_density([1.0] + [0.0] * (spec.d - 1)))
            self.assertTrue(audit["is_isometry"])
            self.assertLess(float(audit["isometry_error"]), 1e-14)

    def test_trace_and_hermiticity_hold_for_pure_and_mixed_inputs(self) -> None:
        for rho in (self.pure_qubit, self.mixed_qubit):
            audit = validate_transition(TransitionCodeSpec(2, 0.37), rho)
            self.assertTrue(audit["is_trace_preserving"])
            self.assertTrue(audit["output_is_hermitian"])
            self.assertTrue(audit["output_is_positive_semidefinite"])
            output_trace = trace(apply_isometry(TransitionCodeSpec(2, 0.37), rho))
            self.assertAlmostEqual(output_trace.real, 1.0, places=14)

    def test_joint_code_range_round_trip_recovers_pure_and_mixed_inputs(self) -> None:
        spec = TransitionCodeSpec(2, 0.37)
        for rho in (self.pure_qubit, self.mixed_qubit):
            recovered = recover_code_state(spec, apply_isometry(spec, rho))
            self.assertLess(max_abs_difference(recovered, rho), 1e-14)
            audit = validate_transition(spec, rho)
            self.assertTrue(audit["is_recoverable_on_code_range"])
            self.assertLess(float(audit["recovery_error"]), 1e-14)

    def test_reduced_states_are_exact_erasure_blocks(self) -> None:
        for probability in (0.0, 0.37, 1.0):
            spec = TransitionCodeSpec(2, probability)
            output = apply_isometry(spec, self.pure_qubit)
            self.assertLess(
                max_abs_difference(
                    partial_trace_child(spec, output),
                    expected_radiation_marginal(spec, self.pure_qubit),
                ),
                1e-14,
            )
            self.assertLess(
                max_abs_difference(
                    partial_trace_radiation(spec, output),
                    expected_child_marginal(spec, self.pure_qubit),
                ),
                1e-14,
            )

    def test_three_dimensional_code_accepts_a_mixed_state(self) -> None:
        spec = TransitionCodeSpec(3, 0.37)
        rho = [
            [0.5 + 0.0j, 0.0j, 0.0j],
            [0.0j, 0.3 + 0.0j, 0.0j],
            [0.0j, 0.0j, 0.2 + 0.0j],
        ]
        output = apply_isometry(spec, rho)
        self.assertAlmostEqual(trace(output).real, 1.0, places=14)
        self.assertLess(
            max_abs_difference(
                partial_trace_child(spec, output),
                expected_radiation_marginal(spec, rho),
            ),
            1e-14,
        )
        self.assertLess(
            max_abs_difference(
                partial_trace_radiation(spec, output),
                expected_child_marginal(spec, rho),
            ),
            1e-14,
        )

    def test_pure_and_mixed_entropy_formulas(self) -> None:
        spec = TransitionCodeSpec(2, 0.37)
        for rho in (self.pure_qubit, self.mixed_qubit):
            input_entropy = von_neumann_entropy(rho)
            output = apply_isometry(spec, rho)
            self.assertAlmostEqual(
                von_neumann_entropy(partial_trace_child(spec, output)),
                analytic_radiation_entropy(spec, input_entropy),
                places=11,
            )
            self.assertAlmostEqual(
                von_neumann_entropy(partial_trace_radiation(spec, output)),
                analytic_child_entropy(spec, input_entropy),
                places=11,
            )

    def test_endpoint_entropy_and_global_fine_grained_entropy(self) -> None:
        self.assertEqual(binary_entropy(0.0), 0.0)
        self.assertEqual(binary_entropy(1.0), 0.0)
        for probability in (0.0, 0.37, 1.0):
            output = apply_isometry(TransitionCodeSpec(2, probability), self.pure_qubit)
            self.assertAlmostEqual(von_neumann_entropy(output), 0.0, places=11)

    def test_no_cloning_contrast_and_no_double_deterministic_copy(self) -> None:
        contrast = no_cloning_overlap_contrast()
        self.assertNotEqual(contrast["input_overlap"], contrast["putative_cloned_output_overlap"])
        spec = TransitionCodeSpec(2, 0.37)
        output = apply_isometry(spec, self.pure_qubit)
        radiation = partial_trace_child(spec, output)
        child = partial_trace_radiation(spec, output)
        fidelity_r = sum(
            (
                self.pure_qubit[row][column] * radiation[column][row]
                for row in range(2)
                for column in range(2)
            ),
            0.0j,
        ).real
        fidelity_c = sum(
            (
                self.pure_qubit[row][column] * child[column][row]
                for row in range(2)
                for column in range(2)
            ),
            0.0j,
        ).real
        self.assertAlmostEqual(fidelity_r, 1.0 - spec.p, places=14)
        self.assertAlmostEqual(fidelity_c, spec.p, places=14)
        self.assertLess(fidelity_r, 1.0)
        self.assertLess(fidelity_c, 1.0)

    def test_invalid_arguments_are_rejected(self) -> None:
        for bad_d in (0, -1, 1.0, True):
            with self.subTest(bad_d=bad_d):
                with self.assertRaises(ValueError):
                    TransitionCodeSpec(bad_d, 0.5)  # type: ignore[arg-type]
        for bad_p in (
            -0.1,
            1.1,
            math.inf,
            math.nan,
            "0.5",
            None,
            True,
            0.5 + 0.0j,
            Decimal("0.5"),
            10**10_000,
        ):
            with self.subTest(bad_p=bad_p):
                with self.assertRaises(ValueError):
                    TransitionCodeSpec(2, bad_p)  # type: ignore[arg-type]
        self.assertEqual(TransitionCodeSpec(2, 1).p, 1.0)
        with self.assertRaises(ValueError):
            apply_isometry(TransitionCodeSpec(2, 0.5), [[1.0]])
        with self.assertRaises(ValueError):
            state_density([2.0])
        for bad_amplitudes in ([math.nan], [math.inf], [complex(0.0, math.nan)]):
            with self.subTest(bad_amplitudes=bad_amplitudes):
                with self.assertRaises(ValueError):
                    state_density(bad_amplitudes)
        for bad_entropy in (math.nan, math.inf, -math.inf):
            with self.subTest(bad_entropy=bad_entropy):
                with self.assertRaises(ValueError):
                    analytic_radiation_entropy(TransitionCodeSpec(2, 0.5), bad_entropy)
                with self.assertRaises(ValueError):
                    analytic_child_entropy(TransitionCodeSpec(2, 0.5), bad_entropy)
        with self.assertRaises(ValueError):
            validate_transition(TransitionCodeSpec(1, 0.5), [[complex(math.nan, 0.0)]])


if __name__ == "__main__":
    unittest.main()
