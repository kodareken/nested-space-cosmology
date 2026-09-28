from __future__ import annotations

from pathlib import Path
import json
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionCheckpoint,
    EvolutionRHS,
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
    accept_step,
    propose_step,
)


class FGCNumericalEngineTests(unittest.TestCase):
    def test_frozen_sbp_operators_satisfy_identity_and_center_parity(self) -> None:
        for order in (4, 2):
            with self.subTest(order=order):
                grid = UniformRadialGrid(0.0, 1.0, 65)
                operator = SBPFirstDerivative(grid, order)
                self.assertLessEqual(operator.sbp_identity_residual(), 1.0e-14)
                radius = grid.coordinates
                values = np.column_stack((1.0 + radius**4, radius + radius**3))
                derivative = operator.differentiate(
                    values, center_parities=(1, -1)
                )
                self.assertEqual(derivative[0, 0], 0.0)
                tolerance = 1.0e-12 if order == 4 else 3.0e-4
                self.assertLessEqual(abs(derivative[0, 1] - 1.0), tolerance)

                alternating = ((-1.0) ** np.arange(grid.point_count))[:, None]
                dissipation = operator.kreiss_oliger_dissipation(
                    alternating, coefficient=1.0 / 64.0
                )
                energy_rate = np.sum(
                    operator.norm_weights[:, None] * alternating * dissipation
                )
                self.assertLess(energy_rate, 0.0)

    def test_both_methods_evaluate_an_independent_candidate_endpoint(self) -> None:
        initial = EvolutionState(
            np.zeros((9, 1)), np.ones((9, 1)), np.zeros((9, 1))
        )

        def rhs(time: float, state: EvolutionState) -> EvolutionRHS:
            shape = state.shape
            return EvolutionRHS(
                state.p,
                np.full(shape, time),
                np.zeros(shape),
                {"evaluation_time": time},
            )

        expected = {
            PRIMARY_METHOD: (5, "rk4_k1", "candidate_endpoint"),
            COMPARATOR_METHOD: (4, "ssprk3_s0", "candidate_endpoint"),
        }
        for method, (count, first, last) in expected.items():
            with self.subTest(method=method):
                proposal = propose_step(
                    method=method,
                    time=0.0,
                    step_size=0.125,
                    state=initial,
                    rhs=rhs,
                )
                self.assertEqual(len(proposal.stages), count)
                self.assertEqual(proposal.stages[0].stage_name, first)
                self.assertEqual(proposal.stages[-1].stage_name, last)
                self.assertIsNot(
                    proposal.stages[-1].state,
                    proposal.stages[-2].state,
                )
                accepted = accept_step(
                    proposal,
                    previous_step_index=0,
                    previous_transaction_serial=0,
                )
                self.assertEqual(accepted.transaction_serial, count)
                self.assertEqual(accepted.time, 0.125)

    def test_guard_failure_exposes_no_accepted_endpoint(self) -> None:
        initial = EvolutionState(
            np.zeros((9, 1)), np.zeros((9, 1)), np.zeros((9, 1))
        )

        def rhs(_time: float, state: EvolutionState) -> EvolutionRHS:
            zero = np.zeros(state.shape)
            return EvolutionRHS(zero, zero, zero, {})

        proposal = propose_step(
            method=PRIMARY_METHOD,
            time=0.0,
            step_size=0.1,
            state=initial,
            rhs=rhs,
        )

        def reject(_proposal):
            raise RuntimeError("injected stop")

        with self.assertRaisesRegex(RuntimeError, "injected stop"):
            accept_step(
                proposal,
                previous_step_index=0,
                previous_transaction_serial=0,
                guards=(reject,),
            )

    def test_checkpoint_is_bitwise_canonical_and_tamper_evident(self) -> None:
        state = EvolutionState(
            np.arange(18.0).reshape(9, 2),
            np.arange(18.0, 36.0).reshape(9, 2),
            np.arange(36.0, 54.0).reshape(9, 2),
        )
        checkpoint = EvolutionCheckpoint(
            time=0.25,
            step_index=4,
            transaction_serial=20,
            state=state,
            metadata={"method": PRIMARY_METHOD, "control": True},
        )
        encoded = checkpoint.canonical_bytes()
        restored = EvolutionCheckpoint.from_bytes(encoded)
        self.assertEqual(encoded, restored.canonical_bytes())
        np.testing.assert_array_equal(state.u, restored.state.u)

        payload = json.loads(encoded)
        payload["state"]["u"]["sha256"] = "0" * 64
        tampered = (
            json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
        ).encode()
        with self.assertRaisesRegex(ValueError, "checksum"):
            EvolutionCheckpoint.from_bytes(tampered)


if __name__ == "__main__":
    unittest.main()
