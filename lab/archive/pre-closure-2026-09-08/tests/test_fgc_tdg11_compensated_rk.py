"""Outcome-blind arithmetic controls; no campaign inputs or accepted state."""

from __future__ import annotations

from fractions import Fraction as Q
import math
import unittest

import numpy as np

from recursive_horizons.fgc.evolution import tdg11_compensated_rk as arithmetic
from recursive_horizons.fgc.evolution.numerical_engine import (
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
)


def _exact(value) -> Q:
    return Q.from_float(float(value))


def _state() -> EvolutionState:
    return EvolutionState(*(np.full((9, 6), value) for value in (1.0, 2.0, 4.0)))


class CompensatedPrimitiveTests(unittest.TestCase):
    def test_two_sum_exact_fraction_oracle(self):
        for left, right in (
            (1.0, 2.0**-54),
            (2.0**53, -(2.0**53) + 1),
            (-1.25, 0.125),
            (0.0, 0.0),
            (2.0**100, -(2.0**-100)),
        ):
            high, low = arithmetic.two_sum(left, right)
            self.assertEqual(_exact(high) + _exact(low), _exact(left) + _exact(right))

    def test_two_product_exact_fraction_oracle(self):
        left = np.array([1.0 + 2**-52, 1.0 / 3, -1.25, 2.0**100, 0.0])
        right = np.array([1.0 - 2**-52, 2.0 / 3, 0.375, 2.0**-80, 17.0])
        high, low = arithmetic.two_product(left, right)
        for a, b, p, e in zip(left, right, high, low, strict=True):
            self.assertEqual(_exact(p) + _exact(e), _exact(a) * _exact(b))

    def test_deterministic_wide_exponent_oracles(self):
        # All data precede any campaign measurement and are exactly specified.
        for exponent in (-180, -40, 0, 35, 180):
            for mantissa in (1.0 + 2**-52, -1.5, 1.0 - 2**-53):
                a, b = math.ldexp(mantissa, exponent), math.ldexp(1.25, -exponent)
                for operation, exact in (
                    (arithmetic.two_sum, _exact(a) + _exact(b)),
                    (arithmetic.two_product, _exact(a) * _exact(b)),
                ):
                    high, low = operation(a, b)
                    self.assertEqual(_exact(high) + _exact(low), exact)

    def test_coefficient_remainder_is_not_erased(self):
        for value in (Q(0), Q(1, 4), Q(1, 3), Q(-1, 6), Q.from_float(0.0001) / 6):
            high, low, rest = arithmetic.twofold_coefficient(value)
            self.assertEqual(_exact(high) + _exact(low) + rest, value)
        self.assertNotEqual(arithmetic.twofold_coefficient(Q(1, 3))[2], 0)

    def test_compensation_recovers_cancellation(self):
        base = np.array([2.0**53])
        terms = ((Q(1), np.array([1.0])), (Q(1), np.array([-(2.0**53)])))
        self.assertEqual(float((base + terms[0][1] + terms[1][1])[0]), 0.0)
        self.assertEqual(float(arithmetic.compensated_update(base, terms)[0]), 1.0)

    def test_no_input_mutation(self):
        base = np.array([1.0, 2.0])
        data = np.array([0.5, -0.25])
        before = (base.tobytes(), data.tobytes())
        answer = arithmetic.compensated_update(base, ((Q(1, 3), data),))
        self.assertEqual(before, (base.tobytes(), data.tobytes()))
        self.assertFalse(np.shares_memory(answer, base))
        self.assertFalse(np.shares_memory(answer, data))

    def test_overflow_and_unsupported_underflow_are_typed(self):
        for operation in (
            lambda: arithmetic.two_sum(float("inf"), 1.0),
            lambda: arithmetic.two_sum(np.finfo(float).max, np.finfo(float).max),
            lambda: arithmetic.two_product(2.0**1000, 1.0),
            lambda: arithmetic.two_product(2.0**-600, 2.0**-600),
            lambda: arithmetic.two_sum(float.fromhex("0x0.0000000000001p-1022"), 0.0),
            lambda: arithmetic.twofold_coefficient(Q(1, 2**1200)),
            lambda: arithmetic.twofold_coefficient(Q(2**1200)),
        ):
            with (
                self.subTest(operation=operation),
                self.assertRaises(arithmetic.CompensatedArithmeticStop),
            ):
                operation()

    def test_invalid_contract_inputs_are_rejected(self):
        for value in (
            True,
            1,
            np.array([1], dtype=np.float32),
            np.array([1], dtype=object),
        ):
            with self.assertRaises(TypeError):
                arithmetic.two_sum(value, 1.0)
        for coefficient in (True, 0.5, 1):
            with self.assertRaises(TypeError):
                arithmetic.twofold_coefficient(coefficient)
        with self.assertRaises(ValueError):
            arithmetic.compensated_update(np.ones(2), ((Q(1), np.ones(3)),))
        for terms in ((), ((Q(1), np.ones(2)),) * 5):
            with self.assertRaises(ValueError):
                arithmetic.compensated_update(np.ones(2), terms)


class CompensatedProposalTests(unittest.TestCase):
    def test_fresh_endpoint_stage_schedule_and_projector(self):
        expected = {
            PRIMARY_METHOD: (
                ("rk4_k1", "rk4_k2", "rk4_k3", "rk4_k4", "candidate_endpoint"),
                (0.0, 0.125, 0.125, 0.25, 0.25),
            ),
            COMPARATOR_METHOD: (
                ("ssprk3_s0", "ssprk3_s1", "ssprk3_s2", "candidate_endpoint"),
                (0.0, 0.25, 0.125, 0.25),
            ),
        }
        for method, (names, times) in expected.items():
            calls, projections = [], []
            state = _state()
            before = tuple(getattr(state, name).tobytes() for name in ("u", "p", "q"))

            def rhs(at, value):
                calls.append((at, value))
                return EvolutionRHS(value.u, value.p, value.q, {})

            def projector(at, value):
                projections.append(at)
                return value

            proposal = arithmetic.propose_compensated_step(
                method=method,
                time=0.0,
                step_size=0.25,
                state=state,
                rhs=rhs,
                projector=projector,
            )
            self.assertEqual(
                tuple(record.stage_name for record in proposal.stages), names
            )
            self.assertEqual(tuple(record.time for record in proposal.stages), times)
            self.assertEqual(len(calls), len(names))
            self.assertEqual(projections[-2:], [0.25, 0.25])
            self.assertIs(calls[-1][1], proposal.candidate_state)
            self.assertIsNot(proposal.stages[-1].rhs, proposal.stages[-2].rhs)
            self.assertEqual(
                before,
                tuple(getattr(state, name).tobytes() for name in ("u", "p", "q")),
            )

    def test_zero_rhs_preserves_state(self):
        state = _state()

        def zero(_time, value):
            return EvolutionRHS(*(np.zeros_like(value.u) for _ in range(3)), {})

        for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
            proposal = arithmetic.propose_compensated_step(
                method=method, time=0.0, step_size=0.125, state=state, rhs=zero
            )
            for block in ("u", "p", "q"):
                np.testing.assert_array_equal(
                    getattr(state, block), getattr(proposal.candidate_state, block)
                )

    def test_exact_linear_polynomial_control(self):
        state = _state()

        def rhs(_time, value):
            return EvolutionRHS(value.u, value.p, value.q, {})

        for method, order in ((PRIMARY_METHOD, 4), (COMPARATOR_METHOD, 3)):
            for width in (Q(1, 8), Q(1, 16), Q(1, 32)):
                factor = sum(
                    (
                        width**power / math.factorial(power)
                        for power in range(order + 1)
                    ),
                    Q(0),
                )
                proposal = arithmetic.propose_compensated_step(
                    method=method,
                    time=0.0,
                    step_size=float(width),
                    state=state,
                    rhs=rhs,
                )
                for block, value in (("u", Q(1)), ("p", Q(2)), ("q", Q(4))):
                    observed = float(getattr(proposal.candidate_state, block)[0, 0])
                    # A control of this finite polynomial, not a universal bound.
                    self.assertLessEqual(
                        abs(_exact(observed) - factor * value),
                        _exact(math.ulp(observed)),
                    )

    def test_wrong_method_rhs_projector_and_time_fail(self):
        state = _state()

        def rhs(_time, value):
            return EvolutionRHS(value.u, value.p, value.q, {})

        values = dict(
            method=PRIMARY_METHOD, time=0.0, step_size=0.125, state=state, rhs=rhs
        )
        for changes in (
            {"method": "RK4"},
            {"step_size": 0},
            {"time": True},
            {"rhs": lambda *_args: None},
            {"projector": lambda *_args: None},
            {"state": object()},
        ):
            with self.assertRaises((TypeError, ValueError)):
                arithmetic.propose_compensated_step(**(values | changes))
        with self.assertRaises(arithmetic.CompensatedArithmeticStop):
            arithmetic.propose_compensated_step(**(values | {"time": 2.0**53}))


if __name__ == "__main__":
    unittest.main()
