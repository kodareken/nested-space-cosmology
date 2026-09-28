"""Prospective, noncommitting compensated RK arithmetic for TDG11.

The mathematical tableaux are unchanged. The binary64 operation schedule is
new and explicitly identified; this module is neither a production authority
nor an error/convergence certificate. Unsupported arithmetic fails closed.
"""

from __future__ import annotations

from fractions import Fraction
import math
from numbers import Real
from typing import Sequence

import numpy as np

from .numerical_engine import (
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    RightHandSide,
    StageRecord,
    StateProjector,
    StepProposal,
)


ARITHMETIC_ID = "tdg11_binary64_twofold_dot2_coherent_rk_v1"
_TINY = np.finfo(np.float64).tiny
_SPLITTER = float(2**27 + 1)
_MAX_TERMS = 4


class CompensatedArithmeticStop(RuntimeError):
    """A premise stop, never a scientific nonpass or a precision fallback."""

    def __init__(self, reason: str) -> None:
        self.reason = reason
        super().__init__(f"TDG11 compensated arithmetic: {reason}")


def _checked(value: np.ndarray) -> np.ndarray:
    if not np.all(np.isfinite(value)):
        raise CompensatedArithmeticStop("nonfinite_or_overflow")
    magnitude = np.abs(value)
    if np.any((magnitude != 0.0) & (magnitude < _TINY)):
        raise CompensatedArithmeticStop("unsupported_underflow")
    return value


def _operand(value: object) -> np.ndarray:
    if type(value) in (float, np.float64):
        return _checked(np.asarray(value, dtype=np.float64))
    if type(value) is not np.ndarray or value.dtype != np.dtype(np.float64):
        raise TypeError("arithmetic operands must be native binary64 arrays or floats")
    return _checked(value)


def _op(operation, left: np.ndarray, right: np.ndarray) -> np.ndarray:
    try:
        with np.errstate(over="raise", under="raise", invalid="raise"):
            answer = operation(left, right)
    except FloatingPointError as error:
        reason = (
            "unsupported_underflow"
            if "underflow" in str(error)
            else "nonfinite_or_overflow"
        )
        raise CompensatedArithmeticStop(reason) from error
    if operation is np.multiply and np.any(
        (left != 0.0) & (right != 0.0) & (answer == 0.0)
    ):
        raise CompensatedArithmeticStop("unsupported_underflow")
    return _checked(answer)


def _add(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return _op(np.add, left, right)


def _sub(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return _op(np.subtract, left, right)


def _mul(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return _op(np.multiply, left, right)


def two_sum(a: object, b: object) -> tuple[np.ndarray, np.ndarray]:
    """Knuth's six-operation error-free sum on the supported normal domain."""

    left, right = np.broadcast_arrays(_operand(a), _operand(b))
    total = _add(left, right)
    b_virtual = _sub(total, left)
    a_virtual = _sub(total, b_virtual)
    b_error = _sub(right, b_virtual)
    a_error = _sub(left, a_virtual)
    error = _add(a_error, b_error)
    return total, error


def _split(value: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    expanded = _mul(np.asarray(_SPLITTER), value)
    high = _sub(expanded, _sub(expanded, value))
    low = _sub(value, high)
    return high, low


def two_product(a: object, b: object) -> tuple[np.ndarray, np.ndarray]:
    """Dekker product; no silent scaling, FMA, or underflow fallback."""

    left, right = np.broadcast_arrays(_operand(a), _operand(b))
    product = _mul(left, right)
    left_high, left_low = _split(left)
    right_high, right_low = _split(right)
    error1 = _sub(product, _mul(left_high, right_high))
    error2 = _sub(error1, _mul(left_low, right_high))
    error3 = _sub(error2, _mul(left_high, right_low))
    error = _sub(_mul(left_low, right_low), error3)
    return product, error


def twofold_coefficient(value: Fraction) -> tuple[float, float, Fraction]:
    """Return two rounded words and their *exact*, possibly nonzero remainder."""

    if type(value) is not Fraction:
        raise TypeError("coefficient must be an exact Fraction")
    try:
        high = float(value)
        if not math.isfinite(high):
            raise CompensatedArithmeticStop("nonfinite_or_overflow")
        remainder = value - Fraction.from_float(high)
        low = float(remainder)
    except OverflowError as error:
        raise CompensatedArithmeticStop("nonfinite_or_overflow") from error
    _operand(high)
    _operand(low)
    if (value and high == 0.0) or (remainder and low == 0.0):
        raise CompensatedArithmeticStop("unsupported_underflow")
    return high, low, remainder - Fraction.from_float(low)


def compensated_update(
    base: np.ndarray,
    terms: Sequence[tuple[Fraction, np.ndarray]],
) -> np.ndarray:
    """A fixed twofold-coefficient Dot2-style base-plus-increments update.

    The rational coefficient remainder beyond two words is not asserted zero.
    Exact-rational controls expose it. The measured output remains binary64.
    """

    if type(base) is not np.ndarray:
        raise TypeError("base must be a binary64 array")
    _operand(base)
    if not isinstance(terms, (tuple, list)) or not 1 <= len(terms) <= _MAX_TERMS:
        raise ValueError("an update requires one to four ordered RHS terms")
    checked_terms: list[tuple[tuple[float, float, Fraction], np.ndarray]] = []
    for term in terms:
        if not isinstance(term, (tuple, list)) or len(term) != 2:
            raise TypeError("each term must contain a Fraction and a binary64 array")
        coefficient, data = term
        if type(data) is not np.ndarray or data.shape != base.shape:
            raise ValueError("RHS term shape differs from the base")
        _operand(data)
        checked_terms.append((twofold_coefficient(coefficient), data))
    accumulator = base.copy()
    correction = np.zeros_like(base)
    for (high, low, _remainder), data in checked_terms:
        for part in (high, low):
            product, product_error = two_product(part, data)
            accumulator, sum_error = two_sum(accumulator, product)
            correction = _add(correction, _add(sum_error, product_error))
    return _add(accumulator, correction)


def _real(value: Real, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real")
    answer = float(value)
    if not math.isfinite(answer):
        raise CompensatedArithmeticStop("nonfinite_or_overflow")
    return answer


def propose_compensated_step(
    *,
    method: str,
    time: Real,
    step_size: Real,
    state: EvolutionState,
    rhs: RightHandSide,
    projector: StateProjector | None = None,
) -> StepProposal:
    """Collect stages and a fresh endpoint; perform no acceptance or I/O.

    ``StepProposal.method`` retains the underlying tableau selector required
    by historical stage guards. Callers must additionally bind ARITHMETIC_ID;
    it is not an assertion of byte-equivalence with the historical proposer.
    """

    if method not in (PRIMARY_METHOD, COMPARATOR_METHOD):
        raise ValueError("unknown evolution method")
    if not isinstance(state, EvolutionState):
        raise TypeError("state must be EvolutionState")
    start, dt = _real(time, "time"), _real(step_size, "step_size")
    if dt <= 0.0:
        raise ValueError("step_size must be positive")
    final, middle = start + dt, start + dt / 2.0
    if not math.isfinite(final) or not start < middle < final:
        raise CompensatedArithmeticStop("unresolved_stage_time_lattice")
    width = Fraction.from_float(dt)
    records: list[StageRecord] = []

    def project(at: float, value: EvolutionState) -> EvolutionState:
        result = value if projector is None else projector(at, value)
        if not isinstance(result, EvolutionState) or result.shape != state.shape:
            raise ValueError("projector returned an incompatible state")
        for array in (result.u, result.p, result.q):
            _operand(array)
        return result

    def evaluate(name: str, at: float, value: EvolutionState) -> EvolutionRHS:
        projected = project(at, value)
        derivative = rhs(at, projected)
        if not isinstance(derivative, EvolutionRHS) or derivative.shape != state.shape:
            raise ValueError("RHS returned an incompatible value")
        for array in (derivative.du, derivative.dp, derivative.dq):
            _operand(array)
        records.append(StageRecord(name, at, projected, derivative))
        return derivative

    def update(terms: tuple[tuple[Fraction, EvolutionRHS], ...]) -> EvolutionState:
        return EvolutionState(
            *(
                compensated_update(
                    getattr(state, block),
                    tuple(
                        (width * coefficient, getattr(value, derivative))
                        for coefficient, value in terms
                    ),
                )
                for block, derivative in (("u", "du"), ("p", "dp"), ("q", "dq"))
            )
        )

    if method == PRIMARY_METHOD:
        k1 = evaluate("rk4_k1", start, state)
        k2 = evaluate(
            "rk4_k2", middle, project(middle, update(((Fraction(1, 2), k1),)))
        )
        k3 = evaluate(
            "rk4_k3", middle, project(middle, update(((Fraction(1, 2), k2),)))
        )
        k4 = evaluate("rk4_k4", final, project(final, update(((Fraction(1), k3),))))
        candidate = update(
            (
                (Fraction(1, 6), k1),
                (Fraction(1, 3), k2),
                (Fraction(1, 3), k3),
                (Fraction(1, 6), k4),
            )
        )
    else:
        k1 = evaluate("ssprk3_s0", start, state)
        k2 = evaluate("ssprk3_s1", final, project(final, update(((Fraction(1), k1),))))
        k3 = evaluate(
            "ssprk3_s2",
            middle,
            project(middle, update(((Fraction(1, 4), k1), (Fraction(1, 4), k2)))),
        )
        candidate = update(
            ((Fraction(1, 6), k1), (Fraction(1, 6), k2), (Fraction(2, 3), k3))
        )
    candidate = project(final, candidate)
    evaluate("candidate_endpoint", final, candidate)
    return StepProposal(method, start, final, state, candidate, tuple(records))


__all__ = [
    "ARITHMETIC_ID",
    "CompensatedArithmeticStop",
    "compensated_update",
    "propose_compensated_step",
    "two_product",
    "two_sum",
    "twofold_coefficient",
]
