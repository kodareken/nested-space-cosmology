"""Audit the RN nodal stencils against a stabilized Fornberg construction.

`nsc_rn_pg._fornberg` builds finite-difference weights with the standard
Fornberg recursion and then verifies them against the polynomial moments
at a `1e-8` relative tolerance. On this runtime that verification fails
for every width of 9 or more, and the widths actually configured are
`_ONE_SIDED_D1_WIDTH = 9` and `_ONE_SIDED_D2_WIDTH = 10`, which act at
the two nodes at each end. Those end nodes are the throat and the outer
boundary, so the failing stencils sit exactly where the regularity
verdict is read.

This module answers one question without changing anything: is the
measured throat steepness a property of the physics, or an artifact of
the same round-off that breaks the moment check?

Nothing here imports, edits, or replaces `nsc_rn_pg`. The stabilized
weights are built here, compared against the recorded ones, and reported.
The existing operators keep their bytes and their provenance.
"""
from __future__ import annotations

import math
from fractions import Fraction
from typing import Dict, List

import numpy as np

API_VERSION = "nsc-rn-stencil-audit-v1"

#: Tolerance the existing guard uses.
MOMENT_TOLERANCE = 1.0e-8
#: Widths the existing chart configures for its one-sided stencils.
CONFIGURED_ONE_SIDED = ((9, 1), (10, 2))
#: Widest stencil worth auditing.
MAX_WIDTH = 12
#: Highest derivative order worth auditing.
MAX_DERIVATIVE = 4


def fornberg_exact(width: int, derivative: int, at: float) -> np.ndarray:
    """Fornberg weights in exact rational arithmetic.

    The nodes are integers, so the exact weights are rational numbers and
    `Fraction` returns them with no round-off at all. Exact rationals are
    the only genuinely stabilized option here: `np.longdouble` is 64-bit
    on Windows, so widening the accumulator changes nothing, which this
    audit measured rather than assumed.
    """
    nodes = [Fraction(i) for i in range(width)]
    order = int(derivative)
    target = Fraction(at).limit_denominator(10 ** 9)
    coeff = [[[Fraction(0) for _ in range(width)] for _ in range(width)]
             for _ in range(order + 1)]
    coeff[0][0][0] = Fraction(1)
    c1 = Fraction(1)
    for column in range(1, width):
        c2 = Fraction(1)
        for row in range(column):
            gap = nodes[column] - nodes[row]
            c2 *= gap
            for level in range(min(order, column) + 1):
                previous = coeff[level][column - 1][row]
                lower = coeff[level - 1][column - 1][row] if level else Fraction(0)
                coeff[level][column][row] = (
                    (nodes[column] - target) * previous - level * lower) / gap
        for level in range(min(order, column) + 1):
            previous = coeff[level][column - 1][column - 1]
            lower = coeff[level - 1][column - 1][column - 1] if level else Fraction(0)
            coeff[level][column][column] = (c1 / c2) * (
                level * lower - (nodes[column - 1] - target) * previous
            )
        c1 = c2
    return np.array([float(v) for v in coeff[order][width - 1]], dtype=np.float64)


def moment_defect(weights: np.ndarray, derivative: int, at: float = 0.0) -> float:
    """Worst relative error against the polynomial moment conditions.

    The moments are taken about the node the stencil differentiates at,
    which is what makes the check meaningful for a one-sided window.
    """
    nodes = np.arange(weights.size, dtype=np.longdouble) - np.longdouble(at)
    order = int(derivative)
    worst = np.longdouble(0.0)
    for power in range(weights.size):
        moment = np.sum(weights * nodes ** power, dtype=np.longdouble)
        expected = (np.longdouble(float(math.factorial(order)))
                    if power == order else np.longdouble(0.0))
        worst = max(worst, abs(moment - expected) / (1.0 + abs(expected)))
    return float(worst)


def weight_comparison(width: int, derivative: int, at: float) -> Dict[str, object]:
    """Compare one stencil: moment defects and the weight difference."""
    extended = fornberg_exact(width, derivative, at).astype(np.longdouble)
    plain = _plain_weights(width, derivative, at).astype(np.longdouble)

    plain_defect = moment_defect(plain, derivative, at)
    extended_defect = moment_defect(extended, derivative, at)
    scale = float(np.max(np.abs(extended))) or 1.0
    difference = float(np.max(np.abs(plain - extended))) / scale
    return {
        "width": width,
        "derivative": derivative,
        "at": float(at),
        "moment_defect_plain": plain_defect,
        "moment_defect_extended": extended_defect,
        "plain_passes": plain_defect <= MOMENT_TOLERANCE,
        "extended_passes": extended_defect <= MOMENT_TOLERANCE,
        "max_relative_weight_difference": difference,
    }


def audit_stencils() -> Dict[str, object]:
    """Every stencil in the configured range, plus a verdict per row."""
    rows: List[Dict[str, object]] = []
    for width in range(5, MAX_WIDTH + 1):
        for derivative in range(0, MAX_DERIVATIVE + 1):
            for at in (0.0, 0.5, 1.0):
                rows.append(weight_comparison(width, derivative, at))

    configured: List[Dict[str, object]] = []
    for width, derivative in CONFIGURED_ONE_SIDED:
        for at in (0.0, 1.0):
            row = weight_comparison(width, derivative, at)
            row["configured"] = True
            configured.append(row)

    diverging = [r for r in configured
                 if r["max_relative_weight_difference"] > 0.0 and not r["extended_passes"]]
    return {
        "api_version": API_VERSION,
        "moment_tolerance": MOMENT_TOLERANCE,
        "max_width": MAX_WIDTH,
        "max_derivative": MAX_DERIVATIVE,
        "rows": rows,
        "configured_one_sided": configured,
        "rows_failing_plain": sum(1 for r in rows if not r["plain_passes"]),
        "rows_failing_extended": sum(1 for r in rows if not r["extended_passes"]),
        "configured_failing_plain": sum(1 for r in configured if not r["plain_passes"]),
        "configured_failing_extended": sum(1 for r in configured if not r["extended_passes"]),
        "stabilized_never_worse": all(
            r["moment_defect_extended"] <= r["moment_defect_plain"] + 1e-18 for r in rows),
    }


def slope_under_both_stencils(radius: np.ndarray, values: np.ndarray) -> Dict[str, float]:
    """Throat logarithmic slope using each stencil's first derivative.

    `values` is sampled on the uniform areal grid. The two end nodes are
    differentiated with the width-9 one-sided stencil, the same rule
    `nsc_rn_pg` uses, so this reproduces the recorded operators exactly.
    The mirrored right end uses the `at = 1` window, as the chart does.
    Comparing the two answers says directly whether the measured
    steepness survives a stabilized operator.
    """
    grid = np.asarray(radius, dtype=np.float64)
    sample = np.asarray(values, dtype=np.float64)
    step = float(grid[1] - grid[0])
    count = sample.size
    if count < 12:
        raise ValueError("need at least twelve samples for a width-9 one-sided window")

    width = min(CONFIGURED_ONE_SIDED[0][0], count)
    plain = np.vstack([_plain_weights(width, 1, 0.0), _plain_weights(width, 1, 1.0)])
    exact = np.vstack([fornberg_exact(width, 1, 0.0), fornberg_exact(width, 1, 1.0)])

    derivative_plain = np.empty(count)
    derivative_exact = np.empty(count)
    for index in range(count):
        if index < 2:
            window = sample[:width]
            derivative_plain[index] = float(np.dot(plain[0], window)) / step
            derivative_exact[index] = float(np.dot(exact[0], window)) / step
        elif index > count - 3:
            mirror = min(index - (count - 2), 1)
            window = sample[-width:][::-1]
            sign = -1.0 if mirror < 2 else 1.0
            derivative_plain[index] = sign * float(np.dot(plain[mirror], window)) / step
            derivative_exact[index] = sign * float(np.dot(exact[mirror], window)) / step
        else:
            derivative_plain[index] = (sample[index + 1] - sample[index - 1]) / (2.0 * step)
            derivative_exact[index] = derivative_plain[index]

    def log_slope(derivative: np.ndarray) -> float:
        return float((np.log(abs(derivative[2])) - np.log(abs(derivative[0])))
                     / (np.log(grid[2]) - np.log(grid[0])))

    slope_plain = log_slope(derivative_plain)
    slope_exact = log_slope(derivative_exact)
    return {
        "width": width,
        "slope_plain": slope_plain,
        "slope_exact": slope_exact,
        "slope_difference": slope_exact - slope_plain,
        "derivative_relative_difference": float(
            np.max(np.abs(derivative_exact - derivative_plain))
            / (np.max(np.abs(derivative_plain)) or 1.0)),
    }


def _plain_weights(width: int, derivative: int, at: float) -> np.ndarray:
    """The double-precision recursion, isolated for comparison."""
    nodes = np.arange(width, dtype=np.float64)
    order = int(derivative)
    coeff = np.zeros((order + 1, width, width), dtype=np.float64)
    coeff[0, 0, 0] = 1.0
    c1 = 1.0
    target = float(at)
    for column in range(1, width):
        c2 = 1.0
        for row in range(column):
            gap = nodes[column] - nodes[row]
            c2 *= gap
            for level in range(min(order, column) + 1):
                previous = coeff[level, column - 1, row]
                lower = coeff[level - 1, column - 1, row] if level else 0.0
                coeff[level, column, row] = ((nodes[column] - target) * previous - level * lower) / gap
        for level in range(min(order, column) + 1):
            previous = coeff[level, column - 1, column - 1]
            lower = coeff[level - 1, column - 1, column - 1] if level else 0.0
            coeff[level, column, column] = (c1 / c2) * (
                level * lower - (nodes[column - 1] - target) * previous
            )
        c1 = c2
    return np.asarray(coeff[order, width - 1, :], dtype=np.float64)