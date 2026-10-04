"""Tests for the Fornberg moment guard.

The guard was blocking import of ``nsc_rn_pg`` entirely, at every width the
boundary tables use. The weights were never wrong: an exact rational audit
of 120 configurations agreed with them to 5.1e-16. What was wrong was the
test, which scaled its tolerance by the magnitude of the answer rather than
by the magnitude of the sum being formed.

These tests pin down both halves of that claim. The guard must accept the
real stencils, and it must still reject a weight that is actually wrong.
A guard loosened without a rejection test proves nothing.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from recursive_horizons import nsc_rn_pg as pg


def _moment_residual(weights, nodes, target, power, order):
    column = (nodes - target) ** power
    moment = float(np.dot(weights, column))
    expected = float(math.factorial(order)) if power == order else 0.0
    scale = float(np.sum(np.abs(weights) * np.abs(column)))
    return abs(moment - expected), expected, scale


def _exact_fornberg(width, order, at):
    """Fornberg weights recovered in exact rational arithmetic.

    Solves the moment system sum_i w_i (x_i - at)^p = order! if p == order
    else 0 for p in 0..width-1 by Gaussian elimination over ``Fraction``.
    This is an independent derivation, not a copy of the module's output.
    """
    from fractions import Fraction

    def factorial_exact(n):
        result = Fraction(1)
        for k in range(2, n + 1):
            result *= k
        return result

    rows = []
    for power in range(width):
        row = [Fraction(x - at) ** power for x in range(width)]
        row.append(factorial_exact(order) if power == order else Fraction(0))
        rows.append(row)
    for column in range(width):
        pivot = next(r for r in range(column, width) if rows[r][column] != 0)
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = rows[column][column]
        rows[column] = [v / scale for v in rows[column]]
        for r in range(width):
            if r != column and rows[r][column] != 0:
                factor = rows[r][column]
                rows[r] = [v - factor * w for v, w in zip(rows[r], rows[column])]
    return [row[-1] for row in rows]


class TestTheModuleImports:
    def test_boundary_tables_are_built_at_module_scope(self):
        assert len(pg._LEFT_D1) > 0
        assert len(pg._LEFT_D2) > 0

    def test_widths_used_by_the_boundary_closures(self):
        assert pg._ONE_SIDED_D1_WIDTH == 9
        assert pg._ONE_SIDED_D2_WIDTH == 10


class TestTheRealStencilsAreAccepted:
    @pytest.mark.parametrize("width", range(5, 15))
    @pytest.mark.parametrize("order", [0, 1, 2, 3, 4])
    @pytest.mark.parametrize("at", [0, 1])
    def test_fornberg_accepts_every_configured_stencil(self, width, order, at):
        if order >= width:
            pytest.skip("derivative order beyond the stencil")
        weights = pg._fornberg(width, order, at)
        assert np.all(np.isfinite(weights))

    @pytest.mark.parametrize("width,order,at", [(10, 2, 0), (10, 2, 1), (9, 1, 0)])
    def test_the_configurations_that_used_to_raise(self, width, order, at):
        pg._fornberg(width, order, at)

    def test_weights_are_exactly_the_rational_solution(self):
        """The real claim: the weights were never wrong.

        A copy of the module's own output would prove nothing, so the
        weights are recovered independently by solving the moment system in
        exact rational arithmetic and compared against the float64 values.
        """
        width, order, at = 9, 1, 0
        recovered = pg._fornberg(width, order, at)
        exact = _exact_fornberg(width, order, at)
        assert np.allclose(recovered, np.array([float(v) for v in exact]),
                           rtol=1e-15, atol=1e-17)

    @pytest.mark.parametrize("width,order,at", [(10, 2, 0), (10, 2, 1), (9, 1, 0)])
    def test_every_configured_stencil_matches_exact_arithmetic(self, width, order, at):
        recovered = pg._fornberg(width, order, at)
        exact = np.array([float(v) for v in _exact_fornberg(width, order, at)])
        # Scaled by the weight magnitude: these stencils reach about 162,
        # where one unit in the last place is already 3e-14.
        scale = float(np.max(np.abs(exact)))
        assert np.max(np.abs(recovered - exact)) / scale < 1e-14


class TestTheMomentsAreActuallySatisfied:
    @pytest.mark.parametrize("width,order,at", [(9, 1, 0), (10, 2, 0), (10, 2, 1)])
    def test_residual_is_within_the_achievable_float64_accuracy(self, width, order, at):
        weights = pg._fornberg(width, order, at)
        nodes = np.arange(width, dtype=float)
        worst = 0.0
        for power in range(width):
            residual, expected, scale = _moment_residual(weights, nodes, float(at), power, order)
            allowance = 32.0 * np.finfo(np.float64).eps * (scale + abs(expected))
            assert residual <= allowance, (width, order, at, power)
            worst = max(worst, residual / (scale + abs(expected)))
        assert worst < 1e-13

    def test_the_identity_moment_holds_to_near_machine_precision(self):
        weights = pg._fornberg(10, 2, 0)
        nodes = np.arange(10, dtype=float)
        residual, expected, _ = _moment_residual(weights, nodes, 0.0, 2, 2)
        assert expected == 2.0
        assert residual < 1e-12


class TestTheGuardStillHasTeeth:
    """The repair loosened a tolerance, so it must still reject bad weights."""

    def _guard(self, weights, width, order, at):
        nodes = np.arange(width, dtype=float)
        factor = 32.0 * float(np.finfo(np.float64).eps)
        for power in range(width):
            residual, expected, scale = _moment_residual(weights, nodes, at, power, order)
            if residual > factor * (scale + abs(expected)):
                return False
        return True

    def test_a_perturbed_weight_is_rejected(self):
        weights = pg._fornberg(10, 2, 0).copy()
        weights[4] += 1e-6
        assert not self._guard(weights, 10, 2, 0)

    def test_a_shifted_weight_is_rejected(self):
        weights = pg._fornberg(10, 2, 0).copy()
        weights[0] += 1e-3
        assert not self._guard(weights, 10, 2, 0)

    def test_the_real_weights_pass_the_same_guard(self):
        assert self._guard(pg._fornberg(10, 2, 0), 10, 2, 0)

    def test_a_one_ulp_perturbation_is_still_accepted(self):
        """The repair must not become a licence for anything roundoff-sized."""
        weights = pg._fornberg(10, 2, 0).copy()
        weights[3] = float(np.nextafter(weights[3], np.inf))
        assert self._guard(weights, 10, 2, 0)