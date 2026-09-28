"""Conditional continuous-residual and assembly arithmetic controls."""
from fractions import Fraction

import numpy as np
import pytest

from recursive_horizons.nsc_ks_continuous_constraint_bound import (
    CELL_COUNT, VERIFICATION_NODE_COUNT, componentwise_continuous_max,
    gamma_n, lipschitz_continuous_max, residual_assembly_roundoff,
    rounded_sum_error, verification_componentwise_continuous_max)


def test_lipschitz_cells_enclose_linear_function_and_report_remainder():
    nodes = (Fraction(0), Fraction(1, 2), Fraction(1))
    values = (Fraction(0), Fraction(1, 2), Fraction(1))
    result = lipschitz_continuous_max(nodes, values, (1, 1))
    assert result['continuous_maximum_upper'] == 1
    assert result['between_node_remainder_upper'] == 0
    assert not result['dense_sampling_used_as_bound']


def test_invalid_derivative_bound_is_rejected():
    with pytest.raises(ValueError, match='endpoint movement'):
        lipschitz_continuous_max((0, 1), (0, 2), (1,))


def test_componentwise_owner_keeps_N_and_beta_separate():
    values = np.array([[0, 1], [1, 0]], dtype=object)
    result = componentwise_continuous_max((0, 1), values,
                                          np.array([[1, 1]], dtype=object))
    assert set(result) == {'N', 'beta'}
    assert result['N']['continuous_maximum_upper'] == 1
    assert result['beta']['continuous_maximum_upper'] == 1


def test_roundoff_uses_directed_gamma_model():
    unit = Fraction(1, 1000)
    assert gamma_n(2, unit) == Fraction(2, 998)
    assert rounded_sum_error((1, 2, 3), unit_roundoff=unit) == Fraction(12, 998)
    result = residual_assembly_roundoff(
        {'N': (1, 2, 3, 4), 'beta': (2, 3)},
        unit_roundoff=unit)
    assert result['N'] == gamma_n(3, unit) * 10
    assert result['beta'] == gamma_n(1, unit) * 5


def test_verification_grid_requires_257_nodes_and_256_cells():
    nodes = tuple(Fraction(index, 256) for index in range(257))
    values = np.array([[Fraction(index, 256), Fraction(1) - Fraction(index, 256)]
                       for index in range(257)], dtype=object)
    derivatives = np.array([[1, 1]] * 256, dtype=object)
    result = verification_componentwise_continuous_max(nodes, values, derivatives)
    assert VERIFICATION_NODE_COUNT == 257
    assert CELL_COUNT == 256
    assert result['N']['continuous_maximum_upper'] == 1
    assert result['beta']['continuous_maximum_upper'] == 1
    with pytest.raises(ValueError, match='257'):
        verification_componentwise_continuous_max(nodes[:-1], values[:-1], derivatives[:-1])
