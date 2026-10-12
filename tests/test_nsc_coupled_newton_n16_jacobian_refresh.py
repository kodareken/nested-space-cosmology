"""The n=16 leftover stall does not authorize EXISTENCE or n=32."""
import importlib
from pathlib import Path

import numpy as np


def _module(name):
    scripts = str(Path(__file__).resolve().parents[1] / 'scripts')
    import sys
    if scripts not in sys.path:
        sys.path.insert(0, scripts)
    return importlib.import_module(name)


def test_unclipped_leftover_stops_the_frozen_walk():
    next5 = _module('derive_nsc_ks_coupled_newton_n16_damped_next5')
    damped = _module('derive_nsc_ks_coupled_newton_n16_damped')
    direction = np.ones(32) / np.sqrt(32)
    gradient = np.full((5, 2), 0.004)
    jacobian = np.zeros((32, 5, 2))
    leftover = next5.unclipped_leftover(gradient, jacobian, direction, 14.8)
    assert leftover[0] > next5.EXISTENCE_TOLERANCE
    assert leftover[1] > next5.EXISTENCE_TOLERANCE
    predicted = damped.predicted_maxima(gradient, jacobian, np.zeros(32))
    assert damped.beats_measured(np.array([0.0038, 0.0038]), predicted)
    assert not np.all(leftover <= next5.EXISTENCE_TOLERANCE)


def test_small_chebyshev_tail_does_not_open_n32():
    leftover = np.array([1.3e-4, 6.3e-4])
    tails = np.array([5.701e-9, 1.0634e-8])
    assert not np.any(tails > 0.1 * leftover)
