#!/usr/bin/env python3
"""Dirac-free checks for the n=64 coefficient layout."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n64_truncated_iterate as N64I
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


def test_blocked_step_lands_in_family_order():
    delta = np.arange(128, dtype=float)
    step = N64I.blocked_to_family(delta)
    assert step.shape == (2, 64)
    assert np.array_equal(step[0, :32], delta[0:32])
    assert np.array_equal(step[0, 32:], delta[64:96])
    assert np.array_equal(step[1, :32], delta[32:64])
    assert np.array_equal(step[1, 32:], delta[96:128])


def test_zero_pad_keeps_the_n32_radius_bound():
    low = LocalIncomingFamily(np.zeros((2, 32)))
    high = N64I.LocalIncomingFamily64(np.zeros((2, 64)))
    assert high.radius_lower_bound() == low.radius_lower_bound()
    assert high.radius_lower_bound() > 0


def test_unit_clip_scales_only_outside_the_ball():
    direction = np.zeros(4)
    direction[0] = 0.25
    _clipped, step, raw_norm, step_norm = N64I.clip_and_scale(direction, 1.0)
    assert raw_norm == 0.25 and step_norm == 0.25
    assert np.array_equal(step, direction)
    large = np.zeros(4)
    large[1] = 4.0
    _clipped, step, raw_norm, step_norm = N64I.clip_and_scale(large, 0.5)
    assert raw_norm == 4.0 and step_norm == 0.5
    assert step[1] == 0.5
