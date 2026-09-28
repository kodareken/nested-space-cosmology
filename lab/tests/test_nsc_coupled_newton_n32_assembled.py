#!/usr/bin/env python3
"""Dirac-free tests for n=32 assembled remap/splice and zero-pad identity."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n32_assembled as N32
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


def test_zero_pad_identity_and_radius():
    trial = json.loads(N32.I4.OUTPUT.read_text())
    family = N32.zero_pad(trial['history']['coefficients'])
    assert profile_identity(family, include_normal_window=True) == N32.EXPECTED_PADDED
    parent = LocalIncomingFamily(np.asarray(trial['history']['coefficients'], float))
    assert family.radius_lower_bound() == parent.radius_lower_bound()
    assert np.array_equal(np.asarray(family.coefficients)[:, :16], parent.coefficients)


def test_remap_and_splice_round_trip():
    rng = np.random.default_rng(0)
    low = rng.normal(size=(32, 47, 2))
    high = rng.normal(size=(32, 47, 2))
    remapped = N32.remap_low_tangent(low)
    assert remapped.shape == (64, 47, 2)
    assert np.array_equal(remapped[0:16], low[0:16])
    assert np.array_equal(remapped[32:48], low[16:32])
    assert np.all(remapped[16:32] == 0) and np.all(remapped[48:64] == 0)
    spliced = N32.splice_high_tangent(remapped, high)
    assert np.array_equal(spliced[16:32], high[0:16])
    assert np.array_equal(spliced[48:64], high[16:32])
    assert np.array_equal(spliced[0:16], low[0:16])
    assert np.array_equal(spliced[32:48], low[16:32])


def test_high_mode_metric_matches_parent_radius_samples():
    trial = json.loads(N32.I4.OUTPUT.read_text())
    family = N32.zero_pad(trial['history']['coefficients'])
    metric = N32.high_mode_metric(family)
    assert len(metric.directions) == N32.FROZEN_PLUS_HIGH
    from recursive_horizons.nsc_ks_source_envelope import (
        computational_z_grid, _sample_axial_profiles)
    z = computational_z_grid(64)
    w_full, u_full = _sample_axial_profiles(family.metric().directions, z)
    w_high, u_high = _sample_axial_profiles(metric.directions, z)
    full = np.asarray(family.metric().amplitudes) @ w_full
    high = np.asarray(metric.amplitudes) @ w_high
    assert np.max(np.abs(full - high)) < 1e-14
    full_u = np.asarray(family.metric().amplitudes) @ u_full
    high_u = np.asarray(metric.amplitudes) @ u_high
    assert np.max(np.abs(full_u - high_u)) < 1e-14
