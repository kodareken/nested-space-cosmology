"""Exact sign reuse and authenticated nonoverlapping batch accounting."""
import importlib.util
import json
from pathlib import Path

import pytest

from recursive_horizons.nsc_incoming_high_band_batch_bound import (
    GROUPS, CachedGeometry, authenticated_batch, sign_symmetry_identity, interval_sign_control,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/'results/development/nsc-incoming-high-band-batch-bound.json'


def test_exact_induction_and_actual_interval_sign_control():
    result = sign_symmetry_identity()
    assert result['base_residual'] == result['norm_squared_residual'] == '0'
    assert result['recurrence_residuals'] == ['0']*23
    assert result['defect_residuals'] == ['0']*25
    assert interval_sign_control()['maximum_endpoint_difference'] == 0.
    assert not result['physical_occupation_symmetry_assumed']


def test_actual28_groups55_windows_and_no_group13_overlap():
    data = authenticated_batch(ROOT)
    assert len(GROUPS) == 28 and set(GROUPS).isdisjoint({12, 14, 22, 32})
    assert sum(len(v) for v in data['windows'].values()) == 55
    assert [r['interval'] for r in data['windows']['13']] == [[16., 160.]]
    assert [r['interval'] for r in data['windows']['6']] == [[32., 40.], [40., 160.]]
    for group, windows in data['windows'].items():
        if group != '13': assert windows[0]['interval'][1] == windows[1]['interval'][0]


def test_cached_geometry_does_not_allow_completed_group_reruns():
    # Rejection occurs before reading any cached geometry or coefficient.
    obj = object.__new__(CachedGeometry)
    for group in (0, 12, 14, 22, 32):
        with pytest.raises(ValueError, match='completed'):
            obj.bound({'index': group})


def test_archived_signed_norm_arrays_and_saved_only_replay(monkeypatch):
    spec = importlib.util.spec_from_file_location('highband', ROOT/'scripts/derive_nsc_incoming_high_band_batch_bound.py')
    script = importlib.util.module_from_spec(spec); spec.loader.exec_module(script)
    assert script.archived_sign_check()['norm_array_difference'] == 0.
    if not OUTPUT.exists(): pytest.skip('post-preparation replay gate')
    def forbidden(*a, **k): raise AssertionError('replay generated interval coefficients')
    monkeypatch.setattr(script, 'interval_sign_control', forbidden)
    monkeypatch.setattr(script, 'local_cross_product_coefficients', forbidden)
    monkeypatch.setattr(CachedGeometry, 'bound', forbidden)
    prior = json.loads(OUTPUT.read_text()); actual = script.replay(prior['manifest'])
    assert actual == prior
    assert actual['group_count'] == 28 and actual['source_window_count'] == 55
    assert actual['combined_intervals_added_once_per_group']
    assert not actual['scope']['occupation_symmetry_assumed']
