"""Actual windows, nonmutating energy rebinding and saved-only replay."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT/'results/development/nsc-incoming-group32-order24-bound.json'


def script():
    spec = importlib.util.spec_from_file_location('group32_bound', ROOT/'scripts/derive_nsc_incoming_group32_order24_bound.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def test_original32_guard_and_window_metadata_authenticate():
    module = script(); owned = module.authenticated_windows()
    assert owned['channel']['index'] == 32
    assert [(r['sign'], r['rows'], r['interval']) for r in owned['low']] == [(1, 32, [32., 40.]), (-1, 32, [32., 40.])]
    assert owned['middle']['left'] == 40. and owned['middle']['right'] == 320.
    assert owned['low_measure_residual'] <= 3e-13


def test_explicit_window_rebinding_preserves_coefficients_and_original():
    module = script(); values = [{'sign': 1, 'coefficient_integral_upper': [1.]*25}]
    radial = {'group': 32, 'lower': 320., 'physical_Riccati_order': 24,
              'intervals': 128, 'centered_remainder_depth': 4, 'per_sign': values}
    before = json.dumps(radial, sort_keys=True)
    for name, window in module.WINDOWS.items():
        rebound = module.bind_window(radial, name)
        assert rebound['per_sign'] is values
        assert rebound['original_radial_endpoint'] == 320.
        assert rebound['lower'] == window[1]
        assert rebound['energy_contraction_window'] == list(window)
    assert json.dumps(radial, sort_keys=True) == before
    with pytest.raises(ValueError, match='authenticated'):
        module.bind_window(radial, 'unowned')
    with pytest.raises(ValueError, match='original group32'):
        module.bind_window({**radial, 'lower': 40.}, 'middle')
    with pytest.raises(ValueError, match='original group32'):
        module.bind_window({**radial, 'physical_Riccati_order': 16}, 'low')


def test_saved_replay_does_not_prepare_radial_or_local_coefficients(monkeypatch):
    if not RESULT.exists(): pytest.skip('post-preparation replay gate')
    module = script()
    def forbidden(*a, **kw): raise AssertionError('replay must not regenerate coefficients')
    monkeypatch.setattr(module, 'retained_order24_bound', forbidden)
    monkeypatch.setattr(module, 'local_cross_product_coefficients', forbidden)
    old = json.loads(RESULT.read_text()); actual = module.replay(ROOT/old['payload']['path'])
    assert actual == old
    assert actual['radial_preparation_passes'] == 1
    assert actual['scope']['matching_order24_source_corrections_required']
    assert not actual['scope']['source_corrections_applied']
    assert actual['scope']['combined_is_alternative_to_sum_of_windows']
    assert actual['radial']['physical_Riccati_order'] == 24
    assert set(actual['windows']) == {'low', 'middle', 'combined'}
