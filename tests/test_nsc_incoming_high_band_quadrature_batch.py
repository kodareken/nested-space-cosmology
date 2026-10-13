"""Complete direct-source coverage, interval aggregation and resumable replay."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import mpmath as mp
import pytest

import recursive_horizons.nsc_incoming_high_band_quadrature_batch as owner
from recursive_horizons.nsc_incoming_source_quadrature_bound import pack
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _range


ROOT = Path(__file__).resolve().parents[1]


from recursive_horizons.nsc_frozen_sources import frozen_source_context
# Authenticate original source bytes explicitly; numerical outputs remain local.
FROZEN_PRODUCER_CONTEXT = "329ef809e466a86ed4a290ee32c19485f85add9d"


@frozen_source_context(FROZEN_PRODUCER_CONTEXT)
def test_original_high_union_rejects_changes_beyond_its_weighted_roundoff():
    current = owner.authenticated_high_inventory(str(ROOT))[6]
    bounds = owner._selected_high_data(str(ROOT))[1][6]
    altered = deepcopy(current)
    altered['original_high_source_per_sign']['1'][0] += float(8*bounds['per_sign']['1'][0])
    with pytest.raises(ValueError, match='arithmetic envelope'):
        owner.require_original_high_union(ROOT, altered, current)
    altered = deepcopy(current)
    altered['channel']['copy_count'] += 1
    with pytest.raises(ValueError, match='provenance changed'):
        owner.require_original_high_union(ROOT, altered, current)


def test_current_high_inventory_cannot_reuse_a_frozen_cache(monkeypatch):
    calls = []
    def checked(root, context):
        calls.append(context)
        return ({}, {})
    monkeypatch.setattr(owner, '_authenticated_high_data', checked)
    # Supply the same uncached entrypoint used by the lru wrapper.
    checked.__wrapped__ = checked
    with frozen_source_context(FROZEN_PRODUCER_CONTEXT):
        owner.authenticated_high_inventory(str(ROOT))
    owner.authenticated_high_inventory(str(ROOT))
    assert calls == [FROZEN_PRODUCER_CONTEXT, None]

@frozen_source_context(FROZEN_PRODUCER_CONTEXT)
def test_original_source_union_is_complete_and_nonoverlapping():
    all_groups = owner.authenticated_high_inventory(str(ROOT))
    assert set(all_groups) == set(range(1, 33))
    for group, data in all_groups.items():
        left = 16 if group in (13, 14) else 32
        right = 320 if group in (10, 11, 12, 31, 32) else 160
        assert data['interval'] == [left, right]
        for sign in data['actual_signs']:
            intervals = sorted(p['interval'] for p in data['selected_original_panels'] if p['sign'] == sign)
            assert intervals == ([[16, 160]] if group in (13, 14) else [[32, 40], [40, right]])
        assert not data['channel']['index'] == 0


@frozen_source_context(FROZEN_PRODUCER_CONTEXT)
def test_correction_or_incomplete_window_cannot_become_full_direct_receipt():
    owned = owner.authenticated_high_inventory(str(ROOT))[6]
    record = json.loads((ROOT/'results/development/nsc-incoming-source-quadrature-bound.json').read_text())
    payload = json.loads((ROOT/record['payload']['path']).read_text())
    prepared = {'purpose': 'direct full source', 'owned_union': owned,
                'direct': payload['windows']['group6_middle_correction']}
    with pytest.raises(ValueError, match='complete direct'):
        owner.replay_direct_group(ROOT, prepared)
    prepared['direct'] = payload['windows']['group6_low']
    with pytest.raises(ValueError, match='complete direct'):
        owner.replay_direct_group(ROOT, prepared)


def test_aggregate_encloses_values_and_does_not_relax_tolerance():
    with _precision(80):
        groups = {str(g): {'source_integral_interval': [pack(_range(mp.mpf('1e-8')-mp.mpf('1e-18'), mp.mpf('1e-8')+mp.mpf('1e-18')))]*4,
                          'original_high_source': [2e-8]*4} for g in owner.GROUPS}
        result = owner.aggregate_groups(groups)
        assert result['numerical_budget_pass']
        assert abs(result['direct_high_vacuum_source'][0]-32e-8) < 1e-22
        assert abs(result['replacement_delta_from_original_high_source'][0]+32e-8) < 1e-22
        assert result['stress_error_upper'][0] >= 32e-18
        del groups['32']
        with pytest.raises(ValueError, match='exactly32'):
            owner.aggregate_groups(groups)
        groups['32'] = {'source_integral_interval': [pack(_range(-mp.mpf('1e-4'), mp.mpf('1e-4')))]*4,
                        'original_high_source': [0.]*4}
        assert not owner.aggregate_groups(groups)['numerical_budget_pass']


@frozen_source_context(FROZEN_PRODUCER_CONTEXT)
def test_complete_batch_replay_never_generates_source(monkeypatch):
    path = ROOT/'results/development/nsc-incoming-high-band-quadrature-batch.json'
    if not path.exists():
        pytest.skip('full directed batch not produced yet')
    record = json.loads(path.read_text())
    assert len(record['group_payloads']) == 32
    assert not record['scope']['old16_baseline_quadrature_carried']
    assert not record['scope']['thermal_or_physical_projector_error_included']
    def forbidden(*args, **kwargs):
        raise AssertionError('source evaluation during replay')
    monkeypatch.setattr(owner, 'certify_source_window', forbidden)
    monkeypatch.setattr(owner, 'cached_rule', forbidden)
    spec = importlib.util.spec_from_file_location('direct_quad_replay', ROOT/'scripts/derive_nsc_incoming_high_band_quadrature_batch.py')
    replay = importlib.util.module_from_spec(spec); spec.loader.exec_module(replay)
    assert replay.make_record(record['group_payloads']) == record
    missing = deepcopy(record['group_payloads'])[:-1]
    with pytest.raises(ValueError, match='exact32'):
        replay.make_record(missing)
