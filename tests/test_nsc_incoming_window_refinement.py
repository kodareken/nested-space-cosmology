"""Explicit window binding, source conventions and reusable no-radial replay."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

import recursive_horizons.nsc_incoming_window_refinement as owner
from recursive_horizons.nsc_incoming_low_high_source import vacuum_source_kernel


ROOT = Path(__file__).resolve().parents[1]


def test_actual_two_window_endpoints_and_signed_inventory():
    for group, kind, left, right in ((12, 'low_vacuum', 32., 40.), (22, 'middle_correction', 40., 160.)):
        source = owner.authenticated_window(ROOT, group, kind, left, right)
        assert source['signs'] == (1, -1)
        assert source['measure_residual'] <= 3e-11
        assert all(row['interval'] == [left, right] for row in source['descriptors'])
        assert all(np.all(v['weights'] > 0) for v in source['selected'].values())
    with pytest.raises(ValueError, match='actual archived cell'):
        owner.authenticated_window(ROOT, 12, 'low_vacuum', 31., 40.)
    with pytest.raises(ValueError, match='actual archived endpoints'):
        owner.authenticated_window(ROOT, 22, 'middle_correction', 32., 160.)


def test_window_kernels_keep_existing_subtraction_and_pauli_vertices():
    with mp.workdps(80):
        coefficients = [mp.mpc(1, j+1)/(j+1) for j in range(24)]
        mass, ell, E = mp.mpf('1.2'), mp.mpf('-2.3'), mp.mpf(35)
        direct = owner.window_kernel(E, mass, ell, coefficients, 'low_vacuum')
        previous = vacuum_source_kernel(E, mass, ell, coefficients)
        assert list(direct) == list(previous)
        a, r = mp.sqrt(3*mp.pi/2-4), mp.sqrt(2)
        def projector(n):
            S = sum(coefficients[j]/(2*E)**(j+1) for j in range(n))
            v = mp.matrix([1, -1j*a*S]); return v*v.H/(v.H*v)[0]
        delta = projector(24)-projector(16)
        s1, s2, s3 = mp.matrix([[0,1],[1,0]]), mp.matrix([[0,-1j],[1j,0]]), mp.matrix([[1,0],[0,-1]])
        vertices = [-mass*s1+ell/r*s2-E/a*s3, -E/a*s3, E/a*mp.eye(2), ell/(2*r)*s2]
        expected = [sum((delta*V)[j,j] for j in range(2)) for V in vertices]
        actual = owner.window_kernel(E, mass, ell, coefficients, 'middle_correction')
        assert max(abs(x-y) for x, y in zip(actual, expected)) < mp.mpf('1e-74')
        assert actual[2] == 0.


def test_radial_interval_view_never_mutates_original_certificate():
    record = json.loads((ROOT/'results/development/nsc-incoming-retained-order24-bound.json').read_text())
    channel = json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels'][12]
    radial = record['radial']; before = deepcopy(radial)
    view = owner.bind_window_certificate(channel, radial, record['local_cross'], 32., 40.)
    assert radial == before
    assert view['original_endpoint_tag'] == 320.
    assert view['consumer_endpoint_tag'] == 40.
    assert view['fresh_interval'] == [32., 40.]
    assert view['coefficient_integrals_unchanged'] and not view['new_radial_recurrence']
    assert view['energy']['lapse_action_error_upper'] < 3e-11
    wrong = dict(radial, physical_Riccati_order=16)
    with pytest.raises(ValueError, match='matching order24'):
        owner.bind_window_certificate(channel, wrong, record['local_cross'], 32., 40.)


def test_two_case_record_has_explicit_before_after_and_thermal_policy():
    record = json.loads((ROOT/'results/development/nsc-incoming-window-reuse.json').read_text())
    assert set(record['windows']) == {'group12_low', 'group22_middle'}
    assert not record['failures']
    assert record['scope']['radial_recurrences_rerun'] == 0
    assert not record['scope']['thermal_source_physically_zero']
    for item in record['windows'].values():
        assert not item['failures'] and item['source_order'] == 24
        assert np.allclose(np.array(item['old_source'])+item['explicit_source_delta'], item['updated_source_approximant'], rtol=0, atol=1e-23)
        assert all(v > 0 for v in item['thermal_scattering_stress_upper'])
        assert item['physical_error_budget']['stationarity_tolerance'] == 3e-11
        assert not item['physical_error_budget']['rigorous_source_quadrature_error_included']
    assert record['windows']['group22_middle']['retained_archived_thermal'] is not None
    assert record['windows']['group12_low']['retained_archived_thermal'] is None


def test_replay_never_prepares_local_coefficients_or_quadrature(monkeypatch):
    def forbidden(*args, **kwargs): raise AssertionError('preparation called during window replay')
    monkeypatch.setattr(owner, '_riccati_at_one', forbidden); monkeypatch.setattr(owner, 'window_quadrature', forbidden)
    spec = importlib.util.spec_from_file_location('window_reuse_replay', ROOT/'scripts/derive_nsc_incoming_window_reuse.py')
    replay = importlib.util.module_from_spec(spec); spec.loader.exec_module(replay)
    prior = json.loads((ROOT/replay.OUTPUT).read_text())
    assert replay.make_record(ROOT/prior['payload']['path']) == prior
