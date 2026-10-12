"""Independent Pauli source conventions and explicit two-window accounting."""
import importlib.util
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from recursive_horizons import nsc_incoming_window_refinement as helper

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/'results/development/nsc-incoming-group32-source.json'


def script():
    spec = importlib.util.spec_from_file_location('source32', ROOT/'scripts/derive_nsc_incoming_group32_source.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def test_frozen_helper_and_actual_group32_windows():
    module = script(); assert module.digest(module.HELPER) == module.HELPER_SHA
    for _, group, kind, left, right, _ in module.CASES:
        owned = helper.authenticated_window(ROOT, group, kind, left, right)
        assert group == 32 and owned['signs'] == (1, -1)
        assert owned['measure_residual'] <= 3e-11
        assert all(p['interval'] == [left, right] for p in owned['descriptors'])


def test_actual_low_and_middle_kernels_match_independent_matrix_vertices():
    channel = json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels'][32]
    with mp.workdps(80):
        a, r = mp.sqrt(3*mp.pi/2-4), mp.sqrt(2)
        m = mp.mpf(channel['compact_mass'])
        s1 = mp.matrix([[0, 1], [1, 0]]); s2 = mp.matrix([[0, -1j], [1j, 0]]); s3 = mp.matrix([[1, 0], [0, -1]])
        for sign in (1, -1):
            ell = sign*mp.mpf(channel['angular_eigenvalue'])
            coefficients = helper._riccati_at_one(m, ell, order=24)
            for kind, E in (('low_vacuum', mp.mpf(35)), ('middle_correction', mp.mpf(64))):
                def projector(n):
                    S = sum((coefficients[j]/(2*E)**(j+1) for j in range(n)), mp.mpc(0))
                    v = mp.matrix([1, -1j*a*S]); return v*v.H/(v.H*v)[0]
                if kind == 'low_vacuum':
                    orders = helper._mp_bloch()(3*mp.pi/4, m, ell, E)
                    b = [sum(orders[j][i, 0] for j in range(5)) for i in range(3)]
                    reference = (mp.eye(2)+b[0]*s1+b[1]*s2+b[2]*s3)/2
                    delta = projector(24)-reference
                else:
                    delta = projector(24)-projector(16)
                vertices = [-m*s1+ell/r*s2-E/a*s3, -E/a*s3, E/a*mp.eye(2), ell/(2*r)*s2]
                direct = [sum((delta*V)[i, i] for i in range(2)) for V in vertices]
                actual = helper.window_kernel(E, m, ell, coefficients, kind)
                assert max(abs(x-y) for x, y in zip(actual, direct)) < mp.mpf('1e-73')
                assert actual[2] == 0


def test_record_has_two_disjoint_explicit_changes_and_nonzero_thermal():
    if not OUTPUT.exists(): pytest.skip('post-preparation source record gate')
    record = json.loads(OUTPUT.read_text())
    assert set(record['windows']) == {'group32_low', 'group32_middle'}
    assert record['scope']['radial_recurrences_rerun'] == 0
    assert record['scope']['modes_or_scattering_solved'] == 0
    assert not record['aggregate']['union_added_as_third_source']
    assert not record['failures']
    deltas = np.zeros(4)
    for row in record['windows'].values():
        np.testing.assert_allclose(np.array(row['old_source'])+row['explicit_source_delta'], row['updated_source_approximant'], rtol=0, atol=3e-22)
        assert row['source_order'] == 24
        assert row['physical_error_budget']['stationarity_tolerance'] == 3e-11
        assert all(v > 0 for v in row['thermal_scattering_stress_upper'])
        assert not row['physical_error_budget']['rigorous_source_quadrature_error_included']
        deltas += row['explicit_source_delta']
    np.testing.assert_array_equal(deltas, record['aggregate']['explicit_source_delta'])
    assert record['windows']['group32_low']['retained_archived_thermal'] is None
    assert record['windows']['group32_middle']['retained_archived_thermal'] is not None


def test_replay_never_prepares_sources_or_radial_coefficients(monkeypatch):
    if not OUTPUT.exists(): pytest.skip('post-preparation source replay gate')
    module = script()
    def forbidden(*a, **kw): raise AssertionError('scientific preparation called by replay')
    monkeypatch.setattr(module, 'prepare_window', forbidden)
    monkeypatch.setattr(helper, '_riccati_at_one', forbidden)
    monkeypatch.setattr(helper, 'window_quadrature', forbidden)
    prior = json.loads(OUTPUT.read_text())
    assert module.make_record(ROOT/prior['payload']['path']) == prior
