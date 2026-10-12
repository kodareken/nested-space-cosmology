"""Independent projector algebra and authenticated additive-correction replay."""
import importlib.util
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

import recursive_horizons.nsc_incoming_middle_order_correction as owner


ROOT = Path(__file__).resolve().parents[1]


def _projector(series, a):
    v = mp.matrix([1, -1j*a*series]); return v*v.H/(v.H*v)[0]


def test_stable_bloch_delta_matches_direct_rank_one_projectors():
    with mp.workdps(90):
        a = mp.sqrt(3*mp.pi/2-4); S = mp.mpc('.3', '.2'); dS = mp.mpc('2e-30', '-3e-30')
        delta = _projector(S+dS, a)-_projector(S, a)
        actual = owner.stable_delta_bloch(S, dS, a)
        expected = [2*mp.re(delta[0, 1]), -2*mp.im(delta[0, 1]), delta[0, 0]-delta[1, 1]]
        assert max(abs(x-y) for x, y in zip(actual, expected)) < mp.mpf('1e-85')
        assert abs(delta[0, 0]+delta[1, 1]) < mp.mpf('1e-85')
        assert list(owner.stable_delta_bloch(S, mp.mpc(0), a)) == [0, 0, 0]


def test_kernel_uses_existing_four_vertices_and_exact_zero_trace():
    with mp.workdps(90):
        E, mass, angular = mp.mpf(20), mp.pi/2, -mp.sqrt(5)
        a, r = mp.sqrt(3*mp.pi/2-4), mp.sqrt(2)
        coefficients = [mp.mpc(1, j+1)/(j+1) for j in range(24)]
        S16 = sum(coefficients[j]/(2*E)**(j+1) for j in range(16))
        S24 = sum(coefficients[j]/(2*E)**(j+1) for j in range(24))
        delta = _projector(S24, a)-_projector(S16, a)
        s1, s2, s3 = mp.matrix([[0, 1], [1, 0]]), mp.matrix([[0, -1j], [1j, 0]]), mp.matrix([[1, 0], [0, -1]])
        vertices = [-mass*s1+angular/r*s2-E/a*s3, -E/a*s3, E/a*mp.eye(2), angular/(2*r)*s2]
        expected = [sum((delta*V)[j, j] for j in range(2)) for V in vertices]
        actual = owner.correction_kernel(E, mass, angular, coefficients)
        assert max(abs(x-y) for x, y in zip(actual, expected)) < mp.mpf('1e-84')
        assert actual[2] == 0
        with pytest.raises(ValueError, match='order24'):
            owner.correction_kernel(E, mass, angular, coefficients[:16])
        with pytest.raises(ValueError, match='interval'):
            owner.correction_kernel(200, mass, angular, coefficients)


def test_authenticated_record_is_explicit_additive_and_physical_bound_open():
    record = json.loads((ROOT/'results/development/nsc-incoming-middle-order-correction.json').read_text())
    assert not record['failures']
    assert record['old_order'] == 16 and record['new_order'] == 24
    assert set(record['per_sign']) == {'1', '-1'}
    assert record['refined48_node_vacuum_correction'][2] == 0.
    assert record['physical_order24_defect_bound'] is None
    assert record['physical_order24_defect_status'] == 'OPEN'
    assert not record['scope']['old_source_artifact_replaced']
    assert not record['scope']['archived_thermal_insertion_changed']
    assert not record['scope']['numerical_order_difference_is_physical_error_bound']
    assert not record['scope']['quadrature_or_precision_indicator_is_rigorous_bound']
    assert np.array_equal(np.array(record['old16_middle_source'])+record['refined48_node_vacuum_correction'],
                          record['additively_corrected_middle_approximation'])


def test_replay_does_not_regenerate_coefficients_or_refined_quadrature(monkeypatch):
    def forbidden(*args, **kwargs): raise AssertionError('scientific preparation entered during replay')
    monkeypatch.setattr(owner, '_riccati_at_one', forbidden)
    monkeypatch.setattr(owner, 'positive_refined_grid', forbidden)
    spec = importlib.util.spec_from_file_location('nsc_order_correction_replay', ROOT/'scripts/derive_nsc_incoming_middle_order_correction.py')
    replay = importlib.util.module_from_spec(spec); spec.loader.exec_module(replay)
    prior = json.loads((ROOT/replay.OUTPUT).read_text())
    assert replay.make_record(ROOT/prior['payload']['path']) == prior
