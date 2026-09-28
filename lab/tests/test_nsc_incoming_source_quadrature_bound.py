"""Exact-rule certification, analytic continuation and honest replay scope."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import mpmath as mp
import pytest

import recursive_horizons.nsc_incoming_source_quadrature_bound as owner
from recursive_horizons.nsc_incoming_source_tail import _riccati_at_one
from recursive_horizons.nsc_incoming_window_refinement import window_kernel
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi, _range


ROOT = Path(__file__).resolve().parents[1]


def channel():
    return json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels'][6]


def test_exact_gauss_root_brackets_and_polynomial_moments():
    rule = owner.gauss_rule(8, 60)
    with _precision(60):
        pairs = owner.certify_gauss_brackets(rule['roots'], 8)
        for degree in range(16):
            value = sum((w*x**degree for x, w in pairs), mp.iv.mpf(0))
            expected = mp.mpf(0) if degree % 2 else mp.mpf(2)/(degree+1)
            assert _lo(value) <= expected <= _hi(value)
        malformed = deepcopy(rule['roots']); malformed[1] = malformed[0]
        with pytest.raises(ValueError, match='disjoint'):
            owner.certify_gauss_brackets(malformed, 8)
        malformed = deepcopy(rule['roots']); malformed[0] = owner.pack(_range(mp.mpf('-.9'), mp.mpf('-.89')))
        with pytest.raises(ArithmeticError, match='endpoint signs'):
            owner.certify_gauss_brackets(malformed, 8)


def test_directed_real_kernel_encloses_existing_source_convention():
    c = channel()
    with _precision(80):
        for sign in (1, -1):
            kernel = owner.IntervalWindowKernel(c, sign)
            mass, ell = mp.mpf(c['compact_mass']), sign*mp.mpf(c['angular_eigenvalue'])
            coefficients = _riccati_at_one(mass, ell, order=24)
            for kind in owner.KINDS:
                for E in (36, 44, 156):
                    actual = kernel(mp.iv.mpf(E), kind)
                    expected = window_kernel(E, mass, ell, coefficients, kind)
                    for a, b in zip(actual, expected):
                        assert _lo(a.real) <= b <= _hi(a.real)
                        assert _lo(a.imag) <= 0 <= _hi(a.imag)
                    assert _lo(actual[2].real) == _hi(actual[2].real) == 0


def test_complex_extension_is_coefficient_conjugate_not_energy_conjugate():
    with _precision(70):
        E = mp.iv.mpc(36, 2)
        cs = [mp.iv.mpc(1, 2), mp.iv.mpc(3, -1)]
        sbar = owner._series(E, [owner._conj_coeff(c) for c in cs])
        conjugated_value = owner._conj_coeff(owner._series(E, cs))
        assert _hi(sbar.imag) < _lo(conjugated_value.imag) or _hi(conjugated_value.imag) < _lo(sbar.imag)


def test_disk_bound_proves_domain_and_rejects_bad_disk():
    with _precision(70):
        kernel = owner.IntervalWindowKernel(channel())
        certificate = owner.disk_bound(kernel, 36, 'low_vacuum')
        assert _lo(owner.unpack(certificate['energy_squared_real_lower'])) > 0
        assert all(_hi(owner.unpack(v)) < 1 for v in certificate['normalization_product_upper'].values())
        assert max(_hi(owner.unpack(v)) for v in certificate['quadrature_error_upper']) < mp.mpf('1e-26')
        with pytest.raises(ValueError, match='disk'):
            owner.disk_bound(kernel, 12, 'low_vacuum')


def test_binary_endpoint_serialization_is_lossless():
    with _precision(60):
        value = mp.iv.pi/7-mp.iv.sqrt(2)
        assert owner.pack(owner.unpack(owner.pack(value))) == owner.pack(value)


def test_pilot_record_scope_and_replay_without_source_generation(monkeypatch):
    record_path = ROOT/'results/development/nsc-incoming-source-quadrature-bound.json'
    if not record_path.exists():
        pytest.skip('pilot record is generated after kernel checks')
    record = json.loads(record_path.read_text())
    assert record['scope']['radial_recurrences_rerun'] == 0
    assert not record['scope']['physical_projector_thermal_errors_included']
    assert 'correction only' in record['scope']['middle_certificate_covers']
    assert set(record['windows']) == {'group6_low', 'group6_middle_correction'}
    assert all(r['lapse_error_upper'] < 3e-11 for r in record['windows'].values())
    def forbidden(*args, **kwargs):
        raise AssertionError('source or coefficient generation during replay')
    monkeypatch.setattr(owner, 'IntervalWindowKernel', forbidden)
    monkeypatch.setattr(owner, 'incoming_riccati_interval_coefficients', forbidden)
    spec = importlib.util.spec_from_file_location('quadrature_replay', ROOT/'scripts/derive_nsc_incoming_source_quadrature_bound.py')
    replay = importlib.util.module_from_spec(spec); spec.loader.exec_module(replay)
    assert replay.make_record(ROOT/record['payload']['path']) == record
