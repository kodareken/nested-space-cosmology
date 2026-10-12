"""Exact endpoint identities, whitelisted branch and local recipe agreement."""
import importlib.util
import inspect
import json
from pathlib import Path

import mpmath as mp
import pytest

from recursive_horizons.nsc_incoming_tail_quadrature_bound import (
    endpoint_certificate, adapt_bloch_ast, analytic_bloch, CompactTailKernel, _gap_root, integer_power,
)
from recursive_horizons.nsc_incoming_source_tail import _mp_bloch, paired_source_function
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/'results/development/nsc-incoming-tail-quadrature-bound.json'
CHANNEL = json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels'][22]


def test_symbolic_cancellation_precedes_directed_arithmetic():
    result = endpoint_certificate()
    assert result['stable_rearrangement_residuals'] == ['0']*3
    assert result['reused_c1_c2_ad1']['vacuum_minus_ad1'] == ['0']*3
    assert result['endpoint_g0'] == 0 and not result['numerical_chopping']
    grades = result['higher_adiabatic_grade_certificate']
    assert [(r['transverse_power'],r['longitudinal_power']) for r in grades['rows']] == [(2,3),(3,4),(4,5),(5,6)]
    assert all(r['grade_residuals'] == [0,0] for r in grades['rows'])
    assert grades['ad2_ad4_have_no_remaining_poles']


def test_generated_AST_whitelist_fails_closed_and_positive_branch_agrees():
    source = inspect.getsource(_mp_bloch())
    _, counts = adapt_bloch_ast(source)
    assert counts == {'square_roots':1,'odd_gap_root_powers':[-3,-5,-7,9]}
    with pytest.raises(ValueError,match='assignment'):
        adapt_bloch_ast(source.replace('x11 = x0 + x10 + x5','x11 = x0 + x10 + 2*x5'))
    with pytest.raises(ValueError,match='whitelist'):
        adapt_bloch_ast(source.replace('x11**(-7/2)','x11**(-9/2)'))
    with pytest.raises(ValueError,match='unowned energy'):
        adapt_bloch_ast(source.replace('    return ', '    x999 = sqrt(x10)\n    return '))
    with _precision(80):
        actual = analytic_bloch()(3*mp.iv.pi/4,mp.iv.mpf('1.25'),mp.iv.mpf('2.5'),mp.iv.mpf(180))
        expected = _mp_bloch()(3*mp.pi/4,mp.mpf('1.25'),mp.mpf('2.5'),mp.mpf(180))
        for row, old in zip(actual,expected):
            for j in range(3): assert _lo(row[j,0].real) <= old[j,0] <= _hi(row[j,0].real)
        negative = _gap_root(mp.iv.mpf(-180),mp.iv.mpf(2),mp.iv.mpf(3))
        assert _hi(negative.real) < 0


def test_two_local_values_match_archived_recipe_without_reintegrating_tail():
    with _precision(80):
        kernels = [CompactTailKernel(CHANNEL,s) for s in (1,-1)]
        original = paired_source_function(CHANNEL)
        for denominator in (160,200):
            x = mp.mpf(1)/denominator; target = original(x)/(x*x)
            actual = [a+b for a,b in zip(kernels[0](mp.iv.mpf(x)),kernels[1](mp.iv.mpf(x)))]
            for j,k in ((0,0),(1,1),(3,2)):
                assert _lo(actual[j].real) <= target[k] <= _hi(actual[j].real)
            assert _lo(actual[2].real) == _hi(actual[2].real) == 0
        for kernel in kernels:
            kernel.margins(mp.iv.mpf(3)/320)


def test_localized_circle_arc_power_wrapping_is_finite_after_algebra_fix():
    with _precision(80):
        kernel = CompactTailKernel(CHANNEL,1); h=mp.iv.mpf(1)/320; radius=2*h
        theta=2*mp.iv.pi*mp.iv.mpf([mp.mpf(1)/32,mp.mpf(2)/32])
        x=mp.iv.mpc(h+radius*mp.iv.cos(theta),radius*mp.iv.sin(theta))
        values=kernel(x)
        assert all(mp.isfinite(endpoint) for v in values for endpoint in (_lo(v.real),_hi(v.real),_lo(v.imag),_hi(v.imag)))
        z=mp.iv.mpc(mp.iv.mpf([-3,-2]),mp.iv.mpf(['-.1','.1']))
        inverse=integer_power(z,-9)
        assert all(mp.isfinite(e) for e in (_lo(inverse.real),_hi(inverse.real),_lo(inverse.imag),_hi(inverse.imag)))


def test_saved_only_replay_and_nonzero_outward_thermal(monkeypatch):
    if not OUTPUT.exists(): pytest.skip('post-preparation certificate gate')
    spec = importlib.util.spec_from_file_location('tail_quad',ROOT/'scripts/derive_nsc_incoming_tail_quadrature_bound.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    def forbidden(*args,**kwargs): raise AssertionError('replay prepared numerical source')
    monkeypatch.setattr(module,'prepare_integral',forbidden)
    monkeypatch.setattr(CompactTailKernel,'__call__',forbidden)
    old = json.loads(OUTPUT.read_text()); current = module.replay(ROOT/old['payload']['path'])
    assert current == old
    assert current['thermal_tail']['exact_zero_components'] == []
    assert current['archived_source_unchanged']
    assert not current['scope']['radial_or_mode_solves']
    assert current['group'] == 22 and current['riccati_order'] == 16
