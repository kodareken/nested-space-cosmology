"""Matched Frobenius coefficient, exact-profile residual and ratio representation."""
import importlib.util
import json
from pathlib import Path

import mpmath as mp
import numpy as np

import recursive_horizons.nsc_incoming_matched_horizon_initializer as owner
from recursive_horizons.nsc_incoming_source_quadrature_bound import unpack
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi


ROOT = Path(__file__).resolve().parents[1]


def inputs():
    c = json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels'][32]
    cfg = json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    return c, cfg


def test_exact_coefficients_cancel_and_affine_limit_is_same_vacuum():
    identity = owner.matched_ratio_identity()
    assert identity['coefficient_residuals'] == ['0', '0']
    assert identity['residual_order'] == 'delta^(5/2)'
    assert identity['vacuum_limit'] == 'P(delta)->diag(1,0)'


def test_chart_Riccati_defect_has_same_projector_norm():
    # The imported scalar/projector identity in this logarithmic chart.
    w, dw = .3+.2j, .07-.04j
    K, D = 1.2-.8j, 2.3
    H = np.array([[-D, K.conjugate()], [K, D]])
    v, dv = np.array([1., w]), np.array([0., dw])
    norm = 1+abs(w)**2
    P = np.outer(v, v.conj())/norm
    dP = (np.outer(dv, v.conj())+np.outer(v, dv.conj()))/norm
    dP -= P*2*np.real(w.conjugate()*dw)/norm
    defect = dP+1j*(H@P-P@H)
    riccati = dw+1j*K+2j*D*w-1j*K.conjugate()*w*w
    assert abs(np.linalg.norm(defect,2)-abs(riccati)/norm) < 3e-15


def test_exact_profile_residual_is_inside_uniform_directed_majorant():
    c, cfg = inputs(); result = owner.matched_initializer_bound(c,cfg)
    assert result['mathematical_endpoint_budget_pass']
    assert result['frozen_source_kappa'] == cfg['surface_gravity']
    assert result['total_lapse_error_upper']['float_enclosure'][1] < 1e-18
    # Independent high-precision evaluation of the exact elementary profile,
    # not the old small-delta fourth-order evaluator and not an ODE solve.
    with mp.workdps(100):
        W = lambda q: 3*(mp.pi-q)+mp.mpf('1.5')*mp.sin(2*q)-mp.sin(q)**2
        qh = mp.findroot(W, mp.mpf('2.657'))
        kg = (3-3*mp.cos(2*qh)+mp.sin(2*qh))/2
        rh = 1/mp.sin(qh); r1 = mp.cos(qh)/mp.sin(qh)**2
        b0 = 2*kg; b1 = (-6*mp.sin(2*qh)-2*mp.cos(2*qh))/2
        m = mp.mpf(c['compact_mass'])
        for row in result['per_sign']:
            ell = row['sign']*mp.mpf(c['angular_eigenvalue'])
            N0 = -m*rh+1j*ell
            K0 = N0/mp.sqrt(b0)
            K1 = -m*r1/mp.sqrt(b0)-N0*b1/(2*b0**mp.mpf('1.5'))
            for E in (24,28,32):
                D0, D1 = E/b0, -E*b1/b0**2
                w0 = -1j*K0/(mp.mpf('.5')+2j*D0)
                w1 = (-1j*K1-2j*D1*w0+1j*K0.conjugate()*w0*w0)/(mp.mpf('1.5')+2j*D0)
                for d in (mp.mpf('5e-11'),mp.mpf('1e-10')):
                    q = qh-d; B = W(q)/d
                    K = (-m/mp.sin(q)+1j*ell)/mp.sqrt(B)
                    D = E/B
                    w = mp.sqrt(d)*(w0+w1*d)
                    wd = mp.sqrt(d)*(w0/2+mp.mpf('1.5')*w1*d)
                    residual = wd+1j*mp.sqrt(d)*K+2j*D*w-1j*mp.sqrt(d)*K.conjugate()*w*w
                    assert abs(residual)/d**mp.mpf('2.5') < row['residual_coefficient_upper']['float_enclosure'][1]


def test_point_binary_ratio_is_enclosed_without_changing_old_columns():
    c,cfg=inputs()
    for sign in(1,-1):
        for E in(24,28,32):
            point=owner.represent_initializer(c,cfg,E,sign)
            with _precision(80):
                re=unpack(point['ratio_interval']['real']['binary_interval'])
                im=unpack(point['ratio_interval']['imag']['binary_interval'])
                z=mp.iv.mpc(*point['binary_ratio'])
                radius=unpack(point['projector_representation_error_upper']['binary_interval'])
                assert _hi(abs(mp.iv.mpc(re,im)-z)) <= _hi(radius)
                assert _hi(radius) < mp.mpf('2e-22')


def test_receipt_replay_does_not_prepare_initializers(monkeypatch):
    path=ROOT/'results/development/nsc-incoming-matched-horizon-initializer.json'
    if not path.exists():
        import pytest
        pytest.skip('receipt follows primitive tests')
    old=json.loads(path.read_text())
    def forbidden(*args,**kwargs):raise AssertionError('initializer generation during replay')
    monkeypatch.setattr(owner,'matched_initializer_bound',forbidden)
    monkeypatch.setattr(owner,'represent_initializer',forbidden)
    spec=importlib.util.spec_from_file_location('matched_replay',ROOT/'scripts/derive_nsc_incoming_matched_horizon_initializer.py')
    replay=importlib.util.module_from_spec(spec);spec.loader.exec_module(replay)
    assert replay.make_record(ROOT/old['payload']['path'])==old
