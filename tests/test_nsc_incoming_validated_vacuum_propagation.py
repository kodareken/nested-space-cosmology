"""Continuous error proofs checked against independent exact evolutions."""
import importlib.util
import json
from pathlib import Path

import mpmath as mp
import pytest

import recursive_horizons.nsc_incoming_validated_vacuum_propagation as owner
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision,_lo,_hi

ROOT=Path(__file__).resolve().parents[1]


def inputs():
    c=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels'][32]
    cfg=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    initial=json.loads((ROOT/'results/development/nsc-incoming-matched-horizon-initializer.json').read_text())
    return c,cfg,initial


def test_directed_Fourier_coefficients_have_correct_sign_and_alias_enclosure():
    class Polynomial:
        def values(self,z):return(1+2*z+3*z*z,mp.iv.mpc(2),mp.iv.mpc(2))
    with _precision(80):
        values=owner.generator_coefficients(Polynomial(),mp.mpf(0),{'generator_norm_upper':20.})
        for j,row in enumerate(values):
            expected=((1,2,3)[j] if j<3 else 0,2 if j==0 else 0,2 if j==0 else 0)
            for value,target in zip(row,expected):
                assert _lo(value.real)<=target<=_hi(value.real)
                assert _lo(value.imag)<=0<=_hi(value.imag)


def test_continuous_residual_encloses_exact_constant_H_step_error():
    class Constant:
        def values(self,z):return(mp.iv.mpc(7),mp.iv.mpc('.3','.2'),mp.iv.mpc('.3','-.2'))
        def disk(self,center):return{'B_real_lower':1.,'sin_real_lower':1.,'generator_norm_upper':10.}
    with _precision(80):
        h=mp.mpf(1)/32
        end,row=owner.taylor_step(Constant(),mp.mpf(0),h,[mp.mpc(1),mp.mpc(0)])
        H=mp.matrix([[-7,mp.mpc('.3','-.2')],[mp.mpc('.3','.2'),7]])
        energy=mp.sqrt(mp.mpf('49.13'))
        exact=(mp.cos(energy*h)*mp.eye(2)-1j*mp.sin(energy*h)/energy*H)*mp.matrix([1,0])
        distance=mp.sqrt(sum(abs(a-b)**2 for a,b in zip(end,exact)))
        assert distance<=row['local_column_error_upper']
        assert row['local_column_error_upper']<1e-22
        assert owner._upper_float(owner._step_error(row)[0])==row['local_column_error_upper']


def test_stable_geometry_matches_exact_profile_and_proves_disk_margins():
    c,cfg,_=inputs()
    with _precision(80):
        g=owner.ExactProfileGenerator(c,cfg,1)
        W=lambda q:3*(mp.pi-q)+mp.mpf('1.5')*mp.sin(2*q)-mp.sin(q)**2
        qh=mp.findroot(W,mp.mpf('2.657'))
        m,ell=mp.mpf(c['compact_mass']),mp.mpf(c['angular_eigenvalue'])
        for y in(mp.mpf(-23),mp.mpf(-10),mp.log(qh-3*mp.pi/4)):
            disk=g.disk(y)
            assert min(disk['B_real_lower'],disk['sin_real_lower'])>0
            value=g.values(mp.iv.mpc(y))
            delta=mp.exp(y);B=W(qh-delta)/delta;r=1/mp.sin(qh-delta)
            expected=(24/B,mp.sqrt(delta/B)*(-m*r+1j*ell),mp.sqrt(delta/B)*(-m*r-1j*ell))
            for v,e in zip(value,expected):
                assert _lo(v.real)<=e.real<=_hi(v.real)
                assert _lo(v.imag)<=e.imag<=_hi(v.imag)


def test_full_receipt_and_completed_checkpoint_resume_without_field_steps(monkeypatch,tmp_path):
    path=ROOT/'results/development/nsc-incoming-validated-vacuum-propagation.json'
    if not path.exists():pytest.skip('fixed field run follows primitive tests')
    record=json.loads(path.read_text())
    assert record['energy']==24 and not record['scope']['source_quadrature_or_window_integral_claimed']
    assert record['numerical_gate_pass']
    artifact=json.loads((ROOT/record['sign_payloads'][0]['path']).read_text())
    p=artifact['propagation'];checkpoint=tmp_path/'checkpoint.json'
    checkpoint.write_text(json.dumps({'binding':p['binding'],'steps':p['steps']}))
    def forbidden(*args,**kwargs):raise AssertionError('new field step in completed replay')
    monkeypatch.setattr(owner,'taylor_step',forbidden)
    c,cfg,initial=inputs()
    resumed=owner.propagate_sign(c,cfg,1,initial,checkpoint=checkpoint,checkpoint_key=p['binding']['key'])
    assert resumed==p
    spec=importlib.util.spec_from_file_location('validated_replay',ROOT/'scripts/derive_nsc_incoming_validated_vacuum_propagation.py')
    replay=importlib.util.module_from_spec(spec);spec.loader.exec_module(replay)
    assert replay.make_record(record['sign_payloads'])==record
