"""Energy-edge retention, initializer interpolation and scoped step replay."""
import importlib.util
import json
from pathlib import Path

import mpmath as mp
import pytest

import recursive_horizons.nsc_incoming_energy_panel as owner
from recursive_horizons.nsc_incoming_matched_horizon_initializer import _geometry,_coefficients
from recursive_horizons.nsc_incoming_validated_vacuum_propagation import point,_ivc
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision,_lo,_hi

ROOT=Path(__file__).resolve().parents[1]


def inputs():
    c=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels'][32]
    cfg=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    return c,cfg


def evaluate(coefficients,x):
    basis=[mp.mpf(1),mp.mpf(x)]
    for n in range(2,len(coefficients)):basis.append(2*x*basis[-1]-basis[-2])
    return [sum((basis[k]*row[i] for k,row in enumerate(coefficients)),mp.mpc(0)) for i in range(2)]


def test_phase_and_energy_product_preserve_full_degree25_residual():
    identity=owner.panel_identities()
    assert identity['scalar_commutator']==['0']*4
    assert identity['Chebyshev_multiplication_residuals']==['0']*25
    assert identity['linear_vertex_kernel_degree']==49
    with _precision(80):
        coefficients=[[mp.mpc(j+1,2-j)/100,mp.mpc(3-j,j+2)/100] for j in range(25)]
        g=(mp.iv.mpc('.3'),mp.iv.mpc('.4'),mp.iv.mpc('.5'))
        product=owner._apply_energy(g,coefficients,1,mp.iv.mpf(2),mp.iv.mpf(3))
        expected_edge=4*g[0]*coefficients[24][1]
        assert _hi(abs(product[25][1]-expected_edge))<mp.mpf('1e-75')
        for x in(map(mp.mpf,('-.7','0','.6','1'))):
            v=evaluate(coefficients,x)
            plus=mp.mpc('-.8','1.5');minus=plus.conjugate()
            expected=[minus*v[1],plus*v[0]+2*(28+4*x)*mp.mpf('.3')*v[1]]
            actual=evaluate(product,x)
            for a,b in zip(actual,expected):
                assert _lo(a.real)<=b.real<=_hi(a.real)
                assert _lo(a.imag)<=b.imag<=_hi(a.imag)


def test_initial_DCT_midpoint_is_inside_componentwise_interpolation_certificate():
    c,cfg=inputs()
    with _precision(80):
        g=_geometry(c,cfg);g['mass']=mp.iv.mpf(c['compact_mass']);g['angular']=mp.iv.mpf(c['angular_eigenvalue'])
        for sign in(1,-1):
            coefficients,record=owner.initial_interpolant(c,cfg,sign)
            assert record['normalization_product_upper']<1e-9
            assert record['column_interpolation_error_upper']<3e-27
            assert record['DCT_coefficient_radius_l1_upper']<1e-35
            for x in(map(mp.mpf,('-1','-.5','0','.5','1'))):
                row=_coefficients(g,mp.iv.mpf(28+4*x),sign)
                w=mp.iv.sqrt(g['delta0'])*(row['w0']+g['delta0']*row['w1'])
                d=mp.iv.sqrt(1+w.real*w.real+w.imag*w.imag)
                exact=[1/d,w/d];actual=evaluate(coefficients,x)
                error=mp.iv.sqrt(sum((abs(v-_ivc(a))**2 for v,a in zip(exact,actual)),mp.iv.mpf(0)))
                assert _hi(error)<record['total_initial_column_error_upper']


def test_replay_charges_degree25_without_hiding_it_in_projected_equation():
    with _precision(80):
        norms=[[0.]*26 for _ in range(73)];norms[0][25]=1.
        row={'step':point(mp.mpf(1)/32),'residual_time_energy_norm_upper':norms,
             'uniform_generator_norm_upper':0.,'column_polynomial_Cheb_l1_upper':1.,
             'centering_coefficient_l1_error_upper':0.}
        result=owner.replay_step(row)
        assert result['energy_degree25_residual_upper']>=1/32
        assert not result['local_gate_pass']
        assert result['projected_low_order_arithmetic_upper']==0.


def test_preflight_receipt_replays_without_another_step(monkeypatch):
    path=ROOT/'results/development/nsc-incoming-energy-panel-preflight.json'
    if not path.exists():pytest.skip('preflight follows primitive tests')
    old=json.loads(path.read_text())
    assert old['scope']['full_panel_run']==False
    assert old['scope']['physical_source_integral']==False
    def forbidden(*args,**kwargs):raise AssertionError('preflight step rerun during replay')
    monkeypatch.setattr(owner,'energy_step',forbidden)
    monkeypatch.setattr(owner,'initial_interpolant',forbidden)
    spec=importlib.util.spec_from_file_location('energy_preflight_replay',ROOT/'scripts/derive_nsc_incoming_energy_panel_preflight.py')
    replay=importlib.util.module_from_spec(spec);spec.loader.exec_module(replay)
    assert replay.make_record(ROOT/old['payload']['path'])==old
