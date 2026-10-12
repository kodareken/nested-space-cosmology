"""Full-window observable conventions and immutable preflight/cache reuse."""
import importlib.util
import json
from pathlib import Path

import mpmath as mp
import pytest

import recursive_horizons.nsc_incoming_validated_energy_window as owner
from recursive_horizons.nsc_incoming_source_tail import _mp_bloch
from recursive_horizons.nsc_incoming_validated_vacuum_propagation import cpoint,unpoint
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision,_lo,_hi
from recursive_horizons.nsc_incoming_window_refinement import authenticated_window

ROOT=Path(__file__).resolve().parents[1]


def inputs():
    c=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels'][32]
    cfg=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    return c,cfg


def test_actual_LOWcell_is_disjoint_from_committed_high_source():
    old=authenticated_window(ROOT,32,'low_vacuum',24.,32.)
    assert old['left']==24 and old['right']==32
    assert old['signs']==(1,-1)
    assert all(p['interval']==[24.,32.] for p in old['descriptors'])


def test_raw_polynomial_moments_match_independent_exact_degree_quadrature():
    c,_=inputs()
    with _precision(80):
        v=[[mp.mpc(0),mp.mpc(0)] for _ in range(25)]
        v[0]=[mp.mpc('1.1','.2'),mp.mpc('.3','-.1')]
        v[1]=[mp.mpc('.02','-.03'),mp.mpc('-.04','.01')]
        v[2]=[mp.mpc('.003'),mp.mpc(0,'.002')]
        serial=[[cpoint(z) for z in row] for row in v]
        a,r=mp.sqrt(3*mp.pi/2-4),mp.sqrt(2);m=mp.mpf(c['compact_mass'])
        x,w=mp.gauss_quadrature(8,'legendre')
        for sign in(1,-1):
            result=owner.raw_polynomial_integral(c,sign,serial);expected=[mp.mpf(0)]*4
            ell=sign*mp.mpf(c['angular_eigenvalue'])
            for xi,wi in zip(x,w):
                E=28+4*xi;col=[v[0][i]+xi*v[1][i]+(2*xi*xi-1)*v[2][i] for i in range(2)]
                norm=abs(col[0])**2+abs(col[1])**2
                b1=2*mp.re(col[0].conjugate()*col[1]);b2=2*mp.im(col[0].conjugate()*col[1]);b3=abs(col[0])**2-abs(col[1])**2
                parallel=-E/a*b3
                values=[-m*b1+ell/r*b2+parallel,parallel,E/a*(norm-1),ell/(2*r)*b2]
                expected=[s+4*wi*t for s,t in zip(expected,values)]
            assert all(_lo(actual)<=reference<=_hi(actual) for actual,reference in zip(result,expected))
            assert abs(expected[2])>1  # Test would catch endpoint normalization.


def test_reference_is_existing_full_ad4_with_frozen_binary_labels():
    c,_=inputs()
    with _precision(80):
        a,r=mp.sqrt(3*mp.pi/2-4),mp.sqrt(2);m=mp.mpf(c['compact_mass'])
        for sign in(1,-1):
            ell=sign*mp.mpf(c['angular_eigenvalue'])
            for E in(24,28,32):
                terms=_mp_bloch()(3*mp.pi/4,m,ell,mp.mpf(E))
                b=[sum(terms[j][i,0] for j in range(5)) for i in range(3)]
                p=-E/a*b[2];expected=[-m*b[0]+ell/r*b[1]+p,p,0,ell/(2*r)*b[1]]
                actual=owner.reference_kernel(c,sign,mp.iv.mpf(E))
                assert all(_lo(v.real)<=e<=_hi(v.real) for v,e in zip(actual,expected))


def test_geometry_cache_reuses_preflight_coefficients_without_generation(monkeypatch,tmp_path):
    record=json.loads((ROOT/'results/development/nsc-incoming-energy-panel-preflight.json').read_text())
    pre=json.loads((ROOT/record['payload']['path']).read_text())['preflight']
    geometry=pre['shared_geometry']
    (tmp_path/'geometry-0000.json').write_text(json.dumps({'key':'test','index':0,'geometry':geometry}))
    def forbidden(*args,**kwargs):raise AssertionError('shared geometry was regenerated')
    monkeypatch.setattr(owner,'shared_geometry',forbidden)
    c,cfg=inputs()
    with _precision(80):
        basis,found=owner.geometry_at(c,cfg,0,unpoint(geometry['center']),tmp_path,'test')
        assert found==geometry and len(basis)==37


def test_complete_window_replays_without_field_or_reference_generation(monkeypatch):
    path=ROOT/'results/development/nsc-incoming-validated-energy-window.json'
    if not path.exists():pytest.skip('full authorized run has not finished')
    old=json.loads(path.read_text())
    assert old['energy_interval']==[24,32] and old['complete_window_budget_pass']
    assert old['source']['exact_vacuum_current']==0
    assert old['positive_thermal_stress_error_upper'][2]>0
    assert not old['scope']['raw_endpoint_polynomial_normalized']
    def forbidden(*args,**kwargs):raise AssertionError('generation during completed window replay')
    monkeypatch.setattr(owner,'energy_step',forbidden)
    monkeypatch.setattr(owner,'reference_integral',forbidden)
    spec=importlib.util.spec_from_file_location('window_replay',ROOT/'scripts/derive_nsc_incoming_validated_energy_window.py')
    replay=importlib.util.module_from_spec(spec);spec.loader.exec_module(replay)
    assert replay.make_record(old['sign_payloads'],old['observable_payload'])==old
