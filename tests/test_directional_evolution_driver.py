"""Driver preflight only: no PDE integration or external cache dependency.

The constant-positive-psi fixture is deliberately off-shell for the
Hamiltonian constraint. It tests the source, momentum identity and ledger
without presenting it as a solved formation state.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

import evolve_nsc_atlas_valve_directional as driver
import check_nsc_atlas_valve_directional as initial_producer


@pytest.fixture
def supplied(tmp_path,monkeypatch):
    monkeypatch.setattr(driver,'ROOT',tmp_path)
    with threadpool_limits(limits=1):system=driver.DirectionalEvolution(16)
    p=system.parent;z=p.mass_squared/p.H**2
    A=math.sqrt(24/(p.kappa*(z+3)));B=p.H**1.5*math.sqrt(12)/math.pi
    C=p.kappa*A*B*math.sqrt(3/5);lam=1.
    coefficients=np.zeros(21);coefficients[0]=math.sqrt(2*math.pi**2)
    cache=tmp_path/'.build/initial.npz';cache.parent.mkdir()
    np.savez_compressed(cache,coefficients=coefficients,degrees=np.arange(0,41,2),
        A=A,B=B,York_C=C,kick_amplitude=lam,H=p.H,vacuum=p.vacuum,mass_squared=p.mass_squared,kappa=p.kappa)
    source=tmp_path/'source.py';source.write_text('# synthetic source identity\n')
    digest=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    pins={'source.py':digest(source)}
    record=dict(schema='nsc-atlas-valve-directional-initial-v1',
        provenance=dict(source_unchanged=True,source_sha256_before=pins,source_sha256_after=pins),
        raw_cache=dict(path='.build/initial.npz',sha256=digest(cache)),
        preparation=dict(A=A,B=B,York_C=C,kick_amplitude=lam))
    path=tmp_path/'initial.json';path.write_text(json.dumps(record))
    driver.load_initial(system,Path('initial.json'))
    return system,coefficients,A,B,C,lam,path,cache


def test_loader_matches_the_directional_source_and_retains_its_time_parity(supplied):
    s,coeff,A,B,C,lam,_,_=supplied
    reference=initial_producer.fields(s.initial_data,coeff,lam,s.u,A,B,C)
    actual=s.pieces(s.initial)
    for key,expected in [('q',reference['psi']**2/s.H),('b',reference['psi']**2/s.H),
                         ('ka',reference['Kr']),('kb',reference['Kt']),('rho',reference['rho']),
                         ('pr',reference['pr']),('pt',reference['pt'])]:
        np.testing.assert_allclose(actual[key],expected,atol=3e-14,rtol=0)
    np.testing.assert_allclose(s.initial[5],A*s.u,atol=0,rtol=0)
    np.testing.assert_allclose(s.initial[7],reference['Pi'],atol=0,rtol=0)
    assert not s.reflection_symmetry
    with threadpool_limits(limits=1):derivative=s.rhs(0.,s.initial.ravel()).reshape(8,s.modes)
    np.testing.assert_allclose(derivative[5],s.initial[7],atol=0,rtol=0)
    assert np.max(abs(derivative[5]))>0
    np.testing.assert_allclose(derivative[5],derivative[5,::-1],atol=2e-16,rtol=0)
    assert np.max(abs(derivative[0]))>0
    np.testing.assert_allclose(derivative[0],-derivative[0,::-1],atol=2e-16,rtol=0)
    np.testing.assert_allclose(derivative[1],-derivative[1,::-1],atol=2e-16,rtol=0)


def test_cap_rates_keep_material_transfer_and_pressure_work_separate(supplied):
    s,_,A,B,C,lam,_,_=supplied;measure=driver.Measurements(s)
    rates=measure.cap_rates(s.initial)
    # psi=1: exact hemisphere work integral int_0^1 u(1-u²)^(5/2)du=1/7.
    expected_flow=4*math.pi*lam*A*B*math.sqrt(3/5)/(6*s.H)
    expected_work=4*math.pi*lam*C*A*A/(42*s.H)
    np.testing.assert_allclose(rates,[expected_flow,expected_work,-expected_work],atol=2e-16,rtol=2e-13)
    assert rates[0]>0 and rates[1]>0 and rates[2]<0
    caps=measure.caps(s.initial)
    assert caps[0]['energy']==pytest.approx(caps[1]['energy'],rel=0,abs=2e-12)
    assert caps[0]['bulk_volume_rate']==pytest.approx(0.,abs=2e-15)
    assert rates[0]+rates[1]==pytest.approx(-(-rates[0]+rates[2]),abs=2e-16)


def test_geometry_checks_physical_momentum_but_does_not_admit_off_shell_H(supplied):
    s,_,_,_,C,lam,_,_=supplied;measure=driver.Measurements(s)
    result=measure.geometry(s.initial)
    assert result['momentum_RMS_over_H2']<1e-12
    assert result['global_max_relative_constraints']['M']<1e-12
    assert result['global_max_relative_constraints']['H']>.02
    assert not result['sampled_admitted']
    u=np.array([-.8,-.2,0.,.3,.75]);v=1-u*u
    plus,minus=s.expansions(s.initial,u)
    Kt=-lam*C*u*v/12
    np.testing.assert_allclose(plus,-2*Kt+2*s.H*u/np.sqrt(v),atol=3e-14,rtol=0)
    np.testing.assert_allclose(minus,-2*Kt-2*s.H*u/np.sqrt(v),atol=3e-14,rtol=0)


def test_loader_refuses_modified_arrays_before_using_them(supplied):
    s,_,_,_,_,_,path,cache=supplied
    with cache.open('ab') as stream:stream.write(b'changed')
    with pytest.raises(ValueError,match='initial array identity'):
        driver.load_initial(s,path)
