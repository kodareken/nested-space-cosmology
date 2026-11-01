"""Adapter/gauge/known-limit controls; no nonlinear formation claim."""
import importlib.util
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_parent_gh import ParentGH
from recursive_horizons.nsc_atlas_gh_pointwise import rhs_nodes

spec=importlib.util.spec_from_file_location('counterrotation_run_driver',Path(__file__).resolve().parents[1]/'scripts/evolve_nsc_atlas_valve_counterrotation.py')
driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)


def test_explicit_projection_never_invokes_old_mixed_source_initializer(monkeypatch):
    # Synthetic constant metric tests packing only, not Einstein admissibility.
    system=ParentGH(4,1.,kick_amplitude=.3,kick_source=2)
    data=SimpleNamespace(parent=system.parent,kick_amplitude=.3,degree=0,A=.1,provenance={},
        values=lambda x:np.ones(len(x)),tensor=lambda x:np.zeros((len(x),3,3)),
        f=lambda x:system.r**1.5*math.sqrt(12)/math.pi*x[:,1]*x[:,3])
    def forbidden(*a,**kw):raise AssertionError('old initial preparation was called')
    monkeypatch.setattr(ParentGH,'initial',forbidden)
    with threadpool_limits(limits=1):
        state,meta=driver.project_initial_state(system,data);J=system.jets(state)
    x=system.basis.points
    np.testing.assert_allclose(J.q[:system.N,1],data.A*x[:,0],atol=3e-16,rtol=0)
    np.testing.assert_allclose(J.p[:system.N,1],.3*data.f(x),atol=3e-16,rtol=0)
    np.testing.assert_allclose(J.q[:,0],0,atol=0,rtol=0)
    np.testing.assert_allclose(J.p[:,0],0,atol=0,rtol=0)
    assert meta['disallowed_character_leakage']<1e-14
    assert meta['observer_role'].startswith('Initially normal inflow-lobe')


def test_first_jet_gauge_constraint_matches_full_rhs_with_shift_and_feedback():
    system=ParentGH(4,1.);rng=np.random.default_rng(42);N=3
    h=np.broadcast_to(np.diag([1.1,1.2,1.3]),(N,3,3)).copy()
    lapse=np.array([.8,1.1,1.3]);beta=rng.normal(size=(N,3))*.1
    g=np.zeros((N,4,4));g[:,1:,1:]=h;g[:,0,1:]=np.einsum('nij,nj->ni',h,beta)
    g[:,1:,0]=g[:,0,1:];g[:,0,0]=-lapse*lapse+np.einsum('ni,nij,nj->n',beta,h,beta)
    symmetric=lambda shape:(lambda z:(z+z.swapaxes(-1,-2))/2)(rng.normal(size=shape)*.03)
    v=symmetric((N,4,4));dg=symmetric((N,3,4,4));dv=np.zeros_like(dg)
    d2g=np.zeros((N,3,3,4,4));q=np.zeros((N,2));p=rng.normal(size=(N,2))*.02
    dq=rng.normal(size=(N,3,2))*.01;d2q=np.zeros((N,3,3,2));dp=np.zeros_like(dq)
    t=.4;H,Hdot=system.gauge_source(t)
    _,_,expected,valid=rhs_nodes(g,v,dg,d2g,dv,q,p,dq,d2q,dp,r=system.r,Hsource0=H,
        HsourceDot0=Hdot,sstar=system.sstar,kappa=system.parent.kappa,gamma0=system.r,lapse_feedback=system.r)
    assert np.all(valid)
    actual=driver.gauge_covector(system,t,g,v,dg,lapse)
    np.testing.assert_allclose(actual,expected,atol=2e-15,rtol=2e-14)


def test_report_recovers_closed_de_sitter_energy_pressure_and_constraints():
    system=ParentGH(4,0.)
    with threadpool_limits(limits=1):
        state=system.initial();reading=driver.diagnostic(system,0.,state)
    volume=2*math.pi**2/system.r**3
    assert abs(reading['proper_slice_volume']/volume-1)<2e-14
    assert abs(reading['total_slice_normal_matter_energy']/(system.parent.vacuum*volume)-1)<2e-14
    assert abs(reading['slice_w']+1)<2e-14
    assert reading['sampled_admitted']
    assert max(reading['global_max_relative_constraints'].values())<2e-14
    assert max(reading['neck_neighborhood_max_relative_constraints'].values())<2e-14
    assert not reading['regular_child_proved']
