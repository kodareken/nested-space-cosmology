"""Arbitrary-point geometry, physical-frame signs and strict cache pins."""
import hashlib
import json
import math
from pathlib import Path
import subprocess

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_finite_parent_data import solve_initial_data
from recursive_horizons.nsc_atlas_gh_constraints import adm_constraints
from recursive_horizons.nsc_atlas_gh_history import DenseGHHistory
from recursive_horizons.nsc_atlas_parent_gh import PAIRS,ParentGH
from recursive_horizons.nsc_atlas_trapped_parent import (
    GHSurfaceGeometry,TimeSymmetricParentGeometry,_first_geometry,compare_surface_neighborhood,
)

ROOT=Path(__file__).resolve().parents[1]


def vacuum_state(tau=.7):
    system=ParentGH(2,0);y=system.initial();C=system.unpack(y)
    a=math.cosh(system.r*tau);adot=system.r*math.sinh(system.r*tau)
    for k,(i,j) in enumerate(PAIRS):
        if i==j and i:
            C[k,0]=a*a*math.sqrt(2*math.pi**2)
            C[k+10,0]=2*a*adot*math.sqrt(2*math.pi**2)
    y[:system.dof]=system.pack(C);y[system.dof+8]=tau
    return system,y,a,adot


def points(seed=80,count=19):
    x=np.random.default_rng(seed).normal(size=(count,4));return x/np.linalg.norm(x,axis=1)[:,None]


def test_vacuum_geometry_and_matter_readouts_at_arbitrary_points():
    with threadpool_limits(limits=1):
        system,y,a,adot=vacuum_state();provider=GHSurfaceGeometry(system,y,coordinate_t=.7,chunk=3)
        X=points();out=provider(X);constraints=provider.adm_sample(X)
    np.testing.assert_allclose(out['h'],np.broadcast_to(a*a*np.eye(3),out['h'].shape),atol=2e-15)
    np.testing.assert_allclose(out['dh'],0,atol=2e-15)
    np.testing.assert_allclose(out['K'],-a*adot*np.broadcast_to(np.eye(3),out['K'].shape),atol=2e-15)
    np.testing.assert_allclose(out['K_eigenvalues'],-adot/a,atol=2e-15)
    np.testing.assert_allclose(out['future_tangential_capacity'],-2*adot/a,atol=2e-15)
    np.testing.assert_allclose(out['Pi'],0,atol=0);np.testing.assert_allclose(out['dq'],0,atol=0)
    np.testing.assert_allclose(out['W'],system.parent.vacuum,atol=2e-15)
    np.testing.assert_allclose(constraints['Kij'],out['K'],atol=3e-15)
    assert max(abs(constraints['H']))<1e-12
    assert np.max(abs(constraints['M']))<1e-12


def test_round_surface_control_has_correct_physical_radius_and_future_signs():
    from recursive_horizons.nsc_atlas_trapped_surfaces import evaluate_surface,round_surface
    with threadpool_limits(limits=1):
        system,y,a,adot=vacuum_state(.7);provider=GHSurfaceGeometry(system,y,chunk=33)
        radius=.9;surface=round_surface([0,0,0,1],radius)
        result=evaluate_surface(surface,provider,frame_scale=provider.r,quadrature_order=(8,17))
    mean=2*system.r/(a*math.tan(radius));trace=-2*adot/a
    np.testing.assert_allclose(result.theta_plus,mean-trace,atol=2e-13,rtol=2e-13)
    np.testing.assert_allclose(result.theta_minus,-mean-trace,atol=2e-13,rtol=2e-13)
    assert result.area==pytest.approx(4*math.pi*(a*math.sin(radius)/system.r)**2,rel=2e-13)
    assert not result.future_trapped_sampled


def test_first_jet_spacetime_K_agrees_with_independent_ADM_koszul():
    rng=np.random.default_rng(9);count=5
    h=np.broadcast_to(np.diag([1.2,.9,1.4]),(count,3,3)).copy()
    beta=np.broadcast_to([.12,-.07,.09],(count,3));N=.83
    g=np.zeros((count,4,4));g[:,1:,1:]=h;g[:,0,1:]=np.einsum('nij,nj->ni',h,beta)
    g[:,1:,0]=g[:,0,1:];g[:,0,0]=-N*N+np.einsum('ni,nij,nj->n',beta,h,beta)
    v=rng.normal(scale=.03,size=g.shape);v=(v+v.swapaxes(1,2))/2
    dg=np.zeros((count,3,4,4));q=np.zeros((count,2));dq=np.zeros((count,3,2));r=.31
    new=_first_geometry(g,v,dg,r)
    old=adm_constraints(g,v,dg,np.zeros((count,3,3,4,4)),np.zeros_like(dg),q,q,dq,r=r,potential=.2)
    np.testing.assert_allclose(new['K'],old['Kij'],atol=3e-16,rtol=2e-14)
    assert np.max(abs(new['K']))>1e-3
    # Physical generalized eigenvalues, independent of an arbitrary basis.
    for i in range(count):
        from scipy.linalg import eigvalsh
        np.testing.assert_allclose(new['K_eigenvalues'][i],eigvalsh(new['K'][i],h[i]),atol=2e-16)


def test_conformal_initial_geometry_conversion_and_analytic_frame_derivatives():
    with threadpool_limits(limits=1):
        data=solve_initial_data(4,epsilon=.125);provider=TimeSymmetricParentGeometry(data,chunk=4)
        X=points(11,13);out=provider(X)
    psi=data.values(X)
    np.testing.assert_allclose(out['h'],psi[:,None,None]**4*np.eye(3),atol=2e-13,rtol=2e-13)
    np.testing.assert_array_equal(out['K'],np.zeros_like(out['K']))
    unit=np.stack((np.c_[X[:,3],X[:,2],-X[:,1],-X[:,0]],
        np.c_[-X[:,2],X[:,3],X[:,0],-X[:,1]],np.c_[X[:,1],-X[:,0],X[:,3],-X[:,2]]),axis=1)
    step=1e-5
    for k in range(3):
        plus=math.cos(provider.r*step)*X+math.sin(provider.r*step)*unit[:,k]
        minus=math.cos(provider.r*step)*X-math.sin(provider.r*step)*unit[:,k]
        expected=(data.values(plus)**4-data.values(minus)**4)/(2*step)
        np.testing.assert_allclose(out['dh'][:,k,0,0],expected,atol=5e-10,rtol=2e-6)
    zero=TimeSymmetricParentGeometry(solve_initial_data(4,epsilon=0),chunk=4)
    assert np.max(abs(zero.adm_sample(X)['H']))<1e-12


def cache_fixture(tmp_path):
    system,y,_,_=vacuum_state(0)
    names=['src/recursive_horizons/nsc_atlas_parent_gh.py',
        'src/recursive_horizons/nsc_atlas_s3_frame.py','src/recursive_horizons/nsc_atlas_s3_hopf.py',
        'src/recursive_horizons/nsc_atlas_parent_circulation.py','src/recursive_horizons/nsc_atlas_parent_scale_selector.py']
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    sources={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}
    metadata=dict(degree=2,epsilon=0.,gauge='lapse_controlled',backend='hopf',source_sha256=sources)
    record=dict(schema='nsc-atlas-parent-gh-v2',degree=2,epsilon=0.,prescribed_GH_gauge='lapse_controlled',
        spatial_transform='hopf',source_sha256=sources,producer_commit=commit,
        producer_source_unchanged_during_run=True,source_bound_to_commit=True,status='incomplete',
        model_inputs=dict(H_reference=system.r,vacuum=system.parent.vacuum,scalar_mass_squared=system.parent.mass_squared,
            volume_kinetic=system.K,sstar=system.sstar),
        snapshots=[dict(coordinate_t=0.,observer_proper_duration=0.)])
    path=tmp_path/'states.npz';np.savez_compressed(path,times=[0.],states=[y],metadata=json.dumps(metadata))
    return path,record,hashlib.sha256(path.read_bytes()).hexdigest()


def test_sparse_cache_requires_array_pin_and_rejects_wrong_source_metadata(tmp_path):
    path,record,sha=cache_fixture(tmp_path)
    with pytest.raises(ValueError,match='explicitly admitted'):
        GHSurfaceGeometry.from_state_cache(path,record)
    loaded=GHSurfaceGeometry.from_state_cache(path,record,expected_sha256=sha)
    assert loaded.provenance['cache_pin_origin']=='explicit_caller_admission'
    assert loaded.proper_duration==0
    with pytest.raises(ValueError,match='SHA256'):
        GHSurfaceGeometry.from_state_cache(path,record,expected_sha256='0'*64)
    wrong=dict(record,source_sha256=dict(record['source_sha256']))
    wrong['source_sha256'][next(iter(wrong['source_sha256']))]='0'*64
    with pytest.raises(ValueError,match='Git source pin'):
        GHSurfaceGeometry.from_state_cache(path,wrong,expected_sha256=sha)
    wrong=dict(record,degree=3)
    with pytest.raises(ValueError,match='metadata'):
        GHSurfaceGeometry.from_state_cache(path,wrong,expected_sha256=sha)


def test_dense_cache_pin_and_physical_endpoint_are_enforced(tmp_path):
    path,record,_=cache_fixture(tmp_path)
    with np.load(path,allow_pickle=False) as data:
        initial=data['states'][0];metadata=json.loads(str(data['metadata']))
    metadata.update(producer_commit=record['producer_commit'],physical_coordinate_endpoint=.1,proper_endpoint=.1)
    class Polynomial:
        t_old=0.;t=.1;y_old=initial
        F=np.zeros((7,len(initial)))
    system=ParentGH(2,0)
    Polynomial.F[0,system.dof+8]=.1
    history=DenseGHHistory();history.record(0.,.1,Polynomial)
    dense_path=tmp_path/'dense.npz';history.save(dense_path,metadata)
    record['snapshots'].append(dict(coordinate_t=.1,observer_proper_duration=.1))
    record['dense_history']=dict(sha256=hashlib.sha256(dense_path.read_bytes()).hexdigest())
    loaded=GHSurfaceGeometry.from_history(dense_path,record,coordinate_t=.05)
    assert loaded.proper_duration==pytest.approx(.05,abs=1e-15)
    with pytest.raises(ValueError,match='outside'):
        GHSurfaceGeometry.from_history(dense_path,record,coordinate_t=.2)
    wrong=dict(record,snapshots=[dict(coordinate_t=.11,observer_proper_duration=.1)])
    with pytest.raises(ValueError,match='endpoint'):
        GHSurfaceGeometry.from_history(dense_path,wrong,coordinate_t=.05)


def test_neighborhood_audit_preserves_physical_tensor_units_and_no_source_mutation():
    with threadpool_limits(limits=1):
        system,y,_,_=vacuum_state(.4);saved=y.copy()
        a=GHSurfaceGeometry(system,y,chunk=2);b=GHSurfaceGeometry(system,y,chunk=7)
        result=compare_surface_neighborhood([a,b],points(31,11))
    row=result['comparisons'][0]
    assert row['metric_relative_strain_RMS']==0
    assert row['extrinsic_curvature_physical_RMS']==0
    assert row['metric_reference_covariant_gradient_RMS_over_reference_scale']==0
    np.testing.assert_array_equal(y,saved)
    with pytest.raises(ValueError):a(np.array([[1,1,0,0]]))
