"""Analytic focusing and regular caustic continuation in the compact frame."""
import math
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_null_transport import (
    STATE_SIZE,DenseNullBundle,pack_null_state,transport_null_bundle,unpack_null_state,
)


EPS=np.zeros((3,3,3))
for i,j,k in ((0,1,2),(1,2,0),(2,0,1)):
    EPS[i,j,k]=1;EPS[j,i,k]=-1


def round_bundle(radii,r,*,focal_affine=None):
    radii=np.asarray(radii);n=len(radii)
    X=np.c_[np.sin(radii),np.zeros((n,2)),np.cos(radii)]
    k=np.broadcast_to([1.,-1.,0.,0.],(n,4)).copy()
    screen=np.broadcast_to([[0.,1.,0.],[0.,0.,1.]],(n,2,3)).copy()
    B=-r/np.tan(radii) if focal_affine is None else np.full(n,-1/focal_affine)
    return dict(X=X,k=k,screen=screen,D=np.broadcast_to(np.eye(2),(n,2,2)).copy(),
        P=B[:,None,None]*np.eye(2))


def static_sphere(r):
    g=np.diag([-1.,1.,1.,1.]);Gamma=np.zeros((4,4,4));Gamma[1:,1:,1:]=r*EPS.transpose(2,0,1)
    R=np.zeros((4,4,4,4));I=np.eye(3)
    R[1:,1:,1:,1:]=r*r*(np.einsum('ac,bd->abcd',I,I)-np.einsum('ad,bc->abcd',I,I))
    Ricci=np.diag([0.,2*r*r,2*r*r,2*r*r])
    return lambda t,X:dict(g=g,Gamma=Gamma,R_cov=R,Ricci=Ricci)


def de_sitter(r):
    def sampler(t,X):
        a=math.cosh(r*t);adot=r*math.sinh(r*t);g=np.diag([-1.,a*a,a*a,a*a]);Gamma=np.zeros((4,4,4))
        Gamma[1:,1:,1:]=r*EPS.transpose(2,0,1)
        Gamma[0,1:,1:]=a*adot*np.eye(3)
        Gamma[1:,0,1:]=Gamma[1:,1:,0]=adot/a*np.eye(3)
        R=r*r*(np.einsum('ac,bd->abcd',g,g)-np.einsum('ad,bc->abcd',g,g))
        return dict(g=g,Gamma=Gamma,R_cov=R,Ricci=3*r*r*g)
    return sampler


def test_static_closed_sphere_focus_and_reexpand_without_Riccati_failure():
    r=.7;radii=np.array([.3,.5]);focus=radii/r;end=focus[-1]+.2
    result=transport_null_bundle(static_sphere(r),round_bundle(radii,r),
        coordinate_interval=(0,end),frame_scale=r,checkpoints=[.1,focus[0],focus[0]+.02,focus[1]],
        max_step=.04,rtol=2e-11,atol=2e-13,caustic_tolerance=1e-8)
    assert result.status=='endpoint_reached' and result.caustic_encountered
    assert result.packed_states.shape==(6,2,STATE_SIZE)
    for t,state,diag in zip(result.times,result.states,result.diagnostics):
        angle=radii-r*t;D=np.sin(angle)/np.sin(radii);P=-r*np.cos(angle)/np.sin(radii)
        expectedX=np.c_[np.sin(angle),np.zeros((2,2)),np.cos(angle)]
        np.testing.assert_allclose(state['X'],expectedX,atol=2e-11,rtol=2e-11)
        np.testing.assert_allclose(state['k'],np.broadcast_to([1.,-1.,0.,0.],(2,4)),atol=2e-12)
        np.testing.assert_allclose(state['D'],D[:,None,None]*np.eye(2),atol=2e-11,rtol=2e-11)
        np.testing.assert_allclose(state['P'],P[:,None,None]*np.eye(2),atol=2e-11,rtol=2e-11)
        np.testing.assert_allclose(state['affine_lambda'],t,atol=2e-12)
        np.testing.assert_allclose(state['integrated_Rkk'],2*r*r*t,atol=2e-12)
        np.testing.assert_allclose(state['integrated_source_Rkk'],state['integrated_Rkk'],atol=2e-12)
        assert diag['matter_Ricci_role']=='geometric_Ricci_fallback_control'
        assert np.max(abs(diag['null_norm_over_normal_energy_squared']))<2e-12
        assert np.max(abs(diag['screen_gram_error']))<2e-11
        assert np.max(abs(diag['screen_k_inner_over_normal_energy']))<2e-12
        assert np.max(abs(diag['wronskian_drift']))<2e-11
        assert np.max(abs(diag['tidal_trace_minus_Ricci_kk']))<3e-11
        assert t<=diag['caustic_prefix_checked_through_t']+2e-15
        assert diag['caustic_prefix_checked_through_t']<=t+.04+2e-15
    for ray,t in enumerate(focus):
        index=int(np.argmin(abs(result.times-t)))
        assert not result.diagnostics[index]['theta_defined'][ray]
        assert math.isnan(result.diagnostics[index]['theta'][ray])
    np.testing.assert_array_equal(result.diagnostics[0]['caustic_rays_so_far'],[False,False])
    np.testing.assert_array_equal(result.diagnostics[1]['caustic_rays_so_far'],[False,False])
    np.testing.assert_array_equal(result.diagnostics[2]['caustic_rays_so_far'],[True,False])
    np.testing.assert_array_equal(result.diagnostics[3]['caustic_rays_so_far'],[True,False])
    assert result.diagnostics[3]['theta_defined'][0]
    np.testing.assert_array_equal(result.diagnostics[-1]['caustic_rays_so_far'],[True,True])
    assert result.diagnostics[0]['caustic_prefix_checked_through_t']==0.
    assert result.summary()['checkpoints'][-1]['caustic_rays_so_far']==2
    assert np.all(result.diagnostics[-1]['theta']>0)
    # Parallel transport in the noncommuting frame has this signed rotation.
    expected_screen=np.array([[0.,0.,math.cos(r*end),math.sin(r*end)],
                              [0.,0.,-math.sin(r*end),math.cos(r*end)]])
    np.testing.assert_allclose(result.final['screen'][0],expected_screen,atol=2e-11)


def test_de_sitter_has_zero_null_tidal_and_affine_linear_jacobi_through_focus():
    r=.3;R=.8;focus=math.asinh(r*R)/r;end=1.7
    initial=round_bundle([.4],r,focal_affine=R)
    result=transport_null_bundle(de_sitter(r),initial,coordinate_interval=(0,end),frame_scale=r,
        checkpoints=[focus],rtol=2e-11,atol=2e-13,max_step=.03,caustic_tolerance=1e-8)
    assert result.status=='endpoint_reached'
    for t,state,diag in zip(result.times,result.states,result.diagnostics):
        affine=math.sinh(r*t)/r;a=math.cosh(r*t);eta=math.atan(math.sinh(r*t))
        np.testing.assert_allclose(state['affine_lambda'],affine,atol=2e-11,rtol=2e-11)
        np.testing.assert_allclose(state['k'][:,0],1/a,atol=2e-11)
        np.testing.assert_allclose(state['D'],(1-affine/R)*np.eye(2)[None,:,:],atol=2e-10)
        np.testing.assert_allclose(state['P'],-np.eye(2)[None,:,:]/R,atol=2e-11)
        np.testing.assert_allclose(state['X'],[[math.sin(.4-eta),0.,0.,math.cos(.4-eta)]],atol=2e-11)
        assert np.max(abs(diag['tidal']))<3e-12
        assert np.max(abs(state['integrated_Rkk']))<3e-12
        assert np.max(abs(diag['screen_gram_error']))<3e-11
    assert not result.diagnostics[1]['theta_defined'][0]
    assert np.all(np.isfinite(result.packed_states))


def test_data_end_is_an_explicit_input_boundary_and_never_a_singularity():
    calls=[];base=static_sphere(.5)
    def sampler(t,X):calls.append(t);return base(t,X)
    initial=round_bundle([.6],.5)
    result=transport_null_bundle(sampler,initial,coordinate_interval=(0,.8),frame_scale=.5,
        checkpoints=[.2,.6],data_end=.4,max_step=.05)
    assert result.status=='data_end_reached'
    assert result.effective_interval==(0.,.4)
    assert result.requested_interval==(0.,.8)
    assert max(calls)<=.4
    assert 'no physical singularity' in result.reason
    assert result.summary()['scope'].startswith('Coordinate-time null')


def test_partial_final_snapshot_retains_its_own_caustic_prefix():
    r=.7;radius=.1;focus=radius/r;base=static_sphere(r)
    def sampler(t,X):
        if t>.35:raise ValueError('synthetic data gap')
        return base(t,X)
    result=transport_null_bundle(sampler,round_bundle([radius],r),
        coordinate_interval=(0,.8),frame_scale=r,checkpoints=[focus],max_step=.025,
        rtol=2e-11,atol=2e-13,caustic_tolerance=1e-8)
    assert result.status=='incomplete' and 'data gap' in result.reason
    assert result.times[-1]>focus
    assert not result.diagnostics[0]['caustic_rays_so_far'][0]
    assert result.diagnostics[-1]['caustic_rays_so_far'][0]
    assert result.diagnostics[-1]['caustic_prefix_checked_through_t']==result.times[-1]


def test_independent_positive_source_exposure_does_not_repair_geometric_focusing():
    r=.5;radius=.6;end=.4;base=static_sphere(r)
    # A canonical scalar-gradient Gram tensor deliberately disagrees with
    # this static geometry's Ricci contraction. The discrepancy is measured,
    # not fed back into the geometric tidal matrix or null trajectory.
    gradient=np.array([r,0.,0.,0.]);source=np.outer(gradient,gradient)
    sampler=lambda t,X:dict(base(t,X),matter_Ricci=source)
    result=transport_null_bundle(sampler,round_bundle([radius],r),
        coordinate_interval=(0,end),frame_scale=r,checkpoints=[.2],
        rtol=2e-11,atol=2e-13,max_step=.03)
    assert result.status=='endpoint_reached'
    for t,state,diag in zip(result.times,result.states,result.diagnostics):
        expected=math.sin(radius-r*t)/math.sin(radius)
        np.testing.assert_allclose(state['D'],expected*np.eye(2)[None,:,:],atol=2e-11)
        np.testing.assert_allclose(diag['Ricci_kk'],2*r*r,atol=2e-12)
        np.testing.assert_allclose(diag['source_Rkk'],r*r,atol=2e-12)
        np.testing.assert_allclose(diag['Ricci_minus_source_Rkk'],r*r,atol=2e-12)
        np.testing.assert_allclose(state['integrated_Rkk'],2*r*r*t,atol=2e-12)
        np.testing.assert_allclose(state['integrated_source_Rkk'],r*r*t,atol=2e-12)
        assert diag['matter_Ricci_role']=='provided_matter_Ricci'
    last=result.summary()['checkpoints'][-1]
    np.testing.assert_allclose(last['integrated_Ricci_minus_source_range'],[r*r*end]*2,atol=2e-12)


def test_anisotropic_geometry_keeps_null_screens_and_jacobi_wronskian():
    from recursive_horizons.nsc_atlas_null_geometry import geometry_from_jets
    r=.5;g=np.diag([-1.,1.,1.5,.8]);zero=np.zeros((4,4));dg=np.zeros((3,4,4))
    geometry=geometry_from_jets(g,zero,dg,np.zeros((3,3,4,4)),np.zeros_like(dg),zero,frame_scale=r)
    sampler=lambda t,X:dict(g=g,**geometry) # Riemann alias exercises the shared helper.
    initial=round_bundle([.4],r)
    initial['screen'][:]=[[[0.,1/math.sqrt(1.5),0.],[0.,0.,1/math.sqrt(.8)]]]
    initial['P'][:]=[[[-.3,.12],[.12,-.2]]]
    result=transport_null_bundle(sampler,initial,coordinate_interval=(0,.8),frame_scale=r,
        rtol=2e-11,atol=2e-13,max_step=.025)
    d=result.diagnostics[-1]
    assert result.status=='endpoint_reached'
    assert np.max(abs(d['null_norm_over_normal_energy_squared']))<2e-11
    assert np.max(abs(d['screen_gram_error']))<2e-11
    assert np.max(abs(d['screen_k_inner_over_normal_energy']))<2e-11
    assert np.max(abs(d['wronskian_drift']))<3e-11
    assert np.max(abs(d['tidal_trace_minus_Ricci_kk']))<2e-11
    assert np.max(abs(result.final['D']-result.final['D'].swapaxes(1,2)))>1e-6


def test_input_copies_ray_order_and_rejects_bad_initial_null_or_screen_data():
    initial=round_bundle([.4,.5],.7);original=initial['X'].copy()
    packed=pack_null_state(initial);packed[0,0]=.9
    np.testing.assert_array_equal(initial['X'],original)
    assert unpack_null_state(pack_null_state(initial))['X'].shape==(2,4)
    bad=dict(initial,k=initial['k'].copy());bad['k'][0,0]=-1
    with pytest.raises(ValueError,match='future'):pack_null_state(bad)
    bad=dict(initial,screen=initial['screen']*2)
    with pytest.raises(ValueError,match='normalization'):
        transport_null_bundle(static_sphere(.7),bad,coordinate_interval=(0,.1),frame_scale=.7)


def test_passive_dense_archive_records_every_step_without_steering_bundle(tmp_path):
    from recursive_horizons.nsc_atlas_gh_history import DenseGHHistory
    r=.5;radius=.6;initial=round_bundle([radius],r);history=DenseGHHistory()
    options=dict(coordinate_interval=(0.,.3),frame_scale=r,checkpoints=[.1,.2],
        rtol=2e-11,atol=2e-13,max_step=.025)
    plain=transport_null_bundle(static_sphere(r),initial,**options)
    recorded=transport_null_bundle(static_sphere(r),initial,dense_sink=history.record,**options)
    assert plain.status==recorded.status=='endpoint_reached'
    assert len(history.bounds)==recorded.steps==plain.steps
    assert recorded.rhs_calls>plain.rhs_calls
    np.testing.assert_array_equal(recorded.times,plain.times)
    np.testing.assert_allclose(recorded.packed_states,plain.packed_states,atol=2e-14,rtol=2e-14)
    path=tmp_path/'null-dense.npz';history.save(path,{'rays':1,'state_size':STATE_SIZE})
    restored=DenseGHHistory.load(path)
    for t,state in zip(recorded.times,recorded.packed_states):
        np.testing.assert_allclose(restored.evaluate(float(t)).reshape(1,STATE_SIZE),state,atol=2e-14,rtol=2e-14)
    # This time was never a stored checkpoint. The actual dense polynomial
    # still follows the analytic Jacobi and affine solution.
    t=.173;state=unpack_null_state(restored.evaluate(t).reshape(1,STATE_SIZE))
    expected=math.sin(radius-r*t)/math.sin(radius)
    np.testing.assert_allclose(state['D'],expected*np.eye(2)[None,:,:],atol=2e-11,rtol=2e-11)
    np.testing.assert_allclose(state['affine_lambda'],t,atol=2e-13)


def dense_reader_fixture(tmp_path,*,metadata_changes=None,focus=False):
    """Synthetic archive-format fixture; no physical trajectory claim."""
    from recursive_horizons.nsc_atlas_gh_history import DenseGHHistory
    root=Path(__file__).resolve().parents[1]
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()
    name='src/recursive_horizons/nsc_atlas_null_transport.py'
    blob=subprocess.check_output(['git','show',commit+':'+name],cwd=root)
    sources={name:hashlib.sha256(blob).hexdigest()}
    parent=dict(producer_commit='a'*40,cache_sha256='b'*64,record_sha256='c'*64)
    metadata=dict(schema='nsc-null-ray-dense-v1',producer_commit=commit,source_sha256=sources,
        parent_history=parent,rays=3,angular_rays=2,protected_probe_rays=1,angular_order=[1,2],
        coordinate_interval=[0.,1.],seed='reader_fixture')
    metadata.update(metadata_changes or {})
    base=np.zeros((3,STATE_SIZE));base[:,3]=1.;base[:,4]=1.
    base[:,16:20]=np.eye(2).ravel();velocity=np.zeros_like(base)
    velocity[:,24]=np.arange(1,4);velocity[:,25]=np.arange(3)
    velocity[:,16]=-.5 if focus else .02
    velocity[:,19]=-2. if focus else .01
    history=DenseGHHistory()
    for left,right in ((0.,.5),(.5,1.)):
        class Polynomial:pass
        dense=Polynomial();dense.t_old=left;dense.t=right
        dense.y_old=(base+left*velocity).ravel();dense.F=np.zeros((7,3*STATE_SIZE));dense.F[0]=(right-left)*velocity.ravel()
        history.record(left,right,dense)
    cache=tmp_path/'ray.npz';history.save(cache,metadata)
    record=dict(schema='nsc-atlas-valve-null-sheet-v1',producer_commit=commit,source_sha256=sources,
        source_bound_to_commit=True,source_unchanged_during_run=True,parent_history=parent,
        angular_order=[1,2],angular_rays=2,protected_probe_rays=1,seed='reader_fixture',
        transport=dict(rays=3,effective_interval=[0.,1.]),
        ray_dense_cache=dict(path=str(cache),schema='nsc-null-ray-dense-v1',intervals=2,
            sha256=hashlib.sha256(cache.read_bytes()).hexdigest()))
    path=tmp_path/'sheet.json';path.write_text(json.dumps(record),encoding='utf-8')
    return path,record,history


def test_dense_reader_events_match_actual_polynomials_per_ray_and_omit_probes(tmp_path):
    path,_,history=dense_reader_fixture(tmp_path);reader=DenseNullBundle.from_sheet_record(path)
    assert reader.angular_rays==2 and reader.rays==3 and reader.order==(1,2)
    assert reader.start==0 and reader.end==1
    times=np.array([.13,.78,.5]);out=reader.events(times)
    for field,value in out.items():
        expected=np.stack([unpack_null_state(history.evaluate(float(t)).reshape(3,STATE_SIZE))[field][ray]
            for ray,t in enumerate(times)])
        np.testing.assert_allclose(value,expected,atol=2e-15,rtol=2e-15)
    # Default angular-only requests omit the trailing protected probe.
    assert reader.events([.1,.9])['X'].shape==(2,4)
    out=reader.events([.8,.2],ray_indices=[2,0])
    np.testing.assert_allclose(out['affine_lambda'],[2.4,.2],atol=2e-15)
    for times,indices in (([-.001],None),([1.001],None),([float('nan')],None),([.2],[3]),([.2],[True])):
        with pytest.raises(ValueError):reader.events(times,ray_indices=indices)


@pytest.mark.parametrize('bad', ['cache_sha','source_sha','schema','parent','bounds','order'])
def test_dense_reader_refuses_broken_pins_and_metadata(tmp_path,bad):
    changes={}
    if bad=='schema':changes['schema']='unknown'
    if bad=='parent':changes['parent_history']=dict(producer_commit='d'*40,cache_sha256='b'*64,record_sha256='c'*64)
    if bad=='bounds':changes['coordinate_interval']=[0.,1.1]
    if bad=='order':changes['angular_order']=[1,3]
    path,record,_=dense_reader_fixture(tmp_path,metadata_changes=changes)
    if bad=='cache_sha':record['ray_dense_cache']['sha256']='0'*64
    if bad=='source_sha':record['source_sha256'][next(iter(record['source_sha256']))]='0'*64
    path.write_text(json.dumps(record),encoding='utf-8')
    with pytest.raises(ValueError):DenseNullBundle.from_sheet_record(path)


def test_stored_polynomial_noncaustic_bound_is_cautious_and_padded(tmp_path):
    path,_,_=dense_reader_fixture(tmp_path);reader=DenseNullBundle.from_sheet_record(path)
    bound=reader.stored_polynomial_noncaustic_bound()
    assert bound['all_sampled_rays_certified']
    assert .6<bound['min_lower_bound']<1/math.sqrt(2)
    assert bound['min_roundoff_padding']>0 and bound['checked_rays']==3
    assert reader.stored_polynomial_noncaustic_bound(ray_indices=[0,1])['checked_rays']==2
    directory=tmp_path/'focus';directory.mkdir()
    path,_,_=dense_reader_fixture(directory,focus=True)
    bound=DenseNullBundle.from_sheet_record(path).stored_polynomial_noncaustic_bound()
    assert not bound['all_sampled_rays_certified'] and bound['min_lower_bound']<0
    assert 'inconclusive' in bound['interpretation']
