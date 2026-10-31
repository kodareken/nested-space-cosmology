"""Common-cut geometry, optical second forms and reconstruction controls."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_null_cuts import evaluate_parametric_cut,initial_null_bundle,_normalize_jets
from recursive_horizons.nsc_atlas_trapped_surfaces import (
    StarShapedSurface,_angular_basis,evaluate_surface,round_surface,
)


NORTH=np.array([0.,0.,0.,1.])


def round_geometry(a,hubble=0.,r=1.):
    def provider(points):
        h=np.broadcast_to((a*r)**2*np.eye(3),(len(points),3,3)).copy()
        return dict(h=h,dh=np.zeros((len(points),3,3,3)),K=-hubble*h)
    return provider


def flat_stereographic(r=.4,L=2.):
    def provider(points):
        factor=(L*r/(1+points[:,3]))**2
        h=factor[:,None,None]*np.eye(3)
        dh=2*r*points[:,:3,None,None]/(1+points[:,3,None,None,None])*h[:,None,:,:]
        return dict(h=h,dh=dh,K=np.zeros_like(h))
    return provider


def test_full_normalization_chain_agrees_with_independent_vector_derivatives():
    # R(u,v)=a+u b+v c+u v d at u=v=0; finite differences test the
    # nontrivial radial-normalization cross terms independently.
    a=np.array([.4,-.2,.7,1.1]);b=np.array([.2,.3,-.1,.4])
    c=np.array([-.3,.1,.2,.1]);d=np.array([.2,-.1,.3,-.2])
    first=np.array([[b,c]]);second=np.zeros((1,2,2,4));second[0,0,1]=second[0,1,0]=d
    X,dX,ddX,norm=_normalize_jets(a[None,:],first,second)
    def f(u,v):
        raw=a+u*b+v*c+u*v*d
        return raw/np.linalg.norm(raw)
    step=1.e-4
    np.testing.assert_allclose(dX[0,0],(f(step,0)-f(-step,0))/(2*step),atol=3.e-10,rtol=0)
    mixed=(f(step,step)-f(step,-step)-f(-step,step)+f(-step,-step))/(4*step**2)
    np.testing.assert_allclose(ddX[0,0,1],mixed,atol=8.e-9,rtol=0)
    np.testing.assert_allclose(np.einsum('ni,nai->na',X,dX),0.,atol=6.e-17,rtol=0)
    np.testing.assert_allclose(np.einsum('ni,nabi->nab',X,ddX),-np.einsum('nai,nbi->nab',dX,dX),atol=2.e-16,rtol=0)


def test_round_closed_de_sitter_cut_and_complete_optical_matrices():
    H=.31;u=.7;r=.23;a=math.cosh(u)/H;h=H*math.tanh(u)
    rho=math.asin(1/math.cosh(u));provider=round_geometry(a,h,r)
    surface=round_surface(NORTH,rho)
    bundle=initial_null_bundle(surface,provider,frame_scale=r,order=(12,25))
    cut=evaluate_parametric_cut(bundle['X'],provider,frame_scale=r,angular_degree=1,order=(12,25),
                                outgoing_k=bundle['outgoing']['k'])
    np.testing.assert_allclose(cut.theta_plus,4*h,atol=2.e-13,rtol=0)
    np.testing.assert_allclose(cut.theta_minus,0.,atol=2.e-13,rtol=0)
    assert cut.area==pytest.approx(4*math.pi/H**2,rel=2.e-14)
    np.testing.assert_allclose(bundle['ingoing']['P'],0.,atol=1.e-15,rtol=0)
    np.testing.assert_allclose(bundle['outgoing']['P'],np.broadcast_to(2*h*np.eye(2),(300,2,2)),atol=1.e-15,rtol=0)
    assert not cut.future_trapped_sampled and not cut.future_mots_sampled
    assert cut.diagnostics['embedding_interpolation_max']<2.e-15
    assert cut.diagnostics['outgoing_normality_max_relative']<2.e-14
    assert cut.diagnostics['outgoing_null_defect_max_relative']<2.e-14


def test_rotated_deformed_graph_matches_existing_exact_graph_evaluator():
    seed=round_surface([.2,-.3,.4,.8],1.05,degree=2)
    c=seed.coefficients.copy();c[2]=.13;c[4]=-.17;c[7]=.09
    surface=StarShapedSurface(seed.center,c,2)
    provider=round_geometry(1.7,.12,r=.6);order=(22,45)
    bundle=initial_null_bundle(surface,provider,frame_scale=.6,order=order)
    exact=evaluate_surface(surface,provider,frame_scale=.6,quadrature_order=order)
    np.testing.assert_allclose(np.trace(bundle['outgoing']['P'],axis1=1,axis2=2),exact.theta_plus,atol=2.e-15,rtol=0)
    np.testing.assert_allclose(np.trace(bundle['ingoing']['P'],axis1=1,axis2=2),exact.theta_minus,atol=2.e-15,rtol=0)
    cut=evaluate_parametric_cut(bundle['X'],provider,frame_scale=.6,angular_degree=12,order=order,
                                outgoing_k=bundle['outgoing']['k'])
    np.testing.assert_allclose(cut.theta_plus,exact.theta_plus,atol=2.e-7,rtol=0)
    np.testing.assert_allclose(cut.theta_minus,exact.theta_minus,atol=2.e-7,rtol=0)
    np.testing.assert_allclose(cut.normal,exact.normal,atol=2.e-8,rtol=0)
    assert cut.area==pytest.approx(exact.area,rel=2.e-9)


def test_initial_null_and_screen_normalization_with_shift_and_full_anisotropic_K():
    h0=np.diag([.8,1.2,1.5]);K0=np.array([[.2,.04,-.03],[.04,-.1,.02],[-.03,.02,.07]])
    N=1.3;beta=np.array([.1,-.07,.04])
    def provider(points):
        count=len(points)
        return dict(h=np.broadcast_to(h0,(count,3,3)),dh=np.zeros((count,3,3,3)),
                    K=np.broadcast_to(K0,(count,3,3)),lapse=np.full(count,N),shift=np.broadcast_to(beta,(count,3)))
    bundle=initial_null_bundle(round_surface([.2,.1,-.2,.9],.8),provider,frame_scale=.7,order=10)
    g=np.zeros((4,4));g[1:,1:]=h0;g[0,1:]=g[1:,0]=h0@beta;g[0,0]=-N*N+beta@h0@beta
    screen=np.concatenate((np.zeros((len(bundle['X']),2,1)),bundle['screen']),axis=2)
    np.testing.assert_allclose(np.einsum('nAi,ij,nBj->nAB',screen,g,screen),np.broadcast_to(np.eye(2),(len(screen),2,2)),atol=8.e-16,rtol=0)
    for sign in ('outgoing','ingoing'):
        k=bundle[sign]['k']
        np.testing.assert_allclose(np.einsum('ni,ij,nj->n',k,g,k),0.,atol=1.e-15,rtol=0)
        np.testing.assert_allclose(np.einsum('nAi,ij,nj->nA',screen,g,k),0.,atol=8.e-16,rtol=0)
    np.testing.assert_allclose(bundle['outgoing']['P']+bundle['ingoing']['P'],-2*bundle['K_screen'],atol=5.e-16,rtol=0)
    assert np.max(abs(bundle['K_screen'][:,0,1]))>.03


def test_flat_round_geometry_reconstruction_converges_under_nonuniform_parameterization():
    order=(20,41);angular=_angular_basis(8,order)
    directions=angular['directions']*np.array([1.15,.9,1.])
    directions/=np.linalg.norm(directions,axis=1)[:,None]
    rho=.9;points=np.c_[math.sin(rho)*directions,np.full(len(directions),math.cos(rho))]
    provider=flat_stereographic();R=2*math.tan(rho/2)
    errors=[];fits=[]
    for degree in (1,3,5,7,9):
        cut=evaluate_parametric_cut(points,provider,frame_scale=.4,angular_degree=degree,order=order)
        errors.append(float(np.max(abs(cut.theta_plus-2/R))))
        fits.append(cut.diagnostics['embedding_interpolation_max'])
    assert all(b<a for a,b in zip(errors,errors[1:]))
    assert all(b<a*.25 for a,b in zip(errors[1:],errors[2:]))
    assert all(b<a*.2 for a,b in zip(fits,fits[1:]))
    assert errors[-1]<.001


def test_known_caustic_immersion_loss_and_mixed_orientation_are_not_admitted():
    provider=round_geometry(1.);order=(11,23);angular=_angular_basis(2,order)
    n=angular['directions'];points=np.c_[n,np.zeros(len(n))]
    with pytest.raises(ValueError,match='upstream caustic'):
        evaluate_parametric_cut(points,provider,frame_scale=1.,angular_degree=2,order=order,caustic_encountered=True)
    with pytest.raises(ValueError,match='optical determinant'):
        evaluate_parametric_cut(points,provider,frame_scale=1.,angular_degree=2,order=order,ray_area_determinant=np.zeros(len(n)))
    folded=np.c_[2*n[:,0]*n[:,2],2*n[:,1]*n[:,2],2*n[:,2]**2-1,np.zeros(len(n))]
    with pytest.raises(ValueError,match='immersion|coincidence'):
        evaluate_parametric_cut(folded,provider,frame_scale=1.,angular_degree=2,order=order)
    bundle=initial_null_bundle(round_surface(NORTH,1.),provider,frame_scale=1.,order=order)
    k=bundle['outgoing']['k'].copy();k[::2,1:]*=-1
    with pytest.raises(ValueError,match='mixed outgoing orientation'):
        evaluate_parametric_cut(bundle['X'],provider,frame_scale=1.,angular_degree=2,order=order,outgoing_k=k)


def test_summary_does_not_promote_a_sampled_cut_to_a_child():
    bundle=initial_null_bundle(round_surface(NORTH,1.),round_geometry(1.),frame_scale=1.,order=8)
    result=evaluate_parametric_cut(bundle['X'],round_geometry(1.),frame_scale=1.,angular_degree=1,order=8)
    summary=result.summary()
    assert summary['immersion_min_gramdet']>0
    assert summary['global_injectivity_certified'] is False
    assert summary['continuum_expansion_signs_certified'] is False
    assert summary['child_formation_proved'] is False
    assert 'common coordinate-time cut' in summary['scope']
