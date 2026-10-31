"""Common-cut geometry, optical second forms and reconstruction controls."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_null_cuts import (
    evaluate_parametric_cut,evaluate_spacetime_cut,initial_null_bundle,_normalize_jets,_spacetime_surface,
)
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


def round_event_geometry(scale,rate,r):
    """Independent physical product/FLRW frame connection, Gamma[a,b,c]."""
    def provider(times,points):
        a=np.asarray(scale(times))*np.ones(len(points))
        hubble=np.asarray(rate(times))*np.ones(len(points))
        g=np.zeros((len(points),4,4));g[:,0,0]=-1
        Gamma=np.zeros((len(points),4,4,4))
        for i in range(3):
            g[:,i+1,i+1]=(a*r)**2
            Gamma[:,0,i+1,i+1]=hubble*(a*r)**2
            Gamma[:,i+1,0,i+1]=Gamma[:,i+1,i+1,0]=hubble
        for i,j,k in ((0,1,2),(1,2,0),(2,0,1)):
            Gamma[:,k+1,i+1,j+1]=r;Gamma[:,k+1,j+1,i+1]=-r
        return dict(g=g,Gamma=Gamma)
    return provider


def test_untilted_spacetime_cut_reduces_to_both_existing_slice_expansions():
    H=.31;t=.7;r=.27;a=math.cosh(H*t)/H;h=H*math.tanh(H*t)
    shape=round_surface([.2,-.1,.4,.8],1.1)
    provider=round_geometry(a,h,r);order=(12,25)
    bundle=initial_null_bundle(shape,provider,frame_scale=r,order=order)
    old=evaluate_parametric_cut(bundle['X'],provider,frame_scale=r,angular_degree=1,order=order,
                                outgoing_k=bundle['outgoing']['k'])
    events=round_event_geometry(lambda times:np.cosh(H*times)/H,lambda times:H*np.tanh(H*times),r)
    cut=evaluate_spacetime_cut(np.full(len(bundle['X']),t),bundle['X'],events,frame_scale=r,
        angular_degree=1,order=order,outgoing_k=bundle['outgoing']['k'])
    np.testing.assert_allclose(cut.theta_plus,old.theta_plus,atol=5.e-15,rtol=0)
    np.testing.assert_allclose(cut.theta_minus,old.theta_minus,atol=5.e-15,rtol=0)
    np.testing.assert_allclose(cut.normal_energy,1.,atol=5.e-16,rtol=0)
    np.testing.assert_array_equal(cut.times,np.full(len(bundle['X']),t))
    assert cut.area==pytest.approx(old.area,rel=2.e-15)


def _flat_null_plane_patch(a,b):
    """Analytic Cartesian graph u=t-z=0, v=(t+z)/2=F(x,y).

    This is a LOCAL patch test. A whole spacelike closed S2 cannot lie in
    this null hyperplane; the public S2 reconstructor is not used here.
    """
    r=.43;L=2.;xy=np.array([[0.,0.],[.2,-.3],[-.1,.4]])
    F=.4+.5*(a*xy[:,0]**2+b*xy[:,1]**2)
    gradient=xy*np.array([a,b]);Hessian=np.broadcast_to(np.diag([a,b]),(len(xy),2,2))
    y=np.c_[xy,F];denominator=L*L+np.sum(y*y,axis=1)
    X=np.c_[2*L*y/denominator[:,None],2*L*L/denominator-1]
    J=np.empty((len(y),4,3));HessX=np.empty((len(y),4,3,3))
    for p in range(3):
        for j in range(3):
            J[:,p,j]=2*L*(p==j)/denominator-4*L*y[:,p]*y[:,j]/denominator**2
            for k in range(3):
                HessX[:,p,j,k]=(-4*L*((p==j)*y[:,k]+(p==k)*y[:,j]+(j==k)*y[:,p])/denominator**2
                                  +16*L*y[:,p]*y[:,j]*y[:,k]/denominator**3)
    for j in range(3):
        J[:,3,j]=-4*L*L*y[:,j]/denominator**2
        for k in range(3):
            HessX[:,3,j,k]=-4*L*L*(j==k)/denominator**2+16*L*L*y[:,j]*y[:,k]/denominator**3
    dy=np.zeros((len(y),2,3));dy[:,0,0]=dy[:,1,1]=1.;dy[:,:,2]=gradient
    ddy=np.zeros((len(y),2,2,3));ddy[:,:,:,2]=Hessian
    dX=np.einsum('nij,nAj->nAi',J,dy)
    ddX=np.einsum('nijk,nAj,nBk->nABi',HessX,dy,dy)+np.einsum('nij,nABj->nABi',J,ddy)
    A=np.array([[[0,0,0,1],[0,0,1,0],[0,-1,0,0],[-1,0,0,0]],
                [[0,0,-1,0],[0,0,0,1],[1,0,0,0],[0,-1,0,0]],
                [[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]]],dtype=float)
    frame=np.einsum('iab,nb->nia',A,X)
    Cartesian_frame=L*r*(frame[:,:,:3].swapaxes(1,2)/(1+X[:,3,None,None])
        -X[:,:3,None]*frame[:,:,3][:,None,:]/(1+X[:,3,None,None])**2)
    k=np.c_[np.ones(len(y)),np.linalg.solve(Cartesian_frame,np.broadcast_to([0.,0.,1.],(len(y),3))[:,:,None])[:,:,0]]
    def events(times,points):
        count=len(points);g=np.zeros((count,4,4));g[:,0,0]=-1
        factor=(L*r/(1+points[:,3]))**2
        g[:,1:,1:]=factor[:,None,None]*np.eye(3)
        Gamma=np.zeros((count,4,4,4));gradient_log=r*points[:,:3]/(1+points[:,3,None])
        for i in range(3):
            for j in range(3):
                for kindex in range(3):
                    sign=(1 if (i,j,kindex) in ((0,1,2),(1,2,0),(2,0,1)) else
                         -1 if (i,j,kindex) in ((1,0,2),(2,1,0),(0,2,1)) else 0)
                    Gamma[:,kindex+1,i+1,j+1]=(r*sign+(kindex==j)*gradient_log[:,i]
                        +(kindex==i)*gradient_log[:,j]-(i==j)*gradient_log[:,kindex])
        return dict(g=g,Gamma=Gamma)
    angular=dict(theta=np.full(len(y),math.pi/2),weights=np.ones(len(y)))
    return _spacetime_surface(F,gradient,Hessian,X,dX,ddX,angular,events,r,k),gradient


@pytest.mark.parametrize('a,b',[(.4,.2),(-.3,-.7)])
def test_tilted_flat_null_plane_has_zero_outgoing_and_positive_twice_laplacian_sign(a,b):
    result,gradient=_flat_null_plane_patch(a,b)
    factor=1/np.sqrt(1+np.sum(gradient**2,axis=1))
    np.testing.assert_allclose(result['theta_plus'],0.,atol=2.e-15,rtol=0)
    np.testing.assert_allclose(result['theta_minus'],2*(a+b)*factor,atol=2.e-15,rtol=0)
    np.testing.assert_allclose(result['normal_energy'],factor,atol=5.e-16,rtol=0)
    # The fixed reciprocal null pair (kout.kin=-2) gives exactly +2 Delta F.
    np.testing.assert_allclose(result['theta_minus']/result['normal_energy'],2*(a+b),atol=2.e-15,rtol=0)
    assert result['outgoing_normality_max_relative']<2.e-15


def test_time_tilted_static_s3_cut_has_laplacian_shift_at_its_time_minimum():
    a=1.8;r=.4;rho=.95;kappa=.12;order=(12,25)
    bundle=initial_null_bundle(round_surface(NORTH,rho),round_geometry(a,r=r),frame_scale=r,order=order)
    angular=_angular_basis(1,order);direction=angular['directions'][71]
    times=.3+kappa*(1-angular['directions']@direction)
    provider=round_event_geometry(lambda t:a,lambda t:0.,r)
    cut=evaluate_spacetime_cut(times,bundle['X'],provider,frame_scale=r,angular_degree=1,order=order,
                               outgoing_k=bundle['outgoing']['k'])
    shift=2*kappa/(a*math.sin(rho))**2
    mean=2/(a*math.tan(rho))
    assert cut.theta_plus[71]==pytest.approx(shift+mean,abs=4.e-14)
    assert cut.theta_minus[71]==pytest.approx(shift-mean,abs=4.e-14)
    assert cut.normal_energy[71]==pytest.approx(1.,abs=5.e-16)
    assert cut.diagnostics['time_interpolation_max']<5.e-16
    assert cut.summary()['child_formation_proved'] is False


def test_nonspacelike_time_tilt_is_rejected():
    order=(12,25);bundle=initial_null_bundle(round_surface(NORTH,1.),round_geometry(1.),frame_scale=1.,order=order)
    angular=_angular_basis(1,order);times=4*angular['directions'][:,2]
    provider=round_event_geometry(lambda t:1.,lambda t:0.,1.)
    with pytest.raises(ValueError,match='spacelike'):
        evaluate_spacetime_cut(times,bundle['X'],provider,frame_scale=1.,angular_degree=1,order=order,
                               outgoing_k=bundle['outgoing']['k'])
