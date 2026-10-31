"""Independent round, cosmological and deformed Schwarzschild horizon controls."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_trapped_surfaces import (
    StarShapedSurface, _angular_basis, evaluate_surface, find_mots, round_surface,
)


NORTH=np.array([0.,0.,0.,1.])


def round_geometry(a, expansion_rate=0., r=1.):
    def provider(points):
        h=np.broadcast_to((a*r)**2*np.eye(3),(len(points),3,3)).copy()
        return h,np.zeros((len(points),3,3,3)),-expansion_rate*h
    return provider


def schwarzschild_pg(mass=1., chart_length=2., r=1.):
    """Euclidean PG slice pulled back by y=L*X_xyz/(1+X4).

    h=(L*r/(1+w))^2 delta and K=(v/R)(h-3/2 s_flat s_flat)
    with infalling PG shift +sqrt(2M/R). No surface code supplies these.
    """
    def provider(points):
        x=points[:,:3];w=points[:,3]
        sine=np.linalg.norm(x,axis=1)
        R=chart_length*sine/(1+w)
        factor=(chart_length*r/(1+w))**2
        h=factor[:,None,None]*np.eye(3)
        dh=(2*r*x/(1+w)[:,None])[:,:,None,None]*h[:,None,:,:]
        n=x/sine[:,None]
        speed=np.sqrt(2*mass/R)
        K=(factor*speed/R)[:,None,None]*(np.eye(3)-1.5*n[:,:,None]*n[:,None,:])
        return dict(h=h,dh=dh,K=K)
    return provider


def test_analytic_s2_harmonic_jets_and_normalization():
    angular=_angular_basis(3,(12,25));Y=angular['values']
    np.testing.assert_allclose(Y.T@(angular['weights'][:,None]*Y),np.eye(16),atol=3.e-15,rtol=0)
    theta,phi=angular['theta'],angular['phi']
    j=angular['labels'].index((1,1,'sin'))
    norm=-math.sqrt(3/(4*math.pi))
    np.testing.assert_allclose(Y[:,j],norm*np.sin(theta)*np.sin(phi),atol=3.e-16,rtol=0)
    np.testing.assert_allclose(angular['first'][:,0,j],norm*np.cos(theta)*np.sin(phi),atol=1.e-15,rtol=0)
    np.testing.assert_allclose(angular['first'][:,1,j],norm*np.sin(theta)*np.cos(phi),atol=3.e-16,rtol=0)
    np.testing.assert_allclose(angular['second'][:,0,1,j],norm*np.cos(theta)*np.cos(phi),atol=1.e-15,rtol=0)
    lap=(angular['second'][:,0,0,:]+np.cos(theta)[:,None]/np.sin(theta)[:,None]*angular['first'][:,0,:]
         +angular['second'][:,1,1,:]/np.sin(theta)[:,None]**2)
    eigen=np.array([ell*(ell+1) for ell,m,kind in angular['labels']])
    np.testing.assert_allclose(lap,-Y*eigen,atol=3.e-14,rtol=0)


@pytest.mark.parametrize('center',[NORTH,np.array([.2,-.5,.3,.7])])
@pytest.mark.parametrize('r',[.27,1.3])
def test_round_s3_mean_curvature_and_area_about_arbitrary_center(center,r):
    a=2.7;rho=.83
    result=evaluate_surface(round_surface(center,rho),round_geometry(a,r=r),frame_scale=r)
    np.testing.assert_allclose(result.mean_curvature,2/(a*math.tan(rho)),atol=3.e-14,rtol=0)
    np.testing.assert_allclose(result.theta_plus,-result.theta_minus,atol=0,rtol=0)
    assert result.area==pytest.approx(4*math.pi*(a*math.sin(rho))**2,rel=2.e-14)
    np.testing.assert_allclose(np.sum(result.points**2,axis=1),1.,atol=5.e-16,rtol=0)
    np.testing.assert_allclose((a*r)**2*np.sum(result.normal**2,axis=1),1.,atol=8.e-16,rtol=0)
    c=np.asarray(center)/np.linalg.norm(center)
    radial=(math.cos(rho)*result.points-c)/math.sin(rho)
    x,y,z,w=result.points.T
    frame=np.stack((np.c_[w,z,-y,-x],np.c_[-z,w,x,-y],np.c_[y,-x,w,-z]),axis=1)
    radial_components=np.einsum('nik,nk->ni',frame,radial)/r
    assert np.min(np.sum(result.normal*radial_components,axis=1))>0


def test_time_symmetric_deformed_surface_cannot_be_future_trapped():
    initial=round_surface([.3,.1,-.2,.9],1.05,degree=3)
    coeff=initial.coefficients.copy();coeff[4]=.24;coeff[8]=-.17;coeff[14]=.09
    surface=StarShapedSurface(initial.center,coeff,initial.degree)
    result=evaluate_surface(surface,round_geometry(1.4),frame_scale=1.)
    assert np.ptp(result.radius)>.08
    np.testing.assert_allclose(result.theta_plus+result.theta_minus,0.,atol=0,rtol=0)
    assert not result.future_trapped_sampled and not result.future_mots_sampled


def test_anisotropic_regular_frame_agrees_with_independent_level_set_hessian():
    # Independent level-set formula H=(h^{-1}-s s):D^2 F/|DF|.
    # This probes nonconformal frame connection terms, absent in round tests.
    r=.38;diagonal=np.array([.7,1.6,1.2]);center=np.array([.2,-.1,.3,.8])
    center/=np.linalg.norm(center);rho=1.13
    def provider(points):
        h=np.broadcast_to(np.diag(diagonal),(len(points),3,3)).copy()
        return h,np.zeros((len(points),3,3,3)),np.zeros_like(h)
    result=evaluate_surface(round_surface(center,rho),provider,frame_scale=r)
    A=np.array([[[0,0,0,1],[0,0,1,0],[0,-1,0,0],[-1,0,0,0]],
                [[0,0,-1,0],[0,0,0,1],[1,0,0,0],[0,-1,0,0]],
                [[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]]],dtype=float)
    grad=-r*np.einsum('a,iab,nb->ni',center,A,result.points)
    hessian=np.empty((len(result.points),3,3))
    for i in range(3):
        for j in range(3):
            hessian[:,i,j]=-r*r*np.einsum('a,ab,nb->n',center,A[j]@A[i],result.points)
            for k in range(3):
                if len({i,j,k})==3:
                    sign=1 if (i,j,k) in ((0,1,2),(1,2,0),(2,0,1)) else -1
                    Gamma=r*sign*(-diagonal[i]+diagonal[j]+diagonal[k])/diagonal[k]
                    hessian[:,i,j]-=Gamma*grad[:,k]
    norm=np.sqrt(np.sum(grad*grad/diagonal,axis=1))
    normal=grad/diagonal/norm[:,None]
    projector=np.diag(1/diagonal)-normal[:,:,None]*normal[:,None,:]
    expected=np.einsum('nij,nij->n',projector,hessian)/norm
    np.testing.assert_allclose(result.mean_curvature,expected,atol=8.e-15,rtol=0)


def test_closed_de_sitter_cosmological_marginal_surface_is_not_future_trapped():
    H=.37;u=.8;r=.23
    a=math.cosh(u)/H;h=H*math.tanh(u)
    rho=math.asin(1/math.cosh(u))
    geometry=round_geometry(a,h,r)
    small=evaluate_surface(round_surface(NORTH,rho),geometry,frame_scale=r)
    np.testing.assert_allclose(small.theta_minus,0.,atol=2.e-15,rtol=0)
    np.testing.assert_allclose(small.theta_plus,4*h,atol=2.e-15,rtol=0)
    assert small.area==pytest.approx(4*math.pi/H**2,rel=2.e-14)
    assert not small.future_trapped_sampled and not small.future_mots_sampled
    # The complementary plus-marginal sphere has POSITIVE ingoing expansion.
    found=find_mots(geometry,center=NORTH,initial_radius=2.1,frame_scale=r,degree=0)
    assert found.converged and not found.future_mots
    np.testing.assert_allclose(found.evaluation.theta_minus,4*h,atol=2.e-9,rtol=0)
    assert found.evaluation.radius.mean()==pytest.approx(math.pi-rho,abs=2.e-9)


@pytest.mark.parametrize('R',[1.,2.,4.])
def test_schwarzschild_pg_expansions_and_trapped_signs(R):
    mass=1.;L=2.;r=.41
    surface=round_surface(NORTH,2*math.atan(R/L))
    result=evaluate_surface(surface,schwarzschild_pg(mass,L,r),frame_scale=r)
    speed=math.sqrt(2*mass/R)
    np.testing.assert_allclose(result.theta_plus,2*(1-speed)/R,atol=4.e-14,rtol=0)
    np.testing.assert_allclose(result.theta_minus,-2*(1+speed)/R,atol=4.e-14,rtol=0)
    assert result.area==pytest.approx(4*math.pi*R**2,rel=2.e-14)
    assert result.future_trapped_sampled==(R<2*mass)
    assert result.future_mots_sampled==(R==2*mass)


def test_newton_finder_solves_a_deformed_closed_surface_not_only_round_radii():
    # The physical horizon is R=2M about NORTH. Describing it about another
    # S3 center forces a genuine nonconstant rho(theta,phi) graph.
    angle=.07;center=np.array([math.sin(angle),0.,0.,math.cos(angle)])
    result=find_mots(schwarzschild_pg(r=.61),center=center,initial_radius=1.48,
                     frame_scale=.61,degree=5,max_nfev=40,expansion_tolerance=2.e-7)
    assert result.optimizer_success and result.converged and result.future_mots
    assert np.ptp(result.evaluation.radius)>.13
    assert np.linalg.norm(result.surface.coefficients[1:])>.1
    # Independent exact Schwarzschild conditions on every denser output node.
    points=result.evaluation.points
    R=2*np.linalg.norm(points[:,:3],axis=1)/(1+points[:,3])
    np.testing.assert_allclose(R,2.,atol=2.e-8,rtol=0)
    np.testing.assert_allclose(result.evaluation.theta_minus,-2.,atol=2.e-7,rtol=0)
    assert result.evaluation.area==pytest.approx(16*math.pi,rel=2.e-9)
    summary=result.summary()
    assert summary['theta_plus_max_abs']<2.e-7
    assert summary['quadrature_order'][0]>14
    assert 'no continuum sign bound' in summary['scope']
    # Reuse every nonconstant shape coefficient on another physical chart
    # scale, rather than discarding it and repeating a coarse round search.
    continued=find_mots(schwarzschild_pg(r=.63),frame_scale=.63,initial_surface=result.surface,
                         expansion_tolerance=2.e-7,max_nfev=4)
    assert continued.converged and continued.future_mots
    assert np.ptp(continued.evaluation.radius)>.13
    assert continued.evaluation.area==pytest.approx(16*math.pi,rel=2.e-9)


def test_shape_chart_and_geometry_failures_are_explicit():
    with pytest.raises(ValueError,match='center'):round_surface([0,0,0,0],1.)
    with pytest.raises(ValueError,match='radius'):round_surface(NORTH,0.)
    with pytest.raises(ValueError,match='coefficients'):StarShapedSurface(NORTH,np.ones(3),1)
    with pytest.raises(ValueError,match='denser'):
        find_mots(round_geometry(1.),center=NORTH,initial_radius=1.,frame_scale=1.,degree=0,
                  quadrature_order=10,validation_order=10)
    def invalid(points):
        return -np.broadcast_to(np.eye(3),(len(points),3,3)),np.zeros((len(points),3,3,3)),np.zeros((len(points),3,3))
    with pytest.raises(ValueError,match='positive definite'):
        evaluate_surface(round_surface(NORTH,1.),invalid,frame_scale=1.)
    seed=round_surface(NORTH,1.,degree=2)
    for conflict in (dict(degree=3),dict(center=[1.,0,0,0]),dict(initial_radius=1.1),dict(radius_bounds=(.1,3.))):
        with pytest.raises(ValueError,match='conflict|combined'):
            find_mots(round_geometry(1.),initial_surface=seed,frame_scale=1.,**conflict)
