"""Exact solutions and nontrivial regular-frame identities for the GH kernel."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_gh_pointwise import rhs_nodes
from recursive_horizons.nsc_atlas_parent_scale_selector import VolumeParent


def data(N=1):
    g=np.repeat(np.diag([-1.,1.,1.,1.])[None],N,axis=0)
    return [g,np.zeros((N,4,4)),np.zeros((N,3,4,4)),np.zeros((N,3,3,4,4)),np.zeros((N,3,4,4)),
            np.zeros((N,2)),np.zeros((N,2)),np.zeros((N,3,2)),np.zeros((N,3,3,2)),np.zeros((N,3,2))]


def parameters():
    minimum=next(x for x in VolumeParent().vacuum_roots() if x['kind']=='minimum')
    return minimum['s'],math.sqrt(minimum['vacuum']/3)


def test_exact_de_sitter_at_multiple_proper_times():
    star,H=parameters()
    for t in (0.,.7,2.):
        jets=data(3);g,v=jets[:2]
        g[:,1:,1:]=math.cosh(H*t)**2*np.eye(3)
        v[:,1:,1:]=H*math.sinh(2*H*t)*np.eye(3)
        ag,aq,C,valid=rhs_nodes(*jets,r=H,Hsource0=3*H*math.tanh(H*t),
                              HsourceDot0=3*H**2/math.cosh(H*t)**2,sstar=star,gamma0=.2)
        expected=np.zeros((4,4));expected[1:,1:]=2*H**2*math.cosh(2*H*t)*np.eye(3)
        assert valid.all()
        np.testing.assert_allclose(ag,np.repeat(expected[None],3,axis=0),atol=3e-13,rtol=3e-13)
        np.testing.assert_allclose(C,0.,atol=5e-15)
        np.testing.assert_allclose(aq,0.,atol=2e-13)


def test_flat_vacuum_and_constant_shift_with_spacelike_coordinate_time():
    star,_=parameters();jets=data(2)
    # ADM lapse is1 despite g_tt>0 at node1. Rejecting by g_tt alone is wrong.
    jets[0][1,0,0]=3.;jets[0][1,0,1]=2.;jets[0][1,1,0]=2.
    ag,aq,C,valid=rhs_nodes(*jets,r=0,Hsource0=0,HsourceDot0=0,sstar=star,kappa=0)
    assert valid.all()
    np.testing.assert_allclose(ag,0.,atol=0)
    np.testing.assert_allclose(C,0.,atol=0)
    np.testing.assert_allclose(aq,0.,atol=2e-13)


def test_homogeneous_flat_lapse_gauge_exact_acceleration():
    star,_=parameters();jets=data()
    N,Nd,Ndd=1.3,.17,-.11
    jets[0][0,0,0]=-N*N;jets[1][0,0,0]=-2*N*Nd
    ag,aq,C,valid=rhs_nodes(*jets,r=0,Hsource0=-Nd/N,HsourceDot0=-Ndd/N+(Nd/N)**2,sstar=star,kappa=0)
    assert valid[0]
    expected=np.zeros((4,4));expected[0,0]=-2*(Nd*Nd+N*Ndd)
    np.testing.assert_allclose(ag[0],expected,atol=1e-14)
    np.testing.assert_allclose(C,0.,atol=2e-16)
    np.testing.assert_allclose(aq,0.,atol=2e-13)


def test_static_flat_metric_in_nonlinear_harmonic_spatial_coordinates():
    # x'=x+alpha(y²-z²) is harmonic in flat space. Its nonconstant metric
    # exercises both spatial second derivatives and the quadratic GH terms.
    star,_=parameters();jets=data();g,v,dg,d2g=jets[:4]
    alpha,y,z=.3,.2,-.1
    g[0,1:,1:]=[[1.,-2*alpha*y,2*alpha*z],[-2*alpha*y,1+4*alpha*alpha*y*y,-4*alpha*alpha*y*z],
                 [2*alpha*z,-4*alpha*alpha*y*z,1+4*alpha*alpha*z*z]]
    dg[0,1,1,2]=dg[0,1,2,1]=-2*alpha
    dg[0,1,2,2]=8*alpha*alpha*y
    dg[0,1,2,3]=dg[0,1,3,2]=-4*alpha*alpha*z
    dg[0,2,1,3]=dg[0,2,3,1]=2*alpha
    dg[0,2,3,3]=8*alpha*alpha*z
    dg[0,2,2,3]=dg[0,2,3,2]=-4*alpha*alpha*y
    d2g[0,1,1,2,2]=8*alpha*alpha;d2g[0,2,2,3,3]=8*alpha*alpha
    d2g[0,1,2,2,3]=d2g[0,1,2,3,2]=-4*alpha*alpha
    d2g[0,2,1,2,3]=d2g[0,2,1,3,2]=-4*alpha*alpha
    ag,_,C,valid=rhs_nodes(*jets,r=0,Hsource0=0,HsourceDot0=0,sstar=star,kappa=0)
    assert valid[0]
    np.testing.assert_allclose(C,0.,atol=2e-15)
    np.testing.assert_allclose(ag,0.,atol=2e-14)


def test_anisotropic_static_s3_jet_matches_left_invariant_ricci():
    # A pointwise curvature control, not constraint-satisfying vacuum data.
    star,_=parameters();jets=data();A=np.array([1.2,.9,1.5]);r=.4
    jets[0][0,1:,1:]=np.diag(A)
    ag,_,C,valid=rhs_nodes(*jets,r=r,Hsource0=0,HsourceDot0=0,sstar=star,kappa=0)
    expected=np.array([-4*r*r*(A[i]**2-(A[(i+1)%3]-A[(i+2)%3])**2)/(A[(i+1)%3]*A[(i+2)%3]) for i in range(3)])
    assert valid[0]
    np.testing.assert_allclose(ag[0,1:,1:],np.diag(expected),atol=2e-14)
    np.testing.assert_allclose(C,0.,atol=2e-15)


def test_canonical_matter_and_exact_scalar_jet_equation():
    star,_=parameters();jets=data();jets[5][0]=[.04,.03];jets[6][0]=[.13,-.17]
    jets[7][0]=[[.02,.01],[-.03,.04],[.01,-.02]]
    jets[8][0,0,0]=[.07,-.09];jets[8][0,1,1]=[-.02,.04]
    kappa=.7;S,chi=jets[5][0];volume=star+S/math.sqrt(84)
    model=VolumeParent();mu2=400*math.exp(-12*volume)
    W=model.potential(volume)+mu2*chi*chi/2
    Wq=np.array([(model.potential(volume,1)-6*mu2*chi*chi)/math.sqrt(84),mu2*chi])
    derivatives=np.r_[jets[6][0][None],jets[7][0]]
    expectedg=2*kappa*(derivatives@derivatives.T+W*jets[0][0])
    expectedq=-Wq+np.einsum('iif->f',jets[8][0])
    ag,aq,C,valid=rhs_nodes(*jets,r=0,Hsource0=0,HsourceDot0=0,sstar=star,kappa=kappa)
    assert valid[0]
    np.testing.assert_allclose(ag[0],expectedg,atol=2e-13)
    np.testing.assert_allclose(aq[0],expectedq,atol=2e-13)
    np.testing.assert_allclose(C,0.,atol=0)


def test_constraint_damping_sign_on_a_flat_off_constraint_jet():
    star,_=parameters();jets=data();H0,gamma=.23,.4
    ag,_,C,valid=rhs_nodes(*jets,r=0,Hsource0=H0,HsourceDot0=0,sstar=star,kappa=0,gamma0=gamma)
    assert valid[0]
    np.testing.assert_allclose(C[0],[H0,0,0,0],atol=0)
    np.testing.assert_allclose(ag[0],gamma*H0*np.eye(4),atol=2e-15)


def test_invalid_metric_and_arithmetic_nodes_are_not_repaired():
    star,_=parameters();jets=data(7)
    jets[0][1,1,1]=-1.  # spatial metric not SPD
    jets[0][2,0,0]=0.   # zero lapse
    jets[0][3,0,0]=np.nan
    jets[5][4,0]=-1e6  # potential exponent is not representable
    jets[4][5,0,0,0]=np.nan
    jets[0][6,0,1]=.1  # nonsymmetric metric is not symmetrized
    ag,aq,C,valid=rhs_nodes(*jets,r=0,Hsource0=0,HsourceDot0=0,sstar=star,kappa=0)
    np.testing.assert_array_equal(valid,[True,False,False,False,False,False,False])
    assert np.isnan(ag[1:]).all() and np.isnan(aq[1:]).all() and np.isnan(C[1:]).all()
    assert jets[0][1,1,1]==-1. and jets[0][2,0,0]==0.


def test_structural_shape_errors_are_explicit():
    star,_=parameters();jets=data();jets[7]=np.zeros((1,2,3))
    with pytest.raises(ValueError,match='shape'):
        rhs_nodes(*jets,r=0,Hsource0=0,HsourceDot0=0,sstar=star)


def test_nonlinear_inhomogeneous_turning_slice_obeys_normal_raychaudhuri():
    """A constrained point jet with nonconstant conformal geometry and shift rate.

    The observer is the initial geodesic normal. Coordinate metric acceleration
    must have the shift-rate Lie derivative removed before reading its physical
    expansion. This is a local identity test, not a new global initial solution.
    """
    from recursive_horizons.nsc_atlas_gh_constraints import adm_constraints

    star,r=parameters();kappa=1.;Psi=.83
    gradient=np.array([.02,-.03,.01]);chi=.04
    dchi=np.array([.003,.005,-.002])
    mu2=400*math.exp(-12*star)
    W=VolumeParent().potential(star)+mu2*chi*chi/2
    rho=W+np.dot(dchi,dchi)/(2*Psi**4)
    laplacian=(6*r*r/Psi**4-2*kappa*rho)*Psi**5/8
    intrinsic_hessian=np.diag([laplacian/3+.01,laplacian/3-.02,laplacian/3+.01])
    intrinsic_hessian[0,1]=intrinsic_hessian[1,0]=.004
    intrinsic_hessian[0,2]=intrinsic_hessian[2,0]=-.006
    intrinsic_hessian[1,2]=intrinsic_hessian[2,1]=.003

    def epsilon(i,j,k):return (i-j)*(j-k)*(k-i)/2
    raw_hessian=intrinsic_hessian.copy()
    for i in range(3):
        for j in range(3):
            raw_hessian[i,j]+=r*sum(epsilon(i,j,k)*gradient[k] for k in range(3))
    jets=data();g,v,dg,d2g,dv,q,p,dq,d2q,dp=jets
    g[0,1:,1:]=Psi**4*np.eye(3)
    v[0,0,1:]=v[0,1:,0]=-2*gradient/Psi
    q[0,1]=chi;dq[0,:,1]=dchi
    for i in range(3):
        dg[0,i,1:,1:]=4*Psi**3*gradient[i]*np.eye(3)
        for j in range(3):
            d2g[0,i,j,1:,1:]=(12*Psi**2*gradient[i]*gradient[j]+4*Psi**3*raw_hessian[i,j])*np.eye(3)
            dv[0,i,0,j+1]=dv[0,i,j+1,0]=-2*(raw_hessian[i,j]/Psi-gradient[i]*gradient[j]/Psi**2)
            d2q[0,i,j,1]=r*sum(epsilon(i,j,k)*dchi[k] for k in range(3))

    # Independent conformal-curvature Hamiltonian check; no ADM/GH equation
    # is used to construct the target expansion acceleration.
    curvature=6*r*r/Psi**4-8*np.trace(intrinsic_hessian)/Psi**5
    assert curvature==pytest.approx(2*kappa*rho,abs=3e-15)
    actual=adm_constraints(g,v,dg,d2g,dv,q,p,dq,r=r,kappa=kappa,potential=np.array([W]))
    np.testing.assert_allclose(actual['H'],0.,atol=8e-14)
    np.testing.assert_allclose(actual['M'],0.,atol=3e-14)
    np.testing.assert_allclose(actual['Kij'],0.,atol=0.)

    connection=np.zeros((3,3,3));pPsi=gradient/Psi
    for k in range(3):
        for i in range(3):
            for j in range(3):
                connection[k,i,j]=(r*epsilon(i,j,k)+2*((1. if k==i else 0.)*pPsi[j]
                    +(1. if k==j else 0.)*pPsi[i]-(1. if i==j else 0.)*pPsi[k]))
    beta_dot_cov=v[0,0,1:]
    lie=np.empty((3,3))
    for i in range(3):
        for j in range(3):
            lie[i,j]=(dv[0,i,0,j+1]+dv[0,j,0,i+1]
                      -sum((connection[k,i,j]+connection[k,j,i])*beta_dot_cov[k] for k in range(3)))
    for Hdot in (0.,3*r*r):
        accg,_,C,valid=rhs_nodes(*jets,r=r,Hsource0=0.,HsourceDot0=Hdot,sstar=star,kappa=kappa,gamma0=r)
        assert valid[0]
        np.testing.assert_allclose(C,0.,atol=3e-16)
        theta_dot=np.trace(accg[0,1:,1:]-lie)/(2*Psi**4)
        assert theta_dot==pytest.approx(kappa*W,abs=5e-13)
