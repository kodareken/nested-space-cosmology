"""Direct physical Ricci/GH identity at general noncoordinate-frame jets.

The oracle uses the Koszul connection of the physical metric, differentiates
it with raw ordered jets, and contracts the Riemann tensor. It does not use
the kernel's reference covariant-DD construction or its Eq. (24) algebra.
The point jets are manufactured local data, not constraint-satisfying initial
data or a global spacetime solution.
"""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_gh_pointwise import rhs_nodes
from recursive_horizons.nsc_atlas_parent_scale_selector import VolumeParent


def bracket(r):
    C=np.zeros((4,4,4))
    for i,j,k in ((1,2,3),(2,3,1),(3,1,2)):
        C[k,i,j]=2*r;C[k,j,i]=-2*r
    return C


def physical_geometry(g,first,second,C):
    """Gamma[a,b,c] defines nabla_(e_b)e_c=Gamma[a,b,c]e_a."""
    inverse=np.linalg.inv(g)
    dinverse=np.array([-inverse@first[m]@inverse for m in range(4)])
    lower=np.empty((4,4,4));dlower=np.empty((4,4,4,4))
    for d in range(4):
        for b in range(4):
            for c in range(4):
                value=first[b,c,d]+first[c,b,d]-first[d,b,c]
                for e in range(4):
                    value+=-g[b,e]*C[e,c,d]+g[c,e]*C[e,d,b]+g[d,e]*C[e,b,c]
                lower[d,b,c]=value/2
                for m in range(4):
                    value=second[m,b,c,d]+second[m,c,b,d]-second[m,d,b,c]
                    for e in range(4):
                        value+=(-first[m,b,e]*C[e,c,d]+first[m,c,e]*C[e,d,b]
                                +first[m,d,e]*C[e,b,c])
                    dlower[m,d,b,c]=value/2
    Gamma=np.einsum('ad,dbc->abc',inverse,lower)
    dGamma=(np.einsum('mad,dbc->mabc',dinverse,lower)
            +np.einsum('ad,mdbc->mabc',inverse,dlower))
    Ricci=np.zeros((4,4))
    for b in range(4):
        for d in range(4):
            for a in range(4):
                Ricci[b,d]+=dGamma[a,a,d,b]-dGamma[d,a,a,b]
                for c in range(4):
                    Ricci[b,d]+=(Gamma[a,a,c]*Gamma[c,d,b]-Gamma[a,d,c]*Gamma[c,a,b]
                                 -C[c,a,d]*Gamma[a,c,b])
    return inverse,dinverse,Gamma,dGamma,Ricci


def gh_constraint(g,first,inverse,dinverse,Gamma,dGamma,Gammabar,H,Hdot,source_gradient=None):
    difference=Gamma-Gammabar
    delta=np.einsum('cd,bcd->b',inverse,difference)
    ddelta=(np.einsum('mcd,bcd->mb',dinverse,difference)
            +np.einsum('cd,mbcd->mb',inverse,dGamma))
    covector=g@delta;covector[0]+=H
    derivative=np.einsum('mbc,c->mb',first,delta)+np.einsum('bc,mc->mb',g,ddelta)
    if source_gradient is None:derivative[0,0]+=Hdot
    else:derivative[:,0]+=source_gradient
    nabla=derivative-np.einsum('cab,c->ab',Gamma,covector)
    return covector,(nabla+nabla.T)/2


def symmetric_metric_jet(rng,scale):
    matrix=rng.normal(size=(4,4))*scale
    return (matrix+matrix.T)/2


def manufactured_jets(seed):
    rng=np.random.default_rng(seed)
    spatial=np.array([[1.2,.08,-.05],[.08,.95,.04],[-.05,.04,1.15]])
    shift=np.array([.17,-.11,.09]);lapse=1.13
    g=np.empty((4,4));g[1:,1:]=spatial
    g[0,1:]=g[1:,0]=spatial@shift
    g[0,0]=-lapse*lapse+shift@spatial@shift
    v=symmetric_metric_jet(rng,.045)
    dg=np.array([symmetric_metric_jet(rng,.035) for _ in range(3)])
    dv=np.array([symmetric_metric_jet(rng,.025) for _ in range(3)])
    d2g=rng.normal(size=(3,3,4,4))*.025
    d2g=(d2g+d2g.swapaxes(0,1))/2
    d2g=(d2g+d2g.swapaxes(2,3))/2
    q=np.array([.012,.021]);p=np.array([.031,-.018])
    dq=rng.normal(size=(3,2))*.017;dp=rng.normal(size=(3,2))*.023
    d2q=rng.normal(size=(3,3,2))*.021;d2q=(d2q+d2q.swapaxes(0,1))/2
    return [g,v,dg,d2g,dv,q,p,dq,d2q,dp]


def test_direct_koszul_curvature_has_the_positive_round_s3_sign():
    r=.27;g=np.diag([-1.,1.,1.,1.])
    _,_,Gamma,_,Ricci=physical_geometry(g,np.zeros((4,4,4)),np.zeros((4,4,4,4)),bracket(r))
    np.testing.assert_allclose(Ricci,np.diag([0.,2*r*r,2*r*r,2*r*r]),atol=1e-15,rtol=0)
    np.testing.assert_allclose(Gamma-Gamma.swapaxes(1,2),bracket(r),atol=1e-16,rtol=0)


@pytest.mark.parametrize('seed,lapse_feedback',[(4723,0.),(5891,0.),(4723,.27)])
def test_general_nonlinear_gh_identity_from_the_physical_ricci_tensor(seed,lapse_feedback):
    r=.27;H=.13;Hdot=.07;gamma0=.27;kappa=1.
    Cframe=bracket(r)
    jets=manufactured_jets(seed)
    g,v,dg,d2g,dv,q,p,dq,d2q,dp=jets
    # Ordered frame jets: [e_i,e_j]f=C^k_ij e_k f for every component f.
    for i in range(3):
        for j in range(3):
            for k in range(3):
                d2g[i,j]+=Cframe[k+1,i+1,j+1]*dg[k]/2
                d2q[i,j]+=Cframe[k+1,i+1,j+1]*dq[k]/2
    star=next(x for x in VolumeParent().vacuum_roots() if x['kind']=='minimum')['s']
    ag,aq,Ckernel,valid=rhs_nodes(*[x[None] for x in jets],r=r,Hsource0=H,HsourceDot0=Hdot,
                                sstar=star,kappa=kappa,gamma0=gamma0,lapse_feedback=lapse_feedback)
    assert valid[0]
    first=np.concatenate([v[None],dg],axis=0)
    second=np.empty((4,4,4,4));second[0,0]=ag[0]
    second[0,1:]=dv;second[1:,0]=dv;second[1:,1:]=d2g
    inverse,dinverse,Gamma,dGamma,Ricci=physical_geometry(g,first,second,Cframe)
    _,_,Gammabar,_,_=physical_geometry(np.diag([-1.,1.,1.,1.]),np.zeros((4,4,4)),
                                     np.zeros((4,4,4,4)),Cframe)
    # Derive the gauge-source chain from raw metric jets and the inverse,
    # independently of the kernel's normal/reference-covariant contraction.
    lapse=1/math.sqrt(-inverse[0,0])
    lapse_gradient=.5*lapse**3*dinverse[:,0,0]
    effective_H=H+lapse_feedback*(lapse*lapse-lapse**-2)
    source_gradient=np.array([Hdot,0.,0.,0.])
    source_gradient+=lapse_feedback*(2*lapse+2*lapse**-3)*lapse_gradient
    C,nablaC=gh_constraint(g,first,inverse,dinverse,Gamma,dGamma,Gammabar,
                         effective_H,Hdot,source_gradient=source_gradient)
    np.testing.assert_allclose(Ricci,Ricci.T,atol=2e-14,rtol=0)
    np.testing.assert_allclose(Gamma-Gamma.swapaxes(1,2),Cframe,atol=1e-14,rtol=0)
    np.testing.assert_allclose(C,Ckernel[0],atol=1e-14,rtol=0)
    assert np.linalg.norm(C)>.05  # This is deliberately off the GH constraint.
    normal=-lapse*inverse[:,0];ncov=g@normal
    damping=gamma0*((np.outer(ncov,C)+np.outer(C,ncov))/2-g*np.dot(normal,C)/2)
    volume=star+q[0]/math.sqrt(84)
    e12,e14,e24=math.exp(-12*volume),math.exp(-14*volume),math.exp(-24*volume)
    mass2=400*e12
    W=7*e12-7.5*e14+(25/32)*e24+mass2*q[1]*q[1]/2
    Dq=np.concatenate([p[None],dq],axis=0)
    stress_reversed=Dq@Dq.T+W*g
    residual=Ricci-nablaC+damping-kappa*stress_reversed
    np.testing.assert_allclose(residual,0.,atol=3e-12,rtol=0,
                               err_msg='direct physical Ricci/GH identity disagrees with returned acceleration')
    if lapse_feedback:
        assert np.linalg.norm(lapse_gradient[1:])>.005
        # This case really requires spatial chain terms. Retaining the time
        # chain alone gives a material failure of the same physical identity.
        missing_spatial=source_gradient.copy();missing_spatial[1:]=0.
        _,wrong_nabla=gh_constraint(g,first,inverse,dinverse,Gamma,dGamma,Gammabar,
                                   effective_H,Hdot,source_gradient=missing_spatial)
        wrong_residual=Ricci-wrong_nabla+damping-kappa*stress_reversed
        assert np.max(abs(wrong_residual))>1e-4

    # Independent canonical scalar equations use the same physical connection.
    second_q=np.empty((4,4,2));second_q[0,0]=aq[0]
    second_q[0,1:]=dp;second_q[1:,0]=dp;second_q[1:,1:]=d2q
    wave=np.einsum('ab,abf->f',inverse,second_q)-np.einsum('ab,cab,cf->f',inverse,Gamma,Dq)
    gradient=np.array([(-84*e12+105*e14-(75/4)*e24-6*mass2*q[1]*q[1])/math.sqrt(84),mass2*q[1]])
    np.testing.assert_allclose(wave,gradient,atol=2e-12,rtol=0)
