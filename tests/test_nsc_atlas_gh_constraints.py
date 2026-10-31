"""Independent exact-geometry controls for the raw-frame ADM diagnostics."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_gh_constraints import adm_constraints

EPS=np.zeros((3,3,3))
for i,j,k in ((0,1,2),(1,2,0),(2,0,1)):
    EPS[i,j,k]=1.;EPS[j,i,k]=-1.


def blank():
    return dict(g=np.diag([-1.,1.,1.,1.]),v=np.zeros((4,4)),dg=np.zeros((3,4,4)),
                d2g=np.zeros((3,3,4,4)),dv=np.zeros((3,4,4)),
                q=np.zeros(2),p=np.zeros(2),dq=np.zeros((3,2)))


def coordinate_function_jet(r,f=.23):
    # Hess_bar f=-r^2 f gamma for one unit-S3 embedding coordinate.
    df=np.array([r*math.sqrt(1-f*f),0.,0.])
    d2f=-r*r*f*np.eye(3)+r*np.einsum('ijk,k->ij',EPS,df)
    np.testing.assert_allclose(d2f-d2f.T,2*r*np.einsum('ijk,k->ij',EPS,df),atol=2e-16)
    return f,df,d2f


def conformal_jet(r,alpha=.13):
    f,df,d2f=coordinate_function_jet(r)
    psi=1+alpha*f*f;dpsi=2*alpha*f*df
    d2psi=2*alpha*(np.outer(df,df)+f*d2f)
    z=blank();z['g'][1:,1:]=psi**4*np.eye(3)
    for i in range(3):
        z['dg'][i,1:,1:]=4*psi**3*dpsi[i]*np.eye(3)
        for j in range(3):
            z['d2g'][i,j,1:,1:]=(12*psi**2*dpsi[i]*dpsi[j]+4*psi**3*d2psi[i,j])*np.eye(3)
    return z,psi,dpsi,d2psi


def test_exact_closed_de_sitter_background_and_canonical_vacuum():
    H=.2714056154003984;kappa=1.7
    for u in (0.,.4,1.2):
        z=blank();c=math.cosh(u);h=H*math.tanh(u)
        z['g'][1:,1:]=c*c*np.eye(3);z['v'][1:,1:]=2*h*c*c*np.eye(3)
        out=adm_constraints(**z,r=H,kappa=kappa,potential=3*H*H/kappa)
        np.testing.assert_allclose(out['Kij'],-h*c*c*np.eye(3),atol=2e-15)
        assert out['R3']==pytest.approx(6*H*H/(c*c),abs=2e-15)
        assert abs(out['H'])<3e-15
        np.testing.assert_allclose(out['M'],0,atol=2e-15)


def test_finite_conformal_metric_matches_its_independent_scalar_curvature():
    r=.31;z,psi,dpsi,d2psi=conformal_jet(r)
    out=adm_constraints(**z,r=r,potential=3*r*r)
    R=psi**-5*(-8*np.trace(d2psi)+6*r*r*psi)
    assert out['R3']==pytest.approx(R,abs=2e-15)
    np.testing.assert_allclose(out['Kij'],0,atol=1e-16)
    assert out['H']==pytest.approx(R-6*r*r,abs=2e-15)
    np.testing.assert_allclose(out['M'],0,atol=1e-16)


def test_initial_GH_time_shift_rates_do_not_change_the_intrinsic_data():
    r=.37;z,psi,dpsi,d2psi=conformal_jet(r)
    baseline=adm_constraints(**z,r=r,potential=.23)
    z['v'][0,1:]=-2*dpsi/psi;z['v'][1:,0]=z['v'][0,1:]
    z['dv'][:,0,1:]=-2*(d2psi/psi-np.outer(dpsi,dpsi)/psi**2)
    z['dv'][:,1:,0]=z['dv'][:,0,1:]
    out=adm_constraints(**z,r=r,potential=.23)
    for key in ('Kij','R3','H','M'):
        np.testing.assert_allclose(out[key],baseline[key],atol=1e-15,rtol=0)


def test_gradient_of_isotropic_extrinsic_curvature_fixes_momentum_sign():
    r=.42;kappa=1.3;f,df,_=coordinate_function_jet(r);z=blank()
    z['v'][1:,1:]=-2*f*np.eye(3)
    for m in range(3):z['dv'][m,1:,1:]=-2*df[m]*np.eye(3)
    z['q'][:]=[.1,.2];z['p'][:]=[.3,.7]
    z['dq'][:]=[[.2,.1],[-.1,.3],[.4,-.2]]
    potential=.51;out=adm_constraints(**z,r=r,kappa=kappa,potential=potential)
    j=-z['dq']@z['p']
    rho=(np.dot(z['p'],z['p'])+np.sum(z['dq']**2))/2+potential
    np.testing.assert_allclose(out['Kij'],f*np.eye(3),atol=2e-16)
    np.testing.assert_allclose(out['j'],j,atol=2e-16)
    np.testing.assert_allclose(out['M'],-2*df-kappa*j,atol=2e-15)
    assert out['rho']==pytest.approx(rho)
    assert out['H']==pytest.approx(6*r*r+6*f*f-2*kappa*rho,abs=2e-15)


def test_nonuniform_lapse_cancels_from_a_prescribed_physical_K():
    r=.25;f,df,d2f=coordinate_function_jet(r);z=blank();b=.4;N=math.exp(b*f)
    dN=b*N*df;d2N=N*(b*d2f+b*b*np.outer(df,df))
    z['g'][0,0]=-N*N;z['dg'][:,0,0]=-2*N*dN
    z['d2g'][:,:,0,0]=-2*(np.outer(dN,dN)+N*d2N)
    z['v'][1:,1:]=-2*N*f*np.eye(3)
    for m in range(3):z['dv'][m,1:,1:]=-2*(N*df[m]+f*dN[m])*np.eye(3)
    out=adm_constraints(**z,r=r,potential=0.)
    assert out['lapse']==pytest.approx(N)
    np.testing.assert_allclose(out['Kij'],f*np.eye(3),atol=3e-16)
    np.testing.assert_allclose(out['M'],-2*df,atol=2e-15)


def test_nonzero_shift_Lie_derivative_gives_K_zero_without_a_coordinate_chart():
    r=.33;z,psi,dpsi,d2psi=conformal_jet(r)
    beta=np.array([2.,.3,-.2]);b2=np.dot(beta,beta);scale=psi**4
    ds=4*psi**3*dpsi;d2s=12*psi**2*np.outer(dpsi,dpsi)+4*psi**3*d2psi
    z['g'][0,1:]=scale*beta;z['g'][1:,0]=z['g'][0,1:];z['g'][0,0]=-1+scale*b2
    for i in range(3):
        z['dg'][i,0,1:]=ds[i]*beta;z['dg'][i,1:,0]=z['dg'][i,0,1:]
        z['dg'][i,0,0]=ds[i]*b2
        for j in range(3):
            z['d2g'][i,j,0,1:]=d2s[i,j]*beta;z['d2g'][i,j,1:,0]=z['d2g'][i,j,0,1:]
            z['d2g'][i,j,0,0]=d2s[i,j]*b2
    velocity=np.dot(beta,ds);dvelocity=d2s@beta
    z['v'][1:,1:]=velocity*np.eye(3);z['v'][0,1:]=velocity*beta
    z['v'][1:,0]=z['v'][0,1:];z['v'][0,0]=velocity*b2
    for i in range(3):
        z['dv'][i,1:,1:]=dvelocity[i]*np.eye(3)
        z['dv'][i,0,1:]=dvelocity[i]*beta;z['dv'][i,1:,0]=z['dv'][i,0,1:]
        z['dv'][i,0,0]=dvelocity[i]*b2
    # A scalar advected by this shift has zero NORMAL momentum.
    z['dq'][:]=[[.1,.2],[.2,-.1],[-.1,.3]];z['p']=beta@z['dq']
    out=adm_constraints(**z,r=r,potential=.2)
    assert z['g'][0,0]>0  # Coordinate time lines may be spacelike; slices remain spacelike.
    assert out['lapse']==pytest.approx(1.,abs=2e-15)
    np.testing.assert_allclose(out['shift'],beta,atol=2e-15)
    np.testing.assert_allclose(out['Kij'],0,atol=3e-16)
    np.testing.assert_allclose(out['j'],0,atol=3e-16)
    np.testing.assert_allclose(out['M'],0,atol=3e-15)
    R=psi**-5*(-8*np.trace(d2psi)+6*r*r*psi)
    assert out['R3']==pytest.approx(R,abs=2e-15)


def test_batch_shapes_and_flat_zero_source_limit():
    z=blank();batch={name:np.repeat(value[None,...],3,axis=0) for name,value in z.items()}
    out=adm_constraints(**batch,r=0.,potential=np.zeros(3))
    assert out['Kij'].shape==(3,3,3)
    np.testing.assert_array_equal(out['H'],0.)
    np.testing.assert_array_equal(out['M'],0.)
    np.testing.assert_array_equal(out['rho'],0.)


def test_anisotropic_SU2_metric_and_constant_frame_rotation():
    # Cartan/Bianchi-IX scalar curvature for h=diag(A,B,C), [E_i,E_j]=2r eps_ijk E_k.
    r=.31;A,B,C=1.2,.8,1.7;k=.19;kappa=1.4;W=.23
    R=2*r*r*(-A*A-B*B-C*C+2*(A*B+B*C+C*A))/(A*B*C)
    angle=.53;c,s=math.cos(angle),math.sin(angle)
    rotation=np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]])
    for O in (np.eye(3),rotation):
        z=blank();h=O@np.diag([A,B,C])@O.T
        z['g'][1:,1:]=h;z['v'][1:,1:]=-2*k*h
        z['p'][:]=[.1,.2];z['dq'][:]=[[.3,-.1],[-.2,.4],[.1,.2]]
        out=adm_constraints(**z,r=r,kappa=kappa,potential=W)
        grad=np.einsum('ij,ia,ja->',np.linalg.inv(h),z['dq'],z['dq'])
        rho=(np.dot(z['p'],z['p'])+grad)/2+W
        assert out['R3']==pytest.approx(R,abs=3e-15)
        np.testing.assert_allclose(out['Kij'],k*h,atol=3e-16)
        assert out['rho']==pytest.approx(rho,abs=3e-16)
        assert out['H']==pytest.approx(R+6*k*k-2*kappa*rho,abs=3e-15)
        np.testing.assert_allclose(out['M'],-kappa*out['j'],atol=3e-15)


def test_bad_shapes_invalid_slices_and_metric_asymmetry_are_rejected():
    z=blank()
    with pytest.raises(ValueError):adm_constraints(**dict(z,d2g=np.zeros((3,4,4))),r=1.,potential=0.)
    bad=blank();bad['g'][0,0]=1.
    with pytest.raises(ValueError):adm_constraints(**bad,r=1.,potential=0.)
    bad=blank();bad['g'][1,2]=.1
    with pytest.raises(ValueError):adm_constraints(**bad,r=1.,potential=0.)
