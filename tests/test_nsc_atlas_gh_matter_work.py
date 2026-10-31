"""Independent stress-work signs, first-jet geometry and local flux controls."""
import math

import numpy as np
import pytest
from numba import set_num_threads

from recursive_horizons.nsc_atlas_gh_matter_work import matter_work


def blank(batch=()):
    return dict(g=np.broadcast_to(np.diag([-1.,1.,1.,1.]),batch+(4,4)).copy(),
                v=np.zeros(batch+(4,4)),dg=np.zeros(batch+(3,4,4)),
                q=np.zeros(batch+(2,)),p=np.zeros(batch+(2,)),dq=np.zeros(batch+(3,2)))


def test_homogeneous_frw_scalar_work_is_negative_pressure_volume_work():
    set_num_threads(1)
    z=blank();a=1.7;H=.29;W=.31
    z['g'][1:,1:]=a*a*np.eye(3);z['v'][1:,1:]=2*H*a*a*np.eye(3)
    z['q'][:]=[.2,.3];z['p'][:]=[.7,-.4]
    out=matter_work(**z,r=.37,potential=W)
    kinetic=.5*z['p']@z['p'];pressure=kinetic-W
    assert out['rho']==pytest.approx(kinetic+W,abs=2e-16)
    assert out['work_integrand']==pytest.approx(-3*H*pressure,abs=2e-16)
    assert out['volume_weight']==pytest.approx(a**3,abs=1e-15)
    assert out['lapse']==pytest.approx(1.,abs=1e-16)
    assert out['extrinsic_trace']==pytest.approx(-3*H,abs=2e-16)
    np.testing.assert_array_equal(out['energy_flux'],0.)
    # A pure full vacuum gives positive normal-energy growth as its volume grows.
    z['p'][:]=0.;vacuum=matter_work(**z,r=.37,potential=W)
    assert vacuum['work_integrand']==pytest.approx(3*H*W,abs=2e-16)


def test_constant_vacuum_subtraction_matches_direct_FRW_excitation_energy_rate():
    z=blank();a=1.6;H=.23;N=1.4;Vstar=.37;masses2=np.array([.7,1.9])
    z['g'][0,0]=-N*N;z['g'][1:,1:]=a*a*np.eye(3)
    z['v'][1:,1:]=2*N*H*a*a*np.eye(3)
    z['q'][:]=[.2,-.3];pi=np.array([.6,-.4]);z['p']=N*pi
    U=.5*np.sum(masses2*z['q']**2);out=matter_work(**z,r=.31,potential=Vstar+U)
    excitation_work=out['work_integrand']+out['lapse']*Vstar*out['extrinsic_trace']
    # Proper-time KG: pi_dot=-3H*pi-dU/dq.  Differentiate a^3*rho_exc
    # in coordinate time, independently of the metric-work contraction.
    pi_coordinate_rate=N*(-3*H*pi-masses2*z['q'])
    rho_exc=.5*pi@pi+U
    density_rate=pi@pi_coordinate_rate+(masses2*z['q'])@z['p']
    direct_energy_rate=a**3*(density_rate+3*N*H*rho_exc)
    assert excitation_work==pytest.approx(-3*N*H*(.5*pi@pi-U),abs=2e-16)
    assert out['volume_weight']*excitation_work==pytest.approx(direct_energy_rate,abs=5e-16)
    assert out['lapse']==pytest.approx(N,abs=3e-16)
    assert out['extrinsic_trace']==pytest.approx(-3*H,abs=2e-16)


def test_flat_massless_travelling_wave_has_zero_work_and_local_flux_balance():
    A=.4;k=1.3;x=.27;t=.19
    def at(x,t):
        z=blank();phase=k*(x-t);z['q'][1]=A*math.cos(phase)
        z['p'][1]=A*k*math.sin(phase);z['dq'][0,1]=-A*k*math.sin(phase)
        return matter_work(**z,r=0.,potential=0.)
    out=at(x,t);expected=A*A*k*k*math.sin(k*(x-t))**2
    assert out['rho']==pytest.approx(expected,abs=1e-17)
    assert out['energy_flux'][0]==pytest.approx(expected,abs=1e-17)
    assert out['work_integrand']==0
    dt=2e-5
    energy_rate=(at(x,t+dt)['rho']-at(x,t-dt)['rho'])/(2*dt)
    flux_div=(at(x+dt,t)['energy_flux'][0]-at(x-dt,t)['energy_flux'][0])/(2*dt)
    assert abs(energy_rate+flux_div)<2e-12


def test_nonuniform_lapse_and_shift_give_prescribed_K_and_signed_lapse_work():
    z=blank();h=np.array([[1.4,.1,.05],[.1,.9,-.07],[.05,-.07,1.2]])
    hi=np.linalg.inv(h);N=1.3;beta=np.array([.3,-.2,.1]);dN=np.array([.2,-.1,.07])
    dbeta=np.array([[.03,.02,-.01],[-.04,.01,.03],[.01,-.02,.02]])
    K=np.array([[.07,.02,-.01],[.02,-.04,.03],[-.01,.03,.05]])
    z['g'][1:,1:]=h;z['g'][0,1:]=h@beta;z['g'][1:,0]=h@beta
    z['g'][0,0]=-N*N+beta@h@beta
    for i in range(3):
        z['dg'][i,0,1:]=h@dbeta[i];z['dg'][i,1:,0]=z['dg'][i,0,1:]
        z['dg'][i,0,0]=-2*N*dN[i]+2*dbeta[i]@h@beta
    z['v'][1:,1:]=dbeta@h+h@dbeta.T-2*N*K
    z['v'][0,1:]=[.7,-.3,.2];z['v'][1:,0]=z['v'][0,1:];z['v'][0,0]=.9
    z['dq'][:]=[[.2,-.3],[.1,.4],[-.2,.15]];pi=np.array([.5,-.7]);z['p']=N*pi+beta@z['dq']
    W=.23;grad2=np.einsum('ij,iA,jA->',hi,z['dq'],z['dq'])
    rho=.5*(pi@pi+grad2)+W
    Scov=z['dq']@z['dq'].T+h*(.5*(pi@pi-grad2)-W)
    jup=hi@(-z['dq']@pi);expected=N*np.sum(K*(hi@Scov@hi))-jup@dN
    out=matter_work(**z,r=0.,potential=W)
    assert out['rho']==pytest.approx(rho,abs=2e-16)
    assert out['work_integrand']==pytest.approx(expected,abs=2e-16)
    np.testing.assert_allclose(out['energy_flux'],math.sqrt(np.linalg.det(h))*(N*jup-beta*rho),atol=4e-16,rtol=0)


def test_regular_s3_first_jet_work_agrees_with_independent_ADM_K_and_source():
    from recursive_horizons.nsc_atlas_gh_constraints import adm_constraints
    rng=np.random.default_rng(123);z=blank((4,));W=np.linspace(.1,.4,4);r=.41
    for n in range(4):
        B=rng.normal(scale=.15,size=(3,3));h=np.eye(3)+B@B.T
        beta=rng.normal(scale=.2,size=3);N=1.1+.1*n
        z['g'][n,1:,1:]=h;z['g'][n,0,1:]=h@beta;z['g'][n,1:,0]=h@beta
        z['g'][n,0,0]=-N*N+beta@h@beta
    for key in ('v','dg'):
        B=rng.normal(scale=.1,size=z[key].shape);z[key][:]=(B+np.swapaxes(B,-1,-2))/2
    z['q'][:]=rng.normal(scale=.1,size=z['q'].shape)
    z['p'][:]=rng.normal(scale=.2,size=z['p'].shape)
    z['dq'][:]=rng.normal(scale=.1,size=z['dq'].shape)
    out=matter_work(**z,r=r,potential=W)
    adm=adm_constraints(**z,d2g=np.zeros((4,3,3,4,4)),dv=np.zeros((4,3,4,4)),r=r,potential=W)
    for n in range(4):
        hi=np.linalg.inv(z['g'][n,1:,1:]);gi=np.linalg.inv(z['g'][n]);N=adm['lapse'][n]
        pi=(z['p'][n]-adm['shift'][n]@z['dq'][n])/N
        grad2=np.einsum('ij,iA,jA->',hi,z['dq'][n],z['dq'][n])
        Scov=z['dq'][n]@z['dq'][n].T+z['g'][n,1:,1:]*(.5*(pi@pi-grad2)-W[n])
        dN=-N**3*np.einsum('a,iab,b->i',gi[0],z['dg'][n],gi[:,0])/2
        expected=N*np.sum(adm['Kij'][n]*(hi@Scov@hi))-(hi@adm['j'][n])@dN
        assert out['rho'][n]==pytest.approx(adm['rho'][n],abs=3e-16)
        assert out['work_integrand'][n]==pytest.approx(expected,abs=2e-16)


def test_varying_lapse_term_sign_can_be_isolated_with_zero_extrinsic_curvature():
    z=blank();N=1.6;dN=np.array([.3,-.1,.2]);z['g'][0,0]=-N*N
    z['dg'][:,0,0]=-2*N*dN;pi=np.array([.4,-.6]);z['p']=N*pi
    z['dq'][:]=[[.2,.1],[-.3,.2],[.1,-.4]]
    out=matter_work(**z,r=.3,potential=.2)
    assert out['work_integrand']==pytest.approx((z['dq']@pi)@dN,abs=3e-17)


def test_multiaxis_batches_and_broadcast_vacuum_potential():
    z=blank((2,3));out=matter_work(**z,r=0.,potential=np.array([[.2],[.7]]))
    np.testing.assert_array_equal(out['rho'],np.broadcast_to(np.array([[.2],[.7]]),(2,3)))
    np.testing.assert_array_equal(out['work_integrand'],0.)
    np.testing.assert_array_equal(out['energy_flux'],0.)
    np.testing.assert_array_equal(out['volume_weight'],1.)


def test_invalid_jets_frame_scale_and_metric_are_rejected():
    z=blank()
    for r in (-1.,float('nan'),True):
        with pytest.raises(ValueError):matter_work(**z,r=r,potential=0.)
    for W in (float('nan'),1j,True,[.1,.2]):
        with pytest.raises(ValueError):matter_work(**z,r=1.,potential=W)
    with pytest.raises(ValueError):matter_work(**dict(z,dg=np.zeros((4,4))),r=1.,potential=0.)
    bad=z['g'].copy();bad[0,0]=1.
    with pytest.raises(ValueError):matter_work(**dict(z,g=bad),r=1.,potential=0.)
    bad=z['g'].copy();bad[1,1]=-1.
    with pytest.raises(ValueError):matter_work(**dict(z,g=bad),r=1.,potential=0.)
