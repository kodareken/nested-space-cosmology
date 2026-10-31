"""Independent quadrature, normal-observer derivative and nodal controls."""
import numpy as np
import pytest

from recursive_horizons.nsc_atlas_gh_wave_diagnostics import wave_phase_diagnostics


def blank(batch=()):
    return dict(g=np.broadcast_to(np.diag([-1.,1.,1.,1.]),batch+(4,4)).copy(),
                dg=np.zeros(batch+(3,4,4)),q=np.zeros(batch+(2,)),p=np.zeros(batch+(2,)),
                dq=np.zeros(batch+(3,2)),d2q=np.zeros(batch+(3,3,2)),dp=np.zeros(batch+(3,2)))


def test_travelling_wave_has_constant_amplitude_and_pure_phase_gradient():
    x=np.linspace(-1.3,2.1,151);A=.31;k=1.7;omega=2.8;t=.19
    theta=k*x-omega*t;z=blank((len(x),))
    z['q'][:,1]=A*np.cos(theta)
    z['p'][:,1]=A*omega*np.sin(theta)  # partial_t cos(kx-omega*t) has this positive sign.
    z['dq'][:,0,1]=-A*k*np.sin(theta)
    z['dp'][:,0,1]=A*omega*k*np.cos(theta)
    z['d2q'][:,0,0,1]=-A*k*k*np.cos(theta)
    out=wave_phase_diagnostics(**z,omega_ref=omega)
    np.testing.assert_allclose(out['amplitude_squared'],A*A,atol=6e-17,rtol=0)
    np.testing.assert_allclose(out['reference_gradient_squared'],A*A*k*k,atol=3e-16,rtol=0)
    np.testing.assert_allclose(out['amplitude_gradient_squared'],0,atol=2e-32,rtol=0)
    np.testing.assert_allclose(out['phase_gradient_contribution'],A*A*k*k,atol=3e-16,rtol=0)
    np.testing.assert_allclose(out['phase_current'][:,0],-A*A*k,atol=1e-16,rtol=0)
    np.testing.assert_allclose(out['phase_current'][:,1:],0,atol=0,rtol=0)
    assert out['undefined_count']==0


def test_real_standing_profile_has_only_amplitude_gradient():
    x=np.linspace(-1.,1.,81);z=blank((len(x),));f=.2*np.exp(-x*x)
    z['q'][:,1]=f;z['dq'][:,0,1]=-2*x*f;z['d2q'][:,0,0,1]=(4*x*x-2)*f
    out=wave_phase_diagnostics(**z,omega_ref=4.)
    expected=4*x*x*f*f
    np.testing.assert_allclose(out['reference_gradient_squared'],expected,atol=1e-17,rtol=0)
    np.testing.assert_allclose(out['amplitude_gradient_squared'],expected,atol=1e-17,rtol=0)
    np.testing.assert_array_equal(out['phase_gradient_contribution'],0)
    np.testing.assert_array_equal(out['phase_current'],0)


def test_normal_momentum_derivative_includes_lapse_shift_and_spatial_metric_rates():
    z=blank();h=np.array([[1.3,.1,0.],[.1,.9,.2],[0.,.2,1.5]])
    dh=np.array([[[.03,.01,0.],[.01,-.02,0.],[0.,0.,.01]],
                 [[.01,0.,.02],[0.,.02,0.],[.02,0.,-.01]],
                 [[-.02,0.,0.],[0.,.01,.01],[0.,.01,.03]]])
    N=1.2;dN=np.array([.03,-.08,.04]);beta=np.array([.2,-.1,.05])
    dbeta=np.array([[.02,-.01,.03],[-.01,.04,.02],[.03,.01,-.02]])
    z['g'][1:,1:]=h;z['g'][0,1:]=h@beta;z['g'][1:,0]=h@beta
    z['g'][0,0]=-N*N+beta@h@beta
    for i in range(3):
        z['dg'][i,1:,1:]=dh[i]
        z['dg'][i,0,1:]=dh[i]@beta+h@dbeta[i]
        z['dg'][i,1:,0]=z['dg'][i,0,1:]
        z['dg'][i,0,0]=-2*N*dN[i]+2*dbeta[i]@h@beta+beta@dh[i]@beta
    chi=.4;pi=.7;dchi=np.array([.2,-.3,.1]);dpi=np.array([.06,-.04,.09])
    second=np.array([[.1,.03,-.02],[.03,-.07,.04],[-.02,.04,.02]])
    z['q'][1]=chi;z['p'][1]=N*pi+beta@dchi
    z['dq'][:,1]=dchi;z['d2q'][:,:,1]=second
    z['dp'][:,1]=dN*pi+N*dpi+dbeta@dchi+second@beta
    omega=2.3;out=wave_phase_diagnostics(**z,omega_ref=omega)
    assert out['normal_momentum']==pytest.approx(pi,abs=3e-16)
    np.testing.assert_allclose(out['normal_momentum_gradient'],dpi,atol=4e-17,rtol=0)
    Z=chi-1j*pi/omega;dZ=dchi-1j*dpi/omega
    amplitude=np.real(np.conj(Z)*dZ);phase=np.imag(np.conj(Z)*dZ)
    np.testing.assert_allclose(out['phase_current'],phase,atol=2e-17,rtol=0)
    assert out['amplitude_gradient_squared']==pytest.approx(amplitude@amplitude/abs(Z)**2)
    assert out['phase_gradient_contribution']==pytest.approx(phase@phase/abs(Z)**2)
    hi=np.linalg.inv(h)
    expected_physical=float(np.real(np.conj(dZ)@hi@dZ))
    assert out['physical_gradient_squared']==pytest.approx(expected_physical,abs=3e-17)
    assert abs(out['identity_residual'])<6e-17
    assert abs(out['physical_identity_residual'])<6e-17


def test_reference_and_physical_contractions_are_not_silently_interchanged():
    z=blank();scale=2.4;z['g'][1:,1:]=scale*scale*np.eye(3)
    z['q'][1]=.3;z['p'][1]=.4;z['dq'][:,1]=[.2,-.1,.3];z['dp'][:,1]=[-.1,.2,.4]
    out=wave_phase_diagnostics(**z,omega_ref=1.7)
    assert out['physical_gradient_squared']==pytest.approx(out['reference_gradient_squared']/scale**2)
    assert out['physical_amplitude_gradient_squared']==pytest.approx(out['amplitude_gradient_squared']/scale**2)
    assert out['physical_phase_gradient_contribution']==pytest.approx(out['phase_gradient_contribution']/scale**2)


def test_noncommuting_frame_uses_ordered_second_scalar_jets():
    r=.37;f=.23;z=blank();beta=np.array([0.,.4,0.]);pi=.19
    eps=np.zeros((3,3,3))
    for i,j,k in ((0,1,2),(1,2,0),(2,0,1)):
        eps[i,j,k]=1.;eps[j,i,k]=-1.
    df=np.array([r*np.sqrt(1-f*f),0.,0.])
    second=-r*r*f*np.eye(3)+r*np.einsum('ijk,k->ij',eps,df)
    assert np.max(abs(second-second.T))>0.1
    z['g'][0,0]=-1+beta@beta;z['g'][0,1:]=beta;z['g'][1:,0]=beta
    z['q'][1]=f;z['dq'][:,1]=df;z['d2q'][:,:,1]=second
    # Constant NORMAL pi: p_coordinate=pi+beta^j E_j chi.
    z['p'][1]=pi+beta@df;z['dp'][:,1]=second@beta
    out=wave_phase_diagnostics(**z,omega_ref=2.)
    assert out['normal_momentum']==pytest.approx(pi)
    np.testing.assert_allclose(out['normal_momentum_gradient'],0,atol=2e-17,rtol=0)


def test_exact_node_is_flagged_without_erasing_its_full_gradient():
    z=blank();z['dq'][:,1]=[1.,.2,.1];z['dp'][:,1]=[.3,.1,.2]
    out=wave_phase_diagnostics(**z,omega_ref=2.)
    assert not bool(out['defined']) and out['undefined_count']==1
    assert out['reference_gradient_squared']>1
    assert out['amplitude_gradient_squared']==0
    assert out['phase_gradient_contribution']==0
    assert out['identity_residual']==0  # Masked: identity is not asserted at this node.


def test_tiny_nonzero_quadrature_is_not_reclassified_as_an_exact_node():
    z=blank();z['q'][1]=1e-80;z['p'][1]=-.3e-80
    z['dq'][:,1]=[.2,-.1,.3];z['dp'][:,1]=[.3,.2,-.1]
    out=wave_phase_diagnostics(**z,omega_ref=2.)
    assert bool(out['defined']) and out['undefined_count']==0
    assert out['amplitude_squared']>0
    assert abs(out['identity_residual'])<1e-16


def test_other_canonical_field_does_not_enter_the_chi_quadrature():
    z=blank();z['q'][1]=.4;z['p'][1]=.7;z['dq'][:,1]=[.1,.2,.3]
    before=wave_phase_diagnostics(**z,omega_ref=2.)
    z['q'][0]=7.;z['p'][0]=-3.;z['dq'][:,0]=10.;z['d2q'][:,:,0]=5.;z['dp'][:,0]=-2.
    after=wave_phase_diagnostics(**z,omega_ref=2.)
    for key in before:np.testing.assert_array_equal(before[key],after[key])


def test_invalid_frequency_shapes_and_foliation_are_rejected():
    z=blank()
    for omega in (0.,-1.,float('nan'),True):
        with pytest.raises(ValueError):wave_phase_diagnostics(**z,omega_ref=omega)
    with pytest.raises(ValueError):wave_phase_diagnostics(**dict(z,dp=np.zeros(2)),omega_ref=2.)
    z['g'][0,0]=1.
    with pytest.raises(ValueError):wave_phase_diagnostics(**z,omega_ref=2.)
