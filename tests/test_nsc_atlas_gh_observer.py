"""Independent stress-tensor, observer-normalization and phase controls."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_gh_observer import observer_readout


def blank(batch=()):
    return dict(g=np.broadcast_to(np.diag([-1.,1.,1.,1.]),batch+(4,4)).copy(),
                q=np.zeros(batch+(2,)),p=np.zeros(batch+(2,)),dq=np.zeros(batch+(3,2)),
                u=np.broadcast_to([1.,0.,0.,0.],batch+(4,)).copy())


def direct_stress(g,p,dq,W):
    gradients=np.vstack([p,dq])
    kinetic=gradients@gradients.T
    return kinetic-g*(.5*np.sum(np.linalg.inv(g)*kinetic)+W)


def test_boosted_observer_measures_both_canonical_fields_not_adm_normal_density():
    z=blank();z['q'][:]=[.7,.3];z['p'][:]=[.4,-.6]
    z['dq'][:]=[[.2,.3],[-.1,.5],[.05,-.2]]
    v=.6;z['u'][:]=np.array([1.,v,0.,0.])/math.sqrt(1-v*v)
    W=.17;T=direct_stress(z['g'],z['p'],z['dq'],W)
    out=observer_readout(**z,potential=W,omega_ref=2.3)
    assert out['density_total']==pytest.approx(z['u']@T@z['u'],abs=2e-16)
    assert abs(out['density_total']-T[0,0])>.1
    assert out['chi_proper_derivative']==pytest.approx(z['u']@np.r_[z['p'][1],z['dq'][:,1]])
    np.testing.assert_array_equal(out['unit_velocity'],z['u'])


def test_nonunit_input_normalizes_density_derivative_and_phase_consistently():
    z=blank();z['q'][1]=-.31;z['p'][:]=[.2,.8];z['dq'][:,1]=[.4,-.2,.1]
    z['u'][:]=np.array([1.,.3,-.2,0.])/math.sqrt(1-.3**2-.2**2)
    expected=observer_readout(**z,potential=.23,omega_ref=1.7)
    factor=2.4;z['u']*=factor
    actual=observer_readout(**z,potential=.23,omega_ref=1.7)
    for key in ('unit_velocity','density_total','chi_proper_derivative','quadrature_amplitude','phase_wrapped'):
        np.testing.assert_allclose(actual[key],expected[key],atol=3e-16,rtol=0)
    assert actual['input_norm_error']==pytest.approx(factor*factor-1)


def test_shifted_metric_readout_agrees_with_direct_tensor_contraction():
    z=blank();h=np.array([[1.7,.2,.1],[.2,.9,0.],[.1,0.,1.3]])
    lapse=1.4;beta=np.array([.3,-.2,.1]);velocity=np.array([.2,.1,-.07])
    gamma=1/math.sqrt(1-velocity@h@velocity)
    z['g'][1:,1:]=h;z['g'][0,1:]=h@beta;z['g'][1:,0]=h@beta
    z['g'][0,0]=-lapse*lapse+beta@h@beta
    z['u'][0]=gamma/lapse;z['u'][1:]=gamma*(velocity-beta/lapse)
    z['q'][1]=.19;z['p'][:]=[.3,-.4];z['dq'][:]=[[.1,.2],[-.3,.4],[.2,-.1]]
    W=.37;T=direct_stress(z['g'],z['p'],z['dq'],W)
    out=observer_readout(**z,potential=W,omega_ref=3.)
    assert out['input_norm_error']<6e-16
    assert out['density_total']==pytest.approx(z['u']@T@z['u'],abs=2e-16)
    normal=np.r_[1.,-beta]/lapse
    assert abs(out['density_total']-normal@T@normal)>1e-3


def test_exact_zero_quadrature_flags_undefined_phase_without_invalidating_geometry():
    z=blank();z['dq'][:,1]=[.2,-.3,.5]
    out=observer_readout(**z,potential=.1,omega_ref=2.)
    assert not bool(out['defined']) and out['quadrature_amplitude']==0
    assert out['phase_wrapped']==0 and out['density_total']>.1
    assert out['input_norm_error']==0


def test_arbitrary_leading_batch_and_broadcast_potential_preserve_wrapped_phase():
    z=blank((2,3));phases=np.array([[-math.pi,-2.1,-.4],[.2,1.7,math.pi]])
    A=.7;omega=2.8
    z['q'][...,1]=A*np.cos(phases);z['p'][...,1]=-A*omega*np.sin(phases)
    W=np.array([[.1],[.2]])
    out=observer_readout(**z,potential=W,omega_ref=omega)
    np.testing.assert_allclose(out['phase_wrapped'],phases,atol=3e-16,rtol=0)
    np.testing.assert_allclose(out['quadrature_amplitude'],A,atol=2e-16,rtol=0)
    assert out['density_total'].shape==(2,3) and out['unit_velocity'].shape==(2,3,4)
    assert np.all(out['defined'])
    assert np.all(out['phase_wrapped']>=-math.pi) and np.all(out['phase_wrapped']<=math.pi)
    expected=.5*z['p'][...,1]**2+W
    np.testing.assert_allclose(out['density_total'],expected,atol=3e-16,rtol=0)


def test_tiny_nonzero_quadrature_is_not_erased_by_an_amplitude_threshold():
    z=blank();z['q'][1]=1e-200;z['p'][1]=3e-200
    out=observer_readout(**z,potential=0.,omega_ref=2.)
    assert bool(out['defined']) and out['quadrature_amplitude']>0
    assert out['phase_wrapped']==pytest.approx(math.atan2(-1.5,1.))


def test_invalid_velocity_metric_frequency_shapes_or_potential_are_rejected():
    z=blank()
    for u in ([0.,0.,0.,0.],[-1.,0.,0.,0.],[1.,1.,0.,0.],[1.,2.,0.,0.],
              [float('nan'),0.,0.,0.]):
        with pytest.raises(ValueError):observer_readout(**dict(z,u=u),potential=0.,omega_ref=2.)
    for omega in (0.,-1.,float('nan'),True):
        with pytest.raises(ValueError):observer_readout(**z,potential=0.,omega_ref=omega)
    for W in (float('inf'),1j,[.2,.3],True):
        with pytest.raises(ValueError):observer_readout(**z,potential=W,omega_ref=2.)
    with pytest.raises(ValueError):observer_readout(**dict(z,dq=np.zeros((2,3))),potential=0.,omega_ref=2.)
    with pytest.raises(ValueError):observer_readout(**dict(z,g=np.eye(4)),potential=0.,omega_ref=2.)
    nonsymmetric=z['g'].copy();nonsymmetric[0,1]=.1
    with pytest.raises(ValueError):observer_readout(**dict(z,g=nonsymmetric),potential=0.,omega_ref=2.)
    bad_spatial=z['g'].copy();bad_spatial[1,1]=-1.
    with pytest.raises(ValueError):observer_readout(**dict(z,g=bad_spatial),potential=0.,omega_ref=2.)
