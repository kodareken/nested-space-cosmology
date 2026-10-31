"""Physical and independent geometric controls for finite neck formation."""
import math

import numpy as np
from numpy.polynomial.legendre import leggauss
import pytest

from recursive_horizons.nsc_atlas_valve_spherical import SphericalValve


def test_exact_closed_de_sitter_recovery_in_proper_time():
    s=SphericalValve(24,epsilon=0.);sol=s.evolve(.3)
    y=sol.y[:,-1].reshape(8,s.modes)
    np.testing.assert_allclose(y[0],math.cosh(s.H*.3)/s.H,rtol=2e-12,atol=1e-13)
    np.testing.assert_allclose(y[2],-s.H*math.tanh(s.H*.3),atol=2e-13)
    np.testing.assert_allclose(y[[1,3,4,5,6,7]],0.,atol=3e-13)
    r=s.diagnostic(.3,y)
    assert not r['closed_neck_future_trapped']
    assert r['Hamiltonian_RMS_over_H_squared']<2e-10
    assert r['mean_w']==pytest.approx(-1.,abs=3e-14)


def test_regular_poles_and_independent_warped_metric_ricci():
    s=SphericalValve(24)
    u=np.linspace(-.9,.9,15)
    y=s.initial.copy();y[1]=.013*(1+s.u**2);y[3]=.02*(1-.1*s.u**2)
    z=s.pieces(y,u);b=z['b'];q=z['q'];bu=z['yu'][0];cu=z['yu'][1]
    _,_,second=s.jets(y,u);qu=bu-2*u*z['y'][1]+(1-u*u)*cu
    R=np.sqrt(1-u*u)*b
    Rchi=u*b-(1-u*u)*bu
    Rchichi=-np.sqrt(1-u*u)*(b+3*u*bu-(1-u*u)*second[0])
    qchi=-np.sqrt(1-u*u)*qu
    Rll=(Rchichi-qchi/q*Rchi)/q**2
    direct_rad=-2*Rll/R
    direct_tan=(1-(Rchi/q)**2-R*Rll)/R**2
    np.testing.assert_allclose(z['Rr'],direct_rad,atol=3e-15)
    np.testing.assert_allclose(z['Rt'],direct_tan,atol=3e-15)
    np.testing.assert_allclose((direct_rad-direct_tan)/(1-u*u),z['Rdiff'],atol=4e-15)
    poles=s.pieces(y,[-1.,1.])
    np.testing.assert_allclose(poles['q'],poles['b'],atol=0)
    np.testing.assert_allclose(poles['ka'],poles['kb'],atol=0)
    np.testing.assert_allclose(poles['Rr'],poles['Rt'],atol=5e-15)


def test_convergent_actual_future_trapped_neck_and_marginal_boundaries():
    rows=[]
    for N in (24,32,48):
        s=SphericalValve(N);sol=s.evolve(.3);y=sol.y[:,-1].reshape(8,N)
        r=s.diagnostic(.3,y);rows.append(r)
        assert r['closed_neck_future_trapped']
        assert r['theta_plus_over_H']<-.1 and r['theta_minus_over_H']<-.1
        assert r['neck_radius_gradient_norm_squared']<0
        assert abs(r['neck_timelike_energy_flux'])<1e-14
        band=s.trapped_band(y)
        assert 0<band['u_edge']<.3
        assert abs(band['marginal_expansion'])<1e-12
        assert band['other_expansion']<0
        assert 0<band['volume']<r['volume']
        assert 0<band['energy']<r['energy']
    assert rows[2]['Hamiltonian_RMS_over_H_squared']<2e-10
    assert rows[1]['Hamiltonian_RMS_over_H_squared']<rows[0]['Hamiltonian_RMS_over_H_squared']/50
    assert rows[2]['Hamiltonian_RMS_over_H_squared']<rows[1]['Hamiltonian_RMS_over_H_squared']/50
    assert abs(rows[2]['theta_plus_over_H']-rows[1]['theta_plus_over_H'])<2e-9


def test_regional_energy_change_is_pressure_work_not_reflection_flux():
    s=SphericalValve(32);sol=s.evolve(.3)
    first=s.diagnostic(0.,s.initial);last=s.diagnostic(.3,sol.y[:,-1].reshape(8,s.modes))
    x,w=leggauss(32);times=.15*(x+1)
    integral=.15*sum(weight*s.diagnostic(float(t),sol.sol(t).reshape(8,s.modes))['pressure_work_rate'] for t,weight in zip(times,w))
    assert abs(last['energy']-first['energy']-integral)<2e-8
    assert last['half_energy']-first['half_energy']>0
    assert abs(last['neck_timelike_energy_flux'])<1e-14


def test_preparation_below_neck_threshold_does_not_form_initial_neck_trap():
    s=SphericalValve(32,epsilon=.75);sol=s.evolve(.02)
    r=s.diagnostic(.02,sol.y[:,-1].reshape(8,s.modes))
    assert r['theta_plus_over_H']>0 and r['theta_minus_over_H']>0
    assert not r['closed_neck_future_trapped']
    assert r['Hamiltonian_RMS_over_H_squared']<1e-8


def test_future_outer_horizons_and_symmetric_receiving_flux():
    s=SphericalValve(48);sol=s.evolve(.3);t=.2
    def state(at):return sol.sol(at).reshape(8,s.modes)
    band=s.trapped_band(state(t));influx=0.
    for sign in (1.,-1.):
        u=sign*band['u_edge'];z=s.pieces(state(t),[u]);R=float(z['radius'][0])
        q=float(z['q'][0]);U=-float(z['kb'][0])*R;radial=float(z['rr'][0])/q
        mass=R/2*(1+U*U-radial*radial)
        kappa_H=mass/R**2-s.kappa*R*(float(z['rho'][0])-float(z['pr'][0]))/4
        sigma=-sign
        def theta_out(at,coord):
            return float(s.expansions(state(at),[coord])[0 if sigma>0 else 1][0])
        h=1e-5
        cross=(theta_out(t+h,u)-theta_out(t-h,u))/(2*h)
        cross+=sigma*math.sqrt(1-u*u)/q*(theta_out(t,u+h)-theta_out(t,u-h))/(2*h)
        assert kappa_H>0 and cross<0
        np.testing.assert_allclose(cross,-4*kappa_H/R,atol=2e-9,rtol=0)
        influx-=4*math.pi*R*R*sigma*float(z['j'][0])
    assert influx>0
    np.testing.assert_allclose(influx,.08561978764,rtol=1e-8)
    assert abs(float(s.pieces(state(t),[0.])['j'][0]))<1e-14
