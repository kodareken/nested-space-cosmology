import numpy as np
import pytest
from recursive_horizons.nsc_atlas_flow_evolution import PeriodicFlow


def test_beltrami_decay_and_energy_balance():
    f=PeriodicFlow(16,length=2*np.pi,viscosity=.1,workers=1)
    z=2*np.pi*np.arange(16)/16
    u=np.zeros((3,16,16,16));u[0]=np.sin(z);u[1]=np.cos(z)
    h=f.project(f.fft(u));rate,D=f.rhs(h)
    np.testing.assert_allclose(rate,-.1*h,atol=2e-12)
    initial=f.diagnostics(h)['energy'];loss=0.
    for _ in range(10):h,inc=f.step(h,.001);loss+=inc
    np.testing.assert_allclose(f.real(h),u*np.exp(-.001),atol=2e-14)
    assert abs(f.diagnostics(h)['energy']-initial+loss)<1e-10


def test_taylor_green_pressure_is_recovered_without_fitting():
    f=PeriodicFlow(16,length=2*np.pi,viscosity=.1,workers=1)
    a=2*np.pi*np.arange(16)/16;x,y,z=np.meshgrid(a,a,a,indexing='ij')
    u=np.array([np.sin(x)*np.cos(y),-np.cos(x)*np.sin(y),np.zeros_like(x)])
    h=f.project(f.fft(u));rate,_=f.rhs(h)
    np.testing.assert_allclose(f.real(rate),-.2*u,atol=3e-15)
    # u.grad(u)=(sin(2x),sin(2y),0)/2, cancelled by this pressure gradient.
    np.testing.assert_allclose(f.pressure(h),(np.cos(2*x)+np.cos(2*y))/4,atol=6e-16)


def test_prepared_flow_is_solenoidal_and_nonlinearity_conserves_energy():
    f=PeriodicFlow(24,workers=1);h=f.prepare();u=f.real(h)
    rate,D=f.rhs(h)
    direct=f.length**3*np.mean(np.sum(u*f.real(rate),axis=0))
    assert direct==pytest.approx(-D,rel=2e-12)
    assert f.diagnostics(h)['divergence_rms']<3e-15
