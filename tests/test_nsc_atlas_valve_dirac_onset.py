"""Actual-parent onset, exact-slice chart and basis-definition controls."""
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_valve_neck import solve_necks
from recursive_horizons.nsc_atlas_valve_spherical import SphericalValve
from recursive_horizons.nsc_atlas_valve_dirac_modes import (
    calculate_dirac_level,initial_parent_mode_pump,finite_parent_mode_pump_check,slice_isotropic_map,
)


@pytest.fixture(scope='module')
def data():
    initial=solve_necks(40,(1.,))[0]
    _,modes=calculate_dirac_level(24,initial=initial)
    return initial,next(m for m in modes if m['sign']>0)


def test_round_and_massless_limits_are_independent_onset_controls(data):
    initial,mode=data;mass=2.
    assert initial_parent_mode_pump(initial,mode,0.)['total_pair_Hessian']==0.
    round_initial=solve_necks(16,(0.,))[0]
    _,modes=calculate_dirac_level(4,initial=round_initial)
    round_mode=next(m for m in modes if m['sign']>0)
    result=initial_parent_mode_pump(round_initial,round_mode,mass)
    assert result['lambda_ddot']==pytest.approx(-round_initial.parent.H**2*round_mode['eigenvalue'],abs=2e-15)
    assert result['projected_mass_ddot']==0.
    expected=-mass*round_mode['eigenvalue']/math.hypot(mass,round_mode['eigenvalue'])*round_initial.parent.H**2
    assert result['total_pair_Hessian']==pytest.approx(expected,abs=2e-15)


def test_actual_finite_geometry_and_radion_tend_to_the_initial_pump(data):
    initial,mode=data;mass=2.5958423951725464
    analytic=initial_parent_mode_pump(initial,mode,mass)
    with threadpool_limits(limits=1):
        system=SphericalValve(80);solution=system.evolve(.01)
        rows=[finite_parent_mode_pump_check(system,solution.sol(t).reshape(8,80),t,initial,mode,mass)
              for t in (.01,.005,.0025)]
    errors=[abs(r['pair_Hessian_finite']-analytic['total_pair_Hessian']) for r in rows]
    assert errors[1]<errors[0]/3.5 and errors[2]<errors[1]/3.5
    assert errors[-1]<6e-7
    for r in rows:
        assert r['mapping_residual']<5e-15
        assert r['angular_metric_identity_relative_error']<2e-12
        assert r['radial_metric_identity_relative_error']<2e-12


def test_exact_slice_map_preserves_round_spheres_and_small_actual_slice(data):
    system=SphericalValve(32,epsilon=0.)
    with threadpool_limits(limits=1):solution=system.evolve(.02)
    u=np.linspace(-.95,.95,51)
    mapped=slice_isotropic_map(system,solution.y[:,-1].reshape(8,32),u)
    np.testing.assert_allclose(mapped['original_u'],u,atol=3e-15,rtol=0)
    np.testing.assert_allclose(mapped['Omega'],math.cosh(system.H*.02)/system.H,atol=4e-13,rtol=0)


def test_fixed_initial_and_instantaneous_final_beta_have_different_leading_orders():
    lam=.7;mass=2.;ld=-.1;md=-.8;energy=math.hypot(lam,mass)
    tau1=np.array([[0,1],[1,0]],complex);tau3=np.diag([1.,-1.]).astype(complex)
    H0=lam*tau1+mass*tau3;Hdd=ld*tau1+md*tau3
    eig,vec=np.linalg.eigh(H0);negative=vec[:,0];positive=vec[:,1]
    M=complex(positive.conj()@Hdd@negative)
    t=.01
    solution=solve_ivp(lambda at,psi:-1j*(H0+.5*at*at*Hdd)@psi,(0.,t),negative,
                       method='DOP853',rtol=1e-13,atol=1e-15)
    final=solution.y[:,-1]
    beta_initial=positive.conj()@final
    _,final_vectors=np.linalg.eigh(H0+.5*t*t*Hdd)
    beta_final=final_vectors[:,1].conj()@final
    assert abs(beta_initial)/t**3==pytest.approx(abs(M)/6,rel=3e-4)
    assert abs(beta_final)/t**2==pytest.approx(abs(M)/(4*energy),rel=3e-4)
