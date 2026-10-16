"""Independent characteristic curves and radial metric constraints."""
import mpmath as mp
import numpy as np
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_wavefront import wavefront_rhs
from recursive_horizons.nsc_atlas_cap_radial_reference import reference_potential
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget


CONTROLS=ControlBudget(maximum_evaluations=50_000,cpu_seconds=15.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('wavefront stage controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def options(budget):
    return dict(orientation=1,wall_radius=.8,alpha=.02,mass=.003,lapse=.04,
                field=.12,clock_rate=1.3,budget=budget)


@pytest.mark.parametrize('sign',[-1,1])
def test_geometric_front_matches_independent_exact_curve_and_amplitude(sign,budget):
    with mp.workdps(70):
        z0,t,j0=map(mp.mpf,('.4','.03','.02'))
        radius=lambda u:mp.tanh(z0+sign*u)
        amplitude=lambda u:j0*mp.cosh(z0+sign*u)**2/mp.cosh(z0)**2
        r,j=float(radius(t)),float(amplitude(t))
        dr,dj=float(mp.diff(radius,t)),float(mp.diff(amplitude,t))
    result=wavefront_rhs(r,j,**(options(budget)|dict(orientation=sign,alpha=0.,mass=0.,lapse=0.,clock_rate=1.)))
    np.testing.assert_allclose([result['radius_rate'],result['jump_rate']],[dr,dj],rtol=3e-14,atol=1e-16)
    assert not result['numerical_error_enclosure']


def test_metric_derivative_cancellation_and_both_physical_jump_equations(budget):
    r,j,pi,phi=.37,-.04,.21,-.13
    opts=options(budget)
    row=wavefront_rhs(r,j,**opts)
    al,n,delta=opts['alpha'],opts['mass'],opts['lapse']
    P=float(reference_potential(np.array([opts['field']]),al,budget=budget)['Pdiff'][0])
    rho=r*r*(P/2+row['A']*(pi*pi+phi*phi)/4)
    delta_r=al*r*(pi*pi+phi*phi)/2
    A_r=-2*r-2*al*(rho/r-n/(r*r))
    independent=np.exp(delta)*(A_r+row['A']*delta_r)
    np.testing.assert_allclose(row['B_r'],independent,rtol=3e-14,atol=1e-16)
    v,tx,B=row['radius_rate'],opts['clock_rate'],row['B']
    dp,dh=row['Pi_jump'],row['Phi_jump']
    np.testing.assert_allclose([v*dp+tx*B*dh,v*dh+tx*B*dp],[0.,0.],atol=2e-17,rtol=0)


def test_clock_scaling_and_zero_amplitude(budget):
    opts=options(budget)
    base=wavefront_rhs(.3,.01,**opts)
    scaled=wavefront_rhs(.3,.01,**(opts|dict(clock_rate=2*opts['clock_rate'])))
    np.testing.assert_allclose([scaled[k] for k in ('radius_rate','jump_rate')],
                               [2*base[k] for k in ('radius_rate','jump_rate')],rtol=2e-15)
    quiet=wavefront_rhs(.3,0.,**opts)
    assert quiet['jump_rate']==quiet['Pi_jump']==quiet['Phi_jump']==0


def test_endpoint_bad_orientation_metric_and_budget_reject(budget):
    opts=options(budget)
    for radius in (0.,opts['wall_radius'],1.):
        with pytest.raises(ValueError): wavefront_rhs(radius,.1,**opts)
    for change in (dict(orientation=0),dict(orientation=True),dict(orientation=1+0j),
                   dict(orientation=np.array([1.])),dict(clock_rate=0.),
                   dict(mass=100.),dict(field=np.nan),dict(lapse=1000.)):
        with pytest.raises(ValueError): wavefront_rhs(.3,.1,**(opts|change))
    with pytest.raises(RuntimeError,match='budget'):
        wavefront_rhs(.3,.1,**(opts|dict(budget=Budget(maximum_evaluations=1,cpu_seconds=1.))))
