"""Independent local KG/energy controls; no parent caches or evolution."""
import math
from types import SimpleNamespace

import numpy as np
import pytest
import sympy as sp

from recursive_horizons.nsc_atlas_null_geometry import first_connection
from scripts.parent_release_response import (
    causal_release_response, instantaneous_release, potential_hessian,
)
from scripts.parent_release_region import (
    analytic_geometry_controls, nodal_radius_control, region_quadrature,
    source_region, normal_boost_flux,
)


K_S = 84.
ETA = np.diag([-1., 1., 1., 1.])


def _potential(s, chi):
    """Independent analytic derivatives, with canonical S=sqrt(84)*s."""
    e12, e14, e24 = np.exp(-12*s), np.exp(-14*s), np.exp(-24*s)
    mu2 = 400*e12
    V = 7*e12-7.5*e14+25/32*e24
    Vs = -84*e12+105*e14-75/4*e24
    return V, Vs, mu2


def _curved_control(H=.2, lapse_gradient=.3):
    """Jets at t=x=0 of ds2=-exp(2gx)dt2+exp(2Ht)dx_i dx_i.

    Scalar Taylor jets obey covariant radion KG at this event.  The
    metric/control state is analytic test data, not scientific preparation.
    """
    g = lapse_gradient
    s, chi = .025, .011
    v, c, d = .017, -.023, .031
    V, Vs, mu2 = _potential(s, chi)
    force = (-Vs+6*mu2*chi**2)/K_S
    s_tt = -3*H*v+d+g*c+force
    metric = ETA.copy()
    dtmetric = np.diag([0., 2*H, 2*H, 2*H])
    dxmetric = np.zeros((3, 4, 4)); dxmetric[0, 0, 0] = -2*g
    Gamma = first_connection(metric, dtmetric, dxmetric, frame_scale=0.)['Gamma']
    ds = np.array([v, c, 0., 0.])
    partial2 = np.diag([s_tt, d, 0., 0.])
    Hess = partial2-np.einsum('cab,c->ab', Gamma, ds)
    return dict(metric=metric, Gamma=Gamma, ds=ds, Hess=Hess, s=s, chi=chi,
                v=v, c=c, d=d, H=H, g=g, s_tt=s_tt, force=force,
                V=V, Vs=Vs, mu2=mu2)


@pytest.mark.parametrize('H,g', [(.2, .3), (-.1, .4), (.15, 0.), (0., 0.)])
def test_covariant_KG_has_correct_expansion_gradient_and_lapse_signs(H, g):
    z = _curved_control(H, g)
    h = z['Hess']; v, c, d = z['v'], z['c'], z['d']
    assert h == pytest.approx(h.T, abs=2e-15)
    assert h[0, 0] == pytest.approx(z['s_tt']-g*c, abs=2e-15)
    assert h[0, 1] == pytest.approx(-g*v-H*c, abs=2e-15)
    assert h[1, 1] == pytest.approx(d-H*v, abs=2e-15)
    assert h[2, 2] == pytest.approx(-H*v, abs=2e-15)
    box = float(np.einsum('ab,ab', ETA, h))
    assert box == pytest.approx(-z['force'], abs=3e-15)
    # n[n(s)] includes the lapse acceleration; Hess_nn does not.
    assert z['s_tt'] == pytest.approx(-3*H*v+d+g*c+z['force'], abs=2e-15)


@pytest.mark.parametrize('speed', [0., .4, -.3])
def test_tilted_geodesic_radion_energy_uses_observer_projected_hessian(speed):
    z = _curved_control()
    gamma = 1/math.sqrt(1-speed**2)
    u = gamma*np.array([1., speed, 0., 0.])
    spatial = ETA+np.outer(u, u)
    rate = float(u@z['ds'])
    acceleration = float(u@z['Hess']@u)
    spatial_hessian = float(np.einsum('ab,ab', spatial, z['Hess']))
    assert acceleration == pytest.approx(spatial_hessian+z['force'], abs=2e-15)
    direct = K_S*rate*acceleration+z['Vs']*rate
    resolved = K_S*rate*spatial_hessian+6*z['mu2']*z['chi']**2*rate
    assert direct == pytest.approx(resolved, abs=2e-15)


def test_geodesic_initially_normal_is_not_accelerated_slice_normal():
    z = _curved_control()
    # Both curves have the same initial unit tangent, but different
    # second derivatives because the slice normals have a^x=g.
    geodesic_second = z['Hess'][0, 0]
    normal_second = z['s_tt']
    assert normal_second-geodesic_second == pytest.approx(z['g']*z['c'], abs=2e-15)
    assert abs(normal_second-geodesic_second) > 1e-3


def test_Minkowski_radion_carrier_energy_exchange_cancels_exactly():
    s, chi, v, pi = .027, -.019, .023, .031
    V, Vs, mu2 = _potential(s, chi)
    sdd = (-Vs+6*mu2*chi**2)/K_S
    chidd = -mu2*chi
    radion_change = K_S*v*sdd+Vs*v
    carrier_change = pi*chidd+mu2*chi*pi-6*mu2*chi**2*v
    assert radion_change == pytest.approx(6*mu2*chi**2*v, abs=2e-15)
    assert radion_change+carrier_change == pytest.approx(0., abs=2e-15)


def test_canonical_potential_hessian_matches_independent_symbolic_action():
    S, chi = sp.symbols('S chi', real=True)
    s = S/sp.sqrt(84)
    action_potential = (7*sp.exp(-12*s)-sp.Rational(15, 2)*sp.exp(-14*s)
        + sp.Rational(25, 32)*sp.exp(-24*s)+200*sp.exp(-12*s)*chi**2)
    symbolic = sp.lambdify((S, chi), sp.hessian(action_potential, (S, chi)), 'numpy')
    volumes = np.array([-.035, .025, .061])
    carriers = np.array([0., .011, -.019])
    expected = np.array([symbolic(math.sqrt(K_S)*v, c) for v, c in zip(volumes, carriers)])
    actual = potential_hessian(volumes, carriers)
    assert actual.shape == (3, 2, 2)
    assert actual == pytest.approx(expected, rel=3e-14, abs=2e-12)


@pytest.mark.parametrize('speed', [0., .4, -.3])
def test_release_helper_closes_actual_curved_metric_KG_with_lapse_and_tilt(speed):
    z = _curved_control()
    gamma = 1/math.sqrt(1-speed**2)
    u = gamma*np.array([1., speed, 0., 0.])
    Huu = float(u@z['Hess']@u)
    rootK = math.sqrt(K_S)
    sstar = -.035646782552940987
    # Repeated local jets exercise the pointwise API, not time
    # reconstruction: no propagated Taylor history is asserted here.
    repeat = lambda x: np.full(2, x)
    values = dict(tau=np.array([0., .01]),
        S=repeat(rootK*(z['s']-sstar)), chi=repeat(z['chi']),
        uS=repeat(rootK*float(u@z['ds'])), uChi=repeat(0.),
        nS=repeat(rootK*z['v']), nChi=repeat(0.),
        HessS_uu=repeat(rootK*Huu), HessS_nn=repeat(rootK*z['Hess'][0, 0]),
        spatial_lapS=repeat(rootK*z['d']),
        lapse_accel_dot_DS=repeat(rootK*z['g']*z['c']), traceK=repeat(-3*z['H']))
    output = causal_release_response(values, sstar=sstar)
    summary, series = output['summary'], output['series']
    assert summary['maximum_abs_kg_normal_residual'] < 1e-13
    assert summary['maximum_abs_kg_observer_residual'] < 1e-13
    labels = series['force_labels']
    forces = np.array(series['measured_forces'])
    assert forces.sum(axis=1) == pytest.approx(rootK*Huu, abs=1e-13)
    assert forces[:, labels.index('lapse_acceleration')] == pytest.approx(
        rootK*z['g']*z['c'], abs=1e-14)
    expected_conversion = rootK*(Huu-z['Hess'][0, 0]-z['g']*z['c'])
    assert forces[:, labels.index('observer_conversion')] == pytest.approx(expected_conversion, abs=1e-14)
    if speed == 0:
        assert forces[:, labels.index('observer_conversion')] == pytest.approx(
            -forces[:, labels.index('lapse_acceleration')], abs=1e-14)


@pytest.mark.parametrize('chi,pi', [(0., 0.), (.003, .001)])
def test_conditional_barriers_are_stationary_points_of_the_declared_comparisons(chi, pi):
    row = instantaneous_release(-.035, chi, pi)
    fixed = row['fixed_chi']
    assert fixed['extrema_exist']
    for kind, sign in (('minimum_s', 1), ('top_s', -1)):
        s = fixed[kind]
        _, Vs, mu2 = _potential(s, chi)
        assert Vs-6*mu2*chi**2 == pytest.approx(0., abs=2e-12)
        assert np.sign(potential_hessian(s, chi)[0, 0]) == sign
    adiabatic = row['fixed_action']
    assert adiabatic['extrema_exist']
    for key in ('minimum_s', 'top_s'):
        s = adiabatic[key]
        _, Vs, _ = _potential(s, chi)
        assert Vs-6*adiabatic['loading']*math.exp(-6*s) == pytest.approx(0., abs=2e-12)
    assert row['scope'].startswith('Conditional')


@pytest.mark.parametrize('z', [0., 7., 1.e6])
def test_source_region_radius_has_feasible_nodal_contact_and_global_support_bound(z):
    # Analytic source-geometry control, not new physical initial data.
    parent = SimpleNamespace(H=1., mass_squared=z, kappa=1.)
    result = source_region(parent)
    anchor = np.array(result['sampler_arguments']['anchor'])
    contact = np.array(result['closest_node_points'])
    a, b, q = result['carrier_a'], result['carrier_b'], result['q']
    radius = result['sampler_arguments']['radius']
    assert contact@anchor == pytest.approx(math.cos(radius), abs=1e-14)
    assert np.sum(contact*contact, axis=1) == pytest.approx(1., abs=1e-14)
    assert a*contact[:, 0]+b*contact[:, 1]*contact[:, 3] == pytest.approx(0., abs=1e-14)
    # Independently parameterize the complete reduced nodal interval.
    # Nonzero X3 only reduces the support function since anchor.X3=0.
    xmax = math.sqrt(1+q*q)-q
    x = np.linspace(-xmax, xmax, 1001)
    p = np.sqrt(np.maximum(0., (1-x*x-2*q*x)/2))
    m = np.sqrt(np.maximum(0., (1-x*x+2*q*x)/2))
    nodes = np.column_stack([x, (p+m)/math.sqrt(2), np.zeros_like(x), (p-m)/math.sqrt(2)])
    assert np.max(nodes@anchor) <= math.cos(radius)+1e-14
    check = nodal_radius_control(parent)
    assert check['root_error'] < 1e-12
    assert check['maximum_dot_error'] < 1e-12


def test_nodal_contact_formula_remains_feasible_at_both_geometric_limits():
    for q in (0., .25, .75, .999999, 1.):
        c = math.sqrt((1-q*q)/2)
        x = (math.sqrt(2)-1)*q
        m = math.sqrt(c*c+2*q*x)
        contact = np.array([x, (c+m)/math.sqrt(2), 0., (c-m)/math.sqrt(2)])
        anchor = np.array([q, c, 0., c])
        assert contact@contact == pytest.approx(1., abs=1e-14)
        assert q*contact[0]+contact[1]*contact[3] == pytest.approx(0., abs=1e-14)
        assert anchor@contact == pytest.approx((1+q*q)/math.sqrt(2)-q*q, abs=1e-14)


def test_round_ball_quadrature_and_boundary_geometry_are_anchor_independent():
    radius, scale = .8522887164833749, .2714056154003984
    for anchor in ([1., 0., 0., 0.], [.5, .5, .5, .5]):
        rule = region_quadrature(anchor, radius, order=8)
        assert np.sum(rule['volume_points']**2, axis=1) == pytest.approx(1., abs=1e-14)
        assert rule['boundary_points']@np.asarray(anchor) == pytest.approx(math.cos(radius), abs=1e-14)
        assert sum(rule['volume_weights']) == pytest.approx(
            2*math.pi*(radius-math.sin(2*radius)/2), abs=1e-13)
        controls = analytic_geometry_controls(anchor, radius, frame_scale=scale, order=8)
        for key in ('round_volume_error', 'round_area_error', 'flat_volume_error', 'flat_area_error'):
            assert controls[key] < 2e-11


@pytest.mark.parametrize('speed', [-.7, 0., .6])
def test_tube_boost_flux_matches_direct_scalar_current_contractions(speed):
    pi = np.array([.4, -.2])
    normal_derivative = np.array([-.3, .5])
    W = .73
    rho = .5*(pi@pi+normal_derivative@normal_derivative)+W
    pr = .5*(pi@pi+normal_derivative@normal_derivative)-W
    j = -float(pi@normal_derivative)
    gamma = 1/math.sqrt(1-speed**2)
    Uphi = gamma*(pi+speed*normal_derivative)
    nBphi = gamma*(normal_derivative+speed*pi)
    direct = -float(Uphi@nBphi)
    assert normal_boost_flux(rho, pr, j, speed) == pytest.approx(direct, abs=2e-15)


@pytest.mark.parametrize('speed', [-.8, .5])
def test_pure_vacuum_has_no_tube_material_flux_despite_shift_energy_transport(speed):
    W, lapse = .23, 1.3
    actual_tube = normal_boost_flux(W, -W, 0., speed)
    coordinate_shift_transport = -lapse*speed*W
    assert actual_tube == pytest.approx(0., abs=1e-15)
    assert abs(coordinate_shift_transport) > .1
    # In contrast, kinetic stress with zero Eulerian j has nonzero
    # observer-boosted flux. The two flux definitions are not interchangeable.
    assert normal_boost_flux(.5, .5, 0., speed) == pytest.approx(
        -speed/(1-speed**2), abs=2e-15)
