"""Independent constraint, wall projection and quadrature controls."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_connection_slice import wall_slice_data
from recursive_horizons.nsc_atlas_geometric_average import (
    KAPPA, clocked_scalar_average_sources, geometric_slice_average,
    periodic_spherical_average, scalar_average_sources, wall_slice_stress,
    partitioned_scalar_average, partition_volume_rates, averaged_gravity_balance,
    expansion_source_inference,
)


def test_partition_composition_from_point_moments_and_unique_surface_source():
    # Weighted points independently supply raw moments; arbitrary regrouping
    # must preserve Q, the interface delta and the integrated constraint.
    theta = np.array([.3, .8, -.4, 1.2, .6, -.2])
    volume = np.array([1., 2., .5, 3., 1.5, .7])
    shear = np.array([.02, .01, .03, .04, .02, .01])
    energy = np.array([.1, .2, .3, .1, .2, .3])
    J = 2*theta**2/3-2*shear
    R = 2*KAPPA*energy-J
    raw = np.column_stack((volume, volume*theta, volume*theta**2,
                           volume*shear, volume*R, volume*energy, -volume*energy))
    names = ['a', 'b', 'c', 'd', 'e', 'f']
    shell_E, shell_S = .23, -.19
    shell = (2*KAPPA*shell_E, shell_E, shell_S)
    direct = partitioned_scalar_average(dict(zip(names, raw)), surface_integrals=shell)
    left = partitioned_scalar_average(dict(zip(names[:3], raw[:3])), surface_integrals=shell)
    right = partitioned_scalar_average(dict(zip(names[3:], raw[3:])))
    merged = partitioned_scalar_average(dict(left=left['moment_vector'], right=right['moment_vector']))
    expected = 2*(np.average(theta**2, weights=volume)
                  -np.average(theta, weights=volume)**2)/3-2*np.average(shear, weights=volume)
    assert direct['sources']['Q'] == pytest.approx(expected)
    assert abs(direct['sources']['hamiltonian_residual']) < 3e-14
    assert direct['between_region_Q'] > 0
    for key, value in direct['sources'].items():
        assert merged['sources'][key] == pytest.approx(value, abs=3e-14)
    assert abs(direct['Q_composition_residual']) < 2e-15
    # Counting the same wall again violates the independently known total E.
    assert direct['moment_vector'][5] == pytest.approx(volume @ energy+shell_E)


def test_moving_partition_reynolds_rate_and_weighted_expansion_derivative():
    import sympy as sp
    t = sp.symbols('t', real=True)
    V1, V2 = 1+t/3+t*t/7, 2+t/5+t*t/11
    H1 = sp.Rational(3, 10)+t/10
    H2 = (sp.diff(V1+V2, t)-3*V1*H1)/(3*V2)
    flux1 = sp.diff(V1, t)-3*V1*H1
    at = {t: sp.Rational(2, 7)}
    value = lambda expression: float(expression.subs(at))
    flux = value(flux1)
    row = partition_volume_rates(dict(a=value(V1), b=value(V2)),
        dict(a=value(H1), b=value(H2)), dict(a=flux, b=-flux))
    total_H = sp.diff(V1+V2, t)/(3*(V1+V2))
    regional_rate = value((V1*sp.diff(H1, t)+V2*sp.diff(H2, t))/(V1+V2))
    assert row['H_normal'] == pytest.approx(value(total_H))
    assert row['regional_volume_rates']['a'] == pytest.approx(value(sp.diff(V1, t)))
    assert row['volume_fraction_rates']['a'] == pytest.approx(value(sp.diff(V1/(V1+V2), t)))
    assert row['total_weight_H_rate'] == pytest.approx(row['direct_weight_H_rate'])
    assert regional_rate+row['total_weight_H_rate'] == pytest.approx(value(sp.diff(total_H, t)))
    assert abs(row['interface_weight_H_rate']) > .01


def test_time_balance_from_independent_friedmann_and_raychaudhuri_with_lapse():
    import sympy as sp
    t = sp.symbols('t', real=True)
    kappa = sp.symbols('kappa', positive=True)
    H, E, S, P = (sp.Function(name)(t) for name in ('H', 'E', 'S', 'P'))
    Q = 3*(sp.diff(H, t)+H*H)+kappa*(E+S)/2-P
    R = 2*kappa*E-6*H*H-Q
    trace, A = (Q-R+2*P)/4, 3*Q+R+2*P
    J = sp.diff(E, t)+3*H*E+H*S
    assert sp.simplify(sp.diff(trace, t)-sp.diff(A, t)/4-H*A+kappa*J) == 0
    substitutions = {H: 1+t/7, E: (1+t*t)/kappa, S: -(1+t)/kappa, P: t/5}
    evaluate = lambda expression: float(expression.subs(substitutions).doit().subs(
        {t: sp.Rational(2, 3), kappa: KAPPA}))
    row = averaged_gravity_balance(H=evaluate(H), Q=evaluate(Q), curvature=evaluate(R),
        lapse_acceleration=evaluate(P), Q_rate=evaluate(sp.diff(Q, t)),
        curvature_rate=evaluate(sp.diff(R, t)), lapse_acceleration_rate=evaluate(sp.diff(P, t)),
        mean_energy=evaluate(E), mean_stress_trace=evaluate(S), mean_energy_rate=evaluate(sp.diff(E, t)))
    inference = expansion_source_inference(H=evaluate(H), H_rate=evaluate(sp.diff(H, t)),
        mean_energy=evaluate(E), mean_stress_trace=evaluate(S))
    assert abs(row['integrability_residual']) < 2e-14
    assert abs(row['trace_balance_residual']) < 2e-14
    assert abs(row['geometric_continuity_residual']) < 2e-14
    assert row['Lambda_trace'] == pytest.approx(inference['Lambda_trace'])
    assert row['vacuum_condition'] == pytest.approx(inference['vacuum_condition'])


def test_nonvacuum_radiation_mixture_can_have_a_constant_vacuum_trace():
    # Manufactured averaged Lambda+radiation history: H=1, Lambda=.6,
    # Hdot=-2(H^2-Lambda/3), explicit source zero. This is a control of the
    # balance theorem, not a derived NSC metric or a fitted vacuum value.
    H, lam = 1., .6
    H_rate = -2*(H*H-lam/3)
    A = 4*lam-12*H*H
    A_rate = -24*H*H_rate
    row = averaged_gravity_balance(H=H, Q=lam+A/4, curvature=A/4-3*lam,
        Q_rate=A_rate/4, curvature_rate=A_rate/4,
        mean_energy=0., mean_stress_trace=0., mean_energy_rate=0.)
    assert row['vacuum_condition'] < 0
    assert row['Lambda_trace'] == pytest.approx(lam)
    assert row['Lambda_trace_rate'] == pytest.approx(0., abs=1e-14)
    assert abs(row['integrability_residual']) < 1e-14
    assert row['kappa_geometric_density']+row['kappa_geometric_pressure'] > 0
    assert row['kappa_geometric_density'] == pytest.approx(3*H*H)


def test_homogeneous_flat_cap_has_no_geometric_correction():
    H, volume = 1.3, 2.7
    energy = 3*H**2*volume/KAPPA
    row = scalar_average_sources(volume, 3*H*volume, 9*H**2*volume,
                                 0., 0., energy, -3*energy)
    assert row['H_normal'] == pytest.approx(H)
    assert row['Q'] == pytest.approx(0., abs=3e-15)
    assert row['kappa_geometric_density'] == pytest.approx(0., abs=3e-15)
    assert row['kappa_geometric_pressure'] == pytest.approx(0., abs=3e-15)
    assert row['total_trace_projection'] == pytest.approx(3*H**2)
    assert row['three_volume_acceleration'] == pytest.approx(3*H**2)


def test_expansion_variance_is_not_automatically_vacuum():
    # Equal volumes with expansion +/-3H have Q=6H^2 and H_D=0.
    row = scalar_average_sources(2., 0., 18., 0., 0., 6/KAPPA, 0.)
    assert row['Q'] == pytest.approx(6.)
    assert row['hamiltonian_residual'] == pytest.approx(0.)
    assert row['kappa_geometric_density'] == pytest.approx(-3.)
    assert row['kappa_geometric_pressure'] == pytest.approx(-3.)
    assert row['geometric_vacuum_condition'] == pytest.approx(18.)


def test_surface_stresses_from_four_tensor_and_restricted_delta():
    wall = wall_slice_data()
    metric = np.diag([-1., 1., 1., 1.])
    normal = np.array([wall['gamma_v'], wall['gamma'], 0., 0.])
    normal_lower = metric@normal
    tensor = -wall['parameters'].sigma*(metric-np.outer(normal_lower, normal_lower))/wall['gamma']
    row = wall_slice_stress(wall)
    assert tensor[0, 0] == pytest.approx(row['energy'])
    assert tensor[1, 1] == pytest.approx(row['radial_stress'])
    assert tensor[2, 2] == tensor[3, 3] == pytest.approx(row['angular_stress'])
    assert np.trace(tensor[1:, 1:]) == pytest.approx(row['spatial_stress_trace'])
    assert -tensor[0, 1] == pytest.approx(wall['momentum_surface'])


def test_averaged_constraint_requires_the_wall_curvature_distribution():
    record = geometric_slice_average()
    row = record['normalized_sources']
    delta = record['integrals_in_ell_units']['curvature_wall_delta']
    assert delta > 100
    assert abs(row['hamiltonian_residual']) < 2e-12
    omitted_delta_residual = row['hamiltonian_residual']-delta/(2*record['volume_over_ell3'])
    assert abs(omitted_delta_residual) > .1
    assert abs(row['scalar_fitting_acceleration_residual']) < 2e-14
    assert max(abs(x) for x in record['checks'].values()) < 2e-10


def test_single_domain_quadratures_agree_and_unit_scaling_cancels():
    reference = geometric_slice_average(quadrature='adaptive')
    for order in ('gauss32', 'gauss64'):
        row = geometric_slice_average(quadrature=order)
        assert row['volume_over_ell3'] == pytest.approx(reference['volume_over_ell3'], rel=2e-11)
        for key, value in row['normalized_sources'].items():
            # The fixed 32-point control resolves the flat taper less closely;
            # the recorded differences, not this acceptance limit, are reported.
            tolerance = 1e-9 if order == 'gauss32' else 2e-12
            assert value == pytest.approx(reference['normalized_sources'][key], abs=tolerance, rel=tolerance)
    scaled = geometric_slice_average(length_scale=3.7)
    for key, value in scaled['normalized_sources'].items():
        assert value == pytest.approx(reference['normalized_sources'][key], abs=5e-11, rel=2e-10)


def test_fixed_domain_reports_scalar_outputs_without_vacuum_identification():
    record = geometric_slice_average()
    row = record['normalized_sources']
    assert record['domain']['Rcut_right_over_ell'] == 4.
    assert 0 < record['cap_volume_fraction'] < 1
    assert abs(row['geometric_vacuum_condition']) > .01
    assert record['spatial_evaluations'] <= 20000
    with pytest.raises(ValueError, match='quadrature'):
        geometric_slice_average(quadrature='unbounded_scan')


def test_reject_inconsistent_or_nonfinite_moments():
    for values in ((0., 0., 0., 0., 0., 0., 0.),
                   (1., math.nan, 1., 0., 0., 0., 0.),
                   (1., 2., 1., 0., 0., 0., 0.)):
        with pytest.raises(ValueError):
            scalar_average_sources(*values)


def test_homogeneous_clock_change_matches_direct_volume_derivatives():
    # a = exp(H*t_proper), dt_proper/dtau_ref = nu. Differentiate a directly.
    H, nu, nu_rate, volume = .7, 1.4, -.3, 2.3
    energy = 3*H**2/KAPPA
    row = clocked_scalar_average_sources(
        volume, lapse_theta_integral=volume*nu*3*H,
        lapse_squared_theta_squared_integral=volume*(nu*3*H)**2,
        lapse_squared_shear_squared_integral=0.,
        lapse_squared_curvature_integral=0.,
        lapse_squared_energy_integral=volume*nu**2*energy,
        lapse_squared_stress_trace_integral=-3*volume*nu**2*energy,
        lapse_rate_theta_integral=volume*nu_rate*3*H,
        lapse_laplacian_integral=0.)
    assert row['H_observer'] == pytest.approx(nu*H)
    assert row['three_volume_acceleration'] == pytest.approx(
        3*(nu_rate*H+(nu*H)**2))
    assert abs(row['hamiltonian_residual']) < 3e-15
    assert abs(row['scalar_fitting_acceleration_residual']) < 3e-15


def test_static_periodic_metric_lapse_cancels_its_actual_ricci_source():
    # ds^2 = -N(x)^2 dt^2 + dx^2 + dy^2 + dz^2, N = 1+eps*cos(x).
    # Direct connection contraction: R_nn=N_xx/N, E=0, kappa*S=2*N_xx/N.
    # Its spatial volume is static despite nonzero spatial lapse variation.
    eps, volume = .4, 2*math.pi
    nodes, weights = np.polynomial.legendre.leggauss(64)
    x, weights = math.pi*(nodes+1), math.pi*weights
    lapse, second = 1+eps*np.cos(x), -eps*np.cos(x)
    laplacian = float(weights@(lapse*second))
    stress_trace = float(weights@(2*lapse*second/KAPPA))
    row = clocked_scalar_average_sources(
        volume, lapse_theta_integral=0.,
        lapse_squared_theta_squared_integral=0.,
        lapse_squared_shear_squared_integral=0.,
        lapse_squared_curvature_integral=0.,
        lapse_squared_energy_integral=0.,
        lapse_squared_stress_trace_integral=stress_trace,
        lapse_rate_theta_integral=0., lapse_laplacian_integral=laplacian)
    assert laplacian/volume == pytest.approx(-eps**2/2, abs=2e-15)
    assert row['three_volume_acceleration'] == pytest.approx(0., abs=2e-15)
    # Doubling the lapse-Laplacian would incorrectly accelerate static volume.
    assert abs(row['three_volume_acceleration']+laplacian/volume) > .07
    total_pressure = row['kappa_mean_pressure']+row['kappa_geometric_pressure']
    assert total_pressure == pytest.approx(0., abs=2e-15)
    assert row['total_trace_projection'] == pytest.approx(0., abs=2e-15)


def _constant(value, count):
    return np.full(count, value)


def test_static_product_circle_has_round_sphere_curvature():
    radius, scale, lapse, length, count = 2.5, 1.7, 1.3, 2*math.pi, 16
    row = periodic_spherical_average(
        _constant(radius, count), _constant(scale, count), _constant(lapse, count),
        _constant(0., count), _constant(0., count), _constant(0., count),
        _constant(0., count), _constant(0., count), length, 3)
    assert row['H_observer'] == pytest.approx(0., abs=1e-12)
    assert row['Q'] == pytest.approx(0., abs=1e-12)
    assert row['lapse_acceleration'] == pytest.approx(0., abs=1e-12)
    assert row['mean_intrinsic_curvature'] == pytest.approx(2/radius**2)
    assert row['geometric_trace_projection'] == pytest.approx(-1/(2*radius**2))
    assert row['geometric_vacuum_condition'] == pytest.approx(2/radius**2)
    assert row['three_volume_acceleration'] == pytest.approx(0., abs=1e-12)
    assert row['curvature_ibp_residual'] == pytest.approx(0., abs=1e-10)
    assert row['lapse_ibp_residual'] == pytest.approx(0., abs=1e-12)


def test_warped_static_sphere_matches_integrated_curvature_identity():
    radius, amplitude, scale, length, count = 2., .2, 1.5, 4., 64
    nodes = np.arange(count)*length/count
    areal = radius*(1+amplitude*np.cos(2*math.pi*nodes/length))
    row = periodic_spherical_average(
        areal, _constant(scale, count), _constant(1., count),
        _constant(0., count), _constant(0., count), _constant(0., count),
        _constant(0., count), _constant(0., count), length, 0)
    expected = 8*math.pi*(scale*length+(radius*amplitude)**2*(2*math.pi/length)**2*length/(2*scale))
    assert row['curvature_proper_integral'] == pytest.approx(expected, rel=1e-12)
    assert row['curvature_ibp_residual'] == pytest.approx(0., abs=1e-9)
    assert row['three_volume_acceleration'] == pytest.approx(0., abs=1e-10)
    varying = 1+.25*np.sin(2*math.pi*nodes/length)
    tilted = periodic_spherical_average(
        areal, _constant(scale, count), varying,
        _constant(0., count), _constant(0., count), _constant(0., count),
        _constant(0., count), _constant(0., count), length, 5)
    assert tilted['lapse_ibp_residual'] == pytest.approx(0., abs=1e-9)
    assert tilted['three_volume_acceleration'] == pytest.approx(
        tilted['direct_volume_acceleration'], abs=1e-9)


def test_homogeneous_anisotropic_volume_and_length_scaling():
    count, length, index = 24, 3., 7
    values = {'radius': 2., 'scale': 3., 'lapse': 1.5,
              'radius_rate': .2, 'scale_rate': -.1, 'lapse_rate': .05,
              'radius_acceleration': .01, 'scale_acceleration': -.02}
    def snapshot(scale):
        return periodic_spherical_average(
            _constant(scale*values['radius'], count), _constant(scale*values['scale'], count),
            _constant(scale*values['lapse'], count), _constant(scale*values['radius_rate'], count),
            _constant(scale*values['scale_rate'], count), _constant(scale*values['lapse_rate'], count),
            _constant(scale*values['radius_acceleration'], count),
            _constant(scale*values['scale_acceleration'], count), length, index)
    row = snapshot(1.)
    volume = 4*math.pi*values['scale']*values['radius']**2*length
    volume_rate = 4*math.pi*(values['scale_rate']*values['radius']**2
                             +2*values['scale']*values['radius']*values['radius_rate'])*length
    volume_acceleration = 4*math.pi*(
        values['scale_acceleration']*values['radius']**2
        +4*values['scale_rate']*values['radius']*values['radius_rate']
        +2*values['scale']*(values['radius_rate']**2
                            +values['radius']*values['radius_acceleration']))*length
    lapse = values['lapse']
    direct = (volume_acceleration/(volume*lapse**2)
              -volume_rate*values['lapse_rate']/(volume*lapse**3)
              -(2/3)*(volume_rate/(volume*lapse))**2)
    expansion = (values['scale_rate']/values['scale']
                 +2*values['radius_rate']/values['radius'])/lapse
    assert row['H_observer'] == pytest.approx(expansion/3)
    assert row['three_volume_acceleration'] == pytest.approx(direct, rel=1e-12)
    assert row['direct_volume_acceleration'] == pytest.approx(direct, rel=1e-12)
    scaled = snapshot(2.5)
    assert scaled['H_observer'] == pytest.approx(row['H_observer']/2.5, rel=1e-12)
    for key in ('Q', 'mean_intrinsic_curvature', 'lapse_acceleration',
                'kappa_mean_energy', 'kappa_geometric_density',
                'three_volume_acceleration', 'geometric_trace_projection',
                'geometric_vacuum_condition'):
        assert scaled[key] == pytest.approx(row[key]/2.5**2, rel=1e-11)


def test_explicit_defaults_preserve_slice_average():
    default = geometric_slice_average()
    explicit = geometric_slice_average(mass_over_ell=.4, taper_stop_over_ell=1.)
    assert explicit['volume_over_ell3'] == default['volume_over_ell3']
    assert explicit['integrals_in_ell_units'] == default['integrals_in_ell_units']
    assert explicit['normalized_sources'] == default['normalized_sources']
    assert explicit['checks'] == default['checks']
    assert explicit['domain']['M_over_ell'] == .4
    assert 'taper_stop_over_ell' not in explicit['domain']
    assert explicit['domain']['H_ell'] == 1. and explicit['domain']['k_ell'] == 2.
    assert explicit['domain']['Rwall_over_ell'] == 1.5
    assert explicit['domain']['Rcut_right_over_ell'] == 4.


def test_unequal_body_moments_stay_finite_and_scale():
    row = geometric_slice_average(mass_over_ell=.6, taper_stop_over_ell=1.3)
    moments = row['integrals_in_ell_units']
    J = 2*moments['theta_squared']/3-2*moments['shear_squared']
    cap_volume = 4*math.pi*1.5**3/3
    assert row['volume_over_ell3'] > cap_volume
    assert math.isfinite(J)
    assert moments['energy_cap'] == pytest.approx(1.6875, abs=1e-12)
    assert moments['stress_trace_cap'] == pytest.approx(-3*1.6875, abs=1e-12)
    assert moments['curvature_wall_delta'] > 0
    assert moments['curvature_bulk']+moments['curvature_wall_delta'] == pytest.approx(
        moments['curvature_total'], abs=1e-12)
    assert abs(row['normalized_sources']['hamiltonian_residual']) < 1e-9
    assert abs(row['checks']['bridge_volume_ell3_difference']) < 1e-8
    assert row['domain']['M_over_ell'] == pytest.approx(.6)
    assert row['domain']['taper_stop_over_ell'] == pytest.approx(1.3)
    assert row['domain']['H_ell'] == 1. and row['domain']['k_ell'] == 2.
    assert row['domain']['Rwall_over_ell'] == 1.5
    assert row['domain']['Rcut_right_over_ell'] == 4.
    assert row['spatial_evaluations'] <= 20000
    assert 'four-tensor' in row['scope']
    scaled = geometric_slice_average(length_scale=3.7, mass_over_ell=.6, taper_stop_over_ell=1.3)
    assert scaled['volume_over_ell3'] == pytest.approx(row['volume_over_ell3'], rel=2e-10, abs=5e-11)
    for key, value in scaled['normalized_sources'].items():
        assert value == pytest.approx(row['normalized_sources'][key], abs=5e-11, rel=2e-10)
    for key in ('theta', 'theta_squared', 'shear_squared', 'curvature_total',
                'energy_cap', 'energy_wall', 'stress_trace_wall'):
        assert scaled['integrals_in_ell_units'][key] == pytest.approx(
            moments[key], rel=2e-10, abs=1e-8)
    with pytest.raises(ValueError):
        geometric_slice_average(mass_over_ell=.6, taper_stop_over_ell=1.2)
