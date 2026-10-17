"""Independent constraint, wall projection and quadrature controls."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_connection_slice import wall_slice_data
from recursive_horizons.nsc_atlas_geometric_average import (
    KAPPA, clocked_scalar_average_sources, geometric_slice_average,
    periodic_spherical_average, scalar_average_sources, wall_slice_stress,
)


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
