"""Independent geometry and correlated-error controls, without observed data."""
import math

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from recursive_horizons.nsc_atlas_cosmological_distance import (
    contrast_weights, distance_contrast, project_distance_contrast,
    quiet_cap_distance_shape, quiet_cap_domain, radial_conformal_lightcone,
)


def example():
    z = np.array([.03, .05, .17, .2])
    hel = z+.001
    mag = 5*np.log10(z*(1+hel))+np.array([17., 17., 17.3, 17.3])
    return z, hel, mag


def test_quiet_distance_from_null_integral_and_finite_wall_domain():
    for z in (.025, .075, .15, .25):
        emission_time = -math.log1p(z)
        comoving, _ = quad(lambda time: math.exp(-time), emission_time, 0.)
        assert float(quiet_cap_distance_shape(z, z)) == pytest.approx((1+z)*comoving, abs=1e-14)
        assert math.exp(emission_time)*comoving == pytest.approx(z/(1+z), abs=1e-14)
    domain = quiet_cap_domain()
    assert domain['emission_lookback_upper'] < domain['wall_horizon_crossing_lookback_lower']
    assert domain['wall_speed_squared_upper'] == 225/64-1+3/10+1/25 < 4


def test_analytic_independent_errors_and_common_covariance_mode():
    z, hel, mag = example()
    for common_variance in (0., .09):
        C = .04*np.eye(4)+common_variance*np.ones((4, 4))
        indices, weights, result = distance_contrast(z, hel, mag, C, np.zeros(4))
        np.testing.assert_array_equal(indices, np.arange(4))
        np.testing.assert_allclose(weights, [-.5, -.5, .5, .5], atol=1e-14)
        assert result['contrast_magnitudes'] == pytest.approx(.3, abs=1e-13)
        assert result['standard_uncertainty_magnitudes'] == pytest.approx(.2, abs=1e-14)
        assert result['exponential_contrast'] == pytest.approx(10**(.3/5))


def test_calibrator_exclusion_row_permutation_and_zero_point_cancellation():
    z, hel, mag = example()
    z = np.r_[z, .04]; hel = np.r_[hel, .041]; mag = np.r_[mag, 999.]
    flags = np.array([0, 0, 0, 0, 1])
    L = np.diag([.1, .2, .3, .15, .4])+np.tril(np.ones((5, 5)), -1)*.01
    C = L @ L.T
    indices, weights, reference = distance_contrast(z, hel, mag, C, flags)
    np.testing.assert_array_equal(indices, np.arange(4))
    order = np.array([3, 4, 0, 2, 1])
    _, _, permuted = distance_contrast(z[order], hel[order], mag[order]+100., C[np.ix_(order, order)], flags[order])
    assert permuted['contrast_magnitudes'] == pytest.approx(reference['contrast_magnitudes'], abs=1e-12)
    assert permuted['standard_uncertainty_magnitudes'] == pytest.approx(reference['standard_uncertainty_magnitudes'], abs=1e-13)
    for scale in (.01, 1., 100.):
        predicted = project_distance_contrast(weights, z[indices], hel[indices], scale*quiet_cap_distance_shape(z[indices], hel[indices]))
        assert predicted == pytest.approx(0., abs=1e-12)


def test_forward_shape_uses_each_redshift_and_same_weights():
    z, hel, _ = example()
    C = np.array([[.1, .02, .01, 0], [.02, .2, 0, .03], [.01, 0, .15, .04], [0, .03, .04, .12]])
    weights, _, checks = contrast_weights(C, [0, 0, 1, 1])
    shift = .3*z*z
    distances = quiet_cap_distance_shape(z, hel)*10**(shift/5)
    assert project_distance_contrast(weights, z, hel, distances) == pytest.approx(float(weights @ shift), abs=1e-13)
    assert checks['KKT_weight_maximum_difference'] < 1e-13
    assert abs(weights.sum()) < 1e-14


def _flat_annulus(length):
    def metric(_time, position):
        radius = length*(3.+position)
        return {'radius': radius, 'Q': 1./(3.+position)}
    return metric


def test_minkowski_annulus_is_euclidean_and_scales_by_length():
    nodes = np.linspace(-2., 2., 5)
    lookbacks = np.array([.2, .5])
    for direction in (1, -1):
        bare = radial_conformal_lightcone(
            lambda time, position: {'radius': 3.+position, 'Q': 1./(3.+position)},
            nodes, reception_time=0., reception_x=0., lookbacks=lookbacks,
            direction=direction)
        scaled = radial_conformal_lightcone(
            _flat_annulus(2.5), nodes, reception_time=0., reception_x=0.,
            lookbacks=lookbacks, direction=direction, quadrature_order=64)
        for row, lookback, scaled_row in zip(bare, lookbacks, scaled):
            assert row['redshift'] == pytest.approx(0., abs=1e-14)
            assert row['angular_diameter_distance'] == pytest.approx(lookback, abs=1e-12)
            assert row['luminosity_distance'] == pytest.approx(lookback, abs=1e-12)
            assert row['observer_lapse'] == pytest.approx(1.)
            assert scaled_row['redshift'] == pytest.approx(0., abs=1e-14)
            assert scaled_row['luminosity_distance'] == pytest.approx(2.5*lookback, abs=1e-11)


def test_off_centre_de_sitter_recovers_quiet_cap_shape():
    hubble, offset = .7, 8.
    nodes = np.linspace(-1/hubble-2., -1/hubble+.2, 6)

    def metric(time, position):
        scale = -1./(hubble*time)
        return {'radius': scale*(offset+position), 'Q': 1./(offset+position)}

    lookbacks = np.array([.15, .4, .9])
    for direction in (1, -1):
        rows = radial_conformal_lightcone(
            metric, nodes, reception_time=-1./hubble, reception_x=.2,
            lookbacks=lookbacks, direction=direction)
        for row, lookback in zip(rows, lookbacks):
            assert row['redshift'] == pytest.approx(hubble*lookback, rel=1e-12)
            assert row['luminosity_distance'] == pytest.approx(
                lookback*(1.+hubble*lookback), rel=1e-12)
            assert row['angular_diameter_distance'] == pytest.approx(
                lookback/(1.+hubble*lookback), rel=1e-12)


def test_jacobi_ode_matches_the_area_integral():
    nodes = np.linspace(-1., .4, 4)
    lookback, direction = .35, 1

    def metric(time, position):
        return {'radius': 2.+.15*time+.1*position,
                'Q': .8+.05*time-.03*position}

    row = radial_conformal_lightcone(
        metric, nodes, reception_time=0., reception_x=.2,
        lookbacks=[lookback], direction=direction, quadrature_order=64)[0]

    def derivatives(chi):
        time, position = -chi, .2+direction*chi
        radius, scale = 2.+.15*time+.1*position, .8+.05*time-.03*position
        lapse = radius*scale
        radius_chi = -.15+direction*.1
        radius_chichi = 0.
        lapse_t = .15*scale+radius*.05
        lapse_x = .1*scale+radius*(-.03)
        lapse_tt, lapse_tx, lapse_xx = 2*.15*.05, .15*(-.03)+.1*.05, 2*.1*(-.03)
        lapse_chi = -lapse_t+direction*lapse_x
        lapse_chichi = -(-lapse_tt+direction*lapse_tx)+direction*(-lapse_tx+direction*lapse_xx)
        return lapse, radius, lapse_chi, lapse_chichi, radius_chi, radius_chichi

    def flow(_chi, state):
        distance, slope = state
        lapse, radius, lapse_chi, _lapse_chichi, radius_chi, radius_chichi = derivatives(_chi)
        coefficient = radius_chichi/radius-2*(lapse_chi/lapse)*(radius_chi/radius)
        return [slope, 2*(lapse_chi/lapse)*slope+coefficient*distance]

    opening = derivatives(0.)[0]
    solved = solve_ivp(flow, (0., lookback), [0., opening], rtol=1e-10, atol=1e-12)
    assert row['angular_diameter_distance'] == pytest.approx(solved.y[0, -1], rel=1e-10)
    assert row['luminosity_distance'] == pytest.approx(
        (row['observer_lapse']/row['emission_lapse'])**2*row['angular_diameter_distance'])


def test_lightcone_refuses_time_extrapolation():
    with pytest.raises(ValueError):
        radial_conformal_lightcone(
            lambda time, position: {'radius': 1.+position, 'Q': 1.},
            [0., 1.], reception_time=.2, reception_x=0., lookbacks=[.3], direction=1)


def test_bad_covariance_and_empty_band_are_rejected():
    with pytest.raises(ValueError, match='symmetric'):
        contrast_weights([[1., .2], [.1, 1.]], [0, 1])
    with pytest.raises(ValueError, match='positive covariance'):
        contrast_weights([[1., 2.], [2., 1.]], [0, 1])
    with pytest.raises(ValueError, match='two populated'):
        contrast_weights(np.eye(2), [0, 0])
    with pytest.raises(ValueError):
        quiet_cap_distance_shape([.1], [-1.])
