"""Geometry checks of the formal f_s owner; no source generators or assembly."""
from dataclasses import replace
import time

import numpy as np
import pytest

from recursive_horizons.nsc_compatible_history_geometry import (
    CompatibleIncomingMetric, CompatibleRadiusDirection,
)
from recursive_horizons.nsc_dirac_source_phase import (
    SOURCE_SIGNS, coefficient_binding, formal_source_phase_coefficient,
    _basis_jets, _characteristic_distance, _normal_coordinate,
)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_source_envelope import (
    RHO_UP_MIN, axial_center, continuum_speed_distance, physical_incoming_interval,
)
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction, LocalIncomingFamily


RHO_UP = 1.03
ANGULAR = np.sqrt(5.)
GAUSS = 8
Z_FD = 1e-6
AMP_FD = 1e-6


def targets():
    left, right = physical_incoming_interval()
    return np.array([left, 0.5 * (left + right), right])


def n8_family():
    coeff = np.zeros((2, 8))
    coeff[0, :4] = [0.03, -0.02, 0.01, 0.005]
    coeff[1, :3] = [1.2, -0.4, 0.25]
    return LocalIncomingFamily(coeff)


def bump(amplitude=0.002, w=(0.25, 0.15, -0.05), u=(0.8, -0.2)):
    center = axial_center()
    direction = CompatibleRadiusDirection(
        LocalAxialFunction(w, center), LocalAxialFunction(u, center), 0.007, 0.03)
    return CompatibleIncomingMetric((float(amplitude),), (direction,))


def phase(family, z=None, **kwargs):
    options = dict(angular=ANGULAR, rho_up=RHO_UP, gauss_nodes=GAUSS)
    options.update(kwargs)
    return formal_source_phase_coefficient(family, targets() if z is None else z, **options)


def test_ell0_is_exact_zero_without_quadrature_and_bounds_stay_none():
    family = n8_family()
    started = time.process_time()
    result = phase(family, angular=0.0)
    assert result.quadrature['skipped'] is True
    assert result.quadrature['n_nodes'] == 0
    assert np.count_nonzero(result.f) == np.count_nonzero(result.f_z) == 0
    assert np.count_nonzero(result.delta_f) == np.count_nonzero(result.delta_f_z) == 0
    assert result.coefficient_error_bound is result.remainder_error_bound is None
    assert result.physical_error_bound is None
    assert result.assembly_correction is False
    assert result.physical_local_gate.startswith('OPEN')
    assert result.signs == SOURCE_SIGNS == (1, -1)
    with pytest.raises(ValueError, match='1.03'):
        phase(family, rho_up=1.02)
    assert time.process_time() - started < 5.


def test_zero_amplitude_keeps_nonzero_tangents_and_frozen_upstream_distance():
    family = bump(0.0)
    result = phase(family)
    assert np.count_nonzero(result.f) == np.count_nonzero(result.f_z) == 0
    assert result.delta_f.shape == (2, 1, 3)
    assert np.max(np.abs(result.delta_f)) > 1e-8
    assert np.max(np.abs(result.delta_f_z)) > 1e-8
    assert result.characteristic_distance == _characteristic_distance(RHO_UP)
    assert result.characteristic_distance == continuum_speed_distance(RHO_UP, 1.0)
    assert result.rho_up == RHO_UP >= RHO_UP_MIN
    other = phase(bump(0.003))
    assert other.characteristic_distance == result.characteristic_distance
    assert other.profile_identity != result.profile_identity


def test_analytic_f_z_and_w_U_tangents_match_finite_differences():
    family = n8_family()
    center = family.center
    z = np.array([center - Z_FD, center, center + Z_FD])
    started = time.process_time()
    result = phase(family, z)
    f_z_fd = (result.f[:, 2] - result.f[:, 0]) / (2. * Z_FD)
    f_z_residual = float(np.max(np.abs(result.f_z[:, 1] - f_z_fd)))
    assert f_z_residual <= 2e-8
    np.testing.assert_allclose(result.delta_n2, np.asarray(SOURCE_SIGNS)[:, None] * result.f_z,
                               rtol=0, atol=0)

    coeff = np.array(family.coefficients, dtype=float)

    def shifted(index, step):
        altered = coeff.copy()
        altered.ravel()[index] += step
        return LocalIncomingFamily(altered)

    residuals = {}
    for name, index in (('w', 1), ('U', 8)):
        plus = phase(shifted(index, AMP_FD), z[1:2])
        minus = phase(shifted(index, -AMP_FD), z[1:2])
        expected = (plus.f - minus.f) / (2. * AMP_FD)
        expected_z = (plus.f_z - minus.f_z) / (2. * AMP_FD)
        residuals[name] = (
            float(np.max(np.abs(result.delta_f[:, index, 1:2] - expected))),
            float(np.max(np.abs(result.delta_f_z[:, index, 1:2] - expected_z))),
        )
        assert residuals[name][0] <= 2e-7
        assert residuals[name][1] <= 5e-6
    cpu = time.process_time() - started
    assert result.delta_f.shape == (2, 16, 3)
    assert np.max(np.abs(result.delta_f[:, 7, :])) > 0.
    assert not np.allclose(result.f[0], result.f[1])
    assert cpu < 25.


def test_reflection_is_odd_to_leading_order_not_finite_alpha_even():
    alpha = 0.002
    plus = phase(bump(alpha))
    minus = phase(bump(-alpha))
    origin = phase(bump(0.0))
    odd = 0.5 * (plus.f - minus.f)
    even = 0.5 * (plus.f + minus.f)
    linear = alpha * origin.delta_f[:, 0, :]
    assert plus.signs == minus.signs == SOURCE_SIGNS
    assert np.max(np.abs(plus.f - minus.f)) > 1e-8
    assert np.max(np.abs(even)) < np.max(np.abs(odd))
    assert np.max(np.abs(plus.f - minus.f)) > np.max(np.abs(plus.f + minus.f))
    np.testing.assert_allclose(odd, linear, rtol=2e-6, atol=2e-10)
    assert plus.characteristic_distance == minus.characteristic_distance == origin.characteristic_distance


def test_distinct_analytic_profiles_bind_windows_and_quadrature():
    center = axial_center()
    first = bump(0.001, w=(0.2, 0.1), u=(0.5,))
    second = bump(0.001, w=(0.2, -0.1), u=(0.5,))
    third = CompatibleIncomingMetric(
        first.amplitudes,
        (replace(first.directions[0], inner_radius=0.006),))
    a, b = phase(first), phase(second)
    assert a.profile_identity != b.profile_identity
    assert a.binding != b.binding
    assert np.max(np.abs(a.f - b.f)) > 1e-12
    assert profile_identity(first) != profile_identity(third)
    digest, description = coefficient_binding(
        first, targets(), angular=ANGULAR, rho_up=RHO_UP, gauss_nodes=GAUSS)
    assert a.binding == digest
    assert description['gauss_nodes'] == GAUSS
    assert description['source_signs'] == [1, -1]
    assert description['profile']['normal_windows_included'] is True
    other = coefficient_binding(
        first, targets(), angular=ANGULAR, rho_up=RHO_UP, gauss_nodes=GAUSS + 1)[0]
    assert other != digest
    assert a.characteristic_distance == b.characteristic_distance
    metric = first
    rho = 1.008
    z = np.array([center])
    basis, basis_z = _basis_jets(metric.directions[0], _normal_coordinate(rho), z)
    np.testing.assert_allclose(basis[0], metric.directions[0].value(_normal_coordinate(rho), float(z[0])),
                               rtol=0, atol=1e-14)
    h = 1e-7
    value_z = (metric.directions[0].value(_normal_coordinate(rho), float(z[0] + h))
               - metric.directions[0].value(_normal_coordinate(rho), float(z[0] - h))) / (2. * h)
    assert abs(float(basis_z[0]) - value_z) <= 2e-8


def test_tiny_history_preserves_the_perturbation_before_radius_addition():
    alpha = 1e-20
    tiny, zero = phase(bump(alpha)), phase(bump(0.))
    expected = alpha*zero.delta_f[:, 0]
    assert np.max(abs(tiny.f)) > 0
    np.testing.assert_allclose(tiny.f, expected, rtol=2e-14, atol=0)
    assert tiny.quadrature['convergence_indicator'] is None
    assert 'panel_absolute_contribution_sums' in tiny.quadrature
    assert tiny.quadrature['radius_lower_bound'] > 0


def test_positive_geometry_and_flat_upstream_are_required_between_nodes():
    coefficient = np.zeros((2, 8)); coefficient[0, 0] = 100.
    with pytest.raises(ValueError, match='radius bound'):
        phase(LocalIncomingFamily(coefficient))
    g = bump()
    nonflat = CompatibleIncomingMetric(g.amplitudes,
        (replace(g.directions[0], outer_radius=.1),))
    with pytest.raises(ValueError, match='flat at the fixed upstream'):
        phase(nonflat)
