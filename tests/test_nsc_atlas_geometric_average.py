"""Independent constraint, wall projection and quadrature controls."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_connection_slice import wall_slice_data
from recursive_horizons.nsc_atlas_geometric_average import (
    KAPPA, geometric_slice_average, scalar_average_sources, wall_slice_stress,
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
