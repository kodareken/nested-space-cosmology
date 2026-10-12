"""Source-preserving operator reconstruction controls; no source campaign."""
import numpy as np
import pytest

from recursive_horizons.nsc_ks_energy_propagator import (
    chebyshev_energy_nodes, evolve_energy_propagator, interpolation_matrix,
)
from recursive_horizons.nsc_ks_difference_envelope import evolve_ks_difference_envelope
from recursive_horizons.nsc_ks_source_envelope import (
    computational_z_grid, physical_incoming_interval, usual_axial_support,
)
from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from test_nsc_ks_batched_constraints import make_panel, bind, MASS, ANGULAR, RHO_UP
from test_nsc_ks_local_constraints import bump_family, coefficients_for


SOLVER = dict(rtol=8e-12, atol=8e-14, max_step=.001)


def test_actual_node_polynomial_and_endpoints():
    interval = (.2, 4.)
    nodes = chebyshev_energy_nodes(interval, 12)
    target = np.r_[interval, nodes, .8, 2.7]
    matrix = interpolation_matrix(nodes, target, interval)
    for degree in (0, 1, 4, 12):
        x = (nodes-2.1)/1.9
        y = (target-2.1)/1.9
        np.testing.assert_allclose(matrix @ (x**degree), y**degree, atol=2e-14, rtol=0)
    np.testing.assert_array_equal(matrix[2:2+len(nodes)], np.eye(len(nodes)))
    with pytest.raises(ValueError, match='outside'):
        interpolation_matrix(nodes, np.array([4.01]), interval)
    with pytest.raises(ValueError, match='integer'):
        chebyshev_energy_nodes(interval, 3.5)


def setup_operator():
    family = bump_family()
    grid = computational_z_grid(32)
    target = np.linspace(*physical_incoming_interval(), 4)
    operator = evolve_energy_propagator(
        (.2, 3.2), 8, family, grid, target, MASS, ANGULAR, RHO_UP,
        axial_support=usual_axial_support(), **SOLVER)
    return family, grid, target, operator


def test_original_coherent_source_and_all_four_fields_match_direct_evolution():
    family, grid, target, operator = setup_operator()
    batch = bind(make_panel(energies=(.2, .53, 1.1, 3.2), weights=(.13, .2, .31, .23)))
    reconstructed = operator.apply(batch.source, batch.initial_columns)
    direct = evolve_ks_difference_envelope(
        batch.source, batch.initial_columns, family, grid, target, MASS, ANGULAR, RHO_UP,
        axial_support=usual_axial_support(), **SOLVER)
    for key in ('columns', 'axial_columns', 'column_tangents', 'axial_tangents'):
        np.testing.assert_allclose(getattr(reconstructed, key), getattr(direct, key), atol=3e-9, rtol=0)
    np.testing.assert_array_equal(reconstructed.source_covariance, batch.source.covariance)
    np.testing.assert_array_equal(reconstructed.source_energies, batch.source.energies)
    np.testing.assert_array_equal(reconstructed.column_weights, batch.source.column_weights)
    np.testing.assert_array_equal(reconstructed.initial_columns, batch.initial_columns)
    assert reconstructed.fixed_preparation_digest == direct.fixed_preparation_digest
    reconstructed.require_history(family)
    assert reconstructed.diagnostics['physical_constraint_status'] == 'OPEN'
    assert reconstructed.diagnostics['interpolation_error_bound'] is None
    coefficients = coefficients_for()
    kwargs = dict(mass=MASS, angular=ANGULAR, axial_scale=coefficients['a'],
                  radius=coefficients['r'], multiplicity=12.)
    def matter(state, covariance):
        return source_column_matter(
            state.weighted_columns, state.weighted_axial_columns, covariance,
            state.weighted_column_tangents, state.weighted_axial_tangents, **kwargs)
    actual, check = matter(reconstructed, batch.source.covariance), matter(direct, batch.source.covariance)
    np.testing.assert_allclose(actual['action_gradient'], check['action_gradient'], atol=3e-9, rtol=0)
    np.testing.assert_allclose(actual['action_gradient_tangent'], check['action_gradient_tangent'], atol=2e-8, rtol=0)
    incoherent = matter(reconstructed, np.diag(np.diag(batch.source.covariance)))
    assert np.max(abs(actual['action_gradient']-incoherent['action_gradient'])) > 1e-5


def test_new_source_changes_state_without_changing_operator_and_extrapolation_rejected():
    _, _, _, operator = setup_operator()
    first = bind(make_panel(energies=(.7,), weights=(.2,), seed=10))
    second = bind(make_panel(energies=(.7,), weights=(.2,), seed=11))
    before = operator.digest
    x, y = operator.apply(first.source, first.initial_columns), operator.apply(second.source, second.initial_columns)
    assert np.max(abs(x.columns-y.columns)) > .1
    assert x.fixed_preparation_digest != y.fixed_preparation_digest
    assert operator.digest == before
    outside = FixedSourcePreparation(first.source.covariance, first.source.column_weights,
                                     np.repeat(3.3, 3))
    with pytest.raises(ValueError, match='outside'):
        operator.apply(outside, first.initial_columns)
    with pytest.raises(TypeError, match='original'):
        operator.apply(np.eye(3), first.initial_columns)


def test_zero_tangent_operator_keeps_empty_direction_axis():
    family = bump_family()
    grid = computational_z_grid(16)
    target = np.linspace(*physical_incoming_interval(), 3)
    operator = evolve_energy_propagator(
        (.2, 3.2), 4, family, grid, target, MASS, ANGULAR, RHO_UP,
        axial_support=usual_axial_support(), tangents='zero', **SOLVER)
    batch = bind(make_panel(energies=(.7,), weights=(.2,)))
    state = operator.apply(batch.source, batch.initial_columns)
    assert state.column_tangents.shape == (0, 3, 2, 3)
    assert state.axial_tangents.shape == (0, 3, 2, 3)
