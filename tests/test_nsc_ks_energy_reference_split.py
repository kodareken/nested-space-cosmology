"""Interpolate the physical change, without an interpolated-reference offset."""
from dataclasses import replace
import numpy as np
import pytest

from recursive_horizons.nsc_ks_energy_propagator import evolve_energy_propagator
from recursive_horizons.nsc_ks_source_envelope import computational_z_grid, physical_incoming_interval, usual_axial_support
from recursive_horizons.nsc_ks_local_constraints import homogeneous_reference_amplitudes
from recursive_horizons.nsc_ks_matter_difference import source_fixed_matter_difference
from test_nsc_ks_energy_propagator import setup_operator, make_panel, bind
from test_nsc_ks_signed_state import bump_family, MASS, ANGULAR, RHO_UP, SOLVER
from test_nsc_ks_local_constraints import coefficients_for


def test_zero_history_reconstructs_the_same_reference_without_a_stress_offset():
    family = bump_family(0.0)
    target = np.linspace(*physical_incoming_interval(), 7)
    op = evolve_energy_propagator(
        (.2, 3.2), 8, family, computational_z_grid(32), target,
        MASS, ANGULAR, RHO_UP, tangents='zero', axial_support=usual_axial_support(), **SOLVER)
    batch = bind(make_panel(energies=(.2, .53, 1.1, 3.2), weights=(.13, .2, .31, .23)))
    state = op.apply(batch.source, batch.initial_columns)
    np.testing.assert_array_equal(state.reference_amplitudes, homogeneous_reference_amplitudes(state))
    assert np.count_nonzero(state.envelope_difference) == 0
    coefficients = coefficients_for()
    matter = source_fixed_matter_difference(
        state, axial_scale=coefficients['a'], radius=coefficients['r'],
        multiplicity=12., reference='cached')
    assert np.count_nonzero(matter['action_gradient_change']) == 0
    assert matter['stress_drift_subtracted'] is False
    assert state.diagnostics['reference_uses_original_energies'] is True
    assert state.diagnostics['fixed_reference_error_bound'] is None


def test_reference_mode_is_bound_and_old_operator_records_keep_their_semantics():
    _family, _grid, _target, direct = setup_operator()
    historical_diagnostics = {key: value for key, value in direct.diagnostics.items()
                              if key != 'reference_mode'}
    historical = replace(direct, diagnostics=historical_diagnostics)
    assert direct.digest != historical.digest
    batch = bind(make_panel(energies=(.7,), weights=(.2,), seed=10))
    current = direct.apply(batch.source, batch.initial_columns)
    old = historical.apply(batch.source, batch.initial_columns)
    assert current.diagnostics['reference_mode'] == 'direct-original-energies'
    assert old.diagnostics['reference_mode'] == 'interpolated'
    np.testing.assert_array_equal(current.envelope_difference, old.envelope_difference)
    np.testing.assert_array_equal(current.column_tangents, old.column_tangents)
    np.testing.assert_array_equal(current.initial_columns, old.initial_columns)
    with pytest.raises(ValueError, match='unknown energy reference'):
        replace(direct, diagnostics={'reference_mode': 'subtract-measured-offset'})
