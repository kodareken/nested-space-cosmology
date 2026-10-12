"""Observable cross/quadratic cancellation and full retarded derivative."""
from dataclasses import replace

import numpy as np

from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter
from recursive_horizons.nsc_ks_matter_difference import source_fixed_matter_difference
from recursive_horizons.nsc_pg_ks_metric_pullback import reference_chart
from test_nsc_ks_difference_envelope import evolve


PARAMETERS = dict(axial_scale=reference_chart(1.)[2], radius=np.sqrt(2), multiplicity=3.)


def explicit_difference(result):
    ref = np.exp(-1j * result.source_energies * result.z[:, None])[:, None, :] * result.reference_amplitudes[None]
    refz = -1j * result.source_energies * ref
    zeros = np.zeros((0, *ref.shape), complex)
    params = dict(mass=result.mass, angular=result.angular, **PARAMETERS)
    current = source_column_matter(result.weighted_columns, result.weighted_axial_columns,
        result.source_covariance, result.weighted_column_tangents, result.weighted_axial_tangents, **params)
    reference = source_column_matter(ref * result.column_weights, refz * result.column_weights,
        result.source_covariance, zeros, zeros, **params)
    return current['action_gradient'] - reference['action_gradient']


def test_difference_preserves_coherence_weights_and_quadratic_terms():
    result = evolve(.4)
    stable = source_fixed_matter_difference(result, **PARAMETERS)
    np.testing.assert_allclose(stable['action_gradient_change'], explicit_difference(result), atol=3e-15, rtol=0)
    zero = source_fixed_matter_difference(evolve(0.), **PARAMETERS)
    assert np.count_nonzero(zero['action_gradient_change']) == 0
    assert stable['physical_constraint_status'] == 'OPEN'


def test_difference_derivative_is_full_fixed_source_tangent():
    h = 1e-4
    base, plus, minus = [source_fixed_matter_difference(evolve(a), **PARAMETERS)
                         for a in (.002, .002 + h, .002 - h)]
    fd = (plus['action_gradient_change'] - minus['action_gradient_change']) / (2 * h)
    np.testing.assert_allclose(fd, base['action_gradient_tangent'][0], atol=3e-10, rtol=0)


def test_pure_common_phase_cancels_in_observables():
    result = evolve(0.)
    # A constant common column phase changes F but neither C nor the momentum.
    phase = np.exp(.07j)
    delta = np.broadcast_to((phase - 1) * result.reference_amplitudes, result.columns.shape).copy()
    total = result.reference_amplitudes[None] + delta
    carrier = np.exp(-1j * result.source_energies * result.z[:, None])[:, None, :]
    changed = replace(result, envelope_difference=delta,
                      columns=carrier * total, axial_columns=carrier * (-1j * result.source_energies * total))
    stable = source_fixed_matter_difference(changed, **PARAMETERS)
    assert np.max(abs(stable['action_gradient_change'])) < 3e-17
    # This is an algebraic norm counterexample, not a solved radius history.
