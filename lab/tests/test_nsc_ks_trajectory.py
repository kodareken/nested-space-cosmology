"""Endpoint continuity and unchanged evolution under reconstruction capture."""
import numpy as np

from recursive_horizons.nsc_ks_trajectory import (
    TrajectorySegment, capture_ks_trajectory, TrajectoryResidualSampler,
)
from recursive_horizons.nsc_ks_source_envelope import computational_z_grid, physical_incoming_interval, usual_axial_support
from recursive_horizons.nsc_ks_difference_envelope import evolve_ks_difference_envelope
from test_nsc_ks_difference_envelope import inputs


def test_anchored_polynomial_endpoints_and_derivative():
    rng = np.random.default_rng(6)
    start, end = rng.normal(size=(2, 4))
    coefficients = rng.normal(size=(6, 4))
    part = TrajectorySegment(1.03, 1.02, start, end, coefficients)
    np.testing.assert_array_equal(part.evaluate(0), start)
    np.testing.assert_array_equal(part.evaluate(1), end)
    h = 1e-5
    numerical = (part.evaluate(.37+h) - part.evaluate(.37-h)) / (2*h*(part.rho_end-part.rho_start))
    np.testing.assert_allclose(part.evaluate(.37, derivative=True), numerical, atol=2e-7, rtol=0)


def test_capture_preserves_the_actual_prepared_state_without_patching_globals():
    args = (*inputs(), computational_z_grid(16), np.linspace(*physical_incoming_interval(), 5), 1.5, np.sqrt(5), 1.03)
    options = dict(axial_support=usual_axial_support(), rtol=2e-12, atol=2e-14, max_step=.001, tangents='zero')
    before = evolve_ks_difference_envelope.__globals__['DOP853']
    original = evolve_ks_difference_envelope(*args, **options)
    prepared, segments = capture_ks_trajectory(*args, **options)
    assert evolve_ks_difference_envelope.__globals__['DOP853'] is before
    for name in ('columns', 'axial_columns', 'reference_amplitudes', 'envelope_difference'):
        np.testing.assert_array_equal(getattr(original, name), getattr(prepared, name))
    assert prepared.fixed_preparation_digest == original.fixed_preparation_digest
    assert len(segments) == original.diagnostics['accepted_steps']
    sampler = TrajectoryResidualSampler(prepared, inputs()[2])
    indicator = sampler.sample(segments[len(segments)//2], .5)
    assert np.isfinite(indicator['weighted_L2_indicators']).all()
    assert indicator['continuous_integral_bounds'] is None
    assert indicator['is_error_certificate'] is False
