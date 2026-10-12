"""Synthetic exact-phase incoming controls; no authenticated 801/3201 runs."""
import numpy as np
import pytest
from unittest.mock import patch
from scipy.linalg import expm
from scipy.sparse import bmat, csr_matrix, diags

from recursive_horizons.nsc_compatible_history_geometry import (
    CompatibleIncomingMetric, CompatibleRadiusDirection,
)
from recursive_horizons.nsc_evolved_incoming_state import (
    AmplitudeOnlyMetric, CachedCompatibleIncomingMetric, FixedSourcePreparation,
)
from recursive_horizons.nsc_exact_phase_prepared_state import (
    ExactPhaseIncoming, prepare_exact_phase_incoming,
)
from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates
from recursive_horizons.nsc_prepared_history_jets import evolve_prepared_field_jets
from recursive_horizons.nsc_transmitting_history_jets import (
    FourthOrderModePropagator, incoming_reference_columns,
)
from recursive_horizons.nsc_transmitting_resolvent import profile


def signed_source():
    return FixedSourcePreparation.from_signed_blocks(
        np.array([0.4, -0.55]),
        np.array([1.2, 0.7]),
        np.array([
            [[0.45, 0.07 + 0.04j], [0.07 - 0.04j, 0.6]],
            [[0.50, 0.02 - 0.03j], [0.02 + 0.03j, 0.35]],
        ]),
    )


def delayed_bump(center, halfwidth):
    def value(z, order):
        if order not in (0, 1, 2, 3):
            raise ValueError('owned incoming derivative order 0..3 required')
        u = (z - center) / halfwidth
        if abs(u) >= 1:
            return 0.
        b, bp = profile(u)
        if order == 0:
            return b
        if order == 1:
            return bp / halfwidth
        d = 1 - u * u
        g1 = -2 * u / d ** 2
        g2 = -2 / d ** 2 - 8 * u * u / d ** 3
        if order == 2:
            return b * (g1 * g1 + g2) / halfwidth ** 2
        g3 = -24 * u / d ** 3 - 48 * u ** 3 / d ** 4
        return b * (g1 ** 3 + 3 * g1 * g2 + g3) / halfwidth ** 3
    return value


def delayed_family(amplitude=0.003, inner=0.05, outer=0.25, shift=0.025, halfwidth=0.018):
    x = np.linspace(-2., 1.2, 17)
    i8 = int(np.argmin(np.abs(x - 0.8)))
    center = chart_coordinates(float(x[i8]))[1] + shift
    directions = (CompatibleRadiusDirection(
        delayed_bump(center, halfwidth), lambda z, n: 0., inner, outer),)
    return CompatibleIncomingMetric((float(amplitude),), directions)


def initial_columns(owner, nsrc):
    x = owner.x
    n = len(x)
    rng = np.linspace(0., 1., n)
    columns = np.zeros((2 * n, nsrc), dtype=complex)
    for k in range(nsrc):
        columns[:n, k] = (0.2 + 0.03 * k) * np.exp(-rng * rng) + 0.04j * np.sin((k + 1) * rng)
        columns[n:, k] = (0.15 - 0.02 * k) * np.exp(-0.3 * rng * rng) + 0.05j * np.cos((k + 2) * rng)
        columns[n - 1, k] += 0.12 + 0.03j * (k + 1)
        columns[2 * n - 1, k] += 0.09 - 0.02j * (k + 1)
    return columns


def setup(amplitude=0.003, ndir=1):
    x = np.linspace(-2., 1.2, 17)
    owner = FourthOrderModePropagator(x, 1.3, 2.)
    source = signed_source()
    phi = initial_columns(owner, source.covariance.shape[0])
    times = np.linspace(0., 0.04, 5)
    family = delayed_family(amplitude)
    provider = family if ndir else AmplitudeOnlyMetric(family)
    return owner, provider, family, times, phi, source


def test_cached_evolution_and_binding_never_rebuild_uncached_chart():
    owner, _, family, times, phi, source = setup()
    cached = CachedCompatibleIncomingMetric(owner.x, family)
    with patch.object(CompatibleIncomingMetric, '_evaluate', side_effect=AssertionError('uncached chart')):
        result = prepare_exact_phase_incoming(owner, times, cached, phi, source)
        result.require_history(cached, times, owner.x)
    result.require_history(family, times, owner.x)


def test_constant_augmented_matrix_matches_independent_expm():
    owner, provider, family, times, phi, source = setup(amplitude=0., ndir=0)
    result = prepare_exact_phase_incoming(owner, times, provider, phi, source)
    rows, nsrc = phi.shape
    L, checks = owner.generator(owner.reference)
    B = incoming_reference_columns(owner, phi)
    A = bmat(
        [[L, csr_matrix(B)], [None, diags(-1j * source.energies)]],
        format='csr',
    ).toarray()
    y0 = np.vstack((phi, np.eye(nsrc, dtype=complex)))
    independent = np.stack([expm((t - times[0]) * A) @ y0 for t in times])
    assert result.state.column_tangents.shape == (0, len(times), 2, nsrc)
    assert result.axial_tangents.shape == (0, len(times), 2, nsrc)
    assert np.max(np.abs(independent[:, :rows] - result.diagnostics['node_fields'])) < 3e-12
    assert np.max(np.abs(independent[:, rows:] - result.diagnostics['phase'])) < 3e-12
    assert result.diagnostics['exact_harmonic_phase_residual'] < 3e-12
    assert result.diagnostics['midpoint_held_forcing'] is False
    assert checks['norm_flux_algebra'] < 3e-11
    assert result.diagnostics['norm_flux_residual'] < 3e-11
    assert result.diagnostics['CF4_time_integrator'] is False
    take = [int(np.flatnonzero(owner.x == 1.)[0]), len(owner.x) + int(np.flatnonzero(owner.x == 1.)[0])]
    R = result.state.history.restriction
    for k, y in enumerate(independent):
        rhs = L @ y[:rows] + B @ y[rows:]
        assert np.max(np.abs(R @ y[take] - result.state.columns[k])) < 3e-12
        assert np.max(np.abs(R @ rhs[take] - result.axial_columns[k])) < 3e-12
    result.validate()
    result.require_history(family, times, owner.x)
    result.require_history(AmplitudeOnlyMetric(family), times, owner.x)


def test_nonconstant_radius_finite_difference_includes_axial_tangent():
    owner, provider, family, times, phi, source = setup(0.003, ndir=1)
    center = prepare_exact_phase_incoming(owner, times, provider, phi, source)
    step = 1e-4
    plus = prepare_exact_phase_incoming(
        owner, times, AmplitudeOnlyMetric(delayed_family(0.003 + step)), phi, source)
    minus = prepare_exact_phase_incoming(
        owner, times, AmplitudeOnlyMetric(delayed_family(0.003 - step)), phi, source)
    field_fd = (plus.diagnostics['node_fields'] - minus.diagnostics['node_fields']) / (2 * step)
    trace_fd = (plus.state.columns - minus.state.columns) / (2 * step)
    axial_fd = (plus.axial_columns - minus.axial_columns) / (2 * step)
    assert center.diagnostics['maximum_radius_change_from_background'] > 0
    assert np.max(np.abs(center.state.column_tangents)) > 1e-8
    assert np.max(np.abs(center.axial_tangents)) > 1e-8
    assert np.max(np.abs(field_fd - center.diagnostics['node_tangents'][0])) < 3e-8
    assert np.max(np.abs(trace_fd - center.state.column_tangents[0])) < 3e-8
    assert np.max(np.abs(axial_fd - center.axial_tangents[0])) < 3e-8
    assert center.diagnostics['norm_flux_residual'] < 3e-11
    assert center.diagnostics['metric_flux_tangent_residual'] < 3e-11
    assert isinstance(center.diagnostics['cache_class'], str)
    assert center.diagnostics['cache_class'] == 'CachedCompatibleIncomingMetric'
    assert plus.state.column_tangents.shape[0] == 0
    energy = source.energies
    frozen = -1j * energy * center.state.columns
    assert np.max(np.abs(center.axial_columns - frozen)) > 1e-6
    assert center.diagnostics['outgoing_momentum_k_equals_minus_E'] is False


def test_exact_source_phase_not_midpoint_held_forcing():
    owner, provider, family, times, phi, source = setup(0., ndir=0)
    exact = prepare_exact_phase_incoming(owner, times, provider, phi, source)
    expected = np.eye(phi.shape[1])[None] * np.exp(
        -1j * (times - times[0])[:, None, None] * source.energies[None, None, :])
    assert np.max(np.abs(exact.diagnostics['phase'] - expected)) < 3e-12
    assert exact.diagnostics['exact_harmonic_phase_residual'] < 3e-12
    B = incoming_reference_columns(owner, phi)
    zeros = np.zeros((0, *phi.shape), complex)

    def incoming(time):
        phase = np.exp(-1j * source.energies * (time - times[0]))
        return B * phase[None, :], zeros

    held = evolve_prepared_field_jets(owner, times, provider, phi, zeros, incoming)
    assert np.max(np.abs(held['evolved_field'] - exact.state.diagnostics['evolved_field'])) > 1e-8
    assert held['auxiliary_identity_residual'] < 3e-14
    assert exact.diagnostics['midpoint_held_forcing'] is False
    owner2, provider2, _, times2, phi2, source2 = setup(0.003, ndir=1)
    moved = prepare_exact_phase_incoming(owner2, times2, provider2, phi2, source2)
    expected2 = np.eye(phi2.shape[1])[None] * np.exp(
        -1j * (times2 - times2[0])[:, None, None] * source2.energies[None, None, :])
    assert np.max(np.abs(moved.diagnostics['phase'] - expected2)) < 3e-12


def test_fixed_preparation_initial_support_and_intrinsic_cache_rejection():
    owner, provider, family, times, phi, source = setup(0.003, ndir=1)
    with pytest.raises(ValueError, match='frozen incoming C0'):
        prepare_exact_phase_incoming(owner, times, provider, phi, source, incoming_C0=np.eye(4))
    with pytest.raises(ValueError, match='frozen incoming C0'):
        prepare_exact_phase_incoming(
            owner, times, provider, phi, {'stress_approximant': np.ones(4), 'parts': {}})
    with pytest.raises(TypeError, match='FixedSourcePreparation'):
        prepare_exact_phase_incoming(owner, times, provider, phi, np.eye(4))
    wave = lambda z, order: float(np.real((1j * 0.7) ** order * np.exp(1j * 0.7 * z)))
    live = CompatibleIncomingMetric(
        (0.003,), (CompatibleRadiusDirection(wave, lambda z, n: 0., 0.05, 0.25),))
    with pytest.raises(ValueError, match='reference-support guard'):
        prepare_exact_phase_incoming(owner, times, live, phi, source)
    rho1 = int(np.flatnonzero(owner.x == 1.)[0])

    class ChangingIntrinsic:
        def values(self, time, x):
            metric_values = owner.reference.copy()
            metric_values[0, rho1] *= 1. + 0.02 * float(time)
            return metric_values

        def log_directions(self, time, x, metric_values):
            directions = np.zeros((1, 4, len(x)))
            directions[0, 0, rho1] = 0.02
            return directions

    with pytest.raises(ValueError, match='unsupported intrinsic-changed family'):
        prepare_exact_phase_incoming(owner, times, ChangingIntrinsic(), phi, source)
    other_grid = np.linspace(-2., 1.2, 21)
    cache = CachedCompatibleIncomingMetric(other_grid, family)
    with pytest.raises(ValueError, match='cached provider grid'):
        prepare_exact_phase_incoming(owner, times, cache, phi, source)
    result = prepare_exact_phase_incoming(owner, times, family, phi, source)
    other = delayed_family(0.003, shift=0.030)
    assert tuple(other.amplitudes) == result.state.history.amplitudes
    assert tuple(d.inner_radius for d in other.directions) == result.state.history.inner_radii
    with pytest.raises(ValueError, match='actual sampled'):
        result.require_history(other, times, owner.x)
    result.require_history(family, times, owner.x)
    with pytest.raises(ValueError, match='sampled-history fingerprint'):
        ExactPhaseIncoming(
            state=result.state,
            axial_columns=result.axial_columns,
            axial_tangents=result.axial_tangents,
            sampled_node_metrics=result.sampled_node_metrics,
            sampled_midpoint_metrics=result.sampled_midpoint_metrics,
            sampled_node_directions=result.sampled_node_directions,
            sampled_midpoint_directions=result.sampled_midpoint_directions,
            sampled_history_fingerprint='0' * 64,
            fixed_preparation_digest=result.fixed_preparation_digest,
            diagnostics={**dict(result.diagnostics), 'sampled_history_fingerprint': '0' * 64},
        )


def test_signed_fibers_weights_and_history_fingerprint():
    owner, provider, family, times, phi, source = setup(0.003, ndir=1)
    result = prepare_exact_phase_incoming(owner, times, provider, phi, source)
    nsrc = 4
    assert result.state.columns.shape == (len(times), 2, nsrc)
    assert result.axial_columns.shape == (len(times), 2, nsrc)
    assert result.state.column_tangents.shape == (1, len(times), 2, nsrc)
    assert result.axial_tangents.shape == (1, len(times), 2, nsrc)
    assert np.array_equal(result.state.source_energies, np.repeat([0.4, -0.55], 2))
    expected_weights = np.repeat(np.sqrt(np.array([1.2, 0.7]) / (2 * np.pi)), 2)
    assert np.max(np.abs(result.state.column_weights - expected_weights)) < 1e-15
    assert np.max(np.abs(result.weighted_axial_columns - result.axial_columns * expected_weights)) < 1e-15
    assert np.max(np.abs(
        result.weighted_axial_tangents - result.axial_tangents * expected_weights)) < 1e-15
    twice = result.axial_columns * expected_weights * expected_weights
    assert np.max(np.abs(result.weighted_axial_columns - twice)) > 1e-4
    kernel = result.state.covariance()
    assert kernel.shape == (len(times), 2, len(times), 2)
    assert np.max(np.abs(kernel[0, :, -1, :])) > 1e-8
    F = result.weighted_columns.reshape(len(times) * 2, nsrc)
    manual = (F @ source.covariance @ F.conj().T).reshape(len(times), 2, len(times), 2)
    assert np.max(np.abs(manual - kernel)) < 1e-14
    assert result.fixed_preparation_digest == result.diagnostics['fixed_preparation_digest']
    changed_phi = phi.copy()
    changed_phi[0, 0] += 0.2
    other = prepare_exact_phase_incoming(owner, times, provider, changed_phi, source)
    assert other.fixed_preparation_digest != result.fixed_preparation_digest
    assert result.state.coverage['physical_history_solution'] == 'OPEN'
    assert result.state.coverage['continuum_error_bound'] is None
    assert result.diagnostics['full_source_error_bound'] is None
    assert result.diagnostics['metric_evolution'] is False
    assert result.diagnostics['optimizer'] is False
    assert result.diagnostics['callback_identity_claimed'] is False
    z = times + chart_coordinates(1.)[1]
    assert np.max(np.abs(result.state.z - z)) < 1e-15
    result.validate()
    result.require_history(CachedCompatibleIncomingMetric(owner.x, family), times, owner.x)
