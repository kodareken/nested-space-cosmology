"""Source-fixed incoming state: actual evolution, frozen-input rejection, replay."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_compatible_history_geometry import (
    CompatibleIncomingMetric, CompatibleRadiusDirection,
)
from recursive_horizons.nsc_evolved_incoming_state import (
    AmplitudeOnlyMetric, CachedCompatibleIncomingMetric, EvolvedIncomingState,
    FixedSourcePreparation, HistoryBinding, evolve_incoming_state,
    restriction_covariance_tangent, rho1_node_index, rho1_restriction_map,
)
import recursive_horizons.nsc_evolved_incoming_state as STATE
from recursive_horizons.nsc_transmitting_history_jets import FourthOrderModePropagator

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'evolved_incoming_control', ROOT / 'scripts/derive_nsc_evolved_incoming_state.py')
CONTROL = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CONTROL)


def wave(frequency, phase=0., scale=1.):
    return lambda z, order: float(
        scale * np.real((1j * frequency) ** order * np.exp(1j * (frequency * z + phase))))


def tiny_source():
    return FixedSourcePreparation.from_signed_blocks(
        np.array([0.4]), np.array([1.2]),
        np.array([[[0.45, 0.07 + 0.04j], [0.07 - 0.04j, 0.6]]]))


def tiny_family(amplitude=0.003):
    x = np.linspace(-2., 1.2, 17)
    owner = FourthOrderModePropagator(x, 1.3, 2.)
    directions = (CompatibleRadiusDirection(wave(0.7), wave(0.4, 0.3), 0.015, 0.08),)
    metric = CompatibleIncomingMetric((float(amplitude),), directions)
    times = np.linspace(0., 0.04, 5)
    n = len(x)
    base = np.vstack((
        np.column_stack((np.exp(-x * x), 0.2 * np.exp(0.2j * x))),
        np.column_stack((0.3 * np.exp(0.1j * x), np.exp(-0.3 * x * x))),
    ))
    right = np.array([n - 1, 2 * n - 1])
    N, beta, q, _ = owner.reference
    speed = np.array([N[-1] / q[-1] - beta[-1], -N[-1] / q[-1] - beta[-1]])
    sat = -speed / owner.weights[-1]
    incident = np.array([[0.4 + 0.1j, -0.2j], [0.2, 0.3 - 0.1j]])

    def incoming(ndir):
        def supplied(time):
            forcing = np.zeros_like(base)
            forcing[right] = sat[:, None] * incident * np.exp(-1j * np.array([0.2, 0.6]) * time)
            return forcing, np.zeros((ndir, *base.shape), complex)
        return supplied

    return owner, metric, times, base, incoming, tiny_source()


def test_cache_matches_owned_provider_at_nonzero_radius():
    guard = CONTROL.RETARDED.support_guard()
    lower, upper = guard['radial_support']
    x = np.array([-2., lower, 0.99, 1., 1.01, upper, 1.2])
    metric = CompatibleIncomingMetric(
        (0.001,), (CompatibleRadiusDirection(
            CONTROL.RETARDED.axial_bump(guard['axial_center']), lambda z, n: 0., .007, .03),))
    cache = CachedCompatibleIncomingMetric(x, metric)
    report = cache.verify_against_owned((0., 0.11, 0.15, 0.21, 0.3))
    assert report['value_residual'] < 3e-12
    assert report['direction_residual'] < 3e-12
    assert report['maximum_radius_change_from_background'] > 0
    zero = CompatibleIncomingMetric((0.,), metric.directions)
    changed = np.max(np.abs(cache.values(0.15, x)[3] - zero.values(0.15, x)[3]))
    assert changed > 0
    assert np.max(np.abs(cache.values(0.15, x)[:3] - zero.values(0.15, x)[:3])) < 3e-14
    assert cache.values(0.15, x)[3, 3] == pytest.approx(np.sqrt(2.), abs=3e-14)


def test_rejects_frozen_incoming_c0_and_legacy_matter(monkeypatch):
    owner, metric, times, base, incoming, source = tiny_family()
    tangents = np.zeros((1, *base.shape), complex)
    with pytest.raises(ValueError, match='frozen incoming C0'):
        evolve_incoming_state(owner, times, metric, base, tangents, incoming(1), source,
                              incoming_C0=np.eye(4))
    with pytest.raises(ValueError, match='frozen incoming C0'):
        evolve_incoming_state(owner, times, metric, base, tangents, incoming(1), source,
                              covariance=np.eye(4))
    with pytest.raises(ValueError, match='matter dictionary'):
        evolve_incoming_state(owner, times, metric, base, tangents, incoming(1),
                              {'stress_approximant': np.ones(4), 'parts': {}})
    with pytest.raises(TypeError, match='FixedSourcePreparation'):
        evolve_incoming_state(owner, times, metric, base, tangents, incoming(1), np.eye(2))

    def forbidden(*args, **kwargs):
        raise AssertionError('evolution must not run on a rejected frozen-C0 interface')
    monkeypatch.setattr(STATE, 'evolve_prepared_field_jets', forbidden)
    with pytest.raises(ValueError, match='frozen incoming C0'):
        evolve_incoming_state(owner, times, metric, base, tangents, incoming(1), source, C0=1)


def test_rejects_unsupported_intrinsic_changed_family():
    owner, metric, times, base, incoming, source = tiny_family()
    rho1 = rho1_node_index(owner.x)

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
        evolve_incoming_state(
            owner, times, ChangingIntrinsic(), base, np.zeros((1, *base.shape), complex),
            incoming(1), source)


def test_factory_executes_owned_prepared_evolution(monkeypatch):
    owner, metric, times, base, incoming, source = tiny_family()
    calls = []
    real = STATE.evolve_prepared_field_jets

    def spy(*args, **kwargs):
        calls.append((args, kwargs))
        return real(*args, **kwargs)

    monkeypatch.setattr(STATE, 'evolve_prepared_field_jets', spy)
    state = evolve_incoming_state(
        owner, times, metric, base, np.zeros((1, *base.shape), complex), incoming(1), source)
    assert len(calls) == 1
    assert calls[0][1].get('sample_rows') is not None or len(calls[0][0]) >= 6
    assert state.columns.shape == (len(times), 2, 2)
    assert state.column_tangents.shape == (1, len(times), 2, 2)
    assert state.weighted_columns.shape == state.columns.shape
    assert state.covariance().shape == (len(times), 2, len(times), 2)
    assert state.covariance_tangent().shape == (1, len(times), 2, len(times), 2)
    assert state.history.restriction_parameter_derivative == 0
    assert state.coverage['continuum_error_bound'] is None
    assert state.coverage['physical_history_solution'] == 'OPEN'
    assert state.diagnostics['outgoing_momentum_k_equals_minus_E'] is False
    state.validate()
    state.require_history(metric, times, owner.x)
    other = CompatibleIncomingMetric((0.01,), metric.directions)
    with pytest.raises(ValueError, match='reusable covariance'):
        state.require_history(other, times, owner.x)
    slim = {k: v for k, v in state.diagnostics.items()
            if k not in ('evolved_field', 'tangent_fields')}
    with pytest.raises(ValueError, match='source binding'):
        EvolvedIncomingState(
            z=state.z, columns=state.columns, column_tangents=state.column_tangents,
            source_covariance=np.eye(2, dtype=complex), column_weights=state.column_weights,
            source_energies=state.source_energies, history=state.history,
            coverage=state.coverage, diagnostics=slim)


def test_cross_z_and_all_source_columns_are_retained():
    z = np.array([0.1, 0.2, 0.35])
    columns = np.array([
        [[1., 0.2j], [0.1, 0.7]],
        [[0.3, 0.4], [0.2j, 0.5]],
        [[0.05, 0.9], [0.6, 0.1j]],
    ], dtype=complex)
    tangents = np.array([columns * 0.1, columns * (-0.07 + 0.02j)])
    source = FixedSourcePreparation.from_signed_blocks(
        np.array([0.3]), np.array([0.8]),
        np.array([[[0.5, 0.1 - 0.05j], [0.1 + 0.05j, 0.4]]]))
    history = HistoryBinding(
        amplitudes=(0.001,), inner_radii=(0.007,), outer_radii=(0.03,),
        times=np.array([0., 0.1, 0.2]), spatial_grid=np.array([-2., 1., 1.2]),
        rho1_index=1, rho1_metric=np.array([1., 2., 1., np.sqrt(2.)]),
        restriction=rho1_restriction_map(1., 1., 2., np.sqrt(2.)),
        restriction_parameter_derivative=0., S1=0.05, source_digest=source.digest)
    state = EvolvedIncomingState(
        z=z, columns=columns, column_tangents=tangents,
        source_covariance=source.covariance, column_weights=source.column_weights,
        source_energies=source.energies, history=history,
        coverage=None, diagnostics={'outgoing_momentum_k_equals_minus_E': False})
    kernel = state.covariance()
    assert kernel.shape == (3, 2, 3, 2)
    assert np.max(np.abs(kernel[0, :, 2, :])) > 1e-4
    truncated = EvolvedIncomingState(
        z=z, columns=columns[:, :, :1], column_tangents=tangents[:, :, :, :1],
        source_covariance=source.covariance[:1, :1],
        column_weights=source.column_weights[:1],
        source_energies=source.energies[:1],
        history=HistoryBinding(
            amplitudes=history.amplitudes, inner_radii=history.inner_radii,
            outer_radii=history.outer_radii, times=history.times,
            spatial_grid=history.spatial_grid, rho1_index=history.rho1_index,
            rho1_metric=history.rho1_metric, restriction=history.restriction,
            restriction_parameter_derivative=0., S1=history.S1,
            source_digest=FixedSourcePreparation(
                source.covariance[:1, :1], source.column_weights[:1],
                source.energies[:1]).digest),
        coverage=None, diagnostics={'outgoing_momentum_k_equals_minus_E': False})
    full = kernel.reshape(6, 6)
    part = truncated.covariance().reshape(6, 6)
    assert np.max(np.abs(full - part)) > 1e-4
    F = state.weighted_columns.reshape(6, 2)
    dF = state.weighted_column_tangents.reshape(2, 6, 2)
    zero = np.zeros((2, 2, 2), complex)
    expected = restriction_covariance_tangent(F, dF, state.source_covariance, zero)
    assert np.max(np.abs(state.covariance_tangent() - expected.reshape(2, 3, 2, 3, 2))) < 1e-15
    fake = zero.copy()
    fake[0, 0, 1] = 0.2
    fake[0, 1, 0] = 0.2
    wrong = restriction_covariance_tangent(F, dF, state.source_covariance, fake)
    assert np.max(np.abs(wrong - expected)) > 1e-4
    frozen = restriction_covariance_tangent(F, np.zeros_like(dF), state.source_covariance, zero)
    assert np.max(np.abs(frozen - expected)) > 1e-4
    assert np.max(np.abs(state.covariance() - state.covariance().transpose(2, 3, 0, 1).conj())) < 3e-15


def test_matching_a_reference_covariance_is_not_a_frozen_input_failure():
    source = tiny_source()
    z = np.array([0.2, 0.4])
    columns = np.zeros((2, 2, 2), complex)
    columns[0, 0, 0] = 1.
    columns[1, 1, 1] = 1.
    tangents = np.zeros((1, 2, 2, 2), complex)
    history = HistoryBinding(
        amplitudes=(0.,), inner_radii=(0.007,), outer_radii=(0.03,),
        times=np.array([0., 0.1]), spatial_grid=np.array([-2., 1., 1.2]),
        rho1_index=1, rho1_metric=np.array([1., 2., 1., np.sqrt(2.)]),
        restriction=rho1_restriction_map(1., 1., 2., np.sqrt(2.)),
        restriction_parameter_derivative=0., S1=0.1, source_digest=source.digest)
    state = EvolvedIncomingState(
        z=z, columns=columns, column_tangents=tangents,
        source_covariance=source.covariance, column_weights=source.column_weights,
        source_energies=source.energies, history=history,
        coverage=None, diagnostics={'outgoing_momentum_k_equals_minus_E': False})
    reference = state.covariance()
    assert np.max(np.abs(reference)) > 0
    assert state.validate() is True


def test_source_energies_are_labels_not_minus_k():
    owner, metric, times, base, incoming, source = tiny_family()
    state = evolve_incoming_state(
        owner, times, metric, base, np.zeros((1, *base.shape), complex), incoming(1), source)
    F = state.weighted_columns.reshape(len(times) * 2, 2)
    manual = (F @ state.source_covariance @ F.conj().T).reshape(len(times), 2, len(times), 2)
    assert np.max(np.abs(manual - state.covariance())) < 1e-14
    assert 'momentum' not in state.diagnostics
    assert state.diagnostics['outgoing_momentum_k_equals_minus_E'] is False
    assert np.array_equal(state.source_energies, source.energies)


def test_tiny_centered_derivatives_and_primal_amplitude_only():
    owner, metric, times, base, incoming, source = tiny_family(0.003)
    step = 1e-4
    center = evolve_incoming_state(
        owner, times, metric, base, np.zeros((1, *base.shape), complex), incoming(1), source)
    plus = evolve_incoming_state(
        owner, times, AmplitudeOnlyMetric(CompatibleIncomingMetric((0.003 + step,), metric.directions)),
        base, np.zeros((0, *base.shape), complex), incoming(0), source)
    minus = evolve_incoming_state(
        owner, times, AmplitudeOnlyMetric(CompatibleIncomingMetric((0.003 - step,), metric.directions)),
        base, np.zeros((0, *base.shape), complex), incoming(0), source)
    field_fd = (plus.diagnostics['evolved_field'] - minus.diagnostics['evolved_field']) / (2 * step)
    trace_fd = (plus.columns - minus.columns) / (2 * step)
    cov_fd = (plus.covariance() - minus.covariance()) / (2 * step)
    assert np.max(np.abs(field_fd - center.diagnostics['tangent_fields'][0])) < 3e-8
    assert np.max(np.abs(trace_fd - center.column_tangents[0])) < 3e-8
    assert np.max(np.abs(cov_fd - center.covariance_tangent()[0])) < 3e-8
    assert center.diagnostics['norm_flux_residual'] < 3e-11


def test_array_only_saved_retarded_restriction_without_old_solver():
    residual = CONTROL.retarded_restriction_replay(
        json.loads((ROOT / 'results/development/nsc-retarded-compatible-response.json').read_text()),
        CONTROL.load_retarded_inputs()[2])
    assert residual['old_solver_rerun'] is False
    assert residual['final_trace_residual'] < 3e-12
    assert residual['final_trace_tangent_residual'] < 3e-12
    assert residual['sampled_kernel_residual'] < 3e-12
    assert residual['restriction_parameter_derivative'] == 0


def test_receipt_replays_arrays_without_propagation(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('scientific producer called during receipt replay')

    monkeypatch.setattr(CONTROL, 'one_state', forbidden)
    monkeypatch.setattr(CONTROL, 'evolve_incoming_state', forbidden)
    monkeypatch.setattr(STATE, 'evolve_prepared_field_jets', forbidden)
    monkeypatch.setattr(CONTROL.RETARDED, 'one_case', forbidden)
    monkeypatch.setattr(CONTROL.RETARDED, 'run', forbidden)
    record = CONTROL.check()
    assert record['control']['field_runs'] == 3
    assert record['control']['field_run_cap'] == 3
    assert record['control']['amplitude'] == 0.001
    assert record['control']['outgoing_momentum_substitution'] is None
    assert record['scope']['continuum_error_bound'] is None
    assert record['scope']['physical_history_solution'] == 'OPEN'
    assert record['scope']['cosmology_EXISTENCE'] == 'OPEN'
    assert record['scope']['frozen_incoming_C0'] is False
    assert record['residuals']['complete_field_derivative'] < record['control']['field_derivative_tolerance']
    assert record['residuals']['complete_coherent_covariance_derivative'] < record['control']['covariance_derivative_tolerance']
    assert record['algebra_maximum'] < record['control']['flux_algebra_tolerance']
    assert min(record['omitted_term_effects'].values()) > 1e-16
    with np.load(ROOT / record['payload']['path'], allow_pickle=False) as arrays:
        C = arrays['source/covariance']
        weights = arrays['source/column_weights']
        F = (arrays['center/columns'] * weights).reshape(-1, C.shape[0])
        dF = (arrays['center/column_tangents'] * weights).reshape(-1, C.shape[0])
        expected = dF @ C @ F.conj().T + F @ C @ dF.conj().T
        saved = arrays['center/covariance_tangent'].reshape(F.shape[0], F.shape[0])
        assert np.max(np.abs(expected - saved)) < 1e-15
        assert arrays['source/blocks'].shape[0] == 4
        assert arrays['center/columns'].shape[-1] == 12
