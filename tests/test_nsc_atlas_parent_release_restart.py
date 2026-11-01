"""Restart/export failure regression controls with no physical evolution."""
from argparse import Namespace
import json
from types import SimpleNamespace

import numpy as np
import pytest

from scripts import continue_nsc_atlas_parent_release as driver
from recursive_horizons.nsc_atlas_gh_history import DenseGHHistory


@pytest.fixture
def fake_run(monkeypatch, tmp_path):
    class System:
        dof, r, epsilon, gauge, backend = 2, 1., 1., 'lapse_controlled', 'hopf'
        calls = 0
        basis = SimpleNamespace(quadrature_shape=(2, 4, 4))

        def initial(self):
            pytest.fail('a restart must never call the initializer')

        def observer_angle(self, state):
            return 0.

        def rhs(self, t, state):
            assert t >= 12.3  # The original absolute coordinate clock.
            self.calls += 1
            return np.zeros_like(state)

    system = System()
    state0 = np.r_[[123., 456.], [0., 0., 0., 1.], [1., 0., 0., 0.], .9, 0.]
    data = dict(system=system, end=12.3,
        history=SimpleNamespace(evaluate=lambda t: state0.copy(), metadata={'quadrature_order': [2, 4]}),
        record={'model_inputs': {'unchanged': True}},
        provenance={'dense_cache': 'synthetic-original-history.npz', 'auxiliary_phase': {'values': [dict(reference_observer_proper_duration=.9,
            sparse_sampled_phase=37.25)]}, 'source_sha256': {}})

    class Dense:
        def __init__(self, old_t, t, old, new):
            self.t_old, self.t, self.y_old = old_t, t, old.copy()
            self.F = np.zeros((7, len(old)))
            self.F[0] = new-old

        def __call__(self, t):
            return self.y_old+(t-self.t_old)/(self.t-self.t_old)*self.F[0]

    class Solver:
        def __init__(self, function, t, state, **kwargs):
            self.function, self.t, self.y = function, t, state.copy()
            self.status, self.t_old = 'running', None
            function(t, self.y)

        def step(self):
            self.function(self.t, self.y)  # Can fail before acceptance.
            self.t_old, self.old = self.t, self.y.copy()
            self.t += .12
            self.y[0] += 7.
            self.y[system.dof+8] += .12

        def dense_output(self):
            return Dense(self.t_old, self.t, self.old, self.y)

    def readout(system, t, state, phase):
        return dict(coordinate_t=float(t), observer_proper_duration=float(state[system.dof+8]),
            s=np.float64(-.035), uS=np.float64(.001), ADM_H_RMS_over_H_squared=.001,
            GH_C_RMS_over_H=.001, numerical_guard_passed=np.bool_(True),
            passive_sampled_phase=phase)

    monkeypatch.setattr(driver, 'ROOT', tmp_path)
    monkeypatch.setattr(driver, 'set_num_threads', lambda n: None)
    monkeypatch.setattr(driver, 'load_original_parent', lambda degree: data)
    monkeypatch.setattr(driver, 'source_identity', lambda: ('synthetic-test-commit', {}))
    monkeypatch.setattr(driver, 'DOP853', Solver)
    monkeypatch.setattr(driver, 'snapshot', readout)
    monkeypatch.setattr(driver, 'anchor', lambda system, state: {'s': -.035, 'uS': .001})
    args = Namespace(degree=20, output='.build/fake', threads=1, budget=1000.,
        proper_endpoint=1.1, rtol=1e-9, atol=1e-11, max_step=.01)
    return SimpleNamespace(args=args, system=system, state0=state0, folder=tmp_path/args.output,
        snapshot=readout)


def test_numpy_checkpoint_flag_and_exact_restart_are_preserved(fake_run):
    result = driver.run(fake_run.args)
    assert result['status'] == 'event_reached'
    assert [row['event'] for row in result['events']] == ['proper_checkpoint', 'proper_endpoint']
    assert all(type(row['terminal']) is bool for row in result['events'])
    with np.load(fake_run.folder/'states.npz', allow_pickle=False) as raw:
        np.testing.assert_array_equal(raw['states'][0], fake_run.state0)
        assert raw['times'][0] == 12.3
        assert raw['accepted_endpoint_coordinate_t'] > result['physical_stop_coordinate_t']
        assert np.all(raw['states'][:, -1] == 0.)  # Phase cannot feed the ODE.
    assert result['auxiliary_phase']['sparse_phase_carry'] == 37.25
    assert result['snapshots'][0]['passive_sampled_phase'] == 37.25
    json.loads((fake_run.folder/'run.json').read_text(encoding='utf-8'))


def test_budget_midstep_keeps_last_accepted_state_and_polynomial(fake_run, monkeypatch):
    original = fake_run.system.rhs
    def budget_failure(t, state):
        if fake_run.system.calls >= 2:
            raise TimeoutError('synthetic budget interruption before next acceptance')
        return original(t, state)
    monkeypatch.setattr(fake_run.system, 'rhs', budget_failure)
    result = driver.run(fake_run.args)
    assert result['status'] == 'budget_stop'
    with np.load(fake_run.folder/'states.npz', allow_pickle=False) as raw:
        assert raw['accepted_endpoint'][0] == fake_run.state0[0]+7
        assert raw['accepted_endpoint_coordinate_t'] == pytest.approx(12.42)
    history = DenseGHHistory.load(fake_run.folder/'history.npz')
    assert len(history.bounds) == 1
    assert history.evaluate(history.bounds[-1][1])[0] == fake_run.state0[0]+7


def test_public_json_error_cannot_prevent_raw_flush(fake_run, monkeypatch):
    dumps = json.dumps
    def failing_public_export(value, *args, **kwargs):
        if (isinstance(value, dict) and value.get('schema') == 'nsc-atlas-parent-release-continuation-v1'
                and value.get('events')):
            raise TypeError('synthetic public record export failure')
        return dumps(value, *args, **kwargs)
    monkeypatch.setattr(driver.json, 'dumps', failing_public_export)
    original_readout = fake_run.snapshot
    def raw_first_readout(system, t, state, phase):
        if state[system.dof+8] > 1.005:
            assert (fake_run.folder/'states.npz').is_file()
            assert (fake_run.folder/'history.npz').is_file()
        return original_readout(system, t, state, phase)
    monkeypatch.setattr(driver, 'snapshot', raw_first_readout)
    result = driver.run(fake_run.args)
    assert result['status'] == 'export_failure'
    with np.load(fake_run.folder/'states.npz', allow_pickle=False) as raw:
        assert raw['accepted_endpoint'][0] == fake_run.state0[0]+7
        np.testing.assert_array_equal(raw['states'][0], fake_run.state0)
    terminal = json.loads((fake_run.folder/'run.json').read_text(encoding='utf-8'))
    assert terminal['status'] == 'export_failure'
    assert terminal['state_cache_sha256'] == driver.sha(fake_run.folder/'states.npz')
    assert terminal['dense_history']['sha256'] == driver.sha(fake_run.folder/'history.npz')


@pytest.mark.parametrize('failure', ['type', 'nonfinite'])
def test_terminal_diagnostic_failures_happen_after_raw_preservation(fake_run, monkeypatch, failure):
    rhs = fake_run.system.rhs
    def budget_failure(t, state):
        if fake_run.system.calls >= 2:
            raise TimeoutError('synthetic budget')
        return rhs(t, state)
    monkeypatch.setattr(fake_run.system, 'rhs', budget_failure)
    original = fake_run.snapshot
    def diagnostic(system, t, state, phase):
        if state[system.dof+8] > 1.005:
            assert (fake_run.folder/'states.npz').is_file()
            assert (fake_run.folder/'history.npz').is_file()
            if failure == 'type':
                raise TypeError('synthetic terminal diagnostic type failure')
            row = original(system, t, state, phase)
            row['s'] = float('nan')
            return row
        return original(system, t, state, phase)
    monkeypatch.setattr(driver, 'snapshot', diagnostic)
    result = driver.run(fake_run.args)
    assert result['status'] == 'budget_stop'
    assert 'accepted_endpoint_diagnostic_failure' in result
    with np.load(fake_run.folder/'states.npz', allow_pickle=False) as raw:
        assert raw['accepted_endpoint'][0] == fake_run.state0[0]+7
    terminal = json.loads((fake_run.folder/'run.json').read_text(encoding='utf-8'))
    assert terminal['status'] == 'budget_stop'


def test_json_numbers_are_native_but_unsupported_data_stays_an_error():
    native = json.loads(json.dumps({'flag': np.bool_(False), 'count': np.int64(2)},
                                   default=driver.json_scalar, allow_nan=False))
    assert native == {'flag': False, 'count': 2}
    with pytest.raises(TypeError, match='unsupported JSON value'):
        json.dumps({'unsupported': object()}, default=driver.json_scalar)
    with pytest.raises(ValueError):
        json.dumps({'nonfinite': np.float64(np.nan)}, default=driver.json_scalar, allow_nan=False)
