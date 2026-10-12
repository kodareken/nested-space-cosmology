"""Endpoint identities and scope; never replace a sampled residual by a bound."""
import importlib.util
import json
from pathlib import Path

import mpmath as mp
import pytest

import recursive_horizons.nsc_incoming_horizon_endpoint_bound as owner
from recursive_horizons.nsc_incoming_source_quadrature_bound import unpack
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi


ROOT = Path(__file__).resolve().parents[1]


def inputs():
    channel = json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels'][32]
    config = json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    return channel, config


def test_finite_frobenius_residual_has_only_owned_last_term():
    assert owner.frobenius_residual_identity()['coefficient_residuals'] == ['0', '0']
    with mp.workdps(80):
        z, nu = mp.mpf('1.2'), mp.mpf('3.4')
        def terms(x, shift):
            term = mp.mpc(1); values = [term]
            for j in range(1, 10):
                term *= -x*x/(4*j*(shift+j-1)); values.append(term)
            return values
        a = lambda x: sum(terms(x, mp.mpf('.5')+1j*nu))
        b = lambda x: -1j*x/(1+2j*nu)*sum(terms(x, mp.mpf('1.5')+1j*nu))
        last = -1j*z/(1+2j*nu)*terms(z, mp.mpf('1.5')+1j*nu)[-1]
        assert abs(z/2*mp.diff(a,z)+1j*z/2*b(z)-1j*z/2*last) < mp.mpf('1e-76')
        assert abs(z/2*mp.diff(b,z)+1j*z/2*a(z)+1j*nu*b(z)) < mp.mpf('1e-76')


def test_exact_profile_coordinate_and_frozen_source_kappa_are_separate():
    c, cfg = inputs(); result = owner.endpoint_bound(c, cfg)
    b = result['bounds']
    with _precision(80):
        assert _lo(unpack(b['q_coordinate_shift']['binary_interval'])) > mp.mpf('1.49e-15')
        assert _hi(unpack(b['geometric_minus_source_kappa']['binary_interval'])) < 0
        frozen = unpack(b['frozen_source_kappa']['binary_interval'])
        assert _lo(frozen) == _hi(frozen) == mp.mpf(cfg['surface_gravity'])
        assert _lo(unpack(b['true_distance_at_stored_endpoint']['binary_interval'])) > _hi(unpack(b['stored_delta']['binary_interval']))
        assert _hi(unpack(b['finite_series_norm_defect_upper']['binary_interval'])) < mp.mpf('2e-101')
    assert not result['mathematical_endpoint_budget_pass']
    assert not result['aligned_coordinate_budget_pass']
    assert result['scope']['archived_floating_initializer_error_bound'] is None
    assert result['scope']['archived_DOP853_propagation_error_bound'] is None
    with pytest.raises(ValueError, match='group32'):
        owner.endpoint_bound(dict(c,index=31),cfg)


def test_covariance_phase_removal_matches_actual_working_frame():
    c, cfg = inputs(); values = owner.archived_formula_controls(c, cfg)
    assert values['samples'] == 6
    assert values['sampled_working_frame_covariance_difference'] < 3e-14
    assert values['sampled_raw_frame_norm_residual'] < 3e-14
    assert 'not uniform' in values['interpretation']


def test_receipt_replay_never_recomputes_horizon_or_endpoint(monkeypatch):
    path = ROOT/'results/development/nsc-incoming-horizon-endpoint-bound.json'
    if not path.exists():
        pytest.skip('endpoint record follows the primitive tests')
    old = json.loads(path.read_text())
    def forbidden(*args, **kwargs):
        raise AssertionError('endpoint generation called during replay')
    monkeypatch.setattr(owner, 'endpoint_bound', forbidden)
    monkeypatch.setattr(owner, 'archived_formula_controls', forbidden)
    spec = importlib.util.spec_from_file_location('endpoint_replay', ROOT/'scripts/derive_nsc_incoming_horizon_endpoint_bound.py')
    replay = importlib.util.module_from_spec(spec); spec.loader.exec_module(replay)
    assert replay.make_record(ROOT/old['payload']['path']) == old
