"""Whole-realization translations, wrapped measurements and AP spin structure."""
from pathlib import Path
import json

import numpy as np
import pytest

import derive_nsc_discovery_translation as cli
from recursive_horizons import nsc_discovery_translation as translation
from recursive_horizons import nsc_discovery_episode as episode
from recursive_horizons import nsc_discovery_backend as backend
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


def test_half_integer_shift_has_antiperiodic_monodromy_and_no_naive_wrap():
    count = 64
    x = np.arange(count) * 8. / count
    values = np.exp(2j * np.pi * 2.5 * x / 8.) + .2 * np.exp(-2j * np.pi * 3.5 * x / 8.)
    delta = .37
    expected = np.exp(2j * np.pi * 2.5 * (x-delta) / 8.) + .2 * np.exp(-2j * np.pi * 3.5 * (x-delta) / 8.)
    assert translation.antiperiodic_shift(values, delta) == pytest.approx(expected, abs=1e-14)
    assert translation.antiperiodic_shift(values, 8.) == pytest.approx(-values, abs=1e-14)
    assert translation.antiperiodic_shift(translation.antiperiodic_shift(values, delta), -delta) == pytest.approx(values, abs=1e-14)
    moved = translation.antiperiodic_shift(values, 5.5)
    assert np.max(np.abs(moved - np.roll(values, 44))) > 1.
    assert np.linalg.norm(moved) == pytest.approx(np.linalg.norm(values), abs=1e-14)


def test_integer_geometry_translation_and_AP_prolongation_commute():
    count = 63
    x = np.arange(count) * 8. / count
    values = 1. + .1*np.cos(6*np.pi*x/8.)
    expected = 1. + .1*np.cos(6*np.pi*(x-.37)/8.)
    assert translation.periodic_shift(values, .37) == pytest.approx(expected, abs=1e-14)
    grid = backend.make_fft_grid(galerkin.build_grid(64, gauge='conformal'))
    rng = np.random.default_rng(21)
    columns = rng.normal(size=(64, 6)) + 1j*rng.normal(size=(64, 6))
    left = galerkin.prolong_columns(grid, translation.antiperiodic_shift(columns, 5.5))
    right = translation.antiperiodic_shift(galerkin.prolong_columns(grid, columns), 5.5)
    assert left == pytest.approx(right, abs=2e-14)


def test_wrapped_window_keeps_span_and_translated_integral():
    grid = backend.make_fft_grid(galerkin.build_grid(64, gauge='conformal'))
    values = 2.+np.cos(2*np.pi*grid.xi_q/8.)
    moved = translation.periodic_shift(values, 5.5)
    assert translation.interval_pieces((1.,3.),5.5) == [(6.5,8.),(0.,.5)]
    expected = translation.translated_integral(grid, values, (1.,3.),0.)
    assert translation.translated_integral(grid,moved,(1.,3.),5.5) == pytest.approx(expected, abs=1e-13)
    assert sum(b-a for a,b in translation.interval_pieces((0.,4.),5.5)) == 4.


@pytest.fixture(scope='module')
def saved_report():
    return translation.assessment()


def test_saved_whole_rate_force_constraint_tidal_and_observer_covariance(saved_report):
    assert saved_report['ok'] is True
    assert saved_report['hashes_unchanged'] is True
    assert saved_report['trajectory_evolved'] is False
    assert saved_report['physical_mechanism_claimed_from_identity'] is False
    assert saved_report['cpu_seconds'] < 30.
    assert {(case['nf'],case['time']) for case in saved_report['cases']} == {(128,.3),(128,1.),(256,.3),(256,1.)}
    for case in saved_report['cases']:
        assert case['input_state_basis_source_preserved'] is True
        assert case['translated']['clock_locations'] == [6.5,7.5,.5]
        assert case['translated']['windows']['child']['pieces'] == [(6.5,8.),(0.,.5)]
        assert len(case['saved_proper_clocks_at_shifted_worldlines']) == 3
        assert all(item['within_roundoff_tolerance'] for item in case['checks'].values())
        assert {'force_Q','force_L','force_beta','constraint_hamilton','constraint_momentum','actual_curvature_R_0101','actual_curvature_R_0202','actual_curvature_R4','actual_curvature_owned_W','diagnostic_reference_amplitudes','diagnostic_source_amplitudes'} <= set(case['checks'])


def test_optional_one_step_and_translated_basis_interpretation():
    case_id = 'nf128_coupled_dt0.0005'
    record, arrays = episode.load_checkpoint(translation.DEFAULT_EPISODE,case_id,0)
    pair = backend.make_fft_pair(episode.pair_from_arrays(arrays,record))
    state = episode.state_from_arrays(arrays,record['momentum_representation'])
    moved, changed = translation.translate_pair_state(pair,state,5.5)
    assert moved.geometry_child_indices is pair.geometry_child_indices
    assert moved.geometry_parent_indices is pair.geometry_parent_indices
    assert not moved.geometry_map.flags.writeable
    for name in ('Q','r','chi','p_Q','p_r','p_chi'):
        assert np.array_equal(getattr(changed,name),getattr(state,name))
    compared = translation.compare_translation(pair,state,step=.0005)
    assert compared['ok'] is True
    assert all(compared['checks']['one_step_'+name]['within_roundoff_tolerance'] for name in ('Q','r','chi','p_Q','p_r','p_chi','phi0','phi1'))
    with pytest.raises(ValueError,match='fermion lattice'):
        translation.translate_pair_state(pair,state,.37)


def test_creation_only_cli_refuses_before_assessment(tmp_path, monkeypatch, saved_report, capsys):
    monkeypatch.setattr(cli,'OUTPUT_DIR',tmp_path)
    monkeypatch.setattr(translation,'assessment',lambda **kwargs: saved_report)
    destination=tmp_path/'new.json'
    assert cli.main(['--write',str(destination)]) == 0
    assert json.loads(destination.read_text())['ok'] is True
    monkeypatch.setattr(translation,'assessment',lambda **kwargs: pytest.fail('assessment started before destination validation'))
    with pytest.raises(FileExistsError):
        cli.main(['--write',str(destination)])
    with pytest.raises(PermissionError):
        cli.main(['--write',str(tmp_path.parent/'outside.json')])
