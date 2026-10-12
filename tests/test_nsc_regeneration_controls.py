"""Controls for changed-source data, frozen geometry, and the saved continuation.

The long episodes are read back. These tests do not start another production run.
"""
import json
import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import numpy as np
import pytest

import derive_nsc_regeneration_controls as regeneration
import derive_nsc_spherical_feedback_episode as episode
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons.nsc_regeneration_controls import (
    BASELINE_OCCUPATIONS,
    CONTENT_EFFECT,
    CRITERION,
    FINE_IDENTITY_IS_ALGEBRAIC_CONTROL,
    REVERSED_OCCUPATIONS,
    UNIFORM_OCCUPATIONS,
    assess_windows,
    bracket_g,
    car_report,
    chart_from_fields,
    continuation_check,
    dynamic_stop_reason,
    frozen_rk4_step,
    galerkin_half_identity_report,
    geometry_bytes,
    load_episode_final,
    prepare_population,
    reversal_amplitude,
    separate_ledgers,
    solve_source_radius,
    source_arrays,
    state_sha256,
    uniform_complement_localization,
    with_occupations,
)
from recursive_horizons.nsc_spherical_coupling import chart_failure

OUT = episode.LAB / "results" / "development" / "nsc-regeneration-controls-v1.json"


def _packet(fermions=64):
    phi0, phi1, _preparation, problems = galerkin.load_physical_columns(fermions)
    assert problems == []
    grid = galerkin.build_grid(fermions)
    state = galerkin.blank_state(grid, phi0, phi1)
    return grid, state


def test_source_radius_matches_the_owner_on_the_baseline_packet():
    grid, state = _packet()
    owned, info = galerkin.solve_initial_radius(grid, state)
    solved, report = solve_source_radius(grid, state)
    assert info["converged"] is True
    assert report["converged"] is True
    assert report["physical_failure"] is False
    assert report["source_mean_subtracted"] is False
    assert np.array_equal(solved.r, owned.r)
    assert report["positive_G"] is True
    assert report["r_min"] > 0.0


def test_changed_occupations_recompute_the_source_and_keep_the_mean():
    grid, state = _packet()
    original = np.array(grid.fine.occupations, copy=True)
    caller_phi = np.array(state.phi0, copy=True)
    uniform_grid, uniform, uniform_report = prepare_population(grid, state, UNIFORM_OCCUPATIONS)
    reversed_grid, reversed_state, reversed_report = prepare_population(grid, state, REVERSED_OCCUPATIONS)
    assert np.array_equal(grid.fine.occupations, original)
    assert np.array_equal(state.phi0, caller_phi)
    assert np.allclose(uniform_grid.fine.occupations, UNIFORM_OCCUPATIONS)
    assert np.allclose(reversed_grid.fine.occupations, REVERSED_OCCUPATIONS)
    for report, solved in ((uniform_report, uniform), (reversed_report, reversed_state)):
        assert report["columns_unchanged"] is True
        assert report["current_unchanged"] is True
        assert report["source_mean_subtracted"] is False
        assert report["mean_retained_in_shift_residual"] is True
        assert report["shift_mean_gap"] < 1e-8
        assert report["physical_failure"] is False
        assert report["positive_G"] is True
        assert report["positive_r"] is True
        assert report["positive_Q"] is True
        assert report["car"]["admissible"] is True
        assert report["projected_hamilton_max"] is not None
        assert report["full_hamilton_max"] is not None
        assert report["evolve"] is False
        assert solved.phi0.shape == state.phi0.shape
    assert uniform_report["full_hamilton_max"] > 1e-2
    assert uniform_report["rho_mean"] != reversed_report["rho_mean"]
    assert car_report(uniform)["gap"] <= 1e-8
    assert np.allclose(original, BASELINE_OCCUPATIONS)


def test_uniform_source_differs_from_the_baseline_source():
    grid, state = _packet()
    from recursive_horizons.nsc_regeneration_controls import source_arrays
    _source, rho_base, current_base = source_arrays(grid, state)
    owned, solved, report = prepare_population(grid, state, UNIFORM_OCCUPATIONS)
    _other, rho_new, current_new = source_arrays(owned, solved)
    assert np.max(np.abs(rho_new - rho_base)) > 0.0
    assert report["current_unchanged"] is True
    assert np.array_equal(np.asarray(current_new), np.asarray(current_new))
    del current_base


def test_frozen_step_keeps_geometry_and_applies_no_fieldwork():
    grid, state = _packet()
    solved, _report = solve_source_radius(grid, state)
    before = geometry_bytes(solved)
    phi = np.array(solved.phi0, copy=True)
    updated, refused = frozen_rk4_step(grid, solved, 5e-4)
    assert geometry_bytes(updated) == before
    assert geometry_bytes(solved) == before
    assert not np.array_equal(updated.phi0, phi)
    assert np.isfinite(refused)


def test_episode_final_loader_is_bitwise():
    state = load_episode_final(episode.NPZ, "nf512_dt_0_0005")
    with np.load(episode.NPZ, allow_pickle=False) as data:
        assert np.array_equal(state.Q, data["nf512_dt_0_0005_final_Q"])
        assert np.array_equal(state.r, data["nf512_dt_0_0005_final_r"])
        assert np.array_equal(state.chi, data["nf512_dt_0_0005_final_chi"])
        assert np.array_equal(state.p_Q, data["nf512_dt_0_0005_final_p_Q"])
        assert np.array_equal(state.p_r, data["nf512_dt_0_0005_final_p_r"])
        assert np.array_equal(state.p_chi, data["nf512_dt_0_0005_final_p_chi"])
        assert np.array_equal(state.phi0, data["nf512_dt_0_0005_final_phi0"])
        assert np.array_equal(state.phi1, data["nf512_dt_0_0005_final_phi1"])
    again = load_episode_final(episode.NPZ, "nf512_dt_0_0005")
    assert state_sha256(state) == state_sha256(again)


def test_two_segments_of_one_drift_are_not_renewal():
    content = np.linspace(6.0, 6.22, 21)
    shares_scale = np.array([0.60, 0.34, 0.03, 0.03])
    samples = content[:, None] * shares_scale[None, :]
    flux = np.tile(np.array([0.05, -3.3, 3.1, 0.15]), (21, 1))
    whole = assess_windows(samples, flux)
    second = assess_windows(samples[10:], flux[10:])
    assert whole["renewed"] is False
    assert whole["one_drift"] is True
    assert whole["maintained"] is True
    assert second["renewed"] is False
    assert reversal_amplitude(content) < CONTENT_EFFECT


def test_a_localized_reversal_is_renewal():
    up = np.linspace(6.0, 6.25, 11)
    down = np.linspace(6.25, 6.05, 11)
    leader = np.concatenate([up, down[1:]])
    rest = np.full(leader.shape, 4.0)
    samples = np.column_stack([leader, rest * 0.8, rest * 0.05, rest * 0.05])
    flux = np.zeros_like(samples)
    verdict = assess_windows(samples, flux)
    assert verdict["renewed"] is True
    assert verdict["localized_end"] is True
    assert verdict["reversal"] >= CONTENT_EFFECT


def test_saved_campaign_record_matches_the_continuation_contract():
    if not OUT.is_file():
        pytest.fail("campaign record is missing")
    record = json.loads(OUT.read_text())
    assert record["schema"] == "NSC-REGENERATION-CONTROLS-v1"
    assert record["mean_source_subtracted"] is False
    assert record["inputs_preserved"] is True
    assert record["episode_driver_modified"] is False
    continuation = record["stages"]["continuation_nf512"]
    state = load_episode_final(episode.NPZ, "nf512_dt_0_0005")
    assert continuation["start_sha256"] == state_sha256(state)
    assert continuation["mean_source_subtracted"] is False
    frozen = record["stages"]["frozen_nf512"]
    assert frozen["geometry_unchanged"] is True
    assert frozen["applied_fieldwork_is_zero"] is True
    uniform = record["initial_data"]["nf512_uniform"]
    reversed_case = record["initial_data"]["nf512_reversed"]
    for report in (uniform, reversed_case):
        assert report["current_unchanged"] is True
        assert report["source_mean_subtracted"] is False
        assert report["shift_mean_gap"] < 1e-8
        assert report["positive_G"] and report["positive_r"] and report["positive_Q"]
        assert report["car"]["admissible"] is True
        assert report["physical_failure"] is False
        assert report["columns_unchanged"] is True
        assert report["rho_gap_from_baseline"] > 0.0
    conclusion = record["conclusion"]
    assert conclusion["splitting_one_drift_is_renewal"] is False
    if not conclusion["candidate_regime"]:
        assert conclusion["perturbations_executed"] is False
    balance = continuation["balance"]
    assert "versus_field_change" in balance
    assert "versus_applied_work" in balance
    assert "field_change_over_saved_exchange" in balance


def _small_packet():
    grid = galerkin.build_grid(14, quadrature=64)
    phi0, phi1, _preparation, problems = galerkin.load_physical_columns(14)
    assert problems == []
    state = galerkin.blank_state(grid, phi0, phi1)
    return grid, state


def test_negative_g_dynamic_chart_continues_and_checkpoints_a_real_exit():
    grid, state = _small_packet()
    fine = galerkin.prolong_state(grid, state)
    _source, rho, _current = source_arrays(grid, state)
    chart = chart_from_fields(grid, {"Q": fine.Q, "rho": rho, "r": fine.r})
    assert chart_failure(grid.fine, fine) is None
    assert chart["positive_G"] is False
    assert chart["G_min"] < 0.0
    assert chart["positive_r"] and chart["positive_Q"] and chart["positive_L"] and chart["positive_lapse"]
    assert chart["g_is_dynamic_stop"] is False
    assert dynamic_stop_reason(None, car_report(state)) == (None, False)
    assert "min G" in CRITERION["chart_stop"]
    assert "when G, r, and Q are positive" in CRITERION["evolution_gate"]
    checked = continuation_check(grid, state, 1.0e-4)
    assert checked["accepted"] is True
    assert checked["chart_stop"] is False
    assert checked["stop_reason"] is None
    assert checked["positive_G"] is False
    assert checked["Q_change"] > 0.0
    assert checked["gram"]["admissible"] is True
    assert checked["positive_r"] and checked["positive_Q"] and checked["positive_L"]
    result = regeneration.integrate(
        grid, state, 1.0e-4, 1.0e-4, 30.0, lambda: 0.0,
        t0=0.0, frozen=False, stop_on_delocalization=False, label="negative-g",
    )
    assert result["completed"] is True
    assert result["stop_reason"] is None
    assert result["chart_stop"] is False
    assert result["g_used_as_dynamic_stop"] is False
    assert result["positive_G"] is False
    assert result["positive_r"] and result["positive_Q"] and result["positive_L"]
    assert result["steps_completed"] == 1
    assert result["Q_change_max"] > 0.0
    assert result["gram_gap_max"] <= 1e-8
    assert result["error"] is None
    assert float(np.min(result["series"]["G_min"])) < 0.0
    assert result["checkpoints"][-1]["positive_G"] is False
    assert result["checkpoints"][-1]["Q_change"] > 0.0
    exited = state.copy()
    exited.Q = np.full_like(exited.Q, -0.1)
    refused = regeneration.integrate(
        grid, exited, 1.0e-4, 1.0e-4, 30.0, lambda: 0.0,
        t0=0.0, frozen=False, stop_on_delocalization=False, label="chart-exit",
    )
    assert refused["chart_stop"] is True
    assert refused["stop_reason"] == "Q_left_positive_chart"
    assert refused["steps_completed"] == 0
    assert refused["error"] is None
    assert refused["checkpoints"][0]["exit_reason"] == "Q_left_positive_chart"
    assert refused["checkpoints"][0]["kept"] is False
    gee = bracket_g(fine.Q, rho, grid)
    assert float(np.min(gee)) < 0.0


def test_galerkin_half_identity_commutes_and_the_fine_embedding_need_not():
    grid, _state = _small_packet()
    report = galerkin_half_identity_report(grid)
    assert report["dimension_H_G"] == 2 * 14
    assert report["dimension_H_fine"] == 2 * 64
    assert report["pullback_gap"] < 1e-8
    assert report["half_identity_commutes_with_H_G"] is True
    assert report["half_identity_commutator_H_G"] < 1e-8
    assert report["embedding_commutes_with_unprojected_H_fine"] is False
    assert report["embedding_commutator_H_fine"] > 1.0
    assert report["fine_identity_is_algebraic_control"] is True
    assert report["fine_identity_is_an_occupation"] is False
    assert FINE_IDENTITY_IS_ALGEBRAIC_CONTROL is True
    owned = with_occupations(grid, np.full(6, 0.5))
    assert np.array_equal(owned.fine.occupations, np.full(6, 0.5))
    with pytest.raises(ValueError, match="six"):
        with_occupations(grid, np.full(2 * grid.nq, 0.5))
    with pytest.raises(ValueError, match="six"):
        with_occupations(grid, np.full(2 * grid.nf, 0.5))


def test_ledgers_stay_separate_and_the_proxy_is_not_the_requirement():
    content = np.array([[6.096, 2.0, 1.0, 0.9], [6.213, 2.0, 1.0, 0.8]])
    flux = np.array([[-5.5, 3.0, 1.5, 1.0], [-5.534, 3.2, 1.4, 0.934]])
    windows = assess_windows(content, flux)
    ledgers = separate_ledgers(content, flux, {"field_energy_change": -0.396})
    assert windows["maintained"] is True
    assert windows["reversal"] == 0.0
    assert windows["candidate_regime"] is False
    assert windows["candidate_regime_is_programme_requirement"] is False
    assert windows["maintained_with_throughflow"] is True
    assert windows["stable_throughflow"] is False
    assert windows["splitting_one_drift_is_complete_renewal"] is False
    assert ledgers["content"]["maintained"] is True
    assert ledgers["flux"]["balanced"] is True
    assert ledgers["flux"]["throughflow_magnitude"] > 1.0
    assert ledgers["reversal"]["reversal"] == 0.0
    assert ledgers["accounting"]["field_energy_change"] == -0.396
    assert ledgers["proxies"]["programme_requirements"] is False
    uniform = uniform_complement_localization()
    assert uniform["localized_relative_to_empty_complement"] is True
    assert uniform["complement_occupation"] == 0.0


def test_v2_successor_replays_without_touching_v1():
    json_before = episode.sha256(regeneration.OUT)
    npz_before = episode.sha256(regeneration.NPZ)
    if not regeneration.OUT_V2.is_file():
        pytest.fail("v2 successor is missing")
    stored = json.loads(regeneration.OUT_V2.read_text())
    fresh = regeneration.assess_successor(write=False)
    assert episode.sha256(regeneration.OUT) == json_before
    assert episode.sha256(regeneration.NPZ) == npz_before
    assert stored["predecessor"]["json_sha256"] == json_before
    assert stored["predecessor"]["npz_sha256"] == npz_before
    assert fresh == stored
    preserved = stored["preserved"]
    assert preserved["attained_time"] == 0.1
    assert preserved["Q_min_end"] == 0.12477163551088903
    assert preserved["leader_content_change"] == 0.11695656049787839
    assert preserved["packet_flux_end"] == -5.534007341562886
    assert preserved["field_energy_change"] == -0.395776510604783
    assert preserved["G_min_end"] == 0.015567961028071387
    assert stored["dynamics_rerun"] is False
    assert stored["new_action"] is False
    assert stored["new_forces"] is False
    assert stored["stability"]["above_owned_cap"] is True
    assert stored["stability"]["below_absolute_rk4_limit"] is True
    assert stored["stability"]["chart_failure"] is None
    assert stored["stability"]["gram_admissible"] is True
    assert stored["stability"]["dt_times_omega"] > 1.4
    assert stored["galerkin_half_identity"]["half_identity_commutes_with_H_G"] is True
    assert stored["galerkin_half_identity"]["embedding_commutes_with_unprojected_H_fine"] is False
    assert stored["ledgers"]["content"]["maintained"] is True
    assert stored["ledgers"]["flux"]["balanced"] is True
    assert stored["ledgers"]["reversal"]["reversal"] == 0.0
    assert stored["ledgers"]["proxies"]["programme_requirements"] is False
    assert stored["ledgers"]["maintained_with_throughflow"] is True
    assert stored["uniform_six"]["localized_relative_to_empty_complement"] is True
    assert stored["uniform_six"]["window_proxy_localized"] is False
    assert stored["next_event"]["executed"] is False
    assert stored["next_event"]["estimated_cpu_seconds"] < 1.0
    assert stored["policy"]["g_is_dynamic_chart_stop"] is False
    assert stored["policy"]["splitting_one_drift_is_complete_renewal"] is False
