"""Replay and accounting tests for the spherical feedback episode.

The long T=0.05 evidence is read back. These tests do not start another
production episode and do not retune its floors.
"""
import json
import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import numpy as np

import derive_nsc_spherical_feedback_episode as episode
from recursive_horizons import nsc_regional_energy_exchange as regional
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin

PAYLOAD_LIMIT = 64 * 1024 * 1024


def _nf512():
    grid, state, metadata = episode.load_v5_state("nf512")
    return grid, state.copy(), metadata


def test_saved_v5_state_is_the_episode_setup_without_a_new_solve():
    grid, state, metadata = _nf512()
    owned_grid, owned_state, radius, digest = regional.load_nf512_initial()
    assert metadata["v5_sha256"] == digest
    assert metadata["occupation_gap"] == 0.0
    assert metadata["quadrature_radius_gap"] == 0.0
    assert np.array_equal(state.r, owned_state.r)
    assert np.array_equal(state.Q, owned_state.Q)
    assert np.array_equal(state.chi, owned_state.chi)
    assert np.array_equal(state.p_Q, owned_state.p_Q)
    assert np.array_equal(state.p_r, owned_state.p_r)
    assert np.array_equal(state.p_chi, owned_state.p_chi)
    assert np.array_equal(state.phi0, owned_state.phi0)
    assert np.array_equal(state.phi1, owned_state.phi1)
    assert np.allclose(galerkin.prolong_geometry(grid, state.r), radius)
    assert grid.nf == owned_grid.nf == 512
    assert grid.nq == owned_grid.nq == 2048
    coarse, _fine = galerkin.load_physical_columns(512)[:2]
    assert np.array_equal(state.phi0, coarse)
    assert np.array_equal(state.phi1, _fine)
    called = {"solve": 0}
    original = galerkin.solve_initial_radius

    def forbidden(*_args, **_kwargs):
        called["solve"] += 1
        raise AssertionError("the episode re-solved the initial radius")

    galerkin.solve_initial_radius = forbidden
    try:
        sample, _rate, fields = episode.observe(grid, state)
    finally:
        galerkin.solve_initial_radius = original
    assert called["solve"] == 0
    assert np.max(np.abs(fields["rho"] - metadata["saved_rho"])) < 1e-10
    assert abs(sample["full_hamilton_max"] - 9.981468739539423e-06) < 1e-12
    assert abs(sample["projected_hamilton_max"] - 8.52043675013457e-09) < 1e-15
    assert sample["chart_proper_max"] == 0.0
    assert abs(sample["lifted_proper_max"] - 1.6655044770048675e-10) < 1e-16
    decision = episode.initial_tolerance_decision(sample["full_hamilton_max"], sample["full_momentum_max"])
    assert decision["historical_initial_1e-8_passed"] is False
    assert decision["abort_evolution"] is False
    assert decision["continue_evolution"] is True
    assert sample["mean_subtracted"] is False
    assert sample["rho_mean_removed_before_constraint"] is False
    assert sample["frozen_source_used"] is False


def test_proper_motion_is_not_coordinate_velocity_or_coordinate_work():
    grid, state, _metadata = _nf512()
    shift = np.array(grid.fine.shift, copy=True)
    length = np.array(grid.fine.length_density, copy=True)
    sample, _rate, fields = episode.observe(grid, state)
    assert sample["proper_shift_identity_gap"] < 1e-12
    assert abs(sample["signed_versus_owner_gap"]) < 1e-12
    assert sample["c2_identity_gap"] < 1e-12
    assert sample["coordinate_over_N_max"] > 0.05
    assert sample["lifted_proper_max"] < 1e-8
    assert sample["coordinate_work_l1"] > 1.0
    assert sample["proper_pressure_work_l1"] < 1e-6
    assert sample["coordinate_work_l1"] > 100.0 * max(sample["proper_pressure_work_l1"], 1e-12)
    assert abs(sum(sample["window_coordinate_work"]) - sample["lifted_fieldwork"]) < 1e-8
    assert abs(sum(sample["window_proper_work"]) - sample["proper_pressure_power"]) < 1e-8
    assert abs(sum(sample["window_proper_flux"])) < 1e-8
    regional_record = json.loads(episode.REGIONAL_JSON.read_text())
    saved_rows = regional_record["nf512"]["initial"]["windows"]
    for saved, row in zip(saved_rows, sample["windows"], strict=True):
        assert saved["name"] == row["name"]
        assert abs(saved["normal_energy"] - row["normal_energy"]) < 1e-9
        assert abs(saved["matter_channels"]["coordinate_metric_work"] - row["coordinate_metric_work"]) < 1e-9
        assert abs(saved["proper_pressure_work"] - row["proper_pressure_work"]) < 1e-9
        assert abs(saved["matter_channels"]["proper_normal_flux"] - row["proper_normal_flux"]) < 1e-9
        occupations = [column["occupation"] for column in row["mode_columns"]]
        assert occupations == [float(value) for value in grid.fine.occupations]
        assert sorted(occupations) == [0.25, 0.25, 0.5, 0.5, 0.75, 0.75]
        column_energy = sum(column["hamiltonian_energy"] for column in row["mode_columns"])
        assert abs(column_energy - row["matter_energy"]) < 1e-8
    projectors = sample["mode_projectors"]
    assert projectors["held_out_complement_hamilton_max"] > 0.99 * projectors["full_quadrature_hamilton_max"]
    assert projectors["geometry_band_projected_hamilton_max"] < 0.01 * projectors["full_quadrature_hamilton_max"]
    assert set(projectors) >= {
        "full_quadrature_hamilton_max",
        "geometry_band_projected_hamilton_max",
        "held_out_complement_hamilton_max",
        "lifted_proper_velocity",
        "unprojected_proper_velocity",
        "chart_proper_velocity",
    }
    mapped = episode._series_map([sample, sample])
    assert mapped["window_normal"].shape == (2, 4)
    assert mapped["window_proper_flux"].shape == (2, 4)
    assert mapped["proper_max"].shape == (2,)
    assert abs(sample["mean_rho"] - float(np.mean(fields["rho"]))) < 1e-12
    assert abs(sample["mean_current"] - float(np.mean(fields["current"]))) < 1e-12
    stepped = galerkin.rk4_step(grid, state, episode.DT)
    after, _rate_after, _fields_after = episode.observe(grid, stepped)
    assert after["proper_shift_identity_gap"] < 1e-12
    assert after["lifted_proper_max"] > 10.0 * sample["lifted_proper_max"]
    assert abs(after["proper_pressure_power"] - after["lifted_fieldwork"]) > 1e-3
    assert np.array_equal(grid.fine.shift, shift)
    assert np.array_equal(grid.fine.length_density, length)


def test_rk4_recomputes_the_source_on_four_distinct_stages_and_replays():
    grid, state, _metadata = _nf512()
    calls = []
    original = galerkin.rates

    def wrapped(grid_argument, state_argument, include_matter_force=True):
        assert include_matter_force is True
        calls.append(float(np.max(np.abs(state_argument.Q - state.Q))))
        return original(grid_argument, state_argument, include_matter_force)

    galerkin.rates = wrapped
    try:
        first = galerkin.rk4_step(grid, state, episode.DT)
    finally:
        galerkin.rates = original
    assert len(calls) == 4
    assert calls[0] == 0.0
    assert max(calls) > 0.0
    assert len(set(np.round(calls, decimals=12))) == 4
    second = galerkin.rk4_step(grid, state.copy(), episode.DT)
    for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1"):
        assert np.array_equal(getattr(first, name), getattr(second, name))
    assert not np.array_equal(first.Q, state.Q)


def test_resolution_rule_preserves_unresolved_and_incomplete_runs():
    assert episode.budget_allows(10.0, 0.5, 300.0) is True
    assert episode.budget_allows(299.5, 0.5, 300.0) is False
    assert episode.budget_allows(300.0, 0.0, 300.0) is False
    assert episode.classify_movement(1.0, 0.005, 1e-8) == "resolved"
    assert episode.classify_movement(1.0, 0.02, 1e-8) == "unresolved"
    assert episode.classify_movement(1e-12, 1e-12, 1e-8) == "below_noise_floor"
    assert episode.classify_movement(1e-12, 0.1, 1e-8) == "unresolved"
    assert episode.summarize_statuses(["resolved", "below_noise_floor"]) == "resolved"
    assert episode.summarize_statuses(["resolved", "unresolved"]) == "unresolved"
    assert episode.summarize_statuses(["below_noise_floor", "below_noise_floor"]) == "no_measurable_evolution"
    assert episode.summarize_statuses(["resolved", "incomplete"]) == "incomplete"
    primary = {"time": np.array([0.0, 0.05]), "proper_max": np.array([0.0, 1.0])}
    agree = {"time": np.array([0.0, 0.05]), "proper_max": np.array([0.0, 1.005])}
    disagree = {"time": np.array([0.0, 0.05]), "proper_max": np.array([0.0, 1.02])}
    short = {"time": np.array([0.0, 0.01]), "proper_max": np.array([0.0, 0.2])}
    resolved = episode.effect_row("proper_velocity_max", "proper_max", 1e-8, primary, agree, 0.05, 0.05)
    unresolved = episode.effect_row("proper_velocity_max", "proper_max", 1e-8, primary, disagree, 0.05, 0.05)
    incomplete = episode.effect_row("proper_velocity_max", "proper_max", 1e-8, primary, short, 0.0005, 0.05)
    assert resolved["status"] == "resolved"
    assert abs(resolved["effect_scale"] - 1.0) < 1e-12
    assert unresolved["status"] == "unresolved"
    assert abs(unresolved["movement"] - 0.02) < 1e-12
    assert incomplete["status"] == "incomplete"
    measured = episode.verdict_from("resolved", "unresolved", True, True, False, False, False)
    assert measured["verdict"] == "MEASURED_FEEDBACK_UNRESOLVED_CONSTRAINT_CONTROL"
    assert measured["effect_resolved"] is True
    assert measured["constraint_consistent_solution"] is False
    assert measured["sensitivity"] is None
    per_row = episode.verdict_from("unresolved", "resolved", True, True, False, False, False)
    assert per_row["verdict"] == "MEASURED_FEEDBACK_PER_OBSERVABLE"
    assert per_row["verdict"] != "EPISODE_MEASURED_EFFECT_UNRESOLVED"
    assert per_row["sensitivity"] is None
    partial = episode.verdict_from("resolved", "resolved", False, True, False, False, True)
    assert partial["verdict"] == "PARTIAL_CPU_BUDGET"
    assert partial["effect_resolved"] is False
    witness = episode.verdict_from("resolved", "resolved", True, False, False, False, False)
    assert witness["verdict"] == "SAME_MODEL_WITNESS_FAILED"
    balance = episode.balance_against_exchange(3.442300311462532e-08, -0.41230177531171325)
    assert balance["within_one_percent_of_exchange"] is True
    assert balance["balance_over_exchange"] < 1e-6
    assert balance["not_one_percent_of_the_balance_error"] is True


def test_saved_episode_replays_and_keeps_the_declared_rule():
    record = episode.verify_saved(replay=True)
    assert record["renewal"] is False
    assert record["effect_resolved"] in (True, False)
    assert record["incoming_gate_prerequisite"] is False
    assert record["inheritance_T_star_recomputed"] is False
    assert record["lambda_cdm_conflict_claimed"] is False
    assert record["model"]["state_reset"] is False
    assert record["model"]["restoring_force"] is False
    assert record["model"]["constraint_projection"] is False
    assert record["model"]["hamiltonian_cached_across_stages"] is False
    assert record["payload_within_64MiB"] is True
    assert record["payload_bytes"] <= PAYLOAD_LIMIT
    assert record["frozen_bytes_unchanged"] is True
    assert record["setup_nf512"]["frozen_rho_used_in_evolution"] is False
    assert record["setup_nf512"]["historical_initial_1e-8"]["abort_evolution"] is False
    names = [item["name"] for item in record["budget"]["intended_runs"]]
    assert names == [spec["name"] for spec in episode.INTENDED_RUNS]
    assert set(record["runs"]) <= set(names)
    v5 = json.loads(episode.V5_JSON.read_text())
    assert v5["time_extension_T_0_05"] is False
    assert episode.sha256(episode.V5_NPZ) == record["hashes_after"]["v5_npz"]
    assert episode.sha256(episode.V5_JSON) == record["hashes_after"]["v5_json"]
    primary = record["runs"]["nf512_dt_0_0005"]
    assert primary["mean_subtracted"] is False
    with np.load(episode.NPZ, allow_pickle=False) as payload:
        attained = float(payload["nf512_dt_0_0005_time"][-1])
        assert abs(attained - primary["attained_T"]) <= 1e-12
        proper = payload["nf512_dt_0_0005_proper_max"]
        assert abs(float(proper[0]) - primary["initial"]["proper_max"]) <= 1e-12
        assert abs(float(proper[-1]) - primary["final"]["proper_max"]) <= 1e-12
    witness = record["same_model_witness"]
    assert witness["matches_saved_v5_window"] is True
    assert witness["max_abs_gap"] <= 1e-8
    for item in record["physical_effects"]:
        for axis in ("time", "space"):
            row = item[axis]
            if row["status"] == "incomplete" or row.get("effect_scale") is None:
                continue
            assert row["status"] == episode.classify_movement(row["effect_scale"], row["movement"], row["floor"])
            assert abs(abs(row["primary_change"]) - row["effect_scale"]) <= 1e-12
    assert record["answer"]["historical_initial_1e-8_aborted_the_run"] is False
    assert record["answer"]["renewal"] is False
    assert record["verdict"] == "MEASURED_FEEDBACK_UNRESOLVED_CONSTRAINT_CONTROL"
    assert record["effect_resolved"] is True
    assert record["constraint_consistent_solution"] is False
    assert record["assessment_revision"]["numerical_values_recomputed"] is False
    assert record["assessment_revision"]["revised_without_rerunning_evolution"] is True
    assert record["constraint_control"]["consistent_solution_validated"] is False
    assert record["constraint_control"]["continuum_error"] == "unknown"
    indicators = record["refinement_indicators"]
    assert len(indicators) == 31
    assert all(item["meets_one_percent"] for item in indicators)
    reported = record["energy_balance"]["nf512_dt_0_0005"]
    assert reported["balance_versus_field_exchange"]["within_one_percent_of_exchange"] is True
    assert reported["balance_versus_coordinate_work"]["within_one_percent_of_exchange"] is True
    assert abs(reported["field_energy_change"]) > 0.4
    assert abs(reported["total_energy_change"]) < 1e-7
    assert reported["balance_versus_field_exchange"]["balance_over_exchange"] < 1e-6
