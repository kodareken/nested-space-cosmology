"""Matched proper-time prediction. Finite differences are the control, not the tangent."""
import os

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

import json
from pathlib import Path

import numpy as np
import pytest

from derive_nsc_discovery_prediction import main as cli_main
from recursive_horizons import nsc_discovery_episode as episode
from recursive_horizons import nsc_discovery_prediction as prediction
from recursive_horizons import nsc_discovery_response as discovery
from recursive_horizons import nsc_nested_parent_child as model
from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


_LAB = Path(__file__).resolve().parents[1]
_REPO = _LAB.parent


def _gap(actual, reference):
    difference = np.asarray(actual) - np.asarray(reference)
    return float(np.max(np.abs(difference)))


def _assert_close(actual, reference, name, rtol=2e-5, atol=2e-8):
    gap = _gap(actual, reference)
    scale = max(1e-8, float(np.max(np.abs(np.asarray(actual)))), float(np.max(np.abs(np.asarray(reference)))))
    assert gap <= rtol * scale + atol, (name, gap, scale, gap / scale)


def _manufactured_pair(nf=48, quadrature=96):
    phi0, phi1 = galerkin.manufactured_columns(nf, seed=4)
    phi1 = np.exp(-0.5j) * phi1
    pair = model.build_pair(
        nf, quadrature=quadrature, source_layout="override",
        columns_override=(phi0, phi1), occupations=coupling.OCCUPATIONS.copy(),
    )
    nodal = galerkin.blank_state(pair.grid, phi0, phi1)
    angle = 2 * np.pi * pair.grid.xi_g / pair.grid.length
    nodal.Q = nodal.Q * (1.0 + 0.03 * np.sin(angle))
    nodal.r = 1.0 + 0.04 * np.cos(angle)
    nodal.chi = 0.02 * np.sin(2.0 * angle)
    nodal.p_Q = 0.01 * np.cos(angle)
    nodal.p_r = 0.015 * np.sin(angle)
    nodal.p_chi = 0.012 * np.cos(2.0 * angle)
    nodal.phi1 = np.exp(-0.5j) * nodal.phi1
    return pair, model.encode_state(pair, nodal)


def _free_tangent(pair, seed=7):
    tangent = discovery.zero_tangent(pair.grid)
    angle = 2 * np.pi * pair.grid.xi_g / pair.grid.length
    generator = np.random.default_rng(seed)
    tangent.Q = 0.02 * np.sin(angle)
    tangent.r = 0.03 * np.cos(angle)
    tangent.chi = 0.01 * np.sin(2.0 * angle)
    tangent.p_r = 0.01 * np.sin(angle)
    raw0 = generator.normal(size=tangent.phi0.shape) + 1j * generator.normal(size=tangent.phi0.shape)
    raw1 = generator.normal(size=tangent.phi1.shape) + 1j * generator.normal(size=tangent.phi1.shape)
    tangent.phi0 = 1e-3 * raw0 / np.linalg.norm(raw0)
    tangent.phi1 = 1e-3 * raw1 / np.linalg.norm(raw1)
    tangent.occupations[:] = prediction.OCCUPATION_CONTRAST
    return tangent


def _displace(state, tangent, step):
    return model.NestedState(*(
        getattr(state, name) + step * getattr(tangent, name) for name in model.STATE_NAMES
    ))


def _solved_separated(nf=64, quadrature=128):
    pair = model.build_pair(nf, quadrature=quadrature)
    state, report = model.initial_state(pair)
    assert report["converged"] is True
    return pair, state


def _input_hashes():
    return prediction.saved_input_hashes()


def test_contract_keeps_stress_open_and_locks_the_contrast():
    contract = prediction.specification()
    assert contract["schema"] == "NSC-DISCOVERY-PREDICTION-v1"
    assert contract["method"] == "rk4_analytic_jacobian_vector"
    assert contract["clock_quadrature"] == "rk4_stages"
    assert contract["formula"] == "deltaO_t - Odot * delta_tau / tau_dot"
    assert contract["primary_readout"] == "child_regional_content"
    assert contract["secondary_readout"] == "child_proper_mean_r"
    assert contract["occupation_contrast"] == [1.0, 1.0, 0.0, 0.0, -1.0, -1.0]
    assert contract["contrast_trace"] == 0.0
    assert contract["derivative_h"] == [0.001, 0.0005]
    assert contract["held_out_alpha"] == 0.01
    assert contract["held_out_distinct_from_derivative_h"] is True
    assert 0.01 not in contract["derivative_h"]
    assert contract["finite_difference_used_as_analytic"] is False
    assert contract["newton_on_analytic_tangent"] is False
    assert contract["t03_source_changed"] is False
    assert contract["t03_constraint_resolved"] is False
    assert contract["effective_stress_evaluated"] is False
    assert contract["effective_stress"] is None
    assert "ctp_future_effective_stress" in contract["missing"]
    source = Path(prediction.__file__).read_text()
    assert "future_effective_stress_specification" not in source
    assert "solve_initial_radius(" not in source


def test_coupled_jv_step_matches_centred_state_map_on_manufactured_state():
    pair, state = _manufactured_pair()
    tangent = _free_tangent(pair)
    limit, _omega, _quadrature = model.stable_timestep(pair, state, 0.01)
    dt = min(0.01, float(limit))
    stepped = prediction.coupled_rk4_step(pair, state, tangent, dt)
    nonlinear = prediction.nonlinear_rk4_step(pair, state, dt)
    reference = model.rk4_step(pair, state, dt)
    for name in model.STATE_NAMES:
        assert _gap(getattr(stepped["state"], name), getattr(reference, name)) == 0.0
        assert _gap(getattr(nonlinear["state"], name), getattr(reference, name)) == 0.0
    assert np.array_equal(stepped["tangent"].occupations, prediction.OCCUPATION_CONTRAST)
    step = 1e-6
    base = np.array(pair.grid.fine.occupations, dtype=float, copy=True)

    def advance(sign):
        weights = base + sign * step * tangent.occupations
        trial = prediction.pair_with_occupations(pair, weights)
        return model.rk4_step(trial, _displace(state, tangent, sign * step), dt)

    plus, minus = advance(1.0), advance(-1.0)
    assert np.array_equal(pair.grid.fine.occupations, base)
    for name in model.STATE_NAMES:
        numeric = (getattr(plus, name) - getattr(minus, name)) / (2.0 * step)
        _assert_close(getattr(stepped["tangent"], name), numeric, "coupled " + name)
    eigenvalues = prediction.assert_admissible_source(stepped["state"].phi0, stepped["state"].phi1, base)
    assert float(np.min(eigenvalues)) >= -1e-8
    assert float(np.max(eigenvalues)) <= 1.0 + 1e-8


def test_proper_clock_integrates_rk_stages_and_tracks_a_finite_variation():
    pair, state = _manufactured_pair()
    tangent = _free_tangent(pair, seed=11)
    limit, _omega, _quadrature = model.stable_timestep(pair, state, 0.01)
    dt = min(0.01, float(limit))
    stepped = prediction.coupled_rk4_step(pair, state, tangent, dt)
    rate1 = model.rates(pair, state)
    jvp1 = discovery.nested_rate_jacobian_vector(pair, state, tangent)
    state2 = model._combine(state, rate1, dt / 2.0)
    tangent2 = prediction._shift_tangent(tangent, jvp1, dt / 2.0)
    rate2 = model.rates(pair, state2)
    jvp2 = discovery.nested_rate_jacobian_vector(pair, state2, tangent2)
    state3 = model._combine(state, rate2, dt / 2.0)
    tangent3 = prediction._shift_tangent(tangent, jvp2, dt / 2.0)
    rate3 = model.rates(pair, state3)
    jvp3 = discovery.nested_rate_jacobian_vector(pair, state3, tangent3)
    state4 = model._combine(state, rate3, dt)
    tangent4 = prediction._shift_tangent(tangent, jvp3, dt)
    tau_samples = []
    delta_samples = []
    for stage_state, stage_tangent in (
        (state, tangent), (state2, tangent2), (state3, tangent3), (state4, tangent4),
    ):
        tau_dot, delta_tau_dot = prediction.clock_rates(pair, stage_state, stage_tangent)
        tau_samples.append(tau_dot)
        delta_samples.append(delta_tau_dot)
    integrated_tau = dt * (tau_samples[0] + 2.0 * tau_samples[1] + 2.0 * tau_samples[2] + tau_samples[3]) / 6.0
    integrated_delta = dt * (
        delta_samples[0] + 2.0 * delta_samples[1] + 2.0 * delta_samples[2] + delta_samples[3]
    ) / 6.0
    assert stepped["tau"] == pytest.approx(integrated_tau)
    assert stepped["delta_tau"] == pytest.approx(integrated_delta)
    frozen_delta = delta_samples[0] * dt
    assert abs(integrated_delta - frozen_delta) > 1e-8
    amplitude = 1e-5
    base = np.array(pair.grid.fine.occupations, dtype=float, copy=True)

    def elapsed(sign):
        weights = base + sign * amplitude * tangent.occupations
        trial = prediction.pair_with_occupations(pair, weights)
        moved = _displace(state, tangent, sign * amplitude)
        return prediction.nonlinear_rk4_step(trial, moved, dt)["tau"]

    numeric_delta = (elapsed(1.0) - elapsed(-1.0)) / (2.0 * amplitude)
    _assert_close(stepped["delta_tau"], numeric_delta, "clock variation", rtol=2e-4, atol=2e-8)
    assert abs(numeric_delta - frozen_delta) > abs(numeric_delta - integrated_delta)


def test_centred_two_h_and_held_out_alpha_on_a_small_solved_source():
    pair, state = _solved_separated()
    stored = np.array(pair.grid.fine.occupations, dtype=float, copy=True)
    report = prediction.finite_difference_control(
        pair, state, duration=0.002, step_cap=0.001, fixed_step=0.001,
    )
    assert np.array_equal(pair.grid.fine.occupations, stored)
    assert report["finite_difference_used_as_analytic"] is False
    assert report["effective_stress"] is None
    assert report["held_out_alpha"] == 0.01
    assert report["derivative_h"] == [0.001, 0.0005]
    assert report["analytic"]["coordinate_time"] == pytest.approx(0.002)
    assert report["analytic"]["step_count"] == 2
    assert report["analytic"]["initial_state_called"] is False
    assert report["analytic"]["newton_used"] is False
    assert report["analytic"]["nonlinear_state_gap"] < 1e-12
    assert report["analytic"]["eigenvalues_admissible"] is True
    rows = {row["h"]: row for row in report["centred"]}
    assert set(rows) == {0.001, 0.0005}
    analytic_primary = report["analytic"]["primary_tau"]
    analytic_secondary = report["analytic"]["secondary_tau"]
    assert abs(analytic_primary) > 1e-8
    for step, row in rows.items():
        assert row["used_as_analytic"] is False
        _assert_close(row["primary"], analytic_primary, f"primary h={step}", rtol=2e-3, atol=2e-6)
        _assert_close(row["secondary"], analytic_secondary, f"secondary h={step}", rtol=5e-3, atol=2e-6)
    assert abs(rows[0.0005]["primary"] - analytic_primary) <= abs(rows[0.001]["primary"] - analytic_primary) + 1e-6
    held = report["held_out"]
    assert held["alpha"] == 0.01
    assert held["used_as_derivative_h"] is False
    assert held["eigenvalues_admissible"] is True
    assert np.sign(held["primary_scaled"]) == np.sign(analytic_primary)
    _assert_close(held["primary_scaled"], analytic_primary, "held-out primary", rtol=5e-2, atol=2e-5)
    held_gap = abs(held["primary_scaled"] - analytic_primary)
    small_gap = abs(rows[0.0005]["primary"] - analytic_primary)
    assert held_gap > small_gap
    for item in report["preparations"]:
        assert item["constructor"] == "nsc_nested_parent_child.initial_state"
        assert item["constructor_grid"] == "dense"
        assert item["converged"] is True
        assert item["constraint_resolved"] is True
        assert item["eigenvalues_admissible"] is True
        assert item["final_eigenvalues_admissible"] is True
        assert abs(sum(item["weights"]) - float(np.sum(stored))) < 1e-12
        assert min(item["weights"]) > 0.0
        assert max(item["weights"]) <= 1.0
    prepared_alphas = sorted(round(item["alpha"], 6) for item in report["preparations"])
    assert prepared_alphas == [-0.001, -0.0005, 0.0005, 0.001, 0.01]


def test_saved_baselines_match_a_short_fft_prefix_and_do_not_resolve_t03(monkeypatch):
    before = _input_hashes()

    def forbid(*_args, **_kwargs):
        raise AssertionError("saved prediction called the initial constructor")

    monkeypatch.setattr(model, "initial_state", forbid)
    monkeypatch.setattr(model, "solve_initial_radius_successor", forbid)
    for case_name in episode.BASELINE_CASES:
        seed = prediction.load_saved_baseline(case_name, 0.0)
        final = prediction.load_saved_baseline(case_name, 0.3)
        assert seed["initial_state_called"] is False
        assert seed["constraint_resolved"] is False
        assert seed["source_changed"] is False
        assert final["comparison_frame"] is True
        assert final["constraint_resolved"] is False
        assert final["source_changed"] is False
        assert seed["pins"]["weights_sha256"] == final["pins"]["weights_sha256"]
        assert seed["state_sha256"] == seed["initial_state_sha256"]
        assert final["state_sha256"] == final["final_state_sha256"]
        assert seed["state_sha256"] != final["state_sha256"]
        identity = prediction.compare_saved_baseline(final["state"], case_name, 0.3)
        assert identity["max_gap"] == 0.0
        assert identity["initial_state_called"] is False
        assert identity["constraint_resolved"] is False
        assert identity["source_changed"] is False
        assert identity["weights"] == [0.75, 0.75, 0.5, 0.5, 0.25, 0.25]
        evolved = prediction.propagate_baseline(
            seed["pair"], seed["state"], duration=0.01, step_cap=seed["step_cap"],
        )
        assert evolved["coordinate_time"] == pytest.approx(0.01)
        assert evolved["coordinate_time"] <= prediction.NONPRODUCTION_DURATION
        assert evolved["initial_state_called"] is False
        assert evolved["evolved_clock_integral"] is True
        assert evolved["frozen_rate_interval"] is False
        assert evolved["eigenvalues_admissible"] is True
        compared = prediction.compare_saved_baseline(evolved["state"], case_name, 0.01)
        assert compared["max_gap"] < 1e-8, (case_name, compared["gaps"])
        assert compared["source_changed"] is False
        assert _gap(evolved["state"].r, seed["state"].r) > 1e-8
    for case_name in ("nf128_baseline_dt0.0005", "nf256_baseline_dt0.0005"):
        seed = prediction.load_saved_baseline(case_name, 0.0)
        forecast = prediction.predict_occupation_response(
            seed["pair"], seed["state"], duration=0.0005, step_cap=seed["step_cap"], fixed_step=0.0005,
        )
        assert forecast["initial_state_called"] is False
        assert forecast["newton_used"] is False
        assert forecast["t03_constraint_resolved"] is False
        assert forecast["t03_source_changed"] is False
        assert forecast["effective_stress"] is None
        assert forecast["evolved_clock_integral"] is True
        assert forecast["frozen_rate_interval"] is False
        assert forecast["eigenvalues_admissible"] is True
        assert forecast["propagated"]["coordinate_time"] == pytest.approx(0.0005)
        assert forecast["radius_tangent"] == "prepared_nested_radius_tangent"
        assert abs(forecast["primary_tau"]) > 0.0
        assert forecast["readout"]["clock_quadrature"] == "rk4_stages"
        assert np.isfinite(forecast["propagated"]["delta_tau"])
        assert np.array_equal(forecast["propagated"]["tangent"].occupations, prediction.OCCUPATION_CONTRAST)
    after = _input_hashes()
    assert before == after
    with pytest.raises(PermissionError, match="production"):
        prediction.propagate_baseline(
            seed["pair"], seed["state"], duration=0.02, step_cap=seed["step_cap"],
        )


def test_cost_counts_and_cli_create_only(tmp_path):
    one = prediction.expected_cost(0.3, 0.0005)
    assert one["steps_at_cap"] == 600
    assert one["analytic_rate_evaluations"] == 2400
    assert one["analytic_jacobian_vector_evaluations"] == 2400
    assert one["nonlinear_baseline_rate_evaluations"] == 2400
    assert one["finite_difference_trajectories"] == 4
    assert one["finite_difference_rate_evaluations"] == 9600
    assert one["finite_difference_dense_initial_solves"] == 4
    assert one["held_out_trajectories"] == 1
    assert one["held_out_rate_evaluations"] == 2400
    assert one["held_out_dense_initial_solves"] == 1
    assert one["effective_stress_evaluations"] == 0
    assert one["t03_resolves"] == 0
    campaign = prediction.production_cost()
    assert campaign["executed"] is False
    assert campaign["comparison_time"] == 0.3
    assert len(campaign["cases"]) == 4
    totals = campaign["totals_if_caps_are_prepared_separately"]
    assert totals["steps_at_cap"] == 1800
    assert totals["analytic_rate_evaluations"] == 7200
    assert totals["analytic_jacobian_vector_evaluations"] == 7200
    assert totals["nonlinear_baseline_rate_evaluations"] == 7200
    assert totals["finite_difference_rate_evaluations"] == 28800
    assert totals["held_out_rate_evaluations"] == 7200
    assert totals["finite_difference_dense_initial_solves"] == 16
    assert totals["finite_difference_dense_initial_solves_if_shared_across_caps"] == 8
    assert totals["held_out_dense_initial_solves_if_shared_across_caps"] == 2
    output = tmp_path / "nsc-discovery-prediction-spec"
    assert cli_main(["--create", "--output", str(output)]) == 0
    manifest = json.loads((output / "prediction-manifest.json").read_text())
    assert manifest["creation_only"] is True
    assert manifest["evolved"] is False
    assert manifest["production_trajectory"] is False
    assert manifest["later_outputs_written"] is False
    assert manifest["settings"]["held_out_alpha"] == 0.01
    assert manifest["settings"]["derivative_h"] == [0.001, 0.0005]
    assert manifest["settings"]["t03_constraint_resolved"] is False
    assert manifest["expected_cost"]["executed"] is False
    assert {frame["case_name"] for frame in manifest["frames"]} == set(episode.BASELINE_CASES)
    for name in prediction.LATER_OUTPUTS:
        assert not (output / name).exists()
    assert cli_main(["--check", "--output", str(output)]) == 0
    refused = _REPO / "lab" / "tests" / "nsc-discovery-prediction-not-a-result"
    with pytest.raises(PermissionError):
        prediction.create_successor(refused)
    assert not refused.exists()
    forecast = prediction.cpu_forecast(256, 0.3, 0.0005, n_cases=2)
    assert forecast["forecast_cpu_seconds"] > 0.0
    assert forecast["hard_gate"] is False
    assert forecast["process_pool"] is False
    assert forecast["single_process"] is True
    assert forecast["rate_evaluations"] > 7200
    short = prediction.accepted_duration(0.002, production=False)
    assert short["target"] == "short_test" and short["production"] is False
    assert prediction.accepted_duration(0.3, production=True)["target"] == "T0.3"
    assert prediction.accepted_duration(1.0, production=True)["target"] == "optional_T1"
    with pytest.raises(PermissionError):
        prediction.accepted_duration(0.3, production=False)
    with pytest.raises(ValueError, match="T=3"):
        prediction.accepted_duration(3.0, production=True)


def test_equal_proper_time_root_is_separate_from_the_taylor_indicator():
    pair, state = _solved_separated()
    report = prediction.equal_proper_time_comparison(
        pair, state, duration=0.002, step_cap=0.001, fixed_step=0.001,
    )
    assert report["hard_acceptance_gate"] is False
    assert report["initial_state_called_on_baseline"] is False
    assert abs(report["tau_target"]) > 0.0
    for row in report["centred"]:
        assert row["reached"] is True
        assert row["method"] == "bracketed_decreasing_dt"
        assert row["uses_linear_clock_correction"] is False
        assert row["taylor_indicator"]["equal_proper_time"] is False
        assert row["taylor_indicator"]["same_formula_as_prediction"] is True
        assert np.isfinite(row["space_error_primary"])
        assert np.isfinite(row["space_error_secondary"])
        assert abs(row["time_error_plus"]["tau_residual"]) < 1e-8
        assert abs(row["time_error_minus"]["tau_residual"]) < 1e-8
        assert row["time_error_plus"]["bisections"] > 0
    held = report["held_out"]
    assert held["alpha"] == 0.01
    assert held["reached"] is True
    assert held["used_as_derivative_h"] is False
    assert held["uses_linear_clock_correction"] is False
    assert held["taylor_indicator"]["uses_linear_clock_correction"] is True
    assert held["taylor_indicator"]["equal_proper_time"] is False
    assert held["space_error"] == pytest.approx(held["measured_change"] - held["predicted_change"])
    assert np.isfinite(held["time_error"]["coordinate_offset"])
    arm = next(item for item in report["arms"] if abs(item["alpha"] - 0.01) <= 1e-15)
    observed = prediction._observables(arm["pair"], arm["event"]["state"])
    assert observed["child_regional_content"] == pytest.approx(held["child_regional_content"])
    assert arm["event"]["dense_interpolation"]["available"] is True
    assert np.isfinite(arm["event"]["dense_interpolation_content_gap"])
    assert np.all(np.diff(arm["history"]["tau"]) > 0.0)


def test_cli_run_check_replays_saved_states_without_evolving(tmp_path, monkeypatch):
    output = tmp_path / "nsc-discovery-prediction-run"
    assert cli_main([
        "--run", "--output", str(output), "--duration", "0.0005", "--nf", "128", "--step-cap", "0.0005",
    ]) == 0
    record = json.loads((output / "prediction-run.json").read_text())
    assert record["status"] == "RUN"
    assert record["evolved"] is True
    assert record["production_trajectory"] is False
    assert record["creation_only"] is True
    assert record["hard_acceptance_gate"] is False
    assert record["equal_proper_time_method"] == "bracketed_decreasing_dt"
    assert record["taylor_indicator_role"] == prediction.TAYLOR_ROLE
    assert record["cpu_forecast"]["forecast_cpu_seconds"] > 0.0
    assert record["cpu_forecast"]["process_pool"] is False
    assert record["late_T3_curvature_unresolved"] is True
    assert record["cases"][0]["nf"] == 128
    assert record["cases"][0]["t03_constraint_resolved"] is False
    assert record["cases"][0]["centred_equal_tau"][0]["uses_linear_clock_correction"] is False
    assert record["cases"][0]["centred_equal_tau"][0]["reached"] is True
    assert record["cases"][0]["centred_equal_tau"][0]["method"] == "bracketed_decreasing_dt"
    assert record["cases"][0]["centred_equal_tau"][1]["reached"] is True
    assert record["cases"][0]["held_out_equal_tau"]["alpha"] == 0.01
    assert record["cases"][0]["held_out_equal_tau"]["reached"] is True
    assert record["cases"][0]["held_out_equal_tau"]["uses_linear_clock_correction"] is False
    assert record["cases"][0]["held_out_equal_tau"]["taylor_indicator"]["equal_proper_time"] is False
    assert (output / "prediction-run.npz").is_file()
    with pytest.raises(FileExistsError):
        prediction.run_saved_comparison((128,), 0.0005, output, step_cap=0.0005, production=False)

    def forbid(*_args, **_kwargs):
        raise AssertionError("readonly replay called an evolution or an initializer")

    monkeypatch.setattr(model, "rk4_step", forbid)
    monkeypatch.setattr(model, "initial_state", forbid)
    monkeypatch.setattr(model, "solve_initial_radius_successor", forbid)
    monkeypatch.setattr(prediction, "propagate_baseline", forbid)
    monkeypatch.setattr(prediction, "propagate_coupled", forbid)
    monkeypatch.setattr(prediction, "equal_proper_time_comparison", forbid)
    assert cli_main(["--check", "--output", str(output)]) == 0
    replay = prediction.replay_run_record(output)
    assert replay["evolution_called"] is False
    assert replay["readonly"] is True
    assert replay["producer_match"] is True
    assert replay["max_observable_gap"] < 1e-8
    assert replay["max_tau_residual"] < 1e-8
