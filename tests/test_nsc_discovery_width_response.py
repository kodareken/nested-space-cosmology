"""Bounded width-frame, implicit preparation, Jv and sealed held-out checks."""
import os
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import json
import numpy as np
import pytest

from derive_nsc_discovery_width_response import main as cli_main
from recursive_horizons import nsc_discovery_width_response as width


@pytest.fixture(scope="module")
def baseline():
    return width.family.prepare_member(1., 1., nf=64, basis="owned_geometry_frame")


def test_continuous_adapter_extends_exact_owned_frames_without_rotating_weights(baseline):
    pair = baseline["pair"]
    for w in (0.8, 1., 1.2):
        phi0, phi1, meta = width.source_columns(pair.grid, w)
        owned0, owned1, _owned_meta = width.family.owned_separated_columns(pair.grid, w)
        assert np.array_equal(phi0, owned0)
        assert np.array_equal(phi1, owned1)
        assert meta["gram_max"] < 1e-9
        assert not meta["unequal_weight_alignment_rotations"]
        assert meta["occupation_trace"] == 3
        assert meta["CAR_min"] >= -1e-8 and meta["CAR_max"] <= 1+1e-8
    fingerprints = []
    for w in (1.001, 0.999, 1.0005, 0.9995, 1.05):
        changed = width.pair_at_width(pair, w)
        fingerprints.append(changed.source_metadata["frame_sha256"])
        assert changed.geometry_map is pair.geometry_map
        assert changed.reference_columns is pair.reference_columns
        assert np.array_equal(changed.weights, pair.weights)
        assert changed.source_metadata["width_hex"] == float(w).hex()
        assert changed.source_metadata["child_support"] == [2-w, 2+w]
        assert changed.source_metadata["child_frame_origin"] == "declared_continuous_width_lobes"
        assert not changed.source_metadata["child_source_equals_reference"]
        assert changed.source_metadata["fine_outside_support_power_fraction"] == changed.source_metadata["tails"]
    assert len(set(fingerprints)) == 5
    assert np.array_equal(pair.source_columns, np.vstack((pair.source_phi0, pair.source_phi1)))


def test_column_direction_retains_implicit_radius_and_cancels_constraint(baseline):
    pair, state = baseline["pair"], baseline["state"]
    tangent, report, _frames = width.frame_direction(pair, state, 0.001)
    assert np.linalg.norm(tangent.phi0) > 0
    assert np.linalg.norm(tangent.r) > 0
    assert report["physical_covariance_change_measured"]
    assert report["implicit_radius_linear_residual"] < 1e-8
    assert report["tangency_max"] < 0.01
    assert not report["newton_used_for_tangent"]
    assert np.array_equal(tangent.occupations, np.zeros(6))
    for name in ("Q", "chi", "p_Q", "p_r", "p_chi"):
        assert np.array_equal(getattr(tangent, name), np.zeros(pair.grid.ng))
    nodal = width.model.reconstruct_state(pair, state)
    delta = width.prediction._nodal_tangent(pair, tangent)
    def residual(sign, change_radius):
        trial = nodal.copy()
        trial.phi0 += sign*1e-6*delta.phi0
        trial.phi1 += sign*1e-6*delta.phi1
        if change_radius:
            trial.r += sign*1e-6*delta.r
        fine = width.galerkin.prolong_state(pair.grid, trial)
        system = width.galerkin.active_fine_system(pair.grid, fine)
        source = width.coupling.source_from_columns(system, fine)
        return width.galerkin.projected_radius_operator(pair.grid, trial.r, source["force_L"]/pair.grid.dx_q)[0]
    source_only = (residual(1, False)-residual(-1, False))/2e-6
    joint = (residual(1, True)-residual(-1, True))/2e-6
    assert np.max(np.abs(joint)) < 2e-5*max(1., np.max(np.abs(source_only)))


def test_source_support_length_derivative_includes_endpoint_term(baseline):
    pair, state = baseline["pair"], baseline["state"]
    tangent, _report, _frames = width.frame_direction(pair, state, 0.0005)
    analytic = width.proper_geometry(pair, state, 1., tangent)["width_derivative_at_fixed_t"]
    epsilon = 1e-6
    displaced = []
    for sign in (-1, 1):
        trial = width.model.NestedState(*(getattr(state, name)+sign*epsilon*getattr(tangent, name)
                                        for name in width.FIELDS))
        displaced.append(width.proper_geometry(pair, trial, 1.+sign*epsilon))
    derivative = (displaced[1]["source_support_proper_length"]-displaced[0]["source_support_proper_length"])/(2*epsilon)
    assert analytic["source_support_endpoint_term"] > 0
    assert derivative == pytest.approx(analytic["source_support_length"], rel=2e-6, abs=1e-6)
    temporal = width.proper_geometry(pair, state, 1., tangent, width_derivative=0.)
    assert temporal["width_derivative_at_fixed_t"]["source_support_endpoint_term"] == 0


def test_predeclared_full_precision_contract_and_budget_limits():
    config = width.inputs()
    assert config["held_out_width"] == 1.05
    assert config["derivative_h"] == [0.001, 0.0005]
    assert config["cpu_budget_seconds"] == 300
    assert config["forecast_factor"] == 1.5
    assert config["universal_scale_power"] is None
    assert config["Q0_varied"] is False
    assert config["inputs_sha256"] != width.inputs(step_cap=0.000500001)["inputs_sha256"]
    with pytest.raises(PermissionError):
        width.inputs(nf=256, duration=0.3)
    with pytest.raises(ValueError):
        width.inputs(cpu_budget=301)
    with pytest.raises(ValueError):
        width.inputs(nf=64)


def test_sealed_forecast_precedes_independent_equal_tau_arm(tmp_path, monkeypatch):
    directory = tmp_path/"width-preview"
    prepared = width.prepare(directory, nf=128, duration=0.003, step_cap=0.003)
    assert prepared["held_out_source_not_constructed"]
    assert not (directory/"measurement.json").exists()
    assert (directory/"prepare.json").stat().st_mode & 0o222 == 0
    assert prepared["directions"][0]["covariance_derivative_frobenius"] > 0
    assert not prepared["directions"][0]["frame_difference_is_analytic_width_derivative"]
    forecast = width.predict(directory)
    assert forecast["status"] == "predicted"
    assert forecast["coordinate_time"] == pytest.approx(0.003)
    assert forecast["initial_delta_phi_retained"]
    assert len(forecast["coefficients"]) == 2
    assert abs(forecast["coefficients"][-1]["coefficient_at_equal_tau"]) > 0
    sealed_digest = width.episode.file_sha256(directory/"prediction.json")
    actual = width.family.solve_prepared
    seen = []
    def independent_constructor(pair):
        assert width.episode.file_sha256(directory/"prediction.json") == sealed_digest
        assert pair.source_metadata["width_scale"] == 1.05
        seen.append(pair.source_metadata["frame_sha256"])
        return actual(pair)
    monkeypatch.setattr(width.family, "solve_prepared", independent_constructor)
    measurement = width.measure(directory)
    assert len(seen) == 1
    assert width.episode.file_sha256(directory/"prediction.json") == sealed_digest
    assert measurement["status"] == "measured_at_equal_tau"
    assert measurement["event"]["uses_linear_clock_correction"] is False
    assert measurement["W_unchanged"] and measurement["observer_unchanged"]
    assert not measurement["initial_radius_copied"]
    assert measurement["universal_percent_gate"] is None
    assert measurement["signed_nonlinear_remainder"] is not None
    assert measurement["aggregate_cpu_seconds"] <= 300
    assert measurement["coordinate_time"] <= 0.01
    assert measurement["gaussian_state"]["actual_covariance_trace"] == pytest.approx(3., abs=1e-5)
    _forecast_record, forecast_arrays = width.load(directory, "prediction")
    _measurement_record, held_arrays = width.load(directory, "measurement")
    assert np.max(np.abs(held_arrays["initial_r"]-forecast_arrays["initial_r"])) > 1e-6
    assert not np.array_equal(held_arrays["source_phi0"], forecast_arrays["source_phi0"])
    assert width.check(directory)["ok"]
    with pytest.raises(FileExistsError):
        width.predict(directory)


def test_budget_stop_seals_last_admissible_arrays_without_launching_heldout(tmp_path, monkeypatch):
    directory = tmp_path/"width-stop"
    width.prepare(directory, nf=128, duration=0.002, step_cap=0.002)
    monkeypatch.setattr(width, "_budget", lambda *args, **kwargs: (True, 300., 301.))
    stopped = width.predict(directory)
    assert stopped["status"] == "budget_stop"
    assert stopped["steps"] == 0
    _record, arrays = width.load(directory, "prediction")
    assert np.array_equal(arrays["baseline_r"], arrays["initial_r"])
    with pytest.raises(ValueError, match="did not reach"):
        width.measure(directory)
    assert not (directory/"measurement.json").exists()
    assert width.check(directory)["ok"]


def test_failed_heldout_constructor_is_retained_without_nonexistence_claim(tmp_path, monkeypatch):
    directory = tmp_path/"width-failed-arm"
    width.prepare(directory, nf=128, duration=0.002, step_cap=0.002)
    width.predict(directory)
    def fail(pair):
        raise ValueError("retained own-source homotopy blocker")
    monkeypatch.setattr(width.family, "solve_prepared", fail)
    result = width.measure(directory)
    assert result["status"] == "initial_preparation_failed"
    assert "homotopy blocker" in result["error"]
    assert result["non_existence_claimed"] is False
    assert width.check(directory)["ok"]


def test_cli_default_is_short_preview_and_check_is_readonly(tmp_path):
    directory = tmp_path/"cli-preview"
    assert cli_main(["--prepare", "--output", str(directory)]) == 0
    record = json.loads((directory/"prepare.json").read_text())
    assert record["inputs"]["duration"] == 0.01
    assert record["inputs"]["nf"] == 128
    assert record["inputs"]["production"] is False
    assert cli_main(["--check", "--output", str(directory)]) == 0
    assert cli_main(["--prepare", "--output", str(directory)]) == 2


def test_failed_baseline_preparation_is_sealed_with_source_and_without_false_claim(tmp_path, monkeypatch):
    directory = tmp_path/"baseline-failed"
    def fail(pair):
        raise ValueError("baseline own-source constructor blocked")
    monkeypatch.setattr(width.family, "solve_prepared", fail)
    record = width.prepare(directory, nf=128, duration=0.002)
    assert record["status"] == "initial_preparation_failed"
    assert record["non_existence_claimed"] is False
    _record, arrays = width.load(directory, "prepare")
    assert arrays["source_phi0"].shape == (128, 6)
    assert "initial_r" not in arrays
    assert width.check(directory)["ok"]
    with pytest.raises(ValueError, match="preparation stopped"):
        width.predict(directory)


def test_sealed_stage_cannot_trigger_duplicate_numerical_work(tmp_path, monkeypatch):
    directory = tmp_path/"no-duplicate"
    width.prepare(directory, nf=128, duration=0.002, step_cap=0.002)
    width.predict(directory)
    def forbidden(*args, **kwargs):
        raise AssertionError("duplicate propagation attempted")
    monkeypatch.setattr(width.prediction, "coupled_rk4_step", forbidden)
    with pytest.raises(FileExistsError):
        width.predict(directory)
    assert not (directory/"measurement.json").exists()


def test_failed_jv_keeps_last_finite_baseline_and_direction(tmp_path, monkeypatch):
    directory = tmp_path/"failed-jv"
    width.prepare(directory, nf=128, duration=0.002, step_cap=0.002)
    actual = width.prediction.coupled_rk4_step
    def invalid(pair, state, tangent, dt):
        result = actual(pair, state, tangent, dt)
        result["tangent"].r[0] = float("nan")
        return result
    monkeypatch.setattr(width.prediction, "coupled_rk4_step", invalid)
    record = width.predict(directory)
    assert record["status"] == "chart_or_numerical_stop"
    assert record["steps"] == 0
    assert "nonfinite" in record["stop"]["error"]
    _record, arrays = width.load(directory, "prediction")
    assert np.array_equal(arrays["baseline_r"], arrays["initial_r"])
    assert np.isfinite(arrays["d0_r"]).all()
    assert width.check(directory)["ok"]


def test_retarded_clock_root_carries_tangent_on_the_same_manufactured_step():
    alpha = .25; w = 1.1; target = .0015
    def stepper(parameter, state, tangent, dt):
        return {"state": state+parameter*dt, "tangent": tangent+dt,
                "tau": (1+alpha*parameter)*dt, "delta_tau": alpha*dt}
    bracket = {"state": np.array([w*w]), "tangent": np.array([2*w]),
               "time": 0., "tau": 0., "delta_tau": 0., "dt_hi": .003}
    event = width.root_retarded_event(w, bracket, target, stepper=stepper)
    assert event["tangent_carried_on_same_rooted_step"]
    assert event["uses_linear_clock_correction_as_event"] is False
    assert abs(event["tau"]-target) < 1e-8
    assert event["state"][0] == pytest.approx(w*w+w*event["time"], abs=1e-14)
    assert event["tangent"][0] == pytest.approx(2*w+event["time"], abs=1e-14)
    assert event["delta_tau"] == pytest.approx(alpha*event["time"], abs=1e-14)
    retarded = event["tangent"][0]-w*event["delta_tau"]/(1+alpha*w)
    exact = 2*w+target/(1+alpha*w)**2
    assert retarded == pytest.approx(exact, abs=1e-8)
    assert abs(retarded-event["tangent"][0]) > 1e-5


def test_manufactured_curvature_is_difference_of_common_clock_slopes():
    alpha, target = .4, .002
    # y(w,t)=w^2+w*t and tau=(1+alpha*w)*t.
    def slope(w):
        return 2*w+target/(1+alpha*w)**2
    exact_second = 2-2*alpha*target/(1+alpha)**3
    second = []
    for h in width.DERIVATIVE_H:
        second.append(width.curvature_from_slopes([slope(1-h), slope(1+h)], h))
    assert second[-1] == pytest.approx(exact_second, abs=1e-9)
    assert abs(second[1]-exact_second) <= abs(second[0]-exact_second)
    delta = width.NONLINEAR_HELD_OUT_WIDTH-1
    assert 1.03 == width.NONLINEAR_HELD_OUT_WIDTH
    quadratic = 1+target/(1+alpha)+delta*slope(1)+.5*delta**2*second[-1]
    actual = (1+delta)**2+(1+delta)*target/(1+alpha*(1+delta))
    assert abs(actual-quadratic) < 1e-8


def test_nonlinear_common_clock_forecast_precedes_new_holdout_and_preserves_old(tmp_path, monkeypatch):
    old = tmp_path/"linear-history"; new = tmp_path/"nonlinear-successor"
    width.prepare(old, nf=128, duration=.002, step_cap=.002)
    width.predict(old)
    old_hashes = {path.name: width.episode.file_sha256(path) for path in old.iterdir()}
    actual_load = width.load; loaded = []
    def load(output, stage):
        loaded.append((str(output), stage))
        return actual_load(output, stage)
    monkeypatch.setattr(width, "load", load)
    record = width.prepare_curvature(new, baseline=old, nf=128, duration=.002, step_cap=.002)
    assert record["status"] == "prepared" and len(record["branches"]) == 4
    assert record["inputs"]["held_out_width"] == 1.03
    assert record["nominal"]["old_measurement_read"] is False
    assert record["old_width1p05_used_for_fit"] is False
    assert all(branch["initial_direction"]["base_width"] == branch["width"] for branch in record["branches"])
    forecast = width.predict_curvature(new)
    assert forecast["status"] == "predicted"
    assert len(forecast["curvature_coefficients"]) == 2
    assert all(row["target_tau"] == forecast["target_tau"] for row in forecast["neighbors"])
    assert all(abs(row["tau"]-forecast["target_tau"]) < 1e-8 for row in forecast["neighbors"])
    assert all(row["event"]["tangent_carried_on_same_rooted_step"] for row in forecast["neighbors"])
    assert any(abs(row["time"]-.002) > 1e-10 for row in forecast["neighbors"])
    assert all(row["time"] <= .01 for row in forecast["neighbors"])
    sealed_hash = width.episode.file_sha256(new/"curvature-prediction.json")
    actual_pair = width.pair_at_width; constructed = []
    def held_pair(pair, center):
        if center == 1.03:
            assert width.episode.file_sha256(new/"curvature-prediction.json") == sealed_hash
            constructed.append(center)
        return actual_pair(pair, center)
    monkeypatch.setattr(width, "pair_at_width", held_pair)
    result = width.measure_curvature(new)
    assert constructed == [1.03]
    assert result["status"] == "measured_at_common_absolute_tau"
    assert result["quadratic_signed_remainder"] is not None
    assert result["all_green_threshold"] is None
    assert result["W_unchanged"] and result["observer_unchanged"]
    assert result["coordinate_time"] <= .01
    assert result["aggregate_cpu_seconds"] < 300
    assert width.check(new)["ok"]
    assert width.episode.file_sha256(new/"curvature-prediction.json") == sealed_hash
    assert {path.name: width.episode.file_sha256(path) for path in old.iterdir()} == old_hashes
    assert not any(output == str(old) and stage == "measurement" for output,stage in loaded)


def test_nonlinear_refuses_old_evidence_destination():
    with pytest.raises(PermissionError, match="immutable"):
        width._nonlinear_output(width.LEGACY_LINEAR_OUTPUT)


def test_nonlinear_redo_nominal_uses_one_shared_actual_clock_target(tmp_path):
    directory = tmp_path/"redo-nonlinear"
    prepared = width.prepare_curvature(directory, nf=128, duration=.001, step_cap=.001, redo_baseline=True)
    assert prepared["nominal"]["cached"] is False
    forecast = width.predict_curvature(directory)
    assert forecast["status"] == "predicted"
    assert forecast["target_tau"] > 0
    assert forecast["nominal"]["status"] == "coordinate_baseline_reached"
    assert all(row["target_tau"] == forecast["target_tau"] for row in forecast["neighbors"])
    assert all(row["time"] <= .01 for row in forecast["neighbors"])
    assert width.check(directory)["ok"]


def _manufactured_confirmation_reference(directory):
    # An ordinary short nominal run supplies a real clock/state; curvature
    # numbers are intentionally manufactured to test immutability, not fitting.
    linear = directory/"linear"
    width.prepare(linear, nf=128, duration=.001, step_cap=.001)
    width.predict(linear)
    old, saved = width.load(linear, "prediction")
    reference = directory/"manufactured-reference"
    config = dict(old["inputs"], experiment="nonlinear-v2", held_out_width=1.03, old_width1p05_used_for_fit=False)
    config["inputs_sha256"] = width._digest({k:v for k,v in config.items() if k != "inputs_sha256"})
    nominal = old["coefficients"][-1]
    record = {"schema": width.NONLINEAR_SCHEMA, "inputs": config, "status": "predicted",
              "producer_identity": width._pin_identity(), "target_tau": old["target_tau"],
              "nominal": {"readout": nominal["readout"], "first_coefficient": nominal["coefficient_at_equal_tau"]},
              "curvature_coefficients": [{"second_proper_clock_coefficient": 9.55,
                  "linear_forecast": 1.234, "quadratic_forecast": 1.238}],
              "quadratic_forecast_h_indicator": .0000025}
    width.seal(reference, width.CURVATURE_STAGES[1], record, saved)
    measured = {"schema": width.NONLINEAR_SCHEMA, "inputs": config,
                "status": "measured_at_common_absolute_tau", "producer_identity": width._pin_identity(),
                "prediction_curvature_json_sha256": width.episode.file_sha256(reference/"curvature-prediction.json"),
                "observed_effect": -.057968}
    width.seal(reference, width.CURVATURE_STAGES[2], measured, {"manufactured": np.array([1.])})
    return reference


def test_one_confirmation_locks_coefficients_and_uses_only_four_state_arms(tmp_path, monkeypatch):
    reference = _manufactured_confirmation_reference(tmp_path)
    hashes = {path.name: width.episode.file_sha256(path) for path in reference.iterdir()}
    output = tmp_path/"numerical-confirmation"
    locked = width.lock_curvature_confirmation(output, reference)
    assert locked["own_control_sources_constructed"] is False
    assert locked["locked_original_forecast"]["second_coefficient"] == 9.55
    assert locked["inputs"]["cases"] == width.confirmation_cases()
    assert locked["inputs"]["cpu_budget_seconds"] == 120
    digest = width.episode.file_sha256(output/"confirmation-lock.json")
    actual = width.family.solve_prepared; calls = []
    def solve(pair):
        assert width.episode.file_sha256(output/"confirmation-lock.json") == digest
        calls.append((pair.grid.nf, pair.source_metadata["width_scale"]))
        return actual(pair)
    def forbidden(*args, **kwargs):
        raise AssertionError("confirmation reran neighbor tangents")
    monkeypatch.setattr(width.family, "solve_prepared", solve)
    monkeypatch.setattr(width.prediction, "coupled_rk4_step", forbidden)
    result = width.measure_confirmation(output)
    assert calls == [(256, 1.), (256, 1.03), (128, 1.), (128, 1.03)]
    assert result["status"] == "confirmation_completed"
    assert len(result["arms"]) == 4 and len(result["resolution_effects"]) == 2
    assert all(arm["time"] <= .01 and arm["target_tau"] == locked["inputs"]["common_absolute_tau_target"] for arm in result["arms"])
    assert all(arm["own_initial_solver"]["converged"] and not arm["initial_radius_copied"] for arm in result["arms"])
    assert result["locked_original_forecast"] == locked["locked_original_forecast"]
    assert result["original_forecast_recalibrated"] is False
    assert result["new_response_coefficients_computed"] is False
    assert result["numerical_indicators"]["additional_refinement_automatically_requested"] is False
    assert result["numerical_indicators"]["space_effect_indicator"] == pytest.approx(abs(result["reference_observed_effect"]-result["resolution_effects"][1]["measured_effect"]))
    assert width.check(output)["ok"]
    assert width.episode.file_sha256(output/"confirmation-lock.json") == digest
    assert {path.name: width.episode.file_sha256(path) for path in reference.iterdir()} == hashes
    with pytest.raises(FileExistsError):
        width.confirm_curvature(output, reference)


def test_confirmation_budget_stop_retains_checkpoint_without_coefficients_update(tmp_path, monkeypatch):
    reference = _manufactured_confirmation_reference(tmp_path)
    output = tmp_path/"stopped-confirmation"
    locked = width.lock_curvature_confirmation(output, reference)
    monkeypatch.setattr(width, "_budget", lambda *args, **kwargs: (True, 120., 121.))
    result = width.measure_confirmation(output)
    assert result["status"] == "budget_stop"
    assert not result["arms"]
    assert result["locked_original_forecast"] == locked["locked_original_forecast"]
    assert result["failures"]
    assert width.check(output)["ok"]
