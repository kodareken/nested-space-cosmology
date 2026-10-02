"""Tiny runner checks. No production evolution and no writes inside the repository."""
import json
import os
from pathlib import Path

import numpy as np
import pytest
from concurrent.futures import Future

from recursive_horizons import nsc_discovery_episode as episode
from recursive_horizons import nsc_nested_parent_child as model
from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from derive_nsc_discovery_episode import main as cli_main


@pytest.fixture(scope="module")
def manufactured():
    pair = model.build_pair(64)
    base = galerkin.blank_state(pair.grid, pair.source_phi0, pair.source_phi1)
    angle = 2 * np.pi * pair.grid.xi_g / pair.grid.length
    base.Q = 1.1 + 0.05 * np.cos(angle)
    base.r = 1.2 + 0.04 * np.sin(2 * angle)
    base.chi = 0.03 * np.cos(2 * angle)
    base.p_Q = 0.07 + 0.03 * np.sin(angle)
    base.p_r = 0.06 * np.cos(2 * angle)
    base.p_chi = 0.04 + 0.02 * np.sin(2 * angle)
    rng = np.random.default_rng(4701)
    columns, _ = np.linalg.qr(rng.normal(size=(2 * pair.grid.nf, 6))
                              + 1j * rng.normal(size=(2 * pair.grid.nf, 6)))
    base.phi0, base.phi1 = columns[:pair.grid.nf], columns[pair.grid.nf:]
    return pair, model.encode_state(pair, base)


def _bytes(state):
    return tuple(np.ascontiguousarray(getattr(state, name)).tobytes() for name in model.STATE_NAMES)


def _write_case(directory, pair, state, *, case_id="manufactured", geometry="evolving"):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    clocks = np.asarray(model.metrics(pair, state)["clock_rates"], dtype=float)
    arrays = {
        "Q": state.Q, "r": state.r, "chi": state.chi,
        "pi_Q": state.p_Q, "pi_r": state.p_r, "pi_chi": state.p_chi,
        "phi0": state.phi0, "phi1": state.phi1,
        "W": np.array(pair.geometry_map, copy=True),
        "source_phi0": np.array(pair.source_phi0, copy=True),
        "source_phi1": np.array(pair.source_phi1, copy=True),
        "observer_columns": np.array(pair.reference_columns, copy=True),
        "source_weights": np.array(pair.weights, copy=True),
        "clock_rates": clocks, "normal_clocks": np.zeros(3),
    }
    record = {
        "case_id": case_id, "nf": pair.grid.nf, "geometry": geometry, "step_cap": 0.001,
        "stations": [1.0, 3.0], "coordinate_time": 0.3, "steps": 0, "stations_reached": [],
        "momentum_representation": episode.CANONICAL_PI, "verify_external_pins": False,
        "initial_state_called": False, "stability_certificate": False,
        "coarse_indices": pair.geometry_coarse_indices.tolist(),
        "child_indices": pair.geometry_child_indices.tolist(),
        "parent_indices": pair.geometry_parent_indices.tolist(),
        "clock_locations": [1.0, 2.0, 3.0],
        "snapshot_kind": "handoff",
        "work_ledger": episode.empty_work_ledger("frozen_geometry" if geometry == "frozen" else "coupled"),
        "dense_propagator_stored": False,
    }
    committed = episode.commit_checkpoint(directory, record, arrays)
    manifest = {"schema": episode.SCHEMA, "stage": 0, "status": "PREPARED", "evolved": False,
                "stations": [1.0, 3.0], "cases": [{"case_id": case_id}],
                "initial_state_called": False, "stability_certificate": False}
    (directory / "manifest.json").write_text(json.dumps(manifest))
    return committed


def test_geometric_restriction_is_stricter_than_field_cfl_and_not_a_certificate(manufactured):
    pair, state = manufactured
    dt, info = episode.step_restriction(pair, state, 1.0)
    measured = np.max(np.abs(np.linalg.eigvals(pair.grid.derivative)))
    assert abs(episode.quadrature_principal_wavenumber(pair.grid) - measured) < 1e-8
    assert info["geometric"]["identity_holds"] is True
    assert info["geometric"]["majorant_used"] is False
    assert info["field_cfl_alone"] is False
    assert info["stability_certificate"] is False
    assert info["zero_residual_gate"] is False
    assert info["max_T_gate"] is False
    assert dt < info["field_only_dt"]
    assert dt == pytest.approx(min(1.0, episode.RK4_HALF_STABILITY / info["combined_omega"]))
    assert episode.residual_blocks_continuation(1.0e6, 1.0e6) is False


def test_diagnostics_keep_curvature_shell_constraints_car_and_energy_separate(manufactured):
    pair, state = manufactured
    report = episode.stability_diagnostics(pair, state, 0.3)
    assert set(("actual_metric_curvature", "chi_shell", "projected_constraints", "CAR",
                "constraint_algebra_residual", "energy", "work", "chart")) <= set(report)
    assert report["collapsed_score"] is None
    assert report["stability_certificate"] is False
    assert report["actual_metric_curvature"]["chi_substituted"] is False
    assert report["chi_shell"]["substituted_into_curvature"] is False
    assert report["projected_constraints"]["gates_continuation"] is False
    assert report["CAR"]["gram_gap"] >= 0
    assert "total" in report["energy"] and "fieldwork_power" in report["work"]
    assert report["chart"]["clamped"] is False
    assert report["chart"]["admissible"] is True
    row = episode.observe(pair, state, 0.3)
    assert row["stability"]["time"] == pytest.approx(0.3)
    assert row["stability"]["stability_certificate"] is False
    if row["observer_module"] == "explicit_adapter":
        assert row["external_observation"] is None
    else:
        assert row["external_observation"]["time"] == pytest.approx(0.3)
        assert row["external_observation"]["chi_proxy_substituted"] is False
    resolved, resolved_state, info = episode.resolve_pair(pair, state, fft=False)
    assert info["backend"] == "dense"
    assert info["fft_applied"] is False
    assert info["initial_state_called"] is False
    assert _bytes(resolved_state) == _bytes(state)
    assert np.array_equal(np.asarray(resolved.geometry_map), np.asarray(pair.geometry_map))
    fft_pair, fft_state, fft_info = episode.resolve_pair(pair, state, fft=True)
    assert fft_info["fft_applied"] is True
    assert fft_info["backend"] == "fft"
    assert fft_info["W_changed"] is False
    assert _bytes(fft_state) == _bytes(state)
    assert np.array_equal(np.asarray(fft_pair.geometry_map), np.asarray(pair.geometry_map))
    _auto_pair, auto_state, auto_info = episode.resolve_pair(pair, state)
    assert auto_info["backend_requested"] == "auto"
    assert auto_info["backend"] == "fft"
    assert auto_info["initializer"] == "dense_galerkin_grid"
    assert auto_info["initial_state_called"] is False
    assert _bytes(auto_state) == _bytes(state)


def test_chart_exit_keeps_last_admissible_state_and_bracket(tmp_path, manufactured):
    pair, state = manufactured
    directory = tmp_path / "chart"
    _write_case(directory, pair, state)

    def fail(pair, current, dt):
        bad = current.copy()
        bad.Q = np.zeros_like(bad.Q)
        raise coupling.PositiveChartExit("Q_left_positive_chart", np.nan, bad)

    result = episode.execute_case(directory, "manufactured", cpu_allowance=1.0e9, max_steps=4,
                                  stepper=fail, cpu_per_step=0.0, backend="dense")
    assert result["status"] == "chart_exit"
    assert result["steps"] == 0
    assert result["state_clamped"] is False
    assert result["physical_instability_claimed"] is False
    record, arrays = episode.load_checkpoint(directory, "manufactured")
    assert record["stop"]["event_bracket"][0] == pytest.approx(0.3)
    assert record["stop"]["event_bracket"][1] > record["stop"]["event_bracket"][0]
    assert record["stop"]["state_clamped"] is False
    assert record["stop"]["physical_instability_claimed"] is False
    assert record["stop"]["failed_state_retained"] is False
    assert np.array_equal(arrays["Q"], state.Q)
    assert float(np.max(np.abs(arrays["Q"]))) > 1.0


def test_budget_and_memory_stop_without_another_step(tmp_path, manufactured):
    pair, state = manufactured
    budget_dir = tmp_path / "budget"
    _write_case(budget_dir, pair, state)
    stopped = episode.execute_case(budget_dir, "manufactured", cpu_allowance=100.0,
                                   forecast_factor=1.5, cpu_per_step=30.0, max_steps=10,
                                   backend="dense")
    assert stopped["status"] == "budget_stop"
    assert stopped["steps"] == 2
    record, arrays = episode.load_checkpoint(budget_dir, "manufactured")
    assert record["steps"] == 2
    assert record["stop"]["stop_class"] == "budget_stop"
    assert record["coordinate_time"] > 0.3
    assert np.array_equal(arrays["W"], np.array(pair.geometry_map))

    memory_dir = tmp_path / "memory"
    _write_case(memory_dir, pair, state)
    memory = episode.execute_case(memory_dir, "manufactured", cpu_allowance=1.0e9, max_steps=5,
                                  memory_limit_bytes=episode.DEFAULT_MEMORY_BYTES,
                                  rss_bytes=lambda: 10 ** 12, backend="dense")
    assert memory["status"] == "memory_stop"
    assert memory["steps"] == 0
    assert list(memory_dir.glob("manufactured-*.json")) == [memory_dir / "manufactured-000000.json"]


def test_restart_matches_an_uninterrupted_manufactured_run(tmp_path, manufactured):
    pair, state = manufactured
    left = tmp_path / "left"
    right = tmp_path / "right"
    _write_case(left, pair, state)
    _write_case(right, pair, state)
    whole = episode.execute_case(left, "manufactured", cpu_allowance=1.0e9, max_steps=3, backend="dense")
    first = episode.execute_case(right, "manufactured", cpu_allowance=1.0e9, max_steps=2, backend="dense")
    second = episode.execute_case(right, "manufactured", cpu_allowance=1.0e9, max_steps=3, backend="dense")
    assert (whole["steps"], first["steps"], second["steps"]) == (3, 2, 3)
    _record, left_arrays = episode.load_checkpoint(left, "manufactured")
    _record, right_arrays = episode.load_checkpoint(right, "manufactured")
    for name in ("Q", "r", "chi", "pi_Q", "pi_r", "pi_chi", "phi0", "phi1", "W"):
        assert np.array_equal(left_arrays[name], right_arrays[name])
    assert right_arrays["pi_Q"].dtype == np.float64
    assert np.array_equal(right_arrays["observer_columns"], pair.reference_columns)


def test_frozen_geometry_keeps_the_metric_and_the_same_observer(tmp_path, manufactured):
    pair, state = manufactured
    directory = tmp_path / "frozen"
    _write_case(directory, pair, state, geometry="frozen")
    result = episode.execute_case(directory, "manufactured", cpu_allowance=1.0e9, max_steps=1,
                                  backend="dense")
    assert result["status"] == "step_limit"
    _record, arrays = episode.load_checkpoint(directory, "manufactured")
    assert np.array_equal(arrays["Q"], state.Q)
    assert np.array_equal(arrays["r"], state.r)
    assert np.array_equal(arrays["chi"], state.chi)
    assert np.array_equal(arrays["pi_Q"], state.p_Q)
    assert np.array_equal(arrays["pi_r"], state.p_r)
    assert np.array_equal(arrays["pi_chi"], state.p_chi)
    assert not np.array_equal(arrays["phi0"], state.phi0)
    assert np.array_equal(arrays["observer_columns"], pair.reference_columns)
    assert np.array_equal(arrays["source_weights"], pair.weights)
    assert not np.shares_memory(arrays["clock_rates"], arrays["source_weights"])
    assert not np.shares_memory(arrays["phi0"], arrays["source_phi0"])
    loaded_pair = episode.pair_from_arrays(arrays, _record)
    loaded_state = episode.state_from_arrays(arrays)
    frozen_report = episode.stability_diagnostics(
        loaded_pair, loaded_state, _record["coordinate_time"], control_mode="frozen_geometry")
    assert frozen_report["coupled_stability_classification"] is False
    assert frozen_report["actual_metric_curvature"]["Q_dot_max"] == 0.0
    assert frozen_report["actual_metric_curvature"]["Q_ddot_max"] == 0.0
    assert frozen_report["actual_metric_curvature"]["L_dot_max"] == 0.0
    assert frozen_report["work"]["metric_coordinate_work"] == 0.0
    assert frozen_report["work"]["field_rate_max"] > 0.0
    assert "spatial_lapse_work_max" in frozen_report["work"]


def test_checkpoint_guards_immutability_size_and_repository(tmp_path, manufactured):
    pair, state = manufactured
    directory = tmp_path / "guards"
    committed = _write_case(directory, pair, state)
    npz = directory / committed["npz"]
    assert not os.access(npz, os.W_OK)
    with pytest.raises(PermissionError):
        npz.open("ab")
    with pytest.raises(FileExistsError):
        episode._commit_bytes(directory, "manufactured", 0, ".npz", b"changed", episode.CHUNK_LIMIT_BYTES)
    record = {"case_id": "oversized", "momentum_representation": episode.CANONICAL_PI,
              "stability_certificate": False}
    _case, arrays = episode.load_checkpoint(directory, "manufactured")
    with pytest.raises(RuntimeError, match="byte limit"):
        episode.commit_checkpoint(directory, record, arrays, limit=1000)
    assert not (directory / "oversized-000000.npz").exists()
    successor = episode.assert_campaign_output(
        episode.LAB / "results" / "development" / "nsc-discovery-episode-v1")
    assert successor.name == "nsc-discovery-episode-v1"
    with pytest.raises(PermissionError):
        episode.assert_campaign_output(episode.V1_NPZ)
    with pytest.raises(PermissionError):
        episode.prepare(episode.LAB / "results" / "development" / "nsc-spherical-feedback-episode-v1")
    with pytest.raises(PermissionError):
        episode.run(episode.LAB / "src" / "recursive_horizons")


def test_check_does_not_mutate_or_evolve(tmp_path, manufactured, monkeypatch):
    pair, state = manufactured
    directory = tmp_path / "check"
    _write_case(directory, pair, state)

    def forbid(*args, **kwargs):
        raise AssertionError("check evolved")

    monkeypatch.setattr(model, "rk4_step", forbid)
    monkeypatch.setattr(model, "initial_state", forbid)
    monkeypatch.setattr(episode, "advance_case", forbid)
    before = {path.name: (path.stat().st_mode, path.read_bytes())
              for path in directory.iterdir() if path.is_file()}
    report = episode.check(directory)
    after = {path.name: (path.stat().st_mode, path.read_bytes())
             for path in directory.iterdir() if path.is_file()}
    assert report["ok"] is True
    assert report["evolved"] is False
    assert report["bytes_written"] == 0
    assert before == after
    assert cli_main(["--check", "--output", str(directory)]) == 0
    assert before == {path.name: (path.stat().st_mode, path.read_bytes())
                      for path in directory.iterdir() if path.is_file()}


def test_stations_defaults_and_one_pool(tmp_path):
    assert episode.normalize_stations(None) == episode.DEFAULT_STATIONS == (1.0, 3.0)
    assert episode.normalize_stations((8, 16, 24)) == (8.0, 16.0, 24.0)
    assert episode.normalize_stations((1, 3, 8, 16, 24)) == (1.0, 3.0, 8.0, 16.0, 24.0)
    with pytest.raises(ValueError):
        episode.normalize_stations((2.0,))
    with pytest.raises(ValueError):
        episode.normalize_stations((3.0, 1.0))

    class Pool:
        created = []

        def __init__(self, max_workers):
            Pool.created.append(max_workers)

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, traceback):
            return False

        def submit(self, function, payload):
            future = Future()
            future.set_result({
                "case_id": payload["case_id"], "status": "budget_stop", "child_cpu_seconds": 1.0,
                "coordinate_time": 0.3, "steps": 0, "stations_reached": [],
                "stop": {"stop_class": "budget_stop"}, "chunk": None,
                "stability_certificate": False, "state_clamped": False,
                "physical_instability_claimed": False,
            })
            return future

    directory = tmp_path / "pool"
    directory.mkdir()
    manifest = {"schema": episode.SCHEMA, "stage": 0, "status": "PREPARED", "evolved": False,
                "cases": [{"case_id": "one"}, {"case_id": "two"}]}
    (directory / "manifest.json").write_text(json.dumps(manifest))
    result = episode.run(directory, executor=Pool)
    assert Pool.created == [episode.DEFAULT_WORKERS]
    assert result["workers"] == 6
    assert result["cpu_budget_seconds"] == 21600.0
    assert result["forecast_factor"] == 1.5
    assert result["memory_limit_bytes"] == 8 * 1024 ** 3
    assert result["executor_pools"] == 1
    assert result["cpu_accounting"] == "sum_of_child_process_time"
    assert result["child_cpu_seconds"] == pytest.approx(2.0)
    assert result["max_T_forced"] is False
    assert result["zero_residual_gate"] is False
    admitted, ledger = episode._ledger_update(tmp_path / "ledger", pilot=1.0, reserve=50.0, budget=40.0)
    assert admitted is False
    assert ledger["spent"] == pytest.approx(1.0)
    assert ledger["reserved"] == pytest.approx(0.0)


def test_canonical_and_nodal_momentum_round_trips(tmp_path, manufactured):
    pair, state = manufactured
    nodal = model.reconstruct_state(pair, state)
    directory = tmp_path / "nodal"
    directory.mkdir()
    arrays = {
        "Q": state.Q, "r": state.r, "chi": state.chi,
        "nodal_p_Q": nodal.p_Q, "nodal_p_r": nodal.p_r, "nodal_p_chi": nodal.p_chi,
        "dx_g": np.array([pair.grid.dx_g]),
        "phi0": state.phi0, "phi1": state.phi1,
        "W": np.array(pair.geometry_map, copy=True),
        "source_phi0": np.array(pair.source_phi0), "source_phi1": np.array(pair.source_phi1),
        "observer_columns": np.array(pair.reference_columns),
        "source_weights": np.array(pair.weights),
        "clock_rates": np.zeros(3), "normal_clocks": np.zeros(3),
    }
    record = {"case_id": "nodal", "momentum_representation": episode.NODAL_MOMENTUM,
              "stability_certificate": False}
    episode.commit_checkpoint(directory, record, arrays)
    _loaded, stored = episode.load_checkpoint(directory, "nodal")
    restored = episode.state_from_arrays(stored, episode.NODAL_MOMENTUM)
    expected = pair.grid.dx_g * (stored["W"].T @ stored["nodal_p_Q"])
    assert np.array_equal(restored.p_Q, expected)
    canonical = tmp_path / "canonical"
    _write_case(canonical, pair, state, case_id="canonical")
    _record, exact = episode.load_checkpoint(canonical, "canonical")
    assert np.array_equal(exact["pi_Q"], state.p_Q)
    assert np.array_equal(episode.state_from_arrays(exact).p_Q, state.p_Q)
    assert _record["momentum_representation"] == episode.CANONICAL_PI


def test_saved_baseline_handoff_uses_exact_bytes(tmp_path, monkeypatch):
    def forbid(*args, **kwargs):
        raise AssertionError("handoff started a new solve or a new basis")

    monkeypatch.setattr(model, "initial_state", forbid)
    monkeypatch.setattr(model, "build_pair", forbid)
    monkeypatch.setattr(model, "solve_initial_radius_successor", forbid)
    calls = {"rk4": 0}

    def count_rk4(*args, **kwargs):
        calls["rk4"] += 1
        raise AssertionError("prepare or check evolved")

    monkeypatch.setattr(model, "rk4_step", count_rk4)
    directory = tmp_path / "handoff"
    manifest = episode.prepare(directory)
    assert manifest["evolved"] is False
    assert manifest["initial_state_called"] is False
    assert [case["case_id"] for case in manifest["cases"]] == list(episode.baseline_case_ids())
    assert manifest["stations"] == [1.0, 3.0]
    record, arrays = episode.load_checkpoint(directory, "nf128_coupled_dt0.0005")
    with np.load(episode.V1_NPZ, allow_pickle=False) as saved, np.load(episode.BASIS_NPZ, allow_pickle=False) as basis:
        assert np.array_equal(arrays["Q"], saved["nf128_baseline_dt0.0005_Q"][-1])
        assert np.array_equal(arrays["pi_Q"], saved["nf128_baseline_dt0.0005_p_Q"][-1])
        assert np.array_equal(arrays["phi0"], saved["nf128_baseline_dt0.0005_phi0"][-1])
        assert np.array_equal(arrays["phi1"], saved["nf128_baseline_dt0.0005_phi1"][-1])
        assert not np.array_equal(arrays["phi0"], saved["nf128_baseline_dt0.0005_phi0"][0])
        assert np.array_equal(arrays["W"], basis["nf128_W"])
        assert np.array_equal(arrays["observer_columns"], basis["nf128_reference_columns"])
    assert record["momentum_representation"] == episode.CANONICAL_PI
    assert record["source_pins"]["final_state_sha256"]
    assert record["coordinate_time"] == pytest.approx(0.3)
    coarse, coarse_arrays = episode.load_checkpoint(directory, "nf128_coupled_dt0.001")
    assert coarse["parent_case"] == "nf128_baseline_dt0.0005"
    assert coarse["step_cap"] == pytest.approx(0.001)
    assert np.array_equal(coarse_arrays["phi0"], arrays["phi0"])
    with np.load(episode.V1_NPZ, allow_pickle=False) as saved:
        assert not np.array_equal(coarse_arrays["phi0"], saved["nf128_baseline_dt0.001_phi0"][-1])
    frozen, frozen_arrays = episode.load_checkpoint(directory, "nf128_frozen_geometry_dt0.0005")
    assert frozen["geometry"] == "frozen"
    assert frozen["control_mode"] == "frozen_geometry"
    assert np.array_equal(frozen_arrays["phi0"], arrays["phi0"])
    assert np.array_equal(frozen_arrays["observer_columns"], arrays["observer_columns"])
    assert record["source_pins_before"]["W_sha256"] == record["source_pins_after"]["W_sha256"]
    assert record["source_pins_after"]["phi_sha256"] != record["source_pins_after"]["source_columns_sha256"]
    before = {path.name: path.read_bytes() for path in directory.iterdir() if path.is_file()}
    report = episode.check(directory)
    after = {path.name: path.read_bytes() for path in directory.iterdir() if path.is_file()}
    assert report["ok"] is True and report["evolved"] is False
    assert before == after
    assert calls["rk4"] == 0


def test_station_midpoint_retention_and_ledger_restart(tmp_path, manufactured):
    pair, base = manufactured
    state = base.copy()
    state.Q = state.Q + 0.05

    def push(_pair, current, dt):
        out = current.copy()
        out.phi0 = current.phi0 * np.exp(-1j * dt)
        return out

    def run(directory, **extra):
        _write_case(directory, pair, state)
        return episode.execute_case(
            directory, "manufactured", cpu_allowance=1.0e9, backend="dense", stepper=push,
            forced_dt=1.0, observation_cadence=1.0, **extra,
        )

    left, right = tmp_path / "left", tmp_path / "right"
    handoff_bytes = None
    run(left)
    run(right, until_time=1.0)
    handoff_bytes = (right / "manufactured-000000.npz").read_bytes()
    run_again = episode.execute_case(
        right, "manufactured", cpu_allowance=1.0e9, backend="dense", stepper=push,
        forced_dt=1.0, observation_cadence=1.0,
    )
    assert run_again["stations_reached"][-1] == pytest.approx(3.0)
    assert (right / "manufactured-000000.npz").read_bytes() == handoff_bytes
    left_stations = episode.read_station_states(left, "manufactured")
    right_stations = episode.read_station_states(right, "manufactured")
    left_times = [item["coordinate_time"] for item in left_stations]
    right_times = [item["coordinate_time"] for item in right_stations]
    for expected in (0.3, 1.0, 3.0):
        assert any(abs(float(value) - expected) <= 1e-12 for value in left_times)
        assert any(abs(float(value) - expected) <= 1e-12 for value in right_times)
    midpoints = [row["time"] for row in episode.read_observations(left, "manufactured")
                 if 1.0 < float(row["time"]) < 3.0]
    assert midpoints
    left_ledger = episode.read_case_ledger(left, "manufactured")
    right_ledger = episode.read_case_ledger(right, "manufactured")
    for name in episode.WORK_CHANNELS:
        assert left_ledger["work_ledger"][name] == pytest.approx(right_ledger["work_ledger"][name])
    assert left_ledger["normal_clocks"] == pytest.approx(right_ledger["normal_clocks"])
    assert left_ledger["work_quadrature"] == "trapezoid_coordinate_time_between_control_rate_samples"
    assert left_ledger["clock_quadrature"] == "trapezoid_coordinate_time_per_accepted_step"
    audit = episode.read_physical_audit(left, "manufactured")
    assert audit["dense_propagator_stored"] is False
    assert "car_gram_gap" in audit["required_evidence_present"]
    assert "aux_discrepancy_max" in audit["required_evidence_present"]
    assert "localization_width" in audit["required_evidence_present"]
    assert audit["overhead"]["observation_cpu_seconds"] > 0.0
    frozen = tmp_path / "frozen"
    _write_case(frozen, pair, state, geometry="frozen")
    episode.execute_case(
        frozen, "manufactured", cpu_allowance=1.0e9, backend="dense", stepper=push,
        forced_dt=1.0, observation_cadence=1.0, until_time=1.0,
    )
    frozen_ledger = episode.read_case_ledger(frozen, "manufactured")["work_ledger"]
    assert frozen_ledger["coordinate_fieldwork"] == pytest.approx(0.0)
    assert frozen_ledger["control_mode"] == "frozen_geometry"
    assert frozen_ledger["dense_propagator_stored"] is False
    started = __import__("time").process_time()
    episode.evolving_step(pair, state, 1.0e-4)
    step_cpu = __import__("time").process_time() - started
    started = __import__("time").process_time()
    row = episode.scalar_observation_row(pair, state, 0.3, "coupled")
    observation_cpu = __import__("time").process_time() - started
    assert step_cpu > 0.0 and observation_cpu > 0.0
    assert row["dense_propagator_stored"] is False
    assert row["R_h_max"] is not None and row["energy_total"] is not None
    print("MEASURED_OVERHEAD", {"step_cpu_seconds": step_cpu,
                                "observation_cpu_seconds": observation_cpu,
                                "observation_over_step": observation_cpu / step_cpu})


def test_saved_state_two_steps_frozen_jets_and_auto_backend(tmp_path, monkeypatch):
    def forbid(*_args, **_kwargs):
        raise AssertionError("handoff or restart started a new solve")

    monkeypatch.setattr(model, "initial_state", forbid)
    monkeypatch.setattr(model, "build_pair", forbid)
    monkeypatch.setattr(model, "solve_initial_radius_successor", forbid)
    from recursive_horizons.nsc_discovery_backend import operator_backend

    handoff = episode.load_saved_handoff("nf128_baseline_dt0.0005")
    spec = episode._case_record(
        handoff, case_id="nf128_coupled_dt0.0005", control_mode="coupled",
        step_cap=0.0005, stations=(1.0, 3.0),
    )
    arrays = episode._arrays_for_handoff(handoff)
    left, right = tmp_path / "left", tmp_path / "right"
    for directory in (left, right):
        directory.mkdir()
        episode.commit_checkpoint(directory, spec, arrays)
    record, stored = episode.load_checkpoint(left, spec["case_id"])
    dense = episode.pair_from_arrays(stored, record)
    state = episode.state_from_arrays(stored)
    assert operator_backend(dense.grid) == "dense"
    resolved, resolved_state, info = episode.resolve_pair(dense, state, backend="auto")
    assert info["backend"] == "fft"
    assert info["initial_state_called"] is False
    assert np.array_equal(resolved_state.phi0, state.phi0)
    assert np.array_equal(np.asarray(resolved.geometry_map), stored["W"])
    original_chunk = (left / "nf128_coupled_dt0.0005-000000.npz").read_bytes()
    whole = episode.execute_case(left, spec["case_id"], cpu_allowance=1.0e9, max_steps=2, backend="auto")
    episode.execute_case(right, spec["case_id"], cpu_allowance=1.0e9, max_steps=1, backend="auto")
    episode.execute_case(right, spec["case_id"], cpu_allowance=1.0e9, max_steps=2, backend="auto")
    assert whole["steps"] == 2
    assert (left / "nf128_coupled_dt0.0005-000000.npz").read_bytes() == original_chunk
    _left_record, left_arrays = episode.load_checkpoint(left, spec["case_id"])
    _right_record, right_arrays = episode.load_checkpoint(right, spec["case_id"])
    for name in ("Q", "r", "chi", "pi_Q", "pi_r", "pi_chi", "phi0", "phi1",
                 "clock_rates", "normal_clocks", "W", "source_phi0", "source_weights"):
        assert np.array_equal(left_arrays[name], right_arrays[name])
        assert not np.shares_memory(left_arrays[name], right_arrays[name])
    resumed_pair = episode.pair_from_arrays(left_arrays, _left_record)
    resumed_state = episode.state_from_arrays(left_arrays)
    assert np.allclose(left_arrays["clock_rates"], model.metrics(resumed_pair, resumed_state)["clock_rates"])
    assert not np.shares_memory(left_arrays["clock_rates"], left_arrays["source_weights"])
    assert not np.array_equal(left_arrays["phi0"], left_arrays["source_phi0"])
    assert _left_record["source_pins_before"]["source_columns_sha256"] == _left_record["source_pins_after"]["source_columns_sha256"]
    assert _left_record["source_pins_before"]["W_sha256"] == _left_record["source_pins_after"]["W_sha256"]
    assert _left_record["source_pins_before"]["phi_sha256"] != _left_record["source_pins_after"]["phi_sha256"]
    frozen = episode.observe(resolved, state, 0.3, control_mode="frozen_geometry")
    report = frozen["stability"]
    assert report["control_mode"] == "frozen_geometry"
    assert report["coupled_stability_classification"] is False
    assert report["actual_metric_curvature"]["Q_dot_max"] == 0.0
    assert report["actual_metric_curvature"]["Q_ddot_max"] == 0.0
    assert report["actual_metric_curvature"]["L_dot_max"] == 0.0
    assert report["work"]["metric_coordinate_work"] == 0.0
    assert report["work"]["field_rate_max"] > 0.0
    assert report["stability_certificate"] is False
    external = frozen["external_observation"]
    assert external is not None, frozen["notes"]
    assert external["control_mode"] == "frozen_geometry"
    assert external["curvature"]["Q_dot_max"] == 0.0
    assert external["curvature"]["Q_ddot_max"] == 0.0
    assert external["coordinate_metric_work_total"] == 0.0
    assert "control_mode" in frozen["observer_arguments"]
    assert "nodal_rate" in frozen["observer_arguments"]
