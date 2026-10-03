"""Parent episode adapter: reconstruction, scaled steps, exact resume, no production run."""
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons import nsc_discovery_episode as discovery
from recursive_horizons import nsc_discovery_backend as backend
from recursive_horizons import nsc_discovery_leading_einstein as leading
from recursive_horizons import nsc_discovery_leading_step_control as step_control
from recursive_horizons import nsc_discovery_parent as prepared_parent
from recursive_horizons import nsc_discovery_parent_episode as parent
from recursive_horizons import nsc_discovery_parent_step_control as parent_step
from recursive_horizons import nsc_spherical_coupling as coupling


def fixture_case(sign=1, population="balanced", rank=2):
    prepared = parent.analytic_fixture(nf=16, rank=rank)
    population_row = next(item for item in prepared["populations"] if item["population_id"] == population)
    pair, state = parent.pair_and_state_from_population(
        population_row, common_k=prepared["common_k"], sign=sign, nf=prepared["nf"], length=prepared["length"],
    )
    return prepared, pair, state


def write_handoff(directory, rank=4):
    prepared, pair, state = fixture_case(rank=rank)
    spec = parent.case_specs(prepared, confirm=False)[0]
    with parent.runtime_adapter():
        record, arrays = parent.build_case_record(
            spec, pair, state, producer_commit=None, source_binding=parent.source_hashes(), input_binding={},
        )
        parent.commit_parent_checkpoint(directory, record, arrays)
    manifest = {
        "schema": parent.SCHEMA,
        "cases": [{"case_id": spec["case_id"]}],
        "numerical_binding": parent.numerical_binding(),
        "source_hashes": parent.source_hashes(),
        "stations": list(parent.STATIONS),
    }
    discovery._write_json(discovery._manifest_path(directory), manifest)
    return spec["case_id"], pair, state


def test_signed_distance_plan_is_pure_and_does_not_hardcode_a_pair_constructor():
    intervals = parent.region_intervals(8)
    assert intervals["centre"] == 4
    assert intervals["child_interval"] == (3.5, 4.5)
    assert intervals["protected_collar"] == (3.0, 5.0)
    assert intervals["parent_interval"] == (1.0, 7.0)
    assert intervals["parent_annulus"] == ((1.0, 2.8), (5.2, 7.0))
    np.testing.assert_allclose(parent.signed_distance([4, 1, 7], 8), [0, -3, 3])
    report = parent.plan()
    assert report["pure"] and report["bytes_written"] == 0 and not report["pool_launched"]
    assert len(report["cases"]) == 6 and len(report["confirmation_cases"]) == 6
    assert {row["nf"] for row in report["cases"]} == {128}
    assert {row["step_cap"] for row in report["cases"]} == {0.001}
    assert {row["nf"] for row in report["confirmation_cases"]} == {256}
    assert {row["step_cap"] for row in report["confirmation_cases"]} == {0.0005}
    assert {row["momentum_sign"] for row in report["cases"]} == {1, -1}
    assert {row["population_id"] for row in report["cases"]} == set(parent.POPULATION_IDS)
    assert report["accepted_parent_record"]["schema"] == "NSC-DISCOVERY-PARENT-v1"
    assert report["accepted_parent_record"]["arrays"] == list(parent.PARENT_RECORD_ARRAYS)
    assert report["stations"] == [1.0, 3.0, 8.0, 16.0, 24.0]
    assert report["numerical_binding"]["raw_coordinate_norm_used"] is False
    assert report["held_out_nonlinear_measurement_before_prediction"] is False
    source = Path(parent.__file__).read_text()
    assert "build_pair(" not in source


def test_momentum_sign_keeps_source_bytes_and_negates_both_momenta():
    prepared = parent.analytic_fixture()
    phi = [item["phi0"].tobytes() for item in prepared["populations"]]
    assert phi[0] == phi[1] == phi[2]
    assert [item["population_id"] for item in prepared["populations"]] == ["child_heavy", "balanced", "parent_heavy"]
    _prepared, plus_pair, plus_state = fixture_case(sign=1)
    _prepared, minus_pair, minus_state = fixture_case(sign=-1)
    assert plus_state.phi0.tobytes() == minus_state.phi0.tobytes()
    assert plus_state.phi1.tobytes() == minus_state.phi1.tobytes()
    assert plus_pair.source_phi0.tobytes() == minus_pair.source_phi0.tobytes()
    np.testing.assert_array_equal(minus_state.p_Q, -plus_state.p_Q)
    np.testing.assert_array_equal(minus_state.p_r, -plus_state.p_r)
    assert np.any(plus_state.p_Q != 0) and np.any(plus_state.p_r != 0)
    signed = dict(prepared["populations"][0], momenta_already_signed=True, pi_Q=-plus_state.p_Q, pi_r=-plus_state.p_r)
    _pair, kept = parent.pair_and_state_from_population(
        signed, common_k=prepared["common_k"], sign=-1, nf=prepared["nf"], length=prepared["length"])
    np.testing.assert_array_equal(kept.p_Q, -plus_state.p_Q)
    assert kept.phi0.tobytes() == plus_state.phi0.tobytes()


def test_reconstruct_parent_pair_restores_cached_arrays(tmp_path):
    case_id, pair, state = write_handoff(tmp_path, rank=4)
    record, arrays = discovery.load_checkpoint(tmp_path, case_id)
    restored = parent.reconstruct_parent_pair(arrays, record)
    revived = parent.state_from_arrays(arrays)
    assert restored.column_rank == restored.covariance_rank == pair.column_rank
    assert restored.weights.shape == (pair.column_rank,)
    assert revived.phi0.shape == (16, pair.column_rank)
    np.testing.assert_array_equal(restored.geometry_map, pair.geometry_map)
    np.testing.assert_array_equal(restored.weights, pair.weights)
    np.testing.assert_array_equal(restored.source_phi0, pair.source_phi0)
    np.testing.assert_array_equal(restored.reference_columns, pair.reference_columns)
    np.testing.assert_array_equal(revived.p_Q, arrays["pi_Q"])
    np.testing.assert_array_equal(revived.p_r, arrays["pi_r"])
    np.testing.assert_array_equal(revived.Q, state.Q)
    assert record["work_ledger"]["coordinate_fieldwork"] == 0.0
    assert np.all(arrays["normal_clocks"] == 0)
    assert record["projector_source_reset"] is False
    with parent.runtime_adapter():
        resolved, same, info = discovery.resolve_pair(restored, revived, backend="fft")
    assert info["fft_applied"] and not info["state_changed"] and not info["W_changed"]
    assert backend.operator_backend(resolved.grid) == "fft"
    np.testing.assert_array_equal(resolved.geometry_map, restored.geometry_map)
    assert all(np.array_equal(getattr(same, name), getattr(revived, name)) for name in leading.FIELDS)


def test_step_restriction_delegates_to_the_shared_owner_and_pins_source():
    _prepared, pair, state = fixture_case(rank=2)
    direct, direct_meta = parent_step.step_restriction(pair, state, 0.001, control_mode="coupled")
    dt, meta = parent.scaled_step_restriction(pair, state, 0.001, control_mode="coupled")
    frozen_dt, frozen_meta = parent.scaled_step_restriction(pair, state, 0.001, control_mode="frozen_geometry")
    assert dt == direct
    assert meta["numerical_mode"] == parent_step.NUMERICAL_MODE == direct_meta["numerical_mode"]
    assert meta["local_real_variables"] == direct_meta["local_real_variables"] == 4 + 4 * pair.column_rank
    assert meta["source_pin"]["column_rank"] == pair.column_rank
    assert meta["source_pin"]["weights_sha256"]
    assert meta["raw_coordinate_norm_used"] is False
    assert frozen_meta["control_mode"] == "frozen_geometry"
    assert 0 < frozen_dt <= 0.001
    assert leading.step_restriction is step_control.RAW_RESTRICTION
    assert "for component in range" not in Path(parent.__file__).read_text()


def test_adapter_restores_and_spawn_initializer_does_not_evolve():
    before = (discovery.SCHEMA, discovery.step_restriction, leading.step_restriction, backend.make_fft_pair,
              discovery.execute_case)
    with pytest.raises(RuntimeError, match="boom"):
        with parent.runtime_adapter():
            assert discovery.SCHEMA == parent.SCHEMA
            assert discovery.step_restriction is parent.scaled_step_restriction
            assert leading.step_restriction is parent.scaled_step_restriction
            assert backend.make_fft_pair is parent.make_fft_pair
            raise RuntimeError("boom")
    assert (discovery.SCHEMA, discovery.step_restriction, leading.step_restriction, backend.make_fft_pair,
            discovery.execute_case) == before
    assert parent.POOL_START_METHOD == "spawn"
    parent.worker_initializer()
    try:
        identity = parent.worker_identity()
    finally:
        parent._worker_scope.__exit__(None, None, None)
        parent._worker_scope = None
    assert identity["schema"] == parent.SCHEMA
    assert identity["scaled_step_installed"] and not identity["raw_step_active"]
    assert identity["fft_pair_installed"] and not identity["trajectory_executed"]
    assert identity["pool_launched_by_initializer"] is False
    assert identity["rates_module"].endswith("nsc_discovery_leading_einstein")
    assert identity["rk4_module"].endswith("nsc_discovery_leading_einstein")
    assert (discovery.SCHEMA, leading.step_restriction, backend.make_fft_pair) == (before[0], before[2], before[3])
    try:
        with ProcessPoolExecutor(max_workers=1, mp_context=multiprocessing.get_context("spawn"),
                                 initializer=parent.worker_initializer) as executor:
            spawned = executor.submit(parent.worker_identity).result(timeout=60)
    except PermissionError as error:
        spawned = error
    if not isinstance(spawned, PermissionError):
        assert spawned["trajectory_executed"] is False and spawned["scaled_step_installed"]
        assert (discovery.SCHEMA, leading.step_restriction) == (before[0], before[2])


def test_rk4_recomputes_source_and_frozen_geometry_jets_are_zero():
    _prepared, pair, state = fixture_case()
    calls = {"n": 0}
    original = coupling.source_from_columns

    def counted(*args, **kwargs):
        calls["n"] += 1
        return original(*args, **kwargs)

    coupling.source_from_columns = counted
    try:
        leading.rk4_step(pair, state.copy(), 1e-4)
    finally:
        coupling.source_from_columns = original
    assert calls["n"] == 4
    assert coupling.source_from_columns is original
    rate = leading.rates(pair, state, control_mode="frozen_geometry")
    for name in leading.REAL_FIELDS:
        assert np.all(getattr(rate, name) == 0)
    assert np.any(rate.phi0 != 0) or np.any(rate.phi1 != 0)
    stepped = leading.frozen_geometry_step(pair, state, 1e-4)
    for name in leading.REAL_FIELDS:
        np.testing.assert_array_equal(getattr(stepped, name), getattr(state, name))
    assert state.phi0.shape[1] == pair.column_rank


def test_source_free_reprepare_keeps_the_caller_field():
    _prepared, pair, state = fixture_case()
    identity = id(state.phi0)
    before = state.phi0.copy()
    new_pair, new_state, report = parent.source_free_reprepare(pair, state)
    assert id(state.phi0) == identity
    np.testing.assert_array_equal(state.phi0, before)
    assert np.all(new_state.phi0 == 0) and np.all(new_state.phi1 == 0)
    np.testing.assert_array_equal(new_pair.weights, pair.weights)
    assert np.any(new_pair.weights > 0)
    assert report["weights_preserved"] and not report["weights_deleted"]
    assert not report["weights_used_as_constraint_solution"]
    assert report["source_reprepared"] and not report["field_deleted_inplace"]
    assert report["intended_empty_source"] and report["constraint_solved"] and report["converged"]
    assert report["leading_constraints"]["raw_C_max"] <= 1e-6
    assert report["leading_constraints"]["D_max"] <= 1e-6
    assert report["leading_constraints"]["source_current_max"] == 0.0
    assert new_pair.source_metadata["population_id"] == pair.source_metadata["population_id"]
    assert new_pair.column_rank == pair.column_rank and new_pair.covariance_rank == 0
    assert report["covariance_rank"] == 0
    assert report["source_column_Gram_distance_from_identity"] == 1
    arrays = parent.arrays_from_parent(new_pair, new_state, clocks=np.zeros(3), clock_rates=np.zeros(3))
    record = {"nf": new_pair.grid.nf, "column_rank": new_pair.column_rank,
              "covariance_rank": 0, "k_common": new_pair.common_k}
    restored = parent.reconstruct_parent_pair(arrays, record)
    assert restored.column_rank == new_pair.column_rank and restored.covariance_rank == 0
    row = parent.scalar_observation(restored, new_state, 0.)
    assert row["covariance_rank"] == 0 and row["source_column_Gram_distance_from_identity"] == 1
    assert new_state is not state and new_pair is not pair


def test_event_checkpoint_resumes_exact_clocks_and_ledger(tmp_path):
    case_id, _pair, _state = write_handoff(tmp_path)
    with parent.runtime_adapter():
        result = parent.execute_case(tmp_path, case_id, cpu_allowance=30.0, max_steps=1)
    assert result["steps"] == 1
    assert result["stop"]["stop_class"] == "step_limit"
    assert result["projector_source_reset"] is False
    first, first_arrays = discovery.load_checkpoint(tmp_path, case_id, 0)
    second, second_arrays = discovery.load_checkpoint(tmp_path, case_id, 1)
    assert second["snapshot_kind"] == "event"
    np.testing.assert_array_equal(second_arrays["source_phi0"], first_arrays["source_phi0"])
    np.testing.assert_array_equal(second_arrays["source_weights"], first_arrays["source_weights"])
    np.testing.assert_array_equal(second_arrays["W"], first_arrays["W"])
    np.testing.assert_array_equal(second_arrays["observer_columns"], first_arrays["observer_columns"])
    assert second["work_ledger"] == discovery.load_checkpoint(tmp_path, case_id, 1)[0]["work_ledger"]
    np.testing.assert_array_equal(second_arrays["normal_clocks"], discovery.load_checkpoint(tmp_path, case_id, 1)[1]["normal_clocks"])
    np.testing.assert_array_equal(second_arrays["pi_Q"], discovery.load_checkpoint(tmp_path, case_id, 1)[1]["pi_Q"])
    assert second["projector_source_reset"] is False
    row = json.loads((tmp_path / f"{case_id}-observations.jsonl").read_text().splitlines()[0])
    assert row["nonlinear_measurement_before_prediction"] is False
    assert row["population_gradient_held_out"] is None
    assert "source_current" in row and "normal_clock_rates" in row
    assert first["coordinate_time"] == 0.0


def test_budget_stop_retains_the_handoff_and_check_is_readonly(tmp_path):
    case_id, _pair, _state = write_handoff(tmp_path)
    payload = {"directory": str(tmp_path), "case_id": case_id, "cpu_budget_seconds": 1e-6,
               "forecast_factor": parent.FORECAST_FACTOR, "memory_limit_bytes": parent.MEMORY_LIMIT_BYTES,
               "max_steps": 1, "fft": None, "backend": "fft"}
    with parent.runtime_adapter():
        result = discovery._case_worker(payload)
    assert result["status"] == "budget_stop"
    assert result["chunk"] is None
    assert discovery._next_ordinal(tmp_path, case_id) == 1
    before = parent._directory_files(tmp_path)
    report = parent.check(tmp_path)
    assert report["ok"] and report["bytes_written"] == 0 and not report["pool_launched"]
    assert parent._directory_files(tmp_path) == before


def test_authenticated_loader_boundary_is_ingested_without_resigning(tmp_path, monkeypatch):
    prepared = parent.analytic_fixture(nf=16, rank=2, common_k=1.0)
    item = prepared["populations"][2]
    phi0 = np.array(item["phi0"], copy=True)
    phi1 = np.array(item["phi1"], copy=True)
    arrays = {
        "Q": np.array(item["Q"], copy=True), "r": np.array(item["r"], copy=True),
        "pi_Q": -np.array(item["pi_Q"], copy=True), "pi_r": -np.array(item["pi_r"], copy=True),
        "phi0": phi0, "phi1": phi1, "W": np.array(item["W"], copy=True),
        "weights": np.array(item["weights"], copy=True),
        "source_columns": np.vstack((phi0, phi1)),
        "reference_columns": np.vstack((phi0, phi1)),
    }
    payload_path = tmp_path / "parent.npz"
    np.savez(payload_path, **arrays)
    payload = payload_path.read_bytes()
    record = {
        "schema": prepared_parent.SCHEMA, "nf": 16, "population": 2, "sign": -1,
        "k": 1.25, "k_common": 1.25, "geometry_map": "identity",
        "intervals": {"child": [3.5, 4.5], "parent": [1.0, 7.0], "annulus": [1.2, 3.0]},
        "clock_locations": list(prepared_parent.CLOCK_LOCATIONS),
        "payload_sha256": __import__("hashlib").sha256(payload).hexdigest(),
        "array_sha256": {key: __import__("hashlib").sha256(np.ascontiguousarray(value).tobytes()).hexdigest()
                         for key, value in arrays.items()},
    }
    (tmp_path / "parent.json").write_text(json.dumps(record))
    # This manufactured record is deliberately not frozen scientific evidence.
    # The actual producer loader must reject its missing authentication pins.
    with pytest.raises(KeyError, match="input_hashes"):
        prepared_parent._load(tmp_path)
    def authenticated_loader_oracle(path):
        assert Path(path) == tmp_path
        return record, {name: value.copy() for name, value in arrays.items()}
    monkeypatch.setattr(prepared_parent, "_load", authenticated_loader_oracle)
    loaded, loaded_arrays = parent.load_parent_record(tmp_path)
    pair, state = parent.adopt_parent_record(loaded, loaded_arrays)
    assert pair.source_metadata.get("population") == 2 or loaded["population"] == 2
    assert state.phi0.tobytes() == phi0.tobytes()
    np.testing.assert_array_equal(state.p_Q, arrays["pi_Q"])
    np.testing.assert_array_equal(state.p_r, arrays["pi_r"])
    assert pair.clock_locations == tuple(prepared_parent.CLOCK_LOCATIONS)
    assert parent.POPULATION_NAME[2] == "parent_heavy"
    with pytest.raises(ValueError, match="6"):
        parent.run(tmp_path, workers=7)
    with pytest.raises(ValueError, match="21600"):
        parent.run(tmp_path, cpu_budget=21601)


def test_explicit_prepare_requires_a_frozen_producer(tmp_path):
    with pytest.raises(ValueError, match="frozen producer"):
        parent.prepare(output=tmp_path / "campaign", execute=True)
    preparation = {
        "schema": parent.PARENT_RECORD_SCHEMA,
        "common_k": 1.0,
        "populations": [{"population_id": name, "population": parent.POPULATION_INDEX[name]}
                        for name in parent.POPULATION_IDS],
    }
    with pytest.raises(ValueError, match="frozen producer"):
        parent.prepare(output=tmp_path / "campaign", execute=True, producer_commit="HEAD", preparation=preparation)
    assert not (tmp_path / "campaign").exists()
    preflight = parent.prepare(output=tmp_path / "campaign", cpu_budget=parent.CPU_BUDGET_SECONDS)
    assert preflight["status"] == "PREFLIGHT" and preflight["bytes_written"] == 0
    assert not (tmp_path / "campaign").exists()


@pytest.mark.parametrize("stations", [(), (0.,1.), (-1.,1.), (1.,float("nan")),
                                      (1.,float("inf")), (3.,1.), (1.,1.)])
def test_station_selection_rejects_invalid_configuration(stations):
    with pytest.raises(ValueError, match="stations"):
        parent.prepare(stations=stations)


def test_selected_stations_are_stored_and_bound_forecast_horizon(tmp_path, monkeypatch):
    _fixture, pair, state = fixture_case()
    source=tmp_path/"source";source.mkdir()
    source_arrays={"Q":state.Q,"r":state.r,"pi_Q":state.p_Q,"pi_r":state.p_r,
                   "phi0":state.phi0,"phi1":state.phi1,"W":pair.geometry_map,
                   "weights":pair.weights,"source_columns":pair.source_columns,
                   "reference_columns":pair.reference_columns}
    np.savez(source/"parent.npz",**source_arrays)
    record={"schema":prepared_parent.SCHEMA,"nf":16,"population":1,"sign":1,
            "k":pair.common_k,"k_common":pair.common_k,"geometry_map":"identity"}
    (source/"parent.json").write_text(json.dumps(record))
    # Narrow authenticated-input boundary oracle, not a fabricated frozen
    # scientific record. Production source authentication is unchanged.
    monkeypatch.setattr(parent,"load_parent_record",lambda path:(record,source_arrays))
    monkeypatch.setattr(leading.preparation,"_git_hashes",lambda commit,pins:commit)
    destination=tmp_path/"campaign"
    report=parent.prepare(source,destination,execute=True,producer_commit="test-oracle",
                          cpu_budget=10.,stations=(1.,3.))
    assert parent.STATIONS==(1.,3.,8.,16.,24.)
    manifest=discovery.read_manifest(destination)
    assert report["stations"]==manifest["stations"]==[1.,3.]
    assert manifest["numerical_binding"]["stations"]==[1.,3.]
    saved,arrays=discovery.load_checkpoint(destination,report["cases"][0]["case_id"])
    assert saved["stations"]==saved["numerical_binding"]["stations"]==[1.,3.]
    np.testing.assert_array_equal(arrays["phi0"],state.phi0)
    np.testing.assert_array_equal(arrays["pi_Q"],state.p_Q)
    estimate=parent.estimate_case_cpu(pair,state,step_cap=.001,stations=saved["stations"],
        current_time=0.,probe_step_cpu=.002,probe_diag_cpu=.001)
    assert estimate["estimated_steps"]==int(np.ceil(3./estimate["dt"]))
    assert estimate["estimated_steps"]<int(np.ceil(24./estimate["dt"]))
    assert parent.plan()["stations"]==[1.,3.,8.,16.,24.]


def test_cli_station_selection_is_prepare_config_not_resume_override(capsys):
    from derive_nsc_discovery_parent_episode import main
    assert main(["--stations","1,3"])==0
    report=json.loads(capsys.readouterr().out)
    assert report["stations"]==[1.,3.]
    assert all(case["stations"]==[1.,3.] for case in report["cases"])
    with pytest.raises(SystemExit):
        main(["--run","--stations","1,3"])
    with pytest.raises(SystemExit):
        main(["--stations","3,1"])
