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
    assert new_pair.common_k == report["control_k"] == report["original_common_k"] == pair.common_k
    assert report["common_k_preserved"] and report["momentum_correction"]["converged"]
    assert new_pair.source_metadata["kinetic_anchor"] == report["kinetic_anchor"]
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


@pytest.mark.parametrize("control_mode", ["coupled", "frozen_geometry", "source_free"])
def test_selected_stations_are_stored_and_bound_forecast_horizon(tmp_path, monkeypatch, control_mode):
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
                          cpu_budget=10.,stations=(1.,3.),control_mode=control_mode,step_cap=.00025)
    assert parent.STATIONS==(1.,3.,8.,16.,24.)
    manifest=discovery.read_manifest(destination)
    assert report["stations"]==manifest["stations"]==[1.,3.]
    assert manifest["numerical_binding"]["stations"]==[1.,3.]
    saved,arrays=discovery.load_checkpoint(destination,report["cases"][0]["case_id"])
    assert saved["stations"]==saved["numerical_binding"]["stations"]==[1.,3.]
    assert manifest["control_mode"]==saved["control_mode"]==saved["numerical_binding"]["control_mode"]==control_mode
    assert manifest["step_cap_override"]==saved["step_cap"]==saved["numerical_binding"]["step_cap"]==.00025
    if control_mode == "source_free":
        assert np.all(arrays["phi0"]==0) and saved["covariance_rank"]==0
        assert saved["control_preparation"]["control_k"]==pair.common_k
        assert saved["control_preparation"]["constraint_solved"]
    else:
        np.testing.assert_array_equal(arrays["phi0"],state.phi0)
        np.testing.assert_array_equal(arrays["pi_Q"],state.p_Q)
        np.testing.assert_array_equal(arrays["r"],state.r)
        assert saved["prescribed_geometry_control"] == (control_mode=="frozen_geometry")
    np.testing.assert_array_equal(source_arrays["phi0"],state.phi0)
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


@pytest.mark.parametrize("cap", [0.,-1.,float("nan"),float("inf")])
def test_cap_override_is_positive_finite(cap):
    with pytest.raises(ValueError,match="positive and finite"):
        parent.prepare(step_cap=cap)


def test_source_free_refuses_inadmissible_common_k_without_new_motion():
    from dataclasses import replace
    _prepared,pair,state=fixture_case()
    before=state.copy()
    with pytest.raises(ValueError,match="OPEN.*refusing"):
        parent.source_free_reprepare(replace(pair,common_k=-1.),state)
    for name in leading.FIELDS:
        np.testing.assert_array_equal(getattr(state,name),getattr(before,name))


def test_tagged_column_regional_content_sums_to_total_without_rebasing():
    _prepared,pair,state=fixture_case()
    row=parent.scalar_observation(pair,state,0.)
    tags=row["tagged_column_probability"]
    assert tags["column_indices"]==list(range(pair.column_rank))
    assert not tags["observer_rebased"] and not tags["independent_particle_or_energy_claim"]
    assert sum(tags["child"])==pytest.approx(row["child_probability"],abs=1e-12)
    assert sum(tags["parent"])==pytest.approx(row["parent_probability"],abs=1e-12)
    fine=leading.fine_state(pair,state)
    density=(abs(fine.phi0[:,0])**2+abs(fine.phi1[:,0])**2)*pair.weights[0]/pair.grid.dx_q
    expected=parent.extent.real_interval_integral(pair.grid,density,pair.child_interval)
    assert tags["child"][0]==pytest.approx(expected,abs=1e-12)
    sample=parent.scalar_observation_row(pair,state,0.)
    assert sample["tagged_column_probability"]==tags


def test_run_metadata_reports_real_pool_and_manifest_bytes(tmp_path,monkeypatch):
    write_handoff(tmp_path)
    def fake_pool(workers):
        assert workers==1
        return object()
    monkeypatch.setattr(parent,"pool",fake_pool)
    def no_trajectory_run(output,**kwargs):
        kwargs["executor"](kwargs["workers"])
        return dict(discovery.read_manifest(output),status="STEP_LIMIT",evolved=True)
    monkeypatch.setattr(discovery,"run",no_trajectory_run)
    result=parent.run(tmp_path,workers=1,cpu_budget=10.)
    assert result["pure"] is False and result["pool_launched"] is True
    assert result["executor_pools"]==1
    assert result["bytes_written"]==discovery._manifest_path(tmp_path).stat().st_size
    assert "final manifest" in result["bytes_written_scope"]
    assert discovery.read_manifest(tmp_path)["bytes_written"]==result["bytes_written"]


def test_historical_check_qualifies_current_mismatch_without_writing(tmp_path,monkeypatch):
    write_handoff(tmp_path)
    manifest=discovery.read_manifest(tmp_path)
    path=next(iter(manifest["source_hashes"]))
    manifest["source_hashes"]={path:"0"*64}
    manifest["producing_commit"]="1"*40
    discovery._write_json(discovery._manifest_path(tmp_path),manifest)
    seen=[]
    monkeypatch.setattr(prepared_parent.provenance,"resolve_pinned_source_bytes",
        lambda root,name,digest,commit:seen.append((name,digest,commit)))
    before=parent._directory_files(tmp_path)
    report=parent.check(tmp_path)
    assert report["ok"] and not report["current_producer_match"]
    assert report["historical_sources_authenticated"] and "not replayed" in report["check_scope"]
    assert seen==[(path,"0"*64,"1"*40)]
    assert parent._directory_files(tmp_path)==before


def test_completed_run_is_noop_without_healing_old_manifest(tmp_path):
    write_handoff(tmp_path)
    manifest=discovery.read_manifest(tmp_path)
    manifest.update(status="STATION_REACHED",evolved=True,pure=True,bytes_written=0)
    discovery._write_json(discovery._manifest_path(tmp_path),manifest)
    before=parent._directory_files(tmp_path)
    result=parent.run(tmp_path,workers=1,cpu_budget=10.)
    assert result["run_noop"] and not result["pool_launched"] and result["bytes_written"]==0
    assert parent._directory_files(tmp_path)==before


def test_cli_control_cap_and_ancestry_station_preview(capsys):
    from derive_nsc_discovery_parent_episode import main
    assert main(["--stations","1,1.25,2.25,3","--control-mode","frozen_geometry","--step-cap","0.0005"])==0
    record=json.loads(capsys.readouterr().out)
    assert record["stations"]==[1.,1.25,2.25,3.]
    assert record["control_mode"]=="frozen_geometry" and record["step_cap_override"]==.0005
    assert all(case["control_mode"]=="frozen_geometry" and case["step_cap"]==.0005 for case in record["cases"])
    with pytest.raises(SystemExit):
        main(["--run","--control-mode","source_free"])
    with pytest.raises(SystemExit):
        main(["--check","--step-cap","0.0005"])


def test_real_magnetic_parent_offset_outranks_obsolete_default_without_state_change():
    source=parent.LAB/"results/development/nsc-discovery-parent-v1/six-magnetic-nf128-v1/nf128_pop1_signplus"
    before={name:discovery.file_sha256(source/name) for name in ("parent.json","parent.npz")}
    record,arrays=parent.load_parent_record(source)
    actual=0.002621213673327129
    assert record["k"]==actual and record["k_common"]<2e-8
    pair,state=parent.adopt_parent_record(dict(record,common_k=None),arrays)
    assert pair.common_k==actual
    fine=leading.fine_state(pair,state)
    measured=float(np.mean(fine.p_Q**2/fine.r**3))
    assert measured==pytest.approx(actual,rel=1e-10,abs=1e-14)
    spec={"case_id":"real-source-readonly-regression","nf":128,"step_cap":.001,
          "stations":[1.,3.],"parent_k":actual,"control_mode":"coupled"}
    with parent.runtime_adapter():
        header,_payload=parent.build_case_record(spec,pair,state,producer_commit=None,
                                                 source_binding={},input_binding={})
    assert header["common_k"]==header["declared_offset_k"]==header["parent_k"]==actual
    assert header["initial_kinetic_anchor_mean"]==pytest.approx(actual,rel=1e-10)
    assert not header["kinetic_anchor_equals_declared_offset_claim"]
    for name,stored in (("p_Q","pi_Q"),("p_r","pi_r"),("phi0","phi0")):
        np.testing.assert_array_equal(getattr(state,name),arrays[stored])
    old_directory=parent.LAB/"results/development/nsc-discovery-parent-episode-magnetic-t3-v1"
    old,old_arrays=discovery.load_checkpoint(old_directory,"nf128_plus_balanced_dt0.001",0)
    assert old["common_k"]<2e-8 and old["parent_k"]==actual
    restored=parent.reconstruct_parent_pair(old_arrays,dict(old,common_k=None))
    assert restored.common_k==actual
    assert old["common_k"]<2e-8  # The historical record itself is not healed.
    assert before=={name:discovery.file_sha256(source/name) for name in before}


def test_declared_offset_is_not_fitted_to_general_collar_kinetic_mean():
    _prepared,pair,state=fixture_case()
    arrays=parent.arrays_from_parent(pair,state,clocks=np.zeros(3),clock_rates=np.zeros(3))
    record={"nf":16,"parent_k":.002621213673327129,"common_k":None,"k_common":1e-8}
    reconstructed=parent.reconstruct_parent_pair(arrays,record)
    assert reconstructed.common_k==record["parent_k"]
    fine=leading.fine_state(reconstructed,state)
    assert not np.isclose(float(np.mean(fine.p_Q**2/fine.r**3)),reconstructed.common_k)
    assert parent.declared_parent_offset({"common_k":None,"k_common":.25})==(.25,"k_common")


def test_real_completed_parent_continuation_preserves_whole_state_and_stocks(tmp_path,monkeypatch):
    source=parent.LAB/"results/development/nsc-discovery-parent-strong-t1-v1"
    before=parent._directory_files(source)
    old=discovery.read_manifest(source)
    def forbidden(*args,**kwargs):raise AssertionError("continuation cannot prepare/adopt/resign")
    monkeypatch.setattr(parent,"adopt_parent_record",forbidden)
    monkeypatch.setattr(parent,"source_free_reprepare",forbidden)
    monkeypatch.setattr(prepared_parent,"prepare_parent",forbidden)
    # Only the new-uncommitted producer gate is an oracle. Actual predecessor
    # records, producer history, source inputs and array hashes authenticate.
    monkeypatch.setattr(leading.preparation,"_git_hashes",lambda commit,pins:commit)
    destination=tmp_path/"continued"
    report=parent.continue_parent(source,destination,stations=(1.25,2.25,3.),
                                  cpu_budget=20.,producer_commit="test-future-freeze")
    assert len(report["cases"])==2 and report["case_selection"]=="all predecessor cases"
    assert report["stations"]==[.25,.5,1.,1.25,2.25,3.]
    assert report["source_hashes"]==parent.source_hashes()
    assert report["numerical_binding"]["step_cap"]==.001
    assert report["continuation_ancestor_CPU_spent"]==16.0031
    assert report["ancestor_CPU_budget"]==20500.
    assert report["continuation_remaining_cpu_budget"]==20.
    for case in old["cases"]:
        previous,arrays=discovery.load_checkpoint(source,case["case_id"])
        saved,copied=discovery.load_checkpoint(destination,case["case_id"])
        assert saved["array_sha256"]==previous["array_sha256"]
        for name,value in arrays.items():np.testing.assert_array_equal(copied[name],value)
        for key in ("work_ledger","channel_sample","channel_time","control_mode","step_cap","parent_k"):
            assert saved[key]==previous[key]
        assert saved["coordinate_time"]==1. and saved["steps"]==previous["steps"]
        assert saved["source_strength"]==previous["source_metadata"]["source_strength"]
        assert saved["numerical_binding"]["stations"]==report["stations"]
        assert saved["continuation_predecessor"]["array_sha256"]==previous["array_sha256"]
        assert not saved["source_reset"] and not saved["prepare_parent_called"] and not saved["momenta_resigned"]
    assert discovery.read_manifest(destination)["producing_commit"]=="test-future-freeze"
    checked=parent.check(destination)
    assert checked["ok"] and checked["bytes_written"]==0
    with pytest.raises(ValueError,match="enlarge"):
        parent.run(destination,cpu_budget=21.)
    assert parent._directory_files(source)==before


def test_continuation_requires_explicit_remaining_budget_and_cli_binding(capsys,tmp_path):
    from derive_nsc_discovery_parent_episode import main
    with pytest.raises(ValueError,match="explicit positive remaining"):
        parent.continue_parent(tmp_path,tmp_path/"new",producer_commit="future")
    with pytest.raises(SystemExit):
        main(["--continue-from",str(tmp_path),"--output",str(tmp_path/"new")])
    with pytest.raises(SystemExit):
        main(["--continue-from",str(tmp_path),"--output",str(tmp_path/"new"),
              "--cpu-budget","20","--producer-commit","future","--step-cap",".0005"])


def _varied_state(state, q_scale, momentum_shift):
    varied = state.copy()
    varied.Q = np.array(state.Q, dtype=float) * float(q_scale)
    varied.p_Q = np.array(state.p_Q, dtype=float) + float(momentum_shift)
    return varied


def _simulated_advance(pair, state, samples, status):
    """Distinct station snapshots. No RK4 campaign and no stored-result replay."""
    snapshots = []
    reached = []
    for index, (instant, q_scale, momentum_shift) in enumerate(samples, start=1):
        varied = _varied_state(state, q_scale, momentum_shift)
        rates = np.array(leading.clock_rates(pair, varied), dtype=float, copy=True)
        reached.append(float(instant))
        snapshots.append({
            "kind": "station",
            "time": float(instant),
            "state": varied,
            "clocks": rates * float(instant),
            "clock_rates": rates,
            "ledger": {
                "coordinate_fieldwork": float(index), "pressure_work": 0.0, "lapse_work": 0.0,
                "boundary_child": 0.0, "boundary_parent": 0.0,
                "quadrature": "trapezoid_coordinate_time", "control_mode": "coupled",
                "dense_propagator_stored": False, "stocks_are_coordinate_trapezoid": True,
            },
            "channel_sample": {"time": float(instant)},
            "channel_time": float(instant),
            "steps": index,
            "stations_reached": list(reached),
        })
    last = snapshots[-1]
    return {
        "time": last["time"], "steps": last["steps"], "status": status,
        "stations_reached": list(last["stations_reached"]),
        "stop": {"stop_class": status}, "snapshots": snapshots, "observations": [],
    }


def _assert_observation_matches_checkpoint(record, arrays):
    observation = record["observation"]
    assert observation["time"] == pytest.approx(float(record["coordinate_time"]), abs=1e-12)
    assert observation["control_mode"] == record["control_mode"]
    assert observation["energy"]["total"] != pytest.approx(123456.0)
    with parent.runtime_adapter():
        pair = parent.reconstruct_parent_pair(arrays, record)
        state = parent.state_from_arrays(arrays)
        pair, state, info = discovery.resolve_pair(pair, state, backend="fft")
        assert info["state_changed"] is False and info["W_changed"] is False
        fine = leading.fine_state(pair, state)
        anchor = float(np.mean(fine.p_Q ** 2 / fine.r ** 3))
        rates = leading.clock_rates(pair, state)
        energy = leading.energy(pair, state)
    assert observation["actual_kinetic_anchor_mean"] == pytest.approx(anchor, rel=1e-12, abs=1e-12)
    np.testing.assert_allclose(observation["normal_clock_rates"], rates, rtol=1e-12, atol=1e-12)
    for key, value in energy.items():
        assert observation["energy"][key] == pytest.approx(value, rel=1e-10, abs=1e-10)
    return anchor


def _chunk_sha256(directory, case_id, ordinal):
    stem = f"{case_id}-{ordinal:06d}"
    directory = Path(directory)
    return (
        discovery.file_sha256(directory / (stem + ".json")),
        discovery.file_sha256(directory / (stem + ".npz")),
    )


def _handoff_with_stale_observation(directory):
    """Commit the predecessor row inside the immutable handoff. Do not patch it later."""
    prepared, pair, state = fixture_case(rank=2)
    spec = parent.case_specs(prepared, confirm=False)[0]
    with parent.runtime_adapter():
        record, arrays = parent.build_case_record(
            spec, pair, state, producer_commit=None, source_binding=parent.source_hashes(), input_binding={},
        )
        record["observation"] = {
            "time": -1.0, "control_mode": record["control_mode"],
            "normal_clock_rates": [0.0, 0.0, 0.0],
            "energy": {"gravity": 0.0, "field": 0.0, "total": 123456.0},
            "actual_kinetic_anchor_mean": -1.0,
        }
        parent.commit_parent_checkpoint(directory, record, arrays)
    manifest = {
        "schema": parent.SCHEMA,
        "cases": [{"case_id": spec["case_id"]}],
        "numerical_binding": parent.numerical_binding(),
        "source_hashes": parent.source_hashes(),
        "stations": list(parent.STATIONS),
        "producing_commit": "manufactured-current-test-identity",
    }
    discovery._write_json(discovery._manifest_path(directory), manifest)
    return spec["case_id"]


def test_snapshot_observations_match_their_own_state_time_and_continuation(tmp_path, monkeypatch):
    source = tmp_path / "predecessor"
    source.mkdir()
    case_id = _handoff_with_stale_observation(source)
    original_chunk = _chunk_sha256(source, case_id, 0)
    original_advance = discovery.advance_case
    first_samples = ((0.25, 1.02, 0.004), (0.5, 0.97, -0.003), (1.0, 1.05, 0.008))

    def first_advance(*args, **kwargs):
        return _simulated_advance(args[0], args[1], first_samples, "station_reached")

    monkeypatch.setattr(discovery, "advance_case", first_advance)
    with parent.runtime_adapter():
        parent.execute_case(source, case_id, cpu_allowance=30.0, diagnostics=False)
    assert _chunk_sha256(source, case_id, 0) == original_chunk
    anchors = []
    for ordinal, sample in enumerate(first_samples, start=1):
        record, arrays = discovery.load_checkpoint(source, case_id, ordinal)
        assert record["coordinate_time"] == sample[0]
        assert record["status"] == ("station_reached" if ordinal == len(first_samples) else "station_retained")
        assert record["observation"]["time"] != -1.0
        anchors.append(_assert_observation_matches_checkpoint(record, arrays))
    assert len(set(np.round(anchors, decimals=12))) == len(anchors)
    checked = parent.check(source)
    assert checked["ok"] and checked["bytes_written"] == 0 and not checked["pool_launched"]

    # Emulate an older producer attaching an earlier summary to a later state,
    # by publishing a separate immutable TEMP checkpoint. No prior chunk changes.
    terminal, terminal_arrays = discovery.load_checkpoint(source, case_id)
    terminal["observation"] = {"time": -1.0, "control_mode": "frozen_geometry",
        "energy": {"gravity": 0., "field": 0., "total": 123456.}}
    with parent.runtime_adapter():
        parent.commit_parent_checkpoint(source, terminal, terminal_arrays)
    assert _chunk_sha256(source, case_id, 0) == original_chunk
    discovery._ledger_update(source, pilot=1.0, budget=100.0)
    before = parent._directory_files(source)
    monkeypatch.setattr(discovery, "advance_case", original_advance)
    monkeypatch.setattr(leading.preparation, "_git_hashes", lambda commit, pins: commit)
    destination = tmp_path / "continued"
    report = parent.continue_parent(
        source, destination, stations=(1.25, 2.25, 3.0), cpu_budget=20.0,
        producer_commit="test-future-freeze")
    assert parent._directory_files(source) == before
    assert report["bytes_written"] > 0 and report["pool_launched"] is False
    copied, copied_arrays = discovery.load_checkpoint(destination, case_id)
    previous, previous_arrays = discovery.load_checkpoint(source, case_id)
    assert copied["snapshot_kind"] == "parent_continuation_handoff"
    assert copied["coordinate_time"] == previous["coordinate_time"] == 1.0
    assert previous["observation"]["time"] == -1.0
    assert copied["observation"]["time"] == pytest.approx(1.0)
    for name, value in previous_arrays.items():
        np.testing.assert_array_equal(copied_arrays[name], value)
    copied_anchor = _assert_observation_matches_checkpoint(copied, copied_arrays)
    assert copied_anchor == pytest.approx(anchors[-1], rel=1e-12, abs=1e-12)
    continuation_ordinal = int(copied["ordinal"])
    continuation_chunk = _chunk_sha256(destination, case_id, continuation_ordinal)
    later_samples = ((1.25, 1.01, 0.002), (1.5, 0.94, -0.006), (2.0, 1.08, 0.01))

    def later_advance(*args, **kwargs):
        return _simulated_advance(args[0], args[1], later_samples, "step_limit")

    monkeypatch.setattr(discovery, "advance_case", later_advance)
    with parent.runtime_adapter():
        parent.execute_case(destination, case_id, cpu_allowance=30.0, diagnostics=True)
    assert _chunk_sha256(destination, case_id, continuation_ordinal) == continuation_chunk
    later_anchors = []
    for ordinal, sample in enumerate(later_samples, start=continuation_ordinal+1):
        record, arrays = discovery.load_checkpoint(destination, case_id, ordinal)
        assert record["coordinate_time"] == sample[0]
        assert abs(float(record["observation"]["time"]) - 1.0) > 1e-8
        later_anchors.append(_assert_observation_matches_checkpoint(record, arrays))
    assert len(set(np.round(later_anchors, decimals=12))) == len(later_anchors)
    assert parent._directory_files(source) == before
