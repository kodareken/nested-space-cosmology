"""Magnetic-radius family, 1e-4 admission and the PG adapter gap."""
import inspect
import json
import os
from pathlib import Path

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import numpy as np
import pytest

import derive_nsc_rn_parent as driver


def test_preview_uses_the_magnetic_radius_and_does_not_write(capsys):
    report = driver.preview()
    assert report["mode"] == "preview" and report["read_only"] is True
    assert report["evolved"] is False and report["bytes_written"] == 0
    assert report["production_trajectory"] is False
    assert report["root_review_overrides"] is False
    assert len(report["cases"]) == 6
    ratios = {row["mass_over_rm"] for row in report["cases"]}
    fractions = {row["energy_fraction"] for row in report["cases"]}
    assert ratios == {1.002, 1.01, 1.04}
    assert fractions == {1e-3, 1e-2}
    for row in report["cases"]:
        assert row["magnetic_r2"] == pytest.approx(row["r_m"] ** 2)
        assert row["magnetic_q"] == 1.0 and row["V4"] == 0.0
        assert row["filled_sea"] is False
        assert row["r_m_identified_with_r_minus"] is False
        assert row["r_minus"] != pytest.approx(row["r_m"], rel=0.0, abs=1e-6)
        assert row["excision"] == pytest.approx(0.5 * (row["r_minus"] + row["r_plus"]))
        assert row["route"]["outer"] == pytest.approx(32.0 * row["r_m"])
        assert row["control"]["outer"] == pytest.approx(64.0 * row["r_m"])
        assert row["route"]["child"] == row["control"]["child"]
        assert row["route"]["points"] >= 32 and row["control"]["points"] >= 64
        assert row["route"]["derivative_ladder"] is False
        assert row["launched"] is False
    assert report["configuration"]["calibration_tolerance"] == 1e-4
    assert report["configuration"]["cases"]["r_m_is_inner_horizon"] is False
    assert report["preserved_nf256"]["opened"] is False
    assert report["preserved_nf256"]["mode"] == "read_only"
    code = driver.main([])
    printed = json.loads(capsys.readouterr().out)
    assert code == 0 and printed["mode"] == "preview" and printed["bytes_written"] == 0


def test_calibration_uses_the_live_pg_reconstruction():
    report = driver.calibrate()
    assert report["injected_analytic_jets"] is False
    assert report["bytes_written"] == 0 and report["production_trajectory"] is False
    assert report["mock_physical_output"] is False
    assert report["root_review_overrides"] is False
    assert report.get("measured_residuals") is None or report["calibration_executed"] is True
    if report["calibration_executed"]:
        residuals = report["measured_residuals"]
        assert set(driver.MEASUREMENT_CHANNELS) <= set(residuals)
        within = all(value is not None and value <= 1e-4 for value in residuals.values())
        endpoint_ok = all(row["endpoint_second_derivative"] is None for row in report["cases"])
        assert report["ok"] is (within and endpoint_ok)
        assert report["passed"] is report["ok"]
    else:
        assert report["ok"] is False
        assert report["measured_residuals"] is None
        assert "require_exterior_killing_frequency" in report["blocker"]
    for row in report["cases"]:
        assert row["r_m_identified_with_r_minus"] is False
    assert "sourcefree_rates" in inspect.getsource(driver._reconstruct_sourcefree)
    assert "rk4_step" in inspect.getsource(driver._measure_case)


def test_boolean_cannot_admit_and_run_does_not_invent_a_state(tmp_path):
    forged = {
        "schema": driver.SCHEMA,
        "mode": "calibration",
        "passed": True,
        "ok": True,
        "root_scientific_calibration_passed": True,
        "calibration_executed": True,
        "source_binding": driver.producers(),
        "measured_residuals": {channel: 0.0 for channel in driver.MEASUREMENT_CHANNELS},
    }
    decision = driver.admission_decision(forged)
    assert decision["admitted"] is False
    assert decision["root_review_overrides"] is False
    assert decision["boolean_pass_ignored"] is True
    assert "immutable measurement payload is missing" in decision["blocker"]
    oversized = dict(forged, measured_residuals={channel: 1.0 for channel in driver.MEASUREMENT_CHANNELS})
    refused = driver.admission_decision(oversized)
    assert refused["admitted"] is False
    assert "exceeds 1e-4" in refused["blocker"]
    with pytest.raises(driver.CampaignBlocker) as blocked:
        driver.prepare(tmp_path / "plan", producer_commit="a" * 40, calibration_record=forged)
    assert blocked.value.payload["evolved"] is False
    assert not (tmp_path / "plan").exists()
    with pytest.raises(driver.CampaignBlocker) as gap:
        driver.run(tmp_path / "run", producer_commit="a" * 40, calibration_record=forged)
    assert gap.value.payload["evolution_called"] is False
    assert gap.value.payload["mock_physical_output"] is False
    assert "state" not in gap.value.payload
    assert "immutable measurement payload is missing" in gap.value.payload["blocker"]
    text = inspect.getsource(driver.run)
    assert "prepare_radial_packet" in text and "rk4_step" in text
    assert "PILOT_RANK" in text
    assert driver.PILOT_RANK == 9
    assert "nsc_discovery_parent_cut_confirmation" not in Path(driver.__file__).read_text()


def test_checkpoint_mechanics_stay_creation_only(tmp_path):
    destination = tmp_path / "stage"
    written = driver.write_stage(
        destination, "prepare", {"note": "mechanic"}, {"marker": np.array([1.0])}, blob_match=True,
    )
    assert written["creation_only"] is True and written["evolved"] is False
    with pytest.raises(driver.CampaignBlocker) as exists:
        driver.write_stage(destination, "prepare", {"note": "again"}, blob_match=True)
    assert exists.value.payload["blocker"] == "checkpoint stage already exists"
    replay = driver.resume(destination)
    assert replay["evolution_called"] is False and replay["stages"][0]["stage"] == "prepare"
    original = driver.CHUNK_BYTES
    driver.CHUNK_BYTES = 8
    try:
        with pytest.raises(driver.CampaignBlocker) as oversized:
            driver.write_stage(
                tmp_path / "big", "dependency", {"note": "big"},
                {"marker": np.arange(8.0)}, blob_match=True,
            )
        assert oversized.value.payload["blocker"] == "checkpoint exceeds 64 MiB"
    finally:
        driver.CHUNK_BYTES = original
    with pytest.raises(driver.CampaignBlocker) as sealed:
        driver.write_stage(
            driver.LAB / "results" / "development" / "nsc-rn-parent-should-not-exist",
            "dependency", {}, blob_match=True,
        )
    assert sealed.value.payload["blocker"] == "existing outputs are sealed"


def test_timing_admission_uses_the_factor_without_shrinking_cfl():
    admitted = driver.admit(0.25, 8, remaining=3.0)
    refused = driver.admit(0.25, 8, remaining=2.0)
    assert admitted["forecast_cpu_seconds"] == pytest.approx(3.0)
    assert admitted["admitted"] is True and admitted["admission_factor"] == 1.5
    assert refused["admitted"] is False
    geometry = driver.case_geometry(1.01)
    route = geometry["route"]
    expected = min(geometry["horizon_gap"] / 16.0, route["lambda_min"] / 24.0, route["sigma"] / 24.0)
    assert route["spacing"] == pytest.approx(expected)
    assert route["pending_scales"] == []
    forecast = driver.forecast_first_pair()
    assert forecast["bytes_written"] == 0 and forecast["production_trajectory"] is False
    assert forecast["admitted"] is False and forecast["root_freeze_required"] is True
    assert forecast["mass_over_rm"] == pytest.approx(1.04)
    assert forecast["energy_fraction"] == pytest.approx(1.0e-3)
    assert forecast["steps_to_first_station"] > 0
    assert forecast["within_cpu_budget"] is True and forecast["within_memory_budget"] is True
    assert forecast["stations_reached"] == ()
    assert forecast["step"]["limited_by"] == "step_cap"
    assert forecast["step"]["dt"] == pytest.approx(0.001 * forecast["r_m"])
    assert forecast["cfl_half_step_not_used"] is True
    assert forecast["fixed_control"] == "initial_sourced_metric_frozen"
    assert forecast["measured_pair_step_seconds"] > forecast["measured_coupled_rk4_seconds"]


def test_resume_continues_the_pair_and_the_fixed_metric_stays_frozen(tmp_path):
    from recursive_horizons import nsc_rn_pg as pg
    from recursive_horizons import nsc_rn_source as source
    grid = pg.build_grid(1.04, r_m=1.0, points=33, r_out=16.0)
    phi = np.zeros((2, grid["points"], 1), dtype=complex)
    phi[0, :, 0] = 0.05 * np.exp(1j * grid["radius"])
    phi[1, :, 0] = 0.02 * np.exp(-0.5j * grid["radius"])
    state = pg.make_state(
        grid, phi, occupations=np.array([0.15]), inner_mass=grid["mass"], killing_frequency=1.0,
    )
    initial = pg.stage_rates(state, grid)
    metric = driver.freeze_sourced_metric(initial)
    fixed = driver._clone_state(state)
    coupled = driver._clone_state(state)
    step = driver.matched_step(initial["cfl_dt"], grid["r_m"])
    assert step["dt"] <= step["cfl_dt"]
    calls = {"stage": 0}
    real = pg.stage_rates
    def wrapped(*args, **kwargs):
        calls["stage"] += 1
        return real(*args, **kwargs)
    pg.stage_rates = wrapped
    try:
        first = driver.march_pair(coupled, fixed, grid, metric, dt=step["dt"], steps=2)
    finally:
        pg.stage_rates = real
    assert calls["stage"] == 8
    assert np.array_equal(metric["lapse"], initial["lapse"])
    assert first["fixed"]["inner_mass"] == state["inner_mass"]
    assert abs(first["coupled"]["inner_mass"] - state["inner_mass"]) > 1.0e-8
    assert abs(first["coupled"]["stocks"]["sat_debit"] - first["fixed"]["stocks"]["sat_debit"]) > 0.0
    assert set(first["coupled"]["probability_stocks"]) == {"remaining", "excision_outflow", "outer_outflow"}
    assert "sat_debit" not in first["coupled"]["probability_stocks"]
    observed = driver.station_observation(first["coupled"], grid, real(first["coupled"], grid))
    assert observed["analytic_curvature_injected"] is False
    assert observed["curvature_status"] == "sourced_time_dependent_jets"
    assert np.all(np.isfinite(observed["R4"]))
    assert observed["required_time_jet_primitive"] is None
    assert observed["trapping_horizon"] is not None
    arrays = driver.state_arrays(first["coupled"], first["fixed"], metric)
    written = driver.write_stage(
        tmp_path / "station", "station-probe",
        {"station_over_rm": 0.0, "fixed_control": driver.FIXED_CONTROL},
        arrays, blob_match=True,
    )
    assert written["payload_bytes"] < driver.CHUNK_BYTES
    loaded = driver.resume(tmp_path / "station")
    stored = loaded["arrays"]["station-probe"]
    restored_coupled = driver.apply_checkpoint(state, stored, "coupled")
    restored_fixed = driver.apply_checkpoint(state, stored, "fixed")
    assert np.array_equal(restored_coupled["phi"], first["coupled"]["phi"])
    assert restored_coupled["occupations"] == pytest.approx(state["occupations"])
    assert restored_coupled["clocks"]["t"] == pytest.approx(first["coupled"]["clocks"]["t"])
    assert restored_coupled["stocks"]["sat_debit"] == pytest.approx(first["coupled"]["stocks"]["sat_debit"])
    assert restored_fixed["probability_stocks"]["remaining"] == pytest.approx(
        first["fixed"]["probability_stocks"]["remaining"],
    )
    continued = driver.march_pair(
        restored_coupled, restored_fixed, grid, metric, dt=step["dt"], steps=1,
    )
    direct = driver.march_pair(coupled, fixed, grid, metric, dt=step["dt"], steps=3)
    assert np.allclose(continued["coupled"]["phi"], direct["coupled"]["phi"])
    assert continued["coupled"]["clocks"]["t"] == pytest.approx(direct["coupled"]["clocks"]["t"])
    assert continued["fixed"]["stocks"]["excision_flux"] == pytest.approx(0.0)
    stopped = driver.march_pair(
        driver._clone_state(state), driver._clone_state(state), grid, metric,
        dt=step["dt"], station_time=1.0, cpu_seconds=0.0,
    )
    assert stopped["stopped"] == "cpu_budget" and stopped["steps"] == 1
    packet = source.prepare_radial_packet(
        mass_over_rm=1.04, energy_target=1.0e-3, rank=3, model=grid["model"],
        radii=grid["radius"], weights=grid["weights"], derivative=grid["derivative"],
    )
    admitted = driver.admit_prepared_packet(packet)
    assert admitted["phi_shape"] == [2, grid["points"], 1]
    assert admitted["car_rank"] == 1
    assert admitted["inward_at_center"] is True
    assert admitted["center_probability_current"] < 0.0
    assert admitted["packet_mode_energy"] > 0.0
    assert admitted["killing_energy"] > 0.0


def test_assess_archive_is_read_only_and_reports_the_full_run_gap(tmp_path, capsys):
    archive = Path("/tmp/nsc-rn-pilot-rank9")
    if not (archive / "t0.npz").is_file() or not (archive / "t10.0000.npz").is_file():
        pytest.skip("rank-9 scratch archive is not on this machine")
    before = {path.name: path.stat().st_mtime_ns for path in archive.iterdir()}
    report = driver.assess_archive(archive)
    after = {path.name: path.stat().st_mtime_ns for path in archive.iterdir()}
    assert before == after
    assert report["bytes_written"] == 0 and report["evolution_called"] is False
    assert report["dynamics_loop"] is False
    normal = report["balances"]["normal_dt_0.25"]
    assert normal["full_delta"] == pytest.approx(-0.00153576756387234, abs=1.0e-15)
    assert normal["physical_integral"] == pytest.approx(-0.00153487098034569, abs=1.0e-15)
    assert normal["full_run_gap"] == pytest.approx(normal["full_delta"] - normal["physical_integral"])
    assert abs(normal["full_run_gap"]) == pytest.approx(8.966e-7, rel=1.0e-3, abs=1.0e-12)
    assert normal["interval_max_is_not_full_run_closure"] is True
    assert report["highlights"]["normal_point_max"]["normal_point_index"] == 868
    assert report["highlights"]["normal_point_max"]["t"] == pytest.approx(0.0)
    coordinate = report["balances"]["coordinate_dt_0.25"]
    assert "full_run_gap" in coordinate
    assert report["frames"] == 41
    assert report["context_hashes"]["summary.json"]
    assert report["context_hashes"]["initial-packet.json"]
    assert report["context_hashes"]["producer-hashes-before.json"]
    assert report["context_hashes"]["original_wrapper"]
    assert report["context_hashes"]["parent_receipt"] is None
    assert report["lineage"]["producing_commit"] is None
    assert report["lineage"]["old89_produced_this_run"] is False
    assert "t0.npz" in report["lineage"]["sidecar_payloads_authenticated"]
    assert report["masks"]["normal_continuity"]["nodes_excluded_total"] == 2
    assert report["masks"]["curvature_bulk"]["boundary_width"] == 10
    assert report["masks"]["normal_continuity"]["boundary_width"] == 1
    assert report["masks"]["cells_forced_to_zero"] is False
    assert "wall_seconds" in report and report["cpu_seconds"] >= 0.0
    t10 = next(row for row in report["stations"] if row["tag"] == "t10.0000")
    assert t10["R4_max_abs"] >= t10["R4_bulk_10_max_abs"]
    assert t10["R4_bulk_10_max_abs"] == pytest.approx(9.175393e-7, rel=1.0e-4)
    assert t10["R4_max_abs"] > t10["R4_bulk_10_max_abs"]
    assert driver.main(["--assess-archive", str(archive)]) == 0
    captured = capsys.readouterr().out
    assert "full_run_gap" in captured
    with pytest.raises(driver.CampaignBlocker):
        driver.write_assessment(tmp_path / "assessment.json", report, producer_commit="a" * 40)
    assert not (tmp_path / "assessment.json").exists()


def test_assess_archive_rejects_a_sidecar_hash_mismatch(tmp_path):
    raw = tmp_path / "t0.npz"
    np.savez(raw, radius=np.array([1.0, 2.0]))
    other = tmp_path / "t1.npz"
    np.savez(other, radius=np.array([1.0, 2.0]))
    (tmp_path / "t0.json").write_text(json.dumps({"payload_sha256": "0" * 64}))
    with pytest.raises(driver.CampaignBlocker, match="payload_sha256"):
        driver.assess_archive(tmp_path)


def test_sat_norm_loss_is_twice_the_raw_coefficient_and_closes_the_scratch():
    import json
    from pathlib import Path
    from recursive_horizons import nsc_rn_pg as pg
    from recursive_horizons import nsc_rn_source as source
    radius = np.linspace(4.0, 8.0, 17)
    pack = source.radial_sbp(radius.size, float(radius[0]), float(radius[-1]))
    grid = {
        "radius": radius, "weights": pack["weights"], "derivative": pack["derivative"],
        "spacing": pack["step"], "points": radius.size,
    }
    phi = np.zeros((2, radius.size, 1), dtype=complex)
    phi[1, -1, 0] = 0.4 + 0.15j
    lapse = np.ones(radius.size)
    shift = np.full(radius.size, 0.2)
    _rates, raw = pg.dirac_rates(grid, phi, lapse, shift, kappa=1)
    incoming = -float(shift[-1] + lapse[-1])
    assert raw == pytest.approx(-incoming * abs(phi[1, -1, 0]) ** 2)
    sat_only = np.zeros_like(phi)
    sat_only[1, -1, 0] = (incoming / pack["weights"][-1]) * phi[1, -1, 0]
    norm_rate = 2.0 * np.real(np.sum(
        pack["weights"] * np.conj(phi[:, :, 0]) * sat_only[:, :, 0]
    ))
    assert norm_rate == pytest.approx(-2.0 * raw, abs=1.0e-12)
    step = 1.0e-8
    moved = phi + step * sat_only
    before = pg._probability_norm(phi, pack["weights"])
    after = pg._probability_norm(moved, pack["weights"])
    assert (after - before) / step == pytest.approx(-2.0 * raw, rel=1.0e-5, abs=1.0e-6)
    scratch = Path("/tmp/nsc-rn-pilot-rank9")
    if not (scratch / "t5.0000.json").is_file() or not (scratch / "t10.0000.json").is_file():
        pytest.skip("scratch T5/T10 ledger is not on this machine")
    t0 = json.loads((scratch / "t0.json").read_text())["accounting"]
    initial = t0["coupled_probability"]["remaining"]
    for name in ("t5.0000", "t10.0000"):
        accounting = json.loads((scratch / f"{name}.json").read_text())["accounting"]
        for arm in ("coupled", "fixed"):
            probability = accounting[arm + "_probability"]
            raw_sat = accounting[arm + "_mass_stocks"]["sat_debit"]
            ledger = driver.sat_norm_ledger(raw_sat, probability, initial)
            gap = initial - (
                probability["remaining"] + probability["excision_outflow"]
                + probability["outer_outflow"]
            )
            assert abs(gap - 2.0 * raw_sat) < 9.6e-13
            assert abs(ledger["probability_closure_residual"]) < 1.0e-11
            assert ledger["included_in_mass_flux"] is False
            assert ledger["coefficient_duplicated"] is False
