"""Manufactured RK4-stage balance checks. No production trajectory is loaded."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

import derive_nsc_discovery_balance as cli
from recursive_horizons import nsc_discovery_balance as balance
from recursive_horizons import nsc_discovery_episode as episode
from recursive_horizons import nsc_discovery_observables as observer
from recursive_horizons import nsc_discovery_regions as regions
from recursive_horizons import nsc_nested_parent_child as model
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


LAB = Path(__file__).resolve().parents[1]
ROOT = LAB.parent
PYTHON = Path("/Users/admin/Documents/BlackHoles-Infinity/.venv/validation/bin/python")
SEALED = tuple(LAB / relative for relative in balance.SEALED_INPUTS)
DT = 0.02


def _sealed_stat():
    return tuple((path.stat().st_mtime_ns, path.stat().st_size, hashlib.sha256(path.read_bytes()).hexdigest()) for path in SEALED)


def _state(pair):
    base = galerkin.blank_state(pair.grid, pair.source_phi0, pair.source_phi1)
    coordinate = pair.grid.xi_g
    angle = 2.0 * np.pi * coordinate / pair.grid.length
    base.Q = 1.1 + 0.05 * np.cos(angle)
    base.r = 1.2 + 0.04 * np.sin(2.0 * angle)
    base.chi = 0.03 * np.cos(2.0 * angle)
    base.p_Q = 0.07 + 0.03 * np.sin(angle)
    base.p_r = 0.06 * np.cos(2.0 * angle)
    base.p_chi = 0.04 + 0.02 * np.sin(2.0 * angle)
    generator = np.random.default_rng(4701)
    columns, _ = np.linalg.qr(
        generator.normal(size=(2 * pair.grid.nf, 6)) + 1j * generator.normal(size=(2 * pair.grid.nf, 6))
    )
    base.phi0, base.phi1 = columns[:pair.grid.nf], columns[pair.grid.nf:]
    return model.encode_state(pair, base)


def _forbid_production(*_args, **_kwargs):
    raise AssertionError("production handoff loaded")


def _forbid_initial(*_args, **_kwargs):
    raise AssertionError("initial state was solved")


@pytest.fixture(scope="module")
def prepared():
    pair = model.build_pair(64)
    state = _state(pair)
    nodal = model.reconstruct_state(pair, state)
    rate, bundle = galerkin.compose_fine_hamiltonian(pair.grid, nodal)
    return {"pair": pair, "state": state, "nodal": nodal, "rate": rate, "bundle": bundle}


def test_interpreter_is_the_pinned_validation_environment():
    assert Path(sys.executable) == PYTHON


def test_endpoint_schedule_uses_the_owned_caps():
    span = balance.TARGET_TIME - balance.HANDOFF_TIME
    for cap, count in ((0.001, 700), (0.0005, 1400)):
        steps = balance.acceptance_steps(balance.HANDOFF_TIME, balance.TARGET_TIME, cap)
        assert len(steps) == count
        assert abs(sum(steps) - span) < 1e-12
        assert steps[0] == cap
        assert steps[-1] <= cap
        assert abs(balance.HANDOFF_TIME + sum(steps) - balance.TARGET_TIME) < 1e-12


def test_nyquist_gap_is_explicit_and_disjoint_fluxes_cancel():
    grid = galerkin.build_grid(10, gauge="conformal")
    nyquist = np.cos(np.pi * np.arange(grid.nq))
    assert abs(balance._nyquist_symbol(grid)) < 1e-12
    zeros = np.zeros(grid.nq)
    fields = {
        "energy": zeros,
        "flux": nyquist,
        "pressure": zeros,
        "lapse": zeros,
        "slope": zeros,
        "coordinate": zeros,
        "probability": zeros,
        "probability_flux": nyquist,
    }
    partition = balance.partition_rates(grid, fields)
    fluxes = []
    for name in balance.DISJOINT_WINDOWS:
        window = partition["windows"][name]
        rates = window["rates"]
        fluxes.append(rates["fixed_boundary_flux"])
        assert window["piston_work_included"] is False
        assert window["surface_work"] == 0.0
        assert abs(rates["nyquist_gap"] - (-rates["derivative_window_flux"] - rates["fixed_boundary_flux"])) < 1e-8
        assert abs(rates["finite_projection_defect"] - (rates["pointwise_projection_defect"] + rates["nyquist_gap"])) < 1e-8
        assert abs(rates["balance_rate"] - (rates["slope_integral"] + rates["reynolds"])) < 1e-8
        assert abs(rates["balance_rate"]) < 1e-8
    assert abs(partition["windows"]["left"]["rates"]["nyquist_gap"]) > 1.0
    assert abs(sum(fluxes)) < 1e-8
    assert abs(sum(row["outward_pair_sum_proper"] for row in partition["cuts"])) < 1e-8
    constant = np.ones(grid.nq)
    uniform = balance.partition_rates(grid, dict(fields, flux=constant, probability_flux=constant))
    assert all(abs(uniform["windows"][name]["rates"]["fixed_boundary_flux"]) < 1e-8 for name in balance.DISJOINT_WINDOWS)
    assert abs(sum(row["outward_pair_sum_proper"] for row in uniform["cuts"])) < 1e-8


def test_manufactured_slope_identity_includes_the_measured_defect(prepared, monkeypatch):
    monkeypatch.setattr(episode, "load_saved_handoff", _forbid_production)
    monkeypatch.setattr(episode, "load_confirmation_handoff", _forbid_production)
    monkeypatch.setattr(model, "initial_state", _forbid_initial)
    monkeypatch.setattr(galerkin, "solve_initial_radius", _forbid_initial)
    pair = prepared["pair"]
    state = prepared["state"]
    bundle = prepared["bundle"]
    rate = prepared["rate"]
    fine = bundle["fine_state"]
    system = bundle["fine_system"]
    lifted = observer._lifted_rate(pair.grid, rate, bundle)
    slope = regions.analytic_normal_energy_slope(system, fine, lifted, bundle["source"])
    row = balance.instantaneous_balance(pair, state, bundle=bundle, nodal_rate=rate, control_mode="coupled")
    owned = regions.analyze(pair, state, 0.0, bundle=bundle, nodal_rate=rate, control_mode="coupled")
    child = row["windows"]["child"]
    assert child["dx_divisions"] == 1
    assert child["piston_work_included"] is False
    assert child["positive_part_used_as_accounted_flux"] is False
    assert child["positive_part_net_flux_is_gross_lower_bound"] is True
    assert row["multiplicity_applied_once"] is True
    assert row["multiplicity"] == 4 * int(system.kappa)
    assert abs(child["balance_rate"] - owned["normal_energy"]["child_fixed"]["rate"]) < 1e-8
    assert abs(child["finite_projection_defect"] - owned["normal_energy"]["child_fixed"]["finite_projection_defect"]) < 1e-8
    assert abs(child["nyquist_gap"] - owned["normal_energy"]["child_fixed"]["endpoint_identity_gap"]) < 1e-8
    assert abs(child["pointwise_projection_defect"] - owned["normal_energy"]["child_fixed"]["pointwise_projection_defect"]) < 1e-8
    integral = model.interval_integral(pair.grid, slope / system.dx, (1.0, 3.0))
    assert child["balance_rate"] == pytest.approx(integral, abs=1e-8)
    assert abs(child["balance_rate"] - (
        child["fixed_boundary_flux"] + child["pressure"] + child["lapse"] + child["finite_projection_defect"] + child["reynolds"]
    )) < 1e-8
    parent = row["windows"]["parent"]
    assert parent["kind"] == "derived_sum"
    assert parent["pressure"] == pytest.approx(sum(row["windows"][name]["pressure"] for name in balance.PARENT_PARTS), abs=1e-10)
    assert row["windows"]["global"]["fixed_boundary_flux"] == pytest.approx(0.0, abs=1e-8)
    assert abs(row["internal_cut_cancellation_residual"]) < 1e-8
    assert row["rejected_overlap"]["used_in_closure"] is False
    assert row["rejected_overlap"]["double_count_pressure"] == pytest.approx(child["pressure"], abs=1e-8)
    assert row["direct_parent_window"]["used_in_closure"] is False
    chiral = balance.sigma2_chiral_summand(fine.phi0, fine.phi1, system.occupations)
    observer_flux = 2.0 * np.sum(np.imag(np.conjugate(fine.phi0) * fine.phi1) * system.occupations[None, :], axis=1)
    assert np.max(np.abs(chiral - observer_flux)) < 1e-12
    assert row["gradient_cut_integrated"] is False
    assert row["piston_work_on_measurement_cut"] is False
    assert row["child_normal_closure_claimed"] is False
    assert row["global_normal_closure_claimed"] is False


def test_owned_rk4_state_matches_bitwise_and_clocks_keep_the_trapezoid(prepared, monkeypatch):
    monkeypatch.setattr(episode, "load_saved_handoff", _forbid_production)
    monkeypatch.setattr(model, "initial_state", _forbid_initial)
    assert DT > 0.01
    pair = prepared["pair"]
    state = prepared["state"].copy()
    source = balance._source_token(pair)
    before = _sealed_stat()
    token = episode.state_sha256(state)
    updated, step = balance.rk4_stage_ledger(pair, state, DT, control_mode="coupled")
    owned = model.rk4_step(pair, prepared["state"].copy(), DT)
    assert episode.state_sha256(state) == token
    assert balance._source_token(pair) == source
    assert _sealed_stat() == before
    for name in model.STATE_NAMES:
        assert np.array_equal(getattr(updated, name), getattr(owned, name))
    assert step["quadrature"] == "rk4_stage_weights"
    assert step["observation_cadence_used"] is False
    weights = DT / 6 * (
        step["stage_fieldwork_power"][0] + 2 * step["stage_fieldwork_power"][1]
        + 2 * step["stage_fieldwork_power"][2] + step["stage_fieldwork_power"][3]
    )
    assert step["coordinate_fieldwork"]["owned_power"] == pytest.approx(weights, abs=0.0, rel=0.0)
    clocks = step["clocks"]
    assert clocks["locations"] == [1.0, 2.0, 3.0]
    assert clocks["protocol"] == "trapezoid_coordinate_time_per_accepted_step"
    assert clocks["rk_clock_role"] == "comparison"
    expected = [0.5 * DT * (clocks["start_rates"][index] + clocks["end_rates"][index]) for index in range(3)]
    assert clocks["trapezoid_increment"] == pytest.approx(expected, abs=0.0, rel=0.0)
    metrics_start = model.metrics(pair, prepared["state"])["clock_rates"]
    metrics_end = model.metrics(pair, updated)["clock_rates"]
    assert clocks["start_rates"] == pytest.approx(metrics_start, abs=1e-12)
    assert clocks["end_rates"] == pytest.approx(metrics_end, abs=1e-12)
    assert clocks["stage_sampler_minus_protocol_start_max"] < 1e-12
    child = step["windows"]["child"]
    assert child["measured_change"] == pytest.approx(child["accounted_change"] + child["quadrature_closure_error"], abs=1e-12)
    assert step["closure"]["child"]["closure_error_over_physical"] is not None
    assert abs(step["closure"]["child"]["quadrature_closure_error"]) <= step["closure"]["child"]["physical_scale"] + 1e-8
    assert step["windows"]["global"]["fixed_boundary_flux"] == pytest.approx(0.0, abs=1e-8)
    assert step["rejected_overlap"]["used_in_closure"] is False
    assert step["piston_work_on_measurement_cut"] is False


def test_stages_reuse_source_and_slope_matches_independent_endpoint_difference(prepared, monkeypatch):
    from recursive_horizons import nsc_regional_energy_exchange as regional
    pair, state = prepared["pair"], prepared["state"]
    calls = []
    original = model.rates

    def counted(*args, **kwargs):
        calls.append(kwargs.get("return_bundle"))
        return original(*args, **kwargs)

    monkeypatch.setattr(model, "rates", counted)
    monkeypatch.setattr(regional, "matter_ledger", _forbid_production)
    balance.rk4_stage_ledger(pair, state, 0.002)
    assert calls == [True] * 4
    rate, _bundle = original(pair, state, return_bundle=True)
    eps = 1e-6
    forward, _ = balance.energy_stocks(pair, model._combine(state, rate, eps))
    backward, _ = balance.energy_stocks(pair, model._combine(state, rate, -eps))
    instant = balance.instantaneous_balance(pair, state)
    for name in balance.DISJOINT_WINDOWS:
        central = (forward[name]["normal_energy"] - backward[name]["normal_energy"]) / (2 * eps)
        assert instant["windows"][name]["slope_integral"] == pytest.approx(central, rel=1e-7, abs=1e-7)


def test_halfstep_reduces_closure_error_and_global_probability_drift(prepared, monkeypatch):
    monkeypatch.setattr(episode, "load_saved_handoff", _forbid_production)
    pair, state = prepared["pair"], prepared["state"]
    coarse_state, coarse = balance.integrate_segment(pair, state, 0., 0.02, 0.01, verify_first_step=False)
    fine_state, fine = balance.integrate_segment(pair, state, 0., 0.02, 0.005, verify_first_step=False)
    initial_trace = sum(balance.energy_stocks(pair, state)[0][name]["probability"] for name in balance.DISJOINT_WINDOWS)
    for record in (coarse, fine):
        global_window = record["windows"]["global"]
        assert abs(global_window["probability_slope_integral"]) < 1e-12
        assert abs(global_window["probability_boundary_flux"]) < 1e-12
        assert abs(global_window["probability_gross_in_lower_bound"] - global_window["probability_gross_out_lower_bound"]) < 1e-12
    coarse_drift = abs(coarse["windows"]["global"]["probability_end"] - initial_trace)
    fine_drift = abs(fine["windows"]["global"]["probability_end"] - initial_trace)
    assert fine_drift < coarse_drift / 10
    assert abs(fine["windows"]["child"]["quadrature_closure_error"]) < abs(coarse["windows"]["child"]["quadrature_closure_error"]) / 8
    assert episode.state_sha256(coarse_state) != episode.state_sha256(fine_state)


def test_frozen_geometry_keeps_lapse_and_has_no_metric_work(prepared, monkeypatch):
    monkeypatch.setattr(episode, "load_saved_handoff", _forbid_production)
    monkeypatch.setattr(model, "initial_state", _forbid_initial)
    pair = prepared["pair"]
    state = prepared["state"]
    coupled = balance.instantaneous_balance(pair, state, control_mode="coupled")
    frozen = balance.instantaneous_balance(pair, state, control_mode="frozen_geometry")
    assert frozen["coordinate_fieldwork"]["owned_power"] == 0.0
    assert frozen["coordinate_fieldwork"]["control_rate_power"] == 0.0
    assert abs(coupled["coordinate_fieldwork"]["owned_power"]) > 1e-8
    assert frozen["windows"]["child"]["lapse"] == pytest.approx(coupled["windows"]["child"]["lapse"], abs=1e-12)
    assert abs(frozen["windows"]["child"]["lapse"]) > 1e-8
    assert frozen["windows"]["child"]["pressure"] == pytest.approx(0.0, abs=1e-8)
    assert frozen["windows"]["global"]["fixed_boundary_flux"] == pytest.approx(0.0, abs=1e-8)
    updated, step = balance.rk4_stage_ledger(pair, state.copy(), DT, control_mode="frozen_geometry")
    owned = episode.frozen_geometry_step(pair, state.copy(), DT)
    for name in model.STATE_NAMES:
        assert np.array_equal(getattr(updated, name), getattr(owned, name))
    assert step["coordinate_fieldwork"]["owned_power"] == 0.0
    assert abs(step["windows"]["child"]["lapse"]) > 1e-10
    assert step["windows"]["child"]["measured_change"] == pytest.approx(
        step["windows"]["child"]["accounted_change"] + step["windows"]["child"]["quadrature_closure_error"], abs=1e-12,
    )


def test_segment_telescopes_the_identity_above_the_production_step(prepared, monkeypatch):
    monkeypatch.setattr(episode, "load_saved_handoff", _forbid_production)
    monkeypatch.setattr(episode, "load_confirmation_handoff", _forbid_production)
    monkeypatch.setattr(model, "initial_state", _forbid_initial)
    pair = prepared["pair"]
    _final, record = balance.integrate_segment(pair, prepared["state"], 0.0, 0.04, DT, verify_first_step=True)
    assert record["step_count"] == 2
    assert all(step > 0.01 for step in record["steps"])
    assert record["first_step_bitwise_owned_rk4"] is True
    assert record["observation_cadence_used"] is False
    assert record["production_trajectory_stored"] is False
    for name in list(balance.DISJOINT_WINDOWS) + ["parent", "global"]:
        window = record["windows"][name]
        assert window["measured_change"] == pytest.approx(window["accounted_change"] + window["quadrature_closure_error"], abs=1e-9)
        assert abs(window["stock_telescoping_gap"]) < 1e-8
        assert window["positive_part_used_as_accounted_flux"] is False
    assert record["windows"]["parent"]["pressure"] == pytest.approx(
        sum(record["windows"][name]["pressure"] for name in balance.PARENT_PARTS), abs=1e-9,
    )
    assert record["rejected_overlap"]["double_count_pressure"] == pytest.approx(record["windows"]["child"]["pressure"], abs=1e-8)
    assert record["rejected_overlap"]["used_in_closure"] is False
    assert abs(record["internal_cut_cancellation_residual"]) < 1e-8
    assert record["clocks"]["protocol"] == "trapezoid_coordinate_time_per_accepted_step"
    assert len(record["clocks"]["rk_minus_trapezoid"]) == 3


def test_cli_check_is_creation_only_and_leaves_inputs_unchanged(tmp_path, monkeypatch, capsys):
    before = _sealed_stat()
    check = subprocess.run(
        [str(PYTHON), str(ROOT / "scripts/lab.py"), "scripts/derive_nsc_discovery_balance.py", "--check"],
        check=True, capture_output=True, text=True,
    )
    record = json.loads(check.stdout)
    assert record["schema"] == balance.CHECK_SCHEMA
    assert record["evolved"] is False
    assert record["arrays_loaded"] is False
    assert record["sealed_npz_parsed"] is False
    assert record["output_written"] is False
    assert record["launch"]["launched_here"] is False
    assert record["identity"]["rk_weight"] == "dt/6 * (f1 + 2*f2 + 2*f3 + f4)"
    assert record["identity"]["piston_work_on_measurement_cut"] is False
    assert record["cost"]["new_trajectory_measured"] is False
    assert record["cost"]["npz_parsed"] is False
    assert record["cost"]["cases"][0]["step_count"] == 700
    assert record["cost"]["cases"][1]["step_count"] == 1400
    assert record["cost"]["two_h_estimated_seconds"] is not None
    monkeypatch.setattr(cli, "OUTPUT_DIR", tmp_path)
    monkeypatch.setattr(balance, "evolve_frozen_segment", _forbid_production)
    output = tmp_path / "balance.json"
    output.write_text("{}\n")
    with pytest.raises(FileExistsError):
        cli.main(["--run", "--write", str(output)])
    assert output.read_text() == "{}\n"
    assert _sealed_stat() == before
    assert not (LAB / "results/development/nsc-discovery-balance-v1.json").exists()
    assert cli.EPISODE_DIR.name == "nsc-discovery-episode-v1"
    new = tmp_path / "forecast.json"
    assert cli.main(["--check", "--target-time", "3", "--write", str(new)]) == 0
    assert json.loads(new.read_text())["settings"]["target_time"] == 3.
    with pytest.raises(FileExistsError):
        cli.main(["--check", "--write", str(new)])


def test_record_replay_authenticates_and_rejects_a_changed_channel(prepared, tmp_path):
    _state, ledger = balance.integrate_segment(prepared["pair"], prepared["state"], 0., 0.004, 0.002)
    inputs = balance.sealed_input_hashes()
    record = {"schema": balance.SCHEMA, "producer_hashes": balance.producer_hashes(),
              "input_hashes_before": inputs, "input_hashes_after": inputs,
              "cases": [{"case_id": "manufactured", "ledger": ledger}]}
    path = tmp_path / "manufactured.json"
    path.write_text(json.dumps(record))
    checked = balance.verify_record(path)
    assert checked["ok"] is True
    assert checked["trajectory_recomputed"] is False
    record["cases"][0]["ledger"]["windows"]["child"]["pressure"] += 1.
    path.write_text(json.dumps(record))
    checked = balance.verify_record(path)
    assert checked["ok"] is False
    assert any("accounted channels disagree" in problem for problem in checked["problems"])
