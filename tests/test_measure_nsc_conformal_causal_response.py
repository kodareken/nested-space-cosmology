"""Manufactured controls for the causal photon response. No archived grid is evolved."""
import time

import numpy as np
from numpy.testing import assert_allclose
import pytest

from recursive_horizons.nsc_atlas_cosmological_distance import radial_conformal_lightcone
from recursive_horizons.nsc_conformal_tangent_evolution import FIELDS, variational_rk4_step
from recursive_horizons.nsc_spherical_coupling import CauchyState
from recursive_horizons.nsc_spherical_galerkin_coupling import (
    blank_state, build_grid, manufactured_columns,
)
from measure_nsc_conformal_causal_response import (
    COMPLETED, PARTIAL, archived_frame, baseline_reports, collect_cases, compare_checkpoint,
    dump_summary, evolve_case, field_scale_gap, forecast_accepts,
    initial_affine_radius_rate,
    log_redshift_derivative, preflight_jvp, production_nodes, prolonged_value,
    quadrature_index, ray_row, rays_from_events, require_fresh, summary_status,
    temporal_envelope, write_raw,
)


ACTIVE_CPU = [0.0]


@pytest.fixture(autouse=True)
def control_budget():
    started = time.process_time()
    yield
    ACTIVE_CPU[0] += time.process_time() - started
    assert ACTIVE_CPU[0] < 30.0, "causal response manufactured controls exhausted 30 CPU seconds"


def analytic_metric(epsilon):
    def metric(time, position):
        time, position = np.broadcast_arrays(np.asarray(time, dtype=float),
                                              np.asarray(position, dtype=float))
        phase = np.sin(2.0 * np.pi * position / 8.0)
        radius = (2.2 + 0.15 * np.cos(2.0 * np.pi * position / 8.0)
                  + 0.05 * np.sin(time) + epsilon * 0.03 * np.cos(time) * phase)
        scale = (0.7 + 0.04 * np.cos(time) * np.cos(4.0 * np.pi * position / 8.0)
                 + epsilon * 0.02 * np.sin(time + position))
        return {"radius": radius, "Q": scale}
    return metric


def analytic_sample(time, position):
    phase = np.sin(2.0 * np.pi * position / 8.0)
    radius = 2.2 + 0.15 * np.cos(2.0 * np.pi * position / 8.0) + 0.05 * np.sin(time)
    scale = 0.7 + 0.04 * np.cos(time) * np.cos(4.0 * np.pi * position / 8.0)
    return {
        "r": float(radius),
        "Q": float(scale),
        "ur": float(0.03 * np.cos(time) * phase),
        "uQ": float(0.02 * np.sin(time + position)),
    }


def test_redshift_derivative_matches_the_saved_lightcone_convention():
    nodes = np.linspace(0.0, 0.3, 5)
    step = 1.0e-6
    for direction in (-1, 1):
        kwargs = dict(time_nodes=nodes, reception_time=0.3, reception_x=1.0,
                      lookbacks=(0.125,), direction=direction, quadrature_order=32)
        plus = radial_conformal_lightcone(analytic_metric(step), **kwargs)[0]
        minus = radial_conformal_lightcone(analytic_metric(-step), **kwargs)[0]
        finite_difference = (np.log1p(plus["redshift"]) - np.log1p(minus["redshift"])) / (2.0 * step)
        emission_x = 1.0 + direction * 0.125
        receiver = analytic_sample(0.3, 1.0)
        emitter = analytic_sample(0.175, emission_x)
        independent = ((receiver["ur"] / receiver["r"] + receiver["uQ"] / receiver["Q"])
                       - (emitter["ur"] / emitter["r"] + emitter["uQ"] / emitter["Q"]))
        assert log_redshift_derivative(receiver, emitter) == pytest.approx(independent, abs=0.0)
        assert log_redshift_derivative(receiver, emitter) == pytest.approx(finite_difference, abs=1.0e-8)
        receiver.update(time=0.3, x=1.0, N=receiver["r"] * receiver["Q"], dlogN=0.0)
        emitter.update(time=0.175, x=emission_x, N=emitter["r"] * emitter["Q"], dlogN=0.0)
        measured = ray_row(direction, receiver, emitter)
        assert measured["redshift_derivative"] == pytest.approx(
            (plus["redshift"] - minus["redshift"]) / (2.0 * step), abs=1.0e-8,
        )
        assert plus["emission_time"] == pytest.approx(0.175)
        assert plus["emission_x"] == pytest.approx(emission_x)


def test_production_events_are_quadrature_nodes_224_256_and_288():
    nodes = np.arange(2048) * (8.0 / 2048)
    found = production_nodes(nodes)
    assert (found["receiver"], found["emitters"][-1], found["emitters"][1]) == (256, 224, 288)
    assert quadrature_index(nodes, 1.0) == 256
    with pytest.raises(ValueError, match="quadrature node"):
        quadrature_index(nodes, 1.01)


def test_prolonged_sample_uses_the_geometry_row_not_the_nearest_coarse_node():
    grid = build_grid(12, quadrature=48, gauge="conformal")
    values = np.sin(2.0 * np.pi * grid.xi_g / grid.length)
    index = quadrature_index(grid.xi_q, 1.0)
    sampled = prolonged_value(grid.A_g[index], values)
    nearest = int(np.argmin(np.abs(grid.xi_g - 1.0)))
    assert sampled == pytest.approx(np.sin(2.0 * np.pi / 8.0), abs=1.0e-12)
    assert abs(float(grid.xi_g[nearest]) - 1.0) > 1.0e-6
    assert sampled != pytest.approx(float(values[nearest]), abs=1.0e-3)


def test_two_tiny_steps_match_the_variational_map_and_a_stopped_run_has_no_photon_result(monkeypatch):
    grid = build_grid(10, quadrature=40, gauge="conformal")
    state = blank_state(grid, *manufactured_columns(grid.nf))
    phi0, phi1 = manufactured_columns(grid.nf, seed=5)
    tangent = CauchyState(
        **{name: 1.0e-3 * np.sin((index + 1) * 2.0 * np.pi * grid.xi_g / grid.length)
           for index, name in enumerate(("Q", "r", "chi", "p_Q", "p_r", "p_chi"))},
        phi0=1.0e-4 * phi0, phi1=1.0e-4 * phi1,
    )
    report, cache = evolve_case(
        grid, state.copy(), tangent.copy(), dt=0.001, target=0.002, cpu_seconds=20.0, name="tiny",
        diagnostic_times=(0.0, 0.002), sample_events=False, checkpoint_times=(),
    )
    reference_state, reference_tangent = state.copy(), tangent
    for _step in range(2):
        reference_state, reference_tangent = variational_rk4_step(
            grid, reference_state, reference_tangent, 0.001,
        )
    for name in FIELDS:
        assert_allclose(cache["state"][name], getattr(reference_state, name), rtol=1.0e-12, atol=1.0e-12)
        assert_allclose(cache["tangent"][name], getattr(reference_tangent, name), rtol=1.0e-12, atol=1.0e-12)
    assert report["target_reached"] is True
    assert report["forecast"]["accepted"] is True
    assert report["endpoint_shortened"] is False
    assert report["photon_result"] is False
    assert [row["time"] for row in report["diagnostics"]] == [0.0, 0.002]
    assert report["energy_variation_is_hard_gate"] is False
    assert report["u_max_abs_fields"] == list(FIELDS)
    assert len(report["diagnostics"][0]["u_max_abs"]) == len(FIELDS)
    # Windows process_time can stay flat across a sub-quantum manufactured step.
    # Advance the producer's clock explicitly so the allocation branch is deterministic.
    import measure_nsc_conformal_causal_response as producer
    original_step = producer.variational_rk4_step
    calls = {"count": 0}

    def fail_second_step(*args, **kwargs):
        calls["count"] += 1
        if calls["count"] == 2:
            raise ValueError("injected second-stage failure")
        return original_step(*args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(producer, "variational_rk4_step", fail_second_step)
        failed, failed_cache = evolve_case(
            grid, state.copy(), tangent.copy(), dt=0.001, target=0.002,
            cpu_seconds=20.0, name="fault", diagnostic_times=(0.0, 0.001, 0.002),
            sample_events=False, checkpoint_times=(0.001,),
        )
    assert failed["time"] == pytest.approx(0.001)
    assert failed["steps"] == 1
    assert failed["reason"] == "ValueError: injected second-stage failure"
    assert failed["target_reached"] is False
    assert failed["photon_result"] is False
    assert failed_cache["time"] == failed["time"]
    assert set(failed_cache["checkpoints"]) == {"0.001"}
    expected_state, _ = original_step(grid, state, tangent, 0.001)
    assert_allclose(failed_cache["state"]["r"], expected_state.r, rtol=0.0, atol=0.0)
    ticks = {"value": 0.0}

    def advance():
        ticks["value"] += 1.0
        return ticks["value"]

    monkeypatch.setattr(producer.time, "process_time", advance)
    stopped, _stopped_cache = evolve_case(
        grid, state, tangent, dt=0.001, target=0.002, cpu_seconds=1.0e-15, name="stopped",
        diagnostic_times=(0.0, 0.002), sample_events=False, checkpoint_times=(),
    )
    assert stopped["target_reached"] is False
    assert stopped["photon_result"] is False
    assert stopped["rays"] == []
    assert stopped["time"] < 0.002
    assert stopped["reason"] == "cpu_allocation_exhausted"


def test_forecast_gate_and_temporal_envelope_are_not_a_convergence_claim():
    assert forecast_accepts(1.0, 0.01, 1, 600, 2000.0)[0] is True
    accepted, projected = forecast_accepts(10.0, 2.0, 1, 600, 1000.0)
    assert accepted is False
    assert projected == pytest.approx(10.0 + 2.0 * 599.0)
    fine = {"name": "nf512_dt_0.0005", "photon_result": True, "target_reached": True,
            "baseline": {"within_tolerance": True}, "rays": [
        {"direction": -1, "log_redshift_derivative": 0.2, "redshift": 0.001,
         "redshift_derivative": 1.001 * 0.2},
        {"direction": 1, "log_redshift_derivative": -0.5, "redshift": -0.007,
         "redshift_derivative": 0.993 * -0.5},
    ]}
    coarse = {"name": "nf512_dt_0.0010", "photon_result": True, "target_reached": True,
              "baseline": {"within_tolerance": True}, "rays": [
        {"direction": -1, "log_redshift_derivative": 0.23, "redshift": 0.001,
         "redshift_derivative": 1.001 * 0.23},
        {"direction": 1, "log_redshift_derivative": -0.47, "redshift": -0.007,
         "redshift_derivative": 0.993 * -0.47},
    ]}
    envelope = temporal_envelope([coarse, fine])
    assert envelope["convergence_claimed"] is False
    assert envelope["continuum_bound"] is False
    assert [row["absolute_difference"] for row in envelope["rows"]] == pytest.approx([0.03, 0.03])
    assert [row["redshift_derivative"]["absolute_difference"] for row in envelope["rows"]] == pytest.approx(
        [1.001 * 0.03, 0.993 * 0.03],
    )
    assert envelope["rows"][0]["fine_centered_interval"] == pytest.approx([0.17, 0.23])
    assert summary_status([
        {**fine, "target_reached": True}, {**coarse, "target_reached": True},
    ]) == COMPLETED
    missing = {**fine, "target_reached": False, "photon_result": False, "rays": []}
    assert summary_status([missing, {**coarse, "target_reached": True}]) == PARTIAL
    assert temporal_envelope([missing, coarse]) is None
    bad_baseline = {**fine, "baseline": {"within_tolerance": False}}
    assert summary_status([bad_baseline, coarse]) == PARTIAL
    assert temporal_envelope([bad_baseline, coarse]) is None
    assert rays_from_events({"r": 2.0, "Q": 0.5, "ur": 0.0, "uQ": 0.0, "N": 1.0,
                             "dlogN": 0.0, "time": 0.3, "x": 1.0}, {}) == []


def test_preflight_reports_eight_blocks_and_a_smooth_directional_difference():
    grid = build_grid(10, quadrature=40, gauge="conformal")
    state = blank_state(grid, *manufactured_columns(grid.nf))
    tangent = CauchyState(
        **{name: 1.0e-4 * np.cos((index + 1) * 2.0 * np.pi * grid.xi_g / grid.length)
           for index, name in enumerate(("Q", "r", "chi", "p_Q", "p_r", "p_chi"))},
        phi0=np.zeros_like(state.phi0), phi1=np.zeros_like(state.phi1),
    )
    report = preflight_jvp(grid, state, tangent)
    assert report["compose_calls"] == 2
    assert report["h"] == 1.0e-4
    assert report["dense_jacobian"] is False
    assert [item["field"] for item in report["blocks"]] == list(FIELDS)
    assert report["accepted"] is True
    assert max(item["relative"] for item in report["blocks"]) <= report["gate"]


def test_initial_affine_radius_oracle_uses_an_independent_fourier_second_derivative():
    from recursive_horizons.nsc_conformal_tangent_evolution import conformal_jvp
    grid = build_grid(12, quadrature=48, gauge="conformal")
    state = blank_state(grid, *manufactured_columns(grid.nf))
    state.Q *= 1 + 0.04 * np.cos(2*np.pi*grid.xi_g/grid.length)
    tangent = CauchyState(**{name: np.zeros_like(getattr(state, name)) for name in FIELDS})
    tangent.r = 0.03*np.sin(6*np.pi*grid.xi_g/grid.length)
    oracle = initial_affine_radius_rate(grid, state, tangent)
    assert_allclose(conformal_jvp(grid, state, tangent).p_r, oracle, rtol=1e-12, atol=1e-12)
    report = preflight_jvp(grid, state, tangent, initial_affine=True)
    assert report["initial_affine_radius"]["accepted"] is True
    assert report["accepted"] is True
    tangent.Q[0] = 1e-6
    with pytest.raises(ValueError, match="delta_Q"):
        initial_affine_radius_rate(grid, state, tangent)


def test_baseline_gap_and_saved_frame_lookup_use_scale_not_bitwise_equality(tmp_path):
    actual = np.array([1.0, 1.0 + 5.0e-9])
    gap = field_scale_gap(actual, np.array([1.0, 1.0]))
    assert gap["within_tolerance"] is True
    assert field_scale_gap(actual, np.array([1.0, 1.0 + 1.0e-6]))["within_tolerance"] is False
    path = tmp_path / "frames.npz"
    frame = {name: np.full(4, float(index + 1)) for index, name in enumerate(
        ("Q", "r", "chi", "p_Q", "p_r", "p_chi"))}
    frame["phi0"] = np.zeros((3, 6), dtype=complex)
    frame["phi1"] = np.ones((3, 6), dtype=complex)
    arrays = {f"nf512_dt_0.0005_frames_{name}": np.stack((value, value))
              for name, value in frame.items()}
    arrays["nf512_dt_0.0005_frame_times"] = np.array([0.05, 0.3])
    np.savez(path, **arrays)
    with np.load(path) as data:
        loaded = archived_frame(data, "nf512_dt_0.0005", 0.3)
    comparison = compare_checkpoint(loaded, frame)
    assert comparison["within_tolerance"] is True
    assert comparison["bitwise_gate"] is False
    reports = baseline_reports(path, [{
        "report": {"name": "nf512_dt_0.0005"},
        "cache": {"checkpoints": {"0.050": frame}},
    }])
    assert reports[0]["within_tolerance"] is False
    assert reports[0]["missing_checkpoints"] == ["0.300"]
    reports = baseline_reports(path, [{
        "report": {"name": "nf512_dt_0.0005"},
        "cache": {"checkpoints": {"0.050": frame, "0.300": frame}},
    }])
    assert reports[0]["within_tolerance"] is True
    cache_path = tmp_path / "reached.npz"
    record = write_raw(cache_path, [{
        "name": "nf512_dt_0.0005", "time": 0.3, "target_reached": True,
        "state": frame, "tangent": frame,
    }])
    assert record["bytes"] < 2_000_000
    with np.load(cache_path) as saved:
        assert saved["nf512_dt_0.0005_time"] == pytest.approx(0.3)
        assert saved["nf512_dt_0.0005_u_r"].shape == frame["r"].shape


def test_fresh_paths_and_summary_budget(tmp_path):
    target = tmp_path / "summary.json"
    target.write_text("already", encoding="utf-8")
    with pytest.raises(FileExistsError):
        require_fresh(target, "summary")
    with pytest.raises(ValueError, match="parent"):
        require_fresh(tmp_path / "missing" / "summary.json", "summary")
    assert len(dump_summary({"schema": "example", "status": PARTIAL}).encode()) < 24_000
    with pytest.raises(ValueError, match="24"):
        dump_summary({"pad": "x" * 24_000})


def test_two_spawned_workers_evolve_one_manufactured_step_each():
    grid = build_grid(10, quadrature=40, gauge="conformal")
    state = blank_state(grid, *manufactured_columns(grid.nf))
    tangent = CauchyState(**{
        name: np.zeros_like(getattr(state, name)) for name in FIELDS
    })
    packed_state = {name: getattr(state, name) for name in (
        "Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1")}
    packed_tangent = {name: getattr(tangent, name) for name in packed_state}
    payloads = []
    for name in ("worker-a", "worker-b"):
        payloads.append({
            "name": name, "dt": 0.001, "target": 0.001, "cpu_seconds": 20.0,
            "nf": 10, "nq": 40, "threads": 4,
            "diagnostic_times": (0.0, 0.001), "sample_events": False,
            "checkpoint_times": (), "state": packed_state, "tangent": packed_tangent,
        })
    results = collect_cases(payloads, workers=2)
    assert {item["report"]["name"] for item in results} == {"worker-a", "worker-b"}
    for item in results:
        report = item["report"]
        assert report["target_reached"] is True, report.get("reason")
        assert report["steps"] == 1
        assert report["photon_result"] is False
        assert report["requested_threads"] == 4
        numeric_threads = [value for value in report["observed_num_threads"] if value is not None]
        assert numeric_threads
        assert max(numeric_threads) <= 4
        assert item["cache"]["state"]["r"].shape == state.r.shape
