"""Full coupled RK4 against the causal split on one common-action geometry."""
import os
import tempfile

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_discovery_backend as backend
from recursive_horizons import nsc_discovery_coupled_memory as memory
from recursive_horizons import nsc_discovery_episode as episode
from recursive_horizons import nsc_nested_parent_child as model
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


@pytest.fixture(scope="session", autouse=True)
def _one_thread():
    with threadpool_limits(limits=1), backend.fft_thread_limit(1):
        yield


def _assert_below(value, tolerance, label):
    assert value <= tolerance, (label, value, tolerance)


def _manufactured(nf=64, seed=4):
    pair = model.build_pair(nf)
    nodal = galerkin.blank_state(pair.grid, pair.source_phi0, pair.source_phi1)
    angle = 2.0 * np.pi * pair.grid.xi_g / pair.grid.length
    nodal.Q = nodal.Q * (1.0 + 0.02 * np.sin(angle))
    nodal.r = 1.0 + 0.03 * np.cos(angle)
    nodal.chi = 0.01 * np.sin(2.0 * angle)
    nodal.p_Q = 0.01 * np.cos(angle)
    nodal.p_r = 0.008 * np.sin(angle)
    nodal.p_chi = 0.006 * np.cos(2.0 * angle)
    phi0, phi1 = galerkin.manufactured_columns(nf, seed=seed)
    nodal.phi0 = phi0
    nodal.phi1 = np.exp(-0.5j) * phi1
    state = model.encode_state(pair, nodal)
    pair, state, info = episode.resolve_pair(pair, state, backend="fft")
    assert info["backend"] == "fft"
    assert info["W_changed"] is False
    return pair, state


def test_exact_split_matches_one_full_rk4_step():
    pair, state = _manufactured()
    step = memory.compare_step(pair, state, 1.0e-4)
    for name in ("field_max_abs", "geometry_max_abs", "force_max_abs", "clock_rate_max_abs",
                 "equation_gap", "generator_gap", "energy_max_abs", "energy_field_closure_abs"):
        _assert_below(step[name], 1e-9, name)
    assert step["cross_retained_in_rate"] is True
    assert step["initial_cross_force_max_abs"] > 1e-6
    assert step["initial_cross_frobenius"] > 1e-6
    assert step["modal_kernel_used_as_force"] is False
    assert step["memory_self_energy_used_as_force"] is False
    assert step["one_physical_metric"] is True
    assert step["prescribed_q"] is False
    assert step["omission_used_as_gate"] is False
    assert step["dense_covariance_stored"] is False
    assert step["dense_exterior_propagator_stored"] is False
    assert step["car_full_admissible"] is True
    assert step["car_split_admissible"] is True
    assert step["multiplicity_applied_once"] is True


def test_cross_force_is_not_the_modal_kernel_or_a_dropped_cross():
    pair, state = _manufactured()
    split = memory.initial_split(pair, state)
    report = memory.cross_force_report(pair, split)
    rate = memory.split_rates(pair, split)
    assert report["modal_kernel_shape"] == [2, 2]
    assert report["thin_row_identity_gap"] < 1e-10
    assert report["modal_kernel_used_as_force"] is False
    assert report["memory_self_energy_used_as_force"] is False
    assert report["dense_covariance_formed"] is False
    for name, force in (("force_L", rate.force_L), ("force_Q", rate.force_Q), ("force_beta", rate.force_beta)):
        assert report["forces"][name]["cross_max_abs"] > 1e-6
        assert np.allclose(force, report["full_" + name])
        assert report["forces"][name]["geometry_uses_full_source"] is True
        assert report["forces"][name]["cross_omitted_from_geometry"] is False


def test_elimination_boundary_term_is_not_a_second_source():
    pair, state = _manufactured()
    identity = memory.elimination_variation(pair, state)
    scale = max(1.0, abs(identity["chart_analytic"]), abs(identity["chart_numeric"]))
    _assert_below(identity["chart_gap"], 1e-6 * scale, "chart")
    _assert_below(identity["lapse_partial_gap"], 1e-6 * scale, "lapse")
    _assert_below(identity["shift_partial_gap"], 1e-6 * scale, "shift")
    _assert_below(identity["integration_by_parts_gap"], 1e-8, "by_parts")
    _assert_below(identity["exterior_residual_integral_abs"], 1e-8, "residual")
    _assert_below(identity["force_source_gap"], 1e-10, "force_source")
    assert identity["boundary_abs"] > 1e-3
    assert identity["force_matches_reconstructed_source"] is True
    assert identity["boundary_inserted_into_force"] is False
    assert identity["second_induced_source"] is False
    assert identity["laboratory_is_coupled_production"] is False
    assert identity["geometry_freeze_is_the_physical_law"] is False
    assert identity["dense_exterior_propagator_stored"] is False


def test_omission_is_not_a_unitary_law_or_a_gate():
    pair, state = _manufactured()
    split = memory.initial_split(pair, state)
    assert np.max(np.abs(split.e_memory)) == 0.0
    for name in ("memory", "drive"):
        diagnostic = memory.omission_diagnostic(pair, split, name)
        assert diagnostic["residual_max_abs"] > 1e-4
        assert diagnostic["unitary_global_law"] is False
        assert diagnostic["autonomous_regeneration"] is False
        assert diagnostic["defines_coupled_geometry"] is False
        assert diagnostic["physical_gate"] is False


def test_rank12_control_is_named_and_does_not_zero_the_force():
    with pytest.raises(memory.ImplementationGap):
        memory._require_rank6(np.zeros((4, 12), dtype=np.complex128), np.linspace(0.1, 0.9, 12))
    assert memory.DEPHASING_GAP["status"] == "implementation_gap"
    assert memory.DEPHASING_GAP["force_zeroed"] is False
    assert memory.DEPHASING_GAP["core_modified"] is False
    assert memory.DEPHASING_GAP["implemented"] is False
    assert memory.DEPHASING_GAP["prerequisite_for_coupled_response"] is False


def test_fixed_child_observer_is_not_the_evolved_source():
    pair, state = _manufactured()
    basis = memory.fixed_child_observer(pair)
    evolved = np.vstack((state.phi0, state.phi1))[:, list(memory.CHILD_COLUMN_INDICES)]
    assert np.max(np.abs(basis - evolved)) > 1e-3
    assert np.allclose(basis, pair.original_columns[:, list(memory.CHILD_COLUMN_INDICES)])


@pytest.mark.parametrize("fermions", (128, 256))
def test_stored_short_window_matches_full_coupled_rk4(fermions):
    report = memory.short_coupled_check(fermions, duration=0.01)
    assert report["nf"] == fermions
    assert report["time0"] == 0.3
    assert report["steps"] == 20
    assert report["dt"] == 0.0005
    assert report["production_evolved"] is False
    assert report["clocks_follow_evolved_metric"] is True
    assert report["stored_handoff_rates_reused_as_future_rates"] is False
    for name in ("field_max_abs", "geometry_max_abs", "force_max_abs", "clock_rate_max_abs",
                 "clock_integral_max_abs", "covariance_spectrum_max_abs", "equation_gap",
                 "generator_gap", "energy_max_abs", "energy_field_closure_abs"):
        _assert_below(report[name], 1e-12, name)
    assert report["geometry_motion_max_abs"] > 1e-8
    assert report["field_motion_max_abs"] > 1e-6
    assert report["clock_advance_max_abs"] > 1e-6
    assert report["cross_retained_in_rate"] is True
    assert report["initial_cross_frobenius"] > 1e-3
    assert report["initial_cross_force_max_abs"] > 1e-4
    assert report["car_full_admissible"] is True
    assert report["car_split_admissible"] is True
    assert report["multiplicity_applied_once"] is True
    assert report["one_physical_metric"] is True
    assert report["omission_used_as_gate"] is False
    assert report["dense_exterior_propagator_stored"] is False
    assert report["dephasing_prerequisite"] is False
    assert report["memory_reset_at_endpoint"] is False
    assert report["endpoint_response_uses_evolved_split"] is True
    assert report["endpoint_memory_norm"] > 1e-6
    _assert_below(report["effective_outside_force_max_abs"], 1e-12, "outside")
    _assert_below(report["complementary_force_max_abs"], 1e-12, "complementary")
    _assert_below(report["endpoint_cross_force_max_abs"], 1e-12, "endpoint_cross")
    _assert_below(report["endpoint_chart_gap"], 1e-6, "endpoint_chart")
    _assert_below(report["endpoint_lapse_partial_gap"], 1e-6, "endpoint_lapse")
    _assert_below(report["endpoint_shift_partial_gap"], 1e-6, "endpoint_shift")
    identity = report["elimination"]
    _assert_below(identity["integration_by_parts_gap"], 1e-6, "stored_by_parts")
    _assert_below(identity["exterior_residual_integral_abs"], 1e-8, "stored_residual")
    _assert_below(identity["force_source_gap"], 1e-10, "stored_force")
    assert identity["boundary_abs"] > 1e-3
    assert identity["second_induced_source"] is False
    assert identity["boundary_inserted_into_force"] is False
    assert identity["laboratory_is_coupled_production"] is False
    for name in ("memory", "drive"):
        assert report["omissions"][name]["residual_max_abs"] > 1e-4
        assert report["omissions"][name]["unitary_global_law"] is False
        assert report["omissions"][name]["autonomous_regeneration"] is False
        assert report["omissions"][name]["physical_gate"] is False


def test_dependency_closure_pins_the_stored_handoff_without_production():
    record = memory.dependency_closure()
    assert record["schema"] == memory.SCHEMA
    assert record["production_segment"] == [0.3, 1.0]
    assert record["segment_duration"] == 0.7
    assert record["allowed_step_caps"] == [0.0005, 0.001, 0.00025]
    assert record["checkpoint_limit_bytes"] == 64 * 1024 * 1024
    assert record["dephasing_prerequisite_for_coupled_response"] is False
    assert record["production_evolved"] is False
    assert record["t3_in_scope"] is False
    assert record["short_check_ran"] is False
    assert record["new_interactions"] is False
    assert record["unphysical_control_gate_required"] is False
    assert record["omission_is_unitary_global_law"] is False
    assert record["equations"]["propagator_stored"] is False
    assert [case["nf"] for case in record["cases"]] == [128, 256]
    for case in record["cases"]:
        assert case["time"] == 0.3
        assert case["dt"] == 0.0005
        assert case["fft_applied"] is True
        assert case["initial_state_called"] is False
        assert case["W_rebuilt"] is False
        assert case["child_columns"] == [2, 3]
        assert len(case["pins"]["W_sha256"]) == 64
        assert len(case["pins"]["phi_sha256"]) == 64
    text = memory.closure_bytes(record)
    assert "nsc_spherical_galerkin_coupling.py" in text
    with pytest.raises(PermissionError):
        memory.assert_output_outside_repository("/Users/admin/Documents/BlackHoles-Infinity/lab/results/closure.json")


def test_full_segment_duration_and_cpu_forecast_admit_only_inside_the_budget():
    extent = memory.segment_extent(0.7, 0.0005)
    assert extent["steps"] == 1400
    assert extent["time0"] == 0.3
    assert extent["time1"] == 1.0
    assert extent["full_segment"] is True
    assert memory.segment_extent(0.7, 0.001)["steps"] == 700
    assert memory.segment_extent(0.7, 0.00025)["steps"] == 2800
    with pytest.raises(ValueError):
        memory.segment_extent(2.7, 0.0005)
    admitted = memory.cpu_forecast(1400, 0.02, cpu_budget_seconds=21600)
    refused = memory.cpu_forecast(1400, 0.02, cpu_budget_seconds=1.0)
    assert admitted["forecast_cpu_seconds"] == pytest.approx(42.0)
    assert admitted["admitted"] is True
    assert admitted["checkpoint_limit_bytes"] == 64 * 1024 * 1024
    assert refused["admitted"] is False


def test_write_and_check_full_segment_plan_without_evolving_it():
    directory = tempfile.mkdtemp(prefix="nsc-coupled-memory-plan-")
    manifest = memory.prepare_campaign(
        directory, fermions=(128,), duration=0.7, step_cap=0.0005, cpu_budget_seconds=21600,
    )
    assert manifest["full_segment"] is True
    assert manifest["duration"] == pytest.approx(0.7)
    assert manifest["time1"] == pytest.approx(1.0)
    assert manifest["steps"] == 1400
    assert manifest["production_evolved"] is False
    assert manifest["admitted"] is True
    assert manifest["forecast_cpu_seconds"] > 0.0
    assert manifest["checkpoint_limit_bytes"] == 64 * 1024 * 1024
    assert manifest["dephasing_prerequisite_for_coupled_response"] is False
    assert manifest["code_sha256"]["src/recursive_horizons/nsc_spherical_galerkin_coupling.py"]
    stopped = memory.run_campaign(directory, cpu_budget_seconds=0.0)
    assert stopped["production_evolved"] is False
    assert stopped["cases"][0]["status"] == "budget_stop"
    assert stopped["cases"][0]["production_steps_taken"] == 0
    report = memory.check_campaign(directory)
    assert report["stepped"] is False
    assert report["chunks"] == 1
    assert report["duration"] == pytest.approx(0.7)
    assert report["replay_force_gap"] == 0.0


def test_short_campaign_replays_force_and_geometry_without_stepping(monkeypatch):
    directory = tempfile.mkdtemp(prefix="nsc-coupled-memory-short-")
    memory.prepare_campaign(
        directory, fermions=(128,), duration=0.01, step_cap=0.0005, cpu_budget_seconds=21600,
    )
    manifest = memory.run_campaign(directory)
    assert manifest["cases"][0]["status"] == "completed"
    assert manifest["cases"][0]["production_steps_taken"] == 20
    assert manifest["cases"][0]["coordinate_time"] == pytest.approx(0.31)
    assert manifest["production_evolved"] is False

    def _refuse_step(*_args, **_kwargs):
        raise AssertionError("check took a coupled step")

    monkeypatch.setattr(memory.model, "rk4_step", _refuse_step)
    monkeypatch.setattr(memory, "rk4_split", _refuse_step)
    report = memory.check_campaign(directory)
    assert report["stepped"] is False
    assert report["chunks"] == 2
    _assert_below(report["replay_force_gap"], 1e-12, "replay_force")
    _assert_below(report["replay_clock_gap"], 1e-12, "replay_clock")
    _assert_below(report["replay_partial_gap"], 1e-8, "replay_partial")
