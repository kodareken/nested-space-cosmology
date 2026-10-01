"""Independent successor handoff and surface-flow interpretation."""
import json
import numpy as np

import derive_nsc_spherical_conformal_episode as episode
import derive_nsc_spherical_conformal_continuation as continuation
import derive_nsc_spherical_conformal_transport as transport
from recursive_horizons.nsc_spherical_coupling import CauchyState


def test_equal_nonzero_entry_and_exit_are_resolved_even_when_net_flow_cancels():
    physical = {}
    for spec in episode.RUN_PLAN:
        physical[spec["name"]] = [{"time": t, "windows": [
            {"normal_left_outward_flux": -1., "normal_right_outward_flux": 1., "normal_boundary_flux": 0.}
            for _ in range(4)]} for t in (0., .005, .01)]
    assessed = transport.assess_surfaces(physical)
    assert all(row["resolved_one_percent_indicators"] for row in assessed if row["channel"] != "normal_boundary_flux")
    assert all(not row["resolved_one_percent_indicators"] for row in assessed if row["channel"] == "normal_boundary_flux")
    assert all(abs(row["absolute_integral_effect"] - .01) < 1e-15 for row in assessed if row["channel"] != "normal_boundary_flux")


def test_saved_successor_handoff_and_original_observer_are_bitwise_unchanged():
    record = json.loads(continuation.OUT.read_text())
    assert record["gauge"] == "conformal"
    assert record["midtrajectory_reset"] is False
    assert record["source_changed"] is False
    assert record["sealed_sources_unchanged_during_run"] is True
    assert record["predecessor_hashes"]["v1_npz"] == episode.sha256(episode.NPZ)
    assert record["predecessor_hashes"]["v1_driver"] == episode.sha256(episode.DRIVER)
    with np.load(episode.NPZ, allow_pickle=False) as predecessor, np.load(continuation.NPZ, allow_pickle=False) as successor:
        for spec in episode.RUN_PLAN:
            label = spec["name"]
            for name in episode.STATE_NAMES:
                assert np.array_equal(successor[label + "_frames_" + name][0], predecessor[label + "_final_" + name])
                assert np.array_equal(successor[label + "_frames_" + name][-1], successor[label + "_final_" + name])
                assert np.array_equal(successor[f"nf{spec['nf']}_original_T0_{name}"], predecessor[f"nf{spec['nf']}_initial_{name}"])
            assert record["handoffs"][label]["v1_final_state_sha256"] == record["handoffs"][label]["v2_initial_state_sha256"]
            assert successor[label + "_frame_times"][0] == .05
            assert abs(successor[label + "_frame_times"][-1] - .3) < 1e-12
            assert successor[label + "_proper_clock_frames"][0] > 0
            assert successor[label + "_frames_L"].shape[1] == 4 * spec["nf"]
            assert np.array_equal(successor[label + "_frames_L_dot"], successor[label + "_frames_rate_Q"])
    for pool in record["cpu_pools"].values():
        assert pool["actual_cpu_seconds"] <= pool["budget_seconds"]
    assert record["combined_v1_v2_payload_bytes"] < episode.PAYLOAD_LIMIT


def test_successor_constraint_budget_and_admissibility_use_actual_step_times():
    record = json.loads(continuation.OUT.read_text())
    for spec in episode.RUN_PLAN:
        label = spec["name"]
        series = record["series"][label]
        result = record["results"][label]
        times = np.array([row["time"] for row in series])
        force = np.array([row["forcing_norm"] for row in series])
        integral = float(np.dot(np.diff(times), .5 * (force[:-1] + force[1:])))
        assert abs(integral - result["integrals"]["forcing_norm"]["trapezoid"]) < 1e-14
        assert abs(result["constraint_end_margin"] - (series[0]["constraint_norm"] + integral - series[-1]["constraint_norm"])) < 1e-14
        assert result["constraint_end_margin"] > 0
        assert result["positive_chart"] is True
        assert result["gram_gap_max"] < continuation.GRAM_LIMIT
        assert result["target_reached"] is True
        assert result["continuum_constraint_certified"] is False
        for dt, row in zip(np.diff(times), series[1:]):
            assert dt * row["omega"] < 2 * np.sqrt(2)
            assert dt <= spec["dt"] + 1e-12


def test_surface_record_preserves_signed_and_absolute_indicators_without_a_net_veto():
    record = json.loads(transport.OUT.read_text())
    assert record["source_unchanged"] is True
    assert record["net_zero_is_a_veto"] is False
    assert record["source_bindings"]["episode_v2"] == episode.sha256(continuation.OUT)
    assert record["source_bindings"]["payload_v2"] == episode.sha256(continuation.NPZ)
    assert len(record["surface_assessments"]) == 12
    for row in record["surface_assessments"]:
        if row["resolved_one_percent_indicators"]:
            assert max(row["time_fraction_of_effect"], row["space_fraction_of_effect"], row["frame_fraction_of_effect"]) < .01
        assert row["certified"] is False
    assert record["observable_error_certified"] is False
