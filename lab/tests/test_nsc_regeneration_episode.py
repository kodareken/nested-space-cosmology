"""Read-only checks for the regeneration-episode continuation.

These tests do not evolve a new trajectory and do not write predecessor records.
"""
import json
import os
import copy

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import numpy as np
import pytest

import derive_nsc_regeneration_episode as episode_run
import derive_nsc_spherical_feedback_episode as episode
from recursive_horizons.nsc_regeneration_controls import load_episode_final, state_sha256


def test_timestep_caps_and_safety_adapter():
    coarse, mark, limiter = episode_run.choose_dt(0.05, 0.15, 100.0, 5e-4)
    assert abs(coarse - 5e-4) < 1e-12
    assert abs(mark - 0.055) < 1e-12
    assert "owned_dt_cap" in limiter
    confirm, _mark, confirm_limiter = episode_run.choose_dt(0.05, 0.15, 100.0, 2.5e-4)
    assert abs(confirm - 2.5e-4) < 1e-12
    assert confirm_limiter == "run_cap"
    adapted, _mark, safety = episode_run.choose_dt(0.05, 0.15, 10000.0, 5e-4)
    assert abs(adapted - 1.4e-4) < 1e-12
    assert safety == "safety_cap"
    assert adapted * 10000.0 <= 1.4 + 1e-12


def test_saved_arc_coordinates_stay_in_their_regions():
    packet = {"x_start": 1.4296875, "x_end": 1.6015625, "wraps": False}
    bridge = {"x_start": 5.09765625, "x_end": 5.8046875, "wraps": False}
    assert episode_run.interval_overlaps(packet, 0.0, 4.0)
    assert not episode_run.interval_overlaps(packet, 4.0, 8.0)
    assert episode_run.interval_overlaps(bridge, 4.0, 8.0)
    assert not episode_run.interval_overlaps(bridge, 0.0, 4.0)


def _quiet_flags():
    return {
        "anti_in_packet": False,
        "trapped_in_bridge": False,
        "large_same_sign_arcs": 2,
    }


def test_declared_events_ignore_the_handoff_and_the_share_proxy():
    start = {"flags": _quiet_flags(), "q_at_areal_max": 1.0}
    assert episode_run.geometric_event([start]) is None
    held = {"flags": _quiet_flags(), "q_at_areal_max": 1.2}
    assert episode_run.geometric_event([start, held]) is None
    entered = {"flags": {**_quiet_flags(), "anti_in_packet": True}, "q_at_areal_max": 1.3}
    assert episode_run.geometric_event([start, held, entered]) == "BRIDGE_ARC_ENTERED_PACKET"
    turned = [
        {"flags": _quiet_flags(), "q_at_areal_max": value}
        for value in (1.0, 1.2, 1.4, 1.3, 1.2)
    ]
    assert episode_run.geometric_event(turned) == "Q_TURNBACK"
    monotone = [
        {"flags": _quiet_flags(), "q_at_areal_max": value}
        for value in (1.0, 1.1, 1.2, 1.3)
    ]
    assert episode_run.geometric_event(monotone) is None


def test_window_identity_is_not_the_flux_sum():
    times = np.array([0.0, 1.0])
    values = np.array([[1.0, 2.0], [3.0, 6.0]])
    assert np.allclose(episode_run.trap_integrate(times, values), [2.0, 4.0])


def test_driver_does_not_call_sealed_writers():
    text = episode_run.DRIVER.read_text()
    for token in ("reassess_draft(", "write_record(", "run_diagnostic("):
        assert token not in text


def test_saved_record_replays_without_writing():
    if not episode_run.OUT.is_file():
        raise AssertionError("episode record was not written")
    before = {name: episode_run.sha256_file(path) for name, path in episode_run.SEALED_FILES.items()}
    owned = {path: episode_run.sha256_file(path) for path in episode_run.OWNED_OUTPUTS}
    record = episode_run.verify_saved(replay=True, source_ref=episode_run.SEALED_SOURCE_REF)
    after = {name: episode_run.sha256_file(path) for name, path in episode_run.SEALED_FILES.items()}
    assert before == after
    assert owned == {path: episode_run.sha256_file(path) for path in episode_run.OWNED_OUTPUTS}
    assert record["state_reset"] is False
    assert record["new_force_added"] is False
    assert record["old_toy_B_matched"] is False
    assert record["proxy_reversal_required"] is False
    assert record["proxy_share_required"] is False
    assert record["safety_cap_is_physical_failure_label"] is False
    assert record["payload_within_64MiB"] is True
    pilot = record["runs"]["nf256_dtmax_0_0005"]
    saved = load_episode_final(episode.NPZ, "nf256_dt_0_0005")
    assert pilot["initial_hash"] == state_sha256(saved)
    assert pilot["source_reset"] is False
    assert pilot["dt_max_realized"] <= 5e-4 + 1e-12
    assert pilot["dt_omega_max"] <= episode_run.ABSOLUTE_RK4 + 1e-8
    for name, balance in record["balances"].items():
        assert balance["flux_sum_is_not_the_window_balance"] is True
        assert "window_normal_residual" in balance
        assert "proper_pressure_work_integral" in balance
        assert "lapse_gradient_work_integral" in balance
        fresh = episode.balance_against_exchange(
            balance["total_energy_change"], balance["coordinate_work_integral"]
        )
        assert fresh["within_one_percent_of_exchange"] == (
            balance["balance_versus_coordinate_work"]["within_one_percent_of_exchange"]
        )
    if record["decision"]["open_remaining_cases"]:
        assert "nf512_dtmax_0_00025" in record["runs"]
        for row in record["effects"]:
            if row["status"] == "incomplete" or row.get("effect_scale") is None:
                continue
            assert row["status"] == episode.classify_movement(
                row["effect_scale"], row["movement"], row["floor"]
            )
    assert record["goal"]["radial_motion_alone_is_sufficient"] is False
    schema = record["payload_schema"]
    assert schema["metric_curvature_computed_here"] is False
    assert schema["chi_plus_2_is_not_the_only_curvature_input"] is True
    with np.load(episode_run.NPZ, allow_pickle=False) as payload:
        assert "nf256_dtmax_0_0005_frame_Q_dot" in payload
        assert "nf256_dtmax_0_0005_frame_p_chi_dot" in payload
        assert "nf256_dtmax_0_0005_frame_indicator_Q_dot" in payload
        assert "nf256_dtmax_0_0005_frame_increment_present" in payload
        present = np.asarray(payload["nf256_dtmax_0_0005_frame_increment_present"])
        assert bool(present[0]) is True
        phi0 = np.asarray(payload["nf256_dtmax_0_0005_final_phi0"])
        assert phi0.shape == saved.phi0.shape
        assert not np.array_equal(phi0, saved.phi0)
    note = episode_run.NOTE.read_text()
    assert "flux sum" in note.lower() or "flux sum" in note
    assert "indicator" in note.lower()


def test_sealed_source_replay_is_explicit_and_rejects_mutations():
    record = json.loads(episode_run.OUT.read_text())
    assert "galerkin" in episode_run._sealed_binding_errors(record)
    assert episode_run._sealed_binding_errors(record, source_ref=episode_run.SEALED_SOURCE_REF) == []
    changed = copy.deepcopy(record)
    changed["sealed_hashes_before"]["galerkin"] = "0" * 64
    assert "galerkin" in episode_run._sealed_binding_errors(changed, source_ref=episode_run.SEALED_SOURCE_REF)
    missing_hash = copy.deepcopy(record)
    del missing_hash["sealed_hashes_before"]["galerkin"]
    assert "galerkin" in episode_run._sealed_binding_errors(missing_hash, source_ref=episode_run.SEALED_SOURCE_REF)
    missing_hash["sealed_hashes_before"]["galerkin"] = None
    assert "galerkin" in episode_run._sealed_binding_errors(missing_hash, source_ref=episode_run.SEALED_SOURCE_REF)
    numerical = copy.deepcopy(record)
    numerical["sealed_hashes_before"]["episode_npz"] = "0" * 64
    assert "episode_npz" in episode_run._sealed_binding_errors(numerical, source_ref=episode_run.SEALED_SOURCE_REF)
    absent_ref = episode_run._sealed_binding_errors(record, source_ref="0" * 40)
    assert "galerkin" in absent_ref
    with pytest.raises(AssertionError, match="galerkin"):
        episode_run.verify_saved(replay=False)
