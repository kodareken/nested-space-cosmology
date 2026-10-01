"""Physical interval integrals and the finalized dynamic-lapse frame ledger."""
import json

import numpy as np

import derive_nsc_spherical_conformal_episode as episode
import derive_nsc_spherical_conformal_frames as frames
from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons import nsc_regional_energy_exchange as regional


def test_interval_weights_integrate_the_retained_trigonometric_density_exactly():
    count, length, start, end = 128, 8., 0., 2.
    coordinate = np.arange(count) * length / count
    window = frames.interval_weights(count, length, start, end)
    for mode in (0, 1, 3, 17, 31):
        wave = mode * 2 * np.pi / length
        density = np.cos(wave * coordinate)
        expected = end - start if mode == 0 else (np.sin(wave * end) - np.sin(wave * start)) / wave
        assert abs(np.dot(window, density) * length / count - expected) < 1e-13
    partition = np.sum([frames.interval_weights(count, length, a, b) for a, b in zip(frames.BOUNDS[:-1], frames.BOUNDS[1:])], axis=0)
    assert np.max(abs(partition - 1)) < 1e-13


def test_dynamic_lapse_directional_normal_energy_matches_the_original_lapse_independent_density():
    grid = galerkin.build_grid(16, quadrature=64, gauge="conformal")
    phi0, phi1 = galerkin.manufactured_columns(16)
    state = galerkin.blank_state(grid, phi0, phi1)
    state.Q += .03 * np.cos(2 * np.pi * grid.xi_g / grid.length)
    state.r += .05 * np.sin(2 * np.pi * grid.xi_g / grid.length)
    state.p_chi += .1
    coarse, bundle = galerkin.compose_fine_hamiltonian(grid, state)
    fine, system = bundle["fine_state"], bundle["fine_system"]
    lifted = episode.legacy.lifted_rate(grid, coarse, bundle)
    actual = frames.stage_normal_slope(system, fine, lifted)
    held_lapse = regional.proper_pointwise_slope(system, fine, lifted)
    assert np.max(abs(actual - held_lapse)) < 1e-8
    source = coupling.source_from_columns(system, fine)
    dynamic = np.sum(source["force_Q"] * lifted.Q + source["force_L"] * lifted.Q)
    assert abs(dynamic - coarse.fieldwork_power) < 1e-12
    assert abs(dynamic - np.sum(source["force_Q"] * lifted.Q)) > 1e-3


def test_saved_profile_and_window_ledger_keep_internal_flux_separate_from_outer_boundary_flux():
    record = json.loads(frames.OUT.read_text())
    assert record["source_trajectory_unchanged"] is True
    assert record["source_trajectory_hashes"]["episode_json"] == episode.sha256(episode.OUT)
    assert record["source_trajectory_hashes"]["episode_npz"] == episode.sha256(episode.NPZ)
    assert record["combined_episode_and_frames_cpu"] < episode.CPU_BUDGET
    with np.load(frames.NPZ, allow_pickle=False) as payload:
        for spec in episode.RUN_PLAN:
            label = spec["name"]
            report = record["results"][label]
            normal = payload[label + "_normal_energy_nodal"]
            current = payload[label + "_normal_flux_nodal"]
            probability = payload[label + "_probability_nodal"]
            count, length = normal.shape[1], 8.
            window = frames.interval_weights(count, length, 0., 2.)
            observed_change = float(np.dot(window, normal[-1] - normal[0]))
            assert abs(observed_change - report["windows"][0]["normal_shell_change"]) < 1e-13
            assert abs(np.sum(probability[-1]) - report["final"]["canonical_probability"]) < 1e-13
            sampled_current_max = float(np.max(abs(current[:, :count // 4])) / (length / count))
            assert abs(sampled_current_max - report["windows"][0]["max_internal_normal_flux"]) < 1e-12
            shell = report["windows"][0]
            assert shell["max_internal_normal_flux"] > 1.
            assert shell["max_abs_normal_boundary_flux"] < 1e-7
            assert report["shell_leader_unchanged"] is True
            assert report["probability_leader_unchanged"] is True
            assert report["continuum_balance_certified"] is False
            assert abs(shell["normal_shell_integrated_gap"]) < .01 * abs(shell["normal_shell_change"])
    size = sum(path.stat().st_size for path in (episode.OUT, episode.NPZ, frames.OUT, frames.NPZ))
    assert size < episode.PAYLOAD_LIMIT
