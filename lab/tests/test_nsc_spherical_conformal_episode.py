"""Independent frozen-control and saved finite forcing-budget checks."""
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

import derive_nsc_spherical_conformal_episode as episode
from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


def test_exact_frozen_control_matches_the_same_dense_dirac_group():
    count, length, kappa = 16, 8., 1
    phi0, phi1 = galerkin.manufactured_columns(count)
    grid = galerkin.build_grid(count, quadrature=64)
    initial = galerkin.blank_state(grid, phi0, phi1)
    p, _metric = coupling.antiperiodic_momentum(count, length)
    mass = kappa * initial.Q[0] * np.eye(count)
    zero = np.zeros_like(p)
    h = np.block([[zero, mass - 1j * p], [mass + 1j * p, zero]])
    expected = expm(-1j * .05 * h) @ np.vstack((phi0, phi1))
    actual0, actual1 = episode.frozen_columns(initial, length, kappa, .05)
    assert np.max(abs(np.vstack((actual0, actual1)) - expected)) < 1e-13


def test_saved_initial_fields_and_final_frames_are_bound_to_the_same_v5_data():
    record = json.loads(episode.OUT.read_text())
    assert record["schema"] == episode.SCHEMA
    assert record["gauge"] == "conformal"
    assert record["source_changed"] is False
    assert record["initial_reset"] is False
    assert record["bound_sources_unchanged_during_run"] is True
    assert record["source_hashes_before"]["v5_npz"] == episode.sha256(episode.legacy.V5_NPZ)
    assert record["source_hashes_before"]["legacy_episode_npz"] == episode.sha256(episode.legacy.NPZ)
    with np.load(episode.NPZ, allow_pickle=False) as saved, np.load(episode.legacy.V5_NPZ, allow_pickle=False) as v5:
        for nf in (256, 512):
            prefix = "" if nf == 256 else "nf512_"
            for name in episode.STATE_NAMES:
                original = "columns_" + name if name.startswith("phi") else "geometry_" + name
                assert np.array_equal(saved[f"nf{nf}_initial_{name}"], v5[prefix + original])
        for spec in episode.RUN_PLAN:
            for name in episode.STATE_NAMES:
                assert np.array_equal(saved[spec["name"] + "_final_" + name], saved[spec["name"] + "_frames_" + name][-1])
    assert record["cpu_seconds"] <= episode.CPU_BUDGET
    assert episode.OUT.stat().st_size + episode.NPZ.stat().st_size <= episode.PAYLOAD_LIMIT


def test_saved_forcing_integrals_and_effects_are_recomputed_from_every_sample():
    record = json.loads(episode.OUT.read_text())
    for spec in episode.RUN_PLAN:
        rows = record["series"][spec["name"]]
        reported = record["results"][spec["name"]]
        times = np.array([row["time"] for row in rows])
        forcing = np.array([row["forcing_norm"] for row in rows])
        dt = spec["dt"]
        trapezoid = np.sum((forcing[:-1] + forcing[1:]) * np.diff(times) / 2)
        simpson = dt / 3 * (forcing[0] + forcing[-1] + 4 * np.sum(forcing[1:-1:2]) + 2 * np.sum(forcing[2:-1:2]))
        assert abs(reported["integrals"]["forcing_norm"]["trapezoid"] - trapezoid) < 1e-14
        assert abs(reported["integrals"]["forcing_norm"]["simpson"] - simpson) < 1e-14
        expected_margin = rows[0]["constraint_norm"] + trapezoid - rows[-1]["constraint_norm"]
        assert abs(reported["constraint_end_margin"] - expected_margin) < 1e-14
        assert reported["constraint_end_margin"] > 0
        assert reported["positive_chart"] is True
        assert reported["continuum_constraint_certified"] is False
        assert reported["forcing_time_integral_certified"] is False
        assert reported["renewal"] is False
        assert rows[-1]["time"] == .05
        field_delta = rows[-1]["field_energy"] - rows[0]["field_energy"]
        gravity_delta = rows[-1]["gravity_energy"] - rows[0]["gravity_energy"]
        assert abs(reported["field_energy_change"] - field_delta) < 1e-15
        assert abs(reported["energy_balance_change"] - field_delta - gravity_delta) < 1e-15
        assert np.max(abs(np.asarray(reported["occupation_change"]) - np.asarray(rows[-1]["occupation"]) + np.asarray(rows[0]["occupation"]))) < 1e-15
        assert reported["sbp_identity_error_max"] < 1e-10
