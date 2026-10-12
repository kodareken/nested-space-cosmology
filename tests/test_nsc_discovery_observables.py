"""Stage-1 observer contract, manufactured window transport, and atlas check.

No saved trajectory is evolved and the production atlas file is not created.
"""
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

import derive_nsc_discovery_atlas as atlas
import derive_nsc_spherical_conformal_curvature as curvature
import derive_nsc_spherical_conformal_frames as frames
from recursive_horizons import nsc_discovery_observables as observer
from recursive_horizons import nsc_nested_parent_child as nested
from recursive_horizons import nsc_spherical_episode_assessment as metric
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


LAB = Path(__file__).resolve().parents[1]
PRODUCTION_ATLAS = LAB / "results/development/nsc-discovery-atlas-v1.json"


def manufactured(pair):
    base = galerkin.blank_state(pair.grid, pair.source_phi0, pair.source_phi1)
    coordinate = pair.grid.xi_g
    angle = 2 * np.pi * coordinate / pair.grid.length
    base.Q = 1.1 + .05 * np.cos(angle)
    base.r = 1.2 + .04 * np.sin(2 * angle)
    base.chi = .03 * np.cos(2 * angle)
    base.p_Q = .07 + .03 * np.sin(angle)
    base.p_r = .06 * np.cos(2 * angle)
    base.p_chi = .04 + .02 * np.sin(2 * angle)
    generator = np.random.default_rng(4701)
    columns, _ = np.linalg.qr(
        generator.normal(size=(2 * pair.grid.nf, 6)) + 1j * generator.normal(size=(2 * pair.grid.nf, 6))
    )
    base.phi0, base.phi1 = columns[:pair.grid.nf], columns[pair.grid.nf:]
    return nested.encode_state(pair, base)


def _evaluate(values, coordinate, length):
    samples = np.asarray(values, dtype=float)
    coefficients = np.fft.fft(samples) / samples.size
    modes = np.fft.fftfreq(samples.size) * samples.size
    phase = np.exp(2j * np.pi * float(coordinate) * modes / float(length))
    return float(np.real(phase @ coefficients))


def _keys(value, found):
    if isinstance(value, dict):
        found.update(value)
        for item in value.values():
            _keys(item, found)
    elif isinstance(value, list):
        for item in value:
            _keys(item, found)


@pytest.fixture(scope="module")
def prepared():
    pair = nested.build_pair(64)
    state = manufactured(pair)
    nodal = nested.reconstruct_state(pair, state)
    rate, bundle = galerkin.compose_fine_hamiltonian(pair.grid, nodal)
    full = observer.observe(pair, state, 0.2, capture_profiles=True, bundle=bundle, nodal_rate=rate)
    bare = observer.observe(pair, state, 0.2, bundle=bundle, nodal_rate=rate)
    return {
        "pair": pair,
        "state": state,
        "nodal": nodal,
        "rate": rate,
        "bundle": bundle,
        "full": full,
        "bare": bare,
    }


def test_observation_is_not_due_on_a_runge_kutta_stage():
    assert observer.OBSERVATION_ON_RK_STAGE is False
    assert observer.observation_due(0, every=5) is True
    assert observer.observation_due(5, every=5) is True
    assert observer.observation_due(4, every=5) is False
    assert observer.observation_due(5, every=5, on_rk_stage=True) is False


def test_compact_row_matches_the_frozen_contract(prepared):
    row = prepared["bare"]
    observer.assert_compact_row(row)
    assert "profiles" not in row
    names = set()
    _keys(row, names)
    assert "dispersal_time" not in names
    child = row["child_window"]
    parent = row["parent_window"]
    assert child["kind"] == "fixed_coordinate_window"
    assert row["mode_projectors"]["kind"] == "canonical_mode_projector"
    assert abs(child["source_probability"] - child["source_probability_per_proper_length"] * row["proper"]["child_proper_length"]) < 1e-9
    assert abs(parent["source_probability"] - parent["source_probability_per_proper_length"] * row["proper"]["parent_proper_length"]) < 1e-9
    assert child["source_probability"] <= parent["source_probability"] + 1e-10
    assert abs(row["source_probability_total"] - row["source_probability_interpolant"]) < 1e-8
    assert abs(row["source_probability_total"] - float(np.sum(prepared["pair"].weights))) < 1e-8
    catalog = {item["x"]: item for item in row["boundary_flux"]}
    assert abs(child["normal_boundary_flux"] - (catalog[1.0]["normal_energy_flux_density"] - catalog[3.0]["normal_energy_flux_density"])) < 1e-12
    assert abs(parent["normal_boundary_flux"] - (catalog[0.0]["normal_energy_flux_density"] - catalog[4.0]["normal_energy_flux_density"])) < 1e-12
    assert catalog[1.0]["window"] == "child" and catalog[0.0]["window"] == "parent"
    assert abs(child["flux_projection_shift"]) < 1e-12
    assert abs(child["flux_projection_shift_pressure"]) < 1e-12
    assert abs(child["source_probability"] - row["mode_projectors"]["source_child_column_power"]) > 1e-4
    for rate in row["proper"]["clock_rates"]:
        assert rate not in row["characteristics"]["crossing_times"].values()
    assert row["characteristics"]["ambiguities"] == []
    assert row["control_mode"] == "coupled"
    assert row["frozen_geometry_control"]["status"] == "unspecified"
    assert row["nyquist_symbol_max"] < 1e-8


def test_profiles_are_separate_from_the_scalar_row(prepared):
    full, bare = prepared["full"], prepared["bare"]
    assert set(full["profiles"]) == set(observer.PROFILE_KEYS)
    for key in observer.REQUIRED_TOP_KEYS:
        if key != "capture_profiles":
            assert full[key] == bare[key]
    peak = int(np.argmax(full["profiles"]["source_probability_per_proper_length"]))
    assert full["profiles"]["x"][peak] == bare["source_probability_per_proper_length_peak_x"]


def test_live_curvature_reuses_the_projected_jet_and_frozen_map(prepared):
    pair, nodal, rate, bundle = prepared["pair"], prepared["nodal"], prepared["rate"], prepared["bundle"]
    row = prepared["bare"]
    q_dot, q_ddot, before = curvature.actual_q_second_rate(pair.grid, nodal, rate, bundle)
    zeros = np.zeros(pair.grid.nq)
    direct = metric.direct_rh_grid(bundle["fine_state"].Q, bundle["fine_state"].Q, zeros, q_dot, q_ddot, pair.grid.length, L_dot=q_dot, beta_dot=zeros)
    fine = bundle["fine_state"]
    assert abs(row["curvature"]["R_h_max"] - float(np.max(direct["R_h"]))) < 1e-12
    assert abs(row["curvature"]["R_h_min"] - float(np.min(direct["R_h"]))) < 1e-12
    assert abs(row["curvature"]["chi_shell_gap_max"] - float(np.max(np.abs(direct["R_h"] - fine.chi - 2.0)))) < 1e-12
    assert row["curvature"]["auxiliary_substituted"] is False
    weyl = metric.weyl_actual(direct["R_h"], fine.r)
    assert abs(row["curvature"]["weyl_C2_max"] - float(np.max(weyl))) < 1e-12
    assert abs(row["curvature"]["Qddot_projection_gap_max"] - float(np.max(np.abs(q_ddot - before)))) < 1e-12
    assert np.max(np.abs(nodal.Q - pair.geometry_map @ prepared["state"].Q)) == 0.0
    assert pair.geometry_map.flags.writeable is False


def test_frozen_geometry_zeros_jets_and_keeps_spatial_lapse_work(prepared):
    pair, state, bundle, rate = prepared["pair"], prepared["state"], prepared["bundle"], prepared["rate"]
    coupled = prepared["bare"]
    frozen = observer.observe(
        pair, state, 0.2, bundle=bundle, nodal_rate=rate, control_mode="frozen_geometry",
    )
    named = observer.observe(
        pair, state, 0.2, bundle=bundle, nodal_rate=rate, preparation="separated_pair_v1",
    )
    original = observer.observe(
        pair, state, 0.2, bundle=bundle, nodal_rate=rate, preparation="original_conformal_v2",
    )
    assert coupled["curvature"]["Q_dot_max"] > 1e-6
    assert frozen["curvature"]["Q_dot_max"] == 0.0
    assert frozen["curvature"]["Q_ddot_max"] == 0.0
    assert frozen["coordinate_metric_work_total"] == 0.0
    assert frozen["child_window"]["proper_pressure_work"] == 0.0
    assert abs(frozen["child_window"]["momentum_lapse_work"] - coupled["child_window"]["momentum_lapse_work"]) < 1e-12
    assert abs(frozen["proper"]["child_r_proper_mean"] - coupled["proper"]["child_r_proper_mean"]) < 1e-12
    assert named["frozen_geometry_control"]["status"] == "saved_control"
    assert original["frozen_geometry_control"]["status"] == "missing"
    assert named["control_mode"] == "coupled"


def test_stale_bundle_and_state_mutation_are_rejected(prepared):
    pair, state, bundle, rate = prepared["pair"], prepared["state"], prepared["bundle"], prepared["rate"]
    before = state.Q.copy()
    observer.observe(pair, state, 0.2, bundle=bundle, nodal_rate=rate)
    assert np.array_equal(state.Q, before)
    shifted = state.copy()
    shifted.r = shifted.r * 1.2
    with pytest.raises(ValueError, match="frozen-W"):
        observer.observe(pair, shifted, 0.0, bundle=bundle, nodal_rate=rate)


def test_characteristic_images_are_ambiguous_and_are_not_a_dispersal_time():
    _tracks, at_two = observer.candidate_boundary_tracks(2.0)
    hits = {(item["track_boundary"], item["coordinate_speed"], item["coincides_with"]) for item in at_two}
    assert (1.0, 1.0, 3.0) in hits
    assert (3.0, -1.0, 1.0) in hits
    _later, at_four = observer.candidate_boundary_tracks(4.0)
    assert any(item["track_boundary"] == 0.0 and item["coincides_with"] == 4.0 for item in at_four)
    width = 0.25
    assert width not in observer.CROSSING_TIMES.values()


def test_constant_window_drops_reynolds_and_contact():
    grid = galerkin.build_grid(32, quadrature=64, gauge="conformal")
    coordinate = grid.xi_q
    length = grid.length
    density = 1.3 + 0.2 * np.cos(2 * np.pi * coordinate / length)
    flux = 0.1 * np.sin(2 * np.pi * coordinate / length)
    source = 0.05 * np.cos(4 * np.pi * coordinate / length)
    contact = 0.4 + 0.1 * np.sin(2 * np.pi * coordinate / length)
    interval = (1.25, 3.5)
    report = observer.moving_window_transport(grid, density, flux, source, interval, (0.0, 0.0))
    expected = nested.interval_integral(grid, source, interval) + _evaluate(flux, interval[0], length) - _evaluate(flux, interval[1], length)
    assert report["fixed_window"] is True
    assert report["reynolds"] == 0.0
    assert report["surface_work"] == 0.0
    assert abs(report["leibniz_rate"] - expected) < 1e-10


def test_moving_window_reynolds_matches_the_seed_manufactured_profile():
    grid = galerkin.build_grid(32, quadrature=64, gauge="conformal")
    coordinate = grid.xi_q
    length = grid.length
    density = 1.3 + 0.2 * np.cos(2 * np.pi * coordinate / length) + 0.05 * np.sin(4 * np.pi * coordinate / length)
    flux = np.zeros_like(density)
    source = np.zeros_like(density)
    interval = (1.25, 3.5)
    velocity = (0.15, -0.25)
    report = observer.moving_window_transport(grid, density, flux, source, interval, velocity)
    expected = _evaluate(density, interval[1], length) * velocity[1] - _evaluate(density, interval[0], length) * velocity[0]

    def shifted(epsilon):
        return nested.interval_integral(
            grid, density, (interval[0] + epsilon * velocity[0], interval[1] + epsilon * velocity[1]),
        )

    finite_difference = (shifted(1e-5) - shifted(-1e-5)) / 2e-5
    assert abs(report["reynolds"] - expected) < 1e-10
    assert abs(finite_difference - expected) < 1e-8
    assert abs(report["leibniz_rate"] - expected) < 1e-10
    assert report["surface_work"] == 0.0


def test_surface_action_is_explicit_and_artificial_cuts_keep_only_reynolds(prepared):
    grid = galerkin.build_grid(32, quadrature=64, gauge="conformal")
    coordinate = grid.xi_q
    length = grid.length
    eta = np.ones_like(coordinate)
    flux = 0.05 * np.cos(2 * np.pi * coordinate / length)
    interior = 0.02 * np.cos(4 * np.pi * coordinate / length)
    traction = 0.4 + 0.1 * np.sin(2 * np.pi * coordinate / length)
    interval = (1.25, 3.5)
    velocity = (0.2, -0.15)
    report = observer.moving_window_transport(grid, eta, flux, interior, interval, velocity)
    surface = observer.surface_action_work(grid, traction, interval, velocity)
    expected_surface = _evaluate(traction, interval[1], length) * velocity[1] - _evaluate(traction, interval[0], length) * velocity[0]
    assert abs(expected_surface) > 1e-4
    assert abs(surface - expected_surface) < 1e-10
    assert report["surface_work"] == 0.0
    assert abs(report["leibniz_rate"] - (report["interior_source"] + report["boundary_flux"] + report["reynolds"])) < 1e-12
    assert abs(report["leibniz_rate"] - (report["leibniz_rate"] + surface)) > 1e-4
    assert prepared["bare"]["child_window"]["surface_work"] == 0.0
    assert prepared["bare"]["child_window"]["energy_reynolds"] == 0.0


def test_atlas_check_distinguishes_preparations_and_does_not_write():
    existed = PRODUCTION_ATLAS.exists()
    report = atlas.check_inputs()
    assert existed is False and PRODUCTION_ATLAS.exists() is False
    assert report["writer_executed"] is False
    assert report["trajectories_loaded"] is False
    assert report["observe_called"] is False
    assert report["metadata_distinct"] is True
    named = {row["name"]: row for row in report["preparations"]}
    assert set(named) == set(observer.PREPARATION_NAMES)
    assert len({row["schema"] for row in report["preparations"] if row["name"] != "new_pair"}) == 4
    for name in ("original_conformal_v2", "new_pair"):
        assert named[name]["frozen_geometry_comparison"]["status"] == "missing"
        assert named[name]["frozen_geometry_comparison"]["saved_control_applied"] is False
    for name in ("separated_pair_v1", "separated_pair_v2", "frozen_basis"):
        comparison = named[name]["frozen_geometry_comparison"]
        assert comparison["status"] == "saved_control"
        assert comparison["hash_matches_saved_control"] is True
        assert comparison["svd_reconstructed"] is False
        assert comparison["saved_validation_svd_reconstructed"] is False
    assert named["separated_pair_v1"]["source_layout"] == "separated"
    assert named["original_conformal_v2"]["gauge"] == "conformal"
    assert named["separated_pair_v1"]["fingerprint"] != named["separated_pair_v2"]["fingerprint"]


def test_saved_nf256_frame_matches_the_reference_window_and_native_momenta():
    report = atlas.saved_frame_reference()
    assert report["time"] == pytest.approx(0.3)
    for channel in ("normal_energy", "coordinate_metric_work", "flux_shift", "flux_proper", "flux_shift_pressure"):
        numbers = report["channels"][channel]
        assert numbers["reported"] == pytest.approx(numbers["raw_sum"], rel=0, abs=1e-8 * max(1.0, abs(numbers["raw_sum"])))
        assert numbers["reference_window"] == pytest.approx(numbers["reported"], rel=0, abs=1e-8 * max(1.0, abs(numbers["raw_sum"])))
    weights = frames.interval_weights(report["count"], 8.0, 0.0, 8.0)
    direct = float(np.dot(weights, report["normal_energy_summand"]))
    assert direct == pytest.approx(report["channels"]["normal_energy"]["reported"], rel=0, abs=1e-8)
    assert report["momentum_conversion_gap"] < 1e-10
    assert report["frozen_Q_dot_max"] == 0.0
    assert report["frozen_coordinate_metric_work_total"] == 0.0
    assert report["frozen_child_lapse_work"] == pytest.approx(report["coupled_child_lapse_work"], rel=0, abs=1e-9)
    assert report["coupled_child_r_mean"] == pytest.approx(report["frozen_child_r_mean"], rel=0, abs=1e-12)


def test_atlas_writer_creates_only_a_new_physical_successor(tmp_path):
    destination = tmp_path / "nsc-discovery-atlas-v1.json"
    case = next(row for row in atlas.BASELINE_CASES if row["case"] == "nf256_dt_0.0005")
    record = atlas.write_atlas(destination, cases=(case,), indices=(-1,))
    assert destination.is_file()
    assert PRODUCTION_ATLAS.exists() is False
    assert record["schema"] == atlas.SCHEMA
    assert record["creation_only"] is True
    assert record["evolution_performed"] is False
    assert record["initial_solve_performed"] is False
    assert record["trajectories_created"] is False
    series = record["series"][0]
    assert series["preparation"] == "original_conformal_v2"
    assert series["rows"][0]["time"] == pytest.approx(0.3)
    assert np.isfinite(series["rows"][0]["R_h_max"])
    assert np.isfinite(series["rows"][0]["localization_peak_x"])
    assert series["final_frozen_geometry"]["Q_dot_max"] == 0.0
    evidence = record["evidence"]
    assert len(evidence["producing_commit"]) == 40
    assert set(evidence["dependencies"]) == {"producer", "observer", "backend", "core", "consumer"}
    assert tuple(evidence["inputs"]) == atlas.bound_input_paths()
    assert evidence["sealed_inputs_regenerated"] is False
    assert evidence["dependencies"]["observer"]["src/recursive_horizons/nsc_discovery_observables.py"] == atlas.sha256(
        LAB / "src/recursive_horizons/nsc_discovery_observables.py"
    )
    payload = destination.with_suffix(".npz")
    stamps = (destination.stat().st_mtime_ns, payload.stat().st_mtime_ns, destination.stat().st_size, payload.stat().st_size)
    with np.load(payload) as stored:
        assert stored["original_conformal_v2_nf256_dt_0.0005_final_R_h"].shape[0] > 0
        assert np.isfinite(stored["original_conformal_v2_nf256_dt_0.0005_normal_energy_total"]).all()
    completed = subprocess.run(
        [sys.executable, str(LAB / "scripts/derive_nsc_discovery_atlas.py"), "--check", "--successor", str(destination)],
        cwd=LAB, env=os_environ(), capture_output=True, text=True, check=False,
    )
    assert completed.returncode == 0, completed.stderr
    checked = json.loads(completed.stdout)
    assert checked["writer_executed"] is False
    assert checked["successor"]["authenticated"] is True
    assert checked["successor"]["rewritten"] is False
    assert (destination.stat().st_mtime_ns, payload.stat().st_mtime_ns, destination.stat().st_size, payload.stat().st_size) == stamps
    assert [row["case"] for row in atlas.BASELINE_CASES] == [
        "nf256_dt_0.0005", "nf512_dt_0.0005", "nf256_baseline_dt0.0005", "nf512_baseline_dt0.0005",
    ]
    with pytest.raises(FileExistsError):
        atlas.write_atlas(destination)
    sealed = LAB / "results/development/nsc-nested-parent-child-v1.json"
    with pytest.raises(FileExistsError):
        atlas.write_atlas(sealed)


def test_atlas_cli_check_is_readonly():
    env = os_environ()
    completed = subprocess.run(
        [sys.executable, str(LAB / "scripts/derive_nsc_discovery_atlas.py"), "--check"],
        cwd=LAB, env=env, capture_output=True, text=True, check=False,
    )
    assert completed.returncode == 0, completed.stderr
    report = json.loads(completed.stdout)
    assert report["writer_executed"] is False
    assert report["saved_frame"]["time"] == pytest.approx(0.3)
    assert report["saved_frame"]["channels"]["normal_energy"]["reported"] == pytest.approx(
        report["saved_frame"]["channels"]["normal_energy"]["raw_sum"], rel=0, abs=1e-6,
    )
    assert report["baseline_shapes"]["original_conformal_v2:nf512_dt_0.0005"]["final_time"] == pytest.approx(0.3)
    assert report["baseline_shapes"]["separated_pair_v2:nf512_baseline_dt0.0005"]["frames"] == 31
    assert PRODUCTION_ATLAS.exists() is False


def os_environ():
    import os
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([str(LAB / "src"), str(LAB / "scripts")])
    return env
