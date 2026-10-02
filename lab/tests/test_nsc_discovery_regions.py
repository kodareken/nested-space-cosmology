"""Manufactured stage-4 region checks. No production trajectory is loaded."""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

import derive_nsc_discovery_regions as cli
import derive_nsc_nested_parent_child as owned
from recursive_horizons import nsc_discovery_episode as episode
from recursive_horizons import nsc_discovery_observables as observer
from recursive_horizons import nsc_discovery_regions as regions
from recursive_horizons import nsc_nested_parent_child as nested
from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


LAB = Path(__file__).resolve().parents[1]
ROOT = LAB.parent
PYTHON = Path("/Users/admin/Documents/BlackHoles-Infinity/.venv/validation/bin/python")
SEALED = tuple(LAB / relative for relative in regions.SEALED_INPUTS)


def _sealed_stat():
    return tuple((path.stat().st_mtime_ns, path.stat().st_size, path.stat().st_mode) for path in SEALED)


def _state(pair, *, source_columns):
    base = galerkin.blank_state(pair.grid, pair.source_phi0, pair.source_phi1)
    coordinate = pair.grid.xi_g
    angle = 2.0 * np.pi * coordinate / pair.grid.length
    base.Q = 1.1 + 0.05 * np.cos(angle)
    base.r = 1.2 + 0.04 * np.sin(2.0 * angle)
    base.chi = 0.03 * np.cos(2.0 * angle)
    base.p_Q = 0.07 + 0.03 * np.sin(angle)
    base.p_r = 0.06 * np.cos(2.0 * angle)
    base.p_chi = 0.04 + 0.02 * np.sin(2.0 * angle)
    if source_columns:
        base.phi0 = np.array(pair.source_phi0, copy=True)
        base.phi1 = np.array(pair.source_phi1, copy=True)
    else:
        generator = np.random.default_rng(4701)
        columns, _ = np.linalg.qr(
            generator.normal(size=(2 * pair.grid.nf, 6)) + 1j * generator.normal(size=(2 * pair.grid.nf, 6))
        )
        base.phi0, base.phi1 = columns[:pair.grid.nf], columns[pair.grid.nf:]
    return nested.encode_state(pair, base)


def _prepare(pair, state):
    nodal = nested.reconstruct_state(pair, state)
    rate, bundle = galerkin.compose_fine_hamiltonian(pair.grid, nodal)
    return nodal, rate, bundle


@pytest.fixture(scope="module")
def prepared():
    pair = nested.build_pair(64)
    source_state = _state(pair, source_columns=True)
    live_state = _state(pair, source_columns=False)
    source = _prepare(pair, source_state)
    live = _prepare(pair, live_state)
    return {
        "pair": pair,
        "source_state": source_state,
        "live_state": live_state,
        "source_bundle": source,
        "live_bundle": live,
    }


def _translation(grid, speed):
    length = float(grid.length)
    coordinate = np.asarray(grid.xi_q, dtype=float)
    wave = 2.0 * math.pi / length
    offset, amplitude = 1.2, 0.8
    density = offset + amplitude * np.cos(wave * (coordinate - 0.5 * length))
    density_rate = -speed * np.asarray(grid.derivative @ density, dtype=float)
    spacing = float(grid.dx_q)
    return {
        "offset": offset,
        "amplitude": amplitude,
        "wave": wave,
        "density": density,
        "density_rate": density_rate,
        "energy": density * spacing,
        "flux": speed * density * spacing,
        "slope": density_rate * spacing,
        "zeros": np.zeros(grid.nq),
    }


def _analytic_rate(profile, left, right, left_dot, right_dot, speed):
    wave = profile["wave"]
    offset = profile["offset"]
    amplitude = profile["amplitude"]
    center = 0.5 * float(2.0 * math.pi / wave)

    def value(coordinate):
        return offset + amplitude * math.cos(wave * (coordinate - center))

    def antiderivative(coordinate):
        # ∫ v A k sin(k(x-c)) dx = -v A cos(k(x-c))
        return -speed * amplitude * math.cos(wave * (coordinate - center))

    accumulated = antiderivative(right) - antiderivative(left)
    return accumulated + value(right) * right_dot - value(left) * left_dot


def test_interpreter_is_the_pinned_validation_environment():
    assert Path(sys.executable) == PYTHON


def test_moving_cut_matches_translating_density_analytic_leibniz():
    grid = galerkin.build_grid(32, gauge="conformal")
    speed = 0.37
    profile = _translation(grid, speed)
    left, right = 1.2, 3.4
    left_dot, right_dot = 0.2, -0.15
    rate = regions.moving_normal_energy_rate(
        grid, profile["energy"], profile["flux"], profile["zeros"], profile["zeros"], profile["slope"],
        (left, right), (left_dot, right_dot),
    )
    expected = _analytic_rate(profile, left, right, left_dot, right_dot, speed)
    assert rate["piston_work_included"] is False
    assert rate["surface_work"] == 0.0
    assert rate["physical_wall"] is False
    assert rate["dx_divisions"] == 1
    assert abs(rate["closure_residual"]) < 1.0e-9
    assert rate["rate"] == pytest.approx(expected, abs=1.0e-9)
    assert rate["rate"] == pytest.approx(
        rate["reynolds"] + rate["fixed_boundary_flux"] + rate["pressure"] + rate["lapse"] + rate["finite_projection_defect"]
    )
    fixed = regions.moving_normal_energy_rate(
        grid, profile["energy"], profile["flux"], profile["zeros"], profile["zeros"], profile["slope"],
        (left, right), (0.0, 0.0),
    )
    assert fixed["reynolds"] == 0.0
    assert fixed["rate"] == pytest.approx(_analytic_rate(profile, left, right, 0.0, 0.0, speed), abs=1.0e-9)
    comoving = regions.moving_normal_energy_rate(
        grid, profile["energy"], profile["flux"], profile["zeros"], profile["zeros"], profile["slope"],
        (left, right), (speed, speed),
    )
    assert comoving["rate"] == pytest.approx(0.0, abs=1.0e-9)
    peak_at_branch = profile["offset"] + profile["amplitude"] * np.cos(profile["wave"] * grid.xi_q)
    branch_rate = -speed * np.asarray(grid.derivative @ peak_at_branch, dtype=float)
    spacing = float(grid.dx_q)
    branch = regions.density_contours(grid, peak_at_branch, branch_rate)
    half = next(level for level in branch["levels"] if level["fraction"] == 0.5)
    assert half["wrapped"] is True
    wrapped = regions.moving_normal_energy_on_cut(
        grid, peak_at_branch * spacing, speed * peak_at_branch * spacing, profile["zeros"], profile["zeros"],
        branch_rate * spacing, half["backward"], half["forward"], speed, speed, wrapped=True,
    )
    assert wrapped["physical_wall"] is False
    assert wrapped["piston_work_included"] is False
    assert wrapped["rate"] == pytest.approx(0.0, abs=1.0e-9)
    assert abs(wrapped["closure_residual"]) < 1.0e-9


def test_fixed_window_normal_energy_matches_owned_slope(prepared):
    pair = prepared["pair"]
    state = prepared["live_state"]
    _nodal, rate, bundle = prepared["live_bundle"]
    fine = bundle["fine_state"]
    system = bundle["fine_system"]
    lifted = observer._lifted_rate(pair.grid, rate, bundle)
    owned_slope = owned.normal_energy_slope(system, fine, lifted, bundle["source"])
    slope = regions.analytic_normal_energy_slope(system, fine, lifted, bundle["source"])
    assert np.max(np.abs(slope - owned_slope)) < 1.0e-12
    row = regions.analyze(pair, state, 0.2, bundle=bundle, nodal_rate=rate, control_mode="coupled")
    for name, interval in (("child_fixed", nested.CHILD_INTERVAL), ("parent_fixed", nested.PARENT_INTERVAL)):
        window = row["normal_energy"][name]
        integral = nested.interval_integral(pair.grid, owned_slope / system.dx, interval)
        assert window["boundary_velocity"] == [0.0, 0.0]
        assert window["reynolds"] == 0.0
        assert window["piston_work_included"] is False
        assert window["rate"] == pytest.approx(integral, abs=1.0e-8)
        assert abs(window["closure_residual"]) < 1.0e-8
    assert row["normal_energy"]["analytical_continuity_closure_residual_max"] < 1.0e-8
    assert row["renewal_asserted"] is False
    assert row["contour_certified"] is False
    def per_proper(sample):
        weights = np.asarray(system.occupations, dtype=float)[None, :]
        mass = np.sum((np.abs(sample.phi0) ** 2 + np.abs(sample.phi1) ** 2) * weights, axis=1)
        return (mass / system.dx) / (sample.r * sample.Q)

    step = 1.0e-6
    finite_difference = (
        per_proper(coupling._combine(fine, lifted, step)) - per_proper(coupling._combine(fine, lifted, -step))
    ) / (2.0 * step)
    actual = regions._probability_rate(fine, system, lifted)
    scale = max(float(np.max(np.abs(actual))), 1.0e-12)
    assert float(np.max(np.abs(finite_difference - actual))) / scale < 1.0e-6


def test_frozen_geometry_zeros_coordinate_rates_and_keeps_the_field_source(prepared):
    pair = prepared["pair"]
    state = prepared["live_state"]
    _nodal, rate, bundle = prepared["live_bundle"]
    coupled = regions.analyze(pair, state, 0.2, bundle=bundle, nodal_rate=rate, control_mode="coupled")
    frozen = regions.analyze(pair, state, 0.2, bundle=bundle, nodal_rate=rate, control_mode="frozen_geometry")
    assert coupled["rates"]["Q_dot_max"] > 1.0e-8
    assert coupled["rates"]["r_dot_max"] > 1.0e-8
    assert frozen["rates"]["Q_dot_max"] == 0.0
    assert frozen["rates"]["r_dot_max"] == 0.0
    assert frozen["rates"]["geometry_jets_zero"] is True
    assert frozen["rates"]["field_column_rate_max"] > 0.0
    assert frozen["rates"]["field_source_kernel_rate_max"] > 0.0
    assert frozen["normal_energy"]["child_fixed"]["pressure"] == pytest.approx(0.0, abs=1.0e-12)
    assert abs(frozen["normal_energy"]["child_fixed"]["lapse"] - coupled["normal_energy"]["child_fixed"]["lapse"]) < 1.0e-12
    assert abs(frozen["normal_energy"]["child_fixed"]["rate"] - coupled["normal_energy"]["child_fixed"]["rate"]) > 1.0e-8
    observed = observer.observe(
        pair, state, 0.2, bundle=bundle, nodal_rate=rate, control_mode="frozen_geometry",
    )
    assert observed["curvature"]["Q_dot_max"] == 0.0
    assert observed["child_window"]["proper_pressure_work"] == pytest.approx(0.0, abs=1.0e-12)


def test_source_resolved_localization_stays_distinct_from_gradient_contours(prepared):
    pair = prepared["pair"]
    state = prepared["source_state"]
    _nodal, rate, bundle = prepared["source_bundle"]
    row = regions.analyze(pair, state, 0.0, bundle=bundle, nodal_rate=rate)
    observed = observer.observe(pair, state, 0.0, bundle=bundle, nodal_rate=rate)
    localization = row["localization"]
    child = localization["originally_child_prepared_pair"]
    assert child["columns"] == [2, 3]
    assert child["ancestry"] == "originally_child_prepared_source_pair"
    assert child["ancestry_is_computational_label"] is True
    assert child["independent_energy_parcel"] is False
    assert child["retained_fraction"] > 0.5
    contents = [item["content"] for item in localization["column_pairs"]]
    assert sum(contents) == pytest.approx(localization["six_column_content"])
    assert localization["six_column_total_width"] == pytest.approx(observed["localization"]["proper_width"])
    assert abs(localization["six_column_total_width"] - child["proper_width"]) > 1.0e-3
    assert localization["modal_complement_is_spatial_exterior"] is False
    assert abs(localization["modal_complement_content"] - localization["spatial_exterior_probability"]) > 1.0e-3
    assert row["positive_probability"] is True
    assert row["probability"]["per_proper_length_minimum"] >= 0.0
    assert row["probability"]["child_per_proper_length"] == pytest.approx(
        row["probability"]["child"] / observed["proper"]["child_proper_length"]
    )
    half = next(level for level in row["contours"]["levels"] if level["fraction"] == 0.5)
    assert half["kind"] == "measurement_gradient_cut"
    assert half["physical_wall"] is False
    assert half["contour_certified"] is False
    assert [half["backward"], half["forward"]] != [1.0, 3.0]
    assert row["normal_energy"]["child_fixed"]["kind"] == "measurement_cut"


def test_contour_velocity_follows_translation_and_refuses_a_shallow_slope():
    grid = galerkin.build_grid(32, gauge="conformal")
    speed = 0.37
    profile = _translation(grid, speed)
    contours = regions.density_contours(grid, profile["density"], profile["density_rate"])
    assert contours["flags"] == {
        "flat": False, "merge": False, "split": False, "multiple_dominant_peaks": False,
    }
    assert contours["contour_certified"] is False
    center = 0.5 * float(grid.length)
    for level in contours["levels"]:
        fraction = level["fraction"]
        height = fraction * (profile["offset"] + profile["amplitude"])
        cosine = (height - profile["offset"]) / profile["amplitude"]
        delta = math.acos(cosine) / profile["wave"]
        assert level["wrapped"] is False
        assert level["backward"] == pytest.approx(center - delta, abs=1.0e-8)
        assert level["forward"] == pytest.approx(center + delta, abs=1.0e-8)
        assert level["velocity_backward"]["status"] == "measured"
        assert level["velocity_forward"]["status"] == "measured"
        assert level["velocity_backward"]["velocity"] == pytest.approx(speed, abs=1.0e-8)
        assert level["velocity_forward"]["velocity"] == pytest.approx(speed, abs=1.0e-8)
        assert level["resolution_status"] == "compared"
        assert level["contour_certified"] is False
    peak_slope = float((grid.derivative @ profile["density"])[int(np.argmax(profile["density"]))])
    ambiguous = regions.implicit_level_velocity(1.0, peak_slope, peak=float(np.max(profile["density"])), spacing=float(grid.dx_q))
    assert ambiguous["status"] == "ambiguous"
    assert ambiguous["velocity"] is None
    assert ambiguous["division"] is False
    shallow = regions.implicit_level_velocity(1.0, 1.0e-8, peak=2.0, spacing=float(grid.dx_q))
    assert shallow["status"] == "ambiguous"
    measured = regions.implicit_level_velocity(-speed * 0.5, 0.5, peak=2.0, spacing=float(grid.dx_q))
    assert measured["velocity"] == pytest.approx(speed)

    length = float(grid.length)
    coordinate = np.asarray(grid.xi_q)
    wave = 2.0 * math.pi / length
    rich = 1.2 + 0.55 * np.cos(wave * coordinate) + 0.08 * np.cos(40.0 * wave * coordinate)
    rich_rate = -speed * np.asarray(grid.derivative @ rich)
    rich_contours = regions.density_contours(grid, rich, rich_rate)
    gaps = [
        level["resolution_gap_backward"] + level["resolution_gap_forward"]
        for level in rich_contours["levels"]
        if level["resolution_status"] == "compared"
    ]
    assert gaps
    assert max(gaps) > 1.0e-4
    assert rich_contours["contour_certified"] is False


def test_contour_flags_distinguish_flat_split_and_merged_peaks():
    grid = galerkin.build_grid(32, gauge="conformal")
    length = float(grid.length)
    coordinate = np.asarray(grid.xi_q, dtype=float)

    def bump(center, width):
        delta = (coordinate - center + 0.5 * length) % length - 0.5 * length
        return np.exp(-0.5 * (delta / width) ** 2)

    flat = regions.density_contours(grid, np.ones(grid.nq))
    assert flat["flags"]["flat"] is True
    assert flat["flags"]["merge"] is False
    assert flat["flags"]["split"] is False
    assert flat["levels"] == []
    assert flat["contour_certified"] is False

    separated = bump(2.0, 0.12) + 0.92 * bump(5.5, 0.12)
    between = (coordinate > 2.4) & (coordinate < 5.1)
    assert float(np.min(separated[between])) < 0.5 * float(np.max(separated))
    split = regions.density_contours(grid, separated)
    assert split["flags"]["multiple_dominant_peaks"] is True
    assert split["flags"]["split"] is True
    assert split["flags"]["merge"] is False

    close = bump(4.0, 0.25) + 0.97 * bump(4.8, 0.25)
    valley = (coordinate > 4.2) & (coordinate < 4.6)
    assert float(np.min(close[valley])) > 0.5 * float(np.max(close))
    merged = regions.density_contours(grid, close)
    assert merged["flags"]["multiple_dominant_peaks"] is True
    assert merged["flags"]["merge"] is True
    assert merged["flags"]["split"] is False


def _station_directory(directory, pair, state):
    clocks = np.asarray(nested.metrics(pair, state)["clock_rates"], dtype=float)
    arrays = {
        "Q": state.Q, "r": state.r, "chi": state.chi,
        "pi_Q": state.p_Q, "pi_r": state.p_r, "pi_chi": state.p_chi,
        "phi0": state.phi0, "phi1": state.phi1,
        "W": np.array(pair.geometry_map, copy=True),
        "source_phi0": np.array(pair.source_phi0, copy=True),
        "source_phi1": np.array(pair.source_phi1, copy=True),
        "observer_columns": np.array(pair.reference_columns, copy=True),
        "source_weights": np.array(pair.weights, copy=True),
        "clock_rates": clocks, "normal_clocks": np.zeros(3),
    }
    record = {
        "case_id": "manufactured", "nf": pair.grid.nf, "control_mode": "coupled",
        "step_cap": 0.001, "stations": [1.0, 3.0], "coordinate_time": 1.0,
        "steps": 4, "stations_reached": [1.0], "snapshot_kind": "station",
        "momentum_representation": episode.CANONICAL_PI, "verify_external_pins": False,
        "initial_state_called": False, "stability_certificate": False,
        "coarse_indices": pair.geometry_coarse_indices.tolist(),
        "child_indices": pair.geometry_child_indices.tolist(),
        "parent_indices": pair.geometry_parent_indices.tolist(),
        "source_metadata": pair.source_metadata, "geometry_metadata": pair.geometry_metadata,
        "clock_locations": [1.0, 2.0, 3.0],
        "work_ledger": episode.empty_work_ledger("coupled"),
        "dense_propagator_stored": False,
    }
    episode.commit_checkpoint(directory, record, arrays)
    return directory


def test_reader_accepts_a_station_chunk_without_reset_or_rewrite(prepared, tmp_path):
    pair = prepared["pair"]
    state = prepared["live_state"]
    directory = _station_directory(tmp_path, pair, state)
    chunk = directory / "manufactured-000000.npz"
    before = hashlib.sha256(chunk.read_bytes()).hexdigest()
    phi_before = state.phi0.tobytes()
    mode = chunk.stat().st_mode & 0o222
    report = regions.read_station_regions(directory, "manufactured")
    after = hashlib.sha256(chunk.read_bytes()).hexdigest()
    assert state.phi0.tobytes() == phi_before
    assert mode == 0
    assert before == after
    assert report["source_reset"] is False
    assert report["historical_rewrite"] is False
    assert report["evolved"] is False
    assert report["sealed_preparations_parsed"] is False
    station = report["stations"][0]
    assert station["snapshot_kind"] == "station"
    assert station["coordinate_time"] == pytest.approx(1.0)
    assert station["state_equals_prepared_source"] is False
    assert station["initial_state_called"] is False
    assert station["chunk_sha256"]["npz"] == before
    direct = regions.analyze(pair, state, 1.0, bundle=prepared["live_bundle"][2], nodal_rate=prepared["live_bundle"][1])
    assert station["regions"]["probability"]["total"] == pytest.approx(direct["probability"]["total"])
    assert station["regions"]["normal_energy"]["child_fixed"]["rate"] == pytest.approx(
        direct["normal_energy"]["child_fixed"]["rate"]
    )


def test_cli_check_and_read_are_creation_only_and_leave_sealed_inputs_unchanged(prepared, tmp_path):
    before = _sealed_stat()
    check = subprocess.run(
        [str(PYTHON), str(ROOT / "scripts/lab.py"), "scripts/derive_nsc_discovery_regions.py", "--check"],
        check=True, capture_output=True, text=True,
    )
    record = json.loads(check.stdout)
    assert record["schema"] == regions.CHECK_SCHEMA
    assert record["arrays_loaded"] is False
    assert record["sealed_preparations_parsed"] is False
    assert record["evolved"] is False
    assert record["output_written"] is False
    assert record["renewal_asserted"] is False
    assert record["contour_certified"] is False
    assert record["settings"]["child_source_columns"] == [2, 3]
    assert record["settings"]["piston_work_on_measurement_cut"] is False
    assert len(record["renewal_classification_requires"]) == 6
    produced = record["producer_hashes"]["src/recursive_horizons/nsc_discovery_regions.py"]
    assert produced == regions.sha256_file(LAB / "src/recursive_horizons/nsc_discovery_regions.py")
    assert {item["path"] for item in record["sealed_inputs"]} == set(regions.SEALED_INPUTS)
    assert all(item["parsed"] is False for item in record["sealed_inputs"])
    directory = _station_directory(tmp_path / "episode", prepared["pair"], prepared["live_state"])
    output = tmp_path / "regions.json"
    read = subprocess.run(
        [str(PYTHON), str(ROOT / "scripts/lab.py"), "scripts/derive_nsc_discovery_regions.py",
         "--read", str(directory), "--case", "manufactured", "--output", str(output)],
        check=True, capture_output=True, text=True,
    )
    written = json.loads(output.read_text())
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    assert written["schema"] == regions.SCHEMA
    assert written["stations"][0]["snapshot_kind"] == "station"
    assert json.loads(read.stdout)["stations"][0]["chunk_sha256"] == written["stations"][0]["chunk_sha256"]
    refused = subprocess.run(
        [str(PYTHON), str(ROOT / "scripts/lab.py"), "scripts/derive_nsc_discovery_regions.py",
         "--read", str(directory), "--case", "manufactured", "--output", str(output)],
        capture_output=True, text=True,
    )
    assert refused.returncode != 0
    assert "refusing to replace" in refused.stderr
    assert hashlib.sha256(output.read_bytes()).hexdigest() == digest
    inside = subprocess.run(
        [str(PYTHON), str(ROOT / "scripts/lab.py"), "scripts/derive_nsc_discovery_regions.py",
         "--read", str(directory), "--case", "manufactured", "--output", str(directory / "extra.json")],
        capture_output=True, text=True,
    )
    assert inside.returncode != 0
    assert _sealed_stat() == before
    assert not (LAB / "results/development/nsc-discovery-regions-v1.json").exists()
