"""Independent metric controls and saved-sample checks. No evolution."""
import hashlib
import json
import os
from pathlib import Path

for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_var, "1")

import numpy as np
import pytest

from recursive_horizons.nsc_spherical_coupling import CALIBRATION, PERIOD, periodic_derivative
from recursive_horizons.nsc_spherical_null_expansion import (
    CPU_BUDGET_S,
    EPISODE_JSON,
    EPISODE_NPZ,
    JSON_BYTE_LIMIT,
    RECORD_PATH,
    boosted_expansions,
    build_record,
    circular_intervals,
    expansions_from_stored_normal,
    null_expansions,
    spectral_derivative,
    static_clock,
)

_MODULE = Path(__file__).resolve().parents[1] / "src" / "recursive_horizons" / "nsc_spherical_null_expansion.py"
_REVIEW = Path(__file__).resolve().parents[1] / "results" / "development" / "nsc-local-boundary-review-v1.json"
_RESPONSE_JSON = Path(__file__).resolve().parents[1] / "results" / "development" / "nsc-coupled-local-response-v1.json"
_RESPONSE_NPZ = Path(__file__).resolve().parents[1] / "results" / "development" / "nsc-coupled-local-response-v1.npz"
_RESPONSE_MODULE = Path(__file__).resolve().parents[1] / "src" / "recursive_horizons" / "nsc_coupled_local_response.py"
_NULL_DOC = Path(__file__).resolve().parents[1] / "docs" / "nsc-spherical-null-expansion.md"
_RESPONSE_DOC = Path(__file__).resolve().parents[1] / "docs" / "nsc-coupled-local-response.md"


@pytest.fixture(scope="module")
def measurement():
    return build_record()


@pytest.fixture(scope="module")
def saved_record():
    return json.loads(RECORD_PATH.read_text())


def _independent_inverse(r, length_density, radial_factor, shift):
    """2x2 inverse from det(g), written apart from the module's inverse helper."""
    area = r * r
    conformal_square = radial_factor * radial_factor
    g_tt = area * (length_density ** 2 - conformal_square * shift ** 2)
    g_tx = -area * conformal_square * shift
    g_xx = -area * conformal_square
    determinant = g_tt * g_xx - g_tx ** 2
    inverse_tt = g_xx / determinant
    inverse_xx = g_tt / determinant
    inverse_tx = -g_tx / determinant
    return g_tt, g_tx, g_xx, inverse_tt, inverse_tx, inverse_xx


def _independent_norms(r, length_density, radial_factor, shift):
    g_tt, g_tx, g_xx, _inverse_tt, _inverse_tx, _inverse_xx = _independent_inverse(
        r, length_density, radial_factor, shift
    )
    lapse = r * length_density
    radial_metric = r * radial_factor
    normal_t = 1.0 / lapse
    normal_x = -shift / lapse
    radial_x = 1.0 / radial_metric
    plus_t, plus_x = normal_t, normal_x + radial_x
    minus_t, minus_x = normal_t, normal_x - radial_x

    def form(left_t, left_x, right_t, right_x):
        return (
            g_tt * left_t * right_t
            + g_tx * (left_t * right_x + left_x * right_t)
            + g_xx * left_x * right_x
        )

    return {
        "nn": form(normal_t, normal_x, normal_t, normal_x),
        "ee": form(0.0, radial_x, 0.0, radial_x),
        "plus": form(plus_t, plus_x, plus_t, plus_x),
        "minus": form(minus_t, minus_x, minus_t, minus_x),
        "cross": form(plus_t, plus_x, minus_t, minus_x),
    }


def test_postprocess_does_not_own_evolution():
    text = _MODULE.read_text()
    for forbidden in ("rk4_step", "evolve_episode", "solve_initial_radius", "full_source_rho"):
        assert forbidden not in text
    assert "frame_chi" not in text
    assert "weyl" not in text


def test_inverse_metric_product_null_norms_and_boost():
    cases = (
        (2.0, 1.0, 1.0, 0.0, -0.5, 0.0, -0.25, -0.25),
        (2.0, 1.0, 1.0, 0.0, 0.5, 0.0, 0.25, 0.25),
        (2.0, 1.0, 1.0, 0.0, 0.0, 0.5, 0.25, -0.25),
        (2.0, 3.0, 4.0, 0.5, -1.25, 0.25, None, None),
    )
    generator = np.random.default_rng(19)
    extra = [
        tuple(float(value) for value in row) + (None, None)
        for row in generator.uniform(
            [0.3, 0.2, 0.2, -1.0, -2.0, -2.0],
            [3.0, 2.0, 2.0, 2.0, 2.0, 2.0],
            size=(6, 6),
        )
    ]
    for radius, lapse_factor, conformal, shift, r_t, r_x, expected_plus, expected_minus in list(cases) + extra:
        theta_plus, theta_minus, product = null_expansions(
            radius, r_t, r_x, lapse_factor, conformal, shift
        )
        _g_tt, _g_tx, _g_xx, inverse_tt, inverse_tx, inverse_xx = _independent_inverse(
            radius, lapse_factor, conformal, shift
        )
        gradient = inverse_tt * r_t ** 2 + 2.0 * inverse_tx * r_t * r_x + inverse_xx * r_x ** 2
        assert abs(product - 4.0 / radius ** 2 * gradient) < 1e-12
        norms = _independent_norms(radius, lapse_factor, conformal, shift)
        assert abs(norms["nn"] - 1.0) < 1e-12
        assert abs(norms["ee"] + 1.0) < 1e-12
        assert abs(norms["plus"]) < 1e-12
        assert abs(norms["minus"]) < 1e-12
        assert abs(norms["cross"] - 2.0) < 1e-12
        for factor in (0.5, 3.0, 8.0):
            scaled_plus, scaled_minus, scaled_product = boosted_expansions(
                theta_plus, theta_minus, factor
            )
            assert abs(scaled_product - product) < 1e-12
            if theta_plus != 0.0:
                assert np.sign(scaled_plus) == np.sign(theta_plus)
            if theta_minus != 0.0:
                assert np.sign(scaled_minus) == np.sign(theta_minus)
        if expected_plus is not None:
            assert abs(float(theta_plus) - expected_plus) < 1e-12
            assert abs(float(theta_minus) - expected_minus) < 1e-12
        flipped_product = (-theta_plus) * (-theta_minus)
        assert abs(flipped_product - product) < 1e-15
        if theta_plus != 0.0:
            assert (-theta_plus) * theta_plus < 0.0


def test_derivative_matches_owner_and_kills_nyquist():
    points = 64
    dense = periodic_derivative(points, PERIOD)
    coordinate = np.arange(points) * (PERIOD / points)
    field = np.sin(4.0 * np.pi * coordinate / PERIOD) + 0.3 * np.cos(np.pi * np.arange(points))
    assert np.max(np.abs(spectral_derivative(field, PERIOD) - dense @ field)) < 1e-11
    nyquist = np.cos(np.pi * np.arange(points))
    assert np.max(np.abs(spectral_derivative(nyquist, PERIOD))) < 1e-12


def test_clock_arc_uses_calibration_literals():
    _coordinate, length_density, shift = static_clock(8)
    assert length_density[0] == pytest.approx(CALIBRATION["b0"])
    assert shift[0] == pytest.approx(CALIBRATION["beta0"])
    assert length_density[2] == pytest.approx(CALIBRATION["b0"] * CALIBRATION["Omega"] ** 2)
    assert shift[2] == pytest.approx(CALIBRATION["beta0"] * CALIBRATION["Omega"] ** 2)


def test_circular_intervals_wrap_and_split():
    assert circular_intervals(np.array([False, False, False])) == []
    split = circular_intervals(np.array([False, True, True, False, True]))
    assert [(item["start_index"], item["end_index"], item["count"]) for item in split] == [
        (1, 2, 2),
        (4, 4, 1),
    ]
    wrapped = circular_intervals(np.array([True, True, False, True]))
    assert wrapped[0]["wraps"] is True
    assert wrapped[0]["count"] == 3


def test_saved_cases_change_character_under_positive_chart(measurement, saved_record):
    assert measurement["evolution_performed"] is False
    assert measurement["inheritance_recomputed"] is False
    assert measurement["chi_used"] is False
    assert measurement["historical_1e-8_used_as_veto"] is False
    assert measurement["global_event_horizon_claimed"] is False
    assert measurement["black_hole_birth_claimed"] is False
    assert measurement["child_region_claimed"] is False
    assert measurement["regeneration_claimed"] is False
    assert measurement["sampling_is_continuous_horizon_proof"] is False
    assert measurement["cpu_budget_exceeded"] is False
    assert measurement["cpu_seconds"] <= CPU_BUDGET_S
    assert measurement["status"] == "MEASURED_SAMPLE_TRAPPING_CHARACTER_CHANGES"
    assert measurement["pattern"]["positive_geometry"] is True
    assert measurement["pattern"]["sample_character_changes"] is True
    assert measurement["pattern"]["arcs_track_areal_extrema"] is True
    assert measurement["pattern"]["strict_margin_agrees"] is True
    agreement = measurement["source_bindings"]["hash_agreement"]
    assert all(agreement.values())
    assert measurement["source_bindings"]["episode_files_unchanged"] is True
    assert measurement["metadata_normal_gap"] < 1e-12
    assert measurement["controls"]["inverse_metric_product_gap"] < 1e-12
    assert measurement["controls"]["null_norm_gap"] < 1e-12
    assert measurement["controls"]["boost_product_gap"] < 1e-12
    assert measurement["controls"]["analytic_sign_gap"] < 1e-12
    assert measurement["controls"]["positive_boost_preserves_each_sign"] is True
    assert measurement["controls"]["derivative_owner"]["field_gap"] < 1e-11
    comparison = measurement["comparisons"]
    assert comparison["dt_raw_counts_identical"] is True
    assert comparison["shared_node_sign_disagreements_max"] == 0
    assert comparison["dt_nf512_theta_plus_max_abs_gap"] < 1e-8
    assert measurement["clearance"]["clearance_ratio"] > 50.0
    by_name = {run["name"]: run for run in measurement["runs"]}
    fine = by_name["nf512_dt_0_0005"]
    coarse = by_name["nf256_dt_0_0005"]
    assert fine["frames"][0]["raw_counts"] == {
        "trapped": 0,
        "anti_trapped": 0,
        "untrapped": 2048,
        "marginal_indicator": 0,
    }
    assert fine["frames"][1]["raw_counts"]["trapped"] == 4
    assert fine["frames"][1]["raw_counts"]["anti_trapped"] == 8
    assert fine["frames"][-1]["raw_counts"]["trapped"] == 45
    assert fine["frames"][-1]["raw_counts"]["anti_trapped"] == 182
    assert fine["frames"][-1]["raw_counts"]["untrapped"] == 1821
    assert coarse["frames"][-1]["raw_counts"]["trapped"] == 23
    assert coarse["frames"][-1]["raw_counts"]["anti_trapped"] == 91
    assert fine["frames"][0]["product_max"] < 0.0
    assert fine["frames"][-1]["product_max"] > 0.0
    assert fine["frames"][0]["critical_point_estimates"]
    assert all(
        point["label"] == "marginal_within_margin"
        for point in fine["frames"][0]["critical_point_estimates"]
    )
    assert fine["frames"][-1]["product_crossing_count"] == 4
    assert fine["positive_r"] is True and fine["positive_Q"] is True
    assert fine["observer_proper_gap"] < 1e-14
    assert fine["high_mode_derivative_bound"] < 1e-8
    final_norm = fine["frames"][-1]["norm_gaps"]["g_nn_minus_one"]
    assert final_norm < 1e-12
    initial = fine["frames"][0]
    final = fine["frames"][-1]
    margin = measurement["margin"]["absolute_expansion_margin"]
    assert margin == pytest.approx(1.8325977887805045e-06)
    saved_fine = next(run for run in saved_record["runs"] if run["name"] == "nf512_dt_0_0005")
    peak = initial["areal_maximum"]
    saved_peak = saved_fine["frames"][0]["areal_maximum"]
    assert peak["index"] == saved_peak["index"] == 395
    assert abs(peak["r_x"]) > 1e-4
    recomputed_plus, recomputed_minus, recomputed_product, _spatial = expansions_from_stored_normal(
        peak["r"],
        peak["normal_velocity"],
        peak["r_x"],
        initial["Q_min"],
    )
    assert float(recomputed_plus) == pytest.approx(saved_peak["theta_plus"], abs=1e-14)
    assert float(recomputed_minus) == pytest.approx(saved_peak["theta_minus"], abs=1e-14)
    assert float(recomputed_plus) == pytest.approx(3.113501917378642e-4, rel=1e-9)
    assert float(recomputed_minus) == pytest.approx(-3.113501506382507e-4, rel=1e-9)
    assert float(recomputed_plus) * float(recomputed_minus) < 0.0
    assert float(recomputed_product) < 0.0
    critical = initial["critical_point_estimates"]
    theta_abs = max(abs(point["theta_both"]) for point in critical)
    assert theta_abs == pytest.approx(1.9636126344543498e-11, rel=1e-9)
    assert abs(theta_abs - 1.96e-11) < 5e-14
    assert theta_abs < margin
    assert abs(initial["normal_max"]) < 1e-9
    assert abs(final["normal_max"]) > 1e-2
    assert abs(final["areal_maximum"]["theta_plus"]) > margin
    assert abs(final["areal_maximum"]["theta_minus"]) > margin
    assert final["areal_maximum"]["theta_plus"] < 0.0 and final["areal_maximum"]["theta_minus"] < 0.0
    assert final["areal_minimum"]["theta_plus"] > margin
    assert final["areal_minimum"]["theta_minus"] > margin
    assert final["same_sign_min_abs_expansion"] > margin
    assert measurement["black_hole_birth_claimed"] is False
    assert measurement["global_event_horizon_claimed"] is False


def test_record_distinguishes_marginal_indicator_from_global_horizon(measurement, saved_record):
    assert saved_record["schema"] == measurement["schema"]
    assert saved_record["status"] == measurement["status"]
    assert saved_record["finding"] == measurement["finding"]
    assert saved_record["global_event_horizon_claimed"] is False
    assert saved_record["sampling_is_continuous_horizon_proof"] is False
    assert "global event horizon" in " ".join(saved_record["open"]).lower() or any(
        "event horizon" in item.lower() for item in saved_record["open"]
    )
    assert RECORD_PATH.stat().st_size < JSON_BYTE_LIMIT
    episode_hash = hashlib.sha256(EPISODE_NPZ.read_bytes()).hexdigest()
    assert episode_hash == saved_record["source_bindings"]["hashes"]["episode_npz"]
    assert episode_hash == saved_record["source_bindings"]["episode_npz_sha256_recorded_in_metadata"]
    json_hash = hashlib.sha256(EPISODE_JSON.read_bytes()).hexdigest()
    assert json_hash == saved_record["source_bindings"]["hashes"]["episode_json"]
    frozen_module = saved_record["source_bindings"]["module_sha256"]
    review = json.loads(_REVIEW.read_text())
    assert review["immutable_v1"]["null_expansion_recorded_module_sha256"] == frozen_module
    assert hashlib.sha256(_MODULE.read_bytes()).hexdigest() == review["current_owner_sha256"]["nsc_spherical_null_expansion.py"]
    assert hashlib.sha256(_MODULE.read_bytes()).hexdigest() != frozen_module


def test_flipped_boost_return_fails_the_sign_flag(monkeypatch):
    import recursive_horizons.nsc_spherical_null_expansion as owner

    real = owner.boosted_expansions

    def flipped(theta_plus, theta_minus, scale):
        plus, minus, product = real(theta_plus, theta_minus, scale)
        return -plus, -minus, product

    monkeypatch.setattr(owner, "boosted_expansions", flipped)
    controls = owner.independent_controls()
    assert controls["positive_boost_preserves_each_sign"] is False
    assert controls["boost_product_gap"] < 1e-12


def test_observer_label_names_initial_mode_columns():
    from recursive_horizons.nsc_coupled_local_response import finite_domain

    label = finite_domain(
        {
            "case": "nf512_dt_0_0005",
            "frame_time": np.zeros(11),
            "interpolation": {"midpoint_linear_scale": 0.0},
            "time_resolution_geometry_gap": 0.0,
        },
        {"t_start": 0.0, "t_end": 0.05, "steps": 10},
        {"nf": 512, "nq": 2048, "length": 8.0},
        0.0,
        True,
    )["observer"]
    assert label == "initial mode columns 0 and 1, weights 0.75 and 0.75, no phase QR"
    assert "plus and minus" not in label


def test_successor_review_binds_immutable_v1(measurement, saved_record):
    review = json.loads(_REVIEW.read_text())
    immutable = review["immutable_v1"]
    assert immutable["episode_json_sha256"] == hashlib.sha256(EPISODE_JSON.read_bytes()).hexdigest()
    assert immutable["episode_npz_sha256"] == hashlib.sha256(EPISODE_NPZ.read_bytes()).hexdigest()
    assert immutable["null_expansion_json_sha256"] == hashlib.sha256(RECORD_PATH.read_bytes()).hexdigest()
    assert immutable["null_expansion_recorded_module_sha256"] == saved_record["source_bindings"]["module_sha256"]
    assert immutable["coupled_response_json_sha256"] == hashlib.sha256(_RESPONSE_JSON.read_bytes()).hexdigest()
    assert immutable["coupled_response_npz_sha256"] == hashlib.sha256(_RESPONSE_NPZ.read_bytes()).hexdigest()
    stored = json.loads(_RESPONSE_JSON.read_text())
    assert stored["finite_approximation_domain"]["observer"] == immutable["stored_observer_label"]
    assert immutable["stored_observer_label"] == "initial region-0 plus and minus source columns, unphased"
    owners = review["current_owner_sha256"]
    assert owners["nsc_spherical_null_expansion.py"] == hashlib.sha256(_MODULE.read_bytes()).hexdigest()
    assert owners["nsc_coupled_local_response.py"] == hashlib.sha256(_RESPONSE_MODULE.read_bytes()).hexdigest()
    assert owners["test_nsc_spherical_null_expansion.py"] == hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    assert owners["nsc-spherical-null-expansion.md"] == hashlib.sha256(_NULL_DOC.read_bytes()).hexdigest()
    assert owners["nsc-coupled-local-response.md"] == hashlib.sha256(_RESPONSE_DOC.read_bytes()).hexdigest()
    saved_fine = next(run for run in saved_record["runs"] if run["name"] == "nf512_dt_0_0005")
    peak = saved_fine["frames"][0]["areal_maximum"]
    assert review["controls"]["positive_boost_preserves_each_sign"] is True
    assert review["controls"]["wrong_return_preserves_each_sign"] is False
    assert review["controls"]["boost_product_gap"] == measurement["controls"]["boost_product_gap"]
    assert review["controls"]["t0_areal_maximum_theta_plus"] == peak["theta_plus"]
    assert review["controls"]["t0_areal_maximum_theta_minus"] == peak["theta_minus"]
    assert review["controls"]["absolute_expansion_margin"] == saved_record["margin"]["absolute_expansion_margin"]
    independent = Path(__file__).resolve().parents[1] / "tests" / "test_nsc_local_boundary_independent.py"
    assert review["reference_unchanged"]["independent_test_sha256"] == hashlib.sha256(independent.read_bytes()).hexdigest()
