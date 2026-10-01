"""Independent checks for the conditional local response and the spherical null samples.

The nf512 Volterra campaign is not run again. Saved occupation series are
re-reduced from the payload, and the autonomous endpoint is re-projected
from the episode spinors. Null expansions are derived here from the line
element. A four-dimensional ODE checks the retarded kernel, the exterior
drive, and the memory omission. Covariance is not treated as a stress.
"""
from __future__ import annotations

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
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from recursive_horizons.nsc_coupled_local_response import (
    hermite_rows,
    interpolate_rows,
    region_observer,
    request_identity,
    resume_disposition,
    source_id,
    spinor_local_moments,
)
from recursive_horizons.nsc_evolving_reduction import (
    causal_memory_kernel,
    evolve_retained_region,
    exterior_kernels,
    fixed_observer_frame,
    prescribed_six_mode_problem,
)
from recursive_horizons.nsc_spherical_coupling import (
    CALIBRATION,
    PERIOD,
    periodic_derivative,
    profile_coordinate,
)
from recursive_horizons.nsc_spherical_null_expansion import static_clock

_LAB = Path(__file__).resolve().parents[1]
_EPISODE = _LAB / "results" / "development" / "nsc-spherical-feedback-episode-v1.npz"
_NULL_RECORD = _LAB / "results" / "development" / "nsc-spherical-null-expansion-v1.json"
_RESPONSE_JSON = _LAB / "results" / "development" / "nsc-coupled-local-response-v1.json"
_RESPONSE_NPZ = _LAB / "results" / "development" / "nsc-coupled-local-response-v1.npz"
_RUNS = (
    "nf512_dt_0_0005",
    "nf512_dt_0_00025",
    "nf256_dt_0_0005",
    "nf256_dt_0_00025",
)


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _clock(count):
    """Calibration lapse factor and shift. L and beta are not evolved."""
    count = int(count)
    coordinate = np.arange(count, dtype=float) * (PERIOD / count)
    profile = profile_coordinate(coordinate, PERIOD)
    scale = np.exp(np.log(CALIBRATION["Omega"]) * profile)
    return coordinate, CALIBRATION["b0"] * scale, CALIBRATION["beta0"] * scale


def _fft_derivative(values, length):
    """Periodic spectral derivative with the Nyquist symbol removed."""
    samples = np.asarray(values, dtype=float)
    count = samples.shape[-1]
    spectrum = np.fft.fft(samples, axis=-1)
    modes = np.fft.fftfreq(count) * count
    wavenumber = np.array(2.0 * np.pi * modes / float(length), dtype=float, copy=True)
    wavenumber[count // 2] = 0.0
    return np.real(np.fft.ifft(spectrum * (1j * wavenumber), axis=-1))


def _inverse_metric(radius, lapse_factor, conformal, shift):
    """Inverse (t, x) block from det of the expanded +--- line element.

    g = r^2 [ L^2 dt^2 - Q^2 (dx + beta dt)^2 - dOmega^2 ].
    Angular derivatives of r are zero, so they do not enter this block.
    """
    area = np.asarray(radius, dtype=float) ** 2
    conformal_square = np.asarray(conformal, dtype=float) ** 2
    shift = np.asarray(shift, dtype=float)
    lapse_factor = np.asarray(lapse_factor, dtype=float)
    g_tt = area * (lapse_factor ** 2 - conformal_square * shift ** 2)
    g_tx = -area * conformal_square * shift
    g_xx = -area * conformal_square
    determinant = g_tt * g_xx - g_tx ** 2
    inverse_tt = g_xx / determinant
    inverse_tx = -g_tx / determinant
    inverse_xx = g_tt / determinant
    return g_tt, g_tx, g_xx, inverse_tt, inverse_tx, inverse_xx


def _quadratic(g_tt, g_tx, g_xx, left_t, left_x, right_t, right_x):
    return (
        g_tt * left_t * right_t
        + g_tx * (left_t * right_x + left_x * right_t)
        + g_xx * left_x * right_x
    )


def _null_data(radius, r_t, r_x, lapse_factor, conformal, shift):
    """theta_pm = 2/r k_pm(r) for k_pm = n ± e1, and the inverse-metric product."""
    radius = np.asarray(radius, dtype=float)
    lapse_factor = np.asarray(lapse_factor, dtype=float)
    conformal = np.asarray(conformal, dtype=float)
    shift = np.asarray(shift, dtype=float)
    r_t = np.asarray(r_t, dtype=float)
    r_x = np.asarray(r_x, dtype=float)
    lapse = radius * lapse_factor
    radial_metric = radius * conformal
    normal_t = 1.0 / lapse
    normal_x = -shift / lapse
    radial_x = 1.0 / radial_metric
    g_tt, g_tx, g_xx, inverse_tt, inverse_tx, inverse_xx = _inverse_metric(
        radius, lapse_factor, conformal, shift
    )
    plus_t, plus_x = normal_t, normal_x + radial_x
    minus_t, minus_x = normal_t, normal_x - radial_x
    normal_velocity = (r_t - shift * r_x) / lapse
    spatial = r_x / radial_metric
    theta_plus = 2.0 / radius * (normal_velocity + spatial)
    theta_minus = 2.0 / radius * (normal_velocity - spatial)
    gradient = (
        inverse_tt * r_t ** 2
        + 2.0 * inverse_tx * r_t * r_x
        + inverse_xx * r_x ** 2
    )
    norms = {
        "nn": _quadratic(g_tt, g_tx, g_xx, normal_t, normal_x, normal_t, normal_x),
        "ee": _quadratic(g_tt, g_tx, g_xx, 0.0, radial_x, 0.0, radial_x),
        "plus": _quadratic(g_tt, g_tx, g_xx, plus_t, plus_x, plus_t, plus_x),
        "minus": _quadratic(g_tt, g_tx, g_xx, minus_t, minus_x, minus_t, minus_x),
        "cross": _quadratic(g_tt, g_tx, g_xx, plus_t, plus_x, minus_t, minus_x),
        "time_component": plus_t,
    }
    return theta_plus, theta_minus, gradient, norms


def _counts(theta_plus, theta_minus):
    trapped = (theta_plus < 0.0) & (theta_minus < 0.0)
    anti = (theta_plus > 0.0) & (theta_minus > 0.0)
    untrapped = ((theta_plus > 0.0) & (theta_minus < 0.0)) | ((theta_plus < 0.0) & (theta_minus > 0.0))
    marginal = ~(trapped | anti | untrapped)
    return {
        "trapped": int(np.sum(trapped)),
        "anti_trapped": int(np.sum(anti)),
        "untrapped": int(np.sum(untrapped)),
        "marginal": int(np.sum(marginal)),
        "trapped_mask": trapped,
        "anti_mask": anti,
    }


def _circular_groups(mask):
    selected = np.flatnonzero(np.asarray(mask, dtype=bool))
    if selected.size == 0:
        return []
    splits = np.where(np.diff(selected) > 1)[0]
    groups = np.split(selected, splits + 1)
    count = int(mask.size)
    if len(groups) > 1 and int(groups[0][0]) == 0 and int(groups[-1][-1]) == count - 1:
        groups = [np.concatenate((groups[-1], groups[0]))] + groups[1:-1]
    return groups


def _critical_thetas(radius, normal, r_x):
    slope = np.asarray(r_x, dtype=float)
    found = []
    for index in range(slope.size):
        neighbor = (index + 1) % slope.size
        left = float(slope[index])
        right = float(slope[neighbor])
        if left == 0.0 or left * right >= 0.0:
            continue
        weight = abs(left) / (abs(left) + abs(right))
        radius_zero = (1.0 - weight) * float(radius[index]) + weight * float(radius[neighbor])
        velocity = (1.0 - weight) * float(normal[index]) + weight * float(normal[neighbor])
        found.append(2.0 * velocity / radius_zero)
    return np.asarray(found, dtype=float)


def _load_run(payload, name):
    radius = np.array(payload[name + "_frame_r"], dtype=float)
    conformal = np.array(payload[name + "_frame_Q"], dtype=float)
    proper = np.array(payload[name + "_frame_proper"], dtype=float)
    k_perp = np.array(payload[name + "_frame_K_perp"], dtype=float)
    times = np.array(payload[name + "_frame_time"], dtype=float)
    _coordinate, lapse_factor, shift = _clock(radius.shape[1])
    r_x = _fft_derivative(radius, PERIOD)
    r_t = proper * (radius * lapse_factor) + shift * r_x
    theta_plus, theta_minus, gradient, norms = _null_data(
        radius, r_t, r_x, lapse_factor, conformal, shift
    )
    return {
        "radius": radius,
        "conformal": conformal,
        "proper": proper,
        "k_perp": k_perp,
        "times": times,
        "lapse_factor": lapse_factor,
        "shift": shift,
        "r_x": r_x,
        "r_t": r_t,
        "theta_plus": theta_plus,
        "theta_minus": theta_minus,
        "gradient": gradient,
        "final_norms": {key: np.asarray(value) for key, value in norms.items()},
    }


@pytest.fixture(scope="module")
def saved():
    before = {
        "episode": _sha256(_EPISODE),
        "null": _sha256(_NULL_RECORD),
        "response": _sha256(_RESPONSE_JSON),
    }
    with np.load(_EPISODE, allow_pickle=False) as payload:
        runs = {name: _load_run(payload, name) for name in _RUNS}
        spinors = {
            "phi0": np.array(payload["nf512_initial_phi0"]),
            "phi1": np.array(payload["nf512_initial_phi1"]),
            "final_phi0": np.array(payload["nf512_dt_0_0005_final_phi0"]),
            "final_phi1": np.array(payload["nf512_dt_0_0005_final_phi1"]),
            "weights": np.array(payload["nf512_occupations"], dtype=float),
            "final_q_dot": np.array(payload["nf512_dt_0_0005_final_rate_Q_dot"], dtype=float),
            "replay_time": float(payload["replay_time"]),
            "files": tuple(payload.files),
        }
    with np.load(_RESPONSE_NPZ, allow_pickle=False) as series:
        arrays = {key: np.array(series[key]) for key in series.files}
    bundle = {
        "runs": runs,
        "spinors": spinors,
        "arrays": arrays,
        "null_record": json.loads(_NULL_RECORD.read_text()),
        "response": json.loads(_RESPONSE_JSON.read_text()),
        "before": before,
    }
    yield bundle
    assert _sha256(_EPISODE) == before["episode"]
    assert _sha256(_NULL_RECORD) == before["null"]
    assert _sha256(_RESPONSE_JSON) == before["response"]


def test_null_formulas_cover_signature_product_norms_and_boost():
    samples = [
        (2.0, 1.0, 1.0, 0.0, -0.5, 0.0, -0.25, -0.25),
        (2.0, 1.0, 1.0, 0.0, 0.5, 0.0, 0.25, 0.25),
        (2.0, 1.0, 1.0, 0.0, 0.0, 0.5, 0.25, -0.25),
        (2.0, 3.0, 4.0, 0.5, -1.25, 0.25, None, None),
    ]
    generator = np.random.default_rng(11)
    for row in generator.uniform([0.4, 0.2, 0.2, -1.0, -1.5, -1.5], [3.0, 2.0, 2.0, 1.5, 1.5, 1.5], size=(8, 6)):
        samples.append(tuple(float(value) for value in row) + (None, None))
    for radius, lapse_factor, conformal, shift, r_t, r_x, expected_plus, expected_minus in samples:
        theta_plus, theta_minus, gradient, norms = _null_data(
            radius, r_t, r_x, lapse_factor, conformal, shift
        )
        product = float(theta_plus * theta_minus)
        assert abs(product - float(4.0 / radius ** 2 * gradient)) < 1e-12
        assert abs(float(norms["nn"]) - 1.0) < 1e-12
        assert abs(float(norms["ee"]) + 1.0) < 1e-12
        assert abs(float(norms["plus"])) < 1e-12
        assert abs(float(norms["minus"])) < 1e-12
        assert abs(float(norms["cross"]) - 2.0) < 1e-12
        assert float(norms["time_component"]) > 0.0
        for factor in (0.5, 2.0, 6.0):
            scaled_product = float((factor * theta_plus) * (theta_minus / factor))
            assert abs(scaled_product - product) < 1e-12
            if float(theta_plus) != 0.0:
                assert np.sign(factor * theta_plus) == np.sign(theta_plus)
            if float(theta_minus) != 0.0:
                assert np.sign(theta_minus / factor) == np.sign(theta_minus)
        if expected_plus is not None:
            assert abs(float(theta_plus) - expected_plus) < 1e-12
            assert abs(float(theta_minus) - expected_minus) < 1e-12
            if expected_plus * expected_minus < 0.0:
                assert product < 0.0
                assert float(gradient) < 0.0
    points = 64
    coordinate = np.arange(points) * (PERIOD / points)
    field = np.sin(2.0 * np.pi * coordinate / PERIOD) + 0.3 * np.cos(np.pi * np.arange(points))
    dense = periodic_derivative(points, PERIOD)
    assert np.max(np.abs(_fft_derivative(field, PERIOD) - dense @ field)) < 1e-12
    assert np.max(np.abs(_fft_derivative(np.cos(np.pi * np.arange(points)), PERIOD))) < 1e-12


def test_saved_samples_are_local_future_characters(saved):
    primary = saved["runs"]["nf512_dt_0_0005"]
    record = saved["null_record"]
    assert record["status"] == "MEASURED_SAMPLE_TRAPPING_CHARACTER_CHANGES"
    assert record["evolution_performed"] is False
    assert record["chi_used"] is False
    assert record["global_event_horizon_claimed"] is False
    assert record["black_hole_birth_claimed"] is False
    assert record["child_region_claimed"] is False
    assert record["regeneration_claimed"] is False
    assert record["sampling_is_continuous_horizon_proof"] is False
    assert record["source_bindings"]["hashes"]["episode_npz"] == saved["before"]["episode"]
    assert all(record["source_bindings"]["hash_agreement"].values())
    assert record["source_bindings"]["episode_files_unchanged"] is True

    expected_times = np.arange(0.0, 0.05 + 1e-12, 0.005)
    product_gap = 0.0
    proper_gap = 0.0
    norm_gap = 0.0
    for name, run in saved["runs"].items():
        assert run["times"].shape == expected_times.shape
        assert np.max(np.abs(run["times"] - expected_times)) < 1e-12
        assert np.min(run["radius"]) > 0.0 and np.min(run["conformal"]) > 0.0
        assert np.min(run["lapse_factor"]) > 0.0
        product = run["theta_plus"] * run["theta_minus"]
        product_gap = max(product_gap, float(np.max(np.abs(product - 4.0 / run["radius"] ** 2 * run["gradient"]))))
        proper_gap = max(proper_gap, float(np.max(np.abs(run["k_perp"] * run["radius"] - run["proper"]))))
        final = run["final_norms"]
        norm_gap = max(
            norm_gap,
            float(np.max(np.abs(final["nn"] - 1.0))),
            float(np.max(np.abs(final["ee"] + 1.0))),
            float(np.max(np.abs(final["plus"]))),
            float(np.max(np.abs(final["minus"]))),
            float(np.max(np.abs(final["cross"] - 2.0))),
        )
        assert float(np.min(final["time_component"])) > 0.0
        owner_clock = static_clock(run["radius"].shape[1])
        assert np.allclose(owner_clock[1], run["lapse_factor"], rtol=0.0, atol=0.0)
        assert np.allclose(owner_clock[2], run["shift"], rtol=0.0, atol=0.0)
        initial = _counts(run["theta_plus"][0], run["theta_minus"][0])
        assert initial["trapped"] == 0 and initial["anti_trapped"] == 0
        assert initial["marginal"] == 0
        assert initial["untrapped"] == run["radius"].shape[1]
        assert float(np.max(run["theta_plus"][0] * run["theta_minus"][0])) < 0.0
        critical = _critical_thetas(run["radius"][0], run["proper"][0], run["r_x"][0])
        assert critical.size == 2
        assert float(np.max(np.abs(critical))) < 1e-8
        for frame in range(1, run["times"].size):
            counts = _counts(run["theta_plus"][frame], run["theta_minus"][frame])
            assert counts["marginal"] == 0
            assert counts["trapped"] > 0 and counts["anti_trapped"] > 0
            maximum = int(np.argmax(run["radius"][frame]))
            minimum = int(np.argmin(run["radius"][frame]))
            trapped_groups = _circular_groups(counts["trapped_mask"])
            anti_groups = _circular_groups(counts["anti_mask"])
            assert len(trapped_groups) == 1 and len(anti_groups) == 1
            assert maximum in set(int(index) for index in trapped_groups[0])
            assert minimum in set(int(index) for index in anti_groups[0])
            assert float(run["theta_plus"][frame, maximum]) < 0.0
            assert float(run["theta_minus"][frame, maximum]) < 0.0
            assert float(run["theta_plus"][frame, minimum]) > 0.0
            assert float(run["theta_minus"][frame, minimum]) > 0.0
    assert product_gap < 1e-15
    assert proper_gap < 1e-15
    assert norm_gap < 1e-14

    fine = primary
    half = saved["runs"]["nf512_dt_0_00025"]
    coarse = saved["runs"]["nf256_dt_0_0005"]
    coarse_half = saved["runs"]["nf256_dt_0_00025"]
    dt_gap = max(
        float(np.max(np.abs(fine["theta_plus"] - half["theta_plus"]))),
        float(np.max(np.abs(coarse["theta_plus"] - coarse_half["theta_plus"]))),
    )
    space_gap = float(np.max(np.abs(fine["theta_plus"][:, ::2] - coarse["theta_plus"])))
    initial_floor = max(
        2.0 * float(np.max(np.abs(run["proper"][0]))) / float(np.min(run["radius"][0]))
        for run in saved["runs"].values()
    )
    margin = 10.0 * max(dt_gap, space_gap, initial_floor)
    assert abs(margin - record["margin"]["absolute_expansion_margin"]) < 1e-18
    same = (fine["theta_plus"] * fine["theta_minus"]) > 0.0
    smallest = float(np.min(np.minimum(
        np.abs(fine["theta_plus"][same]), np.abs(fine["theta_minus"][same])
    )))
    assert smallest > margin
    assert dt_gap < 1e-8
    assert space_gap < 1e-6
    for run, other in ((fine, half), (coarse, coarse_half)):
        for frame in range(run["times"].size):
            left = _counts(run["theta_plus"][frame], run["theta_minus"][frame])
            right = _counts(other["theta_plus"][frame], other["theta_minus"][frame])
            assert left["trapped"] == right["trapped"]
            assert left["anti_trapped"] == right["anti_trapped"]
    disagreements = 0
    for frame in range(fine["times"].size):
        fine_counts = _counts(fine["theta_plus"][frame], fine["theta_minus"][frame])
        coarse_counts = _counts(coarse["theta_plus"][frame], coarse["theta_minus"][frame])
        fine_label = np.where(fine_counts["trapped_mask"], -1, np.where(fine_counts["anti_mask"], 1, 2))
        coarse_label = np.where(coarse_counts["trapped_mask"], -1, np.where(coarse_counts["anti_mask"], 1, 2))
        disagreements += int(np.sum(fine_label[::2] != coarse_label))
    assert disagreements == 0
    final_counts = _counts(fine["theta_plus"][-1], fine["theta_minus"][-1])
    assert final_counts["trapped"] == 45
    assert final_counts["anti_trapped"] == 182
    assert final_counts["untrapped"] == 1821
    coarse_final = _counts(coarse["theta_plus"][-1], coarse["theta_minus"][-1])
    assert coarse_final["trapped"] == 23
    assert coarse_final["anti_trapped"] == 91
    assert coarse_final["untrapped"] == 910
    product = fine["theta_plus"][-1] * fine["theta_minus"][-1]
    crossings = 0
    for index in range(product.size):
        neighbor = (index + 1) % product.size
        if product[index] == 0.0 or product[index] * product[neighbor] < 0.0:
            crossings += 1
    assert crossings == 4
    assert not np.any(fine["theta_plus"] == 0.0) and not np.any(fine["theta_minus"] == 0.0)
    # The quadrature node at the initial areal maximum is still opposite-sign.
    # The marginal object is the interpolated zero of r_x, inside the normal-velocity noise.
    initial_maximum = int(np.argmax(fine["radius"][0]))
    assert fine["theta_plus"][0, initial_maximum] * fine["theta_minus"][0, initial_maximum] < 0.0
    critical = _critical_thetas(fine["radius"][0], fine["proper"][0], fine["r_x"][0])
    assert float(np.max(np.abs(critical))) < margin
    index_at_two = int(round(2.0 / (PERIOD / fine["radius"].shape[1])))
    assert fine["lapse_factor"][0] == pytest.approx(CALIBRATION["b0"])
    assert fine["shift"][0] == pytest.approx(CALIBRATION["beta0"])
    assert fine["lapse_factor"][index_at_two] == pytest.approx(CALIBRATION["b0"] * CALIBRATION["Omega"] ** 2)
    dense = periodic_derivative(fine["radius"].shape[1], PERIOD)
    dense_slope = np.stack([dense @ row for row in fine["radius"]])
    derivative_gap = float(np.max(np.abs(dense_slope - fine["r_x"])))
    # The explicit derivative matrix and the FFT disagree by dense-matrix roundoff.
    # That gap stays far below the sign margin, and it does not move the sample counts.
    assert derivative_gap < 1e-8
    assert derivative_gap < 0.01 * margin
    dense_r_t = fine["proper"] * (fine["radius"] * fine["lapse_factor"]) + fine["shift"] * dense_slope
    dense_plus, dense_minus, _dense_gradient, _dense_norms = _null_data(
        fine["radius"], dense_r_t, dense_slope, fine["lapse_factor"], fine["conformal"], fine["shift"]
    )
    for frame in range(fine["times"].size):
        fft_counts = _counts(fine["theta_plus"][frame], fine["theta_minus"][frame])
        dense_counts = _counts(dense_plus[frame], dense_minus[frame])
        assert fft_counts["trapped"] == dense_counts["trapped"]
        assert fft_counts["anti_trapped"] == dense_counts["anti_trapped"]
        assert fft_counts["untrapped"] == dense_counts["untrapped"]

    # Replacing the instantaneous rate by the 0.005 frame secant moves theta by more
    # than the smallest same-sign sample. Edge signs are not resolved by that secant.
    step = np.diff(fine["times"])
    secant = (fine["radius"][1:] - fine["radius"][:-1]) / step[:, None]
    averaged = 0.5 * (fine["r_t"][:-1] + fine["r_t"][1:])
    theta_from_secant = 2.0 / fine["radius"][:-1] * np.abs(secant - averaged) / (
        fine["radius"][:-1] * fine["lapse_factor"]
    )
    assert float(np.max(theta_from_secant)) > smallest


def test_saved_response_binds_source_observer_and_separates_gaps(saved):
    response = saved["response"]
    arrays = saved["arrays"]
    spinors = saved["spinors"]
    assert response["stress_claimed"] is False
    assert response["regeneration"] is False
    assert "metric variation" in response["stress_gap"]
    assert response["headline"]["t_end"] == 0.035
    assert response["headline"]["t_start"] == 0.0
    phase = response["phases"]["full_window"]
    assert phase["complete"] is True
    assert phase["pilot_window_unchanged"] == [0.0, 0.035]
    assert phase["t_start"] == 0.0 and phase["t_end"] == 0.05
    assert phase["steps"] == 10 and phase["substeps"] == 2
    assert phase["declared_geometry"] == "piecewise-linear-stored-Q"
    assert phase["field_cross_active"] is False
    assert "memory" in phase["omissions"]["active"]
    assert "outside_drive" in phase["omissions"]["active"]
    assert "cross" not in phase["omissions"]["active"]

    source = source_id(spinors["phi0"], spinors["phi1"], spinors["weights"])
    assert response["requests"]["stored-frames"]["source"] == source
    assert response["requests"]["pilot"]["source"] == source
    assert response["requests"]["stored-frames"]["case"] == "nf512_dt_0_0005"
    assert response["requests"]["stored-frames"]["substeps"] == 2
    assert response["requests"]["stored-frames"]["memory"] is True
    assert response["requests"]["stored-frames"]["outside_drive"] is True
    assert response["requests"]["stored-frames"]["drop_cross_covariance"] is False
    assert response["bindings"]["episode_npz_sha256"] == saved["before"]["episode"]
    assert response["bindings"]["episode_verdict_used_as_acceptance"] is False
    assert response["bindings"]["geometry_rate_series_stored"] is False
    assert resume_disposition(
        response, request_identity("stored-frames", "nf512_dt_0_0005", 2, source)
    ) == "return"
    assert resume_disposition(
        response, request_identity("stored-frames", "nf512_dt_0_0005", 1, source)
    ) == "reject"
    assert resume_disposition(
        response, request_identity("pilot", "nf512_dt_0_0005", 2, source)
    ) == "reject"

    observer, columns, gram = region_observer(spinors["phi0"], spinors["phi1"], region=0)
    assert np.array_equal(observer, columns[:, :2])
    assert np.max(np.abs(gram - np.eye(6))) < 1e-15
    np.testing.assert_allclose(spinors["weights"], [0.75, 0.75, 0.5, 0.5, 0.25, 0.25], atol=0.0)
    initial = spinor_local_moments(observer, spinors["phi0"], spinors["phi1"], spinors["weights"])
    final = spinor_local_moments(
        observer, spinors["final_phi0"], spinors["final_phi1"], spinors["weights"]
    )
    np.testing.assert_allclose(initial["occupation"], spinors["weights"][:2], atol=1e-12)
    np.testing.assert_allclose(final["occupation"], phase["saved_autonomous"]["final_occupation"], atol=0.0)
    frame = fixed_observer_frame(observer)["frame"]
    assert np.array_equal(frame[:, :2], observer)
    coefficients = frame.conj().T @ columns
    weighted = coefficients * spinors["weights"]
    cross = weighted[:2] @ coefficients[2:].conj().T
    cross_norm = float(np.linalg.norm(cross, ord="fro"))
    assert cross_norm < 1e-14
    assert abs(cross_norm - phase["field_cross_norm"]) < 1e-18
    parent = weighted[:2] @ coefficients[:2].conj().T
    np.testing.assert_allclose(np.real(np.diag(parent)), spinors["weights"][:2], atol=1e-12)

    full = arrays["full_occupation_full"]
    retained = arrays["full_occupation_retained"]
    assert arrays["full_times"].shape == (11,)
    assert arrays["times"].shape == (8,)
    assert float(arrays["times"][-1]) == pytest.approx(0.035)
    assert float(arrays["full_times"][-1]) == pytest.approx(0.05)
    reduction_gap = float(np.max(np.abs(retained - full)))
    occupation_change = float(np.max(np.abs(full - full[0])))
    endpoint_gap = float(np.max(np.abs(full[-1] - np.asarray(final["occupation"]))))
    memory_gap = float(np.max(np.abs(arrays["full_occupation_memory_off"] - full)))
    drive_gap = float(np.max(np.abs(arrays["full_occupation_outside_drive_off"] - full)))
    assert reduction_gap == pytest.approx(7.353273811995242e-4)
    assert occupation_change == pytest.approx(0.3305919037025379)
    assert endpoint_gap == pytest.approx(3.4123144243558556e-7)
    assert memory_gap == pytest.approx(0.36709975607817125)
    assert drive_gap == pytest.approx(0.028764901996227665)
    assert memory_gap > occupation_change > drive_gap > reduction_gap > endpoint_gap
    coherence_gap = float(np.max(np.abs(
        arrays["full_coherence_retained"] - arrays["full_coherence_full"]
    )))
    assert coherence_gap == pytest.approx(5.5090520756831094e-6)
    assert coherence_gap < float(np.max(np.abs(arrays["full_coherence_full"])))

    frame_q = saved["runs"]["nf512_dt_0_0005"]["conformal"]
    frame_t = saved["runs"]["nf512_dt_0_0005"]["times"]
    assert spinors["final_q_dot"].shape == (frame_q.shape[1] // 4 - 1,)
    assert not any("frame" in name and "phi" in name for name in spinors["files"])
    assert spinors["replay_time"] == pytest.approx(0.0005)
    assert spinors["replay_time"] not in set(np.round(frame_t, decimals=12))
    variants = phase["interpolation_variants"]
    assert variants["saved_nodal_Q_rate_series"] is False
    midpoints = 0.5 * (frame_t[:-1] + frame_t[1:])
    linear = interpolate_rows(frame_t, frame_q, midpoints)
    curved = hermite_rows(frame_t, frame_q, midpoints)
    midpoint_gap = float(np.max(np.abs(linear - curved)))
    assert midpoint_gap == pytest.approx(variants["midpoint_max_abs_linear_versus_hermite"])
    assert midpoint_gap == pytest.approx(8.693988543928555e-5)
    assert float(np.max(np.abs(hermite_rows(frame_t, frame_q, frame_t) - frame_q))) == 0.0
    # The recorded Hermite spinor endpoint was not re-integrated in this test.
    hermite_occupation = phase["hermite_indicator"]["occupation_error_against_declared_full"]["max_abs"]
    assert hermite_occupation < reduction_gap
    assert hermite_occupation < drive_gap

    problem = prescribed_six_mode_problem()
    six_frame = fixed_observer_frame(problem["local_basis"])["frame"]
    transformed = six_frame.conj().T @ problem["covariance"] @ six_frame
    six_cross = float(np.linalg.norm(transformed[:2, 2:], ord="fro"))
    fresh = evolve_retained_region(
        problem["hamiltonian"],
        problem["local_basis"],
        np.linspace(0.0, 1.0, 17),
        covariance=problem["covariance"],
        backend="streamed",
    )
    dropped = evolve_retained_region(
        problem["hamiltonian"],
        problem["local_basis"],
        np.linspace(0.0, 1.0, 17),
        covariance=problem["covariance"],
        drop_cross_covariance=True,
        backend="streamed",
    )
    separation = float(np.linalg.norm(fresh["covariance_total"] - dropped["covariance_total"]))
    recorded = phase["correlated_cross_control"]
    assert recorded["control"] == "prescribed_six_mode_problem"
    assert recorded["active"] is True
    assert six_cross == pytest.approx(recorded["cross_norm"])
    assert separation == pytest.approx(recorded["covariance_separation"])
    assert six_cross > 1e-3 and separation > 1e-3
    assert six_cross > 1e8 * cross_norm


def _hermitian_pulse(time, scale=1.0):
    t = float(time)
    local = np.array(
        [[0.35 + 0.2 * np.sin(t), 0.12 - 0.18j], [0.12 + 0.18j, -0.22 + 0.1 * np.cos(t)]],
        dtype=np.complex128,
    )
    exterior = np.array(
        [[0.5 + 0.15 * np.cos(1.3 * t), 0.08j], [-0.08j, -0.16 + 0.05 * np.sin(t)]],
        dtype=np.complex128,
    )
    coupling = scale * np.array(
        [[0.24 + 0.05j, 0.11], [0.04 - 0.07j, -0.18 + 0.06j]],
        dtype=np.complex128,
    ) * (1.0 + 0.25 * np.sin(0.8 * t))
    matrix = np.zeros((4, 4), dtype=np.complex128)
    matrix[:2, :2] = 0.5 * (local + local.conj().T)
    matrix[2:, 2:] = 0.5 * (exterior + exterior.conj().T)
    matrix[:2, 2:] = coupling
    matrix[2:, :2] = coupling.conj().T
    return matrix


def _integrate_schrodinger(hamiltonian, times, columns):
    initial = np.asarray(columns, dtype=np.complex128)
    dimension, width = initial.shape

    def pack(state):
        return np.concatenate((state.real.ravel(), state.imag.ravel()))

    def unpack(flat):
        half = dimension * width
        return (flat[:half] + 1j * flat[half:]).reshape(dimension, width)

    def derivative(_time, flat):
        return pack(-1j * (hamiltonian(_time) @ unpack(flat)))

    solution = solve_ivp(
        derivative,
        (float(times[0]), float(times[-1])),
        pack(initial),
        t_eval=np.asarray(times, dtype=float),
        method="DOP853",
        rtol=1e-8,
        atol=1e-10,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return np.stack([unpack(solution.y[:, index]) for index in range(solution.y.shape[1])])


def _midpoint_columns(hamiltonian, times, columns):
    state = np.array(columns, dtype=np.complex128, copy=True)
    history = np.empty((times.size,) + state.shape, dtype=np.complex128)
    history[0] = state
    for index in range(times.size - 1):
        width = float(times[index + 1] - times[index])
        state = expm(-1j * width * hamiltonian(times[index] + 0.5 * width)) @ state
        history[index + 1] = state
    return history


def _moments(projected, weights):
    scale = np.asarray(weights, dtype=float)
    occupation = np.sum((np.abs(projected) ** 2) * scale, axis=-1)
    coherence = np.sum(projected[:, 0, :] * np.conjugate(projected[:, 1, :]) * scale, axis=-1)
    return occupation, coherence


def test_small_ode_kernel_is_causal_and_covariance_is_not_a_stress():
    basis = np.zeros((4, 2), dtype=np.complex128)
    basis[0, 0] = np.exp(0.4j)
    basis[1, 1] = np.exp(-0.9j)
    completed = fixed_observer_frame(basis)
    assert np.array_equal(completed["frame"][:, :2], basis)
    generator = np.random.default_rng(4)
    draw = generator.normal(size=(4, 3)) + 1j * generator.normal(size=(4, 3))
    columns, _triangle = np.linalg.qr(draw)
    weights = np.array([0.8, 0.45, 0.2])
    errors = []
    for nodes in (5, 9):
        times = np.linspace(0.0, 0.4, nodes)
        reduced = evolve_retained_region(
            _hermitian_pulse,
            basis,
            times,
            columns,
            weights=weights,
            backend="streamed",
        )
        reference = _integrate_schrodinger(_hermitian_pulse, times, columns)
        projected = np.einsum("ij,tjk->tik", basis.conj().T, reference)
        reduced_occupation, reduced_coherence = _moments(reduced["amplitudes"], weights)
        full_occupation, full_coherence = _moments(projected, weights)
        errors.append(float(np.max(np.abs(reduced_occupation - full_occupation))))
        assert np.array_equal(reduced["local_basis"], basis)
        assert int(reduced["allocation"]["time_indexed_exterior_propagator_bytes"]) == 0
        assert reduced["history"].shape[-1] == 3
        assert reduced["history"].shape[-2] == 2
        np.testing.assert_allclose(reduced_occupation[0], full_occupation[0], atol=1e-12)
        np.testing.assert_allclose(reduced_coherence[0], full_coherence[0], atol=1e-12)
        midpoint = _midpoint_columns(_hermitian_pulse, times, columns)
        midpoint_occupation, _midpoint_coherence = _moments(
            np.einsum("ij,tjk->tik", basis.conj().T, midpoint),
            weights,
        )
        assert float(np.max(np.abs(midpoint_occupation - full_occupation))) < errors[-1]
    assert errors[1] < errors[0]
    assert errors[1] < 1e-4

    times = np.linspace(0.0, 0.4, 9)
    baseline = evolve_retained_region(
        _hermitian_pulse, basis, times, columns, weights=weights, backend="streamed"
    )
    reference = _integrate_schrodinger(_hermitian_pulse, times, columns)
    full_occupation, _full_coherence = _moments(
        np.einsum("ij,tjk->tik", basis.conj().T, reference), weights
    )
    reduction_gap = float(np.max(np.abs(_moments(baseline["amplitudes"], weights)[0] - full_occupation)))
    for name, flags in (
        ("memory", {"memory": False}),
        ("drive", {"outside_drive": False}),
        ("coupling", {"coupling": False}),
    ):
        omitted = evolve_retained_region(
            _hermitian_pulse,
            basis,
            times,
            columns,
            weights=weights,
            backend="streamed",
            **flags,
        )
        occupation, _coherence = _moments(omitted["amplitudes"], weights)
        assert float(np.max(np.abs(occupation - full_occupation))) > 10.0 * reduction_gap
        if name == "memory":
            assert np.linalg.norm(omitted["history"]) == 0.0
        if name == "drive":
            assert np.isclose(
                np.linalg.norm(omitted["initial_exterior"]),
                np.linalg.norm(baseline["initial_exterior"]),
            )
            assert np.linalg.norm(baseline["initial_exterior"]) > 0.0

    def future_changed(time):
        matrix = _hermitian_pulse(time)
        if float(time) > 0.35:
            matrix = matrix + 0.4 * np.eye(4)
        return matrix

    coarse_times = np.linspace(0.0, 0.4, 3)
    original = evolve_retained_region(
        _hermitian_pulse, basis, coarse_times, columns, weights=weights, backend="streamed"
    )
    changed = evolve_retained_region(
        future_changed, basis, coarse_times, columns, weights=weights, backend="streamed"
    )
    assert np.array_equal(original["amplitudes"][1], changed["amplitudes"][1])
    assert np.linalg.norm(original["amplitudes"][2] - changed["amplitudes"][2]) > 1e-6

    dense = evolve_retained_region(
        _hermitian_pulse,
        basis,
        times,
        columns,
        weights=weights,
        covariance=None,
        backend="dense",
    )
    later, earlier = times.size - 1, 1
    future = exterior_kernels(
        float(times[later]),
        float(times[earlier]),
        coupling_t=dense["coupling_block"][later],
        propagator_t=dense["exterior_propagator"][later],
        coupling_s=dense["coupling_block"][earlier],
        propagator_s=dense["exterior_propagator"][earlier],
        exterior_covariance=np.eye(2),
        initial_cross=np.zeros((2, 2)),
    )
    past = exterior_kernels(
        float(times[earlier]),
        float(times[later]),
        coupling_t=dense["coupling_block"][earlier],
        propagator_t=dense["exterior_propagator"][earlier],
        coupling_s=dense["coupling_block"][later],
        propagator_s=dense["exterior_propagator"][later],
        exterior_covariance=np.eye(2),
        initial_cross=np.zeros((2, 2)),
    )
    equal = exterior_kernels(
        float(times[later]),
        float(times[later]),
        coupling_t=dense["coupling_block"][later],
        propagator_t=dense["exterior_propagator"][later],
        coupling_s=dense["coupling_block"][later],
        propagator_s=dense["exterior_propagator"][later],
        exterior_covariance=np.eye(2),
        initial_cross=np.zeros((2, 2)),
    )
    source_product = dense["coupling_block"][later] @ dense["exterior_propagator"][later]
    source_product = source_product @ (
        dense["exterior_propagator"][earlier].conj().T @ dense["coupling_block"][earlier].conj().T
    )
    assert np.linalg.norm(past["retarded"]) == 0.0
    assert np.linalg.norm(future["retarded"]) > 1e-6
    np.testing.assert_allclose(future["retarded"], -1j * source_product, atol=1e-12)
    np.testing.assert_allclose(equal["retarded"], -0.5j * (
        dense["coupling_block"][later]
        @ dense["exterior_propagator"][later]
        @ dense["exterior_propagator"][later].conj().T
        @ dense["coupling_block"][later].conj().T
    ), atol=1e-12)
    integrand = causal_memory_kernel(
        _hermitian_pulse, dense["frame"], 2, times, earlier, later, substeps=1
    )
    assert np.linalg.norm(integrand["kernel"]) > 0.0

    draw = generator.normal(size=(4, 4)) + 1j * generator.normal(size=(4, 4))
    covariance = draw @ draw.conj().T
    covariance = 0.5 * (covariance + covariance.conj().T)
    covariance /= np.trace(covariance)
    kept = evolve_retained_region(
        _hermitian_pulse, basis, times, covariance=covariance, backend="streamed"
    )
    without_cross = evolve_retained_region(
        _hermitian_pulse,
        basis,
        times,
        covariance=covariance,
        drop_cross_covariance=True,
        backend="streamed",
    )
    assert np.linalg.norm(kept["covariance_total"] - without_cross["covariance_total"]) > 1e-4
    assert np.allclose(kept["retained_map"], without_cross["retained_map"], rtol=0.0, atol=1e-12)
    assert kept["amplitudes"] is None and without_cross["amplitudes"] is None
    assert "metric" not in evolve_retained_region.__code__.co_varnames


_V1_JSON = _LAB / "results" / "development" / "nsc-regeneration-controls-v1.json"
_V1_NPZ = _LAB / "results" / "development" / "nsc-regeneration-controls-v1.npz"
_V2_JSON = _LAB / "results" / "development" / "nsc-regeneration-controls-v2.json"
_CONTROLS_NOTE = _LAB / "docs" / "nsc-regeneration-controls.md"
_V1_JSON_SHA = "3abe72de1be7f5090569086a3134664727e25849eecada6f12b89040bdddc922"
_V1_NPZ_SHA = "8509a1d0a8151c189ce3346d7d47f7c4a6453065d39ec01fda35515af037c96b"


def test_control_policy_keeps_g_initial_and_proxies_off_the_gate():
    """G > 0 admits the convex initial bracket. The dynamic chart is r, Q, L."""
    from recursive_horizons.nsc_regeneration_controls import (
        CAR_GAP_MAX,
        CRITERION,
        PERTURBATIONS_FROM_PROXY,
        bracket_g,
        chart_from_fields,
        dynamic_stop_reason,
    )
    from recursive_horizons.nsc_spherical_coupling import chart_failure
    from recursive_horizons.nsc_spherical_galerkin_coupling import build_grid

    assert "gee" not in dynamic_stop_reason.__code__.co_varnames
    assert dynamic_stop_reason(None, {"admissible": True}) == (None, False)
    assert dynamic_stop_reason("Q_left_positive_chart", {"admissible": True}) == (
        "Q_left_positive_chart",
        True,
    )
    assert dynamic_stop_reason(None, {"admissible": False}) == ("GRAM_LEFT_IDENTITY", False)
    assert PERTURBATIONS_FROM_PROXY is False
    assert "not a dynamic stop" in CRITERION["chart_stop"]
    assert "not programme requirements" in CRITERION["proxies"]
    assert "does not by itself open the fine budget" in CRITERION["fine_only_if"]

    grid = build_grid(14, quadrature=64)
    radius = np.full(grid.nq, 2.0)
    conformal = np.full(grid.nq, 0.25)
    rho = np.full(grid.nq, -80.0)
    gee = bracket_g(conformal, rho, grid)
    chart = chart_from_fields(grid, {"Q": conformal, "rho": rho, "r": radius})
    assert float(np.min(gee)) < 0.0
    assert chart["positive_G"] is False
    assert chart["g_is_dynamic_stop"] is False
    assert chart["positive_r"] and chart["positive_Q"] and chart["positive_L"] and chart["positive_lapse"]
    assert chart_failure(grid.fine, type("S", (), {
        "Q": conformal, "r": radius, "chi": np.zeros(grid.nq),
        "p_Q": np.zeros(grid.nq), "p_r": np.zeros(grid.nq), "p_chi": np.zeros(grid.nq),
        "phi0": np.zeros((grid.nq, 1)), "phi1": np.zeros((grid.nq, 1)),
    })()) is None
    broken = conformal.copy()
    broken[0] = -0.1
    assert chart_failure(grid.fine, type("S", (), {
        "Q": broken, "r": radius, "chi": np.zeros(grid.nq),
        "p_Q": np.zeros(grid.nq), "p_r": np.zeros(grid.nq), "p_chi": np.zeros(grid.nq),
        "phi0": np.zeros((grid.nq, 1)), "phi1": np.zeros((grid.nq, 1)),
    })()) == "Q_left_positive_chart"
    assert CAR_GAP_MAX == 1e-8


def test_v2_reads_stored_endpoint_without_touching_v1():
    """Stability and Gram come from the stored T=0.10 state. No RK4 step is taken."""
    from recursive_horizons.nsc_conformal_adm_source import direct_hamiltonian
    from recursive_horizons.nsc_covariant_operator import CovariantStaticMetric
    from recursive_horizons.nsc_spherical_coupling import CauchyState, apply_dirac, chart_failure
    from recursive_horizons.nsc_spherical_galerkin_coupling import (
        build_grid,
        load_physical_columns,
        prolong_state,
        stable_timestep,
        subspace_frequency,
    )
    from recursive_horizons.nsc_spherical_galerkin_coupling import blank_state

    before_json = _sha256(_V1_JSON)
    before_npz = _sha256(_V1_NPZ)
    note = _CONTROLS_NOTE.read_text()
    assert before_json == _V1_JSON_SHA and before_json in note
    assert before_npz == _V1_NPZ_SHA and before_npz in note
    v1 = json.loads(_V1_JSON.read_text())
    v2 = json.loads(_V2_JSON.read_text())
    assert v1["status"] == "MEASURED"
    assert v1["cpu_seconds"] < 600.0
    assert v1["conclusion"]["fine_cpu_seconds"] == 0.0
    assert v1["conclusion"]["perturbations_executed"] is False
    assert v1["conclusion"]["candidate_regime"] is False
    assert v1["conclusion"]["maintained_structure"] is True
    assert v1["conclusion"]["renewed_structure"] is False
    fine = v1["stages"]["continuation_nf512"]
    assert fine["stop_reason"] is None and fine["chart_stop"] is False
    assert fine["dt"] == 0.0005 and fine["steps_completed"] == 100
    assert fine["attained_time"] == 0.1
    assert "positive_L" not in fine["checkpoints"][-1]
    assert "gram_gap" not in fine["checkpoints"][-1]
    assert v2["dynamics_rerun"] is False and v2["new_action"] is False and v2["new_forces"] is False
    assert v2["predecessor"]["json_sha256"] == before_json
    assert v2["predecessor"]["npz_sha256"] == before_npz
    assert v2["predecessor"]["immutable"] is True
    assert v2["next_event"]["executed"] is False
    assert v2["policy"]["g_is_dynamic_chart_stop"] is False
    assert v2["policy"]["candidate_regime_is_programme_requirement"] is False
    assert v2["ledgers"]["proxies"]["programme_requirements"] is False
    assert v2["ledgers"]["proxies"]["candidate_regime"] is False
    assert v2["ledgers"]["maintained_with_throughflow"] is True
    assert v2["ledgers"]["reversal"]["renewed"] is False

    with np.load(_V1_NPZ, allow_pickle=False) as stored:
        final = CauchyState(
            Q=np.array(stored["continuation_nf512_final_Q"], dtype=float),
            r=np.array(stored["continuation_nf512_final_r"], dtype=float),
            chi=np.array(stored["continuation_nf512_final_chi"], dtype=float),
            p_Q=np.array(stored["continuation_nf512_final_p_Q"], dtype=float),
            p_r=np.array(stored["continuation_nf512_final_p_r"], dtype=float),
            p_chi=np.array(stored["continuation_nf512_final_p_chi"], dtype=float),
            phi0=np.array(stored["continuation_nf512_final_phi0"]),
            phi1=np.array(stored["continuation_nf512_final_phi1"]),
        )
        unitarity = np.array(stored["continuation_nf512_unitarity"], dtype=float)
        q_min = np.array(stored["continuation_nf512_Q_min"], dtype=float)
        times = np.array(stored["continuation_nf512_time"], dtype=float)
        content = np.array(stored["continuation_nf512_content"], dtype=float)
        flux = np.array(stored["continuation_nf512_flux"], dtype=float)
    gram = final.phi0.conj().T @ final.phi0 + final.phi1.conj().T @ final.phi1
    gram_gap = float(np.max(np.abs(gram - np.eye(6))))
    assert gram_gap == pytest.approx(1.3069428872469757e-9, rel=0.0, abs=0.0)
    assert gram_gap <= 1e-8
    assert float(np.max(unitarity)) == pytest.approx(1.3069464399606545e-9, rel=0.0, abs=0.0)
    assert abs(gram_gap - float(np.max(unitarity))) < 1e-14

    grid = build_grid(512)
    omega, quadrature_omega = subspace_frequency(grid, final)
    step, omega_again, _quadrature_again = stable_timestep(grid, final)
    prolonged = prolong_state(grid, final)
    assert chart_failure(grid.fine, prolonged) is None
    assert float(np.min(prolonged.r)) > 0.0
    assert float(np.min(prolonged.Q)) > 0.0
    assert float(np.min(grid.fine.length_density)) > 0.0
    assert omega == pytest.approx(2998.824498749641, rel=0.0, abs=0.0)
    assert omega_again == omega
    assert quadrature_omega == pytest.approx(804.247719318987, rel=0.0, abs=0.0)
    dt_omega = 0.0005 * omega
    assert dt_omega == pytest.approx(1.4994122493748205, rel=0.0, abs=0.0)
    assert 1.4 < dt_omega < 2.0 * np.sqrt(2.0)
    assert step == pytest.approx(min(5.0e-4, 1.4 / omega), rel=0.0, abs=0.0)
    assert step == pytest.approx(4.668495940938623e-4, rel=0.0, abs=0.0)
    lapse_peak = float(np.max(grid.fine.length_density))
    shift_peak = float(np.max(np.abs(grid.fine.shift)))
    wavenumber = float(np.max(np.abs(grid.modes_f)) * 2.0 * np.pi / grid.length)
    mass_peak = float(np.max(np.abs(grid.fine.kappa * grid.fine.length_density)))
    diagnostic = (lapse_peak / q_min + shift_peak) * wavenumber + mass_peak
    assert abs(diagnostic[-1] - omega) / omega == pytest.approx(6.122469423745533e-4, rel=1e-12, abs=0.0)
    crossed = np.flatnonzero(0.0005 * diagnostic >= 1.4)
    assert int(crossed[0]) == 93
    assert float(times[93]) == pytest.approx(0.0965)
    assert float(q_min[93]) == pytest.approx(0.1345054206385445, rel=0.0, abs=0.0)
    assert float(content[93, 0] - content[0, 0]) == pytest.approx(0.10855338655846847, rel=0.0, abs=0.0)
    shares = content / np.sum(content, axis=1)[:, None]
    packet = shares[-1] >= 0.30
    assert float(np.sum(flux[-1, packet])) == pytest.approx(-5.534007341562886, rel=0.0, abs=1e-12)
    assert abs(float(np.sum(flux[-1]))) < 1e-12
    assert float(content[-1, 0] - content[0, 0]) == pytest.approx(0.11695656049787839, rel=0.0, abs=0.0)
    assert v2["stability"]["above_owned_cap"] is True
    assert v2["stability"]["below_absolute_rk4_limit"] is True
    assert v2["stability"]["gram_admissible"] is True
    assert v2["stability"]["chart_failure"] is None
    assert v2["stability"]["positive_L"] is True

    small = build_grid(14, quadrature=64)
    phi0, phi1, _preparation, problems = load_physical_columns(14)
    assert problems == []
    fine_state = prolong_state(small, blank_state(small, phi0, phi1))
    system = small.fine
    metric = CovariantStaticMetric(
        system.length,
        system.length_density * fine_state.r,
        fine_state.Q * fine_state.r,
        fine_state.r,
        eta=0.5,
    )
    hamiltonian = direct_hamiltonian(metric, system.shift, system.kappa)
    embedding = np.zeros((2 * system.points, 28), dtype=np.complex128)
    embedding[: system.points, :14] = small.U_f
    embedding[system.points :, 14:] = small.U_f
    projected = embedding.conj().T @ hamiltonian @ embedding
    half = 0.5 * np.eye(28, dtype=np.complex128)
    assert np.linalg.norm(projected @ half - half @ projected) == 0.0
    random = np.random.default_rng(0).normal(size=(28, 28))
    assert np.linalg.norm(random @ half - half @ random) == 0.0
    embedded = embedding @ half @ embedding.conj().T
    embedding_gap = float(np.linalg.norm(hamiltonian @ embedded - embedded @ hamiltonian))
    assert embedding_gap == pytest.approx(11.084247729335008, rel=1e-12, abs=0.0)
    assert embedding_gap > 1.0
    eye = np.eye(14, dtype=np.complex128)
    zero = np.zeros((14, 14), dtype=np.complex128)
    image0, image1 = apply_dirac(
        small.U_f @ np.concatenate((eye, zero), axis=1),
        small.U_f @ np.concatenate((zero, eye), axis=1),
        system.length_density,
        fine_state.Q,
        system.shift,
        system.kappa,
        system.momentum,
    )
    pulled = np.vstack((small.U_f.conj().T @ image0, small.U_f.conj().T @ image1))
    assert float(np.max(np.abs(pulled - projected))) == pytest.approx(7.105427357601002e-15, rel=0.0, abs=0.0)
    assert _sha256(_V1_JSON) == before_json
    assert _sha256(_V1_NPZ) == before_npz
