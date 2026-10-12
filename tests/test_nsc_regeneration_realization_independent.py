"""Independent realization audit. Reads stored payloads. Does not evolve a state.

The Gram gap is not the covariance. C = Phi W Phi†, with W the six Gaussian
weights. A passing label in a JSON file is not evidence.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_spherical_cauchy_weak import (
    saved_lapse_existence,
    saved_shift_constraint_bound,
)
from recursive_horizons.nsc_spherical_galerkin_coupling import (
    build_grid,
    prolong_geometry,
    pull_geometry,
)
from recursive_horizons.nsc_spherical_null_expansion import (
    spectral_derivative,
    static_clock,
    summarize_frame,
)

LAB = Path(__file__).resolve().parents[1]
DEV = LAB / "results" / "development"
EPISODE_NPZ = DEV / "nsc-regeneration-episode-v1.npz"
EPISODE_JSON = DEV / "nsc-regeneration-episode-v1.json"
PREDECESSOR_NPZ = DEV / "nsc-spherical-feedback-episode-v1.npz"
RESPONSE_NPZ = DEV / "nsc-coupled-local-response-v2.npz"
RESPONSE_JSON = DEV / "nsc-coupled-local-response-v2.json"
CONTROLS_NPZ = DEV / "nsc-regeneration-controls-v1.npz"
ERROR_JSON = DEV / "nsc-spherical-cauchy-error-v1.json"
NULL_JSON = DEV / "nsc-spherical-null-expansion-v1.json"
ASSESSMENT = LAB / "src" / "recursive_horizons" / "nsc_spherical_episode_assessment.py"

RUNS = (
    "nf256_dtmax_0_0005",
    "nf256_dtmax_0_00025",
    "nf512_dtmax_0_0005",
    "nf512_dtmax_0_00025",
)
HANDOFF = {
    "nf256_dtmax_0_0005": "nf256_dt_0_0005",
    "nf256_dtmax_0_00025": "nf256_dt_0_0005",
    "nf512_dtmax_0_0005": "nf512_dt_0_0005",
    "nf512_dtmax_0_00025": "nf512_dt_0_0005",
}
SIBLING = {
    "nf256_dtmax_0_0005": "nf256_dt_0_00025",
    "nf512_dtmax_0_0005": "nf512_dt_0_00025",
}
FIELDS = ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1")
WEIGHTS = np.array([0.75, 0.75, 0.5, 0.5, 0.25, 0.25])
PRIMARY = "nf512_dtmax_0_00025"
PACKET_END = 4.0


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _state_hash(arrays):
    digest = hashlib.sha256()
    for array in arrays:
        digest.update(np.ascontiguousarray(array).tobytes())
    return digest.hexdigest()


def _trap(times, values):
    dt = np.diff(times)
    midpoint = 0.5 * (values[:-1] + values[1:])
    scale = dt.reshape((-1,) + (1,) * (values.ndim - 1))
    return np.sum(midpoint * scale, axis=0)


def _frame(payload, run, name, index=0):
    if name in ("phi0", "phi1"):
        return payload[f"{run}_frame_{name}"][index]
    return payload[f"{run}_frame_coarse_{name}"][index]


def _covariance_spectrum(phi0, phi1, weights=WEIGHTS):
    """Nonzero eigenvalues of C = Phi W Phi†. They are the spectrum of W^{1/2} Gram W^{1/2}."""
    phi = np.vstack((phi0, phi1))
    gram = phi.conj().T @ phi
    scale = np.sqrt(weights)
    reduced = (scale[:, None] * gram) * scale[None, :]
    reduced = 0.5 * (reduced + reduced.conj().T)
    return gram, np.linalg.eigvalsh(reduced).real


def _cross_frobenius(observer, phi0, phi1, weights=WEIGHTS):
    """||V† C P_perp||_F without forming the dense exterior block."""
    phi = np.vstack((phi0, phi1))
    coefficients = observer.conj().T @ phi
    parent = (coefficients * weights) @ coefficients.conj().T
    parent = 0.5 * (parent + parent.conj().T)
    gram = phi.conj().T @ phi
    left = coefficients * weights
    product = float(np.real(np.trace(left @ gram @ left.conj().T)))
    gram_v = observer.conj().T @ observer
    parent_piece = float(np.real(np.trace(parent @ np.linalg.inv(gram_v) @ parent.conj().T)))
    return float(np.sqrt(max(0.0, product - parent_piece))), parent, gram_v


@pytest.fixture(scope="module")
def payload():
    return {
        "episode": np.load(EPISODE_NPZ),
        "episode_json": json.loads(EPISODE_JSON.read_text()),
        "predecessor": np.load(PREDECESSOR_NPZ),
        "response": np.load(RESPONSE_NPZ),
        "response_json": json.loads(RESPONSE_JSON.read_text()),
        "controls": np.load(CONTROLS_NPZ),
        "error": json.loads(ERROR_JSON.read_text()),
        "margin": float(json.loads(NULL_JSON.read_text())["margin"]["absolute_expansion_margin"]),
    }


def test_bitwise_handoff(payload):
    episode = payload["episode"]
    predecessor = payload["predecessor"]
    record = payload["episode_json"]
    assert _sha256(EPISODE_NPZ) == "17364df24efead5814c3ab2a2ab5f82d9e09bf412c8c82e68d57a3b38498d7a1"
    assert _sha256(PREDECESSOR_NPZ) == "2864d3a8d3ba413e840c395961f14a3d71eabee702c27fd9f47c7d333eb24b77"
    for run, origin in HANDOFF.items():
        arrays = [_frame(episode, run, name) for name in FIELDS]
        reference = [predecessor[f"{origin}_final_{name}"] for name in FIELDS]
        for got, ref in zip(arrays, reference):
            assert np.array_equal(got, ref)
        assert _state_hash(arrays) == record["runs"][run]["initial_hash"]
        assert abs(float(episode[f"{run}_time"][0]) - 0.05) < 1e-15
        assert abs(float(episode[f"{run}_frame_time"][0]) - 0.05) < 1e-15
        assert abs(float(episode[f"{run}_time"][-1]) - 0.085) < 1e-12
    assert record["runs"]["nf512_dtmax_0_00025"]["initial_hash"] == (
        "cd6b1f88a12635b137aaddc47444c06cc681c7568271ef989312f48f43a0649e"
    )
    for run, sibling in SIBLING.items():
        assert not np.array_equal(
            episode[f"{run}_frame_phi0"][0],
            predecessor[f"{sibling}_final_phi0"],
        )
    assert record["state_reset"] is False
    assert record["new_force_added"] is False


def test_covariance_car_not_the_gram_flag(payload):
    episode = payload["episode"]
    for run in RUNS:
        for phi0, phi1 in (
            (episode[f"{run}_frame_phi0"][0], episode[f"{run}_frame_phi1"][0]),
            (episode[f"{run}_final_phi0"], episode[f"{run}_final_phi1"]),
        ):
            gram, spectrum = _covariance_spectrum(phi0, phi1)
            gap = float(np.max(np.abs(gram - np.eye(6))))
            assert gap < 1.0e-8
            assert spectrum.shape == (6,)
            order = np.sort(spectrum)
            assert np.allclose(order, np.sort(WEIGHTS), rtol=0.0, atol=2.0e-9)
            assert float(order.min()) > 0.2
            assert float(order.max()) < 0.76
            trace = float(np.real(np.trace(gram * WEIGHTS)))
            assert abs(trace - 3.0) < 1.0e-8
    # The controls flag names the mode Gram. It is not this spectrum.
    assert abs(float(np.sort(_covariance_spectrum(
        episode[f"{PRIMARY}_frame_phi0"][0],
        episode[f"{PRIMARY}_frame_phi1"][0],
    )[1])[0]) - 0.25) < 2.0e-9


def test_chart_timestep_and_four_runs(payload):
    episode = payload["episode"]
    record = payload["episode_json"]
    effects = {}
    for run in RUNS:
        radius = episode[f"{run}_frame_quad_r"]
        conformal = episode[f"{run}_frame_quad_Q"]
        assert float(radius.min()) > 4.2
        assert float(conformal.min()) > 0.16
        assert bool(episode[f"{run}_positive_r"].all())
        assert bool(episode[f"{run}_positive_Q"].all())
        assert bool(episode[f"{run}_positive_L"].all())
        assert float(episode[f"{run}_G_min"].min()) > 0.02
        product = episode[f"{run}_dt_omega_values"]
        assert float(np.max(np.abs(episode[f"{run}_dt_values"] * episode[f"{run}_omega_values"] - product))) == 0.0
        assert float(product.max()) < 1.4
        assert float(product.max()) < 2.0 * np.sqrt(2.0)
        assert int(record["runs"][run]["rejected_count"]) == 0
        assert record["runs"][run]["stop_reason"] == "BRIDGE_ARC_ENTERED_PACKET"
        effects[run] = {
            "field": float(episode[f"{run}_field_energy"][-1] - episode[f"{run}_field_energy"][0]),
            "content": float(episode[f"{run}_leader_content"][-1] - episode[f"{run}_leader_content"][0]),
            "chi": float(episode[f"{run}_chi_max"][-1] - episode[f"{run}_chi_max"][0]),
            "proper": float(episode[f"{run}_proper_max"][-1] - episode[f"{run}_proper_max"][0]),
            "flux": float(episode[f"{run}_packet_flux"][-1]),
        }
    primary = effects[PRIMARY]
    assert primary["field"] < -0.2
    assert primary["content"] > 0.05
    assert primary["flux"] < -4.0
    for run, row in effects.items():
        for key, value in row.items():
            assert abs(value - primary[key]) <= 0.01 * max(abs(primary[key]), 1.0e-12)
    _, lapse, shift = static_clock(int(episode[f"{PRIMARY}_frame_quad_r"].shape[-1]), 8.0)
    assert float(np.min(lapse)) > 0.1
    assert np.isfinite(shift).all()


def test_window_ledger_and_projectors(payload):
    episode = payload["episode"]
    for run in RUNS:
        times = episode[f"{run}_time"]
        normal = episode[f"{run}_window_normal"]
        rate = (
            episode[f"{run}_observer_boundary"]
            + episode[f"{run}_observer_pressure"]
            + episode[f"{run}_observer_lapse"]
        )
        residual = (normal[-1] - normal[0]) - _trap(times, rate)
        field = float(episode[f"{run}_field_energy"][-1] - episode[f"{run}_field_energy"][0])
        gravity = float(episode[f"{run}_gravity_energy"][-1] - episode[f"{run}_gravity_energy"][0])
        assert float(np.max(np.abs(residual))) < 2.0e-7
        assert float(np.max(np.abs(residual))) < 1.0e-5 * abs(field)
        measured = episode[f"{run}_window_coordinate"]
        split = episode[f"{run}_window_matter"] + episode[f"{run}_window_gravity"]
        assert float(np.max(np.abs(measured - split))) < 1.0e-12
        assert float(np.max(np.abs(episode[f"{run}_window_proper_flux"].sum(axis=1)))) < 1.0e-10
        assert abs(float(episode[f"{run}_packet_flux"][-1] + episode[f"{run}_reservoir_flux"][-1])) < 1.0e-10
        work = float(episode[f"{run}_coordinate_work_integral"][-1] - episode[f"{run}_coordinate_work_integral"][0])
        pressure = float(
            episode[f"{run}_proper_pressure_work_integral"][-1]
            - episode[f"{run}_proper_pressure_work_integral"][0]
        )
        assert abs(work - field) < 1.0e-4 * abs(field)
        assert abs(pressure) > 1.0e-4
        assert abs(pressure) < 0.01 * abs(field)
        assert abs(field + gravity) < 1.0e-4 * abs(field)
        assert int(episode[f"{run}_leader"][0]) == 0
        assert int(episode[f"{run}_leader"][-1]) == 0
        assert float(episode[f"{run}_leader_share"][-1]) > float(episode[f"{run}_leader_share"][0]) > 0.6
    full = float(episode[f"{PRIMARY}_full_hamilton_max"].max())
    projected = float(episode[f"{PRIMARY}_projected_hamilton_max"].max())
    held = float(episode[f"{PRIMARY}_held_out_hamilton_max"].max())
    assert full > 1.0e-4
    assert projected > 0.5 * full
    assert 0.0 < held < full
    coarse = float(episode["nf256_dtmax_0_00025_full_hamilton_max"].max())
    assert coarse > 50.0 * full


def test_realized_jet_projection_and_proxy_normalization(payload):
    episode = payload["episode"]
    for run in RUNS:
        step = episode[f"{run}_frame_increment_dt"]
        increment = episode[f"{run}_frame_increment_Q_dot"]
        indicator = episode[f"{run}_frame_indicator_Q_dot"]
        assert bool(episode[f"{run}_frame_increment_present"].all())
        assert float(np.max(np.abs(increment / step[:, None] - indicator))) == 0.0
        qdot = episode[f"{run}_frame_Q_dot"]
        secant = (qdot[1:] - qdot[:-1]) / np.diff(episode[f"{run}_frame_time"])[:, None]
        assert float(np.max(np.abs(secant - indicator[:-1]))) > 1.0
        chi = episode[f"{run}_frame_quad_chi"]
        radius = episode[f"{run}_frame_quad_r"]
        proxy = chi ** 2 / (3.0 * radius ** 4)
        assert float(np.max(np.abs(episode[f"{run}_frame_quad_weyl_C2"] - proxy))) < 1.0e-14
    for fermions, run in ((256, "nf256_dtmax_0_00025"), (512, PRIMARY)):
        grid = build_grid(fermions)
        rate = episode[f"{run}_frame_Q_dot"]
        gap = max(
            float(np.max(np.abs(pull_geometry(grid, prolong_geometry(grid, row)) - row)))
            for row in rate
        )
        assert gap < 1.0e-12
    assert payload["episode_json"]["payload_schema"]["metric_curvature_computed_here"] is False


def test_arc_crossing_is_not_a_time_split(payload):
    episode = payload["episode"]
    record = payload["episode_json"]
    margin = payload["margin"]
    radius = episode[f"{PRIMARY}_frame_quad_r"]
    _, lapse, shift = static_clock(radius.shape[-1], 8.0)
    edges = []
    for index in (0, -1):
        sample = radius[index]
        frame = summarize_frame(
            float(episode[f"{PRIMARY}_frame_time"][index]),
            sample,
            episode[f"{PRIMARY}_frame_quad_Q"][index],
            episode[f"{PRIMARY}_frame_quad_proper"][index],
            episode[f"{PRIMARY}_frame_quad_K_perp"][index],
            spectral_derivative(sample, 8.0),
            lapse,
            shift,
            margin,
        )
        anti = frame["anti_trapped_intervals"]
        assert len(anti) == 1
        edges.append((float(anti[0]["x_start"]), float(anti[0]["x_end"])))
    assert edges[0][0] > PACKET_END
    assert edges[1][0] < PACKET_END < edges[1][1]
    assert abs(edges[1][0] - 3.96875) < 1.0e-12
    for run in RUNS:
        stored = record["runs"][run]["expansion_frames"]
        start = stored[0]["anti_trapped_intervals"][0]["x_start"]
        final = stored[-1]["anti_trapped_intervals"][0]
        assert start > PACKET_END
        assert final["x_start"] < PACKET_END < final["x_end"]
        assert abs(final["x_start"] - 3.96875) < 1.0e-12
    reversal = record["structures"][PRIMARY]["joined_ledgers"]["reversal"]
    assert reversal["renewed"] is False
    assert reversal["one_drift"] is True
    assert reversal["reversal"] == 0.0


def test_frozen_geometry_and_changed_sources(payload):
    controls = payload["controls"]
    predecessor = payload["predecessor"]
    coupled = float(
        predecessor["nf512_dt_0_0005_field_energy"][-1]
        - predecessor["nf512_dt_0_0005_field_energy"][0]
    )
    frozen_field = float(controls["frozen_nf512_field_energy"][-1] - controls["frozen_nf512_field_energy"][0])
    assert coupled < -0.4
    assert abs(frozen_field) < 1.0e-6
    assert float(controls["frozen_nf512_chi_max"][-1] - controls["frozen_nf512_chi_max"][0]) == 0.0
    assert float(controls["frozen_nf512_proper_max"][-1] - controls["frozen_nf512_proper_max"][0]) == 0.0
    assert float(controls["frozen_nf512_leader_content"][-1] - controls["frozen_nf512_leader_content"][0]) > 0.1
    assert float(controls["frozen_nf512_full_h"].max()) > 20.0
    assert float(np.max(np.abs(
        controls["frozen_nf512_final_phi0"] - predecessor["nf512_initial_phi0"]
    ))) > 0.1
    assert abs(float(controls["frozen_nf512_field_energy"][0] - predecessor["nf512_dt_0_0005_field_energy"][0])) < 1.0e-9
    for name, floor in (("uniform_nf512", -0.5), ("reversed_nf512", -0.6)):
        change = float(controls[f"{name}_field_energy"][-1] - controls[f"{name}_field_energy"][0])
        assert change < floor
        assert float(controls[f"{name}_full_h"][0]) < 1.0e-4
        assert float(controls[f"{name}_G_min"].min()) > 0.04
        assert float(controls[f"{name}_Q_min"].min()) > 0.2
        assert float(controls[f"{name}_r_min"].min()) > 4.0
    assert float(controls["uniform_nf256_full_h"][0]) > 1.0e-3
    assert float(controls["reversed_nf256_full_h"][0]) > 1.0e-3


def test_local_response_omissions_from_stored_series(payload):
    response = payload["response"]
    record = payload["response_json"]
    predecessor = payload["predecessor"]
    episode = payload["episode"]
    observer = np.vstack((
        predecessor["nf512_initial_phi0"],
        predecessor["nf512_initial_phi1"],
    ))[:, :2].copy()
    transported = np.vstack((
        episode[f"{PRIMARY}_frame_phi0"][0],
        episode[f"{PRIMARY}_frame_phi1"][0],
    ))
    assert not np.array_equal(observer, transported[:, :2])
    cross, parent, gram_v = _cross_frobenius(
        observer,
        episode[f"{PRIMARY}_frame_phi0"][0],
        episode[f"{PRIMARY}_frame_phi1"][0],
    )
    assert float(np.max(np.abs(gram_v - np.eye(2)))) < 1.0e-12
    assert cross == pytest.approx(0.4422387105894591, rel=0.0, abs=1.0e-12)
    parent_ev = np.sort(np.linalg.eigvalsh(parent).real)
    assert parent_ev[0] == pytest.approx(0.4193978304297541, rel=0.0, abs=1.0e-12)
    assert parent_ev[1] == pytest.approx(0.6474479155567984, rel=0.0, abs=1.0e-12)
    full = response["occupation_full"]
    autonomous = response["occupation_autonomous"]
    retained = response["occupation_retained"]
    assert float(np.max(np.abs(autonomous[-1] - autonomous[0]))) == pytest.approx(0.2674889273028506)
    assert float(np.max(np.abs(full - autonomous))) == pytest.approx(1.1666190363746054e-07)
    assert float(np.max(np.abs(retained - full))) == pytest.approx(0.0002570304093015008)
    memory = float(np.max(np.abs(response["occupation_memory_off"] - full)))
    drive = float(np.max(np.abs(response["occupation_outside_drive_off"] - full)))
    assert memory == pytest.approx(0.06952775462900784)
    assert drive == pytest.approx(0.15526034047924991)
    assert min(memory, drive, cross) > float(np.max(np.abs(retained - full)))
    names = set(response.files)
    assert "occupation_cross_off" not in names
    assert "probe_cross_occupation_movement" in names
    probe_cross = float(np.max(np.abs(response["probe_cross_occupation_movement"])))
    assert probe_cross == pytest.approx(0.08797142267734728)
    recorded_cross = record["phases"]["full_controls"]["omissions"]["initial_cross"]["occupation_diagonal_movement"]
    assert recorded_cross == pytest.approx(0.2043417632796783)
    assert record["stress_claimed"] is False
    assert record["force_claimed"] is False
    assert record["renewal"] is False
    assert np.allclose(predecessor["nf512_occupations"], WEIGHTS)


def test_initial_constraint_scope_is_not_an_evolution_bound(payload):
    error = payload["error"]
    predecessor = payload["predecessor"]
    episode = payload["episode"]
    assert error["evolution_error_bound"] is None
    assert error["minimizer_existence_proved"] is None
    assert error["total_certified"] is False
    assert error["mean_current_removed"] is False
    assert error["pointwise_radius_error_upper_bound"] == pytest.approx(4.060849694930819e-06)
    theorem = saved_lapse_existence()
    assert theorem["hypotheses_hold"] is True
    assert theorem["unique_positive_critical_point"] is True
    assert theorem["smooth_critical_point"] is True
    assert theorem["q_squared"] == pytest.approx(0.0631642220827373)
    assert theorem["g_lower"] == pytest.approx(0.053733032411766975)
    assert theorem["evolution_error_bound"] is None
    assert theorem["historical_field_left_null"] is True
    shift = saved_shift_constraint_bound()
    assert shift["mean_contains_zero"] is False
    assert shift["mean_exactly_zero"] is False
    assert shift["source_modified"] is False
    assert shift["mean_deleted"] is False
    assert shift["mean_j_lower"] > 4.0e-17
    assert shift["symmetry_defect_max"] == 0.0
    correction = shift["correction"]
    assert correction["lambda_upper"] < 0.0
    assert correction["lambda_lower"] == pytest.approx(-1.2038809379699564e-15, rel=0.0, abs=1.0e-28)
    assert correction["critical_point_mean_leftover_upper"] <= 3.857333633418152e-21
    assert correction["source_modified"] is False
    assert shift["evolution_error_bound"] is None
    assert float(np.max(np.abs(predecessor["nf512_initial_p_r"]))) == 0.0
    assert float(np.max(np.abs(predecessor["nf512_initial_chi"]))) == 0.0
    assert float(np.max(np.abs(episode[f"{PRIMARY}_frame_coarse_p_r"][0]))) > 1.0
    assert _sha256(ERROR_JSON) == "e0a37e1177b2a483eca0eb3e44243893bbdc04667496b732dd7e63467a97bdd2"
    # The assessment helper is owned elsewhere. Record its bytes; do not treat a later edit as this audit.
    assert ASSESSMENT.is_file()
    assert len(_sha256(ASSESSMENT)) == 64
