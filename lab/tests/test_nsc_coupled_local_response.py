"""Conditional Galerkin consumer. The prescribed six-mode Schur identity is not retested.

The small band uses the rank-6 packet and an independent geometry pulse.
The stored-episode test uses the saved nf256 frames, columns, and weights.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from recursive_horizons.nsc_coupled_local_response import (
    EPISODE_NPZ,
    DECLARED_INTERPOLATION,
    DEFAULT_JSON,
    DEFAULT_NPZ,
    MeasuredBudget,
    action_residual,
    column_cross_control,
    correlated_cross_separation,
    compare_retained,
    evolve_streamed,
    fourier_galerkin_factors,
    frame_source_blocks,
    galerkin_matrix,
    hermite_with_slopes,
    independent_control,
    independent_radial_density,
    interpolation_scale,
    hermite_rows,
    interpolate_rows,
    load_episode_case,
    load_successor_inputs,
    multiplication_matrix,
    occupations_and_coherence,
    one_body_blocks,
    plan_remaining,
    prolonged_nodal_rates,
    quadratic_estimate,
    request_identity,
    region_observer,
    resume_disposition,
    reject_linear_operator_hamiltonian,
    sha256_file,
    source_covariance,
    source_id,
    spinor_local_moments,
    successor_identity,
    successor_resume_disposition,
    verify_stored_operator,
    build_case_grid,
    hamiltonian_from_radial,
    stored_radial,
)
from recursive_horizons.nsc_evolving_reduction import evolve_retained_region
from recursive_horizons.nsc_spherical_coupling import OCCUPATIONS


@pytest.fixture(scope="module")
def control():
    prepared = independent_control(24)
    assert prepared["stored_episode"] is False
    assert prepared["regeneration"] is False
    return prepared


def _dop853(hamiltonian, times, columns):
    initial = np.asarray(columns, dtype=np.complex128)
    dimension, width = initial.shape

    def pack(state):
        return np.concatenate((state.real.ravel(), state.imag.ravel()))

    def unpack(flat):
        half = dimension * width
        return (flat[:half] + 1j * flat[half:]).reshape(dimension, width)

    def derivative(_time, flat):
        matrix = np.asarray(hamiltonian(_time), dtype=np.complex128)
        return pack(-1j * (matrix @ unpack(flat)))

    solution = solve_ivp(
        derivative,
        (float(times[0]), float(times[-1])),
        pack(initial),
        t_eval=np.asarray(times, dtype=float),
        method="DOP853",
        rtol=1e-9,
        atol=1e-11,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return np.stack([unpack(solution.y[:, index]) for index in range(solution.y.shape[1])])


def test_fourier_matrix_matches_apply_dirac_on_independent_controls(control):
    grid = control["grid"]
    factors = control["factors"]
    assert factors["fine_identity_image"] is False
    assert factors["fourier_band"].shape == (grid.nq, grid.nf)
    radial = independent_radial_density(grid, 0.2)
    residual = action_residual(grid, factors, radial, control["phi0"], control["phi1"])
    assert residual["relative"] < 1e-12
    matrix = galerkin_matrix(factors, radial)
    assert matrix.shape == (2 * grid.nf, 2 * grid.nf)
    multiplied = multiplication_matrix(factors, radial)
    sandwich = grid.U_f.conj().T @ (radial[:, None] * grid.U_f)
    sandwich = 0.5 * (sandwich + sandwich.conj().T)
    assert np.allclose(multiplied, sandwich, rtol=0.0, atol=1e-10)
    leaked = grid.fine.momentum @ grid.U_f
    projected = grid.U_f @ (grid.U_f.conj().T @ leaked)
    invariance = np.linalg.norm(leaked - projected) / np.linalg.norm(leaked)
    assert invariance < 1e-10


def test_time_dependent_operator_and_source_covariance_do_not_commute(control):
    matrix_0 = galerkin_matrix(control["factors"], control["radial"](0.0))
    matrix_1 = galerkin_matrix(control["factors"], control["radial"](0.4))
    assert np.linalg.norm(matrix_0 @ matrix_1 - matrix_1 @ matrix_0) > 1e-6
    covariance = source_covariance(control["columns"], control["weights"])
    assert np.linalg.norm(matrix_0 @ covariance - covariance @ matrix_0) > 1e-6
    assert np.max(control["weights"]) > np.min(control["weights"])
    np.testing.assert_allclose(control["weights"], OCCUPATIONS)


def test_observer_is_the_initial_region_pair_and_coupling_is_galerkin(control):
    observer, columns, gram = region_observer(control["phi0"], control["phi1"], region=0)
    assert np.array_equal(observer, columns[:, :2])
    assert np.max(np.abs(gram - np.eye(6))) < 1e-10
    times = np.linspace(0.0, 0.04, 5)
    comparison = compare_retained(
        control["hamiltonian"],
        control["observer"],
        times,
        control["columns"],
        control["weights"],
    )
    assert np.array_equal(comparison["result"]["local_basis"], control["observer"])
    assert comparison["allocation"]["time_indexed_exterior_propagator_bytes"] == 0
    assert comparison["allocation"]["history_shape"][-1] == 6
    assert comparison["allocation"]["history_shape"][-2] == comparison["allocation"]["exterior_dimension"]
    assert comparison["allocation"]["avoided_dense_propagator_bytes"] > comparison["allocation"]["history_bytes"]
    np.testing.assert_allclose(comparison["occupation_retained"][0], control["weights"][:2], atol=1e-12)
    assert np.linalg.norm(comparison["result"]["coupling_block"][0]) > 0.0


def test_streamed_occupation_converges_to_the_independent_field(control):
    errors = []
    midpoint_errors = []
    for nodes in (5, 9, 17):
        times = np.linspace(0.0, 0.05, nodes)
        comparison = compare_retained(
            control["hamiltonian"],
            control["observer"],
            times,
            control["columns"],
            control["weights"],
        )
        samples = _dop853(control["hamiltonian"], times, control["columns"])
        projected = np.einsum("ij,tjk->tik", control["observer"].conj().T, samples)
        occupation, _coherence = occupations_and_coherence(projected, control["weights"])
        errors.append(float(np.max(np.abs(comparison["occupation_retained"] - occupation))))
        midpoint_errors.append(float(np.max(np.abs(comparison["occupation_full"] - occupation))))
    assert errors[2] < errors[1] < errors[0]
    assert errors[2] < errors[0] / 3.0
    assert midpoint_errors[2] < errors[0]


def test_omissions_are_active_and_history_is_not_a_dense_propagator(control):
    times = np.linspace(0.0, 0.05, 9)
    baseline = compare_retained(
        control["hamiltonian"],
        control["observer"],
        times,
        control["columns"],
        control["weights"],
    )
    reduction_error = baseline["occupation_error"]["max_abs"]
    assert baseline["history_norm"] > 0.0
    assert baseline["initial_exterior_norm"] > 0.0
    active = []
    for name, flags in (
        ("memory", {"memory": False}),
        ("outside_drive", {"outside_drive": False}),
        ("coupling", {"coupling": False}),
    ):
        result, occupation, _coherence = evolve_streamed(
            control["hamiltonian"],
            control["observer"],
            times,
            control["columns"],
            control["weights"],
            **flags,
        )
        separation = float(np.max(np.abs(occupation - baseline["occupation_full"])))
        if separation > 10.0 * reduction_error:
            active.append(name)
        if name == "memory":
            assert np.linalg.norm(result["history"]) == 0.0
        if name == "outside_drive":
            assert np.isclose(np.linalg.norm(result["initial_exterior"]), baseline["initial_exterior_norm"])
        assert result["history"].shape[-1] == control["columns"].shape[1]
        assert int(result["allocation"]["time_indexed_exterior_propagator_bytes"]) == 0
    assert {"memory", "outside_drive", "coupling"} <= set(active)


def test_source_cross_block_is_preserved_rather_than_replaced(control):
    covariance = source_covariance(control["columns"], control["weights"])
    times = np.linspace(0.0, 0.03, 4)
    kept = evolve_retained_region(
        control["hamiltonian"],
        control["observer"],
        times,
        covariance=covariance,
        backend="streamed",
    )
    dropped = evolve_retained_region(
        control["hamiltonian"],
        control["observer"],
        times,
        covariance=covariance,
        drop_cross_covariance=True,
        backend="streamed",
    )
    blocks = frame_source_blocks(kept["frame"], control["columns"], control["weights"])
    delta = np.linalg.norm(kept["covariance_total"] - dropped["covariance_total"])
    assert blocks["cross_norm"] < 1e-8
    assert delta < 1e-8
    np.testing.assert_allclose(np.real(np.diag(blocks["parent"])), control["weights"][:2], atol=1e-10)
    assert int(kept["allocation"]["time_indexed_exterior_propagator_bytes"]) == 0
    # The covariance path evolves one column per mode, so its history is
    # (time, N_outside, N) rather than a square propagator.
    assert kept["history_of_frame"].shape[-1] == kept["frame"].shape[0]
    assert kept["history_of_frame"].shape[-2] == kept["frame"].shape[0] - 2


def test_dense_callback_rejects_a_linear_operator():
    name = reject_linear_operator_hamiltonian()
    assert name in {"TypeError", "ValueError"}


def test_linear_interpolation_keeps_nodes_and_names_the_quadratic_scale():
    times = np.linspace(0.0, 1.0, 5)
    values = (times ** 2)[:, None]
    nodes = interpolate_rows(times, values, times)
    assert np.allclose(nodes, values, rtol=0.0, atol=0.0)
    midpoint = interpolate_rows(times, values, [0.125])[0, 0]
    assert midpoint == pytest.approx(0.03125)
    scale = interpolation_scale(times, values)
    assert scale["midpoint_linear_scale"] == pytest.approx(0.015625)
    assert abs(midpoint - 0.125 ** 2) == pytest.approx(scale["midpoint_linear_scale"])


def test_budget_estimate_shrinks_the_prefix():
    assert quadratic_estimate(2.0, 2, 4) == pytest.approx(8.0)
    assert quadratic_estimate(2.0, 2, 4, substeps=2) == pytest.approx(16.0)
    wide = plan_remaining(1.0e9, 1.0, 2, 10)
    assert wide["main_steps"] == 10
    assert wide["refine"] is True
    assert wide["omissions"] is True
    tight = plan_remaining(0.01, 10.0, 2, 10)
    assert tight["main_steps"] == 2
    assert tight["omissions"] is False
    assert tight["refine"] is False
    budget = MeasuredBudget(0.05)
    assert budget.allows(1.0) is False


def test_loader_preserves_stored_frames_bindings_and_partner_gap():
    case = load_episode_case("nf256_dt_0_0005")
    with np.load(EPISODE_NPZ, allow_pickle=False) as data:
        assert np.array_equal(case["frame_Q"], data["nf256_dt_0_0005_frame_Q"])
        assert np.array_equal(case["phi0"], data["nf256_initial_phi0"])
        assert np.array_equal(case["phi1"], data["nf256_initial_phi1"])
        assert np.array_equal(case["weights"], data["nf256_occupations"])
        assert np.array_equal(case["binding_values"], data["nf256_binding_values"])
        gap = float(np.max(np.abs(data["nf256_dt_0_0005_frame_Q"] - data["nf256_dt_0_00025_frame_Q"])))
    assert case["time_resolution_geometry_gap"] == gap
    assert case["episode_verdict_used_as_acceptance"] is False
    assert case["geometry_rate_series_stored"] is False
    assert case["interpolation"]["midpoint_linear_scale"] > 0.0


def test_stored_nf256_prefix_uses_the_same_preparation():
    case = load_episode_case("nf256_dt_0_0005")
    grid, factors = build_case_grid(case)
    verified = verify_stored_operator(grid, factors, case)
    assert verified["initial_action_relative"] < 1e-12
    assert verified["final_rate_relative"] < 1e-12
    assert verified["initial_prolongation_gap"] < 1e-8
    observer, columns, gram = region_observer(case["phi0"], case["phi1"], region=0)
    assert np.array_equal(observer, columns[:, :2])
    assert np.max(np.abs(gram - np.eye(6))) < 1e-12
    np.testing.assert_allclose(case["weights"], OCCUPATIONS)
    hamiltonian = hamiltonian_from_radial(factors, stored_radial(case))
    comparison = compare_retained(
        hamiltonian,
        observer,
        case["frame_time"][:3],
        columns,
        case["weights"],
    )
    assert np.array_equal(comparison["result"]["local_basis"], observer)
    assert comparison["allocation"]["time_indexed_exterior_propagator_bytes"] == 0
    assert comparison["allocation"]["history_shape"] == [3, 2 * grid.nf - 2, 6]
    assert comparison["occupation_error"]["relative_to_signal"] < 1e-3
    np.testing.assert_allclose(comparison["occupation_retained"][0], case["weights"][:2], atol=1e-12)


def test_hermite_variant_keeps_nodes_and_moves_midpoints():
    times = np.linspace(0.0, 1.0, 5)
    values = (times ** 3)[:, None]
    assert np.allclose(hermite_rows(times, values, times), values, rtol=0.0, atol=1e-15)
    midpoint = 0.125
    linear = interpolate_rows(times, values, [midpoint])[0, 0]
    curved = hermite_rows(times, values, [midpoint])[0, 0]
    assert linear == pytest.approx(0.0078125)
    assert abs(curved - linear) > 1e-4


def test_rank6_cross_is_inactive_and_six_mode_cross_is_active(control):
    covariance = source_covariance(control["columns"], control["weights"])
    times = np.linspace(0.0, 0.02, 3)
    kept = evolve_retained_region(
        control["hamiltonian"],
        control["observer"],
        times,
        covariance=covariance,
        backend="streamed",
    )
    blocks = frame_source_blocks(kept["frame"], control["columns"], control["weights"])
    assert blocks["cross_norm"] < 1e-8
    measured = correlated_cross_separation()
    assert measured["active"] is True
    assert measured["cross_norm"] > 1e-4
    assert measured["covariance_separation"] > 1e-4
    assert measured["stores_dense_exterior_propagator"] is False


def _finished_pilot(source):
    identity = request_identity("pilot", "nf512_dt_0_0005", 1, source)
    return {
        "status": "MEASURED_CONDITIONAL_REDUCTION",
        "requested_domain": "pilot",
        "requested_case": "nf512_dt_0_0005",
        "requested_substeps": 1,
        "headline": {"complete": True, "t_start": 0.0, "t_end": 0.035},
        "phases": {"comparison": {"complete": True, "t_end": 0.035, "substeps": 1}},
        "requests": {"pilot": identity},
    }


def test_resume_rejects_a_different_domain_substeps_source_or_settings():
    source = source_id(np.array([[1.0]]), np.array([[0.0]]), np.array([0.75, 0.75]))
    finished = _finished_pilot(source)
    same = request_identity("pilot", "nf512_dt_0_0005", 1, source)
    assert resume_disposition(finished, same) == "return"
    assert resume_disposition(finished, request_identity("stored-frames", "nf512_dt_0_0005", 2, source)) == "reject"
    assert resume_disposition(finished, request_identity("pilot", "nf512_dt_0_0005", 2, source)) == "reject"
    assert resume_disposition(finished, request_identity("pilot", "nf512_dt_0_00025", 1, source)) == "reject"
    assert resume_disposition(finished, request_identity("pilot", "nf512_dt_0_0005", 1, "other-source")) == "reject"
    changed = request_identity(
        "pilot",
        "nf512_dt_0_0005",
        1,
        source,
        settings={
            "backend": "streamed",
            "interpolation": "cubic-hermite",
            "memory": True,
            "outside_drive": True,
            "coupling": True,
            "drop_cross_covariance": False,
        },
    )
    assert changed["interpolation"] != DECLARED_INTERPOLATION
    assert resume_disposition(finished, changed) == "reject"
    legacy = {
        "status": "MEASURED_CONDITIONAL_REDUCTION",
        "requested_domain": "pilot",
        "requested_case": "nf512_dt_0_0005",
        "requested_substeps": 1,
        "headline": {"complete": True, "t_end": 0.035},
    }
    assert resume_disposition(legacy, request_identity("stored-frames", "nf512_dt_0_0005", 2, source)) == "reject"
    assert resume_disposition(legacy, same) == "return"


def test_saved_record_resume_rejects_a_changed_request_and_keeps_the_pilot():
    record_path = Path(__file__).resolve().parents[1] / "results" / "development" / "nsc-coupled-local-response-v1.json"
    record = json.loads(record_path.read_text())
    case = load_episode_case("nf512_dt_0_0005")
    source = source_id(case["phi0"], case["phi1"], case["weights"])
    assert record["headline"]["t_end"] == 0.035
    assert record["cpu_seconds_this_process"] == 18.504784
    assert resume_disposition(record, request_identity("stored-frames", "nf512_dt_0_0005", 2, source)) == "return"
    assert resume_disposition(record, request_identity("stored-frames", "nf512_dt_0_0005", 1, source)) == "reject"
    assert resume_disposition(record, request_identity("pilot", "nf512_dt_0_0005", 2, source)) == "reject"
    assert resume_disposition(record, request_identity("pilot", "nf256_dt_0_0005", 1, source)) == "reject"
    assert resume_disposition(record, request_identity("pilot", "nf512_dt_0_0005", 1, source)) == "return"


def test_saved_endpoint_spinor_projects_onto_the_initial_observer():
    case = load_episode_case("nf512_dt_0_0005")
    observer, _columns, _gram = region_observer(case["phi0"], case["phi1"], region=0)
    initial = spinor_local_moments(observer, case["phi0"], case["phi1"], case["weights"])
    final = spinor_local_moments(observer, case["final_phi0"], case["final_phi1"], case["weights"])
    np.testing.assert_allclose(initial["occupation"], case["weights"][:2], atol=1e-12)
    assert np.max(np.abs(np.asarray(final["occupation"]) - np.asarray(initial["occupation"]))) > 1e-3


def test_v1_payload_hashes_stay_immutable():
    assert sha256_file(DEFAULT_JSON) == "33fe043b30b037432801f12905b3545a20ca2070b9edc46bd3a11aa0d76201a1"
    assert sha256_file(DEFAULT_NPZ) == "6f4bfce40a4e34a958f0ac1e3b6f763484505a9800b3e34c7362125cbe160ad1"


def test_rate_hermite_uses_supplied_slopes_and_keeps_nodes():
    times = np.linspace(0.0, 1.0, 5)
    values = (times ** 2)[:, None]
    slopes = (2.0 * times)[:, None]
    assert np.allclose(hermite_with_slopes(times, values, slopes, times), values, rtol=0.0, atol=1e-12)
    midpoint = hermite_with_slopes(times, values, slopes, [0.125])[0, 0]
    assert midpoint == pytest.approx(0.125 ** 2, abs=1e-12)
    finite = hermite_rows(times, values, [0.125])[0, 0]
    assert abs(midpoint - finite) > 1e-4


def test_hamiltonian_cache_returns_the_same_matrix(control):
    radial = control["radial"]
    hamiltonian = hamiltonian_from_radial(control["factors"], radial)
    first = hamiltonian(0.2)
    second = hamiltonian(0.2)
    assert first is second
    assert hamiltonian.stats["builds"] == 1
    assert hamiltonian.stats["hits"] == 1
    other = hamiltonian(0.4)
    assert other is not first
    assert hamiltonian.stats["builds"] == 2


def test_successor_observer_is_the_t0_pair_and_the_cross_is_actual():
    before = sha256_file(EPISODE_NPZ)
    loaded = load_successor_inputs()
    assert sha256_file(EPISODE_NPZ) == before
    assert np.array_equal(loaded["observer"], loaded["origin_columns"][:, :2])
    assert not np.array_equal(loaded["observer"], loaded["transported_columns"][:, :2])
    assert np.array_equal(loaded["phi0_frames"][0], loaded["transported_columns"][: loaded["phi0"].shape[0]])
    np.testing.assert_allclose(loaded["weights"], OCCUPATIONS)
    blocks = loaded["initial_blocks"]
    assert blocks["admissible"] is True
    assert blocks["cross_frobenius"] > 1e-2
    assert float(np.min(blocks["covariance_eigenvalues"])) >= -1e-8
    assert float(np.max(blocks["covariance_eigenvalues"])) <= 1.0 + 1e-8
    assert float(np.max(blocks["parent_eigenvalues"])) <= 1.0 + 1e-8
    assert loaded["arc_x_start"] == pytest.approx(3.96875)
    assert loaded["structure_renewed"] is False
    assert loaded["old_toy_B_not_used"] is True
    assert loaded["state_reset"] is False
    identity = successor_identity(loaded)
    finished = {
        "schema": "NSC-COUPLED-LOCAL-RESPONSE-v2",
        "status": "MEASURED_SUCCESSOR_RESPONSE",
        "request": identity,
        "phases": {"full_window": {"complete": True, "substeps": identity["substeps"]}},
    }
    assert successor_resume_disposition(finished, identity) == "return"
    changed = dict(identity)
    changed["observer"] = "0" * 64
    assert successor_resume_disposition(finished, changed) == "reject"
    assert successor_resume_disposition(finished, successor_identity(loaded, substeps=1)) == "reject"


def test_prolonged_rate_reproduces_the_saved_coarse_sample():
    loaded = load_successor_inputs()
    grid, _factors = build_case_grid({
        "phi0": loaded["phi0"],
        "frame_Q": loaded["frame_Q"],
        "coarse_Q": loaded["coarse_Q"][0],
    })
    _fine, report = prolonged_nodal_rates(grid, loaded["coarse_Q_dot"])
    assert report["usable"] is True
    assert report["pullback_max_abs"] < 1e-10


def test_actual_column_cross_matches_the_covariance_omission(control):
    columns = np.array(control["columns"], copy=True)
    cosine = float(np.cos(0.4))
    sine = float(np.sin(0.4))
    rotated = np.array(columns, copy=True)
    rotated[:, 0] = cosine * columns[:, 0] + sine * columns[:, 2]
    rotated[:, 2] = -sine * columns[:, 0] + cosine * columns[:, 2]
    weights = control["weights"]
    blocks = one_body_blocks(control["observer"], rotated, weights)
    assert blocks["cross_frobenius"] > 1e-3
    assert blocks["admissible"] is True
    times = np.linspace(0.0, 0.04, 3)
    split = column_cross_control(
        control["hamiltonian"],
        control["observer"],
        times,
        rotated,
        weights,
        substeps=1,
    )
    assert split["superposition_max_abs"] < 1e-8
    assert split["active"] is True
    assert split["synthetic_six_mode_substituted"] is False
    covariance = source_covariance(rotated, weights)
    kept = evolve_retained_region(
        control["hamiltonian"],
        control["observer"],
        times,
        covariance=covariance,
        backend="streamed",
    )
    dropped = evolve_retained_region(
        control["hamiltonian"],
        control["observer"],
        times,
        covariance=covariance,
        drop_cross_covariance=True,
        backend="streamed",
    )
    assert np.allclose(split["covariance_total"], kept["covariance_total"], rtol=0.0, atol=1e-8)
    assert np.allclose(split["covariance_without_cross"], dropped["covariance_total"], rtol=0.0, atol=1e-8)
    assert int(kept["allocation"]["time_indexed_exterior_propagator_bytes"]) == 0
