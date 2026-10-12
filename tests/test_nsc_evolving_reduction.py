"""Independent checks for the retarded retained-region reduction.

The full-system samples below use DOP853. The reducer module does not.
Prescribed pulses are numerical controls, not regeneration. The incoming
gate is not closed.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from recursive_horizons.nsc_boundary_state import GaussianBoundaryState
from recursive_horizons.nsc_evolving_reduction import (
    covariance_terms,
    evolve_retained_region,
    exterior_kernels,
    fixed_observer_frame,
    free_exterior_propagator,
    hamiltonian_blocks,
    prescribed_six_mode_problem,
)
import recursive_horizons.nsc_evolving_reduction as reduction


def integrate_full_system(hamiltonian, times, columns):
    """Sample ``i psi_dot = H(t) psi`` with DOP853.

    This oracle is intentionally outside the reducer. Agreement is an
    observed discrepancy, not a certified bound.
    """
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
        rtol=1e-11,
        atol=1e-13,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return np.stack([unpack(solution.y[:, index]) for index in range(solution.y.shape[1])])


def _amplitude_errors(reduced, samples, basis):
    projected = np.einsum("ij,tjk->tik", basis.conj().T, samples)
    difference = reduced["amplitudes"] - projected
    norms = np.linalg.norm(difference, axis=(1, 2))
    return float(norms.max()), float(norms[-1])


def _signal(left, right):
    difference = left["amplitudes"] - right["amplitudes"]
    return float(np.linalg.norm(difference, axis=(1, 2)).max())


def measure_prescribed_controls():
    """Observed six-mode errors on uniform grids 32, 64 and 128."""
    problem = prescribed_six_mode_problem()
    hamiltonian = problem["hamiltonian"]
    basis = problem["local_basis"]
    covariance = problem["covariance"]
    identity = np.eye(6, dtype=np.complex128)
    amplitude_max = {}
    amplitude_final = {}
    covariance_final = {}
    car = {}
    trapezoid = {}
    assembly = {}
    fine = None
    fine_samples = None
    fine_times = None
    for steps in (32, 64, 128):
        times = np.linspace(0.0, 1.0, steps + 1)
        reduced = evolve_retained_region(
            hamiltonian,
            basis,
            times,
            identity,
            covariance=covariance,
        )
        samples = integrate_full_system(hamiltonian, times, identity)
        worst, final = _amplitude_errors(reduced, samples, basis)
        observed = basis.conj().T @ samples[-1]
        covariance_error = observed @ covariance @ observed.conj().T
        amplitude_max[steps] = worst
        amplitude_final[steps] = final
        covariance_final[steps] = float(
            np.linalg.norm(reduced["covariance_total"][-1] - covariance_error, ord="fro")
        )
        car[steps] = float(reduced["car_deviation_max"])
        trapezoid[steps] = float(reduced["trapezoid_residual_max"])
        assembly[steps] = float(reduced["covariance_assembly_residual_max"])
        if steps == 128:
            fine = reduced
            fine_samples = samples
            fine_times = times
    omitted_memory = evolve_retained_region(
        hamiltonian, basis, fine_times, identity, memory=False
    )
    omitted_drive = evolve_retained_region(
        hamiltonian, basis, fine_times, identity, outside_drive=False
    )
    coupling_off = evolve_retained_region(
        hamiltonian, basis, fine_times, identity, coupling=False
    )
    cross_off = evolve_retained_region(
        hamiltonian,
        basis,
        fine_times,
        covariance=covariance,
        drop_cross_covariance=True,
    )
    memory_error, _memory_final = _amplitude_errors(omitted_memory, fine_samples, basis)
    drive_error, _drive_final = _amplitude_errors(omitted_drive, fine_samples, basis)
    coupled_coupling_error, _ignored = _amplitude_errors(coupling_off, fine_samples, basis)
    frame = fine["frame"]

    def decoupled(time):
        local, _coupling, exterior = hamiltonian_blocks(hamiltonian(time), frame, 2)
        blocked = np.zeros((6, 6), dtype=np.complex128)
        blocked[:2, :2] = local
        blocked[2:, 2:] = exterior
        return frame @ blocked @ frame.conj().T

    decoupled_samples = integrate_full_system(decoupled, fine_times, identity)
    decoupled_error, _decoupled_final = _amplitude_errors(coupling_off, decoupled_samples, basis)
    cross_difference = fine["covariance_total"] - cross_off["covariance_total"]
    terms = fine["covariance_terms_final"]
    midpoint = len(fine_times) // 2
    kernels = exterior_kernels(
        float(fine_times[-1]),
        float(fine_times[midpoint]),
        coupling_t=fine["coupling_block"][-1],
        propagator_t=fine["exterior_propagator"][-1],
        coupling_s=fine["coupling_block"][midpoint],
        propagator_s=fine["exterior_propagator"][midpoint],
        exterior_covariance=fine["frame_blocks"]["child"],
        initial_cross=fine["frame_blocks"]["child_parent"],
    )
    return {
        "amplitude_error_max": amplitude_max,
        "amplitude_error_final": amplitude_final,
        "covariance_error_final": covariance_final,
        "car_deviation_max": car,
        "trapezoid_residual_max": trapezoid,
        "assembly_residual_max": assembly,
        "signal_omit_memory": _signal(fine, omitted_memory),
        "signal_omit_drive": _signal(fine, omitted_drive),
        "signal_coupling_zero": _signal(fine, coupling_off),
        "signal_drop_cross_final": float(np.linalg.norm(cross_difference[-1], ord="fro")),
        "error_omit_memory": memory_error,
        "error_omit_drive": drive_error,
        "error_coupling_zero_against_coupled_oracle": coupled_coupling_error,
        "error_coupling_zero_against_decoupled_oracle": decoupled_error,
        "term_parent": float(np.linalg.norm(terms["parent"], ord="fro")),
        "term_child": float(np.linalg.norm(terms["child"], ord="fro")),
        "term_initial_parent_child": float(np.linalg.norm(terms["initial_parent_child"], ord="fro")),
        "term_initial_child_parent": float(np.linalg.norm(terms["initial_child_parent"], ord="fro")),
        "term_total": float(np.linalg.norm(terms["total"], ord="fro")),
        "kernel_occupied": float(np.linalg.norm(kernels["occupied"], ord="fro")),
        "kernel_empty": float(np.linalg.norm(kernels["empty"], ord="fro")),
        "kernel_retarded": float(np.linalg.norm(kernels["retarded"], ord="fro")),
        "kernel_initial": float(np.linalg.norm(kernels["initial_occupied"], ord="fro")),
        "kernel_occupied_empty_gap": float(
            np.linalg.norm(kernels["occupied"] - kernels["empty"], ord="fro")
        ),
        "history_final": float(np.linalg.norm(fine["history"][-1], ord="fro")),
        "history_initial": float(np.linalg.norm(fine["history"][0], ord="fro")),
        "initial_exterior": float(np.linalg.norm(fine["initial_exterior"], ord="fro")),
        "omitted_memory_history": float(np.linalg.norm(omitted_memory["history"])),
        "exterior_propagator_shape": list(fine["exterior_propagator"].shape),
        "amplitude_shape": list(fine["amplitudes"].shape),
    }


@pytest.fixture(scope="module")
def controls():
    measured = measure_prescribed_controls()
    json.dumps(measured)
    return measured


def test_reducer_does_not_call_a_full_system_integrator():
    source = Path(reduction.__file__).read_text(encoding="utf-8")
    assert "solve_ivp" not in source
    assert "DOP853" not in source


def test_moving_projector_and_bad_inputs_are_rejected():
    problem = prescribed_six_mode_problem()
    basis = problem["local_basis"]
    times = np.linspace(0.0, 1.0, 5)
    identity = np.eye(6, dtype=np.complex128)
    with pytest.raises(ValueError, match="Berry"):
        evolve_retained_region(problem["hamiltonian"], lambda time: basis, times, identity)
    with pytest.raises(ValueError, match="uniform"):
        evolve_retained_region(problem["hamiltonian"], basis, np.array([0.0, 0.1, 0.4]), identity)

    def nonhermitian(_time):
        matrix = np.zeros((6, 6), dtype=np.complex128)
        matrix[0, 1] = 1.0
        return matrix

    with pytest.raises(ValueError, match="Hermitian"):
        evolve_retained_region(nonhermitian, basis, times, identity)
    with pytest.raises(ValueError, match="0<=C<=I"):
        evolve_retained_region(problem["hamiltonian"], basis, times, covariance=2.0 * identity)
    with pytest.raises(ValueError, match="weights"):
        evolve_retained_region(problem["hamiltonian"], basis, times, weights=np.ones(2))
    with pytest.raises(ValueError, match="cross"):
        evolve_retained_region(
            problem["hamiltonian"],
            basis,
            times,
            identity,
            drop_cross_covariance=True,
        )
    skewed = basis.copy()
    skewed[:, 1] *= 0.2
    with pytest.raises(ValueError, match="QR"):
        fixed_observer_frame(skewed)


def test_frame_keeps_observer_phases_that_qr_would_remove():
    generator = np.random.default_rng(1)
    draw = generator.normal(size=(6, 2)) + 1j * generator.normal(size=(6, 2))
    orthogonal, _triangle = np.linalg.qr(draw)
    phased = orthogonal * np.exp(1j * np.array([0.4, -1.1]))
    completed = fixed_observer_frame(phased)
    assert np.array_equal(completed["frame"][:, :2], phased)
    assert np.array_equal(completed["local_basis"], phased)
    qr_factor, _raw = np.linalg.qr(phased)
    assert not np.allclose(qr_factor, phased)


def test_linear_unknown_is_the_retained_block_only(monkeypatch):
    problem = prescribed_six_mode_problem()
    seen = []
    original = reduction.np.linalg.solve

    def spy(matrix, right_hand):
        seen.append(matrix.shape)
        return original(matrix, right_hand)

    monkeypatch.setattr(reduction.np.linalg, "solve", spy)
    times = np.linspace(0.0, 0.4, 5)
    evolve_retained_region(problem["hamiltonian"], problem["local_basis"], times, np.eye(6))
    assert seen
    assert all(shape == (2, 2) for shape in seen)


def test_constant_pair_matches_matrix_exponential_and_uses_memory():
    hamiltonian = np.array([[0.3, 0.4], [0.4, 0.0]], dtype=np.complex128)
    basis = np.zeros((2, 1), dtype=np.complex128)
    basis[0, 0] = 1.0
    identity = np.eye(2, dtype=np.complex128)
    errors = []
    for steps in (32, 64, 128):
        times = np.linspace(0.0, 1.0, steps + 1)
        reduced = evolve_retained_region(lambda _time: hamiltonian, basis, times, identity)
        reference = np.stack(
            [basis.conj().T @ expm(-1j * hamiltonian * time) for time in times]
        )
        errors.append(float(np.linalg.norm(reduced["amplitudes"] - reference, axis=(1, 2)).max()))
    assert errors[0] > errors[1] > errors[2]
    assert errors[0] / errors[1] > 3.5
    assert errors[1] / errors[2] > 3.5
    assert errors[2] < 2e-6
    times = np.linspace(0.0, 1.0, 65)
    full = evolve_retained_region(lambda _time: hamiltonian, basis, times, identity)
    forgotten = evolve_retained_region(
        lambda _time: hamiltonian, basis, times, identity, memory=False
    )
    signal = _signal(full, forgotten)
    assert signal > 100.0 * errors[1]


def test_independent_oracle_matches_matrix_exponential():
    hamiltonian = np.array([[0.3, 0.4], [0.4, 0.0]], dtype=np.complex128)
    times = np.linspace(0.0, 1.0, 9)
    samples = integrate_full_system(lambda _time: hamiltonian, times, np.eye(2, dtype=np.complex128))
    for index, time in enumerate(times):
        assert np.allclose(samples[index], expm(-1j * hamiltonian * time), atol=1e-9, rtol=0)


def test_time_dependent_blocks_converge():
    def hamiltonian(time):
        pulse = np.sin(2.0 * time)
        return np.array(
            [[pulse, 0.35], [0.35, 0.2 * np.cos(time)]],
            dtype=np.complex128,
        )

    basis = np.zeros((2, 1), dtype=np.complex128)
    basis[0, 0] = 1.0
    identity = np.eye(2, dtype=np.complex128)
    errors = []
    for steps in (32, 64):
        times = np.linspace(0.0, 1.0, steps + 1)
        reduced = evolve_retained_region(hamiltonian, basis, times, identity)
        samples = integrate_full_system(hamiltonian, times, identity)
        worst, _final = _amplitude_errors(reduced, samples, basis)
        errors.append(worst)
    assert errors[0] / errors[1] > 3.5
    assert errors[1] < 2e-4


def test_generic_columns_match_the_retained_map_without_a_second_ode():
    problem = prescribed_six_mode_problem()
    times = np.linspace(0.0, 0.7, 17)
    generator = np.random.default_rng(4)
    columns = generator.normal(size=(6, 6)) + 1j * generator.normal(size=(6, 6))
    mapped = evolve_retained_region(
        problem["hamiltonian"], problem["local_basis"], times, np.eye(6, dtype=np.complex128)
    )
    direct = evolve_retained_region(
        problem["hamiltonian"], problem["local_basis"], times, columns
    )
    predicted = np.einsum("tij,jk->tik", mapped["retained_map"], columns)
    assert np.allclose(direct["amplitudes"], predicted, atol=1e-9, rtol=0)
    short = columns[:, :2]
    direct_short = evolve_retained_region(
        problem["hamiltonian"], problem["local_basis"], times, short
    )
    predicted_short = np.einsum("tij,jk->tik", mapped["retained_map"], short)
    assert np.allclose(direct_short["amplitudes"], predicted_short, atol=1e-9, rtol=0)
    orthogonal, _triangle = np.linalg.qr(columns[:, :3])
    weights = np.array([0.15, 0.4, 0.8])
    weighted = evolve_retained_region(
        problem["hamiltonian"],
        problem["local_basis"],
        times,
        orthogonal,
        weights=weights,
    )
    column_state = (orthogonal * weights) @ orthogonal.conj().T
    from_covariance = evolve_retained_region(
        problem["hamiltonian"],
        problem["local_basis"],
        times,
        covariance=column_state,
    )
    assert np.allclose(
        weighted["column_weight_covariance"],
        from_covariance["covariance_total"],
        atol=1e-9,
        rtol=0,
    )


def test_history_accumulator_is_the_factored_trapezoid():
    problem = prescribed_six_mode_problem()
    times = np.linspace(0.0, 1.0, 21)
    reduced = evolve_retained_region(
        problem["hamiltonian"],
        problem["local_basis"],
        times,
        np.eye(6, dtype=np.complex128),
    )
    assert reduced["exterior_propagator"].shape[-1] == 4
    assert reduced["amplitudes"].shape[1:] == (2, 6)
    assert np.linalg.norm(reduced["history"][0]) == 0.0
    assert np.linalg.norm(reduced["initial_exterior"]) > 0.1
    step = float(times[1] - times[0])
    accumulator = np.zeros_like(reduced["history"][0])
    coupling = reduced["coupling_block"]
    propagator = reduced["exterior_propagator"]
    amplitudes = reduced["amplitudes"]
    for index in range(len(times) - 1):
        kernel_now = propagator[index].conj().T @ coupling[index].conj().T
        kernel_new = propagator[index + 1].conj().T @ coupling[index + 1].conj().T
        accumulator = accumulator + (step / 2.0) * (
            kernel_now @ amplitudes[index] + kernel_new @ amplitudes[index + 1]
        )
        assert np.allclose(accumulator, reduced["history"][index + 1], atol=1e-12, rtol=0)
    frame = reduced["frame"]

    def exterior(time):
        return hamiltonian_blocks(problem["hamiltonian"](time), frame, 2)[2]

    recomputed = free_exterior_propagator(exterior, times, substeps=1)
    assert np.allclose(recomputed, propagator, atol=1e-12, rtol=0)


def test_kernels_match_the_fixed_frame_boundary_state():
    problem = prescribed_six_mode_problem()
    times = np.linspace(0.0, 1.0, 25)
    reduced = evolve_retained_region(
        problem["hamiltonian"],
        problem["local_basis"],
        times,
        covariance=problem["covariance"],
    )
    mask = np.array([True, True, False, False, False, False])
    frame_covariance = reduced["frame"].conj().T @ problem["covariance"] @ reduced["frame"]
    state = GaussianBoundaryState(frame_covariance, mask)
    later, earlier = 18, 5
    produced = exterior_kernels(
        float(times[later]),
        float(times[earlier]),
        coupling_t=reduced["coupling_block"][later],
        propagator_t=reduced["exterior_propagator"][later],
        coupling_s=reduced["coupling_block"][earlier],
        propagator_s=reduced["exterior_propagator"][earlier],
        exterior_covariance=reduced["frame_blocks"]["child"],
        initial_cross=reduced["frame_blocks"]["child_parent"],
    )
    reference = state.kernels(
        float(times[later]),
        float(times[earlier]),
        link_t=reduced["coupling_block"][later],
        child_t=reduced["exterior_propagator"][later],
        link_s=reduced["coupling_block"][earlier],
        child_s=reduced["exterior_propagator"][earlier],
    )
    for name in produced:
        assert np.allclose(produced[name], reference[name], atol=1e-10, rtol=0)
    assert np.linalg.norm(produced["occupied"] - produced["empty"]) > 1e-4
    assert np.linalg.norm(produced["initial_occupied"]) > 0.0
    same_time = exterior_kernels(
        float(times[later]),
        float(times[later]),
        coupling_t=reduced["coupling_block"][later],
        propagator_t=reduced["exterior_propagator"][later],
        coupling_s=reduced["coupling_block"][later],
        propagator_s=reduced["exterior_propagator"][later],
        exterior_covariance=reduced["frame_blocks"]["child"],
        initial_cross=reduced["frame_blocks"]["child_parent"],
    )
    source = reduced["coupling_block"][later] @ reduced["exterior_propagator"][later]
    assert np.allclose(same_time["retarded"], -0.5j * (source @ source.conj().T), atol=1e-12, rtol=0)
    blocks = (
        reduced["frame_blocks"]["parent"],
        reduced["frame_blocks"]["parent_child"],
        reduced["frame_blocks"]["child_parent"],
        reduced["frame_blocks"]["child"],
    )
    rebuilt = covariance_terms(reduced["parent_response"][0], reduced["child_response"][0], blocks)
    assert np.allclose(rebuilt["total"], blocks[0], atol=1e-12, rtol=0)
    assert np.linalg.norm(rebuilt["initial_parent_child"]) == 0.0


def test_rotated_fixed_basis_matches_the_oracle():
    problem = prescribed_six_mode_problem()
    generator = np.random.default_rng(2)
    draw = generator.normal(size=(6, 2)) + 1j * generator.normal(size=(6, 2))
    orthogonal, _triangle = np.linalg.qr(draw)
    basis = orthogonal * np.exp(1j * np.array([0.4, -1.1]))
    times = np.linspace(0.0, 1.0, 65)
    identity = np.eye(6, dtype=np.complex128)
    reduced = evolve_retained_region(problem["hamiltonian"], basis, times, identity)
    samples = integrate_full_system(problem["hamiltonian"], times, identity)
    worst, _final = _amplitude_errors(reduced, samples, basis)
    assert worst < 2e-3
    assert np.array_equal(reduced["local_basis"], basis)
    assert reduced["car_deviation_max"] < 1e-3


def test_prescribed_window_is_a_control_not_a_trajectory():
    problem = prescribed_six_mode_problem()
    assert problem["regeneration"] is False
    assert problem["coupled_trajectory"] is False
    assert problem["prescribed_history"] is True
    assert np.allclose(problem["hamiltonian"](0.0), problem["window"], atol=1e-12, rtol=0)
    perturbation = problem["perturbation"]
    assert perturbation[2, 2] == 1.0
    assert perturbation[0, 2] == 0.1
    assert perturbation[2, 0] == 0.1
    assert np.count_nonzero(perturbation) == 3


def test_six_mode_grids_converge_and_negative_controls_separate(controls):
    amplitude = controls["amplitude_error_max"]
    covariance = controls["covariance_error_final"]
    car = controls["car_deviation_max"]
    assert amplitude[32] > amplitude[64] > amplitude[128]
    assert amplitude[32] / amplitude[64] > 3.5
    assert amplitude[64] / amplitude[128] > 3.5
    assert amplitude[128] < 1e-4
    assert covariance[32] > covariance[64] > covariance[128]
    assert covariance[32] / covariance[64] > 3.5
    assert car[32] > car[64] > car[128]
    assert car[32] / car[64] > 3.5
    assert controls["trapezoid_residual_max"][128] < 1e-12
    assert controls["assembly_residual_max"][128] < 1e-12
    fine = amplitude[128]
    assert controls["signal_omit_memory"] > 100.0 * fine
    assert controls["signal_omit_drive"] > 100.0 * fine
    assert controls["signal_coupling_zero"] > 100.0 * fine
    assert controls["error_omit_memory"] > 20.0 * fine
    assert controls["error_omit_drive"] > 20.0 * fine
    assert controls["error_coupling_zero_against_coupled_oracle"] > 20.0 * fine
    assert controls["error_coupling_zero_against_decoupled_oracle"] < 3.0 * fine
    assert controls["signal_drop_cross_final"] > 50.0 * covariance[128]
    assert controls["term_child"] > 1e-2
    assert controls["term_initial_parent_child"] > 1e-4
    assert controls["term_initial_child_parent"] > 1e-4
    assert controls["kernel_occupied_empty_gap"] > 1e-3
    assert controls["kernel_initial"] > 0.0
    assert controls["history_initial"] == 0.0
    assert controls["history_final"] > 1e-2
    assert controls["initial_exterior"] > 0.1
    assert controls["omitted_memory_history"] == 0.0
    assert controls["exterior_propagator_shape"][-1] == 4
    assert controls["amplitude_shape"][1:] == [2, 6]


if __name__ == "__main__":
    print(json.dumps(measure_prescribed_controls(), indent=2, sort_keys=True))
