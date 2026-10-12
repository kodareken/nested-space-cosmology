"""Independent checks for the retained-region reducer.

The generator is one prescribed four-dimensional Hermitian matrix function
of time, in a fixed complex observer basis. No coupled H(g(t)) is supplied
and none is built. Full-system samples use DOP853 outside the reducer.
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from recursive_horizons.nsc_boundary_state import GaussianBoundaryState
from recursive_horizons.nsc_evolving_reduction import (
    covariance_terms,
    evolve_retained_region,
    exterior_kernels,
    fixed_observer_frame,
    hamiltonian_blocks,
)
import recursive_horizons.nsc_evolving_reduction as reduction


_DRAW = np.array(
    [
        [0.7 + 0.4j, -0.2 + 0.5j],
        [-0.3 + 0.8j, 0.6 - 0.1j],
        [0.5 - 0.6j, 0.4 + 0.7j],
        [0.2 + 0.3j, -0.8 + 0.2j],
    ],
    dtype=np.complex128,
)
_RETAINED = 2
_DIMENSION = 4
_FINAL = 0.5
_GRIDS = (16, 32, 64)


def _local_basis():
    orthogonal, _triangle = np.linalg.qr(_DRAW)
    phases = np.exp(1j * np.array([0.55, -1.15]))
    return orthogonal * phases


def _prescribed_blocks(time):
    """Explicit frame blocks. Time dependence is an input, not a solved trajectory."""
    t = float(time)
    sine = np.sin(1.3 * t)
    cosine = np.cos(0.6 * t)
    local = np.array(
        [
            [0.4 + 0.15 * sine, (0.12 - 0.07j) + 0.05 * cosine],
            [(0.12 + 0.07j) + 0.05 * cosine, -0.25 + 0.1 * cosine],
        ],
        dtype=np.complex128,
    )
    exterior = np.array(
        [
            [0.55 + 0.08 * cosine, -0.1 + 0.14j],
            [-0.1 - 0.14j, 0.2 - 0.12 * sine],
        ],
        dtype=np.complex128,
    )
    coupling = np.array(
        [
            [0.3 + 0.1 * sine, (0.15 - 0.2j) + 0.04j * cosine],
            [0.08 + 0.15j, -0.22 + 0.06 * sine],
        ],
        dtype=np.complex128,
    )
    local = 0.5 * (local + local.conj().T)
    exterior = 0.5 * (exterior + exterior.conj().T)
    return local, coupling, exterior


def _pack_blocks(local, coupling, exterior):
    full = np.zeros((_DIMENSION, _DIMENSION), dtype=np.complex128)
    full[:_RETAINED, :_RETAINED] = local
    full[:_RETAINED, _RETAINED:] = coupling
    full[_RETAINED:, :_RETAINED] = coupling.conj().T
    full[_RETAINED:, _RETAINED:] = exterior
    return full


def _plane_rotation(first, second, angle, phase):
    unitary = np.eye(_DIMENSION, dtype=np.complex128)
    cosine, sine = np.cos(angle), np.sin(angle)
    factor = np.exp(1j * phase)
    unitary[first, first] = cosine
    unitary[second, second] = cosine
    unitary[first, second] = -sine * factor
    unitary[second, first] = sine * np.conjugate(factor)
    return unitary


def _frame_covariance():
    mixing = _plane_rotation(0, 2, 0.45, 0.7) @ _plane_rotation(1, 3, -0.35, -1.1)
    weights = np.diag(np.array([0.2, 0.45, 0.7, 0.9], dtype=float))
    covariance = mixing @ weights @ mixing.conj().T
    return 0.5 * (covariance + covariance.conj().T)


@lru_cache(maxsize=1)
def example():
    basis = _local_basis()
    completed = fixed_observer_frame(basis)
    frame = completed["frame"]

    def hamiltonian(time):
        return frame @ _pack_blocks(*_prescribed_blocks(time)) @ frame.conj().T

    covariance = frame @ _frame_covariance() @ frame.conj().T
    covariance = 0.5 * (covariance + covariance.conj().T)
    return {
        "basis": basis,
        "frame": frame,
        "hamiltonian": hamiltonian,
        "covariance": covariance,
        "frame_covariance": _frame_covariance(),
    }


def integrate_full_system(hamiltonian, times, columns):
    """Sample ``i psi_dot = H(t) psi`` with DOP853, outside the reducer."""
    initial = np.asarray(columns, dtype=np.complex128)
    if initial.ndim == 1:
        initial = initial.reshape(-1, 1)
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


def integrate_memory_free(hamiltonian, frame, times, column):
    """Retained drive without the history integral, integrated outside the reducer.

    State is ``(X, W)`` with ``i W_dot = E W`` and
    ``X_dot = -i A X - i B W Y0``. No accumulator ``Z``.
    """
    coefficients = frame.conj().T @ np.asarray(column, dtype=np.complex128)
    if coefficients.ndim == 1:
        coefficients = coefficients.reshape(-1, 1)
    retained = _RETAINED
    initial_retained = coefficients[:retained]
    initial_exterior = coefficients[retained:]
    exterior = frame.shape[0] - retained
    initial_propagator = np.eye(exterior, dtype=np.complex128)
    retained_size = retained * initial_retained.shape[1]
    propagator_size = exterior * exterior

    def pack(amplitude, propagator):
        return np.concatenate(
            (
                amplitude.real.ravel(),
                amplitude.imag.ravel(),
                propagator.real.ravel(),
                propagator.imag.ravel(),
            )
        )

    def unpack(flat):
        amplitude = (
            flat[:retained_size] + 1j * flat[retained_size : 2 * retained_size]
        ).reshape(retained, initial_retained.shape[1])
        rest = flat[2 * retained_size :]
        propagator = (rest[:propagator_size] + 1j * rest[propagator_size:]).reshape(exterior, exterior)
        return amplitude, propagator

    def derivative(time, flat):
        amplitude, propagator = unpack(flat)
        local, coupling, exterior_block = hamiltonian_blocks(hamiltonian(time), frame, retained)
        amplitude_dot = -1j * (local @ amplitude) - 1j * (coupling @ (propagator @ initial_exterior))
        propagator_dot = -1j * (exterior_block @ propagator)
        return pack(amplitude_dot, propagator_dot)

    solution = solve_ivp(
        derivative,
        (float(times[0]), float(times[-1])),
        pack(initial_retained, initial_propagator),
        t_eval=np.asarray(times, dtype=float),
        method="DOP853",
        rtol=1e-11,
        atol=1e-13,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return np.stack([unpack(solution.y[:, index])[0] for index in range(solution.y.shape[1])])


def _project(samples, basis):
    return np.einsum("ij,tjk->tik", basis.conj().T, samples)


def _max_frobenius(difference):
    return float(np.linalg.norm(difference, axis=(1, 2)).max())


def _history_deviation(reduced):
    step = float(reduced["times"][1] - reduced["times"][0])
    accumulator = np.zeros_like(reduced["history"][0])
    deviation = 0.0
    coupling = reduced["coupling_block"]
    propagator = reduced["exterior_propagator"]
    amplitudes = reduced["amplitudes"]
    for index in range(reduced["times"].size - 1):
        kernel_now = propagator[index].conj().T @ coupling[index].conj().T
        kernel_new = propagator[index + 1].conj().T @ coupling[index + 1].conj().T
        accumulator = accumulator + (step / 2.0) * (
            kernel_now @ amplitudes[index] + kernel_new @ amplitudes[index + 1]
        )
        deviation = max(
            deviation,
            float(np.linalg.norm(accumulator - reduced["history"][index + 1], ord="fro")),
        )
    return deviation


def _richardson_coefficient(delta_half, delta_final, half_time, final_time):
    earlier = delta_half / (half_time ** 2)
    later = delta_final / (final_time ** 2)
    return 2.0 * earlier - later


@lru_cache(maxsize=1)
def measurements():
    built = example()
    hamiltonian = built["hamiltonian"]
    basis = built["basis"]
    frame = built["frame"]
    covariance = built["covariance"]
    identity = np.eye(_DIMENSION, dtype=np.complex128)
    amplitude_error = {}
    amplitude_final = {}
    covariance_error = {}
    covariance_final = {}
    car = {}
    trapezoid = {}
    assembly = {}
    fine = None
    fine_projected = None
    fine_oracle = None
    for steps in _GRIDS:
        times = np.linspace(0.0, _FINAL, steps + 1)
        reduced = evolve_retained_region(
            hamiltonian,
            basis,
            times,
            identity,
            covariance=covariance,
        )
        samples = integrate_full_system(hamiltonian, times, identity)
        projected = _project(samples, basis)
        oracle = np.einsum("tia,ab,tjb->tij", projected, covariance, projected.conj())
        amplitude_error[steps] = _max_frobenius(reduced["amplitudes"] - projected)
        amplitude_final[steps] = float(np.linalg.norm(reduced["amplitudes"][-1] - projected[-1], ord="fro"))
        covariance_error[steps] = _max_frobenius(reduced["covariance_total"] - oracle)
        covariance_final[steps] = float(
            np.linalg.norm(reduced["covariance_total"][-1] - oracle[-1], ord="fro")
        )
        car[steps] = float(reduced["car_deviation_max"])
        trapezoid[steps] = float(reduced["trapezoid_residual_max"])
        assembly[steps] = float(reduced["covariance_assembly_residual_max"])
        if steps == _GRIDS[-1]:
            fine = reduced
            fine_projected = projected
            fine_oracle = oracle

    retained_column = basis[:, :1]
    exterior_column = frame[:, _RETAINED : _RETAINED + 1]
    fine_times = fine["times"]
    driven = evolve_retained_region(hamiltonian, basis, fine_times, exterior_column)
    driven_samples = integrate_full_system(hamiltonian, fine_times, exterior_column)
    driven_projected = _project(driven_samples, basis)
    drive_error = _max_frobenius(driven["amplitudes"] - driven_projected)
    silent = evolve_retained_region(
        hamiltonian,
        basis,
        fine_times,
        exterior_column,
        outside_drive=False,
    )
    forgotten = evolve_retained_region(
        hamiltonian,
        basis,
        fine_times,
        identity,
        memory=False,
    )
    cross_off = evolve_retained_region(
        hamiltonian,
        basis,
        fine_times,
        identity,
        covariance=covariance,
        drop_cross_covariance=True,
    )
    blocks = (
        fine["frame_blocks"]["parent"],
        fine["frame_blocks"]["parent_child"],
        fine["frame_blocks"]["child_parent"],
        fine["frame_blocks"]["child"],
    )
    initial_terms = covariance_terms(fine["parent_response"][0], fine["child_response"][0], blocks)
    final_terms = covariance_terms(fine["parent_response"][-1], fine["child_response"][-1], blocks)
    transformed = frame.conj().T @ covariance @ frame
    child = transformed[_RETAINED:, _RETAINED:]
    child_parent = transformed[_RETAINED:, :_RETAINED]
    later = fine_times.size - 1
    earlier = fine_times.size // 2
    produced = exterior_kernels(
        float(fine_times[later]),
        float(fine_times[earlier]),
        coupling_t=fine["coupling_block"][later],
        propagator_t=fine["exterior_propagator"][later],
        coupling_s=fine["coupling_block"][earlier],
        propagator_s=fine["exterior_propagator"][earlier],
        exterior_covariance=child,
        initial_cross=child_parent,
    )
    source_t = fine["coupling_block"][later] @ fine["exterior_propagator"][later]
    source_s = fine["coupling_block"][earlier] @ fine["exterior_propagator"][earlier]
    manual_occupied = source_t @ child @ source_s.conj().T
    manual_empty = source_t @ (np.eye(child.shape[0], dtype=np.complex128) - child) @ source_s.conj().T
    mask = np.array([True, True, False, False])
    reference = GaussianBoundaryState(transformed, mask).kernels(
        float(fine_times[later]),
        float(fine_times[earlier]),
        link_t=fine["coupling_block"][later],
        child_t=fine["exterior_propagator"][later],
        link_s=fine["coupling_block"][earlier],
        child_s=fine["exterior_propagator"][earlier],
    )
    same = exterior_kernels(
        float(fine_times[later]),
        float(fine_times[later]),
        coupling_t=fine["coupling_block"][later],
        propagator_t=fine["exterior_propagator"][later],
        coupling_s=fine["coupling_block"][later],
        propagator_s=fine["exterior_propagator"][later],
        exterior_covariance=child,
        initial_cross=child_parent,
    )
    reversed_kernel = exterior_kernels(
        float(fine_times[earlier]),
        float(fine_times[later]),
        coupling_t=fine["coupling_block"][earlier],
        propagator_t=fine["exterior_propagator"][earlier],
        coupling_s=fine["coupling_block"][later],
        propagator_s=fine["exterior_propagator"][later],
        exterior_covariance=child,
        initial_cross=child_parent,
    )

    short = np.linspace(0.0, 0.01, 81)
    half_index = short.size // 2
    local, coupling, _exterior = hamiltonian_blocks(hamiltonian(0.0), frame, _RETAINED)
    step = 1e-6
    local_later = hamiltonian_blocks(hamiltonian(step), frame, _RETAINED)[0]
    local_earlier = hamiltonian_blocks(hamiltonian(-step), frame, _RETAINED)[0]
    local_dot = (local_later - local_earlier) / (2.0 * step)
    coefficients = frame.conj().T @ retained_column
    initial_retained = coefficients[:_RETAINED]
    initial_exterior = coefficients[_RETAINED:]
    initial_velocity = -1j * (local @ initial_retained)
    memory_term = -(coupling @ (coupling.conj().T @ initial_retained))
    without_memory = -1j * (local_dot @ initial_retained) - 1j * (local @ initial_velocity)
    second_minus = without_memory + memory_term
    second_plus = without_memory - memory_term
    reduced_short = evolve_retained_region(hamiltonian, basis, short, retained_column)
    reduced_silent_memory = evolve_retained_region(
        hamiltonian,
        basis,
        short,
        retained_column,
        memory=False,
    )
    full_short = _project(integrate_full_system(hamiltonian, short, retained_column), basis)
    free_short = integrate_memory_free(hamiltonian, frame, short, retained_column)
    delta_full = full_short - free_short
    delta_reducer = reduced_short["amplitudes"] - reduced_silent_memory["amplitudes"]
    coefficient_full = _richardson_coefficient(
        delta_full[half_index],
        delta_full[-1],
        float(short[half_index]),
        float(short[-1]),
    )
    coefficient_reducer = _richardson_coefficient(
        delta_reducer[half_index],
        delta_reducer[-1],
        float(short[half_index]),
        float(short[-1]),
    )
    predicted_coefficient = 0.5 * memory_term

    def coefficient_errors(extracted):
        minus = float(np.linalg.norm(extracted - predicted_coefficient, ord="fro"))
        plus = float(np.linalg.norm(extracted + predicted_coefficient, ord="fro"))
        return minus, plus

    full_minus, full_plus = coefficient_errors(coefficient_full)
    reducer_minus, reducer_plus = coefficient_errors(coefficient_reducer)
    taylor = (
        initial_retained
        + float(short[-1]) * initial_velocity
        + (float(short[-1]) ** 2) / 2.0 * second_minus
    )
    taylor_flipped = (
        initial_retained
        + float(short[-1]) * initial_velocity
        + (float(short[-1]) ** 2) / 2.0 * second_plus
    )
    drive_short_times = np.linspace(0.0, 0.01, 81)
    drive_short = evolve_retained_region(hamiltonian, basis, drive_short_times, exterior_column)
    drive_short_full = _project(
        integrate_full_system(hamiltonian, drive_short_times, exterior_column),
        basis,
    )
    exterior_coefficients = frame.conj().T @ exterior_column
    drive_velocity = -1j * (coupling @ exterior_coefficients[_RETAINED:])
    drive_scale = drive_short["amplitudes"][-1] / float(drive_short_times[-1])
    propagator_defect = 0.0
    identity_exterior = np.eye(fine["exterior_propagator"].shape[-1], dtype=np.complex128)
    for propagator in fine["exterior_propagator"]:
        propagator_defect = max(
            propagator_defect,
            float(np.linalg.norm(propagator.conj().T @ propagator - identity_exterior, ord="fro")),
        )
    kernel_gap = {
        name: float(np.linalg.norm(produced[name] - reference[name], ord="fro"))
        for name in produced
    }
    return {
        "amplitude_error_max": amplitude_error,
        "amplitude_error_final": amplitude_final,
        "covariance_error_max": covariance_error,
        "covariance_error_final": covariance_final,
        "car_deviation_max": car,
        "trapezoid_residual_max": trapezoid,
        "assembly_residual_max": assembly,
        "history_trapezoid_deviation": _history_deviation(fine),
        "history_initial": float(np.linalg.norm(fine["history"][0], ord="fro")),
        "history_final": float(np.linalg.norm(fine["history"][-1], ord="fro")),
        "omitted_memory_history": float(np.linalg.norm(forgotten["history"])),
        "omitted_memory_amplitude": _max_frobenius(forgotten["amplitudes"] - fine_projected),
        "memory_separation": _max_frobenius(fine["amplitudes"] - forgotten["amplitudes"]),
        "cross_initial": float(np.linalg.norm(initial_terms["initial_parent_child"], ord="fro")),
        "cross_initial_other": float(np.linalg.norm(initial_terms["initial_child_parent"], ord="fro")),
        "cross_final": float(np.linalg.norm(final_terms["initial_parent_child"], ord="fro")),
        "cross_final_other": float(np.linalg.norm(final_terms["initial_child_parent"], ord="fro")),
        "cross_drop_amplitude": _max_frobenius(fine["amplitudes"] - cross_off["amplitudes"]),
        "cross_drop_covariance": float(
            np.linalg.norm(fine["covariance_total"][-1] - cross_off["covariance_total"][-1], ord="fro")
        ),
        "cross_matches_drop": float(
            np.linalg.norm(
                (fine["covariance_total"][-1] - cross_off["covariance_total"][-1])
                - (final_terms["initial_parent_child"] + final_terms["initial_child_parent"]),
                ord="fro",
            )
        ),
        "initial_covariance_defect": float(
            np.linalg.norm(fine["covariance_total"][0] - blocks[0], ord="fro")
        ),
        "kernel_occupied": float(np.linalg.norm(produced["occupied"], ord="fro")),
        "kernel_empty": float(np.linalg.norm(produced["empty"], ord="fro")),
        "kernel_gap": float(np.linalg.norm(produced["occupied"] - produced["empty"], ord="fro")),
        "kernel_occupied_formula": float(np.linalg.norm(produced["occupied"] - manual_occupied, ord="fro")),
        "kernel_empty_formula": float(np.linalg.norm(produced["empty"] - manual_empty, ord="fro")),
        "kernel_boundary_state": kernel_gap,
        "kernel_equal_time_retarded": float(
            np.linalg.norm(same["retarded"] + 0.5j * (source_t @ source_t.conj().T), ord="fro")
        ),
        "kernel_reversed_retarded": float(np.linalg.norm(reversed_kernel["retarded"], ord="fro")),
        "kernel_not_added": float(
            np.linalg.norm(fine["covariance_total"][-1] + produced["occupied"] - fine_oracle[-1], ord="fro")
        ),
        "drive_error": drive_error,
        "drive_omitted_norm": float(np.linalg.norm(silent["amplitudes"])),
        "drive_omitted_history": float(np.linalg.norm(silent["history"])),
        "drive_recorded_exterior": float(np.linalg.norm(silent["initial_exterior"], ord="fro")),
        "drive_amplitude": float(np.linalg.norm(driven["amplitudes"][-1], ord="fro")),
        "drive_short_error": _max_frobenius(drive_short["amplitudes"] - drive_short_full),
        "drive_leading_relative": float(
            np.linalg.norm(drive_scale - drive_velocity, ord="fro")
            / np.linalg.norm(drive_velocity, ord="fro")
        ),
        "drive_flipped_relative": float(
            np.linalg.norm(drive_scale + drive_velocity, ord="fro")
            / np.linalg.norm(drive_velocity, ord="fro")
        ),
        "memory_term_norm": float(np.linalg.norm(memory_term, ord="fro")),
        "memory_coefficient_full_minus": full_minus,
        "memory_coefficient_full_plus": full_plus,
        "memory_coefficient_reducer_minus": reducer_minus,
        "memory_coefficient_reducer_plus": reducer_plus,
        "memory_coefficient_reducer_vs_full": float(
            np.linalg.norm(coefficient_reducer - coefficient_full, ord="fro")
        ),
        "taylor_error": float(np.linalg.norm(full_short[-1] - taylor, ord="fro")),
        "taylor_flipped_error": float(np.linalg.norm(full_short[-1] - taylor_flipped, ord="fro")),
        "exterior_propagator_defect": propagator_defect,
        "exterior_dimension": int(fine["exterior_propagator"].shape[-1]),
        "amplitude_shape": list(fine["amplitudes"].shape),
        "initial_exterior_norm_for_memory_column": float(np.linalg.norm(initial_exterior, ord="fro")),
        "basis_imaginary_norm": float(np.linalg.norm(basis.imag)),
        "basis_preserved": bool(
            np.array_equal(fine["local_basis"], basis) and np.array_equal(fine["frame"][:, :_RETAINED], basis)
        ),
    }


def test_fixed_complex_basis_is_preserved():
    built = example()
    basis = built["basis"]
    data = measurements()
    assert basis.shape == (_DIMENSION, _RETAINED)
    assert data["basis_imaginary_norm"] > 0.1
    for column in range(basis.shape[1]):
        assert np.count_nonzero(np.abs(basis[:, column]) > 1e-2) >= 3
    # A unitary matrix is already a QR factor, so qr(basis) cannot show rephasing.
    # Forcing each column's first large entry onto the positive real axis does.
    rephased = np.array(basis, copy=True)
    for column in range(rephased.shape[1]):
        pivot = int(np.flatnonzero(np.abs(rephased[:, column]) > 1e-8)[0])
        rephased[:, column] *= np.exp(-1j * np.angle(rephased[pivot, column]))
    assert not np.allclose(rephased, basis)
    assert "qr(" not in Path(reduction.__file__).read_text(encoding="utf-8")
    assert data["basis_preserved"]
    assert np.array_equal(built["frame"][:, :_RETAINED], basis)
    completed = fixed_observer_frame(basis)
    assert np.array_equal(completed["local_basis"], basis)
    assert not np.allclose(completed["local_basis"], rephased)
    assert np.array_equal(completed["frame"][:, :_RETAINED], basis)


def test_reducer_uses_exterior_history_not_full_ode(monkeypatch):
    source = Path(reduction.__file__).read_text(encoding="utf-8")
    assert "solve_ivp" not in source
    assert "DOP853" not in source
    built = example()
    times = np.linspace(0.0, 0.4, 5)
    seen_times = []
    shapes = []
    original = reduction.np.linalg.solve

    def wrapped(time):
        seen_times.append(float(time))
        return built["hamiltonian"](time)

    def spy(matrix, right_hand):
        shapes.append(matrix.shape)
        return original(matrix, right_hand)

    monkeypatch.setattr(reduction.np.linalg, "solve", spy)
    reduced = evolve_retained_region(wrapped, built["basis"], times, built["basis"][:, :1])
    assert shapes
    assert all(shape == (_RETAINED, _RETAINED) for shape in shapes)
    assert len(shapes) == times.size - 1
    grid = np.asarray(times, dtype=float)
    stencil = list(grid)
    stencil.extend(0.5 * (grid[:-1] + grid[1:]))
    assert len(seen_times) == 2 * grid.size
    for sample in seen_times:
        assert min(abs(sample - point) for point in stencil) < 1e-9
    for point in stencil:
        assert min(abs(sample - point) for sample in seen_times) < 1e-9
    data = measurements()
    assert data["exterior_dimension"] == _DIMENSION - _RETAINED
    assert data["amplitude_shape"][1] == _RETAINED
    assert data["history_initial"] == 0.0
    assert data["history_final"] > 1e-3
    assert data["history_trapezoid_deviation"] < 1e-12
    assert data["omitted_memory_history"] == 0.0
    assert data["memory_separation"] > 20.0 * data["amplitude_error_max"][_GRIDS[-1]]
    assert reduced["flags"]["history_integrator"] == "uniform-trapezoid-semiseparable"
    assert reduced["flags"]["exterior_integrator"] == "midpoint-expm"
    assert data["exterior_propagator_defect"] < 1e-8


def test_retained_amplitudes_and_covariance_match_independent_evolution():
    data = measurements()
    finest = _GRIDS[-1]
    assert data["amplitude_error_max"][finest] < 1e-3
    assert data["covariance_error_max"][finest] < 1e-3
    assert data["amplitude_error_final"][finest] <= data["amplitude_error_max"][finest]
    assert data["covariance_error_final"][finest] <= data["covariance_error_max"][finest]
    assert data["trapezoid_residual_max"][finest] < 1e-12
    assert data["assembly_residual_max"][finest] < 1e-12
    assert data["memory_separation"] > 20.0 * data["amplitude_error_max"][finest]
    assert data["omitted_memory_amplitude"] > 20.0 * data["amplitude_error_max"][finest]


def test_refinement_converges():
    data = measurements()
    for label, series in (
        ("amplitude", data["amplitude_error_max"]),
        ("covariance", data["covariance_error_max"]),
        ("car", data["car_deviation_max"]),
    ):
        coarse, middle, fine = (series[steps] for steps in _GRIDS)
        assert coarse > middle > fine, (label, coarse, middle, fine)
        assert coarse / middle > 3.5, (label, coarse / middle)
        assert middle / fine > 3.5, (label, middle / fine)


def test_initial_cross_covariance():
    data = measurements()
    assert data["cross_initial"] < 1e-12
    assert data["cross_initial_other"] < 1e-12
    assert data["initial_covariance_defect"] < 1e-12
    assert data["cross_final"] > 1e-4
    assert data["cross_final_other"] > 1e-4
    assert data["cross_drop_amplitude"] < 1e-12
    assert data["cross_matches_drop"] < 1e-12
    assert data["cross_drop_covariance"] > 50.0 * data["covariance_error_final"][_GRIDS[-1]]


def test_occupied_and_empty_exterior_kernels():
    data = measurements()
    assert data["kernel_occupied"] > 1e-3
    assert data["kernel_empty"] > 1e-3
    assert data["kernel_gap"] > 1e-3
    assert data["kernel_occupied_formula"] < 1e-12
    assert data["kernel_empty_formula"] < 1e-12
    assert max(data["kernel_boundary_state"].values()) < 1e-10
    assert data["kernel_equal_time_retarded"] < 1e-12
    assert data["kernel_reversed_retarded"] < 1e-15
    assert data["kernel_not_added"] > 1e-3
    assert data["assembly_residual_max"][_GRIDS[-1]] < 1e-12


def test_nonzero_exterior_drive():
    data = measurements()
    assert data["drive_recorded_exterior"] > 0.5
    assert data["drive_omitted_norm"] < 1e-9
    assert data["drive_omitted_history"] < 1e-9
    assert data["drive_amplitude"] > 1e-2
    assert data["drive_error"] < 1e-3
    assert data["drive_error"] * 20.0 < data["drive_amplitude"]
    assert data["drive_short_error"] < 1e-6
    assert data["drive_leading_relative"] < 0.2
    assert data["drive_flipped_relative"] > 1.5


def test_omitting_the_coupling_keeps_preparation_and_matches_decoupled_evolution():
    """Zeroing B is not the same as dropping memory, the exterior drive, or C_AE.

    The same frame, covariance, and initial columns are passed. An independent
    full system with that block removed is the comparison.
    """
    built = example()
    basis = built["basis"]
    frame = built["frame"]
    hamiltonian = built["hamiltonian"]
    times = np.linspace(0.0, _FINAL, _GRIDS[1] + 1)
    identity = np.eye(_DIMENSION, dtype=np.complex128)
    coupled = evolve_retained_region(
        hamiltonian,
        basis,
        times,
        identity,
        covariance=built["covariance"],
    )
    uncoupled = evolve_retained_region(
        hamiltonian,
        basis,
        times,
        identity,
        covariance=built["covariance"],
        coupling=False,
    )
    forgotten = evolve_retained_region(
        hamiltonian,
        basis,
        times,
        identity,
        covariance=built["covariance"],
        memory=False,
    )
    assert np.array_equal(uncoupled["local_basis"], basis)
    assert np.allclose(uncoupled["initial_exterior"], coupled["initial_exterior"])
    assert np.allclose(uncoupled["frame_blocks"]["parent"], coupled["frame_blocks"]["parent"])
    assert np.allclose(
        uncoupled["frame_blocks"]["parent_child"],
        coupled["frame_blocks"]["parent_child"],
    )
    assert float(np.linalg.norm(uncoupled["frame_blocks"]["parent_child"])) > 1e-3
    assert float(np.linalg.norm(uncoupled["coupling_block"])) == 0.0
    assert float(np.linalg.norm(coupled["coupling_block"])) > 0.1
    assert float(np.linalg.norm(uncoupled["initial_exterior"])) > 0.1
    assert float(np.linalg.norm(uncoupled["history"])) == 0.0

    def decoupled(time):
        local, _coupling, exterior = hamiltonian_blocks(hamiltonian(time), frame, _RETAINED)
        zeros = np.zeros((_RETAINED, _DIMENSION - _RETAINED), dtype=np.complex128)
        return frame @ _pack_blocks(local, zeros, exterior) @ frame.conj().T

    projected = _project(integrate_full_system(decoupled, times, identity), basis)
    error = _max_frobenius(uncoupled["amplitudes"] - projected)
    separation = _max_frobenius(coupled["amplitudes"] - uncoupled["amplitudes"])
    memory_gap = _max_frobenius(forgotten["amplitudes"] - uncoupled["amplitudes"])
    assert error < 1e-3
    assert separation > 20.0 * error
    assert memory_gap > 20.0 * error


def test_memory_sign_from_short_time_expansion():
    """Leading memory piece of ``X(t)`` is ``-(t**2)/2 B(0) B(0)† X(0)``.

    That coefficient is read from the independent full evolution minus a
    memory-free retained integration, then compared with the reducer.
    """
    data = measurements()
    assert data["initial_exterior_norm_for_memory_column"] < 1e-10
    assert data["memory_term_norm"] > 1e-2
    assert data["memory_coefficient_full_minus"] * 20.0 < data["memory_coefficient_full_plus"]
    assert data["memory_coefficient_reducer_minus"] * 20.0 < data["memory_coefficient_reducer_plus"]
    assert data["memory_coefficient_reducer_vs_full"] < 0.05 * data["memory_term_norm"]
    assert data["taylor_error"] * 5.0 < data["taylor_flipped_error"]


if __name__ == "__main__":
    print(json.dumps(measurements(), indent=2, sort_keys=True, default=float))
