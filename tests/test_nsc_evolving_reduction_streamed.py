"""Streamed retained-region backend against the dense Volterra control.

The generator is Hermitian. Distinct times do not commute. Full-system
samples use DOP853 outside the reducer. No coupled H(g(t)) is built, and
the incoming-gate campaign is not run.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.sparse.linalg import LinearOperator

from recursive_horizons.nsc_evolving_reduction import (
    causal_memory_kernel,
    evolve_retained_region,
    exterior_column_transport,
    exterior_kernels,
    exterior_kernels_from_cohorts,
    hamiltonian_blocks,
    midpoint_exponential_action,
    prescribed_six_mode_problem,
)
import recursive_horizons.nsc_evolving_reduction as reduction


_GRIDS = (16, 32, 64)
_FINAL = 1.0
_SUBSTEPS = 2


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


def _project(samples, basis):
    return np.einsum("ij,tjk->tik", basis.conj().T, samples)


def _max_norm(difference):
    return float(np.linalg.norm(difference, axis=(1, 2)).max())


def _pair(hamiltonian, basis, times, columns, *, substeps=_SUBSTEPS, **kwargs):
    dense = evolve_retained_region(
        hamiltonian,
        basis,
        times,
        columns,
        propagator_substeps=substeps,
        **kwargs,
    )
    streamed = evolve_retained_region(
        hamiltonian,
        basis,
        times,
        columns,
        propagator_substeps=substeps,
        backend="streamed",
        **kwargs,
    )
    return dense, streamed


@pytest.fixture(scope="module")
def prescribed():
    problem = prescribed_six_mode_problem()
    frame = evolve_retained_region(
        problem["hamiltonian"],
        problem["local_basis"],
        np.linspace(0.0, 0.1, 3),
        np.eye(6, dtype=np.complex128),
    )["frame"]
    exterior_samples = [
        hamiltonian_blocks(problem["hamiltonian"](time), frame, 2)[2]
        for time in (0.0, 0.7)
    ]
    commutator = exterior_samples[0] @ exterior_samples[1] - exterior_samples[1] @ exterior_samples[0]
    problem = dict(problem)
    problem["frame_probe"] = frame
    problem["exterior_commutator"] = float(np.linalg.norm(commutator, ord="fro"))
    problem["exterior_hermitian"] = float(
        np.linalg.norm(exterior_samples[0] - exterior_samples[0].conj().T, ord="fro")
    )
    return problem


@pytest.fixture(scope="module")
def comparison(prescribed):
    hamiltonian = prescribed["hamiltonian"]
    basis = prescribed["local_basis"]
    covariance = prescribed["covariance"]
    identity = np.eye(6, dtype=np.complex128)
    amplitude_error = {}
    covariance_error = {}
    dense_gap = {}
    fine = None
    fine_stream = None
    fine_projected = None
    for steps in _GRIDS:
        times = np.linspace(0.0, _FINAL, steps + 1)
        dense, streamed = _pair(
            hamiltonian,
            basis,
            times,
            identity,
            covariance=covariance,
        )
        samples = integrate_full_system(hamiltonian, times, identity)
        projected = _project(samples, basis)
        oracle = np.einsum("tia,ab,tjb->tij", projected, covariance, projected.conj())
        amplitude_error[steps] = _max_norm(streamed["amplitudes"] - projected)
        covariance_error[steps] = _max_norm(streamed["covariance_total"] - oracle)
        dense_gap[steps] = max(
            _max_norm(streamed["amplitudes"] - dense["amplitudes"]),
            _max_norm(streamed["covariance_total"] - dense["covariance_total"]),
            _max_norm(streamed["history"] - dense["history"]),
        )
        if steps == _GRIDS[-1]:
            fine = dense
            fine_stream = streamed
            fine_projected = projected
    times = fine["times"]
    omitted = {}
    for name, flags in (
        ("memory", {"memory": False}),
        ("drive", {"outside_drive": False}),
        ("coupling", {"coupling": False}),
    ):
        dense_flag, stream_flag = _pair(hamiltonian, basis, times, identity, **flags)
        omitted[name] = {
            "gap": _max_norm(dense_flag["amplitudes"] - stream_flag["amplitudes"]),
            "dense_separation": _max_norm(fine["amplitudes"] - dense_flag["amplitudes"]),
            "stream_separation": _max_norm(fine_stream["amplitudes"] - stream_flag["amplitudes"]),
            "history": float(np.linalg.norm(stream_flag["history"])),
            "initial_exterior": float(np.linalg.norm(stream_flag["initial_exterior"], ord="fro")),
            "required_columns": int(stream_flag["allocation"]["required_exterior_columns"]),
        }
    dense_cross, stream_cross = _pair(
        hamiltonian,
        basis,
        times,
        identity,
        covariance=covariance,
        drop_cross_covariance=True,
    )
    cross_dense = fine["covariance_total"] - dense_cross["covariance_total"]
    cross_stream = fine_stream["covariance_total"] - stream_cross["covariance_total"]
    later = times.size - 1
    earlier = times.size // 3
    assembled = causal_memory_kernel(
        hamiltonian,
        fine["frame"],
        2,
        times,
        later,
        earlier,
        substeps=_SUBSTEPS,
    )
    dense_kernel = (
        fine["coupling_block"][later]
        @ fine["exterior_propagator"][later]
        @ fine["exterior_propagator"][earlier].conj().T
        @ fine["coupling_block"][earlier].conj().T
    )
    cohort_kernels = exterior_kernels_from_cohorts(
        float(times[later]),
        float(times[earlier]),
        assembled["cohort_later"],
        assembled["cohort_earlier"],
        fine["frame_blocks"]["child"],
        fine["frame_blocks"]["child_parent"],
    )
    dense_kernels = exterior_kernels(
        float(times[later]),
        float(times[earlier]),
        coupling_t=fine["coupling_block"][later],
        propagator_t=fine["exterior_propagator"][later],
        coupling_s=fine["coupling_block"][earlier],
        propagator_s=fine["exterior_propagator"][earlier],
        exterior_covariance=fine["frame_blocks"]["child"],
        initial_cross=fine["frame_blocks"]["child_parent"],
    )

    def decoupled(time):
        local, _coupling, exterior = hamiltonian_blocks(hamiltonian(time), fine["frame"], 2)
        blocked = np.zeros((6, 6), dtype=np.complex128)
        blocked[:2, :2] = local
        blocked[2:, 2:] = exterior
        return fine["frame"] @ blocked @ fine["frame"].conj().T

    coupling_off = omitted["coupling"]
    decoupled_samples = _project(integrate_full_system(decoupled, times, identity), basis)
    # The coupling-off streamed amplitudes are recovered by a fresh call so the
    # decoupled oracle uses the same retained samples as the omission pair.
    _dense_off, stream_off = _pair(hamiltonian, basis, times, identity, coupling=False)
    return {
        "amplitude_error": amplitude_error,
        "covariance_error": covariance_error,
        "dense_gap": dense_gap,
        "fine": fine,
        "fine_stream": fine_stream,
        "fine_projected": fine_projected,
        "omitted": omitted,
        "cross_amplitude": _max_norm(stream_cross["amplitudes"] - fine_stream["amplitudes"]),
        "cross_resolution": float(np.linalg.norm(cross_dense - cross_stream)),
        "cross_signal": float(np.linalg.norm(cross_stream[-1], ord="fro")),
        "kernel_gap": float(np.linalg.norm(assembled["kernel"] - dense_kernel, ord="fro")),
        "kernel_norm": float(np.linalg.norm(dense_kernel, ord="fro")),
        "kernel_blocks": {
            name: float(np.linalg.norm(cohort_kernels[name] - dense_kernels[name], ord="fro"))
            for name in dense_kernels
        },
        "occupied_empty_gap": float(
            np.linalg.norm(cohort_kernels["occupied"] - cohort_kernels["empty"], ord="fro")
        ),
        "initial_cross_norm": float(np.linalg.norm(cohort_kernels["initial_occupied"], ord="fro")),
        "decoupled_error": _max_norm(stream_off["amplitudes"] - decoupled_samples),
        "coupling_separation": coupling_off["stream_separation"],
        "stores_kernel_propagator": assembled["stores_exterior_propagator"],
    }


def test_module_does_not_integrate_the_full_system():
    source = Path(reduction.__file__).read_text(encoding="utf-8")
    assert "solve_ivp" not in source
    assert "DOP853" not in source
    assert "expm_multiply" in source


def test_default_backend_remains_the_dense_propagator(prescribed):
    times = np.linspace(0.0, 0.3, 5)
    identity = np.eye(6, dtype=np.complex128)
    default = evolve_retained_region(
        prescribed["hamiltonian"],
        prescribed["local_basis"],
        times,
        identity,
    )
    named = evolve_retained_region(
        prescribed["hamiltonian"],
        prescribed["local_basis"],
        times,
        identity,
        backend="dense",
    )
    assert default["flags"]["backend"] == "dense"
    assert default["flags"]["history_integrator"] == "uniform-trapezoid-semiseparable"
    assert default["flags"]["exterior_integrator"] == "midpoint-expm"
    assert default["flags"]["stores_dense_exterior_propagator"] is True
    assert default["exterior_propagator"].shape == (5, 4, 4)
    assert default["allocation"]["time_indexed_exterior_propagator_bytes"] == (
        5 * 4 * 4 * 16
    )
    assert np.array_equal(default["amplitudes"], named["amplitudes"])
    assert np.array_equal(default["exterior_propagator"], named["exterior_propagator"])
    with pytest.raises(ValueError, match="backend"):
        evolve_retained_region(
            prescribed["hamiltonian"],
            prescribed["local_basis"],
            times,
            identity,
            backend="full",
        )


def test_generator_is_hermitian_and_noncommuting(prescribed):
    covariance = prescribed["covariance"]
    off_diagonal = covariance - np.diag(np.diag(covariance))
    eigenvalues = np.linalg.eigvalsh(covariance)
    assert prescribed["exterior_hermitian"] < 1e-12
    assert prescribed["exterior_commutator"] > 1e-2
    assert float(np.linalg.norm(off_diagonal)) > 1e-3
    assert float(eigenvalues.min()) > 0.0
    assert float(eigenvalues.max()) < 1.0
    link = hamiltonian_blocks(prescribed["hamiltonian"](0.4), prescribed["frame_probe"], 2)[1]
    assert float(np.linalg.norm(link.imag)) > 1e-4


def test_observer_phases_and_column_api_match(prescribed):
    generator = np.random.default_rng(3)
    draw = generator.normal(size=(6, 2)) + 1j * generator.normal(size=(6, 2))
    orthogonal, _triangle = np.linalg.qr(draw)
    basis = orthogonal * np.exp(1j * np.array([0.4, -1.1]))
    times = np.linspace(0.0, 0.5, 11)
    dense, streamed = _pair(
        prescribed["hamiltonian"],
        basis,
        times,
        np.eye(6, dtype=np.complex128),
        covariance=prescribed["covariance"],
        substeps=3,
    )
    assert np.array_equal(streamed["local_basis"], basis)
    assert np.array_equal(streamed["frame"][:, :2], basis)
    assert streamed["exterior_propagator"] is None
    assert np.allclose(streamed["amplitudes"], dense["amplitudes"], atol=1e-12, rtol=0)
    assert np.allclose(streamed["history"], dense["history"], atol=1e-12, rtol=0)
    assert np.allclose(streamed["covariance_total"], dense["covariance_total"], atol=1e-12, rtol=0)
    assert streamed["flags"]["memory_kernel"] == "B(t) U_E(t,s) B(s)†"
    assert streamed["flags"]["exterior_integrator"] == "midpoint-expm-multiply"
    raw_columns = generator.normal(size=(6, 3)) + 1j * generator.normal(size=(6, 3))
    short_columns, _triangle = np.linalg.qr(raw_columns)
    weights = np.array([0.2, 0.5, 0.8])
    weighted_dense, weighted_stream = _pair(
        prescribed["hamiltonian"],
        basis,
        times,
        short_columns,
        weights=weights,
        substeps=1,
    )
    assert np.allclose(
        weighted_stream["column_weight_covariance"],
        weighted_dense["column_weight_covariance"],
        atol=1e-12,
        rtol=0,
    )


def test_dense_streamed_and_dop853_controls(comparison):
    gap = comparison["dense_gap"]
    amplitude = comparison["amplitude_error"]
    covariance = comparison["covariance_error"]
    assert max(gap.values()) < 1e-12
    assert amplitude[16] > amplitude[32] > amplitude[64]
    assert covariance[16] > covariance[32] > covariance[64]
    assert amplitude[16] / amplitude[32] > 3.0
    assert amplitude[32] / amplitude[64] > 3.0
    assert covariance[16] / covariance[32] > 3.0
    assert amplitude[64] < 1e-3
    assert covariance[64] < 1e-4
    assert gap[64] * 1e6 < amplitude[64]
    fine = comparison["fine_stream"]
    assert fine["trapezoid_residual_max"] < 1e-12
    assert fine["covariance_assembly_residual_max"] < 1e-12
    assert float(np.linalg.norm(fine["history"][0])) == 0.0
    assert float(np.linalg.norm(fine["history"][-1])) > 1e-2
    memory = comparison["omitted"]["memory"]
    drive = comparison["omitted"]["drive"]
    coupling = comparison["omitted"]["coupling"]
    assert memory["gap"] < 1e-12
    assert drive["gap"] < 1e-12
    assert coupling["gap"] < 1e-12
    assert memory["history"] == 0.0
    assert memory["initial_exterior"] > 0.1
    assert drive["initial_exterior"] > 0.1
    assert drive["history"] > 1e-2
    assert coupling["history"] == 0.0
    assert coupling["required_columns"] == 0
    assert memory["required_columns"] > 0
    finest = amplitude[64]
    for item in (memory, drive, coupling):
        assert item["stream_separation"] > 20.0 * finest
        assert item["dense_separation"] > 20.0 * finest
        assert abs(item["stream_separation"] - item["dense_separation"]) < 1e-12
    assert comparison["cross_amplitude"] < 1e-12
    assert comparison["cross_resolution"] < 1e-12
    assert comparison["cross_signal"] > 50.0 * covariance[64]
    assert comparison["decoupled_error"] < 3.0 * finest
    assert comparison["coupling_separation"] > 20.0 * comparison["decoupled_error"]


def test_causal_kernel_uses_source_cohorts_not_a_surrogate(prescribed, comparison):
    assert comparison["kernel_gap"] < 1e-12
    assert comparison["kernel_norm"] > 1e-3
    assert max(comparison["kernel_blocks"].values()) < 1e-12
    assert comparison["occupied_empty_gap"] > 1e-3
    assert comparison["initial_cross_norm"] > 0.0
    assert comparison["stores_kernel_propagator"] is False
    hamiltonian = prescribed["hamiltonian"]
    basis = prescribed["local_basis"]
    times = np.linspace(0.0, 1.0, 9)
    dense, streamed = _pair(hamiltonian, basis, times, np.eye(6, dtype=np.complex128), substeps=1)
    later, earlier = 8, 0
    assembled = causal_memory_kernel(hamiltonian, dense["frame"], 2, times, later, earlier, substeps=1)
    direct = (
        dense["coupling_block"][later]
        @ dense["exterior_propagator"][later]
        @ dense["exterior_propagator"][earlier].conj().T
        @ dense["coupling_block"][earlier].conj().T
    )
    assert np.allclose(assembled["kernel"], direct, atol=1e-12, rtol=0)
    frame = dense["frame"]

    def exterior_at(time):
        return hamiltonian_blocks(hamiltonian(time), frame, 2)[2]

    frozen = np.array(dense["coupling_block"][earlier].conj().T, copy=True)
    for left in range(earlier, later):
        width = float(times[left + 1] - times[left])
        frozen = expm(-1j * width * exterior_at(times[left])) @ frozen
    surrogate = dense["coupling_block"][later] @ frozen
    surrogate_gap = float(np.linalg.norm(surrogate - direct, ord="fro"))
    assert surrogate_gap > 1e3 * comparison["kernel_gap"]
    assert surrogate_gap > 1e-4
    zero = causal_memory_kernel(
        hamiltonian,
        frame,
        2,
        times,
        later,
        earlier,
        substeps=1,
        coupling=False,
    )
    assert np.linalg.norm(zero["kernel"]) == 0.0
    transformed = frame.conj().T @ prescribed["covariance"] @ frame
    reversed_kernel = exterior_kernels_from_cohorts(
        float(times[earlier]),
        float(times[later]),
        assembled["cohort_earlier"],
        assembled["cohort_later"],
        transformed[2:, 2:],
        transformed[2:, :2],
    )
    cohort_kernels = exterior_kernels_from_cohorts(
        float(times[later]),
        float(times[later]),
        assembled["cohort_later"],
        assembled["cohort_later"],
        transformed[2:, 2:],
        transformed[2:, :2],
    )
    source = dense["coupling_block"][later] @ dense["exterior_propagator"][later]
    assert np.allclose(
        cohort_kernels["retarded"],
        -0.5j * (source @ source.conj().T),
        atol=1e-12,
        rtol=0,
    )
    assert np.linalg.norm(reversed_kernel["retarded"]) < 1e-15


def test_history_binds_the_conjugate_source_cohort(prescribed):
    hamiltonian = prescribed["hamiltonian"]
    basis = prescribed["local_basis"]
    times = np.linspace(0.0, 0.4, 9)
    streamed = evolve_retained_region(
        hamiltonian,
        basis,
        times,
        np.eye(6, dtype=np.complex128),
        backend="streamed",
        propagator_substeps=3,
    )
    frame = streamed["frame"]

    def generator(time):
        return hamiltonian_blocks(hamiltonian(time), frame, 2)[2]

    step = float(times[1] - times[0])
    accumulator = np.zeros_like(streamed["history"][0])
    wrong = np.zeros_like(accumulator)
    for index in range(times.size - 1):
        source_now = streamed["coupling_block"][index].conj().T @ streamed["amplitudes"][index]
        source_new = streamed["coupling_block"][index + 1].conj().T @ streamed["amplitudes"][index + 1]
        plain_now = streamed["coupling_block"][index].T @ streamed["amplitudes"][index]
        plain_new = streamed["coupling_block"][index + 1].T @ streamed["amplitudes"][index + 1]
        left = exterior_column_transport(
            generator, times, index, source_now, substeps=3, adjoint=True
        )
        right = exterior_column_transport(
            generator, times, index + 1, source_new, substeps=3, adjoint=True
        )
        wrong_left = exterior_column_transport(
            generator, times, index, plain_now, substeps=3, adjoint=True
        )
        wrong_right = exterior_column_transport(
            generator, times, index + 1, plain_new, substeps=3, adjoint=True
        )
        accumulator = accumulator + (step / 2.0) * (left + right)
        wrong = wrong + (step / 2.0) * (wrong_left + wrong_right)
        assert np.allclose(accumulator, streamed["history"][index + 1], atol=1e-12, rtol=0)
    assert float(np.linalg.norm(accumulator - wrong)) > 1e-6


def test_linear_operator_action_matches_the_dense_block(prescribed, monkeypatch):
    hamiltonian = prescribed["hamiltonian"]
    basis = prescribed["local_basis"]
    times = np.linspace(0.0, 0.5, 7)
    prepared = evolve_retained_region(
        hamiltonian,
        basis,
        times,
        np.eye(6, dtype=np.complex128),
        covariance=prescribed["covariance"],
        backend="streamed",
        propagator_substeps=2,
    )
    frame = prepared["frame"]

    def matrix_generator(time):
        return hamiltonian_blocks(hamiltonian(time), frame, 2)[2]

    def operator_generator(time):
        matrix = np.asarray(matrix_generator(time), dtype=np.complex128)

        def matvec(vector, matrix=matrix):
            return matrix @ np.asarray(vector, dtype=np.complex128)

        return LinearOperator(matrix.shape, matvec=matvec, rmatvec=matvec, dtype=np.complex128)

    generator_rng = np.random.default_rng(5)
    columns = generator_rng.normal(size=(4, 3)) + 1j * generator_rng.normal(size=(4, 3))
    dense_action = exterior_column_transport(
        matrix_generator, times, times.size - 1, columns, substeps=3, adjoint=False
    )
    operator_action = exterior_column_transport(
        operator_generator, times, times.size - 1, columns, substeps=3, adjoint=False
    )
    adjoint_action = exterior_column_transport(
        operator_generator, times, times.size - 1, columns, substeps=3, adjoint=True
    )
    reference = np.array(columns, dtype=np.complex128, copy=True)
    for left in range(times.size - 1):
        width = float(times[left + 1] - times[left])
        sample = width / 3
        for part in range(3):
            midpoint = float(times[left] + (part + 0.5) * sample)
            reference = expm(-1j * sample * matrix_generator(midpoint)) @ reference
    assert np.allclose(operator_action, dense_action, atol=1e-12, rtol=0)
    assert np.allclose(dense_action, reference, atol=1e-12, rtol=0)
    identity = np.eye(4, dtype=np.complex128)
    propagated = exterior_column_transport(
        matrix_generator, times, times.size - 1, identity, substeps=3, adjoint=False
    )
    assert np.allclose(adjoint_action, propagated.conj().T @ columns, atol=1e-12, rtol=0)

    original = reduction._expm_multiply_columns

    def as_operator(operator, coefficient, state, workspace):
        if isinstance(operator, LinearOperator):
            return original(operator, coefficient, state, workspace)
        matrix = np.asarray(operator, dtype=np.complex128)
        wrapped = LinearOperator(
            matrix.shape,
            matvec=lambda vector, matrix=matrix: matrix @ np.asarray(vector, dtype=np.complex128),
            rmatvec=lambda vector, matrix=matrix: matrix.conj().T @ np.asarray(vector, dtype=np.complex128),
            dtype=np.complex128,
        )
        return original(wrapped, coefficient, state, workspace)

    monkeypatch.setattr(reduction, "_expm_multiply_columns", as_operator)
    wrapped_run = evolve_retained_region(
        hamiltonian,
        basis,
        times,
        np.eye(6, dtype=np.complex128),
        covariance=prescribed["covariance"],
        backend="streamed",
        propagator_substeps=2,
    )
    assert np.allclose(wrapped_run["amplitudes"], prepared["amplitudes"], atol=1e-11, rtol=0)
    assert np.allclose(wrapped_run["covariance_total"], prepared["covariance_total"], atol=1e-11, rtol=0)
    interval = midpoint_exponential_action(
        operator_generator,
        0.0,
        0.2,
        columns[:, :1],
        substeps=4,
        adjoint=True,
    )
    assert interval.shape == (4, 1)


def test_streamed_solve_is_the_retained_block_only(prescribed, monkeypatch):
    shapes = []
    original = reduction.np.linalg.solve

    def spy(matrix, right_hand):
        shapes.append(matrix.shape)
        return original(matrix, right_hand)

    def reject_dense_propagator(*_args, **_kwargs):
        raise AssertionError("streamed evolution requested the dense propagator")

    def reject_dense_history(*_args, **_kwargs):
        raise AssertionError("streamed evolution used the dense history step")

    monkeypatch.setattr(reduction.np.linalg, "solve", spy)
    monkeypatch.setattr(reduction, "free_exterior_propagator", reject_dense_propagator)
    monkeypatch.setattr(reduction, "_volterra", reject_dense_history)
    times = np.linspace(0.0, 0.4, 6)
    evolved = evolve_retained_region(
        prescribed["hamiltonian"],
        prescribed["local_basis"],
        times,
        prescribed["local_basis"][:, :1],
        backend="streamed",
        propagator_substeps=2,
    )
    assert shapes
    assert all(shape == (2, 2) for shape in shapes)
    assert len(shapes) == times.size - 1
    assert evolved["exterior_propagator"] is None
    assert evolved["allocation"]["dense_exterior_propagator_frames"] == 0


def test_input_validation_covers_the_streamed_backend(prescribed):
    basis = prescribed["local_basis"]
    times = np.linspace(0.0, 0.4, 5)
    identity = np.eye(6, dtype=np.complex128)
    with pytest.raises(ValueError, match="Berry"):
        evolve_retained_region(
            prescribed["hamiltonian"],
            lambda time: basis,
            times,
            identity,
            backend="streamed",
        )
    with pytest.raises(ValueError, match="uniform"):
        evolve_retained_region(
            prescribed["hamiltonian"],
            basis,
            np.array([0.0, 0.1, 0.4]),
            identity,
            backend="streamed",
        )

    def nonhermitian(_time):
        matrix = np.zeros((6, 6), dtype=np.complex128)
        matrix[0, 1] = 1.0
        return matrix

    with pytest.raises(ValueError, match="Hermitian"):
        evolve_retained_region(nonhermitian, basis, times, identity, backend="streamed")
    with pytest.raises(ValueError, match="0<=C<=I"):
        evolve_retained_region(
            prescribed["hamiltonian"],
            basis,
            times,
            covariance=2.0 * identity,
            backend="streamed",
        )
    with pytest.raises(ValueError, match="weights"):
        evolve_retained_region(
            prescribed["hamiltonian"],
            basis,
            times,
            weights=np.ones(2),
            backend="streamed",
        )
    with pytest.raises(ValueError, match="cross"):
        evolve_retained_region(
            prescribed["hamiltonian"],
            basis,
            times,
            identity,
            drop_cross_covariance=True,
            backend="streamed",
        )
    with pytest.raises(ValueError, match="outside"):
        causal_memory_kernel(
            prescribed["hamiltonian"],
            prescribed["frame_probe"],
            2,
            times,
            5,
            0,
        )
    with pytest.raises(ValueError, match="index"):
        causal_memory_kernel(
            prescribed["hamiltonian"],
            prescribed["frame_probe"],
            2,
            times,
            True,
            0,
        )
    with pytest.raises(ValueError, match="ordered"):
        midpoint_exponential_action(lambda time: np.eye(2), 0.4, 0.1, np.eye(2))
    with pytest.raises(TypeError, match="callable"):
        exterior_column_transport(np.eye(2), times, 1, np.eye(2))


def _allocation_problem(dimension, retained, columns):
    generator = np.random.default_rng(19 + dimension)
    def herm():
        raw = generator.normal(size=(dimension, dimension)) + 1j * generator.normal(
            size=(dimension, dimension)
        )
        return 0.5 * (raw + raw.conj().T)

    base, first, second = herm(), herm(), herm()
    draw = generator.normal(size=(dimension, retained)) + 1j * generator.normal(
        size=(dimension, retained)
    )
    orthogonal, _triangle = np.linalg.qr(draw)
    basis = orthogonal * np.exp(1j * generator.normal(size=retained))
    initial = generator.normal(size=(dimension, columns)) + 1j * generator.normal(
        size=(dimension, columns)
    )

    def hamiltonian(time):
        pulse = float(time)
        return base + np.sin(2.0 * pulse) * first + np.cos(1.1 * pulse) * second

    return hamiltonian, basis, initial


def _watch_shapes(monkeypatch, real_empty, real_zeros):
    shapes = []

    def empty(shape, *args, **kwargs):
        shapes.append(tuple(shape) if isinstance(shape, tuple) else (int(shape),))
        return real_empty(shape, *args, **kwargs)

    def zeros(shape, *args, **kwargs):
        shapes.append(tuple(shape) if isinstance(shape, tuple) else (int(shape),))
        return real_zeros(shape, *args, **kwargs)

    monkeypatch.setattr(reduction.np, "empty", empty)
    monkeypatch.setattr(reduction.np, "zeros", zeros)
    return shapes


def _square_exterior_trajectories(shapes, times, exterior_dimension):
    forbidden = []
    for shape in shapes:
        if len(shape) >= 3 and shape[-1] == exterior_dimension and shape[-2] == exterior_dimension:
            if shape[0] >= times:
                forbidden.append(shape)
    return forbidden


def test_allocation_has_no_time_indexed_exterior_propagator(monkeypatch):
    dimension = 12
    retained = 2
    columns = 2
    exterior = dimension - retained
    hamiltonian, basis, initial = _allocation_problem(dimension, retained, columns)
    probe_times = np.linspace(0.0, 0.2, 4)
    probe = evolve_retained_region(hamiltonian, basis, probe_times, initial, backend="dense")
    first_exterior = hamiltonian_blocks(hamiltonian(0.0), probe["frame"], retained)[2]
    second_exterior = hamiltonian_blocks(hamiltonian(0.3), probe["frame"], retained)[2]
    commutator = first_exterior @ second_exterior - second_exterior @ first_exterior
    assert float(np.linalg.norm(first_exterior - first_exterior.conj().T)) < 1e-10
    assert float(np.linalg.norm(commutator)) > 1e-3

    counts = {}
    histories = {}
    frames = {}
    transported = {}
    real_empty = reduction.np.empty
    real_zeros = reduction.np.zeros
    for steps in (8, 16):
        times = np.linspace(0.0, 0.5, steps + 1)
        shapes = _watch_shapes(monkeypatch, real_empty, real_zeros)
        streamed = evolve_retained_region(
            hamiltonian,
            basis,
            times,
            initial,
            backend="streamed",
            propagator_substeps=1,
        )
        assert _square_exterior_trajectories(shapes, times.size, exterior) == []
        allocation = streamed["allocation"]
        assert allocation["time_indexed_exterior_propagator_bytes"] == 0
        assert allocation["dense_exterior_propagator_frames"] == 0
        assert allocation["propagator_storage"] == "none"
        assert allocation["frame_completion"] == "dense-null-space-once"
        assert streamed["exterior_propagator"] is None
        assert streamed["history"].shape == (times.size, exterior, columns)
        assert allocation["history_bytes"] == times.size * exterior * columns * 16
        assert allocation["current_exterior_block_bytes"] == exterior * exterior * 16
        assert allocation["coupling_bytes"] == times.size * retained * exterior * 16
        assert allocation["frame_bytes"] == dimension * dimension * 16
        assert allocation["max_transported_columns"] == columns + retained + columns
        assert allocation["max_transported_columns"] < exterior
        assert allocation["source_cohort_width"] == retained
        counts[steps] = allocation["time_indexed_exterior_propagator_bytes"]
        histories[steps] = allocation["history_bytes"]
        frames[steps] = allocation["frame_bytes"]
        transported[steps] = allocation["max_transported_columns"]
        dense_shapes = _watch_shapes(monkeypatch, real_empty, real_zeros)
        dense = evolve_retained_region(
            hamiltonian,
            basis,
            times,
            initial,
            backend="dense",
            propagator_substeps=1,
        )
        assert (times.size, exterior, exterior) in dense_shapes
        assert dense["allocation"]["time_indexed_exterior_propagator_bytes"] == (
            times.size * exterior * exterior * 16
        )
        assert np.allclose(streamed["amplitudes"], dense["amplitudes"], atol=1e-12, rtol=0)
    assert counts[8] == 0 and counts[16] == 0
    assert histories[16] / histories[8] == pytest.approx(17 / 9)
    assert frames[8] == frames[16]
    assert transported[8] == transported[16]
    dense_short = (8 + 1) * exterior * exterior * 16
    dense_long = (16 + 1) * exterior * exterior * 16
    assert dense_long > dense_short
    assert histories[16] < dense_long
    assert histories[16] / dense_long == pytest.approx(columns / exterior)
