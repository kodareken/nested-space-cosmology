"""Retained-region reduction by a retarded, semiseparable history integral.

The unknown is the retained amplitude. Its exterior is the free propagator
of the exterior block plus one running integral, not a full-system
Schrödinger equation rewritten in retained labels. A later geometry enters
only by supplying the callable ``H(t)``. This module does not build that
trajectory, and the prescribed six-mode pulse below is not regeneration.

``backend="dense"`` stores the free exterior propagator. ``backend="streamed"``
is opt-in: it applies the same midpoint product, through ``expm_multiply``,
only to the source columns of ``B(s)†``, the drive columns, and the running
memory image. It does not store ``W(t)``.

Units are ``hbar = 1``. The local basis is fixed. A moving projector would
need a Berry term and is rejected. Occupied and empty kernels use ``C_EE``
and ``I-C_EE``. Initial cross blocks stay in the covariance assembly.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import sympy as sp
from scipy.linalg import expm, null_space
from scipy.sparse.linalg import LinearOperator, expm_multiply

from .nsc_influence import _covariance
from .nsc_nested_qualities import finite_window

_HERMITIAN_TOLERANCE = 1e-8
_ORTHONORMAL_TOLERANCE = 1e-9
_COMPLETION_TOLERANCE = 1e-8


def _reject_moving_basis(local_basis):
    if callable(local_basis):
        raise ValueError(
            "time-dependent local projectors are unsupported without a Berry term"
        )


def _uniform_times(times):
    grid = np.asarray(times, dtype=float)
    if grid.ndim != 1 or grid.size < 2 or not np.isfinite(grid).all():
        raise ValueError("times must be a finite grid with at least two points")
    steps = np.diff(grid)
    if np.any(steps <= 0):
        raise ValueError("times must strictly increase")
    scale = max(1.0, float(abs(steps[0])))
    if not np.allclose(steps, steps[0], rtol=0, atol=1e-12 * scale):
        raise ValueError(
            "the history trapezoid requires a uniform step; "
            "refine that grid instead of integrating the full system"
        )
    return grid


def _hermitian_matrix(matrix, label):
    values = np.asarray(matrix, dtype=np.complex128)
    if values.ndim != 2 or values.shape[0] != values.shape[1] or not np.isfinite(values).all():
        raise ValueError(f"{label} must be a finite square matrix")
    defect = float(np.linalg.norm(values - values.conj().T, ord="fro"))
    scale = max(1.0, float(np.linalg.norm(values, ord="fro")))
    if defect > _HERMITIAN_TOLERANCE * scale:
        raise ValueError(f"{label} must be Hermitian")
    return 0.5 * (values + values.conj().T)


def _as_complex_matrix(matrix):
    converted = np.empty((matrix.rows, matrix.cols), dtype=np.complex128)
    for row in range(matrix.rows):
        for col in range(matrix.cols):
            converted[row, col] = complex(matrix[row, col])
    return converted


def _positive_substeps(substeps):
    if isinstance(substeps, bool) or not isinstance(substeps, int) or substeps < 1:
        raise ValueError("propagator refinement must be a positive integer")
    return substeps


def fixed_observer_frame(local_basis):
    """Complete ``T = [V, null_space(V†)]`` without rephasing ``V``.

    QR of ``V`` is not used. Its triangular factor would change observer phases.
    """
    _reject_moving_basis(local_basis)
    supplied = np.array(local_basis, dtype=np.complex128, copy=True)
    if supplied.ndim != 2:
        raise ValueError("local basis must be one fixed matrix of orthonormal columns")
    dimension, retained = supplied.shape
    if retained < 1 or retained >= dimension:
        raise ValueError("retained region must be a proper nonempty subspace")
    if not np.isfinite(supplied).all():
        raise ValueError("local basis must be finite")
    gram = supplied.conj().T @ supplied
    if not np.allclose(gram, np.eye(retained), rtol=0, atol=_ORTHONORMAL_TOLERANCE):
        raise ValueError(
            "local basis must already be orthonormal; it is not replaced by a QR factor"
        )
    exterior = null_space(supplied.conj().T)
    if exterior.shape != (dimension, dimension - retained):
        raise ValueError("null_space(V†) did not complete the observer frame")
    exterior_gram = exterior.conj().T @ exterior
    if not np.allclose(exterior_gram, np.eye(dimension - retained), rtol=0, atol=_COMPLETION_TOLERANCE):
        raise ValueError("exterior completion is not orthonormal")
    leakage = float(np.linalg.norm(supplied.conj().T @ exterior, ord="fro"))
    if leakage > _COMPLETION_TOLERANCE:
        raise ValueError("exterior completion is not orthogonal to the observer")
    frame = np.concatenate((supplied, exterior), axis=1)
    if not np.array_equal(frame[:, :retained], supplied):
        raise RuntimeError("observer columns were not preserved exactly")
    if not np.allclose(frame.conj().T @ frame, np.eye(dimension), rtol=0, atol=_COMPLETION_TOLERANCE):
        raise ValueError("completed frame is not unitary")
    return {
        "local_basis": supplied,
        "exterior_basis": exterior,
        "frame": frame,
    }


def hamiltonian_blocks(matrix, frame, retained_dimension):
    """Return ``(H_AA, H_AE, H_EE)`` in the supplied frame."""
    basis = np.asarray(frame, dtype=np.complex128)
    dimension = basis.shape[0]
    if basis.shape != (dimension, dimension):
        raise ValueError("frame must be square")
    if not isinstance(retained_dimension, int) or isinstance(retained_dimension, bool):
        raise ValueError("retained dimension must be an integer")
    if not 0 < retained_dimension < dimension:
        raise ValueError("retained dimension must describe a proper subspace")
    hamiltonian = _hermitian_matrix(matrix, "Hamiltonian")
    if hamiltonian.shape != (dimension, dimension):
        raise ValueError("Hamiltonian and frame dimensions differ")
    transformed = basis.conj().T @ hamiltonian @ basis
    transformed = 0.5 * (transformed + transformed.conj().T)
    retained = transformed[:retained_dimension, :retained_dimension]
    coupling = transformed[:retained_dimension, retained_dimension:]
    exterior = transformed[retained_dimension:, retained_dimension:]
    return retained, coupling, exterior


def free_exterior_propagator(exterior_hamiltonian, times, substeps=1):
    """Advance ``i Wdot = E(t) W`` by midpoint exponentials.

    ``substeps`` refines each output interval. ``W(times[0]) = I``.
    """
    if not callable(exterior_hamiltonian):
        raise TypeError("exterior Hamiltonian must be a callable E(t)")
    grid = _uniform_times(times)
    refinement = _positive_substeps(substeps)
    initial = _hermitian_matrix(exterior_hamiltonian(float(grid[0])), "exterior Hamiltonian")
    exterior_dimension = initial.shape[0]
    propagator = np.empty((grid.size, exterior_dimension, exterior_dimension), dtype=np.complex128)
    current = np.eye(exterior_dimension, dtype=np.complex128)
    propagator[0] = current
    for left in range(grid.size - 1):
        width = float(grid[left + 1] - grid[left])
        sample = width / refinement
        for part in range(refinement):
            midpoint = float(grid[left] + (part + 0.5) * sample)
            generator = _hermitian_matrix(
                exterior_hamiltonian(midpoint),
                "exterior Hamiltonian",
            )
            if generator.shape != (exterior_dimension, exterior_dimension):
                raise ValueError("exterior Hamiltonian changed dimension")
            current = expm(-1j * sample * generator) @ current
        propagator[left + 1] = current
    return propagator


def _backend_name(backend):
    if not isinstance(backend, str) or backend not in ("dense", "streamed"):
        raise ValueError("backend must be 'dense' or 'streamed'")
    return backend


def _grid_index(index, size, label):
    if isinstance(index, bool) or not isinstance(index, int):
        raise ValueError(f"{label} must be an output-node index")
    if index < 0 or index >= size:
        raise ValueError(f"{label} is outside the time grid")
    return index


def _interval_midpoints(t_left, t_right, substeps):
    width = float(t_right - t_left)
    sample = width / substeps
    return [
        (float(t_left + (part + 0.5) * sample), sample)
        for part in range(substeps)
    ]


def _grid_midpoints(times, index, substeps):
    pairs = []
    for left in range(index):
        pairs.extend(_interval_midpoints(float(times[left]), float(times[left + 1]), substeps))
    return pairs


def _linear_operator_trace(operator):
    dimension = operator.shape[0]
    total = 0j
    basis_vector = np.zeros(dimension, dtype=np.complex128)
    for index in range(dimension):
        basis_vector[index] = 1.0
        total += complex(np.vdot(basis_vector, operator.matvec(basis_vector)))
        basis_vector[index] = 0.0
    return total


def _note_generator(workspace, matrix):
    if workspace is None:
        return
    workspace["current_exterior_block_bytes"] = max(
        workspace["current_exterior_block_bytes"],
        int(np.asarray(matrix).nbytes),
    )


def _note_columns(workspace, *arrays):
    if workspace is None:
        return
    total = sum(int(np.asarray(array).nbytes) for array in arrays)
    workspace["peak_column_workspace_bytes"] = max(
        workspace["peak_column_workspace_bytes"],
        total,
    )
    width = max(int(np.asarray(array).shape[1]) for array in arrays)
    workspace["max_transported_columns"] = max(
        workspace["max_transported_columns"],
        width,
    )


def _expm_multiply_columns(operator, coefficient, columns, workspace):
    """Action of ``exp(coefficient * E)`` on columns.

    A dense block is checked Hermitian and discarded after the multiply.
    A ``LinearOperator`` is the same action without a formed block. The
    trace passed to ``expm_multiply`` is the exact trace of the scaled
    operator; it fixes the scaling input and is not a surrogate evolution.
    """
    if isinstance(operator, LinearOperator):
        if operator.shape != (columns.shape[0], columns.shape[0]):
            raise ValueError("exterior Hamiltonian changed dimension")
        scaled = coefficient * operator
        produced = expm_multiply(scaled, columns, traceA=_linear_operator_trace(scaled))
    else:
        matrix = _hermitian_matrix(operator, "exterior Hamiltonian")
        if matrix.shape != (columns.shape[0], columns.shape[0]):
            raise ValueError("exterior Hamiltonian changed dimension")
        _note_generator(workspace, matrix)
        scaled = coefficient * matrix
        produced = expm_multiply(scaled, columns, traceA=np.trace(scaled))
    produced = np.asarray(produced, dtype=np.complex128)
    if produced.shape != columns.shape:
        raise RuntimeError("exterior column action changed shape")
    return np.ascontiguousarray(produced)


def _act_on_midpoints(generator, samples, columns, *, adjoint, workspace):
    ordered = list(samples)
    if adjoint:
        ordered.reverse()
    current = np.array(columns, dtype=np.complex128, copy=True)
    if current.ndim != 2 or current.shape[1] < 1 or not np.isfinite(current).all():
        raise ValueError("exterior columns must be a finite matrix with one or more columns")
    _note_columns(workspace, current)
    for time, width in ordered:
        coefficient = (1j * width) if adjoint else (-1j * width)
        current = _expm_multiply_columns(generator(float(time)), coefficient, current, workspace)
    return current


def midpoint_exponential_action(
    generator,
    t_left,
    t_right,
    state,
    *,
    substeps=1,
    adjoint=False,
    workspace=None,
):
    """Apply one midpoint interval of ``i Wdot = E(t) W`` to columns.

    ``generator(t)`` returns the exterior block as a Hermitian matrix or as
    a ``LinearOperator``. ``adjoint=False`` applies the forward product used
    by ``free_exterior_propagator``. ``adjoint=True`` applies the conjugate
    transpose of that same product. No propagator matrix is stored.
    ``workspace`` is the reducer's allocation record; callers omit it.
    """
    if not callable(generator):
        raise TypeError("exterior generator must be a callable E(t)")
    left = float(t_left)
    right = float(t_right)
    if not np.isfinite(left) or not np.isfinite(right) or right < left:
        raise ValueError("midpoint interval must be finite and ordered")
    refinement = _positive_substeps(substeps)
    if right == left:
        current = np.array(state, dtype=np.complex128, copy=True)
        if current.ndim != 2 or not np.isfinite(current).all():
            raise ValueError("exterior columns must be a finite matrix")
        return current
    return _act_on_midpoints(
        generator,
        _interval_midpoints(left, right, refinement),
        state,
        adjoint=bool(adjoint),
        workspace=workspace,
    )


def exterior_column_transport(
    generator,
    times,
    index,
    state,
    *,
    substeps=1,
    adjoint=False,
    workspace=None,
):
    """Apply ``W(times[index], times[0])`` or its adjoint to columns.

    The product is the same midpoint composition as
    ``free_exterior_propagator``. ``index`` is an output node. The action
    is recomputed from ``generator`` and is not read from a stored ``W``.
    """
    if not callable(generator):
        raise TypeError("exterior generator must be a callable E(t)")
    grid = _uniform_times(times)
    refinement = _positive_substeps(substeps)
    node = _grid_index(index, grid.size, "transport index")
    return _act_on_midpoints(
        generator,
        _grid_midpoints(grid, node, refinement),
        state,
        adjoint=bool(adjoint),
        workspace=workspace,
    )


def _exterior_generator(hamiltonian, frame, retained_dimension):
    dimension = frame.shape[0]

    def generator(time):
        matrix = _hermitian_matrix(hamiltonian(float(time)), "H(t)")
        if matrix.shape != (dimension, dimension):
            raise ValueError("H(t) changed dimension")
        _local, _coupling, exterior = hamiltonian_blocks(matrix, frame, retained_dimension)
        return exterior

    return generator


def _coupling_at(hamiltonian, frame, retained_dimension, time, *, coupling):
    matrix = _hermitian_matrix(hamiltonian(float(time)), "H(t)")
    if matrix.shape[0] != frame.shape[0]:
        raise ValueError("H(t) changed dimension")
    _local, link, _exterior = hamiltonian_blocks(matrix, frame, retained_dimension)
    if not coupling:
        link = np.zeros_like(link)
    return link


def causal_memory_kernel(
    hamiltonian,
    frame,
    retained_dimension,
    times,
    later,
    earlier,
    *,
    substeps=1,
    coupling=True,
):
    """Assemble ``B(t) U_E(t, s) B(s)†`` from source-column cohorts.

    ``later`` and ``earlier`` are indices into ``times``. Each cohort is
    ``B†`` adjoint-propagated back to ``times[0]`` by the midpoint product.
    The kernel is ``cohort_later† @ cohort_earlier``. This is the integrand
    factor in the retained history, without the retarded ``-i θ`` of the
    boundary-state kernel. No exterior propagator is stored.

    The frame must already be the fixed completion of the observer. This
    function does not rebuild or rephase ``V``.
    """
    _require_callable_hamiltonian(hamiltonian)
    grid = _uniform_times(times)
    refinement = _positive_substeps(substeps)
    later_index = _grid_index(later, grid.size, "later")
    earlier_index = _grid_index(earlier, grid.size, "earlier")
    basis = np.asarray(frame, dtype=np.complex128)
    link_later = _coupling_at(
        hamiltonian,
        basis,
        retained_dimension,
        grid[later_index],
        coupling=bool(coupling),
    )
    link_earlier = _coupling_at(
        hamiltonian,
        basis,
        retained_dimension,
        grid[earlier_index],
        coupling=bool(coupling),
    )
    generator = _exterior_generator(hamiltonian, basis, retained_dimension)
    cohort_later = exterior_column_transport(
        generator,
        grid,
        later_index,
        link_later.conj().T,
        substeps=refinement,
        adjoint=True,
    )
    cohort_earlier = exterior_column_transport(
        generator,
        grid,
        earlier_index,
        link_earlier.conj().T,
        substeps=refinement,
        adjoint=True,
    )
    return {
        "kernel": cohort_later.conj().T @ cohort_earlier,
        "cohort_later": cohort_later,
        "cohort_earlier": cohort_earlier,
        "time_later": float(grid[later_index]),
        "time_earlier": float(grid[earlier_index]),
        "stores_exterior_propagator": False,
    }


def _causal_theta(t, s):
    if t > s:
        return 1.0
    if t == s:
        return 0.5
    return 0.0


def _stream_workspace(frame, exterior_dimension, retained_dimension, output_times):
    return {
        "backend": "streamed",
        "dense_exterior_propagator_frames": 0,
        "time_indexed_exterior_propagator_bytes": 0,
        "history_bytes": 0,
        "amplitude_bytes": 0,
        "coupling_bytes": 0,
        "retained_block_bytes": 0,
        "current_exterior_block_bytes": 0,
        "peak_column_workspace_bytes": 0,
        "max_transported_columns": 0,
        "required_exterior_columns": 0,
        "frame_bytes": int(np.asarray(frame).nbytes),
        "frame_completion": "dense-null-space-once",
        "exterior_dimension": int(exterior_dimension),
        "retained_dimension": int(retained_dimension),
        "column_count": 0,
        "output_times": int(output_times),
        "source_cohort_width": int(retained_dimension),
        "propagator_storage": "none",
        "history_replay": "adjoint-midpoint-on-source-cohorts",
    }


def _column_derivative(
    local,
    amplitude,
    coupling,
    memory_columns,
    drive_columns,
    *,
    memory,
    outside_drive,
):
    """Retained derivative with the causal column images already applied.

    ``memory_columns`` is ``S = U Z``, so ``B S`` is the history integral
    against ``B(t) U_E(t, s) B(s)†``. ``drive_columns`` is ``U Y_0``.
    """
    derivative = -1j * (local @ amplitude)
    if memory:
        derivative = derivative - coupling @ memory_columns
    if outside_drive:
        derivative = derivative - 1j * (coupling @ drive_columns)
    return derivative


def _split_transport(moved, pieces):
    advanced = {}
    cursor = 0
    for label, part in pieces:
        width = part.shape[1]
        advanced[label] = moved[:, cursor : cursor + width]
        cursor += width
    return advanced


def _volterra_streamed(
    local_blocks,
    coupling_blocks,
    exterior_at,
    times,
    initial_retained,
    initial_exterior,
    step,
    *,
    substeps,
    memory,
    outside_drive,
    coupling_active,
    workspace,
):
    """Trapezoid step whose exterior actions are column transports.

    The linear unknown at each step is only the retained matrix. Each
    column of ``X`` contributes one drive column and one memory column.
    The source cohort is ``B†``, so its width is the retained dimension.
    One midpoint interval advances those columns. The interaction-picture
    history reapplies the adjoint of that product to the new source cohort
    and does not read a stored ``W``.
    """
    count = local_blocks.shape[0]
    retained_dimension = local_blocks.shape[1]
    exterior_dimension = initial_exterior.shape[0]
    column_count = initial_retained.shape[1]
    amplitudes = np.empty((count, retained_dimension, column_count), dtype=np.complex128)
    history = np.zeros((count, exterior_dimension, column_count), dtype=np.complex128)
    amplitude = np.array(initial_retained, dtype=np.complex128, copy=True)
    accumulator = np.zeros((exterior_dimension, column_count), dtype=np.complex128)
    memory_columns = np.zeros((exterior_dimension, column_count), dtype=np.complex128)
    drive_columns = np.array(initial_exterior, dtype=np.complex128, copy=True)
    amplitudes[0] = amplitude
    residual = 0.0
    cached_cohort = None
    transport = bool(coupling_active) and (bool(memory) or bool(outside_drive))
    if transport:
        workspace["required_exterior_columns"] = max(
            int(workspace["required_exterior_columns"]),
            int(column_count),
        )
    workspace["column_count"] = max(int(workspace["column_count"]), int(column_count))
    for index in range(count - 1):
        local_now = local_blocks[index]
        local_new = local_blocks[index + 1]
        coupling_now = coupling_blocks[index]
        coupling_new = coupling_blocks[index + 1]
        current_derivative = _column_derivative(
            local_now,
            amplitude,
            coupling_now,
            memory_columns,
            drive_columns,
            memory=memory,
            outside_drive=outside_drive,
        )
        step_matrix = np.eye(retained_dimension, dtype=np.complex128) + 1j * (step / 2.0) * local_new
        bracket = current_derivative
        advanced = {}
        if transport and memory:
            cohort_now = coupling_now.conj().T
            pieces = [("memory", memory_columns), ("cohort", cohort_now)]
            if outside_drive:
                pieces.append(("drive", drive_columns))
            bundle = np.concatenate([part for _label, part in pieces], axis=1)
            _note_columns(workspace, bundle)
            moved = midpoint_exponential_action(
                exterior_at,
                times[index],
                times[index + 1],
                bundle,
                substeps=substeps,
                workspace=workspace,
            )
            advanced = _split_transport(moved, pieces)
            step_matrix = step_matrix + (step * step / 4.0) * (
                coupling_new @ coupling_new.conj().T
            )
            remembered_image = advanced["memory"] + (step / 2.0) * (advanced["cohort"] @ amplitude)
            bracket = bracket - coupling_new @ remembered_image
        elif transport and outside_drive:
            _note_columns(workspace, drive_columns)
            advanced["drive"] = midpoint_exponential_action(
                exterior_at,
                times[index],
                times[index + 1],
                drive_columns,
                substeps=substeps,
                workspace=workspace,
            )
        if transport and outside_drive:
            bracket = bracket - 1j * (coupling_new @ advanced["drive"])
        right_hand = amplitude + (step / 2.0) * bracket
        try:
            updated = np.linalg.solve(step_matrix, right_hand)
        except np.linalg.LinAlgError as error:
            raise np.linalg.LinAlgError(
                "retained trapezoid matrix is singular; refine the uniform step"
            ) from error
        if transport and memory:
            if index == 0:
                interaction_now = coupling_now.conj().T @ amplitude
            else:
                interaction_now = cached_cohort @ amplitude
            cohort_new = exterior_column_transport(
                exterior_at,
                times,
                index + 1,
                coupling_new.conj().T,
                substeps=substeps,
                adjoint=True,
                workspace=workspace,
            )
            interaction_new = cohort_new @ updated
            updated_history = accumulator + (step / 2.0) * (interaction_now + interaction_new)
            updated_memory = advanced["memory"] + (step / 2.0) * (
                advanced["cohort"] @ amplitude + coupling_new.conj().T @ updated
            )
            cached_cohort = cohort_new
        else:
            updated_history = accumulator
            updated_memory = memory_columns
        updated_drive = advanced["drive"] if transport and outside_drive else drive_columns
        updated_derivative = _column_derivative(
            local_new,
            updated,
            coupling_new,
            updated_memory,
            updated_drive,
            memory=memory,
            outside_drive=outside_drive,
        )
        discrepancy = updated - amplitude - (step / 2.0) * (current_derivative + updated_derivative)
        residual = max(residual, float(np.linalg.norm(discrepancy, ord="fro")))
        amplitude = updated
        accumulator = updated_history
        memory_columns = updated_memory
        drive_columns = updated_drive
        amplitudes[index + 1] = amplitude
        history[index + 1] = accumulator
    workspace["history_bytes"] += int(history.nbytes)
    workspace["amplitude_bytes"] += int(amplitudes.nbytes)
    return amplitudes, history, residual


def _require_callable_hamiltonian(hamiltonian):
    if not callable(hamiltonian):
        raise TypeError(
            "H(t) must be a callable; a later geometry is passed as H(g(t)) "
            "and is not constructed here"
        )


def _sample_blocks(hamiltonian, frame, retained_dimension, times, *, coupling):
    dimension = frame.shape[0]
    count = times.size
    exterior_dimension = dimension - retained_dimension
    retained = np.empty((count, retained_dimension, retained_dimension), dtype=np.complex128)
    link = np.empty((count, retained_dimension, exterior_dimension), dtype=np.complex128)
    for index, time in enumerate(times):
        matrix = _hermitian_matrix(hamiltonian(float(time)), "H(t)")
        if matrix.shape != (dimension, dimension):
            raise ValueError("H(t) changed dimension")
        local, coupling_block, _exterior = hamiltonian_blocks(matrix, frame, retained_dimension)
        retained[index] = local
        if coupling:
            link[index] = coupling_block
        else:
            link[index] = 0
    def exterior_at(time):
        matrix = _hermitian_matrix(hamiltonian(float(time)), "H(t)")
        if matrix.shape != (dimension, dimension):
            raise ValueError("H(t) changed dimension")
        _local, _coupling, exterior = hamiltonian_blocks(matrix, frame, retained_dimension)
        return exterior
    return retained, link, exterior_at


def _initial_coefficients(frame, initial_columns, dimension):
    columns = np.asarray(initial_columns, dtype=np.complex128)
    if columns.ndim == 1:
        columns = columns.reshape(dimension, 1)
    if columns.ndim != 2 or columns.shape[0] != dimension or not np.isfinite(columns).all():
        raise ValueError("initial columns must be finite vectors in the Hamiltonian basis")
    return columns, frame.conj().T @ columns


def _split_coefficients(coefficients, retained_dimension):
    return coefficients[:retained_dimension], coefficients[retained_dimension:]


def _xdot(local, source, amplitude, history, initial_exterior, *, memory, outside_drive):
    derivative = -1j * (local @ amplitude)
    if memory:
        derivative = derivative - source @ history
    if outside_drive:
        derivative = derivative - 1j * (source @ initial_exterior)
    return derivative


def _volterra(local_blocks, coupling_blocks, propagator, initial_retained, initial_exterior, step, *, memory, outside_drive):
    """Trapezoid step for ``X`` and ``Z`` on one uniform grid.

    ``Z`` stores ``∫ W(s)† B(s)† X(s) ds``. The linear unknown at each step
    is only the retained matrix. ``K = W† B†``.
    """
    count = local_blocks.shape[0]
    retained_dimension = local_blocks.shape[1]
    exterior_dimension = initial_exterior.shape[0]
    column_count = initial_retained.shape[1]
    amplitudes = np.empty((count, retained_dimension, column_count), dtype=np.complex128)
    history = np.empty((count, exterior_dimension, column_count), dtype=np.complex128)
    amplitude = np.array(initial_retained, dtype=np.complex128, copy=True)
    accumulator = np.zeros((exterior_dimension, column_count), dtype=np.complex128)
    amplitudes[0] = amplitude
    history[0] = accumulator
    residual = 0.0
    for index in range(count - 1):
        local_now = local_blocks[index]
        local_new = local_blocks[index + 1]
        coupling_now = coupling_blocks[index]
        coupling_new = coupling_blocks[index + 1]
        propagator_now = propagator[index]
        propagator_new = propagator[index + 1]
        kernel_now = propagator_now.conj().T @ coupling_now.conj().T
        kernel_new = propagator_new.conj().T @ coupling_new.conj().T
        source_now = coupling_now @ propagator_now
        source_new = coupling_new @ propagator_new
        current_derivative = _xdot(
            local_now,
            source_now,
            amplitude,
            accumulator,
            initial_exterior,
            memory=memory,
            outside_drive=outside_drive,
        )
        step_matrix = np.eye(retained_dimension, dtype=np.complex128) + 1j * (step / 2.0) * local_new
        bracket = current_derivative
        if memory:
            step_matrix = step_matrix + (step * step / 4.0) * (source_new @ kernel_new)
            remembered = accumulator + (step / 2.0) * (kernel_now @ amplitude)
            bracket = bracket - source_new @ remembered
        if outside_drive:
            bracket = bracket - 1j * (source_new @ initial_exterior)
        right_hand = amplitude + (step / 2.0) * bracket
        try:
            updated = np.linalg.solve(step_matrix, right_hand)
        except np.linalg.LinAlgError as error:
            raise np.linalg.LinAlgError(
                "retained trapezoid matrix is singular; refine the uniform step"
            ) from error
        if memory:
            updated_history = accumulator + (step / 2.0) * (
                kernel_now @ amplitude + kernel_new @ updated
            )
        else:
            updated_history = accumulator
        updated_derivative = _xdot(
            local_new,
            source_new,
            updated,
            updated_history,
            initial_exterior,
            memory=memory,
            outside_drive=outside_drive,
        )
        discrepancy = updated - amplitude - (step / 2.0) * (current_derivative + updated_derivative)
        residual = max(residual, float(np.linalg.norm(discrepancy, ord="fro")))
        amplitude = updated
        accumulator = updated_history
        amplitudes[index + 1] = amplitude
        history[index + 1] = accumulator
    return amplitudes, history, residual


def _frame_blocks(covariance, frame, retained_dimension):
    transformed = frame.conj().T @ covariance @ frame
    transformed = 0.5 * (transformed + transformed.conj().T)
    parent = transformed[:retained_dimension, :retained_dimension]
    parent_child = transformed[:retained_dimension, retained_dimension:]
    child_parent = transformed[retained_dimension:, :retained_dimension]
    child = transformed[retained_dimension:, retained_dimension:]
    return parent, parent_child, child_parent, child


def covariance_terms(parent_response, child_response, blocks, *, drop_cross=False):
    """Four blocks of ``R C R†`` at one time, matching the boundary-state split.

    ``parent_response`` multiplies the initial retained amplitude.
    ``child_response`` multiplies the initial exterior amplitude.
    """
    parent_block, parent_child, child_parent, child_block = blocks
    parent_term = parent_response @ parent_block @ parent_response.conj().T
    child_term = child_response @ child_block @ child_response.conj().T
    if drop_cross:
        initial_parent_child = np.zeros_like(parent_term)
        initial_child_parent = np.zeros_like(parent_term)
    else:
        initial_parent_child = parent_response @ parent_child @ child_response.conj().T
        initial_child_parent = child_response @ child_parent @ parent_response.conj().T
    total = parent_term + child_term + initial_parent_child + initial_child_parent
    return {
        "parent": parent_term,
        "child": child_term,
        "initial_parent_child": initial_parent_child,
        "initial_child_parent": initial_child_parent,
        "total": total,
    }


def exterior_kernels(
    t,
    s,
    *,
    coupling_t,
    propagator_t,
    coupling_s,
    propagator_s,
    exterior_covariance,
    initial_cross,
):
    """Occupied, empty, and retarded kernels at a pair of times.

    ``F(t) = B(t) W(t)``. ``θ(0) = 1/2``, as in the fixed-frame boundary
    state. The equal-time value does not enter the history integral.
    ``initial_cross`` is ``C_EA``.
    """
    if not np.isfinite(t) or not np.isfinite(s):
        raise ValueError("kernel times must be finite")
    source_t = np.asarray(coupling_t, dtype=np.complex128) @ np.asarray(propagator_t, dtype=np.complex128)
    source_s = np.asarray(coupling_s, dtype=np.complex128) @ np.asarray(propagator_s, dtype=np.complex128)
    exterior = _hermitian_matrix(exterior_covariance, "exterior covariance")
    if source_t.shape[1] != exterior.shape[0] or source_s.shape != source_t.shape:
        raise ValueError("kernel factors must match the exterior covariance")
    cross = np.asarray(initial_cross, dtype=np.complex128)
    if cross.shape != (exterior.shape[0], source_t.shape[0]):
        raise ValueError("initial cross block must be C_EA")
    identity = np.eye(exterior.shape[0], dtype=np.complex128)
    occupied = source_t @ exterior @ source_s.conj().T
    empty = source_t @ (identity - exterior) @ source_s.conj().T
    theta = _causal_theta(t, s)
    retarded = -1j * theta * (source_t @ source_s.conj().T)
    initial = source_t @ cross
    return {
        "occupied": occupied,
        "empty": empty,
        "retarded": retarded,
        "lesser": 1j * occupied,
        "greater": -1j * empty,
        "keldysh": -1j * (empty - occupied),
        "initial_occupied": initial,
        "initial_empty": -initial,
    }


def exterior_kernels_from_cohorts(
    t,
    s,
    cohort_t,
    cohort_s,
    exterior_covariance,
    initial_cross,
):
    """Occupied, empty, and retarded kernels from adjoint source cohorts.

    ``cohort(t) = U_E(t, t_0)† B(t)†``. Then ``B(t) U_E(t, s) B(s)†`` is
    ``cohort(t)† cohort(s)``, and ``F(t) C F(s)†`` is ``cohort(t)† C cohort(s)``.
    ``θ(0) = 1/2``, as in ``exterior_kernels``. ``initial_cross`` is ``C_EA``.
    """
    if not np.isfinite(t) or not np.isfinite(s):
        raise ValueError("kernel times must be finite")
    later = np.asarray(cohort_t, dtype=np.complex128)
    earlier = np.asarray(cohort_s, dtype=np.complex128)
    if later.ndim != 2 or later.shape != earlier.shape:
        raise ValueError("source cohorts must be matching exterior-by-retained matrices")
    exterior = _hermitian_matrix(exterior_covariance, "exterior covariance")
    if later.shape[0] != exterior.shape[0]:
        raise ValueError("kernel factors must match the exterior covariance")
    cross = np.asarray(initial_cross, dtype=np.complex128)
    if cross.shape != (exterior.shape[0], later.shape[1]):
        raise ValueError("initial cross block must be C_EA")
    identity = np.eye(exterior.shape[0], dtype=np.complex128)
    gram = later.conj().T @ earlier
    occupied = later.conj().T @ exterior @ earlier
    empty = later.conj().T @ (identity - exterior) @ earlier
    theta = _causal_theta(t, s)
    retarded = -1j * theta * gram
    initial = later.conj().T @ cross
    return {
        "occupied": occupied,
        "empty": empty,
        "retarded": retarded,
        "lesser": 1j * occupied,
        "greater": -1j * empty,
        "keldysh": -1j * (empty - occupied),
        "initial_occupied": initial,
        "initial_empty": -initial,
    }


def _car_deviation(amplitudes):
    deviation = 0.0
    identity = np.eye(amplitudes.shape[1], dtype=np.complex128)
    for sample in amplitudes:
        discrepancy = sample @ sample.conj().T - identity
        deviation = max(deviation, float(np.linalg.norm(discrepancy, ord="fro")))
    return deviation


def _is_full_orthonormal(columns, dimension):
    return (
        columns.shape == (dimension, dimension)
        and np.allclose(columns.conj().T @ columns, np.eye(dimension), rtol=0, atol=1e-8)
    )


def _column_covariance(columns, weights):
    scale = np.asarray(weights, dtype=float)
    if scale.ndim != 1 or scale.shape[0] != columns.shape[1] or not np.isfinite(scale).all():
        raise ValueError("weights must be one finite value for each initial column")
    if scale.min() < -1e-12 or scale.max() > 1.0 + 1e-12:
        raise ValueError("column weights must lie in [0, 1]")
    covariance = (columns * scale) @ columns.conj().T
    return _covariance(covariance), scale


def evolve_retained_region(
    hamiltonian,
    local_basis,
    times,
    initial_columns=None,
    *,
    weights=None,
    covariance=None,
    memory=True,
    outside_drive=True,
    coupling=True,
    drop_cross_covariance=False,
    propagator_substeps=1,
    backend="dense",
):
    """Evolve retained amplitudes, and optionally the retained covariance.

    ``initial_columns`` are vectors in the same basis as ``H(t)``. A later
    rank-6 source is just those columns; this path does not build the full
    map. ``covariance`` is the one-body matrix ``C``. Identity columns, or
    this covariance path, expose the retained map ``R`` and the same-time
    check ``R R†`` against ``I``.

    ``memory=False`` drops the history accumulator. ``outside_drive=False``
    drops the initial exterior amplitude from the derivative but still
    records it. ``coupling=False`` zeroes ``H_AE`` only.
    ``drop_cross_covariance=True`` drops ``C_AE`` and ``C_EA`` in the
    covariance assembly and does not change the amplitudes.

    ``backend="dense"`` stores ``W(t)``. ``backend="streamed"`` keeps that
    map's retained equation but transports only source, drive, and memory
    columns. A later trajectory can pass ``lambda t: operator(g(t))`` and
    ``backend="streamed"``; this function does not build ``g(t)``.
    """
    _require_callable_hamiltonian(hamiltonian)
    backend_name = _backend_name(backend)
    if weights is not None and initial_columns is None:
        raise ValueError("weights require initial columns")
    if initial_columns is None and covariance is None:
        raise ValueError("pass initial columns or a full covariance")
    if drop_cross_covariance and covariance is None:
        raise ValueError("dropping the initial cross block requires a covariance")
    completed = fixed_observer_frame(local_basis)
    frame = completed["frame"]
    retained_dimension = completed["local_basis"].shape[1]
    dimension = frame.shape[0]
    grid = _uniform_times(times)
    local_blocks, coupling_blocks, exterior_at = _sample_blocks(
        hamiltonian,
        frame,
        retained_dimension,
        grid,
        coupling=bool(coupling),
    )
    refinement = _positive_substeps(propagator_substeps)
    if backend_name == "dense":
        propagator = free_exterior_propagator(
            exterior_at,
            grid,
            substeps=refinement,
        )
        allocation = {
            "backend": "dense",
            "dense_exterior_propagator_frames": int(grid.size),
            "time_indexed_exterior_propagator_bytes": int(propagator.nbytes),
            "frame_bytes": int(frame.nbytes),
            "frame_completion": "dense-null-space-once",
            "exterior_dimension": int(propagator.shape[-1]),
            "output_times": int(grid.size),
            "propagator_storage": "time-indexed-dense-W",
        }
    else:
        propagator = None
        allocation = _stream_workspace(
            frame,
            dimension - retained_dimension,
            retained_dimension,
            grid.size,
        )
        allocation["coupling_bytes"] = int(coupling_blocks.nbytes)
        allocation["retained_block_bytes"] = int(local_blocks.nbytes)
    step = float(grid[1] - grid[0])
    result = {
        "times": grid,
        "local_basis": completed["local_basis"],
        "exterior_basis": completed["exterior_basis"],
        "frame": frame,
        "retained_block": local_blocks,
        "coupling_block": coupling_blocks,
        "exterior_propagator": propagator,
        "amplitudes": None,
        "history": None,
        "initial_retained": None,
        "initial_exterior": None,
        "column_weight_covariance": None,
        "parent_response": None,
        "child_response": None,
        "retained_map": None,
        "car_deviation_max": None,
        "frame_blocks": None,
        "covariance_total": None,
        "covariance_terms_final": None,
        "covariance_assembly_residual_max": None,
        "trapezoid_residual_max": 0.0,
        "flags": {
            "memory": bool(memory),
            "outside_drive": bool(outside_drive),
            "coupling": bool(coupling),
            "drop_cross_covariance": bool(drop_cross_covariance),
            "propagator_substeps": refinement,
            "backend": backend_name,
            "history_integrator": (
                "uniform-trapezoid-causal-columns"
                if backend_name == "streamed"
                else "uniform-trapezoid-semiseparable"
            ),
            "exterior_integrator": (
                "midpoint-expm-multiply" if backend_name == "streamed" else "midpoint-expm"
            ),
            "stores_dense_exterior_propagator": backend_name == "dense",
        },
    }
    if backend_name == "streamed":
        result["flags"]["memory_kernel"] = "B(t) U_E(t,s) B(s)†"

    def integrate(initial_retained, initial_exterior):
        if backend_name == "dense":
            return _volterra(
                local_blocks,
                coupling_blocks,
                propagator,
                initial_retained,
                initial_exterior,
                step,
                memory=bool(memory),
                outside_drive=bool(outside_drive),
            )
        return _volterra_streamed(
            local_blocks,
            coupling_blocks,
            exterior_at,
            grid,
            initial_retained,
            initial_exterior,
            step,
            substeps=refinement,
            memory=bool(memory),
            outside_drive=bool(outside_drive),
            coupling_active=bool(coupling),
            workspace=allocation,
        )

    residuals = []
    if initial_columns is not None:
        columns, coefficients = _initial_coefficients(frame, initial_columns, dimension)
        initial_retained, initial_exterior = _split_coefficients(coefficients, retained_dimension)
        weight_scale = None
        if weights is not None:
            weight_scale = _column_covariance(columns, weights)[1]
        amplitudes, history, residual = integrate(initial_retained, initial_exterior)
        residuals.append(residual)
        result["amplitudes"] = amplitudes
        result["history"] = history
        result["initial_retained"] = initial_retained
        result["initial_exterior"] = initial_exterior
        if weight_scale is not None:
            weighted = np.empty(
                (amplitudes.shape[0], retained_dimension, retained_dimension),
                dtype=np.complex128,
            )
            for index, sample in enumerate(amplitudes):
                weighted[index] = (sample * weight_scale) @ sample.conj().T
            result["column_weight_covariance"] = weighted
        if _is_full_orthonormal(columns, dimension):
            result["retained_map"] = amplitudes
            result["car_deviation_max"] = _car_deviation(amplitudes)
    if covariance is not None:
        validated = _covariance(covariance)
        if validated.shape != (dimension, dimension):
            raise ValueError("covariance must use the Hamiltonian basis")
        blocks = _frame_blocks(validated, frame, retained_dimension)
        coefficients = frame.conj().T @ frame
        initial_retained, initial_exterior = _split_coefficients(coefficients, retained_dimension)
        amplitudes, history, residual = integrate(initial_retained, initial_exterior)
        residuals.append(residual)
        parent = amplitudes[:, :, :retained_dimension]
        child = amplitudes[:, :, retained_dimension:]
        retained_map = np.einsum("tij,jk->tik", amplitudes, frame.conj().T)
        totals = np.empty(
            (grid.size, retained_dimension, retained_dimension),
            dtype=np.complex128,
        )
        assembly = 0.0
        final_terms = None
        for index in range(grid.size):
            terms = covariance_terms(
                parent[index],
                child[index],
                blocks,
                drop_cross=bool(drop_cross_covariance),
            )
            totals[index] = terms["total"]
            if not drop_cross_covariance:
                direct = retained_map[index] @ validated @ retained_map[index].conj().T
                assembly = max(
                    assembly,
                    float(np.linalg.norm(terms["total"] - direct, ord="fro")),
                )
            if index == grid.size - 1:
                final_terms = terms
        result["parent_response"] = parent
        result["child_response"] = child
        result["history_of_frame"] = history
        if result["retained_map"] is None:
            result["retained_map"] = retained_map
        elif not np.allclose(result["retained_map"], retained_map, rtol=0, atol=1e-8):
            raise RuntimeError("identity columns and the frame map are not the same linear solution")
        result["car_deviation_max"] = _car_deviation(retained_map)
        result["frame_blocks"] = {
            "parent": blocks[0],
            "parent_child": blocks[1],
            "child_parent": blocks[2],
            "child": blocks[3],
        }
        result["covariance_total"] = totals
        result["covariance_terms_final"] = final_terms
        result["covariance_assembly_residual_max"] = (
            None if drop_cross_covariance else assembly
        )
    if residuals:
        result["trapezoid_residual_max"] = max(residuals)
    result["allocation"] = allocation
    return result


@lru_cache(maxsize=1)
def _inherited_window_and_state():
    local = sp.Matrix([[1, sp.I / 5], [-sp.I / 5, 2]])
    link = sp.Matrix(
        [
            [sp.Rational(1, 4), sp.I / 7],
            [sp.Rational(1, 9), sp.Rational(1, 6)],
        ]
    )
    window = finite_window(local, link, sp.Rational(3, 2), 0, 3)
    covariance = sp.eye(window.rows) / 2 - window**3 / 2048
    operator = _hermitian_matrix(_as_complex_matrix(window), "inherited window")
    state = _covariance(_as_complex_matrix(covariance))
    commutator = operator @ state - state @ operator
    if float(np.linalg.norm(commutator, ord="fro")) > 1e-8:
        raise RuntimeError("inherited covariance is not stationary on the frozen window")
    return operator, state


def prescribed_six_mode_problem():
    """Fixed six-mode window plus a prescribed Hermitian pulse.

    ``H(t) = J + 0.2 sin(2t) W`` with ``W[2, 2] = 1`` and
    ``W[0, 2] = W[2, 0] = 0.1``. The pulse is a numerical control. It is
    not a regenerated history and not a coupled geometry. ``J`` and ``C``
    are the inherited finite-turnover window and correlated Gaussian state.
    """
    window, state = _inherited_window_and_state()
    operator = np.array(window, dtype=np.complex128, copy=True)
    covariance = np.array(state, dtype=np.complex128, copy=True)
    perturbation = np.zeros((6, 6), dtype=np.complex128)
    perturbation[2, 2] = 1.0
    perturbation[0, 2] = 0.1
    perturbation[2, 0] = 0.1
    local_basis = np.zeros((6, 2), dtype=np.complex128)
    local_basis[0, 0] = 1.0
    local_basis[1, 1] = 1.0

    def hamiltonian(time):
        pulse = 0.2 * np.sin(2.0 * float(time))
        return operator + pulse * perturbation

    return {
        "window": operator.copy(),
        "covariance": covariance.copy(),
        "perturbation": perturbation.copy(),
        "pulse_amplitude": 0.2,
        "pulse_frequency": 2.0,
        "hamiltonian": hamiltonian,
        "local_basis": local_basis.copy(),
        "dimension": 6,
        "retained_dimension": 2,
        "prescribed_history": True,
        "regeneration": False,
        "coupled_trajectory": False,
    }
