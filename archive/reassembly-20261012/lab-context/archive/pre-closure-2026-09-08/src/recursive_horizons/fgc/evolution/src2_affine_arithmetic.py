"""Outcome-neutral arithmetic attacks on the evolved GR-0 affine source wall.

SRC2 does not change the REF1 equations, the PROTO11 derivative map, or the
strict raw residual gate.  It exposes the binary64 affine system used by the
existing source evaluator and compares algebraically equivalent solution
routes.  Every candidate is judged by a fresh call to the complete
unredefined residual authority, never by the reconstructed affine system
alone.

The exact-rational operations in this module operate only on captured
binary64 coefficients.  They diagnose finite-precision linear algebra; they
do not claim exact continuum coefficients or introduce a new physical model.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isfinite
from numbers import Real
from typing import Iterable, Sequence

import numpy as np

from .cal2_source_diagnosis import diagnose_gr0_grid_accelerations
from .gr0_calibration import FIELD_COUNT
from .gr0_direct_source import gr0_ref1_residual_batch


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _array(name: str, value: object, *, ndim: int) -> np.ndarray:
    answer = np.asarray(value, dtype=np.float64)
    if answer.ndim != ndim or not np.all(np.isfinite(answer)):
        raise ValueError(f"{name} must be a finite rank-{ndim} binary64 array")
    return np.ascontiguousarray(answer)


def _immutable(name: str, value: object, *, ndim: int) -> np.ndarray:
    answer = _array(name, value, ndim=ndim).copy()
    answer.setflags(write=False)
    return answer


def _validated_jet(
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radii: object,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    base = _array("u", u, ndim=2)
    velocity = _array("p", p, ndim=2)
    gradient = _array("q", q, ndim=2)
    mixed = _array("p_r", p_r, ndim=2)
    radial_second = _array("q_r", q_r, ndim=2)
    radius = _array("radii", radii, ndim=1)
    points, fields = base.shape
    if fields != FIELD_COUNT:
        raise ValueError("u must contain the six canonical ADM fields")
    if any(
        item.shape != base.shape
        for item in (velocity, gradient, mixed, radial_second)
    ) or radius.shape != (points,):
        raise ValueError("SRC2 source-jet shapes differ")
    if np.any(radius <= 0.0):
        raise ValueError("SRC2 radii must be strictly positive")
    return base, velocity, gradient, mixed, radial_second, radius


@dataclass(frozen=True, slots=True)
class CapturedAffineSystem:
    """One binary64 representation of the pointwise six-field affine source."""

    label: str
    seed_scale: float
    constant: np.ndarray
    jacobian: np.ndarray

    def __post_init__(self) -> None:
        if not isinstance(self.label, str) or not self.label:
            raise ValueError("affine-system label must be nonempty")
        scale = _finite("seed_scale", self.seed_scale, positive=True)
        constant = _immutable("constant", self.constant, ndim=2)
        jacobian = _immutable("jacobian", self.jacobian, ndim=3)
        if constant.shape[1] != FIELD_COUNT or jacobian.shape != (
            constant.shape[0],
            FIELD_COUNT,
            FIELD_COUNT,
        ):
            raise ValueError("captured affine-system shapes differ")
        object.__setattr__(self, "seed_scale", scale)
        object.__setattr__(self, "constant", constant)
        object.__setattr__(self, "jacobian", jacobian)

    @property
    def point_count(self) -> int:
        return int(self.constant.shape[0])


@dataclass(frozen=True, slots=True)
class ArithmeticCandidate:
    """Complete-residual assessment of one algebraically equivalent route."""

    label: str
    affine_system_label: str
    acceleration: np.ndarray
    complete_residual: np.ndarray
    complete_residual_infinity: float
    reconstructed_residual_infinity: float
    strict_raw_gate_passed: bool
    maximum_residual_point_index: int
    maximum_residual_row_index: int
    maximum_residual_radius: float
    complete_refinement_history: tuple[float, ...]
    exact_binary_refinement_steps: int
    selected_exact_binary_points: tuple[int, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.label, str) or not self.label:
            raise ValueError("candidate label must be nonempty")
        if not isinstance(self.affine_system_label, str) or not self.affine_system_label:
            raise ValueError("candidate affine-system label must be nonempty")
        acceleration = _immutable("acceleration", self.acceleration, ndim=2)
        residual = _immutable("complete_residual", self.complete_residual, ndim=2)
        if acceleration.shape != residual.shape or acceleration.shape[1] != FIELD_COUNT:
            raise ValueError("candidate acceleration and residual shapes differ")
        object.__setattr__(self, "acceleration", acceleration)
        object.__setattr__(self, "complete_residual", residual)
        for name in (
            "complete_residual_infinity",
            "reconstructed_residual_infinity",
            "maximum_residual_radius",
        ):
            value = _finite(name, getattr(self, name))
            if value < 0.0:
                raise ValueError(f"{name} must be nonnegative")
            object.__setattr__(self, name, value)
        if not isinstance(self.strict_raw_gate_passed, bool):
            raise TypeError("strict_raw_gate_passed must be bool")
        if (
            isinstance(self.maximum_residual_point_index, bool)
            or not isinstance(self.maximum_residual_point_index, int)
            or self.maximum_residual_point_index < 0
            or self.maximum_residual_point_index >= acceleration.shape[0]
        ):
            raise ValueError("maximum_residual_point_index lies outside the fixture")
        if (
            isinstance(self.maximum_residual_row_index, bool)
            or not isinstance(self.maximum_residual_row_index, int)
            or not 0 <= self.maximum_residual_row_index < FIELD_COUNT
        ):
            raise ValueError("maximum_residual_row_index lies outside the source rows")
        if not self.complete_refinement_history:
            raise ValueError("candidate refinement history must be nonempty")
        if any(
            not isfinite(value) or value < 0.0
            for value in self.complete_refinement_history
        ):
            raise ValueError("candidate refinement history must be finite and nonnegative")
        if (
            isinstance(self.exact_binary_refinement_steps, bool)
            or not isinstance(self.exact_binary_refinement_steps, int)
            or self.exact_binary_refinement_steps < 0
        ):
            raise ValueError("exact_binary_refinement_steps must be nonnegative")
        if tuple(sorted(set(self.selected_exact_binary_points))) != (
            self.selected_exact_binary_points
        ):
            raise ValueError("selected exact-binary points must be sorted and unique")


@dataclass(frozen=True, slots=True)
class AffineArithmeticComparison:
    raw_tolerance: float
    baseline: ArithmeticCandidate
    candidates: tuple[ArithmeticCandidate, ...]

    def __post_init__(self) -> None:
        tolerance = _finite("raw_tolerance", self.raw_tolerance, positive=True)
        if not isinstance(self.baseline, ArithmeticCandidate):
            raise TypeError("baseline must be ArithmeticCandidate")
        if not self.candidates or any(
            not isinstance(item, ArithmeticCandidate) for item in self.candidates
        ):
            raise ValueError("at least one arithmetic candidate is required")
        labels = (self.baseline.label,) + tuple(item.label for item in self.candidates)
        if len(set(labels)) != len(labels):
            raise ValueError("arithmetic-candidate labels must be unique")
        object.__setattr__(self, "raw_tolerance", tolerance)

    @property
    def best(self) -> ArithmeticCandidate:
        return min(
            (self.baseline, *self.candidates),
            key=lambda item: (item.complete_residual_infinity, item.label),
        )


def assemble_forward_affine_system(
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radii: object,
) -> CapturedAffineSystem:
    """Reproduce the exact zero-plus-unit affine extraction used by PROTO11."""

    base, velocity, gradient, mixed, radial_second, radius = _validated_jet(
        u, p, q, p_r, q_r, radii
    )
    points = base.shape[0]
    seeds = np.zeros((FIELD_COUNT + 1, points, FIELD_COUNT), dtype=np.float64)
    for field in range(FIELD_COUNT):
        seeds[field + 1, :, field] = 1.0
    evaluated = gr0_ref1_residual_batch(
        base, velocity, gradient, seeds, mixed, radial_second, radius
    ).full_residual
    constant = evaluated[0]
    jacobian = np.moveaxis(evaluated[1:] - constant[None, :, :], 0, 2)
    return CapturedAffineSystem(
        label="forward_zero_plus_unit",
        seed_scale=1.0,
        constant=constant,
        jacobian=jacobian,
    )


def assemble_symmetric_affine_system(
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radii: object,
    *,
    seed_scale: Real,
) -> CapturedAffineSystem:
    """Extract the same affine rows from symmetric power-of-two probes."""

    scale = _finite("seed_scale", seed_scale, positive=True)
    mantissa, _exponent = np.frexp(scale)
    if mantissa != 0.5:
        raise ValueError("symmetric affine seed_scale must be a power of two")
    base, velocity, gradient, mixed, radial_second, radius = _validated_jet(
        u, p, q, p_r, q_r, radii
    )
    points = base.shape[0]
    seeds = np.zeros((1 + 2 * FIELD_COUNT, points, FIELD_COUNT), dtype=np.float64)
    for field in range(FIELD_COUNT):
        seeds[1 + 2 * field, :, field] = scale
        seeds[2 + 2 * field, :, field] = -scale
    evaluated = gr0_ref1_residual_batch(
        base, velocity, gradient, seeds, mixed, radial_second, radius
    ).full_residual
    jacobian = np.empty((points, FIELD_COUNT, FIELD_COUNT), dtype=np.float64)
    denominator = 2.0 * scale
    for field in range(FIELD_COUNT):
        jacobian[:, :, field] = (
            evaluated[1 + 2 * field] - evaluated[2 + 2 * field]
        ) / denominator
    return CapturedAffineSystem(
        label=f"symmetric_power_two_{scale.hex()}",
        seed_scale=scale,
        constant=evaluated[0],
        jacobian=jacobian,
    )


def complete_residual(
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radii: object,
    acceleration: object,
) -> np.ndarray:
    """Evaluate the complete unredefined REF1 residual for one acceleration."""

    base, velocity, gradient, mixed, radial_second, radius = _validated_jet(
        u, p, q, p_r, q_r, radii
    )
    candidate = _array("acceleration", acceleration, ndim=2)
    if candidate.shape != base.shape:
        raise ValueError("SRC2 acceleration shape differs from the captured jet")
    return np.ascontiguousarray(
        gr0_ref1_residual_batch(
            base,
            velocity,
            gradient,
            candidate[None, :, :],
            mixed,
            radial_second,
            radius,
        ).full_residual[0]
    )


def _power_two_inverse_scales(maxima: np.ndarray) -> np.ndarray:
    values = _array("maxima", maxima, ndim=2)
    if values.shape[1] != FIELD_COUNT or np.any(values < 0.0):
        raise ValueError("equilibration maxima have invalid shape or sign")
    scales = np.ones_like(values)
    nonzero = values > 0.0
    _mantissa, exponents = np.frexp(values[nonzero])
    scales[nonzero] = np.ldexp(np.ones(np.count_nonzero(nonzero)), -exponents)
    if not np.all(np.isfinite(scales)) or np.any(scales <= 0.0):
        raise ValueError("power-of-two equilibration scale became invalid")
    return scales


def solve_power_two_equilibrated(
    system: CapturedAffineSystem,
    right_hand_side: object | None = None,
) -> np.ndarray:
    """Solve with exact-in-binary row/column power-of-two equilibration."""

    if not isinstance(system, CapturedAffineSystem):
        raise TypeError("system must be CapturedAffineSystem")
    rhs = (
        -system.constant
        if right_hand_side is None
        else _array("right_hand_side", right_hand_side, ndim=2)
    )
    if rhs.shape != system.constant.shape:
        raise ValueError("equilibrated right-hand side shape differs")
    row_max = np.max(np.abs(system.jacobian), axis=2)
    if np.any(row_max == 0.0):
        raise ValueError("captured affine system has a zero row")
    row_scale = _power_two_inverse_scales(row_max)
    row_scaled = row_scale[:, :, None] * system.jacobian
    scaled_rhs = row_scale * rhs
    column_max = np.max(np.abs(row_scaled), axis=1)
    if np.any(column_max == 0.0):
        raise ValueError("captured affine system has a zero column")
    column_scale = _power_two_inverse_scales(column_max)
    equilibrated = row_scaled * column_scale[:, None, :]
    try:
        scaled_solution = np.linalg.solve(
            equilibrated, scaled_rhs[..., None]
        )[..., 0]
    except np.linalg.LinAlgError as error:
        raise ValueError("equilibrated affine source is singular") from error
    solution = column_scale * scaled_solution
    if not np.all(np.isfinite(solution)):
        raise ValueError("equilibrated affine source became nonfinite")
    return np.ascontiguousarray(solution)


def _fraction_solve(matrix: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """Solve one 6x6 binary64 system exactly as rational coefficients."""

    coefficients = _array("matrix", matrix, ndim=2)
    vector = _array("rhs", rhs, ndim=1)
    if coefficients.shape != (FIELD_COUNT, FIELD_COUNT) or vector.shape != (
        FIELD_COUNT,
    ):
        raise ValueError("exact-binary solve requires one six-field system")
    augmented = [
        [Fraction.from_float(float(value)) for value in coefficients[row]]
        + [Fraction.from_float(float(vector[row]))]
        for row in range(FIELD_COUNT)
    ]
    for column in range(FIELD_COUNT):
        pivot = max(
            range(column, FIELD_COUNT),
            key=lambda row: abs(augmented[row][column]),
        )
        if augmented[pivot][column] == 0:
            raise ValueError("exact-binary affine source is singular")
        if pivot != column:
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(FIELD_COUNT):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor == 0:
                continue
            augmented[row] = [
                left - factor * right
                for left, right in zip(
                    augmented[row], augmented[column], strict=True
                )
            ]
    answer = np.asarray(
        [float(augmented[row][-1]) for row in range(FIELD_COUNT)],
        dtype=np.float64,
    )
    if not np.all(np.isfinite(answer)):
        raise ValueError("exact-binary affine root became nonfinite")
    return answer


def exact_binary_solve(
    system: CapturedAffineSystem,
    point_indices: Iterable[int],
    *,
    right_hand_side: object | None = None,
    base_solution: object | None = None,
) -> np.ndarray:
    """Replace selected point roots by exact-rational binary64 solves."""

    if not isinstance(system, CapturedAffineSystem):
        raise TypeError("system must be CapturedAffineSystem")
    rhs = (
        -system.constant
        if right_hand_side is None
        else _array("right_hand_side", right_hand_side, ndim=2)
    )
    if rhs.shape != system.constant.shape:
        raise ValueError("exact-binary right-hand side shape differs")
    if base_solution is None:
        answer = solve_power_two_equilibrated(system, rhs)
    else:
        answer = _array("base_solution", base_solution, ndim=2).copy()
        if answer.shape != rhs.shape:
            raise ValueError("exact-binary base solution shape differs")
    indices = tuple(point_indices)
    if any(
        isinstance(index, bool)
        or not isinstance(index, int)
        or not 0 <= index < system.point_count
        for index in indices
    ):
        raise ValueError("exact-binary point index lies outside the fixture")
    if len(set(indices)) != len(indices):
        raise ValueError("exact-binary point indices must be unique")
    for index in indices:
        answer[index] = _fraction_solve(system.jacobian[index], rhs[index])
    return np.ascontiguousarray(answer)


def _residual_location(
    residual: np.ndarray, radii: np.ndarray
) -> tuple[float, int, int, float]:
    absolute = np.abs(residual)
    flat = int(np.argmax(absolute))
    point, row = np.unravel_index(flat, residual.shape)
    return float(absolute[point, row]), int(point), int(row), float(radii[point])


def _reconstructed_residual_infinity(
    system: CapturedAffineSystem, acceleration: np.ndarray
) -> float:
    # Long-double accumulation diagnoses the captured affine representation;
    # it never substitutes for the complete binary64 REF1 gate.
    matrix = np.asarray(system.jacobian, dtype=np.longdouble)
    vector = np.asarray(acceleration, dtype=np.longdouble)
    constant = np.asarray(system.constant, dtype=np.longdouble)
    reconstructed = constant + np.einsum(
        "nij,nj->ni", matrix, vector, optimize=False
    )
    return float(np.max(np.abs(reconstructed), initial=np.longdouble(0.0)))


def _candidate(
    *,
    label: str,
    system: CapturedAffineSystem,
    initial_acceleration: np.ndarray,
    jet: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    raw_tolerance: float,
    maximum_refinement_iterations: int,
    exact_binary_point_limit: int,
    exact_binary_selection_ratio: float,
) -> ArithmeticCandidate:
    base, velocity, gradient, mixed, radial_second, radius = jet
    acceleration = _array("initial_acceleration", initial_acceleration, ndim=2).copy()
    if acceleration.shape != system.constant.shape:
        raise ValueError("candidate acceleration and affine system differ")
    history: list[float] = []
    residual = complete_residual(
        base,
        velocity,
        gradient,
        mixed,
        radial_second,
        radius,
        acceleration,
    )
    residual_norm, _point, _row, _radius = _residual_location(residual, radius)
    history.append(residual_norm)
    for _iteration in range(maximum_refinement_iterations):
        correction = solve_power_two_equilibrated(system, -residual)
        proposal = acceleration + correction
        proposal_residual = complete_residual(
            base,
            velocity,
            gradient,
            mixed,
            radial_second,
            radius,
            proposal,
        )
        proposal_norm, _point, _row, _radius = _residual_location(
            proposal_residual, radius
        )
        history.append(proposal_norm)
        if proposal_norm >= residual_norm:
            break
        acceleration = proposal
        residual = proposal_residual
        residual_norm = proposal_norm
        if residual_norm < raw_tolerance:
            break

    point_norms = np.max(np.abs(residual), axis=1)
    selected = tuple(
        sorted(
            int(index)
            for index in np.argsort(point_norms)[::-1][:exact_binary_point_limit]
            if point_norms[index] >= exact_binary_selection_ratio * raw_tolerance
        )
    )
    exact_steps = 0
    if selected:
        for _iteration in range(2):
            proposal = exact_binary_solve(
                system,
                selected,
                right_hand_side=-residual,
                base_solution=np.zeros_like(acceleration),
            )
            proposal = acceleration + proposal
            proposal_residual = complete_residual(
                base,
                velocity,
                gradient,
                mixed,
                radial_second,
                radius,
                proposal,
            )
            proposal_norm, _point, _row, _radius = _residual_location(
                proposal_residual, radius
            )
            history.append(proposal_norm)
            if proposal_norm >= residual_norm:
                break
            acceleration = proposal
            residual = proposal_residual
            residual_norm = proposal_norm
            exact_steps += 1
            if residual_norm < raw_tolerance:
                break

    residual_norm, point, row, maximum_radius = _residual_location(residual, radius)
    return ArithmeticCandidate(
        label=label,
        affine_system_label=system.label,
        acceleration=acceleration,
        complete_residual=residual,
        complete_residual_infinity=residual_norm,
        reconstructed_residual_infinity=_reconstructed_residual_infinity(
            system, acceleration
        ),
        strict_raw_gate_passed=residual_norm < raw_tolerance,
        maximum_residual_point_index=point,
        maximum_residual_row_index=row,
        maximum_residual_radius=maximum_radius,
        complete_refinement_history=tuple(history),
        exact_binary_refinement_steps=exact_steps,
        selected_exact_binary_points=selected,
    )


def compare_affine_arithmetic(
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radii: object,
    *,
    raw_tolerance: Real,
    condition_number_maximum: Real,
    symmetric_seed_scales: Sequence[Real] = (1.0, 16.0, 256.0),
    maximum_refinement_iterations: int = 16,
    exact_binary_point_limit: int = 8,
    exact_binary_selection_ratio: Real = 0.5,
) -> tuple[AffineArithmeticComparison, CapturedAffineSystem, tuple[CapturedAffineSystem, ...]]:
    """Compare stable routes under one unchanged complete raw residual gate."""

    jet = _validated_jet(u, p, q, p_r, q_r, radii)
    tolerance = _finite("raw_tolerance", raw_tolerance, positive=True)
    condition_limit = _finite(
        "condition_number_maximum", condition_number_maximum, positive=True
    )
    if (
        isinstance(maximum_refinement_iterations, bool)
        or not isinstance(maximum_refinement_iterations, int)
        or not 0 <= maximum_refinement_iterations <= 16
    ):
        raise ValueError("maximum_refinement_iterations must lie in [0,16]")
    if (
        isinstance(exact_binary_point_limit, bool)
        or not isinstance(exact_binary_point_limit, int)
        or not 1 <= exact_binary_point_limit <= 64
    ):
        raise ValueError("exact_binary_point_limit must lie in [1,64]")
    selection_ratio = _finite(
        "exact_binary_selection_ratio",
        exact_binary_selection_ratio,
        positive=True,
    )
    if selection_ratio > 1.0:
        raise ValueError("exact_binary_selection_ratio must not exceed one")
    scales = tuple(_finite("symmetric_seed_scale", item, positive=True) for item in symmetric_seed_scales)
    if not scales or len(set(scales)) != len(scales):
        raise ValueError("symmetric seed scales must be nonempty and unique")

    forward = assemble_forward_affine_system(*jet)
    diagnosis = diagnose_gr0_grid_accelerations(
        *jet,
        raw_tolerance=tolerance,
        condition_number_maximum=condition_limit,
        maximum_refinement_iterations=maximum_refinement_iterations,
    )
    # The baseline record is the exact PROTO11 authority.  Construct it
    # directly so no SRC2 refinement route is even evaluated on its behalf.
    baseline_residual = diagnosis.residuals
    baseline_norm, point, row, maximum_radius = _residual_location(
        baseline_residual, jet[-1]
    )
    baseline = ArithmeticCandidate(
        label="PROTO11_binary64_baseline",
        affine_system_label=forward.label,
        acceleration=diagnosis.accelerations,
        complete_residual=baseline_residual,
        complete_residual_infinity=baseline_norm,
        reconstructed_residual_infinity=_reconstructed_residual_infinity(
            forward, diagnosis.accelerations
        ),
        strict_raw_gate_passed=baseline_norm < tolerance,
        maximum_residual_point_index=point,
        maximum_residual_row_index=row,
        maximum_residual_radius=maximum_radius,
        complete_refinement_history=tuple(
            item.residual_infinity for item in diagnosis.iterations
        ),
        exact_binary_refinement_steps=0,
        selected_exact_binary_points=(),
    )

    candidates: list[ArithmeticCandidate] = []
    equilibrated = solve_power_two_equilibrated(forward)
    candidates.append(
        _candidate(
            label="forward_power_two_equilibrated",
            system=forward,
            initial_acceleration=equilibrated,
            jet=jet,
            raw_tolerance=tolerance,
            maximum_refinement_iterations=maximum_refinement_iterations,
            exact_binary_point_limit=exact_binary_point_limit,
            exact_binary_selection_ratio=selection_ratio,
        )
    )
    symmetric_systems = tuple(
        assemble_symmetric_affine_system(*jet, seed_scale=scale)
        for scale in scales
    )
    for system in symmetric_systems:
        candidates.append(
            _candidate(
                label=f"{system.label}_equilibrated",
                system=system,
                initial_acceleration=solve_power_two_equilibrated(system),
                jet=jet,
                raw_tolerance=tolerance,
                maximum_refinement_iterations=maximum_refinement_iterations,
                exact_binary_point_limit=exact_binary_point_limit,
                exact_binary_selection_ratio=selection_ratio,
            )
        )
    return (
        AffineArithmeticComparison(
            raw_tolerance=tolerance,
            baseline=baseline,
            candidates=tuple(candidates),
        ),
        forward,
        symmetric_systems,
    )


__all__ = [
    "AffineArithmeticComparison",
    "ArithmeticCandidate",
    "CapturedAffineSystem",
    "assemble_forward_affine_system",
    "assemble_symmetric_affine_system",
    "compare_affine_arithmetic",
    "complete_residual",
    "exact_binary_solve",
    "solve_power_two_equilibrated",
]
