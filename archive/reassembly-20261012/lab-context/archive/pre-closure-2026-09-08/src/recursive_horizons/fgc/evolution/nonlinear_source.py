"""Safeguarded floating solve of the complete six-row REF1 residual.

The physical residual is always evaluated by
``modified_harmonic_full_residuals``.  ``FloatJet2`` supplies binary64 scalar
arithmetic to that same tensor code, and ``FloatFirstTangent`` differentiates
the six acceleration columns by forward automatic differentiation.  This file
therefore adds a solver and diagnostics, not a second set of field equations.

The current certified use is the compact QIFT1 local branch box.  It is not a
regular-centre formulation, a continuum evolution system, a collapse result,
or retained-EFT authorization.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isfinite
from numbers import Real
from typing import Any, Mapping, Sequence

import numpy as np

from ..modified_harmonic_reference import modified_harmonic_full_residuals
from ..reference_connection import ReferenceConnection
from ..spherical_reduction import BASE_FIELD_ORDER, Jet2, SphericalState
from .floating_jet import (
    FloatFirstTangent,
    FloatJet2,
    float_jet,
    scalar_primal,
    scalar_primal_and_tangent,
)


FIELD_COUNT = len(BASE_FIELD_ORDER)
PARAMETER_GROUPS = ("u", "p", "q", "p_r", "q_r")
PARAMETER_SLOTS = {
    "u": "value",
    "p": "dt",
    "q": "dr",
    "p_r": "dtr",
    "q_r": "drr",
}
PARAMETER_ORDER = tuple(
    f"{group}.{field}" for group in PARAMETER_GROUPS for field in BASE_FIELD_ORDER
)


STOP_OUTCOME_LABELS = {
    "nonfinite_input": "invalid_implementation_or_nonconverged_run",
    "invalid_state": "invalid_implementation_or_nonconverged_run",
    "parameter_box_violation": "stopped_branch_loss",
    "initial_acceleration_outside_authorized_box": "stopped_branch_loss",
    "kinetic_condition_limit": "stopped_branch_loss",
    "linear_solve_failure": "invalid_implementation_or_nonconverged_run",
    "nonfinite_residual_or_jacobian": "invalid_implementation_or_nonconverged_run",
    "residual_failed_to_decrease": "invalid_implementation_or_nonconverged_run",
    "iteration_limit": "invalid_implementation_or_nonconverged_run",
    "branch_displacement_limit": "stopped_branch_loss",
    "final_acceleration_outside_authorized_box": "stopped_branch_loss",
}


def _finite(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    return answer


def _vector(name: str, value: Sequence[Real], *, length: int = FIELD_COUNT) -> np.ndarray:
    if isinstance(value, (str, bytes)) or len(value) != length:
        raise ValueError(f"{name} must contain {length} entries")
    try:
        array = np.asarray([_finite(f"{name}[{index}]", item) for index, item in enumerate(value)], dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain only finite real entries") from exc
    if array.shape != (length,) or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a finite length-{length} vector")
    return array


@dataclass(frozen=True, slots=True)
class NonlinearSolverConfig:
    residual_tolerance: float = 1.0e-12
    maximum_iterations: int = 16
    condition_number_maximum: float = 1.0e10
    normalized_branch_displacement_maximum: float = 1.0 / 16.0
    parameter_half_width: float = 1.0 / 65536.0
    acceleration_half_width: float = 1.0 / 128.0
    maximum_backtracks: int = 24
    minimum_step_fraction: float = 1.0 / (2**24)

    def __post_init__(self) -> None:
        for name in (
            "residual_tolerance",
            "condition_number_maximum",
            "normalized_branch_displacement_maximum",
            "parameter_half_width",
            "acceleration_half_width",
            "minimum_step_fraction",
        ):
            value = _finite(name, getattr(self, name))
            if value <= 0.0:
                raise ValueError(f"{name} must be positive")
            object.__setattr__(self, name, value)
        for name in ("maximum_iterations", "maximum_backtracks"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ValueError(f"{name} must be a positive integer")


@dataclass(frozen=True, slots=True)
class NewtonIteration:
    iteration: int
    accepted_step_fraction: float
    residual_infinity_before: float
    residual_infinity_after: float
    jacobian_condition_infinity: float
    acceleration_infinity: float
    branch_displacement_infinity: float


@dataclass(frozen=True, slots=True)
class AccelerationSolveResult:
    acceleration_order: tuple[str, ...]
    accelerations: tuple[float, ...]
    residual_vector: tuple[float, ...]
    residual_infinity: float
    iterations: int
    jacobian_condition_infinity: float
    branch_displacement_infinity: float
    iteration_history: tuple[NewtonIteration, ...]
    converged: bool
    branch_continuity_preserved: bool
    physical_residual: str = "complete_unredefined_ACT1_VAR1_REF1_rows"
    differentiation_method: str = "binary64_first_tangent_forward_AD_through_same_REF1_evaluator"


class AccelerationSolveStop(RuntimeError):
    """Typed fail-closed termination before a solver premise is crossed."""

    def __init__(self, reason: str, message: str, diagnostics: Mapping[str, Any] | None = None):
        if reason not in STOP_OUTCOME_LABELS:
            raise ValueError("unknown acceleration stop reason")
        super().__init__(message)
        self.reason = reason
        self.outcome_label = STOP_OUTCOME_LABELS[reason]
        self.diagnostics = dict(diagnostics or {})


def parameter_vector(state: SphericalState) -> tuple[float, ...]:
    if not isinstance(state, SphericalState):
        raise TypeError("state must be a SphericalState")
    return tuple(
        _finite(
            f"{group}.{field}",
            scalar_primal(getattr(getattr(state, field), PARAMETER_SLOTS[group])),
        )
        for group in PARAMETER_GROUPS
        for field in BASE_FIELD_ORDER
    )


def state_from_parameter_vector(
    center: SphericalState,
    parameters: Sequence[Real],
    *,
    accelerations: Sequence[Real] | None = None,
) -> SphericalState:
    """Build one floating state in the frozen QIFT1 parameter order."""

    if not isinstance(center, SphericalState):
        raise TypeError("center must be a SphericalState")
    values = _vector("parameters", parameters, length=len(PARAMETER_ORDER))
    acceleration = (
        np.asarray([_finite(f"center acceleration {field}", getattr(center, field).dtt) for field in BASE_FIELD_ORDER], dtype=np.float64)
        if accelerations is None
        else _vector("accelerations", accelerations)
    )
    fields: dict[str, FloatJet2] = {}
    for field_index, field in enumerate(BASE_FIELD_ORDER):
        by_group = {
            group: values[group_index * FIELD_COUNT + field_index]
            for group_index, group in enumerate(PARAMETER_GROUPS)
        }
        fields[field] = FloatJet2(
            value=by_group["u"],
            dt=by_group["p"],
            dr=by_group["q"],
            dtt=acceleration[field_index],
            dtr=by_group["p_r"],
            drr=by_group["q_r"],
        )
    try:
        return SphericalState(
            **fields,
            planck_mass=center.planck_mass,
            beta=center.beta,
            mu=center.mu,
            g4=center.g4,
            alpha=center.alpha,
            eta=center.eta,
            branch=center.branch,
        )
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise AccelerationSolveStop("invalid_state", "floating state violates the REF1 domain", {"error": str(exc)}) from exc


def _state_with_accelerations(
    state: SphericalState,
    accelerations: Sequence[Real | FloatFirstTangent],
) -> SphericalState:
    if len(accelerations) != FIELD_COUNT:
        raise ValueError("accelerations must contain six entries")
    fields = {
        field: float_jet(getattr(state, field), dtt=accelerations[index])
        for index, field in enumerate(BASE_FIELD_ORDER)
    }
    return SphericalState(
        **fields,
        planck_mass=state.planck_mass,
        beta=state.beta,
        mu=state.mu,
        g4=state.g4,
        alpha=state.alpha,
        eta=state.eta,
        branch=state.branch,
    )


def _residual_objects(
    state: SphericalState,
    accelerations: Sequence[Real | FloatFirstTangent],
    *,
    reference: ReferenceConnection,
    coordinate_radius: Fraction | int,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
) -> tuple[object, ...]:
    numeric_state = _state_with_accelerations(state, accelerations)
    return tuple(
        modified_harmonic_full_residuals(
            numeric_state,
            reference=reference,
            coordinate_radius=coordinate_radius,
            tilde_normal_factor=tilde_normal_factor,
            hat_normal_factor=hat_normal_factor,
        )["full_residual_vector"]
    )


def float_ref1_residual(
    state: SphericalState,
    accelerations: Sequence[Real],
    *,
    reference: ReferenceConnection,
    coordinate_radius: Fraction | int,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
) -> np.ndarray:
    acceleration = _vector("accelerations", accelerations)
    try:
        objects = _residual_objects(
            state,
            acceleration,
            reference=reference,
            coordinate_radius=coordinate_radius,
            tilde_normal_factor=tilde_normal_factor,
            hat_normal_factor=hat_normal_factor,
        )
        residual = np.asarray(
            [scalar_primal_and_tangent(value)[0] for value in objects],
            dtype=np.float64,
        )
    except (ArithmeticError, TypeError, ValueError) as exc:
        raise AccelerationSolveStop(
            "nonfinite_residual_or_jacobian",
            "floating REF1 residual evaluation failed",
            {"error": str(exc)},
        ) from exc
    if residual.shape != (FIELD_COUNT,) or not np.all(np.isfinite(residual)):
        raise AccelerationSolveStop(
            "nonfinite_residual_or_jacobian",
            "floating REF1 residual is nonfinite or has the wrong shape",
        )
    return residual


def acceleration_jacobian(
    state: SphericalState,
    accelerations: Sequence[Real],
    *,
    reference: ReferenceConnection,
    coordinate_radius: Fraction | int,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return residual, analytic acceleration Jacobian, and infinity condition."""

    acceleration = _vector("accelerations", accelerations)
    primals: list[np.ndarray] = []
    columns: list[np.ndarray] = []
    try:
        for column in range(FIELD_COUNT):
            seeded: list[float | FloatFirstTangent] = list(acceleration)
            seeded[column] = FloatFirstTangent.seed(acceleration[column])
            parts = tuple(
                scalar_primal_and_tangent(value)
                for value in _residual_objects(
                    state,
                    seeded,
                    reference=reference,
                    coordinate_radius=coordinate_radius,
                    tilde_normal_factor=tilde_normal_factor,
                    hat_normal_factor=hat_normal_factor,
                )
            )
            primals.append(np.asarray([item[0] for item in parts], dtype=np.float64))
            columns.append(np.asarray([item[1] for item in parts], dtype=np.float64))
    except (ArithmeticError, TypeError, ValueError) as exc:
        raise AccelerationSolveStop(
            "nonfinite_residual_or_jacobian",
            "floating REF1 acceleration differentiation failed",
            {"error": str(exc)},
        ) from exc
    residual = primals[0]
    primal_scale = max(1.0, float(np.linalg.norm(residual, ord=np.inf)))
    primal_tolerance = 64.0 * np.finfo(np.float64).eps * primal_scale
    if any(
        float(np.linalg.norm(residual - other, ord=np.inf)) > primal_tolerance
        for other in primals[1:]
    ):
        raise AccelerationSolveStop(
            "nonfinite_residual_or_jacobian",
            "seeded acceleration columns do not reproduce one residual primal within the binary64 roundoff contract",
            {"primal_tolerance": primal_tolerance},
        )
    jacobian = np.column_stack(columns)
    if (
        residual.shape != (FIELD_COUNT,)
        or jacobian.shape != (FIELD_COUNT, FIELD_COUNT)
        or not np.all(np.isfinite(residual))
        or not np.all(np.isfinite(jacobian))
    ):
        raise AccelerationSolveStop(
            "nonfinite_residual_or_jacobian",
            "floating REF1 residual or acceleration Jacobian is nonfinite",
        )
    try:
        condition = float(np.linalg.cond(jacobian, p=np.inf))
    except np.linalg.LinAlgError as exc:
        raise AccelerationSolveStop(
            "linear_solve_failure", "kinetic condition estimate failed"
        ) from exc
    if not isfinite(condition):
        raise AccelerationSolveStop(
            "kinetic_condition_limit", "kinetic block is singular or nonfinite"
        )
    return residual, jacobian, condition


def finite_difference_acceleration_jacobian(
    state: SphericalState,
    accelerations: Sequence[Real],
    *,
    step: Real,
    reference: ReferenceConnection,
    coordinate_radius: Fraction | int,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
) -> np.ndarray:
    acceleration = _vector("accelerations", accelerations)
    width = _finite("finite-difference step", step)
    if width <= 0.0:
        raise ValueError("finite-difference step must be positive")
    columns = []
    for column in range(FIELD_COUNT):
        plus = acceleration.copy()
        minus = acceleration.copy()
        plus[column] += width
        minus[column] -= width
        columns.append(
            (
                float_ref1_residual(
                    state,
                    plus,
                    reference=reference,
                    coordinate_radius=coordinate_radius,
                    tilde_normal_factor=tilde_normal_factor,
                    hat_normal_factor=hat_normal_factor,
                )
                - float_ref1_residual(
                    state,
                    minus,
                    reference=reference,
                    coordinate_radius=coordinate_radius,
                    tilde_normal_factor=tilde_normal_factor,
                    hat_normal_factor=hat_normal_factor,
                )
            )
            / (2.0 * width)
        )
    answer = np.column_stack(columns)
    if not np.all(np.isfinite(answer)):
        raise AccelerationSolveStop(
            "nonfinite_residual_or_jacobian",
            "finite-difference acceleration Jacobian is nonfinite",
        )
    return answer


def _check_parameter_box(
    state: SphericalState,
    center: SphericalState,
    half_width: float,
) -> float:
    difference = np.abs(
        np.asarray(parameter_vector(state)) - np.asarray(parameter_vector(center))
    )
    maximum = float(np.max(difference))
    if maximum > half_width:
        raise AccelerationSolveStop(
            "parameter_box_violation",
            "state lies outside the authorized QIFT1 parameter box",
            {"maximum_parameter_displacement": maximum, "half_width": half_width},
        )
    return maximum


def solve_accelerations(
    state: SphericalState,
    *,
    reference: ReferenceConnection,
    coordinate_radius: Fraction | int,
    tilde_normal_factor: Fraction | int,
    hat_normal_factor: Fraction | int,
    branch_center: SphericalState,
    warm_start: Sequence[Real] | None = None,
    config: NonlinearSolverConfig | None = None,
) -> AccelerationSolveResult:
    """Safeguarded Newton solve on the already-proved QIFT1 local branch."""

    settings = config or NonlinearSolverConfig()
    if not isinstance(state, SphericalState) or not isinstance(branch_center, SphericalState):
        raise TypeError("state and branch_center must be SphericalState instances")
    if state.branch != "FGC-QR" or branch_center.branch != "FGC-QR":
        raise AccelerationSolveStop(
            "invalid_state", "SRC1 currently certifies the FGC-QR branch only"
        )
    _check_parameter_box(state, branch_center, settings.parameter_half_width)
    try:
        initial = (
            np.zeros(FIELD_COUNT, dtype=np.float64)
            if warm_start is None
            else _vector("warm_start", warm_start)
        )
    except (TypeError, ValueError) as exc:
        raise AccelerationSolveStop(
            "nonfinite_input",
            "warm start is not a finite six-component vector",
            {"error": str(exc)},
        ) from exc
    if float(np.max(np.abs(initial))) >= settings.acceleration_half_width:
        raise AccelerationSolveStop(
            "initial_acceleration_outside_authorized_box",
            "warm start lies outside the open authorized acceleration box",
            {
                "warm_start_infinity": float(np.max(np.abs(initial))),
                "acceleration_half_width": settings.acceleration_half_width,
            },
        )
    acceleration = initial.copy()
    history: list[NewtonIteration] = []
    last_condition = 0.0

    for iteration in range(settings.maximum_iterations + 1):
        residual, jacobian, condition = acceleration_jacobian(
            state,
            acceleration,
            reference=reference,
            coordinate_radius=coordinate_radius,
            tilde_normal_factor=tilde_normal_factor,
            hat_normal_factor=hat_normal_factor,
        )
        last_condition = condition
        residual_norm = float(np.linalg.norm(residual, ord=np.inf))
        if condition > settings.condition_number_maximum:
            raise AccelerationSolveStop(
                "kinetic_condition_limit",
                "kinetic condition number exceeds the declared maximum",
                {
                    "condition_infinity": condition,
                    "condition_maximum": settings.condition_number_maximum,
                    "iteration": iteration,
                },
            )
        displacement = float(np.linalg.norm(acceleration - initial, ord=np.inf))
        if residual_norm <= settings.residual_tolerance:
            if displacement > settings.normalized_branch_displacement_maximum:
                raise AccelerationSolveStop(
                    "branch_displacement_limit",
                    "converged root exceeds the branch-displacement limit",
                    {"displacement": displacement},
                )
            if float(np.max(np.abs(acceleration))) >= settings.acceleration_half_width:
                raise AccelerationSolveStop(
                    "final_acceleration_outside_authorized_box",
                    "converged root lies outside the open QIFT1 acceleration box",
                )
            return AccelerationSolveResult(
                acceleration_order=tuple(f"p_t.{field}" for field in BASE_FIELD_ORDER),
                accelerations=tuple(float(value) for value in acceleration),
                residual_vector=tuple(float(value) for value in residual),
                residual_infinity=residual_norm,
                iterations=iteration,
                jacobian_condition_infinity=last_condition,
                branch_displacement_infinity=displacement,
                iteration_history=tuple(history),
                converged=True,
                branch_continuity_preserved=True,
            )
        if iteration == settings.maximum_iterations:
            raise AccelerationSolveStop(
                "iteration_limit",
                "nonlinear acceleration solve exhausted the iteration limit",
                {"residual_infinity": residual_norm, "iteration": iteration},
            )
        try:
            direction = np.linalg.solve(jacobian, -residual)
        except np.linalg.LinAlgError as exc:
            raise AccelerationSolveStop(
                "linear_solve_failure", "Newton kinetic solve failed"
            ) from exc
        if not np.all(np.isfinite(direction)):
            raise AccelerationSolveStop(
                "linear_solve_failure", "Newton direction is nonfinite"
            )

        accepted: tuple[np.ndarray, float, float] | None = None
        step_fraction = 1.0
        for _ in range(settings.maximum_backtracks + 1):
            if step_fraction < settings.minimum_step_fraction:
                break
            trial = acceleration + step_fraction * direction
            if float(np.max(np.abs(trial))) >= settings.acceleration_half_width:
                step_fraction /= 2.0
                continue
            trial_residual = float_ref1_residual(
                state,
                trial,
                reference=reference,
                coordinate_radius=coordinate_radius,
                tilde_normal_factor=tilde_normal_factor,
                hat_normal_factor=hat_normal_factor,
            )
            trial_norm = float(np.linalg.norm(trial_residual, ord=np.inf))
            if trial_norm < residual_norm:
                accepted = trial, trial_norm, step_fraction
                break
            step_fraction /= 2.0
        if accepted is None:
            raise AccelerationSolveStop(
                "residual_failed_to_decrease",
                "no safeguarded Newton step strictly reduced the complete REF1 residual",
                {"iteration": iteration, "residual_infinity": residual_norm},
            )
        trial, trial_norm, accepted_fraction = accepted
        trial_displacement = float(np.linalg.norm(trial - initial, ord=np.inf))
        if trial_displacement > settings.normalized_branch_displacement_maximum:
            raise AccelerationSolveStop(
                "branch_displacement_limit",
                "accepted Newton step would cross the branch-displacement limit",
                {"displacement": trial_displacement},
            )
        history.append(
            NewtonIteration(
                iteration=iteration,
                accepted_step_fraction=accepted_fraction,
                residual_infinity_before=residual_norm,
                residual_infinity_after=trial_norm,
                jacobian_condition_infinity=condition,
                acceleration_infinity=float(np.linalg.norm(trial, ord=np.inf)),
                branch_displacement_infinity=trial_displacement,
            )
        )
        acceleration = trial

    raise AssertionError("unreachable nonlinear solver state")
