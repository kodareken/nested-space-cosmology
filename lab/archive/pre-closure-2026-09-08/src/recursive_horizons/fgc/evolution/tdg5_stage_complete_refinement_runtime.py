"""Runtime implementation of the frozen TDG5 temporal-refinement pair.

TDG5 compares one full Runge--Kutta step with two half steps on the same
spatial grid.  Both paths start from one bitwise-identical accepted state and
run through isolated copies of the inherited PROTO7 transaction.  The coarse
path is evidence only.  A prepared pair can expose the fine endpoint only
through :func:`commit_tdg5_gr0_refinement_pair`, which first verifies that the
accepted state and the real monitor/causal ledgers have not changed.

For every owned component of ``(u,p,q)``, the continuous numerical record is
the declared cubic Hermite extension formed from endpoint values and fresh
endpoint right-hand sides.  The implementation evaluates both endpoints and
all numerically real interior roots of the cubic derivative.  A separately
computed Bernstein convex-hull envelope, inflated by binary64 operation and
coefficient-construction debits, provides a fail-closed outward upper bound
even for a nearly degenerate derivative discriminant.

This module implements an instrument.  It defines no acceptance threshold,
does not reopen CAL10, and authorizes no PROTO14, GR-0 eligibility, SGB-L,
FGC-QR, DEF1, retained-EFT, transition, or physical claim.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import inf, isfinite, nextafter
from numbers import Real
from typing import Callable, Mapping

import numpy as np

from .boundary_domain import CausalBudgetState
from .numerical_engine import (
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    StepProposal,
    array_content_sha256,
)
from .proto5_runtime import (
    GR0RuntimeMonitorState,
    GR0RuntimeStageTransaction,
)
from .proto7_runtime import (
    ProposalFailureEvidence,
    Proto7StepAttempt,
    attempt_proto7_step,
)


TDG5_RUNTIME_METHODS: Mapping[str, tuple[str, int, int]] = {
    PRIMARY_METHOD: ("RK4", 4, 15),
    COMPARATOR_METHOD: ("SSPRK3", 3, 7),
}
TDG5_OWNED_INNER_ROW = 1
TDG5_PROJECTOR_OWNED_OUTER_ROWS = 4
_EPSILON = np.finfo(np.float64).eps
_SMALLEST_SUBNORMAL = np.nextafter(np.float64(0.0), np.float64(np.inf))


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


def _nonnegative_integer(name: str, value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _gamma(operation_count: int) -> float:
    if (
        isinstance(operation_count, bool)
        or not isinstance(operation_count, int)
        or operation_count < 1
    ):
        raise ValueError("operation_count must be a positive integer")
    product = operation_count * _EPSILON
    if product >= 1.0:
        raise ValueError("binary64 operation budget does not define gamma_n")
    return product / (1.0 - product)


def _same_binary64(left: float, right: float) -> bool:
    return np.float64(left).tobytes() == np.float64(right).tobytes()


def _state_hash(state: EvolutionState) -> str:
    if not isinstance(state, EvolutionState):
        raise TypeError("state must be EvolutionState")
    return array_content_sha256(state.u, state.p, state.q)


def _rhs_hash(rhs: EvolutionRHS) -> str:
    if not isinstance(rhs, EvolutionRHS):
        raise TypeError("rhs must be EvolutionRHS")
    return array_content_sha256(rhs.du, rhs.dp, rhs.dq)


def _owned_state(state: EvolutionState) -> np.ndarray:
    if state.shape[0] <= TDG5_OWNED_INNER_ROW + TDG5_PROJECTOR_OWNED_OUTER_ROWS:
        raise ValueError("state has no TDG5-owned interior rows")
    owned = slice(TDG5_OWNED_INNER_ROW, -TDG5_PROJECTOR_OWNED_OUTER_ROWS)
    return np.stack((state.u[owned], state.p[owned], state.q[owned]), axis=0)


def _owned_rhs(rhs: EvolutionRHS) -> np.ndarray:
    if rhs.shape[0] <= TDG5_OWNED_INNER_ROW + TDG5_PROJECTOR_OWNED_OUTER_ROWS:
        raise ValueError("RHS has no TDG5-owned interior rows")
    owned = slice(TDG5_OWNED_INNER_ROW, -TDG5_PROJECTOR_OWNED_OUTER_ROWS)
    return np.stack((rhs.du[owned], rhs.dp[owned], rhs.dq[owned]), axis=0)


def _linear_combination_with_radius(
    terms: tuple[tuple[float, np.ndarray, np.ndarray], ...],
) -> tuple[np.ndarray, np.ndarray]:
    if not terms:
        raise ValueError("linear combination requires at least one term")
    shape = terms[0][1].shape
    nominal = np.zeros(shape, dtype=np.float64)
    inherited = np.zeros(shape, dtype=np.float64)
    magnitude = np.zeros(shape, dtype=np.float64)
    for coefficient, value, radius in terms:
        scale = _finite("linear coefficient", coefficient)
        data = np.asarray(value, dtype=np.float64)
        uncertainty = np.asarray(radius, dtype=np.float64)
        if data.shape != shape or uncertainty.shape != shape:
            raise ValueError("linear-combination term shapes differ")
        if np.any(~np.isfinite(data)) or np.any(~np.isfinite(uncertainty)):
            raise ValueError("linear-combination terms must be finite")
        if np.any(uncertainty < 0.0):
            raise ValueError("coefficient radii must be nonnegative")
        nominal = nominal + scale * data
        inherited = inherited + abs(scale) * uncertainty
        magnitude = magnitude + abs(scale) * np.abs(data)
    operation_count = max(1, 2 * len(terms) - 1)
    operation_debit = _gamma(operation_count) * magnitude
    operation_debit = operation_debit + np.where(
        magnitude > 0.0,
        operation_count * _SMALLEST_SUBNORMAL,
        0.0,
    )
    radius = inherited + operation_debit
    if np.any(~np.isfinite(nominal)) or np.any(~np.isfinite(radius)):
        raise FloatingPointError("TDG5 coefficient construction overflowed")
    return nominal, radius


def _scaled_slope(
    width: float, rhs: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    value = width * rhs
    radius = _gamma(1) * np.abs(value)
    radius = radius + np.where(
        (width != 0.0) & (rhs != 0.0), _SMALLEST_SUBNORMAL, 0.0
    )
    if np.any(~np.isfinite(value)) or np.any(~np.isfinite(radius)):
        raise FloatingPointError("TDG5 scaled endpoint slope overflowed")
    return value, radius


def _hermite_coefficients_with_radius(
    left: np.ndarray,
    left_rhs: np.ndarray,
    right: np.ndarray,
    right_rhs: np.ndarray,
    *,
    width: float,
) -> tuple[np.ndarray, np.ndarray]:
    arrays = tuple(np.asarray(item, dtype=np.float64) for item in (
        left,
        left_rhs,
        right,
        right_rhs,
    ))
    if len({item.shape for item in arrays}) != 1:
        raise ValueError("Hermite endpoint arrays have different shapes")
    if any(np.any(~np.isfinite(item)) for item in arrays):
        raise ValueError("Hermite endpoint arrays must be finite")
    zero = np.zeros_like(arrays[0])
    left_slope, left_slope_radius = _scaled_slope(width, arrays[1])
    right_slope, right_slope_radius = _scaled_slope(width, arrays[3])
    coefficient_2, radius_2 = _linear_combination_with_radius(
        (
            (-3.0, arrays[0], zero),
            (-2.0, left_slope, left_slope_radius),
            (3.0, arrays[2], zero),
            (-1.0, right_slope, right_slope_radius),
        )
    )
    coefficient_3, radius_3 = _linear_combination_with_radius(
        (
            (2.0, arrays[0], zero),
            (1.0, left_slope, left_slope_radius),
            (-2.0, arrays[2], zero),
            (1.0, right_slope, right_slope_radius),
        )
    )
    coefficients = np.stack(
        (arrays[0], left_slope, coefficient_2, coefficient_3), axis=-1
    )
    radii = np.stack(
        (zero, left_slope_radius, radius_2, radius_3), axis=-1
    )
    return coefficients, radii


def _restrict_to_half_with_radius(
    coefficients: np.ndarray,
    radii: np.ndarray,
    *,
    half: int,
) -> tuple[np.ndarray, np.ndarray]:
    if isinstance(half, bool) or half not in (0, 1):
        raise ValueError("half must be zero or one")
    if coefficients.shape != radii.shape or coefficients.shape[-1] != 4:
        raise ValueError("cubic coefficient/radius shapes differ")
    terms = tuple((coefficients[..., index], radii[..., index]) for index in range(4))
    if half == 0:
        rows = (
            ((1.0, 0),),
            ((0.5, 1),),
            ((0.25, 2),),
            ((0.125, 3),),
        )
    else:
        rows = (
            ((1.0, 0), (0.5, 1), (0.25, 2), (0.125, 3)),
            ((0.5, 1), (0.5, 2), (0.375, 3)),
            ((0.25, 2), (0.375, 3)),
            ((0.125, 3),),
        )
    output = []
    output_radii = []
    for row in rows:
        value, radius = _linear_combination_with_radius(
            tuple((weight, terms[index][0], terms[index][1]) for weight, index in row)
        )
        output.append(value)
        output_radii.append(radius)
    return np.stack(output, axis=-1), np.stack(output_radii, axis=-1)


def _difference_with_radius(
    left: np.ndarray,
    left_radius: np.ndarray,
    right: np.ndarray,
    right_radius: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if not (
        left.shape == left_radius.shape == right.shape == right_radius.shape
        and left.shape[-1] == 4
    ):
        raise ValueError("difference coefficient/radius shapes differ")
    difference = left - right
    radius = (
        left_radius
        + right_radius
        + _gamma(1) * (np.abs(left) + np.abs(right))
    )
    if np.any(~np.isfinite(difference)) or np.any(~np.isfinite(radius)):
        raise FloatingPointError("TDG5 cubic difference overflowed")
    return difference, radius


@dataclass(frozen=True, slots=True)
class Binary64CubicEnvelope:
    polynomial_count: int
    endpoint_candidate_count: int
    real_interior_root_count: int
    ambiguous_discriminant_count: int
    raw_candidate_maximum: float
    coefficient_construction_debit: float
    outward_arithmetic_debit: float
    bernstein_certification_slack: float
    certified_continuous_upper_bound: float

    def __post_init__(self) -> None:
        for name in (
            "polynomial_count",
            "endpoint_candidate_count",
            "real_interior_root_count",
            "ambiguous_discriminant_count",
        ):
            _nonnegative_integer(name, getattr(self, name))
        for name in (
            "raw_candidate_maximum",
            "coefficient_construction_debit",
            "outward_arithmetic_debit",
            "bernstein_certification_slack",
            "certified_continuous_upper_bound",
        ):
            value = _finite(name, getattr(self, name))
            if value < 0.0:
                raise ValueError(f"{name} must be nonnegative")
        if self.endpoint_candidate_count != 2 * self.polynomial_count:
            raise ValueError("every cubic must retain both endpoints")
        if self.certified_continuous_upper_bound < self.raw_candidate_maximum:
            raise ValueError("continuous upper bound is below the candidate maximum")


def binary64_cubic_absolute_envelope(
    coefficients: object,
    coefficient_radii: object | None = None,
) -> Binary64CubicEnvelope:
    """Evaluate cubic extrema and return a fail-closed outward upper bound.

    ``coefficients[..., 4]`` defines the nominal binary64 cubics.  Optional
    nonnegative radii enclose coefficient-construction roundoff.  Stationary
    points are evaluated with a scaled, cancellation-resistant quadratic
    formula.  The Bernstein control polygon supplies an independent global
    enclosure; its gap from the stationary-point estimate remains public.
    """

    values = np.asarray(coefficients, dtype=np.float64)
    if values.ndim < 1 or values.shape[-1] != 4 or values.size == 0:
        raise ValueError("coefficients must have a nonempty final axis of four")
    if np.any(~np.isfinite(values)):
        raise ValueError("coefficients must be finite binary64 values")
    if coefficient_radii is None:
        radii = np.zeros_like(values)
    else:
        radii = np.asarray(coefficient_radii, dtype=np.float64)
        if radii.shape != values.shape:
            raise ValueError("coefficient radii must match coefficient shape")
        if np.any(~np.isfinite(radii)) or np.any(radii < 0.0):
            raise ValueError("coefficient radii must be finite and nonnegative")

    flat = values.reshape((-1, 4))
    flat_radii = radii.reshape((-1, 4))
    a, b, c, d = (flat[:, index] for index in range(4))
    candidate = np.maximum(np.abs(a), np.abs(a + b + c + d))

    # Scale the monomial coefficients by an exact power of two *before*
    # forming 3d and 2c.  This keeps every finite input away from derivative
    # overflow without changing its stationary points.
    derivative_base = np.maximum.reduce((np.abs(b), np.abs(c), np.abs(d)))
    nonconstant = derivative_base > 0.0
    _, exponent = np.frexp(derivative_base)
    scaled_b = np.zeros_like(b)
    scaled_c = np.zeros_like(c)
    scaled_d = np.zeros_like(d)
    scaled_b[nonconstant] = np.ldexp(b[nonconstant], -exponent[nonconstant])
    scaled_c[nonconstant] = np.ldexp(c[nonconstant], -exponent[nonconstant])
    scaled_d[nonconstant] = np.ldexp(d[nonconstant], -exponent[nonconstant])
    normalized_a = 3.0 * scaled_d
    normalized_b = 2.0 * scaled_c
    normalized_c = scaled_b

    linear = nonconstant & (normalized_a == 0.0) & (normalized_b != 0.0)
    roots = [np.full(flat.shape[0], np.nan), np.full(flat.shape[0], np.nan)]
    roots[0][linear] = -normalized_c[linear] / normalized_b[linear]

    quadratic = nonconstant & (normalized_a != 0.0)
    discriminant = normalized_b * normalized_b - 4.0 * normalized_a * normalized_c
    discriminant_debit = _gamma(5) * (
        normalized_b * normalized_b
        + 4.0 * np.abs(normalized_a * normalized_c)
    )
    ambiguous = quadratic & (np.abs(discriminant) <= discriminant_debit)
    real_quadratic = quadratic & (discriminant >= 0.0)
    square_root = np.zeros_like(discriminant)
    square_root[real_quadratic] = np.sqrt(discriminant[real_quadratic])
    q_value = np.zeros_like(discriminant)
    q_value[real_quadratic] = -0.5 * (
        normalized_b[real_quadratic]
        + np.copysign(
            square_root[real_quadratic], normalized_b[real_quadratic]
        )
    )
    nonzero_q = real_quadratic & (q_value != 0.0)
    roots[0][nonzero_q] = q_value[nonzero_q] / normalized_a[nonzero_q]
    roots[1][nonzero_q] = normalized_c[nonzero_q] / q_value[nonzero_q]
    repeated = real_quadratic & ~nonzero_q
    roots[0][repeated] = (
        -normalized_b[repeated] / (2.0 * normalized_a[repeated])
    )

    interior_root_count = 0
    with np.errstate(over="raise", invalid="raise"):
        for root in roots:
            interior = np.isfinite(root) & (root > 0.0) & (root < 1.0)
            interior_root_count += int(np.count_nonzero(interior))
            if np.any(interior):
                point = root[interior]
                evaluated = np.abs(
                    a[interior]
                    + point
                    * (
                        b[interior]
                        + point * (c[interior] + point * d[interior])
                    )
                )
                candidate[interior] = np.maximum(candidate[interior], evaluated)

    scale = np.sum(np.abs(flat), axis=1)
    coefficient_debit = np.sum(flat_radii, axis=1)
    arithmetic_debit = _gamma(16) * scale
    arithmetic_debit = arithmetic_debit + np.where(
        scale > 0.0, 16.0 * _SMALLEST_SUBNORMAL, 0.0
    )
    nonzero_candidate = candidate > 0.0
    ulp = np.zeros_like(candidate)
    ulp[nonzero_candidate] = (
        np.nextafter(candidate[nonzero_candidate], np.inf)
        - candidate[nonzero_candidate]
    )
    arithmetic_debit = arithmetic_debit + ulp

    bernstein = np.column_stack(
        (
            a,
            a + b / 3.0,
            a + 2.0 * b / 3.0 + c / 3.0,
            a + b + c + d,
        )
    )
    bernstein_nominal = np.max(np.abs(bernstein), axis=1)
    bernstein_roundoff = _gamma(16) * scale
    bernstein_roundoff = bernstein_roundoff + np.where(
        scale > 0.0, 16.0 * _SMALLEST_SUBNORMAL, 0.0
    )
    candidate_upper = candidate + arithmetic_debit + coefficient_debit
    bernstein_upper = bernstein_nominal + bernstein_roundoff + coefficient_debit
    certified = np.maximum(candidate_upper, bernstein_upper)
    positive = certified > 0.0
    certified[positive] = np.nextafter(certified[positive], np.inf)
    if np.any(~np.isfinite(certified)):
        raise FloatingPointError("TDG5 continuous envelope overflowed")

    raw_maximum = float(np.max(candidate, initial=0.0))
    coefficient_maximum = float(np.max(coefficient_debit, initial=0.0))
    arithmetic_maximum = float(np.max(arithmetic_debit, initial=0.0))
    certified_maximum = float(np.max(certified, initial=0.0))
    slack = max(
        0.0,
        certified_maximum
        - raw_maximum
        - coefficient_maximum
        - arithmetic_maximum,
    )
    return Binary64CubicEnvelope(
        polynomial_count=flat.shape[0],
        endpoint_candidate_count=2 * flat.shape[0],
        real_interior_root_count=interior_root_count,
        ambiguous_discriminant_count=int(np.count_nonzero(ambiguous)),
        raw_candidate_maximum=raw_maximum,
        coefficient_construction_debit=coefficient_maximum,
        outward_arithmetic_debit=arithmetic_maximum,
        bernstein_certification_slack=slack,
        certified_continuous_upper_bound=certified_maximum,
    )


@dataclass(frozen=True, slots=True)
class TDG5ContinuousDifferenceEvidence:
    method: str
    method_label: str
    formal_order: int
    richardson_denominator: int
    owned_row_count: int
    field_count: int
    left_half: Binary64CubicEnvelope
    right_half: Binary64CubicEnvelope
    raw_candidate_maximum: float
    certified_continuous_upper_bound: float
    conditional_fine_path_richardson_debit: float
    declared_cubic_is_exact_PDE_history: bool = False
    trajectory_asymptotic_regime_proved: bool = False
    rigorous_global_PDE_error_bound_proved: bool = False

    def __post_init__(self) -> None:
        if self.method not in TDG5_RUNTIME_METHODS:
            raise ValueError("unknown TDG5 runtime method")
        expected_label, expected_order, expected_denominator = TDG5_RUNTIME_METHODS[
            self.method
        ]
        if (
            self.method_label != expected_label
            or self.formal_order != expected_order
            or self.richardson_denominator != expected_denominator
        ):
            raise ValueError("TDG5 method-owned error semantics differ")
        if self.owned_row_count < 1 or self.field_count < 1:
            raise ValueError("TDG5 complete-state dimensions must be positive")
        if self.declared_cubic_is_exact_PDE_history is not False:
            raise ValueError("TDG5 cubic may not be promoted to exact PDE history")
        if (
            self.trajectory_asymptotic_regime_proved is not False
            or self.rigorous_global_PDE_error_bound_proved is not False
        ):
            raise ValueError("TDG5 numerical debit may not become a PDE theorem")
        expected = self.certified_continuous_upper_bound / self.richardson_denominator
        if self.conditional_fine_path_richardson_debit != expected:
            raise ValueError("TDG5 Richardson debit differs from its method owner")


def _proposal_endpoint_records(
    proposal: StepProposal,
) -> tuple[EvolutionRHS, EvolutionRHS]:
    if not isinstance(proposal, StepProposal) or not proposal.stages:
        raise TypeError("proposal must retain its stage records")
    first = proposal.stages[0]
    last = proposal.stages[-1]
    if last.stage_name != "candidate_endpoint":
        raise ValueError("proposal omits the fresh candidate-endpoint RHS")
    if not _same_binary64(first.time, proposal.initial_time):
        raise ValueError("proposal first RHS is not at the initial time")
    if not _same_binary64(last.time, proposal.final_time):
        raise ValueError("proposal endpoint RHS is not at the final time")
    if _state_hash(first.state) != _state_hash(proposal.initial_state):
        raise ValueError("proposal first-stage state differs from its initial state")
    if _state_hash(last.state) != _state_hash(proposal.candidate_state):
        raise ValueError("proposal endpoint state differs from its candidate state")
    return first.rhs, last.rhs


def assess_tdg5_continuous_difference(
    *,
    coarse: StepProposal,
    fine_left: StepProposal,
    fine_right: StepProposal,
) -> TDG5ContinuousDifferenceEvidence:
    """Compare complete coarse/fine Hermite records on both half intervals."""

    if not all(isinstance(item, StepProposal) for item in (coarse, fine_left, fine_right)):
        raise TypeError("TDG5 comparison requires three StepProposal values")
    if coarse.method not in TDG5_RUNTIME_METHODS:
        raise ValueError("unknown TDG5 method")
    if fine_left.method != coarse.method or fine_right.method != coarse.method:
        raise ValueError("TDG5 coarse and fine methods differ")
    if not _same_binary64(coarse.initial_time, fine_left.initial_time):
        raise ValueError("TDG5 paths do not share one initial time")
    if _state_hash(coarse.initial_state) != _state_hash(fine_left.initial_state):
        raise ValueError("TDG5 paths do not share one bitwise initial state")
    if not _same_binary64(fine_left.final_time, fine_right.initial_time):
        raise ValueError("TDG5 fine half steps are not contiguous")
    if _state_hash(fine_left.candidate_state) != _state_hash(fine_right.initial_state):
        raise ValueError("TDG5 fine midpoint state is not bitwise continuous")
    if not _same_binary64(coarse.final_time, fine_right.final_time):
        raise ValueError("TDG5 coarse and fine paths do not share one endpoint")
    coarse_width = coarse.final_time - coarse.initial_time
    left_width = fine_left.final_time - fine_left.initial_time
    right_width = fine_right.final_time - fine_right.initial_time
    if not (
        _same_binary64(left_width, coarse_width / 2.0)
        and _same_binary64(right_width, coarse_width / 2.0)
    ):
        raise ValueError("TDG5 fine paths are not two equal half steps")

    coarse_start_rhs, coarse_end_rhs = _proposal_endpoint_records(coarse)
    left_start_rhs, left_end_rhs = _proposal_endpoint_records(fine_left)
    right_start_rhs, right_end_rhs = _proposal_endpoint_records(fine_right)
    if _rhs_hash(coarse_start_rhs) != _rhs_hash(left_start_rhs):
        raise ValueError("TDG5 shared initial state produced different RHS fields")
    if _rhs_hash(left_end_rhs) != _rhs_hash(right_start_rhs):
        raise ValueError("TDG5 fine midpoint produced different endpoint/start RHS fields")

    coarse_coefficients, coarse_radii = _hermite_coefficients_with_radius(
        _owned_state(coarse.initial_state),
        _owned_rhs(coarse_start_rhs),
        _owned_state(coarse.candidate_state),
        _owned_rhs(coarse_end_rhs),
        width=coarse_width,
    )
    half_evidence: list[Binary64CubicEnvelope] = []
    for half, proposal, start_rhs, end_rhs in (
        (0, fine_left, left_start_rhs, left_end_rhs),
        (1, fine_right, right_start_rhs, right_end_rhs),
    ):
        restricted, restricted_radius = _restrict_to_half_with_radius(
            coarse_coefficients, coarse_radii, half=half
        )
        fine_coefficients, fine_radii = _hermite_coefficients_with_radius(
            _owned_state(proposal.initial_state),
            _owned_rhs(start_rhs),
            _owned_state(proposal.candidate_state),
            _owned_rhs(end_rhs),
            width=proposal.final_time - proposal.initial_time,
        )
        difference, difference_radius = _difference_with_radius(
            restricted,
            restricted_radius,
            fine_coefficients,
            fine_radii,
        )
        half_evidence.append(
            binary64_cubic_absolute_envelope(difference, difference_radius)
        )

    method_label, order, denominator = TDG5_RUNTIME_METHODS[coarse.method]
    raw = max(item.raw_candidate_maximum for item in half_evidence)
    certified = max(
        item.certified_continuous_upper_bound for item in half_evidence
    )
    return TDG5ContinuousDifferenceEvidence(
        method=coarse.method,
        method_label=method_label,
        formal_order=order,
        richardson_denominator=denominator,
        owned_row_count=coarse.initial_state.shape[0]
        - TDG5_OWNED_INNER_ROW
        - TDG5_PROJECTOR_OWNED_OUTER_ROWS,
        field_count=coarse.initial_state.shape[1],
        left_half=half_evidence[0],
        right_half=half_evidence[1],
        raw_candidate_maximum=raw,
        certified_continuous_upper_bound=certified,
        conditional_fine_path_richardson_debit=certified / denominator,
    )


class TDG5RefinementPathStop(RuntimeError):
    """One shadow path failed or requested a retry; no real state committed."""

    def __init__(
        self,
        *,
        path: str,
        cause: BaseException | None = None,
        retry: ProposalFailureEvidence | None = None,
    ) -> None:
        if path not in {"coarse", "fine_left", "fine_right"}:
            raise ValueError("unknown TDG5 refinement path")
        if (cause is None) == (retry is None):
            raise ValueError("TDG5 path stop requires exactly one cause or retry")
        self.path = path
        self.cause = cause
        self.retry = retry
        self.real_transaction_advanced = False
        self.accepted_state_advanced = False
        kind = "retry" if retry is not None else type(cause).__name__
        super().__init__(f"TDG5 {path} shadow path stopped: {kind}")


def _clone_gr0_transaction(
    transaction: GR0RuntimeStageTransaction,
) -> GR0RuntimeStageTransaction:
    if not isinstance(transaction, GR0RuntimeStageTransaction):
        raise TypeError("transaction must be GR0RuntimeStageTransaction")
    clone = GR0RuntimeStageTransaction(
        thresholds=transaction.thresholds,
        causal_state=transaction.causal_state,
        boundary_geometry=transaction.boundary_geometry,
        grid_spacing=transaction.grid_spacing,
        cfl_maximum=transaction.cfl_maximum,
        hat_normal_factor=transaction.hat_normal_factor,
    )
    clone.state = transaction.state
    return clone


def _shadow_attempt(
    *,
    path: str,
    method: str,
    time: float,
    step_size: float,
    state: EvolutionState,
    rhs: Callable[[float, EvolutionState], EvolutionRHS],
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    previous_step_index: int,
    previous_transaction_serial: int,
) -> Proto7StepAttempt:
    try:
        attempt = attempt_proto7_step(
            method=method,
            time=time,
            step_size=step_size,
            state=state,
            rhs=rhs,
            projector=projector,
            transaction=transaction,
            previous_step_index=previous_step_index,
            previous_transaction_serial=previous_transaction_serial,
            preaccept=None,
        )
    except BaseException as error:
        if isinstance(error, (KeyboardInterrupt, SystemExit)):
            raise
        raise TDG5RefinementPathStop(path=path, cause=error) from error
    if attempt.retry is not None:
        raise TDG5RefinementPathStop(path=path, retry=attempt.retry)
    if attempt.accepted is None:
        raise RuntimeError("TDG5 shadow attempt returned no accepted step")
    return attempt


@dataclass(frozen=True, slots=True)
class TDG5PreparedGR0Refinement:
    method: str
    initial_time: float
    final_time: float
    initial_state_sha256: str
    previous_step_index: int
    previous_transaction_serial: int
    original_monitor_state: GR0RuntimeMonitorState
    original_causal_state: CausalBudgetState
    coarse: Proto7StepAttempt
    fine_left: Proto7StepAttempt
    fine_right: Proto7StepAttempt
    fine_monitor_state: GR0RuntimeMonitorState
    fine_causal_state: CausalBudgetState
    continuous_difference: TDG5ContinuousDifferenceEvidence
    coarse_path_is_evidence_only: bool = True
    only_fine_path_is_committable: bool = True
    real_transaction_advanced_during_prepare: bool = False
    accepted_state_advanced_during_prepare: bool = False

    def __post_init__(self) -> None:
        if self.method not in TDG5_RUNTIME_METHODS:
            raise ValueError("unknown TDG5 prepared method")
        if len(self.initial_state_sha256) != 64:
            raise ValueError("initial_state_sha256 must be a SHA-256 digest")
        _nonnegative_integer("previous_step_index", self.previous_step_index)
        _nonnegative_integer(
            "previous_transaction_serial", self.previous_transaction_serial
        )
        if not _same_binary64(self.coarse.proposal.initial_time, self.initial_time):
            raise ValueError("prepared coarse path starts at the wrong time")
        if not _same_binary64(self.fine_right.proposal.final_time, self.final_time):
            raise ValueError("prepared fine path ends at the wrong time")
        if (
            self.coarse_path_is_evidence_only is not True
            or self.only_fine_path_is_committable is not True
            or self.real_transaction_advanced_during_prepare is not False
            or self.accepted_state_advanced_during_prepare is not False
        ):
            raise ValueError("TDG5 prepared pair violates shadow ownership")


def prepare_tdg5_gr0_refinement_pair(
    *,
    method: str,
    time: Real,
    step_size: Real,
    state: EvolutionState,
    rhs: Callable[[float, EvolutionState], EvolutionRHS],
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    previous_step_index: int,
    previous_transaction_serial: int,
) -> TDG5PreparedGR0Refinement:
    """Prepare coarse/fine shadow paths while leaving real state untouched."""

    if method not in TDG5_RUNTIME_METHODS:
        raise ValueError("unknown TDG5 runtime method")
    if not isinstance(state, EvolutionState):
        raise TypeError("state must be EvolutionState")
    if not callable(rhs):
        raise TypeError("rhs must be callable")
    if projector is not None and not callable(projector):
        raise TypeError("projector must be callable or None")
    if not isinstance(transaction, GR0RuntimeStageTransaction):
        raise TypeError("transaction must be GR0RuntimeStageTransaction")
    start = _finite("time", time)
    width = _finite("step_size", step_size, positive=True)
    step_index = _nonnegative_integer("previous_step_index", previous_step_index)
    serial = _nonnegative_integer(
        "previous_transaction_serial", previous_transaction_serial
    )
    if serial != transaction.state.last_transaction_serial + 1:
        raise ValueError("member serial does not follow the transaction ledger")
    if not _same_binary64(transaction.causal_state.accepted_time, start):
        raise ValueError("causal ledger time differs from the accepted state time")
    half = width / 2.0
    midpoint = start + half
    final = start + width
    if not _same_binary64(midpoint + half, final):
        raise ValueError("step is not bitwise divisible into two contiguous halves")

    original_monitor = transaction.state
    original_causal = transaction.causal_state
    initial_hash = _state_hash(state)

    coarse_transaction = _clone_gr0_transaction(transaction)
    coarse = _shadow_attempt(
        path="coarse",
        method=method,
        time=start,
        step_size=width,
        state=state,
        rhs=rhs,
        projector=projector,
        transaction=coarse_transaction,
        previous_step_index=step_index,
        previous_transaction_serial=serial,
    )

    fine_transaction = _clone_gr0_transaction(transaction)
    fine_left = _shadow_attempt(
        path="fine_left",
        method=method,
        time=start,
        step_size=half,
        state=state,
        rhs=rhs,
        projector=projector,
        transaction=fine_transaction,
        previous_step_index=step_index,
        previous_transaction_serial=serial,
    )
    left_accepted = fine_left.accepted
    if left_accepted is None:
        raise RuntimeError("TDG5 fine-left shadow did not accept")
    fine_right = _shadow_attempt(
        path="fine_right",
        method=method,
        time=midpoint,
        step_size=half,
        state=left_accepted.state,
        rhs=rhs,
        projector=projector,
        transaction=fine_transaction,
        previous_step_index=left_accepted.step_index,
        previous_transaction_serial=left_accepted.transaction_serial,
    )
    if transaction.state != original_monitor or transaction.causal_state != original_causal:
        raise RuntimeError("TDG5 shadow preparation mutated the real transaction")
    if _state_hash(state) != initial_hash:
        raise RuntimeError("TDG5 shadow preparation mutated the accepted state")

    continuous = assess_tdg5_continuous_difference(
        coarse=coarse.proposal,
        fine_left=fine_left.proposal,
        fine_right=fine_right.proposal,
    )
    return TDG5PreparedGR0Refinement(
        method=method,
        initial_time=start,
        final_time=final,
        initial_state_sha256=initial_hash,
        previous_step_index=step_index,
        previous_transaction_serial=serial,
        original_monitor_state=original_monitor,
        original_causal_state=original_causal,
        coarse=coarse,
        fine_left=fine_left,
        fine_right=fine_right,
        fine_monitor_state=fine_transaction.state,
        fine_causal_state=fine_transaction.causal_state,
        continuous_difference=continuous,
    )


@dataclass(frozen=True, slots=True)
class TDG5CommittedGR0Refinement:
    method: str
    time: float
    step_index: int
    transaction_serial: int
    state: EvolutionState
    continuous_difference: TDG5ContinuousDifferenceEvidence
    committed_path: str = "fine_two_half_steps"
    coarse_path_committed: bool = False
    fine_half_steps_committed_atomically: bool = True

    def __post_init__(self) -> None:
        if self.method not in TDG5_RUNTIME_METHODS:
            raise ValueError("unknown TDG5 committed method")
        _nonnegative_integer("step_index", self.step_index)
        _nonnegative_integer("transaction_serial", self.transaction_serial)
        if not isinstance(self.state, EvolutionState):
            raise TypeError("committed state must be EvolutionState")
        if (
            self.committed_path != "fine_two_half_steps"
            or self.coarse_path_committed is not False
            or self.fine_half_steps_committed_atomically is not True
        ):
            raise ValueError("TDG5 commit exposed the wrong path")


def commit_tdg5_gr0_refinement_pair(
    prepared: TDG5PreparedGR0Refinement,
    *,
    transaction: GR0RuntimeStageTransaction,
    current_time: Real,
    current_state: EvolutionState,
    current_step_index: int,
    current_transaction_serial: int,
) -> TDG5CommittedGR0Refinement:
    """Atomically adopt only the prepared fine transaction and endpoint."""

    if not isinstance(prepared, TDG5PreparedGR0Refinement):
        raise TypeError("prepared must be TDG5PreparedGR0Refinement")
    if not isinstance(transaction, GR0RuntimeStageTransaction):
        raise TypeError("transaction must be GR0RuntimeStageTransaction")
    time = _finite("current_time", current_time)
    step_index = _nonnegative_integer("current_step_index", current_step_index)
    serial = _nonnegative_integer(
        "current_transaction_serial", current_transaction_serial
    )
    if (
        not _same_binary64(time, prepared.initial_time)
        or _state_hash(current_state) != prepared.initial_state_sha256
        or step_index != prepared.previous_step_index
        or serial != prepared.previous_transaction_serial
        or transaction.state != prepared.original_monitor_state
        or transaction.causal_state != prepared.original_causal_state
    ):
        raise RuntimeError("TDG5 commit boundary drifted after shadow preparation")
    fine_accepted = prepared.fine_right.accepted
    if fine_accepted is None:
        raise RuntimeError("TDG5 prepared pair has no fine endpoint")
    prior_monitor = transaction.state
    prior_causal = transaction.causal_state
    try:
        transaction.state = prepared.fine_monitor_state
        transaction.causal_state = prepared.fine_causal_state
    except BaseException:
        transaction.state = prior_monitor
        transaction.causal_state = prior_causal
        raise
    return TDG5CommittedGR0Refinement(
        method=prepared.method,
        time=fine_accepted.time,
        step_index=fine_accepted.step_index,
        transaction_serial=fine_accepted.transaction_serial,
        state=fine_accepted.state,
        continuous_difference=prepared.continuous_difference,
    )


__all__ = [
    "Binary64CubicEnvelope",
    "TDG5CommittedGR0Refinement",
    "TDG5ContinuousDifferenceEvidence",
    "TDG5PreparedGR0Refinement",
    "TDG5RefinementPathStop",
    "TDG5_RUNTIME_METHODS",
    "assess_tdg5_continuous_difference",
    "binary64_cubic_absolute_envelope",
    "commit_tdg5_gr0_refinement_pair",
    "prepare_tdg5_gr0_refinement_pair",
]
