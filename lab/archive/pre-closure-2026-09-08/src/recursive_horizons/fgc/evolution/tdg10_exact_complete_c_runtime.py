"""Pure TDG10 adapter from a prepared TDG6 compositor to exact all-18 complete-C.

This module consumes an already-prepared :class:`TDG6PreparedGR0Compositor`.
It never calls ``prepare_tdg6_gr0_compositor`` or
``commit_tdg6_gr0_compositor``, never constructs shadows, and never mutates a
transaction, ledger, tracer, monitor, causal, or PDE state.  The only
authorized tableau is SSPRK3 (``COMPARATOR_METHOD``) with exactly one outer,
two medium, and four fine accepted attempts.  Production owned-row count is
2044; synthetic tests may pass ``expected_owned_row_count``.

For each channel in ``TDG6_COMPLETE_STATE_CHANNELS`` order the adapter builds
one refinement-row stream from the actual finite endpoint values and the
fresh endpoint RHS records of those seven attempts.  A segment is the plain
five-tuple ``(left_value, left_rhs, right_value, right_rhs, width)`` of
built-in finite binary64 floats.  Owned-domain slicing is the TDG6/TDG5
interior slice, not the full grid: the centre row and the four
projector-owned outer rows are excluded.

Each of the 18 channels is assessed independently by the public
``assess_exact_complete_c_rows`` core with explicit candidate ceilings and
refinement depth.  All-of admission is exactly
``decision.admission_passed`` on every channel.  Compact evidence retains
channel identity, per-channel exact evidence, the failed-channel list, and
the finest-pair exact ``Fraction`` upper vector.  It does not retain prepared
shadows or candidate lists.

Resource exhaustion and route disagreement stay typed exact-core failures.
They are never rewritten as a physical classification or as a channel
admission decision.

Historical private helpers
--------------------------
TDG6's owned-domain slice and fresh endpoint-RHS extraction are underscore
helpers on the TDG5 runtime; no public equivalent exists.  This adapter
isolates exactly three of them:

* ``tdg5._owned_state``
* ``tdg5._owned_rhs``
* ``tdg5._proposal_endpoint_records``

They are the same helpers TDG6 already uses when enclosing its binary64
envelope.  This module does not import TDG7/8/9 binders, does not read the
prepared envelope decision, and does not authorize a successor commit.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from typing import Final

import numpy as np

from .numerical_engine import (
    COMPARATOR_METHOD,
    EvolutionRHS,
    EvolutionState,
    StepProposal,
    array_content_sha256,
)
from . import tdg5_stage_complete_refinement_runtime as tdg5
from .tdg6_temporal_admission_design import (
    TDG6_COMPLETE_STATE_CHANNELS,
    TDG6_LEVEL_STEP_COUNTS,
)
from .tdg6_temporal_admission_runtime import (
    TDG6PreparedGR0Compositor,
    TDG6ShadowPath,
)
from .tdg10_exact_complete_c_admission import (
    ExactCompleteCAdmissionEvidence,
    ExactCompleteCResourceExhausted,
    ExactCompleteCRouteDisagreement,
    assess_exact_complete_c_rows,
)


TDG10_PRODUCTION_OWNED_ROW_COUNT: Final[int] = 2044
TDG10_SSPRK3_STAGE_RECORD_COUNT: Final[int] = 4
TDG10_PRODUCTION_MAXIMUM_CANDIDATES_D01: Final[int] = 16352
TDG10_PRODUCTION_MAXIMUM_CANDIDATES_D12: Final[int] = 32704
TDG10_PRODUCTION_REFINEMENT_DEPTH: Final[int] = 160
_BLOCK_NAMES: Final[tuple[str, ...]] = ("u", "p", "q")
_FIELD_NAMES: Final[tuple[str, ...]] = ("alpha", "v", "lambda", "R", "phi", "chi")
_COMPLETE_STATE_SHAPE: Final[tuple[int, int]] = (3, 6)


def _positive_integer(value: object, *, name: str, allow_zero: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if value < (0 if allow_zero else 1):
        qualifier = "nonnegative" if allow_zero else "positive"
        raise ValueError(f"{name} must be {qualifier}")
    return value


def _builtin_finite_float(value: object, *, name: str) -> float:
    if type(value) is float:
        answer = value
    elif type(value) is np.float64:
        answer = float(value)
    else:
        raise TypeError(f"{name} must be a built-in binary64 float")
    if not math.isfinite(answer):
        raise ValueError(f"{name} must be finite")
    return answer


def _same_binary64(left: object, right: object) -> bool:
    left64 = _builtin_finite_float(left, name="left binary64 value")
    right64 = _builtin_finite_float(right, name="right binary64 value")
    return np.float64(left64).tobytes() == np.float64(right64).tobytes()


def _state_digest(state: EvolutionState) -> str:
    if not isinstance(state, EvolutionState):
        raise TypeError("state must be EvolutionState")
    return array_content_sha256(state.u, state.p, state.q)


def _rhs_digest(rhs: EvolutionRHS) -> str:
    if not isinstance(rhs, EvolutionRHS):
        raise TypeError("rhs must be EvolutionRHS")
    return array_content_sha256(rhs.du, rhs.dp, rhs.dq)


def _owned_state(state: EvolutionState) -> np.ndarray:
    """TDG6 owned-domain state slice; see the module docstring."""

    return tdg5._owned_state(state)


def _owned_rhs(rhs: EvolutionRHS) -> np.ndarray:
    """TDG6 owned-domain RHS slice; see the module docstring."""

    return tdg5._owned_rhs(rhs)


def _proposal_endpoint_records(
    proposal: StepProposal,
) -> tuple[EvolutionRHS, EvolutionRHS]:
    """Fresh first/last stage RHS records; see the module docstring."""

    return tdg5._proposal_endpoint_records(proposal)


def _subdivision_boundaries(
    start: float, width: float, count: int
) -> tuple[float, ...]:
    step = width / count
    boundaries = tuple(start + index * step for index in range(count + 1))
    expected_final = start + width
    if not _same_binary64(boundaries[-1], expected_final):
        raise ValueError("TDG10 subdivision does not share the macro endpoint")
    for left, right in zip(boundaries, boundaries[1:]):
        if not _same_binary64(right - left, step):
            raise ValueError("TDG10 subdivision is not bitwise uniform")
    return boundaries


def _require_owned_array(
    value: object, *, name: str, owned_row_count: int
) -> np.ndarray:
    if not isinstance(value, np.ndarray):
        raise TypeError(f"{name} must be a numpy array")
    if value.dtype != np.dtype("<f8") or not value.flags.c_contiguous:
        raise ValueError(f"{name} must be C-contiguous little-endian binary64")
    if value.shape != (_COMPLETE_STATE_SHAPE[0], owned_row_count, _COMPLETE_STATE_SHAPE[1]):
        raise ValueError(
            f"{name} component/order shape must be "
            f"(3, {owned_row_count}, 6), not {value.shape}"
        )
    if np.any(~np.isfinite(value)):
        raise ValueError(f"{name} must be finite")
    return value


def _segment_width(proposal: StepProposal) -> float:
    if not isinstance(proposal, StepProposal):
        raise TypeError("proposal must be StepProposal")
    start = _builtin_finite_float(proposal.initial_time, name="proposal.initial_time")
    end = _builtin_finite_float(proposal.final_time, name="proposal.final_time")
    width = end - start
    if type(width) is not float or not math.isfinite(width) or width <= 0.0:
        raise ValueError("proposal width must be a positive finite binary64 float")
    return width


def _validate_path_boundaries(
    path: TDG6ShadowPath,
    *,
    level: str,
    method: str,
    start: float,
    final: float,
    width: float,
    initial_state_sha256: str,
    expected_attempts: int,
) -> str:
    if not isinstance(path, TDG6ShadowPath):
        raise TypeError(f"{level} must be TDG6ShadowPath")
    if path.level != level:
        raise ValueError(f"{level} path has the wrong level name")
    if path.method != method:
        raise ValueError(f"{level} path is not SSPRK3")
    if len(path.attempts) != expected_attempts:
        raise ValueError(
            f"{level} path must contain exactly {expected_attempts} attempts"
        )
    if path.initial_state_sha256 != initial_state_sha256:
        raise ValueError("TDG10 paths do not share one bitwise initial state")
    if not (
        _same_binary64(path.initial_time, start)
        and _same_binary64(path.final_time, final)
    ):
        raise ValueError("TDG10 paths do not span one macro interval")
    boundaries = _subdivision_boundaries(start, width, expected_attempts)
    previous_digest: str | None = None
    previous_rhs_digest: str | None = None
    initial_rhs_digest: str | None = None
    for index, (attempt, left, right) in enumerate(
        zip(path.attempts, boundaries[:-1], boundaries[1:], strict=True)
    ):
        if attempt.accepted is None or attempt.retry is not None:
            raise ValueError("TDG10 path retained an unaccepted proposal")
        proposal = attempt.proposal
        if proposal.method != method:
            raise ValueError("TDG10 proposal is not SSPRK3")
        if len(proposal.stages) != TDG10_SSPRK3_STAGE_RECORD_COUNT:
            raise ValueError("TDG10 SSPRK3 proposal must retain 4 stage records")
        start_rhs, end_rhs = _proposal_endpoint_records(proposal)
        start_rhs_digest = _rhs_digest(start_rhs)
        end_rhs_digest = _rhs_digest(end_rhs)
        if index == 0:
            initial_rhs_digest = start_rhs_digest
        elif start_rhs_digest != previous_rhs_digest:
            raise ValueError(
                "TDG10 shadow proposal shared state produced different RHS fields"
            )
        if not (
            _same_binary64(proposal.initial_time, left)
            and _same_binary64(proposal.final_time, right)
        ):
            raise ValueError("TDG10 proposal boundary differs from the 1/2/4 subdivision")
        digest = _state_digest(proposal.initial_state)
        if index == 0 and digest != initial_state_sha256:
            raise ValueError("TDG10 path initial state differs from the compositor")
        if previous_digest is not None and digest != previous_digest:
            raise ValueError("TDG10 shadow proposal states are not contiguous")
        previous_digest = _state_digest(proposal.candidate_state)
        previous_rhs_digest = end_rhs_digest
        if _state_digest(attempt.accepted.state) != previous_digest:
            raise ValueError("TDG10 accepted state differs from the candidate")
    if initial_rhs_digest is None:
        raise ValueError("TDG10 path retained no initial RHS")
    return initial_rhs_digest


def _validate_prepared(
    prepared: TDG6PreparedGR0Compositor,
    *,
    expected_owned_row_count: int,
) -> tuple[float, float]:
    if not isinstance(prepared, TDG6PreparedGR0Compositor):
        raise TypeError("prepared must be TDG6PreparedGR0Compositor")
    if prepared.method != COMPARATOR_METHOD:
        raise ValueError("TDG10 exact complete-C runtime requires SSPRK3")
    start = _builtin_finite_float(prepared.initial_time, name="initial_time")
    final = _builtin_finite_float(prepared.final_time, name="final_time")
    width = final - start
    if type(width) is not float or not math.isfinite(width) or width <= 0.0:
        raise ValueError("macro width must be a positive finite binary64 float")
    if (
        len(prepared.outer.attempts),
        len(prepared.medium.attempts),
        len(prepared.fine.attempts),
    ) != TDG6_LEVEL_STEP_COUNTS:
        raise ValueError("TDG10 prepared paths must have 1/2/4 attempts")
    initial_rhs_digests: list[str] = []
    for path, level, expected_attempts in (
        (prepared.outer, "outer", 1),
        (prepared.medium, "medium", 2),
        (prepared.fine, "fine", 4),
    ):
        initial_rhs_digests.append(
            _validate_path_boundaries(
                path,
                level=level,
                method=prepared.method,
                start=start,
                final=final,
                width=width,
                initial_state_sha256=prepared.initial_state_sha256,
                expected_attempts=expected_attempts,
            )
        )
    if len(set(initial_rhs_digests)) != 1:
        raise ValueError("TDG10 paths' shared initial state produced different RHS fields")
    sample = _owned_state(prepared.outer.attempts[0].proposal.initial_state)
    _require_owned_array(
        sample, name="owned initial state", owned_row_count=expected_owned_row_count
    )
    return start, width


SegmentArrays = tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]


def _proposal_segment(
    proposal: StepProposal, *, owned_row_count: int, name: str
) -> SegmentArrays:
    start_rhs, end_rhs = _proposal_endpoint_records(proposal)
    left = _require_owned_array(
        _owned_state(proposal.initial_state),
        name=f"{name}.left",
        owned_row_count=owned_row_count,
    )
    left_rhs = _require_owned_array(
        _owned_rhs(start_rhs),
        name=f"{name}.left_rhs",
        owned_row_count=owned_row_count,
    )
    right = _require_owned_array(
        _owned_state(proposal.candidate_state),
        name=f"{name}.right",
        owned_row_count=owned_row_count,
    )
    right_rhs = _require_owned_array(
        _owned_rhs(end_rhs),
        name=f"{name}.right_rhs",
        owned_row_count=owned_row_count,
    )
    return left, left_rhs, right, right_rhs, _segment_width(proposal)


def _path_segments(
    path: TDG6ShadowPath, *, owned_row_count: int
) -> tuple[SegmentArrays, ...]:
    return tuple(
        _proposal_segment(
            attempt.proposal,
            owned_row_count=owned_row_count,
            name=f"{path.level}[{index}]",
        )
        for index, attempt in enumerate(path.attempts)
    )


def _scalar_segment(
    item: SegmentArrays,
    *,
    block: int,
    field: int,
    row: int,
    name: str,
) -> tuple[float, float, float, float, float]:
    left, left_rhs, right, right_rhs, width = item
    values = (
        _builtin_finite_float(left[block, row, field], name=f"{name}.left_value"),
        _builtin_finite_float(left_rhs[block, row, field], name=f"{name}.left_rhs"),
        _builtin_finite_float(right[block, row, field], name=f"{name}.right_value"),
        _builtin_finite_float(right_rhs[block, row, field], name=f"{name}.right_rhs"),
        _builtin_finite_float(width, name=f"{name}.width"),
    )
    if any(type(value) is not float for value in values):
        raise TypeError(f"{name} must contain built-in binary64 floats")
    return values


def _channel_rows(
    *,
    outer: tuple[SegmentArrays, ...],
    medium: tuple[SegmentArrays, ...],
    fine: tuple[SegmentArrays, ...],
    channel: str,
    owned_row_count: int,
    outer_width: float,
) -> tuple[tuple[object, object, object], ...]:
    if channel not in TDG6_COMPLETE_STATE_CHANNELS:
        raise ValueError("unknown TDG10 channel")
    if len(outer) != 1 or len(medium) != 2 or len(fine) != 4:
        raise ValueError("TDG10 channel stream must be segmented 1/2/4")
    if any(not _same_binary64(item[4], outer_width / 2.0) for item in medium):
        raise ValueError("medium widths must equal half the outer width")
    if any(not _same_binary64(item[4], outer_width / 4.0) for item in fine):
        raise ValueError("fine widths must equal one quarter of the outer width")
    if not _same_binary64(outer[0][4], outer_width):
        raise ValueError("outer width differs from the macro interval")
    block_name, field_name = channel.split(":", 1)
    block = _BLOCK_NAMES.index(block_name)
    field = _FIELD_NAMES.index(field_name)
    rows: list[tuple[object, object, object]] = []
    for row in range(owned_row_count):
        rows.append(
            (
                _scalar_segment(
                    outer[0],
                    block=block,
                    field=field,
                    row=row,
                    name=f"{channel}[{row}].outer",
                ),
                tuple(
                    _scalar_segment(
                        item,
                        block=block,
                        field=field,
                        row=row,
                        name=f"{channel}[{row}].medium[{index}]",
                    )
                    for index, item in enumerate(medium)
                ),
                tuple(
                    _scalar_segment(
                        item,
                        block=block,
                        field=field,
                        row=row,
                        name=f"{channel}[{row}].fine[{index}]",
                    )
                    for index, item in enumerate(fine)
                ),
            )
        )
    if len(rows) != owned_row_count:
        raise ValueError("TDG10 channel omitted an owned row")
    return tuple(rows)


class TDG10ExactCompleteCRuntimeClosed(RuntimeError):
    """Fail-closed numerical stop; never a physical classification."""

    physical_classification = False

    def __init__(self, *, channel: str, cause: BaseException) -> None:
        if channel not in TDG6_COMPLETE_STATE_CHANNELS:
            raise ValueError("unknown TDG10 closed channel")
        if not isinstance(
            cause,
            (ExactCompleteCResourceExhausted, ExactCompleteCRouteDisagreement),
        ):
            raise TypeError("runtime closed cause must be an exact-core exception")
        self.channel = channel
        self.cause = cause
        self.evidence = cause.evidence
        super().__init__(
            f"tdg10_exact_complete_c_runtime_closed:{channel}:{type(cause).__name__}"
        )


@dataclass(frozen=True, slots=True)
class TDG10ExactCompleteCChannelEvidence:
    """Compact exact complete-C evidence for one complete-state channel."""

    channel: str
    evidence: ExactCompleteCAdmissionEvidence

    def __post_init__(self) -> None:
        if self.channel not in TDG6_COMPLETE_STATE_CHANNELS:
            raise ValueError("unknown TDG10 channel")
        if not isinstance(self.evidence, ExactCompleteCAdmissionEvidence):
            raise TypeError("evidence must be ExactCompleteCAdmissionEvidence")


@dataclass(frozen=True, slots=True)
class TDG10ExactCompleteCRuntimeEvidence:
    """All-18 exact complete-C diagnostic; never a commit or method award."""

    method: str
    owned_row_count: int
    channel_evidence: tuple[TDG10ExactCompleteCChannelEvidence, ...]
    complete_admission_passed: bool
    failed_channels: tuple[str, ...]
    finest_pair_exact_upper_vector: tuple[Fraction, ...]
    maximum_candidates_D01: int
    maximum_candidates_D12: int
    refinement_depth: int
    admission_is_all_of: bool = True
    state_advance_authorized: bool = False
    production_method_earned: bool = False
    PDE_state_committed: bool = False
    GR0_calibration_completed: bool = False
    candidate_branch_opened: bool = False

    def __post_init__(self) -> None:
        if self.method != COMPARATOR_METHOD:
            raise ValueError("TDG10 exact complete-C runtime requires SSPRK3")
        owned_row_count = _positive_integer(
            self.owned_row_count, name="owned_row_count"
        )
        _positive_integer(self.maximum_candidates_D01, name="maximum_candidates_D01")
        _positive_integer(self.maximum_candidates_D12, name="maximum_candidates_D12")
        _positive_integer(
            self.refinement_depth, name="refinement_depth", allow_zero=True
        )
        if len(self.channel_evidence) != len(TDG6_COMPLETE_STATE_CHANNELS):
            raise ValueError("TDG10 channel evidence is incomplete or reordered")
        channels = tuple(item.channel for item in self.channel_evidence)
        if channels != TDG6_COMPLETE_STATE_CHANNELS:
            raise ValueError("TDG10 channel evidence is incomplete or reordered")
        for item in self.channel_evidence:
            if not isinstance(item, TDG10ExactCompleteCChannelEvidence):
                raise TypeError(
                    "channel_evidence must contain TDG10ExactCompleteCChannelEvidence"
                )
            if item.evidence.row_count != owned_row_count:
                raise ValueError(
                    "exact complete-C row count differs from owned-row count"
                )
            if (
                item.evidence.maximum_candidates_D01 != self.maximum_candidates_D01
                or item.evidence.maximum_candidates_D12 != self.maximum_candidates_D12
                or item.evidence.refinement_depth != self.refinement_depth
            ):
                raise ValueError(
                    "exact complete-C ceilings differ from the runtime request"
                )
        expected_failed = tuple(
            item.channel
            for item in self.channel_evidence
            if not item.evidence.decision.admission_passed
        )
        if expected_failed != self.failed_channels:
            raise ValueError("TDG10 failed-channel ledger differs from channel evidence")
        if self.complete_admission_passed is not (not expected_failed):
            raise ValueError("TDG10 all-of admission differs from channel evidence")
        expected_upper = tuple(
            item.evidence.d12.upper for item in self.channel_evidence
        )
        if self.finest_pair_exact_upper_vector != expected_upper:
            raise ValueError(
                "finest-pair exact upper vector differs from channel evidence"
            )
        if len(self.finest_pair_exact_upper_vector) != len(TDG6_COMPLETE_STATE_CHANNELS):
            raise ValueError("finest-pair exact upper vector must have 18 channels")
        if (
            self.admission_is_all_of is not True
            or self.state_advance_authorized is not False
            or self.production_method_earned is not False
            or self.PDE_state_committed is not False
            or self.GR0_calibration_completed is not False
            or self.candidate_branch_opened is not False
        ):
            raise ValueError(
                "TDG10 exact complete-C runtime crossed its nonclaim boundary"
            )


def assess_tdg10_exact_complete_c(
    prepared: TDG6PreparedGR0Compositor,
    *,
    maximum_candidates_D01: int,
    maximum_candidates_D12: int,
    refinement_depth: int = 160,
    expected_owned_row_count: int = TDG10_PRODUCTION_OWNED_ROW_COUNT,
) -> TDG10ExactCompleteCRuntimeEvidence:
    """Convert a prepared SSPRK3 compositor into the exact all-18 diagnostic."""

    d01_ceiling = _positive_integer(
        maximum_candidates_D01, name="maximum_candidates_D01"
    )
    d12_ceiling = _positive_integer(
        maximum_candidates_D12, name="maximum_candidates_D12"
    )
    depth = _positive_integer(
        refinement_depth, name="refinement_depth", allow_zero=True
    )
    owned_row_count = _positive_integer(
        expected_owned_row_count, name="expected_owned_row_count"
    )
    if owned_row_count == TDG10_PRODUCTION_OWNED_ROW_COUNT and (
        d01_ceiling != TDG10_PRODUCTION_MAXIMUM_CANDIDATES_D01
        or d12_ceiling != TDG10_PRODUCTION_MAXIMUM_CANDIDATES_D12
        or depth != TDG10_PRODUCTION_REFINEMENT_DEPTH
    ):
        raise ValueError(
            "TDG10 production owned-row assessment requires the frozen "
            "D01/D12 candidate ceilings and refinement depth"
        )
    _, width = _validate_prepared(
        prepared, expected_owned_row_count=owned_row_count
    )
    outer = _path_segments(prepared.outer, owned_row_count=owned_row_count)
    medium = _path_segments(prepared.medium, owned_row_count=owned_row_count)
    fine = _path_segments(prepared.fine, owned_row_count=owned_row_count)

    channel_evidence: list[TDG10ExactCompleteCChannelEvidence] = []
    for channel in TDG6_COMPLETE_STATE_CHANNELS:
        rows = _channel_rows(
            outer=outer,
            medium=medium,
            fine=fine,
            channel=channel,
            owned_row_count=owned_row_count,
            outer_width=width,
        )
        try:
            evidence = assess_exact_complete_c_rows(
                rows,
                maximum_candidates_D01=d01_ceiling,
                maximum_candidates_D12=d12_ceiling,
                refinement_depth=depth,
            )
        except (
            ExactCompleteCResourceExhausted,
            ExactCompleteCRouteDisagreement,
        ) as exc:
            raise TDG10ExactCompleteCRuntimeClosed(
                channel=channel, cause=exc
            ) from exc
        channel_evidence.append(
            TDG10ExactCompleteCChannelEvidence(channel=channel, evidence=evidence)
        )

    failed = tuple(
        item.channel
        for item in channel_evidence
        if not item.evidence.decision.admission_passed
    )
    return TDG10ExactCompleteCRuntimeEvidence(
        method=prepared.method,
        owned_row_count=owned_row_count,
        channel_evidence=tuple(channel_evidence),
        complete_admission_passed=not failed,
        failed_channels=failed,
        finest_pair_exact_upper_vector=tuple(
            item.evidence.d12.upper for item in channel_evidence
        ),
        maximum_candidates_D01=d01_ceiling,
        maximum_candidates_D12=d12_ceiling,
        refinement_depth=depth,
    )


__all__ = [
    "TDG10ExactCompleteCChannelEvidence",
    "TDG10ExactCompleteCRuntimeClosed",
    "TDG10ExactCompleteCRuntimeEvidence",
    "TDG10_PRODUCTION_OWNED_ROW_COUNT",
    "TDG10_PRODUCTION_MAXIMUM_CANDIDATES_D01",
    "TDG10_PRODUCTION_MAXIMUM_CANDIDATES_D12",
    "TDG10_PRODUCTION_REFINEMENT_DEPTH",
    "TDG10_SSPRK3_STAGE_RECORD_COUNT",
    "assess_tdg10_exact_complete_c",
]
