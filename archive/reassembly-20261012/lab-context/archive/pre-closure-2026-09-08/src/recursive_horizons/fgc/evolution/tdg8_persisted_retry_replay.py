"""Recovery-only replay of one exact persisted TDG7 retry predecessor.

The immutable TDG7 owner correctly rejects an active temporal ledger at its
fresh entrypoint.  This overlay supplies the narrower missing operation without
changing that certified module: rebuild TDG6's seven shadow paths through its
public compositor, validate them against an already reconstructed exact TDG7
plan, and expose only admission and retry replan.  There is deliberately no
commit surface.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from . import tdg6_temporal_admission_runtime as tdg6
from . import tdg7_stage_safe_runtime as tdg7
from .numerical_engine import EvolutionRHS, EvolutionState
from .proto5_runtime import GR0RuntimeStageTransaction
from .tdg6_temporal_admission_design import (
    TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP,
)
from .tdg7_binary64_subdivision_lattice import (
    TDG7Binary64SubdivisionPlan,
    validate_tdg7_binary64_subdivision_plan,
)


@dataclass(frozen=True, slots=True)
class TDG8PersistedRetryReplay:
    """Opaque, non-committable handle for one reconstructed shadow proposal."""

    _prepared: tdg7.TDG7PreparedStageSafeRuntime
    expected_prior_retry_count: int

    def __post_init__(self) -> None:
        tdg7.validate_tdg7_stage_safe_prepared(self._prepared)
        ledger = self._prepared.prepared_tdg6.original_temporal_ledger
        if (
            isinstance(self.expected_prior_retry_count, bool)
            or not isinstance(self.expected_prior_retry_count, int)
            or self.expected_prior_retry_count < 1
            or self.expected_prior_retry_count
            >= TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
            or ledger.current_macro_step_temporal_retry_count
            != self.expected_prior_retry_count
            or len(ledger.serialized_temporal_rejections)
            < self.expected_prior_retry_count
        ):
            raise ValueError("TDG8 persisted replay ledger count differs")

    @property
    def plan(self) -> TDG7Binary64SubdivisionPlan:
        return self._prepared.plan


def _checked(replay: TDG8PersistedRetryReplay) -> TDG8PersistedRetryReplay:
    if not isinstance(replay, TDG8PersistedRetryReplay):
        raise TypeError("replay must be TDG8PersistedRetryReplay")
    replay.__post_init__()
    return replay


def prepare_tdg8_persisted_retry_replay(
    *,
    replay_plan: TDG7Binary64SubdivisionPlan,
    expected_prior_retry_count: int,
    method: str,
    state: EvolutionState,
    rhs: Callable[[float, EvolutionState], EvolutionRHS],
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    coordinates: object,
    temporal_ledger: tdg6.TDG6TemporalLedger,
    previous_step_index: int,
    previous_transaction_serial: int,
) -> TDG8PersistedRetryReplay:
    """Prepare one exact active-ledger predecessor without a commit capability."""
    if not isinstance(replay_plan, TDG7Binary64SubdivisionPlan):
        raise TypeError("replay_plan must be a TDG7 plan")
    validate_tdg7_binary64_subdivision_plan(replay_plan)
    if (
        isinstance(expected_prior_retry_count, bool)
        or not isinstance(expected_prior_retry_count, int)
        or expected_prior_retry_count < 1
        or expected_prior_retry_count >= TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
    ):
        raise ValueError("TDG8 persisted replay retry count is outside its bound")
    if not isinstance(temporal_ledger, tdg6.TDG6TemporalLedger):
        raise TypeError("temporal_ledger must be TDG6TemporalLedger")
    if (
        temporal_ledger.current_macro_step_temporal_retry_count
        != expected_prior_retry_count
        or len(temporal_ledger.serialized_temporal_rejections)
        < expected_prior_retry_count
    ):
        raise ValueError("TDG8 persisted replay ledger count differs")

    inner = tdg6.prepare_tdg6_gr0_compositor(
        method=method,
        time=replay_plan.current,
        step_size=replay_plan.macro_width,
        state=state,
        rhs=rhs,
        projector=projector,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinates,
        temporal_ledger=temporal_ledger,
        previous_step_index=previous_step_index,
        previous_transaction_serial=previous_transaction_serial,
    )
    prepared = tdg7.TDG7PreparedStageSafeRuntime(replay_plan, inner)
    return TDG8PersistedRetryReplay(prepared, expected_prior_retry_count)


def require_tdg8_persisted_retry_replay_admission(
    replay: TDG8PersistedRetryReplay,
    **kwargs: object,
) -> None:
    """Apply unchanged TDG7/TDG6 admission to the replay-only preparation."""
    checked = _checked(replay)
    tdg7.require_tdg7_stage_safe_admission(checked._prepared, **kwargs)


def replan_tdg8_after_persisted_retry_rejection(
    replay: TDG8PersistedRetryReplay,
    retry: tdg6.TDG6TemporalRetryRequired,
) -> tdg7.TDG7RetrySuccessor:
    """Return the ordinary exact successor after the replayed rejection."""
    checked = _checked(replay)
    return tdg7.replan_tdg7_after_tdg6_retry(checked._prepared, retry)


__all__ = (
    "TDG8PersistedRetryReplay",
    "prepare_tdg8_persisted_retry_replay",
    "replan_tdg8_after_persisted_retry_rejection",
    "require_tdg8_persisted_retry_replay_admission",
)
