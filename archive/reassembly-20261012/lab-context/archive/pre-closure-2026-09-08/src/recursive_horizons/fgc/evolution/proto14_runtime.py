"""PROTO14 temporal-aware successor member runtime.

This module owns the smallest production adapter between the immutable
PROTO7 member semantics and TDG6's one/two/four same-grid temporal admission.
It neither changes the equations, source evaluator, projector, CFL policy,
nor the PROTO7 source-only retry rule.  It only replaces one accepted macro
step with an atomic sequence:

``prepare seven shadows -> durably classify/retry -> commit fine path``.

The runner owns campaign-level provenance, output namespaces, common events,
and checkpoint files.  This module owns the per-member temporal ledger and
the exact state required to snapshot or restore it.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import isfinite
from numbers import Real
from typing import Any, Callable, Mapping

import numpy as np

from .numerical_engine import EvolutionState, array_content_sha256
from .proto5_runtime import CFLRetryRequired, GR0RuntimeStop, GR0RuntimeStageTransaction
from .proto7_runtime import Proto7TerminalStop
from .tdg6_temporal_admission_runtime import (
    TDG6RefinementPathStop,
    TDG6TemporalLedger,
    TDG6TemporalRetryExhausted,
    TDG6TemporalRetryRequired,
    commit_tdg6_gr0_compositor,
    prepare_tdg6_gr0_compositor,
    require_tdg6_temporal_admission,
    restore_tdg6_checkpoint_extension,
    tdg6_checkpoint_extension,
)


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    if positive and answer <= 0.0:
        raise ValueError(f"{name} must be positive")
    return answer


class Proto14SourceRetryExhausted(RuntimeError):
    """The inherited PROTO7 source-retry recovery has exhausted its budget."""

    physical_classification = False

    def __init__(self, message: str, evidence: Mapping[str, object]) -> None:
        super().__init__(message)
        self.evidence = dict(evidence)


class Proto14InvalidTemporalRuntime(RuntimeError):
    """A TDG6 shadow failed outside a frozen recoverable ownership route."""

    physical_classification = False


@dataclass
class Proto14RunMember:
    """A PROTO7 member with an independently owned TDG6 temporal ledger.

    ``initial``, ``operator``, ``projector``, ``transaction``, and ``tracers``
    are deliberately preserved as opaque inherited objects.  The adapter
    consumes their existing PROTO7 semantics without reconstructing fields or
    altering the declared GR-0 equations.
    """

    amplitude: str
    method_label: str
    integrator_id: str
    spatial_order: int
    point_count: int
    input_hash: str
    initial: Any
    state: EvolutionState
    operator: Callable[[float, EvolutionState], Any]
    projector: Callable[[float, EvolutionState], EvolutionState] | None
    transaction: GR0RuntimeStageTransaction
    tracers: Any
    time: float = 0.0
    step_index: int = 0
    transaction_serial: int = 0
    CFL_retry_count: int = 0
    source_retry_count: int = 0
    temporal_ledger: TDG6TemporalLedger | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.state, EvolutionState):
            raise TypeError("PROTO14 member state must be EvolutionState")
        if not isinstance(self.transaction, GR0RuntimeStageTransaction):
            raise TypeError("PROTO14 member transaction must be GR0RuntimeStageTransaction")
        self.time = _finite("member time", self.time)
        if self.step_index < 0 or self.transaction_serial < 0:
            raise ValueError("PROTO14 member indices must be nonnegative")
        if self.CFL_retry_count < 0 or self.source_retry_count < 0:
            raise ValueError("PROTO14 retry counters must be nonnegative")
        if self.temporal_ledger is None:
            self.temporal_ledger = TDG6TemporalLedger.zero(initial_time=self.time)
        if not isinstance(self.temporal_ledger, TDG6TemporalLedger):
            raise TypeError("temporal_ledger must be TDG6TemporalLedger")
        if np.float64(self.temporal_ledger.last_accepted_time).tobytes() != np.float64(
            self.time
        ).tobytes():
            raise ValueError("PROTO14 temporal ledger time differs from member time")

    @property
    def key(self) -> str:
        return f"{self.method_label}-{self.point_count}"

    @classmethod
    def from_proto7(
        cls,
        member: object,
        *,
        temporal_ledger: TDG6TemporalLedger | None = None,
    ) -> "Proto14RunMember":
        """Copy one historical PROTO7 member without rebuilding its internals."""

        required = (
            "amplitude", "method_label", "integrator_id", "spatial_order",
            "point_count", "input_hash", "initial", "state", "operator",
            "projector", "transaction", "tracers", "time", "step_index",
            "transaction_serial", "CFL_retry_count", "source_retry_count",
        )
        missing = tuple(name for name in required if not hasattr(member, name))
        if missing:
            raise TypeError(f"PROTO7 member omits fields: {', '.join(missing)}")
        return cls(
            **{name: getattr(member, name) for name in required},
            temporal_ledger=temporal_ledger,
        )

    def snapshot(self) -> dict[str, Any]:
        """Capture inherited accepted state plus the TDG6-owned ledger."""

        return {
            "state": self.state,
            "time": self.time,
            "step_index": self.step_index,
            "transaction_serial": self.transaction_serial,
            "CFL_retry_count": self.CFL_retry_count,
            "source_retry_count": self.source_retry_count,
            "monitor_state": self.transaction.state,
            "causal_state": self.transaction.causal_state,
            "tracer_positions": self.tracers.positions.copy(),
            "tracer_proper_times": self.tracers.proper_times.copy(),
            "event_proper_times": [item.copy() for item in self.tracers.event_proper_times],
            "event_fields": [item.copy() for item in self.tracers.event_fields],
            "temporal_ledger": self.temporal_ledger,
        }

    def restore(self, snapshot: Mapping[str, Any]) -> None:
        """Restore an exact accepted boundary, including temporal evidence."""

        required = {
            "state", "time", "step_index", "transaction_serial", "CFL_retry_count",
            "source_retry_count", "monitor_state", "causal_state", "tracer_positions",
            "tracer_proper_times", "event_proper_times", "event_fields", "temporal_ledger",
        }
        if set(snapshot) != required:
            raise ValueError("PROTO14 snapshot keys differ")
        ledger = snapshot["temporal_ledger"]
        if not isinstance(ledger, TDG6TemporalLedger):
            raise TypeError("PROTO14 snapshot omits a TDG6 temporal ledger")
        self.state = snapshot["state"]
        self.time = _finite("snapshot time", snapshot["time"])
        self.step_index = int(snapshot["step_index"])
        self.transaction_serial = int(snapshot["transaction_serial"])
        self.CFL_retry_count = int(snapshot["CFL_retry_count"])
        self.source_retry_count = int(snapshot["source_retry_count"])
        self.transaction.state = snapshot["monitor_state"]
        self.transaction.causal_state = snapshot["causal_state"]
        self.tracers.positions = snapshot["tracer_positions"].copy()
        self.tracers.proper_times = snapshot["tracer_proper_times"].copy()
        self.tracers.event_proper_times = [item.copy() for item in snapshot["event_proper_times"]]
        self.tracers.event_fields = [item.copy() for item in snapshot["event_fields"]]
        self.temporal_ledger = ledger
        self.__post_init__()

    def _accepted_boundary_preserved(self, snapshot: Mapping[str, Any]) -> bool:
        return (
            array_content_sha256(self.state.u, self.state.p, self.state.q)
            == array_content_sha256(snapshot["state"].u, snapshot["state"].p, snapshot["state"].q)
            and self.time == snapshot["time"]
            and self.step_index == snapshot["step_index"]
            and self.transaction_serial == snapshot["transaction_serial"]
            and self.CFL_retry_count == snapshot["CFL_retry_count"]
            and self.source_retry_count == snapshot["source_retry_count"]
            and self.transaction.state == snapshot["monitor_state"]
            and self.transaction.causal_state == snapshot["causal_state"]
            and np.array_equal(self.tracers.positions, snapshot["tracer_positions"])
            and np.array_equal(self.tracers.proper_times, snapshot["tracer_proper_times"])
            and len(self.tracers.event_proper_times)
            == len(snapshot["event_proper_times"])
            and len(self.tracers.event_fields) == len(snapshot["event_fields"])
            and all(
                np.array_equal(actual, expected)
                for actual, expected in zip(
                    self.tracers.event_proper_times,
                    snapshot["event_proper_times"],
                    strict=True,
                )
            )
            and all(
                np.array_equal(actual, expected)
                for actual, expected in zip(
                    self.tracers.event_fields,
                    snapshot["event_fields"],
                    strict=True,
                )
            )
            and self.temporal_ledger == snapshot["temporal_ledger"]
        )

    def advance_to(
        self,
        target_time: float,
        *,
        retry_factor: float,
        maximum_CFL_retries: int,
        maximum_source_retries: int,
        minimum_step_size: float,
        rejected_trial_sink: Callable[[Mapping[str, Any]], None],
        temporal_rejection_sink: Callable[[Mapping[str, object]], None],
    ) -> None:
        """Advance through TDG6-admitted fine macro steps only.

        Source and CFL recovery retain their inherited ownership.  A TDG6
        rejection is durably serialized by the TDG6 sink before its ledger is
        adopted and before the next half-size proposal is considered.
        """

        if not callable(rejected_trial_sink) or not callable(temporal_rejection_sink):
            raise TypeError("PROTO14 rejection sinks must be callable")
        target = _finite("target_time", target_time)
        factor = _finite("retry_factor", retry_factor, positive=True)
        minimum = _finite("minimum_step_size", minimum_step_size, positive=True)
        if factor >= 1.0:
            raise ValueError("retry_factor must reduce retry step sizes")
        tolerance = 4096.0 * np.finfo(np.float64).eps * max(1.0, target)
        while self.time < target - tolerance:
            previous_speed = max(
                self.transaction.causal_state.previous_speed_upper, np.finfo(np.float64).tiny
            )
            proposed_size = min(
                target - self.time,
                self.transaction.cfl_maximum * self.initial.grid.spacing / previous_speed,
            )
            initial_proposed_size = proposed_size
            cfl_retries = 0
            source_retries = 0
            last_source_evidence: dict[str, Any] | None = None
            while True:
                if proposed_size < minimum:
                    if last_source_evidence is not None:
                        raise Proto14SourceRetryExhausted(
                            "PROTO14 inherited source retry fell below the frozen minimum step",
                            {
                                "exhaustion_kind": "minimum_step_size",
                                "minimum_step_size": minimum,
                                "attempted_step_size": proposed_size,
                                "source_retry_count_for_accepted_step": source_retries,
                                "last_rejected_proposal": last_source_evidence,
                            },
                        )
                    raise Proto14InvalidTemporalRuntime(
                        "adaptive macro step fell below the frozen minimum"
                    )
                boundary = self.snapshot()
                try:
                    prepared = prepare_tdg6_gr0_compositor(
                        method=self.integrator_id,
                        time=self.time,
                        step_size=proposed_size,
                        state=self.state,
                        rhs=self.operator,
                        projector=self.projector,
                        transaction=self.transaction,
                        tracers=self.tracers,
                        coordinates=self.initial.grid.coordinates,
                        temporal_ledger=self.temporal_ledger,
                        previous_step_index=self.step_index,
                        previous_transaction_serial=self.transaction_serial,
                    )
                except TDG6RefinementPathStop as stopped:
                    if not self._accepted_boundary_preserved(boundary):
                        raise RuntimeError("PROTO14 TDG6 shadow stop changed accepted boundary")
                    if stopped.source_retry is not None:
                        source_retries += 1
                        self.source_retry_count += 1
                        last_source_evidence = {
                            **asdict(stopped.source_retry),
                            "event_type": "rejected_unaccepted_source_only_proposal",
                            "amplitude": self.amplitude,
                            "member": self.key,
                            "method": self.method_label,
                            "point_count": self.point_count,
                            "TDG6_path": stopped.path,
                            "initial_step_size_for_accepted_step": initial_proposed_size,
                            "attempted_step_size": proposed_size,
                            "retry_count_for_accepted_step": source_retries,
                            "cumulative_member_source_retry_count": self.source_retry_count,
                            "target_common_event_time": target,
                            "external_transaction_state_preserved": True,
                            "complete_failure_evidence_serialized_before_retry": True,
                            "temporal_retry_owned_by_TDG6": False,
                        }
                        rejected_trial_sink(last_source_evidence)
                        if source_retries > maximum_source_retries:
                            raise Proto14SourceRetryExhausted(
                                "PROTO14 inherited source retry exceeded the frozen retry count",
                                {
                                    "exhaustion_kind": "maximum_source_retries",
                                    "maximum_source_retries": maximum_source_retries,
                                    "attempted_step_size": proposed_size,
                                    "source_retry_count_for_accepted_step": source_retries,
                                    "last_rejected_proposal": last_source_evidence,
                                },
                            )
                        proposed_size *= factor
                        continue
                    cause = stopped.cause
                    if isinstance(cause, CFLRetryRequired):
                        cfl_retries += 1
                        self.CFL_retry_count += 1
                        if cfl_retries > maximum_CFL_retries:
                            raise Proto14InvalidTemporalRuntime(
                                "frozen maximum CFL retries exceeded"
                            ) from cause
                        proposed_size *= factor
                        continue
                    if isinstance(cause, (Proto7TerminalStop, GR0RuntimeStop)):
                        raise cause
                    raise Proto14InvalidTemporalRuntime(
                        f"TDG6 shadow stopped outside a recoverable route: {type(cause).__name__}"
                    ) from cause

                try:
                    require_tdg6_temporal_admission(
                        prepared,
                        transaction=self.transaction,
                        tracers=self.tracers,
                        temporal_ledger=self.temporal_ledger,
                        current_time=self.time,
                        current_state=self.state,
                        current_step_index=self.step_index,
                        current_transaction_serial=self.transaction_serial,
                        durable_rejection_sink=temporal_rejection_sink,
                    )
                except TDG6TemporalRetryRequired as retry:
                    if not self._accepted_boundary_preserved(boundary):
                        raise RuntimeError("PROTO14 temporal retry changed accepted boundary")
                    self.temporal_ledger = retry.updated_ledger
                    proposed_size = retry.retry_step_size
                    continue
                except TDG6TemporalRetryExhausted as exhausted:
                    if not self._accepted_boundary_preserved(boundary):
                        raise RuntimeError("PROTO14 temporal exhaustion changed accepted boundary")
                    self.temporal_ledger = exhausted.updated_ledger
                    raise

                committed = commit_tdg6_gr0_compositor(
                    prepared,
                    transaction=self.transaction,
                    tracers=self.tracers,
                    temporal_ledger=self.temporal_ledger,
                    current_time=self.time,
                    current_state=self.state,
                    current_step_index=self.step_index,
                    current_transaction_serial=self.transaction_serial,
                )
                self.state = committed.state
                self.time = committed.time
                self.step_index = committed.step_index
                self.transaction_serial = committed.transaction_serial
                self.temporal_ledger = committed.temporal_ledger
                break
        if np.float64(self.time).tobytes() != np.float64(target).tobytes():
            if abs(self.time - target) > tolerance:
                raise RuntimeError("PROTO14 member did not land on the common event")
            self.time = target
            # Landing normalization is only allowed at an already-admitted
            # binary64-near target.  The ledger remains bitwise at its true
            # accepted endpoint, so a mismatch is a hard failure.
            if np.float64(self.temporal_ledger.last_accepted_time).tobytes() != np.float64(
                target
            ).tobytes():
                raise RuntimeError("PROTO14 landing normalization drifted TDG6 ledger time")


def proto14_checkpoint_extension(
    member: Proto14RunMember,
) -> tuple[dict[str, object], np.ndarray]:
    """Return the TDG6 checkpoint extension for one exact member boundary."""

    if not isinstance(member, Proto14RunMember):
        raise TypeError("member must be Proto14RunMember")
    return tdg6_checkpoint_extension(member.temporal_ledger)


def restore_proto14_checkpoint_extension(
    member: Proto14RunMember,
    metadata: Mapping[str, object],
    accumulated_debit: object,
) -> None:
    """Restore a TDG6 extension only if its accepted time matches the member."""

    if not isinstance(member, Proto14RunMember):
        raise TypeError("member must be Proto14RunMember")
    ledger = restore_tdg6_checkpoint_extension(metadata, accumulated_debit)
    if np.float64(ledger.last_accepted_time).tobytes() != np.float64(member.time).tobytes():
        raise ValueError("PROTO14 checkpoint temporal ledger time differs from member time")
    member.temporal_ledger = ledger


__all__ = [
    "Proto14InvalidTemporalRuntime",
    "Proto14RunMember",
    "Proto14SourceRetryExhausted",
    "proto14_checkpoint_extension",
    "restore_proto14_checkpoint_extension",
]
