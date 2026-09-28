"""Atomic composition of HLT1 stage health and BND2 causal bookkeeping.

Runge--Kutta internal stages are provisional.  A later stage or the explicit
candidate endpoint may fail even when earlier stage diagnostics pass.  This
adapter therefore evaluates HLT1 on a temporary monitor and BND2 on its pure
candidate ledger, then commits both only after the entire step is admissible.
On failure the physical endpoint and causal ledger remain unchanged while the
first HLT1 failure receipt is made permanent.

The transaction is reusable by NUM1 controls and the later COL1 engine.  It
does not manufacture any diagnostic values; callers must supply the complete
``StageHealthSnapshot`` from the actual trial state.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Mapping

from .boundary_domain import (
    BoundaryGeometry,
    CausalBudgetState,
    assess_causal_step,
)
from .health_monitor import (
    ClassicalHealthMonitor,
    ClassicalHealthStop,
    HealthMonitorState,
    StageHealthSnapshot,
)
from .numerical_engine import StageRecord, StepProposal


SpeedEvaluator = Callable[[StageRecord], float]
SnapshotFactory = Callable[[StageRecord, int, float], StageHealthSnapshot]


@dataclass(frozen=True, slots=True)
class RuntimeTransactionReceipt:
    accepted_stage_count: int
    first_transaction_serial: int
    last_transaction_serial: int
    interval_speed_upper: float
    remaining_causal_buffer: float
    strict_margin_over_required_buffer: float

    def as_mapping(self) -> Mapping[str, object]:
        return {
            "accepted_stage_count": self.accepted_stage_count,
            "first_transaction_serial": self.first_transaction_serial,
            "last_transaction_serial": self.last_transaction_serial,
            "interval_speed_upper": self.interval_speed_upper,
            "remaining_causal_buffer": self.remaining_causal_buffer,
            "strict_margin_over_required_buffer": self.strict_margin_over_required_buffer,
        }


class RuntimeStageTransaction:
    """Stateful all-or-nothing HLT1/BND2 proposal guard."""

    def __init__(
        self,
        *,
        health_monitor: ClassicalHealthMonitor,
        causal_state: CausalBudgetState,
        boundary_geometry: BoundaryGeometry,
        speed_evaluator: SpeedEvaluator,
        snapshot_factory: SnapshotFactory,
    ) -> None:
        if not isinstance(health_monitor, ClassicalHealthMonitor):
            raise TypeError("health_monitor must be ClassicalHealthMonitor")
        if not isinstance(causal_state, CausalBudgetState):
            raise TypeError("causal_state must be CausalBudgetState")
        if not isinstance(boundary_geometry, BoundaryGeometry):
            raise TypeError("boundary_geometry must be BoundaryGeometry")
        if not callable(speed_evaluator) or not callable(snapshot_factory):
            raise TypeError("runtime evaluators must be callable")
        self.health_monitor = health_monitor
        self.causal_state = causal_state
        self.boundary_geometry = boundary_geometry
        self.speed_evaluator = speed_evaluator
        self.snapshot_factory = snapshot_factory

    @staticmethod
    def _failed_state_from_original(
        original: HealthMonitorState,
        failure: ClassicalHealthStop,
    ) -> HealthMonitorState:
        return HealthMonitorState(
            accepted_stage_count=original.accepted_stage_count,
            last_accepted_time=original.last_accepted_time,
            last_accepted_stage_index=original.last_accepted_stage_index,
            first_failed_premise=failure.state.first_failed_premise,
            first_failed_stage_index=failure.state.first_failed_stage_index,
            first_failed_time=failure.state.first_failed_time,
        )

    def __call__(self, proposal: StepProposal) -> Mapping[str, object]:
        if not isinstance(proposal, StepProposal):
            raise TypeError("proposal must be StepProposal")
        if not proposal.stages or proposal.stages[-1].stage_name != "candidate_endpoint":
            raise ValueError("runtime transaction requires an explicit candidate endpoint")
        original_health = self.health_monitor.state
        if original_health.first_failed_premise is not None:
            raise ClassicalHealthStop(
                original_health, (original_health.first_failed_premise,)
            )

        speeds = tuple(float(self.speed_evaluator(stage)) for stage in proposal.stages)
        boundary = assess_causal_step(
            self.causal_state,
            self.boundary_geometry,
            trial_time=proposal.final_time,
            stage_coordinate_speed_uppers=speeds[:-1],
            candidate_endpoint_coordinate_speed_upper=speeds[-1],
        )

        temporary = ClassicalHealthMonitor(self.health_monitor.thresholds)
        temporary.state = original_health
        first_serial = original_health.last_accepted_stage_index + 1
        try:
            for offset, stage in enumerate(proposal.stages):
                snapshot = self.snapshot_factory(
                    stage,
                    first_serial + offset,
                    boundary.strict_margin_over_required_buffer,
                )
                temporary.assess_and_accept(snapshot)
        except ClassicalHealthStop as failure:
            # Roll back internal-stage accepts but retain the first failure.
            self.health_monitor.state = self._failed_state_from_original(
                original_health, failure
            )
            raise ClassicalHealthStop(self.health_monitor.state, failure.failures) from failure

        if not boundary.admissible:
            # A correct snapshot factory must pass the BND2 margin to HLT1,
            # whose final-priority boundary stop should already have fired.
            raise RuntimeError(
                "BND2 rejected a step that the composed HLT1 snapshot accepted"
            )

        self.health_monitor.state = temporary.state
        self.causal_state = boundary.candidate_state
        receipt = RuntimeTransactionReceipt(
            accepted_stage_count=len(proposal.stages),
            first_transaction_serial=first_serial,
            last_transaction_serial=first_serial + len(proposal.stages) - 1,
            interval_speed_upper=boundary.interval_speed_upper,
            remaining_causal_buffer=boundary.remaining_causal_buffer,
            strict_margin_over_required_buffer=boundary.strict_margin_over_required_buffer,
        )
        return receipt.as_mapping()


__all__ = [
    "RuntimeStageTransaction",
    "RuntimeTransactionReceipt",
    "SnapshotFactory",
    "SpeedEvaluator",
]
