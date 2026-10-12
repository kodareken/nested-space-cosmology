"""Prospective C1 compositor under synthetic qualification; no store authority.

The numerical proposals are the unchanged binary64 engine's proposals. Exact
accumulation-corrected Hermite paths are admission evidence only. The only
state returned by commit is the actual four-substep binary64 endpoint.

The enclosing future PROTO19/HLT17 protocol must authenticate source/operator
identity, the HLT15 origin, durable cursor and writer. This in-memory kernel
does not infer any of those facts from a caller-supplied callback or hash.
"""

from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass
from hashlib import sha256
import json
from math import isfinite
from numbers import Real
from typing import Callable, Mapping

import numpy as np

from recursive_horizons.evidence_io import canonical_json_bytes

from .numerical_engine import (
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
)
from .proto5_runtime import GR0RuntimeStageTransaction
from .tdg11_imp1_enclosure import (
    IMP1EnclosureLimits,
    IMP1FamilyAssessment,
    assess_imp1_family,
)
from .tdg11_imp1_ledger import (
    TDG11IMP1Ledger,
    TDG11IMP1TemporalRejection,
    accept_imp1_step,
    append_imp1_rejection,
    imp1_checkpoint_extension,
)
from .tdg11_msel1_reconstruction import (
    ORIGINAL_ARITHMETIC_ID,
    ValidatedRecordedFamily,
    validate_recorded_family,
)
from .tdg6_temporal_admission_design import (
    TDG6_COMPLETE_STATE_CHANNELS,
    TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP,
    TDG6_MINIMUM_MACRO_STEP,
)
from . import tdg6_temporal_admission_runtime as tdg6
from .tdg7_binary64_subdivision_lattice import (
    CoordinateLatticeLimitReached,
    TDG7Binary64SubdivisionPlan,
    plan_forward_proto14_subdivision,
    validate_tdg7_binary64_subdivision_plan,
)


ARTIFACT_ID = "FGC-1-TDG11-IMP1"
RUNTIME_SCHEMA_VERSION = 1
REJECTION_EVENT = "rejected_TDG11_IMP1_temporal_admission"
FRESH_VS_RETRY_DURABLE_PROTOCOL_IMPLEMENTED = False
CAMPAIGN_EXECUTION_AUTHORIZED = False
_STAGE_COUNTS = {PRIMARY_METHOD: 5, COMPARATOR_METHOD: 4}
_STAGE_NAMES = {
    PRIMARY_METHOD: ("rk4_k1", "rk4_k2", "rk4_k3", "rk4_k4", "candidate_endpoint"),
    COMPARATOR_METHOD: ("ssprk3_s0", "ssprk3_s1", "ssprk3_s2", "candidate_endpoint"),
}
_PREPARED_DOMAIN = b"TDG11-IMP1-PREPARED-BOUNDARY-v1\n"

# These helpers have no equivalent public clone/restore/shadow interfaces.
# Their defining files, and the transitive TDG5 transaction cloner, are pinned
# by the IMP1 certificate. No old module is patched or relabeled.
_clone_tracers = tdg6._clone_tracers
_restore_tracers = tdg6._restore_tracers
_tracer_snapshot = tdg6._tracer_snapshot
_shadow_path = tdg6._shadow_path


def _finite(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _same_bits(left: Real, right: Real) -> bool:
    return _finite("left", left).hex() == _finite("right", right).hex()


def _index(name: str, value: object) -> int:
    if type(value) is not int or value < 0:
        raise ValueError(f"{name} must be a nonnegative built-in integer")
    return value


def _state_hash(state: EvolutionState) -> str:
    if not isinstance(state, EvolutionState) or state.shape[1] != 6:
        raise ValueError("IMP1 state must contain all six frozen fields")
    return array_content_sha256(state.u, state.p, state.q)


def _coordinates(coordinates: object, state: EvolutionState) -> np.ndarray:
    result = np.ascontiguousarray(np.asarray(coordinates, dtype=np.float64))
    if (
        result.ndim != 1
        or result.size != state.shape[0]
        or result.size <= 5
        or not np.all(np.isfinite(result))
        or not np.all(result[1:] > result[:-1])
    ):
        raise ValueError("IMP1 coordinates must match the ordered finite state grid")
    return result


def _wire(value: object) -> object:
    """Bit-preserving receipt data, not an endpoint serialization facility."""

    if is_dataclass(value) and not isinstance(value, type):
        return {item.name: _wire(getattr(value, item.name)) for item in fields(value)}
    if isinstance(value, tuple):
        return [_wire(item) for item in value]
    if isinstance(value, list):
        return [_wire(item) for item in value]
    if isinstance(value, Mapping):
        return {key: _wire(item) for key, item in value.items()}
    if isinstance(value, np.generic):
        return {"numpy_scalar_dtype": str(value.dtype), "value": _wire(value.item())}
    if type(value) is float:
        return {"binary64_hex": _finite("receipt float", value).hex()}
    if value is None or type(value) in (bool, int, str):
        return value
    raise TypeError(f"unsupported IMP1 receipt value {type(value).__name__}")


def _digest(value: object) -> str:
    return sha256(canonical_json_bytes(_wire(value))).hexdigest()


def _guard_configuration(transaction: GR0RuntimeStageTransaction) -> tuple[object, ...]:
    return (
        transaction.thresholds,
        transaction.boundary_geometry,
        float(transaction.grid_spacing),
        float(transaction.cfl_maximum),
        float(transaction.hat_normal_factor),
    )


def _ledger_digest(ledger: TDG11IMP1Ledger) -> str:
    return sha256(canonical_json_bytes(imp1_checkpoint_extension(ledger))).hexdigest()


class IMP1CoordinateLatticeStop(CoordinateLatticeLimitReached):
    """A coordinate stop before any source/state/tracer hook is invoked."""

    physical_classification = False


def _plan(time: Real, event_target: Real, cap: Real) -> TDG7Binary64SubdivisionPlan:
    try:
        return plan_forward_proto14_subdivision(
            time, event_target, cap, minimum_width=float(TDG6_MINIMUM_MACRO_STEP)
        )
    except CoordinateLatticeLimitReached as error:
        raise IMP1CoordinateLatticeStop(error.reason) from error


class IMP1RefinementPathStop(RuntimeError):
    """Inherited source/health failure; not a C1 temporal rejection."""

    temporal_retry = False
    accepted_state_advanced = False
    real_transaction_advanced = False
    real_tracer_advanced = False

    def __init__(self, cause: tdg6.TDG6RefinementPathStop) -> None:
        self.path = cause.path
        self.cause = cause.cause
        self.source_retry = cause.source_retry
        self.source_retry_owned_by_PROTO7 = cause.source_retry is not None
        super().__init__(f"IMP1 inherited shadow path stop: {self.path}: {cause}")


class IMP1RuntimeResourceStop(RuntimeError):
    """Prospective kernel resource limit, not a temporal or physical nonpass."""

    physical_classification = False
    temporal_retry = False
    accepted_state_advanced = False

    def __init__(self, reason: str) -> None:
        self.reason = reason
        super().__init__(f"IMP1 runtime resource stop: {reason}")


@dataclass(frozen=True, slots=True, init=False)
class TDG11IMP1PreparedRuntime:
    """Factory-only receipt for the actual guarded 1/2/4 proposal family."""

    plan: TDG7Binary64SubdivisionPlan
    method: str
    original_ledger: TDG11IMP1Ledger
    original_monitor_state: object
    original_causal_state: object
    original_tracer_snapshot: object
    guard_configuration: tuple[object, ...]
    coordinates_sha256: str
    initial_state_sha256: str
    previous_step_index: int
    previous_transaction_serial: int
    outer: tdg6.TDG6ShadowPath
    medium: tdg6.TDG6ShadowPath
    fine: tdg6.TDG6ShadowPath
    recorded_family: ValidatedRecordedFamily
    assessment: IMP1FamilyAssessment
    limits: IMP1EnclosureLimits
    receipt_sha256: str
    original_rhs: Callable = field(repr=False, compare=False)
    original_projector: Callable | None = field(repr=False, compare=False)

    def __init__(self, *args: object, **kwargs: object) -> None:
        raise TypeError("use an IMP1 fresh or retry preparation entrypoint")

    @property
    def initial_time(self) -> float:
        return self.plan.current

    @property
    def final_time(self) -> float:
        return self.plan.endpoint

    @property
    def stage_record_count(self) -> int:
        return sum(
            path.stage_record_count for path in (self.outer, self.medium, self.fine)
        )


def _prepared_mapping(prepared: TDG11IMP1PreparedRuntime) -> dict[str, object]:
    return {
        "schema_version": RUNTIME_SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "method": prepared.method,
        "plan": _wire(prepared.plan),
        "original_ledger_sha256": _ledger_digest(prepared.original_ledger),
        "original_monitor_state": _wire(prepared.original_monitor_state),
        "original_causal_state": _wire(prepared.original_causal_state),
        "original_tracer_snapshot": _wire(prepared.original_tracer_snapshot),
        "guard_configuration": _wire(prepared.guard_configuration),
        "coordinates_sha256": prepared.coordinates_sha256,
        "initial_state_sha256": prepared.initial_state_sha256,
        "previous_step_index": prepared.previous_step_index,
        "previous_transaction_serial": prepared.previous_transaction_serial,
        "recorded_family_sha256": prepared.recorded_family.family_sha256,
        "assessment": prepared.assessment.as_mapping(),
        "limits": _wire(prepared.limits),
        "paths": [
            {
                "level": path.level,
                "final_monitor_state": _wire(path.final_monitor_state),
                "final_causal_state": _wire(path.final_causal_state),
                "final_tracer_snapshot": _wire(path.final_tracer_snapshot),
                "final_positions_sha256": array_content_sha256(
                    path.final_tracer_positions
                ),
                "final_proper_times_sha256": array_content_sha256(
                    path.final_tracer_proper_times
                ),
                "guard_receipts": [
                    _wire(attempt.accepted.guard_receipts) for attempt in path.attempts
                ],
                "stage_diagnostics": [
                    [_wire(stage.rhs.diagnostics) for stage in attempt.proposal.stages]
                    for attempt in path.attempts
                ],
            }
            for path in (prepared.outer, prepared.medium, prepared.fine)
        ],
        "corrected_endpoint_committable": False,
        "outer_or_medium_committable": False,
        "campaign_execution_authorized": False,
        "live_retry_callback_replacement_permitted": False,
    }


def _receipt_digest(prepared: TDG11IMP1PreparedRuntime) -> str:
    return sha256(
        _PREPARED_DOMAIN + canonical_json_bytes(_prepared_mapping(prepared))
    ).hexdigest()


def _validate_paths(prepared: TDG11IMP1PreparedRuntime) -> None:
    expected_records = _STAGE_COUNTS[prepared.method]
    for path, boundaries, count in zip(
        (prepared.outer, prepared.medium, prepared.fine),
        (prepared.plan.one_full, prepared.plan.two_half, prepared.plan.four_quarter),
        (1, 2, 4),
        strict=True,
    ):
        if not isinstance(path, tdg6.TDG6ShadowPath):
            raise TypeError("IMP1 path must be a guarded TDG6ShadowPath data record")
        path.__post_init__()
        if (
            path.method != prepared.method
            or len(path.attempts) != count
            or path.initial_state_sha256 != prepared.initial_state_sha256
            or path.stage_record_count != count * expected_records
            or not _same_bits(path.initial_time, prepared.plan.current)
            or not _same_bits(path.final_time, prepared.plan.endpoint)
        ):
            raise ValueError("IMP1 path differs from its method, origin or stage plan")
        for index, (attempt, left, right) in enumerate(
            zip(path.attempts, boundaries[:-1], boundaries[1:], strict=True)
        ):
            proposal = attempt.proposal
            midpoint = left + (right - left) / 2.0
            times = (
                (left, midpoint, midpoint, right, right)
                if prepared.method == PRIMARY_METHOD
                else (left, right, midpoint, right)
            )
            accepted = attempt.accepted
            if (
                accepted is None
                or not _same_bits(proposal.initial_time, left)
                or not _same_bits(proposal.final_time, right)
                or tuple(stage.stage_name for stage in proposal.stages)
                != _STAGE_NAMES[prepared.method]
                or any(
                    not _same_bits(stage.time, expected)
                    for stage, expected in zip(proposal.stages, times, strict=True)
                )
                or accepted.step_index != prepared.previous_step_index + index + 1
                or accepted.transaction_serial
                != prepared.previous_transaction_serial + (index + 1) * expected_records
                or not _same_bits(accepted.time, right)
                or accepted.stage_records is not proposal.stages
                or _state_hash(accepted.state) != _state_hash(proposal.candidate_state)
            ):
                raise ValueError(
                    "IMP1 accepted path record differs from its actual proposal"
                )
        if (
            path.final_monitor_state.accepted_stage_count
            != prepared.original_monitor_state.accepted_stage_count
            + count * expected_records
            or path.final_monitor_state.last_transaction_serial
            != prepared.previous_transaction_serial + count * expected_records - 1
            or not _same_bits(
                path.final_causal_state.accepted_time, prepared.plan.endpoint
            )
            or path.final_monitor_state.first_failed_premise is not None
        ):
            raise ValueError(
                "IMP1 shadow ledger differs from the guarded accepted stages"
            )


def validate_imp1_prepared(prepared: TDG11IMP1PreparedRuntime) -> None:
    """Recheck mutable proposal arrays and every retained adoption receipt."""

    if not isinstance(prepared, TDG11IMP1PreparedRuntime):
        raise TypeError("prepared must be TDG11IMP1PreparedRuntime")
    validate_tdg7_binary64_subdivision_plan(prepared.plan)
    if prepared.method not in _STAGE_COUNTS:
        raise ValueError("unknown IMP1 method")
    if not isinstance(prepared.original_ledger, TDG11IMP1Ledger):
        raise TypeError("IMP1 cannot adopt a TDG6 ledger as its current ledger")
    if not isinstance(prepared.assessment, IMP1FamilyAssessment):
        raise TypeError("IMP1 requires its own full-family assessment")
    _validate_paths(prepared)
    checked = validate_recorded_family(
        tuple(
            tuple(attempt.proposal for attempt in path.attempts)
            for path in (prepared.outer, prepared.medium, prepared.fine)
        ),
        method=prepared.method,
        arithmetic_id=ORIGINAL_ARITHMETIC_ID,
        expected_owned_row_count=prepared.recorded_family.owned_row_count,
    )
    if (
        checked.family_sha256 != prepared.recorded_family.family_sha256
        or prepared.assessment.family_sha256 != checked.family_sha256
        or prepared.assessment.method != prepared.method
        or _receipt_digest(prepared) != prepared.receipt_sha256
    ):
        raise ValueError("IMP1 proposal family or admission/adoption receipt changed")


def _check_current_boundary(
    *,
    method: str,
    time: Real,
    state: EvolutionState,
    transaction: GR0RuntimeStageTransaction,
    temporal_ledger: TDG11IMP1Ledger,
    previous_step_index: int,
    previous_transaction_serial: int,
) -> None:
    if not isinstance(transaction, GR0RuntimeStageTransaction):
        raise TypeError("transaction must be GR0RuntimeStageTransaction")
    if not isinstance(temporal_ledger, TDG11IMP1Ledger):
        raise TypeError("temporal_ledger must be TDG11IMP1Ledger")
    _index("previous_step_index", previous_step_index)
    _index("previous_transaction_serial", previous_transaction_serial)
    monitor = transaction.state
    count = _index("monitor accepted_stage_count", monitor.accepted_stage_count)
    if (
        type(monitor.last_transaction_serial) is not int
        or count != previous_transaction_serial
        or count != monitor.last_transaction_serial + 1
        or (count > 0 and not _same_bits(monitor.last_accepted_time, time))
    ):
        raise ValueError("IMP1 monitor count, serial or accepted time is inconsistent")
    failure = (
        monitor.first_failed_premise,
        monitor.first_failed_transaction_serial,
        monitor.first_failed_time,
    )
    if any(item is None for item in failure) and any(
        item is not None for item in failure
    ):
        raise ValueError("IMP1 monitor contains incomplete failure metadata")
    if (
        temporal_ledger.method != method
        or not _same_bits(temporal_ledger.last_accepted_time, time)
        or not _same_bits(transaction.causal_state.accepted_time, time)
        or temporal_ledger.current_state_sha256 != _state_hash(state)
        or temporal_ledger.current_step_index != previous_step_index
        or temporal_ledger.current_transaction_serial != previous_transaction_serial
        or previous_transaction_serial != transaction.state.last_transaction_serial + 1
    ):
        raise ValueError("IMP1 accepted runtime and ledger boundaries disagree")


def _revalidate_boundary(
    prepared: TDG11IMP1PreparedRuntime,
    *,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    temporal_ledger: TDG11IMP1Ledger,
    current_time: Real,
    current_state: EvolutionState,
    current_step_index: int,
    current_transaction_serial: int,
    coordinates: object,
) -> None:
    validate_imp1_prepared(prepared)
    _check_current_boundary(
        method=prepared.method,
        time=current_time,
        state=current_state,
        transaction=transaction,
        temporal_ledger=temporal_ledger,
        previous_step_index=current_step_index,
        previous_transaction_serial=current_transaction_serial,
    )
    if (
        not _same_bits(current_time, prepared.plan.current)
        or _state_hash(current_state) != prepared.initial_state_sha256
        or current_step_index != prepared.previous_step_index
        or current_transaction_serial != prepared.previous_transaction_serial
        or _ledger_digest(temporal_ledger) != _ledger_digest(prepared.original_ledger)
        or _digest(transaction.state) != _digest(prepared.original_monitor_state)
        or _digest(transaction.causal_state) != _digest(prepared.original_causal_state)
        or _digest(_guard_configuration(transaction))
        != _digest(prepared.guard_configuration)
        or _tracer_snapshot(tracers) != prepared.original_tracer_snapshot
        or array_content_sha256(_coordinates(coordinates, current_state))
        != prepared.coordinates_sha256
    ):
        raise ValueError("IMP1 prepared boundary is stale")


def _prepare(
    *,
    plan: TDG7Binary64SubdivisionPlan,
    method: str,
    state: EvolutionState,
    rhs: Callable[[float, EvolutionState], EvolutionRHS],
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    coordinates: object,
    temporal_ledger: TDG11IMP1Ledger,
    previous_step_index: int,
    previous_transaction_serial: int,
    limits: IMP1EnclosureLimits,
) -> TDG11IMP1PreparedRuntime:
    if method not in _STAGE_COUNTS:
        raise ValueError("unknown IMP1 method")
    if not callable(rhs) or (projector is not None and not callable(projector)):
        raise TypeError("invalid IMP1 source or projector callback")
    if not isinstance(limits, IMP1EnclosureLimits):
        raise TypeError("limits must be IMP1EnclosureLimits")
    _check_current_boundary(
        method=method,
        time=plan.current,
        state=state,
        transaction=transaction,
        temporal_ledger=temporal_ledger,
        previous_step_index=previous_step_index,
        previous_transaction_serial=previous_transaction_serial,
    )
    coordinate_array = _coordinates(coordinates, state).copy()
    coordinate_array.setflags(write=False)
    if state.shape[0] - 5 > limits.maximum_owned_rows:
        raise IMP1RuntimeResourceStop("maximum_owned_rows")
    state_hash = _state_hash(state)
    coordinate_hash = array_content_sha256(coordinate_array)
    monitor, causal = transaction.state, transaction.causal_state
    tracer_snapshot = _tracer_snapshot(tracers)
    guard_configuration = _guard_configuration(transaction)
    ledger_hash = _ledger_digest(temporal_ledger)

    def unchanged() -> bool:
        return (
            _state_hash(state) == state_hash
            and array_content_sha256(_coordinates(coordinates, state))
            == coordinate_hash
            and _digest(transaction.state) == _digest(monitor)
            and _digest(transaction.causal_state) == _digest(causal)
            and _digest(_guard_configuration(transaction))
            == _digest(guard_configuration)
            and _tracer_snapshot(tracers) == tracer_snapshot
            and _ledger_digest(temporal_ledger) == ledger_hash
        )

    paths: list[tdg6.TDG6ShadowPath] = []
    try:
        for level, boundaries in zip(
            ("outer", "medium", "fine"),
            (plan.one_full, plan.two_half, plan.four_quarter),
            strict=True,
        ):
            paths.append(
                _shadow_path(
                    level=level,
                    method=method,
                    boundaries=boundaries,
                    state=EvolutionState(state.u, state.p, state.q),
                    rhs=rhs,
                    projector=projector,
                    transaction=transaction,
                    tracers=tracers,
                    coordinates=coordinate_array,
                    previous_step_index=previous_step_index,
                    previous_transaction_serial=previous_transaction_serial,
                )
            )
        family = validate_recorded_family(
            tuple(
                tuple(attempt.proposal for attempt in path.attempts) for path in paths
            ),
            method=method,
            arithmetic_id=ORIGINAL_ARITHMETIC_ID,
            expected_owned_row_count=state.shape[0] - 5,
        )
        assessment = assess_imp1_family(family, limits=limits)
    except tdg6.TDG6RefinementPathStop as error:
        raise IMP1RefinementPathStop(error) from error
    finally:
        if not unchanged():
            raise RuntimeError("IMP1 preparation changed a real accepted boundary")

    values = dict(
        plan=plan,
        method=method,
        original_ledger=temporal_ledger,
        original_monitor_state=monitor,
        original_causal_state=causal,
        original_tracer_snapshot=tracer_snapshot,
        guard_configuration=guard_configuration,
        coordinates_sha256=coordinate_hash,
        initial_state_sha256=state_hash,
        previous_step_index=previous_step_index,
        previous_transaction_serial=previous_transaction_serial,
        outer=paths[0],
        medium=paths[1],
        fine=paths[2],
        recorded_family=family,
        assessment=assessment,
        limits=limits,
        original_rhs=rhs,
        original_projector=projector,
    )
    prepared = object.__new__(TDG11IMP1PreparedRuntime)
    for name, value in values.items():
        object.__setattr__(prepared, name, value)
    object.__setattr__(prepared, "receipt_sha256", _receipt_digest(prepared))
    validate_imp1_prepared(prepared)
    return prepared


def prepare_imp1_initial_runtime(
    *,
    method: str,
    time: Real,
    event_target: Real,
    requested_cap: Real,
    state: EvolutionState,
    rhs: Callable[[float, EvolutionState], EvolutionRHS],
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    coordinates: object,
    temporal_ledger: TDG11IMP1Ledger,
    previous_step_index: int,
    previous_transaction_serial: int,
    limits: IMP1EnclosureLimits = IMP1EnclosureLimits(),
) -> TDG11IMP1PreparedRuntime:
    """Prepare a fresh macro attempt; lattice validation is the first action."""

    plan = _plan(time, event_target, requested_cap)
    if not isinstance(temporal_ledger, TDG11IMP1Ledger):
        raise TypeError("temporal_ledger must be TDG11IMP1Ledger")
    if temporal_ledger.current_macro_step_temporal_retry_count != 0:
        raise ValueError("IMP1 fresh entrypoint cannot consume an active retry ledger")
    return _prepare(
        plan=plan,
        method=method,
        state=state,
        rhs=rhs,
        projector=projector,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinates,
        temporal_ledger=temporal_ledger,
        previous_step_index=previous_step_index,
        previous_transaction_serial=previous_transaction_serial,
        limits=limits,
    )


def _rejection_record(prepared: TDG11IMP1PreparedRuntime) -> dict[str, object]:
    """Identity/counter record; full classification is hash-bound separately."""

    ledger = prepared.original_ledger
    return TDG11IMP1TemporalRejection(
        event_type=REJECTION_EVENT,
        method=prepared.method,
        time_hex=prepared.plan.current.hex(),
        attempted_width_hex=prepared.plan.macro_width.hex(),
        next_cap_hex=(prepared.plan.macro_width / 2.0).hex(),
        accepted_state_sha256=prepared.initial_state_sha256,
        accepted_step_index=prepared.previous_step_index,
        accepted_transaction_serial=prepared.previous_transaction_serial,
        assessment_sha256=prepared.assessment.assessment_sha256,
        preparation_sha256=prepared.receipt_sha256,
        event_target_hex=prepared.plan.event_target.hex(),
        channel_order=TDG6_COMPLETE_STATE_CHANNELS,
        failed_channels=tuple(
            item.channel
            for item in prepared.assessment.channels
            if item.corrected.decision.admission_passed is not True
        ),
        cumulative_temporal_retry_count=ledger.cumulative_temporal_retry_count + 1,
        retry_count_for_current_macro_step=ledger.current_macro_step_temporal_retry_count
        + 1,
    ).as_mapping()


class IMP1TemporalRetryRequired(RuntimeError):
    """One durably acknowledged nonadmission, not a scientific verdict."""

    physical_classification = False
    accepted_state_advanced = False

    def __init__(
        self, record: Mapping[str, object], updated_ledger: TDG11IMP1Ledger
    ) -> None:
        self.record_bytes = canonical_json_bytes(dict(record))
        self.updated_ledger = updated_ledger
        self.retry_step_size = float.fromhex(str(record["next_cap_hex"]))
        super().__init__("IMP1 temporal admission requests a stage-safe half-cap retry")

    @property
    def evidence(self) -> dict[str, object]:
        return json.loads(self.record_bytes)


class IMP1TemporalRetryExhausted(RuntimeError):
    """Retain the final rejection; do not continue below the frozen limits."""

    physical_classification = False
    accepted_state_advanced = False

    def __init__(
        self,
        *,
        reason: str,
        record: Mapping[str, object],
        updated_ledger: TDG11IMP1Ledger,
    ) -> None:
        if reason not in {
            "maximum_temporal_retries",
            "minimum_macro_step",
            "coordinate_lattice",
        }:
            raise ValueError("unknown IMP1 retry exhaustion reason")
        self.reason = reason
        self.record_bytes = canonical_json_bytes(dict(record))
        self.updated_ledger = updated_ledger
        super().__init__(f"IMP1 temporal retry exhausted: {reason}")

    @property
    def evidence(self) -> dict[str, object]:
        return json.loads(self.record_bytes)


def require_imp1_admission(
    prepared: TDG11IMP1PreparedRuntime,
    *,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    temporal_ledger: TDG11IMP1Ledger,
    current_time: Real,
    current_state: EvolutionState,
    current_step_index: int,
    current_transaction_serial: int,
    coordinates: object,
    durable_rejection_sink: Callable[[Mapping[str, object]], object] | None = None,
) -> None:
    """Require all channels; emit a retry only after the durable sink succeeds.

    The sink is an explicit enclosing-protocol obligation. This kernel checks
    order and immutable evidence, not disk durability by introspecting an
    arbitrary callable. A resource or source stop never reaches this function
    as a fabricated temporal classification.
    """

    current = dict(
        transaction=transaction,
        tracers=tracers,
        temporal_ledger=temporal_ledger,
        current_time=current_time,
        current_state=current_state,
        current_step_index=current_step_index,
        current_transaction_serial=current_transaction_serial,
        coordinates=coordinates,
    )
    _revalidate_boundary(prepared, **current)
    if prepared.assessment.admission_passed is True:
        return
    if not callable(durable_rejection_sink):
        raise TypeError("IMP1 nonadmission requires a durable rejection sink")
    record = _rejection_record(prepared)
    expected_bytes = canonical_json_bytes(record)
    # Preflight the complete new checkpoint before asking the sink to write.
    # Constructing an immutable value does not install it into a real member.
    updated = append_imp1_rejection(temporal_ledger, record)
    _ledger_digest(updated)
    durable_rejection_sink(record)
    if canonical_json_bytes(record) != expected_bytes:
        raise RuntimeError("IMP1 rejection sink altered the acknowledged evidence")
    _revalidate_boundary(prepared, **current)
    if (
        updated.current_macro_step_temporal_retry_count
        > TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
    ):
        raise IMP1TemporalRetryExhausted(
            reason="maximum_temporal_retries", record=record, updated_ledger=updated
        )
    if prepared.plan.macro_width / 2.0 < float(TDG6_MINIMUM_MACRO_STEP):
        raise IMP1TemporalRetryExhausted(
            reason="minimum_macro_step", record=record, updated_ledger=updated
        )
    raise IMP1TemporalRetryRequired(record, updated)


@dataclass(frozen=True, slots=True)
class IMP1RetrySuccessor:
    """Immutable preceding nonadmission and its exact next TDG7-safe plan."""

    predecessor: TDG11IMP1PreparedRuntime
    plan: TDG7Binary64SubdivisionPlan
    updated_ledger: TDG11IMP1Ledger
    rejection_bytes: bytes

    def __post_init__(self) -> None:
        validate_imp1_prepared(self.predecessor)
        validate_tdg7_binary64_subdivision_plan(self.plan)
        if self.predecessor.assessment.admission_passed is not False:
            raise ValueError("IMP1 cannot turn an admitted family into a retry")
        expected_record = _rejection_record(self.predecessor)
        if self.rejection_bytes != canonical_json_bytes(expected_record):
            raise ValueError("IMP1 retry evidence differs from its predecessor")
        expected_ledger = append_imp1_rejection(
            self.predecessor.original_ledger, expected_record
        )
        if _ledger_digest(self.updated_ledger) != _ledger_digest(expected_ledger):
            raise ValueError(
                "IMP1 retry ledger omits or changes the preceding rejection"
            )
        if (
            self.updated_ledger.current_macro_step_temporal_retry_count
            > TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
        ):
            raise ValueError("IMP1 retry exceeds the frozen retry count")
        expected_plan = _plan(
            self.predecessor.plan.current,
            self.predecessor.plan.event_target,
            self.predecessor.plan.macro_width / 2.0,
        )
        if (
            self.plan != expected_plan
            or self.plan.macro_width >= self.predecessor.plan.macro_width
        ):
            raise ValueError("IMP1 retry does not use the certified halved cap")


def replan_imp1_retry(
    prepared: TDG11IMP1PreparedRuntime, retry: IMP1TemporalRetryRequired
) -> IMP1RetrySuccessor:
    validate_imp1_prepared(prepared)
    if not isinstance(retry, IMP1TemporalRetryRequired):
        raise TypeError("retry must be IMP1TemporalRetryRequired")
    if retry.record_bytes != canonical_json_bytes(_rejection_record(prepared)):
        raise ValueError("IMP1 retry object differs from the prepared predecessor")
    try:
        plan = _plan(
            prepared.plan.current,
            prepared.plan.event_target,
            prepared.plan.macro_width / 2.0,
        )
    except IMP1CoordinateLatticeStop as error:
        raise IMP1TemporalRetryExhausted(
            reason="coordinate_lattice",
            record=retry.evidence,
            updated_ledger=retry.updated_ledger,
        ) from error
    return IMP1RetrySuccessor(prepared, plan, retry.updated_ledger, retry.record_bytes)


def prepare_imp1_retry_runtime(
    successor: IMP1RetrySuccessor,
    *,
    state: EvolutionState,
    rhs: Callable[[float, EvolutionState], EvolutionRHS],
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    coordinates: object,
) -> TDG11IMP1PreparedRuntime:
    """No caller replacement time, cap, ledger, method or indices are accepted."""

    if not isinstance(successor, IMP1RetrySuccessor):
        raise TypeError("successor must be IMP1RetrySuccessor")
    successor.__post_init__()
    prior = successor.predecessor
    if rhs is not prior.original_rhs or projector is not prior.original_projector:
        raise ValueError(
            "IMP1 live retry cannot replace its source or projector callback"
        )
    _revalidate_boundary(
        prior,
        transaction=transaction,
        tracers=tracers,
        temporal_ledger=prior.original_ledger,
        current_time=prior.plan.current,
        current_state=state,
        current_step_index=prior.previous_step_index,
        current_transaction_serial=prior.previous_transaction_serial,
        coordinates=coordinates,
    )
    return _prepare(
        plan=successor.plan,
        method=prior.method,
        state=state,
        rhs=rhs,
        projector=projector,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinates,
        temporal_ledger=successor.updated_ledger,
        previous_step_index=prior.previous_step_index,
        previous_transaction_serial=prior.previous_transaction_serial,
        limits=prior.limits,
    )


@dataclass(frozen=True, slots=True)
class TDG11IMP1CommittedRuntime:
    """In-memory fine endpoint; intentionally not an endpoint codec."""

    state: EvolutionState
    time: float
    step_index: int
    transaction_serial: int
    temporal_ledger: TDG11IMP1Ledger
    method: str
    assessment_sha256: str
    preparation_sha256: str
    committed_proposal_count: int = 4
    corrected_endpoint_adopted: bool = False
    outer_path_adopted: bool = False
    medium_path_adopted: bool = False
    campaign_execution_authorized: bool = False

    def __post_init__(self) -> None:
        if self.method not in _STAGE_COUNTS:
            raise ValueError("unknown IMP1 committed method")
        if (
            type(self.committed_proposal_count) is not int
            or self.committed_proposal_count != 4
            or self.corrected_endpoint_adopted is not False
            or self.outer_path_adopted is not False
            or self.medium_path_adopted is not False
            or self.campaign_execution_authorized is not False
        ):
            raise ValueError("IMP1 committed record crosses its fine-only kernel scope")
        if not isinstance(self.temporal_ledger, TDG11IMP1Ledger):
            raise TypeError("IMP1 commit requires its own ledger")
        if (
            self.temporal_ledger.method != self.method
            or not _same_bits(self.temporal_ledger.last_accepted_time, self.time)
            or self.temporal_ledger.current_state_sha256 != _state_hash(self.state)
            or self.temporal_ledger.current_step_index != self.step_index
            or self.temporal_ledger.current_transaction_serial
            != self.transaction_serial
        ):
            raise ValueError("IMP1 committed endpoint and ledger disagree")
        for digest in (self.assessment_sha256, self.preparation_sha256):
            if (
                type(digest) is not str
                or len(digest) != 64
                or any(char not in "0123456789abcdef" for char in digest)
            ):
                raise ValueError("invalid IMP1 commit receipt digest")


def commit_imp1_runtime(
    prepared: TDG11IMP1PreparedRuntime,
    *,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    temporal_ledger: TDG11IMP1Ledger,
    current_time: Real,
    current_state: EvolutionState,
    current_step_index: int,
    current_transaction_serial: int,
    coordinates: object,
) -> TDG11IMP1CommittedRuntime:
    """Adopt only the fine binary64 endpoint after complete revalidation.

    New ledger construction precedes real mutations. Late failure restores
    all monitor, causal and tracer state, including tracer event histories.
    The caller receives a new immutable state only after this commit succeeds.
    """

    _revalidate_boundary(
        prepared,
        transaction=transaction,
        tracers=tracers,
        temporal_ledger=temporal_ledger,
        current_time=current_time,
        current_state=current_state,
        current_step_index=current_step_index,
        current_transaction_serial=current_transaction_serial,
        coordinates=coordinates,
    )
    if prepared.assessment.admission_passed is not True:
        raise ValueError("IMP1 cannot commit an unadmitted family")
    accepted = prepared.fine.final_accepted
    updated = accept_imp1_step(
        temporal_ledger,
        final_time=accepted.time,
        state_sha256=_state_hash(accepted.state),
        step_index=accepted.step_index,
        transaction_serial=accepted.transaction_serial,
        debit_vector=prepared.assessment.public_fine_debits,
        assessment_sha256=prepared.assessment.assessment_sha256,
    )
    committed = TDG11IMP1CommittedRuntime(
        state=accepted.state,
        time=accepted.time,
        step_index=accepted.step_index,
        transaction_serial=accepted.transaction_serial,
        temporal_ledger=updated,
        method=prepared.method,
        assessment_sha256=prepared.assessment.assessment_sha256,
        preparation_sha256=prepared.receipt_sha256,
    )
    monitor, causal, prior_tracers = (
        transaction.state,
        transaction.causal_state,
        _clone_tracers(tracers),
    )
    guard_configuration = _guard_configuration(transaction)
    try:
        transaction.state = prepared.fine.final_monitor_state
        transaction.causal_state = prepared.fine.final_causal_state
        tracers.commit_advance(
            prepared.fine.final_tracer_positions,
            prepared.fine.final_tracer_proper_times,
        )
        # The tracer hook has run arbitrary callback code. Revalidate the
        # actual proposal arrays, not only the previously captured hash label.
        validate_imp1_prepared(prepared)
        if (
            _tracer_snapshot(tracers) != prepared.fine.final_tracer_snapshot
            or transaction.state != prepared.fine.final_monitor_state
            or transaction.causal_state != prepared.fine.final_causal_state
            or _digest(_guard_configuration(transaction))
            != _digest(guard_configuration)
            or _receipt_digest(prepared) != prepared.receipt_sha256
            or _state_hash(current_state) != prepared.initial_state_sha256
            or _state_hash(accepted.state) != updated.current_state_sha256
            or array_content_sha256(_coordinates(coordinates, current_state))
            != prepared.coordinates_sha256
        ):
            raise RuntimeError(
                "IMP1 adopted boundary differs from the guarded fine receipt"
            )
    except BaseException:
        transaction.state = monitor
        transaction.causal_state = causal
        (
            transaction.thresholds,
            transaction.boundary_geometry,
            transaction.grid_spacing,
            transaction.cfl_maximum,
            transaction.hat_normal_factor,
        ) = guard_configuration
        try:
            _restore_tracers(tracers, prior_tracers)
            if (
                _tracer_snapshot(tracers) != prepared.original_tracer_snapshot
                or transaction.state != monitor
                or transaction.causal_state != causal
            ):
                raise RuntimeError("restored IMP1 boundary differs")
        except BaseException as error:
            raise RuntimeError("IMP1 complete tracer/ledger rollback failed") from error
        raise
    return committed


__all__ = [
    "ARTIFACT_ID",
    "CAMPAIGN_EXECUTION_AUTHORIZED",
    "FRESH_VS_RETRY_DURABLE_PROTOCOL_IMPLEMENTED",
    "IMP1CoordinateLatticeStop",
    "IMP1RefinementPathStop",
    "IMP1RuntimeResourceStop",
    "IMP1TemporalRetryRequired",
    "IMP1TemporalRetryExhausted",
    "IMP1RetrySuccessor",
    "TDG11IMP1PreparedRuntime",
    "TDG11IMP1CommittedRuntime",
    "commit_imp1_runtime",
    "prepare_imp1_initial_runtime",
    "prepare_imp1_retry_runtime",
    "replan_imp1_retry",
    "require_imp1_admission",
    "validate_imp1_prepared",
]
