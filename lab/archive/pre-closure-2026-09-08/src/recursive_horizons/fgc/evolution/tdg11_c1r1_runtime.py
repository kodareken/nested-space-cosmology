"""C1R1 runtime adapter; not production method selection or authority.

This kernel prepares, admits, retries and commits the same binary64 1/2/4
shadows as IMP1, but assesses them with the C1R1 exact-representation
enclosure.  It does not call the sealed IMP1 ``_prepare`` and swap the
assessment: that path already spent the expensive work.

Returned assessments and ledger rejection records use the IMP1
reference-wire format.  Prepared receipts are C1R1-named and record
implementation provenance separately from that wire.  The only committable
endpoint remains the original binary64 fine state.

Private ``_imp1_reference_*`` aliases are compatibility dependencies on
sealed shadow, guard, ledger and wire helpers.  They are not a claim that
IMP1 assessed the family.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
import json
from numbers import Real
from typing import Callable, Mapping

from recursive_horizons.evidence_io import canonical_json_bytes

from .numerical_engine import (
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
)
from .proto5_runtime import GR0RuntimeStageTransaction
from .tdg11_c1r1_enclosure import (
    C1R1_IMPLEMENTATION_ID,
    C1R1_MATHEMATICAL_OBJECT,
    C1R1_REFERENCE_WIRE_EVALUATOR_ID,
    assess_c1r1_family,
)
from .tdg11_imp1_enclosure import (
    IMP1EnclosureLimits,
    IMP1FamilyAssessment,
)
from .tdg11_imp1_ledger import (
    IMP1_REJECTION_EVENT_TYPE,
    TDG11IMP1Ledger,
    TDG11IMP1TemporalRejection,
    accept_imp1_step,
    append_imp1_rejection,
)
from . import tdg11_imp1_runtime as _imp1_reference_runtime
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
from . import tdg6_temporal_admission_runtime as _tdg6_shadow_runtime
from .tdg7_binary64_subdivision_lattice import (
    CoordinateLatticeLimitReached,
    TDG7Binary64SubdivisionPlan,
    validate_tdg7_binary64_subdivision_plan,
)


ARTIFACT_ID = "FGC-1-TDG11-C1R1"
RUNTIME_SCHEMA_VERSION = 1
REFERENCE_WIRE_ARTIFACT_ID = _imp1_reference_runtime.ARTIFACT_ID
FRESH_VS_RETRY_DURABLE_PROTOCOL_IMPLEMENTED = False
CAMPAIGN_EXECUTION_AUTHORIZED = False
_PREPARED_DOMAIN = b"TDG11-C1R1-PREPARED-BOUNDARY-v1\n"

# Private compatibility dependencies.  These helpers have no C1R1-owned
# clone; the shadows, tracers, coordinate lattice and IMP1 ledger codec
# remain the sealed instruments.  Names record that fact.
_imp1_reference_stage_counts = _imp1_reference_runtime._STAGE_COUNTS
_imp1_reference_clone_tracers = _imp1_reference_runtime._clone_tracers
_imp1_reference_restore_tracers = _imp1_reference_runtime._restore_tracers
_imp1_reference_tracer_snapshot = _imp1_reference_runtime._tracer_snapshot
_imp1_reference_shadow_path = _imp1_reference_runtime._shadow_path
_imp1_reference_wire = _imp1_reference_runtime._wire
_imp1_reference_digest = _imp1_reference_runtime._digest
_imp1_reference_state_hash = _imp1_reference_runtime._state_hash
_imp1_reference_coordinates = _imp1_reference_runtime._coordinates
_imp1_reference_same_bits = _imp1_reference_runtime._same_bits
_imp1_reference_guard_configuration = _imp1_reference_runtime._guard_configuration
_imp1_reference_ledger_digest = _imp1_reference_runtime._ledger_digest
_imp1_reference_plan = _imp1_reference_runtime._plan
_imp1_reference_check_current_boundary = (
    _imp1_reference_runtime._check_current_boundary
)
_imp1_reference_validate_shadow_paths = _imp1_reference_runtime._validate_paths


class C1R1CoordinateLatticeStop(CoordinateLatticeLimitReached):
    """A coordinate stop before any source/state/tracer hook is invoked."""

    physical_classification = False


class C1R1RefinementPathStop(RuntimeError):
    """Inherited source/health failure; not a C1 temporal rejection."""

    temporal_retry = False
    accepted_state_advanced = False
    real_transaction_advanced = False
    real_tracer_advanced = False

    def __init__(self, cause: _tdg6_shadow_runtime.TDG6RefinementPathStop) -> None:
        self.path = cause.path
        self.cause = cause.cause
        self.source_retry = cause.source_retry
        self.source_retry_owned_by_PROTO7 = cause.source_retry is not None
        super().__init__(f"C1R1 inherited shadow path stop: {self.path}: {cause}")


class C1R1RuntimeResourceStop(RuntimeError):
    """Prospective kernel resource limit, not a temporal or physical nonpass."""

    physical_classification = False
    temporal_retry = False
    accepted_state_advanced = False

    def __init__(self, reason: str) -> None:
        self.reason = reason
        super().__init__(f"C1R1 runtime resource stop: {reason}")


def _plan(time: Real, event_target: Real, cap: Real) -> TDG7Binary64SubdivisionPlan:
    try:
        return _imp1_reference_plan(time, event_target, cap)
    except _imp1_reference_runtime.IMP1CoordinateLatticeStop as error:
        raise C1R1CoordinateLatticeStop(error.reason) from error


@dataclass(frozen=True, slots=True, init=False)
class TDG11C1R1PreparedRuntime:
    """Factory-only C1R1 receipt; assessment payload is IMP1 reference-wire."""

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
    outer: _tdg6_shadow_runtime.TDG6ShadowPath
    medium: _tdg6_shadow_runtime.TDG6ShadowPath
    fine: _tdg6_shadow_runtime.TDG6ShadowPath
    recorded_family: ValidatedRecordedFamily
    assessment: IMP1FamilyAssessment
    limits: IMP1EnclosureLimits
    implementation_id: str
    mathematical_object: str
    reference_wire_evaluator_id: str
    receipt_sha256: str
    original_rhs: Callable = field(repr=False, compare=False)
    original_projector: Callable | None = field(repr=False, compare=False)

    def __init__(self, *args: object, **kwargs: object) -> None:
        raise TypeError("use a C1R1 fresh or retry preparation entrypoint")

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


def _prepared_mapping(prepared: TDG11C1R1PreparedRuntime) -> dict[str, object]:
    return {
        "schema_version": RUNTIME_SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "implementation_id": prepared.implementation_id,
        "mathematical_object": prepared.mathematical_object,
        "reference_wire_artifact_id": REFERENCE_WIRE_ARTIFACT_ID,
        "reference_wire_evaluator_id": prepared.reference_wire_evaluator_id,
        "method": prepared.method,
        "plan": _imp1_reference_wire(prepared.plan),
        "original_ledger_sha256": _imp1_reference_ledger_digest(
            prepared.original_ledger
        ),
        "original_monitor_state": _imp1_reference_wire(
            prepared.original_monitor_state
        ),
        "original_causal_state": _imp1_reference_wire(prepared.original_causal_state),
        "original_tracer_snapshot": _imp1_reference_wire(
            prepared.original_tracer_snapshot
        ),
        "guard_configuration": _imp1_reference_wire(prepared.guard_configuration),
        "coordinates_sha256": prepared.coordinates_sha256,
        "initial_state_sha256": prepared.initial_state_sha256,
        "previous_step_index": prepared.previous_step_index,
        "previous_transaction_serial": prepared.previous_transaction_serial,
        "recorded_family_sha256": prepared.recorded_family.family_sha256,
        "assessment_reference_wire": prepared.assessment.as_mapping(),
        "limits": _imp1_reference_wire(prepared.limits),
        "paths": [
            {
                "level": path.level,
                "final_monitor_state": _imp1_reference_wire(path.final_monitor_state),
                "final_causal_state": _imp1_reference_wire(path.final_causal_state),
                "final_tracer_snapshot": _imp1_reference_wire(
                    path.final_tracer_snapshot
                ),
                "final_positions_sha256": array_content_sha256(
                    path.final_tracer_positions
                ),
                "final_proper_times_sha256": array_content_sha256(
                    path.final_tracer_proper_times
                ),
                "guard_receipts": [
                    _imp1_reference_wire(attempt.accepted.guard_receipts)
                    for attempt in path.attempts
                ],
                "stage_diagnostics": [
                    [
                        _imp1_reference_wire(stage.rhs.diagnostics)
                        for stage in attempt.proposal.stages
                    ]
                    for attempt in path.attempts
                ],
            }
            for path in (prepared.outer, prepared.medium, prepared.fine)
        ],
        "corrected_endpoint_committable": False,
        "outer_or_medium_committable": False,
        "campaign_execution_authorized": False,
        "live_retry_callback_replacement_permitted": False,
        "production_authority": False,
    }


def _receipt_digest(prepared: TDG11C1R1PreparedRuntime) -> str:
    return sha256(
        _PREPARED_DOMAIN + canonical_json_bytes(_prepared_mapping(prepared))
    ).hexdigest()


def validate_c1r1_prepared(prepared: TDG11C1R1PreparedRuntime) -> None:
    """Recheck mutable proposal arrays and the C1R1 adoption receipt."""

    if not isinstance(prepared, TDG11C1R1PreparedRuntime):
        raise TypeError("prepared must be TDG11C1R1PreparedRuntime")
    validate_tdg7_binary64_subdivision_plan(prepared.plan)
    if prepared.method not in _imp1_reference_stage_counts:
        raise ValueError("unknown C1R1 method")
    if not isinstance(prepared.original_ledger, TDG11IMP1Ledger):
        raise TypeError("C1R1 cannot adopt a TDG6 ledger as its current ledger")
    if not isinstance(prepared.assessment, IMP1FamilyAssessment):
        raise TypeError("C1R1 requires the IMP1 reference-wire family assessment")
    if (
        prepared.implementation_id != C1R1_IMPLEMENTATION_ID
        or prepared.mathematical_object != C1R1_MATHEMATICAL_OBJECT
        or prepared.reference_wire_evaluator_id != C1R1_REFERENCE_WIRE_EVALUATOR_ID
        or prepared.assessment.evaluator_id != C1R1_REFERENCE_WIRE_EVALUATOR_ID
    ):
        raise ValueError("C1R1 implementation or reference-wire identity changed")
    _imp1_reference_validate_shadow_paths(prepared)
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
        raise ValueError("C1R1 proposal family or admission/adoption receipt changed")


def _revalidate_boundary(
    prepared: TDG11C1R1PreparedRuntime,
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
    validate_c1r1_prepared(prepared)
    _imp1_reference_check_current_boundary(
        method=prepared.method,
        time=current_time,
        state=current_state,
        transaction=transaction,
        temporal_ledger=temporal_ledger,
        previous_step_index=current_step_index,
        previous_transaction_serial=current_transaction_serial,
    )
    if (
        not _imp1_reference_same_bits(current_time, prepared.plan.current)
        or _imp1_reference_state_hash(current_state) != prepared.initial_state_sha256
        or current_step_index != prepared.previous_step_index
        or current_transaction_serial != prepared.previous_transaction_serial
        or _imp1_reference_ledger_digest(temporal_ledger)
        != _imp1_reference_ledger_digest(prepared.original_ledger)
        or _imp1_reference_digest(transaction.state)
        != _imp1_reference_digest(prepared.original_monitor_state)
        or _imp1_reference_digest(transaction.causal_state)
        != _imp1_reference_digest(prepared.original_causal_state)
        or _imp1_reference_digest(
            _imp1_reference_guard_configuration(transaction)
        )
        != _imp1_reference_digest(prepared.guard_configuration)
        or _imp1_reference_tracer_snapshot(tracers)
        != prepared.original_tracer_snapshot
        or array_content_sha256(_imp1_reference_coordinates(coordinates, current_state))
        != prepared.coordinates_sha256
    ):
        raise ValueError("C1R1 prepared boundary is stale")


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
) -> TDG11C1R1PreparedRuntime:
    if method not in _imp1_reference_stage_counts:
        raise ValueError("unknown C1R1 method")
    if not callable(rhs) or (projector is not None and not callable(projector)):
        raise TypeError("invalid C1R1 source or projector callback")
    if not isinstance(limits, IMP1EnclosureLimits):
        raise TypeError("limits must be IMP1EnclosureLimits")
    _imp1_reference_check_current_boundary(
        method=method,
        time=plan.current,
        state=state,
        transaction=transaction,
        temporal_ledger=temporal_ledger,
        previous_step_index=previous_step_index,
        previous_transaction_serial=previous_transaction_serial,
    )
    coordinate_array = _imp1_reference_coordinates(coordinates, state).copy()
    coordinate_array.setflags(write=False)
    if state.shape[0] - 5 > limits.maximum_owned_rows:
        raise C1R1RuntimeResourceStop("maximum_owned_rows")
    state_hash = _imp1_reference_state_hash(state)
    coordinate_hash = array_content_sha256(coordinate_array)
    monitor, causal = transaction.state, transaction.causal_state
    tracer_snapshot = _imp1_reference_tracer_snapshot(tracers)
    guard_configuration = _imp1_reference_guard_configuration(transaction)
    ledger_hash = _imp1_reference_ledger_digest(temporal_ledger)

    def unchanged() -> bool:
        return (
            _imp1_reference_state_hash(state) == state_hash
            and array_content_sha256(
                _imp1_reference_coordinates(coordinates, state)
            )
            == coordinate_hash
            and _imp1_reference_digest(transaction.state)
            == _imp1_reference_digest(monitor)
            and _imp1_reference_digest(transaction.causal_state)
            == _imp1_reference_digest(causal)
            and _imp1_reference_digest(
                _imp1_reference_guard_configuration(transaction)
            )
            == _imp1_reference_digest(guard_configuration)
            and _imp1_reference_tracer_snapshot(tracers) == tracer_snapshot
            and _imp1_reference_ledger_digest(temporal_ledger) == ledger_hash
        )

    paths: list[_tdg6_shadow_runtime.TDG6ShadowPath] = []
    try:
        for level, boundaries in zip(
            ("outer", "medium", "fine"),
            (plan.one_full, plan.two_half, plan.four_quarter),
            strict=True,
        ):
            paths.append(
                _imp1_reference_shadow_path(
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
        assessment = assess_c1r1_family(family, limits=limits)
    except _tdg6_shadow_runtime.TDG6RefinementPathStop as error:
        raise C1R1RefinementPathStop(error) from error
    finally:
        if not unchanged():
            raise RuntimeError("C1R1 preparation changed a real accepted boundary")

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
        implementation_id=C1R1_IMPLEMENTATION_ID,
        mathematical_object=C1R1_MATHEMATICAL_OBJECT,
        reference_wire_evaluator_id=C1R1_REFERENCE_WIRE_EVALUATOR_ID,
        original_rhs=rhs,
        original_projector=projector,
    )
    prepared = object.__new__(TDG11C1R1PreparedRuntime)
    for name, value in values.items():
        object.__setattr__(prepared, name, value)
    object.__setattr__(prepared, "receipt_sha256", _receipt_digest(prepared))
    validate_c1r1_prepared(prepared)
    return prepared


def prepare_c1r1_initial_runtime(
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
) -> TDG11C1R1PreparedRuntime:
    """Prepare a fresh macro attempt; lattice validation is the first action."""

    plan = _plan(time, event_target, requested_cap)
    if not isinstance(temporal_ledger, TDG11IMP1Ledger):
        raise TypeError("temporal_ledger must be TDG11IMP1Ledger")
    if temporal_ledger.current_macro_step_temporal_retry_count != 0:
        raise ValueError("C1R1 fresh entrypoint cannot consume an active retry ledger")
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


def _rejection_record(prepared: TDG11C1R1PreparedRuntime) -> dict[str, object]:
    """IMP1 ledger-wire rejection; preparation hash is the C1R1 receipt."""

    ledger = prepared.original_ledger
    return TDG11IMP1TemporalRejection(
        event_type=IMP1_REJECTION_EVENT_TYPE,
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


class C1R1TemporalRetryRequired(RuntimeError):
    """One durably acknowledged nonadmission, not a scientific verdict."""

    physical_classification = False
    accepted_state_advanced = False

    def __init__(
        self, record: Mapping[str, object], updated_ledger: TDG11IMP1Ledger
    ) -> None:
        self.record_bytes = canonical_json_bytes(dict(record))
        self.updated_ledger = updated_ledger
        self.retry_step_size = float.fromhex(str(record["next_cap_hex"]))
        super().__init__("C1R1 temporal admission requests a stage-safe half-cap retry")

    @property
    def evidence(self) -> dict[str, object]:
        return json.loads(self.record_bytes)


class C1R1TemporalRetryExhausted(RuntimeError):
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
            raise ValueError("unknown C1R1 retry exhaustion reason")
        self.reason = reason
        self.record_bytes = canonical_json_bytes(dict(record))
        self.updated_ledger = updated_ledger
        super().__init__(f"C1R1 temporal retry exhausted: {reason}")

    @property
    def evidence(self) -> dict[str, object]:
        return json.loads(self.record_bytes)


def require_c1r1_admission(
    prepared: TDG11C1R1PreparedRuntime,
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
    """Require all channels; emit a retry only after the durable sink succeeds."""

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
        raise TypeError("C1R1 nonadmission requires a durable rejection sink")
    record = _rejection_record(prepared)
    expected_bytes = canonical_json_bytes(record)
    updated = append_imp1_rejection(temporal_ledger, record)
    _imp1_reference_ledger_digest(updated)
    durable_rejection_sink(record)
    if canonical_json_bytes(record) != expected_bytes:
        raise RuntimeError("C1R1 rejection sink altered the acknowledged evidence")
    _revalidate_boundary(prepared, **current)
    if (
        updated.current_macro_step_temporal_retry_count
        > TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
    ):
        raise C1R1TemporalRetryExhausted(
            reason="maximum_temporal_retries", record=record, updated_ledger=updated
        )
    if prepared.plan.macro_width / 2.0 < float(TDG6_MINIMUM_MACRO_STEP):
        raise C1R1TemporalRetryExhausted(
            reason="minimum_macro_step", record=record, updated_ledger=updated
        )
    raise C1R1TemporalRetryRequired(record, updated)


@dataclass(frozen=True, slots=True)
class C1R1RetrySuccessor:
    """Immutable preceding nonadmission and its exact next TDG7-safe plan."""

    predecessor: TDG11C1R1PreparedRuntime
    plan: TDG7Binary64SubdivisionPlan
    updated_ledger: TDG11IMP1Ledger
    rejection_bytes: bytes

    def __post_init__(self) -> None:
        validate_c1r1_prepared(self.predecessor)
        validate_tdg7_binary64_subdivision_plan(self.plan)
        if self.predecessor.assessment.admission_passed is not False:
            raise ValueError("C1R1 cannot turn an admitted family into a retry")
        expected_record = _rejection_record(self.predecessor)
        if self.rejection_bytes != canonical_json_bytes(expected_record):
            raise ValueError("C1R1 retry evidence differs from its predecessor")
        expected_ledger = append_imp1_rejection(
            self.predecessor.original_ledger, expected_record
        )
        if _imp1_reference_ledger_digest(
            self.updated_ledger
        ) != _imp1_reference_ledger_digest(expected_ledger):
            raise ValueError(
                "C1R1 retry ledger omits or changes the preceding rejection"
            )
        if (
            self.updated_ledger.current_macro_step_temporal_retry_count
            > TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
        ):
            raise ValueError("C1R1 retry exceeds the frozen retry count")
        expected_plan = _plan(
            self.predecessor.plan.current,
            self.predecessor.plan.event_target,
            self.predecessor.plan.macro_width / 2.0,
        )
        if (
            self.plan != expected_plan
            or self.plan.macro_width >= self.predecessor.plan.macro_width
        ):
            raise ValueError("C1R1 retry does not use the certified halved cap")


def replan_c1r1_retry(
    prepared: TDG11C1R1PreparedRuntime, retry: C1R1TemporalRetryRequired
) -> C1R1RetrySuccessor:
    validate_c1r1_prepared(prepared)
    if not isinstance(retry, C1R1TemporalRetryRequired):
        raise TypeError("retry must be C1R1TemporalRetryRequired")
    if retry.record_bytes != canonical_json_bytes(_rejection_record(prepared)):
        raise ValueError("C1R1 retry object differs from the prepared predecessor")
    try:
        plan = _plan(
            prepared.plan.current,
            prepared.plan.event_target,
            prepared.plan.macro_width / 2.0,
        )
    except C1R1CoordinateLatticeStop as error:
        raise C1R1TemporalRetryExhausted(
            reason="coordinate_lattice",
            record=retry.evidence,
            updated_ledger=retry.updated_ledger,
        ) from error
    return C1R1RetrySuccessor(prepared, plan, retry.updated_ledger, retry.record_bytes)


def prepare_c1r1_retry_runtime(
    successor: C1R1RetrySuccessor,
    *,
    state: EvolutionState,
    rhs: Callable[[float, EvolutionState], EvolutionRHS],
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    coordinates: object,
) -> TDG11C1R1PreparedRuntime:
    """No caller replacement time, cap, ledger, method or indices are accepted."""

    if not isinstance(successor, C1R1RetrySuccessor):
        raise TypeError("successor must be C1R1RetrySuccessor")
    successor.__post_init__()
    prior = successor.predecessor
    if rhs is not prior.original_rhs or projector is not prior.original_projector:
        raise ValueError(
            "C1R1 live retry cannot replace its source or projector callback"
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
class TDG11C1R1CommittedRuntime:
    """In-memory fine endpoint; intentionally not an endpoint codec."""

    state: EvolutionState
    time: float
    step_index: int
    transaction_serial: int
    temporal_ledger: TDG11IMP1Ledger
    method: str
    assessment_sha256: str
    preparation_sha256: str
    implementation_id: str = C1R1_IMPLEMENTATION_ID
    committed_proposal_count: int = 4
    corrected_endpoint_adopted: bool = False
    outer_path_adopted: bool = False
    medium_path_adopted: bool = False
    campaign_execution_authorized: bool = False

    def __post_init__(self) -> None:
        if self.method not in _imp1_reference_stage_counts:
            raise ValueError("unknown C1R1 committed method")
        if (
            type(self.committed_proposal_count) is not int
            or self.committed_proposal_count != 4
            or self.corrected_endpoint_adopted is not False
            or self.outer_path_adopted is not False
            or self.medium_path_adopted is not False
            or self.campaign_execution_authorized is not False
            or self.implementation_id != C1R1_IMPLEMENTATION_ID
        ):
            raise ValueError("C1R1 committed record crosses its fine-only kernel scope")
        if not isinstance(self.temporal_ledger, TDG11IMP1Ledger):
            raise TypeError("C1R1 commit requires the IMP1 ledger wire")
        if (
            self.temporal_ledger.method != self.method
            or not _imp1_reference_same_bits(
                self.temporal_ledger.last_accepted_time, self.time
            )
            or self.temporal_ledger.current_state_sha256
            != _imp1_reference_state_hash(self.state)
            or self.temporal_ledger.current_step_index != self.step_index
            or self.temporal_ledger.current_transaction_serial
            != self.transaction_serial
        ):
            raise ValueError("C1R1 committed endpoint and ledger disagree")
        for digest in (self.assessment_sha256, self.preparation_sha256):
            if (
                type(digest) is not str
                or len(digest) != 64
                or any(char not in "0123456789abcdef" for char in digest)
            ):
                raise ValueError("invalid C1R1 commit receipt digest")


def commit_c1r1_runtime(
    prepared: TDG11C1R1PreparedRuntime,
    *,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    temporal_ledger: TDG11IMP1Ledger,
    current_time: Real,
    current_state: EvolutionState,
    current_step_index: int,
    current_transaction_serial: int,
    coordinates: object,
) -> TDG11C1R1CommittedRuntime:
    """Adopt only the fine binary64 endpoint after complete revalidation."""

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
        raise ValueError("C1R1 cannot commit an unadmitted family")
    accepted = prepared.fine.final_accepted
    updated = accept_imp1_step(
        temporal_ledger,
        final_time=accepted.time,
        state_sha256=_imp1_reference_state_hash(accepted.state),
        step_index=accepted.step_index,
        transaction_serial=accepted.transaction_serial,
        debit_vector=prepared.assessment.public_fine_debits,
        assessment_sha256=prepared.assessment.assessment_sha256,
    )
    committed = TDG11C1R1CommittedRuntime(
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
        _imp1_reference_clone_tracers(tracers),
    )
    guard_configuration = _imp1_reference_guard_configuration(transaction)
    try:
        transaction.state = prepared.fine.final_monitor_state
        transaction.causal_state = prepared.fine.final_causal_state
        tracers.commit_advance(
            prepared.fine.final_tracer_positions,
            prepared.fine.final_tracer_proper_times,
        )
        validate_c1r1_prepared(prepared)
        if (
            _imp1_reference_tracer_snapshot(tracers)
            != prepared.fine.final_tracer_snapshot
            or transaction.state != prepared.fine.final_monitor_state
            or transaction.causal_state != prepared.fine.final_causal_state
            or _imp1_reference_digest(
                _imp1_reference_guard_configuration(transaction)
            )
            != _imp1_reference_digest(guard_configuration)
            or _receipt_digest(prepared) != prepared.receipt_sha256
            or _imp1_reference_state_hash(current_state)
            != prepared.initial_state_sha256
            or _imp1_reference_state_hash(accepted.state)
            != updated.current_state_sha256
            or array_content_sha256(
                _imp1_reference_coordinates(coordinates, current_state)
            )
            != prepared.coordinates_sha256
        ):
            raise RuntimeError(
                "C1R1 adopted boundary differs from the guarded fine receipt"
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
            _imp1_reference_restore_tracers(tracers, prior_tracers)
            if (
                _imp1_reference_tracer_snapshot(tracers)
                != prepared.original_tracer_snapshot
                or transaction.state != monitor
                or transaction.causal_state != causal
            ):
                raise RuntimeError("restored C1R1 boundary differs")
        except BaseException as error:
            raise RuntimeError("C1R1 complete tracer/ledger rollback failed") from error
        raise
    return committed


__all__ = [
    "ARTIFACT_ID",
    "CAMPAIGN_EXECUTION_AUTHORIZED",
    "C1R1CoordinateLatticeStop",
    "C1R1RefinementPathStop",
    "C1R1RetrySuccessor",
    "C1R1RuntimeResourceStop",
    "C1R1TemporalRetryExhausted",
    "C1R1TemporalRetryRequired",
    "FRESH_VS_RETRY_DURABLE_PROTOCOL_IMPLEMENTED",
    "REFERENCE_WIRE_ARTIFACT_ID",
    "TDG11C1R1CommittedRuntime",
    "TDG11C1R1PreparedRuntime",
    "commit_c1r1_runtime",
    "prepare_c1r1_initial_runtime",
    "prepare_c1r1_retry_runtime",
    "replan_c1r1_retry",
    "require_c1r1_admission",
    "validate_c1r1_prepared",
]
