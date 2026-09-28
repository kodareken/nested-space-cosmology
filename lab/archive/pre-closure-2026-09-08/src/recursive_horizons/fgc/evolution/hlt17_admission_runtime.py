"""Fixed C1R1 admission runtime for prospective in-memory HLT17.

This is not a backend registry and not a production authority. The HLT17
in-memory path uses one named C1R1 implementation. IMP1 remains the sealed
mathematical and ledger-wire reference. Prepared objects stay C1R1-named;
an IMP1 prepared object cannot be adopted as C1R1 execution, and IMP1
``_prepare`` is not run and then replaced.
"""

from __future__ import annotations

from typing import Callable, Final, Mapping
from numbers import Real

from . import tdg11_c1r1_runtime as c1r1
from . import tdg11_imp1_runtime as imp1
from .numerical_engine import EvolutionRHS, EvolutionState
from .proto5_runtime import GR0RuntimeStageTransaction
from .tdg11_c1r1_enclosure import (
    C1R1_IMPLEMENTATION_ID,
    C1R1_MATHEMATICAL_OBJECT,
    C1R1_REFERENCE_WIRE_EVALUATOR_ID,
)
from .tdg11_imp1_enclosure import IMP1EnclosureLimits
from .tdg11_imp1_ledger import TDG11IMP1Ledger
from .tdg7_binary64_subdivision_lattice import TDG7Binary64SubdivisionPlan


HLT17_ADMISSION_RUNTIME_ID: Final[str] = "hlt17_c1r1_admission_runtime_v1"
HLT17_ADMISSION_ARTIFACT_ID: Final[str] = "FGC-1-HLT17-C1R1-ADMISSION-v1"
HLT17_IMPLEMENTATION_ID: Final[str] = C1R1_IMPLEMENTATION_ID
HLT17_MATHEMATICAL_OBJECT: Final[str] = C1R1_MATHEMATICAL_OBJECT
HLT17_REFERENCE_WIRE_EVALUATOR_ID: Final[str] = C1R1_REFERENCE_WIRE_EVALUATOR_ID
HLT17_REFERENCE_WIRE_ARTIFACT_ID: Final[str] = c1r1.REFERENCE_WIRE_ARTIFACT_ID
HLT17_C1R1_ARTIFACT_ID: Final[str] = c1r1.ARTIFACT_ID
CAMPAIGN_EXECUTION_AUTHORIZED: Final[bool] = False

IMPLEMENTATION_IDENTITY_KEYS: Final[tuple[str, ...]] = (
    "implementation_id",
    "mathematical_object",
    "reference_wire_artifact_id",
    "reference_wire_evaluator_id",
    "admission_runtime_id",
)

# C1R1 private instruments used by HLT17 replay and reduced-plan prepare.
# None of these is a new public historical API. IMP1._prepare is not pinned.
PINNED_C1R1_PREPARE: Final = c1r1._prepare
PINNED_C1R1_PLAN: Final = c1r1._plan
PINNED_C1R1_REVALIDATE_BOUNDARY: Final = c1r1._revalidate_boundary
PINNED_C1R1_REJECTION_RECORD: Final = c1r1._rejection_record
PINNED_C1R1_RECEIPT_DIGEST: Final = c1r1._receipt_digest
PINNED_C1R1_PREPARED_MAPPING: Final = c1r1._prepared_mapping
PINNED_C1R1_VALIDATE_PREPARED: Final = c1r1.validate_c1r1_prepared

PINNED_C1R1_PRIVATE_DEPENDENCIES: Final[tuple[tuple[str, object], ...]] = (
    ("tdg11_c1r1_runtime._prepare", PINNED_C1R1_PREPARE),
    ("tdg11_c1r1_runtime._plan", PINNED_C1R1_PLAN),
    ("tdg11_c1r1_runtime._revalidate_boundary", PINNED_C1R1_REVALIDATE_BOUNDARY),
    ("tdg11_c1r1_runtime._rejection_record", PINNED_C1R1_REJECTION_RECORD),
    ("tdg11_c1r1_runtime._receipt_digest", PINNED_C1R1_RECEIPT_DIGEST),
    ("tdg11_c1r1_runtime._prepared_mapping", PINNED_C1R1_PREPARED_MAPPING),
)

C1R1CoordinateLatticeStop = c1r1.C1R1CoordinateLatticeStop
C1R1RefinementPathStop = c1r1.C1R1RefinementPathStop
C1R1RuntimeResourceStop = c1r1.C1R1RuntimeResourceStop
C1R1TemporalRetryRequired = c1r1.C1R1TemporalRetryRequired
C1R1TemporalRetryExhausted = c1r1.C1R1TemporalRetryExhausted
C1R1RetrySuccessor = c1r1.C1R1RetrySuccessor
TDG11C1R1PreparedRuntime = c1r1.TDG11C1R1PreparedRuntime
TDG11C1R1CommittedRuntime = c1r1.TDG11C1R1CommittedRuntime


def hlt17_implementation_identity() -> dict[str, str]:
    """Closed identity of the fixed C1R1 admission implementation."""

    mapping = {
        "implementation_id": HLT17_IMPLEMENTATION_ID,
        "mathematical_object": HLT17_MATHEMATICAL_OBJECT,
        "reference_wire_artifact_id": HLT17_REFERENCE_WIRE_ARTIFACT_ID,
        "reference_wire_evaluator_id": HLT17_REFERENCE_WIRE_EVALUATOR_ID,
        "admission_runtime_id": HLT17_ADMISSION_RUNTIME_ID,
    }
    if tuple(mapping) != IMPLEMENTATION_IDENTITY_KEYS:
        raise RuntimeError("HLT17 implementation identity keys differ")
    if mapping["implementation_id"] != C1R1_IMPLEMENTATION_ID:
        raise RuntimeError("HLT17 implementation_id is not the named C1R1 object")
    if mapping["reference_wire_artifact_id"] != imp1.ARTIFACT_ID:
        raise RuntimeError("HLT17 reference-wire artifact is not the sealed IMP1 wire")
    return mapping


def require_hlt17_implementation_identity(value: object, *, label: str) -> dict[str, str]:
    """Require the fixed C1R1 identity; refuse any other backend marker."""

    expected = hlt17_implementation_identity()
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a mapping")
    item = dict(value)
    if set(item) != set(IMPLEMENTATION_IDENTITY_KEYS):
        raise ValueError(f"{label} implementation identity is incomplete or broadened")
    for name in IMPLEMENTATION_IDENTITY_KEYS:
        found = item[name]
        if type(found) is not str or not found:
            raise TypeError(f"{label} {name} must be nonempty text")
        if found != expected[name]:
            raise ValueError(f"{label} {name} is not the fixed C1R1 admission identity")
    return expected


def _reject_imp1_prepared(prepared: object) -> None:
    if isinstance(prepared, imp1.TDG11IMP1PreparedRuntime):
        raise TypeError(
            "HLT17 cannot adopt an IMP1 prepared runtime as C1R1 execution"
        )


def _reject_imp1_committed(committed: object) -> None:
    if type(committed) is imp1.TDG11IMP1CommittedRuntime:
        raise TypeError(
            "HLT17 cannot adopt an IMP1 committed runtime as C1R1 execution"
        )


def _reject_imp1_successor(successor: object) -> None:
    if isinstance(successor, imp1.IMP1RetrySuccessor):
        raise TypeError(
            "HLT17 cannot adopt an IMP1 retry successor as C1R1 execution"
        )


def require_hlt17_prepared(prepared: object) -> TDG11C1R1PreparedRuntime:
    """Require a live C1R1 prepared family. IMP1 prepared objects are foreign."""

    _reject_imp1_prepared(prepared)
    if not isinstance(prepared, TDG11C1R1PreparedRuntime):
        raise TypeError("prepared must be TDG11C1R1PreparedRuntime")
    if prepared.implementation_id != HLT17_IMPLEMENTATION_ID:
        raise ValueError("prepared implementation_id is not the fixed C1R1 object")
    if prepared.mathematical_object != HLT17_MATHEMATICAL_OBJECT:
        raise ValueError("prepared mathematical object is not the C1 accumulation object")
    if prepared.reference_wire_evaluator_id != HLT17_REFERENCE_WIRE_EVALUATOR_ID:
        raise ValueError("prepared reference-wire evaluator identity differs")
    PINNED_C1R1_VALIDATE_PREPARED(prepared)
    return prepared


def require_hlt17_committed(committed: object) -> TDG11C1R1CommittedRuntime:
    """Require a C1R1 fine commit. IMP1 committed objects are foreign."""

    _reject_imp1_committed(committed)
    if type(committed) is not TDG11C1R1CommittedRuntime:
        raise TypeError("committed must be TDG11C1R1CommittedRuntime")
    committed.__post_init__()
    if committed.implementation_id != HLT17_IMPLEMENTATION_ID:
        raise ValueError("committed implementation_id is not the fixed C1R1 object")
    if committed.corrected_endpoint_adopted is not False:
        raise ValueError("cannot scrape an unqualified corrected endpoint")
    return committed


def require_hlt17_retry_successor(successor: object) -> C1R1RetrySuccessor:
    """Require a C1R1 half-cap successor. IMP1 successors are foreign."""

    _reject_imp1_successor(successor)
    if not isinstance(successor, C1R1RetrySuccessor):
        raise TypeError("successor must be C1R1RetrySuccessor")
    successor.__post_init__()
    require_hlt17_prepared(successor.predecessor)
    return successor


def prepare_hlt17_admission_fresh(
    *,
    method: str,
    time: Real,
    event_target: Real,
    requested_cap: object,
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
    """Public C1R1 fresh entrypoint. Does not call sealed IMP1 ``_prepare``."""

    prepared = c1r1.prepare_c1r1_initial_runtime(
        method=method,
        time=time,
        event_target=event_target,
        requested_cap=requested_cap,  # type: ignore[arg-type]
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
    return require_hlt17_prepared(prepared)


def prepare_hlt17_admission_with_plan(
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
    """Reduced-plan C1R1 prepare. Does not call sealed IMP1 ``_prepare``."""

    prepared = PINNED_C1R1_PREPARE(
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
    return require_hlt17_prepared(prepared)


def prepare_hlt17_admission_retry(
    successor: object,
    *,
    state: EvolutionState,
    rhs: Callable[[float, EvolutionState], EvolutionRHS],
    projector: Callable[[float, EvolutionState], EvolutionState] | None,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    coordinates: object,
) -> TDG11C1R1PreparedRuntime:
    checked = require_hlt17_retry_successor(successor)
    prepared = c1r1.prepare_c1r1_retry_runtime(
        checked,
        state=state,
        rhs=rhs,
        projector=projector,
        transaction=transaction,
        tracers=tracers,
        coordinates=coordinates,
    )
    return require_hlt17_prepared(prepared)


def require_hlt17_admission_family(
    prepared: object,
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
    bound = require_hlt17_prepared(prepared)
    c1r1.require_c1r1_admission(
        bound,
        transaction=transaction,
        tracers=tracers,
        temporal_ledger=temporal_ledger,
        current_time=current_time,
        current_state=current_state,
        current_step_index=current_step_index,
        current_transaction_serial=current_transaction_serial,
        coordinates=coordinates,
        durable_rejection_sink=durable_rejection_sink,
    )


def commit_hlt17_admission(
    prepared: object,
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
    bound = require_hlt17_prepared(prepared)
    committed = c1r1.commit_c1r1_runtime(
        bound,
        transaction=transaction,
        tracers=tracers,
        temporal_ledger=temporal_ledger,
        current_time=current_time,
        current_state=current_state,
        current_step_index=current_step_index,
        current_transaction_serial=current_transaction_serial,
        coordinates=coordinates,
    )
    return require_hlt17_committed(committed)


def revalidate_hlt17_admission_boundary(
    prepared: object,
    *,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    temporal_ledger: TDG11IMP1Ledger,
    current_time: Real,
    current_state: EvolutionState,
    current_step_index: int,
    current_transaction_serial: int,
    coordinates: object,
) -> TDG11C1R1PreparedRuntime:
    bound = require_hlt17_prepared(prepared)
    PINNED_C1R1_REVALIDATE_BOUNDARY(
        bound,
        transaction=transaction,
        tracers=tracers,
        temporal_ledger=temporal_ledger,
        current_time=current_time,
        current_state=current_state,
        current_step_index=current_step_index,
        current_transaction_serial=current_transaction_serial,
        coordinates=coordinates,
    )
    return bound


def hlt17_admission_plan(
    time: Real, event_target: Real, cap: Real
) -> TDG7Binary64SubdivisionPlan:
    return PINNED_C1R1_PLAN(time, event_target, cap)


def hlt17_rejection_record(prepared: object) -> dict[str, object]:
    return dict(PINNED_C1R1_REJECTION_RECORD(require_hlt17_prepared(prepared)))


def hlt17_receipt_digest(prepared: object) -> str:
    bound = require_hlt17_prepared(prepared)
    digest = PINNED_C1R1_RECEIPT_DIGEST(bound)
    if digest != bound.receipt_sha256:
        raise ValueError("C1R1 receipt digest is not redundant")
    return digest


def make_hlt17_retry_successor(
    reconstructed: object,
    plan: TDG7Binary64SubdivisionPlan,
    updated_ledger: TDG11IMP1Ledger,
    rejection_bytes: bytes,
) -> C1R1RetrySuccessor:
    predecessor = require_hlt17_prepared(reconstructed)
    successor = C1R1RetrySuccessor(
        predecessor, plan, updated_ledger, rejection_bytes
    )
    return require_hlt17_retry_successor(successor)


FUTURE_STORE_IMPLEMENTATION_SEAM: Final[str] = (
    "protocol_v19 and the HLT17 campaign store must persist this cursor/"
    "codec implementation identity and refuse IMP1-cursor-v1, member-state-v1, "
    "foreign IMP1 prepared objects, and silent C1R1 backend changes; those "
    "files are owned by a separate store worker"
)


__all__ = [
    "CAMPAIGN_EXECUTION_AUTHORIZED",
    "C1R1CoordinateLatticeStop",
    "C1R1RefinementPathStop",
    "C1R1RetrySuccessor",
    "C1R1RuntimeResourceStop",
    "C1R1TemporalRetryExhausted",
    "C1R1TemporalRetryRequired",
    "FUTURE_STORE_IMPLEMENTATION_SEAM",
    "HLT17_ADMISSION_ARTIFACT_ID",
    "HLT17_ADMISSION_RUNTIME_ID",
    "HLT17_C1R1_ARTIFACT_ID",
    "HLT17_IMPLEMENTATION_ID",
    "HLT17_MATHEMATICAL_OBJECT",
    "HLT17_REFERENCE_WIRE_ARTIFACT_ID",
    "HLT17_REFERENCE_WIRE_EVALUATOR_ID",
    "IMPLEMENTATION_IDENTITY_KEYS",
    "PINNED_C1R1_PLAN",
    "PINNED_C1R1_PREPARE",
    "PINNED_C1R1_PREPARED_MAPPING",
    "PINNED_C1R1_PRIVATE_DEPENDENCIES",
    "PINNED_C1R1_RECEIPT_DIGEST",
    "PINNED_C1R1_REJECTION_RECORD",
    "PINNED_C1R1_REVALIDATE_BOUNDARY",
    "PINNED_C1R1_VALIDATE_PREPARED",
    "TDG11C1R1CommittedRuntime",
    "TDG11C1R1PreparedRuntime",
    "commit_hlt17_admission",
    "hlt17_admission_plan",
    "hlt17_implementation_identity",
    "hlt17_receipt_digest",
    "hlt17_rejection_record",
    "make_hlt17_retry_successor",
    "prepare_hlt17_admission_fresh",
    "prepare_hlt17_admission_retry",
    "prepare_hlt17_admission_with_plan",
    "require_hlt17_admission_family",
    "require_hlt17_committed",
    "require_hlt17_implementation_identity",
    "require_hlt17_prepared",
    "require_hlt17_retry_successor",
    "revalidate_hlt17_admission_boundary",
]
