"""Narrow HLT17 adapter for C1R1 replay, reduced-plan prepare, and mixed retries.

Public C1R1 fresh is used only at zero active retry depth with no overlay.
Replay of a persisted rejection and any source/CFL reduced plan use the
explicitly pinned private C1R1 ``_prepare`` instrument. Sealed IMP1
``_prepare`` is not run and then replaced. IMP1 remains the ledger and
rejection-wire reference. This adapter does not read a store or authenticate
an origin, physical source, or campaign writer.
"""

from __future__ import annotations

from dataclasses import dataclass, fields, field, is_dataclass
from hashlib import sha256
import math
from numbers import Real
import types
from typing import Callable, Final, Mapping

from recursive_horizons.evidence_io import CanonicalJSONError, canonical_json_bytes

from . import tdg6_temporal_admission_runtime as tdg6
from . import tdg11_imp1_runtime as imp1
from .hlt17_admission_runtime import (
    C1R1CoordinateLatticeStop,
    C1R1RefinementPathStop,
    C1R1RetrySuccessor,
    C1R1RuntimeResourceStop,
    C1R1TemporalRetryExhausted,
    C1R1TemporalRetryRequired,
    FUTURE_STORE_IMPLEMENTATION_SEAM,
    HLT17_IMPLEMENTATION_ID,
    PINNED_C1R1_PLAN,
    PINNED_C1R1_PREPARE,
    PINNED_C1R1_PREPARED_MAPPING,
    PINNED_C1R1_PRIVATE_DEPENDENCIES,
    PINNED_C1R1_RECEIPT_DIGEST,
    PINNED_C1R1_REJECTION_RECORD,
    PINNED_C1R1_REVALIDATE_BOUNDARY,
    TDG11C1R1CommittedRuntime,
    TDG11C1R1PreparedRuntime,
    commit_hlt17_admission,
    hlt17_implementation_identity,
    hlt17_receipt_digest,
    hlt17_rejection_record,
    make_hlt17_retry_successor,
    prepare_hlt17_admission_fresh,
    prepare_hlt17_admission_retry,
    require_hlt17_admission_family,
    require_hlt17_committed,
    require_hlt17_prepared,
)
from .hlt17_imp1_cursor import (
    CFL_OWNER,
    FRESH_READY,
    FUTURE_AUTHENTICATED_ORIGIN_SEAM,
    FUTURE_CAMPAIGN_TRANSACTION_SEAM,
    FUTURE_CHECKPOINT_JOURNAL_SEAM,
    FUTURE_CONTENT_ADDRESSED_STORE_SEAM,
    FUTURE_MEMBER_CODEC_SEAM,
    FUTURE_SOURCE_MANIFEST_SEAM,
    HLT17IMP1Cursor,
    HLT17IMP1Error,
    HLT17InvalidPremise,
    HLT17OwnerRetryExhausted,
    HLT17RetryOverlay,
    IDENTITY_SCOPE_LIVE,
    SOURCE_OWNER,
    hlt17_member_spec,
    close_fine_adoption,
    close_source_cfl_overlay,
    close_temporal_rejection,
    encode_tdg7_plan,
    imp1_enclosure_limits_for_member,
    ledger_sha256,
    plan_imp1_lattice,
    revalidate_hlt17_cursor,
    seed_hlt17_cursor,
)
from .numerical_engine import EvolutionRHS, EvolutionState, array_content_sha256
from .proto5_runtime import CFLRetryRequired, GR0RuntimeStageTransaction
from .tdg11_imp1_enclosure import IMP1EnclosureLimits
from .tdg11_imp1_ledger import (
    TDG11IMP1Ledger,
    TDG11IMP1TemporalRejection,
    append_imp1_rejection,
    imp1_checkpoint_extension,
)
from .tdg7_binary64_subdivision_lattice import (
    TDG7Binary64SubdivisionPlan,
    validate_tdg7_binary64_subdivision_plan,
)


# Sealed IMP1 ledger/wire helpers remain reference-wire compatibility.
# They are not a claim that IMP1 assessed the family. Actual prepare is C1R1.
PINNED_IMP1_CHECK_CURRENT_BOUNDARY: Final = imp1._check_current_boundary
PINNED_IMP1_LEDGER_DIGEST: Final = imp1._ledger_digest
PINNED_IMP1_STATE_HASH: Final = imp1._state_hash
PINNED_IMP1_DIGEST: Final = imp1._digest
PINNED_IMP1_WIRE: Final = imp1._wire
PINNED_IMP1_GUARD_CONFIGURATION: Final = imp1._guard_configuration
PINNED_IMP1_COORDINATES: Final = imp1._coordinates
PINNED_TDG6_TRACER_SNAPSHOT: Final = tdg6._tracer_snapshot
PINNED_TDG6_CLONE_TRACERS: Final = tdg6._clone_tracers
PINNED_TDG6_RESTORE_TRACERS: Final = tdg6._restore_tracers

PINNED_IMP1_PRIVATE_DEPENDENCIES: Final[tuple[tuple[str, object], ...]] = (
    ("tdg11_imp1_runtime._check_current_boundary", PINNED_IMP1_CHECK_CURRENT_BOUNDARY),
    ("tdg11_imp1_runtime._ledger_digest", PINNED_IMP1_LEDGER_DIGEST),
    ("tdg11_imp1_runtime._state_hash", PINNED_IMP1_STATE_HASH),
    ("tdg11_imp1_runtime._digest", PINNED_IMP1_DIGEST),
    ("tdg11_imp1_runtime._wire", PINNED_IMP1_WIRE),
    ("tdg11_imp1_runtime._guard_configuration", PINNED_IMP1_GUARD_CONFIGURATION),
    ("tdg11_imp1_runtime._coordinates", PINNED_IMP1_COORDINATES),
    ("tdg6_temporal_admission_runtime._tracer_snapshot", PINNED_TDG6_TRACER_SNAPSHOT),
    ("tdg6_temporal_admission_runtime._clone_tracers", PINNED_TDG6_CLONE_TRACERS),
    ("tdg6_temporal_admission_runtime._restore_tracers", PINNED_TDG6_RESTORE_TRACERS),
)
PINNED_HLT17_PRIVATE_DEPENDENCIES: Final[tuple[tuple[str, object], ...]] = (
    PINNED_C1R1_PRIVATE_DEPENDENCIES + PINNED_IMP1_PRIVATE_DEPENDENCIES
)

FUTURE_INTEGRATION_SEAMS: Final[dict[str, str]] = {
    "member_codec": FUTURE_MEMBER_CODEC_SEAM,
    "content_addressed_store": FUTURE_CONTENT_ADDRESSED_STORE_SEAM,
    "checkpoint_journal": FUTURE_CHECKPOINT_JOURNAL_SEAM,
    "authenticated_origin": FUTURE_AUTHENTICATED_ORIGIN_SEAM,
    "campaign_transaction": FUTURE_CAMPAIGN_TRANSACTION_SEAM,
    "source_manifest": FUTURE_SOURCE_MANIFEST_SEAM,
    "store_implementation_identity": FUTURE_STORE_IMPLEMENTATION_SEAM,
}

DISPOSITIONS: Final[frozenset[str]] = frozenset(
    {
        "prepared_fresh",
        "replay_authenticated",
        "prepared_successor",
        "prepared_overlay",
        "temporal_retry_required",
        "temporal_retry_exhausted",
        "source_retry_required",
        "source_retry_exhausted",
        "cfl_retry_required",
        "cfl_retry_exhausted",
        "accepted_fine",
        "invalid_premise",
    }
)


MAX_BINDING_DEPTH: Final[int] = 16
EXPLICIT_BINDING_SEAM: Final[str] = "hlt17_synthetic_callable_binding"


class HLT17UnsupportedCallableBinding(HLT17IMP1Error):
    """The live callable is opaque, mutable, or otherwise not inspectable."""


def _fail_premise(message: str) -> None:
    raise HLT17InvalidPremise(message)


def _fail_binding(message: str) -> None:
    raise HLT17UnsupportedCallableBinding(message)


def _code_material(code: types.CodeType, *, depth: int, seen: set[int]) -> dict[str, object]:
    return {
        "co_code_sha256": sha256(code.co_code).hexdigest(),
        "co_consts": [
            _binding_atom(item, depth=depth + 1, seen=seen) for item in code.co_consts
        ],
        "co_names": list(code.co_names),
        "co_varnames": list(code.co_varnames),
        "co_freevars": list(code.co_freevars),
        "co_cellvars": list(code.co_cellvars),
        "co_argcount": code.co_argcount,
        "co_posonlyargcount": code.co_posonlyargcount,
        "co_kwonlyargcount": code.co_kwonlyargcount,
        "co_nlocals": code.co_nlocals,
        "co_flags": code.co_flags,
    }


def _function_material(
    func: types.FunctionType, *, depth: int, seen: set[int]
) -> dict[str, object]:
    marker = id(func)
    if marker in seen:
        _fail_binding("callable binding encountered a cyclic function")
    seen = set(seen)
    seen.add(marker)
    defaults = func.__defaults__
    kwdefaults = func.__kwdefaults__
    closure = func.__closure__
    freevars = func.__code__.co_freevars
    cells: list[object] = []
    if closure is not None:
        if len(closure) != len(freevars):
            _fail_binding("python function free variables differ from closure cells")
        for name, cell in zip(freevars, closure, strict=True):
            try:
                contents = cell.cell_contents
            except ValueError:
                _fail_binding(f"python function closure cell {name!r} is empty")
            cells.append(
                {
                    "name": name,
                    "value": _binding_atom(contents, depth=depth + 1, seen=seen),
                }
            )
    return {
        "kind": "python_function",
        "module": getattr(func, "__module__", ""),
        "qualname": getattr(func, "__qualname__", getattr(func, "__name__", "")),
        "code": _code_material(func.__code__, depth=depth, seen=seen),
        "defaults": (
            None
            if defaults is None
            else [_binding_atom(item, depth=depth + 1, seen=seen) for item in defaults]
        ),
        "kwdefaults": (
            None
            if kwdefaults is None
            else {
                str(name): _binding_atom(item, depth=depth + 1, seen=seen)
                for name, item in kwdefaults.items()
            }
        ),
        "closure": cells,
        "physical_source_authenticated": False,
    }


def _binding_atom(value: object, *, depth: int, seen: set[int]) -> object:
    if depth > MAX_BINDING_DEPTH:
        _fail_binding("callable binding exceeded its inspection depth")
    if value is None or type(value) is bool or type(value) is str:
        return value
    if type(value) is int:
        return value
    if type(value) is float:
        if not math.isfinite(value):
            _fail_binding("callable binding contains a nonfinite float")
        return {"binary64_hex": value.hex()}
    if type(value) is bytes:
        return {"bytes_hex": value.hex()}
    if value is Ellipsis:
        return {"ellipsis": True}
    if type(value) is tuple:
        return {
            "tuple": [
                _binding_atom(item, depth=depth + 1, seen=seen) for item in value
            ]
        }
    if type(value) is frozenset:
        encoded = [
            _binding_atom(item, depth=depth + 1, seen=seen) for item in value
        ]
        return {"frozenset": sorted(encoded, key=lambda item: repr(item))}
    if type(value) is types.CodeType:
        return {"code": _code_material(value, depth=depth, seen=seen)}
    if type(value) is types.FunctionType:
        return _function_material(value, depth=depth, seen=seen)
    if isinstance(value, types.MethodType):
        return {
            "kind": "bound_method",
            "func": _function_material(value.__func__, depth=depth, seen=seen),
            "self": _callable_object_material(value.__self__, depth=depth, seen=seen),
        }
    if callable(value):
        return _callable_object_material(value, depth=depth, seen=seen)
    _fail_binding(
        f"unsupported callable-binding value {type(value).__module__}."
        f"{type(value).__qualname__}"
    )
    raise AssertionError("unreachable")


def _explicit_binding_material(value: object) -> dict[str, object]:
    seam = getattr(value, EXPLICIT_BINDING_SEAM, None)
    if not callable(seam):
        _fail_binding("callable object omitted the explicit synthetic binding seam")
    try:
        material = seam()
    except HLT17UnsupportedCallableBinding:
        raise
    except Exception as error:
        raise HLT17UnsupportedCallableBinding(
            "explicit synthetic callable binding failed"
        ) from error
    if type(material) is not dict:
        _fail_binding("explicit synthetic callable binding must return a dict")
    if material.get("physical_source_authenticated") is True:
        _fail_binding("synthetic callable binding cannot claim physical authentication")
    try:
        canonical_json_bytes(material)
    except (CanonicalJSONError, TypeError, ValueError) as error:
        raise HLT17UnsupportedCallableBinding(
            "explicit synthetic callable binding is not canonical JSON"
        ) from error
    return material


def _frozen_dataclass_callable_material(
    value: object, *, depth: int, seen: set[int]
) -> dict[str, object]:
    params = getattr(value, "__dataclass_params__", None)
    if params is None or params.frozen is not True:
        _fail_binding("callable dataclass is not frozen")
    payload = {
        name: _binding_atom(getattr(value, name), depth=depth + 1, seen=seen)
        for name in (item.name for item in fields(value))
    }
    call = getattr(type(value), "__call__", None)
    func = getattr(call, "__func__", call)
    if type(func) is not types.FunctionType:
        _fail_binding("frozen dataclass callable has no python __call__")
    return {
        "kind": "frozen_dataclass_callable",
        "type": f"{type(value).__module__}.{type(value).__qualname__}",
        "fields": payload,
        "call": _function_material(func, depth=depth, seen=seen),
        "physical_source_authenticated": False,
    }


def _callable_object_material(
    value: object, *, depth: int, seen: set[int]
) -> dict[str, object]:
    marker = id(value)
    if marker in seen:
        _fail_binding("callable binding encountered a cyclic object")
    seen = set(seen)
    seen.add(marker)
    if hasattr(value, EXPLICIT_BINDING_SEAM):
        return {
            "kind": "explicit_synthetic_binding",
            "type": f"{type(value).__module__}.{type(value).__qualname__}",
            "material": _explicit_binding_material(value),
            "physical_source_authenticated": False,
        }
    if is_dataclass(value) and not isinstance(value, type):
        return _frozen_dataclass_callable_material(value, depth=depth, seen=seen)
    _fail_binding(
        "opaque or mutable callable object has no inspectable immutable "
        "configuration; supply hlt17_synthetic_callable_binding or a frozen "
        "dataclass"
    )
    raise AssertionError("unreachable")


def synthetic_callable_binding_material(
    rhs: object, projector: object | None
) -> dict[str, object]:
    """Inspect immutable python callable identity. Not a physical source proof."""

    if not callable(rhs):
        _fail_binding("rhs callback is not callable")
    if projector is not None and not callable(projector):
        _fail_binding("projector callback is not callable")
    return {
        "rhs": _binding_atom(rhs, depth=0, seen=set()),
        "projector": (
            None if projector is None else _binding_atom(projector, depth=0, seen=set())
        ),
        "physical_source_authenticated": False,
        "source_manifest_authenticated": False,
    }


def synthetic_callable_binding_digest(
    rhs: Callable[..., object],
    projector: Callable[..., object] | None,
) -> str:
    """Digest of inspectable synthetic callable binding, not source authentication."""

    return sha256(
        canonical_json_bytes(synthetic_callable_binding_material(rhs, projector))
    ).hexdigest()


def live_boundary_identities(
    *,
    state: EvolutionState,
    rhs: Callable[..., object],
    projector: Callable[..., object] | None,
    transaction: GR0RuntimeStageTransaction,
    tracers: object,
    coordinates: object,
) -> dict[str, str]:
    return {
        "physical_state_sha256": PINNED_IMP1_STATE_HASH(state),
        "coordinates_sha256": array_content_sha256(
            PINNED_IMP1_COORDINATES(coordinates, state)
        ),
        "synthetic_callable_binding_sha256": synthetic_callable_binding_digest(
            rhs, projector
        ),
        "monitor_sha256": PINNED_IMP1_DIGEST(transaction.state),
        "causal_state_sha256": PINNED_IMP1_DIGEST(transaction.causal_state),
        "tracer_sha256": PINNED_IMP1_DIGEST(PINNED_TDG6_TRACER_SNAPSHOT(tracers)),
        "guard_configuration_sha256": PINNED_IMP1_DIGEST(
            PINNED_IMP1_GUARD_CONFIGURATION(transaction)
        ),
    }


@dataclass(frozen=True, slots=True)
class HLT17LiveBoundary:
    """In-memory accepted boundary. Callbacks are not a durable store identity."""

    member_key: str
    state: EvolutionState
    transaction: GR0RuntimeStageTransaction
    tracers: object
    coordinates: object
    temporal_ledger: TDG11IMP1Ledger
    rhs: Callable[[float, EvolutionState], EvolutionRHS] = field(repr=False, compare=False)
    projector: Callable[[float, EvolutionState], EvolutionState] | None = field(
        default=None, repr=False, compare=False
    )
    previous_step_index: int = 0
    previous_transaction_serial: int = 0
    descriptor_sha256: str = ""
    event_target: float = 0.0
    limits: IMP1EnclosureLimits | None = None

    def __post_init__(self) -> None:
        spec = hlt17_member_spec(self.member_key)
        if not isinstance(self.state, EvolutionState):
            raise TypeError("state must be EvolutionState")
        if self.state.shape[0] != spec.point_count:
            _fail_premise(
                "live grid point_count differs from the declared HLT17 member identity"
            )
        try:
            coords = PINNED_IMP1_COORDINATES(self.coordinates, self.state)
        except (TypeError, ValueError) as error:
            raise HLT17InvalidPremise(
                "live coordinates are not an ordered finite member grid"
            ) from error
        if coords.size != spec.point_count:
            _fail_premise(
                "live coordinates differ from the declared HLT17 member point_count"
            )
        if spec.identity_scope == IDENTITY_SCOPE_LIVE and coords.size != spec.point_count:
            _fail_premise("production member must use its exact live point_count")
        if type(self.descriptor_sha256) is not str or len(self.descriptor_sha256) != 64:
            _fail_premise("descriptor_sha256 must be a SHA-256 digest")
        if isinstance(self.event_target, bool) or not isinstance(self.event_target, Real):
            raise TypeError("event_target must be a finite real scalar")
        if not math.isfinite(float(self.event_target)):
            _fail_premise("event_target must be finite")
        if type(self.previous_step_index) is not int or self.previous_step_index < 0:
            raise TypeError("previous_step_index must be a nonnegative built-in integer")
        if (
            type(self.previous_transaction_serial) is not int
            or self.previous_transaction_serial < 0
        ):
            raise TypeError(
                "previous_transaction_serial must be a nonnegative built-in integer"
            )
        if self.limits is not None:
            expected = spec.enclosure_limits()
            if self.limits.as_mapping() != expected.as_mapping():
                _fail_premise("live IMP1 limits differ from the declared member spec")
        synthetic_callable_binding_digest(self.rhs, self.projector)

    def identities(self) -> dict[str, str]:
        return live_boundary_identities(
            state=self.state,
            rhs=self.rhs,
            projector=self.projector,
            transaction=self.transaction,
            tracers=self.tracers,
            coordinates=self.coordinates,
        )

    def resolved_limits(self) -> IMP1EnclosureLimits:
        expected = imp1_enclosure_limits_for_member(self.member_key)
        if self.limits is None:
            return expected
        if self.limits.as_mapping() != expected.as_mapping():
            _fail_premise("live IMP1 limits differ from the declared member spec")
        return expected


@dataclass(frozen=True, slots=True)
class HLT17ReplayHandle:
    """Authenticated reconstruction of one persisted C1R1 rejection.

    ``temporal_successor.plan`` is the exact IMP1-wire half-cap. Overlay
    reductions must not replace that plan. Reconstruction is C1R1-named.
    """

    cursor: HLT17IMP1Cursor
    reconstructed: TDG11C1R1PreparedRuntime
    temporal_successor: C1R1RetrySuccessor
    sink_writes: int = 0


@dataclass(frozen=True, slots=True)
class HLT17BridgeResult:
    disposition: str
    cursor: HLT17IMP1Cursor
    prepared: TDG11C1R1PreparedRuntime | None = None
    replay: HLT17ReplayHandle | None = None
    temporal_successor: C1R1RetrySuccessor | None = None
    attempted_plan: TDG7Binary64SubdivisionPlan | None = None
    next_plan: TDG7Binary64SubdivisionPlan | None = None
    committed: TDG11C1R1CommittedRuntime | None = None
    evidence: Mapping[str, object] = field(default_factory=dict)
    sink_writes: int = 0
    accepted_state_advanced: bool = False

    def __post_init__(self) -> None:
        if self.disposition not in DISPOSITIONS:
            raise HLT17IMP1Error("unknown HLT17 bridge disposition")
        if self.sink_writes < 0 or type(self.sink_writes) is not int:
            raise TypeError("sink_writes must be a nonnegative built-in integer")
        if type(self.accepted_state_advanced) is not bool:
            raise TypeError("accepted_state_advanced must be a built-in bool")
        if self.accepted_state_advanced is True and self.disposition != "accepted_fine":
            raise HLT17IMP1Error(
                "accepted_state_advanced is only valid for accepted_fine"
            )
        if self.disposition == "accepted_fine":
            if (
                self.prepared is None
                or self.committed is None
                or self.accepted_state_advanced is not True
            ):
                raise HLT17IMP1Error(
                    "accepted_fine requires prepared and committed runtimes and the accepted-state flag"
                )
            if self.sink_writes != 0:
                raise HLT17IMP1Error("accepted_fine must not write a rejection sink")
            if type(self.committed) is not TDG11C1R1CommittedRuntime:
                raise TypeError("committed must be TDG11C1R1CommittedRuntime")
            if self.committed.implementation_id != HLT17_IMPLEMENTATION_ID:
                raise TypeError("committed implementation_id is not the fixed C1R1 object")
        elif self.committed is not None:
            raise HLT17IMP1Error(
                f"{self.disposition} must not carry a committed runtime"
            )
        if self.disposition in {
            "prepared_fresh",
            "prepared_successor",
            "prepared_overlay",
            "temporal_retry_required",
        }:
            if self.prepared is None:
                raise HLT17IMP1Error(f"{self.disposition} requires a prepared family")
            if self.accepted_state_advanced:
                raise HLT17IMP1Error(
                    f"{self.disposition} must not claim accepted-state advance"
                )
        if self.disposition == "replay_authenticated":
            if self.replay is None or self.sink_writes != 0:
                raise HLT17IMP1Error(
                    "replay reconstruction requires a handle and must not write a sink"
                )
            if self.accepted_state_advanced:
                raise HLT17IMP1Error("replay must not adopt a real state")
        if self.disposition == "temporal_retry_required":
            if self.sink_writes != 1 or self.next_plan is None:
                raise HLT17IMP1Error(
                    "temporal_retry_required requires one sink write and a successor plan"
                )
        if self.disposition == "temporal_retry_exhausted":
            if self.accepted_state_advanced:
                raise HLT17IMP1Error("temporal exhaustion must not adopt a real state")
            if "updated_imp1_ledger" not in self.evidence:
                raise HLT17IMP1Error(
                    "temporal exhaustion must retain the updated IMP1 ledger as evidence"
                )
            if self.evidence.get("executable_retry_cursor") is not False:
                raise HLT17IMP1Error(
                    "temporal exhaustion must not advertise an executable retry cursor"
                )
        if self.disposition in {
            "source_retry_required",
            "cfl_retry_required",
            "source_retry_exhausted",
            "cfl_retry_exhausted",
            "invalid_premise",
        }:
            if self.accepted_state_advanced:
                raise HLT17IMP1Error(
                    f"{self.disposition} must not claim accepted-state advance"
                )
        if self.disposition == "invalid_premise" and self.prepared is not None:
            raise HLT17IMP1Error(
                "invalid premise must not carry an adoptable prepared family"
            )


def authenticate_live_boundary(
    cursor: HLT17IMP1Cursor, boundary: HLT17LiveBoundary
) -> None:
    """Refuse stale, swapped, or retargeted live inputs before any C1R1 prepare."""

    revalidate_hlt17_cursor(cursor)
    if not isinstance(boundary, HLT17LiveBoundary):
        raise TypeError("boundary must be HLT17LiveBoundary")
    expected_identity = hlt17_implementation_identity()
    if (
        cursor.implementation_id != expected_identity["implementation_id"]
        or cursor.admission_runtime_id != expected_identity["admission_runtime_id"]
    ):
        _fail_premise("live cursor implementation identity is not the fixed C1R1 runtime")
    if boundary.member_key != cursor.member_key:
        _fail_premise("live member differs from the HLT17 cursor")
    if boundary.temporal_ledger.method != cursor.method:
        _fail_premise("live IMP1 method differs from the HLT17 cursor")
    if ledger_sha256(boundary.temporal_ledger) != ledger_sha256(cursor.ledger):
        _fail_premise("live IMP1 ledger differs from the HLT17 cursor")
    if PINNED_IMP1_LEDGER_DIGEST(boundary.temporal_ledger) != PINNED_IMP1_LEDGER_DIGEST(
        cursor.ledger
    ):
        _fail_premise("live IMP1 ledger digest differs from the pinned IMP1 formula")
    if not (
        boundary.event_target.hex() == cursor.event_target.hex()
    ):
        _fail_premise("live event target differs from the HLT17 cursor")
    identities = boundary.identities()
    if identities["physical_state_sha256"] != cursor.physical_state_sha256:
        _fail_premise("live physical u/p/q hash differs from the HLT17 cursor")
    if boundary.descriptor_sha256 != cursor.descriptor_sha256:
        _fail_premise("live serialized descriptor hash differs from the HLT17 cursor")
    if identities["physical_state_sha256"] == boundary.descriptor_sha256:
        _fail_premise("descriptor and physical u/p/q identities are not interchangeable")
    for name in (
        "coordinates_sha256",
        "synthetic_callable_binding_sha256",
        "monitor_sha256",
        "causal_state_sha256",
        "tracer_sha256",
        "guard_configuration_sha256",
    ):
        if identities[name] != getattr(cursor, name):
            _fail_premise(f"live {name} differs from the HLT17 cursor")
    if boundary.previous_step_index != cursor.accepted_step_index:
        _fail_premise("live step index differs from the HLT17 cursor")
    if boundary.previous_transaction_serial != cursor.accepted_transaction_serial:
        _fail_premise("live transaction serial differs from the HLT17 cursor")
    boundary.resolved_limits()
    try:
        PINNED_IMP1_CHECK_CURRENT_BOUNDARY(
            method=cursor.method,
            time=cursor.accepted_time,
            state=boundary.state,
            transaction=boundary.transaction,
            temporal_ledger=boundary.temporal_ledger,
            previous_step_index=boundary.previous_step_index,
            previous_transaction_serial=boundary.previous_transaction_serial,
        )
    except (TypeError, ValueError) as error:
        raise HLT17InvalidPremise(
            "IMP1 accepted boundary check failed before preparation"
        ) from error


def _runtime_fingerprint(
    transaction: GR0RuntimeStageTransaction, tracers: object
) -> tuple[object, ...]:
    snapshot = PINNED_TDG6_TRACER_SNAPSHOT(tracers)
    return (
        transaction.state,
        transaction.causal_state,
        PINNED_IMP1_GUARD_CONFIGURATION(transaction),
        snapshot,
        transaction.state.accepted_stage_count,
        transaction.state.last_transaction_serial,
        transaction.state.last_accepted_time,
        snapshot.labels_sha256,
        snapshot.positions_sha256,
        snapshot.proper_times_sha256,
        snapshot.event_history_sha256,
        snapshot.event_count,
    )


def _capture_runtime_boundary(boundary: HLT17LiveBoundary) -> dict[str, object]:
    transaction = boundary.transaction
    return {
        "monitor": transaction.state,
        "causal": transaction.causal_state,
        "guard": PINNED_IMP1_GUARD_CONFIGURATION(transaction),
        "tracers": PINNED_TDG6_CLONE_TRACERS(boundary.tracers),
        "fingerprint": _runtime_fingerprint(transaction, boundary.tracers),
    }


def _restore_runtime_boundary(
    boundary: HLT17LiveBoundary, captured: Mapping[str, object]
) -> None:
    transaction = boundary.transaction
    transaction.state = captured["monitor"]  # type: ignore[assignment]
    transaction.causal_state = captured["causal"]  # type: ignore[assignment]
    (
        transaction.thresholds,
        transaction.boundary_geometry,
        transaction.grid_spacing,
        transaction.cfl_maximum,
        transaction.hat_normal_factor,
    ) = captured["guard"]  # type: ignore[misc]
    PINNED_TDG6_RESTORE_TRACERS(boundary.tracers, captured["tracers"])
    if _runtime_fingerprint(transaction, boundary.tracers) != captured["fingerprint"]:
        raise RuntimeError("HLT17 complete accepted-boundary rollback failed")


def _bind_prepared_to_current(
    cursor: HLT17IMP1Cursor,
    boundary: HLT17LiveBoundary,
    prepared: object,
) -> TDG11C1R1PreparedRuntime:
    """Bind prepared evidence to the CURRENT cursor/boundary, not its own claims."""

    authenticate_live_boundary(cursor, boundary)
    try:
        prepared = require_hlt17_prepared(prepared)
    except TypeError:
        raise
    except (ValueError,) as error:
        raise HLT17InvalidPremise(
            "prepared family failed C1R1 validation before adoption"
        ) from error
    if prepared.original_rhs is not boundary.rhs:
        _fail_premise("prepared rhs differs from the current live boundary")
    if prepared.original_projector is not boundary.projector:
        _fail_premise("prepared projector differs from the current live boundary")
    if prepared.implementation_id != cursor.implementation_id:
        _fail_premise("prepared implementation_id differs from the current cursor")
    if prepared.method != cursor.method:
        _fail_premise("prepared method differs from the current cursor")
    if prepared.plan.current.hex() != cursor.accepted_time.hex():
        _fail_premise("prepared current time differs from the cursor accepted time")
    if prepared.plan.event_target.hex() != cursor.event_target.hex():
        _fail_premise("prepared event target differs from the current cursor")
    if prepared.previous_step_index != boundary.previous_step_index:
        _fail_premise("prepared step index differs from the current live boundary")
    if prepared.previous_transaction_serial != boundary.previous_transaction_serial:
        _fail_premise(
            "prepared transaction serial differs from the current live boundary"
        )
    prepared_ledger = PINNED_IMP1_LEDGER_DIGEST(prepared.original_ledger)
    current_ledger = PINNED_IMP1_LEDGER_DIGEST(boundary.temporal_ledger)
    if prepared_ledger != current_ledger:
        _fail_premise(
            "stale predecessor-ledger preparation cannot be admitted or committed"
        )
    next_plan = cursor.next_attempt_plan()
    if next_plan is not None and prepared.plan != next_plan:
        _fail_premise(
            "prepared TDG7 plan differs from the exact cursor next-attempt plan"
        )
    try:
        PINNED_C1R1_REVALIDATE_BOUNDARY(
            prepared,
            transaction=boundary.transaction,
            tracers=boundary.tracers,
            temporal_ledger=boundary.temporal_ledger,
            current_time=cursor.accepted_time,
            current_state=boundary.state,
            current_step_index=boundary.previous_step_index,
            current_transaction_serial=boundary.previous_transaction_serial,
            coordinates=boundary.coordinates,
        )
    except (TypeError, ValueError) as error:
        raise HLT17InvalidPremise(
            "prepared family is not bound to the current cursor and boundary"
        ) from error
    return prepared


def _validate_committed_against_cursor(
    committed: TDG11C1R1CommittedRuntime, cursor: HLT17IMP1Cursor
) -> None:
    committed = require_hlt17_committed(committed)
    if committed.method != cursor.method:
        _fail_premise("committed method differs from the returned cursor")
    if committed.time.hex() != cursor.accepted_time.hex():
        _fail_premise("committed time differs from the returned cursor")
    if committed.step_index != cursor.accepted_step_index:
        _fail_premise("committed step index differs from the returned cursor")
    if committed.transaction_serial != cursor.accepted_transaction_serial:
        _fail_premise("committed transaction serial differs from the returned cursor")
    if PINNED_IMP1_STATE_HASH(committed.state) != cursor.physical_state_sha256:
        _fail_premise("committed physical state differs from the returned cursor")
    if PINNED_IMP1_LEDGER_DIGEST(committed.temporal_ledger) != PINNED_IMP1_LEDGER_DIGEST(
        cursor.ledger
    ):
        _fail_premise("committed IMP1 ledger differs from the returned cursor")


def with_ledger(
    boundary: HLT17LiveBoundary, ledger: TDG11IMP1Ledger
) -> HLT17LiveBoundary:
    """Rebind the live IMP1 ledger after a cursor transition; no store write."""

    return HLT17LiveBoundary(
        member_key=boundary.member_key,
        state=boundary.state,
        rhs=boundary.rhs,
        projector=boundary.projector,
        transaction=boundary.transaction,
        tracers=boundary.tracers,
        coordinates=boundary.coordinates,
        temporal_ledger=ledger,
        previous_step_index=boundary.previous_step_index,
        previous_transaction_serial=boundary.previous_transaction_serial,
        descriptor_sha256=boundary.descriptor_sha256,
        event_target=boundary.event_target,
        limits=boundary.limits,
    )


def seed_cursor_from_boundary(boundary: HLT17LiveBoundary) -> HLT17IMP1Cursor:
    identities = boundary.identities()
    return seed_hlt17_cursor(
        member_key=boundary.member_key,
        event_target=boundary.event_target,
        ledger=boundary.temporal_ledger,
        descriptor_sha256=boundary.descriptor_sha256,
        physical_state_sha256=identities["physical_state_sha256"],
        coordinates_sha256=identities["coordinates_sha256"],
        synthetic_callable_binding_sha256=identities["synthetic_callable_binding_sha256"],
        monitor_sha256=identities["monitor_sha256"],
        causal_state_sha256=identities["causal_state_sha256"],
        tracer_sha256=identities["tracer_sha256"],
        guard_configuration_sha256=identities["guard_configuration_sha256"],
    )


def _prepare_with_plan(
    *,
    plan: TDG7Binary64SubdivisionPlan,
    boundary: HLT17LiveBoundary,
    ledger: TDG11IMP1Ledger,
) -> TDG11C1R1PreparedRuntime:
    validate_tdg7_binary64_subdivision_plan(plan)
    limits = boundary.resolved_limits()
    try:
        PINNED_IMP1_CHECK_CURRENT_BOUNDARY(
            method=boundary.temporal_ledger.method,
            time=plan.current,
            state=boundary.state,
            transaction=boundary.transaction,
            temporal_ledger=ledger,
            previous_step_index=boundary.previous_step_index,
            previous_transaction_serial=boundary.previous_transaction_serial,
        )
        prepared = PINNED_C1R1_PREPARE(
            plan=plan,
            method=boundary.temporal_ledger.method,
            state=boundary.state,
            rhs=boundary.rhs,
            projector=boundary.projector,
            transaction=boundary.transaction,
            tracers=boundary.tracers,
            coordinates=boundary.coordinates,
            temporal_ledger=ledger,
            previous_step_index=boundary.previous_step_index,
            previous_transaction_serial=boundary.previous_transaction_serial,
            limits=limits,
        )
        return require_hlt17_prepared(prepared)
    except (
        C1R1RefinementPathStop,
        C1R1RuntimeResourceStop,
        C1R1CoordinateLatticeStop,
    ):
        raise
    except (TypeError, ValueError) as error:
        raise HLT17InvalidPremise(
            "C1R1 reduced-plan preparation failed before a successor existed"
        ) from error


def _canonical_rejection(record: Mapping[str, object] | TDG11IMP1TemporalRejection) -> bytes:
    mapping = record.as_mapping() if isinstance(record, TDG11IMP1TemporalRejection) else dict(record)
    return canonical_json_bytes(mapping)


def _path_owner(stop: C1R1RefinementPathStop) -> str | None:
    if stop.source_retry is not None:
        return SOURCE_OWNER
    if isinstance(stop.cause, CFLRetryRequired):
        return CFL_OWNER
    return None


def _overlay_evidence(
    stop: C1R1RefinementPathStop, plan: TDG7Binary64SubdivisionPlan, owner: str
) -> dict[str, object]:
    common = {
        "path": stop.path,
        "attempted_requested_cap_hex": plan.requested_cap.hex(),
        "attempted_macro_width_hex": plan.macro_width.hex(),
        "accepted_boundary_preserved": True,
        "physical_classification": False,
    }
    if owner == SOURCE_OWNER:
        if stop.source_retry is None:
            _fail_premise("source overlay omitted PROTO7 source-retry evidence")
        return {
            "kind": "source_retry",
            "source_retry": PINNED_IMP1_WIRE(stop.source_retry),
            **common,
        }
    if not isinstance(stop.cause, CFLRetryRequired):
        _fail_premise("CFL overlay omitted CFLRetryRequired evidence")
    return {
        "kind": "cfl_retry",
        "observed_ratio_hex": float(stop.cause.observed_ratio).hex(),
        "maximum_ratio_hex": float(stop.cause.maximum_ratio).hex(),
        **common,
    }


def _overlay_result(
    cursor: HLT17IMP1Cursor,
    *,
    owner: str,
    attempted_plan: TDG7Binary64SubdivisionPlan,
    stop: C1R1RefinementPathStop,
) -> HLT17BridgeResult:
    evidence = _overlay_evidence(stop, attempted_plan, owner)
    updated = close_source_cfl_overlay(
        cursor,
        owner=owner,
        attempted_plan=attempted_plan,
        evidence=evidence,
    )
    overlay = updated.overlay
    assert isinstance(overlay, HLT17RetryOverlay)
    if owner == SOURCE_OWNER:
        disposition = (
            "source_retry_exhausted" if overlay.exhausted else "source_retry_required"
        )
    else:
        disposition = (
            "cfl_retry_exhausted" if overlay.exhausted else "cfl_retry_required"
        )
    return HLT17BridgeResult(
        disposition=disposition,
        cursor=updated,
        attempted_plan=attempted_plan,
        next_plan=overlay.next_plan,
        evidence=dict(overlay.evidence),
        sink_writes=0,
        accepted_state_advanced=False,
    )


def _invalid_prepare_result(
    cursor: HLT17IMP1Cursor,
    *,
    attempted_plan: TDG7Binary64SubdivisionPlan | None,
    error: BaseException,
) -> HLT17BridgeResult:
    return HLT17BridgeResult(
        disposition="invalid_premise",
        cursor=cursor,
        attempted_plan=attempted_plan,
        evidence={
            "kind": "invalid_premise",
            "exception_type": type(error).__name__,
            "reason": str(error),
            "spends_retry": False,
            "spends_debit": False,
        },
        sink_writes=0,
        accepted_state_advanced=False,
    )


def _prepared_evidence(kind: str, prepared: TDG11C1R1PreparedRuntime) -> dict[str, object]:
    identity = hlt17_implementation_identity()
    return {
        "kind": kind,
        "requested_cap_hex": prepared.plan.requested_cap.hex(),
        "macro_width_hex": prepared.plan.macro_width.hex(),
        "preparation_sha256": prepared.receipt_sha256,
        "assessment_sha256": prepared.assessment.assessment_sha256,
        "admission_passed": prepared.assessment.admission_passed,
        "implementation_id": prepared.implementation_id,
        "admission_runtime_id": identity["admission_runtime_id"],
        "reference_wire_artifact_id": identity["reference_wire_artifact_id"],
        "reference_wire_evaluator_id": prepared.reference_wire_evaluator_id,
    }


def prepare_hlt17_fresh(
    cursor: HLT17IMP1Cursor,
    boundary: HLT17LiveBoundary,
    *,
    requested_cap: object,
) -> HLT17BridgeResult:
    """Public C1R1 fresh entrypoint; refused at any active retry depth or overlay."""

    authenticate_live_boundary(cursor, boundary)
    if not cursor.public_fresh_permitted():
        _fail_premise("public C1R1 fresh is permitted only at zero active retry depth")
    try:
        plan = PINNED_C1R1_PLAN(
            cursor.accepted_time, cursor.event_target, requested_cap
        )
    except C1R1CoordinateLatticeStop as error:
        return _invalid_prepare_result(cursor, attempted_plan=None, error=error)
    if plan.requested_cap.hex() != float(requested_cap).hex():
        # requested_cap is an identity, even when the lattice floors the width.
        if plan.requested_cap != float(requested_cap):
            _fail_premise("fresh TDG7 plan did not retain the requested cap identity")
    try:
        prepared = prepare_hlt17_admission_fresh(
            method=cursor.method,
            time=cursor.accepted_time,
            event_target=cursor.event_target,
            requested_cap=requested_cap,
            state=boundary.state,
            rhs=boundary.rhs,
            projector=boundary.projector,
            transaction=boundary.transaction,
            tracers=boundary.tracers,
            coordinates=boundary.coordinates,
            temporal_ledger=boundary.temporal_ledger,
            previous_step_index=boundary.previous_step_index,
            previous_transaction_serial=boundary.previous_transaction_serial,
            limits=boundary.resolved_limits(),
        )
    except C1R1RefinementPathStop as stop:
        owner = _path_owner(stop)
        if owner is None:
            return _invalid_prepare_result(cursor, attempted_plan=plan, error=stop)
        return _overlay_result(cursor, owner=owner, attempted_plan=plan, stop=stop)
    except (C1R1RuntimeResourceStop, C1R1CoordinateLatticeStop) as error:
        return _invalid_prepare_result(cursor, attempted_plan=plan, error=error)
    if prepared.plan != plan:
        _fail_premise("public C1R1 fresh plan differs from the certified TDG7 plan")
    return HLT17BridgeResult(
        disposition="prepared_fresh",
        cursor=cursor,
        prepared=prepared,
        attempted_plan=prepared.plan,
        next_plan=None,
        evidence=_prepared_evidence("prepared_fresh", prepared),
        sink_writes=0,
    )


def replay_hlt17_rejection(
    cursor: HLT17IMP1Cursor, boundary: HLT17LiveBoundary
) -> HLT17ReplayHandle:
    """Reconstruct the exact saved rejected proposal; never a new retry or sink write.

    Reconstruction uses the predecessor IMP1 ledger and the complete TDG7 plan
    (requested_cap, not executed macro_width). The entire canonical
    rejection/preparation receipt is compared.
    """

    authenticate_live_boundary(cursor, boundary)
    if cursor.mode != "RETRY_PENDING" or cursor.temporal is None:
        _fail_premise("replay requires a RETRY_PENDING cursor with temporal identity")
    temporal = cursor.temporal
    saved_plan = temporal.predecessor_plan
    try:
        replanned = PINNED_C1R1_PLAN(
            saved_plan.current, saved_plan.event_target, saved_plan.requested_cap
        )
    except C1R1CoordinateLatticeStop as error:
        raise HLT17InvalidPremise(
            "persisted predecessor plan no longer replans from its requested cap"
        ) from error
    if replanned != saved_plan:
        _fail_premise("persisted predecessor plan is not the certified requested-cap lattice")
    if saved_plan.requested_cap.hex() == saved_plan.macro_width.hex():
        pass
    else:
        width_only = PINNED_C1R1_PLAN(
            saved_plan.current, saved_plan.event_target, saved_plan.macro_width
        )
        if width_only == saved_plan:
            _fail_premise(
                "requested_cap identity was lost; width-only reconstruction is insufficient"
            )
    try:
        reconstructed = _prepare_with_plan(
            plan=saved_plan,
            boundary=boundary,
            ledger=temporal.predecessor_ledger,
        )
    except C1R1RefinementPathStop as error:
        raise HLT17InvalidPremise(
            "replay source/health failure is an invalid premise, not a new attempt"
        ) from error
    except C1R1RuntimeResourceStop as error:
        raise HLT17InvalidPremise(
            "replay resource failure is an invalid premise, not a new attempt"
        ) from error
    except C1R1CoordinateLatticeStop as error:
        raise HLT17InvalidPremise(
            "replay lattice failure is an invalid premise, not a new attempt"
        ) from error
    if reconstructed.assessment.admission_passed is True:
        _fail_premise("persisted rejection now replays as an admitted predecessor")
    if reconstructed.receipt_sha256 != temporal.preparation_sha256:
        _fail_premise("replay preparation receipt differs from the persisted cursor")
    if reconstructed.assessment.assessment_sha256 != temporal.assessment_sha256:
        _fail_premise("replay assessment receipt differs from the persisted cursor")
    if reconstructed.plan != saved_plan:
        _fail_premise("replay executed a different TDG7 plan than the persisted predecessor")
    reconstructed_record = hlt17_rejection_record(reconstructed)
    if _canonical_rejection(reconstructed_record) != _canonical_rejection(temporal.rejection):
        _fail_premise("replay IMP1-wire rejection receipt differs from the persisted cursor")
    expected_updated = append_imp1_rejection(
        temporal.predecessor_ledger, reconstructed_record
    )
    if PINNED_IMP1_LEDGER_DIGEST(expected_updated) != PINNED_IMP1_LEDGER_DIGEST(
        cursor.ledger
    ):
        _fail_premise("replay updated IMP1 ledger differs from the persisted cursor")
    if hlt17_receipt_digest(reconstructed) != reconstructed.receipt_sha256:
        _fail_premise("replay receipt digest is not redundant")
    try:
        successor = make_hlt17_retry_successor(
            reconstructed,
            temporal.successor_plan,
            cursor.ledger,
            _canonical_rejection(temporal.rejection),
        )
    except (TypeError, ValueError) as error:
        raise HLT17InvalidPremise(
            "replay cannot reconstruct the exact C1R1 half-cap successor"
        ) from error
    if successor.plan != temporal.successor_plan:
        _fail_premise("C1R1RetrySuccessor.plan is not the persisted half-cap identity")
    return HLT17ReplayHandle(
        cursor=cursor,
        reconstructed=reconstructed,
        temporal_successor=successor,
        sink_writes=0,
    )


def prepare_hlt17_next(
    cursor: HLT17IMP1Cursor,
    boundary: HLT17LiveBoundary,
    *,
    replay: HLT17ReplayHandle | None = None,
) -> HLT17BridgeResult:
    """Prepare the temporal successor or overlay plan after authenticated replay."""

    authenticate_live_boundary(cursor, boundary)
    if cursor.public_fresh_permitted():
        _fail_premise("zero-depth fresh attempts must use prepare_hlt17_fresh")
    handle = replay
    if cursor.mode != FRESH_READY:
        if handle is None:
            handle = replay_hlt17_rejection(cursor, boundary)
        elif handle.cursor.kernel_transition_digest != cursor.kernel_transition_digest:
            _fail_premise("replay handle is stale relative to the HLT17 cursor")
        if handle.temporal_successor.plan != cursor.temporal.successor_plan:  # type: ignore[union-attr]
            _fail_premise("replay handle replaced C1R1RetrySuccessor.plan")
    plan = cursor.next_attempt_plan()
    if plan is None:
        owner = cursor.latest_owner or SOURCE_OWNER
        raise HLT17OwnerRetryExhausted(
            owner,
            "no_executable_next_plan",
            attempted_plan=None if cursor.overlay is None else cursor.overlay.attempted_plan,
            ledger=cursor.ledger,
        )
    attempted = plan
    try:
        if cursor.overlay is None and handle is not None:
            prepared = prepare_hlt17_admission_retry(
                handle.temporal_successor,
                state=boundary.state,
                rhs=boundary.rhs,
                projector=boundary.projector,
                transaction=boundary.transaction,
                tracers=boundary.tracers,
                coordinates=boundary.coordinates,
            )
            disposition = "prepared_successor"
            successor = handle.temporal_successor
        else:
            if cursor.overlay is not None and handle is not None:
                if handle.temporal_successor.plan == cursor.overlay.next_plan:
                    _fail_premise("overlay plan must not replace C1R1RetrySuccessor.plan")
            prepared = _prepare_with_plan(
                plan=plan,
                boundary=boundary,
                ledger=cursor.ledger,
            )
            disposition = "prepared_overlay" if cursor.overlay is not None else "prepared_successor"
            successor = handle.temporal_successor if handle is not None else None
    except C1R1RefinementPathStop as stop:
        owner = _path_owner(stop)
        if owner is None:
            return _invalid_prepare_result(cursor, attempted_plan=attempted, error=stop)
        return _overlay_result(cursor, owner=owner, attempted_plan=attempted, stop=stop)
    except (C1R1RuntimeResourceStop, C1R1CoordinateLatticeStop) as error:
        return _invalid_prepare_result(cursor, attempted_plan=attempted, error=error)
    if prepared.plan != attempted:
        _fail_premise("next attempt executed a different TDG7 plan than the cursor")
    if successor is not None and cursor.overlay is not None:
        if successor.plan == prepared.plan and prepared.plan != cursor.temporal.successor_plan:  # type: ignore[union-attr]
            _fail_premise("overlay preparation replaced the IMP1 half-cap successor")
    return HLT17BridgeResult(
        disposition=disposition,
        cursor=cursor,
        prepared=prepared,
        replay=handle,
        temporal_successor=successor,
        attempted_plan=prepared.plan,
        next_plan=None if cursor.overlay is None else cursor.overlay.next_plan,
        evidence={
            **_prepared_evidence(disposition, prepared),
            "imp1_half_cap_hex": (
                None
                if cursor.temporal is None
                else cursor.temporal.imp1_half_cap.hex()
            ),
        },
        sink_writes=0,
    )


class _MemoryRejectionSink:
    """In-process IMP1 sink callback. Not durable publication or a store write."""

    def __init__(self) -> None:
        self.records: list[dict[str, object]] = []

    def __call__(self, record: Mapping[str, object]) -> None:
        self.records.append(dict(record))


def require_hlt17_admission(
    cursor: HLT17IMP1Cursor,
    boundary: HLT17LiveBoundary,
    prepared: TDG11C1R1PreparedRuntime,
) -> HLT17BridgeResult:
    """Require C1R1 admission against the current cursor/boundary.

    A typed nonadmission is recorded only in the in-process IMP1-wire sink
    callback. That callback is not durable publication. The rejection event
    type remains the IMP1 ledger wire.
    """

    bound = _bind_prepared_to_current(cursor, boundary, prepared)
    sink = _MemoryRejectionSink()
    try:
        require_hlt17_admission_family(
            bound,
            transaction=boundary.transaction,
            tracers=boundary.tracers,
            temporal_ledger=boundary.temporal_ledger,
            current_time=cursor.accepted_time,
            current_state=boundary.state,
            current_step_index=boundary.previous_step_index,
            current_transaction_serial=boundary.previous_transaction_serial,
            coordinates=boundary.coordinates,
            durable_rejection_sink=sink,
        )
    except C1R1TemporalRetryRequired as retry:
        if len(sink.records) != 1:
            _fail_premise("C1R1 nonadmission sink write count differs")
        if _canonical_rejection(sink.records[0]) != retry.record_bytes:
            _fail_premise("C1R1 nonadmission sink receipt differs")
        updated = close_temporal_rejection(
            cursor,
            predecessor_plan=bound.plan,
            predecessor_ledger=boundary.temporal_ledger,
            updated_ledger=retry.updated_ledger,
            rejection=retry.evidence,
            preparation_sha256=bound.receipt_sha256,
            assessment_sha256=bound.assessment.assessment_sha256,
        )
        return HLT17BridgeResult(
            disposition="temporal_retry_required",
            cursor=updated,
            prepared=bound,
            temporal_successor=None,
            attempted_plan=bound.plan,
            next_plan=updated.next_attempt_plan(),
            evidence={
                "kind": "temporal_retry_required",
                "rejection": dict(retry.evidence),
                "implementation_id": bound.implementation_id,
                "imp1_half_cap_hex": updated.temporal.imp1_half_cap.hex()  # type: ignore[union-attr]
                if updated.temporal is not None
                else None,
            },
            sink_writes=1,
        )
    except C1R1TemporalRetryExhausted as retry:
        return HLT17BridgeResult(
            disposition="temporal_retry_exhausted",
            cursor=cursor,
            prepared=bound,
            attempted_plan=bound.plan,
            evidence={
                "kind": "temporal_retry_exhausted",
                "reason": retry.reason,
                "rejection": dict(retry.evidence),
                "implementation_id": bound.implementation_id,
                "updated_imp1_ledger": imp1_checkpoint_extension(retry.updated_ledger),
                "updated_ledger_sha256": ledger_sha256(retry.updated_ledger),
                "executable_retry_cursor": False,
            },
            sink_writes=len(sink.records),
        )
    if bound.assessment.admission_passed is not True:
        _fail_premise("C1R1 require returned without admission or a typed retry")
    if sink.records:
        _fail_premise("admitted C1R1 family must not write a rejection sink")
    return HLT17BridgeResult(
        disposition="prepared_fresh" if cursor.mode == FRESH_READY and cursor.overlay is None else "prepared_successor",
        cursor=cursor,
        prepared=bound,
        attempted_plan=bound.plan,
        evidence={
            **_prepared_evidence("admitted", bound),
            "admission_passed": True,
        },
        sink_writes=0,
    )


def commit_hlt17_kernel(
    cursor: HLT17IMP1Cursor,
    boundary: HLT17LiveBoundary,
    prepared: TDG11C1R1PreparedRuntime,
) -> TDG11C1R1CommittedRuntime:
    """One C1R1 fine commit on the supplied boundary. Does not close the cursor.

    The live transaction is committed at most once. Encoding and descriptor
    finalization belong to the joint checkpoint layer so the new cursor can
    name the actual new descriptor. Late kernel failure restores this boundary.
    """

    captured = _capture_runtime_boundary(boundary)
    try:
        bound = _bind_prepared_to_current(cursor, boundary, prepared)
        if bound.assessment.admission_passed is not True:
            _fail_premise("cannot commit an unadmitted C1R1 family")
        committed = commit_hlt17_admission(
            bound,
            transaction=boundary.transaction,
            tracers=boundary.tracers,
            temporal_ledger=boundary.temporal_ledger,
            current_time=cursor.accepted_time,
            current_state=boundary.state,
            current_step_index=boundary.previous_step_index,
            current_transaction_serial=boundary.previous_transaction_serial,
            coordinates=boundary.coordinates,
        )
        committed = require_hlt17_committed(committed)
        if committed.corrected_endpoint_adopted is not False:
            _fail_premise("cannot scrape an unqualified corrected endpoint")
        return committed
    except BaseException:
        _restore_runtime_boundary(boundary, captured)
        raise


def finalize_hlt17_fine(
    cursor: HLT17IMP1Cursor,
    boundary: HLT17LiveBoundary,
    prepared: TDG11C1R1PreparedRuntime,
    committed: TDG11C1R1CommittedRuntime,
    *,
    descriptor_sha256: str,
) -> HLT17BridgeResult:
    """Close the cursor against an already-committed runtime and new descriptor.

    ``descriptor_sha256`` is the actual new accepted-encoding digest. This does
    not call C1R1 commit again. The supplied boundary must already carry the
    committed monitor/tracer/causal identities; re-binding would compare the
    old cursor to the advanced live objects.
    """

    prepared = require_hlt17_prepared(prepared)
    if prepared.assessment.admission_passed is not True:
        _fail_premise("cannot finalize an unadmitted C1R1 family")
    committed = require_hlt17_committed(committed)
    if prepared.original_rhs is not boundary.rhs:
        _fail_premise("prepared rhs differs from the current live boundary")
    if prepared.original_projector is not boundary.projector:
        _fail_premise("prepared projector differs from the current live boundary")
    identities = live_boundary_identities(
        state=committed.state,
        rhs=boundary.rhs,
        projector=boundary.projector,
        transaction=boundary.transaction,
        tracers=boundary.tracers,
        coordinates=boundary.coordinates,
    )
    updated = close_fine_adoption(
        cursor,
        committed_ledger=committed.temporal_ledger,
        descriptor_sha256=descriptor_sha256,
        physical_state_sha256=identities["physical_state_sha256"],
        coordinates_sha256=identities["coordinates_sha256"],
        monitor_sha256=identities["monitor_sha256"],
        causal_state_sha256=identities["causal_state_sha256"],
        tracer_sha256=identities["tracer_sha256"],
        guard_configuration_sha256=identities["guard_configuration_sha256"],
    )
    _validate_committed_against_cursor(committed, updated)
    if updated.descriptor_sha256 != descriptor_sha256:
        _fail_premise("finalized cursor does not name the supplied new descriptor")
    if updated.descriptor_sha256 == cursor.descriptor_sha256:
        _fail_premise("finalized cursor retained the old descriptor")
    return HLT17BridgeResult(
        disposition="accepted_fine",
        cursor=updated,
        prepared=prepared,
        attempted_plan=prepared.plan,
        committed=committed,
        evidence={
            "kind": "accepted_fine",
            "preparation_sha256": prepared.receipt_sha256,
            "assessment_sha256": prepared.assessment.assessment_sha256,
            "implementation_id": prepared.implementation_id,
            "admission_runtime_id": hlt17_implementation_identity()["admission_runtime_id"],
            "accepted_time_hex": committed.time.hex(),
            "physical_classification": False,
            "corrected_endpoint_adopted": committed.corrected_endpoint_adopted,
        },
        sink_writes=0,
        accepted_state_advanced=True,
    )


def commit_hlt17_fine(
    cursor: HLT17IMP1Cursor,
    boundary: HLT17LiveBoundary,
    prepared: TDG11C1R1PreparedRuntime,
    *,
    descriptor_sha256: str,
) -> HLT17BridgeResult:
    """Kernel-commit then finalize the cursor with the actual new descriptor.

    Callers that still need a one-shot bridge result must supply the new
    descriptor digest. The joint checkpoint layer encodes first from a
    provisional copy, then finalizes, so this convenience path is not the
    production descriptor owner.
    """

    captured = _capture_runtime_boundary(boundary)
    try:
        committed = commit_hlt17_kernel(cursor, boundary, prepared)
        return finalize_hlt17_fine(
            cursor,
            boundary,
            prepared,
            committed,
            descriptor_sha256=descriptor_sha256,
        )
    except BaseException:
        _restore_runtime_boundary(boundary, captured)
        raise


def plan_from_requested_cap_not_width(
    *, current: object, event_target: object, requested_cap: object
) -> TDG7Binary64SubdivisionPlan:
    """Lattice helper that keeps requested_cap distinct from executed width."""

    plan = plan_imp1_lattice(current, event_target, requested_cap)
    width_only = plan_imp1_lattice(current, event_target, plan.macro_width)
    if (
        plan.requested_cap.hex() != plan.macro_width.hex()
        and encode_tdg7_plan(width_only) == encode_tdg7_plan(plan)
    ):
        raise HLT17IMP1Error("requested_cap identity collapsed to executed width")
    return plan


def checkpoint_bytes(ledger: TDG11IMP1Ledger) -> bytes:
    return canonical_json_bytes(imp1_checkpoint_extension(ledger))


__all__ = [
    "DISPOSITIONS",
    "FUTURE_INTEGRATION_SEAMS",
    "EXPLICIT_BINDING_SEAM",
    "HLT17BridgeResult",
    "HLT17LiveBoundary",
    "HLT17ReplayHandle",
    "HLT17UnsupportedCallableBinding",
    "PINNED_C1R1_PLAN",
    "PINNED_C1R1_PREPARE",
    "PINNED_C1R1_PREPARED_MAPPING",
    "PINNED_C1R1_PRIVATE_DEPENDENCIES",
    "PINNED_C1R1_RECEIPT_DIGEST",
    "PINNED_C1R1_REJECTION_RECORD",
    "PINNED_C1R1_REVALIDATE_BOUNDARY",
    "PINNED_HLT17_PRIVATE_DEPENDENCIES",
    "PINNED_IMP1_CHECK_CURRENT_BOUNDARY",
    "PINNED_IMP1_COORDINATES",
    "PINNED_IMP1_DIGEST",
    "PINNED_IMP1_GUARD_CONFIGURATION",
    "PINNED_IMP1_LEDGER_DIGEST",
    "PINNED_IMP1_PRIVATE_DEPENDENCIES",
    "PINNED_IMP1_STATE_HASH",
    "PINNED_IMP1_WIRE",
    "PINNED_TDG6_CLONE_TRACERS",
    "PINNED_TDG6_RESTORE_TRACERS",
    "PINNED_TDG6_TRACER_SNAPSHOT",
    "authenticate_live_boundary",
    "commit_hlt17_fine",
    "commit_hlt17_kernel",
    "finalize_hlt17_fine",
    "live_boundary_identities",
    "plan_from_requested_cap_not_width",
    "prepare_hlt17_fresh",
    "prepare_hlt17_next",
    "replay_hlt17_rejection",
    "require_hlt17_admission",
    "seed_cursor_from_boundary",
    "synthetic_callable_binding_digest",
    "synthetic_callable_binding_material",
    "with_ledger",
]
