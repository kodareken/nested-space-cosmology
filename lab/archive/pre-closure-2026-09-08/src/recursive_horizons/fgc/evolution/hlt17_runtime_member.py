"""Independent HLT17 in-memory member carrier for prospective PROTO19.

This is a new runtime object.  It is not a Proto14 subclass, not a relabeled
HLT16 member, and not an execution authority.  Construction consumes already
provided static initial/operator/projector/transaction/tracer objects and a
separately typed IMP1 ledger.  Fresh/retry/commit use the fixed C1R1
admission runtime; they do not call sealed IMP1 ``_prepare``.  It does not
reconstruct a source or grid, and it does not authenticate an HLT15 origin,
source/operator closure, cursor, or durable store.

Pinned private dependencies, with no public clone/restore/snapshot API:

* ``tdg6_temporal_admission_runtime._clone_tracers``
* ``tdg6_temporal_admission_runtime._restore_tracers``
* ``tdg6_temporal_admission_runtime._tracer_snapshot``

Those sealed helpers are referenced, not patched.  Historical
``Proto14RunMember.advance_to`` is never called.

Integration seams owned elsewhere:

* HLT17 cursor/bridge authenticates the complete predecessor cursor.  This
  carrier never hashes a small identity projection and calls it the cursor.
* Durable journal/store binds closure, origin, and full-cursor references.
  ``cursor_generation``, checkpoint/journal revisions, and retry-plan overlays
  are not accepted physical-state generation.
* The sealed HLT15 origin protocol is ``FGC-2-SF1-PROTO17``.
* Live production authority admits only RK4 2049/4097/8193 and SSPRK3
  4097/8193/16385.  This carrier can hold a true small synthetic identity.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
import math
from numbers import Real
import re
from typing import Any, Callable, Mapping

import numpy as np

from recursive_horizons.evidence_io import CanonicalJSONError, canonical_json_bytes

from .boundary_domain import CausalBudgetState
from .numerical_engine import (
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionState,
    array_content_sha256,
)
from .proto5_runtime import GR0RuntimeMonitorState, GR0RuntimeStageTransaction
from .tdg11_imp1_enclosure import IMP1EnclosureLimits
from .tdg11_imp1_ledger import (
    TDG11IMP1Ledger,
    inherited_tdg6_ledger,
)
from .hlt17_admission_runtime import (
    C1R1RetrySuccessor,
    C1R1TemporalRetryRequired,
    TDG11C1R1CommittedRuntime,
    TDG11C1R1PreparedRuntime,
    commit_hlt17_admission,
    prepare_hlt17_admission_fresh,
    prepare_hlt17_admission_retry,
    require_hlt17_admission_family,
    require_hlt17_committed,
    require_hlt17_prepared,
    require_hlt17_retry_successor,
)
from .tdg6_temporal_admission_runtime import TDG6TemporalLedger
from . import tdg6_temporal_admission_runtime as tdg6
from .tdg7_binary64_subdivision_lattice import (
    TDG7Binary64SubdivisionPlan,
    validate_tdg7_binary64_subdivision_plan,
)


RUNTIME_PROTOCOL = "FGC-2-SF1-PROTO19"
IDENTITY_SCOPE_SYNTHETIC = "synthetic"
IDENTITY_SCOPE_LIVE = "live"
IDENTITY_SCOPES = (IDENTITY_SCOPE_SYNTHETIC, IDENTITY_SCOPE_LIVE)
METHOD_LABELS = ("RK4", "SSPRK3")
INTEGRATOR_BY_LABEL = {
    "RK4": PRIMARY_METHOD,
    "SSPRK3": COMPARATOR_METHOD,
}
SPATIAL_ORDER_BY_LABEL = {"RK4": 4, "SSPRK3": 2}
HLT17_LIVE_MEMBER_KEYS = (
    "RK4-2049",
    "RK4-4097",
    "RK4-8193",
    "SSPRK3-4097",
    "SSPRK3-8193",
    "SSPRK3-16385",
)
LIVE_MEMBER_PAIRS = (
    ("RK4", 2049),
    ("RK4", 4097),
    ("RK4", 8193),
    ("SSPRK3", 4097),
    ("SSPRK3", 8193),
    ("SSPRK3", 16385),
)
MIN_POINT_COUNT = 9
MAX_POINT_COUNT = 16385
CAMPAIGN_EXECUTION_AUTHORIZED = False
_SHA256 = re.compile(r"\A[0-9a-f]{64}\Z")
REFERENCE_KINDS = (
    "source_operator_closure",
    "origin_receipt",
)

# Sealed TDG6 helpers with no public equivalent.  Do not patch the owner.
_clone_tracers = tdg6._clone_tracers
_restore_tracers = tdg6._restore_tracers
_tracer_snapshot = tdg6._tracer_snapshot

PINNED_PRIVATE_DEPENDENCIES = (
    "recursive_horizons.fgc.evolution.tdg6_temporal_admission_runtime._clone_tracers",
    "recursive_horizons.fgc.evolution.tdg6_temporal_admission_runtime._restore_tracers",
    "recursive_horizons.fgc.evolution.tdg6_temporal_admission_runtime._tracer_snapshot",
)

HLT17_CURSOR_BRIDGE_SEAMS = (
    "predecessor_cursor_identity",
    "predecessor_cursor_identity_sha256",
    "predecessor_cursor_sha256",
    "predecessor_descriptor_sha256",
    "source_operator_closure",
    "origin_reference",
    "durable_store_bound",
    "accepted_generation_is_not_cursor_revision",
    "adopt_committed_does_not_call_kernel_commit",
)


class HLT17RuntimeMemberError(ValueError):
    """Malformed HLT17 identity, capture, restore, or adoption data."""


class HLT17RollbackError(RuntimeError):
    """A late failure could not restore the complete prior member boundary."""


def _fail(message: str) -> None:
    raise HLT17RuntimeMemberError(message)


def _text(value: object, label: str) -> str:
    if type(value) is not str or not value:
        raise TypeError(f"{label} must be nonempty text")
    return value


def _digest_text(value: object, label: str) -> str:
    text = _text(value, label)
    if _SHA256.fullmatch(text) is None:
        _fail(f"{label} must be a SHA-256 digest")
    return text


def _bool_false(value: object, label: str) -> bool:
    if type(value) is not bool or value is not False:
        _fail(f"{label} must be false")
    return False


def _nonnegative_int(value: object, label: str) -> int:
    if type(value) is not int or value < 0:
        raise TypeError(f"{label} must be a nonnegative built-in integer")
    return value


def _finite(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _same_bits(left: float, right: float) -> bool:
    return _finite("left", left).hex() == _finite("right", right).hex()


def _type_identity(value: object) -> str | None:
    if value is None:
        return None
    if callable(value) and hasattr(value, "__module__") and hasattr(value, "__qualname__"):
        return f"{value.__module__}.{value.__qualname__}"
    cls = type(value)
    return f"{cls.__module__}.{cls.__qualname__}"


def _hex_number(value: object, label: str) -> str:
    return _finite(label, value).hex()


def _canonical_mapping(value: Mapping[str, object], label: str) -> dict[str, object]:
    payload = dict(value) if type(value) is dict else {str(key): item for key, item in value.items()}
    try:
        canonical_json_bytes(payload)
    except (CanonicalJSONError, TypeError, ValueError) as error:
        raise HLT17RuntimeMemberError(f"{label} is not canonical JSON") from error
    return payload


def _normalized_static(value: object, label: str) -> object:
    if value is None or type(value) is bool or type(value) is str:
        return value
    if type(value) is int:
        return value
    if isinstance(value, bool) or not isinstance(value, Real):
        if hasattr(value, "__dict__") or hasattr(value, "__dataclass_fields__"):
            names = getattr(value, "__dataclass_fields__", None)
            if names is not None:
                return {
                    name: _normalized_static(getattr(value, name), f"{label}.{name}")
                    for name in names
                }
            raise HLT17RuntimeMemberError(f"unsupported static identity type at {label}")
        raise HLT17RuntimeMemberError(f"unsupported static identity type at {label}")
    number = _finite(label, value)
    if type(value) is int:
        return int(value)
    return {"binary64_hex": number.hex()}


def _operator_identity(operator: object) -> dict[str, object]:
    fields: dict[str, object] = {"type": _type_identity(operator)}
    for name in (
        "ko_dissipation",
        "residual_tolerance",
        "raw_tolerance",
        "kinetic_condition_maximum",
    ):
        if hasattr(operator, name):
            fields[name] = _normalized_static(getattr(operator, name), name)
    derivative = getattr(operator, "derivative", None)
    if derivative is not None and hasattr(derivative, "order"):
        fields["derivative_order"] = _nonnegative_int(derivative.order, "derivative order")
    return fields


def _state_hash(state: EvolutionState) -> str:
    if not isinstance(state, EvolutionState) or state.shape[1] != 6:
        _fail("HLT17 state must contain all six frozen fields")
    return array_content_sha256(state.u, state.p, state.q)


def _grid(initial: object) -> object:
    grid = getattr(initial, "grid", None)
    if grid is None:
        raise TypeError(
            "HLT17 initial must provide a grid; the carrier does not reconstruct one"
        )
    for name in ("minimum", "maximum", "spacing", "point_count", "coordinates"):
        if not hasattr(grid, name):
            raise TypeError(f"HLT17 grid omits {name}")
    return grid


def _coordinates_from_grid(grid: object, point_count: int) -> np.ndarray:
    raw = np.ascontiguousarray(np.asarray(grid.coordinates, dtype=np.float64))
    if (
        raw.ndim != 1
        or raw.size != point_count
        or raw.size < MIN_POINT_COUNT
        or raw.size > MAX_POINT_COUNT
        or not np.all(np.isfinite(raw))
        or not np.all(raw[1:] > raw[:-1])
    ):
        _fail("HLT17 grid coordinates are not a strictly ordered finite identity")
    count = getattr(grid, "point_count")
    if type(count) is not int or count != point_count or count != raw.size:
        _fail("HLT17 grid point count differs from the member identity")
    copied = np.array(raw, dtype=np.dtype("<f8"), order="C", copy=True)
    copied.setflags(write=False)
    return copied


def _copy_array(value: object, label: str) -> np.ndarray:
    array = np.ascontiguousarray(np.asarray(value, dtype=np.float64))
    if array.size == 0 or not np.all(np.isfinite(array)):
        _fail(f"{label} is not a nonempty finite array")
    copied = np.array(array, dtype=np.dtype("<f8"), order="C", copy=True)
    copied.setflags(write=False)
    return copied


@dataclass(frozen=True, slots=True)
class HLT17ExternalReference:
    """A supplied digest.  Matching it does not authenticate the origin."""

    kind: str
    sha256: str
    authenticated: bool = False

    def __post_init__(self) -> None:
        kind = _text(self.kind, "reference kind")
        if kind not in REFERENCE_KINDS:
            _fail("unknown HLT17 external reference kind")
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "sha256", _digest_text(self.sha256, "reference sha256"))
        object.__setattr__(
            self, "authenticated", _bool_false(self.authenticated, "authenticated")
        )

    def as_mapping(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "sha256": self.sha256,
            "authenticated": False,
        }


@dataclass(frozen=True, slots=True)
class HLT17MemberIdentity:
    """Declared carrier identity.  Live keys are never inferred from array size."""

    scope: str
    method_label: str
    point_count: int
    integrator_id: str
    spatial_order: int
    campaign_id: str
    amplitude: str

    def __post_init__(self) -> None:
        scope = _text(self.scope, "scope")
        if scope not in IDENTITY_SCOPES:
            _fail("HLT17 identity scope differs")
        method = _text(self.method_label, "method_label")
        if method not in METHOD_LABELS:
            _fail("HLT17 method label differs")
        count = _nonnegative_int(self.point_count, "point_count")
        if count < MIN_POINT_COUNT or count > MAX_POINT_COUNT:
            _fail("HLT17 point count is outside the bounded identity range")
        integrator = _text(self.integrator_id, "integrator_id")
        if integrator != INTEGRATOR_BY_LABEL[method]:
            _fail("HLT17 integrator identity disagrees with the method label")
        order = _nonnegative_int(self.spatial_order, "spatial_order")
        if order != SPATIAL_ORDER_BY_LABEL[method]:
            _fail("HLT17 spatial order disagrees with the method label")
        campaign = _text(self.campaign_id, "campaign_id")
        amplitude = _text(self.amplitude, "amplitude")
        key = f"{method}-{count}"
        live_pair = (method, count)
        if scope == IDENTITY_SCOPE_LIVE:
            if key not in HLT17_LIVE_MEMBER_KEYS or live_pair not in LIVE_MEMBER_PAIRS:
                _fail("live HLT17 identity is outside the future production set")
            if amplitude != "3":
                _fail("live HLT17 amplitude differs")
        else:
            if key in HLT17_LIVE_MEMBER_KEYS or live_pair in LIVE_MEMBER_PAIRS:
                _fail("synthetic HLT17 identity cannot use a live production key")
            if amplitude != "synthetic":
                _fail("synthetic HLT17 amplitude differs")
        object.__setattr__(self, "scope", scope)
        object.__setattr__(self, "method_label", method)
        object.__setattr__(self, "point_count", count)
        object.__setattr__(self, "integrator_id", integrator)
        object.__setattr__(self, "spatial_order", order)
        object.__setattr__(self, "campaign_id", campaign)
        object.__setattr__(self, "amplitude", amplitude)

    @property
    def member_key(self) -> str:
        return f"{self.method_label}-{self.point_count}"

    def as_mapping(self) -> dict[str, object]:
        return {
            "scope": self.scope,
            "method_label": self.method_label,
            "point_count": self.point_count,
            "integrator_id": self.integrator_id,
            "spatial_order": self.spatial_order,
            "campaign_id": self.campaign_id,
            "amplitude": self.amplitude,
            "member_key": self.member_key,
        }


@dataclass(eq=False)
class HLT17MemberCapture:
    """Exact dynamic boundary copy.  Static source/grid objects stay on the member."""

    state: EvolutionState
    time: float
    step_index: int
    transaction_serial: int
    source_retry_count: int
    CFL_retry_count: int
    monitor_state: GR0RuntimeMonitorState
    causal_state: CausalBudgetState
    tracer_labels: np.ndarray
    tracer_positions: np.ndarray
    tracer_proper_times: np.ndarray
    event_proper_times: tuple[np.ndarray, ...]
    event_fields: tuple[np.ndarray, ...]
    tracer_cutoff: float
    tracer_outer_radius: float
    temporal_ledger: TDG11IMP1Ledger
    executed_plan: TDG7Binary64SubdivisionPlan | None
    coordinates_sha256: str
    inherited_snapshot: str
    operator: object
    projector: object | None
    initial: object


@dataclass(slots=True)
class HLT17RuntimeMember:
    """In-memory HLT17 member carrier with capture, restore, and C1R1 adoption."""

    identity: HLT17MemberIdentity
    initial: object
    operator: Callable[..., Any]
    projector: Callable[..., Any] | None
    transaction: GR0RuntimeStageTransaction
    tracers: object
    state: EvolutionState
    temporal_ledger: TDG11IMP1Ledger
    source_operator_closure: HLT17ExternalReference
    origin_reference: HLT17ExternalReference
    time: float
    step_index: int = 0
    transaction_serial: int = 0
    CFL_retry_count: int = 0
    source_retry_count: int = 0
    executed_plan: TDG7Binary64SubdivisionPlan | None = None
    _coordinates: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if type(self) is not HLT17RuntimeMember:
            raise TypeError("HLT17RuntimeMember must not be subclassed as a historical record")
        if not isinstance(self.identity, HLT17MemberIdentity):
            raise TypeError("identity must be HLT17MemberIdentity")
        if not callable(self.operator):
            raise TypeError("operator must be the provided source callable")
        if self.projector is not None and not callable(self.projector):
            raise TypeError("projector must be callable or None")
        if not isinstance(self.transaction, GR0RuntimeStageTransaction):
            raise TypeError("transaction must be GR0RuntimeStageTransaction")
        if not isinstance(self.state, EvolutionState):
            raise TypeError("state must be EvolutionState")
        if not isinstance(self.temporal_ledger, TDG11IMP1Ledger):
            raise TypeError("temporal_ledger must be TDG11IMP1Ledger")
        if isinstance(self.temporal_ledger, TDG6TemporalLedger):
            raise TypeError("HLT17 cannot lawfully carry a TDG6 ledger as its live ledger")
        if not isinstance(self.source_operator_closure, HLT17ExternalReference):
            raise TypeError("source_operator_closure must be HLT17ExternalReference")
        if not isinstance(self.origin_reference, HLT17ExternalReference):
            raise TypeError("origin_reference must be HLT17ExternalReference")
        if self.source_operator_closure.kind != "source_operator_closure":
            _fail("source_operator_closure kind differs")
        if self.origin_reference.kind != "origin_receipt":
            _fail("origin_reference kind differs")
        if self.origin_reference.authenticated is not False:
            _fail("origin construction does not authenticate the external origin")
        if self.source_operator_closure.authenticated is not False:
            _fail("closure construction does not authenticate the external origin")
        if self.temporal_ledger.historical_origin_authenticated is not False:
            _fail("IMP1 origin authentication is not claimed by carrier construction")
        if CAMPAIGN_EXECUTION_AUTHORIZED is not False:
            _fail("HLT17 carrier is not a production execution authority")
        self.time = _finite("member time", self.time)
        self.step_index = _nonnegative_int(self.step_index, "step_index")
        self.transaction_serial = _nonnegative_int(
            self.transaction_serial, "transaction_serial"
        )
        self.CFL_retry_count = _nonnegative_int(self.CFL_retry_count, "CFL_retry_count")
        self.source_retry_count = _nonnegative_int(
            self.source_retry_count, "source_retry_count"
        )
        if self.state.shape[0] != self.identity.point_count or self.state.shape[1] != 6:
            _fail("HLT17 state shape disagrees with the declared identity")
        grid = _grid(self.initial)
        coordinates = _coordinates_from_grid(grid, self.identity.point_count)
        try:
            existing = object.__getattribute__(self, "_coordinates")
        except AttributeError:
            existing = None
        if existing is None:
            object.__setattr__(self, "_coordinates", coordinates)
        elif array_content_sha256(existing) != array_content_sha256(coordinates):
            _fail("provided grid coordinates became stale")
        if self.temporal_ledger.method != self.identity.integrator_id:
            _fail("IMP1 ledger method differs from the member integrator")
        if not _same_bits(self.temporal_ledger.last_accepted_time, self.time):
            _fail("IMP1 ledger time differs from member time")
        if self.temporal_ledger.current_state_sha256 != _state_hash(self.state):
            _fail("IMP1 ledger state hash differs from member state")
        if self.temporal_ledger.current_step_index != self.step_index:
            _fail("IMP1 ledger step index differs from member step index")
        if self.temporal_ledger.current_transaction_serial != self.transaction_serial:
            _fail("IMP1 ledger serial differs from member serial")
        if self.origin_reference.sha256 != self.temporal_ledger.origin_receipt_sha256:
            _fail("origin reference digest differs from the IMP1 origin receipt reference")
        monitor = self.transaction.state
        if not isinstance(monitor, GR0RuntimeMonitorState):
            raise TypeError("monitor must be GR0RuntimeMonitorState")
        count = monitor.accepted_stage_count
        if type(count) is not int or count < 0:
            _fail("monitor accepted_stage_count differs")
        serial = monitor.last_transaction_serial
        if type(serial) is not int or serial < -1:
            _fail("monitor last_transaction_serial differs")
        if count != self.transaction_serial or count != serial + 1:
            _fail("monitor count/serial disagree with the member transaction serial")
        if count > 0 and not _same_bits(monitor.last_accepted_time, self.time):
            _fail("nonempty monitor accepted time differs from member time")
        causal = self.transaction.causal_state
        if not isinstance(causal, CausalBudgetState):
            raise TypeError("causal_state must be CausalBudgetState")
        if not _same_bits(causal.accepted_time, self.time):
            _fail("causal accepted time differs from member time")
        _tracer_snapshot(self.tracers)
        if self.executed_plan is not None:
            validate_tdg7_binary64_subdivision_plan(self.executed_plan)
        inherited_tdg6_ledger(self.temporal_ledger)

    @property
    def key(self) -> str:
        return self.identity.member_key

    @property
    def method_label(self) -> str:
        return self.identity.method_label

    @property
    def point_count(self) -> int:
        return self.identity.point_count

    @property
    def integrator_id(self) -> str:
        return self.identity.integrator_id

    @property
    def coordinates(self) -> np.ndarray:
        return self._coordinates

    @property
    def coordinates_sha256(self) -> str:
        return array_content_sha256(self._coordinates)

    def inherited_tdg6_sibling(self) -> TDG6TemporalLedger:
        """Frozen inherited TDG6 snapshot.  It is never the live IMP1 ledger."""

        return inherited_tdg6_ledger(self.temporal_ledger)

    def runtime_template(self) -> dict[str, object]:
        """Static shell identity.  Digest agreement is not origin authentication."""

        grid = _grid(self.initial)
        transaction = self.transaction
        template = {
            "scope": self.identity.scope,
            "amplitude": self.identity.amplitude,
            "method": self.identity.method_label,
            "point_count": self.identity.point_count,
            "integrator_id": self.identity.integrator_id,
            "spatial_order": self.identity.spatial_order,
            "member_type": _type_identity(self),
            "initial_type": _type_identity(self.initial),
            "operator": _operator_identity(self.operator),
            "projector": {"callable": _type_identity(self.projector)},
            "transaction_type": _type_identity(transaction),
            "tracer_type": _type_identity(self.tracers),
            "grid": {
                "minimum_hex": _hex_number(grid.minimum, "grid minimum"),
                "maximum_hex": _hex_number(grid.maximum, "grid maximum"),
                "spacing_hex": _hex_number(grid.spacing, "grid spacing"),
                "point_count": self.identity.point_count,
                "coordinates_sha256": self.coordinates_sha256,
            },
            "transaction": {
                "thresholds": _normalized_static(transaction.thresholds, "thresholds"),
                "boundary_geometry": _normalized_static(
                    transaction.boundary_geometry, "boundary_geometry"
                ),
                "grid_spacing_hex": _hex_number(transaction.grid_spacing, "grid_spacing"),
                "cfl_maximum_hex": _hex_number(transaction.cfl_maximum, "cfl_maximum"),
                "hat_normal_factor_hex": _hex_number(
                    transaction.hat_normal_factor, "hat_normal_factor"
                ),
            },
        }
        _canonical_mapping(template, "runtime template")
        return template

    def runtime_template_sha256(self) -> str:
        return sha256(canonical_json_bytes(self.runtime_template())).hexdigest()

    def physical_state_sha256(self) -> str:
        return _state_hash(self.state)

    def capture(self) -> HLT17MemberCapture:
        """Copy every dynamic field, including tracer histories and the IMP1 ledger."""

        self.__post_init__()
        cloned = _clone_tracers(self.tracers)
        snapshot = _tracer_snapshot(cloned)
        if snapshot != _tracer_snapshot(self.tracers):
            _fail("tracer capture is not a complete identity copy")
        plan = self.executed_plan
        if plan is not None:
            validate_tdg7_binary64_subdivision_plan(plan)
        return HLT17MemberCapture(
            state=self.state,
            time=self.time,
            step_index=self.step_index,
            transaction_serial=self.transaction_serial,
            source_retry_count=self.source_retry_count,
            CFL_retry_count=self.CFL_retry_count,
            monitor_state=self.transaction.state,
            causal_state=self.transaction.causal_state,
            tracer_labels=_copy_array(cloned.labels, "tracer labels"),
            tracer_positions=_copy_array(cloned.positions, "tracer positions"),
            tracer_proper_times=_copy_array(cloned.proper_times, "tracer proper times"),
            event_proper_times=tuple(
                _copy_array(item, "event proper times") for item in cloned.event_proper_times
            ),
            event_fields=tuple(
                _copy_array(item, "event fields") for item in cloned.event_fields
            ),
            tracer_cutoff=_finite("tracer cutoff", cloned.cutoff),
            tracer_outer_radius=_finite("tracer outer radius", cloned.outer_radius),
            temporal_ledger=self.temporal_ledger,
            executed_plan=plan,
            coordinates_sha256=self.coordinates_sha256,
            inherited_snapshot=self.temporal_ledger.inherited_snapshot,
            operator=self.operator,
            projector=self.projector,
            initial=self.initial,
        )

    def _prevalidate_capture(self, snapshot: object) -> HLT17MemberCapture:
        if not isinstance(snapshot, HLT17MemberCapture):
            raise TypeError("snapshot must be HLT17MemberCapture")
        if snapshot.operator is not self.operator:
            _fail("capture/restore keeps the original source identity")
        if snapshot.projector is not self.projector:
            _fail("capture/restore keeps the original projector identity")
        if snapshot.initial is not self.initial:
            _fail("capture/restore keeps the original initial/grid identity")
        if snapshot.coordinates_sha256 != self.coordinates_sha256:
            _fail("capture grid identity differs from the live member")
        live_grid_hash = array_content_sha256(
            np.ascontiguousarray(
                np.asarray(_grid(self.initial).coordinates, dtype=np.float64)
            )
        )
        if live_grid_hash != self.coordinates_sha256:
            _fail("provided grid coordinates became stale")
        if not isinstance(snapshot.state, EvolutionState):
            raise TypeError("capture state must be EvolutionState")
        if snapshot.state.shape != (self.identity.point_count, 6):
            _fail("capture state shape differs")
        time = _finite("capture time", snapshot.time)
        step = _nonnegative_int(snapshot.step_index, "capture step_index")
        serial = _nonnegative_int(
            snapshot.transaction_serial, "capture transaction_serial"
        )
        source_retries = _nonnegative_int(
            snapshot.source_retry_count, "capture source_retry_count"
        )
        cfl_retries = _nonnegative_int(
            snapshot.CFL_retry_count, "capture CFL_retry_count"
        )
        if not isinstance(snapshot.monitor_state, GR0RuntimeMonitorState):
            raise TypeError("capture monitor must be GR0RuntimeMonitorState")
        if not isinstance(snapshot.causal_state, CausalBudgetState):
            raise TypeError("capture causal state must be CausalBudgetState")
        if not isinstance(snapshot.temporal_ledger, TDG11IMP1Ledger):
            raise TypeError("capture ledger must be TDG11IMP1Ledger")
        if snapshot.inherited_snapshot != self.temporal_ledger.inherited_snapshot:
            _fail("restore/adoption cannot replace the immutable inherited TDG6 sibling")
        if snapshot.temporal_ledger.inherited_snapshot != snapshot.inherited_snapshot:
            _fail("capture inherited snapshot disagrees with its IMP1 ledger")
        if snapshot.temporal_ledger.method != self.identity.integrator_id:
            _fail("capture ledger method differs")
        if not _same_bits(snapshot.temporal_ledger.last_accepted_time, time):
            _fail("capture ledger time differs")
        if snapshot.temporal_ledger.current_state_sha256 != _state_hash(snapshot.state):
            _fail("capture ledger state hash differs")
        if snapshot.temporal_ledger.current_step_index != step:
            _fail("capture ledger step differs")
        if snapshot.temporal_ledger.current_transaction_serial != serial:
            _fail("capture ledger serial differs")
        plan = snapshot.executed_plan
        if plan is not None:
            validate_tdg7_binary64_subdivision_plan(plan)
        labels = _copy_array(snapshot.tracer_labels, "tracer labels")
        positions = _copy_array(snapshot.tracer_positions, "tracer positions")
        proper = _copy_array(snapshot.tracer_proper_times, "tracer proper times")
        if labels.shape != positions.shape or labels.shape != proper.shape:
            _fail("capture tracer shapes differ")
        events_t = tuple(
            _copy_array(item, "event proper times") for item in snapshot.event_proper_times
        )
        events_f = tuple(_copy_array(item, "event fields") for item in snapshot.event_fields)
        if len(events_t) != len(events_f) or not events_t:
            _fail("capture tracer event history is empty or inconsistent")
        cutoff = _finite("tracer cutoff", snapshot.tracer_cutoff)
        outer = _finite("tracer outer radius", snapshot.tracer_outer_radius)
        if not (cutoff > 0.0 and outer > 0.0):
            _fail("tracer cutoff and outer radius must be positive")
        return HLT17MemberCapture(
            state=snapshot.state,
            time=time,
            step_index=step,
            transaction_serial=serial,
            source_retry_count=source_retries,
            CFL_retry_count=cfl_retries,
            monitor_state=snapshot.monitor_state,
            causal_state=snapshot.causal_state,
            tracer_labels=labels,
            tracer_positions=positions,
            tracer_proper_times=proper,
            event_proper_times=events_t,
            event_fields=events_f,
            tracer_cutoff=cutoff,
            tracer_outer_radius=outer,
            temporal_ledger=snapshot.temporal_ledger,
            executed_plan=plan,
            coordinates_sha256=snapshot.coordinates_sha256,
            inherited_snapshot=snapshot.inherited_snapshot,
            operator=snapshot.operator,
            projector=snapshot.projector,
            initial=snapshot.initial,
        )

    def _apply_capture(self, snapshot: HLT17MemberCapture) -> None:
        self.state = snapshot.state
        self.time = snapshot.time
        self.step_index = snapshot.step_index
        self.transaction_serial = snapshot.transaction_serial
        self.source_retry_count = snapshot.source_retry_count
        self.CFL_retry_count = snapshot.CFL_retry_count
        self.transaction.state = snapshot.monitor_state
        self.transaction.causal_state = snapshot.causal_state
        self.tracers.labels = np.array(snapshot.tracer_labels, dtype=np.float64, copy=True)
        self.tracers.positions = np.array(
            snapshot.tracer_positions, dtype=np.float64, copy=True
        )
        self.tracers.proper_times = np.array(
            snapshot.tracer_proper_times, dtype=np.float64, copy=True
        )
        self.tracers.event_proper_times = [
            np.array(item, dtype=np.float64, copy=True)
            for item in snapshot.event_proper_times
        ]
        self.tracers.event_fields = [
            np.array(item, dtype=np.float64, copy=True) for item in snapshot.event_fields
        ]
        self.tracers.cutoff = snapshot.tracer_cutoff
        self.tracers.outer_radius = snapshot.tracer_outer_radius
        self.temporal_ledger = snapshot.temporal_ledger
        self.executed_plan = snapshot.executed_plan

    def restore(
        self,
        snapshot: HLT17MemberCapture,
        *,
        fault_hook: Callable[[str], None] | None = None,
    ) -> None:
        """Restore only after complete prevalidation; late failure rolls back."""

        checked = self._prevalidate_capture(snapshot)
        prior = self.capture()
        try:
            self._apply_capture(checked)
            if fault_hook is not None:
                fault_hook("after_restore_mutation")
            self.__post_init__()
            if self.temporal_ledger.inherited_snapshot != prior.inherited_snapshot:
                _fail("restore rewrote the inherited TDG6 sibling")
            if self.operator is not prior.operator or self.projector is not prior.projector:
                _fail("restore replaced source or projector identity")
        except BaseException as error:
            try:
                self._apply_capture(prior)
                self.__post_init__()
            except BaseException as rollback_error:
                raise HLT17RollbackError(
                    "HLT17 complete restore rollback failed"
                ) from rollback_error
            raise error

    def prepare_fresh(
        self,
        *,
        event_target: Real,
        requested_cap: Real,
        limits: IMP1EnclosureLimits = IMP1EnclosureLimits(),
    ) -> TDG11C1R1PreparedRuntime:
        """Plan the C1R1 lattice before any source/projector hook."""

        self.__post_init__()
        return prepare_hlt17_admission_fresh(
            method=self.identity.integrator_id,
            time=self.time,
            event_target=event_target,
            requested_cap=requested_cap,
            state=self.state,
            rhs=self.operator,
            projector=self.projector,
            transaction=self.transaction,
            tracers=self.tracers,
            coordinates=self.coordinates,
            temporal_ledger=self.temporal_ledger,
            previous_step_index=self.step_index,
            previous_transaction_serial=self.transaction_serial,
            limits=limits,
        )

    def require_admission(
        self,
        prepared: TDG11C1R1PreparedRuntime,
        *,
        durable_rejection_sink: Callable[[Mapping[str, object]], object] | None = None,
    ) -> None:
        self.__post_init__()
        bound = require_hlt17_prepared(prepared)
        if bound.original_rhs is not self.operator:
            _fail("prepared source identity differs from the live member")
        if bound.original_projector is not self.projector:
            _fail("prepared projector identity differs from the live member")
        require_hlt17_admission_family(
            bound,
            transaction=self.transaction,
            tracers=self.tracers,
            temporal_ledger=self.temporal_ledger,
            current_time=self.time,
            current_state=self.state,
            current_step_index=self.step_index,
            current_transaction_serial=self.transaction_serial,
            coordinates=self.coordinates,
            durable_rejection_sink=durable_rejection_sink,
        )

    def adopt_retry_ledger(
        self,
        retry: C1R1TemporalRetryRequired,
        *,
        fault_hook: Callable[[str], None] | None = None,
    ) -> None:
        """Install the durably acknowledged retry ledger without advancing state."""

        if not isinstance(retry, C1R1TemporalRetryRequired):
            raise TypeError("retry must be C1R1TemporalRetryRequired")
        self.__post_init__()
        updated = retry.updated_ledger
        if not isinstance(updated, TDG11IMP1Ledger):
            raise TypeError("retry ledger must be TDG11IMP1Ledger")
        if updated.inherited_snapshot != self.temporal_ledger.inherited_snapshot:
            _fail("retry adoption cannot rewrite the inherited TDG6 sibling")
        if updated.current_state_sha256 != self.physical_state_sha256():
            _fail("retry ledger cannot replace the accepted physical state")
        if not _same_bits(updated.last_accepted_time, self.time):
            _fail("retry ledger cannot advance accepted time")
        if updated.current_step_index != self.step_index:
            _fail("retry ledger cannot advance the step index")
        if updated.current_transaction_serial != self.transaction_serial:
            _fail("retry ledger cannot advance the transaction serial")
        if updated.current_macro_step_temporal_retry_count < 1:
            _fail("retry ledger does not record a temporal retry")
        prior = self.capture()
        try:
            self.temporal_ledger = updated
            if fault_hook is not None:
                fault_hook("after_retry_install")
            self.__post_init__()
            if _state_hash(self.state) != _state_hash(prior.state):
                _fail("retry adoption mutated physical state")
            if self.temporal_ledger.inherited_snapshot != prior.inherited_snapshot:
                _fail("retry adoption rewrote the inherited TDG6 sibling")
        except BaseException as error:
            try:
                self._apply_capture(prior)
                self.__post_init__()
            except BaseException as rollback_error:
                raise HLT17RollbackError(
                    "HLT17 complete retry-ledger rollback failed"
                ) from rollback_error
            raise error

    def prepare_retry(self, successor: C1R1RetrySuccessor) -> TDG11C1R1PreparedRuntime:
        self.__post_init__()
        checked = require_hlt17_retry_successor(successor)
        if checked.updated_ledger is not self.temporal_ledger:
            if (
                checked.updated_ledger.inherited_snapshot
                != self.temporal_ledger.inherited_snapshot
                or checked.updated_ledger.current_state_sha256
                != self.physical_state_sha256()
            ):
                _fail("retry successor ledger does not match the live member")
        return prepare_hlt17_admission_retry(
            checked,
            state=self.state,
            rhs=self.operator,
            projector=self.projector,
            transaction=self.transaction,
            tracers=self.tracers,
            coordinates=self.coordinates,
        )

    def commit_prepared(
        self,
        prepared: TDG11C1R1PreparedRuntime,
        *,
        fault_hook: Callable[[str], None] | None = None,
    ) -> TDG11C1R1CommittedRuntime:
        """Adopt the C1R1 fine endpoint onto this carrier with complete rollback."""

        self.__post_init__()
        bound = require_hlt17_prepared(prepared)
        if bound.original_rhs is not self.operator:
            _fail("commit keeps the original source identity")
        if bound.original_projector is not self.projector:
            _fail("commit keeps the original projector identity")
        if bound.method != self.identity.integrator_id:
            _fail("prepared method differs from the live member")
        prior = self.capture()
        try:
            committed = commit_hlt17_admission(
                bound,
                transaction=self.transaction,
                tracers=self.tracers,
                temporal_ledger=self.temporal_ledger,
                current_time=self.time,
                current_state=self.state,
                current_step_index=self.step_index,
                current_transaction_serial=self.transaction_serial,
                coordinates=self.coordinates,
            )
            if fault_hook is not None:
                fault_hook("after_c1r1_commit")
            committed = require_hlt17_committed(committed)
            if committed.temporal_ledger.inherited_snapshot != prior.inherited_snapshot:
                _fail("C1R1 commit rewrote the inherited TDG6 sibling")
            self.state = committed.state
            self.time = committed.time
            self.step_index = committed.step_index
            self.transaction_serial = committed.transaction_serial
            self.temporal_ledger = committed.temporal_ledger
            self.executed_plan = bound.plan
            if fault_hook is not None:
                fault_hook("after_member_install")
            self.__post_init__()
            if self.operator is not prior.operator or self.projector is not prior.projector:
                _fail("commit replaced source or projector identity")
            return committed
        except BaseException as error:
            try:
                self._apply_capture(prior)
                self.__post_init__()
            except BaseException as rollback_error:
                raise HLT17RollbackError(
                    "HLT17 complete adoption rollback failed"
                ) from rollback_error
            raise error

    def adopt_committed(
        self,
        committed: TDG11C1R1CommittedRuntime,
        *,
        executed_plan: TDG7Binary64SubdivisionPlan,
        source_transaction: GR0RuntimeStageTransaction | None = None,
        source_tracers: object | None = None,
        source_retry_count: int | None = None,
        CFL_retry_count: int | None = None,
        rollback_to: HLT17MemberCapture | None = None,
        fault_hook: Callable[[str], None] | None = None,
    ) -> TDG11C1R1CommittedRuntime:
        """Install an already-committed C1R1 endpoint. Never commits the kernel again.

        The kernel may already have advanced this carrier's transaction/tracers.
        Pass ``rollback_to`` captured before that commit; do not revalidate the
        stale member serial against the advanced monitor first.
        """

        committed = require_hlt17_committed(committed)
        if committed.method != self.identity.integrator_id:
            _fail("committed method differs from the live member")
        if committed.campaign_execution_authorized is not False:
            _fail("committed record cannot authorize a campaign")
        if not isinstance(executed_plan, TDG7Binary64SubdivisionPlan):
            raise TypeError("executed_plan must be TDG7Binary64SubdivisionPlan")
        validate_tdg7_binary64_subdivision_plan(executed_plan)
        if not _same_bits(executed_plan.boundaries[4], committed.time):
            _fail("executed plan endpoint is not the committed time")
        if committed.temporal_ledger.inherited_snapshot != self.temporal_ledger.inherited_snapshot:
            _fail("committed ledger cannot rewrite the inherited TDG6 sibling")
        if source_transaction is None:
            source_transaction = self.transaction
        if source_tracers is None:
            source_tracers = self.tracers
        if not isinstance(source_transaction, GR0RuntimeStageTransaction):
            raise TypeError("source_transaction must be GR0RuntimeStageTransaction")
        next_source = (
            self.source_retry_count
            if source_retry_count is None
            else _nonnegative_int(source_retry_count, "source_retry_count")
        )
        next_cfl = (
            self.CFL_retry_count
            if CFL_retry_count is None
            else _nonnegative_int(CFL_retry_count, "CFL_retry_count")
        )
        prior = rollback_to if rollback_to is not None else self.capture()
        if not isinstance(prior, HLT17MemberCapture):
            raise TypeError("rollback_to must be HLT17MemberCapture")
        try:
            if source_transaction is not self.transaction:
                self.transaction.state = source_transaction.state
                self.transaction.causal_state = source_transaction.causal_state
            if source_tracers is not self.tracers:
                _restore_tracers(self.tracers, source_tracers)
            self.state = committed.state
            self.time = committed.time
            self.step_index = committed.step_index
            self.transaction_serial = committed.transaction_serial
            self.temporal_ledger = committed.temporal_ledger
            self.executed_plan = executed_plan
            self.source_retry_count = next_source
            self.CFL_retry_count = next_cfl
            if fault_hook is not None:
                fault_hook("after_committed_install")
            self.__post_init__()
            if self.operator is not prior.operator or self.projector is not prior.projector:
                _fail("committed adoption replaced source or projector identity")
            if self.temporal_ledger.inherited_snapshot != prior.inherited_snapshot:
                _fail("committed adoption rewrote the inherited TDG6 sibling")
            if _state_hash(self.state) != committed.temporal_ledger.current_state_sha256:
                _fail("committed adoption disagrees with the committed ledger")
            return committed
        except BaseException as error:
            try:
                self._apply_capture(prior)
                self.__post_init__()
            except BaseException as rollback_error:
                raise HLT17RollbackError(
                    "HLT17 complete committed-adoption rollback failed"
                ) from rollback_error
            raise error

    def install_retry_overlay_state(
        self,
        ledger: TDG11IMP1Ledger,
        *,
        source_retry_count: int,
        CFL_retry_count: int,
        executed_plan: TDG7Binary64SubdivisionPlan | None,
        fault_hook: Callable[[str], None] | None = None,
    ) -> None:
        """Install a checked retry ledger and totals without advancing state or time.

        ``executed_plan`` is the last accepted plan. Retry-overlay encodings omit
        it; this installer must not silently erase it.
        """

        if not isinstance(ledger, TDG11IMP1Ledger):
            raise TypeError("ledger must be TDG11IMP1Ledger")
        self.__post_init__()
        if ledger.inherited_snapshot != self.temporal_ledger.inherited_snapshot:
            _fail("retry overlay cannot rewrite the inherited TDG6 sibling")
        if ledger.current_state_sha256 != self.physical_state_sha256():
            _fail("retry overlay cannot replace the accepted physical state")
        if not _same_bits(ledger.last_accepted_time, self.time):
            _fail("retry overlay cannot advance accepted time")
        if ledger.current_step_index != self.step_index:
            _fail("retry overlay cannot advance the step index")
        if ledger.current_transaction_serial != self.transaction_serial:
            _fail("retry overlay cannot advance the transaction serial")
        if executed_plan is not None:
            validate_tdg7_binary64_subdivision_plan(executed_plan)
        next_source = _nonnegative_int(source_retry_count, "source_retry_count")
        next_cfl = _nonnegative_int(CFL_retry_count, "CFL_retry_count")
        prior = self.capture()
        try:
            self.temporal_ledger = ledger
            self.source_retry_count = next_source
            self.CFL_retry_count = next_cfl
            self.executed_plan = executed_plan
            if fault_hook is not None:
                fault_hook("after_retry_overlay_install")
            self.__post_init__()
            if _state_hash(self.state) != _state_hash(prior.state):
                _fail("retry overlay mutated physical state")
            if not _same_bits(self.time, prior.time):
                _fail("retry overlay mutated accepted time")
            if self.temporal_ledger.inherited_snapshot != prior.inherited_snapshot:
                _fail("retry overlay rewrote the inherited TDG6 sibling")
            if executed_plan is not None and self.executed_plan != executed_plan:
                _fail("retry overlay dropped the last accepted executed plan")
        except BaseException as error:
            try:
                self._apply_capture(prior)
                self.__post_init__()
            except BaseException as rollback_error:
                raise HLT17RollbackError(
                    "HLT17 complete retry-overlay rollback failed"
                ) from rollback_error
            raise error


__all__ = [
    "CAMPAIGN_EXECUTION_AUTHORIZED",
    "HLT17_CURSOR_BRIDGE_SEAMS",
    "HLT17_LIVE_MEMBER_KEYS",
    "HLT17ExternalReference",
    "HLT17MemberCapture",
    "HLT17MemberIdentity",
    "HLT17RollbackError",
    "HLT17RuntimeMember",
    "HLT17RuntimeMemberError",
    "IDENTITY_SCOPE_LIVE",
    "IDENTITY_SCOPE_SYNTHETIC",
    "PINNED_PRIVATE_DEPENDENCIES",
    "RUNTIME_PROTOCOL",
]
