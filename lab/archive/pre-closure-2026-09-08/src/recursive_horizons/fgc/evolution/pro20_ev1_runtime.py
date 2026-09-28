"""One-shot PRO20-EV1 first-event runner core.

This is an implementation parent, not an execution freeze.  It reconstructs
the authenticated HLT15 t=23/16 origin through ``pro20_origin`` and
``pro20_source_factory``, seeds the production store atomically, acquires
one writer, and advances the canonical unfinished or pending member one
attempt at a time to t=3/2.  Runner labels are not evidence.  Outcomes are
mapped from immutable bundle bytes, cursor ownership, and typed stops.

There is no resume or takeover path.  A visible production namespace
consumes the one-shot.  Temporal exhaustion publishes the unchanged accepted
bundle plus non-executable terminal evidence.  Typed resource or premise
stops after seed publish immutable terminal evidence.  Unexpected,
programming, or post-publication uncertainty abandons the writer and leaves
the durable session unclosed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
import math
import os
from pathlib import Path
import resource
import sys
import time
from types import MappingProxyType
from typing import Callable, Mapping

from recursive_horizons.evidence_io import load_canonical_json, read_regular_file

from .hlt17_admission_runtime import (
    C1R1CoordinateLatticeStop,
    C1R1RefinementPathStop,
    C1R1RuntimeResourceStop,
    C1R1TemporalRetryExhausted,
    hlt17_implementation_identity,
)
from .hlt17_imp1_cursor import (
    CFL_OWNER,
    FRESH_READY,
    HLT17_CFL_RETRY_CAP,
    HLT17InvalidPremise,
    HLT17OwnerRetryExhausted,
    HLT17_SOURCE_RETRY_CAP,
    HLT17IMP1Cursor,
    RETRY_PENDING,
    SOURCE_OWNER,
    ledger_sha256,
    revalidate_hlt17_cursor,
    restore_hlt17_cursor,
)
from .hlt17_member_checkpoint import (
    ACCEPTED_FINE,
    HLT17BridgeResult,
    HLT17CheckpointResult,
    HLT17GenerationBundle,
    HLT17InMemoryCheckpoint,
)
from .hlt17_member_codec import HLT17MemberEncoding
from .pro20_ev1_protocol import (
    EVENT_ORIGIN_HEX,
    EVENT_TARGET,
    EVENT_TARGET_HEX,
    PRODUCTION_NAMESPACE,
    PRO20EV1AuthorityReceipt,
    PRO20EV1ValidatedGeneration,
    build_authority_receipt,
    next_required_member_key,
)
from .pro20_ev1_store import (
    PRODUCTION_CONTAINER,
    PRO20EV1CampaignStore,
    PRO20EV1CampaignStoreError,
    PRO20EV1PostpublicationUncertainty,
    PRO20EV1TerminalStore,
    PRO20EV1UnclosedSession,
    PRO20EV1WriterCapability,
)
from .pro20_origin import MEMBER_KEYS, ORIGIN_TIME_HEX, capture_pro20_historical_origin
from .pro20_source_factory import (
    CAMPAIGN_ID,
    PROTOCOL,
    build_pro20_runtime_origin,
)
from .proto19_gr0_static_factory import STATIC_INPUT_PATHS
from .protocol_v19 import (
    BLOB_UNCHANGED_KINDS,
    CURSOR_CHANGED_RETRY_KINDS,
    RECORD_KINDS,
    TEMPORAL_EXHAUSTION_REASONS,
    TERMINAL_KINDS,
)
from .tdg11_imp1_ledger import imp1_checkpoint_extension


ARTIFACT_ID = "FGC-1-PRO20-EV1-FRZ1"
SCHEMA = "FGC-1-PRO20-EV1-runner-v1"
IMPLEMENTATION_RELATIVE = "src/recursive_horizons/fgc/evolution/pro20_ev1_runtime.py"
AUTHORITY_RELATIVE = "src/recursive_horizons/fgc/evolution/pro20_ev1_authority.py"
RUNNER_RELATIVE = "scripts/run_fgc_pro20_ev1.py"
AUTHORITY_VALIDATOR_MODULE = "recursive_horizons.fgc.evolution.pro20_ev1_authority"
AUTHORITY_VALIDATOR_NAME = "validate_authority_delta"
DEFAULT_HOST = "local"
OWNER_TOKEN = "pro20-ev1-one-shot"
CFL_MAXIMUM_HEX = "0x1.0000000000000p-3"
MAX_TOTAL_WALL_SECONDS = 86400.0
MAX_RSS_BYTES = 4294967296
MAX_NAMESPACE_BYTES = 17179869184
MIN_FREE_DISK_BYTES = 34359738368
MAX_PUBLISHED_ATTEMPTS = 4096
PLANNING_ACCEPTED_STEP_ESTIMATE = 224
SOURCE_RETRY_CAP = HLT17_SOURCE_RETRY_CAP
CFL_RETRY_CAP = HLT17_CFL_RETRY_CAP
FRESH_CAP_FORMULA = "min(remaining, cfl_maximum*grid_spacing/previous_speed_upper)"
AUTHORITY_SEAM = (
    "authority/delta validation is the later exact-delta freeze seam "
    f"{AUTHORITY_VALIDATOR_MODULE}.{AUTHORITY_VALIDATOR_NAME}; "
    "this parent does not accept a boolean skip/force flag"
)
INDEPENDENT_BINDER_SEAM = (
    "an independent PREF1 binder must reconstruct evidence without trusting "
    "runner labels; this runner is not that binder"
)
STATE_MACHINE = MappingProxyType(
    {
        "default_or_status": "read-only plan record; no source, store, or namespace",
        "authority": "fail closed before source construction or namespace creation",
        "origin": "pro20_origin then pro20_source_factory; one construction",
        "seed": "atomic production container publication",
        "writer": "one exclusive session; no resume or takeover",
        "schedule": "canonical unfinished member, or the single pending owner",
        "fresh_cap": FRESH_CAP_FORMULA,
        "pending_plan": "cursor.next_attempt_plan only",
        "publish": "map structural outcomes onto production protocol kinds",
        "temporal_exhaustion": (
            "unchanged accepted bundle plus non-executable terminal evidence"
        ),
        "typed_stop_after_seed": "immutable terminal evidence; close the writer",
        "unexpected": "abandon writer; leave the durable session unclosed",
        "success": "store terminal first_event_complete; calibration_eligible false",
    }
)
PROTOCOL_KIND_BY_STRUCTURAL_CLASS = MappingProxyType(
    {
        "accepted_physical_advance": "accepted_fine",
        "source_overlay_retry": "source_retry",
        "source_overlay_exhausted": "source_retry_exhausted",
        "cfl_overlay_retry": "cfl_retry",
        "cfl_overlay_exhausted": "cfl_retry_exhausted",
        "temporal_pending_successor": "temporal_retry",
        "temporal_exhausted_unchanged": "temporal_retry_exhausted",
        "typed_resource_stop": "resource_exhausted",
        "typed_premise_stop": "invalid_premise",
    }
)


class PRO20EV1RuntimeError(ValueError):
    """Typed PRO20-EV1 runner stop. Not a scientific classification."""

    outcome = "invalid"
    published = False
    phase = "invalid"
    abandon_writer = False
    physical_source_constructed = False


class PRO20EV1AuthoritySeamError(PRO20EV1RuntimeError):
    outcome = "invalid"
    phase = "authority_seam"


class PRO20EV1PremiseError(PRO20EV1RuntimeError):
    outcome = "invalid"
    phase = "premise"


class PRO20EV1ResourceError(PRO20EV1RuntimeError):
    outcome = "resource"
    phase = "resource"


class PRO20EV1UnexpectedError(PRO20EV1RuntimeError):
    outcome = "inconclusive"
    phase = "unexpected"
    abandon_writer = True


class PRO20EV1PostpublicationError(PRO20EV1RuntimeError):
    outcome = "inconclusive"
    phase = "postpublication"
    published = True
    abandon_writer = True


def _fail(message: str, cls: type[PRO20EV1RuntimeError] = PRO20EV1RuntimeError) -> None:
    raise cls(message)


def _sha256_hex(raw: bytes) -> str:
    if type(raw) is not bytes:
        raise PRO20EV1PremiseError("payload is not immutable bytes")
    return sha256(raw).hexdigest()


def _same_bits(left: float, right: float) -> bool:
    return float(left).hex() == float(right).hex()


def current_rss_bytes() -> int:
    usage = resource.getrusage(resource.RUSAGE_SELF)
    rss = int(usage.ru_maxrss)
    if rss < 0:
        _fail("RSS observation is negative", PRO20EV1ResourceError)
    if sys.platform == "darwin":
        return rss
    return rss * 1024


def free_disk_bytes(path: Path) -> int:
    try:
        stats = os.statvfs(os.fspath(path))
    except OSError as error:
        raise PRO20EV1ResourceError("free disk cannot be observed") from error
    available = int(stats.f_bavail) * int(stats.f_frsize)
    if available < 0:
        _fail("free disk observation is negative", PRO20EV1ResourceError)
    return available


def namespace_byte_count(root: Path) -> int:
    target = root.joinpath(*PRODUCTION_CONTAINER.split("/"))
    try:
        info = target.lstat()
    except FileNotFoundError:
        return 0
    except OSError as error:
        raise PRO20EV1ResourceError("namespace cannot be observed") from error
    if not (info.st_mode & 0o170000 == 0o040000):
        _fail("namespace is not a directory", PRO20EV1ResourceError)
    total = 0
    pending = [target]
    while pending:
        current = pending.pop()
        try:
            with os.scandir(current) as iterator:
                children = list(iterator)
        except OSError as error:
            raise PRO20EV1ResourceError("namespace cannot be scanned") from error
        for child in children:
            try:
                child_info = child.stat(follow_symlinks=False)
            except OSError as error:
                raise PRO20EV1ResourceError("namespace entry cannot be read") from error
            mode = child_info.st_mode
            if mode & 0o170000 == 0o120000:
                _fail("namespace contains a symlink", PRO20EV1ResourceError)
            if mode & 0o170000 == 0o040000:
                pending.append(Path(child.path))
                continue
            if mode & 0o170000 != 0o100000:
                _fail("namespace contains a nonregular node", PRO20EV1ResourceError)
            total += int(child_info.st_size)
    return total


def recommended_resource_ceilings() -> dict[str, object]:
    return {
        "max_rss_bytes": MAX_RSS_BYTES,
        "max_total_wall_seconds": MAX_TOTAL_WALL_SECONDS,
        "max_wall_seconds_per_member": {key: MAX_TOTAL_WALL_SECONDS for key in MEMBER_KEYS},
        "max_namespace_bytes": MAX_NAMESPACE_BYTES,
        "min_free_disk_bytes": MIN_FREE_DISK_BYTES,
        "max_published_attempts": MAX_PUBLISHED_ATTEMPTS,
        "source_retry_cap": SOURCE_RETRY_CAP,
        "cfl_retry_cap": CFL_RETRY_CAP,
        "scientific_threshold_fitted": False,
        "planning_accepted_step_estimate": PLANNING_ACCEPTED_STEP_ESTIMATE,
        "planning_accepted_step_estimate_is_not_a_gate": True,
        "coordinator_freeze_required": True,
    }


def plan_status() -> dict[str, object]:
    """Default CLI/status record. Performs no I/O, source call, or write."""

    identity = hlt17_implementation_identity()
    return {
        "artifact_id": ARTIFACT_ID,
        "schema": SCHEMA,
        "role": "implementation_parent_one_shot_runner",
        "protocol": PROTOCOL,
        "campaign_id": CAMPAIGN_ID,
        "namespace": PRODUCTION_NAMESPACE,
        "member_keys": list(MEMBER_KEYS),
        "event_origin_time_hex": ORIGIN_TIME_HEX,
        "event_target_hex": EVENT_TARGET_HEX,
        "fresh_cap_formula": FRESH_CAP_FORMULA,
        "pending_plans_from_cursor_only": True,
        "state_machine": dict(STATE_MACHINE),
        "resource_ceilings": recommended_resource_ceilings(),
        "implementation_identity": dict(identity),
        "authority_delta_validation": AUTHORITY_SEAM,
        "independent_binder_seam": INDEPENDENT_BINDER_SEAM,
        "physical_source_constructed": False,
        "campaign_execution_authorized": False,
        "calibration_eligible": False,
        "resume_authorized": False,
        "takeover_authorized": False,
        "boolean_flag_accepted": False,
    }


def implementation_identity(repository_root: Path | None = None) -> dict[str, object]:
    identity = hlt17_implementation_identity()
    record: dict[str, object] = {
        "artifact_id": ARTIFACT_ID,
        "schema": SCHEMA,
        "path": IMPLEMENTATION_RELATIVE,
        "protocol": PROTOCOL,
        "campaign_id": CAMPAIGN_ID,
        "c1r1_actual": dict(identity),
        "campaign_execution_authorized": False,
        "calibration_eligible": False,
    }
    if repository_root is not None:
        raw = read_regular_file(Path(repository_root), IMPLEMENTATION_RELATIVE)
        record["sha256"] = _sha256_hex(raw)
    return record


def _cursor_from_bundle(bundle: HLT17GenerationBundle) -> HLT17IMP1Cursor:
    mapping = load_canonical_json(bundle.cursor_bytes)
    if type(mapping) is not dict:
        _fail("bundle cursor is not a JSON object", PRO20EV1PremiseError)
    return restore_hlt17_cursor(mapping)


def _bundles_identical(left: HLT17GenerationBundle, right: HLT17GenerationBundle) -> bool:
    return (
        left.descriptor == right.descriptor
        and left.payload == right.payload
        and left.cursor_bytes == right.cursor_bytes
    )


def _encoding_changed(left: HLT17GenerationBundle, right: HLT17GenerationBundle) -> bool:
    return left.descriptor != right.descriptor or left.payload != right.payload


def fresh_requested_cap(checkpoint: HLT17InMemoryCheckpoint) -> float:
    """``min(remaining, cfl_maximum * grid_spacing / previous_speed_upper)``."""

    if not isinstance(checkpoint, HLT17InMemoryCheckpoint):
        raise TypeError("checkpoint must be HLT17InMemoryCheckpoint")
    checkpoint.agree()
    cursor = checkpoint.cursor
    if not cursor.public_fresh_permitted():
        _fail("fresh cap is not available on a pending cursor", PRO20EV1PremiseError)
    remaining = float(cursor.event_target) - float(cursor.accepted_time)
    if not math.isfinite(remaining) or remaining <= 0.0:
        _fail("fresh cap requires positive remaining time", PRO20EV1PremiseError)
    transaction = checkpoint.member.transaction
    spacing = float(transaction.grid_spacing)
    cfl_maximum = float(transaction.cfl_maximum)
    speed = float(transaction.causal_state.previous_speed_upper)
    if not all(
        math.isfinite(value) and value > 0.0
        for value in (spacing, cfl_maximum, speed)
    ):
        _fail("fresh cap inputs are not positive", PRO20EV1PremiseError)
    cfl_cap = cfl_maximum * spacing / speed
    if not math.isfinite(cfl_cap) or cfl_cap <= 0.0:
        _fail("CFL cap is not a positive finite width", PRO20EV1PremiseError)
    return remaining if remaining < cfl_cap else cfl_cap


def pending_plan_from_cursor(cursor: HLT17IMP1Cursor):
    """Return the exact cursor-owned next plan. Fresh cursors have none."""

    if not isinstance(cursor, HLT17IMP1Cursor):
        raise TypeError("cursor must be HLT17IMP1Cursor")
    cursor = revalidate_hlt17_cursor(cursor)
    if cursor.public_fresh_permitted():
        return None
    plan = cursor.next_attempt_plan()
    if plan is None:
        _fail(
            "pending retry cursor has no executable next plan",
            PRO20EV1PremiseError,
        )
    return plan


def schedule_next_member(generation: PRO20EV1ValidatedGeneration) -> str:
    """Canonical unfinished member, or the single pending retry owner."""

    return next_required_member_key(generation)


@dataclass(frozen=True, slots=True)
class AttemptClassification:
    """Protocol kind derived from bytes and typed stops, not runner labels."""

    kind: str
    member_key: str
    successor_bundle: HLT17GenerationBundle
    terminal_evidence: Mapping[str, object] | None
    accepted_state_advanced: bool
    executable_retry_cursor: bool
    physical_state_preserved: bool
    structural_class: str

    def __post_init__(self) -> None:
        if self.kind not in RECORD_KINDS or self.kind == "seed":
            _fail("classified kind is not a production attempt kind", PRO20EV1PremiseError)
        if self.member_key not in MEMBER_KEYS:
            _fail("classified member is outside the live cohort", PRO20EV1PremiseError)
        if type(self.accepted_state_advanced) is not bool:
            raise TypeError("accepted_state_advanced must be a built-in bool")
        if type(self.executable_retry_cursor) is not bool:
            raise TypeError("executable_retry_cursor must be a built-in bool")
        if type(self.physical_state_preserved) is not bool:
            raise TypeError("physical_state_preserved must be a built-in bool")
        if self.kind in BLOB_UNCHANGED_KINDS:
            if self.accepted_state_advanced or self.executable_retry_cursor:
                _fail(
                    f"{self.kind} cannot advance state or advertise an executable cursor",
                    PRO20EV1PremiseError,
                )
            if type(self.terminal_evidence) is not dict:
                _fail(f"{self.kind} omitted terminal evidence", PRO20EV1PremiseError)
            if self.terminal_evidence.get("kind") != self.kind:
                _fail(f"{self.kind} evidence kind differs", PRO20EV1PremiseError)
            if self.terminal_evidence.get("executable_retry_cursor") is not False:
                _fail(
                    f"{self.kind} evidence advertised an executable retry cursor",
                    PRO20EV1PremiseError,
                )
        if self.kind == "accepted_fine" and self.accepted_state_advanced is not True:
            _fail("accepted_fine requires accepted-state advance", PRO20EV1PremiseError)
        if self.kind in TERMINAL_KINDS and self.executable_retry_cursor:
            _fail("typed terminal cannot remain executable", PRO20EV1PremiseError)
        object.__setattr__(
            self,
            "terminal_evidence",
            None if self.terminal_evidence is None else dict(self.terminal_evidence),
        )


def _overlay_exhausted(cursor: HLT17IMP1Cursor) -> bool:
    overlay = cursor.overlay
    if overlay is None:
        return False
    if overlay.exhausted is True or overlay.next_plan is None:
        return True
    cap = SOURCE_RETRY_CAP if overlay.owner == SOURCE_OWNER else CFL_RETRY_CAP
    return overlay.owner_current_count > cap


def classify_attempt(
    *,
    member_key: str,
    predecessor_bundle: HLT17GenerationBundle,
    successor_bundle: HLT17GenerationBundle,
    accepted_state_advanced: bool,
    resource_stop: Mapping[str, object] | None = None,
    premise_stop: Mapping[str, object] | None = None,
    temporal_exhaustion_evidence: Mapping[str, object] | None = None,
) -> AttemptClassification:
    """Map one attempt from bytes and typed stops. Labels are not consulted."""

    if member_key not in MEMBER_KEYS:
        _fail("member_key is outside the live cohort", PRO20EV1PremiseError)
    if type(accepted_state_advanced) is not bool:
        raise TypeError("accepted_state_advanced must be a built-in bool")
    identical = _bundles_identical(predecessor_bundle, successor_bundle)
    encoding_changed = _encoding_changed(predecessor_bundle, successor_bundle)
    predecessor_cursor = _cursor_from_bundle(predecessor_bundle)
    successor_cursor = _cursor_from_bundle(successor_bundle)
    physical_preserved = (
        predecessor_cursor.physical_state_sha256 == successor_cursor.physical_state_sha256
        and _same_bits(predecessor_cursor.accepted_time, successor_cursor.accepted_time)
    )
    if resource_stop is not None:
        if not identical:
            _fail(
                "resource stop must publish the unchanged accepted bundle",
                PRO20EV1PremiseError,
            )
        extras = {
            key: value
            for key, value in dict(resource_stop).items()
            if key not in {"kind", "reason", "executable_retry_cursor", "spends_retry"}
        }
        evidence = {
            **extras,
            "kind": "resource_exhausted",
            "reason": str(resource_stop.get("reason") or "resource_ceiling"),
            "executable_retry_cursor": False,
            "spends_retry": False,
        }
        return AttemptClassification(
            kind="resource_exhausted",
            member_key=member_key,
            successor_bundle=predecessor_bundle,
            terminal_evidence=evidence,
            accepted_state_advanced=False,
            executable_retry_cursor=False,
            physical_state_preserved=True,
            structural_class="typed_resource_stop",
        )
    if premise_stop is not None:
        if not identical:
            _fail(
                "premise stop must publish the unchanged accepted bundle",
                PRO20EV1PremiseError,
            )
        evidence = {
            "kind": "invalid_premise",
            "reason": str(premise_stop.get("reason") or "invalid_premise"),
            "executable_retry_cursor": False,
            "spends_retry": False,
        }
        return AttemptClassification(
            kind="invalid_premise",
            member_key=member_key,
            successor_bundle=predecessor_bundle,
            terminal_evidence=evidence,
            accepted_state_advanced=False,
            executable_retry_cursor=False,
            physical_state_preserved=True,
            structural_class="typed_premise_stop",
        )
    if accepted_state_advanced:
        if identical or not encoding_changed or physical_preserved:
            _fail(
                "accepted advance did not change the accepted physical encoding",
                PRO20EV1PremiseError,
            )
        if successor_cursor.mode != FRESH_READY:
            _fail("accepted fine did not return a fresh cursor", PRO20EV1PremiseError)
        return AttemptClassification(
            kind="accepted_fine",
            member_key=member_key,
            successor_bundle=successor_bundle,
            terminal_evidence=None,
            accepted_state_advanced=True,
            executable_retry_cursor=False,
            physical_state_preserved=False,
            structural_class="accepted_physical_advance",
        )
    if encoding_changed:
        _fail(
            "non-fine attempt changed accepted encoding without an advance flag",
            PRO20EV1UnexpectedError,
        )
    if identical:
        if temporal_exhaustion_evidence is None:
            _fail(
                "unchanged bundle is not a typed terminal without evidence",
                PRO20EV1UnexpectedError,
            )
        evidence = dict(temporal_exhaustion_evidence)
        evidence["kind"] = "temporal_retry_exhausted"
        evidence["executable_retry_cursor"] = False
        if evidence.get("reason") not in TEMPORAL_EXHAUSTION_REASONS:
            _fail("temporal exhaustion reason differs", PRO20EV1PremiseError)
        return AttemptClassification(
            kind="temporal_retry_exhausted",
            member_key=member_key,
            successor_bundle=predecessor_bundle,
            terminal_evidence=evidence,
            accepted_state_advanced=False,
            executable_retry_cursor=False,
            physical_state_preserved=True,
            structural_class="temporal_exhausted_unchanged",
        )
    overlay = successor_cursor.overlay
    if overlay is not None:
        if not physical_preserved:
            _fail("source/CFL retry changed accepted physical state", PRO20EV1UnexpectedError)
        exhausted = _overlay_exhausted(successor_cursor)
        if overlay.owner == SOURCE_OWNER:
            kind = "source_retry_exhausted" if exhausted else "source_retry"
            structural = (
                "source_overlay_exhausted" if exhausted else "source_overlay_retry"
            )
        elif overlay.owner == CFL_OWNER:
            kind = "cfl_retry_exhausted" if exhausted else "cfl_retry"
            structural = "cfl_overlay_exhausted" if exhausted else "cfl_overlay_retry"
        else:
            _fail("overlay owner is not source or CFL", PRO20EV1UnexpectedError)
            raise AssertionError("unreachable")
        evidence = None
        if kind in TERMINAL_KINDS:
            evidence = {
                "kind": kind,
                "exhausted": True,
                "executable_retry_cursor": False,
                "next_plan": None,
            }
        return AttemptClassification(
            kind=kind,
            member_key=member_key,
            successor_bundle=successor_bundle,
            terminal_evidence=evidence,
            accepted_state_advanced=False,
            executable_retry_cursor=kind in CURSOR_CHANGED_RETRY_KINDS
            and kind not in TERMINAL_KINDS,
            physical_state_preserved=True,
            structural_class=structural,
        )
    if (
        successor_cursor.mode == RETRY_PENDING
        and successor_cursor.temporal is not None
        and successor_cursor.next_attempt_plan() is not None
    ):
        if not physical_preserved:
            _fail("temporal retry changed accepted physical state", PRO20EV1UnexpectedError)
        return AttemptClassification(
            kind="temporal_retry",
            member_key=member_key,
            successor_bundle=successor_bundle,
            terminal_evidence=None,
            accepted_state_advanced=False,
            executable_retry_cursor=True,
            physical_state_preserved=True,
            structural_class="temporal_pending_successor",
        )
    _fail(
        "cursor changed without an executable retry plan; temporal exhaustion "
        "must not be absorbed into an executable cursor",
        PRO20EV1UnexpectedError,
    )
    raise AssertionError("unreachable")


def restore_checkpoint_from_bundle(
    checkpoint: HLT17InMemoryCheckpoint,
    bundle: HLT17GenerationBundle,
) -> None:
    encoding = HLT17MemberEncoding(bundle.descriptor, bundle.payload)
    cursor = _cursor_from_bundle(bundle)
    checkpoint.restore(
        encoding,
        cursor,
        last_accepted_executed_plan=encoding.executed_plan,
    )


def _resource_evidence(reason: str, **fields: object) -> dict[str, object]:
    payload = {"reason": reason, **fields}
    return payload


def check_resource_ceilings(
    *,
    repository_root: Path,
    started: float,
    now: Callable[[], float],
    rss: Callable[[], int],
    disk_free: Callable[[Path], int],
    namespace_bytes: Callable[[Path], int],
    published_attempts: int,
) -> None:
    elapsed = float(now()) - float(started)
    if elapsed > MAX_TOTAL_WALL_SECONDS:
        _fail("max total wall exceeded", PRO20EV1ResourceError)
    if int(rss()) > MAX_RSS_BYTES:
        _fail("max RSS exceeded", PRO20EV1ResourceError)
    if int(disk_free(repository_root)) < MIN_FREE_DISK_BYTES:
        _fail("free disk is below the operational ceiling", PRO20EV1ResourceError)
    if int(namespace_bytes(repository_root)) > MAX_NAMESPACE_BYTES:
        _fail("namespace byte ceiling exceeded", PRO20EV1ResourceError)
    if published_attempts > MAX_PUBLISHED_ATTEMPTS:
        _fail("published-attempt ceiling exceeded", PRO20EV1ResourceError)


def execute_scheduled_attempt(
    checkpoint: HLT17InMemoryCheckpoint,
    *,
    requested_cap: float | None = None,
) -> AttemptClassification:
    """One prepare/admit/adopt path with temporal-exhaustion intercept."""

    if not isinstance(checkpoint, HLT17InMemoryCheckpoint):
        raise TypeError("checkpoint must be HLT17InMemoryCheckpoint")
    checkpoint.agree()
    predecessor = checkpoint.generation_bundle()
    member_key = checkpoint.member.key
    cursor = checkpoint.cursor
    if checkpoint.at_event_target():
        _fail("parked member cannot be advanced", PRO20EV1PremiseError)
    try:
        if cursor.public_fresh_permitted():
            cap = fresh_requested_cap(checkpoint) if requested_cap is None else float(
                requested_cap
            )
            if requested_cap is not None and cap.hex() != float(requested_cap).hex():
                cap = float(requested_cap)
            prepared = checkpoint.prepare(requested_cap=cap)
        else:
            if pending_plan_from_cursor(cursor) is None:
                _fail("pending cursor omitted its next plan", PRO20EV1PremiseError)
            prepared = checkpoint.prepare(requested_cap=requested_cap)
    except HLT17OwnerRetryExhausted as error:
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            premise_stop={
                "reason": (
                    f"persisted_{error.owner}_exhaustion_without_terminal:"
                    f"{error.reason}"
                )
            },
        )
    except HLT17InvalidPremise as error:
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            premise_stop={"reason": str(error) or "invalid_premise"},
        )
    except C1R1RuntimeResourceStop as error:
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            resource_stop=_resource_evidence(str(error) or "c1r1_runtime_resource"),
        )
    except (C1R1CoordinateLatticeStop, C1R1RefinementPathStop) as error:
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            premise_stop={"reason": str(error) or type(error).__name__},
        )
    if isinstance(prepared, HLT17CheckpointResult):
        _fail("prepare returned a checkpoint result instead of a bridge", PRO20EV1PremiseError)
    if not isinstance(prepared, HLT17BridgeResult):
        raise TypeError("prepare must return HLT17BridgeResult")
    if prepared.accepted_state_advanced:
        _fail("prepare advanced accepted state", PRO20EV1UnexpectedError)
    overlay_dispositions = {
        "source_retry_required",
        "cfl_retry_required",
        "source_retry_exhausted",
        "cfl_retry_exhausted",
    }
    if prepared.disposition in overlay_dispositions:
        absorbed = checkpoint.absorb_nonfine(prepared)
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=absorbed.bundle,
            accepted_state_advanced=False,
        )
    if prepared.disposition == "temporal_retry_required":
        absorbed = checkpoint.absorb_nonfine(prepared)
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=absorbed.bundle,
            accepted_state_advanced=False,
        )
    if prepared.disposition == "temporal_retry_exhausted":
        evidence = dict(prepared.evidence)
        evidence["kind"] = "temporal_retry_exhausted"
        evidence["executable_retry_cursor"] = False
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            temporal_exhaustion_evidence=evidence,
        )
    if prepared.disposition == "invalid_premise":
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            premise_stop={"reason": str(prepared.evidence.get("reason") or "invalid_premise")},
        )
    if prepared.prepared is None:
        _fail("prepare omitted an admissible family", PRO20EV1PremiseError)
    try:
        admitted = checkpoint.require_admission(prepared.prepared)
    except C1R1RuntimeResourceStop as error:
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            resource_stop=_resource_evidence(str(error) or "c1r1_runtime_resource"),
        )
    except C1R1TemporalRetryExhausted as error:
        mapping = imp1_checkpoint_extension(error.updated_ledger)
        evidence = {
            "kind": "temporal_retry_exhausted",
            "reason": error.reason,
            "updated_imp1_ledger": mapping,
            "updated_ledger_sha256": ledger_sha256(error.updated_ledger),
            "executable_retry_cursor": False,
        }
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            temporal_exhaustion_evidence=evidence,
        )
    except (HLT17InvalidPremise, HLT17OwnerRetryExhausted) as error:
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            premise_stop={"reason": str(error) or "invalid_premise"},
        )
    if isinstance(admitted, HLT17CheckpointResult):
        _fail("admission returned a checkpoint result", PRO20EV1PremiseError)
    if admitted.disposition in overlay_dispositions:
        absorbed = checkpoint.absorb_nonfine(admitted)
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=absorbed.bundle,
            accepted_state_advanced=False,
        )
    if admitted.disposition == "temporal_retry_required":
        absorbed = checkpoint.absorb_nonfine(admitted)
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=absorbed.bundle,
            accepted_state_advanced=False,
        )
    if admitted.disposition == "temporal_retry_exhausted":
        evidence = dict(admitted.evidence)
        evidence["kind"] = "temporal_retry_exhausted"
        evidence["executable_retry_cursor"] = False
        if evidence.get("reason") not in TEMPORAL_EXHAUSTION_REASONS:
            evidence["reason"] = str(evidence.get("reason") or "maximum_temporal_retries")
            if evidence["reason"] not in TEMPORAL_EXHAUSTION_REASONS:
                evidence["reason"] = "maximum_temporal_retries"
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            temporal_exhaustion_evidence=evidence,
        )
    if admitted.disposition == "invalid_premise":
        return classify_attempt(
            member_key=member_key,
            predecessor_bundle=predecessor,
            successor_bundle=predecessor,
            accepted_state_advanced=False,
            premise_stop={"reason": str(admitted.evidence.get("reason") or "invalid_premise")},
        )
    adopted = checkpoint.adopt_fine(admitted.prepared)
    if adopted.disposition != ACCEPTED_FINE or adopted.accepted_state_advanced is not True:
        _fail("fine adoption did not advance accepted state", PRO20EV1UnexpectedError)
    return classify_attempt(
        member_key=member_key,
        predecessor_bundle=predecessor,
        successor_bundle=adopted.bundle,
        accepted_state_advanced=True,
    )


@dataclass
class RuntimeHooks:
    """Injectable origin, step, store, and resource seams for focused tests."""

    origin_loader: Callable[[Path], object] | None = None
    runtime_origin_builder: Callable[..., object] | None = None
    static_input_loader: Callable[[Path], Mapping[str, bytes]] | None = None
    step_fn: Callable[[HLT17InMemoryCheckpoint], AttemptClassification] | None = None
    store_cls: type[PRO20EV1CampaignStore] = PRO20EV1CampaignStore
    now: Callable[[], float] = time.monotonic
    rss: Callable[[], int] = current_rss_bytes
    disk_free: Callable[[Path], int] = free_disk_bytes
    namespace_bytes: Callable[[Path], int] = namespace_byte_count
    authority_validator: Callable[..., Mapping[str, object]] | None = None
    implementation_identity_fn: Callable[[Path | None], Mapping[str, object]] = (
        implementation_identity
    )
    source_constructions: list[str] = field(default_factory=list)
    adoptions: list[str] = field(default_factory=list)


def _load_static_inputs(root: Path) -> dict[str, bytes]:
    payloads: dict[str, bytes] = {}
    for relative in STATIC_INPUT_PATHS:
        payloads[relative] = read_regular_file(root, relative)
    return payloads


def _require_namespace_absent(root: Path) -> None:
    container = root.joinpath(*PRODUCTION_CONTAINER.split("/"))
    try:
        container.lstat()
    except FileNotFoundError:
        parent = container.parent
        if not parent.exists():
            return
        try:
            names = os.listdir(parent)
        except OSError as error:
            raise PRO20EV1PremiseError("namespace parent cannot be listed") from error
        if any(
            name.startswith(".pro20-event1.evidence-io-stage-")
            or name.startswith(".pro20-event1.")
            for name in names
        ):
            _fail("staging namespace already exists", PRO20EV1PremiseError)
        return
    except OSError as error:
        raise PRO20EV1PremiseError("namespace cannot be read") from error
    _fail("production namespace already exists; one-shot authority is consumed")


def _receipt_from_authority(mapping: Mapping[str, object]) -> PRO20EV1AuthorityReceipt:
    hashes = mapping.get("receipt_identities")
    if not isinstance(hashes, Mapping):
        hashes = {
            "authority_sha256": mapping.get("authority_sha256"),
            "implementation_sha256": mapping.get("implementation_sha256"),
            "config_sha256": mapping.get("config_sha256"),
            "source_sha256": mapping.get("source_closure_sha256")
            or mapping.get("source_sha256"),
            "origin_sha256": mapping.get("origin_capture_sha256")
            or mapping.get("origin_sha256"),
            "environment_sha256": mapping.get("environment_sha256"),
        }
    return build_authority_receipt(
        campaign_id=str(mapping.get("campaign_id") or CAMPAIGN_ID),
        authority_sha256=str(hashes["authority_sha256"]),
        implementation_sha256=str(hashes["implementation_sha256"]),
        config_sha256=str(hashes["config_sha256"]),
        source_sha256=str(hashes["source_sha256"]),
        origin_sha256=str(hashes["origin_sha256"]),
        environment_sha256=str(hashes["environment_sha256"]),
    )


def _working_checkpoints(runtime_origin: object) -> dict[str, HLT17InMemoryCheckpoint]:
    members = getattr(runtime_origin, "members", None)
    if not isinstance(members, Mapping) or tuple(members) != MEMBER_KEYS:
        _fail("runtime origin member order differs", PRO20EV1PremiseError)
    working: dict[str, HLT17InMemoryCheckpoint] = {}
    for key in MEMBER_KEYS:
        record = members[key]
        checkpoint = getattr(record, "checkpoint", None)
        if not isinstance(checkpoint, HLT17InMemoryCheckpoint):
            _fail(f"{key} omitted a live checkpoint", PRO20EV1PremiseError)
        working[key] = checkpoint
    return working


def _require_origin_runtime_identities(
    origin: object,
    runtime_origin: object,
    *,
    receipt: PRO20EV1AuthorityReceipt,
    authority: Mapping[str, object],
) -> dict[str, object]:
    expected_origin = str(receipt.mapping["origin_sha256"])
    if getattr(origin, "capture_sha256", None) != expected_origin:
        _fail("live origin capture identity differs", PRO20EV1PremiseError)
    historical = authority.get("historical_source_tree")
    if not isinstance(historical, Mapping):
        _fail("authority omitted historical source-tree proof", PRO20EV1PremiseError)
    expected_bridge = historical.get("generation1_bridge_content_id")
    if getattr(origin, "generation1_bridge_content_id", None) != expected_bridge:
        _fail("live generation-one bridge identity differs", PRO20EV1PremiseError)
    if getattr(runtime_origin, "captured_origin_sha256", None) != expected_origin:
        _fail("runtime origin capture identity differs", PRO20EV1PremiseError)
    expected_source = str(receipt.mapping["source_sha256"])
    if getattr(runtime_origin, "source_closure_sha256", None) != expected_source:
        _fail("runtime source-closure identity differs", PRO20EV1PremiseError)
    members = getattr(runtime_origin, "members", None)
    if not isinstance(members, Mapping) or tuple(members) != MEMBER_KEYS:
        _fail("runtime origin member order differs", PRO20EV1PremiseError)
    expected_configurations = authority.get("source_configuration")
    if not isinstance(expected_configurations, Mapping):
        _fail("authority omitted source-configuration identities", PRO20EV1PremiseError)
    observed_configurations: dict[str, str] = {}
    for key in MEMBER_KEYS:
        record = members[key]
        configuration = getattr(record, "source_configuration_sha256", None)
        if configuration != expected_configurations.get(key):
            _fail(f"{key} source-configuration identity differs", PRO20EV1PremiseError)
        if getattr(record, "captured_origin_sha256", None) != expected_origin:
            _fail(f"{key} captured-origin identity differs", PRO20EV1PremiseError)
        if getattr(record, "source_closure_sha256", None) != expected_source:
            _fail(f"{key} source-closure identity differs", PRO20EV1PremiseError)
        observed_configurations[key] = str(configuration)
    construction = getattr(runtime_origin, "construction_sha256", None)
    if (
        type(construction) is not str
        or len(construction) != 64
        or construction.lower() != construction
        or any(character not in "0123456789abcdef" for character in construction)
    ):
        _fail("runtime construction identity is absent", PRO20EV1PremiseError)
    return {
        "origin_capture_sha256": expected_origin,
        "generation1_bridge_content_id": expected_bridge,
        "source_closure_sha256": expected_source,
        "source_configuration": observed_configurations,
        "runtime_construction_sha256": construction,
    }


def _seed_bundles(
    runtime_origin: object,
) -> dict[str, HLT17GenerationBundle]:
    members = getattr(runtime_origin, "members")
    return {key: members[key].bundle for key in MEMBER_KEYS}


def _load_authority_validator() -> Callable[..., Mapping[str, object]]:
    from importlib import import_module

    try:
        module = import_module(AUTHORITY_VALIDATOR_MODULE)
        validator = getattr(module, AUTHORITY_VALIDATOR_NAME)
    except (ImportError, AttributeError) as error:
        raise PRO20EV1AuthoritySeamError(AUTHORITY_SEAM) from error
    if not callable(validator):
        _fail(AUTHORITY_SEAM, PRO20EV1AuthoritySeamError)
    return validator


def _require_authority(
    *,
    repository_root: Path,
    authority_commit: str,
    implementation_sha256: str,
    validator: Callable[..., Mapping[str, object]] | None,
) -> Mapping[str, object]:
    loaded = validator if validator is not None else _load_authority_validator()
    try:
        receipt = loaded(
            repository_root=repository_root,
            authority_commit=authority_commit,
            implementation_sha256=implementation_sha256,
        )
    except PRO20EV1RuntimeError:
        raise
    except Exception as error:
        raise PRO20EV1AuthoritySeamError(
            "PRO20 exact-delta authority validation failed"
        ) from error
    if type(receipt) is bool:
        _fail(
            "authority/delta validation must not be a boolean flag",
            PRO20EV1AuthoritySeamError,
        )
    if not isinstance(receipt, Mapping):
        _fail("authority receipt is not a mapping", PRO20EV1AuthoritySeamError)
    return dict(receipt)


def _publish_typed_terminal(
    store: PRO20EV1CampaignStore,
    capability: PRO20EV1WriterCapability,
    classification: AttemptClassification,
) -> PRO20EV1ValidatedGeneration:
    return store.publish_attempt(
        capability,
        kind=classification.kind,
        member_key=classification.member_key,
        successor_bundle=classification.successor_bundle,
        terminal_evidence=classification.terminal_evidence,
    )


def run_first_event(
    repository_root: Path,
    *,
    authority_commit: str | None = None,
    hooks: RuntimeHooks | None = None,
    host: str = DEFAULT_HOST,
) -> dict[str, object]:
    """One-shot first-event run. Status/default never call this."""

    hooks = RuntimeHooks() if hooks is None else hooks
    root = Path(repository_root)
    if not root.is_absolute():
        _fail("repository root is not canonical", PRO20EV1PremiseError)
    started = hooks.now()
    published = False
    physical = False
    capability: PRO20EV1WriterCapability | None = None
    store: PRO20EV1CampaignStore | None = None
    record: dict[str, object] = {
        "artifact_id": ARTIFACT_ID,
        "schema": SCHEMA,
        "campaign_id": CAMPAIGN_ID,
        "namespace": PRODUCTION_NAMESPACE,
        "event_origin_time_hex": EVENT_ORIGIN_HEX,
        "event_target_hex": EVENT_TARGET.hex(),
        "physical_source_constructed": False,
        "campaign_execution_authorized": False,
        "calibration_eligible": False,
        "resume_authorized": False,
        "takeover_authorized": False,
        "planning_accepted_step_estimate": PLANNING_ACCEPTED_STEP_ESTIMATE,
        "planning_accepted_step_estimate_is_not_a_gate": True,
        "attempts": [],
    }
    try:
        _require_namespace_absent(root)
        check_resource_ceilings(
            repository_root=root,
            started=started,
            now=hooks.now,
            rss=hooks.rss,
            disk_free=hooks.disk_free,
            namespace_bytes=hooks.namespace_bytes,
            published_attempts=0,
        )
        implementation = hooks.implementation_identity_fn(root)
        if authority_commit is None:
            _fail("run refuses absent authority", PRO20EV1AuthoritySeamError)
        authority = _require_authority(
            repository_root=root,
            authority_commit=authority_commit,
            implementation_sha256=str(implementation["sha256"]),
            validator=hooks.authority_validator,
        )
        record["authority"] = {
            "authority_commit": authority.get("authority_commit"),
            "implementation_commit": authority.get("implementation_commit"),
            "physical_source_qualification_executed": False,
        }
        receipt = _receipt_from_authority(authority)
        origin_loader = hooks.origin_loader or capture_pro20_historical_origin
        static_loader = hooks.static_input_loader or _load_static_inputs
        builder = hooks.runtime_origin_builder or build_pro20_runtime_origin
        origin = origin_loader(root)
        hooks.source_constructions.append("origin")
        static_inputs = static_loader(root)
        source_closure = (
            getattr(origin, "source_closure_sha256", None)
            if hasattr(origin, "source_closure_sha256")
            else None
        )
        if source_closure is None and hasattr(receipt, "mapping"):
            source_closure = receipt.mapping.get("source_sha256")
        runtime_origin = builder(
            origin,
            static_input_bytes=static_inputs,
            source_closure_sha256=str(source_closure or receipt.mapping["source_sha256"]),
        )
        hooks.source_constructions.append("factory")
        physical = True
        record["physical_source_constructed"] = True
        runtime_identities = _require_origin_runtime_identities(
            origin,
            runtime_origin,
            receipt=receipt,
            authority=authority,
        )
        record["runtime_identities"] = runtime_identities
        if len([item for item in hooks.source_constructions if item == "factory"]) != 1:
            _fail("physical source factory was constructed more than once", PRO20EV1UnexpectedError)
        bundles = _seed_bundles(runtime_origin)
        working = _working_checkpoints(runtime_origin)
        store = hooks.store_cls.publish_seed_store(
            root, receipt=receipt, bundles=bundles
        )
        published = True
        capability = store.acquire_writer(owner_token=OWNER_TOKEN, host=host)
        attempts: list[dict[str, object]] = []
        while True:
            view = store.authenticate()
            last = view.generations[-1]
            if last.terminal or last.first_event_complete:
                break
            try:
                check_resource_ceilings(
                    repository_root=root,
                    started=started,
                    now=hooks.now,
                    rss=hooks.rss,
                    disk_free=hooks.disk_free,
                    namespace_bytes=hooks.namespace_bytes,
                    published_attempts=len(attempts),
                )
            except PRO20EV1ResourceError as error:
                member_key = schedule_next_member(last)
                classification = classify_attempt(
                    member_key=member_key,
                    predecessor_bundle=last.bundles[member_key],
                    successor_bundle=last.bundles[member_key],
                    accepted_state_advanced=False,
                    resource_stop=_resource_evidence(str(error)),
                )
                last = _publish_typed_terminal(store, capability, classification)
                attempts.append(
                    {
                        "member_key": member_key,
                        "kind": classification.kind,
                        "structural_class": classification.structural_class,
                        "store_generation": last.store_generation,
                        "terminal": True,
                    }
                )
                break
            if len(attempts) >= MAX_PUBLISHED_ATTEMPTS:
                member_key = schedule_next_member(last)
                classification = classify_attempt(
                    member_key=member_key,
                    predecessor_bundle=last.bundles[member_key],
                    successor_bundle=last.bundles[member_key],
                    accepted_state_advanced=False,
                    resource_stop=_resource_evidence("published-attempt ceiling"),
                )
                last = _publish_typed_terminal(store, capability, classification)
                attempts.append(
                    {
                        "member_key": member_key,
                        "kind": classification.kind,
                        "structural_class": classification.structural_class,
                    }
                )
                break
            member_key = schedule_next_member(last)
            step = hooks.step_fn or execute_scheduled_attempt
            try:
                restore_checkpoint_from_bundle(
                    working[member_key], last.bundles[member_key]
                )
                classification = step(working[member_key])
            except PRO20EV1ResourceError as error:
                classification = classify_attempt(
                    member_key=member_key,
                    predecessor_bundle=last.bundles[member_key],
                    successor_bundle=last.bundles[member_key],
                    accepted_state_advanced=False,
                    resource_stop=_resource_evidence(str(error)),
                )
            except PRO20EV1PremiseError as error:
                classification = classify_attempt(
                    member_key=member_key,
                    predecessor_bundle=last.bundles[member_key],
                    successor_bundle=last.bundles[member_key],
                    accepted_state_advanced=False,
                    premise_stop={"reason": str(error)},
                )
            if classification.kind == "accepted_fine":
                hooks.adoptions.append(member_key)
            last = store.publish_attempt(
                capability,
                kind=classification.kind,
                member_key=classification.member_key,
                successor_bundle=classification.successor_bundle,
                terminal_evidence=classification.terminal_evidence,
            )
            attempts.append(
                {
                    "member_key": classification.member_key,
                    "kind": classification.kind,
                    "structural_class": classification.structural_class,
                    "accepted_state_advanced": classification.accepted_state_advanced,
                    "executable_retry_cursor": classification.executable_retry_cursor,
                    "physical_state_preserved": classification.physical_state_preserved,
                    "store_generation": last.store_generation,
                    "disposition": last.disposition,
                    "terminal": last.terminal,
                    "first_event_complete": last.first_event_complete,
                    "calibration_eligible": False,
                }
            )
            if last.terminal or last.first_event_complete:
                break
        store.close_writer(capability)
        capability = None
        view = store.authenticate()
        last = view.generations[-1]
        record.update(
            {
                "attempts": attempts,
                "published_attempts": len(attempts),
                "terminal": last.terminal,
                "first_event_complete": last.first_event_complete,
                "disposition": last.disposition,
                "checkpoint_sha256": last.checkpoint_sha256,
                "journal_sha256": last.journal_sha256,
                "unclosed_session": view.unclosed_session,
                "calibration_eligible": False,
                "calibration_claimed": False,
                "physics_claimed": False,
                "wall_seconds": hooks.now() - started,
            }
        )
        return record
    except PRO20EV1PostpublicationUncertainty as error:
        if capability is not None and store is not None:
            store.abandon_writer(capability)
        wrapped = PRO20EV1PostpublicationError(str(error))
        wrapped.physical_source_constructed = physical
        raise wrapped from error
    except PRO20EV1UnexpectedError:
        if capability is not None and store is not None:
            store.abandon_writer(capability)
        raise
    except (PRO20EV1TerminalStore, PRO20EV1UnclosedSession, PRO20EV1CampaignStoreError) as error:
        if capability is not None and store is not None:
            store.abandon_writer(capability)
        wrapped = PRO20EV1UnexpectedError(str(error))
        wrapped.physical_source_constructed = physical
        raise wrapped from error
    except PRO20EV1RuntimeError as error:
        if capability is not None and store is not None:
            if published:
                try:
                    store.abandon_writer(capability)
                except Exception:
                    pass
        error.physical_source_constructed = physical
        raise
    except Exception as error:
        if capability is not None and store is not None:
            store.abandon_writer(capability)
        wrapped = PRO20EV1UnexpectedError(str(error))
        wrapped.physical_source_constructed = physical
        raise wrapped from error


def require_authority_delta(
    *,
    repository_root: Path,
    authority_commit: str,
    implementation_sha256: str,
    validator: Callable[..., Mapping[str, object]] | None = None,
) -> Mapping[str, object]:
    """Exact-delta freeze receipt. Never constructs source or a store."""

    return _require_authority(
        repository_root=Path(repository_root),
        authority_commit=authority_commit,
        implementation_sha256=implementation_sha256,
        validator=validator,
    )


__all__ = [
    "ARTIFACT_ID",
    "AttemptClassification",
    "AUTHORITY_SEAM",
    "CFL_MAXIMUM_HEX",
    "CAMPAIGN_ID",
    "DEFAULT_HOST",
    "FRESH_CAP_FORMULA",
    "INDEPENDENT_BINDER_SEAM",
    "MAX_NAMESPACE_BYTES",
    "MAX_PUBLISHED_ATTEMPTS",
    "MAX_RSS_BYTES",
    "MAX_TOTAL_WALL_SECONDS",
    "MIN_FREE_DISK_BYTES",
    "PLANNING_ACCEPTED_STEP_ESTIMATE",
    "PRODUCTION_NAMESPACE",
    "PRO20EV1AuthoritySeamError",
    "PRO20EV1PostpublicationError",
    "PRO20EV1PremiseError",
    "PRO20EV1ResourceError",
    "PRO20EV1RuntimeError",
    "PRO20EV1UnexpectedError",
    "PROTOCOL_KIND_BY_STRUCTURAL_CLASS",
    "RuntimeHooks",
    "SCHEMA",
    "STATE_MACHINE",
    "classify_attempt",
    "current_rss_bytes",
    "execute_scheduled_attempt",
    "fresh_requested_cap",
    "free_disk_bytes",
    "implementation_identity",
    "namespace_byte_count",
    "pending_plan_from_cursor",
    "plan_status",
    "recommended_resource_ceilings",
    "require_authority_delta",
    "restore_checkpoint_from_bundle",
    "run_first_event",
    "schedule_next_member",
]
