"""Prospective authority for a bounded TDG8 retry chain.

RCV1 repaired one historical metadata defect.  RCV2 binds the resulting
campaign once and authorizes the *existing* HLT16/TDG6 transition relation for
the remainder of the original GR-0 event.  It does not enumerate future
checkpoint hashes and never treats a retry, exhaustion, or event completion as
physics.

Only :func:`build_prelaunch` observes the exact mutable entry while producing
the compact result.  :func:`validate_compact` is deliberately store-blind.
The runner separately invokes :func:`inspect_descendant` immediately before a
status report or mutation.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tomllib
from types import MappingProxyType
from typing import Any, Iterable, Mapping

from . import hlt16_campaign_recovery as campaign_recovery
from . import hlt16_lifecycle as lifecycle
from . import proto15_runtime as p15
from . import tdg8_retry_recovery_authority as rcv1
from .hlt16_campaign_schema import (
    CampaignCheckpoint,
    validate_checkpoint_successor,
    validate_journal,
)
from .hlt16_campaign_store import HLT16CampaignStore, HLT16CampaignStoreError
from .numerical_engine import array_content_sha256
from .proto17_pure_construction import MEMBER_KEYS


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG8-RCV2-AUTH1"
PROJECT_VERSION = "0.11.0"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CLASSIFICATION = "premise_only_bounded_retry_chain_and_same_event_resume_authority"
CONFIG_PATH = "configs/fgc/fgc-1-tdg8-rcv2-auth1.toml"
RESULT_PATH = "results/fgc-1-tdg8-rcv2-auth1.json"
OWNER_DOCUMENT = "docs/fgc-tdg8-rcv2-auth1.md"

RCV1_AUTHORITY_COMMIT = "6036403e580474821a417c60e8dd8ddf73c52a17"
RCV1_RESULT_PATH = "results/fgc-1-tdg8-rcv1-auth1.json"
RCV1_RESULT_SHA256 = "5248f03f64c0a37d8cccc844d198ac6aa3ae7cde3c074c516142bcd7cc83b046"
ORIGINAL_EXECUTION_COMMIT = "6df67967e6e0bc63eca4dc6951a60d3787368f26"
ORIGINAL_PLAN_SHA256 = "e649e477b98411a3051f9946f2a9248b875d2b598818b8983ce06caa969d2126"
DESTINATION_PATH = "runs/fgc-2-sf1/proto19/calibration"
CAMPAIGN_ID = "FGC-2-SF1-PROTO18-GR0-A3-TDG8-SUCCESSOR-1"

ENTRY_CHECKPOINT_GENERATION = 8
ENTRY_CHECKPOINT_SHA256 = "a6401ff970aea002f8991c745a4cf7fe91c6d0a556907b04280d87e2c04d8936"
ENTRY_CHECKPOINT_JOURNAL_SEQUENCE = 8
ENTRY_CHECKPOINT_JOURNAL_SHA256 = "63af61720215223eb3a660b055a07fa897a317b9580128d1b8f79a7022bab86f"
ENTRY_REJECTION_SEQUENCE = 9
ENTRY_REJECTION_SHA256 = "f1d005bb01cf4ea473551ac9fdface8b582d9e5f20cdf63169193d9177ddabed"
ENTRY_TRANSITION_SEQUENCE = 10
ENTRY_TRANSITION_SHA256 = "5b533eb7009a9c7c3f353d26813f9cbbde8574afae77492929b435b4c701c5ff"
ENTRY_MEMBER_KEY = "RK4-2049"
ENTRY_DESCRIPTOR_SHA256 = "77847126340e78c8ac300fac2795bceda724c15e0894e34444b0da42a84166b4"
ENTRY_PHYSICAL_STATE_SHA256 = "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a"
ENTRY_ACCEPTED_TIME_HEX = "0x1.78554de5a30e0p+0"
ENTRY_RETRY_COUNT = 2
ENTRY_PENDING_CAP_HEX = "0x1.aaa9612df8000p-11"

FIRST_RECOVERY_GENERATION = 9
FIRST_RECOVERY_CHECKPOINT_SHA256 = "eb6fddc480c94aa7ed15c80399fc2b26637693efef02399f267075fc6c258e56"

EVENT = 23
TARGET = {"rational": "3/2", "binary64_hex": "0x1.8000000000000p+0"}
SUCCESSOR_EVENT = 24
SUCCESSOR_TARGET = {"rational": "25/16", "binary64_hex": "0x1.9000000000000p+0"}
TEMPORAL_RETRY_CAP = 32
MINIMUM_MACRO_STEP_HEX = "0x1.0000000000000p-30"

_MAX_CONFIG_BYTES = 1024 * 1024
_MAX_RESULT_BYTES = 4 * 1024 * 1024
_MAX_BOUND_BYTES = 32 * 1024 * 1024
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")

IMPLEMENTATION_INVENTORY = (
    ("authority", "src/recursive_horizons/fgc/evolution/tdg8_bounded_retry_authority.py"),
    ("runner", "scripts/run_fgc_tdg8_bounded_retry.py"),
    ("reproducer", "scripts/reproduce_fgc_tdg8_rcv2_auth1.py"),
    ("predecessor_authority", "src/recursive_horizons/fgc/evolution/tdg8_retry_recovery_authority.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_campaign_recovery.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_campaign_runtime.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_campaign_schema.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_campaign_store.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_lifecycle.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_progression_attempt.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/proto19_progression_contract.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/proto19_gr0_static_factory.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/tdg8_successor_runtime.py"),
    ("orchestrator", "scripts/run_fgc_pro19_event1.py"),
    ("owner_document", OWNER_DOCUMENT),
)

# A single RCV2 commit may contain the generic pending-retry recovery repair
# because the authority must bind that corrected implementation.  Nothing
# outside this exact task surface may hitchhike in the executable image.
AUTHORITY_DELTA_PATHS = tuple(sorted({
    "Makefile",
    CONFIG_PATH,
    RESULT_PATH,
    OWNER_DOCUMENT,
    "scripts/check_repo.py",
    "scripts/reproduce_fgc_tdg8_rcv2_auth1.py",
    "scripts/run_fgc_tdg8_bounded_retry.py",
    "src/recursive_horizons/fgc/evolution/hlt16_campaign_recovery.py",
    "src/recursive_horizons/fgc/evolution/tdg8_bounded_retry_authority.py",
    "tests/test_fgc_hlt16_campaign_recovery.py",
    "tests/test_fgc_tdg8_bounded_retry_authority.py",
    "tests/test_fgc_tdg8_bounded_retry_runner.py",
}))

NONCLAIMS = (
    "RCV2 authorizes one bounded GR-0 event continuation, not a common-event, calibration, trapping, candidate, mechanism, or physical result.",
    "A TDG6 rejection and its metadata recovery preserve the last accepted physical state; a later independently admitted fine step may advance that state under the original plan.",
    "Temporal exhaustion remains invalid_implementation_or_nonconverged_run and is never a physical obstruction or FGC-QR result.",
    "RCV2 never reimports generation zero, rechecks destination absence, opens SGB-L or FGC-QR, changes a threshold, or authorizes retained-EFT or transition claims.",
)
SCOPE = {
    "status_read_only": True,
    "one_prelaunch_live_observation": True,
    "compact_verifier_store_blind": True,
    "generic_metadata_reconciliation": True,
    "same_event_gr0_resume": True,
    "destination_absence_recheck": False,
    "generation_zero_reimport": False,
    "candidate_branches_forbidden": True,
    "threshold_or_retry_rule_change": False,
    "calibration_result_earned": False,
    "physical_result_earned": False,
}
CLAIMS = {
    "bounded_retry_chain_authorized": True,
    "original_plan_preserved": True,
    "physical_state_advanced_by_recovery": False,
    "common_event_completed": False,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "physical_result_earned": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
}

_ALLOWED_SUFFIXES = frozenset({
    (),
    ("accepted_fine",),
    ("source_rejection",),
    ("cfl_rejection",),
    ("tdg6_rejection",),
    ("tdg6_rejection", "cursor_transition"),
    ("terminal_lock",),
    ("source_rejection", "terminal_lock"),
    ("cfl_rejection", "terminal_lock"),
    ("tdg6_rejection", "terminal_lock"),
    ("common_event_commit",),
})


class TDG8BoundedRetryAuthorityError(RuntimeError):
    """The compact authority, entry, or authenticated descendant differs."""


@dataclass(frozen=True, slots=True)
class TDG8BoundedRetryExecutionAuthority:
    authority_commit: str
    predecessor_authority_commit: str
    original_execution_commit: str
    original_plan_sha256: str
    config_sha256: str
    result_sha256: str
    campaign_id: str
    destination_path: str
    entry_checkpoint_sha256: str
    entry_rejection_sha256: str
    entry_transition_sha256: str
    first_recovery_checkpoint_sha256: str
    temporal_retry_cap: int
    minimum_macro_step_hex: str
    implementation_inventory: tuple[tuple[str, str, str], ...]
    environment: Mapping[str, str]
    candidate_branches_forbidden: bool = True


@dataclass(frozen=True, slots=True)
class BoundedRetryBoundary:
    checkpoint_generation: int
    checkpoint_sha256: str
    journal_sequence: int
    journal_tip_sha256: str
    event: int
    target: Mapping[str, str]
    disposition: str
    suffix_kinds: tuple[str, ...]
    suffix_classification: str
    active_write: bool
    writer_state: str | None
    entry_ancestor_present: bool
    first_recovery_ancestor_present: bool
    temporal_pending_members: tuple[str, ...]
    terminal: bool
    terminal_lock_present: bool
    safe_to_restart: bool
    recoverable_under_rcv2: bool


def canonical(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _decode_json(raw: bytes, label: str) -> dict[str, Any]:
    def no_duplicates(items: Iterable[tuple[str, object]]) -> dict[str, object]:
        answer: dict[str, object] = {}
        for key, value in items:
            if key in answer:
                raise ValueError(key)
            answer[key] = value
        return answer

    try:
        value = json.loads(
            raw.decode("utf-8"), object_pairs_hook=no_duplicates,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TDG8BoundedRetryAuthorityError(f"{label} is malformed") from exc
    if not isinstance(value, dict):
        raise TDG8BoundedRetryAuthorityError(f"{label} is not object-valued")
    return value


def _expected_predecessor() -> dict[str, object]:
    return {
        "rcv1_authority_commit": RCV1_AUTHORITY_COMMIT,
        "rcv1_result_path": RCV1_RESULT_PATH,
        "rcv1_result_sha256": RCV1_RESULT_SHA256,
        "original_execution_commit": ORIGINAL_EXECUTION_COMMIT,
        "original_plan_sha256": ORIGINAL_PLAN_SHA256,
    }


def _expected_entry() -> dict[str, object]:
    return {
        "destination_path": DESTINATION_PATH,
        "campaign_id": CAMPAIGN_ID,
        "branch": "GR-0",
        "amplitude": "3",
        "event": EVENT,
        "target_rational": TARGET["rational"],
        "target_binary64_hex": TARGET["binary64_hex"],
        "checkpoint_generation": ENTRY_CHECKPOINT_GENERATION,
        "checkpoint_sha256": ENTRY_CHECKPOINT_SHA256,
        "checkpoint_journal_sequence": ENTRY_CHECKPOINT_JOURNAL_SEQUENCE,
        "checkpoint_journal_sha256": ENTRY_CHECKPOINT_JOURNAL_SHA256,
        "rejection_sequence": ENTRY_REJECTION_SEQUENCE,
        "rejection_sha256": ENTRY_REJECTION_SHA256,
        "transition_sequence": ENTRY_TRANSITION_SEQUENCE,
        "transition_sha256": ENTRY_TRANSITION_SHA256,
        "member_key": ENTRY_MEMBER_KEY,
        "descriptor_sha256": ENTRY_DESCRIPTOR_SHA256,
        "physical_state_sha256": ENTRY_PHYSICAL_STATE_SHA256,
        "accepted_time_binary64_hex": ENTRY_ACCEPTED_TIME_HEX,
        "retry_count": ENTRY_RETRY_COUNT,
        "pending_cap_binary64_hex": ENTRY_PENDING_CAP_HEX,
        "writer_state": "stale_verified_lock",
    }


def _expected_first_reconciliation() -> dict[str, object]:
    return {
        "checkpoint_generation": FIRST_RECOVERY_GENERATION,
        "checkpoint_sha256": FIRST_RECOVERY_CHECKPOINT_SHA256,
        "journal_sequence": ENTRY_TRANSITION_SEQUENCE,
        "journal_sha256": ENTRY_TRANSITION_SHA256,
        "member_mode": "RETRY_PENDING",
        "pending_owner": "temporal",
        "pending_cap_binary64_hex": ENTRY_PENDING_CAP_HEX,
        "retry_count": ENTRY_RETRY_COUNT,
        "physical_state_advanced": False,
    }


def _expected_bounded_policy() -> dict[str, object]:
    return {
        "temporal_retry_cap_per_macro_step": TEMPORAL_RETRY_CAP,
        "minimum_macro_step_binary64_hex": MINIMUM_MACRO_STEP_HEX,
        "retry_factor_numerator": 1,
        "retry_factor_denominator": 2,
        "same_event_only": True,
        "event_successor": SUCCESSOR_EVENT,
        "event_successor_target_rational": SUCCESSOR_TARGET["rational"],
        "event_successor_target_binary64_hex": SUCCESSOR_TARGET["binary64_hex"],
        "authenticated_descendants_allowed": True,
        "recognized_publication_cuts_allowed": True,
        "rejection_preserves_accepted_descriptor": True,
        "rejection_preserves_accepted_time": True,
    }


def _parse_config(raw: bytes) -> dict[str, Any]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise TDG8BoundedRetryAuthorityError("RCV2 config is malformed") from exc
    if not isinstance(value, dict) or set(value) != {
        "schema_version", "artifact_id", "project_version", "target_protocol",
        "classification", "nonclaims", "predecessor", "entry",
        "first_reconciliation", "bounded_policy", "scope", "claims",
        "implementation_inventory",
    }:
        raise TDG8BoundedRetryAuthorityError("RCV2 config fields differ")
    if (
        value["schema_version"] != SCHEMA_VERSION
        or value["artifact_id"] != ARTIFACT_ID
        or value["project_version"] != PROJECT_VERSION
        or value["target_protocol"] != TARGET_PROTOCOL
        or value["classification"] != CLASSIFICATION
        or value["nonclaims"] != list(NONCLAIMS)
        or value["predecessor"] != _expected_predecessor()
        or value["entry"] != _expected_entry()
        or value["first_reconciliation"] != _expected_first_reconciliation()
        or value["bounded_policy"] != _expected_bounded_policy()
        or value["scope"] != SCOPE
        or value["claims"] != CLAIMS
    ):
        raise TDG8BoundedRetryAuthorityError("RCV2 config contract differs")
    rows = value["implementation_inventory"]
    if (
        not isinstance(rows, list)
        or len(rows) != len(IMPLEMENTATION_INVENTORY)
        or tuple((item.get("role"), item.get("path")) for item in rows if isinstance(item, Mapping))
        != IMPLEMENTATION_INVENTORY
        or any(not isinstance(item, Mapping) or set(item) != {"role", "path"} for item in rows)
    ):
        raise TDG8BoundedRetryAuthorityError("RCV2 implementation inventory differs")
    return value


def _require_frozen_policy() -> None:
    if (
        p15.TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP != TEMPORAL_RETRY_CAP
        or float(p15.TDG6_MINIMUM_MACRO_STEP).hex() != MINIMUM_MACRO_STEP_HEX
    ):
        raise TDG8BoundedRetryAuthorityError("RCV2 frozen TDG6 policy differs")


def _derive_first_recovery(previous: CampaignCheckpoint,
                           records: tuple[Mapping[str, Any], ...]) -> CampaignCheckpoint:
    """Derive the checkpoint already fixed by the exact two-record entry."""
    if len(records) != 2:
        raise TDG8BoundedRetryAuthorityError("RCV2 entry record count differs")
    rejection, transition = records
    if (rejection.get("kind"), transition.get("kind")) != (
        "tdg6_rejection", "cursor_transition",
    ):
        raise TDG8BoundedRetryAuthorityError("RCV2 entry suffix differs")
    key = str(transition["payload"]["member_key"])
    before = previous.members[key]
    try:
        ledger = p15._updated_retry_ledger(
            p15._ledger_from_mapping(before.ledger),
            rejection["payload"]["evidence"],
        )
        successor_cursor = dict(transition["payload"]["successor_cursor"])
        pending = successor_cursor["retry_successor_payload_or_none"]
        members = dict(previous.members)
        members[key] = replace(
            before,
            cursor=successor_cursor,
            ledger=p15._ledger_mapping(ledger),
            pending_owner="temporal",
            pending_cap_hex=str(pending["half_cap"]),
        )
        current = CampaignCheckpoint(
            plan_sha256=previous.plan_sha256,
            authorization_commit=previous.authorization_commit,
            protocol=previous.protocol,
            campaign_id=previous.campaign_id,
            event=previous.event,
            target=dict(previous.target),
            members=members,
            journal_tip_sha256=str(transition["record_sha256"]),
            parent_sha256=previous.sha256,
            generation=previous.generation + 1,
            journal_sequence=int(transition["sequence"]),
            disposition="nonterminal",
            terminal=None,
        )
        validate_checkpoint_successor(previous, current, records)
    except Exception as exc:
        raise TDG8BoundedRetryAuthorityError(
            "RCV2 exact entry does not determine generation nine"
        ) from exc
    return current


def _inventory(root: Path) -> list[dict[str, str]]:
    answer = []
    for role, relative in IMPLEMENTATION_INVENTORY:
        raw = rcv1._read_leaf(root, relative, f"RCV2 {role}", _MAX_BOUND_BYTES)
        answer.append({"role": role, "path": relative, "sha256": _sha(raw)})
    return answer


def inspect_live_entry(root: Path) -> dict[str, Any]:
    """Read the one prospective entry without publishing or taking a lease."""
    repository = rcv1._repository(root)
    store = HLT16CampaignStore(repository / DESTINATION_PATH)
    try:
        snapshot = store.authenticated_snapshot()
        status = store.inspect_recovery()
    except HLT16CampaignStoreError as exc:
        raise TDG8BoundedRetryAuthorityError("RCV2 live entry does not validate") from exc
    previous = snapshot.checkpoint
    records = tuple(dict(item) for item in snapshot.suffix_records)
    if (
        previous.generation != ENTRY_CHECKPOINT_GENERATION
        or previous.sha256 != ENTRY_CHECKPOINT_SHA256
        or previous.journal_sequence != ENTRY_CHECKPOINT_JOURNAL_SEQUENCE
        or previous.journal_tip_sha256 != ENTRY_CHECKPOINT_JOURNAL_SHA256
        or previous.authorization_commit != ORIGINAL_EXECUTION_COMMIT
        or previous.plan_sha256 != ORIGINAL_PLAN_SHA256
        or previous.campaign_id != CAMPAIGN_ID
        or previous.protocol != TARGET_PROTOCOL
        or previous.event != EVENT
        or previous.target != TARGET
        or previous.disposition != "nonterminal"
        or snapshot.suffix_kinds != ("tdg6_rejection", "cursor_transition")
        or snapshot.staging_paths
        or snapshot.orphan_descriptor_sha256 is not None
        or snapshot.orphan_payload_semantic_sha256 is not None
        or status.writer_state != "stale_verified_lock"
    ):
        raise TDG8BoundedRetryAuthorityError("RCV2 prospective entry differs")
    rejection, transition = records
    if (
        rejection.get("sequence") != ENTRY_REJECTION_SEQUENCE
        or rejection.get("record_sha256") != ENTRY_REJECTION_SHA256
        or transition.get("sequence") != ENTRY_TRANSITION_SEQUENCE
        or transition.get("record_sha256") != ENTRY_TRANSITION_SHA256
        or transition.get("previous_record_sha256") != ENTRY_REJECTION_SHA256
    ):
        raise TDG8BoundedRetryAuthorityError("RCV2 entry journal identity differs")
    try:
        persisted = store.load_state(ENTRY_DESCRIPTOR_SHA256)
        physical = array_content_sha256(
            persisted.arrays["u"], persisted.arrays["p"], persisted.arrays["q"],
        )
        preview = campaign_recovery._resolve_tdg6_rejection(
            previous, rejection, persisted.arrays,
        )
    except Exception as exc:
        raise TDG8BoundedRetryAuthorityError("RCV2 entry physical state cannot be restored") from exc
    recovered = preview.checkpoint
    member = recovered.members[ENTRY_MEMBER_KEY]
    if (
        tuple(preview.records) != records
        or preview.evolution_state_sha256 != physical
        or preview.physical_state_advanced
        or recovered.generation != FIRST_RECOVERY_GENERATION
        or recovered.sha256 != FIRST_RECOVERY_CHECKPOINT_SHA256
        or recovered.journal_sequence != ENTRY_TRANSITION_SEQUENCE
        or recovered.journal_tip_sha256 != ENTRY_TRANSITION_SHA256
        or member.descriptor_sha256 != ENTRY_DESCRIPTOR_SHA256
        or physical != ENTRY_PHYSICAL_STATE_SHA256
        or member.cursor["accepted_boundary_time"]["binary64_hex"] != ENTRY_ACCEPTED_TIME_HEX
        or member.cursor["mode"] != "RETRY_PENDING"
        or member.pending_owner != "temporal"
        or member.pending_cap_hex != ENTRY_PENDING_CAP_HEX
        or member.ledger["current_macro_step_temporal_retry_count"] != ENTRY_RETRY_COUNT
    ):
        raise TDG8BoundedRetryAuthorityError("RCV2 first reconciliation differs")
    return {
        "entry": _expected_entry(),
        "first_reconciliation": _expected_first_reconciliation(),
        "entry_suffix_kinds": list(snapshot.suffix_kinds),
        "entry_physical_state_advanced": False,
    }


def _payload(config: Mapping[str, Any], inventory: list[dict[str, str]],
             live: Mapping[str, Any], environment: Mapping[str, str]) -> dict[str, Any]:
    return {
        "predecessor": dict(config["predecessor"]),
        "entry": dict(live["entry"]),
        "first_reconciliation": dict(live["first_reconciliation"]),
        "entry_suffix_kinds": list(live["entry_suffix_kinds"]),
        "entry_physical_state_advanced": live["entry_physical_state_advanced"],
        "bounded_policy": dict(config["bounded_policy"]),
        "implementation_inventory": inventory,
        "environment": dict(environment),
        "scope": dict(SCOPE),
        "claims": dict(CLAIMS),
        "nonclaims": list(NONCLAIMS),
    }


def build_prelaunch(config_raw: bytes, root: Path) -> dict[str, Any]:
    repository = rcv1._repository(root)
    config = _parse_config(config_raw)
    _require_frozen_policy()
    if rcv1._read_leaf(repository, CONFIG_PATH, "RCV2 config", _MAX_CONFIG_BYTES) != config_raw:
        raise TDG8BoundedRetryAuthorityError("live RCV2 config bytes differ")
    try:
        (repository / RESULT_PATH).lstat()
    except FileNotFoundError:
        pass
    else:
        raise TDG8BoundedRetryAuthorityError("RCV2 result already exists")
    predecessor = rcv1._read_leaf(
        repository, RCV1_RESULT_PATH, "RCV1 result", _MAX_RESULT_BYTES,
    )
    if _sha(predecessor) != RCV1_RESULT_SHA256:
        raise TDG8BoundedRetryAuthorityError("RCV1 result changed")
    live = inspect_live_entry(repository)
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL,
        "classification": CLASSIFICATION,
        "gate_status": "pass",
        "source_config_sha256": _sha(config_raw),
        "artifact_payload": _payload(
            config, _inventory(repository), live, rcv1.observed_environment(),
        ),
    }


def _parse_result(raw: bytes) -> dict[str, Any]:
    result = _decode_json(raw, "RCV2 result")
    payload = result.get("artifact_payload")
    if not isinstance(payload, Mapping):
        raise TDG8BoundedRetryAuthorityError("RCV2 result payload differs")
    rows = payload.get("implementation_inventory")
    if not isinstance(rows, list) or len(rows) != len(IMPLEMENTATION_INVENTORY):
        raise TDG8BoundedRetryAuthorityError("RCV2 result inventory differs")
    for row, expected in zip(rows, IMPLEMENTATION_INVENTORY, strict=True):
        if (
            not isinstance(row, Mapping)
            or set(row) != {"role", "path", "sha256"}
            or (row.get("role"), row.get("path")) != expected
            or not isinstance(row.get("sha256"), str)
            or not _SHA256.fullmatch(str(row["sha256"]))
        ):
            raise TDG8BoundedRetryAuthorityError("RCV2 result inventory row differs")
    return result


def validate_compact(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    """Validate tracked bytes without reading Git, a run store, or arrays."""
    config = _parse_config(config_raw)
    result = _parse_result(result_raw)
    payload = result["artifact_payload"]
    fixed_live = {
        "entry": _expected_entry(),
        "first_reconciliation": _expected_first_reconciliation(),
        "entry_suffix_kinds": ["tdg6_rejection", "cursor_transition"],
        "entry_physical_state_advanced": False,
    }
    inventory = [dict(item) for item in payload["implementation_inventory"]]
    expected = {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL,
        "classification": CLASSIFICATION,
        "gate_status": "pass",
        "source_config_sha256": _sha(config_raw),
        "artifact_payload": _payload(
            config, inventory, fixed_live, payload.get("environment", {}),
        ),
    }
    if canonical(result) != result_raw or canonical(expected) != result_raw:
        raise TDG8BoundedRetryAuthorityError("RCV2 compact result differs")
    return result


def _require_git_image(root: Path, commit: str) -> None:
    if not _COMMIT.fullmatch(commit):
        raise TDG8BoundedRetryAuthorityError("RCV2 authority commit is malformed")
    head = rcv1._git(root, "rev-parse", "--verify", "HEAD").stdout.decode("ascii").strip()
    resolved = rcv1._git(root, "rev-parse", "--verify", f"{commit}^{{commit}}").stdout.decode("ascii").strip()
    if head != commit or resolved != commit:
        raise TDG8BoundedRetryAuthorityError("RCV2 authority commit differs from HEAD")
    ancestor = subprocess.run(
        (
            "git", "--no-replace-objects", "--no-optional-locks", "-C", str(root),
            "merge-base", "--is-ancestor", RCV1_AUTHORITY_COMMIT, commit,
        ),
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=False,
        env={
            **{key: value for key, value in os.environ.items() if not key.startswith("GIT_")},
            "GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C",
        },
    )
    if ancestor.returncode != 0:
        raise TDG8BoundedRetryAuthorityError("RCV1 authority is not an ancestor")
    changed = rcv1._git(
        root, "diff", "--name-only", "-z", "--no-renames",
        f"{RCV1_AUTHORITY_COMMIT}..{commit}", "--",
    ).stdout
    try:
        paths = tuple(sorted(item.decode("utf-8") for item in changed.split(b"\0") if item))
    except UnicodeDecodeError as exc:
        raise TDG8BoundedRetryAuthorityError("RCV2 Git delta is not UTF-8") from exc
    if paths != AUTHORITY_DELTA_PATHS:
        raise TDG8BoundedRetryAuthorityError("RCV2 authority Git delta differs")
    dirty = subprocess.run(
        (
            "git", "--no-replace-objects", "--no-optional-locks",
            "-c", "core.fsmonitor=false", "-c", "core.untrackedCache=false",
            "-C", str(root), "diff", "--quiet", "--no-ext-diff",
            "--ignore-submodules=all", "HEAD", "--",
        ),
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=False,
        env={
            **{key: value for key, value in os.environ.items() if not key.startswith("GIT_")},
            "GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C",
        },
    )
    if dirty.returncode != 0:
        raise TDG8BoundedRetryAuthorityError("tracked tree differs from RCV2 authority commit")


def authorize_execution(root: Path, config_raw: bytes, result_raw: bytes,
                        authority_commit: str) -> TDG8BoundedRetryExecutionAuthority:
    repository = rcv1._repository(root)
    result = validate_compact(config_raw, result_raw)
    _require_frozen_policy()
    _require_git_image(repository, authority_commit)
    for path, expected in ((CONFIG_PATH, config_raw), (RESULT_PATH, result_raw)):
        live = rcv1._read_leaf(repository, path, f"RCV2 committed {path}", _MAX_RESULT_BYTES)
        if live != expected or rcv1._committed_bytes(repository, authority_commit, path) != expected:
            raise TDG8BoundedRetryAuthorityError(f"RCV2 committed bytes differ: {path}")
    predecessor = rcv1._read_leaf(repository, RCV1_RESULT_PATH, "RCV1 result", _MAX_RESULT_BYTES)
    if (
        _sha(predecessor) != RCV1_RESULT_SHA256
        or rcv1._committed_bytes(repository, RCV1_AUTHORITY_COMMIT, RCV1_RESULT_PATH) != predecessor
        or rcv1._committed_bytes(repository, authority_commit, RCV1_RESULT_PATH) != predecessor
    ):
        raise TDG8BoundedRetryAuthorityError("RCV1 result binding differs")
    rows: list[tuple[str, str, str]] = []
    for item in result["artifact_payload"]["implementation_inventory"]:
        role, path, digest = str(item["role"]), str(item["path"]), str(item["sha256"])
        live = rcv1._read_leaf(repository, path, f"RCV2 implementation {path}", _MAX_BOUND_BYTES)
        if _sha(live) != digest or rcv1._committed_bytes(repository, authority_commit, path) != live:
            raise TDG8BoundedRetryAuthorityError(f"RCV2 implementation bytes differ: {path}")
        rows.append((role, path, digest))
    environment = result["artifact_payload"].get("environment")
    if environment != rcv1.observed_environment():
        raise TDG8BoundedRetryAuthorityError("RCV2 execution environment differs")
    return TDG8BoundedRetryExecutionAuthority(
        authority_commit=authority_commit,
        predecessor_authority_commit=RCV1_AUTHORITY_COMMIT,
        original_execution_commit=ORIGINAL_EXECUTION_COMMIT,
        original_plan_sha256=ORIGINAL_PLAN_SHA256,
        config_sha256=_sha(config_raw),
        result_sha256=_sha(result_raw),
        campaign_id=CAMPAIGN_ID,
        destination_path=DESTINATION_PATH,
        entry_checkpoint_sha256=ENTRY_CHECKPOINT_SHA256,
        entry_rejection_sha256=ENTRY_REJECTION_SHA256,
        entry_transition_sha256=ENTRY_TRANSITION_SHA256,
        first_recovery_checkpoint_sha256=FIRST_RECOVERY_CHECKPOINT_SHA256,
        temporal_retry_cap=TEMPORAL_RETRY_CAP,
        minimum_macro_step_hex=MINIMUM_MACRO_STEP_HEX,
        implementation_inventory=tuple(rows),
        environment=MappingProxyType(dict(environment)),
    )


def _read_historical_json(store: HLT16CampaignStore, relative: str,
                          label: str) -> dict[str, Any]:
    try:
        return _decode_json(
            rcv1._read_leaf(store.root, relative, label, _MAX_RESULT_BYTES), label,
        )
    except (
        OSError,
        TDG8BoundedRetryAuthorityError,
        rcv1.TDG8RetryRecoveryAuthorityError,
    ) as exc:
        raise TDG8BoundedRetryAuthorityError(f"{label} is unavailable") from exc


def _require_entry_ancestors(store: HLT16CampaignStore, latest: CampaignCheckpoint) -> bool:
    checkpoint_path = (
        f"checkpoints/{ENTRY_CHECKPOINT_GENERATION:020d}-{ENTRY_CHECKPOINT_SHA256}.json"
    )
    checkpoint = _read_historical_json(store, checkpoint_path, "RCV2 entry checkpoint")
    if (
        checkpoint.get("checkpoint_sha256") != ENTRY_CHECKPOINT_SHA256
        or checkpoint.get("generation") != ENTRY_CHECKPOINT_GENERATION
        or checkpoint.get("journal_sequence") != ENTRY_CHECKPOINT_JOURNAL_SEQUENCE
        or checkpoint.get("journal_tip_sha256") != ENTRY_CHECKPOINT_JOURNAL_SHA256
        or checkpoint.get("authorization_commit") != ORIGINAL_EXECUTION_COMMIT
        or checkpoint.get("plan_sha256") != ORIGINAL_PLAN_SHA256
        or checkpoint.get("campaign_id") != CAMPAIGN_ID
    ):
        raise TDG8BoundedRetryAuthorityError("RCV2 entry checkpoint ancestor differs")
    for sequence, digest, kind in (
        (ENTRY_REJECTION_SEQUENCE, ENTRY_REJECTION_SHA256, "tdg6_rejection"),
        (ENTRY_TRANSITION_SEQUENCE, ENTRY_TRANSITION_SHA256, "cursor_transition"),
    ):
        relative = f"journal/{sequence:020d}-{digest}.journal"
        record = validate_journal(
            _read_historical_json(store, relative, f"RCV2 entry {kind}")
        )
        if record["record_sha256"] != digest or record["sequence"] != sequence or record["kind"] != kind:
            raise TDG8BoundedRetryAuthorityError("RCV2 entry journal ancestor differs")
    if latest.generation >= FIRST_RECOVERY_GENERATION:
        path = (
            f"checkpoints/{FIRST_RECOVERY_GENERATION:020d}-"
            f"{FIRST_RECOVERY_CHECKPOINT_SHA256}.json"
        )
        recovered = _read_historical_json(store, path, "RCV2 generation-nine checkpoint")
        if (
            recovered.get("checkpoint_sha256") != FIRST_RECOVERY_CHECKPOINT_SHA256
            or recovered.get("generation") != FIRST_RECOVERY_GENERATION
            or recovered.get("journal_sequence") != ENTRY_TRANSITION_SEQUENCE
            or recovered.get("journal_tip_sha256") != ENTRY_TRANSITION_SHA256
            or recovered.get("parent_sha256") != ENTRY_CHECKPOINT_SHA256
        ):
            raise TDG8BoundedRetryAuthorityError("RCV2 generation-nine ancestor differs")
        return True
    return False


def _require_same_event(checkpoint: CampaignCheckpoint) -> None:
    if checkpoint.event == EVENT:
        if checkpoint.target != TARGET:
            raise TDG8BoundedRetryAuthorityError("RCV2 event-23 target differs")
        return
    if checkpoint.event == SUCCESSOR_EVENT:
        if checkpoint.disposition != "event_complete" or checkpoint.target != SUCCESSOR_TARGET:
            raise TDG8BoundedRetryAuthorityError("RCV2 event-24 boundary differs")
        return
    raise TDG8BoundedRetryAuthorityError("RCV2 lineage escaped the one-event scope")


def _require_temporal_pending(store: HLT16CampaignStore,
                              checkpoint: CampaignCheckpoint) -> tuple[str, ...]:
    pending_keys: list[str] = []
    for key in MEMBER_KEYS:
        state = checkpoint.members[key]
        cursor = p15.Proto15Cursor(dict(state.cursor))
        cursor.validate()
        if cursor.mode != "RETRY_PENDING":
            continue
        pending_keys.append(key)
        if state.pending_owner != "temporal" or state.pending_cap_hex is None:
            raise TDG8BoundedRetryAuthorityError("RCV2 pending temporal owner differs")
        cap = float.fromhex(state.pending_cap_hex)
        count = int(state.ledger["current_macro_step_temporal_retry_count"])
        if (
            not math.isfinite(cap)
            or cap < float(p15.TDG6_MINIMUM_MACRO_STEP)
            or count < 1
            or count > TEMPORAL_RETRY_CAP
        ):
            raise TDG8BoundedRetryAuthorityError("RCV2 pending temporal bound differs")
        try:
            persisted = store.load_state(state.descriptor_sha256)
            physical = array_content_sha256(
                persisted.arrays["u"], persisted.arrays["p"], persisted.arrays["q"],
            )
            lifecycle.validate_retry_payload(
                cursor, p15._ledger_from_mapping(state.ledger),
                evolution_state_sha256=physical,
            )
        except Exception as exc:
            raise TDG8BoundedRetryAuthorityError("RCV2 pending temporal payload differs") from exc
    return tuple(pending_keys)


def inspect_descendant(root: Path,
                       receipt: TDG8BoundedRetryExecutionAuthority) -> BoundedRetryBoundary:
    """Authenticate any lawful descendant without mutating the store."""
    repository = rcv1._repository(root)
    if (
        not isinstance(receipt, TDG8BoundedRetryExecutionAuthority)
        or receipt.original_execution_commit != ORIGINAL_EXECUTION_COMMIT
        or receipt.original_plan_sha256 != ORIGINAL_PLAN_SHA256
        or receipt.campaign_id != CAMPAIGN_ID
        or receipt.temporal_retry_cap != TEMPORAL_RETRY_CAP
        or receipt.minimum_macro_step_hex != MINIMUM_MACRO_STEP_HEX
        or not receipt.candidate_branches_forbidden
    ):
        raise TDG8BoundedRetryAuthorityError("RCV2 typed execution authority differs")
    _require_frozen_policy()
    store = HLT16CampaignStore(repository / DESTINATION_PATH)
    try:
        snapshot = store.authenticated_snapshot()
        status = store.inspect_recovery()
    except HLT16CampaignStoreError as exc:
        raise TDG8BoundedRetryAuthorityError("RCV2 descendant store does not validate") from exc
    checkpoint = snapshot.checkpoint
    if (
        checkpoint.generation < ENTRY_CHECKPOINT_GENERATION
        or checkpoint.authorization_commit != ORIGINAL_EXECUTION_COMMIT
        or checkpoint.plan_sha256 != ORIGINAL_PLAN_SHA256
        or checkpoint.campaign_id != CAMPAIGN_ID
        or checkpoint.protocol != TARGET_PROTOCOL
        or tuple(checkpoint.members) != MEMBER_KEYS
    ):
        raise TDG8BoundedRetryAuthorityError("RCV2 descendant authority differs")
    _require_same_event(checkpoint)
    if tuple(snapshot.suffix_kinds) not in _ALLOWED_SUFFIXES:
        raise TDG8BoundedRetryAuthorityError("RCV2 descendant suffix differs")
    if checkpoint.generation == ENTRY_CHECKPOINT_GENERATION:
        if (
            checkpoint.sha256 != ENTRY_CHECKPOINT_SHA256
            or tuple(snapshot.suffix_kinds) != ("tdg6_rejection", "cursor_transition")
            or tuple(record["record_sha256"] for record in snapshot.suffix_records)
            != (ENTRY_REJECTION_SHA256, ENTRY_TRANSITION_SHA256)
        ):
            raise TDG8BoundedRetryAuthorityError("RCV2 entry descendant differs")
    first_recovery = _require_entry_ancestors(store, checkpoint)
    try:
        entry = store.load_state(ENTRY_DESCRIPTOR_SHA256)
        if array_content_sha256(
            entry.arrays["u"], entry.arrays["p"], entry.arrays["q"],
        ) != ENTRY_PHYSICAL_STATE_SHA256:
            raise TDG8BoundedRetryAuthorityError("RCV2 entry physical ancestor differs")
    except TDG8BoundedRetryAuthorityError:
        raise
    except Exception as exc:
        raise TDG8BoundedRetryAuthorityError("RCV2 entry physical ancestor cannot be restored") from exc
    pending = _require_temporal_pending(store, checkpoint)
    recoverable_writer = (
        status.writer_state == "stale_verified_lock"
        or (
            not status.active_write
            and snapshot.suffix_classification == "clean"
        )
    )
    return BoundedRetryBoundary(
        checkpoint_generation=checkpoint.generation,
        checkpoint_sha256=checkpoint.sha256,
        journal_sequence=checkpoint.journal_sequence,
        journal_tip_sha256=checkpoint.journal_tip_sha256,
        event=checkpoint.event,
        target=dict(checkpoint.target),
        disposition=checkpoint.disposition,
        suffix_kinds=tuple(snapshot.suffix_kinds),
        suffix_classification=snapshot.suffix_classification,
        active_write=status.active_write,
        writer_state=status.writer_state,
        entry_ancestor_present=True,
        first_recovery_ancestor_present=first_recovery,
        temporal_pending_members=pending,
        terminal=status.terminal or checkpoint.disposition in {
            "scientific_terminal", "invalid_terminal",
        },
        terminal_lock_present=snapshot.terminal_lock_present,
        safe_to_restart=status.safe_to_restart,
        recoverable_under_rcv2=(
            recoverable_writer
            and not snapshot.terminal_lock_present
            and checkpoint.event == EVENT
        ),
    )


__all__ = [
    "ARTIFACT_ID", "AUTHORITY_DELTA_PATHS", "BoundedRetryBoundary",
    "CAMPAIGN_ID", "CONFIG_PATH", "DESTINATION_PATH",
    "FIRST_RECOVERY_CHECKPOINT_SHA256", "IMPLEMENTATION_INVENTORY",
    "MINIMUM_MACRO_STEP_HEX", "ORIGINAL_EXECUTION_COMMIT",
    "ORIGINAL_PLAN_SHA256", "RCV1_AUTHORITY_COMMIT", "RESULT_PATH",
    "TDG8BoundedRetryAuthorityError", "TDG8BoundedRetryExecutionAuthority",
    "TEMPORAL_RETRY_CAP", "authorize_execution", "build_prelaunch",
    "canonical", "inspect_descendant", "inspect_live_entry",
    "validate_compact",
]
