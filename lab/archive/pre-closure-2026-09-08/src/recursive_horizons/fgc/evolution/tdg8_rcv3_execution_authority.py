"""Prospective authority for one bounded event on the RCV3 projection.

PREF1 authenticated the installed generation-nine projection without exposing
an evolution entry point.  AUTH1 consumes that immutable binder, preserves the
original internal HLT16 campaign identity, and authorizes only the existing
GR-0 event-23 transition relation through the event-24 boundary.

The compact verifier is deliberately store-blind.  A runner must separately
call :func:`inspect_descendant` immediately before reporting status or taking
the writer lease.  This module has no generation-zero, bootstrap, SGB-L,
FGC-QR, calibration-promotion, or physical-result entry point.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import re
import stat
import subprocess
import tomllib
from types import MappingProxyType
from typing import Any, Iterable, Mapping

from . import proto15_runtime as p15
from . import tdg8_bounded_retry_authority as rcv2
from . import tdg8_rcv3_pref1_binder as pref1
from .hlt16_campaign_store import HLT16CampaignStore, HLT16CampaignStoreError
from .numerical_engine import array_content_sha256
from .proto17_pure_construction import MEMBER_KEYS


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG8-RCV3-AUTH1"
PROJECT_VERSION = "0.11.0"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CLASSIFICATION = "premise_only_authenticated_generation9_same_event_resume_authority"
CONFIG_PATH = "configs/fgc/fgc-1-tdg8-rcv3-auth1.toml"
RESULT_PATH = "results/fgc-1-tdg8-rcv3-auth1.json"
OWNER_DOCUMENT = "docs/fgc-tdg8-rcv3-auth1.md"

PREF1_AUTHORITY_COMMIT = "46abf80753a53e13cf7206fb76d90dc28e4be242"
PREF1_CONFIG_PATH = pref1.CONFIG_PATH
PREF1_CONFIG_SHA256 = "2cdbf7e4ae97bd8cdba2e4a3c47cd798007fd2212c6267c169164232d93f3626"
PREF1_RESULT_PATH = pref1.RESULT_PATH
PREF1_RESULT_SHA256 = "c61bcb6ad65fec73cc43e6ecc92c7516234b77560c21b0dac1bb72b50160fa0a"

PROJECTION_ID = pref1.PROJECTION_ID
DESTINATION_WRAPPER = pref1.DESTINATION_WRAPPER
DESTINATION_PATH = pref1.DESTINATION_STORE
RECEIPT_PATH = pref1.RECEIPT_PATH
RECEIPT_SHA256 = pref1.RECEIPT_SHA256
RECEIPT_RAW_SHA256 = pref1.RECEIPT_RAW_SHA256

CAMPAIGN_ID = pref1.CAMPAIGN_ID
ORIGINAL_EXECUTION_COMMIT = pref1.INTERNAL_AUTHORIZATION_COMMIT
ORIGINAL_PLAN_SHA256 = pref1.PLAN_SHA256
ENTRY_CHECKPOINT_GENERATION = pref1.ANCHOR_GENERATION
ENTRY_CHECKPOINT_SHA256 = pref1.ANCHOR_CHECKPOINT_SHA256
ENTRY_CHECKPOINT_RAW_SHA256 = pref1.ANCHOR_CHECKPOINT_RAW_SHA256
ENTRY_JOURNAL_SEQUENCE = pref1.ANCHOR_JOURNAL_SEQUENCE
ENTRY_JOURNAL_SHA256 = pref1.ANCHOR_JOURNAL_SHA256
ENTRY_JOURNAL_RAW_SHA256 = pref1.ANCHOR_JOURNAL_RAW_SHA256
ENTRY_DESCRIPTOR_SHA256 = str(pref1.EXPECTED_MEMBERS[0]["descriptor_sha256"])
ENTRY_PHYSICAL_STATE_SHA256 = str(pref1.EXPECTED_MEMBERS[0]["physical_state_sha256"])
ENTRY_PENDING_CAP_HEX = str(pref1.EXPECTED_MEMBERS[0]["pending_cap_hex"])
ENTRY_RETRY_COUNT = int(pref1.EXPECTED_MEMBERS[0]["retry_count"])

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
    ("authority", "src/recursive_horizons/fgc/evolution/tdg8_rcv3_execution_authority.py"),
    ("runner", "scripts/run_fgc_tdg8_rcv3_event.py"),
    ("reproducer", "scripts/reproduce_fgc_tdg8_rcv3_auth1.py"),
    ("predecessor_authority", "src/recursive_horizons/fgc/evolution/tdg8_bounded_retry_authority.py"),
    ("predecessor_authority", "src/recursive_horizons/fgc/evolution/tdg8_retry_recovery_authority.py"),
    ("predecessor_binder", pref1.RUNTIME_PATH.replace("fork_runtime", "pref1_binder")),
    ("projection_runtime", pref1.RUNTIME_PATH),
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

# AUTH1 is committed as a single narrow checkpoint immediately after PREF1.
# Exact-delta verification prevents unrelated executable changes from joining
# the authority image unnoticed.
AUTHORITY_DELTA_PATHS = tuple(sorted({
    "Makefile",
    "README.md",
    CONFIG_PATH,
    "docs/claim-ledger.md",
    "docs/fgc-runtime-matrix.md",
    OWNER_DOCUMENT,
    "docs/research-roadmap.md",
    "paper/fgc-local-defocusing/README.md",
    "results/README.md",
    RESULT_PATH,
    "scripts/check_repo.py",
    "scripts/reproduce_fgc_tdg8_rcv3_auth1.py",
    "scripts/run_fgc_tdg8_rcv3_event.py",
    "src/recursive_horizons/fgc/evolution/tdg8_rcv3_execution_authority.py",
    "tests/test_check_repo_tdg8_rcv3_auth1.py",
    "tests/test_fgc_tdg8_rcv3_execution_authority.py",
    "tests/test_fgc_tdg8_rcv3_event_runner.py",
}))

NONCLAIMS = (
    "AUTH1 authorizes only the installed RCV3 GR-0 generation-nine continuation through event 23; it does not complete an event, trapping calibration, candidate trajectory, mechanism test, or physical result.",
    "The internal HLT16 campaign ID, original execution commit, plan bytes, generation-nine checkpoint, and sequence-ten journal remain unchanged; projection identity is external and is bound by the receipt, PREF1, and AUTH1.",
    "The existing depth-two temporal retry, halved cap, retry ceiling 32, minimum macro-step, TDG6/TDG7 admission, source/CFL ownership, and scientific stops are preserved without relaxation.",
    "A retry or metadata reconciliation never advances physical state; only an independently admitted accepted-fine transition may do so.",
    "Retry exhaustion, a provenance failure, or an invalid runtime is not a GR-0, FGC-QR, gradient, mechanism, or physical obstruction.",
    "AUTH1 never reimports generation zero, reruns bootstrap, resumes the terminal source, rechecks destination absence, opens SGB-L or FGC-QR, or authorizes retained-EFT or transition claims.",
)
SCOPE = {
    "status_read_only": True,
    "compact_verifier_store_blind": True,
    "live_projection_reauthenticated_before_mutation": True,
    "same_event_gr0_resume": True,
    "source_terminal_store_read_only": True,
    "source_campaign_resume_authorized": False,
    "bootstrap_reexecution": False,
    "destination_absence_recheck": False,
    "generation_zero_reimport": False,
    "candidate_branches_forbidden": True,
    "threshold_or_retry_rule_change": False,
    "physical_input_change": False,
    "calibration_result_earned": False,
    "physical_result_earned": False,
}
CLAIMS = {
    "PREF1_projection_authentication_bound": True,
    "bounded_retry_chain_authorized": True,
    "original_internal_plan_preserved": True,
    "physical_state_advanced_by_projection": False,
    "execution_authorized_by_artifact_alone": False,
    "execution_started": False,
    "common_event_completed": False,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "physical_result_earned": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
}


class TDG8RCV3ExecutionAuthorityError(RuntimeError):
    """The prospective contract, immutable binder, or live descendant differs."""


@dataclass(frozen=True, slots=True)
class TDG8RCV3ExecutionAuthority:
    authority_commit: str
    pref1_authority_commit: str
    pref1_result_sha256: str
    projection_id: str
    receipt_sha256: str
    original_execution_commit: str
    original_plan_sha256: str
    config_sha256: str
    result_sha256: str
    campaign_id: str
    destination_path: str
    entry_checkpoint_sha256: str
    entry_journal_sha256: str
    temporal_retry_cap: int
    minimum_macro_step_hex: str
    implementation_inventory: tuple[tuple[str, str, str], ...]
    environment: Mapping[str, str]
    candidate_branches_forbidden: bool = True


@dataclass(frozen=True, slots=True)
class RCV3BoundedBoundary:
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
    external_receipt_authenticated: bool
    entry_ancestor_present: bool
    temporal_pending_members: tuple[str, ...]
    terminal: bool
    terminal_lock_present: bool
    safe_to_restart: bool


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
        value = json.loads(raw.decode("utf-8"), object_pairs_hook=no_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TDG8RCV3ExecutionAuthorityError(f"{label} is malformed") from exc
    if not isinstance(value, dict):
        raise TDG8RCV3ExecutionAuthorityError(f"{label} is not object-valued")
    return value


def _expected_predecessor() -> dict[str, object]:
    return {
        "pref1_authority_commit": PREF1_AUTHORITY_COMMIT,
        "pref1_config_path": PREF1_CONFIG_PATH,
        "pref1_config_sha256": PREF1_CONFIG_SHA256,
        "pref1_result_path": PREF1_RESULT_PATH,
        "pref1_result_sha256": PREF1_RESULT_SHA256,
    }


def _expected_projection() -> dict[str, object]:
    return {
        "projection_id": PROJECTION_ID,
        "wrapper_path": DESTINATION_WRAPPER,
        "store_path": DESTINATION_PATH,
        "receipt_path": RECEIPT_PATH,
        "receipt_sha256": RECEIPT_SHA256,
        "receipt_raw_sha256": RECEIPT_RAW_SHA256,
    }


def _expected_internal_boundary() -> dict[str, object]:
    return {
        "campaign_id": CAMPAIGN_ID,
        "authorization_commit": ORIGINAL_EXECUTION_COMMIT,
        "plan_sha256": ORIGINAL_PLAN_SHA256,
        "branch": "GR-0",
        "amplitude": "3",
        "event": EVENT,
        "target_rational": TARGET["rational"],
        "target_binary64_hex": TARGET["binary64_hex"],
        "checkpoint_generation": ENTRY_CHECKPOINT_GENERATION,
        "checkpoint_sha256": ENTRY_CHECKPOINT_SHA256,
        "checkpoint_raw_sha256": ENTRY_CHECKPOINT_RAW_SHA256,
        "journal_sequence": ENTRY_JOURNAL_SEQUENCE,
        "journal_sha256": ENTRY_JOURNAL_SHA256,
        "journal_raw_sha256": ENTRY_JOURNAL_RAW_SHA256,
        "retry_member_key": "RK4-2049",
        "retry_descriptor_sha256": ENTRY_DESCRIPTOR_SHA256,
        "retry_physical_state_sha256": ENTRY_PHYSICAL_STATE_SHA256,
        "retry_depth": ENTRY_RETRY_COUNT,
        "pending_cap_binary64_hex": ENTRY_PENDING_CAP_HEX,
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
        "recognized_publication_cuts_allowed": True,
        "rejection_preserves_accepted_descriptor": True,
        "rejection_preserves_accepted_time": True,
    }


def _parse_config(raw: bytes) -> dict[str, Any]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 AUTH1 config is malformed") from exc
    expected_fields = {
        "schema_version", "artifact_id", "project_version", "target_protocol",
        "classification", "nonclaims", "predecessor", "projection",
        "internal_boundary", "bounded_policy", "scope", "claims",
        "implementation_inventory",
    }
    if not isinstance(value, dict) or set(value) != expected_fields:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 AUTH1 config fields differ")
    if (
        value["schema_version"] != SCHEMA_VERSION
        or value["artifact_id"] != ARTIFACT_ID
        or value["project_version"] != PROJECT_VERSION
        or value["target_protocol"] != TARGET_PROTOCOL
        or value["classification"] != CLASSIFICATION
        or value["nonclaims"] != list(NONCLAIMS)
        or value["predecessor"] != _expected_predecessor()
        or value["projection"] != _expected_projection()
        or value["internal_boundary"] != _expected_internal_boundary()
        or value["bounded_policy"] != _expected_bounded_policy()
        or value["scope"] != SCOPE
        or value["claims"] != CLAIMS
    ):
        raise TDG8RCV3ExecutionAuthorityError("RCV3 AUTH1 config contract differs")
    rows = value["implementation_inventory"]
    if (
        not isinstance(rows, list)
        or len(rows) != len(IMPLEMENTATION_INVENTORY)
        or any(not isinstance(item, Mapping) or set(item) != {"role", "path"} for item in rows)
        or tuple((item["role"], item["path"]) for item in rows) != IMPLEMENTATION_INVENTORY
    ):
        raise TDG8RCV3ExecutionAuthorityError("RCV3 AUTH1 implementation inventory differs")
    return value


def _require_frozen_policy() -> None:
    try:
        rcv2._require_frozen_policy()
    except Exception as exc:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 AUTH1 frozen policy differs") from exc


def _inventory(root: Path) -> list[dict[str, str]]:
    answer: list[dict[str, str]] = []
    for role, relative in IMPLEMENTATION_INVENTORY:
        raw = rcv2.rcv1._read_leaf(root, relative, f"RCV3 AUTH1 {role}", _MAX_BOUND_BYTES)
        answer.append({"role": role, "path": relative, "sha256": _sha(raw)})
    return answer


def _validate_pref1_bytes(repository: Path) -> tuple[bytes, bytes]:
    config_raw = rcv2.rcv1._read_leaf(
        repository, PREF1_CONFIG_PATH, "RCV3 PREF1 config", _MAX_CONFIG_BYTES,
    )
    result_raw = rcv2.rcv1._read_leaf(
        repository, PREF1_RESULT_PATH, "RCV3 PREF1 result", _MAX_RESULT_BYTES,
    )
    if _sha(config_raw) != PREF1_CONFIG_SHA256 or _sha(result_raw) != PREF1_RESULT_SHA256:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 PREF1 tracked identities differ")
    try:
        pref1.validate_compact_result(config_raw, result_raw)
    except Exception as exc:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 PREF1 compact result is invalid") from exc
    return config_raw, result_raw


def _payload(config: Mapping[str, Any], inventory: list[dict[str, str]],
             environment: Mapping[str, str]) -> dict[str, Any]:
    return {
        "predecessor": dict(config["predecessor"]),
        "projection": dict(config["projection"]),
        "internal_boundary": dict(config["internal_boundary"]),
        "bounded_policy": dict(config["bounded_policy"]),
        "implementation_inventory": inventory,
        "environment": dict(environment),
        "scope": dict(SCOPE),
        "claims": dict(CLAIMS),
        "nonclaims": list(NONCLAIMS),
    }


def build_prelaunch(config_raw: bytes, root: Path) -> dict[str, Any]:
    """Build compact authority evidence without opening either run store."""
    repository = rcv2.rcv1._repository(root)
    config = _parse_config(config_raw)
    _require_frozen_policy()
    if rcv2.rcv1._read_leaf(
        repository, CONFIG_PATH, "RCV3 AUTH1 config", _MAX_CONFIG_BYTES,
    ) != config_raw:
        raise TDG8RCV3ExecutionAuthorityError("live RCV3 AUTH1 config bytes differ")
    try:
        (repository / RESULT_PATH).lstat()
    except FileNotFoundError:
        pass
    else:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 AUTH1 result already exists")
    pref1_config, pref1_result = _validate_pref1_bytes(repository)
    if (
        rcv2.rcv1._committed_bytes(repository, PREF1_AUTHORITY_COMMIT, PREF1_CONFIG_PATH)
        != pref1_config
        or rcv2.rcv1._committed_bytes(repository, PREF1_AUTHORITY_COMMIT, PREF1_RESULT_PATH)
        != pref1_result
    ):
        raise TDG8RCV3ExecutionAuthorityError("RCV3 PREF1 immutable commit binding differs")
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL,
        "classification": CLASSIFICATION,
        "gate_status": "pass",
        "source_config_sha256": _sha(config_raw),
        "artifact_payload": _payload(
            config, _inventory(repository), rcv2.rcv1.observed_environment(),
        ),
    }


def _parse_result(raw: bytes) -> dict[str, Any]:
    result = _decode_json(raw, "RCV3 AUTH1 result")
    payload = result.get("artifact_payload")
    if not isinstance(payload, Mapping):
        raise TDG8RCV3ExecutionAuthorityError("RCV3 AUTH1 result payload differs")
    rows = payload.get("implementation_inventory")
    if not isinstance(rows, list) or len(rows) != len(IMPLEMENTATION_INVENTORY):
        raise TDG8RCV3ExecutionAuthorityError("RCV3 AUTH1 result inventory differs")
    for row, expected in zip(rows, IMPLEMENTATION_INVENTORY, strict=True):
        if (
            not isinstance(row, Mapping)
            or set(row) != {"role", "path", "sha256"}
            or (row.get("role"), row.get("path")) != expected
            or not isinstance(row.get("sha256"), str)
            or not _SHA256.fullmatch(str(row["sha256"]))
        ):
            raise TDG8RCV3ExecutionAuthorityError("RCV3 AUTH1 result inventory row differs")
    return result


def validate_compact(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    """Validate tracked AUTH1 bytes without reading Git, stores, or arrays."""
    config = _parse_config(config_raw)
    result = _parse_result(result_raw)
    payload = result["artifact_payload"]
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
            config, inventory, payload.get("environment", {}),
        ),
    }
    if canonical(result) != result_raw or canonical(expected) != result_raw:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 AUTH1 compact result differs")
    return result


def _require_git_image(root: Path, commit: str) -> None:
    if not _COMMIT.fullmatch(commit):
        raise TDG8RCV3ExecutionAuthorityError("RCV3 AUTH1 authority commit is malformed")
    head = rcv2.rcv1._git(root, "rev-parse", "--verify", "HEAD").stdout.decode("ascii").strip()
    resolved = rcv2.rcv1._git(
        root, "rev-parse", "--verify", f"{commit}^{{commit}}",
    ).stdout.decode("ascii").strip()
    if head != commit or resolved != commit:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 AUTH1 authority commit differs from HEAD")
    ancestor = subprocess.run(
        (
            "git", "--no-replace-objects", "--no-optional-locks", "-C", str(root),
            "merge-base", "--is-ancestor", PREF1_AUTHORITY_COMMIT, commit,
        ),
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=False,
        env={
            **{key: value for key, value in os.environ.items() if not key.startswith("GIT_")},
            "GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C",
        },
    )
    if ancestor.returncode != 0:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 PREF1 commit is not an ancestor")
    changed = rcv2.rcv1._git(
        root, "diff", "--name-only", "-z", "--no-renames",
        f"{PREF1_AUTHORITY_COMMIT}..{commit}", "--",
    ).stdout
    try:
        paths = tuple(sorted(item.decode("utf-8") for item in changed.split(b"\0") if item))
    except UnicodeDecodeError as exc:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 AUTH1 Git delta is not UTF-8") from exc
    if paths != AUTHORITY_DELTA_PATHS:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 AUTH1 authority Git delta differs")
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
        raise TDG8RCV3ExecutionAuthorityError("tracked tree differs from RCV3 AUTH1 commit")


def authorize_execution(root: Path, config_raw: bytes, result_raw: bytes,
                        authority_commit: str) -> TDG8RCV3ExecutionAuthority:
    repository = rcv2.rcv1._repository(root)
    result = validate_compact(config_raw, result_raw)
    _require_frozen_policy()
    _require_git_image(repository, authority_commit)
    for path, expected in ((CONFIG_PATH, config_raw), (RESULT_PATH, result_raw)):
        live = rcv2.rcv1._read_leaf(repository, path, f"RCV3 AUTH1 committed {path}", _MAX_RESULT_BYTES)
        if live != expected or rcv2.rcv1._committed_bytes(repository, authority_commit, path) != expected:
            raise TDG8RCV3ExecutionAuthorityError(f"RCV3 AUTH1 committed bytes differ: {path}")
    pref1_config, pref1_result = _validate_pref1_bytes(repository)
    for commit in (PREF1_AUTHORITY_COMMIT, authority_commit):
        if (
            rcv2.rcv1._committed_bytes(repository, commit, PREF1_CONFIG_PATH) != pref1_config
            or rcv2.rcv1._committed_bytes(repository, commit, PREF1_RESULT_PATH) != pref1_result
        ):
            raise TDG8RCV3ExecutionAuthorityError("RCV3 PREF1 lineage binding differs")
    rows: list[tuple[str, str, str]] = []
    for item in result["artifact_payload"]["implementation_inventory"]:
        role, path, digest = str(item["role"]), str(item["path"]), str(item["sha256"])
        live = rcv2.rcv1._read_leaf(repository, path, f"RCV3 AUTH1 implementation {path}", _MAX_BOUND_BYTES)
        if _sha(live) != digest or rcv2.rcv1._committed_bytes(repository, authority_commit, path) != live:
            raise TDG8RCV3ExecutionAuthorityError(f"RCV3 AUTH1 implementation bytes differ: {path}")
        rows.append((role, path, digest))
    environment = result["artifact_payload"].get("environment")
    if environment != rcv2.rcv1.observed_environment():
        raise TDG8RCV3ExecutionAuthorityError("RCV3 AUTH1 execution environment differs")
    return TDG8RCV3ExecutionAuthority(
        authority_commit=authority_commit,
        pref1_authority_commit=PREF1_AUTHORITY_COMMIT,
        pref1_result_sha256=PREF1_RESULT_SHA256,
        projection_id=PROJECTION_ID,
        receipt_sha256=RECEIPT_SHA256,
        original_execution_commit=ORIGINAL_EXECUTION_COMMIT,
        original_plan_sha256=ORIGINAL_PLAN_SHA256,
        config_sha256=_sha(config_raw),
        result_sha256=_sha(result_raw),
        campaign_id=CAMPAIGN_ID,
        destination_path=DESTINATION_PATH,
        entry_checkpoint_sha256=ENTRY_CHECKPOINT_SHA256,
        entry_journal_sha256=ENTRY_JOURNAL_SHA256,
        temporal_retry_cap=TEMPORAL_RETRY_CAP,
        minimum_macro_step_hex=MINIMUM_MACRO_STEP_HEX,
        implementation_inventory=tuple(rows),
        environment=MappingProxyType(dict(environment)),
    )


def _real_directory(repository: Path, relative: str, label: str) -> Path:
    current = repository
    for part in Path(relative).parts:
        current = current / part
        try:
            metadata = current.lstat()
        except OSError as exc:
            raise TDG8RCV3ExecutionAuthorityError(f"{label} is unavailable") from exc
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
            raise TDG8RCV3ExecutionAuthorityError(f"{label} is redirected or unsafe")
    return current


def _authenticate_external_receipt(repository: Path) -> None:
    raw = rcv2.rcv1._read_leaf(
        repository, RECEIPT_PATH, "RCV3 external receipt", _MAX_RESULT_BYTES,
    )
    if len(raw) != 3070 or _sha(raw) != RECEIPT_RAW_SHA256:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 external receipt raw identity differs")
    value = _decode_json(raw, "RCV3 external receipt")
    body = dict(value)
    observed = body.pop("receipt_sha256", None)
    if (
        observed != RECEIPT_SHA256
        or _sha(pref1._canonical(body)) != RECEIPT_SHA256
        or value != pref1._expected_receipt()
    ):
        raise TDG8RCV3ExecutionAuthorityError("RCV3 external receipt contract differs")


def _require_typed_receipt(receipt: TDG8RCV3ExecutionAuthority) -> None:
    if (
        not isinstance(receipt, TDG8RCV3ExecutionAuthority)
        or receipt.pref1_authority_commit != PREF1_AUTHORITY_COMMIT
        or receipt.pref1_result_sha256 != PREF1_RESULT_SHA256
        or receipt.projection_id != PROJECTION_ID
        or receipt.receipt_sha256 != RECEIPT_SHA256
        or receipt.original_execution_commit != ORIGINAL_EXECUTION_COMMIT
        or receipt.original_plan_sha256 != ORIGINAL_PLAN_SHA256
        or receipt.campaign_id != CAMPAIGN_ID
        or receipt.destination_path != DESTINATION_PATH
        or receipt.entry_checkpoint_sha256 != ENTRY_CHECKPOINT_SHA256
        or receipt.entry_journal_sha256 != ENTRY_JOURNAL_SHA256
        or receipt.temporal_retry_cap != TEMPORAL_RETRY_CAP
        or receipt.minimum_macro_step_hex != MINIMUM_MACRO_STEP_HEX
        or not receipt.candidate_branches_forbidden
    ):
        raise TDG8RCV3ExecutionAuthorityError("RCV3 AUTH1 typed execution authority differs")


def inspect_descendant(root: Path,
                       receipt: TDG8RCV3ExecutionAuthority) -> RCV3BoundedBoundary:
    """Authenticate the projected lineage without mutating it."""
    repository = rcv2.rcv1._repository(root)
    _require_typed_receipt(receipt)
    _require_frozen_policy()
    _real_directory(repository, DESTINATION_WRAPPER, "RCV3 projection wrapper")
    destination = _real_directory(repository, DESTINATION_PATH, "RCV3 campaign store")
    _authenticate_external_receipt(repository)
    store = HLT16CampaignStore(destination)
    try:
        snapshot = store.authenticated_snapshot()
        status = store.inspect_recovery()
        closure = store.authenticated_snapshot()
        closure_status = store.inspect_recovery()
    except HLT16CampaignStoreError as exc:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 descendant store does not validate") from exc
    if (
        closure.checkpoint.sha256 != snapshot.checkpoint.sha256
        or closure.suffix_kinds != snapshot.suffix_kinds
        or closure.suffix_records != snapshot.suffix_records
        or closure.staging_paths != snapshot.staging_paths
        or closure.orphan_descriptor_sha256 != snapshot.orphan_descriptor_sha256
        or closure.orphan_payload_semantic_sha256 != snapshot.orphan_payload_semantic_sha256
        or closure.terminal_lock_present != snapshot.terminal_lock_present
        or closure_status != status
    ):
        raise TDG8RCV3ExecutionAuthorityError(
            "RCV3 descendant changed during read-only inspection"
        )
    checkpoint = snapshot.checkpoint
    if (
        checkpoint.generation < ENTRY_CHECKPOINT_GENERATION
        or checkpoint.authorization_commit != ORIGINAL_EXECUTION_COMMIT
        or checkpoint.plan_sha256 != ORIGINAL_PLAN_SHA256
        or checkpoint.campaign_id != CAMPAIGN_ID
        or checkpoint.protocol != TARGET_PROTOCOL
        or tuple(checkpoint.members) != MEMBER_KEYS
    ):
        raise TDG8RCV3ExecutionAuthorityError("RCV3 descendant internal authority differs")
    try:
        rcv2._require_same_event(checkpoint)
    except Exception as exc:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 descendant escaped one-event scope") from exc
    if tuple(snapshot.suffix_kinds) not in rcv2._ALLOWED_SUFFIXES:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 descendant suffix differs")
    try:
        entry_ancestor = rcv2._require_entry_ancestors(store, checkpoint)
    except Exception as exc:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 generation-nine ancestry differs") from exc
    if not entry_ancestor:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 generation-nine ancestor is absent")
    try:
        entry = store.load_state(ENTRY_DESCRIPTOR_SHA256)
        if array_content_sha256(
            entry.arrays["u"], entry.arrays["p"], entry.arrays["q"],
        ) != ENTRY_PHYSICAL_STATE_SHA256:
            raise TDG8RCV3ExecutionAuthorityError("RCV3 entry physical ancestor differs")
        pending = rcv2._require_temporal_pending(store, checkpoint)
    except TDG8RCV3ExecutionAuthorityError:
        raise
    except Exception as exc:
        raise TDG8RCV3ExecutionAuthorityError("RCV3 descendant state cannot be restored") from exc
    if checkpoint.generation == ENTRY_CHECKPOINT_GENERATION and (
        checkpoint.sha256 != ENTRY_CHECKPOINT_SHA256
        or checkpoint.journal_sequence != ENTRY_JOURNAL_SEQUENCE
        or checkpoint.journal_tip_sha256 != ENTRY_JOURNAL_SHA256
        or snapshot.suffix_kinds
    ):
        raise TDG8RCV3ExecutionAuthorityError("RCV3 exact generation-nine boundary differs")
    terminal = status.terminal or checkpoint.disposition in {
        "scientific_terminal", "invalid_terminal",
    }
    safe = bool(
        status.safe_to_restart
        and not status.active_write
        and not snapshot.terminal_lock_present
        and checkpoint.event == EVENT
    )
    return RCV3BoundedBoundary(
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
        external_receipt_authenticated=True,
        entry_ancestor_present=True,
        temporal_pending_members=pending,
        terminal=terminal,
        terminal_lock_present=snapshot.terminal_lock_present,
        safe_to_restart=safe,
    )


__all__ = [
    "ARTIFACT_ID", "AUTHORITY_DELTA_PATHS", "CAMPAIGN_ID", "CLASSIFICATION",
    "CONFIG_PATH", "DESTINATION_PATH", "DESTINATION_WRAPPER", "EVENT",
    "IMPLEMENTATION_INVENTORY", "MINIMUM_MACRO_STEP_HEX",
    "ORIGINAL_EXECUTION_COMMIT", "ORIGINAL_PLAN_SHA256", "PREF1_AUTHORITY_COMMIT",
    "PREF1_RESULT_SHA256", "PROJECTION_ID", "RCV3BoundedBoundary",
    "RECEIPT_SHA256", "RESULT_PATH", "SUCCESSOR_EVENT", "SUCCESSOR_TARGET",
    "TARGET", "TDG8RCV3ExecutionAuthority", "TDG8RCV3ExecutionAuthorityError",
    "TEMPORAL_RETRY_CAP", "authorize_execution", "build_prelaunch", "canonical",
    "inspect_descendant", "validate_compact",
]
