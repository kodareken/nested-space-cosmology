"""Prospective authority for the exact TDG8 retry recovery and resume.

The first fresh TDG8 attempt is immutable evidence.  This module binds its
complete generation-seven checkpoint and sole sequence-seven TDG6 rejection,
reproduces the representation defect without mutating the campaign, and
authorizes only the unique metadata recovery plus continuation of the original
GR-0 event under a repaired committed image.

Compact verification is deliberately store-blind.  A mutating runner must
re-authenticate the live anchor immediately before taking the writer lease.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version
import json
import os
from pathlib import Path
import platform
import re
import stat
import subprocess
import sys
import tomllib
from types import MappingProxyType
from typing import Any, Iterable, Mapping

from . import proto15_runtime as p15
from . import hlt16_campaign_recovery as campaign_recovery
from .hlt16_campaign_recovery import preview_tdg6_recovery
from .hlt16_campaign_store import HLT16CampaignStore, HLT16CampaignStoreError
from .hlt16_state_store import HLT16StateStoreError, read_nofollow
from .numerical_engine import array_content_sha256


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG8-RCV1-AUTH1"
PROJECT_VERSION = "0.11.0"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CLASSIFICATION = (
    "premise_only_exact_retry_recovery_and_same_event_resume_authority"
)
CONFIG_PATH = "configs/fgc/fgc-1-tdg8-rcv1-auth1.toml"
RESULT_PATH = "results/fgc-1-tdg8-rcv1-auth1.json"

ORIGINAL_EXECUTION_COMMIT = "6df67967e6e0bc63eca4dc6951a60d3787368f26"
ORIGINAL_PLAN_SHA256 = (
    "e649e477b98411a3051f9946f2a9248b875d2b598818b8983ce06caa969d2126"
)
ORIGINAL_AUTHORITY_RESULT_PATH = "results/fgc-1-tdg8-run1-auth1.json"
ORIGINAL_AUTHORITY_RESULT_SHA256 = (
    "6ee22000bb76ac9f69ed3fe1056a7be319d24dcbfe37b6160eb8e9d0e3ec40de"
)
DESTINATION_PATH = "runs/fgc-2-sf1/proto19/calibration"
CAMPAIGN_ID = "FGC-2-SF1-PROTO18-GR0-A3-TDG8-SUCCESSOR-1"
MEMBER_KEY = "RK4-2049"

ANCHOR_CHECKPOINT_GENERATION = 7
ANCHOR_CHECKPOINT_SHA256 = (
    "51aec0cdc38d8fbfa6c5f81dc4f2f20e7332cc392a8621c2126c706fd9e3e3fd"
)
ANCHOR_JOURNAL_SEQUENCE = 7
ANCHOR_JOURNAL_SHA256 = (
    "6e8b55b6f9a4920046f15f3591968144664d0647aecfd4b962d85ad0e96ff981"
)
ANCHOR_CHECKPOINT_JOURNAL_TIP_SHA256 = (
    "ec4fd6e65657f1a7b354bdc6c2d90392bffb67f69d78c1161a440b3eeabd19db"
)
ANCHOR_DESCRIPTOR_SHA256 = (
    "77847126340e78c8ac300fac2795bceda724c15e0894e34444b0da42a84166b4"
)
ANCHOR_PHYSICAL_STATE_SHA256 = (
    "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a"
)
ANCHOR_ACCEPTED_TIME = {
    "rational": "206891375579527/140737488355328",
    "binary64_hex": "0x1.78554de5a30e0p+0",
}

EXPECTED_CHECKPOINT_GENERATION = 8
EXPECTED_CHECKPOINT_SHA256 = (
    "a6401ff970aea002f8991c745a4cf7fe91c6d0a556907b04280d87e2c04d8936"
)
EXPECTED_JOURNAL_SEQUENCE = 8
EXPECTED_JOURNAL_SHA256 = (
    "63af61720215223eb3a660b055a07fa897a317b9580128d1b8f79a7022bab86f"
)
EXPECTED_PENDING_CAP_HEX = "0x1.aaa9612df9000p-10"
FROZEN_MINIMUM_WIDTH_HEX = "0x1.0000000000000p-30"

_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")
_MAX_CONFIG_BYTES = 1024 * 1024
_MAX_RESULT_BYTES = 4 * 1024 * 1024
_MAX_BOUND_BYTES = 32 * 1024 * 1024

IMPLEMENTATION_INVENTORY = (
    ("authority", "src/recursive_horizons/fgc/evolution/tdg8_retry_recovery_authority.py"),
    ("runner", "scripts/run_fgc_tdg8_retry_recovery.py"),
    ("reproducer", "scripts/reproduce_fgc_tdg8_rcv1_auth1.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_campaign_runtime.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_campaign_recovery.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_campaign_store.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_lifecycle.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_progression_attempt.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/tdg8_successor_runtime.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/proto19_progression_contract.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/proto19_gr0_static_factory.py"),
    ("orchestrator", "scripts/run_fgc_pro19_event1.py"),
)

# Every other execution dependency is inherited byte-for-byte from the
# original committed image.  The result is materialized only after the single
# campaign-runtime repair has stabilized.
AUTHORITY_DELTA_PATHS = tuple(sorted((
    "Makefile",
    CONFIG_PATH,
    RESULT_PATH,
    "scripts/check_repo.py",
    "scripts/reproduce_fgc_tdg8_rcv1_auth1.py",
    "scripts/run_fgc_tdg8_retry_recovery.py",
    "src/recursive_horizons/fgc/evolution/hlt16_campaign_runtime.py",
    "src/recursive_horizons/fgc/evolution/tdg8_retry_recovery_authority.py",
    "tests/test_fgc_hlt16_campaign_runtime.py",
    "tests/test_fgc_tdg8_retry_recovery_authority.py",
    "tests/test_fgc_tdg8_retry_recovery_runner.py",
)))

SCOPE = {
    "status_read_only": True,
    "recovery_metadata_only": True,
    "same_event_resume_only": True,
    "destination_absence_recheck": False,
    "generation_zero_reimport": False,
    "accepted_step_replay": False,
    "threshold_or_retry_rule_change": False,
    "candidate_branches_forbidden": True,
}
CLAIMS = {
    "original_attempt_preserved": True,
    "physical_state_advanced_by_recovery": False,
    "common_event_completed": False,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "physical_result_earned": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
}
NONCLAIMS = (
    "This authority preserves the original append-only TDG8 attempt and does not reclassify its numerical stop as physics.",
    "Recovery may publish only the uniquely derived metadata successor; it cannot advance or replace an accepted state.",
    "Resume remains GR-0-only, uses the original event plan, and opens no SGB-L or FGC-QR branch.",
    "This authority records no common-event, calibration, trapping, DEF1, retained-EFT, transition, mechanism, or physical result.",
)


class TDG8RetryRecoveryAuthorityError(RuntimeError):
    """The compact authority, committed image, or live anchor differs."""


@dataclass(frozen=True, slots=True)
class TDG8RetryRecoveryExecutionAuthority:
    authority_commit: str
    original_execution_commit: str
    original_plan_sha256: str
    config_sha256: str
    result_sha256: str
    campaign_id: str
    destination_path: str
    anchor_checkpoint_sha256: str
    anchor_journal_sha256: str
    expected_checkpoint_sha256: str
    expected_journal_sha256: str
    implementation_inventory: tuple[tuple[str, str, str], ...]
    environment: Mapping[str, str]
    candidate_branches_forbidden: bool = True


def canonical(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _safe_relative(value: object, label: str) -> str:
    if not isinstance(value, str):
        raise TDG8RetryRecoveryAuthorityError(f"{label} is not a path")
    path = Path(value)
    if path.is_absolute() or not path.parts or any(
        part in {"", ".", ".."} for part in path.parts
    ) or path.as_posix() != value:
        raise TDG8RetryRecoveryAuthorityError(f"{label} is unsafe")
    return value


def _read_leaf(root: Path, relative: str, label: str, maximum: int) -> bytes:
    try:
        raw = read_nofollow(root, _safe_relative(relative, label), label)
    except HLT16StateStoreError as exc:
        raise TDG8RetryRecoveryAuthorityError(f"{label} cannot be read safely") from exc
    if len(raw) > maximum:
        raise TDG8RetryRecoveryAuthorityError(f"{label} exceeds its size bound")
    return raw


def _decode_json(raw: bytes, label: str) -> dict[str, Any]:
    def pairs(items: Iterable[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in items:
            if key in result:
                raise ValueError(key)
            result[key] = value
        return result

    try:
        value = json.loads(raw.decode("utf-8"), object_pairs_hook=pairs)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TDG8RetryRecoveryAuthorityError(f"{label} is malformed") from exc
    if not isinstance(value, dict):
        raise TDG8RetryRecoveryAuthorityError(f"{label} is not object-valued")
    return value


def _git(root: Path, *arguments: str) -> subprocess.CompletedProcess[bytes]:
    environment = {
        key: value for key, value in os.environ.items() if not key.startswith("GIT_")
    }
    environment.update({"GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C"})
    try:
        result = subprocess.run(
            (
                "git", "--no-replace-objects", "--no-optional-locks",
                "-c", "core.fsmonitor=false", "-c", "core.untrackedCache=false",
                "-C", str(root), *arguments,
            ),
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False, env=environment,
        )
    except OSError as exc:
        raise TDG8RetryRecoveryAuthorityError("Git query could not start") from exc
    if result.returncode != 0:
        raise TDG8RetryRecoveryAuthorityError("Git authority query failed")
    return result


def _repository(root: Path) -> Path:
    repository = Path(os.path.abspath(os.fspath(root)))
    try:
        metadata = repository.lstat()
    except OSError as exc:
        raise TDG8RetryRecoveryAuthorityError("repository root is absent") from exc
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        raise TDG8RetryRecoveryAuthorityError("repository root is unsafe")
    return repository


def _parse_config(raw: bytes) -> dict[str, Any]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise TDG8RetryRecoveryAuthorityError("RCV1 config is malformed") from exc
    if not isinstance(value, dict):
        raise TDG8RetryRecoveryAuthorityError("RCV1 config is not object-valued")
    expected_top = {
        "schema_version", "artifact_id", "project_version", "target_protocol",
        "classification", "nonclaims", "predecessor", "live_anchor",
        "diagnosis", "expected_recovery", "implementation_inventory", "scope",
        "claims",
    }
    if set(value) != expected_top:
        raise TDG8RetryRecoveryAuthorityError("RCV1 config keys differ")
    if (
        value["schema_version"] != SCHEMA_VERSION
        or value["artifact_id"] != ARTIFACT_ID
        or value["project_version"] != PROJECT_VERSION
        or value["target_protocol"] != TARGET_PROTOCOL
        or value["classification"] != CLASSIFICATION
        or value["nonclaims"] != list(NONCLAIMS)
        or value["scope"] != SCOPE
        or value["claims"] != CLAIMS
    ):
        raise TDG8RetryRecoveryAuthorityError("RCV1 config contract differs")
    predecessor = value["predecessor"]
    anchor = value["live_anchor"]
    diagnosis = value["diagnosis"]
    recovery = value["expected_recovery"]
    if predecessor != {
        "original_execution_commit": ORIGINAL_EXECUTION_COMMIT,
        "original_plan_sha256": ORIGINAL_PLAN_SHA256,
        "original_authority_result_path": ORIGINAL_AUTHORITY_RESULT_PATH,
        "original_authority_result_sha256": ORIGINAL_AUTHORITY_RESULT_SHA256,
    }:
        raise TDG8RetryRecoveryAuthorityError("RCV1 predecessor differs")
    if anchor != {
        "destination_path": DESTINATION_PATH, "campaign_id": CAMPAIGN_ID,
        "branch": "GR-0", "amplitude": "3", "event": 23,
        "target_rational": "3/2",
        "checkpoint_generation": ANCHOR_CHECKPOINT_GENERATION,
        "checkpoint_sha256": ANCHOR_CHECKPOINT_SHA256,
        "journal_sequence": ANCHOR_JOURNAL_SEQUENCE,
        "journal_sha256": ANCHOR_JOURNAL_SHA256,
        "member_key": MEMBER_KEY,
        "physical_state_sha256": ANCHOR_PHYSICAL_STATE_SHA256,
        "writer_state": "stale_verified_lock",
    }:
        raise TDG8RetryRecoveryAuthorityError("RCV1 live anchor differs")
    if diagnosis != {
        "classification": "omitted_frozen_minimum_width_metadata_only",
        "field": "minimum_width_hex_or_none",
        "observed_buggy_value": "none",
        "required_value": FROZEN_MINIMUM_WIDTH_HEX,
        "all_other_plan_fields_equal": True,
        "threshold_or_rule_change": False,
    }:
        raise TDG8RetryRecoveryAuthorityError("RCV1 diagnosis differs")
    if recovery != {
        "checkpoint_generation": EXPECTED_CHECKPOINT_GENERATION,
        "checkpoint_sha256": EXPECTED_CHECKPOINT_SHA256,
        "journal_sequence": EXPECTED_JOURNAL_SEQUENCE,
        "journal_sha256": EXPECTED_JOURNAL_SHA256,
        "member_mode": "RETRY_PENDING", "pending_owner": "temporal",
        "pending_cap_hex": EXPECTED_PENDING_CAP_HEX,
        "physical_state_advanced": False,
    }:
        raise TDG8RetryRecoveryAuthorityError("RCV1 recovery contract differs")
    rows = value["implementation_inventory"]
    if not isinstance(rows, list) or tuple(
        (row.get("role"), row.get("path")) if isinstance(row, Mapping) else (None, None)
        for row in rows
    ) != IMPLEMENTATION_INVENTORY or any(
        not isinstance(row, Mapping) or set(row) != {"role", "path"}
        for row in rows
    ):
        raise TDG8RetryRecoveryAuthorityError("RCV1 implementation inventory differs")
    for _role, path in IMPLEMENTATION_INVENTORY:
        _safe_relative(path, "RCV1 implementation path")
    return value


def observed_environment() -> dict[str, str]:
    try:
        numpy_version = version("numpy")
    except PackageNotFoundError as exc:
        raise TDG8RetryRecoveryAuthorityError("NumPy is absent") from exc
    return {
        "implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "numpy_version": numpy_version,
        "system": platform.system(),
        "machine": platform.machine(),
        "platform": sys.platform,
        "byteorder": sys.byteorder,
    }


def _inventory(root: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for role, path in IMPLEMENTATION_INVENTORY:
        raw = _read_leaf(root, path, f"RCV1 implementation {path}", _MAX_BOUND_BYTES)
        rows.append({"role": role, "path": path, "sha256": _sha(raw)})
    return rows


def _diagnose(record: Mapping[str, Any], target_hex: str) -> dict[str, Any]:
    try:
        evidence = record["payload"]["evidence"]
        initial = float(evidence["initial_time"])
        retry = float(evidence["retry_step_size"])
        target = float.fromhex(target_hex)
        buggy = p15._plan_mapping(
            p15.plan_forward_proto14_subdivision(initial, target, retry)
        )
        correct = p15._plan_mapping(
            p15.plan_forward_proto14_subdivision(
                initial, target, retry,
                minimum_width=float(p15.TDG6_MINIMUM_MACRO_STEP),
            )
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise TDG8RetryRecoveryAuthorityError("RCV1 retry evidence cannot derive") from exc
    differences = {
        key: {"buggy": buggy[key], "required": correct[key]}
        for key in buggy if buggy[key] != correct[key]
    }
    expected = {
        "minimum_width_hex_or_none": {
            "buggy": None, "required": FROZEN_MINIMUM_WIDTH_HEX,
        }
    }
    if differences != expected:
        raise TDG8RetryRecoveryAuthorityError("RCV1 defect is not minimum metadata only")
    return {
        "classification": "omitted_frozen_minimum_width_metadata_only",
        "differing_fields": differences,
        "all_other_plan_fields_equal": True,
        "numerical_boundaries_equal": buggy["boundaries_hex"] == correct["boundaries_hex"],
        "threshold_or_rule_change": False,
    }


def inspect_live_anchor(
    root: Path,
    *,
    permitted_writer_states: frozenset[str] = frozenset({"stale_verified_lock"}),
) -> dict[str, Any]:
    """Read and derive the unique live edge without acquiring a lease."""

    repository = _repository(root)
    destination = repository / DESTINATION_PATH
    store = HLT16CampaignStore(destination)
    try:
        snapshot = store.authenticated_snapshot()
        status = store.inspect_recovery()
        preview = preview_tdg6_recovery(store)
    except (HLT16CampaignStoreError, ValueError) as exc:
        raise TDG8RetryRecoveryAuthorityError("RCV1 live store does not validate") from exc
    checkpoint = snapshot.checkpoint
    if (
        checkpoint.authorization_commit != ORIGINAL_EXECUTION_COMMIT
        or checkpoint.plan_sha256 != ORIGINAL_PLAN_SHA256
        or checkpoint.campaign_id != CAMPAIGN_ID
        or checkpoint.protocol != TARGET_PROTOCOL
        or checkpoint.event != 23
        or checkpoint.target != {"rational": "3/2", "binary64_hex": "0x1.8000000000000p+0"}
        or checkpoint.generation != ANCHOR_CHECKPOINT_GENERATION
        or checkpoint.sha256 != ANCHOR_CHECKPOINT_SHA256
        or checkpoint.journal_sequence != 6
        or snapshot.suffix_kinds != ("tdg6_rejection",)
        or len(snapshot.suffix_records) != 1
        or snapshot.staging_paths
        or snapshot.orphan_descriptor_sha256 is not None
        or snapshot.orphan_payload_semantic_sha256 is not None
        or snapshot.terminal_lock_present
        or status.state != "recognized_uncheckpointed_journal_suffix"
        or not status.active_write
        or status.writer_state not in permitted_writer_states
        or status.safe_to_restart
    ):
        raise TDG8RetryRecoveryAuthorityError("RCV1 generation-seven boundary differs")
    record = snapshot.suffix_records[0]
    if (
        record.get("sequence") != ANCHOR_JOURNAL_SEQUENCE
        or record.get("record_sha256") != ANCHOR_JOURNAL_SHA256
        or record.get("campaign_id") != CAMPAIGN_ID
        or record.get("generation") != EXPECTED_CHECKPOINT_GENERATION
        or record.get("event") != 23
        or record.get("kind") != "tdg6_rejection"
        or record.get("payload", {}).get("member_key") != MEMBER_KEY
    ):
        raise TDG8RetryRecoveryAuthorityError("RCV1 rejection anchor differs")
    member = checkpoint.members[MEMBER_KEY]
    if preview.evolution_state_sha256 != ANCHOR_PHYSICAL_STATE_SHA256:
        raise TDG8RetryRecoveryAuthorityError("RCV1 physical-state identity differs")
    recovered = preview.checkpoint
    recovered_member = recovered.members[MEMBER_KEY]
    if (
        preview.predecessor_checkpoint_sha256 != ANCHOR_CHECKPOINT_SHA256
        or preview.physical_state_advanced
        or recovered.generation != EXPECTED_CHECKPOINT_GENERATION
        or recovered.sha256 != EXPECTED_CHECKPOINT_SHA256
        or recovered.journal_sequence != EXPECTED_JOURNAL_SEQUENCE
        or recovered.journal_tip_sha256 != EXPECTED_JOURNAL_SHA256
        or tuple((item["kind"], item["sequence"], item["record_sha256"])
                 for item in preview.records) != (
            ("tdg6_rejection", 7, ANCHOR_JOURNAL_SHA256),
            ("cursor_transition", 8, EXPECTED_JOURNAL_SHA256),
        )
        or recovered_member.descriptor_sha256 != member.descriptor_sha256
        or recovered_member.cursor["accepted_boundary_time"] != member.cursor["accepted_boundary_time"]
        or recovered_member.cursor["mode"] != "RETRY_PENDING"
        or recovered_member.pending_owner != "temporal"
        or recovered_member.pending_cap_hex != EXPECTED_PENDING_CAP_HEX
        or recovered_member.cfl_current != 1
        or recovered_member.cfl_total != 1
        or recovered_member.ledger["cumulative_temporal_retry_count"] != 1
    ):
        raise TDG8RetryRecoveryAuthorityError("RCV1 unique recovery preview differs")
    diagnosis = _diagnose(record, checkpoint.target["binary64_hex"])
    return {
        "anchor": {
            "authorization_commit": checkpoint.authorization_commit,
            "plan_sha256": checkpoint.plan_sha256,
            "checkpoint_generation": checkpoint.generation,
            "checkpoint_sha256": checkpoint.sha256,
            "checkpoint_journal_sequence": checkpoint.journal_sequence,
            "checkpoint_journal_tip_sha256": checkpoint.journal_tip_sha256,
            "suffix_sequence": record["sequence"],
            "suffix_sha256": record["record_sha256"],
            "member_key": MEMBER_KEY,
            "descriptor_sha256": member.descriptor_sha256,
            "physical_state_sha256": preview.evolution_state_sha256,
            "accepted_time": dict(member.cursor["accepted_boundary_time"]),
            "writer_state": status.writer_state,
            "physical_state_advanced": False,
        },
        "diagnosis": diagnosis,
        "expected_recovery": {
            "checkpoint_generation": recovered.generation,
            "checkpoint_sha256": recovered.sha256,
            "journal_sequence": recovered.journal_sequence,
            "journal_sha256": recovered.journal_tip_sha256,
            "member_mode": recovered_member.cursor["mode"],
            "pending_owner": recovered_member.pending_owner,
            "pending_cap_hex": recovered_member.pending_cap_hex,
            "descriptor_sha256": recovered_member.descriptor_sha256,
            "physical_state_sha256": preview.evolution_state_sha256,
            "physical_state_advanced": False,
        },
    }


def inspect_recovery_boundary(
    root: Path,
    *,
    permitted_writer_states: frozenset[str] = frozenset({"stale_verified_lock"}),
) -> dict[str, Any]:
    """Classify the three exact recoverable publication cuts.

    This reader admits only the original rejection suffix, that suffix plus
    its byte-exact cursor transition, or the resulting clean generation-eight
    checkpoint.  It never publishes, removes, adopts, or replays anything.
    """

    repository = _repository(root)
    store = HLT16CampaignStore(repository / DESTINATION_PATH)
    try:
        snapshot = store.authenticated_snapshot()
        status = store.inspect_recovery()
    except HLT16CampaignStoreError as exc:
        raise TDG8RetryRecoveryAuthorityError(
            "RCV1 recovery boundary does not validate"
        ) from exc
    checkpoint = snapshot.checkpoint
    if (
        checkpoint.authorization_commit != ORIGINAL_EXECUTION_COMMIT
        or checkpoint.plan_sha256 != ORIGINAL_PLAN_SHA256
        or checkpoint.campaign_id != CAMPAIGN_ID
        or checkpoint.protocol != TARGET_PROTOCOL
        or checkpoint.event != 23
        or checkpoint.target
        != {"rational": "3/2", "binary64_hex": "0x1.8000000000000p+0"}
        or snapshot.staging_paths
        or snapshot.orphan_descriptor_sha256 is not None
        or snapshot.orphan_payload_semantic_sha256 is not None
        or snapshot.terminal_lock_present
        or not status.active_write
        or status.writer_state not in permitted_writer_states
        or status.safe_to_restart
    ):
        raise TDG8RetryRecoveryAuthorityError(
            "RCV1 recoverable campaign boundary differs"
        )

    if (
        checkpoint.generation == ANCHOR_CHECKPOINT_GENERATION
        and checkpoint.sha256 == ANCHOR_CHECKPOINT_SHA256
        and snapshot.suffix_kinds == ("tdg6_rejection",)
    ):
        live = inspect_live_anchor(
            repository, permitted_writer_states=permitted_writer_states
        )
        return {
            "state": "generation7_rejection_only",
            "checkpoint_sha256": checkpoint.sha256,
            "journal_tip_sha256": live["anchor"]["suffix_sha256"],
            "physical_state_sha256": live["anchor"]["physical_state_sha256"],
        }

    if (
        checkpoint.generation == ANCHOR_CHECKPOINT_GENERATION
        and checkpoint.sha256 == ANCHOR_CHECKPOINT_SHA256
        and snapshot.suffix_kinds == ("tdg6_rejection", "cursor_transition")
        and len(snapshot.suffix_records) == 2
    ):
        rejection, transition = snapshot.suffix_records
        member = checkpoint.members[MEMBER_KEY]
        if (
            checkpoint.journal_sequence != 6
            or checkpoint.journal_tip_sha256
            != ANCHOR_CHECKPOINT_JOURNAL_TIP_SHA256
            or member.descriptor_sha256 != ANCHOR_DESCRIPTOR_SHA256
            or member.cursor["accepted_boundary_time"] != ANCHOR_ACCEPTED_TIME
            or rejection.get("sequence") != ANCHOR_JOURNAL_SEQUENCE
            or rejection.get("record_sha256") != ANCHOR_JOURNAL_SHA256
            or rejection.get("kind") != "tdg6_rejection"
            or rejection.get("payload", {}).get("member_key") != MEMBER_KEY
            or transition.get("sequence") != EXPECTED_JOURNAL_SEQUENCE
            or transition.get("record_sha256") != EXPECTED_JOURNAL_SHA256
            or transition.get("kind") != "cursor_transition"
            or transition.get("previous_record_sha256") != ANCHOR_JOURNAL_SHA256
            or transition.get("payload", {}).get("member_key") != MEMBER_KEY
        ):
            raise TDG8RetryRecoveryAuthorityError(
                "RCV1 intermediate transition identity differs"
            )
        try:
            persisted = store.load_state(member.descriptor_sha256)
            preview = campaign_recovery._resolve_tdg6_rejection(
                checkpoint, rejection, persisted.arrays
            )
        except Exception as exc:
            raise TDG8RetryRecoveryAuthorityError(
                "RCV1 intermediate transition cannot derive"
            ) from exc
        if (
            tuple(dict(item) for item in preview.records)
            != tuple(dict(item) for item in snapshot.suffix_records)
            or preview.predecessor_checkpoint_sha256 != ANCHOR_CHECKPOINT_SHA256
            or preview.evolution_state_sha256 != ANCHOR_PHYSICAL_STATE_SHA256
            or preview.physical_state_advanced
            or preview.checkpoint.generation != EXPECTED_CHECKPOINT_GENERATION
            or preview.checkpoint.sha256 != EXPECTED_CHECKPOINT_SHA256
            or preview.checkpoint.journal_sequence != EXPECTED_JOURNAL_SEQUENCE
            or preview.checkpoint.journal_tip_sha256 != EXPECTED_JOURNAL_SHA256
        ):
            raise TDG8RetryRecoveryAuthorityError(
                "RCV1 intermediate transition does not imply the frozen checkpoint"
            )
        return {
            "state": "generation7_transition_published",
            "checkpoint_sha256": checkpoint.sha256,
            "journal_tip_sha256": EXPECTED_JOURNAL_SHA256,
            "physical_state_sha256": preview.evolution_state_sha256,
        }

    if (
        checkpoint.generation == EXPECTED_CHECKPOINT_GENERATION
        and checkpoint.sha256 == EXPECTED_CHECKPOINT_SHA256
        and checkpoint.parent_sha256 == ANCHOR_CHECKPOINT_SHA256
        and checkpoint.journal_sequence == EXPECTED_JOURNAL_SEQUENCE
        and checkpoint.journal_tip_sha256 == EXPECTED_JOURNAL_SHA256
        and checkpoint.disposition == "nonterminal"
        and snapshot.suffix_classification == "clean"
        and not snapshot.suffix_records
    ):
        member = checkpoint.members[MEMBER_KEY]
        try:
            persisted = store.load_state(member.descriptor_sha256)
            physical = array_content_sha256(
                persisted.arrays["u"], persisted.arrays["p"], persisted.arrays["q"]
            )
        except Exception as exc:
            raise TDG8RetryRecoveryAuthorityError(
                "RCV1 recovered physical state cannot be restored"
            ) from exc
        if (
            member.descriptor_sha256 != ANCHOR_DESCRIPTOR_SHA256
            or member.cursor["accepted_boundary_time"] != ANCHOR_ACCEPTED_TIME
            or member.cursor["mode"] != "RETRY_PENDING"
            or member.pending_owner != "temporal"
            or member.pending_cap_hex != EXPECTED_PENDING_CAP_HEX
            or member.cfl_current != 1
            or member.cfl_total != 1
            or member.ledger["cumulative_temporal_retry_count"] != 1
            or physical != ANCHOR_PHYSICAL_STATE_SHA256
        ):
            raise TDG8RetryRecoveryAuthorityError(
                "RCV1 recovered generation-eight state differs"
            )
        return {
            "state": "generation8_checkpoint_published",
            "checkpoint_sha256": checkpoint.sha256,
            "journal_tip_sha256": checkpoint.journal_tip_sha256,
            "physical_state_sha256": physical,
        }

    raise TDG8RetryRecoveryAuthorityError(
        "RCV1 campaign is outside the three recoverable publication cuts"
    )


def _payload(config: Mapping[str, Any], inventory: list[dict[str, str]],
             live: Mapping[str, Any], environment: Mapping[str, str]) -> dict[str, Any]:
    return {
        "predecessor": dict(config["predecessor"]),
        "live_anchor": dict(live["anchor"]),
        "diagnosis": dict(live["diagnosis"]),
        "expected_recovery": dict(live["expected_recovery"]),
        "implementation_inventory": inventory,
        "environment": dict(environment),
        "scope": dict(SCOPE),
        "claims": dict(CLAIMS),
        "nonclaims": list(NONCLAIMS),
    }


def build_prelaunch(config_raw: bytes, root: Path) -> dict[str, Any]:
    repository = _repository(root)
    config = _parse_config(config_raw)
    if _read_leaf(repository, CONFIG_PATH, "RCV1 config", _MAX_CONFIG_BYTES) != config_raw:
        raise TDG8RetryRecoveryAuthorityError("live RCV1 config bytes differ")
    result_path = repository / RESULT_PATH
    try:
        result_path.lstat()
    except FileNotFoundError:
        pass
    else:
        raise TDG8RetryRecoveryAuthorityError("RCV1 result already exists")
    original = _read_leaf(
        repository, ORIGINAL_AUTHORITY_RESULT_PATH, "original TDG8 authority result",
        _MAX_RESULT_BYTES,
    )
    if _sha(original) != ORIGINAL_AUTHORITY_RESULT_SHA256:
        raise TDG8RetryRecoveryAuthorityError("original TDG8 authority result changed")
    live = inspect_live_anchor(repository)
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL,
        "classification": CLASSIFICATION,
        "gate_status": "pass",
        "source_config_sha256": _sha(config_raw),
        "artifact_payload": _payload(
            config, _inventory(repository), live, observed_environment()
        ),
    }


def _parse_result(raw: bytes) -> dict[str, Any]:
    result = _decode_json(raw, "RCV1 result")
    payload = result.get("artifact_payload")
    if not isinstance(payload, Mapping):
        raise TDG8RetryRecoveryAuthorityError("RCV1 result payload differs")
    rows = payload.get("implementation_inventory")
    if not isinstance(rows, list) or len(rows) != len(IMPLEMENTATION_INVENTORY):
        raise TDG8RetryRecoveryAuthorityError("RCV1 result inventory differs")
    for row, expected in zip(rows, IMPLEMENTATION_INVENTORY, strict=True):
        if (
            not isinstance(row, Mapping)
            or set(row) != {"role", "path", "sha256"}
            or (row.get("role"), row.get("path")) != expected
            or not isinstance(row.get("sha256"), str)
            or not _SHA256.fullmatch(str(row["sha256"]))
        ):
            raise TDG8RetryRecoveryAuthorityError("RCV1 result inventory row differs")
    return result


def _fixed_live_contract() -> dict[str, Any]:
    """Return the store facts frozen before result publication."""

    return {
        "anchor": {
            "authorization_commit": ORIGINAL_EXECUTION_COMMIT,
            "plan_sha256": ORIGINAL_PLAN_SHA256,
            "checkpoint_generation": ANCHOR_CHECKPOINT_GENERATION,
            "checkpoint_sha256": ANCHOR_CHECKPOINT_SHA256,
            "checkpoint_journal_sequence": 6,
            "checkpoint_journal_tip_sha256": ANCHOR_CHECKPOINT_JOURNAL_TIP_SHA256,
            "suffix_sequence": ANCHOR_JOURNAL_SEQUENCE,
            "suffix_sha256": ANCHOR_JOURNAL_SHA256,
            "member_key": MEMBER_KEY,
            "descriptor_sha256": ANCHOR_DESCRIPTOR_SHA256,
            "physical_state_sha256": ANCHOR_PHYSICAL_STATE_SHA256,
            "accepted_time": dict(ANCHOR_ACCEPTED_TIME),
            "writer_state": "stale_verified_lock",
            "physical_state_advanced": False,
        },
        "diagnosis": {
            "classification": "omitted_frozen_minimum_width_metadata_only",
            "differing_fields": {"minimum_width_hex_or_none": {"buggy": None, "required": FROZEN_MINIMUM_WIDTH_HEX}},
            "all_other_plan_fields_equal": True,
            "numerical_boundaries_equal": True,
            "threshold_or_rule_change": False,
        },
        "expected_recovery": {
            "checkpoint_generation": EXPECTED_CHECKPOINT_GENERATION,
            "checkpoint_sha256": EXPECTED_CHECKPOINT_SHA256,
            "journal_sequence": EXPECTED_JOURNAL_SEQUENCE,
            "journal_sha256": EXPECTED_JOURNAL_SHA256,
            "member_mode": "RETRY_PENDING", "pending_owner": "temporal",
            "pending_cap_hex": EXPECTED_PENDING_CAP_HEX,
            "descriptor_sha256": ANCHOR_DESCRIPTOR_SHA256,
            "physical_state_sha256": ANCHOR_PHYSICAL_STATE_SHA256,
            "physical_state_advanced": False,
        },
    }


def validate_compact(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    config = _parse_config(config_raw)
    result = _parse_result(result_raw)
    payload = result["artifact_payload"]
    inventory = [dict(item) for item in payload["implementation_inventory"]]
    expected_live = {
        "anchor": payload["live_anchor"],
        "diagnosis": payload["diagnosis"],
        "expected_recovery": payload["expected_recovery"],
    }
    expected = {
        "schema_version": SCHEMA_VERSION, "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION, "target_protocol": TARGET_PROTOCOL,
        "classification": CLASSIFICATION, "gate_status": "pass",
        "source_config_sha256": _sha(config_raw),
        "artifact_payload": _payload(
            config, inventory, expected_live, payload.get("environment", {})
        ),
    }
    fixed_live = _fixed_live_contract()
    fixed_payload = {
        "live_anchor": fixed_live["anchor"],
        "diagnosis": fixed_live["diagnosis"],
        "expected_recovery": fixed_live["expected_recovery"],
    }
    for name, value in fixed_payload.items():
        if payload.get(name) != value:
            raise TDG8RetryRecoveryAuthorityError(f"RCV1 compact {name} differs")
    if canonical(result) != result_raw or canonical(expected) != result_raw:
        raise TDG8RetryRecoveryAuthorityError("RCV1 compact result differs")
    return result


def _require_git_image(root: Path, commit: str) -> None:
    if not _COMMIT.fullmatch(commit):
        raise TDG8RetryRecoveryAuthorityError("authority commit is malformed")
    head = _git(root, "rev-parse", "--verify", "HEAD").stdout.decode("ascii").strip()
    resolved = _git(root, "rev-parse", "--verify", f"{commit}^{{commit}}").stdout.decode("ascii").strip()
    if head != commit or resolved != commit:
        raise TDG8RetryRecoveryAuthorityError("authority commit differs from HEAD")
    ancestor = subprocess.run(
        (
            "git", "--no-replace-objects", "--no-optional-locks", "-C", str(root),
            "merge-base", "--is-ancestor", ORIGINAL_EXECUTION_COMMIT, commit,
        ),
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=False, env={**{k: v for k, v in os.environ.items() if not k.startswith("GIT_")}, "GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C"},
    )
    if ancestor.returncode != 0:
        raise TDG8RetryRecoveryAuthorityError("original execution commit is not an ancestor")
    changed = _git(
        root, "diff", "--name-only", "-z", "--no-renames",
        f"{ORIGINAL_EXECUTION_COMMIT}..{commit}", "--",
    ).stdout
    try:
        paths = tuple(sorted(item.decode("utf-8") for item in changed.split(b"\0") if item))
    except UnicodeDecodeError as exc:
        raise TDG8RetryRecoveryAuthorityError("authority delta is not UTF-8") from exc
    if paths != AUTHORITY_DELTA_PATHS:
        raise TDG8RetryRecoveryAuthorityError("authority Git delta differs")
    dirty = subprocess.run(
        (
            "git", "--no-replace-objects", "--no-optional-locks",
            "-c", "core.fsmonitor=false", "-c", "core.untrackedCache=false",
            "-C", str(root), "diff", "--quiet", "--no-ext-diff",
            "--ignore-submodules=all", "HEAD", "--",
        ),
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=False, env={**{k: v for k, v in os.environ.items() if not k.startswith("GIT_")}, "GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C"},
    )
    if dirty.returncode != 0:
        raise TDG8RetryRecoveryAuthorityError("tracked tree differs from authority commit")


def _committed_bytes(root: Path, commit: str, path: str) -> bytes:
    tracked = _git(root, "ls-tree", "-r", "--name-only", commit, "--", path).stdout.decode("utf-8").splitlines()
    if tracked != [path]:
        raise TDG8RetryRecoveryAuthorityError(f"authority path is not tracked: {path}")
    return _git(root, "show", f"{commit}:{path}").stdout


def authorize_execution(root: Path, config_raw: bytes, result_raw: bytes,
                        authority_commit: str) -> TDG8RetryRecoveryExecutionAuthority:
    repository = _repository(root)
    config = _parse_config(config_raw)
    result = validate_compact(config_raw, result_raw)
    _require_git_image(repository, authority_commit)
    for path, expected in ((CONFIG_PATH, config_raw), (RESULT_PATH, result_raw)):
        live = _read_leaf(repository, path, f"RCV1 committed {path}", _MAX_RESULT_BYTES)
        if live != expected or _committed_bytes(repository, authority_commit, path) != expected:
            raise TDG8RetryRecoveryAuthorityError(f"RCV1 committed bytes differ: {path}")
    original = _read_leaf(repository, ORIGINAL_AUTHORITY_RESULT_PATH, "original authority result", _MAX_RESULT_BYTES)
    if (
        _sha(original) != ORIGINAL_AUTHORITY_RESULT_SHA256
        or _committed_bytes(repository, ORIGINAL_EXECUTION_COMMIT, ORIGINAL_AUTHORITY_RESULT_PATH) != original
        or _committed_bytes(repository, authority_commit, ORIGINAL_AUTHORITY_RESULT_PATH) != original
    ):
        raise TDG8RetryRecoveryAuthorityError("original one-time authority result was replaced")
    rows: list[tuple[str, str, str]] = []
    for item in result["artifact_payload"]["implementation_inventory"]:
        role, path, digest = str(item["role"]), str(item["path"]), str(item["sha256"])
        live = _read_leaf(repository, path, f"RCV1 implementation {path}", _MAX_BOUND_BYTES)
        if _sha(live) != digest or _committed_bytes(repository, authority_commit, path) != live:
            raise TDG8RetryRecoveryAuthorityError(f"RCV1 implementation bytes differ: {path}")
        rows.append((role, path, digest))
    environment = result["artifact_payload"].get("environment")
    if environment != observed_environment():
        raise TDG8RetryRecoveryAuthorityError("RCV1 execution environment differs")
    return TDG8RetryRecoveryExecutionAuthority(
        authority_commit=authority_commit,
        original_execution_commit=ORIGINAL_EXECUTION_COMMIT,
        original_plan_sha256=ORIGINAL_PLAN_SHA256,
        config_sha256=_sha(config_raw), result_sha256=_sha(result_raw),
        campaign_id=CAMPAIGN_ID, destination_path=DESTINATION_PATH,
        anchor_checkpoint_sha256=ANCHOR_CHECKPOINT_SHA256,
        anchor_journal_sha256=ANCHOR_JOURNAL_SHA256,
        expected_checkpoint_sha256=EXPECTED_CHECKPOINT_SHA256,
        expected_journal_sha256=EXPECTED_JOURNAL_SHA256,
        implementation_inventory=tuple(rows),
        environment=MappingProxyType(dict(environment)),
    )


__all__ = [
    "ANCHOR_CHECKPOINT_SHA256", "ANCHOR_JOURNAL_SHA256", "ARTIFACT_ID",
    "AUTHORITY_DELTA_PATHS", "CONFIG_PATH", "DESTINATION_PATH",
    "EXPECTED_CHECKPOINT_SHA256", "EXPECTED_JOURNAL_SHA256",
    "IMPLEMENTATION_INVENTORY", "ORIGINAL_EXECUTION_COMMIT",
    "ORIGINAL_PLAN_SHA256", "RESULT_PATH", "TDG8RetryRecoveryAuthorityError",
    "TDG8RetryRecoveryExecutionAuthority", "authorize_execution",
    "build_prelaunch", "canonical", "inspect_live_anchor",
    "inspect_recovery_boundary", "validate_compact",
]
