"""Read-only SID1 authority for the persisted PRO19 generation-seven suffix.

SID1 freezes only the identity and recursive decoding of one existing TDG6
rejection.  It neither acquires a writer lease nor reconciles, advances, or
changes a campaign store.  Keeping this reader independent of the campaign
store's recovery coordinator lets it authenticate the historical suffix while
that coordinator's nested-object decoder is the object under repair.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import os
from pathlib import Path
import stat
import tomllib
from typing import Any, Mapping

import numpy as np

from . import tdg5_stage_complete_refinement_runtime as tdg5
from . import tdg6_temporal_admission_runtime as tdg6
from .hlt16_state_store import canonical, decode_payload, read_nofollow
from .numerical_engine import array_content_sha256
from .proto19_launch_authority import (
    LaunchAuthorityReceipt,
    Proto19LaunchAuthorityError,
    _nofollow_regular_bytes,
    authorize_first_event,
)
from .proto19_progression_contract import ProgressionPlan, construct_first_event


ARTIFACT_ID = "FGC-1-PRO19-SID1-FRZ1"
SCHEMA_VERSION = 1
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CLASSIFICATION = "premise_only_exact_generation7_suffix_identity_and_decoder_contract_frozen"
CONFIG_PATH = "configs/fgc/fgc-1-pro19-sid1-frz1.toml"
RESULT_PATH = "results/fgc-1-pro19-sid1-frz1.json"
STORE_PATH = "runs/fgc-2-sf1/proto17/calibration"

HISTORICAL_COMMIT = "5d913f84260a4ed0e3058d53ef7dcb18cab79d58"
ORIGINAL_AUTHORIZATION_COMMIT = "c11ba422ce49ddd6f6175de1a9d2da675b97f2de"
PLAN_SHA256 = "f3b365978aeb2ae6063742807a3161bf5d6aa0c8392cd283fdec2a79d3e46060"
CAMPAIGN_ID = "FGC-2-SF1-PROTO17-GR0-A3-CAL-290a65bcd6a2ba820683"

PARENT_CHECKPOINT_SHA256 = "6277b3be4ffe1bd79b858e2b4770c244b88b9971cbff67ffb1a40d40f14e2e63"
CHECKPOINT_SHA256 = "c0ebfe09ef4058e117e46d34424edc8139419bfb4680cec979a1f5ab64d37be9"
CFL_RECORD_SHA256 = "78e65c08380a4fc979269e21db314232fd3fdc0e428c7e259f2031e66bd30bf7"
TDG6_SUFFIX_SHA256 = "e23d81705f1038f8240249658d0160cd5ca8ea8e33969cb4b32de3f03d426473"
DESCRIPTOR_SHA256 = "6b0dd985ab95ac52950348df5079a6b8fc4f5e7a78cf95031237cb0155c49c56"
SEMANTIC_PAYLOAD_SHA256 = "0e48e0e2396ea7d9ef79d1641e6944a2dfdaa3716bb65e0b882b97156ff9372f"
RAW_ARCHIVE_SHA256 = "f47a4f533eb8e0ed9f8ea337d520fa6f4e0bb6c1b27a4ef8023ee78bc830f5a4"
EVOLUTION_STATE_SHA256 = "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a"

MEMBER_KEY = "RK4-2049"
MEMBER_TIME_HEX = "0x1.78554de5a30e0p+0"
TDG6_LEDGER_SHA256 = "2763399f9824fcb476dc2e617d35cf1510ef08f11e65757ed5ef7e7840fd4984"
UPQ_BYTES_SHA256 = {
    "u": "856fc327d3163d038f9fdbfc660968509dd7360faa0bb9d158721eb3049c2fe4",
    "p": "32939115b86938cf8a813bbf9a50a8c90bbaceaee161fa65695deb435d1243b7",
    "q": "48ec07bcace3935173c06e7b41915cacb17f14e1c9d3500522d4eb59a205c4a3",
}
CHECKPOINT_INVENTORY_COUNT = 8
CHECKPOINT_INVENTORY_SHA256 = (
    "301b1ec030b820565bbd7856811160d613a1cc1301ea141a12d9bbf13088966c"
)
JOURNAL_INVENTORY_COUNT = 8
JOURNAL_INVENTORY_SHA256 = (
    "7c60adf94e78b76bc62856fd59c1b8b8cfbd099a696f2a74d4b226f8f94604da"
)
RECOVERY_TRANSITION_SHA256 = (
    "949b29e46081bf5a42c34461b3894419db86aeca6c733b7a2636c91c6cb2ee4a"
)
RECOVERY_CHECKPOINT_SHA256 = (
    "6d2b6e37e7c7cfd4ce64504353d886099003f959ea2aaee1a14fc778ce14fa46"
)
RECOVERY_CURSOR_SHA256 = (
    "6d594999f6aaebb5d56f0fd6ff12499e2b22a47e3daf74ac8b13b15d273fa177"
)
RECOVERY_PENDING_CAP_HEX = "0x1.aaa9612df9000p-10"

EXPECTED_RECOVERY_CONTRACT: dict[str, object] = {
    "generation": 8,
    "journal_sequence": 8,
    "transition_record_sha256": RECOVERY_TRANSITION_SHA256,
    "checkpoint_sha256": RECOVERY_CHECKPOINT_SHA256,
    "member_cursor_sha256": RECOVERY_CURSOR_SHA256,
    "member_descriptor_sha256": DESCRIPTOR_SHA256,
    "evolution_state_sha256": EVOLUTION_STATE_SHA256,
    "accepted_time_hex": MEMBER_TIME_HEX,
    "cfl_current": 1,
    "cfl_total": 1,
    "temporal_retry_total": 1,
    "pending_owner": "temporal",
    "pending_cap_hex": RECOVERY_PENDING_CAP_HEX,
    "physical_state_advanced": False,
    "rejected_proposal_replayed": False,
    "candidate_branch_opened": False,
}


class SID1AuthorityError(RuntimeError):
    """A SID1 static identity, schema, or recursive decoder check failed."""


@dataclass(frozen=True, slots=True)
class SID1RecoveryReceipt:
    """One committed corrected image bound to the exact pre-recovery suffix."""

    continuation_commit: str
    manifest_sha256: str
    config_sha256: str
    result_sha256: str
    original_authorization_commit: str
    original_plan_sha256: str
    campaign_id: str
    recovery_checkpoint_generation: int
    recovery_checkpoint_sha256: str
    checkpoint_journal_sequence: int
    checkpoint_journal_tip_sha256: str
    suffix_sequence: int
    suffix_sha256: str
    member_key: str
    descriptor_sha256: str
    evolution_state_sha256: str
    expected_transition_sha256: str
    expected_checkpoint_sha256: str
    expected_cursor_sha256: str
    expected_pending_cap_hex: str
    authority_paths: tuple[tuple[str, str, str], ...]
    environment: Mapping[str, str]
    progression_plan: ProgressionPlan


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _canonical_result(value: object) -> bytes:
    try:
        return (
            json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False)
            + "\n"
        ).encode("ascii")
    except (TypeError, ValueError) as exc:
        raise SID1AuthorityError("result is not canonical JSON") from exc


def _json(raw: bytes, label: str) -> dict[str, Any]:
    def pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
        answer: dict[str, Any] = {}
        for key, value in items:
            if key in answer:
                raise ValueError(key)
            answer[key] = value
        return answer
    try:
        value = json.loads(raw.decode("ascii"), object_pairs_hook=pairs)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise SID1AuthorityError(f"{label} is malformed") from exc
    if not isinstance(value, dict) or canonical(value) != raw:
        raise SID1AuthorityError(f"{label} is noncanonical")
    return value


def _require_hash(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or value.lower() != value:
        raise SID1AuthorityError(f"{label} is not SHA-256")
    try:
        int(value, 16)
    except ValueError as exc:
        raise SID1AuthorityError(f"{label} is not SHA-256") from exc
    return value


def _object_hash(value: Mapping[str, Any], key: str, label: str) -> None:
    supplied = _require_hash(value.get(key), label)
    body = dict(value)
    body.pop(key, None)
    if _sha(canonical(body)) != supplied:
        raise SID1AuthorityError(f"{label} differs")


def _read(root: Path, relative: str, label: str) -> bytes:
    try:
        return read_nofollow(root, relative, label)
    except Exception as exc:
        raise SID1AuthorityError(f"cannot safely read {label}") from exc


def _inventory_sha(names: tuple[str, ...]) -> str:
    raw = (
        json.dumps(
            list(names),
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("ascii")
    return _sha(raw)


def _safe_inventory(root: Path, relative: str, label: str) -> tuple[str, ...]:
    """List one directory through no-follow descriptors without mutation."""
    root_fd = directory_fd = -1
    try:
        root_fd = os.open(
            root,
            os.O_RDONLY
            | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_NOFOLLOW", 0),
        )
        if not stat.S_ISDIR(os.fstat(root_fd).st_mode):
            raise SID1AuthorityError(f"{label} root is unsafe")
        directory_fd = os.open(
            relative,
            os.O_RDONLY
            | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=root_fd,
        )
        if not stat.S_ISDIR(os.fstat(directory_fd).st_mode):
            raise SID1AuthorityError(f"{label} directory is unsafe")
        names: list[str] = []
        with os.scandir(directory_fd) as entries:
            for entry in entries:
                metadata = entry.stat(follow_symlinks=False)
                if (
                    stat.S_ISLNK(metadata.st_mode)
                    or not stat.S_ISREG(metadata.st_mode)
                    or metadata.st_nlink != 1
                ):
                    raise SID1AuthorityError(f"{label} leaf is unsafe")
                names.append(entry.name)
        return tuple(sorted(names))
    except SID1AuthorityError:
        raise
    except OSError as exc:
        raise SID1AuthorityError(f"cannot safely inventory {label}") from exc
    finally:
        if directory_fd != -1:
            os.close(directory_fd)
        if root_fd != -1:
            os.close(root_fd)


def _interval(value: object, label: str) -> tdg6.TDG6Binary64MagnitudeInterval:
    if not isinstance(value, Mapping):
        raise SID1AuthorityError(f"{label} is not a mapping")
    raw = dict(value)
    envelope = raw.get("envelope")
    if not isinstance(envelope, Mapping):
        raise SID1AuthorityError(f"{label} omits its cubic envelope")
    try:
        raw["envelope"] = tdg5.Binary64CubicEnvelope(**dict(envelope))
        return tdg6.TDG6Binary64MagnitudeInterval(**raw)
    except (TypeError, ValueError) as exc:
        raise SID1AuthorityError(f"{label} differs from its frozen nested schema") from exc


def decode_tdg6_evidence(value: Mapping[str, Any]) -> tdg6.TDG6TemporalRetryEvidence:
    """Recursively decode the stored evidence without executing a PDE proposal."""
    raw = dict(value)
    admissions = raw.get("channel_admissions")
    if not isinstance(admissions, list):
        raise SID1AuthorityError("TDG6 evidence omits ordered channel admissions")
    decoded = []
    for index, item in enumerate(admissions):
        if not isinstance(item, Mapping):
            raise SID1AuthorityError("TDG6 channel admission is not a mapping")
        channel = dict(item)
        channel["outer_difference"] = _interval(channel.get("outer_difference"), f"outer interval {index}")
        channel["finest_difference"] = _interval(channel.get("finest_difference"), f"finest interval {index}")
        try:
            decoded.append(tdg6.TDG6RuntimeChannelAdmission(**channel))
        except (TypeError, ValueError) as exc:
            raise SID1AuthorityError(f"TDG6 channel admission {index} differs") from exc
    raw["channel_admissions"] = tuple(decoded)
    if isinstance(raw.get("failed_channels"), list):
        raw["failed_channels"] = tuple(raw["failed_channels"])
    try:
        return tdg6.TDG6TemporalRetryEvidence(**raw)
    except (TypeError, ValueError) as exc:
        raise SID1AuthorityError("TDG6 retry evidence differs") from exc


def expected_anchor() -> dict[str, object]:
    return {
        "historical_commit": HISTORICAL_COMMIT,
        "authorization_commit": ORIGINAL_AUTHORIZATION_COMMIT,
        "plan_sha256": PLAN_SHA256,
        "campaign_id": CAMPAIGN_ID,
        "protocol": TARGET_PROTOCOL,
        "branch": "GR-0",
        "amplitude": "3",
        "event": 23,
        "target_rational": "3/2",
        "parent_checkpoint_sha256": PARENT_CHECKPOINT_SHA256,
        "checkpoint_generation": 7,
        "checkpoint_sha256": CHECKPOINT_SHA256,
        "checkpoint_journal_sequence": 6,
        "checkpoint_journal_tip_sha256": CFL_RECORD_SHA256,
        "cfl_record_sequence": 6,
        "cfl_record_sha256": CFL_RECORD_SHA256,
        "tdg6_suffix_sequence": 7,
        "tdg6_suffix_sha256": TDG6_SUFFIX_SHA256,
        "member_key": MEMBER_KEY,
        "descriptor_sha256": DESCRIPTOR_SHA256,
        "semantic_payload_sha256": SEMANTIC_PAYLOAD_SHA256,
        "raw_archive_sha256": RAW_ARCHIVE_SHA256,
        "evolution_state_sha256": EVOLUTION_STATE_SHA256,
        "member_accepted_time_hex": MEMBER_TIME_HEX,
        "tdg6_ledger_sha256": TDG6_LEDGER_SHA256,
        "upq_bytes_sha256": dict(UPQ_BYTES_SHA256),
        "checkpoint_inventory_count": CHECKPOINT_INVENTORY_COUNT,
        "checkpoint_inventory_sha256": CHECKPOINT_INVENTORY_SHA256,
        "journal_inventory_count": JOURNAL_INVENTORY_COUNT,
        "journal_inventory_sha256": JOURNAL_INVENTORY_SHA256,
        "checkpoint_disposition": "nonterminal",
        "checkpoint_suffix_kinds": ["tdg6_rejection"],
        "state_advanced": False,
        "candidate_branch_opened": False,
        "physical_result_earned": False,
    }


def derive_live_anchor(repository: Path) -> dict[str, object]:
    """Authenticate exactly the SID1 raw checkpoint/suffix/payload facts.

    This intentionally reads only the explicit leaves required by SID1.  It
    does not use the campaign store's recovery decoder and has no writer or
    publication capability.
    """
    root = Path(repository) / STORE_PATH
    checkpoint_inventory = _safe_inventory(root, "checkpoints", "checkpoint")
    journal_inventory = _safe_inventory(root, "journal", "journal")
    if (
        len(checkpoint_inventory) != CHECKPOINT_INVENTORY_COUNT
        or _inventory_sha(checkpoint_inventory) != CHECKPOINT_INVENTORY_SHA256
        or len(journal_inventory) != JOURNAL_INVENTORY_COUNT
        or _inventory_sha(journal_inventory) != JOURNAL_INVENTORY_SHA256
    ):
        raise SID1AuthorityError("SID1 store moved beyond its frozen inventory")
    checkpoint_rel = f"checkpoints/{7:020d}-{CHECKPOINT_SHA256}.json"
    cfl_rel = f"journal/{6:020d}-{CFL_RECORD_SHA256}.journal"
    suffix_rel = f"journal/{7:020d}-{TDG6_SUFFIX_SHA256}.journal"
    descriptor_rel = f"states/{DESCRIPTOR_SHA256}.json"
    checkpoint = _json(_read(root, checkpoint_rel, "generation-seven checkpoint"), "generation-seven checkpoint")
    cfl = _json(_read(root, cfl_rel, "CFL rejection record"), "CFL rejection record")
    suffix = _json(_read(root, suffix_rel, "TDG6 suffix record"), "TDG6 suffix record")
    descriptor = _json(_read(root, descriptor_rel, "persisted member descriptor"), "persisted member descriptor")
    _object_hash(checkpoint, "checkpoint_sha256", "generation-seven checkpoint hash")
    _object_hash(cfl, "record_sha256", "CFL rejection record hash")
    _object_hash(suffix, "record_sha256", "TDG6 suffix record hash")
    _object_hash(descriptor, "descriptor_sha256", "member descriptor hash")
    if checkpoint.get("checkpoint_sha256") != CHECKPOINT_SHA256 or checkpoint.get("parent_sha256") != PARENT_CHECKPOINT_SHA256:
        raise SID1AuthorityError("generation-seven checkpoint identity differs")
    if checkpoint.get("generation") != 7 or checkpoint.get("journal_sequence") != 6 or checkpoint.get("journal_tip_sha256") != CFL_RECORD_SHA256:
        raise SID1AuthorityError("generation-seven checkpoint journal identity differs")
    if cfl.get("record_sha256") != CFL_RECORD_SHA256 or suffix.get("record_sha256") != TDG6_SUFFIX_SHA256 or descriptor.get("descriptor_sha256") != DESCRIPTOR_SHA256:
        raise SID1AuthorityError("SID1 leaf address identity differs")
    if checkpoint.get("authorization_commit") != ORIGINAL_AUTHORIZATION_COMMIT or checkpoint.get("plan_sha256") != PLAN_SHA256 or checkpoint.get("campaign_id") != CAMPAIGN_ID or checkpoint.get("protocol") != TARGET_PROTOCOL:
        raise SID1AuthorityError("generation-seven campaign identity differs")
    if checkpoint.get("target", {}).get("rational") != "3/2" or checkpoint.get("disposition") != "nonterminal":
        raise SID1AuthorityError("generation-seven target/disposition differs")
    if cfl.get("sequence") != 6 or cfl.get("kind") != "cfl_rejection" or cfl.get("previous_record_sha256") != "bf92cbbab28fc55259f1d6cd1249fa73d5e740f067e4e4204e0dd3d1cfb854e5":
        raise SID1AuthorityError("CFL record chain differs")
    if suffix.get("sequence") != 7 or suffix.get("kind") != "tdg6_rejection" or suffix.get("previous_record_sha256") != CFL_RECORD_SHA256:
        raise SID1AuthorityError("TDG6 suffix chain differs")
    if (
        cfl.get("generation") != 7
        or suffix.get("generation") != 8
        or any(item.get("campaign_id") != CAMPAIGN_ID or item.get("event") != 23 for item in (cfl, suffix))
    ):
        raise SID1AuthorityError("journal campaign identity differs")
    member = checkpoint.get("members", {}).get(MEMBER_KEY)
    if not isinstance(member, Mapping) or member.get("descriptor_sha256") != DESCRIPTOR_SHA256:
        raise SID1AuthorityError("checkpoint member descriptor differs")
    cursor = member.get("cursor", {})
    if cursor.get("accepted_state_sha256") != DESCRIPTOR_SHA256 or cursor.get("accepted_boundary_time", {}).get("binary64_hex") != MEMBER_TIME_HEX or cursor.get("mode") != "FRESH_READY":
        raise SID1AuthorityError("checkpoint member cursor differs")
    if member.get("cfl_current") != 1 or member.get("cfl_total") != 1 or member.get("pending_owner") != "cfl" or member.get("pending_cap_hex") != "0x1.aaa9612df9800p-9":
        raise SID1AuthorityError("checkpoint CFL retry state differs")
    if cursor.get("TDG6_ledger_sha256") != TDG6_LEDGER_SHA256:
        raise SID1AuthorityError("checkpoint TDG6 ledger identity differs")
    payload = suffix.get("payload")
    if not isinstance(payload, Mapping) or payload.get("member_key") != MEMBER_KEY or payload.get("predecessor_descriptor_sha256") != DESCRIPTOR_SHA256:
        raise SID1AuthorityError("TDG6 suffix member identity differs")
    evidence = payload.get("evidence")
    if not isinstance(evidence, Mapping):
        raise SID1AuthorityError("TDG6 suffix omits evidence")
    decoded = decode_tdg6_evidence(evidence)
    if decoded.initial_state_sha256 != EVOLUTION_STATE_SHA256 or decoded.failed_channels != ("u:R",) or decoded.retry_count_for_current_macro_step != 1:
        raise SID1AuthorityError("decoded TDG6 evidence identity differs")
    if payload.get("evidence", {}).get("accepted_state_bitwise_preserved") is not True or payload.get("evidence", {}).get("accumulated_debit_unchanged") is not True:
        raise SID1AuthorityError("TDG6 suffix crossed its no-commit boundary")
    if descriptor.get("semantic_sha256") != SEMANTIC_PAYLOAD_SHA256 or descriptor.get("raw_archive_sha256") != RAW_ARCHIVE_SHA256:
        raise SID1AuthorityError("member payload identity differs")
    arrays_manifest = descriptor.get("arrays")
    if not isinstance(arrays_manifest, list):
        raise SID1AuthorityError("member array manifest differs")
    manifest = {item.get("name"): item for item in arrays_manifest if isinstance(item, Mapping)}
    if set(UPQ_BYTES_SHA256) - set(manifest):
        raise SID1AuthorityError("member u/p/q manifest differs")
    if any(manifest[name].get("bytes_sha256") != digest for name, digest in UPQ_BYTES_SHA256.items()):
        raise SID1AuthorityError("member u/p/q byte identity differs")
    payload_rel = f"payloads/{SEMANTIC_PAYLOAD_SHA256}.npz"
    raw_payload = _read(root, payload_rel, "persisted member payload")
    if _sha(raw_payload) != RAW_ARCHIVE_SHA256:
        raise SID1AuthorityError("member raw archive hash differs")
    try:
        arrays, semantic, observed_manifest = decode_payload(raw_payload, expected_semantic_sha256=SEMANTIC_PAYLOAD_SHA256, expected_manifest=arrays_manifest)
    except Exception as exc:
        raise SID1AuthorityError("member payload cannot be safely decoded") from exc
    if semantic != SEMANTIC_PAYLOAD_SHA256 or list(observed_manifest) != arrays_manifest:
        raise SID1AuthorityError("member payload semantic manifest differs")
    if array_content_sha256(arrays["u"], arrays["p"], arrays["q"]) != EVOLUTION_STATE_SHA256:
        raise SID1AuthorityError("member u/p/q evolution identity differs")
    observed_upq = {name: _sha(np.asarray(arrays[name], dtype="<f8", order="C").tobytes(order="C")) for name in UPQ_BYTES_SHA256}
    if observed_upq != UPQ_BYTES_SHA256:
        raise SID1AuthorityError("member u/p/q byte hashes differ")
    return expected_anchor()


def _validate_config(value: Mapping[str, Any]) -> dict[str, Any]:
    required = {"schema_version", "artifact_id", "project_version", "target_protocol", "classification", "live_anchor", "decoder_contract", "expected_recovery", "scope", "claims", "nonclaims"}
    if set(value) != required:
        raise SID1AuthorityError("SID1 config fields differ")
    if value["schema_version"] != SCHEMA_VERSION or value["artifact_id"] != ARTIFACT_ID or value["target_protocol"] != TARGET_PROTOCOL or value["classification"] != CLASSIFICATION or value["live_anchor"] != expected_anchor():
        raise SID1AuthorityError("SID1 config identity differs")
    if value["decoder_contract"] != {"recursive_envelope_decode_required": True, "in_memory_only": True, "PDE_proposal_execution": False, "recovery_publication": False}:
        raise SID1AuthorityError("SID1 decoder contract differs")
    if value["expected_recovery"] != EXPECTED_RECOVERY_CONTRACT:
        raise SID1AuthorityError("SID1 expected recovery differs")
    if value["scope"] != {"raw_store_mutation": False, "writer_lease_acquisition": False, "checkpoint_reconciliation": False, "candidate_branch_opened": False, "GEN0_reimport": False}:
        raise SID1AuthorityError("SID1 scope differs")
    if value["claims"] != {"SID1_contract_frozen": True, "corrected_decoder_image_authenticated": False, "runtime_recovery_receipt_required": True, "GR0_calibration_completed": False, "candidate_execution_authorized": False, "physical_result_earned": False}:
        raise SID1AuthorityError("SID1 claims differ")
    if not isinstance(value["project_version"], str) or not value["project_version"] or not isinstance(value["nonclaims"], list):
        raise SID1AuthorityError("SID1 presentation fields differ")
    return dict(value)


def build_sid1_result(config_raw: bytes, anchor: Mapping[str, Any]) -> dict[str, Any]:
    try:
        config = tomllib.loads(config_raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise SID1AuthorityError("SID1 config is malformed") from exc
    if not isinstance(config, Mapping):
        raise SID1AuthorityError("SID1 config is not a table")
    checked = _validate_config(config)
    if dict(anchor) != expected_anchor():
        raise SID1AuthorityError("SID1 live anchor differs")
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "project_version": checked["project_version"],
        "target_protocol": TARGET_PROTOCOL,
        "classification": CLASSIFICATION,
        "gate_status": "pass",
        "source_config_sha256": _sha(config_raw),
        "artifact_payload": {
            "live_anchor": expected_anchor(),
            "decoder_contract": dict(checked["decoder_contract"]),
            "expected_recovery": dict(checked["expected_recovery"]),
            "scope": dict(checked["scope"]),
            "claims": dict(checked["claims"]),
            "nonclaims": list(checked["nonclaims"]),
        },
    }


def _parse_config(raw: bytes) -> dict[str, Any]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise SID1AuthorityError("SID1 config is malformed") from exc
    if not isinstance(value, Mapping):
        raise SID1AuthorityError("SID1 config is not a table")
    return _validate_config(value)


def _parse_result(raw: bytes) -> dict[str, Any]:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, item in pairs:
            if key in value:
                raise ValueError(key)
            value[key] = item
        return value

    try:
        result = json.loads(
            raw.decode("ascii"), object_pairs_hook=reject_duplicates
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise SID1AuthorityError("SID1 result is malformed") from exc
    if not isinstance(result, dict) or _canonical_result(result) != raw:
        raise SID1AuthorityError("SID1 result is noncanonical")
    return result


def _captured_bytes(
    root: Path, receipt: LaunchAuthorityReceipt, relative: str
) -> bytes:
    captured = {path: digest for _role, path, digest in receipt.authority_paths}
    if relative not in captured:
        raise SID1AuthorityError("SID1 predecessor is outside the launch manifest")
    try:
        raw = _nofollow_regular_bytes(root, relative)
    except Proto19LaunchAuthorityError as exc:
        raise SID1AuthorityError("SID1 predecessor cannot be read safely") from exc
    if _sha(raw) != captured[relative]:
        raise SID1AuthorityError("SID1 predecessor differs from committed image")
    return raw


def authorize_sid1_recovery(
    root: Path,
    *,
    continuation_commit: str,
    launch_manifest_path: str,
    store_anchor: Mapping[str, Any],
) -> SID1RecoveryReceipt:
    """Authorize only the exact decoded sequence-seven recovery edge."""
    repository = Path(root)
    try:
        current = authorize_first_event(
            repository,
            authorization_commit=continuation_commit,
            manifest_path=launch_manifest_path,
        )
    except Proto19LaunchAuthorityError as exc:
        raise SID1AuthorityError("corrected SID1 execution image differs") from exc
    if current.progression_plan.branch != "GR-0":
        raise SID1AuthorityError("corrected SID1 image opened a candidate branch")

    config_raw = _captured_bytes(repository, current, CONFIG_PATH)
    result_raw = _captured_bytes(repository, current, RESULT_PATH)
    _parse_config(config_raw)
    result = _parse_result(result_raw)
    if result != build_sid1_result(config_raw, store_anchor):
        raise SID1AuthorityError("SID1 result/config/store cross-binding differs")

    progression_paths = (
        "results/fgc-1-pro19-frz1.json",
        "results/fgc-1-pro18-auth1.json",
    )
    captured = {
        path: _captured_bytes(repository, current, path)
        for path in progression_paths
    }
    original_plan = construct_first_event(
        repository,
        authorization_commit=ORIGINAL_AUTHORIZATION_COMMIT,
        evidence_bytes=captured,
    )
    if (
        original_plan.sha256 != PLAN_SHA256
        or original_plan.campaign_id != CAMPAIGN_ID
        or original_plan.branch != "GR-0"
        or original_plan.amplitude != "3"
    ):
        raise SID1AuthorityError("SID1 original progression plan differs")

    return SID1RecoveryReceipt(
        continuation_commit=continuation_commit,
        manifest_sha256=current.manifest_sha256,
        config_sha256=_sha(config_raw),
        result_sha256=_sha(result_raw),
        original_authorization_commit=ORIGINAL_AUTHORIZATION_COMMIT,
        original_plan_sha256=PLAN_SHA256,
        campaign_id=CAMPAIGN_ID,
        recovery_checkpoint_generation=7,
        recovery_checkpoint_sha256=CHECKPOINT_SHA256,
        checkpoint_journal_sequence=6,
        checkpoint_journal_tip_sha256=CFL_RECORD_SHA256,
        suffix_sequence=7,
        suffix_sha256=TDG6_SUFFIX_SHA256,
        member_key=MEMBER_KEY,
        descriptor_sha256=DESCRIPTOR_SHA256,
        evolution_state_sha256=EVOLUTION_STATE_SHA256,
        expected_transition_sha256=RECOVERY_TRANSITION_SHA256,
        expected_checkpoint_sha256=RECOVERY_CHECKPOINT_SHA256,
        expected_cursor_sha256=RECOVERY_CURSOR_SHA256,
        expected_pending_cap_hex=RECOVERY_PENDING_CAP_HEX,
        authority_paths=current.authority_paths,
        environment=current.environment,
        progression_plan=original_plan,
    )


__all__ = [
    "ARTIFACT_ID", "CHECKPOINT_SHA256", "CFL_RECORD_SHA256", "CONFIG_PATH",
    "DESCRIPTOR_SHA256", "EXPECTED_RECOVERY_CONTRACT", "RESULT_PATH",
    "SID1AuthorityError",
    "SID1RecoveryReceipt", "STORE_PATH", "TDG6_SUFFIX_SHA256",
    "authorize_sid1_recovery", "build_sid1_result", "decode_tdg6_evidence",
    "derive_live_anchor", "expected_anchor",
]
