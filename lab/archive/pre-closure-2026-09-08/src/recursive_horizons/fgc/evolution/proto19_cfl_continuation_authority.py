"""One-use authority for the corrected PRO19 CFL-retry continuation.

The first production attempt stopped before publishing its proposed CFL
rejection because the coordinator mixed an outer requested cap with the
smaller TDG7 macro width that was actually attempted.  This module does not
rewrite that campaign.  It authenticates one later committed execution image
and permits it to resume only the exact original plan from the exact last
healthy generation-six checkpoint.

This module is deliberately store-agnostic and read-only: the runner supplies
an independently derived store anchor, and the existing campaign store remains
the sole owner of stale-writer takeover and durable mutation.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import tomllib
from typing import Any, Mapping

from .proto19_launch_authority import (
    LaunchAuthorityReceipt,
    Proto19LaunchAuthorityError,
    _nofollow_regular_bytes,
    authorize_first_event,
)
from .proto19_progression_contract import ProgressionPlan, construct_first_event


ARTIFACT_ID = "FGC-1-PRO19-CFL1-FRZ1"
SCHEMA_VERSION = 1
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CLASSIFICATION = "premise_only_exact_checkpoint_continuation_frozen"
CONFIG_PATH = "configs/fgc/fgc-1-pro19-cfl1-frz1.toml"
RESULT_PATH = "results/fgc-1-pro19-cfl1-frz1.json"

ORIGINAL_COMMIT = "c11ba422ce49ddd6f6175de1a9d2da675b97f2de"
ORIGINAL_PLAN = "f3b365978aeb2ae6063742807a3161bf5d6aa0c8392cd283fdec2a79d3e46060"
ORIGINAL_CAMPAIGN = "FGC-2-SF1-PROTO17-GR0-A3-CAL-290a65bcd6a2ba820683"
GENERATION = 6
CHECKPOINT_SHA256 = "6277b3be4ffe1bd79b858e2b4770c244b88b9971cbff67ffb1a40d40f14e2e63"
JOURNAL_SEQUENCE = 5
JOURNAL_TIP_SHA256 = "bf92cbbab28fc55259f1d6cd1249fa73d5e740f067e4e4204e0dd3d1cfb854e5"
MEMBER_KEY = "RK4-2049"
MEMBER_DESCRIPTOR_SHA256 = "6b0dd985ab95ac52950348df5079a6b8fc4f5e7a78cf95031237cb0155c49c56"
MEMBER_TIME_HEX = "0x1.78554de5a30e0p+0"

ORIGINAL_CAMPAIGN_CONTRACT: dict[str, object] = {
    "authorization_commit": ORIGINAL_COMMIT,
    "plan_sha256": ORIGINAL_PLAN,
    "campaign_id": ORIGINAL_CAMPAIGN,
    "branch": "GR-0",
    "amplitude": "3",
    "event": 23,
    "target_rational": "3/2",
}

SAFE_RECOVERY_CONTRACT: dict[str, object] = {
    "checkpoint_generation": GENERATION,
    "checkpoint_sha256": CHECKPOINT_SHA256,
    "journal_sequence": JOURNAL_SEQUENCE,
    "journal_tip_sha256": JOURNAL_TIP_SHA256,
    "disposition": "nonterminal",
    "suffix_classification": "clean",
    "suffix_kinds": [],
    "staging_paths": [],
    "terminal_lock_present": False,
    "writer_state": "stale_verified_lock",
    "active_write": True,
    "member_key": MEMBER_KEY,
    "member_descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
    "member_accepted_time_hex": MEMBER_TIME_HEX,
    "member_mode": "FRESH_READY",
    "source_retry_total": 0,
    "cfl_retry_total": 0,
    "temporal_retry_total": 0,
    "pending_owner": "none",
    "pending_cap_hex": "none",
    "durable_cfl_rejection_exists": False,
    "durable_terminal_exists": False,
}

DIAGNOSED_STOP_CONTRACT: dict[str, object] = {
    "classification": "invalid_implementation_or_nonconverged_run",
    "detail": "cfl_rejection cap reduction differs",
    "owner_layer": "HLT16_source_CFL_rejection_serialization",
    "accepted_fine_steps_before_stop": 5,
    "outer_requested_cap_hex": "0x1.aaa9612df9ba5p-8",
    "actual_macro_width_hex": "0x1.aaa9612df9800p-8",
    "correct_successor_cap_hex": "0x1.aaa9612df9800p-9",
    "observed_CFL_ratio_hex": "0x1.000004e8e3345p-3",
    "maximum_CFL_ratio_hex": "0x1.0000000000000p-3",
    "accepted_boundary_restored": True,
    "durable_rejection_published": False,
    "physical_inference_earned": False,
}

SCOPE_CONTRACT: dict[str, object] = {
    "permitted_action": "repair_retry_record_identity_and_resume_exact_generation_6_edge",
    "target_or_threshold_change": False,
    "retry_cap_change": False,
    "physical_input_change": False,
    "plan_change": False,
    "checkpoint_rollback": False,
    "candidate_branch_opened": False,
    "raw_GEN0_reimport_permitted": False,
    "manual_stale_lock_removal_permitted": False,
}

CLAIMS_CONTRACT: dict[str, object] = {
    "continuation_contract_frozen": True,
    "corrected_execution_image_authenticated": False,
    "corrected_execution_image_authentication_required": True,
    "exact_generation_6_continuation_authorized_by_artifact_alone": False,
    "original_campaign_identity_preserved": True,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "physical_result_earned": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
}

NONCLAIMS = (
    "The earlier stop is not a GR-0, trapping, FGC-QR, or physical result.",
    "This authority does not adopt another plan, checkpoint, campaign, or branch.",
    "This authority does not authorize SGB-L, FGC-QR, DEF1, or retained-EFT evolution.",
)


class CFLContinuationAuthorityError(RuntimeError):
    """The one-use continuation image or exact recovery anchor differs."""


@dataclass(frozen=True, slots=True)
class CFLContinuationReceipt:
    continuation_commit: str
    continuation_manifest_sha256: str
    continuation_config_sha256: str
    continuation_result_sha256: str
    original_authorization_commit: str
    original_plan_sha256: str
    campaign_id: str
    recovery_checkpoint_generation: int
    recovery_checkpoint_sha256: str
    recovery_journal_sequence: int
    recovery_journal_tip_sha256: str
    recovery_member_key: str
    recovery_member_descriptor_sha256: str
    authority_paths: tuple[tuple[str, str, str], ...]
    environment: Mapping[str, str]
    progression_plan: ProgressionPlan


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _canonical(value: object) -> bytes:
    try:
        return (
            json.dumps(
                value,
                sort_keys=True,
                indent=2,
                ensure_ascii=True,
                allow_nan=False,
            )
            + "\n"
        ).encode("ascii")
    except (TypeError, ValueError) as exc:
        raise CFLContinuationAuthorityError(
            "continuation result is not canonical JSON"
        ) from exc


def _validate_config_mapping(value: Mapping[str, Any]) -> dict[str, Any]:
    required = {
        "schema_version",
        "artifact_id",
        "project_version",
        "target_protocol",
        "classification",
        "original_campaign",
        "safe_recovery_ancestor",
        "diagnosed_stop",
        "scope",
        "claims",
        "nonclaims",
    }
    if set(value) != required:
        raise CFLContinuationAuthorityError("continuation config fields differ")
    if (
        value["schema_version"] != SCHEMA_VERSION
        or value["artifact_id"] != ARTIFACT_ID
        or value["target_protocol"] != TARGET_PROTOCOL
        or value["classification"] != CLASSIFICATION
        or not isinstance(value["project_version"], str)
        or not value["project_version"]
        or value["original_campaign"] != ORIGINAL_CAMPAIGN_CONTRACT
        or value["safe_recovery_ancestor"] != SAFE_RECOVERY_CONTRACT
        or value["diagnosed_stop"] != DIAGNOSED_STOP_CONTRACT
        or value["scope"] != SCOPE_CONTRACT
        or value["claims"] != CLAIMS_CONTRACT
        or value["nonclaims"] != list(NONCLAIMS)
    ):
        raise CFLContinuationAuthorityError("continuation config contract differs")
    return dict(value)


def _parse_config(raw: bytes) -> dict[str, Any]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise CFLContinuationAuthorityError("continuation config is malformed") from exc
    if not isinstance(value, Mapping):
        raise CFLContinuationAuthorityError("continuation config is not a table")
    return _validate_config_mapping(value)


def _parse_result(raw: bytes) -> dict[str, Any]:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(key)
            result[key] = value
        return result

    try:
        value = json.loads(raw.decode("ascii"), object_pairs_hook=reject_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise CFLContinuationAuthorityError("continuation result is malformed") from exc
    if not isinstance(value, dict) or _canonical(value) != raw:
        raise CFLContinuationAuthorityError("continuation result is noncanonical")
    return value


def expected_store_anchor() -> dict[str, object]:
    """Return the exact read-only recovery facts the runner must derive."""
    return {**ORIGINAL_CAMPAIGN_CONTRACT, **SAFE_RECOVERY_CONTRACT}


def build_cfl1_freeze(
    config: Mapping[str, Any], store_anchor: Mapping[str, Any]
) -> dict[str, Any]:
    """Build compact premise-only evidence without opening or mutating a store."""
    checked = _validate_config_mapping(config)
    if dict(store_anchor) != expected_store_anchor():
        raise CFLContinuationAuthorityError("continuation store anchor differs")
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "project_version": checked["project_version"],
        "target_protocol": TARGET_PROTOCOL,
        "classification": CLASSIFICATION,
        "gate_status": "pass",
        "source_config_sha256": "",
        "artifact_payload": {
            "original_campaign": dict(ORIGINAL_CAMPAIGN_CONTRACT),
            "safe_recovery_ancestor": dict(SAFE_RECOVERY_CONTRACT),
            "diagnosed_stop": dict(DIAGNOSED_STOP_CONTRACT),
            "scope": dict(SCOPE_CONTRACT),
            "claims": dict(CLAIMS_CONTRACT),
            "nonclaims": list(NONCLAIMS),
        },
    }


def build_cfl1_result(
    config_raw: bytes, store_anchor: Mapping[str, Any]
) -> dict[str, Any]:
    config = _parse_config(config_raw)
    result = build_cfl1_freeze(config, store_anchor)
    result["source_config_sha256"] = _sha(config_raw)
    return result


def _captured_bytes(
    root: Path,
    receipt: LaunchAuthorityReceipt,
    path: str,
) -> bytes:
    hashes = {
        relative: digest for _role, relative, digest in receipt.authority_paths
    }
    if path not in hashes:
        raise CFLContinuationAuthorityError(
            "continuation predecessor is outside the closed launch manifest"
        )
    try:
        raw = _nofollow_regular_bytes(root, path)
    except Proto19LaunchAuthorityError as exc:
        raise CFLContinuationAuthorityError(
            "continuation predecessor cannot be read safely"
        ) from exc
    if _sha(raw) != hashes[path]:
        raise CFLContinuationAuthorityError(
            "continuation predecessor differs from the authenticated image"
        )
    return raw


def authorize_cfl_continuation(
    root: Path,
    *,
    continuation_commit: str,
    launch_manifest_path: str,
    store_anchor: Mapping[str, Any],
) -> CFLContinuationReceipt:
    """Authenticate one corrected image for one exact persisted campaign edge."""
    repository = Path(root)
    try:
        current = authorize_first_event(
            repository,
            authorization_commit=continuation_commit,
            manifest_path=launch_manifest_path,
        )
    except Proto19LaunchAuthorityError as exc:
        raise CFLContinuationAuthorityError(
            "corrected execution image authority differs"
        ) from exc
    if current.progression_plan.branch != "GR-0":
        raise CFLContinuationAuthorityError("corrected image opened a candidate branch")

    config_raw = _captured_bytes(repository, current, CONFIG_PATH)
    result_raw = _captured_bytes(repository, current, RESULT_PATH)
    config = _parse_config(config_raw)
    result = _parse_result(result_raw)
    expected = build_cfl1_result(config_raw, store_anchor)
    if result != expected:
        raise CFLContinuationAuthorityError(
            "continuation result/config/store cross-binding differs"
        )

    progression_paths = (
        "results/fgc-1-pro19-frz1.json",
        "results/fgc-1-pro18-auth1.json",
    )
    captured = {
        path: _captured_bytes(repository, current, path)
        for path in progression_paths
    }
    plan = construct_first_event(
        repository,
        authorization_commit=ORIGINAL_COMMIT,
        evidence_bytes=captured,
    )
    if (
        plan.sha256 != ORIGINAL_PLAN
        or plan.campaign_id != ORIGINAL_CAMPAIGN
        or plan.branch != "GR-0"
        or plan.amplitude != "3"
    ):
        raise CFLContinuationAuthorityError("original progression plan differs")

    return CFLContinuationReceipt(
        continuation_commit=continuation_commit,
        continuation_manifest_sha256=current.manifest_sha256,
        continuation_config_sha256=_sha(config_raw),
        continuation_result_sha256=_sha(result_raw),
        original_authorization_commit=ORIGINAL_COMMIT,
        original_plan_sha256=ORIGINAL_PLAN,
        campaign_id=ORIGINAL_CAMPAIGN,
        recovery_checkpoint_generation=GENERATION,
        recovery_checkpoint_sha256=CHECKPOINT_SHA256,
        recovery_journal_sequence=JOURNAL_SEQUENCE,
        recovery_journal_tip_sha256=JOURNAL_TIP_SHA256,
        recovery_member_key=MEMBER_KEY,
        recovery_member_descriptor_sha256=MEMBER_DESCRIPTOR_SHA256,
        authority_paths=current.authority_paths,
        environment=current.environment,
        progression_plan=plan,
    )


__all__ = [
    "ARTIFACT_ID",
    "CHECKPOINT_SHA256",
    "CFLContinuationAuthorityError",
    "CFLContinuationReceipt",
    "CONFIG_PATH",
    "JOURNAL_TIP_SHA256",
    "ORIGINAL_CAMPAIGN",
    "ORIGINAL_COMMIT",
    "ORIGINAL_PLAN",
    "RESULT_PATH",
    "authorize_cfl_continuation",
    "build_cfl1_result",
    "expected_store_anchor",
]
