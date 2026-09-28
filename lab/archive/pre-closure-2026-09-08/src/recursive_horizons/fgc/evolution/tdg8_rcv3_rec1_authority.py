"""Exact recovery authority for the interrupted RCV3 event-23 continuation.

REC1 owns one narrow transition: authenticate the immutable generation-nine
checkpoint plus its exact TDG6/cursor suffix, reconcile that suffix to one
predicted generation-ten checkpoint, close the recovery lease, and only then
hand the unchanged event-23 plan to the already sealed event engine.  Compact
verification is store-blind; live inspection is explicit and read-only.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import tomllib
from types import MappingProxyType
from typing import Any, Iterable, Mapping

from . import proto15_runtime as p15
from . import tdg8_rcv3_execution_authority as auth1
from . import tdg8_rcv3_pref1_binder as pref1
from .hlt16_campaign_store import HLT16CampaignStore, HLT16CampaignStoreError
from .numerical_engine import array_content_sha256
from .proto17_pure_construction import MEMBER_KEYS


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG8-RCV3-REC1-AUTH1"
PROJECT_VERSION = "0.11.0"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CLASSIFICATION = "premise_only_exact_recovery_then_same_event_continuation_authority"
CONFIG_PATH = "configs/fgc/fgc-1-tdg8-rcv3-rec1-auth1.toml"
RESULT_PATH = "results/fgc-1-tdg8-rcv3-rec1-auth1.json"
OWNER_DOCUMENT = "docs/fgc-tdg8-rcv3-rec1-auth1.md"

PRIOR_AUTH1_COMMIT = "afe855be31fe983e940479bffa669e9e60a26344"
PRIOR_AUTH1_CONFIG_PATH = auth1.CONFIG_PATH
PRIOR_AUTH1_CONFIG_SHA256 = "0a453fe3b2ea61a6bfdeea1140228e7028cef668e80e34a2618caeef7d69cbf2"
PRIOR_AUTH1_RESULT_PATH = auth1.RESULT_PATH
PRIOR_AUTH1_RESULT_SHA256 = "6796474ba6f39bd5b425fe564cc758587bbe75eb2ed8a5566a5696af828f3618"
NORMALIZATION_FIX_COMMIT = "ee888156b226c8210a902ed733fd144e0ba24598"
FIXED_SCHEMA_PATH = "src/recursive_horizons/fgc/evolution/hlt16_campaign_schema.py"
FIXED_SCHEMA_SHA256 = "2b6884aaedecfeeaf82a633f92c3418ef0b5889fc620b90faa84bdc714694d48"

PROJECTION_ID = auth1.PROJECTION_ID
DESTINATION_WRAPPER = auth1.DESTINATION_WRAPPER
DESTINATION_PATH = auth1.DESTINATION_PATH
RECEIPT_PATH = auth1.RECEIPT_PATH
RECEIPT_SHA256 = auth1.RECEIPT_SHA256
RECEIPT_RAW_SHA256 = auth1.RECEIPT_RAW_SHA256
CAMPAIGN_ID = auth1.CAMPAIGN_ID
ORIGINAL_EXECUTION_COMMIT = auth1.ORIGINAL_EXECUTION_COMMIT
ORIGINAL_PLAN_SHA256 = auth1.ORIGINAL_PLAN_SHA256

ENTRY_LEAF_COUNT = 54
ENTRY_BYTE_COUNT = 6_328_040
ENTRY_MANIFEST_SHA256 = "774ea5aad22a52c031abf8e6ee8462c8066d0fd9570638639cad77036752043a"
ENTRY_DIRECTORY_COUNTS = {
    "checkpoints": 10, "journal": 13, "locks": 2,
    "payloads": 11, "receipts": 1, "states": 17,
}
ENTRY_CHECKPOINT_GENERATION = 9
ENTRY_CHECKPOINT_SHA256 = "eb6fddc480c94aa7ed15c80399fc2b26637693efef02399f267075fc6c258e56"
ENTRY_CHECKPOINT_RAW_SHA256 = "07572d4450829ff2e2bad65e217a917b0d31ba677a3b8fc20a3be43e849ba846"
ENTRY_JOURNAL_SEQUENCE = 10
ENTRY_JOURNAL_SHA256 = "5b533eb7009a9c7c3f353d26813f9cbbde8574afae77492929b435b4c701c5ff"
ENTRY_MEMBER_KEY = "RK4-2049"
ENTRY_DESCRIPTOR_SHA256 = "77847126340e78c8ac300fac2795bceda724c15e0894e34444b0da42a84166b4"
ENTRY_PHYSICAL_STATE_SHA256 = "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a"

SEQ11_PATH = "journal/00000000000000000011-341cd8cd328434d85774bcb452889ac1886c520a72e945a92337b72cf564161a.journal"
SEQ11_SHA256 = "341cd8cd328434d85774bcb452889ac1886c520a72e945a92337b72cf564161a"
SEQ11_RAW_SHA256 = "e4f9ca42f2d014d4bdbd875a0f24a9bbb9c04159fb0d9e0e7de8d0d72db61fdf"
SEQ11_BYTES = 23_839
SEQ12_PATH = "journal/00000000000000000012-d25d371b67cd227edc735479f719970b382443d9097ee718a4dd3888c0b6cca4.journal"
SEQ12_SHA256 = "d25d371b67cd227edc735479f719970b382443d9097ee718a4dd3888c0b6cca4"
SEQ12_RAW_SHA256 = "dd02233dd5e4cf28becbc4fe29259d0879005d806b70493aa2281255081eb219"
SEQ12_BYTES = 101_869

ENTRY_WRITER_PATH = "locks/active-write.lock"
ENTRY_WRITER_RAW_SHA256 = "278f0c5aa4cfd86661e8bbaef44e340e39592c9b0206cc8b4dc216f387f29605"
ENTRY_WRITER_BYTES = 394
ENTRY_WRITER_HOST = "DouglasMac.local"
ENTRY_WRITER_PID = 57_763
ENTRY_WRITER_SCHEMA = "FGC-1-HLT16-active-writer-v1"
ENTRY_RETIRED_WRITER_PATH = f"locks/.active-write.lock.hlt16-quarantine-{ENTRY_WRITER_RAW_SHA256}"

PREDICTED_GENERATION = 10
PREDICTED_CHECKPOINT_SHA256 = "6dc263e7719c9a422e9fed1b81ac573f3127607a2285f3c55b095d9019f60187"
PREDICTED_CHECKPOINT_RAW_SHA256 = "9abb59999809d22354d6b5f810bc592cd2673e5e4d354d0ec58e9ad53522dc4b"
PREDICTED_CHECKPOINT_BYTES = 186_497
PREDICTED_RETRY_DEPTH = 3
PREDICTED_PENDING_CAP_HEX = "0x1.aaa9612df8000p-12"
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
    ("authority", "src/recursive_horizons/fgc/evolution/tdg8_rcv3_rec1_authority.py"),
    ("runner", "scripts/run_fgc_tdg8_rcv3_rec1_event.py"),
    ("reproducer", "scripts/reproduce_fgc_tdg8_rcv3_rec1_auth1.py"),
    ("predecessor_authority", "src/recursive_horizons/fgc/evolution/tdg8_rcv3_execution_authority.py"),
    ("predecessor_authority", "src/recursive_horizons/fgc/evolution/tdg8_bounded_retry_authority.py"),
    ("predecessor_authority", "src/recursive_horizons/fgc/evolution/tdg8_retry_recovery_authority.py"),
    ("predecessor_binder", "src/recursive_horizons/fgc/evolution/tdg8_rcv3_pref1_binder.py"),
    ("projection_runtime", "src/recursive_horizons/fgc/evolution/tdg8_rcv3_fork_runtime.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_campaign_recovery.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_campaign_runtime.py"),
    ("runtime", FIXED_SCHEMA_PATH),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_campaign_store.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_lifecycle.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/hlt16_progression_attempt.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/proto19_progression_contract.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/proto19_gr0_static_factory.py"),
    ("runtime", "src/recursive_horizons/fgc/evolution/tdg8_successor_runtime.py"),
    ("orchestrator", "scripts/run_fgc_pro19_event1.py"),
    ("owner_document", OWNER_DOCUMENT),
)

AUTHORITY_DELTA_PATHS = tuple(sorted({
    "Makefile", "README.md", CONFIG_PATH, "docs/claim-ledger.md",
    "docs/fgc-runtime-matrix.md", OWNER_DOCUMENT, "docs/research-roadmap.md",
    "paper/fgc-local-defocusing/README.md", "results/README.md", RESULT_PATH,
    "scripts/check_repo.py", "scripts/reproduce_fgc_tdg8_rcv3_rec1_auth1.py",
    "scripts/run_fgc_tdg8_rcv3_rec1_event.py",
    "src/recursive_horizons/fgc/evolution/tdg8_rcv3_rec1_authority.py",
    "tests/test_check_repo_tdg8_rcv3_rec1_auth1.py",
    "tests/test_fgc_tdg8_rcv3_rec1_authority.py",
    "tests/test_fgc_tdg8_rcv3_rec1_event_runner.py",
}))

NONCLAIMS = (
    "REC1 authorizes only exact metadata reconciliation followed by the unchanged GR-0 event-23 continuation; it does not itself complete recovery, an event, trapping calibration, or a physical result.",
    "The generation-nine physical arrays, campaign identity, internal execution commit, plan, thresholds, retry ceiling, minimum width, source/CFL ownership, and scientific stops remain unchanged.",
    "The predicted generation-ten checkpoint advances only retry metadata: no accepted physical state, descriptor, payload, or accepted time changes.",
    "A recovery, writer, provenance, retry, or runtime failure is not a GR-0, FGC-QR, gradient, mechanism, or physical obstruction.",
    "REC1 has no generation-zero, bootstrap, SGB-L, FGC-QR, calibration-promotion, retained-EFT, or physical-transition entry point.",
)
SCOPE = {
    "status_read_only": True,
    "compact_verifier_store_blind": True,
    "exact_seq11_seq12_entry_required": True,
    "exact_predicted_generation10_required_before_PDE": True,
    "same_campaign_no_fork": True,
    "same_event_gr0_resume": True,
    "generation_zero_reimport": False,
    "bootstrap_reexecution": False,
    "candidate_branches_forbidden": True,
    "threshold_or_retry_rule_change": False,
    "physical_input_change": False,
}
CLAIMS = {
    "prior_AUTH1_and_normalization_fix_bound": True,
    "exact_recovery_reconciliation_authorized": True,
    "bounded_event23_continuation_authorized": True,
    "unchanged_internal_plan_preserved": True,
    "candidate_path_excluded": True,
    "recovery_completed": False,
    "event24_completed": False,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "mechanism_result_earned": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
}


class TDG8RCV3REC1AuthorityError(RuntimeError):
    """The immutable recovery contract or its exact live lineage differs."""


class REC1Phase(str, Enum):
    EXACT_PRE_RECOVERY = "exact_pre_recovery"
    GEN9_RECONCILE_PENDING = "generation9_reconcile_pending"
    GEN10_HANDOFF_PENDING = "generation10_handoff_pending"
    EXACT_CLEAN_GEN10 = "exact_clean_generation10"
    LAWFUL_EVENT23_DESCENDANT = "lawful_event23_descendant"
    EVENT24_COMPLETE = "event24_complete"
    TYPED_TERMINAL = "typed_terminal"


@dataclass(frozen=True, slots=True)
class TDG8RCV3REC1Authority:
    authority_commit: str
    prior_auth1_commit: str
    normalization_fix_commit: str
    projection_id: str
    original_execution_commit: str
    original_plan_sha256: str
    campaign_id: str
    destination_path: str
    entry_checkpoint_sha256: str
    predicted_checkpoint_sha256: str
    config_sha256: str
    result_sha256: str
    implementation_inventory: tuple[tuple[str, str, str], ...]
    environment: Mapping[str, str]
    candidate_branches_forbidden: bool = True


@dataclass(frozen=True, slots=True)
class REC1Boundary:
    phase: REC1Phase
    checkpoint_generation: int
    checkpoint_sha256: str
    journal_sequence: int
    journal_tip_sha256: str
    event: int
    target: Mapping[str, str]
    disposition: str
    suffix_kinds: tuple[str, ...]
    active_write: bool
    writer_state: str | None
    writer_checkpoint_sha256: str | None
    retired_writer_count: int
    exact_generation10_ancestor: bool
    physical_state_unchanged_at_recovery: bool
    safe_to_continue: bool
    terminal_lock_present: bool


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
                      allow_nan=False).encode("ascii")


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _decode_json(raw: bytes, label: str) -> dict[str, Any]:
    def unique(items: Iterable[tuple[str, object]]) -> dict[str, object]:
        answer: dict[str, object] = {}
        for key, value in items:
            if key in answer:
                raise ValueError(key)
            answer[key] = value
        return answer
    try:
        value = json.loads(raw.decode("utf-8"), object_pairs_hook=unique)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TDG8RCV3REC1AuthorityError(f"{label} is malformed") from exc
    if not isinstance(value, dict):
        raise TDG8RCV3REC1AuthorityError(f"{label} is not object-valued")
    return value


def _expected_predecessor() -> dict[str, object]:
    return {
        "auth1_authority_commit": PRIOR_AUTH1_COMMIT,
        "auth1_config_path": PRIOR_AUTH1_CONFIG_PATH,
        "auth1_config_sha256": PRIOR_AUTH1_CONFIG_SHA256,
        "auth1_result_path": PRIOR_AUTH1_RESULT_PATH,
        "auth1_result_sha256": PRIOR_AUTH1_RESULT_SHA256,
        "normalization_fix_commit": NORMALIZATION_FIX_COMMIT,
        "fixed_schema_path": FIXED_SCHEMA_PATH,
        "fixed_schema_sha256": FIXED_SCHEMA_SHA256,
    }


def _expected_projection() -> dict[str, object]:
    return {
        "projection_id": PROJECTION_ID, "wrapper_path": DESTINATION_WRAPPER,
        "store_path": DESTINATION_PATH, "receipt_path": RECEIPT_PATH,
        "receipt_sha256": RECEIPT_SHA256, "receipt_raw_sha256": RECEIPT_RAW_SHA256,
    }


def _expected_recovery_entry() -> dict[str, object]:
    return {
        "leaf_count": ENTRY_LEAF_COUNT, "byte_count": ENTRY_BYTE_COUNT,
        "manifest_sha256": ENTRY_MANIFEST_SHA256,
        "directory_leaf_counts": dict(ENTRY_DIRECTORY_COUNTS),
        "checkpoint_generation": ENTRY_CHECKPOINT_GENERATION,
        "checkpoint_sha256": ENTRY_CHECKPOINT_SHA256,
        "checkpoint_raw_sha256": ENTRY_CHECKPOINT_RAW_SHA256,
        "journal_sequence": ENTRY_JOURNAL_SEQUENCE,
        "journal_sha256": ENTRY_JOURNAL_SHA256,
        "suffix_kinds": ["tdg6_rejection", "cursor_transition"],
        "seq11_path": SEQ11_PATH, "seq11_sha256": SEQ11_SHA256,
        "seq11_raw_sha256": SEQ11_RAW_SHA256, "seq11_byte_count": SEQ11_BYTES,
        "seq12_path": SEQ12_PATH, "seq12_sha256": SEQ12_SHA256,
        "seq12_raw_sha256": SEQ12_RAW_SHA256, "seq12_byte_count": SEQ12_BYTES,
        "writer_path": ENTRY_WRITER_PATH,
        "writer_raw_sha256": ENTRY_WRITER_RAW_SHA256,
        "writer_byte_count": ENTRY_WRITER_BYTES,
        "writer_schema": ENTRY_WRITER_SCHEMA,
        "writer_authorization_commit": ORIGINAL_EXECUTION_COMMIT,
        "writer_plan_sha256": ORIGINAL_PLAN_SHA256,
        "writer_checkpoint_sha256": ENTRY_CHECKPOINT_SHA256,
        "writer_host": ENTRY_WRITER_HOST, "writer_pid": ENTRY_WRITER_PID,
    }


def _expected_predicted_recovery() -> dict[str, object]:
    return {
        "generation": PREDICTED_GENERATION,
        "checkpoint_sha256": PREDICTED_CHECKPOINT_SHA256,
        "checkpoint_raw_sha256": PREDICTED_CHECKPOINT_RAW_SHA256,
        "checkpoint_byte_count": PREDICTED_CHECKPOINT_BYTES,
        "parent_sha256": ENTRY_CHECKPOINT_SHA256,
        "journal_sequence": 12, "journal_tip_sha256": SEQ12_SHA256,
        "event": EVENT, "target_rational": TARGET["rational"],
        "target_binary64_hex": TARGET["binary64_hex"],
        "disposition": "nonterminal", "retry_member_key": ENTRY_MEMBER_KEY,
        "descriptor_sha256": ENTRY_DESCRIPTOR_SHA256,
        "physical_state_sha256": ENTRY_PHYSICAL_STATE_SHA256,
        "retry_depth": PREDICTED_RETRY_DEPTH,
        "pending_cap_binary64_hex": PREDICTED_PENDING_CAP_HEX,
        "new_state_count": 0, "new_payload_count": 0,
    }


def _expected_bounded_policy() -> dict[str, object]:
    return {
        "two_phase_reconcile_then_continue": True,
        "exact_generation10_barrier_before_PDE": True,
        "temporal_retry_cap_per_macro_step": TEMPORAL_RETRY_CAP,
        "minimum_macro_step_binary64_hex": MINIMUM_MACRO_STEP_HEX,
        "same_event_only": True, "event_successor": SUCCESSOR_EVENT,
        "event_successor_target_rational": SUCCESSOR_TARGET["rational"],
        "event_successor_target_binary64_hex": SUCCESSOR_TARGET["binary64_hex"],
        "source_CFL_and_scientific_stops_unchanged": True,
        "owner_token_tracked": False,
    }


def _parse_config(raw: bytes) -> dict[str, Any]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise TDG8RCV3REC1AuthorityError("REC1 config is malformed") from exc
    expected_fields = {
        "schema_version", "artifact_id", "project_version", "target_protocol",
        "classification", "nonclaims", "predecessor", "projection",
        "recovery_entry", "predicted_recovery", "bounded_policy", "scope",
        "claims", "implementation_inventory",
    }
    if not isinstance(value, dict) or set(value) != expected_fields:
        raise TDG8RCV3REC1AuthorityError("REC1 config fields differ")
    if (
        value["schema_version"] != SCHEMA_VERSION or value["artifact_id"] != ARTIFACT_ID
        or value["project_version"] != PROJECT_VERSION
        or value["target_protocol"] != TARGET_PROTOCOL
        or value["classification"] != CLASSIFICATION
        or value["nonclaims"] != list(NONCLAIMS)
        or value["predecessor"] != _expected_predecessor()
        or value["projection"] != _expected_projection()
        or value["recovery_entry"] != _expected_recovery_entry()
        or value["predicted_recovery"] != _expected_predicted_recovery()
        or value["bounded_policy"] != _expected_bounded_policy()
        or value["scope"] != SCOPE or value["claims"] != CLAIMS
    ):
        raise TDG8RCV3REC1AuthorityError("REC1 config contract differs")
    rows = value["implementation_inventory"]
    if (
        not isinstance(rows, list) or len(rows) != len(IMPLEMENTATION_INVENTORY)
        or any(not isinstance(row, Mapping) or set(row) != {"role", "path"} for row in rows)
        or tuple((row["role"], row["path"]) for row in rows) != IMPLEMENTATION_INVENTORY
    ):
        raise TDG8RCV3REC1AuthorityError("REC1 implementation inventory differs")
    return value


def _inventory(root: Path) -> list[dict[str, str]]:
    answer: list[dict[str, str]] = []
    for role, relative in IMPLEMENTATION_INVENTORY:
        raw = auth1.rcv2.rcv1._read_leaf(root, relative, f"REC1 {role}", _MAX_BOUND_BYTES)
        answer.append({"role": role, "path": relative, "sha256": _sha(raw)})
    return answer


def _payload(config: Mapping[str, Any], inventory: list[dict[str, str]],
             environment: Mapping[str, str]) -> dict[str, Any]:
    return {
        "predecessor": dict(config["predecessor"]),
        "projection": dict(config["projection"]),
        "recovery_entry": dict(config["recovery_entry"]),
        "predicted_recovery": dict(config["predicted_recovery"]),
        "bounded_policy": dict(config["bounded_policy"]),
        "implementation_inventory": inventory, "environment": dict(environment),
        "scope": dict(SCOPE), "claims": dict(CLAIMS), "nonclaims": list(NONCLAIMS),
    }


def _validate_predecessors(repository: Path) -> None:
    for commit, path, expected in (
        (PRIOR_AUTH1_COMMIT, PRIOR_AUTH1_CONFIG_PATH, PRIOR_AUTH1_CONFIG_SHA256),
        (PRIOR_AUTH1_COMMIT, PRIOR_AUTH1_RESULT_PATH, PRIOR_AUTH1_RESULT_SHA256),
        (NORMALIZATION_FIX_COMMIT, FIXED_SCHEMA_PATH, FIXED_SCHEMA_SHA256),
    ):
        live = auth1.rcv2.rcv1._read_leaf(repository, path, f"REC1 predecessor {path}", _MAX_BOUND_BYTES)
        if _sha(live) != expected or auth1.rcv2.rcv1._committed_bytes(repository, commit, path) != live:
            raise TDG8RCV3REC1AuthorityError(f"REC1 predecessor bytes differ: {path}")


def build_prelaunch(config_raw: bytes, root: Path) -> dict[str, Any]:
    """Build premise-only evidence after authenticating the exact live entry."""
    repository = auth1.rcv2.rcv1._repository(root)
    config = _parse_config(config_raw)
    if auth1.rcv2.rcv1._read_leaf(repository, CONFIG_PATH, "REC1 config", _MAX_CONFIG_BYTES) != config_raw:
        raise TDG8RCV3REC1AuthorityError("live REC1 config bytes differ")
    try:
        (repository / RESULT_PATH).lstat()
    except FileNotFoundError:
        pass
    else:
        raise TDG8RCV3REC1AuthorityError("REC1 result already exists")
    _validate_predecessors(repository)
    auth1._authenticate_external_receipt(repository)
    # A synthetic receipt is enough here because inspect_store validates only
    # the fixed typed identities, never Git or compact result bytes.
    receipt = TDG8RCV3REC1Authority(
        authority_commit="0" * 40, prior_auth1_commit=PRIOR_AUTH1_COMMIT,
        normalization_fix_commit=NORMALIZATION_FIX_COMMIT, projection_id=PROJECTION_ID,
        original_execution_commit=ORIGINAL_EXECUTION_COMMIT,
        original_plan_sha256=ORIGINAL_PLAN_SHA256, campaign_id=CAMPAIGN_ID,
        destination_path=DESTINATION_PATH, entry_checkpoint_sha256=ENTRY_CHECKPOINT_SHA256,
        predicted_checkpoint_sha256=PREDICTED_CHECKPOINT_SHA256,
        config_sha256=_sha(config_raw), result_sha256="0" * 64,
        implementation_inventory=(), environment=MappingProxyType({}),
    )
    boundary = inspect_store(repository / DESTINATION_PATH, receipt)
    if boundary.phase is not REC1Phase.EXACT_PRE_RECOVERY:
        raise TDG8RCV3REC1AuthorityError("REC1 live entry is not exact pre-recovery")
    return {
        "schema_version": SCHEMA_VERSION, "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION, "target_protocol": TARGET_PROTOCOL,
        "classification": CLASSIFICATION, "gate_status": "pass",
        "source_config_sha256": _sha(config_raw),
        "artifact_payload": _payload(config, _inventory(repository),
                                     auth1.rcv2.rcv1.observed_environment()),
    }


def _parse_result(raw: bytes) -> dict[str, Any]:
    result = _decode_json(raw, "REC1 result")
    payload = result.get("artifact_payload")
    rows = payload.get("implementation_inventory") if isinstance(payload, Mapping) else None
    if not isinstance(rows, list) or len(rows) != len(IMPLEMENTATION_INVENTORY):
        raise TDG8RCV3REC1AuthorityError("REC1 result inventory differs")
    for row, expected in zip(rows, IMPLEMENTATION_INVENTORY, strict=True):
        if (
            not isinstance(row, Mapping) or set(row) != {"role", "path", "sha256"}
            or (row.get("role"), row.get("path")) != expected
            or not isinstance(row.get("sha256"), str) or not _SHA256.fullmatch(str(row["sha256"]))
        ):
            raise TDG8RCV3REC1AuthorityError("REC1 result inventory row differs")
    return result


def validate_compact(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    config = _parse_config(config_raw)
    result = _parse_result(result_raw)
    payload = result["artifact_payload"]
    expected = {
        "schema_version": SCHEMA_VERSION, "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION, "target_protocol": TARGET_PROTOCOL,
        "classification": CLASSIFICATION, "gate_status": "pass",
        "source_config_sha256": _sha(config_raw),
        "artifact_payload": _payload(config, [dict(row) for row in payload["implementation_inventory"]],
                                     payload.get("environment", {})),
    }
    if canonical(result) != result_raw or canonical(expected) != result_raw:
        raise TDG8RCV3REC1AuthorityError("REC1 compact result differs")
    return result


def _require_git_image(root: Path, commit: str) -> None:
    if not _COMMIT.fullmatch(commit):
        raise TDG8RCV3REC1AuthorityError("REC1 authority commit is malformed")
    head = auth1.rcv2.rcv1._git(root, "rev-parse", "--verify", "HEAD").stdout.decode("ascii").strip()
    resolved = auth1.rcv2.rcv1._git(root, "rev-parse", "--verify", f"{commit}^{{commit}}").stdout.decode("ascii").strip()
    if head != commit or resolved != commit:
        raise TDG8RCV3REC1AuthorityError("REC1 authority commit differs from HEAD")
    for ancestor in (PRIOR_AUTH1_COMMIT, NORMALIZATION_FIX_COMMIT):
        check = subprocess.run(
            ("git", "--no-replace-objects", "--no-optional-locks", "-C", str(root),
             "merge-base", "--is-ancestor", ancestor, commit),
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            check=False, env={**{k: v for k, v in os.environ.items() if not k.startswith("GIT_")},
                              "GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C"},
        )
        if check.returncode != 0:
            raise TDG8RCV3REC1AuthorityError("REC1 predecessor commit is not an ancestor")
    changed = auth1.rcv2.rcv1._git(
        root, "diff", "--name-only", "-z", "--no-renames",
        f"{NORMALIZATION_FIX_COMMIT}..{commit}", "--",
    ).stdout
    try:
        paths = tuple(sorted(item.decode("utf-8") for item in changed.split(b"\0") if item))
    except UnicodeDecodeError as exc:
        raise TDG8RCV3REC1AuthorityError("REC1 Git delta is not UTF-8") from exc
    if paths != AUTHORITY_DELTA_PATHS:
        raise TDG8RCV3REC1AuthorityError("REC1 authority Git delta differs")


def authorize_execution(root: Path, config_raw: bytes, result_raw: bytes,
                        authority_commit: str) -> TDG8RCV3REC1Authority:
    repository = auth1.rcv2.rcv1._repository(root)
    result = validate_compact(config_raw, result_raw)
    _require_git_image(repository, authority_commit)
    _validate_predecessors(repository)
    for path, expected in ((CONFIG_PATH, config_raw), (RESULT_PATH, result_raw)):
        live = auth1.rcv2.rcv1._read_leaf(repository, path, f"REC1 committed {path}", _MAX_RESULT_BYTES)
        if live != expected or auth1.rcv2.rcv1._committed_bytes(repository, authority_commit, path) != live:
            raise TDG8RCV3REC1AuthorityError(f"REC1 committed bytes differ: {path}")
    rows: list[tuple[str, str, str]] = []
    for item in result["artifact_payload"]["implementation_inventory"]:
        role, path, digest = str(item["role"]), str(item["path"]), str(item["sha256"])
        live = auth1.rcv2.rcv1._read_leaf(repository, path, f"REC1 implementation {path}", _MAX_BOUND_BYTES)
        if _sha(live) != digest or auth1.rcv2.rcv1._committed_bytes(repository, authority_commit, path) != live:
            raise TDG8RCV3REC1AuthorityError(f"REC1 implementation bytes differ: {path}")
        rows.append((role, path, digest))
    environment = result["artifact_payload"].get("environment")
    if environment != auth1.rcv2.rcv1.observed_environment():
        raise TDG8RCV3REC1AuthorityError("REC1 execution environment differs")
    return TDG8RCV3REC1Authority(
        authority_commit=authority_commit, prior_auth1_commit=PRIOR_AUTH1_COMMIT,
        normalization_fix_commit=NORMALIZATION_FIX_COMMIT, projection_id=PROJECTION_ID,
        original_execution_commit=ORIGINAL_EXECUTION_COMMIT,
        original_plan_sha256=ORIGINAL_PLAN_SHA256, campaign_id=CAMPAIGN_ID,
        destination_path=DESTINATION_PATH, entry_checkpoint_sha256=ENTRY_CHECKPOINT_SHA256,
        predicted_checkpoint_sha256=PREDICTED_CHECKPOINT_SHA256,
        config_sha256=_sha(config_raw), result_sha256=_sha(result_raw),
        implementation_inventory=tuple(rows), environment=MappingProxyType(dict(environment)),
    )


def _require_receipt(receipt: TDG8RCV3REC1Authority) -> None:
    if (
        not isinstance(receipt, TDG8RCV3REC1Authority)
        or receipt.prior_auth1_commit != PRIOR_AUTH1_COMMIT
        or receipt.normalization_fix_commit != NORMALIZATION_FIX_COMMIT
        or receipt.projection_id != PROJECTION_ID
        or receipt.original_execution_commit != ORIGINAL_EXECUTION_COMMIT
        or receipt.original_plan_sha256 != ORIGINAL_PLAN_SHA256
        or receipt.campaign_id != CAMPAIGN_ID or receipt.destination_path != DESTINATION_PATH
        or receipt.entry_checkpoint_sha256 != ENTRY_CHECKPOINT_SHA256
        or receipt.predicted_checkpoint_sha256 != PREDICTED_CHECKPOINT_SHA256
        or not receipt.candidate_branches_forbidden
    ):
        raise TDG8RCV3REC1AuthorityError("REC1 typed authority differs")


def _tree(path: Path):
    try:
        return pref1._snapshot(path.parent, path.name)
    except Exception as exc:
        raise TDG8RCV3REC1AuthorityError("REC1 store tree is unsafe") from exc


def _leaf(tree: Any, relative: str, raw_sha256: str, byte_count: int) -> None:
    item = tree.leaves.get(relative)
    if item is None or item.raw_sha256 != raw_sha256 or item.byte_count != byte_count:
        raise TDG8RCV3REC1AuthorityError(f"REC1 exact leaf differs: {relative}")


def _directory_counts(tree: Any) -> dict[str, int]:
    return {
        name: sum(1 for path in tree.leaves if path.startswith(f"{name}/"))
        for name in ENTRY_DIRECTORY_COUNTS
    }


def _only_recovery_locks(tree: Any, *, active_required: bool) -> bool:
    paths = {path for path in tree.leaves if path.startswith("locks/")}
    fixed = {"locks/writer.guard"}
    if active_required:
        fixed.add(ENTRY_WRITER_PATH)
    quarantines = {
        path for path in paths
        if path.startswith("locks/.active-write.lock.hlt16-quarantine-")
        and _SHA256.fullmatch(path.rsplit("-", 1)[-1])
    }
    return paths == fixed | quarantines


def _writer(store: HLT16CampaignStore, status: Any) -> Mapping[str, Any] | None:
    try:
        value = store._active_writer()
    except Exception as exc:
        raise TDG8RCV3REC1AuthorityError("REC1 writer lease cannot be decoded") from exc
    if (value is None) != (not status.active_write):
        raise TDG8RCV3REC1AuthorityError("REC1 writer status differs")
    return value


def _writer_authority(value: Mapping[str, Any], trusted: set[str]) -> None:
    if (
        value.get("schema") != ENTRY_WRITER_SCHEMA
        or value.get("authorization_commit") != ORIGINAL_EXECUTION_COMMIT
        or value.get("plan_sha256") != ORIGINAL_PLAN_SHA256
        or value.get("checkpoint_sha256") not in trusted
        or not isinstance(value.get("host"), str) or not value.get("host")
        or isinstance(value.get("pid"), bool) or not isinstance(value.get("pid"), int)
    ):
        raise TDG8RCV3REC1AuthorityError("REC1 writer authority differs")


def _physical_entry(store: HLT16CampaignStore) -> None:
    try:
        state = store.load_state(ENTRY_DESCRIPTOR_SHA256)
        digest = array_content_sha256(state.arrays["u"], state.arrays["p"], state.arrays["q"])
    except Exception as exc:
        raise TDG8RCV3REC1AuthorityError("REC1 entry physical state cannot be restored") from exc
    if digest != ENTRY_PHYSICAL_STATE_SHA256:
        raise TDG8RCV3REC1AuthorityError("REC1 entry physical state differs")


def _predicted_checkpoint(tree: Any, checkpoint: Any) -> None:
    relative = f"checkpoints/{PREDICTED_GENERATION:020d}-{PREDICTED_CHECKPOINT_SHA256}.json"
    _leaf(tree, relative, PREDICTED_CHECKPOINT_RAW_SHA256, PREDICTED_CHECKPOINT_BYTES)
    member = checkpoint.members.get(ENTRY_MEMBER_KEY)
    retry = member.cursor.get("retry_successor_payload_or_none") if member is not None else None
    if (
        checkpoint.generation != PREDICTED_GENERATION
        or checkpoint.sha256 != PREDICTED_CHECKPOINT_SHA256
        or checkpoint.parent_sha256 != ENTRY_CHECKPOINT_SHA256
        or checkpoint.journal_sequence != 12 or checkpoint.journal_tip_sha256 != SEQ12_SHA256
        or checkpoint.event != EVENT or dict(checkpoint.target) != TARGET
        or checkpoint.disposition != "nonterminal" or member is None
        or member.descriptor_sha256 != ENTRY_DESCRIPTOR_SHA256
        or member.pending_owner != "temporal" or member.pending_cap_hex != PREDICTED_PENDING_CAP_HEX
        or not isinstance(retry, Mapping) or retry.get("retry_count") != PREDICTED_RETRY_DEPTH
    ):
        raise TDG8RCV3REC1AuthorityError("REC1 predicted generation-ten checkpoint differs")


def inspect_store(destination: Path, receipt: TDG8RCV3REC1Authority, *,
                  owned_writer_pid: int | None = None) -> REC1Boundary:
    """Classify only the closed REC1 recovery/continuation state machine."""
    _require_receipt(receipt)
    destination = Path(destination)
    before = _tree(destination)
    store = HLT16CampaignStore(destination)
    try:
        snapshot = store.authenticated_snapshot()
        status = store.inspect_recovery()
        _records, checkpoints, _suffix = store._validate_full_store(repair_stages=False)
        closure = store.authenticated_snapshot()
        closure_status = store.inspect_recovery()
    except HLT16CampaignStoreError as exc:
        raise TDG8RCV3REC1AuthorityError("REC1 store does not validate") from exc
    after = _tree(destination)
    if (
        before != after or snapshot != closure or status != closure_status
        or not checkpoints or checkpoints[-1].sha256 != snapshot.checkpoint.sha256
    ):
        raise TDG8RCV3REC1AuthorityError("REC1 store changed during read-only inspection")
    checkpoint = snapshot.checkpoint
    if (
        checkpoint.authorization_commit != ORIGINAL_EXECUTION_COMMIT
        or checkpoint.plan_sha256 != ORIGINAL_PLAN_SHA256
        or checkpoint.campaign_id != CAMPAIGN_ID or checkpoint.protocol != TARGET_PROTOCOL
        or tuple(checkpoint.members) != MEMBER_KEYS
    ):
        raise TDG8RCV3REC1AuthorityError("REC1 internal authority differs")
    trusted = {item.sha256 for item in checkpoints}
    writer = _writer(store, status)
    if writer is not None:
        _writer_authority(writer, trusted)
    retired = tuple(path for path in before.leaves if path.startswith("locks/.active-write.lock.hlt16-quarantine-"))
    old_retired = ENTRY_RETIRED_WRITER_PATH in before.leaves
    if old_retired:
        _leaf(before, ENTRY_RETIRED_WRITER_PATH, ENTRY_WRITER_RAW_SHA256, ENTRY_WRITER_BYTES)
    _physical_entry(store)

    phase: REC1Phase
    exact_gen10 = any(item.generation == PREDICTED_GENERATION and item.sha256 == PREDICTED_CHECKPOINT_SHA256
                      for item in checkpoints)
    physical_unchanged = checkpoint.generation >= PREDICTED_GENERATION
    if checkpoint.generation == ENTRY_CHECKPOINT_GENERATION:
        if (
            checkpoint.sha256 != ENTRY_CHECKPOINT_SHA256
            or checkpoint.journal_sequence != ENTRY_JOURNAL_SEQUENCE
            or checkpoint.journal_tip_sha256 != ENTRY_JOURNAL_SHA256
            or snapshot.suffix_kinds != ("tdg6_rejection", "cursor_transition")
            or snapshot.staging_paths or snapshot.orphan_descriptor_sha256 is not None
            or snapshot.orphan_payload_semantic_sha256 is not None
        ):
            raise TDG8RCV3REC1AuthorityError("REC1 exact generation-nine suffix differs")
        _leaf(before, SEQ11_PATH, SEQ11_RAW_SHA256, SEQ11_BYTES)
        _leaf(before, SEQ12_PATH, SEQ12_RAW_SHA256, SEQ12_BYTES)
        if writer is None:
            raise TDG8RCV3REC1AuthorityError("REC1 generation-nine recovery lease is absent")
        if writer["checkpoint_sha256"] != ENTRY_CHECKPOINT_SHA256:
            raise TDG8RCV3REC1AuthorityError("REC1 generation-nine writer anchor differs")
        exact_entry = (
            len(before.leaves) == ENTRY_LEAF_COUNT and before.byte_count() == ENTRY_BYTE_COUNT
            and before.manifest_sha256() == ENTRY_MANIFEST_SHA256
            and _directory_counts(before) == ENTRY_DIRECTORY_COUNTS
        )
        if exact_entry:
            _leaf(before, ENTRY_WRITER_PATH, ENTRY_WRITER_RAW_SHA256, ENTRY_WRITER_BYTES)
            if (
                writer.get("host") != ENTRY_WRITER_HOST or writer.get("pid") != ENTRY_WRITER_PID
                or status.writer_state != "stale_verified_lock" or old_retired
            ):
                raise TDG8RCV3REC1AuthorityError("REC1 exact stale entry writer differs")
            phase = REC1Phase.EXACT_PRE_RECOVERY
        else:
            # Repeated process death between takeover and reconciliation may
            # append further content-addressed retired leases.  The campaign
            # validator has authenticated every body and authority; accept the
            # monotone evidence chain, not one accidental quarantine count.
            counts = _directory_counts(before)
            if (
                not old_retired or len(retired) < 1
                or any(counts[name] != ENTRY_DIRECTORY_COUNTS[name] for name in counts if name != "locks")
                or counts["locks"] != len(retired) + 2
                or not _only_recovery_locks(before, active_required=True)
            ):
                raise TDG8RCV3REC1AuthorityError("REC1 generation-nine takeover cut differs")
            allowed = status.writer_state == "stale_verified_lock" or (
                owned_writer_pid is not None and writer.get("host") == os.uname().nodename
                and writer.get("pid") == owned_writer_pid and status.writer_state == "live_local_writer"
            )
            if not allowed:
                raise TDG8RCV3REC1AuthorityError("REC1 generation-nine takeover owner differs")
            phase = REC1Phase.GEN9_RECONCILE_PENDING
        physical_unchanged = True
    elif checkpoint.generation >= PREDICTED_GENERATION:
        if not exact_gen10:
            raise TDG8RCV3REC1AuthorityError("REC1 predicted generation-ten ancestor is absent")
        generation10 = next(item for item in checkpoints if item.generation == PREDICTED_GENERATION)
        _predicted_checkpoint(before, generation10)
        if checkpoint.generation == PREDICTED_GENERATION and snapshot.suffix_kinds:
            # A suffix written by the existing event engine is a descendant,
            # never part of metadata reconciliation.
            phase = REC1Phase.LAWFUL_EVENT23_DESCENDANT
        elif checkpoint.generation == PREDICTED_GENERATION and writer is not None and writer["checkpoint_sha256"] == ENTRY_CHECKPOINT_SHA256:
            allowed = status.writer_state == "stale_verified_lock" or (
                owned_writer_pid is not None and writer.get("host") == os.uname().nodename
                and writer.get("pid") == owned_writer_pid and status.writer_state == "live_local_writer"
            )
            if not old_retired or not allowed:
                raise TDG8RCV3REC1AuthorityError("REC1 generation-ten handoff lease differs")
            phase = REC1Phase.GEN10_HANDOFF_PENDING
        elif checkpoint.event == SUCCESSOR_EVENT:
            if dict(checkpoint.target) != SUCCESSOR_TARGET or checkpoint.disposition != "event_complete" or snapshot.suffix_kinds:
                raise TDG8RCV3REC1AuthorityError("REC1 event-24 completion differs")
            phase = REC1Phase.EVENT24_COMPLETE
        elif checkpoint.disposition in {"scientific_terminal", "invalid_terminal"}:
            if checkpoint.event != EVENT or dict(checkpoint.target) != TARGET or not snapshot.terminal_lock_present:
                raise TDG8RCV3REC1AuthorityError("REC1 typed terminal differs")
            phase = REC1Phase.TYPED_TERMINAL
        elif checkpoint.event == EVENT and dict(checkpoint.target) == TARGET:
            if checkpoint.generation == PREDICTED_GENERATION and not snapshot.suffix_kinds and writer is None:
                if (
                    not old_retired or len(retired) < 2
                    or not _only_recovery_locks(before, active_required=False)
                ):
                    raise TDG8RCV3REC1AuthorityError("REC1 clean generation-ten lease ledger differs")
                phase = REC1Phase.EXACT_CLEAN_GEN10
            else:
                phase = REC1Phase.LAWFUL_EVENT23_DESCENDANT
        else:
            raise TDG8RCV3REC1AuthorityError("REC1 descendant escaped one-event scope")
    else:
        raise TDG8RCV3REC1AuthorityError("REC1 checkpoint predates its exact entry")

    safe = phase in {REC1Phase.EXACT_PRE_RECOVERY, REC1Phase.EXACT_CLEAN_GEN10,
                     REC1Phase.LAWFUL_EVENT23_DESCENDANT}
    if writer is not None and status.writer_state != "stale_verified_lock":
        safe = False
    return REC1Boundary(
        phase=phase, checkpoint_generation=checkpoint.generation,
        checkpoint_sha256=checkpoint.sha256, journal_sequence=checkpoint.journal_sequence,
        journal_tip_sha256=checkpoint.journal_tip_sha256, event=checkpoint.event,
        target=dict(checkpoint.target), disposition=checkpoint.disposition,
        suffix_kinds=tuple(snapshot.suffix_kinds), active_write=status.active_write,
        writer_state=status.writer_state,
        writer_checkpoint_sha256=None if writer is None else str(writer["checkpoint_sha256"]),
        retired_writer_count=len(retired), exact_generation10_ancestor=exact_gen10,
        physical_state_unchanged_at_recovery=physical_unchanged,
        safe_to_continue=safe, terminal_lock_present=snapshot.terminal_lock_present,
    )


def inspect_boundary(root: Path, receipt: TDG8RCV3REC1Authority) -> REC1Boundary:
    repository = auth1.rcv2.rcv1._repository(root)
    _require_receipt(receipt)
    auth1._real_directory(repository, DESTINATION_WRAPPER, "REC1 projection wrapper")
    destination = auth1._real_directory(repository, DESTINATION_PATH, "REC1 campaign store")
    auth1._authenticate_external_receipt(repository)
    return inspect_store(destination, receipt)


__all__ = [
    "ARTIFACT_ID", "AUTHORITY_DELTA_PATHS", "CAMPAIGN_ID", "CLASSIFICATION",
    "CONFIG_PATH", "DESTINATION_PATH", "ENTRY_CHECKPOINT_SHA256", "EVENT",
    "IMPLEMENTATION_INVENTORY", "MINIMUM_MACRO_STEP_HEX", "NORMALIZATION_FIX_COMMIT",
    "ORIGINAL_EXECUTION_COMMIT", "ORIGINAL_PLAN_SHA256", "PREDICTED_CHECKPOINT_SHA256",
    "PRIOR_AUTH1_COMMIT", "PROJECT_VERSION", "REC1Boundary", "REC1Phase", "RESULT_PATH",
    "SUCCESSOR_EVENT", "SUCCESSOR_TARGET", "TARGET", "TDG8RCV3REC1Authority",
    "TDG8RCV3REC1AuthorityError", "TEMPORAL_RETRY_CAP", "authorize_execution",
    "build_prelaunch", "canonical", "inspect_boundary", "inspect_store", "validate_compact",
]
