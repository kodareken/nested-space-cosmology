"""Prospective HLT17 SRCQ1-REC1 recovery core.

Implementation A is a deterministic no-write harness. It authenticates the
sealed Attempt2 compact/raw identities, rebuilds the byte-bound origin and
runtime, skips RK member qualification as predecessor evidence, and resumes
at SSPRK3-4097. Source/CFL overlays are absorbed only into an isolated
in-memory checkpoint. The core never adopts an endpoint or writes a campaign
store. Ordinary status inspection does not call the physical source.
"""

from __future__ import annotations

from hashlib import sha256
import importlib
import json
import math
from pathlib import Path
import re
import time
from types import MappingProxyType
from typing import Any, Callable, Mapping

from recursive_horizons.evidence_io import (
    UnsafePathError,
    canonical_json_bytes,
    read_regular_file,
)

from .hlt17_admission_runtime import (
    C1R1CoordinateLatticeStop,
    C1R1RefinementPathStop,
    C1R1RuntimeResourceStop,
    hlt17_implementation_identity,
    require_hlt17_implementation_identity,
)
from .hlt17_imp1_cursor import (
    HLT17_CFL_RETRY_CAP,
    HLT17_SOURCE_RETRY_CAP,
    HLT17OwnerRetryExhausted,
    encode_tdg7_plan,
)
from .hlt17_member_checkpoint import HLT17InMemoryCheckpoint
from .hlt17_member_checkpoint import _isolated_member as isolate_hlt17_runtime_member
from .hlt17_srcq1 import (
    FORBIDDEN_ENDPOINT_KEYS,
    HLT17SRCQ1AuthoritySeamError,
    HLT17SRCQ1CapError,
    HLT17SRCQ1Error,
    HLT17SRCQ1IdentityError,
    HLT17SRCQ1PublicationError,
    HLT17SRCQ1ResourceError,
    HLT17SRCQ1SourceError,
    ResourceCeilings,
    STATIC_INPUT_PATHS,
    STATIC_INPUT_SHA256,
    c1r1_shadow_record,
    check_inherited_counters,
    check_method_owned_identity,
    current_rss_bytes,
    derive_prospective_requested_caps,
    evaluate_bound_source_once,
    member_fingerprints,
    observe_environment,
    publish_qualification_result as publish_srcq1_qualification_result,
    recommended_resource_ceilings,
    refuse_endpoint_payload,
    static_input_identities,
)
from .hlt17_srcq1 import _check_resources as check_srcq1_resources
from .proto5_runtime import CFLRetryRequired
from .pro20_origin import (
    MEMBER_KEYS,
    ORIGIN_TIME_HEX,
    Pro20OriginCapture,
    Pro20OriginError,
    revalidate_pro20_origin_capture,
)
from .pro20_source_factory import (
    CAMPAIGN_ID,
    EVENT_TARGET,
    PROTOCOL,
    Pro20SourceFactoryError,
    build_pro20_runtime_origin,
)


SCHEMA = "FGC-1-HLT17-SRCQ1-REC1-qualification-v1"
ARTIFACT_ID = "FGC-1-HLT17-SRCQ1-REC1"
IMPLEMENTATION_RELATIVE = "src/recursive_horizons/fgc/evolution/hlt17_srcq1_rec1.py"
RESULT_LEAF = "qualification.json"
AUTHORITY_VALIDATOR_MODULE = "recursive_horizons.fgc.evolution.hlt17_srcq1_rec1_auth1"
AUTHORITY_VALIDATOR_NAME = "validate_authority_delta"
AUTHORITY_SEAM = (
    "authority/delta validation is the coordinator freeze-B seam "
    f"{AUTHORITY_VALIDATOR_MODULE}.{AUTHORITY_VALIDATOR_NAME}; "
    "implementation A does not accept a boolean skip/force flag"
)
PREDECESSOR_ARTIFACT_ID = "FGC-1-HLT17-SRCQ1"
ATTEMPT2_ARTIFACT_ID = "FGC-1-HLT17-SRCQ1-ATTEMPT2"
ATTEMPT2_COMPACT_RELATIVE = "results/fgc-1-hlt17-srcq1-attempt2.json"
ATTEMPT2_RAW_LEAF_RELATIVE = (
    "runs/fgc-2-sf1/hlt17-srcq1-source-origin-qualification/qualification.json"
)
ATTEMPT2_SEALED_COMPACT_COMMIT = "efb53b839a36a7f7f52018a44ef4eeaee42c1752"
ATTEMPT2_SEALED_COMPACT_SHA256 = (
    "7860491360fe058d6c5019f780103e6bf13b2c243d3c1d59ec3f74109636ab5f"
)
ATTEMPT2_METADATA_LINKED_COMPACT_COMMIT = (
    "4891feab58922e449a06585a848eb6a11df994ea"
)
ATTEMPT2_METADATA_LINKED_COMPACT_SHA256 = (
    "ba656a63343a87b389e388e39e142e8cfc94194c96c33326f6da0963f4504c48"
)
ATTEMPT2_DERIVATION_DOCUMENT = "docs/fgc-hlt17-srcq1-frz1.md"
ATTEMPT2_COMPACT_SHA256 = ATTEMPT2_SEALED_COMPACT_SHA256
ATTEMPT2_RAW_LEAF_SHA256 = (
    "889a8d17551a1eea6ed56f43f602fefb3f61ab78637bbdd0c995905e22259d89"
)
ATTEMPT2_RAW_DIRECTORY_PAYLOAD_SHA256 = (
    "7b06deb68c989734f8dd669f2f1a8f4531aa90b82cc0e9ff34330db363123009"
)
ATTEMPT2_RAW_LEAF_BYTES = 45589
ATTEMPT2_CLASSIFICATION = (
    "completed_terminal_reports_cfl_evidence_incomplete_no_state_advance"
)
ATTEMPT2_AUTHORITY_COMMIT = "f3a31dc6fce35baa90d305a61e89cee9706fed38"
ATTEMPT2_IMPLEMENTATION_COMMIT = "7123aa385b3f4165186acd241f532b58ebbb3a9b"
ATTEMPT2_IMPLEMENTATION_SHA256 = (
    "e801532e01061122b03bd2e5b8b1359926a65445b676d695fcdc6da371784bab"
)
ORIGIN_CAPTURE_SHA256 = (
    "2824af06c674336d479ec53faa4269f02fb598055860ec22f2ad84602a83fe72"
)
PREDECESSOR_MEMBER_KEYS = ("RK4-2049", "RK4-4097", "RK4-8193")
RESUME_MEMBER_KEYS = ("SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385")
FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER = MappingProxyType(
    {
        "RK4-2049": "0x1.aaa90b0fb5c26p-8",
        "RK4-4097": "0x1.aaaaaa9fefc9cp-9",
        "RK4-8193": "0x1.aaaaaaaaaa654p-10",
        "SSPRK3-4097": "0x1.aaaaaa9908e70p-9",
        "SSPRK3-8193": "0x1.aaaaaaaaa8b26p-10",
        "SSPRK3-16385": "0x1.aaaaaaaaaaa91p-11",
    }
)
PREPARED_DISPOSITIONS = frozenset(
    {"prepared_fresh", "prepared_overlay", "prepared_successor"}
)
CFL_REQUIRED = "cfl_retry_required"
CFL_EXHAUSTED = "cfl_retry_exhausted"
SOURCE_REQUIRED = "source_retry_required"
SOURCE_EXHAUSTED = "source_retry_exhausted"
CFL_DISPOSITIONS = frozenset({CFL_REQUIRED, CFL_EXHAUSTED})
SOURCE_DISPOSITIONS = frozenset({SOURCE_REQUIRED, SOURCE_EXHAUSTED})
TEMPORAL_DISPOSITIONS = frozenset(
    {"temporal_retry_required", "temporal_retry_exhausted"}
)
OVERLAY_CONTINUE = frozenset({CFL_REQUIRED, SOURCE_REQUIRED})
MAX_PREPARE_ATTEMPTS = (HLT17_SOURCE_RETRY_CAP + 1) + (HLT17_CFL_RETRY_CAP + 1)
_COMMIT = re.compile(r"\A[0-9a-f]{40}(?:[0-9a-f]{24})?\Z")
_SHA = re.compile(r"\A[0-9a-f]{64}\Z")
RESOURCE_CEILINGS_ARE_COORDINATOR_FREEZE_VALUES = True


class HLT17SRCQ1REC1Error(HLT17SRCQ1Error):
    """Typed SRCQ1-REC1 stop; not a campaign classification."""

    outcome = "invalid"
    published = False
    phase = "invalid"
    physical_work_started = False


class HLT17SRCQ1REC1CapError(HLT17SRCQ1CapError, HLT17SRCQ1REC1Error):
    outcome = "invalid"


class HLT17SRCQ1REC1SourceError(HLT17SRCQ1SourceError, HLT17SRCQ1REC1Error):
    outcome = "source"


class HLT17SRCQ1REC1ResourceError(HLT17SRCQ1ResourceError, HLT17SRCQ1REC1Error):
    outcome = "resource"


class HLT17SRCQ1REC1IdentityError(HLT17SRCQ1IdentityError, HLT17SRCQ1REC1Error):
    outcome = "invalid"


class HLT17SRCQ1REC1AuthoritySeamError(
    HLT17SRCQ1AuthoritySeamError, HLT17SRCQ1REC1Error
):
    outcome = "invalid"
    phase = "authority_seam"
    physical_work_started = False


class HLT17SRCQ1REC1InconclusiveError(HLT17SRCQ1REC1Error):
    outcome = "inconclusive"


class HLT17SRCQ1REC1PublicationError(HLT17SRCQ1PublicationError, HLT17SRCQ1REC1Error):
    outcome = "invalid"
    phase = "publication"

    def __init__(
        self,
        message: str,
        *,
        published: bool = False,
        physical_work_started: bool = False,
    ) -> None:
        HLT17SRCQ1Error.__init__(self, message)
        self.published = published
        self.physical_work_started = physical_work_started
        if published:
            self.outcome = "inconclusive"


def _sha256_hex(raw: bytes) -> str:
    if type(raw) is not bytes:
        raise HLT17SRCQ1REC1Error("payload is not immutable bytes")
    return sha256(raw).hexdigest()


def _require_sha(value: object, label: str) -> str:
    if type(value) is not str or _SHA.fullmatch(value) is None:
        raise HLT17SRCQ1REC1Error(f"{label} must be a lowercase SHA-256 digest")
    return value


def _require_commit(value: object, label: str) -> str:
    if type(value) is not str or _COMMIT.fullmatch(value) is None:
        raise HLT17SRCQ1REC1AuthoritySeamError(
            f"{label} must be an exact lowercase 40- or 64-digit commit"
        )
    return value


def _mapping(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping) or type(value) is bool:
        raise HLT17SRCQ1REC1Error(f"{label} must be a mapping")
    return dict(value)


def _json_mapping(raw: bytes, label: str) -> dict[str, Any]:
    if type(raw) is not bytes:
        raise HLT17SRCQ1REC1Error(f"{label} is not immutable bytes")
    try:
        parsed = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise HLT17SRCQ1REC1IdentityError(f"{label} is not JSON") from error
    return _mapping(parsed, label)


def _require_bool(value: object, expected: bool, label: str) -> bool:
    if type(value) is not bool or value is not expected:
        raise HLT17SRCQ1REC1IdentityError(f"{label} differs")
    return value


def _plan_mapping(plan: object) -> dict[str, object] | None:
    if plan is None:
        return None
    if isinstance(plan, Mapping):
        return dict(plan)
    if hasattr(plan, "requested_cap"):
        return encode_tdg7_plan(plan)
    raise HLT17SRCQ1REC1Error("attempted or next plan cannot be encoded")


def implementation_identity(repository_root: Path | None = None) -> dict[str, object]:
    identity = hlt17_implementation_identity()
    record: dict[str, object] = {
        "artifact_id": ARTIFACT_ID,
        "schema": SCHEMA,
        "path": IMPLEMENTATION_RELATIVE,
        "predecessor_artifact_id": PREDECESSOR_ARTIFACT_ID,
        "protocol": PROTOCOL,
        "campaign_id": CAMPAIGN_ID,
        "c1r1_actual": dict(identity),
        "imp1_reference_wire": {
            "reference_wire_artifact_id": identity["reference_wire_artifact_id"],
            "reference_wire_evaluator_id": identity["reference_wire_evaluator_id"],
        },
        "c1r1_is_not_imp1_execution": True,
        "campaign_execution_authorized": False,
        "rk_members_remeasured": False,
    }
    if identity["implementation_id"] == identity["reference_wire_artifact_id"]:
        raise HLT17SRCQ1REC1IdentityError(
            "C1R1 actual identity collapsed into the IMP1 wire"
        )
    if repository_root is not None:
        raw = read_regular_file(Path(repository_root), IMPLEMENTATION_RELATIVE)
        record["sha256"] = _sha256_hex(raw)
    return record


def plan_status() -> dict[str, object]:
    """Default CLI/status record. Performs no I/O, source call, or write."""

    identity = hlt17_implementation_identity()
    ceilings = recommended_resource_ceilings()
    return {
        "artifact_id": ARTIFACT_ID,
        "schema": SCHEMA,
        "role": "implementation_A_no_write_recovery_core",
        "physical_source_qualification_executed": False,
        "rk_members_remeasured": False,
        "campaign_execution_authorized": False,
        "store_publication_authorized": False,
        "endpoint_adoption_authorized": False,
        "authority_delta_validation": AUTHORITY_SEAM,
        "protocol": PROTOCOL,
        "campaign_id": CAMPAIGN_ID,
        "event_target_hex": EVENT_TARGET.hex(),
        "origin_time_hex": ORIGIN_TIME_HEX,
        "resume_member_keys": list(RESUME_MEMBER_KEYS),
        "predecessor_member_keys": list(PREDECESSOR_MEMBER_KEYS),
        "member_keys": list(MEMBER_KEYS),
        "frozen_prospective_requested_cap_hex_by_member": dict(
            FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER
        ),
        "prospective_requested_cap_formula": (
            "cfl_maximum*grid_spacing/inherited_previous_speed_upper"
        ),
        "attempt2": {
            "artifact_id": ATTEMPT2_ARTIFACT_ID,
            "classification": ATTEMPT2_CLASSIFICATION,
            "compact_relative": ATTEMPT2_COMPACT_RELATIVE,
            "compact_identities": {
                "sealed_historical_blob": {
                    "commit": ATTEMPT2_SEALED_COMPACT_COMMIT,
                    "sha256": ATTEMPT2_SEALED_COMPACT_SHA256,
                    "derivation_document": None,
                },
                "metadata_linked_blob": {
                    "commit": ATTEMPT2_METADATA_LINKED_COMPACT_COMMIT,
                    "sha256": ATTEMPT2_METADATA_LINKED_COMPACT_SHA256,
                    "derivation_document": ATTEMPT2_DERIVATION_DOCUMENT,
                },
            },
            "sealed_compact_sha256": ATTEMPT2_SEALED_COMPACT_SHA256,
            "metadata_linked_compact_sha256": ATTEMPT2_METADATA_LINKED_COMPACT_SHA256,
            "compact_sha256": ATTEMPT2_SEALED_COMPACT_SHA256,
            "raw_leaf_relative": ATTEMPT2_RAW_LEAF_RELATIVE,
            "raw_leaf_sha256": ATTEMPT2_RAW_LEAF_SHA256,
            "raw_directory_payload_sha256": ATTEMPT2_RAW_DIRECTORY_PAYLOAD_SHA256,
            "raw_leaf_bytes": ATTEMPT2_RAW_LEAF_BYTES,
            "retry_authorized": False,
        },
        "static_input_paths": list(STATIC_INPUT_PATHS),
        "static_input_sha256": dict(STATIC_INPUT_SHA256),
        "c1r1_actual": dict(identity),
        "imp1_reference_wire": {
            "reference_wire_artifact_id": identity["reference_wire_artifact_id"],
            "reference_wire_evaluator_id": identity["reference_wire_evaluator_id"],
        },
        "resource_ceilings": ceilings.as_mapping(),
        "source_retry_cap": HLT17_SOURCE_RETRY_CAP,
        "cfl_retry_cap": HLT17_CFL_RETRY_CAP,
        "one_real_run": {
            "requires": ["exact --authority-commit", "absent --output-directory"],
            "forbids": [
                "boolean skip/force authority flag",
                "existing namespace",
                "cap retune after a result",
                "RK remeasurement",
                "informal retry",
                "endpoint serialization",
                "campaign store publication",
            ],
            "authority_module": AUTHORITY_VALIDATOR_MODULE,
            "authority_symbol": AUTHORITY_VALIDATOR_NAME,
        },
    }


def load_coordinator_authority_delta_validator() -> Callable[..., Mapping[str, object]]:
    try:
        module = importlib.import_module(AUTHORITY_VALIDATOR_MODULE)
    except ImportError as error:
        raise HLT17SRCQ1REC1AuthoritySeamError(AUTHORITY_SEAM) from error
    validator = getattr(module, AUTHORITY_VALIDATOR_NAME, None)
    if not callable(validator):
        raise HLT17SRCQ1REC1AuthoritySeamError(AUTHORITY_SEAM)
    return validator


def require_authority_delta(
    *,
    repository_root: Path,
    authority_commit: str,
    implementation_sha256: str,
    validator: Callable[..., object] | None = None,
) -> dict[str, object]:
    commit = _require_commit(authority_commit, "authority_commit")
    digest = _require_sha(implementation_sha256, "implementation_sha256")
    loaded = validator
    if loaded is None:
        loaded = load_coordinator_authority_delta_validator()
    if not callable(loaded):
        raise HLT17SRCQ1REC1AuthoritySeamError(AUTHORITY_SEAM)
    receipt = loaded(
        repository_root=Path(repository_root),
        authority_commit=commit,
        implementation_sha256=digest,
    )
    if type(receipt) is bool:
        raise HLT17SRCQ1REC1AuthoritySeamError(
            "authority/delta validation must not be a boolean flag"
        )
    mapping = _mapping(receipt, "authority delta receipt")
    if mapping.get("authority_commit") != commit:
        raise HLT17SRCQ1REC1AuthoritySeamError(
            "authority receipt authority_commit differs"
        )
    if mapping.get("implementation_sha256") != digest:
        raise HLT17SRCQ1REC1AuthoritySeamError(
            "authority receipt implementation_sha256 differs"
        )
    closure = mapping.get("source_closure_sha256")
    if type(closure) is not str or _SHA.fullmatch(closure) is None:
        raise HLT17SRCQ1REC1AuthoritySeamError(
            "authority receipt omits exact source_closure_sha256"
        )
    mapping["boolean_flag_accepted"] = False
    return mapping


def require_frozen_srcq1_caps(value: object) -> dict[str, str]:
    mapping = _mapping(value, "frozen SRCQ1 requested caps")
    if tuple(mapping) != MEMBER_KEYS:
        raise HLT17SRCQ1REC1CapError("frozen requested-cap order/inventory differs")
    result: dict[str, str] = {}
    for key in MEMBER_KEYS:
        text = mapping[key]
        expected = FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER[key]
        if type(text) is not str or text.lower() != text or text != expected:
            raise HLT17SRCQ1REC1CapError(
                f"{key} requested cap is not the frozen SRCQ1 hex"
            )
        result[key] = text
    return result


def authenticate_attempt2_predecessor(
    *,
    compact_bytes: bytes,
    raw_bytes: bytes | None = None,
    require_raw: bool = True,
) -> dict[str, object]:
    """Authenticate Attempt2 compact/raw hashes and state-safe flags."""

    compact_digest = _sha256_hex(compact_bytes)
    compact = _json_mapping(compact_bytes, "Attempt2 compact")
    if compact_digest == ATTEMPT2_SEALED_COMPACT_SHA256:
        if "derivation_document" in compact:
            raise HLT17SRCQ1REC1IdentityError(
                "sealed Attempt2 compact must not carry derivation_document"
            )
        compact_identity = "sealed_historical_blob"
        derivation_document = None
    elif compact_digest == ATTEMPT2_METADATA_LINKED_COMPACT_SHA256:
        if compact.get("derivation_document") != ATTEMPT2_DERIVATION_DOCUMENT:
            raise HLT17SRCQ1REC1IdentityError(
                "metadata-linked Attempt2 compact derivation_document differs"
            )
        compact_identity = "metadata_linked_blob"
        derivation_document = ATTEMPT2_DERIVATION_DOCUMENT
    else:
        raise HLT17SRCQ1REC1IdentityError("Attempt2 compact SHA-256 differs")
    if compact.get("artifact_id") != ATTEMPT2_ARTIFACT_ID:
        raise HLT17SRCQ1REC1IdentityError("Attempt2 compact artifact differs")
    if compact.get("classification") != ATTEMPT2_CLASSIFICATION:
        raise HLT17SRCQ1REC1IdentityError("Attempt2 classification differs")
    if compact.get("authority_commit") != ATTEMPT2_AUTHORITY_COMMIT:
        raise HLT17SRCQ1REC1IdentityError("Attempt2 authority commit differs")
    if compact.get("implementation_commit") != ATTEMPT2_IMPLEMENTATION_COMMIT:
        raise HLT17SRCQ1REC1IdentityError("Attempt2 implementation commit differs")
    if compact.get("implementation_sha256") != ATTEMPT2_IMPLEMENTATION_SHA256:
        raise HLT17SRCQ1REC1IdentityError("Attempt2 implementation SHA-256 differs")
    if compact.get("origin_capture_sha256") != ORIGIN_CAPTURE_SHA256:
        raise HLT17SRCQ1REC1IdentityError("Attempt2 origin capture differs")
    _require_bool(compact.get("retry_authorized"), False, "Attempt2 retry_authorized")
    if compact.get("completed_member_count") != 4:
        raise HLT17SRCQ1REC1IdentityError("Attempt2 completed-member count differs")
    state = _mapping(compact.get("state_boundary"), "Attempt2 state_boundary")
    _require_bool(state.get("accepted_state_advanced"), False, "accepted_state_advanced")
    _require_bool(state.get("campaign_store_written"), False, "campaign_store_written")
    _require_bool(
        state.get("endpoint_adopted_or_serialized"),
        False,
        "endpoint_adopted_or_serialized",
    )
    _require_bool(state.get("pro20_namespace_absent"), True, "pro20_namespace_absent")
    raw_terminal = _mapping(compact.get("raw_terminal"), "Attempt2 raw_terminal")
    if raw_terminal.get("leaf_path") != ATTEMPT2_RAW_LEAF_RELATIVE:
        raise HLT17SRCQ1REC1IdentityError("Attempt2 raw leaf path differs")
    if raw_terminal.get("leaf_sha256") != ATTEMPT2_RAW_LEAF_SHA256:
        raise HLT17SRCQ1REC1IdentityError("Attempt2 recorded raw leaf SHA-256 differs")
    if raw_terminal.get("directory_payload_sha256") != ATTEMPT2_RAW_DIRECTORY_PAYLOAD_SHA256:
        raise HLT17SRCQ1REC1IdentityError(
            "Attempt2 recorded directory payload SHA-256 differs"
        )
    if raw_terminal.get("leaf_bytes") != ATTEMPT2_RAW_LEAF_BYTES:
        raise HLT17SRCQ1REC1IdentityError("Attempt2 recorded raw leaf bytes differ")
    outcomes = compact.get("member_outcomes")
    if type(outcomes) is not list or len(outcomes) != 4:
        raise HLT17SRCQ1REC1IdentityError("Attempt2 member outcomes differ")
    by_key = {}
    for item in outcomes:
        row = _mapping(item, "Attempt2 member outcome")
        key = row.get("member_key")
        if type(key) is not str or key in by_key:
            raise HLT17SRCQ1REC1IdentityError("Attempt2 member outcome key differs")
        by_key[key] = row
    if tuple(by_key) != PREDECESSOR_MEMBER_KEYS + ("SSPRK3-4097",):
        raise HLT17SRCQ1REC1IdentityError("Attempt2 member outcome order differs")
    for key in PREDECESSOR_MEMBER_KEYS:
        row = by_key[key]
        if row.get("outcome") != "prepared":
            raise HLT17SRCQ1REC1IdentityError(f"Attempt2 {key} was not prepared")
        _require_bool(
            row.get("source_evaluation_recorded"),
            True,
            f"{key} source_evaluation_recorded",
        )
    ssp = by_key["SSPRK3-4097"]
    if ssp.get("outcome") != "runner_reported_cfl":
        raise HLT17SRCQ1REC1IdentityError("Attempt2 SSPRK3-4097 outcome differs")
    _require_bool(
        ssp.get("source_evaluation_recorded"),
        True,
        "SSPRK3-4097 source_evaluation_recorded",
    )
    _require_bool(
        ssp.get("cfl_ratio_evidence_recorded"),
        False,
        "SSPRK3-4097 cfl_ratio_evidence_recorded",
    )
    raw_mapping: dict[str, Any] | None = None
    raw_digest = None
    if require_raw or raw_bytes is not None:
        if raw_bytes is None:
            raise HLT17SRCQ1REC1IdentityError("Attempt2 raw terminal bytes are absent")
        if len(raw_bytes) != ATTEMPT2_RAW_LEAF_BYTES:
            raise HLT17SRCQ1REC1IdentityError("Attempt2 raw leaf byte length differs")
        raw_digest = _sha256_hex(raw_bytes)
        if raw_digest != ATTEMPT2_RAW_LEAF_SHA256:
            raise HLT17SRCQ1REC1IdentityError("Attempt2 raw leaf SHA-256 differs")
        raw_mapping = _json_mapping(raw_bytes, "Attempt2 raw terminal")
        _require_bool(
            raw_mapping.get("accepted_state_advanced"),
            False,
            "raw accepted_state_advanced",
        )
        _require_bool(raw_mapping.get("endpoint_adopted"), False, "raw endpoint_adopted")
        _require_bool(raw_mapping.get("store_published"), False, "raw store_published")
        _require_bool(
            raw_mapping.get("campaign_execution_authorized"),
            False,
            "raw campaign_execution_authorized",
        )
        if raw_mapping.get("halted_outcome") != "cfl":
            raise HLT17SRCQ1REC1IdentityError("Attempt2 raw halted_outcome differs")
        raw_members = raw_mapping.get("members")
        if type(raw_members) is not list or len(raw_members) != 4:
            raise HLT17SRCQ1REC1IdentityError("Attempt2 raw member inventory differs")
        ssp_raw = _mapping(raw_members[3], "Attempt2 raw SSPRK3-4097")
        if ssp_raw.get("member_key") != "SSPRK3-4097":
            raise HLT17SRCQ1REC1IdentityError("Attempt2 raw resume member differs")
        if ssp_raw.get("prepare_disposition") != CFL_REQUIRED:
            raise HLT17SRCQ1REC1IdentityError(
                "Attempt2 raw SSPRK3-4097 disposition differs"
            )
        if ssp_raw.get("c1r1") is not None:
            raise HLT17SRCQ1REC1IdentityError("Attempt2 raw SSPRK3-4097 carried C1R1")
        encoded = json.dumps(ssp_raw, sort_keys=True)
        if "observed_ratio_hex" in encoded or "maximum_ratio_hex" in encoded:
            raise HLT17SRCQ1REC1IdentityError(
                "Attempt2 raw unexpectedly contains CFL ratio evidence"
            )
    record = {
        "artifact_id": ATTEMPT2_ARTIFACT_ID,
        "classification": ATTEMPT2_CLASSIFICATION,
        "compact_identity": compact_identity,
        "compact_sha256": compact_digest,
        "sealed_compact_sha256": ATTEMPT2_SEALED_COMPACT_SHA256,
        "metadata_linked_compact_sha256": ATTEMPT2_METADATA_LINKED_COMPACT_SHA256,
        "derivation_document": derivation_document,
        "raw_leaf_sha256": raw_digest if raw_digest is not None else ATTEMPT2_RAW_LEAF_SHA256,
        "raw_directory_payload_sha256": ATTEMPT2_RAW_DIRECTORY_PAYLOAD_SHA256,
        "raw_leaf_bytes": ATTEMPT2_RAW_LEAF_BYTES,
        "retry_authorized": False,
        "cfl_terminal_independently_reproduced": False,
        "runner_label_sufficient": False,
        "accepted_state_advanced": False,
        "campaign_store_written": False,
        "endpoint_adopted_or_serialized": False,
        "pro20_namespace_absent": True,
        "compact": compact,
        "raw": raw_mapping,
        "member_outcomes": [dict(by_key[key]) for key in by_key],
    }
    refuse_endpoint_payload(record)
    return record


def predecessor_rk_evidence_references(
    attempt2: Mapping[str, object],
) -> list[dict[str, object]]:
    compact = _mapping(attempt2, "Attempt2 predecessor")
    outcomes = compact.get("member_outcomes")
    if type(outcomes) is not list:
        outcomes = _mapping(compact.get("compact"), "Attempt2 compact").get(
            "member_outcomes", []
        )
    by_key: dict[str, dict[str, object]] = {}
    for item in outcomes:
        row = _mapping(item, "predecessor outcome")
        key = row.get("member_key")
        if type(key) is str:
            by_key[key] = row
    raw_members = []
    raw = compact.get("raw")
    if isinstance(raw, Mapping):
        members = raw.get("members")
        if type(members) is list:
            raw_members = members
    raw_by_key = {}
    for item in raw_members:
        if isinstance(item, Mapping) and type(item.get("member_key")) is str:
            raw_by_key[str(item["member_key"])] = dict(item)
    records: list[dict[str, object]] = []
    for key in PREDECESSOR_MEMBER_KEYS:
        row = by_key.get(key)
        if row is None:
            raise HLT17SRCQ1REC1IdentityError(f"predecessor evidence omitted {key}")
        raw_row = raw_by_key.get(key, {})
        c1r1 = raw_row.get("c1r1") if isinstance(raw_row, Mapping) else None
        receipt = None
        disposition = raw_row.get("prepare_disposition") if raw_row else None
        if isinstance(c1r1, Mapping):
            receipt = c1r1.get("receipt_sha256")
        record = {
            "member_key": key,
            "role": "predecessor_evidence",
            "remeasured": False,
            "source_evaluated_by_rec1": False,
            "outcome": row.get("outcome"),
            "source_evaluation_recorded": row.get("source_evaluation_recorded"),
            "attempt2_prepare_disposition": disposition,
            "attempt2_c1r1_receipt_sha256": receipt,
            "attempt2_compact_sha256": compact.get("compact_sha256", ATTEMPT2_SEALED_COMPACT_SHA256),
            "attempt2_sealed_compact_sha256": ATTEMPT2_SEALED_COMPACT_SHA256,
            "attempt2_metadata_linked_compact_sha256": ATTEMPT2_METADATA_LINKED_COMPACT_SHA256,
            "attempt2_raw_leaf_sha256": ATTEMPT2_RAW_LEAF_SHA256,
            "accepted_state_advanced": False,
            "endpoint_adopted": False,
        }
        refuse_endpoint_payload(record)
        records.append(record)
    return records


def _require_positive_finite_hex(value: object, label: str) -> float:
    if type(value) is not str or value.lower() != value:
        raise HLT17SRCQ1REC1Error(f"{label} is not canonical lowercase binary64 hex")
    try:
        number = float.fromhex(value)
    except ValueError as error:
        raise HLT17SRCQ1REC1Error(f"{label} is malformed") from error
    if not math.isfinite(number) or number <= 0.0:
        raise HLT17SRCQ1REC1Error(f"{label} is not a finite positive binary64")
    if number.hex() != value:
        raise HLT17SRCQ1REC1Error(f"{label} is not a round-trip binary64")
    return number


def require_independent_cfl_ratio(evidence: object) -> dict[str, object]:
    """Require kind cfl_retry and observed_ratio>maximum_ratio. Labels never prove CFL."""

    mapping = _mapping(evidence, "CFL evidence")
    observed_hex = mapping.get("observed_ratio_hex")
    maximum_hex = mapping.get("maximum_ratio_hex")
    if type(observed_hex) is not str or type(maximum_hex) is not str:
        raise HLT17SRCQ1REC1InconclusiveError(
            "CFL stop omitted observed/max ratio evidence"
        )
    if mapping.get("kind") != "cfl_retry":
        raise HLT17SRCQ1REC1Error("CFL evidence kind is not cfl_retry")
    observed = _require_positive_finite_hex(observed_hex, "observed_ratio_hex")
    maximum = _require_positive_finite_hex(maximum_hex, "maximum_ratio_hex")
    if not (observed > maximum):
        raise HLT17SRCQ1REC1Error(
            "CFL stop does not independently satisfy observed_ratio>maximum_ratio"
        )
    return {
        "kind": "cfl_retry",
        "observed_ratio_hex": observed.hex(),
        "maximum_ratio_hex": maximum.hex(),
        "observed_exceeds_maximum": True,
        "runner_label_sufficient": False,
    }


def require_source_overlay_evidence(evidence: object) -> dict[str, object]:
    mapping = _mapping(evidence, "source overlay evidence")
    if mapping.get("kind") != "source_retry":
        raise HLT17SRCQ1REC1InconclusiveError("source overlay omitted source_retry kind")
    if not isinstance(mapping.get("source_retry"), Mapping):
        raise HLT17SRCQ1REC1InconclusiveError(
            "source overlay omitted wired source-retry evidence"
        )
    return dict(mapping)


def retain_bridge_result(result: object) -> dict[str, object]:
    """Retain complete HLT17BridgeResult evidence without endpoints."""

    if getattr(result, "accepted_state_advanced", False) is True:
        raise HLT17SRCQ1REC1Error("bridge result advanced accepted state")
    if getattr(result, "committed", None) is not None:
        raise HLT17SRCQ1REC1Error("bridge result serialized a committed endpoint")
    disposition = getattr(result, "disposition", None)
    if type(disposition) is not str:
        raise HLT17SRCQ1REC1Error("bridge result omitted disposition")
    evidence = getattr(result, "evidence", {}) or {}
    evidence_map = _mapping(evidence, "bridge evidence")
    if not hasattr(result, "sink_writes"):
        raise HLT17SRCQ1REC1Error("bridge result omitted sink_writes")
    sink_writes = result.sink_writes
    if type(sink_writes) is not int or sink_writes != 0:
        raise HLT17SRCQ1REC1Error("bridge sink_writes must be the integer zero")
    cursor = getattr(result, "cursor", None)
    overlay = getattr(cursor, "overlay", None) if cursor is not None else None
    overlay_map = None
    if overlay is not None and hasattr(overlay, "as_mapping"):
        overlay_map = overlay.as_mapping()
    elif isinstance(overlay, Mapping):
        overlay_map = dict(overlay)
    record = {
        "disposition": disposition,
        "accepted_state_advanced": False,
        "sink_writes": sink_writes,
        "evidence": dict(evidence_map),
        "attempted_plan": _plan_mapping(getattr(result, "attempted_plan", None)),
        "next_plan": _plan_mapping(getattr(result, "next_plan", None)),
        "prepared_present": getattr(result, "prepared", None) is not None,
        "committed_present": False,
        "replay_present": getattr(result, "replay", None) is not None,
        "cursor": {
            "source_current": getattr(cursor, "source_current", None),
            "source_total": getattr(cursor, "source_total", None),
            "cfl_current": getattr(cursor, "cfl_current", None),
            "cfl_total": getattr(cursor, "cfl_total", None),
            "overlay": overlay_map,
        },
    }
    for key in FORBIDDEN_ENDPOINT_KEYS:
        if key in record or key in evidence_map:
            raise HLT17SRCQ1REC1Error("bridge result serializes a forbidden endpoint")
    refuse_endpoint_payload(record)
    canonical_json_bytes(record)
    return record


def isolate_working_checkpoint(checkpoint: object) -> object:
    """Copy the in-memory checkpoint. Overlay absorption must not touch live."""

    isolate = getattr(checkpoint, "isolate", None)
    if callable(isolate):
        working = isolate()
        if working is checkpoint:
            raise HLT17SRCQ1REC1Error("isolation returned the live checkpoint")
        return working
    if isinstance(checkpoint, HLT17InMemoryCheckpoint):
        isolated_member = isolate_hlt17_runtime_member(checkpoint.member)
        working = HLT17InMemoryCheckpoint(
            isolated_member,
            checkpoint.cursor,
            checkpoint.accepted_encoding,
        )
        if working is checkpoint or working.member is checkpoint.member:
            raise HLT17SRCQ1REC1Error("HLT17 isolation failed to copy the member")
        return working
    raise HLT17SRCQ1REC1Error("checkpoint cannot be isolated in memory")


def classify_prepare_outcome(
    result: object, prepared_record: Mapping[str, object] | None
) -> str:
    disposition = getattr(result, "disposition", None)
    if disposition in SOURCE_DISPOSITIONS:
        return "source"
    if disposition in CFL_DISPOSITIONS:
        return "cfl"
    if disposition in TEMPORAL_DISPOSITIONS:
        return "temporal"
    if disposition == "invalid_premise":
        return "invalid"
    if disposition == "accepted_fine":
        raise HLT17SRCQ1REC1Error("qualification adopted an endpoint")
    if disposition not in PREPARED_DISPOSITIONS:
        return "inconclusive"
    if not prepared_record:
        return "inconclusive"
    decisions = prepared_record.get("eighteen_channel_decisions")
    if type(decisions) is not list:
        return "inconclusive"
    classifications = [item.get("classification") for item in decisions]
    if any(item == "order_inconclusive" for item in classifications):
        return "inconclusive"
    if prepared_record.get("admission_passed") is True:
        return "prepared"
    return "temporal"


def _working_counters(working: object, result: object) -> dict[str, object]:
    member = getattr(working, "member", None)
    cursor = getattr(result, "cursor", None)
    if cursor is None:
        cursor = getattr(working, "cursor", None)
    return {
        "working_source_retry_count": getattr(
            member, "source_retry_count", getattr(working, "source_retry_count", None)
        ),
        "working_CFL_retry_count": getattr(
            member, "CFL_retry_count", getattr(working, "CFL_retry_count", None)
        ),
        "cursor_source_current": getattr(cursor, "source_current", None),
        "cursor_source_total": getattr(cursor, "source_total", None),
        "cursor_cfl_current": getattr(cursor, "cfl_current", None),
        "cursor_cfl_total": getattr(cursor, "cfl_total", None),
    }


def _overlay_exhausted(result: object) -> bool:
    disposition = getattr(result, "disposition", None)
    if disposition in {CFL_EXHAUSTED, SOURCE_EXHAUSTED}:
        return True
    cursor = getattr(result, "cursor", None)
    overlay = getattr(cursor, "overlay", None) if cursor is not None else None
    if overlay is not None and getattr(overlay, "exhausted", False) is True:
        return True
    if isinstance(overlay, Mapping) and overlay.get("exhausted") is True:
        return True
    return getattr(result, "next_plan", None) is None


def _absorb_overlay(working: object, result: object) -> object:
    absorb = getattr(working, "absorb_nonfine", None)
    if not callable(absorb):
        raise HLT17SRCQ1REC1Error("isolated checkpoint cannot absorb a source/CFL overlay")
    absorbed = absorb(result)
    if getattr(absorbed, "accepted_state_advanced", False) is True:
        raise HLT17SRCQ1REC1Error("overlay absorption advanced accepted state")
    return absorbed


def _captured_member(origin: object, member_key: str):
    members = getattr(origin, "members", ())
    matches = [item for item in members if getattr(item, "member_key", None) == member_key]
    if len(matches) != 1:
        raise HLT17SRCQ1REC1IdentityError(f"captured origin member differs: {member_key}")
    return matches[0]


def _require_expected_identity(
    actual: Mapping[str, object],
    expected: Mapping[str, object] | None,
    label: str,
) -> None:
    if expected is None:
        return
    wanted = _mapping(expected, label)
    for key, value in wanted.items():
        if actual.get(key) != value:
            raise HLT17SRCQ1REC1IdentityError(f"foreign {label}: {key}")


def _shadow_record(prepared: object) -> dict[str, object] | None:
    if prepared is None:
        return None
    return c1r1_shadow_record(prepared)


def qualify_member_with_overlays(
    *,
    member_key: str,
    view: object,
    captured: object,
    initial_cap: float,
    evaluate: Callable[[object], Mapping[str, object]],
    now: Callable[[], float],
    rss: Callable[[], int],
    ceilings: ResourceCeilings,
    started: float,
    expected_config: str | None,
) -> dict[str, object]:
    """Evaluate source once, then follow isolated overlay.next_plan."""

    member = view.member
    binding = view.source_binding
    owned = check_method_owned_identity(member_key, member)
    check_inherited_counters(member_key, member, captured)
    configuration = _require_sha(
        getattr(binding, "configuration_sha256", None),
        f"{member_key} source configuration",
    )
    if expected_config is not None and expected_config != configuration:
        raise HLT17SRCQ1REC1IdentityError(f"foreign source configuration: {member_key}")
    if configuration != getattr(view, "source_configuration_sha256", configuration):
        raise HLT17SRCQ1REC1IdentityError(
            f"{member_key} source configuration identity differs"
        )
    if initial_cap.hex() != FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER[member_key]:
        raise HLT17SRCQ1REC1CapError(f"{member_key} initial cap is not frozen SRCQ1 hex")
    before = member_fingerprints(member, binding)
    member_started = now()
    rss_before = int(rss())
    check_srcq1_resources(
        member_key=member_key,
        ceilings=ceilings,
        wall_seconds=0.0,
        total_wall_seconds=now() - started,
        rss=rss_before,
    )
    source_record = dict(evaluate(view))
    live_checkpoint = view.checkpoint
    if hasattr(live_checkpoint, "agree"):
        live_checkpoint.agree()
    working = isolate_working_checkpoint(live_checkpoint)
    if hasattr(working, "agree"):
        working.agree()
    attempts: list[dict[str, object]] = []
    expected_plan: dict[str, object] | None = None
    halted: str | None = None
    c1r1_record = None
    outcome = "inconclusive"
    last_disposition = None
    prepare_count = 0
    overlay_count = 0
    source_owner_count = 0
    cfl_owner_count = 0
    stop_reason: str | None = None
    for prepare_index in range(MAX_PREPARE_ATTEMPTS):
        prepared_result = working.prepare(requested_cap=initial_cap)
        prepare_count += 1
        if getattr(prepared_result, "accepted_state_advanced", False) is True:
            raise HLT17SRCQ1REC1Error("C1R1 shadow prepare advanced accepted state")
        if getattr(prepared_result, "committed", None) is not None:
            raise HLT17SRCQ1REC1Error("C1R1 shadow prepare carried a committed runtime")
        retained = retain_bridge_result(prepared_result)
        last_disposition = retained["disposition"]
        if expected_plan is not None:
            if retained["attempted_plan"] != expected_plan:
                raise HLT17SRCQ1REC1Error(
                    f"{member_key} prepare did not follow overlay.next_plan"
                )
        elif prepare_index == 0:
            attempted = retained["attempted_plan"]
            if (
                isinstance(attempted, Mapping)
                and attempted.get("requested_cap_hex")
                != FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER[member_key]
            ):
                raise HLT17SRCQ1REC1CapError(
                    f"{member_key} first prepare is not the frozen initial cap"
                )
        prepared = getattr(prepared_result, "prepared", None)
        member_c1r1 = _shadow_record(prepared) if prepared is not None else None
        classified = classify_prepare_outcome(prepared_result, member_c1r1)
        attempt_record = {
            "prepare_index": prepare_index,
            "requested_initial_cap_hex": initial_cap.hex(),
            "bridge": retained,
            "outcome": classified,
            "counters": _working_counters(working, prepared_result),
            "overlay_absorbed": False,
            "source_owner_count": source_owner_count,
            "cfl_owner_count": cfl_owner_count,
        }
        try:
            if last_disposition in CFL_DISPOSITIONS:
                retained["independent_cfl"] = require_independent_cfl_ratio(
                    retained["evidence"]
                )
            if last_disposition in SOURCE_DISPOSITIONS:
                retained["independent_source"] = require_source_overlay_evidence(
                    retained["evidence"]
                )
        except HLT17SRCQ1Error as error:
            attempt_record["outcome"] = error.outcome
            attempt_record["reason"] = str(error)
            refuse_endpoint_payload(attempt_record)
            attempts.append(attempt_record)
            outcome = error.outcome
            halted = error.outcome
            stop_reason = str(error)
            break
        owner_stop: str | None = None
        if last_disposition in CFL_DISPOSITIONS:
            cfl_owner_count += 1
            attempt_record["cfl_owner_count"] = cfl_owner_count
            if cfl_owner_count > HLT17_CFL_RETRY_CAP + 1:
                owner_stop = "CFL owner count exceeded the exhaustion edge"
            elif (
                cfl_owner_count == HLT17_CFL_RETRY_CAP + 1
                and last_disposition != CFL_EXHAUSTED
                and not _overlay_exhausted(prepared_result)
            ):
                owner_stop = "CFL owner count 33 is not exhausted"
        if last_disposition in SOURCE_DISPOSITIONS:
            source_owner_count += 1
            attempt_record["source_owner_count"] = source_owner_count
            if source_owner_count > HLT17_SOURCE_RETRY_CAP + 1:
                owner_stop = "source owner count exceeded the exhaustion edge"
            elif (
                source_owner_count == HLT17_SOURCE_RETRY_CAP + 1
                and last_disposition != SOURCE_EXHAUSTED
                and not _overlay_exhausted(prepared_result)
            ):
                owner_stop = "source owner count 33 is not exhausted"
        if owner_stop is not None:
            attempt_record["outcome"] = "invalid"
            attempt_record["reason"] = owner_stop
            refuse_endpoint_payload(attempt_record)
            attempts.append(attempt_record)
            outcome = "invalid"
            halted = "invalid"
            stop_reason = owner_stop
            break
        refuse_endpoint_payload(attempt_record)
        attempts.append(attempt_record)
        if classified == "prepared":
            c1r1_record = member_c1r1
            outcome = "prepared"
            halted = None
            break
        if last_disposition in OVERLAY_CONTINUE or last_disposition in {
            CFL_EXHAUSTED,
            SOURCE_EXHAUSTED,
        }:
            _absorb_overlay(working, prepared_result)
            attempt_record["overlay_absorbed"] = True
            overlay_count += 1
            attempt_record["counters"] = _working_counters(working, prepared_result)
            if last_disposition in {CFL_EXHAUSTED, SOURCE_EXHAUSTED} or _overlay_exhausted(
                prepared_result
            ):
                outcome = classified
                halted = classified
                break
            next_plan = retained["next_plan"]
            if not isinstance(next_plan, Mapping):
                attempt_record["outcome"] = "invalid"
                attempt_record["reason"] = f"{member_key} overlay omitted next_plan"
                outcome = "invalid"
                halted = "invalid"
                stop_reason = attempt_record["reason"]
                break
            expected_plan = next_plan
            continue
        outcome = classified
        halted = classified
        break
    else:
        outcome = "invalid"
        halted = "invalid"
        stop_reason = "prepare loop exceeded independent source/CFL exhaustion edges"
        attempts.append(
            {
                "prepare_index": prepare_count,
                "outcome": "invalid",
                "reason": stop_reason,
                "overlay_absorbed": False,
                "source_owner_count": source_owner_count,
                "cfl_owner_count": cfl_owner_count,
            }
        )
    wall = now() - member_started
    rss_after = int(rss())
    try:
        check_srcq1_resources(
            member_key=member_key,
            ceilings=ceilings,
            wall_seconds=wall,
            total_wall_seconds=now() - started,
            rss=rss_after,
        )
    except HLT17SRCQ1ResourceError as error:
        outcome = "resource"
        halted = "resource"
        stop_reason = str(error)
    after = member_fingerprints(member, binding)
    if after != before:
        outcome = "source"
        halted = "source"
        stop_reason = (
            f"{member_key} live fingerprints changed after isolated overlay work"
        )
    record = {
        "member_key": member_key,
        "role": "newly_owned",
        "remeasured": True,
        "source_evaluated_by_rec1": True,
        "outcome": outcome,
        "reason": stop_reason,
        "requested_cap_hex": initial_cap.hex(),
        "owned": owned,
        "source_evaluation": source_record,
        "c1r1": c1r1_record,
        "prepare_disposition": last_disposition,
        "attempts": attempts,
        "prepare_count": prepare_count,
        "overlay_count": overlay_count,
        "source_owner_count": source_owner_count,
        "cfl_owner_count": cfl_owner_count,
        "fingerprints_before": before,
        "fingerprints_after": after,
        "inherited_counters_preserved_on_live_member": True,
        "wall_seconds": wall,
        "rss_bytes_before": rss_before,
        "rss_bytes_after": rss_after,
        "accepted_state_advanced": False,
        "endpoint_adopted": False,
        "halted_outcome": halted,
    }
    refuse_endpoint_payload(record)
    return record


def qualify_srcq1_rec1(
    origin: object,
    *,
    static_input_bytes: Mapping[str, bytes],
    source_closure_sha256: str,
    environment: Mapping[str, str],
    attempt2_compact_bytes: bytes,
    attempt2_raw_bytes: bytes,
    resource_ceilings: ResourceCeilings | None = None,
    runtime_origin: object | None = None,
    runtime_origin_builder: Callable[..., object] | None = None,
    source_evaluator: Callable[..., Mapping[str, object]] | None = None,
    clock: Callable[[], float] | None = None,
    rss_bytes: Callable[[], int] | None = None,
    expected_origin_sha256: str | None = None,
    expected_environment: Mapping[str, str] | None = None,
    expected_implementation: Mapping[str, str] | None = None,
    expected_source_configurations: Mapping[str, str] | None = None,
    prospective_requested_cap_hex_by_member: Mapping[str, str] | None = None,
) -> dict[str, object]:
    """Recover SRCQ1 from SSPRK3-4097. No write, adopt, RK remeasurement, or store."""

    attempt2 = authenticate_attempt2_predecessor(
        compact_bytes=attempt2_compact_bytes,
        raw_bytes=attempt2_raw_bytes,
        require_raw=True,
    )
    predecessor = predecessor_rk_evidence_references(attempt2)
    if isinstance(origin, Pro20OriginCapture):
        try:
            revalidate_pro20_origin_capture(origin)
        except Pro20OriginError as error:
            raise HLT17SRCQ1REC1IdentityError(
                "origin capture is an invalid premise"
            ) from error
    elif not hasattr(origin, "capture_sha256") or not hasattr(origin, "members"):
        raise HLT17SRCQ1REC1IdentityError("origin capture type differs")
    expected_origin = expected_origin_sha256 or ORIGIN_CAPTURE_SHA256
    if origin.capture_sha256 != _require_sha(expected_origin, "expected origin"):
        raise HLT17SRCQ1REC1IdentityError("foreign origin identity")
    closure = _require_sha(source_closure_sha256, "source closure")
    env = _mapping(environment, "environment")
    _require_expected_identity(env, expected_environment, "environment")
    implementation = hlt17_implementation_identity()
    if expected_implementation is not None:
        require_hlt17_implementation_identity(
            expected_implementation, label="expected implementation"
        )
        if dict(expected_implementation) != implementation:
            raise HLT17SRCQ1REC1IdentityError("foreign implementation identity")
    static_identities = static_input_identities(static_input_bytes)
    ceilings = (
        recommended_resource_ceilings() if resource_ceilings is None else resource_ceilings
    )
    if type(ceilings) is not ResourceCeilings:
        raise HLT17SRCQ1REC1ResourceError("resource ceilings type differs")
    frozen = require_frozen_srcq1_caps(
        prospective_requested_cap_hex_by_member
        if prospective_requested_cap_hex_by_member is not None
        else FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER
    )
    builder = runtime_origin_builder or build_pro20_runtime_origin
    now = clock or time.perf_counter
    rss = rss_bytes or current_rss_bytes
    started = now()
    initial_rss = int(rss())
    runtime = runtime_origin
    if runtime is None:
        try:
            runtime = builder(
                origin,
                static_input_bytes=static_input_bytes,
                source_closure_sha256=closure,
            )
        except Pro20SourceFactoryError as error:
            raise HLT17SRCQ1REC1SourceError(
                "physical static-factory/origin construction failed"
            ) from error
    members = getattr(runtime, "members", None)
    if not isinstance(members, Mapping) or tuple(members) != MEMBER_KEYS:
        raise HLT17SRCQ1REC1IdentityError("runtime origin member order differs")
    if getattr(runtime, "captured_origin_sha256", None) != origin.capture_sha256:
        raise HLT17SRCQ1REC1IdentityError("runtime origin capture identity differs")
    if getattr(runtime, "source_closure_sha256", None) != closure:
        raise HLT17SRCQ1REC1IdentityError("runtime origin source-closure identity differs")
    requested_caps, requested_cap_facts = derive_prospective_requested_caps(
        origin,
        runtime,
        expected_cap_hex_by_member=frozen,
    )
    for key in MEMBER_KEYS:
        if requested_caps[key].hex() != frozen[key]:
            raise HLT17SRCQ1REC1CapError(
                f"{key} derived cap is not the frozen SRCQ1 hex"
            )
    evaluate = source_evaluator or (
        lambda view: evaluate_bound_source_once(
            binding=view.source_binding,
            time=float(view.member.time),
            state=view.member.state,
            transaction=view.member.transaction,
            tracers=view.member.tracers,
        )
    )
    expected_configs = (
        None
        if expected_source_configurations is None
        else _mapping(expected_source_configurations, "expected source configurations")
    )
    newly_owned: list[dict[str, object]] = []
    halted: str | None = None
    evaluated: list[str] = []
    try:
        for key in RESUME_MEMBER_KEYS:
            view = members[key]
            if getattr(view, "member_key", None) != key:
                raise HLT17SRCQ1REC1IdentityError(f"{key} runtime member identity differs")
            captured = _captured_member(origin, key)
            expected_config = None if expected_configs is None else expected_configs.get(key)
            record = qualify_member_with_overlays(
                member_key=key,
                view=view,
                captured=captured,
                initial_cap=requested_caps[key],
                evaluate=evaluate,
                now=now,
                rss=rss,
                ceilings=ceilings,
                started=started,
                expected_config=expected_config if type(expected_config) is str else None,
            )
            evaluated.append(key)
            newly_owned.append(record)
            if record.get("outcome") != "prepared":
                halted = str(record.get("halted_outcome") or record.get("outcome"))
                break
    except HLT17SRCQ1Error as error:
        key = RESUME_MEMBER_KEYS[len(newly_owned)] if len(newly_owned) < 3 else "unknown"
        newly_owned.append(
            {
                "member_key": key,
                "role": "newly_owned",
                "remeasured": True,
                "outcome": error.outcome,
                "reason": str(error),
                "requested_cap_hex": frozen.get(key),
                "accepted_state_advanced": False,
                "endpoint_adopted": False,
            }
        )
        halted = error.outcome
    except (C1R1RefinementPathStop, C1R1RuntimeResourceStop, C1R1CoordinateLatticeStop) as error:
        if isinstance(error, C1R1RuntimeResourceStop):
            outcome = "resource"
        elif isinstance(error, C1R1CoordinateLatticeStop):
            outcome = "invalid"
        elif isinstance(error, C1R1RefinementPathStop) and isinstance(
            error.cause, CFLRetryRequired
        ):
            try:
                require_independent_cfl_ratio(
                    {
                        "kind": "cfl_retry",
                        "observed_ratio_hex": float(error.cause.observed_ratio).hex(),
                        "maximum_ratio_hex": float(error.cause.maximum_ratio).hex(),
                    }
                )
                outcome = "cfl"
            except HLT17SRCQ1REC1InconclusiveError:
                outcome = "inconclusive"
            except HLT17SRCQ1Error:
                outcome = "invalid"
        else:
            outcome = "source"
        key = RESUME_MEMBER_KEYS[len(newly_owned)] if len(newly_owned) < 3 else "unknown"
        newly_owned.append(
            {
                "member_key": key,
                "role": "newly_owned",
                "remeasured": True,
                "outcome": outcome,
                "reason": str(error),
                "requested_cap_hex": frozen.get(key),
                "accepted_state_advanced": False,
                "endpoint_adopted": False,
            }
        )
        halted = outcome
    except HLT17OwnerRetryExhausted as error:
        outcome = error.owner if error.owner in {"source", "cfl", "temporal"} else "invalid"
        key = RESUME_MEMBER_KEYS[len(newly_owned)] if len(newly_owned) < 3 else "unknown"
        newly_owned.append(
            {
                "member_key": key,
                "role": "newly_owned",
                "outcome": outcome,
                "reason": str(error),
                "accepted_state_advanced": False,
                "endpoint_adopted": False,
            }
        )
        halted = outcome
    combined: list[dict[str, object]] = []
    predecessor_by_key = {item["member_key"]: item for item in predecessor}
    newly_by_key = {item["member_key"]: item for item in newly_owned}
    for key in MEMBER_KEYS:
        if key in predecessor_by_key:
            combined.append(predecessor_by_key[key])
        elif key in newly_by_key:
            combined.append(newly_by_key[key])
        else:
            combined.append(
                {
                    "member_key": key,
                    "role": "not_reached",
                    "remeasured": False,
                    "source_evaluated_by_rec1": False,
                    "outcome": None,
                    "accepted_state_advanced": False,
                    "endpoint_adopted": False,
                }
            )
    result = {
        "schema": SCHEMA,
        "artifact_id": ARTIFACT_ID,
        "protocol": PROTOCOL,
        "campaign_id": CAMPAIGN_ID,
        "implementation": {
            **implementation_identity(),
            **implementation,
        },
        "environment": env,
        "input": {
            "static_input_sha256": static_identities,
            "source_closure_sha256": closure,
        },
        "origin": {
            "capture_sha256": origin.capture_sha256,
            "generation1_bridge_content_id": getattr(
                origin, "generation1_bridge_content_id", None
            ),
            "origin_time_hex": ORIGIN_TIME_HEX,
        },
        "attempt2": {
            "artifact_id": ATTEMPT2_ARTIFACT_ID,
            "classification": ATTEMPT2_CLASSIFICATION,
            "compact_identity": attempt2["compact_identity"],
            "compact_sha256": attempt2["compact_sha256"],
            "sealed_compact_sha256": ATTEMPT2_SEALED_COMPACT_SHA256,
            "metadata_linked_compact_sha256": ATTEMPT2_METADATA_LINKED_COMPACT_SHA256,
            "derivation_document": attempt2["derivation_document"],
            "raw_leaf_sha256": attempt2["raw_leaf_sha256"],
            "raw_directory_payload_sha256": attempt2["raw_directory_payload_sha256"],
            "raw_leaf_bytes": attempt2["raw_leaf_bytes"],
            "retry_authorized": False,
            "cfl_terminal_independently_reproduced": False,
            "runner_label_sufficient": False,
        },
        "source_configuration": {
            key: members[key].source_configuration_sha256
            for key in MEMBER_KEYS
            if key in members and hasattr(members[key], "source_configuration_sha256")
        },
        "prospective_initial_requested_caps": requested_cap_facts,
        "resource_ceilings": ceilings.as_mapping(),
        "resource_facts": {
            "total_wall_seconds": now() - started,
            "rss_bytes_before_factory": initial_rss,
            "rss_bytes": int(rss()),
        },
        "predecessor_rk_members": predecessor,
        "members": combined,
        "newly_owned_members": newly_owned,
        "member_order": list(MEMBER_KEYS),
        "resume_member_keys": list(RESUME_MEMBER_KEYS),
        "evaluated_member_keys": evaluated,
        "rk_members_remeasured": False,
        "halted_outcome": halted,
        "completed_member_count": len(newly_owned),
        "phase": "qualification",
        "physical_source_qualification_attempted": True,
        "accepted_state_advanced": False,
        "endpoint_adopted": False,
        "store_published": False,
        "campaign_execution_authorized": False,
        "c1r1_actual_versus_imp1_reference_wire": {
            "c1r1_implementation_id": implementation["implementation_id"],
            "imp1_reference_wire_artifact_id": implementation["reference_wire_artifact_id"],
            "distinct": True,
        },
    }
    refuse_endpoint_payload(result)
    canonical_json_bytes(result)
    return result


def publish_qualification_result(
    output_directory: Path,
    result: Mapping[str, object],
    *,
    fault_hook: Callable[[str], None] | None = None,
    physical_work_started: bool | None = None,
) -> dict[str, object]:
    refuse_endpoint_payload(result)
    attempted = physical_work_started
    if attempted is None:
        attempted = result.get("physical_source_qualification_attempted", False)
    if type(attempted) is not bool:
        raise HLT17SRCQ1REC1Error(
            "physical_source_qualification_attempted must be a built-in bool"
        )
    try:
        return publish_srcq1_qualification_result(
            output_directory, result, fault_hook=fault_hook
        )
    except HLT17SRCQ1PublicationError as error:
        raise HLT17SRCQ1REC1PublicationError(
            str(error),
            published=error.published,
            physical_work_started=attempted,
        ) from error


def _call_loader(value: object, loader: Callable | None, label: str) -> object:
    if value is not None:
        return value
    if loader is None:
        raise HLT17SRCQ1REC1IdentityError(f"{label} is absent")
    return loader()


def qualify_and_publish(
    origin: object | None,
    *,
    repository_root: Path,
    output_directory: Path,
    authority_commit: str,
    static_input_bytes: Mapping[str, bytes] | None,
    source_closure_sha256: str,
    environment: Mapping[str, str] | None,
    attempt2_compact_bytes: bytes | None,
    attempt2_raw_bytes: bytes | None,
    authority_delta_validator: Callable[..., object] | None = None,
    resource_ceilings: ResourceCeilings | None = None,
    runtime_origin: object | None = None,
    runtime_origin_builder: Callable[..., object] | None = None,
    source_evaluator: Callable[..., Mapping[str, object]] | None = None,
    clock: Callable[[], float] | None = None,
    rss_bytes: Callable[[], int] | None = None,
    expected_origin_sha256: str | None = None,
    expected_environment: Mapping[str, str] | None = None,
    expected_implementation: Mapping[str, str] | None = None,
    expected_source_configurations: Mapping[str, str] | None = None,
    fault_hook: Callable[[str], None] | None = None,
    origin_loader: Callable[[], object] | None = None,
    static_input_loader: Callable[[], Mapping[str, bytes]] | None = None,
    environment_loader: Callable[[], Mapping[str, str]] | None = None,
    attempt2_compact_loader: Callable[[], bytes] | None = None,
    attempt2_raw_loader: Callable[[], bytes] | None = None,
) -> dict[str, object]:
    root = Path(repository_root)
    implementation = implementation_identity(root)
    digest = _require_sha(implementation.get("sha256"), "implementation sha256")
    authority = require_authority_delta(
        repository_root=root,
        authority_commit=authority_commit,
        implementation_sha256=digest,
        validator=authority_delta_validator,
    )
    receipt_closure = _require_sha(
        authority.get("source_closure_sha256"), "source_closure_sha256"
    )
    if source_closure_sha256 != receipt_closure:
        raise HLT17SRCQ1REC1AuthoritySeamError(
            "source_closure_sha256 differs from authority receipt"
        )
    requested_caps = authority.get("prospective_requested_cap_hex_by_member")
    if not isinstance(requested_caps, Mapping):
        raise HLT17SRCQ1REC1AuthoritySeamError(
            "authority receipt omits prospective requested-cap mapping"
        )
    require_frozen_srcq1_caps(requested_caps)
    expected_output = authority.get("output_namespace")
    if type(expected_output) is not str or not expected_output:
        raise HLT17SRCQ1REC1AuthoritySeamError(
            "authority receipt omits the qualification output namespace"
        )
    output = Path(output_directory)
    expected_path = Path(expected_output)
    expected_path = expected_path if expected_path.is_absolute() else root / expected_path
    if not output.is_absolute() or output != expected_path:
        raise HLT17SRCQ1REC1AuthoritySeamError(
            "qualification output directory differs from authority receipt"
        )
    phase = "repository_input"
    physical_work_started = False
    loaded_origin = origin
    loaded_environment: Mapping[str, str] | dict[str, object] = environment or {}
    try:
        compact = _call_loader(
            attempt2_compact_bytes, attempt2_compact_loader, "Attempt2 compact"
        )
        raw = _call_loader(attempt2_raw_bytes, attempt2_raw_loader, "Attempt2 raw")
        loaded_environment = _call_loader(
            environment, environment_loader, "environment"
        )
        loaded_origin = _call_loader(origin, origin_loader, "origin capture")
        static_inputs = _call_loader(
            static_input_bytes, static_input_loader, "static factory inputs"
        )
        phase = "qualification"
        physical_work_started = True
        result = qualify_srcq1_rec1(
            loaded_origin,
            static_input_bytes=static_inputs,  # type: ignore[arg-type]
            source_closure_sha256=receipt_closure,
            environment=loaded_environment,  # type: ignore[arg-type]
            attempt2_compact_bytes=compact,  # type: ignore[arg-type]
            attempt2_raw_bytes=raw,  # type: ignore[arg-type]
            resource_ceilings=resource_ceilings,
            runtime_origin=runtime_origin,
            runtime_origin_builder=runtime_origin_builder,
            source_evaluator=source_evaluator,
            clock=clock,
            rss_bytes=rss_bytes,
            expected_origin_sha256=expected_origin_sha256,
            expected_environment=expected_environment,
            expected_implementation=expected_implementation,
            expected_source_configurations=expected_source_configurations,
            prospective_requested_cap_hex_by_member=requested_caps,
        )
    except UnsafePathError:
        typed = HLT17SRCQ1REC1IdentityError(f"{phase} path is absent or unsafe")
        typed.phase = phase
        typed.physical_work_started = physical_work_started
        result = {
            "schema": SCHEMA,
            "artifact_id": ARTIFACT_ID,
            "protocol": PROTOCOL,
            "campaign_id": CAMPAIGN_ID,
            "halted_outcome": typed.outcome,
            "reason": str(typed),
            "implementation": dict(implementation),
            "environment": dict(loaded_environment),
            "input": {"source_closure_sha256": receipt_closure},
            "origin": {"capture_sha256": getattr(loaded_origin, "capture_sha256", None)},
            "resource_ceilings": (
                recommended_resource_ceilings()
                if resource_ceilings is None
                else resource_ceilings
            ).as_mapping(),
            "completed_member_count": 0,
            "rk_members_remeasured": False,
            "phase": phase,
            "physical_source_qualification_attempted": physical_work_started,
            "accepted_state_advanced": False,
            "endpoint_adopted": False,
            "store_published": False,
            "campaign_execution_authorized": False,
        }
    except Pro20OriginError as error:
        typed = HLT17SRCQ1REC1IdentityError("origin capture is an invalid premise")
        typed.phase = "repository_input"
        typed.physical_work_started = False
        result = {
            "schema": SCHEMA,
            "artifact_id": ARTIFACT_ID,
            "protocol": PROTOCOL,
            "campaign_id": CAMPAIGN_ID,
            "halted_outcome": typed.outcome,
            "reason": str(error),
            "implementation": dict(implementation),
            "environment": dict(loaded_environment),
            "input": {"source_closure_sha256": receipt_closure},
            "origin": {"capture_sha256": getattr(loaded_origin, "capture_sha256", None)},
            "resource_ceilings": (
                recommended_resource_ceilings()
                if resource_ceilings is None
                else resource_ceilings
            ).as_mapping(),
            "completed_member_count": 0,
            "rk_members_remeasured": False,
            "phase": "repository_input",
            "physical_source_qualification_attempted": False,
            "accepted_state_advanced": False,
            "endpoint_adopted": False,
            "store_published": False,
            "campaign_execution_authorized": False,
        }
    except HLT17SRCQ1Error as error:
        if not getattr(error, "phase", None) or error.phase == "invalid":
            error.phase = phase
        if not getattr(error, "physical_work_started", False):
            error.physical_work_started = physical_work_started
        result = {
            "schema": SCHEMA,
            "artifact_id": ARTIFACT_ID,
            "protocol": PROTOCOL,
            "campaign_id": CAMPAIGN_ID,
            "halted_outcome": error.outcome,
            "reason": str(error),
            "implementation": dict(implementation),
            "environment": dict(loaded_environment),
            "input": {"source_closure_sha256": receipt_closure},
            "origin": {"capture_sha256": getattr(loaded_origin, "capture_sha256", None)},
            "resource_ceilings": (
                recommended_resource_ceilings()
                if resource_ceilings is None
                else resource_ceilings
            ).as_mapping(),
            "completed_member_count": 0,
            "rk_members_remeasured": False,
            "phase": error.phase,
            "physical_source_qualification_attempted": error.physical_work_started,
            "accepted_state_advanced": False,
            "endpoint_adopted": False,
            "store_published": False,
            "campaign_execution_authorized": False,
        }
    result["authority_commit"] = _require_commit(authority_commit, "authority_commit")
    result["authority_delta"] = authority
    result["implementation"]["sha256"] = digest
    result.setdefault("phase", phase)
    result.setdefault(
        "physical_source_qualification_attempted", physical_work_started
    )
    tracked_attempt = result.get(
        "physical_source_qualification_attempted", physical_work_started
    )
    if type(tracked_attempt) is not bool:
        raise HLT17SRCQ1REC1Error(
            "physical_source_qualification_attempted must be a built-in bool"
        )
    result["physical_source_qualification_attempted"] = tracked_attempt
    try:
        publication = publish_qualification_result(
            Path(output_directory),
            result,
            fault_hook=fault_hook,
            physical_work_started=tracked_attempt,
        )
    except HLT17SRCQ1REC1PublicationError as error:
        error.phase = "publication"
        error.physical_work_started = tracked_attempt
        raise
    result["publication"] = publication
    return result


__all__ = [
    "ARTIFACT_ID",
    "ATTEMPT2_CLASSIFICATION",
    "ATTEMPT2_COMPACT_SHA256",
    "ATTEMPT2_DERIVATION_DOCUMENT",
    "ATTEMPT2_METADATA_LINKED_COMPACT_COMMIT",
    "ATTEMPT2_METADATA_LINKED_COMPACT_SHA256",
    "ATTEMPT2_RAW_LEAF_SHA256",
    "ATTEMPT2_SEALED_COMPACT_COMMIT",
    "ATTEMPT2_SEALED_COMPACT_SHA256",
    "AUTHORITY_SEAM",
    "AUTHORITY_VALIDATOR_MODULE",
    "AUTHORITY_VALIDATOR_NAME",
    "FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER",
    "HLT17SRCQ1REC1AuthoritySeamError",
    "HLT17SRCQ1REC1CapError",
    "HLT17SRCQ1REC1Error",
    "HLT17SRCQ1REC1IdentityError",
    "HLT17SRCQ1REC1InconclusiveError",
    "HLT17SRCQ1REC1PublicationError",
    "HLT17SRCQ1REC1ResourceError",
    "HLT17SRCQ1REC1SourceError",
    "PREDECESSOR_MEMBER_KEYS",
    "RESUME_MEMBER_KEYS",
    "SCHEMA",
    "authenticate_attempt2_predecessor",
    "classify_prepare_outcome",
    "implementation_identity",
    "isolate_working_checkpoint",
    "load_coordinator_authority_delta_validator",
    "observe_environment",
    "plan_status",
    "predecessor_rk_evidence_references",
    "publish_qualification_result",
    "qualify_and_publish",
    "qualify_srcq1_rec1",
    "recommended_resource_ceilings",
    "require_authority_delta",
    "require_frozen_srcq1_caps",
    "require_independent_cfl_ratio",
    "retain_bridge_result",
]
