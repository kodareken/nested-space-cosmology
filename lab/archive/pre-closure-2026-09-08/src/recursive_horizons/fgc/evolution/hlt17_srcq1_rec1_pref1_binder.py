"""Independent no-write binder for the completed SRCQ1-REC1 qualification.

The live path authenticates authority ``be6a179`` and parent A2, the freeze
image, Attempt2 identities, the REC1 raw leaf, the historical runs baseline,
Planck hashes, and PRO20/state absence. It then rebuilds origin/runtime from
pro20/H17/C1R1 owners, derives the six caps, binds Attempt2 RK receipts
without remeasurement, and reconstructs the three SSP source/C1R1 paths in
memory. Comparison uses independent facts, never runner
outcome/disposition/admission labels.

Ordinary compact verification is raw-, runs-, source-, Git-, and shadow-blind.
This module does not import REC1 runner, authority classification, or
publication decisions, and it never serializes an endpoint or writes state.
"""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from hashlib import sha256
import json
import math
import os
from pathlib import Path, PurePosixPath
import platform
import re
import stat
import sys
import tomllib
from types import MappingProxyType
from typing import Any, Callable, Mapping

from recursive_horizons.evidence_io import (
    CanonicalJSONError,
    CommitParents,
    EvidenceIOError,
    InspectDelta,
    InspectTree,
    ReadBlob,
    UnsafePathError,
    canonical_json_bytes,
    git_read,
    load_canonical_json,
    read_regular_file,
)

from .hlt17_admission_runtime import (
    C1R1CoordinateLatticeStop,
    C1R1RefinementPathStop,
    C1R1RuntimeResourceStop,
    hlt17_implementation_identity,
    require_hlt17_prepared,
)
from .hlt17_imp1_cursor import (
    HLT17_CFL_RETRY_CAP,
    HLT17_SOURCE_RETRY_CAP,
    encode_tdg7_plan,
)
from .hlt17_member_checkpoint import HLT17InMemoryCheckpoint
from .hlt17_member_checkpoint import _isolated_member as isolate_hlt17_runtime_member
from .hlt17_srcq1 import (
    FORBIDDEN_ENDPOINT_KEYS,
    ResourceCeilings,
    check_inherited_counters,
    check_method_owned_identity,
    current_rss_bytes,
    member_fingerprints,
    observe_environment,
    recommended_resource_ceilings,
    refuse_endpoint_payload,
    static_input_identities,
)
from .numerical_engine import EvolutionRHS, EvolutionState, array_content_sha256
from .pro20_origin import (
    MEMBER_KEYS,
    Pro20OriginCapture,
    Pro20OriginError,
    revalidate_pro20_origin_capture,
)
from .pro20_source_factory import (
    Pro20SourceFactoryError,
    build_pro20_runtime_origin,
)
from .proto5_runtime import CFLRetryRequired
from .tdg6_temporal_admission_design import (
    TDG6_COMPLETE_STATE_CHANNELS,
    TDG6_PASSING_CLASSES,
)


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-HLT17-SRCQ1-REC1-PREF1"
CLASSIFICATION = (
    "independently_bound_all_six_source_c1r1_qualified_no_state_advance"
)
CONFIG_PATH = "configs/fgc/fgc-1-hlt17-srcq1-rec1-pref1.toml"
RESULT_PATH = "results/fgc-1-hlt17-srcq1-rec1-pref1.json"
OWNER_DOCUMENT = "docs/fgc-hlt17-srcq1-rec1-pref1.md"
FREEZE_ARTIFACT_ID = "FGC-1-HLT17-SRCQ1-REC1-FRZ1"
FREEZE_CONFIG_PATH = "configs/fgc/fgc-1-hlt17-srcq1-rec1-frz1.toml"
FREEZE_OWNER_PATH = "docs/fgc-hlt17-srcq1-rec1-frz1.md"
FREEZE_MAKE_PATH = "mk/current-foundation.mk"
AUTHORITY_COMMIT = "be6a1792a30c39fc48ea7435368cbc8c17718d4c"
AUTHORITY_PARENT = "1fc17e5e2327cb3d8776386fc12d92a48af07f9c"
ALLOWED_AUTHORITY_DELTA = (
    FREEZE_CONFIG_PATH,
    FREEZE_OWNER_PATH,
    FREEZE_MAKE_PATH,
)
FREEZE_CONFIG_SHA256 = (
    "6a11a65dd7de04b36dcf02c5b1774e031efe7cb1638b0445c465045086989525"
)
FREEZE_OWNER_SHA256 = (
    "b4955efca7b82e8add4eb714d495ff317e2a7849f5b7a3250c26f59e4a9a5be5"
)
FREEZE_MAKE_SHA256 = (
    "f2780359f29880af949891a5b790a00ff09de11895d020c9f1cbc62a2f126440"
)
IMPLEMENTATION_SHA256 = (
    "faf15958643be946d91be235690e5ed5a2b8e89a0b6d0845288956b70064c2dc"
)
RUNNER_SHA256 = (
    "fd5947022a8322190449a8dfa4f7422c02db2438ce8d656c4b67916fe853da3c"
)
REC1_AUTHORITY_SHA256 = (
    "294ae85eb6f72a10da7e759be69b5111124208b4a54e3c77ad1bf0941f7775b4"
)
SOURCE_CLOSURE_SHA256 = (
    "48e507ca674c7d164b905bb47d9c6778902886979ddb3dbe906e7d5cfd258ce0"
)
SOURCE_CLOSURE_FILE_COUNT = 248
ORIGIN_CAPTURE_SHA256 = (
    "2824af06c674336d479ec53faa4269f02fb598055860ec22f2ad84602a83fe72"
)
SOURCE_CLOSURE_DOMAIN = b"FGC-1-HLT17-SRCQ1-REC1-SOURCE-CLOSURE-v1\n"
SOURCE_TREE = "src/recursive_horizons"
REC1_HARNESS_PATH = "src/recursive_horizons/fgc/evolution/hlt17_srcq1_rec1.py"
REC1_RUNNER_PATH = "scripts/qualify_fgc_hlt17_srcq1_rec1.py"
REC1_AUTHORITY_PATH = (
    "src/recursive_horizons/fgc/evolution/hlt17_srcq1_rec1_auth1.py"
)

ATTEMPT2_ARTIFACT_ID = "FGC-1-HLT17-SRCQ1-ATTEMPT2"
ATTEMPT2_COMPACT_RELATIVE = "results/fgc-1-hlt17-srcq1-attempt2.json"
ATTEMPT2_RAW_DIRECTORY = (
    "runs/fgc-2-sf1/hlt17-srcq1-source-origin-qualification"
)
ATTEMPT2_RAW_LEAF_RELATIVE = f"{ATTEMPT2_RAW_DIRECTORY}/qualification.json"
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

REC1_OUTPUT_NAMESPACE = (
    "runs/fgc-2-sf1/hlt17-srcq1-rec1-source-origin-qualification"
)
REC1_RAW_LEAF_RELATIVE = f"{REC1_OUTPUT_NAMESPACE}/qualification.json"
REC1_RAW_LEAF_SHA256 = (
    "ec9cd18b48ebfae93508182444ed5b689c566a7bc4cadcdc6a002d3550b17f16"
)
REC1_RAW_LEAF_BYTES = 105375
REC1_RAW_DIRECTORY_PAYLOAD_SHA256 = (
    "567d14095cc8631ed04715b332a4a35968d5a43dee1983bf062dc7647aa823ea"
)
PRO20_EVENT_NAMESPACE = "runs/fgc-2-sf1/pro20-event1/calibration"
HISTORICAL_RUNS_ROOT = "runs"
HISTORICAL_RUNS_FILE_COUNT = 300
HISTORICAL_RUNS_BYTE_COUNT = 150437927
HISTORICAL_RUNS_DIGEST = (
    "77aa79f0300bd6866fcb491fe38aabf813a0064cbc7d511ef067364b4998ed88"
)
PLANCK_FILES = MappingProxyType(
    {
        "collected-data/smica_2048.fits": (
            2013312960,
            "60952c645eb33d151905ddf5837477e15ca02a3261feca4ceca3c3feece4f9ac",
        ),
        "collected-data/mask_common.fits": (
            201335040,
            "23f49a3479073408b83a7f1bb2cca490d9be7f934b2e0e891524b6a1d06a5fc4",
        ),
    }
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
SOURCE_DIAGNOSTIC_KEYS = (
    "source_residual_infinity",
    "source_raw_gate_passed",
    "source_refinement_iterations",
    "source_residual_decreased_monotonically",
    "kinetic_condition_infinity",
    "reduction_constraint_infinity",
    "coordinate_speed_upper",
    "minimum_lapse",
    "minimum_radial_metric",
    "minimum_areal_radius_away_from_center",
)
MEMBER_CAP_INPUTS = (
    (
        "RK4-2049",
        "0x1.0000000000000p-4",
        "0x1.33345e70b5188p+0",
        "0x1.aaa90b0fb5c26p-8",
    ),
    (
        "RK4-4097",
        "0x1.0000000000000p-5",
        "0x1.3333333aecf3ep+0",
        "0x1.aaaaaa9fefc9cp-9",
    ),
    (
        "RK4-8193",
        "0x1.0000000000000p-6",
        "0x1.3333333333653p+0",
        "0x1.aaaaaaaaaa654p-10",
    ),
    (
        "SSPRK3-4097",
        "0x1.0000000000000p-5",
        "0x1.3333333fe51c4p+0",
        "0x1.aaaaaa9908e70p-9",
    ),
    (
        "SSPRK3-8193",
        "0x1.0000000000000p-6",
        "0x1.33333333349e5p+0",
        "0x1.aaaaaaaaa8b26p-10",
    ),
    (
        "SSPRK3-16385",
        "0x1.0000000000000p-7",
        "0x1.3333333333346p+0",
        "0x1.aaaaaaaaaaa91p-11",
    ),
)
SOURCE_CONFIGURATION_SHA256_BY_MEMBER = MappingProxyType(
    {
        "RK4-2049": (
            "76e0cae19347d29f3958decf115df9a20953b0d831cda904b41aa14f017e7565"
        ),
        "RK4-4097": (
            "8abbdc5978e989262855281ba37e783642adac354d85fca57aa931107dbbd860"
        ),
        "RK4-8193": (
            "55b8247845791d5aa5c5131a521f668fa830628583059fb5ea55e2fa1aa82eb5"
        ),
        "SSPRK3-4097": (
            "d39b7bbb33dd3c8ffa8985d76d9b1deb5249a725f73ed158ca46e01b105c1dd8"
        ),
        "SSPRK3-8193": (
            "25671dbf22da99f1662e4fc209d2680dc0adc8f2e58b0f3d983999d705ba5f9c"
        ),
        "SSPRK3-16385": (
            "3f97c57ec9514fa7bef695a8f967110bd1dfa8141d1be8acdbc07190d3ccae61"
        ),
    }
)
ATTEMPT2_RK_RECEIPTS = MappingProxyType(
    {
        "RK4-2049": (
            "25e6e1de623d5f2b75ca8c5610cf32112dfd7795996db2dbc8289e394211d8cb"
        ),
        "RK4-4097": (
            "005577c7dddf10bac0793fbf26c02b766bd931064df9e04726149c551bc83523"
        ),
        "RK4-8193": (
            "f3461dbd9c5075ac316269b89ddc9f34bb48593897db580dbc99d82f9fc936f2"
        ),
    }
)
EXPECTED_SSP4097_CFL = MappingProxyType(
    {
        "kind": "cfl_retry",
        "observed_ratio_hex": "0x1.000000003e80dp-3",
        "maximum_ratio_hex": "0x1.0000000000000p-3",
        "observed_exceeds_maximum": True,
        "runner_label_sufficient": False,
    }
)
EXPECTED_SSP4097_NEXT_PLAN_CAP_HEX = "0x1.aaaaaa9908000p-10"
EXPECTED_C1R1_RECEIPTS = MappingProxyType(
    {
        "SSPRK3-4097": (
            "1a97dab121e9c3f3134c50286bd46505d519dc8728dee9d23963629f4fca56e2"
        ),
        "SSPRK3-8193": (
            "523a990cd60dabb71ba8df55c53c634c9ad21fa9a2631bcf754b0fae0166e8a6"
        ),
        "SSPRK3-16385": (
            "b8b53246970d4f38016488304bf90c4c00d0c7e9529540d23ac5cefaa40d3159"
        ),
    }
)
EXPECTED_C1R1_ASSESSMENTS = MappingProxyType(
    {
        "SSPRK3-4097": (
            "479377cf5d50966739e161fba636a6e55c9d79b65da7e3becea45fbb55cf1e80"
        ),
        "SSPRK3-8193": (
            "752c90bc17244344bbebacbb4a7cf911e5f6b22fbd1659d71e1fff8bf0367e2c"
        ),
        "SSPRK3-16385": (
            "724a37a6c1c9c6c24e0103519b0955aa29d332f28d2c6d4a1763aec69e1e3709"
        ),
    }
)
EXPECTED_C1R1_FAMILIES = MappingProxyType(
    {
        "SSPRK3-4097": (
            "3ef7750c302c3c73580a44ab0b6d7c546d8e8e6c829eb3f68531dfd07087f6e7"
        ),
        "SSPRK3-8193": (
            "6d6b7d7c99a37b05432d9382dd23c0a32db0b4d6dcb74eae15b9a23ae158afa9"
        ),
        "SSPRK3-16385": (
            "1912556f35887fede2f16bc556e547fe1823d1801cc1655c790c96b122412b11"
        ),
    }
)
RECORDED_RAW_WALL_SECONDS = 371.6832666249975
RECORDED_RAW_RSS_BYTES = 1133412352
CFL_MAXIMUM_HEX = "0x1.0000000000000p-3"
PREPARED_DISPOSITIONS = frozenset(
    {"prepared_fresh", "prepared_overlay", "prepared_successor"}
)
CFL_REQUIRED = "cfl_retry_required"
CFL_EXHAUSTED = "cfl_retry_exhausted"
SOURCE_REQUIRED = "source_retry_required"
SOURCE_EXHAUSTED = "source_retry_exhausted"
CFL_DISPOSITIONS = frozenset({CFL_REQUIRED, CFL_EXHAUSTED})
SOURCE_DISPOSITIONS = frozenset({SOURCE_REQUIRED, SOURCE_EXHAUSTED})
OVERLAY_CONTINUE = frozenset({CFL_REQUIRED, SOURCE_REQUIRED})
MAX_PREPARE_ATTEMPTS = (HLT17_SOURCE_RETRY_CAP + 1) + (HLT17_CFL_RETRY_CAP + 1)
ENVIRONMENT_KEYS = (
    "python_implementation",
    "python_version",
    "numpy_version",
    "blas_name",
    "blas_version",
    "system",
    "machine",
    "byteorder",
    "executable",
    "executable_realpath",
    "executable_sha256",
    "numpy_extension_realpath",
    "numpy_extension_sha256",
    "kernel_release",
)
_COMMIT_RE = re.compile(r"\A(?:[0-9a-f]{40}|[0-9a-f]{64})\Z")
_SHA_RE = re.compile(r"\A[0-9a-f]{64}\Z")
_REGULAR_GIT_MODES = frozenset({"100644", "100755"})
COMPACT_RESULT_SHA256: str | None = (
    "ea69bc1fd7ee31c41982f7eabfc3bf3e93c7e237765e71e5b112dadf69b34e75"
)
ABSENT_COMPACT_MESSAGE = (
    "FGC-1-HLT17-SRCQ1-REC1-PREF1 compact result is absent at "
    f"{RESULT_PATH}. Default --check validates tracked compact bytes only "
    "and never reconstructs, reads raw/runs, or inspects Git."
)
NONCLAIMS = (
    "PREF1 licenses only a separately frozen PRO20 first-event authority.",
    "PREF1 does not authorize PRO20 event execution.",
    "PREF1 does not authorize calibration, SGB-L health, candidate, "
    "mechanism, or physics claims.",
    "PREF1 does not advance accepted state, write a campaign store, or "
    "serialize an endpoint.",
    "Runner outcome, disposition, and admission labels are never sufficient.",
)
EXPECTED_CLAIMS = MappingProxyType(
    {
        "independently_bound_all_six_source_c1r1_qualified": True,
        "accepted_state_advanced": False,
        "campaign_store_written": False,
        "endpoint_adopted_or_serialized": False,
        "pro20_namespace_absent": True,
        "rk_members_remeasured": False,
        "runner_label_sufficient": False,
        "pro20_first_event_authority_licensed_separately": True,
        "pro20_event_execution_authorized": False,
        "calibration_authorized": False,
        "sgbl_health_established": False,
        "candidate_execution_authorized": False,
        "mechanism_claim_authorized": False,
        "physical_claim_authorized": False,
    }
)


class HLT17SRCQ1REC1PREF1Error(RuntimeError):
    """The independent REC1-PREF1 certificate could not be established."""


class HLT17SRCQ1REC1PREF1ResourceError(HLT17SRCQ1REC1PREF1Error):
    """Live reconstruction exceeded a freeze resource ceiling."""


def _fail(message: str) -> None:
    raise HLT17SRCQ1REC1PREF1Error(message)


def _sha256_hex(raw: bytes) -> str:
    if type(raw) is not bytes:
        _fail("payload is not immutable bytes")
    return sha256(raw).hexdigest()


def _require_sha(value: object, label: str) -> str:
    if type(value) is not str or _SHA_RE.fullmatch(value) is None:
        _fail(f"{label} must be a lowercase SHA-256 digest")
    return value


def _require_commit(value: object, label: str) -> str:
    if type(value) is not str or _COMMIT_RE.fullmatch(value) is None:
        _fail(f"{label} must be an exact lowercase commit")
    return value


def _mapping(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping) or type(value) is bool:
        _fail(f"{label} must be a mapping")
    return dict(value)


def _json_mapping(raw: bytes, label: str) -> dict[str, Any]:
    if type(raw) is not bytes:
        _fail(f"{label} is not immutable bytes")
    try:
        parsed = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise HLT17SRCQ1REC1PREF1Error(f"{label} is not JSON") from error
    return _mapping(parsed, label)


def _require_bool(value: object, expected: bool, label: str) -> bool:
    if type(value) is not bool or value is not expected:
        _fail(f"{label} differs")
    return value


def _root(value: Path) -> Path:
    if not isinstance(value, Path):
        value = Path(value)
    if not value.is_absolute():
        _fail("repository root is not canonical")
    try:
        resolved = value.resolve(strict=True)
    except OSError as error:
        raise HLT17SRCQ1REC1PREF1Error(
            "repository root cannot be resolved"
        ) from error
    if resolved != value:
        _fail("repository root traverses a symlink")
    return value


def _same(actual: object, expected: object, label: str) -> None:
    if actual != expected:
        _fail(f"{label} differs")


def canonical_result(value: object) -> bytes:
    _refuse(value)
    return canonical_json_bytes(value) + b"\n"


def _refuse(value: object) -> None:
    try:
        refuse_endpoint_payload(value)
    except Exception as error:
        raise HLT17SRCQ1REC1PREF1Error(
            "result serializes a forbidden endpoint"
        ) from error


def _positive_finite_hex(value: object, label: str) -> float:
    if type(value) is not str or value.lower() != value:
        _fail(f"{label} is not canonical lowercase binary64 hex")
    try:
        number = float.fromhex(value)
    except ValueError as error:
        raise HLT17SRCQ1REC1PREF1Error(f"{label} is malformed") from error
    if not math.isfinite(number) or number <= 0.0:
        _fail(f"{label} is not a finite positive binary64")
    if number.hex() != value:
        _fail(f"{label} is not a round-trip binary64")
    return number


def independent_cfl_ratio(evidence: object) -> dict[str, object]:
    """Require kind cfl_retry and observed_ratio>maximum_ratio. Labels never prove CFL."""

    mapping = _mapping(evidence, "CFL evidence")
    observed_hex = mapping.get("observed_ratio_hex")
    maximum_hex = mapping.get("maximum_ratio_hex")
    if type(observed_hex) is not str or type(maximum_hex) is not str:
        _fail("CFL stop omitted observed/max ratio evidence")
    if mapping.get("kind") != "cfl_retry":
        _fail("CFL evidence kind is not cfl_retry")
    observed = _positive_finite_hex(observed_hex, "observed_ratio_hex")
    maximum = _positive_finite_hex(maximum_hex, "maximum_ratio_hex")
    if not (observed > maximum):
        _fail("CFL stop does not independently satisfy observed_ratio>maximum_ratio")
    record = {
        "kind": "cfl_retry",
        "observed_ratio_hex": observed.hex(),
        "maximum_ratio_hex": maximum.hex(),
        "observed_exceeds_maximum": True,
        "runner_label_sufficient": False,
    }
    _refuse(record)
    return record


def independent_assessment_decisions(assessment: Mapping[str, object]) -> list[dict[str, object]]:
    channels = assessment.get("channels")
    if type(channels) is not list or len(channels) != len(
        TDG6_COMPLETE_STATE_CHANNELS
    ):
        _fail("C1R1 assessment is not eighteen complete-state channels")
    records: list[dict[str, object]] = []
    for expected, item in zip(TDG6_COMPLETE_STATE_CHANNELS, channels, strict=True):
        channel = _mapping(item, "channel assessment")
        if channel.get("channel") != expected:
            _fail("eighteen-channel assessment order differs")
        corrected = _mapping(channel.get("corrected"), "corrected path")
        decision = _mapping(corrected.get("decision"), "channel decision")
        record = {
            "channel": expected,
            "admission_passed": channel.get("admission_passed"),
            "classification": decision.get("classification"),
            "public_fine_debit": channel.get("public_fine_debit"),
            "gate_debit": channel.get("gate_debit"),
            "extra_enclosure_debit": channel.get("extra_enclosure_debit"),
        }
        records.append(record)
    return records


def independent_c1r1_record(prepared: object) -> dict[str, object]:
    bound = require_hlt17_prepared(prepared)
    identity = hlt17_implementation_identity()
    if bound.implementation_id != identity["implementation_id"]:
        _fail("prepared implementation is not the fixed C1R1 object")
    if bound.reference_wire_evaluator_id != identity["reference_wire_evaluator_id"]:
        _fail("prepared IMP1 reference wire differs")
    assessment = bound.assessment.as_mapping()
    decisions = independent_assessment_decisions(assessment)
    independent_channel_admission(decisions)
    record = {
        "prepared": True,
        "implementation_identity": dict(identity),
        "imp1_reference_wire": {
            "reference_wire_artifact_id": identity["reference_wire_artifact_id"],
            "reference_wire_evaluator_id": identity["reference_wire_evaluator_id"],
        },
        "c1r1_is_not_imp1_execution": True,
        "plan": encode_tdg7_plan(bound.plan),
        "eighteen_channel_decisions": decisions,
        "admission_passed": bound.assessment.admission_passed,
        "public_fine_debits": list(assessment["public_fine_debits"]),
        "receipt_sha256": bound.receipt_sha256,
        "assessment_sha256": bound.assessment.assessment_sha256,
        "family_sha256": bound.assessment.family_sha256,
        "corrected_endpoint_committable": False,
        "outer_or_medium_committable": False,
        "endpoint_adopted": False,
        "campaign_execution_authorized": False,
    }
    _refuse(record)
    return record


def _identity_hash(value: object, label: str) -> str:
    if hasattr(value, "as_mapping") and callable(value.as_mapping):
        payload = value.as_mapping()
    elif is_dataclass(value) and not isinstance(value, type):
        payload = asdict(value)
    elif isinstance(value, Mapping):
        payload = dict(value)
    else:
        _fail(f"{label} cannot be fingerprinted")
    return _sha256_hex(canonical_json_bytes(payload))


def _tracer_fingerprint(tracers: object) -> str:
    positions = getattr(tracers, "positions", None)
    proper_times = getattr(tracers, "proper_times", None)
    if positions is None or proper_times is None:
        _fail("tracer fingerprint requires positions and proper times")
    return array_content_sha256(positions, proper_times)


def independent_source_evaluation(
    *,
    binding: object,
    time: float,
    state: EvolutionState,
    transaction: object,
    tracers: object,
) -> dict[str, object]:
    if hasattr(binding, "validate"):
        binding.validate()
    if type(state) is not EvolutionState:
        _fail("source evaluation state type differs")
    before_state = array_content_sha256(state.u, state.p, state.q)
    before_monitor = _identity_hash(transaction.state, "monitor")
    before_causal = _identity_hash(transaction.causal_state, "causal")
    before_tracer = _tracer_fingerprint(tracers)
    before_config = getattr(binding, "configuration_sha256", None)
    rhs = binding.rhs(time, state)
    if type(rhs) is not EvolutionRHS:
        _fail("bound source did not return EvolutionRHS")
    if array_content_sha256(state.u, state.p, state.q) != before_state:
        _fail("bound source mutated its input state")
    if _identity_hash(transaction.state, "monitor") != before_monitor:
        _fail("bound source mutated monitor state")
    if _identity_hash(transaction.causal_state, "causal") != before_causal:
        _fail("bound source mutated causal state")
    if _tracer_fingerprint(tracers) != before_tracer:
        _fail("bound source mutated tracer state")
    if before_config is not None and binding.configuration_sha256 != before_config:
        _fail("source configuration changed during evaluation")
    if hasattr(binding, "validate"):
        binding.validate()
    diagnostics = dict(rhs.diagnostics)
    missing = [name for name in SOURCE_DIAGNOSTIC_KEYS if name not in diagnostics]
    if missing:
        _fail("source diagnostics are incomplete: " + ",".join(missing))
    if diagnostics["source_raw_gate_passed"] is not True:
        _fail("source raw gate failed")
    if diagnostics["source_residual_decreased_monotonically"] is not True:
        _fail("source residual monotonicity failed")
    for name in (
        "source_residual_infinity",
        "kinetic_condition_infinity",
        "reduction_constraint_infinity",
        "coordinate_speed_upper",
        "minimum_lapse",
        "minimum_radial_metric",
        "minimum_areal_radius_away_from_center",
    ):
        value = diagnostics[name]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            _fail(f"{name} is not a finite scalar")
        number = float(value)
        if not math.isfinite(number):
            _fail(f"{name} is not finite")
        diagnostics[name] = number
    for name in (
        "minimum_lapse",
        "minimum_radial_metric",
        "minimum_areal_radius_away_from_center",
    ):
        if diagnostics[name] <= 0.0:
            _fail(f"{name} is not positive")
    encoded = {
        key: (
            value.hex()
            if type(value) is float
            else value
            if type(value) in (bool, int, str)
            else str(value)
        )
        for key, value in diagnostics.items()
        if type(key) is str
    }
    return {
        "diagnostics": encoded,
        "source_health": "recorded",
        "raw_health": "recorded",
        "kinetic_health": "recorded",
        "causal_sha256": before_causal,
        "constraint_health": "recorded",
        "input_state_sha256": before_state,
        "monitor_sha256": before_monitor,
        "tracer_sha256": before_tracer,
        "mutated_input": False,
        "mutated_monitor": False,
        "mutated_causal": False,
        "mutated_tracer": False,
    }


def independent_channel_admission(decisions: object) -> dict[str, object]:
    if type(decisions) is not list or len(decisions) != len(
        TDG6_COMPLETE_STATE_CHANNELS
    ):
        _fail("C1R1 decisions are not eighteen complete-state channels")
    classifications: list[str] = []
    admissions: list[bool] = []
    for expected, item in zip(TDG6_COMPLETE_STATE_CHANNELS, decisions, strict=True):
        row = _mapping(item, "channel decision")
        if row.get("channel") != expected:
            _fail("eighteen-channel order differs")
        classification = row.get("classification")
        if type(classification) is not str:
            _fail(f"{expected} classification is absent")
        passed = row.get("admission_passed")
        if type(passed) is not bool:
            _fail(f"{expected} admission flag is absent")
        if classification not in TDG6_PASSING_CLASSES or passed is not True:
            _fail(f"{expected} is not independently admitted")
        classifications.append(classification)
        admissions.append(passed)
    record = {
        "eighteen_channel_classifications": classifications,
        "eighteen_channel_admission": admissions,
        "independent_admission_passed": True,
        "runner_label_sufficient": False,
    }
    _refuse(record)
    return record


def _plan_mapping(plan: object) -> dict[str, object] | None:
    if plan is None:
        return None
    if isinstance(plan, Mapping):
        return dict(plan)
    if hasattr(plan, "requested_cap"):
        return encode_tdg7_plan(plan)
    _fail("attempted or next plan cannot be encoded")
    raise AssertionError


def retain_bridge_result(result: object) -> dict[str, object]:
    if getattr(result, "accepted_state_advanced", False) is True:
        _fail("bridge result advanced accepted state")
    if getattr(result, "committed", None) is not None:
        _fail("bridge result serialized a committed endpoint")
    disposition = getattr(result, "disposition", None)
    if type(disposition) is not str:
        _fail("bridge result omitted disposition")
    evidence = _mapping(getattr(result, "evidence", {}) or {}, "bridge evidence")
    if not hasattr(result, "sink_writes"):
        _fail("bridge result omitted sink_writes")
    sink_writes = result.sink_writes
    if type(sink_writes) is not int or sink_writes != 0:
        _fail("bridge sink_writes must be the integer zero")
    record = {
        "disposition": disposition,
        "accepted_state_advanced": False,
        "sink_writes": sink_writes,
        "evidence": dict(evidence),
        "attempted_plan": _plan_mapping(getattr(result, "attempted_plan", None)),
        "next_plan": _plan_mapping(getattr(result, "next_plan", None)),
        "prepared_present": getattr(result, "prepared", None) is not None,
        "committed_present": False,
    }
    for key in FORBIDDEN_ENDPOINT_KEYS:
        if key in record or key in evidence:
            _fail("bridge result serializes a forbidden endpoint")
    _refuse(record)
    canonical_json_bytes(record)
    return record


def isolate_working_checkpoint(checkpoint: object) -> object:
    isolate = getattr(checkpoint, "isolate", None)
    if callable(isolate):
        working = isolate()
        if working is checkpoint:
            _fail("isolation returned the live checkpoint")
        return working
    if isinstance(checkpoint, HLT17InMemoryCheckpoint):
        isolated_member = isolate_hlt17_runtime_member(checkpoint.member)
        working = HLT17InMemoryCheckpoint(
            isolated_member,
            checkpoint.cursor,
            checkpoint.accepted_encoding,
        )
        if working is checkpoint or working.member is checkpoint.member:
            _fail("HLT17 isolation failed to copy the member")
        return working
    _fail("checkpoint cannot be isolated in memory")
    raise AssertionError


def classify_prepare_outcome(
    result: object, prepared_record: Mapping[str, object] | None
) -> str:
    disposition = getattr(result, "disposition", None)
    if disposition in SOURCE_DISPOSITIONS:
        return "source"
    if disposition in CFL_DISPOSITIONS:
        return "cfl"
    if disposition in {"temporal_retry_required", "temporal_retry_exhausted"}:
        return "temporal"
    if disposition == "invalid_premise":
        return "invalid"
    if disposition == "accepted_fine":
        _fail("qualification adopted an endpoint")
    if disposition not in PREPARED_DISPOSITIONS:
        return "inconclusive"
    if not prepared_record:
        return "inconclusive"
    try:
        independent_channel_admission(
            prepared_record.get("eighteen_channel_decisions")
        )
    except HLT17SRCQ1REC1PREF1Error:
        decisions = prepared_record.get("eighteen_channel_decisions")
        if type(decisions) is list and any(
            item.get("classification") == "order_inconclusive"
            for item in decisions
            if isinstance(item, Mapping)
        ):
            return "inconclusive"
        return "temporal"
    if prepared_record.get("admission_passed") is True:
        return "prepared"
    return "temporal"


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
        _fail("isolated checkpoint cannot absorb a source/CFL overlay")
    absorbed = absorb(result)
    if getattr(absorbed, "accepted_state_advanced", False) is True:
        _fail("overlay absorption advanced accepted state")
    return absorbed


def _directory_payload_digest(files: Mapping[str, bytes]) -> str:
    if type(files) is not dict or len(files) != 1:
        _fail("raw directory must contain exactly one leaf")
    inventory = {path: _sha256_hex(payload) for path, payload in files.items()}
    return _sha256_hex(canonical_json_bytes(inventory))


def historical_runs_digest(
    records: Mapping[str, tuple[int, str]] | list[tuple[str, int, str]]
) -> str:
    """SHA-256 of sha256sum lines ``digest`` two spaces ``path`` newline."""

    rows: list[tuple[str, int, str]]
    if isinstance(records, Mapping):
        rows = [
            (path, int(item[0]), _require_sha(item[1], "runs leaf"))
            for path, item in records.items()
        ]
    else:
        rows = list(records)
    lines = "".join(
        f"{digest}  {path}\n" for path, _size, digest in sorted(rows, key=lambda item: item[0])
    )
    return _sha256_hex(lines.encode("ascii"))


def expected_claims() -> dict[str, object]:
    return dict(EXPECTED_CLAIMS)


def expected_compact_payload() -> dict[str, object]:
    payload = {
        "artifact_id": ARTIFACT_ID,
        "schema_version": SCHEMA_VERSION,
        "classification": CLASSIFICATION,
        "authority": {
            "authority_commit": AUTHORITY_COMMIT,
            "authority_parent": AUTHORITY_PARENT,
            "exact_delta_paths": list(ALLOWED_AUTHORITY_DELTA),
            "freeze_config_sha256": FREEZE_CONFIG_SHA256,
            "freeze_owner_sha256": FREEZE_OWNER_SHA256,
            "freeze_make_sha256": FREEZE_MAKE_SHA256,
            "implementation_sha256": IMPLEMENTATION_SHA256,
            "source_closure_sha256": SOURCE_CLOSURE_SHA256,
            "origin_capture_sha256": ORIGIN_CAPTURE_SHA256,
        },
        "attempt2": {
            "sealed_compact_sha256": ATTEMPT2_SEALED_COMPACT_SHA256,
            "metadata_linked_compact_sha256": ATTEMPT2_METADATA_LINKED_COMPACT_SHA256,
            "raw_leaf_sha256": ATTEMPT2_RAW_LEAF_SHA256,
            "raw_directory_payload_sha256": ATTEMPT2_RAW_DIRECTORY_PAYLOAD_SHA256,
            "rk_receipts": dict(ATTEMPT2_RK_RECEIPTS),
            "rk_members_remeasured": False,
        },
        "raw_terminal": {
            "leaf_relative": REC1_RAW_LEAF_RELATIVE,
            "leaf_sha256": REC1_RAW_LEAF_SHA256,
            "leaf_bytes": REC1_RAW_LEAF_BYTES,
            "directory_payload_sha256": REC1_RAW_DIRECTORY_PAYLOAD_SHA256,
        },
        "historical_runs_baseline": {
            "file_count": HISTORICAL_RUNS_FILE_COUNT,
            "byte_count": HISTORICAL_RUNS_BYTE_COUNT,
            "digest": HISTORICAL_RUNS_DIGEST,
            "excluded_namespace": REC1_OUTPUT_NAMESPACE,
        },
        "planck": {
            path: {"byte_count": size, "sha256": digest}
            for path, (size, digest) in PLANCK_FILES.items()
        },
        "caps": dict(FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER),
        "ssp4097": {
            "initial_cap_hex": FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER[
                "SSPRK3-4097"
            ],
            "independent_cfl": dict(EXPECTED_SSP4097_CFL),
            "next_plan_cap_hex": EXPECTED_SSP4097_NEXT_PLAN_CAP_HEX,
            "reconstructed_prepare_disposition": "prepared_overlay",
            "overlay_count": 1,
            "prepare_count": 2,
            "receipt_sha256": EXPECTED_C1R1_RECEIPTS["SSPRK3-4097"],
            "assessment_sha256": EXPECTED_C1R1_ASSESSMENTS["SSPRK3-4097"],
            "family_sha256": EXPECTED_C1R1_FAMILIES["SSPRK3-4097"],
        },
        "ssp8193": {
            "initial_cap_hex": FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER[
                "SSPRK3-8193"
            ],
            "reconstructed_prepare_disposition": "prepared_fresh",
            "overlay_count": 0,
            "prepare_count": 1,
            "receipt_sha256": EXPECTED_C1R1_RECEIPTS["SSPRK3-8193"],
            "assessment_sha256": EXPECTED_C1R1_ASSESSMENTS["SSPRK3-8193"],
            "family_sha256": EXPECTED_C1R1_FAMILIES["SSPRK3-8193"],
        },
        "ssp16385": {
            "initial_cap_hex": FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER[
                "SSPRK3-16385"
            ],
            "reconstructed_prepare_disposition": "prepared_fresh",
            "overlay_count": 0,
            "prepare_count": 1,
            "receipt_sha256": EXPECTED_C1R1_RECEIPTS["SSPRK3-16385"],
            "assessment_sha256": EXPECTED_C1R1_ASSESSMENTS["SSPRK3-16385"],
            "family_sha256": EXPECTED_C1R1_FAMILIES["SSPRK3-16385"],
        },
        "recorded_raw_resources": {
            "total_wall_seconds": RECORDED_RAW_WALL_SECONDS,
            "rss_bytes": RECORDED_RAW_RSS_BYTES,
            "within_freeze_ceilings": True,
        },
        "claims": expected_claims(),
        "nonclaims": list(NONCLAIMS),
        "written": False,
        "endpoint_serialized": False,
        "campaign_store_written": False,
        "pro20_event_execution_authorized": False,
    }
    _refuse(payload)
    return payload


def expected_compact() -> dict[str, object]:
    return {
        "artifact_id": ARTIFACT_ID,
        "artifact_payload": expected_compact_payload(),
    }


def expected_config() -> dict[str, object]:
    return {
        "artifact_id": ARTIFACT_ID,
        "schema_version": SCHEMA_VERSION,
        "classification": CLASSIFICATION,
        "authority_commit": AUTHORITY_COMMIT,
        "authority_parent": AUTHORITY_PARENT,
        "freeze_config_path": FREEZE_CONFIG_PATH,
        "freeze_config_sha256": FREEZE_CONFIG_SHA256,
        "result_path": RESULT_PATH,
        "owner_document": OWNER_DOCUMENT,
        "accepted_state_advanced": False,
        "campaign_store_written": False,
        "endpoint_adopted_or_serialized": False,
        "pro20_event_execution_authorized": False,
        "written": False,
    }


def render_config(config: Mapping[str, object]) -> bytes:
    mapping = _mapping(config, "PREF1 config")
    expected = expected_config()
    if mapping != expected:
        _fail("PREF1 config fields differ")
    lines = []
    for key, value in expected.items():
        if type(value) is bool:
            rendered = "true" if value else "false"
        elif type(value) is int and not isinstance(value, bool):
            rendered = str(value)
        elif type(value) is str:
            rendered = json.dumps(value, ensure_ascii=True)
        else:
            _fail("unsupported PREF1 config scalar")
        lines.append(f"{key} = {rendered}")
    return ("\n".join(lines) + "\n").encode("ascii")


def emit_config_bytes() -> bytes:
    return render_config(expected_config())


def _read_tracked(root: Path, relative: str) -> bytes:
    try:
        return read_regular_file(root, relative)
    except (EvidenceIOError, UnsafePathError, OSError) as error:
        raise HLT17SRCQ1REC1PREF1Error(
            f"tracked PREF1 leaf is absent: {relative}"
        ) from error


def verify_compact(repository: Path) -> dict[str, object]:
    """Validate tracked compact bytes only. No raw, runs, source, Git, or reconstruction."""

    root = Path(repository)
    try:
        config_raw = _read_tracked(root, CONFIG_PATH)
        result_raw = _read_tracked(root, RESULT_PATH)
    except HLT17SRCQ1REC1PREF1Error as error:
        raise HLT17SRCQ1REC1PREF1Error(ABSENT_COMPACT_MESSAGE) from error
    return validate_compact_result(config_raw, result_raw)


def validate_compact_result(config_raw: bytes, result_raw: bytes) -> dict[str, object]:
    if type(config_raw) is not bytes or type(result_raw) is not bytes:
        _fail("compact bytes differ")
    if render_config(expected_config()) != config_raw:
        try:
            parsed = tomllib.loads(config_raw.decode("utf-8"))
        except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
            raise HLT17SRCQ1REC1PREF1Error("PREF1 config is malformed") from error
        if parsed != expected_config():
            _fail("PREF1 config differs")
        if render_config(parsed) != config_raw:
            _fail("PREF1 config bytes differ")
    if (
        COMPACT_RESULT_SHA256 is not None
        and _sha256_hex(result_raw) != COMPACT_RESULT_SHA256
    ):
        _fail("PREF1 compact SHA-256 differs")
    result = _json_mapping(result_raw, "PREF1 compact")
    if canonical_result(result) != result_raw:
        _fail("PREF1 compact bytes are not canonical")
    expected = expected_compact()
    if result.get("artifact_id") != ARTIFACT_ID:
        _fail("PREF1 compact artifact differs")
    payload = _mapping(result.get("artifact_payload"), "PREF1 compact payload")
    if payload.get("classification") != CLASSIFICATION:
        _fail("PREF1 classification differs")
    if payload != expected["artifact_payload"]:
        _fail("PREF1 compact payload differs")
    if payload["claims"] != expected_claims():
        _fail("PREF1 claims differ")
    _require_bool(payload.get("written"), False, "compact written")
    _require_bool(payload.get("endpoint_serialized"), False, "endpoint_serialized")
    _refuse(result)
    return result


def classification_if_established(facts: Mapping[str, object]) -> str:
    mapping = _mapping(facts, "independent facts")
    required = (
        "authority_bound",
        "freeze_bound",
        "attempt2_bound",
        "rec1_raw_bound",
        "historical_runs_bound",
        "planck_bound",
        "pro20_absent",
        "caps_derived",
        "rk_receipts_bound",
        "ssp4097_cfl_then_prepared",
        "ssp8193_fresh",
        "ssp16385_fresh",
        "eighteen_channels_bound",
        "no_state_advance",
        "resources_within_ceilings",
        "endpoints_absent",
    )
    for name in required:
        if mapping.get(name) is not True:
            _fail(f"independent predicate {name} is not established")
    if mapping.get("runner_label_sufficient") is True:
        _fail("runner labels must not be sufficient")
    return CLASSIFICATION


def _git(root: Path, operation: object) -> object:
    try:
        return git_read(root, operation)
    except EvidenceIOError as error:
        raise HLT17SRCQ1REC1PREF1Error("committed image cannot be inspected") from error


def _blob(root: Path, commit: str, path: str) -> bytes:
    raw = _git(root, ReadBlob(commit=commit, path=path))
    if type(raw) is not bytes:
        _fail(f"{path} blob is not bytes")
    return raw


def _stat_identity(metadata: os.stat_result) -> tuple[int, ...]:
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_mode,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
        metadata.st_nlink,
    )


def _lexists(path: Path) -> bool:
    try:
        path.lstat()
    except FileNotFoundError:
        return False
    except OSError as error:
        raise HLT17SRCQ1REC1PREF1Error("namespace path cannot be read") from error
    return True


def _require_relative(relative: str, label: str) -> Path:
    value = Path(relative)
    if value.is_absolute() or any(part in {"", ".", ".."} for part in value.parts):
        _fail(f"{label} is unsafe")
    return value


def _hash_software_file(path: str) -> str:
    if type(path) is not str or not path or "\0" in path or len(path) > 1024:
        _fail("software path is unsafe")
    try:
        resolved = os.path.realpath(path)
        before = os.lstat(resolved)
    except OSError as error:
        raise HLT17SRCQ1REC1PREF1Error("software path cannot be read") from error
    if (
        type(resolved) is not str
        or not resolved
        or "\0" in resolved
        or stat.S_ISLNK(before.st_mode)
        or not stat.S_ISREG(before.st_mode)
        or before.st_size < 0
    ):
        _fail("software path is not a regular file")
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(resolved, flags)
    except OSError as error:
        raise HLT17SRCQ1REC1PREF1Error("software path cannot be opened") from error
    try:
        opened = os.fstat(descriptor)
        if _stat_identity(opened) != _stat_identity(before):
            _fail("software path changed before read")
        digest = sha256()
        remaining = int(before.st_size)
        while remaining > 0:
            chunk = os.read(descriptor, min(remaining, 1 << 20))
            if not chunk:
                break
            digest.update(chunk)
            remaining -= len(chunk)
        after = os.fstat(descriptor)
        path_after = os.lstat(resolved)
        if (
            remaining != 0
            or _stat_identity(after) != _stat_identity(before)
            or _stat_identity(path_after) != _stat_identity(before)
            or not stat.S_ISREG(after.st_mode)
        ):
            _fail("software path changed during read")
        return digest.hexdigest()
    finally:
        os.close(descriptor)


def _runtime_regular_path(path: str) -> str:
    if type(path) is not str or not path or "\0" in path:
        _fail("runtime path is unsafe")
    try:
        resolved = os.path.realpath(path)
        metadata = os.lstat(resolved)
    except OSError as error:
        raise HLT17SRCQ1REC1PREF1Error("runtime path cannot be resolved") from error
    if (
        type(resolved) is not str
        or not resolved
        or "\0" in resolved
        or len(resolved) > 1024
        or stat.S_ISLNK(metadata.st_mode)
        or not stat.S_ISREG(metadata.st_mode)
    ):
        _fail("runtime path is not a regular file")
    return resolved


def environment_identity() -> dict[str, str]:
    observed = observe_environment()
    executable = _runtime_regular_path(sys.executable)
    import numpy._core._multiarray_umath as numpy_extension

    extension_path = getattr(numpy_extension, "__file__", None)
    extension = _runtime_regular_path(str(extension_path))
    record = {
        **observed,
        "executable_realpath": executable,
        "executable_sha256": _hash_software_file(executable),
        "numpy_extension_realpath": extension,
        "numpy_extension_sha256": _hash_software_file(extension),
        "kernel_release": platform.release(),
    }
    if tuple(record) != ENVIRONMENT_KEYS:
        _fail("expanded environment keys differ")
    for key in ENVIRONMENT_KEYS:
        value = record[key]
        if type(value) is not str or not value:
            _fail(f"environment {key} is empty")
    return record


def _tree_python_paths(root: Path, commit: str) -> tuple[str, ...]:
    tree = _git(root, InspectTree(commit=commit, path=SOURCE_TREE))
    paths: list[str] = []
    seen: set[str] = set()
    for entry in tree:
        path = getattr(entry, "path", None)
        kind = getattr(entry, "kind", None)
        mode = getattr(entry, "mode", None)
        if type(path) is not str or type(kind) is not str or type(mode) is not str:
            _fail("source closure tree entry differs")
        if not path.endswith((".py", ".pyi")):
            continue
        if kind != "blob" or mode not in _REGULAR_GIT_MODES:
            _fail("source closure contains a nonregular Python path")
        pure = PurePosixPath(path)
        if pure.is_absolute() or ".." in pure.parts:
            _fail("source closure path is unsafe")
        if path in seen:
            _fail("source closure repeats a path")
        seen.add(path)
        paths.append(path)
    if REC1_RUNNER_PATH not in seen:
        paths.append(REC1_RUNNER_PATH)
        seen.add(REC1_RUNNER_PATH)
    return tuple(sorted(paths))


def source_closure_identity(
    root: Path,
    implementation_commit: str,
    *,
    environment: Mapping[str, str],
) -> dict[str, object]:
    repository = _root(root)
    commit = _require_commit(implementation_commit, "implementation commit")
    records: list[tuple[str, str]] = []
    digest = sha256(SOURCE_CLOSURE_DOMAIN)
    for path in _tree_python_paths(repository, commit):
        raw = _blob(repository, commit, path)
        value = _sha256_hex(raw)
        records.append((path, value))
        digest.update(f"{path}\0{value}\n".encode("ascii"))
    observed = dict(environment)
    if tuple(observed) != ENVIRONMENT_KEYS:
        _fail("source-closure environment keys differ")
    for name in ("executable_sha256", "numpy_extension_sha256"):
        digest.update(
            f"environment:{name}\0{_require_sha(observed[name], name)}\n".encode(
                "ascii"
            )
        )
    return {
        "sha256": digest.hexdigest(),
        "file_count": len(records),
        "files": tuple(records),
    }


def authenticate_authority(root: Path) -> dict[str, object]:
    repository = _root(root)
    authority = _require_commit(AUTHORITY_COMMIT, "authority commit")
    parents = _git(repository, CommitParents(authority))
    if type(parents) is not tuple or len(parents) != 1:
        _fail("REC1 authority must have one implementation parent")
    parent = parents[0]
    if parent != AUTHORITY_PARENT:
        _fail("REC1 authority parent is not A2")
    delta = _git(repository, InspectDelta(commit=authority, against=parent))
    paths = tuple(sorted(item.path for item in delta))
    if paths != tuple(sorted(ALLOWED_AUTHORITY_DELTA)):
        _fail("REC1 authority delta paths differ")
    freeze_hashes = {
        FREEZE_CONFIG_PATH: FREEZE_CONFIG_SHA256,
        FREEZE_OWNER_PATH: FREEZE_OWNER_SHA256,
        FREEZE_MAKE_PATH: FREEZE_MAKE_SHA256,
    }
    for path, digest in freeze_hashes.items():
        committed = _blob(repository, authority, path)
        if _sha256_hex(committed) != digest:
            _fail(f"authority freeze blob differs: {path}")
        try:
            live = read_regular_file(repository, path)
        except EvidenceIOError as error:
            raise HLT17SRCQ1REC1PREF1Error(
                f"live freeze path cannot be read: {path}"
            ) from error
        if live != committed:
            _fail(f"live freeze path differs from authority commit: {path}")
    implementation = _sha256_hex(_blob(repository, parent, REC1_HARNESS_PATH))
    runner = _sha256_hex(_blob(repository, parent, REC1_RUNNER_PATH))
    rec1_authority = _sha256_hex(_blob(repository, parent, REC1_AUTHORITY_PATH))
    if (
        implementation != IMPLEMENTATION_SHA256
        or runner != RUNNER_SHA256
        or rec1_authority != REC1_AUTHORITY_SHA256
    ):
        _fail("A2 implementation/runner/authority hash differs")
    return {
        "authority_commit": authority,
        "authority_parent": parent,
        "exact_delta_paths": list(ALLOWED_AUTHORITY_DELTA),
        "implementation_sha256": implementation,
        "runner_sha256": runner,
        "authority_sha256": rec1_authority,
    }


def authenticate_freeze_config(
    root: Path, *, environment: Mapping[str, str]
) -> dict[str, object]:
    repository = _root(root)
    raw = read_regular_file(repository, FREEZE_CONFIG_PATH)
    if _sha256_hex(raw) != FREEZE_CONFIG_SHA256:
        _fail("live freeze config SHA-256 differs")
    try:
        config = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise HLT17SRCQ1REC1PREF1Error("freeze config is malformed") from error
    mapping = _mapping(config, "freeze config")
    if mapping.get("artifact_id") != FREEZE_ARTIFACT_ID:
        _fail("freeze artifact differs")
    if mapping.get("implementation_commit") != AUTHORITY_PARENT:
        _fail("freeze implementation parent differs")
    if mapping.get("source_closure_sha256") != SOURCE_CLOSURE_SHA256:
        _fail("freeze source closure differs")
    if mapping.get("origin_capture_sha256") != ORIGIN_CAPTURE_SHA256:
        _fail("freeze origin capture differs")
    observed = dict(environment)
    if mapping.get("environment") != observed:
        _fail("freeze environment differs")
    identity = hlt17_implementation_identity()
    if mapping.get("implementation_identity") != identity:
        _fail("freeze C1R1 implementation identity differs")
    _require_bool(mapping.get("accepted_state_advanced"), False, "freeze state")
    _require_bool(mapping.get("campaign_store_written"), False, "freeze store")
    _require_bool(
        mapping.get("endpoint_adopted_or_serialized"), False, "freeze endpoint"
    )
    _require_bool(mapping.get("pro20_namespace_absent"), True, "freeze PRO20")
    _require_bool(mapping.get("campaign_execution_authorized"), False, "freeze campaign")
    return mapping


def authenticate_attempt2_compacts(root: Path) -> dict[str, object]:
    repository = _root(root)
    sealed = _blob(
        repository, ATTEMPT2_SEALED_COMPACT_COMMIT, ATTEMPT2_COMPACT_RELATIVE
    )
    linked = _blob(
        repository,
        ATTEMPT2_METADATA_LINKED_COMPACT_COMMIT,
        ATTEMPT2_COMPACT_RELATIVE,
    )
    try:
        live = read_regular_file(repository, ATTEMPT2_COMPACT_RELATIVE)
    except EvidenceIOError as error:
        raise HLT17SRCQ1REC1PREF1Error("live Attempt2 compact cannot be read") from error
    if _sha256_hex(sealed) != ATTEMPT2_SEALED_COMPACT_SHA256:
        _fail("sealed Attempt2 compact SHA-256 differs")
    if _sha256_hex(linked) != ATTEMPT2_METADATA_LINKED_COMPACT_SHA256:
        _fail("metadata-linked Attempt2 compact SHA-256 differs")
    if live != linked:
        _fail("live Attempt2 compact is not the metadata-linked blob")
    sealed_obj = _json_mapping(sealed, "sealed Attempt2 compact")
    linked_obj = _json_mapping(linked, "metadata-linked Attempt2 compact")
    if "derivation_document" in sealed_obj:
        _fail("sealed Attempt2 compact must not carry derivation_document")
    expected = dict(sealed_obj)
    expected["derivation_document"] = ATTEMPT2_DERIVATION_DOCUMENT
    if linked_obj != expected:
        _fail("metadata-linked compact is not sealed payload plus derivation_document")
    if linked_obj.get("artifact_id") != ATTEMPT2_ARTIFACT_ID:
        _fail("Attempt2 compact artifact differs")
    if linked_obj.get("classification") != ATTEMPT2_CLASSIFICATION:
        _fail("Attempt2 classification is a sealed predecessor claim, not a PREF1 class")
    if linked_obj.get("authority_commit") != ATTEMPT2_AUTHORITY_COMMIT:
        _fail("Attempt2 authority commit differs")
    if linked_obj.get("implementation_commit") != ATTEMPT2_IMPLEMENTATION_COMMIT:
        _fail("Attempt2 implementation commit differs")
    if linked_obj.get("implementation_sha256") != ATTEMPT2_IMPLEMENTATION_SHA256:
        _fail("Attempt2 implementation SHA-256 differs")
    state = _mapping(linked_obj.get("state_boundary"), "Attempt2 state_boundary")
    _require_bool(state.get("accepted_state_advanced"), False, "Attempt2 state")
    _require_bool(state.get("campaign_store_written"), False, "Attempt2 store")
    _require_bool(
        state.get("endpoint_adopted_or_serialized"), False, "Attempt2 endpoint"
    )
    _require_bool(state.get("pro20_namespace_absent"), True, "Attempt2 PRO20")
    return {
        "sealed": sealed_obj,
        "linked": linked_obj,
        "compact_bytes": live,
    }


def _read_single_leaf_directory(
    root: Path,
    directory_relative: str,
    *,
    leaf_name: str,
    expected_bytes: int,
    expected_leaf_sha256: str,
    expected_directory_sha256: str,
    label: str,
) -> tuple[bytes, dict[str, Any]]:
    repository = _root(root)
    relative = _require_relative(directory_relative, f"{label} directory")
    directory = repository / relative
    try:
        metadata = directory.lstat()
    except OSError as error:
        raise HLT17SRCQ1REC1PREF1Error(f"{label} directory is absent") from error
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        _fail(f"{label} directory is unsafe")
    flags = (
        os.O_RDONLY
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    try:
        descriptor = os.open(directory, flags)
    except OSError as error:
        raise HLT17SRCQ1REC1PREF1Error(f"{label} directory cannot be opened") from error
    try:
        opened = os.fstat(descriptor)
        if _stat_identity(opened) != _stat_identity(metadata):
            _fail(f"{label} directory changed before read")
        names = tuple(os.listdir(descriptor))
        if names != (leaf_name,):
            _fail(f"{label} directory leaf count or name differs")
        leaf_flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
        try:
            leaf = os.open(leaf_name, leaf_flags, dir_fd=descriptor)
        except OSError as error:
            raise HLT17SRCQ1REC1PREF1Error(f"{label} leaf cannot be opened") from error
        try:
            leaf_before = os.fstat(leaf)
            if (
                not stat.S_ISREG(leaf_before.st_mode)
                or leaf_before.st_size != expected_bytes
            ):
                _fail(f"{label} leaf is not the expected regular file")
            chunks: list[bytes] = []
            remaining = expected_bytes
            while remaining > 0:
                chunk = os.read(leaf, min(remaining, 1 << 20))
                if not chunk:
                    break
                chunks.append(chunk)
                remaining -= len(chunk)
            raw = b"".join(chunks)
            leaf_after = os.fstat(leaf)
            if remaining != 0 or _stat_identity(leaf_after) != _stat_identity(
                leaf_before
            ):
                _fail(f"{label} leaf changed during read")
        finally:
            os.close(leaf)
        if tuple(os.listdir(descriptor)) != (leaf_name,):
            _fail(f"{label} directory changed during read")
        directory_after = os.fstat(descriptor)
        path_after = directory.lstat()
        if (
            _stat_identity(directory_after) != _stat_identity(opened)
            or _stat_identity(path_after) != _stat_identity(metadata)
        ):
            _fail(f"{label} directory changed during read")
    finally:
        os.close(descriptor)
    if len(raw) != expected_bytes or _sha256_hex(raw) != expected_leaf_sha256:
        _fail(f"{label} leaf SHA-256 differs")
    digest = _directory_payload_digest({leaf_name: raw})
    if digest != expected_directory_sha256:
        _fail(f"{label} directory payload SHA-256 differs")
    payload = _json_mapping(raw, f"{label} terminal")
    _require_bool(payload.get("accepted_state_advanced"), False, f"{label} state")
    _require_bool(payload.get("endpoint_adopted"), False, f"{label} endpoint")
    _require_bool(payload.get("store_published"), False, f"{label} store")
    _require_bool(
        payload.get("campaign_execution_authorized"), False, f"{label} campaign"
    )
    _refuse(payload)
    return raw, payload


def authenticate_attempt2_raw(root: Path) -> dict[str, object]:
    raw, payload = _read_single_leaf_directory(
        root,
        ATTEMPT2_RAW_DIRECTORY,
        leaf_name="qualification.json",
        expected_bytes=ATTEMPT2_RAW_LEAF_BYTES,
        expected_leaf_sha256=ATTEMPT2_RAW_LEAF_SHA256,
        expected_directory_sha256=ATTEMPT2_RAW_DIRECTORY_PAYLOAD_SHA256,
        label="Attempt2 raw",
    )
    members = payload.get("members")
    if type(members) is not list or len(members) != 4:
        _fail("Attempt2 raw member inventory differs")
    receipts: dict[str, str] = {}
    for key, item in zip(PREDECESSOR_MEMBER_KEYS, members[:3], strict=True):
        row = _mapping(item, f"Attempt2 {key}")
        if row.get("member_key") != key:
            _fail(f"Attempt2 raw {key} order differs")
        c1r1 = _mapping(row.get("c1r1"), f"Attempt2 {key} C1R1")
        receipt = _require_sha(c1r1.get("receipt_sha256"), f"{key} receipt")
        if receipt != ATTEMPT2_RK_RECEIPTS[key]:
            _fail(f"Attempt2 {key} receipt differs")
        receipts[key] = receipt
    ssp = _mapping(members[3], "Attempt2 raw SSPRK3-4097")
    if ssp.get("member_key") != "SSPRK3-4097":
        _fail("Attempt2 raw resume member differs")
    if ssp.get("c1r1") is not None:
        _fail("Attempt2 raw SSPRK3-4097 carried C1R1")
    return {"raw_bytes": raw, "payload": payload, "rk_receipts": receipts}


def authenticate_rec1_raw(root: Path) -> dict[str, object]:
    raw, payload = _read_single_leaf_directory(
        root,
        REC1_OUTPUT_NAMESPACE,
        leaf_name="qualification.json",
        expected_bytes=REC1_RAW_LEAF_BYTES,
        expected_leaf_sha256=REC1_RAW_LEAF_SHA256,
        expected_directory_sha256=REC1_RAW_DIRECTORY_PAYLOAD_SHA256,
        label="REC1 raw",
    )
    if "halted_outcome" not in payload or payload.get("halted_outcome") is not None:
        _fail("REC1 raw halt claim differs from the no-halt terminal")
    if payload.get("rk_members_remeasured") is not False:
        _fail("REC1 raw reports RK remeasurement")
    members = payload.get("members")
    if type(members) is not list or len(members) != 6:
        _fail("REC1 raw member inventory differs")
    resources = _mapping(payload.get("resource_facts"), "REC1 raw resources")
    if resources.get("total_wall_seconds") != RECORDED_RAW_WALL_SECONDS:
        _fail("REC1 raw wall claim differs")
    if resources.get("rss_bytes") != RECORDED_RAW_RSS_BYTES:
        _fail("REC1 raw RSS claim differs")
    return {"raw_bytes": raw, "payload": payload}


def require_pro20_namespace_absent(root: Path) -> None:
    repository = _root(root)
    relative = _require_relative(PRO20_EVENT_NAMESPACE, "PRO20 event namespace")
    target = repository / relative
    if _lexists(target):
        _fail("PRO20 event namespace already exists")
    current = repository
    for part in relative.parts:
        current = current / part
        if not _lexists(current):
            return
        try:
            metadata = current.lstat()
        except OSError as error:
            raise HLT17SRCQ1REC1PREF1Error("PRO20 namespace cannot be read") from error
        if stat.S_ISLNK(metadata.st_mode):
            _fail("PRO20 event namespace is unsafe")


def authenticate_planck(root: Path) -> dict[str, object]:
    repository = _root(root)
    records: dict[str, object] = {}
    for relative, (size, digest) in PLANCK_FILES.items():
        path = repository / relative
        try:
            metadata = path.lstat()
        except FileNotFoundError:
            records[relative] = {
                "present": False,
                "byte_count": size,
                "sha256": digest,
            }
            continue
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
            _fail(f"Planck path is unsafe: {relative}")
        if metadata.st_size != size:
            _fail(f"Planck byte count differs: {relative}")
        if _hash_software_file(str(path)) != digest:
            _fail(f"Planck SHA-256 differs: {relative}")
        records[relative] = {
            "present": True,
            "byte_count": size,
            "sha256": digest,
        }
    return records


def inventory_historical_runs(root: Path) -> dict[str, object]:
    repository = _root(root)
    relative = _require_relative(HISTORICAL_RUNS_ROOT, "runs root")
    base = repository / relative
    try:
        info = base.lstat()
    except OSError as error:
        raise HLT17SRCQ1REC1PREF1Error("historical runs tree is absent") from error
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
        _fail("historical runs tree is unsafe")
    pending = [str(relative)]
    files: list[str] = []
    excluded_prefix = REC1_OUTPUT_NAMESPACE + "/"
    while pending:
        current_relative = pending.pop()
        current = repository / current_relative
        try:
            with os.scandir(current) as iterator:
                children = sorted(iterator, key=lambda item: item.name)
        except OSError as error:
            raise HLT17SRCQ1REC1PREF1Error(
                "historical runs tree cannot be scanned"
            ) from error
        for child in children:
            child_relative = f"{current_relative}/{child.name}"
            try:
                child_info = child.stat(follow_symlinks=False)
            except OSError as error:
                raise HLT17SRCQ1REC1PREF1Error(
                    f"historical runs path is unsafe: {child_relative}"
                ) from error
            if stat.S_ISLNK(child_info.st_mode):
                _fail(f"historical runs path is unsafe: {child_relative}")
            if stat.S_ISDIR(child_info.st_mode):
                pending.append(child_relative)
                continue
            if not stat.S_ISREG(child_info.st_mode) or child_info.st_nlink != 1:
                _fail(f"historical runs path is unsafe: {child_relative}")
            if (
                child_relative == REC1_RAW_LEAF_RELATIVE
                or child_relative.startswith(excluded_prefix)
            ):
                continue
            files.append(child_relative)
    files = sorted(files)
    rows: list[tuple[str, int, str]] = []
    total = 0
    for path in files:
        absolute = repository / path
        try:
            metadata = absolute.lstat()
        except OSError as error:
            raise HLT17SRCQ1REC1PREF1Error(
                f"historical runs path vanished: {path}"
            ) from error
        if (
            stat.S_ISLNK(metadata.st_mode)
            or not stat.S_ISREG(metadata.st_mode)
            or metadata.st_nlink != 1
        ):
            _fail(f"historical runs path is unsafe: {path}")
        digest = _hash_software_file(str(absolute))
        rows.append((path, metadata.st_size, digest))
        total += metadata.st_size
    digest = historical_runs_digest(rows)
    if (
        len(rows) != HISTORICAL_RUNS_FILE_COUNT
        or total != HISTORICAL_RUNS_BYTE_COUNT
        or digest != HISTORICAL_RUNS_DIGEST
    ):
        _fail("historical runs baseline differs")
    return {
        "file_count": len(rows),
        "byte_count": total,
        "digest": digest,
        "excluded_namespace": REC1_OUTPUT_NAMESPACE,
    }


def independent_six_caps(
    origin: object, runtime: object
) -> tuple[dict[str, float], dict[str, dict[str, str]]]:
    members = getattr(runtime, "members", None)
    captured_members = getattr(origin, "members", None)
    if not isinstance(members, Mapping) or tuple(members) != MEMBER_KEYS:
        _fail("runtime origin member order differs for cap derivation")
    if not isinstance(captured_members, (tuple, list)):
        _fail("captured origin members differ for cap derivation")
    captured = {item.member_key: item for item in captured_members}
    if tuple(captured) != MEMBER_KEYS:
        _fail("captured origin member order differs for cap derivation")
    caps: dict[str, float] = {}
    facts: dict[str, dict[str, str]] = {}
    for key in MEMBER_KEYS:
        try:
            causal = load_canonical_json(captured[key].causal_bytes)
        except (CanonicalJSONError, TypeError, ValueError) as error:
            raise HLT17SRCQ1REC1PREF1Error(
                f"{key} retained causal bytes differ"
            ) from error
        causal_map = _mapping(causal, f"{key} retained causal state")
        speed = causal_map.get("previous_speed_upper")
        if isinstance(speed, bool) or not isinstance(speed, (int, float)):
            _fail(f"{key} inherited causal speed differs")
        previous_speed = float(speed)
        transaction = members[key].member.transaction
        spacing = float(transaction.grid_spacing)
        cfl = float(transaction.cfl_maximum)
        if not all(
            math.isfinite(item) and item > 0.0
            for item in (previous_speed, spacing, cfl)
        ):
            _fail(f"{key} cap inputs are outside their positive finite domain")
        cap = cfl * spacing / previous_speed
        expected = FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER[key]
        if not math.isfinite(cap) or cap <= 0.0 or cap.hex() != expected:
            _fail(f"{key} derived requested cap differs")
        caps[key] = cap
        facts[key] = {
            "formula": "cfl_maximum*grid_spacing/inherited_previous_speed_upper",
            "cfl_maximum_hex": cfl.hex(),
            "grid_spacing_hex": spacing.hex(),
            "inherited_previous_speed_upper_hex": previous_speed.hex(),
            "requested_cap_hex": cap.hex(),
        }
    return caps, facts


def derive_six_caps(origin: object, runtime: object) -> dict[str, dict[str, str]]:
    caps, facts = independent_six_caps(origin, runtime)
    for key, _spacing, _speed, requested in MEMBER_CAP_INPUTS:
        derived = caps[key]
        spacing = float.fromhex(facts[key]["grid_spacing_hex"])
        speed = float.fromhex(facts[key]["inherited_previous_speed_upper_hex"])
        cfl = float.fromhex(facts[key]["cfl_maximum_hex"])
        if (cfl * spacing / speed).hex() != derived.hex():
            _fail(f"{key} cap formula differs")
        if derived.hex() != requested:
            _fail(f"{key} derived cap is not the frozen SRCQ1 hex")
        if facts[key]["cfl_maximum_hex"] != CFL_MAXIMUM_HEX:
            _fail(f"{key} CFL maximum differs")
    return facts


def bind_attempt2_rk_receipts(
    attempt2_raw: Mapping[str, object], rec1_raw: Mapping[str, object]
) -> list[dict[str, object]]:
    raw_members = attempt2_raw.get("members")
    rec1_members = rec1_raw.get("members")
    if type(raw_members) is not list or type(rec1_members) is not list:
        _fail("RK predecessor inventories differ")
    records: list[dict[str, object]] = []
    rec1_by_key = {}
    for item in rec1_members:
        row = _mapping(item, "REC1 member")
        key = row.get("member_key")
        if type(key) is str:
            rec1_by_key[key] = row
    for key, item in zip(PREDECESSOR_MEMBER_KEYS, raw_members[:3], strict=True):
        row = _mapping(item, f"Attempt2 {key}")
        c1r1 = _mapping(row.get("c1r1"), f"Attempt2 {key} C1R1")
        receipt = _require_sha(c1r1.get("receipt_sha256"), f"{key} receipt")
        if receipt != ATTEMPT2_RK_RECEIPTS[key]:
            _fail(f"Attempt2 {key} receipt differs")
        rec1_row = rec1_by_key.get(key)
        if rec1_row is None:
            _fail(f"REC1 omitted predecessor {key}")
        if rec1_row.get("remeasured") is True:
            _fail(f"REC1 remeasured predecessor {key}")
        if rec1_row.get("attempt2_c1r1_receipt_sha256") != receipt:
            _fail(f"REC1 {key} receipt reference differs from Attempt2 raw")
        record = {
            "member_key": key,
            "role": "predecessor_evidence",
            "remeasured": False,
            "attempt2_c1r1_receipt_sha256": receipt,
            "independent_receipt_bound": True,
            "accepted_state_advanced": False,
            "endpoint_adopted": False,
        }
        _refuse(record)
        records.append(record)
    return records


def reproduce_ssp_member(
    *,
    member_key: str,
    view: object,
    captured: object,
    initial_cap: float,
    evaluate: Callable[[object], Mapping[str, object]],
    expected_config: str | None,
) -> dict[str, object]:
    member = view.member
    binding = view.source_binding
    check_method_owned_identity(member_key, member)
    check_inherited_counters(member_key, member, captured)
    configuration = _require_sha(
        getattr(binding, "configuration_sha256", None),
        f"{member_key} source configuration",
    )
    if expected_config is not None and expected_config != configuration:
        _fail(f"foreign source configuration: {member_key}")
    if initial_cap.hex() != FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER[member_key]:
        _fail(f"{member_key} initial cap is not frozen SRCQ1 hex")
    before = member_fingerprints(member, binding)
    source_record = dict(evaluate(view))
    live_checkpoint = view.checkpoint
    if hasattr(live_checkpoint, "agree"):
        live_checkpoint.agree()
    working = isolate_working_checkpoint(live_checkpoint)
    if hasattr(working, "agree"):
        working.agree()
    attempts: list[dict[str, object]] = []
    expected_plan: dict[str, object] | None = None
    c1r1_record = None
    last_disposition = None
    independent_cfl = None
    overlay_count = 0
    source_owner_count = 0
    cfl_owner_count = 0
    prepared_index = None
    for prepare_index in range(MAX_PREPARE_ATTEMPTS):
        try:
            prepared_result = working.prepare(requested_cap=initial_cap)
        except C1R1RuntimeResourceStop as error:
            raise HLT17SRCQ1REC1PREF1ResourceError(str(error)) from error
        except C1R1CoordinateLatticeStop as error:
            raise HLT17SRCQ1REC1PREF1Error(str(error)) from error
        except C1R1RefinementPathStop as error:
            cause = getattr(error, "cause", None)
            if isinstance(cause, CFLRetryRequired):
                independent_cfl_ratio(
                    {
                        "kind": "cfl_retry",
                        "observed_ratio_hex": float(cause.observed_ratio).hex(),
                        "maximum_ratio_hex": float(cause.maximum_ratio).hex(),
                    }
                )
            raise HLT17SRCQ1REC1PREF1Error(str(error)) from error
        if getattr(prepared_result, "accepted_state_advanced", False) is True:
            _fail("C1R1 shadow prepare advanced accepted state")
        if getattr(prepared_result, "committed", None) is not None:
            _fail("C1R1 shadow prepare carried a committed runtime")
        retained = retain_bridge_result(prepared_result)
        last_disposition = retained["disposition"]
        if expected_plan is not None:
            if retained["attempted_plan"] != expected_plan:
                _fail(f"{member_key} prepare did not follow overlay.next_plan")
        elif prepare_index == 0:
            attempted = retained["attempted_plan"]
            if (
                isinstance(attempted, Mapping)
                and attempted.get("requested_cap_hex")
                != FROZEN_PROSPECTIVE_REQUESTED_CAP_HEX_BY_MEMBER[member_key]
            ):
                _fail(f"{member_key} first prepare is not the frozen initial cap")
        prepared = getattr(prepared_result, "prepared", None)
        member_c1r1 = (
            independent_c1r1_record(prepared) if prepared is not None else None
        )
        classified = classify_prepare_outcome(prepared_result, member_c1r1)
        attempt_record = {
            "prepare_index": prepare_index,
            "requested_initial_cap_hex": initial_cap.hex(),
            "bridge": retained,
            "independent_class": classified,
            "overlay_absorbed": False,
            "source_owner_count": source_owner_count,
            "cfl_owner_count": cfl_owner_count,
        }
        if last_disposition in CFL_DISPOSITIONS:
            independent_cfl = independent_cfl_ratio(retained["evidence"])
            retained["independent_cfl"] = independent_cfl
            cfl_owner_count += 1
            attempt_record["cfl_owner_count"] = cfl_owner_count
        if last_disposition in SOURCE_DISPOSITIONS:
            source_evidence = _mapping(retained["evidence"], "source overlay")
            if source_evidence.get("kind") != "source_retry":
                _fail("source overlay omitted source_retry kind")
            source_owner_count += 1
            attempt_record["source_owner_count"] = source_owner_count
        _refuse(attempt_record)
        attempts.append(attempt_record)
        if classified == "prepared":
            c1r1_record = member_c1r1
            prepared_index = prepare_index
            break
        if last_disposition in OVERLAY_CONTINUE or last_disposition in {
            CFL_EXHAUSTED,
            SOURCE_EXHAUSTED,
        }:
            _absorb_overlay(working, prepared_result)
            attempt_record["overlay_absorbed"] = True
            overlay_count += 1
            if last_disposition in {CFL_EXHAUSTED, SOURCE_EXHAUSTED} or _overlay_exhausted(
                prepared_result
            ):
                break
            next_plan = retained["next_plan"]
            if not isinstance(next_plan, Mapping):
                _fail(f"{member_key} overlay omitted next_plan")
            expected_plan = next_plan
            continue
        break
    after = member_fingerprints(member, binding)
    if after != before:
        _fail(f"{member_key} live fingerprints changed after isolated overlay work")
    if c1r1_record is None:
        _fail(f"{member_key} did not independently reconstruct a C1R1 prepared path")
    admission = independent_channel_admission(
        c1r1_record.get("eighteen_channel_decisions")
    )
    record = {
        "member_key": member_key,
        "role": "newly_owned",
        "remeasured": False,
        "source_evaluated_once": True,
        "source_evaluation": source_record,
        "c1r1": c1r1_record,
        "reconstructed_prepare_disposition": last_disposition,
        "independent_admission": admission,
        "independent_cfl": independent_cfl,
        "attempts": attempts,
        "prepare_count": len(attempts),
        "overlay_count": overlay_count,
        "source_owner_count": source_owner_count,
        "cfl_owner_count": cfl_owner_count,
        "prepared_index": prepared_index,
        "next_plan_cap_hex": (
            None
            if expected_plan is None
            else expected_plan.get("requested_cap_hex")
        ),
        "accepted_state_advanced": False,
        "endpoint_adopted": False,
        "rk_remeasured": False,
    }
    _refuse(record)
    return record


def _captured_member(origin: object, member_key: str):
    members = getattr(origin, "members", ())
    matches = [
        item for item in members if getattr(item, "member_key", None) == member_key
    ]
    if len(matches) != 1:
        _fail(f"captured origin member differs: {member_key}")
    return matches[0]


def compare_ssp_to_raw(
    reconstructed: Mapping[str, object], raw_member: Mapping[str, object]
) -> None:
    key = reconstructed["member_key"]
    raw = _mapping(raw_member, f"raw {key}")
    if raw.get("member_key") != key:
        _fail(f"raw {key} identity differs")
    raw_c1r1 = _mapping(raw.get("c1r1"), f"raw {key} C1R1")
    rec_c1r1 = _mapping(reconstructed.get("c1r1"), f"reconstructed {key} C1R1")
    for name in ("receipt_sha256", "assessment_sha256", "family_sha256"):
        if rec_c1r1.get(name) != raw_c1r1.get(name):
            _fail(f"{key} independent {name} differs from raw")
        if rec_c1r1.get(name) != {
            "receipt_sha256": EXPECTED_C1R1_RECEIPTS,
            "assessment_sha256": EXPECTED_C1R1_ASSESSMENTS,
            "family_sha256": EXPECTED_C1R1_FAMILIES,
        }[name][key]:
            _fail(f"{key} independent {name} is not the bound identity")
    raw_decisions = raw_c1r1.get("eighteen_channel_decisions")
    rec_decisions = rec_c1r1.get("eighteen_channel_decisions")
    if rec_decisions != raw_decisions:
        _fail(f"{key} independent eighteen-channel facts differ from raw")
    independent_channel_admission(rec_decisions)
    if reconstructed.get("reconstructed_prepare_disposition") != raw.get(
        "prepare_disposition"
    ):
        _fail(
            f"{key} reconstructed disposition differs from raw; "
            "the raw label is compared as a fact, not trusted as classification"
        )
    if reconstructed.get("prepare_count") != raw.get("prepare_count"):
        _fail(f"{key} prepare count differs")
    if reconstructed.get("overlay_count") != raw.get("overlay_count"):
        _fail(f"{key} overlay count differs")
    if key == "SSPRK3-4097":
        if reconstructed.get("independent_cfl") != dict(EXPECTED_SSP4097_CFL):
            _fail("SSPRK3-4097 independent CFL ratio differs")
        if reconstructed.get("next_plan_cap_hex") != EXPECTED_SSP4097_NEXT_PLAN_CAP_HEX:
            _fail("SSPRK3-4097 did not follow the bound next_plan cap")
        if reconstructed.get("overlay_count") != 1:
            _fail("SSPRK3-4097 overlay count is not one")
        if reconstructed.get("reconstructed_prepare_disposition") != "prepared_overlay":
            _fail("SSPRK3-4097 reconstructed disposition is not prepared_overlay")
    else:
        if reconstructed.get("overlay_count") != 0:
            _fail(f"{key} was not a fresh reconstruction")
        if reconstructed.get("reconstructed_prepare_disposition") != "prepared_fresh":
            _fail(f"{key} reconstructed disposition is not prepared_fresh")
        if reconstructed.get("independent_cfl") is not None:
            _fail(f"{key} unexpectedly reconstructed CFL evidence")


def _check_resources(ceilings: ResourceCeilings, wall: float, rss: int) -> None:
    if type(wall) is not float or not math.isfinite(wall) or wall < 0.0:
        raise HLT17SRCQ1REC1PREF1ResourceError("wall observation is not finite")
    if type(rss) is not int or rss < 0:
        raise HLT17SRCQ1REC1PREF1ResourceError("RSS observation is negative")
    if rss > ceilings.max_rss_bytes:
        raise HLT17SRCQ1REC1PREF1ResourceError("RSS exceeded the freeze ceiling")
    if wall > ceilings.max_total_wall_seconds:
        raise HLT17SRCQ1REC1PREF1ResourceError("wall exceeded the freeze ceiling")
    if RECORDED_RAW_RSS_BYTES > ceilings.max_rss_bytes:
        raise HLT17SRCQ1REC1PREF1ResourceError("recorded raw RSS exceeds the freeze ceiling")
    if RECORDED_RAW_WALL_SECONDS > ceilings.max_total_wall_seconds:
        raise HLT17SRCQ1REC1PREF1ResourceError(
            "recorded raw wall exceeds the freeze ceiling"
        )


def reproduce_three_ssp_paths(
    origin: object,
    runtime: object,
    *,
    source_evaluator: Callable[..., Mapping[str, object]] | None = None,
) -> list[dict[str, object]]:
    members = getattr(runtime, "members", None)
    if not isinstance(members, Mapping) or tuple(members) != MEMBER_KEYS:
        _fail("runtime origin member order differs")
    evaluate = source_evaluator or (
        lambda view: independent_source_evaluation(
            binding=view.source_binding,
            time=float(view.member.time),
            state=view.member.state,
            transaction=view.member.transaction,
            tracers=view.member.tracers,
        )
    )
    caps, _facts = independent_six_caps(origin, runtime)
    newly: list[dict[str, object]] = []
    for key in RESUME_MEMBER_KEYS:
        view = members[key]
        captured = _captured_member(origin, key)
        record = reproduce_ssp_member(
            member_key=key,
            view=view,
            captured=captured,
            initial_cap=caps[key],
            evaluate=evaluate,
            expected_config=SOURCE_CONFIGURATION_SHA256_BY_MEMBER[key],
        )
        newly.append(record)
    return newly


def build_independent_compact(
    *,
    ssp: Mapping[str, Mapping[str, object]],
) -> dict[str, object]:
    payload = expected_compact_payload()
    for key, section in (
        ("SSPRK3-4097", "ssp4097"),
        ("SSPRK3-8193", "ssp8193"),
        ("SSPRK3-16385", "ssp16385"),
    ):
        record = _mapping(ssp.get(key), key)
        payload[section]["receipt_sha256"] = record["c1r1"]["receipt_sha256"]
        payload[section]["assessment_sha256"] = record["c1r1"]["assessment_sha256"]
        payload[section]["family_sha256"] = record["c1r1"]["family_sha256"]
        payload[section]["reconstructed_prepare_disposition"] = record[
            "reconstructed_prepare_disposition"
        ]
        payload[section]["overlay_count"] = record["overlay_count"]
        payload[section]["prepare_count"] = record["prepare_count"]
        if key == "SSPRK3-4097":
            payload[section]["independent_cfl"] = dict(record["independent_cfl"])
            payload[section]["next_plan_cap_hex"] = record["next_plan_cap_hex"]
    if payload != expected_compact_payload():
        _fail("independent compact facts differ from the bound certificate")
    result = {"artifact_id": ARTIFACT_ID, "artifact_payload": payload}
    _refuse(result)
    return result


def bind_live(
    repository: Path,
    *,
    origin: object | None = None,
    runtime: object | None = None,
    runtime_origin_builder: Callable[..., object] | None = None,
    source_evaluator: Callable[..., Mapping[str, object]] | None = None,
    static_input_bytes: Mapping[str, bytes] | None = None,
    environment: Mapping[str, str] | None = None,
    skip_git: bool = False,
    skip_raw: bool = False,
    rec1_raw_payload: Mapping[str, object] | None = None,
    attempt2_raw_payload: Mapping[str, object] | None = None,
    historical_runs: Mapping[str, object] | None = None,
    planck: Mapping[str, object] | None = None,
    clock: Callable[[], float] | None = None,
    rss_bytes: Callable[[], int] | None = None,
) -> dict[str, object]:
    """Independent live reconstruction. Prints only at the CLI; never writes."""

    import time as time_module

    root = _root(Path(repository))
    now = clock or time_module.perf_counter
    rss = rss_bytes or current_rss_bytes
    started = now()
    initial_rss = int(rss())
    env = dict(environment) if environment is not None else environment_identity()
    if skip_git:
        authority = {
            "authority_commit": AUTHORITY_COMMIT,
            "authority_parent": AUTHORITY_PARENT,
            "exact_delta_paths": list(ALLOWED_AUTHORITY_DELTA),
            "implementation_sha256": IMPLEMENTATION_SHA256,
        }
        freeze = {"environment": env, "source_closure_sha256": SOURCE_CLOSURE_SHA256}
        closure = SOURCE_CLOSURE_SHA256
    else:
        authority = authenticate_authority(root)
        freeze = authenticate_freeze_config(root, environment=env)
        source = source_closure_identity(
            root, AUTHORITY_PARENT, environment=env
        )
        if source["sha256"] != SOURCE_CLOSURE_SHA256:
            _fail("independent source closure differs")
        if source["file_count"] != SOURCE_CLOSURE_FILE_COUNT:
            _fail("independent source-closure file count differs")
        closure = str(source["sha256"])
        attempt2_compacts = authenticate_attempt2_compacts(root)
        del attempt2_compacts
    if skip_raw:
        if rec1_raw_payload is None or attempt2_raw_payload is None:
            _fail("synthetic live path omitted raw payloads")
        rec1_raw = dict(rec1_raw_payload)
        attempt2_raw = dict(attempt2_raw_payload)
    else:
        attempt2 = authenticate_attempt2_raw(root)
        rec1 = authenticate_rec1_raw(root)
        rec1_raw = rec1["payload"]
        attempt2_raw = attempt2["payload"]
    if historical_runs is None:
        runs = inventory_historical_runs(root)
    else:
        runs = dict(historical_runs)
        if (
            runs.get("file_count") != HISTORICAL_RUNS_FILE_COUNT
            or runs.get("byte_count") != HISTORICAL_RUNS_BYTE_COUNT
            or runs.get("digest") != HISTORICAL_RUNS_DIGEST
        ):
            _fail("historical runs baseline differs")
    if planck is None:
        planck_record = authenticate_planck(root)
    else:
        planck_record = dict(planck)
    if not skip_raw:
        require_pro20_namespace_absent(root)
    static = static_input_bytes
    if static is None:
        from .hlt17_srcq1 import read_static_factory_inputs

        static = read_static_factory_inputs(root)
    static_input_identities(static)
    loaded_origin = origin
    if loaded_origin is None:
        from .pro20_origin import capture_pro20_historical_origin

        loaded_origin = capture_pro20_historical_origin(root)
    if isinstance(loaded_origin, Pro20OriginCapture):
        try:
            revalidate_pro20_origin_capture(loaded_origin)
        except Pro20OriginError as error:
            raise HLT17SRCQ1REC1PREF1Error("origin capture is invalid") from error
    if getattr(loaded_origin, "capture_sha256", None) != ORIGIN_CAPTURE_SHA256:
        _fail("foreign origin identity")
    loaded_runtime = runtime
    if loaded_runtime is None:
        builder = runtime_origin_builder or build_pro20_runtime_origin
        try:
            loaded_runtime = builder(
                loaded_origin,
                static_input_bytes=static,
                source_closure_sha256=closure,
            )
        except Pro20SourceFactoryError as error:
            raise HLT17SRCQ1REC1PREF1Error(
                "physical static-factory/origin construction failed"
            ) from error
    cap_facts = derive_six_caps(loaded_origin, loaded_runtime)
    del cap_facts
    rk = bind_attempt2_rk_receipts(attempt2_raw, rec1_raw)
    ssp_records = reproduce_three_ssp_paths(
        loaded_origin,
        loaded_runtime,
        source_evaluator=source_evaluator,
    )
    raw_members = {
        _mapping(item, "raw member").get("member_key"): item
        for item in rec1_raw.get("members", [])
    }
    ssp_by_key = {item["member_key"]: item for item in ssp_records}
    for key in RESUME_MEMBER_KEYS:
        compare_ssp_to_raw(ssp_by_key[key], raw_members[key])
    wall = float(now() - started)
    observed_rss = int(rss())
    ceilings = recommended_resource_ceilings()
    _check_resources(ceilings, wall, observed_rss)
    _check_resources(ceilings, wall, initial_rss)
    facts = {
        "authority_bound": True,
        "freeze_bound": True,
        "attempt2_bound": True,
        "rec1_raw_bound": True,
        "historical_runs_bound": True,
        "planck_bound": True,
        "pro20_absent": True,
        "caps_derived": True,
        "rk_receipts_bound": True,
        "ssp4097_cfl_then_prepared": True,
        "ssp8193_fresh": True,
        "ssp16385_fresh": True,
        "eighteen_channels_bound": True,
        "no_state_advance": True,
        "resources_within_ceilings": True,
        "endpoints_absent": True,
        "runner_label_sufficient": False,
    }
    classification = classification_if_established(facts)
    compact = build_independent_compact(ssp=ssp_by_key)
    if compact["artifact_payload"]["classification"] != classification:
        _fail("compact classification drifted")
    _refuse(compact)
    return {
        "config": expected_config(),
        "result": compact,
        "classification": classification,
        "written": False,
        "authority": authority,
        "freeze_source_closure_sha256": freeze.get("source_closure_sha256"),
        "historical_runs": runs,
        "planck": planck_record,
        "rk": rk,
        "live_wall_seconds": wall,
        "live_rss_bytes": observed_rss,
    }


__all__ = (
    "ABSENT_COMPACT_MESSAGE",
    "ARTIFACT_ID",
    "AUTHORITY_COMMIT",
    "AUTHORITY_PARENT",
    "CLASSIFICATION",
    "CONFIG_PATH",
    "HLT17SRCQ1REC1PREF1Error",
    "OWNER_DOCUMENT",
    "RESULT_PATH",
    "bind_live",
    "canonical_result",
    "classification_if_established",
    "emit_config_bytes",
    "expected_compact",
    "expected_config",
    "historical_runs_digest",
    "independent_cfl_ratio",
    "render_config",
    "reproduce_ssp_member",
    "validate_compact_result",
    "verify_compact",
)
