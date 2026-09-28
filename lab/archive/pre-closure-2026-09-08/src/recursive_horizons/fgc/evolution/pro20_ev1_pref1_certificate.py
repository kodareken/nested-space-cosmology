"""Independent post-terminal certificate for the one authorized PRO20 event.

This module consumes the import-independent PREF1 binder core.  It does not
import PRO20 runner, protocol, store, or authority decision code.  Ordinary
compact verification is Git-, raw-, store-, source-, Planck-, and replay-blind.
Live binding authenticates Git/freeze/root/store/terminal/historical/Planck/
environment/source identities and independently rebuilds the six seed bundles.
The authorized event created the production campaign store and persisted 146
accepted fine states.  The binder and certificate write nothing, serialize or
adopt no additional endpoint, resume nothing, and do not rerun the event.
"""

from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path, PurePosixPath
import platform
import re
import stat
import tomllib
from types import MappingProxyType
from typing import Any, Mapping

from recursive_horizons.evidence_io import (
    CommitParents,
    EvidenceIOError,
    InspectDelta,
    InspectTree,
    ReadBlob,
    UnsafePathError,
    canonical_json_bytes,
    git_read,
    read_regular_file,
)

from .hlt17_imp1_cursor import HLT17_MEMBER_KEYS
from .hlt17_member_checkpoint import HLT17GenerationBundle
from .hlt17_srcq1 import observe_environment, read_static_factory_inputs
from .pro20_ev1_pref1_binder import (
    ARTIFACT_ID,
    PRODUCTION_NAMESPACE,
    RECORD_KINDS,
    bind_independent_seed,
    bind_terminal_store,
    independent_outcome,
    independent_progression,
    snapshot_store,
)
from .pro20_ev1_pref1_binder import PRO20EV1PREF1Error as _BinderError
from .pro20_origin import (
    Pro20OriginError,
    capture_pro20_historical_origin,
    revalidate_pro20_origin_capture,
)
from .pro20_source_factory import Pro20SourceFactoryError, build_pro20_runtime_origin


SCHEMA_VERSION = 1
CLASSIFICATION = "independently_bound_resource_exhausted_no_common_event"
CONFIG_PATH = "configs/fgc/fgc-1-pro20-ev1-pref1.toml"
RESULT_PATH = "results/fgc-1-pro20-ev1-pref1.json"
OWNER_DOCUMENT = "docs/fgc-pro20-ev1-pref1.md"
FREEZE_ARTIFACT_ID = "FGC-1-PRO20-EV1-FRZ1"
FREEZE_CONFIG_PATH = "configs/fgc/fgc-1-pro20-ev1-frz1.toml"
FREEZE_OWNER_PATH = "docs/fgc-pro20-ev1-frz1.md"
FREEZE_MAKE_PATH = "mk/current-foundation.mk"
AUTHORITY_COMMIT = "40c870045d04dc732be153699edf11a409ab8737"
AUTHORITY_PARENT = "6e24358319940d31a4b914663d698e88e6c73905"
ALLOWED_AUTHORITY_DELTA = (
    "CHANGELOG.md",
    "README.md",
    "configs/fgc/artifact-catalog.json",
    "configs/fgc/fgc-1-pro20-ev1-frz1.toml",
    "docs/active-code-map.md",
    "docs/claim-ledger.md",
    "docs/fgc-pro20-ev1-frz1.md",
    "docs/fgc-runtime-matrix.md",
    "docs/research-roadmap.md",
    "mk/current-foundation.mk",
    "results/README.md",
    "scripts/build_artifact_catalog.py",
    "scripts/repo_checks/catalog.py",
    "tests/test_phase_minus1_artifact_catalog.py",
    "tests/test_phase_minus1_make_routing.py",
)
FREEZE_CONFIG_SHA256 = (
    "5dc5dff9a2a52accd3916fb74effa355e73207ea578293efa0c0888c3e63c0c8"
)
FREEZE_OWNER_SHA256 = (
    "affb6cfaf9dae901e157ab4d688dc1daa08906f7f72b7e6d72c085c397dfc6ff"
)
FREEZE_MAKE_SHA256 = (
    "066ca6d4a5d8566d8912416faeb212fdb3675041440146da7677a6b8bb8bce33"
)
IMPLEMENTATION_SHA256 = (
    "aa1dd7b6ef041f8bb17ac8dd4c26fa52d83c319d5d2ba12b2f157d0620103864"
)
AUTHORITY_MODULE_SHA256 = (
    "24eec1e99bfa400c50f32ae308b491baf51bc302d3a12f01036b7f414e56f900"
)
RUNNER_SHA256 = (
    "c03583ae57b648dc79537ef354b9be61b022b164f88bc3cb9dd69e922f00c434"
)
SOURCE_CLOSURE_SHA256 = (
    "5ad9f84dd9f38c7f47f591e4a55ce525861abcb20265354cdea1b4b564f67922"
)
SOURCE_CLOSURE_FILE_COUNT = 253
ORIGIN_CAPTURE_SHA256 = (
    "2824af06c674336d479ec53faa4269f02fb598055860ec22f2ad84602a83fe72"
)
ENVIRONMENT_SHA256 = (
    "d9638dba30833371e67d151db62e77268f92878f14eeea8a3031804b113cb73a"
)
AUTHORITY_RECEIPT_SHA256 = (
    "77b5978bdead70c7ca193728e6d8895208e396314f09cf70906fe6f92e222242"
)
SOURCE_CLOSURE_DOMAIN = b"FGC-1-PRO20-EV1-SOURCE-CLOSURE-v1\n"
SOURCE_TREE = "src/recursive_horizons"
HARNESS_PATH = "src/recursive_horizons/fgc/evolution/pro20_ev1_runtime.py"
RUNNER_PATH = "scripts/run_fgc_pro20_ev1.py"
AUTHORITY_PATH = "src/recursive_horizons/fgc/evolution/pro20_ev1_authority.py"
CAMPAIGN_ID = "FGC-2-SF1-PRO20-EVENT1"
STORE_LEAF_COUNT = 1081
STORE_BYTE_COUNT = 155341592
STORE_LEXICAL_SHA256 = (
    "3104609e029f686221a6ff66e48e098ee394b600fd5a9cac40c76f5a9fd28a1c"
)
CHECKPOINT_SHA256 = (
    "de809e456fafe74364c29165f963734722361621b878ef2d11fabdc5358d76b4"
)
JOURNAL_SHA256 = (
    "a391d38614b065c7ec3bed187446efb5734645e61e00bd942d502b5f77e7e553"
)
TERMINAL_LOCK_SHA256 = (
    "c9503cf85348f24c2251e83c9479c7ce594bc28b16ae4087982c290f2c5ba607"
)
SESSION_COUNT = 1
GENERATION_COUNT = 192
ATTEMPT_COUNT = 191
KIND_COUNTS = MappingProxyType(
    {
        "seed": 1,
        "accepted_fine": 146,
        "source_retry": 0,
        "cfl_retry": 44,
        "temporal_retry": 0,
        "source_retry_exhausted": 0,
        "cfl_retry_exhausted": 0,
        "temporal_retry_exhausted": 0,
        "invalid_premise": 0,
        "resource_exhausted": 1,
    }
)
PARKED_MEMBER_KEYS = (
    "RK4-2049",
    "RK4-4097",
    "RK4-8193",
    "SSPRK3-4097",
)
ACCEPTED_GENERATIONS = MappingProxyType(
    {
        "RK4-2049": 14,
        "RK4-4097": 28,
        "RK4-8193": 39,
        "SSPRK3-4097": 29,
        "SSPRK3-8193": 36,
        "SSPRK3-16385": 0,
    }
)
ACCEPTED_TIMES = MappingProxyType(
    {
        "RK4-2049": "0x1.8000000000000p+0",
        "RK4-4097": "0x1.8000000000000p+0",
        "RK4-8193": "0x1.8000000000000p+0",
        "SSPRK3-4097": "0x1.8000000000000p+0",
        "SSPRK3-8193": "0x1.7effffffffe80p+0",
        "SSPRK3-16385": "0x1.7000000000000p+0",
    }
)
TERMINAL_MEMBER_KEY = "SSPRK3-8193"
TERMINAL_REASON = "max RSS exceeded"
ORIGIN_TIME_HEX = "0x1.7000000000000p+0"
SEED_IDENTITY_STREAM_SHA256 = (
    "cde71e2585f82649625283cc916b4a85ab8743918cb029a90a074c0fa29fb57f"
)
HISTORICAL_RUNS_ROOT = "runs"
HISTORICAL_RUNS_FILE_COUNT = 300
HISTORICAL_RUNS_BYTE_COUNT = 150437927
HISTORICAL_RUNS_DIGEST = (
    "77aa79f0300bd6866fcb491fe38aabf813a0064cbc7d511ef067364b4998ed88"
)
REC1_OUTPUT_NAMESPACE = (
    "runs/fgc-2-sf1/hlt17-srcq1-rec1-source-origin-qualification"
)
PRODUCTION_CONTAINER = "runs/fgc-2-sf1/pro20-event1"
EXCLUDED_HISTORICAL_NAMESPACES = (
    REC1_OUTPUT_NAMESPACE,
    PRODUCTION_CONTAINER,
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
    "3ec879fdfd309864b913b036f6b8c08c4717cb40269a18749cc19b8c0ca7b09b"
)
ABSENT_COMPACT_MESSAGE = (
    "FGC-1-PRO20-EV1-PREF1 compact result is absent at "
    f"{RESULT_PATH}. Default --check validates tracked compact bytes only "
    "and never reconstructs, reads raw/runs/store, or inspects Git."
)
NONCLAIMS = (
    "PREF1 binds a typed resource_exhausted stop with partial accepted "
    "state and no common event.",
    "The authorized event created the production campaign store and "
    "persisted 146 accepted fine states; the independent binder and "
    "certificate wrote nothing and serialized or adopted no additional "
    "endpoint.",
    "PREF1 does not authorize PRO21, Wave 2, calibration, SGB-L, DEF1, "
    "holdout, candidate, mechanism, or physics.",
    "PREF1 does not authorize resume, takeover, a new namespace, or another "
    "production run.",
    "Runner kind and disposition labels are never sufficient.",
    "RSRC1 is a planning owner only; any new measurement needs a separately "
    "committed RSRC1-FRZ1.",
)
FORBIDDEN_UNQUALIFIED_KEYS = (
    "campaign_store_written",
    "endpoint_adopted_or_serialized",
    "endpoint_serialized",
)
EXPECTED_CLAIMS = MappingProxyType(
    {
        "accepted_fine_states_persisted": 146,
        "accepted_state_advanced": True,
        "all_six_parked": False,
        "binder_mutated_store": False,
        "binder_serialized_or_adopted_endpoint": False,
        "binder_wrote_campaign_store": False,
        "calibration_authorized": False,
        "calibration_eligible": False,
        "candidate_execution_authorized": False,
        "common_event_completed": False,
        "def1_error_map_passed": False,
        "first_event_complete": False,
        "holdout_authorized": False,
        "independently_bound_resource_exhausted": True,
        "mechanism_claimed": False,
        "partial_accepted_state": True,
        "physics_claimed": False,
        "pro21_authorized": False,
        "production_campaign_store_present": True,
        "production_store_written_by_event": True,
        "resume_authorized": False,
        "runner_label_sufficient": False,
        "sgbl_health_established": False,
        "ssprk3_16385_untouched": True,
        "takeover_authorized": False,
        "wave2_authorized": False,
        "written": False,
    }
)


class PRO20EV1PREF1CertificateError(RuntimeError):
    """The independent PRO20 PREF1 certificate could not be established."""


def _fail(message: str) -> None:
    raise PRO20EV1PREF1CertificateError(message)


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
        raise PRO20EV1PREF1CertificateError(f"{label} is not JSON") from error
    return _mapping(parsed, label)


def _require_bool(value: object, expected: bool, label: str) -> bool:
    if type(value) is not bool or value is not expected:
        _fail(f"{label} differs")
    return value


def _refuse_unqualified_store_fields(mapping: Mapping[str, object]) -> None:
    for key in FORBIDDEN_UNQUALIFIED_KEYS:
        if key in mapping:
            _fail(f"unqualified store field {key} is present")


def _root(value: Path) -> Path:
    if not isinstance(value, Path):
        value = Path(value)
    if not value.is_absolute():
        _fail("repository root is not canonical")
    try:
        resolved = value.resolve(strict=True)
    except OSError as error:
        raise PRO20EV1PREF1CertificateError(
            "repository root cannot be resolved"
        ) from error
    if resolved != value:
        _fail("repository root traverses a symlink")
    return value


def _same(actual: object, expected: object, label: str) -> None:
    if actual != expected:
        _fail(f"{label} differs")


def canonical_result(value: object) -> bytes:
    return canonical_json_bytes(value) + b"\n"


def expected_claims() -> dict[str, object]:
    return dict(EXPECTED_CLAIMS)


def expected_kind_counts() -> dict[str, int]:
    counts = {kind: 0 for kind in RECORD_KINDS}
    counts.update(dict(KIND_COUNTS))
    return counts


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
        "partial_accepted_state": True,
        "accepted_state_advanced": True,
        "common_event_completed": False,
        "production_campaign_store_present": True,
        "production_store_written_by_event": True,
        "accepted_fine_states_persisted": 146,
        "binder_mutated_store": False,
        "binder_wrote_campaign_store": False,
        "binder_serialized_or_adopted_endpoint": False,
        "pro21_wave2_authorized": False,
        "calibration_authorized": False,
        "written": False,
    }


def render_config(config: Mapping[str, object]) -> bytes:
    mapping = _mapping(config, "PREF1 config")
    _refuse_unqualified_store_fields(mapping)
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


def expected_compact_payload() -> dict[str, object]:
    return {
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
            "authority_module_sha256": AUTHORITY_MODULE_SHA256,
            "runner_sha256": RUNNER_SHA256,
            "source_closure_sha256": SOURCE_CLOSURE_SHA256,
            "origin_capture_sha256": ORIGIN_CAPTURE_SHA256,
            "environment_sha256": ENVIRONMENT_SHA256,
            "authority_receipt_sha256": AUTHORITY_RECEIPT_SHA256,
        },
        "store": {
            "namespace": PRODUCTION_NAMESPACE,
            "campaign_id": CAMPAIGN_ID,
            "leaf_count": STORE_LEAF_COUNT,
            "byte_count": STORE_BYTE_COUNT,
            "lexical_sha256": STORE_LEXICAL_SHA256,
            "checkpoint_sha256": CHECKPOINT_SHA256,
            "journal_sha256": JOURNAL_SHA256,
            "terminal_lock_sha256": TERMINAL_LOCK_SHA256,
            "session_count": SESSION_COUNT,
        },
        "progression": {
            "generation_count_including_seed": GENERATION_COUNT,
            "attempt_count": ATTEMPT_COUNT,
            "kind_counts": expected_kind_counts(),
            "accepted_fine_count": 146,
            "accepted_state_advanced": True,
            "partial_accepted_state": True,
            "parked_member_keys": list(PARKED_MEMBER_KEYS),
            "parked_member_count": len(PARKED_MEMBER_KEYS),
            "member_keys": list(HLT17_MEMBER_KEYS),
            "accepted_generation_by_member": dict(ACCEPTED_GENERATIONS),
            "accepted_generation_in_canonical_member_order": [
                ACCEPTED_GENERATIONS[key] for key in HLT17_MEMBER_KEYS
            ],
            "accepted_time_hex_by_member": dict(ACCEPTED_TIMES),
            "terminal_member_key": TERMINAL_MEMBER_KEY,
            "terminal_kind": "resource_exhausted",
            "terminal_reason": TERMINAL_REASON,
            "ssprk3_8193_accepted_time_hex": ACCEPTED_TIMES["SSPRK3-8193"],
            "ssprk3_16385_accepted_time_hex": ORIGIN_TIME_HEX,
            "all_six_parked": False,
            "common_event_completed": False,
        },
        "seed": {
            "identity_stream_sha256": SEED_IDENTITY_STREAM_SHA256,
            "independently_reconstructed_origin_matches_raw_seed": True,
            "root_manifest_origin_label_sufficient": False,
            "state_advanced_by_seed_binding": False,
        },
        "historical_runs_baseline": {
            "file_count": HISTORICAL_RUNS_FILE_COUNT,
            "byte_count": HISTORICAL_RUNS_BYTE_COUNT,
            "digest": HISTORICAL_RUNS_DIGEST,
            "excluded_namespaces": list(EXCLUDED_HISTORICAL_NAMESPACES),
        },
        "planck": {
            path: {"byte_count": size, "sha256": digest}
            for path, (size, digest) in PLANCK_FILES.items()
        },
        "claims": expected_claims(),
        "nonclaims": list(NONCLAIMS),
        "written": False,
        "production_campaign_store_present": True,
        "production_store_written_by_event": True,
        "accepted_fine_states_persisted": 146,
        "binder_mutated_store": False,
        "binder_wrote_campaign_store": False,
        "binder_serialized_or_adopted_endpoint": False,
        "pro21_wave2_authorized": False,
    }


def expected_compact() -> dict[str, object]:
    return {
        "artifact_id": ARTIFACT_ID,
        "artifact_payload": expected_compact_payload(),
    }


def _read_tracked(root: Path, relative: str) -> bytes:
    try:
        return read_regular_file(root, relative)
    except (EvidenceIOError, UnsafePathError, OSError) as error:
        raise PRO20EV1PREF1CertificateError(
            f"tracked PREF1 leaf is absent: {relative}"
        ) from error


def verify_compact(repository: Path) -> dict[str, object]:
    """Validate tracked compact bytes only. No Git, raw, store, source, Planck, or replay."""

    root = Path(repository)
    try:
        config_raw = _read_tracked(root, CONFIG_PATH)
        result_raw = _read_tracked(root, RESULT_PATH)
    except PRO20EV1PREF1CertificateError as error:
        raise PRO20EV1PREF1CertificateError(ABSENT_COMPACT_MESSAGE) from error
    return validate_compact_result(config_raw, result_raw)


def validate_compact_result(config_raw: bytes, result_raw: bytes) -> dict[str, object]:
    if type(config_raw) is not bytes or type(result_raw) is not bytes:
        _fail("compact bytes differ")
    if render_config(expected_config()) != config_raw:
        try:
            parsed = tomllib.loads(config_raw.decode("utf-8"))
        except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
            raise PRO20EV1PREF1CertificateError("PREF1 config is malformed") from error
        _refuse_unqualified_store_fields(parsed)
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
    _refuse_unqualified_store_fields(payload)
    claims = _mapping(payload.get("claims"), "PREF1 claims")
    _refuse_unqualified_store_fields(claims)
    if payload.get("classification") != CLASSIFICATION:
        _fail("PREF1 classification differs")
    if payload != expected["artifact_payload"]:
        _fail("PREF1 compact payload differs")
    if claims != expected_claims():
        _fail("PREF1 claims differ")
    _require_bool(payload.get("written"), False, "compact written")
    _require_bool(
        payload.get("production_campaign_store_present"),
        True,
        "production_campaign_store_present",
    )
    _require_bool(
        payload.get("production_store_written_by_event"),
        True,
        "production_store_written_by_event",
    )
    if payload.get("accepted_fine_states_persisted") != 146:
        _fail("accepted_fine_states_persisted differs")
    _require_bool(payload.get("binder_mutated_store"), False, "binder_mutated_store")
    _require_bool(
        payload.get("binder_wrote_campaign_store"),
        False,
        "binder_wrote_campaign_store",
    )
    _require_bool(
        payload.get("binder_serialized_or_adopted_endpoint"),
        False,
        "binder_serialized_or_adopted_endpoint",
    )
    _require_bool(
        payload.get("pro21_wave2_authorized"), False, "pro21_wave2_authorized"
    )
    return result


def classification_if_established(facts: Mapping[str, object]) -> str:
    mapping = _mapping(facts, "independent facts")
    required = (
        "authority_bound",
        "freeze_bound",
        "root_bound",
        "store_bound",
        "terminal_bound",
        "historical_runs_bound",
        "planck_bound",
        "environment_bound",
        "source_bound",
        "seed_independently_rebuilt",
        "resource_exhausted_terminal",
        "no_common_event",
        "partial_accepted_state",
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
        raise PRO20EV1PREF1CertificateError(
            "committed image cannot be inspected"
        ) from error


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


def _hash_software_file(path: str) -> str:
    if type(path) is not str or not path or "\0" in path or len(path) > 1024:
        _fail("software path is unsafe")
    try:
        resolved = os.path.realpath(path)
        before = os.lstat(resolved)
    except OSError as error:
        raise PRO20EV1PREF1CertificateError("software path cannot be read") from error
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
        raise PRO20EV1PREF1CertificateError(
            "software path cannot be opened"
        ) from error
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
        raise PRO20EV1PREF1CertificateError(
            "runtime path cannot be resolved"
        ) from error
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
    executable = _runtime_regular_path(str(observed["executable"]))
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


def environment_sha256(environment: Mapping[str, str]) -> str:
    if tuple(environment) != ENVIRONMENT_KEYS:
        _fail("environment keys differ")
    return _sha256_hex(canonical_json_bytes(dict(environment)))


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
    for path in (RUNNER_PATH, AUTHORITY_PATH, HARNESS_PATH):
        if path not in seen:
            paths.append(path)
            seen.add(path)
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
        _fail("PRO20 authority must have one implementation parent")
    parent = parents[0]
    if parent != AUTHORITY_PARENT:
        _fail("PRO20 authority parent is not the implementation commit")
    delta = _git(repository, InspectDelta(commit=authority, against=parent))
    paths = tuple(sorted(item.path for item in delta))
    if paths != tuple(sorted(ALLOWED_AUTHORITY_DELTA)):
        _fail("PRO20 authority delta paths differ")
    freeze_hashes = {
        FREEZE_CONFIG_PATH: FREEZE_CONFIG_SHA256,
        FREEZE_OWNER_PATH: FREEZE_OWNER_SHA256,
        FREEZE_MAKE_PATH: FREEZE_MAKE_SHA256,
    }
    for path, digest in freeze_hashes.items():
        committed = _blob(repository, authority, path)
        if _sha256_hex(committed) != digest:
            _fail(f"authority freeze blob differs: {path}")
    live_config = read_regular_file(repository, FREEZE_CONFIG_PATH)
    if _sha256_hex(live_config) != FREEZE_CONFIG_SHA256:
        _fail("live freeze config SHA-256 differs from the authority blob")
    if live_config != _blob(repository, authority, FREEZE_CONFIG_PATH):
        _fail("live freeze config differs from the authority commit")
    implementation = _sha256_hex(_blob(repository, parent, HARNESS_PATH))
    runner = _sha256_hex(_blob(repository, parent, RUNNER_PATH))
    authority_module = _sha256_hex(_blob(repository, parent, AUTHORITY_PATH))
    if (
        implementation != IMPLEMENTATION_SHA256
        or runner != RUNNER_SHA256
        or authority_module != AUTHORITY_MODULE_SHA256
    ):
        _fail("implementation/runner/authority hash differs")
    return {
        "authority_commit": authority,
        "authority_parent": parent,
        "exact_delta_paths": list(ALLOWED_AUTHORITY_DELTA),
        "implementation_sha256": implementation,
        "runner_sha256": runner,
        "authority_module_sha256": authority_module,
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
        raise PRO20EV1PREF1CertificateError("freeze config is malformed") from error
    mapping = _mapping(config, "freeze config")
    if mapping.get("artifact_id") != FREEZE_ARTIFACT_ID:
        _fail("freeze artifact differs")
    if mapping.get("implementation_commit") != AUTHORITY_PARENT:
        _fail("freeze implementation parent differs")
    if mapping.get("implementation_sha256") != IMPLEMENTATION_SHA256:
        _fail("freeze implementation SHA-256 differs")
    if mapping.get("authority_sha256") != AUTHORITY_MODULE_SHA256:
        _fail("freeze authority SHA-256 differs")
    if mapping.get("runner_sha256") != RUNNER_SHA256:
        _fail("freeze runner SHA-256 differs")
    if mapping.get("source_closure_sha256") != SOURCE_CLOSURE_SHA256:
        _fail("freeze source closure differs")
    if mapping.get("origin_capture_sha256") != ORIGIN_CAPTURE_SHA256:
        _fail("freeze origin capture differs")
    observed = dict(environment)
    if mapping.get("environment") != observed:
        _fail("freeze environment differs")
    if environment_sha256(observed) != ENVIRONMENT_SHA256:
        _fail("freeze environment SHA-256 differs")
    return mapping


def historical_runs_digest(
    records: Mapping[str, tuple[int, str]] | list[tuple[str, int, str]]
) -> str:
    """SHA-256 of sha256sum lines ``digest`` two spaces ``path`` newline."""

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


def _require_relative(relative: str, label: str) -> Path:
    value = Path(relative)
    if value.is_absolute() or any(part in {"", ".", ".."} for part in value.parts):
        _fail(f"{label} is unsafe")
    return value


def _excluded_historical(relative: str) -> bool:
    for prefix in EXCLUDED_HISTORICAL_NAMESPACES:
        if relative == prefix or relative.startswith(prefix + "/"):
            return True
    return False


def inventory_historical_runs(root: Path) -> dict[str, object]:
    repository = _root(root)
    relative = _require_relative(HISTORICAL_RUNS_ROOT, "runs root")
    base = repository / relative
    try:
        info = base.lstat()
    except OSError as error:
        raise PRO20EV1PREF1CertificateError("historical runs tree is absent") from error
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
        _fail("historical runs tree is unsafe")
    pending = [str(relative)]
    files: list[str] = []
    while pending:
        current_relative = pending.pop()
        current = repository / current_relative
        try:
            with os.scandir(current) as iterator:
                children = sorted(iterator, key=lambda item: item.name)
        except OSError as error:
            raise PRO20EV1PREF1CertificateError(
                "historical runs tree cannot be scanned"
            ) from error
        for child in children:
            child_relative = f"{current_relative}/{child.name}"
            try:
                child_info = child.stat(follow_symlinks=False)
            except OSError as error:
                raise PRO20EV1PREF1CertificateError(
                    f"historical runs path is unsafe: {child_relative}"
                ) from error
            if stat.S_ISLNK(child_info.st_mode):
                _fail(f"historical runs path is unsafe: {child_relative}")
            if stat.S_ISDIR(child_info.st_mode):
                pending.append(child_relative)
                continue
            if not stat.S_ISREG(child_info.st_mode) or child_info.st_nlink != 1:
                _fail(f"historical runs path is unsafe: {child_relative}")
            if _excluded_historical(child_relative):
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
            raise PRO20EV1PREF1CertificateError(
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
        "excluded_namespaces": list(EXCLUDED_HISTORICAL_NAMESPACES),
    }


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


def store_identity(repository: Path) -> dict[str, object]:
    """Read-only store identity. Never writes."""

    try:
        snapshot = snapshot_store(Path(repository))
    except _BinderError as error:
        raise PRO20EV1PREF1CertificateError(str(error)) from error
    return {
        "lexical_sha256": snapshot.lexical_sha256,
        "byte_count": snapshot.byte_count,
        "leaf_count": len(snapshot.leaves),
        "leaf_sha256s": {
            relative: leaf.sha256 for relative, leaf in snapshot.leaves.items()
        },
    }


def rebuild_independent_seed_bundles(
    repository: Path,
    *,
    source_closure_sha256: str,
    origin: object | None = None,
    runtime: object | None = None,
    static_input_bytes: Mapping[str, bytes] | None = None,
) -> dict[str, HLT17GenerationBundle]:
    """Rebuild six origin bundles from lower origin/factory owners. No PDE replay."""

    loaded_origin = origin
    if loaded_origin is None:
        try:
            loaded_origin = capture_pro20_historical_origin(repository)
        except Pro20OriginError as error:
            raise PRO20EV1PREF1CertificateError(
                "historical origin capture failed"
            ) from error
    try:
        revalidate_pro20_origin_capture(loaded_origin)
    except Pro20OriginError as error:
        raise PRO20EV1PREF1CertificateError("origin capture is invalid") from error
    if getattr(loaded_origin, "capture_sha256", None) != ORIGIN_CAPTURE_SHA256:
        _fail("foreign origin identity")
    loaded_runtime = runtime
    if loaded_runtime is None:
        static = static_input_bytes
        if static is None:
            static = read_static_factory_inputs(repository)
        try:
            loaded_runtime = build_pro20_runtime_origin(
                loaded_origin,
                static_input_bytes=static,
                source_closure_sha256=source_closure_sha256,
            )
        except Pro20SourceFactoryError as error:
            raise PRO20EV1PREF1CertificateError(
                "independent seed reconstruction failed"
            ) from error
    if getattr(loaded_runtime, "captured_origin_sha256", None) != ORIGIN_CAPTURE_SHA256:
        _fail("runtime origin capture identity differs")
    if getattr(loaded_runtime, "source_closure_sha256", None) != source_closure_sha256:
        _fail("runtime origin source identity differs")
    members = getattr(loaded_runtime, "members", None)
    if not isinstance(members, Mapping) or tuple(members) != tuple(HLT17_MEMBER_KEYS):
        _fail("independent seed member order differs")
    bundles: dict[str, HLT17GenerationBundle] = {}
    for key in HLT17_MEMBER_KEYS:
        bundle = members[key].bundle
        if not isinstance(bundle, HLT17GenerationBundle):
            _fail(f"independent seed bundle type differs for {key}")
        bundles[key] = bundle
    return bundles


def _bind_store(repository: Path):
    try:
        return bind_terminal_store(repository)
    except _BinderError as error:
        raise PRO20EV1PREF1CertificateError(str(error)) from error


def _compare_terminal(terminal: object) -> dict[str, object]:
    try:
        outcome = independent_outcome(terminal)
        progression = independent_progression(terminal)
    except _BinderError as error:
        raise PRO20EV1PREF1CertificateError(str(error)) from error
    if outcome.get("classification") != CLASSIFICATION:
        _fail("independent outcome classification differs")
    if outcome.get("terminal_kind") != "resource_exhausted":
        _fail("independent terminal kind is not resource_exhausted")
    if outcome.get("common_event_completed") is not False:
        _fail("common event must remain incomplete")
    if getattr(terminal, "checkpoint_sha256", None) != CHECKPOINT_SHA256:
        _fail("terminal checkpoint SHA-256 differs")
    if getattr(terminal, "journal_sha256", None) != JOURNAL_SHA256:
        _fail("terminal journal SHA-256 differs")
    if getattr(terminal, "terminal_lock_sha256", None) != TERMINAL_LOCK_SHA256:
        _fail("terminal lock SHA-256 differs")
    if getattr(terminal, "session_count", None) != SESSION_COUNT:
        _fail("closed session count differs")
    if getattr(terminal, "store_leaf_count", None) != STORE_LEAF_COUNT:
        _fail("store leaf count differs")
    if getattr(terminal, "store_byte_count", None) != STORE_BYTE_COUNT:
        _fail("store byte count differs")
    if getattr(terminal, "store_lexical_sha256", None) != STORE_LEXICAL_SHA256:
        _fail("store lexical SHA-256 differs")
    if progression.get("generation_count_including_seed") != GENERATION_COUNT:
        _fail("generation count differs")
    if progression.get("attempt_count") != ATTEMPT_COUNT:
        _fail("attempt count differs")
    if progression.get("kind_counts") != expected_kind_counts():
        _fail("kind counts differ")
    if tuple(progression.get("parked_member_keys") or ()) != PARKED_MEMBER_KEYS:
        _fail("parked member keys differ")
    if progression.get("accepted_generation_by_member") != dict(ACCEPTED_GENERATIONS):
        _fail("accepted generations differ")
    if progression.get("accepted_time_hex_by_member") != dict(ACCEPTED_TIMES):
        _fail("accepted times differ")
    if progression.get("accepted_state_advanced") is not True:
        _fail("partial accepted state must remain explicit")
    if progression.get("accepted_fine_count") != 146:
        _fail("accepted fine count differs")
    if progression.get("binder_mutated_store") is not False:
        _fail("binder mutated the production store")
    if progression.get("common_event_completed") is not False:
        _fail("progression reported a common event")
    last = terminal.generations[-1]
    if last.changed_member_key != TERMINAL_MEMBER_KEY:
        _fail("terminal member differs")
    if last.derived_kind != "resource_exhausted":
        _fail("derived terminal kind differs")
    if ACCEPTED_TIMES["SSPRK3-16385"] != ORIGIN_TIME_HEX:
        _fail("SSPRK3-16385 origin time constant drifted")
    if last.member_states["SSPRK3-16385"]["accepted_time_hex"] != ORIGIN_TIME_HEX:
        _fail("SSPRK3-16385 is not the untouched origin time")
    if last.member_states["SSPRK3-16385"]["accepted_generation"] != 0:
        _fail("SSPRK3-16385 accepted generation is not zero")
    if last.member_states["SSPRK3-8193"]["accepted_time_hex"] != ACCEPTED_TIMES[
        "SSPRK3-8193"
    ]:
        _fail("SSPRK3-8193 accepted time differs")
    return {"outcome": outcome, "progression": progression}


def bind_live(
    repository: Path,
    *,
    origin: object | None = None,
    runtime: object | None = None,
    static_input_bytes: Mapping[str, bytes] | None = None,
    environment: Mapping[str, str] | None = None,
    skip_git: bool = False,
    skip_store: bool = False,
    skip_source: bool = False,
    terminal: object | None = None,
    seed_bundles: Mapping[str, HLT17GenerationBundle] | None = None,
    historical_runs: Mapping[str, object] | None = None,
    planck: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """Independent no-write live binding. Never publishes or reruns the event.

    The production campaign store is the event's closed terminal. This path
    authenticates it and must leave it unchanged.
    """

    root = _root(Path(repository))
    before = None if skip_store else store_identity(root)
    env = dict(environment) if environment is not None else environment_identity()
    if environment_sha256(env) != ENVIRONMENT_SHA256:
        _fail("live environment SHA-256 differs")
    if skip_git:
        authority = {
            "authority_commit": AUTHORITY_COMMIT,
            "authority_parent": AUTHORITY_PARENT,
            "exact_delta_paths": list(ALLOWED_AUTHORITY_DELTA),
            "implementation_sha256": IMPLEMENTATION_SHA256,
            "runner_sha256": RUNNER_SHA256,
            "authority_module_sha256": AUTHORITY_MODULE_SHA256,
        }
        freeze = {"source_closure_sha256": SOURCE_CLOSURE_SHA256}
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
    loaded_terminal = terminal
    if loaded_terminal is None:
        if skip_store:
            _fail("synthetic live path omitted the bound terminal")
        loaded_terminal = _bind_store(root)
    compared = _compare_terminal(loaded_terminal)
    if (
        not skip_store
        and getattr(loaded_terminal, "campaign_id", None) != CAMPAIGN_ID
    ):
        _fail("campaign identity differs")
    rebuilt = seed_bundles
    if rebuilt is None:
        if skip_source:
            _fail("synthetic live path omitted independent seed bundles")
        rebuilt = rebuild_independent_seed_bundles(
            root,
            source_closure_sha256=closure,
            origin=origin,
            runtime=runtime,
            static_input_bytes=static_input_bytes,
        )
    try:
        seed = bind_independent_seed(loaded_terminal, rebuilt)
    except _BinderError as error:
        raise PRO20EV1PREF1CertificateError(str(error)) from error
    if seed.get("identity_stream_sha256") != SEED_IDENTITY_STREAM_SHA256:
        _fail("independent seed identity stream differs")
    if seed.get("independently_reconstructed_origin_matches_raw_seed") is not True:
        _fail("independent seed does not match raw generation zero")
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
    after = None if skip_store else store_identity(root)
    if before is not None and after != before:
        _fail("live PREF1 mutated the production store")
    facts = {
        "authority_bound": True,
        "freeze_bound": True,
        "root_bound": True,
        "store_bound": True,
        "terminal_bound": True,
        "historical_runs_bound": True,
        "planck_bound": True,
        "environment_bound": True,
        "source_bound": True,
        "seed_independently_rebuilt": True,
        "resource_exhausted_terminal": True,
        "no_common_event": True,
        "partial_accepted_state": True,
        "runner_label_sufficient": False,
    }
    classification = classification_if_established(facts)
    compact = expected_compact()
    if compact["artifact_payload"]["classification"] != classification:
        _fail("compact classification drifted")
    return {
        "config": expected_config(),
        "result": compact,
        "classification": classification,
        "written": False,
        "authority": authority,
        "freeze_source_closure_sha256": freeze.get("source_closure_sha256"),
        "historical_runs": runs,
        "planck": planck_record,
        "seed": seed,
        "outcome": compared["outcome"],
        "progression": compared["progression"],
        "store_identity": after if after is not None else before,
    }


__all__ = (
    "ABSENT_COMPACT_MESSAGE",
    "ARTIFACT_ID",
    "AUTHORITY_COMMIT",
    "AUTHORITY_PARENT",
    "CLASSIFICATION",
    "COMPACT_RESULT_SHA256",
    "CONFIG_PATH",
    "OWNER_DOCUMENT",
    "PRO20EV1PREF1CertificateError",
    "RESULT_PATH",
    "bind_live",
    "canonical_result",
    "classification_if_established",
    "emit_config_bytes",
    "expected_compact",
    "expected_config",
    "historical_runs_digest",
    "render_config",
    "store_identity",
    "validate_compact_result",
    "verify_compact",
)
