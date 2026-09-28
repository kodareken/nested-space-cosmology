"""Prospective recovery authority for the bounded TDG9 TI2 counterfactual."""

from __future__ import annotations

import ast
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import re
import stat
from typing import NoReturn

from . import tdg9_ti1_authority as ti1
from .hlt16_campaign_store import HLT16CampaignStore
from .proto19_gr0_static_factory import build_static_gr0_shells


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG9-TI2-FRZ1"
CLASSIFICATION = "premise_only_same_state_ssprk3_tableau_recovery_authority"
PROJECT_VERSION = ti1.PROJECT_VERSION
TARGET_PROTOCOL = ti1.TARGET_PROTOCOL
CONFIG_PATH = "configs/fgc/fgc-1-tdg9-ti2-frz1.toml"
RESULT_PATH = "results/fgc-1-tdg9-ti2-frz1.json"
OWNER_DOCUMENT = "docs/fgc-tdg9-ti2-frz1.md"
OUTPUT_NAMESPACE = (
    "runs/fgc-2-sf1/tdg9-ti2/ssprk3-tableau-on-rk4-2049-retries-3-5"
)
STAGING_PREFIX = ".ssprk3-tableau-on-rk4-2049-retries-3-5-ti2.stage-"
BASE_COMMIT = "8a32636654dbd6479fcb91a88d4557a8090041d3"

TI1_CONFIG_PATH = ti1.CONFIG_PATH
TI1_CONFIG_SHA256 = (
    "e36615549fa8e818285b826442a3ab6865abde0a2deef8b192f70c2457c14fa7"
)
TI1_RESULT_PATH = ti1.RESULT_PATH
TI1_RESULT_SHA256 = (
    "811942089ec67c5a1b0877c08c1fb051267f4800477084fb86c9ffc44f8aed7b"
)
TI1_RUNNER_PATH = ti1.RUNNER_PATH
TI1_RUNNER_SHA256 = (
    "abee4a9c4aa25dc8c634eaf5f6f2c54e45d74a754740fea874042eef9db8d7bb"
)
TI1_AUTHORITY_PATH = (
    "src/recursive_horizons/fgc/evolution/tdg9_ti1_authority.py"
)
TI1_AUTHORITY_SHA256 = (
    "bf657c95ad20975550edeebbdcd4b3ece462c4a8e21274f563b18d8087f3b01c"
)
TI1_OWNER_DOCUMENT = ti1.OWNER_DOCUMENT
TI1_OWNER_DOCUMENT_SHA256 = (
    "872a4b15ac3a56075f3252b91b7e084b42f7c2a7742e2ab1ab45447cd26d38d9"
)
TI1_REPRODUCER_PATH = "scripts/reproduce_fgc_tdg9_ti1_frz1.py"
TI1_REPRODUCER_SHA256 = (
    "3c79d3b11a0f79d987c33fbe9c45a27dd53374374326361fa17ee18811ae0be7"
)
TI1_INVALID_STDOUT_SHA256 = (
    "b5b2241f86e637ebf4f341d5223797b16938631cb86401f9c77d2596a8048555"
)
TI1_INVALID_STDOUT_PROVENANCE = (
    "operator_observed_not_raw_namespace_bound_or_timestamp_authenticated"
)

RUNNER_PATH = "scripts/run_fgc_tdg9_ti2.py"
RUNNER_SHA256 = (
    "488f05b65b71dd2617d1d77c2efa49496aa1a62698b37866e07bbb9c5792cd57"
)
IMPLEMENTATION_BINDINGS = ti1.IMPLEMENTATION_BINDINGS

LOC2_PREF2_CONFIG_PATH = ti1.LOC2_PREF2_CONFIG_PATH
LOC2_PREF2_CONFIG_SHA256 = ti1.LOC2_PREF2_CONFIG_SHA256
LOC2_PREF2_RESULT_PATH = ti1.LOC2_PREF2_RESULT_PATH
LOC2_PREF2_RESULT_SHA256 = ti1.LOC2_PREF2_RESULT_SHA256
LOC2_PREF2_BINDER_PATH = ti1.LOC2_PREF2_BINDER_PATH
LOC2_PREF2_BINDER_SHA256 = ti1.LOC2_PREF2_BINDER_SHA256
SEALED_STORE_LEAF_COUNT = ti1.SEALED_STORE_LEAF_COUNT
SEALED_STORE_SNAPSHOT_SHA256 = ti1.SEALED_STORE_SNAPSHOT_SHA256
MEMBER_KEY = ti1.MEMBER_KEY
MEMBER_DESCRIPTOR_SHA256 = ti1.MEMBER_DESCRIPTOR_SHA256
PHYSICAL_STATE_SHA256 = ti1.PHYSICAL_STATE_SHA256
ACCEPTED_TIME_HEX = ti1.ACCEPTED_TIME_HEX
POINT_COUNT = ti1.POINT_COUNT
OWNED_ROW_COUNT = ti1.OWNED_ROW_COUNT
COORDINATES_SHA256 = ti1.COORDINATES_SHA256
GRID_SPACING_HEX = ti1.GRID_SPACING_HEX
OUTER_RADIUS_HEX = ti1.OUTER_RADIUS_HEX
FAILED_OCCURRENCES = ti1.FAILED_OCCURRENCES
BASE_CUBIC_COUNT = ti1.BASE_CUBIC_COUNT
COMPONENT_CUBIC_COUNT = ti1.COMPONENT_CUBIC_COUNT
FAMILY_CANDIDATE_CEILING = ti1.FAMILY_CANDIDATE_CEILING
AGGREGATE_CANDIDATE_CEILING = ti1.AGGREGATE_CANDIDATE_CEILING
PRIMARY_REFINEMENT_DEPTH = ti1.PRIMARY_REFINEMENT_DEPTH
GLOBAL_PROOF_BIT_CEILING = ti1.GLOBAL_PROOF_BIT_CEILING
SEMANTICS = dict(ti1.SEMANTICS)
WORK_BUDGET = dict(ti1.WORK_BUDGET)
DECISION = dict(ti1.DECISION)
ENVIRONMENT = dict(ti1.ENVIRONMENT)
SCOPE = dict(ti1.SCOPE)

TRANSACTION_SHA256_BY_RETRY = {
    3: "7a90163d2eb4252fa7a1bbdf55d12f92a9127ed7cba5a37c28456e72d732b22b",
    4: "abab435c8222de0a86368e9562a2949578de38fe822c6032319410b49e4c6f6c",
    5: "69fd069034fdd3e2a80a5f3fa989b639b6659105a8ff3e8742e32e78739fd461",
}
REPLAYS = tuple(
    {
        **dict(item),
        "transaction_sha256": TRANSACTION_SHA256_BY_RETRY[int(item["retry"])],
    }
    for item in ti1.REPLAYS
)
FINGERPRINT_ENCODING = {
    "schema": "finite_builtin_binary64_hex_v1",
    "float_representation": "mapping_with_single_binary64_hex_key",
    "finite_builtin_float_only": True,
    "nan_and_infinity_rejected": True,
    "numpy_foreign_and_extended_numeric_types_rejected": True,
    "general_scientific_json_exact_unchanged": True,
}
RECOVERY = {
    "TI1_authority_commit": BASE_COMMIT,
    "TI1_config_path": TI1_CONFIG_PATH,
    "TI1_config_sha256": TI1_CONFIG_SHA256,
    "TI1_result_path": TI1_RESULT_PATH,
    "TI1_result_sha256": TI1_RESULT_SHA256,
    "TI1_runner_path": TI1_RUNNER_PATH,
    "TI1_runner_sha256": TI1_RUNNER_SHA256,
    "TI1_authority_path": TI1_AUTHORITY_PATH,
    "TI1_authority_sha256": TI1_AUTHORITY_SHA256,
    "TI1_owner_document": TI1_OWNER_DOCUMENT,
    "TI1_owner_document_sha256": TI1_OWNER_DOCUMENT_SHA256,
    "TI1_reproducer_path": TI1_REPRODUCER_PATH,
    "TI1_reproducer_sha256": TI1_REPRODUCER_SHA256,
    "operator_observed_invalid_stdout_sha256": TI1_INVALID_STDOUT_SHA256,
    "operator_observed_stdout_provenance": TI1_INVALID_STDOUT_PROVENANCE,
    "invalid_owner": "serialization",
    "invalid_code": "nonexact_value",
    "invalid_detail": "float",
    "first_live_float_path": "monitor.last_accepted_time",
    "first_live_float_binary64_hex": ACCEPTED_TIME_HEX,
    "failure_occurs_before_prepare_shadow": True,
    "zero_shadow_proposals_derived_from_source_order": True,
    "historical_member_fingerprint_serializer": "_json_exact",
    "historical_run_restore_line": 1591,
    "historical_run_prepare_shadow_line": 1592,
    "historical_restore_precedes_shadow": True,
    "old_output_namespace": ti1.OUTPUT_NAMESPACE,
    "old_result_not_reclassified": True,
}
CLAIMS = {
    "TI1_invalid_attempt_preserved": True,
    "TI2_execution_authorized": True,
    "TI2_result_earned": False,
    "same_state_tableau_sensitivity_measured": False,
    "state_advance_authorized": False,
    "successor_remedy_selected": False,
    "common_event_completed": False,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "mechanism_result_earned": False,
    "physical_result_earned": False,
}

DELTA_PATHS = tuple(
    sorted(
        (
            "Makefile",
            "README.md",
            CONFIG_PATH,
            "docs/claim-ledger.md",
            "docs/fgc-runtime-matrix.md",
            OWNER_DOCUMENT,
            "docs/research-roadmap.md",
            "results/README.md",
            RESULT_PATH,
            "scripts/check_repo.py",
            "scripts/reproduce_fgc_tdg9_ti2_frz1.py",
            TI1_RUNNER_PATH,
            RUNNER_PATH,
            "src/recursive_horizons/fgc/evolution/tdg9_ti2_authority.py",
            "tests/test_fgc_tdg9_ti1_runner.py",
            "tests/test_fgc_tdg9_ti2_authority.py",
            "tests/test_fgc_tdg9_ti2_runner.py",
            "tests/test_check_repo_tdg9_ti2_frz1.py",
        )
    )
)
ALLOWED_UNTRACKED = ti1.ALLOWED_UNTRACKED
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")


class TI2AuthorityError(RuntimeError):
    """The prospective TI2 image differs from its frozen recovery authority."""


def _fail(message: str) -> NoReturn:
    raise TI2AuthorityError(message)


canonical_pretty = ti1.canonical_pretty
read_leaf = ti1.read_leaf


def _git(root: Path, *arguments: str) -> bytes:
    try:
        return ti1._git(root, *arguments)  # type: ignore[attr-defined]
    except ti1.TI1AuthorityError as exc:
        raise TI2AuthorityError("Git image differs") from exc


def _paths(root: Path, *arguments: str) -> tuple[str, ...]:
    raw = _git(root, *arguments)
    if not raw:
        return ()
    chunks = raw.split(b"\0")
    if chunks[-1] != b"":
        _fail("Git path stream differs")
    return tuple(item.decode() for item in chunks[:-1])


def _working(root: Path) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    return (
        _paths(root, "diff", "--name-only", "-z", "--"),
        _paths(root, "diff", "--cached", "--name-only", "-z", "--"),
        _paths(root, "ls-files", "--others", "--exclude-standard", "-z", "--"),
    )


def _expected_config() -> dict[str, object]:
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL,
        "owner_document": OWNER_DOCUMENT,
        "output_namespace": OUTPUT_NAMESPACE,
        "recovery": RECOVERY,
        "implementation": {
            "runner_path": RUNNER_PATH,
            "runner_sha256": RUNNER_SHA256,
            **{
                key: value
                for binding_key, path, digest in IMPLEMENTATION_BINDINGS
                for key, value in (
                    (binding_key, path),
                    (binding_key.replace("_path", "_sha256"), digest),
                )
            },
        },
        "predecessor": {
            "LOC2_PREF2_config_path": LOC2_PREF2_CONFIG_PATH,
            "LOC2_PREF2_config_sha256": LOC2_PREF2_CONFIG_SHA256,
            "LOC2_PREF2_result_path": LOC2_PREF2_RESULT_PATH,
            "LOC2_PREF2_result_sha256": LOC2_PREF2_RESULT_SHA256,
            "LOC2_PREF2_binder_path": LOC2_PREF2_BINDER_PATH,
            "LOC2_PREF2_binder_sha256": LOC2_PREF2_BINDER_SHA256,
            "sealed_store_leaf_count": SEALED_STORE_LEAF_COUNT,
            "sealed_store_snapshot_sha256": SEALED_STORE_SNAPSHOT_SHA256,
        },
        "fingerprint_encoding": FINGERPRINT_ENCODING,
        "semantics": SEMANTICS,
        "selection": {
            "member_key": MEMBER_KEY,
            "descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
            "physical_state_sha256": PHYSICAL_STATE_SHA256,
            "accepted_time_hex": ACCEPTED_TIME_HEX,
            "point_count": POINT_COUNT,
            "owned_row_count": OWNED_ROW_COUNT,
            "coordinates_sha256": COORDINATES_SHA256,
            "grid_spacing_hex": GRID_SPACING_HEX,
            "outer_radius_hex": OUTER_RADIUS_HEX,
            "retries": [3, 4, 5],
            "failed_occurrences": [
                {"retry": retry, "channel": channel}
                for retry, channel in FAILED_OCCURRENCES
            ],
        },
        "replays": [dict(item) for item in REPLAYS],
        "work_budget": WORK_BUDGET,
        "decision": DECISION,
        "environment": ENVIRONMENT,
        "scope": SCOPE,
        "claims": CLAIMS,
    }


def parse_config(raw: bytes) -> dict[str, object]:
    import tomllib

    try:
        value = tomllib.loads(raw.decode())
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise TI2AuthorityError("invalid TI2 TOML") from exc
    if value != _expected_config():
        _fail("TI2 config differs from frozen contract")
    return value


def require_output_absent(root: Path) -> None:
    target = root / OUTPUT_NAMESPACE
    current = root
    for part in target.relative_to(root).parent.parts:
        current /= part
        try:
            metadata = current.lstat()
        except FileNotFoundError:
            return
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
            _fail("TI2 output parent is unsafe")
    try:
        target.lstat()
    except FileNotFoundError:
        pass
    else:
        _fail("TI2 output namespace already exists")
    if any(item.name.startswith(STAGING_PREFIX) for item in target.parent.iterdir()):
        _fail("TI2 staging namespace already exists")


def _verify_historical_bindings(root: Path) -> None:
    for path, digest in (
        (TI1_CONFIG_PATH, TI1_CONFIG_SHA256),
        (TI1_RESULT_PATH, TI1_RESULT_SHA256),
        (TI1_RUNNER_PATH, TI1_RUNNER_SHA256),
        (TI1_AUTHORITY_PATH, TI1_AUTHORITY_SHA256),
        (TI1_OWNER_DOCUMENT, TI1_OWNER_DOCUMENT_SHA256),
        (TI1_REPRODUCER_PATH, TI1_REPRODUCER_SHA256),
    ):
        if sha256(_git(root, "show", f"{BASE_COMMIT}:{path}")).hexdigest() != digest:
            _fail(f"historical TI1 binding differs: {path}")


def _verify_historical_source_order(root: Path) -> dict[str, object]:
    """Prove the bound TI1 failure site precedes any shadow proposal."""

    raw = _git(root, "show", f"{BASE_COMMIT}:{TI1_RUNNER_PATH}")
    try:
        tree = ast.parse(raw.decode())
    except (UnicodeDecodeError, SyntaxError) as exc:
        raise TI2AuthorityError("historical TI1 runner is not parseable") from exc
    functions = {
        node.name: node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    try:
        fingerprint = functions["_member_fingerprint"]
        restore = functions["_restore_replay"]
        run = functions["run"]
    except KeyError as exc:
        raise TI2AuthorityError("historical TI1 source owner differs") from exc

    def call_name(node: ast.Call) -> str | None:
        return node.func.id if isinstance(node.func, ast.Name) else None

    fingerprint_calls = {
        call_name(node)
        for node in ast.walk(fingerprint)
        if isinstance(node, ast.Call)
    }
    restore_calls = {
        call_name(node)
        for node in ast.walk(restore)
        if isinstance(node, ast.Call)
    }
    run_calls = sorted(
        (
            (node.lineno, call_name(node))
            for node in ast.walk(run)
            if isinstance(node, ast.Call)
            and call_name(node) in {"_restore_replay", "_prepare_shadow"}
        ),
        key=lambda item: item[0],
    )
    expected_calls = [(1591, "_restore_replay"), (1592, "_prepare_shadow")]
    if (
        "_json_exact" not in fingerprint_calls
        or "_fingerprint_exact" in fingerprint_calls
        or "_member_fingerprint" not in restore_calls
        or run_calls != expected_calls
    ):
        _fail("historical TI1 source-order proof differs")
    return {
        "member_fingerprint_serializer": "_json_exact",
        "restore_calls_member_fingerprint": True,
        "run_restore_line": 1591,
        "run_prepare_shadow_line": 1592,
        "restore_precedes_shadow": True,
        "shadow_proposals_before_serialization_failure": 0,
    }


def _verify_live_bindings(root: Path, *, commit: str | None) -> None:
    bindings = (
        (RUNNER_PATH, RUNNER_SHA256),
        *((path, digest) for _key, path, digest in IMPLEMENTATION_BINDINGS),
    )
    for path, digest in bindings:
        raw = read_leaf(root, path) if commit is None else _git(root, "show", f"{commit}:{path}")
        if sha256(raw).hexdigest() != digest:
            _fail(f"TI2 implementation binding differs: {path}")


def _snapshot_store(root: Path) -> tuple[int, str]:
    from scripts import run_fgc_tdg9_loc1 as loc1

    try:
        return loc1._snapshot_store(root)
    except Exception as exc:
        raise TI2AuthorityError("sealed store cannot be authenticated") from exc


def real_predecessor_fingerprints(root: Path) -> tuple[dict[str, object], ...]:
    """Restore all three exact predecessors and fingerprint without proposals."""

    from scripts import run_fgc_tdg9_ti2 as runner

    from . import tdg9_ar1_authority as ar1

    repository = root.resolve()
    store = HLT16CampaignStore(repository / ar1.PREF2_STORE_PATH)
    shells = build_static_gr0_shells(repository)
    receipts: list[dict[str, object]] = []
    for replay in REPLAYS:
        try:
            restored = runner._restore_replay(
                repository, store, shells, replay
            )
            fingerprint = restored.fingerprint
        except Exception as exc:
            raise TI2AuthorityError(
                f"retry-{replay['retry']} predecessor fingerprint failed"
            ) from exc
        retry = int(replay["retry"])
        if (
            fingerprint["transaction_sha256"]
            != TRANSACTION_SHA256_BY_RETRY[retry]
            or fingerprint["state_sha256"] != PHYSICAL_STATE_SHA256
            or fingerprint["accepted_time_hex"] != ACCEPTED_TIME_HEX
            or fingerprint["descriptor_sha256"] != MEMBER_DESCRIPTOR_SHA256
        ):
            _fail(f"retry-{retry} predecessor fingerprint differs")
        receipts.append(
            {
                "retry": retry,
                "transaction_sha256": fingerprint["transaction_sha256"],
                "state_sha256": fingerprint["state_sha256"],
                "accepted_time_hex": fingerprint["accepted_time_hex"],
                "descriptor_sha256": fingerprint["descriptor_sha256"],
                "shadow_proposal_constructed": False,
            }
        )
    return tuple(receipts)


def _expected_fingerprint_receipts() -> list[dict[str, object]]:
    return [
        {
            "retry": retry,
            "transaction_sha256": TRANSACTION_SHA256_BY_RETRY[retry],
            "state_sha256": PHYSICAL_STATE_SHA256,
            "accepted_time_hex": ACCEPTED_TIME_HEX,
            "descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
            "shadow_proposal_constructed": False,
        }
        for retry in (3, 4, 5)
    ]


def build_prelaunch(
    config_raw: bytes,
    root: Path,
    *,
    store_snapshot: tuple[int, str],
) -> dict[str, object]:
    repository = root.resolve()
    parse_config(config_raw)
    ti1._environment()  # type: ignore[attr-defined]
    head = _git(repository, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    if head != BASE_COMMIT:
        _fail("prelaunch HEAD is not the sealed invalid TI1 authority")
    unstaged, staged, untracked = _working(repository)
    observed = tuple(
        sorted((*unstaged, *(item for item in untracked if item not in ALLOWED_UNTRACKED)))
    )
    if staged or observed != DELTA_PATHS:
        _fail("prospective TI2 delta differs")
    for relative in DELTA_PATHS:
        read_leaf(repository, relative)
    _verify_historical_bindings(repository)
    source_order = _verify_historical_source_order(repository)
    _verify_live_bindings(repository, commit=None)
    if store_snapshot != (SEALED_STORE_LEAF_COUNT, SEALED_STORE_SNAPSHOT_SHA256):
        _fail("sealed store differs at TI2 prelaunch")
    ti1.require_output_absent(repository)
    require_output_absent(repository)
    fingerprints = list(real_predecessor_fingerprints(repository))
    if fingerprints != _expected_fingerprint_receipts():
        _fail("restored predecessor fingerprints differ")
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "gate_status": "pass",
        "source_config_sha256": sha256(config_raw).hexdigest(),
        "artifact_payload": {
            **_expected_config(),
            "store_identity_observed_at_prelaunch": {
                "leaf_count": store_snapshot[0],
                "sha256": store_snapshot[1],
            },
            "TI1_output_namespace_absent_after_invalid_attempt": True,
            "TI2_output_namespace_absent_at_prelaunch": True,
            "TI1_source_order_proof": source_order,
            "restored_predecessor_fingerprints": fingerprints,
            "shadow_executed": False,
        },
    }


def _unique(items: list[tuple[str, object]]) -> dict[str, object]:
    answer: dict[str, object] = {}
    for key, value in items:
        if key in answer:
            _fail("duplicate compact JSON key")
        answer[key] = value
    return answer


def validate_compact(config_raw: bytes, result_raw: bytes) -> dict[str, object]:
    parse_config(config_raw)
    try:
        result = json.loads(result_raw, object_pairs_hook=_unique)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TI2AuthorityError("invalid TI2 compact JSON") from exc
    expected = {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "gate_status": "pass",
        "source_config_sha256": sha256(config_raw).hexdigest(),
        "artifact_payload": {
            **_expected_config(),
            "store_identity_observed_at_prelaunch": {
                "leaf_count": SEALED_STORE_LEAF_COUNT,
                "sha256": SEALED_STORE_SNAPSHOT_SHA256,
            },
            "TI1_output_namespace_absent_after_invalid_attempt": True,
            "TI2_output_namespace_absent_at_prelaunch": True,
            "TI1_source_order_proof": {
                "member_fingerprint_serializer": "_json_exact",
                "restore_calls_member_fingerprint": True,
                "run_restore_line": 1591,
                "run_prepare_shadow_line": 1592,
                "restore_precedes_shadow": True,
                "shadow_proposals_before_serialization_failure": 0,
            },
            "restored_predecessor_fingerprints": _expected_fingerprint_receipts(),
            "shadow_executed": False,
        },
    }
    if result != expected or result_raw != canonical_pretty(result):
        _fail("TI2 compact result differs")
    return result


@dataclass(frozen=True, slots=True)
class TI2Authority:
    authority_commit: str
    execution_authorized: bool = True
    tableau_selector: str = "SSPRK3"
    actual_spatial_operator: str = "inherited_RK4_2049_SBP4"
    state_advance_authorized: bool = False
    production_comparator_authorized: bool = False


def authorize(
    root: Path,
    config_raw: bytes,
    result_raw: bytes,
    authority_commit: str,
) -> TI2Authority:
    repository = root.resolve()
    if not _COMMIT.fullmatch(authority_commit):
        _fail("TI2 authority commit is malformed")
    validate_compact(config_raw, result_raw)
    ti1._environment()  # type: ignore[attr-defined]
    head = _git(repository, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    parents = _git(repository, "rev-list", "--parents", "-n", "1", authority_commit).decode().split()[1:]
    if head != authority_commit or parents != [BASE_COMMIT]:
        _fail("TI2 authority is not the exact direct successor")
    unstaged, staged, untracked = _working(repository)
    if unstaged or staged or any(item not in ALLOWED_UNTRACKED for item in untracked):
        _fail("TI2 authority working image is not clean")
    delta = tuple(
        sorted(
            _paths(
                repository,
                "diff",
                "--name-only",
                "-z",
                BASE_COMMIT,
                authority_commit,
                "--",
            )
        )
    )
    if delta != DELTA_PATHS:
        _fail("committed TI2 delta differs")
    if _git(repository, "show", f"{authority_commit}:{CONFIG_PATH}") != config_raw:
        _fail("committed TI2 config differs")
    if _git(repository, "show", f"{authority_commit}:{RESULT_PATH}") != result_raw:
        _fail("committed TI2 compact result differs")
    _verify_historical_bindings(repository)
    _verify_historical_source_order(repository)
    _verify_live_bindings(repository, commit=authority_commit)
    ti1.require_output_absent(repository)
    if _snapshot_store(repository) != (
        SEALED_STORE_LEAF_COUNT,
        SEALED_STORE_SNAPSHOT_SHA256,
    ):
        _fail("sealed store differs before TI2 authorization")
    if list(real_predecessor_fingerprints(repository)) != _expected_fingerprint_receipts():
        _fail("restored predecessor fingerprints differ before authorization")
    return TI2Authority(authority_commit=authority_commit)


__all__ = [
    "ARTIFACT_ID",
    "CONFIG_PATH",
    "DELTA_PATHS",
    "OUTPUT_NAMESPACE",
    "RESULT_PATH",
    "TI2Authority",
    "TI2AuthorityError",
    "authorize",
    "build_prelaunch",
    "canonical_pretty",
    "parse_config",
    "read_leaf",
    "real_predecessor_fingerprints",
    "require_output_absent",
    "validate_compact",
]
