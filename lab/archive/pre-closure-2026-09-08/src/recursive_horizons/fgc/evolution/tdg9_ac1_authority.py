"""Prospective authority for one retry-3 SSPRK3-on-SBP4 all-channel audit."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import re
import stat
from typing import NoReturn

from . import tdg9_ti2_authority as ti2
from . import tdg9_ti2_pref1_binder as pref1
from .hlt16_campaign_store import HLT16CampaignStore
from .proto19_gr0_static_factory import build_static_gr0_shells
from .tdg6_temporal_admission_design import TDG6_COMPLETE_STATE_CHANNELS


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG9-AC1-FRZ1"
CLASSIFICATION = (
    "premise_only_retry3_ssprk3_on_sbp4_all_channel_production_classifier_audit"
)
PROJECT_VERSION = ti2.PROJECT_VERSION
TARGET_PROTOCOL = ti2.TARGET_PROTOCOL
CONFIG_PATH = "configs/fgc/fgc-1-tdg9-ac1-frz1.toml"
RESULT_PATH = "results/fgc-1-tdg9-ac1-frz1.json"
OWNER_DOCUMENT = "docs/fgc-tdg9-ac1-frz1.md"
OUTPUT_NAMESPACE = (
    "runs/fgc-2-sf1/tdg9-ac1/ssprk3-sbp4-retry3-all-channels"
)
STAGING_PREFIX = ".ssprk3-sbp4-retry3-all-channels-ac1.stage-"
BASE_COMMIT = "660f369365ab8417f4c68eded46acb9bb2b05a5b"

TI2_AUTHORITY_COMMIT = pref1.AUTHORITY_COMMIT
TI2_CONFIG_PATH = ti2.CONFIG_PATH
TI2_CONFIG_SHA256 = (
    "8631ce20e22a56b00f4a9580df1ec2e2e98176b14c4761dc78e316cb4a771546"
)
TI2_RESULT_PATH = ti2.RESULT_PATH
TI2_RESULT_SHA256 = (
    "2f2f2f844e9611b1023cee869e9f01ce33cfdef2df01cf488302eaaa28a1e3f5"
)
TI2_AUTHORITY_PATH = (
    "src/recursive_horizons/fgc/evolution/tdg9_ti2_authority.py"
)
TI2_AUTHORITY_SHA256 = (
    "7715685100e56434a36edbec73ca684fe2fd88010af3eaeae1a59761884512b5"
)
TI2_RUNNER_PATH = ti2.RUNNER_PATH
TI2_RUNNER_SHA256 = (
    "488f05b65b71dd2617d1d77c2efa49496aa1a62698b37866e07bbb9c5792cd57"
)
TI2_OWNER_DOCUMENT = ti2.OWNER_DOCUMENT
TI2_OWNER_DOCUMENT_SHA256 = (
    "91f0f9fbf410e6e896d12c5bf01d0c7105f6681114caded73ad21c0a1f668080"
)
TI2_REPRODUCER_PATH = "scripts/reproduce_fgc_tdg9_ti2_frz1.py"
TI2_REPRODUCER_SHA256 = (
    "42bef89656720404668e3d4daab60295d7b06b9a85d79410bc197de41ece97c4"
)
TI2_RAW_NAMESPACE = pref1.RAW_NAMESPACE
TI2_RAW_SCHEMA = pref1.RAW_SCHEMA
TI2_RAW_MANIFEST_SHA256 = pref1.RAW_MANIFEST_SHA256
TI2_RAW_TERMINAL_SHA256 = pref1.RAW_TERMINAL_SHA256
TI2_RAW_CLASSIFICATION = pref1.RAW_CLASSIFICATION

PREF1_CONFIG_PATH = pref1.CONFIG_PATH
PREF1_CONFIG_SHA256 = (
    "7023d7d1f528e71914d0cacfee6f8074d434fc1210378f8aeb94b0b4f85613c2"
)
PREF1_RESULT_PATH = pref1.RESULT_PATH
PREF1_RESULT_SHA256 = pref1.COMPACT_RESULT_SHA256
PREF1_BINDER_PATH = (
    "src/recursive_horizons/fgc/evolution/tdg9_ti2_pref1_binder.py"
)
PREF1_BINDER_SHA256 = (
    "89e8242ed94ff93619a3d3530cb3e8d6d282fe038090be02cec07c9ad501b653"
)
PREF1_OWNER_DOCUMENT = pref1.OWNER_DOCUMENT
PREF1_OWNER_DOCUMENT_SHA256 = (
    "60718f49e7d0b50b4b9ca24f6a239e9b6d7fe99c3d84e8daabde8f4982b2a37f"
)
PREF1_REPRODUCER_PATH = "scripts/reproduce_fgc_tdg9_ti2_pref1.py"
PREF1_REPRODUCER_SHA256 = (
    "b78abc49f38181eb8da08199c2aeedb0cbca1522258f7d643068fae5e7de75c7"
)

RUNNER_PATH = "scripts/run_fgc_tdg9_ac1.py"
RUNNER_SHA256 = (
    "bc2f28e99a747caabac3f419b05574bfe06fcddc5e7b5f1d21879e3639e77605"
)
IMPLEMENTATION_BINDINGS = (
    (
        "tableau_runtime_path",
        "src/recursive_horizons/fgc/evolution/numerical_engine.py",
        "8ea99604a85b4e1acbf869ed2462ea4aa8301c06e068330a1d0ce66a51a0bebf",
    ),
    (
        "compositor_path",
        "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_runtime.py",
        "039de0bd008a2536bea6f7187afc3856bfc841d14790cbae5ca2fe27f7eb8e2d",
    ),
    (
        "classifier_path",
        "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_design.py",
        "87d9a90e272a89479690c0d432f2d5cd979d88841615aaaeb9b49a2c0cca651e",
    ),
    (
        "TI2_runner_path",
        TI2_RUNNER_PATH,
        TI2_RUNNER_SHA256,
    ),
    (
        "TI2_authority_path",
        TI2_AUTHORITY_PATH,
        TI2_AUTHORITY_SHA256,
    ),
    (
        "PREF1_binder_path",
        PREF1_BINDER_PATH,
        PREF1_BINDER_SHA256,
    ),
)

SEALED_STORE_LEAF_COUNT = ti2.SEALED_STORE_LEAF_COUNT
SEALED_STORE_SNAPSHOT_SHA256 = ti2.SEALED_STORE_SNAPSHOT_SHA256
MEMBER_KEY = ti2.MEMBER_KEY
MEMBER_DESCRIPTOR_SHA256 = ti2.MEMBER_DESCRIPTOR_SHA256
PHYSICAL_STATE_SHA256 = ti2.PHYSICAL_STATE_SHA256
ACCEPTED_TIME_HEX = ti2.ACCEPTED_TIME_HEX
POINT_COUNT = ti2.POINT_COUNT
OWNED_ROW_COUNT = ti2.OWNED_ROW_COUNT
COORDINATES_SHA256 = ti2.COORDINATES_SHA256
GRID_SPACING_HEX = ti2.GRID_SPACING_HEX
OUTER_RADIUS_HEX = ti2.OUTER_RADIUS_HEX
CHANNEL_ORDER = TDG6_COMPLETE_STATE_CHANNELS
RETRY = 3
WIDTH_HEX = "0x1.aaa9612df8000p-11"
TRANSACTION_SHA256 = ti2.TRANSACTION_SHA256_BY_RETRY[RETRY]
REPLAY = dict(ti2.REPLAYS[0])
HISTORICAL_JOURNAL_SHA256 = (
    "0d21f650ccc46db778fd81fbebf5acab6394e1c2f7972b7940e24c76af9e1f59"
)
TI2_SAMPLED_RETRY3_CHANNELS = ("u:alpha", "u:R")
ENVIRONMENT = dict(ti2.ENVIRONMENT)

SEMANTICS = {
    "experiment_label": (
        "retry3_SSPRK3_tableau_on_inherited_SBP4_all_18_channel_"
        "production_classifier_audit"
    ),
    "tableau_selector": "SSPRK3",
    "tableau_runtime_selector": ti2.SEMANTICS["tableau_runtime_selector"],
    "source_member_method": "RK4",
    "source_member_integrator": ti2.SEMANTICS["source_member_integrator"],
    "actual_spatial_operator": "inherited_RK4_2049_SBP4",
    "RHS_owner": ti2.SEMANTICS["RHS_owner"],
    "projector_owner": ti2.SEMANTICS["projector_owner"],
    "transaction_owner": ti2.SEMANTICS["transaction_owner"],
    "tracer_owner": ti2.SEMANTICS["tracer_owner"],
    "ledger_owner": ti2.SEMANTICS["ledger_owner"],
    "only_changed_variable": "explicit_time_tableau_selector",
    "channel_order_owner": "TDG6_COMPLETE_STATE_CHANNELS",
    "classifier_owner": "unchanged_binary64_TDG6_production_classifier",
    "production_SSPRK3_comparator": False,
    "independent_method_agreement": False,
}
WORK_BUDGET = {
    "retry_count": 1,
    "retry": RETRY,
    "shadow_paths": 7,
    "shadow_proposals": 7,
    "SSPRK3_records_per_proposal": 4,
    "maximum_stage_and_endpoint_RHS_records": 28,
    "channel_count": len(CHANNEL_ORDER),
    "fourth_width_authorized": False,
    "retry_4_authorized": False,
    "retry_5_authorized": False,
    "FAILED_OCCURRENCES_filter_applied": False,
    "LOC1_exact_localization_authorized": False,
    "LOC2_exact_localization_authorized": False,
    "resource_escalation_authorized": False,
}
DECISION = {
    "unchanged_contraction_predicate": "D12_squared_times_8_le_D01_squared",
    "minimum_observed_order": "3/2",
    "complete_admission_is_all_of": True,
    "order_inconclusive_is_channel_failure": True,
    "complete_pass_iff_failed_channels_empty": True,
    "complete_pass_iff_all_channel_admission_flags_true": True,
    "debit_equals_finest_upper": True,
    "successor_remedy_selected": False,
    "production_method_earned_on_all_pass": False,
}
SCOPE = {
    "live_store_read_only": True,
    "compact_verifier_store_blind": True,
    "compact_verifier_raw_blind": True,
    "compact_verifier_shadow_blind": True,
    "prelaunch_may_restore_retry_3": True,
    "prelaunch_constructs_zero_shadow_proposals": True,
    "same_state_tableau_counterfactual_only": True,
    "historical_threshold_changed": False,
    "ULP_used_as_tolerance": False,
    "production_SSPRK3_comparator": False,
    "independent_method_agreement": False,
    "campaign_store_mutation_authorized": False,
    "PDE_state_commit_authorized": False,
    "temporal_retry_admission_authorized": False,
    "fine_path_commit_authorized": False,
    "continuation_authorized": False,
    "GR0_calibration_authorized": False,
    "candidate_branches_authorized": False,
    "mechanism_result_authorized": False,
    "physical_result_authorized": False,
}
CLAIMS = {
    "TI2_PREF1_bound": True,
    "AC1_execution_authorized": True,
    "AC1_result_earned": False,
    "production_method_earned": False,
    "production_SSPRK3_comparator": False,
    "independent_method_agreement": False,
    "state_advance_authorized": False,
    "successor_remedy_selected": False,
    "common_event_completed": False,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "mechanism_result_earned": False,
    "physical_result_earned": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
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
            "scripts/reproduce_fgc_tdg9_ac1_frz1.py",
            RUNNER_PATH,
            "src/recursive_horizons/fgc/evolution/tdg9_ac1_authority.py",
            "tests/test_check_repo_tdg9_ac1_frz1.py",
            "tests/test_fgc_tdg9_ac1_authority.py",
            "tests/test_fgc_tdg9_ac1_runner.py",
        )
    )
)
ALLOWED_UNTRACKED = ti2.ALLOWED_UNTRACKED
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")


class AC1AuthorityError(RuntimeError):
    """The prospective AC1 image differs from its frozen audit authority."""


def _fail(message: str) -> NoReturn:
    raise AC1AuthorityError(message)


canonical_pretty = ti2.canonical_pretty
read_leaf = ti2.read_leaf


def _git(root: Path, *arguments: str) -> bytes:
    try:
        return ti2._git(root, *arguments)  # type: ignore[attr-defined]
    except ti2.TI2AuthorityError as exc:
        raise AC1AuthorityError("Git image differs") from exc


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


def _require_no_replace_refs(root: Path) -> None:
    raw = _git(root, "for-each-ref", "--format=%(refname)", "refs/replace")
    if raw.strip():
        _fail("Git replace refs are present")


def _implementation() -> dict[str, object]:
    answer: dict[str, object] = {
        "runner_path": RUNNER_PATH,
        "runner_sha256": RUNNER_SHA256,
    }
    for key, path, digest in IMPLEMENTATION_BINDINGS:
        answer[key] = path
        answer[key.replace("_path", "_sha256")] = digest
    return answer


def _expected_config() -> dict[str, object]:
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL,
        "owner_document": OWNER_DOCUMENT,
        "output_namespace": OUTPUT_NAMESPACE,
        "predecessor": {
            "commit": BASE_COMMIT,
            "TI2_authority_commit": TI2_AUTHORITY_COMMIT,
            "TI2_config_path": TI2_CONFIG_PATH,
            "TI2_config_sha256": TI2_CONFIG_SHA256,
            "TI2_result_path": TI2_RESULT_PATH,
            "TI2_result_sha256": TI2_RESULT_SHA256,
            "TI2_authority_path": TI2_AUTHORITY_PATH,
            "TI2_authority_sha256": TI2_AUTHORITY_SHA256,
            "TI2_runner_path": TI2_RUNNER_PATH,
            "TI2_runner_sha256": TI2_RUNNER_SHA256,
            "TI2_owner_document": TI2_OWNER_DOCUMENT,
            "TI2_owner_document_sha256": TI2_OWNER_DOCUMENT_SHA256,
            "TI2_reproducer_path": TI2_REPRODUCER_PATH,
            "TI2_reproducer_sha256": TI2_REPRODUCER_SHA256,
            "TI2_raw_namespace": TI2_RAW_NAMESPACE,
            "TI2_raw_schema": TI2_RAW_SCHEMA,
            "TI2_raw_manifest_sha256": TI2_RAW_MANIFEST_SHA256,
            "TI2_raw_terminal_sha256": TI2_RAW_TERMINAL_SHA256,
            "TI2_raw_classification": TI2_RAW_CLASSIFICATION,
            "PREF1_config_path": PREF1_CONFIG_PATH,
            "PREF1_config_sha256": PREF1_CONFIG_SHA256,
            "PREF1_result_path": PREF1_RESULT_PATH,
            "PREF1_result_sha256": PREF1_RESULT_SHA256,
            "PREF1_binder_path": PREF1_BINDER_PATH,
            "PREF1_binder_sha256": PREF1_BINDER_SHA256,
            "PREF1_owner_document": PREF1_OWNER_DOCUMENT,
            "PREF1_owner_document_sha256": PREF1_OWNER_DOCUMENT_SHA256,
            "PREF1_reproducer_path": PREF1_REPRODUCER_PATH,
            "PREF1_reproducer_sha256": PREF1_REPRODUCER_SHA256,
            "sealed_store_leaf_count": SEALED_STORE_LEAF_COUNT,
            "sealed_store_snapshot_sha256": SEALED_STORE_SNAPSHOT_SHA256,
        },
        "implementation": _implementation(),
        "fingerprint_encoding": dict(ti2.FINGERPRINT_ENCODING),
        "semantics": dict(SEMANTICS),
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
            "retry": RETRY,
            "attempted_width_hex": WIDTH_HEX,
            "channel_order": list(CHANNEL_ORDER),
            "TI2_sampled_retry3_channels": list(TI2_SAMPLED_RETRY3_CHANNELS),
        },
        "replay": dict(REPLAY),
        "work_budget": dict(WORK_BUDGET),
        "decision": dict(DECISION),
        "environment": dict(ENVIRONMENT),
        "scope": dict(SCOPE),
        "claims": dict(CLAIMS),
    }


def parse_config(raw: bytes) -> dict[str, object]:
    import tomllib

    try:
        value = tomllib.loads(raw.decode())
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise AC1AuthorityError("invalid AC1 TOML") from exc
    if value != _expected_config():
        _fail("AC1 config differs from frozen contract")
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
            _fail("AC1 output parent is unsafe")
    try:
        target.lstat()
    except FileNotFoundError:
        pass
    else:
        _fail("AC1 output namespace already exists")
    if any(item.name.startswith(STAGING_PREFIX) for item in target.parent.iterdir()):
        _fail("AC1 staging namespace already exists")


def _verify_predecessor_bindings(root: Path) -> None:
    for path, digest in (
        (TI2_CONFIG_PATH, TI2_CONFIG_SHA256),
        (TI2_RESULT_PATH, TI2_RESULT_SHA256),
        (TI2_AUTHORITY_PATH, TI2_AUTHORITY_SHA256),
        (TI2_RUNNER_PATH, TI2_RUNNER_SHA256),
        (TI2_OWNER_DOCUMENT, TI2_OWNER_DOCUMENT_SHA256),
        (TI2_REPRODUCER_PATH, TI2_REPRODUCER_SHA256),
        (PREF1_CONFIG_PATH, PREF1_CONFIG_SHA256),
        (PREF1_RESULT_PATH, PREF1_RESULT_SHA256),
        (PREF1_BINDER_PATH, PREF1_BINDER_SHA256),
        (PREF1_OWNER_DOCUMENT, PREF1_OWNER_DOCUMENT_SHA256),
        (PREF1_REPRODUCER_PATH, PREF1_REPRODUCER_SHA256),
    ):
        if sha256(_git(root, "show", f"{BASE_COMMIT}:{path}")).hexdigest() != digest:
            _fail(f"historical predecessor binding differs: {path}")


def _verify_live_bindings(root: Path, *, commit: str | None) -> None:
    bindings = (
        (RUNNER_PATH, RUNNER_SHA256),
        *((path, digest) for _key, path, digest in IMPLEMENTATION_BINDINGS),
    )
    for path, digest in bindings:
        raw = read_leaf(root, path) if commit is None else _git(root, "show", f"{commit}:{path}")
        if sha256(raw).hexdigest() != digest:
            _fail(f"AC1 implementation binding differs: {path}")


def _snapshot_store(root: Path) -> tuple[int, str]:
    try:
        return ti2._snapshot_store(root)  # type: ignore[attr-defined]
    except ti2.TI2AuthorityError as exc:
        raise AC1AuthorityError("sealed store cannot be authenticated") from exc


def _expected_fingerprint_receipt() -> dict[str, object]:
    return {
        "retry": RETRY,
        "transaction_sha256": TRANSACTION_SHA256,
        "state_sha256": PHYSICAL_STATE_SHA256,
        "accepted_time_hex": ACCEPTED_TIME_HEX,
        "descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
        "shadow_proposal_constructed": False,
    }


def expected_raw_replay_receipt() -> dict[str, object]:
    return {
        "member_key": MEMBER_KEY,
        "retry": RETRY,
        "attempted_width_hex": WIDTH_HEX,
        "accepted_time_hex": ACCEPTED_TIME_HEX,
        "state_sha256": PHYSICAL_STATE_SHA256,
        "descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
        "transaction_sha256": TRANSACTION_SHA256,
        "historical_journal_sha256": HISTORICAL_JOURNAL_SHA256,
    }


def real_predecessor_fingerprints(root: Path) -> tuple[dict[str, object], ...]:
    """Restore only retry 3 and fingerprint it without constructing a proposal."""

    from scripts import run_fgc_tdg9_ti2 as runner

    from . import tdg9_ar1_authority as ar1

    repository = root.resolve()
    store = HLT16CampaignStore(repository / ar1.PREF2_STORE_PATH)
    shells = build_static_gr0_shells(repository)
    try:
        restored = runner._restore_replay(repository, store, shells, REPLAY)
        fingerprint = restored.fingerprint
    except Exception as exc:
        raise AC1AuthorityError("retry-3 predecessor fingerprint failed") from exc
    if (
        fingerprint["transaction_sha256"] != TRANSACTION_SHA256
        or fingerprint["state_sha256"] != PHYSICAL_STATE_SHA256
        or fingerprint["accepted_time_hex"] != ACCEPTED_TIME_HEX
        or fingerprint["descriptor_sha256"] != MEMBER_DESCRIPTOR_SHA256
    ):
        _fail("retry-3 predecessor fingerprint differs")
    return (
        {
            "retry": RETRY,
            "transaction_sha256": fingerprint["transaction_sha256"],
            "state_sha256": fingerprint["state_sha256"],
            "accepted_time_hex": fingerprint["accepted_time_hex"],
            "descriptor_sha256": fingerprint["descriptor_sha256"],
            "shadow_proposal_constructed": False,
        },
    )


def expected_compact(config_raw: bytes) -> dict[str, object]:
    parse_config(config_raw)
    return {
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
            "AC1_output_namespace_absent_at_prelaunch": True,
            "replace_refs_absent": True,
            "restored_predecessor_fingerprints": [_expected_fingerprint_receipt()],
            "shadow_executed": False,
            "shadow_proposals_constructed": 0,
        },
    }


def build_prelaunch(
    config_raw: bytes,
    root: Path,
    *,
    store_snapshot: tuple[int, str],
) -> dict[str, object]:
    repository = root.resolve()
    parse_config(config_raw)
    ti2.ti1._environment()  # type: ignore[attr-defined]
    head = _git(repository, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    if head != BASE_COMMIT:
        _fail("prelaunch HEAD is not the sealed TI2-PREF1 commit")
    _require_no_replace_refs(repository)
    unstaged, staged, untracked = _working(repository)
    observed = tuple(
        sorted((*unstaged, *(item for item in untracked if item not in ALLOWED_UNTRACKED)))
    )
    if staged or observed != DELTA_PATHS:
        _fail("prospective AC1 delta differs")
    for relative in DELTA_PATHS:
        read_leaf(repository, relative)
    _verify_predecessor_bindings(repository)
    _verify_live_bindings(repository, commit=None)
    if store_snapshot != (SEALED_STORE_LEAF_COUNT, SEALED_STORE_SNAPSHOT_SHA256):
        _fail("sealed store differs at AC1 prelaunch")
    require_output_absent(repository)
    fingerprints = list(real_predecessor_fingerprints(repository))
    if fingerprints != [_expected_fingerprint_receipt()]:
        _fail("restored predecessor fingerprints differ")
    return {
        **expected_compact(config_raw),
        "artifact_payload": {
            **expected_compact(config_raw)["artifact_payload"],
            "store_identity_observed_at_prelaunch": {
                "leaf_count": store_snapshot[0],
                "sha256": store_snapshot[1],
            },
            "restored_predecessor_fingerprints": fingerprints,
            "shadow_executed": False,
            "shadow_proposals_constructed": 0,
        },
    }


def _unique(items: list[tuple[str, object]]) -> dict[str, object]:
    answer: dict[str, object] = {}
    for key, value in items:
        if key in answer:
            _fail("duplicate compact JSON key")
        answer[key] = value
    return answer


def _reject_compact_json_constant(token: str) -> NoReturn:
    _fail(f"nonfinite AC1 compact JSON constant: {token}")


def validate_compact(config_raw: bytes, result_raw: bytes) -> dict[str, object]:
    expected = expected_compact(config_raw)
    try:
        result = json.loads(
            result_raw,
            object_pairs_hook=_unique,
            parse_constant=_reject_compact_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise AC1AuthorityError("invalid AC1 compact JSON") from exc
    if result != expected or result_raw != canonical_pretty(result):
        _fail("AC1 compact result differs")
    payload = result["artifact_payload"]
    if payload["claims"]["AC1_result_earned"] is not False:
        _fail("AC1 compact result claimed an audit outcome")
    if "complete_admission_passed" in payload:
        _fail("AC1 compact result claimed a classifier outcome")
    return result


@dataclass(frozen=True, slots=True)
class AC1Authority:
    authority_commit: str
    execution_authorized: bool = True
    tableau_selector: str = "SSPRK3"
    actual_spatial_operator: str = "inherited_RK4_2049_SBP4"
    retry: int = RETRY
    channel_count: int = 18
    state_advance_authorized: bool = False
    production_comparator_authorized: bool = False


def authorize(
    root: Path,
    config_raw: bytes,
    result_raw: bytes,
    authority_commit: str,
) -> AC1Authority:
    repository = root.resolve()
    if not _COMMIT.fullmatch(authority_commit):
        _fail("AC1 authority commit is malformed")
    validate_compact(config_raw, result_raw)
    ti2.ti1._environment()  # type: ignore[attr-defined]
    head = _git(repository, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    parents = (
        _git(repository, "rev-list", "--parents", "-n", "1", authority_commit)
        .decode()
        .split()[1:]
    )
    if head != authority_commit or parents != [BASE_COMMIT]:
        _fail("AC1 authority is not the exact direct successor")
    _require_no_replace_refs(repository)
    unstaged, staged, untracked = _working(repository)
    if unstaged or staged or any(item not in ALLOWED_UNTRACKED for item in untracked):
        _fail("AC1 authority working image is not clean")
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
        _fail("committed AC1 delta differs")
    if _git(repository, "show", f"{authority_commit}:{CONFIG_PATH}") != config_raw:
        _fail("committed AC1 config differs")
    if _git(repository, "show", f"{authority_commit}:{RESULT_PATH}") != result_raw:
        _fail("committed AC1 compact result differs")
    _verify_predecessor_bindings(repository)
    _verify_live_bindings(repository, commit=authority_commit)
    if _snapshot_store(repository) != (
        SEALED_STORE_LEAF_COUNT,
        SEALED_STORE_SNAPSHOT_SHA256,
    ):
        _fail("sealed store differs before AC1 authorization")
    if list(real_predecessor_fingerprints(repository)) != [
        _expected_fingerprint_receipt()
    ]:
        _fail("restored predecessor fingerprints differ before authorization")
    return AC1Authority(authority_commit=authority_commit)


__all__ = [
    "ARTIFACT_ID",
    "CHANNEL_ORDER",
    "CONFIG_PATH",
    "DELTA_PATHS",
    "HISTORICAL_JOURNAL_SHA256",
    "OUTPUT_NAMESPACE",
    "RESULT_PATH",
    "AC1Authority",
    "AC1AuthorityError",
    "authorize",
    "build_prelaunch",
    "canonical_pretty",
    "expected_compact",
    "expected_raw_replay_receipt",
    "parse_config",
    "read_leaf",
    "real_predecessor_fingerprints",
    "require_output_absent",
    "validate_compact",
]
