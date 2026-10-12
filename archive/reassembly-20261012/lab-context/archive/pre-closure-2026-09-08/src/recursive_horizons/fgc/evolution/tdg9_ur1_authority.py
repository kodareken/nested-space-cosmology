"""Prospective authority for one retry-3 u:R envelope-owner diagnostic."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import re
import stat
from typing import NoReturn

from . import tdg9_ac1_authority as ac1
from . import tdg9_ac1_pref1_binder as pref1


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG9-UR1-FRZ1"
CLASSIFICATION = (
    "premise_only_retry3_u_R_row_hash_and_zero_lower_envelope_owner_diagnostic"
)
PROJECT_VERSION = ac1.PROJECT_VERSION
TARGET_PROTOCOL = ac1.TARGET_PROTOCOL
CONFIG_PATH = "configs/fgc/fgc-1-tdg9-ur1-frz1.toml"
RESULT_PATH = "results/fgc-1-tdg9-ur1-frz1.json"
OWNER_DOCUMENT = "docs/fgc-tdg9-ur1-frz1.md"
OUTPUT_NAMESPACE = "runs/fgc-2-sf1/tdg9-ur1/retry3-ur-envelope-owner"
STAGING_PREFIX = ".retry3-ur-envelope-owner-ur1.stage-"
BASE_COMMIT = "22f7dc5e181af45eb43f6f19dee23fbc29f89c0c"

AC1_AUTHORITY_COMMIT = pref1.AUTHORITY_COMMIT
AC1_ARTIFACT_ID = pref1.AC1_ARTIFACT_ID
AC1_PREF1_ARTIFACT_ID = pref1.ARTIFACT_ID
AC1_PREF1_RESULT_PATH = pref1.RESULT_PATH
AC1_PREF1_RESULT_SHA256 = (
    "173183871b2d750f98e9758cdc5f702be40b0bbb73d10e707672e0e68752935d"
)
AC1_RAW_NAMESPACE = pref1.RAW_NAMESPACE
AC1_RAW_SCHEMA = pref1.RAW_SCHEMA
AC1_RAW_CLASSIFICATION = pref1.RAW_CLASSIFICATION
AC1_RAW_MANIFEST_SHA256 = (
    "5cd74aa33555c6253353311c246515680b50c8c908a4ff295949a98a7e8f0912"
)
AC1_RAW_TERMINAL_SHA256 = (
    "515ddf31d2fe919e4228668bea8701f282931d45e462542c26675382c95b834c"
)

RUNNER_PATH = "scripts/run_fgc_tdg9_ur1.py"
RUNNER_SHA256 = (
    "f660ba40df882d6ded87d8b0c020f8328fa50fa26a0a61f7f5686febe32c5680"
)
IMPLEMENTATION_BINDINGS = (
    (
        "TI2_runner_path",
        "scripts/run_fgc_tdg9_ti2.py",
        "488f05b65b71dd2617d1d77c2efa49496aa1a62698b37866e07bbb9c5792cd57",
    ),
    (
        "LOC1_runner_path",
        "scripts/run_fgc_tdg9_loc1.py",
        "19d12813b7ee3fde6b044b224124d62909170fdfaffb728414a378984bac0996",
    ),
    (
        "envelope_runtime_path",
        "src/recursive_horizons/fgc/evolution/tdg5_stage_complete_refinement_runtime.py",
        "75fc11943f459ca7549f75a812fe961e55ae04bbdb13450c1fdd509b5b53a080",
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
        "tableau_runtime_path",
        "src/recursive_horizons/fgc/evolution/numerical_engine.py",
        "8ea99604a85b4e1acbf869ed2462ea4aa8301c06e068330a1d0ce66a51a0bebf",
    ),
    (
        "AC1_authority_path",
        "src/recursive_horizons/fgc/evolution/tdg9_ac1_authority.py",
        "3ccd293ff67413466af9390065357e68226ddd46ee23df7bceda10aec6dc4c19",
    ),
    (
        "AC1_PREF1_binder_path",
        "src/recursive_horizons/fgc/evolution/tdg9_ac1_pref1_binder.py",
        "d53f6c7501629f375b477999d79e5e894241c3899e767405e8ee032c8c60c17d",
    ),
)

SEALED_STORE_LEAF_COUNT = 115
SEALED_STORE_SNAPSHOT_SHA256 = (
    "5917d70de3080dcc6b7abeb7fdb7cc3370c6140cbeb37bde0c142d6a2d36a445"
)
MEMBER_KEY = pref1.MEMBER_KEY
MEMBER_DESCRIPTOR_SHA256 = pref1.MEMBER_DESCRIPTOR_SHA256
PHYSICAL_STATE_SHA256 = pref1.PHYSICAL_STATE_SHA256
ACCEPTED_TIME_HEX = pref1.ACCEPTED_TIME_HEX
POINT_COUNT = pref1.POINT_COUNT
OWNED_ROW_COUNT = pref1.OWNED_ROW_COUNT
RETRY = pref1.RETRY
WIDTH_HEX = pref1.WIDTH_HEX
TRANSACTION_SHA256 = pref1.TRANSACTION_SHA256
HISTORICAL_JOURNAL_SHA256 = pref1.HISTORICAL_JOURNAL_SHA256
REPLAY = dict(ac1.REPLAY)
ENVIRONMENT = dict(ac1.ENVIRONMENT)

PUBLISHED_CHANNEL = "u:R"
ROW_HASH_DOMAIN = "TDG9-AR1-BINARY64-HERMITE-ROWS-v1"
TI2_RETRY3_U_R_ROW_SHA256 = (
    "b4a48c7bf72e015f3f3d9a340ea550d608514fb8cd886562265158c33c0b3969"
)
TI2_RETRY3_U_R_RELATED_COMPLETE_C_CLASS = "sufficient_contraction_pass"
AC1_U_R_D01_LOWER_HEX = "0x0.0p+0"
AC1_U_R_D01_UPPER_HEX = "0x1.2d198e246e459p-38"
AC1_U_R_D12_LOWER_HEX = "0x0.0p+0"
AC1_U_R_D12_UPPER_HEX = "0x1.2d1c31bb91376p-38"
AC1_U_R_CLASS = "order_inconclusive"
ZERO_LOWER_OWNERS = (
    "raw_candidate_maximum_is_zero",
    "coefficient_plus_arithmetic_debit_clips_positive_raw_maximum",
    "downward_binary64_rounding_of_positive_exact_lower",
)

SEMANTICS = {
    "experiment_label": (
        "retry3_SSPRK3_on_inherited_SBP4_u_R_row_hash_and_"
        "zero_lower_envelope_owner_diagnostic"
    ),
    "tableau_selector": "SSPRK3",
    "tableau_runtime_selector": ac1.SEMANTICS["tableau_runtime_selector"],
    "source_member_method": "RK4",
    "actual_spatial_operator": "inherited_RK4_2049_SBP4",
    "published_channel": PUBLISHED_CHANNEL,
    "row_hash_domain": ROW_HASH_DOMAIN,
    "envelope_owner": "Binary64CubicEnvelope",
    "classifier_owner": "classify_tdg6_channel",
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
    "published_channel_count": 1,
    "published_channel": PUBLISHED_CHANNEL,
    "published_row_count": OWNED_ROW_COUNT,
    "fourth_width_authorized": False,
    "retry_4_authorized": False,
    "retry_5_authorized": False,
    "resource_escalation_authorized": False,
}
DECISION = {
    "exact_lower_formula": (
        "max(0,Fraction(raw_candidate_maximum)-"
        "Fraction(coefficient_construction_debit)-"
        "Fraction(outward_arithmetic_debit))"
    ),
    "bernstein_certification_slack_recorded": True,
    "bernstein_certification_slack_part_of_lower_clip": False,
    "stored_lower_reproduced_by_downward_binary64_rounding": True,
    "zero_lower_owner_enum": list(ZERO_LOWER_OWNERS),
    "AC1_D01_interval_hex": [AC1_U_R_D01_LOWER_HEX, AC1_U_R_D01_UPPER_HEX],
    "AC1_D12_interval_hex": [AC1_U_R_D12_LOWER_HEX, AC1_U_R_D12_UPPER_HEX],
    "required_design_classifier": AC1_U_R_CLASS,
    "related_TI2_complete_C_class_recorded_only_after_row_hash_match": True,
    "related_TI2_complete_C_class_is_not_replacement_admission": True,
    "successor_remedy_selected": False,
}
SCOPE = {
    "live_store_read_only": True,
    "AC1_raw_read_only": True,
    "compact_verifier_store_blind": True,
    "compact_verifier_raw_blind": True,
    "compact_verifier_shadow_blind": True,
    "prelaunch_may_restore_retry_3": True,
    "prelaunch_constructs_zero_shadow_proposals": True,
    "run_publishes_only_u_R": True,
    "campaign_store_mutation_authorized": False,
    "PDE_state_commit_authorized": False,
    "temporal_retry_admission_authorized": False,
    "fine_path_commit_authorized": False,
    "continuation_authorized": False,
    "envelope_or_admission_runtime_change_authorized": False,
    "GR0_calibration_authorized": False,
    "candidate_branches_authorized": False,
    "mechanism_result_authorized": False,
    "physical_result_authorized": False,
}
CLAIMS = {
    "AC1_PREF1_bound": True,
    "UR1_execution_authorized": True,
    "UR1_result_earned": False,
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
            RESULT_PATH,
            "results/README.md",
            "scripts/check_repo.py",
            "scripts/reproduce_fgc_tdg9_ur1_frz1.py",
            RUNNER_PATH,
            "src/recursive_horizons/fgc/evolution/tdg9_ur1_authority.py",
            "tests/test_check_repo_tdg9_ur1_frz1.py",
            "tests/test_fgc_tdg9_ur1_authority.py",
            "tests/test_fgc_tdg9_ur1_runner.py",
        )
    )
)
ALLOWED_UNTRACKED = ac1.ALLOWED_UNTRACKED
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")


class UR1AuthorityError(RuntimeError):
    """The prospective UR1 image differs from its frozen authority."""


def _fail(message: str) -> NoReturn:
    raise UR1AuthorityError(message)


canonical_pretty = ac1.canonical_pretty
read_leaf = ac1.read_leaf


def _git(root: Path, *arguments: str) -> bytes:
    try:
        return ac1._git(root, *arguments)  # type: ignore[attr-defined]
    except ac1.AC1AuthorityError as exc:
        raise UR1AuthorityError("Git image differs") from exc


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
    if _git(root, "for-each-ref", "--format=%(refname)", "refs/replace").strip():
        _fail("Git replace refs are present")


def _implementation() -> dict[str, object]:
    value: dict[str, object] = {
        "runner_path": RUNNER_PATH,
        "runner_sha256": RUNNER_SHA256,
    }
    for key, path, digest in IMPLEMENTATION_BINDINGS:
        value[key] = path
        value[key.replace("_path", "_sha256")] = digest
    return value


def expected_raw_replay_receipt() -> dict[str, object]:
    return {
        "accepted_time_hex": ACCEPTED_TIME_HEX,
        "attempted_width_hex": WIDTH_HEX,
        "descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
        "historical_journal_sha256": HISTORICAL_JOURNAL_SHA256,
        "member_key": MEMBER_KEY,
        "retry": RETRY,
        "state_sha256": PHYSICAL_STATE_SHA256,
        "transaction_sha256": TRANSACTION_SHA256,
    }


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
            "AC1_authority_commit": AC1_AUTHORITY_COMMIT,
            "AC1_artifact_id": AC1_ARTIFACT_ID,
            "AC1_PREF1_artifact_id": AC1_PREF1_ARTIFACT_ID,
            "AC1_PREF1_result_path": AC1_PREF1_RESULT_PATH,
            "AC1_PREF1_result_sha256": AC1_PREF1_RESULT_SHA256,
            "AC1_raw_namespace": AC1_RAW_NAMESPACE,
            "AC1_raw_schema": AC1_RAW_SCHEMA,
            "AC1_raw_classification": AC1_RAW_CLASSIFICATION,
            "AC1_raw_manifest_sha256": AC1_RAW_MANIFEST_SHA256,
            "AC1_raw_terminal_sha256": AC1_RAW_TERMINAL_SHA256,
            "sealed_store_leaf_count": SEALED_STORE_LEAF_COUNT,
            "sealed_store_snapshot_sha256": SEALED_STORE_SNAPSHOT_SHA256,
        },
        "implementation": _implementation(),
        "fingerprint_encoding": dict(ac1.ti2.FINGERPRINT_ENCODING),
        "semantics": dict(SEMANTICS),
        "selection": {
            "member_key": MEMBER_KEY,
            "descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
            "physical_state_sha256": PHYSICAL_STATE_SHA256,
            "accepted_time_hex": ACCEPTED_TIME_HEX,
            "point_count": POINT_COUNT,
            "owned_row_count": OWNED_ROW_COUNT,
            "retry": RETRY,
            "attempted_width_hex": WIDTH_HEX,
            "published_channel": PUBLISHED_CHANNEL,
            "row_hash_domain": ROW_HASH_DOMAIN,
            "TI2_retry3_u_R_row_sha256": TI2_RETRY3_U_R_ROW_SHA256,
            "TI2_retry3_u_R_related_complete_C_class": (
                TI2_RETRY3_U_R_RELATED_COMPLETE_C_CLASS
            ),
        },
        "replay": dict(REPLAY),
        "replay_receipt": expected_raw_replay_receipt(),
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
        raise UR1AuthorityError("invalid UR1 TOML") from exc
    if value != _expected_config():
        _fail("UR1 config differs from frozen contract")
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
            _fail("UR1 output parent is unsafe")
    try:
        target.lstat()
    except FileNotFoundError:
        pass
    else:
        _fail("UR1 output namespace already exists")
    if any(item.name.startswith(STAGING_PREFIX) for item in target.parent.iterdir()):
        _fail("UR1 staging namespace already exists")


def _verify_predecessor_bindings(root: Path) -> None:
    raw = _git(root, "show", f"{BASE_COMMIT}:{AC1_PREF1_RESULT_PATH}")
    if sha256(raw).hexdigest() != AC1_PREF1_RESULT_SHA256:
        _fail("historical AC1-PREF1 compact binding differs")


def _verify_live_bindings(root: Path, *, commit: str | None) -> None:
    bindings = (
        (RUNNER_PATH, RUNNER_SHA256),
        *((path, digest) for _key, path, digest in IMPLEMENTATION_BINDINGS),
    )
    for path, digest in bindings:
        raw = read_leaf(root, path) if commit is None else _git(root, "show", f"{commit}:{path}")
        if sha256(raw).hexdigest() != digest:
            _fail(f"UR1 implementation binding differs: {path}")


def _verify_ac1_raw(root: Path) -> dict[str, object]:
    try:
        manifest_raw, terminal_raw = pref1._raw_pair(root)  # type: ignore[attr-defined]
    except Exception as exc:
        raise UR1AuthorityError("sealed AC1 raw pair cannot be authenticated") from exc
    observed = {
        "namespace": AC1_RAW_NAMESPACE,
        "leaf_count": 2,
        "manifest_sha256": sha256(manifest_raw).hexdigest(),
        "terminal_sha256": sha256(terminal_raw).hexdigest(),
    }
    expected = {
        "namespace": AC1_RAW_NAMESPACE,
        "leaf_count": 2,
        "manifest_sha256": AC1_RAW_MANIFEST_SHA256,
        "terminal_sha256": AC1_RAW_TERMINAL_SHA256,
    }
    if observed != expected:
        _fail("sealed AC1 raw pair differs")
    return observed


def _snapshot_store(root: Path) -> tuple[int, str]:
    try:
        return ac1._snapshot_store(root)  # type: ignore[attr-defined]
    except ac1.AC1AuthorityError as exc:
        raise UR1AuthorityError("sealed store cannot be authenticated") from exc


def _expected_fingerprint_receipt() -> dict[str, object]:
    return {
        "retry": RETRY,
        "transaction_sha256": TRANSACTION_SHA256,
        "state_sha256": PHYSICAL_STATE_SHA256,
        "accepted_time_hex": ACCEPTED_TIME_HEX,
        "descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
        "shadow_proposal_constructed": False,
    }


def real_predecessor_fingerprints(root: Path) -> tuple[dict[str, object], ...]:
    """Restore only retry 3 and fingerprint it without constructing a proposal."""

    try:
        receipts = ac1.real_predecessor_fingerprints(root.resolve())
    except ac1.AC1AuthorityError as exc:
        raise UR1AuthorityError("retry-3 predecessor fingerprint failed") from exc
    expected = (_expected_fingerprint_receipt(),)
    if receipts != expected:
        _fail("retry-3 predecessor fingerprint differs")
    return receipts


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
            "AC1_raw_identity_observed_at_prelaunch": {
                "namespace": AC1_RAW_NAMESPACE,
                "leaf_count": 2,
                "manifest_sha256": AC1_RAW_MANIFEST_SHA256,
                "terminal_sha256": AC1_RAW_TERMINAL_SHA256,
            },
            "UR1_output_namespace_absent_at_prelaunch": True,
            "replace_refs_absent": True,
            "restored_predecessor_fingerprints": [_expected_fingerprint_receipt()],
            "shadow_executed": False,
            "shadow_proposals_constructed": 0,
            "UR1_result_earned": False,
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
    ac1.ti2.ti1._environment()  # type: ignore[attr-defined]
    head = _git(repository, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    if head != BASE_COMMIT:
        _fail("prelaunch HEAD is not the sealed AC1-PREF1 commit")
    _require_no_replace_refs(repository)
    unstaged, staged, untracked = _working(repository)
    observed = tuple(
        sorted((*unstaged, *(item for item in untracked if item not in ALLOWED_UNTRACKED)))
    )
    if staged or observed != DELTA_PATHS:
        _fail("prospective UR1 delta differs")
    for relative in DELTA_PATHS:
        read_leaf(repository, relative)
    _verify_predecessor_bindings(repository)
    _verify_live_bindings(repository, commit=None)
    sealed_store = (SEALED_STORE_LEAF_COUNT, SEALED_STORE_SNAPSHOT_SHA256)
    if store_snapshot != sealed_store:
        _fail("sealed store differs at UR1 prelaunch")
    ac1_raw = _verify_ac1_raw(repository)
    require_output_absent(repository)
    fingerprints = list(real_predecessor_fingerprints(repository))
    if fingerprints != [_expected_fingerprint_receipt()]:
        _fail("restored predecessor fingerprints differ")
    value = expected_compact(config_raw)
    return {
        **value,
        "artifact_payload": {
            **value["artifact_payload"],
            "store_identity_observed_at_prelaunch": {
                "leaf_count": store_snapshot[0],
                "sha256": store_snapshot[1],
            },
            "AC1_raw_identity_observed_at_prelaunch": ac1_raw,
            "restored_predecessor_fingerprints": fingerprints,
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
    _fail(f"nonfinite UR1 compact JSON constant: {token}")


def validate_compact(config_raw: bytes, result_raw: bytes) -> dict[str, object]:
    expected = expected_compact(config_raw)
    try:
        result = json.loads(
            result_raw,
            object_pairs_hook=_unique,
            parse_constant=_reject_compact_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise UR1AuthorityError("invalid UR1 compact JSON") from exc
    if result != expected or result_raw != canonical_pretty(result):
        _fail("UR1 compact result differs")
    payload = result["artifact_payload"]
    if payload["claims"]["UR1_result_earned"] is not False:
        _fail("UR1 compact result claimed a diagnostic outcome")
    if payload["shadow_executed"] is not False or payload["shadow_proposals_constructed"] != 0:
        _fail("UR1 compact result claimed shadow work")
    if "zero_lower_owner" in payload or "u_R" in payload:
        _fail("UR1 compact result claimed a diagnostic outcome")
    return result


@dataclass(frozen=True, slots=True)
class UR1Authority:
    authority_commit: str
    execution_authorized: bool = True
    tableau_selector: str = "SSPRK3"
    actual_spatial_operator: str = "inherited_RK4_2049_SBP4"
    retry: int = RETRY
    published_channel: str = PUBLISHED_CHANNEL
    state_advance_authorized: bool = False
    temporal_admission_authorized: bool = False


def authorize(
    root: Path,
    config_raw: bytes,
    result_raw: bytes,
    authority_commit: str,
) -> UR1Authority:
    repository = root.resolve()
    if not _COMMIT.fullmatch(authority_commit):
        _fail("UR1 authority commit is malformed")
    validate_compact(config_raw, result_raw)
    ac1.ti2.ti1._environment()  # type: ignore[attr-defined]
    head = _git(repository, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    parents = (
        _git(repository, "rev-list", "--parents", "-n", "1", authority_commit)
        .decode()
        .split()[1:]
    )
    if head != authority_commit or parents != [BASE_COMMIT]:
        _fail("UR1 authority is not the exact direct successor")
    _require_no_replace_refs(repository)
    unstaged, staged, untracked = _working(repository)
    if unstaged or staged or any(item not in ALLOWED_UNTRACKED for item in untracked):
        _fail("UR1 authority working image is not clean")
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
        _fail("committed UR1 delta differs")
    if _git(repository, "show", f"{authority_commit}:{CONFIG_PATH}") != config_raw:
        _fail("committed UR1 config differs")
    if _git(repository, "show", f"{authority_commit}:{RESULT_PATH}") != result_raw:
        _fail("committed UR1 compact result differs")
    _verify_predecessor_bindings(repository)
    _verify_live_bindings(repository, commit=authority_commit)
    if _snapshot_store(repository) != (
        SEALED_STORE_LEAF_COUNT,
        SEALED_STORE_SNAPSHOT_SHA256,
    ):
        _fail("sealed store differs before UR1 authorization")
    _verify_ac1_raw(repository)
    if real_predecessor_fingerprints(repository) != (_expected_fingerprint_receipt(),):
        _fail("restored predecessor fingerprints differ before authorization")
    return UR1Authority(authority_commit=authority_commit)


__all__ = [
    "AC1_U_R_CLASS",
    "AC1_U_R_D01_LOWER_HEX",
    "AC1_U_R_D01_UPPER_HEX",
    "AC1_U_R_D12_LOWER_HEX",
    "AC1_U_R_D12_UPPER_HEX",
    "ARTIFACT_ID",
    "CONFIG_PATH",
    "DELTA_PATHS",
    "OUTPUT_NAMESPACE",
    "PUBLISHED_CHANNEL",
    "RESULT_PATH",
    "ROW_HASH_DOMAIN",
    "STAGING_PREFIX",
    "TI2_RETRY3_U_R_RELATED_COMPLETE_C_CLASS",
    "TI2_RETRY3_U_R_ROW_SHA256",
    "UR1Authority",
    "UR1AuthorityError",
    "ZERO_LOWER_OWNERS",
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
