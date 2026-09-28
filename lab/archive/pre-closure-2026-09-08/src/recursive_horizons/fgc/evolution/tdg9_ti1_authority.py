"""Prospective authority for the bounded TDG9 TI1 tableau counterfactual."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import os
import platform
from pathlib import Path
import re
import stat
import subprocess
import sys
import tomllib
from typing import NoReturn

import numpy as np


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG9-TI1-FRZ1"
CLASSIFICATION = "premise_only_same_state_ssprk3_tableau_counterfactual_authority"
PROJECT_VERSION = "0.11.0"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CONFIG_PATH = "configs/fgc/fgc-1-tdg9-ti1-frz1.toml"
RESULT_PATH = "results/fgc-1-tdg9-ti1-frz1.json"
OWNER_DOCUMENT = "docs/fgc-tdg9-ti1-frz1.md"
OUTPUT_NAMESPACE = (
    "runs/fgc-2-sf1/tdg9-ti1/ssprk3-tableau-on-rk4-2049-retries-3-5"
)
STAGING_PREFIX = ".ssprk3-tableau-on-rk4-2049-retries-3-5.stage-"
BASE_COMMIT = "de7fb1dd52f85d6f618bd91e5e5b20953baa463f"

LOC2_PREF2_CONFIG_PATH = "configs/fgc/fgc-1-tdg9-loc2-pref2.toml"
LOC2_PREF2_CONFIG_SHA256 = (
    "eb67e45c79995ab642aa7945fa4d75c5a5987d6008cd0ffaa5ce8b7905b142a7"
)
LOC2_PREF2_RESULT_PATH = "results/fgc-1-tdg9-loc2-pref2.json"
LOC2_PREF2_RESULT_SHA256 = (
    "2f1630f449aa4d942cd0b0e45cf6ebb0f43862f2e387da51172584680aa03232"
)
LOC2_PREF2_BINDER_PATH = (
    "src/recursive_horizons/fgc/evolution/tdg9_loc2_pref2_binder.py"
)
LOC2_PREF2_BINDER_SHA256 = (
    "f509a35245b9bc9d6aca38b0eb51f1bb358a4e0945e78c4ab8cf3234f5e1b74d"
)
SEALED_STORE_LEAF_COUNT = 115
SEALED_STORE_SNAPSHOT_SHA256 = (
    "5917d70de3080dcc6b7abeb7fdb7cc3370c6140cbeb37bde0c142d6a2d36a445"
)

RUNNER_PATH = "scripts/run_fgc_tdg9_ti1.py"
RUNNER_SHA256 = "abee4a9c4aa25dc8c634eaf5f6f2c54e45d74a754740fea874042eef9db8d7bb"
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
        "AR1_authority_path",
        "src/recursive_horizons/fgc/evolution/tdg9_ar1_authority.py",
        "154c7dbe675bbf788cd199c6224f74e439b469cad0f6bf6b0f015a785001647d",
    ),
    (
        "AR1_reconstruction_path",
        "scripts/run_fgc_tdg9_loc1.py",
        "19d12813b7ee3fde6b044b224124d62909170fdfaffb728414a378984bac0996",
    ),
    (
        "LOC2_route_path",
        "scripts/run_fgc_tdg9_loc2.py",
        "e3ad84f818fe33afce707bcbec1e21a99e0892b3d924536549e0ead740351122",
    ),
    (
        "primary_localizer_path",
        "src/recursive_horizons/fgc/evolution/tdg9_local_extrema.py",
        "b97bca8c49313ad28322a1236143e03c5a7f0d3b6bcdb2940965e531a122b28a",
    ),
    (
        "independent_v2_localizer_path",
        "src/recursive_horizons/fgc/evolution/tdg9_local_extrema_independent_v2.py",
        "750c38bb6bad4b37d946f66708055f3dc0835c807244e9c25c6e87e4e87e0511",
    ),
)

MEMBER_KEY = "RK4-2049"
MEMBER_DESCRIPTOR_SHA256 = (
    "77847126340e78c8ac300fac2795bceda724c15e0894e34444b0da42a84166b4"
)
PHYSICAL_STATE_SHA256 = (
    "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a"
)
ACCEPTED_TIME_HEX = "0x1.78554de5a30e0p+0"
POINT_COUNT = 2049
OWNED_ROW_COUNT = 2044
COORDINATES_SHA256 = (
    "2a347b314de5e1e94ac0c999db283e03498070ece6762c53252fa41c1e04a646"
)
GRID_SPACING_HEX = "0x1.0000000000000p-4"
OUTER_RADIUS_HEX = "0x1.0000000000000p+7"
FAILED_OCCURRENCES = (
    (3, "u:alpha"),
    (3, "u:R"),
    (4, "u:alpha"),
    (4, "u:lambda"),
    (4, "u:R"),
    (5, "u:alpha"),
    (5, "u:v"),
    (5, "u:lambda"),
    (5, "u:R"),
    (5, "q:R"),
)
REPLAYS = (
    {
        "retry": 3,
        "prior_retry_count": 2,
        "predecessor_generation": 9,
        "checkpoint_sha256": "eb6fddc480c94aa7ed15c80399fc2b26637693efef02399f267075fc6c258e56",
        "checkpoint_raw_sha256": "07572d4450829ff2e2bad65e217a917b0d31ba677a3b8fc20a3be43e849ba846",
        "journal_sequence": 11,
        "journal_sha256": "341cd8cd328434d85774bcb452889ac1886c520a72e945a92337b72cf564161a",
        "journal_raw_sha256": "e4f9ca42f2d014d4bdbd875a0f24a9bbb9c04159fb0d9e0e7de8d0d72db61fdf",
        "attempted_width_hex": "0x1.aaa9612df8000p-11",
    },
    {
        "retry": 4,
        "prior_retry_count": 3,
        "predecessor_generation": 10,
        "checkpoint_sha256": "6dc263e7719c9a422e9fed1b81ac573f3127607a2285f3c55b095d9019f60187",
        "checkpoint_raw_sha256": "9abb59999809d22354d6b5f810bc592cd2673e5e4d354d0ec58e9ad53522dc4b",
        "journal_sequence": 13,
        "journal_sha256": "1bd21553f3b6e0af613997da094a46d4dcd966a0e4346692bae62b05046862ca",
        "journal_raw_sha256": "26542403a17b05664d5a447b88860186d5f554b1f4d555624493ed40b941694c",
        "attempted_width_hex": "0x1.aaa9612df8000p-12",
    },
    {
        "retry": 5,
        "prior_retry_count": 4,
        "predecessor_generation": 11,
        "checkpoint_sha256": "11ec8a80d300de72075661a97b616aecead8ec78208bb541b93e62ca4c8124c8",
        "checkpoint_raw_sha256": "786bdcb835e043c65ddd33fad3a2af0e3d33e768dc26a9f9f00d2757b0da3ea5",
        "journal_sequence": 15,
        "journal_sha256": "475014d918952adef032d71b81073773364a37fa795491e739aae609f77e8537",
        "journal_raw_sha256": "1997342415f74f330ce6b45ca6fd550ad0e5e674b85148fea5047ff82ac3f152",
        "attempted_width_hex": "0x1.aaa9612df0000p-13",
    },
)

BASE_CUBIC_COUNT = 122_640
COMPONENT_CUBIC_COUNT = 367_920
FAMILY_CANDIDATE_CEILING = 490_560
AGGREGATE_CANDIDATE_CEILING = 1_471_680
PRIMARY_REFINEMENT_DEPTH = 160
GLOBAL_PROOF_BIT_CEILING = 9216

SEMANTICS = {
    "experiment_label": (
        "same_state_SSPRK3_tableau_counterfactual_on_RK4_2049_SBP4_operator"
    ),
    "tableau_selector": "SSPRK3",
    "tableau_runtime_selector": (
        "second_order_diagonal_norm_SBP_plus_SSPRK3"
    ),
    "source_member_method": "RK4",
    "source_member_integrator": (
        "fourth_order_diagonal_norm_SBP_with_second_order_boundary_closure_plus_RK4"
    ),
    "actual_spatial_operator": "inherited_RK4_2049_SBP4",
    "RHS_owner": "Proto12GR0EvolutionOperator",
    "projector_owner": "make_gr0_center_boundary_projector.<locals>.projector",
    "transaction_owner": "GR0RuntimeStageTransaction",
    "tracer_owner": "NormalFlowTracers",
    "ledger_owner": "TDG6TemporalLedger",
    "only_changed_variable": "explicit_time_tableau_selector",
    "production_SSPRK3_comparator": False,
    "independent_method_agreement": False,
}
WORK_BUDGET = {
    "retry_count": 3,
    "shadow_paths_per_retry": 7,
    "maximum_shadow_proposals": 21,
    "SSPRK3_records_per_proposal": 4,
    "maximum_stage_and_endpoint_RHS_records": 84,
    "occurrence_count": 10,
    "level_count": 20,
    "component_level_count": 60,
    "base_cubic_count": BASE_CUBIC_COUNT,
    "component_cubic_count": COMPONENT_CUBIC_COUNT,
    "family_candidate_ceiling": FAMILY_CANDIDATE_CEILING,
    "aggregate_candidate_ceiling": AGGREGATE_CANDIDATE_CEILING,
    "primary_refinement_depth": PRIMARY_REFINEMENT_DEPTH,
    "global_proof_bit_ceiling": GLOBAL_PROOF_BIT_CEILING,
    "fourth_width_authorized": False,
    "resource_escalation_authorized": False,
}
DECISION = {
    "unchanged_contraction_predicate": "D12_squared_times_8_le_D01_squared",
    "minimum_observed_order": "3/2",
    "original_complete_failure_required": True,
    "changed_ownership_alone_counts_as_clearance": False,
    "clearance_requires_shadow_complete_pass_or_exact_zero": True,
    "inconclusive_preempts_clearance_or_persistence_terminal": True,
    "successor_remedy_selected": False,
}
ENVIRONMENT = {
    "python_implementation": "CPython",
    "python_version": "3.14.3",
    "numpy_version": "2.5.1",
    "system": "Darwin",
    "sys_platform": "darwin",
    "machine": "arm64",
    "byteorder": "little",
}
SCOPE = {
    "live_store_read_only": True,
    "compact_verifier_store_blind": True,
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
    "LOC2_PREF2_bound": True,
    "TI1_execution_authorized": True,
    "TI1_result_earned": False,
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
            "scripts/reproduce_fgc_tdg9_ti1_frz1.py",
            RUNNER_PATH,
            "src/recursive_horizons/fgc/evolution/tdg9_ti1_authority.py",
            "tests/test_check_repo_tdg9_ti1_frz1.py",
            "tests/test_fgc_tdg9_ti1_authority.py",
            "tests/test_fgc_tdg9_ti1_runner.py",
        )
    )
)
ALLOWED_UNTRACKED = (".qdrant-initialized",)
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")


class TI1AuthorityError(RuntimeError):
    """The prospective TI1 image differs from its frozen authority."""


def _fail(message: str) -> NoReturn:
    raise TI1AuthorityError(message)


def canonical_pretty(value: object) -> bytes:
    return (
        json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"
    ).encode()


def _unique(items: list[tuple[str, object]]) -> dict[str, object]:
    answer: dict[str, object] = {}
    for key, value in items:
        if key in answer:
            _fail("duplicate compact JSON key")
        answer[key] = value
    return answer


def _git(root: Path, *arguments: str) -> bytes:
    environment = dict(
        os.environ,
        LC_ALL="C",
        LANG="C",
        GIT_NO_REPLACE_OBJECTS="1",
        GIT_CONFIG_NOSYSTEM="1",
    )
    try:
        return subprocess.run(
            (
                "git",
                "--no-replace-objects",
                "--no-optional-locks",
                "-c",
                "core.fsmonitor=false",
                "-c",
                "core.untrackedCache=false",
                *arguments,
            ),
            cwd=root,
            env=environment,
            check=True,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
    except subprocess.CalledProcessError as exc:
        raise TI1AuthorityError("Git image differs") from exc


def _paths(root: Path, *arguments: str) -> tuple[str, ...]:
    raw = _git(root, *arguments)
    if not raw:
        return ()
    chunks = raw.split(b"\0")
    if chunks[-1] != b"":
        _fail("Git path stream differs")
    return tuple(item.decode() for item in chunks[:-1])


def read_leaf(root: Path, relative: str) -> bytes:
    """Read one regular leaf without following a link or accepting a race."""

    candidate = Path(relative)
    if candidate.is_absolute() or any(
        part in {"", ".", ".."} for part in candidate.parts
    ):
        _fail("unsafe authority path")
    path = root / candidate
    before = path.lstat()
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        _fail("unsafe authority leaf")
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        active = os.fstat(descriptor)
        blocks: list[bytes] = []
        while block := os.read(descriptor, 1 << 20):
            blocks.append(block)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    identity = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    if (
        identity
        != (active.st_dev, active.st_ino, active.st_size, active.st_mtime_ns)
        or identity
        != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    ):
        _fail("authority leaf changed")
    raw = b"".join(blocks)
    if len(raw) != before.st_size:
        _fail("authority leaf size changed")
    return raw


def _working(root: Path) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    return (
        _paths(root, "diff", "--name-only", "-z", "--"),
        _paths(root, "diff", "--cached", "--name-only", "-z", "--"),
        _paths(root, "ls-files", "--others", "--exclude-standard", "-z", "--"),
    )


def _environment() -> dict[str, str]:
    observed = {
        "python_implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "system": platform.system(),
        "sys_platform": sys.platform,
        "machine": platform.machine(),
        "byteorder": sys.byteorder,
    }
    if observed != ENVIRONMENT:
        _fail("execution environment differs")
    return observed


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
            "LOC2_PREF2_config_path": LOC2_PREF2_CONFIG_PATH,
            "LOC2_PREF2_config_sha256": LOC2_PREF2_CONFIG_SHA256,
            "LOC2_PREF2_result_path": LOC2_PREF2_RESULT_PATH,
            "LOC2_PREF2_result_sha256": LOC2_PREF2_RESULT_SHA256,
            "LOC2_PREF2_binder_path": LOC2_PREF2_BINDER_PATH,
            "LOC2_PREF2_binder_sha256": LOC2_PREF2_BINDER_SHA256,
            "sealed_store_leaf_count": SEALED_STORE_LEAF_COUNT,
            "sealed_store_snapshot_sha256": SEALED_STORE_SNAPSHOT_SHA256,
        },
        "implementation": _implementation(),
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
    try:
        value = tomllib.loads(raw.decode())
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise TI1AuthorityError("invalid TI1 TOML") from exc
    if value != _expected_config():
        _fail("TI1 config differs from frozen contract")
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
            _fail("TI1 output parent is unsafe")
    try:
        target.lstat()
    except FileNotFoundError:
        pass
    else:
        _fail("TI1 output namespace already exists")
    if any(item.name.startswith(STAGING_PREFIX) for item in os.scandir(target.parent)):
        _fail("TI1 staging namespace already exists")


def _verify_hash_bindings(root: Path, *, commit: str | None) -> None:
    bindings = (
        (LOC2_PREF2_CONFIG_PATH, LOC2_PREF2_CONFIG_SHA256),
        (LOC2_PREF2_RESULT_PATH, LOC2_PREF2_RESULT_SHA256),
        (LOC2_PREF2_BINDER_PATH, LOC2_PREF2_BINDER_SHA256),
        (RUNNER_PATH, RUNNER_SHA256),
        *((path, digest) for _key, path, digest in IMPLEMENTATION_BINDINGS),
    )
    for path, digest in bindings:
        raw = read_leaf(root, path) if commit is None else _git(root, "show", f"{commit}:{path}")
        if sha256(raw).hexdigest() != digest:
            _fail(f"TI1 implementation binding differs: {path}")


def build_prelaunch(
    config_raw: bytes,
    root: Path,
    *,
    store_snapshot: tuple[int, str],
) -> dict[str, object]:
    repository = root.resolve()
    parse_config(config_raw)
    _environment()
    head = _git(repository, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    if head != BASE_COMMIT:
        _fail("prelaunch HEAD is not sealed LOC2 PREF2")
    unstaged, staged, untracked = _working(repository)
    observed = tuple(
        sorted((*unstaged, *(item for item in untracked if item not in ALLOWED_UNTRACKED)))
    )
    if staged or observed != DELTA_PATHS:
        _fail("prospective TI1 delta differs")
    for relative in DELTA_PATHS:
        read_leaf(repository, relative)
    _verify_hash_bindings(repository, commit=None)
    if store_snapshot != (SEALED_STORE_LEAF_COUNT, SEALED_STORE_SNAPSHOT_SHA256):
        _fail("sealed campaign store differs at TI1 prelaunch")
    require_output_absent(repository)
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
            "output_namespace_absent_at_prelaunch": True,
            "shadow_executed": False,
        },
    }


def validate_compact(config_raw: bytes, result_raw: bytes) -> dict[str, object]:
    parse_config(config_raw)
    try:
        result = json.loads(result_raw, object_pairs_hook=_unique)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TI1AuthorityError("invalid TI1 compact JSON") from exc
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
            "output_namespace_absent_at_prelaunch": True,
            "shadow_executed": False,
        },
    }
    if result != expected or result_raw != canonical_pretty(result):
        _fail("TI1 compact result differs")
    return result


@dataclass(frozen=True, slots=True)
class TI1Authority:
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
) -> TI1Authority:
    repository = root.resolve()
    if not _COMMIT.fullmatch(authority_commit):
        _fail("TI1 authority commit is malformed")
    validate_compact(config_raw, result_raw)
    _environment()
    head = _git(repository, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    parents = _git(
        repository, "rev-list", "--parents", "-n", "1", authority_commit
    ).decode().split()[1:]
    if head != authority_commit or parents != [BASE_COMMIT]:
        _fail("TI1 authority is not the exact direct successor")
    unstaged, staged, untracked = _working(repository)
    if unstaged or staged or any(item not in ALLOWED_UNTRACKED for item in untracked):
        _fail("TI1 authority working image is not clean")
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
        _fail("committed TI1 delta differs")
    if _git(repository, "show", f"{authority_commit}:{CONFIG_PATH}") != config_raw:
        _fail("committed TI1 config differs")
    if _git(repository, "show", f"{authority_commit}:{RESULT_PATH}") != result_raw:
        _fail("committed TI1 result differs")
    _verify_hash_bindings(repository, commit=authority_commit)
    return TI1Authority(authority_commit)


__all__ = (
    "ACCEPTED_TIME_HEX",
    "AGGREGATE_CANDIDATE_CEILING",
    "ARTIFACT_ID",
    "BASE_COMMIT",
    "BASE_CUBIC_COUNT",
    "CLASSIFICATION",
    "CONFIG_PATH",
    "COORDINATES_SHA256",
    "DELTA_PATHS",
    "FAILED_OCCURRENCES",
    "FAMILY_CANDIDATE_CEILING",
    "LOC2_PREF2_CONFIG_PATH",
    "LOC2_PREF2_RESULT_PATH",
    "MEMBER_KEY",
    "MEMBER_DESCRIPTOR_SHA256",
    "OUTPUT_NAMESPACE",
    "OWNED_ROW_COUNT",
    "PHYSICAL_STATE_SHA256",
    "PRIMARY_REFINEMENT_DEPTH",
    "REPLAYS",
    "RESULT_PATH",
    "RUNNER_PATH",
    "SEMANTICS",
    "SEALED_STORE_LEAF_COUNT",
    "SEALED_STORE_SNAPSHOT_SHA256",
    "STAGING_PREFIX",
    "TI1Authority",
    "TI1AuthorityError",
    "WORK_BUDGET",
    "authorize",
    "build_prelaunch",
    "canonical_pretty",
    "parse_config",
    "read_leaf",
    "require_output_absent",
    "validate_compact",
)
