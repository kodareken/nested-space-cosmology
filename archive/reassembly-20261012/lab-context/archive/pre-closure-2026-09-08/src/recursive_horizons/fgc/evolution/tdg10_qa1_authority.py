"""Prospective authority for one retry-3 exact complete-C qualification.

TDG10-QA1 may authorize one diagnostic, non-accepted measurement of the
immutable retry-3 GR-0 proposal under SSPRK3-on-inherited-SBP4 with the
exact radius-free complete-C interval owner.  It does not authorize a
campaign-state advance, a new protocol ID, or restoration of generation 10.
"""

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
from typing import Mapping, NoReturn

from . import tdg9_ti2_authority as ti2
from .tdg6_temporal_admission_design import TDG6_COMPLETE_STATE_CHANNELS


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG10-QA1-FRZ1"
CLASSIFICATION = (
    "premise_only_retry3_ssprk3_sbp4_exact_complete_C_all_channel_qualification"
)
PROJECT_VERSION = "0.11.0"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CONFIG_PATH = "configs/fgc/fgc-1-tdg10-qa1-frz1.toml"
RESULT_PATH = "results/fgc-1-tdg10-qa1-frz1.json"
OWNER_DOCUMENT = "docs/fgc-tdg10-qa1-frz1.md"
AUTHORITY_PATH = (
    "src/recursive_horizons/fgc/evolution/tdg10_qa1_authority.py"
)
RUNNER_PATH = "scripts/run_fgc_tdg10_qa1.py"
EXACT_RUNTIME_PATH = (
    "src/recursive_horizons/fgc/evolution/tdg10_exact_complete_c_runtime.py"
)
EXACT_ADMISSION_PATH = (
    "src/recursive_horizons/fgc/evolution/tdg10_exact_complete_c_admission.py"
)
TI2_RUNNER_PATH = "scripts/run_fgc_tdg9_ti2.py"
UR1_PREF1_ARTIFACT_ID = "FGC-1-TDG9-UR1-PREF1"
UR1_PREF1_RESULT_PATH = "results/fgc-1-tdg9-ur1-pref1.json"
UR1_PREF1_RESULT_SHA256 = (
    "27c2aef9b033092285b868f583441842de99f61941ff57aa7e62749cc8d91300"
)
OUTPUT_NAMESPACE = (
    "runs/fgc-2-sf1/tdg10-qa1/retry3-ssprk3-sbp4-exact-complete-c"
)
STAGING_PREFIX = ".retry3-ssprk3-sbp4-exact-complete-c-qa1.stage-"
BASE_COMMIT = "49c514ac283ae3ec8091a6638442da6fdaf58db0"

STORE_PATH = "runs/fgc-2-sf1/tdg8-rcv3/calibration"
SEALED_STORE_LEAF_COUNT = 115
SEALED_STORE_SNAPSHOT_SHA256 = (
    "5917d70de3080dcc6b7abeb7fdb7cc3370c6140cbeb37bde0c142d6a2d36a445"
)
PREDECESSOR_GENERATION = 9
FORBIDDEN_RETRY3_GENERATION = 10
CHECKPOINT_SHA256 = (
    "eb6fddc480c94aa7ed15c80399fc2b26637693efef02399f267075fc6c258e56"
)
CHECKPOINT_RAW_SHA256 = (
    "07572d4450829ff2e2bad65e217a917b0d31ba677a3b8fc20a3be43e849ba846"
)
CHECKPOINT_JOURNAL_SEQUENCE = 10
CHECKPOINT_JOURNAL_SHA256 = (
    "5b533eb7009a9c7c3f353d26813f9cbbde8574afae77492929b435b4c701c5ff"
)
CHECKPOINT_JOURNAL_RAW_SHA256 = (
    "b50dc39ebd8bbdd3729d40d9e9ae4b223ff2ef19389d4ca74f9a316bac72202f"
)
HISTORICAL_RETRY3_REJECTION_SEQUENCE = 11
HISTORICAL_RETRY3_REJECTION_SHA256 = (
    "341cd8cd328434d85774bcb452889ac1886c520a72e945a92337b72cf564161a"
)
HISTORICAL_RETRY3_REJECTION_RAW_SHA256 = (
    "e4f9ca42f2d014d4bdbd875a0f24a9bbb9c04159fb0d9e0e7de8d0d72db61fdf"
)
MEMBER_KEY = "RK4-2049"
POINT_COUNT = 2049
OWNED_ROW_COUNT = 2044
RETRY = 3
PRIOR_RETRY_COUNT = 2
ACCEPTED_TIME_HEX = "0x1.78554de5a30e0p+0"
WIDTH_HEX = "0x1.aaa9612df8000p-11"
PHYSICAL_STATE_SHA256 = (
    "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a"
)
MEMBER_DESCRIPTOR_SHA256 = (
    "77847126340e78c8ac300fac2795bceda724c15e0894e34444b0da42a84166b4"
)
COORDINATES_SHA256 = (
    "2a347b314de5e1e94ac0c999db283e03498070ece6762c53252fa41c1e04a646"
)
TRANSACTION_SHA256 = (
    "7a90163d2eb4252fa7a1bbdf55d12f92a9127ed7cba5a37c28456e72d732b22b"
)
HISTORICAL_JOURNAL_SHA256 = (
    "0d21f650ccc46db778fd81fbebf5acab6394e1c2f7972b7940e24c76af9e1f59"
)
CURSOR_MODE = "RETRY_PENDING"
PENDING_OWNER = "temporal"
CHANNEL_ORDER = TDG6_COMPLETE_STATE_CHANNELS
CHANNEL_COUNT = 18
TABLEAU_SELECTOR = "SSPRK3"
TABLEAU_RUNTIME_SELECTOR = (
    "second_order_diagonal_norm_SBP_plus_SSPRK3"
)
ACTUAL_SPATIAL_OPERATOR = "inherited_RK4_2049_SBP4"
INTERVAL_OWNER = "exact_radius_free_complete_C_dual_rational_localizer"
PRIMARY_REFINEMENT_DEPTH = 160
PER_CHANNEL_MAXIMUM_CANDIDATES_D01 = 16352
PER_CHANNEL_MAXIMUM_CANDIDATES_D12 = 32704

REPLAY = dict(ti2.REPLAYS[0])
ENVIRONMENT = {
    "python_implementation": "CPython",
    "python_version": "3.14.3",
    "numpy_version": "2.5.1",
    "system": "Darwin",
    "machine": "arm64",
    "byteorder": "little",
}
IMPLEMENTATION_INVENTORY = (
    ("runner_path", RUNNER_PATH),
    ("authority_path", AUTHORITY_PATH),
    ("exact_runtime_path", EXACT_RUNTIME_PATH),
    ("exact_admission_path", EXACT_ADMISSION_PATH),
    ("TI2_runner_path", TI2_RUNNER_PATH),
)
IMPLEMENTATION_DIGEST_KEYS = tuple(
    key.replace("_path", "_sha256") for key, _path in IMPLEMENTATION_INVENTORY
)
RUNTIME_BYTE_PATHS = (
    CONFIG_PATH,
    RESULT_PATH,
    RUNNER_PATH,
    AUTHORITY_PATH,
)

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
            "scripts/reproduce_fgc_tdg10_qa1_frz1.py",
            RUNNER_PATH,
            AUTHORITY_PATH,
            EXACT_ADMISSION_PATH,
            EXACT_RUNTIME_PATH,
            "tests/test_check_repo_tdg10_qa1_frz1.py",
            "tests/test_fgc_tdg10_exact_complete_c_admission.py",
            "tests/test_fgc_tdg10_exact_complete_c_runtime.py",
            "tests/test_fgc_tdg10_qa1_authority.py",
            "tests/test_fgc_tdg10_qa1_runner.py",
        )
    )
)
ALLOWED_UNTRACKED = (".qdrant-initialized",)
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")


class QA1AuthorityError(RuntimeError):
    """The prospective QA1 image differs from its frozen qualification authority."""


def _fail(message: str) -> NoReturn:
    raise QA1AuthorityError(message)


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
        raise QA1AuthorityError("Git image differs") from exc


def _paths(root: Path, *arguments: str) -> tuple[str, ...]:
    raw = _git(root, *arguments)
    if not raw:
        return ()
    chunks = raw.split(b"\0")
    if chunks[-1] != b"":
        _fail("Git path stream differs")
    return tuple(item.decode() for item in chunks[:-1])


def _stat_identity(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_nlink)


def read_leaf(root: Path, relative: str) -> bytes:
    """Read one regular leaf without following any parent or leaf symlink."""

    candidate = Path(relative)
    if candidate.is_absolute() or any(
        part in {"", ".", ".."} for part in candidate.parts
    ):
        _fail("unsafe authority path")
    parts = candidate.parts
    flags_dir = (
        os.O_RDONLY
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    flags_leaf = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    directory_fd = -1
    leaf_fd = -1
    try:
        try:
            directory_fd = os.open(root, flags_dir)
        except OSError as exc:
            raise QA1AuthorityError("unsafe authority root") from exc
        root_stat = os.fstat(directory_fd)
        if stat.S_ISLNK(root_stat.st_mode) or not stat.S_ISDIR(root_stat.st_mode):
            _fail("unsafe authority root")
        for component in parts[:-1]:
            try:
                child_fd = os.open(component, flags_dir, dir_fd=directory_fd)
            except OSError as exc:
                raise QA1AuthorityError("authority parent is unsafe") from exc
            child_stat = os.fstat(child_fd)
            if stat.S_ISLNK(child_stat.st_mode) or not stat.S_ISDIR(child_stat.st_mode):
                os.close(child_fd)
                _fail("authority parent is unsafe")
            os.close(directory_fd)
            directory_fd = child_fd
        try:
            leaf_fd = os.open(parts[-1], flags_leaf, dir_fd=directory_fd)
        except FileNotFoundError as exc:
            raise QA1AuthorityError("authority leaf is absent") from exc
        except OSError as exc:
            raise QA1AuthorityError("unsafe authority leaf") from exc
        before = os.fstat(leaf_fd)
        if (
            stat.S_ISLNK(before.st_mode)
            or not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_size < 0
        ):
            _fail("unsafe authority leaf")
        remaining = before.st_size
        blocks: list[bytes] = []
        while remaining:
            block = os.read(leaf_fd, min(1 << 20, remaining))
            if not block:
                _fail("authority leaf changed")
            blocks.append(block)
            remaining -= len(block)
        if os.read(leaf_fd, 1):
            _fail("authority leaf changed")
        after = os.fstat(leaf_fd)
        if _stat_identity(before) != _stat_identity(after):
            _fail("authority leaf changed")
        raw = b"".join(blocks)
        if len(raw) != before.st_size:
            _fail("authority leaf size changed")
        return raw
    finally:
        if leaf_fd != -1:
            os.close(leaf_fd)
        if directory_fd != -1:
            os.close(directory_fd)


def _working(root: Path) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    return (
        _paths(root, "diff", "--name-only", "-z", "--"),
        _paths(root, "diff", "--cached", "--name-only", "-z", "--"),
        _paths(root, "ls-files", "--others", "--exclude-standard", "-z", "--"),
    )


def _require_no_replace_refs(root: Path) -> None:
    if _git(root, "for-each-ref", "--format=%(refname)", "refs/replace").strip():
        _fail("Git replace refs are present")


def _environment() -> dict[str, str]:
    import numpy as np

    observed = {
        "python_implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "system": platform.system(),
        "machine": platform.machine(),
        "byteorder": sys.byteorder,
    }
    if observed != ENVIRONMENT:
        _fail("execution environment differs")
    return observed


def _require_frozen_retry3_replay() -> None:
    if (
        int(REPLAY["retry"]) != RETRY
        or int(REPLAY["prior_retry_count"]) != PRIOR_RETRY_COUNT
        or int(REPLAY["predecessor_generation"]) != PREDECESSOR_GENERATION
        or int(REPLAY["predecessor_generation"]) == FORBIDDEN_RETRY3_GENERATION
        or str(REPLAY["checkpoint_sha256"]) != CHECKPOINT_SHA256
        or str(REPLAY["checkpoint_raw_sha256"]) != CHECKPOINT_RAW_SHA256
        or int(REPLAY["journal_sequence"]) != HISTORICAL_RETRY3_REJECTION_SEQUENCE
        or str(REPLAY["journal_sha256"]) != HISTORICAL_RETRY3_REJECTION_SHA256
        or str(REPLAY["journal_raw_sha256"]) != HISTORICAL_RETRY3_REJECTION_RAW_SHA256
        or str(REPLAY["attempted_width_hex"]) != WIDTH_HEX
        or str(REPLAY.get("transaction_sha256", TRANSACTION_SHA256))
        != TRANSACTION_SHA256
    ):
        _fail("frozen generation-9 retry-3 replay differs")


def _predecessor() -> dict[str, object]:
    return {
        "store": STORE_PATH,
        "store_leaf_count": SEALED_STORE_LEAF_COUNT,
        "store_snapshot_sha256": SEALED_STORE_SNAPSHOT_SHA256,
        "generation": PREDECESSOR_GENERATION,
        "checkpoint_sha256": CHECKPOINT_SHA256,
        "checkpoint_raw_sha256": CHECKPOINT_RAW_SHA256,
        "checkpoint_journal_sequence": CHECKPOINT_JOURNAL_SEQUENCE,
        "checkpoint_journal_sha256": CHECKPOINT_JOURNAL_SHA256,
        "checkpoint_journal_raw_sha256": CHECKPOINT_JOURNAL_RAW_SHA256,
        "historical_retry3_rejection_sequence": (
            HISTORICAL_RETRY3_REJECTION_SEQUENCE
        ),
        "historical_retry3_rejection_sha256": HISTORICAL_RETRY3_REJECTION_SHA256,
        "historical_retry3_rejection_raw_sha256": (
            HISTORICAL_RETRY3_REJECTION_RAW_SHA256
        ),
        "member_key": MEMBER_KEY,
        "point_count": POINT_COUNT,
        "owned_row_count": OWNED_ROW_COUNT,
        "retry": RETRY,
        "prior_retry_count": PRIOR_RETRY_COUNT,
        "accepted_time_binary64_hex": ACCEPTED_TIME_HEX,
        "attempted_width_binary64_hex": WIDTH_HEX,
        "physical_state_sha256": PHYSICAL_STATE_SHA256,
        "member_descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
        "coordinates_sha256": COORDINATES_SHA256,
        "transaction_sha256": TRANSACTION_SHA256,
        "historical_journal_sha256": HISTORICAL_JOURNAL_SHA256,
        "cursor_mode": CURSOR_MODE,
        "pending_owner": PENDING_OWNER,
        "immediate_predecessor_artifact_id": UR1_PREF1_ARTIFACT_ID,
        "immediate_predecessor_commit": BASE_COMMIT,
        "immediate_predecessor_result_path": UR1_PREF1_RESULT_PATH,
        "immediate_predecessor_result_sha256": UR1_PREF1_RESULT_SHA256,
    }


def _formulation() -> dict[str, object]:
    return {
        "source_member_method": "RK4",
        "source_member_integrator": (
            "fourth_order_diagonal_norm_SBP_with_second_order_boundary_closure_plus_RK4"
        ),
        "tableau_selector": TABLEAU_SELECTOR,
        "tableau_runtime_selector": TABLEAU_RUNTIME_SELECTOR,
        "actual_spatial_operator": ACTUAL_SPATIAL_OPERATOR,
        "successor_formulation_if_later_adopted": (
            "SSPRK3_on_inherited_RK4_2049_SBP4"
        ),
        "interval_owner": INTERVAL_OWNER,
        "channel_order_owner": "TDG6_COMPLETE_STATE_CHANNELS",
        "channel_count": CHANNEL_COUNT,
        "complete_admission_is_all_of": True,
        "minimum_observed_order": "3/2",
        "order_squared_multiplier": 8,
        "equality_passes": True,
        "outer_subintervals": 2,
        "finest_subintervals": 4,
        "primary_refinement_depth": PRIMARY_REFINEMENT_DEPTH,
        "per_channel_maximum_candidates_D01": PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
        "per_channel_maximum_candidates_D12": PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
        "independent_route_required": True,
        "absolute_tolerance_used": False,
        "physical_signal_used_for_normalization": False,
        "declared_cubic_is_exact_PDE_history": False,
        "pde_trajectory_order_certified": False,
    }


def _transaction() -> dict[str, object]:
    return {
        "shadow_path_count": 7,
        "shadow_proposal_count": 7,
        "SSPRK3_records_per_proposal": 4,
        "maximum_stage_and_endpoint_records": 28,
        "historical_store_read_only": True,
        "output_namespace_must_be_absent": True,
        "overwrite_forbidden": True,
        "diagnostic_fine_endpoint_serialized_iff_all_channels_pass": True,
        "diagnostic_fine_endpoint_is_accepted_state": False,
        "campaign_state_write_authorized": False,
        "PDE_state_commit_authorized": False,
        "fine_path_commit_authorized": False,
        "independent_postrun_binder_required": True,
        "later_state_adoption_requires_separate_freeze": True,
    }


def _exclusions() -> dict[str, object]:
    return {
        "retry_4_authorized": False,
        "retry_5_authorized": False,
        "fourth_width_authorized": False,
        "new_grid_authorized": False,
        "threshold_change_authorized": False,
        "historical_module_mutation_authorized": False,
        "common_event_authorized": False,
        "GR0_calibration_authorized": False,
        "production_SSPRK3_comparator_earned": False,
        "independent_method_agreement_earned": False,
        "SGBL_authorized": False,
        "FGCQR_authorized": False,
        "candidate_execution_authorized": False,
        "mechanism_result_authorized": False,
        "physical_result_authorized": False,
        "retained_EFT_evolution_authorized": False,
        "physical_transition_claim_authorized": False,
        "external_publication_authorized": False,
        "push_authorized": False,
    }


def _require_implementation_block(value: object) -> dict[str, object]:
    """Validate the exact path schema and config-supplied hashes. No self-hash constants."""

    if not isinstance(value, Mapping):
        _fail("QA1 implementation block is required")
    expected_keys = []
    for key, path in IMPLEMENTATION_INVENTORY:
        expected_keys.append(key)
        expected_keys.append(key.replace("_path", "_sha256"))
    if set(value) != set(expected_keys):
        _fail("QA1 implementation inventory keys differ")
    answer: dict[str, object] = {}
    for key, path in IMPLEMENTATION_INVENTORY:
        digest_key = key.replace("_path", "_sha256")
        declared = value.get(key)
        digest = value.get(digest_key)
        if declared != path:
            _fail(f"QA1 implementation path differs: {key}")
        if not isinstance(digest, str) or not _SHA256.fullmatch(digest):
            _fail(f"QA1 implementation hash is malformed: {digest_key}")
        answer[key] = path
        answer[digest_key] = digest
    return answer


def _expected_config(implementation: Mapping[str, object]) -> dict[str, object]:
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL,
        "owner_document": OWNER_DOCUMENT,
        "base_commit": BASE_COMMIT,
        "output_namespace": OUTPUT_NAMESPACE,
        "predecessor": _predecessor(),
        "formulation": _formulation(),
        "transaction": _transaction(),
        "environment": dict(ENVIRONMENT),
        "exclusions": _exclusions(),
        "implementation": dict(implementation),
    }


def parse_config(raw: bytes) -> dict[str, object]:
    import tomllib

    try:
        value = tomllib.loads(raw.decode())
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise QA1AuthorityError("invalid QA1 TOML") from exc
    _require_frozen_retry3_replay()
    implementation = _require_implementation_block(value.get("implementation"))
    expected = _expected_config(implementation)
    if value != expected:
        _fail("QA1 config differs from frozen contract")
    predecessor = value["predecessor"]
    if (
        predecessor["generation"] != PREDECESSOR_GENERATION
        or predecessor["generation"] == FORBIDDEN_RETRY3_GENERATION
        or predecessor["checkpoint_journal_sequence"]
        != CHECKPOINT_JOURNAL_SEQUENCE
        or predecessor["historical_retry3_rejection_sequence"]
        != HISTORICAL_RETRY3_REJECTION_SEQUENCE
        or predecessor["retry"] != RETRY
        or predecessor["prior_retry_count"] != PRIOR_RETRY_COUNT
        or predecessor["cursor_mode"] != CURSOR_MODE
        or predecessor["pending_owner"] != PENDING_OWNER
        or predecessor["immediate_predecessor_artifact_id"] != UR1_PREF1_ARTIFACT_ID
        or predecessor["immediate_predecessor_commit"] != BASE_COMMIT
        or predecessor["immediate_predecessor_result_path"] != UR1_PREF1_RESULT_PATH
        or predecessor["immediate_predecessor_result_sha256"] != UR1_PREF1_RESULT_SHA256
    ):
        _fail("QA1 retry-3 predecessor boundary differs")
    if value["target_protocol"] != TARGET_PROTOCOL:
        _fail("QA1 target protocol differs")
    if value["project_version"] != PROJECT_VERSION:
        _fail("QA1 project version differs")
    if any(value["exclusions"].values()) or any(
        value["transaction"][name] is not False
        for name in (
            "diagnostic_fine_endpoint_is_accepted_state",
            "campaign_state_write_authorized",
            "PDE_state_commit_authorized",
            "fine_path_commit_authorized",
        )
    ):
        _fail("QA1 authorized a state, campaign, candidate, or physics flag")
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
            _fail("QA1 output parent is unsafe")
    try:
        target.lstat()
    except FileNotFoundError:
        pass
    else:
        _fail("QA1 output namespace already exists")
    if any(item.name.startswith(STAGING_PREFIX) for item in target.parent.iterdir()):
        _fail("QA1 staging namespace already exists")


def filled_implementation_pairs(
    config: Mapping[str, object],
) -> tuple[tuple[str, str], ...]:
    """Return the required (path, sha256) pairs from the parsed implementation block."""

    block = _require_implementation_block(config.get("implementation"))
    return tuple(
        (path, str(block[key.replace("_path", "_sha256")]))
        for key, path in IMPLEMENTATION_INVENTORY
    )


def verify_implementation_inventory(
    root: Path,
    config: Mapping[str, object],
    *,
    commit: str | None,
) -> tuple[tuple[str, str], ...]:
    pairs = filled_implementation_pairs(config)
    if len(pairs) != len(IMPLEMENTATION_INVENTORY):
        _fail("QA1 implementation inventory is incomplete")
    for path, digest in pairs:
        raw = (
            read_leaf(root, path)
            if commit is None
            else _git(root, "show", f"{commit}:{path}")
        )
        if sha256(raw).hexdigest() != digest:
            _fail(f"QA1 implementation binding differs: {path}")
    return pairs


def _verify_ur1_pref1_result(root: Path) -> None:
    live = read_leaf(root, UR1_PREF1_RESULT_PATH)
    if sha256(live).hexdigest() != UR1_PREF1_RESULT_SHA256:
        _fail("live UR1-PREF1 compact result differs")
    committed = _git(root, "show", f"{BASE_COMMIT}:{UR1_PREF1_RESULT_PATH}")
    if sha256(committed).hexdigest() != UR1_PREF1_RESULT_SHA256:
        _fail("committed UR1-PREF1 compact result differs")
    if committed != live:
        _fail("UR1-PREF1 compact result drifted from the base commit")


def _verify_predecessor_leaves(root: Path) -> None:
    checkpoint = (
        f"{STORE_PATH}/checkpoints/"
        f"{PREDECESSOR_GENERATION:020d}-{CHECKPOINT_SHA256}.json"
    )
    journal_tip = (
        f"{STORE_PATH}/journal/"
        f"{CHECKPOINT_JOURNAL_SEQUENCE:020d}-{CHECKPOINT_JOURNAL_SHA256}.journal"
    )
    rejection = (
        f"{STORE_PATH}/journal/"
        f"{HISTORICAL_RETRY3_REJECTION_SEQUENCE:020d}-"
        f"{HISTORICAL_RETRY3_REJECTION_SHA256}.journal"
    )
    for relative, digest in (
        (checkpoint, CHECKPOINT_RAW_SHA256),
        (journal_tip, CHECKPOINT_JOURNAL_RAW_SHA256),
        (rejection, HISTORICAL_RETRY3_REJECTION_RAW_SHA256),
    ):
        if sha256(read_leaf(root, relative)).hexdigest() != digest:
            _fail(f"predecessor identity differs: {relative}")
    checkpoint_raw = json.loads(read_leaf(root, checkpoint))
    if (
        checkpoint_raw.get("generation") != PREDECESSOR_GENERATION
        or checkpoint_raw.get("checkpoint_sha256") != CHECKPOINT_SHA256
        or checkpoint_raw.get("journal_sequence") != CHECKPOINT_JOURNAL_SEQUENCE
        or checkpoint_raw.get("journal_tip_sha256") != CHECKPOINT_JOURNAL_SHA256
    ):
        _fail("generation-9 checkpoint tip is not journal sequence 10")
    journal_raw = json.loads(read_leaf(root, journal_tip))
    if (
        journal_raw.get("sequence") != CHECKPOINT_JOURNAL_SEQUENCE
        or journal_raw.get("record_sha256") != CHECKPOINT_JOURNAL_SHA256
        or journal_raw.get("kind") == "tdg6_rejection"
    ):
        _fail("checkpoint journal sequence 10 is not the generation-9 tip")
    rejection_raw = json.loads(read_leaf(root, rejection))
    if (
        rejection_raw.get("sequence") != HISTORICAL_RETRY3_REJECTION_SEQUENCE
        or rejection_raw.get("record_sha256") != HISTORICAL_RETRY3_REJECTION_SHA256
        or rejection_raw.get("kind") != "tdg6_rejection"
    ):
        _fail("sequence 11 is not the historical retry-3 rejection")


def _snapshot_store(root: Path) -> tuple[int, str]:
    from scripts import run_fgc_tdg9_ti2 as runner

    try:
        return runner._snapshot_store(root)
    except Exception as exc:
        raise QA1AuthorityError("sealed store cannot be authenticated") from exc


def _expected_fingerprint_receipt() -> dict[str, object]:
    return {
        "retry": RETRY,
        "predecessor_generation": PREDECESSOR_GENERATION,
        "checkpoint_journal_sequence": CHECKPOINT_JOURNAL_SEQUENCE,
        "historical_retry3_rejection_sequence": (
            HISTORICAL_RETRY3_REJECTION_SEQUENCE
        ),
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
        "predecessor_generation": PREDECESSOR_GENERATION,
        "checkpoint_journal_sequence": CHECKPOINT_JOURNAL_SEQUENCE,
        "historical_retry3_rejection_sequence": (
            HISTORICAL_RETRY3_REJECTION_SEQUENCE
        ),
        "attempted_width_hex": WIDTH_HEX,
        "accepted_time_hex": ACCEPTED_TIME_HEX,
        "state_sha256": PHYSICAL_STATE_SHA256,
        "descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
        "transaction_sha256": TRANSACTION_SHA256,
        "historical_journal_sha256": HISTORICAL_JOURNAL_SHA256,
    }


def real_predecessor_fingerprints(root: Path) -> tuple[dict[str, object], ...]:
    """Restore only generation 9 / retry 3 and fingerprint without proposals."""

    from scripts import run_fgc_tdg9_ti2 as runner

    from . import tdg9_ar1_authority as ar1
    from .hlt16_campaign_store import HLT16CampaignStore
    from .proto19_gr0_static_factory import build_static_gr0_shells

    _require_frozen_retry3_replay()
    if int(REPLAY["predecessor_generation"]) == FORBIDDEN_RETRY3_GENERATION:
        _fail("retry-3 must not restore generation 10")
    repository = root.resolve()
    store = HLT16CampaignStore(repository / ar1.PREF2_STORE_PATH)
    shells = build_static_gr0_shells(repository)
    try:
        restored = runner._restore_replay(repository, store, shells, REPLAY)
        fingerprint = restored.fingerprint
    except Exception as exc:
        raise QA1AuthorityError("retry-3 predecessor fingerprint failed") from exc
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
            "predecessor_generation": PREDECESSOR_GENERATION,
            "checkpoint_journal_sequence": CHECKPOINT_JOURNAL_SEQUENCE,
            "historical_retry3_rejection_sequence": (
                HISTORICAL_RETRY3_REJECTION_SEQUENCE
            ),
            "transaction_sha256": fingerprint["transaction_sha256"],
            "state_sha256": fingerprint["state_sha256"],
            "accepted_time_hex": fingerprint["accepted_time_hex"],
            "descriptor_sha256": fingerprint["descriptor_sha256"],
            "shadow_proposal_constructed": False,
        },
    )


def expected_compact(config_raw: bytes) -> dict[str, object]:
    parsed = parse_config(config_raw)
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "gate_status": "pass",
        "source_config_sha256": sha256(config_raw).hexdigest(),
        "artifact_payload": {
            **_expected_config(parsed["implementation"]),
            "store_identity_observed_at_prelaunch": {
                "leaf_count": SEALED_STORE_LEAF_COUNT,
                "sha256": SEALED_STORE_SNAPSHOT_SHA256,
            },
            "QA1_output_namespace_absent_at_prelaunch": True,
            "replace_refs_absent": True,
            "restored_predecessor_fingerprints": [
                _expected_fingerprint_receipt()
            ],
            "shadow_executed": False,
            "shadow_proposals_constructed": 0,
            "diagnostic_qualification_only": True,
            "new_protocol_id_authorized": False,
            "target_protocol": TARGET_PROTOCOL,
            "retry3_restores_generation": PREDECESSOR_GENERATION,
            "retry3_never_restores_generation_10": True,
            "checkpoint_journal_is_sequence_10": True,
            "historical_retry3_rejection_is_sequence_11": True,
            "state_advance_authorized": False,
            "campaign_state_write_authorized": False,
            "PDE_state_commit_authorized": False,
            "fine_path_commit_authorized": False,
            "candidate_execution_authorized": False,
            "mechanism_result_authorized": False,
            "physical_result_authorized": False,
        },
    }


def build_prelaunch(
    config_raw: bytes,
    root: Path,
    *,
    store_snapshot: tuple[int, str],
) -> dict[str, object]:
    repository = root.resolve()
    config = parse_config(config_raw)
    _environment()
    head = _git(repository, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    if head != BASE_COMMIT:
        _fail("prelaunch HEAD is not the sealed QA1 base commit")
    _require_no_replace_refs(repository)
    unstaged, staged, untracked = _working(repository)
    observed = tuple(
        sorted(
            (
                *unstaged,
                *(item for item in untracked if item not in ALLOWED_UNTRACKED),
            )
        )
    )
    if staged or observed != DELTA_PATHS:
        _fail("prospective QA1 delta differs")
    for relative in DELTA_PATHS:
        read_leaf(repository, relative)
    verify_implementation_inventory(repository, config, commit=None)
    _verify_ur1_pref1_result(repository)
    _verify_predecessor_leaves(repository)
    sealed = (SEALED_STORE_LEAF_COUNT, SEALED_STORE_SNAPSHOT_SHA256)
    if store_snapshot != sealed:
        _fail("sealed store differs at QA1 prelaunch")
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
            "restored_predecessor_fingerprints": fingerprints,
        },
    }


def _reject_compact_json_constant(token: str) -> NoReturn:
    _fail(f"nonfinite QA1 compact JSON constant: {token}")


def validate_compact(config_raw: bytes, result_raw: bytes) -> dict[str, object]:
    expected = expected_compact(config_raw)
    try:
        result = json.loads(
            result_raw,
            object_pairs_hook=_unique,
            parse_constant=_reject_compact_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise QA1AuthorityError("invalid QA1 compact JSON") from exc
    if result != expected or result_raw != canonical_pretty(result):
        _fail("QA1 compact result differs")
    payload = result["artifact_payload"]
    if payload["shadow_executed"] is not False or payload[
        "shadow_proposals_constructed"
    ] != 0:
        _fail("QA1 compact result claimed shadow work")
    if payload["state_advance_authorized"] is not False:
        _fail("QA1 compact result authorized a state advance")
    for forbidden in (
        "complete_admission_passed",
        "failed_channels",
        "diagnostic_payload_sha256",
        "fine_endpoint_descriptor_sha256",
    ):
        if forbidden in payload:
            _fail("QA1 compact result claimed a diagnostic outcome")
    return result


def _require_bytes_match_commit(root: Path, commit: str, relative: str) -> None:
    if _git(root, "show", f"{commit}:{relative}") != read_leaf(root, relative):
        _fail(f"working {relative} differs from the authority commit")


@dataclass(frozen=True, slots=True)
class QA1Authority:
    authority_commit: str
    execution_authorized: bool = True
    diagnostic_qualification_only: bool = True
    tableau_selector: str = TABLEAU_SELECTOR
    actual_spatial_operator: str = ACTUAL_SPATIAL_OPERATOR
    retry: int = RETRY
    predecessor_generation: int = PREDECESSOR_GENERATION
    checkpoint_journal_sequence: int = CHECKPOINT_JOURNAL_SEQUENCE
    historical_retry3_rejection_sequence: int = (
        HISTORICAL_RETRY3_REJECTION_SEQUENCE
    )
    target_protocol: str = TARGET_PROTOCOL
    new_protocol_id_authorized: bool = False
    state_advance_authorized: bool = False
    campaign_state_write_authorized: bool = False
    PDE_state_commit_authorized: bool = False
    fine_path_commit_authorized: bool = False
    candidate_execution_authorized: bool = False
    mechanism_result_authorized: bool = False
    physical_result_authorized: bool = False
    production_SSPRK3_comparator_earned: bool = False
    independent_method_agreement_earned: bool = False


def authorize(
    root: Path,
    config_raw: bytes,
    result_raw: bytes,
    authority_commit: str,
) -> QA1Authority:
    repository = root.resolve()
    if not _COMMIT.fullmatch(authority_commit):
        _fail("QA1 authority commit is malformed")
    config = validate_compact(config_raw, result_raw)
    _environment()
    head = _git(repository, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    parents = (
        _git(repository, "rev-list", "--parents", "-n", "1", authority_commit)
        .decode()
        .split()[1:]
    )
    if head != authority_commit or parents != [BASE_COMMIT]:
        _fail("QA1 authority is not the exact direct successor")
    _require_no_replace_refs(repository)
    unstaged, staged, untracked = _working(repository)
    if unstaged or staged or any(item not in ALLOWED_UNTRACKED for item in untracked):
        _fail("QA1 authority working image is not clean")
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
        _fail("committed QA1 delta differs")
    if _git(repository, "show", f"{authority_commit}:{CONFIG_PATH}") != config_raw:
        _fail("committed QA1 config differs")
    if _git(repository, "show", f"{authority_commit}:{RESULT_PATH}") != result_raw:
        _fail("committed QA1 compact result differs")
    for relative in RUNTIME_BYTE_PATHS:
        _require_bytes_match_commit(repository, authority_commit, relative)
    verify_implementation_inventory(
        repository, config["artifact_payload"], commit=authority_commit
    )
    _verify_ur1_pref1_result(repository)
    _verify_predecessor_leaves(repository)
    if _snapshot_store(repository) != (
        SEALED_STORE_LEAF_COUNT,
        SEALED_STORE_SNAPSHOT_SHA256,
    ):
        _fail("sealed store differs before QA1 authorization")
    if list(real_predecessor_fingerprints(repository)) != [
        _expected_fingerprint_receipt()
    ]:
        _fail("restored predecessor fingerprints differ before authorization")
    return QA1Authority(authority_commit=authority_commit)


__all__ = [
    "ARTIFACT_ID",
    "AUTHORITY_PATH",
    "BASE_COMMIT",
    "CHANNEL_ORDER",
    "CHECKPOINT_JOURNAL_SEQUENCE",
    "CHECKPOINT_SHA256",
    "CLASSIFICATION",
    "CONFIG_PATH",
    "DELTA_PATHS",
    "FORBIDDEN_RETRY3_GENERATION",
    "HISTORICAL_JOURNAL_SHA256",
    "HISTORICAL_RETRY3_REJECTION_SEQUENCE",
    "IMPLEMENTATION_INVENTORY",
    "MEMBER_KEY",
    "PROJECT_VERSION",
    "UR1_PREF1_ARTIFACT_ID",
    "UR1_PREF1_RESULT_PATH",
    "UR1_PREF1_RESULT_SHA256",
    "OUTPUT_NAMESPACE",
    "PREDECESSOR_GENERATION",
    "QA1Authority",
    "QA1AuthorityError",
    "REPLAY",
    "RESULT_PATH",
    "RETRY",
    "RUNNER_PATH",
    "STAGING_PREFIX",
    "STORE_PATH",
    "TARGET_PROTOCOL",
    "authorize",
    "build_prelaunch",
    "canonical_pretty",
    "expected_compact",
    "expected_raw_replay_receipt",
    "filled_implementation_pairs",
    "parse_config",
    "read_leaf",
    "real_predecessor_fingerprints",
    "require_output_absent",
    "validate_compact",
    "verify_implementation_inventory",
]
