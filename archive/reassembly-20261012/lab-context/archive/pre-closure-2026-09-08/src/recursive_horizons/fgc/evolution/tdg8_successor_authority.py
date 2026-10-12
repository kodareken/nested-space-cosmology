"""Committed-image authority for one fresh TDG8 successor event.

The sealed PRO19 campaign is evidence, not a restart target.  This module
binds one prospective, GR-0-only successor image to the immutable PREF28,
TDG8-DIAG1, PREF27, and AUTH1 records and to the exact implementation bytes
that may create a new PROTO19 namespace.

Destination absence is a one-time prelaunch observation.  Compact validation
and execution authorization deliberately do not inspect either campaign
store, so a successfully created destination cannot invalidate the authority
that permitted its creation.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version
import json
import os
from pathlib import Path
import platform
import stat
import subprocess
import sys
import tomllib
from types import MappingProxyType
from typing import Any, Mapping


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG8-RUN1-AUTH1"
PROJECT_VERSION = "0.11.0"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CLASSIFICATION = "premise_only_committed_fresh_gr0_successor_authority"
CONFIG_PATH = "configs/fgc/fgc-1-tdg8-run1-auth1.toml"
RESULT_PATH = "results/fgc-1-tdg8-run1-auth1.json"

PREDECESSOR_REPAIR_COMMIT = "f1e96c5e3406705bbf3354a18e5458334e923119"
CAMPAIGN_ID = "FGC-2-SF1-PROTO18-GR0-A3-TDG8-SUCCESSOR-1"
SOURCE_CAMPAIGN_PATH = "runs/fgc-2-sf1/proto17/calibration"
DESTINATION_CAMPAIGN_PATH = "runs/fgc-2-sf1/proto19/calibration"

_MAX_CONFIG_BYTES = 1024 * 1024
_MAX_RESULT_BYTES = 4 * 1024 * 1024
_MAX_BOUND_FILE_BYTES = 32 * 1024 * 1024


@dataclass(frozen=True, slots=True)
class ImmutableBinding:
    """One historically sealed result consumed by the successor authority."""

    name: str
    artifact_id: str
    path: str
    sha256: str


@dataclass(frozen=True, slots=True)
class ImplementationBinding:
    """One exact implementation leaf in the authorized execution image."""

    role: str
    path: str
    sha256: str


IMMUTABLE_BINDINGS = (
    ImmutableBinding(
        name="pref28",
        artifact_id="FGC-1-PRO19-PREF28",
        path="results/fgc-1-pro19-pref28.json",
        sha256="c53482e7965e5e3ed79bad5abcd7323a28868801ef82392360916ae8979659e2",
    ),
    ImmutableBinding(
        name="tdg8_diag1",
        artifact_id="FGC-1-TDG8-DIAG1",
        path="results/fgc-1-tdg8-diag1.json",
        sha256="9605fd52383cd5722ad8d48e41b2ebca4a3ce6415eaa91cd6b5e5715106c896f",
    ),
    ImmutableBinding(
        name="pref27",
        artifact_id="FGC-1-PRO18-PREF27",
        path="results/fgc-1-pro18-pref27.json",
        sha256="102b385d058066d668e9e8541c974c736c583f6d5c4afaf9b99e95a16aabeb83",
    ),
    ImmutableBinding(
        name="auth1",
        artifact_id="FGC-1-PRO18-AUTH1",
        path="results/fgc-1-pro18-auth1.json",
        sha256="dd967568fbae0bd19e063af94ea8e2d53480c98afec5c78bd86601c20faf4efa",
    ),
)

# The hashes are supplied by the materialized config, but the closure's path
# and role set is code-owned.  A self-consistent config cannot omit a runtime
# owner or replace it with an adjacent implementation.
REQUIRED_IMPLEMENTATION_INVENTORY = (
    ("reproducer", "scripts/reproduce_fgc_tdg8_run1_auth1.py"),
    ("runner", "scripts/run_fgc_tdg8_successor_event.py"),
    (
        "runtime",
        "src/recursive_horizons/fgc/evolution/hlt16_campaign_runtime.py",
    ),
    (
        "runtime",
        "src/recursive_horizons/fgc/evolution/hlt16_campaign_store.py",
    ),
    (
        "runtime",
        "src/recursive_horizons/fgc/evolution/hlt16_lifecycle.py",
    ),
    (
        "runtime",
        "src/recursive_horizons/fgc/evolution/hlt16_member_codec.py",
    ),
    (
        "runtime",
        "src/recursive_horizons/fgc/evolution/hlt16_progression_attempt.py",
    ),
    (
        "authority",
        "src/recursive_horizons/fgc/evolution/tdg8_successor_authority.py",
    ),
    (
        "runtime",
        "src/recursive_horizons/fgc/evolution/tdg8_successor_runtime.py",
    ),
)

# The successor authority is intentionally much smaller than another complete
# import-closure capture.  Every execution dependency not listed here must be
# byte-for-byte inherited from PREDECESSOR_REPAIR_COMMIT.  Requiring the exact
# Git delta makes that inheritance a machine-checked statement: an adjacent
# source edit cannot hide behind a clean, self-consistent later commit, and a
# missing promised owner cannot silently disappear from the authority image.
AUTHORITY_DELTA_PATHS = tuple(
    sorted(
        (
            "Makefile",
            CONFIG_PATH,
            RESULT_PATH,
            "scripts/check_repo.py",
            "scripts/reproduce_fgc_tdg8_run1_auth1.py",
            "scripts/run_fgc_tdg8_successor_event.py",
            "src/recursive_horizons/fgc/evolution/hlt16_campaign_runtime.py",
            "src/recursive_horizons/fgc/evolution/hlt16_campaign_store.py",
            "src/recursive_horizons/fgc/evolution/hlt16_lifecycle.py",
            "src/recursive_horizons/fgc/evolution/hlt16_member_codec.py",
            "src/recursive_horizons/fgc/evolution/tdg8_successor_authority.py",
            "src/recursive_horizons/fgc/evolution/tdg8_successor_runtime.py",
            "tests/test_fgc_hlt16_campaign_store.py",
            "tests/test_fgc_tdg8_successor_authority.py",
            "tests/test_fgc_tdg8_successor_runner.py",
            "tests/test_fgc_tdg8_successor_runtime.py",
            "tests/test_hlt16_lifecycle.py",
        )
    )
)

PINNED_ENVIRONMENT: dict[str, str] = {
    "implementation": "CPython",
    "python_version": "3.14.3",
    "numpy_version": "2.5.1",
    "system": "Darwin",
    "machine": "arm64",
    "platform": "darwin",
    "byteorder": "little",
}

PROGRESSION_CONTRACT: dict[str, object] = {
    "protocol": TARGET_PROTOCOL,
    "campaign_id": CAMPAIGN_ID,
    "branch": "GR-0",
    "amplitude": "3",
    "event": 23,
    "start_rational": "23/16",
    "start_binary64_hex": "0x1.7000000000000p+0",
    "target_rational": "3/2",
    "target_binary64_hex": "0x1.8000000000000p+0",
    "event_count": 1,
    "candidate_branches_forbidden": True,
}

STORES_CONTRACT: dict[str, object] = {
    "source_path": SOURCE_CAMPAIGN_PATH,
    "source_disposition": "terminal_read_only",
    "source_resume_forbidden": True,
    "destination_path": DESTINATION_CAMPAIGN_PATH,
    "destination_initially_absent_required": True,
    "destination_creation_mode": "fresh_authenticated_generation_zero",
}

SCOPE_CONTRACT: dict[str, object] = {
    "one_common_event_only": True,
    "candidate_branches_forbidden": True,
    "locked_source_mutation_authorized": False,
    "locked_source_resume_authorized": False,
    "postexecution_destination_absence_reobservation": False,
    "threshold_or_retry_rule_change_authorized": False,
    "TDG6_or_TDG7_rule_change_authorized": False,
    "physical_input_change_authorized": False,
}

CLAIMS_CONTRACT: dict[str, object] = {
    "predecessor_repair_bound": True,
    "immutable_evidence_bound": True,
    "committed_execution_image_required": True,
    "successor_execution_authorized_by_artifact_alone": False,
    "source_terminal_retroactively_reclassified": False,
    "common_event_completed": False,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "physical_result_earned": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
}

NONCLAIMS = (
    "This authority does not mutate, resume, roll back, or reclassify the "
    "locked PROTO17 campaign.",
    "This authority permits at most one fresh GR-0 common event and opens no "
    "SGB-L or FGC-QR candidate.",
    "This authority records no common event, calibration, trapping, DEF1, retained-EFT, transition, mechanism, or physical result.",
)


class TDG8SuccessorAuthorityError(RuntimeError):
    """The prelaunch record, committed image, or fixed successor scope differs."""


@dataclass(frozen=True, slots=True)
class TDG8SuccessorExecutionAuthority:
    """Typed permission for the runner's single fresh common event."""

    authority_commit: str
    predecessor_repair_commit: str
    config_sha256: str
    result_sha256: str
    campaign_id: str
    target_protocol: str
    branch: str
    amplitude: str
    event: int
    start_rational: str
    target_rational: str
    event_count: int
    source_campaign_path: str
    destination_campaign_path: str
    source_terminal_read_only: bool
    candidate_branches_forbidden: bool
    immutable_bindings: tuple[ImmutableBinding, ...]
    implementation_inventory: tuple[ImplementationBinding, ...]
    environment: Mapping[str, str]


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def canonical_result(value: object) -> bytes:
    """Encode an authority result as deterministic repository JSON."""

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
        raise TDG8SuccessorAuthorityError(
            "TDG8 successor result is not canonical JSON"
        ) from exc


def _digest(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or value != value.lower():
        raise TDG8SuccessorAuthorityError(f"{label} digest differs")
    try:
        int(value, 16)
    except ValueError as exc:
        raise TDG8SuccessorAuthorityError(f"{label} digest differs") from exc
    return value


def _commit(value: object) -> str:
    if not isinstance(value, str) or len(value) != 40 or value != value.lower():
        raise TDG8SuccessorAuthorityError(
            "authority commit is not a lowercase 40-character Git id"
        )
    try:
        int(value, 16)
    except ValueError as exc:
        raise TDG8SuccessorAuthorityError(
            "authority commit is not hexadecimal"
        ) from exc
    return value


def _relative_path(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or "\\" in value or "\x00" in value:
        raise TDG8SuccessorAuthorityError(f"{label} is not a portable relative path")
    path = Path(value)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise TDG8SuccessorAuthorityError(f"{label} escapes the repository")
    return path.as_posix()


def _exact(value: object, expected: object) -> bool:
    """Type-sensitive recursive equality for TOML-derived contracts."""

    if isinstance(expected, Mapping):
        return (
            isinstance(value, Mapping)
            and set(value) == set(expected)
            and all(_exact(value[key], expected[key]) for key in expected)
        )
    if isinstance(expected, list):
        return (
            isinstance(value, list)
            and len(value) == len(expected)
            and all(_exact(actual, wanted) for actual, wanted in zip(value, expected))
        )
    return type(value) is type(expected) and value == expected


def _binding_rows() -> list[dict[str, str]]:
    return [
        {
            "name": item.name,
            "artifact_id": item.artifact_id,
            "path": item.path,
            "sha256": item.sha256,
        }
        for item in IMMUTABLE_BINDINGS
    ]


def _predecessor_contract() -> dict[str, object]:
    return {
        "repair_commit": PREDECESSOR_REPAIR_COMMIT,
        "must_be_ancestor_of_authority_commit": True,
    }


def _validate_inventory(value: object) -> list[dict[str, str]]:
    if not isinstance(value, list) or len(value) != len(
        REQUIRED_IMPLEMENTATION_INVENTORY
    ):
        raise TDG8SuccessorAuthorityError("implementation inventory closure differs")
    expected_roles = {path: role for role, path in REQUIRED_IMPLEMENTATION_INVENTORY}
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for item in value:
        if not isinstance(item, Mapping) or set(item) != {"role", "path", "sha256"}:
            raise TDG8SuccessorAuthorityError("implementation inventory entry differs")
        path = _relative_path(item.get("path"), "implementation inventory path")
        role = item.get("role")
        digest = _digest(item.get("sha256"), "implementation inventory")
        if path in seen or expected_roles.get(path) != role:
            raise TDG8SuccessorAuthorityError(
                "implementation inventory path or role differs"
            )
        seen.add(path)
        rows.append({"role": str(role), "path": path, "sha256": digest})
    if set(seen) != set(expected_roles) or rows != sorted(
        rows, key=lambda row: row["path"]
    ):
        raise TDG8SuccessorAuthorityError("implementation inventory closure differs")
    return rows


def _validate_config_mapping(value: Mapping[str, Any]) -> dict[str, Any]:
    required = {
        "schema_version",
        "artifact_id",
        "project_version",
        "target_protocol",
        "classification",
        "nonclaims",
        "predecessor",
        "immutable_bindings",
        "implementation_inventory",
        "environment",
        "progression",
        "stores",
        "scope",
        "claims",
    }
    if set(value) != required:
        raise TDG8SuccessorAuthorityError("TDG8 successor config fields differ")
    fixed: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL,
        "classification": CLASSIFICATION,
        "nonclaims": list(NONCLAIMS),
        "predecessor": _predecessor_contract(),
        "immutable_bindings": _binding_rows(),
        "environment": PINNED_ENVIRONMENT,
        "progression": PROGRESSION_CONTRACT,
        "stores": STORES_CONTRACT,
        "scope": SCOPE_CONTRACT,
        "claims": CLAIMS_CONTRACT,
    }
    for key, expected in fixed.items():
        if not _exact(value.get(key), expected):
            raise TDG8SuccessorAuthorityError(f"TDG8 successor {key} contract differs")
    checked = dict(value)
    checked["implementation_inventory"] = _validate_inventory(
        value.get("implementation_inventory")
    )
    return checked


def _parse_config(raw: bytes) -> dict[str, Any]:
    if len(raw) > _MAX_CONFIG_BYTES:
        raise TDG8SuccessorAuthorityError(
            "TDG8 successor config exceeds its byte budget"
        )
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise TDG8SuccessorAuthorityError("TDG8 successor config is malformed") from exc
    if not isinstance(value, Mapping):
        raise TDG8SuccessorAuthorityError("TDG8 successor config is not a table")
    return _validate_config_mapping(value)


def _json_object(raw: bytes, label: str) -> dict[str, Any]:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(key)
            result[key] = value
        return result

    try:
        value = json.loads(raw.decode("utf-8"), object_pairs_hook=reject_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TDG8SuccessorAuthorityError(f"{label} is malformed") from exc
    if not isinstance(value, dict):
        raise TDG8SuccessorAuthorityError(f"{label} is not object-valued")
    return value


def _parse_result(raw: bytes) -> dict[str, Any]:
    if len(raw) > _MAX_RESULT_BYTES:
        raise TDG8SuccessorAuthorityError(
            "TDG8 successor result exceeds its byte budget"
        )
    value = _json_object(raw, "TDG8 successor result")
    if canonical_result(value) != raw:
        raise TDG8SuccessorAuthorityError("TDG8 successor result is noncanonical")
    return value


def _read_repository_leaf(
    root: Path,
    relative: str,
    *,
    label: str,
    maximum_bytes: int = _MAX_BOUND_FILE_BYTES,
) -> bytes:
    """Read one bounded unique regular file through a no-follow FD walk."""

    relative = _relative_path(relative, label)
    components = Path(relative).parts
    directory_fd = -1
    leaf_fd = -1
    try:
        directory_fd = os.open(
            root,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
        root_stat = os.fstat(directory_fd)
        if not stat.S_ISDIR(root_stat.st_mode):
            raise TDG8SuccessorAuthorityError("repository root is unsafe")
        for component in components[:-1]:
            child_fd = os.open(
                component,
                os.O_RDONLY
                | getattr(os, "O_DIRECTORY", 0)
                | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=directory_fd,
            )
            child_stat = os.fstat(child_fd)
            if not stat.S_ISDIR(child_stat.st_mode):
                os.close(child_fd)
                raise TDG8SuccessorAuthorityError(f"{label} has an unsafe ancestor")
            os.close(directory_fd)
            directory_fd = child_fd
        leaf_fd = os.open(
            components[-1],
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=directory_fd,
        )
        before = os.fstat(leaf_fd)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_size < 0
            or before.st_size > maximum_bytes
        ):
            raise TDG8SuccessorAuthorityError(
                f"{label} is not a bounded unique regular file"
            )
        remaining = before.st_size
        chunks: list[bytes] = []
        while remaining:
            block = os.read(leaf_fd, min(1024 * 1024, remaining))
            if not block:
                raise TDG8SuccessorAuthorityError(f"{label} ended during read")
            chunks.append(block)
            remaining -= len(block)
        if os.read(leaf_fd, 1):
            raise TDG8SuccessorAuthorityError(f"{label} grew during read")
        after = os.fstat(leaf_fd)

        def identity(item: os.stat_result) -> tuple[int, int, int, int, int, int]:
            return (
                item.st_dev,
                item.st_ino,
                item.st_mode,
                item.st_nlink,
                item.st_size,
                item.st_mtime_ns,
            )

        if identity(before) != identity(after):
            raise TDG8SuccessorAuthorityError(f"{label} changed during read")
        return b"".join(chunks)
    except TDG8SuccessorAuthorityError:
        raise
    except FileNotFoundError as exc:
        raise TDG8SuccessorAuthorityError(f"{label} is absent") from exc
    except OSError as exc:
        raise TDG8SuccessorAuthorityError(f"{label} cannot be read safely") from exc
    finally:
        if leaf_fd != -1:
            os.close(leaf_fd)
        if directory_fd != -1:
            os.close(directory_fd)


def _require_path_absent(root: Path, relative: str, *, label: str) -> None:
    """Prove absence without following any existing ancestor or leaf link."""

    relative = _relative_path(relative, label)
    components = Path(relative).parts
    directory_fd = -1
    try:
        directory_fd = os.open(
            root,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
        if not stat.S_ISDIR(os.fstat(directory_fd).st_mode):
            raise TDG8SuccessorAuthorityError("repository root is unsafe")
        for component in components[:-1]:
            try:
                child_fd = os.open(
                    component,
                    os.O_RDONLY
                    | getattr(os, "O_DIRECTORY", 0)
                    | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=directory_fd,
                )
            except FileNotFoundError:
                return
            if not stat.S_ISDIR(os.fstat(child_fd).st_mode):
                os.close(child_fd)
                raise TDG8SuccessorAuthorityError(f"{label} has an unsafe ancestor")
            os.close(directory_fd)
            directory_fd = child_fd
        try:
            leaf_stat = os.stat(
                components[-1], dir_fd=directory_fd, follow_symlinks=False
            )
        except FileNotFoundError:
            return
        if stat.S_ISLNK(leaf_stat.st_mode):
            raise TDG8SuccessorAuthorityError(f"{label} is an unsafe link")
        raise TDG8SuccessorAuthorityError(f"{label} is already present")
    except TDG8SuccessorAuthorityError:
        raise
    except OSError as exc:
        raise TDG8SuccessorAuthorityError(f"{label} cannot be probed safely") from exc
    finally:
        if directory_fd != -1:
            os.close(directory_fd)


def _git(
    root: Path, *args: str, check: bool = True
) -> subprocess.CompletedProcess[bytes]:
    environment = {
        key: value for key, value in os.environ.items() if not key.startswith("GIT_")
    }
    environment.update({"GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C"})
    result = subprocess.run(
        (
            "git",
            "--no-replace-objects",
            "--no-optional-locks",
            "-c",
            "core.fsmonitor=false",
            "-c",
            "core.untrackedCache=false",
            "-C",
            str(root),
            *args,
        ),
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        env=environment,
    )
    if check and result.returncode:
        raise TDG8SuccessorAuthorityError(
            f"Git authority command failed: {' '.join(args)}"
        )
    return result


def _verify_repository_root(root: Path) -> None:
    top = _git(root, "rev-parse", "--show-toplevel").stdout.decode("utf-8").strip()
    if os.path.realpath(top) != os.path.realpath(root):
        raise TDG8SuccessorAuthorityError("authority root is not the Git worktree root")


def _require_predecessor_ancestry(root: Path, descendant: str) -> None:
    if _git(
        root,
        "merge-base",
        "--is-ancestor",
        PREDECESSOR_REPAIR_COMMIT,
        descendant,
        check=False,
    ).returncode:
        raise TDG8SuccessorAuthorityError(
            "predecessor repair commit is not an ancestor of the authority image"
        )


def _require_exact_authority_delta(root: Path, descendant: str) -> None:
    """Require the committed successor image to change exactly its owners."""

    raw = _git(
        root,
        "diff",
        "--no-renames",
        "--name-only",
        "-z",
        PREDECESSOR_REPAIR_COMMIT,
        descendant,
        "--",
    ).stdout
    if raw and not raw.endswith(b"\0"):
        raise TDG8SuccessorAuthorityError("authority Git delta is malformed")
    encoded = raw[:-1].split(b"\0") if raw else []
    observed: list[str] = []
    try:
        for item in encoded:
            if not item:
                raise TDG8SuccessorAuthorityError(
                    "authority Git delta contains an empty path"
                )
            observed.append(
                _relative_path(item.decode("utf-8"), "authority Git delta path")
            )
    except UnicodeDecodeError as exc:
        raise TDG8SuccessorAuthorityError(
            "authority Git delta path is not UTF-8"
        ) from exc
    if len(set(observed)) != len(observed):
        raise TDG8SuccessorAuthorityError("authority Git delta repeats a path")
    if tuple(sorted(observed)) != AUTHORITY_DELTA_PATHS:
        raise TDG8SuccessorAuthorityError(
            "authority Git delta differs from the exact successor owner set"
        )


def observed_environment() -> dict[str, str]:
    """Return the standard-library-visible numerical process identity."""

    try:
        numpy_version = version("numpy")
    except PackageNotFoundError as exc:
        raise TDG8SuccessorAuthorityError(
            "NumPy distribution metadata is absent"
        ) from exc
    return {
        "implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "numpy_version": numpy_version,
        "system": platform.system(),
        "machine": platform.machine(),
        "platform": sys.platform,
        "byteorder": sys.byteorder,
    }


def _verify_binding_semantics(binding: ImmutableBinding, raw: bytes) -> None:
    value = _json_object(raw, binding.name)
    if value.get("artifact_id") != binding.artifact_id:
        raise TDG8SuccessorAuthorityError(
            f"{binding.name} immutable artifact identity differs"
        )
    if binding.name == "auth1":
        payload = value.get("artifact_payload")
        numerical = (
            payload.get("numerical_environment")
            if isinstance(payload, Mapping)
            else None
        )
        contract = numerical.get("contract") if isinstance(numerical, Mapping) else None
        if not isinstance(contract, Mapping):
            raise TDG8SuccessorAuthorityError("AUTH1 numerical environment is absent")
        observed = {key: contract.get(key) for key in PINNED_ENVIRONMENT}
        if not _exact(observed, PINNED_ENVIRONMENT):
            raise TDG8SuccessorAuthorityError("AUTH1 numerical environment differs")


def _verify_live_contract(
    root: Path,
    config: Mapping[str, Any],
) -> tuple[tuple[ImmutableBinding, ...], tuple[ImplementationBinding, ...]]:
    for binding in IMMUTABLE_BINDINGS:
        raw = _read_repository_leaf(root, binding.path, label=binding.name)
        if _sha(raw) != binding.sha256:
            raise TDG8SuccessorAuthorityError(
                f"{binding.name} immutable binding bytes differ"
            )
        _verify_binding_semantics(binding, raw)

    inventory: list[ImplementationBinding] = []
    for row in config["implementation_inventory"]:
        raw = _read_repository_leaf(
            root,
            row["path"],
            label=f"implementation {row['path']}",
        )
        if _sha(raw) != row["sha256"]:
            raise TDG8SuccessorAuthorityError(
                f"implementation binding differs: {row['path']}"
            )
        inventory.append(
            ImplementationBinding(
                role=row["role"], path=row["path"], sha256=row["sha256"]
            )
        )
    if observed_environment() != PINNED_ENVIRONMENT:
        raise TDG8SuccessorAuthorityError("pinned execution environment differs")
    return IMMUTABLE_BINDINGS, tuple(inventory)


def _result_payload(config: Mapping[str, Any]) -> dict[str, Any]:
    bindings = [
        {
            **row,
            "observed_sha256": row["sha256"],
            "binding_status": "sealed",
        }
        for row in _binding_rows()
    ]
    inventory = [
        {
            **row,
            "observed_sha256": row["sha256"],
            "binding_status": "sealed",
        }
        for row in config["implementation_inventory"]
    ]
    return {
        "predecessor": _predecessor_contract(),
        "immutable_bindings": bindings,
        "implementation_inventory": inventory,
        "environment": {
            "contract": dict(PINNED_ENVIRONMENT),
            "observed": dict(PINNED_ENVIRONMENT),
            "binding_status": "sealed",
        },
        "progression": dict(PROGRESSION_CONTRACT),
        "stores": dict(STORES_CONTRACT),
        "prelaunch_observation": {
            "destination_path": DESTINATION_CAMPAIGN_PATH,
            "destination_initially_absent": True,
            "observation_phase": "before_first_creation",
            "postexecution_reobservation_required": False,
        },
        "scope": dict(SCOPE_CONTRACT),
        "claims": dict(CLAIMS_CONTRACT),
        "nonclaims": list(NONCLAIMS),
    }


def _expected_result(config_raw: bytes, config: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL,
        "classification": CLASSIFICATION,
        "gate_status": "pass",
        "source_config_sha256": _sha(config_raw),
        "artifact_payload": _result_payload(config),
    }


def build_prelaunch(
    config_raw: bytes,
    root: Path,
    observe_destination: bool = True,
) -> dict[str, Any]:
    """Observe and bind the one-time prelaunch boundary.

    ``observe_destination=False`` is rejected rather than manufacturing an
    absence claim.  Call :func:`validate_compact` for every post-write or
    post-execution check.
    """

    if not observe_destination:
        raise TDG8SuccessorAuthorityError(
            "prelaunch cannot be built without the one-time destination observation"
        )
    repository = Path(root)
    config = _parse_config(config_raw)
    _verify_repository_root(repository)
    _require_predecessor_ancestry(repository, "HEAD")
    live_config = _read_repository_leaf(
        repository,
        CONFIG_PATH,
        label="TDG8 successor config",
        maximum_bytes=_MAX_CONFIG_BYTES,
    )
    if live_config != config_raw:
        raise TDG8SuccessorAuthorityError("live TDG8 successor config bytes differ")
    _verify_live_contract(repository, config)
    # Existing compact authority blocks re-observation, and the destination
    # check is deliberately last to keep the observation close to publication.
    _require_path_absent(repository, RESULT_PATH, label="TDG8 successor result")
    _require_path_absent(
        repository,
        DESTINATION_CAMPAIGN_PATH,
        label="TDG8 successor destination",
    )
    return _expected_result(config_raw, config)


def validate_compact(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    """Validate the immutable compact pair without observing either store."""

    config = _parse_config(config_raw)
    result = _parse_result(result_raw)
    expected = _expected_result(config_raw, config)
    if result_raw != canonical_result(expected):
        raise TDG8SuccessorAuthorityError("TDG8 successor compact result differs")
    return result


def _verify_committed_leaf(
    root: Path,
    commit: str,
    relative: str,
    expected_raw: bytes,
) -> None:
    relative = _relative_path(relative, "committed authority path")
    tracked = (
        _git(root, "ls-tree", "-r", "--name-only", commit, "--", relative)
        .stdout.decode("utf-8")
        .splitlines()
    )
    if tracked != [relative]:
        raise TDG8SuccessorAuthorityError(
            f"authority path is not tracked at the authority commit: {relative}"
        )
    try:
        committed_size = int(
            _git(root, "cat-file", "-s", f"{commit}:{relative}").stdout.strip()
        )
    except ValueError as exc:
        raise TDG8SuccessorAuthorityError(
            f"committed authority size is malformed: {relative}"
        ) from exc
    if committed_size != len(expected_raw):
        raise TDG8SuccessorAuthorityError(
            f"committed authority size differs: {relative}"
        )
    committed = _git(root, "show", f"{commit}:{relative}").stdout
    live = _read_repository_leaf(root, relative, label=f"authority path {relative}")
    if committed != expected_raw or live != expected_raw:
        raise TDG8SuccessorAuthorityError(
            f"committed authority bytes differ: {relative}"
        )


def authorize_execution(
    root: Path,
    config_raw: bytes,
    result_raw: bytes,
    authority_commit: str,
) -> TDG8SuccessorExecutionAuthority:
    """Authenticate the exact clean committed image for one fresh event.

    This function is intentionally destination- and source-store blind.  The
    runner owns exclusive fresh-namespace creation after this receipt exists.
    """

    config = _parse_config(config_raw)
    validate_compact(config_raw, result_raw)
    repository = Path(root)
    commit = _commit(authority_commit)
    _verify_repository_root(repository)
    head = _git(repository, "rev-parse", "HEAD").stdout.decode("ascii").strip()
    if head != commit:
        raise TDG8SuccessorAuthorityError(
            "current HEAD differs from the TDG8 successor authority commit"
        )
    resolved = (
        _git(repository, "rev-parse", "--verify", f"{commit}^{{commit}}")
        .stdout.decode("ascii")
        .strip()
    )
    if resolved != commit:
        raise TDG8SuccessorAuthorityError("authority commit does not resolve exactly")
    _require_predecessor_ancestry(repository, commit)
    if _git(repository, "status", "--porcelain=v1", "--untracked-files=no").stdout:
        raise TDG8SuccessorAuthorityError("tracked worktree changes are present")

    _verify_committed_leaf(repository, commit, CONFIG_PATH, config_raw)
    _verify_committed_leaf(repository, commit, RESULT_PATH, result_raw)
    bindings, inventory = _verify_live_contract(repository, config)
    for binding in bindings:
        raw = _read_repository_leaf(repository, binding.path, label=binding.name)
        _verify_committed_leaf(repository, commit, binding.path, raw)
    for binding in inventory:
        raw = _read_repository_leaf(
            repository,
            binding.path,
            label=f"implementation {binding.path}",
        )
        _verify_committed_leaf(repository, commit, binding.path, raw)
    _require_exact_authority_delta(repository, commit)

    # Close the narrow race around the initial clean-tree and HEAD checks.
    if _git(repository, "rev-parse", "HEAD").stdout.decode("ascii").strip() != commit:
        raise TDG8SuccessorAuthorityError("HEAD changed during authority verification")
    if _git(repository, "status", "--porcelain=v1", "--untracked-files=no").stdout:
        raise TDG8SuccessorAuthorityError(
            "tracked worktree changed during authority verification"
        )

    return TDG8SuccessorExecutionAuthority(
        authority_commit=commit,
        predecessor_repair_commit=PREDECESSOR_REPAIR_COMMIT,
        config_sha256=_sha(config_raw),
        result_sha256=_sha(result_raw),
        campaign_id=CAMPAIGN_ID,
        target_protocol=TARGET_PROTOCOL,
        branch="GR-0",
        amplitude="3",
        event=23,
        start_rational="23/16",
        target_rational="3/2",
        event_count=1,
        source_campaign_path=SOURCE_CAMPAIGN_PATH,
        destination_campaign_path=DESTINATION_CAMPAIGN_PATH,
        source_terminal_read_only=True,
        candidate_branches_forbidden=True,
        immutable_bindings=bindings,
        implementation_inventory=inventory,
        environment=MappingProxyType(dict(PINNED_ENVIRONMENT)),
    )


__all__ = [
    "ARTIFACT_ID",
    "AUTHORITY_DELTA_PATHS",
    "CAMPAIGN_ID",
    "CLASSIFICATION",
    "CONFIG_PATH",
    "DESTINATION_CAMPAIGN_PATH",
    "IMMUTABLE_BINDINGS",
    "ImplementationBinding",
    "ImmutableBinding",
    "PINNED_ENVIRONMENT",
    "PREDECESSOR_REPAIR_COMMIT",
    "PROJECT_VERSION",
    "REQUIRED_IMPLEMENTATION_INVENTORY",
    "RESULT_PATH",
    "SCHEMA_VERSION",
    "SOURCE_CAMPAIGN_PATH",
    "TARGET_PROTOCOL",
    "TDG8SuccessorAuthorityError",
    "TDG8SuccessorExecutionAuthority",
    "authorize_execution",
    "build_prelaunch",
    "canonical_result",
    "observed_environment",
    "validate_compact",
]
