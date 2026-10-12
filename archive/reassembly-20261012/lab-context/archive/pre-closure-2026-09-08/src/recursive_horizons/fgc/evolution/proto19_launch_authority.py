"""External Git/image authority for the first HLT16 progression event.

The authority observes only committed source/configuration evidence and the
interpreter environment.  It never opens a run store, imports the numerical
engine, or advances a trajectory.  The external Git commit is deliberately
an argument: this module cannot authenticate a commit which contains its own
new authorization result.
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
from typing import Any, Mapping

from .proto19_progression_contract import ProgressionPlan, construct_first_event


SEALED_PRO19_FREEZE_COMMIT = "6be5a1b88ef5e1448599beb111b96b5ddb8a79c8"
MANIFEST_SCHEMA = "FGC-1-PRO19-launch-authority-v1"
REQUIRED_ROLES = frozenset({"runtime", "config", "result", "documentation", "reproducer", "runner", "test"})
_MON16_ARTIFACT_ID = "FGC-1-HLT16-MON16"
_MON16_CONFIG_PATH = "configs/fgc/fgc-1-hlt16-mon16.toml"
_MON16_RESULT_PATH = "results/fgc-1-hlt16-mon16.json"
_MON16_SEALED_CLASSIFICATION = "sealed_durable_runtime_qualified_nonexecuting"
_MON16_CLAIMS = {
    "durable_runtime_implemented": True,
    "durable_runtime_qualified": True,
    "one_event_preflight_authorized": False,
    "state_advance_authorized": False,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "physical_transition_claim_authorized": False,
}
_MON16_NONEXECUTION = {
    "arrays_opened": False,
    "live_run_store_opened": False,
    "runner_invoked": False,
    "one_event_preflight_completed": False,
}


class Proto19LaunchAuthorityError(RuntimeError):
    """The external launch image differs from the prospectively frozen authority."""


@dataclass(frozen=True, slots=True)
class LaunchAuthorityReceipt:
    authorization_commit: str
    manifest_path: str
    manifest_sha256: str
    freeze_commit: str
    freeze_is_ancestor_of_authorization: bool
    authority_paths: tuple[tuple[str, str, str], ...]
    environment: Mapping[str, str]
    progression_plan: ProgressionPlan


def _sha(value: bytes) -> str:
    return sha256(value).hexdigest()


def _commit(value: object) -> str:
    if not isinstance(value, str) or len(value) != 40:
        raise Proto19LaunchAuthorityError("authorization commit is not a 40-character Git id")
    try:
        int(value, 16)
    except ValueError as error:
        raise Proto19LaunchAuthorityError("authorization commit is not hexadecimal") from error
    if value != value.lower():
        raise Proto19LaunchAuthorityError("authorization commit must be lowercase")
    return value


def _git(root: Path, *args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    result = subprocess.run(("git", "-C", str(root), *args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if check and result.returncode:
        raise Proto19LaunchAuthorityError(f"Git authority command failed: {' '.join(args)}")
    return result


def _relative_path(value: object) -> str:
    if not isinstance(value, str) or not value or "\\" in value:
        raise Proto19LaunchAuthorityError("authority path is not a portable relative path")
    path = Path(value)
    if path.is_absolute() or ".." in path.parts or "." in path.parts:
        raise Proto19LaunchAuthorityError("authority path escapes the repository")
    return path.as_posix()


def _nofollow_regular_bytes(root: Path, relative: str) -> bytes:
    """Read one regular repository file through an FD-relative no-follow walk."""
    parts = Path(relative).parts
    if not parts:
        raise Proto19LaunchAuthorityError("authority path is empty")
    directory_fd = -1
    leaf_fd = -1
    try:
        directory_fd = os.open(
            root,
            os.O_RDONLY
            | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_NOFOLLOW", 0),
        )
        if not stat.S_ISDIR(os.fstat(directory_fd).st_mode):
            raise Proto19LaunchAuthorityError("repository root is not a directory")
        for component in parts[:-1]:
            child_fd = os.open(
                component,
                os.O_RDONLY
                | getattr(os, "O_DIRECTORY", 0)
                | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=directory_fd,
            )
            if not stat.S_ISDIR(os.fstat(child_fd).st_mode):
                os.close(child_fd)
                raise Proto19LaunchAuthorityError(
                    "authority path contains a non-directory component"
                )
            os.close(directory_fd)
            directory_fd = child_fd
        leaf_fd = os.open(
            parts[-1],
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=directory_fd,
        )
        leaf = os.fstat(leaf_fd)
        if not stat.S_ISREG(leaf.st_mode) or leaf.st_nlink != 1:
            raise Proto19LaunchAuthorityError(
                "authority path is not a unique regular file"
            )
        chunks: list[bytes] = []
        while True:
            block = os.read(leaf_fd, 1024 * 1024)
            if not block:
                break
            chunks.append(block)
        return b"".join(chunks)
    except FileNotFoundError as error:
        raise Proto19LaunchAuthorityError("authority path is absent") from error
    except OSError as error:
        raise Proto19LaunchAuthorityError("authority path cannot be read safely") from error
    finally:
        if leaf_fd != -1:
            os.close(leaf_fd)
        if directory_fd != -1:
            os.close(directory_fd)


def _canonical_json(raw: bytes, label: str) -> Mapping[str, Any]:
    def duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        output: dict[str, Any] = {}
        for key, value in pairs:
            if key in output:
                raise ValueError(key)
            output[key] = value
        return output
    try:
        value = json.loads(raw.decode("utf-8"), object_pairs_hook=duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise Proto19LaunchAuthorityError(f"{label} is not valid JSON") from error
    if not isinstance(value, Mapping):
        raise Proto19LaunchAuthorityError(f"{label} is not object-valued")
    return value


def observed_environment() -> dict[str, str]:
    """The standard-library-observable numerical execution identity."""
    try:
        numpy_version = version("numpy")
    except PackageNotFoundError as error:
        raise Proto19LaunchAuthorityError("NumPy distribution metadata is absent") from error
    return {
        "implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "numpy_version": numpy_version,
        "system": platform.system(),
        "machine": platform.machine(),
        "platform": sys.platform,
        "byteorder": sys.byteorder,
    }


def _environment_from_auth1(root: Path, path: str) -> dict[str, str]:
    payload = _canonical_json(_nofollow_regular_bytes(root, path), "AUTH1 result")
    artifact = payload.get("artifact_payload")
    if not isinstance(artifact, Mapping):
        raise Proto19LaunchAuthorityError("AUTH1 result omits its artifact payload")
    numerical = artifact.get("numerical_environment")
    if not isinstance(numerical, Mapping) or not isinstance(numerical.get("contract"), Mapping):
        raise Proto19LaunchAuthorityError("AUTH1 result omits its numerical environment")
    contract = numerical["contract"]
    source = {
        "implementation": contract.get("implementation"),
        "python_version": contract.get("python_version"),
        "numpy_version": contract.get("numpy_version"),
        "system": contract.get("system"),
        "machine": contract.get("machine"),
        "platform": contract.get("platform"),
        "byteorder": contract.get("byteorder"),
    }
    if any(not isinstance(value, str) or not value for value in source.values()):
        raise Proto19LaunchAuthorityError("AUTH1 environment has a non-text required field")
    return dict(source)  # type: ignore[arg-type]


def _load_manifest(root: Path, manifest_path: str, commit: str) -> tuple[dict[str, Any], bytes]:
    relative = _relative_path(manifest_path)
    tracked = _git(root, "ls-tree", "-r", "--name-only", commit, "--", relative).stdout.decode("utf-8").splitlines()
    if tracked != [relative]:
        raise Proto19LaunchAuthorityError("launch manifest is not tracked at the authorization commit")
    committed = _git(root, "show", f"{commit}:{relative}").stdout
    live = _nofollow_regular_bytes(root, relative)
    if live != committed:
        raise Proto19LaunchAuthorityError("launch manifest live bytes differ from authorization commit")
    try:
        manifest = tomllib.loads(live.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise Proto19LaunchAuthorityError("launch manifest is malformed") from error
    if not isinstance(manifest, dict) or manifest.get("schema") != MANIFEST_SCHEMA:
        raise Proto19LaunchAuthorityError("launch manifest schema differs")
    return manifest, live


def _manifest_paths(manifest: Mapping[str, Any]) -> tuple[tuple[str, str, str], ...]:
    entries = manifest.get("authority_path")
    if not isinstance(entries, list) or not entries:
        raise Proto19LaunchAuthorityError("launch manifest omits authority paths")
    result: list[tuple[str, str, str]] = []
    seen: set[str] = set()
    roles: set[str] = set()
    for item in entries:
        if not isinstance(item, Mapping) or set(item) != {"role", "path", "sha256"}:
            raise Proto19LaunchAuthorityError("authority path entry fields differ")
        role = item["role"]
        path = _relative_path(item["path"])
        digest = item["sha256"]
        if not isinstance(role, str) or role not in REQUIRED_ROLES or path in seen:
            raise Proto19LaunchAuthorityError("authority path role or uniqueness differs")
        if not isinstance(digest, str) or len(digest) != 64 or digest.lower() != digest:
            raise Proto19LaunchAuthorityError("authority path digest differs")
        try:
            int(digest, 16)
        except ValueError as error:
            raise Proto19LaunchAuthorityError("authority path digest differs") from error
        seen.add(path); roles.add(role); result.append((role, path, digest))
    if roles != REQUIRED_ROLES:
        raise Proto19LaunchAuthorityError("launch manifest does not close every required authority role")
    return tuple(sorted(result, key=lambda item: item[1]))


def _verify_paths(
    root: Path,
    commit: str,
    entries: tuple[tuple[str, str, str], ...],
) -> dict[str, bytes]:
    """Capture the exact live bytes proved equal to the committed image."""
    tracked = set(_git(root, "ls-tree", "-r", "--name-only", commit).stdout.decode("utf-8").splitlines())
    captured: dict[str, bytes] = {}
    for _role, path, expected in entries:
        if path not in tracked:
            raise Proto19LaunchAuthorityError("authority path is not tracked at authorization commit")
        committed = _git(root, "show", f"{commit}:{path}").stdout
        live = _nofollow_regular_bytes(root, path)
        if live != committed or _sha(committed) != expected:
            raise Proto19LaunchAuthorityError("authority path bytes or digest differ")
        captured[path] = live
    return captured


def _verify_environment(
    root: Path,
    manifest: Mapping[str, Any],
    entries: tuple[tuple[str, str, str], ...],
) -> dict[str, str]:
    environment = manifest.get("environment")
    if not isinstance(environment, Mapping) or set(environment) not in ({"source", "auth1_result_path"}, {"source", "contract"}):
        raise Proto19LaunchAuthorityError("launch environment contract fields differ")
    source = environment.get("source")
    if source == "AUTH1":
        auth1_path = _relative_path(environment["auth1_result_path"])
        if auth1_path not in {path for _role, path, _digest in entries}:
            raise Proto19LaunchAuthorityError("AUTH1 environment result is outside the closed authority manifest")
        expected = _environment_from_auth1(root, auth1_path)
    elif source == "HLT16":
        contract = environment.get("contract")
        if not isinstance(contract, Mapping) or set(contract) != set(observed_environment()):
            raise Proto19LaunchAuthorityError("HLT16 environment contract differs")
        expected = dict(contract)
        if any(not isinstance(value, str) or not value for value in expected.values()):
            raise Proto19LaunchAuthorityError("HLT16 environment contract is malformed")
    else:
        raise Proto19LaunchAuthorityError("launch environment source differs")
    observed = observed_environment()
    if observed != expected:
        raise Proto19LaunchAuthorityError("current Python/NumPy/platform environment differs from the launch contract")
    return observed


def _mon16_config_bindings(raw: bytes) -> tuple[tuple[str, str, str], ...]:
    """Parse only MON16's declared predecessor/inventory binding surface.

    This stays deliberately independent from the MON16 artifact builder: the
    launch gate needs a second parser with no ability to refresh or construct
    a certificate.
    """
    try:
        config = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise Proto19LaunchAuthorityError("MON16 config is malformed") from exc
    if not isinstance(config, Mapping) or config.get("artifact_id") != _MON16_ARTIFACT_ID:
        raise Proto19LaunchAuthorityError("MON16 config identity differs")
    bindings: list[tuple[str, str, str]] = []
    seen: set[str] = set()
    for section, required_keys in (
        ("predecessors", {"artifact_id", "path", "sha256"}),
        ("inventory", {"role", "path", "sha256"}),
    ):
        rows = config.get(section)
        if not isinstance(rows, list) or not rows:
            raise Proto19LaunchAuthorityError(f"MON16 config {section} differs")
        seen_in_section: set[str] = set()
        for row in rows:
            if not isinstance(row, Mapping) or set(row) != required_keys:
                raise Proto19LaunchAuthorityError(f"MON16 config {section} entry differs")
            path = _relative_path(row.get("path"))
            digest = row.get("sha256")
            role = "predecessor" if section == "predecessors" else row.get("role")
            if (not isinstance(role, str) or not isinstance(digest, str)
                    or len(digest) != 64 or digest.lower() != digest or path in seen_in_section):
                raise Proto19LaunchAuthorityError(f"MON16 config {section} binding differs")
            try:
                int(digest, 16)
            except ValueError as exc:
                raise Proto19LaunchAuthorityError(f"MON16 config {section} binding differs") from exc
            seen_in_section.add(path)
            bindings.append((str(role), path, digest))
    return tuple(sorted(bindings, key=lambda item: item[1]))


def _verify_sealed_mon16(
    config_raw: bytes,
    result_raw: bytes,
    captured: Mapping[str, bytes],
) -> None:
    """Validate the MON16 result as a sealed, non-executing predecessor.

    This intentionally parses the compact result independently of the MON16
    builder.  A launch image cannot convert a pending inventory or a result
    with expanded claims into authority merely because its own file hash was
    listed in the manifest.
    """
    bindings = _mon16_config_bindings(config_raw)
    expected_paths = {_MON16_CONFIG_PATH, _MON16_RESULT_PATH, *(path for _role, path, _sha256 in bindings)}
    if not expected_paths.issubset(captured):
        raise Proto19LaunchAuthorityError("launch manifest omits MON16 configured bindings")
    config_sha256 = _sha(config_raw)
    result = _canonical_json(result_raw, "HLT16 MON16 result")
    if result.get("artifact_id") != _MON16_ARTIFACT_ID:
        raise Proto19LaunchAuthorityError("MON16 result identity differs")
    if result.get("classification") != _MON16_SEALED_CLASSIFICATION:
        raise Proto19LaunchAuthorityError("MON16 result is not sealed")
    payload = result.get("artifact_payload")
    if not isinstance(payload, Mapping):
        raise Proto19LaunchAuthorityError("MON16 result omits artifact payload")
    if payload.get("claims") != _MON16_CLAIMS:
        raise Proto19LaunchAuthorityError("MON16 durable-runtime claim boundary differs")
    if payload.get("explicit_nonexecution") != _MON16_NONEXECUTION:
        raise Proto19LaunchAuthorityError("MON16 nonexecution boundary differs")
    pending = payload.get("coordinator_refresh_required_paths")
    if not isinstance(pending, list) or pending:
        raise Proto19LaunchAuthorityError("MON16 result retains pending coordinator bindings")
    if result.get("source_config_sha256") != config_sha256:
        raise Proto19LaunchAuthorityError("MON16 result/source-config binding differs")

    expected_predecessors = {
        path: digest for role, path, digest in bindings if role == "predecessor"
    }
    expected_inventory = {
        path: (role, digest) for role, path, digest in bindings if role != "predecessor"
    }
    predecessors = payload.get("predecessors")
    inventory = payload.get("inventory")
    if not isinstance(predecessors, list) or not isinstance(inventory, list):
        raise Proto19LaunchAuthorityError("MON16 sealed binding records differ")

    actual_predecessors: dict[str, Mapping[str, Any]] = {}
    actual_inventory: dict[str, Mapping[str, Any]] = {}
    for label, rows, target in (
        ("predecessor", predecessors, actual_predecessors),
        ("inventory", inventory, actual_inventory),
    ):
        for item in rows:
            if not isinstance(item, Mapping):
                raise Proto19LaunchAuthorityError(f"MON16 {label} record differs")
            path = _relative_path(item.get("path"))
            if path in target or item.get("binding_status") != "sealed":
                raise Proto19LaunchAuthorityError(f"MON16 {label} record differs")
            target[path] = item
    if set(actual_predecessors) != set(expected_predecessors) or set(actual_inventory) != set(expected_inventory):
        raise Proto19LaunchAuthorityError("MON16 result/config binding set differs")

    for path, digest in expected_predecessors.items():
        item = actual_predecessors[path]
        observed = _sha(captured[path])
        if item.get("sha256") != digest or item.get("observed_sha256") != digest or observed != digest:
            raise Proto19LaunchAuthorityError("MON16 predecessor actual-byte binding differs")
    for path, (role, digest) in expected_inventory.items():
        item = actual_inventory[path]
        observed = _sha(captured[path])
        if (item.get("role") != role or item.get("sha256") != digest
                or item.get("observed_sha256") != digest or observed != digest):
            raise Proto19LaunchAuthorityError("MON16 inventory actual-byte binding differs")


def authorize_first_event(
    root: Path,
    *,
    authorization_commit: str,
    manifest_path: str,
) -> LaunchAuthorityReceipt:
    """Authenticate the external image, then construct the static first-event plan."""
    repository = Path(root)
    commit = _commit(authorization_commit)
    head = _git(repository, "rev-parse", "HEAD").stdout.decode("ascii").strip()
    if head != commit:
        raise Proto19LaunchAuthorityError("current HEAD differs from external authorization commit")
    resolved = _git(repository, "rev-parse", "--verify", f"{commit}^{{commit}}")
    if resolved.stdout.decode("ascii").strip() != commit:
        raise Proto19LaunchAuthorityError("authorization commit does not resolve exactly")
    if _git(repository, "status", "--porcelain", "--untracked-files=no").stdout:
        raise Proto19LaunchAuthorityError("tracked worktree changes are present")
    if _git(
        repository,
        "merge-base",
        "--is-ancestor",
        SEALED_PRO19_FREEZE_COMMIT,
        commit,
        check=False,
    ).returncode != 0:
        raise Proto19LaunchAuthorityError(
            "sealed PRO19 freeze is not an ancestor of the authorization commit"
        )
    manifest, bytes_ = _load_manifest(repository, manifest_path, commit)
    entries = _manifest_paths(manifest)
    captured = _verify_paths(repository, commit, entries)
    environment = _verify_environment(repository, manifest, entries)
    progression_evidence = (
        "results/fgc-1-pro19-frz1.json",
        "results/fgc-1-pro18-auth1.json",
    )
    required_evidence = {
        *progression_evidence,
        _MON16_RESULT_PATH,
        _MON16_CONFIG_PATH,
    }
    if not required_evidence.issubset(captured):
        raise Proto19LaunchAuthorityError(
            "launch manifest omits compact progression evidence"
        )
    _verify_sealed_mon16(
        captured[_MON16_CONFIG_PATH],
        captured[_MON16_RESULT_PATH],
        captured,
    )
    plan = construct_first_event(
        repository,
        authorization_commit=commit,
        # The pure PRO19 plan owns its frozen two-result contract.  MON16's
        # result and config are authenticated here for the independent
        # cross-bind; neither is smuggled into that older constructor API.
        evidence_bytes={path: captured[path] for path in progression_evidence},
    )
    return LaunchAuthorityReceipt(
        authorization_commit=commit,
        manifest_path=_relative_path(manifest_path),
        manifest_sha256=_sha(bytes_),
        freeze_commit=SEALED_PRO19_FREEZE_COMMIT,
        freeze_is_ancestor_of_authorization=True,
        authority_paths=entries,
        environment=environment,
        progression_plan=plan,
    )


__all__ = [
    "LaunchAuthorityReceipt", "MANIFEST_SCHEMA", "Proto19LaunchAuthorityError",
    "REQUIRED_ROLES", "SEALED_PRO19_FREEZE_COMMIT", "authorize_first_event", "observed_environment",
]
