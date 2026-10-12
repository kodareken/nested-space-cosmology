"""Prospective committed-image authority for TDG11-MSEL1.

A valid compact freeze is not execution authority. Live authorization inspects
the committed successor image, environment, process table, output namespace,
and the sealed store. Compact verification reads only tracked config, result,
implementation, and pinned-reference bytes.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass, is_dataclass
from hashlib import sha1, sha256
import json
import math
import os
from pathlib import Path
import re
import stat
from typing import Mapping, NoReturn

from recursive_horizons.evidence_io import (
    DEFAULT_MAX_BYTES,
    CommitParents,
    EvidenceIOError,
    InspectDelta,
    InspectTree,
    InspectWorktree,
    ResolveCommit,
    TreeEntry,
    canonical_json_bytes,
    git_read,
    load_canonical_json,
    read_regular_file,
)

from . import tdg11_msel1_contract as contract
from .tdg11_msel1_contract import MSEL1ContractError


_COMMIT = re.compile(r"[0-9a-f]{40}\Z")
_STAGING_PREFIX = ".tdg11-msel1.stage-"
_SNAPSHOT_DOMAIN = b"TDG9-LOC1-STORE-SNAPSHOT-v1\n"
_PS_MAX_BYTES = DEFAULT_MAX_BYTES
_FORBIDDEN_PROCESS = (
    "scripts/run_fgc_",
    "run_fgc_test_partition",
    "-m unittest",
    "unittest.main",
    "-m pytest",
    "pytest",
    "codex review",
    "review-agent",
)
_REGULAR_GIT_MODES = frozenset({"100644", "100755"})


class MSEL1AuthorityError(RuntimeError):
    """The prospective TDG11-MSEL1 image is not authorized."""


@dataclass(frozen=True, slots=True)
class MSEL1Authority:
    authority_commit: str
    config_sha256: str
    result_sha256: str
    environment: dict[str, str]
    store_snapshot: tuple[int, str]
    config: dict[str, object]


@dataclass(frozen=True, slots=True)
class PredecessorInspection:
    retry: int
    generation: int
    checkpoint_sha256: str
    journal_tip_sha256: str
    descriptor_sha256: str
    pending_owner: str
    accepted_time_hex: str
    prior_retry_count: int
    checkpoint_raw_sha256: str
    journal_tip_raw_sha256: str
    rejection_raw_sha256: str


@dataclass(frozen=True, slots=True)
class RestoredPredecessor:
    retry: int
    member: object
    checkpoint_sha256: str
    descriptor_sha256: str
    fingerprint: dict[str, object]


def _fail(message: str) -> NoReturn:
    raise MSEL1AuthorityError(message)


def _absolute_root(root: Path) -> Path:
    if not isinstance(root, Path):
        root = Path(root)
    supplied = Path(os.fspath(root))
    if not supplied.is_absolute():
        _fail("root must be an absolute canonical path")
    try:
        resolved = supplied.resolve(strict=True)
    except OSError as error:
        raise MSEL1AuthorityError("root cannot be resolved") from error
    if resolved != supplied:
        _fail("root path traverses a symlink")
    return supplied


def _wrap(error: Exception, message: str) -> NoReturn:
    raise MSEL1AuthorityError(message) from error


def _git_blob_id(payload: bytes) -> str:
    return sha1(
        b"blob " + str(len(payload)).encode("ascii") + b"\0" + payload,
        usedforsecurity=False,
    ).hexdigest()


def _worktree_git_mode(metadata: os.stat_result) -> str:
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
        _fail("worktree path is not a regular file")
    return "100755" if metadata.st_mode & 0o111 else "100644"


def _open_nofollow_regular(path: str) -> int:
    if not hasattr(os, "O_NOFOLLOW"):
        _fail("no-follow reads are unavailable")
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
    try:
        return os.open(path, flags)
    except OSError as error:
        raise MSEL1AuthorityError("regular file cannot be opened") from error


def _read_fd_regular(
    descriptor: int,
    *,
    before: os.stat_result,
    allow_hardlink: bool,
    max_bytes: int,
) -> bytes:
    if (
        stat.S_ISLNK(before.st_mode)
        or not stat.S_ISREG(before.st_mode)
        or before.st_size < 0
        or before.st_size > max_bytes
        or (not allow_hardlink and before.st_nlink != 1)
    ):
        _fail("regular file is unsafe")
    active = os.fstat(descriptor)
    identity = _file_identity(before)
    if identity != _file_identity(active) or not stat.S_ISREG(active.st_mode):
        _fail("regular file raced")
    chunks: list[bytes] = []
    remaining = int(before.st_size)
    while remaining > 0:
        chunk = os.read(descriptor, min(remaining, 1 << 20))
        if not chunk:
            break
        chunks.append(chunk)
        remaining -= len(chunk)
    raw = b"".join(chunks)
    after = os.fstat(descriptor)
    if (
        len(raw) != before.st_size
        or _file_identity(after) != identity
        or not stat.S_ISREG(after.st_mode)
        or (not allow_hardlink and after.st_nlink != 1)
    ):
        _fail("regular file changed during read")
    return raw


def _file_identity(metadata: os.stat_result) -> tuple[int, ...]:
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_mode,
        metadata.st_nlink,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
    )


def _hash_software_file(path: str) -> str:
    if type(path) is not str or not path or "\0" in path or len(path) > 1024:
        _fail("software path is unsafe")
    try:
        before = os.lstat(path)
    except OSError as error:
        raise MSEL1AuthorityError("software path cannot be read") from error
    descriptor = _open_nofollow_regular(path)
    try:
        raw = _read_fd_regular(
            descriptor,
            before=before,
            allow_hardlink=True,
            max_bytes=DEFAULT_MAX_BYTES,
        )
        if _file_identity(os.lstat(path)) != _file_identity(before):
            _fail("software path changed during read")
        return sha256(raw).hexdigest()
    finally:
        os.close(descriptor)


def _runtime_regular_path(path: str) -> str:
    if type(path) is not str or not path or "\0" in path:
        _fail("runtime path is unsafe")
    try:
        resolved = os.path.realpath(path)
    except OSError as error:
        raise MSEL1AuthorityError("runtime path cannot be resolved") from error
    if type(resolved) is not str or not resolved or "\0" in resolved:
        _fail("runtime path is unsafe")
    if len(resolved) > 1024:
        _fail("runtime path exceeds the environment bound")
    try:
        metadata = os.lstat(resolved)
    except OSError as error:
        raise MSEL1AuthorityError("runtime path cannot be read") from error
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
        _fail("runtime path is not a regular file")
    return resolved


def _numpy_extension_path() -> str:
    import importlib

    for name in (
        "numpy._core._multiarray_umath",
        "numpy.core._multiarray_umath",
    ):
        try:
            module = importlib.import_module(name)
        except ImportError:
            continue
        path = getattr(module, "__file__", None)
        if type(path) is str and path:
            return _runtime_regular_path(path)
    _fail("NumPy compiled extension is unavailable")


def observe_environment() -> dict[str, str]:
    """Observe the live arithmetic image. Compact verification never calls this."""

    import platform
    import sys

    import numpy as np

    configuration = np.show_config(mode="dicts")
    if type(configuration) is not dict:
        _fail("NumPy build configuration is unavailable")
    dependencies = configuration.get("Build Dependencies", {})
    blas = dependencies.get("blas", {}) if type(dependencies) is dict else {}
    if type(blas) is not dict:
        _fail("NumPy BLAS configuration is unavailable")
    python_executable = _runtime_regular_path(sys.executable)
    numpy_extension = _numpy_extension_path()
    observed = {
        "python_implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "numpy_version": str(np.__version__),
        "blas_name": str(blas.get("name") or ""),
        "blas_version": str(blas.get("version") or "unknown"),
        "system": platform.system(),
        "machine": platform.machine(),
        "byteorder": sys.byteorder,
        "kernel_release": platform.release(),
        "kernel_version": platform.version(),
        "macos_version": platform.mac_ver()[0],
        "python_executable": python_executable,
        "python_executable_sha256": _hash_software_file(python_executable),
        "numpy_extension": numpy_extension,
        "numpy_extension_sha256": _hash_software_file(numpy_extension),
    }
    try:
        return contract._environment_record(observed)
    except MSEL1ContractError as error:
        _wrap(error, "environment image differs")


def _store_leaf_bytes(root: Path, relative: str) -> bytes:
    try:
        return read_regular_file(root, relative)
    except EvidenceIOError as error:
        _wrap(error, "store leaf cannot be read safely")


def _hash_store_tree(base: Path) -> tuple[int, str]:
    base = _absolute_root(base)
    try:
        metadata = base.lstat()
    except OSError as error:
        raise MSEL1AuthorityError("store root cannot be read") from error
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        _fail("store root is unsafe")
    digest = sha256(_SNAPSHOT_DOMAIN)
    count = 0
    stack = [base]
    directory_identities = {}
    while stack:
        directory = stack.pop()
        flags = (
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0)
        )
        try:
            directory_fd = os.open(directory, flags)
        except OSError as error:
            raise MSEL1AuthorityError("store directory cannot be opened") from error
        try:
            directory_identities[directory] = _file_identity(os.fstat(directory_fd))
            with os.scandir(directory_fd) as entries:
                items = sorted(
                    (
                        (
                            item.name,
                            item.is_symlink(),
                            item.stat(follow_symlinks=False),
                        )
                        for item in entries
                    ),
                    key=lambda row: row[0],
                )
            prefix = "" if directory == base else directory.relative_to(base).as_posix()
            for name, is_link, metadata in items:
                if is_link:
                    _fail("store contains a symlink")
                relative = f"{prefix}/{name}" if prefix else name
                if stat.S_ISDIR(metadata.st_mode):
                    stack.append(directory / name)
                elif stat.S_ISREG(metadata.st_mode):
                    raw = _store_leaf_bytes(base, relative)
                    digest.update(
                        (relative + "\0" + sha256(raw).hexdigest() + "\n").encode(
                            "ascii"
                        )
                    )
                    count += 1
                else:
                    _fail("store contains a special file")
        finally:
            os.close(directory_fd)
    for directory, identity in directory_identities.items():
        _absolute_root(directory)
        if _file_identity(directory.lstat()) != identity:
            _fail("store directory changed during snapshot")
    return count, digest.hexdigest()


def snapshot_store(root: Path) -> tuple[int, str]:
    """Hash the sealed store with the historical TDG9-LOC1 snapshot domain."""

    root = _absolute_root(root)
    pair = _hash_store_tree(root / contract.STORE_PATH)
    if pair != (contract.STORE_LEAF_COUNT, contract.STORE_SHA256):
        _fail("store snapshot differs")
    return pair


def require_output_absent(root: Path) -> None:
    """Require the output namespace parent chain and a missing final leaf."""

    root = _absolute_root(root)
    relative = Path(contract.OUTPUT_NAMESPACE)
    if relative.is_absolute() or any(
        part in {"", ".", ".."} for part in relative.parts
    ):
        _fail("output namespace is unsafe")
    current = root
    for part in relative.parent.parts:
        current = current / part
        try:
            metadata = current.lstat()
        except FileNotFoundError:
            _fail("output parent is absent")
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
            _fail("output parent is unsafe")
    target = root / relative
    if os.path.lexists(target):
        _fail("output namespace already exists")
    try:
        names = os.listdir(current)
    except OSError as error:
        raise MSEL1AuthorityError("output parent cannot be listed") from error
    if any(name.startswith(_STAGING_PREFIX) for name in names):
        _fail("staging namespace already exists")


def _implementation_records(root: Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    for path in contract.IMPLEMENTATION_PATHS:
        try:
            digest = sha256(read_regular_file(root, path)).hexdigest()
        except EvidenceIOError as error:
            _wrap(error, f"implementation path cannot be read: {path}")
        records.append({"path": path, "sha256": digest})
    return records


def _require_pinned_references(root: Path) -> None:
    for path, digest in contract.PINNED_REFERENCES.items():
        try:
            observed = sha256(read_regular_file(root, path)).hexdigest()
        except EvidenceIOError as error:
            _wrap(error, f"pinned reference cannot be read: {path}")
        if observed != digest:
            _fail(f"pinned reference differs: {path}")


def _tracked_payload(root: Path) -> tuple[bytes, bytes, dict[str, object]]:
    try:
        config_raw = read_regular_file(root, contract.CONFIG_PATH)
        result_raw = read_regular_file(root, contract.RESULT_PATH)
        result = load_canonical_json(result_raw.removesuffix(b"\n"))
        expected = contract.validate_compact(config_raw, result)
    except (MSEL1ContractError, EvidenceIOError) as error:
        _wrap(error, "tracked compact bundle differs")
    config = expected["config"]
    if type(config) is not dict:
        _fail("compact config is missing")
    declared = config.get("implementation")
    if type(declared) is not list:
        _fail("implementation inventory differs")
    for item, path in zip(declared, contract.IMPLEMENTATION_PATHS, strict=True):
        if type(item) is not dict:
            _fail("implementation inventory differs")
        try:
            observed = sha256(read_regular_file(root, item["path"])).hexdigest()
        except EvidenceIOError as error:
            _wrap(error, f"implementation path cannot be read: {path}")
        if item["path"] != path or observed != item["sha256"]:
            _fail(f"implementation binding differs: {path}")
    _require_pinned_references(root)
    return config_raw, result_raw, expected


def validate_compact_bundle(root: Path) -> dict[str, object]:
    """Git/store/shadow/environment-blind compact freeze verification."""

    root = _absolute_root(root)
    _config_raw, _result_raw, expected = _tracked_payload(root)
    return expected


def emit_config_bytes(root: Path) -> bytes:
    """Deterministic config from the live environment and current hashes."""

    root = _absolute_root(root)
    environment = observe_environment()
    implementation = _implementation_records(root)
    _require_pinned_references(root)
    try:
        return contract.render_config(
            contract.expected_config(
                environment=environment, implementation=implementation
            )
        )
    except MSEL1ContractError as error:
        _wrap(error, "emitted config differs")


def emit_result_object(root: Path) -> dict[str, object]:
    """Compact freeze of the existing tracked config; no file writes."""

    root = _absolute_root(root)
    try:
        config_raw = read_regular_file(root, contract.CONFIG_PATH)
        return contract.compact_freeze(config_raw)
    except (MSEL1ContractError, EvidenceIOError) as error:
        _wrap(error, "emitted compact result differs")


def _require_authority_commit(authority_commit: str) -> str:
    if type(authority_commit) is not str or not _COMMIT.fullmatch(authority_commit):
        _fail("authority commit is not a 40-hex SHA-1 identity")
    return authority_commit


def _git_query(root: Path, operation: object) -> object:
    try:
        return git_read(root, operation)
    except EvidenceIOError as error:
        _wrap(error, "committed image cannot be inspected")


def _require_committed_image(root: Path, authority_commit: str) -> None:
    head = _git_query(root, ResolveCommit("HEAD"))
    if head != authority_commit:
        _fail("authority commit is not HEAD")
    parents = _git_query(root, CommitParents(authority_commit))
    if parents != (contract.BASE_COMMIT,):
        _fail("authority commit is not the exact direct successor")
    delta = _git_query(
        root, InspectDelta(commit=authority_commit, against=contract.BASE_COMMIT)
    )
    paths = tuple(sorted(item.path for item in delta))
    if paths != contract.DELTA_PATHS:
        _fail("committed delta differs")
    worktree = _git_query(root, InspectWorktree())
    if not worktree.clean:
        _fail("authority working image is not clean")
    tree = _git_query(root, InspectTree(commit=authority_commit))
    ceiling = contract.RESOURCES["git_listing_entry_ceiling"]
    if type(ceiling) is not int or len(tree) < 1 or len(tree) > ceiling:
        _fail("committed tree listing differs")
    seen: set[str] = set()
    for entry in tree:
        if not isinstance(entry, TreeEntry):
            _fail("committed tree entry differs")
        if (
            entry.kind != "blob"
            or entry.mode not in _REGULAR_GIT_MODES
            or len(entry.object_id) != 40
            or not _COMMIT.fullmatch(entry.object_id)
        ):
            _fail("committed tree contains a nonregular path")
        if entry.path in seen:
            _fail("committed tree repeats a path")
        seen.add(entry.path)
        try:
            payload = read_regular_file(root, entry.path)
            metadata = (root / entry.path).lstat()
        except (EvidenceIOError, OSError) as error:
            _wrap(error, f"committed path cannot be read: {entry.path}")
        if _worktree_git_mode(metadata) != entry.mode:
            _fail(f"committed mode differs: {entry.path}")
        if _git_blob_id(payload) != entry.object_id:
            _fail(f"committed bytes differ: {entry.path}")


def _process_rows() -> tuple[tuple[int, int, str], ...]:
    import subprocess
    from recursive_horizons.evidence_io import _communicate_bounded

    try:
        proc = subprocess.Popen(
            ["/bin/ps", "-axo", "pid=,ppid=,command="],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            close_fds=True,
            start_new_session=True,
        )
    except OSError as error:
        raise MSEL1AuthorityError("process table cannot be observed") from error
    assert proc.stdout is not None
    try:
        raw = _communicate_bounded(proc, max_bytes=_PS_MAX_BYTES, timeout=5)
        code = proc.returncode
    except Exception as error:
        proc.kill()
        proc.wait()
        raise MSEL1AuthorityError("process table cannot be observed") from error
    if len(raw) > _PS_MAX_BYTES or code != 0:
        _fail("process table is unbounded or unavailable")
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise MSEL1AuthorityError("process table is not text") from error
    rows: list[tuple[int, int, str]] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split(None, 2)
        if len(parts) < 2:
            _fail("process table is malformed")
        try:
            pid = int(parts[0])
            ppid = int(parts[1])
        except ValueError:
            _fail("process table is malformed")
        command = parts[2] if len(parts) == 3 else ""
        rows.append((pid, ppid, command))
    return tuple(rows)


def _protected_pids(rows: tuple[tuple[int, int, str], ...]) -> set[int]:
    parents = {pid: ppid for pid, ppid, _command in rows}
    protected = {os.getpid(), os.getppid()}
    current = os.getpid()
    seen: set[int] = set()
    while current and current not in seen:
        seen.add(current)
        protected.add(current)
        current = parents.get(current, 0)
    return protected


def _require_process_preflight() -> None:
    rows = _process_rows()
    protected = _protected_pids(rows)
    for pid, _ppid, command in rows:
        if pid in protected:
            continue
        if any(marker in command for marker in _FORBIDDEN_PROCESS):
            _fail("another scientific runner, test partition, or review is active")


def _replay_spec(retry: int) -> dict[str, object]:
    if type(retry) is not int:
        _fail("retry must be an integer selector")
    for spec in contract.replay_specs():
        if spec["retry"] == retry:
            return dict(spec)
    _fail(f"retry is not a selected predecessor: {retry}")


def _checkpoint_relative(generation: int, digest: str) -> str:
    return f"{contract.STORE_PATH}/checkpoints/{generation:020d}-{digest}.json"


def _journal_relative(sequence: int, digest: str) -> str:
    return f"{contract.STORE_PATH}/journal/{sequence:020d}-{digest}.journal"


def _mapping_field(value: object, key: str) -> object:
    if isinstance(value, Mapping):
        return value.get(key)
    return getattr(value, key, None)


def _accepted_time_hex(cursor: object) -> str:
    payload = cursor
    nested = _mapping_field(payload, "payload")
    if isinstance(nested, Mapping) and "accepted_boundary_time" in nested:
        payload = nested
    accepted = _mapping_field(payload, "accepted_boundary_time")
    if isinstance(accepted, Mapping):
        hex_value = accepted.get("binary64_hex")
        if type(hex_value) is str:
            return hex_value
    _fail("predecessor accepted time differs")


def _pending_owner(state: object) -> str:
    owner = _mapping_field(state, "pending_owner")
    if type(owner) is str:
        return owner
    _fail("predecessor pending owner differs")


def _descriptor_sha256(state: object) -> str:
    digest = _mapping_field(state, "descriptor_sha256")
    if type(digest) is str:
        return digest
    _fail("predecessor descriptor differs")


def _prior_retry_count(state: object) -> int:
    ledger = _mapping_field(state, "ledger")
    count = _mapping_field(ledger, "current_macro_step_temporal_retry_count")
    if type(count) is int and not isinstance(count, bool):
        return count
    _fail("predecessor temporal retry count differs")


def _json_object(raw: bytes) -> dict[str, object]:
    try:
        value = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise MSEL1AuthorityError("predecessor leaf is not JSON") from error
    if type(value) is not dict:
        _fail("predecessor leaf is not an object")
    return value


def _inspect_replay(
    root: Path,
    store: object,
    spec: Mapping[str, object],
) -> PredecessorInspection:
    generation = spec["generation"]
    retry = spec["retry"]
    if type(generation) is not int or type(retry) is not int:
        _fail("replay identity differs")
    try:
        checkpoint = store.authenticated_checkpoint_at_generation(generation)
    except Exception as error:
        raise MSEL1AuthorityError(
            f"generation-{generation} checkpoint cannot be authenticated"
        ) from error
    if (
        checkpoint.sha256 != spec["checkpoint_sha256"]
        or checkpoint.generation != generation
        or checkpoint.journal_sequence != spec["journal_tip_sequence"]
        or checkpoint.journal_tip_sha256 != spec["journal_tip_sha256"]
    ):
        _fail(f"generation-{generation} checkpoint identity differs")
    members = checkpoint.members
    if not isinstance(members, Mapping) or contract.MEMBER_KEY not in members:
        _fail("predecessor member is absent")
    state = members[contract.MEMBER_KEY]
    descriptor = _descriptor_sha256(state)
    owner = _pending_owner(state)
    accepted = _accepted_time_hex(_mapping_field(state, "cursor"))
    prior = _prior_retry_count(state)
    cursor_mode = _mapping_field(_mapping_field(state, "cursor"), "mode")
    if (
        descriptor != contract.DESCRIPTOR_SHA256
        or owner != "temporal"
        or accepted != contract.ACCEPTED_TIME_HEX
        or prior != spec["prior_retry_count"]
        or cursor_mode != "RETRY_PENDING"
    ):
        _fail(f"retry-{retry} member identity differs")
    checkpoint_relative = _checkpoint_relative(
        generation, str(spec["checkpoint_sha256"])
    )
    journal_relative = _journal_relative(
        int(spec["journal_tip_sequence"]), str(spec["journal_tip_sha256"])
    )
    rejection_relative = _journal_relative(
        int(spec["rejection_sequence"]), str(spec["rejection_sha256"])
    )
    checkpoint_raw = _store_leaf_bytes(root, checkpoint_relative)
    journal_raw = _store_leaf_bytes(root, journal_relative)
    rejection_raw = _store_leaf_bytes(root, rejection_relative)
    if (
        sha256(checkpoint_raw).hexdigest() != spec["checkpoint_raw_sha256"]
        or sha256(journal_raw).hexdigest() != spec["journal_tip_raw_sha256"]
        or sha256(rejection_raw).hexdigest() != spec["rejection_raw_sha256"]
    ):
        _fail(f"retry-{retry} predecessor leaf hash differs")
    checkpoint_json = _json_object(checkpoint_raw)
    journal_json = _json_object(journal_raw)
    rejection_json = _json_object(rejection_raw)
    if (
        checkpoint_json.get("generation") != generation
        or checkpoint_json.get("checkpoint_sha256") != spec["checkpoint_sha256"]
        or checkpoint_json.get("journal_sequence") != spec["journal_tip_sequence"]
        or checkpoint_json.get("journal_tip_sha256") != spec["journal_tip_sha256"]
    ):
        _fail(f"retry-{retry} checkpoint leaf identity differs")
    if (
        journal_json.get("sequence") != spec["journal_tip_sequence"]
        or journal_json.get("record_sha256") != spec["journal_tip_sha256"]
        or journal_json.get("kind") == "tdg6_rejection"
    ):
        _fail(f"retry-{retry} journal tip is not the accepted generation tip")
    if (
        rejection_json.get("sequence") != spec["rejection_sequence"]
        or rejection_json.get("record_sha256") != spec["rejection_sha256"]
        or rejection_json.get("kind") != "tdg6_rejection"
    ):
        _fail(f"retry-{retry} historical rejection leaf differs")
    return PredecessorInspection(
        retry=retry,
        generation=generation,
        checkpoint_sha256=str(spec["checkpoint_sha256"]),
        journal_tip_sha256=str(spec["journal_tip_sha256"]),
        descriptor_sha256=descriptor,
        pending_owner=owner,
        accepted_time_hex=accepted,
        prior_retry_count=prior,
        checkpoint_raw_sha256=str(spec["checkpoint_raw_sha256"]),
        journal_tip_raw_sha256=str(spec["journal_tip_raw_sha256"]),
        rejection_raw_sha256=str(spec["rejection_raw_sha256"]),
    )


def inspect_predecessors(root: Path) -> tuple[PredecessorInspection, ...]:
    """Authenticate generations 9/10/11 read-only. Rejection leaves stay unrestored."""

    from .hlt16_campaign_store import HLT16CampaignStore, HLT16CampaignStoreError

    root = _absolute_root(root)
    try:
        store = HLT16CampaignStore(root / contract.STORE_PATH)
    except HLT16CampaignStoreError as error:
        _wrap(error, "campaign store cannot be opened read-only")
    inspections = []
    for spec in contract.replay_specs():
        inspections.append(_inspect_replay(root, store, spec))
    if tuple(item.generation for item in inspections) != (9, 10, 11):
        _fail("predecessor generation set differs")
    return tuple(inspections)


def _jsonable(value: object) -> object:
    if value is None or type(value) in (str, int, bool):
        return value
    if type(value) is float:
        if not math.isfinite(value):
            _fail("fingerprint contains a nonfinite binary64")
        return value.hex()
    if type(value) is dict:
        return {str(key): _jsonable(item) for key, item in value.items()}
    if type(value) is list or type(value) is tuple:
        return [_jsonable(item) for item in value]
    if is_dataclass(value) and not isinstance(value, type):
        return _jsonable(asdict(value))
    _fail(f"fingerprint cannot encode {type(value).__name__}")


def _tracer_history_sha256(member: object) -> str:
    from .numerical_engine import array_content_sha256

    tracers = member.tracers
    arrays = [
        tracers.labels,
        tracers.positions,
        tracers.proper_times,
        *tuple(tracers.event_proper_times),
        *tuple(tracers.event_fields),
    ]
    return array_content_sha256(*arrays)


def _member_fingerprint(member: object, descriptor_sha256: str) -> dict[str, object]:
    from .numerical_engine import array_content_sha256

    ledger = member.temporal_ledger
    if ledger is None:
        _fail("restored temporal ledger is absent")
    transaction = member.transaction
    payload = {
        "monitor": asdict(transaction.state),
        "causal": asdict(transaction.causal_state),
        "ledger": asdict(ledger),
        "step_index": member.step_index,
        "transaction_serial": member.transaction_serial,
        "source_retry_count": member.source_retry_count,
        "CFL_retry_count": member.CFL_retry_count,
    }
    return {
        "member_key": member.key,
        "method_label": member.method_label,
        "integrator_id": member.integrator_id,
        "spatial_order": member.spatial_order,
        "point_count": member.point_count,
        "accepted_time_hex": float(member.time).hex(),
        "state_sha256": array_content_sha256(
            member.state.u, member.state.p, member.state.q
        ),
        "coordinates_sha256": array_content_sha256(member.initial.grid.coordinates),
        "descriptor_sha256": descriptor_sha256,
        "tracer_history_sha256": _tracer_history_sha256(member),
        "transaction_sha256": sha256(
            canonical_json_bytes(_jsonable(payload))
        ).hexdigest(),
        "source_retry_count": member.source_retry_count,
        "CFL_retry_count": member.CFL_retry_count,
        "current_macro_step_temporal_retry_count": (
            ledger.current_macro_step_temporal_retry_count
        ),
    }


def restore_predecessor(
    root: Path,
    retry: int,
    *,
    templates: Mapping[str, object],
) -> RestoredPredecessor:
    """Restore one RK4-2049 predecessor onto a caller-built template."""

    from .hlt16_campaign_runtime import restore_member_with_overlay
    from .hlt16_campaign_store import HLT16CampaignStore, HLT16CampaignStoreError
    from .numerical_engine import PRIMARY_METHOD, array_content_sha256

    root = _absolute_root(root)
    if not isinstance(templates, Mapping) or len(templates) != 6:
        _fail("static GR-0 templates were not supplied")
    if contract.MEMBER_KEY not in templates:
        _fail("RK4-2049 template is absent")
    spec = _replay_spec(retry)
    generation = spec["generation"]
    if type(generation) is not int:
        _fail("replay generation differs")
    try:
        store = HLT16CampaignStore(root / contract.STORE_PATH)
        checkpoint = store.authenticated_checkpoint_at_generation(generation)
    except HLT16CampaignStoreError as error:
        _wrap(error, "predecessor checkpoint cannot be authenticated")
    if checkpoint.sha256 != spec["checkpoint_sha256"]:
        _fail(f"retry-{retry} checkpoint identity differs")
    members = checkpoint.members
    if not isinstance(members, Mapping) or contract.MEMBER_KEY not in members:
        _fail("predecessor member is absent")
    state = members[contract.MEMBER_KEY]
    descriptor = _descriptor_sha256(state)
    member = deepcopy(templates[contract.MEMBER_KEY])
    try:
        restored = restore_member_with_overlay(
            store, checkpoint, member, key=contract.MEMBER_KEY
        )
    except Exception as error:
        raise MSEL1AuthorityError(
            f"retry-{retry} predecessor cannot be restored"
        ) from error
    if restored is not None:
        member = restored
    ledger = member.temporal_ledger
    try:
        state_digest = array_content_sha256(
            member.state.u, member.state.p, member.state.q
        )
        coordinates_digest = array_content_sha256(member.initial.grid.coordinates)
    except Exception as error:
        raise MSEL1AuthorityError("restored arrays cannot be hashed") from error
    if (
        member.key != contract.MEMBER_KEY
        or member.point_count != contract.POINT_COUNT
        or member.integrator_id != PRIMARY_METHOD
        or member.spatial_order != 4
        or float(member.time).hex() != contract.ACCEPTED_TIME_HEX
        or ledger is None
        or ledger.current_macro_step_temporal_retry_count != spec["prior_retry_count"]
        or state_digest != contract.PHYSICAL_STATE_SHA256
        or coordinates_digest != contract.COORDINATES_SHA256
        or descriptor != contract.DESCRIPTOR_SHA256
    ):
        _fail(f"retry-{retry} restored identity differs")
    fingerprint = _member_fingerprint(member, descriptor)
    if (
        fingerprint["state_sha256"] != contract.PHYSICAL_STATE_SHA256
        or fingerprint["coordinates_sha256"] != contract.COORDINATES_SHA256
        or fingerprint["accepted_time_hex"] != contract.ACCEPTED_TIME_HEX
        or fingerprint["descriptor_sha256"] != contract.DESCRIPTOR_SHA256
    ):
        _fail(f"retry-{retry} predecessor fingerprint differs")
    return RestoredPredecessor(
        retry=retry,
        member=member,
        checkpoint_sha256=str(spec["checkpoint_sha256"]),
        descriptor_sha256=descriptor,
        fingerprint=fingerprint,
    )


def authorize_execution(root: Path, *, authority_commit: str) -> MSEL1Authority:
    """Authorize one prospective diagnostic against the committed successor image."""

    root = _absolute_root(root)
    commit = _require_authority_commit(authority_commit)
    _require_committed_image(root, commit)
    config_raw, result_raw, compact = _tracked_payload(root)
    config = compact["config"]
    if type(config) is not dict:
        _fail("authorized config is missing")
    environment = observe_environment()
    if environment != config["environment"]:
        _fail("observed environment differs from the config snapshot")
    _require_process_preflight()
    require_output_absent(root)
    store = snapshot_store(root)
    inspect_predecessors(root)
    return MSEL1Authority(
        authority_commit=commit,
        config_sha256=sha256(config_raw).hexdigest(),
        result_sha256=sha256(result_raw).hexdigest(),
        environment=dict(environment),
        store_snapshot=store,
        config=deepcopy(config),
    )


__all__ = (
    "MSEL1Authority",
    "MSEL1AuthorityError",
    "PredecessorInspection",
    "RestoredPredecessor",
    "authorize_execution",
    "emit_config_bytes",
    "emit_result_object",
    "inspect_predecessors",
    "observe_environment",
    "require_output_absent",
    "restore_predecessor",
    "snapshot_store",
    "validate_compact_bundle",
)
