"""Prospective freeze for the TDG8 RCV3 generation-nine projection.

RCV3 does not edit or resume the terminal source campaign.  This module binds
the complete terminal tree, selects the exact accepted prefix through
generation nine / journal sequence ten, and records that a new projection
namespace was absent before any bootstrap.  The selected prefix deliberately
omits terminal history and every lock: guards and quarantined writer leases are
coordination state, not append-only campaign lineage.

The live builder is a one-time, read-only observation.  Compact verification
uses only the tracked config/result bytes and never re-observes destination
absence after a later bootstrap.
"""

from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path, PurePosixPath
import stat
import subprocess
import tomllib
from typing import Any, Iterable, Mapping, Sequence


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG8-RCV3-FRZ1"
PROJECT_VERSION = "0.11.0"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CLASSIFICATION = "premise_only_generation9_prefix_projection_freeze"
CONFIG_PATH = "configs/fgc/fgc-1-tdg8-rcv3-frz1.toml"
RESULT_PATH = "results/fgc-1-tdg8-rcv3-frz1.json"
OWNER_DOCUMENT = "docs/fgc-tdg8-rcv3.md"
BOOTSTRAP_RUNTIME_PATH = (
    "src/recursive_horizons/fgc/evolution/tdg8_rcv3_fork_runtime.py"
)
BOOTSTRAP_RUNTIME_SHA256 = (
    "9c7e40dd76d126fa4ce708025483d0c15a65e9f36a80c5919335bd53c73dca36"
)
BOOTSTRAP_SCRIPT_PATH = "scripts/bootstrap_fgc_tdg8_rcv3.py"
BOOTSTRAP_SCRIPT_SHA256 = (
    "0e3bbfc4f23508c8a5a03764c54c56d60a3e6e9c7f2eeec53bc9e1a9f78b614a"
)

SOURCE_STORE_PATH = "runs/fgc-2-sf1/proto19/calibration"
DESTINATION_WRAPPER_PATH = "runs/fgc-2-sf1/tdg8-rcv3"
DESTINATION_STORE_PATH = "runs/fgc-2-sf1/tdg8-rcv3/calibration"
PROJECTION_ID = "FGC-1-TDG8-RCV3-PROJECTION-1"

REPAIR_COMMIT = "cde091d018d4196f8ba212ce59d606b0fb8db992"
REPAIR_SOURCES = (
    (
        "src/recursive_horizons/fgc/evolution/hlt16_progression_attempt.py",
        "d43f51ce92beda8843ea6d8bbef19f20af230d996611e3890cda43292f6e46fa",
        "b9fae332e3e48d004f6edbf1ee6f61fad45c8bfc",
    ),
    (
        "src/recursive_horizons/fgc/evolution/tdg8_persisted_retry_replay.py",
        "13b9de4dfee683afb405cbe17b007ffd8679b36bb946c43f0fdf469ee2dce9ee",
        "7661293ac583f4434d8e1535735bb2a842ccc5db",
    ),
)

SOURCE_CAMPAIGN_ID = "FGC-2-SF1-PROTO18-GR0-A3-TDG8-SUCCESSOR-1"
SOURCE_AUTHORIZATION_COMMIT = "6df67967e6e0bc63eca4dc6951a60d3787368f26"
SOURCE_PLAN_SHA256 = "e649e477b98411a3051f9946f2a9248b875d2b598818b8983ce06caa969d2126"

GENERATION9_CHECKPOINT_SHA256 = (
    "eb6fddc480c94aa7ed15c80399fc2b26637693efef02399f267075fc6c258e56"
)
GENERATION9_RAW_SHA256 = (
    "07572d4450829ff2e2bad65e217a917b0d31ba677a3b8fc20a3be43e849ba846"
)
SEQUENCE10_JOURNAL_SHA256 = (
    "5b533eb7009a9c7c3f353d26813f9cbbde8574afae77492929b435b4c701c5ff"
)
SEQUENCE10_RAW_SHA256 = (
    "b50dc39ebd8bbdd3729d40d9e9ae4b223ff2ef19389d4ca74f9a316bac72202f"
)
GENERATION10_TERMINAL_SHA256 = (
    "291eb466d5cb566b1f03822bf109c0c9e0f6ff930e9be496d26309b53b5005ea"
)
SEQUENCE11_TERMINAL_SHA256 = (
    "394c803df15f0b73827425dec4c914a1dd998aad8de23d8a8ecf3d637e1cf340"
)
TERMINAL_LOCK_SHA256 = (
    "3bfb54c0e23d0940d52418402cdd24e8444a843398f3bbfb194d32eac05e8606"
)

MEMBER_KEY = "RK4-2049"
MEMBER_DESCRIPTOR_SHA256 = (
    "77847126340e78c8ac300fac2795bceda724c15e0894e34444b0da42a84166b4"
)
MEMBER_PHYSICAL_STATE_SHA256 = (
    "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a"
)
MEMBER_ACCEPTED_TIME_HEX = "0x1.78554de5a30e0p+0"
MEMBER_CURSOR_SHA256 = (
    "8f17e3efd76df310c97229672363609881c32cd063cd2386d6086387181a95c4"
)
RETRY_DEPTH = 2
PENDING_CAP_HEX = "0x1.aaa9612df8000p-11"

EXPECTED_DIRECTORIES = (
    "checkpoints",
    "journal",
    "locks",
    "payloads",
    "receipts",
    "states",
)
FULL_MANIFEST_COUNT = 59
FULL_MANIFEST_BYTES = 6_342_040
FULL_MANIFEST_SHA256 = (
    "b076b5233e5e2ee9898444f6cade495e2e0ed8d24697348f640485104d2cf145"
)
PREFIX_MANIFEST_COUNT = 50
PREFIX_MANIFEST_BYTES = 6_201_938
PREFIX_MANIFEST_SHA256 = (
    "d84e85c3dd29ac1554111ed00b37289d2b89273e7a7c2c47e3d21f5b9fda9e0b"
)

GENERATION9_PATH = (
    "checkpoints/00000000000000000009-"
    f"{GENERATION9_CHECKPOINT_SHA256}.json"
)
SEQUENCE10_PATH = (
    "journal/00000000000000000010-"
    f"{SEQUENCE10_JOURNAL_SHA256}.journal"
)
GENERATION10_PATH = (
    "checkpoints/00000000000000000010-"
    f"{GENERATION10_TERMINAL_SHA256}.json"
)
SEQUENCE11_PATH = (
    "journal/00000000000000000011-"
    f"{SEQUENCE11_TERMINAL_SHA256}.journal"
)
EXCLUDED_PATHS = (
    GENERATION10_PATH,
    SEQUENCE11_PATH,
    "locks/.active-write.lock.hlt16-quarantine-1d1ebcfbcb57ecf91adb191815acd57c1137609cd50f35fa98e7d18d2c91086c",
    "locks/.active-write.lock.hlt16-quarantine-2f76eeac85746761d38617804da02f0d520676a86a243a5e29b2a9f0fdfa71e9",
    "locks/.active-write.lock.hlt16-quarantine-6c409d46d6c049cced6a26d9ad733911070dad9d53e966674a7fdce4ba372946",
    "locks/.active-write.lock.hlt16-quarantine-f474863d36b1c0312ff677560f1385c2913add7353bde694f08e69161bcef5ec",
    "locks/bootstrap.guard",
    "locks/terminal.lock",
    "locks/writer.guard",
)

NONCLAIMS = (
    "FRZ1 freezes a byte-for-byte recovery projection input; it does not create, authenticate, or execute that projection.",
    "The projected prefix deliberately preserves the source campaign identity internally; PROJECTION_ID and the fixed destination path provide its external fork identity.",
    "Generation ten, sequence eleven, the terminal lock, four quarantined writer leases, and two coordination guards remain bound as source-terminal evidence but are excluded from the projected lineage.",
    "FRZ1 binds one projection-only runtime and bootstrap script; no evolution or candidate entry point is authorized by those bytes.",
    "FRZ1 records no common-event, calibration, trapping, activation, DEF1, retained-EFT, transition, mechanism, or physical result.",
)

SCOPE = {
    "source_store_read_only": True,
    "source_terminal_tree_bound": True,
    "projection_prefix_selected": True,
    "destination_absence_one_time": True,
    "bootstrap_projection_authorized": True,
    "bootstrap_completed": False,
    "destination_store_created": False,
    "bounded_execution_authorized": False,
    "source_campaign_resume_authorized": False,
    "threshold_or_retry_rule_change": False,
    "candidate_branches_forbidden": True,
}

CLAIMS = {
    "source_terminal_bound": True,
    "last_valid_generation9_bound": True,
    "prefix_manifest_frozen": True,
    "destination_absence_observed": True,
    "bootstrap_authorized": True,
    "bootstrap_completed": False,
    "projection_authenticated": False,
    "bounded_execution_authorized": False,
    "common_event_completed": False,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "physical_result_earned": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
}


class TDG8RCV3FreezeError(RuntimeError):
    """The RCV3 freeze config, source, namespace, or compact result differs."""


def canonical(value: object) -> bytes:
    """Canonical compact bytes used for hashes and embedded store objects."""
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("ascii")
    except (TypeError, ValueError) as exc:
        raise TDG8RCV3FreezeError("value is not canonicalizable") from exc


def canonical_result(value: object) -> bytes:
    """Canonical tracked result representation."""
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
        raise TDG8RCV3FreezeError("result is not canonicalizable") from exc


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _require_sha(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or value.lower() != value:
        raise TDG8RCV3FreezeError(f"{label} is not SHA-256")
    try:
        int(value, 16)
    except ValueError as exc:
        raise TDG8RCV3FreezeError(f"{label} is not SHA-256") from exc
    return value


def _safe_relative(relative: str) -> tuple[str, ...]:
    path = PurePosixPath(relative)
    if (
        not relative
        or path.is_absolute()
        or str(path) != relative
        or any(part in ("", ".", "..") for part in path.parts)
    ):
        raise TDG8RCV3FreezeError("unsafe relative path")
    return path.parts


def _read_fd_all(descriptor: int) -> bytes:
    chunks: list[bytes] = []
    while True:
        chunk = os.read(descriptor, 1024 * 1024)
        if not chunk:
            return b"".join(chunks)
        chunks.append(chunk)


def _open_directory(parent_fd: int, name: str, label: str) -> int:
    try:
        descriptor = os.open(
            name,
            os.O_RDONLY
            | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent_fd,
        )
    except OSError as exc:
        raise TDG8RCV3FreezeError(f"cannot safely open {label}") from exc
    if not stat.S_ISDIR(os.fstat(descriptor).st_mode):
        os.close(descriptor)
        raise TDG8RCV3FreezeError(f"{label} is not a directory")
    return descriptor


def read_nofollow(root: Path, relative: str, label: str) -> bytes:
    """Read one regular, single-link leaf without following any component."""
    parts = _safe_relative(relative)
    descriptors: list[int] = []
    try:
        root_fd = os.open(
            Path(root),
            os.O_RDONLY
            | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_NOFOLLOW", 0),
        )
        descriptors.append(root_fd)
        directory_fd = root_fd
        for part in parts[:-1]:
            directory_fd = _open_directory(directory_fd, part, label)
            descriptors.append(directory_fd)
        leaf_fd = os.open(
            parts[-1],
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=directory_fd,
        )
        descriptors.append(leaf_fd)
        metadata = os.fstat(leaf_fd)
        if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
            raise TDG8RCV3FreezeError(f"{label} is not a safe regular leaf")
        raw = _read_fd_all(leaf_fd)
        if len(raw) != metadata.st_size:
            raise TDG8RCV3FreezeError(f"{label} changed while read")
        return raw
    except TDG8RCV3FreezeError:
        raise
    except OSError as exc:
        raise TDG8RCV3FreezeError(f"cannot safely read {label}") from exc
    finally:
        for descriptor in reversed(descriptors):
            os.close(descriptor)


def _walk_tree(root: Path) -> tuple[tuple[str, ...], tuple[dict[str, object], ...]]:
    """Hash a two-level campaign tree using descriptor-relative no-follow I/O."""
    root_fd = -1
    directories: list[str] = []
    rows: list[dict[str, object]] = []
    try:
        root_fd = os.open(
            Path(root),
            os.O_RDONLY
            | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_NOFOLLOW", 0),
        )
        if not stat.S_ISDIR(os.fstat(root_fd).st_mode):
            raise TDG8RCV3FreezeError("source store is not a directory")
        with os.scandir(root_fd) as iterator:
            top_names = sorted(entry.name for entry in iterator)
        if tuple(top_names) != EXPECTED_DIRECTORIES:
            raise TDG8RCV3FreezeError("source store directory set differs")
        for directory_name in top_names:
            directory_fd = _open_directory(root_fd, directory_name, directory_name)
            try:
                directories.append(directory_name)
                with os.scandir(directory_fd) as iterator:
                    names = sorted(entry.name for entry in iterator)
                for name in names:
                    metadata = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
                    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
                        raise TDG8RCV3FreezeError("source tree contains an unsafe leaf")
                    if metadata.st_nlink != 1:
                        raise TDG8RCV3FreezeError("source tree contains a multiply linked leaf")
                    leaf_fd = os.open(
                        name,
                        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
                        dir_fd=directory_fd,
                    )
                    try:
                        opened = os.fstat(leaf_fd)
                        if (
                            not stat.S_ISREG(opened.st_mode)
                            or opened.st_nlink != 1
                            or opened.st_dev != metadata.st_dev
                            or opened.st_ino != metadata.st_ino
                        ):
                            raise TDG8RCV3FreezeError("source leaf changed during open")
                        raw = _read_fd_all(leaf_fd)
                        if len(raw) != opened.st_size:
                            raise TDG8RCV3FreezeError("source leaf changed during read")
                    finally:
                        os.close(leaf_fd)
                    rows.append({
                        "path": f"{directory_name}/{name}",
                        "byte_count": len(raw),
                        "sha256": _sha(raw),
                    })
            finally:
                os.close(directory_fd)
    except TDG8RCV3FreezeError:
        raise
    except OSError as exc:
        raise TDG8RCV3FreezeError("cannot safely inventory source tree") from exc
    finally:
        if root_fd != -1:
            os.close(root_fd)
    rows.sort(key=lambda row: str(row["path"]))
    return tuple(directories), tuple(rows)


def _json(raw: bytes, label: str) -> dict[str, Any]:
    def reject_duplicates(items: Iterable[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, item in items:
            if key in value:
                raise ValueError(key)
            value[key] = item
        return value
    try:
        value = json.loads(raw.decode("utf-8"), object_pairs_hook=reject_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TDG8RCV3FreezeError(f"{label} is malformed") from exc
    if not isinstance(value, dict):
        raise TDG8RCV3FreezeError(f"{label} is not object-valued")
    return value


def _verify_object_hash(value: Mapping[str, Any], key: str, expected: str, label: str) -> None:
    if value.get(key) != expected:
        raise TDG8RCV3FreezeError(f"{label} identity differs")
    body = dict(value)
    body.pop(key, None)
    if _sha(canonical(body)) != expected:
        raise TDG8RCV3FreezeError(f"{label} content hash differs")


def _manifest_summary(rows: Sequence[Mapping[str, object]]) -> dict[str, object]:
    normalized = [dict(row) for row in rows]
    return {
        "leaf_count": len(normalized),
        "byte_count": sum(int(row["byte_count"]) for row in normalized),
        "manifest_sha256": _sha(canonical(normalized)),
        "leaves": normalized,
    }


def _validate_manifest_rows(rows: object, label: str) -> list[dict[str, object]]:
    if not isinstance(rows, list):
        raise TDG8RCV3FreezeError(f"{label} leaves are not a list")
    normalized: list[dict[str, object]] = []
    previous = ""
    for row in rows:
        if not isinstance(row, Mapping) or set(row) != {"path", "byte_count", "sha256"}:
            raise TDG8RCV3FreezeError(f"{label} leaf schema differs")
        path = row.get("path")
        byte_count = row.get("byte_count")
        if not isinstance(path, str):
            raise TDG8RCV3FreezeError(f"{label} leaf path differs")
        _safe_relative(path)
        if path <= previous:
            raise TDG8RCV3FreezeError(f"{label} leaf ordering differs")
        if not isinstance(byte_count, int) or isinstance(byte_count, bool) or byte_count < 0:
            raise TDG8RCV3FreezeError(f"{label} byte count differs")
        digest = _require_sha(row.get("sha256"), f"{label} leaf hash")
        normalized.append({"path": path, "byte_count": byte_count, "sha256": digest})
        previous = path
    return normalized


def _validate_summary(
    value: object,
    *,
    label: str,
    count: int,
    byte_count: int,
    digest: str,
) -> list[dict[str, object]]:
    if not isinstance(value, Mapping) or set(value) != {
        "leaf_count", "byte_count", "manifest_sha256", "leaves"
    }:
        raise TDG8RCV3FreezeError(f"{label} manifest schema differs")
    rows = _validate_manifest_rows(value.get("leaves"), label)
    observed = _manifest_summary(rows)
    if (
        value.get("leaf_count") != count
        or value.get("byte_count") != byte_count
        or value.get("manifest_sha256") != digest
        or observed["leaf_count"] != count
        or observed["byte_count"] != byte_count
        or observed["manifest_sha256"] != digest
    ):
        raise TDG8RCV3FreezeError(f"{label} manifest identity differs")
    return rows


def expected_source_anchor() -> dict[str, object]:
    return {
        "campaign_id": SOURCE_CAMPAIGN_ID,
        "authorization_commit": SOURCE_AUTHORIZATION_COMMIT,
        "plan_sha256": SOURCE_PLAN_SHA256,
        "protocol": TARGET_PROTOCOL,
        "branch": "GR-0",
        "amplitude": "3",
        "event": 23,
        "target_rational": "3/2",
        "target_binary64_hex": "0x1.8000000000000p+0",
        "generation9_checkpoint_sha256": GENERATION9_CHECKPOINT_SHA256,
        "generation9_raw_sha256": GENERATION9_RAW_SHA256,
        "sequence10_journal_sha256": SEQUENCE10_JOURNAL_SHA256,
        "sequence10_raw_sha256": SEQUENCE10_RAW_SHA256,
        "generation10_terminal_sha256": GENERATION10_TERMINAL_SHA256,
        "sequence11_terminal_sha256": SEQUENCE11_TERMINAL_SHA256,
        "terminal_lock_sha256": TERMINAL_LOCK_SHA256,
        "member_key": MEMBER_KEY,
        "member_descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
        "member_physical_state_sha256": MEMBER_PHYSICAL_STATE_SHA256,
        "member_accepted_time_hex": MEMBER_ACCEPTED_TIME_HEX,
        "member_cursor_sha256": MEMBER_CURSOR_SHA256,
        "member_mode": "RETRY_PENDING",
        "pending_owner": "temporal",
        "retry_depth": RETRY_DEPTH,
        "pending_cap_binary64_hex": PENDING_CAP_HEX,
        "generation9_disposition": "nonterminal",
        "generation10_disposition": "invalid_terminal",
    }


def _expected_repair_sources() -> list[dict[str, str]]:
    return [
        {"path": path, "sha256": digest, "git_blob": blob}
        for path, digest, blob in REPAIR_SOURCES
    ]


def expected_bootstrap_contract() -> dict[str, object]:
    return {
        "runtime_path": BOOTSTRAP_RUNTIME_PATH,
        "runtime_sha256": BOOTSTRAP_RUNTIME_SHA256,
        "script_path": BOOTSTRAP_SCRIPT_PATH,
        "script_sha256": BOOTSTRAP_SCRIPT_SHA256,
        "production_install_function": "tdg8_rcv3_fork_runtime.install_recovery_fork",
        "production_spec": "tdg8_rcv3_fork_runtime.PRODUCTION_SPEC",
        "candidate_or_evolution_entrypoints": False,
    }


def _expected_exclusion() -> dict[str, object]:
    return {
        "excluded_paths": list(EXCLUDED_PATHS),
        "excluded_leaf_count": len(EXCLUDED_PATHS),
        "projected_lock_leaf_count": 0,
        "destination_locks_start_empty": True,
        "terminal_lineage_excluded": True,
        "retired_or_coordination_locks_excluded": True,
    }


def expected_destination_observation() -> dict[str, object]:
    return {
        "projection_id": PROJECTION_ID,
        "wrapper_path": DESTINATION_WRAPPER_PATH,
        "store_path": DESTINATION_STORE_PATH,
        "wrapper_absent": True,
        "store_absent": True,
        "observed_before_bootstrap": True,
        "compact_reobservation_required": False,
    }


def _config_expected() -> dict[str, object]:
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL,
        "classification": CLASSIFICATION,
        "nonclaims": list(NONCLAIMS),
        "repair": {
            "commit": REPAIR_COMMIT,
            "sources": _expected_repair_sources(),
        },
        "bootstrap": expected_bootstrap_contract(),
        "source": {
            "store_path": SOURCE_STORE_PATH,
            "expected_directories": list(EXPECTED_DIRECTORIES),
            "full_manifest_leaf_count": FULL_MANIFEST_COUNT,
            "full_manifest_byte_count": FULL_MANIFEST_BYTES,
            "full_manifest_sha256": FULL_MANIFEST_SHA256,
            "anchor": expected_source_anchor(),
        },
        "projection": {
            "projection_id": PROJECTION_ID,
            "destination_wrapper_path": DESTINATION_WRAPPER_PATH,
            "destination_store_path": DESTINATION_STORE_PATH,
            "prefix_manifest_leaf_count": PREFIX_MANIFEST_COUNT,
            "prefix_manifest_byte_count": PREFIX_MANIFEST_BYTES,
            "prefix_manifest_sha256": PREFIX_MANIFEST_SHA256,
            "exclusion": _expected_exclusion(),
        },
        "scope": dict(SCOPE),
        "claims": dict(CLAIMS),
    }


def parse_config(raw: bytes) -> dict[str, Any]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise TDG8RCV3FreezeError("RCV3 FRZ1 config is malformed") from exc
    if not isinstance(value, dict) or value != _config_expected():
        raise TDG8RCV3FreezeError("RCV3 FRZ1 config identity differs")
    return value


def _git(repository: Path, *arguments: str) -> bytes:
    process = subprocess.run(
        ["git", *arguments],
        cwd=repository,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if process.returncode != 0:
        raise TDG8RCV3FreezeError("cannot authenticate repair commit")
    return process.stdout


def derive_repair_evidence(repository: Path) -> dict[str, object]:
    """Bind the two execution sources repaired at immutable commit cde091d."""
    head = _git(repository, "rev-parse", "HEAD").decode("ascii").strip()
    ancestry = subprocess.run(
        ["git", "merge-base", "--is-ancestor", REPAIR_COMMIT, head],
        cwd=repository,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    if ancestry.returncode != 0:
        raise TDG8RCV3FreezeError("repair commit is not an ancestor of HEAD")
    observed: list[dict[str, str]] = []
    for path, expected_sha, expected_blob in REPAIR_SOURCES:
        committed = _git(repository, "show", f"{REPAIR_COMMIT}:{path}")
        blob = _git(repository, "rev-parse", f"{REPAIR_COMMIT}:{path}").decode("ascii").strip()
        live = read_nofollow(repository, path, f"repair source {path}")
        if _sha(committed) != expected_sha or blob != expected_blob or live != committed:
            raise TDG8RCV3FreezeError("repair source bytes differ")
        observed.append({"path": path, "sha256": expected_sha, "git_blob": expected_blob})
    return {
        "commit": REPAIR_COMMIT,
        "head_at_observation": head,
        "repair_commit_is_ancestor": True,
        "live_sources_equal_commit": True,
        "sources": observed,
    }


def derive_bootstrap_evidence(repository: Path) -> dict[str, object]:
    """Bind the exact projection-only production runtime and bootstrap bytes."""
    expected = expected_bootstrap_contract()
    for kind in ("runtime", "script"):
        path = str(expected[f"{kind}_path"])
        digest = str(expected[f"{kind}_sha256"])
        raw = read_nofollow(repository, path, f"RCV3 {kind}")
        if _sha(raw) != digest:
            raise TDG8RCV3FreezeError(f"RCV3 {kind} bytes differ")
    return expected


def _derive_anchor(store: Path, rows: Sequence[Mapping[str, object]]) -> dict[str, object]:
    row_map = {str(row["path"]): row for row in rows}
    if row_map.get(GENERATION9_PATH, {}).get("sha256") != GENERATION9_RAW_SHA256:
        raise TDG8RCV3FreezeError("generation-nine raw bytes differ")
    if row_map.get(SEQUENCE10_PATH, {}).get("sha256") != SEQUENCE10_RAW_SHA256:
        raise TDG8RCV3FreezeError("sequence-ten raw bytes differ")
    generation9 = _json(
        read_nofollow(store, GENERATION9_PATH, "generation-nine checkpoint"),
        "generation-nine checkpoint",
    )
    sequence10 = _json(
        read_nofollow(store, SEQUENCE10_PATH, "sequence-ten journal"),
        "sequence-ten journal",
    )
    generation10 = _json(
        read_nofollow(store, GENERATION10_PATH, "generation-ten checkpoint"),
        "generation-ten checkpoint",
    )
    sequence11 = _json(
        read_nofollow(store, SEQUENCE11_PATH, "sequence-eleven journal"),
        "sequence-eleven journal",
    )
    terminal_lock = _json(
        read_nofollow(store, "locks/terminal.lock", "terminal lock"),
        "terminal lock",
    )
    _verify_object_hash(
        generation9,
        "checkpoint_sha256",
        GENERATION9_CHECKPOINT_SHA256,
        "generation-nine checkpoint",
    )
    _verify_object_hash(
        sequence10,
        "record_sha256",
        SEQUENCE10_JOURNAL_SHA256,
        "sequence-ten journal",
    )
    _verify_object_hash(
        generation10,
        "checkpoint_sha256",
        GENERATION10_TERMINAL_SHA256,
        "generation-ten checkpoint",
    )
    _verify_object_hash(
        sequence11,
        "record_sha256",
        SEQUENCE11_TERMINAL_SHA256,
        "sequence-eleven journal",
    )
    _verify_object_hash(
        terminal_lock,
        "lock_sha256",
        TERMINAL_LOCK_SHA256,
        "terminal lock",
    )
    expected_common = {
        "campaign_id": SOURCE_CAMPAIGN_ID,
        "event": 23,
    }
    if any(
        value.get("campaign_id") != expected_common["campaign_id"]
        or value.get("event") != expected_common["event"]
        for value in (generation9, generation10, sequence10, sequence11)
    ):
        raise TDG8RCV3FreezeError("campaign/event identity differs")
    if (
        generation9.get("generation") != 9
        or generation9.get("journal_sequence") != 10
        or generation9.get("journal_tip_sha256") != SEQUENCE10_JOURNAL_SHA256
        or generation9.get("disposition") != "nonterminal"
        or generation9.get("authorization_commit") != SOURCE_AUTHORIZATION_COMMIT
        or generation9.get("plan_sha256") != SOURCE_PLAN_SHA256
        or generation9.get("protocol") != TARGET_PROTOCOL
        or generation9.get("target") != {
            "binary64_hex": "0x1.8000000000000p+0", "rational": "3/2"
        }
    ):
        raise TDG8RCV3FreezeError("generation-nine boundary differs")
    if (
        sequence10.get("sequence") != 10
        or sequence10.get("generation") != 9
        or sequence10.get("kind") != "cursor_transition"
        or sequence10.get("previous_record_sha256")
        != "f1d005bb01cf4ea473551ac9fdface8b582d9e5f20cdf63169193d9177ddabed"
    ):
        raise TDG8RCV3FreezeError("sequence-ten boundary differs")
    if (
        generation10.get("generation") != 10
        or generation10.get("parent_sha256") != GENERATION9_CHECKPOINT_SHA256
        or generation10.get("journal_sequence") != 11
        or generation10.get("journal_tip_sha256") != SEQUENCE11_TERMINAL_SHA256
        or generation10.get("disposition") != "invalid_terminal"
        or sequence11.get("sequence") != 11
        or sequence11.get("generation") != 10
        or sequence11.get("kind") != "terminal_lock"
        or sequence11.get("previous_record_sha256") != SEQUENCE10_JOURNAL_SHA256
        or terminal_lock.get("checkpoint_sha256") != GENERATION10_TERMINAL_SHA256
        or terminal_lock.get("journal_tip_sha256") != SEQUENCE11_TERMINAL_SHA256
    ):
        raise TDG8RCV3FreezeError("source terminal boundary differs")
    member = generation9.get("members", {}).get(MEMBER_KEY)
    if not isinstance(member, Mapping):
        raise TDG8RCV3FreezeError("generation-nine member is absent")
    cursor = member.get("cursor")
    successor = cursor.get("retry_successor_payload_or_none") if isinstance(cursor, Mapping) else None
    if (
        member.get("descriptor_sha256") != MEMBER_DESCRIPTOR_SHA256
        or member.get("pending_owner") != "temporal"
        or member.get("pending_cap_hex") != PENDING_CAP_HEX
        or not isinstance(cursor, Mapping)
        or cursor.get("mode") != "RETRY_PENDING"
        or cursor.get("accepted_state_sha256") != MEMBER_DESCRIPTOR_SHA256
        or cursor.get("accepted_boundary_time", {}).get("binary64_hex")
        != MEMBER_ACCEPTED_TIME_HEX
        or cursor.get("cursor_chain_sha256") != MEMBER_CURSOR_SHA256
        or not isinstance(successor, Mapping)
        or successor.get("retry_count") != RETRY_DEPTH
        or successor.get("half_cap") != PENDING_CAP_HEX
        or successor.get("TDG6_rejection_evidence", {}).get("initial_state_sha256")
        != MEMBER_PHYSICAL_STATE_SHA256
    ):
        raise TDG8RCV3FreezeError("generation-nine retry state differs")
    return expected_source_anchor()


def derive_source_evidence(repository: Path) -> dict[str, object]:
    """Independently bind the full source tree and its exact projected prefix."""
    store = Path(repository) / SOURCE_STORE_PATH
    first_directories, first_rows = _walk_tree(store)
    second_directories, second_rows = _walk_tree(store)
    if first_directories != second_directories or first_rows != second_rows:
        raise TDG8RCV3FreezeError("source tree changed during observation")
    full = _manifest_summary(first_rows)
    if (
        first_directories != EXPECTED_DIRECTORIES
        or full["leaf_count"] != FULL_MANIFEST_COUNT
        or full["byte_count"] != FULL_MANIFEST_BYTES
        or full["manifest_sha256"] != FULL_MANIFEST_SHA256
    ):
        raise TDG8RCV3FreezeError("full source manifest differs")
    full_paths = {str(row["path"]) for row in first_rows}
    if not set(EXCLUDED_PATHS).issubset(full_paths):
        raise TDG8RCV3FreezeError("source exclusion set is incomplete")
    prefix_rows = tuple(row for row in first_rows if row["path"] not in EXCLUDED_PATHS)
    prefix = _manifest_summary(prefix_rows)
    if (
        prefix["leaf_count"] != PREFIX_MANIFEST_COUNT
        or prefix["byte_count"] != PREFIX_MANIFEST_BYTES
        or prefix["manifest_sha256"] != PREFIX_MANIFEST_SHA256
        or any(str(row["path"]).startswith("locks/") for row in prefix_rows)
    ):
        raise TDG8RCV3FreezeError("projected prefix manifest differs")
    return {
        "store_path": SOURCE_STORE_PATH,
        "directories": list(first_directories),
        "full_manifest": full,
        "prefix_manifest": prefix,
        "exclusion": _expected_exclusion(),
        "anchor": _derive_anchor(store, first_rows),
    }


def _path_absent_nofollow(repository: Path, relative: str) -> bool:
    """Return absence without following a symlink in any existing component."""
    parts = _safe_relative(relative)
    descriptors: list[int] = []
    try:
        root_fd = os.open(
            Path(repository),
            os.O_RDONLY
            | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_NOFOLLOW", 0),
        )
        descriptors.append(root_fd)
        directory_fd = root_fd
        for index, part in enumerate(parts):
            try:
                metadata = os.stat(part, dir_fd=directory_fd, follow_symlinks=False)
            except FileNotFoundError:
                return True
            if stat.S_ISLNK(metadata.st_mode):
                raise TDG8RCV3FreezeError("destination path contains a symlink")
            if index == len(parts) - 1:
                return False
            if not stat.S_ISDIR(metadata.st_mode):
                raise TDG8RCV3FreezeError("destination ancestor is not a directory")
            directory_fd = _open_directory(directory_fd, part, "destination ancestor")
            descriptors.append(directory_fd)
        return False
    except TDG8RCV3FreezeError:
        raise
    except OSError as exc:
        raise TDG8RCV3FreezeError("cannot safely inspect destination") from exc
    finally:
        for descriptor in reversed(descriptors):
            os.close(descriptor)


def destination_absence_status(repository: Path) -> tuple[bool, bool]:
    """Return wrapper/store absence while rejecting unsafe path components."""
    wrapper_absent = _path_absent_nofollow(repository, DESTINATION_WRAPPER_PATH)
    store_absent = _path_absent_nofollow(repository, DESTINATION_STORE_PATH)
    return wrapper_absent, store_absent


def observe_destination_absence(repository: Path) -> dict[str, object]:
    wrapper_absent, store_absent = destination_absence_status(repository)
    if not wrapper_absent or not store_absent:
        raise TDG8RCV3FreezeError("RCV3 destination was not absent")
    return expected_destination_observation()


def derive_live_evidence(repository: Path) -> dict[str, object]:
    """Perform the sole live, read-only FRZ1 observation."""
    return {
        "repair": derive_repair_evidence(repository),
        "bootstrap": derive_bootstrap_evidence(repository),
        "source": derive_source_evidence(repository),
        "destination_absence": observe_destination_absence(repository),
        "scope": dict(SCOPE),
        "claims": dict(CLAIMS),
        "nonclaims": list(NONCLAIMS),
    }


def _validate_repair(value: object) -> dict[str, object]:
    if not isinstance(value, Mapping) or set(value) != {
        "commit", "head_at_observation", "repair_commit_is_ancestor",
        "live_sources_equal_commit", "sources",
    }:
        raise TDG8RCV3FreezeError("repair evidence schema differs")
    head = value.get("head_at_observation")
    if head != REPAIR_COMMIT:
        raise TDG8RCV3FreezeError("repair observation HEAD differs")
    if (
        value.get("commit") != REPAIR_COMMIT
        or value.get("repair_commit_is_ancestor") is not True
        or value.get("live_sources_equal_commit") is not True
        or value.get("sources") != _expected_repair_sources()
    ):
        raise TDG8RCV3FreezeError("repair evidence differs")
    return dict(value)


def _validate_source(value: object) -> dict[str, object]:
    if not isinstance(value, Mapping) or set(value) != {
        "store_path", "directories", "full_manifest", "prefix_manifest",
        "exclusion", "anchor",
    }:
        raise TDG8RCV3FreezeError("source evidence schema differs")
    if (
        value.get("store_path") != SOURCE_STORE_PATH
        or value.get("directories") != list(EXPECTED_DIRECTORIES)
        or value.get("exclusion") != _expected_exclusion()
        or value.get("anchor") != expected_source_anchor()
    ):
        raise TDG8RCV3FreezeError("source evidence identity differs")
    full_rows = _validate_summary(
        value.get("full_manifest"),
        label="full source",
        count=FULL_MANIFEST_COUNT,
        byte_count=FULL_MANIFEST_BYTES,
        digest=FULL_MANIFEST_SHA256,
    )
    prefix_rows = _validate_summary(
        value.get("prefix_manifest"),
        label="projected prefix",
        count=PREFIX_MANIFEST_COUNT,
        byte_count=PREFIX_MANIFEST_BYTES,
        digest=PREFIX_MANIFEST_SHA256,
    )
    full_map = {str(row["path"]): row for row in full_rows}
    expected_prefix = [row for row in full_rows if row["path"] not in EXCLUDED_PATHS]
    if prefix_rows != expected_prefix or set(EXCLUDED_PATHS) != set(full_map) - {
        str(row["path"]) for row in prefix_rows
    }:
        raise TDG8RCV3FreezeError("prefix is not the exact full-manifest projection")
    if any(str(row["path"]).startswith("locks/") for row in prefix_rows):
        raise TDG8RCV3FreezeError("projected prefix contains a lock")
    return dict(value)


def validate_evidence(value: object) -> dict[str, object]:
    if not isinstance(value, Mapping) or set(value) != {
        "repair", "bootstrap", "source", "destination_absence", "scope", "claims", "nonclaims"
    }:
        raise TDG8RCV3FreezeError("FRZ1 evidence schema differs")
    _validate_repair(value.get("repair"))
    if value.get("bootstrap") != expected_bootstrap_contract():
        raise TDG8RCV3FreezeError("bootstrap implementation evidence differs")
    _validate_source(value.get("source"))
    if value.get("destination_absence") != expected_destination_observation():
        raise TDG8RCV3FreezeError("destination-absence evidence differs")
    if (
        value.get("scope") != SCOPE
        or value.get("claims") != CLAIMS
        or value.get("nonclaims") != list(NONCLAIMS)
    ):
        raise TDG8RCV3FreezeError("FRZ1 scope or claims differ")
    return dict(value)


def build_result(config_raw: bytes, evidence: Mapping[str, object]) -> dict[str, object]:
    config = parse_config(config_raw)
    checked = validate_evidence(evidence)
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "project_version": config["project_version"],
        "target_protocol": TARGET_PROTOCOL,
        "classification": CLASSIFICATION,
        "gate_status": "pass",
        "source_config_sha256": _sha(config_raw),
        "artifact_payload": checked,
    }


def parse_result(raw: bytes) -> dict[str, Any]:
    value = _json(raw, "RCV3 FRZ1 result")
    if canonical_result(value) != raw:
        raise TDG8RCV3FreezeError("RCV3 FRZ1 result is noncanonical")
    return value


def validate_compact(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    config = parse_config(config_raw)
    result = parse_result(result_raw)
    if set(result) != {
        "schema_version", "artifact_id", "project_version", "target_protocol",
        "classification", "gate_status", "source_config_sha256", "artifact_payload",
    }:
        raise TDG8RCV3FreezeError("RCV3 FRZ1 result fields differ")
    if (
        result.get("schema_version") != SCHEMA_VERSION
        or result.get("artifact_id") != ARTIFACT_ID
        or result.get("project_version") != config["project_version"]
        or result.get("target_protocol") != TARGET_PROTOCOL
        or result.get("classification") != CLASSIFICATION
        or result.get("gate_status") != "pass"
        or result.get("source_config_sha256") != _sha(config_raw)
    ):
        raise TDG8RCV3FreezeError("RCV3 FRZ1 result identity differs")
    validate_evidence(result.get("artifact_payload"))
    return result


__all__ = [
    "ARTIFACT_ID", "BOOTSTRAP_RUNTIME_PATH", "BOOTSTRAP_RUNTIME_SHA256",
    "BOOTSTRAP_SCRIPT_PATH", "BOOTSTRAP_SCRIPT_SHA256", "CLAIMS", "CLASSIFICATION", "CONFIG_PATH",
    "DESTINATION_STORE_PATH", "DESTINATION_WRAPPER_PATH", "EXCLUDED_PATHS",
    "FULL_MANIFEST_COUNT", "FULL_MANIFEST_SHA256", "GENERATION9_CHECKPOINT_SHA256",
    "NONCLAIMS", "PREFIX_MANIFEST_COUNT", "PREFIX_MANIFEST_SHA256", "PROJECTION_ID",
    "REPAIR_COMMIT", "RESULT_PATH", "SCOPE", "SOURCE_STORE_PATH",
    "TDG8RCV3FreezeError", "build_result", "canonical_result",
    "derive_bootstrap_evidence", "derive_live_evidence", "derive_repair_evidence", "derive_source_evidence",
    "destination_absence_status", "expected_bootstrap_contract",
    "expected_destination_observation", "observe_destination_absence", "parse_config",
    "parse_result", "read_nofollow", "validate_compact", "validate_evidence",
]
