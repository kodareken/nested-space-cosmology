"""Pure, no-follow tooling for the PRO19 first-event launch manifest.

The manifest is deliberately only an *external-image inventory*.  It neither
opens a campaign store nor imports a numerical operator.  ``--refresh`` is
the coordinator-only step which replaces visible placeholders with hashes;
``--verify`` rejects every placeholder, so a scaffold cannot become launch
authority by accident.
"""
from __future__ import annotations

from hashlib import sha256
import os
from pathlib import Path
import stat
import tomllib
from typing import Any, Mapping

MANIFEST_SCHEMA = "FGC-1-PRO19-launch-authority-v1"
REQUIRED_ROLES = frozenset({"runtime", "config", "result", "documentation", "reproducer", "runner", "test"})
PENDING_SHA256 = "PENDING_COORDINATOR_REFRESH"
AUTH1_PATH = "results/fgc-1-pro18-auth1.json"
MON16_CONFIG_PATH = "configs/fgc/fgc-1-hlt16-mon16.toml"
MON16_RESULT_PATH = "results/fgc-1-hlt16-mon16.json"


class Proto19LaunchManifestError(ValueError):
    """The manifest cannot serve as a complete, safe launch inventory."""


def _relative_path(value: object) -> str:
    if not isinstance(value, str) or not value or "\\" in value:
        raise Proto19LaunchManifestError("manifest path is not portable")
    path = Path(value)
    if path.is_absolute() or "." in path.parts or ".." in path.parts:
        raise Proto19LaunchManifestError("manifest path escapes the repository")
    if path.parts[0] not in {"configs", "docs", "results", "scripts", "src", "tests"}:
        raise Proto19LaunchManifestError("manifest path has an unauthorized root")
    return path.as_posix()


def nofollow_regular_bytes(root: Path, relative: str) -> bytes:
    """Read one unique regular file without following any path component."""
    parts = Path(relative).parts
    directory_fd = leaf_fd = -1
    try:
        directory_fd = os.open(root, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0))
        if not stat.S_ISDIR(os.fstat(directory_fd).st_mode):
            raise Proto19LaunchManifestError("repository root is not a directory")
        for component in parts[:-1]:
            child_fd = os.open(component, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0), dir_fd=directory_fd)
            if not stat.S_ISDIR(os.fstat(child_fd).st_mode):
                os.close(child_fd)
                raise Proto19LaunchManifestError("manifest path contains a non-directory")
            os.close(directory_fd)
            directory_fd = child_fd
        leaf_fd = os.open(parts[-1], os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0), dir_fd=directory_fd)
        metadata = os.fstat(leaf_fd)
        if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
            raise Proto19LaunchManifestError("manifest binding is not a unique regular file")
        chunks: list[bytes] = []
        while block := os.read(leaf_fd, 1024 * 1024):
            chunks.append(block)
        return b"".join(chunks)
    except FileNotFoundError as error:
        raise Proto19LaunchManifestError("manifest binding is absent") from error
    except OSError as error:
        raise Proto19LaunchManifestError("manifest binding cannot be read safely") from error
    finally:
        if leaf_fd != -1:
            os.close(leaf_fd)
        if directory_fd != -1:
            os.close(directory_fd)


def _hash(value: object, *, allow_pending: bool) -> str:
    if allow_pending and value == PENDING_SHA256:
        return PENDING_SHA256
    if not isinstance(value, str) or len(value) != 64 or value.lower() != value:
        raise Proto19LaunchManifestError("manifest digest differs")
    try:
        int(value, 16)
    except ValueError as error:
        raise Proto19LaunchManifestError("manifest digest differs") from error
    return value


def parse_manifest(raw: bytes, *, manifest_path: str, allow_pending: bool) -> dict[str, Any]:
    """Validate and normalize the exact manifest contract without reading files."""
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise Proto19LaunchManifestError("manifest is malformed TOML") from error
    if not isinstance(value, dict) or set(value) != {"schema", "environment", "authority_path"}:
        raise Proto19LaunchManifestError("manifest top-level fields differ")
    if value["schema"] != MANIFEST_SCHEMA:
        raise Proto19LaunchManifestError("manifest schema differs")
    environment = value["environment"]
    if not isinstance(environment, Mapping) or dict(environment) != {"source": "AUTH1", "auth1_result_path": AUTH1_PATH}:
        raise Proto19LaunchManifestError("manifest environment differs")
    own = _relative_path(manifest_path)
    entries = value["authority_path"]
    if not isinstance(entries, list) or not entries:
        raise Proto19LaunchManifestError("manifest authority inventory is absent")
    normalized: list[dict[str, str]] = []
    paths: set[str] = set()
    roles: set[str] = set()
    for entry in entries:
        if not isinstance(entry, Mapping) or set(entry) != {"role", "path", "sha256"}:
            raise Proto19LaunchManifestError("manifest authority entry fields differ")
        role = entry["role"]
        path = _relative_path(entry["path"])
        if not isinstance(role, str) or role not in REQUIRED_ROLES:
            raise Proto19LaunchManifestError("manifest authority role differs")
        if path == own:
            raise Proto19LaunchManifestError("manifest must not self-bind")
        if path in paths:
            raise Proto19LaunchManifestError("manifest authority path is duplicated")
        paths.add(path); roles.add(role)
        normalized.append({"role": role, "path": path, "sha256": _hash(entry["sha256"], allow_pending=allow_pending)})
    if roles != set(REQUIRED_ROLES):
        raise Proto19LaunchManifestError("manifest does not close every authority role")
    required = {
        AUTH1_PATH,
        "results/fgc-1-pro19-frz1.json",
        MON16_CONFIG_PATH,
        MON16_RESULT_PATH,
    }
    if not required.issubset(paths):
        raise Proto19LaunchManifestError("manifest omits required compact progression evidence")
    return {"schema": MANIFEST_SCHEMA, "environment": dict(environment), "authority_path": sorted(normalized, key=lambda item: item["path"])}


def _require_mon16_config_bindings(root: Path, parsed: Mapping[str, Any]) -> None:
    """Require the manifest to close every path MON16 itself configured.

    The manifest layer deliberately validates only the inventory shape here.
    Launch authority independently compares these entries and result records
    against the captured bytes before an event can be authorized.
    """
    try:
        config = tomllib.loads(nofollow_regular_bytes(root, MON16_CONFIG_PATH).decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise Proto19LaunchManifestError("MON16 config is malformed") from exc
    if not isinstance(config, Mapping) or config.get("artifact_id") != "FGC-1-HLT16-MON16":
        raise Proto19LaunchManifestError("MON16 config identity differs")
    paths = {entry["path"] for entry in parsed["authority_path"]}
    required = {MON16_CONFIG_PATH, MON16_RESULT_PATH}
    for section in ("predecessors", "inventory"):
        rows = config.get(section)
        if not isinstance(rows, list) or not rows:
            raise Proto19LaunchManifestError(f"MON16 config {section} differs")
        for row in rows:
            if not isinstance(row, Mapping):
                raise Proto19LaunchManifestError(f"MON16 config {section} entry differs")
            try:
                required.add(_relative_path(row.get("path")))
            except Proto19LaunchManifestError as exc:
                raise Proto19LaunchManifestError(f"MON16 config {section} entry differs") from exc
    if not required.issubset(paths):
        raise Proto19LaunchManifestError("manifest omits MON16 configured bindings")


def refresh_manifest(root: Path, raw: bytes, *, manifest_path: str) -> dict[str, Any]:
    """Replace every digest using no-follow reads; all listed paths must exist."""
    parsed = parse_manifest(raw, manifest_path=manifest_path, allow_pending=True)
    _require_mon16_config_bindings(root, parsed)
    refreshed: list[dict[str, str]] = []
    for entry in parsed["authority_path"]:
        refreshed.append({**entry, "sha256": sha256(nofollow_regular_bytes(root, entry["path"])).hexdigest()})
    return {"schema": parsed["schema"], "environment": parsed["environment"], "authority_path": refreshed}


def render_manifest(manifest: Mapping[str, Any]) -> bytes:
    """Render the constrained TOML form deterministically."""
    parsed = parse_manifest(
        _render_unchecked(manifest), manifest_path="configs/fgc/fgc-1-pro19-launch-authority.toml", allow_pending=True,
    )
    return _render_unchecked(parsed)


def _render_unchecked(manifest: Mapping[str, Any]) -> bytes:
    environment = manifest["environment"]
    lines = [f'schema = "{manifest["schema"]}"', "", "[environment]", f'source = "{environment["source"]}"', f'auth1_result_path = "{environment["auth1_result_path"]}"', ""]
    for entry in manifest["authority_path"]:
        lines.extend(("[[authority_path]]", f'role = "{entry["role"]}"', f'path = "{entry["path"]}"', f'sha256 = "{entry["sha256"]}"', ""))
    return ("\n".join(lines)).encode("utf-8")


def verify_manifest(root: Path, raw: bytes, *, manifest_path: str) -> dict[str, Any]:
    """Verify a fully refreshed manifest against safe live reads."""
    parsed = parse_manifest(raw, manifest_path=manifest_path, allow_pending=False)
    _require_mon16_config_bindings(root, parsed)
    for entry in parsed["authority_path"]:
        observed = sha256(nofollow_regular_bytes(root, entry["path"])).hexdigest()
        if observed != entry["sha256"]:
            raise Proto19LaunchManifestError(f"manifest binding differs: {entry['path']}")
    return parsed


__all__ = [
    "AUTH1_PATH", "MON16_CONFIG_PATH", "MON16_RESULT_PATH", "PENDING_SHA256", "Proto19LaunchManifestError", "nofollow_regular_bytes",
    "parse_manifest", "refresh_manifest", "render_manifest", "verify_manifest",
]
