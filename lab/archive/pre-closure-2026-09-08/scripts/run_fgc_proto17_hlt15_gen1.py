#!/usr/bin/env python3
"""Authenticate AUTH1, reimport actual inputs, and create PROTO17 generation zero.

This file is a standard-library-only bootstrap until the external tuple,
committed AUTH1 blob, complete live source closure, and numerical environment
have all been checked.  It then loads the pinned production modules through
minimal namespace packages, bypassing the repository's broad ``__init__``
surfaces.  The only caller-controlled value is the external AUTH1 tuple.
"""
from __future__ import annotations

import argparse
import ast
from hashlib import sha256
import importlib.machinery
import importlib.util
import json
import os
from pathlib import Path
import platform
import re
import stat
import subprocess
import sys
import tomllib
from types import ModuleType
from typing import Any, Mapping, TYPE_CHECKING

if TYPE_CHECKING:  # Static closure declaration; never executed by HLT15.
    from recursive_horizons.fgc.evolution.proto17_hlt15_runtime import Proto17HLT15Runtime
    from recursive_horizons.fgc.evolution.proto18_auth1_inputs import Auth1GenesisInputs
    from recursive_horizons.fgc.evolution.proto18_authority import Proto18AuthorityError


ROOT = Path(__file__).resolve().parents[1]


class Proto17HLT15LaunchError(ValueError):
    """HLT15 authority, environment, import, or input admission failed."""


_HEX40 = re.compile(r"^[0-9a-f]{40}$")
_HEX64 = re.compile(r"^[0-9a-f]{64}$")
_TUPLE_FIELDS = frozenset({
    "authorization_commit", "authorization_result_path",
    "authorization_result_sha256", "genesis_spec_sha256",
})
_PAYLOAD_FIELDS = frozenset({
    "pref26_binding", "authority_inputs", "auth1_input_evidence",
    "auth1_input_evidence_sha256", "source_hashes", "campaign_id",
    "genesis_spec", "genesis_spec_sha256", "source_closure",
    "source_closure_sha256", "numerical_environment", "claims",
})
_TRUE_CLAIMS = frozenset({
    "canonical_payload_constructed", "PREF26_evidence_bound",
    "complete_calibration_GenesisSpec_constructed",
    "calibration_launch_closure_pinned",
})
_FALSE_CLAIMS = frozenset({
    "committed_external_Git_live_import_validation_completed", "HLT15_authorized",
    "future_output_roots_observed", "namespace_created",
    "pretrajectory_operation_authorized", "fresh_GR0_dynamic_calibration_completed",
    "candidate_execution_authorized", "physical_transition_claim_authorized",
})
_DIRECT_PATHS = {
    "runtime_module": "src/recursive_horizons/fgc/evolution/proto17_hlt15_runtime.py",
    "adapter_module": "src/recursive_horizons/fgc/evolution/proto18_auth1_inputs.py",
    "runner": "scripts/run_fgc_proto17_hlt15_gen1.py",
}
_MODULE_LOAD_ORDER = (
    ("recursive_horizons.fgc.evolution.protocol_v17",
     "src/recursive_horizons/fgc/evolution/protocol_v17.py"),
    ("recursive_horizons.fgc.evolution.proto17_pure_construction",
     "src/recursive_horizons/fgc/evolution/proto17_pure_construction.py"),
    ("recursive_horizons.fgc.evolution.proto17_hlt15_runtime",
     "src/recursive_horizons/fgc/evolution/proto17_hlt15_runtime.py"),
    ("recursive_horizons.fgc.evolution.numerical_engine",
     "src/recursive_horizons/fgc/evolution/numerical_engine.py"),
    ("recursive_horizons.fgc.evolution.proto18_production_inputs",
     "src/recursive_horizons/fgc/evolution/proto18_production_inputs.py"),
    ("recursive_horizons.fgc.evolution.proto18_auth1_inputs",
     "src/recursive_horizons/fgc/evolution/proto18_auth1_inputs.py"),
    ("recursive_horizons.fgc.evolution.proto18_authority",
     "src/recursive_horizons/fgc/evolution/proto18_authority.py"),
)
_EXPECTED_CLOSURE_PATHS = frozenset(path for _module, path in _MODULE_LOAD_ORDER) | {
    _DIRECT_PATHS["runner"],
}
_AUTHORITY_INPUT_PATHS = {
    "authority_config": "configs/fgc/fgc-1-pro18-auth1.toml",
    "raw_source_config": "configs/fgc/fgc-1-pro18-frz1.toml",
}


def _compact(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")


def _digest(value: object) -> str:
    return sha256(_compact(value)).hexdigest()


def _indented(value: object) -> bytes:
    return (json.dumps(
        value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False,
    ) + "\n").encode("utf-8")


def _decode(raw: bytes, *, label: str, style: str) -> dict[str, Any]:
    def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(key)
            result[key] = value
        return result

    try:
        value = json.loads(raw.decode("utf-8"), object_pairs_hook=no_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise Proto17HLT15LaunchError(f"{label} is malformed") from exc
    expected = _indented(value) if style == "indented" else _compact(value)
    if not isinstance(value, dict) or expected != raw:
        raise Proto17HLT15LaunchError(f"{label} is noncanonical")
    return value


def _sha(value: object, *, label: str) -> str:
    if not isinstance(value, str) or not _HEX64.fullmatch(value):
        raise Proto17HLT15LaunchError(f"{label} must be lowercase SHA-256")
    return value


def _relative(value: object, *, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise Proto17HLT15LaunchError(f"{label} path is absent")
    path = Path(value)
    if path.is_absolute() or "\\" in value or any(part in {"", ".", ".."} for part in path.parts):
        raise Proto17HLT15LaunchError(f"{label} path is unsafe")
    return value


def _git(root: Path, *args: str) -> bytes:
    try:
        return subprocess.run(
            ["git", *args], cwd=root, check=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
    except (OSError, subprocess.CalledProcessError) as exc:
        raise Proto17HLT15LaunchError("HLT15 Git authority operation failed") from exc


def _require_isolated_detached_image(root: Path, commit: str) -> None:
    """Require the exact detached image declared by PRO18.

    This check intentionally rejects an ordinary descendant checkout.  HLT15
    is a launch operation, so its executable image is the authenticated AUTH1
    commit itself rather than merely a tree that happens to contain it.
    """
    if _git(root, "rev-parse", "HEAD").decode().strip() != commit:
        raise Proto17HLT15LaunchError("HLT15 HEAD is not the exact AUTH1 commit")
    detached = subprocess.run(
        ["git", "symbolic-ref", "-q", "HEAD"], cwd=root,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )
    if detached.returncode == 0:
        raise Proto17HLT15LaunchError("HLT15 launch image is not detached")
    if detached.returncode not in {1}:
        raise Proto17HLT15LaunchError("HLT15 detached-image check failed")


def _require_tracked_clean_image(root: Path) -> None:
    """Reject index/worktree drift without enumerating untracked paths.

    diff checks intentionally ignore untracked paths: raw predecessor archives
    are untracked inputs, and enumerating untracked paths here could observe
    the future PROTO17 root before its explicit runtime boundary.
    """
    for arguments in (
        ("diff", "--quiet", "--no-ext-diff", "--ignore-submodules=none", "--"),
        ("diff", "--cached", "--quiet", "--no-ext-diff", "--ignore-submodules=none", "--"),
    ):
        observed = subprocess.run(
            ["git", *arguments], cwd=root,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        if observed.returncode == 1:
            raise Proto17HLT15LaunchError(
                "HLT15 detached launch image is not tracked-clean"
            )
        if observed.returncode != 0:
            raise Proto17HLT15LaunchError("HLT15 tracked-image check failed")


def _tracked_regular_blob(root: Path, commit: str, relative: str) -> None:
    entries = [entry for entry in _git(
        root, "ls-tree", "-r", "-z", commit, "--", relative,
    ).split(b"\0") if entry]
    if len(entries) != 1:
        raise Proto17HLT15LaunchError("authority path is not one tracked blob")
    try:
        mode_type_hash, observed_path = entries[0].split(b"\t", 1)
        mode, kind, _object = mode_type_hash.split(b" ", 2)
    except ValueError as exc:
        raise Proto17HLT15LaunchError("authority tree entry is malformed") from exc
    if (mode not in {b"100644", b"100755"} or kind != b"blob"
            or observed_path.decode("utf-8") != relative):
        raise Proto17HLT15LaunchError("authority path is not a tracked regular blob")


def _live_regular_bytes(root: Path, relative: str) -> bytes:
    """Read one repository-relative regular file without following symlinks."""
    parts = Path(_relative(relative, label="live authority input")).parts
    root_stat = root.lstat()
    if stat.S_ISLNK(root_stat.st_mode) or not stat.S_ISDIR(root_stat.st_mode):
        raise Proto17HLT15LaunchError("repository root is unsafe")
    flags_dir = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_DIRECTORY", 0)
    flags_file = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    opened: list[int] = []
    descriptor: int | None = None
    try:
        parent = os.open(root, flags_dir)
        opened.append(parent)
        for part in parts[:-1]:
            before = os.stat(part, dir_fd=parent, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                raise Proto17HLT15LaunchError("authority path has an unsafe ancestor")
            child = os.open(part, flags_dir, dir_fd=parent)
            after = os.fstat(child)
            if (after.st_dev, after.st_ino) != (before.st_dev, before.st_ino):
                raise Proto17HLT15LaunchError("authority path changed during traversal")
            opened.append(child)
            parent = child
        before = os.stat(parts[-1], dir_fd=parent, follow_symlinks=False)
        if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
            raise Proto17HLT15LaunchError("authority path is not a live regular file")
        descriptor = os.open(parts[-1], flags_file, dir_fd=parent)
        actual = os.fstat(descriptor)
        if (actual.st_dev, actual.st_ino) != (before.st_dev, before.st_ino):
            raise Proto17HLT15LaunchError("authority file changed before read")
        chunks: list[bytes] = []
        while True:
            chunk = os.read(descriptor, 1 << 20)
            if not chunk:
                break
            chunks.append(chunk)
        payload = b"".join(chunks)
        after = os.fstat(descriptor)
        if (len(payload) != before.st_size
                or (after.st_size, after.st_mtime_ns, after.st_ctime_ns)
                   != (before.st_size, before.st_mtime_ns, before.st_ctime_ns)):
            raise Proto17HLT15LaunchError("authority file changed during read")
        return payload
    except Proto17HLT15LaunchError:
        raise
    except OSError as exc:
        raise Proto17HLT15LaunchError("authority path cannot be read safely") from exc
    finally:
        if descriptor is not None:
            os.close(descriptor)
        for item in reversed(opened):
            os.close(item)


def _committed_and_live(
    root: Path, commit: str, relative: str, expected_sha256: str,
) -> bytes:
    relative = _relative(relative, label="authority blob")
    expected = _sha(expected_sha256, label=f"authority blob {relative}")
    _tracked_regular_blob(root, commit, relative)
    historical = _git(root, "show", f"{commit}:{relative}")
    if sha256(historical).hexdigest() != expected:
        raise Proto17HLT15LaunchError("committed authority blob hash differs")
    if _live_regular_bytes(root, relative) != historical:
        raise Proto17HLT15LaunchError("live authority blob differs from its commit")
    return historical


def _contains_forbidden(value: object) -> bool:
    forbidden = {
        "authorization_commit", "authorization_result_path",
        "authorization_result_sha256",
    }
    if isinstance(value, Mapping):
        return any(key in forbidden or _contains_forbidden(item) for key, item in value.items())
    return isinstance(value, list) and any(_contains_forbidden(item) for item in value)


def _module_path(module: str) -> str | None:
    if not module.startswith("recursive_horizons."):
        return None
    return f"src/{module.replace('.', '/')}.py"


def _static_local_imports(relative: str, source: bytes) -> set[str]:
    try:
        tree = ast.parse(source, filename=relative)
    except SyntaxError as exc:
        raise Proto17HLT15LaunchError("pinned source does not parse") from exc
    current = None
    if relative.startswith("src/") and relative.endswith(".py"):
        current = relative.removeprefix("src/").removesuffix(".py").replace("/", ".")
    package = current.rsplit(".", 1)[0] if current and "." in current else ""
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                target = _module_path(alias.name)
                if target:
                    found.add(target)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                if not package:
                    continue
                pieces = package.split(".")
                base = ".".join(pieces[: len(pieces) - node.level + 1])
                module = f"{base}.{node.module}" if node.module else base
            else:
                module = node.module or ""
            target = _module_path(module)
            if target:
                found.add(target)
            if node.module is None:
                for alias in node.names:
                    child = _module_path(f"{module}.{alias.name}")
                    if child:
                        found.add(child)
    return found


def observed_environment() -> dict[str, object]:
    """Observe the evolution image without importing repository-local code."""
    try:
        import numpy as np
    except ImportError as exc:
        raise Proto17HLT15LaunchError("NumPy evolution environment is unavailable") from exc
    configuration = np.show_config(mode="dicts")
    encoded = json.dumps(configuration, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return {
        "contract_id": "FGC-1-PRO18-AUTH1-evolution-environment-v1",
        "implementation": "CPython" if sys.implementation.name == "cpython" else sys.implementation.name,
        "python_version": platform.python_version(),
        "python_build": " | ".join(platform.python_build()),
        "python_compiler": platform.python_compiler(),
        "numpy_version": np.__version__,
        "numpy_build_configuration": encoded,
        "blas_lapack": "Accelerate" if "accelerate" in encoded.lower() else "other",
        "operating_system": platform.platform(),
        "system": platform.system(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "platform": sys.platform,
        "byteorder": sys.byteorder,
    }


def _validate_environment(payload: Mapping[str, Any]) -> None:
    environment = payload.get("numerical_environment")
    if (not isinstance(environment, Mapping)
            or set(environment) != {"contract", "canonical_sha256"}
            or not isinstance(environment["contract"], Mapping)
            or _digest(environment["contract"]) != environment["canonical_sha256"]):
        raise Proto17HLT15LaunchError("AUTH1 environment contract differs")
    if observed_environment() != dict(environment["contract"]):
        raise Proto17HLT15LaunchError("live numerical environment differs from AUTH1")


def verify_authority_tuple(
    repository_root: Path,
    authority_tuple: Mapping[str, Any],
    *,
    check_environment: bool = True,
) -> dict[str, Any]:
    """Verify external Git/blob/live authority without observing an output root."""
    if not isinstance(authority_tuple, Mapping) or set(authority_tuple) != _TUPLE_FIELDS:
        raise Proto17HLT15LaunchError("external AUTH1 tuple fields differ")
    root = Path(repository_root).resolve()
    commit = authority_tuple["authorization_commit"]
    result_path = _relative(authority_tuple["authorization_result_path"], label="AUTH1 result")
    result_sha = _sha(authority_tuple["authorization_result_sha256"], label="AUTH1 result")
    tuple_genesis_sha = _sha(authority_tuple["genesis_spec_sha256"], label="GenesisSpec")
    if not isinstance(commit, str) or not _HEX40.fullmatch(commit):
        raise Proto17HLT15LaunchError("AUTH1 commit identity differs")
    if _git(root, "rev-parse", "--verify", f"{commit}^{{commit}}").decode().strip() != commit:
        raise Proto17HLT15LaunchError("AUTH1 commit does not resolve exactly")
    _require_isolated_detached_image(root, commit)
    raw = _committed_and_live(root, commit, result_path, result_sha)
    record = _decode(raw, label="AUTH1 result", style="indented")
    if set(record) != {"artifact_id", "artifact_payload"} or record["artifact_id"] != "FGC-1-PRO18-AUTH1":
        raise Proto17HLT15LaunchError("AUTH1 result outer schema differs")
    payload = record["artifact_payload"]
    if not isinstance(payload, Mapping) or set(payload) != _PAYLOAD_FIELDS:
        raise Proto17HLT15LaunchError("AUTH1 result payload schema differs")
    if _contains_forbidden(record):
        raise Proto17HLT15LaunchError("AUTH1 result self-references its external authority")
    claims = payload["claims"]
    if (not isinstance(claims, Mapping)
            or set(claims) != _TRUE_CLAIMS | _FALSE_CLAIMS
            or any(claims[name] is not True for name in _TRUE_CLAIMS)
            or any(claims[name] is not False for name in _FALSE_CLAIMS)):
        raise Proto17HLT15LaunchError("AUTH1 claims were promoted or weakened")

    authority_inputs = payload["authority_inputs"]
    if not isinstance(authority_inputs, Mapping) or set(authority_inputs) != set(_AUTHORITY_INPUT_PATHS):
        raise Proto17HLT15LaunchError("AUTH1 input binding inventory differs")
    bound_inputs: dict[str, bytes] = {}
    for name, expected_path in _AUTHORITY_INPUT_PATHS.items():
        item = authority_inputs[name]
        if (not isinstance(item, Mapping) or set(item) != {"path", "sha256"}
                or item["path"] != expected_path):
            raise Proto17HLT15LaunchError("AUTH1 input binding differs")
        bound_inputs[name] = _committed_and_live(root, commit, item["path"], item["sha256"])

    pref26 = payload["pref26_binding"]
    required_pref26 = {
        "historical_commit", "compact_result", "compact_result_sha256",
        "AUTH1_input_evidence_sha256",
    }
    if not isinstance(pref26, Mapping) or set(pref26) != required_pref26:
        raise Proto17HLT15LaunchError("PREF26 binding differs")
    historical_commit = pref26["historical_commit"]
    if not isinstance(historical_commit, str) or not _HEX40.fullmatch(historical_commit):
        raise Proto17HLT15LaunchError("PREF26 historical commit differs")
    try:
        subprocess.run(
            ["git", "merge-base", "--is-ancestor", historical_commit, commit], cwd=root,
            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise Proto17HLT15LaunchError("PREF26 commit is not an AUTH1 ancestor") from exc
    pref26_path = _relative(pref26["compact_result"], label="PREF26 compact result")
    pref26_sha = _sha(pref26["compact_result_sha256"], label="PREF26 compact result")
    pref26_historical = _committed_and_live(root, historical_commit, pref26_path, pref26_sha)
    if _git(root, "show", f"{commit}:{pref26_path}") != pref26_historical:
        raise Proto17HLT15LaunchError("AUTH1 commit does not preserve PREF26 bytes")
    pref26_record = _decode(pref26_historical, label="PREF26 compact result", style="indented")
    pref26_payload = pref26_record.get("artifact_payload")
    if (not isinstance(pref26_payload, Mapping)
            or pref26_payload.get("AUTH1_input_evidence") != payload["auth1_input_evidence"]):
        raise Proto17HLT15LaunchError("PREF26 embedded AUTH1 evidence differs")

    closure = payload["source_closure"]
    if not isinstance(closure, list) or _digest(closure) != payload["source_closure_sha256"]:
        raise Proto17HLT15LaunchError("AUTH1 source closure digest differs")
    seen: set[str] = set()
    role_map: dict[str, Mapping[str, Any]] = {}
    closure_bytes: dict[str, bytes] = {}
    for entry in closure:
        if not isinstance(entry, Mapping) or set(entry) != {"path", "sha256", "role"}:
            raise Proto17HLT15LaunchError("AUTH1 source closure record differs")
        path = _relative(entry["path"], label="AUTH1 source closure")
        role = entry["role"]
        if (path in seen or role not in {"runtime_module", "adapter_module", "runner", "dependency"}
                or not isinstance(role, str)):
            raise Proto17HLT15LaunchError("AUTH1 source closure identity differs")
        if role != "dependency":
            if role in role_map or _DIRECT_PATHS.get(role) != path:
                raise Proto17HLT15LaunchError("AUTH1 direct source role differs")
            role_map[role] = entry
        seen.add(path)
        closure_bytes[path] = _committed_and_live(root, commit, path, entry["sha256"])
    if seen != _EXPECTED_CLOSURE_PATHS or set(role_map) != set(_DIRECT_PATHS):
        raise Proto17HLT15LaunchError("AUTH1 complete source closure differs")
    for path, source in closure_bytes.items():
        if not _static_local_imports(path, source) <= seen:
            raise Proto17HLT15LaunchError("AUTH1 source closure omits a local import")

    source_hashes = payload["source_hashes"]
    expected_hash_fields = {
        "protocol_config_sha256", "protocol_freeze_result_sha256", "hlt13_result_sha256",
        "run_plan_sha256", "runtime_module_sha256", "adapter_module_sha256", "runner_sha256",
    }
    genesis = payload["genesis_spec"]
    if (not isinstance(source_hashes, Mapping) or set(source_hashes) != expected_hash_fields
            or not isinstance(genesis, Mapping)
            or _digest(genesis) != payload["genesis_spec_sha256"]
            or payload["genesis_spec_sha256"] != tuple_genesis_sha):
        raise Proto17HLT15LaunchError("AUTH1 GenesisSpec or source-hash identity differs")
    for role, field in (
        ("runtime_module", "runtime_module_sha256"),
        ("adapter_module", "adapter_module_sha256"),
        ("runner", "runner_sha256"),
    ):
        if (role_map[role]["sha256"] != source_hashes[field]
                or genesis.get(field) != source_hashes[field]):
            raise Proto17HLT15LaunchError("AUTH1 direct source hash differs")
    evidence = payload["auth1_input_evidence"]
    if (not isinstance(evidence, Mapping)
            or evidence.get("AUTH1_input_evidence_sha256") != payload["auth1_input_evidence_sha256"]
            or payload["auth1_input_evidence_sha256"] != pref26["AUTH1_input_evidence_sha256"]):
        raise Proto17HLT15LaunchError("AUTH1 compact PREF26 evidence pin differs")
    bare_evidence = dict(evidence)
    bare_evidence.pop("AUTH1_input_evidence_sha256", None)
    if _digest(bare_evidence) != payload["auth1_input_evidence_sha256"]:
        raise Proto17HLT15LaunchError("AUTH1 compact PREF26 evidence digest differs")
    if check_environment:
        _validate_environment(payload)
    _require_tracked_clean_image(root)
    return {
        "record": record,
        "payload": dict(payload),
        "authority_config_bytes": bound_inputs["authority_config"],
        "raw_source_config_bytes": bound_inputs["raw_source_config"],
        "pref26_record": pref26_record,
        "closure_bytes": closure_bytes,
        "authority_tuple": dict(authority_tuple),
    }


def _install_production_modules(root: Path, closure: Mapping[str, bytes]) -> dict[str, ModuleType]:
    contamination = sorted(
        name for name in sys.modules
        if name == "recursive_horizons" or name.startswith("recursive_horizons.")
    )
    if contamination:
        raise Proto17HLT15LaunchError("repository-local modules were imported before authority")
    packages = (
        ("recursive_horizons", root / "src/recursive_horizons"),
        ("recursive_horizons.fgc", root / "src/recursive_horizons/fgc"),
        ("recursive_horizons.fgc.evolution", root / "src/recursive_horizons/fgc/evolution"),
    )
    for name, path in packages:
        module = ModuleType(name)
        module.__path__ = [str(path)]  # type: ignore[attr-defined]
        module.__package__ = name
        module.__spec__ = importlib.machinery.ModuleSpec(name, loader=None, is_package=True)
        sys.modules[name] = module
    loaded: dict[str, ModuleType] = {}
    try:
        for name, relative in _MODULE_LOAD_ORDER:
            if relative not in closure:
                raise Proto17HLT15LaunchError("pinned module is absent from verified closure")
            path = root / relative
            spec = importlib.util.spec_from_file_location(name, path)
            if spec is None or spec.loader is None or spec.origin is None:
                raise Proto17HLT15LaunchError("pinned module import specification differs")
            if Path(spec.origin).resolve() != path.resolve():
                raise Proto17HLT15LaunchError("pinned module import origin differs")
            module = importlib.util.module_from_spec(spec)
            sys.modules[name] = module
            try:
                code = compile(closure[relative], str(path), "exec", dont_inherit=True)
                exec(code, module.__dict__)
            except BaseException:
                sys.modules.pop(name, None)
                raise
            if Path(module.__file__).resolve() != path.resolve():
                raise Proto17HLT15LaunchError("loaded module origin differs")
            loaded[name] = module
        allowed = {name for name, _path in packages} | {name for name, _path in _MODULE_LOAD_ORDER}
        observed = {
            name for name in sys.modules
            if name == "recursive_horizons" or name.startswith("recursive_horizons.")
        }
        if observed != allowed:
            raise Proto17HLT15LaunchError("undeclared repository-local module was imported")
        return loaded
    except BaseException:
        for name in list(sys.modules):
            if name == "recursive_horizons" or name.startswith("recursive_horizons."):
                sys.modules.pop(name, None)
        raise


def run_hlt15_gen1(repository_root: Path, authority_tuple: Mapping[str, Any]) -> dict[str, Any]:
    """Run all gates in one invocation; materialization is the final operation."""
    root = Path(repository_root).resolve()
    verified = verify_authority_tuple(root, authority_tuple, check_environment=True)
    modules = _install_production_modules(root, verified["closure_bytes"])
    authority = modules["recursive_horizons.fgc.evolution.proto18_authority"]
    adapter = modules["recursive_horizons.fgc.evolution.proto18_auth1_inputs"]
    runtime_module = modules["recursive_horizons.fgc.evolution.proto17_hlt15_runtime"]
    try:
        authority_config = tomllib.loads(verified["authority_config_bytes"].decode("utf-8"))
        raw_config = tomllib.loads(verified["raw_source_config_bytes"].decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise Proto17HLT15LaunchError("AUTH1 committed configuration is unreadable") from exc
    try:
        reconstructed = authority.verify_authority_payload(
            root, authority_config, verified["record"],
        )
        if reconstructed != verified["record"]:
            raise Proto17HLT15LaunchError("AUTH1 full content reconstruction differs")
        payload = verified["payload"]
        genesis = adapter.build_calibration_genesis(
            payload["auth1_input_evidence"], campaign_id=payload["campaign_id"],
            namespace="runs/fgc-2-sf1/proto17/calibration",
            source_hashes=payload["source_hashes"],
            numerical_environment=payload["numerical_environment"]["contract"],
        )
        if (genesis.genesis_spec != payload["genesis_spec"]
                or _digest(genesis.genesis_spec) != payload["genesis_spec_sha256"]):
            raise Proto17HLT15LaunchError("AUTH1 GenesisSpec rederivation differs")
        source_identities = adapter.source_identities_from_pref26_result(
            verified["pref26_record"], payload["auth1_input_evidence"],
        )
        adapter.reimport_matches_pref26_evidence(
            root, raw_config, payload["auth1_input_evidence"],
            expected_source_identities=source_identities,
        )
    except Proto17HLT15LaunchError:
        raise
    except Exception as exc:
        raise Proto17HLT15LaunchError("HLT15 authority reconstruction or raw reimport failed") from exc
    launch = verified["authority_tuple"]
    binding = {
        "artifact_id": "FGC-1-PRO18-AUTH1",
        **launch,
        "source_closure_sha256": payload["source_closure_sha256"],
        "numerical_environment_sha256": payload["numerical_environment"]["canonical_sha256"],
        "pref26_auth1_input_evidence_sha256": payload["auth1_input_evidence_sha256"],
    }
    runtime = runtime_module.Proto17HLT15Runtime(
        root, payload["genesis_spec"], authority_binding=binding,
    )
    return dict(runtime.materialize_generation_zero())


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("authority_tuple", type=Path, help="canonical external AUTH1 tuple JSON")
    args = parser.parse_args()
    authority_tuple = _decode(
        args.authority_tuple.read_bytes(), label="external AUTH1 tuple", style="compact",
    )
    result = run_hlt15_gen1(ROOT, authority_tuple)
    print(json.dumps({
        "artifact_id": result["receipt"]["artifact_id"],
        "checkpoint_sha256": result["checkpoint"]["checkpoint_sha256"],
        "receipt_sha256": result["receipt"]["receipt_sha256"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
