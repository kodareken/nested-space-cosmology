"""Fail-closed HLT14 authority and raw-bundle admission primitives.

This module deliberately stops before campaign materialization: it authenticates
an externally committed authority result and verifies synthetic restart bundles
only.  It never creates a namespace, imports an evolution implementation, or
advances a numerical state.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import re
import struct
import subprocess
from typing import Any, Mapping, Sequence
import zipfile

import numpy as np
from .proto17_pure_construction import Proto17ConstructionError, build_genesis


class HLT14InputError(ValueError):
    """Base class for a non-authorized HLT14 input or provenance boundary."""


class HLT14AuthorityError(HLT14InputError):
    """The external authority tuple cannot authenticate its declared input."""


class HLT14BundleError(HLT14InputError):
    """A raw restart bundle differs from its declared, finite input contract."""


class HLT14NamespaceForbiddenError(HLT14InputError):
    """This admission layer was asked to touch a future production namespace."""


_HEX40 = re.compile(r"^[0-9a-f]{40}$")
_HEX64 = re.compile(r"^[0-9a-f]{64}$")
_ARRAY_NAMES = ("u", "p", "q", "tracer_positions", "tracer_proper_times", "event_proper_times", "event_fields")
_STATE_FIELDS = {"schema_id", "member_key", "method", "point_count", "coordinate_time", "accepted_boundary_time", "step_index", "transaction_serial", "input_hash", "source_retry_count", "CFL_retry_count", "runtime_monitor_state", "causal_state", "tracer_state", "event_history_state", "physical_arrays", "physical_state_sha256", "restart_payload_sha256"}
_DESCRIPTOR_FIELDS = {"member_key", "method", "point_count", "source_checkpoint", "source_bundle_key", "accepted_boundary_time", "step_index", "transaction_serial", "state_object", "state_object_canonical_sha256", "initial_TDG6_ledger", "initial_TDG6_ledger_sha256"}
_BUNDLE_FIELDS = {"source_bundle_key", "source_checkpoint", "relative_input_path", "raw_file_sha256", "required_npz_member_keys"}
_METADATA_FIELDS = {"member_key", "method", "point_count", "coordinate_time", "accepted_boundary_time", "step_index", "transaction_serial", "input_hash", "source_retry_count", "CFL_retry_count", "runtime_monitor_state", "causal_state", "tracer_state", "event_history_state", "initial_TDG6_ledger", "initial_TDG6_ledger_sha256"}
_MAX_METADATA_BYTES = 1 << 20


def canonical_bytes(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
                      allow_nan=False).encode("utf-8")


def canonical_sha256(value: object) -> str:
    return sha256(canonical_bytes(value)).hexdigest()


def decode_canonical_json(raw: bytes, *, context: str) -> dict[str, Any]:
    def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise HLT14AuthorityError(f"{context}: duplicate JSON key")
            result[key] = value
        return result
    try:
        value = json.loads(raw.decode("utf-8"), object_pairs_hook=no_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise HLT14AuthorityError(f"{context}: malformed JSON") from error
    if not isinstance(value, dict) or canonical_bytes(value) != raw:
        raise HLT14AuthorityError(f"{context}: JSON is not canonical")
    return value


def _sha(value: object, *, context: str, error: type[HLT14InputError] = HLT14InputError) -> str:
    if not isinstance(value, str) or not _HEX64.fullmatch(value):
        raise error(f"{context}: expected lowercase SHA-256")
    return value


def safe_relative_path(value: object, *, context: str) -> str:
    if not isinstance(value, str) or not value:
        raise HLT14InputError(f"{context}: path must be nonempty text")
    path = Path(value)
    if path.is_absolute() or "\\" in value or any(part in {"", ".", ".."} for part in value.split("/")):
        raise HLT14InputError(f"{context}: path is not a safe relative path")
    return value


def _under(root: Path, relative: str, *, context: str) -> Path:
    lexical = root / safe_relative_path(relative, context=context)
    if lexical.is_symlink() or any(parent.is_symlink() for parent in lexical.parents if parent != root.parent):
        raise HLT14InputError(f"{context}: symlinked path is forbidden")
    candidate = lexical.resolve(strict=False)
    try:
        candidate.relative_to(root.resolve())
    except ValueError as error:
        raise HLT14InputError(f"{context}: path escapes root") from error
    return candidate


@dataclass(frozen=True, slots=True)
class LaunchAuthorityTuple:
    authorization_commit: str
    authorization_result_path: str
    authorization_result_sha256: str
    genesis_spec_sha256: str

    def __post_init__(self) -> None:
        if not isinstance(self.authorization_commit, str) or not _HEX40.fullmatch(self.authorization_commit):
            raise HLT14AuthorityError("authority commit must be a full lowercase commit")
        try:
            safe_relative_path(self.authorization_result_path, context="authority result")
        except HLT14InputError as error:
            raise HLT14AuthorityError("authority result path is unsafe") from error
        _sha(self.authorization_result_sha256, context="authority result", error=HLT14AuthorityError)
        _sha(self.genesis_spec_sha256, context="GenesisSpec", error=HLT14AuthorityError)


@dataclass(frozen=True, slots=True)
class VerifiedAuthority:
    launch: LaunchAuthorityTuple
    result: Mapping[str, Any]
    genesis_spec: Mapping[str, Any]
    pinned_sources: tuple[Mapping[str, Any], ...]


class LaunchAuthorityVerifier:
    """Authenticate a canonical authority result strictly from a Git commit."""

    def __init__(self, repository_root: Path) -> None:
        self.repository_root = Path(repository_root).resolve()

    def _git(self, *args: str) -> bytes:
        try:
            return subprocess.run(["git", *args], cwd=self.repository_root, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, check=True).stdout
        except (OSError, subprocess.CalledProcessError) as error:
            raise HLT14AuthorityError("authority Git operation failed") from error

    def _tracked_regular_blob(self, commit: str, relative: str) -> None:
        listing = self._git("ls-tree", "-r", "-z", commit, "--", relative).split(b"\0")
        entries = [entry for entry in listing if entry]
        if len(entries) != 1:
            raise HLT14AuthorityError("authority result is not one tracked blob")
        try:
            mode_type_hash, path = entries[0].split(b"\t", 1)
            mode, kind, _object = mode_type_hash.split(b" ", 2)
        except ValueError as error:
            raise HLT14AuthorityError("authority tree entry is malformed") from error
        if mode not in {b"100644", b"100755"} or kind != b"blob" or path.decode("utf-8") != relative:
            raise HLT14AuthorityError("authority result is not a tracked regular blob")

    def verify(self, launch: LaunchAuthorityTuple, *, check_import_origins: bool = True) -> VerifiedAuthority:
        self._git("cat-file", "-e", f"{launch.authorization_commit}^{{commit}}")
        if subprocess.run(["git", "merge-base", "--is-ancestor", launch.authorization_commit, "HEAD"], cwd=self.repository_root).returncode:
            raise HLT14AuthorityError("authority commit is not an ancestor of HEAD")
        self._tracked_regular_blob(launch.authorization_commit, launch.authorization_result_path)
        raw = self._git("show", f"{launch.authorization_commit}:{launch.authorization_result_path}")
        if sha256(raw).hexdigest() != launch.authorization_result_sha256:
            raise HLT14AuthorityError("authority result raw hash differs")
        result = decode_canonical_json(raw, context="authority result")
        required = {"artifact_id", "genesis_spec", "pinned_sources"}
        if (set(result) != required or result.get("artifact_id") != "FGC-1-HLT14-MON14"
                or not isinstance(result["genesis_spec"], dict) or not isinstance(result["pinned_sources"], list)):
            raise HLT14AuthorityError("authority result schema differs")
        # The external tuple, not the authority payload, authenticates these values.
        forbidden = {"authorization_commit", "authorization_result_path", "authorization_result_sha256", "genesis_spec_sha256"}
        if self._contains_forbidden_key(result, forbidden):
            raise HLT14AuthorityError("authority result self-references its external tuple")
        if canonical_sha256(result["genesis_spec"]) != launch.genesis_spec_sha256:
            raise HLT14AuthorityError("GenesisSpec hash differs")
        try:
            derived = build_genesis(result["genesis_spec"])
        except (Proto17ConstructionError, TypeError, ValueError) as error:
            raise HLT14AuthorityError("GenesisSpec is not a complete legal PROTO17 construction") from error
        if derived["genesis_spec"] != result["genesis_spec"]:
            raise HLT14AuthorityError("GenesisSpec differs from its exact normalized derivation")
        pinned = tuple(self._verify_pinned_sources(launch.authorization_commit, result["pinned_sources"], result["genesis_spec"], check_import_origins))
        return VerifiedAuthority(launch=launch, result=result, genesis_spec=derived["genesis_spec"], pinned_sources=pinned)

    @staticmethod
    def _contains_forbidden_key(value: object, forbidden: set[str]) -> bool:
        if isinstance(value, Mapping):
            return any(key in forbidden or LaunchAuthorityVerifier._contains_forbidden_key(item, forbidden) for key, item in value.items())
        if isinstance(value, list):
            return any(LaunchAuthorityVerifier._contains_forbidden_key(item, forbidden) for item in value)
        return False

    def _verify_pinned_sources(self, commit: str, items: Sequence[object], genesis: Mapping[str, object], check_origins: bool) -> list[Mapping[str, Any]]:
        if not items:
            raise HLT14AuthorityError("authority result pins no source modules")
        seen: set[str] = set(); seen_roles: set[str] = set(); checked: list[Mapping[str, Any]] = []
        role_hashes = {"runtime_module": genesis["runtime_module_sha256"], "adapter_module": genesis["adapter_module_sha256"], "runner": genesis["runner_sha256"]}
        for supplied in items:
            if not isinstance(supplied, Mapping) or set(supplied) != {"path", "sha256", "module", "role"}:
                raise HLT14AuthorityError("pinned source schema differs")
            relative = safe_relative_path(supplied["path"], context="pinned source")
            expected = _sha(supplied["sha256"], context="pinned source", error=HLT14AuthorityError)
            module, role = supplied["module"], supplied["role"]
            if relative in seen or role not in role_hashes or role in seen_roles or not isinstance(module, str) or not module or expected != role_hashes[role]:
                raise HLT14AuthorityError("pinned source identity differs")
            seen.add(relative); seen_roles.add(role); self._tracked_regular_blob(commit, relative)
            historical = self._git("show", f"{commit}:{relative}")
            live = _under(self.repository_root, relative, context="pinned source")
            if live.is_symlink() or any(parent.is_symlink() for parent in live.parents if parent != self.repository_root.parent):
                raise HLT14AuthorityError("pinned source symlink is forbidden")
            if not live.is_file() or sha256(historical).hexdigest() != expected or live.read_bytes() != historical:
                raise HLT14AuthorityError("pinned source bytes differ from authority commit")
            if check_origins:
                importlib.invalidate_caches()
                spec = importlib.util.find_spec(module)
                origin = None if spec is None else spec.origin
                if origin is None or Path(origin).resolve() != live.resolve():
                    raise HLT14AuthorityError("pinned module import-origin preflight differs")
            checked.append(dict(supplied))
        if seen_roles != set(role_hashes):
            raise HLT14AuthorityError("authority does not pin all required source roles")
        return checked


def _array_hash(array: np.ndarray) -> str:
    normalized = np.ascontiguousarray(array.astype(np.dtype("<f8"), copy=False))
    return sha256(normalized.tobytes(order="C")).hexdigest()


def shape_framed_array_content_sha256(*arrays: object) -> str:
    """Local exact equivalent of numerical_engine.array_content_sha256.

    Each array contributes little-endian uint32 ndim, little-endian uint64
    shape entries, then C-order little-endian float64 bytes.  This is not a
    hash-of-hashes and deliberately matches the established restart contract.
    """
    digest = sha256()
    for index, value in enumerate(arrays):
        array = np.asarray(value)
        if array.dtype.kind not in {"f"} or not np.isfinite(array).all():
            raise HLT14BundleError(f"aggregate arrays[{index}] is nonfinite or nonfloating")
        normalized = np.asarray(array, dtype="<f8", order="C")
        digest.update(struct.pack("<I", normalized.ndim))
        digest.update(struct.pack(f"<{normalized.ndim}Q", *normalized.shape))
        digest.update(normalized.tobytes(order="C"))
    return digest.hexdigest()


@dataclass(frozen=True, slots=True)
class VerifiedBundle:
    raw_file_sha256: str
    array_hashes: Mapping[str, str]
    physical_state_sha256: str
    restart_payload_sha256: str
    state_object: Mapping[str, Any]
    state_object_canonical_sha256: str


class VerifiedBundleAdapter:
    """Verify a synthetic NPZ against one complete visible descriptor mapping."""

    def __init__(self, input_root: Path) -> None:
        self.input_root = Path(input_root).resolve()

    @classmethod
    def from_manifest_root(
        cls, authority_root: Path, genesis_spec: Mapping[str, object],
    ) -> "VerifiedBundleAdapter":
        """Bind the adapter root to the manifest's authenticated relative root."""
        manifest = genesis_spec.get("physical_input_manifest")
        if not isinstance(manifest, Mapping):
            raise HLT14BundleError("GenesisSpec physical-input manifest is absent")
        try:
            input_root = _under(
                Path(authority_root).resolve(),
                manifest["relative_input_root"],
                context="physical input root",
            )
        except (KeyError, HLT14InputError) as error:
            raise HLT14BundleError("physical input root is not safely bound") from error
        cls.forbid_production_namespace(input_root)
        return cls(input_root)

    @staticmethod
    def forbid_production_namespace(path: Path) -> None:
        if "runs" in Path(path).parts:
            raise HLT14NamespaceForbiddenError("HLT14 input adapter must not touch a production namespace")

    def verify_genesis(self, genesis_spec: Mapping[str, object]) -> Mapping[str, VerifiedBundle]:
        """Admit all six bundles only after complete GenesisSpec derivation."""
        try:
            derived = build_genesis(genesis_spec)
        except (Proto17ConstructionError, TypeError, ValueError) as error:
            raise HLT14BundleError("GenesisSpec is not complete/legal before bundle admission") from error
        if derived["genesis_spec"] != genesis_spec:
            raise HLT14BundleError("GenesisSpec is not its exact normalized derivation")
        entries = {item["source_bundle_key"]: item for item in genesis_spec["physical_input_manifest"]["bundle_members"]}
        return {item["member_key"]: self.verify(item, entries[item["source_bundle_key"]]) for item in genesis_spec["member_descriptors"]}

    def verify(self, member_descriptor: Mapping[str, Any], bundle_entry: Mapping[str, Any]) -> VerifiedBundle:
        if not isinstance(member_descriptor, Mapping) or set(member_descriptor) != _DESCRIPTOR_FIELDS:
            raise HLT14BundleError("PROTO17 member descriptor schema differs")
        if not isinstance(bundle_entry, Mapping) or set(bundle_entry) != _BUNDLE_FIELDS:
            raise HLT14BundleError("PROTO17 physical-manifest bundle schema differs")
        if member_descriptor["source_bundle_key"] != bundle_entry["source_bundle_key"] or member_descriptor["source_checkpoint"] != bundle_entry["source_checkpoint"]:
            raise HLT14BundleError("member descriptor and manifest source ownership differ")
        path = _under(self.input_root, bundle_entry["relative_input_path"], context="raw bundle")
        self.forbid_production_namespace(path)
        if not path.is_file() or path.is_symlink():
            raise HLT14BundleError("raw bundle is not a regular file")
        state = member_descriptor["state_object"]
        if not isinstance(state, Mapping) or set(state) != _STATE_FIELDS:
            raise HLT14BundleError("PROTO17 state object schema differs")
        if sha256(path.read_bytes()).hexdigest() != _sha(bundle_entry["raw_file_sha256"], context="raw bundle", error=HLT14BundleError):
            raise HLT14BundleError("raw bundle file hash differs")
        keys = bundle_entry["required_npz_member_keys"]
        arrays = state["physical_arrays"]
        if not isinstance(keys, list) or not isinstance(arrays, list) or len(arrays) != len(_ARRAY_NAMES):
            raise HLT14BundleError("raw bundle descriptor members differ")
        expected_keys = [item.get("npz_storage_key") for item in arrays] + ["metadata_utf8"]
        if keys != expected_keys or len(set(keys)) != len(keys):
            raise HLT14BundleError("raw bundle required key ordering differs")
        self._check_zip_members(path, keys, arrays)
        hashes: dict[str, str] = {}
        loaded: dict[str, np.ndarray] = {}
        try:
            with np.load(path, allow_pickle=False) as archive:
                if set(archive.files) != set(keys):
                    raise HLT14BundleError("NPZ member set differs")
                raw_metadata = archive["metadata_utf8"]
                if raw_metadata.dtype != np.dtype("uint8") or raw_metadata.ndim != 1 or not raw_metadata.flags.c_contiguous or raw_metadata.nbytes > _MAX_METADATA_BYTES:
                    raise HLT14BundleError("NPZ metadata dtype/rank/layout/size differs")
                metadata = decode_canonical_json(raw_metadata.tobytes(), context="NPZ metadata")
                expected_metadata = {
                    "member_key": member_descriptor["member_key"], "method": member_descriptor["method"],
                    "point_count": member_descriptor["point_count"], "coordinate_time": state["coordinate_time"],
                    "accepted_boundary_time": state["accepted_boundary_time"], "step_index": state["step_index"],
                    "transaction_serial": state["transaction_serial"], "input_hash": state["input_hash"],
                    "source_retry_count": state["source_retry_count"], "CFL_retry_count": state["CFL_retry_count"],
                    "runtime_monitor_state": state["runtime_monitor_state"], "causal_state": state["causal_state"],
                    "tracer_state": state["tracer_state"], "event_history_state": state["event_history_state"],
                    "initial_TDG6_ledger": member_descriptor["initial_TDG6_ledger"],
                    "initial_TDG6_ledger_sha256": member_descriptor["initial_TDG6_ledger_sha256"],
                }
                if metadata != expected_metadata or set(metadata) != _METADATA_FIELDS:
                    raise HLT14BundleError("NPZ metadata does not bind complete descriptor state")
                for name, item in zip(_ARRAY_NAMES, arrays, strict=True):
                    if not isinstance(item, Mapping) or set(item) != {"logical_name", "source_bundle_key", "npz_storage_key", "dtype", "layout", "shape", "little_endian_c_bytes_sha256"}:
                        raise HLT14BundleError("array descriptor schema differs")
                    expected_shape = {
                        "u": [member_descriptor["point_count"], 6],
                        "p": [member_descriptor["point_count"], 6],
                        "q": [member_descriptor["point_count"], 6],
                        "tracer_positions": [48],
                        "tracer_proper_times": [48],
                        "event_proper_times": [24, 48],
                        "event_fields": [24, 48, 6],
                    }[name]
                    array = archive[item["npz_storage_key"]]
                    if item["logical_name"] != name or item["source_bundle_key"] != member_descriptor["source_bundle_key"] or item["dtype"] != "<f8" or item["layout"] != "C":
                        raise HLT14BundleError("array descriptor identity differs")
                    if item["shape"] != expected_shape:
                        raise HLT14BundleError("array descriptor shape differs from the frozen member shape")
                    if array.dtype.str != "<f8" or not array.flags.c_contiguous or list(array.shape) != expected_shape or not np.isfinite(array).all():
                        raise HLT14BundleError("array dtype/layout/shape/finiteness differs")
                    observed = _array_hash(array)
                    if observed != _sha(item["little_endian_c_bytes_sha256"], context="array bytes", error=HLT14BundleError):
                        raise HLT14BundleError("array raw byte hash differs")
                    hashes[name] = observed
                    loaded[name] = array
        except (OSError, ValueError, zipfile.BadZipFile) as error:
            if isinstance(error, HLT14BundleError):
                raise
            raise HLT14BundleError("NPZ bundle is unreadable") from error
        physical = shape_framed_array_content_sha256(loaded["u"], loaded["p"], loaded["q"])
        restart = shape_framed_array_content_sha256(*(loaded[name] for name in _ARRAY_NAMES))
        if physical != _sha(state["physical_state_sha256"], context="physical state", error=HLT14BundleError) or restart != _sha(state["restart_payload_sha256"], context="restart payload", error=HLT14BundleError):
            raise HLT14BundleError("aggregate state hash differs")
        if (state["member_key"] != member_descriptor["member_key"] or state["method"] != member_descriptor["method"] or state["point_count"] != member_descriptor["point_count"] or state["physical_state_sha256"] != physical or state["restart_payload_sha256"] != restart):
            raise HLT14BundleError("rich PROTO17 state object does not bind verified inputs")
        state_hash = canonical_sha256(state)
        if state_hash != _sha(member_descriptor["state_object_canonical_sha256"], context="state object", error=HLT14BundleError):
            raise HLT14BundleError("state object canonical hash differs")
        return VerifiedBundle(raw_file_sha256=bundle_entry["raw_file_sha256"], array_hashes=hashes,
                              physical_state_sha256=physical, restart_payload_sha256=restart,
                              state_object=dict(state), state_object_canonical_sha256=state_hash)

    @staticmethod
    def _check_zip_members(
        path: Path, keys: Sequence[str], array_descriptors: Sequence[Mapping[str, object]],
    ) -> None:
        try:
            with zipfile.ZipFile(path) as archive:
                infos = archive.infolist(); names = [info.filename for info in infos]
        except (OSError, zipfile.BadZipFile) as error:
            raise HLT14BundleError("NPZ ZIP structure is unreadable") from error
        expected = [f"{key}.npy" for key in keys]
        byte_limits = {"metadata_utf8.npy": _MAX_METADATA_BYTES}
        try:
            for descriptor in array_descriptors:
                shape = descriptor["shape"]
                count = 1
                for dimension in shape:
                    count *= int(dimension)
                # NPY headers are tiny relative to the arrays.  The fixed
                # allowance remains generous while preventing a declared
                # small array from expanding into an unbounded allocation.
                byte_limits[f"{descriptor['npz_storage_key']}.npy"] = count * 8 + 65536
        except (KeyError, TypeError, ValueError) as error:
            raise HLT14BundleError("array size contract is malformed") from error
        if (len(names) != len(set(names)) or set(names) != set(expected)
                or any(
                    "/" in info.filename
                    or "\\" in info.filename
                    or info.filename.startswith(".")
                    or info.flag_bits & 1
                    or info.file_size > byte_limits.get(info.filename, 0)
                    for info in infos
                )):
            raise HLT14BundleError("NPZ ZIP member names differ")


__all__ = ["HLT14AuthorityError", "HLT14BundleError", "HLT14InputError", "HLT14NamespaceForbiddenError", "LaunchAuthorityTuple", "LaunchAuthorityVerifier", "VerifiedAuthority", "VerifiedBundle", "VerifiedBundleAdapter", "canonical_bytes", "canonical_sha256", "decode_canonical_json", "safe_relative_path", "shape_framed_array_content_sha256"]
