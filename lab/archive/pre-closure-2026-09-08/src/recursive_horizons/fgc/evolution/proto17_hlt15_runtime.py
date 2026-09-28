"""Generation-zero-only production store materializer for HLT15.

Unlike the synthetic HLT14 runtime this class has no accepted-record replay,
no recovery/resume operation, and no common-event commit API. It materializes
only an externally authenticated AUTH1 generation-zero checkpoint and receipt,
then independently reopens that immutable generation-zero store.
"""
from __future__ import annotations

from copy import deepcopy
import ctypes
import errno
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import stat
import tempfile
from typing import Any, Mapping

from .proto17_pure_construction import (
    Proto17ConstructionError,
    build_genesis,
    canonical,
    digest,
    validate_checkpoint,
)


class Proto17HLT15Error(ValueError):
    """HLT15 generation-zero materialization is invalid or unsafe."""


_AUTHORITY_FIELDS = frozenset({
    "artifact_id", "authorization_commit", "authorization_result_path",
    "authorization_result_sha256", "genesis_spec_sha256",
    "source_closure_sha256", "numerical_environment_sha256",
    "pref26_auth1_input_evidence_sha256",
})


def _hex(value: object, length: int, label: str) -> str:
    if not isinstance(value, str) or len(value) != length or value.lower() != value:
        raise Proto17HLT15Error(f"{label} identity differs")
    try:
        int(value, 16)
    except ValueError as exc:
        raise Proto17HLT15Error(f"{label} identity differs") from exc
    return value


def _safe_relative(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise Proto17HLT15Error(f"{label} path differs")
    path = Path(value)
    if path.is_absolute() or "\\" in value or any(part in {"", ".", ".."} for part in path.parts):
        raise Proto17HLT15Error(f"{label} path differs")
    return value


def _atomic(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("xb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        descriptor = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    finally:
        temporary.unlink(missing_ok=True)


def _canonical_object(path: Path, label: str) -> dict[str, Any]:
    try:
        raw = path.read_bytes()
        value = json.loads(raw.decode("utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise Proto17HLT15Error(f"{label} is unreadable") from exc
    if not isinstance(value, dict) or canonical(value) != raw:
        raise Proto17HLT15Error(f"{label} is noncanonical")
    return value


def _lstat_absent(path: Path) -> None:
    try:
        os.lstat(path)
    except FileNotFoundError:
        return
    except OSError as exc:
        raise Proto17HLT15Error("HLT15 namespace precondition cannot be inspected") from exc
    raise Proto17HLT15Error("HLT15 refuses overwrite or resume of a namespace")


def _rename_exclusive(source: Path, destination: Path) -> None:
    """Atomically adopt a directory without replacing an arriving target.

    AUTH1 pins Darwin/arm64. Darwin's ``renamex_np(RENAME_EXCL)`` is the
    required no-clobber primitive; failing to expose it is an environment stop,
    never permission to fall back to a racy check-then-rename operation.
    """
    if os.uname().sysname != "Darwin":
        raise Proto17HLT15Error("exclusive namespace adoption is unavailable")
    libc = ctypes.CDLL(None, use_errno=True)
    renamex = getattr(libc, "renamex_np", None)
    if renamex is None:
        raise Proto17HLT15Error("exclusive namespace adoption is unavailable")
    renamex.argtypes = (ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint)
    renamex.restype = ctypes.c_int
    if renamex(os.fsencode(source), os.fsencode(destination), 0x00000004) != 0:
        error = ctypes.get_errno()
        if error in {errno.EEXIST, errno.ENOTEMPTY}:
            raise Proto17HLT15Error("HLT15 namespace appeared during staging")
        raise Proto17HLT15Error(
            f"exclusive namespace adoption failed with errno {error}"
        )


def _safe_parent_chain(repository_root: Path, namespace: str) -> Path:
    root = Path(repository_root)
    before = root.lstat()
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
        raise Proto17HLT15Error("repository root is unsafe")
    current = root
    parts = Path(_safe_relative(namespace, "PROTO17 namespace")).parts
    for part in parts[:-1]:
        candidate = current / part
        try:
            observed = candidate.lstat()
        except FileNotFoundError:
            try:
                candidate.mkdir(mode=0o755)
                observed = candidate.lstat()
            except OSError as exc:
                raise Proto17HLT15Error("namespace parent cannot be created safely") from exc
        if stat.S_ISLNK(observed.st_mode) or not stat.S_ISDIR(observed.st_mode):
            raise Proto17HLT15Error("namespace parent is symlinked or non-directory")
        current = candidate
    return root.joinpath(*parts)


class Proto17HLT15Runtime:
    """Atomic, no-overwrite production GEN1 materializer.

    The repository root is supplied by the already authenticated runner. The
    output root is derived only from the GenesisSpec; callers cannot override
    it. Construction performs no output-root observation.
    """

    def __init__(
        self,
        repository_root: Path,
        genesis_spec: Mapping[str, object],
        *,
        authority_binding: Mapping[str, object],
    ) -> None:
        try:
            derived = build_genesis(genesis_spec)
        except (Proto17ConstructionError, TypeError, ValueError) as exc:
            raise Proto17HLT15Error("AUTH1 GenesisSpec does not derive") from exc
        if not isinstance(authority_binding, Mapping) or set(authority_binding) != _AUTHORITY_FIELDS:
            raise Proto17HLT15Error("external AUTH1 binding fields differ")
        binding = deepcopy(dict(authority_binding))
        if (binding["artifact_id"] != "FGC-1-PRO18-AUTH1"
                or _hex(binding["authorization_commit"], 40, "AUTH1 commit")
                   != binding["authorization_commit"]
                or _safe_relative(binding["authorization_result_path"], "AUTH1 result")
                   != binding["authorization_result_path"]):
            raise Proto17HLT15Error("external AUTH1 binding identity differs")
        for name in _AUTHORITY_FIELDS - {
            "artifact_id", "authorization_commit", "authorization_result_path",
        }:
            _hex(binding[name], 64, name)
        if binding["genesis_spec_sha256"] != digest(derived["genesis_spec"]):
            raise Proto17HLT15Error("external AUTH1 GenesisSpec binding differs")
        self.repository_root = Path(repository_root)
        self.spec = deepcopy(derived["genesis_spec"])
        self.genesis = deepcopy(derived["checkpoint"])
        self.store_plan = deepcopy(derived["store_plan"])
        self.authority_binding = binding

    @staticmethod
    def _state_path(root: Path, state_hash: str) -> Path:
        return root / "states" / f"{state_hash}.json"

    @staticmethod
    def _checkpoint_path(root: Path, checkpoint: Mapping[str, Any]) -> Path:
        return root / "checkpoints" / f"00000000000000000000-{checkpoint['checkpoint_sha256']}.json"

    def _receipt(self, checkpoint: Mapping[str, Any]) -> dict[str, Any]:
        bare = {
            "artifact_id": "FGC-1-HLT15-GEN1",
            "kind": "generation_zero_materialized",
            "campaign_id": checkpoint["campaign_id"],
            "campaign_generation": 0,
            "checkpoint_sha256": checkpoint["checkpoint_sha256"],
            "state_object_set_sha256": checkpoint["accepted_state_set_sha256"],
            "namespace_relative_path": self.spec["namespace"],
            "external_authority": deepcopy(self.authority_binding),
            "claims": {
                "generation_zero_materialized": True,
                "state_advanced": False,
                "calibration_completed": False,
                "candidate_execution_authorized": False,
                "physical_transition_claim_authorized": False,
            },
        }
        return {**bare, "receipt_sha256": digest(bare)}

    def materialize_generation_zero(self) -> Mapping[str, Any]:
        root = _safe_parent_chain(self.repository_root, self.spec["namespace"])
        # This is the explicit HLT15 output-root boundary. Every authority,
        # environment, import, and raw-reimport check must already have passed.
        _lstat_absent(root)
        staging = Path(tempfile.mkdtemp(prefix=f".{root.name}.hlt15-stage-", dir=root.parent))
        adopted = False
        try:
            for state_hash, state in self.store_plan["state_objects"].items():
                payload = canonical(state)
                if sha256(payload).hexdigest() != state_hash:
                    raise Proto17HLT15Error("derived state object content address differs")
                _atomic(self._state_path(staging, state_hash), payload)
            _atomic(self._checkpoint_path(staging, self.genesis), canonical(self.genesis))
            receipt = self._receipt(self.genesis)
            _atomic(staging / "receipts" / "generation-zero.json", canonical(receipt))
            descriptor = os.open(staging, os.O_RDONLY)
            try:
                os.fsync(descriptor)
            finally:
                os.close(descriptor)
            _rename_exclusive(staging, root)
            adopted = True
            parent = os.open(root.parent, os.O_RDONLY)
            try:
                os.fsync(parent)
            finally:
                os.close(parent)
            return self.reopen_verify()
        except BaseException:
            # Before adoption, remove only this invocation's private staging
            # directory. After adoption preserve the failed store verbatim as
            # forensic evidence; HLT15 has no overwrite or resume path.
            if not adopted and staging.exists():
                shutil.rmtree(staging)
            raise

    def reopen_verify(self) -> Mapping[str, Any]:
        root = self.repository_root / self.spec["namespace"]
        checkpoint_path = self._checkpoint_path(root, self.genesis)
        observed = _canonical_object(checkpoint_path, "generation-zero checkpoint")
        try:
            if validate_checkpoint(observed) != self.genesis:
                raise Proto17HLT15Error("generation-zero checkpoint differs from AUTH1 derivation")
        except Proto17ConstructionError as exc:
            raise Proto17HLT15Error("generation-zero checkpoint construction differs") from exc
        for state_hash, state in self.store_plan["state_objects"].items():
            path = self._state_path(root, state_hash)
            if _canonical_object(path, "generation-zero state") != state:
                raise Proto17HLT15Error("generation-zero state differs from AUTH1 derivation")
        receipt = _canonical_object(root / "receipts" / "generation-zero.json", "generation-zero receipt")
        bare = dict(receipt)
        observed_digest = bare.pop("receipt_sha256", None)
        if observed_digest != digest(bare) or receipt != self._receipt(self.genesis):
            raise Proto17HLT15Error("generation-zero receipt differs from materialized store")
        return {"checkpoint": deepcopy(self.genesis), "receipt": deepcopy(receipt)}


__all__ = ["Proto17HLT15Error", "Proto17HLT15Runtime"]
