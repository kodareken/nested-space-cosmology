"""Pure contract for the first authenticated PROTO18 GR-0 event.

The constructor reads compact, tracked evidence only. It has no array, run
namespace, numerical-engine, or branch-selection capability. The external
authorization commit is supplied by the eventual runner so the plan can bind
the immutable HLT16 prelaunch checkpoint without a self-hash.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
import os
from pathlib import Path
import stat
from typing import Any, Iterable, Mapping


MEMBERS = (
    "RK4-2049",
    "RK4-4097",
    "RK4-8193",
    "SSPRK3-4097",
    "SSPRK3-8193",
    "SSPRK3-16385",
)
FIRST_EDGE_PROTOCOL = "FGC-2-SF1-PROTO18"
FIRST_EDGE_BRANCH = "GR-0"
FIRST_EDGE_AMPLITUDE = "3"
STATE_SCHEMA_ID = "FGC-1-HLT16-member-state-v1"


class Proto19ProgressionContractError(ValueError):
    """Compact evidence cannot construct the frozen first-event plan."""


def _canonical(value: object) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")


def _digest(value: object) -> str:
    return sha256(_canonical(value)).hexdigest()


def _sha256(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64:
        raise Proto19ProgressionContractError(f"{label} is not SHA-256")
    try:
        int(value, 16)
    except ValueError as error:
        raise Proto19ProgressionContractError(f"{label} is not SHA-256") from error
    if value.lower() != value:
        raise Proto19ProgressionContractError(f"{label} is not SHA-256")
    return value


def _commit(value: object) -> str:
    if not isinstance(value, str) or len(value) != 40:
        raise Proto19ProgressionContractError(
            "authorization commit is not a Git object id"
        )
    try:
        int(value, 16)
    except ValueError as error:
        raise Proto19ProgressionContractError(
            "authorization commit is not a Git object id"
        ) from error
    if value.lower() != value:
        raise Proto19ProgressionContractError(
            "authorization commit is not a Git object id"
        )
    return value


@dataclass(frozen=True, slots=True)
class TimeIdentity:
    rational: str
    binary64_hex: str

    def __post_init__(self) -> None:
        try:
            value = Fraction(self.rational)
        except (TypeError, ValueError, ZeroDivisionError) as error:
            raise Proto19ProgressionContractError("invalid rational time") from error
        if str(value) != self.rational or float(value).hex() != self.binary64_hex:
            raise Proto19ProgressionContractError("time identity differs")


@dataclass(frozen=True, slots=True)
class MemberPlan:
    member_key: str
    source_state_object_sha256: str
    method: str
    point_count: int

    def __post_init__(self) -> None:
        if (
            self.member_key not in MEMBERS
            or self.member_key != f"{self.method}-{self.point_count}"
            or isinstance(self.point_count, bool)
            or self.point_count <= 0
        ):
            raise Proto19ProgressionContractError("member plan differs")
        _sha256(self.source_state_object_sha256, "source state object")


@dataclass(frozen=True, slots=True)
class StopPolicy:
    source_retry_cap: int
    cfl_retry_cap: int
    temporal_retry_cap: int
    scientific_stops: tuple[str, ...]

    def __post_init__(self) -> None:
        required = (
            "branch_loss",
            "kinetic_loss",
            "hyperbolicity_loss",
            "health_stop",
            "constraint_stop",
            "scale_stop",
            "aliasing_stop",
            "boundary_stop",
            "mass_flux_inconsistency",
        )
        if (
            (self.source_retry_cap, self.cfl_retry_cap, self.temporal_retry_cap)
            != (32, 32, 32)
            or self.scientific_stops != required
        ):
            raise Proto19ProgressionContractError("stop policy differs")


@dataclass(frozen=True, slots=True)
class ProgressionPlan:
    authorization_commit: str
    protocol_artifact_id: str
    campaign_id: str
    source_receipt_sha256: str
    source_checkpoint_sha256: str
    tdg6_contract_sha256: str
    tdg7_contract_sha256: str
    state_schema_id: str
    start_time: TimeIdentity
    target_time: TimeIdentity
    common_event_index: int
    branch: str
    amplitude: str
    members: tuple[MemberPlan, ...]
    stop_policy: StopPolicy

    def __post_init__(self) -> None:
        _commit(self.authorization_commit)
        for label, value in (
            ("source receipt", self.source_receipt_sha256),
            ("source checkpoint", self.source_checkpoint_sha256),
            ("TDG6 contract", self.tdg6_contract_sha256),
            ("TDG7 contract", self.tdg7_contract_sha256),
        ):
            _sha256(value, label)
        if (
            self.protocol_artifact_id != FIRST_EDGE_PROTOCOL
            or self.branch != FIRST_EDGE_BRANCH
            or self.amplitude != FIRST_EDGE_AMPLITUDE
            or self.state_schema_id != STATE_SCHEMA_ID
            or self.common_event_index != 23
            or tuple(member.member_key for member in self.members) != MEMBERS
            or self.start_time != TimeIdentity("23/16", "0x1.7000000000000p+0")
            or self.target_time != TimeIdentity("3/2", "0x1.8000000000000p+0")
        ):
            raise Proto19ProgressionContractError("first progression plan differs")

    def mapping(self) -> dict[str, Any]:
        return asdict(self)

    @property
    def sha256(self) -> str:
        return _digest(self.mapping())


def _safe_relative(value: str) -> tuple[str, ...]:
    path = Path(value)
    if (
        path.is_absolute()
        or "\\" in value
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise Proto19ProgressionContractError("compact evidence path is unsafe")
    return path.parts


def _read_nofollow(root: Path, relative: str, label: str) -> bytes:
    """Read one unique regular file through an all-no-follow descriptor chain."""
    parts = _safe_relative(relative)
    directory_flags = (
        os.O_RDONLY
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    leaf_flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    directories: list[int] = []
    leaf = -1
    try:
        parent = os.open(root, directory_flags)
        directories.append(parent)
        for component in parts[:-1]:
            before = os.stat(component, dir_fd=parent, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                raise Proto19ProgressionContractError(f"{label} parent is unsafe")
            child = os.open(component, directory_flags, dir_fd=parent)
            after = os.fstat(child)
            if (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino):
                os.close(child)
                raise Proto19ProgressionContractError(f"{label} parent raced")
            directories.append(child)
            parent = child
        before = os.stat(parts[-1], dir_fd=parent, follow_symlinks=False)
        if (
            stat.S_ISLNK(before.st_mode)
            or not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
        ):
            raise Proto19ProgressionContractError(f"{label} is unsafe")
        leaf = os.open(parts[-1], leaf_flags, dir_fd=parent)
        active = os.fstat(leaf)
        identity = (
            before.st_dev,
            before.st_ino,
            before.st_size,
            before.st_mtime_ns,
            before.st_ctime_ns,
        )
        if identity != (
            active.st_dev,
            active.st_ino,
            active.st_size,
            active.st_mtime_ns,
            active.st_ctime_ns,
        ):
            raise Proto19ProgressionContractError(f"{label} raced")
        chunks: list[bytes] = []
        while True:
            block = os.read(leaf, 1 << 20)
            if not block:
                break
            chunks.append(block)
        raw = b"".join(chunks)
        after = os.fstat(leaf)
        if len(raw) != before.st_size or identity != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
            after.st_ctime_ns,
        ):
            raise Proto19ProgressionContractError(f"{label} changed during read")
        return raw
    except Proto19ProgressionContractError:
        raise
    except OSError as error:
        raise Proto19ProgressionContractError(f"{label} is unavailable") from error
    finally:
        if leaf != -1:
            os.close(leaf)
        for descriptor in reversed(directories):
            os.close(descriptor)


def _load_object(raw: bytes, label: str) -> Mapping[str, Any]:
    def reject_duplicates(items: Iterable[tuple[str, object]]) -> dict[str, object]:
        answer: dict[str, object] = {}
        for key, value in items:
            if key in answer:
                raise ValueError(key)
            answer[key] = value
        return answer

    try:
        value = json.loads(
            raw.decode("utf-8"), object_pairs_hook=reject_duplicates
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise Proto19ProgressionContractError(f"{label} is malformed") from error
    if not isinstance(value, Mapping):
        raise Proto19ProgressionContractError(f"{label} is not object-valued")
    return value


def construct_first_event(
    root: Path,
    *,
    authorization_commit: str,
    evidence_bytes: Mapping[str, bytes] | None = None,
) -> ProgressionPlan:
    """Construct one event from exactly the compact bytes already authenticated.

    ``evidence_bytes`` is the launch authority's race-free handoff.  The
    direct path remains useful for focused inspection, but it uses the same
    all-no-follow reader and therefore cannot silently substitute a symlink.
    """

    repository = Path(root)
    paths = (
        "results/fgc-1-pro19-frz1.json",
        "results/fgc-1-pro18-auth1.json",
    )
    if evidence_bytes is None:
        captured = {
            path: _read_nofollow(repository, path, Path(path).stem)
            for path in paths
        }
    else:
        if set(evidence_bytes) != set(paths) or any(
            not isinstance(value, bytes) for value in evidence_bytes.values()
        ):
            raise Proto19ProgressionContractError(
                "compact evidence byte inventory differs"
            )
        captured = dict(evidence_bytes)
    freeze = _load_object(captured[paths[0]], "PRO19")
    authority = _load_object(captured[paths[1]], "AUTH1")
    if (
        freeze.get("artifact_id") != "FGC-1-PRO19-FRZ1"
        or authority.get("artifact_id") != "FGC-1-PRO18-AUTH1"
    ):
        raise Proto19ProgressionContractError("compact artifact identity differs")
    artifact_payload = freeze.get("artifact_payload")
    auth_payload = authority.get("artifact_payload")
    if not isinstance(artifact_payload, Mapping) or not isinstance(auth_payload, Mapping):
        raise Proto19ProgressionContractError("compact payload is absent")
    bindings = artifact_payload.get("bindings")
    if not isinstance(bindings, list) or any(
        not isinstance(item, Mapping) for item in bindings
    ):
        raise Proto19ProgressionContractError("PRO19 binding inventory differs")
    binding_map = {
        str(item.get("path")): str(item.get("sha256"))
        for item in bindings
    }
    if (
        len(binding_map) != len(bindings)
        or binding_map.get("results/fgc-1-tdg6-imp2.json")
        != "6a2cac24f65e067aed7ff7de299f9c1f0f3c3319e9aee692763b8df1db4fd5b5"
        or binding_map.get("results/fgc-1-tdg7-imp3.json")
        != "364223f0f0f93c446b9d5b6cf09d3aa3cbe5b48bd1c175faf96a34459fd5b779"
        or binding_map.get("results/fgc-1-pro18-pref27.json")
        != "102b385d058066d668e9e8541c974c736c583f6d5c4afaf9b99e95a16aabeb83"
    ):
        raise Proto19ProgressionContractError("PRO19 predecessor binding differs")
    edge = artifact_payload.get("first_edge")
    if not isinstance(edge, Mapping) or (
        edge.get("members") != list(MEMBERS)
        or edge.get("branch") != FIRST_EDGE_BRANCH
        or edge.get("amplitude") != FIRST_EDGE_AMPLITUDE
        or edge.get("start_time") != "23/16"
        or edge.get("target_time") != "3/2"
        or edge.get("common_event_index") != 23
        or edge.get("event_count") != 1
        or edge.get("raw_import_permitted") is not False
        or edge.get("run_permitted") is not False
    ):
        raise Proto19ProgressionContractError("frozen edge differs")
    genesis = auth_payload.get("genesis_spec")
    if not isinstance(genesis, Mapping):
        raise Proto19ProgressionContractError("AUTH1 genesis is absent")
    descriptors = genesis.get("member_descriptors")
    if (
        not isinstance(descriptors, list)
        or [item.get("member_key") for item in descriptors if isinstance(item, Mapping)]
        != list(MEMBERS)
        or len(descriptors) != len(MEMBERS)
    ):
        raise Proto19ProgressionContractError("descriptor inventory differs")
    members: list[MemberPlan] = []
    for descriptor in descriptors:
        if not isinstance(descriptor, Mapping):
            raise Proto19ProgressionContractError("member descriptor differs")
        key = str(descriptor["member_key"])
        method, point_count = key.split("-", 1)
        members.append(
            MemberPlan(
                member_key=key,
                source_state_object_sha256=str(
                    descriptor["state_object_canonical_sha256"]
                ),
                method=method,
                point_count=int(point_count),
            )
        )
    return ProgressionPlan(
        authorization_commit=_commit(authorization_commit),
        protocol_artifact_id=FIRST_EDGE_PROTOCOL,
        campaign_id=str(genesis["campaign_id"]),
        source_receipt_sha256=(
            "78ae89a9f74193658003cf51f8567c4088db80fd6fb86867c2d581a0c16d5baf"
        ),
        source_checkpoint_sha256=(
            "7c059dcc2197a9f57baf8012f4113c240917090171ee06cc639e9e8734603d2d"
        ),
        tdg6_contract_sha256=(
            "6a2cac24f65e067aed7ff7de299f9c1f0f3c3319e9aee692763b8df1db4fd5b5"
        ),
        tdg7_contract_sha256=(
            "364223f0f0f93c446b9d5b6cf09d3aa3cbe5b48bd1c175faf96a34459fd5b779"
        ),
        state_schema_id=STATE_SCHEMA_ID,
        start_time=TimeIdentity("23/16", "0x1.7000000000000p+0"),
        target_time=TimeIdentity("3/2", "0x1.8000000000000p+0"),
        common_event_index=23,
        branch=FIRST_EDGE_BRANCH,
        amplitude=FIRST_EDGE_AMPLITUDE,
        members=tuple(members),
        stop_policy=StopPolicy(
            source_retry_cap=32,
            cfl_retry_cap=32,
            temporal_retry_cap=32,
            scientific_stops=(
                "branch_loss",
                "kinetic_loss",
                "hyperbolicity_loss",
                "health_stop",
                "constraint_stop",
                "scale_stop",
                "aliasing_stop",
                "boundary_stop",
                "mass_flux_inconsistency",
            ),
        ),
    )


__all__ = [
    "FIRST_EDGE_AMPLITUDE",
    "FIRST_EDGE_BRANCH",
    "FIRST_EDGE_PROTOCOL",
    "MEMBERS",
    "STATE_SCHEMA_ID",
    "MemberPlan",
    "ProgressionPlan",
    "Proto19ProgressionContractError",
    "StopPolicy",
    "TimeIdentity",
    "construct_first_event",
]
