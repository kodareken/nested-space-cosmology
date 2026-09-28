"""Independent post-result binder for the completed TDG9 LOC2 diagnostic.

The live construction path authenticates the sealed authority and two raw
leaves, reconstructs the ten selected finite-dimensional occurrences from the
read-only campaign store, and independently applies the primary and v2 exact
localizers.  The ordinary compact path validates only tracked evidence and
never opens the raw namespace or campaign store.  This module deliberately
does not import the LOC2 runner.
"""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import stat
import subprocess
import tomllib
from typing import Any, Iterable, Iterator, Mapping, NoReturn, Sequence

import numpy as np

from . import tdg9_ar1_authority as ar1
from . import tdg9_ar1_pref1_binder as ar1_pref1
from . import tdg9_loc1_authority as loc1_authority
from . import tdg9_loc2_authority as loc2_authority
from . import tdg9_local_extrema as primary
from . import tdg9_local_extrema_independent_v2 as independent_v2
from .hlt16_campaign_store import HLT16CampaignStore
from .proto19_gr0_static_factory import build_static_gr0_shells


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG9-LOC2-PREF2"
CLASSIFICATION = "independently_bound_scale_invariant_localization_terminal"
CONFIG_PATH = "configs/fgc/fgc-1-tdg9-loc2-pref2.toml"
RESULT_PATH = "results/fgc-1-tdg9-loc2-pref2.json"
OWNER_DOCUMENT = "docs/fgc-tdg9-loc2-pref2.md"
TARGET_PROTOCOL = loc2_authority.TARGET_PROTOCOL

AUTHORITY_COMMIT = "1cd78ee20100ae4ae6495161146c6f2fcf7de1bf"
AUTHORITY_PARENT = loc2_authority.BASE_COMMIT
RAW_NAMESPACE = loc2_authority.OUTPUT_NAMESPACE
RAW_SCHEMA = "FGC-1-TDG9-LOC2-raw-v1"
RAW_CLASSIFICATION = "bounded_exact_localization_completed"
RUNNER_ID = "FGC-1-TDG9-LOC2-RUN1"
RAW_MANIFEST_SHA256 = "8189b030180444fdeaa8c7df7238550316ae0613efb7f1f8a1466d49b327858f"
RAW_TERMINAL_SHA256 = "571bb2979bd73229fbefe07d5a2d137a8e68b68661e94b7dce67f63ac8f04232"
COMPACT_RESULT_SHA256 = "2f1630f449aa4d942cd0b0e45cf6ebb0f43862f2e387da51172584680aa03232"

SEALED_STORE_LEAF_COUNT = 115
SEALED_STORE_SHA256 = "5917d70de3080dcc6b7abeb7fdb7cc3370c6140cbeb37bde0c142d6a2d36a445"
EXPECTED_COMPONENT_CANDIDATES = {
    "complete_C": 326_164,
    "slope_S": 391_791,
    "value_V": 245_280,
}
EXPECTED_TOTAL_CANDIDATES = 963_235
EXPECTED_ROOT_METHODS = {
    "adaptive_nonsquare": 21,
    "exact_endpoint": 3_013,
    "exact_endpoint_deflation": 28,
}
EXPECTED_ROUTE_CLASSIFICATIONS = {
    "nonunique_or_interval_inconclusive": 31,
    "unique_maximum": 29,
}
EXPECTED_SURVIVOR_COUNT = 3_062
_BLOCKS = ("u", "p", "q")
_FIELDS = ("alpha", "v", "lambda", "R", "phi", "chi")
_MAX_RAW_LEAF_BYTES = 64 * 1024 * 1024
_MAX_STORE_LEAF_BYTES = 128 * 1024 * 1024

_AUTHORITY_BLOBS = (
    (
        loc2_authority.CONFIG_PATH,
        "b22b54584d07894d7fe70badce012c8fc7f9d75d50f16fd05ad97d507c15451d",
    ),
    (
        loc2_authority.RESULT_PATH,
        "62e15c731b53054d32a6bcb42dcbf95c0b879a1ce559edacbf1571e69fe22718",
    ),
    (
        loc2_authority.OWNER_DOCUMENT,
        "5ea8ed1ccf9c410de24ccb78091f8c53a118999e93ac1899f2805e7bf58d4a08",
    ),
    (
        "src/recursive_horizons/fgc/evolution/tdg9_loc2_authority.py",
        "2da55fb0126091aaad501c1b7d63759dc62e47615bd72ea898f65099a9bcbdf7",
    ),
    (loc2_authority.EVALUATOR_PATH, loc2_authority.EVALUATOR_SHA256),
    (loc2_authority.RUNNER_PATH, loc2_authority.RUNNER_SHA256),
)


class TDG9LOC2PREF2Error(ValueError):
    """The compact contract, raw terminal, store, or independent replay differs."""

    def __init__(self, stop_id: str, detail: object) -> None:
        self.stop_id = str(stop_id)
        self.detail = " ".join(str(detail).split())[:640]
        super().__init__(f"{self.stop_id}: {self.detail}")


def _stop(stop_id: str, detail: object) -> NoReturn:
    raise TDG9LOC2PREF2Error(stop_id, detail)


def _canonical(value: object) -> bytes:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("ascii")
    except (TypeError, ValueError) as exc:
        raise TDG9LOC2PREF2Error("PREF2_CANONICAL_DRIFT", exc) from exc


def canonical_result(value: object) -> bytes:
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
        raise TDG9LOC2PREF2Error("PREF2_COMPACT_DRIFT", exc) from exc


def _pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in items:
        if key in answer:
            raise ValueError(f"duplicate key {key}")
        answer[key] = value
    return answer


def _json(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(
            raw.decode("ascii"),
            object_pairs_hook=_pairs,
            parse_constant=lambda item: (_ for _ in ()).throw(ValueError(item)),
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TDG9LOC2PREF2Error("PREF2_JSON_DRIFT", label) from exc
    if not isinstance(value, dict) or canonical_result(value) != raw:
        _stop("PREF2_JSON_DRIFT", f"{label} is not canonical pretty JSON")
    return value


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
        raise TDG9LOC2PREF2Error("PREF2_GIT_DRIFT", arguments) from exc


def _identity(metadata: os.stat_result) -> tuple[int, int, int, int]:
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
    )


def _open_directory(root: Path, relative: str) -> int:
    candidate = Path(relative)
    if candidate.is_absolute() or any(
        part in {"", ".", ".."} for part in candidate.parts
    ):
        _stop("PREF2_PATH_UNSAFE", relative)
    descriptor = os.open(
        root,
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        for part in candidate.parts:
            child = os.open(
                part,
                os.O_RDONLY
                | getattr(os, "O_DIRECTORY", 0)
                | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=descriptor,
            )
            metadata = os.fstat(child)
            if not stat.S_ISDIR(metadata.st_mode):
                os.close(child)
                _stop("PREF2_PATH_UNSAFE", relative)
            os.close(descriptor)
            descriptor = child
        return descriptor
    except Exception:
        os.close(descriptor)
        raise


def _read_leaf_at(
    directory_fd: int,
    name: str,
    *,
    maximum: int,
    label: str,
) -> bytes:
    if not name or "/" in name or name in {".", ".."}:
        _stop("PREF2_PATH_UNSAFE", label)
    before = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        _stop("PREF2_LEAF_UNSAFE", label)
    if before.st_size > maximum:
        _stop("PREF2_LEAF_SIZE", label)
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=directory_fd,
    )
    try:
        active = os.fstat(descriptor)
        if _identity(before) != _identity(active):
            _stop("PREF2_LEAF_RACED", label)
        chunks: list[bytes] = []
        total = 0
        while total <= maximum:
            block = os.read(descriptor, min(1 << 20, maximum + 1 - total))
            if not block:
                break
            chunks.append(block)
            total += len(block)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    final = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    if (
        _identity(before) != _identity(after)
        or _identity(before) != _identity(final)
        or total != before.st_size
        or total > maximum
    ):
        _stop("PREF2_LEAF_CHANGED", label)
    return b"".join(chunks)


def _read_repo_leaf(root: Path, relative: str, maximum: int) -> bytes:
    candidate = Path(relative)
    parent = candidate.parent.as_posix()
    if parent == ".":
        descriptor = os.open(
            root,
            os.O_RDONLY
            | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_NOFOLLOW", 0),
        )
    else:
        descriptor = _open_directory(root, parent)
    try:
        return _read_leaf_at(
            descriptor,
            candidate.name,
            maximum=maximum,
            label=relative,
        )
    finally:
        os.close(descriptor)


def _raw_pair(root: Path) -> tuple[bytes, bytes]:
    descriptor = _open_directory(root, RAW_NAMESPACE)
    try:
        before = os.fstat(descriptor)
        with os.scandir(descriptor) as entries:
            names = tuple(sorted(item.name for item in entries))
        if names != ("manifest.json", "terminal.json"):
            _stop("PREF2_RAW_TREE_DRIFT", names)
        manifest_raw = _read_leaf_at(
            descriptor,
            "manifest.json",
            maximum=1 << 20,
            label=f"{RAW_NAMESPACE}/manifest.json",
        )
        terminal_raw = _read_leaf_at(
            descriptor,
            "terminal.json",
            maximum=_MAX_RAW_LEAF_BYTES,
            label=f"{RAW_NAMESPACE}/terminal.json",
        )
        with os.scandir(descriptor) as entries:
            final_names = tuple(sorted(item.name for item in entries))
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if names != final_names or _identity(before) != _identity(after):
        _stop("PREF2_RAW_TREE_CHANGED", RAW_NAMESPACE)
    if (
        sha256(manifest_raw).hexdigest() != RAW_MANIFEST_SHA256
        or sha256(terminal_raw).hexdigest() != RAW_TERMINAL_SHA256
    ):
        _stop("PREF2_RAW_HASH_DRIFT", RAW_NAMESPACE)
    return manifest_raw, terminal_raw


def _raw_snapshot(
    root: Path,
) -> tuple[bytes, bytes, dict[str, Any], dict[str, Any]]:
    manifest_raw, terminal_raw = _raw_pair(root)
    return (
        manifest_raw,
        terminal_raw,
        _json(manifest_raw, "raw manifest"),
        _json(terminal_raw, "raw terminal"),
    )


def _snapshot_store(root: Path) -> tuple[int, str]:
    """Reproduce the sealed stack-order store digest with no-following reads."""

    base_fd = _open_directory(root, ar1.PREF2_STORE_PATH)
    digest = sha256(b"TDG9-LOC1-STORE-SNAPSHOT-v1\n")
    count = 0
    stack: list[tuple[int, str]] = [(base_fd, "")]
    try:
        while stack:
            directory_fd, prefix = stack.pop()
            try:
                before = os.fstat(directory_fd)
                with os.scandir(directory_fd) as entries:
                    items = sorted(entries, key=lambda item: item.name)
                initial_names = tuple(item.name for item in items)
                children: list[tuple[int, str]] = []
                for item in items:
                    metadata = item.stat(follow_symlinks=False)
                    relative = f"{prefix}/{item.name}" if prefix else item.name
                    if item.is_symlink():
                        _stop("PREF2_STORE_SYMLINK", relative)
                    if stat.S_ISDIR(metadata.st_mode):
                        child_fd = os.open(
                            item.name,
                            os.O_RDONLY
                            | getattr(os, "O_DIRECTORY", 0)
                            | getattr(os, "O_NOFOLLOW", 0),
                            dir_fd=directory_fd,
                        )
                        active = os.fstat(child_fd)
                        if _identity(metadata) != _identity(active):
                            os.close(child_fd)
                            _stop("PREF2_STORE_DIRECTORY_RACED", relative)
                        children.append((child_fd, relative))
                    elif stat.S_ISREG(metadata.st_mode):
                        raw = _read_leaf_at(
                            directory_fd,
                            item.name,
                            maximum=_MAX_STORE_LEAF_BYTES,
                            label=relative,
                        )
                        digest.update(
                            (
                                relative
                                + "\0"
                                + sha256(raw).hexdigest()
                                + "\n"
                            ).encode("ascii")
                        )
                        count += 1
                    else:
                        _stop("PREF2_STORE_FOREIGN_TYPE", relative)
                with os.scandir(directory_fd) as entries:
                    final_names = tuple(sorted(item.name for item in entries))
                after = os.fstat(directory_fd)
                if initial_names != final_names or _identity(before) != _identity(after):
                    for child_fd, _ in children:
                        os.close(child_fd)
                    _stop("PREF2_STORE_DIRECTORY_CHANGED", prefix or ".")
                stack.extend(children)
            finally:
                os.close(directory_fd)
    except Exception:
        for directory_fd, _ in stack:
            try:
                os.close(directory_fd)
            except OSError:
                pass
        raise
    return count, digest.hexdigest()


def _expected_config() -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "target_protocol": TARGET_PROTOCOL,
        "owner_document": OWNER_DOCUMENT,
        "raw": {
            "namespace": RAW_NAMESPACE,
            "schema": RAW_SCHEMA,
            "classification": RAW_CLASSIFICATION,
            "runner_id": RUNNER_ID,
            "authority_commit": AUTHORITY_COMMIT,
            "manifest_sha256": RAW_MANIFEST_SHA256,
            "terminal_sha256": RAW_TERMINAL_SHA256,
            "leaf_count": 2,
        },
        "authority": {
            "parent_commit": AUTHORITY_PARENT,
            "blob_count": len(_AUTHORITY_BLOBS),
            "blobs": [
                {"path": path, "sha256": digest}
                for path, digest in _AUTHORITY_BLOBS
            ],
        },
        "store": {
            "path": ar1.PREF2_STORE_PATH,
            "leaf_count": SEALED_STORE_LEAF_COUNT,
            "snapshot_sha256": SEALED_STORE_SHA256,
            "snapshot_algorithm": "TDG9-LOC1-STORE-SNAPSHOT-v1",
        },
        "selection": {
            "method": loc2_authority.METHOD,
            "member_key": loc2_authority.MEMBER_KEY,
            "point_count": loc2_authority.POINT_COUNT,
            "owned_row_count": loc2_authority.OWNED_ROW_COUNT,
            "grid_spacing_hex": loc2_authority.GRID_SPACING_HEX,
            "outer_radius_hex": loc2_authority.OUTER_RADIUS_HEX,
            "retries": [3, 4, 5],
            "failed_occurrences": [
                {"retry": retry, "channel": channel}
                for retry, channel in loc2_authority.FAILED_OCCURRENCES
            ],
        },
        "work_budget": {
            "base_cubic_count": loc2_authority.BASE_CUBIC_COUNT,
            "component_cubic_count": loc2_authority.COMPONENT_CUBIC_COUNT,
            "family_candidate_ceiling": loc2_authority.FAMILY_CANDIDATE_CEILING,
            "aggregate_candidate_ceiling": (
                loc2_authority.AGGREGATE_CANDIDATE_CEILING
            ),
            "primary_refinement_depth": loc2_authority.PRIMARY_REFINEMENT_DEPTH,
            "v2_initial_refinement_bits": loc2_authority.INITIAL_REFINEMENT_BITS,
            "v2_global_proof_bit_ceiling": (
                loc2_authority.GLOBAL_PROOF_BIT_CEILING
            ),
        },
        "scope": {
            "one_time_live_raw_and_store_reconstruction": True,
            "ordinary_verifier_raw_blind": True,
            "ordinary_verifier_store_blind": True,
            "LOC2_runner_decision_code_imported": False,
            "campaign_store_mutation_authorized": False,
            "raw_namespace_mutation_authorized": False,
            "successor_remedy_selected": False,
            "fourth_width_authorized": False,
            "stage2_source_decomposition_authorized": False,
            "SSPRK3_comparator_authorized": False,
            "PDE_state_commit_authorized": False,
            "continuation_authorized": False,
            "candidate_branches_authorized": False,
            "mechanism_result_authorized": False,
            "physical_result_authorized": False,
        },
        "claims": {
            "LOC2_raw_result_independently_bound": True,
            "LOC2_exact_localization_completed": True,
            "primary_v2_route_agreement_completed": True,
            "all_ten_occurrences_both_components_independently_fail": True,
            "temporal_instrument_remedy_selected": False,
            "common_event_completed": False,
            "GR0_calibration_completed": False,
            "candidate_execution_authorized": False,
            "mechanism_result_earned": False,
            "physical_result_earned": False,
            "retained_EFT_evolution_authorized": False,
            "physical_transition_claim_authorized": False,
        },
    }


def _config(raw: bytes) -> dict[str, Any]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise TDG9LOC2PREF2Error("PREF2_CONFIG_DRIFT", exc) from exc
    if value != _expected_config():
        _stop("PREF2_CONFIG_DRIFT", "typed PREF2 config differs")
    return value


def _authenticate_authority(root: Path, manifest: Mapping[str, Any]) -> None:
    if manifest.get("authority_commit") != AUTHORITY_COMMIT:
        _stop("PREF2_AUTHORITY_DRIFT", manifest.get("authority_commit"))
    if _git(root, "cat-file", "-t", AUTHORITY_COMMIT).strip() != b"commit":
        _stop("PREF2_AUTHORITY_DRIFT", "authority object is not a commit")
    parents = _git(root, "show", "-s", "--format=%P", AUTHORITY_COMMIT).split()
    if parents != [AUTHORITY_PARENT.encode("ascii")]:
        _stop("PREF2_AUTHORITY_LINEAGE_DRIFT", parents)
    for relative, expected in _AUTHORITY_BLOBS:
        observed = sha256(
            _git(root, "show", f"{AUTHORITY_COMMIT}:{relative}")
        ).hexdigest()
        if observed != expected:
            _stop("PREF2_AUTHORITY_BLOB_DRIFT", relative)


def _expected_manifest() -> dict[str, Any]:
    return {
        "ERR1_result_sha256": loc2_authority.ERR1_RESULT_SHA256,
        "artifact_id": loc2_authority.ARTIFACT_ID,
        "authority_commit": AUTHORITY_COMMIT,
        "base_cubic_count": loc2_authority.BASE_CUBIC_COUNT,
        "component_cubic_count": loc2_authority.COMPONENT_CUBIC_COUNT,
        "output_leaves": ["manifest.json", "terminal.json"],
        "runner_id": RUNNER_ID,
        "schema": RAW_SCHEMA,
    }


def _pref1_failure_map(root: Path) -> dict[tuple[int, str], Mapping[str, Any]]:
    raw = _read_repo_leaf(
        root,
        loc1_authority.PREF1_RESULT_PATH,
        maximum=2 * 1024 * 1024,
    )
    if sha256(raw).hexdigest() != loc1_authority.PREF1_RESULT_SHA256:
        _stop("PREF2_PREF1_HASH_DRIFT", loc1_authority.PREF1_RESULT_PATH)
    value = _json(raw, "tracked TDG9 AR1 PREF1 result")
    try:
        retries = value["artifact_payload"]["independent_replay"]["retries"]
    except (KeyError, TypeError) as exc:
        raise TDG9LOC2PREF2Error("PREF2_PREF1_SHAPE_DRIFT", exc) from exc
    answer: dict[tuple[int, str], Mapping[str, Any]] = {}
    for retry in retries:
        for channel in retry["channels"]:
            if channel["classification"] == "sufficient_contraction_failure":
                key = (int(retry["retry"]), str(channel["channel"]))
                if key in answer:
                    _stop("PREF2_PREF1_SHAPE_DRIFT", f"duplicate {key}")
                answer[key] = channel
    if tuple(answer) != loc2_authority.FAILED_OCCURRENCES:
        _stop("PREF2_OCCURRENCE_SELECTION_DRIFT", tuple(answer))
    return answer


def _row_values(segment: Any, block: int, row: int, field: int) -> tuple[float, ...]:
    values = (
        float(segment.left[block, row, field]),
        float(segment.left_rhs[block, row, field]),
        float(segment.right[block, row, field]),
        float(segment.right_rhs[block, row, field]),
        segment.width,
    )
    if any(type(value) is not float or not math.isfinite(value) for value in values):
        _stop("PREF2_BINARY64_ROW_DRIFT", (row, block, field))
    return values


def _rows(surface: Sequence[Sequence[Any]], channel: str) -> Iterator[object]:
    try:
        block_name, field_name = channel.split(":", 1)
        block = _BLOCKS.index(block_name)
        field = _FIELDS.index(field_name)
    except ValueError as exc:
        raise TDG9LOC2PREF2Error("PREF2_CHANNEL_DRIFT", channel) from exc
    if tuple(map(len, surface)) != (1, 2, 4):
        _stop("PREF2_SURFACE_SHAPE_DRIFT", tuple(map(len, surface)))
    outer, medium, fine = surface
    for row in range(loc2_authority.OWNED_ROW_COUNT):
        yield (
            _row_values(outer[0], block, row, field),
            tuple(_row_values(item, block, row, field) for item in medium),
            tuple(_row_values(item, block, row, field) for item in fine),
        )


def _row_hash(rows: Iterable[object]) -> str:
    digest = sha256(b"TDG9-AR1-BINARY64-HERMITE-ROWS-v1\n")
    count = 0
    for ordinal, row in enumerate(rows):
        outer, medium, fine = row  # type: ignore[misc]
        values = [
            value
            for segment in (outer, *medium, *fine)
            for value in segment
        ]
        if len(values) != 35 or any(type(value) is not float for value in values):
            _stop("PREF2_ROW_STREAM_DRIFT", ordinal)
        digest.update(
            (
                f"{ordinal}|"
                + "|".join(value.hex() for value in values)
                + "\n"
            ).encode("ascii")
        )
        count += 1
    if count != loc2_authority.OWNED_ROW_COUNT:
        _stop("PREF2_ROW_COUNT_DRIFT", count)
    return digest.hexdigest()


def _exact_hermite(
    segment: Sequence[float],
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    if len(segment) != 5:
        _stop("PREF2_HERMITE_INPUT_DRIFT", len(segment))
    left, left_rhs, right, right_rhs, width = (
        Fraction(*value.as_integer_ratio()) for value in segment
    )
    slope_left = width * left_rhs
    slope_right = width * right_rhs
    return (
        left,
        slope_left,
        -3 * left - 2 * slope_left + 3 * right - slope_right,
        2 * left + slope_left - 2 * right + slope_right,
    )


def _restrict(
    coefficients: tuple[Fraction, Fraction, Fraction, Fraction], half: int
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    a0, a1, a2, a3 = coefficients
    if half == 0:
        return a0, a1 / 2, a2 / 4, a3 / 8
    if half == 1:
        return (
            a0 + a1 / 2 + a2 / 4 + a3 / 8,
            a1 / 2 + a2 / 2 + 3 * a3 / 8,
            a2 / 4 + 3 * a3 / 8,
            a3 / 8,
        )
    _stop("PREF2_HALF_DRIFT", half)


def _components(
    parent: Sequence[float], child: Sequence[float], half: int
) -> dict[str, tuple[Fraction, Fraction, Fraction, Fraction]]:
    restricted = _restrict(_exact_hermite(parent), half)
    child_coefficients = _exact_hermite(child)
    complete = tuple(
        left - right
        for left, right in zip(restricted, child_coefficients, strict=True)
    )
    parent_left = restricted[0]
    parent_right = sum(restricted, Fraction(0))
    parent_left_slope = restricted[1]
    parent_right_slope = restricted[1] + 2 * restricted[2] + 3 * restricted[3]
    child_left = Fraction(*child[0].as_integer_ratio())
    child_right = Fraction(*child[2].as_integer_ratio())
    child_width = Fraction(*child[4].as_integer_ratio())
    state_left = parent_left - child_left
    state_right = parent_right - child_right
    rhs_left = parent_left_slope - child_width * Fraction(
        *child[1].as_integer_ratio()
    )
    rhs_right = parent_right_slope - child_width * Fraction(
        *child[3].as_integer_ratio()
    )
    value = (
        state_left,
        Fraction(0),
        -3 * state_left + 3 * state_right,
        2 * state_left - 2 * state_right,
    )
    slope = (
        Fraction(0),
        rhs_left,
        -2 * rhs_left - rhs_right,
        rhs_left + rhs_right,
    )
    if tuple(
        left + right for left, right in zip(value, slope, strict=True)
    ) != complete:
        _stop("PREF2_COMPONENT_IDENTITY_DRIFT", half)
    return {"value_V": value, "slope_S": slope, "complete_C": complete}


def _region(row: int) -> str:
    for start, end, name in loc1_authority.ROW_REGIONS:
        if start <= row <= end:
            return name
    _stop("PREF2_ROW_REGION_DRIFT", row)


def _cubics(
    rows: Sequence[object], level: str, component: str
) -> tuple[list[primary.LocalCubic], list[independent_v2.IndependentLocalCubicV2]]:
    if component not in {"value_V", "slope_S", "complete_C"}:
        _stop("PREF2_COMPONENT_DRIFT", component)
    first: list[primary.LocalCubic] = []
    second: list[independent_v2.IndependentLocalCubicV2] = []
    ordinal = 0
    for row, raw in enumerate(rows):
        outer, medium, fine = raw  # type: ignore[misc]
        if level == "D01":
            pairs = (
                (outer, medium[index], index, 0, index) for index in range(2)
            )
        elif level == "D12":
            pairs = (
                (medium[index // 2], fine[index], index, index // 2, index)
                for index in range(4)
            )
        else:
            _stop("PREF2_LEVEL_DRIFT", level)
        for parent, child, subinterval, parent_index, child_index in pairs:
            half = subinterval if level == "D01" else subinterval % 2
            coefficients = _components(parent, child, half)[component]
            metadata = {
                "component": component,
                "level": level,
                "subinterval": subinterval,
                "parent_index": parent_index,
                "child_index": child_index,
                "owned_row": row,
                "global_index": row + 1,
                "radius_hex": float((row + 1) / 16).hex(),
                "row_region": _region(row),
                "ordinal": ordinal,
            }
            first.append(primary.LocalCubic(coefficients, metadata))
            second.append(
                independent_v2.IndependentLocalCubicV2(coefficients, metadata)
            )
            ordinal += 1
    return first, second


def _fraction(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def _exact_json(value: object) -> object:
    if isinstance(value, Fraction):
        return _fraction(value)
    if is_dataclass(value) and not isinstance(value, type):
        return _exact_json(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _exact_json(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_exact_json(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    _stop("PREF2_EXACT_SERIALIZATION_DRIFT", type(value).__name__)


def _candidate_key(item: object) -> tuple[object, ...]:
    return (
        item.polynomial_ordinal,  # type: ignore[attr-defined]
        item.location,  # type: ignore[attr-defined]
        item.location_ordinal,  # type: ignore[attr-defined]
        tuple(item.metadata),  # type: ignore[attr-defined]
    )


def _primary_stationary_count_digest(
    cubics: Sequence[primary.LocalCubic],
) -> str:
    digest = sha256(b"TDG9-LOC2-STATIONARY-COUNT-STREAM-v1\n")
    for ordinal, cubic in enumerate(cubics):
        roots = primary._stationary_intervals(  # type: ignore[attr-defined]
            cubic.coefficients,
            depth=loc2_authority.PRIMARY_REFINEMENT_DEPTH,
        )
        digest.update(f"{ordinal}|{len(roots)}\n".encode("ascii"))
    return digest.hexdigest()


def _route_agreement(
    first: object,
    second: object,
    primary_stationary_digest: str,
    detail: object,
) -> dict[str, Any]:
    first_items = {
        _candidate_key(item): item for item in first.candidates  # type: ignore[attr-defined]
    }
    second_items = {
        _candidate_key(item): item for item in second.candidates  # type: ignore[attr-defined]
    }
    same_keys = (
        len(first_items) == len(first.candidates)  # type: ignore[attr-defined]
        and len(second_items) == len(second.candidates)  # type: ignore[attr-defined]
        and tuple(sorted(first_items)) == tuple(sorted(second_items))
    )
    survivor_overlaps = same_keys and all(
        max(
            first_items[key].parameter_lower,
            second_items[key].parameter_lower,
        )
        <= min(
            first_items[key].parameter_upper,
            second_items[key].parameter_upper,
        )
        and max(
            first_items[key].absolute_lower,
            second_items[key].absolute_lower,
        )
        <= min(
            first_items[key].absolute_upper,
            second_items[key].absolute_upper,
        )
        for key in first_items
    )
    global_overlap = (
        max(first.global_absolute_lower, second.global_absolute_lower)  # type: ignore[attr-defined]
        <= min(first.global_absolute_upper, second.global_absolute_upper)  # type: ignore[attr-defined]
    )
    facts = {
        "polynomial_count_agrees": (
            first.polynomial_count == second.polynomial_count  # type: ignore[attr-defined]
        ),
        "candidate_count_agrees": (
            first.candidate_count == second.candidate_count  # type: ignore[attr-defined]
        ),
        "classification_agrees": (
            first.classification == second.classification  # type: ignore[attr-defined]
        ),
        "stationary_count_digest_agrees": (
            primary_stationary_digest
            == second.stationary_count_stream_sha256  # type: ignore[attr-defined]
        ),
        "survivor_keys_agree": same_keys,
        "survivor_intervals_overlap": survivor_overlaps,
        "global_intervals_overlap": global_overlap,
    }
    if not all(facts.values()):
        _stop("PREF2_ROUTE_DISAGREEMENT", (detail, facts))
    return facts


def _ulp_record(value: float) -> dict[str, object]:
    exact = Fraction(*value.as_integer_ratio())
    lower = Fraction(*math.nextafter(value, -math.inf).as_integer_ratio())
    upper = Fraction(*math.nextafter(value, math.inf).as_integer_ratio())
    return {
        "hex": value.hex(),
        "exact": _fraction(exact),
        "downward": _fraction(exact - lower),
        "upward": _fraction(upper - exact),
    }


def _hermite_inputs(
    rows: Sequence[object], level: str, candidate: object
) -> dict[str, object]:
    metadata = dict(candidate.metadata)  # type: ignore[attr-defined]
    row = int(metadata["owned_row"])
    subinterval = int(metadata["subinterval"])
    outer, medium, fine = rows[row]  # type: ignore[misc]
    if level == "D01":
        parent, child, half = outer, medium[subinterval], subinterval
    elif level == "D12":
        parent = medium[subinterval // 2]
        child = fine[subinterval]
        half = subinterval % 2
    else:
        _stop("PREF2_LEVEL_DRIFT", level)
    restricted = _restrict(_exact_hermite(parent), half)
    components = _components(parent, child, half)
    parent_left = restricted[0]
    parent_right = sum(restricted, Fraction(0))
    parent_left_slope = restricted[1]
    parent_right_slope = restricted[1] + 2 * restricted[2] + 3 * restricted[3]
    child_left = Fraction(*child[0].as_integer_ratio())
    child_right = Fraction(*child[2].as_integer_ratio())
    child_width = Fraction(*child[4].as_integer_ratio())
    child_left_slope = child_width * Fraction(*child[1].as_integer_ratio())
    child_right_slope = child_width * Fraction(*child[3].as_integer_ratio())
    state_left = parent_left - child_left
    state_right = parent_right - child_right
    rhs_left = parent_left_slope - child_left_slope
    rhs_right = parent_right_slope - child_right_slope
    parent_update = parent_right - parent_left
    child_update = Fraction(*child[2].as_integer_ratio()) - Fraction(
        *child[0].as_integer_ratio()
    )
    return {
        "restricted_parent_endpoint_data": {
            "left_value": _fraction(parent_left),
            "right_value": _fraction(parent_right),
            "left_normalized_slope": _fraction(parent_left_slope),
            "right_normalized_slope": _fraction(parent_right_slope),
        },
        "child_endpoint_data": {
            "left_value": _fraction(child_left),
            "right_value": _fraction(child_right),
            "left_normalized_slope": _fraction(child_left_slope),
            "right_normalized_slope": _fraction(child_right_slope),
        },
        "endpoint_state_difference": {
            "left": _fraction(state_left),
            "right": _fraction(state_right),
        },
        "width_scaled_endpoint_RHS_difference": {
            "left": _fraction(rhs_left),
            "right": _fraction(rhs_right),
        },
        "component_polynomials": {
            name: [_fraction(item) for item in coefficients]
            for name, coefficients in components.items()
        },
        "complete_difference_polynomial": [
            _fraction(item) for item in components["complete_C"]
        ],
        "complete_equals_value_plus_slope": True,
        "updates": {
            "restricted_parent": _fraction(parent_update),
            "child": _fraction(child_update),
        },
        "binary64_contributors": {
            "parent_left": _ulp_record(parent[0]),
            "parent_left_RHS": _ulp_record(parent[1]),
            "parent_right": _ulp_record(parent[2]),
            "parent_right_RHS": _ulp_record(parent[3]),
            "parent_width": _ulp_record(parent[4]),
            "child_left": _ulp_record(child[0]),
            "child_left_RHS": _ulp_record(child[1]),
            "child_right": _ulp_record(child[2]),
            "child_right_RHS": _ulp_record(child[3]),
            "child_width": _ulp_record(child[4]),
        },
        "ULP_used_as_tolerance": False,
    }


def _contraction_assessment(outer: object, finest: object) -> dict[str, object]:
    left_lower = outer.global_absolute_lower**2  # type: ignore[attr-defined]
    left_upper = outer.global_absolute_upper**2  # type: ignore[attr-defined]
    right_lower = 8 * finest.global_absolute_lower**2  # type: ignore[attr-defined]
    right_upper = 8 * finest.global_absolute_upper**2  # type: ignore[attr-defined]
    exact_zero = (
        outer.global_absolute_upper == 0  # type: ignore[attr-defined]
        and finest.global_absolute_upper == 0  # type: ignore[attr-defined]
    )
    if exact_zero:
        classification = "exact_zero"
    elif right_upper <= left_lower:
        classification = "sufficient_contraction_pass"
    elif right_lower > left_upper:
        classification = "sufficient_contraction_failure"
    else:
        classification = "threshold_inconclusive"
    return {
        "classification": classification,
        "sufficient_contraction_pass_certified": (
            classification == "sufficient_contraction_pass"
        ),
        "sufficient_contraction_failure_certified": (
            classification == "sufficient_contraction_failure"
        ),
        "threshold_inconclusive": classification == "threshold_inconclusive",
        "left_D01_squared": {
            "lower": _fraction(left_lower),
            "upper": _fraction(left_upper),
        },
        "right_8_D12_squared": {
            "lower": _fraction(right_lower),
            "upper": _fraction(right_upper),
        },
    }


def _component_ownership(
    components: Mapping[str, Mapping[str, object]],
) -> dict[str, object]:
    classes = {
        name: str(value["contraction"]["classification"])  # type: ignore[index]
        for name, value in components.items()
    }
    value_failure = classes["value_V"] == "sufficient_contraction_failure"
    slope_failure = classes["slope_S"] == "sufficient_contraction_failure"
    complete_failure = classes["complete_C"] == "sufficient_contraction_failure"
    value_pass = classes["value_V"] == "sufficient_contraction_pass"
    slope_pass = classes["slope_S"] == "sufficient_contraction_pass"
    value_clear = value_pass or classes["value_V"] == "exact_zero"
    slope_clear = slope_pass or classes["slope_S"] == "exact_zero"
    if complete_failure and value_failure and slope_clear:
        classification = "endpoint_state_owned_failure"
    elif complete_failure and slope_failure and value_clear:
        classification = "width_scaled_RHS_owned_failure"
    elif complete_failure and value_failure and slope_failure:
        classification = "both_components_independently_fail"
    elif complete_failure and value_pass and slope_pass:
        classification = "mixed_only_failure"
    elif complete_failure:
        classification = "component_ownership_inconclusive"
    else:
        classification = "complete_not_failure"
    cancellation: dict[str, bool] = {}
    reinforcement: dict[str, bool] = {}
    for level in ("D01", "D12"):
        value = components["value_V"][level]
        slope = components["slope_S"][level]
        complete = components["complete_C"][level]
        cancellation[level] = complete.global_absolute_upper < max(  # type: ignore[attr-defined]
            value.global_absolute_lower, slope.global_absolute_lower  # type: ignore[attr-defined]
        )
        reinforcement[level] = complete.global_absolute_lower > max(  # type: ignore[attr-defined]
            value.global_absolute_upper, slope.global_absolute_upper  # type: ignore[attr-defined]
        )
    return {
        "classification": classification,
        "value_only_failure": complete_failure and value_failure and slope_clear,
        "slope_only_failure": complete_failure and slope_failure and value_clear,
        "both_components_independently_fail": value_failure and slope_failure,
        "mixed_only_failure": complete_failure and value_pass and slope_pass,
        "cancellation_certified_by_level": cancellation,
        "reinforcement_certified_by_level": reinforcement,
    }


def _interval_record(item: object) -> dict[str, object]:
    return {
        "polynomial_ordinal": item.polynomial_ordinal,  # type: ignore[attr-defined]
        "location": item.location,  # type: ignore[attr-defined]
        "location_ordinal": item.location_ordinal,  # type: ignore[attr-defined]
        "metadata": _exact_json(item.metadata),  # type: ignore[attr-defined]
        "parameter_lower": _fraction(item.parameter_lower),  # type: ignore[attr-defined]
        "parameter_upper": _fraction(item.parameter_upper),  # type: ignore[attr-defined]
        "absolute_lower": _fraction(item.absolute_lower),  # type: ignore[attr-defined]
        "absolute_upper": _fraction(item.absolute_upper),  # type: ignore[attr-defined]
    }


def _component_summary(
    *,
    component: str,
    level: str,
    rows: Sequence[object],
    raw_component: Mapping[str, Any],
    component_evidence: dict[str, dict[str, object]],
) -> tuple[dict[str, Any], dict[str, int]]:
    expected = (
        2 * loc2_authority.OWNED_ROW_COUNT
        if level == "D01"
        else 4 * loc2_authority.OWNED_ROW_COUNT
    )
    first_cubics, second_cubics = _cubics(rows, level, component)
    if len(first_cubics) != expected or len(second_cubics) != expected:
        _stop("PREF2_POLYNOMIAL_COUNT_DRIFT", (component, level))
    first = primary.localize_absolute_maximum(
        first_cubics,
        maximum_candidates=4 * expected,
        refinement_depth=loc2_authority.PRIMARY_REFINEMENT_DEPTH,
    )
    try:
        second = independent_v2.localize_absolute_maximum_independently_v2(
            second_cubics,
            maximum_candidates=4 * expected,
        )
    except independent_v2.RootIsolationInconclusive as exc:
        raise TDG9LOC2PREF2Error(
            "PREF2_V2_LOCALIZATION_INCONCLUSIVE",
            (component, level, exc.reason, exc.detail),
        ) from exc
    stationary_digest = _primary_stationary_count_digest(first_cubics)
    agreement = _route_agreement(
        first,
        second,
        stationary_digest,
        (component, level),
    )
    primary_json = _exact_json(first)
    second_json = _exact_json(second)
    decompositions = [_hermite_inputs(rows, level, item) for item in first.candidates]
    component_evidence[component][level] = first
    expected_raw = {
        "polynomial_count": expected,
        "primary_stationary_count_stream_sha256": stationary_digest,
        "primary": primary_json,
        "independent_v2": second_json,
        "co_maximizer_count": len(first.candidates),
        "maximizer_input_decomposition": decompositions,
    }
    # The contraction is appended after both levels have been evaluated.
    if any(raw_component.get(key) != value for key, value in expected_raw.items()):
        _stop("PREF2_RAW_COMPONENT_DRIFT", (component, level))
    first_items = {_candidate_key(item): item for item in first.candidates}
    second_items = {_candidate_key(item): item for item in second.candidates}
    survivor_keys = [
        {
            "polynomial_ordinal": key[0],
            "location": key[1],
            "location_ordinal": key[2],
            "metadata": _exact_json(key[3]),
        }
        for key in sorted(first_items)
    ]
    interval_pairs = [
        {
            "key": survivor_keys[index],
            "primary": _interval_record(first_items[key]),
            "independent_v2": _interval_record(second_items[key]),
        }
        for index, key in enumerate(sorted(first_items))
    ]
    root_methods: dict[str, int] = {}
    for item in second.candidates:
        root_methods[item.root_method] = root_methods.get(item.root_method, 0) + 1
    summary = {
        "component": component,
        "level": level,
        "polynomial_count": expected,
        "candidate_count": first.candidate_count,
        "co_maximizer_count": len(first.candidates),
        "classification": first.classification,
        "primary_stationary_count_stream_sha256": stationary_digest,
        "independent_v2_stationary_count_stream_sha256": (
            second.stationary_count_stream_sha256
        ),
        "route_agreement": agreement,
        "primary_evidence_sha256": sha256(_canonical(primary_json)).hexdigest(),
        "independent_v2_evidence_sha256": sha256(
            _canonical(second_json)
        ).hexdigest(),
        "survivor_key_stream_sha256": sha256(
            _canonical(survivor_keys)
        ).hexdigest(),
        "survivor_interval_pair_stream_sha256": sha256(
            _canonical(interval_pairs)
        ).hexdigest(),
        "maximizer_input_decomposition_sha256": sha256(
            _canonical(decompositions)
        ).hexdigest(),
        "v2_root_method_counts": root_methods,
    }
    return summary, root_methods


def _occurrence_summary(
    *,
    rows: Sequence[object],
    retry: int,
    channel: str,
    pref1: Mapping[str, Any],
    raw_occurrence: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, int]]:
    row_hash = _row_hash(rows)
    if row_hash != pref1.get("row_stream_sha256"):
        _stop("PREF2_PREF1_ROW_DRIFT", (retry, channel))
    for raw_key, pref1_key in (
        ("PREF1_row_stream_sha256", "row_stream_sha256"),
        ("PREF1_primary_evidence_sha256", "primary_evidence_sha256"),
        ("PREF1_independent_evidence_sha256", "independent_evidence_sha256"),
    ):
        if raw_occurrence.get(raw_key) != pref1.get(pref1_key):
            _stop("PREF2_PREF1_BINDING_DRIFT", (retry, channel, raw_key))
    if (
        raw_occurrence.get("retry") != retry
        or raw_occurrence.get("channel") != channel
    ):
        _stop("PREF2_OCCURRENCE_ORDER_DRIFT", (retry, channel))
    raw_levels = raw_occurrence.get("levels")
    if not isinstance(raw_levels, list) or len(raw_levels) != 2:
        _stop("PREF2_RAW_LEVEL_SHAPE_DRIFT", (retry, channel))
    component_evidence: dict[str, dict[str, object]] = {
        "value_V": {},
        "slope_S": {},
        "complete_C": {},
    }
    component_summaries: list[dict[str, Any]] = []
    root_methods: dict[str, int] = {}
    recomputed_levels: list[dict[str, Any]] = []
    for level, raw_level in zip(("D01", "D12"), raw_levels, strict=True):
        base_count = (
            2 * loc2_authority.OWNED_ROW_COUNT
            if level == "D01"
            else 4 * loc2_authority.OWNED_ROW_COUNT
        )
        if (
            not isinstance(raw_level, Mapping)
            or raw_level.get("level") != level
            or raw_level.get("base_cubic_count") != base_count
            or not isinstance(raw_level.get("components"), Mapping)
            or set(raw_level["components"])  # type: ignore[index]
            != {"value_V", "slope_S", "complete_C"}
        ):
            _stop("PREF2_RAW_LEVEL_DRIFT", (retry, channel, level))
        components: dict[str, Any] = {}
        for component in ("value_V", "slope_S", "complete_C"):
            raw_component = raw_level["components"][component]  # type: ignore[index]
            if not isinstance(raw_component, Mapping):
                _stop("PREF2_RAW_COMPONENT_DRIFT", (retry, channel, component))
            summary, methods = _component_summary(
                component=component,
                level=level,
                rows=rows,
                raw_component=raw_component,
                component_evidence=component_evidence,
            )
            component_summaries.append(summary)
            for name, count in methods.items():
                root_methods[name] = root_methods.get(name, 0) + count
            components[component] = {
                key: raw_component[key]
                for key in (
                    "polynomial_count",
                    "primary_stationary_count_stream_sha256",
                    "primary",
                    "independent_v2",
                    "co_maximizer_count",
                    "maximizer_input_decomposition",
                )
            }
        recomputed_levels.append(
            {
                "level": level,
                "base_cubic_count": base_count,
                "components": components,
            }
        )
    for component in ("value_V", "slope_S", "complete_C"):
        contraction = _contraction_assessment(
            component_evidence[component]["D01"],
            component_evidence[component]["D12"],
        )
        component_evidence[component]["contraction"] = contraction
        for level_index in (0, 1):
            recomputed_levels[level_index]["components"][component][
                "contraction"
            ] = contraction
            raw_contraction = raw_levels[level_index]["components"][component][  # type: ignore[index]
                "contraction"
            ]
            if raw_contraction != contraction:
                _stop(
                    "PREF2_CONTRACTION_DRIFT",
                    (retry, channel, component, level_index),
                )
        for summary in component_summaries:
            if summary["component"] == component:
                summary["contraction"] = contraction["classification"]
                summary["contraction_evidence_sha256"] = sha256(
                    _canonical(contraction)
                ).hexdigest()
    ownership = _component_ownership(component_evidence)
    expected_occurrence = {
        "retry": retry,
        "channel": channel,
        "PREF1_row_stream_sha256": pref1["row_stream_sha256"],
        "PREF1_primary_evidence_sha256": pref1["primary_evidence_sha256"],
        "PREF1_independent_evidence_sha256": pref1[
            "independent_evidence_sha256"
        ],
        "levels": recomputed_levels,
        "component_ownership": ownership,
    }
    if expected_occurrence != raw_occurrence:
        _stop("PREF2_RAW_OCCURRENCE_DRIFT", (retry, channel))
    return (
        {
            "retry": retry,
            "channel": channel,
            "PREF1_row_stream_sha256": pref1["row_stream_sha256"],
            "PREF1_primary_evidence_sha256": pref1["primary_evidence_sha256"],
            "PREF1_independent_evidence_sha256": pref1[
                "independent_evidence_sha256"
            ],
            "raw_occurrence_sha256": sha256(
                _canonical(raw_occurrence)
            ).hexdigest(),
            "component_ownership": ownership,
            "component_ownership_sha256": sha256(
                _canonical(ownership)
            ).hexdigest(),
            "components": component_summaries,
        },
        root_methods,
    )


def expected_evidence(config_raw: bytes) -> dict[str, Any]:
    config = _config(config_raw)
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "target_protocol": TARGET_PROTOCOL,
        "raw": dict(config["raw"]),
        "authority": dict(config["authority"]),
        "store": dict(config["store"]),
        "selection": dict(config["selection"]),
        "work_budget": dict(config["work_budget"]),
        "scope": dict(config["scope"]),
        "claims": dict(config["claims"]),
    }


def bind_raw_result(config_raw: bytes, repository: Path) -> dict[str, Any]:
    expected = expected_evidence(config_raw)
    root = Path(os.path.abspath(os.fspath(repository)))
    try:
        root_metadata = root.lstat()
    except OSError as exc:
        raise TDG9LOC2PREF2Error("PREF2_REPOSITORY_ABSENT", root) from exc
    if stat.S_ISLNK(root_metadata.st_mode) or not stat.S_ISDIR(root_metadata.st_mode):
        _stop("PREF2_REPOSITORY_UNSAFE", root)
    manifest_raw, terminal_raw, manifest, terminal = _raw_snapshot(root)
    if manifest != _expected_manifest():
        _stop("PREF2_RAW_MANIFEST_DRIFT", "manifest schema differs")
    _authenticate_authority(root, manifest)
    before_store = _snapshot_store(root)
    sealed_store = (SEALED_STORE_LEAF_COUNT, SEALED_STORE_SHA256)
    if before_store != sealed_store:
        _stop("PREF2_STORE_IDENTITY_DRIFT", before_store)
    failures = _pref1_failure_map(root)
    raw_occurrences = terminal.get("occurrences")
    if not isinstance(raw_occurrences, list) or len(raw_occurrences) != 10:
        _stop("PREF2_RAW_OCCURRENCE_SHAPE_DRIFT", type(raw_occurrences).__name__)
    store = HLT16CampaignStore(root / ar1.PREF2_STORE_PATH)
    shells = build_static_gr0_shells(root)
    coordinates = shells[loc2_authority.MEMBER_KEY].initial.grid.coordinates
    expected_coordinates = (
        np.arange(loc2_authority.POINT_COUNT, dtype=np.float64)
        * float.fromhex(loc2_authority.GRID_SPACING_HEX)
    )
    if (
        type(coordinates) is not np.ndarray
        or coordinates.dtype.str != "<f8"
        or coordinates.shape != (loc2_authority.POINT_COUNT,)
        or not np.array_equal(coordinates, expected_coordinates)
        or float(coordinates[-1]).hex() != loc2_authority.OUTER_RADIUS_HEX
    ):
        _stop("PREF2_GRID_IDENTITY_DRIFT", getattr(coordinates, "shape", None))
    surfaces: dict[int, tuple[tuple[Any, ...], ...]] = {}
    historical_hashes: dict[str, str] = {}
    for replay in ar1.REPLAYS:
        try:
            prepared, historical = ar1_pref1._prepare(root, store, shells, replay)
            surface = ar1_pref1._surface(prepared)
        except Exception as exc:
            raise TDG9LOC2PREF2Error(
                "PREF2_REPLAY_RECONSTRUCTION_DRIFT", replay["retry"]
            ) from exc
        retry = int(replay["retry"])
        surfaces[retry] = surface
        historical_hashes[str(retry)] = sha256(
            _canonical(historical)
        ).hexdigest()
    occurrence_summaries: list[dict[str, Any]] = []
    component_candidate_counts = {
        "complete_C": 0,
        "slope_S": 0,
        "value_V": 0,
    }
    route_classifications: dict[str, int] = {}
    root_methods: dict[str, int] = {}
    base_cubic_count = 0
    route_agreement_count = 0
    survivor_count = 0
    component_contraction_failure_count = 0
    ownership_both_count = 0
    for (retry, channel), raw_occurrence in zip(
        loc2_authority.FAILED_OCCURRENCES,
        raw_occurrences,
        strict=True,
    ):
        if not isinstance(raw_occurrence, Mapping):
            _stop("PREF2_RAW_OCCURRENCE_SHAPE_DRIFT", (retry, channel))
        rows = tuple(_rows(surfaces[retry], channel))
        summary, methods = _occurrence_summary(
            rows=rows,
            retry=retry,
            channel=channel,
            pref1=failures[(retry, channel)],
            raw_occurrence=raw_occurrence,
        )
        occurrence_summaries.append(summary)
        if summary["component_ownership"]["classification"] == (
            "both_components_independently_fail"
        ):
            ownership_both_count += 1
        for item in summary["components"]:
            base_cubic_count += (
                item["polynomial_count"] if item["component"] == "complete_C" else 0
            )
            component_candidate_counts[item["component"]] += item["candidate_count"]
            route_classifications[item["classification"]] = (
                route_classifications.get(item["classification"], 0) + 1
            )
            route_agreement_count += int(all(item["route_agreement"].values()))
            survivor_count += item["co_maximizer_count"]
            component_contraction_failure_count += int(
                item["level"] == "D01"
                and item["contraction"] == "sufficient_contraction_failure"
            )
        for name, count in methods.items():
            root_methods[name] = root_methods.get(name, 0) + count
    total_candidates = sum(component_candidate_counts.values())
    if (
        base_cubic_count != loc2_authority.BASE_CUBIC_COUNT
        or component_candidate_counts != EXPECTED_COMPONENT_CANDIDATES
        or total_candidates != EXPECTED_TOTAL_CANDIDATES
        or any(
            count > loc2_authority.FAMILY_CANDIDATE_CEILING
            for count in component_candidate_counts.values()
        )
        or total_candidates > loc2_authority.AGGREGATE_CANDIDATE_CEILING
        or route_classifications != EXPECTED_ROUTE_CLASSIFICATIONS
        or root_methods != EXPECTED_ROOT_METHODS
        or route_agreement_count != 60
        or survivor_count != EXPECTED_SURVIVOR_COUNT
        or component_contraction_failure_count != 30
        or ownership_both_count != 10
    ):
        _stop(
            "PREF2_AGGREGATE_REDUCTION_DRIFT",
            {
                "base": base_cubic_count,
                "candidates": component_candidate_counts,
                "routes": route_classifications,
                "roots": root_methods,
                "survivors": survivor_count,
                "contractions": component_contraction_failure_count,
                "ownership": ownership_both_count,
            },
        )
    expected_terminal = {
        "PDE_state_committed": False,
        "SSPRK3_comparator_executed": False,
        "ULP_used_as_tolerance": False,
        "artifact_id": loc2_authority.ARTIFACT_ID,
        "authority_commit": AUTHORITY_COMMIT,
        "base_cubic_count": base_cubic_count,
        "candidate_branch_opened": False,
        "classification": RAW_CLASSIFICATION,
        "component_candidate_counts": component_candidate_counts,
        "evaluator_inconclusive": False,
        "fourth_width_executed": False,
        "historical_TDG6_evidence_sha256": historical_hashes,
        "occurrences": raw_occurrences,
        "physical_result_earned": False,
        "runner_id": RUNNER_ID,
        "schema": RAW_SCHEMA,
        "stage2_source_decomposition_executed": False,
        "store_snapshot_after": {
            "leaf_count": SEALED_STORE_LEAF_COUNT,
            "sha256": SEALED_STORE_SHA256,
        },
        "store_snapshot_before": {
            "leaf_count": SEALED_STORE_LEAF_COUNT,
            "sha256": SEALED_STORE_SHA256,
        },
        "store_unchanged": True,
        "total_VSC_candidate_count": total_candidates,
    }
    if terminal != expected_terminal:
        _stop("PREF2_RAW_TERMINAL_DRIFT", "terminal schema or nonclaim differs")
    after_store = _snapshot_store(root)
    manifest_after, terminal_after = _raw_pair(root)
    if (
        after_store != before_store
        or manifest_after != manifest_raw
        or terminal_after != terminal_raw
    ):
        _stop("PREF2_INPUT_MUTATED_DURING_BIND", "raw/store input changed")
    return {
        **expected,
        "raw_binding": {
            "manifest_sha256": sha256(manifest_raw).hexdigest(),
            "terminal_sha256": sha256(terminal_raw).hexdigest(),
            "leaf_count": 2,
            "raw_namespace_unchanged": True,
            "canonical_duplicate_free_schema_verified": True,
        },
        "authority_binding": {
            "commit": AUTHORITY_COMMIT,
            "parent_commit": AUTHORITY_PARENT,
            "blob_count": len(_AUTHORITY_BLOBS),
            "all_committed_blob_hashes_verified": True,
        },
        "store_binding": {
            "binder_snapshot_before": {
                "leaf_count": before_store[0],
                "sha256": before_store[1],
            },
            "binder_snapshot_after": {
                "leaf_count": after_store[0],
                "sha256": after_store[1],
            },
            "raw_snapshot_before": terminal["store_snapshot_before"],
            "raw_snapshot_after": terminal["store_snapshot_after"],
            "store_unchanged": True,
        },
        "independent_replay": {
            "occurrence_count": 10,
            "level_count": 20,
            "component_level_count": 60,
            "base_cubic_count": base_cubic_count,
            "component_cubic_count": base_cubic_count * 3,
            "component_candidate_counts": component_candidate_counts,
            "total_VSC_candidate_count": total_candidates,
            "route_agreement_count": route_agreement_count,
            "route_classification_counts": route_classifications,
            "survivor_count": survivor_count,
            "v2_root_method_counts": root_methods,
            "component_contraction_failure_count": (
                component_contraction_failure_count
            ),
            "both_components_independently_fail_count": ownership_both_count,
            "historical_TDG6_evidence_sha256": historical_hashes,
            "occurrences": occurrence_summaries,
        },
        "conclusion": {
            "LOC2_bounded_exact_localization_completed": True,
            "primary_v2_agree_on_all_60_component_levels": True,
            "stationary_count_digests_agree_on_all_60_component_levels": True,
            "survivor_keys_and_intervals_agree_on_all_60_component_levels": True,
            "all_10_failures_have_endpoint_state_and_width_scaled_RHS_failures": True,
            "mixed_only_failure_count": 0,
            "successor_remedy_selected": False,
            "PDE_or_candidate_state_opened": False,
            "physics_inference_permitted": False,
        },
    }


def build_pref2_result(
    config_raw: bytes,
    repository: Path,
    *,
    live: bool = True,
) -> dict[str, Any]:
    payload = (
        bind_raw_result(config_raw, repository)
        if live
        else expected_evidence(config_raw)
    )
    return {"artifact_id": ARTIFACT_ID, "artifact_payload": payload}


def validate_compact_result(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    expected = expected_evidence(config_raw)
    if sha256(result_raw).hexdigest() != COMPACT_RESULT_SHA256:
        _stop("PREF2_COMPACT_DRIFT", "compact result hash differs")
    result = _json(result_raw, "compact result")
    if (
        result.get("artifact_id") != ARTIFACT_ID
        or not isinstance(result.get("artifact_payload"), Mapping)
    ):
        _stop("PREF2_COMPACT_DRIFT", "compact identity differs")
    payload = result["artifact_payload"]
    for key, value in expected.items():
        if payload.get(key) != value:
            _stop("PREF2_COMPACT_DRIFT", key)
    raw_binding = payload.get("raw_binding")
    if raw_binding != {
        "canonical_duplicate_free_schema_verified": True,
        "leaf_count": 2,
        "manifest_sha256": RAW_MANIFEST_SHA256,
        "raw_namespace_unchanged": True,
        "terminal_sha256": RAW_TERMINAL_SHA256,
    }:
        _stop("PREF2_COMPACT_DRIFT", "raw binding")
    replay = payload.get("independent_replay")
    if not isinstance(replay, Mapping):
        _stop("PREF2_COMPACT_DRIFT", "independent replay absent")
    expected_aggregate = {
        "occurrence_count": 10,
        "level_count": 20,
        "component_level_count": 60,
        "base_cubic_count": loc2_authority.BASE_CUBIC_COUNT,
        "component_cubic_count": loc2_authority.COMPONENT_CUBIC_COUNT,
        "component_candidate_counts": EXPECTED_COMPONENT_CANDIDATES,
        "total_VSC_candidate_count": EXPECTED_TOTAL_CANDIDATES,
        "route_agreement_count": 60,
        "route_classification_counts": EXPECTED_ROUTE_CLASSIFICATIONS,
        "survivor_count": EXPECTED_SURVIVOR_COUNT,
        "v2_root_method_counts": EXPECTED_ROOT_METHODS,
        "component_contraction_failure_count": 30,
        "both_components_independently_fail_count": 10,
        "historical_TDG6_evidence_sha256": {
            "3": "0d21f650ccc46db778fd81fbebf5acab6394e1c2f7972b7940e24c76af9e1f59",
            "4": "fc68ce01980cb7a66f82de59bbd59d5e53ace08537b6b0dbfe0be1e070e6912f",
            "5": "b602d12ab00757090734ac63042f172a372154d51b41c7d286c7c0a758b62221",
        },
    }
    for key, value in expected_aggregate.items():
        if replay.get(key) != value:
            _stop("PREF2_COMPACT_DRIFT", f"aggregate {key}")
    occurrences = replay.get("occurrences")
    if (
        not isinstance(occurrences, list)
        or len(occurrences) != 10
        or [
            (item.get("retry"), item.get("channel")) for item in occurrences
        ]
        != list(loc2_authority.FAILED_OCCURRENCES)
        or any(len(item.get("components", ())) != 6 for item in occurrences)
    ):
        _stop("PREF2_COMPACT_DRIFT", "occurrence summaries")
    conclusion = payload.get("conclusion")
    if conclusion != {
        "LOC2_bounded_exact_localization_completed": True,
        "PDE_or_candidate_state_opened": False,
        "all_10_failures_have_endpoint_state_and_width_scaled_RHS_failures": True,
        "mixed_only_failure_count": 0,
        "physics_inference_permitted": False,
        "primary_v2_agree_on_all_60_component_levels": True,
        "stationary_count_digests_agree_on_all_60_component_levels": True,
        "successor_remedy_selected": False,
        "survivor_keys_and_intervals_agree_on_all_60_component_levels": True,
    }:
        _stop("PREF2_COMPACT_DRIFT", "conclusion")
    return result


__all__ = (
    "ARTIFACT_ID",
    "CONFIG_PATH",
    "OWNER_DOCUMENT",
    "RESULT_PATH",
    "TDG9LOC2PREF2Error",
    "bind_raw_result",
    "build_pref2_result",
    "canonical_result",
    "expected_evidence",
    "validate_compact_result",
)
