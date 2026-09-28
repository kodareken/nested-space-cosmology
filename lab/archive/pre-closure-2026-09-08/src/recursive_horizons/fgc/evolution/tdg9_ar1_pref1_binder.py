"""Independent, outcome-neutral binder for the completed TDG9 AR1 diagnostic.

The live path is intentionally expensive and store-reading.  It authenticates
the two raw leaves, reconstructs the three selected GR-0 TDG6 compositors from
their persisted predecessors, and invokes both exact arithmetic evaluators.
The compact path validates only tracked evidence and never opens either raw
namespace.  This module imports neither the AR1 runner nor TDG8's persisted
retry replay wrapper.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass, is_dataclass
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

from . import proto15_runtime as p15
from . import tdg5_stage_complete_refinement_runtime as tdg5
from . import tdg6_temporal_admission_runtime as tdg6
from . import tdg7_binary64_subdivision_lattice as lattice
from . import tdg8_rcv3_pref2_binder as pref2
from . import tdg9_ar1_authority as authority
from . import tdg9_exact_temporal_arithmetic as exact_primary
from . import tdg9_exact_temporal_arithmetic_independent as exact_independent
from .hlt16_campaign_store import HLT16CampaignStore
from .hlt16_member_codec import HLT16MemberSnapshot, restore_member as restore_snapshot
from .proto19_gr0_static_factory import build_static_gr0_shells


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG9-AR1-PREF1"
CLASSIFICATION = "independently_bound_exact_arithmetic_discriminator_terminal"
CONFIG_PATH = "configs/fgc/fgc-1-tdg9-ar1-pref1.toml"
RESULT_PATH = "results/fgc-1-tdg9-ar1-pref1.json"
OWNER_DOCUMENT = "docs/fgc-tdg9-ar1-pref1.md"
RAW_NAMESPACE = authority.OUTPUT_NAMESPACE
RAW_MANIFEST_SHA256 = "13f2bf0c89a2303da4c2b5b6a376688676b3b956abcfa52533c0c19fb6a0a77d"
RAW_TERMINAL_SHA256 = "c864b2bffc61451b2f3e183e481425d07e19a1409ee81d327da958ab14994432"
COMPACT_RESULT_SHA256 = "e061bcbcebf0919ef237cb1e757b0fd9c40c0b37348b3bfb3409354c97b84bdf"
AUTHORITY_COMMIT = "fd475aee1870ad60fc8cbaeacb3e33d7b97f831e"
RAW_CLASSIFICATION = "legacy_enclosure_not_sole_owner_on_frozen_samples"
RAW_SCHEMA = "FGC-1-TDG9-AR1-raw-result-v1"
RUNNER_ID = "FGC-1-TDG9-AR1-RUN1"
_MAX_LEAF_BYTES = 16 * 1024 * 1024
_FIELD_NAMES = ("alpha", "v", "lambda", "R", "phi", "chi")
_BLOCK_NAMES = ("u", "p", "q")
_EXPECTED_FAILURES = {
    3: ("u:alpha", "u:R"),
    4: ("u:alpha", "u:lambda", "u:R"),
    5: ("u:alpha", "u:v", "u:lambda", "u:R", "q:R"),
}


class TDG9AR1PREF1Error(ValueError):
    """The compact contract, raw result, store, or independent replay differs."""

    def __init__(self, stop_id: str, detail: object) -> None:
        self.stop_id = str(stop_id)
        self.detail = " ".join(str(detail).split())[:640]
        super().__init__(f"{self.stop_id}: {self.detail}")


def _stop(stop_id: str, detail: object) -> NoReturn:
    raise TDG9AR1PREF1Error(stop_id, detail)


def _canonical(value: object) -> bytes:
    try:
        return json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
            allow_nan=False,
        ).encode("ascii")
    except (TypeError, ValueError) as exc:
        raise TDG9AR1PREF1Error("PREF1_CANONICAL_DRIFT", exc) from exc


def canonical_result(value: object) -> bytes:
    try:
        return (
            json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False)
            + "\n"
        ).encode("ascii")
    except (TypeError, ValueError) as exc:
        raise TDG9AR1PREF1Error("PREF1_COMPACT_DRIFT", exc) from exc


def _pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in items:
        if key in answer:
            raise ValueError(f"duplicate key {key}")
        answer[key] = value
    return answer


def _json(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(raw.decode("ascii"), object_pairs_hook=_pairs,
                           parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TDG9AR1PREF1Error("PREF1_RAW_JSON_DRIFT", label) from exc
    if not isinstance(value, dict) or canonical_result(value) != raw:
        _stop("PREF1_RAW_JSON_DRIFT", f"{label} is not canonical pretty JSON")
    return value


def _config(raw: bytes) -> dict[str, Any]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise TDG9AR1PREF1Error("PREF1_CONFIG_DRIFT", exc) from exc
    expected = {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "target_protocol": authority.TARGET_PROTOCOL,
        "owner_document": OWNER_DOCUMENT,
        "raw": {
            "namespace": RAW_NAMESPACE,
            "manifest_sha256": RAW_MANIFEST_SHA256,
            "terminal_sha256": RAW_TERMINAL_SHA256,
            "classification": RAW_CLASSIFICATION,
            "authority_commit": AUTHORITY_COMMIT,
            "leaf_count": 2,
        },
        "replay": {
            "retries": [3, 4, 5],
            "channel_count_per_retry": 18,
            "owned_row_count": authority.OWNED_ROW_COUNT,
            "evaluator_calls": authority.MAX_EVALUATOR_CALLS,
            "max_depth": authority.MAX_DEPTH,
            "max_nodes_per_channel_per_evaluator": authority.MAX_NODES_PER_CHANNEL_PER_EVALUATOR,
        },
        "scope": {
            "live_store_read_only": True,
            "raw_namespace_read_only": True,
            "compact_verifier_store_blind": True,
            "runner_decision_code_imported": False,
            "persisted_retry_wrapper_imported": False,
            "campaign_store_mutation_authorized": False,
            "PDE_state_commit_authorized": False,
            "continuation_authorized": False,
            "candidate_branches_authorized": False,
        },
        "claims": {
            "raw_AR1_result_independently_bound": True,
            "AR1_binder_completed": True,
            "legacy_enclosure_sole_owner_rejected_on_frozen_samples": True,
            "common_event_completed": False,
            "GR0_calibration_completed": False,
            "candidate_execution_authorized": False,
            "mechanism_result_earned": False,
            "physical_result_earned": False,
            "retained_EFT_evolution_authorized": False,
            "physical_transition_claim_authorized": False,
        },
    }
    if value != expected:
        _stop("PREF1_CONFIG_DRIFT", "typed PREF1 config differs")
    return value


def _read_leaf(root: Path, relative: str, maximum: int = _MAX_LEAF_BYTES) -> bytes:
    candidate = Path(relative)
    if candidate.is_absolute() or any(part in {"", ".", ".."} for part in candidate.parts):
        _stop("PREF1_PATH_UNSAFE", relative)
    current = root
    for part in candidate.parts[:-1]:
        current = current / part
        try:
            metadata = current.lstat()
        except OSError as exc:
            raise TDG9AR1PREF1Error("PREF1_PATH_ABSENT", relative) from exc
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
            _stop("PREF1_PATH_UNSAFE", relative)
    path = root / candidate
    try:
        before = path.lstat()
    except OSError as exc:
        raise TDG9AR1PREF1Error("PREF1_LEAF_ABSENT", relative) from exc
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode) or before.st_size > maximum:
        _stop("PREF1_LEAF_UNSAFE", relative)
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        active = os.fstat(descriptor)
        identity = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
        if identity != (active.st_dev, active.st_ino, active.st_size, active.st_mtime_ns):
            _stop("PREF1_LEAF_RACED", relative)
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
    if identity != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
        _stop("PREF1_LEAF_CHANGED", relative)
    raw = b"".join(chunks)
    if len(raw) != before.st_size or len(raw) > maximum:
        _stop("PREF1_LEAF_SIZE", relative)
    return raw


def _raw_snapshot(root: Path) -> tuple[bytes, bytes, dict[str, Any], dict[str, Any]]:
    target = root / RAW_NAMESPACE
    try:
        metadata = target.lstat()
        names = tuple(sorted(item.name for item in os.scandir(target)))
    except OSError as exc:
        raise TDG9AR1PREF1Error("PREF1_RAW_NAMESPACE_ABSENT", RAW_NAMESPACE) from exc
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode) or names != ("manifest.json", "terminal.json"):
        _stop("PREF1_RAW_NAMESPACE_PARTIAL_OR_UNSAFE", names)
    manifest = _read_leaf(root, f"{RAW_NAMESPACE}/manifest.json")
    terminal = _read_leaf(root, f"{RAW_NAMESPACE}/terminal.json")
    if sha256(manifest).hexdigest() != RAW_MANIFEST_SHA256 or sha256(terminal).hexdigest() != RAW_TERMINAL_SHA256:
        _stop("PREF1_RAW_HASH_DRIFT", "raw two-leaf hash differs")
    return manifest, terminal, _json(manifest, "manifest"), _json(terminal, "terminal")


def _git(root: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ["git", "--no-replace-objects", "--no-optional-locks", "-C", str(root), *arguments],
        check=True, capture_output=True, text=True,
        env={"PATH": os.environ.get("PATH", ""), "LC_ALL": "C"},
    )
    return completed.stdout.strip()


def _git_bytes(root: Path, *arguments: str) -> bytes:
    completed = subprocess.run(
        ["git", "--no-replace-objects", "--no-optional-locks", "-C", str(root), *arguments],
        check=True, capture_output=True,
        env={"PATH": os.environ.get("PATH", ""), "LC_ALL": "C"},
    )
    return completed.stdout


def _authenticate_authority(root: Path, manifest: Mapping[str, Any]) -> None:
    if manifest.get("authority_commit") != AUTHORITY_COMMIT or _git(root, "cat-file", "-t", AUTHORITY_COMMIT) != "commit":
        _stop("PREF1_AUTHORITY_COMMIT_DRIFT", manifest.get("authority_commit"))
    parents = _git(root, "show", "-s", "--format=%P", AUTHORITY_COMMIT).split()
    if parents != [authority.EVALUATOR_COMMIT]:
        _stop("PREF1_AUTHORITY_LINEAGE_DRIFT", parents)
    for relative, digest in (
        (authority.CONFIG_PATH, "db44602e79dfd79458292fd6482b4520b9369f751f9f73540e20dd636fea13f5"),
        (authority.RESULT_PATH, "2ce42ad7f571ade72cb705ed4577affa9b6b86eebc3b730959c765c77e25cfaa"),
    ):
        raw = _git_bytes(root, "show", f"{AUTHORITY_COMMIT}:{relative}")
        if sha256(raw).hexdigest() != digest:
            _stop("PREF1_AUTHORITY_BLOB_DRIFT", relative)


def _authenticate_pref2(root: Path) -> str:
    config = _read_leaf(root, authority.PREF2_CONFIG_PATH, 2 * 1024 * 1024)
    result = _read_leaf(root, authority.PREF2_RESULT_PATH, 8 * 1024 * 1024)
    if sha256(config).hexdigest() != authority.PREF2_CONFIG_SHA256 or sha256(result).hexdigest() != authority.PREF2_RESULT_SHA256:
        _stop("PREF1_PREF2_BLOB_DRIFT", "tracked PREF2 differs")
    regenerated = pref2.build_pref2_result(config, root)
    if pref2.canonical_result(regenerated) != result:
        _stop("PREF1_PREF2_LIVE_DRIFT", "live store differs from PREF2")
    store = regenerated["artifact_payload"]["terminal_evidence"]["store"]
    if store["manifest_sha256"] != authority.PREF2_STORE_MANIFEST_SHA256:
        _stop("PREF1_PREF2_MANIFEST_DRIFT", store["manifest_sha256"])
    return str(store["manifest_sha256"])


@dataclass(frozen=True, slots=True)
class _Segment:
    left: np.ndarray
    left_rhs: np.ndarray
    right: np.ndarray
    right_rhs: np.ndarray
    width: float


def _journal_evidence(root: Path, replay: Mapping[str, Any]) -> Mapping[str, Any]:
    relative = f"{authority.PREF2_STORE_PATH}/journal/{int(replay['journal_sequence']):020d}-{replay['journal_sha256']}.journal"
    raw = _read_leaf(root, relative)
    if sha256(raw).hexdigest() != replay["journal_raw_sha256"]:
        _stop("PREF1_JOURNAL_HASH_DRIFT", replay["retry"])
    value = json.loads(raw, object_pairs_hook=_pairs)
    evidence = value.get("payload", {}).get("evidence")
    if not isinstance(evidence, Mapping):
        _stop("PREF1_JOURNAL_EVIDENCE_DRIFT", replay["retry"])
    return evidence


def _plan(cursor: p15.Proto15Cursor, replay: Mapping[str, Any]) -> Any:
    pending = cursor.payload.get("retry_successor_payload_or_none")
    if not isinstance(pending, Mapping):
        _stop("PREF1_RETRY_PLAN_ABSENT", replay["retry"])
    persisted = p15._plan_mapping(pending["successor_plan"])
    answer = lattice.plan_forward_proto14_subdivision(
        float.fromhex(str(persisted["current_hex"])),
        float.fromhex(str(persisted["event_target_hex"])),
        float.fromhex(str(persisted["requested_cap_hex"])),
        minimum_width=float.fromhex(str(persisted["minimum_width_hex_or_none"])),
    )
    if p15._plan_mapping(answer) != persisted or answer.macro_width.hex() != replay["attempted_width_hex"]:
        _stop("PREF1_RETRY_PLAN_DRIFT", replay["retry"])
    return answer


def _prepare(root: Path, store: HLT16CampaignStore, shells: Mapping[str, object], replay: Mapping[str, Any]) -> tuple[tdg6.TDG6PreparedGR0Compositor, Mapping[str, Any]]:
    generation = int(replay["predecessor_generation"])
    checkpoint = store.authenticated_checkpoint_at_generation(generation)
    if checkpoint.sha256 != replay["checkpoint_sha256"]:
        _stop("PREF1_CHECKPOINT_DRIFT", generation)
    relative = f"{authority.PREF2_STORE_PATH}/checkpoints/{generation:020d}-{checkpoint.sha256}.json"
    if sha256(_read_leaf(root, relative)).hexdigest() != replay["checkpoint_raw_sha256"]:
        _stop("PREF1_CHECKPOINT_RAW_DRIFT", generation)
    member_state = checkpoint.members[authority.MEMBER_KEY]
    member = deepcopy(shells[authority.MEMBER_KEY])
    try:
        persisted = store.load_state(member_state.descriptor_sha256)
        snapshot = HLT16MemberSnapshot(
            persisted.arrays,
            persisted.descriptor["metadata"],
        )
        restore_snapshot(member, snapshot)
        generation_identity = persisted.descriptor["metadata"]["runtime_identity"]
        if generation_identity["generation"] == 0:
            member.source_retry_count += member_state.source_total
            member.CFL_retry_count += member_state.cfl_total
        elif (
            member.source_retry_count < member_state.source_total
            or member.CFL_retry_count < member_state.cfl_total
        ):
            _stop("PREF1_RETRY_OVERLAY_DRIFT", replay["retry"])
        member.temporal_ledger = p15._ledger_from_mapping(member_state.ledger)
    except TDG9AR1PREF1Error:
        raise
    except Exception as exc:
        raise TDG9AR1PREF1Error(
            "PREF1_MEMBER_RESTORE_DRIFT", replay["retry"]
        ) from exc
    cursor = p15.Proto15Cursor(dict(member_state.cursor))
    cursor.validate()
    ledger = member.temporal_ledger
    if ledger is None or cursor.mode != "RETRY_PENDING" or member_state.pending_owner != "temporal" or ledger.current_macro_step_temporal_retry_count != replay["prior_retry_count"]:
        _stop("PREF1_PREDECESSOR_SELECTION_DRIFT", replay["retry"])
    plan = _plan(cursor, replay)
    prepared = tdg6.prepare_tdg6_gr0_compositor(
        method=member.integrator_id, time=plan.current, step_size=plan.macro_width,
        state=member.state, rhs=member.operator, projector=member.projector,
        transaction=member.transaction, tracers=member.tracers,
        coordinates=member.initial.grid.coordinates, temporal_ledger=ledger,
        previous_step_index=member.step_index,
        previous_transaction_serial=member.transaction_serial,
    )
    historical = _journal_evidence(root, replay)
    seen: list[Mapping[str, Any]] = []
    try:
        tdg6.require_tdg6_temporal_admission(
            prepared, transaction=member.transaction, tracers=member.tracers,
            temporal_ledger=ledger, current_time=member.time, current_state=member.state,
            current_step_index=member.step_index,
            current_transaction_serial=member.transaction_serial,
            durable_rejection_sink=lambda item: seen.append(dict(item)),
        )
    except tdg6.TDG6TemporalRetryRequired as exc:
        if exc.evidence.retry_count_for_current_macro_step != replay["retry"]:
            _stop("PREF1_RETRY_COUNT_DRIFT", replay["retry"])
    except tdg6.TDG6TemporalRetryExhausted as exc:
        _stop("PREF1_UNEXPECTED_EXHAUSTION", exc.reason)
    else:
        _stop("PREF1_HISTORICAL_REJECTION_ADMITTED", replay["retry"])
    if len(seen) != 1 or _canonical(p15._json_safe(seen[0])) != _canonical(historical):
        _stop("PREF1_HISTORICAL_EVIDENCE_DRIFT", replay["retry"])
    return prepared, historical


def _owned(value: object, label: str) -> np.ndarray:
    if type(value) is not np.ndarray:
        _stop("PREF1_ARRAY_TYPE_DRIFT", label)
    if value.dtype.str != "<f8" or value.shape != (3, authority.OWNED_ROW_COUNT, 6) or not value.flags.c_contiguous or not bool(np.isfinite(value).all()):
        _stop("PREF1_ARRAY_SHAPE_OR_VALUE_DRIFT", label)
    return value


def _segment(proposal: Any) -> _Segment:
    start_rhs, end_rhs = tdg5._proposal_endpoint_records(proposal)
    width = proposal.final_time - proposal.initial_time
    if type(width) is not float or not math.isfinite(width) or width <= 0.0:
        _stop("PREF1_WIDTH_DRIFT", width)
    return _Segment(
        _owned(tdg5._owned_state(proposal.initial_state), "left"),
        _owned(tdg5._owned_rhs(start_rhs), "left_rhs"),
        _owned(tdg5._owned_state(proposal.candidate_state), "right"),
        _owned(tdg5._owned_rhs(end_rhs), "right_rhs"), width,
    )


def _surface(prepared: tdg6.TDG6PreparedGR0Compositor) -> tuple[tuple[_Segment, ...], ...]:
    paths = tuple(tuple(attempt.proposal for attempt in path.attempts) for path in (prepared.outer, prepared.medium, prepared.fine))
    if tuple(map(len, paths)) != (1, 2, 4):
        _stop("PREF1_PROPOSAL_SHAPE_DRIFT", tuple(map(len, paths)))
    return tuple(tuple(_segment(item) for item in path) for path in paths)


def _rows(surface: tuple[tuple[_Segment, ...], ...], channel: str) -> Iterator[object]:
    block_name, field_name = channel.split(":", 1)
    block = _BLOCK_NAMES.index(block_name)
    field = _FIELD_NAMES.index(field_name)
    def segment(item: _Segment, row: int) -> tuple[float, ...]:
        values = (float(item.left[block,row,field]), float(item.left_rhs[block,row,field]), float(item.right[block,row,field]), float(item.right_rhs[block,row,field]), item.width)
        if any(type(value) is not float or not math.isfinite(value) for value in values):
            _stop("PREF1_SEGMENT_BINARY64_DRIFT", channel)
        return values
    outer, medium, fine = surface
    for row in range(authority.OWNED_ROW_COUNT):
        yield (segment(outer[0], row), tuple(segment(item,row) for item in medium), tuple(segment(item,row) for item in fine))


def _row_hash(rows: Iterable[object]) -> str:
    digest = sha256(b"TDG9-AR1-BINARY64-HERMITE-ROWS-v1\n")
    for ordinal, row in enumerate(rows):
        outer, medium, fine = row  # type: ignore[misc]
        flattened = [value for segment in (outer, *medium, *fine) for value in segment]
        if len(flattened) != 35 or any(type(value) is not float for value in flattened):
            _stop("PREF1_ROW_STREAM_DRIFT", ordinal)
        digest.update((f"{ordinal}|" + "|".join(value.hex() for value in flattened) + "\n").encode())
    return digest.hexdigest()


def _fraction(value: Fraction) -> dict[str, str]:
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


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
    _stop("PREF1_EXACT_EVIDENCE_DRIFT", type(value).__name__)


def _overlap(first: Any, second: Any) -> bool:
    return max(first.lower, second.lower) <= min(first.upper, second.upper)


def _assess(rows_factory: Any) -> dict[str, Any]:
    row_hash = _row_hash(rows_factory())
    primary = independent = primary_exhaustion = independent_exhaustion = None
    try:
        primary = exact_primary.assess_exact_temporal_refinement(rows_factory(), max_depth=authority.MAX_DEPTH, max_nodes=authority.MAX_NODES_PER_CHANNEL_PER_EVALUATOR)
    except exact_primary.ExactTemporalArithmeticResourceExhausted as exc:
        primary_exhaustion = exc.evidence
    try:
        independent = exact_independent.assess_exact_temporal_refinement_independently(rows_factory(), max_depth=authority.MAX_DEPTH, max_nodes=authority.MAX_NODES_PER_CHANNEL_PER_EVALUATOR)
    except exact_independent.IndependentExactTemporalArithmeticResourceExhausted as exc:
        independent_exhaustion = exc.evidence
    if primary is None or independent is None:
        return {"classification":"resource_inconclusive","row_stream_sha256":row_hash,"primary":None if primary is None else _exact_json(primary),"independent":None if independent is None else _exact_json(independent),"primary_resource_exhaustion":_exact_json(primary_exhaustion),"independent_resource_exhaustion":_exact_json(independent_exhaustion)}
    if primary.classification != independent.classification or primary.combined_coefficient_stream_sha256 != independent.combined_coefficient_stream_sha256 or not _overlap(primary.outer_difference.interval, independent.outer_difference.interval) or not _overlap(primary.finest_difference.interval, independent.finest_difference.interval):
        _stop("PREF1_EVALUATOR_DISAGREEMENT", row_hash)
    return {"classification":primary.classification,"row_stream_sha256":row_hash,"primary":_exact_json(primary),"independent":_exact_json(independent),"primary_resource_exhaustion":None,"independent_resource_exhaustion":None}


def _reduce(retries: Sequence[Mapping[str, Any]]) -> str:
    classes = [channel["classification"] for retry in retries for channel in retry["channels"]]
    if any(item == "sufficient_contraction_failure" for item in classes):
        return "legacy_enclosure_not_sole_owner_on_frozen_samples"
    if any(item == "resource_inconclusive" for item in classes):
        return "arithmetic_discriminator_inconclusive"
    if classes and all(item in {"exact_zero","sufficient_contraction_pass"} for item in classes):
        return "legacy_enclosure_owned_nonadmission_on_all_frozen_samples"
    _stop("PREF1_CLASSIFICATION_DRIFT", classes)


def expected_evidence(config_raw: bytes) -> dict[str, Any]:
    config = _config(config_raw)
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "target_protocol": authority.TARGET_PROTOCOL,
        "raw": dict(config["raw"]),
        "replay": dict(config["replay"]),
        "scope": dict(config["scope"]),
        "claims": dict(config["claims"]),
    }


def bind_raw_diagnostic(config_raw: bytes, repository: Path) -> dict[str, Any]:
    expected = expected_evidence(config_raw)
    root = Path(os.path.abspath(os.fspath(repository)))
    try:
        root_metadata = root.lstat()
    except OSError as exc:
        raise TDG9AR1PREF1Error("PREF1_REPOSITORY_ABSENT", root) from exc
    if stat.S_ISLNK(root_metadata.st_mode) or not stat.S_ISDIR(root_metadata.st_mode):
        _stop("PREF1_REPOSITORY_UNSAFE", root)
    manifest_raw, terminal_raw, manifest, terminal = _raw_snapshot(root)
    _authenticate_authority(root, manifest)
    if manifest.get("schema") != RAW_SCHEMA or manifest.get("artifact_id") != authority.ARTIFACT_ID or manifest.get("runner_id") != RUNNER_ID:
        _stop("PREF1_RAW_MANIFEST_DRIFT", "identity differs")
    before_store = _authenticate_pref2(root)
    if terminal.get("schema") != RAW_SCHEMA or terminal.get("authority_commit") != AUTHORITY_COMMIT or terminal.get("evaluator_calls") != 108:
        _stop("PREF1_RAW_TERMINAL_DRIFT", "identity/call count differs")
    if any(terminal.get(key) is not value for key,value in {"store_unchanged":True,"campaign_store_mutated":False,"PDE_state_committed":False,"continuation_authorized":False,"fourth_width_executed":False,"automatic_resource_escalation_used":False,"candidate_branch_opened":False,"physical_result_earned":False}.items()):
        _stop("PREF1_CLAIM_PROMOTION", "raw terminal crossed scope")
    store = HLT16CampaignStore(root / authority.PREF2_STORE_PATH)
    shells = build_static_gr0_shells(root)
    replayed: list[dict[str, Any]] = []
    calls = 0
    raw_retries = terminal.get("retries")
    if not isinstance(raw_retries, list) or len(raw_retries) != 3:
        _stop("PREF1_RAW_RETRY_SHAPE", type(raw_retries).__name__)
    for replay, raw_retry in zip(authority.REPLAYS, raw_retries, strict=True):
        prepared, historical = _prepare(root, store, shells, replay)
        surface = _surface(prepared)
        channels: list[dict[str, Any]] = []
        raw_channels = raw_retry.get("channels")
        if not isinstance(raw_channels, list) or len(raw_channels) != 18:
            _stop("PREF1_RAW_CHANNEL_SHAPE", replay["retry"])
        for channel, raw_channel in zip(authority.CHANNEL_ORDER, raw_channels, strict=True):
            result = _assess(lambda s=surface, c=channel: _rows(s, c))
            calls += 2
            full = {"channel":channel,"owned_row_count":authority.OWNED_ROW_COUNT,**result}
            if full != raw_channel:
                _stop("PREF1_RAW_CHANNEL_EVIDENCE_DRIFT", f"retry {replay['retry']} {channel}")
            channels.append(full)
        replay_result = {"retry":replay["retry"],"attempted_width_hex":replay["attempted_width_hex"],"historical_TDG6_evidence_sha256":sha256(_canonical(historical)).hexdigest(),"historical_TDG6_evidence_reproduced_byte_identically":True,"channels":channels}
        if replay_result != raw_retry:
            _stop("PREF1_RAW_RETRY_EVIDENCE_DRIFT", replay["retry"])
        replayed.append(replay_result)
    classification = _reduce(replayed)
    if calls != authority.MAX_EVALUATOR_CALLS or classification != terminal.get("classification") or classification != RAW_CLASSIFICATION:
        _stop("PREF1_AGGREGATE_DRIFT", {"calls":calls,"classification":classification})
    after_store = _authenticate_pref2(root)
    manifest_raw_after, terminal_raw_after, _, _ = _raw_snapshot(root)
    if before_store != after_store or manifest_raw != manifest_raw_after or terminal_raw != terminal_raw_after:
        _stop("PREF1_INPUT_MUTATED_DURING_BIND", "store/raw changed")
    compact_retries = []
    for retry in replayed:
        compact_channels = []
        for channel in retry["channels"]:
            compact_channels.append({
                "channel": channel["channel"],
                "classification": channel["classification"],
                "row_stream_sha256": channel["row_stream_sha256"],
                "primary_evidence_sha256": sha256(_canonical(channel["primary"] if channel["primary"] is not None else channel["primary_resource_exhaustion"])).hexdigest(),
                "independent_evidence_sha256": sha256(_canonical(channel["independent"] if channel["independent"] is not None else channel["independent_resource_exhaustion"])).hexdigest(),
            })
        compact_retries.append({
            "retry": retry["retry"], "attempted_width_hex": retry["attempted_width_hex"],
            "historical_TDG6_evidence_sha256": retry["historical_TDG6_evidence_sha256"],
            "failure_channels": [item["channel"] for item in retry["channels"] if item["classification"] == "sufficient_contraction_failure"],
            "channels": compact_channels,
        })
    return {
        **expected,
        "raw_binding": {"manifest_sha256":sha256(manifest_raw).hexdigest(),"terminal_sha256":sha256(terminal_raw).hexdigest(),"leaf_count":2,"raw_namespace_unchanged":True},
        "store_binding": {"manifest_before":before_store,"manifest_after":after_store,"store_unchanged":True},
        "independent_replay": {"evaluator_calls":calls,"classification":classification,"retries":compact_retries},
        "conclusion": {
            "legacy_binary64_enclosure_is_not_the_sole_owner_on_the_three_frozen_samples": True,
            "exact_discrete_failures_exist": True,
            "failed_channel_count": sum(len(item["failure_channels"]) for item in compact_retries),
            "successor_selected": False,
            "physics_inference_permitted": False,
        },
    }


def build_pref1_result(config_raw: bytes, repository: Path, *, live: bool = True) -> dict[str, Any]:
    payload = bind_raw_diagnostic(config_raw, repository) if live else expected_evidence(config_raw)
    return {"artifact_id": ARTIFACT_ID, "artifact_payload": payload}


def validate_compact_result(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    expected = expected_evidence(config_raw)
    if sha256(result_raw).hexdigest() != COMPACT_RESULT_SHA256:
        _stop("PREF1_COMPACT_DRIFT", "compact result hash differs")
    result = _json(result_raw, "compact result")
    if result.get("artifact_id") != ARTIFACT_ID or not isinstance(result.get("artifact_payload"), Mapping):
        _stop("PREF1_COMPACT_DRIFT", "compact identity differs")
    payload = result["artifact_payload"]
    for key,value in expected.items():
        if payload.get(key) != value:
            _stop("PREF1_COMPACT_DRIFT", key)
    replay = payload.get("independent_replay")
    if not isinstance(replay, Mapping) or replay.get("evaluator_calls") != 108 or replay.get("classification") != RAW_CLASSIFICATION:
        _stop("PREF1_COMPACT_DRIFT", "replay summary differs")
    retries = replay.get("retries")
    if not isinstance(retries, list) or [item.get("retry") for item in retries] != [3,4,5]:
        _stop("PREF1_COMPACT_DRIFT", "retry summary differs")
    for item in retries:
        if tuple(item.get("failure_channels", ())) != _EXPECTED_FAILURES[item["retry"]] or len(item.get("channels", ())) != 18:
            _stop("PREF1_COMPACT_DRIFT", f"retry {item.get('retry')}")
    conclusion = payload.get("conclusion")
    if not isinstance(conclusion, Mapping) or conclusion != {
        "exact_discrete_failures_exist":True,
        "failed_channel_count":10,
        "legacy_binary64_enclosure_is_not_the_sole_owner_on_the_three_frozen_samples":True,
        "physics_inference_permitted":False,
        "successor_selected":False,
    }:
        _stop("PREF1_COMPACT_DRIFT", "conclusion differs")
    return result


__all__ = (
    "TDG9AR1PREF1Error", "bind_raw_diagnostic", "build_pref1_result",
    "canonical_result", "expected_evidence", "validate_compact_result",
)
