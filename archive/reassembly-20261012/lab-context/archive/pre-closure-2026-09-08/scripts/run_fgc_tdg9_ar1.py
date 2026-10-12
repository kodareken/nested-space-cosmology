#!/usr/bin/env python3
"""Run exactly one read-only TDG9 audit over RCV3 retries 3, 4, and 5.

The RCV3 store is authenticated before and after deterministic replay and is
never opened for writing.  The only durable mutation is an ignored AR1 result
namespace, installed from a private sibling directory with exclusive rename
after all three widths have terminated.  No PDE state is committed.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import asdict, dataclass, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import stat
import sys
import tempfile
from typing import Any, Iterable, Iterator, Mapping, NoReturn, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
_SOURCE_ROOT = str(ROOT / "src")
sys.path[:] = [entry for entry in sys.path if entry != _SOURCE_ROOT]
sys.path.insert(0, _SOURCE_ROOT)

from recursive_horizons.fgc.evolution import hlt16_campaign_runtime as campaign_runtime  # noqa: E402
from recursive_horizons.fgc.evolution import proto15_runtime as p15  # noqa: E402
from recursive_horizons.fgc.evolution import tdg5_stage_complete_refinement_runtime as tdg5  # noqa: E402
from recursive_horizons.fgc.evolution import tdg6_temporal_admission_runtime as tdg6  # noqa: E402
from recursive_horizons.fgc.evolution import tdg7_binary64_subdivision_lattice as lattice  # noqa: E402
from recursive_horizons.fgc.evolution import tdg8_persisted_retry_replay as tdg8_replay  # noqa: E402
from recursive_horizons.fgc.evolution import tdg8_rcv3_pref2_binder as pref2  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ar1_authority as authority  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_exact_temporal_arithmetic as exact_primary  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_exact_temporal_arithmetic_independent as exact_independent  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_store import HLT16CampaignStore  # noqa: E402
from recursive_horizons.fgc.evolution.proto17_hlt15_runtime import _rename_exclusive  # noqa: E402
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import build_static_gr0_shells  # noqa: E402


RUNNER_ID = "FGC-1-TDG9-AR1-RUN1"
RUNNER_SCHEMA = "FGC-1-TDG9-AR1-read-only-runner-v1"
RAW_RESULT_SCHEMA = "FGC-1-TDG9-AR1-raw-result-v1"
_MAX_LEAF_BYTES = 16 * 1024 * 1024
_MAX_DETAIL = 640
_FIELD_NAMES = ("alpha", "v", "lambda", "R", "phi", "chi")
_BLOCK_NAMES = ("u", "p", "q")


@dataclass(frozen=True, slots=True)
class _HermiteSegmentSurface:
    left: Any
    left_rhs: Any
    right: Any
    right_rhs: Any
    width: float


class TDG9AR1RunnerError(RuntimeError):
    """An authority, replay, arithmetic, or output premise differs."""

    def __init__(self, owner: str, code: str, detail: object) -> None:
        self.owner = _safe_text(owner)
        self.code = _safe_text(code)
        self.detail = _safe_text(detail)
        super().__init__(f"{self.owner}/{self.code}: {self.detail}")

    def as_record(self) -> dict[str, object]:
        return {
            "schema": RUNNER_SCHEMA,
            "runner_id": RUNNER_ID,
            "classification": "invalid_provenance_or_implementation",
            "error_owner": self.owner,
            "error_code": self.code,
            "error_detail": self.detail,
            "campaign_store_mutated": False,
            "PDE_state_committed": False,
            "continuation_authorized": False,
            "candidate_branch_opened": False,
            "physical_result_earned": False,
        }


def _safe_text(value: object) -> str:
    text = " ".join(str(value).replace(str(ROOT), "<repo>").split())
    return text.encode("ascii", "backslashreplace").decode("ascii")[:_MAX_DETAIL] or "none"


def _fail(owner: str, code: str, detail: object) -> NoReturn:
    raise TDG9AR1RunnerError(owner, code, detail)


def _canonical(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def _pretty(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def _unique(pairs: list[tuple[str, object]]) -> dict[str, object]:
    answer: dict[str, object] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError("duplicate JSON key")
        answer[key] = value
    return answer


def _fixed_root(root: Path) -> Path:
    supplied = Path(os.path.abspath(os.fspath(root)))
    try:
        metadata = supplied.lstat()
    except OSError as exc:
        raise TDG9AR1RunnerError("provenance", "repository_absent", exc) from exc
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        _fail("provenance", "repository_unsafe", "repository is not a real directory")
    if supplied.resolve() != ROOT.resolve():
        _fail("authority", "repository_scope", "runner requires its fixed repository")
    return supplied


def _read_leaf(root: Path, relative: str, maximum: int = _MAX_LEAF_BYTES) -> bytes:
    candidate = Path(relative)
    if candidate.is_absolute() or any(part in {"", ".", ".."} for part in candidate.parts):
        _fail("provenance", "leaf_path_unsafe", relative)
    path = root / candidate
    try:
        before = path.lstat()
    except OSError as exc:
        raise TDG9AR1RunnerError("provenance", "leaf_absent", relative) from exc
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode) or before.st_size > maximum:
        _fail("provenance", "leaf_unsafe", relative)
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        active = os.fstat(descriptor)
        identity = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
        if identity != (active.st_dev, active.st_ino, active.st_size, active.st_mtime_ns):
            _fail("provenance", "leaf_raced", relative)
        chunks: list[bytes] = []
        size = 0
        while size <= maximum:
            block = os.read(descriptor, min(1 << 20, maximum + 1 - size))
            if not block:
                break
            chunks.append(block)
            size += len(block)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if identity != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
        _fail("provenance", "leaf_changed", relative)
    raw = b"".join(chunks)
    if len(raw) != before.st_size or len(raw) > maximum:
        _fail("provenance", "leaf_size", relative)
    return raw


def _execution_authority(root: Path, commit: str) -> authority.TDG9AR1Authority:
    try:
        config = _read_leaf(root, authority.CONFIG_PATH, 2 * 1024 * 1024)
        result = _read_leaf(root, authority.RESULT_PATH, 8 * 1024 * 1024)
        receipt = authority.authorize_diagnostic(root, config, result, commit)
    except TDG9AR1RunnerError:
        raise
    except Exception as exc:
        raise TDG9AR1RunnerError("authority", "committed_image_rejected", exc) from exc
    if (
        receipt.pref2_commit != authority.PREF2_COMMIT
        or receipt.evaluator_commit != authority.EVALUATOR_COMMIT
        or receipt.retries != (3, 4, 5)
        or receipt.retry_width_hex != tuple(
            item["attempted_width_hex"] for item in authority.REPLAYS
        )
        or receipt.channel_order != authority.CHANNEL_ORDER
        or receipt.environment != tuple(authority.ENVIRONMENT.items())
        or receipt.selection != tuple(authority.SELECTION.items())
        or receipt.authority_delta_paths != authority.AUTHORITY_DELTA_PATHS
        or receipt.max_evaluator_calls != 108
        or not receipt.diagnostic_execution_authorized
        or receipt.campaign_state_mutation_authorized
        or receipt.continuation_authorized
        or receipt.candidate_branches_authorized
    ):
        _fail("authority", "typed_scope_rejected", "AR1 typed authority differs")
    return receipt


def _authenticate_pref2(root: Path) -> str:
    config = _read_leaf(root, authority.PREF2_CONFIG_PATH, 2 * 1024 * 1024)
    result = _read_leaf(root, authority.PREF2_RESULT_PATH, 8 * 1024 * 1024)
    if sha256(config).hexdigest() != authority.PREF2_CONFIG_SHA256:
        _fail("provenance", "pref2_config_drift", "PREF2 config differs")
    if sha256(result).hexdigest() != authority.PREF2_RESULT_SHA256:
        _fail("provenance", "pref2_result_drift", "PREF2 result differs")
    try:
        regenerated = pref2.build_pref2_result(config, root)
    except Exception as exc:
        raise TDG9AR1RunnerError("provenance", "pref2_live_rejected", exc) from exc
    if pref2.canonical_result(regenerated) != result:
        _fail("provenance", "pref2_live_drift", "live terminal differs")
    payload = regenerated.get("artifact_payload")
    if not isinstance(payload, Mapping):
        _fail("provenance", "pref2_payload", "PREF2 payload differs")
    terminal_evidence = payload.get("terminal_evidence")
    store = (
        terminal_evidence.get("store")
        if isinstance(terminal_evidence, Mapping)
        else None
    )
    if not isinstance(store, Mapping) or store.get("manifest_sha256") != authority.PREF2_STORE_MANIFEST_SHA256:
        _fail("provenance", "pref2_manifest", "PREF2 store manifest differs")
    return authority.PREF2_STORE_MANIFEST_SHA256


def _output_path(root: Path) -> Path:
    return root / authority.OUTPUT_NAMESPACE


def _inspect_output(root: Path) -> dict[str, object]:
    target = _output_path(root)
    current = root
    for part in Path(authority.OUTPUT_NAMESPACE).parent.parts:
        current = current / part
        try:
            metadata = current.lstat()
        except FileNotFoundError:
            return {"state": "absent", "output_namespace": authority.OUTPUT_NAMESPACE}
        except OSError as exc:
            raise TDG9AR1RunnerError("output", "namespace_unreadable", exc) from exc
        if stat.S_ISLNK(metadata.st_mode):
            _fail("output", "namespace_symlink", current)
        if not stat.S_ISDIR(metadata.st_mode):
            _fail("output", "namespace_parent_not_directory", current)
    try:
        with os.scandir(target.parent) as entries:
            parent_names = tuple(sorted(item.name for item in entries))
    except OSError as exc:
        raise TDG9AR1RunnerError("output", "namespace_parent_unreadable", exc) from exc
    stages = tuple(
        name for name in parent_names if name.startswith(authority.STAGING_PREFIX)
    )
    if stages:
        _fail("output", "private_stage_present", stages)
    try:
        target_metadata = target.lstat()
    except FileNotFoundError:
        return {"state": "absent", "output_namespace": authority.OUTPUT_NAMESPACE}
    except OSError as exc:
        raise TDG9AR1RunnerError("output", "namespace_unreadable", exc) from exc
    if stat.S_ISLNK(target_metadata.st_mode):
        _fail("output", "namespace_symlink", target)
    if not stat.S_ISDIR(target_metadata.st_mode):
        _fail("output", "namespace_not_directory", target)
    try:
        with os.scandir(target) as entries:
            observed = tuple(sorted(entries, key=lambda item: item.name))
    except OSError as exc:
        raise TDG9AR1RunnerError("output", "namespace_unreadable", exc) from exc
    names = tuple(item.name for item in observed)
    if names != ("manifest.json", "terminal.json"):
        _fail("output", "namespace_partial_or_foreign", names)
    if any(
        item.is_symlink() or not item.is_file(follow_symlinks=False)
        for item in observed
    ):
        _fail("output", "namespace_leaf_unsafe", names)
    manifest = _read_leaf(root, f"{authority.OUTPUT_NAMESPACE}/manifest.json")
    terminal = _read_leaf(root, f"{authority.OUTPUT_NAMESPACE}/terminal.json")
    try:
        manifest_value = json.loads(manifest, object_pairs_hook=_unique)
        terminal_value = json.loads(terminal, object_pairs_hook=_unique)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TDG9AR1RunnerError("output", "namespace_invalid_json", exc) from exc
    if _pretty(manifest_value) != manifest or _pretty(terminal_value) != terminal:
        _fail("output", "namespace_noncanonical", "AR1 output differs")
    return {
        "state": "terminal",
        "output_namespace": authority.OUTPUT_NAMESPACE,
        "manifest_sha256": sha256(manifest).hexdigest(),
        "terminal_sha256": sha256(terminal).hexdigest(),
        "classification": terminal_value.get("classification"),
    }


def status(root: Path, *, authority_commit: str) -> dict[str, object]:
    repository = _fixed_root(root)
    receipt = _execution_authority(repository, authority_commit)
    store_manifest = _authenticate_pref2(repository)
    output = _inspect_output(repository)
    return {
        "schema": RUNNER_SCHEMA,
        "runner_id": RUNNER_ID,
        "authority_commit": receipt.authority_commit,
        "store_manifest_sha256": store_manifest,
        "PREF2_terminal_authenticated": True,
        "status_read_only": True,
        "campaign_store_mutated": False,
        "PDE_state_committed": False,
        "candidate_branch_opened": False,
        **output,
    }


def _journal_evidence(root: Path, replay: Mapping[str, object]) -> Mapping[str, object]:
    sequence = int(replay["journal_sequence"])
    digest = str(replay["journal_sha256"])
    relative = (
        f"{authority.PREF2_STORE_PATH}/journal/"
        f"{sequence:020d}-{digest}.journal"
    )
    raw = _read_leaf(root, relative)
    if sha256(raw).hexdigest() != replay["journal_raw_sha256"]:
        _fail("replay", "journal_raw_hash", f"retry {replay['retry']}")
    try:
        value = json.loads(raw, object_pairs_hook=_unique)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TDG9AR1RunnerError("replay", "journal_decode", exc) from exc
    if (
        value.get("record_sha256") != digest
        or value.get("sequence") != sequence
        or value.get("kind") != "tdg6_rejection"
    ):
        _fail("replay", "journal_identity", f"retry {replay['retry']}")
    payload = value.get("payload")
    evidence = payload.get("evidence") if isinstance(payload, Mapping) else None
    if not isinstance(evidence, Mapping):
        _fail("replay", "journal_evidence_absent", f"retry {replay['retry']}")
    if (
        evidence.get("retry_count_for_current_macro_step") != replay["retry"]
        or float(evidence.get("attempted_macro_step_size", 0.0)).hex()
        != replay["attempted_width_hex"]
        or evidence.get("initial_state_sha256") != authority.PREF2_PHYSICAL_STATE_SHA256
    ):
        _fail("replay", "journal_evidence_identity", f"retry {replay['retry']}")
    return evidence


def _typed_plan(cursor: p15.Proto15Cursor, replay: Mapping[str, object]) -> Any:
    pending = cursor.payload.get("retry_successor_payload_or_none")
    if not isinstance(pending, Mapping):
        _fail("replay", "pending_plan_absent", f"retry {replay['retry']}")
    try:
        persisted = p15._plan_mapping(pending["successor_plan"])
        plan = lattice.plan_forward_proto14_subdivision(
            float.fromhex(str(persisted["current_hex"])),
            float.fromhex(str(persisted["event_target_hex"])),
            float.fromhex(str(persisted["requested_cap_hex"])),
            minimum_width=float.fromhex(str(persisted["minimum_width_hex_or_none"])),
        )
    except Exception as exc:
        raise TDG9AR1RunnerError("replay", "pending_plan_invalid", exc) from exc
    if p15._plan_mapping(plan) != persisted:
        _fail("replay", "pending_plan_reconstruction", f"retry {replay['retry']}")
    if plan.macro_width.hex() != replay["attempted_width_hex"]:
        _fail("replay", "attempted_width_drift", f"retry {replay['retry']}")
    return plan


def _require_checkpoint_selection(
    store: HLT16CampaignStore,
    checkpoint: Any,
    member_state: Any,
) -> Mapping[str, object]:
    target = {
        "binary64_hex": authority.TARGET_BINARY64_HEX,
        "rational": authority.TARGET_RATIONAL,
    }
    if (
        checkpoint.protocol != authority.TARGET_PROTOCOL
        or checkpoint.campaign_id != authority.CAMPAIGN_ID
        or checkpoint.plan_sha256 != authority.CAMPAIGN_PLAN_SHA256
        or checkpoint.authorization_commit != authority.CAMPAIGN_AUTHORIZATION_COMMIT
        or checkpoint.event != authority.EVENT
        or checkpoint.target != target
        or authority.MEMBER_KEY not in checkpoint.members
    ):
        _fail("selection", "checkpoint_selection", checkpoint.generation)
    cursor = member_state.cursor
    if not isinstance(cursor, Mapping) or (
        cursor.get("protocol_artifact_id") != authority.TARGET_PROTOCOL
        or cursor.get("campaign_id") != authority.CAMPAIGN_ID
        or cursor.get("member_key") != authority.MEMBER_KEY
        or cursor.get("method") != authority.METHOD
        or cursor.get("point_count") != authority.POINT_COUNT
        or cursor.get("committed_common_event_index") != authority.EVENT
        or cursor.get("event_target_time") != target
        or cursor.get("accepted_state_sha256") != member_state.descriptor_sha256
        or cursor.get("mode") != "RETRY_PENDING"
    ):
        _fail("selection", "cursor_selection", checkpoint.generation)
    try:
        persisted = store.load_state(member_state.descriptor_sha256)
        descriptor = persisted.descriptor
        metadata = descriptor["metadata"]
        runtime_identity = metadata["runtime_identity"]
        template = metadata["template"]
    except Exception as exc:
        raise TDG9AR1RunnerError("selection", "descriptor_selection", exc) from exc
    if not all(
        isinstance(item, Mapping)
        for item in (descriptor, metadata, runtime_identity, template)
    ) or (
        runtime_identity.get("protocol_artifact_id") != authority.TARGET_PROTOCOL
        or runtime_identity.get("campaign_id") != authority.CAMPAIGN_ID
        or runtime_identity.get("branch") != authority.BRANCH
        or runtime_identity.get("amplitude") != authority.AMPLITUDE
        or runtime_identity.get("member_key") != authority.MEMBER_KEY
        or runtime_identity.get("method") != authority.METHOD
        or runtime_identity.get("point_count") != authority.POINT_COUNT
        or template.get("amplitude") != authority.AMPLITUDE
        or template.get("method") != authority.METHOD
        or template.get("point_count") != authority.POINT_COUNT
    ):
        _fail("selection", "descriptor_selection", checkpoint.generation)
    return runtime_identity


def _require_replayed_evidence(
    replayed: object,
    historical: Mapping[str, object],
    replay: Mapping[str, object],
    *,
    rejected_retry_count: object,
) -> Mapping[str, object]:
    try:
        normalized = p15._json_safe(replayed)
    except Exception as exc:
        raise TDG9AR1RunnerError("replay", "replayed_evidence_invalid", exc) from exc
    if not isinstance(normalized, Mapping):
        _fail("replay", "replayed_evidence_invalid", f"retry {replay['retry']}")
    replayed_bytes = _canonical(normalized)
    historical_bytes = _canonical(historical)
    if (
        normalized != historical
        or replayed_bytes != historical_bytes
        or rejected_retry_count != replay["retry"]
    ):
        _fail("replay", "historical_evidence_mismatch", f"retry {replay['retry']}")
    return normalized


def _prepare_replay(
    root: Path,
    store: HLT16CampaignStore,
    shells: Mapping[str, object],
    replay: Mapping[str, object],
) -> tuple[Any, Mapping[str, object]]:
    generation = int(replay["predecessor_generation"])
    try:
        checkpoint = store.authenticated_checkpoint_at_generation(generation)
    except Exception as exc:
        raise TDG9AR1RunnerError("replay", "checkpoint_unavailable", exc) from exc
    if checkpoint.sha256 != replay["checkpoint_sha256"]:
        _fail("replay", "checkpoint_content_hash", f"generation {generation}")
    checkpoint_path = (
        f"{authority.PREF2_STORE_PATH}/checkpoints/"
        f"{generation:020d}-{checkpoint.sha256}.json"
    )
    if sha256(_read_leaf(root, checkpoint_path)).hexdigest() != replay["checkpoint_raw_sha256"]:
        _fail("replay", "checkpoint_raw_hash", f"generation {generation}")
    member_state = checkpoint.members.get(authority.MEMBER_KEY)
    if member_state is None:
        _fail("replay", "member_absent", authority.MEMBER_KEY)
    _require_checkpoint_selection(store, checkpoint, member_state)
    member = deepcopy(shells[authority.MEMBER_KEY])
    try:
        campaign_runtime.restore_member_with_overlay(
            store, checkpoint, member, key=authority.MEMBER_KEY
        )
        cursor = p15.Proto15Cursor(dict(member_state.cursor))
        cursor.validate()
    except Exception as exc:
        raise TDG9AR1RunnerError("replay", "member_restore", exc) from exc
    ledger = member.temporal_ledger
    if (
        ledger is None
        or member.key != authority.MEMBER_KEY
        or member.amplitude != authority.AMPLITUDE
        or member.method_label != authority.METHOD
        or member.point_count != authority.POINT_COUNT
        or ledger.current_macro_step_temporal_retry_count != replay["prior_retry_count"]
        or cursor.mode != "RETRY_PENDING"
        or member_state.pending_owner != "temporal"
    ):
        _fail("replay", "predecessor_retry_depth", f"retry {replay['retry']}")
    plan = _typed_plan(cursor, replay)
    try:
        prepared = tdg8_replay.prepare_tdg8_persisted_retry_replay(
            replay_plan=plan,
            expected_prior_retry_count=int(replay["prior_retry_count"]),
            method=member.integrator_id,
            state=member.state,
            rhs=member.operator,
            projector=member.projector,
            transaction=member.transaction,
            tracers=member.tracers,
            coordinates=member.initial.grid.coordinates,
            temporal_ledger=ledger,
            previous_step_index=member.step_index,
            previous_transaction_serial=member.transaction_serial,
        )
    except Exception as exc:
        raise TDG9AR1RunnerError("replay", "proposal_reconstruction", exc) from exc
    historical = _journal_evidence(root, replay)
    seen: list[Mapping[str, object]] = []
    try:
        tdg8_replay.require_tdg8_persisted_retry_replay_admission(
            prepared,
            transaction=member.transaction,
            tracers=member.tracers,
            temporal_ledger=ledger,
            current_time=member.time,
            current_state=member.state,
            current_step_index=member.step_index,
            current_transaction_serial=member.transaction_serial,
            durable_rejection_sink=lambda item: seen.append(dict(item)),
        )
    except tdg6.TDG6TemporalRetryRequired as rejected:
        if len(seen) != 1:
            _fail("replay", "rejection_count", f"retry {replay['retry']}")
        _require_replayed_evidence(
            seen[0], historical, replay,
            rejected_retry_count=rejected.evidence.retry_count_for_current_macro_step,
        )
    except tdg6.TDG6TemporalRetryExhausted as exc:
        raise TDG9AR1RunnerError("replay", "unexpected_exhaustion", exc) from exc
    else:
        _fail("replay", "historical_rejection_admitted", f"retry {replay['retry']}")
    return prepared, historical


def _proposals(prepared: Any) -> tuple[tuple[Any, ...], tuple[Any, ...], tuple[Any, ...]]:
    try:
        inner = prepared._prepared.prepared_tdg6
        paths = tuple(
            tuple(attempt.proposal for attempt in path.attempts)
            for path in (inner.outer, inner.medium, inner.fine)
        )
    except Exception as exc:
        raise TDG9AR1RunnerError("replay", "prepared_surface", exc) from exc
    if tuple(len(path) for path in paths) != (1, 2, 4):
        _fail("replay", "prepared_path_count", tuple(len(path) for path in paths))
    return paths  # type: ignore[return-value]


def _require_owned_array(value: object, *, label: str) -> np.ndarray:
    if type(value) is not np.ndarray:
        _fail("replay", "Hermite_owned_array_type", label)
    array = value
    if (
        array.dtype.str != "<f8"
        or array.shape != (3, authority.OWNED_ROW_COUNT, 6)
        or not array.flags.c_contiguous
    ):
        _fail(
            "replay",
            "Hermite_owned_array_layout",
            {
                "label": label,
                "dtype": array.dtype.str,
                "shape": array.shape,
                "C_contiguous": bool(array.flags.c_contiguous),
            },
        )
    if not bool(np.isfinite(array).all()):
        _fail("replay", "Hermite_owned_array_nonfinite", label)
    return array


def _segment_surface(proposal: Any) -> _HermiteSegmentSurface:
    try:
        start_rhs, end_rhs = tdg5._proposal_endpoint_records(proposal)
        left = _require_owned_array(
            tdg5._owned_state(proposal.initial_state), label="left_state"
        )
        left_rhs = _require_owned_array(
            tdg5._owned_rhs(start_rhs), label="left_rhs"
        )
        right = _require_owned_array(
            tdg5._owned_state(proposal.candidate_state), label="right_state"
        )
        right_rhs = _require_owned_array(
            tdg5._owned_rhs(end_rhs), label="right_rhs"
        )
        width = proposal.final_time - proposal.initial_time
    except Exception as exc:
        if isinstance(exc, TDG9AR1RunnerError):
            raise
        raise TDG9AR1RunnerError("replay", "Hermite_endpoint_record", exc) from exc
    if (
        type(width) is not float
        or not math.isfinite(width)
        or width <= 0.0
    ):
        _fail("replay", "Hermite_width", {"type": type(width).__name__, "width": width})
    return _HermiteSegmentSurface(left, left_rhs, right, right_rhs, width)


def _hermite_surface(
    prepared: Any,
) -> tuple[
    tuple[_HermiteSegmentSurface, ...],
    tuple[_HermiteSegmentSurface, ...],
    tuple[_HermiteSegmentSurface, ...],
]:
    paths = _proposals(prepared)
    return tuple(
        tuple(_segment_surface(proposal) for proposal in path) for path in paths
    )  # type: ignore[return-value]


def _segment(
    surface: _HermiteSegmentSurface,
    block: int,
    row: int,
    field: int,
) -> tuple[float, ...]:
    answer = (
        float(surface.left[block, row, field]),
        float(surface.left_rhs[block, row, field]),
        float(surface.right[block, row, field]),
        float(surface.right_rhs[block, row, field]),
        surface.width,
    )
    if any(type(value) is not float or not math.isfinite(value) for value in answer):
        _fail("replay", "Hermite_binary64_type", "segment differs")
    return answer


def _rows_factory(
    surface: tuple[
        tuple[_HermiteSegmentSurface, ...],
        tuple[_HermiteSegmentSurface, ...],
        tuple[_HermiteSegmentSurface, ...],
    ],
    channel: str,
) -> tuple[int, Any]:
    if channel not in authority.CHANNEL_ORDER:
        _fail("arithmetic", "channel_unknown", channel)
    block_name, field_name = channel.split(":", 1)
    block = _BLOCK_NAMES.index(block_name)
    field = _FIELD_NAMES.index(field_name)
    outer, medium, fine = surface
    row_count = int(outer[0].left.shape[1])
    if row_count != authority.OWNED_ROW_COUNT:
        _fail("replay", "owned_row_count", row_count)

    def rows() -> Iterator[object]:
        for row in range(row_count):
            yield (
                _segment(outer[0], block, row, field),
                tuple(_segment(item, block, row, field) for item in medium),
                tuple(_segment(item, block, row, field) for item in fine),
            )

    return row_count, rows


def _row_stream_hash(rows: Iterable[object]) -> str:
    digest = sha256(b"TDG9-AR1-BINARY64-HERMITE-ROWS-v1\n")
    for ordinal, row in enumerate(rows):
        outer, medium, fine = row  # type: ignore[misc]
        flattened: list[float] = []
        for segment in (outer, *medium, *fine):
            if not isinstance(segment, Sequence) or len(segment) != 5:
                _fail("arithmetic", "row_shape", ordinal)
            for value in segment:
                if type(value) is not float:
                    _fail("arithmetic", "row_not_binary64", ordinal)
                flattened.append(value)
        digest.update(
            (f"{ordinal}|" + "|".join(value.hex() for value in flattened) + "\n").encode("ascii")
        )
    return digest.hexdigest()


def _fraction(value: Fraction) -> dict[str, str]:
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def _json_exact(value: object) -> object:
    if isinstance(value, Fraction):
        return _fraction(value)
    if is_dataclass(value) and not isinstance(value, type):
        return _json_exact(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _json_exact(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_exact(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    _fail("arithmetic", "evidence_not_exact_json", type(value).__name__)


def _overlaps(first: Any, second: Any) -> bool:
    return max(first.lower, second.lower) <= min(first.upper, second.upper)


def _assess_channel(rows_factory: Any) -> dict[str, object]:
    try:
        row_hash = _row_stream_hash(rows_factory())
    except TDG9AR1RunnerError:
        raise
    except Exception as exc:
        raise TDG9AR1RunnerError("arithmetic", "row_stream_invalid", exc) from exc
    primary: Any = None
    independent: Any = None
    primary_exhaustion: Any = None
    independent_exhaustion: Any = None
    try:
        primary = exact_primary.assess_exact_temporal_refinement(
            rows_factory(),
            max_depth=authority.MAX_DEPTH,
            max_nodes=authority.MAX_NODES_PER_CHANNEL_PER_EVALUATOR,
        )
    except exact_primary.ExactTemporalArithmeticResourceExhausted as exc:
        primary_exhaustion = exc.evidence
    except Exception as exc:
        raise TDG9AR1RunnerError("arithmetic", "primary_evaluator_invalid", exc) from exc
    try:
        independent = exact_independent.assess_exact_temporal_refinement_independently(
            rows_factory(),
            max_depth=authority.MAX_DEPTH,
            max_nodes=authority.MAX_NODES_PER_CHANNEL_PER_EVALUATOR,
        )
    except exact_independent.IndependentExactTemporalArithmeticResourceExhausted as exc:
        independent_exhaustion = exc.evidence
    except Exception as exc:
        raise TDG9AR1RunnerError("arithmetic", "independent_evaluator_invalid", exc) from exc

    if primary is None or independent is None:
        return {
            "classification": "resource_inconclusive",
            "row_stream_sha256": row_hash,
            "primary": None if primary is None else _json_exact(primary),
            "independent": None if independent is None else _json_exact(independent),
            "primary_resource_exhaustion": _json_exact(primary_exhaustion),
            "independent_resource_exhaustion": _json_exact(independent_exhaustion),
        }
    if (
        primary.row_count != independent.row_count
        or primary.classification != independent.classification
        or primary.outer_difference.coefficient_stream_sha256
        != independent.outer_difference.coefficient_stream_sha256
        or primary.finest_difference.coefficient_stream_sha256
        != independent.finest_difference.coefficient_stream_sha256
        or primary.combined_coefficient_stream_sha256
        != independent.combined_coefficient_stream_sha256
        or not _overlaps(primary.outer_difference.interval, independent.outer_difference.interval)
        or not _overlaps(primary.finest_difference.interval, independent.finest_difference.interval)
    ):
        _fail("arithmetic", "evaluator_disagreement", "completed evaluators differ")
    return {
        "classification": primary.classification,
        "row_stream_sha256": row_hash,
        "primary": _json_exact(primary),
        "independent": _json_exact(independent),
        "primary_resource_exhaustion": None,
        "independent_resource_exhaustion": None,
    }


def _terminal_classification(retries: Sequence[Mapping[str, object]]) -> str:
    if len(retries) != len(authority.REPLAYS) or any(
        not isinstance(retry.get("channels"), Sequence)
        or isinstance(retry.get("channels"), (str, bytes))
        or len(retry["channels"]) != len(authority.CHANNEL_ORDER)  # type: ignore[arg-type]
        for retry in retries
    ):
        _fail("arithmetic", "classification_shape", "retry/channel count differs")
    channels = [
        channel
        for retry in retries
        for channel in retry["channels"]  # type: ignore[index]
    ]
    classifications = [item["classification"] for item in channels]
    if any(value == "sufficient_contraction_failure" for value in classifications):
        return "legacy_enclosure_not_sole_owner_on_frozen_samples"
    if any(value == "resource_inconclusive" for value in classifications):
        return "arithmetic_discriminator_inconclusive"
    if classifications and all(value == "exact_zero" for value in classifications):
        # Exact zero is not called a p>=3/2 pass.  It is mapped directly to the
        # same aggregate owner result because a zero discrete difference cannot
        # be the nonadmission obstruction on these finite replay cubics.
        return "legacy_enclosure_owned_nonadmission_on_all_frozen_samples"
    if all(value in {"exact_zero", "sufficient_contraction_pass"} for value in classifications):
        return "legacy_enclosure_owned_nonadmission_on_all_frozen_samples"
    _fail("arithmetic", "classification_reduction", classifications)


def _manifest(receipt: authority.TDG9AR1Authority) -> dict[str, object]:
    return {
        "schema": RAW_RESULT_SCHEMA,
        "artifact_id": authority.ARTIFACT_ID,
        "runner_id": RUNNER_ID,
        "authority_commit": receipt.authority_commit,
        "PREF2_commit": receipt.pref2_commit,
        "evaluator_commit": receipt.evaluator_commit,
        "environment": dict(receipt.environment),
        "selection": dict(receipt.selection),
        "authority_delta_paths": list(receipt.authority_delta_paths),
        "store_manifest_sha256": receipt.store_manifest_sha256,
        "output_namespace": receipt.output_namespace,
        "replays": [dict(item) for item in authority.REPLAYS],
        "channel_order": list(receipt.channel_order),
        "max_depth": receipt.max_depth,
        "max_nodes_per_channel_per_evaluator": receipt.max_nodes_per_channel_per_evaluator,
        "max_evaluator_calls": receipt.max_evaluator_calls,
        "PDE_state_commit_authorized": False,
        "campaign_store_mutation_authorized": False,
        "continuation_authorized": False,
        "candidate_branches_authorized": False,
    }


def _ensure_parent(root: Path) -> Path:
    relative = Path(authority.OUTPUT_NAMESPACE).parent
    current = root
    for part in relative.parts:
        current = current / part
        try:
            metadata = current.lstat()
        except FileNotFoundError:
            current.mkdir(mode=0o755)
            metadata = current.lstat()
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
            _fail("output", "parent_unsafe", current)
    return current


def _write_fsynced(path: Path, raw: bytes) -> None:
    descriptor = os.open(
        path,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
        0o644,
    )
    try:
        view = memoryview(raw)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                _fail("output", "write_failed", path.name)
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _publish(root: Path, manifest: Mapping[str, object], terminal: Mapping[str, object]) -> None:
    parent = _ensure_parent(root)
    target = _output_path(root)
    if target.exists() or target.is_symlink():
        _fail("output", "namespace_preexisting", authority.OUTPUT_NAMESPACE)
    staging = Path(tempfile.mkdtemp(prefix=authority.STAGING_PREFIX, dir=parent))
    _write_fsynced(staging / "manifest.json", _pretty(manifest))
    _write_fsynced(staging / "terminal.json", _pretty(terminal))
    directory_fd = os.open(staging, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
    try:
        _rename_exclusive(staging, target)
    except Exception as exc:
        raise TDG9AR1RunnerError("output", "exclusive_publication_failed", exc) from exc
    parent_fd = os.open(parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(parent_fd)
    finally:
        os.close(parent_fd)


def run(root: Path, *, authority_commit: str) -> dict[str, object]:
    repository = _fixed_root(root)
    receipt = _execution_authority(repository, authority_commit)
    if _inspect_output(repository)["state"] != "absent":
        _fail("output", "namespace_preexisting", authority.OUTPUT_NAMESPACE)
    before_manifest = _authenticate_pref2(repository)
    store = HLT16CampaignStore(repository / authority.PREF2_STORE_PATH)
    try:
        shells = build_static_gr0_shells(repository)
    except Exception as exc:
        raise TDG9AR1RunnerError("replay", "static_shells", exc) from exc
    if tuple(shells) != (
        "RK4-2049", "RK4-4097", "RK4-8193",
        "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385",
    ):
        _fail("replay", "static_shell_order", tuple(shells))
    retry_results: list[dict[str, object]] = []
    evaluator_calls = 0
    for replay in authority.REPLAYS:
        prepared, historical = _prepare_replay(repository, store, shells, replay)
        surface = _hermite_surface(prepared)
        channels: list[dict[str, object]] = []
        for channel in authority.CHANNEL_ORDER:
            row_count, rows_factory = _rows_factory(surface, channel)
            result = _assess_channel(rows_factory)
            evaluator_calls += 2
            channels.append({"channel": channel, "owned_row_count": row_count, **result})
        retry_results.append(
            {
                "retry": replay["retry"],
                "attempted_width_hex": replay["attempted_width_hex"],
                "historical_TDG6_evidence_sha256": sha256(_canonical(historical)).hexdigest(),
                "historical_TDG6_evidence_reproduced_byte_identically": True,
                "channels": channels,
            }
        )
    if evaluator_calls != authority.MAX_EVALUATOR_CALLS:
        _fail("arithmetic", "evaluator_call_count", evaluator_calls)
    after_manifest = _authenticate_pref2(repository)
    if before_manifest != after_manifest:
        _fail("provenance", "store_changed_during_replay", "PREF2 manifest differs")
    classification = _terminal_classification(retry_results)
    manifest = _manifest(receipt)
    terminal = {
        "schema": RAW_RESULT_SCHEMA,
        "artifact_id": authority.ARTIFACT_ID,
        "runner_id": RUNNER_ID,
        "authority_commit": receipt.authority_commit,
        "classification": classification,
        "store_manifest_before": before_manifest,
        "store_manifest_after": after_manifest,
        "store_unchanged": True,
        "evaluator_calls": evaluator_calls,
        "retries": retry_results,
        "campaign_store_mutated": False,
        "PDE_state_committed": False,
        "continuation_authorized": False,
        "fourth_width_executed": False,
        "automatic_resource_escalation_used": False,
        "candidate_branch_opened": False,
        "physical_result_earned": False,
    }
    # A foreign namespace appearing during the read-only calculation cannot be
    # adopted or overwritten.  The exclusive rename owns the final decision.
    _publish(repository, manifest, terminal)
    return {
        "schema": RUNNER_SCHEMA,
        "runner_id": RUNNER_ID,
        "state": "terminal",
        "classification": classification,
        "output_namespace": authority.OUTPUT_NAMESPACE,
        "evaluator_calls": evaluator_calls,
        "campaign_store_mutated": False,
        "PDE_state_committed": False,
        "continuation_authorized": False,
        "candidate_branch_opened": False,
        "physical_result_earned": False,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authority-commit", required=True)
    parser.add_argument("--status", action="store_true")
    arguments = parser.parse_args(argv)
    try:
        result = (
            status(ROOT, authority_commit=arguments.authority_commit)
            if arguments.status
            else run(ROOT, authority_commit=arguments.authority_commit)
        )
    except TDG9AR1RunnerError as exc:
        print(json.dumps(exc.as_record(), sort_keys=True))
        return 2
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
