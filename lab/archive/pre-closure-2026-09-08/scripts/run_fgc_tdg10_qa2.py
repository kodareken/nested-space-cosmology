#!/usr/bin/env python3
"""Run one two-width SSPRK3 exact complete-C robustness diagnostic."""

from __future__ import annotations

import argparse
import ctypes
import errno
from hashlib import sha256
import importlib
import json
import os
from pathlib import Path
import secrets
import stat
import sys
from typing import Mapping, NoReturn, Sequence

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from scripts import run_fgc_tdg9_ti2 as ti2_runner  # noqa: E402
from scripts import run_fgc_tdg10_qa1 as qa1_runner  # noqa: E402
from recursive_horizons.fgc.evolution import tdg6_temporal_admission_runtime as tdg6  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ar1_authority as ar1  # noqa: E402
from recursive_horizons.fgc.evolution import tdg10_qa2_authority as authority  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_store import (  # noqa: E402
    HLT16CampaignStore,
)
from recursive_horizons.fgc.evolution.numerical_engine import COMPARATOR_METHOD  # noqa: E402
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import (  # noqa: E402
    build_static_gr0_shells,
)
from recursive_horizons.fgc.evolution.tdg10_exact_complete_c_admission import (  # noqa: E402
    ExactCompleteCResourceExhausted,
    ExactCompleteCRouteDisagreement,
)
from recursive_horizons.fgc.evolution.tdg10_exact_complete_c_runtime import (  # noqa: E402
    TDG10ExactCompleteCRuntimeClosed,
)


RUNNER_ID = "FGC-1-TDG10-QA2-RUN1"
RAW_SCHEMA = "FGC-1-TDG10-QA2-raw-v1"
EXACT_COMPLETE_C_RUNTIME_MODULE = (
    "recursive_horizons.fgc.evolution.tdg10_exact_complete_c_runtime"
)
EXACT_COMPLETE_C_ASSESS_NAME = "assess_tdg10_exact_complete_c"
INCONCLUSIVE_EXCEPTIONS = (
    ExactCompleteCResourceExhausted,
    ExactCompleteCRouteDisagreement,
    TDG10ExactCompleteCRuntimeClosed,
)
BOTH_PASS_CLASS = "completed_both_width_all_channel_pass_no_state_advance"
NONPASS_CLASS = "completed_one_or_more_width_all_channel_nonpass_no_state_advance"
INCONCLUSIVE_CLASS = "exact_route_or_resource_inconclusive_no_state_advance"
PREMISE_STOP_CLASS = "shadow_or_scientific_premise_stop_no_state_advance"
INVALID_CLASS = "invalid_provenance_or_implementation"
WIDTH_PASS_CLASS = "all_18_channel_pass"
WIDTH_NONPASS_CLASS = "one_or_more_channel_nonpass"
TERMINAL_CLASSES = frozenset(
    {
        BOTH_PASS_CLASS,
        NONPASS_CLASS,
        INCONCLUSIVE_CLASS,
        PREMISE_STOP_CLASS,
        INVALID_CLASS,
    }
)
PUBLISHABLE_TERMINAL_CLASSES = TERMINAL_CLASSES - {INVALID_CLASS}
_OUTPUT_LEAF_LIMIT = 512 * 1024 * 1024
_RENAME_EXCL = 0x00000004
_FORBIDDEN_OUTPUT_NAMES = frozenset(
    {
        "checkpoints",
        "journal",
        "members",
        "cursor",
        "ledger",
        "campaign",
        "common-event",
        "common_event",
        "payloads",
        "fine-endpoint",
    }
)


class QA2RunnerError(RuntimeError):
    """Typed fail-closed QA2 runner error."""

    def __init__(self, owner: str, code: str, detail: object) -> None:
        self.owner = str(owner)
        self.code = str(code)
        self.detail = " ".join(str(detail).replace(str(ROOT), "<repo>").split())[:640]
        super().__init__(f"{self.owner}/{self.code}: {self.detail}")


def _fail(owner: str, code: str, detail: object) -> NoReturn:
    raise QA2RunnerError(owner, code, detail)


def _wrap_qa1(exc: Exception) -> NoReturn:
    if isinstance(exc, qa1_runner.QA1RunnerError):
        raise QA2RunnerError(exc.owner, exc.code, exc.detail) from exc
    raise QA2RunnerError("qa1", "adapter", exc) from exc


def _canonical(value: object) -> bytes:
    return qa1_runner._canonical(value)


def _pretty(value: object) -> bytes:
    return qa1_runner._pretty(value)


def _unique(items: list[tuple[str, object]]) -> dict[str, object]:
    try:
        return qa1_runner._unique(items)
    except qa1_runner.QA1RunnerError as exc:
        _wrap_qa1(exc)


def _require_keys(
    value: object, expected: set[str], *, label: str
) -> Mapping[str, object]:
    try:
        return qa1_runner._require_keys(value, expected, label=label)
    except qa1_runner.QA1RunnerError as exc:
        _wrap_qa1(exc)


def _nonclaims() -> dict[str, object]:
    return {
        "diagnostic_qualification_only": True,
        "diagnostic_fine_endpoint_serialized": False,
        "diagnostic_fine_endpoint_is_accepted_state": False,
        "new_protocol_id_authorized": False,
        "target_protocol": authority.TARGET_PROTOCOL,
        "production_SSPRK3_comparator_earned": False,
        "independent_method_agreement_earned": False,
        "production_method_earned": False,
        "state_advance_authorized": False,
        "campaign_state_write_authorized": False,
        "PDE_state_commit_authorized": False,
        "PDE_state_committed": False,
        "fine_path_commit_authorized": False,
        "fine_path_committed": False,
        "temporal_retry_admission_called": False,
        "candidate_execution_authorized": False,
        "candidate_branch_opened": False,
        "mechanism_result_authorized": False,
        "physical_result_authorized": False,
        "common_event_authorized": False,
        "GR0_calibration_authorized": False,
        "SGBL_authorized": False,
        "FGCQR_authorized": False,
        "retained_EFT_evolution_authorized": False,
        "physical_transition_claim_authorized": False,
        "external_publication_authorized": False,
        "push_authorized": False,
        "retry_3_authorized": False,
        "fourth_width_authorized": False,
        "old_member_adoption_authorized": False,
        "RA1_authorized": False,
        "successor_remedy_selected": False,
        "common_event_completed": False,
        "GR0_calibration_completed": False,
        "mechanism_result_earned": False,
        "physical_result_earned": False,
        "QA1_licenses_QA2_method_design_never_old_member_adoption": True,
    }


def _execution_authority(root: Path, commit: str) -> authority.QA2Authority:
    try:
        return authority.authorize(
            root,
            authority.read_leaf(root, authority.CONFIG_PATH),
            authority.read_leaf(root, authority.RESULT_PATH),
            commit,
        )
    except Exception as exc:
        raise QA2RunnerError("authority", "rejected", exc) from exc


def _snapshot_store(root: Path) -> tuple[int, str]:
    try:
        return ti2_runner._snapshot_store(root)
    except Exception as exc:
        raise QA2RunnerError("provenance", "store_snapshot", exc) from exc


def _restore_replay(
    root: Path, replay: Mapping[str, object]
) -> ti2_runner.RestoredReplay:
    width = authority.width_spec(int(replay["retry"]))
    if int(replay["predecessor_generation"]) == authority.FORBIDDEN_RETRY3_GENERATION:
        _fail("replay", "must_not_restore_generation_9", 9)
    if int(replay["predecessor_generation"]) == width[
        "forbidden_predecessor_generation"
    ]:
        _fail(
            "replay",
            "forbidden_predecessor_generation",
            width["forbidden_predecessor_generation"],
        )
    if int(replay["predecessor_generation"]) != width["generation"]:
        _fail("replay", "generation_mismatch", replay["predecessor_generation"])
    if int(replay["journal_sequence"]) != width["historical_rejection_sequence"]:
        _fail(
            "replay",
            "historical_rejection_mismatch",
            replay["journal_sequence"],
        )
    store = HLT16CampaignStore(root / ar1.PREF2_STORE_PATH)
    shells = build_static_gr0_shells(root)
    try:
        return ti2_runner._restore_replay(root, store, shells, replay)
    except ti2_runner.TI1RunnerError as exc:
        raise QA2RunnerError(exc.owner, exc.code, exc.detail) from exc


def _prepare_shadow(
    restored: ti2_runner.RestoredReplay,
) -> tdg6.TDG6PreparedGR0Compositor:
    try:
        return ti2_runner._prepare_shadow(restored)
    except tdg6.TDG6RefinementPathStop:
        raise
    except ti2_runner.TI1RunnerError as exc:
        raise QA2RunnerError(exc.owner, exc.code, exc.detail) from exc


def _assess_prepared_exact_complete_c(prepared: object) -> object:
    try:
        module = importlib.import_module(EXACT_COMPLETE_C_RUNTIME_MODULE)
        assess = getattr(module, EXACT_COMPLETE_C_ASSESS_NAME)
    except (ModuleNotFoundError, AttributeError) as exc:
        raise QA2RunnerError("exact", "adapter_unavailable", exc) from exc
    return assess(
        prepared,
        maximum_candidates_D01=authority.PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
        maximum_candidates_D12=authority.PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
        refinement_depth=authority.PRIMARY_REFINEMENT_DEPTH,
        expected_owned_row_count=authority.OWNED_ROW_COUNT,
    )


def serialize_exact_assessment(assessment: object) -> dict[str, object]:
    try:
        serialized = qa1_runner.serialize_exact_assessment(assessment)
    except qa1_runner.QA1RunnerError as exc:
        _wrap_qa1(exc)
    if (
        serialized.get("channel_order") != list(authority.CHANNEL_ORDER)
        or serialized.get("channel_count") != authority.CHANNEL_COUNT
        or serialized.get("interval_owner") != authority.INTERVAL_OWNER
        or serialized.get("admission_is_all_of") is not True
        or serialized.get("independent_route_required") is not True
    ):
        _fail("exact", "serialized_identity", serialized.get("channel_count"))
    return serialized


def _replay_receipt(restored: ti2_runner.RestoredReplay) -> dict[str, object]:
    fingerprint = restored.fingerprint
    if not isinstance(fingerprint, Mapping):
        _fail("replay", "fingerprint_type", type(fingerprint).__name__)
    retry = int(restored.replay["retry"])
    width = authority.width_spec(retry)
    receipt = {
        "member_key": fingerprint.get("member_key"),
        "retry": width["retry"],
        "predecessor_generation": width["generation"],
        "checkpoint_journal_sequence": width["checkpoint_journal_sequence"],
        "historical_rejection_sequence": width["historical_rejection_sequence"],
        "attempted_width_hex": width["attempted_width_hex"],
        "accepted_time_hex": fingerprint.get("accepted_time_hex"),
        "state_sha256": fingerprint.get("state_sha256"),
        "descriptor_sha256": fingerprint.get("descriptor_sha256"),
        "transaction_sha256": fingerprint.get("transaction_sha256"),
        "historical_journal_sha256": restored.historical_journal_sha256,
        "ledger_owner": authority.LEDGER_OWNER,
        "tracer_owner": authority.TRACER_OWNER,
        "transaction_owner": authority.TRANSACTION_OWNER,
    }
    expected = authority.expected_raw_replay_receipt(retry)
    if receipt != expected:
        _fail("replay", "receipt_identity", receipt)
    return receipt


def reduce_qa2_terminal(width_records: Sequence[Mapping[str, object]]) -> str:
    if len(width_records) != 2:
        _fail("reduction", "width_count", len(width_records))
    retries = tuple(item.get("retry") for item in width_records)
    if retries != authority.SELECTED_RETRIES:
        _fail("reduction", "selected_retries", retries)
    classes = tuple(item.get("width_class") for item in width_records)
    allowed = {
        WIDTH_PASS_CLASS,
        WIDTH_NONPASS_CLASS,
        "inconclusive",
        "premise_stop",
    }
    unknown = tuple(item for item in classes if item not in allowed)
    if unknown:
        _fail("reduction", "unknown_width_class", unknown)
    if any(item == "premise_stop" for item in classes):
        return PREMISE_STOP_CLASS
    if any(item == "inconclusive" for item in classes):
        return INCONCLUSIVE_CLASS
    if any(item != WIDTH_PASS_CLASS for item in classes):
        return NONPASS_CLASS
    return BOTH_PASS_CLASS


def _terminal_base(
    authority_commit: str,
    classification: str,
    before: tuple[int, str],
    after: tuple[int, str],
    *,
    width_robustness_passed: bool,
) -> dict[str, object]:
    if classification not in TERMINAL_CLASSES:
        _fail("reduction", "terminal_class", classification)
    if width_robustness_passed is True and classification != BOTH_PASS_CLASS:
        _fail("reduction", "robustness_without_both_pass", classification)
    if classification == BOTH_PASS_CLASS and width_robustness_passed is not True:
        _fail("reduction", "both_pass_without_robustness", classification)
    return {
        "schema": RAW_SCHEMA,
        "artifact_id": authority.ARTIFACT_ID,
        "runner_id": RUNNER_ID,
        "classification": classification,
        "authority_commit": authority_commit,
        "tableau_selector": authority.TABLEAU_SELECTOR,
        "tableau_runtime_selector": COMPARATOR_METHOD,
        "actual_spatial_operator": authority.ACTUAL_SPATIAL_OPERATOR,
        "selected_retries": [4, 5],
        "channel_order": list(authority.CHANNEL_ORDER),
        "width_robustness_passed": width_robustness_passed,
        "store_snapshot_before": {"leaf_count": before[0], "sha256": before[1]},
        "store_snapshot_after": {"leaf_count": after[0], "sha256": after[1]},
        "store_unchanged": before == after,
        **_nonclaims(),
    }


def _manifest(authority_commit: str) -> dict[str, object]:
    return {
        "schema": RAW_SCHEMA,
        "artifact_id": authority.ARTIFACT_ID,
        "runner_id": RUNNER_ID,
        "authority_commit": authority_commit,
        "tableau_selector": authority.TABLEAU_SELECTOR,
        "actual_spatial_operator": authority.ACTUAL_SPATIAL_OPERATOR,
        "selected_retries": [4, 5],
        "target_protocol": authority.TARGET_PROTOCOL,
        "diagnostic_qualification_only": True,
        "diagnostic_fine_endpoint_serialized": False,
        "state_advance_authorized": False,
        "width_robustness_passed_requires_both_widths_all_18_pass": True,
        "output_leaves": ["manifest.json", "terminal.json"],
    }


def _shadow_budget() -> dict[str, int]:
    return {
        "shadow_path_count": 7,
        "shadow_proposal_count": 7,
        "SSPRK3_stage_and_endpoint_record_count": 28,
    }


def _open_output_parent(root: Path) -> tuple[int, str]:
    relative = Path(authority.OUTPUT_NAMESPACE)
    descriptor = os.open(
        root,
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        for part in relative.parent.parts:
            try:
                child = os.open(
                    part,
                    os.O_RDONLY
                    | getattr(os, "O_DIRECTORY", 0)
                    | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=descriptor,
                )
            except FileNotFoundError:
                os.mkdir(part, mode=0o755, dir_fd=descriptor)
                child = os.open(
                    part,
                    os.O_RDONLY
                    | getattr(os, "O_DIRECTORY", 0)
                    | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=descriptor,
                )
            os.close(descriptor)
            descriptor = child
        return descriptor, relative.name
    except Exception:
        os.close(descriptor)
        raise


def _write_leaf(directory_fd: int, name: str, raw: bytes) -> None:
    leaf = os.open(
        name,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
        0o644,
        dir_fd=directory_fd,
    )
    try:
        pending = memoryview(raw)
        while pending:
            written = os.write(leaf, pending)
            if written <= 0:
                _fail("output", "short_write", name)
            pending = pending[written:]
        os.fsync(leaf)
    finally:
        os.close(leaf)


def _cleanup_stage(
    stage_fd: int | None,
    stage_name: str,
    parent_fd: int,
) -> None:
    if stage_fd is None:
        return
    for name in ("manifest.json", "terminal.json"):
        try:
            os.unlink(name, dir_fd=stage_fd)
        except FileNotFoundError:
            pass
    os.close(stage_fd)
    try:
        os.rmdir(stage_name, dir_fd=parent_fd)
    except FileNotFoundError:
        pass


def _publish(
    root: Path,
    manifest: Mapping[str, object],
    terminal: Mapping[str, object],
) -> None:
    parent_fd, target_name = _open_output_parent(root)
    stage_name = f"{authority.STAGING_PREFIX}{os.getpid()}-{secrets.token_hex(8)}"
    stage_fd: int | None = None
    try:
        try:
            os.stat(target_name, dir_fd=parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            _fail("output", "namespace_exists", authority.OUTPUT_NAMESPACE)
        os.mkdir(stage_name, mode=0o700, dir_fd=parent_fd)
        stage_fd = os.open(
            stage_name,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent_fd,
        )
        _write_leaf(stage_fd, "manifest.json", _pretty(manifest))
        _write_leaf(stage_fd, "terminal.json", _pretty(terminal))
        os.fsync(stage_fd)
        if os.uname().sysname != "Darwin":
            _fail("output", "exclusive_adoption_unavailable", os.uname().sysname)
        rename = getattr(ctypes.CDLL(None, use_errno=True), "renameatx_np", None)
        if rename is None:
            _fail("output", "exclusive_adoption_unavailable", "renameatx_np")
        rename.argtypes = (
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        )
        rename.restype = ctypes.c_int
        if rename(
            parent_fd,
            os.fsencode(stage_name),
            parent_fd,
            os.fsencode(target_name),
            _RENAME_EXCL,
        ):
            error = ctypes.get_errno()
            if error in {errno.EEXIST, errno.ENOTEMPTY}:
                _fail("output", "namespace_arrived", authority.OUTPUT_NAMESPACE)
            _fail("output", "exclusive_adoption_failed", error)
        os.close(stage_fd)
        stage_fd = None
        os.fsync(parent_fd)
    except Exception:
        _cleanup_stage(stage_fd, stage_name, parent_fd)
        stage_fd = None
        raise
    finally:
        os.close(parent_fd)


def _publish_terminal(
    root: Path,
    manifest: Mapping[str, object],
    terminal: Mapping[str, object],
    *,
    store_before: tuple[int, str],
    store_after: tuple[int, str],
) -> None:
    if store_after != store_before:
        _fail("provenance", "campaign_store_mutated", (store_before, store_after))
    classification = terminal.get("classification")
    if classification not in PUBLISHABLE_TERMINAL_CLASSES:
        _fail("output", "invalid_not_published", classification)
    if terminal.get("diagnostic_fine_endpoint_serialized") is not False:
        _fail("output", "diagnostic_endpoint_serialized", classification)
    _publish(root, manifest, terminal)


def _identity(value: os.stat_result) -> tuple[int, int, int, int, int, int]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
        value.st_nlink,
    )


def _read_output_leaf(directory_fd: int, name: str, before: os.stat_result) -> bytes:
    if (
        stat.S_ISLNK(before.st_mode)
        or not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or before.st_size < 0
        or before.st_size > _OUTPUT_LEAF_LIMIT
    ):
        _fail("status", "unsafe_output_leaf", name)
    descriptor = -1
    try:
        descriptor = os.open(
            name,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=directory_fd,
        )
        active = os.fstat(descriptor)
        expected = _identity(before)
        if _identity(active) != expected:
            _fail("status", "output_leaf_substituted", name)
        chunks: list[bytes] = []
        remaining = before.st_size
        while remaining:
            block = os.read(descriptor, min(1 << 20, remaining))
            if not block:
                _fail("status", "output_leaf_short_read", name)
            chunks.append(block)
            remaining -= len(block)
        if os.read(descriptor, 1):
            _fail("status", "output_leaf_grew", name)
        if _identity(os.fstat(descriptor)) != expected:
            _fail("status", "output_leaf_changed", name)
        after = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        if _identity(after) != expected:
            _fail("status", "output_leaf_substituted", name)
        return b"".join(chunks)
    except OSError as exc:
        raise QA2RunnerError("status", "unsafe_output_leaf", name) from exc
    finally:
        if descriptor != -1:
            os.close(descriptor)


def _parse_output_json(raw: bytes, *, name: str) -> dict[str, object]:
    def reject_constant(token: str) -> NoReturn:
        _fail("status", "nonfinite_json_constant", (name, token))

    try:
        value = json.loads(
            raw,
            object_pairs_hook=_unique,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise QA2RunnerError("status", "noncanonical_output", name) from exc
    if not isinstance(value, dict) or raw != _pretty(value):
        _fail("status", "noncanonical_output", name)
    return value


def _snapshot_output(root: Path) -> tuple[dict[str, object], dict[str, bytes], dict[str, str]] | None:
    flags = (
        os.O_RDONLY
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    root_before = root.lstat()
    if stat.S_ISLNK(root_before.st_mode) or not stat.S_ISDIR(root_before.st_mode):
        _fail("status", "unsafe_repository_root", root)
    descriptors: list[int] = []
    identities: list[tuple[int, int, int, int, int, int]] = []
    try:
        current = os.open(root, flags)
        descriptors.append(current)
        identities.append(_identity(os.fstat(current)))
        if identities[0] != _identity(root_before):
            _fail("status", "repository_root_substituted", root)
        for part in Path(authority.OUTPUT_NAMESPACE).parts:
            try:
                before = os.stat(part, dir_fd=current, follow_symlinks=False)
            except FileNotFoundError:
                return None
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                _fail("status", "unsafe_namespace", part)
            child = os.open(part, flags, dir_fd=current)
            if _identity(os.fstat(child)) != _identity(before):
                os.close(child)
                _fail("status", "output_directory_substituted", part)
            descriptors.append(child)
            identities.append(_identity(before))
            current = child
        output_fd = descriptors[-1]
        names = tuple(sorted(os.listdir(output_fd)))
        if any(name in _FORBIDDEN_OUTPUT_NAMES for name in names):
            _fail("status", "campaign_artifact_present", names)
        values: dict[str, object] = {}
        blobs: dict[str, bytes] = {}
        hashes: dict[str, str] = {}
        json_names = ("manifest.json", "terminal.json")
        if tuple(names) != json_names:
            _fail("status", "partial_or_foreign_namespace", names)
        for name in names:
            before = os.stat(name, dir_fd=output_fd, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
                _fail("status", "unsafe_namespace", name)
            raw = _read_output_leaf(output_fd, name, before)
            blobs[name] = raw
            hashes[name] = sha256(raw).hexdigest()
            values[name] = _parse_output_json(raw, name=name)
        for index, (descriptor, expected) in enumerate(
            zip(descriptors, identities, strict=True)
        ):
            if _identity(os.fstat(descriptor)) != expected:
                _fail("status", "output_directory_changed", index)
        if _identity(root.lstat()) != identities[0]:
            _fail("status", "repository_root_substituted", root)
        return values, blobs, hashes
    except OSError as exc:
        raise QA2RunnerError(
            "status", "output_snapshot_unreadable", authority.OUTPUT_NAMESPACE
        ) from exc
    finally:
        for descriptor in reversed(descriptors):
            os.close(descriptor)


def _validate_replay_receipt(value: object, *, retry: int) -> dict[str, object]:
    expected = authority.expected_raw_replay_receipt(retry)
    if not isinstance(value, Mapping) or set(value) != set(expected):
        _fail("status", "replay_receipt_schema", (sorted(expected), value))
    receipt = {key: value[key] for key in expected}
    if receipt != expected:
        _fail("status", "replay_receipt_identity", receipt)
    return receipt


def _validate_exact_summary(value: object) -> Mapping[str, object]:
    try:
        return qa1_runner._validate_exact_summary(value)
    except qa1_runner.QA1RunnerError as exc:
        _wrap_qa1(exc)


def _validate_width_record(value: object) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        _fail("status", "width_record_type", type(value).__name__)
    retry = value.get("retry")
    if retry not in authority.SELECTED_RETRIES:
        _fail("status", "width_retry", retry)
    width_class = value.get("width_class")
    base_keys = {
        "retry",
        "width_class",
        "predecessor_generation",
        "attempted_width_hex",
        "replay_receipt",
    }
    if width_class == "premise_stop":
        _require_keys(value, base_keys | {"typed_stop"}, label="premise width")
    elif width_class == "inconclusive":
        _require_keys(
            value,
            base_keys
            | {
                "typed_inconclusive",
                "shadow_path_count",
                "shadow_proposal_count",
                "SSPRK3_stage_and_endpoint_record_count",
            },
            label="inconclusive width",
        )
        if (
            value.get("shadow_path_count") != 7
            or value.get("shadow_proposal_count") != 7
            or value.get("SSPRK3_stage_and_endpoint_record_count") != 28
        ):
            _fail("status", "width_shadow_budget", retry)
    elif width_class in {WIDTH_PASS_CLASS, WIDTH_NONPASS_CLASS}:
        _require_keys(
            value,
            base_keys
            | {
                "shadow_path_count",
                "shadow_proposal_count",
                "SSPRK3_stage_and_endpoint_record_count",
                "exact_complete_C",
                "complete_admission_passed",
            },
            label="completed width",
        )
        exact = _validate_exact_summary(value["exact_complete_C"])
        passed = exact["complete_admission_passed"] is True and tuple(
            exact["failed_channels"]
        ) == ()
        expected_class = WIDTH_PASS_CLASS if passed else WIDTH_NONPASS_CLASS
        if (
            width_class != expected_class
            or value.get("complete_admission_passed") is not passed
            or value.get("shadow_path_count") != 7
            or value.get("shadow_proposal_count") != 7
            or value.get("SSPRK3_stage_and_endpoint_record_count") != 28
        ):
            _fail("status", "width_class_identity", retry)
    else:
        _fail("status", "width_class", width_class)
    spec = authority.width_spec(int(retry))
    if (
        value.get("predecessor_generation") != spec["generation"]
        or value.get("attempted_width_hex") != spec["attempted_width_hex"]
    ):
        _fail("status", "width_identity", retry)
    _validate_replay_receipt(value["replay_receipt"], retry=int(retry))
    return value


def _validate_terminal(
    value: object, *, authority_commit: str
) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        _fail("status", "terminal_schema", "terminal is not a mapping")
    classification = value.get("classification")
    if classification not in PUBLISHABLE_TERMINAL_CLASSES:
        _fail("status", "terminal_class_not_publishable", classification)
    sealed = (
        authority.SEALED_STORE_LEAF_COUNT,
        authority.SEALED_STORE_SNAPSHOT_SHA256,
    )
    robustness = value.get("width_robustness_passed") is True
    expected_base = _terminal_base(
        authority_commit,
        str(classification),
        sealed,
        sealed,
        width_robustness_passed=robustness,
    )
    if any(value.get(name) != expected for name, expected in expected_base.items()):
        _fail("status", "terminal_identity", classification)
    widths = value.get("widths")
    if not isinstance(widths, list) or len(widths) != 2:
        _fail("status", "width_records", type(widths).__name__)
    records = [_validate_width_record(item) for item in widths]
    expected_class = reduce_qa2_terminal(records)
    if classification != expected_class:
        _fail("status", "terminal_reduction", (classification, expected_class))
    extra = {"widths"}
    if classification == PREMISE_STOP_CLASS:
        extra.add("typed_stop")
    if classification == INCONCLUSIVE_CLASS:
        extra.add("typed_inconclusive")
    _require_keys(value, set(expected_base) | extra, label="qa2 terminal")
    if classification == PREMISE_STOP_CLASS:
        stop = _require_keys(
            value["typed_stop"], {"type", "detail", "retry"}, label="premise stop"
        )
        if (
            not isinstance(stop["type"], str)
            or not stop["type"]
            or not isinstance(stop["detail"], str)
            or not stop["detail"]
            or stop["retry"] not in authority.SELECTED_RETRIES
        ):
            _fail("status", "premise_stop_terminal", stop)
    if classification == INCONCLUSIVE_CLASS:
        closed = _require_keys(
            value["typed_inconclusive"],
            {"type", "reason", "detail", "retry"},
            label="typed inconclusive",
        )
        if (
            closed["type"]
            not in {
                "ExactCompleteCResourceExhausted",
                "ExactCompleteCRouteDisagreement",
                "TDG10ExactCompleteCRuntimeClosed",
            }
            or closed["retry"] not in authority.SELECTED_RETRIES
            or not isinstance(closed["reason"], str)
            or not closed["reason"]
        ):
            _fail("status", "inconclusive_terminal", closed)
    return value


def _inspect_output(root: Path, *, authority_commit: str) -> dict[str, object]:
    snapshot = _snapshot_output(root)
    if snapshot is None:
        return {"state": "absent"}
    values, _blobs, hashes = snapshot
    terminal = _validate_terminal(
        values["terminal.json"], authority_commit=authority_commit
    )
    if values["manifest.json"] != _manifest(authority_commit):
        _fail("status", "manifest_identity", values["manifest.json"])
    return {
        "state": "terminal",
        "hashes": hashes,
        "classification": str(terminal["classification"]),
        "width_robustness_passed": terminal["width_robustness_passed"] is True,
    }


def status(root: Path, *, authority_commit: str) -> dict[str, object]:
    _execution_authority(root.resolve(), authority_commit)
    before = _snapshot_store(root.resolve())
    if before != (
        authority.SEALED_STORE_LEAF_COUNT,
        authority.SEALED_STORE_SNAPSHOT_SHA256,
    ):
        _fail("status", "sealed_store_identity", before)
    output = _inspect_output(root.resolve(), authority_commit=authority_commit)
    return {
        "artifact_id": authority.ARTIFACT_ID,
        "authority_commit": authority_commit,
        "output": output,
        "store": {"leaf_count": before[0], "sha256": before[1]},
        "safe_to_run": output["state"] == "absent",
        "state_advance_authorized": False,
        "status_read_only": True,
    }


def _width_record_from_assessment(
    restored: ti2_runner.RestoredReplay,
    serialized: Mapping[str, object],
) -> dict[str, object]:
    width = authority.width_spec(int(restored.replay["retry"]))
    passed = serialized["complete_admission_passed"] is True and tuple(
        serialized["failed_channels"]
    ) == ()
    return {
        "retry": width["retry"],
        "width_class": WIDTH_PASS_CLASS if passed else WIDTH_NONPASS_CLASS,
        "predecessor_generation": width["generation"],
        "attempted_width_hex": width["attempted_width_hex"],
        "replay_receipt": _replay_receipt(restored),
        **_shadow_budget(),
        "exact_complete_C": serialized,
        "complete_admission_passed": passed,
    }


def run(root: Path, *, authority_commit: str) -> dict[str, object]:
    repository = root.resolve()
    receipt = _execution_authority(repository, authority_commit)
    if (
        receipt.tableau_selector != authority.TABLEAU_SELECTOR
        or receipt.selected_retries != authority.SELECTED_RETRIES
        or receipt.state_advance_authorized
        or receipt.new_protocol_id_authorized
        or receipt.diagnostic_fine_endpoint_serialized
        or receipt.production_method_earned
        or receipt.width_robustness_passed
        or not receipt.diagnostic_qualification_only
    ):
        _fail("authority", "receipt_semantics", receipt)
    authority.require_output_absent(repository)
    before = _snapshot_store(repository)
    sealed = (
        authority.SEALED_STORE_LEAF_COUNT,
        authority.SEALED_STORE_SNAPSHOT_SHA256,
    )
    if before != sealed:
        _fail("provenance", "sealed_store_identity", (before, sealed))

    width_records: list[dict[str, object]] = []
    typed_stop: dict[str, object] | None = None
    typed_inconclusive: dict[str, object] | None = None
    for replay in authority.REPLAYS:
        width = authority.width_spec(int(replay["retry"]))
        restored = _restore_replay(repository, replay)
        replay_receipt = _replay_receipt(restored)
        try:
            try:
                prepared = _prepare_shadow(restored)
            except tdg6.TDG6RefinementPathStop as exc:
                record = {
                    "retry": width["retry"],
                    "width_class": "premise_stop",
                    "predecessor_generation": width["generation"],
                    "attempted_width_hex": width["attempted_width_hex"],
                    "replay_receipt": replay_receipt,
                    "typed_stop": {
                        "type": type(exc).__name__,
                        "detail": " ".join(str(exc).split())[:640],
                    },
                }
                width_records.append(record)
                if typed_stop is None:
                    typed_stop = {**record["typed_stop"], "retry": width["retry"]}
                continue
            try:
                assessment = _assess_prepared_exact_complete_c(prepared)
            except INCONCLUSIVE_EXCEPTIONS as exc:
                evidence = getattr(exc, "evidence", None)
                reason = getattr(evidence, "reason", type(exc).__name__)
                record = {
                    "retry": width["retry"],
                    "width_class": "inconclusive",
                    "predecessor_generation": width["generation"],
                    "attempted_width_hex": width["attempted_width_hex"],
                    "replay_receipt": replay_receipt,
                    **_shadow_budget(),
                    "typed_inconclusive": {
                        "type": type(exc).__name__,
                        "reason": str(reason),
                        "detail": " ".join(str(exc).split())[:640],
                    },
                }
                width_records.append(record)
                if typed_inconclusive is None:
                    typed_inconclusive = {
                        **record["typed_inconclusive"],
                        "retry": width["retry"],
                    }
                continue
            serialized = serialize_exact_assessment(assessment)
            width_records.append(_width_record_from_assessment(restored, serialized))
        finally:
            after_fingerprint = ti2_runner._member_fingerprint(
                restored.member, str(restored.fingerprint["descriptor_sha256"])
            )
            if after_fingerprint != restored.fingerprint:
                _fail("provenance", "member_mutated", after_fingerprint)

    after = _snapshot_store(repository)
    if after != before:
        _fail("provenance", "campaign_store_mutated", (before, after))
    classification = reduce_qa2_terminal(width_records)
    robustness = classification == BOTH_PASS_CLASS
    terminal = {
        **_terminal_base(
            authority_commit,
            classification,
            before,
            after,
            width_robustness_passed=robustness,
        ),
        "widths": width_records,
    }
    if classification == PREMISE_STOP_CLASS:
        if typed_stop is None:
            _fail("reduction", "premise_stop_missing", classification)
        terminal["typed_stop"] = typed_stop
    if classification == INCONCLUSIVE_CLASS:
        if typed_inconclusive is None:
            _fail("reduction", "inconclusive_missing", classification)
        terminal["typed_inconclusive"] = typed_inconclusive
    _publish_terminal(
        repository,
        _manifest(authority_commit),
        terminal,
        store_before=before,
        store_after=after,
    )
    return terminal


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
    except QA2RunnerError as exc:
        print(
            json.dumps(
                {
                    "artifact_id": authority.ARTIFACT_ID,
                    "classification": INVALID_CLASS,
                    "owner": exc.owner,
                    "code": exc.code,
                    "detail": exc.detail,
                },
                sort_keys=True,
                separators=(",", ":"),
            )
        )
        return 2
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
