#!/usr/bin/env python3
"""Run one retry-3 SSPRK3-on-SBP4 all-channel production-classifier audit."""

from __future__ import annotations

import argparse
import ctypes
import errno
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite
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
from recursive_horizons.fgc.evolution import tdg6_temporal_admission_runtime as tdg6  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ac1_authority as authority  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ar1_authority as ar1  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_store import (  # noqa: E402
    HLT16CampaignStore,
)
from recursive_horizons.fgc.evolution.numerical_engine import COMPARATOR_METHOD  # noqa: E402
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import (  # noqa: E402
    build_static_gr0_shells,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    CertifiedMagnitudeInterval,
    classify_tdg6_channel,
)


RUNNER_ID = "FGC-1-TDG9-AC1-RUN1"
RAW_SCHEMA = "FGC-1-TDG9-AC1-raw-v1"
TERMINAL_CLASSES = frozenset(
    {
        "completed_retry3_all_18_channels_pass",
        "completed_retry3_one_or_more_channels_fail",
        "shadow_proposal_premise_stop",
        "invalid_provenance_or_implementation",
    }
)
PUBLISHABLE_TERMINAL_CLASSES = TERMINAL_CLASSES - {
    "invalid_provenance_or_implementation"
}
_OUTPUT_LEAF_LIMIT = 512 * 1024 * 1024
_HEX = frozenset("0123456789abcdef")


class AC1RunnerError(RuntimeError):
    """Typed fail-closed AC1 runner error."""

    def __init__(self, owner: str, code: str, detail: object) -> None:
        self.owner = str(owner)
        self.code = str(code)
        self.detail = " ".join(str(detail).replace(str(ROOT), "<repo>").split())[:640]
        super().__init__(f"{self.owner}/{self.code}: {self.detail}")


def _fail(owner: str, code: str, detail: object) -> NoReturn:
    raise AC1RunnerError(owner, code, detail)


def _pretty(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def _unique(items: list[tuple[str, object]]) -> dict[str, object]:
    answer: dict[str, object] = {}
    for key, value in items:
        if key in answer:
            _fail("serialization", "duplicate_JSON_key", key)
        answer[key] = value
    return answer


def _exact_binary64(value: object, *, label: str) -> dict[str, str]:
    if type(value) is not float:
        _fail("serialization", "nonexact_binary64", label)
    if not isfinite(value):
        _fail("serialization", "nonfinite_binary64", label)
    fraction = Fraction.from_float(value)
    encoded = {
        "binary64_hex": value.hex(),
        "numerator": str(fraction.numerator),
        "denominator": str(fraction.denominator),
    }
    rebuilt = float.fromhex(encoded["binary64_hex"])
    if rebuilt.hex() != value.hex() or Fraction.from_float(rebuilt) != fraction:
        _fail("serialization", "binary64_roundtrip", label)
    return encoded


def _fraction_from_exact(value: object, *, label: str) -> Fraction:
    if not isinstance(value, Mapping):
        _fail("serialization", "exact_mapping", label)
    if set(value) != {"binary64_hex", "numerator", "denominator"}:
        _fail("serialization", "exact_keys", label)
    hex_value = value["binary64_hex"]
    numerator = value["numerator"]
    denominator = value["denominator"]
    if (
        not isinstance(hex_value, str)
        or not isinstance(numerator, str)
        or not isinstance(denominator, str)
    ):
        _fail("serialization", "exact_types", label)
    try:
        rebuilt = float.fromhex(hex_value)
        fraction = Fraction(int(numerator), int(denominator))
    except (ValueError, ZeroDivisionError) as exc:
        raise AC1RunnerError("serialization", "exact_parse", label) from exc
    if (
        not isfinite(rebuilt)
        or rebuilt.hex() != hex_value
        or str(fraction.numerator) != numerator
        or str(fraction.denominator) != denominator
        or Fraction.from_float(rebuilt) != fraction
    ):
        _fail("serialization", "exact_roundtrip", label)
    return fraction


def _interval_from_exact(value: object, *, label: str) -> CertifiedMagnitudeInterval:
    if not isinstance(value, Mapping) or set(value) != {"lower", "upper"}:
        _fail("serialization", "interval_keys", label)
    return CertifiedMagnitudeInterval(
        _fraction_from_exact(value["lower"], label=f"{label}.lower"),
        _fraction_from_exact(value["upper"], label=f"{label}.upper"),
    )


def serialize_continuous_admission(
    admission: tdg6.TDG6ContinuousAdmissionEvidence,
) -> dict[str, object]:
    if not isinstance(admission, tdg6.TDG6ContinuousAdmissionEvidence):
        _fail("admission", "type", type(admission).__name__)
    channels: list[dict[str, object]] = []
    for item in admission.channel_admissions:
        channels.append(
            {
                "channel": item.channel,
                "classification": item.classification,
                "admission_passed": item.admission_passed,
                "temporal_retry_permitted": item.temporal_retry_permitted,
                "order_threshold_resolved": item.order_threshold_resolved,
                "order_threshold_passed": item.order_threshold_passed,
                "outer": {
                    "lower": _exact_binary64(
                        item.outer_difference.lower_bound, label=f"{item.channel}.outer.lower"
                    ),
                    "upper": _exact_binary64(
                        item.outer_difference.upper_bound, label=f"{item.channel}.outer.upper"
                    ),
                },
                "finest": {
                    "lower": _exact_binary64(
                        item.finest_difference.lower_bound,
                        label=f"{item.channel}.finest.lower",
                    ),
                    "upper": _exact_binary64(
                        item.finest_difference.upper_bound,
                        label=f"{item.channel}.finest.upper",
                    ),
                },
                "finest_pair_debit": _exact_binary64(
                    item.finest_pair_debit, label=f"{item.channel}.debit"
                ),
            }
        )
    serialized = {
        "method": admission.method,
        "owned_row_count": admission.owned_row_count,
        "channel_order": list(authority.CHANNEL_ORDER),
        "channels": channels,
        "complete_admission_passed": admission.complete_admission_passed,
        "failed_channels": list(admission.failed_channels),
        "finest_pair_debit_vector": [
            _exact_binary64(value, label=f"debit[{index}]")
            for index, value in enumerate(admission.finest_pair_debit_vector)
        ],
        "admission_is_all_of": True,
    }
    validate_all_of_identity(serialized)
    return serialized


def validate_all_of_identity(serialized: Mapping[str, object]) -> None:
    channels = serialized.get("channels")
    if not isinstance(channels, list) or len(channels) != len(authority.CHANNEL_ORDER):
        _fail("admission", "channel_count", type(channels).__name__)
    names = tuple(item.get("channel") for item in channels if isinstance(item, Mapping))
    if names != authority.CHANNEL_ORDER:
        _fail("admission", "channel_order", names)
    failed = tuple(
        item["channel"]
        for item in channels
        if isinstance(item, Mapping) and item.get("admission_passed") is not True
    )
    stored_failed = serialized.get("failed_channels")
    if not isinstance(stored_failed, list):
        _fail("admission", "failed_channels_type", type(stored_failed).__name__)
    if tuple(stored_failed) != failed:
        _fail("admission", "failed_channels_identity", stored_failed)
    complete = serialized.get("complete_admission_passed") is True
    if complete is not (not failed):
        _fail("admission", "complete_pass_identity", complete)
    if complete is not all(
        isinstance(item, Mapping) and item.get("admission_passed") is True
        for item in channels
    ):
        _fail("admission", "all_flags_identity", complete)
    for item in channels:
        if not isinstance(item, Mapping):
            _fail("admission", "channel_type", type(item).__name__)
        outer = _interval_from_exact(item.get("outer"), label=f"{item['channel']}.outer")
        finest = _interval_from_exact(
            item.get("finest"), label=f"{item['channel']}.finest"
        )
        decision = classify_tdg6_channel(outer, finest)
        debit = _fraction_from_exact(
            item.get("finest_pair_debit"), label=f"{item['channel']}.debit"
        )
        if (
            decision.classification != item.get("classification")
            or decision.admission_passed is not item.get("admission_passed")
            or decision.temporal_retry_permitted
            is not item.get("temporal_retry_permitted")
            or decision.order_threshold_resolved
            is not item.get("order_threshold_resolved")
            or decision.order_threshold_passed is not item.get("order_threshold_passed")
            or decision.finest_pair_debit != finest.upper
            or debit != finest.upper
        ):
            _fail("admission", "reclassification", item["channel"])
        if item.get("classification") == "order_inconclusive" and (
            item.get("admission_passed") is not False
        ):
            _fail("admission", "order_inconclusive_not_failure", item["channel"])
    vector = tuple(
        _fraction_from_exact(item, label=f"debit[{index}]")
        for index, item in enumerate(serialized.get("finest_pair_debit_vector", ()))
    )
    expected_vector = tuple(
        _fraction_from_exact(item["finest_pair_debit"], label=f"{item['channel']}.debit")
        for item in channels
        if isinstance(item, Mapping)
    )
    if vector != expected_vector:
        _fail("admission", "debit_vector", len(vector))
    if serialized.get("admission_is_all_of") is not True:
        _fail("admission", "all_of_flag", serialized.get("admission_is_all_of"))


def reduce_retry3_terminal(
    *, complete_admission_passed: bool, failed_channels: Sequence[str]
) -> str:
    if complete_admission_passed is True and tuple(failed_channels) == ():
        return "completed_retry3_all_18_channels_pass"
    return "completed_retry3_one_or_more_channels_fail"


def _execution_authority(root: Path, commit: str) -> authority.AC1Authority:
    try:
        return authority.authorize(
            root,
            authority.read_leaf(root, authority.CONFIG_PATH),
            authority.read_leaf(root, authority.RESULT_PATH),
            commit,
        )
    except Exception as exc:
        raise AC1RunnerError("authority", "rejected", exc) from exc


def _snapshot_store(root: Path) -> tuple[int, str]:
    try:
        return ti2_runner._snapshot_store(root)
    except Exception as exc:
        raise AC1RunnerError("provenance", "store_snapshot", exc) from exc


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


def _publish(
    root: Path, manifest: Mapping[str, object], terminal: Mapping[str, object]
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
        for name, value in (("manifest.json", manifest), ("terminal.json", terminal)):
            leaf = os.open(
                name,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
                0o644,
                dir_fd=stage_fd,
            )
            try:
                pending = memoryview(_pretty(value))
                while pending:
                    written = os.write(leaf, pending)
                    if written <= 0:
                        _fail("output", "short_write", name)
                    pending = pending[written:]
                os.fsync(leaf)
            finally:
                os.close(leaf)
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
            0x00000004,
        ):
            error = ctypes.get_errno()
            if error in {errno.EEXIST, errno.ENOTEMPTY}:
                _fail("output", "namespace_arrived", authority.OUTPUT_NAMESPACE)
            _fail("output", "exclusive_adoption_failed", error)
        os.close(stage_fd)
        stage_fd = None
        os.fsync(parent_fd)
    except Exception:
        if stage_fd is not None:
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
        raise AC1RunnerError("status", "unsafe_output_leaf", name) from exc
    finally:
        if descriptor != -1:
            os.close(descriptor)


def _snapshot_output(root: Path) -> tuple[dict[str, object], dict[str, str]] | None:
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
    edge_names: list[str] = []
    leaf_identities: dict[str, tuple[int, int, int, int, int, int]] = {}
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
            edge_names.append(part)
            current = child
        output_fd = descriptors[-1]
        names = tuple(sorted(os.listdir(output_fd)))
        expected_names = ("manifest.json", "terminal.json")
        if names != expected_names:
            _fail("status", "partial_or_foreign_namespace", names)
        values: dict[str, object] = {}
        hashes: dict[str, str] = {}
        for name in names:
            before = os.stat(name, dir_fd=output_fd, follow_symlinks=False)
            leaf_identities[name] = _identity(before)
            raw = _read_output_leaf(output_fd, name, before)
            value = _parse_output_json(raw, name=name)
            values[name] = value
            hashes[name] = sha256(raw).hexdigest()
        if tuple(sorted(os.listdir(output_fd))) != expected_names:
            _fail("status", "output_changed_during_snapshot", names)
        for name, expected in leaf_identities.items():
            after = os.stat(name, dir_fd=output_fd, follow_symlinks=False)
            if _identity(after) != expected:
                _fail("status", "output_leaf_substituted", name)
        for index, (descriptor, expected) in enumerate(
            zip(descriptors, identities, strict=True)
        ):
            if _identity(os.fstat(descriptor)) != expected:
                _fail("status", "output_directory_changed", index)
        for index, name in enumerate(edge_names, start=1):
            after = os.stat(name, dir_fd=descriptors[index - 1], follow_symlinks=False)
            if _identity(after) != identities[index]:
                _fail("status", "output_directory_substituted", name)
        if _identity(root.lstat()) != identities[0]:
            _fail("status", "repository_root_substituted", root)
        return values, hashes
    except OSError as exc:
        raise AC1RunnerError(
            "status", "output_snapshot_unreadable", authority.OUTPUT_NAMESPACE
        ) from exc
    finally:
        for descriptor in reversed(descriptors):
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
        raise AC1RunnerError("status", "noncanonical_output", name) from exc
    if not isinstance(value, dict) or raw != _pretty(value):
        _fail("status", "noncanonical_output", name)
    return value


def _require_keys(
    value: object, expected: set[str], *, label: str
) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != expected:
        _fail("status", "terminal_schema", (label, sorted(expected)))
    return value


def expected_replay_receipt() -> dict[str, object]:
    return authority.expected_raw_replay_receipt()


def _validate_replay_receipt(value: object) -> dict[str, object]:
    expected = expected_replay_receipt()
    if not isinstance(value, Mapping) or set(value) != set(expected):
        _fail("status", "replay_receipt_schema", (sorted(expected), value))
    receipt = {key: value[key] for key in expected}
    if receipt != expected:
        _fail("status", "replay_receipt_identity", receipt)
    return receipt


def _terminal_base(
    authority_commit: str,
    classification: str,
    before: tuple[int, str],
    after: tuple[int, str],
) -> dict[str, object]:
    if classification not in TERMINAL_CLASSES:
        _fail("reduction", "terminal_class", classification)
    return {
        "schema": RAW_SCHEMA,
        "artifact_id": authority.ARTIFACT_ID,
        "runner_id": RUNNER_ID,
        "classification": classification,
        "authority_commit": authority_commit,
        "experiment_label": authority.SEMANTICS["experiment_label"],
        "tableau_selector": "SSPRK3",
        "tableau_runtime_selector": COMPARATOR_METHOD,
        "actual_spatial_operator": "inherited_RK4_2049_SBP4",
        "retry": authority.RETRY,
        "attempted_width_hex": authority.WIDTH_HEX,
        "channel_order": list(authority.CHANNEL_ORDER),
        "production_SSPRK3_comparator": False,
        "independent_method_agreement": False,
        "production_method_earned": False,
        "store_snapshot_before": {"leaf_count": before[0], "sha256": before[1]},
        "store_snapshot_after": {"leaf_count": after[0], "sha256": after[1]},
        "store_unchanged": before == after,
        "PDE_state_committed": False,
        "temporal_retry_admission_called": False,
        "fine_path_committed": False,
        "successor_remedy_selected": False,
        "common_event_completed": False,
        "GR0_calibration_completed": False,
        "candidate_branch_opened": False,
        "mechanism_result_earned": False,
        "physical_result_earned": False,
        "all_pass_licenses_only_this_frozen_state_and_width": True,
        "fail_licenses_only_diagnosis_of_failed_channel_set": True,
    }


def _manifest(authority_commit: str) -> dict[str, object]:
    return {
        "schema": RAW_SCHEMA,
        "artifact_id": authority.ARTIFACT_ID,
        "runner_id": RUNNER_ID,
        "authority_commit": authority_commit,
        "TI2_PREF1_result_sha256": authority.PREF1_RESULT_SHA256,
        "experiment_label": authority.SEMANTICS["experiment_label"],
        "tableau_selector": "SSPRK3",
        "actual_spatial_operator": "inherited_RK4_2049_SBP4",
        "retry": authority.RETRY,
        "production_SSPRK3_comparator": False,
        "output_leaves": ["manifest.json", "terminal.json"],
    }


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
    expected_base = _terminal_base(
        authority_commit, str(classification), sealed, sealed
    )
    if any(value.get(name) != expected for name, expected in expected_base.items()):
        _fail("status", "terminal_identity", classification)
    base_keys = set(expected_base)
    if classification == "shadow_proposal_premise_stop":
        _require_keys(
            value,
            base_keys | {"typed_stop", "replay_receipt"},
            label="premise stop terminal",
        )
        stop = _require_keys(
            value["typed_stop"], {"type", "detail"}, label="premise stop"
        )
        if (
            stop["type"] != "TDG6RefinementPathStop"
            or not isinstance(stop["detail"], str)
            or not stop["detail"]
            or len(stop["detail"]) > 640
        ):
            _fail("status", "premise_stop_terminal", stop)
        _validate_replay_receipt(value["replay_receipt"])
        return value
    _require_keys(
        value,
        base_keys
        | {
            "replay_receipt",
            "shadow_path_count",
            "shadow_proposal_count",
            "SSPRK3_stage_and_endpoint_record_count",
            "continuous_admission",
        },
        label="completed terminal",
    )
    _validate_replay_receipt(value["replay_receipt"])
    admission = value.get("continuous_admission")
    if not isinstance(admission, Mapping):
        _fail("status", "admission_schema", type(admission).__name__)
    validate_all_of_identity(admission)
    failed = tuple(admission["failed_channels"])
    complete = admission["complete_admission_passed"] is True
    expected = reduce_retry3_terminal(
        complete_admission_passed=complete, failed_channels=failed
    )
    if (
        classification != expected
        or value.get("shadow_path_count") != 7
        or value.get("shadow_proposal_count") != 7
        or value.get("SSPRK3_stage_and_endpoint_record_count") != 28
        or value.get("production_method_earned") is not False
    ):
        _fail("status", "completed_terminal_counts", classification)
    return value


def _inspect_output(root: Path, *, authority_commit: str) -> dict[str, object]:
    snapshot = _snapshot_output(root)
    if snapshot is None:
        return {"state": "absent"}
    values, hashes = snapshot
    manifest = values["manifest.json"]
    if manifest != _manifest(authority_commit):
        _fail("status", "manifest_identity", manifest)
    terminal = _validate_terminal(
        values["terminal.json"], authority_commit=authority_commit
    )
    return {
        "state": "terminal",
        "hashes": hashes,
        "classification": terminal["classification"],
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
    }


def _restore_replay(root: Path) -> ti2_runner.RestoredReplay:
    store = HLT16CampaignStore(root / ar1.PREF2_STORE_PATH)
    shells = build_static_gr0_shells(root)
    try:
        return ti2_runner._restore_replay(root, store, shells, authority.REPLAY)
    except ti2_runner.TI1RunnerError as exc:
        raise AC1RunnerError(exc.owner, exc.code, exc.detail) from exc


def _prepare_shadow(
    restored: ti2_runner.RestoredReplay,
) -> tdg6.TDG6PreparedGR0Compositor:
    try:
        return ti2_runner._prepare_shadow(restored)
    except tdg6.TDG6RefinementPathStop:
        raise
    except ti2_runner.TI1RunnerError as exc:
        raise AC1RunnerError(exc.owner, exc.code, exc.detail) from exc


def _replay_receipt(restored: ti2_runner.RestoredReplay) -> dict[str, object]:
    fingerprint = restored.fingerprint
    if not isinstance(fingerprint, Mapping):
        _fail("replay", "fingerprint_type", type(fingerprint).__name__)
    return _validate_replay_receipt(
        {
            "member_key": fingerprint.get("member_key"),
            "retry": authority.RETRY,
            "attempted_width_hex": authority.WIDTH_HEX,
            "accepted_time_hex": fingerprint.get("accepted_time_hex"),
            "state_sha256": fingerprint.get("state_sha256"),
            "descriptor_sha256": fingerprint.get("descriptor_sha256"),
            "transaction_sha256": fingerprint.get("transaction_sha256"),
            "historical_journal_sha256": restored.historical_journal_sha256,
        }
    )


def run(root: Path, *, authority_commit: str) -> dict[str, object]:
    repository = root.resolve()
    receipt = _execution_authority(repository, authority_commit)
    if (
        receipt.tableau_selector != "SSPRK3"
        or receipt.retry != authority.RETRY
        or receipt.state_advance_authorized
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
    replay_receipt: dict[str, object] = {}
    try:
        restored = _restore_replay(repository)
        replay_receipt = _replay_receipt(restored)
        prepared = _prepare_shadow(restored)
    except tdg6.TDG6RefinementPathStop as exc:
        after = _snapshot_store(repository)
        if after != before:
            _fail("provenance", "campaign_store_mutated", (before, after))
        terminal = {
            **_terminal_base(
                authority_commit, "shadow_proposal_premise_stop", before, after
            ),
            "typed_stop": {
                "type": type(exc).__name__,
                "detail": " ".join(str(exc).split())[:640],
            },
            "replay_receipt": replay_receipt,
        }
        _publish_terminal(
            repository,
            _manifest(authority_commit),
            terminal,
            store_before=before,
            store_after=after,
        )
        return terminal

    serialized = serialize_continuous_admission(prepared.continuous_admission)
    classification = reduce_retry3_terminal(
        complete_admission_passed=bool(serialized["complete_admission_passed"]),
        failed_channels=tuple(serialized["failed_channels"]),
    )
    after = _snapshot_store(repository)
    if after != before:
        _fail("provenance", "campaign_store_mutated", (before, after))
    terminal = {
        **_terminal_base(authority_commit, classification, before, after),
        "replay_receipt": replay_receipt,
        "shadow_path_count": 7,
        "shadow_proposal_count": 7,
        "SSPRK3_stage_and_endpoint_record_count": 28,
        "continuous_admission": serialized,
    }
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
    except AC1RunnerError as exc:
        print(
            json.dumps(
                {
                    "artifact_id": authority.ARTIFACT_ID,
                    "classification": "invalid_provenance_or_implementation",
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
