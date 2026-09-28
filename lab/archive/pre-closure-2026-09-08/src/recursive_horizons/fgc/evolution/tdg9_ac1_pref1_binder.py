"""Independent post-result binder for the completed TDG9 AC1 terminal.

The one-time live path authenticates the immutable AC1 authority and the two
canonical raw leaves, independently decodes every exact binary64/rational
interval, and reclassifies all eighteen channels through the design-only
TDG6 classifier.  It deliberately does not import the AC1 runner, does not
import TI2 shadow restore or prepare helpers, and does not reconstruct the
seven SSPRK3 proposals.  Ordinary verification consumes only the tracked
compact certificate and is raw/store blind.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import stat
import subprocess
import tomllib
from typing import Any, Mapping, NoReturn, Sequence

from .tdg6_temporal_admission_design import (
    CertifiedMagnitudeInterval,
    TDG6_COMPLETE_STATE_CHANNELS,
    classify_tdg6_channel,
)


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG9-AC1-PREF1"
CLASSIFICATION = (
    "independently_bound_retry3_all_channel_production_classifier_terminal"
)
CONFIG_PATH = "configs/fgc/fgc-1-tdg9-ac1-pref1.toml"
RESULT_PATH = "results/fgc-1-tdg9-ac1-pref1.json"
OWNER_DOCUMENT = "docs/fgc-tdg9-ac1-pref1.md"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"

AUTHORITY_COMMIT = "ec4195d5f09788bad3e08936b9506716a9eafa3d"
AUTHORITY_PARENT = "660f369365ab8417f4c68eded46acb9bb2b05a5b"
RAW_NAMESPACE = "runs/fgc-2-sf1/tdg9-ac1/ssprk3-sbp4-retry3-all-channels"
RAW_SCHEMA = "FGC-1-TDG9-AC1-raw-v1"
RAW_CLASSIFICATION = "completed_retry3_one_or_more_channels_fail"
RUNNER_ID = "FGC-1-TDG9-AC1-RUN1"
RAW_MANIFEST_SHA256 = (
    "5cd74aa33555c6253353311c246515680b50c8c908a4ff295949a98a7e8f0912"
)
RAW_TERMINAL_SHA256 = (
    "515ddf31d2fe919e4228668bea8701f282931d45e462542c26675382c95b834c"
)
COMPACT_RESULT_SHA256 = (
    "173183871b2d750f98e9758cdc5f702be40b0bbb73d10e707672e0e68752935d"
)

STORE_PATH = "runs/fgc-2-sf1/tdg8-rcv3/calibration"
SEALED_STORE_LEAF_COUNT = 115
SEALED_STORE_SHA256 = (
    "5917d70de3080dcc6b7abeb7fdb7cc3370c6140cbeb37bde0c142d6a2d36a445"
)
AC1_ARTIFACT_ID = "FGC-1-TDG9-AC1-FRZ1"
TI2_PREF1_RESULT_SHA256 = (
    "4a7a4bb15c35766b9a73dd734f657aeae2c99692c37b00d7a9697dc369ccf546"
)
EXPERIMENT_LABEL = (
    "retry3_SSPRK3_tableau_on_inherited_SBP4_all_18_channel_"
    "production_classifier_audit"
)
TABLEAU_RUNTIME_SELECTOR = "second_order_diagonal_norm_SBP_plus_SSPRK3"
MEMBER_KEY = "RK4-2049"
PHYSICAL_STATE_SHA256 = (
    "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a"
)
ACCEPTED_TIME_HEX = "0x1.78554de5a30e0p+0"
POINT_COUNT = 2049
OWNED_ROW_COUNT = 2044
RETRY = 3
WIDTH_HEX = "0x1.aaa9612df8000p-11"
MEMBER_DESCRIPTOR_SHA256 = (
    "77847126340e78c8ac300fac2795bceda724c15e0894e34444b0da42a84166b4"
)
TRANSACTION_SHA256 = (
    "7a90163d2eb4252fa7a1bbdf55d12f92a9127ed7cba5a37c28456e72d732b22b"
)
HISTORICAL_JOURNAL_SHA256 = (
    "0d21f650ccc46db778fd81fbebf5acab6394e1c2f7972b7940e24c76af9e1f59"
)
WORK_BUDGET = {
    "FAILED_OCCURRENCES_filter_applied": False,
    "LOC1_exact_localization_authorized": False,
    "LOC2_exact_localization_authorized": False,
    "SSPRK3_records_per_proposal": 4,
    "channel_count": 18,
    "fourth_width_authorized": False,
    "maximum_stage_and_endpoint_RHS_records": 28,
    "resource_escalation_authorized": False,
    "retry": 3,
    "retry_4_authorized": False,
    "retry_5_authorized": False,
    "retry_count": 1,
    "shadow_paths": 7,
    "shadow_proposals": 7,
}
CHANNEL_ORDER = TDG6_COMPLETE_STATE_CHANNELS
EXPECTED_ADMITTED_COUNT = 17
EXPECTED_FAILED_CHANNELS = ("u:R",)
EXPECTED_REPLAY_RECEIPT = {
    "accepted_time_hex": ACCEPTED_TIME_HEX,
    "attempted_width_hex": WIDTH_HEX,
    "descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
    "historical_journal_sha256": HISTORICAL_JOURNAL_SHA256,
    "member_key": MEMBER_KEY,
    "retry": RETRY,
    "state_sha256": PHYSICAL_STATE_SHA256,
    "transaction_sha256": TRANSACTION_SHA256,
}
EXPECTED_CHANNEL_CLASSES: tuple[
    tuple[str, str, bool, bool, bool, bool | None], ...
] = (
    ("u:alpha", "enclosure_dominated_debit_only", True, False, False, None),
    ("u:v", "resolved_order_pass", True, False, True, True),
    ("u:lambda", "enclosure_dominated_debit_only", True, False, False, None),
    ("u:R", "order_inconclusive", False, True, False, None),
    ("u:phi", "resolved_order_pass", True, False, True, True),
    ("u:chi", "resolved_order_pass", True, False, True, True),
    ("p:alpha", "resolved_order_pass", True, False, True, True),
    ("p:v", "resolved_order_pass", True, False, True, True),
    ("p:lambda", "resolved_order_pass", True, False, True, True),
    ("p:R", "resolved_order_pass", True, False, True, True),
    ("p:phi", "resolved_order_pass", True, False, True, True),
    ("p:chi", "resolved_order_pass", True, False, True, True),
    ("q:alpha", "resolved_order_pass", True, False, True, True),
    ("q:v", "resolved_order_pass", True, False, True, True),
    ("q:lambda", "resolved_order_pass", True, False, True, True),
    ("q:R", "resolved_order_pass", True, False, True, True),
    ("q:phi", "resolved_order_pass", True, False, True, True),
    ("q:chi", "resolved_order_pass", True, False, True, True),
)
EXPECTED_CLASSIFICATION_COUNTS = {
    "enclosure_dominated_debit_only": 2,
    "order_inconclusive": 1,
    "resolved_order_pass": 15,
}
_EXACT_KEYS = frozenset({"binary64_hex", "numerator", "denominator"})
_INTERVAL_KEYS = frozenset({"lower", "upper"})
_CHANNEL_KEYS = frozenset(
    {
        "admission_passed",
        "channel",
        "classification",
        "finest",
        "finest_pair_debit",
        "order_threshold_passed",
        "order_threshold_resolved",
        "outer",
        "temporal_retry_permitted",
    }
)
_ADMISSION_KEYS = frozenset(
    {
        "admission_is_all_of",
        "channel_order",
        "channels",
        "complete_admission_passed",
        "failed_channels",
        "finest_pair_debit_vector",
        "method",
        "owned_row_count",
    }
)
_COMPLETED_EXTRA_KEYS = frozenset(
    {
        "SSPRK3_stage_and_endpoint_record_count",
        "continuous_admission",
        "replay_receipt",
        "shadow_path_count",
        "shadow_proposal_count",
    }
)
_MAX_RAW_LEAF_BYTES = 8 * 1024 * 1024
_MAX_STORE_LEAF_BYTES = 128 * 1024 * 1024

_AUTHORITY_BLOBS = (
    (
        "configs/fgc/fgc-1-tdg9-ac1-frz1.toml",
        "aab7b3931a43131f57af035063033ae285f0f196a34c32d81ba7f4c99fe5eaf0",
    ),
    (
        "results/fgc-1-tdg9-ac1-frz1.json",
        "e0019fb6e22d5213a6b170d866af3965c45dbb1b83a1c9d612f1218df520b8af",
    ),
    (
        "docs/fgc-tdg9-ac1-frz1.md",
        "9dc2a2f7bb2f9df89a93acd25fb7d22b41eba490449e06914b6b33f0b0fbebfb",
    ),
    (
        "scripts/run_fgc_tdg9_ac1.py",
        "bc2f28e99a747caabac3f419b05574bfe06fcddc5e7b5f1d21879e3639e77605",
    ),
    (
        "src/recursive_horizons/fgc/evolution/tdg9_ac1_authority.py",
        "3ccd293ff67413466af9390065357e68226ddd46ee23df7bceda10aec6dc4c19",
    ),
    (
        "scripts/reproduce_fgc_tdg9_ac1_frz1.py",
        "d7911135924eaddb25a09d29534b1f121f229b2df58df16656bca96778334f99",
    ),
)


class TDG9AC1PREF1Error(ValueError):
    """The compact contract, raw terminal, store, or independent replay differs."""

    def __init__(self, stop_id: str, detail: object) -> None:
        self.stop_id = str(stop_id)
        self.detail = " ".join(str(detail).split())[:640]
        super().__init__(f"{self.stop_id}: {self.detail}")


def _stop(stop_id: str, detail: object) -> NoReturn:
    raise TDG9AC1PREF1Error(stop_id, detail)


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
        raise TDG9AC1PREF1Error("PREF1_COMPACT_DRIFT", exc) from exc


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
        raise TDG9AC1PREF1Error("PREF1_JSON_DRIFT", label) from exc
    if not isinstance(value, dict) or canonical_result(value) != raw:
        _stop("PREF1_JSON_DRIFT", f"{label} is not canonical pretty JSON")
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
        raise TDG9AC1PREF1Error("PREF1_GIT_DRIFT", arguments) from exc


def _identity(metadata: os.stat_result) -> tuple[int, int, int, int, int]:
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
    )


def _open_directory(root: Path, relative: str) -> int:
    candidate = Path(relative)
    if candidate.is_absolute() or any(
        part in {"", ".", ".."} for part in candidate.parts
    ):
        _stop("PREF1_PATH_UNSAFE", relative)
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
            if not stat.S_ISDIR(os.fstat(child).st_mode):
                os.close(child)
                _stop("PREF1_PATH_UNSAFE", relative)
            os.close(descriptor)
            descriptor = child
        return descriptor
    except Exception:
        os.close(descriptor)
        raise


def _read_leaf_at(
    directory_fd: int, name: str, *, maximum: int, label: str
) -> bytes:
    if not name or "/" in name or name in {".", ".."}:
        _stop("PREF1_PATH_UNSAFE", label)
    before = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    if (
        stat.S_ISLNK(before.st_mode)
        or not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or before.st_size < 0
        or before.st_size > maximum
    ):
        _stop("PREF1_LEAF_UNSAFE", label)
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=directory_fd,
    )
    try:
        active = os.fstat(descriptor)
        if _identity(active) != _identity(before):
            _stop("PREF1_LEAF_RACED", label)
        chunks: list[bytes] = []
        remaining = before.st_size
        while remaining:
            block = os.read(descriptor, min(1 << 20, remaining))
            if not block:
                _stop("PREF1_LEAF_SHORT_READ", label)
            chunks.append(block)
            remaining -= len(block)
        if os.read(descriptor, 1):
            _stop("PREF1_LEAF_GREW", label)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    final = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    if _identity(before) != _identity(after) or _identity(before) != _identity(final):
        _stop("PREF1_LEAF_CHANGED", label)
    return b"".join(chunks)


def _raw_pair(root: Path) -> tuple[bytes, bytes]:
    descriptor = _open_directory(root, RAW_NAMESPACE)
    try:
        before = os.fstat(descriptor)
        with os.scandir(descriptor) as entries:
            names = tuple(sorted(item.name for item in entries))
        if names != ("manifest.json", "terminal.json"):
            _stop("PREF1_RAW_TREE_DRIFT", names)
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
        _stop("PREF1_RAW_TREE_CHANGED", RAW_NAMESPACE)
    if (
        sha256(manifest_raw).hexdigest() != RAW_MANIFEST_SHA256
        or sha256(terminal_raw).hexdigest() != RAW_TERMINAL_SHA256
    ):
        _stop("PREF1_RAW_HASH_DRIFT", RAW_NAMESPACE)
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

    base_fd = _open_directory(root, STORE_PATH)
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
                        _stop("PREF1_STORE_SYMLINK", relative)
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
                            _stop("PREF1_STORE_DIRECTORY_RACED", relative)
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
                        _stop("PREF1_STORE_FOREIGN_TYPE", relative)
                with os.scandir(directory_fd) as entries:
                    final_names = tuple(sorted(item.name for item in entries))
                after = os.fstat(directory_fd)
                if initial_names != final_names or _identity(before) != _identity(
                    after
                ):
                    for child_fd, _ in children:
                        os.close(child_fd)
                    _stop("PREF1_STORE_DIRECTORY_CHANGED", prefix or ".")
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


def _mapping(value: object, keys: frozenset[str], label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or set(value) != keys:
        _stop("PREF1_SCHEMA_DRIFT", (label, sorted(keys)))
    return value


def _fraction_from_exact(value: object, *, label: str) -> Fraction:
    encoded = _mapping(value, _EXACT_KEYS, label)
    hex_value = encoded["binary64_hex"]
    numerator = encoded["numerator"]
    denominator = encoded["denominator"]
    if (
        not isinstance(hex_value, str)
        or not isinstance(numerator, str)
        or not isinstance(denominator, str)
    ):
        _stop("PREF1_EXACT_TYPE_DRIFT", label)
    try:
        rebuilt = float.fromhex(hex_value)
        fraction = Fraction(int(numerator), int(denominator))
    except (ValueError, ZeroDivisionError) as exc:
        raise TDG9AC1PREF1Error("PREF1_EXACT_PARSE_DRIFT", label) from exc
    if (
        type(rebuilt) is not float
        or not math.isfinite(rebuilt)
        or rebuilt.hex() != hex_value
        or str(fraction.numerator) != numerator
        or str(fraction.denominator) != denominator
        or Fraction.from_float(rebuilt) != fraction
    ):
        _stop("PREF1_EXACT_ROUNDTRIP_DRIFT", label)
    return fraction


def _interval_from_exact(value: object, *, label: str) -> CertifiedMagnitudeInterval:
    encoded = _mapping(value, _INTERVAL_KEYS, f"{label} interval")
    try:
        return CertifiedMagnitudeInterval(
            _fraction_from_exact(encoded["lower"], label=f"{label}.lower"),
            _fraction_from_exact(encoded["upper"], label=f"{label}.upper"),
        )
    except (TypeError, ValueError) as exc:
        raise TDG9AC1PREF1Error("PREF1_INTERVAL_DRIFT", label) from exc


def reduce_retry3_terminal(
    *, complete_admission_passed: bool, failed_channels: Sequence[str]
) -> str:
    if complete_admission_passed is True and tuple(failed_channels) == ():
        return "completed_retry3_all_18_channels_pass"
    return "completed_retry3_one_or_more_channels_fail"


def _channel_summary(item: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "admission_passed": item["admission_passed"],
        "channel": item["channel"],
        "classification": item["classification"],
        "debit_equals_finest_upper": True,
        "order_threshold_passed": item["order_threshold_passed"],
        "order_threshold_resolved": item["order_threshold_resolved"],
        "temporal_retry_permitted": item["temporal_retry_permitted"],
    }


def _reclassify_admission(admission: Mapping[str, Any]) -> dict[str, Any]:
    payload = _mapping(admission, _ADMISSION_KEYS, "continuous_admission")
    if (
        payload.get("admission_is_all_of") is not True
        or payload.get("method") != TABLEAU_RUNTIME_SELECTOR
        or payload.get("owned_row_count") != OWNED_ROW_COUNT
        or payload.get("channel_order") != list(CHANNEL_ORDER)
    ):
        _stop("PREF1_ADMISSION_CONTRACT_DRIFT", "all-of identity header")
    channels = payload.get("channels")
    vector = payload.get("finest_pair_debit_vector")
    if not isinstance(channels, list) or len(channels) != len(CHANNEL_ORDER):
        _stop("PREF1_CHANNEL_COUNT_DRIFT", type(channels).__name__)
    if not isinstance(vector, list) or len(vector) != len(CHANNEL_ORDER):
        _stop("PREF1_DEBIT_VECTOR_DRIFT", type(vector).__name__)
    summaries: list[dict[str, Any]] = []
    failed: list[str] = []
    counts: dict[str, int] = {}
    observed_classes: list[tuple[str, str, bool, bool, bool, bool | None]] = []
    for index, (raw_channel, name) in enumerate(
        zip(channels, CHANNEL_ORDER, strict=True)
    ):
        if not isinstance(raw_channel, Mapping):
            _stop("PREF1_CHANNEL_TYPE_DRIFT", name)
        item = _mapping(raw_channel, _CHANNEL_KEYS, name)
        if item.get("channel") != name:
            _stop("PREF1_CHANNEL_ORDER_DRIFT", (name, item.get("channel")))
        outer = _interval_from_exact(item.get("outer"), label=f"{name}.outer")
        finest = _interval_from_exact(item.get("finest"), label=f"{name}.finest")
        try:
            decision = classify_tdg6_channel(outer, finest)
        except (TypeError, ValueError) as exc:
            raise TDG9AC1PREF1Error("PREF1_RECLASSIFY_DRIFT", name) from exc
        debit = _fraction_from_exact(
            item.get("finest_pair_debit"), label=f"{name}.debit"
        )
        vector_debit = _fraction_from_exact(
            vector[index], label=f"debit[{index}]"
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
            or vector_debit != finest.upper
        ):
            _stop("PREF1_RECLASSIFY_DRIFT", name)
        if item.get("classification") == "order_inconclusive" and (
            item.get("admission_passed") is not False
        ):
            _stop("PREF1_ORDER_INCONCLUSIVE_NOT_FAILURE", name)
        if item.get("admission_passed") is not True:
            failed.append(name)
        classification = str(item["classification"])
        counts[classification] = counts.get(classification, 0) + 1
        observed_classes.append(
            (
                name,
                classification,
                bool(item["admission_passed"]),
                bool(item["temporal_retry_permitted"]),
                bool(item["order_threshold_resolved"]),
                item["order_threshold_passed"]
                if item["order_threshold_passed"] is None
                else bool(item["order_threshold_passed"]),
            )
        )
        summaries.append(_channel_summary(item))
    stored_failed = payload.get("failed_channels")
    complete = payload.get("complete_admission_passed") is True
    if not isinstance(stored_failed, list) or tuple(stored_failed) != tuple(failed):
        _stop("PREF1_FAILED_CHANNELS_DRIFT", stored_failed)
    if complete is not (not failed):
        _stop("PREF1_COMPLETE_PASS_DRIFT", complete)
    if complete is not all(item["admission_passed"] is True for item in summaries):
        _stop("PREF1_ALL_FLAGS_DRIFT", complete)
    if tuple(observed_classes) != EXPECTED_CHANNEL_CLASSES:
        _stop("PREF1_CHANNEL_CLASS_DRIFT", observed_classes)
    if dict(sorted(counts.items())) != EXPECTED_CLASSIFICATION_COUNTS:
        _stop("PREF1_CLASSIFICATION_COUNT_DRIFT", counts)
    if tuple(failed) != EXPECTED_FAILED_CHANNELS:
        _stop("PREF1_FAILED_SET_DRIFT", failed)
    if len(summaries) - len(failed) != EXPECTED_ADMITTED_COUNT:
        _stop("PREF1_ADMITTED_COUNT_DRIFT", len(summaries) - len(failed))
    return {
        "all_of_identity_verified": True,
        "admitted_count": EXPECTED_ADMITTED_COUNT,
        "channel_count": len(CHANNEL_ORDER),
        "channels": summaries,
        "classification_counts": dict(EXPECTED_CLASSIFICATION_COUNTS),
        "complete_admission_passed": False,
        "debit_equals_finest_upper_count": len(CHANNEL_ORDER),
        "failed_channels": list(EXPECTED_FAILED_CHANNELS),
        "failed_channels_match_admission_flags": True,
        "failed_count": len(EXPECTED_FAILED_CHANNELS),
        "raw_terminal_class_matches_independent_reduction": True,
    }


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
            "path": STORE_PATH,
            "leaf_count": SEALED_STORE_LEAF_COUNT,
            "snapshot_sha256": SEALED_STORE_SHA256,
            "snapshot_algorithm": "TDG9-LOC1-STORE-SNAPSHOT-v1",
        },
        "selection": {
            "experiment_label": EXPERIMENT_LABEL,
            "tableau_selector": "SSPRK3",
            "tableau_runtime_selector": TABLEAU_RUNTIME_SELECTOR,
            "actual_spatial_operator": "inherited_RK4_2049_SBP4",
            "production_SSPRK3_comparator": False,
            "independent_method_agreement": False,
            "member_key": MEMBER_KEY,
            "physical_state_sha256": PHYSICAL_STATE_SHA256,
            "accepted_time_hex": ACCEPTED_TIME_HEX,
            "point_count": POINT_COUNT,
            "owned_row_count": OWNED_ROW_COUNT,
            "retry": RETRY,
            "attempted_width_hex": WIDTH_HEX,
            "channel_order": list(CHANNEL_ORDER),
        },
        "work_budget": dict(WORK_BUDGET),
        "scope": {
            "one_time_live_raw_and_store_authentication": True,
            "ordinary_verifier_raw_blind": True,
            "ordinary_verifier_store_blind": True,
            "AC1_runner_decision_code_imported": False,
            "TI2_shadow_restore_or_prepare_imported": False,
            "seven_proposals_reconstructed": False,
            "same_state_tableau_counterfactual_only": True,
            "production_SSPRK3_comparator": False,
            "independent_method_agreement": False,
            "campaign_store_mutation_authorized": False,
            "raw_namespace_mutation_authorized": False,
            "temporal_retry_admission_authorized": False,
            "PDE_state_commit_authorized": False,
            "successor_remedy_selected": False,
            "common_event_authorized": False,
            "GR0_calibration_authorized": False,
            "candidate_branches_authorized": False,
            "mechanism_result_authorized": False,
            "physical_result_authorized": False,
        },
        "claims": {
            "AC1_raw_result_independently_bound": True,
            "seventeen_of_eighteen_channels_admitted": True,
            "only_u_R_nonadmitted": True,
            "u_R_class_order_inconclusive": True,
            "production_method_earned": False,
            "production_SSPRK3_comparator": False,
            "independent_method_agreement": False,
            "state_advance_authorized": False,
            "successor_remedy_selected": False,
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
        raise TDG9AC1PREF1Error("PREF1_CONFIG_DRIFT", exc) from exc
    if value != _expected_config():
        _stop("PREF1_CONFIG_DRIFT", "typed PREF1 config differs")
    return value


def _authenticate_authority(root: Path, manifest: Mapping[str, Any]) -> None:
    if manifest.get("authority_commit") != AUTHORITY_COMMIT:
        _stop("PREF1_AUTHORITY_DRIFT", manifest.get("authority_commit"))
    if _git(root, "cat-file", "-t", AUTHORITY_COMMIT).strip() != b"commit":
        _stop("PREF1_AUTHORITY_DRIFT", "authority object is not a commit")
    parents = _git(root, "show", "-s", "--format=%P", AUTHORITY_COMMIT).split()
    if parents != [AUTHORITY_PARENT.encode("ascii")]:
        _stop("PREF1_AUTHORITY_LINEAGE_DRIFT", parents)
    for relative, expected in _AUTHORITY_BLOBS:
        observed = sha256(
            _git(root, "show", f"{AUTHORITY_COMMIT}:{relative}")
        ).hexdigest()
        if observed != expected:
            _stop("PREF1_AUTHORITY_BLOB_DRIFT", relative)


def _expected_manifest() -> dict[str, Any]:
    return {
        "TI2_PREF1_result_sha256": TI2_PREF1_RESULT_SHA256,
        "actual_spatial_operator": "inherited_RK4_2049_SBP4",
        "artifact_id": AC1_ARTIFACT_ID,
        "authority_commit": AUTHORITY_COMMIT,
        "experiment_label": EXPERIMENT_LABEL,
        "output_leaves": ["manifest.json", "terminal.json"],
        "production_SSPRK3_comparator": False,
        "retry": RETRY,
        "runner_id": RUNNER_ID,
        "schema": RAW_SCHEMA,
        "tableau_selector": "SSPRK3",
    }


def _expected_terminal_base() -> dict[str, Any]:
    sealed = {
        "leaf_count": SEALED_STORE_LEAF_COUNT,
        "sha256": SEALED_STORE_SHA256,
    }
    return {
        "GR0_calibration_completed": False,
        "PDE_state_committed": False,
        "actual_spatial_operator": "inherited_RK4_2049_SBP4",
        "all_pass_licenses_only_this_frozen_state_and_width": True,
        "artifact_id": AC1_ARTIFACT_ID,
        "attempted_width_hex": WIDTH_HEX,
        "authority_commit": AUTHORITY_COMMIT,
        "candidate_branch_opened": False,
        "channel_order": list(CHANNEL_ORDER),
        "classification": RAW_CLASSIFICATION,
        "common_event_completed": False,
        "experiment_label": EXPERIMENT_LABEL,
        "fail_licenses_only_diagnosis_of_failed_channel_set": True,
        "fine_path_committed": False,
        "independent_method_agreement": False,
        "mechanism_result_earned": False,
        "physical_result_earned": False,
        "production_SSPRK3_comparator": False,
        "production_method_earned": False,
        "retry": RETRY,
        "runner_id": RUNNER_ID,
        "schema": RAW_SCHEMA,
        "store_snapshot_after": sealed,
        "store_snapshot_before": sealed,
        "store_unchanged": True,
        "successor_remedy_selected": False,
        "tableau_runtime_selector": TABLEAU_RUNTIME_SELECTOR,
        "tableau_selector": "SSPRK3",
        "temporal_retry_admission_called": False,
    }


def _validate_replay_receipt(value: object) -> dict[str, Any]:
    receipt = _mapping(value, frozenset(EXPECTED_REPLAY_RECEIPT), "replay_receipt")
    observed = {key: receipt[key] for key in EXPECTED_REPLAY_RECEIPT}
    if observed != EXPECTED_REPLAY_RECEIPT:
        _stop("PREF1_REPLAY_RECEIPT_DRIFT", observed)
    return dict(EXPECTED_REPLAY_RECEIPT)


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
    metadata = root.lstat()
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        _stop("PREF1_REPOSITORY_UNSAFE", root)
    manifest_raw, terminal_raw, manifest, terminal = _raw_snapshot(root)
    if manifest != _expected_manifest():
        _stop("PREF1_RAW_MANIFEST_DRIFT", "manifest schema differs")
    _authenticate_authority(root, manifest)
    before_store = _snapshot_store(root)
    sealed_store = (SEALED_STORE_LEAF_COUNT, SEALED_STORE_SHA256)
    if before_store != sealed_store:
        _stop("PREF1_STORE_IDENTITY_DRIFT", before_store)
    base = _expected_terminal_base()
    if set(terminal) != set(base) | _COMPLETED_EXTRA_KEYS:
        _stop("PREF1_TERMINAL_SCHEMA_DRIFT", sorted(terminal))
    if any(terminal.get(name) != value for name, value in base.items()):
        _stop("PREF1_TERMINAL_NONCLAIM_DRIFT", terminal.get("classification"))
    receipt = _validate_replay_receipt(terminal.get("replay_receipt"))
    if (
        terminal.get("shadow_path_count") != 7
        or terminal.get("shadow_proposal_count") != 7
        or terminal.get("SSPRK3_stage_and_endpoint_record_count") != 28
    ):
        _stop("PREF1_SHADOW_COUNT_DRIFT", "completed terminal counts")
    admission = terminal.get("continuous_admission")
    if not isinstance(admission, Mapping):
        _stop("PREF1_ADMISSION_TYPE_DRIFT", type(admission).__name__)
    replay = _reclassify_admission(admission)
    independent_class = reduce_retry3_terminal(
        complete_admission_passed=bool(admission["complete_admission_passed"]),
        failed_channels=tuple(admission["failed_channels"]),
    )
    if (
        independent_class != RAW_CLASSIFICATION
        or terminal.get("classification") != independent_class
    ):
        _stop("PREF1_TERMINAL_CLASS_DRIFT", independent_class)
    after_store = _snapshot_store(root)
    manifest_after, terminal_after = _raw_pair(root)
    if (
        after_store != before_store
        or manifest_after != manifest_raw
        or terminal_after != terminal_raw
    ):
        _stop("PREF1_INPUT_MUTATED_DURING_BIND", "raw/store input changed")
    return {
        **expected,
        "raw_binding": {
            "canonical_duplicate_free_schema_verified": True,
            "leaf_count": 2,
            "manifest_sha256": sha256(manifest_raw).hexdigest(),
            "raw_namespace_unchanged": True,
            "terminal_sha256": sha256(terminal_raw).hexdigest(),
        },
        "authority_binding": {
            "all_committed_blob_hashes_verified": True,
            "blob_count": len(_AUTHORITY_BLOBS),
            "commit": AUTHORITY_COMMIT,
            "parent_commit": AUTHORITY_PARENT,
        },
        "store_binding": {
            "binder_snapshot_after": {
                "leaf_count": after_store[0],
                "sha256": after_store[1],
            },
            "binder_snapshot_before": {
                "leaf_count": before_store[0],
                "sha256": before_store[1],
            },
            "raw_snapshot_after": terminal["store_snapshot_after"],
            "raw_snapshot_before": terminal["store_snapshot_before"],
            "store_unchanged": True,
        },
        "independent_replay": {
            **replay,
            "SSPRK3_stage_and_endpoint_record_count": 28,
            "replay_receipt": receipt,
            "shadow_path_count": 7,
            "shadow_proposal_count": 7,
        },
        "conclusion": {
            "PDE_or_candidate_state_opened": False,
            "frozen_SSPRK3_on_inherited_SBP4_state_width_has_one_unresolved_u_R_channel": True,
            "independent_method_agreement": False,
            "licenses_only_separately_frozen_diagnosis_of_u_R_or_sharper_discriminator": True,
            "physics_inference_permitted": False,
            "production_SSPRK3_comparator": False,
            "production_method_earned": False,
            "state_advance_authorized": False,
            "successor_remedy_selected": False,
        },
    }


def build_pref1_result(
    config_raw: bytes, repository: Path, *, live: bool = True
) -> dict[str, Any]:
    payload = (
        bind_raw_result(config_raw, repository)
        if live
        else expected_evidence(config_raw)
    )
    return {"artifact_id": ARTIFACT_ID, "artifact_payload": payload}


def _expected_channel_summaries() -> list[dict[str, Any]]:
    return [
        {
            "admission_passed": admitted,
            "channel": channel,
            "classification": classification,
            "debit_equals_finest_upper": True,
            "order_threshold_passed": order_passed,
            "order_threshold_resolved": order_resolved,
            "temporal_retry_permitted": retry_permitted,
        }
        for (
            channel,
            classification,
            admitted,
            retry_permitted,
            order_resolved,
            order_passed,
        ) in EXPECTED_CHANNEL_CLASSES
    ]


def validate_compact_result(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    expected = expected_evidence(config_raw)
    if sha256(result_raw).hexdigest() != COMPACT_RESULT_SHA256:
        _stop("PREF1_COMPACT_DRIFT", "compact result hash differs")
    result = _json(result_raw, "compact result")
    if (
        result.get("artifact_id") != ARTIFACT_ID
        or not isinstance(result.get("artifact_payload"), Mapping)
    ):
        _stop("PREF1_COMPACT_DRIFT", "compact identity differs")
    payload = result["artifact_payload"]
    for key, value in expected.items():
        if payload.get(key) != value:
            _stop("PREF1_COMPACT_DRIFT", key)
    if payload.get("raw_binding") != {
        "canonical_duplicate_free_schema_verified": True,
        "leaf_count": 2,
        "manifest_sha256": RAW_MANIFEST_SHA256,
        "raw_namespace_unchanged": True,
        "terminal_sha256": RAW_TERMINAL_SHA256,
    }:
        _stop("PREF1_COMPACT_DRIFT", "raw binding")
    if payload.get("authority_binding") != {
        "all_committed_blob_hashes_verified": True,
        "blob_count": len(_AUTHORITY_BLOBS),
        "commit": AUTHORITY_COMMIT,
        "parent_commit": AUTHORITY_PARENT,
    }:
        _stop("PREF1_COMPACT_DRIFT", "authority binding")
    expected_snapshot = {
        "leaf_count": SEALED_STORE_LEAF_COUNT,
        "sha256": SEALED_STORE_SHA256,
    }
    if payload.get("store_binding") != {
        "binder_snapshot_after": expected_snapshot,
        "binder_snapshot_before": expected_snapshot,
        "raw_snapshot_after": expected_snapshot,
        "raw_snapshot_before": expected_snapshot,
        "store_unchanged": True,
    }:
        _stop("PREF1_COMPACT_DRIFT", "store binding")
    replay = payload.get("independent_replay")
    if not isinstance(replay, Mapping):
        _stop("PREF1_COMPACT_DRIFT", "independent replay absent")
    expected_replay = {
        "SSPRK3_stage_and_endpoint_record_count": 28,
        "all_of_identity_verified": True,
        "admitted_count": EXPECTED_ADMITTED_COUNT,
        "channel_count": len(CHANNEL_ORDER),
        "channels": _expected_channel_summaries(),
        "classification_counts": dict(EXPECTED_CLASSIFICATION_COUNTS),
        "complete_admission_passed": False,
        "debit_equals_finest_upper_count": len(CHANNEL_ORDER),
        "failed_channels": list(EXPECTED_FAILED_CHANNELS),
        "failed_channels_match_admission_flags": True,
        "failed_count": len(EXPECTED_FAILED_CHANNELS),
        "raw_terminal_class_matches_independent_reduction": True,
        "replay_receipt": dict(EXPECTED_REPLAY_RECEIPT),
        "shadow_path_count": 7,
        "shadow_proposal_count": 7,
    }
    if replay != expected_replay:
        _stop("PREF1_COMPACT_DRIFT", "independent replay")
    if payload.get("conclusion") != {
        "PDE_or_candidate_state_opened": False,
        "frozen_SSPRK3_on_inherited_SBP4_state_width_has_one_unresolved_u_R_channel": True,
        "independent_method_agreement": False,
        "licenses_only_separately_frozen_diagnosis_of_u_R_or_sharper_discriminator": True,
        "physics_inference_permitted": False,
        "production_SSPRK3_comparator": False,
        "production_method_earned": False,
        "state_advance_authorized": False,
        "successor_remedy_selected": False,
    }:
        _stop("PREF1_COMPACT_DRIFT", "conclusion")
    return result


__all__ = (
    "ARTIFACT_ID",
    "CONFIG_PATH",
    "OWNER_DOCUMENT",
    "RESULT_PATH",
    "TDG9AC1PREF1Error",
    "bind_raw_result",
    "build_pref1_result",
    "canonical_result",
    "expected_evidence",
    "reduce_retry3_terminal",
    "validate_compact_result",
)
