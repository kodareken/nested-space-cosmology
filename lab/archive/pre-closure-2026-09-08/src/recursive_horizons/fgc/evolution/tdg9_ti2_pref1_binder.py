"""Independent post-result binder for the completed TDG9 TI2 shadow.

The one-time live path authenticates the immutable TI2 authority and raw
terminal, restores the three sealed RK4-2049 retry predecessors, constructs
the same-state SSPRK3-tableau shadows through the lower-level TDG6 API, and
independently repeats all 60 primary/v2 exact localization decisions.  It
deliberately does not import the TI2 runner.  Ordinary verification consumes
only the tracked compact certificate and is raw/store blind.
"""

from __future__ import annotations

from copy import deepcopy
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
from typing import Any, Mapping, NoReturn, Sequence

import numpy as np

from . import hlt16_campaign_runtime as campaign_runtime
from . import proto15_runtime as p15
from . import tdg6_temporal_admission_runtime as tdg6
from . import tdg9_ar1_authority as ar1
from . import tdg9_ar1_pref1_binder as ar1_pref1
from . import tdg9_loc2_pref2_binder as loc2_pref2
from . import tdg9_local_extrema as primary
from . import tdg9_local_extrema_independent_v2 as independent_v2
from . import tdg9_ti2_authority as ti2
from .hlt16_campaign_store import HLT16CampaignStore
from .numerical_engine import COMPARATOR_METHOD, PRIMARY_METHOD, array_content_sha256
from .proto19_gr0_static_factory import build_static_gr0_shells


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG9-TI2-PREF1"
CLASSIFICATION = "independently_bound_same_state_tableau_shadow_terminal"
CONFIG_PATH = "configs/fgc/fgc-1-tdg9-ti2-pref1.toml"
RESULT_PATH = "results/fgc-1-tdg9-ti2-pref1.json"
OWNER_DOCUMENT = "docs/fgc-tdg9-ti2-pref1.md"
TARGET_PROTOCOL = ti2.TARGET_PROTOCOL

AUTHORITY_COMMIT = "963f772c7fd5de676c37c3ee385a930701df6341"
AUTHORITY_PARENT = ti2.BASE_COMMIT
RAW_NAMESPACE = ti2.OUTPUT_NAMESPACE
RAW_SCHEMA = "FGC-1-TDG9-TI2-raw-v1"
RAW_CLASSIFICATION = (
    "completed_one_or_more_complete_failures_persist_under_tableau_shadow"
)
RUNNER_ID = "FGC-1-TDG9-TI2-RUN1"
RAW_MANIFEST_SHA256 = "afa3d26975c116dfdd350e409c39462769f8b8b5e96e22902496a7704b43e883"
RAW_TERMINAL_SHA256 = "e5e8dd4e8d0dee3d32fefc170a5ed0e4f07c0c249200d74d4296ee9a3a601311"
COMPACT_RESULT_SHA256 = "4a7a4bb15c35766b9a73dd734f657aeae2c99692c37b00d7a9697dc369ccf546"

SEALED_STORE_LEAF_COUNT = 115
SEALED_STORE_SHA256 = "5917d70de3080dcc6b7abeb7fdb7cc3370c6140cbeb37bde0c142d6a2d36a445"
EXPECTED_COMPONENT_CANDIDATES = {
    "complete_C": 332_438,
    "slope_S": 397_721,
    "value_V": 245_280,
}
EXPECTED_TOTAL_CANDIDATES = 975_439
EXPECTED_ROUTE_CLASSIFICATIONS = {
    "nonunique_or_interval_inconclusive": 10,
    "unique_maximum": 50,
}
EXPECTED_SURVIVOR_COUNT = 93
EXPECTED_ROOT_METHODS = {
    "adaptive_nonsquare": 16,
    "exact_endpoint": 60,
    "exact_endpoint_deflation": 17,
}
EXPECTED_OWNERSHIP_COUNTS = {
    "both_components_independently_fail": 4,
    "complete_not_failure": 4,
    "endpoint_state_owned_failure": 2,
}
EXPECTED_TRANSITIONS = {
    "both_components_independently_fail->both_components_independently_fail": 4,
    "both_components_independently_fail->complete_not_failure": 4,
    "both_components_independently_fail->endpoint_state_owned_failure": 2,
}
EXPECTED_COMPLETE_COUNTS = {"failure": 6, "inconclusive": 0, "pass": 4}
EXPECTED_OCCURRENCE_CLASSES = (
    (3, "u:alpha", "complete_not_failure", "sufficient_contraction_pass", True),
    (3, "u:R", "complete_not_failure", "sufficient_contraction_pass", True),
    (
        4,
        "u:alpha",
        "endpoint_state_owned_failure",
        "sufficient_contraction_failure",
        False,
    ),
    (
        4,
        "u:lambda",
        "endpoint_state_owned_failure",
        "sufficient_contraction_failure",
        False,
    ),
    (
        4,
        "u:R",
        "both_components_independently_fail",
        "sufficient_contraction_failure",
        False,
    ),
    (
        5,
        "u:alpha",
        "both_components_independently_fail",
        "sufficient_contraction_failure",
        False,
    ),
    (5, "u:v", "complete_not_failure", "sufficient_contraction_pass", True),
    (
        5,
        "u:lambda",
        "both_components_independently_fail",
        "sufficient_contraction_failure",
        False,
    ),
    (
        5,
        "u:R",
        "both_components_independently_fail",
        "sufficient_contraction_failure",
        False,
    ),
    (5, "q:R", "complete_not_failure", "sufficient_contraction_pass", True),
)
_MAX_RAW_LEAF_BYTES = 8 * 1024 * 1024
_MAX_REPO_LEAF_BYTES = 4 * 1024 * 1024

_AUTHORITY_BLOBS = (
    (ti2.CONFIG_PATH, "8631ce20e22a56b00f4a9580df1ec2e2e98176b14c4761dc78e316cb4a771546"),
    (ti2.RESULT_PATH, "2f2f2f844e9611b1023cee869e9f01ce33cfdef2df01cf488302eaaa28a1e3f5"),
    (ti2.OWNER_DOCUMENT, "91f0f9fbf410e6e896d12c5bf01d0c7105f6681114caded73ad21c0a1f668080"),
    (ti2.RUNNER_PATH, "488f05b65b71dd2617d1d77c2efa49496aa1a62698b37866e07bbb9c5792cd57"),
    (
        "src/recursive_horizons/fgc/evolution/tdg9_ti2_authority.py",
        "7715685100e56434a36edbec73ca684fe2fd88010af3eaeae1a59761884512b5",
    ),
    (
        "scripts/reproduce_fgc_tdg9_ti2_frz1.py",
        "42bef89656720404668e3d4daab60295d7b06b9a85d79410bc197de41ece97c4",
    ),
)


class TDG9TI2PREF1Error(ValueError):
    """The compact contract, raw terminal, store, or independent replay differs."""

    def __init__(self, stop_id: str, detail: object) -> None:
        self.stop_id = str(stop_id)
        self.detail = " ".join(str(detail).split())[:640]
        super().__init__(f"{self.stop_id}: {self.detail}")


def _stop(stop_id: str, detail: object) -> NoReturn:
    raise TDG9TI2PREF1Error(stop_id, detail)


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
        raise TDG9TI2PREF1Error("PREF1_CANONICAL_DRIFT", exc) from exc


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
        raise TDG9TI2PREF1Error("PREF1_COMPACT_DRIFT", exc) from exc


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
        raise TDG9TI2PREF1Error("PREF1_JSON_DRIFT", label) from exc
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
        raise TDG9TI2PREF1Error("PREF1_GIT_DRIFT", arguments) from exc


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


def _read_repo_leaf(root: Path, relative: str, maximum: int) -> bytes:
    candidate = Path(relative)
    parent = candidate.parent.as_posix()
    descriptor = (
        os.open(
            root,
            os.O_RDONLY
            | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_NOFOLLOW", 0),
        )
        if parent == "."
        else _open_directory(root, parent)
    )
    try:
        return _read_leaf_at(
            descriptor, candidate.name, maximum=maximum, label=relative
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
    try:
        return loc2_pref2._snapshot_store(root)  # type: ignore[attr-defined]
    except Exception as exc:
        raise TDG9TI2PREF1Error("PREF1_STORE_SNAPSHOT_DRIFT", exc) from exc


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
            "experiment_label": ti2.SEMANTICS["experiment_label"],
            "tableau_selector": "SSPRK3",
            "tableau_runtime_selector": COMPARATOR_METHOD,
            "actual_spatial_operator": "inherited_RK4_2049_SBP4",
            "production_SSPRK3_comparator": False,
            "independent_method_agreement": False,
            "member_key": ti2.MEMBER_KEY,
            "physical_state_sha256": ti2.PHYSICAL_STATE_SHA256,
            "accepted_time_hex": ti2.ACCEPTED_TIME_HEX,
            "point_count": ti2.POINT_COUNT,
            "owned_row_count": ti2.OWNED_ROW_COUNT,
            "retries": [3, 4, 5],
            "failed_occurrences": [
                {"retry": retry, "channel": channel}
                for retry, channel in ti2.FAILED_OCCURRENCES
            ],
        },
        "work_budget": {
            "base_cubic_count": ti2.BASE_CUBIC_COUNT,
            "component_cubic_count": ti2.COMPONENT_CUBIC_COUNT,
            "family_candidate_ceiling": ti2.FAMILY_CANDIDATE_CEILING,
            "aggregate_candidate_ceiling": ti2.AGGREGATE_CANDIDATE_CEILING,
            "primary_refinement_depth": ti2.PRIMARY_REFINEMENT_DEPTH,
            "v2_global_proof_bit_ceiling": ti2.GLOBAL_PROOF_BIT_CEILING,
            "shadow_proposal_count": 21,
            "SSPRK3_stage_and_endpoint_record_count": 84,
            "route_count": 60,
        },
        "scope": {
            "one_time_live_raw_and_store_reconstruction": True,
            "ordinary_verifier_raw_blind": True,
            "ordinary_verifier_store_blind": True,
            "TI2_runner_decision_code_imported": False,
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
            "TI2_raw_result_independently_bound": True,
            "same_state_tableau_sensitivity_measured": True,
            "all_ten_complete_failures_cleared": False,
            "one_or_more_complete_failures_persist": True,
            "shadow_complete_failure_count": 6,
            "shadow_complete_pass_count": 4,
            "shadow_complete_inconclusive_count": 0,
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
        raise TDG9TI2PREF1Error("PREF1_CONFIG_DRIFT", exc) from exc
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
        "LOC2_PREF2_result_sha256": ti2.LOC2_PREF2_RESULT_SHA256,
        "actual_spatial_operator": "inherited_RK4_2049_SBP4",
        "artifact_id": ti2.ARTIFACT_ID,
        "authority_commit": AUTHORITY_COMMIT,
        "experiment_label": ti2.SEMANTICS["experiment_label"],
        "output_leaves": ["manifest.json", "terminal.json"],
        "production_SSPRK3_comparator": False,
        "runner_id": RUNNER_ID,
        "schema": RAW_SCHEMA,
        "tableau_selector": "SSPRK3",
    }


def _original_occurrences(root: Path) -> dict[tuple[int, str], dict[str, object]]:
    config_raw = _read_repo_leaf(root, ti2.LOC2_PREF2_CONFIG_PATH, 1 << 20)
    result_raw = _read_repo_leaf(root, ti2.LOC2_PREF2_RESULT_PATH, 8 << 20)
    if (
        sha256(config_raw).hexdigest() != ti2.LOC2_PREF2_CONFIG_SHA256
        or sha256(result_raw).hexdigest() != ti2.LOC2_PREF2_RESULT_SHA256
    ):
        _stop("PREF1_LOC2_BINDING_DRIFT", "tracked compact predecessor differs")
    result = loc2_pref2.validate_compact_result(config_raw, result_raw)
    answer: dict[tuple[int, str], dict[str, object]] = {}
    for item in result["artifact_payload"]["independent_replay"]["occurrences"]:
        key = (int(item["retry"]), str(item["channel"]))
        complete = {
            str(component["contraction"])
            for component in item["components"]
            if component["component"] == "complete_C"
        }
        if complete != {"sufficient_contraction_failure"}:
            _stop("PREF1_ORIGINAL_FAILURE_DRIFT", key)
        answer[key] = {
            "original_complete_classification": "sufficient_contraction_failure",
            "original_complete_failure": True,
            "original_ownership_class": item["component_ownership"]["classification"],
            "original_ownership_sha256": item["component_ownership_sha256"],
            "original_row_stream_sha256": item["PREF1_row_stream_sha256"],
        }
    if tuple(answer) != ti2.FAILED_OCCURRENCES:
        _stop("PREF1_OCCURRENCE_SELECTION_DRIFT", tuple(answer))
    return answer


def _class_name(value: object) -> str:
    kind = type(value)
    return f"{kind.__module__}.{kind.__qualname__}"


def _callable_name(value: object) -> str:
    return (
        f"{getattr(value, '__module__', type(value).__module__)}."
        f"{getattr(value, '__qualname__', type(value).__qualname__)}"
    )


def _fingerprint_exact(value: object) -> object:
    if isinstance(value, Fraction):
        return {
            "numerator": str(value.numerator),
            "denominator": str(value.denominator),
        }
    if type(value) is float:
        if not math.isfinite(value):
            _stop("PREF1_FINGERPRINT_DRIFT", "nonfinite binary64")
        return {"binary64_hex": value.hex()}
    if is_dataclass(value) and not isinstance(value, type):
        return _fingerprint_exact(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _fingerprint_exact(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_fingerprint_exact(item) for item in value]
    if value is None or type(value) in {str, int, bool}:
        return value
    _stop("PREF1_FINGERPRINT_DRIFT", type(value).__name__)


def _tracer_history_hash(member: object) -> str:
    arrays = [
        np.asarray(member.tracers.labels, dtype=np.float64),
        np.asarray(member.tracers.positions, dtype=np.float64),
        np.asarray(member.tracers.proper_times, dtype=np.float64),
    ]
    arrays.extend(
        np.asarray(row, dtype=np.float64) for row in member.tracers.event_proper_times
    )
    arrays.extend(
        np.asarray(row, dtype=np.float64) for row in member.tracers.event_fields
    )
    return array_content_sha256(*arrays)


def _member_fingerprint(member: object, descriptor_sha256: str) -> dict[str, object]:
    ledger = member.temporal_ledger
    if ledger is None:
        _stop("PREF1_REPLAY_DRIFT", "temporal ledger absent")
    return {
        "member_key": member.key,
        "method_label": member.method_label,
        "source_integrator": member.integrator_id,
        "point_count": member.point_count,
        "accepted_time_hex": float(member.time).hex(),
        "state_sha256": array_content_sha256(
            member.state.u, member.state.p, member.state.q
        ),
        "coordinates_sha256": array_content_sha256(member.initial.grid.coordinates),
        "grid_spacing_hex": float(member.initial.grid.spacing).hex(),
        "outer_radius_hex": float(member.initial.grid.maximum).hex(),
        "descriptor_sha256": descriptor_sha256,
        "operator_class": _class_name(member.operator),
        "projector_callable": _callable_name(member.projector),
        "transaction_class": _class_name(member.transaction),
        "tracer_class": _class_name(member.tracers),
        "ledger_class": _class_name(ledger),
        "transaction_sha256": sha256(
            _canonical(
                _fingerprint_exact(
                    {
                        "monitor": asdict(member.transaction.state),
                        "causal": asdict(member.transaction.causal_state),
                        "ledger": asdict(ledger),
                        "step_index": member.step_index,
                        "transaction_serial": member.transaction_serial,
                    }
                )
            )
        ).hexdigest(),
        "tracer_history_sha256": _tracer_history_hash(member),
    }


def _restore_shadow(
    root: Path,
    store: HLT16CampaignStore,
    shells: Mapping[str, object],
    replay: Mapping[str, object],
) -> tuple[tdg6.TDG6PreparedGR0Compositor, dict[str, object]]:
    generation = int(replay["predecessor_generation"])
    checkpoint = store.authenticated_checkpoint_at_generation(generation)
    if checkpoint.sha256 != replay["checkpoint_sha256"]:
        _stop("PREF1_CHECKPOINT_DRIFT", generation)
    relative = (
        f"{ar1.PREF2_STORE_PATH}/checkpoints/"
        f"{generation:020d}-{checkpoint.sha256}.json"
    )
    if sha256(_read_repo_leaf(root, relative, 8 << 20)).hexdigest() != replay[
        "checkpoint_raw_sha256"
    ]:
        _stop("PREF1_CHECKPOINT_RAW_DRIFT", generation)
    state = checkpoint.members.get(ti2.MEMBER_KEY)
    if state is None:
        _stop("PREF1_MEMBER_ABSENT", ti2.MEMBER_KEY)
    member = deepcopy(shells[ti2.MEMBER_KEY])
    try:
        campaign_runtime.restore_member_with_overlay(
            store, checkpoint, member, key=ti2.MEMBER_KEY
        )
        cursor = p15.Proto15Cursor(dict(state.cursor))
        cursor.validate()
    except Exception as exc:
        raise TDG9TI2PREF1Error("PREF1_RESTORE_DRIFT", replay["retry"]) from exc
    ledger = member.temporal_ledger
    if (
        ledger is None
        or member.key != ti2.MEMBER_KEY
        or member.method_label != "RK4"
        or member.integrator_id != PRIMARY_METHOD
        or member.point_count != ti2.POINT_COUNT
        or float(member.time).hex() != ti2.ACCEPTED_TIME_HEX
        or ledger.current_macro_step_temporal_retry_count
        != replay["prior_retry_count"]
        or cursor.mode != "RETRY_PENDING"
        or state.pending_owner != "temporal"
    ):
        _stop("PREF1_PREDECESSOR_DRIFT", replay["retry"])
    fingerprint = _member_fingerprint(member, state.descriptor_sha256)
    retry = int(replay["retry"])
    fixed = {
        "member_key": ti2.MEMBER_KEY,
        "method_label": "RK4",
        "source_integrator": PRIMARY_METHOD,
        "point_count": ti2.POINT_COUNT,
        "accepted_time_hex": ti2.ACCEPTED_TIME_HEX,
        "state_sha256": ti2.PHYSICAL_STATE_SHA256,
        "coordinates_sha256": ti2.COORDINATES_SHA256,
        "grid_spacing_hex": ti2.GRID_SPACING_HEX,
        "outer_radius_hex": ti2.OUTER_RADIUS_HEX,
        "descriptor_sha256": ti2.MEMBER_DESCRIPTOR_SHA256,
        "transaction_sha256": ti2.TRANSACTION_SHA256_BY_RETRY[retry],
    }
    if any(fingerprint.get(name) != expected for name, expected in fixed.items()):
        _stop("PREF1_FINGERPRINT_DRIFT", retry)
    before = dict(fingerprint)
    try:
        prepared = tdg6.prepare_tdg6_gr0_compositor(
            method=COMPARATOR_METHOD,
            time=member.time,
            step_size=float.fromhex(str(replay["attempted_width_hex"])),
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
        raise TDG9TI2PREF1Error("PREF1_SHADOW_PREMISE_STOP", retry) from exc
    after = _member_fingerprint(member, state.descriptor_sha256)
    paths = (prepared.outer, prepared.medium, prepared.fine)
    proposals = tuple(attempt.proposal for path in paths for attempt in path.attempts)
    if (
        after != before
        or prepared.method != COMPARATOR_METHOD
        or prepared.initial_state_sha256 != ti2.PHYSICAL_STATE_SHA256
        or tuple(len(path.attempts) for path in paths) != (1, 2, 4)
        or len(proposals) != 7
        or any(proposal.method != COMPARATOR_METHOD for proposal in proposals)
        or any(len(proposal.stages) != 4 for proposal in proposals)
    ):
        _stop("PREF1_ONE_VARIABLE_IDENTITY_DRIFT", retry)
    historical = ar1_pref1._journal_evidence(root, replay)  # type: ignore[attr-defined]
    receipt = {
        "predecessor": fingerprint,
        "historical_journal_sha256": sha256(_canonical(historical)).hexdigest(),
        "attempted_width_hex": replay["attempted_width_hex"],
        "shadow_path_count": 7,
        "shadow_proposal_count": 7,
        "SSPRK3_stage_and_endpoint_record_count": 28,
    }
    return prepared, receipt


def _component_summary(
    *,
    rows: Sequence[object],
    component: str,
    level: str,
    raw: Mapping[str, Any],
    evidence: dict[str, dict[str, object]],
) -> tuple[dict[str, Any], dict[str, int]]:
    expected = (2 if level == "D01" else 4) * ti2.OWNED_ROW_COUNT
    first_cubics, second_cubics = loc2_pref2._cubics(  # type: ignore[attr-defined]
        rows, level, component
    )
    if len(first_cubics) != expected or len(second_cubics) != expected:
        _stop("PREF1_POLYNOMIAL_COUNT_DRIFT", (component, level))
    first = primary.localize_absolute_maximum(
        first_cubics,
        maximum_candidates=4 * expected,
        refinement_depth=ti2.PRIMARY_REFINEMENT_DEPTH,
    )
    try:
        second = independent_v2.localize_absolute_maximum_independently_v2(
            second_cubics, maximum_candidates=4 * expected
        )
    except independent_v2.RootIsolationInconclusive as exc:
        raise TDG9TI2PREF1Error(
            "PREF1_V2_LOCALIZATION_INCONCLUSIVE", (component, level, exc.reason)
        ) from exc
    stationary = loc2_pref2._primary_stationary_count_digest(  # type: ignore[attr-defined]
        first_cubics
    )
    agreement = loc2_pref2._route_agreement(  # type: ignore[attr-defined]
        first, second, stationary, (component, level)
    )
    first_json = loc2_pref2._exact_json(first)  # type: ignore[attr-defined]
    second_json = loc2_pref2._exact_json(second)  # type: ignore[attr-defined]
    expected_raw = {
        "component": component,
        "level": level,
        "polynomial_count": expected,
        "candidate_count": first.candidate_count,
        "primary_stationary_count_stream_sha256": stationary,
        "independent_stationary_count_stream_sha256": (
            second.stationary_count_stream_sha256
        ),
        "primary": first_json,
        "independent_v2": second_json,
        "survivor_keys_and_intervals_compared": True,
    }
    if any(raw.get(key) != value for key, value in expected_raw.items()):
        _stop("PREF1_RAW_COMPONENT_DRIFT", (component, level))
    evidence[component][level] = first
    root_methods: dict[str, int] = {}
    for item in second.candidates:
        root_methods[item.root_method] = root_methods.get(item.root_method, 0) + 1
    return (
        {
            "component": component,
            "level": level,
            "polynomial_count": expected,
            "candidate_count": first.candidate_count,
            "classification": first.classification,
            "co_maximizer_count": len(first.candidates),
            "stationary_count_stream_sha256": stationary,
            "route_agreement": agreement,
            "primary_evidence_sha256": sha256(_canonical(first_json)).hexdigest(),
            "independent_v2_evidence_sha256": sha256(
                _canonical(second_json)
            ).hexdigest(),
        },
        root_methods,
    )


def _occurrence_summary(
    *,
    rows: Sequence[object],
    retry: int,
    channel: str,
    original: Mapping[str, object],
    raw: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, int]]:
    if raw.get("retry") != retry or raw.get("channel") != channel:
        _stop("PREF1_OCCURRENCE_ORDER_DRIFT", (retry, channel))
    if any(raw.get(name) != value for name, value in original.items()):
        _stop("PREF1_ORIGINAL_OCCURRENCE_DRIFT", (retry, channel))
    row_hash = loc2_pref2._row_hash(rows)  # type: ignore[attr-defined]
    if raw.get("row_stream_sha256") != row_hash:
        _stop("PREF1_ROW_STREAM_DRIFT", (retry, channel))
    raw_components = raw.get("components")
    expected_order = tuple(
        (component, level)
        for component in ("value_V", "slope_S", "complete_C")
        for level in ("D01", "D12")
    )
    if not isinstance(raw_components, list) or len(raw_components) != 6:
        _stop("PREF1_COMPONENT_SHAPE_DRIFT", (retry, channel))
    evidence: dict[str, dict[str, object]] = {
        name: {} for name in ("value_V", "slope_S", "complete_C")
    }
    summaries: list[dict[str, Any]] = []
    methods: dict[str, int] = {}
    for raw_component, (component, level) in zip(
        raw_components, expected_order, strict=True
    ):
        if not isinstance(raw_component, Mapping):
            _stop("PREF1_COMPONENT_SHAPE_DRIFT", (retry, channel, component, level))
        summary, counts = _component_summary(
            rows=rows,
            component=component,
            level=level,
            raw=raw_component,
            evidence=evidence,
        )
        summaries.append(summary)
        for name, count in counts.items():
            methods[name] = methods.get(name, 0) + count
    for component in ("value_V", "slope_S", "complete_C"):
        contraction = loc2_pref2._contraction_assessment(  # type: ignore[attr-defined]
            evidence[component]["D01"], evidence[component]["D12"]
        )
        evidence[component]["contraction"] = contraction
        for summary, raw_component in zip(summaries, raw_components, strict=True):
            if summary["component"] == component:
                if raw_component.get("contraction") != contraction:
                    _stop("PREF1_CONTRACTION_DRIFT", (retry, channel, component))
                summary["contraction"] = contraction["classification"]
                summary["contraction_evidence_sha256"] = sha256(
                    _canonical(contraction)
                ).hexdigest()
    ownership = loc2_pref2._component_ownership(evidence)  # type: ignore[attr-defined]
    complete = str(evidence["complete_C"]["contraction"]["classification"])  # type: ignore[index]
    failed = complete == "sufficient_contraction_failure"
    passed = complete in {"sufficient_contraction_pass", "exact_zero"}
    inconclusive = complete == "threshold_inconclusive"
    expected_tail = {
        "shadow_complete_classification": complete,
        "shadow_complete_failure": failed,
        "shadow_complete_pass": passed,
        "shadow_complete_inconclusive": inconclusive,
        "shadow_ownership": ownership,
        "shadow_ownership_class": ownership["classification"],
        "complete_failure_cleared": passed,
        "ownership_class_changed": (
            original["original_ownership_class"] != ownership["classification"]
        ),
    }
    if any(raw.get(name) != value for name, value in expected_tail.items()):
        _stop("PREF1_OCCURRENCE_REDUCTION_DRIFT", (retry, channel))
    return (
        {
            "retry": retry,
            "channel": channel,
            "row_stream_sha256": row_hash,
            "raw_occurrence_sha256": sha256(_canonical(raw)).hexdigest(),
            "original_ownership_class": original["original_ownership_class"],
            "shadow_ownership_class": ownership["classification"],
            "shadow_ownership_sha256": sha256(_canonical(ownership)).hexdigest(),
            "shadow_complete_classification": complete,
            "complete_failure_cleared": passed,
            "components": summaries,
        },
        methods,
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
    originals = _original_occurrences(root)
    raw_occurrences = terminal.get("occurrences")
    if not isinstance(raw_occurrences, list) or len(raw_occurrences) != 10:
        _stop("PREF1_RAW_OCCURRENCE_SHAPE_DRIFT", type(raw_occurrences).__name__)
    store = HLT16CampaignStore(root / ar1.PREF2_STORE_PATH)
    shells = build_static_gr0_shells(root)
    surfaces: dict[int, tuple[tuple[Any, ...], ...]] = {}
    replay_receipts: dict[str, object] = {}
    for replay in ti2.REPLAYS:
        prepared, receipt = _restore_shadow(root, store, shells, replay)
        retry = int(replay["retry"])
        surfaces[retry] = ar1_pref1._surface(prepared)  # type: ignore[attr-defined]
        replay_receipts[str(retry)] = receipt
    if terminal.get("replay_receipts") != replay_receipts:
        _stop("PREF1_REPLAY_RECEIPT_DRIFT", "raw replay receipts differ")
    occurrences: list[dict[str, Any]] = []
    candidate_counts = {"complete_C": 0, "slope_S": 0, "value_V": 0}
    route_classes: dict[str, int] = {}
    root_methods: dict[str, int] = {}
    survivor_count = 0
    ownership_counts: dict[str, int] = {}
    complete_counts = {"failure": 0, "inconclusive": 0, "pass": 0}
    transition: dict[str, int] = {}
    route_count = 0
    for (retry, channel), raw_occurrence in zip(
        ti2.FAILED_OCCURRENCES, raw_occurrences, strict=True
    ):
        if not isinstance(raw_occurrence, Mapping):
            _stop("PREF1_RAW_OCCURRENCE_SHAPE_DRIFT", (retry, channel))
        rows = tuple(loc2_pref2._rows(surfaces[retry], channel))  # type: ignore[attr-defined]
        summary, methods = _occurrence_summary(
            rows=rows,
            retry=retry,
            channel=channel,
            original=originals[(retry, channel)],
            raw=raw_occurrence,
        )
        occurrences.append(summary)
        owner = str(summary["shadow_ownership_class"])
        ownership_counts[owner] = ownership_counts.get(owner, 0) + 1
        key = f"{summary['original_ownership_class']}->{owner}"
        transition[key] = transition.get(key, 0) + 1
        complete = str(summary["shadow_complete_classification"])
        if complete == "sufficient_contraction_failure":
            complete_counts["failure"] += 1
        elif complete in {"sufficient_contraction_pass", "exact_zero"}:
            complete_counts["pass"] += 1
        else:
            complete_counts["inconclusive"] += 1
        for item in summary["components"]:
            route_count += 1
            component = str(item["component"])
            candidate_counts[component] += int(item["candidate_count"])
            route_class = str(item["classification"])
            route_classes[route_class] = route_classes.get(route_class, 0) + 1
            survivor_count += int(item["co_maximizer_count"])
        for name, count in methods.items():
            root_methods[name] = root_methods.get(name, 0) + count
    transition = dict(sorted(transition.items()))
    total_candidates = sum(candidate_counts.values())
    aggregate = {
        "occurrence_count": len(occurrences),
        "component_level_count": route_count,
        "component_candidate_counts": candidate_counts,
        "total_VSC_candidate_count": total_candidates,
        "route_classification_counts": dict(sorted(route_classes.items())),
        "survivor_count": survivor_count,
        "v2_root_method_counts": dict(sorted(root_methods.items())),
        "shadow_complete_counts": complete_counts,
        "shadow_ownership_counts": dict(sorted(ownership_counts.items())),
        "ownership_transition_matrix": transition,
    }
    expected_aggregate = {
        "occurrence_count": 10,
        "component_level_count": 60,
        "component_candidate_counts": EXPECTED_COMPONENT_CANDIDATES,
        "total_VSC_candidate_count": EXPECTED_TOTAL_CANDIDATES,
        "route_classification_counts": EXPECTED_ROUTE_CLASSIFICATIONS,
        "survivor_count": EXPECTED_SURVIVOR_COUNT,
        "v2_root_method_counts": EXPECTED_ROOT_METHODS,
        "shadow_complete_counts": EXPECTED_COMPLETE_COUNTS,
        "shadow_ownership_counts": EXPECTED_OWNERSHIP_COUNTS,
        "ownership_transition_matrix": EXPECTED_TRANSITIONS,
    }
    if aggregate != expected_aggregate:
        _stop("PREF1_AGGREGATE_REDUCTION_DRIFT", aggregate)
    terminal_base = {
        "schema": RAW_SCHEMA,
        "artifact_id": ti2.ARTIFACT_ID,
        "runner_id": RUNNER_ID,
        "classification": RAW_CLASSIFICATION,
        "authority_commit": AUTHORITY_COMMIT,
        "experiment_label": ti2.SEMANTICS["experiment_label"],
        "tableau_selector": "SSPRK3",
        "tableau_runtime_selector": COMPARATOR_METHOD,
        "actual_spatial_operator": "inherited_RK4_2049_SBP4",
        "production_SSPRK3_comparator": False,
        "independent_method_agreement": False,
        "store_snapshot_before": {
            "leaf_count": SEALED_STORE_LEAF_COUNT,
            "sha256": SEALED_STORE_SHA256,
        },
        "store_snapshot_after": {
            "leaf_count": SEALED_STORE_LEAF_COUNT,
            "sha256": SEALED_STORE_SHA256,
        },
        "store_unchanged": True,
        "PDE_state_committed": False,
        "temporal_retry_admission_called": False,
        "fine_path_committed": False,
        "successor_remedy_selected": False,
        "common_event_completed": False,
        "GR0_calibration_completed": False,
        "candidate_branch_opened": False,
        "mechanism_result_earned": False,
        "physical_result_earned": False,
        "replay_receipts": replay_receipts,
        "occurrence_count": 10,
        "occurrences": raw_occurrences,
        "ownership_transition_matrix": transition,
        "original_complete_failure_count": 10,
        "shadow_complete_failure_count": complete_counts["failure"],
        "shadow_complete_pass_count": complete_counts["pass"],
        "shadow_complete_inconclusive_count": complete_counts["inconclusive"],
        "component_candidate_counts": candidate_counts,
        "total_VSC_candidate_count": total_candidates,
        "shadow_proposal_count": 21,
        "SSPRK3_stage_and_endpoint_record_count": 84,
        "primary_v2_route_agreement_count": 60,
    }
    if terminal != terminal_base:
        _stop("PREF1_RAW_TERMINAL_DRIFT", "terminal schema/nonclaim differs")
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
            **aggregate,
            "shadow_proposal_count": 21,
            "SSPRK3_stage_and_endpoint_record_count": 84,
            "primary_v2_route_agreement_count": 60,
            "replay_receipts_sha256": sha256(
                _canonical(replay_receipts)
            ).hexdigest(),
            "occurrences": occurrences,
        },
        "conclusion": {
            "same_state_tableau_sensitivity_measured": True,
            "four_of_ten_complete_failures_clear_under_shadow": True,
            "six_of_ten_complete_failures_persist_under_shadow": True,
            "two_persistent_failures_become_endpoint_state_owned": True,
            "production_SSPRK3_comparator": False,
            "independent_method_agreement": False,
            "state_advance_authorized": False,
            "successor_remedy_selected": False,
            "PDE_or_candidate_state_opened": False,
            "physics_inference_permitted": False,
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
    expected_aggregate = {
        "occurrence_count": 10,
        "component_level_count": 60,
        "component_candidate_counts": EXPECTED_COMPONENT_CANDIDATES,
        "total_VSC_candidate_count": EXPECTED_TOTAL_CANDIDATES,
        "route_classification_counts": EXPECTED_ROUTE_CLASSIFICATIONS,
        "survivor_count": EXPECTED_SURVIVOR_COUNT,
        "v2_root_method_counts": EXPECTED_ROOT_METHODS,
        "shadow_complete_counts": EXPECTED_COMPLETE_COUNTS,
        "shadow_ownership_counts": EXPECTED_OWNERSHIP_COUNTS,
        "ownership_transition_matrix": EXPECTED_TRANSITIONS,
        "shadow_proposal_count": 21,
        "SSPRK3_stage_and_endpoint_record_count": 84,
        "primary_v2_route_agreement_count": 60,
    }
    for key, value in expected_aggregate.items():
        if replay.get(key) != value:
            _stop("PREF1_COMPACT_DRIFT", f"aggregate {key}")
    occurrences = replay.get("occurrences")
    if (
        not isinstance(occurrences, list)
        or len(occurrences) != 10
        or [
            (
                item.get("retry"),
                item.get("channel"),
                item.get("shadow_ownership_class"),
                item.get("shadow_complete_classification"),
                item.get("complete_failure_cleared"),
            )
            for item in occurrences
        ]
        != list(EXPECTED_OCCURRENCE_CLASSES)
        or any(len(item.get("components", ())) != 6 for item in occurrences)
    ):
        _stop("PREF1_COMPACT_DRIFT", "occurrence summaries")
    conclusion = payload.get("conclusion")
    if conclusion != {
        "PDE_or_candidate_state_opened": False,
        "four_of_ten_complete_failures_clear_under_shadow": True,
        "independent_method_agreement": False,
        "physics_inference_permitted": False,
        "production_SSPRK3_comparator": False,
        "same_state_tableau_sensitivity_measured": True,
        "six_of_ten_complete_failures_persist_under_shadow": True,
        "state_advance_authorized": False,
        "successor_remedy_selected": False,
        "two_persistent_failures_become_endpoint_state_owned": True,
    }:
        _stop("PREF1_COMPACT_DRIFT", "conclusion")
    return result


__all__ = (
    "ARTIFACT_ID",
    "CONFIG_PATH",
    "OWNER_DOCUMENT",
    "RESULT_PATH",
    "TDG9TI2PREF1Error",
    "bind_raw_result",
    "build_pref1_result",
    "canonical_result",
    "expected_evidence",
    "validate_compact_result",
)
