"""Static SID3 authority for one guarded PRO19 generation-eight resume image.

This module consumes bytes that a standard-library bootstrap has already
captured under :mod:`proto19_execution_closure`.  It does not read Git, inspect
the mutable campaign store, acquire a writer, or advance a trajectory.  Its
job is deliberately narrower:

* bind the exact A-owned Python/input closure and the later C-owned authority
  artifacts;
* reconstruct the original GR-0 :class:`ProgressionPlan` from captured bytes;
* validate SID2's compact generation-eight certificate without reopening the
  ignored run store; and
* return immutable bytes to the runner for its own read-only live preflight.

Consequently a passing SID3 artifact never says that a trajectory is safe to
resume.  Only the runner's later store-aware preflight may make that transient
statement, and the runner must repeat it immediately before mutation.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import importlib.metadata
import json
import os
from pathlib import Path
import stat
import tomllib
from types import MappingProxyType
from typing import Any, Iterable, Mapping

import numpy as np
from numpy._core import _multiarray_umath

from .proto19_execution_closure import (
    ExecutionClosureError,
    ExecutionClosureRecord,
    parse_execution_closure_record,
)
from .proto19_progression_contract import ProgressionPlan, construct_first_event
from .proto19_sid2_binder import validate_compact_result as validate_sid2_result


ARTIFACT_ID = "FGC-1-PRO19-SID3-AUTH1"
SCHEMA_VERSION = 1
PROJECT_VERSION = "0.11.0"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
CLASSIFICATION = "premise_only_committed_gr0_resume_image_authority"

CONFIG_PATH = "configs/fgc/fgc-1-pro19-sid3-auth1.toml"
CLOSURE_PATH = "configs/fgc/fgc-1-pro19-sid3-execution-closure.json"
RESULT_PATH = "results/fgc-1-pro19-sid3-auth1.json"
STORE_PATH = "runs/fgc-2-sf1/proto17/calibration"

SID2_COMMIT = "c0ac458517025d44a020b0b0c3417ca1b844429b"
SID2_CONFIG_PATH = "configs/fgc/fgc-1-pro19-sid2-pref1.toml"
SID2_CONFIG_SHA256 = "e00ae01660ae9467c121b9c7190554dcdaf154ff292a0e7f013a614401e2c58d"
SID2_RESULT_PATH = "results/fgc-1-pro19-sid2-pref1.json"
SID2_RESULT_SHA256 = "2c70d87f6f5dd022330492562fd395004c7ec1eb9ec1d9e714ee0f0ceceb7738"

ORIGINAL_AUTHORIZATION_COMMIT = "c11ba422ce49ddd6f6175de1a9d2da675b97f2de"
ORIGINAL_PLAN_SHA256 = (
    "f3b365978aeb2ae6063742807a3161bf5d6aa0c8392cd283fdec2a79d3e46060"
)
CAMPAIGN_ID = "FGC-2-SF1-PROTO17-GR0-A3-CAL-290a65bcd6a2ba820683"

GENERATION = 8
GENERATION8_CHECKPOINT_SHA256 = (
    "6d2b6e37e7c7cfd4ce64504353d886099003f959ea2aaee1a14fc778ce14fa46"
)
GENERATION8_JOURNAL_SEQUENCE = 8
GENERATION8_JOURNAL_SHA256 = (
    "949b29e46081bf5a42c34461b3894419db86aeca6c733b7a2636c91c6cb2ee4a"
)
GENERATION8_CURSOR_SHA256 = (
    "6d594999f6aaebb5d56f0fd6ff12499e2b22a47e3daf74ac8b13b15d273fa177"
)
GENERATION8_PENDING_CAP_HEX = "0x1.aaa9612df9000p-10"
GENERATION8_ACCEPTED_TIME_HEX = "0x1.78554de5a30e0p+0"
GENERATION8_MEMBER_KEY = "RK4-2049"

PROGRESSION_EVIDENCE_PATHS = (
    "results/fgc-1-pro18-auth1.json",
    "results/fgc-1-pro19-frz1.json",
)
FACTORY_STATIC_INPUT_PATHS = (
    "configs/fgc/fgc-1-cal9-run1.toml",
    "configs/fgc/fgc-1-rsp2-run1.toml",
    "results/fgc-1-hlt10-mon10.json",
    "results/fgc-1-rsp2-frz1.json",
)
MANDATORY_STATIC_INPUT_PATHS = tuple(
    sorted(
        {
            *PROGRESSION_EVIDENCE_PATHS,
            *FACTORY_STATIC_INPUT_PATHS,
            SID2_CONFIG_PATH,
            SID2_RESULT_PATH,
        }
    )
)

AUTHORITY_MODULE = "recursive_horizons.fgc.evolution.proto19_resume_authority"
RUNNER_MODULE = "scripts.run_fgc_pro19_event1"
# Commit A owns this complete name-to-origin permission surface.  Commit C may
# bind hashes for these files, but it cannot widen the executable closure by
# adding another A-owned module and regenerating matching authority metadata.
APPROVED_IMPORT_MODULE_PATHS = MappingProxyType(
    {
        "recursive_horizons.fgc.action": "src/recursive_horizons/fgc/action.py",
        "recursive_horizons.fgc.evolution.boundary_domain": "src/recursive_horizons/fgc/evolution/boundary_domain.py",
        "recursive_horizons.fgc.evolution.cal2_source_diagnosis": "src/recursive_horizons/fgc/evolution/cal2_source_diagnosis.py",
        "recursive_horizons.fgc.evolution.cal4_common_event_diagnosis": "src/recursive_horizons/fgc/evolution/cal4_common_event_diagnosis.py",
        "recursive_horizons.fgc.evolution.calibration_runtime": "src/recursive_horizons/fgc/evolution/calibration_runtime.py",
        "recursive_horizons.fgc.evolution.floating_jet": "src/recursive_horizons/fgc/evolution/floating_jet.py",
        "recursive_horizons.fgc.evolution.gr0_calibration": "src/recursive_horizons/fgc/evolution/gr0_calibration.py",
        "recursive_horizons.fgc.evolution.gr0_direct_source": "src/recursive_horizons/fgc/evolution/gr0_direct_source.py",
        "recursive_horizons.fgc.evolution.health_monitor": "src/recursive_horizons/fgc/evolution/health_monitor.py",
        "recursive_horizons.fgc.evolution.hlt16_campaign_recovery": "src/recursive_horizons/fgc/evolution/hlt16_campaign_recovery.py",
        "recursive_horizons.fgc.evolution.hlt16_campaign_runtime": "src/recursive_horizons/fgc/evolution/hlt16_campaign_runtime.py",
        "recursive_horizons.fgc.evolution.hlt16_campaign_schema": "src/recursive_horizons/fgc/evolution/hlt16_campaign_schema.py",
        "recursive_horizons.fgc.evolution.hlt16_campaign_store": "src/recursive_horizons/fgc/evolution/hlt16_campaign_store.py",
        "recursive_horizons.fgc.evolution.hlt16_lifecycle": "src/recursive_horizons/fgc/evolution/hlt16_lifecycle.py",
        "recursive_horizons.fgc.evolution.hlt16_member_codec": "src/recursive_horizons/fgc/evolution/hlt16_member_codec.py",
        "recursive_horizons.fgc.evolution.hlt16_progression_attempt": "src/recursive_horizons/fgc/evolution/hlt16_progression_attempt.py",
        "recursive_horizons.fgc.evolution.hlt16_state_store": "src/recursive_horizons/fgc/evolution/hlt16_state_store.py",
        "recursive_horizons.fgc.evolution.numerical_engine": "src/recursive_horizons/fgc/evolution/numerical_engine.py",
        "recursive_horizons.fgc.evolution.proto10_runtime": "src/recursive_horizons/fgc/evolution/proto10_runtime.py",
        "recursive_horizons.fgc.evolution.proto11_runtime": "src/recursive_horizons/fgc/evolution/proto11_runtime.py",
        "recursive_horizons.fgc.evolution.proto12_runtime": "src/recursive_horizons/fgc/evolution/proto12_runtime.py",
        "recursive_horizons.fgc.evolution.proto14_runtime": "src/recursive_horizons/fgc/evolution/proto14_runtime.py",
        "recursive_horizons.fgc.evolution.proto15_runtime": "src/recursive_horizons/fgc/evolution/proto15_runtime.py",
        "recursive_horizons.fgc.evolution.proto17_pure_construction": "src/recursive_horizons/fgc/evolution/proto17_pure_construction.py",
        "recursive_horizons.fgc.evolution.proto19_cfl_continuation_authority": "src/recursive_horizons/fgc/evolution/proto19_cfl_continuation_authority.py",
        "recursive_horizons.fgc.evolution.proto19_execution_closure": "src/recursive_horizons/fgc/evolution/proto19_execution_closure.py",
        "recursive_horizons.fgc.evolution.proto19_gr0_static_factory": "src/recursive_horizons/fgc/evolution/proto19_gr0_static_factory.py",
        "recursive_horizons.fgc.evolution.proto19_launch_authority": "src/recursive_horizons/fgc/evolution/proto19_launch_authority.py",
        "recursive_horizons.fgc.evolution.proto19_progression_contract": "src/recursive_horizons/fgc/evolution/proto19_progression_contract.py",
        "recursive_horizons.fgc.evolution.proto19_resume_authority": "src/recursive_horizons/fgc/evolution/proto19_resume_authority.py",
        "recursive_horizons.fgc.evolution.proto19_sid1_authority": "src/recursive_horizons/fgc/evolution/proto19_sid1_authority.py",
        "recursive_horizons.fgc.evolution.proto19_sid2_binder": "src/recursive_horizons/fgc/evolution/proto19_sid2_binder.py",
        "recursive_horizons.fgc.evolution.proto4_admission": "src/recursive_horizons/fgc/evolution/proto4_admission.py",
        "recursive_horizons.fgc.evolution.proto5_runtime": "src/recursive_horizons/fgc/evolution/proto5_runtime.py",
        "recursive_horizons.fgc.evolution.proto6_runtime": "src/recursive_horizons/fgc/evolution/proto6_runtime.py",
        "recursive_horizons.fgc.evolution.proto7_runtime": "src/recursive_horizons/fgc/evolution/proto7_runtime.py",
        "recursive_horizons.fgc.evolution.proto8_runtime": "src/recursive_horizons/fgc/evolution/proto8_runtime.py",
        "recursive_horizons.fgc.evolution.proto9_runtime": "src/recursive_horizons/fgc/evolution/proto9_runtime.py",
        "recursive_horizons.fgc.evolution.protocol_v17": "src/recursive_horizons/fgc/evolution/protocol_v17.py",
        "recursive_horizons.fgc.evolution.protocol_v8": "src/recursive_horizons/fgc/evolution/protocol_v8.py",
        "recursive_horizons.fgc.evolution.spectral_sensitivity": "src/recursive_horizons/fgc/evolution/spectral_sensitivity.py",
        "recursive_horizons.fgc.evolution.spectral_sensitivity_v12": "src/recursive_horizons/fgc/evolution/spectral_sensitivity_v12.py",
        "recursive_horizons.fgc.evolution.src3_reference_balanced_source": "src/recursive_horizons/fgc/evolution/src3_reference_balanced_source.py",
        "recursive_horizons.fgc.evolution.src4_vectorized_reference_source": "src/recursive_horizons/fgc/evolution/src4_vectorized_reference_source.py",
        "recursive_horizons.fgc.evolution.tdg5_stage_complete_refinement_runtime": "src/recursive_horizons/fgc/evolution/tdg5_stage_complete_refinement_runtime.py",
        "recursive_horizons.fgc.evolution.tdg6_temporal_admission_design": "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_design.py",
        "recursive_horizons.fgc.evolution.tdg6_temporal_admission_runtime": "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_runtime.py",
        "recursive_horizons.fgc.evolution.tdg7_binary64_subdivision_lattice": "src/recursive_horizons/fgc/evolution/tdg7_binary64_subdivision_lattice.py",
        "recursive_horizons.fgc.evolution.tdg7_stage_safe_runtime": "src/recursive_horizons/fgc/evolution/tdg7_stage_safe_runtime.py",
        "recursive_horizons.fgc.evolution.vectorized_source": "src/recursive_horizons/fgc/evolution/vectorized_source.py",
        "recursive_horizons.fgc.exact_interval": "src/recursive_horizons/fgc/exact_interval.py",
        "recursive_horizons.fgc.exact_linear_algebra": "src/recursive_horizons/fgc/exact_linear_algebra.py",
        "recursive_horizons.fgc.exact_tangent": "src/recursive_horizons/fgc/exact_tangent.py",
        "recursive_horizons.fgc.initial_data_family": "src/recursive_horizons/fgc/initial_data_family.py",
        "recursive_horizons.fgc.initial_data_preflight": "src/recursive_horizons/fgc/initial_data_preflight.py",
        "recursive_horizons.fgc.interval_tangent": "src/recursive_horizons/fgc/interval_tangent.py",
        "recursive_horizons.fgc.modified_harmonic": "src/recursive_horizons/fgc/modified_harmonic.py",
        "recursive_horizons.fgc.modified_harmonic_constraints": "src/recursive_horizons/fgc/modified_harmonic_constraints.py",
        "recursive_horizons.fgc.modified_harmonic_first_order": "src/recursive_horizons/fgc/modified_harmonic_first_order.py",
        "recursive_horizons.fgc.modified_harmonic_implicit": "src/recursive_horizons/fgc/modified_harmonic_implicit.py",
        "recursive_horizons.fgc.modified_harmonic_reference": "src/recursive_horizons/fgc/modified_harmonic_reference.py",
        "recursive_horizons.fgc.reference_connection": "src/recursive_horizons/fgc/reference_connection.py",
        "recursive_horizons.fgc.scoped_run_authorization": "src/recursive_horizons/fgc/scoped_run_authorization.py",
        "recursive_horizons.fgc.spherical_reduction": "src/recursive_horizons/fgc/spherical_reduction.py",
        "recursive_horizons.fgc.spherical_symbol": "src/recursive_horizons/fgc/spherical_symbol.py",
        "scripts.run_fgc_pro19_event1": "scripts/run_fgc_pro19_event1.py",
    }
)
MANDATORY_IMPORT_MODULES = frozenset(
    {
        AUTHORITY_MODULE,
        RUNNER_MODULE,
        "recursive_horizons.fgc.evolution.proto19_execution_closure",
        "recursive_horizons.fgc.evolution.proto19_gr0_static_factory",
        "recursive_horizons.fgc.evolution.proto19_progression_contract",
        "recursive_horizons.fgc.evolution.proto19_sid2_binder",
    }
)
MANDATORY_NAMESPACE_PATHS = {
    "recursive_horizons": "src/recursive_horizons",
    "recursive_horizons.fgc": "src/recursive_horizons/fgc",
    "recursive_horizons.fgc.evolution": "src/recursive_horizons/fgc/evolution",
    "scripts": "scripts",
}
MANDATORY_FORBIDDEN_MODULE_PREFIXES = frozenset(
    {
        "recursive_horizons.fgc.evolution.covariant_principal_health",
        "recursive_horizons.fgc.evolution.initial_state_bridge",
        "recursive_horizons.fgc.evolution.nonlinear_source",
        "recursive_horizons.fgc.evolution.proto19_progression_inputs",
        "recursive_horizons.fgc.evolution.static_initial_admission",
        "scripts.run_fgc_gr0_calibration",
    }
)
MANDATORY_FORBIDDEN_SCRIPT_PREFIXES = frozenset(
    {
        "scripts/run_fgc_gr0_calibration.py",
        "src/recursive_horizons/fgc/evolution/proto19_progression_inputs.py",
        "src/recursive_horizons/fgc/evolution/static_initial_admission.py",
    }
)

EXPECTED_NONCLAIMS = [
    "SID3 binds a committed GR-0 execution image and immutable inputs; it does not inspect or authorize the live mutable store by artifact alone.",
    "A passing SID3 artifact leaves safe_to_resume false until a separate guarded read-only preflight authenticates generation eight and every runtime premise.",
    "SID3 does not complete a common event or GR-0 calibration and does not authorize SGB-L, FGC-QR, DEF1, retained-EFT, transition, or physical claims.",
    "The committed GR-0 image loads shared analytic modules that also contain dormant candidate-capable definitions and may instantiate their default parameter records; SID3 proves that no candidate-only module, runtime configuration, initial/evolution state, branch selection, or outcome is opened, not that every loaded source file contains only GR-0 definitions.",
]
EXPECTED_SCOPE = {
    "raw_store_mutation": False,
    "writer_lease_acquisition": False,
    "PDE_proposal_execution": False,
    "GEN0_reimport": False,
    "trajectory_resume": False,
    "candidate_branch_opened": False,
    "shared_candidate_capable_definitions_loaded": True,
    "read_only_runtime_preflight_required": True,
    "runtime_recheck_before_mutation_required": True,
}
EXPECTED_CLAIMS = {
    "committed_resume_contract_frozen": True,
    "committed_resume_authority_present": True,
    "committed_resume_image_authenticated_by_artifact_alone": False,
    "exact_generation8_inputs_bound": True,
    "safe_to_resume_trajectory": False,
    "trajectory_resume_authorized_by_artifact_alone": False,
    "GR0_calibration_completed": False,
    "candidate_execution_authorized": False,
    "candidate_runtime_configuration_state_or_outcome_opened": False,
    "physical_result_earned": False,
}


class ResumeAuthorityError(RuntimeError):
    """The prospective SID3 contract or captured committed image differs."""


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _require_sha(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or value.lower() != value:
        raise ResumeAuthorityError(f"{label} is not SHA-256")
    try:
        int(value, 16)
    except ValueError as exc:
        raise ResumeAuthorityError(f"{label} is not SHA-256") from exc
    return value


def _require_oid(value: object, length: int, label: str) -> str:
    if not isinstance(value, str) or len(value) != length or value.lower() != value:
        raise ResumeAuthorityError(f"{label} is not a Git object id")
    try:
        int(value, 16)
    except ValueError as exc:
        raise ResumeAuthorityError(f"{label} is not a Git object id") from exc
    return value


def _canonical_json(value: object) -> bytes:
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
        raise ResumeAuthorityError("SID3 JSON is not canonicalizable") from exc


def canonical_result(value: object) -> bytes:
    """Return the sole accepted tracked JSON representation."""

    return _canonical_json(value)


def _reject_duplicates(items: Iterable[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in items:
        if key in answer:
            raise ValueError(key)
        answer[key] = value
    return answer


def _json(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(raw.decode("ascii"), object_pairs_hook=_reject_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise ResumeAuthorityError(f"{label} is malformed") from exc
    if not isinstance(value, dict) or _canonical_json(value) != raw:
        raise ResumeAuthorityError(f"{label} is noncanonical")
    return value


def _absolute_regular_identity(path: Path, label: str) -> tuple[str, str]:
    real = Path(os.path.realpath(path))
    try:
        before = real.lstat()
    except OSError as exc:
        raise ResumeAuthorityError(f"{label} is unavailable") from exc
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        raise ResumeAuthorityError(f"{label} is not a regular file")
    descriptor = -1
    try:
        descriptor = os.open(real, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
        active = os.fstat(descriptor)
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
            raise ResumeAuthorityError(f"{label} raced before read")
        chunks: list[bytes] = []
        while True:
            chunk = os.read(descriptor, 1 << 20)
            if not chunk:
                break
            chunks.append(chunk)
        after = os.fstat(descriptor)
        if identity != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
            after.st_ctime_ns,
        ):
            raise ResumeAuthorityError(f"{label} changed during read")
    except OSError as exc:
        raise ResumeAuthorityError(f"{label} cannot be read safely") from exc
    finally:
        if descriptor != -1:
            os.close(descriptor)
    return str(real), _sha(b"".join(chunks))


@dataclass(frozen=True, slots=True)
class NumpyRuntimeIdentity:
    version: str
    package_origin_realpath: str
    package_origin_sha256: str
    core_extension_origin_realpath: str
    core_extension_sha256: str
    distribution_metadata_realpath: str
    distribution_metadata_sha256: str

    def mapping(self) -> dict[str, str]:
        return {
            "version": self.version,
            "package_origin_realpath": self.package_origin_realpath,
            "package_origin_sha256": self.package_origin_sha256,
            "core_extension_origin_realpath": self.core_extension_origin_realpath,
            "core_extension_sha256": self.core_extension_sha256,
            "distribution_metadata_realpath": self.distribution_metadata_realpath,
            "distribution_metadata_sha256": self.distribution_metadata_sha256,
        }


def observe_numpy_runtime_identity() -> NumpyRuntimeIdentity:
    """Bind NumPy's Python entry point, native core, and package metadata."""

    package_file = getattr(np, "__file__", None)
    core_file = getattr(_multiarray_umath, "__file__", None)
    if not isinstance(package_file, str) or not isinstance(core_file, str):
        raise ResumeAuthorityError("NumPy package/native origin is absent")
    package_path, package_hash = _absolute_regular_identity(
        Path(package_file), "NumPy package origin"
    )
    core_path, core_hash = _absolute_regular_identity(
        Path(core_file), "NumPy core extension"
    )
    try:
        distribution = importlib.metadata.distribution("numpy")
        matches = [
            distribution.locate_file(item)
            for item in (distribution.files or ())
            if str(item).endswith(".dist-info/METADATA")
        ]
    except importlib.metadata.PackageNotFoundError as exc:
        raise ResumeAuthorityError("NumPy distribution metadata is absent") from exc
    if len(matches) != 1:
        raise ResumeAuthorityError("NumPy distribution METADATA is ambiguous")
    metadata_path, metadata_hash = _absolute_regular_identity(
        Path(matches[0]), "NumPy distribution METADATA"
    )
    if not isinstance(np.__version__, str) or not np.__version__:
        raise ResumeAuthorityError("NumPy version is absent")
    return NumpyRuntimeIdentity(
        version=np.__version__,
        package_origin_realpath=package_path,
        package_origin_sha256=package_hash,
        core_extension_origin_realpath=core_path,
        core_extension_sha256=core_hash,
        distribution_metadata_realpath=metadata_path,
        distribution_metadata_sha256=metadata_hash,
    )


def _numpy_identity(value: object) -> NumpyRuntimeIdentity:
    fields = {
        "version",
        "package_origin_realpath",
        "package_origin_sha256",
        "core_extension_origin_realpath",
        "core_extension_sha256",
        "distribution_metadata_realpath",
        "distribution_metadata_sha256",
    }
    if not isinstance(value, Mapping) or set(value) != fields:
        raise ResumeAuthorityError("NumPy runtime identity fields differ")
    for field in fields:
        if not isinstance(value[field], str) or not value[field]:
            raise ResumeAuthorityError("NumPy runtime identity value differs")
    for field in fields:
        if field.endswith("_sha256"):
            _require_sha(value[field], f"NumPy {field}")
        elif field.endswith("_realpath") and (
            not os.path.isabs(value[field])
            or os.path.realpath(value[field]) != value[field]
        ):
            raise ResumeAuthorityError(f"NumPy {field} is not an absolute realpath")
    return NumpyRuntimeIdentity(**dict(value))


def _predecessor() -> dict[str, object]:
    return {
        "artifact_id": "FGC-1-PRO19-SID2-PREF1",
        "sealed_commit": SID2_COMMIT,
        "config_path": SID2_CONFIG_PATH,
        "config_sha256": SID2_CONFIG_SHA256,
        "result_path": SID2_RESULT_PATH,
        "result_sha256": SID2_RESULT_SHA256,
    }


def _progression() -> dict[str, object]:
    return {
        "original_authorization_commit": ORIGINAL_AUTHORIZATION_COMMIT,
        "original_plan_sha256": ORIGINAL_PLAN_SHA256,
        "campaign_id": CAMPAIGN_ID,
        "branch": "GR-0",
        "amplitude": "3",
        "event": 23,
        "start_rational": "23/16",
        "target_rational": "3/2",
        "candidate_branch_opened": False,
    }


def _recovered_generation() -> dict[str, object]:
    return {
        "generation": GENERATION,
        "checkpoint_sha256": GENERATION8_CHECKPOINT_SHA256,
        "journal_sequence": GENERATION8_JOURNAL_SEQUENCE,
        "journal_sha256": GENERATION8_JOURNAL_SHA256,
        "member_key": GENERATION8_MEMBER_KEY,
        "member_accepted_time_hex": GENERATION8_ACCEPTED_TIME_HEX,
        "member_cursor_sha256": GENERATION8_CURSOR_SHA256,
        "member_mode": "RETRY_PENDING",
        "pending_owner": "temporal",
        "pending_cap_hex": GENERATION8_PENDING_CAP_HEX,
        "disposition": "nonterminal",
        "candidate_branch_opened": False,
        "physical_state_advanced_by_recovery": False,
    }


def _string_list(value: object, label: str) -> list[str]:
    if (
        not isinstance(value, list)
        or any(not isinstance(item, str) or not item for item in value)
        or value != sorted(set(value))
    ):
        raise ResumeAuthorityError(f"{label} is not canonical")
    return list(value)


def _closure_mapping(raw: bytes) -> tuple[dict[str, Any], ExecutionClosureRecord]:
    value = _json(raw, "SID3 execution closure")
    try:
        record = parse_execution_closure_record(value)
    except ExecutionClosureError as exc:
        raise ResumeAuthorityError("SID3 execution closure differs") from exc
    return value, record


def _execution_contract(
    value: object,
    *,
    closure_raw: bytes,
    closure: ExecutionClosureRecord,
) -> dict[str, Any]:
    fields = {
        "path",
        "raw_sha256",
        "canonical_sha256",
        "implementation_commit",
        "authority_module",
        "runner_module",
        "required_import_modules",
        "required_namespaces",
        "required_forbidden_module_prefixes",
        "required_forbidden_script_prefixes",
        "authority_delta_paths",
        "static_input_paths",
    }
    if not isinstance(value, Mapping) or set(value) != fields:
        raise ResumeAuthorityError("SID3 execution contract fields differ")
    imports = _string_list(value["required_import_modules"], "required imports")
    namespaces = _string_list(value["required_namespaces"], "required namespaces")
    forbidden_modules = _string_list(
        value["required_forbidden_module_prefixes"], "forbidden modules"
    )
    forbidden_scripts = _string_list(
        value["required_forbidden_script_prefixes"], "forbidden scripts"
    )
    delta = _string_list(value["authority_delta_paths"], "authority delta")
    inputs = _string_list(value["static_input_paths"], "static inputs")
    file_inputs = sorted(item.path for item in closure.files if item.kind == "input")
    expected = {
        "path": CLOSURE_PATH,
        "raw_sha256": _sha(closure_raw),
        "canonical_sha256": closure.canonical_sha256,
        "implementation_commit": closure.implementation_commit,
        "authority_module": AUTHORITY_MODULE,
        "runner_module": RUNNER_MODULE,
        "required_import_modules": sorted(item.module for item in closure.imports),
        "required_namespaces": sorted(item.module for item in closure.namespaces),
        "required_forbidden_module_prefixes": list(closure.forbidden_module_prefixes),
        "required_forbidden_script_prefixes": list(closure.forbidden_script_prefixes),
        "authority_delta_paths": list(closure.authority_delta_paths),
        "static_input_paths": file_inputs,
    }
    normalized = {
        **dict(value),
        "required_import_modules": imports,
        "required_namespaces": namespaces,
        "required_forbidden_module_prefixes": forbidden_modules,
        "required_forbidden_script_prefixes": forbidden_scripts,
        "authority_delta_paths": delta,
        "static_input_paths": inputs,
    }
    if normalized != expected:
        raise ResumeAuthorityError("SID3 execution contract/closure differs")
    observed_import_paths = {item.module: item.path for item in closure.imports}
    if observed_import_paths != APPROVED_IMPORT_MODULE_PATHS:
        raise ResumeAuthorityError("SID3 import permission surface differs")
    if not MANDATORY_IMPORT_MODULES.issubset(imports):
        raise ResumeAuthorityError("SID3 closure omits a mandatory runtime module")
    if namespaces != sorted(MANDATORY_NAMESPACE_PATHS):
        raise ResumeAuthorityError("SID3 closure namespace set differs")
    observed_namespaces = {item.module: item.path for item in closure.namespaces}
    if observed_namespaces != MANDATORY_NAMESPACE_PATHS:
        raise ResumeAuthorityError("SID3 closure namespace path differs")
    if not MANDATORY_FORBIDDEN_MODULE_PREFIXES.issubset(forbidden_modules):
        raise ResumeAuthorityError("SID3 closure omits a forbidden legacy module")
    if not MANDATORY_FORBIDDEN_SCRIPT_PREFIXES.issubset(forbidden_scripts):
        raise ResumeAuthorityError("SID3 closure omits a forbidden legacy script")
    if not set(MANDATORY_STATIC_INPUT_PATHS).issubset(inputs):
        raise ResumeAuthorityError("SID3 closure omits a mandatory static input")
    if not {CONFIG_PATH, CLOSURE_PATH, RESULT_PATH}.issubset(delta):
        raise ResumeAuthorityError("SID3 authority delta omits its own artifacts")
    disallowed_tokens = ("static_initial_admission", "run_fgc_gr0_calibration")
    if any(
        any(token in item.module for token in disallowed_tokens)
        for item in closure.imports
    ):
        raise ResumeAuthorityError("SID3 closure imports a forbidden legacy path")
    return normalized


def _parse_config(
    config_raw: bytes,
    closure_raw: bytes,
) -> tuple[dict[str, Any], ExecutionClosureRecord]:
    try:
        value = tomllib.loads(config_raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise ResumeAuthorityError("SID3 config is malformed") from exc
    fields = {
        "schema_version",
        "artifact_id",
        "project_version",
        "target_protocol",
        "classification",
        "nonclaims",
        "predecessor",
        "execution_closure",
        "progression",
        "recovered_generation",
        "numpy_runtime",
        "scope",
        "claims",
    }
    if not isinstance(value, dict) or set(value) != fields:
        raise ResumeAuthorityError("SID3 config fields differ")
    closure_value, closure = _closure_mapping(closure_raw)
    if (
        value["schema_version"] != SCHEMA_VERSION
        or value["artifact_id"] != ARTIFACT_ID
        or value["project_version"] != PROJECT_VERSION
        or value["target_protocol"] != TARGET_PROTOCOL
        or value["classification"] != CLASSIFICATION
        or value["nonclaims"] != EXPECTED_NONCLAIMS
        or value["predecessor"] != _predecessor()
        or value["progression"] != _progression()
        or value["recovered_generation"] != _recovered_generation()
        or value["scope"] != EXPECTED_SCOPE
        or value["claims"] != EXPECTED_CLAIMS
    ):
        raise ResumeAuthorityError("SID3 exact contract differs")
    _execution_contract(
        value["execution_closure"], closure_raw=closure_raw, closure=closure
    )
    configured_numpy = _numpy_identity(value["numpy_runtime"])
    observed_numpy = observe_numpy_runtime_identity()
    if configured_numpy != observed_numpy:
        raise ResumeAuthorityError("SID3 NumPy runtime identity differs")
    # Keep the duplicate-rejecting parsed value alive through result building.
    if closure_value != closure.to_mapping():
        raise ResumeAuthorityError("SID3 execution closure normalization differs")
    return value, closure


def _result_from_config(
    config_raw: bytes,
    closure_raw: bytes,
    config: Mapping[str, Any],
    closure: ExecutionClosureRecord,
) -> dict[str, Any]:
    inputs = {item.path: item.sha256 for item in closure.files if item.kind == "input"}
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "target_protocol": TARGET_PROTOCOL,
        "classification": CLASSIFICATION,
        "gate_status": "pass",
        "source_config_sha256": _sha(config_raw),
        "execution_closure_raw_sha256": _sha(closure_raw),
        "execution_closure_canonical_sha256": closure.canonical_sha256,
        "artifact_payload": {
            "predecessor": dict(config["predecessor"]),
            "execution_image": {
                "implementation_commit": closure.implementation_commit,
                "file_pin_count": len(closure.files),
                "import_modules": [item.module for item in closure.imports],
                "namespaces": [item.module for item in closure.namespaces],
                "forbidden_module_prefixes": list(closure.forbidden_module_prefixes),
                "forbidden_script_prefixes": list(closure.forbidden_script_prefixes),
                "authority_delta_paths": list(closure.authority_delta_paths),
                "static_input_sha256": inputs,
                "numpy_runtime": dict(config["numpy_runtime"]),
            },
            "progression": dict(config["progression"]),
            "recovered_generation": dict(config["recovered_generation"]),
            "scope": dict(config["scope"]),
            "claims": dict(config["claims"]),
            "nonclaims": list(config["nonclaims"]),
        },
    }


def build_sid3_result(config_raw: bytes, closure_raw: bytes) -> dict[str, Any]:
    """Build the deterministic artifact; no store or execution is opened."""

    config, closure = _parse_config(config_raw, closure_raw)
    return _result_from_config(config_raw, closure_raw, config, closure)


def validate_sid3_result(
    config_raw: bytes,
    closure_raw: bytes,
    result_raw: bytes,
) -> dict[str, Any]:
    """Validate exact tracked SID3 bytes without granting runtime permission."""

    config, closure = _parse_config(config_raw, closure_raw)
    value = _json(result_raw, "SID3 result")
    expected = _result_from_config(config_raw, closure_raw, config, closure)
    if value != expected:
        raise ResumeAuthorityError("SID3 result/config/closure binding differs")
    return value


@dataclass(frozen=True, slots=True)
class ResumeAuthorityReceipt:
    """Immutable committed image handed to the store-aware runner."""

    authority_commit: str
    implementation_commit: str
    closure_sha256: str
    closure_canonical_sha256: str
    config_sha256: str
    result_sha256: str
    sid2_commit: str
    sid2_config_sha256: str
    sid2_result_sha256: str
    original_authorization_commit: str
    original_plan_sha256: str
    campaign_id: str
    store_path: str
    checkpoint_generation: int
    checkpoint_sha256: str
    journal_sequence: int
    journal_tip_sha256: str
    recovery_member_key: str
    recovery_cursor_sha256: str
    recovery_pending_cap_hex: str
    numpy_runtime: NumpyRuntimeIdentity
    progression_plan: ProgressionPlan
    static_input_bytes: Mapping[str, bytes]

    def __post_init__(self) -> None:
        copied = {
            key: bytes(value) for key, value in sorted(self.static_input_bytes.items())
        }
        object.__setattr__(self, "static_input_bytes", MappingProxyType(copied))

    @property
    def factory_static_input_bytes(self) -> Mapping[str, bytes]:
        """The exact four-byte inventory accepted by the GR-0 shell factory."""

        return MappingProxyType(
            {path: self.static_input_bytes[path] for path in FACTORY_STATIC_INPUT_PATHS}
        )

    @property
    def progression_evidence_bytes(self) -> Mapping[str, bytes]:
        """The two compact predecessors used to rebuild the original plan."""

        return MappingProxyType(
            {path: self.static_input_bytes[path] for path in PROGRESSION_EVIDENCE_PATHS}
        )


def authorize_resume_image(
    config_raw: bytes,
    closure_raw: bytes,
    result_raw: bytes,
    *,
    captured_files: Mapping[str, bytes],
    authority_commit: str,
) -> ResumeAuthorityReceipt:
    """Bind one exact committed image, still without inspecting live state.

    ``captured_files`` must contain every file pin returned by the execution
    closure verifier, not a caller-selected subset.  The returned
    ``static_input_bytes`` are therefore a race-free subset of that complete
    image and can be passed directly to the GR-0 factory.
    """

    config, closure = _parse_config(config_raw, closure_raw)
    result = validate_sid3_result(config_raw, closure_raw, result_raw)
    oid_length = 40 if closure.git_object_format == "sha1" else 64
    commit = _require_oid(authority_commit, oid_length, "authority commit")
    if not isinstance(captured_files, Mapping) or any(
        not isinstance(path, str) or not isinstance(raw, bytes)
        for path, raw in captured_files.items()
    ):
        raise ResumeAuthorityError("captured closure bytes differ")
    file_map = {item.path: item for item in closure.files}
    if set(captured_files) != set(file_map):
        raise ResumeAuthorityError("captured closure inventory differs")
    for path, pin in file_map.items():
        if _sha(captured_files[path]) != pin.sha256:
            raise ResumeAuthorityError(f"captured closure byte differs: {path}")
    static_inputs = {
        path: captured_files[path]
        for path, pin in file_map.items()
        if pin.kind == "input"
    }
    expected_inputs = set(config["execution_closure"]["static_input_paths"])
    if set(static_inputs) != expected_inputs:
        raise ResumeAuthorityError("captured static-input inventory differs")

    if (
        _sha(static_inputs[SID2_CONFIG_PATH]) != SID2_CONFIG_SHA256
        or _sha(static_inputs[SID2_RESULT_PATH]) != SID2_RESULT_SHA256
    ):
        raise ResumeAuthorityError("captured SID2 compact bytes differ")
    try:
        sid2 = validate_sid2_result(
            static_inputs[SID2_CONFIG_PATH], static_inputs[SID2_RESULT_PATH]
        )
    except Exception as exc:
        raise ResumeAuthorityError("captured SID2 certificate is invalid") from exc
    anchor = sid2.get("artifact_payload", {}).get("live_anchor", {})
    generation = anchor.get("generation8", {}) if isinstance(anchor, Mapping) else {}
    retry = anchor.get("retry_state", {}) if isinstance(anchor, Mapping) else {}
    if (
        generation.get("generation") != GENERATION
        or generation.get("checkpoint_sha256") != GENERATION8_CHECKPOINT_SHA256
        or generation.get("journal_sequence") != GENERATION8_JOURNAL_SEQUENCE
        or generation.get("journal_sha256") != GENERATION8_JOURNAL_SHA256
        or retry.get("member_key") != GENERATION8_MEMBER_KEY
        or retry.get("cursor_sha256") != GENERATION8_CURSOR_SHA256
        or retry.get("mode") != "RETRY_PENDING"
        or retry.get("pending_cap_hex") != GENERATION8_PENDING_CAP_HEX
    ):
        raise ResumeAuthorityError("SID2 generation-eight anchor differs")

    evidence = {path: static_inputs[path] for path in PROGRESSION_EVIDENCE_PATHS}
    try:
        plan = construct_first_event(
            Path("."),
            authorization_commit=ORIGINAL_AUTHORIZATION_COMMIT,
            evidence_bytes=evidence,
        )
    except Exception as exc:
        raise ResumeAuthorityError(
            "original progression plan cannot be rebuilt"
        ) from exc
    if (
        plan.sha256 != ORIGINAL_PLAN_SHA256
        or plan.campaign_id != CAMPAIGN_ID
        or plan.branch != "GR-0"
        or plan.amplitude != "3"
        or plan.common_event_index != 23
        or plan.target_time.rational != "3/2"
    ):
        raise ResumeAuthorityError("original progression plan differs")
    numpy_runtime = _numpy_identity(config["numpy_runtime"])
    if result["artifact_payload"]["claims"]["safe_to_resume_trajectory"] is not False:
        raise ResumeAuthorityError("SID3 artifact attempted runtime promotion")
    return ResumeAuthorityReceipt(
        authority_commit=commit,
        implementation_commit=closure.implementation_commit,
        closure_sha256=_sha(closure_raw),
        closure_canonical_sha256=closure.canonical_sha256,
        config_sha256=_sha(config_raw),
        result_sha256=_sha(result_raw),
        sid2_commit=SID2_COMMIT,
        sid2_config_sha256=SID2_CONFIG_SHA256,
        sid2_result_sha256=SID2_RESULT_SHA256,
        original_authorization_commit=ORIGINAL_AUTHORIZATION_COMMIT,
        original_plan_sha256=ORIGINAL_PLAN_SHA256,
        campaign_id=CAMPAIGN_ID,
        store_path=STORE_PATH,
        checkpoint_generation=GENERATION,
        checkpoint_sha256=GENERATION8_CHECKPOINT_SHA256,
        journal_sequence=GENERATION8_JOURNAL_SEQUENCE,
        journal_tip_sha256=GENERATION8_JOURNAL_SHA256,
        recovery_member_key=GENERATION8_MEMBER_KEY,
        recovery_cursor_sha256=GENERATION8_CURSOR_SHA256,
        recovery_pending_cap_hex=GENERATION8_PENDING_CAP_HEX,
        numpy_runtime=numpy_runtime,
        progression_plan=plan,
        static_input_bytes=static_inputs,
    )


__all__ = [
    "APPROVED_IMPORT_MODULE_PATHS",
    "ARTIFACT_ID",
    "AUTHORITY_MODULE",
    "CLASSIFICATION",
    "CLOSURE_PATH",
    "CONFIG_PATH",
    "FACTORY_STATIC_INPUT_PATHS",
    "MANDATORY_FORBIDDEN_MODULE_PREFIXES",
    "MANDATORY_FORBIDDEN_SCRIPT_PREFIXES",
    "MANDATORY_IMPORT_MODULES",
    "MANDATORY_STATIC_INPUT_PATHS",
    "NumpyRuntimeIdentity",
    "RESULT_PATH",
    "RUNNER_MODULE",
    "STORE_PATH",
    "ResumeAuthorityError",
    "ResumeAuthorityReceipt",
    "authorize_resume_image",
    "build_sid3_result",
    "canonical_result",
    "observe_numpy_runtime_identity",
    "validate_sid3_result",
]
