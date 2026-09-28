#!/usr/bin/env python3
"""Replay and capture the immutable PROTO11 GR-0 affine source wall.

This runner consumes CAL8/PREF10 and the original raw PROTO11 campaign.  It
replays only the implicated RK4-8193 member, requires the first evolved event
and all 164 rejection records to match the immutable ledger exactly, captures
the last candidate-endpoint source state, and compares stable arithmetic
routes under the unchanged complete REF1 residual gate.

The output is raw ignored run evidence.  It authorizes no successor protocol,
candidate trajectory, mechanism claim, or retained-EFT claim.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite
import os
from pathlib import Path
import subprocess
import sys
import time
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_cal8_pref10 as cal8  # noqa: E402
from scripts import run_fgc_gr0_calibration as inherited  # noqa: E402
from scripts import run_fgc_gr0_calibration_v7 as proto7_runner  # noqa: E402
from scripts import run_fgc_gr0_calibration_v11 as proto11_runner  # noqa: E402
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto11_runtime import (  # noqa: E402
    Proto11GR0EvolutionOperator,
    reference_balanced_spatial_derivatives,
)
from recursive_horizons.fgc.evolution.src2_affine_arithmetic import (  # noqa: E402
    AffineArithmeticComparison,
    ArithmeticCandidate,
    CapturedAffineSystem,
    compare_affine_arithmetic,
)


Q = Fraction
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-src2-cap1.toml"
RUNNER_ID = "FGC-1-SRC2-CAP1-RUNNER"


def _serial(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return _serial(value.tolist())
    if isinstance(value, np.generic):
        return _serial(value.item())
    if isinstance(value, (tuple, list)):
        return [_serial(item) for item in value]
    if hasattr(value, "__dataclass_fields__"):
        return _serial(asdict(value))
    if isinstance(value, float):
        if not isfinite(value):
            raise ValueError("nonfinite value cannot enter SRC2 evidence")
        return value
    if value is None or isinstance(value, (str, bool, int)):
        return value
    raise TypeError(f"unsupported SRC2 evidence type {type(value).__name__}")


def _canonical(value: Any) -> bytes:
    return (
        json.dumps(
            _serial(value),
            sort_keys=True,
            indent=2,
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _canonical_line(value: Any) -> bytes:
    return (
        json.dumps(
            _serial(value),
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _git_head() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPOSITORY,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _write_exclusive(path: Path, payload: bytes) -> None:
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        try:
            path.unlink(missing_ok=True)
        finally:
            raise


def _write_npz_exclusive(
    path: Path,
    *,
    metadata: Mapping[str, Any],
    arrays: Mapping[str, np.ndarray],
) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            np.savez_compressed(
                handle,
                metadata_utf8=np.frombuffer(_canonical(metadata), dtype=np.uint8),
                **arrays,
            )
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _load_canonical_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict) or path.read_bytes() != _canonical(value):
        raise ValueError(f"{path} is not canonical JSON")
    return value


def _load_canonical_jsonl(path: Path) -> list[dict[str, Any]]:
    payload = path.read_bytes()
    if not payload or not payload.endswith(b"\n"):
        raise ValueError("PROTO11 event ledger is absent or incomplete")
    records: list[dict[str, Any]] = []
    offset = 0
    for line in payload.splitlines(keepends=True):
        value = json.loads(line)
        if not isinstance(value, dict) or line != _canonical_line(value):
            raise ValueError(f"PROTO11 event ledger line at byte {offset} is noncanonical")
        records.append(value)
        offset += len(line)
    return records


def _bits_equal(left: float, right: float) -> bool:
    return np.float64(left).tobytes() == np.float64(right).tobytes()


def _required_tables(config: Mapping[str, Any]) -> None:
    if set(config) != {
        "schema_version",
        "artifact_id",
        "project_version",
        "metric_signature",
        "riemann_convention",
        "scope",
        "lineage",
        "replay",
        "arithmetic",
        "output",
        "claims",
    }:
        raise ValueError("SRC2 capture config top-level contract differs")
    if (
        config["schema_version"] != 1
        or config["artifact_id"] != "FGC-1-SRC2-CAP1"
        or config["metric_signature"] != "-+++"
        or config["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("SRC2 capture config identity differs")
    scope = config["scope"]
    if (
        scope.get("branch") != "GR-0"
        or scope.get("amplitude") != "5/2"
        or scope.get("method") != "RK4"
        or scope.get("point_count") != 8193
        or scope.get("SGBL_trajectory_read") is not False
        or scope.get("FGCQR_trajectory_read") is not False
        or scope.get("mechanism_question_answered") is not False
    ):
        raise ValueError("SRC2 capture scope is absent or over-broad")
    claims = config["claims"]
    if not claims or any(value is not False for value in claims.values()):
        raise ValueError("SRC2 capture config prematurely promotes a claim")
    arithmetic = config["arithmetic"]
    if (
        Q(arithmetic.get("raw_complete_residual_maximum", "0")) != Q(1, 10**12)
        or Q(arithmetic.get("kinetic_condition_number_maximum", "0")) != Q(10**10)
        or arithmetic.get("maximum_complete_refinement_iterations") != 16
        or arithmetic.get("symmetric_power_two_seed_scales") != ["1", "16", "256"]
        or arithmetic.get("exact_binary_point_limit") != 8
        or Q(arithmetic.get("exact_binary_selection_ratio", "0")) != Q(1, 2)
        or arithmetic.get("complete_unredefined_binary64_residual_is_the_only_gate")
        is not True
        or arithmetic.get("reconstructed_affine_residual_is_diagnostic_only")
        is not True
        or arithmetic.get("threshold_change_forbidden") is not True
        or arithmetic.get("continuum_equation_change_forbidden") is not True
        or arithmetic.get("reference_derivative_map_change_forbidden") is not True
    ):
        raise ValueError("SRC2 arithmetic contract differs")
    output = config["output"]
    if (
        output.get("fresh_namespace_required") is not True
        or output.get("clean_tracked_worktree_required") is not True
        or output.get("overwrite_forbidden") is not True
    ):
        raise ValueError("SRC2 output contract is not fail-closed")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    _required_tables(config)
    return config


def _validate_lineage(config: Mapping[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    lineage = config["lineage"]
    for path_key, hash_key in (
        ("cal8_result", "cal8_result_sha256"),
        ("hlt9_authorization", "hlt9_authorization_sha256"),
        ("cal8_run_plan", "cal8_run_plan_sha256"),
        ("campaign_manifest", "campaign_manifest_sha256"),
        ("campaign_event_log", "campaign_event_log_sha256"),
        ("campaign_checkpoint", "campaign_checkpoint_sha256"),
        ("campaign_result", "campaign_result_sha256"),
    ):
        path = REPOSITORY / lineage[path_key]
        if not path.is_file() or _sha(path) != lineage[hash_key]:
            raise ValueError(f"SRC2 lineage drifted: {path_key}")
    stored_cal8 = _load_canonical_json(REPOSITORY / lineage["cal8_result"])
    if _serial(cal8.reproduce(cal8.DEFAULT_CONFIG)) != stored_cal8:
        raise ValueError("CAL8/PREF10 no longer reproduces before SRC2 capture")
    if (
        stored_cal8.get("gate_status", {}).get(
            "SRC2_affine_source_arithmetic_preflight_required"
        )
        is not True
        or stored_cal8.get("gate_status", {}).get("PROTO12_frozen") is not False
        or stored_cal8.get("gate_status", {}).get(
            "FGCQR_holdout_execution_authorized"
        )
        is not False
    ):
        raise ValueError("CAL8 does not authorize only the SRC2 diagnostic boundary")
    events = _load_canonical_jsonl(REPOSITORY / lineage["campaign_event_log"])
    return stored_cal8, events


class TargetedSourceWallCapture:
    """Delegate unchanged PROTO11 RHS calls while retaining one target state."""

    def __init__(
        self,
        delegate: Proto11GR0EvolutionOperator,
        *,
        target_time: float,
        target_residual: float,
    ) -> None:
        if not isinstance(delegate, Proto11GR0EvolutionOperator):
            raise TypeError("SRC2 capture delegate must be the PROTO11 operator")
        self.delegate = delegate
        self.target_time = float(target_time)
        self.target_residual = float(target_residual)
        self.total_calls = 0
        self.failed_calls = 0
        self.target_matches = 0
        self.captured_state: EvolutionState | None = None
        self.captured_rhs: EvolutionRHS | None = None

    def __call__(self, coordinate_time: float, state: EvolutionState) -> EvolutionRHS:
        rhs = self.delegate(coordinate_time, state)
        self.total_calls += 1
        residual = float(rhs.diagnostics["source_residual_infinity"])
        if residual >= self.delegate.raw_tolerance:
            self.failed_calls += 1
        if _bits_equal(coordinate_time, self.target_time) and _bits_equal(
            residual, self.target_residual
        ):
            self.target_matches += 1
            self.captured_state = EvolutionState(state.u, state.p, state.q)
            self.captured_rhs = rhs
        return rhs


def _expected_records(
    events: list[dict[str, Any]], config: Mapping[str, Any]
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    replay = config["replay"]
    common = [
        item
        for item in events
        if item.get("event_type") == "common_event"
        and item.get("amplitude") == "5/2"
        and item.get("event_index") == replay["first_common_event_index"]
    ]
    rejections = [
        item
        for item in events
        if item.get("event_type") == "rejected_unaccepted_source_only_proposal"
        and item.get("amplitude") == "5/2"
        and item.get("member") == "RK4-8193"
    ]
    if len(common) != 1 or len(rejections) != replay["expected_source_only_rejection_count"]:
        raise ValueError("immutable SRC2 replay records are incomplete")
    return common[0], rejections


def _verify_first_event(member: Any, common: Mapping[str, Any], replay: Mapping[str, Any]) -> None:
    diagnostic = common["assessment"]["member_diagnostics"]["RK4-8193"]
    state_hash = array_content_sha256(member.state.u, member.state.p, member.state.q)
    if (
        state_hash != replay["first_common_event_state_sha256"]
        or state_hash != diagnostic["state_sha256"]
        or member.step_index != replay["first_common_event_step_index"]
        or member.step_index != diagnostic["step_index"]
        or member.transaction_serial != replay["first_common_event_transaction_serial"]
        or member.transaction_serial != diagnostic["transaction_serial"]
        or member.CFL_retry_count != replay["first_common_event_CFL_retry_count"]
        or member.CFL_retry_count != diagnostic["CFL_retry_count"]
        or member.source_retry_count != replay["first_common_event_source_retry_count"]
        or member.source_retry_count != diagnostic["source_retry_count"]
    ):
        raise ValueError("SRC2 replay did not reproduce the first evolved common event")


def _verify_terminal(member: Any, replay: Mapping[str, Any]) -> None:
    state_hash = array_content_sha256(member.state.u, member.state.p, member.state.q)
    if (
        state_hash != replay["expected_terminal_accepted_state_sha256"]
        or not _bits_equal(member.time, float(replay["expected_terminal_accepted_time"]))
        or member.step_index != replay["expected_terminal_step_index"]
        or member.transaction_serial != replay["expected_terminal_transaction_serial"]
        or member.source_retry_count != replay["expected_terminal_source_retry_count"]
    ):
        raise ValueError("SRC2 replay did not reproduce the terminal accepted state")


def _candidate_summary(
    candidate: ArithmeticCandidate,
    baseline: ArithmeticCandidate,
) -> dict[str, Any]:
    return {
        "label": candidate.label,
        "affine_system_label": candidate.affine_system_label,
        "acceleration_sha256": array_content_sha256(candidate.acceleration),
        "complete_residual_sha256": array_content_sha256(candidate.complete_residual),
        "complete_residual_infinity": candidate.complete_residual_infinity,
        "reconstructed_residual_infinity": candidate.reconstructed_residual_infinity,
        "strict_raw_gate_passed": candidate.strict_raw_gate_passed,
        "maximum_residual_point_index": candidate.maximum_residual_point_index,
        "maximum_residual_row_index": candidate.maximum_residual_row_index,
        "maximum_residual_radius": candidate.maximum_residual_radius,
        "complete_refinement_history": list(candidate.complete_refinement_history),
        "exact_binary_refinement_steps": candidate.exact_binary_refinement_steps,
        "selected_exact_binary_points": list(candidate.selected_exact_binary_points),
        "acceleration_difference_from_baseline_infinity": float(
            np.max(
                np.abs(candidate.acceleration - baseline.acceleration),
                initial=0.0,
            )
        ),
    }


def _system_summary(system: CapturedAffineSystem) -> dict[str, Any]:
    row_norm = np.max(np.sum(np.abs(system.jacobian), axis=2), axis=1)
    inverse = np.linalg.inv(system.jacobian)
    inverse_norm = np.max(np.sum(np.abs(inverse), axis=2), axis=1)
    return {
        "label": system.label,
        "seed_scale": system.seed_scale,
        "constant_sha256": array_content_sha256(system.constant),
        "jacobian_sha256": array_content_sha256(system.jacobian),
        "constant_infinity": float(np.max(np.abs(system.constant), initial=0.0)),
        "jacobian_infinity": float(np.max(np.abs(system.jacobian), initial=0.0)),
        "kinetic_condition_infinity_maximum": float(
            np.max(row_norm * inverse_norm, initial=0.0)
        ),
    }


def _fixture_metadata(
    *,
    config_path: Path,
    config: Mapping[str, Any],
    observer: TargetedSourceWallCapture,
    comparison: AffineArithmeticComparison,
    forward: CapturedAffineSystem,
    symmetric: tuple[CapturedAffineSystem, ...],
    arrays: Mapping[str, np.ndarray],
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "runner_id": RUNNER_ID,
        "git_commit": _git_head(),
        "config_sha256": _sha(config_path),
        "branch": "GR-0",
        "amplitude": "5/2",
        "member": "RK4-8193",
        "captured_stage_time": observer.target_time,
        "captured_source_residual": observer.target_residual,
        "target_matches": observer.target_matches,
        "complete_residual_gate": comparison.raw_tolerance,
        "array_sha256": {
            name: array_content_sha256(value) for name, value in arrays.items()
        },
        "forward_system": _system_summary(forward),
        "symmetric_systems": [_system_summary(item) for item in symmetric],
        "continuum_equations_changed": False,
        "PROTO11_reference_derivative_map_changed": False,
        "raw_residual_threshold_changed": False,
        "mechanism_question_answered": False,
        "FGCQR_outcome_read": False,
        "config_claims": dict(config["claims"]),
    }


def run(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    inherited._require_clean_tracked_worktree()
    config = load_config(config_path)
    if config_path.resolve() != DEFAULT_CONFIG.resolve():
        raise ValueError("SRC2 capture accepts only the canonical config path")
    stored_cal8, events = _validate_lineage(config)
    common, expected_rejections = _expected_records(events, config)
    replay = config["replay"]
    arithmetic = config["arithmetic"]
    output = config["output"]
    output_root = REPOSITORY / output["root"]
    if output_root.exists():
        raise ValueError("fresh SRC2 capture namespace already exists")

    implementation_paths = (
        Path("scripts/run_fgc_src2_affine_capture.py"),
        Path("src/recursive_horizons/fgc/evolution/src2_affine_arithmetic.py"),
        Path("src/recursive_horizons/fgc/evolution/gr0_direct_source.py"),
        Path("src/recursive_horizons/fgc/evolution/proto11_runtime.py"),
        Path("scripts/run_fgc_gr0_calibration_v11.py"),
    )
    manifest = {
        "schema_version": 1,
        "runner_id": RUNNER_ID,
        "git_commit": _git_head(),
        "config_path": config_path.relative_to(REPOSITORY).as_posix(),
        "config_sha256": _sha(config_path),
        "cal8_result_sha256": config["lineage"]["cal8_result_sha256"],
        "PROTO11_campaign_id": stored_cal8["artifact_payload"]["immutable_campaign"][
            "campaign_id"
        ],
        "implementation_sha256": {
            path.as_posix(): _sha(REPOSITORY / path) for path in implementation_paths
        },
        "created_unix_time": time.time(),
        "outcome_fields_present_at_creation": False,
        "FGCQR_outcome_read": False,
    }
    plan, authorization = proto11_runner._validate_authorization(
        REPOSITORY / config["lineage"]["cal8_run_plan"],
        REPOSITORY / config["lineage"]["hlt9_authorization"],
        require_fresh_namespace=False,
    )
    members = proto11_runner._build_members(plan, authorization, "5/2")
    member = members.get("RK4-8193")
    if member is None or not isinstance(member.operator, Proto11GR0EvolutionOperator):
        raise ValueError("SRC2 implicated PROTO11 member is absent")
    observer = TargetedSourceWallCapture(
        member.operator,
        target_time=float(replay["expected_terminal_failure_stage_time"]),
        target_residual=float(replay["expected_terminal_failure_residual"]),
    )
    member.operator = observer
    rejections: list[dict[str, Any]] = []
    advance_arguments = {
        "retry_factor": float(Q(replay["retry_factor"])),
        "maximum_CFL_retries": replay["maximum_CFL_retries_per_step"],
        "maximum_source_retries": replay["maximum_source_retries_per_step"],
        "minimum_step_size": float(Q(replay["minimum_step_size"])),
        "rejected_trial_sink": lambda record: rejections.append(dict(record)),
    }
    started = time.monotonic()
    member.advance_to(float(Q(replay["first_common_event_time"])), **advance_arguments)
    if rejections:
        raise ValueError("SRC2 replay encountered a source rejection before the first event")
    _verify_first_event(member, common, replay)
    first_event_elapsed = time.monotonic() - started

    terminal_error: proto7_runner.Proto7SourceRetryExhausted | None = None
    try:
        member.advance_to(float(Q(replay["target_common_event_time"])), **advance_arguments)
    except proto7_runner.Proto7SourceRetryExhausted as error:
        terminal_error = error
    if terminal_error is None:
        raise ValueError("SRC2 replay did not reproduce source-retry exhaustion")
    replay_elapsed = time.monotonic() - started
    if len(rejections) != len(expected_rejections) or any(
        _canonical_line(observed) != _canonical_line(expected)
        for observed, expected in zip(rejections, expected_rejections, strict=True)
    ):
        raise ValueError("SRC2 replay rejection trace differs from the immutable ledger")
    if _canonical_line(terminal_error.evidence["last_rejected_proposal"]) != _canonical_line(
        expected_rejections[-1]
    ):
        raise ValueError("SRC2 terminal exception differs from the immutable ledger")
    _verify_terminal(member, replay)
    if (
        observer.target_matches != replay["expected_terminal_matching_source_calls"]
        or observer.captured_state is None
        or observer.captured_rhs is None
    ):
        raise ValueError("SRC2 did not capture the exact terminal source calls")

    derivative = observer.delegate.derivative
    p_r, q_r, _differentiated_u = reference_balanced_spatial_derivatives(
        observer.captured_state, derivative
    )
    captured = observer.captured_state
    jet = (
        captured.u[1:],
        captured.p[1:],
        captured.q[1:],
        p_r[1:],
        q_r[1:],
        derivative.grid.coordinates[1:],
    )
    comparison, forward, symmetric = compare_affine_arithmetic(
        *jet,
        raw_tolerance=float(Q(arithmetic["raw_complete_residual_maximum"])),
        condition_number_maximum=float(
            Q(arithmetic["kinetic_condition_number_maximum"])
        ),
        symmetric_seed_scales=tuple(
            float(Q(item)) for item in arithmetic["symmetric_power_two_seed_scales"]
        ),
        maximum_refinement_iterations=arithmetic[
            "maximum_complete_refinement_iterations"
        ],
        exact_binary_point_limit=arithmetic["exact_binary_point_limit"],
        exact_binary_selection_ratio=float(
            Q(arithmetic["exact_binary_selection_ratio"])
        ),
    )
    if not _bits_equal(
        comparison.baseline.complete_residual_infinity,
        float(replay["expected_terminal_failure_residual"]),
    ):
        raise ValueError("captured baseline residual differs from the terminal ledger")

    arrays = {
        "u": np.ascontiguousarray(jet[0]),
        "p": np.ascontiguousarray(jet[1]),
        "q": np.ascontiguousarray(jet[2]),
        "p_r": np.ascontiguousarray(jet[3]),
        "q_r": np.ascontiguousarray(jet[4]),
        "radii": np.ascontiguousarray(jet[5]),
        "forward_constant": np.ascontiguousarray(forward.constant),
        "forward_jacobian": np.ascontiguousarray(forward.jacobian),
        "baseline_acceleration": np.ascontiguousarray(comparison.baseline.acceleration),
        "baseline_complete_residual": np.ascontiguousarray(
            comparison.baseline.complete_residual
        ),
        "best_acceleration": np.ascontiguousarray(comparison.best.acceleration),
        "best_complete_residual": np.ascontiguousarray(
            comparison.best.complete_residual
        ),
    }
    fixture_metadata = _fixture_metadata(
        config_path=config_path,
        config=config,
        observer=observer,
        comparison=comparison,
        forward=forward,
        symmetric=symmetric,
        arrays=arrays,
    )
    final_paths = {
        name: REPOSITORY / output[name]
        for name in ("manifest", "replay_result", "fixture", "arithmetic_result")
    }
    for name, path in final_paths.items():
        try:
            path.relative_to(output_root)
        except ValueError as error:
            raise ValueError(f"SRC2 {name} path escapes the output root") from error
    output_root.parent.mkdir(parents=True, exist_ok=True)
    staging_root = output_root.with_name(
        f".{output_root.name}.partial-{os.getpid()}-{time.time_ns()}"
    )
    staging_root.mkdir(parents=False, exist_ok=False)

    def staged(name: str) -> Path:
        relative = final_paths[name].relative_to(output_root)
        path = staging_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    manifest_path = staged("manifest")
    _write_exclusive(manifest_path, _canonical(manifest))
    fixture_path = staged("fixture")
    _write_npz_exclusive(fixture_path, metadata=fixture_metadata, arrays=arrays)

    replay_result = {
        "schema_version": 1,
        "runner_id": RUNNER_ID,
        "git_commit": _git_head(),
        "config_sha256": _sha(config_path),
        "first_common_event_reproduced": True,
        "complete_rejection_trace_reproduced": True,
        "source_only_rejection_count": len(rejections),
        "terminal_accepted_state_reproduced": True,
        "terminal_exception_evidence_reproduced": True,
        "target_source_call_matches": observer.target_matches,
        "operator_call_count": observer.total_calls,
        "failed_operator_call_count": observer.failed_calls,
        "first_common_event_elapsed_seconds": first_event_elapsed,
        "complete_replay_elapsed_seconds": replay_elapsed,
        "fixture_path": final_paths["fixture"].relative_to(REPOSITORY).as_posix(),
        "fixture_sha256": _sha(fixture_path),
        "FGCQR_outcome_read": False,
        "mechanism_question_answered": False,
    }
    replay_path = staged("replay_result")
    _write_exclusive(replay_path, _canonical(replay_result))

    arithmetic_result = {
        "schema_version": 1,
        "runner_id": RUNNER_ID,
        "git_commit": _git_head(),
        "config_sha256": _sha(config_path),
        "fixture_sha256": _sha(fixture_path),
        "raw_complete_residual_maximum": comparison.raw_tolerance,
        "baseline": _candidate_summary(comparison.baseline, comparison.baseline),
        "candidates": [
            _candidate_summary(item, comparison.baseline)
            for item in comparison.candidates
        ],
        "best_candidate_label": comparison.best.label,
        "best_complete_residual_infinity": comparison.best.complete_residual_infinity,
        "best_candidate_strict_raw_gate_passed": (
            comparison.best.strict_raw_gate_passed
        ),
        "forward_system": _system_summary(forward),
        "symmetric_systems": [_system_summary(item) for item in symmetric],
        "complete_unredefined_binary64_residual_was_the_only_gate": True,
        "reconstructed_affine_residual_used_as_gate": False,
        "continuum_equations_changed": False,
        "PROTO11_reference_derivative_map_changed": False,
        "raw_residual_threshold_changed": False,
        "PROTO12_frozen": False,
        "GR0_calibration_completed": False,
        "FGCQR_holdout_execution_authorized": False,
        "mechanism_question_answered": False,
    }
    arithmetic_path = staged("arithmetic_result")
    _write_exclusive(arithmetic_path, _canonical(arithmetic_result))
    if output_root.exists():
        raise ValueError("fresh SRC2 capture namespace appeared during replay")
    os.replace(staging_root, output_root)
    manifest_path = final_paths["manifest"]
    replay_path = final_paths["replay_result"]
    fixture_path = final_paths["fixture"]
    arithmetic_path = final_paths["arithmetic_result"]
    return {
        "manifest": manifest_path.relative_to(REPOSITORY).as_posix(),
        "manifest_sha256": _sha(manifest_path),
        "replay_result": replay_path.relative_to(REPOSITORY).as_posix(),
        "replay_result_sha256": _sha(replay_path),
        "fixture": fixture_path.relative_to(REPOSITORY).as_posix(),
        "fixture_sha256": _sha(fixture_path),
        "arithmetic_result": arithmetic_path.relative_to(REPOSITORY).as_posix(),
        "arithmetic_result_sha256": _sha(arithmetic_path),
        "best_candidate_label": comparison.best.label,
        "best_complete_residual_infinity": comparison.best.complete_residual_infinity,
        "best_candidate_strict_raw_gate_passed": (
            comparison.best.strict_raw_gate_passed
        ),
        "PROTO12_frozen": False,
        "FGCQR_holdout_execution_authorized": False,
    }


def check(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    _stored_cal8, events = _validate_lineage(config)
    common, rejections = _expected_records(events, config)
    replay = config["replay"]
    return {
        "authorized_to_capture": True,
        "artifact_id": config["artifact_id"],
        "first_common_event_state_sha256": common["assessment"][
            "member_diagnostics"
        ]["RK4-8193"]["state_sha256"],
        "source_only_rejection_count": len(rejections),
        "terminal_accepted_state_sha256": rejections[-1][
            "accepted_state_sha256"
        ],
        "terminal_failure_residual": rejections[-1]["failed_evaluations"][-1][
            "source_residual_infinity"
        ],
        "expected_source_only_rejection_count": replay[
            "expected_source_only_rejection_count"
        ],
        "output_created_by_check": False,
        "output_namespace_exists": (REPOSITORY / config["output"]["root"]).exists(),
        "PROTO12_frozen": False,
        "FGCQR_holdout_execution_authorized": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    config_path = args.config.resolve()
    if args.check:
        print(json.dumps(check(config_path), sort_keys=True))
        return
    print(json.dumps(run(config_path), sort_keys=True))


if __name__ == "__main__":
    main()
