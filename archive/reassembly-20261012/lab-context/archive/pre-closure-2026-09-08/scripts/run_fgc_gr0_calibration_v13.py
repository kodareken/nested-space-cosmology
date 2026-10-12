#!/usr/bin/env python3
"""Run the HLT11-authorized PROTO13 GR-0 restart calibration.

The runner restores six complete members at ``t=23/16`` and continues the
unchanged PROTO12 equations.  RK4 and SSPRK3 own separate resolution ladders;
all diagnostics are recomputed per method, while the trapped sign is compared
only on exact shared physical nodes through independent Richardson intervals.
No SGB-L or FGC-QR state can be constructed by this runner.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import time as wall_time
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_hlt11_mon11 as hlt11  # noqa: E402
from scripts import reproduce_fgc_hlt10_mon10 as hlt10  # noqa: E402
from scripts import reproduce_fgc_rsp2_frz1 as rsp2_freeze  # noqa: E402
from scripts import run_fgc_gr0_calibration as inherited  # noqa: E402
from scripts import run_fgc_gr0_calibration_v6 as proto6_runner  # noqa: E402
from scripts import run_fgc_gr0_calibration_v7 as proto7_runner  # noqa: E402
from scripts import run_fgc_gr0_calibration_v12 as proto12_runner  # noqa: E402
from scripts import run_fgc_rsp2_constraint_study as rsp2_runner  # noqa: E402
from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    PROTO4_SPECTRAL_FIELD_ORDER,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    GR0RuntimeMonitorState,
    GR0RuntimeStop,
    proto5_temporal_spectral_admission,
)
from recursive_horizons.fgc.evolution.proto6_runtime import (  # noqa: E402
    restore_common_event_snapshots,
)
from recursive_horizons.fgc.evolution.proto7_runtime import (  # noqa: E402
    Proto7TerminalStop,
)
from recursive_horizons.fgc.evolution.proto13_runtime import (  # noqa: E402
    PROTO13_COMPARATOR_POINT_COUNTS,
    PROTO13_PRIMARY_POINT_COUNTS,
    proto13_gr0_common_event,
    proto13_method_owned_trapped_assessment,
)


DEFAULT_PLAN = REPOSITORY / "configs/fgc/fgc-1-pro13-run1.toml"
DEFAULT_AUTHORIZATION = REPOSITORY / "results/fgc-1-hlt11-mon11.json"
RUNNER_ID = "FGC-1-CAL10-RUN1-RUNNER"
RESTART_TIME = 23.0 / 16.0


def _validate_authorization(
    plan_path: Path,
    authorization_path: Path,
    *,
    require_fresh_namespace: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    plan = hlt11.validate_run_plan(plan_path)
    runtime_environment = hlt11.validate_numerical_runtime(plan)
    authorization = hlt11.load_canonical_result(authorization_path)
    gates = authorization.get("gate_status", {})
    payload = authorization.get("artifact_payload", {})
    frozen = payload.get("frozen_run_plan", {})
    if (
        authorization.get("artifact_id") != hlt11.ARTIFACT_ID
        or gates.get("PROTO13_successor_runtime_implemented") is not True
        or gates.get("PROTO13_fresh_GR0_dynamic_calibration_authorized") is not True
        or gates.get("classical_spherical_diagnostic_authorized") is not False
        or gates.get("SGBL_execution_authorized") is not False
        or gates.get("FGCQR_holdout_execution_authorized") is not False
        or frozen.get("run_plan_sha256") != hlt11._sha(plan_path)
        or frozen.get("primary_point_counts") != [2049, 4097, 8193]
        or frozen.get("comparator_point_counts") != [4097, 8193, 16385]
        or frozen.get("restart_coordinate_time") != RESTART_TIME
        or frozen.get("numerical_runtime_contract") != plan["numerical_runtime"]
        or payload.get("runtime_environment_observed") != runtime_environment
    ):
        raise ValueError("HLT11 authorization is absent, incomplete, or over-broad")
    for ledger_name in ("source_config_sha256", "predecessor_sha256", "implementation_sha256"):
        ledger = authorization.get(ledger_name)
        if not isinstance(ledger, dict) or not ledger:
            raise ValueError(f"HLT11 {ledger_name} is absent")
        for relative, expected in ledger.items():
            path = REPOSITORY / relative
            if not path.is_file() or hlt11._sha(path) != expected:
                raise ValueError(f"HLT11-authorized file drifted: {relative}")
    owner = REPOSITORY / authorization.get("derivation_document", "")
    if not owner.is_file() or hlt11._sha(owner) != authorization.get("derivation_document_sha256"):
        raise ValueError("HLT11 owner document drifted")
    if require_fresh_namespace:
        reproduced = hlt11.record(hlt11.DEFAULT_CONFIG)
        if inherited._serial(reproduced) != authorization:
            raise ValueError("HLT11 authorization does not reproduce before launch")
    return plan, authorization


def _build_frozen_member_shells() -> dict[str, proto7_runner.Proto7RunMember]:
    base_plan = hlt10.validate_run_plan(REPOSITORY / "configs/fgc/fgc-1-cal9-run1.toml")
    base_authorization = hlt10.load_canonical_result(
        REPOSITORY / "results/fgc-1-hlt10-mon10.json"
    )
    base = proto12_runner._build_members(base_plan, base_authorization, "3")
    desired = {
        "RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193"
    }
    members = {key: value for key, value in base.items() if key in desired}
    rsp2_plan = rsp2_freeze.validate_run_plan(
        REPOSITORY / "configs/fgc/fgc-1-rsp2-run1.toml"
    )
    rsp2_authorization = rsp2_freeze.load_canonical_result(
        REPOSITORY / "results/fgc-1-rsp2-frz1.json"
    )
    newest = rsp2_runner._build_member(rsp2_plan, rsp2_authorization)
    members[newest.key] = newest
    if set(members) != {
        "RK4-2049", "RK4-4097", "RK4-8193",
        "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385",
    }:
        raise RuntimeError("PROTO13 member shells differ from the frozen ladders")
    return members


def _restore_frozen_members(
    authorization: Mapping[str, Any],
) -> dict[str, proto7_runner.Proto7RunMember]:
    members = _build_frozen_member_shells()
    freeze = json.loads(
        (REPOSITORY / "results/fgc-1-pro13-frz1.json").read_text(encoding="utf-8")
    )
    expected = {
        item["key"]: item
        for item in freeze["artifact_payload"]["restart_members"]
    }
    checkpoint_paths = {
        "PROTO12": REPOSITORY / "runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz",
        "RSP2": REPOSITORY / "runs/fgc-2-sf1/rsp2/amplitude3-ssprk3-momentum/latest-checkpoint.npz",
    }
    handles = {
        name: np.load(path, allow_pickle=False) for name, path in checkpoint_paths.items()
    }
    try:
        metadata = {
            name: json.loads(bytes(handle["metadata_utf8"]).decode("utf-8"))
            for name, handle in handles.items()
        }
        for key, member in members.items():
            frozen = expected[key]
            source = frozen["source_checkpoint"]
            handle = handles[source]
            stored = metadata[source]["members"][key]
            prefix = key.replace("-", "_")
            arrays = {
                suffix: handle[f"{prefix}_{suffix}"].copy()
                for suffix in (
                    "u", "p", "q", "tracer_positions", "tracer_proper_times",
                    "event_proper_times", "event_fields",
                )
            }
            state = EvolutionState(arrays["u"], arrays["p"], arrays["q"])
            if (
                array_content_sha256(state.u, state.p, state.q) != frozen["state_sha256"]
                or array_content_sha256(*arrays.values()) != frozen["restart_payload_sha256"]
                or stored["input_hash"] != member.input_hash
                or stored["input_hash"] != frozen["input_hash"]
            ):
                raise ValueError(f"PROTO13 frozen restart differs: {key}")
            member.state = state
            member.time = float(stored["time"])
            member.step_index = int(stored["step_index"])
            member.transaction_serial = int(stored["transaction_serial"])
            member.CFL_retry_count = int(stored["CFL_retry_count"])
            member.source_retry_count = int(stored["source_retry_count"])
            member.transaction.state = GR0RuntimeMonitorState(
                **stored["runtime_monitor_state"]
            )
            member.transaction.causal_state = CausalBudgetState(**stored["causal_state"])
            member.tracers.positions = arrays["tracer_positions"]
            member.tracers.proper_times = arrays["tracer_proper_times"]
            member.tracers.event_proper_times = [
                row.copy() for row in arrays["event_proper_times"]
            ]
            member.tracers.event_fields = [
                row.copy() for row in arrays["event_fields"]
            ]
            before = array_content_sha256(state.u, state.p, state.q)
            rhs = member.operator(RESTART_TIME, member.state)
            if (
                rhs.diagnostics.get("source_raw_gate_passed") is not True
                or rhs.diagnostics.get("SRC4_tensor_contracted_reference_source") is not True
                or array_content_sha256(state.u, state.p, state.q) != before
                or member.time != RESTART_TIME
                or len(member.tracers.event_fields) != 24
            ):
                raise ValueError(f"PROTO13 restart source premise failed: {key}")
    finally:
        for handle in handles.values():
            handle.close()
    return members


def _method_members(
    members: Mapping[str, proto7_runner.Proto7RunMember], method: str
) -> tuple[proto7_runner.Proto7RunMember, ...]:
    expected = (
        PROTO13_PRIMARY_POINT_COUNTS if method == "RK4" else PROTO13_COMPARATOR_POINT_COUNTS
    )
    records = tuple(
        sorted(
            (item for item in members.values() if item.method_label == method),
            key=lambda item: item.point_count,
        )
    )
    if tuple(item.point_count for item in records) != expected:
        raise ValueError(f"PROTO13 {method} member ownership differs")
    return records


def _event_assessment(
    plan: Mapping[str, Any],
    members: Mapping[str, proto7_runner.Proto7RunMember],
    coordinate_time: float,
) -> dict[str, Any]:
    for member in members.values():
        if np.float64(member.time).tobytes() != np.float64(coordinate_time).tobytes():
            raise ValueError("PROTO13 common-event member is not time aligned")
    primary = _method_members(members, "RK4")
    comparator = _method_members(members, "SSPRK3")

    def common(records: tuple[proto7_runner.Proto7RunMember, ...], method: str):
        return proto13_gr0_common_event(
            [item.state for item in records],
            [item.initial.grid for item in records],
            accepted_stage_counts=[item.transaction.state.accepted_stage_count for item in records],
            method=method,
            coordinate_time=coordinate_time,
            cutoff=16.0,
            measurement_radius_maximum=24.0,
            taper_fraction=1.0 / 8.0,
            fixed_outer_rows=4,
        )

    primary_common = common(primary, "RK4")
    comparator_common = common(comparator, "SSPRK3")
    trapped = proto13_method_owned_trapped_assessment(
        primary_states=[item.state for item in primary],
        primary_grids=[item.initial.grid for item in primary],
        comparator_states=[item.state for item in comparator],
        comparator_grids=[item.initial.grid for item in comparator],
        primary_common_event=primary_common,
        comparator_common_event=comparator_common,
        measurement_radius_maximum=24.0,
        minimum_observed_order=1.5,
        positive_margin_factor=4.0,
    )
    sample_count = len(primary[0].tracers.event_fields)
    temporal: dict[str, Any] = {
        "available": False,
        "required_minimum_samples": 64,
        "sample_count": sample_count,
        "admission_passed": False,
    }
    if sample_count >= 64:
        methods = {}
        for method, records in (("RK4", primary), ("SSPRK3", comparator)):
            methods[method] = proto5_temporal_spectral_admission(
                point_counts=[item.point_count for item in records],
                proper_times_by_resolution=[
                    np.asarray(item.tracers.event_proper_times) for item in records
                ],
                field_histories_by_resolution=[
                    np.asarray(item.tracers.event_fields) for item in records
                ],
                field_names=PROTO4_SPECTRAL_FIELD_ORDER,
                cutoff=16.0,
                taper_fraction=1.0 / 8.0,
            )
        temporal = {
            "available": True,
            "sample_count": sample_count,
            "methods": methods,
            "admission_passed": all(item.admission_passed for item in methods.values()),
        }
    qualified = trapped.trapped_sign_passed and temporal["admission_passed"]
    diagnostics = {}
    for item in members.values():
        boundary = item.transaction.boundary_geometry
        causal = item.transaction.causal_state
        diagnostics[item.key] = {
            "step_index": item.step_index,
            "transaction_serial": item.transaction_serial,
            "accepted_stage_count": item.transaction.state.accepted_stage_count,
            "CFL_retry_count": item.CFL_retry_count,
            "source_retry_count": item.source_retry_count,
            "accumulated_characteristic_distance": causal.accumulated_characteristic_distance,
            "remaining_boundary_margin": (
                boundary.outer_radius - boundary.measurement_radius
                - boundary.sbp_stencil_reach
                - causal.accumulated_characteristic_distance
                - boundary.minimum_causal_buffer
            ),
            "state_sha256": array_content_sha256(item.state.u, item.state.p, item.state.q),
            "minimum_tracer_radius": float(np.min(item.tracers.positions)),
            "maximum_tracer_radius": float(np.max(item.tracers.positions)),
        }
    return {
        "coordinate_time": coordinate_time,
        "primary_common_event": primary_common,
        "comparator_common_event": comparator_common,
        "temporal_spectral_admission": temporal,
        "trapped_assessment": trapped,
        "qualified_trapped_common_event": qualified,
        "method_owned_ladders_enforced": True,
        "common_physical_node_comparison_enforced": True,
        "member_diagnostics": diagnostics,
    }


def _manifest(
    plan_path: Path,
    authorization_path: Path,
    authorization: Mapping[str, Any],
) -> dict[str, Any]:
    frozen = authorization["artifact_payload"]["frozen_run_plan"]
    return {
        "schema_version": 1,
        "runner_id": RUNNER_ID,
        "campaign_id": frozen["campaign_id"],
        "authorization_checkpoint_commit": inherited._git_head(),
        "plan_path": _relative(plan_path),
        "plan_sha256": hlt11._sha(plan_path),
        "authorization_path": _relative(authorization_path),
        "authorization_sha256": hlt11._sha(authorization_path),
        "implementation_sha256": authorization["implementation_sha256"],
        "numerical_runtime_contract": dict(
            authorization["artifact_payload"]["frozen_run_plan"][
                "numerical_runtime_contract"
            ]
        ),
        "runtime_environment": hlt11.validate_numerical_runtime(
            hlt11.validate_run_plan(plan_path)
        ),
        "restart_coordinate_time": RESTART_TIME,
        "created_unix_time": wall_time.time(),
        "outcome_fields_present_at_creation": False,
    }


def _relative(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _member_metadata(member: proto7_runner.Proto7RunMember) -> dict[str, Any]:
    value = proto6_runner._member_metadata(member)
    value["state_sha256"] = array_content_sha256(member.state.u, member.state.p, member.state.q)
    return value


def _checkpoint_metadata(
    *,
    manifest: Mapping[str, Any],
    event_index: int,
    consecutive: int,
    event_log: Path,
    members: Mapping[str, proto7_runner.Proto7RunMember],
    terminal: bool,
    terminal_result: Mapping[str, Any] | None,
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "campaign_id": manifest["campaign_id"],
        "manifest_sha256": sha256(inherited._canonical(manifest)).hexdigest(),
        "amplitude": "3",
        "amplitude_index": 0,
        "completed_common_event_index": event_index,
        "consecutive_qualified_events": consecutive,
        "completed_amplitude_records": [],
        "event_log": inherited._event_log_identity(inherited._event_log_bytes(event_log)),
        "terminal": terminal,
        "terminal_result": terminal_result,
        "members": {key: _member_metadata(item) for key, item in members.items()},
    }


def _restore_campaign_checkpoint(
    path: Path,
    *,
    manifest: Mapping[str, Any],
    authorization: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, proto7_runner.Proto7RunMember], bytes]:
    members = _build_frozen_member_shells()
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(bytes(archive["metadata_utf8"]).decode("utf-8"))
        event_log = bytes(archive["event_log_utf8"])
        if (
            metadata.get("campaign_id") != manifest["campaign_id"]
            or metadata.get("manifest_sha256") != sha256(inherited._canonical(manifest)).hexdigest()
            or metadata.get("event_log") != inherited._event_log_identity(event_log)
            or set(metadata.get("members", {})) != set(members)
        ):
            raise ValueError("PROTO13 checkpoint identity differs")
        for key, member in members.items():
            prefix = key.replace("-", "_")
            stored = metadata["members"][key]
            if stored["input_hash"] != member.input_hash:
                raise ValueError(f"PROTO13 checkpoint input hash differs: {key}")
            member.state = EvolutionState(
                archive[f"{prefix}_u"].copy(),
                archive[f"{prefix}_p"].copy(),
                archive[f"{prefix}_q"].copy(),
            )
            if array_content_sha256(member.state.u, member.state.p, member.state.q) != stored["state_sha256"]:
                raise ValueError(f"PROTO13 checkpoint state hash differs: {key}")
            member.time = float(stored["time"])
            member.step_index = int(stored["step_index"])
            member.transaction_serial = int(stored["transaction_serial"])
            member.CFL_retry_count = int(stored["CFL_retry_count"])
            member.source_retry_count = int(stored["source_retry_count"])
            member.transaction.state = GR0RuntimeMonitorState(**stored["runtime_monitor_state"])
            member.transaction.causal_state = CausalBudgetState(**stored["causal_state"])
            member.tracers.positions = archive[f"{prefix}_tracer_positions"].copy()
            member.tracers.proper_times = archive[f"{prefix}_tracer_proper_times"].copy()
            member.tracers.event_proper_times = [row.copy() for row in archive[f"{prefix}_event_proper_times"]]
            member.tracers.event_fields = [row.copy() for row in archive[f"{prefix}_event_fields"]]
    return metadata, members, event_log


def _validate_resume_manifest(
    manifest: Mapping[str, Any],
    plan_path: Path,
    authorization_path: Path,
    authorization: Mapping[str, Any],
) -> None:
    frozen = authorization["artifact_payload"]["frozen_run_plan"]
    if (
        manifest.get("runner_id") != RUNNER_ID
        or manifest.get("campaign_id") != frozen["campaign_id"]
        or manifest.get("plan_sha256") != hlt11._sha(plan_path)
        or manifest.get("authorization_sha256") != hlt11._sha(authorization_path)
        or manifest.get("implementation_sha256") != authorization["implementation_sha256"]
        or manifest.get("authorization_checkpoint_commit") != inherited._git_head()
        or manifest.get("numerical_runtime_contract") != frozen["numerical_runtime_contract"]
        or manifest.get("runtime_environment") != hlt11.validate_numerical_runtime(
            hlt11.validate_run_plan(plan_path)
        )
    ):
        raise ValueError("PROTO13 resume manifest differs from authorization")


def run_campaign(
    *, plan_path: Path, authorization_path: Path, resume: bool
) -> dict[str, Any]:
    inherited._require_clean_tracked_worktree()
    plan, authorization = _validate_authorization(
        plan_path, authorization_path, require_fresh_namespace=not resume
    )
    provenance = plan["provenance"]
    output_root = REPOSITORY / provenance["output_root"]
    manifest_path = REPOSITORY / provenance["campaign_manifest"]
    event_log = output_root / provenance["append_only_event_log_name"]
    checkpoint_path = output_root / provenance["checkpoint_name"]
    result_path = output_root / provenance["result_name"]
    resume_recovery = None

    if resume:
        if not output_root.is_dir() or not manifest_path.is_file() or result_path.exists():
            raise ValueError("PROTO13 resume requires one nonterminal campaign")
        manifest = inherited._load_manifest(manifest_path)
        _validate_resume_manifest(manifest, plan_path, authorization_path, authorization)
        metadata, members, saved_log = _restore_campaign_checkpoint(
            checkpoint_path, manifest=manifest, authorization=authorization
        )
        if metadata.get("terminal") is True:
            raise ValueError("PROTO13 terminal checkpoint cannot resume")
        event_index = int(metadata["completed_common_event_index"])
        consecutive = int(metadata["consecutive_qualified_events"])
        resume_recovery = inherited._recover_event_log_from_checkpoint(event_log, saved_log)
    else:
        if output_root.exists():
            raise ValueError("fresh PROTO13 output root already exists")
        # Restore and assess every frozen premise before creating the namespace.
        members = _restore_frozen_members(authorization)
        initial = _event_assessment(plan, members, RESTART_TIME)
        if (
            not initial["primary_common_event"].admission_passed
            or not initial["comparator_common_event"].admission_passed
            or initial["temporal_spectral_admission"]["available"]
        ):
            raise ValueError("PROTO13 restart admission changed after HLT11")
        output_root.mkdir(parents=True, exist_ok=False)
        manifest = _manifest(plan_path, authorization_path, authorization)
        inherited._write_exclusive(manifest_path, inherited._canonical(manifest))
        event_index = 23
        consecutive = 0
        inherited._append_fsync(
            event_log,
            inherited._canonical_line({
                "event_type": "restart_common_event",
                "event_index": event_index,
                "assessment": initial,
                "trajectory_advanced_by_PROTO13": False,
                "counts_toward_new_consecutive_trapped_events": False,
            }),
        )
        checkpoint = _checkpoint_metadata(
            manifest=manifest, event_index=event_index, consecutive=consecutive,
            event_log=event_log, members=members, terminal=False, terminal_result=None,
        )
        proto6_runner._write_checkpoint(
            checkpoint_path, metadata=checkpoint, members=members, event_log_path=event_log
        )

    output_interval = 1.0 / 16.0
    checkpoint_every = 4
    final_event_index = 512
    required_consecutive = 8
    retry_factor = 0.5
    minimum_step = 1.0 / 1073741824.0
    candidate_stop = None
    selected = False
    started = wall_time.monotonic()

    while event_index < final_event_index:
        target_index = event_index + 1
        target_time = target_index * output_interval
        snapshots = {key: item.snapshot() for key, item in members.items()}
        active = None

        def rejected_trial_sink(evidence: Mapping[str, Any]) -> None:
            inherited._append_fsync(event_log, inherited._canonical_line(evidence))

        try:
            for active in members.values():
                active.advance_to(
                    target_time,
                    retry_factor=retry_factor,
                    maximum_CFL_retries=32,
                    maximum_source_retries=32,
                    minimum_step_size=minimum_step,
                    rejected_trial_sink=rejected_trial_sink,
                )
            for member in members.values():
                member.append_common_event()
            assessment = _event_assessment(plan, members, target_time)
        except Proto7TerminalStop as stop:
            restore_common_event_snapshots(members, snapshots)
            candidate_stop = {
                "classification": "scientific_universal_runtime_stop",
                "reason": stop.reason,
                "failures": list(stop.failures),
                "complete_failure_evidence": asdict(stop.evidence),
                "member": None if active is None else active.key,
                "target_common_event_index": target_index,
                "last_fully_completed_common_event_index": event_index,
                "cross_member_state_rolled_back": True,
            }
            break
        except proto7_runner.Proto7SourceRetryExhausted as stop:
            restore_common_event_snapshots(members, snapshots)
            candidate_stop = {
                "classification": "scientific_source_retry_exhausted",
                "reason": stop.reason,
                "complete_failure_evidence": stop.evidence,
                "member": None if active is None else active.key,
                "target_common_event_index": target_index,
                "last_fully_completed_common_event_index": event_index,
                "cross_member_state_rolled_back": True,
            }
            break
        except GR0RuntimeStop as stop:
            restore_common_event_snapshots(members, snapshots)
            candidate_stop = {
                "classification": "invalid_implementation_or_nonconverged_run",
                "reason": "terminal_stop_without_complete_PROTO7_evidence",
                "inherited_reason": stop.reason,
                "failures": list(stop.failures),
                "target_common_event_index": target_index,
                "last_fully_completed_common_event_index": event_index,
            }
            break
        except BaseException as error:
            if isinstance(error, KeyboardInterrupt):
                raise
            restore_common_event_snapshots(members, snapshots)
            candidate_stop = {
                "classification": "invalid_implementation_or_nonconverged_run",
                "reason": inherited._source_failure_reason(error),
                "error_type": type(error).__name__,
                "error_message": str(error),
                "target_common_event_index": target_index,
                "last_fully_completed_common_event_index": event_index,
            }
            break

        event_index = target_index
        consecutive = consecutive + 1 if assessment["qualified_trapped_common_event"] else 0
        inherited._append_fsync(
            event_log,
            inherited._canonical_line({
                "event_type": "common_event", "amplitude": "3",
                "event_index": event_index, "assessment": assessment,
                "consecutive_qualified_trapped_common_events": consecutive,
            }),
        )
        if event_index % checkpoint_every == 0 or consecutive > 0:
            checkpoint = _checkpoint_metadata(
                manifest=manifest, event_index=event_index, consecutive=consecutive,
                event_log=event_log, members=members, terminal=False, terminal_result=None,
            )
            proto6_runner._write_checkpoint(
                checkpoint_path, metadata=checkpoint, members=members, event_log_path=event_log
            )
        if event_index % 8 == 0 or consecutive > 0:
            trapped = assessment["trapped_assessment"]
            print(
                f"PROTO13 event={event_index}/{final_event_index} t={target_time:.6g} "
                f"trapped_score={trapped.primary_trapped_scores[-1]:.6g}/"
                f"{trapped.comparator_trapped_scores[-1]:.6g} "
                f"qualified_run={consecutive}/{required_consecutive}",
                flush=True,
            )
        if consecutive >= required_consecutive:
            selected = True
            break
        common_pass = (
            assessment["primary_common_event"].admission_passed
            and assessment["comparator_common_event"].admission_passed
        )
        temporal = assessment["temporal_spectral_admission"]
        if not common_pass or (temporal["available"] and not temporal["admission_passed"]):
            candidate_stop = {
                "classification": "common_event_constraint_or_spectral_stop",
                "reason": (
                    "common_event_constraint_or_spatial_spectral_admission"
                    if not common_pass else "causal_past_temporal_spectral_admission"
                ),
                "event_index": event_index,
                "coordinate_time": target_time,
            }
            break

    if selected:
        classification = "GR0_calibration_completed_selected_amplitude"
    elif candidate_stop and candidate_stop["classification"] == "invalid_implementation_or_nonconverged_run":
        classification = "invalid_implementation_or_nonconverged_run"
    else:
        classification = "calibration_failed_no_eligible_GR0_case"
    record = {
        "amplitude": "3",
        "last_completed_common_event_index": event_index,
        "last_completed_coordinate_time": event_index * output_interval,
        "consecutive_qualified_trapped_common_events": consecutive,
        "eligible": selected,
        "stop": candidate_stop,
        "member_source_retry_counts": {key: item.source_retry_count for key, item in members.items()},
        "member_CFL_retry_counts": {key: item.CFL_retry_count for key, item in members.items()},
        "member_final_state_hashes": {
            key: array_content_sha256(item.state.u, item.state.p, item.state.q)
            for key, item in members.items()
        },
    }
    result = {
        "schema_version": 1,
        "runner_id": RUNNER_ID,
        "campaign_id": manifest["campaign_id"],
        "classification": classification,
        "selected_amplitude": "3" if selected else None,
        "amplitude_records": [record],
        "event_log": _relative(event_log),
        "event_log_sha256": hlt11._sha(event_log),
        "elapsed_wall_seconds": wall_time.monotonic() - started,
        "resume_recovery": resume_recovery,
        "restart_coordinate_time": RESTART_TIME,
        "numerical_runtime_contract": manifest["numerical_runtime_contract"],
        "runtime_environment": manifest["runtime_environment"],
        "method_owned_ladders_enabled": True,
        "common_physical_node_comparison_enabled": True,
        "FGCQR_outcome_read": False,
        "SGBL_outcome_read": False,
        "holdout_execution_authorized": False,
        "retained_EFT_evolution_authorized": False,
        "mechanism_question_answered": False,
    }
    inherited._write_exclusive(result_path, inherited._canonical(result))
    checkpoint = _checkpoint_metadata(
        manifest=manifest, event_index=event_index, consecutive=consecutive,
        event_log=event_log, members=members, terminal=True, terminal_result=result,
    )
    proto6_runner._write_checkpoint(
        checkpoint_path, metadata=checkpoint, members=members, event_log_path=event_log
    )
    return inherited._serial(result)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--authorization", type=Path, default=DEFAULT_AUTHORIZATION)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--check-authorization-only", action="store_true")
    args = parser.parse_args()
    if args.check_authorization_only:
        inherited._require_clean_tracked_worktree()
        plan, authorization = _validate_authorization(
            args.plan.resolve(), args.authorization.resolve(), require_fresh_namespace=True
        )
        print(json.dumps({
            "authorized": True,
            "campaign_id": authorization["artifact_payload"]["frozen_run_plan"]["campaign_id"],
            "amplitude": plan["scope"]["amplitude"],
            "primary_point_counts": plan["method_owned_ladders"]["primary_point_counts"],
            "comparator_point_counts": plan["method_owned_ladders"]["comparator_point_counts"],
            "output_created": False,
        }, sort_keys=True))
        return
    result = run_campaign(
        plan_path=args.plan.resolve(),
        authorization_path=args.authorization.resolve(),
        resume=args.resume,
    )
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
