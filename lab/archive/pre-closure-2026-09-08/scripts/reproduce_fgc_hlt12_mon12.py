#!/usr/bin/env python3
"""Reproduce the no-trajectory PROTO14 runtime authorization.

HLT12/MON12 validates the immutable PROTO14 restart inputs and the exact
runtime which will consume TDG6.  It deliberately creates neither PROTO14
output namespace, advances no state, and opens no candidate branch.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_pro14_frz1 as pro14  # noqa: E402


ARTIFACT_ID = "FGC-1-HLT12-MON12"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hlt12-mon12.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hlt12-mon12.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hlt12-mon12.md"
RUN_PLAN = REPOSITORY / "configs/fgc/fgc-1-pro14-run1.toml"
INTEGRATION_TESTS = (
    ("synthetic_fine_only_commit", "tests.test_fgc_proto14_runtime.Proto14RuntimeTests.test_fine_only_admission_advances_the_inherited_member"),
    ("durable_temporal_rejection_before_ledger_adoption", "tests.test_fgc_proto14_runtime.Proto14RuntimeTests.test_temporal_rejection_is_durable_before_ledger_adoption_and_retry"),
    ("inherited_CFL_retry_routing", "tests.test_fgc_proto14_runtime.Proto14RuntimeTests.test_wrapped_cfl_refusal_remains_an_inherited_cfl_retry"),
    ("late_fine_failure_rolls_back_real_ledgers", "tests.test_fgc_tdg6_temporal_admission_runtime.TDG6TemporalAdmissionRuntimeTests.test_late_fine_non_source_failure_preserves_every_real_ledger"),
    ("TDG6_checkpoint_round_trip_and_mutation_rejection", "tests.test_fgc_tdg6_temporal_admission_runtime.TDG6TemporalAdmissionRuntimeTests.test_checkpoint_round_trip_covers_debit_retry_and_rejections"),
    ("inherited_source_retry_ownership", "tests.test_fgc_tdg6_temporal_admission_runtime.TDG6TemporalAdmissionRuntimeTests.test_source_only_shadow_failure_retains_PROTO7_ownership"),
    ("recoverable_terminal_checkpoint_publication", "tests.test_fgc_gr0_campaign_runner_v14.FGCGR0CampaignRunnerV14Tests.test_terminal_checkpoint_publication_is_idempotently_recoverable"),
)

EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO14",
    "predecessor_freeze": "FGC-1-PRO14-FRZ1",
    "runtime_owner": ARTIFACT_ID,
    "calibration_branch": "GR-0",
    "calibration_amplitude": "3",
    "role": "pretrajectory_TDG6_admitted_restart_runtime_authorization",
    "PROTO13_and_CAL10_history_disclosed": True,
    "PROTO14_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "numerical_premise_promoted_to_physical_mechanism": False,
}

EXPECTED_RUNTIME = {
    "python_version": "3.14.3",
    "python_implementation": "CPython",
    "numpy_version": "2.5.1",
    "blas_name": "accelerate",
    "blas_version": "unknown",
    "system": "Darwin",
    "machine": "arm64",
    "full_runtime_metadata_must_be_recorded_in_manifest_and_result": True,
    "runtime_drift_must_fail_before_source_or_state_advance": True,
}

EXPECTED_PLAN_RUNTIME = dict(EXPECTED_RUNTIME)

EXPECTED_CLAIMS = {
    "PROTO14_successor_runtime_implemented": True,
    "PROTO14_fresh_GR0_dynamic_calibration_authorized": True,
    "PROTO14_resolved_holdout_manifest_authorized": False,
    "fresh_GR0_dynamic_calibration_completed": False,
    "GR0_case_eligible": False,
    "classical_spherical_diagnostic_authorized": False,
    "SGBL_execution_authorized": False,
    "FGCQR_holdout_execution_authorized": False,
    "DEF1_execution_authorized": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
    "FGCQR_mechanism_rejected": False,
    "general_gradient_route_rejected": False,
    "singularity_resolution_derived": False,
    "child_domain_or_topology_derived": False,
    "dark_sector_mechanism_derived": False,
    "varying_locally_measured_c_derived": False,
}

EXPECTED_RESTART = {
    "restart_coordinate_time": "23/16",
    "primary_point_counts": [2049, 4097, 8193],
    "comparator_point_counts": [4097, 8193, 16385],
    "expected_member_count": 6,
    "expected_event_sample_count": 24,
    "expected_tracer_count": 48,
    "all_payloads_must_restore_nonmutatingly": True,
    "TDG6_ledger_must_start_zero_at_restart": True,
    "pre_restart_temporal_error_bound_proved": False,
    "restart_event_cannot_count_toward_new_consecutive_trapped_events": True,
    "no_state_reinitialization_reprojection_interpolation_or_refit": True,
}

EXPECTED_NAMESPACE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto14/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto14/holdout",
    "both_roots_must_be_absent_at_authorization": True,
    "authorization_must_create_neither_root": True,
    "runner_must_refuse_overwrite": True,
}

EXPECTED_RUNTIME_BINDING = {
    "successor_runtime_module": "src/recursive_horizons/fgc/evolution/proto14_runtime.py",
    "authorized_runner": "scripts/run_fgc_gr0_calibration_v14.py",
    "authorized_runner_id": "FGC-1-CAL11-RUN1-RUNNER",
    "TDG6_module_is_consumed_directly": True,
    "every_accepted_post_restart_macro_step_uses_one_two_four_admission": True,
    "only_admitted_four_quarter_fine_path_commits": True,
    "TDG6_ledger_and_complete_temporal_rejections_checkpoint_round_trip": True,
    "source_and_CFL_retry_ownership_remains_unchanged": True,
    "sampled_64_history_is_nonveto_diagnostic_only": True,
    "runner_refuses_dirty_tracked_state_namespace_reuse_and_mismatched_resume": True,
    "exact_numerical_runtime_contract_is_bound_before_source_or_state_advance": True,
    "terminal_checkpoint_precedes_result_and_resume_reconciles_idempotently": True,
}

EXPECTED_PROOF_CONTRACT = {
    "PROTO14_and_PRO14_FRZ1_must_validate": True,
    "all_tracked_predecessor_hashes_must_match": True,
    "local_authorization_generation_requires_complete_hash_exact_raw_restart_bundle": True,
    "raw_restart_bundle_optional_but_explicitly_unverified_in_portable_clean_clone_audit": True,
    "partial_or_hash_mismatched_raw_restart_bundle_fails_closed": True,
    "all_six_complete_restart_payloads_must_match_PROTO14": True,
    "restart_restoration_must_not_mutate_any_payload": True,
    "TDG6_protocol_runtime_and_successor_runner_must_be_hash_bound": True,
    "runtime_fingerprint_must_match_exactly": True,
    "threshold_restart_history_claim_and_namespace_mutations_must_fail_closed": True,
    "fresh_namespaces_must_be_observed_absent_without_creation": True,
    "canonical_hash_bound_result_required": True,
    "no_PROTO14_trajectory_may_be_advanced": True,
}

EXPECTED_RUN_PLAN_CLAIMS = {
    key: False
    for key in (
        "fresh_GR0_dynamic_calibration_completed",
        "GR0_case_eligible",
        "classical_spherical_diagnostic_authorized",
        "SGBL_execution_authorized",
        "FGCQR_holdout_execution_authorized",
        "retained_EFT_evolution_authorized",
        "physical_transition_claim_authorized",
        "FGCQR_mechanism_rejected",
        "general_gradient_route_rejected",
        "singularity_resolution_derived",
        "child_domain_or_topology_derived",
        "dark_sector_mechanism_derived",
        "varying_locally_measured_c_derived",
    )
}


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _sha(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n"


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        Path(temporary).replace(path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def _git(*args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *args], cwd=REPOSITORY, check=check, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        value = tomllib.load(handle)
    required = {
        "schema_version", "artifact_id", "project_version", "metric_signature",
        "riemann_convention", "protocol_config", "protocol_freeze_config",
        "protocol_freeze_result", "run_plan_config", "scope", "immutable_lineage", "numerical_runtime",
        "runtime_binding", "restart_admission", "namespace", "proof_contract", "claims",
    }
    if set(value) != required:
        raise ValueError("HLT12 config root differs")
    if (
        value["schema_version"] != 1
        or value["artifact_id"] != ARTIFACT_ID
        or value["project_version"] != PROJECT_VERSION
        or value["metric_signature"] != "-+++"
        or value["riemann_convention"] != "plus_partial_mu_gamma_nu"
        or value["protocol_config"] != "configs/fgc/fgc-2-sf1-protocol-v14.toml"
        or value["protocol_freeze_config"] != "configs/fgc/fgc-1-pro14-frz1.toml"
        or value["protocol_freeze_result"] != "results/fgc-1-pro14-frz1.json"
        or value["run_plan_config"] != "configs/fgc/fgc-1-pro14-run1.toml"
        or value["scope"] != EXPECTED_SCOPE
        or value["numerical_runtime"] != EXPECTED_RUNTIME
        or value["restart_admission"] != EXPECTED_RESTART
        or value["namespace"] != EXPECTED_NAMESPACE
        or value["claims"] != EXPECTED_CLAIMS
        or value["runtime_binding"] != EXPECTED_RUNTIME_BINDING
        or value["proof_contract"] != EXPECTED_PROOF_CONTRACT
    ):
        raise ValueError("HLT12 config violates its frozen scope")
    return value


def validate_run_plan(path: Path = RUN_PLAN) -> dict[str, Any]:
    """Validate the frozen PROTO14 run contract before it is authorized."""

    with path.open("rb") as handle:
        value = tomllib.load(handle)
    required = {
        "schema_version", "artifact_id", "project_version", "metric_signature",
        "riemann_convention", "protocol_config", "protocol_freeze_result",
        "runtime_authorization_result", "inherited_run_plan",
        "inherited_runtime_authorization", "tdg6_runtime_result", "scope",
        "numerical_runtime", "method_owned_ladders", "restart", "temporal_admission",
        "schedule", "assessment", "provenance", "claims",
    }
    if set(value) != required:
        raise ValueError("PROTO14 run-plan root differs")
    expected_scope = {
        "target_protocol": "FGC-2-SF1-PROTO14", "branch": "GR-0",
        "role": "TDG6_all_step_temporal_admission_restart_calibration",
        "physical_equations": "unredefined_GR0_specialization_of_REF1", "amplitude": "3",
        "candidate_action_fields_or_health_stops_forbidden": True,
        "retained_EFT_interpretation": False,
    }
    expected_ladders = {
        "primary_method": "RK4", "primary_point_counts": [2049, 4097, 8193],
        "comparator_method": "SSPRK3", "comparator_point_counts": [4097, 8193, 16385],
        "common_physical_point_count": 2049,
        "each_method_must_pass_its_own_complete_admission": True,
        "raw_unequal_finest_grid_collocation_forbidden": True,
    }
    expected_restart = {
        "coordinate_time": "23/16",
        "primary_and_comparator_4097_8193_checkpoint": "runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz",
        "comparator_16385_checkpoint": "runs/fgc-2-sf1/rsp2/amplitude3-ssprk3-momentum/latest-checkpoint.npz",
        "complete_state_tracer_history_monitor_and_causal_payload_required": True,
        "continue_step_stage_transaction_source_retry_CFL_retry_and_causal_ledgers": True,
        "TDG6_ledger_starts_at_restart": True,
        "reinitialization_reprojection_interpolation_and_refit_forbidden": True,
        "restart_source_and_complete_common_event_prechecks_required": True,
    }
    expected_temporal = {
        "instrument_artifact_id": "FGC-1-TDG6-IMP2", "level_step_counts": [1, 2, 4],
        "channel_count": 18, "minimum_observed_order": "3/2",
        "exact_outward_order_test": "8*U12^2<=L01^2",
        "every_accepted_post_restart_macro_step_must_pass": True,
        "only_fine_path_may_commit": True, "temporal_retry_factor": "1/2",
        "maximum_temporal_retries_per_macro_step": 32,
        "minimum_macro_step": "1/1073741824",
        "temporal_rejections_durable_before_retry": True,
        "sampled_64_history_diagnostic_only": True,
    }
    expected_schedule = {
        "first_new_common_event_coordinate_time": "3/2", "final_coordinate_time": "32",
        "common_output_interval": "1/16", "checkpoint_interval": "1/4",
        "minimum_consecutive_trapped_common_events": 8,
        "trapped_sign_margin_over_combined_error_factor": 4,
        "temporal_history_minimum_samples": 64,
    }
    expected_assessment = {
        "minimum_constraint_finest_pair_order": "3/2", "maximum_nested_tail_ratio": "1/4",
        "trapped_scores_use_exact_common_physical_nodes": True,
        "each_method_supplies_own_finest_pair_Richardson_interval": True,
        "cross_method_Richardson_intervals_must_overlap": True,
        "constraint_and_spatial_admissions_are_hard_vetoes": True,
        "sampled_temporal_history_is_nonveto_diagnostic": True,
        "no_threshold_reduction_fitted_asymptote_or_post_outcome_change": True,
    }
    expected_provenance = {
        "output_root": "runs/fgc-2-sf1/proto14/calibration",
        "campaign_manifest": "runs/fgc-2-sf1/proto14/calibration/manifest.json",
        "append_only_event_log_name": "events.jsonl", "checkpoint_name": "latest-checkpoint.npz",
        "result_name": "campaign-result.json", "fresh_output_root_required_at_campaign_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_authorization_manifest_and_checkpoint": True,
        "authorization_commit_must_equal_manifest_git_commit": True,
        "terminal_checkpoint_must_precede_result_publication": True,
        "resume_must_idempotently_publish_or_verify_the_checkpoint_owned_terminal_result": True,
        "canonical_reduced_result_is_tracked_only_after_campaign_classification": True,
    }
    if (
        value["schema_version"] != 1 or value["artifact_id"] != "FGC-1-CAL11-RUN1-PLAN"
        or value["project_version"] != PROJECT_VERSION or value["metric_signature"] != "-+++"
        or value["riemann_convention"] != "plus_partial_mu_gamma_nu"
        or value["protocol_config"] != "configs/fgc/fgc-2-sf1-protocol-v14.toml"
        or value["protocol_freeze_result"] != "results/fgc-1-pro14-frz1.json"
        or value["runtime_authorization_result"] != "results/fgc-1-hlt12-mon12.json"
        or value["inherited_run_plan"] != "configs/fgc/fgc-1-pro13-run1.toml"
        or value["inherited_runtime_authorization"] != "results/fgc-1-hlt11-mon11.json"
        or value["tdg6_runtime_result"] != "results/fgc-1-tdg6-imp2.json"
        or value["scope"] != expected_scope or value["numerical_runtime"] != EXPECTED_PLAN_RUNTIME
        or value["method_owned_ladders"] != expected_ladders or value["restart"] != expected_restart
        or value["temporal_admission"] != expected_temporal or value["schedule"] != expected_schedule
        or value["assessment"] != expected_assessment or value["provenance"] != expected_provenance
        or value["claims"] != EXPECTED_RUN_PLAN_CLAIMS
    ):
        raise ValueError("PROTO14 run-plan violates its frozen scope")
    return value


def numerical_runtime_environment() -> dict[str, Any]:
    """Observe the exact arithmetic backend before source or state work."""

    configuration = np.show_config(mode="dicts")
    if not isinstance(configuration, dict):
        raise ValueError("NumPy build configuration is unavailable")
    dependencies = configuration.get("Build Dependencies", {})
    blas = dependencies.get("blas", {}) if isinstance(dependencies, dict) else {}
    lapack = dependencies.get("lapack", {}) if isinstance(dependencies, dict) else {}
    if not isinstance(blas, dict) or not isinstance(lapack, dict):
        raise ValueError("NumPy BLAS/LAPACK configuration is unavailable")
    return {
        "python_executable": str(Path(sys.executable).resolve()),
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "python_compiler": platform.python_compiler(),
        "python_build": list(platform.python_build()),
        "numpy_version": np.__version__,
        "blas": dict(blas),
        "lapack": dict(lapack),
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "version": platform.version(),
            "machine": platform.machine(),
            "processor": platform.processor(),
        },
        "numpy_build_configuration": configuration,
    }


def validate_numerical_runtime(contract: Mapping[str, Any]) -> dict[str, Any]:
    if dict(contract) != EXPECTED_RUNTIME:
        raise ValueError("HLT12 numerical runtime contract differs")
    observed = numerical_runtime_environment()
    actual = {
        "python_version": observed["python_version"],
        "python_implementation": observed["python_implementation"],
        "numpy_version": observed["numpy_version"],
        "blas_name": observed["blas"].get("name"),
        "blas_version": observed["blas"].get("version"),
        "system": observed["platform"]["system"],
        "machine": observed["platform"]["machine"],
        "full_runtime_metadata_must_be_recorded_in_manifest_and_result": True,
        "runtime_drift_must_fail_before_source_or_state_advance": True,
    }
    if actual != EXPECTED_RUNTIME:
        raise ValueError("HLT12 numerical runtime differs from its frozen contract")
    return observed


def _lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("HLT12 checkpoint is not an ancestor of HEAD")
    pairs = (
        (config["protocol_config"], lineage["protocol_config_sha256"]),
        (config["protocol_freeze_config"], lineage["protocol_freeze_config_sha256"]),
        (config["protocol_freeze_result"], lineage["protocol_freeze_result_sha256"]),
        (lineage["protocol_module"], lineage["protocol_module_sha256"]),
        (lineage["TDG6_runtime_module"], lineage["TDG6_runtime_module_sha256"]),
    )
    tracked = []
    for relative, expected in pairs:
        blob = _git("show", f"{commit}:{relative}").stdout
        if sha256(blob).hexdigest() != expected or _sha(REPOSITORY / relative) != expected:
            raise ValueError(f"HLT12 tracked lineage differs: {relative}")
        tracked.append({"path": relative, "sha256": expected})
    raw = {}
    for name in ("PROTO12", "RSP2"):
        relative = lineage[f"{name}_checkpoint"]
        expected = lineage[f"{name}_checkpoint_sha256"]
        if not (REPOSITORY / relative).is_file() or _sha(REPOSITORY / relative) != expected:
            raise ValueError(f"HLT12 {name} checkpoint differs")
        raw[relative] = expected
    run_plan = REPOSITORY / config["run_plan_config"]
    expected_plan_sha = lineage["run_plan_sha256"]
    if not run_plan.is_file() or _sha(run_plan) != expected_plan_sha:
        raise ValueError("HLT12 run plan differs from its hash-bound contract")
    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "tracked_blobs": tracked,
        "raw_restart_checkpoints": raw,
        "run_plan": {"path": config["run_plan_config"], "sha256": expected_plan_sha, "created_after_predecessor_checkpoint": True},
    }


def _restart_payloads(config: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Use PRO14's hash checks while guarding all raw inputs from mutation."""

    lineage = config["immutable_lineage"]
    checkpoints = [REPOSITORY / lineage[f"{name}_checkpoint"] for name in ("PROTO12", "RSP2")]
    before = {_rel(path): _sha(path) for path in checkpoints}
    # PRO14-FRZ1 is the owner of the historical PROTO13 manifest, so reuse
    # its exact restoration contract rather than constructing a second one.
    records = pro14._restore_members(
        pro14.load_config(REPOSITORY / config["protocol_freeze_config"])
    )  # type: ignore[attr-defined]
    after = {_rel(path): _sha(path) for path in checkpoints}
    if before != after or len(records) != 6:
        raise ValueError("HLT12 restart restoration mutated or failed to restore payloads")
    expected_ladders = ([2049, 4097, 8193], [4097, 8193, 16385])
    if ([item["point_count"] for item in records[:3]], [item["point_count"] for item in records[3:]]) != expected_ladders:
        raise ValueError("HLT12 restart ladder differs")
    return [
        {
            "key": item["key"], "source_checkpoint": item["source_checkpoint"],
            "method": item["method"], "point_count": item["point_count"],
            "coordinate_time": item["coordinate_time"], "state_sha256": item["state_sha256"],
            "restart_payload_sha256": item["restart_payload_sha256"], "input_hash": item["input_hash"],
            "step_index": item["step_index"], "transaction_serial": item["transaction_serial"],
            "accepted_stage_count": item["accepted_stage_count"], "source_retry_count": item["source_retry_count"],
            "CFL_retry_count": item["CFL_retry_count"], "event_sample_count": item["event_sample_count"],
            "tracer_count": item["tracer_count"], "state_shape": item["state_shape"],
            "TDG6_ledger_initialization": {"accepted_macro_step_count": 0, "cumulative_temporal_retry_count": 0, "serialized_temporal_rejection_count": 0, "accumulated_debit_channel_count": 18, "all_debits_exact_zero": True},
        }
        for item in records
    ]


def _runtime_implementation(config: Mapping[str, Any]) -> dict[str, str]:
    bindings = config["runtime_binding"]
    paths = [
        REPOSITORY / config["immutable_lineage"]["TDG6_runtime_module"],
        REPOSITORY / bindings["successor_runtime_module"],
        REPOSITORY / bindings["authorized_runner"],
        Path(__file__).resolve(),
    ]
    result = {}
    for path in paths:
        if not path.is_file():
            raise ValueError(f"HLT12 implementation is absent: {_rel(path)}")
        result[_rel(path)] = _sha(path)
    runner_source = (REPOSITORY / bindings["authorized_runner"]).read_text(encoding="utf-8")
    if f'RUNNER_ID = "{bindings["authorized_runner_id"]}"' not in runner_source:
        raise ValueError("HLT12 authorized runner identity differs")
    return result


def _integration_controls() -> dict[str, Any]:
    """Execute bounded synthetic adapter controls without opening a namespace."""

    sources = (
        REPOSITORY / "tests/test_fgc_proto14_runtime.py",
        REPOSITORY / "tests/test_fgc_tdg6_temporal_admission_runtime.py",
        REPOSITORY / "tests/test_fgc_gr0_campaign_runner_v14.py",
    )
    for path in sources:
        if not path.is_file():
            raise ValueError(f"HLT12 integration control source is absent: {_rel(path)}")
    evidence = {}
    for name, test_id in INTEGRATION_TESTS:
        completed = subprocess.run(
            [sys.executable, "-m", "unittest", test_id],
            cwd=REPOSITORY,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"HLT12 integration control failed: {name}: {completed.stderr.strip()}"
            )
        evidence[name] = True
    evidence["sampled_64_history_is_nonveto_diagnostic_only"] = True
    evidence["production_namespace_created"] = False
    evidence["sources_sha256"] = {_rel(path): _sha(path) for path in sources}
    return evidence


def _namespace(config: Mapping[str, Any], historical: Mapping[str, Any] | None) -> dict[str, Any]:
    if historical is not None:
        records = historical.get("records")
        if not (
            historical.get("both_new_namespaces_absent") is True
            and historical.get("authorization_created_no_namespace") is True
            and isinstance(records, list) and len(records) == 2
        ):
            raise ValueError("HLT12 historical namespace evidence differs")
        return deepcopy(dict(historical))
    records = []
    for name in ("calibration_output_root", "holdout_output_root"):
        relative = config["namespace"][name]
        if (REPOSITORY / relative).exists():
            raise ValueError(f"HLT12 prospective namespace exists: {relative}")
        records.append({"path": relative, "absent_before_authorization": True})
    return {"both_new_namespaces_absent": True, "authorization_created_no_namespace": True, "records": records}


def _mutation_controls(config: Mapping[str, Any]) -> dict[str, bool]:
    controls: dict[str, bool] = {}
    config_source = DEFAULT_CONFIG.read_text(encoding="utf-8")
    attacks = {
        "threshold": ('restart_coordinate_time = "23/16"', 'restart_coordinate_time = "3/2"'),
        "runtime": ('python_version = "3.14.3"', 'python_version = "3.14.6"'),
        "claim": (
            "FGCQR_holdout_execution_authorized = false",
            "FGCQR_holdout_execution_authorized = true",
        ),
        "namespace": (
            'calibration_output_root = "runs/fgc-2-sf1/proto14/calibration"',
            'calibration_output_root = "runs/fgc-2-sf1/proto13/calibration"',
        ),
        "runner": (
            'authorized_runner = "scripts/run_fgc_gr0_calibration_v14.py"',
            'authorized_runner = "scripts/run_fgc_gr0_calibration_v13.py"',
        ),
    }
    for name, (before, after) in attacks.items():
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "attacked.toml"
            target.write_text(config_source.replace(before, after, 1), encoding="utf-8")
            try:
                load_config(target)
            except ValueError:
                controls[name] = True
            else:
                controls[name] = False
    plan_source = (REPOSITORY / config["run_plan_config"]).read_text(encoding="utf-8")
    plan_attacks = {
        "plan_threshold": (
            'minimum_observed_order = "3/2"',
            'minimum_observed_order = "149/100"',
        ),
        "plan_temporal_veto": (
            "sampled_64_history_diagnostic_only = true",
            "sampled_64_history_diagnostic_only = false",
        ),
        "plan_candidate_scope": (
            "candidate_action_fields_or_health_stops_forbidden = true",
            "candidate_action_fields_or_health_stops_forbidden = false",
        ),
    }
    for name, (before, after) in plan_attacks.items():
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "attacked.toml"
            target.write_text(plan_source.replace(before, after, 1), encoding="utf-8")
            try:
                validate_run_plan(target)
            except ValueError:
                controls[name] = True
            else:
                controls[name] = False
    if not all(controls.values()):
        raise RuntimeError("HLT12 mutation control failed")
    return controls


def record(config_path: Path = DEFAULT_CONFIG, *, historical_namespace_evidence: Mapping[str, Any] | None = None) -> dict[str, Any]:
    config = load_config(config_path)
    plan = validate_run_plan(REPOSITORY / config["run_plan_config"])
    runtime = validate_numerical_runtime(config["numerical_runtime"])
    if plan["numerical_runtime"] != EXPECTED_PLAN_RUNTIME:
        raise ValueError("HLT12 run plan runtime differs from authorization runtime")
    pro14.verify_canonical(REPOSITORY / config["protocol_freeze_result"], REPOSITORY / config["protocol_freeze_config"])
    lineage = _lineage(config)
    restart_payloads = _restart_payloads(config)
    implementation = _runtime_implementation(config)
    integration_controls = _integration_controls()
    namespace = _namespace(config, historical_namespace_evidence)
    controls = _mutation_controls(config)
    campaign_seed = {
        "scope": config["scope"], "restart_payloads": [item["restart_payload_sha256"] for item in restart_payloads],
        "implementation_sha256": implementation, "runtime_contract": config["numerical_runtime"],
        "run_plan_sha256": _sha(REPOSITORY / config["run_plan_config"]),
    }
    gates = dict(config["claims"])
    return {
        "schema_version": 1, "artifact_id": ARTIFACT_ID, "project_version": PROJECT_VERSION,
        "classification": "pretrajectory_PROTO14_TDG6_restart_runtime_authorization",
        "generated_by": "scripts/reproduce_fgc_hlt12_mon12.py",
        "derivation_document": _rel(OWNER_DOCUMENT), "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): _sha(config_path), config["protocol_config"]: _sha(REPOSITORY / config["protocol_config"]), config["protocol_freeze_config"]: _sha(REPOSITORY / config["protocol_freeze_config"]), config["run_plan_config"]: _sha(REPOSITORY / config["run_plan_config"])},
        "predecessor_sha256": {config["protocol_freeze_result"]: config["immutable_lineage"]["protocol_freeze_result_sha256"], config["run_plan_config"]: config["immutable_lineage"]["run_plan_sha256"], config["immutable_lineage"]["PROTO12_checkpoint"]: config["immutable_lineage"]["PROTO12_checkpoint_sha256"], config["immutable_lineage"]["RSP2_checkpoint"]: config["immutable_lineage"]["RSP2_checkpoint_sha256"]},
        "implementation_sha256": implementation, "scope_bindings": dict(config["scope"]),
        "artifact_payload": {
            "immutable_lineage": lineage,
            "frozen_run_plan": {"campaign_id": sha256(_canonical(campaign_seed).encode("utf-8")).hexdigest(), "run_plan_sha256": _sha(REPOSITORY / config["run_plan_config"]), "numerical_runtime_contract": dict(plan["numerical_runtime"]), "restart_coordinate_time": 1.4375, "first_new_common_event_coordinate_time": 1.5, "final_coordinate_time": 32.0, "primary_point_counts": [2049, 4097, 8193], "comparator_point_counts": [4097, 8193, 16385], "TDG6_every_accepted_macro_step_required": True, "sampled_64_history_is_nonveto_diagnostic_only": True},
            "runtime_environment_observed": runtime,
            "restart_payloads": restart_payloads,
            "executed_integration_controls": integration_controls,
            "namespace_precondition": namespace,
            "mutation_controls": controls,
            "decision": {"PROTO14_successor_runtime_implemented": True, "PROTO14_fresh_GR0_dynamic_calibration_authorized": True, "candidate_execution_authorized": False, "retained_EFT_evolution_authorized": False},
            "epistemic_boundary": {"trajectory_advanced": False, "fresh_GR0_dynamic_calibration_completed": False, "eligible_GR0_case_selected": False, "candidate_or_mechanism_outcome_read": False, "TDG6_debit_is_global_PDE_or_Raychaudhuri_error_bound": False},
        },
        "gate_status": gates,
        "nonclaims": {key: value for key, value in gates.items() if value is False},
    }


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    value = json.loads(source, object_pairs_hook=_reject_duplicate_pairs)
    if not isinstance(value, dict) or value.get("artifact_id") != ARTIFACT_ID:
        raise ValueError("HLT12 canonical result identity differs")
    if source != _canonical(value):
        raise ValueError("HLT12 result is not sorted canonical JSON")
    return value


def verify_canonical(path: Path = DEFAULT_OUTPUT, config_path: Path = DEFAULT_CONFIG) -> None:
    observed = load_canonical_result(path)
    expected = record(config_path, historical_namespace_evidence=observed["artifact_payload"]["namespace_precondition"])
    if observed != expected:
        raise ValueError("HLT12 output differs from a fresh reproduction")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        verify_canonical(args.output, args.config)
        print(f"verified {args.output}")
        return
    _atomic_write(args.output, _canonical(record(args.config)).encode("utf-8"))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
