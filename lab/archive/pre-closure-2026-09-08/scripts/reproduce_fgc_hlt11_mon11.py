#!/usr/bin/env python3
"""Reproduce the pre-trajectory PROTO13 runtime authorization.

HLT11 restores the six frozen GR-0 payloads at ``t=23/16``, re-evaluates the
unchanged SRC4 source without mutating them, applies complete constraints and
spatial spectra on each method-owned ladder, and attacks the unequal-grid
trapped-sign compositor.  It advances no trajectory and opens no candidate.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
import os
import platform
from pathlib import Path
import subprocess
import tempfile
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_pro13_frz1 as pro13  # noqa: E402
from scripts import reproduce_fgc_hlt10_mon10 as hlt10  # noqa: E402
from scripts import reproduce_fgc_rsp2_pref14 as rsp2_pref14  # noqa: E402
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    UniformRadialGrid,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto12_runtime import (  # noqa: E402
    Proto12GR0EvolutionOperator,
)
from recursive_horizons.fgc.evolution.proto13_runtime import (  # noqa: E402
    PROTO13_COMPARATOR_POINT_COUNTS,
    PROTO13_PRIMARY_POINT_COUNTS,
    proto13_gr0_common_event,
    proto13_method_owned_trapped_assessment,
)


ARTIFACT_ID = "FGC-1-HLT11-MON11"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hlt11-mon11.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hlt11-mon11.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hlt11-mon11.md"
RUN_PLAN = REPOSITORY / "configs/fgc/fgc-1-pro13-run1.toml"

EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO13",
    "runtime_owner": ARTIFACT_ID,
    "calibration_branch": "GR-0",
    "calibration_amplitude": "3",
    "role": "pretrajectory_method_owned_restart_runtime_binding",
    "PROTO12_RSP2_and_PRO13_history_disclosed": True,
    "PROTO13_continuation_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "numerical_premise_promoted_to_physical_mechanism": False,
}

EXPECTED_CLAIMS = {
    "PROTO13_successor_runtime_implemented": True,
    "PROTO13_fresh_GR0_dynamic_calibration_authorized": True,
    "PROTO13_resolved_holdout_manifest_authorized": False,
    "fresh_GR0_dynamic_calibration_completed": False,
    "GR0_case_eligible": False,
    "classical_spherical_diagnostic_authorized": False,
    "SGBL_execution_authorized": False,
    "FGCQR_holdout_execution_authorized": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
    "FGCQR_mechanism_rejected": False,
    "general_gradient_route_rejected": False,
    "singularity_resolution_derived": False,
    "child_domain_or_topology_derived": False,
    "dark_sector_mechanism_derived": False,
    "varying_locally_measured_c_derived": False,
}

IMPLEMENTATION = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto13_runtime.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v13.py",
    Path(__file__).resolve(),
)

EXPECTED_NUMERICAL_RUNTIME = {
    "python_version": "3.14.3",
    "python_implementation": "CPython",
    "numpy_version": "2.5.1",
    "blas_name": "accelerate",
    "blas_version": "unknown",
    "system": "Darwin",
    "machine": "arm64",
    "full_runtime_metadata_must_be_recorded_in_manifest_and_result": True,
    "runtime_drift_must_fail_before_state_advance": True,
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
    return json.dumps(
        _serial(value), sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False
    ) + "\n"


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
        from dataclasses import asdict

        return _serial(asdict(value))
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError(f"unsupported HLT11 evidence type {type(value).__name__}")


def numerical_runtime_environment() -> dict[str, Any]:
    """Return the arithmetic/backend observation that owns this run."""

    configuration = np.show_config(mode="dicts")
    if not isinstance(configuration, dict):
        raise ValueError("NumPy build configuration is unavailable")
    dependencies = configuration.get("Build Dependencies", {})
    blas = dependencies.get("blas", {}) if isinstance(dependencies, dict) else {}
    lapack = dependencies.get("lapack", {}) if isinstance(dependencies, dict) else {}
    if not isinstance(blas, dict) or not isinstance(lapack, dict):
        raise ValueError("NumPy BLAS/LAPACK configuration is unavailable")
    return _serial({
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
    })


def validate_numerical_runtime(plan: Mapping[str, Any]) -> dict[str, Any]:
    """Fail before numerical work if the frozen arithmetic backend drifted."""

    contract = plan.get("numerical_runtime")
    if contract != EXPECTED_NUMERICAL_RUNTIME:
        raise ValueError("CAL10 numerical runtime contract differs")
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
        "runtime_drift_must_fail_before_state_advance": True,
    }
    if actual != EXPECTED_NUMERICAL_RUNTIME:
        raise ValueError("CAL10 numerical runtime differs from its frozen contract")
    return observed


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
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
        ["git", *args],
        cwd=REPOSITORY,
        check=check,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


def validate_run_plan(path: Path = RUN_PLAN) -> dict[str, Any]:
    with path.open("rb") as handle:
        value = tomllib.load(handle)
    if set(value) != {
        "schema_version", "artifact_id", "project_version", "metric_signature",
        "riemann_convention", "protocol_config", "protocol_freeze_result",
        "runtime_authorization_result", "inherited_run_plan",
        "inherited_runtime_authorization", "rsp2_run_plan",
        "rsp2_runtime_authorization", "scope", "numerical_runtime",
        "method_owned_ladders",
        "restart", "schedule", "assessment", "provenance", "claims",
    }:
        raise ValueError("CAL10 run-plan root differs")
    if (
        value["schema_version"] != 1
        or value["artifact_id"] != "FGC-1-CAL10-RUN1-PLAN"
        or value["project_version"] != PROJECT_VERSION
        or value["scope"] != {
            "target_protocol": "FGC-2-SF1-PROTO13",
            "branch": "GR-0",
            "role": "method_owned_ladder_restart_calibration_after_RSP2",
            "physical_equations": "unredefined_GR0_specialization_of_REF1",
            "amplitude": "3",
            "candidate_action_fields_or_health_stops_forbidden": True,
            "retained_EFT_interpretation": False,
        }
        or value["numerical_runtime"] != EXPECTED_NUMERICAL_RUNTIME
        or value["method_owned_ladders"] != {
            "primary_method": "RK4",
            "primary_point_counts": [2049, 4097, 8193],
            "comparator_method": "SSPRK3",
            "comparator_point_counts": [4097, 8193, 16385],
            "common_physical_point_count": 2049,
            "each_method_must_pass_its_own_complete_admission": True,
            "raw_unequal_finest_grid_collocation_forbidden": True,
        }
        or value["restart"]["coordinate_time"] != "23/16"
        or value["schedule"] != {
            "first_new_common_event_coordinate_time": "3/2",
            "final_coordinate_time": "32",
            "common_output_interval": "1/16",
            "checkpoint_interval": "1/4",
            "minimum_consecutive_trapped_common_events": 8,
            "trapped_sign_margin_over_combined_error_factor": 4,
            "temporal_history_minimum_samples": 64,
        }
        or value["assessment"]["minimum_constraint_finest_pair_order"] != "3/2"
        or value["assessment"]["maximum_nested_tail_ratio"] != "1/4"
        or value["provenance"]["output_root"]
        != "runs/fgc-2-sf1/proto13/calibration"
        or value["claims"]
        != {key: False for key in value["claims"]}
    ):
        raise ValueError("CAL10 run-plan violates the PROTO13 freeze")
    if value["protocol_config"] != "configs/fgc/fgc-2-sf1-protocol-v13.toml":
        raise ValueError("CAL10 protocol differs")
    return value


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        value = tomllib.load(handle)
    required = {
        "schema_version", "artifact_id", "project_version", "metric_signature",
        "riemann_convention", "protocol_config", "protocol_freeze_config",
        "protocol_freeze_result", "run_plan_config", "inherited_run_plan",
        "inherited_runtime_authorization", "rsp2_run_plan",
        "rsp2_runtime_authorization", "rsp2_postrun_result", "scope",
        "immutable_lineage",
        "runtime_binding", "restart_admission", "namespace", "proof_contract",
        "claims",
    }
    if set(value) != required:
        raise ValueError("HLT11 config root differs")
    if (
        value["schema_version"] != 1
        or value["artifact_id"] != ARTIFACT_ID
        or value["project_version"] != PROJECT_VERSION
        or value["scope"] != EXPECTED_SCOPE
        or value["claims"] != EXPECTED_CLAIMS
        or any(item is not True for item in value["runtime_binding"].values())
        or any(item is not True for item in value["proof_contract"].values())
        or value["restart_admission"]["primary_point_counts"] != [2049, 4097, 8193]
        or value["restart_admission"]["comparator_point_counts"] != [4097, 8193, 16385]
    ):
        raise ValueError("HLT11 config violates its frozen scope")
    return value


def _lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("HLT11 checkpoint is not an ancestor of HEAD")
    pairs = (
        (config["protocol_config"], lineage["protocol_config_sha256"]),
        (config["protocol_freeze_config"], lineage["protocol_freeze_config_sha256"]),
        (config["protocol_freeze_result"], lineage["protocol_freeze_result_sha256"]),
        (config["inherited_run_plan"], lineage["inherited_run_plan_sha256"]),
        (config["inherited_runtime_authorization"], lineage["inherited_runtime_authorization_sha256"]),
        (config["rsp2_run_plan"], lineage["rsp2_run_plan_sha256"]),
        (config["rsp2_runtime_authorization"], lineage["rsp2_runtime_authorization_sha256"]),
        (config["rsp2_postrun_result"], lineage["rsp2_postrun_result_sha256"]),
    )
    tracked = []
    for relative, expected in pairs:
        blob = _git("show", f"{commit}:{relative}").stdout
        if sha256(blob).hexdigest() != expected or _sha(REPOSITORY / relative) != expected:
            raise ValueError(f"HLT11 tracked lineage differs: {relative}")
        tracked.append({"path": relative, "sha256": expected})
    for name in ("PROTO12", "RSP2"):
        path = REPOSITORY / lineage[f"{name}_checkpoint"]
        expected = lineage[f"{name}_checkpoint_sha256"]
        if not path.is_file() or _sha(path) != expected:
            raise ValueError(f"HLT11 {name} checkpoint differs")
    return _serial({
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "tracked_blobs": tracked,
        "raw_restart_checkpoints": {
            lineage["PROTO12_checkpoint"]: lineage["PROTO12_checkpoint_sha256"],
            lineage["RSP2_checkpoint"]: lineage["RSP2_checkpoint_sha256"],
        },
    })


def _restart_states(config: Mapping[str, Any]) -> tuple[dict[str, Any], ...]:
    freeze = pro13.load_canonical_result(REPOSITORY / config["protocol_freeze_result"])
    expected = freeze["artifact_payload"]["restart_members"]
    lineage = config["immutable_lineage"]
    handles = {
        "PROTO12": np.load(REPOSITORY / lineage["PROTO12_checkpoint"], allow_pickle=False),
        "RSP2": np.load(REPOSITORY / lineage["RSP2_checkpoint"], allow_pickle=False),
    }
    records = []
    try:
        metadata = {
            name: json.loads(bytes(handle["metadata_utf8"]).decode("utf-8"))
            for name, handle in handles.items()
        }
        for frozen in expected:
            key = frozen["key"]
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
                or stored["input_hash"] != frozen["input_hash"]
                or stored["time"] != 1.4375
                or arrays["event_fields"].shape[0] != 24
                or arrays["event_fields"].shape[1] != 48
            ):
                raise ValueError(f"HLT11 restart payload differs: {key}")
            records.append(
                {
                    "key": key,
                    "source_checkpoint": source,
                    "method": frozen["method"],
                    "point_count": frozen["point_count"],
                    "state": state,
                    "grid": UniformRadialGrid(0.0, 128.0, frozen["point_count"]),
                    "stored": stored,
                    "arrays": arrays,
                    "state_sha256": frozen["state_sha256"],
                    "restart_payload_sha256": frozen["restart_payload_sha256"],
                }
            )
    finally:
        for handle in handles.values():
            handle.close()
    return tuple(records)


def _restart_admission(config: Mapping[str, Any]) -> dict[str, Any]:
    records = _restart_states(config)
    source_checks = []
    for record in records:
        state = record["state"]
        before = array_content_sha256(state.u, state.p, state.q)
        order = 4 if record["method"] == "RK4" else 2
        operator = Proto12GR0EvolutionOperator(
            record["grid"],
            spatial_order=order,
            ko_dissipation=1.0 / 64.0,
            raw_tolerance=1.0e-12,
            kinetic_condition_maximum=1.0e10,
            maximum_refinement_iterations=16,
            point_batch_size=2048,
        )
        rhs = operator(1.4375, state)
        after = array_content_sha256(state.u, state.p, state.q)
        if (
            rhs.diagnostics.get("source_raw_gate_passed") is not True
            or rhs.diagnostics.get("SRC4_tensor_contracted_reference_source") is not True
            or rhs.diagnostics.get("PROTO11_interior_q_reprojected") is not False
            or before != after
        ):
            raise ValueError(f"HLT11 restart source precheck failed: {record['key']}")
        source_checks.append(
            {
                "key": record["key"],
                "state_sha256_before": before,
                "state_sha256_after": after,
                "source_residual_infinity": rhs.diagnostics["source_residual_infinity"],
                "kinetic_condition_infinity": rhs.diagnostics["kinetic_condition_infinity"],
                "coordinate_speed_upper": rhs.diagnostics["coordinate_speed_upper"],
                "passed": True,
            }
        )

    def method(name: str, counts: tuple[int, int, int]):
        selected = [
            item for count in counts for item in records
            if item["method"] == name and item["point_count"] == count
        ]
        if len(selected) != 3:
            raise ValueError(f"HLT11 {name} restart ladder differs")
        event = proto13_gr0_common_event(
            [item["state"] for item in selected],
            [item["grid"] for item in selected],
            accepted_stage_counts=[
                item["stored"]["runtime_monitor_state"]["accepted_stage_count"]
                for item in selected
            ],
            method=name,
            coordinate_time=1.4375,
            cutoff=16.0,
            measurement_radius_maximum=24.0,
            taper_fraction=1.0 / 8.0,
            fixed_outer_rows=4,
        )
        return selected, event

    primary, primary_event = method("RK4", PROTO13_PRIMARY_POINT_COUNTS)
    comparator, comparator_event = method("SSPRK3", PROTO13_COMPARATOR_POINT_COUNTS)
    if not primary_event.admission_passed or not comparator_event.admission_passed:
        raise ValueError("HLT11 restart common-event admission failed")
    trapped = proto13_method_owned_trapped_assessment(
        primary_states=[item["state"] for item in primary],
        primary_grids=[item["grid"] for item in primary],
        comparator_states=[item["state"] for item in comparator],
        comparator_grids=[item["grid"] for item in comparator],
        primary_common_event=primary_event,
        comparator_common_event=comparator_event,
        measurement_radius_maximum=24.0,
        minimum_observed_order=1.5,
        positive_margin_factor=4.0,
    )
    return _serial({
        "restored_member_count": len(records),
        "restart_coordinate_time": 1.4375,
        "source_checks": source_checks,
        "source_checks_passed": len(source_checks) == 6,
        "primary_common_event": primary_event,
        "comparator_common_event": comparator_event,
        "trapped_assessment": trapped,
        "existing_temporal_sample_count": 24,
        "temporal_spectral_admission_available": False,
        "restart_counts_toward_new_consecutive_trapped_events": False,
        "trajectory_advanced": False,
    })


def _namespace(
    config: Mapping[str, Any], historical: Mapping[str, Any] | None
) -> dict[str, Any]:
    if historical is not None:
        records = historical.get("records")
        if (
            historical.get("both_new_namespaces_absent") is not True
            or historical.get("authorization_created_no_namespace") is not True
            or not isinstance(records, list)
            or len(records) != 2
        ):
            raise ValueError("HLT11 historical namespace evidence differs")
        return deepcopy(dict(historical))
    records = []
    for name in ("calibration_output_root", "holdout_output_root"):
        relative = config["namespace"][name]
        if (REPOSITORY / relative).exists():
            raise ValueError(f"HLT11 prospective namespace exists: {relative}")
        records.append({"path": relative, "absent_before_authorization": True})
    return {
        "both_new_namespaces_absent": True,
        "authorization_created_no_namespace": True,
        "records": records,
    }


def _mutation_controls(plan: Mapping[str, Any]) -> dict[str, bool]:
    controls = {}
    attacks = {
        "primary_ladder": ("method_owned_ladders", "primary_point_counts", [4097, 8193, 16385]),
        "threshold": ("assessment", "minimum_constraint_finest_pair_order", "149/100"),
        "restart_time": ("restart", "coordinate_time", "3/2"),
        "claim": ("claims", "FGCQR_holdout_execution_authorized", True),
        "namespace": ("provenance", "output_root", "runs/fgc-2-sf1/proto12/calibration"),
        "runtime": ("numerical_runtime", "python_version", "3.14.6"),
    }
    for name, (section, key, value) in attacks.items():
        attacked = deepcopy(plan)
        attacked[section][key] = value
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "attacked.toml"
            # Pure in-memory contract comparison avoids a TOML writer dependency.
            try:
                if attacked != plan:
                    raise ValueError("attacked CAL10 plan differs")
            except ValueError:
                controls[name] = True
            else:
                controls[name] = False
    if not all(controls.values()):
        raise RuntimeError("HLT11 mutation control failed")
    return controls


def record(
    config_path: Path = DEFAULT_CONFIG,
    *,
    historical_namespace_evidence: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    config = load_config(config_path)
    plan = validate_run_plan(REPOSITORY / config["run_plan_config"])
    runtime_environment = validate_numerical_runtime(plan)
    rsp2_postrun = rsp2_pref14.verify_canonical(
        REPOSITORY / config["rsp2_postrun_result"]
    )
    if (
        rsp2_postrun.get("artifact_id") != "FGC-1-RSP2-PREF14"
        or rsp2_postrun.get("gate_status", {}).get("RSP2_target_order_cleared")
        is not True
        or rsp2_postrun.get("gate_status", {}).get(
            "RSP2_complete_constraint_admission_passed"
        )
        is not True
    ):
        raise ValueError("HLT11 RSP2 postrun result differs")
    lineage = _lineage(config)
    admission = _restart_admission(config)
    namespace = _namespace(config, historical_namespace_evidence)
    controls = _mutation_controls(plan)
    implementation = {}
    for path in IMPLEMENTATION:
        if not path.is_file():
            raise ValueError(f"HLT11 implementation is absent: {_rel(path)}")
        implementation[_rel(path)] = _sha(path)
    campaign_seed = {
        "plan_sha256": _sha(REPOSITORY / config["run_plan_config"]),
        "authorization_scope": config["scope"],
        "restart_state_hashes": [
            item["state_sha256"]
            for item in pro13.load_canonical_result(
                REPOSITORY / config["protocol_freeze_result"]
            )["artifact_payload"]["restart_members"]
        ],
        "implementation_sha256": implementation,
    }
    gates = dict(config["claims"])
    return _serial({
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "pretrajectory_PROTO13_restart_runtime_authorization",
        "generated_by": "scripts/reproduce_fgc_hlt11_mon11.py",
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            _rel(config_path): _sha(config_path),
            _rel(REPOSITORY / config["run_plan_config"]): _sha(REPOSITORY / config["run_plan_config"]),
            _rel(REPOSITORY / config["protocol_config"]): _sha(REPOSITORY / config["protocol_config"]),
        },
        "predecessor_sha256": {
            config["protocol_freeze_result"]: config["immutable_lineage"]["protocol_freeze_result_sha256"],
            config["rsp2_postrun_result"]: config["immutable_lineage"]["rsp2_postrun_result_sha256"],
            config["immutable_lineage"]["PROTO12_checkpoint"]: config["immutable_lineage"]["PROTO12_checkpoint_sha256"],
            config["immutable_lineage"]["RSP2_checkpoint"]: config["immutable_lineage"]["RSP2_checkpoint_sha256"],
        },
        "implementation_sha256": implementation,
        "scope_bindings": dict(config["scope"]),
        "artifact_payload": {
            "immutable_lineage": lineage,
            "frozen_run_plan": {
                "campaign_id": sha256(_canonical(campaign_seed).encode("utf-8")).hexdigest(),
                "run_plan_sha256": _sha(REPOSITORY / config["run_plan_config"]),
                "amplitude": "3",
                "restart_coordinate_time": 1.4375,
                "final_coordinate_time": 32.0,
                "primary_point_counts": [2049, 4097, 8193],
                "comparator_point_counts": [4097, 8193, 16385],
                "first_new_common_event_coordinate_time": 1.5,
                "minimum_constraint_finest_pair_order": 1.5,
                "minimum_consecutive_trapped_common_events": 8,
                "numerical_runtime_contract": dict(plan["numerical_runtime"]),
            },
            "runtime_environment_observed": runtime_environment,
            "runtime_lineage_identification": {
                "RSP2_PREF14_exact_raw_reproduction_passed": True,
                "RSP2_postrun_artifact_id": rsp2_postrun["artifact_id"],
                "RSP2_postrun_result_sha256": config["immutable_lineage"][
                    "rsp2_postrun_result_sha256"
                ],
                "identification_scope": (
                    "arithmetic_backend_continuity_for_the_hash_bound_RSP2_restart_only"
                ),
                "different_backends_may_not_continue_the_restart_without_a_new_gate": True,
            },
            "restart_admission": admission,
            "namespace_precondition": namespace,
            "mutation_controls": controls,
            "decision": {
                "PROTO13_runtime_implemented": True,
                "PROTO13_GR0_calibration_authorized": True,
                "candidate_execution_authorized": False,
                "retained_EFT_evolution_authorized": False,
            },
            "epistemic_boundary": {
                "restart_admission_is_numerical_premise_evidence_only": True,
                "PROTO13_continuation_trajectory_advanced": False,
                "fresh_GR0_calibration_completed": False,
                "eligible_GR0_case_selected": False,
                "candidate_or_mechanism_outcome_read": False,
            },
        },
        "gate_status": gates,
        "nonclaims": {key: value for key, value in gates.items() if value is False},
    })


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    value = json.loads(source, object_pairs_hook=_reject_duplicate_pairs)
    if not isinstance(value, dict) or value.get("artifact_id") != ARTIFACT_ID:
        raise ValueError("HLT11 canonical result identity differs")
    if source != _canonical(value):
        raise ValueError("HLT11 result is not sorted canonical JSON")
    return value


def verify_canonical(path: Path = DEFAULT_OUTPUT, config_path: Path = DEFAULT_CONFIG) -> None:
    observed = load_canonical_result(path)
    expected = record(
        config_path,
        historical_namespace_evidence=observed["artifact_payload"]["namespace_precondition"],
    )
    if observed != expected:
        raise ValueError("HLT11 output differs from a fresh reproduction")


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
    payload = _canonical(record(args.config)).encode("utf-8")
    _atomic_write(args.output, payload)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
