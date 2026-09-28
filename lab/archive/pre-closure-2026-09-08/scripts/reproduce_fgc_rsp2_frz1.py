#!/usr/bin/env python3
"""Reproduce the prospective FGC-1-RSP2-FRZ1 authorization record.

RSP2 asks one late-event numerical question left open by CAL9/PREF13.  It
hash-binds the two terminal SSPRK3 predecessor states, constructs one new
16,385-point input, verifies the unchanged SRC4 source at that input, and
freezes the strict 3/2 constraint-order decision before any new trajectory is
read.  It cannot select a GR-0 case or open SGB-L/FGC-QR execution.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite, log
from pathlib import Path
import subprocess
import sys
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.cal4_common_event_diagnosis import (  # noqa: E402
    owned_semidiscrete_constraint_snapshot,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto11_runtime import (  # noqa: E402
    project_gr0_reference_balanced_state,
)
from recursive_horizons.fgc.evolution.proto12_runtime import (  # noqa: E402
    Proto12GR0EvolutionOperator,
)
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    PROTO4_CONSTRAINT_ORDER,
)
from recursive_horizons.fgc.evolution.rsp2_checkpoint import (  # noqa: E402
    RSP2_CAL9_STAGE_COUNTS,
    RSP2_CAL9_STATE_HASHES,
    RSP2_CAL9_TARGET_TIME,
    file_sha256,
    load_rsp2_cal9_predecessors,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-RSP2-FRZ1"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-rsp2-frz1.toml"
DEFAULT_PLAN = REPOSITORY / "configs/fgc/fgc-1-rsp2-run1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-rsp2-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-rsp2-frz1.md"
CAL9_RESULT = REPOSITORY / "results/fgc-1-cal9-pref13.json"
CAL9_CHECKPOINT = (
    REPOSITORY / "runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz"
)
CHECKPOINT_COMMIT = "da1750ffe4139e806e151846fc9945c07b30b449"
TARGET_COMPONENT = "radial_momentum"
TARGET_INDEX = PROTO4_CONSTRAINT_ORDER.index(TARGET_COMPONENT)

EXPECTED_SCOPE = {
    "role": "pre_trajectory_late_event_one_finer_grid_constraint_order_freeze_and_runtime_authorization",
    "branch": "GR-0",
    "amplitude": "3",
    "method": "SSPRK3",
    "trajectory_advanced": False,
    "CAL9_outcome_reinterpreted": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "mechanism_question_answered": False,
}

EXPECTED_CLAIMS = {
    "RSP2_runtime_implemented": True,
    "RSP2_execution_authorized": True,
    "RSP2_outcome_read": False,
    "RSP2_target_order_cleared": False,
    "RSP2_complete_constraint_admission_passed": False,
    "PROTO13_frozen": False,
    "fresh_GR0_dynamic_calibration_completed": False,
    "classical_spherical_diagnostic_authorized": False,
    "SGBL_execution_authorized": False,
    "FGCQR_holdout_execution_authorized": False,
    "FGCQR_mechanism_rejected": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
    "general_gradient_route_rejected": False,
    "singularity_resolution_derived": False,
    "child_domain_or_topology_derived": False,
    "dark_sector_mechanism_derived": False,
    "varying_locally_measured_c_derived": False,
}

IMPLEMENTATION = (
    REPOSITORY / "scripts/reproduce_fgc_rsp2_frz1.py",
    REPOSITORY / "scripts/run_fgc_rsp2_constraint_study.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/rsp2_checkpoint.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/rsp2_constraint_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto12_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/cal4_common_event_diagnosis.py",
)


def _rel(path: Path) -> str:
    return path.relative_to(REPOSITORY).as_posix()


def _sha(path: Path) -> str:
    return file_sha256(path)


def _serial(value: Any) -> Any:
    if is_dataclass(value) and not isinstance(value, type):
        return _serial(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_serial(item) for item in value]
    if isinstance(value, np.ndarray):
        return _serial(value.tolist())
    if isinstance(value, np.generic):
        return _serial(value.item())
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, float):
        if not isfinite(value):
            raise ValueError("canonical RSP2 evidence cannot contain nonfinite data")
        return value
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    raise TypeError(f"unsupported RSP2 evidence type: {type(value)!r}")


def _canonical(value: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(_serial(value), indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def _digest(value: Any) -> str:
    payload = json.dumps(
        _serial(value), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def _strict_keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(f"{name} keys differ: {sorted(set(value) ^ expected)}")


def _git(*args: str) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *args],
        cwd=REPOSITORY,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def _load_toml(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    if not isinstance(raw, dict):
        raise ValueError(f"{_rel(path)} must contain a TOML table")
    return raw


def _load_canonical_json(path: Path, artifact_id: str) -> dict[str, Any]:
    raw = path.read_bytes()
    value = json.loads(raw)
    if not isinstance(value, dict) or raw != _canonical(value):
        raise ValueError(f"{_rel(path)} is not canonical JSON")
    if value.get("artifact_id") != artifact_id:
        raise ValueError(f"{_rel(path)} artifact identity differs")
    return value


def validate_run_plan_data(plan: Mapping[str, Any]) -> None:
    _strict_keys(
        "RSP2 run plan",
        plan,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "freeze_result",
            "scope",
            "question",
            "physical_inputs",
            "method",
            "numerics",
            "resource_scope",
            "universal_thresholds",
            "endpoint_classification",
            "provenance",
            "claims",
        },
    )
    scope = plan["scope"]
    question = plan["question"]
    method = plan["method"]
    numerics = plan["numerics"]
    physical = plan["physical_inputs"]
    claims = plan["claims"]
    if (
        plan["schema_version"] != 1
        or plan["artifact_id"] != "FGC-1-RSP2-RUN1-PLAN"
        or plan["project_version"] != PROJECT_VERSION
        or plan["metric_signature"] != "-+++"
        or plan["riemann_convention"] != "plus_partial_mu_gamma_nu"
        or plan.get("freeze_result") != "results/fgc-1-rsp2-frz1.json"
        or scope.get("branch") != "GR-0"
        or scope.get("amplitude") != "3"
        or scope.get("method") != "SSPRK3"
        or scope.get("SGBL_or_FGCQR_state_read") is not False
        or question.get("target_component") != TARGET_COMPONENT
        or question.get("combined_point_counts") != [4097, 8193, 16385]
        or question.get("tested_finest_pair") != "8193_to_16385"
        or question.get("unchanged_minimum_finest_pair_order") != "3/2"
        or question.get("strict_greater_than_or_equal_pass_rule") is not True
        or method.get("method_label") != "SSPRK3"
        or method.get("integrator_id")
        != "second_order_diagonal_norm_SBP_plus_SSPRK3"
        or method.get("spatial_order") != 2
        or numerics.get("new_resolution") != 16385
        or numerics.get("final_coordinate_time") != "23/16"
        or numerics.get("minimum_constraint_finest_pair_order") != "3/2"
        or numerics.get("source_point_batch_size") != 2048
        or physical.get("chi_amplitude") != "3"
        or physical.get("outer_radius") != "128"
        or claims.get("RSP2_execution_authorized") is not False
        or any(value is not False for value in claims.values())
    ):
        raise ValueError("RSP2 run plan changed its frozen scientific question")


def validate_run_plan(path: Path = DEFAULT_PLAN) -> dict[str, Any]:
    plan = _load_toml(path)
    validate_run_plan_data(plan)
    return plan


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = _load_toml(path)
    _strict_keys(
        "RSP2 freeze config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "run_plan",
            "owner_document",
            "scope",
            "lineage",
            "study_freeze",
            "predecessor_terminal",
            "runtime_binding",
            "proof_contract",
            "namespace",
            "claims",
        },
    )
    if (
        config.get("schema_version") != 1
        or config.get("artifact_id") != ARTIFACT_ID
        or config.get("project_version") != PROJECT_VERSION
        or config.get("run_plan") != _rel(DEFAULT_PLAN)
        or config.get("owner_document") != _rel(OWNER_DOCUMENT)
        or config.get("scope") != EXPECTED_SCOPE
        or config.get("claims") != EXPECTED_CLAIMS
        or config["study_freeze"].get("combined_point_counts")
        != [4097, 8193, 16385]
        or config["study_freeze"].get("unchanged_minimum_finest_pair_order")
        != "3/2"
        or config["study_freeze"].get("only_the_new_16385_member_is_advanced")
        is not True
        or config["proof_contract"].get("PROTO13_requires_a_separate_post_result_decision")
        is not True
    ):
        raise ValueError("RSP2 freeze config differs from the prospective contract")
    return config


def _parameters(plan: Mapping[str, Any], *, point_amplitude: str = "3") -> PulseParameters:
    physical = plan["physical_inputs"]
    return PulseParameters(
        chi_amplitude=float(Q(point_amplitude)),
        center=float(Q(physical["chi_center"])),
        half_width=float(Q(physical["chi_half_width"])),
        phi_amplitude=float(Q(physical["phi_seed_amplitude"])),
        planck_mass=float(Q(physical["planck_mass"])),
        scalar_mass=float(Q(physical["scalar_mass"])),
        quartic_coupling=float(Q(physical["quartic_coupling"])),
    )


def build_rsp2_initial_state(plan: Mapping[str, Any]):
    physical = plan["physical_inputs"]
    method = plan["method"]
    initial = construct_gr0_grid_initial_data(
        _parameters(plan),
        point_count=plan["numerics"]["new_resolution"],
        outer_radius=float(Q(physical["outer_radius"])),
        constraint_method=method["constraint_solve_method"],
        diagnostic_spatial_order=method["spatial_order"],
    )
    state = project_gr0_reference_balanced_state(
        initial, spatial_order=method["spatial_order"]
    )
    return initial, state


def rsp2_operator(plan: Mapping[str, Any], grid):
    numerics = plan["numerics"]
    thresholds = plan["universal_thresholds"]
    return Proto12GR0EvolutionOperator(
        grid,
        spatial_order=plan["method"]["spatial_order"],
        ko_dissipation=float(Q(numerics["ko_dissipation"])),
        raw_tolerance=float(Q(thresholds["source_residual_infinity_max"])),
        kinetic_condition_maximum=float(
            Q(thresholds["kinetic_condition_number_max"])
        ),
        maximum_refinement_iterations=thresholds["source_iteration_max"],
        point_batch_size=numerics["source_point_batch_size"],
    )


def _validate_lineage(config: Mapping[str, Any]) -> tuple[dict[str, Any], Any]:
    lineage = config["lineage"]
    if (
        lineage.get("checkpoint_commit") != CHECKPOINT_COMMIT
        or _git("merge-base", "--is-ancestor", CHECKPOINT_COMMIT, "HEAD").returncode
        != 0
        or _sha(CAL9_RESULT) != lineage.get("cal9_result_sha256")
        or _sha(REPOSITORY / lineage["cal9_config"])
        != lineage.get("cal9_config_sha256")
        or _sha(CAL9_CHECKPOINT) != lineage.get("cal9_checkpoint_sha256")
    ):
        raise ValueError("RSP2 immutable CAL9 lineage differs")
    cal9 = _load_canonical_json(CAL9_RESULT, "FGC-1-CAL9-PREF13")
    gates = cal9.get("gate_status", {})
    if (
        gates.get("PROTO12_campaign_terminated_normally") is not True
        or gates.get("PROTO12_GR0_case_eligible") is not False
        or gates.get("RSP2_high_ladder_constraint_preflight_required") is not True
        or gates.get("FGCQR_holdout_execution_authorized") is not False
        or any(value is not False for value in cal9.get("nonclaims", {}).values())
    ):
        raise ValueError("CAL9 does not require the narrow RSP2 successor")
    predecessors = load_rsp2_cal9_predecessors(
        CAL9_CHECKPOINT,
        expected_sha256=lineage["cal9_checkpoint_sha256"],
        outer_radius=128.0,
    )
    return cal9, predecessors


def _predecessor_evidence(plan: Mapping[str, Any], predecessors: Any) -> list[dict[str, Any]]:
    records = []
    effective = []
    for state, grid, stage_count, input_hash, state_hash in zip(
        predecessors.states,
        predecessors.grids,
        predecessors.accepted_stage_counts,
        predecessors.input_hashes,
        predecessors.state_hashes,
        strict=True,
    ):
        snapshot = owned_semidiscrete_constraint_snapshot(
            state,
            grid,
            coordinate_time=predecessors.coordinate_time,
            diagnostic_spatial_order=2,
            fixed_outer_rows=plan["numerics"]["fixed_outer_rows"],
            accepted_stage_count=stage_count,
            planck_mass=float(Q(plan["physical_inputs"]["planck_mass"])),
            scalar_mass=float(Q(plan["physical_inputs"]["scalar_mass"])),
            quartic_coupling=float(Q(plan["physical_inputs"]["quartic_coupling"])),
            length_unit=float(Q(plan["physical_inputs"]["length_unit_L0"])),
        )
        raw = float(snapshot.owned_domain_norms.component_infinity[TARGET_INDEX])
        enclosed = float(
            snapshot.effective_owned_norms_for_order_only.component_infinity[
                TARGET_INDEX
            ]
        )
        effective.append(enclosed)
        records.append(
            {
                "method": "SSPRK3",
                "point_count": grid.point_count,
                "coordinate_time": predecessors.coordinate_time,
                "accepted_stage_count": stage_count,
                "input_hash": input_hash,
                "state_sha256": state_hash,
                "radial_momentum_raw_owned_norm": raw,
                "accumulated_roundoff_enclosure": snapshot.accumulated_roundoff_enclosure,
                "radial_momentum_effective_norm_for_order_only": enclosed,
                "radial_momentum_owned_maximum_radius": float(
                    snapshot.owned_component_maximum_radii[TARGET_INDEX]
                ),
            }
        )
    observed_order = log(effective[0] / effective[1]) / log(2.0)
    expected = float(Q(plan["question"]["predecessor_pair_order"]))
    if observed_order != expected:
        raise ValueError("RSP2 predecessor radial-momentum order does not recompute")
    frozen = (
        plan["question"]["predecessor_pair_shortfall"],
        plan["question"]["predecessor_pair_order"],
    )
    if frozen != ("0.0009052615636708783", "1.4990947384363291"):
        raise ValueError("RSP2 predecessor question differs")
    return records


def _new_input_evidence(plan: Mapping[str, Any]) -> tuple[dict[str, Any], Any, Any]:
    initial, state = build_rsp2_initial_state(plan)
    before = array_content_sha256(state.u, state.p, state.q)
    operator = rsp2_operator(plan, initial.grid)
    rhs = operator(0.0, state)
    after = array_content_sha256(state.u, state.p, state.q)
    diagnostics = rhs.diagnostics
    source_maximum = float(
        Q(plan["universal_thresholds"]["source_residual_infinity_max"])
    )
    source_precheck = _source_precheck_certificate(
        diagnostics["source_residual_infinity"], source_maximum
    )
    if (
        before != after
        or diagnostics.get("source_raw_gate_passed") is not True
        or diagnostics.get("SRC4_tensor_contracted_reference_source") is not True
        or diagnostics.get("PROTO11_interior_q_reprojected") is not False
        or diagnostics.get("PROTO11_exact_reference_equilibrium_applied") is not False
    ):
        raise ValueError("RSP2 16385-point initial source precheck failed")
    record = {
        "amplitude": "3",
        "method": "SSPRK3",
        "integrator_id": plan["method"]["integrator_id"],
        "spatial_order": 2,
        "point_count": 16385,
        "grid_spacing": initial.grid.spacing,
        "projected_state_sha256": before,
        "source_instrument": "FGC-1-SRC4-VEC1",
        "reference_state_map": "PROTO11",
        "source_residual_precheck": source_precheck,
        "kinetic_condition_infinity": diagnostics["kinetic_condition_infinity"],
        "source_point_batch_size": plan["numerics"]["source_point_batch_size"],
        "accepted_state_source_gate_passed": True,
        "operator_input_state_unchanged": True,
        "q_remains_independently_evolved": True,
        "interior_q_reprojected": False,
        "trajectory_advanced": False,
    }
    record["expanded_run_config_sha256"] = _digest(record)
    return record, initial, state


def _source_precheck_certificate(
    observed_residual: float, maximum_residual: float
) -> dict[str, Any]:
    """Return a stable authorization fact for a thread-sensitive reduction.

    The binary64 infinity norm is still emitted by the trajectory diagnostics,
    but its last bits can depend on the platform linear-algebra reduction.  It
    is therefore evidence for the source gate, not an input to the immutable
    study identity.  The exact projected-state hash remains the input binding.
    """

    observed = float(observed_residual)
    maximum = float(maximum_residual)
    if not isfinite(observed) or observed < 0.0:
        raise ValueError("RSP2 source residual must be finite and nonnegative")
    if not isfinite(maximum) or maximum <= 0.0:
        raise ValueError("RSP2 source residual maximum must be finite and positive")
    if observed > maximum:
        raise ValueError("RSP2 source residual exceeds the frozen maximum")
    return {
        "maximum": maximum,
        "passed": True,
        "raw_binary64_reduction_used_in_study_identity": False,
        "raw_value_preserved_by_runtime_diagnostics": True,
    }


def _runtime_controls(plan: Mapping[str, Any]) -> dict[str, Any]:
    # Exercise one independent small state so the mutation control is cheap and
    # cannot read or advance the frozen 16,385-point trajectory.
    small = construct_gr0_grid_initial_data(
        _parameters(plan),
        point_count=129,
        outer_radius=float(Q(plan["physical_inputs"]["outer_radius"])),
        constraint_method="SSPRK3",
        diagnostic_spatial_order=2,
    )
    baseline = project_gr0_reference_balanced_state(small, spatial_order=2)
    operator = rsp2_operator(plan, small.grid)
    base_rhs = operator(0.0, baseline)
    changed_q = baseline.q.copy()
    changed_bits = changed_q.view(np.uint64)
    changed_bits[11, 5] ^= np.uint64(1 << 40)
    changed = EvolutionState(baseline.u, baseline.p, changed_q)
    changed_rhs = operator(0.0, changed)
    visible = array_content_sha256(base_rhs.dp) != array_content_sha256(
        changed_rhs.dp
    )
    if not visible:
        raise ValueError("RSP2 one-bit SRC4 mutation is not visible")

    mutation_controls = {}
    for name, mutate in {
        "grid": lambda value: value["numerics"].__setitem__("new_resolution", 8193),
        "threshold": lambda value: value["question"].__setitem__(
            "unchanged_minimum_finest_pair_order", "149/100"
        ),
        "time": lambda value: value["numerics"].__setitem__(
            "final_coordinate_time", "3/2"
        ),
        "method": lambda value: value["method"].__setitem__(
            "method_label", "RK4"
        ),
        "component": lambda value: value["question"].__setitem__(
            "target_component", "gauge_r"
        ),
        "claim": lambda value: value["claims"].__setitem__(
            "FGCQR_holdout_execution_authorized", True
        ),
    }.items():
        candidate = json.loads(json.dumps(plan))
        mutate(candidate)
        try:
            validate_run_plan_data(candidate)
        except (TypeError, ValueError, KeyError):
            mutation_controls[name] = True
        else:
            mutation_controls[name] = False
    if not all(mutation_controls.values()):
        raise ValueError("RSP2 plan mutation control failed closed")
    return {
        "one_bit_SRC4_input_mutation_visible": True,
        "mutation_controls": mutation_controls,
        "invalid_runtime_stop_cannot_be_persistent": True,
        "unexpected_exception_classifies_invalid": True,
        "all_controls_passed": True,
    }


def _namespace_evidence(
    config: Mapping[str, Any], historical: Mapping[str, Any] | None
) -> dict[str, Any]:
    relative = config["namespace"]["output_root"]
    if historical is not None:
        expected = {
            "path": relative,
            "empty_before_authorization": True,
            "authorization_created_no_namespace": True,
        }
        if historical != expected:
            raise ValueError("historical RSP2 namespace evidence differs")
        return dict(historical)
    path = REPOSITORY / relative
    if path.exists():
        raise ValueError("RSP2 output namespace must be absent before authorization")
    return {
        "path": relative,
        "empty_before_authorization": True,
        "authorization_created_no_namespace": True,
    }


def reproduce(
    config_path: Path = DEFAULT_CONFIG,
    *,
    historical_namespace_evidence: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    config = load_config(config_path)
    plan_path = REPOSITORY / config["run_plan"]
    plan = validate_run_plan(plan_path)
    cal9, predecessors = _validate_lineage(config)
    predecessor_records = _predecessor_evidence(plan, predecessors)
    new_input, _initial, _state = _new_input_evidence(plan)
    controls = _runtime_controls(plan)
    namespace = _namespace_evidence(config, historical_namespace_evidence)
    implementation = {_rel(path): _sha(path) for path in IMPLEMENTATION}
    plan_hash = _sha(plan_path)
    input_manifest = {
        "predecessor_state_sha256": list(predecessors.state_hashes),
        "new_input": new_input,
    }
    study_id = _digest(
        {
            "artifact_id": plan["artifact_id"],
            "plan_sha256": plan_hash,
            "input_manifest_sha256": _digest(input_manifest),
            "checkpoint_sha256": predecessors.checkpoint_sha256,
            "target_coordinate_time": "23/16",
            "point_counts": [4097, 8193, 16385],
            "output_root": config["namespace"]["output_root"],
        }
    )
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "pre_trajectory_late_event_one_finer_grid_constraint_order_runtime_authorization",
        "generated_by": _rel(Path(__file__)),
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            _rel(config_path): _sha(config_path),
            _rel(plan_path): plan_hash,
        },
        "predecessor_sha256": {
            _rel(CAL9_RESULT): _sha(CAL9_RESULT),
            _rel(CAL9_CHECKPOINT): predecessors.checkpoint_sha256,
        },
        "implementation_sha256": implementation,
        "scope_bindings": dict(config["scope"]),
        "gate_status": dict(EXPECTED_CLAIMS),
        "artifact_payload": {
            "immutable_lineage": {
                "checkpoint_commit": CHECKPOINT_COMMIT,
                "checkpoint_is_ancestor_of_HEAD": True,
                "CAL9_canonical_hash_matched": True,
                "CAL9_checkpoint_hash_matched": True,
                "CAL9_campaign_id": predecessors.campaign_id,
                "CAL9_event_log_sha256": predecessors.event_log_sha256,
                "CAL9_event_log_line_count": predecessors.event_log_line_count,
                "CAL9_no_GR0_case_eligible": cal9["gate_status"][
                    "PROTO12_GR0_case_eligible"
                ]
                is False,
            },
            "frozen_run_plan": {
                "artifact_id": plan["artifact_id"],
                "run_plan_sha256": plan_hash,
                "study_id": study_id,
                "input_manifest_sha256": _digest(input_manifest),
                "amplitude": "3",
                "method": "SSPRK3",
                "predecessor_point_counts": [4097, 8193],
                "new_point_count": 16385,
                "combined_point_counts": [4097, 8193, 16385],
                "target_coordinate_time": "23/16",
                "target_component": TARGET_COMPONENT,
                "minimum_finest_pair_order": 1.5,
                "strict_greater_than_or_equal_pass_rule": True,
                "only_the_new_16385_member_is_advanced": True,
            },
            "predecessor_terminal_evidence": predecessor_records,
            "new_run_input": new_input,
            "runtime_controls": controls,
            "namespace_precondition": namespace,
            "decision": {
                "RSP2_execution_authorized": True,
                "RSP2_outcome_read": False,
                "PROTO13_requires_separate_post_result_decision": True,
            },
            "epistemic_boundary": {
                "passing_target_means_only_preasymptotic_on_tested_ladder": True,
                "failing_target_means_only_persistent_through_tested_pair": True,
                "runtime_stop_or_invalid_run_is_not_persistence": True,
                "threshold_reduction_fit_or_post_outcome_tuning_authorized": False,
                "candidate_or_mechanism_outcome_read": False,
            },
        },
        "nonclaims": {
            key: False
            for key in EXPECTED_CLAIMS
            if key not in {"RSP2_runtime_implemented", "RSP2_execution_authorized"}
        },
    }


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    return _load_canonical_json(path, ARTIFACT_ID)


def verify_canonical(
    path: Path = DEFAULT_OUTPUT, config_path: Path = DEFAULT_CONFIG
) -> dict[str, Any]:
    observed = load_canonical_result(path)
    historical = observed["artifact_payload"]["namespace_precondition"]
    reproduced = reproduce(
        config_path, historical_namespace_evidence=historical
    )
    if _canonical(observed) != _canonical(reproduced):
        raise ValueError("canonical FGC-1-RSP2-FRZ1 record differs")
    return observed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    config_path = args.config.resolve()
    output_path = args.output.resolve()
    if args.check:
        verify_canonical(output_path, config_path)
        print(f"verified {output_path}")
        return
    record = reproduce(config_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_bytes(_canonical(record))
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
