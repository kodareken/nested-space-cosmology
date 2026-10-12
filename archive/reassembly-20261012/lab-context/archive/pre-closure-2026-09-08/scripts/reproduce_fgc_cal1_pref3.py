#!/usr/bin/env python3
"""Regenerate the FGC-1-CAL1-PREF3 semidiscrete composition diagnosis."""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
import sys
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-cal1-pref3.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-cal1-pref3.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-cal1-pref3.md"

from recursive_horizons.fgc.evolution.calibration_runtime import (  # noqa: E402
    CONSTRAINT_ROUNDOFF_OPERATION_BUDGET,
    gr0_semidiscrete_constraint_snapshot,
    project_gr0_semidiscrete_state,
    proto5_semidiscrete_constraint_admission,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    ConstraintEventSample,
    common_event_constraint_admission,
)
from recursive_horizons.fgc.initial_data_preflight import PulseParameters  # noqa: E402
from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    validate_sf1_protocol,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-CAL1-PREF3"
PROJECT_VERSION = "0.11.0"
REQUIRED_GATE = "PROTO4_semidiscrete_common_event_contract_obstructed"
EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v4.toml",
    "protocol_freeze_result": "results/fgc-1-pro4-frz1.json",
    "successor_monitor_result": "results/fgc-1-hlt2-mon2.json",
    "static_ledger_result": "results/fgc-1-id2-all1.json",
}
EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO4",
    "role": "pre_trajectory_semidiscrete_common_event_composition_audit",
    "calibration_branch": "GR-0",
    "eligible_amplitudes": ["5/2", "3"],
    "methods": ["RK4", "SSPRK3"],
    "resolutions": [1025, 2049, 4097],
    "calibration_trajectory_read": False,
    "calibration_output_namespace_read_or_written": False,
    "holdout_output_namespace_read_or_written": False,
    "FGCQR_evolution_outcome_read": False,
    "collapse_or_trapped_outcome_classified": False,
}
EXPECTED_PROTO4 = {
    "continuum_initial_q_source": "ID2_complete_grid_analytic_first_derivatives",
    "kinematic_reduction_constraint": "q_minus_D_h_u",
    "coarsest_normalized_guard_max": "1/50",
    "finest_normalized_guard_max": "1/1000",
    "minimum_every_adjacent_pair_order": "3/2",
    "exact_zero_pair_special_case_only": True,
    "same_guards_and_pair_rule_for_both_methods": True,
}
EXPECTED_PREVIEW = {
    "physical_u_and_p_are_bitwise_preserved": True,
    "semidiscrete_q_projection": "q_equals_native_SBP_D_h_u_before_first_stage",
    "raw_normalized_residuals_decide_magnitude_guards": True,
    "roundoff_operation_budget": CONSTRAINT_ROUNDOFF_OPERATION_BUDGET,
    "roundoff_enclosure_use": "convergence_zero_classification_only_and_public_nonnegative_error_component",
    "primary_coarsest_guard_max": "1/50",
    "primary_finest_guard_max": "1/1000",
    "comparator_coarsest_guard_max": "1/10",
    "comparator_finest_guard_max": "1/200",
    "minimum_finest_adjacent_pair_order": "3/2",
    "earlier_pair_must_be_monotone_but_is_diagnostic_for_asymptotic_order": True,
    "new_protocol_version_required": True,
}
EXPECTED_PROOF = {
    "PROTO4_HLT2_and_ID2_must_be_canonical_and_passed",
    "only_ID2_eligible_amplitudes_may_be_composed",
    "both_methods_and_all_three_grids_must_be_evaluated",
    "original_ID2_state_and_projected_semidiscrete_state_must_both_be_hash_bound",
    "original_PROTO4_common_event_decision_must_be_executed_without_override",
    "raw_physical_gauge_and_reduction_norms_must_be_serialized",
    "roundoff_enclosure_must_not_change_raw_magnitude_guards",
    "prospective_repair_must_pass_at_t0_for_every_eligible_input",
    "no_dynamic_output_or_candidate_outcome_may_be_read",
    "protocol_obstruction_is_not_a_mechanism_obstruction",
    "canonical_hash_bound_result_required",
}
EXPECTED_CLAIMS = {
    "PROTO4_fresh_GR0_dynamic_calibration_authorized": False,
    "PROTO5_premise_revision_required": True,
    "fresh_GR0_dynamic_calibration_completed": False,
    "FGCQR_holdout_execution_authorized": False,
    "classical_spherical_diagnostic_authorized": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
    "FGCQR_mechanism_rejected": False,
    "general_gradient_route_rejected": False,
    "singularity_resolution_derived": False,
    "child_domain_or_topology_derived": False,
    "dark_sector_mechanism_derived": False,
    "varying_locally_measured_c_derived": False,
}
IMPLEMENTATION = tuple(
    REPOSITORY / name
    for name in (
        "src/recursive_horizons/fgc/evolution/calibration_runtime.py",
        "src/recursive_horizons/fgc/evolution/gr0_calibration.py",
        "src/recursive_horizons/fgc/evolution/gr0_direct_source.py",
        "src/recursive_horizons/fgc/evolution/proto4_admission.py",
        "scripts/reproduce_fgc_cal1_pref3.py",
    )
)


def _pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


def _serial(value: Any) -> Any:
    if isinstance(value, Fraction):
        return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"
    if is_dataclass(value):
        return _serial(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return _serial(value.tolist())
    if isinstance(value, np.generic):
        return _serial(value.item())
    if isinstance(value, (tuple, list)):
        return [_serial(item) for item in value]
    if isinstance(value, float):
        if not isfinite(value):
            raise ValueError("nonfinite value cannot enter canonical evidence")
        return value
    if value is None or isinstance(value, (str, bool, int)):
        return value
    raise TypeError(f"unsupported evidence value {type(value).__name__}")


def _canonical(value: Any) -> str:
    return json.dumps(_serial(value), indent=2, sort_keys=True, allow_nan=False) + "\n"


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _load_json(path: Path, name: str) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    value = json.loads(source, object_pairs_hook=_pairs)
    if not isinstance(value, dict) or source != _canonical(value):
        raise ValueError(f"{name} must be canonical unique-key JSON")
    return value


def _path(name: str, value: Any, expected: str) -> Path:
    if value != expected:
        raise ValueError(f"{name} must name {expected}")
    path = (REPOSITORY / expected).resolve()
    if not path.is_file() or _rel(path) != expected:
        raise ValueError(f"{name} must be a canonical repository file")
    return path


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    expected_root = {
        "schema_version", "artifact_id", "project_version", "metric_signature",
        "riemann_convention", *EXPECTED_PATHS, "scope", "PROTO4_contract",
        "prospective_repair_preview", "proof_contract", "claims",
    }
    if set(raw) != expected_root:
        raise ValueError("CAL1 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("CAL1 identity or convention differs")
    if raw["scope"] != EXPECTED_SCOPE:
        raise ValueError("CAL1 scope differs")
    if raw["PROTO4_contract"] != EXPECTED_PROTO4:
        raise ValueError("CAL1 PROTO4 contract differs")
    if raw["prospective_repair_preview"] != EXPECTED_PREVIEW:
        raise ValueError("CAL1 prospective repair differs")
    if set(raw["proof_contract"]) != EXPECTED_PROOF or any(
        value is not True for value in raw["proof_contract"].values()
    ):
        raise ValueError("CAL1 proof contract differs")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("CAL1 claims differ")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    with paths["protocol_config"].open("rb") as handle:
        protocol = validate_sf1_protocol(tomllib.load(handle))
    if protocol["artifact_id"] != "FGC-2-SF1-PROTO4":
        raise ValueError("CAL1 protocol identity differs")
    requirements = (
        ("protocol_freeze_result", "FGC-1-PRO4-FRZ1", "PROTO4_outcome_neutral_protocol_frozen"),
        ("successor_monitor_result", "FGC-1-HLT2-MON2", "PROTO4_successor_admission_monitor_implemented"),
        ("static_ledger_result", "FGC-1-ID2-ALL1", "all_amplitude_all_case_static_initial_admission_completed"),
    )
    predecessors = {}
    for key, artifact, gate in requirements:
        result = _load_json(paths[key], artifact)
        if result.get("artifact_id") != artifact or result.get("gate_status", {}).get(gate) is not True:
            raise ValueError(f"{artifact} required gate differs")
        predecessors[key] = result
    eligible = [
        item["amplitude"]
        for item in predecessors["static_ledger_result"]["artifact_payload"]["dynamic_calibration_eligibility"]
        if item["eligible_for_fresh_dynamic_GR0_calibration"]
    ]
    if eligible != EXPECTED_SCOPE["eligible_amplitudes"]:
        raise ValueError("ID2 eligible amplitude order differs")
    return {"raw": raw, "paths": paths, "predecessors": predecessors}


def _method_record(amplitude_text: str, method: str, order: int) -> dict[str, Any]:
    original_samples = []
    projected_samples = []
    resolutions = []
    for point_count in EXPECTED_SCOPE["resolutions"]:
        initial = construct_gr0_grid_initial_data(
            PulseParameters(chi_amplitude=float(Q(amplitude_text))),
            point_count=point_count,
            constraint_method=method,
            diagnostic_spatial_order=order,
        )
        original = gr0_semidiscrete_constraint_snapshot(
            initial.state,
            initial.grid,
            coordinate_time=0.0,
            diagnostic_spatial_order=order,
        )
        projected_state = project_gr0_semidiscrete_state(
            initial,
            spatial_order=order,
        )
        projected = gr0_semidiscrete_constraint_snapshot(
            projected_state,
            initial.grid,
            coordinate_time=0.0,
            diagnostic_spatial_order=order,
        )
        original_samples.append(original)
        projected_samples.append(projected)
        resolutions.append({
            "point_count": point_count,
            "original_ID2_state_sha256": array_content_sha256(
                initial.state.u, initial.state.p, initial.state.q
            ),
            "projected_semidiscrete_state_sha256": array_content_sha256(
                projected_state.u, projected_state.p, projected_state.q
            ),
            "u_bitwise_preserved": bool(np.array_equal(initial.state.u, projected_state.u)),
            "p_bitwise_preserved": bool(np.array_equal(initial.state.p, projected_state.p)),
            "q_projection_change_infinity": float(
                np.max(np.abs(initial.state.q - projected_state.q), initial=0.0)
            ),
            "original_raw_global_constraint": original.raw_norms.global_infinity,
            "original_raw_component_constraints": dict(zip(
                original.raw_norms.component_names,
                original.raw_norms.component_infinity,
                strict=True,
            )),
            "projected_raw_global_constraint": projected.raw_norms.global_infinity,
            "projected_raw_component_constraints": dict(zip(
                projected.raw_norms.component_names,
                projected.raw_norms.component_infinity,
                strict=True,
            )),
            "normalized_roundoff_zero_enclosure": dict(zip(
                projected.raw_norms.component_names,
                projected.normalized_roundoff_zero_enclosure,
                strict=True,
            )),
        })

    def proto4(records):
        return common_event_constraint_admission(tuple(
            ConstraintEventSample(item.point_count, item.coordinate_time, item.raw_norms)
            for item in records
        ))

    original_decision = proto4(original_samples)
    projected_decision = proto4(projected_samples)
    successor = proto5_semidiscrete_constraint_admission(
        projected_samples,
        method=method,
        coarsest_guard_maximum=(1.0 / 50.0 if order == 4 else 1.0 / 10.0),
        finest_guard_maximum=(1.0 / 1000.0 if order == 4 else 1.0 / 200.0),
    )
    return {
        "method": method,
        "spatial_order": order,
        "resolutions": resolutions,
        "original_ID2_state_under_literal_PROTO4": original_decision,
        "projected_state_under_literal_PROTO4": projected_decision,
        "projected_state_under_prospective_repair": successor,
    }


def record(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    loaded = load_config(path)
    methods = (("RK4", 4), ("SSPRK3", 2))
    candidates = []
    for amplitude in EXPECTED_SCOPE["eligible_amplitudes"]:
        candidates.append({
            "amplitude": amplitude,
            "methods": [
                _method_record(amplitude, method, order)
                for method, order in methods
            ],
        })
    literal_original_passes = [
        item["original_ID2_state_under_literal_PROTO4"].admission_passed
        for candidate in candidates for item in candidate["methods"]
    ]
    literal_projected_passes = [
        item["projected_state_under_literal_PROTO4"].admission_passed
        for candidate in candidates for item in candidate["methods"]
    ]
    preview_passes = [
        item["projected_state_under_prospective_repair"].admission_passed
        for candidate in candidates for item in candidate["methods"]
    ]
    if any(literal_original_passes) or any(literal_projected_passes) or not all(preview_passes):
        raise RuntimeError("CAL1 composition controls do not establish the frozen diagnosis")
    raw = loaded["raw"]
    paths = loaded["paths"]
    gate_status = {
        REQUIRED_GATE: True,
        "PROTO4_original_ID2_states_fail_literal_common_event_admission": True,
        "PROTO4_projected_states_still_fail_literal_common_event_admission": True,
        "prospective_semidiscrete_repair_passes_every_t0_composition": True,
        **raw["claims"],
    }
    return _serial({
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "pre_trajectory_semidiscrete_test_contract_obstruction_not_mechanism_result",
        "generated_by": "scripts/reproduce_fgc_cal1_pref3.py",
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(path): _sha(path)},
        "predecessor_sha256": {
            raw[name]: _sha(paths[name]) for name in EXPECTED_PATHS
        },
        "implementation_sha256": {
            _rel(item): _sha(item) for item in IMPLEMENTATION
        },
        "scope_bindings": {
            "target_protocol": "FGC-2-SF1-PROTO4",
            "calibration_branch": "GR-0",
            "eligible_amplitudes": EXPECTED_SCOPE["eligible_amplitudes"],
            "methods": EXPECTED_SCOPE["methods"],
            "resolutions": EXPECTED_SCOPE["resolutions"],
            "dynamic_output_namespace_accessed": False,
        },
        "certificate_contract": sorted(EXPECTED_PROOF),
        "artifact_payload": {
            "decision": "PROTO4_cannot_authorize_fresh_calibration_until_a_versioned_semidiscrete_premise_repair_is_frozen",
            "candidate_compositions": candidates,
            "diagnosis": {
                "missing_continuum_to_semidiscrete_q_map": True,
                "single_method_independent_constraint_guards_reject_declared_second_order_comparator_at_t0": True,
                "exact_zero_only_order_rule_turns_binary64_roundoff_into_fake_nonconvergence": True,
                "physical_equations_or_candidate_outcome_implicated": False,
                "new_protocol_version_required": True,
            },
            "epistemic_boundary": {
                "trajectory_read": False,
                "collapse_or_trapping_classified": False,
                "FGCQR_outcome_read": False,
                "mechanism_question_answered": False,
            },
        },
        "gate_status": gate_status,
        "nonclaims": {
            "fresh_GR0_calibration_completed": False,
            "trapped_sphere_observed": False,
            "FGCQR_holdout_evolved": False,
            "regulator_activation_observed": False,
            "positive_Raychaudhuri_margin_observed": False,
            "FGCQR_mechanism_rejected": False,
            "general_gradient_route_rejected": False,
            "retained_EFT_validity": False,
            "physical_transition": False,
        },
    })


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    return _load_json(path, "CAL1 result")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    payload = record(arguments.config.resolve())
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(_canonical(payload), encoding="utf-8")
    print(
        f"wrote {arguments.output} ({REQUIRED_GATE}=true; trajectory_read=false; mechanism_result=false)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
