#!/usr/bin/env python3
"""Regenerate the FGC-1-BND2-CP1 causal-isolation certificate."""

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


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-bnd2-cp1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-bnd2-cp1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-bnd2-cp1.md"

from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    ALL_CONE_LOCAL_SPEED_BOUND,
    BoundaryControlStop,
    BoundaryGeometry,
    CausalBudgetState,
    LOCAL_CONE_SPEED_BOUNDS,
    accept_causal_step,
    aligned_outer_boundary_grids,
    exact_coordinate_speed_upper_bound,
    exact_reference_budget,
)
from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    FUTURE_COMMON_NONCLAIMS,
    FUTURE_EVIDENCE_CLASSIFICATIONS,
    RUN1_SYM1_ARTIFACT_ID,
    SF1_PROTOCOL_ARTIFACT_ID,
    validate_sf1_protocol,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-BND2-CP1"
PROJECT_VERSION = "0.11.0"
REQUIRED_GATE = "spherical_boundary_or_domain_of_dependence_control_passed"
EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v3.toml",
    "run1_config": "configs/fgc/fgc-1-run1-sym1.toml",
    "action_config": "configs/fgc/fgc-1-action-gate.toml",
    "action_result": "results/fgc-1-action-gate.json",
    "variation_config": "configs/fgc/fgc-1-metric-variation.toml",
    "variation_result": "results/fgc-1-metric-variation.json",
    "frozen_boundary_config": "configs/fgc/fgc-1-hyp1-bnd1-md1.toml",
    "frozen_boundary_result": "results/fgc-1-hyp1-bnd1-md1.json",
    "physical_constraints_config": "configs/fgc/fgc-1-con4-phy1.toml",
    "physical_constraints_result": "results/fgc-1-con4-phy1.json",
    "initial_data_config": "configs/fgc/fgc-1-id1-fam1.toml",
    "initial_data_result": "results/fgc-1-id1-fam1.json",
    "run_domain_config": "configs/fgc/fgc-1-dom4-run1.toml",
    "run_domain_result": "results/fgc-1-dom4-run1.json",
    "multidirectional_health_config": "configs/fgc/fgc-1-hyp2-md1.toml",
    "multidirectional_health_result": "results/fgc-1-hyp2-md1.json",
}
RESULT_REQUIREMENTS = {
    "action_result": ("FGC-1-ACT1", "declared_action_and_scalar_algebra_certificate_passed"),
    "variation_result": ("FGC-1-VAR1", "var1_covariant_metric_equation_derived_and_cross_checked"),
    "frozen_boundary_result": ("FGC-1-HYP1-BND1-MD1", "uniform_frozen_radial_main_system_boundary_dissipation_passed"),
    "physical_constraints_result": ("FGC-1-CON4-PHY1", "physical_gauge_reduction_constraint_system_closed"),
    "initial_data_result": ("FGC-1-ID1-FAM1", "nonzero_width_finite_mass_constraint_compatible_family_constructed"),
    "run_domain_result": ("FGC-1-DOM4-RUN1", "nonzero_classical_spherical_run_envelope_passed"),
    "multidirectional_health_result": ("FGC-1-HYP2-MD1", "quantitative_all_covector_weak_coupling_health_envelope_passed"),
}
IMPLEMENTATION = tuple(
    REPOSITORY / path
    for path in (
        "src/recursive_horizons/fgc/evolution/boundary_domain.py",
        "scripts/reproduce_fgc_bnd2_cp1.py",
    )
)


def _pairs(pairs):
    answer = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON object key: {key}")
        answer[key] = value
    return answer


def _mapping(name: str, value: Any) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a table")
    return value


def _keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(
            f"{name} keys differ: missing={sorted(expected - set(value))}, "
            f"extra={sorted(set(value) - expected)}"
        )


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _fraction(name: str, value: Any, *, positive: bool = False, nonnegative: bool = False) -> Fraction:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a canonical rational string")
    try:
        answer = Q(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} must be a canonical rational string") from exc
    if _fraction_text(answer) != value:
        raise ValueError(f"{name} must use canonical rational encoding")
    if positive and answer <= 0:
        raise ValueError(f"{name} must be positive")
    if nonnegative and answer < 0:
        raise ValueError(f"{name} must be nonnegative")
    return answer


def _path(name: str, value: Any, expected: str) -> Path:
    if value != expected:
        raise ValueError(f"{name} must name {expected}")
    path = (REPOSITORY / expected).resolve()
    if not path.is_file() or path.relative_to(REPOSITORY).as_posix() != expected:
        raise ValueError(f"{name} must be canonical, traversal free, and existing")
    return path


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _load_json(path: Path, name: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs)
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{name} must be valid unique-key JSON") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a JSON object")
    return value


def _serial(value: Any) -> Any:
    if isinstance(value, Fraction):
        return _fraction_text(value)
    if is_dataclass(value):
        return _serial(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
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


def _load_protocol(path: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    return raw, validate_sf1_protocol(raw)


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    top = {
        "schema_version", "artifact_id", "project_version", "metric_signature",
        "riemann_convention", *EXPECTED_PATHS, "scope", "cone_envelope",
        "causal_budget", "discretization", "outer_boundary_move",
        "constraint_control", "proof_contract", "claims",
    }
    _keys("config", raw, top)
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("BND2 identity or convention differs")
    paths = {name: _path(name, raw[name], expected) for name, expected in EXPECTED_PATHS.items()}

    scope = _mapping("scope", raw["scope"])
    _keys("scope", scope, {"target_branch", "physical_equations", "control_route", "complete_nonlinear_constraint_preserving_IBVP_proved", "holdout_execution_performed"})
    if scope != {
        "target_branch": "FGC-QR",
        "physical_equations": "unredefined_ACT1_VAR1",
        "control_route": "boundary_free_measurement_domain_by_dynamic_all_cone_causal_exclusion",
        "complete_nonlinear_constraint_preserving_IBVP_proved": False,
        "holdout_execution_performed": False,
    }:
        raise ValueError("BND2 scope differs")

    cones = _mapping("cone_envelope", raw["cone_envelope"])
    _keys("cone_envelope", cones, {"frame", "ADM_metric", "coordinate_speed_formula", "physical_local_speed_absolute_maximum", "tilde_local_speed_absolute_maximum", "hat_local_speed_absolute_maximum", "all_cone_local_speed_absolute_maximum", "evolving_lapse_shift_and_radial_metric_required"})
    cone_values = {
        name: _fraction(f"cone_envelope.{name}", cones[name], positive=True)
        for name in (
            "physical_local_speed_absolute_maximum",
            "tilde_local_speed_absolute_maximum",
            "hat_local_speed_absolute_maximum",
            "all_cone_local_speed_absolute_maximum",
        )
    }
    if (
        cones["frame"] != "physical_orthonormal_frame"
        or cones["ADM_metric"] != "ds2=-N2dt2+Lambda2(dr+shift_dt)2+R2dOmega2"
        or cones["coordinate_speed_formula"] != "abs_shift_plus_lapse_over_Lambda_times_abs_z"
        or cones["evolving_lapse_shift_and_radial_metric_required"] is not True
        or cone_values != {
            "physical_local_speed_absolute_maximum": LOCAL_CONE_SPEED_BOUNDS["physical"],
            "tilde_local_speed_absolute_maximum": LOCAL_CONE_SPEED_BOUNDS["tilde"],
            "hat_local_speed_absolute_maximum": LOCAL_CONE_SPEED_BOUNDS["hat"],
            "all_cone_local_speed_absolute_maximum": ALL_CONE_LOCAL_SPEED_BOUND,
        }
    ):
        raise ValueError("BND2 cone envelope differs")

    budget = _mapping("causal_budget", raw["causal_budget"])
    _keys("causal_budget", budget, {"length_unit_L0", "nominal_outer_radius", "measurement_radius_maximum", "final_time", "minimum_remaining_buffer", "trial_stage_rule", "strict_acceptance", "failed_trial_does_not_mutate_last_accepted_ledger"})
    budget_values = {name: _fraction(f"causal_budget.{name}", budget[name], positive=name != "measurement_radius_maximum", nonnegative=name == "measurement_radius_maximum") for name in ("length_unit_L0", "nominal_outer_radius", "measurement_radius_maximum", "final_time", "minimum_remaining_buffer")}
    if (
        budget_values != {"length_unit_L0": Q(4), "nominal_outer_radius": Q(128), "measurement_radius_maximum": Q(24), "final_time": Q(32), "minimum_remaining_buffer": Q(16)}
        or budget["trial_stage_rule"] != "debit_maximum_of_previous_endpoint_all_internal_RK_stages_and_explicit_candidate_endpoint_speed_envelopes"
        or budget["strict_acceptance"] != "remaining_buffer_greater_than_minimum"
        or budget["failed_trial_does_not_mutate_last_accepted_ledger"] is not True
    ):
        raise ValueError("BND2 causal budget differs")

    discretization = _mapping("discretization", raw["discretization"])
    _keys("discretization", discretization, {"nominal_point_counts", "primary_SBP_stencil_reach_intervals", "comparator_SBP_stencil_reach_intervals", "stencil_padding_is_not_a_strict_semidiscrete_graph_domain_theorem"})
    if discretization != {"nominal_point_counts": [1025, 2049, 4097], "primary_SBP_stencil_reach_intervals": 3, "comparator_SBP_stencil_reach_intervals": 1, "stencil_padding_is_not_a_strict_semidiscrete_graph_domain_theorem": True}:
        raise ValueError("BND2 discretization differs")

    move = _mapping("outer_boundary_move", raw["outer_boundary_move"])
    _keys("outer_boundary_move", move, {"outer_radii", "preserve_nominal_interior_grid_spacing", "reference_lapse", "reference_radial_metric", "reference_shift", "actual_solution_invariance_check_deferred_to_NUM1_and_ROB1"})
    radii = tuple(_fraction(f"outer_boundary_move.outer_radii[{i}]", value, positive=True) for i, value in enumerate(move["outer_radii"]))
    reference = tuple(_fraction(f"outer_boundary_move.{name}", move[name], positive=name != "reference_shift", nonnegative=name == "reference_shift") for name in ("reference_lapse", "reference_radial_metric", "reference_shift"))
    if radii != (Q(96), Q(128), Q(160)) or reference != (Q(1), Q(1), Q(0)) or move["preserve_nominal_interior_grid_spacing"] is not True or move["actual_solution_invariance_check_deferred_to_NUM1_and_ROB1"] is not True:
        raise ValueError("BND2 outer-boundary move contract differs")

    controls = _mapping("constraint_control", raw["constraint_control"])
    _keys("constraint_control", controls, {"BND1_incoming_main_modes_each_annular_end", "CON4_incoming_gauge_components_each_annular_end", "CON4_reduction_constraints_have_zero_speed", "measurement_region_control_is_by_boundary_free_domain_of_dependence", "nonlinear_constraint_preserving_boundary_map_claimed"})
    if controls != {"BND1_incoming_main_modes_each_annular_end": 6, "CON4_incoming_gauge_components_each_annular_end": 2, "CON4_reduction_constraints_have_zero_speed": True, "measurement_region_control_is_by_boundary_free_domain_of_dependence": True, "nonlinear_constraint_preserving_boundary_map_claimed": False}:
        raise ValueError("BND2 constraint-control contract differs")
    proof = _mapping("proof_contract", raw["proof_contract"])
    _keys(
        "proof_contract",
        proof,
        {
            "HYP2_orthonormal_roots_must_be_converted_to_coordinate_speeds",
            "physical_and_both_auxiliary_cones_must_be_included",
            "evolving_coordinate_speed_envelope_must_be_executable",
            "candidate_endpoint_speed_must_not_be_substituted_by_last_RK_stage",
            "SBP_stencil_reach_must_be_debited",
            "outer_boundary_variants_must_align_at_fixed_interior_spacing",
            "reference_budget_must_pass_at_every_grid_and_boundary_variant",
            "typed_injected_buffer_failure_must_precede_acceptance",
            "incoming_constraints_controlled_only_inside_the_causally_isolated_measurement_region",
            "canonical_hash_bound_result_required",
            "no_holdout_execution",
        },
    )
    if not proof or any(value is not True for value in proof.values()):
        raise ValueError("every BND2 proof contract flag must be true")
    claims = _mapping("claims", raw["claims"])
    if claims != FUTURE_COMMON_NONCLAIMS:
        raise ValueError("BND2 claims must remain fail-closed")

    protocol_raw, protocol = _load_protocol(paths["protocol_config"])
    if protocol["artifact_id"] != SF1_PROTOCOL_ARTIFACT_ID:
        raise ValueError("BND2 protocol differs")
    results: dict[str, dict[str, Any]] = {}
    for name, (artifact_id, gate) in RESULT_REQUIREMENTS.items():
        result = _load_json(paths[name], name)
        if result.get("artifact_id") != artifact_id or result.get("gate_status", {}).get(gate) is not True:
            raise ValueError(f"{name} identity or required gate differs")
        results[name] = result
    return {"raw": raw, "paths": paths, "protocol_raw": protocol_raw, "protocol": protocol, "results": results, "budget": budget_values, "outer_radii": radii, "cone_values": cone_values}


def _scope_bindings(loaded: Mapping[str, Any]) -> dict[str, Any]:
    paths = loaded["paths"]
    protocol = loaded["protocol"]
    return {
        "run1_artifact_id": RUN1_SYM1_ARTIFACT_ID,
        "run1_config_sha256": _sha(paths["run1_config"]),
        "protocol_artifact_id": SF1_PROTOCOL_ARTIFACT_ID,
        "protocol_config_sha256": _sha(paths["protocol_config"]),
        "protocol_semantic_holdout_contract_sha256": protocol["semantic_holdout_contract_sha256"],
        "action_artifact_id": "FGC-1-ACT1",
        "action_config_sha256": _sha(paths["action_config"]),
        "action_result_sha256": _sha(paths["action_result"]),
        "variation_artifact_id": "FGC-1-VAR1",
        "variation_config_sha256": _sha(paths["variation_config"]),
        "variation_result_sha256": _sha(paths["variation_result"]),
        "target_branch": "FGC-QR",
        "physical_equations": "unredefined_ACT1_VAR1",
        "declared_run_envelope_sha256": protocol["semantic_holdout_contract_sha256"],
    }


def _exact_grid_records(loaded: Mapping[str, Any]) -> list[dict[str, Any]]:
    budget = loaded["budget"]
    raw = loaded["raw"]
    grids = aligned_outer_boundary_grids(
        nominal_outer_radius=budget["nominal_outer_radius"],
        nominal_point_counts=raw["discretization"]["nominal_point_counts"],
        outer_radii=loaded["outer_radii"],
        measurement_radius=budget["measurement_radius_maximum"],
        stencil_reach_intervals=raw["discretization"]["primary_SBP_stencil_reach_intervals"],
    )
    records = []
    for grid in grids:
        reference = exact_reference_budget(
            grid,
            final_time=budget["final_time"],
            minimum_causal_buffer=budget["minimum_remaining_buffer"],
            coordinate_speed_bound=exact_coordinate_speed_upper_bound(lapse=Q(1), radial_metric=Q(1), shift=Q(0)),
        )
        records.append({**asdict(grid), "reference_budget": reference})
    if not all(record["reference_budget"]["passed"] for record in records):
        raise ValueError("a BND2 outer-boundary reference budget failed")
    return records


def _typed_runtime_controls() -> dict[str, Any]:
    geometry = BoundaryGeometry(outer_radius=96.0, measurement_radius=24.0, minimum_causal_buffer=16.0, sbp_stencil_reach=0.375)
    original = CausalBudgetState(previous_speed_upper=1.2)
    accepted = accept_causal_step(
        original,
        geometry,
        trial_time=1.0,
        stage_coordinate_speed_uppers=(1.1, 1.2, 1.15, 1.2),
        candidate_endpoint_coordinate_speed_upper=1.3,
    )
    before_failure = accepted
    try:
        accept_causal_step(
            accepted,
            geometry,
            trial_time=32.0,
            stage_coordinate_speed_uppers=(4.0, 4.0, 4.0, 4.0),
            candidate_endpoint_coordinate_speed_upper=5.0,
        )
    except BoundaryControlStop as stop:
        failure = {
            "typed_reason": stop.reason,
            "trial_time": stop.assessment.candidate_state.accepted_time,
            "candidate_remaining_buffer": stop.assessment.remaining_causal_buffer,
            "candidate_margin": stop.assessment.strict_margin_over_required_buffer,
            "trial_rejected": True,
        }
    else:
        raise ValueError("injected BND2 buffer failure was accepted")
    return {
        "accepted_control": {
            "accepted_time": accepted.accepted_time,
            "accumulated_characteristic_distance": accepted.accumulated_characteristic_distance,
            "accepted_endpoint_speed_upper": accepted.previous_speed_upper,
        },
        "injected_failure": failure,
        "last_accepted_state_unchanged_by_failed_pure_transaction": accepted == before_failure,
    }


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    loaded = load_config(config_path)
    paths = loaded["paths"]
    results = loaded["results"]
    bnd1 = results["frozen_boundary_result"]["uniform_frozen_boundary_certificate"]
    con4_quantitative = results["physical_constraints_result"]["artifact_payload"]["quantitative_evidence"]
    bnd1_inner = bnd1["boundary_operators"]["inner"]["rank"]
    bnd1_outer = bnd1["boundary_operators"]["outer"]["rank"]
    gauge_records = [item["gauge_characteristics"] for item in con4_quantitative["state_controls"]]
    if bnd1_inner != 6 or bnd1_outer != 6 or any(item["inner_incoming_component_count"] != 2 or item["outer_incoming_component_count"] != 2 for item in gauge_records):
        raise ValueError("BND1/CON4 incoming characteristic counts differ")
    if not bnd1["characteristics"]["incoming_counts_verified_from_strict_intervals"]:
        raise ValueError("BND1 incoming classifications are not strict")
    # Spell out the ledger and bind every scientific predecessor config/result.
    predecessor_names = tuple(name for name in EXPECTED_PATHS if name not in {"protocol_config", "run1_config"})
    grids = _exact_grid_records(loaded)
    runtime = _typed_runtime_controls()
    smallest = min(grids, key=lambda item: (item["outer_radius"], -item["spacing"]))
    quantitative = {
        "orthonormal_to_coordinate_cone_map": {
            "formula": "dr/dt=-shift+(lapse/Lambda)z",
            "absolute_envelope": "abs(shift)+(lapse/Lambda)*6/5",
            "physical_local_absolute_bound": LOCAL_CONE_SPEED_BOUNDS["physical"],
            "tilde_local_absolute_bound": LOCAL_CONE_SPEED_BOUNDS["tilde"],
            "hat_local_absolute_bound": LOCAL_CONE_SPEED_BOUNDS["hat"],
            "evolving_metric_values_required_at_runtime": True,
        },
        "causal_budget_contract": {
            "formula": "r_outer-r_measure-stencil_reach-integral_v_coordinate_max_dt",
            "minimum_remaining_buffer": loaded["budget"]["minimum_remaining_buffer"],
            "strict_inequality_required": True,
            "trial_stage_maximum_is_debited_before_acceptance": True,
            "candidate_endpoint_speed_is_evaluated_explicitly": True,
            "runtime_controls": runtime,
        },
        "outer_boundary_move_control": {
            "fixed_interior_spacing_at_each_resolution": True,
            "grid_records": grids,
            "smallest_reference_remaining_buffer": smallest["reference_budget"]["remaining_causal_buffer"],
            "smallest_reference_margin_over_required_buffer": smallest["reference_budget"]["strict_margin_over_required_buffer"],
            "actual_evolved_observable_move_agreement_deferred_to_NUM1_and_ROB1": True,
        },
        "incoming_constraint_control": {
            "BND1_main_incoming_rank_inner": bnd1_inner,
            "BND1_main_incoming_rank_outer": bnd1_outer,
            "CON4_gauge_incoming_components_inner": 2,
            "CON4_gauge_incoming_components_outer": 2,
            "CON4_reduction_constraints_zero_speed": True,
            "control_mechanism": "exclude_boundary_domain_of_dependence_from_retained_measurement_interval",
            "nonlinear_constraint_preserving_boundary_map_proved": False,
        },
        "epistemic_boundary": {
            "stencil_pad_is_not_exact_semidiscrete_compact_support": True,
            "complete_nonlinear_constraint_preserving_IBVP_proved": False,
            "trajectory_speed_envelope_is_deferred_to_HLT1": True,
            "actual_boundary_move_solution_agreement_is_deferred_to_NUM1_and_ROB1": True,
            "time_evolution_performed": False,
            "FGCQR_holdout_outcome_inspected": False,
        },
    }
    scope = _scope_bindings(loaded)
    return _serial({
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": FUTURE_EVIDENCE_CLASSIFICATIONS["boundary_control"],
        "generated_by": "scripts/reproduce_fgc_bnd2_cp1.py",
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): _sha(config_path)},
        "predecessor_sha256": {_rel(paths[name]): _sha(paths[name]) for name in predecessor_names},
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "scope_bindings": scope,
        "certificate_contract": {
            "artifact_specific_payload_validated": True,
            "canonical_reproduction_passed": True,
            "scope_bindings_verified": True,
            "outcome_data_not_used_to_select_certificate_contract": True,
        },
        "artifact_payload": {
            "run_envelope_semantic_sha256": scope["declared_run_envelope_sha256"],
            "evolving_physical_and_auxiliary_cones_included": True,
            "SBP_stencil_reach_included": True,
            "measured_affine_interval_outside_boundary_domain_of_dependence": True,
            "outer_boundary_move_check_passed": True,
            "incoming_constraint_characteristics_controlled": True,
            "quantitative_evidence": quantitative,
        },
        "gate_status": {REQUIRED_GATE: True},
        "nonclaims": dict(FUTURE_COMMON_NONCLAIMS),
    })


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    value = _load_json(path, "BND2 result")
    if path.read_text(encoding="utf-8") != _canonical(value):
        raise ValueError("BND2 result must use canonical sorted JSON")
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    arguments.output.write_text(_canonical(record(arguments.config)), encoding="utf-8")


if __name__ == "__main__":
    main()
