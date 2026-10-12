#!/usr/bin/env python3
"""Reproduce the outcome-blind HLT7 spectral-conditioning diagnosis."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-cal6-pref8.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-cal6-pref8.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-cal6-pref8.md"

from scripts.reproduce_fgc_hlt7_mon7 import Q  # noqa: E402
from recursive_horizons.fgc.evolution.cal4_common_event_diagnosis import (  # noqa: E402
    finest_pair_nested_spectral_admission,
)
from recursive_horizons.fgc.evolution.calibration_runtime import (  # noqa: E402
    project_gr0_semidiscrete_state,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.evolution.health_monitor import (  # noqa: E402
    compact_vacuum_buffer_window,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    PROTO4_SPECTRAL_FIELD_ORDER,
    SpectralThresholds,
    proper_radial_profile,
    spectral_field_budgets,
)
from recursive_horizons.fgc.evolution.spectral_sensitivity import (  # noqa: E402
    resolved_or_saturated_admission,
    spectral_tail_sensitivity,
    three_grid_profile_convergence,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


ARTIFACT_ID = "FGC-1-CAL6-PREF8"
PROJECT_VERSION = "0.11.0"
EXPECTED_PATHS = {
    "runtime_authorization_result": "results/fgc-1-hlt7-mon7.json",
    "runtime_authorization_config": "configs/fgc/fgc-1-hlt7-mon7.toml",
    "run_plan_config": "configs/fgc/fgc-1-cal6-run1.toml",
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v9.toml",
}
EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO9",
    "role": "pre_trajectory_proper_spectral_ratio_conditioning_diagnosis",
    "calibration_branch": "GR-0",
    "PROTO9_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "mechanism_question_answered": False,
}
EXPECTED_IMMUTABLE = {
    "checkpoint_commit": "234c5a5cc2c0b6bdc16d22146cdb9a9c6ae5b98d",
    "runtime_authorization_result_sha256": "e6eb5832302eea39f7adf7fee2996cf659979ebbc58317b8f6d3b7e092f789c5",
    "runtime_authorization_config_sha256": "9cf667038cbddcad2c92b1728cf439e7131ec2577d939a13d32868538eef7bfd",
    "run_plan_config_sha256": "ed0daf5be8353274131ad20ea13534fed491aa238a4234f07648e2a3bd6cafb3",
    "protocol_config_sha256": "81ce8426b006518d58accae50e43fcbfe942d013dfd010a716ff265b57afb48e",
    "all_predecessor_blobs_must_be_read_from_checkpoint": True,
    "checkpoint_must_be_ancestor_of_HEAD": True,
}
EXPECTED_DIAGNOSTIC = {
    "point_counts": [2049, 4097, 8193],
    "field_order": list(PROTO4_SPECTRAL_FIELD_ORDER),
    "unchanged_maximum_nested_tail_ratio": "1/4",
    "unchanged_maximum_top_band_field_power_fraction": "1/1048576",
    "unchanged_maximum_top_band_derivative_power_fraction": "1/1024",
    "coarse_to_medium_field_and_derivative_ratios_must_pass_directly": True,
    "medium_and_fine_individual_absolute_budgets_must_pass": True,
    "failed_finest_pair_ratio_may_only_be_classified_as_diagnostically_saturated": True,
    "medium_and_fine_top_band_erasure_must_fit_inside_each_measured_round_trip_scale": True,
    "all_three_round_trip_residuals_must_contract_strictly_or_be_exact_zero": True,
    "whole_profile_infinity_and_RMS_differences_must_contract_strictly_or_be_exact_zero": True,
    "all_raw_power_fractions_ratios_residuals_and_erasure_witnesses_must_be_serialized": True,
    "tail_erasure_is_an_orthogonal_FFT_diagnostic_input_witness": True,
    "round_trip_interpolation_diagnostic_is_not_a_continuum_error_bound": True,
    "diagnostic_saturation_is_not_physical_resolution_or_mechanism_evidence": True,
    "no_epsilon_floor_or_new_numeric_tolerance_is_permitted": True,
}
EXPECTED_REVISION = {
    "new_protocol_required": "FGC-2-SF1-PROTO10",
    "only_finest_pair_nested_tail_interpretation_may_change": True,
    "resolution_ladder_unchanged": True,
    "physical_inputs_equations_methods_CFL_and_stops_unchanged": True,
    "source_solver_raw_residual_threshold_and_retry_budget_unchanged": True,
    "constraint_ownership_magnitude_and_order_rules_unchanged": True,
    "absolute_spectral_budgets_and_direct_ratio_ceiling_unchanged": True,
    "coarse_to_medium_direct_ratio_veto_unchanged": True,
    "finest_pair_direct_ratio_remains_the_primary_route": True,
    "saturation_route_must_add_profile_convergence_and_map_conditioning_guards": True,
    "fresh_PROTO10_namespace_required": True,
    "separate_freeze_and_runtime_certificate_required": True,
}
EXPECTED_PROOF_KEYS = {
    "immutable_HLT7_must_remain_canonical_and_pre_trajectory",
    "all_twelve_projected_state_hashes_must_match_HLT7",
    "all_four_t0_common_events_must_reconstruct_from_direct_states",
    "HLT7_raw_ratios_and_failed_admissions_must_reproduce_exactly",
    "all_medium_and_fine_absolute_budgets_must_pass_unchanged",
    "all_coarse_to_medium_field_and_derivative_tail_ratios_must_pass_directly",
    "every_nonzero_complete_profile_must_contract_in_infinity_and_RMS_norm",
    "every_nonzero_round_trip_residual_must_contract_on_both_adjacent_pairs",
    "every_failed_finest_pair_ratio_must_have_a_top_band_erasure_witness_below_the_measured_map_non_idempotence",
    "prospective_resolved_or_saturated_rule_must_pass_all_four_t0_cases",
    "no_trajectory_namespace_may_be_created_or_read",
    "canonical_hash_bound_result_required",
}
EXPECTED_CLAIMS = {
    "PROTO9_preflight_conditioning_diagnosis_completed": True,
    "PROTO9_ratio_only_failure_is_candidate_or_gradient_evidence": False,
    "PROTO10_spectral_contract_revision_required": True,
    "PROTO10_frozen": False,
    "PROTO10_successor_runtime_implemented": False,
    "PROTO10_fresh_GR0_dynamic_calibration_authorized": False,
    "classical_spherical_diagnostic_authorized": False,
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
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto4_admission.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/cal4_common_event_diagnosis.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto8_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto9_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/spectral_sensitivity.py",
    Path(__file__).resolve(),
)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _digest(value: object) -> str:
    return sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()


def _git_blob(commit: str, relative: str) -> bytes:
    completed = subprocess.run(
        ["git", "show", f"{commit}:{relative}"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise ValueError(f"immutable CAL6/PREF8 predecessor is unavailable: {relative}")
    return completed.stdout


def _json_bytes(blob: bytes, name: str) -> dict[str, Any]:
    duplicates: list[str] = []

    def hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        answer: dict[str, Any] = {}
        for key, value in pairs:
            if key in answer:
                duplicates.append(key)
            answer[key] = value
        return answer

    value = json.loads(blob.decode("utf-8"), object_pairs_hook=hook)
    if duplicates or not isinstance(value, dict):
        raise ValueError(f"{name} is noncanonical JSON")
    return value


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    expected_root = {
        "schema_version",
        "artifact_id",
        "project_version",
        "metric_signature",
        "riemann_convention",
        *EXPECTED_PATHS,
        "scope",
        "immutable_preflight",
        "diagnostic_definition",
        "prospective_revision",
        "proof_contract",
        "claims",
    }
    if set(raw) != expected_root:
        raise ValueError("CAL6/PREF8 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("CAL6/PREF8 identity or convention differs")
    for name, expected in EXPECTED_PATHS.items():
        if raw[name] != expected or not (REPOSITORY / expected).is_file():
            raise ValueError(f"CAL6/PREF8 {name} differs")
    if raw["scope"] != EXPECTED_SCOPE:
        raise ValueError("CAL6/PREF8 scope differs")
    if raw["immutable_preflight"] != EXPECTED_IMMUTABLE:
        raise ValueError("CAL6/PREF8 immutable preflight differs")
    if raw["diagnostic_definition"] != EXPECTED_DIAGNOSTIC:
        raise ValueError("CAL6/PREF8 diagnostic definition differs")
    if raw["prospective_revision"] != EXPECTED_REVISION:
        raise ValueError("CAL6/PREF8 prospective revision differs")
    if set(raw["proof_contract"]) != EXPECTED_PROOF_KEYS or any(
        value is not True for value in raw["proof_contract"].values()
    ):
        raise ValueError("CAL6/PREF8 proof contract differs")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("CAL6/PREF8 claims differ")
    return raw


def _immutable_preflight(raw: Mapping[str, Any]) -> dict[str, Any]:
    immutable = raw["immutable_preflight"]
    commit = immutable["checkpoint_commit"]
    if subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    ).returncode != 0:
        raise ValueError("CAL6/PREF8 checkpoint is not an ancestor of HEAD")
    blobs: dict[str, bytes] = {}
    for name, relative in EXPECTED_PATHS.items():
        blob = _git_blob(commit, relative)
        if sha256(blob).hexdigest() != immutable[f"{name}_sha256"]:
            raise ValueError(f"immutable CAL6/PREF8 {name} hash differs")
        if _sha(REPOSITORY / relative) != immutable[f"{name}_sha256"]:
            raise ValueError(f"current CAL6/PREF8 {name} drifted")
        blobs[name] = blob
    hlt7 = _json_bytes(blobs["runtime_authorization_result"], "immutable HLT7")
    payload = hlt7.get("artifact_payload", {})
    events = payload.get("initial_common_events", [])
    inputs = payload.get("frozen_run_inputs", [])
    if (
        hlt7.get("artifact_id") != "FGC-1-HLT7-MON7"
        or hlt7.get("classification")
        != "pre_trajectory_PROTO9_resolution_ladder_runtime_with_t0_nested_tail_obstruction"
        or hlt7.get("gate_status", {}).get("PROTO9_successor_runtime_implemented")
        is not True
        or hlt7.get("gate_status", {}).get(
            "PROTO9_fresh_GR0_dynamic_calibration_authorized"
        )
        is not False
        or len(inputs) != 12
        or len(events) != 4
        or any(event.get("trajectory_advanced") is not False for event in events)
        or any(event.get("admission_passed") is not False for event in events)
        or any(value is not False for value in hlt7.get("nonclaims", {}).values())
    ):
        raise ValueError("immutable HLT7 preflight boundary differs")
    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "predecessor_sha256": {
            EXPECTED_PATHS[name]: immutable[f"{name}_sha256"]
            for name in EXPECTED_PATHS
        },
        "hlt7": hlt7,
        "run_plan": tomllib.loads(blobs["run_plan_config"].decode("utf-8")),
    }


def _budget_payload(budget: object) -> dict[str, Any]:
    return asdict(budget)  # type: ignore[arg-type]


def _sensitivity_payload(witness: object) -> dict[str, Any]:
    return asdict(witness)  # type: ignore[arg-type]


def _case_diagnosis(
    amplitude: str,
    method: str,
    plan: Mapping[str, Any],
    hlt7_inputs: Mapping[tuple[str, str, int], Mapping[str, Any]],
    hlt7_event: Mapping[str, Any],
) -> dict[str, Any]:
    physical = plan["physical_inputs"]
    numerics = plan["numerics"]
    order = 4 if method == "RK4" else 2
    point_counts = (2049, 4097, 8193)
    profiles = []
    budgets = []
    sensitivities = []
    grid_payloads = []
    for point_count in point_counts:
        initial = construct_gr0_grid_initial_data(
            PulseParameters(
                chi_amplitude=float(Q(amplitude)),
                center=float(Q(physical["chi_center"])),
                half_width=float(Q(physical["chi_half_width"])),
                phi_amplitude=float(Q(physical["phi_seed_amplitude"])),
                planck_mass=float(Q(physical["planck_mass"])),
                scalar_mass=float(Q(physical["scalar_mass"])),
                quartic_coupling=float(Q(physical["quartic_coupling"])),
            ),
            point_count=point_count,
            outer_radius=float(Q(physical["outer_radius"])),
            constraint_method=method,
            diagnostic_spatial_order=order,
        )
        state = project_gr0_semidiscrete_state(initial, spatial_order=order)
        state_hash = array_content_sha256(state.u, state.p, state.q)
        expected_hash = hlt7_inputs[(amplitude, method, point_count)][
            "projected_state_sha256"
        ]
        if state_hash != expected_hash:
            raise ValueError("CAL6/PREF8 reconstructed state differs from HLT7")
        profile = proper_radial_profile(
            state.u,
            state.q,
            initial.grid.coordinates,
            cutoff=float(Q(physical["cutoff_Lambda"])),
            measurement_radius_maximum=float(Q(physical["measurement_radius_maximum"])),
        )
        proper = profile.proper_coordinates
        window = compact_vacuum_buffer_window(
            proper,
            support_minimum=float(proper[0]),
            support_maximum=float(proper[-1]),
            taper_width=float(
                (proper[-1] - proper[0])
                * float(Q(numerics["proper_spectral_taper_fraction"]))
            ),
        )
        budget = dict(
            spectral_field_budgets(
                profile.deviations,
                proper,
                cutoff=float(Q(physical["cutoff_Lambda"])),
                window=window,
            )
        )
        sensitivity = {
            name: spectral_tail_sensitivity(
                profile.deviations[:, index],
                proper,
                window=window,
                round_trip_interpolation_infinity=(
                    profile.round_trip_interpolation_infinity[index]
                ),
                budget=budget[name],
            )
            for index, name in enumerate(PROTO4_SPECTRAL_FIELD_ORDER)
        }
        profiles.append(profile)
        budgets.append(budget)
        sensitivities.append(sensitivity)
        grid_payloads.append(
            {
                "point_count": point_count,
                "spectral_sample_count": int(proper.size),
                "projected_state_sha256": state_hash,
                "round_trip_interpolation_infinity_by_field": {
                    name: float(profile.round_trip_interpolation_infinity[index])
                    for index, name in enumerate(PROTO4_SPECTRAL_FIELD_ORDER)
                },
                "spectral_budgets_by_field": {
                    name: _budget_payload(budget[name])
                    for name in PROTO4_SPECTRAL_FIELD_ORDER
                },
                "tail_sensitivity_by_field": {
                    name: _sensitivity_payload(sensitivity[name])
                    for name in PROTO4_SPECTRAL_FIELD_ORDER
                },
            }
        )

    raw_admission = finest_pair_nested_spectral_admission(point_counts, budgets)
    convergence = three_grid_profile_convergence(point_counts, profiles)
    successor = resolved_or_saturated_admission(
        point_counts,
        budgets,
        sensitivities,
        convergence,
    )
    raw_field = {name: list(values) for name, values in raw_admission.field_power_tail_ratios.items()}
    raw_derivative = {
        name: list(values)
        for name, values in raw_admission.derivative_power_tail_ratios.items()
    }
    if (
        raw_field != hlt7_event["field_power_tail_ratios"]
        or raw_derivative != hlt7_event["derivative_power_tail_ratios"]
        or raw_admission.admission_passed is not False
        or hlt7_event["admission_passed"] is not False
    ):
        raise ValueError("CAL6/PREF8 did not reproduce the HLT7 raw obstruction")

    adjacent = [asdict(item) for item in convergence.adjacent_differences]
    return {
        "amplitude": amplitude,
        "method": method,
        "coordinate_time": 0.0,
        "point_counts": list(point_counts),
        "grid_diagnostics": grid_payloads,
        "raw_PROTO9_admission": {
            "field_power_tail_ratios": raw_field,
            "derivative_power_tail_ratios": raw_derivative,
            "finest_pair_individual_budgets_passed": (
                raw_admission.finest_pair_individual_budgets_passed
            ),
            "every_nested_tail_ratio_passed": (
                raw_admission.every_nested_tail_ratio_passed
            ),
            "admission_passed": raw_admission.admission_passed,
        },
        "whole_profile_convergence": {
            "adjacent_differences": adjacent,
            "exact_zero_by_field": dict(convergence.exact_zero_by_field),
            "infinity_contracted_by_field": dict(
                convergence.infinity_contracted_by_field
            ),
            "rms_contracted_by_field": dict(convergence.rms_contracted_by_field),
            "round_trip_contracted_by_field": dict(
                convergence.round_trip_contracted_by_field
            ),
            "complete_profile_contraction_by_field": dict(
                convergence.complete_profile_contraction_by_field
            ),
            "every_field_contracted": convergence.every_field_contracted,
        },
        "prospective_resolved_or_saturated_admission": {
            "field_power_tail_ratios": {
                name: list(values)
                for name, values in successor.field_power_tail_ratios.items()
            },
            "derivative_power_tail_ratios": {
                name: list(values)
                for name, values in successor.derivative_power_tail_ratios.items()
            },
            "coarse_to_medium_direct_passed_by_field": dict(
                successor.coarse_to_medium_direct_passed_by_field
            ),
            "coarse_to_medium_direct_passed_by_derivative": dict(
                successor.coarse_to_medium_direct_passed_by_derivative
            ),
            "finest_pair_field_classification": dict(
                successor.finest_pair_field_classification
            ),
            "finest_pair_derivative_classification": dict(
                successor.finest_pair_derivative_classification
            ),
            "finest_pair_individual_budgets_passed": (
                successor.finest_pair_individual_budgets_passed
            ),
            "every_profile_contracted": successor.every_profile_contracted,
            "saturation_used": successor.saturation_used,
            "every_tail_resolved_or_saturated": (
                successor.every_tail_resolved_or_saturated
            ),
            "admission_passed": successor.admission_passed,
        },
        "trajectory_read": False,
    }


def _aggregate(cases: list[Mapping[str, Any]]) -> dict[str, Any]:
    failed_raw = 0
    saturated = 0
    direct_finest = 0
    erasure_ratios = []
    field_headrooms = []
    derivative_headrooms = []
    profile_contraction_ratios = []
    round_trip_contraction_ratios = []
    limits = SpectralThresholds()
    for case in cases:
        raw = case["raw_PROTO9_admission"]
        successor = case["prospective_resolved_or_saturated_admission"]
        for kind in ("field", "derivative"):
            ratios = raw[f"{kind}_power_tail_ratios"]
            classes = successor[f"finest_pair_{kind}_classification"]
            for name in PROTO4_SPECTRAL_FIELD_ORDER:
                if ratios[name][1] is None or ratios[name][1] >= limits.maximum_nested_tail_ratio:
                    failed_raw += 1
                if classes[name] == "diagnostically_saturated":
                    saturated += 1
                elif classes[name] == "directly_resolved":
                    direct_finest += 1
        grids = case["grid_diagnostics"]
        for grid in grids[1:]:
            for name in PROTO4_SPECTRAL_FIELD_ORDER:
                budget = grid["spectral_budgets_by_field"][name]
                if budget["top_band_field_power_fraction"] > 0.0:
                    field_headrooms.append(
                        limits.maximum_top_band_field_power_fraction
                        / budget["top_band_field_power_fraction"]
                    )
                if budget["top_band_derivative_weighted_power_fraction"] > 0.0:
                    derivative_headrooms.append(
                        limits.maximum_top_band_derivative_power_fraction
                        / budget["top_band_derivative_weighted_power_fraction"]
                    )
                witness = grid["tail_sensitivity_by_field"][name]
                if witness["erasure_over_round_trip"] is not None:
                    erasure_ratios.append(witness["erasure_over_round_trip"])
        adjacent = case["whole_profile_convergence"]["adjacent_differences"]
        for name in PROTO4_SPECTRAL_FIELD_ORDER:
            for norm in ("infinity_by_field", "rms_by_field"):
                coarse = adjacent[0][norm][name]
                fine = adjacent[1][norm][name]
                if coarse > 0.0:
                    profile_contraction_ratios.append(fine / coarse)
            rt = [
                grid["round_trip_interpolation_infinity_by_field"][name]
                for grid in grids
            ]
            for denominator, numerator in zip(rt, rt[1:]):
                if denominator > 0.0:
                    round_trip_contraction_ratios.append(numerator / denominator)
    return {
        "case_count": len(cases),
        "raw_failed_finest_pair_metric_count": failed_raw,
        "prospective_diagnostically_saturated_metric_count": saturated,
        "prospective_directly_resolved_finest_pair_metric_count": direct_finest,
        "minimum_medium_or_fine_field_absolute_budget_headroom": min(field_headrooms),
        "minimum_medium_or_fine_derivative_absolute_budget_headroom": min(
            derivative_headrooms
        ),
        "maximum_top_band_erasure_over_round_trip": max(erasure_ratios),
        "maximum_complete_profile_finest_pair_difference_ratio": max(
            profile_contraction_ratios
        ),
        "maximum_round_trip_adjacent_ratio": max(round_trip_contraction_ratios),
        "all_raw_PROTO9_admissions_failed": all(
            case["raw_PROTO9_admission"]["admission_passed"] is False
            for case in cases
        ),
        "all_prospective_resolved_or_saturated_admissions_passed": all(
            case["prospective_resolved_or_saturated_admission"]["admission_passed"]
            is True
            for case in cases
        ),
    }


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    raw = load_config(config_path)
    immutable = _immutable_preflight(raw)
    hlt7 = immutable["hlt7"]
    inputs = {
        (item["amplitude"], item["method"], item["point_count"]): item
        for item in hlt7["artifact_payload"]["frozen_run_inputs"]
    }
    events = {
        (item["amplitude"], item["method"]): item
        for item in hlt7["artifact_payload"]["initial_common_events"]
    }
    cases = [
        _case_diagnosis(
            amplitude,
            method,
            immutable["run_plan"],
            inputs,
            events[(amplitude, method)],
        )
        for amplitude in ("5/2", "3")
        for method in ("RK4", "SSPRK3")
    ]
    aggregate = _aggregate(cases)
    if (
        aggregate["raw_failed_finest_pair_metric_count"] != 16
        or aggregate["prospective_diagnostically_saturated_metric_count"] != 16
        or aggregate["prospective_directly_resolved_finest_pair_metric_count"] != 32
        or aggregate["all_raw_PROTO9_admissions_failed"] is not True
        or aggregate["all_prospective_resolved_or_saturated_admissions_passed"]
        is not True
    ):
        raise ValueError("CAL6/PREF8 conditioning diagnosis no longer has its frozen shape")
    for root in (
        "runs/fgc-2-sf1/proto9/calibration",
        "runs/fgc-2-sf1/proto9/holdout",
    ):
        path = REPOSITORY / root
        if path.exists() and (not path.is_dir() or any(path.iterdir())):
            raise ValueError("CAL6/PREF8 would cross a trajectory namespace boundary")

    payload = {
        "immutable_preflight": {
            "checkpoint_commit": immutable["checkpoint_commit"],
            "checkpoint_is_ancestor_of_HEAD": True,
            "artifact_id": hlt7["artifact_id"],
            "classification": hlt7["classification"],
            "PROTO9_successor_runtime_implemented": True,
            "PROTO9_fresh_GR0_dynamic_calibration_authorized": False,
            "trajectory_read": False,
        },
        "diagnostic_contract": dict(raw["diagnostic_definition"]),
        "case_diagnostics": cases,
        "aggregate": aggregate,
        "decision": {
            "PROTO9_ratio_only_failure_is_candidate_or_gradient_evidence": False,
            "HLT7_finest_pair_ratio_is_conditioned_below_the_diagnostic_map_scale": True,
            "smallest_discriminating_successor": (
                "retain_the_direct_coarse_ratio_and_all_absolute_budgets_while_adding_"
                "a_profile_convergence_guarded_finest_pair_saturation_route"
            ),
            "new_protocol_required": raw["prospective_revision"][
                "new_protocol_required"
            ],
            "current_campaign_authorized": False,
        },
        "prospective_revision": dict(raw["prospective_revision"]),
        "epistemic_boundary": {
            "round_trip_is_a_continuum_error_bound": False,
            "tail_erasure_is_a_physical_field_perturbation": False,
            "diagnostic_saturation_proves_physical_resolution": False,
            "trajectory_or_candidate_outcome_read": False,
            "mechanism_question_answered": False,
        },
    }
    public_immutable = dict(immutable["predecessor_sha256"])
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "generated_by": "scripts/reproduce_fgc_cal6_pref8.py",
        "classification": "pre_trajectory_proper_spectral_ratio_conditioning_diagnosis",
        "source_config_sha256": {
            str(config_path.relative_to(REPOSITORY)): _sha(config_path)
        },
        "derivation_document": str(OWNER_DOCUMENT.relative_to(REPOSITORY)),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "predecessor_sha256": public_immutable,
        "implementation_sha256": {
            str(path.relative_to(REPOSITORY)): _sha(path) for path in IMPLEMENTATION
        },
        "artifact_payload": payload,
        "gate_status": dict(raw["claims"]),
        "nonclaims": {
            "fresh_GR0_recalibration_completed": False,
            "PROTO9_trajectory_opened": False,
            "SGBL_comparison_completed": False,
            "FGCQR_trajectory_opened": False,
            "trapped_sphere_formed_dynamically": False,
            "positive_Raychaudhuri_margin_measured": False,
            "retained_EFT_validity": False,
            "physical_transition": False,
            "singularity_resolution": False,
            "child_domain_or_topology": False,
            "dark_sector_mechanism": False,
            "varying_locally_measured_c": False,
        },
    }


def verify_canonical(path: Path, config_path: Path = DEFAULT_CONFIG) -> None:
    canonical = _json_bytes(path.read_bytes(), "canonical CAL6/PREF8")
    fresh = record(config_path)
    if canonical != fresh:
        raise ValueError("CAL6/PREF8 output differs from a fresh reproduction")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    config = arguments.config.resolve()
    output = arguments.output.resolve()
    if arguments.check:
        verify_canonical(output, config)
        print(f"verified {output}")
        return
    result = record(config)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(result, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )
    print(
        f"wrote {output} "
        f"(raw failures={result['artifact_payload']['aggregate']['raw_failed_finest_pair_metric_count']}; "
        "PROTO10 frozen=false; trajectory_read=false)"
    )


if __name__ == "__main__":
    main()
