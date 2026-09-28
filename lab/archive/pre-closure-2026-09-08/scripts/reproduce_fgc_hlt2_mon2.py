#!/usr/bin/env python3
"""Regenerate the outcome-neutral FGC-1-HLT2-MON2 certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite, pi
from pathlib import Path
import sys
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hlt2-mon2.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hlt2-mon2.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hlt2-mon2.md"

from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.evolution.health_monitor import (  # noqa: E402
    compact_vacuum_buffer_window,
)
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    FGCQR_BRANCH_DEFINITION_OWNER,
    PROTO4_BRANCH_ORDER,
    PROTO4_CONSTRAINT_ORDER,
    PROTO4_SPECTRAL_FIELD_ORDER,
    ConstraintEventSample,
    SpectralPowerBudget,
    branch_stop_applicability,
    causal_past_window,
    common_event_constraint_admission,
    conservative_observable_error_sum,
    conservative_richardson_error,
    nested_spectral_admission,
    normalized_constraint_norms,
    proper_radial_profile,
    spectral_field_budgets,
    windowed_spectral_power_budget,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)
from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    PROTO4_CANDIDATE_BRANCH_STOP_IDS,
    PROTO4_REPLACED_STOP_IDS,
    PROTO4_UNIVERSAL_STOP_IDS,
    SF1_PROTOCOL_V4_ARTIFACT_ID,
    validate_sf1_protocol,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-HLT2-MON2"
PROJECT_VERSION = "0.11.0"
REQUIRED_GATE = "PROTO4_successor_admission_monitor_implemented"
EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v4.toml",
    "protocol_freeze_config": "configs/fgc/fgc-1-pro4-frz1.toml",
    "protocol_freeze_result": "results/fgc-1-pro4-frz1.json",
    "historical_monitor_config": "configs/fgc/fgc-1-hlt1-mon1.toml",
    "historical_monitor_result": "results/fgc-1-hlt1-mon1.json",
    "diagnosis_config": "configs/fgc/fgc-1-cal0-pref2.toml",
    "diagnosis_result": "results/fgc-1-cal0-pref2.json",
    "nonlinear_source_result": "results/fgc-1-src1-nl1.json",
    "run_domain_result": "results/fgc-1-dom4-run1.json",
    "FGCQR_health_result": "results/fgc-1-hyp2-md1.json",
}
EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO4",
    "monitor_role": "outcome_neutral_branch_applicability_and_convergent_numerical_admission",
    "calibration_branch": "GR-0",
    "comparison_branch": "SGB-L",
    "holdout_branch": "FGC-QR",
    "historical_HLT1_result_preserved": True,
    "PROTO4_physical_study_unchanged": True,
    "fresh_calibration_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_defocusing_outcome_classified": False,
}
EXPECTED_BRANCH = {
    "branch_order": list(PROTO4_BRANCH_ORDER),
    "universal_runtime_stop_ids": list(PROTO4_UNIVERSAL_STOP_IDS),
    "candidate_branch_only_stop_ids": list(PROTO4_CANDIDATE_BRANCH_STOP_IDS),
    "replaced_stop_ids": list(PROTO4_REPLACED_STOP_IDS),
    "GR0_candidate_stops_forbidden": True,
    "FGCQR_definition_owner": FGCQR_BRANCH_DEFINITION_OWNER,
    "FGCQR_definition_available": True,
    "SGBL_definition_owner": "",
    "SGBL_definition_available": False,
    "SGBL_cannot_borrow_FGCQR_definition": True,
    "missing_branch_definition_forbids_that_branch_run": True,
}
EXPECTED_SPECTRAL = {
    "field_order": list(PROTO4_SPECTRAL_FIELD_ORDER),
    "center_limits": "shift_over_r=partial_r_shift_and_R_over_r=partial_r_R_at_exact_r_zero",
    "proper_radial_coordinate": "cumulative_trapezoid_of_positive_lambda_dr",
    "uniform_resampling": "deterministic_piecewise_linear_with_round_trip_diagnostic",
    "spatial_window": "C_infinity_two_edge_compact_window_inside_measurement_interval",
    "temporal_window": "C_infinity_two_edge_window_using_causal_past_samples_only",
    "minimum_spatial_samples": 16,
    "minimum_temporal_samples": 64,
    "unresolved_top_fraction": "1/8",
    "maximum_top_band_field_power_fraction": "1/1048576",
    "maximum_top_band_derivative_weighted_power_fraction": "1/1024",
    "maximum_nested_refinement_tail_ratio": "1/4",
    "maximum_RMS_scale_over_Lambda": "1/4",
    "minimum_nested_resolutions": 3,
    "strict_upper_bounds": True,
    "one_sided_real_FFT_multiplicities_included": True,
    "last_bin_alias_boolean_is_diagnostic_only": True,
    "zero_signal_has_zero_tail_and_RMS_scale": True,
    "interpolation_diagnostic_is_not_a_continuum_error_bound": True,
}
EXPECTED_CONSTRAINT = {
    "field_order": list(PROTO4_CONSTRAINT_ORDER),
    "normalization": "max_sum_absolute_unredefined_terms_M_Pl_squared_over_L0_squared_1e-30_per_component_and_point",
    "coarsest_normalized_guard_max": "1/50",
    "finest_normalized_guard_max": "1/1000",
    "minimum_observed_convergence_order": "3/2",
    "minimum_nested_resolutions": 3,
    "strict_upper_guards": True,
    "common_GR0_event": "bitwise_equal_coordinate_time",
    "common_DEF1_event": "equal_coordinate_time_affine_label_and_affine_origin",
    "exact_zero_pair_passes_without_invented_infinite_order": True,
    "nonzero_after_exact_zero_fails": True,
    "single_resolution_smallness_forbidden": True,
    "numerical_smallness_is_not_a_constraint_propagation_theorem": True,
}
EXPECTED_ERROR = {
    "Richardson_formula": "max_abs_fine_minus_aligned_coarse_divided_by_refinement_ratio_power_order_minus_1",
    "Richardson_requires_common_event_and_common_observable_units": True,
    "every_component_must_be_nonnegative": True,
    "constraint_component_required": True,
    "components_are_summed_without_cancellation": True,
    "constraint_to_Raychaudhuri_observable_stability_map_supplied": False,
    "constraint_error_cannot_yet_be_converted_to_DEF1_margin": True,
}
EXPECTED_PROOF_KEYS = {
    "PROTO4_must_validate_exactly",
    "PROTO4_freeze_must_be_canonical_and_passed",
    "HLT1_and_CAL0_must_remain_historical_inputs",
    "all_26_legacy_stops_must_be_partitioned_exactly_once",
    "GR0_must_activate_only_universal_runtime_stops",
    "candidate_stops_must_require_exact_branch_owner",
    "FGCQR_owner_artifacts_must_have_FGCQR_scope",
    "SGBL_missing_owner_must_fail_closed",
    "proper_reference_deviations_and_center_limits_must_be_executable",
    "weighted_field_derivative_and_RMS_spectra_must_be_executable",
    "last_bin_diagnostic_must_not_decide_admission",
    "nested_spectral_nonconvergence_must_fail",
    "term_sum_constraint_normalization_must_be_executable",
    "common_event_constraint_nonconvergence_must_fail",
    "Richardson_and_error_sum_must_be_additive_only",
    "declared_compact_initial_profile_must_pass_weighted_spectral_control_on_three_grids",
    "canonical_hash_bound_result_required",
    "no_new_output_namespace_may_be_read_or_written",
}
EXPECTED_CLAIMS = {
    "PROTO4_resolved_holdout_manifest_authorized": False,
    "fresh_GR0_calibration_authorized": False,
    "SGBL_comparison_execution_authorized": False,
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
IMPLEMENTATION = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto4_admission.py",
    REPOSITORY / "scripts/reproduce_fgc_hlt2_mon2.py",
)


def _pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON object key: {key}")
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
    if isinstance(value, (tuple, list)):
        return [_serial(item) for item in value]
    if isinstance(value, np.generic):
        return _serial(value.item())
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


def _path(name: str, value: Any, expected: str) -> Path:
    if value != expected:
        raise ValueError(f"{name} must name {expected}")
    candidate = (REPOSITORY / expected).resolve()
    if not candidate.is_file() or _rel(candidate) != expected:
        raise ValueError(f"{name} must be a canonical existing repository path")
    return candidate


def _load_json(path: Path, name: str, *, canonical: bool = True) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    try:
        value = json.loads(source, object_pairs_hook=_pairs)
    except (json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{name} must be valid unique-key JSON") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a JSON object")
    if canonical and source != _canonical(value):
        raise ValueError(f"{name} must be canonical sorted JSON")
    return value


def _check_result(
    path: Path,
    *,
    artifact_id: str,
    required_gate: str,
    target_branch: str | None = None,
) -> dict[str, Any]:
    value = _load_json(path, artifact_id)
    if value.get("artifact_id") != artifact_id or value.get("project_version") != PROJECT_VERSION:
        raise ValueError(f"{artifact_id} identity differs")
    if value.get("gate_status", {}).get(required_gate) is not True:
        raise ValueError(f"{artifact_id} required gate is not true")
    if target_branch is not None and value.get("scope_bindings", {}).get("target_branch") != target_branch:
        raise ValueError(f"{artifact_id} branch scope differs")
    return value


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    expected_root = {
        "schema_version", "artifact_id", "project_version", "metric_signature",
        "riemann_convention", *EXPECTED_PATHS, "scope", "branch_applicability",
        "spectral_admission", "constraint_admission", "error_budget",
        "proof_contract", "claims",
    }
    if set(raw) != expected_root:
        raise ValueError("HLT2 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("HLT2 identity or convention differs")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    if raw["scope"] != EXPECTED_SCOPE:
        raise ValueError("HLT2 scope differs")
    if raw["branch_applicability"] != EXPECTED_BRANCH:
        raise ValueError("HLT2 branch applicability differs")
    if raw["spectral_admission"] != EXPECTED_SPECTRAL:
        raise ValueError("HLT2 spectral admission differs")
    if raw["constraint_admission"] != EXPECTED_CONSTRAINT:
        raise ValueError("HLT2 constraint admission differs")
    if raw["error_budget"] != EXPECTED_ERROR:
        raise ValueError("HLT2 error-budget contract differs")
    proof = raw["proof_contract"]
    if set(proof) != EXPECTED_PROOF_KEYS or any(value is not True for value in proof.values()):
        raise ValueError("HLT2 proof contract differs")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("HLT2 claims must remain fail-closed")

    with paths["protocol_config"].open("rb") as handle:
        protocol = validate_sf1_protocol(tomllib.load(handle))
    if protocol.get("artifact_id") != SF1_PROTOCOL_V4_ARTIFACT_ID:
        raise ValueError("HLT2 protocol binding differs")
    freeze = _check_result(
        paths["protocol_freeze_result"],
        artifact_id="FGC-1-PRO4-FRZ1",
        required_gate="PROTO4_outcome_neutral_protocol_frozen",
    )
    if freeze.get("gate_status", {}).get("PROTO4_resolved_holdout_manifest_authorized") is not False:
        raise ValueError("PROTO4 freeze no longer keeps the holdout closed")
    _check_result(
        paths["historical_monitor_result"],
        artifact_id="FGC-1-HLT1-MON1",
        required_gate="classical_health_and_typed_stop_monitoring_verified",
    )
    _check_result(
        paths["diagnosis_result"],
        artifact_id="FGC-1-CAL0-PREF2",
        required_gate="PROTO3_pre_holdout_numerical_contract_obstruction_verified",
    )
    _check_result(
        paths["nonlinear_source_result"],
        artifact_id="FGC-1-SRC1-NL1",
        required_gate="nonlinear_REF1_source_and_branch_solver_verified",
        target_branch="FGC-QR",
    )
    _check_result(
        paths["run_domain_result"],
        artifact_id="FGC-1-DOM4-RUN1",
        required_gate="nonzero_classical_spherical_run_envelope_passed",
        target_branch="FGC-QR",
    )
    _check_result(
        paths["FGCQR_health_result"],
        artifact_id="FGC-1-HYP2-MD1",
        required_gate="quantitative_all_covector_weak_coupling_health_envelope_passed",
        target_branch="FGC-QR",
    )
    return {"raw": raw, "paths": paths, "protocol": protocol}


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    return _load_json(path, "HLT2 result", canonical=True)


def _constraint_sample(point_count: int, residual_scale: float, *, time: float = 0.0) -> ConstraintEventSample:
    residual = np.full((point_count, len(PROTO4_CONSTRAINT_ORDER)), residual_scale)
    terms = np.zeros((point_count, len(PROTO4_CONSTRAINT_ORDER), 3))
    return ConstraintEventSample(
        point_count,
        time,
        normalized_constraint_norms(
            residual,
            terms,
            planck_mass=2.0,
            length_unit=4.0,
        ),
    )


def _synthetic_spectral_controls() -> dict[str, Any]:
    coordinate = np.arange(256, dtype=np.float64) * (2.0 * pi / 256.0)
    low = windowed_spectral_power_budget(
        np.sin(2.0 * coordinate), coordinate, cutoff=16.0
    )
    high = windowed_spectral_power_budget(
        np.sin(112.0 * coordinate), coordinate, cutoff=16.0
    )
    compact = np.zeros_like(coordinate)
    inside = np.abs(coordinate - pi) < 1.0
    x = coordinate[inside] - pi
    compact[inside] = np.exp(1.0 - 1.0 / (1.0 - x * x))
    compact_budget = windowed_spectral_power_budget(
        compact, coordinate, cutoff=16.0
    )
    times = np.linspace(-4.0, 0.0, 64)
    temporal = windowed_spectral_power_budget(
        np.sin(2.0 * times),
        times,
        cutoff=16.0,
        window=causal_past_window(times),
    )
    if (
        not low.individual_admission_passed
        or high.individual_admission_passed
        or not compact_budget.individual_admission_passed
        or not compact_budget.last_bin_alias_diagnostic
        or not temporal.individual_admission_passed
    ):
        raise RuntimeError("HLT2 synthetic spectral controls failed")

    def fixture(field: float, derivative: float) -> SpectralPowerBudget:
        return SpectralPowerBudget(
            sample_count=129,
            top_band_first_bin=56,
            nyquist_bin=64,
            total_field_power=1.0,
            top_band_field_power_fraction=field,
            total_derivative_weighted_power=1.0,
            top_band_derivative_weighted_power_fraction=derivative,
            rms_angular_scale=1.0,
            rms_scale_over_cutoff=1.0 / 16.0,
            last_occupied_bin=64,
            last_bin_alias_diagnostic=True,
            individual_admission_passed=True,
        )

    nested_records = [
        {name: fixture(1.0e-8 * scale, 1.0e-6 * scale) for name in PROTO4_SPECTRAL_FIELD_ORDER}
        for scale in (1.0, 0.125, 0.015625)
    ]
    nested_pass = nested_spectral_admission((129, 257, 513), nested_records)
    nested_fail_records = [dict(record) for record in nested_records]
    nested_fail_records[-1]["chi_over_Lambda"] = fixture(1.0e-8, 1.0e-6)
    nested_fail = nested_spectral_admission((129, 257, 513), nested_fail_records)
    if not nested_pass.admission_passed or nested_fail.admission_passed:
        raise RuntimeError("HLT2 nested spectral controls failed")
    return {
        "low_mode": low,
        "injected_high_mode": high,
        "compact_last_bin_control": compact_budget,
        "causal_64_sample_temporal_control": temporal,
        "nested_convergent_control": nested_pass,
        "nested_nonconvergent_injection": nested_fail,
    }


def _declared_compact_profile_control() -> dict[str, Any]:
    counts = (1025, 2049, 4097)
    budgets = []
    records = []
    for point_count in counts:
        initial = construct_gr0_grid_initial_data(
            PulseParameters(chi_amplitude=2.0), point_count=point_count
        )
        profile = proper_radial_profile(
            initial.state.u,
            initial.state.q,
            initial.grid.coordinates,
            cutoff=16.0,
            measurement_radius_maximum=24.0,
        )
        proper = profile.proper_coordinates
        window = compact_vacuum_buffer_window(
            proper,
            support_minimum=float(proper[0]),
            support_maximum=float(proper[-1]),
            taper_width=float((proper[-1] - proper[0]) / 8.0),
        )
        field_budgets = spectral_field_budgets(
            profile.deviations,
            proper,
            cutoff=16.0,
            window=window,
        )
        budgets.append(field_budgets)
        records.append(
            {
                "point_count": point_count,
                "proper_measurement_radius": float(proper[-1]),
                "maximum_round_trip_interpolation_infinity": float(
                    np.max(profile.round_trip_interpolation_infinity, initial=0.0)
                ),
                "round_trip_interpolation_infinity_by_field": dict(
                    zip(PROTO4_SPECTRAL_FIELD_ORDER, profile.round_trip_interpolation_infinity)
                ),
                "field_budgets": field_budgets,
                "every_individual_budget_passed": all(
                    value.individual_admission_passed for value in field_budgets.values()
                ),
            }
        )
    nested = nested_spectral_admission(counts, budgets)
    at_least_one_last_bin = any(
        value.last_bin_alias_diagnostic for record in budgets for value in record.values()
    )
    if not nested.admission_passed or not at_least_one_last_bin:
        raise RuntimeError("declared compact profile did not pass the PROTO4 weighted control")
    return {
        "branch": "GR-0",
        "amplitude": "2",
        "event": "initial_constraint_solved_slice_t=0",
        "dynamic_outcome_used": False,
        "records": records,
        "nested_admission": nested,
        "at_least_one_last_bin_diagnostic_true": at_least_one_last_bin,
        "last_bin_diagnostic_ignored_by_decision": True,
    }


def _constraint_controls() -> dict[str, Any]:
    convergent = tuple(
        _constraint_sample(n, (1.0 / (n - 1)) ** 2) for n in (33, 65, 129)
    )
    admitted = common_event_constraint_admission(convergent)
    nonconvergent = common_event_constraint_admission(
        tuple(_constraint_sample(n, 1.0e-4) for n in (33, 65, 129))
    )
    misaligned_records = list(convergent)
    misaligned_records[-1] = _constraint_sample(
        129, (1.0 / 128.0) ** 2, time=1.0 / 8.0
    )
    misaligned = common_event_constraint_admission(misaligned_records)
    DEF1_missing_affine = common_event_constraint_admission(
        convergent,
        event_alignment_contract="DEF1_coordinate_time_affine_label_and_origin",
    )
    DEF1_affine_records = tuple(
        ConstraintEventSample(
            item.point_count,
            item.coordinate_time,
            item.norms,
            affine_label=1.5,
            affine_origin=0.0,
        )
        for item in convergent
    )
    DEF1_affine_aligned = common_event_constraint_admission(
        DEF1_affine_records,
        event_alignment_contract="DEF1_coordinate_time_affine_label_and_origin",
    )
    exact_zero = common_event_constraint_admission(
        tuple(_constraint_sample(n, 0.0) for n in (33, 65, 129))
    )
    residuals = np.asarray([[1.0, 2.0], [3.0, 4.0]])
    terms = np.asarray(
        [
            [[1.0, -2.0], [0.0, 0.0]],
            [[0.0, 0.0], [1.0, -3.0]],
        ]
    )
    normalization = normalized_constraint_norms(
        residuals,
        terms,
        planck_mass=2.0,
        length_unit=4.0,
        component_names=("control_a", "control_b"),
    )
    richardson = conservative_richardson_error(
        np.asarray((1.0, 2.0)),
        np.asarray((1.25, 1.5)),
        refinement_ratio=2.0,
        assumed_order=2.0,
        event_alignment=admitted,
    )
    misaligned_Richardson_rejected = False
    try:
        conservative_richardson_error(
            np.asarray((1.0, 2.0)),
            np.asarray((1.25, 1.5)),
            refinement_ratio=2.0,
            assumed_order=2.0,
            event_alignment=misaligned,
        )
    except ValueError:
        misaligned_Richardson_rejected = True
    error_sum = conservative_observable_error_sum(
        {"constraint": richardson, "interpolation": 0.25, "affine": 0.125}
    )
    if (
        not admitted.admission_passed
        or nonconvergent.admission_passed
        or misaligned.admission_passed
        or DEF1_missing_affine.common_event_alignment_passed
        or not DEF1_affine_aligned.common_event_alignment_passed
        or not exact_zero.admission_passed
        or tuple(normalization.component_infinity) != (12.0, 8.0)
        or abs(richardson - 1.0 / 6.0) > 1.0e-15
        or not misaligned_Richardson_rejected
        or abs(error_sum - (richardson + 0.375)) > 1.0e-15
    ):
        raise RuntimeError("HLT2 constraint/error controls failed")
    return {
        "manufactured_second_order_common_event": admitted,
        "nonconvergent_injection": nonconvergent,
        "common_event_misalignment_injection": misaligned,
        "DEF1_missing_affine_alignment_injection": DEF1_missing_affine,
        "DEF1_affine_alignment_control": DEF1_affine_aligned,
        "exact_zero_control": exact_zero,
        "term_sum_normalization_control": normalization,
        "Richardson_additive_control": {
            "estimate": richardson,
            "conservative_sum": error_sum,
            "passed_common_event_required": True,
            "misaligned_common_event_rejected": misaligned_Richardson_rejected,
            "subtraction_available": False,
        },
    }


def record(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    loaded = load_config(path)
    raw = loaded["raw"]
    paths = loaded["paths"]
    protocol = loaded["protocol"]

    gr0 = branch_stop_applicability("GR-0")
    fgc_missing = branch_stop_applicability("FGC-QR")
    fgc = branch_stop_applicability(
        "FGC-QR", candidate_definition_owner=FGCQR_BRANCH_DEFINITION_OWNER
    )
    sgbl = branch_stop_applicability("SGB-L")
    rejected_cross_owner = False
    try:
        branch_stop_applicability(
            "SGB-L", candidate_definition_owner=FGCQR_BRANCH_DEFINITION_OWNER
        )
    except ValueError:
        rejected_cross_owner = True
    if (
        not gr0.runtime_monitor_complete
        or gr0.candidate_branch_stop_ids
        or fgc_missing.runtime_monitor_complete
        or not fgc.runtime_monitor_complete
        or sgbl.runtime_monitor_complete
        or not rejected_cross_owner
    ):
        raise RuntimeError("HLT2 branch applicability controls failed")

    spectral = _synthetic_spectral_controls()
    compact = _declared_compact_profile_control()
    constraints = _constraint_controls()

    gate_status = {
        REQUIRED_GATE: True,
        "PROTO4_branch_applicability_executable": True,
        "PROTO4_weighted_spectral_admission_executable": True,
        "PROTO4_common_event_constraint_admission_executable": True,
        "declared_compact_initial_profile_weighted_spectral_control_passed": True,
        "FGCQR_candidate_health_definition_bound": True,
        "SGBL_candidate_health_definition_bound": False,
        "constraint_to_DEF1_observable_stability_map_supplied": False,
        **raw["claims"],
    }
    nonclaims = {
        "fresh_GR0_calibration_completed": False,
        "SGBL_control_evolved": False,
        "FGCQR_holdout_evolved": False,
        "collapse_or_trapped_interval_observed": False,
        "regulator_activation_observed": False,
        "positive_Raychaudhuri_margin_observed": False,
        "candidate_branch_rejected": False,
        "general_gradient_route_rejected": False,
        "retained_EFT_validity": False,
        "physical_transition": False,
        "singularity_resolution": False,
        "child_domain_or_topology": False,
        "dark_sector_mechanism": False,
        "varying_locally_measured_c": False,
    }
    predecessor = {raw[name]: _sha(paths[name]) for name in EXPECTED_PATHS}
    source_config = {_rel(path): _sha(path)}
    implementation = {_rel(item): _sha(item) for item in IMPLEMENTATION}
    return _serial({
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "outcome_neutral_PROTO4_successor_admission_monitor_certificate",
        "generated_by": "scripts/reproduce_fgc_hlt2_mon2.py",
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": source_config,
        "predecessor_sha256": predecessor,
        "implementation_sha256": implementation,
        "scope_bindings": {
            "protocol_artifact_id": protocol["artifact_id"],
            "protocol_config_sha256": _sha(paths["protocol_config"]),
            "protocol_semantic_holdout_contract_sha256": protocol[
                "semantic_holdout_contract_sha256"
            ],
            "physical_study_inherited_unchanged": True,
            "calibration_branch": "GR-0",
            "comparison_branch": "SGB-L",
            "holdout_branch": "FGC-QR",
            "FGCQR_definition_owner": FGCQR_BRANCH_DEFINITION_OWNER,
            "SGBL_definition_owner": None,
        },
        "certificate_contract": sorted(EXPECTED_PROOF_KEYS),
        "artifact_payload": {
            "decision": "PROTO4_branch_and_convergent_numerical_admission_formulas_are_executable_but_no_run_is_authorized",
            "branch_applicability": {
                "GR0": gr0,
                "FGCQR_without_owner": fgc_missing,
                "FGCQR_with_exact_owner": fgc,
                "SGBL_without_owner": sgbl,
                "SGBL_rejected_FGCQR_cross_owner": rejected_cross_owner,
                "all_26_historical_stops_partitioned_exactly_once": True,
            },
            "quantitative_evidence": {
                "synthetic_spectral_controls": spectral,
                "declared_compact_initial_profile_control": compact,
                "constraint_and_error_controls": constraints,
            },
            "epistemic_boundary": {
                "historical_HLT1_and_CAL0_results_rewritten": False,
                "fresh_calibration_output_namespace_accessed": False,
                "holdout_output_namespace_accessed": False,
                "dynamic_trajectory_used": False,
                "SGBL_health_definition_missing": True,
                "constraint_to_Raychaudhuri_stability_map_missing": True,
                "interpolation_round_trip_is_not_a_continuum_bound": True,
                "constraint_smallness_is_not_a_propagation_theorem": True,
                "last_bin_occupation_is_diagnostic_only": True,
                "mechanism_question_answered": False,
            },
        },
        "gate_status": gate_status,
        "nonclaims": nonclaims,
    })


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    payload = record(arguments.config.resolve())
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(_canonical(payload), encoding="utf-8")
    print(
        f"wrote {arguments.output} ({REQUIRED_GATE}=true; "
        "SGBL_definition=false; evolution_authorized=false)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
