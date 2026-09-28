#!/usr/bin/env python3
"""Regenerate the FGC-1-CAL0-PREF2 pre-holdout protocol diagnosis."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-cal0-pref2.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-cal0-pref2.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-cal0-pref2.md"

from recursive_horizons.fgc.evolution.health_monitor import (  # noqa: E402
    compact_vacuum_buffer_window,
    windowed_spectral_support,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
    pulse_fields,
)
from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    SF1_PROTOCOL_ARTIFACT_ID,
    validate_sf1_protocol,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-CAL0-PREF2"
PROJECT_VERSION = "0.11.0"
REQUIRED_TRUE_GATE = "PROTO3_pre_holdout_numerical_contract_obstruction_verified"
EXPECTED_CLAIMS = {
    "classical_spherical_diagnostic_authorized": False,
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
EXPECTED_PROOF = {
    "active_PROTO3_must_validate_unchanged",
    "declared_profile_must_be_evaluated_by_HLT1_estimator",
    "every_frozen_resolution_must_be_checked",
    "nonzero_amplitude_scaling_invariance_must_be_checked",
    "top_band_energy_fractions_must_be_reported_not_used_to_override_PROTO3",
    "HLT1_DOM4_HYP2_branch_scopes_must_be_compared_to_calibration_branch",
    "PROTO3_file_preserved_unchanged",
    "no_FGCQR_or_SGBL_evolution_outcome_access",
    "protocol_failure_is_not_FGCQR_mechanism_failure",
    "canonical_hash_bound_result_required",
}
IMPLEMENTATION = tuple(
    REPOSITORY / path
    for path in (
        "src/recursive_horizons/fgc/initial_data_preflight.py",
        "src/recursive_horizons/fgc/evolution/health_monitor.py",
        "src/recursive_horizons/fgc/evolution/gr0_direct_source.py",
        "src/recursive_horizons/fgc/evolution/gr0_calibration.py",
        "src/recursive_horizons/fgc/scoped_run_authorization.py",
        "scripts/reproduce_fgc_cal0_pref2.py",
    )
)


def _canonical(value: object) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        indent=2,
        ensure_ascii=True,
        allow_nan=False,
    ) + "\n"


def _unique_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    try:
        value = json.loads(source, object_pairs_hook=_unique_pairs)
    except json.JSONDecodeError as exc:
        raise ValueError("CAL0 result is not valid unique-key JSON") from exc
    if not isinstance(value, dict) or source != _canonical(value):
        raise ValueError("CAL0 result is not canonical sorted JSON")
    return value


def _keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(f"{name} keys differ")


def _fraction(name: str, value: object, *, positive: bool = False) -> Fraction:
    if not isinstance(value, str):
        raise TypeError(f"{name} must be a rational string")
    try:
        answer = Q(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} is not rational") from exc
    if positive and answer <= 0:
        raise ValueError(f"{name} must be positive")
    return answer


def _repo_path(name: str, value: object) -> Path:
    if not isinstance(value, str) or not value:
        raise TypeError(f"{name} must be a repository-relative path")
    candidate = (REPOSITORY / value).resolve()
    try:
        candidate.relative_to(REPOSITORY.resolve())
    except ValueError as exc:
        raise ValueError(f"{name} escapes the repository") from exc
    if not candidate.is_file():
        raise ValueError(f"{name} does not exist")
    return candidate


def _load_toml(path: Path) -> tuple[dict[str, Any], str]:
    source = path.read_text(encoding="utf-8")
    try:
        value = tomllib.loads(source)
    except tomllib.TOMLDecodeError as exc:
        raise ValueError(f"{path.name} is not valid TOML") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{path.name} root must be a table")
    return value, source


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    raw, source = _load_toml(path)
    _keys(
        "CAL0 root",
        raw,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "protocol_config",
            "health_monitor_config",
            "run_domain_config",
            "numerical_validation_config",
            "scope",
            "spectral_input",
            "applicability_diagnosis",
            "protocol_diagnosis",
            "proof_contract",
            "claims",
        },
    )
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("CAL0 identity or convention differs")
    scope = raw["scope"]
    if not isinstance(scope, Mapping) or dict(scope) != {
        "target_protocol": SF1_PROTOCOL_ARTIFACT_ID,
        "calibration_branch": "GR-0",
        "holdout_branch": "FGC-QR",
        "physical_equations": "unredefined_ACT1_VAR1",
        "diagnosis_role": "pre_holdout_input_and_gate_applicability_audit",
        "GR0_outcomes_may_be_used": True,
        "FGCQR_evolution_outcomes_inspected": False,
        "SGBL_evolution_outcomes_inspected": False,
    }:
        raise ValueError("CAL0 scope differs")
    spectral = raw["spectral_input"]
    if not isinstance(spectral, Mapping):
        raise TypeError("spectral_input must be a table")
    _keys(
        "spectral_input",
        spectral,
        {
            "field",
            "chi_amplitude",
            "center",
            "half_width",
            "measurement_radius",
            "resolutions",
            "relative_amplitude_floor",
            "absolute_amplitude_floor",
            "unresolved_top_fraction",
            "cutoff_Lambda",
            "proper_wavenumber_over_cutoff_max",
            "window_taper_width",
            "same_relative_spectrum_for_all_nonzero_chi_amplitudes",
        },
    )
    parsed = {
        key: _fraction(f"spectral_input.{key}", spectral[key], positive=True)
        for key in (
            "chi_amplitude",
            "center",
            "half_width",
            "measurement_radius",
            "relative_amplitude_floor",
            "absolute_amplitude_floor",
            "unresolved_top_fraction",
            "cutoff_Lambda",
            "proper_wavenumber_over_cutoff_max",
            "window_taper_width",
        )
    }
    if (
        spectral["field"] != "chi"
        or spectral["resolutions"] != [1025, 2049, 4097]
        or parsed
        != {
            "chi_amplitude": Q(2),
            "center": Q(12),
            "half_width": Q(2),
            "measurement_radius": Q(24),
            "relative_amplitude_floor": Q(1, 1099511627776),
            "absolute_amplitude_floor": Q(1, 10**30),
            "unresolved_top_fraction": Q(1, 8),
            "cutoff_Lambda": Q(16),
            "proper_wavenumber_over_cutoff_max": Q(1, 4),
            "window_taper_width": Q(2),
        }
        or spectral["same_relative_spectrum_for_all_nonzero_chi_amplitudes"]
        is not True
    ):
        raise ValueError("CAL0 spectral input differs")
    applicability = raw["applicability_diagnosis"]
    if not isinstance(applicability, Mapping) or dict(applicability) != {
        "PROTO3_requires_all_RUN1_health_and_constraint_stops_for_GR0_calibration": True,
        "HLT1_target_branch": "FGC-QR",
        "DOM4_target_branch": "FGC-QR",
        "HYP2_target_branch": "FGC-QR",
        "GR0_has_zero_FGC_coefficient_deformation": True,
        "GR0_acceleration_root_is_unique_affine_root": True,
        "coordinate_acceleration_cap_is_not_an_invariant_GR0_collapse_criterion": True,
        "branch_specific_calibration_applicability_map_present": False,
    }:
        raise ValueError("CAL0 applicability diagnosis differs")
    diagnosis = raw["protocol_diagnosis"]
    if not isinstance(diagnosis, Mapping) or dict(diagnosis) != {
        "declared_compact_pulse_fails_frozen_spectral_gate_at_t0": True,
        "NUM1_validated_estimator_mechanics_but_not_declared_pulse_admission": True,
        "GR0_calibration_and_FGCQR_target_health_are_conflated": True,
        "resolved_PRO3_holdout_manifest_cannot_be_honestly_sealed": True,
        "premise_based_successor_protocol_required": True,
        "outcome_based_FGCQR_parameter_change": False,
    }:
        raise ValueError("CAL0 protocol diagnosis differs")
    proof = raw["proof_contract"]
    if not isinstance(proof, Mapping):
        raise TypeError("proof_contract must be a table")
    _keys("proof_contract", proof, EXPECTED_PROOF)
    if any(value is not True for value in proof.values()):
        raise ValueError("every CAL0 proof-contract premise must remain true")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("CAL0 claims must remain fail-closed")

    loaded: dict[str, tuple[Path, dict[str, Any], str]] = {}
    for key in (
        "protocol_config",
        "health_monitor_config",
        "run_domain_config",
        "numerical_validation_config",
    ):
        item_path = _repo_path(key, raw[key])
        item_raw, item_source = _load_toml(item_path)
        loaded[key] = (item_path, item_raw, item_source)
    raw_protocol = loaded["protocol_config"][1]
    protocol = validate_sf1_protocol(raw_protocol)
    if protocol["artifact_id"] != SF1_PROTOCOL_ARTIFACT_ID:
        raise ValueError("CAL0 did not load active PROTO3")

    # Bind the diagnosis to the *declared* PROTO3 pulse, grids, measurement
    # region, and HLT1 estimator rather than to an independently chosen
    # look-alike fixture.  L0=4 is frozen by PROTO3's units declaration.
    protocol_chi = raw_protocol["initial_data"]["chi"]
    protocol_numerics = raw_protocol["numerics"]
    protocol_parameters = raw_protocol["dimensionless_parameters"]
    protocol_health = raw_protocol["health_stops"]
    health_monitor = loaded["health_monitor_config"][1]
    health_spectral = health_monitor["spectral_estimators"]
    if (
        protocol_chi["profile"] != "compact_C_infinity_bump_of_r_times_chi"
        or protocol_chi["unit_bump_formula"]
        != "B(x)=exp(1-1/(1-x^2))_for_abs(x)<1_and_0_otherwise"
        or protocol_chi["profile_formula"]
        != "r_times_chi(r)=amplitude*B((r-center)/half_width)"
        or Q(protocol_chi["center_over_L0"]) * 4 != parsed["center"]
        or Q(protocol_chi["half_width_over_L0"]) * 4 != parsed["half_width"]
        or str(parsed["chi_amplitude"])
        not in protocol_chi["amplitude_candidates"]
        or str(Q(5)) not in protocol_chi["amplitude_candidates"]
    ):
        raise ValueError("CAL0 pulse does not equal the declared PROTO3 chi input")
    if (
        tuple(protocol_numerics["resolutions"]) != tuple(spectral["resolutions"])
        or Q(protocol_numerics["measurement_radius_max_over_L0"]) * 4
        != parsed["measurement_radius"]
        or Q(protocol_parameters["cutoff_Lambda"]) != parsed["cutoff_Lambda"]
        or Q(protocol_health["proper_wavenumber_over_Lambda_max"])
        != parsed["proper_wavenumber_over_cutoff_max"]
    ):
        raise ValueError("CAL0 grids, measurement region, or cutoff differ from PROTO3")
    if (
        Q(health_spectral["relative_amplitude_floor"])
        != parsed["relative_amplitude_floor"]
        or Q(health_spectral["absolute_amplitude_floor"])
        != parsed["absolute_amplitude_floor"]
        or Q(health_spectral["unresolved_top_fraction"])
        != parsed["unresolved_top_fraction"]
        or health_spectral["support_frequency"]
        != "largest_rfft_bin_above_max_absolute_and_relative_floor"
        or health_spectral["alias_rule"]
        != "support_in_top_one_eighth_is_a_stop"
    ):
        raise ValueError("CAL0 estimator does not equal the frozen HLT1 contract")
    return {
        "raw": raw,
        "source_sha256": sha256(source.encode("utf-8")).hexdigest(),
        "spectral": parsed,
        "resolutions": tuple(spectral["resolutions"]),
        "loaded": loaded,
        "protocol": protocol,
        "raw_protocol": raw_protocol,
    }


def _spectral_scan(configuration: Mapping[str, Any]) -> dict[str, Any]:
    values = configuration["spectral"]
    records: list[dict[str, Any]] = []
    support_controls: dict[int, tuple[int, int, bool]] = {}
    for amplitude in (values["chi_amplitude"], Q(5)):
        for point_count in configuration["resolutions"]:
            full_coordinates = np.linspace(0.0, 128.0, point_count)
            coordinates = full_coordinates[
                full_coordinates <= float(values["measurement_radius"])
            ]
            parameters = PulseParameters(
                chi_amplitude=float(amplitude),
                center=float(values["center"]),
                half_width=float(values["half_width"]),
            )
            samples = np.fromiter(
                (
                    0.0
                    if radius == 0.0
                    else pulse_fields(float(radius), parameters)["chi"]
                    for radius in coordinates
                ),
                dtype=np.float64,
                count=coordinates.size,
            )
            window = compact_vacuum_buffer_window(
                coordinates,
                support_minimum=0.0,
                support_maximum=float(values["measurement_radius"]),
                taper_width=float(values["window_taper_width"]),
            )
            support = windowed_spectral_support(
                samples,
                coordinates,
                window=window,
                relative_amplitude_floor=float(values["relative_amplitude_floor"]),
                absolute_amplitude_floor=float(values["absolute_amplitude_floor"]),
                unresolved_top_fraction=float(values["unresolved_top_fraction"]),
            )
            transformed = np.fft.rfft(samples * window)
            energy = np.abs(transformed) ** 2
            bins = np.arange(energy.size, dtype=np.float64)
            top_energy = float(
                np.sum(energy[support.first_unresolved_bin :])
                / np.sum(energy)
            )
            gradient_energy = bins**2 * energy
            top_gradient_energy = float(
                np.sum(gradient_energy[support.first_unresolved_bin :])
                / np.sum(gradient_energy)
            )
            key = (support.support_bin, support.first_unresolved_bin, support.alias_band_occupied)
            if amplitude == values["chi_amplitude"]:
                support_controls[point_count] = key
            elif key != support_controls[point_count]:
                raise ValueError("nonzero amplitude scaling changed the relative support decision")
            records.append(
                {
                    "amplitude": str(amplitude),
                    "point_count": point_count,
                    "measurement_sample_count": int(coordinates.size),
                    "spacing": float(full_coordinates[1] - full_coordinates[0]),
                    "support_bin": support.support_bin,
                    "nyquist_bin": support.nyquist_bin,
                    "first_unresolved_bin": support.first_unresolved_bin,
                    "angular_support": support.angular_support,
                    "proper_wavenumber_over_cutoff": (
                        support.angular_support / float(values["cutoff_Lambda"])
                    ),
                    "alias_band_occupied": support.alias_band_occupied,
                    "top_eighth_field_energy_fraction": top_energy,
                    "top_eighth_gradient_energy_fraction": top_gradient_energy,
                }
            )
    if not records or any(not item["alias_band_occupied"] for item in records):
        raise ValueError("declared PROTO3 pulse unexpectedly passed its frozen alias rule")
    return {
        "estimator": "HLT1_windowed_spectral_support",
        "records": records,
        "all_frozen_resolutions_fail_at_initial_time": True,
        "nonzero_amplitude_scaling_invariance_verified": True,
        "tail_energy_report_does_not_override_PROTO3_boolean": True,
    }


def _applicability(configuration: Mapping[str, Any]) -> dict[str, Any]:
    loaded = configuration["loaded"]
    protocol = loaded["protocol_config"][1]
    hlt = loaded["health_monitor_config"][1]
    dom4 = loaded["run_domain_config"][1]
    num1 = loaded["numerical_validation_config"][1]
    hyp_path = _repo_path(
        "HLT1 multidirectional_health_config",
        hlt["multidirectional_health_config"],
    )
    hyp2, hyp_source = _load_toml(hyp_path)
    if (
        protocol["calibration"]["branch"] != "GR-0"
        or protocol["calibration"]["all_RUN1_health_and_constraint_stops_apply"]
        is not True
        or hlt["scope"]["target_branch"] != "FGC-QR"
        or dom4["scope"]["target_branch"] != "FGC-QR"
        or hyp2["scope"]["target_branch"] != "FGC-QR"
    ):
        raise ValueError("CAL0 branch-scope premises differ")
    if "branch_applicability" in protocol or "branch_applicability" in hlt:
        raise ValueError("PROTO3 unexpectedly contains a branch-applicability map")
    if num1["scope"]["FGCQR_trajectory_evaluated"] is not False:
        raise ValueError("NUM1 scope unexpectedly includes an FGC-QR trajectory")
    return {
        "calibration_branch": "GR-0",
        "PROTO3_all_RUN1_stops_apply": True,
        "target_branch_owners": {
            "HLT1": hlt["scope"]["target_branch"],
            "DOM4": dom4["scope"]["target_branch"],
            "HYP2": hyp2["scope"]["target_branch"],
        },
        "DOM4_coordinate_acceleration_limit": dom4["source_branch"][
            "acceleration_infinity_maximum"
        ],
        "branch_applicability_map_present": False,
        "NUM1_declared_profile_admission_tested": False,
        "scope_conflation_verified": True,
        "hyp2_config_path": str(hyp_path.relative_to(REPOSITORY)),
        "hyp2_config_sha256": sha256(hyp_source.encode("utf-8")).hexdigest(),
    }


def record(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    configuration = load_config(path)
    loaded = configuration["loaded"]
    spectral = _spectral_scan(configuration)
    applicability = _applicability(configuration)
    predecessors = {
        str(item_path.relative_to(REPOSITORY)): sha256(source.encode("utf-8")).hexdigest()
        for item_path, _raw, source in loaded.values()
    }
    implementation = {
        str(item.relative_to(REPOSITORY)): sha256(item.read_bytes()).hexdigest()
        for item in IMPLEMENTATION
    }
    source_relative = str(path.resolve().relative_to(REPOSITORY.resolve()))
    payload = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "pre_holdout_protocol_numerical_contract_obstruction",
        "generated_by": "scripts/reproduce_fgc_cal0_pref2.py",
        "derivation_document": str(OWNER_DOCUMENT.relative_to(REPOSITORY)),
        "derivation_document_sha256": sha256(OWNER_DOCUMENT.read_bytes()).hexdigest(),
        "source_config_sha256": {
            source_relative: configuration["source_sha256"]
        },
        "predecessor_sha256": predecessors,
        "implementation_sha256": implementation,
        "certificate_contract": {
            "artifact_specific_payload_validated": True,
            "canonical_reproduction_passed": True,
            "scope_bindings_verified": True,
            "outcome_data_not_used_to_select_certificate_contract": True,
        },
        "artifact_payload": {
            "spectral_input_obstruction": spectral,
            "branch_applicability_obstruction": applicability,
            "successor_requirements": {
                "branch_specific_stop_applicability_required": True,
                "spectral_tail_weight_and_convergence_required": True,
                "deterministic_constraint_scale_and_common_event_convergence_required": True,
                "action_profiles_amplitude_order_and_holdout_cases_must_remain_unchanged": True,
                "new_output_namespace_required": True,
            },
            "epistemic_boundary": {
                "protocol_design_route_closed": "PROTO3_as_executable_resolved_holdout_contract",
                "FGCQR_evolution_outcome_inspected": False,
                "SGBL_evolution_outcome_inspected": False,
                "FGCQR_action_or_mechanism_tested": False,
                "FGCQR_action_or_mechanism_rejected": False,
                "general_gradient_mechanism_rejected": False,
            },
        },
        "gate_status": {
            REQUIRED_TRUE_GATE: True,
            "PROTO3_resolved_holdout_manifest_authorized": False,
            "PROTO4_premise_revision_required": True,
            "classical_spherical_diagnostic_authorized": False,
            "FGCQR_holdout_execution_authorized": False,
            "retained_EFT_evolution_authorized": False,
            "physical_transition_claim_authorized": False,
        },
        "nonclaims": EXPECTED_CLAIMS,
    }
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    value = record(arguments.config)
    arguments.output.write_text(_canonical(value), encoding="utf-8")


if __name__ == "__main__":
    main()
