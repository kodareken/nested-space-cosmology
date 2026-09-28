#!/usr/bin/env python3
"""Regenerate the FGC-1-ID0-PREF1 protocol-preflight obstruction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
import subprocess
import sys
import tomllib
from typing import Any, Mapping


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-id0-pref1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-id0-pref1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-id0-pref1.md"

from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
    exact_gr0_constraint_crosscheck,
    observed_convergence_order,
    proto1_compactness_upper_bound,
    solve_gr0_initial_slice,
)
from recursive_horizons.fgc.scoped_run_authorization import (  # noqa: E402
    RUN1_SYM1_ARTIFACT_ID,
    SF1_PROTOCOL_V1_ARTIFACT_ID,
    validate_sf1_protocol,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-ID0-PREF1"
PROJECT_VERSION = "0.11.0"
REQUIRED_TRUE_GATE = "PROTO1_initial_data_preflight_obstruction_verified"
EXPECTED_ROOT_KEYS = {
    "schema_version",
    "artifact_id",
    "project_version",
    "metric_signature",
    "riemann_convention",
    "protocol_config",
    "historical_preflight_checkpoint",
    "scope",
    "protocol_diagnosis",
    "exact_controls",
    "analytic_bound",
    "numerical_preflight",
    "proof_contract",
    "claims",
}
HISTORICAL_CHECKPOINT = {
    "commit": "d4f0cc8f58408ef4ee12fb32fe3619916e231795",
    "classification": "immutable_pre_PROTO2_ID0_checkpoint",
    "id0_config_path": "configs/fgc/fgc-1-id0-pref1.toml",
    "id0_config_sha256": "305338f3b940f046b59072448dfb010b7fa2701e290573c08703a2cac59a7086",
    "id0_result_path": "results/fgc-1-id0-pref1.json",
    "id0_result_sha256": "ac9f0726a1e336dc41dd100612d2ca3c64a1913449ac2d3844f42ee31d86ba22",
    "run1_config_path": "configs/fgc/fgc-1-run1-sym1.toml",
    "run1_config_sha256": "7b4341bf5b9712648daff4d247de0a786e1b19fc0ea1aeac761a1f4da39424f0",
    "run1_result_path": "results/fgc-1-run1-sym1.json",
    "run1_result_sha256": "269d0d4b594c87b6d54a84ef860056322d4c3174653af32ab06f9d02a84deca8",
    "con4_config_path": "configs/fgc/fgc-1-con4-phy1.toml",
    "con4_config_sha256": "f438db34cf4895fa990f4574a6e3e32dfdf59a3bc864d3505e40a37de2c22b03",
    "con4_result_path": "results/fgc-1-con4-phy1.json",
    "con4_result_sha256": "8fb15d3d7298b27a1607084bfc5acd8a88d67a427e9e4583ef6560000240d1e7",
    "ctr1_config_path": "configs/fgc/fgc-1-ctr1-reg1.toml",
    "ctr1_config_sha256": "c4eb5a6d646ea6bb1d8bd43d3c3911ce4c2a5182c8dde6e25a5221e07f7415a8",
    "ctr1_result_path": "results/fgc-1-ctr1-reg1.json",
    "ctr1_result_sha256": "998b3577556a79290da07eaeec90224406239a55560f4f9459b283d57f15979d",
}
EXPECTED_PROOF_CONTRACT = {
    "specialized_constraints_derived_from_declared_slice",
    "specialized_constraints_cross_checked_against_complete_unredefined_evaluator",
    "two_nontrivial_exact_rational_controls_required",
    "entire_original_amplitude_box_bounded_analytically",
    "inverse_metric_bootstrap_closed",
    "independent_RK4_and_SSPRK3_preflight_required",
    "all_original_candidates_checked_at_nested_resolutions",
    "successor_candidate_scan_uses_GR0_only",
    "PROTO1_file_preserved_unchanged",
    "no_FGCQR_holdout_execution",
    "protocol_failure_is_not_FGCQR_mechanism_failure",
    "canonical_hash_bound_result_required",
}
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
EXACT_CONTROL_KEYS = {
    "radius",
    "radial_metric",
    "angular_extrinsic_curvature",
    "radial_metric_derivative",
    "angular_extrinsic_curvature_derivative",
    "phi",
    "phi_r",
    "phi_pi",
    "chi",
    "chi_r",
    "chi_pi",
}
IMPLEMENTATION = tuple(
    REPOSITORY / path
    for path in (
        "src/recursive_horizons/fgc/initial_data_preflight.py",
        "src/recursive_horizons/fgc/constraint_system.py",
        "src/recursive_horizons/fgc/modified_harmonic_constraints.py",
        "src/recursive_horizons/fgc/spherical_reduction.py",
        "src/recursive_horizons/fgc/scoped_run_authorization.py",
        "scripts/reproduce_fgc_id0_pref1.py",
    )
)


def _keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        missing = sorted(expected - set(value))
        extra = sorted(set(value) - expected)
        raise ValueError(f"{name} keys differ; missing={missing!r}, extra={extra!r}")


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _path(name: str, value: Any, expected: str) -> Path:
    if not isinstance(value, str) or value != expected:
        raise ValueError(f"{name} must name {expected}")
    path = (REPOSITORY / value).resolve()
    try:
        relative = path.relative_to(REPOSITORY).as_posix()
    except ValueError as exc:
        raise ValueError(f"{name} must stay inside repository") from exc
    if relative != value or not path.is_file():
        raise ValueError(f"{name} must be canonical, traversal free, and existing")
    return path


def _text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _fraction(name: str, value: Any, *, positive: bool = False) -> Fraction:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a canonical rational string")
    try:
        answer = Q(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{name} must be a canonical rational string") from exc
    if _text(answer) != value:
        raise ValueError(f"{name} must be a canonical rational string")
    if positive and answer <= 0:
        raise ValueError(f"{name} must be positive")
    return answer


def _rational_list(name: str, value: Any, *, positive: bool = False) -> tuple[Fraction, ...]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"{name} must be a nonempty rational-string list")
    result = tuple(_fraction(f"{name}[{index}]", item, positive=positive) for index, item in enumerate(value))
    if len(set(result)) != len(result):
        raise ValueError(f"{name} must contain unique values")
    return result


def _pairs(pairs):
    output = {}
    for key, value in pairs:
        if key in output:
            raise ValueError(f"duplicate JSON object key: {key}")
        output[key] = value
    return output


def _load_json(path: Path, name: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_pairs)
    except (json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"{name} must be valid unique-key JSON") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a JSON object")
    return value


def _serial(value: Any) -> Any:
    if isinstance(value, Fraction):
        return _text(value)
    if is_dataclass(value):
        return _serial(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_serial(item) for item in value]
    if value is None or isinstance(value, (str, bool, int, float)):
        if isinstance(value, float) and not isfinite(value):
            raise ValueError("nonfinite floating value cannot enter canonical evidence")
        return value
    raise TypeError(f"unsupported evidence value {type(value).__name__}")


def _canonical(value: Any) -> str:
    return json.dumps(_serial(value), indent=2, sort_keys=True, allow_nan=False) + "\n"


def _git(*arguments: str) -> bytes:
    completed = subprocess.run(
        ["git", *arguments],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    )
    if completed.returncode != 0:
        detail = completed.stderr.decode("utf-8", errors="replace").strip()
        raise ValueError(f"historical Git checkpoint is unavailable: {detail}")
    return completed.stdout


def _git_blob(commit: str, path: str) -> bytes:
    if not path or path.startswith("/") or ".." in Path(path).parts:
        raise ValueError("historical checkpoint path must be canonical and relative")
    return _git("show", f"{commit}:{path}")


def _historical_json(commit: str, path: str, name: str) -> dict[str, Any]:
    source = _git_blob(commit, path)
    try:
        value = json.loads(source.decode("utf-8"), object_pairs_hook=_pairs)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise ValueError(f"historical {name} must be valid unique-key JSON") from exc
    if not isinstance(value, dict) or source.decode("utf-8") != _canonical(value):
        raise ValueError(f"historical {name} must use canonical sorted JSON")
    return value


def _validate_historical_result_ledgers(
    result: Mapping[str, Any], *, commit: str, name: str
) -> None:
    for ledger_name in (
        "source_config_sha256",
        "predecessor_sha256",
        "implementation_sha256",
    ):
        ledger = result.get(ledger_name)
        if not isinstance(ledger, Mapping) or not ledger:
            raise ValueError(f"historical {name} {ledger_name} must be nonempty")
        for path, digest in ledger.items():
            if not isinstance(path, str) or not isinstance(digest, str):
                raise ValueError(f"historical {name} {ledger_name} entry differs")
            if sha256(_git_blob(commit, path)).hexdigest() != digest:
                raise ValueError(f"historical {name} ledger is stale for {path}")
    document = result.get("derivation_document")
    document_digest = result.get("derivation_document_sha256")
    if not isinstance(document, str) or not isinstance(document_digest, str):
        raise ValueError(f"historical {name} derivation binding is absent")
    if sha256(_git_blob(commit, document)).hexdigest() != document_digest:
        raise ValueError(f"historical {name} derivation binding is stale")


def _validate_run1(record: Mapping[str, Any]) -> None:
    if record.get("artifact_id") != RUN1_SYM1_ARTIFACT_ID or record.get("project_version") != PROJECT_VERSION:
        raise ValueError("RUN1 predecessor identity differs")
    status = record.get("gate_status")
    if not isinstance(status, Mapping):
        raise ValueError("RUN1 gate status is absent")
    for gate in (
        "classical_spherical_diagnostic_authorized",
        "FGCQR_holdout_execution_authorized",
        "retained_EFT_evolution_authorized",
        "physical_transition_claim_authorized",
    ):
        if status.get(gate) is not False:
            raise ValueError("RUN1 must remain fail-closed before the preflight")
    audit = record.get("scoped_run_authorization_audit")
    if not isinstance(audit, Mapping):
        raise ValueError("RUN1 scoped audit is absent")
    audits = audit.get("authorization_audits")
    if not isinstance(audits, Mapping):
        raise ValueError("RUN1 authorization audits are absent")
    classical = audits.get("classical_spherical_diagnostic")
    if not isinstance(classical, Mapping) or classical.get("passed_predicate_count") != 2:
        raise ValueError("ID0-PREF1 requires the sealed pre-ID1 RUN1 2/8 stop")


def _validate_historical_checkpoint(value: Any) -> dict[str, Any]:
    checkpoint = value if isinstance(value, Mapping) else None
    if checkpoint is None or dict(checkpoint) != HISTORICAL_CHECKPOINT:
        raise ValueError("ID0-PREF1 historical pre-PROTO2 checkpoint differs")
    commit = HISTORICAL_CHECKPOINT["commit"]
    resolved = _git("rev-parse", "--verify", f"{commit}^{{commit}}").decode("ascii").strip()
    if resolved != commit:
        raise ValueError("ID0-PREF1 historical checkpoint does not resolve exactly")
    ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    )
    if ancestor.returncode != 0:
        raise ValueError("ID0-PREF1 historical checkpoint is not an ancestor of HEAD")

    file_hashes: dict[str, str] = {}
    for key, path in checkpoint.items():
        if not key.endswith("_path"):
            continue
        digest_key = key.removesuffix("_path") + "_sha256"
        digest = checkpoint.get(digest_key)
        if not isinstance(path, str) or not isinstance(digest, str):
            raise ValueError("ID0-PREF1 historical file binding differs")
        observed = sha256(_git_blob(commit, path)).hexdigest()
        if observed != digest:
            raise ValueError(f"ID0-PREF1 historical blob hash differs for {path}")
        file_hashes[path] = observed

    historical_id0 = _historical_json(
        commit, checkpoint["id0_result_path"], "ID0 result"
    )
    _validate_historical_result_ledgers(
        historical_id0, commit=commit, name="ID0 result"
    )
    id0_status = historical_id0.get("gate_status")
    id0_scope = historical_id0.get("scope_bindings")
    if (
        historical_id0.get("artifact_id") != ARTIFACT_ID
        or not isinstance(id0_status, Mapping)
        or id0_status.get(REQUIRED_TRUE_GATE) is not True
        or not isinstance(id0_scope, Mapping)
        or id0_scope.get("target_protocol_artifact_id")
        != SF1_PROTOCOL_V1_ARTIFACT_ID
        or id0_scope.get("FGCQR_holdout_outcomes_inspected") is not False
    ):
        raise ValueError("historical ID0 result does not preserve the pre-holdout obstruction")

    historical_run1 = _historical_json(
        commit, checkpoint["run1_result_path"], "RUN1 result"
    )
    _validate_historical_result_ledgers(
        historical_run1, commit=commit, name="RUN1 result"
    )
    _validate_run1(historical_run1)

    historical_con4 = _historical_json(
        commit, checkpoint["con4_result_path"], "CON4 result"
    )
    historical_ctr1 = _historical_json(
        commit, checkpoint["ctr1_result_path"], "CTR1 result"
    )
    for name, result, artifact_id, gate in (
        (
            "CON4 result",
            historical_con4,
            "FGC-1-CON4-PHY1",
            "physical_gauge_reduction_constraint_system_closed",
        ),
        (
            "CTR1 result",
            historical_ctr1,
            "FGC-1-CTR1-REG1",
            "regular_center_formulation_verified",
        ),
    ):
        _validate_historical_result_ledgers(result, commit=commit, name=name)
        status = result.get("gate_status")
        if (
            result.get("artifact_id") != artifact_id
            or not isinstance(status, Mapping)
            or status.get(gate) is not True
        ):
            raise ValueError(f"historical {name} positive gate differs")

    return {
        "commit": commit,
        "classification": checkpoint["classification"],
        "commit_is_ancestor_of_HEAD": True,
        "file_sha256": file_hashes,
        "original_ID0_canonical_ledgers_verified_at_commit": True,
        "original_RUN1_two_of_eight_stop_verified_at_commit": True,
        "original_CON4_and_CTR1_positive_gates_verified_at_commit": True,
        "FGCQR_holdout_outcomes_inspected_before_ID0": False,
    }


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    source = path.resolve().read_bytes()
    raw = tomllib.loads(source.decode("utf-8"))
    _keys("root", raw, EXPECTED_ROOT_KEYS)
    if {
        "schema_version": raw["schema_version"],
        "artifact_id": raw["artifact_id"],
        "project_version": raw["project_version"],
        "metric_signature": raw["metric_signature"],
        "riemann_convention": raw["riemann_convention"],
    } != {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
    }:
        raise ValueError("ID0-PREF1 identity or convention differs")

    protocol_path = _path(
        "protocol_config",
        raw["protocol_config"],
        "configs/fgc/fgc-2-sf1-protocol.toml",
    )
    protocol_raw = tomllib.loads(protocol_path.read_text(encoding="utf-8"))
    protocol = validate_sf1_protocol(protocol_raw)
    if (
        protocol["artifact_id"] != SF1_PROTOCOL_V1_ARTIFACT_ID
        or protocol["protocol_version"] != 1
    ):
        raise ValueError("ID0-PREF1 requires the original frozen PROTO1")
    historical_checkpoint = _validate_historical_checkpoint(
        raw["historical_preflight_checkpoint"]
    )

    scope = raw["scope"]
    if not isinstance(scope, Mapping):
        raise ValueError("scope must be a table")
    _keys(
        "scope",
        scope,
        {
            "target_protocol",
            "target_branch",
            "physical_equations",
            "initial_slice",
            "areal_radius",
            "extrinsic_curvature_eigenvalues",
            "regular_center_vacuum_interior",
            "asymptotically_finite_mass_vacuum_exterior",
            "FGCQR_holdout_outcomes_inspected",
            "SGBL_outcomes_inspected",
        },
    )
    if dict(scope) != {
        "target_protocol": SF1_PROTOCOL_V1_ARTIFACT_ID,
        "target_branch": "GR-0",
        "physical_equations": "unredefined_ACT1_VAR1",
        "initial_slice": "unit_lapse_zero_shift_maximal_polar_areal",
        "areal_radius": "R=r",
        "extrinsic_curvature_eigenvalues": [
            "K^r_r=-2*k",
            "K^theta_theta=k",
            "K^phi_phi=k",
        ],
        "regular_center_vacuum_interior": True,
        "asymptotically_finite_mass_vacuum_exterior": True,
        "FGCQR_holdout_outcomes_inspected": False,
        "SGBL_outcomes_inspected": False,
    }:
        raise ValueError("ID0-PREF1 scope differs")

    diagnosis = raw["protocol_diagnosis"]
    if not isinstance(diagnosis, Mapping):
        raise ValueError("protocol_diagnosis must be a table")
    expected_diagnosis = {
        "declared_chi_profile": "r_times_chi=amplitude*B((r-center)/half_width)",
        "declared_chi_momentum": "partial_t(r_chi)=partial_r(r_chi)",
        "declared_phi_profile": "phi=amplitude*B((r-center)/half_width)",
        "phi_unit_normal_momentum_frozen_by_PROTO1": False,
        "minimal_preflight_completion": "Pi_phi=0",
        "missing_phi_momentum_alone_prevents_executable_initial_data_contract": True,
        "original_chi_amplitudes_are_too_small_under_minimal_completion": True,
        "premise_based_protocol_revision_required": True,
        "outcome_based_parameter_change": False,
    }
    if dict(diagnosis) != expected_diagnosis:
        raise ValueError("ID0-PREF1 protocol diagnosis differs")

    exact_controls = raw["exact_controls"]
    if not isinstance(exact_controls, Mapping) or set(exact_controls) != {"fixture_A", "fixture_B"}:
        raise ValueError("ID0-PREF1 exact controls differ")
    parsed_controls: dict[str, dict[str, Fraction]] = {}
    for fixture_id, fixture in exact_controls.items():
        if not isinstance(fixture, Mapping):
            raise ValueError(f"{fixture_id} must be a table")
        _keys(f"exact_controls.{fixture_id}", fixture, EXACT_CONTROL_KEYS)
        parsed_controls[fixture_id] = {
            key: _fraction(f"exact_controls.{fixture_id}.{key}", value)
            for key, value in fixture.items()
        }
        if parsed_controls[fixture_id]["radius"] <= 0 or parsed_controls[fixture_id]["radial_metric"] <= 0:
            raise ValueError(f"{fixture_id} lies outside the spacelike positive-radius domain")

    analytic = raw["analytic_bound"]
    if not isinstance(analytic, Mapping):
        raise ValueError("analytic_bound must be a table")
    _keys(
        "analytic_bound",
        analytic,
        {
            "amplitude_maximum",
            "phi_amplitude",
            "support_minimum",
            "support_maximum",
            "planck_mass",
            "assumed_inverse_radial_metric_squared_upper_bound",
            "bump_absolute_bound",
            "bump_x_derivative_absolute_bound",
            "bump_r_derivative_absolute_bound",
            "required_compactness_upper_separator",
            "protocol_initial_compactness_minimum",
            "bootstrap_must_close_strictly",
            "regular_center_forbids_hidden_vacuum_mass_or_r_minus_3_momentum_seed",
        },
    )
    parsed_analytic = {
        key: _fraction(f"analytic_bound.{key}", value, positive=True)
        for key, value in analytic.items()
        if isinstance(value, str)
    }
    if parsed_analytic != {
        "amplitude_maximum": Q(1, 8),
        "phi_amplitude": Q(1, 131072),
        "support_minimum": Q(10),
        "support_maximum": Q(14),
        "planck_mass": Q(2),
        "assumed_inverse_radial_metric_squared_upper_bound": Q(4),
        "bump_absolute_bound": Q(1),
        "bump_x_derivative_absolute_bound": Q(3),
        "bump_r_derivative_absolute_bound": Q(3, 2),
        "required_compactness_upper_separator": Q(1, 64),
        "protocol_initial_compactness_minimum": Q(1, 10),
    }:
        raise ValueError("ID0-PREF1 analytic-bound constants differ")
    if analytic["bootstrap_must_close_strictly"] is not True or analytic[
        "regular_center_forbids_hidden_vacuum_mass_or_r_minus_3_momentum_seed"
    ] is not True:
        raise ValueError("ID0-PREF1 analytic-bound premises must remain true")

    numerical = raw["numerical_preflight"]
    if not isinstance(numerical, Mapping):
        raise ValueError("numerical_preflight must be a table")
    _keys(
        "numerical_preflight",
        numerical,
        {
            "original_amplitude_candidates",
            "candidate_repair_scan",
            "methods",
            "nested_step_counts",
            "report_step_count",
            "minimum_observed_outer_mass_convergence_order",
            "maximum_cross_method_peak_compactness_difference",
            "original_box_peak_compactness_numerical_upper_bound",
            "repair_scan_peak_compactness_minimum",
            "repair_scan_peak_compactness_maximum",
            "repair_scan_is_GR0_premise_design_only",
            "repair_scan_does_not_select_a_FGCQR_outcome",
        },
    )
    original = _rational_list(
        "numerical_preflight.original_amplitude_candidates",
        numerical["original_amplitude_candidates"],
        positive=True,
    )
    if tuple(protocol["amplitude_candidates"]) != tuple(_text(value) for value in original):
        raise ValueError("ID0-PREF1 original amplitudes differ from frozen PROTO1")
    repair = _rational_list(
        "numerical_preflight.candidate_repair_scan",
        numerical["candidate_repair_scan"],
        positive=True,
    )
    if list(numerical["methods"]) != ["RK4", "SSPRK3"]:
        raise ValueError("ID0-PREF1 requires independent RK4 and SSPRK3 preflights")
    steps = numerical["nested_step_counts"]
    if steps != [32, 64, 128] or any(isinstance(value, bool) for value in steps):
        raise ValueError("ID0-PREF1 nested step counts differ")
    report_steps = numerical["report_step_count"]
    if isinstance(report_steps, bool) or not isinstance(report_steps, int) or report_steps != 4096:
        raise ValueError("ID0-PREF1 report step count differs")
    numerical_thresholds = {
        key: _fraction(f"numerical_preflight.{key}", numerical[key], positive=True)
        for key in (
            "minimum_observed_outer_mass_convergence_order",
            "maximum_cross_method_peak_compactness_difference",
            "original_box_peak_compactness_numerical_upper_bound",
            "repair_scan_peak_compactness_minimum",
            "repair_scan_peak_compactness_maximum",
        )
    }
    if numerical_thresholds != {
        "minimum_observed_outer_mass_convergence_order": Q(5, 2),
        "maximum_cross_method_peak_compactness_difference": Q(1, 100000000),
        "original_box_peak_compactness_numerical_upper_bound": Q(1, 2000),
        "repair_scan_peak_compactness_minimum": Q(1, 10),
        "repair_scan_peak_compactness_maximum": Q(3, 4),
    }:
        raise ValueError("ID0-PREF1 numerical thresholds differ")
    if numerical["repair_scan_is_GR0_premise_design_only"] is not True or numerical[
        "repair_scan_does_not_select_a_FGCQR_outcome"
    ] is not True:
        raise ValueError("ID0-PREF1 repair scan scope differs")

    proof = raw["proof_contract"]
    if not isinstance(proof, Mapping):
        raise ValueError("proof_contract must be a table")
    _keys("proof_contract", proof, EXPECTED_PROOF_CONTRACT)
    if any(value is not True for value in proof.values()):
        raise ValueError("every ID0-PREF1 proof-contract premise must remain true")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("ID0-PREF1 claims must remain fail-closed")

    return {
        "raw": raw,
        "source_sha256": sha256(source).hexdigest(),
        "protocol_path": protocol_path,
        "protocol_raw": protocol_raw,
        "protocol": protocol,
        "historical_checkpoint": historical_checkpoint,
        "exact_controls": parsed_controls,
        "analytic": parsed_analytic,
        "original_amplitudes": original,
        "repair_amplitudes": repair,
        "methods": tuple(numerical["methods"]),
        "nested_steps": tuple(steps),
        "report_steps": report_steps,
        "numerical_thresholds": numerical_thresholds,
    }


def _exact_controls(configuration: Mapping[str, Any]) -> dict[str, Any]:
    records: dict[str, Any] = {}
    for fixture_id, fixture in configuration["exact_controls"].items():
        control = exact_gr0_constraint_crosscheck(**fixture)
        if not control["exact_pair_equality"] or not control["source_is_unredefined_ACT1_VAR1"]:
            raise ValueError(f"{fixture_id} exact unredefined cross-check failed")
        records[fixture_id] = control
    return {
        "constraint_formulae": {
            "Hamiltonian": "M_Pl^2*((1-lambda^-2)/r^2+2*lambda_r/(r*lambda^3)-3*k^2)-rho=0",
            "radial_momentum": "2*M_Pl^2*(k_r+3*k/r)-P=0",
            "rho": "((Pi_phi)^2+(Pi_chi)^2+lambda^-2*((phi_r)^2+(chi_r)^2))/2+V(phi)",
            "P": "Pi_phi*phi_r+Pi_chi*chi_r",
        },
        "constraint_ODE": {
            "lambda_r": "r*lambda^3*(rho/M_Pl^2-(1-lambda^-2)/r^2+3*k^2)/2",
            "k_r": "P/(2*M_Pl^2)-3*k/r",
        },
        "fixtures": records,
        "fixture_count": len(records),
        "all_exact_pair_equalities_passed": True,
    }


def _analytic_bound(configuration: Mapping[str, Any]) -> dict[str, Any]:
    result = proto1_compactness_upper_bound()
    expected = configuration["analytic"]
    comparisons = {
        "amplitude_maximum_matches": expected["amplitude_maximum"] == Q(1, 8),
        "phi_amplitude_matches": expected["phi_amplitude"] == Q(1, 131072),
        "support_matches": (
            expected["support_minimum"],
            expected["support_maximum"],
        ) == (Q(10), Q(14)),
        "bump_bound_matches": result["bump_absolute_bound"] == expected["bump_absolute_bound"],
        "bump_derivative_bounds_match": (
            result["bump_x_derivative_absolute_bound"],
            result["bump_r_derivative_absolute_bound"],
        ) == (
            expected["bump_x_derivative_absolute_bound"],
            expected["bump_r_derivative_absolute_bound"],
        ),
        "inverse_metric_bootstrap_assumption_matches": (
            result["assumed_inverse_radial_metric_squared_upper_bound"]
            == expected["assumed_inverse_radial_metric_squared_upper_bound"]
        ),
        "compactness_strictly_below_separator": (
            result["compactness_absolute_bound"]
            < expected["required_compactness_upper_separator"]
        ),
        "separator_strictly_below_protocol_floor": (
            expected["required_compactness_upper_separator"]
            < expected["protocol_initial_compactness_minimum"]
        ),
        "inverse_metric_bootstrap_closed": result["inverse_metric_bootstrap_closed"],
    }
    if not all(comparisons.values()):
        raise ValueError("ID0-PREF1 analytic compactness bound failed")
    return {
        "derivation": {
            "bump_derivative_bound": "for_y=1/(1-x^2)>=1,_abs(B_x)<=2*y^2*exp(1-y)<=8/e<3",
            "momentum_integral": "abs(r^3*k)<=integral(abs(P)*r^3/(2*M_Pl^2),dr)",
            "mass_integral": "abs(m)<=integral((r^2*rho+abs(r^3*k*P))/(2*M_Pl^2),dr)",
            "compactness_identity": "2*m/r=1+r^2*k^2-lambda^-2",
            "bootstrap_logic": "a_first_exit_from_lambda^-2<=4_is_excluded_by_the_stronger_derived_bound",
        },
        "bound": result,
        "comparisons": comparisons,
        "applies_to_every_original_amplitude_under_Pi_phi_zero_completion": True,
    }


def _solution_summary(solution) -> dict[str, Any]:
    return {
        "method": solution.method,
        "step_count": solution.step_count,
        "peak_compactness": solution.peak_compactness,
        "peak_radius": solution.peak_radius,
        "outer_mass": solution.outer_mass,
        "outer_k_times_r_cubed": solution.outer_k_times_r_cubed,
        "finite_mass": solution.finite_mass,
        "no_initial_trapped_sphere": solution.no_initial_trapped_sphere,
    }


def _numerical_preflight(configuration: Mapping[str, Any]) -> dict[str, Any]:
    thresholds = configuration["numerical_thresholds"]
    convergence_minimum = float(thresholds["minimum_observed_outer_mass_convergence_order"])
    cross_method_maximum = float(thresholds["maximum_cross_method_peak_compactness_difference"])
    original_peak_maximum = float(thresholds["original_box_peak_compactness_numerical_upper_bound"])
    repair_minimum = float(thresholds["repair_scan_peak_compactness_minimum"])
    repair_maximum = float(thresholds["repair_scan_peak_compactness_maximum"])

    original_records: list[dict[str, Any]] = []
    for amplitude in configuration["original_amplitudes"]:
        method_records: dict[str, Any] = {}
        report_solutions = {}
        for method in configuration["methods"]:
            nested = tuple(
                solve_gr0_initial_slice(
                    PulseParameters(float(amplitude)),
                    step_count=steps,
                    method=method,
                )
                for steps in configuration["nested_steps"]
            )
            order = observed_convergence_order(*(item.outer_mass for item in nested))
            if order is None or not isfinite(order) or order < convergence_minimum:
                raise ValueError(
                    f"original amplitude {_text(amplitude)} {method} outer-mass convergence failed"
                )
            report = solve_gr0_initial_slice(
                PulseParameters(float(amplitude)),
                step_count=configuration["report_steps"],
                method=method,
            )
            if report.peak_compactness >= original_peak_maximum:
                raise ValueError("PROTO1 original-box numerical separation failed")
            method_records[method] = {
                "nested_outer_masses": [item.outer_mass for item in nested],
                "observed_outer_mass_convergence_order": order,
                "minimum_required_order": convergence_minimum,
                "report_solution": _solution_summary(report),
            }
            report_solutions[method] = report
        difference = abs(
            report_solutions["RK4"].peak_compactness
            - report_solutions["SSPRK3"].peak_compactness
        )
        if difference > cross_method_maximum:
            raise ValueError("PROTO1 cross-method compactness comparison failed")
        original_records.append(
            {
                "amplitude": amplitude,
                "methods": method_records,
                "cross_method_peak_compactness_difference": difference,
                "both_methods_below_one_over_2000": True,
            }
        )

    repair_records: list[dict[str, Any]] = []
    for amplitude in configuration["repair_amplitudes"]:
        solutions = {
            method: solve_gr0_initial_slice(
                PulseParameters(float(amplitude)),
                step_count=configuration["report_steps"],
                method=method,
            )
            for method in configuration["methods"]
        }
        difference = abs(
            solutions["RK4"].peak_compactness
            - solutions["SSPRK3"].peak_compactness
        )
        if difference > cross_method_maximum:
            raise ValueError("repair-scan cross-method compactness comparison failed")
        if not all(
            repair_minimum <= item.peak_compactness <= repair_maximum
            and item.no_initial_trapped_sphere
            and item.finite_mass
            for item in solutions.values()
        ):
            raise ValueError("repair-scan amplitude missed the frozen initial compactness window")
        repair_records.append(
            {
                "amplitude": amplitude,
                "methods": {
                    method: _solution_summary(solution)
                    for method, solution in solutions.items()
                },
                "cross_method_peak_compactness_difference": difference,
                "inside_original_initial_compactness_window": True,
            }
        )

    maximum_original = max(
        item["methods"][method]["report_solution"]["peak_compactness"]
        for item in original_records
        for method in configuration["methods"]
    )
    return {
        "interpretation": "floating_GR0_preflight_supporting_an_independent_exact_analytic_obstruction",
        "nested_step_counts": configuration["nested_steps"],
        "report_step_count": configuration["report_steps"],
        "methods": configuration["methods"],
        "original_box": {
            "records": original_records,
            "maximum_peak_compactness": maximum_original,
            "declared_numerical_upper_bound": original_peak_maximum,
            "all_candidates_below_declared_upper_bound": True,
            "all_methods_pass_outer_mass_convergence_threshold": True,
        },
        "premise_repair_scan": {
            "records": repair_records,
            "all_candidates_inside_original_initial_compactness_window": True,
            "uses_only_GR0_constraint_solutions": True,
            "FGCQR_outcomes_used": False,
            "does_not_freeze_or_authorize_a_successor_protocol": True,
        },
        "numerics_are_not_needed_for_the_entire_box_rejection": True,
    }


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    configuration = load_config(config_path)
    exact = _exact_controls(configuration)
    analytic = _analytic_bound(configuration)
    numerical = _numerical_preflight(configuration)
    protocol_hash = configuration["protocol"]["semantic_holdout_contract_sha256"]
    history = configuration["historical_checkpoint"]
    return _serial(
        {
            "schema_version": 1,
            "artifact_id": ARTIFACT_ID,
            "project_version": PROJECT_VERSION,
            "classification": "pre_holdout_GR0_protocol_obstruction_not_an_FGCQR_mechanism_result",
            "generated_by": _rel(Path(__file__)),
            "derivation_document": _rel(OWNER_DOCUMENT),
            "derivation_document_sha256": _sha(OWNER_DOCUMENT),
            "source_config_sha256": {_rel(config_path): configuration["source_sha256"]},
            "predecessor_sha256": {
                _rel(configuration["protocol_path"]): _sha(configuration["protocol_path"]),
            },
            "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
            "historical_preflight_checkpoint": history,
            "scope_bindings": {
                "target_protocol_artifact_id": SF1_PROTOCOL_V1_ARTIFACT_ID,
                "target_protocol_version": 1,
                "protocol_config_sha256": _sha(configuration["protocol_path"]),
                "protocol_semantic_holdout_contract_sha256": protocol_hash,
                "run1_artifact_id": RUN1_SYM1_ARTIFACT_ID,
                "historical_checkpoint_commit": history["commit"],
                "run1_config_sha256": history["file_sha256"][
                    HISTORICAL_CHECKPOINT["run1_config_path"]
                ],
                "run1_result_sha256": history["file_sha256"][
                    HISTORICAL_CHECKPOINT["run1_result_path"]
                ],
                "target_branch": "GR-0",
                "physical_equations": "unredefined_ACT1_VAR1",
                "FGCQR_holdout_outcomes_inspected": False,
            },
            "certificate_contract": {
                "artifact_specific_payload_validated": True,
                "canonical_reproduction_passed": True,
                "scope_bindings_verified": True,
                "historical_pre_PROTO2_checkpoint_verified": True,
                "frozen_PROTO1_preserved": True,
                "premise_failure_found_before_FGCQR_holdout": True,
                "outcome_data_not_used_to_revise_protocol": True,
            },
            "artifact_payload": {
                "protocol_diagnosis": dict(configuration["raw"]["protocol_diagnosis"]),
                "exact_constraint_specialization": exact,
                "analytic_entire_box_obstruction": analytic,
                "independent_numerical_preflight": numerical,
                "decision": {
                    "PROTO1_initial_data_contract_executable": False,
                    "PROTO1_original_amplitude_box_can_reach_declared_initial_compactness_floor_under_minimal_completion": False,
                    "PROTO1_calibration_or_holdout_may_execute": False,
                    "new_versioned_protocol_required": True,
                    "repair_must_freeze_Pi_phi": True,
                    "repair_may_use_only_the_reported_GR0_premise_scan_before_refreezing": True,
                },
                "epistemic_boundary": {
                    "protocol_design_route_closed": "PROTO1_as_written",
                    "FGCQR_action_or_mechanism_tested": False,
                    "FGCQR_action_or_mechanism_rejected": False,
                    "general_gradient_mechanism_tested": False,
                    "general_gradient_mechanism_rejected": False,
                    "scientific_role": "prevent_a_false_negative_or_non_event_from_being_misread_as_physics",
                },
            },
            "gate_status": {
                REQUIRED_TRUE_GATE: True,
                "PROTO1_initial_data_calibration_authorized": False,
                "PROTO1_FGCQR_holdout_execution_authorized": False,
                "PROTO2_premise_revision_required": True,
                "classical_spherical_diagnostic_authorized": False,
                "retained_EFT_evolution_authorized": False,
                "physical_transition_claim_authorized": False,
            },
            "nonclaims": dict(EXPECTED_CLAIMS),
        }
    )


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    try:
        value = json.loads(source, object_pairs_hook=_pairs)
    except (json.JSONDecodeError, ValueError) as exc:
        raise ValueError("ID0-PREF1 result must be valid unique-key JSON") from exc
    if not isinstance(value, dict) or source != _canonical(value):
        raise ValueError("ID0-PREF1 result must use canonical sorted JSON")
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    arguments.output.write_text(_canonical(record(arguments.config)), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
