#!/usr/bin/env python3
"""Reproduce PROTO12's immutable, outcome-neutral protocol freeze.

The certificate consumes only numerical-premise evidence.  It recomputes the
pairwise spectral-conditioning classification from RSP1's hash-bound terminal
checkpoint, freezes SRC4 as the unchanged GR-0 source instrument, and refuses
every fresh trajectory or physical claim.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_rsp1_pref12 as rsp1_reproduction  # noqa: E402
from scripts import run_fgc_gr0_calibration as inherited  # noqa: E402
from scripts import run_fgc_rsp1_resolution_study as rsp1_runner  # noqa: E402
from recursive_horizons.fgc.evolution.health_monitor import (  # noqa: E402
    compact_vacuum_buffer_window,
)
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    PROTO4_SPECTRAL_FIELD_ORDER,
    SpectralThresholds,
    proper_radial_profile,
    spectral_field_budgets,
)
from recursive_horizons.fgc.evolution.protocol_v11 import (  # noqa: E402
    validate_sf1_protocol_v11,
)
from recursive_horizons.fgc.evolution.protocol_v12 import (  # noqa: E402
    PROTO12_AMENDMENT,
    PROTO12_CLAIMS,
    validate_sf1_protocol_v12,
)
from recursive_horizons.fgc.evolution.rsp1_resolution_runtime import (  # noqa: E402
    RSP1_POINT_COUNTS,
    rsp1_gr0_common_event,
)
from recursive_horizons.fgc.evolution.spectral_sensitivity import (  # noqa: E402
    spectral_tail_sensitivity,
    three_grid_profile_convergence,
)
from recursive_horizons.fgc.evolution.spectral_sensitivity_v12 import (  # noqa: E402
    pairwise_resolved_or_saturated_admission,
)


ARTIFACT_ID = "FGC-1-PRO12-FRZ1"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-pro12-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-pro12-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-pro12-frz1.md"
PROTOCOL_MODULE = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/protocol_v12.py"
)
SPECTRAL_MODULE = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/spectral_sensitivity_v12.py"
)
LEGACY_SPECTRAL_MODULE = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/spectral_sensitivity.py"
)
IMPLEMENTATION = (
    Path(__file__).resolve(),
    PROTOCOL_MODULE,
    SPECTRAL_MODULE,
    LEGACY_SPECTRAL_MODULE,
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto4_admission.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/health_monitor.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/rsp1_resolution_runtime.py",
    REPOSITORY / "scripts/run_fgc_rsp1_resolution_study.py",
)
Q = Fraction


EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO12",
    "predecessor_protocol": "FGC-2-SF1-PROTO11",
    "calibration_diagnosis_artifact": "FGC-1-CAL8-PREF10",
    "resolution_evidence_artifact": "FGC-1-RSP1-PREF12",
    "source_instrument_artifact": "FGC-1-SRC4-VEC1",
    "calibration_branch": "GR-0",
    "holdout_branch": "FGC-QR",
    "freeze_role": (
        "premise_only_source_instrument_and_pairwise_spectral_protocol_certificate"
    ),
    "PROTO11_GR0_calibration_trajectory_inspected": True,
    "RSP1_GR0_resolution_trajectory_inspected": True,
    "PROTO12_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "mechanism_question_answered": False,
}

EXPECTED_PROOF_KEYS = {
    "PROTO12_overlay_must_validate_exactly",
    "checkpoint_must_be_ancestor_of_HEAD",
    "all_tracked_lineage_must_be_read_from_the_immutable_checkpoint",
    "immutable_blob_hashes_must_match_PROTO12_amendment",
    "immutable_PROTO11_must_validate_unchanged",
    "CAL8_RSP1_and_SRC4_must_remain_numerical_premise_evidence",
    "SRC4_equations_branch_thresholds_and_reference_map_must_remain_unchanged",
    "RSP1_terminal_checkpoint_must_recompute_before_conditioning",
    "RSP1_generic_all_field_raw_failure_must_remain_public",
    "every_raw_power_fraction_and_adjacent_ratio_must_be_serialized",
    "direct_ratio_below_one_quarter_must_remain_primary_on_each_pair",
    "pairwise_saturation_requires_both_pair_budgets_and_map_witnesses",
    "pairwise_saturation_requires_complete_three_grid_profile_contraction",
    "round_trip_residual_must_not_be_called_a_continuum_error_bound",
    "no_epsilon_floor_new_tolerance_asymptote_or_hidden_ratio_is_permitted",
    "RSP1_high_ladder_must_not_replace_the_PROTO12_calibration_ladder",
    "original_amplitudes_ladder_methods_schedule_stops_observables_and_claims_must_remain_unchanged",
    "new_output_namespaces_must_be_declared_but_not_accessed",
    "canonical_hash_bound_result_required",
}

EXPECTED_CLAIMS = {
    "PROTO12_frozen": True,
    "PROTO12_SRC4_source_backend_contract_frozen": True,
    "PROTO12_pairwise_spectral_classifier_derived": True,
    "RSP1_pairwise_conditioning_recomputed": True,
    "RSP1_generic_all_field_raw_spectral_admission_passed": False,
    **PROTO12_CLAIMS,
    "fresh_GR0_dynamic_calibration_completed": False,
}


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bytes_sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _canonical(value: Mapping[str, Any]) -> str:
    return json.dumps(
        inherited._serial(value),
        sort_keys=True,
        indent=2,
        ensure_ascii=True,
        allow_nan=False,
    ) + "\n"


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


def _load_canonical_json_bytes(payload: bytes, artifact_id: str) -> dict[str, Any]:
    value = json.loads(payload)
    if not isinstance(value, dict) or value.get("artifact_id") != artifact_id:
        raise ValueError(f"{artifact_id} canonical object differs")
    if payload != _canonical(value).encode("utf-8"):
        raise ValueError(f"{artifact_id} JSON is not canonical")
    return value


def _load_canonical_json(path: Path, artifact_id: str) -> dict[str, Any]:
    return _load_canonical_json_bytes(path.read_bytes(), artifact_id)


def _git(*args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *args],
        cwd=REPOSITORY,
        check=check,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def _checkpoint_blob(commit: str, relative: str, expected_sha: str) -> bytes:
    payload = _git("show", f"{commit}:{relative}").stdout
    if _bytes_sha(payload) != expected_sha:
        raise ValueError(f"immutable checkpoint blob differs: {relative}")
    return payload


def load_protocol(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        protocol = tomllib.load(handle)
    validate_sf1_protocol_v12(protocol)
    return protocol


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    if set(config) != {
        "schema_version",
        "artifact_id",
        "project_version",
        "metric_signature",
        "riemann_convention",
        "protocol_config",
        "scope",
        "immutable_lineage",
        "proof_contract",
        "claims",
    }:
        raise ValueError("PRO12-FRZ1 config root differs")
    if (
        config["schema_version"] != 1
        or config["artifact_id"] != ARTIFACT_ID
        or config["project_version"] != "0.11.0"
        or config["metric_signature"] != "-+++"
        or config["riemann_convention"] != "plus_partial_mu_gamma_nu"
        or config["protocol_config"]
        != "configs/fgc/fgc-2-sf1-protocol-v12.toml"
        or config["scope"] != EXPECTED_SCOPE
        or set(config["proof_contract"]) != EXPECTED_PROOF_KEYS
        or any(value is not True for value in config["proof_contract"].values())
        or config["claims"] != EXPECTED_CLAIMS
    ):
        raise ValueError("PRO12-FRZ1 config violates its fail-closed scope")

    amendment = PROTO12_AMENDMENT
    expected_lineage = {
        "checkpoint_commit": amendment["evidence_checkpoint_commit"],
        "predecessor_protocol_config": amendment["predecessor_protocol_config"],
        "predecessor_protocol_sha256": amendment["predecessor_protocol_sha256"],
        "calibration_diagnosis_config": amendment[
            "calibration_diagnosis_config"
        ],
        "calibration_diagnosis_config_sha256": amendment[
            "calibration_diagnosis_config_sha256"
        ],
        "calibration_diagnosis_result": amendment[
            "calibration_diagnosis_result"
        ],
        "calibration_diagnosis_result_sha256": amendment[
            "calibration_diagnosis_result_sha256"
        ],
        "resolution_evidence_config": amendment["resolution_evidence_config"],
        "resolution_evidence_config_sha256": amendment[
            "resolution_evidence_config_sha256"
        ],
        "resolution_evidence_result": amendment["resolution_evidence_result"],
        "resolution_evidence_result_sha256": amendment[
            "resolution_evidence_result_sha256"
        ],
        "resolution_checkpoint": (
            "runs/fgc-2-sf1/rsp1/amplitude3-spectrum/latest-checkpoint.npz"
        ),
        "resolution_checkpoint_sha256": (
            "0d5d613b6067b17dda35eb2c180a4540cd99b020edec7e4fd56f24f1ef4123f9"
        ),
        "source_instrument_config": amendment["source_instrument_config"],
        "source_instrument_config_sha256": amendment[
            "source_instrument_config_sha256"
        ],
        "source_instrument_result": amendment["source_instrument_result"],
        "source_instrument_result_sha256": amendment[
            "source_instrument_result_sha256"
        ],
    }
    if config["immutable_lineage"] != expected_lineage:
        raise ValueError("PRO12 immutable lineage differs")
    return config


def _validate_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("PRO12 evidence checkpoint is not an ancestor of HEAD")

    predecessor_bytes = _checkpoint_blob(
        commit,
        lineage["predecessor_protocol_config"],
        lineage["predecessor_protocol_sha256"],
    )
    predecessor = tomllib.loads(predecessor_bytes.decode("utf-8"))
    predecessor_validation = validate_sf1_protocol_v11(predecessor)

    calibration_config_bytes = _checkpoint_blob(
        commit,
        lineage["calibration_diagnosis_config"],
        lineage["calibration_diagnosis_config_sha256"],
    )
    calibration_config = tomllib.loads(calibration_config_bytes.decode("utf-8"))
    calibration_bytes = _checkpoint_blob(
        commit,
        lineage["calibration_diagnosis_result"],
        lineage["calibration_diagnosis_result_sha256"],
    )
    calibration = _load_canonical_json_bytes(
        calibration_bytes, "FGC-1-CAL8-PREF10"
    )
    if (
        calibration_config.get("artifact_id") != "FGC-1-CAL8-PREF10"
        or calibration.get("classification")
        != "post_calibration_PROTO11_evolved_source_arithmetic_and_direct_spectral_obstruction_diagnosis"
        or calibration.get("gate_status", {}).get(
            "PROTO11_campaign_terminated_normally"
        )
        is not True
        or calibration.get("gate_status", {}).get("PROTO11_GR0_case_eligible")
        is not False
        or calibration.get("gate_status", {}).get("PROTO12_frozen") is not False
        or calibration.get("gate_status", {}).get(
            "FGCQR_holdout_execution_authorized"
        )
        is not False
    ):
        raise ValueError("immutable CAL8 premise boundary differs")

    resolution_config_bytes = _checkpoint_blob(
        commit,
        lineage["resolution_evidence_config"],
        lineage["resolution_evidence_config_sha256"],
    )
    resolution_config = tomllib.loads(resolution_config_bytes.decode("utf-8"))
    resolution_bytes = _checkpoint_blob(
        commit,
        lineage["resolution_evidence_result"],
        lineage["resolution_evidence_result_sha256"],
    )
    resolution = _load_canonical_json_bytes(
        resolution_bytes, "FGC-1-RSP1-PREF12"
    )
    if (
        resolution_config.get("artifact_id") != "FGC-1-RSP1-PREF12"
        or resolution.get("gate_status", {}).get(
            "RSP1_terminal_checkpoint_recomputed"
        )
        is not True
        or resolution.get("gate_status", {}).get(
            "amplitude_three_spectral_veto_cleared_for_successor_design"
        )
        is not True
        or resolution.get("gate_status", {}).get(
            "RSP1_generic_all_field_raw_spectral_admission_passed"
        )
        is not False
        or resolution.get("gate_status", {}).get("PROTO12_frozen") is not False
        or resolution.get("gate_status", {}).get("GR0_case_eligible") is not False
    ):
        raise ValueError("immutable RSP1 premise boundary differs")

    source_config_bytes = _checkpoint_blob(
        commit,
        lineage["source_instrument_config"],
        lineage["source_instrument_config_sha256"],
    )
    source_config = tomllib.loads(source_config_bytes.decode("utf-8"))
    source_bytes = _checkpoint_blob(
        commit,
        lineage["source_instrument_result"],
        lineage["source_instrument_result_sha256"],
    )
    source = _load_canonical_json_bytes(source_bytes, "FGC-1-SRC4-VEC1")
    source_gates = source.get("gate_status", {})
    if (
        source_config.get("artifact_id") != "FGC-1-SRC4-VEC1"
        or any(
            source_gates.get(name) is not True
            for name in (
                "SRC4_tensor_contraction_evaluator_derived",
                "SRC4_exact_reference_and_exact_oracle_controls_passed",
                "SRC4_independent_SRC3_differential_controls_passed",
                "SRC4_CAP1_full_grid_strict_source_gate_passed",
                "SRC4_runtime_source_instrument_authorized_for_a_separately_frozen_GR0_protocol",
            )
        )
        or source_gates.get("PROTO12_frozen") is not False
        or source_gates.get("FGCQR_holdout_execution_authorized") is not False
    ):
        raise ValueError("immutable SRC4 premise boundary differs")

    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "predecessor_protocol_validation": predecessor_validation,
        "calibration_diagnosis": {
            "artifact_id": calibration["artifact_id"],
            "classification": calibration["classification"],
            "config_sha256": lineage["calibration_diagnosis_config_sha256"],
            "result_sha256": lineage["calibration_diagnosis_result_sha256"],
            "GR0_case_eligible": False,
            "is_mechanism_evidence": False,
        },
        "resolution_evidence": {
            "artifact_id": resolution["artifact_id"],
            "classification": resolution["classification"],
            "config_sha256": lineage["resolution_evidence_config_sha256"],
            "result_sha256": lineage["resolution_evidence_result_sha256"],
            "generic_all_field_raw_admission": False,
            "is_mechanism_evidence": False,
        },
        "source_instrument": {
            "artifact_id": source["artifact_id"],
            "classification": source["classification"],
            "config_sha256": lineage["source_instrument_config_sha256"],
            "result_sha256": lineage["source_instrument_result_sha256"],
            "equations_or_thresholds_changed": False,
            "is_mechanism_evidence": False,
        },
    }


def _conditioning_for_method(
    plan: Mapping[str, Any],
    members: Mapping[str, Any],
    method: str,
) -> dict[str, Any]:
    records = rsp1_runner._method_members(members, method)
    physical = plan["physical_inputs"]
    numerics = plan["numerics"]
    cutoff = float(Q(physical["cutoff_Lambda"]))
    measurement_maximum = float(Q(physical["measurement_radius_maximum"]))
    taper_fraction = float(Q(numerics["proper_spectral_taper_fraction"]))
    thresholds = SpectralThresholds()
    coordinate_time = float(Q(plan["numerics"]["final_coordinate_time"]))

    raw = rsp1_gr0_common_event(
        [item.state for item in records],
        [item.initial.grid for item in records],
        accepted_stage_counts=[
            item.transaction.state.accepted_stage_count for item in records
        ],
        method=method,
        coordinate_time=coordinate_time,
        cutoff=cutoff,
        measurement_radius_maximum=measurement_maximum,
        taper_fraction=taper_fraction,
        fixed_outer_rows=numerics["fixed_outer_rows"],
        point_batch_size=numerics["source_point_batch_size"],
        spectral_thresholds=thresholds,
    )

    profiles = []
    budgets = []
    sensitivities = []
    for member in records:
        profile = proper_radial_profile(
            member.state.u,
            member.state.q,
            member.initial.grid.coordinates,
            cutoff=cutoff,
            measurement_radius_maximum=measurement_maximum,
        )
        proper = profile.proper_coordinates
        window = compact_vacuum_buffer_window(
            proper,
            support_minimum=float(proper[0]),
            support_maximum=float(proper[-1]),
            taper_width=float((proper[-1] - proper[0]) * taper_fraction),
        )
        budget = dict(
            spectral_field_budgets(
                profile.deviations,
                proper,
                cutoff=cutoff,
                window=window,
                thresholds=thresholds,
            )
        )
        witness = {
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
        sensitivities.append(witness)

    convergence = three_grid_profile_convergence(RSP1_POINT_COUNTS, profiles)
    guarded = pairwise_resolved_or_saturated_admission(
        RSP1_POINT_COUNTS,
        budgets,
        sensitivities,
        convergence,
        thresholds=thresholds,
    )
    if (
        dict(raw.raw_spatial_spectral_admission.field_power_tail_ratios)
        != dict(guarded.field_power_tail_ratios)
        or dict(raw.raw_spatial_spectral_admission.derivative_power_tail_ratios)
        != dict(guarded.derivative_power_tail_ratios)
    ):
        raise ValueError(f"PROTO12 changed RSP1 raw spectral ratios: {method}")
    if (
        raw.constraint_admission.admission_passed is not True
        or raw.all_medium_and_fine_absolute_budgets_passed is not True
        or raw.profile_convergence.every_field_contracted is not True
        or raw.raw_spatial_spectral_admission.admission_passed is not False
        or guarded.medium_and_fine_individual_budgets_passed is not True
        or guarded.every_profile_contracted is not True
        or guarded.saturation_used is not True
        or guarded.every_tail_resolved_or_saturated is not True
        or guarded.admission_passed is not True
    ):
        raise ValueError(f"PROTO12 pairwise qualification failed: {method}")

    erasure_ratios = [
        witness.erasure_over_round_trip
        for by_grid in sensitivities
        for witness in by_grid.values()
        if witness.erasure_over_round_trip is not None
    ]
    if not erasure_ratios or not all(np.isfinite(value) for value in erasure_ratios):
        raise ValueError(f"PROTO12 erasure evidence is incomplete: {method}")
    return inherited._serial(
        {
            "method": method,
            "coordinate_time": coordinate_time,
            "point_counts": RSP1_POINT_COUNTS,
            "raw_generic_all_field_admission_passed": False,
            "constraint_admission_passed": True,
            "raw_spatial_spectral_admission": raw.raw_spatial_spectral_admission,
            "spectral_budgets": budgets,
            "spectral_tail_sensitivities": sensitivities,
            "profile_convergence": convergence,
            "pairwise_admission": guarded,
            "maximum_erasure_over_round_trip": max(erasure_ratios),
            "round_trip_is_a_continuum_error_bound": False,
            "diagnostic_saturation_is_physical_resolution_evidence": False,
        }
    )


def _recompute_rsp1_conditioning(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    checkpoint_path = REPOSITORY / lineage["resolution_checkpoint"]
    if not checkpoint_path.is_file() or _sha(checkpoint_path) != lineage[
        "resolution_checkpoint_sha256"
    ]:
        raise ValueError("RSP1 terminal checkpoint hash differs")

    rsp1_config_path = REPOSITORY / lineage["resolution_evidence_config"]
    rsp1_result_path = REPOSITORY / lineage["resolution_evidence_result"]
    if (
        _sha(rsp1_config_path) != lineage["resolution_evidence_config_sha256"]
        or _sha(rsp1_result_path) != lineage["resolution_evidence_result_sha256"]
    ):
        raise ValueError("current RSP1 certificate differs from immutable lineage")
    rsp1_config = rsp1_reproduction.load_config(rsp1_config_path)
    rsp1_result = rsp1_reproduction._load_json(
        rsp1_result_path, "FGC-1-RSP1-PREF12"
    )
    if rsp1_result != rsp1_reproduction.reproduce(rsp1_config_path):
        raise ValueError("stored RSP1/PREF12 result differs from reconstruction")

    plan_path = REPOSITORY / rsp1_config["run_plan_config"]
    authorization_path = REPOSITORY / rsp1_config["runtime_authorization_result"]
    plan, authorization = rsp1_runner._validate_authorization(
        plan_path,
        authorization_path,
        require_fresh_namespace=False,
    )
    metadata, members, _ = rsp1_runner._restore_checkpoint(
        checkpoint_path,
        plan=plan,
        authorization=authorization,
    )
    if metadata.get("terminal") is not True or metadata.get("amplitude") != "3":
        raise ValueError("RSP1 checkpoint is not the frozen terminal amplitude")
    expected_members = rsp1_result["artifact_payload"]["study_diagnosis"]["members"]
    for key, member in members.items():
        observed = inherited.array_content_sha256(
            member.state.u, member.state.p, member.state.q
        )
        if observed != expected_members[key]["state_sha256"]:
            raise ValueError(f"RSP1 restored terminal state differs: {key}")

    methods = {
        method: _conditioning_for_method(plan, members, method)
        for method in ("RK4", "SSPRK3")
    }
    maximum = max(
        method["maximum_erasure_over_round_trip"] for method in methods.values()
    )
    return {
        "source_checkpoint": _rel(checkpoint_path),
        "source_checkpoint_sha256": lineage["resolution_checkpoint_sha256"],
        "source_study_result_sha256": lineage["resolution_evidence_result_sha256"],
        "terminal_states_restored_and_hash_matched": True,
        "RSP1_generic_all_field_raw_admission_remains_public_and_false": True,
        "all_pairwise_conditioning_admissions_passed": True,
        "maximum_erasure_over_round_trip": maximum,
        "methods": methods,
        "classification": (
            "high_ladder_pairwise_diagnostic_conditioning_qualified_for_protocol_design"
        ),
        "is_calibration_or_mechanism_evidence": False,
    }


def record(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(path)
    protocol_path = (REPOSITORY / config["protocol_config"]).resolve()
    protocol = load_protocol(protocol_path)
    validation = validate_sf1_protocol_v12(protocol)
    lineage = _validate_lineage(config)
    conditioning = _recompute_rsp1_conditioning(config)
    roots = (
        REPOSITORY / validation["calibration_output_root"],
        REPOSITORY / validation["holdout_output_root"],
    )
    namespace = [{"path": _rel(root), "exists": root.exists()} for root in roots]
    if any(item["exists"] for item in namespace):
        raise ValueError("PROTO12 output namespace is not fresh at freeze time")

    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": config["project_version"],
        "metric_signature": config["metric_signature"],
        "riemann_convention": config["riemann_convention"],
        "classification": (
            "outcome_neutral_PROTO12_source_instrument_and_pairwise_spectral_freeze"
        ),
        "generated_by": _rel(Path(__file__)),
        "scope_bindings": dict(config["scope"]),
        "artifact_payload": {
            "immutable_lineage": lineage,
            "protocol_validation": validation,
            "RSP1_pairwise_conditioning_recomputation": conditioning,
            "namespace_precondition": {
                "records": namespace,
                "both_fresh_namespaces_absent": True,
                "freeze_created_no_namespace": True,
            },
            "decision": (
                "PROTO12_is_frozen_but_HLT10_runtime_fresh_calibration_and_every_"
                "candidate_or_physical_outcome_remain_closed"
            ),
            "epistemic_boundary": {
                "PROTO12_runtime_implemented": False,
                "PROTO12_trajectory_read": False,
                "fresh_state_constructed_under_PROTO12": False,
                "GR0_calibration_authorized": False,
                "GR0_calibration_completed": False,
                "SGBL_or_FGCQR_trajectory_read": False,
                "collapse_or_trapped_outcome_classified": False,
                "mechanism_question_answered": False,
                "round_trip_residual_used_as_continuum_error_bound": False,
            },
        },
        "gate_status": dict(config["claims"]),
        "implementation_sha256": {_rel(item): _sha(item) for item in IMPLEMENTATION},
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            _rel(path): _sha(path),
            _rel(protocol_path): _sha(protocol_path),
        },
        "nonclaims": {
            "PROTO12_runtime_implemented": False,
            "fresh_GR0_calibration_authorized": False,
            "fresh_GR0_calibration_completed": False,
            "trapped_sphere_observed": False,
            "SGBL_control_completed": False,
            "FGCQR_holdout_evolved": False,
            "regulator_activation_observed": False,
            "positive_Raychaudhuri_margin_observed": False,
            "FGCQR_mechanism_rejected": False,
            "general_gradient_route_rejected": False,
            "retained_EFT_validity": False,
            "physical_transition": False,
        },
    }


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    return _load_canonical_json(path, ARTIFACT_ID)


def verify_canonical(
    path: Path = DEFAULT_OUTPUT,
    config_path: Path = DEFAULT_CONFIG,
) -> None:
    if load_canonical_result(path) != record(config_path):
        raise ValueError("stored PRO12-FRZ1 result differs from reconstruction")


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
        print(f"verified {ARTIFACT_ID}: {output_path}")
        return
    payload = _canonical(record(config_path)).encode("utf-8")
    _atomic_write(output_path, payload)
    print(f"wrote {ARTIFACT_ID}: {output_path}")


if __name__ == "__main__":
    main()
