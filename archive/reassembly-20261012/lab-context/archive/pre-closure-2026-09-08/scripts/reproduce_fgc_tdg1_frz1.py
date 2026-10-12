#!/usr/bin/env python3
"""Reproduce the prospective TDG1 temporal-gate diagnosis freeze."""

from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import tomllib
from typing import Any, Callable, Mapping


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.tdg1_temporal_diagnosis import (  # noqa: E402
    TDG1_ABSOLUTE_AMPLITUDE_FLOOR,
    TDG1_CUTOFF,
    TDG1_DOWNSAMPLE_COUNTS,
    TDG1_FIELD_NAMES,
    TDG1_MAXIMUM_NESTED_TAIL_RATIO,
    TDG1_MINIMUM_HISTORY_DIFFERENCE_ORDER,
    TDG1_SAMPLE_COUNT,
    TDG1_SYNTHETIC_CONTROL_NAMES,
    TDG1_TAPER_FRACTION,
    TDG1_TRACER_COUNT,
    TDG1_VIEW_NAMES,
    synthetic_temporal_gate_controls,
)


ARTIFACT_ID = "FGC-1-TDG1-FRZ1"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg1-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-tdg1-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-tdg1-frz1.md"
DIAGNOSIS_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg1_temporal_diagnosis.py"
)
PROTO5_RUNTIME = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto5_runtime.py"
)
IMPLEMENTATION = (Path(__file__).resolve(), DIAGNOSIS_MODULE, PROTO5_RUNTIME)


EXPECTED_CLAIMS = {
    "TDG1_diagnostic_contract_frozen": True,
    "TDG1_synthetic_counterexample_confirmed": True,
    "inherited_normalized_raw_temporal_tail_ratio_is_valid_general_spatial_convergence_test": False,
    "TDG1_terminal_histories_diagnosed": False,
    "replacement_temporal_admission_defined": False,
    "PROTO14_frozen": False,
    "fresh_GR0_dynamic_calibration_completed": False,
    "GR0_case_eligible": False,
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


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _sha(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bytes_sha(payload: bytes) -> str:
    return sha256(payload).hexdigest()


def _verify_optional_raw_hash(path: Path, expected_hash: str) -> None:
    """Verify an ignored raw file when present, while supporting clean clones."""

    if not path.exists():
        return
    if not path.is_file() or _sha(path) != expected_hash:
        raise ValueError("TDG1 source checkpoint hash differs")


def _canonical(value: Mapping[str, Any]) -> str:
    return json.dumps(
        value, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False
    ) + "\n"


def _canonical_bytes(value: Mapping[str, Any]) -> bytes:
    return _canonical(value).encode("utf-8")


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


def _strict_keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(f"{name} keys differ")


def _fraction(value: str) -> float:
    if "/" in value:
        numerator, denominator = value.split("/", 1)
        return float(int(numerator) / int(denominator))
    return float(value)


def validate_config_data(config: Mapping[str, Any]) -> None:
    _strict_keys(
        "TDG1 freeze config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "source_result",
            "source_config",
            "source_checkpoint",
            "scope",
            "immutable_lineage",
            "frozen_diagnostic",
            "synthetic_controls",
            "outcome_classes",
            "proof_contract",
            "successor_boundary",
            "claims",
        },
    )
    if (
        config.get("schema_version") != 1
        or config.get("artifact_id") != ARTIFACT_ID
        or config.get("project_version") != PROJECT_VERSION
        or config.get("metric_signature") != "-+++"
        or config.get("riemann_convention") != "plus_partial_mu_gamma_nu"
        or config.get("source_result") != "results/fgc-1-cal10-pref15.json"
        or config.get("source_config") != "configs/fgc/fgc-1-cal10-pref15.toml"
        or config.get("source_checkpoint")
        != "runs/fgc-2-sf1/proto13/calibration/latest-checkpoint.npz"
    ):
        raise ValueError("TDG1 freeze top-level contract differs")

    scope = config.get("scope", {})
    lineage = config.get("immutable_lineage", {})
    diagnostic = config.get("frozen_diagnostic", {})
    controls = config.get("synthetic_controls", {})
    outcomes = config.get("outcome_classes", {})
    proof = config.get("proof_contract", {})
    successor = config.get("successor_boundary", {})
    claims = config.get("claims", {})
    if scope != {
        "target": "FGC-2-SF1-TDG1",
        "role": "prospective_no_trajectory_temporal_gate_diagnosis_freeze",
        "calibration_branch": "GR-0",
        "source_terminal_event": 63,
        "source_terminal_coordinate_time": "63/16",
        "source_terminal_outcome_already_bound_by_CAL10": True,
        "synthetic_controls_consume_no_campaign_histories": True,
        "TDG1_terminal_histories_read": False,
        "SGBL_trajectory_read": False,
        "FGCQR_trajectory_read": False,
        "replacement_temporal_admission_defined": False,
        "physical_or_candidate_question_answered": False,
    }:
        raise ValueError("TDG1 freeze scope differs")
    if lineage != {
        "checkpoint_commit": "fa56acf4d8c6297eb6f9b5dc8afdb5e2b5ce9020",
        "source_result_sha256": "1ebbdc7222baee8134a6d466589770ee445aee594f10630fc09aa6264cb4ac18",
        "source_config_sha256": "3661ba8cb6b933ab47daf8897c577b4631544605175dbe15d0c084542cf6c9eb",
        "source_checkpoint_sha256": "0432b266081b1ec1af3ac2a9efb99bcdd67dc8afef837241ec27fe5d37432f97",
        "historical_proto5_runtime_sha256": "607c9d5c353e37f8c68c032ec70df773c326335ba723a7db7317191bc7260545",
        "checkpoint_commit_must_be_ancestor_of_HEAD": True,
        "tracked_source_blobs_must_match_checkpoint_commit_and_worktree": True,
        "source_checkpoint_is_read_only": True,
        "source_checkpoint_may_not_resume": True,
    }:
        raise ValueError("TDG1 freeze immutable lineage differs")
    if diagnostic != {
        "sample_count": TDG1_SAMPLE_COUNT,
        "tracer_count": TDG1_TRACER_COUNT,
        "field_names": list(TDG1_FIELD_NAMES),
        "primary_method": "RK4",
        "primary_point_counts": [2049, 4097, 8193],
        "comparator_method": "SSPRK3",
        "comparator_point_counts": [4097, 8193, 16385],
        "cutoff": "16",
        "taper_fraction": "1/8",
        "historical_maximum_nested_tail_ratio": "1/4",
        "historical_threshold_must_not_change": True,
        "declared_absolute_FFT_amplitude_floor": "1e-30",
        "minimum_raw_history_difference_order": "3/2",
        "fixed_signal_views": list(TDG1_VIEW_NAMES),
        "fixed_resampling_sensitivity_counts": list(TDG1_DOWNSAMPLE_COUNTS),
        "cross_resolution_histories_use_common_proper_time_intersection": True,
        "round_trip_interpolation_debit_must_remain_public": True,
        "resampling_sensitivity_is_not_independent_temporal_convergence_evidence": True,
        "raw_history_difference_order_is_diagnostic_not_replacement_admission": True,
    }:
        raise ValueError("TDG1 frozen diagnostic differs")
    if (
        _fraction(diagnostic["cutoff"]) != TDG1_CUTOFF
        or _fraction(diagnostic["taper_fraction"]) != TDG1_TAPER_FRACTION
        or _fraction(diagnostic["historical_maximum_nested_tail_ratio"])
        != TDG1_MAXIMUM_NESTED_TAIL_RATIO
        or _fraction(diagnostic["declared_absolute_FFT_amplitude_floor"])
        != TDG1_ABSOLUTE_AMPLITUDE_FLOOR
        or _fraction(diagnostic["minimum_raw_history_difference_order"])
        != TDG1_MINIMUM_HISTORY_DIFFERENCE_ORDER
    ):
        raise ValueError("TDG1 numerical constants differ")
    if controls != {
        "control_names": list(TDG1_SYNTHETIC_CONTROL_NAMES),
        "zero_history_must_pass": True,
        "every_nonzero_smooth_history_must_pass_individual_budgets": True,
        "every_identical_nonzero_history_is_expected_to_expose_unit_nested_ratios": True,
        "amplitude_scaled_control_contracts_by_four_per_spatial_refinement": True,
        "amplitude_scaled_power_contracts_by_sixteen_per_spatial_refinement": True,
        "normalization_cancellation_must_remain_visible": True,
        "structural_counterexample_is_independent_of_PROTO13_history_values": True,
    }:
        raise ValueError("TDG1 synthetic control contract differs")
    if outcomes != {
        "structural_counterexample_confirmed": "inherited_normalized_raw_tail_ratio_is_not_a_general_spatial_convergence_observable",
        "structural_counterexample_not_confirmed": "diagnostic_implementation_or_premise_failure",
        "actual_history_diagnosis_pending": "no_raw_history_conclusion_at_freeze",
        "individual_budget_cause_must_remain_separate": True,
        "arithmetic_floor_window_sampling_interpolation_and_continuum_causes_may_not_be_merged": True,
        "no_history_view_may_retroactively_pass_PROTO13": True,
        "no_diagnostic_class_may_authorize_a_new_trajectory": True,
    }:
        raise ValueError("TDG1 outcome classes differ")
    if set(proof.values()) != {True}:
        raise ValueError("TDG1 proof contract must be all-of")
    if successor != {
        "TDG1_diagnostic_contract_frozen": True,
        "TDG1_actual_history_execution_authorized": True,
        "TDG1_actual_history_execution_completed": False,
        "replacement_temporal_gate_design_authorized": False,
        "PROTO14_frozen": False,
        "fresh_GR0_dynamic_calibration_completed": False,
        "GR0_case_eligible": False,
        "SGBL_execution_authorized": False,
        "FGCQR_holdout_execution_authorized": False,
        "DEF1_execution_authorized": False,
    }:
        raise ValueError("TDG1 successor boundary differs")
    if claims != EXPECTED_CLAIMS:
        raise ValueError("TDG1 claims differ")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    validate_config_data(config)
    return config


def _verify_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("TDG1 checkpoint commit is not an ancestor of HEAD")
    tracked = {
        config["source_result"]: lineage["source_result_sha256"],
        config["source_config"]: lineage["source_config_sha256"],
        _rel(PROTO5_RUNTIME): lineage["historical_proto5_runtime_sha256"],
    }
    blobs: dict[str, str] = {}
    for relative, expected_hash in tracked.items():
        committed = _git("show", f"{commit}:{relative}").stdout
        current = (REPOSITORY / relative).read_bytes()
        if committed != current or _bytes_sha(committed) != expected_hash:
            raise ValueError(f"TDG1 immutable tracked source differs: {relative}")
        blobs[relative] = expected_hash
    checkpoint = REPOSITORY / config["source_checkpoint"]
    _verify_optional_raw_hash(
        checkpoint, lineage["source_checkpoint_sha256"]
    )
    return {
        "checkpoint_commit": commit,
        "checkpoint_commit_is_ancestor_of_HEAD": True,
        "tracked_blob_sha256": blobs,
        "source_checkpoint_sha256": lineage["source_checkpoint_sha256"],
        "source_checkpoint_history_arrays_loaded": False,
        "source_checkpoint_resumed": False,
    }


def _expect_config_rejection(
    config: Mapping[str, Any],
    mutate: Callable[[dict[str, Any]], None],
) -> bool:
    attacked = deepcopy(config)
    mutate(attacked)
    try:
        validate_config_data(attacked)
    except (TypeError, ValueError):
        return True
    return False


def _mutation_controls(config: Mapping[str, Any]) -> dict[str, bool]:
    controls = {
        "historical_threshold_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["frozen_diagnostic"].__setitem__(
                "historical_maximum_nested_tail_ratio", "1/3"
            ),
        ),
        "synthetic_control_removal_rejected": _expect_config_rejection(
            config,
            lambda value: value["synthetic_controls"]["control_names"].pop(),
        ),
        "source_checkpoint_hash_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "source_checkpoint_sha256", "0" * 64
            ),
        ),
        "raw_history_read_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["scope"].__setitem__(
                "TDG1_terminal_histories_read", True
            ),
        ),
        "replacement_gate_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "replacement_temporal_admission_defined", True
            ),
        ),
        "GR0_eligibility_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "GR0_case_eligible", True
            ),
        ),
    }
    if not all(controls.values()):
        raise RuntimeError("TDG1 mutation controls did not all fail closed")
    return controls


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    lineage = _verify_lineage(config)
    synthetic = synthetic_temporal_gate_controls()
    if synthetic["structural_counterexample_confirmed"] is not True:
        raise RuntimeError("TDG1 structural counterexample did not reproduce")
    mutations = _mutation_controls(config)
    payload = {
        "immutable_lineage": lineage,
        "frozen_diagnostic": {
            **config["frozen_diagnostic"],
            "diagnosis_module_sha256": _sha(DIAGNOSIS_MODULE),
        },
        "synthetic_preflight": synthetic,
        "mutation_controls": mutations,
        "claim_boundary": {
            "raw_terminal_histories_diagnosed": False,
            "historical_PROTO13_result_reclassified": False,
            "historical_threshold_changed": False,
            "replacement_temporal_admission_defined": False,
            "new_trajectory_authorized": False,
            "candidate_branch_opened": False,
        },
        "successor_boundary": config["successor_boundary"],
    }
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
        "classification": "prospective_TDG1_temporal_gate_diagnosis_freeze",
        "artifact_payload": payload,
        "gate_status": config["claims"],
        "nonclaims": {
            key: value
            for key, value in config["claims"].items()
            if value is False
        },
        "source_config_sha256": {_rel(config_path): _sha(config_path)},
        "implementation_sha256": {
            _rel(path): _sha(path) for path in IMPLEMENTATION
        },
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": (
            _sha(OWNER_DOCUMENT) if OWNER_DOCUMENT.is_file() else None
        ),
        "generated_by": _rel(Path(__file__)),
    }


def _validate_stored_record(
    value: Mapping[str, Any], config: Mapping[str, Any]
) -> None:
    if (
        value.get("artifact_id") != ARTIFACT_ID
        or value.get("classification")
        != "prospective_TDG1_temporal_gate_diagnosis_freeze"
        or value.get("gate_status") != EXPECTED_CLAIMS
        or value.get("artifact_payload", {})
        .get("synthetic_preflight", {})
        .get("structural_counterexample_confirmed")
        is not True
        or value.get("artifact_payload", {})
        .get("immutable_lineage", {})
        .get("source_checkpoint_history_arrays_loaded")
        is not False
        or value.get("artifact_payload", {})
        .get("claim_boundary", {})
        .get("new_trajectory_authorized")
        is not False
        or value.get("source_config_sha256")
        != {_rel(DEFAULT_CONFIG): _sha(DEFAULT_CONFIG)}
    ):
        raise ValueError("TDG1 stored record claim boundary differs")
    if config["claims"] != EXPECTED_CLAIMS:
        raise ValueError("TDG1 stored record config claims differ")


def verify_canonical(
    config_path: Path = DEFAULT_CONFIG,
    output_path: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    config = load_config(config_path)
    observed = record(config_path)
    raw = output_path.read_bytes()
    stored = json.loads(raw)
    if not isinstance(stored, dict) or raw != _canonical_bytes(stored):
        raise ValueError("TDG1 stored result is not canonical JSON")
    _validate_stored_record(stored, config)
    if stored != observed:
        raise ValueError("TDG1 stored result differs from reproduction")
    return observed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write and args.verify:
        raise SystemExit("choose either --write or --verify")
    value = record(args.config)
    if args.write:
        _atomic_write(args.output, _canonical_bytes(value))
        print(f"wrote {_rel(args.output)}")
    elif args.verify:
        verify_canonical(args.config, args.output)
        print(f"PASS {ARTIFACT_ID}")
    else:
        print(_canonical(value), end="")


if __name__ == "__main__":
    main()
