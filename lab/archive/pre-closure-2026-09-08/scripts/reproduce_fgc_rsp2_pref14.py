#!/usr/bin/env python3
"""Reproduce PREF14's independent post-run RSP2 certificate."""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import tomllib
from typing import Any, Mapping, Sequence

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    UniformRadialGrid,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.rsp2_checkpoint import (  # noqa: E402
    load_rsp2_cal9_predecessors,
)
from recursive_horizons.fgc.evolution.rsp2_pref14_diagnosis import (  # noqa: E402
    PREF14_POINT_COUNTS,
    PREF14_TARGET_TIME,
    diagnose_rsp2_pref14_endpoint,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-RSP2-PREF14"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-rsp2-pref14.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-rsp2-pref14.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-rsp2-pref14.md"
DIAGNOSIS_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/rsp2_pref14_diagnosis.py"
)
IMPLEMENTATION = (
    Path(__file__).resolve(),
    DIAGNOSIS_MODULE,
    REPOSITORY / "src/recursive_horizons/fgc/evolution/rsp2_checkpoint.py",
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/cal4_common_event_diagnosis.py",
)


EXPECTED_CLAIMS = {
    "RSP2_study_terminated_normally": True,
    "RSP2_new_member_reached_endpoint": True,
    "RSP2_terminal_checkpoint_recomputed": True,
    "RSP2_source_retry_total_zero": True,
    "RSP2_CFL_retry_total_zero": True,
    "RSP2_outcome_read": True,
    "RSP2_target_order_cleared": True,
    "RSP2_complete_constraint_admission_passed": True,
    "RSP2_target_preasymptotic_on_tested_ladder": True,
    "RSP2_target_persistent_through_tested_ladder": False,
    "PROTO13_design_may_begin": True,
    "PROTO13_frozen": False,
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


def _canonical(value: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(
            value,
            indent=2,
            sort_keys=True,
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _canonical_line(value: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _serial(value: Any) -> Any:
    if is_dataclass(value):
        return _serial(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_serial(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


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


def _load_json(path: Path, artifact_id: str | None = None) -> dict[str, Any]:
    raw = path.read_bytes()
    value = json.loads(raw)
    if not isinstance(value, dict) or raw != _canonical(value):
        raise ValueError(f"{_rel(path)} must contain canonical JSON")
    if artifact_id is not None and value.get("artifact_id") != artifact_id:
        raise ValueError(f"{_rel(path)} artifact identifier differs")
    return value


def _load_events(path: Path) -> tuple[list[dict[str, Any]], bytes]:
    raw = path.read_bytes()
    if not raw.endswith(b"\n"):
        raise ValueError("RSP2 event log must end in a newline")
    records: list[dict[str, Any]] = []
    rebuilt = bytearray()
    for line in raw.splitlines():
        value = json.loads(line)
        if not isinstance(value, dict):
            raise ValueError("RSP2 event records must be objects")
        records.append(value)
        rebuilt.extend(_canonical_line(value))
    if bytes(rebuilt) != raw:
        raise ValueError("RSP2 event log is not canonical JSONL")
    return records, raw


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


def _as_float(value: str) -> float:
    return float(Q(value))


def validate_config_data(config: Mapping[str, Any]) -> None:
    _strict_keys(
        "PREF14 config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "runtime_authorization_result",
            "run_plan_config",
            "predecessor_checkpoint",
            "study_manifest",
            "study_event_log",
            "study_checkpoint",
            "study_result",
            "scope",
            "immutable_study",
            "terminal_member",
            "constraint_reduction",
            "proof_contract",
            "successor_boundary",
            "claims",
        },
    )
    expected_paths = {
        "runtime_authorization_result": "results/fgc-1-rsp2-frz1.json",
        "run_plan_config": "configs/fgc/fgc-1-rsp2-run1.toml",
        "predecessor_checkpoint": (
            "runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz"
        ),
        "study_manifest": (
            "runs/fgc-2-sf1/rsp2/amplitude3-ssprk3-momentum/manifest.json"
        ),
        "study_event_log": (
            "runs/fgc-2-sf1/rsp2/amplitude3-ssprk3-momentum/events.jsonl"
        ),
        "study_checkpoint": (
            "runs/fgc-2-sf1/rsp2/amplitude3-ssprk3-momentum/latest-checkpoint.npz"
        ),
        "study_result": (
            "runs/fgc-2-sf1/rsp2/amplitude3-ssprk3-momentum/study-result.json"
        ),
    }
    if (
        config.get("schema_version") != 1
        or config.get("artifact_id") != ARTIFACT_ID
        or config.get("project_version") != PROJECT_VERSION
        or config.get("metric_signature") != "-+++"
        or config.get("riemann_convention") != "plus_partial_mu_gamma_nu"
        or any(config.get(key) != value for key, value in expected_paths.items())
    ):
        raise ValueError("PREF14 top-level contract differs")

    scope = config["scope"]
    _strict_keys(
        "PREF14 scope",
        scope,
        {
            "target_protocol",
            "role",
            "calibration_branch",
            "amplitude",
            "method",
            "GR0_numerical_trajectory_read",
            "SGBL_trajectory_read",
            "FGCQR_trajectory_read",
            "collapse_or_trapped_outcome_classified",
            "mechanism_question_answered",
        },
    )
    if scope != {
        "target_protocol": "FGC-1-RSP2-FRZ1",
        "role": "post_run_checkpoint_recomputed_one_finer_grid_constraint_order_result",
        "calibration_branch": "GR-0",
        "amplitude": "3",
        "method": "SSPRK3",
        "GR0_numerical_trajectory_read": True,
        "SGBL_trajectory_read": False,
        "FGCQR_trajectory_read": False,
        "collapse_or_trapped_outcome_classified": False,
        "mechanism_question_answered": False,
    }:
        raise ValueError("PREF14 scope differs")

    immutable = config["immutable_study"]
    _strict_keys(
        "PREF14 immutable study",
        immutable,
        {
            "authorization_commit",
            "authorization_sha256",
            "study_id",
            "manifest_sha256",
            "event_log_sha256",
            "checkpoint_sha256",
            "study_result_sha256",
            "event_log_line_count",
            "event_log_byte_count",
            "accepted_boundary_count",
            "terminal_classification",
            "elapsed_wall_seconds",
        },
    )
    if (
        immutable.get("authorization_commit")
        != "81412c1b6b86648a289e865d8ab40b9bbab26370"
        or immutable.get("authorization_sha256")
        != "a6c147d7b369387e70318ee42eb56c97332702a227fd93d7ad20a1cfe2c1aa22"
        or immutable.get("study_id")
        != "a3a7b331d4fa72ff2e37fd674256a61b216c809539c7a0c4b154c76961355257"
        or immutable.get("event_log_line_count") != 25
        or immutable.get("event_log_byte_count") != 15060
        or immutable.get("accepted_boundary_count") != 23
        or immutable.get("terminal_classification")
        != "completed_target_and_complete_constraint_pass"
        or _as_float(immutable["elapsed_wall_seconds"]) != 8078.001129082986
    ):
        raise ValueError("PREF14 immutable study contract differs")

    member = config["terminal_member"]
    _strict_keys(
        "PREF14 terminal member",
        member,
        {
            "key",
            "point_count",
            "coordinate_time",
            "step_index",
            "transaction_serial",
            "accepted_stage_count",
            "source_retry_count",
            "CFL_retry_count",
            "state_sha256",
            "accumulated_characteristic_distance",
            "remaining_boundary_margin",
        },
    )
    if (
        member.get("key") != "SSPRK3-16385"
        or member.get("point_count") != 16385
        or _as_float(member["coordinate_time"]) != PREF14_TARGET_TIME
        or member.get("step_index") != 1771
        or member.get("transaction_serial") != 7084
        or member.get("accepted_stage_count") != 7084
        or member.get("source_retry_count") != 0
        or member.get("CFL_retry_count") != 0
        or member.get("state_sha256")
        != "91e64d0864390cb4ac63cfdbdbbcb8ca5ad6c2f47cdd4d6741e05ef560adad0d"
    ):
        raise ValueError("PREF14 terminal member differs")

    reduction = config["constraint_reduction"]
    _strict_keys(
        "PREF14 constraint reduction",
        reduction,
        {
            "method",
            "point_counts",
            "accepted_stage_counts",
            "state_sha256",
            "coordinate_time",
            "fixed_outer_rows",
            "coarsest_guard_maximum",
            "finest_guard_maximum",
            "target_component",
            "minimum_finest_pair_order",
            "target_raw_owned_norms",
            "target_accumulated_roundoff_enclosures",
            "target_effective_norms_for_order_only",
            "target_adjacent_pair_orders",
            "target_adjacent_pair_margins",
            "minimum_finite_finest_pair_order",
            "component_finest_pair_orders",
        },
    )
    expected_hashes = [
        "ab96602c0f05633375710efee2be01be84e24c32a31d07d6ef23390788a72162",
        "c46abc3f0f72fa2a13c1fee7bb329b94698cc2c26003406ae3eb1592f11b648c",
        "91e64d0864390cb4ac63cfdbdbbcb8ca5ad6c2f47cdd4d6741e05ef560adad0d",
    ]
    if (
        reduction.get("method") != "SSPRK3"
        or reduction.get("point_counts") != list(PREF14_POINT_COUNTS)
        or reduction.get("accepted_stage_counts") != [2524, 3948, 7084]
        or reduction.get("state_sha256") != expected_hashes
        or _as_float(reduction["coordinate_time"]) != PREF14_TARGET_TIME
        or reduction.get("fixed_outer_rows") != 4
        or _as_float(reduction["coarsest_guard_maximum"]) != 0.1
        or _as_float(reduction["finest_guard_maximum"]) != 0.005
        or reduction.get("target_component") != "radial_momentum"
        or _as_float(reduction["minimum_finest_pair_order"]) != 1.5
    ):
        raise ValueError("PREF14 constraint reduction contract differs")
    expected_vectors = {
        "target_raw_owned_norms": [
            0.0060603980047813765,
            0.0021440219481765766,
            0.000542836914675905,
        ],
        "target_accumulated_roundoff_enclosures": [
            2.296474121976644e-09,
            3.5915945773012936e-09,
            6.4437699620611966e-09,
        ],
        "target_effective_norms_for_order_only": [
            0.0060603957083072545,
            0.0021440183565819993,
            0.0005428304709059429,
        ],
        "target_adjacent_pair_orders": [
            1.4990947384363291,
            1.9817436463819333,
        ],
        "target_adjacent_pair_margins": [
            -0.0009052615636708783,
            0.48174364638193334,
        ],
    }
    for key, expected in expected_vectors.items():
        if [_as_float(value) for value in reduction[key]] != expected:
            raise ValueError(f"PREF14 {key} differs")
    if (
        _as_float(reduction["minimum_finite_finest_pair_order"])
        != 1.8428133958883794
        or {
            key: _as_float(value)
            for key, value in reduction["component_finest_pair_orders"].items()
        }
        != {
            "Hamiltonian": 2.1246424153634873,
            "radial_momentum": 1.9817436463819333,
            "gauge_t": 1.8516167425060845,
            "gauge_r": 1.8428133958883794,
        }
    ):
        raise ValueError("PREF14 component-order contract differs")

    proof = config["proof_contract"]
    if not proof or any(value is not True for value in proof.values()):
        raise ValueError("PREF14 proof contract differs")
    successor = config["successor_boundary"]
    if successor != {
        "PROTO13_design_may_begin": True,
        "PROTO13_frozen": False,
        "PROTO13_must_use_method_owned_resolution_ladders": True,
        "amplitude_three_is_the_only_prospectively_designated_calibration_case": True,
        "fresh_GR0_dynamic_calibration_completed": False,
        "GR0_case_eligible": False,
        "SGBL_execution_authorized": False,
        "FGCQR_holdout_execution_authorized": False,
        "DEF1_execution_authorized": False,
    }:
        raise ValueError("PREF14 successor boundary differs")
    if config["claims"] != EXPECTED_CLAIMS:
        raise ValueError("PREF14 claim boundary differs")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    validate_config_data(config)
    return config


def _validate_commit_lineage(
    config: Mapping[str, Any],
    authorization_path: Path,
    authorization: Mapping[str, Any],
    manifest: Mapping[str, Any],
) -> None:
    immutable = config["immutable_study"]
    commit = immutable["authorization_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("RSP2 authorization commit is not an ancestor of HEAD")
    authorization_blob = _git("show", f"{commit}:{_rel(authorization_path)}").stdout
    if (
        _bytes_sha(authorization_blob) != immutable["authorization_sha256"]
        or _sha(authorization_path) != immutable["authorization_sha256"]
        or authorization_blob != _canonical(authorization)
    ):
        raise ValueError("immutable RSP2 authorization differs")
    implementation = manifest.get("implementation_sha256")
    if not isinstance(implementation, dict) or not implementation:
        raise ValueError("RSP2 launch implementation ledger is absent")
    for relative, expected in implementation.items():
        blob = _git("show", f"{commit}:{relative}").stdout
        path = REPOSITORY / relative
        if (
            _bytes_sha(blob) != expected
            or not path.is_file()
            or _sha(path) != expected
        ):
            raise ValueError(f"immutable RSP2 implementation differs: {relative}")


def _validate_manifest(
    config: Mapping[str, Any],
    authorization: Mapping[str, Any],
    manifest: Mapping[str, Any],
) -> None:
    immutable = config["immutable_study"]
    frozen = authorization["artifact_payload"]["frozen_run_plan"]
    if (
        manifest.get("schema_version") != 1
        or manifest.get("runner_id") != "FGC-1-RSP2-RUN1-RUNNER"
        or manifest.get("study_id") != immutable["study_id"]
        or manifest.get("git_commit") != immutable["authorization_commit"]
        or manifest.get("authorization_path")
        != config["runtime_authorization_result"]
        or manifest.get("authorization_sha256")
        != immutable["authorization_sha256"]
        or manifest.get("plan_path") != config["run_plan_config"]
        or manifest.get("plan_sha256") != frozen["run_plan_sha256"]
        or manifest.get("input_manifest_sha256")
        != frozen["input_manifest_sha256"]
        or manifest.get("predecessor_checkpoint_sha256")
        != authorization["predecessor_sha256"][config["predecessor_checkpoint"]]
        or manifest.get("amplitude") != "3"
        or manifest.get("method") != "SSPRK3"
        or manifest.get("new_point_count") != 16385
        or manifest.get("outcome_fields_present_at_creation") is not False
    ):
        raise ValueError("RSP2 launch manifest differs")


def _validate_events(
    config: Mapping[str, Any],
    records: Sequence[Mapping[str, Any]],
    payload: bytes,
) -> None:
    immutable = config["immutable_study"]
    if (
        len(records) != 25
        or len(payload) != immutable["event_log_byte_count"]
        or _bytes_sha(payload) != immutable["event_log_sha256"]
    ):
        raise ValueError("RSP2 event-log identity differs")
    initial = records[0]
    if (
        initial.get("event_type") != "initial_premise_event"
        or initial.get("boundary_index") != 0
        or initial.get("coordinate_time") != 0.0
        or initial.get("trajectory_advanced") is not False
    ):
        raise ValueError("RSP2 initial event differs")
    for index, record in enumerate(records[1:24], start=1):
        member = record.get("member", {})
        if (
            record.get("event_type") != "accepted_progress_boundary"
            or record.get("boundary_index") != index
            or np.float64(record.get("coordinate_time")).tobytes()
            != np.float64(index / 16.0).tobytes()
            or member.get("source_retry_count") != 0
            or member.get("CFL_retry_count") != 0
        ):
            raise ValueError(f"RSP2 accepted boundary differs: {index}")
    terminal = records[-1]
    if (
        terminal.get("event_type") != "evolved_RSP2_endpoint"
        or terminal.get("boundary_index") != 23
        or np.float64(terminal.get("coordinate_time")).tobytes()
        != np.float64(PREF14_TARGET_TIME).tobytes()
        or terminal.get("classification")
        != "completed_target_and_complete_constraint_pass"
    ):
        raise ValueError("RSP2 terminal event differs")


def _readonly_state_array(value: object, *, point_count: int) -> np.ndarray:
    answer = np.ascontiguousarray(np.asarray(value, dtype=np.float64)).copy()
    if answer.shape != (point_count, 6) or not np.all(np.isfinite(answer)):
        raise ValueError("RSP2 terminal state array differs")
    answer.setflags(write=False)
    return answer


def _restore_terminal_checkpoint(
    config: Mapping[str, Any],
    manifest: Mapping[str, Any],
    event_payload: bytes,
    result: Mapping[str, Any],
) -> tuple[dict[str, Any], EvolutionState]:
    checkpoint_path = REPOSITORY / config["study_checkpoint"]
    expected_arrays = {
        "metadata_utf8",
        "event_log_utf8",
        "SSPRK3_16385_u",
        "SSPRK3_16385_p",
        "SSPRK3_16385_q",
        "SSPRK3_16385_tracer_positions",
        "SSPRK3_16385_tracer_proper_times",
        "SSPRK3_16385_event_proper_times",
        "SSPRK3_16385_event_fields",
    }
    with np.load(checkpoint_path, allow_pickle=False) as archive:
        if set(archive.files) != expected_arrays:
            raise ValueError("RSP2 checkpoint array inventory differs")
        metadata_payload = bytes(archive["metadata_utf8"])
        metadata = json.loads(metadata_payload)
        if not isinstance(metadata, dict) or metadata_payload != _canonical(metadata):
            raise ValueError("RSP2 checkpoint metadata are not canonical")
        embedded_events = bytes(archive["event_log_utf8"])
        state = EvolutionState(
            _readonly_state_array(archive["SSPRK3_16385_u"], point_count=16385),
            _readonly_state_array(archive["SSPRK3_16385_p"], point_count=16385),
            _readonly_state_array(archive["SSPRK3_16385_q"], point_count=16385),
        )
        auxiliary_shapes = {
            "SSPRK3_16385_tracer_positions": (48,),
            "SSPRK3_16385_tracer_proper_times": (48,),
            "SSPRK3_16385_event_proper_times": (24, 48),
            "SSPRK3_16385_event_fields": (24, 48, 6),
        }
        if any(
            archive[key].shape != shape or not np.all(np.isfinite(archive[key]))
            for key, shape in auxiliary_shapes.items()
        ):
            raise ValueError("RSP2 checkpoint tracer evidence differs")
    member = metadata.get("members", {}).get("SSPRK3-16385", {})
    if (
        metadata.get("schema_version") != 1
        or metadata.get("campaign_id") != manifest["study_id"]
        or metadata.get("study_id") != manifest["study_id"]
        or metadata.get("manifest_sha256") != _bytes_sha(_canonical(manifest))
        or metadata.get("amplitude") != "3"
        or metadata.get("amplitude_index") != 0
        or metadata.get("completed_common_event_index") != 23
        or metadata.get("terminal") is not True
        or metadata.get("terminal_result") != result
        or metadata.get("event_log")
        != {
            "byte_count": len(event_payload),
            "line_count": event_payload.count(b"\n"),
            "sha256": _bytes_sha(event_payload),
        }
        or embedded_events != event_payload
        or member.get("state_sha256") != config["terminal_member"]["state_sha256"]
        or member.get("time") != PREF14_TARGET_TIME
        or member.get("step_index") != 1771
        or member.get("transaction_serial") != 7084
        or member.get("source_retry_count") != 0
        or member.get("CFL_retry_count") != 0
    ):
        raise ValueError("RSP2 terminal checkpoint metadata differ")
    if array_content_sha256(state.u, state.p, state.q) != member["state_sha256"]:
        raise ValueError("RSP2 terminal checkpoint state hash differs")
    return metadata, state


def _mutation_controls(
    config: Mapping[str, Any],
    states: tuple[EvolutionState, EvolutionState, EvolutionState],
    grids: tuple[UniformRadialGrid, UniformRadialGrid, UniformRadialGrid],
    stage_counts: tuple[int, int, int],
    state_hashes: tuple[str, str, str],
) -> dict[str, bool]:
    reduction = config["constraint_reduction"]
    kwargs = {
        "accepted_stage_counts": stage_counts,
        "expected_state_sha256": state_hashes,
        "coordinate_time": PREF14_TARGET_TIME,
        "target_component": "radial_momentum",
        "minimum_finest_pair_order": 1.5,
        "fixed_outer_rows": reduction["fixed_outer_rows"],
        "coarsest_guard_maximum": _as_float(reduction["coarsest_guard_maximum"]),
        "finest_guard_maximum": _as_float(reduction["finest_guard_maximum"]),
    }

    threshold = diagnose_rsp2_pref14_endpoint(
        states,
        grids,
        **{**kwargs, "minimum_finest_pair_order": 2.0},
    )
    threshold_visible = not threshold.target_finest_pair_passed

    def rejected(callable_object: Any) -> bool:
        try:
            callable_object()
        except (TypeError, ValueError):
            return True
        return False

    time_visible = rejected(
        lambda: diagnose_rsp2_pref14_endpoint(
            states, grids, **{**kwargs, "coordinate_time": 1.5}
        )
    )
    component_visible = rejected(
        lambda: diagnose_rsp2_pref14_endpoint(
            states, grids, **{**kwargs, "target_component": "Hamiltonian"}
        )
    )
    wrong_grids = (grids[0], grids[1], UniformRadialGrid(0.0, 128.0, 8193))
    grid_visible = rejected(
        lambda: diagnose_rsp2_pref14_endpoint(states, wrong_grids, **kwargs)
    )

    changed_q = states[2].q.copy()
    changed_q.view(np.uint64)[31, 4] ^= np.uint64(1)
    changed_state = EvolutionState(states[2].u, states[2].p, changed_q)
    state_visible = rejected(
        lambda: diagnose_rsp2_pref14_endpoint(
            (states[0], states[1], changed_state), grids, **kwargs
        )
    )
    stage_visible = rejected(
        lambda: diagnose_rsp2_pref14_endpoint(
            states,
            grids,
            **{**kwargs, "accepted_stage_counts": (2524, 3948, 7085)},
        )
    )
    claim_mutation = deepcopy(config)
    claim_mutation["claims"]["FGCQR_holdout_execution_authorized"] = True
    claim_visible = rejected(lambda: validate_config_data(claim_mutation))
    return {
        "threshold": threshold_visible,
        "grid": grid_visible,
        "time": time_visible,
        "component": component_visible,
        "claim": claim_visible,
        "one_bit_state": state_visible,
        "accepted_stage_count": stage_visible,
    }


def reproduce(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    immutable = config["immutable_study"]
    raw_paths = {
        "manifest": REPOSITORY / config["study_manifest"],
        "event_log": REPOSITORY / config["study_event_log"],
        "checkpoint": REPOSITORY / config["study_checkpoint"],
        "study_result": REPOSITORY / config["study_result"],
    }
    if any(not path.is_file() for path in raw_paths.values()):
        raise FileNotFoundError("complete RSP2 raw bundle is required to reproduce PREF14")
    expected_raw_hashes = {
        "manifest": immutable["manifest_sha256"],
        "event_log": immutable["event_log_sha256"],
        "checkpoint": immutable["checkpoint_sha256"],
        "study_result": immutable["study_result_sha256"],
    }
    observed_raw_hashes = {key: _sha(path) for key, path in raw_paths.items()}
    if observed_raw_hashes != expected_raw_hashes:
        raise ValueError("RSP2 raw bundle hash differs")

    authorization_path = REPOSITORY / config["runtime_authorization_result"]
    authorization = _load_json(authorization_path, "FGC-1-RSP2-FRZ1")
    manifest = _load_json(raw_paths["manifest"])
    result = _load_json(raw_paths["study_result"])
    event_records, event_payload = _load_events(raw_paths["event_log"])
    _validate_commit_lineage(config, authorization_path, authorization, manifest)
    _validate_manifest(config, authorization, manifest)
    _validate_events(config, event_records, event_payload)
    metadata, fine_state = _restore_terminal_checkpoint(
        config, manifest, event_payload, result
    )

    predecessor_hash = authorization["predecessor_sha256"][
        config["predecessor_checkpoint"]
    ]
    predecessors = load_rsp2_cal9_predecessors(
        REPOSITORY / config["predecessor_checkpoint"],
        expected_sha256=predecessor_hash,
        outer_radius=128.0,
    )
    states = (predecessors.states[0], predecessors.states[1], fine_state)
    grids = (
        predecessors.grids[0],
        predecessors.grids[1],
        UniformRadialGrid(0.0, 128.0, 16385),
    )
    stage_counts = (2524, 3948, 7084)
    state_hashes = tuple(config["constraint_reduction"]["state_sha256"])
    reduction = config["constraint_reduction"]
    diagnosis = diagnose_rsp2_pref14_endpoint(
        states,
        grids,
        accepted_stage_counts=stage_counts,
        expected_state_sha256=state_hashes,
        coordinate_time=_as_float(reduction["coordinate_time"]),
        target_component=reduction["target_component"],
        minimum_finest_pair_order=_as_float(reduction["minimum_finest_pair_order"]),
        fixed_outer_rows=reduction["fixed_outer_rows"],
        coarsest_guard_maximum=_as_float(reduction["coarsest_guard_maximum"]),
        finest_guard_maximum=_as_float(reduction["finest_guard_maximum"]),
    )
    independent_endpoint = _serial(diagnosis.endpoint_assessment)
    terminal_event = event_records[-1]
    if (
        independent_endpoint != result.get("endpoint_assessment")
        or independent_endpoint != terminal_event.get("assessment")
        or metadata.get("terminal_result") != result
        or result.get("classification")
        != "completed_target_and_complete_constraint_pass"
        or result.get("RSP2_target_order_cleared") is not True
        or result.get("RSP2_complete_constraint_admission_passed") is not True
        or result.get("PROTO13_frozen") is not False
        or result.get("SGBL_outcome_read") is not False
        or result.get("FGCQR_outcome_read") is not False
        or result.get("mechanism_question_answered") is not False
        or result.get("stop") is not None
    ):
        raise ValueError("independent PREF14 endpoint differs from the runner record")

    target = independent_endpoint
    expected_vectors = {
        "target_raw_owned_norms": [
            _as_float(value) for value in reduction["target_raw_owned_norms"]
        ],
        "target_accumulated_roundoff_enclosures": [
            _as_float(value)
            for value in reduction["target_accumulated_roundoff_enclosures"]
        ],
        "target_effective_norms_for_order_only": [
            _as_float(value)
            for value in reduction["target_effective_norms_for_order_only"]
        ],
    }
    if any(target[key] != expected for key, expected in expected_vectors.items()):
        raise ValueError("PREF14 target vector differs from the frozen result")
    if (
        list(diagnosis.target_adjacent_pair_orders)
        != [_as_float(value) for value in reduction["target_adjacent_pair_orders"]]
        or list(diagnosis.target_adjacent_pair_margins)
        != [_as_float(value) for value in reduction["target_adjacent_pair_margins"]]
        or diagnosis.minimum_finite_finest_pair_order
        != _as_float(reduction["minimum_finite_finest_pair_order"])
        or diagnosis.target_previous_pair_failed is not True
        or diagnosis.target_finest_pair_passed is not True
        or diagnosis.target_preasymptotic_on_tested_ladder is not True
        or diagnosis.complete_constraint_admission_passed is not True
    ):
        raise ValueError("PREF14 adjacent-pair or complete-admission result differs")

    mutations = _mutation_controls(
        config, states, grids, stage_counts, state_hashes
    )
    if not all(mutations.values()):
        raise ValueError("PREF14 mutation control failed")
    implementation = {_rel(path): _sha(path) for path in IMPLEMENTATION}
    nonclaims = {
        key: False for key, value in EXPECTED_CLAIMS.items() if value is False
    }
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "post_run_checkpoint_recomputed_target_and_complete_constraint_pass",
        "generated_by": _rel(Path(__file__)),
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): _sha(config_path)},
        "predecessor_sha256": {
            _rel(authorization_path): _sha(authorization_path),
            config["predecessor_checkpoint"]: predecessor_hash,
        },
        "raw_evidence_sha256": {
            config["study_manifest"]: observed_raw_hashes["manifest"],
            config["study_event_log"]: observed_raw_hashes["event_log"],
            config["study_checkpoint"]: observed_raw_hashes["checkpoint"],
            config["study_result"]: observed_raw_hashes["study_result"],
        },
        "implementation_sha256": implementation,
        "scope_bindings": dict(config["scope"]),
        "gate_status": dict(EXPECTED_CLAIMS),
        "artifact_payload": {
            "immutable_study": {
                **dict(config["immutable_study"]),
                "authorization_commit_is_ancestor_of_HEAD": True,
                "authorization_blob_matches_commit": True,
                "launch_manifest_bound_before_outcome": True,
            },
            "terminal_member": dict(config["terminal_member"]),
            "independent_constraint_reduction": {
                "method": "SSPRK3",
                "point_counts": list(PREF14_POINT_COUNTS),
                "accepted_stage_counts": list(stage_counts),
                "state_sha256": list(diagnosis.state_sha256),
                "coordinate_time": PREF14_TARGET_TIME,
                "target_component": "radial_momentum",
                "minimum_finest_pair_order": 1.5,
                "target_raw_owned_norms": target["target_raw_owned_norms"],
                "target_accumulated_roundoff_enclosures": target[
                    "target_accumulated_roundoff_enclosures"
                ],
                "target_effective_norms_for_order_only": target[
                    "target_effective_norms_for_order_only"
                ],
                "target_adjacent_pair_orders": list(
                    diagnosis.target_adjacent_pair_orders
                ),
                "target_adjacent_pair_margins": list(
                    diagnosis.target_adjacent_pair_margins
                ),
                "target_previous_pair_failed": diagnosis.target_previous_pair_failed,
                "target_finest_pair_passed": diagnosis.target_finest_pair_passed,
                "target_preasymptotic_on_tested_ladder": (
                    diagnosis.target_preasymptotic_on_tested_ladder
                ),
                "component_finest_pair_orders": target["constraint_admission"][
                    "component_finest_pair_orders"
                ],
                "component_status": target["constraint_admission"][
                    "component_status"
                ],
                "minimum_finite_finest_pair_order": (
                    diagnosis.minimum_finite_finest_pair_order
                ),
                "complete_constraint_admission_passed": (
                    diagnosis.complete_constraint_admission_passed
                ),
                "classification": target["classification"],
            },
            "exact_crosschecks": {
                "checkpoint_event_log_matches_external": True,
                "checkpoint_result_matches_external": True,
                "terminal_event_matches_external_result": True,
                "independent_endpoint_matches_terminal_event": True,
                "independent_endpoint_matches_external_result": True,
                "raw_norms_preserved_before_roundoff_classification": True,
            },
            "runtime_provenance": {
                "event_log_line_count": 25,
                "accepted_boundary_count": 23,
                "source_retry_count": 0,
                "CFL_retry_count": 0,
                "resume_recovery": result["resume_recovery"],
                "stop": result["stop"],
                "elapsed_wall_seconds": result["elapsed_wall_seconds"],
            },
            "mutation_controls": mutations,
            "successor_boundary": dict(config["successor_boundary"]),
            "epistemic_boundary": {
                "CAL9_shortfall_was_preasymptotic_on_this_tested_ladder": True,
                "RSP2_is_not_a_fresh_GR0_calibration": True,
                "RSP2_does_not_select_an_eligible_GR0_case": True,
                "PROTO13_requires_a_new_prospective_freeze": True,
                "candidate_or_mechanism_outcome_read": False,
                "retained_EFT_or_physical_transition_promoted": False,
            },
        },
        "nonclaims": nonclaims,
    }


def _validate_tracked_record(
    record: Mapping[str, Any], config_path: Path = DEFAULT_CONFIG
) -> None:
    config = load_config(config_path)
    if (
        record.get("schema_version") != 1
        or record.get("artifact_id") != ARTIFACT_ID
        or record.get("project_version") != PROJECT_VERSION
        or record.get("classification")
        != "post_run_checkpoint_recomputed_target_and_complete_constraint_pass"
        or record.get("generated_by") != _rel(Path(__file__))
        or record.get("scope_bindings") != config["scope"]
        or record.get("gate_status") != EXPECTED_CLAIMS
        or record.get("nonclaims")
        != {key: False for key, value in EXPECTED_CLAIMS.items() if value is False}
        or record.get("source_config_sha256") != {_rel(config_path): _sha(config_path)}
        or record.get("derivation_document") != _rel(OWNER_DOCUMENT)
        or record.get("derivation_document_sha256") != _sha(OWNER_DOCUMENT)
    ):
        raise ValueError("tracked PREF14 record boundary differs")
    for ledger_name in (
        "predecessor_sha256",
        "implementation_sha256",
    ):
        ledger = record.get(ledger_name)
        if not isinstance(ledger, dict) or not ledger:
            raise ValueError(f"tracked PREF14 {ledger_name} is absent")
        for relative, expected in ledger.items():
            path = REPOSITORY / relative
            if not path.is_file() or _sha(path) != expected:
                raise ValueError(f"tracked PREF14 file drifted: {relative}")
    raw = record.get("raw_evidence_sha256")
    if not isinstance(raw, dict) or set(raw) != {
        config["study_manifest"],
        config["study_event_log"],
        config["study_checkpoint"],
        config["study_result"],
    }:
        raise ValueError("tracked PREF14 raw ledger differs")
    expected_raw = {
        config["study_manifest"]: config["immutable_study"]["manifest_sha256"],
        config["study_event_log"]: config["immutable_study"]["event_log_sha256"],
        config["study_checkpoint"]: config["immutable_study"]["checkpoint_sha256"],
        config["study_result"]: config["immutable_study"]["study_result_sha256"],
    }
    if raw != expected_raw:
        raise ValueError("tracked PREF14 raw hashes differ")
    payload = record.get("artifact_payload", {})
    reduction = payload.get("independent_constraint_reduction", {})
    successor = payload.get("successor_boundary", {})
    if (
        reduction.get("target_adjacent_pair_orders")
        != [1.4990947384363291, 1.9817436463819333]
        or reduction.get("target_adjacent_pair_margins")
        != [-0.0009052615636708783, 0.48174364638193334]
        or reduction.get("complete_constraint_admission_passed") is not True
        or reduction.get("classification")
        != "completed_target_and_complete_constraint_pass"
        or successor != config["successor_boundary"]
        or not all(payload.get("mutation_controls", {}).values())
    ):
        raise ValueError("tracked PREF14 numerical or successor evidence differs")


def _raw_bundle_is_complete_or_absent(paths: Sequence[Path]) -> bool:
    """Return whether all raw files exist; reject a misleading partial bundle."""

    records = tuple(Path(path) for path in paths)
    if not records:
        raise ValueError("RSP2 raw bundle path set is empty")
    present = tuple(path.is_file() for path in records)
    if any(present) and not all(present):
        raise ValueError("partial RSP2 raw bundle is not admissible")
    return all(present)


def verify_canonical(
    output_path: Path = DEFAULT_OUTPUT,
    config_path: Path = DEFAULT_CONFIG,
) -> dict[str, Any]:
    stored = _load_json(output_path, ARTIFACT_ID)
    _validate_tracked_record(stored, config_path)
    config = load_config(config_path)
    raw_paths = tuple(
        REPOSITORY / config[key]
        for key in (
            "study_manifest",
            "study_event_log",
            "study_checkpoint",
            "study_result",
        )
    )
    if _raw_bundle_is_complete_or_absent(raw_paths):
        reproduced = reproduce(config_path)
        if reproduced != stored:
            raise ValueError("canonical PREF14 record differs from raw replay")
    return stored


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
    _atomic_write(output_path, _canonical(record))
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
