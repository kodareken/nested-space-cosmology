#!/usr/bin/env python3
"""Reproduce SRC4/VEC1's tensor-contraction source certificate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess
import tempfile
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.gr0_direct_source import (  # noqa: E402
    _metric_jet_from_adm,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.src2_exact_oracle import (  # noqa: E402
    exact_affine_system,
    exact_gaussian_solve,
    exact_ref1_residual,
)
from recursive_horizons.fgc.evolution.src3_reference_balanced_source import (  # noqa: E402
    _inverse_and_derivative,
    _reference_covariant_connection_difference,
    _reference_covariant_ricci,
    diagnose_gr0_reference_balanced_accelerations,
    gr0_reference_balanced_ref1_residual_batch,
)
from recursive_horizons.fgc.evolution.src4_vectorized_reference_source import (  # noqa: E402
    _vectorized_reference_covariant_connection_difference,
    _vectorized_reference_covariant_ricci,
    diagnose_gr0_vectorized_reference_accelerations,
    gr0_vectorized_reference_ref1_residual_batch,
)


ARTIFACT_ID = "FGC-1-SRC4-VEC1"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-src4-vec1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-src4-vec1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-src4-vec1.md"
SOURCE_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/src4_vectorized_reference_source.py"
)
IMPLEMENTATION = (
    Path(__file__).resolve(),
    SOURCE_MODULE,
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/src3_reference_balanced_source.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/src2_exact_oracle.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/gr0_direct_source.py",
)
FIELD_NAMES = ("u", "p", "q", "p_r", "q_r")
EXPECTED_CLAIMS = {
    "SRC4_tensor_contraction_evaluator_derived": True,
    "SRC4_exact_reference_and_exact_oracle_controls_passed": True,
    "SRC4_independent_SRC3_differential_controls_passed": True,
    "SRC4_CAP1_full_grid_strict_source_gate_passed": True,
    "SRC4_runtime_source_instrument_authorized_for_a_separately_frozen_GR0_protocol": True,
    "wall_clock_speedup_is_a_scientific_result": False,
    "PROTO12_frozen": False,
    "fresh_GR0_calibration_completed": False,
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
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _optional_raw_exists(path: Path) -> bool:
    """Is an optional untracked raw fixture available for extra validation?"""

    return path.is_file()


def _bytes_sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _canonical(value: Mapping[str, Any]) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        indent=2,
        ensure_ascii=True,
        allow_nan=False,
    ) + "\n"


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
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


def _strict_keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(f"{name} keys differ")


def _fraction_float(value: object) -> float:
    if not isinstance(value, str):
        raise TypeError("frozen numerical values must be strings")
    return float(Fraction(value))


def _hex_float(value: object) -> float:
    if not isinstance(value, str):
        raise TypeError("binary64 values must be hex strings")
    answer = float.fromhex(value)
    if not np.isfinite(answer) or answer.hex() != value:
        raise ValueError("binary64 value is not canonical")
    return answer


def _hex_row(value: object, *, name: str) -> np.ndarray:
    if not isinstance(value, list) or len(value) != 6:
        raise ValueError(f"{name} must contain six values")
    return np.asarray([_hex_float(item) for item in value], dtype=np.float64)


def _ordered_float_bits(value: float) -> int:
    signed = struct.unpack(">q", struct.pack(">d", value))[0]
    return 0x8000000000000000 - signed if signed < 0 else signed


def _ulp_distance(left: float, right: float) -> int:
    if left == right:
        return 0
    return abs(_ordered_float_bits(left) - _ordered_float_bits(right))


def _git(*args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *args],
        cwd=REPOSITORY,
        check=check,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    _strict_keys(
        "SRC4 config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "scope",
            "immutable_lineage",
            "evaluator",
            "proof_contract",
            "expected",
            "successor_boundary",
            "claims",
        },
    )
    if (
        config["schema_version"] != 1
        or config["artifact_id"] != ARTIFACT_ID
        or config["metric_signature"] != "-+++"
        or config["riemann_convention"] != "plus_partial_mu_gamma_nu"
        or config["scope"]
        != {
            "role": "prospectively_frozen_tensor_contraction_backend_for_unchanged_SRC3_GR0_REF1_evaluator",
            "branch": "GR-0",
            "existing_SRC3_and_RSP1_outputs_used_as_design_evidence": True,
            "fresh_GR0_trajectory_read": False,
            "SGBL_trajectory_read": False,
            "FGCQR_trajectory_read": False,
            "collapse_or_trapped_outcome_classified": False,
            "mechanism_question_answered": False,
        }
        or any(
            value is not True
            for value in config["proof_contract"].values()
            if isinstance(value, bool)
        )
        or config["claims"] != EXPECTED_CLAIMS
        or config["successor_boundary"].get("next_gate") != "FGC-1-PRO12-FRZ1"
        or any(
            value is not True
            for key, value in config["successor_boundary"].items()
            if key != "next_gate"
        )
    ):
        raise ValueError("SRC4 top-level contract differs")
    expected = config["expected"]
    if set(expected) != {"PREF11", "DYADIC_A", "DYADIC_B", "CAP1"}:
        raise ValueError("SRC4 expected-control tables differ")
    return config


def _validate_lineage(config: Mapping[str, Any]) -> None:
    lineage = config["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("SRC4 checkpoint is not an ancestor of HEAD")
    bindings = (
        ("SRC3_result", "SRC3_result_sha256"),
        ("RSP1_result", "RSP1_result_sha256"),
        ("SRC3_source", "SRC3_source_sha256"),
        ("PREF11_point_fixture", "PREF11_point_fixture_sha256"),
        ("SRC3_control_fixture", "SRC3_control_fixture_sha256"),
    )
    for path_key, hash_key in bindings:
        path = REPOSITORY / lineage[path_key]
        if not path.is_file() or _sha(path) != lineage[hash_key]:
            raise ValueError(f"SRC4 immutable lineage differs: {lineage[path_key]}")
        committed = _git("show", f"{commit}:{lineage[path_key]}").stdout
        if _bytes_sha(committed) != lineage[hash_key]:
            raise ValueError(f"SRC4 checkpoint blob differs: {lineage[path_key]}")
    src3 = json.loads((REPOSITORY / lineage["SRC3_result"]).read_text())
    rsp1 = json.loads((REPOSITORY / lineage["RSP1_result"]).read_text())
    if (
        src3.get("artifact_id") != "FGC-1-SRC3"
        or src3.get("gate_status", {}).get(
            "SRC3_runtime_source_instrument_authorized_for_a_separately_frozen_GR0_protocol"
        )
        is not True
        or rsp1.get("artifact_id") != "FGC-1-RSP1-PREF12"
        or rsp1.get("gate_status", {}).get(
            "amplitude_three_spectral_veto_cleared_for_successor_design"
        )
        is not True
        or rsp1.get("gate_status", {}).get("PROTO12_frozen") is not False
    ):
        raise ValueError("SRC4 predecessor scope differs")


def _load_fixtures(
    config: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    lineage = config["immutable_lineage"]
    point = json.loads(
        (REPOSITORY / lineage["PREF11_point_fixture"]).read_text(encoding="utf-8")
    )
    controls = json.loads(
        (REPOSITORY / lineage["SRC3_control_fixture"]).read_text(encoding="utf-8")
    )
    if not isinstance(point, dict) or not isinstance(controls, dict):
        raise ValueError("SRC4 fixtures must be mappings")
    return point, controls


def _lower_jet(value: Mapping[str, Any]) -> tuple[np.ndarray, ...]:
    return tuple(_hex_row(value[name], name=name) for name in FIELD_NAMES)


def _expected_record(
    expected: Mapping[str, Any],
    solved: Any,
) -> None:
    if (
        solved.residual_infinity != _fraction_float(expected["residual_infinity"])
        or solved.kinetic_condition_infinity_maximum
        != _fraction_float(expected["kinetic_condition_infinity_maximum"])
        or array_content_sha256(solved.accelerations)
        != expected["acceleration_sha256"]
        or array_content_sha256(solved.residuals) != expected["residual_sha256"]
    ):
        raise ValueError("SRC4 frozen root record differs")


def _reference_record(controls: Mapping[str, Any]) -> dict[str, Any]:
    radii = np.asarray(
        [_hex_float(item) for item in controls["exact_reference_radii_binary64_hex"]],
        dtype=np.float64,
    )
    points = radii.size
    u = np.zeros((points, 6), dtype=np.float64)
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    p_r = np.zeros_like(u)
    q_r = np.zeros_like(u)
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = radii
    q[:, 3] = 1.0
    residual = gr0_vectorized_reference_ref1_residual_batch(
        u,
        p,
        q,
        np.zeros((1, points, 6), dtype=np.float64),
        p_r,
        q_r,
        radii,
    )
    arrays = (
        residual.full_residual,
        residual.unredefined_metric_residual,
        residual.scalar_residual,
        residual.gauge_constraint,
        residual.hamiltonian_constraint,
        residual.momentum_constraint,
        residual.ricci_scalar,
        residual.ricci_squared,
    )
    solved = diagnose_gr0_vectorized_reference_accelerations(
        u, p, q, p_r, q_r, radii
    )
    if (
        not all(np.array_equal(value, np.zeros_like(value)) for value in arrays)
        or not solved.raw_gate_passed
        or solved.residual_infinity != 0.0
        or not np.array_equal(
            solved.accelerations,
            np.zeros_like(solved.accelerations),
        )
    ):
        raise ValueError("SRC4 exact reference is not bitwise stationary")
    return {
        "radius_count": points,
        "complete_residual_bitwise_zero": True,
        "acceleration_bitwise_zero": True,
    }


def _point_record(
    config: Mapping[str, Any],
    point: Mapping[str, Any],
) -> dict[str, Any]:
    lower = _lower_jet(point["lower_jet"])
    radius = _hex_float(point["coordinate_radius"])
    arrays = tuple(value[None, :] for value in lower)
    solved = diagnose_gr0_vectorized_reference_accelerations(
        *arrays,
        np.asarray((radius,), dtype=np.float64),
    )
    src3 = diagnose_gr0_reference_balanced_accelerations(
        *arrays,
        np.asarray((radius,), dtype=np.float64),
    )
    _expected_record(config["expected"]["PREF11"], solved)
    exact_constant, exact_matrix = exact_affine_system(*lower, radius)
    exact_root = exact_gaussian_solve(
        exact_matrix,
        tuple(-value for value in exact_constant),
    )
    rounded = np.asarray([float(value) for value in exact_root], dtype=np.float64)
    exact_residual = exact_ref1_residual(*lower, radius, solved.accelerations[0])
    exact_infinity = float(
        max((abs(value) for value in exact_residual), default=Fraction(0))
    )
    ulp_maximum = max(
        _ulp_distance(left, right)
        for left, right in zip(solved.accelerations[0], rounded, strict=True)
    )
    delta = float(np.max(np.abs(solved.accelerations - src3.accelerations)))
    proof = config["proof_contract"]
    if (
        not solved.raw_gate_passed
        or exact_infinity
        > _fraction_float(proof["PREF11_exact_dyadic_root_residual_maximum"])
        or ulp_maximum > proof["PREF11_maximum_ULP_distance_from_rounded_exact_root"]
        or delta
        > _fraction_float(proof["compact_root_acceleration_delta_from_SRC3_maximum"])
    ):
        raise ValueError("SRC4 PREF11 point control failed")
    return {
        "radius_binary64_hex": radius.hex(),
        "residual_infinity": solved.residual_infinity,
        "exact_dyadic_residual_infinity": exact_infinity,
        "maximum_ULP_distance_from_rounded_exact_root": ulp_maximum,
        "SRC3_acceleration_delta_infinity": delta,
        "acceleration_sha256": array_content_sha256(solved.accelerations),
        "residual_sha256": array_content_sha256(solved.residuals),
        "strict_raw_gate_passed": True,
    }


def _control_record(
    config: Mapping[str, Any],
    control: Mapping[str, Any],
) -> dict[str, Any]:
    lower = _lower_jet(control["lower_jet"])
    radius = _hex_float(control["radius_binary64_hex"])
    arrays = tuple(value[None, :] for value in lower)
    seeds = np.zeros((7, 1, 6), dtype=np.float64)
    for field in range(6):
        seeds[field + 1, 0, field] = 1.0
    vector = gr0_vectorized_reference_ref1_residual_batch(
        *arrays[:3],
        seeds,
        *arrays[3:],
        np.asarray((radius,), dtype=np.float64),
    ).full_residual[:, 0]
    src3_seed = gr0_reference_balanced_ref1_residual_batch(
        *arrays[:3],
        seeds,
        *arrays[3:],
        np.asarray((radius,), dtype=np.float64),
    ).full_residual[:, 0]
    exact = np.asarray(
        [
            [
                float(value)
                for value in exact_ref1_residual(*lower, radius, acceleration)
            ]
            for acceleration in seeds[:, 0]
        ],
        dtype=np.float64,
    )
    scale = max(1.0, float(np.max(np.abs(exact))))
    exact_error = float(np.max(np.abs(vector - exact))) / scale
    src3_error = float(np.max(np.abs(vector - src3_seed))) / scale
    solved = diagnose_gr0_vectorized_reference_accelerations(
        *arrays,
        np.asarray((radius,), dtype=np.float64),
    )
    src3_root = diagnose_gr0_reference_balanced_accelerations(
        *arrays,
        np.asarray((radius,), dtype=np.float64),
    )
    table = "DYADIC_A" if control["control_id"] == "DYADIC-A" else "DYADIC_B"
    expected = config["expected"][table]
    if expected["control_id"] != control["control_id"]:
        raise ValueError("SRC4 control identity differs")
    _expected_record(expected, solved)
    limit = _fraction_float(
        config["proof_contract"][
            "independent_control_normalized_residual_error_maximum"
        ]
    )
    root_delta = float(np.max(np.abs(solved.accelerations - src3_root.accelerations)))
    if (
        exact_error > limit
        or src3_error > limit
        or root_delta
        > _fraction_float(
            config["proof_contract"][
                "compact_root_acceleration_delta_from_SRC3_maximum"
            ]
        )
        or not solved.raw_gate_passed
    ):
        raise ValueError("SRC4 independent control failed")
    return {
        "control_id": control["control_id"],
        "radius_binary64_hex": radius.hex(),
        "exact_seed_normalized_error": exact_error,
        "SRC3_seed_normalized_delta": src3_error,
        "SRC3_root_acceleration_delta_infinity": root_delta,
        "residual_infinity": solved.residual_infinity,
        "kinetic_condition_infinity_maximum": (
            solved.kinetic_condition_infinity_maximum
        ),
        "acceleration_sha256": array_content_sha256(solved.accelerations),
        "residual_sha256": array_content_sha256(solved.residuals),
        "strict_raw_gate_passed": True,
    }


def _contraction_record(
    config: Mapping[str, Any],
    controls: Mapping[str, Any],
) -> dict[str, Any]:
    records = controls["independent_nontrivial_controls"]
    lower = {
        name: np.stack(
            [_hex_row(item["lower_jet"][name], name=name) for item in records],
            axis=0,
        )
        for name in FIELD_NAMES
    }
    radii = np.asarray(
        [_hex_float(item["radius_binary64_hex"]) for item in records],
        dtype=np.float64,
    )
    seeds = np.zeros((7, len(records), 6), dtype=np.float64)
    for field in range(6):
        seeds[field + 1, :, field] = 1.0
    metric, first, second = _metric_jet_from_adm(
        lower["u"],
        lower["p"],
        lower["q"],
        seeds,
        lower["p_r"],
        lower["q_r"],
    )
    inverse, inverse_derivative = _inverse_and_derivative(metric, first)
    src3 = _reference_covariant_connection_difference(
        metric,
        first,
        second,
        inverse,
        inverse_derivative,
        radii,
    )
    src4 = _vectorized_reference_covariant_connection_difference(
        metric,
        first,
        second,
        inverse,
        inverse_derivative,
        radii,
    )
    connection_delta = float(np.max(np.abs(src3[2] - src4[2])))
    derivative_delta = float(np.max(np.abs(src3[3] - src4[3])))
    src3_ricci = _reference_covariant_ricci(src3[0], src3[2], src3[3])
    src4_ricci = _vectorized_reference_covariant_ricci(
        src4[0], src4[2], src4[3]
    )
    ricci_scale = max(1.0, float(np.max(np.abs(src3_ricci))))
    ricci_delta = float(np.max(np.abs(src3_ricci - src4_ricci))) / ricci_scale
    proof = config["proof_contract"]
    if (
        not np.array_equal(src3[0], src4[0])
        or not np.array_equal(src3[1], src4[1])
        or connection_delta
        > _fraction_float(
            proof["vectorized_connection_difference_maximum_absolute_delta_from_SRC3"]
        )
        or derivative_delta
        > _fraction_float(
            proof[
                "vectorized_connection_derivative_maximum_absolute_delta_from_SRC3"
            ]
        )
        or ricci_delta
        > _fraction_float(proof["vectorized_Ricci_normalized_delta_maximum"])
    ):
        raise ValueError("SRC4 tensor-contraction differential control failed")
    return {
        "reference_connection_bitwise_equal": True,
        "reference_connection_derivative_bitwise_equal": True,
        "connection_difference_delta_infinity": connection_delta,
        "connection_difference_derivative_delta_infinity": derivative_delta,
        "Ricci_normalized_delta_infinity": ricci_delta,
        "point_and_seed_axes_contracted": False,
    }


def _attacks(point: Mapping[str, Any]) -> dict[str, Any]:
    lower = list(_lower_jet(point["lower_jet"]))
    radius = _hex_float(point["coordinate_radius"])
    arrays = tuple(value[None, :] for value in lower)
    baseline = diagnose_gr0_vectorized_reference_accelerations(
        *arrays,
        np.asarray((radius,), dtype=np.float64),
    )
    mutation = lower[0].copy()
    mutation[0] = np.nextafter(mutation[0], np.inf)
    changed = diagnose_gr0_vectorized_reference_accelerations(
        mutation[None, :],
        *(value[None, :] for value in lower[1:]),
        np.asarray((radius,), dtype=np.float64),
    )
    failures: dict[str, bool] = {}
    try:
        diagnose_gr0_vectorized_reference_accelerations(
            *arrays,
            np.asarray((0.0,), dtype=np.float64),
        )
    except ValueError:
        failures["nonpositive_radius"] = True
    try:
        diagnose_gr0_vectorized_reference_accelerations(
            *arrays,
            np.asarray((radius,), dtype=np.float64),
            condition_number_maximum=1.0,
        )
    except ValueError:
        failures["condition_limit"] = True
    try:
        gr0_vectorized_reference_ref1_residual_batch(
            *arrays[:3],
            np.zeros((1, 1, 6), dtype=np.float64),
            *arrays[3:],
            np.asarray((radius,), dtype=np.float64),
            tilde_normal_factor=4.0,
            hat_normal_factor=4.0,
        )
    except ValueError:
        failures["invalid_auxiliary_factors"] = True
    visible = not np.array_equal(baseline.accelerations, changed.accelerations)
    if not visible or failures != {
        "nonpositive_radius": True,
        "condition_limit": True,
        "invalid_auxiliary_factors": True,
    }:
        raise ValueError("SRC4 mutation or typed-failure controls failed")
    return {
        "one_bit_alpha_mutation_visible": visible,
        "baseline_acceleration_sha256": array_content_sha256(
            baseline.accelerations
        ),
        "mutated_acceleration_sha256": array_content_sha256(changed.accelerations),
        "typed_failure_controls": failures,
    }


def _capture_record(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    expected = config["expected"]["CAP1"]
    raw_path = REPOSITORY / lineage["CAP1_raw_fixture"]
    expected_acceleration_delta = _fraction_float(
        expected["SRC4_acceleration_delta_from_SRC3"]
    )
    expected_condition_delta = _fraction_float(
        expected["SRC4_condition_delta_from_SRC3"]
    )
    if _optional_raw_exists(raw_path):
        if _sha(raw_path) != lineage["CAP1_raw_fixture_sha256"]:
            raise ValueError("SRC4 CAP1 raw fixture hash differs")
        with np.load(raw_path, allow_pickle=False) as raw:
            required = {
                "metadata_utf8",
                "u",
                "p",
                "q",
                "p_r",
                "q_r",
                "radii",
                "forward_constant",
                "forward_jacobian",
                "baseline_acceleration",
                "baseline_complete_residual",
                "best_acceleration",
                "best_complete_residual",
            }
            if set(raw.files) != required:
                raise ValueError("SRC4 CAP1 raw fixture inventory differs")
            arguments = tuple(raw[name] for name in (*FIELD_NAMES, "radii"))
            src3 = diagnose_gr0_reference_balanced_accelerations(*arguments)
            src4 = diagnose_gr0_vectorized_reference_accelerations(*arguments)
        last = src4.iterations[-1]
        _expected_record(expected, src4)
        acceleration_delta = float(
            np.max(np.abs(src4.accelerations - src3.accelerations))
        )
        condition_delta = abs(
            src4.kinetic_condition_infinity_maximum
            - src3.kinetic_condition_infinity_maximum
        )
        if (
            src4.accelerations.shape[0] != expected["point_count"]
            or not src3.raw_gate_passed
            or not src4.raw_gate_passed
            or src3.residual_infinity
            != _fraction_float(expected["SRC3_residual_infinity"])
            or src3.kinetic_condition_infinity_maximum
            != _fraction_float(expected["SRC3_kinetic_condition_infinity_maximum"])
            or last.maximum_residual_point_index
            != expected["maximum_residual_point_index"]
            or last.maximum_residual_row_index
            != expected["maximum_residual_row_index"]
            or last.maximum_residual_radius.hex()
            != expected["maximum_residual_radius_binary64_hex"]
            or acceleration_delta
            > _fraction_float(
                config["proof_contract"][
                    "CAP1_acceleration_delta_from_SRC3_maximum"
                ]
            )
            or condition_delta
            > _fraction_float(
                config["proof_contract"]["CAP1_condition_delta_from_SRC3_maximum"]
            )
            or acceleration_delta != expected_acceleration_delta
            or condition_delta != expected_condition_delta
        ):
            raise ValueError("SRC4 CAP1 full-grid control differs")
    return {
        "point_count": expected["point_count"],
        "SRC4_residual_infinity": _fraction_float(expected["residual_infinity"]),
        "SRC4_kinetic_condition_infinity_maximum": _fraction_float(
            expected["kinetic_condition_infinity_maximum"]
        ),
        "SRC4_acceleration_sha256": expected["acceleration_sha256"],
        "SRC4_residual_sha256": expected["residual_sha256"],
        "SRC3_residual_infinity": _fraction_float(
            expected["SRC3_residual_infinity"]
        ),
        "SRC3_kinetic_condition_infinity_maximum": _fraction_float(
            expected["SRC3_kinetic_condition_infinity_maximum"]
        ),
        "both_strict_raw_gates_passed_in_frozen_CAP1_evidence": True,
        "SRC4_acceleration_delta_from_SRC3": expected_acceleration_delta,
        "SRC4_condition_delta_from_SRC3": expected_condition_delta,
        "raw_fixture_validation_contract": (
            "recomputed_when_present_otherwise_bound_by_frozen_hash"
        ),
        "raw_fixture_optional_for_clean_clone_reproduction": True,
        "hash_mismatched_raw_fixture_fails_closed": True,
    }


def reproduce(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    _validate_lineage(config)
    point, controls = _load_fixtures(config)
    reference = _reference_record(controls)
    point_result = _point_record(config, point)
    independent = [
        _control_record(config, item)
        for item in controls["independent_nontrivial_controls"]
    ]
    contractions = _contraction_record(config, controls)
    attacks = _attacks(point)
    capture = _capture_record(config)
    lineage = config["immutable_lineage"]
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": config["project_version"],
        "classification": "prospectively_frozen_tensor_contraction_backend_for_unchanged_SRC3_GR0_REF1",
        "generated_by": _rel(Path(__file__)),
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): _sha(config_path)},
        "predecessor_sha256": {
            lineage["SRC3_result"]: lineage["SRC3_result_sha256"],
            lineage["RSP1_result"]: lineage["RSP1_result_sha256"],
            lineage["SRC3_source"]: lineage["SRC3_source_sha256"],
            lineage["PREF11_point_fixture"]: lineage[
                "PREF11_point_fixture_sha256"
            ],
            lineage["SRC3_control_fixture"]: lineage[
                "SRC3_control_fixture_sha256"
            ],
        },
        "implementation_sha256": {
            _rel(path): _sha(path) for path in IMPLEMENTATION
        },
        "scope_bindings": dict(config["scope"]),
        "gate_status": dict(EXPECTED_CLAIMS),
        "nonclaims": {
            key: value for key, value in EXPECTED_CLAIMS.items() if value is False
        },
        "artifact_payload": {
            "immutable_checkpoint_commit": lineage["checkpoint_commit"],
            "evaluator_contract": dict(config["evaluator"]),
            "exact_reference_control": reference,
            "PREF11_exact_dyadic_control": point_result,
            "independent_controls": independent,
            "tensor_contraction_differential": contractions,
            "CAP1_full_grid_control": capture,
            "adversarial_controls": attacks,
            "successor_boundary": dict(config["successor_boundary"]),
            "epistemic_boundary": {
                "observed_results": [
                    "the tensor contractions preserve exact spherical Minkowski bitwise",
                    "PREF11 and two independent dyadic controls agree with the exact oracle and independent SRC3 backend within frozen binary64 bounds",
                    "the captured 8192-point CAP1 state preserves the unchanged strict source and kinetic gates",
                    "a one-bit mutation remains visible and typed invalid domains fail closed",
                ],
                "disallowed_implications": [
                    "wall-clock timing is not a scientific result or source of authorization",
                    "SRC4 does not change an equation, threshold, root branch, field, coupling, or reference geometry",
                    "RSP1 generic all-field raw spectral admission remains false",
                    "PROTO12 is not frozen and no fresh calibration or candidate trajectory is authorized",
                    "no collapse, defocusing, retained-EFT, singularity-resolution, child-domain, dark-sector, or varying-c claim follows",
                ],
            },
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    record = reproduce(args.config.resolve())
    payload = _canonical(record).encode("utf-8")
    if args.check:
        if not args.output.is_file() or args.output.read_bytes() != payload:
            raise SystemExit("stored SRC4/VEC1 result differs")
        print(json.dumps({"verified": True, "artifact_id": ARTIFACT_ID}))
        return
    _atomic_write(args.output.resolve(), payload)
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
