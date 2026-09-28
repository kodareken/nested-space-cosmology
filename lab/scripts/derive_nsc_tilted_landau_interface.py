#!/usr/bin/env python3
"""Build/check the tilted transmitting weighted-isometry family record."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.nsc_mode_resolved_cauchy_state import (  # noqa: E402
    ModeResolvedCauchyState,
    array_digest,
    deterministic_npz_bytes,
)
from recursive_horizons.nsc_tilted_landau_interface import (  # noqa: E402
    TransmittingTiltedLandauInterface,
    frequency_fourier_unitary,
    landau_spin_factors,
)


OUTPUT = ROOT / "results/development/nsc-tilted-landau-interface.json"
ARTIFACT_DIRECTORY = ROOT / "results/development/artifacts"
INPUTS = (
    "results/development/nsc-mode-resolved-cauchy-state.json",
    "results/development/nsc-constraint-complete-neck.json",
    "results/nsc-3-boundary-response.json",
    "results/nsc-8-chiral-boundary.json",
    "results/development/curvature-eft.json",
)
SOURCES = (
    "src/recursive_horizons/nsc_tilted_landau_interface.py",
    "scripts/derive_nsc_tilted_landau_interface.py",
    "docs/nsc-tilted-landau-interface.md",
)
TOLERANCE = 3e-11


def hashes(paths):
    return {
        path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
        for path in paths
    }


def compare(expected, actual, path="$", atol=3e-12, rtol=3e-12):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys():
            raise AssertionError(f"keys differ at {path}")
        for key in expected:
            compare(expected[key], actual[key], f"{path}/{key}", atol, rtol)
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise AssertionError(f"list differs at {path}")
        for index, (left, right) in enumerate(zip(expected, actual)):
            compare(left, right, f"{path}/{index}", atol, rtol)
    elif isinstance(expected, float):
        if isinstance(actual, bool) or not isinstance(actual, (int, float)):
            raise AssertionError(f"numeric type differs at {path}")
        if abs(expected - actual) > atol + rtol * abs(expected):
            raise AssertionError(f"number differs at {path}: {expected} != {actual}")
    elif type(expected) is not type(actual) or expected != actual:
        raise AssertionError(f"value differs at {path}: {expected!r} != {actual!r}")


def load_records():
    return tuple(json.loads((ROOT / path).read_text()) for path in INPUTS)


def load_state(manifest):
    path = ROOT / manifest["payload"]["path"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != manifest["payload"]["sha256"]:
        raise ValueError("mode-state payload hash mismatch")
    return ModeResolvedCauchyState.load(path, manifest["channels"])


def build_witness_arrays(state, rapidity):
    interface = TransmittingTiltedLandauInterface(rapidity, TOLERANCE)
    offsets = np.asarray(state.arrays["sample_offsets"], dtype=np.int64)
    dimensions = []
    value_offsets = [0]
    mixed_values = []
    kernel_values = []
    rows = []
    identity_frobenius = 0.0
    mixing_frobenius = 0.0
    for channel_index, channel in enumerate(state.channels):
        start, end = map(int, offsets[channel_index:channel_index + 2])
        weights = state.arrays["quadrature_weight"][start:end]
        covariance = state.arrays["covariance_seed"][start:end]
        dimension = 2 * (end - start)
        identity = np.eye(dimension, dtype=complex)
        mixing = frequency_fourier_unitary(end - start)
        coordinate = interface.coordinate_kernel(weights, mixing)
        diagonal = interface.residuals(weights, covariance, identity)
        mixed = interface.residuals(weights, covariance, mixing)
        dimensions.append(dimension)
        value_offsets.append(value_offsets[-1] + dimension * dimension)
        mixed_values.append(mixing.reshape(-1))
        kernel_values.append(coordinate.reshape(-1))
        delta = mixing - identity
        identity_frobenius += float(np.linalg.norm(identity) ** 2)
        mixing_frobenius += float(np.linalg.norm(delta) ** 2)
        rows.append({
            "channel_index": channel_index,
            "family": int(channel["family"]),
            "compact_level": int(channel["compact_level"]),
            "angular_level": int(channel["angular_level"]),
            "frequency_nodes": end - start,
            "dimension": dimension,
            "identity_witness": diagonal,
            "frequency_mixing_witness": mixed,
        })
    arrays = {
        "channel_dimensions": np.asarray(dimensions, dtype="<i8"),
        "block_value_offsets": np.asarray(value_offsets, dtype="<i8"),
        "U_frequency_mixing_values": np.asarray(np.concatenate(mixed_values), dtype="<c16"),
        "K_coordinate_frequency_mixing_values": np.asarray(np.concatenate(kernel_values), dtype="<c16"),
    }
    return arrays, rows, identity_frobenius ** 0.5, mixing_frobenius ** 0.5


def verify_artifact_arrays(stored, calculated):
    if stored.keys() != calculated.keys():
        raise ValueError("interface artifact array names differ")
    for name in stored:
        if stored[name].dtype != calculated[name].dtype or stored[name].shape != calculated[name].shape:
            raise ValueError(f"interface artifact metadata differs: {name}")
        if np.max(np.abs(stored[name] - calculated[name]), initial=0.0) > 3e-13:
            raise ValueError(f"interface artifact values differ: {name}")


def payload_manifest(path, arrays):
    raw = path.read_bytes()
    return {
        "path": str(path.relative_to(ROOT)),
        "bytes": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "format": "deterministic uncompressed NPZ; complete channel-block witness matrices",
        "allow_pickle": False,
        "arrays": {
            name: {
                "dtype": value.dtype.str,
                "shape": list(value.shape),
                "sha256": array_digest(name, value),
            }
            for name, value in sorted(arrays.items())
        },
    }


def calculate(path, stored_arrays, records):
    state_record, selected, boundary, chiral, curvature = records
    state = load_state(state_record)
    rapidity = float(selected["source_selected_frame"]["rapidity"])
    calculated_arrays, rows, identity_size, witness_distance = build_witness_arrays(state, rapidity)
    verify_artifact_arrays(stored_arrays, calculated_arrays)

    identity_rows = [row["identity_witness"] for row in rows]
    mixing_rows = [row["frequency_mixing_witness"] for row in rows]
    maximum_interface = max(
        max(row["interface_variation_max_abs"] for row in identity_rows),
        max(row["interface_variation_max_abs"] for row in mixing_rows),
    )
    maximum_whitened = max(
        max(row["interface_variation_whitened_norm"] for row in identity_rows),
        max(row["interface_variation_whitened_norm"] for row in mixing_rows),
    )
    maximum_unitarity = max(
        max(row["canonical_unitarity_max_abs"] for row in identity_rows),
        max(row["canonical_unitarity_max_abs"] for row in mixing_rows),
    )
    maximum_car = max(
        max(row["CAR_spectrum_transport_max_abs"] for row in identity_rows),
        max(row["CAR_spectrum_transport_max_abs"] for row in mixing_rows),
    )
    minimum_car = min(row["target_minimum_eigenvalue"] for row in mixing_rows)
    maximum_car_eigenvalue = max(row["target_maximum_eigenvalue"] for row in mixing_rows)
    minimum_rank = min(row["kernel_rank"] for row in mixing_rows)
    total_dimension = sum(row["dimension"] for row in rows)
    sum_squared_dimensions = sum(row["dimension"] ** 2 for row in rows)
    offblock = min(row["frequency_offblock_frobenius"] for row in mixing_rows)
    pass_residuals = (
        maximum_interface < TOLERANCE
        and maximum_whitened < TOLERANCE
        and maximum_unitarity < TOLERANCE
        and maximum_car < TOLERANCE
        and minimum_car >= -TOLERANCE
        and maximum_car_eigenvalue <= 1 + TOLERANCE
        and minimum_rank == min(row["dimension"] for row in rows)
        and all(row["rank_defect"] == 0 for row in mixing_rows)
        and offblock > 1.0
        and witness_distance > 1.0
    )
    if not pass_residuals:
        raise ArithmeticError("weighted interface family witness failed")

    spin_metric, spin_half, spin_inverse_half = landau_spin_factors(rapidity)
    locked = selected["locked_inputs"]
    return {
        "schema": "NSC-TILTED-LANDAU-INTERFACE-v1",
        "owner": "TransmittingTiltedLandauInterface",
        "status": (
            "CONSTRUCTION PASS / PHYSICAL SELECTOR OPEN: the complete retained "
            "channel-block family is weighted-isometric and CAR preserving, "
            "including frequency-mixing witnesses, but the boundary/endpoint "
            "data do not select one physical kernel"
        ),
        "source_hashes": hashes(SOURCES),
        "input_hashes": hashes(INPUTS),
        "locked_inputs": {
            "A": locked["A"],
            "magnetic_flux": locked["q"],
            "Omega": locked["Omega"],
            "zeta": locked["zeta"],
            "V_full": locked["V_full"],
            "occupations_or_parameters_refitted": False,
        },
        "domain": {
            "imported_spatial_boundary_domain": chiral["algebra_and_domain"]["transmission_domain"],
            "imported_spatial_kernel_scope": chiral["algebra_and_domain"]["kernel_scope"],
            "imported_oriented_identity_transmission_residual": chiral["algebra_and_domain"]["identity_transmission_oriented_flux_residual"],
            "boundary_parent_child_orientation": boundary["exact_relations"]["orientation"],
            "curvature_boundary_allocation": "Euler, box-R and induced boundary variations remain separate in the imported first-order four-dimensional EFT owner",
            "in_scope": "33 retained symmetry channels; arbitrary dense frequency/spin unitary V_c inside each channel; quadrature-weighted tilted Cauchy norm",
            "out_of_scope": "mixing different angular/compact channels, a selected hypersurface embedding/history, new boundary shell, finite stress evaluation, and metric stepping",
        },
        "weighted_interface_family": {
            "coordinate_metric_seed": "G0=diag(quadrature_weight) tensor I2",
            "coordinate_metric_Landau": "GL=diag(quadrature_weight) tensor exp(-eta sigma2)",
            "coordinate_kernel": "K_c=GL_c^(-1/2) V_c G0_c^(1/2)",
            "canonical_Cauchy_map": "U_c=GL_c^(1/2) K_c G0_c^(-1/2)=V_c",
            "admissible_group": "direct_product_c U(2*n_c)",
            "B_L_half_is_U": False,
            "rapidity": rapidity,
            "spin_metric_B_L": [[{"real": float(x.real), "imag": float(x.imag)} for x in row] for row in spin_metric],
            "spin_metric_positive_sqrt": [[{"real": float(x.real), "imag": float(x.imag)} for x in row] for row in spin_half],
            "spin_metric_inverse_sqrt": [[{"real": float(x.real), "imag": float(x.imag)} for x in row] for row in spin_inverse_half],
        },
        "rank_nonuniqueness_certificate": {
            "channels": len(rows),
            "frequency_nodes": int(state.arrays["sample_offsets"][-1]),
            "total_canonical_dimension": total_dimension,
            "full_block_map_rank": total_dimension,
            "full_block_map_rank_defect": 0,
            "real_kernel_unknowns": 2 * sum_squared_dimensions,
            "weighted_isometry_constraint_jacobian_rank": sum_squared_dimensions,
            "admissible_family_real_dimension_and_linearized_nullity": sum_squared_dimensions,
            "identity_family_element_frobenius": identity_size,
            "identity_to_mixing_witness_frobenius": witness_distance,
            "endpoint_selects_unique_kernel": False,
            "reason": "positive G0 and GL reduce the endpoint equation to unrestricted V_c dagger V_c=I in every retained channel",
        },
        "residuals": {
            "declared_tolerance": TOLERANCE,
            "maximum_interface_variation_max_abs": maximum_interface,
            "maximum_whitened_interface_norm": maximum_whitened,
            "maximum_canonical_unitarity_max_abs": maximum_unitarity,
            "maximum_CAR_spectrum_transport_max_abs": maximum_car,
            "minimum_target_covariance_eigenvalue": minimum_car,
            "maximum_target_covariance_eigenvalue": maximum_car_eigenvalue,
            "minimum_frequency_mixing_offblock_frobenius": offblock,
            "all_channel_kernel_rank_defects_zero": True,
            "channel_residuals": rows,
        },
        "payload": payload_manifest(path, stored_arrays),
        "physical_selector": {
            "status": "OPEN",
            "selected": False,
            "finite_stress_evaluated": False,
            "updated_nulls_or_constraints": None,
            "interface_T01_assigned": None,
            "reason": "the current finite spatial boundary response and Landau endpoint constrain the weighted isometry but provide no mode-resolved frequency kernel or intervening embedding",
            "exact_next_owner": "same-action mode-resolved transmitting boundary scattering/embedding operator that selects one V_c in every retained channel",
        },
        "gate": {
            "weighted_interface_variation_pass": True,
            "channel_block_rank_pass": True,
            "frequency_mixing_CAR_witness_pass": True,
            "physical_frequency_mixing_kernel_selected": False,
            "extended_history_gate_ready_from_this_owner": True,
            "extended_existence_or_nonexistence_claimed": False,
            "coupled_evolution_reopened": False,
        },
        "scope": {
            "diagnostic_duration_or_metric_profile_chosen": False,
            "new_counterflow_or_physical_term_added": False,
            "A_q_Omega_zeta_or_V_full_changed": False,
            "old_unit_radius_stress_copied": False,
            "finite_stress_or_nulls_fabricated": False,
            "metric_timestep_started": False,
            "GR_or_heat_kernel_rederived": False,
        },
        "comparison": {
            "fields": "all JSON fields and artifact arrays",
            "exact": "keys, strings, booleans, integer ranks, hashes, dtypes and shapes",
            "float_atol": 3e-12,
            "float_rtol": 3e-12,
            "exceptions": [],
        },
    }


def load_artifact(record):
    path = ROOT / record["payload"]["path"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != record["payload"]["sha256"]:
        raise ValueError("interface artifact hash mismatch")
    with np.load(path, allow_pickle=False) as payload:
        arrays = {name: np.array(payload[name], copy=True) for name in payload.files}
    return path, arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = load_records()
    if args.check:
        expected = json.loads(OUTPUT.read_text())
        path, arrays = load_artifact(expected)
        actual = calculate(path, arrays, records)
        compare(expected, actual)
        print("tilted weighted-isometry family verified; physical frequency kernel remains OPEN")
        return

    output = args.output or OUTPUT
    if output.exists():
        raise FileExistsError("refusing to overwrite recorded evidence")
    state = load_state(records[0])
    rapidity = float(records[1]["source_selected_frame"]["rapidity"])
    arrays, _rows, _identity_size, _distance = build_witness_arrays(state, rapidity)
    raw = deterministic_npz_bytes(arrays)
    digest = hashlib.sha256(raw).hexdigest()
    ARTIFACT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    path = ARTIFACT_DIRECTORY / f"nsc-tilted-landau-interface.{digest}.npz"
    if path.exists() and path.read_bytes() != raw:
        raise FileExistsError("content-addressed artifact collision")
    path.write_bytes(raw)
    record = calculate(path, arrays, records)
    output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(output)
    print(path)


if __name__ == "__main__":
    main()
