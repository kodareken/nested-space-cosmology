"""Record sourced moments for unequal cap-wall bodies A and B.

A is the preserved closed one-body average. B is a new sourced cap, wall
and bridge with the same action and a different mass and taper. No connector
metric, evolution, Lambda scan or endpoint solve is evaluated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np

from recursive_horizons.nsc_atlas_connection_slice import compact_body_summary
from recursive_horizons.nsc_atlas_geometric_average import geometric_slice_average


ROOT = Path(__file__).resolve().parents[1]
PRESERVED_A = "results/summaries/nsc-atlas-geometric-average-v1.json"
SOURCE_PATHS = (
    "scripts/record_nsc_atlas_unequal_children.py",
    "src/recursive_horizons/nsc_atlas_connection_slice.py",
    "src/recursive_horizons/nsc_atlas_geometric_average.py",
    "src/recursive_horizons/nsc_atlas_moving_boundaries.py",
    "tests/test_nsc_atlas_connection_slice.py",
    "tests/test_nsc_atlas_geometric_average.py",
    PRESERVED_A,
)
CPU_SECONDS = 30
SAMPLE_LIMIT = 100_000
IDENTITY_LIMIT = 1e-8
FINITE_DIFFERENCE_LIMIT = 1e-7
MOMENT_COMPARISON_RTOL = 1e-7
BODY_A = {"id": "A", "mass_over_ell": .4, "taper_stop_over_ell": 1.}
BODY_B = {"id": "B", "mass_over_ell": .6, "taper_stop_over_ell": 1.3}
MOMENT_KEYS = (
    "V_over_ell3", "I_theta_over_ell2", "J_over_ell",
    "intrinsic_curvature_over_ell", "curvature_bulk_over_ell",
    "curvature_wall_delta_over_ell", "cap_energy_over_ell", "wall_energy_over_ell",
    "cap_stress_trace_over_ell", "wall_stress_trace_over_ell",
)


def frozen_sources():
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    hashes = {}
    for name in SOURCE_PATHS:
        raw = (ROOT/name).read_bytes()
        try:
            committed = subprocess.check_output(["git", "show", f"{commit}:{name}"], cwd=ROOT)
        except subprocess.CalledProcessError as error:
            raise ValueError("unequal-children source is not frozen: "+name) from error
        if raw != committed:
            raise ValueError("unequal-children source is not frozen: "+name)
        hashes[name] = hashlib.sha256(raw).hexdigest()
    return commit, hashes


def moment_block(row):
    integrals = row["integrals_in_ell_units"]
    moments = {
        "V_over_ell3": row["volume_over_ell3"],
        "I_theta_over_ell2": integrals["theta"],
        "J_over_ell": 2*integrals["theta_squared"]/3-2*integrals["shear_squared"],
        "intrinsic_curvature_over_ell": integrals["curvature_total"],
        "curvature_bulk_over_ell": integrals["curvature_bulk"],
        "curvature_wall_delta_over_ell": integrals["curvature_wall_delta"],
        "cap_energy_over_ell": integrals["energy_cap"],
        "wall_energy_over_ell": integrals["energy_wall"],
        "cap_stress_trace_over_ell": integrals["stress_trace_cap"],
        "wall_stress_trace_over_ell": integrals["stress_trace_wall"],
    }
    total = moments["curvature_bulk_over_ell"]+moments["curvature_wall_delta_over_ell"]
    if not math.isclose(total, moments["intrinsic_curvature_over_ell"], rel_tol=0., abs_tol=1e-9):
        raise ValueError("curvature total does not include the wall delta")
    return moments


def preserved_body_a():
    payload = json.loads((ROOT/PRESERVED_A).read_text(encoding="utf-8"))
    if payload.get("schema") != "nsc-atlas-geometric-average-summary-v1":
        raise ValueError("preserved A summary schema mismatch")
    reference = payload["reference"]
    if reference.get("schema") != "nsc-atlas-geometric-slice-average-v1":
        raise ValueError("preserved A average schema mismatch")
    domain = reference["domain"]
    expected = {"H_ell": 1., "k_ell": 2., "M_over_ell": .4,
                "Rwall_over_ell": 1.5, "Rcut_right_over_ell": 4.}
    for key, value in expected.items():
        if domain.get(key) != value:
            raise ValueError("preserved A domain is not the closed reference: "+key)
    moments = moment_block(reference)
    stored = payload["time_symmetric_vacuum_extension_condition"]["J_over_ell"]
    if not math.isclose(moments["J_over_ell"], stored, rel_tol=0., abs_tol=1e-12):
        raise ValueError("preserved A moment J does not match its stored extension value")
    return {
        "id": "A", "mass_over_ell": BODY_A["mass_over_ell"],
        "taper_stop_over_ell": BODY_A["taper_stop_over_ell"],
        "taper_source": "historical default of the closed average; its saved domain dict does not repeat the taper",
        "quadrature": reference["quadrature"], "clock": reference["clock"],
        "recomputed": False, "moments": moments,
    }


def require_identity(name, value):
    if not math.isfinite(value) or abs(value) > IDENTITY_LIMIT:
        raise ValueError(name+" failed its construction identity")


def require_body(summary, body):
    inputs = summary["inputs"]
    for key, value in (("mass_over_ell", body["mass_over_ell"]),
                       ("taper_stop_over_ell", body["taper_stop_over_ell"]),
                       ("H_ell", 1.), ("k_ell", 2.),
                       ("Rwall_over_ell", 1.5), ("Rcut_right_over_ell", 4.)):
        if not math.isclose(inputs[key], value, rel_tol=0., abs_tol=1e-15):
            raise ValueError(body["id"]+" input drifted: "+key)
    if inputs["xi"] != "1/6" or inputs["quiet_scalar_amplitude"] != 0. or inputs["scalar_field_solves"] != 0:
        raise ValueError(body["id"]+" changed the shared quiet scalar action")
    checks = summary["checks"]
    for key in (
        "boost_unit_max_error", "cap_reconstruction_max_error", "exterior_reconstruction_max_error",
        "surface_energy_projection_ell_error", "surface_momentum_projection_ell_error",
        "shell_hamiltonian_delta_ell_error", "shell_momentum_delta_ell_error",
        "vacuum_hamiltonian_max_ell2", "vacuum_momentum_max_ell2",
        "transition_mass_max_error_over_ell", "bridge_sample_mass_max_error_over_ell",
        "embedding_spatial_norm_max_error", "embedding_future_normal_norm_max_error",
        "taper_bridge_overlap_max_error", "throat_radius_error_over_ell", "throat_Rprime_abs",
        "throat_Rsecond_error_ell", "throat_scalar_ell2", "right_end_scalar_ell2",
        "right_end_mass_error_over_ell", "transition_32_64_distance_max_error_over_ell",
    ):
        require_identity(body["id"]+" "+key, checks[key])
    if not checks["right_end_finite"] or checks["right_end_Rprime"] <= 0 or not checks["time_symmetric_join"]:
        raise ValueError(body["id"]+" right end is not a positive time-symmetric vacuum end")
    if not math.isclose(checks["right_end_R_over_ell"], 6., rel_tol=0., abs_tol=1e-12):
        raise ValueError(body["id"]+" right display sphere drifted")
    topology = summary["topology"]
    if not topology["regular_centre"] or topology["centre_is_ae_end"]:
        raise ValueError(body["id"]+" does not have a regular non-asymptotic centre")
    if not topology["vacuum_chart_has_two_asymptotic_directions"]:
        raise ValueError(body["id"]+" vacuum chart did not show two asymptotic directions")
    if topology["left_vacuum_end_included"] or topology["ae_ends"] != 1 or not topology["selected_right_end"]:
        raise ValueError(body["id"]+" does not have exactly one selected right AE end")
    if topology["other_ae_ends_after_selecting_right"] != 0 or not topology["selected_end_is_time_symmetric_vacuum"]:
        raise ValueError(body["id"]+" still retains another AE end")
    for row in checks["finite_difference_constraints"]:
        for key in ("hamiltonian_max_ell2", "momentum_max_ell2", "Rsecond_max_ell_error"):
            if not math.isfinite(row[key]) or abs(row[key]) > FINITE_DIFFERENCE_LIMIT:
                raise ValueError(body["id"]+" independent radial check failed: "+key)


def plain(value):
    if isinstance(value, dict):
        return {str(key): plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        if len(value) > 8:
            raise ValueError("compact record retained a sample array")
        return [plain(item) for item in value]
    if isinstance(value, np.floating):
        value = float(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("non-finite compact result")
        return value
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise ValueError("result is not a compact scalar record")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("use a new output filename; preserve earlier results")
    started = time.process_time()
    commit, hashes = frozen_sources()
    reference = preserved_body_a()
    summary_a = compact_body_summary(**{key: BODY_A[key] for key in ("mass_over_ell", "taper_stop_over_ell")})
    summary_b = compact_body_summary(**{key: BODY_B[key] for key in ("mass_over_ell", "taper_stop_over_ell")})
    require_body(summary_a, BODY_A)
    require_body(summary_b, BODY_B)
    sourced_b = geometric_slice_average(mass_over_ell=BODY_B["mass_over_ell"],
                                        taper_stop_over_ell=BODY_B["taper_stop_over_ell"])
    if sourced_b["quadrature"] != "gauss64" or sourced_b["clock"] != reference["clock"]:
        raise ValueError("B moments left the declared normal clock or default quadrature")
    domain = sourced_b["domain"]
    expected_domain = {"H_ell": 1., "k_ell": 2., "Rwall_over_ell": 1.5,
                       "Rcut_right_over_ell": 4., "M_over_ell": .6, "taper_stop_over_ell": 1.3}
    if any(not math.isclose(domain[key], value, rel_tol=0., abs_tol=1e-15)
           for key, value in expected_domain.items()):
        raise ValueError("B average domain is not the fixed sourced body")
    for key in ("hamiltonian_ell2", "bridge_volume_ell3_difference", "wall_trace_ell",
                "scalar_fitting_acceleration_ell2"):
        require_identity("B average "+key, sourced_b["checks"][key])
    moments_b = moment_block(sourced_b)
    comparison_b = geometric_slice_average(
        mass_over_ell=BODY_B["mass_over_ell"],
        taper_stop_over_ell=BODY_B["taper_stop_over_ell"], quadrature="gauss32")
    moments32 = moment_block(comparison_b)
    moment_differences = {key: abs(moments_b[key]-moments32[key]) for key in MOMENT_KEYS}
    if any(moment_differences[key] > MOMENT_COMPARISON_RTOL*max(1.,abs(moments_b[key]))
           for key in MOMENT_KEYS):
        raise ValueError("B 32/64 moment comparison failed")
    samples = (summary_a["geometric_samples"]+summary_b["geometric_samples"]
               +sourced_b["spatial_evaluations"]+comparison_b["spatial_evaluations"])
    if samples > SAMPLE_LIMIT:
        raise RuntimeError("unequal-children sample budget exhausted")
    if time.process_time()-started > CPU_SECONDS:
        raise RuntimeError("unequal-children CPU budget exhausted")
    summed = {key: reference["moments"][key]+moments_b[key] for key in MOMENT_KEYS}
    record = {
        "schema": "nsc-atlas-unequal-children-v1",
        "admission_owner": "docs/prior-art-and-open-claim.md#unequal-sourced-inputs-for-the-common-parent-18-october-2026",
        "method_commit": commit,
        "source_sha256": hashes,
        "command": "python scripts/run_python.py scripts/record_nsc_atlas_unequal_children.py --output <path>",
        "inputs": {
            "shared_action": {
                "H_ell": 1., "k_ell": 2., "xi": "1/6", "quiet_scalar_amplitude": 0.,
                "wall_law": "pressure=-sigma, fixed positive sigma",
                "note": "Same action and couplings. Distinct masses and taper stops are states and slicing, not new couplings.",
            },
            "bodies": [BODY_A, BODY_B],
            "Rwall_over_ell": 1.5, "Rcut_right_over_ell": 4.,
            "clock": reference["clock"],
        },
        "bounds": {
            "cpu_seconds_max": CPU_SECONDS,
            "sampled_geometric_points_max": SAMPLE_LIMIT,
            "sampled_geometric_points": samples,
            "charge": "One pointwise geometric jet or quadrature sample, counting repeated evaluations.",
            "field_solves": 0, "physical_time_steps": 0,
            "mass_or_cutoff_scan": False,
            "identity_abs_gate": IDENTITY_LIMIT,
            "finite_difference_abs_gate": FINITE_DIFFERENCE_LIMIT,
            "moment_comparison_relative_gate": MOMENT_COMPARISON_RTOL,
            "finite_difference_and_quadrature_are_rigorous_bounds": False,
        },
        "bodies": {
            "A": reference,
            "B": {
                "id": "B", "mass_over_ell": .6, "taper_stop_over_ell": 1.3,
                "quadrature": "gauss64", "clock": sourced_b["clock"],
                "recomputed": True, "moments": moments_b,
                "wall": summary_b["wall"], "topology": summary_b["topology"],
                "rate_convention": "Rdot=dR/d(wall proper time), in geometric units; not a local three-speed.",
                "moment_32_64_absolute_differences": moment_differences,
                "average_checks": {
                    "hamiltonian_ell2": sourced_b["checks"]["hamiltonian_ell2"],
                    "bridge_volume_ell3_difference": sourced_b["checks"]["bridge_volume_ell3_difference"],
                    "wall_trace_ell": sourced_b["checks"]["wall_trace_ell"],
                    "scalar_fit_is_a_four_tensor": False,
                    "lambda_mechanism_scan": False,
                },
            },
        },
        "checks": {"A": summary_a["checks"], "B": summary_b["checks"]},
        "sum_without_connector": {
            "quadratures": "preserved adaptive A plus gauss64 B, compared with gauss32 B",
            "moments": summed,
        },
        "connector": {
            "V_C_over_ell3": None, "I_C_over_ell2": None, "J_C_over_ell": None,
            "status": "absent",
            "reason": "No common connector metric is evaluated, so its volume and extrinsic-curvature moments stay absent.",
        },
        "conditional_end_gluing": {
            "theorem": "1004.1355v1 Theorem 1 applied only to the selected smooth vacuum AE ends",
            "each_body_regular_centre": True,
            "each_body_ae_ends": 1,
            "vacuum_chart_ae_directions": 2,
            "sourced_body_includes_only_the_right_direction": True,
            "retained_other_ae_ends": 0,
            "designated_common_outer_ae_ends": 1,
            "total_parent_ae_ends": 1,
            "complete_vacuum_bridge_other_ae_end_after_selecting_one": 1,
            "exact_sum_adm_charge": False,
            "numeric_separation": None,
            "shared_future": False,
            "sibling_signal": False,
            "connector_metric_evaluated": False,
            "reason": (
                "Each sourced body has a finite-distance regular centre and a bridge parameter "
                "that is finite and negative at the taper, so the vacuum chart's left asymptotic "
                "direction is not part of the body. Theorem 1 combines the two selected right ends "
                "into one exterior and would retain any unselected end. None remain. A complete "
                "vacuum bridge has both chart directions, so selecting one would have retained the other."
            ),
        },
        "scope": (
            "Fixed two-body source moments in the declared normal clock, with connector moments absent. "
            "Not a connector metric, formation or evolution, a generated Lambda, or an observational fit. "
            "The isolated scalar fit is not a four-tensor."
        ),
    }
    record["runtime"] = {"cpu_seconds": time.process_time()-started}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(plain(record), indent=2, sort_keys=True, allow_nan=False)+"\n",
                           encoding="utf-8", newline="\n")
    print(f"Wrote {args.output}; {time.process_time()-started:.3f} CPU seconds; no evolution.")


if __name__ == "__main__":
    main()
