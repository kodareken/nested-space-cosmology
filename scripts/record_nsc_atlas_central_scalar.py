#!/usr/bin/env python3
"""Record a reflected incoming-ray refinement of the formal central scalar.

Only the leading field and incoming integrals are newly evaluated. Saved
wall, cap-clock and exterior-correction records are authenticated explicitly.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

import record_nsc_atlas_cap_response as cap_record
import record_nsc_atlas_lift_clock as clock_record
import record_nsc_atlas_scalar_evolution as exterior_record
from record_nsc_atlas_scalar_inputs import BASE, _plain, _real, _count, _equal, digest, encoded, runtime_versions
from record_nsc_atlas_forcing_cache import _same

ROOT = Path(__file__).resolve().parents[1]
CAP = "results/development/nsc-atlas-cap-response-v1/radial_fine"
CLOCK = "results/development/nsc-atlas-lift-clock-v1"
EXTERIOR = "results/development/nsc-atlas-scalar-evolution-v1"
SCHEMA = "nsc-atlas-central-scalar-comparison-v1"
STATUS = "FORMAL_CENTRAL_SCALAR_INCOMING_REFINEMENT"
SCOPE = "formal central cubic scalar coefficient; sampled refinement, no physical amplitude or uniform remainder"
TIMES = (.125, .375, .5, .625, .875)
RADIAL_ORDERS = (16, 32)
SETTINGS = dict(field_nodes=129, geometry_knots=257, derivative_order=2,
                small_radius_order=8, centre_radius=1e-4, ray_order=16)
INTEGRAL_KINDS = ("forcing", "optical_time", "lambda_plus")
INTEGRAL_SCALARS = ("q", "x", "u_wall", "reflected_profile_defect", "central_profile_defect",
    "regular_center_K1", "first_arrival_margin", "last_arrival_margin", "no_return_margin")
SOURCES = tuple(dict.fromkeys((*exterior_record.SOURCES, *cap_record.SOURCES,
    CAP+"/response.json", CAP+"/manifest.json", EXTERIOR+"/evolution.json", EXTERIOR+"/manifest.json",
    "docs/nsc-atlas-central-scalar-response.md", "src/recursive_horizons/nsc_atlas_central_scalar.py",
    "tests/test_nsc_atlas_central_scalar.py", "scripts/record_nsc_atlas_central_scalar.py",
    "tests/test_record_nsc_atlas_central_scalar.py")))


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def _rows(values):
    lookup = {row["t"]: row for row in values}
    try:
        return [lookup[t] for t in TIMES]
    except KeyError as error:
        raise ValueError("central reconstruction requires five common retained cuts") from error


def input_snapshot():
    cap_record.check(ROOT/CAP, source_context="current")
    exterior_record.check(ROOT/EXTERIOR)  # Also authenticates the selected clock record.
    cap = json.loads((ROOT/CAP/"response.json").read_bytes())
    clock = json.loads((ROOT/CLOCK/"clock.json").read_bytes())
    exterior = json.loads((ROOT/EXTERIOR/"evolution.json").read_bytes())
    expected_cap = {**BASE, "radial_order": 32, "time_panels": 16}
    if (cap["comparison"] != "radial_fine" or cap["settings"][1] != expected_cap
            or clock["settings"][1] != {**BASE, "wall_time_panels": 32}
            or exterior["clock_time_panels"] != 32 or exterior["solutions"][3]["cache_times"] != 33
            or exterior["solutions"][3]["settings"] != {"grid_nodes": 129, **exterior_record.ODE_CASES[3]}):
        raise ValueError("central reconstruction fixed input settings disagree")
    domain = cap["domain"]
    if (domain != clock["domain"] or any(domain[key] != exterior["domain"][key] for key in exterior["domain"])):
        raise ValueError("central reconstruction retained domains disagree")
    caps = _rows(cap["solutions"][1]["cap_samples"])
    walls = _rows(clock["solutions"][1]["wall_samples"])
    boundary = _rows(clock["solutions"][1]["boundary_samples"])
    exterior_rows = [_rows(solution["wall_samples"]) for solution in exterior["solutions"]]
    for i, (c, w, b) in enumerate(zip(caps, walls, boundary)):
        pairs = [(c[key], b[key]) for key in ("V", "x", "e", "e_s")]
        pairs += [(c["V0_s"], w["Vdot"]), (c["q0"], w["reference_q_plus"]),
                  (c["q0"], b["q_plus"]), (w["w"], b["wall_profile"])]
        pairs += [(row[i]["V"], w["V"]) for row in exterior_rows]
        pairs += [(row[i]["U_wall"], w["U"]) for row in exterior_rows]
        if any(not math.isclose(a, v, rel_tol=2e-12, abs_tol=2e-13) for a, v in pairs):
            raise ValueError("central reconstruction changed a common reference or geometric datum")
    files = (CAP+"/response.json", CAP+"/manifest.json", CLOCK+"/clock.json", CLOCK+"/manifest.json",
             EXTERIOR+"/evolution.json", EXTERIOR+"/manifest.json")
    receipt = {"files": {name: digest((ROOT/name).read_bytes()) for name in files},
        "cap_source_context": "current", "cap_solution_index": 1, "clock_solution_index": 1,
        "exterior_solution_index": 3, "cap_radial_order": 32, "cap_ray_order": 16,
        "cap_time_panels": 16, "cap_saved_wall_time_panels": 8, "current_wall_time_panels": 32,
        "new_wall_solves": 0, "new_cap_clock_solves": 0, "new_scalar_propagations": 0}
    return dict(domain=domain, caps=caps, walls=walls, boundary=boundary,
                exterior=exterior_rows, receipt=receipt)


def reconstruct_rows(inputs, i, profile, integrals, qV):
    from recursive_horizons.nsc_atlas_central_scalar import reconstruct_central_scalar
    wall, cap, boundary = inputs["walls"][i], inputs["caps"][i], inputs["boundary"][i]
    result = []
    for rows in inputs["exterior"]:
        row = rows[i]
        result.append(reconstruct_central_scalar(wall, cap, profile, integrals,
            y_wall=row["y_wall"], y_wall_s=wall["Vdot"]*row["y_wall_V"],
            J_minus=boundary["J_minus"], arrival_spline_derivative=qV))
    return result


def leading_compatibility(inputs, i, profile, qV):
    w = inputs["walls"][i]
    return {"field_value_defect": profile["f"]-w["w"],
            "field_spline_tangent_defect": profile["f_prime"]-w["j"]/(qV*w["Vdot"])}


def diagnostics(solutions):
    rows = [solution["samples"] for solution in solutions]
    quantities = ("p1_center", "Delta_C_fixed_s", "Delta_C_fixed_s_s", "Delta_C_fixed_V", "Delta_C_fixed_V_s",
                  "incoming_wall_derivative", "incoming_forcing_term")
    differences = {key: max(abs(a["reconstruction"][key]-b["reconstruction"][key]) for a, b in zip(*rows))
                   for key in quantities}
    axes = []
    for name, a, b in (("cache_time", 0, 1), ("correction_grid", 2, 1), ("ODE_tolerance", 1, 3)):
        axes.append({"axis": name, "exterior_solution_indices": [a, b],
            "maximum_absolute_p1_center_difference": max(abs(row["exterior_controls"][a]["p1_center"]-
                row["exterior_controls"][b]["p1_center"]) for row in rows[1])})
    return {"incoming_radial_differences": differences, "exterior_axes_at_radial32": axes,
        "maximum_absolute_p0_center": max(abs(row["reconstruction"]["p0_center"]) for row in rows[1]),
        "maximum_absolute_p1_center": max(abs(row["reconstruction"]["p1_center"]) for row in rows[1]),
        "maximum_absolute_reused_ray_integral_defect": [max(abs(row["reconstruction"]["reused_ray_integral_defect"])
            for row in samples) for samples in rows],
        "maximum_absolute_saved_Delta_C_fixed_s_difference": [max(abs(row["saved_cap_arrival_difference"])
            for row in samples) for samples in rows],
        "static_null_frame_check": "algebraic identity only; not an independent second solution",
        "uniform_error_enclosure": False}


def compare():
    from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget, solve_prepared_field
    from recursive_horizons.nsc_atlas_bulk_response import BulkBudget
    from recursive_horizons.nsc_atlas_cap_field import CapLeadingField
    from recursive_horizons.nsc_atlas_central_scalar import incoming_response_integrals
    shared = PreparedFieldBudget()
    metrics = BulkBudget(started=shared.started, external=shared)
    inputs = input_snapshot()
    field = solve_prepared_field(129, 257, budget=shared)
    cap = CapLeadingField(field, **{key: SETTINGS[key] for key in ("derivative_order", "small_radius_order", "centre_radius")})
    if any(not _equal(cap.domain[key], value) for key, value in inputs["domain"].items()):
        raise ValueError("new leading field leaves the authenticated cap domain")
    profiles, derivatives, defects = [], [], []
    for i, c in enumerate(inputs["caps"]):
        profile = {key: float(value[0]) for key, value in cap.profile(np.array([c["q0"]]), second=True).items()}
        qV = float(field.geometry["arrival"](c["V"], 1))
        difference = leading_compatibility(inputs, i, profile, qV)
        if (abs(float(field.geometry["wall"](c["V"]))-inputs["walls"][i]["U"]) > 2e-13
                or abs(float(field.geometry["arrival"](c["V"]))-c["q0"]) > 2e-13
                or abs(difference["field_value_defect"]) > 2e-9*max(1., abs(inputs["walls"][i]["w"]))
                or abs(difference["field_spline_tangent_defect"]) > 2e-9*max(1., abs(profile["f_prime"]))):
            raise ValueError("new leading profile differs from the retained wall response")
        profiles.append(profile)
        derivatives.append(qV)
        defects.append(difference)
    solutions = []
    for order in RADIAL_ORDERS:
        samples = []
        for i, t in enumerate(TIMES):
            c = inputs["caps"][i]
            integral = incoming_response_integrals(cap, c["q0"], c["x"],
                radial_order=order, ray_order=16, budget=metrics)
            controls = reconstruct_rows(inputs, i, profiles[i], integral, derivatives[i])
            samples.append(dict(t=t, profile=profiles[i], arrival_spline_derivative=derivatives[i],
                leading_compatibility=defects[i], integrals=integral, reconstruction=controls[3],
                exterior_controls=controls,
                saved_cap_arrival_difference=controls[3]["Delta_C_fixed_s"]-c["Delta_C_fixed_s"]))
        solutions.append(_plain(dict(radial_order=order, ray_order=16, samples=samples)))
    shared.check()
    metrics.check()
    return _plain(dict(schema=SCHEMA, status=STATUS, scope=SCOPE, runtime=runtime_versions(), settings=SETTINGS,
        source_scale=1., domain=inputs["domain"], inputs=inputs["receipt"], solutions=solutions,
        diagnostics=diagnostics(solutions), refinement_is_rigorous_error_enclosure=False,
        profile_evaluation=dict(counts=dict(cap.evaluation_counts), maximum_inverse_bracket_width=cap.max_inverse_bracket_width),
        execution=dict(cpu_seconds=time.process_time()-shared.started, cpu_seconds_max=120.,
            scalar_kernel_values=shared.kernel_values, scalar_kernel_values_max=30_000_000,
            metric_values=metrics.metric_values, source_quadrature_nodes=metrics.source_nodes,
            field_rhs_calls=field.rhs_calls, leading_field_solves=1, wall_solves=0, cap_clock_solves=0,
            scalar_propagations=0, incoming_ray_evaluations=10,
            incoming_ray_nodes=sum(row["integrals"]["ray_nodes"] for s in solutions for row in s["samples"]),
            optical_radial_nodes=sum(row["integrals"]["optical_radial_nodes"] for s in solutions for row in s["samples"]))))


def create(directory):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only record: choose a fresh directory or --check")
    bindings, runtime = source_bindings(), runtime_versions()
    record = compare()
    if bindings != source_bindings() or runtime != runtime_versions() or record["runtime"] != runtime:
        raise RuntimeError("calculation inputs or runtime changed during comparison")
    validate_record(record)
    raw = encoded(record)
    manifest = {"schema": SCHEMA+"-manifest", "source_bindings": bindings,
        "outputs": {"central.json": {"sha256": digest(raw), "bytes": len(raw)}}}
    directory.mkdir(parents=True, exist_ok=False)
    (directory/"central.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return {"created": str(directory), "execution": record["execution"], "diagnostics": record["diagnostics"]}


def check(directory):
    directory = Path(directory)
    path = directory/"manifest.json"
    manifest = json.loads(path.read_bytes())
    if (path.is_symlink() or manifest.get("schema") != SCHEMA+"-manifest"
            or manifest.get("source_bindings") != source_bindings()):
        raise ValueError("central scalar source binding or schema mismatch")
    if ({p.name for p in directory.iterdir()} != {"manifest.json", "central.json"}
            or set(manifest.get("outputs", {})) != {"central.json"}):
        raise ValueError("central scalar inventory mismatch")
    path = directory/"central.json"
    raw = path.read_bytes()
    if path.is_symlink() or manifest["outputs"]["central.json"] != {"sha256": digest(raw), "bytes": len(raw)}:
        raise ValueError("central scalar output mismatch")
    validate_record(json.loads(raw))
    return {"checked": str(directory), "physical_solves": 0, "incoming_integrals_recomputed": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    print(json.dumps(check(args.output) if args.check else create(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
