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
from validation_source_successors import verify_saved_bindings, verify_saved_source

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
    "analytic_regular_center_K1", "first_arrival_margin", "last_arrival_margin", "no_return_margin")
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
    domain = cap["solutions"][1]["domain"]
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
    edge_rows = {row["t"]: row for row in clock["solutions"][1]["boundary_samples"]}
    from recursive_horizons.nsc_atlas_prepared_response import _geometry
    reference = _geometry(257)
    return dict(domain=domain, caps=caps, walls=walls, boundary=boundary,
                exterior=exterior_rows, receipt=receipt,
                profile_edges=[edge_rows[t]["q_plus"] for t in (0., .25, .75, 1.)],
                arrival_spline_derivatives=[float(reference["arrival"](row["V"], 1)) for row in caps])


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
        "maximum_absolute_quiet_recovery_defect": max(abs(value) for samples in rows
            for value in [*samples[0]["profile"].values(), samples[0]["reconstruction"]["p0_center"],
                samples[0]["reconstruction"]["p1_center"],
                *(samples[0]["integrals"][kind+"_integral"] for kind in INTEGRAL_KINDS)]),
        "maximum_absolute_central_frame_algebra_defect": max(abs(row["reconstruction"]["central_frame_algebra_defect"])
            for samples in rows for row in samples),
        "center_sample_scope": "two positive-radius samples per ray; analytic centre value is a separate identity",
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


def _quadrature_shape(q, x, edges, order):
    cuts = [0., x]
    for edge in edges:
        if q > edge:
            for factor in (1., .5):
                r = math.tanh(factor*(q-edge))
                if 0 < r < x:
                    cuts.append(r)
    cuts = np.unique(cuts)
    nodes = np.polynomial.legendre.leggauss(16)[0]
    radii = np.concatenate([a+(b-a)*(1+nodes)/2 for a, b in zip(cuts, cuts[1:])])
    optical = 0
    for r in radii:
        T = q-math.atanh(r)
        radial_cuts = {0., float(r)}
        radial_cuts.update(math.tanh(abs(edge-T)) for edge in edges if 0 < math.tanh(abs(edge-T)) < r)
        optical += order*(len(radial_cuts)-1)
    return radii, optical


def _validate_integrals(integral, inputs, i, order):
    names = {*INTEGRAL_SCALARS, "ray_nodes", "optical_radial_nodes", "center_samples"}
    names.update(kind+suffix for kind in INTEGRAL_KINDS for suffix in
                 ("_integral", "_positive", "_negative", "_sample_min", "_sample_max"))
    if set(integral) != names:
        raise ValueError("central incoming integral inventory mismatch")
    for key in names-{"ray_nodes", "optical_radial_nodes", "center_samples"}:
        _real(integral[key], key)
    q, x, d = inputs["caps"][i]["q0"], inputs["caps"][i]["x"], inputs["domain"]
    expected = dict(q=q, x=x, u_wall=q-2*math.atanh(x), analytic_regular_center_K1=0.,
        first_arrival_margin=q-d["qa"], last_arrival_margin=d["qb"]-q,
        no_return_margin=d["qa"]-q+2*math.atanh(x))
    if any(not _equal(integral[key], value) for key, value in expected.items()):
        raise ValueError("central incoming coordinates or margins differ")
    if min(integral[key] for key in ("first_arrival_margin", "last_arrival_margin", "no_return_margin")) <= 0:
        raise ValueError("central scalar requires strictly interior retained arrivals")
    if integral["reflected_profile_defect"] < 0 or integral["central_profile_defect"] < 0:
        raise ValueError("central profile defect magnitude is negative")
    radii, optical = _quadrature_shape(q, x, inputs["profile_edges"], order)
    if not isinstance(integral["center_samples"], list) or len(integral["center_samples"]) != 2:
        raise ValueError("central positive-radius sample inventory mismatch")
    for sample, r in zip(integral["center_samples"], radii[:2]):
        if set(sample) != {"r", "T", "W1", "K1", "K1_over_r_squared"}:
            raise ValueError("central positive-radius sample fields mismatch")
        for key, value in sample.items():
            _real(value, key)
        if (sample["r"] <= 0 or not _equal(sample["r"], float(r))
                or not _equal(sample["T"], q-math.atanh(r))
                or not _equal(sample["K1"], (1-r*r)*sample["W1"]/4)
                or not _equal(sample["K1_over_r_squared"], sample["K1"]/(r*r))):
            raise ValueError("central positive-radius sample identity mismatch")
    for key, value in (("ray_nodes", len(radii)), ("optical_radial_nodes", optical)):
        _count(integral[key], key, 1_000_000)
        if integral[key] != value:
            raise ValueError("central incoming quadrature count mismatch")
    for kind in INTEGRAL_KINDS:
        value, positive, negative, low, high = [integral[kind+suffix] for suffix in
            ("_integral", "_positive", "_negative", "_sample_min", "_sample_max")]
        if (positive < 0 or negative > 0 or low > high or value != positive+negative
                or value/x < low-1e-12*max(1., abs(low)) or value/x > high+1e-12*max(1., abs(high))):
            raise ValueError("central signed integral or sampled extrema disagree")


def validate_record(record):
    try:
        if (record.get("schema") != SCHEMA or record.get("status") != STATUS or record.get("scope") != SCOPE
                or record.get("refinement_is_rigorous_error_enclosure") is not False
                or record.get("source_scale") != 1. or isinstance(record.get("source_scale"), bool)):
            raise ValueError("central scalar declaration mismatch")
        _same(record["settings"], SETTINGS, "central settings")
        runtime = record["runtime"]
        if (set(runtime) != {"python", "numpy", "scipy"}
                or any(not isinstance(v, str) or not v or not v[0].isdigit() or "." not in v for v in runtime.values())):
            raise ValueError("central scalar runtime mismatch")
        inputs = input_snapshot()
        _same(record["domain"], inputs["domain"], "central domain")
        _same(record["inputs"], inputs["receipt"], "central input authentication")
        if not isinstance(record["solutions"], list) or len(record["solutions"]) != 2:
            raise ValueError("central scalar requires two incoming radial settings")
        for solution, order in zip(record["solutions"], RADIAL_ORDERS):
            if (type(solution["radial_order"]) is not int or solution["radial_order"] != order
                    or type(solution["ray_order"]) is not int or solution["ray_order"] != 16):
                raise ValueError("central scalar incoming settings mismatch")
            samples = solution["samples"]
            if len(samples) != 5 or [row["t"] for row in samples] != list(TIMES):
                raise ValueError("central scalar five-cut inventory mismatch")
            for i, row in enumerate(samples):
                if set(row["profile"]) != {"f", "f_prime", "f_second"}:
                    raise ValueError("central leading profile inventory mismatch")
                for key, value in row["profile"].items():
                    _real(value, key)
                qV = row["arrival_spline_derivative"]
                _real(qV, "arrival spline derivative")
                if qV <= 0:
                    raise ValueError("central arrival map must increase")
                if not _equal(qV, inputs["arrival_spline_derivatives"][i]):
                    raise ValueError("central arrival spline differs from the reference geometry")
                leading = leading_compatibility(inputs, i, row["profile"], qV)
                _same(row["leading_compatibility"], leading, "leading profile compatibility")
                if (abs(leading["field_value_defect"]) > 2e-9*max(1., abs(inputs["walls"][i]["w"]))
                        or abs(leading["field_spline_tangent_defect"]) > 2e-9*max(1., abs(row["profile"]["f_prime"]))):
                    raise ValueError("central leading profile disagrees with saved field")
                _validate_integrals(row["integrals"], inputs, i, order)
                controls = reconstruct_rows(inputs, i, row["profile"], row["integrals"], qV)
                _same(row["exterior_controls"], controls, "central exterior-axis reconstruction")
                _same(row["reconstruction"], controls[3], "central main reconstruction")
                if i == 0 and max(abs(value) for value in [*row["profile"].values(),
                        *(row["integrals"][kind+"_integral"] for kind in INTEGRAL_KINDS),
                        *(control[key] for control in controls for key in ("p0_center", "p1_center"))]) > 1e-24:
                    raise ValueError("central scalar lost the quiet recovery control")
                _real(row["saved_cap_arrival_difference"], "saved cap arrival difference")
                if not _equal(row["saved_cap_arrival_difference"], controls[3]["Delta_C_fixed_s"]-inputs["caps"][i]["Delta_C_fixed_s"]):
                    raise ValueError("central saved cap arrival comparison differs")
        for left, right in zip(*(solution["samples"] for solution in record["solutions"])):
            for key in ("profile", "arrival_spline_derivative", "leading_compatibility"):
                if left[key] != right[key]:
                    raise ValueError("central radial comparison changed the leading profile")
        _same(record["diagnostics"], diagnostics(record["solutions"]), "central diagnostics")
        profile = record["profile_evaluation"]
        if set(profile["counts"]) != cap_record.PROFILE_COUNTS:
            raise ValueError("central profile count inventory mismatch")
        for key, value in profile["counts"].items():
            _count(value, key)
        c = profile["counts"]
        if not c["second_points"] <= c["wall_points"] <= c["profile_requested"] or c["cache_hits"] > c["profile_requested"]:
            raise ValueError("central profile count mismatch")
        _real(profile["maximum_inverse_bracket_width"], "inverse bracket", True)
        if profile["maximum_inverse_bracket_width"] > 4*np.finfo(float).eps*max(abs(record["domain"]["Va"]), abs(record["domain"]["Vb"])):
            raise ValueError("central inverse bracket mismatch")
        cost = record["execution"]
        _real(cost["cpu_seconds"], "CPU time", True)
        _real(cost["cpu_seconds_max"], "CPU ceiling")
        _count(cost["scalar_kernel_values_max"], "kernel ceiling")
        if cost["cpu_seconds"] > 120 or cost["cpu_seconds_max"] != 120 or cost["scalar_kernel_values_max"] != 30_000_000:
            raise ValueError("central scalar budget changed or exceeded")
        for key in ("scalar_kernel_values", "metric_values", "source_quadrature_nodes"):
            _count(cost[key], key, 30_000_000)
        if cost["metric_values"]+cost["source_quadrature_nodes"] > 30_000_000:
            raise ValueError("central metric budget exceeded")
        _count(cost["field_rhs_calls"], "field RHS calls", 6000)
        for key, value in (("leading_field_solves", 1), ("wall_solves", 0), ("cap_clock_solves", 0),
                            ("scalar_propagations", 0), ("incoming_ray_evaluations", 10)):
            if type(cost[key]) is not int or cost[key] != value:
                raise ValueError("central reuse accounting mismatch")
        for key, samplekey in (("incoming_ray_nodes", "ray_nodes"), ("optical_radial_nodes", "optical_radial_nodes")):
            _count(cost[key], key, 10_000_000)
            if cost[key] != sum(row["integrals"][samplekey] for s in record["solutions"] for row in s["samples"]):
                raise ValueError("central incoming node accounting mismatch")
    except (KeyError, TypeError, IndexError) as error:
        raise ValueError("incomplete central scalar record") from error
    return True


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
            or not verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=SOURCES)):
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
