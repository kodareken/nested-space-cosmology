#!/usr/bin/env python3
"""Refine the incoming ray rule against the authenticated central record."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

import record_nsc_atlas_central_scalar as base
from record_nsc_atlas_scalar_inputs import _plain, _real, _count, _equal, digest, encoded, runtime_versions
from record_nsc_atlas_forcing_cache import _same

ROOT = Path(__file__).resolve().parents[1]
PREDECESSOR = "results/development/nsc-atlas-central-scalar-v1"
SCHEMA = "nsc-atlas-central-ray-comparison-v1"
STATUS = "FORMAL_CENTRAL_SCALAR_RAY_REFINEMENT"
SCOPE = "formal central cubic scalar coefficient; incoming ray refinement, no physical amplitude or uniform remainder"
SETTINGS = {**base.SETTINGS, "radial_order": 32, "ray_order": 32}
SOURCES = (*base.SOURCES, PREDECESSOR+"/central.json", PREDECESSOR+"/manifest.json",
           "scripts/record_nsc_atlas_central_ray.py", "tests/test_record_nsc_atlas_central_ray.py")


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def input_snapshot():
    base.check(ROOT/PREDECESSOR)
    previous = json.loads((ROOT/PREDECESSOR/"central.json").read_bytes())
    inputs = base.input_snapshot()
    if (previous["settings"] != base.SETTINGS or previous["domain"] != inputs["domain"]
            or previous["solutions"][1]["radial_order"] != 32
            or previous["solutions"][1]["ray_order"] != 16):
        raise ValueError("central ray predecessor changed the fixed preparation")
    receipt = dict(predecessor_files={name: digest((ROOT/PREDECESSOR/name).read_bytes())
                                     for name in ("central.json", "manifest.json")},
        predecessor_solution_index=1, predecessor_execution_reused=False,
        retained_inputs=inputs["receipt"])
    return inputs, previous["solutions"][1], receipt


def diagnostics(solutions):
    coarse, fine = (s["samples"] for s in solutions)
    keys = ("p1_center", "Delta_C_fixed_s", "Delta_C_fixed_s_s", "Delta_C_fixed_V", "Delta_C_fixed_V_s",
            "exterior_trace_term", "arrival_rate_term", "arrival_shift_term", "outgoing_wall_term", "incoming_forcing_term")
    result = base.diagnostics(solutions)
    result.pop("incoming_radial_differences")
    result["incoming_ray_differences"] = {
        key: max(abs(a["reconstruction"][key]-b["reconstruction"][key]) for a, b in zip(coarse, fine)) for key in keys}
    result["incoming_integral_differences"] = {
        key: max(abs(a["integrals"][key+"_integral"]-b["integrals"][key+"_integral"])
                 for a, b in zip(coarse, fine)) for key in base.INTEGRAL_KINDS}
    result["coarse_side"] = "authenticated saved radial32/ray16; its earlier costs are not counted as new work"
    return result


def compare():
    from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget, solve_prepared_field
    from recursive_horizons.nsc_atlas_bulk_response import BulkBudget
    from recursive_horizons.nsc_atlas_cap_field import CapLeadingField
    from recursive_horizons.nsc_atlas_central_scalar import incoming_response_integrals
    shared = PreparedFieldBudget()
    metrics = BulkBudget(started=shared.started, external=shared)
    inputs, coarse, receipt = input_snapshot()
    field = solve_prepared_field(129, 257, budget=shared)
    cap = CapLeadingField(field, **{key: SETTINGS[key] for key in ("derivative_order", "small_radius_order", "centre_radius")})
    _same(cap.domain, inputs["domain"], "new cap domain")
    samples, profiles = [], []
    for i, old in enumerate(coarse["samples"]):
        c = inputs["caps"][i]
        profile = {key: float(value[0]) for key, value in cap.profile(np.array([c["q0"]]), second=True).items()}
        _same(profile, old["profile"], "recomputed leading profile")
        _same(float(field.geometry["arrival"](c["V"], 1)), old["arrival_spline_derivative"], "arrival tangent")
        _same(float(field.geometry["arrival"](c["V"])), c["q0"], "arrival coordinate")
        _same(float(field.geometry["wall"](c["V"])), inputs["walls"][i]["U"], "wall coordinate")
        profiles.append(profile)
        integral = incoming_response_integrals(cap, c["q0"], c["x"],
            radial_order=32, ray_order=32, budget=metrics)
        controls = base.reconstruct_rows(inputs, i, old["profile"], integral, old["arrival_spline_derivative"])
        samples.append({**{key: old[key] for key in ("t", "profile", "arrival_spline_derivative", "leading_compatibility")},
            "integrals": integral, "reconstruction": controls[3], "exterior_controls": controls,
            "saved_cap_arrival_difference": controls[3]["Delta_C_fixed_s"]-c["Delta_C_fixed_s"]})
    solutions = _plain([coarse, dict(radial_order=32, ray_order=32, samples=samples)])
    shared.check()
    metrics.check()
    return _plain(dict(schema=SCHEMA, status=STATUS, scope=SCOPE, runtime=runtime_versions(), settings=SETTINGS,
        source_scale=1., domain=inputs["domain"], inputs=receipt, solutions=solutions,
        recomputed_profiles=profiles, diagnostics=diagnostics(solutions), refinement_is_rigorous_error_enclosure=False,
        profile_evaluation=dict(counts=dict(cap.evaluation_counts), maximum_inverse_bracket_width=cap.max_inverse_bracket_width),
        execution=dict(cpu_seconds=time.process_time()-shared.started, cpu_seconds_max=120.,
            scalar_kernel_values=shared.kernel_values, scalar_kernel_values_max=30_000_000,
            metric_values=metrics.metric_values, source_quadrature_nodes=metrics.source_nodes,
            field_rhs_calls=field.rhs_calls, leading_field_solves=1, wall_solves=0, cap_clock_solves=0,
            scalar_propagations=0, incoming_ray_evaluations=5, reused_incoming_ray_evaluations=5,
            incoming_ray_nodes=sum(row["integrals"]["ray_nodes"] for row in samples),
            optical_radial_nodes=sum(row["integrals"]["optical_radial_nodes"] for row in samples))))


def quadrature_shape(q, x, edges):
    """Count the new 32-point rules without querying a physical field."""
    cuts = [0., x]
    for edge in edges:
        if q > edge:
            for factor in (1., .5):
                r = math.tanh(factor*(q-edge))
                if 0 < r < x:
                    cuts.append(r)
    cuts = np.unique(cuts)
    nodes = np.polynomial.legendre.leggauss(32)[0]
    radii = np.concatenate([a+(b-a)*(1+nodes)/2 for a, b in zip(cuts, cuts[1:])])
    optical = 0
    for r in radii:
        T = q-float(np.arctanh(r))
        radial_cuts = {0., float(r)}
        radial_cuts.update(math.tanh(abs(float(edge)-T)) for edge in edges
                           if 0 < math.tanh(abs(float(edge)-T)) < r)
        optical += 32*(len(radial_cuts)-1)
    return radii, optical


def validate_integrals(integral, inputs, i):
    names = {*base.INTEGRAL_SCALARS, "ray_nodes", "optical_radial_nodes", "center_samples"}
    names.update(kind+suffix for kind in base.INTEGRAL_KINDS for suffix in
                 ("_integral", "_positive", "_negative", "_sample_min", "_sample_max"))
    if set(integral) != names:
        raise ValueError("central ray integral inventory mismatch")
    for key in names-{"ray_nodes", "optical_radial_nodes", "center_samples"}:
        _real(integral[key], key)
    q, x, d = inputs["caps"][i]["q0"], inputs["caps"][i]["x"], inputs["domain"]
    expected = dict(q=q, x=x, u_wall=q-2*math.atanh(x), analytic_regular_center_K1=0.,
        first_arrival_margin=q-d["qa"], last_arrival_margin=d["qb"]-q,
        no_return_margin=d["qa"]-q+2*math.atanh(x))
    if any(not _equal(integral[key], value) for key, value in expected.items()):
        raise ValueError("central ray coordinates or margins differ")
    if (min(integral[key] for key in ("first_arrival_margin", "last_arrival_margin", "no_return_margin")) <= 0
            or min(integral[key] for key in ("reflected_profile_defect", "central_profile_defect")) < 0):
        raise ValueError("central ray leaves its retained domain or has a negative defect")
    radii, optical = quadrature_shape(q, x, inputs["profile_edges"])
    if not isinstance(integral["center_samples"], list) or len(integral["center_samples"]) != 2:
        raise ValueError("central ray positive-radius inventory mismatch")
    for sample, r in zip(integral["center_samples"], radii[:2]):
        if set(sample) != {"r", "T", "W1", "K1", "K1_over_r_squared"}:
            raise ValueError("central ray positive-radius fields mismatch")
        for key, value in sample.items():
            _real(value, key)
        if (sample["r"] <= 0 or not _equal(sample["r"], float(r))
                or not _equal(sample["T"], q-float(np.arctanh(r)))
                or not _equal(sample["K1"], (1-r*r)*sample["W1"]/4)
                or not _equal(sample["K1_over_r_squared"], sample["K1"]/(r*r))):
            raise ValueError("central ray positive-radius identity mismatch")
    for key, value in (("ray_nodes", len(radii)), ("optical_radial_nodes", optical)):
        _count(integral[key], key, 1_000_000)
        if integral[key] != value:
            raise ValueError("central ray quadrature count mismatch")
    for kind in base.INTEGRAL_KINDS:
        value, positive, negative, low, high = [integral[kind+suffix] for suffix in
            ("_integral", "_positive", "_negative", "_sample_min", "_sample_max")]
        if (positive < 0 or negative > 0 or low > high or value != positive+negative
                or value/x < low-1e-12*max(1., abs(low)) or value/x > high+1e-12*max(1., abs(high))):
            raise ValueError("central ray signed integral or sampled extrema disagree")


def validate_record(record):
    try:
        if (record.get("schema") != SCHEMA or record.get("status") != STATUS or record.get("scope") != SCOPE
                or record.get("refinement_is_rigorous_error_enclosure") is not False
                or record.get("source_scale") != 1. or isinstance(record.get("source_scale"), bool)):
            raise ValueError("central ray declaration mismatch")
        _same(record["settings"], SETTINGS, "central ray settings")
        runtime = record["runtime"]
        if (set(runtime) != {"python", "numpy", "scipy"}
                or any(not isinstance(v, str) or not v or not v[0].isdigit() or "." not in v for v in runtime.values())):
            raise ValueError("central ray runtime mismatch")
        inputs, coarse, receipt = input_snapshot()
        _same(record["domain"], inputs["domain"], "central ray domain")
        _same(record["inputs"], receipt, "central ray input authentication")
        solutions = record["solutions"]
        if not isinstance(solutions, list) or len(solutions) != 2 or solutions[0] != coarse:
            raise ValueError("central ray requires the unchanged authenticated coarse solution")
        fine = solutions[1]
        if (type(fine["radial_order"]) is not int or fine["radial_order"] != 32
                or type(fine["ray_order"]) is not int or fine["ray_order"] != 32
                or len(fine["samples"]) != 5 or [row["t"] for row in fine["samples"]] != list(base.TIMES)):
            raise ValueError("central ray refinement inventory mismatch")
        if not isinstance(record["recomputed_profiles"], list) or len(record["recomputed_profiles"]) != 5:
            raise ValueError("central ray recomputed profile inventory mismatch")
        for i, (old, row) in enumerate(zip(coarse["samples"], fine["samples"])):
            for key in ("profile", "arrival_spline_derivative", "leading_compatibility"):
                if row[key] != old[key]:
                    raise ValueError("central ray refinement changed a retained leading input")
            _same(record["recomputed_profiles"][i], old["profile"], "recomputed leading profile")
            validate_integrals(row["integrals"], inputs, i)
            controls = base.reconstruct_rows(inputs, i, old["profile"], row["integrals"], old["arrival_spline_derivative"])
            _same(row["exterior_controls"], controls, "central ray exterior axes")
            _same(row["reconstruction"], controls[3], "central ray reconstruction")
            _real(row["saved_cap_arrival_difference"], "saved cap arrival difference")
            if not _equal(row["saved_cap_arrival_difference"], controls[3]["Delta_C_fixed_s"]-inputs["caps"][i]["Delta_C_fixed_s"]):
                raise ValueError("central ray saved cap comparison differs")
            if i == 0 and max(abs(value) for value in [*row["profile"].values(),
                    *(row["integrals"][kind+"_integral"] for kind in base.INTEGRAL_KINDS),
                    *(control[key] for control in controls for key in ("p0_center", "p1_center"))]) > 1e-24:
                raise ValueError("central ray lost the quiet recovery control")
        _same(record["diagnostics"], diagnostics(solutions), "central ray diagnostics")
        profile = record["profile_evaluation"]
        if set(profile["counts"]) != base.cap_record.PROFILE_COUNTS:
            raise ValueError("central ray profile inventory mismatch")
        for key, value in profile["counts"].items():
            _count(value, key)
        c = profile["counts"]
        if not c["second_points"] <= c["wall_points"] <= c["profile_requested"] or c["cache_hits"] > c["profile_requested"]:
            raise ValueError("central ray profile count mismatch")
        _real(profile["maximum_inverse_bracket_width"], "inverse bracket", True)
        if profile["maximum_inverse_bracket_width"] > 4*np.finfo(float).eps*max(abs(record["domain"]["Va"]), abs(record["domain"]["Vb"])):
            raise ValueError("central ray inverse bracket mismatch")
        cost = record["execution"]
        _real(cost["cpu_seconds"], "CPU time", True)
        _real(cost["cpu_seconds_max"], "CPU ceiling")
        _count(cost["scalar_kernel_values_max"], "kernel ceiling")
        if cost["cpu_seconds"] > 120 or cost["cpu_seconds_max"] != 120 or cost["scalar_kernel_values_max"] != 30_000_000:
            raise ValueError("central ray budget changed or exceeded")
        for key in ("scalar_kernel_values", "metric_values", "source_quadrature_nodes"):
            _count(cost[key], key, 30_000_000)
        if cost["metric_values"]+cost["source_quadrature_nodes"] > 30_000_000:
            raise ValueError("central ray metric budget exceeded")
        _count(cost["field_rhs_calls"], "field RHS calls", 6000)
        for key, value in (("leading_field_solves", 1), ("wall_solves", 0), ("cap_clock_solves", 0),
                           ("scalar_propagations", 0), ("incoming_ray_evaluations", 5), ("reused_incoming_ray_evaluations", 5)):
            if type(cost[key]) is not int or cost[key] != value:
                raise ValueError("central ray reuse accounting mismatch")
        for key, samplekey in (("incoming_ray_nodes", "ray_nodes"), ("optical_radial_nodes", "optical_radial_nodes")):
            _count(cost[key], key, 10_000_000)
            if cost[key] != sum(row["integrals"][samplekey] for row in fine["samples"]):
                raise ValueError("central ray new-node accounting mismatch")
    except (KeyError, TypeError, IndexError) as error:
        raise ValueError("incomplete central ray record") from error
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
        "outputs": {"ray.json": {"sha256": digest(raw), "bytes": len(raw)}}}
    directory.mkdir(parents=True, exist_ok=False)
    (directory/"ray.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return {"created": str(directory), "execution": record["execution"], "diagnostics": record["diagnostics"]}


def check(directory):
    directory = Path(directory)
    path = directory/"manifest.json"
    manifest = json.loads(path.read_bytes())
    if (path.is_symlink() or manifest.get("schema") != SCHEMA+"-manifest"
            or manifest.get("source_bindings") != source_bindings()):
        raise ValueError("central ray source binding or schema mismatch")
    if ({p.name for p in directory.iterdir()} != {"manifest.json", "ray.json"}
            or set(manifest.get("outputs", {})) != {"ray.json"}):
        raise ValueError("central ray inventory mismatch")
    path = directory/"ray.json"
    raw = path.read_bytes()
    if path.is_symlink() or manifest["outputs"]["ray.json"] != {"sha256": digest(raw), "bytes": len(raw)}:
        raise ValueError("central ray output mismatch")
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
