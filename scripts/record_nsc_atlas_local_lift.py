#!/usr/bin/env python3
"""Compare local boundary lifts on 16/32 active cells with shared heldouts.

The same 81 direct moment jets serve both fits. Authentication reconstructs
only polynomials from saved jets, never a physical field or metric solution.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np
from numpy.polynomial import Polynomial
from validation_source_successors import verify_saved_bindings

import record_nsc_atlas_profile_moment as predecessor
from record_nsc_atlas_scalar_inputs import (
    BASE, JET_KEYS, JET_DEFECTS, MARGINS, PROFILE_COUNTS, _prepare, _plain,
    _real, _count, _equal, interpolation_target, digest, encoded, runtime_versions,
)

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-local-lift-comparison-v1"
STATUS = "FORMAL_LOCAL_BOUNDARY_LIFT_REFINEMENT"
SCOPE = "local polynomial boundary-input refinement; no evolved scalar correction, central waveform or uniform remainder"
INPUT = "results/development/nsc-atlas-profile-moment-v1"
ACTIVE_CELLS = (16, 32)
SETTINGS = {**BASE, "moment_order": 8, "outer_cells": 4, "stencil_nodes": 5}
SOURCES = tuple(dict.fromkeys((*predecessor.SOURCES,
    INPUT+"/moment.json", INPUT+"/manifest.json",
    "src/recursive_horizons/nsc_atlas_adaptive_lift.py", "tests/test_nsc_atlas_adaptive_lift.py",
    "scripts/record_nsc_atlas_local_lift.py", "tests/test_record_nsc_atlas_local_lift.py")))


def knot_times(active):
    return np.unique(np.r_[np.linspace(0., .25, 5), np.linspace(.25, .75, active+1), np.linspace(.75, 1., 5)])


FINE_TIMES = knot_times(32)
WITHHELD = tuple(float(t) for t in (FINE_TIMES[:-1]+FINE_TIMES[1:])/2)
DIRECT_TIMES = tuple(sorted(set(FINE_TIMES) | set(WITHHELD)))


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def input_snapshot():
    path = ROOT/INPUT
    predecessor.check(path)
    saved = json.loads((path/"moment.json").read_bytes())
    return {"path": INPUT, "files": {name: digest((path/name).read_bytes()) for name in ("moment.json", "manifest.json")},
            "domain": saved["domain"], "physical_solves": 0}


def solution_diagnostics(solution):
    keys = ("b_defect", "scaled_slope_defect", "lift_derivative_defect", "wall_transport_defect")
    rows, cells = solution["validation_samples"], solution["cells"]
    per_cell = [{"cell": cell["index"], "block": cell["block"],
        "withheld_samples": sum(row["cell"] == cell["index"] for row in rows),
        **{key: max(abs(row[key]) for row in rows if row["cell"] == cell["index"]) for key in keys}} for cell in cells]
    per_block = [{"block": block, **{key: max(row[key] for row in per_cell if row["block"] == block) for key in keys}}
                 for block in range(3)]
    value_jumps, slope_jumps = [], []
    for left, right in zip(cells[:-1], cells[1:]):
        p, q = Polynomial(left["coefficients_in_cell_coordinate"]), Polynomial(right["coefficients_in_cell_coordinate"])
        value_jumps.append(abs(float(p(1.)-q(0.))))
        slope_jumps.append(min(left["width"], right["width"])*abs(float(p.deriv()(1.)/left["width"]-q.deriv()(0.)/right["width"])))
    return {"per_cell": per_cell, "per_block": per_block,
            "maximum_endpoint_value_jump": max(value_jumps),
            "maximum_scaled_endpoint_slope_jump": max(slope_jumps)}


def comparison_diagnostics(solutions):
    target = interpolation_target()
    left, right = solutions
    return {"interpolation_target": target,
        "sampled_b_target_met": [s["maximum_absolute_defects"]["b_defect"] <= target["absolute_b_target"] for s in solutions],
        "maximum_lift_difference": max(abs(left["terminal_lift"]-right["terminal_lift"]),
            *(abs(a["B_lift"]-b["B_lift"]) for a, b in zip(left["validation_samples"], right["validation_samples"]))),
        "solutions": [solution_diagnostics(solution) for solution in solutions]}


def compare():
    from recursive_horizons.nsc_atlas_scalar_boundary import ScalarBoundary
    from recursive_horizons.nsc_atlas_profile_moment import ProfileMoment, MomentBoundaryJet
    from recursive_horizons.nsc_atlas_adaptive_lift import AdaptiveBoundaryLift
    inputs = input_snapshot()
    shared, metrics, _, _, wall, cap = _prepare()
    boundary = ScalarBoundary(cap, wall, radial_order=16, ray_order=16, budget=metrics)
    moment = ProfileMoment(cap, order=8, budget=metrics)
    jets = MomentBoundaryJet(boundary, moment_order=8, moment=moment)
    fits = [AdaptiveBoundaryLift(jets, knot_times=knot_times(active), nodes=5, validation_times=WITHHELD)
            for active in ACTIVE_CELLS]
    solutions = [fit.summary() for fit in fits]
    direct = [dict(jets.cache[t]) for t in sorted(jets.cache)]
    shared.check()
    metrics.check()
    return _plain({"schema": SCHEMA, "status": STATUS, "scope": SCOPE, "runtime": runtime_versions(),
        "input": inputs, "comparison_base": SETTINGS, "varied_quantity": "active_cells",
        "settings": [{**SETTINGS, "active_cells": count} for count in ACTIVE_CELLS], "source_scale": 1.,
        "domain": dict(cap.domain), "direct_samples": direct, "solutions": solutions,
        "diagnostics": comparison_diagnostics(solutions), "refinement_is_rigorous_error_enclosure": False,
        "unique_boundary_jet_evaluations": len(direct), "ray_nodes": jets.ray_nodes, "moment": moment.summary(),
        "maximum_absolute_jet_defects": {key: max(abs(row[key]) for row in direct) for key in JET_DEFECTS},
        "minimum_domain_margins": {key: min(row[key] for row in direct) for key in MARGINS},
        "profile_evaluation": {"counts": dict(cap.evaluation_counts), "maximum_inverse_bracket_width": cap.max_inverse_bracket_width},
        "execution": {"cpu_seconds": time.process_time()-shared.started, "cpu_seconds_max": shared.cpu_seconds,
            "scalar_kernel_values": shared.kernel_values, "scalar_kernel_values_max": shared.kernel_values_max,
            "metric_values": metrics.metric_values, "source_quadrature_nodes": metrics.source_nodes,
            "rhs_calls": metrics.rhs_calls, "field_solves": 1, "wall_solves": 1, "moment_evaluators": 1}})


def _same_saved(actual, expected):
    """Finite saved-data consistency, allowing only floating polynomial roundoff."""
    if isinstance(expected, dict):
        return isinstance(actual, dict) and set(actual) == set(expected) and all(_same_saved(actual[k], v) for k, v in expected.items())
    if isinstance(expected, list):
        return isinstance(actual, list) and len(actual) == len(expected) and all(_same_saved(a, b) for a, b in zip(actual, expected))
    if isinstance(expected, bool):
        return type(actual) is bool and actual == expected
    if isinstance(expected, int):
        return type(actual) is int and actual == expected
    if isinstance(expected, float):
        return type(actual) in (int, float) and math.isfinite(actual) and math.isfinite(expected) and _equal(actual, expected)
    return actual == expected


def _saved_fits(record):
    """Rebuild only declared polynomial fits from the authenticated direct jets."""
    from recursive_horizons.nsc_atlas_adaptive_lift import AdaptiveBoundaryLift
    direct = record["direct_samples"]
    by_time = {row["t"]: row for row in direct}
    by_V = {row["V"]: row["U"] for row in direct}
    def wall(values):
        values = np.asarray(values, dtype=float)
        return np.array([by_V[float(v)] for v in values.ravel()]).reshape(values.shape)
    jets = SimpleNamespace(domain=record["domain"], budget=SimpleNamespace(check=lambda: None),
        geometry={"wall": wall}, cache=by_time, sample=lambda t: dict(by_time[t]))
    return [_plain(AdaptiveBoundaryLift(jets, knot_times=knot_times(count), nodes=5, validation_times=WITHHELD).summary())
            for count in ACTIVE_CELLS]


def validate_record(record):
    try:
        if any(record.get(key) != value for key, value in {"schema": SCHEMA, "status": STATUS, "scope": SCOPE}.items()):
            raise ValueError("local lift scope or schema mismatch")
        if (record.get("refinement_is_rigorous_error_enclosure") is not False
                or record.get("source_scale") != 1. or isinstance(record.get("source_scale"), bool)):
            raise ValueError("local lift source or non-enclosure mismatch")
        if (record["comparison_base"] != SETTINGS or record["varied_quantity"] != "active_cells"
                or record["settings"] != [{**SETTINGS, "active_cells": count} for count in ACTIVE_CELLS]):
            raise ValueError("local lift comparison settings mismatch")
        for setting in record["settings"]:
            for key, value in setting.items():
                if key in ("wall_tolerance", "centre_radius"):
                    _real(value, key)
                elif type(value) is not int:
                    raise ValueError("local lift integer setting mismatch")
        runtime = record["runtime"]
        if (set(runtime) != {"python", "numpy", "scipy"}
                or any(not isinstance(v, str) or not v or not v[0].isdigit() or "." not in v for v in runtime.values())):
            raise ValueError("local lift runtime mismatch")
        inputs = input_snapshot()
        if record["input"] != inputs or record["domain"] != inputs["domain"]:
            raise ValueError("local lift predecessor identity or domain mismatch")
        domain, direct = record["domain"], record["direct_samples"]
        if (not isinstance(direct, list) or len(direct) != 81 or [row["t"] for row in direct] != list(DIRECT_TIMES)
                or type(record["unique_boundary_jet_evaluations"]) is not int or record["unique_boundary_jet_evaluations"] != 81):
            raise ValueError("local lift direct-query inventory mismatch")
        for row in direct:
            if not set(predecessor.SAMPLE_KEYS) <= set(row):
                raise ValueError("local lift required direct jets are missing")
            for key, value in row.items():
                if key == "boundary_cache_checked":
                    if type(value) is not bool:
                        raise ValueError("local lift cache status must be boolean")
                elif key == "ray_nodes":
                    _count(value, key, 17*SETTINGS["ray_order"])
                    if value % SETTINGS["ray_order"]:
                        raise ValueError("local lift ray accounting mismatch")
                else:
                    _real(value, "jet "+key)
            if (not domain["Ua"]-1e-13 <= row["U"] <= domain["Ub"]+1e-13
                    or not _equal(row["V"], domain["Va"]+(domain["Vb"]-domain["Va"])*row["t"])
                    or row["U_V"] <= 0 or row["r_U"] <= 0 or not .45 <= row["x"] <= .55
                    or not domain["qa"]-1e-13 <= row["q_plus"] <= domain["qb"]+1e-13
                    or not row["q_minus"] < domain["qa"] or any(row[key] < -1e-13 for key in MARGINS)
                    or row["positive_integral"] < 0 or row["negative_integral"] > 0):
                raise ValueError("local lift direct jets exit the declared domain")
            expected = {"b": row["r_U"]*row["c1_r"], "b_V": row["b_U"]*row["U_V"],
                "q_minus": row["q_plus"]-2*math.atanh(row["x"]),
                "J_minus": row["positive_integral"]+row["negative_integral"], "J_minus_u": -row["integral_Ku_f"],
                "J_minus_s": -row["q_plus_s"]*row["K1_wall"]*row["wall_profile"]-row["q_minus_s"]*row["integral_Ku_f"]}
            if any(not _equal(row[key], value) for key, value in expected.items()):
                raise ValueError("local lift saved derivative identity disagrees")
        if (any(b["U"] <= a["U"] for a, b in zip(direct[:-1], direct[1:]))
                or not _equal(direct[0]["U"], domain["Ua"]) or not _equal(direct[-1]["U"], domain["Ub"])
                or any(abs(direct[0][key]) > 1e-24 for key in ("b", "b_U", "b_V", "J_minus", "J_minus_s"))):
            raise ValueError("local lift ordering or quiet anchor lost")
        _count(record["ray_nodes"], "total ray nodes")
        if record["ray_nodes"] != sum(row["ray_nodes"] for row in direct):
            raise ValueError("local lift cumulative ray count mismatch")
        for key, expected in (("maximum_absolute_jet_defects", {k: max(abs(row[k]) for row in direct) for k in JET_DEFECTS}),
                              ("minimum_domain_margins", {k: min(row[k] for row in direct) for k in MARGINS})):
            if record[key] != expected:
                raise ValueError("local lift saved direct diagnostics disagree")
        rebuilt = _saved_fits(record)
        if not _same_saved(record["solutions"], rebuilt):
            raise ValueError("local lift saved polynomial, stencil or heldout result disagrees")
        if not _same_saved(record["diagnostics"], comparison_diagnostics(rebuilt)):
            raise ValueError("local lift saved comparison or target disagrees")
        moment = record["moment"]
        count_keys = ("cells", "prefix_nodes", "partial_nodes", "profile_points", "profile_batches",
                      "value_queries", "cache_hits", "cached_values", "mass_points")
        if (set(moment) != {"schema", "domain", "settings", "terminal_C", *count_keys}
                or moment["schema"] != "nsc-atlas-profile-moment-v1" or moment["domain"] != domain
                or moment["settings"] != {"moment_order": 8} or type(moment["settings"]["moment_order"]) is not int):
            raise ValueError("local lift moment declaration mismatch")
        for key in count_keys:
            _count(moment[key], key, 30_000_000)
        _real(moment["terminal_C"], "terminal moment", True)
        if (not 1 <= moment["cells"] <= 1024 or moment["prefix_nodes"] != 8*moment["cells"]
                or moment["partial_nodes"] % 8 or moment["partial_nodes"] > 8*moment["value_queries"]
                or moment["mass_points"] != record["ray_nodes"]+81 or moment["value_queries"] != moment["mass_points"]
                or moment["cached_values"] > min(32768, moment["value_queries"]) or moment["cache_hits"] > moment["value_queries"]
                or moment["profile_points"] > moment["prefix_nodes"]+moment["partial_nodes"]+moment["mass_points"]
                or not (moment["profile_points"]+512)//513 <= moment["profile_batches"] <= moment["profile_points"]):
            raise ValueError("local lift moment accounting mismatch")
        profile = record["profile_evaluation"]
        if set(profile["counts"]) != PROFILE_COUNTS:
            raise ValueError("local lift profile count inventory mismatch")
        for key, value in profile["counts"].items():
            _count(value, key)
        counts = profile["counts"]
        if counts["second_points"] != 0 or counts["wall_points"] > counts["profile_requested"] or counts["cache_hits"] > counts["profile_requested"]:
            raise ValueError("local lift second derivatives or profile counts mismatch")
        _real(profile["maximum_inverse_bracket_width"], "inverse bracket", True)
        if profile["maximum_inverse_bracket_width"] > 4*np.finfo(float).eps*max(abs(domain["Va"]), abs(domain["Vb"])):
            raise ValueError("local lift inverse bracket exceeds its numerical tolerance")
        cost = record["execution"]
        _real(cost["cpu_seconds_max"], "CPU ceiling")
        _count(cost["scalar_kernel_values_max"], "kernel ceiling")
        if cost["cpu_seconds_max"] != 120 or cost["scalar_kernel_values_max"] != 30_000_000:
            raise ValueError("local lift budget ceilings changed")
        _real(cost["cpu_seconds"], "CPU time", True)
        if cost["cpu_seconds"] > 120:
            raise ValueError("local lift CPU budget exceeded")
        for key in ("scalar_kernel_values", "metric_values", "source_quadrature_nodes"):
            _count(cost[key], key, 30_000_000)
        if cost["metric_values"]+cost["source_quadrature_nodes"] > 30_000_000:
            raise ValueError("local lift metric budget exceeded")
        _count(cost["rhs_calls"], "RHS calls", 12000)
        if any(type(cost[key]) is not int or cost[key] != 1 for key in ("field_solves", "wall_solves", "moment_evaluators")):
            raise ValueError("local lift evaluator count mismatch")
    except (KeyError, TypeError, IndexError) as error:
        raise ValueError("incomplete local lift record") from error
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
        "outputs": {"lift.json": {"sha256": digest(raw), "bytes": len(raw)}}}
    directory.mkdir(parents=True, exist_ok=False)
    (directory/"lift.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return {"created": str(directory), "execution": record["execution"],
        "withheld_maxima": [solution["maximum_absolute_defects"] for solution in record["solutions"]],
        "sampled_b_target_met": record["diagnostics"]["sampled_b_target_met"],
        "maximum_lift_difference": record["diagnostics"]["maximum_lift_difference"]}


def check(directory):
    directory = Path(directory)
    manifest_path = directory/"manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    if (manifest_path.is_symlink() or manifest.get("schema") != SCHEMA+"-manifest"
            or not verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=SOURCES)):
        raise ValueError("local lift source binding or schema mismatch")
    if ({p.name for p in directory.iterdir()} != {"manifest.json", "lift.json"}
            or set(manifest.get("outputs", {})) != {"lift.json"}):
        raise ValueError("local lift inventory mismatch")
    path = directory/"lift.json"
    raw = path.read_bytes()
    if path.is_symlink() or manifest["outputs"]["lift.json"] != {"sha256": digest(raw), "bytes": len(raw)}:
        raise ValueError("local lift output mismatch")
    validate_record(json.loads(raw))
    return {"checked": str(directory), "physical_solves": 0, "saved_polynomial_reconstructions": 2}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    print(json.dumps(check(args.output) if args.check else create(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
