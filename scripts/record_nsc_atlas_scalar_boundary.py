#!/usr/bin/env python3
"""Create or authenticate a bounded cubic scalar boundary-input comparison.

Creation shares one budget across preparation, wall and both boundary queries.
Authentication checks saved declarations and bytes without solving equations.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

from record_nsc_atlas_cap_response import (
    SOURCES as CAP_SOURCES, BASE as CAP_BASE, DOMAIN_KEYS, PROFILE_COUNTS,
    TIMES, digest, encoded, runtime_versions,
)

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-scalar-boundary-comparison-v1"
STATUS = "FORMAL_CUBIC_SCALAR_BOUNDARY_REFINEMENT"
SCOPE = "formal cubic scalar boundary input; no evolved waveform, boundary interpolant or finite-amplitude remainder"
ORDER = "chi/epsilon=c0+alpha*c1+O(alpha^2); alpha=kappa4*epsilon^2/ell^2"
SOURCES = tuple(dict.fromkeys((*CAP_SOURCES,
    "docs/nsc-atlas-scalar-metric-response.md",
    "scripts/record_nsc_atlas_scalar_boundary.py",
    "src/recursive_horizons/nsc_atlas_scalar_boundary.py",
    "tests/test_nsc_atlas_scalar_boundary.py",
    "tests/test_nsc_atlas_scalar_metric_response.py",
    "tests/test_nsc_atlas_scalar_flux.py",
    "tests/test_record_nsc_atlas_scalar_boundary.py")))
BASE = {key: value for key, value in CAP_BASE.items() if key != "time_panels"}
VARIATIONS = {
    "field": ("field_nodes", 129, 257),
    "radial": ("radial_order", 8, 16),
    "ray": ("ray_order", 8, 16),
    "radial_fine": ("radial_order", 16, 32),
    "ray_fine": ("ray_order", 16, 32),
}
COMPARISON_BASE_OVERRIDES = {kind: {"field_nodes": 129}
                            for kind in ("radial_fine", "ray_fine")}
SOLUTION_SCHEMA = "nsc-atlas-scalar-boundary-v1"
SOLUTION_SCOPE = "formal cubic outgoing boundary input; no propagated scalar correction or finite-amplitude remainder"
QUANTITIES = ("J_minus", "J_minus_r", "c1_r", "backscatter_term", "displacement_term")
DEFECTS = {"coordinate_integral_defect", "potential_identity_defect", "wall_flux_defect", "wall_profile_defect"}
MARGINS = {"no_return_margin", "quiet_centre_time_margin", "initial_time_margin",
           "last_profile_margin", "first_profile_margin", "active_radius_lower"}
SAMPLE_KEYS = {"t", "V", "x", "beta_plus", "beta_minus", "v", "V0_s", "q_minus", "q_plus",
    "T_wall", "rho", "V1", "w", "p_wall", "J_minus", "J_minus_r", "positive_integral",
    "negative_integral", "absolute_integral", "K1_sample_min", "K1_sample_max",
    "source_sample_min", "source_sample_max", "ray_nodes", "ray_nodes_r", "c1_r",
    "backscatter_term", "displacement_term", "e", "displacement_at_fixed_V", "c0_rr",
    "q_minus_s", "d", *DEFECTS, *MARGINS}


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def differences(summaries):
    rows = [summary["boundary_samples"] for summary in summaries]
    result = {}
    for key in QUANTITIES:
        left, right = (np.asarray([row[key] for row in values]) for values in rows)
        error = float(np.max(abs(left-right)))
        scale = float(max(np.max(abs(left)), np.max(abs(right))))
        result[key] = {"maximum_absolute_difference": error,
                       "maximum_coefficient_magnitude": scale,
                       "fraction_of_maximum_magnitude": error/scale if scale else 0.}
    return result


def compare(kind):
    from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget, solve_prepared_field
    from recursive_horizons.nsc_atlas_wedge_field import WedgeField
    from recursive_horizons.nsc_atlas_wall_evolution import solve_wall_evolution
    from recursive_horizons.nsc_atlas_bulk_response import BulkBudget
    from recursive_horizons.nsc_atlas_characteristic_initialization import source_geometry
    from recursive_horizons.nsc_atlas_cap_field import CapLeadingField
    from recursive_horizons.nsc_atlas_scalar_boundary import ScalarBoundary

    varied, low, high = VARIATIONS[kind]
    base = {**BASE, **COMPARISON_BASE_OVERRIDES.get(kind, {})}
    shared = PreparedFieldBudget()
    metrics = BulkBudget(started=shared.started, external=shared)
    geometry = source_geometry()
    fields, walls, cap_fields, summaries, settings = {}, {}, {}, [], []
    for value in (low, high):
        config = {**base, varied: value}
        field_key = config["field_nodes"], config["geometry_knots"]
        if field_key not in fields:
            fields[field_key] = solve_prepared_field(*field_key, budget=shared)
        wall_settings = {key.removeprefix("wall_"): val for key, val in config.items()
                         if key.startswith("wall_")}
        wall_key = field_key, tuple(sorted(wall_settings.items()))
        if wall_key not in walls:
            walls[wall_key] = solve_wall_evolution(WedgeField(fields[field_key]), geometry,
                                                  budget=metrics, source_scale=1., **wall_settings)
        cap_key = field_key, config["derivative_order"], config["small_radius_order"], config["centre_radius"]
        if cap_key not in cap_fields:
            cap_fields[cap_key] = CapLeadingField(fields[field_key],
                derivative_order=config["derivative_order"], small_radius_order=config["small_radius_order"],
                centre_radius=config["centre_radius"])
        boundary = ScalarBoundary(cap_fields[cap_key], walls[wall_key], budget=metrics,
            radial_order=config["radial_order"], ray_order=config["ray_order"])
        summary = boundary.summary(TIMES)
        for key in ("rhs_calls", "metric_values", "source_nodes", "process_cpu_seconds",
                    "scalar_kernel_values", "profile_counts"):
            summary.pop(key, None)
        summaries.append(summary)
        settings.append(config)
    shared.check()
    metrics.check()
    return {"schema": SCHEMA, "comparison": kind, "varied_quantity": varied,
        "status": STATUS, "scope": SCOPE, "order": ORDER, "source_scale": 1., "runtime": runtime_versions(),
        "comparison_base": base, "settings": settings, "solutions": summaries,
        "boundary_differences": differences(summaries),
        "refinement_is_rigorous_error_enclosure": False,
        "profile_evaluation": {
            "counts": {key: sum(field.evaluation_counts[key] for field in cap_fields.values())
                       for key in PROFILE_COUNTS},
            "maximum_inverse_bracket_width": max(field.max_inverse_bracket_width for field in cap_fields.values())},
        "execution": {"cpu_seconds": time.process_time()-shared.started,
            "cpu_seconds_max": shared.cpu_seconds, "scalar_kernel_values": shared.kernel_values,
            "scalar_kernel_values_max": shared.kernel_values_max,
            "metric_values": metrics.metric_values, "source_quadrature_nodes": metrics.source_nodes,
            "rhs_calls": metrics.rhs_calls, "field_solves": len(fields),
            "wall_solves": len(walls), "boundary_evaluators": 2}}


def _real(value, label, nonnegative=False):
    if (isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or (nonnegative and value < 0)):
        raise ValueError("invalid scalar boundary "+label)


def _count(value, label, maximum=None):
    if type(value) is not int or value < 0 or (maximum is not None and value > maximum):
        raise ValueError("invalid scalar boundary "+label)


def _equal(left, right):
    return math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-13)


def _validate_solution(solution, setting):
    if (solution["schema"] != SOLUTION_SCHEMA or solution["scope"] != SOLUTION_SCOPE
            or any(type(solution[key]) is not int or solution[key] != setting[key]
                   for key in ("radial_order", "ray_order"))):
        raise ValueError("scalar boundary solution settings or scope mismatch")
    domain = solution["domain"]
    if not isinstance(domain, dict) or set(domain) != DOMAIN_KEYS:
        raise ValueError("scalar boundary domain inventory mismatch")
    for value in domain.values():
        _real(value, "domain")
    if not (domain["U0"] < domain["Ua"] < domain["Ub"] < 0
            and domain["Va"] < domain["Vb"] < 0 and domain["U0"]*domain["Va"] < 1
            and domain["Ta"] < domain["Tb"] < domain["qa"] < domain["qb"]
            and _equal(domain["xa"], .45) and _equal(domain["xb"], .55)
            and _equal(domain["Ta"], domain["qa"]-math.atanh(domain["xa"]))
            and _equal(domain["Tb"], domain["qb"]-math.atanh(domain["xb"]))
            and _equal(domain["qmin"], domain["Ta"]-math.atanh(domain["xa"]))):
        raise ValueError("scalar boundary domain exits the original no-return cap")
    defects, margins = solution["maximum_absolute_defects"], solution["minimum_domain_margins"]
    if set(defects) != DEFECTS or set(margins) != MARGINS:
        raise ValueError("scalar boundary diagnostic inventory mismatch")
    for value in defects.values():
        _real(value, "defect", True)
    for value in margins.values():
        _real(value, "domain margin")
        if value < -1e-13:
            raise ValueError("scalar boundary negative domain margin")
    rows = solution["boundary_samples"]
    if not isinstance(rows, list) or len(rows) != len(TIMES):
        raise ValueError("scalar boundary sample coordinates mismatch")
    for row, t in zip(rows, TIMES):
        if not isinstance(row, dict) or set(row) != SAMPLE_KEYS or row["t"] != t:
            raise ValueError("scalar boundary sample inventory or coordinates mismatch")
        for key, value in row.items():
            if key in ("ray_nodes", "ray_nodes_r"):
                _count(value, key, 17*setting["ray_order"])
            else:
                _real(value, "sample "+key)
        x, q, u = row["x"], row["q_plus"], row["q_minus"]
        if (not .45-1e-13 <= x <= .55+1e-13 or not domain["qa"]-1e-13 <= q <= domain["qb"]+1e-13
                or not _equal(row["V"], domain["Va"]+(domain["Vb"]-domain["Va"])*t)
                or not _equal(u, q-2*math.atanh(x)) or not _equal(row["T_wall"], (u+q)/2)
                or row["V0_s"] <= 0 or row["d"] <= 0 or row["q_minus_s"] <= 0):
            raise ValueError("scalar boundary sample exits the declared model")
        expected_margins = dict(no_return_margin=domain["qa"]-u,
            quiet_centre_time_margin=domain["qa"]-row["T_wall"],
            initial_time_margin=(u+domain["qa"])/2-domain["Ta"],
            last_profile_margin=domain["qb"]-q, first_profile_margin=u-domain["qmin"],
            active_radius_lower=math.tanh((domain["qa"]-u)/2))
        if any(not _equal(row[key], expected_margins[key]) or row[key] < -1e-13 for key in MARGINS):
            raise ValueError("scalar boundary saved domain margins disagree")
        if (row["positive_integral"] < 0 or row["negative_integral"] > 0
                or row["J_minus"] != row["positive_integral"]+row["negative_integral"]
                or row["absolute_integral"] != row["positive_integral"]-row["negative_integral"]
                or row["K1_sample_min"] > row["K1_sample_max"]
                or row["source_sample_min"] > row["source_sample_max"]
                or row["ray_nodes"] != row["ray_nodes_r"]
                or row["ray_nodes"] % setting["ray_order"]
                or (q > domain["qa"]+1e-13 and row["ray_nodes"] < setting["ray_order"])):
            raise ValueError("scalar boundary quadrature diagnostics disagree")
        e = row["rho"]-x*row["p_wall"]**2/12
        offset = e-row["v"]*row["V1"]/row["V0_s"]
        hessian = .8*row["w"]/(x**3*row["d"]**2)
        backscatter = 2*row["q_minus_s"]*row["J_minus"]/row["d"]
        flux = row["d"]*(row["c1_r"]+offset*hessian)/2-row["q_minus_s"]*row["J_minus"]
        identities = {"e": e, "displacement_at_fixed_V": offset, "c0_rr": hessian,
            "d": row["v"]-row["beta_plus"], "q_minus_s": (row["beta_minus"]-row["v"])/(1-x*x),
            "backscatter_term": backscatter, "displacement_term": -offset*hessian,
            "c1_r": row["backscatter_term"]+row["displacement_term"],
            "wall_flux_defect": flux, "coordinate_integral_defect": row["J_minus"]-row["J_minus_r"]}
        if any(not _equal(row[key], value) for key, value in identities.items()):
            raise ValueError("scalar boundary saved wall or coordinate identity disagrees")
    if (not _equal(rows[0]["q_plus"], domain["qa"]) or not _equal(rows[-1]["q_plus"], domain["qb"])
            or any(b["q_plus"] <= a["q_plus"] for a, b in zip(rows[:-1], rows[1:]))
            or any(abs(rows[0][key]) > 1e-24 for key in QUANTITIES)):
        raise ValueError("scalar boundary ordered cuts or quiet anchor lost")
    if any(defects[key] != max(abs(row[key]) for row in rows) for key in DEFECTS):
        raise ValueError("scalar boundary saved defect maxima disagree")
    if any(margins[key] != min(row[key] for row in rows) for key in MARGINS):
        raise ValueError("scalar boundary saved domain minima disagree")
    for key, expected in (("boundary_evaluations", len(TIMES)),
                          ("ray_nodes_q", sum(row["ray_nodes"] for row in rows)),
                          ("ray_nodes_r", sum(row["ray_nodes_r"] for row in rows))):
        if type(solution[key]) is not int or solution[key] != expected:
            raise ValueError("scalar boundary query count mismatch")


def validate_record(record, kind):
    """Authenticate declarations and saved algebra, without asserting a waveform."""
    try:
        if kind not in VARIATIONS or not isinstance(record, dict):
            raise ValueError("scalar boundary comparison declaration mismatch")
        expected = {"schema": SCHEMA, "comparison": kind, "status": STATUS, "scope": SCOPE, "order": ORDER}
        if any(record.get(key) != value for key, value in expected.items()):
            raise ValueError("scalar boundary comparison declaration mismatch")
        if (record.get("refinement_is_rigorous_error_enclosure") is not False
                or record.get("source_scale") != 1. or isinstance(record.get("source_scale"), bool)):
            raise ValueError("scalar boundary source or non-enclosure declaration mismatch")
        varied, low, high = VARIATIONS[kind]
        base = {**BASE, **COMPARISON_BASE_OVERRIDES.get(kind, {})}
        settings = [{**base, varied: value} for value in (low, high)]
        if (record["varied_quantity"] != varied or record["settings"] != settings
                or record["comparison_base"] != base):
            raise ValueError("scalar boundary comparison settings mismatch")
        for setting in record["settings"]:
            if any(type(setting[key]) is not int for key in BASE if key not in ("wall_tolerance", "centre_radius")):
                raise ValueError("scalar boundary integer setting mismatch")
            _real(setting["wall_tolerance"], "wall tolerance")
            _real(setting["centre_radius"], "centre radius")
        runtime = record["runtime"]
        if (not isinstance(runtime, dict) or set(runtime) != {"python", "numpy", "scipy"}
                or any(not isinstance(v, str) or not v or not v[0].isdigit() or "." not in v for v in runtime.values())):
            raise ValueError("scalar boundary runtime mismatch")
        solutions = record["solutions"]
        if not isinstance(solutions, list) or len(solutions) != 2:
            raise ValueError("scalar boundary comparison requires two solutions")
        for solution, setting in zip(solutions, settings):
            _validate_solution(solution, setting)
            if any(not math.isclose(solution["domain"][key], solutions[0]["domain"][key], rel_tol=0, abs_tol=1e-13)
                   for key in DOMAIN_KEYS):
                raise ValueError("scalar boundary domains differ beyond endpoint roundoff")
        if record["boundary_differences"] != differences(solutions):
            raise ValueError("scalar boundary saved differences disagree")
        profile = record["profile_evaluation"]
        if set(profile["counts"]) != PROFILE_COUNTS:
            raise ValueError("scalar boundary profile count inventory mismatch")
        for key, value in profile["counts"].items():
            _count(value, key)
        counts = profile["counts"]
        if (not counts["second_points"] <= counts["wall_points"] <= counts["profile_requested"]
                or counts["cache_hits"] > counts["profile_requested"]):
            raise ValueError("scalar boundary profile count mismatch")
        _real(profile["maximum_inverse_bracket_width"], "inverse bracket", True)
        domain = solutions[0]["domain"]
        if profile["maximum_inverse_bracket_width"] > 4*np.finfo(float).eps*max(abs(domain["Va"]), abs(domain["Vb"])):
            raise ValueError("scalar boundary inverse bracket exceeds its numerical tolerance")
        execution = record["execution"]
        _real(execution["cpu_seconds_max"], "CPU budget")
        _count(execution["scalar_kernel_values_max"], "kernel budget")
        if execution["cpu_seconds_max"] != 120 or execution["scalar_kernel_values_max"] != 30_000_000:
            raise ValueError("scalar boundary budget ceilings changed")
        _real(execution["cpu_seconds"], "CPU time", True)
        if execution["cpu_seconds"] > 120:
            raise ValueError("scalar boundary CPU budget exceeded")
        for key in ("scalar_kernel_values", "metric_values", "source_quadrature_nodes"):
            _count(execution[key], key, 30_000_000)
        if execution["metric_values"]+execution["source_quadrature_nodes"] > 30_000_000:
            raise ValueError("scalar boundary metric budget exceeded")
        _count(execution["rhs_calls"], "RHS calls", 12000)
        for key, expected in (("boundary_evaluators", 2), ("wall_solves", 2 if kind == "field" else 1),
                              ("field_solves", 2 if kind == "field" else 1)):
            if type(execution[key]) is not int or execution[key] != expected:
                raise ValueError("scalar boundary solve count mismatch")
    except (KeyError, TypeError, IndexError) as error:
        raise ValueError("incomplete scalar boundary comparison record") from error
    return True


def create(directory, kind):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only record: choose a fresh directory or --check")
    if kind not in VARIATIONS:
        raise ValueError("scalar boundary comparison is not admitted")
    bindings, runtime = source_bindings(), runtime_versions()
    record = compare(kind)
    if bindings != source_bindings() or runtime != runtime_versions() or record["runtime"] != runtime:
        raise RuntimeError("calculation inputs or runtime changed during comparison")
    validate_record(record, kind)
    raw = encoded(record)
    manifest = {"schema": SCHEMA+"-manifest", "source_bindings": bindings,
        "comparison": kind, "outputs": {"boundary.json": {"sha256": digest(raw), "bytes": len(raw)}}}
    directory.mkdir(parents=True, exist_ok=False)
    (directory/"boundary.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return {"created": str(directory), "comparison": kind,
            "execution": record["execution"], "differences": record["boundary_differences"]}


def check(directory):
    from validation_source_successors import verify_saved_bindings
    directory = Path(directory)
    manifest_path = directory/"manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    if (manifest_path.is_symlink() or manifest.get("schema") != SCHEMA+"-manifest"
            or manifest.get("comparison") not in VARIATIONS):
        raise ValueError("scalar boundary source binding or schema mismatch")
    verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=SOURCES)
    if ({p.name for p in directory.iterdir()} != {"manifest.json", "boundary.json"}
            or set(manifest.get("outputs", {})) != {"boundary.json"}):
        raise ValueError("scalar boundary inventory mismatch")
    path = directory/"boundary.json"
    raw = path.read_bytes()
    if path.is_symlink() or manifest["outputs"]["boundary.json"] != {"sha256": digest(raw), "bytes": len(raw)}:
        raise ValueError("scalar boundary output mismatch")
    validate_record(json.loads(raw), manifest["comparison"])
    context="current" if manifest["source_bindings"] == source_bindings() else "reviewed-successors"
    return {"checked": str(directory), "comparison": manifest["comparison"],
            "source_context": context, "physical_solves": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", choices=VARIATIONS, default="field")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    print(json.dumps(check(args.output) if args.check else create(args.output, args.comparison), sort_keys=True))


if __name__ == "__main__":
    main()
