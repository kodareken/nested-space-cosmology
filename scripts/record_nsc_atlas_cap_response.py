#!/usr/bin/env python3
"""Create or authenticate a finite cap-metric and central-null-clock comparison.

Creation uses one cumulative budget for preparation, wall, cap and diagnostics.
Authentication checks declarations and bytes without solving those equations.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

from record_nsc_atlas_wall_evolution import (
    SOURCES as WALL_SOURCES, digest, encoded, runtime_versions,
)

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-cap-response-comparison-v1"
STATUS = "FORMAL_CAP_METRIC_NULL_CLOCK_REFINEMENT"
SCOPE = "quadratic cap metric and central null clock; no higher-order scalar amplitude or finite-amplitude enclosure"
ORDER = "T_C=q0+alpha*Delta_C+O(alpha^2); alpha=kappa4*epsilon^2/ell^2"
SOLUTION_SCHEMA = "nsc-atlas-cap-response-v1"
SOLUTION_SCOPE = "formal cap metric and central null-arrival coefficients; no corrected scalar amplitude or finite-amplitude remainder"
SOURCES = tuple(dict.fromkeys((*WALL_SOURCES,
    "docs/nsc-atlas-cap-response.md", "scripts/record_nsc_atlas_cap_response.py",
    "src/recursive_horizons/nsc_atlas_cap_field.py",
    "src/recursive_horizons/nsc_atlas_cap_response.py",
    "tests/test_nsc_atlas_cap_field.py", "tests/test_nsc_atlas_cap_response.py",
    "tests/test_record_nsc_atlas_cap_response.py")))
# The primary field comparison uses the finer common cap quadratures. Selected
# auxiliary comparisons hold a 129-node field fixed to respect the same budget.
BASE = dict(field_nodes=257, geometry_knots=257, wall_radial_order=2,
            wall_field_order=32, wall_time_panels=8, wall_source_order=64,
            wall_tolerance=2e-9, derivative_order=2, small_radius_order=8,
            centre_radius=1e-4, radial_order=16, ray_order=16, time_panels=16)
VARIATIONS = {
    "field": ("field_nodes", 129, 257),
    "radial": ("radial_order", 8, 16),
    "ray": ("ray_order", 8, 16),
    "time": ("time_panels", 8, 16),
    "geometry": ("geometry_knots", 129, 257),
    "derivative": ("derivative_order", 2, 4),
    "centre_order": ("small_radius_order", 8, 16),
    "centre_radius": ("centre_radius", 1e-4, 1e-3),
    "ray_fine": ("ray_order", 16, 32),
}
COMPARISON_BASE_OVERRIDES = {kind: {"field_nodes": 129} for kind in
                            ("geometry", "derivative", "centre_order", "centre_radius", "ray_fine")}
TIMES = tuple(float(t) for t in np.linspace(0, 1, 9))
QUANTITIES = ("nu_E_wall", "nu_J_wall", "delta_E_wall", "T1", "T1_s", "ray_integral",
              "Delta_C_fixed_s", "Delta_C_fixed_V")
DOMAIN_KEYS = {"U0", "Ua", "Ub", "Va", "Vb", "Ta", "Tb", "qa", "qb", "qmin", "xa", "xb"}
DEFECTS = {"wall_einstein_mass_defect", "wall_jordan_mass_defect",
           "wall_frame_mass_defect", "wall_scalar_value_defect", "wall_scalar_tangent_defect",
           "wall_scalar_normal_defect", "reference_arrival_rate_defect", "wall_normalization_defect",
           "central_clock_normalization_defect", "fixed_V_conversion_defect"}
SAMPLE_KEYS = {"t", "V", "x", "T0_wall", "q0", "nu_E_wall", "nu_J_wall", "delta_E_wall",
               "nu_E_potential_wall", "integral_L_minus", "T1", "T1_s", "e", "e_s", "A1",
               "q0_s", "V0_s", "ray_integral", "Delta_C_fixed_s", "Delta_C_fixed_V",
               "exterior_label_shift", "V1", "wall_einstein_mass_defect",
               "wall_jordan_mass_defect", "wall_normalization_defect",
               "ray_positive_integral", "ray_negative_integral",
               "ray_integrand_sample_min", "ray_integrand_sample_max"}
PROFILE_COUNTS = {"profile_requested", "wall_points", "second_points", "cache_hits"}


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def differences(summaries):
    """Compare the same coefficients at the same reference wall cuts."""
    rows = [summary["cap_samples"] for summary in summaries]
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
    from recursive_horizons.nsc_atlas_cap_response import solve_cap_response

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
        cap_field = cap_fields[cap_key]
        cap = solve_cap_response(cap_field, walls[wall_key], budget=metrics,
            **{key: config[key] for key in ("radial_order", "ray_order", "time_panels")})
        summary = cap.summary(TIMES)
        # Costs belong to the whole comparison; the shared counters are cumulative.
        for key in ("rhs_calls", "metric_values", "source_nodes", "process_cpu_seconds"):
            summary.pop(key, None)
        summaries.append(summary)
        settings.append(config)
    shared.check()
    metrics.check()
    return {"schema": SCHEMA, "comparison": kind, "varied_quantity": varied,
        "status": STATUS, "scope": SCOPE, "order": ORDER, "source_scale": 1., "runtime": runtime_versions(),
        "comparison_base": base, "settings": settings, "solutions": summaries,
        "cap_differences": differences(summaries),
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
            "wall_solves": len(walls), "cap_solves": 2}}


def validate_record(record, kind):
    """Check the declared finite experiment, without verifying an unsaved trajectory."""
    def real(value, label, nonnegative=False):
        if (isinstance(value, bool) or not isinstance(value, (int, float))
                or not math.isfinite(value) or (nonnegative and value < 0)):
            raise ValueError("invalid cap comparison "+label)

    def count(value, label, maximum=None):
        if type(value) is not int or value < 0 or (maximum is not None and value > maximum):
            raise ValueError("invalid cap comparison "+label)

    def equal(left, right):
        return math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-13)

    try:
        if kind not in VARIATIONS or not isinstance(record, dict):
            raise ValueError("cap comparison declaration mismatch")
        expected = {"schema": SCHEMA, "comparison": kind, "status": STATUS, "scope": SCOPE, "order": ORDER}
        if any(record.get(key) != value for key, value in expected.items()):
            raise ValueError("cap comparison declaration mismatch")
        if (record.get("refinement_is_rigorous_error_enclosure") is not False
                or record.get("source_scale") != 1. or isinstance(record.get("source_scale"), bool)):
            raise ValueError("cap comparison source or non-enclosure declaration mismatch")
        varied, low, high = VARIATIONS[kind]
        base = {**BASE, **COMPARISON_BASE_OVERRIDES.get(kind, {})}
        settings = [{**base, varied: value} for value in (low, high)]
        if (record["varied_quantity"] != varied or record["settings"] != settings
                or record["comparison_base"] != base):
            raise ValueError("cap comparison settings mismatch")
        for setting in record["settings"]:
            if any(type(setting[key]) is not int for key in BASE if key not in ("wall_tolerance", "centre_radius")):
                raise ValueError("cap comparison integer setting mismatch")
            real(setting["wall_tolerance"], "wall tolerance")
            real(setting["centre_radius"], "centre radius")
        runtime = record["runtime"]
        if (not isinstance(runtime, dict) or set(runtime) != {"python", "numpy", "scipy"}
                or any(not isinstance(v, str) or not v or not v[0].isdigit() or "." not in v for v in runtime.values())):
            raise ValueError("cap comparison runtime mismatch")
        solutions = record["solutions"]
        if not isinstance(solutions, list) or len(solutions) != 2:
            raise ValueError("cap comparison requires two solutions")
        for solution, setting in zip(solutions, settings):
            if (solution["schema"] != SOLUTION_SCHEMA or solution["scope"] != SOLUTION_SCOPE
                    or type(solution["clock_time_order"]) is not int or solution["clock_time_order"] != 4
                    or any(type(solution[key]) is not int or solution[key] != setting[key]
                           for key in ("radial_order", "ray_order", "time_panels"))):
                raise ValueError("cap comparison solution settings or scope mismatch")
            domain = solution["domain"]
            if not isinstance(domain, dict) or set(domain) != DOMAIN_KEYS or domain != solutions[0]["domain"]:
                raise ValueError("cap comparison domains differ")
            for value in domain.values():
                real(value, "domain")
            if not (domain["U0"] < domain["Ua"] < domain["Ub"] < 0
                    and domain["Va"] < domain["Vb"] < 0 and domain["U0"]*domain["Va"] < 1
                    and domain["Ta"] < domain["Tb"] < domain["qb"] and domain["qa"] < domain["qb"]
                    and equal(domain["xa"], .45) and equal(domain["xb"], .55)
                    and equal(domain["Ta"], domain["qa"]-math.atanh(domain["xa"]))
                    and equal(domain["Tb"], domain["qb"]-math.atanh(domain["xb"]))
                    and equal(domain["qmin"], domain["Ta"]-math.atanh(domain["xa"]))):
                raise ValueError("cap comparison domain exits the original finite cap")
            defects = solution["maximum_absolute_defects"]
            if not isinstance(defects, dict) or set(defects) != DEFECTS:
                raise ValueError("cap comparison defect inventory mismatch")
            for value in defects.values():
                real(value, "defect", True)
            rows = solution["cap_samples"]
            if not isinstance(rows, list) or len(rows) != len(TIMES):
                raise ValueError("cap comparison sample coordinates mismatch")
            for row, t in zip(rows, TIMES):
                if not isinstance(row, dict) or set(row) != SAMPLE_KEYS or row["t"] != t:
                    raise ValueError("cap comparison sample inventory or coordinates mismatch")
                for key, value in row.items():
                    real(value, "sample "+key)
                if (not .45-1e-13 <= row["x"] <= .55+1e-13
                        or not domain["qa"]-1e-13 <= row["q0"] <= domain["qb"]+1e-13
                        or not equal(row["V"], domain["Va"]+(domain["Vb"]-domain["Va"])*t)
                        or not equal(row["q0"], row["T0_wall"]+math.atanh(row["x"]))
                        or row["V0_s"] <= 0 or row["q0_s"] <= 0
                        or not 0 <= row["nu_E_potential_wall"] <= row["nu_E_wall"]):
                    raise ValueError("cap comparison samples leave the declared model or clock")
                if (not equal(row["wall_einstein_mass_defect"], row["nu_E_wall"]-row["integral_L_minus"])
                        or not equal(row["Delta_C_fixed_s"], row["T1"]+row["e"]/(1-row["x"]**2)+row["ray_integral"])
                        or not equal(row["exterior_label_shift"], row["q0_s"]*row["V1"]/row["V0_s"])
                        or not equal(row["Delta_C_fixed_V"], row["Delta_C_fixed_s"]-row["exterior_label_shift"])):
                    raise ValueError("cap comparison saved mass or clock conversion disagrees")
                if (row["ray_positive_integral"] < 0 or row["ray_negative_integral"] > 0
                        or row["ray_integrand_sample_min"] > row["ray_integrand_sample_max"]
                        or row["ray_integral"] != row["ray_positive_integral"]+row["ray_negative_integral"]):
                    raise ValueError("cap comparison ray quadrature diagnostics disagree")
                if any(abs(row[key]) > defects[key] for key in
                       ("wall_einstein_mass_defect", "wall_jordan_mass_defect", "wall_normalization_defect")):
                    raise ValueError("cap comparison displayed defect exceeds saved maximum")
            if (not equal(rows[0]["q0"], domain["qa"]) or not equal(rows[-1]["q0"], domain["qb"])
                    or any(b["q0"] <= a["q0"] for a, b in zip(rows[:-1], rows[1:]))):
                raise ValueError("cap comparison arrival cuts are not the ordered finite interval")
            if any(abs(rows[0][key]) > 1e-24 for key in QUANTITIES):
                raise ValueError("cap comparison lost the quiet clock anchor")
            for key, expected in (("wall_metric_evaluations", 4*setting["time_panels"]+len(TIMES)),
                                  ("ray_evaluations", len(TIMES))):
                if type(solution[key]) is not int or solution[key] != expected:
                    raise ValueError("cap comparison query count mismatch")
            count(solution["ray_nodes"], "ray nodes", 9*len(TIMES)*setting["ray_order"])
            if solution["ray_nodes"] < len(TIMES)*setting["ray_order"]:
                raise ValueError("cap comparison ray quadrature inventory mismatch")
        if record["cap_differences"] != differences(solutions):
            raise ValueError("cap comparison saved differences disagree")
        profile = record["profile_evaluation"]
        if set(profile["counts"]) != PROFILE_COUNTS:
            raise ValueError("cap comparison profile count inventory mismatch")
        for key, value in profile["counts"].items():
            count(value, key)
        counts = profile["counts"]
        if (not counts["second_points"] <= counts["wall_points"] <= counts["profile_requested"]
                or counts["cache_hits"] > counts["profile_requested"]):
            raise ValueError("cap comparison profile count mismatch")
        real(profile["maximum_inverse_bracket_width"], "inverse bracket", True)
        if profile["maximum_inverse_bracket_width"] > 4*np.finfo(float).eps*max(abs(domain["Va"]), abs(domain["Vb"])):
            raise ValueError("cap comparison inverse bracket exceeds its declared numerical tolerance")
        execution = record["execution"]
        if execution["cpu_seconds_max"] != 120 or execution["scalar_kernel_values_max"] != 30_000_000:
            raise ValueError("cap comparison budget ceilings changed")
        real(execution["cpu_seconds"], "CPU time", True)
        if execution["cpu_seconds"] > 120:
            raise ValueError("cap comparison CPU budget exceeded")
        for key in ("scalar_kernel_values", "metric_values", "source_quadrature_nodes"):
            count(execution[key], key, 30_000_000)
        if execution["metric_values"]+execution["source_quadrature_nodes"] > 30_000_000:
            raise ValueError("cap comparison metric budget exceeded")
        count(execution["rhs_calls"], "RHS calls", 12000)
        varied_field = kind in ("field", "geometry")
        for key, expected in (("cap_solves", 2), ("wall_solves", 2 if varied_field else 1),
                              ("field_solves", 2 if varied_field else 1)):
            if type(execution[key]) is not int or execution[key] != expected:
                raise ValueError("cap comparison solve count mismatch")
    except (KeyError, TypeError, IndexError) as error:
        raise ValueError("incomplete cap comparison record") from error
    return True


def create(directory, kind):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only record: choose a fresh directory or --check")
    if kind not in VARIATIONS:
        raise ValueError("cap comparison is not admitted")
    bindings, runtime = source_bindings(), runtime_versions()
    record = compare(kind)
    if bindings != source_bindings() or runtime != runtime_versions() or record["runtime"] != runtime:
        raise RuntimeError("calculation inputs or runtime changed during comparison")
    validate_record(record, kind)
    raw = encoded(record)
    manifest = {"schema": SCHEMA+"-manifest", "source_bindings": bindings,
        "comparison": kind, "outputs": {"response.json": {"sha256": digest(raw), "bytes": len(raw)}}}
    directory.mkdir(parents=True, exist_ok=False)
    (directory/"response.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return {"created": str(directory), "comparison": kind,
            "execution": record["execution"], "differences": record["cap_differences"]}


def check(directory):
    directory = Path(directory)
    manifest = json.loads((directory/"manifest.json").read_bytes())
    if (manifest.get("schema") != SCHEMA+"-manifest"
            or manifest.get("source_bindings") != source_bindings()
            or manifest.get("comparison") not in VARIATIONS):
        raise ValueError("cap comparison source binding or schema mismatch")
    if ({p.name for p in directory.iterdir()} != {"manifest.json", "response.json"}
            or set(manifest.get("outputs", {})) != {"response.json"}):
        raise ValueError("cap comparison inventory mismatch")
    path = directory/"response.json"
    raw = path.read_bytes()
    if path.is_symlink() or manifest["outputs"]["response.json"] != {"sha256": digest(raw), "bytes": len(raw)}:
        raise ValueError("cap comparison output mismatch")
    validate_record(json.loads(raw), manifest["comparison"])
    return {"checked": str(directory), "comparison": manifest["comparison"], "cap_solves": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", choices=VARIATIONS, default="field")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    print(json.dumps(check(args.output) if args.check else create(args.output, args.comparison), sort_keys=True))


if __name__ == "__main__":
    main()
