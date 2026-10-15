#!/usr/bin/env python3
"""Save one bounded comparison of the same prepared field and moving wall.

Creation is explicit and preserves existing directories. Authentication checks
saved declarations and bytes without performing another physical calculation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-wall-evolution-comparison-v1"
SOURCES = (
    "docs/nsc-atlas-wall-evolution.md", "scripts/record_nsc_atlas_wall_evolution.py",
    "pyproject.toml", "scripts/run_python.py",
    *("src/recursive_horizons/"+name+".py" for name in (
        "nsc_atlas_wall_evolution", "nsc_atlas_wedge_field", "nsc_atlas_prepared_response",
        "nsc_atlas_bulk_response", "nsc_atlas_cap_wall_response",
        "nsc_atlas_characteristic_initialization", "nsc_atlas_field_preparation",
        "nsc_atlas_field_transfer", "nsc_atlas_causal_reference",
        "nsc_atlas_causal_observers", "nsc_atlas_connection_slice",
        "nsc_atlas_moving_boundaries")),
    "tests/test_nsc_atlas_wall_evolution.py", "tests/test_nsc_atlas_wedge_field.py",
    "tests/test_nsc_atlas_prepared_response.py", "tests/test_nsc_atlas_bulk_response.py",
    "tests/test_nsc_atlas_cap_wall_response.py",
    "tests/test_record_nsc_atlas_wall_evolution.py",
)
BASE = dict(field_nodes=257, geometry_knots=257, radial_order=4,
            field_order=32, time_panels=16, source_order=64, tolerance=2e-9)
VARIATIONS = {
    "field": ("field_nodes", 129, 257),
    "radial": ("radial_order", 2, 4),
    "derivative": ("field_order", 32, 64),
    "time": ("time_panels", 8, 16),
    "geometry": ("geometry_knots", 129, 257),
    "source": ("source_order", 32, 64),
    "ode": ("tolerance", 2e-9, 2e-10),
}
TIMES = tuple(float(t) for t in np.linspace(0, 1, 9))
STATUS = "FORMAL_MOVING_WALL_REFINEMENT"
ORDER = "r_W=x+alpha*rho; V_W=V0*(1+alpha*eta); alpha=kappa4*epsilon^2/ell^2"
INTERPRETATION = "dimensionless quadratic coefficients; no finite epsilon selected"
SOLUTION_SCHEMA = "nsc-atlas-wall-evolution-v1"
SCOPE = "formal moving-wall and exterior V-clock coefficients; cap metric and central observer remain open"
QUANTITIES = ("mu_plus", "nu_minus", "rho", "rho_s", "delta_beta_plus", "zeta_J", "eta", "V1")
SAMPLE_KEYS = {"t", "U", "V", "x", "beta_plus", "mu_plus", "nu_minus",
               "cap_total_mass", "rho", "rho_s", "displacement_rate", "delta_beta_plus",
               "delta_normal_acceleration", "radial_mu_E", "radial_mu_J", "zeta_E", "zeta_J",
               "eta", "eta_s", "V1", "V1_s", "mass_path_defect"}
DEFECTS = {"mass_path_defect", "frame_mass_defect", "wall_lapse_frame_defect",
           "wall_physical_p_U_defect", "reference_clock_normalization_defect",
           "angular_jump_defect", "two_displacement_equations_defect",
           "clock_normalization_defect", "chi_U_lower_limit_bound"}


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+"\n").encode()


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def runtime_versions():
    from scipy import __version__ as scipy_version
    return {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy_version}


def differences(summaries):
    """Compare like coefficients at the declared common time cuts."""
    left, right = (summary["wall_samples"] for summary in summaries)
    result = {}
    for key in QUANTITIES:
        a, b = (np.asarray([row[key] for row in rows]) for rows in (left, right))
        error = float(np.max(abs(a-b)))
        scale = float(max(np.max(abs(a)), np.max(abs(b))))
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

    varied, low, high = VARIATIONS[kind]
    shared = PreparedFieldBudget()
    metrics = BulkBudget(started=shared.started, external=shared)
    geometry = source_geometry()
    fields, summaries, settings = {}, [], []
    for value in (low, high):
        config = {**BASE, varied: value}
        key = config["field_nodes"], config["geometry_knots"]
        if key not in fields:
            fields[key] = WedgeField(solve_prepared_field(*key, budget=shared))
        solution = solve_wall_evolution(fields[key], geometry, budget=metrics,
            **{name: config[name] for name in BASE if name not in ("field_nodes", "geometry_knots")})
        summary = solution.summary(TIMES)
        # These are cumulative shared counters, not costs of this solution.
        for name in ("rhs_calls", "metric_values", "source_nodes", "process_cpu_seconds"):
            summary.pop(name)
        summaries.append(summary)
        settings.append(config)
    shared.check()
    metrics.check()
    return {"schema": SCHEMA, "comparison": kind, "varied_quantity": varied,
            "status": STATUS, "order": ORDER, "interpretation": INTERPRETATION,
            "runtime": runtime_versions(), "settings": settings, "solutions": summaries,
            "wall_differences": differences(summaries),
            "refinement_is_rigorous_error_enclosure": False,
            "execution": {"cpu_seconds": time.process_time()-shared.started,
                "cpu_seconds_max": shared.cpu_seconds,
                "scalar_kernel_values": shared.kernel_values,
                "scalar_kernel_values_max": shared.kernel_values_max,
                "metric_values": metrics.metric_values, "source_quadrature_nodes": metrics.source_nodes,
                "rhs_calls": metrics.rhs_calls, "field_solves": len(fields), "wall_solves": 2}}


def validate_record(record, kind):
    """Authenticate the stated finite experiment, without numerical replay.

    Saved defects and differences are diagnostics, not a continuum enclosure
    or verification of the unsaved trajectory. Runtime describes its producer.
    """
    def real(value, label, *, nonnegative=False):
        if (isinstance(value, bool) or not isinstance(value, (int, float))
                or not math.isfinite(value) or (nonnegative and value < 0)):
            raise ValueError("invalid wall comparison "+label)

    def count(value, label, maximum):
        if type(value) is not int or not 0 <= value <= maximum:
            raise ValueError("invalid wall comparison "+label)

    try:
        if kind not in VARIATIONS or not isinstance(record, dict):
            raise ValueError("wall comparison declaration mismatch")
        for key, expected in {"schema": SCHEMA, "comparison": kind, "status": STATUS,
                              "order": ORDER, "interpretation": INTERPRETATION}.items():
            if record.get(key) != expected:
                raise ValueError("wall comparison declaration mismatch")
        if record.get("refinement_is_rigorous_error_enclosure") is not False:
            raise ValueError("wall comparison non-enclosure declaration mismatch")
        varied, low, high = VARIATIONS[kind]
        settings = [{**BASE, varied: value} for value in (low, high)]
        if record["varied_quantity"] != varied or record["settings"] != settings:
            raise ValueError("wall comparison settings mismatch")
        for setting in record["settings"]:
            if any(type(setting[key]) is not int for key in BASE if key != "tolerance"):
                raise ValueError("wall comparison integer setting mismatch")
        runtime = record["runtime"]
        if (not isinstance(runtime, dict) or set(runtime) != {"python", "numpy", "scipy"}
                or any(not isinstance(v, str) or not v or not v[0].isdigit() or "." not in v
                       for v in runtime.values())):
            raise ValueError("wall comparison runtime mismatch")
        solutions = record["solutions"]
        if not isinstance(solutions, list) or len(solutions) != 2:
            raise ValueError("wall comparison requires two solutions")
        for solution, setting in zip(solutions, settings):
            if (solution["schema"] != SOLUTION_SCHEMA or solution["scope"] != SCOPE
                    or solution["source_scale"] != 1. or isinstance(solution["source_scale"], bool)
                    or type(solution["clock_time_order"]) is not int or solution["clock_time_order"] != 4
                    or any(solution[key] != setting[key] for key in
                           ("radial_order", "field_order", "time_panels", "source_order", "tolerance"))):
                raise ValueError("wall comparison solution settings or scope mismatch")
            domain = solution["domain"]
            if (not isinstance(domain, dict) or set(domain) != {"U0", "Ua", "Ub", "Va", "Vb"}
                    or domain != solutions[0]["domain"]):
                raise ValueError("wall comparison domains differ")
            for value in domain.values():
                real(value, "domain")
            if not (domain["U0"] < domain["Ua"] < domain["Ub"] < 0
                    and domain["Va"] < domain["Vb"] < 0 and domain["U0"]*domain["Va"] < 1):
                raise ValueError("wall comparison domain exits the retained IV wedge")
            defects = solution["maximum_absolute_defects"]
            if not isinstance(defects, dict) or set(defects) != DEFECTS:
                raise ValueError("wall comparison defect inventory mismatch")
            for value in defects.values():
                real(value, "defect", nonnegative=True)
            rows = solution["wall_samples"]
            if not isinstance(rows, list) or len(rows) != len(TIMES):
                raise ValueError("wall comparison sample coordinates mismatch")
            for row, t in zip(rows, TIMES):
                if not isinstance(row, dict) or set(row) != SAMPLE_KEYS or row["t"] != t:
                    raise ValueError("wall comparison sample inventory or coordinates mismatch")
                for key, value in row.items():
                    real(value, "sample "+key)
                if (not domain["Ua"]-1e-13 <= row["U"] <= domain["Ub"]+1e-13
                        or not .45 <= row["x"] <= .55
                        or abs(row["V"]-(domain["Va"]+(domain["Vb"]-domain["Va"])*t)) > 1e-13):
                    raise ValueError("wall comparison samples leave the declared strip")
                if row["rho_s"] != row["displacement_rate"]:
                    raise ValueError("wall comparison displacement derivative aliases disagree")
                if abs(row["mass_path_defect"]) > defects["mass_path_defect"]:
                    raise ValueError("wall comparison mass diagnostic understates its displayed defect")
            if (abs(rows[0]["U"]-domain["Ua"]) > 1e-13
                    or abs(rows[-1]["U"]-domain["Ub"]) > 1e-13
                    or any(b["U"] <= a["U"] for a, b in zip(rows[:-1], rows[1:]))):
                raise ValueError("wall comparison wall trace is not the ordered finite interval")
            expected_queries = 4*setting["time_panels"]+len(TIMES)
            if type(solution["radial_evaluations"]) is not int or solution["radial_evaluations"] != expected_queries:
                raise ValueError("wall comparison radial query count mismatch")
        if record["wall_differences"] != differences(solutions):
            raise ValueError("wall comparison saved differences disagree")
        execution = record["execution"]
        if execution["cpu_seconds_max"] != 120 or execution["scalar_kernel_values_max"] != 30_000_000:
            raise ValueError("wall comparison budget ceilings changed")
        real(execution["cpu_seconds"], "CPU time", nonnegative=True)
        if execution["cpu_seconds"] > 120:
            raise ValueError("wall comparison CPU budget exceeded")
        for key in ("scalar_kernel_values", "metric_values", "source_quadrature_nodes"):
            count(execution[key], key, 30_000_000)
        if execution["metric_values"]+execution["source_quadrature_nodes"] > 30_000_000:
            raise ValueError("wall comparison metric budget exceeded")
        count(execution["rhs_calls"], "RHS calls", 12000)
        for key, expected in (("wall_solves", 2), ("field_solves", 2 if kind in ("field", "geometry") else 1)):
            if type(execution[key]) is not int or execution[key] != expected:
                raise ValueError("wall comparison solve count mismatch")
    except (KeyError, TypeError, IndexError) as error:
        raise ValueError("incomplete wall comparison record") from error
    return True


def create(directory, kind):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only record: choose a fresh directory or --check")
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
            "execution": record["execution"], "differences": record["wall_differences"]}


def check(directory):
    directory = Path(directory)
    manifest = json.loads((directory/"manifest.json").read_bytes())
    if (manifest.get("schema") != SCHEMA+"-manifest"
            or manifest.get("source_bindings") != source_bindings()
            or manifest.get("comparison") not in VARIATIONS):
        raise ValueError("wall comparison source binding or schema mismatch")
    if ({p.name for p in directory.iterdir()} != {"manifest.json", "response.json"}
            or set(manifest.get("outputs", {})) != {"response.json"}):
        raise ValueError("wall comparison inventory mismatch")
    path = directory/"response.json"
    raw = path.read_bytes()
    if path.is_symlink() or manifest["outputs"]["response.json"] != {"sha256": digest(raw), "bytes": len(raw)}:
        raise ValueError("wall comparison output mismatch")
    validate_record(json.loads(raw), manifest["comparison"])
    return {"checked": str(directory), "comparison": manifest["comparison"], "wall_solves": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", choices=VARIATIONS, default="field")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    print(json.dumps(check(args.output) if args.check else create(args.output, args.comparison), sort_keys=True))


if __name__ == "__main__":
    main()
