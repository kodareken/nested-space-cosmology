#!/usr/bin/env python3
"""Record one bounded, independently varied exterior-metric comparison.

Each comparison starts from the same physical preparation. Creation preserves
existing records; --check authenticates their inputs and outputs without a solve.
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
SCHEMA = "nsc-atlas-bulk-response-comparison-v1"
OWNER = "docs/nsc-atlas-bulk-response.md"
SOURCES = (
    OWNER, "scripts/record_nsc_atlas_bulk_response.py", "pyproject.toml",
    "scripts/run_python.py",
    "src/recursive_horizons/nsc_atlas_bulk_response.py",
    "src/recursive_horizons/nsc_atlas_prepared_response.py",
    "src/recursive_horizons/nsc_atlas_characteristic_initialization.py",
    "src/recursive_horizons/nsc_atlas_field_preparation.py",
    "src/recursive_horizons/nsc_atlas_field_transfer.py",
    "src/recursive_horizons/nsc_atlas_causal_reference.py",
    "src/recursive_horizons/nsc_atlas_causal_observers.py",
    "src/recursive_horizons/nsc_atlas_connection_slice.py",
    "src/recursive_horizons/nsc_atlas_moving_boundaries.py",
    "tests/test_nsc_atlas_bulk_response.py",
    "tests/test_nsc_atlas_prepared_response.py",
    "tests/test_record_nsc_atlas_bulk_response.py",
)
BASE = dict(metric_nodes=65, field_nodes=257, geometry_knots=257,
            field_order=32, source_order=64, tolerance=2e-9)
VARIATIONS = {
    "metric": ("metric_nodes", 33, 65),
    "field": ("field_nodes", 129, 257),
    "geometry": ("geometry_knots", 129, 257),
    "derivative": ("field_order", 32, 64),
    "source": ("source_order", 32, 64),
    "time": ("tolerance", 2e-9, 2e-10),
}
QUANTITIES = ("h_J", "c_J", "h_J_U", "h_J_V", "c_J_U", "c_J_V", "mu_E", "mu_J")
TIMES = tuple(float(t) for t in np.linspace(0, 1, 9))
STATUS = "FORMAL_EXTERIOR_RECTANGLE_REFINEMENT"
ORDER = "r_J=r0+alpha*h_J; log C_J=log C0+alpha*c_J; alpha=kappa4*epsilon^2/ell^2"
INTERPRETATION = "dimensionless quadratic coefficients; no finite epsilon selected"
SCOPE = "formal exterior rectangle; U=Ua data for subsequent wall matching"
DEFECTS = ("mass_defect_U", "mass_defect_V", "null_U_integral_defect",
           "null_V_integral_defect", "areal_null_lapse_defect", "frame_mass_defect")


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+"\n").encode()


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def runtime_versions():
    from scipy import __version__ as scipy_version
    return {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy_version}


def leg_differences(summaries):
    result = {}
    for key in QUANTITIES:
        left, right = (np.array([row[key] for row in summary["transverse_leg_samples"]])
                       for summary in summaries)
        difference = float(np.max(abs(left-right)))
        scale = float(max(np.max(abs(left)), np.max(abs(right))))
        result[key] = {"maximum_absolute_difference": difference,
                       "maximum_coefficient_magnitude": scale,
                       "fraction_of_maximum_magnitude": difference/scale if scale else 0.}
    return result


def validate_record(record, kind):
    """Check the declared experiment and saved summaries, without a PDE replay.

    This cannot independently reconstruct unsaved bulk arrays or establish a
    continuum bound. Runtime versions describe the producing environment;
    authenticating a record need not use the same installed packages.
    """
    def real(value, label, *, nonnegative=False):
        if (isinstance(value, bool) or not isinstance(value, (int, float))
                or not math.isfinite(value) or (nonnegative and value < 0)):
            raise ValueError("invalid bulk comparison "+label)

    def count(value, label, *, maximum):
        if type(value) is not int or not 0 <= value <= maximum:
            raise ValueError("invalid bulk comparison "+label)

    try:
        if kind not in VARIATIONS or not isinstance(record, dict):
            raise ValueError("bulk comparison declaration mismatch")
        declarations = {"schema": SCHEMA, "comparison": kind, "status": STATUS,
                        "order": ORDER, "interpretation": INTERPRETATION}
        if any(record.get(key) != value for key, value in declarations.items()):
            raise ValueError("bulk comparison declaration mismatch")
        if record.get("refinement_is_rigorous_error_enclosure") is not False:
            raise ValueError("bulk comparison non-enclosure declaration mismatch")
        varied, low, high = VARIATIONS[kind]
        expected = [{**BASE, varied: value} for value in (low, high)]
        if record["varied_quantity"] != varied or record["settings"] != expected:
            raise ValueError("bulk comparison settings mismatch")
        for setting in record["settings"]:
            for key in BASE:
                if key != "tolerance" and type(setting[key]) is not int:
                    raise ValueError("bulk comparison settings must use integer grids and orders")
        runtime = record["runtime"]
        if (set(runtime) != {"python", "numpy", "scipy"}
                or any(not isinstance(value, str) or not value or not value[0].isdigit()
                       or "." not in value for value in runtime.values())):
            raise ValueError("bulk comparison runtime declaration mismatch")
        solutions = record["solutions"]
        if not isinstance(solutions, list) or len(solutions) != 2:
            raise ValueError("bulk comparison needs two solutions")
        for solution, setting in zip(solutions, expected):
            if (solution["schema"] != "nsc-atlas-bulk-response-v1" or solution["scope"] != SCOPE
                    or any(solution[key] != setting[key] for key in ("metric_nodes", "source_order", "field_order"))):
                raise ValueError("bulk comparison solution scope or settings mismatch")
            domain = solution["domain"]
            if set(domain) != {"U0", "Ua", "Va", "Vb"} or domain != solutions[0]["domain"]:
                raise ValueError("bulk comparison domains differ")
            for value in domain.values():
                real(value, "domain")
            if not (domain["U0"] < domain["Ua"] < 0 and domain["Va"] < domain["Vb"] < 0
                    and domain["U0"]*domain["Va"] < 1):
                raise ValueError("bulk comparison is outside the retained IV rectangle")
            defects = solution["maximum_absolute_defects"]
            if set(defects) != set(DEFECTS):
                raise ValueError("bulk comparison defect inventory mismatch")
            for value in defects.values():
                real(value, "defect", nonnegative=True)
            rows = solution["transverse_leg_samples"]
            if not isinstance(rows, list) or len(rows) != len(TIMES):
                raise ValueError("bulk comparison sample inventory mismatch")
            for row, t in zip(rows, TIMES):
                if set(row) != {"t", "mass_defect_U_gauss", *QUANTITIES} or row["t"] != t:
                    raise ValueError("bulk comparison sample coordinates mismatch")
                for value in row.values():
                    real(value, "sample")
        if record["transverse_leg_differences"] != leg_differences(solutions):
            raise ValueError("bulk comparison saved endpoint differences disagree")
        common = record["common_grid_maximum_absolute_differences"]
        if set(common) != set(QUANTITIES):
            raise ValueError("bulk comparison common-grid inventory mismatch")
        for key, value in common.items():
            real(value, "common-grid difference", nonnegative=True)
            if value < record["transverse_leg_differences"][key]["maximum_absolute_difference"]:
                raise ValueError("bulk comparison common grid omits its endpoint difference")
        execution = record["execution"]
        if execution["cpu_seconds_max"] != 120 or execution["scalar_kernel_values_max"] != 30_000_000:
            raise ValueError("bulk comparison budget ceilings changed")
        real(execution["cpu_seconds"], "CPU time", nonnegative=True)
        if execution["cpu_seconds"] > execution["cpu_seconds_max"]:
            raise ValueError("bulk comparison CPU budget exceeded")
        for key in ("scalar_kernel_values", "metric_values", "source_quadrature_nodes"):
            count(execution[key], "budget count", maximum=30_000_000)
        if execution["metric_values"]+execution["source_quadrature_nodes"] > 30_000_000:
            raise ValueError("bulk comparison metric budget exceeded")
        count(execution["metric_rhs_calls"], "RHS count", maximum=12_000)
        for key, value in (("metric_solves", 2), ("field_solves", 2 if kind in ("field", "geometry") else 1)):
            if type(execution[key]) is not int or execution[key] != value:
                raise ValueError("bulk comparison solve count mismatch")
    except (KeyError, TypeError, IndexError) as error:
        raise ValueError("incomplete bulk comparison record") from error
    return True


def compare(kind):
    from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget, solve_prepared_field
    from recursive_horizons.nsc_atlas_bulk_response import BulkBudget, solve_bulk_rectangle
    from recursive_horizons.nsc_atlas_characteristic_initialization import source_geometry

    varied, low, high = VARIATIONS[kind]
    shared = PreparedFieldBudget()
    geometry = source_geometry()
    metrics = BulkBudget(started=shared.started, external=shared)
    fields = {}
    solutions, summaries, settings = [], [], []
    for value in (low, high):
        config = {**BASE, varied: value}
        key = config["field_nodes"], config["geometry_knots"]
        if key not in fields:
            fields[key] = solve_prepared_field(*key, budget=shared)
        solution = solve_bulk_rectangle(
            fields[key], geometry, nodes=config["metric_nodes"],
            source_order=config["source_order"], field_order=config["field_order"],
            tolerance=config["tolerance"], budget=metrics)
        summary = solution.summary(TIMES, independent_order=64)
        # Execution counters describe the whole comparison below.
        for name in ("rhs_calls", "metric_values", "source_nodes", "process_cpu_seconds"):
            summary.pop(name)
        summaries.append(summary)
        solutions.append(solution)
        settings.append(config)
    differences = leg_differences(summaries)
    # Compare metric fields on common nodes, not just the next initial leg.
    common_differences = dict.fromkeys(QUANTITIES, 0.)
    for t in TIMES:
        left, right = (solution.sample(t) for solution in solutions)
        if len(right["U"]) % (len(left["U"])-1) != 1:
            raise ValueError("comparison grids do not share exact nested nodes")
        stride = (len(right["U"])-1)//(len(left["U"])-1)
        if not np.array_equal(left["U"], right["U"][::stride]):
            raise ValueError("comparison coordinates do not agree")
        for key in QUANTITIES:
            common_differences[key] = max(common_differences[key],
                                         float(np.max(abs(left[key]-right[key][::stride]))))
    shared.check()
    metrics.check()
    return {
        "schema": SCHEMA, "comparison": kind, "varied_quantity": varied,
        "status": STATUS, "order": ORDER, "interpretation": INTERPRETATION,
        "runtime": runtime_versions(),
        "settings": settings, "solutions": summaries,
        "transverse_leg_differences": differences,
        "common_grid_maximum_absolute_differences": common_differences,
        "refinement_is_rigorous_error_enclosure": False,
        "execution": {"cpu_seconds": time.process_time()-shared.started,
                      "cpu_seconds_max": shared.cpu_seconds,
                      "scalar_kernel_values": shared.kernel_values,
                      "scalar_kernel_values_max": shared.kernel_values_max,
                      "metric_values": metrics.metric_values,
                      "source_quadrature_nodes": metrics.source_nodes,
                      "metric_rhs_calls": metrics.rhs_calls,
                      "field_solves": len(fields), "metric_solves": 2},
    }


def create(directory, kind):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only record: choose a fresh directory or --check")
    bindings = source_bindings()
    runtime = runtime_versions()
    record = compare(kind)
    if bindings != source_bindings():
        raise RuntimeError("calculation inputs changed during comparison")
    validate_record(record, kind)
    if record["runtime"] != runtime or runtime_versions() != runtime:
        raise RuntimeError("calculation runtime changed during comparison")
    raw = encoded(record)
    manifest = {"schema": SCHEMA+"-manifest", "source_bindings": bindings,
                "comparison": kind,
                "outputs": {"response.json": {"sha256": digest(raw), "bytes": len(raw)}}}
    directory.mkdir(parents=True, exist_ok=False)
    (directory/"response.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return {"created": str(directory), "comparison": kind, "execution": record["execution"],
            "differences": record["common_grid_maximum_absolute_differences"]}


def check(directory):
    directory = Path(directory)
    manifest = json.loads((directory/"manifest.json").read_bytes())
    if (manifest.get("schema") != SCHEMA+"-manifest"
            or manifest.get("source_bindings") != source_bindings()
            or manifest.get("comparison") not in VARIATIONS):
        raise ValueError("bulk comparison source binding or schema mismatch")
    if ({p.name for p in directory.iterdir()} != {"manifest.json", "response.json"}
            or set(manifest.get("outputs", {})) != {"response.json"}):
        raise ValueError("bulk comparison inventory mismatch")
    path = directory/"response.json"
    raw = path.read_bytes()
    if path.is_symlink() or manifest["outputs"]["response.json"] != {"sha256": digest(raw), "bytes": len(raw)}:
        raise ValueError("bulk comparison output mismatch")
    record = json.loads(raw)
    if record.get("schema") != SCHEMA or record.get("comparison") != manifest["comparison"]:
        raise ValueError("bulk comparison declaration mismatch")
    validate_record(record, manifest["comparison"])
    return {"checked": str(directory), "comparison": record["comparison"], "metric_solves": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", choices=VARIATIONS, default="metric")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    print(json.dumps(check(args.output) if args.check else create(args.output, args.comparison), sort_keys=True))


if __name__ == "__main__":
    main()
