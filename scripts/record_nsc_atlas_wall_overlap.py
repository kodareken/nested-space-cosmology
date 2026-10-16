#!/usr/bin/env python3
"""Record the independent radial-constraint/rectangle-metric overlap.

Creation runs one fixed finite comparison. Checking authenticates source and
output bytes and the experiment declaration, without a field or metric solve.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

from record_nsc_atlas_wall_evolution import (
    SOURCES as WALL_SOURCES, encoded, digest, runtime_versions,
)


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-wall-rectangle-overlap-v1"
STATUS = "FORMAL_RADIAL_RECTANGLE_OVERLAP"
SCOPE = "quadratic exterior coefficients at U=Ua; no full spacetime or finite-amplitude error enclosure"
SOURCES = tuple(dict.fromkeys((*WALL_SOURCES,
    "scripts/record_nsc_atlas_wall_overlap.py", "tests/test_record_nsc_atlas_wall_overlap.py")))
SETTINGS = dict(field_nodes=257, geometry_knots=257, metric_nodes=65,
                radial_order=4, field_order=32, source_order=64, tolerance=2e-9)
TIMES = tuple(float(t) for t in np.linspace(0, 1, 9))
DIAGNOSTICS = ("mu_E_difference", "mu_J_difference", "zeta_E_difference",
               "rectangle_mass_defect_U", "rectangle_mass_defect_V", "rectangle_lapse_defect")
SAMPLE_KEYS = {"t", "U", "V", *DIAGNOSTICS}


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def maxima(rows):
    return {key: max(abs(row[key]) for row in rows) for key in DIAGNOSTICS}


def calculate():
    from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget, solve_prepared_field
    from recursive_horizons.nsc_atlas_characteristic_initialization import source_geometry
    from recursive_horizons.nsc_atlas_wedge_field import WedgeField
    from recursive_horizons.nsc_atlas_bulk_response import BulkBudget, solve_bulk_rectangle
    from recursive_horizons.nsc_atlas_wall_evolution import radial_response

    shared = PreparedFieldBudget()
    budget = BulkBudget(started=shared.started, external=shared)
    field = solve_prepared_field(SETTINGS["field_nodes"], SETTINGS["geometry_knots"], budget=shared)
    geometry = source_geometry()
    rectangle = solve_bulk_rectangle(field, geometry, nodes=SETTINGS["metric_nodes"],
        source_order=SETTINGS["source_order"], field_order=SETTINGS["field_order"],
        tolerance=SETTINGS["tolerance"], budget=budget)
    wedge = WedgeField(field)
    rows = []
    for t in TIMES:
        metric = rectangle.sample(t)
        radial = radial_response(wedge, geometry, t, upper_U=geometry["Ua"],
            radial_order=SETTINGS["radial_order"], field_order=SETTINGS["field_order"],
            source_order=SETTINGS["source_order"], budget=budget)
        rows.append({"t": t, "U": geometry["Ua"], "V": radial["V"],
            "mu_E_difference": radial["mu_E"]-float(metric["mu_E"][-1]),
            "mu_J_difference": radial["mu_J"]-float(metric["mu_J"][-1]),
            "zeta_E_difference": radial["zeta_E"]-float(metric["zeta_E"][-1]),
            "rectangle_mass_defect_U": float(metric["mass_defect_U"][-1]),
            "rectangle_mass_defect_V": float(metric["mass_defect_V"][-1]),
            "rectangle_lapse_defect": float(metric["areal_null_lapse_defect"][-1])})
    shared.check()
    budget.check()
    return {"schema": SCHEMA, "status": STATUS, "scope": SCOPE,
        "runtime": runtime_versions(), "settings": dict(SETTINGS),
        "domain": {key: geometry[key] for key in ("U0", "Ua", "Ub", "Va", "Vb")},
        "difference_convention": "radial constraints minus propagated rectangle metric at U=Ua",
        "refinement_is_rigorous_error_enclosure": False,
        "samples": rows, "maximum_absolute": maxima(rows),
        "execution": {"cpu_seconds": time.process_time()-shared.started,
            "cpu_seconds_max": shared.cpu_seconds, "scalar_kernel_values": shared.kernel_values,
            "scalar_kernel_values_max": shared.kernel_values_max,
            "metric_values": budget.metric_values, "source_quadrature_nodes": budget.source_nodes,
            "rhs_calls": budget.rhs_calls, "field_solves": 1, "metric_solves": 1,
            "wall_solves": 0, "radial_evaluations": len(TIMES)}}


def validate_record(record):
    def real(value, label, nonnegative=False):
        if (isinstance(value, bool) or not isinstance(value, (int, float))
                or not math.isfinite(value) or (nonnegative and value < 0)):
            raise ValueError("invalid overlap "+label)

    def count(value, label, maximum):
        if type(value) is not int or not 0 <= value <= maximum:
            raise ValueError("invalid overlap "+label)

    try:
        expected = {"schema": SCHEMA, "status": STATUS, "scope": SCOPE,
            "difference_convention": "radial constraints minus propagated rectangle metric at U=Ua"}
        if not isinstance(record, dict) or any(record.get(key) != value for key, value in expected.items()):
            raise ValueError("overlap declaration mismatch")
        if record.get("refinement_is_rigorous_error_enclosure") is not False:
            raise ValueError("overlap non-enclosure declaration mismatch")
        if (record["settings"] != SETTINGS
                or any(type(record["settings"][key]) is not int for key in SETTINGS if key != "tolerance")):
            raise ValueError("overlap settings mismatch")
        runtime = record["runtime"]
        if (not isinstance(runtime, dict) or set(runtime) != {"python", "numpy", "scipy"}
                or any(not isinstance(value, str) or not value or not value[0].isdigit()
                       or "." not in value for value in runtime.values())):
            raise ValueError("overlap runtime mismatch")
        domain = record["domain"]
        if not isinstance(domain, dict) or set(domain) != {"U0", "Ua", "Ub", "Va", "Vb"}:
            raise ValueError("overlap domain inventory mismatch")
        for value in domain.values():
            real(value, "domain")
        if not (domain["U0"] < domain["Ua"] < domain["Ub"] < 0
                and domain["Va"] < domain["Vb"] < 0 and domain["U0"]*domain["Va"] < 1):
            raise ValueError("overlap domain exits the finite IV wedge")
        rows = record["samples"]
        if not isinstance(rows, list) or len(rows) != len(TIMES):
            raise ValueError("overlap requires the nine fixed time cuts")
        for row, t in zip(rows, TIMES):
            if not isinstance(row, dict) or set(row) != SAMPLE_KEYS or row["t"] != t:
                raise ValueError("overlap sample inventory or coordinates mismatch")
            for key, value in row.items():
                real(value, "sample "+key)
            if (row["U"] != domain["Ua"]
                    or abs(row["V"]-domain["Va"]-(domain["Vb"]-domain["Va"])*t) > 1e-13):
                raise ValueError("overlap samples exit the stated transverse leg")
        if record["maximum_absolute"] != maxima(rows):
            raise ValueError("overlap diagnostic maxima disagree with samples")
        execution = record["execution"]
        if execution["cpu_seconds_max"] != 120 or execution["scalar_kernel_values_max"] != 30_000_000:
            raise ValueError("overlap budget ceilings changed")
        real(execution["cpu_seconds"], "CPU time", nonnegative=True)
        if execution["cpu_seconds"] > 120:
            raise ValueError("overlap CPU budget exceeded")
        for key in ("scalar_kernel_values", "metric_values", "source_quadrature_nodes"):
            count(execution[key], key, 30_000_000)
        if execution["metric_values"]+execution["source_quadrature_nodes"] > 30_000_000:
            raise ValueError("overlap metric budget exceeded")
        count(execution["rhs_calls"], "RHS calls", 12000)
        for key, expected_count in (("field_solves", 1), ("metric_solves", 1),
                                    ("wall_solves", 0), ("radial_evaluations", len(TIMES))):
            if type(execution[key]) is not int or execution[key] != expected_count:
                raise ValueError("overlap solve/query count mismatch")
    except (KeyError, TypeError, IndexError) as error:
        raise ValueError("incomplete overlap record") from error
    return True


def create(directory):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only record: choose a fresh directory or --check")
    bindings, runtime = source_bindings(), runtime_versions()
    record = calculate()
    if bindings != source_bindings() or runtime != runtime_versions() or record["runtime"] != runtime:
        raise RuntimeError("overlap inputs or runtime changed during calculation")
    validate_record(record)
    raw = encoded(record)
    manifest = {"schema": SCHEMA+"-manifest", "source_bindings": bindings,
        "outputs": {"response.json": {"sha256": digest(raw), "bytes": len(raw)}}}
    directory.mkdir(parents=True, exist_ok=False)
    (directory/"response.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return {"created": str(directory), "execution": record["execution"],
            "maximum_absolute": record["maximum_absolute"]}


def check(directory):
    from validation_source_successors import verify_saved_bindings
    directory = Path(directory)
    if ({path.name for path in directory.iterdir()} != {"manifest.json", "response.json"}
            or any((directory/name).is_symlink() for name in ("manifest.json", "response.json"))):
        raise ValueError("overlap output inventory mismatch")
    manifest = json.loads((directory/"manifest.json").read_bytes())
    if manifest.get("schema") != SCHEMA+"-manifest":
        raise ValueError("overlap source binding or schema mismatch")
    verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=SOURCES)
    raw = (directory/"response.json").read_bytes()
    if manifest.get("outputs") != {"response.json": {"sha256": digest(raw), "bytes": len(raw)}}:
        raise ValueError("overlap output mismatch")
    validate_record(json.loads(raw))
    context="current" if manifest["source_bindings"] == source_bindings() else "reviewed-successors"
    return {"checked": str(directory), "field_solves": 0, "metric_solves": 0, "wall_solves": 0,
            "source_context": context}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    print(json.dumps(check(args.output) if args.check else create(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
