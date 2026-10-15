#!/usr/bin/env python3
"""Assess local Hermite fits to authenticated saved boundary jets; no solves.

Only the original 17 knots train each fit. The original 16 withheld midpoints
remain validation data. Stencils cannot cross either source edge, and the
reported errors are sampled diagnostics, not uniform error enclosures.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.interpolate import KroghInterpolator

import record_nsc_atlas_scalar_inputs as predecessor
from record_nsc_atlas_scalar_inputs import digest, encoded, runtime_versions

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-scalar-lift-assessment-v1"
SCOPE = "saved-data interpolation sensitivity; no new samples, physical solve or uniform enclosure"
INPUT = "results/development/nsc-atlas-scalar-inputs-v1/lift"
SOURCES = ("scripts/assess_nsc_atlas_scalar_lift.py",
           "tests/test_assess_nsc_atlas_scalar_lift.py")
BLOCKS = ((0., .25), (.25, .75), (.75, 1.))
NODE_COUNTS = (4, 5, 7, 9)
TRAINING_TIMES = tuple(i/16 for i in range(17))
WITHHELD_TIMES = tuple((2*i+1)/32 for i in range(16))
CPU_LIMIT = 30.


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def input_snapshot():
    """Reuse the full predecessor authentication without preparing any field."""
    checked = predecessor.check(ROOT/INPUT)
    if checked["comparison"] != "lift" or checked["physical_solves"] != 0:
        raise ValueError("assessment requires authenticated saved lift data")
    raw = {name: (ROOT/INPUT/name).read_bytes() for name in ("inputs.json", "manifest.json")}
    manifest = json.loads(raw["manifest.json"])
    files = {name: {"sha256": digest(value), "bytes": len(value)} for name, value in raw.items()}
    if manifest["outputs"]["inputs.json"] != files["inputs.json"]:
        raise ValueError("input changed after source authentication")
    binding = {"directory": INPUT, "files": files,
        "source_validation": {"checker": "scripts/record_nsc_atlas_scalar_inputs.py:check",
            "comparison": "lift", "physical_solves": 0,
            "source_bindings_sha256": digest(encoded(manifest["source_bindings"]))}}
    return json.loads(raw["inputs.json"]), binding


def _finite(value):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
        raise ValueError("sample coordinates and jets must be finite numbers")
    return float(value)


def local_prediction(training, query, requested_nodes):
    """Normalized confluent interpolation using a nearest local training stencil.

    The slope error is later scaled by the original bracketing mesh-cell width,
    independently of the width used to normalize the interpolation polynomial.
    """
    if type(requested_nodes) is not int or requested_nodes not in NODE_COUNTS:
        raise ValueError("only the admitted 4/5/7/9-node stencils are supported")
    if [row["t"] for row in training] != list(TRAINING_TIMES):
        raise ValueError("training must be the original 17-knot mesh")
    coordinates = np.array([_finite(row["U"]) for row in training])
    if np.any(np.diff(coordinates) <= 0):
        raise ValueError("training U coordinates must increase strictly")
    for row in training:
        _finite(row["b"])
        _finite(row["b_U"])
    t, U = _finite(query["t"]), _finite(query["U"])
    if t not in WITHHELD_TIMES or t in TRAINING_TIMES:
        raise ValueError("queries must be original withheld midpoints, never training knots")
    block_index = next(i for i, (lo, hi) in enumerate(BLOCKS) if lo < t < hi)
    lo, hi = BLOCKS[block_index]
    available = [i for i, row in enumerate(training) if lo <= row["t"] <= hi]
    chosen = sorted(sorted(available, key=lambda i: (abs(coordinates[i]-U), i))
                    [:min(requested_nodes, len(available))])
    # These checks prohibit extrapolation even for nonuniform saved coordinates.
    if not coordinates[chosen[0]] < U < coordinates[chosen[-1]]:
        raise ValueError("nearest source-block stencil would extrapolate")
    right = int(np.searchsorted(coordinates, U))
    if not (0 < right < len(training)
            and training[right-1]["t"] < t < training[right]["t"]
            and coordinates[right-1] < U < coordinates[right]):
        raise ValueError("withheld coordinate does not lie in its original training cell")
    width = float(coordinates[chosen[-1]]-coordinates[chosen[0]])
    centre = float(coordinates[chosen[0]]+width/2)
    z = np.repeat((coordinates[chosen]-centre)/width, 2)
    values = np.array([[training[i]["b"], width*training[i]["b_U"]] for i in chosen]).ravel()
    fit = KroghInterpolator(z, values)
    z_query = (U-centre)/width
    value, slope = float(fit(z_query)), float(fit.derivative(z_query)/width)
    if not math.isfinite(value) or not math.isfinite(slope):
        raise ValueError("nonfinite local Hermite prediction")
    return {"t": t, "U": U, "block": block_index,
        "training_indices": chosen, "training_times": [training[i]["t"] for i in chosen],
        "nodes": len(chosen), "normalization_centre": centre, "normalization_width": width,
        "cell_width": float(coordinates[right]-coordinates[right-1]),
        "b_interpolated": value, "b_U_interpolated": slope,
        "extrapolated": False}


def assess(record):
    """Deterministic saved-data postprocessing; callers authenticate the record."""
    if record.get("comparison") != "lift":
        raise ValueError("assessment only accepts the lift comparison")
    solutions = [solution for solution in record["solutions"] if solution["knot_count"] == 17]
    if len(solutions) != 1:
        raise ValueError("one original 17-knot solution is required")
    solution = solutions[0]
    training, withheld = solution["knot_samples"], solution["validation_samples"]
    if [row["t"] for row in withheld] != list(WITHHELD_TIMES):
        raise ValueError("validation must retain the original 16 withheld midpoints")
    target = _finite(record["diagnostics"]["interpolation_target"]["absolute_b_target"])
    if target <= 0:
        raise ValueError("sampled interpolation target must be positive")
    methods = []
    for count in NODE_COUNTS:
        rows = []
        for query in withheld:
            row = local_prediction(training, query, count)
            value, slope = _finite(query["b_reference"]), _finite(query["b_U_reference"])
            row.update(b_reference=value, b_U_reference=slope,
                b_defect=row["b_interpolated"]-value,
                scaled_slope_defect=row["cell_width"]*(row["b_U_interpolated"]-slope))
            rows.append(row)
        blocks = []
        for index, interval in enumerate(BLOCKS):
            selected = [row for row in rows if row["block"] == index]
            blocks.append({"interval": list(interval), "withheld_points": len(selected),
                "available_training_nodes": sum(interval[0] <= row["t"] <= interval[1] for row in training),
                "used_training_nodes": min(count, sum(interval[0] <= row["t"] <= interval[1] for row in training)),
                "maximum_absolute_b_defect": max(abs(row["b_defect"]) for row in selected),
                "maximum_absolute_scaled_slope_defect": max(abs(row["scaled_slope_defect"]) for row in selected)})
        maximum = max(abs(row["b_defect"]) for row in rows)
        methods.append({"requested_nodes": count, "blocks": blocks, "withheld_samples": rows,
            "maximum_absolute_b_defect": maximum,
            "maximum_absolute_scaled_slope_defect": max(abs(row["scaled_slope_defect"]) for row in rows),
            "sampled_b_target_met": maximum <= target})
    return {"interpolation_variable": "U_W", "method": "locally normalized confluent Krogh Hermite",
        "training_times": list(TRAINING_TIMES), "withheld_times": list(WITHHELD_TIMES),
        "source_blocks": [list(pair) for pair in BLOCKS], "heldout_used_for_training": False,
        "extrapolated_queries": 0, "slope_error_scale": "original bracketing 17-knot U-cell width",
        "absolute_b_target": target, "is_uniform_error_bound": False,
        "original_cubic_maximum_absolute_b_defect": max(abs(row["b_defect"]) for row in withheld),
        "methods": methods}


def _budget(started):
    elapsed = time.process_time()-started
    if not 0 <= elapsed <= CPU_LIMIT:
        raise RuntimeError("saved-data assessment exceeded its 30-process-CPU-second admission")
    return elapsed


def _summary(record):
    return {"absolute_b_target": record["assessment"]["absolute_b_target"],
        "methods": [{key: method[key] for key in ("requested_nodes", "maximum_absolute_b_defect",
            "maximum_absolute_scaled_slope_defect", "sampled_b_target_met", "blocks")}
            for method in record["assessment"]["methods"]]}


def create(directory):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only assessment: choose a fresh directory or --check")
    started = time.process_time()
    sources, runtime = source_bindings(), runtime_versions()
    saved, binding = input_snapshot()
    assessment = assess(saved)
    if sources != source_bindings() or binding != input_snapshot()[1] or runtime != runtime_versions():
        raise RuntimeError("assessment inputs or runtime changed during postprocessing")
    record = {"schema": SCHEMA, "scope": SCOPE,
        "input": binding, "runtime": runtime, "assessment": assessment,
        "execution": {"cpu_seconds": _budget(started), "cpu_seconds_max": CPU_LIMIT,
            "field_solves": 0, "wall_solves": 0, "scalar_kernel_values": 0}}
    raw = encoded(record)
    manifest = {"schema": SCHEMA+"-manifest", "source_bindings": sources, "input": binding,
        "outputs": {"assessment.json": {"sha256": digest(raw), "bytes": len(raw)}}}
    directory.mkdir(parents=True, exist_ok=False)
    (directory/"assessment.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return {"created": str(directory), "execution": record["execution"], "summary": _summary(record)}


def check(directory):
    started = time.process_time()
    directory = Path(directory)
    if {p.name for p in directory.iterdir()} != {"assessment.json", "manifest.json"}:
        raise ValueError("assessment file inventory mismatch")
    if any((directory/name).is_symlink() for name in ("assessment.json", "manifest.json")):
        raise ValueError("assessment files must not be symlinks")
    manifest = json.loads((directory/"manifest.json").read_bytes())
    raw = (directory/"assessment.json").read_bytes()
    if (manifest.get("schema") != SCHEMA+"-manifest" or manifest.get("source_bindings") != source_bindings()
            or manifest.get("outputs") != {"assessment.json": {"sha256": digest(raw), "bytes": len(raw)}}):
        raise ValueError("assessment output or source binding mismatch")
    saved, binding = input_snapshot()
    record = json.loads(raw)
    if (manifest.get("input") != binding or record.get("input") != binding
            or record.get("schema") != SCHEMA or record.get("scope") != SCOPE
            or record.get("assessment") != assess(saved)):
        raise ValueError("assessment input or deterministic postprocessing mismatch")
    runtime = record.get("runtime")
    if (not isinstance(runtime, dict) or set(runtime) != {"python", "numpy", "scipy"}
            or any(not isinstance(value, str) or not value or not value[0].isdigit()
                   or "." not in value for value in runtime.values())):
        raise ValueError("assessment runtime metadata mismatch")
    cost = record["execution"]
    if (cost.get("cpu_seconds_max") != CPU_LIMIT or not 0 <= _finite(cost["cpu_seconds"]) <= CPU_LIMIT
            or any(type(cost.get(key)) is not int or cost[key] != 0
                   for key in ("field_solves", "wall_solves", "scalar_kernel_values"))):
        raise ValueError("assessment execution budget mismatch")
    _budget(started)
    return {"checked": str(directory), "physical_solves": 0, "summary": _summary(record)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    print(json.dumps(check(args.output) if args.check else create(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
