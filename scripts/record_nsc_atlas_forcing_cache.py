#!/usr/bin/env python3
"""Create or authenticate a bounded forcing-cache time comparison.

Authentication reconstructs interpolants from retained slices and never
evaluates the prepared field or evolves the scalar correction.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np
from validation_source_successors import verify_saved_bindings

from record_nsc_atlas_scalar_inputs import (
    BASE, SOURCES as INPUT_SOURCES, DOMAIN_KEYS, FORCING_ARRAYS, FORCING_COUNTS,
    FORCING_SCALARS, PROFILE_COUNTS, _prepare, _plain, _real, _count, _equal,
    digest, encoded, runtime_versions,
)

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-forcing-cache-comparison-v1"
STATUS = "FORMAL_FORCING_CACHE_REFINEMENT"
SCOPE = "sampled forcing-cache interpolation sensitivity; no evolved scalar correction or uniform error enclosure"
TRAINING_COUNTS = (17, 33)
TRAINING_TIMES = tuple(i/32 for i in range(33))
WITHHELD = tuple((2*i+1)/64 for i in range(32))
DIRECT_TIMES = tuple(i/64 for i in range(65))
SETTINGS = {**BASE, "eta_nodes": 129, "withheld_eta_nodes": 257,
            "forcing_radial_order": 2, "forcing_field_order": 32,
            "forcing_source_order": 64}
DEFECT_KEYS = ("wall_mass_defect", "source_p_defect", "source_p_V_defect",
               "wall_chi_U", "chi_U_lower_limit_bound")
SOURCES = tuple(dict.fromkeys((*INPUT_SOURCES,
    "docs/nsc-atlas-scalar-evolution.md",
    "src/recursive_horizons/nsc_atlas_forcing_cache.py",
    "tests/test_nsc_atlas_forcing_cache.py",
    "scripts/record_nsc_atlas_forcing_cache.py",
    "tests/test_record_nsc_atlas_forcing_cache.py")))


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def _same(actual, expected, label="saved interpolation"):
    """Compare retained numerical diagnostics without replaying physical data."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(actual) != set(expected):
            raise ValueError(label+" inventory mismatch")
        for key in expected:
            _same(actual[key], expected[key], label+" "+key)
    elif isinstance(expected, (list, tuple)):
        if not isinstance(actual, list) or len(actual) != len(expected):
            raise ValueError(label+" shape mismatch")
        for a, b in zip(actual, expected):
            _same(a, b, label)
    elif isinstance(expected, bool) or expected is None or isinstance(expected, str):
        if type(actual) is not type(expected) or actual != expected:
            raise ValueError(label+" mismatch")
    elif type(expected) is int:
        if type(actual) is not int or actual != expected:
            raise ValueError(label+" count mismatch")
    else:
        _real(actual, label)
        if not _equal(actual, expected):
            raise ValueError(label+" numerical mismatch")


def _solutions(provider, wedge, geometry, budget):
    from recursive_horizons.nsc_atlas_forcing_cache import ForcingTimeCache
    fits = [ForcingTimeCache(wedge, geometry, np.linspace(0, 1, 129),
            np.linspace(0, 1, count), WITHHELD, budget=budget,
            radial_order=2, field_order=32, source_order=64, slice_provider=provider)
            for count in TRAINING_COUNTS]
    results = []
    for fit in fits:
        summary = fit.summary()
        summary.pop("slice_provider")
        summary["validation_samples"] = [
            {key: value for key, value in row.items() if key != "reference"}
            for row in fit.validation_rows]
        results.append(_plain(summary))
    return results


def comparison_diagnostics(solutions):
    from recursive_horizons.nsc_atlas_forcing_cache import DEFECT_KEYS as CACHE_DEFECTS
    keys = (*CACHE_DEFECTS, "scaled_g1_V_defect")
    return {"maximum_absolute_defects": [
            {key: max(abs(row[key]) for row in solution["validation_samples"]) for key in keys}
            for solution in solutions],
        "training_times": list(TRAINING_COUNTS), "withheld_times": 32,
        "temporal_scope": "withheld times at the 129 radial training knots",
        "combined_scope": "withheld times at the 128 radial midpoints",
        "scaled_source_slope": "time-cell V width times g1_V discrepancy",
        "uniform_error_enclosure": False}


def _saved_solutions(record):
    """Use exact retained wall values and slices; missing inputs fail closed."""
    from recursive_horizons.nsc_atlas_bulk_response import BulkBudget
    from recursive_horizons.nsc_atlas_forcing_cache import ForcingSliceProvider
    d = record["domain"]
    by_V = {row["V"]: row["U_wall"] for row in record["direct_samples"]}
    def wall(V):
        values = np.asarray(V)
        try:
            return np.array([by_V[float(value)] for value in values.ravel()]).reshape(values.shape)
        except KeyError as error:
            raise ValueError("authentication requested an unrecorded wall coordinate") from error
    shared = SimpleNamespace(started=time.process_time(), kernel_values=0, check=lambda *_: None)
    metrics = BulkBudget(started=shared.started, external=shared)
    wedge = SimpleNamespace(domain=dict(d), geometry={"wall": wall}, budget=shared)
    rows = [{key: value for key, value in row.items()
             if key not in ("wall_integral_L_plus", "wall_mass_defect")}
            for row in record["direct_samples"]]
    provider = ForcingSliceProvider.from_records(wedge, record["source_geometry"], rows,
        budget=metrics, radial_order=2, field_order=32, source_order=64)
    return _solutions(provider, wedge, record["source_geometry"], metrics)


def compare():
    from recursive_horizons.nsc_atlas_forcing_cache import ForcingSliceProvider
    shared, metrics, geometry, wedge, wall, cap = _prepare()
    provider = ForcingSliceProvider(wedge, geometry, budget=metrics,
        radial_order=2, field_order=32, source_order=64)
    solutions = _solutions(provider, wedge, geometry, metrics)
    rows = []
    for raw in sorted(provider.cache.values(), key=lambda row: row["t"]):
        row = dict(raw)
        row["wall_integral_L_plus"] = wall._wall(row["t"])["integral_L_plus"]
        row["wall_mass_defect"] = float(row["mu_E"][-1])-row["wall_integral_L_plus"]
        rows.append(_plain(row))
    provider_summary = _plain(provider.summary())
    shared.check()
    metrics.check()
    record = {"schema": SCHEMA, "status": STATUS, "scope": SCOPE,
        "source_scale": 1., "runtime": runtime_versions(), "settings": dict(SETTINGS),
        "domain": {key: float(geometry[key]) for key in DOMAIN_KEYS},
        "source_geometry": _plain(geometry), "direct_samples": rows, "solutions": solutions,
        "diagnostics": comparison_diagnostics(solutions), "slice_provider": provider_summary,
        "refinement_is_rigorous_error_enclosure": False, "heldout_used_for_training": False,
        "profile_evaluation": {"counts": dict(cap.evaluation_counts),
            "maximum_inverse_bracket_width": cap.max_inverse_bracket_width},
        "execution": {"cpu_seconds": time.process_time()-shared.started, "cpu_seconds_max": shared.cpu_seconds,
            "scalar_kernel_values": shared.kernel_values, "scalar_kernel_values_max": shared.kernel_values_max,
            "metric_values": metrics.metric_values, "source_quadrature_nodes": metrics.source_nodes,
            "rhs_calls": metrics.rhs_calls, "field_solves": 1, "wall_solves": 1,
            "direct_forcing_slices": len(provider.cache)}}
    record["direct_diagnostics"] = _validate_direct(record)
    return _plain(record)


def _validate_direct(record):
    d, rows = record["domain"], record["direct_samples"]
    if (not isinstance(rows, list) or len(rows) != 65
            or [row["t"] for row in rows] != list(DIRECT_TIMES)):
        raise ValueError("forcing cache direct time inventory mismatch")
    native = np.linspace(d["U0"], d["Ub"], BASE["field_nodes"])
    scalar_keys = set(FORCING_SCALARS)
    for i, row in enumerate(rows):
        n = 129 if i % 2 == 0 else 257
        if set(row) != set(FORCING_ARRAYS+FORCING_SCALARS+FORCING_COUNTS) | {"eta", "slice_costs"}:
            raise ValueError("forcing cache direct slice inventory mismatch")
        _same(row["eta"], np.linspace(0, 1, n).tolist(), "direct eta grid")
        if set(row["slice_costs"]) != {"kernel_values", "metric_values", "source_nodes"}:
            raise ValueError("forcing cache slice cost inventory mismatch")
        for key, value in row["slice_costs"].items():
            _count(value, key, 30_000_000)
        for key in FORCING_ARRAYS:
            if not isinstance(row[key], list) or len(row[key]) != n:
                raise ValueError("forcing cache direct array shape mismatch")
            for value in row[key]:
                _real(value, key)
        for key in scalar_keys:
            _real(row[key], key)
        for key in FORCING_COUNTS:
            _count(row[key], key, 513*5)
        U = np.asarray(row["U"])
        expected = d["U0"]+(row["U_wall"]-d["U0"])*np.linspace(0, 1, n)
        expected[-1] = row["U_wall"]
        if (not d["Ua"]-1e-13 <= row["U_wall"] <= d["Ub"]+1e-13
                or not _equal(row["V"], d["Va"]+(d["Vb"]-d["Va"])*row["t"])
                or not np.allclose(U, expected, rtol=0, atol=2e-15)
                or U[0] != d["U0"] or U[-1] != row["U_wall"]
                or np.any(np.diff(U) <= 0) or any(value <= 0 for value in row["K0"])
                or row["source_r_U"] <= 0 or row["chi_U_lower_limit_bound"] < 0
                or row["wall_chi_U"] != 0 or row["F"][-1] != 0):
            raise ValueError("forcing cache domain or wall identity mismatch")
        cuts = np.unique(np.r_[d["U0"], native[(native > d["U0"]) & (native < row["U_wall"])], U])
        gauss = cuts[:-1, None]+(np.diff(cuts)/2)[:, None]*(1+np.polynomial.legendre.leggauss(2)[0])
        points = len(np.unique(np.r_[cuts, gauss.ravel()]))
        if (row["radial_panels"] != len(cuts)-1 or row["radial_nodes"] != gauss.size
                or row["field_points"] != points or row["field_batches"] != (points+512)//513):
            raise ValueError("forcing cache radial quadrature accounting mismatch")
        g1_V = -(row["source_H_E_V"]+row["source_H_E"]*row["source_P0_r"]/2)*row["source_chi_r"]
        g1_V += row["source_H_E"]*row["K0"][0]*row["source_chi"]/row["source_r_U"]
        if (not _equal(row["g1"], -row["source_H_E"]*row["source_chi_r"])
                or not _equal(row["g1_V"], g1_V)
                or row["mu_E"][0] != row["source_mu_E"] or row["zeta_E"][0] != row["source_zeta_E"]
                or not _equal(row["wall_mass_defect"], row["mu_E"][-1]-row["wall_integral_L_plus"])):
            raise ValueError("forcing cache source or independent mass identity mismatch")
        if row["t"] <= .25 and (any(abs(row[key]) > 1e-24 for key in ("g1", "g1_V", "wall_integral_L_plus"))
                or any(abs(value) > 1e-24 for key in ("F", "G", "mu_E", "zeta_E") for value in row[key])):
            raise ValueError("forcing cache quiet interval mismatch")
    return {"unique_direct_slices": 65, "training_slices": 33, "withheld_slices": 32,
        "field_points": sum(row["field_points"] for row in rows),
        "field_batches": sum(row["field_batches"] for row in rows),
        "radial_nodes": sum(row["radial_nodes"] for row in rows),
        "maximum_absolute_direct_defects": {key: max(abs(row[key]) for row in rows) for key in DEFECT_KEYS}}


def validate_record(record):
    try:
        if (not isinstance(record, dict) or record.get("schema") != SCHEMA
                or record.get("status") != STATUS or record.get("scope") != SCOPE
                or record.get("refinement_is_rigorous_error_enclosure") is not False
                or record.get("heldout_used_for_training") is not False
                or record.get("source_scale") != 1. or isinstance(record.get("source_scale"), bool)):
            raise ValueError("forcing cache declaration mismatch")
        _same(record["settings"], SETTINGS, "fixed settings")
        runtime = record["runtime"]
        if (not isinstance(runtime, dict) or set(runtime) != {"python", "numpy", "scipy"}
                or any(not isinstance(v, str) or not v or not v[0].isdigit() or "." not in v for v in runtime.values())):
            raise ValueError("forcing cache runtime mismatch")
        d = record["domain"]
        if set(d) != DOMAIN_KEYS:
            raise ValueError("forcing cache domain inventory mismatch")
        for value in d.values():
            _real(value, "domain")
        if not (d["U0"] < d["Ua"] < d["Ub"] < 0 and d["Va"] < d["Vb"] < 0 and d["U0"]*d["Va"] < 1):
            raise ValueError("forcing cache exits the retained domain")
        from recursive_horizons.nsc_atlas_characteristic_initialization import validate_geometry
        validate_geometry(record["source_geometry"])
        if any(record["source_geometry"][key] != d[key] for key in d):
            raise ValueError("forcing cache source geometry differs from domain")
        _same(record["direct_diagnostics"], _validate_direct(record), "direct diagnostics")
        solutions = _saved_solutions(record)
        _same(record["solutions"], solutions)
        _same(record["diagnostics"], comparison_diagnostics(solutions), "comparison diagnostics")
        profile = record["profile_evaluation"]
        if set(profile["counts"]) != PROFILE_COUNTS:
            raise ValueError("forcing cache profile counter inventory mismatch")
        for key, value in profile["counts"].items():
            _count(value, key)
        counts = profile["counts"]
        if not counts["second_points"] <= counts["wall_points"] <= counts["profile_requested"] or counts["cache_hits"] > counts["profile_requested"]:
            raise ValueError("forcing cache profile counter mismatch")
        _real(profile["maximum_inverse_bracket_width"], "inverse bracket", True)
        if profile["maximum_inverse_bracket_width"] > 4*np.finfo(float).eps*max(abs(d["Va"]), abs(d["Vb"])):
            raise ValueError("forcing cache inverse bracket mismatch")
        cost = record["execution"]
        _real(cost["cpu_seconds_max"], "CPU ceiling")
        _count(cost["scalar_kernel_values_max"], "kernel ceiling")
        if cost["cpu_seconds_max"] != 120 or cost["scalar_kernel_values_max"] != 30_000_000:
            raise ValueError("forcing cache budget ceilings changed")
        _real(cost["cpu_seconds"], "CPU time", True)
        if cost["cpu_seconds"] > 120:
            raise ValueError("forcing cache CPU budget exceeded")
        for key in ("scalar_kernel_values", "metric_values", "source_quadrature_nodes"):
            _count(cost[key], key, 30_000_000)
        if cost["metric_values"]+cost["source_quadrature_nodes"] > 30_000_000:
            raise ValueError("forcing cache metric budget exceeded")
        _count(cost["rhs_calls"], "RHS calls", 12000)
        for key, value in (("field_solves", 1), ("wall_solves", 1), ("direct_forcing_slices", 65)):
            if type(cost[key]) is not int or cost[key] != value:
                raise ValueError("forcing cache solve or reuse count mismatch")
        provider = record["slice_provider"]
        expected_counts = {key: sum(row[key] for row in record["direct_samples"])
                           for key in ("radial_nodes", "field_points", "field_batches")}
        expected_counts.update({key: sum(row["slice_costs"][key] for row in record["direct_samples"])
                                for key in ("kernel_values", "metric_values", "source_nodes")})
        expected_counts.update(physical_slices=65, cache_hits=49)
        expected_counts["maximum_absolute_input_defects"] = {
            "quiet_input_defect": max(max(abs(row[key]) for key in ("g1", "g1_V"))
                                     for row in record["direct_samples"] if row["t"] <= .25),
            "wall_flux_defect": max(abs(row["F"][-1]) for row in record["direct_samples"]),
            **{key: max(abs(row[key]) for row in record["direct_samples"])
               for key in ("source_p_defect", "source_p_V_defect")}}
        expected_counts["maximum_absolute_input_defects"]["quiet_input_defect"] = max(
            expected_counts["maximum_absolute_input_defects"]["quiet_input_defect"],
            *(abs(value) for row in record["direct_samples"] if row["t"] <= .25 for key in ("F", "G") for value in row[key]))
        _same(provider, expected_counts, "slice reuse accounting")
        if (provider["kernel_values"] > cost["scalar_kernel_values"]
                or provider["metric_values"] > cost["metric_values"]
                or provider["source_nodes"] > cost["source_quadrature_nodes"]):
            raise ValueError("forcing cache cumulative cost is smaller than its slices")
    except (KeyError, TypeError, IndexError) as error:
        raise ValueError("incomplete forcing cache record") from error
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
        "outputs": {"cache.json": {"sha256": digest(raw), "bytes": len(raw)}}}
    directory.mkdir(parents=True, exist_ok=False)
    (directory/"cache.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return {"created": str(directory), "execution": record["execution"],
            "summary": record["diagnostics"]}


def check(directory):
    directory = Path(directory)
    path = directory/"manifest.json"
    manifest = json.loads(path.read_bytes())
    if (path.is_symlink() or manifest.get("schema") != SCHEMA+"-manifest"
            or not verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=SOURCES)):
        raise ValueError("forcing cache source binding or schema mismatch")
    if ({p.name for p in directory.iterdir()} != {"manifest.json", "cache.json"}
            or set(manifest.get("outputs", {})) != {"cache.json"}):
        raise ValueError("forcing cache inventory mismatch")
    path = directory/"cache.json"
    raw = path.read_bytes()
    if path.is_symlink() or manifest["outputs"]["cache.json"] != {"sha256": digest(raw), "bytes": len(raw)}:
        raise ValueError("forcing cache output mismatch")
    validate_record(json.loads(raw))
    return {"checked": str(directory), "physical_solves": 0, "saved_interpolant_reconstructions": 2}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    print(json.dumps(check(args.output) if args.check else create(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
