#!/usr/bin/env python3
"""Locate the boundary interpolation discrepancy using authenticated saved jets."""
import argparse
import hashlib
import json
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np

import record_nsc_atlas_local_lift as owner
from recursive_horizons.nsc_atlas_adaptive_lift import AdaptiveBoundaryLift

ROOT = Path(__file__).resolve().parents[1]
INPUT = "results/development/nsc-atlas-local-lift-v1"
SCHEMA = "nsc-atlas-lift-component-assessment-v1"
SOURCES = (INPUT+"/lift.json", INPUT+"/manifest.json",
           "scripts/assess_nsc_atlas_lift_components.py")
NAMES = ("backscatter", "material_radius", "clock")


def bindings():
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES}


def component_jets(saved):
    """Product-rule decomposition; preserve the supplied continuum first jets."""
    r, rs, e, es, E, Es, C, Cs = [saved[key] for key in (
        "r_U", "r_U_s", "e", "e_s", "displacement_at_fixed_V",
        "displacement_at_fixed_V_s", "c0_rr", "c0_rr_s")]
    qm, qmss, d, ds, J, Js, Vs, UV = [saved[key] for key in (
        "q_minus_s", "q_minus_ss", "d", "d_s", "J_minus", "J_minus_s", "V0_s", "U_V")]
    a, a_s = 2*qm/d, 2*qmss/d-2*qm*ds/d**2
    values = {
        "backscatter": (r*a*J, rs*a*J+r*(a_s*J+a*Js)),
        "material_radius": (-r*e*C, -rs*e*C-r*es*C-r*e*Cs),
        "clock": (r*(e-E)*C, rs*(e-E)*C+r*(es-Es)*C+r*(e-E)*Cs),
        "total": (saved["b"], saved["b_s"]),
    }
    return {name: dict(saved, b=b, b_s=bs, b_V=bs/Vs, b_U=bs/(Vs*UV))
            for name, (b, bs) in values.items()}


def assess():
    started = time.process_time()
    def budget_check():
        if time.process_time()-started > 30:
            raise RuntimeError("saved component assessment exceeded 30 process-CPU seconds")
    owner.check(ROOT/INPUT)
    record = json.loads((ROOT/INPUT/"lift.json").read_bytes())
    components = {key: {} for key in (*NAMES, "total")}
    closure = {key: 0. for key in ("b", "b_s", "b_U")}
    for saved in record["direct_samples"]:
        rows = component_jets(saved)
        for name, row in rows.items():
            components[name][saved["t"]] = row
        for key in closure:
            closure[key] = max(closure[key], abs(sum(rows[name][key] for name in NAMES)-saved[key]))
    byV = {row["V"]: row["U"] for row in record["direct_samples"]}
    def wall(V):
        values = np.asarray(V)
        return np.array([byV[float(v)] for v in values.ravel()]).reshape(values.shape)
    fits, summaries = {}, []
    for active in owner.ACTIVE_CELLS:
        for name, data in components.items():
            budget_check()
            jets = SimpleNamespace(domain=record["domain"], budget=SimpleNamespace(check=budget_check),
                geometry={"wall": wall}, cache=data, sample=lambda t, data=data: dict(data[t]))
            fit = AdaptiveBoundaryLift(jets, knot_times=owner.knot_times(active), nodes=5,
                                      validation_times=owner.WITHHELD)
            fits[active, name] = fit
            worst = max(fit.validation_rows, key=lambda row: abs(row["b_defect"]))
            slope = max(fit.validation_rows, key=lambda row: abs(row["scaled_slope_defect"]))
            summaries.append({"active_cells": active, "component": name,
                "maximum_absolute_value_defect": abs(worst["b_defect"]), "value_peak_t": worst["t"],
                "peak_signed_value_defect": worst["b_defect"],
                "maximum_absolute_scaled_slope_defect": abs(slope["scaled_slope_defect"]),
                "slope_peak_t": slope["t"],
                "sampled_b_target_met": abs(worst["b_defect"]) <= record["diagnostics"]["interpolation_target"]["absolute_b_target"]})
    fit_closure = {str(active): max(abs(sum(fits[active, name].validation_rows[i]["b_defect"]
        for name in NAMES)-fits[active, "total"].validation_rows[i]["b_defect"])
        for i in range(len(owner.WITHHELD))) for active in owner.ACTIVE_CELLS}
    if max(closure.values()) > 1e-10 or max(fit_closure.values()) > 1e-11:
        raise ValueError("component decomposition failed value/derivative sum closure")
    budget_check()
    return {"schema": SCHEMA, "scope": "saved-jet interpolation diagnosis; no continuum clock error bound",
        "input": INPUT, "source_bindings": bindings(), "direct_jet_closure": closure,
        "fit_error_sum_closure": fit_closure, "components": summaries,
        "physical_solves": 0, "scalar_kernel_values": 0, "cpu_seconds_max": 30,
        "cpu_seconds": time.process_time()-started}


def create(path):
    path = Path(path)
    if path.exists():
        raise FileExistsError("component assessment is creation-only")
    before = bindings()
    record = assess()
    if before != bindings() or record["source_bindings"] != before:
        raise RuntimeError("component assessment sources changed during calculation")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write((json.dumps(record, indent=2, sort_keys=True, allow_nan=False)+"\n").encode())
    return {"created": str(path), "cpu_seconds": record["cpu_seconds"], "components": record["components"]}


def check(path):
    path = Path(path)
    if path.is_symlink():
        raise ValueError("component assessment cannot be a symlink")
    saved = json.loads(path.read_bytes())
    elapsed = saved.pop("cpu_seconds")
    if isinstance(elapsed, bool) or not isinstance(elapsed, (int, float)) or not 0 <= elapsed <= 30:
        raise ValueError("invalid component assessment budget")
    rebuilt = assess()
    rebuilt.pop("cpu_seconds")
    if not owner._same_saved(saved, rebuilt):
        raise ValueError("saved component assessment differs from its authenticated recomputation")
    return {"checked": str(path), "physical_solves": 0, "saved_component_fits": 8}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    print(json.dumps(check(args.output) if args.check else create(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
