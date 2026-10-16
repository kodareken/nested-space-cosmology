#!/usr/bin/env python3
"""Bounded saved-state transport audit; no source or time evolution replay.

The cap oracle is a pullback of the saved zero-amplitude wall trace. Crossed
exterior evaluations hold the saved state and cap metric fixed. Neither is a
new solution or an error enclosure.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.interpolate import PchipInterpolator

import record_nsc_atlas_compatible_evolution as producer
from validation_source_successors import verify_saved_source
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_nonlinear_boundary import hybrid_wall_boundary
from recursive_horizons.nsc_atlas_reduced_cap import reduced_cap_metric
from recursive_horizons.nsc_atlas_vacuum_departure import vacuum_slice, departure_slice

ROOT = Path(__file__).resolve().parents[1]
INPUT = producer.RELATIVE_OUTPUT + "/compatible.json"
INPUT_SHA256 = "a7e27b25ce171951415f97b64e055fd93062e3c21cd709c19c99c7543167c434"
INPUT_MANIFEST_SHA256 = "7a8b5b6215dacb7a554c2c9b60e14f698217d944859604312775b9aeb48d3ff7"
SCHEMA = "nsc-atlas-compatible-transport-audit-v1"
LIMITS = dict(cpu_seconds=30., evaluated_points=100000)
CUTS = (16, 69, 79, 128)
FIELDS = ("r", "L", "c")
OWN_SOURCES = ("scripts/diagnose_nsc_atlas_compatible_transport.py",
               "tests/test_diagnose_nsc_atlas_compatible_transport.py")
plain, encoded, digest, statistics = (producer.plain, producer.encoded,
                                    producer.digest, producer.statistics)


def authenticate(root=ROOT):
    """Check fixed payload, manifest and all producer bindings; no validators."""
    root = Path(root).resolve()
    path = root / INPUT
    manifest_path = path.with_name("manifest.json")
    raw, manifest_raw = path.read_bytes(), manifest_path.read_bytes()
    manifest = json.loads(manifest_raw)
    bindings = manifest.get("source_bindings", {})
    if (digest(raw) != INPUT_SHA256
            or digest(manifest_raw) != INPUT_MANIFEST_SHA256
            or manifest.get("schema") != producer.SCHEMA + "-manifest"
            or manifest.get("outputs") != {"compatible.json": dict(sha256=digest(raw), bytes=len(raw))}
            or set(bindings) != set(producer.SOURCES) or len(bindings) != 51):
        raise ValueError("fixed compatible input/manifest mismatch")
    retained = {}
    total_bytes = len(raw) + len(manifest_raw)
    for name, expected in bindings.items():
        candidate = (root / name).resolve()
        if not candidate.is_relative_to(root):
            raise ValueError("input binding exits the repository")
        data = verify_saved_source(name,expected,root)
        total_bytes += len(data)
        if name == producer.previous.RELATIVE_OUTPUT + "/coupled.json":
            retained[name] = data
    record = json.loads(raw)
    parent_name = producer.previous.RELATIVE_OUTPUT + "/coupled.json"
    parent = json.loads(retained[parent_name])
    if (record.get("schema") != producer.SCHEMA
            or record.get("settings") != producer.SETTINGS
            or record["parent_input"]["sha256"] != bindings[parent_name]
            or record["parent_input"]["path"] != parent_name):
        raise ValueError("fixed input identity/settings mismatch")
    identity = dict(path=INPUT, sha256=digest(raw), manifest_sha256=digest(manifest_raw),
                    authenticated_source_files=len(bindings), hashed_bytes=total_bytes,
                    recursive_validation=False)
    return record, parent, identity


def departure_jets(values, U, method):
    """Differentiate the stored departures, preserving the analytic vacuum."""
    U = np.asarray(U, dtype=float)
    if (U.ndim != 1 or len(U) < 3 or not np.all(np.isfinite(U))
            or np.any(np.diff(U) <= 0)):
        raise ValueError("require at least three increasing finite coordinates")
    if method not in ("first", "second"):
        raise ValueError("derivative method must be first or second")
    output = {}
    for key, name in (("r", "r_Ubar"), ("L", "log_C_Ubar"), ("c", "c_Ubar")):
        z = np.asarray(values[key], dtype=float)
        if z.shape != U.shape or not np.all(np.isfinite(z)):
            raise ValueError("departure shape/finite mismatch")
        output[name] = (np.r_[np.diff(z)/np.diff(U), (z[-1]-z[-2])/(U[-1]-U[-2])]
                        if method == "first" else np.gradient(z, U, edge_order=2))
    return output


def cap_pullback(case, budget):
    """Compare the full retained cap primitive to its incoming-wave oracle."""
    if case["alpha"] != 0 or case["amplitude"] != "zero":
        raise ValueError("incoming oracle requires the saved zero-amplitude cap")
    samples, n = case["snapshots"], case["cells"]
    if not 3 <= len(samples) <= 129 or not 3 <= n <= 128:
        raise ValueError("cap audit exceeds the saved finite grid")
    T = np.array([s["physical"]["wall"]["T"] for s in samples])
    R = np.array([s["physical"]["wall"]["r"] for s in samples])
    cW = np.array([s["physical"]["wall"]["c"] for s in samples])
    if (not np.all(np.isfinite(np.r_[T, R, cW])) or np.any(R <= 0)
            or np.any(R >= 1)):
        raise ValueError("cap oracle exits the static finite-radius domain")
    qW = T + np.arctanh(R)
    if np.any(np.diff(qW) <= 0) or cW[0] != 0 or np.max(T) >= qW[0]:
        raise ValueError("require increasing arrival labels and quiet reflected branch")
    budget.check(3*len(samples))
    pchip = PchipInterpolator(qW, cW, extrapolate=False)
    maxima = {k: dict(absolute=-1.) for k in ("linear", "pchip", "interpolation_sensitivity")}
    rows = []
    for i, sample in enumerate(samples):
        radius = R[i]*np.linspace(0., 1., n+1)
        q = T[i] + np.arctanh(radius)
        primitive = np.r_[0., np.cumsum(np.asarray(sample["state"]["cap"]["Iq"], dtype=float))]
        if primitive.shape != q.shape or not np.all(np.isfinite(primitive)):
            raise ValueError("invalid saved cap primitive")
        if not np.isclose(primitive[-1], cW[i], rtol=1e-13, atol=1e-14):
            raise ValueError("zero cap primitive does not match its wall trace")
        if q[-1] != qW[i] or np.max(q) > qW[-1]:
            raise ValueError("incoming pullback needs an unrecorded future wall value")
        active = q >= qW[0]
        linear, cubic = np.zeros(n+1), np.zeros(n+1)
        linear[active] = np.interp(q[active], qW, cW)
        cubic[active] = pchip(q[active])
        errors = dict(linear=primitive-linear, pchip=primitive-cubic,
                      interpolation_sensitivity=cubic-linear)
        budget.check(6*(n+1))
        row = dict(snapshot=i, x=sample["x"], T=float(T[i]), R=float(R[i]),
                   q_wall=float(qW[i]), quiet_points=int(np.count_nonzero(~active)))
        for key, difference in errors.items():
            j = int(np.argmax(abs(difference)))
            row[key] = statistics(difference)
            if abs(difference[j]) > maxima[key]["absolute"]:
                maxima[key] = dict(absolute=float(abs(difference[j])), signed=float(difference[j]),
                    snapshot=i, face=j, x=sample["x"], y=j/n, r=float(radius[j]), q=float(q[j]))
        rows.append(row)
    return dict(cells=n, steps=case["steps"], q_wall=[float(qW[0]), float(qW[-1])],
        reflected_branch_upper_label=float(np.max(T)), rows=rows, maxima=maxima,
        interpretation="saved zero wall trace pulled back with the exact de Sitter incoming characteristics; "
                       "linear/PCHIP sensitivity is sampled, not a uniform interpolation bound")


def wall_readout(prefix, wall, gauge, cap_wall, alpha, budget):
    """Same algebraic wall map for each frozen-state prefix; no wall ODE."""
    outgoing = {"r_V": float(prefix["r_x"][-1]), "c_V": float(prefix["c_x"][-1]),
                "log_C_V": float(prefix["log_C_x"][-1])}
    common = dict(r=wall["r"], c=wall["c"], eta_plus=gauge["eta_plus"])
    C = float(np.exp(wall["L"]))
    first = hybrid_wall_boundary(common, outgoing, cap_wall,
                                 dict(C_wall=C, log_C_s=0.), alpha=alpha)
    xs = float(first["Vdot_plus"])
    Lx = 2*gauge["eta_plus_x"]-2*float(first["k_plus"])/xs+2*outgoing["log_C_V"]
    second = hybrid_wall_boundary(common, outgoing, cap_wall,
                                  dict(C_wall=C, log_C_s=Lx*xs), alpha=alpha)
    rates = dict(r=float(second["r_s"])/xs, L=float(Lx), c=float(second["c_s"])/xs)
    Q = dict(r=outgoing["r_V"], L=outgoing["log_C_V"], c=outgoing["c_V"])
    P = {key: float(prefix[name][-1]) for key, name in
         (("r", "r_Ubar"), ("L", "log_C_Ubar"), ("c", "c_Ubar"))}
    required = {key: (rates[key]-Q[key])/gauge["W_x"] for key in FIELDS}
    residual = {key: rates[key]-Q[key]-gauge["W_x"]*P[key] for key in FIELDS}
    budget.check(2)
    return dict(longitudinal=Q, transverse=P, required_transverse=required,
                wall_rates=rates, wall_kinematic_defect=residual,
                incoming_w_plus=float(second["w_plus"]), x_s=xs)


def crossed_stage(case, parent, index, budget):
    """One fixed state, one cap radial reconstruction, four exterior prefixes."""
    sample, n = case["snapshots"][index], case["cells"]
    alpha, x = float(case["alpha"]), float(sample["x"])
    gauge = producer.subset(parent["gauge"]["rows"], 4*index)
    source = producer.subset(parent["initial_data"][case["amplitude"]]["snapshot_rows"], 4*index)
    if gauge["x"] != x or source["x"] != x:
        raise ValueError("selected stage is absent from the exact input cache")
    eta, state = np.linspace(0., 1., n+1), sample["state"]
    wall = {key: float(gauge[key])+float(state["wall"][key]) for key in ("r", "L", "T", "s")}
    metric = reduced_cap_metric(eta, wall["r"], state["cap"], alpha=alpha,
                                quadrature_order=4, budget=budget)
    wall["c"] = float(metric["wall"]["c"])
    U = eta*float(gauge["W"])
    vacuum = vacuum_slice(U, x, budget=budget)
    values, Z = {}, {}
    for key in FIELDS:
        base = vacuum[key] if key != "c" else np.zeros_like(U)
        src = float(source["log_C" if key == "L" else key])
        offset = src-base[0]
        Z[key] = np.r_[0., state["departure"][key], wall[key]-base[-1]-offset]
        values[key] = base+offset+Z[key]
        values[key][0], values[key][-1] = src, wall[key]
        np.testing.assert_allclose(values[key][1:-1]-src,
            sample["physical"]["increment"][key], rtol=1e-12, atol=1e-13)
    for key in wall:
        if not np.isclose(wall[key], sample["physical"]["wall"][key], rtol=1e-12, atol=1e-13):
            raise ValueError("selected physical wall does not reconstruct from its raw state")
    cap_wall = {key: metric["wall"][key] for key in ("A", "delta", "w_minus")}
    results, prefixes = {}, {}
    for method in ("first", "second"):
        jets = departure_jets(Z, U, method)
        budget.check(3*(n+1))
        for order in (2, 4):
            name = f"{method}_{order}"
            prefix = departure_slice(U, values, jets, source, x=x, alpha=alpha,
                                      tube=producer.TUBE, budget=budget, quadrature_order=order)
            results[name] = wall_readout(prefix, wall, gauge, cap_wall, alpha, budget)
            prefixes[name] = {key: prefix[name] for key, name in
                             (("r", "r_x"), ("L", "log_C_x"), ("c", "c_x"))}
    comparisons = {}
    for label, left, right in (("quadrature_first", "first_4", "first_2"),
                               ("quadrature_second", "second_4", "second_2"),
                               ("derivative_order2", "second_2", "first_2"),
                               ("derivative_order4", "second_4", "first_4")):
        comparisons[label] = dict(
            longitudinal={key: statistics(prefixes[left][key]-prefixes[right][key]) for key in FIELDS},
            wall_kinematic_change={key: results[left]["wall_kinematic_defect"][key]
                -results[right]["wall_kinematic_defect"][key] for key in FIELDS})
        budget.check(3*(n+1))
    radial_identity = {}
    for order in (2, 4):
        a, b = results[f"first_{order}"], results[f"second_{order}"]
        radial_identity[str(order)] = (b["wall_kinematic_defect"]["r"]-a["wall_kinematic_defect"]["r"]
            +float(gauge["W_x"])*(b["transverse"]["r"]-a["transverse"]["r"]))
    return plain(dict(snapshot=index, x=x, cells=n, alpha=alpha, wall=wall,
        cap_wall=cap_wall, cap_radial_reconstructions=1, results=results,
        comparisons=comparisons, radial_jet_change_identity=radial_identity,
        interpretation="fixed saved physical state and cap reconstruction; endpoint slope sensitivity "
                       "and prefix sensitivity do not establish the accuracy of the cap state or a new evolution"))


def audit(root=ROOT, *, budget=None):
    if budget is None:
        budget = Budget(maximum_evaluations=LIMITS["evaluated_points"], cpu_seconds=LIMITS["cpu_seconds"])
    if (not isinstance(budget, Budget)
            or type(budget.maximum_evaluations) is not int
            or not 0 < budget.maximum_evaluations <= LIMITS["evaluated_points"]
            or isinstance(budget.cpu_seconds, bool)
            or not math.isfinite(budget.cpu_seconds)
            or not 0 < budget.cpu_seconds <= LIMITS["cpu_seconds"]
            or type(budget.evaluations) is not int or budget.evaluations < 0
            or not math.isfinite(budget.started) or budget.started > time.process_time()):
        raise ValueError("audit caller allowance exceeds the fixed limits")
    budget.check()
    record, parent, identity = authenticate(root)
    own = {name: digest((Path(root)/name).read_bytes()) for name in OWN_SOURCES}
    before = budget.evaluations
    caps = [cap_pullback(record["cases"][i], budget) for i in (0, 1)]
    cap_points = budget.evaluations-before
    before = budget.evaluations
    stages = [crossed_stage(record["cases"][3], parent, i, budget) for i in CUTS]
    stage_points = budget.evaluations-before
    # Catch any change during calculation before returning or writing a result.
    _, _, after = authenticate(root)
    if after != identity or own != {name: digest((Path(root)/name).read_bytes()) for name in OWN_SOURCES}:
        raise ValueError("audit input/source changed during calculation")
    budget.check()
    return plain(dict(schema=SCHEMA, input=identity, audit_source_bindings=own,
        runtime=producer.runtime_versions(), limits=LIMITS,
        settings=dict(zero_cases=[0, 1], exterior_case=3, exterior_cuts=CUTS,
                      cap_quadrature_order=4, prefix_orders=[2, 4],
                      derivative_methods=["first", "second"],
                      second_derivative_method="numpy.gradient(edge_order=2) on departure values"),
        cap_pullback=caps, crossed_stages=stages,
        execution=dict(evaluated_points=budget.evaluations, cap_pullback_points=cap_points,
            crossed_stage_points=stage_points, cpu_seconds=time.process_time()-budget.started,
            hash_authentication_passes=2,
            source_integrations=0, time_steps=0, recursive_validations=0,
            cap_radial_reconstructions=len(CUTS), exterior_prefix_evaluations=4*len(CUTS)),
        scope="saved-state mechanism audit; no new physical trajectory, convergence result or uniform error bound"))


def create(path=None, *, root=ROOT):
    if path is not None and Path(path).exists():
        raise FileExistsError("audit output is creation-only")
    result = audit(root)
    if path is not None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(encoded(result))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="create one JSON file; otherwise print JSON")
    args = parser.parse_args()
    print(json.dumps(create(args.output), sort_keys=True, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
