#!/usr/bin/env python3
"""Constraint transport from frozen arrays only; no physical operator replay."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time
import uuid

import numpy as np

import record_nsc_atlas_ssprk3_refinement as producer
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget

ROOT = Path(__file__).resolve().parents[1]
INPUT = producer.RELATIVE_OUTPUT + "/refinement.json"
INPUT_SHA256 = "ab95bdcd43b6df78e776f5282387e3f7398d8996a25d4e3de4e5a8226fc25e02"
MANIFEST_SHA256 = "f85a1491eeba5c77514936e48f69a4335d5a92edda2a0db2a0dc624a3701f041"
BASE_INPUT = "results/development/nsc-atlas-coupled-evolution-v1/coupled.json"
SCHEMA = "nsc-atlas-constraint-transport-audit-v1"
LIMITS = dict(cpu_seconds=30., evaluated_points=1000000)
OWN_SOURCES = ("scripts/diagnose_nsc_atlas_constraint_transport.py",
               "tests/test_diagnose_nsc_atlas_constraint_transport.py")
FIELDS = ("E_x", "inherited", "residual", "transport", "E_Ubar")
plain, encoded, digest = producer.plain, producer.encoded, producer.digest
runtime_versions = producer.runtime_versions
SCOPE = ("saved right-trace source constraint and centred departure constraints; "
         "weighted transport residual and diagnostic differentiation sensitivity, "
         "not continuous evolution error, a convergence order or an error enclosure")


def shared_budget(value=None):
    if value is None:
        value = Budget(maximum_evaluations=LIMITS["evaluated_points"], cpu_seconds=30.)
    if (not isinstance(value, Budget) or type(value.maximum_evaluations) is not int
            or not 0 < value.maximum_evaluations <= LIMITS["evaluated_points"]
            or type(value.evaluations) is not int or value.evaluations < 0
            or isinstance(value.cpu_seconds, bool) or not math.isfinite(value.cpu_seconds)
            or not 0 < value.cpu_seconds <= 30. or isinstance(value.started,bool)
            or not math.isfinite(value.started) or value.started < 0
            or value.started > time.process_time()):
        raise ValueError("require a live shared allowance within30CPU/1M points")
    value.check()
    return value


def authenticate(root, budget):
    """Hash the fixed payload, manifest and86 bindings; call no validator."""
    budget.check()
    root = Path(root).resolve()
    path = root / INPUT
    if ({p.name for p in path.parent.iterdir()} != {"refinement.json", "manifest.json"}
            or any(p.is_symlink() for p in path.parent.iterdir())):
        raise ValueError("fixed refinement inventory mismatch")
    raw, mr = path.read_bytes(), path.with_name("manifest.json").read_bytes()
    manifest = json.loads(mr)
    bindings = manifest.get("source_bindings", {})
    if (digest(raw) != INPUT_SHA256 or digest(mr) != MANIFEST_SHA256
            or manifest.get("schema") != producer.SCHEMA + "-manifest"
            or manifest.get("outputs") != {"refinement.json": dict(sha256=digest(raw), bytes=len(raw))}
            or set(bindings) != set(producer.SOURCES) or len(bindings) != 86):
        raise ValueError("pinned refinement payload/manifest mismatch")
    base_raw, hashed_bytes = None, len(raw) + len(mr)
    for name, expected in bindings.items():
        candidate = (root / name).resolve()
        if not candidate.is_relative_to(root):
            raise ValueError("bound input exits repository")
        data = candidate.read_bytes()
        hashed_bytes += len(data)
        if digest(data) != expected:
            raise ValueError("bound input changed: " + name)
        if name == BASE_INPUT:
            base_raw = data
    record = json.loads(raw)
    if (record.get("schema") != producer.SCHEMA or record.get("settings") != producer.SETTINGS
            or record.get("runtime") != runtime_versions() or base_raw is None):
        raise ValueError("fixed refinement declaration or source ancestor mismatch")
    base = json.loads(base_raw)
    knots = {label: np.asarray(base["initial_data"][label]["exact_rows"]["x"])[::2]
             for label in ("zero", "finite")}
    budget.check()
    return record, knots, dict(path=INPUT, sha256=digest(raw), bytes=len(raw),
        manifest_sha256=digest(mr), authenticated_input_files=86, hashed_bytes=hashed_bytes,
        recursive_scientific_validation=False)


def array(value, shape, name, budget):
    raw = np.asarray(value)
    if (raw.shape != shape or raw.dtype.kind not in "iuf"
            or (isinstance(value, (list, tuple)) and any(
                isinstance(v, (bool, np.bool_)) for v in np.asarray(value, dtype=object).flat))):
        raise ValueError(name + " has invalid real array shape")
    budget.check(raw.size)
    result = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(result)):
        raise ValueError(name + " contains nonfinite values")
    return result


def same(a, b, name):
    if not np.array_equal(np.asarray(a), np.asarray(b)):
        raise ValueError(name + " does not match exactly")


def close(a, b, name):
    if not np.allclose(a, b, rtol=2e-13, atol=3e-14):
        raise ValueError(name + " exceeds arithmetic roundoff tolerance")


def source_constraint(source, alpha, budget):
    """Independent source arithmetic, using actual cached curvature, not its ODE target."""
    if isinstance(alpha, bool) or not math.isfinite(alpha) or alpha < 0:
        raise ValueError("alpha must be finite and nonnegative")
    n = len(source["r"])
    values = {k: array(source[k], (n,), "source." + k, budget)
              for k in ("r", "r_x", "r_xx", "log_C_x", "c", "c_x")}
    if np.any(values["r"] <= 0):
        raise ValueError("source radius must be positive")
    budget.check(5*n)
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        r, rx, c = values["r"], values["r_x"], values["c"]
        u = c/r
        ux = (values["c_x"]-u*rx)/r
        result = values["r_xx"]-values["log_C_x"]*rx + alpha*r*ux*ux/(2*(1+alpha*u*u/6))
    if not np.all(np.isfinite(result)):
        raise ValueError("source constraint arithmetic is nonfinite")
    return result


def summary(values, mask, x, eta, U, budget):
    budget.check(values.size)
    indices = np.flatnonzero(mask)
    if not len(indices):
        return dict(count=0, minimum=None, maximum=None, maximum_absolute=None, rms=None, location=None)
    take = values.ravel()[indices]
    flat = int(indices[np.argmax(abs(take))])
    i, j = np.unravel_index(flat, values.shape)
    return dict(count=len(indices), minimum=float(np.min(take)), maximum=float(np.max(take)),
        maximum_absolute=float(np.max(abs(take))), rms=float(np.sqrt(np.mean(take*take))),
        location=dict(interior_index=[int(i), int(j)], time_grid_index=int(i+1),
            eta_grid_index=int(j+1), x=float(x[i]), eta=float(eta[j]), Ubar=float(U[i,j]),
            signed=float(values[i,j])))


def summaries(fields, masks, x, eta, U, budget):
    return {key: {name: summary(value, mask, x, eta, U, budget) for name, mask in masks.items()}
            for key, value in fields.items()}


def decompose_grid(x, eta, W, radius, E_x, E_Ubar, source_r, source_E, *, fronts, W0, budget):
    """Pure finite-array connection; no wall-born incoming constraint is supplied."""
    nt, ne = len(x), len(eta)
    if not 4 <= nt <= 129 or not 4 <= ne <= 33:
        raise ValueError("diagnostic exceeds retained grid bounds")
    x, eta = array(x, (nt,), "x", budget), array(eta, (ne,), "eta", budget)
    for values, name in ((x, "x"), (eta, "eta")):
        if (np.any(np.diff(values) <= 0)
                or not np.allclose(np.diff(values), np.diff(values)[0], rtol=1e-10, atol=1e-14)):
            raise ValueError(name + " must be uniform and increasing")
    if eta[0] != 0 or eta[-1] != 1 or not math.isfinite(W0) or W0 <= 0:
        raise ValueError("invalid original eta interval or separator")
    W = array(W, (nt,), "W", budget)
    source_r, source_E = (array(v, (nt,), k, budget) for v,k in ((source_r,"source_r"),(source_E,"source_E")))
    radius = array(radius, (nt, ne-2), "radius", budget)
    shape = (nt-2, ne-2)
    Ex, EU = array(E_x, shape, "E_x", budget), array(E_Ubar, shape, "E_Ubar", budget)
    if np.any(W <= 0) or np.any(source_r <= 0) or np.any(radius <= 0):
        raise ValueError("radius and gauge width must be positive")
    if np.any(np.diff(W) <= 0) or W[0] != W0:
        raise ValueError("require the original strictly expanding wall gauge and separator")
    front = array(fronts, (2,), "source fronts", budget)
    if not x[0] <= front[0] < front[1] <= x[-1]:
        raise ValueError("source fronts must be ordered in the original interval")
    r = radius[1:-1]
    weighted = source_r[1:-1]*source_E[1:-1]
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        inherited = weighted[:,None]/r
        fields = dict(E_x=Ex, inherited=inherited, residual=Ex-inherited,
                      transport=r*Ex-weighted[:,None], E_Ubar=EU)
    if not all(np.all(np.isfinite(v)) for v in fields.values()):
        raise ValueError("nonfinite transport arithmetic")
    budget.check(5*Ex.size+nt)
    U = W[1:-1,None]*eta[None,1:-1]
    bottom = W[1:-1,None]*eta[None,2:] <= W0
    wall = W[1:-1,None]*eta[None,:-2] > W0
    separator = ~(bottom | wall)
    near = np.zeros(nt-2, dtype=bool)
    for at in front:
        near |= (x[:-2] <= at) & (at <= x[2:])
    near = np.broadcast_to(near[:,None], shape).copy()
    masks = dict(all_interior=np.ones(shape,dtype=bool), front_near=near,
                 away_from_front=~near, bottom_origin=bottom, wall_origin=wall, separator=separator)
    budget.check(5*Ex.size)
    return dict(time_nodes=nt, eta_nodes=ne, x=x, eta=eta, W=W,
        interior_x=x[1:-1], interior_eta=eta[1:-1], interior_Ubar=U,
        radius=r, source_r=source_r, source_E=source_E, source_weighted=source_r*source_E,
        arrays=fields, masks=masks,
        summaries=summaries(fields,masks,x[1:-1],eta[1:-1],U,budget),
        location_indices="time_grid_index and eta_grid_index refer to this full saved grid",
        incoming_constraint=dict(bottom_expected=0.,wall_expected="UNKNOWN",
            stencil="all three eta nodes at this time; equality at lower separator is not wall-origin"))


def compare_sampling(fine, coarse, budget):
    """Compare identical points with the UNION of the two time-front masks."""
    same(fine["eta"],coarse["eta"],"matched eta nodes")
    target, source = coarse["interior_x"], fine["interior_x"]
    indices = np.searchsorted(source,target)
    if np.any(indices >= len(source)):
        raise ValueError("coarse interior time not present in fine history")
    same(source[indices],target,"matched interior times")
    same(fine["W"][1:-1][indices],coarse["W"][1:-1],"matched wall gauge")
    same(fine["interior_Ubar"][indices],coarse["interior_Ubar"],"matched Ubar coordinates")
    same(fine["radius"][indices],coarse["radius"],"matched physical radius")
    same(fine["source_weighted"][1:-1][indices],coarse["source_weighted"][1:-1],"matched source baseline")
    for key in ("bottom_origin","wall_origin","separator"):
        same(fine["masks"][key][indices],coarse["masks"][key],"matched origin mask")
    near = fine["masks"]["front_near"][indices] | coarse["masks"]["front_near"]
    masks = dict(coarse["masks"],front_near=near,away_from_front=~near)
    difference = {k:coarse["arrays"][k]-fine["arrays"][k][indices] for k in FIELDS}
    if not np.array_equal(difference["E_Ubar"],np.zeros_like(difference["E_Ubar"])):
        raise ValueError("E_Ubar changed despite identical saved spatial stencil")
    budget.check(7*coarse["arrays"]["E_x"].size)
    return dict(fine_interior_indices=indices,coarse_interior_x=target,
        difference_coarse_minus_fine=difference,union_masks=masks,
        summaries=summaries(difference,masks,target,coarse["interior_eta"],coarse["interior_Ubar"],budget),
        E_Ubar_same_spatial_stencil_exact=bool(np.array_equal(difference["E_Ubar"],np.zeros_like(difference["E_Ubar"]))),
        interpretation="diagnostic differentiation sensitivity; E_Ubar equality is not spatial refinement or a zero wall constraint")


def analyze(record, knots, budget):
    cases = record["cases"]
    if len(cases) != 2:
        raise ValueError("both zero/finite cases required")
    xcache = array(record["cache"]["x_nodes"],(1025,),"cache x",budget)
    gauge = record["cache"]["gauge_rows"]
    same(gauge["x"],xcache,"gauge cache coordinates")
    W = array(gauge["W"],(1025,),"cached W",budget)[::8]
    x, eta = xcache[::8], np.linspace(0,1,33)
    result = []
    for ci,(case,label) in enumerate(zip(cases,("zero","finite"))):
        if (case["amplitude"],case["cells"],case["steps"]) != (label,32,512):
            raise ValueError("require original ordered32-cell512-step cases")
        if len(case["snapshots"]) != 129:
            raise ValueError("require129 saved snapshots")
        same([s["x"] for s in case["snapshots"]],x,"saved cut coordinates")
        same(knots[label],x,"saved cuts on the original129 Hermite knots")
        cached = record["cache"]["source_rows"][label]
        same(cached["x"],xcache,"source cache coordinates")
        source = {k:array(cached[k],(1025,),"cached source."+k,budget)[::8]
                  for k in ("r","r_x","r_xx","log_C_x","c","c_x","source_constraint")}
        derived = source_constraint(source,case["alpha"],budget)
        close(derived,source["source_constraint"],"independent source constraint")
        history = record["diagnostics"]["histories"][ci]
        Es = array(history["departure129"]["source_constraint"],(129,),"saved source E",budget)
        close(Es,derived,"history/source jet constraint")
        increments = array([s["physical"]["increment"]["r"] for s in case["snapshots"]],(129,31),"physical radius increments",budget)
        radius = source["r"][:,None]+increments
        budget.check(radius.size)
        grids = {}
        for nt,stride in ((129,1),(65,2)):
            saved = history["departure"+str(nt)]
            if (saved["settings"]["time_nodes"],saved["settings"]["eta_nodes"],saved["settings"]["analytic_vacuum_jets"]) != (nt,33,True):
                raise ValueError("saved diagnostic grid/method mismatch")
            original_E = array(saved["source_constraint"],(nt,),"grid source E",budget)
            same(original_E,Es[::stride],"same source residual under time coarsening")
            grids[str(nt)] = decompose_grid(x[::stride],eta,W[::stride],radius[::stride],
                saved["E_x"],saved["E_Ubar"],source["r"][::stride],Es[::stride],
                fronts=record["geometry"]["support_x"],W0=W[0],budget=budget)
        result.append(dict(amplitude=label,alpha=case["alpha"],source=dict(x=x,
            retained_constraint=source["source_constraint"],jet_constraint=derived,history_constraint=Es,
            jet_minus_retained=derived-source["source_constraint"],history_minus_jet=Es-derived,
            max_jet_identity_difference=float(np.max(abs(derived-source["source_constraint"]))),
            all_saved_cuts_are_original_Hermite_knots=True,
            curvature_convention="right trace at interior knots; left trace only at final endpoint; not a centred curvature"),
            grids=grids,sampling_comparison=compare_sampling(grids["129"],grids["65"],budget)))
    budget.check()
    return result


def audit(root=ROOT, *, budget=None):
    budget = shared_budget(budget)
    initial, phase = budget.evaluations, "input authentication"
    try:
        record, knots, identity = authenticate(root,budget)
        bindings = {name:digest((Path(root)/name).read_bytes()) for name in OWN_SOURCES}
        phase = "saved constraint decomposition"
        cases = analyze(record,knots,budget)
        phase = "final input authentication"
        _,_,after = authenticate(root,budget)
        if identity != after or bindings != {name:digest((Path(root)/name).read_bytes()) for name in OWN_SOURCES}:
            raise ValueError("audit inputs or producer changed during evaluation")
        budget.check()
        return plain(dict(schema=SCHEMA,input=identity,audit_source_bindings=bindings,
            runtime=runtime_versions(),limits=LIMITS,scope=SCOPE,cases=cases,
            execution=dict(evaluated_points=budget.evaluations-initial,cpu_seconds=time.process_time()-budget.started,
                hash_authentication_passes=2,source_evaluations=0,profile_evaluations=0,
                cache_reconstructions=0,stage_calls=0,time_steps=0,scientific_validator_calls=0),
            interpretation=dict(inherited_source_is_exact_solution=False,wall_incoming_constraint="UNKNOWN",
                front_masks="union for matched-grid comparisons",stencil_comparison_is_evolution_refinement=False,
                residual_is_continuum_error=False,uniform_error_enclosure=False)))
    except Exception as error:
        error.constraint_failure = dict(phase=phase,exception=type(error).__name__,message=str(error),
            evaluated_points=budget.evaluations-initial,cpu_seconds=time.process_time()-budget.started)
        raise


def create(path, *, root=ROOT):
    path = Path(path)
    if path.exists():
        raise FileExistsError("constraint audit is creation-only")
    budget = shared_budget()
    try:
        result = audit(root,budget=budget)
        raw = encoded(result)
        budget.check()
    except Exception as error:
        receipt = Path(root)/".build/constraint-transport-failures"/(uuid.uuid4().hex+".json")
        receipt.parent.mkdir(parents=True,exist_ok=True)
        with receipt.open("xb") as stream:
            stream.write(encoded(dict(schema=SCHEMA+"-failed-attempt",status="FAILED_NOT_COMPLETED",limits=LIMITS,
                failure=getattr(error,"constraint_failure",dict(exception=type(error).__name__,message=str(error))))))
        error.failure_receipt = str(receipt)
        raise
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("xb") as stream:
        stream.write(raw)
    return result


def check(path, *, root=ROOT):
    budget = shared_budget()
    saved = json.loads(Path(path).read_bytes())
    current = audit(root,budget=budget)
    for row in (saved,current):
        cpu = row["execution"].pop("cpu_seconds")
        if isinstance(cpu,bool) or not math.isfinite(cpu) or not 0 <= cpu <= 30.:
            raise ValueError("invalid audit CPU receipt")
    if encoded(saved) != encoded(current):
        raise ValueError("constraint transport recomputation mismatch")
    budget.check()
    return dict(checked=str(path),evaluated_points=current["execution"]["evaluated_points"],
        source_evaluations=0,cache_reconstructions=0,stage_calls=0,time_steps=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    try:
        result = check(args.output) if args.check else dict(created=str(args.output),execution=create(args.output)["execution"])
    except Exception as error:
        if hasattr(error,"failure_receipt"):
            print(json.dumps(dict(failed=str(error),diagnostic=error.failure_receipt)))
        raise
    print(json.dumps(result,sort_keys=True,allow_nan=False))


if __name__ == "__main__":
    main()
