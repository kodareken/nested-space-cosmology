#!/usr/bin/env python3
"""Estimate and transport the incoming wall constraint from frozen histories."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time
import uuid

import numpy as np

import diagnose_nsc_atlas_constraint_transport as prior
from validation_source_successors import verify_saved_bindings, verify_saved_source
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_vacuum_departure import vacuum_slice

ROOT = prior.ROOT
TRANSPORT_INPUT = "results/development/nsc-atlas-constraint-transport-audit-v1.json"
TRANSPORT_SHA256 = "72e0205dbe8d2854307fc6ef3d5dd761327c75c3451fc5bd11ed6952ce2db9e0"
SCHEMA = "nsc-atlas-wall-constraint-audit-v1"
LIMITS = dict(cpu_seconds=30., evaluated_points=2000000)
OWN_SOURCES = ("scripts/diagnose_nsc_atlas_wall_constraint.py",
               "tests/test_diagnose_nsc_atlas_wall_constraint.py")
VARIANTS = (("backward3_stride1",3,1), ("backward3_stride2",3,2),
            ("backward5_stride1",5,1), ("backward5_stride2",5,2))
encoded, plain, digest = prior.encoded, prior.plain, prior.digest
array, same, close = prior.array, prior.same, prior.close
runtime_versions = prior.runtime_versions


def shared_budget(value=None):
    if value is None:
        value = Budget(maximum_evaluations=LIMITS["evaluated_points"],cpu_seconds=30.)
    if (not isinstance(value,Budget) or type(value.maximum_evaluations) is not int
            or not 0 < value.maximum_evaluations <= LIMITS["evaluated_points"]
            or type(value.evaluations) is not int or value.evaluations < 0
            or isinstance(value.cpu_seconds,bool) or not math.isfinite(value.cpu_seconds)
            or not 0 < value.cpu_seconds <= 30. or isinstance(value.started,bool)
            or not math.isfinite(value.started) or not 0 <= value.started <= time.process_time()):
        raise ValueError("require a live shared allowance within30CPU/2M points")
    value.check()
    return value


def authenticate(root,budget):
    record,_,identity = prior.authenticate(root,budget)
    raw = (Path(root)/TRANSPORT_INPUT).read_bytes()
    if digest(raw) != TRANSPORT_SHA256:
        raise ValueError("pinned source-transport result changed")
    retained = json.loads(raw)
    bindings = retained.get("audit_source_bindings",{})
    if (retained.get("schema") != prior.SCHEMA or retained.get("input") != identity
            or retained.get("runtime") != runtime_versions() or set(bindings) != set(prior.OWN_SOURCES)):
        raise ValueError("source-transport provenance mismatch")
    for name,expected in bindings.items():
        verify_saved_source(name,expected,root)
    budget.check()
    return record,retained,dict(refinement=identity,source_transport=dict(path=TRANSPORT_INPUT,
        sha256=digest(raw),bytes=len(raw),producer_bindings=bindings),scientific_validator_calls=0)


def backward_jets(values,*,points,stride,budget):
    """Differentiate the supplied uniform eta-node values at eta=1."""
    if type(points) is not int or points not in (3,5) or type(stride) is not int or stride not in (1,2):
        raise ValueError("require declared backward stencil and stride")
    nt,ne = np.shape(values)
    if not 2 <= nt <= 129 or not 9 <= ne <= 33:
        raise ValueError("backward jets exceed finite saved grid")
    z = array(values,(nt,ne),"departure values",budget)
    h = stride/(ne-1)
    if points == 3:
        d1,w1 = np.array([3.,-4.,1.]),2*h
        d2,w2 = np.array([2.,-5.,4.,-1.]),h*h
    else:
        d1,w1 = np.array([25.,-48.,36.,-16.,3.]),12*h
        d2,w2 = np.array([35.,-104.,114.,-56.,11.]),12*h*h
    i1 = ne-1-stride*np.arange(len(d1)); i2 = ne-1-stride*np.arange(len(d2))
    if min(i1[-1],i2[-1]) < 0:
        raise ValueError("backward stencil leaves eta interval")
    budget.check(nt*(len(d1)+len(d2)))
    return dict(first=z[:,i1]@d1/w1,second=z[:,i2]@d2/w2,
                first_indices=i1,second_indices=i2,eta_spacing=h)


def wall_estimate(departure,wall,vacuum,W,alpha,*,points,stride,budget):
    nt = len(W)
    W = array(W,(nt,),"wall width",budget)
    if np.any(W <= 0) or np.any(np.diff(W) <= 0) or isinstance(alpha,bool) or not math.isfinite(alpha) or alpha < 0:
        raise ValueError("invalid expanding wall or alpha")
    traces = {k:array(wall[k],(nt,),"wall."+k,budget) for k in ("r","L","c")}
    if np.any(traces["r"] <= 0):
        raise ValueError("wall radius must be positive")
    jets = {k:backward_jets(departure[k],points=points,stride=stride,budget=budget) for k in ("r","L","c")}
    rU = vacuum["r_Ubar"]+jets["r"]["first"]/W
    rUU = vacuum["r_Ubar_Ubar"]+jets["r"]["second"]/(W*W)
    LU = vacuum["log_C_Ubar"]+jets["L"]["first"]/W
    cU = jets["c"]["first"]/W
    r,c = traces["r"],traces["c"]
    with np.errstate(over="raise",invalid="raise",divide="raise"):
        u = c/r; uU = (cU-u*rU)/r
        terms = dict(curvature=rUU,lapse_product=-LU*rU,scalar=alpha*r*uU*uU/(2*(1+alpha*u*u/6)))
        EU = terms["curvature"]+terms["lapse_product"]+terms["scalar"]
        B = r*EU
    if not all(np.all(np.isfinite(v)) for v in (rU,rUU,LU,cU,uU,EU,B)):
        raise ValueError("nonfinite wall constraint arithmetic")
    budget.check(9*nt)
    return dict(points=points,stride=stride,eta_spacing=jets["r"]["eta_spacing"],
        first_indices=jets["r"]["first_indices"],second_indices=jets["r"]["second_indices"],
        jets=dict(r_Ubar=rU,r_Ubar_Ubar=rUU,log_C_Ubar=LU,c_Ubar=cU,u_Ubar=uU),
        terms=terms,E_Ubar=EU,B=B)


def causal_plan(W,U,observations,wall_mask,*,sampling,budget):
    """A birth bracket is usable only when both wall samples already exist."""
    W = array(W,(len(W),),"birth wall widths",budget)
    U = array(U,np.shape(U),"target Ubar",budget)
    obs = np.asarray(observations)
    mask = np.asarray(wall_mask)
    if (sampling not in (1,2) or type(sampling) is not int or len(W)%2 != 1
            or np.any(np.diff(W) <= 0) or W[0] <= 0 or U.ndim != 2
            or obs.dtype.kind not in "iu" or obs.shape != (U.shape[0],)
            or np.any((obs < 0)|(obs >= len(W))) or mask.dtype.kind != "b" or mask.shape != U.shape):
        raise ValueError("invalid causal birth coordinates or mask")
    ix = np.arange(0,len(W),sampling); train = W[ix]
    high = np.searchsorted(train,U,side="left")
    available = high < len(train)
    safe = np.minimum(high,len(train)-1)
    exact = available & (train[safe] == U)
    low = np.where(exact,safe,safe-1)
    available &= low >= 0
    lo = ix[np.maximum(low,0)]; hi = ix[safe]
    available &= hi <= obs[:,None]
    available &= mask
    fraction = np.zeros_like(U)
    distinct = available & (hi != lo)
    fraction[distinct] = (U[distinct]-W[lo[distinct]])/(W[hi[distinct]]-W[lo[distinct]])
    if np.any((fraction[available] < 0)|(fraction[available] > 1)):
        raise ValueError("birth interpolation would extrapolate")
    budget.check(6*U.size)
    return dict(wall_samples=len(ix),sampling=sampling,supported=available,
        unsupported_wall_origin=mask & ~available,
        lower_index=np.where(available,lo,-1),upper_index=np.where(available,hi,-1),fraction=fraction)


def transport_values(B,plan,budget):
    B = array(B,(len(B),),"wall B",budget)
    supported = plan["supported"]
    result = np.zeros(supported.shape)
    lo,hi = plan["lower_index"][supported],plan["upper_index"][supported]
    f = plan["fraction"][supported]
    result[supported] = (1-f)*B[lo]+f*B[hi]
    budget.check(int(np.count_nonzero(supported)))
    return result


def subset_masks(wall_mask,common):
    if wall_mask.shape != common.shape or wall_mask.shape[1] != 31:
        raise ValueError("require original32-cell interior masks")
    # The wall stencil envelope is indices24..32; a centred interior stencil
    # is disjoint iff its highest node index is strictly below24.
    disjoint = np.arange(2,33)[None,:] < 24
    return dict(all_wall_origin=wall_mask & common,
        disjoint_wall_stencils=wall_mask & common & disjoint,
        overlapping_wall_stencils=wall_mask & common & ~disjoint)


def field_summaries(fields,masks,x,eta,U,budget):
    return {k:{name:prior.summary(v,mask,x,eta,U,budget) for name,mask in masks.items()} for k,v in fields.items()}


def masked(values,support):
    result = np.asarray(values,dtype=object)
    result[~support] = None
    return result


def coverage(plan,x,eta,U,budget):
    unavailable = plan["unsupported_wall_origin"]
    budget.check(unavailable.size)
    where = np.argwhere(unavailable)
    first = None
    if len(where):
        i,j = where[0]
        first = dict(interior_index=[int(i),int(j)],time_grid_index=int(i+1),eta_grid_index=int(j+1),
                     x=float(x[i]),eta=float(eta[j]),Ubar=float(U[i,j]))
    return dict(wall_samples=plan["wall_samples"],supported_count=int(np.count_nonzero(plan["supported"])),
        unsupported_count=len(where),first_unsupported=first,
        supported=plan["supported"],unsupported_wall_origin=unavailable,
        lower_index=plan["lower_index"],upper_index=plan["upper_index"],fraction=plan["fraction"])


def mask_domains(masks,x,eta,U,budget):
    result = {}
    for name,mask in masks.items():
        budget.check(mask.size)
        indices = np.argwhere(mask)
        if not len(indices):
            result[name] = dict(count=0,x_interval=None,eta_interval=None,Ubar_interval=None)
        else:
            times,places = indices[:,0],indices[:,1]
            result[name] = dict(count=len(indices),
                x_interval=[float(x[np.min(times)]),float(x[np.max(times)])],
                eta_interval=[float(eta[np.min(places)]),float(eta[np.max(places)])],
                Ubar_interval=[float(np.min(U[mask])),float(np.max(U[mask]))])
    return result


def case_transport(estimates,grid,W,budget):
    nt = len(W)
    if not 5 <= nt <= 129 or nt%2 != 1:
        raise ValueError("require bounded odd wall history")
    shape = (nt-2,31)
    x = array(grid["interior_x"],(nt-2,),"interior x",budget)
    eta = array(grid["interior_eta"],(31,),"interior eta",budget)
    U = array(grid["interior_Ubar"],shape,"interior Ubar",budget)
    radius = array(grid["radius"],shape,"interior radius",budget)
    EU = array(grid["arrays"]["E_Ubar"],shape,"interior E_Ubar",budget)
    wall_mask = np.asarray(grid["masks"]["wall_origin"],dtype=bool)
    same(grid["W"],W,"retained wall gauge")
    same(U,W[1:-1,None]*eta[None,:],"retained interior null labels")
    same(wall_mask,W[1:-1,None]*np.linspace(0,1,33)[None,:-2] > W[0],"wall-origin stencil mask")
    plans = {("129" if stride==1 else "65"):causal_plan(W,U,np.arange(1,nt-1),wall_mask,sampling=stride,budget=budget)
             for stride in (1,2)}
    common = plans["129"]["supported"] & plans["65"]["supported"]
    masks = subset_masks(wall_mask,common)
    weighted = radius*EU
    budget.check(weighted.size)
    transported,records = {},{}
    for name,estimate in estimates.items():
        transported[name],records[name] = {},{}
        for label,plan in plans.items():
            B = transport_values(estimate["B"],plan,budget)
            residual = weighted-B
            budget.check(residual.size)
            transported[name][label] = B
            records[name][label] = dict(B=masked(B,plan["supported"]),residual=masked(residual,plan["supported"]),
                summaries=field_summaries(dict(B=B,residual=residual),masks,x,eta,U,budget))
    sampling = {}
    for name in estimates:
        difference = transported[name]["65"]-transported[name]["129"]
        budget.check(difference.size)
        sampling[name] = dict(B_65_minus_129=masked(difference,common),
            residual_65_minus_129=masked(-difference,common),
            summaries=field_summaries(dict(B=difference,residual=-difference),masks,x,eta,U,budget))
    spatial = {}
    baseline = VARIANTS[0][0]
    for label in plans:
        spatial[label] = {}
        for name in estimates:
            if name == baseline:
                continue
            difference = transported[name][label]-transported[baseline][label]
            budget.check(difference.size)
            spatial[label][name] = dict(B_minus_baseline=masked(difference,common),
                residual_minus_baseline=masked(-difference,common),
                summaries=field_summaries(dict(B=difference,residual=-difference),masks,x,eta,U,budget))
    return dict(raw_r_E_Ubar=weighted,
        raw_summaries=field_summaries(dict(raw=weighted),masks,x,eta,U,budget),
        coverage={k:coverage(v,x,eta,U,budget) for k,v in plans.items()},common_supported=common,
        comparison_masks=masks,comparison_mask_coverage=mask_domains(masks,x,eta,U,budget),
        transported=records,sampling_comparisons=sampling,
        spatial_comparisons=spatial,baseline_variant=baseline,
        unsupported_array_values="null; never interpreted as zero",summary_domain="common causally supported wall-origin subset",
        stencil_independence="interior centred stencil highest index<24; centres1..22; wall stencil envelope24..32",
        coordinate_reference="129 grid in pinned source-transport input; indices refer to that grid")


def analyze(record,retained,budget):
    x = array(record["cache"]["x_nodes"],(1025,),"cache x",budget)[::8]
    gauge = record["cache"]["gauge_rows"]
    W = array(gauge["W"],(1025,),"cached W",budget)[::8]
    if len(record["cases"]) != 2 or len(retained["cases"]) != 2:
        raise ValueError("require both fixed amplitude cases")
    vacuum = vacuum_slice(W,x,budget=budget)  # exactly129 shared analytic points
    for actual,saved in (("r","r"),("L","L"),("r_Ubar","r_Ubar"),("log_C_Ubar","log_C_Ubar")):
        close(vacuum[actual],np.asarray(gauge[saved])[::8],"vacuum/gauge "+saved)
    output = []
    for ci,label in enumerate(("zero","finite")):
        case,old = record["cases"][ci],retained["cases"][ci]
        if (case["amplitude"],case["cells"],case["steps"]) != (label,32,512) or old["amplitude"] != label:
            raise ValueError("fixed case order/settings mismatch")
        same(case["alpha"],old["alpha"],"case alpha")
        snaps = case["snapshots"]
        same([s["x"] for s in snaps],x,"saved wall cuts")
        source = record["cache"]["source_rows"][label]
        same(source["x"],record["cache"]["x_nodes"],"cached source coordinates")
        offsets = {k:array(source["log_C" if k=="L" else k],(1025,),"source."+k,budget)[::8] for k in ("r","L","c")}
        offsets["r"] = offsets["r"]-x
        wall = {k:array([s["physical"]["wall"][k] for s in snaps],(129,),"physical wall."+k,budget) for k in ("r","L","c")}
        departures = {}
        assembly_defects = {}
        for k in ("r","L","c"):
            middle = array([s["state"]["departure"][k] for s in snaps],(129,31),"saved departure."+k,budget)
            base = vacuum[k] if k != "c" else np.zeros(129)
            zW = wall[k]-base-offsets[k]
            departures[k] = np.column_stack((np.zeros(129),middle,zW))
            assembly = base+offsets[k]+zW
            close(assembly,wall[k],"source/vacuum/wall assembly")
            assembly_defects[k] = float(np.max(abs(assembly-wall[k])))
            budget.check(3*129)
        estimates = {name:wall_estimate(departures,wall,vacuum,W,case["alpha"],points=points,stride=stride,budget=budget)
                     for name,points,stride in VARIANTS}
        grid = old["grids"]["129"]
        same(grid["x"],x,"retained129 cuts")
        close(grid["radius"],np.asarray(source["r"])[::8][1:-1,None]+
            np.asarray([s["physical"]["increment"]["r"] for s in snaps])[1:-1],"retained physical radius")
        output.append(dict(amplitude=label,alpha=case["alpha"],wall_estimates=estimates,
            wall_assembly_max_difference=assembly_defects,
            transport=case_transport(estimates,grid,W,budget)))
    return dict(x=x,W=W,cases=output,vacuum_wall_points=129)


def audit(root=ROOT,*,budget=None):
    budget = shared_budget(budget)
    initial,phase = budget.evaluations,"input authentication"
    try:
        record,retained,identity = authenticate(root,budget)
        bindings = {name:digest((Path(root)/name).read_bytes()) for name in OWN_SOURCES}
        runtime = runtime_versions()
        phase = "saved wall curvature and causal transport"
        result = analyze(record,retained,budget)
        phase = "final input authentication"
        _,_,after = authenticate(root,budget)
        if (identity != after or bindings != {name:digest((Path(root)/name).read_bytes()) for name in OWN_SOURCES}
                or runtime != runtime_versions()):
            raise ValueError("audit input or producer changed during evaluation")
        budget.check()
        return plain(dict(schema=SCHEMA,input=identity,audit_source_bindings=bindings,runtime=runtime,limits=LIMITS,
            result=result,execution=dict(evaluated_points=budget.evaluations-initial,cpu_seconds=time.process_time()-budget.started,
                analytic_vacuum_points=129,source_evaluations=0,cache_reconstructions=0,profile_evaluations=0,
                stage_calls=0,time_steps=0,scientific_validator_calls=0,hash_authentication_passes=2),
            interpretation=dict(scope="saved wall-curvature/transport consistency and stencil sensitivity",
                boundary_constraint_is_imposed=False,birth_interpolation="causal piecewise linear; no extrapolation",
                stencil_comparison_is_evolution_refinement=False,uniform_error_enclosure=False,
                disjoint_support_proves_error_independence=False,
                nonzero_wall_constraint_alone_proves_junction_failure=False,
                source_curvature="cached source right traces and centred histories are different diagnostics")))
    except Exception as error:
        error.wall_constraint_failure = dict(phase=phase,exception=type(error).__name__,message=str(error),
            evaluated_points=budget.evaluations-initial,cpu_seconds=time.process_time()-budget.started)
        raise


def create(path,*,root=ROOT):
    path = Path(path)
    if path.exists():
        raise FileExistsError("wall constraint audit is creation-only")
    budget = shared_budget()
    try:
        result = audit(root,budget=budget)
        raw = encoded(result)
        budget.check()
    except Exception as error:
        receipt = Path(root)/".build/wall-constraint-failures"/(uuid.uuid4().hex+".json")
        receipt.parent.mkdir(parents=True,exist_ok=True)
        with receipt.open("xb") as stream:
            stream.write(encoded(dict(schema=SCHEMA+"-failed-attempt",status="FAILED_NOT_COMPLETED",limits=LIMITS,
                failure=getattr(error,"wall_constraint_failure",dict(exception=type(error).__name__,message=str(error))))))
        error.failure_receipt = str(receipt)
        raise
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("xb") as stream:
        stream.write(raw)
    return result


def check(path,*,root=ROOT):
    budget = shared_budget()
    saved = json.loads(Path(path).read_bytes())
    current = audit(root,budget=budget)
    if saved.get("audit_source_bindings") != current.get("audit_source_bindings"):
        verify_saved_bindings(saved.get("audit_source_bindings"),root,expected_paths=OWN_SOURCES)
        current["audit_source_bindings"] = saved["audit_source_bindings"]
    for row in (saved,current):
        cpu = row["execution"].pop("cpu_seconds")
        if isinstance(cpu,bool) or not math.isfinite(cpu) or not 0 <= cpu <= 30.:
            raise ValueError("invalid audit CPU receipt")
    if encoded(saved) != encoded(current):
        raise ValueError("wall constraint recomputation mismatch")
    budget.check()
    return dict(checked=str(path),evaluated_points=current["execution"]["evaluated_points"],
        analytic_vacuum_points=129,source_evaluations=0,stage_calls=0,time_steps=0)


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
