#!/usr/bin/env python3
"""Fixed saved-row audit of exterior boundary derivative sensitivity.

No stage, source, profile, gauge quadrature, history reconstruction or time
evolution is called. Localizing a numerical shape does not identify its cause.
"""
from __future__ import annotations

import argparse
import json
import math
from numbers import Real
from pathlib import Path
import time
import uuid

import numpy as np

import record_nsc_atlas_operator_comparison as producer
import analyze_nsc_atlas_spatial_refinement as spatial
from validation_source_successors import verify_saved_bindings, verify_saved_source
from record_nsc_atlas_nonlinear_rectangle import _same_recomputed as same_recomputed
from recursive_horizons.nsc_atlas_departure_gradient import (
    quadratic_departure_gradient, first_order_departure_gradient,
)

ROOT = Path(__file__).resolve().parents[1]
INPUT = producer.RELATIVE_OUTPUT+"/operator.json"
INPUT_SHA256 = "312bc770f7baf63b2aff3d02c26642b32db212febdd34a7e3e6f40fded539d09"
MANIFEST_SHA256 = "9bdf673c83cc0992f9484a3c0f014337e65cdeb0470d953bc0204028676ec1cd"
SCHEMA = "nsc-atlas-boundary-rows-audit-v1"
LIMITS = dict(cpu_seconds=30., evaluated_points=150000)
CUTS = (63, 68, 128)
PATCHES = dict(wall=np.arange(56,65), control=np.arange(28,37))
OWN_SOURCES = ("scripts/diagnose_nsc_atlas_boundary_rows.py",
               "tests/test_diagnose_nsc_atlas_boundary_rows.py")
SCOPE = ("fixed cuts63,68,128 and wall56:64/control28:36; signed saved-array "
         "stencil and method sensitivities, not a constraint certificate or causal diagnosis")
Budget, plain, encoded, digest = producer.Budget, producer.plain, producer.encoded, producer.digest
runtime_versions = producer.runtime_versions
array, same = spatial.array, spatial.same
FIELDS = ("r", "L", "c")


def shared_budget(value=None):
    if value is None:
        value = Budget(maximum_evaluations=150000, cpu_seconds=30.)
    if (not isinstance(value,Budget) or type(value.maximum_evaluations) is not int
            or not 0 < value.maximum_evaluations <= 150000
            or type(value.evaluations) is not int or value.evaluations < 0):
        raise ValueError("require a shared bounded budget")
    for name in ("started","cpu_seconds"):
        v=getattr(value,name)
        if (isinstance(v,(bool,np.bool_)) or not isinstance(v,Real)
                or not math.isfinite(v) or v < 0
                or (name=="cpu_seconds" and not 0 < v <= 30.)):
            raise ValueError("invalid shared budget clock")
    if value.started > time.process_time():
        raise ValueError("budget cannot start in the future")
    value.check()
    return value


def authenticate(root,budget):
    """One hash pass over latest payload/manifest and its104 frozen bindings."""
    root=Path(root).resolve(); path=root/INPUT
    if ({p.name for p in path.parent.iterdir()} != {"operator.json","manifest.json"}
            or any(p.is_symlink() for p in path.parent.iterdir())):
        raise ValueError("operator input inventory mismatch")
    raw=path.read_bytes(); mr=path.with_name("manifest.json").read_bytes()
    manifest=json.loads(mr); bindings=manifest.get("source_bindings",{})
    if (digest(raw)!=INPUT_SHA256 or digest(mr)!=MANIFEST_SHA256
            or manifest.get("schema")!=producer.SCHEMA+"-manifest"
            or manifest.get("outputs")!={"operator.json":dict(sha256=digest(raw),bytes=len(raw))}
            or set(bindings)!=set(producer.SOURCES) or len(bindings)!=104):
        raise ValueError("pinned operator input mismatch")
    retained={}; hashed_bytes=len(raw)+len(mr)
    for name,expected in bindings.items():
        p=(root/name).resolve()
        if not p.is_relative_to(root):
            raise ValueError("input binding exits repository")
        data=verify_saved_source(name,expected,root); hashed_bytes+=len(data)
        if name in (producer.INPUT,producer.CACHE_INPUT):
            retained[name]=data
    if set(retained)!={producer.INPUT,producer.CACHE_INPUT}:
        raise ValueError("required saved baseline/cache bindings missing")
    current,baseline,ancestor=json.loads(raw),json.loads(retained[producer.INPUT]),json.loads(retained[producer.CACHE_INPUT])
    if (current.get("schema")!=producer.SCHEMA or current.get("settings")!=producer.SETTINGS
            or current.get("runtime")!=runtime_versions()):
        raise ValueError("operator declaration/runtime mismatch")
    for key in ("geometry","profile_inputs","cache_identity"):
        if encoded(current[key])!=encoded(baseline[key]):
            raise ValueError("changed common "+key)
    if encoded(current["cache_identity"])!=encoded(producer.cache_identity(ancestor)):
        raise ValueError("ancestor cache identity mismatch")
    budget.check()
    identity=dict(path=INPUT,sha256=digest(raw),manifest_sha256=digest(mr),
        authenticated_input_files=104,source_bindings=bindings,hashed_bytes=hashed_bytes,
        baseline_path=producer.INPUT,baseline_sha256=bindings[producer.INPUT],
        cache_path=producer.CACHE_INPUT,cache_sha256=bindings[producer.CACHE_INPUT],
        recursive_validation=False)
    return current,baseline,ancestor["cache"],identity


def constraint_terms(physical,vacuum,first,second,W,alpha,budget):
    """Five vacuum-subtracted terms using independently supplied eta jets.

    Physical source offsets enter r,c but have no U derivative. No radial or
    temporal constraint is imposed to obtain these derivatives.
    """
    r=np.asarray(physical["r"]); shape=r.shape
    if r.ndim not in (1,2) or r.size>54:
        raise ValueError("bounded selected constraint nodes required")
    p={k:array(physical[k],shape,"physical."+k,budget) for k in ("r","c")}
    v={k:array(vacuum[k],shape,"vacuum."+k,budget)
       for k in ("r_Ubar","r_Ubar_Ubar","log_C_Ubar")}
    g={k:array(first[k],shape,"eta first."+k,budget) for k in FIELDS}
    rr=array(second,shape,"eta radius curvature",budget)
    W=np.asarray(W,dtype=float)
    try: W=np.broadcast_to(W,shape)
    except ValueError as error: raise ValueError("width cannot broadcast to selected nodes") from error
    if (not np.all(np.isfinite(W)) or np.any(W<=0) or np.any(p["r"]<=0)
            or isinstance(alpha,bool) or not isinstance(alpha,Real) or not math.isfinite(alpha) or alpha<0):
        raise ValueError("positive radius/width and finite nonnegative alpha required")
    with np.errstate(over="raise",invalid="raise",divide="raise"):
        dr,gL,dc=g["r"]/W,g["L"]/W,g["c"]/W
        rU=v["r_Ubar"]+dr; LU=v["log_C_Ubar"]+gL
        u=p["c"]/p["r"]; uU=(dc-u*rU)/p["r"]
        terms=dict(departure_curvature=rr/W**2,
            background_lapse=-v["log_C_Ubar"]*dr,
            background_radius=-v["r_Ubar"]*gL,
            departure_product=-dr*gL,
            scalar=alpha*p["r"]*uU*uU/(2*(1+alpha*u*u/6)))
        total=sum(terms.values())
        full=dict(curvature=v["r_Ubar_Ubar"]+terms["departure_curvature"],
                  lapse_product=-LU*rU,scalar=terms["scalar"])
        vac=v["r_Ubar_Ubar"]-v["log_C_Ubar"]*v["r_Ubar"]
        full_total=sum(full.values())
    if not all(np.all(np.isfinite(a)) for a in (*terms.values(),full_total,total)):
        raise ValueError("nonfinite selected constraint")
    budget.check(12*r.size)
    return dict(terms=terms,E_Ubar=total,full_terms=full,full_E_Ubar=full_total,
        analytic_vacuum_residual=vac,algebraic_identity_defect=full_total-total-vac,
        jets=dict(r_Ubar=rU,log_C_Ubar=LU,c_Ubar=dc,u_Ubar=uU))


def gradient_choices(Z,budget):
    """Change operators on IDENTICAL saved arrays, without advancing them."""
    z=array(Z,(3,65),"one cut channel departures",budget); h=1/64
    eta=np.linspace(0,1,65)
    q=quadratic_departure_gradient(eta,z,budget=budget)
    old=first_order_departure_gradient(eta,z,budget=budget)
    third=z[:,64]-3*z[:,63]+3*z[:,62]-z[:,61]
    back63=(3*z[:,63]-4*z[:,62]+z[:,61])/(2*h)
    center62=(z[:,63]-z[:,61])/(2*h)
    differences=np.stack((q[:,62]-center62,q[:,63]-back63),axis=1)
    target=np.stack((-third/(2*h),third/(2*h)),axis=1)
    budget.check(8*z.size+12)
    return dict(quadratic=q,first_order=old,quadratic_minus_first=q-old,
        transition_rows=[62,63],alternative_names=["centered","backward_quadratic"],
        alternatives=np.stack((center62,back63),axis=1),quadratic_minus_alternative=differences,
        third_difference=third,third_difference_target=target,identity_defect=differences-target)


def _selected_cache(cache,budget):
    x=array(cache["x_nodes"],(1025,),"cache x",budget)
    if np.any(np.diff(x)<=0): raise ValueError("cache coordinates must increase")
    ix=np.asarray(CUTS)*8
    keys=("x","W","W_x","W_xx")
    gauge={k:array(np.asarray(cache["gauge_rows"][k])[ix],(3,),"gauge."+k,budget) for k in keys}
    same(gauge["x"],x[ix],"cache gauge keys")
    if np.any(gauge["W"]<=0) or np.any(gauge["W_x"]<=0):
        raise ValueError("positive expanding gauge required")
    sources={}
    for label in ("zero","finite"):
        sources[label]={k:array(np.asarray(cache["source_rows"][label][k])[ix],(3,),"source."+k,budget)
                        for k in ("x","r","log_C","c")}
        same(sources[label]["x"],x[ix],"cache source keys")
    return x[ix],gauge,sources


def _retained(record,method,ci,patch,cut):
    if cut==128: return None
    key="129" if method=="quadratic" else "native64_129"
    data=record["diagnostics"]["constraint_transport"][ci]["grids"][key]
    same(data["eta"],np.linspace(0,1,65)[1:-1],"retained eta")
    return np.asarray(data["arrays"]["E_Ubar"])[cut-1,patch[1:-1]-1]


def analyze(current,baseline,cache,budget):
    """Only three selected cuts, two patches and four fixed saved cases."""
    initial=budget.evaluations
    x,gauge,sources=_selected_cache(cache,budget); W=gauge["W"]
    all_indices=np.r_[PATCHES["wall"],PATCHES["control"]]
    U=W[:,None]*all_indices[None,:]/64
    #54 shared analytic points; no reference-wall construction or quadrature.
    vv=spatial.vacuum_slice(U.ravel(),np.repeat(x,18),budget=budget)
    vacuum={k:np.asarray(v).reshape(3,18) for k,v in vv.items()}
    wall_vacuum={k:v[:,8] for k,v in vacuum.items()}
    cases=[]
    for method,record in (("first_order",baseline),("quadratic",current)):
        if len(record["cases"])!=2: raise ValueError("ordered zero/finite cases required")
        for ci,label in enumerate(("zero","finite")):
            case=record["cases"][ci]; other=baseline["cases"][ci]
            if ((case["amplitude"],case["cells"],case["steps"])!=(label,64,512)
                    or case["alpha"]!=other["alpha"] or len(case["snapshots"])!=129):
                raise ValueError("fixed same-grid cases required")
            if not isinstance(case["alpha"],Real) or case["alpha"]<0 or (label=="zero" and case["alpha"]!=0):
                raise ValueError("invalid fixed amplitude")
            selected=dict(case,snapshots=[case["snapshots"][i] for i in CUTS])
            same([row["x"] for row in selected["snapshots"]],x,"selected saved keys")
            source=sources[label]
            Z,wall=spatial.wall_departures(selected,source,wall_vacuum,budget)
            choices=[gradient_choices(np.stack([Z[k][i] for k in FIELDS]),budget) for i in range(3)]
            patches={}
            for pi,(name,indices) in enumerate(PATCHES.items()):
                v={k:a[:,pi*9:(pi+1)*9] for k,a in vacuum.items()}
                z={k:Z[k][:,indices] for k in FIELDS}
                offsets=dict(r=source["r"]-x,L=source["log_C"],c=source["c"])
                physical={k:(v[k] if k!="c" else np.zeros((3,9)))+offsets[k][:,None]+z[k] for k in FIELDS}
                # Independent saved physical views must reconstruct the same nodes.
                views={}
                for k in FIELDS:
                    vals=[]
                    for row in selected["snapshots"]:
                        full=np.r_[source[k if k!="L" else "log_C"][len(vals)],
                                   np.asarray(row["physical"]["increment"][k])+source[k if k!="L" else "log_C"][len(vals)],
                                   row["physical"]["wall"][k]]
                        vals.append(full[indices])
                    views[k]=array(vals,(3,9),"saved physical patch."+k,budget)
                    if not np.allclose(physical[k],views[k],rtol=3e-13,atol=3e-14):
                        raise ValueError("raw departure/physical patch reconstruction mismatch")
                first={k:(z[k][:,2:]-z[k][:,:-2])*32 for k in FIELDS}
                second=(z["r"][:,2:]-2*z["r"][:,1:-1]+z["r"][:,:-2])*64**2
                terms=constraint_terms({k:a[:,1:-1] for k,a in physical.items()},
                    {k:a[:,1:-1] for k,a in v.items()},first,second,W[:,None],case["alpha"],budget)
                retained=[]; defects=[]
                for i,cut in enumerate(CUTS):
                    old=_retained(record,method,ci,indices,cut)
                    retained.append(None if old is None else old)
                    defect=None if old is None else terms["E_Ubar"][i]-old
                    if defect is not None and not np.allclose(terms["E_Ubar"][i],old,rtol=2e-12,atol=5e-11):
                        raise ValueError("retained incoming constraint does not match selected stencil")
                    defects.append(defect)
                budget.check(12*27+14)
                patches[name]=dict(indices=indices,interior_indices=indices[1:-1],departures=z,
                    physical=physical,physical_reconstruction_defect={k:physical[k]-views[k] for k in FIELDS},
                    centered=terms,retained_E_Ubar=retained,retained_difference=defects,
                    same_data_gradients={key:np.stack([c[key][:,indices] for c in choices])
                                         for key in ("quadratic","first_order","quadratic_minus_first")})
            wall_estimates={}
            for points in (3,5):
                for stride in (1,2):
                    key=f"p{points}_s{stride}"
                    row=spatial.wall_estimate(Z,wall,wall_vacuum,W,case["alpha"],points=points,stride=stride,budget=budget)
                    # The identical wall stencils are already retained in each input.
                    old=record["diagnostics"]["wall_curvature"]["cases"][ci]["estimates"]
                    if method=="first_order": old=old["fine64"]
                    differences={k:row["arrays"][k]-np.asarray(old[key]["arrays"][k])[list(CUTS)] for k in row["arrays"]}
                    if any(np.max(abs(v))>5e-11 for v in differences.values()):
                        raise ValueError("retained wall stencil reconstruction mismatch")
                    row["retained_differences"]=differences; wall_estimates[key]=row
            cases.append(dict(method=method,amplitude=label,alpha=case["alpha"],patches=patches,
                closure_identity=dict(transition_rows=[62,63],alternative_names=["centered","backward_quadratic"],
                    **{k:np.stack([c[k] for c in choices]) for k in
                    ("third_difference","third_difference_target","identity_defect","alternatives","quadratic_minus_alternative")}),
                wall_estimates=wall_estimates))
    comparison=[]
    for ci,label in enumerate(("zero","finite")):
        old,new=cases[ci],cases[ci+2]; patches={}
        for name in PATCHES:
            a,b=new["patches"][name],old["patches"][name]
            terms={k:a["centered"]["terms"][k]-b["centered"]["terms"][k] for k in a["centered"]["terms"]}
            budget.check(12*27)
            patches[name]=dict(departures={k:a["departures"][k]-b["departures"][k] for k in FIELDS},
                physical={k:a["physical"][k]-b["physical"][k] for k in FIELDS},
                centered_terms=terms,E_Ubar=a["centered"]["E_Ubar"]-b["centered"]["E_Ubar"])
        comparison.append(dict(amplitude=label,sign="quadratic evolution minus first-order evolution",patches=patches))
    budget.check()
    return plain(dict(cuts=CUTS,x=x,W=W,eta_spacing=1/64,patch_indices=PATCHES,
        sources=sources,gauge=gauge,cases=cases,same_stencil_method_difference=comparison,
        channel_order=FIELDS,analytic_vacuum_points=54,
        interpretation=dict(same_data_gradient_change_is_evolution_change=False,
            wall_derivative_is_prescribed=False,terminal_cut_has_retained_interior_E=False,
            curvature_is_converged=False,causation_identified=False),
        evaluated_points=budget.evaluations-initial))


def _own_bindings(root):
    return {name:digest((Path(root)/name).read_bytes()) for name in OWN_SOURCES}


def _validation(saved,inputs,identity,bindings,budget,*,root=ROOT):
    if (saved.get("schema")!=SCHEMA or saved.get("scope")!=SCOPE or saved.get("limits")!=LIMITS
            or saved.get("input")!=identity
            or saved.get("runtime")!=runtime_versions()):
        raise ValueError("boundary-row audit declaration mismatch")
    if saved.get("audit_source_bindings") != bindings:
        try:
            verify_saved_bindings(saved.get("audit_source_bindings"),root,expected_paths=bindings)
        except ValueError as error:
            raise ValueError("boundary-row audit declaration mismatch: " + str(error)) from error
    expected=analyze(*inputs,budget)
    if not same_recomputed(saved.get("analysis"),plain(expected)):
        raise ValueError("boundary-row saved arithmetic mismatch")
    e=saved["execution"]; cost=expected["evaluated_points"]
    if (set(e)!={"analysis_evaluations","creation_validation_evaluations","evaluated_points","initial_evaluations",
                 "cpu_seconds","hash_authentication_passes","stage_calls","time_steps","source_integrations",
                 "gauge_quadratures","profile_evaluations","history_reconstructions","scientific_validator_calls"}
            or e["analysis_evaluations"]!=cost or e["creation_validation_evaluations"]!=cost
            or e["evaluated_points"]!=2*cost or type(e["initial_evaluations"])is not int or e["initial_evaluations"]<0
            or e["initial_evaluations"]+2*cost>150000 or e["hash_authentication_passes"]!=1
            or any(e[k]!=0 for k in ("stage_calls","time_steps","source_integrations","gauge_quadratures",
                                       "profile_evaluations","history_reconstructions","scientific_validator_calls"))
            or isinstance(e["cpu_seconds"],bool) or not isinstance(e["cpu_seconds"],Real)
            or not math.isfinite(e["cpu_seconds"]) or not 0<=e["cpu_seconds"]<=30):
        raise ValueError("invalid boundary-row execution receipt")
    return cost


def create(path,*,root=ROOT,budget=None):
    path=Path(path)
    if path.exists(): raise FileExistsError("boundary-row output is creation-only")
    budget=shared_budget(budget); initial=budget.evaluations; phase="input authentication"
    try:
        current,old,cache,identity=authenticate(root,budget); bindings=_own_bindings(root)
        phase="selected saved arithmetic"; result=analyze(current,old,cache,budget)
        cost=result["evaluated_points"]
        payload=plain(dict(schema=SCHEMA,scope=SCOPE,limits=LIMITS,input=identity,audit_source_bindings=bindings,
            runtime=runtime_versions(),analysis=result,execution=dict(analysis_evaluations=cost,
                creation_validation_evaluations=cost,evaluated_points=2*cost,initial_evaluations=initial,
                cpu_seconds=time.process_time()-budget.started,hash_authentication_passes=1,
                stage_calls=0,time_steps=0,source_integrations=0,gauge_quadratures=0,profile_evaluations=0,
                history_reconstructions=0,scientific_validator_calls=0)))
        phase="creation validation"
        _validation(payload,(current,old,cache),identity,bindings,budget,root=root)
        if bindings!=_own_bindings(root): raise ValueError("audit producer changed during calculation")
        payload["execution"]["cpu_seconds"]=time.process_time()-budget.started
        budget.check(); raw=encoded(payload); budget.check()
    except Exception as error:
        receipt=Path(root)/".build/boundary-row-failures"/(uuid.uuid4().hex+".json")
        receipt.parent.mkdir(parents=True,exist_ok=True)
        with receipt.open("xb") as stream:
            stream.write(encoded(dict(schema=SCHEMA+"-failed-attempt",phase=phase,
                reason=str(error),exception=type(error).__name__,evaluated_points=budget.evaluations-initial,
                cpu_seconds=time.process_time()-budget.started,retries=0)))
        error.failure_receipt=str(receipt)
        raise
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("xb") as stream: stream.write(raw)
    return payload


def check(path,*,root=ROOT,budget=None):
    budget=shared_budget(budget); start=budget.evaluations
    saved=json.loads(Path(path).read_bytes())
    current,old,cache,identity=authenticate(root,budget); bindings=_own_bindings(root)
    _validation(saved,(current,old,cache),identity,bindings,budget,root=root)
    if bindings!=_own_bindings(root): raise ValueError("audit producer changed during authentication")
    budget.check()
    return dict(checked=str(path),evaluated_points=budget.evaluations-start,
        cpu_seconds=time.process_time()-budget.started,stage_calls=0,time_steps=0,source_integrations=0)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True); parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    try:
        result=check(args.output) if args.check else dict(created=str(args.output),execution=create(args.output)["execution"])
    except Exception as error:
        if hasattr(error,"failure_receipt"): print(json.dumps(dict(failed=str(error),diagnostic=error.failure_receipt)))
        raise
    print(json.dumps(result,sort_keys=True,allow_nan=False))


if __name__=="__main__": main()
