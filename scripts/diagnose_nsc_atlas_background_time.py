#!/usr/bin/env python3
"""Saved-data profile-chain and four-frozen-stage scalar time audit.

No source integration, physical preparation or time evolution is replayed.
The profile quadrature residual is a known chain contribution, not by itself
a bug or an explanation of the observed trajectory error.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq

import record_nsc_atlas_background_evolution as producer
from validation_source_successors import verify_saved_bindings, verify_saved_source

ROOT = Path(__file__).resolve().parents[1]
INPUT = producer.RELATIVE_OUTPUT+"/background.json"
BASE = producer.previous.previous.RELATIVE_OUTPUT+"/coupled.json"
INPUT_SHA256 = "4c98aea3b71f911522f45ca092bca183effec1f5a5c0faf53b79513a746b66ca"
MANIFEST_SHA256 = "b9a657079dbef081d9a343e1814a6f8271582bb9370443dbfd91305401c4584c"
SCHEMA = "nsc-atlas-background-time-audit-v1"
LIMITS = dict(cpu_seconds=30., evaluated_points=150000)
OWN_SOURCES = ("scripts/diagnose_nsc_atlas_background_time.py",
               "tests/test_diagnose_nsc_atlas_background_time.py")
CASES, CUT = (1, 2, 4, 5), 82
plain, encoded, digest = producer.plain, producer.encoded, producer.digest
statistics, Budget = producer.statistics, producer.Budget


def authenticate(root=ROOT):
    """Hash the pinned output, manifest and 66 bound inputs; no validators."""
    root = Path(root).resolve()
    path = root/INPUT
    if ({p.name for p in path.parent.iterdir()} != {"background.json", "manifest.json"}
            or any(p.is_symlink() for p in path.parent.iterdir())):
        raise ValueError("background inventory mismatch")
    raw, mr = path.read_bytes(), path.with_name("manifest.json").read_bytes()
    manifest = json.loads(mr)
    bindings = manifest.get("source_bindings", {})
    if (digest(raw) != INPUT_SHA256 or digest(mr) != MANIFEST_SHA256
            or manifest.get("schema") != producer.SCHEMA+"-manifest"
            or manifest.get("outputs") != {"background.json":dict(sha256=digest(raw),bytes=len(raw))}
            or set(bindings) != set(producer.SOURCES) or len(bindings) != 66):
        raise ValueError("pinned background payload/manifest mismatch")
    base_raw = None
    hashed_bytes = len(raw)+len(mr)
    for name, expected in bindings.items():
        candidate = (root/name).resolve()
        if not candidate.is_relative_to(root):
            raise ValueError("bound input exits the repository")
        data = verify_saved_source(name,expected,root)
        hashed_bytes += len(data)
        if name == BASE:
            base_raw = data
    if base_raw is None:
        raise ValueError("original513-row base is not bound")
    record = json.loads(raw)
    if record.get("schema") != producer.SCHEMA or record.get("settings") != producer.SETTINGS:
        raise ValueError("pinned background declaration mismatch")
    return record, json.loads(base_raw), dict(path=INPUT,sha256=digest(raw),
        manifest_sha256=digest(mr),authenticated_input_files=66,hashed_bytes=hashed_bytes,
        recursive_validation=False)


def _vector(value, name, minimum=2):
    if isinstance(value,(list,tuple)) and any(isinstance(v,(bool,np.bool_)) for v in value):
        raise ValueError(name+" must not contain booleans")
    raw = np.asarray(value)
    if (raw.dtype.kind not in "iuf" or raw.ndim != 1 or not minimum <= raw.size <= 513
            or not np.all(np.isfinite(raw))):
        raise ValueError(name+" needs bounded finite real vectors")
    return np.asarray(raw,dtype=float)


def _crossed(profile, low, high, budget):
    """Strictly interior knots; endpoints do not cross a polynomial piece."""
    knots = profile.breakpoints
    lo, hi = np.minimum(low,high), np.maximum(low,high)
    budget.check(2*np.size(lo))
    return np.searchsorted(knots,hi,side="left") > np.searchsorted(knots,lo,side="right")


def reference_chain(profile, gauge, actual_zero_time, budget):
    """Known reference B-B0 minus trapezoid(B_x), onexisting256/512 grids."""
    values = {k:_vector(gauge[k],k) for k in ("x","T","r","T_x","r_x_total")}
    x,T,R,Tx,Rx = (values[k] for k in ("x","T","r","T_x","r_x_total"))
    if (any(a.shape!=x.shape for a in values.values()) or len(x)!=513
            or np.any(np.diff(x)<=0) or not np.allclose(np.diff(x),np.diff(x)[0],rtol=1e-10,atol=1e-14)):
        raise ValueError("reference audit requires513 uniform rows")
    actual = _vector(actual_zero_time,"actual zero time",129)
    if actual.size!=129:
        raise ValueError("require129 actual saved differences")
    bg = profile.background(T,R,budget=budget)
    B, Bx = bg["chi"], Tx*bg["A0"]*bg["p"]+Rx*bg["q"]
    plus, minus = T+np.arctanh(R), T-np.arctanh(R)
    if np.any(np.diff(plus)<=0) or np.any(np.diff(minus)<0):
        raise ValueError("reference null labels must be monotone")
    budget.check(6*len(x))
    grids = {}
    for n in (256,512):
        stride = 512//n
        xx,b,bx,qp,qm = (a[::stride] for a in (x,B,Bx,plus,minus))
        defect = np.diff(b)-np.diff(xx)*(bx[:-1]+bx[1:])/2
        crossing = (_crossed(profile,qp[:-1],qp[1:],budget)
                    | _crossed(profile,qm[:-1],qm[1:],budget))
        cumulative = lambda a: np.r_[0.,np.cumsum(a)][::n//128]
        total = cumulative(defect)
        crossed, within = cumulative(np.where(crossing,defect,0.)), cumulative(np.where(crossing,0.,defect))
        budget.check(4*n)
        grids[str(n)] = dict(steps=n,step_residual_statistics=statistics(defect),
            R0=total,crossing_cumulative=crossed,within_piece_cumulative=within,
            crossing_steps=np.flatnonzero(crossing),crossing_count=int(np.sum(crossing)),
            within_piece_count=int(np.sum(~crossing)),signed_endpoint=float(total[-1]),
            partition_defect=float(np.max(abs(total-crossed-within))))
    change = grids["512"]["R0"]-grids["256"]["R0"]
    residual = actual-change
    return plain(dict(x=x[::4],B=B[::4],B_x=Bx[::4],grids=grids,
        known_chain_time_difference=change,actual_zero_wall_time_difference=actual,
        actual_minus_known_chain=residual,
        statistics={"known_chain":statistics(change),"actual_zero":statistics(actual),
                    "actual_minus_known_chain":statistics(residual)},
        scope="reference geometry only; incoming field feedback is not included; matching curves do not alone establish causality"))


def linear_profile_path(profile, *, T,R,T_x,R_x,h,budget):
    """Profile-only straight geometry path, not an accelerated SSPRK2 step."""
    data = (T,R,T_x,R_x,h)
    if any(isinstance(v,(bool,np.bool_)) or np.asarray(v).ndim!=0
           or np.asarray(v).dtype.kind not in "iuf" or not math.isfinite(float(v)) for v in data):
        raise ValueError("linear path data must be finite real scalars")
    T,R,Tx,Rx,h = map(float,data)
    if h<=0 or Tx<=0 or min(R,R+h*Rx)<=0 or max(R,R+h*Rx)>=1:
        raise ValueError("linear path exits positive finite radius/time")
    if abs(Rx)>=Tx*(1-max(R,R+h*Rx)**2):
        raise ValueError("linear path must keep both null labels increasing")
    tau = np.array([0.,h/2,h])
    times,radii = T+Tx*tau,R+Rx*tau
    bg = profile.background(times,radii,budget=budget)
    B = bg["chi"]
    rate = Tx*bg["A0"]*bg["p"]+Rx*bg["q"]
    increment = float(B[-1]-B[0])
    coarse = increment-h*(rate[0]+rate[-1])/2
    fine = increment-h*(rate[0]+2*rate[1]+rate[-1])/4
    crossings=[]
    for name,sign in (("incoming",1),("reflected",-1)):
        labels=times+sign*np.arctanh(radii)
        knots=profile.breakpoints
        for knot in knots[(knots>labels[0])&(knots<labels[-1])]:
            def residual(t):
                budget.check(1)
                return T+Tx*t+sign*np.arctanh(R+Rx*t)-knot
            root=brentq(residual,0.,h,xtol=max(h*1e-12,np.finfo(float).tiny),
                         rtol=8*np.finfo(float).eps,maxiter=64)
            crossings.append(dict(branch=name,label=float(knot),time=float(root),fraction=float(root/h)))
    budget.check(9)
    return plain(dict(h=h,T=T,R=R,T_x=Tx,R_x=Rx,times=tau,primitive=B,rate=rate,
        exact_endpoint_increment=increment,one_step_defect=float(coarse),
        two_half_steps_defect=float(fine),fine_minus_coarse=float(fine-coarse),
        crossings=crossings,upper_label_margin=float(profile.domain["last_label"]-max(times+np.arctanh(radii))),
        scope="fixed actual initial geometry rates; straight geometry path only, not actual accelerated trajectory or a propagated error"))


def wall_chain(stage,budget):
    """Independent invariant and canonical wall-rate decompositions."""
    wall,rates,hybrid,metric,cap = (stage[k] for k in
        ("physical_wall","wall_rates","hybrid","cap_metric","cap_rhs"))
    R,f = float(wall["r"]),float(metric["wall"]["f"])
    xs = float(hybrid["Vdot_plus"])
    D = float(hybrid["D"])
    if not all(math.isfinite(v) and v>0 for v in (R,xs,D)):
        raise ValueError("positive wall radius and parameter rate required")
    u = float(wall["c"])/R
    Qc = float(stage["exterior_longitudinal"]["c"][-1])
    returned = float(hybrid["C_minus"])/(2*xs)
    Jc = float(rates["c"])
    Jchi = float(rates["r"])*f+R*float(hybrid["f_s"])/xs
    frame = float(rates["r"])*f+(Jc-u*float(rates["r"]))/math.sqrt(D)
    bg = metric["background_faces"]
    Bx = float(rates["T"])*(1-R*R)*float(bg["p"][-1])+float(rates["r"])*float(bg["q"][-1])
    Zx = float(np.sum(cap["rates"]["Iq"]))
    values=dict(Q_c=Qc,cap_c_return=returned,Jc=Jc,Jchi=Jchi,B_x=Bx,departure_chi_rate=Zx)
    defects=dict(invariant_split_defect=Jc-Qc-returned,
        canonical_frame_defect=Jchi-frame,departure_chain_defect=Jchi-Bx-Zx,
        canonical_flux_defect=Jchi-float(cap["wall_chi_x"]))
    if not all(math.isfinite(v) for v in (*values.values(),*defects.values())):
        raise ValueError("nonfinite wall-chain decomposition")
    budget.check(len(cap["rates"]["Iq"])+12)
    return dict(**values,**defects)


def frozen_stage(case,base,profile,index,budget):
    sample,n = case["snapshots"][index],case["cells"]
    if n!=32 or case["steps"] not in (256,512) or case["amplitude"] not in ("zero","finite"):
        raise ValueError("audit requires the selected32-cell cases")
    gauge=producer.subset(base["gauge"]["rows"],4*index)
    source=producer.subset(base["initial_data"][case["amplitude"]]["snapshot_rows"],4*index)
    if gauge["x"]!=sample["x"] or source["x"]!=sample["x"]:
        raise ValueError("frozen stage is not an exact cached coordinate")
    eta=np.linspace(0.,1.,n+1)
    before=budget.evaluations
    stage=producer.background_stage(sample["state"],x=sample["x"],eta=eta,y_faces=eta,
        source=source,gauge=gauge,alpha=case["alpha"],tube=producer.TUBE,budget=budget,profile=profile)
    producer._equal(stage["physical_wall"],sample["physical"]["wall"],"frozen wall reconstruction")
    producer._equal(stage["physical_cap"],sample["physical"]["cap"],"frozen physical cap")
    row=wall_chain(stage,budget)
    x=np.asarray(base["gauge"]["rows"]["x"])
    h=float(x[2]-x[0])
    rates,wall=stage["wall_rates"],stage["physical_wall"]
    path=linear_profile_path(profile,T=wall["T"],R=wall["r"],T_x=rates["T"],R_x=rates["r"],h=h,budget=budget)
    return plain(dict(amplitude=case["amplitude"],steps=case["steps"],cells=n,alpha=case["alpha"],
        snapshot=index,x=sample["x"],wall=wall,wall_rates=rates,chain=row,linear_path=path,
        evaluated_points=budget.evaluations-before))


def _shared(budget):
    if budget is None:
        return Budget(maximum_evaluations=LIMITS["evaluated_points"],cpu_seconds=LIMITS["cpu_seconds"])
    if (not isinstance(budget,Budget) or type(budget.maximum_evaluations) is not int
            or not 0<budget.maximum_evaluations<=LIMITS["evaluated_points"]
            or isinstance(budget.cpu_seconds,bool) or not math.isfinite(budget.cpu_seconds)
            or not 0<budget.cpu_seconds<=30 or type(budget.evaluations) is not int
            or budget.evaluations<0 or not math.isfinite(budget.started) or budget.started>time.process_time()):
        raise ValueError("invalid shared scalar-time audit allowance")
    return budget


def audit(root=ROOT,*,budget=None):
    budget=_shared(budget)
    budget.check()
    initial=budget.evaluations
    def phase(name,call):
        try:
            result=call()
            budget.check()
            return result
        except Exception as error:
            raise RuntimeError(f"scalar-time audit phase={name}; cumulative_points={budget.evaluations}; "
                               f"cpu_seconds={time.process_time()-budget.started:.9g}; {error}") from error
    record,base,identity=phase("input authentication",lambda: authenticate(root))
    own={name:digest((Path(root)/name).read_bytes()) for name in OWN_SOURCES}
    start=budget.evaluations
    profile,summary=phase("retained profile fit",lambda: producer.fit_profile(record["profile_inputs"],budget))
    producer.profile_replay_equal(summary,record["profile"],record["profile_inputs"])
    profile_count=budget.evaluations-start
    start=budget.evaluations
    reference=phase("reference profile chain",lambda: reference_chain(profile,base["gauge"]["rows"],
        record["diagnostics"]["time_zero"]["wall_arrays"]["c"],budget))
    reference_count=budget.evaluations-start
    stages=[]
    for i in CASES:
        case=record["cases"][i]
        name=f"frozen stage case={i} amplitude={case['amplitude']} steps={case['steps']} cut={CUT}"
        stages.append(phase(name,lambda: frozen_stage(case,base,profile,CUT,budget)))
    keys=("Q_c","cap_c_return","Jc","Jchi","B_x","departure_chi_rate")
    paired={k:(stages[3]["chain"][k]-stages[2]["chain"][k])
               -(stages[1]["chain"][k]-stages[0]["chain"][k]) for k in keys}
    _,_,after=phase("final input authentication",lambda: authenticate(root))
    if after!=identity or own!={name:digest((Path(root)/name).read_bytes()) for name in OWN_SOURCES}:
        raise ValueError("audit inputs or producer changed during evaluation")
    budget.check()
    return plain(dict(schema=SCHEMA,input=identity,audit_source_bindings=own,
        runtime=producer.runtime_versions(),limits=LIMITS,
        settings=dict(cases=CASES,cut=CUT,reference_rows=513,retained_cuts=129,steps=[256,512],
            knot_partition="strict interior null-label crossings; knot endpoints stay within one piece"),
        reference_chain=reference,frozen_stages=stages,paired_time_rate_change=paired,
        execution=dict(evaluated_points=budget.evaluations,initial_evaluations=initial,profile_fit_points=profile_count,
            reference_points=reference_count,frozen_stage_points=sum(v["evaluated_points"] for v in stages),
            cpu_seconds=time.process_time()-budget.started,hash_authentication_passes=2,
            profile_stage_extractions=0,source_integrations=0,time_evolutions=0,frozen_stage_reconstructions=4),
        scope="known profile chain contribution and frozen-stage rates; neither establishes a causal error attribution or a uniform bound"))


def create(path,*,root=ROOT):
    path=Path(path)
    if path.exists():
        raise FileExistsError("audit output is creation-only")
    result=audit(root)
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("xb") as stream:
        stream.write(encoded(result))
    return result


def check(path,*,root=ROOT):
    saved=json.loads(Path(path).read_bytes())
    current=audit(root)
    # Provenance is authenticated separately; scientific fields remain exact.
    if saved.get("audit_source_bindings") != current.get("audit_source_bindings"):
        verify_saved_bindings(saved.get("audit_source_bindings"),root,expected_paths=OWN_SOURCES)
        current["audit_source_bindings"] = saved["audit_source_bindings"]
    for row in (saved,current):
        value=row["execution"].pop("cpu_seconds")
        if isinstance(value,bool) or not math.isfinite(value) or not 0<=value<=30:
            raise ValueError("invalid saved audit CPU")
    if encoded(saved)!=encoded(current):
        raise ValueError("scalar-time audit recomputation mismatch")
    return dict(checked=str(path),source_integrations=0,time_evolutions=0,
                frozen_stage_reconstructions=4,evaluated_points=current["execution"]["evaluated_points"])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    if args.check:
        result=check(args.output)
    else:
        row=create(args.output)
        result=dict(created=str(args.output),execution=row["execution"])
    print(json.dumps(result,sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()
