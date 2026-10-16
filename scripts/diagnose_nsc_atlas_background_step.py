#!/usr/bin/env python3
"""One-interval SSPRK2 scalar ledger from pinned saved starting states.

Two amplitudes each take one coarse step and two half steps from the same
start. Exactly twelve coupled stages are evaluated. No source integration,
profile-stage extraction or full trajectory is replayed. The ledger separates
flux quadrature, nonlinear frame, reference geometry and background increments;
it does not assign the whole retained temporal error to any single term.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

import diagnose_nsc_atlas_background_time as previous
from recursive_horizons.nsc_atlas_background_evolution import background_stage, _combine, _state
from recursive_horizons.nsc_atlas_reduced_cap import _wall_conversion

ROOT=Path(__file__).resolve().parents[1]
SCHEMA="nsc-atlas-background-step-audit-v1"
LIMITS=dict(cpu_seconds=30.,evaluated_points=250000,stage_calls=12)
SOURCES=(*previous.OWN_SOURCES,"scripts/diagnose_nsc_atlas_background_step.py",
         "tests/test_diagnose_nsc_atlas_background_step.py")
CASES,CUT=(1,4),82
COMPONENTS=("exterior_flux","cap_return_flux","frame_remainder",
            "reference_geometry_remainder","background_remainder")
producer=previous.producer
Budget,plain,encoded,digest=previous.Budget,previous.plain,previous.encoded,previous.digest


def frame(radius,chi,alpha):
    """Canonical primitive to invariant scalar, with its smooth alpha=0 limit."""
    values=(radius,chi,alpha)
    if any(isinstance(v,(bool,np.bool_)) or np.asarray(v).ndim!=0
           or np.asarray(v).dtype.kind not in "iuf" or not math.isfinite(float(v)) for v in values):
        raise ValueError("frame needs finite real scalars")
    radius,chi,alpha=map(float,values)
    if radius<=0 or alpha<0:
        raise ValueError("frame requires positive radius and nonnegative alpha")
    answer=radius*_wall_conversion(chi/radius,alpha)
    if not math.isfinite(answer):
        raise ValueError("nonfinite invariant scalar")
    return answer


def compact_state(state,gauge,profile,alpha,budget):
    """Accepted/predictor readout without an extra metric or coupled stage."""
    wall={k:float(gauge[k])+float(state["wall"][k]) for k in ("r","L","T","s")}
    bg=profile.background(wall["T"],wall["r"],budget=budget)
    B=float(bg["chi"])
    # Match the frozen cap's ordered primitive exactly. Pairwise np.sum can
    # differ by roundoff and would fail the bitwise saved-wall comparison.
    Z=float(np.cumsum(state["cap"]["Iq"])[-1])
    chi=B+Z
    wall["c"]=frame(wall["r"],chi,alpha)
    budget.check(len(state["cap"]["Iq"])+8)
    return dict(wall=wall,background_chi=B,departure_chi=Z,canonical_chi=chi,
                canonical_f=chi/wall["r"],raw_state_sha256=digest(encoded(plain(state))))


def _cache(rows,index):
    if type(index) is not int or index<0 or index>=len(rows["x"]) or len(rows["x"])>513:
        raise ValueError("cached row index outside bounded inputs")
    row=producer.subset(rows,index)
    if not math.isfinite(float(row["x"])):
        raise ValueError("nonfinite cached coordinate")
    return row


def accepted_step(initial,left,right,*,source_rows,gauge_rows,profile,alpha,budget,counter,label):
    """Advance exactly one frozen SSPRK2 step and record its scalar identity.

    Reference geometry is restored at each cached coordinate, including the
    predictor and accepted endpoint. The accepted endpoint needs only a scalar
    profile readout; there is no thirteenth/final coupled-stage evaluation.
    """
    start_count=budget.evaluations
    g0,g1=_cache(gauge_rows,left),_cache(gauge_rows,right)
    s0,s1=_cache(source_rows,left),_cache(source_rows,right)
    x0,x1=float(g0["x"]),float(g1["x"])
    h=x1-x0
    if right<=left or h<=0 or s0["x"]!=x0 or s1["x"]!=x1:
        raise ValueError("step needs increasing exactly matched cached source/gauge rows")
    ne,nc=len(initial["departure"]["c"])+1,len(initial["cap"]["Iq"])
    state=_state(initial,ne,nc)
    eta,y=np.linspace(0.,1.,ne+1),np.linspace(0.,1.,nc+1)
    phase="initial readout"
    def evaluate(value,gauge,source,name):
        if counter["stage_calls"]>=LIMITS["stage_calls"]:
            raise ValueError("twelve-stage ceiling exhausted")
        counter["stage_calls"]+=1
        row=background_stage(value,x=float(gauge["x"]),eta=eta,y_faces=y,
            source=source,gauge=gauge,alpha=alpha,tube=producer.TUBE,budget=budget,profile=profile)
        courant=h*float(row["courant_rate"])
        if not math.isfinite(courant) or not 0<=courant<=.4:
            raise ValueError("declared characteristic Courant limit exceeded")
        counter["max_courant"]=max(counter["max_courant"],courant)
        return row,courant
    try:
        initial_view=compact_state(state,g0,profile,alpha,budget)
        phase="initial stage"
        first,cfl0=evaluate(state,g0,s0,phase)
        producer._equal(first["physical_wall"],initial_view["wall"],"initial independent scalar readout")
        chain0=previous.wall_chain(first,budget)
        predictor=_combine(state,first["rates"],1.,h)
        phase="predictor readout"
        predictor_view=compact_state(predictor,g1,profile,alpha,budget)
        phase="predictor stage"
        second,cfl1=evaluate(predictor,g1,s1,phase)
        producer._equal(second["physical_wall"],predictor_view["wall"],"predictor independent scalar readout")
        chainp=previous.wall_chain(second,budget)
        corrected=_combine(predictor,second["rates"],1.,h)
        accepted=_combine(state,corrected,.5,.5)
        phase="accepted readout and ledger"
        final=compact_state(accepted,g1,profile,alpha,budget)
        trap=lambda key:h*(chain0[key]+chainp[key])/2
        w0,wp,w1=(v["wall"] for v in (initial_view,predictor_view,final))
        Rrk=w0["r"]+h*(first["wall_rates"]["r"]+second["wall_rates"]["r"])/2
        Trk=w0["T"]+h*(first["wall_rates"]["T"]+second["wall_rates"]["T"])/2
        chi_rk=initial_view["canonical_chi"]+trap("Jchi")
        RB=final["background_chi"]-initial_view["background_chi"]-trap("B_x")
        c_rk=frame(Rrk,chi_rk,alpha)
        c_at_actual_radius=frame(w1["r"],chi_rk,alpha)
        c_with_background=frame(w1["r"],chi_rk+RB,alpha)
        parts=dict(exterior_flux=trap("Q_c"),cap_return_flux=trap("cap_c_return"),
            frame_remainder=c_rk-w0["c"]-trap("Jc"),
            reference_geometry_remainder=c_at_actual_radius-c_rk,
            background_remainder=c_with_background-c_at_actual_radius)
        increment=w1["c"]-w0["c"]
        geometry={}
        for key,rate,trial in (("r","r_x_total",Rrk),("T","T_x",Trk)):
            predicted=float(g1[key])-float(g0[key])-h*(float(g0[rate])+float(g1[rate]))/2
            geometry[key]=dict(physical_rk=trial,accepted=w1[key],actual_defect=w1[key]-trial,
                               reference_defect=predicted,closure=w1[key]-trial-predicted)
        predictor_geometry={k:wp[k]-w0[k]-h*first["wall_rates"][k] for k in ("r","T")}
        closures=dict(canonical=final["canonical_chi"]-chi_rk-RB,
            invariant=increment-sum(parts.values()),
            invariant_flux=trap("Jc")-trap("Q_c")-trap("cap_c_return"),
            departure=final["departure_chi"]-initial_view["departure_chi"]-trap("departure_chi_rate"))
        budget.check(40)
        if not all(math.isfinite(v) for v in (*parts.values(),*closures.values())):
            raise ValueError("nonfinite accepted-step ledger")
        row=dict(label=label,left_index=left,right_index=right,x_left=x0,x_right=x1,h=h,alpha=alpha,
            initial=initial_view,predictor=predictor_view,accepted=final,
            stages=[dict(chain=chain0,wall_rates=first["wall_rates"],courant=cfl0),
                    dict(chain=chainp,wall_rates=second["wall_rates"],courant=cfl1)],
            rk=dict(radius=Rrk,T=Trk,chi=chi_rk,c=c_rk),background_chain=RB,
            reference_geometry=geometry,predictor_geometry_defect=predictor_geometry,
            increment=increment,components=parts,closures=closures,
            evaluated_points=budget.evaluations-start_count)
        return accepted,plain(row)
    except Exception as error:
        raise RuntimeError(f"local scalar ledger {label}; phase={phase}; indices={left},{right}; "
            f"stage_calls={counter['stage_calls']}; cumulative_points={budget.evaluations}; "
            f"cpu_seconds={time.process_time()-budget.started:.9g}; {error}") from error


def _shared(budget):
    if budget is None:
        return Budget(maximum_evaluations=250000,cpu_seconds=30.)
    if (not isinstance(budget,Budget) or type(budget.maximum_evaluations) is not int
            or not 0<budget.maximum_evaluations<=250000 or type(budget.evaluations) is not int
            or budget.evaluations<0 or isinstance(budget.cpu_seconds,bool)
            or not math.isfinite(budget.cpu_seconds) or not 0<budget.cpu_seconds<=30
            or not math.isfinite(budget.started) or budget.started>time.process_time()):
        raise ValueError("invalid shared local-step allowance")
    budget.check()
    return budget


def audit(root=ROOT,*,budget=None):
    budget=_shared(budget)
    initial_count=budget.evaluations
    counter=dict(stage_calls=0,max_courant=0.)
    phase="input authentication"
    try:
        record,base,identity=previous.authenticate(root)
        bindings={name:digest((Path(root)/name).read_bytes()) for name in SOURCES}
        phase="retained profile fit"
        before=budget.evaluations
        profile,summary=producer.fit_profile(record["profile_inputs"],budget)
        producer._equal(summary,record["profile"],"retained computational profile")
        profile_count=budget.evaluations-before
        gauge=base["gauge"]["rows"]
        x=np.asarray(gauge["x"],dtype=float)
        if (x.shape!=(513,) or np.any(np.diff(x)<=0) or not np.all(np.isfinite(x))
                or not np.allclose(np.diff(x),np.diff(x)[0],rtol=1e-10,atol=1e-14)):
            raise ValueError("require the original uniform513 cached times")
        left,mid,right=4*CUT,4*CUT+1,4*CUT+2
        cases=[]
        for case_id in CASES:
            case=record["cases"][case_id]
            label=case["amplitude"]
            phase=f"case={case_id} amplitude={label}"
            if (case["cells"],case["steps"],label)!=(32,256,"zero" if case_id==1 else "finite"):
                raise ValueError("require saved zero/finite32/256 starts")
            sample=case["snapshots"][CUT]
            if sample["x"]!=x[left]:
                raise ValueError("starting state is not the exact cached coordinate")
            source=base["initial_data"][label]["snapshot_rows"]
            options=dict(source_rows=source,gauge_rows=gauge,profile=profile,alpha=case["alpha"],
                         budget=budget,counter=counter)
            initial=sample["state"]
            _,coarse=accepted_step(initial,left,right,label=label+" coarse",**options)
            producer._equal(coarse["initial"]["wall"],sample["physical"]["wall"],"saved starting wall")
            middle,half1=accepted_step(initial,left,mid,label=label+" first half",**options)
            _,half2=accepted_step(middle,mid,right,label=label+" second half",**options)
            fine={k:half1["components"][k]+half2["components"][k] for k in COMPONENTS}
            difference={k:fine[k]-coarse["components"][k] for k in COMPONENTS}
            delta=half2["accepted"]["wall"]["c"]-coarse["accepted"]["wall"]["c"]
            cases.append(dict(case=case_id,amplitude=label,alpha=case["alpha"],coarse=coarse,
                halves=[half1,half2],fine_components=fine,fine_minus_coarse=difference,
                fine_minus_coarse_increment=delta,comparison_closure=delta-sum(difference.values())))
        if counter["stage_calls"]!=12:
            raise ValueError("local ledger must contain exactly twelve stage calls")
        paired={k:cases[1]["fine_minus_coarse"][k]-cases[0]["fine_minus_coarse"][k] for k in COMPONENTS}
        paired_value=cases[1]["fine_minus_coarse_increment"]-cases[0]["fine_minus_coarse_increment"]
        phase="final input authentication"
        _,_,after=previous.authenticate(root)
        if identity!=after or bindings!={name:digest((Path(root)/name).read_bytes()) for name in SOURCES}:
            raise ValueError("ledger input or producer changed during evaluation")
        budget.check()
        return plain(dict(schema=SCHEMA,input=identity,audit_source_bindings=bindings,limits=LIMITS,
            runtime=producer.runtime_versions(),settings=dict(cases=CASES,cut=CUT,cache_indices=[left,mid,right],
                time_scheme="frozen SSPRK2 _combine; one coarse and two halves from identical raw start",
                exterior_quadrature_order=2,cap_quadrature_order=4,maximum_courant=.4),
            cases=cases,paired_fine_minus_coarse=paired,paired_increment=paired_value,
            paired_closure=paired_value-sum(paired.values()),
            execution=dict(initial_evaluations=initial_count,evaluated_points=budget.evaluations,
                profile_fit_points=profile_count,local_step_points=sum(c["coarse"]["evaluated_points"]
                    +sum(h["evaluated_points"] for h in c["halves"]) for c in cases),
                cpu_seconds=time.process_time()-budget.started,stage_calls=12,local_ssprk2_steps=6,
                full_trajectory_replays=0,source_integrations=0,profile_stage_extractions=0,
                accepted_endpoint_extra_stage_calls=0,hash_authentication_passes=2,
                maximum_courant=counter["max_courant"],retries=0),
            scope="one actual interval and exact increment decomposition; local coarse/fine sensitivity, not a global error bound or causal attribution for the full saved trajectory"))
    except Exception as error:
        raise RuntimeError(f"local scalar audit phase={phase}; stage_calls={counter['stage_calls']}; "
            f"cumulative_points={budget.evaluations}; cpu_seconds={time.process_time()-budget.started:.9g}; {error}") from error


def create(path,*,root=ROOT):
    path=Path(path)
    if path.exists():
        raise FileExistsError("local scalar audit output is creation-only")
    result=audit(root)
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("xb") as stream:
        stream.write(encoded(result))
    return result


def check(path,*,root=ROOT):
    saved=json.loads(Path(path).read_bytes())
    current=audit(root)
    for row in (saved,current):
        cpu=row["execution"].pop("cpu_seconds")
        if isinstance(cpu,bool) or not math.isfinite(cpu) or not 0<=cpu<=30:
            raise ValueError("invalid ledger CPU")
    if encoded(saved)!=encoded(current):
        raise ValueError("local scalar ledger recomputation mismatch")
    return dict(checked=str(path),stage_calls=12,local_ssprk2_steps=6,full_trajectory_replays=0,
                source_integrations=0,evaluated_points=current["execution"]["evaluated_points"])


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
