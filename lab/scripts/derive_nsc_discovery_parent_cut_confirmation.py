#!/usr/bin/env python3
"""One NF256 confirmation at the frozen NF128 absolute proper-clock target."""
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import argparse
import json
import os
for name in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[name]="1"
import numpy as np
from threadpoolctl import threadpool_limits
from recursive_horizons import nsc_discovery_parent_cut_response as cut
from recursive_horizons import nsc_discovery_parent_episode as episode
from recursive_horizons import nsc_discovery_backend as backend

SCHEMA="NSC-DISCOVERY-PARENT-CUT-CONFIRMATION-v1"
PARENT_COMMIT="7004a7915a9375c2d5a420853158c4cf19bcdca1"
TARGET=2.8568786296681035
CAP=.0005
BUDGET=600.
OUTPUT=cut.LAB/"results/development/nsc-discovery-parent-cut-confirmation-v1"
OWN=("lab/scripts/derive_nsc_discovery_parent_cut_confirmation.py",
     "lab/tests/test_nsc_discovery_parent_cut_confirmation.py","lab/docs/nsc-discovery-parent-cut-confirmation.md")
NF128_EFFECT={"original_column0_child":.0004292078493,"child_content":.0028420021143,
              "child_proper_length":.0047672394727}


def authenticate(source=cut.OUTPUT):
    source=Path(source).resolve();records={};inputs={}
    for stage in ("prepare","prediction","measurement-a0.05-T2.25"):
        record,a=cut.response._load_stage(source,stage)
        if record["producing_commit"]!=PARENT_COMMIT:raise ValueError("requires frozen7004 cut producer")
        for path,digest in record["producers"].items():
            cut.provenance.resolve_pinned_source_bytes(cut.ROOT,path,digest,commit=PARENT_COMMIT)
            if path.startswith("lab/src/") and cut.parent._sha_file(cut.ROOT/path)!=digest:
                raise ValueError("executing frozen owner changed: "+path)
        for path,digest in record["inputs"].items():
            if cut.parent._sha_file(path)!=digest:raise ValueError("frozen input changed")
        inputs.update(record["inputs"])
        inputs.update({str(source/(stage+suffix)):cut.parent._sha_file(source/(stage+suffix)) for suffix in (".json",".npz")})
        records[stage]=(record,a)
    prepared,a=records["prepare"];forecast,_=records["prediction"]
    measured,_=records["measurement-a0.05-T2.25"]
    target=next(row for row in forecast["forecast"]["rows"] if abs(row["time"]-2.25)<1e-11)
    if not forecast["forecast_locked"] or target["tau"]!=TARGET or forecast["prepare_json_sha256"]!=inputs[str(source/"prepare.json")]:
        raise ValueError("absolute proper-clock lock changed")
    if not measured["matched_proper_clock"] or abs(measured["measurement"]["tau"]-TARGET)>1e-9:
        raise ValueError("completed NF128 reference does not match the locked clock")
    return prepared,a,inputs


def periodic_extension(values,new_size,length):
    return cut.extent.real_periodic_values(SimpleNamespace(nq=len(values),length=length),values,
                                          np.arange(new_size)*length/new_size)


def extend(pair,state,nf=256):
    """Physical periodic interpolation and canonical AP isometry; no normalization."""
    grid=cut.parent.galerkin.build_grid(nf,gauge="conformal",length=pair.grid.length)
    carrier=backend.SpinorCarrier(pair.grid.nf,nf,pair.grid.length,canonical=True)
    columns=lambda frame:np.vstack((carrier@frame[:pair.grid.nf],carrier@frame[pair.grid.nf:]))
    source,reference=columns(pair.source_columns),columns(pair.reference_columns)
    holder=cut.response.make_holder(grid,weights=pair.weights,source_columns=source,reference_columns=reference,
        geometry_map=episode.full_frame(grid.ng),source_metadata=dict(pair.source_metadata,band_extension_only=True))
    nodal=cut.leading.decode(pair,state)
    lifted=cut.leading.State(*(periodic_extension(getattr(nodal,n),grid.ng,grid.length) for n in cut.leading.REAL_FIELDS),
                             carrier@state.phi0,carrier@state.phi1)
    result=cut.leading.encode(holder,lifted)
    gram=lambda s:s.phi0.conj().T@s.phi0+s.phi1.conj().T@s.phi1
    report={"from_nf":pair.grid.nf,"nf":nf,"new_W":"existing full identity frame",
        "source_Gram_gap":float(np.max(abs(source.conj().T@source-pair.source_columns.conj().T@pair.source_columns))),
        "reference_Gram_gap":float(np.max(abs(reference.conj().T@reference-pair.reference_columns.conj().T@pair.reference_columns))),
        "live_Gram_gap":float(np.max(abs(gram(result)-gram(state)))),"source_reselected":False,
        "source_normalized":False,"geometry_relaxed":False,"canonical_AP_density_units_preserved":True}
    if max(report[n] for n in ("source_Gram_gap","reference_Gram_gap","live_Gram_gap"))>1e-10:
        raise ValueError("band extension changed actual frame Gram")
    cut.leading.check_chart(holder,result);cut.car(holder,holder.weights,result)
    return holder,result,report


def reprepare(pair,state,k,weights,remaining):
    held,solved,report=cut.reprepare(pair,state,k,-1,weights,cpu_limit=min(30.,remaining))
    for name in ("Q","r","phi0","phi1"):
        if not np.array_equal(getattr(solved,name),getattr(state,name)):raise ValueError("hidden initial geometry/source change")
    fine=cut.leading.check_chart(held,solved);constraints=cut.leading.constraints(held,solved)
    if np.max(fine.p_Q)>=0 or max(constraints["raw_C_max"],constraints["D_max"])>1e-6:
        raise ValueError("finite initial constraints/sign unresolved")
    report=dict(report,actual_constraints=constraints,Q_r_Phi_fixed=True,k_fixed=k,
        initial_momentum_correction={n:float(np.max(abs(getattr(solved,n)-getattr(state,n)))) for n in ("p_Q","p_r")},
        momentum_correction_representation="canonical pi=dx_g W.T p")
    return held,solved,report


def evolve(pair,state,target,deadline,cap=CAP):
    time=tau=0.;steps=0;stop="MATCHED_PROPER_CLOCK"
    try:
        while tau<target-1e-11:
            if cut.parent._cpu_time()>=deadline-1.:stop="CPU_BUDGET_STOP";break
            dt,_=cut.response._admitted_step(pair,state,cap)
            increment=cut.response._nonlinear_centre_increment(pair,state,dt)
            if not np.isfinite(increment) or increment<=0:raise ValueError("nonpositive proper-clock increment")
            if tau+increment>target:
                lo=0.;hi=dt
                for _ in range(32):
                    if cut.parent._cpu_time()>=deadline-1.:stop="CPU_BUDGET_STOP";break
                    mid=(lo+hi)/2
                    if tau+cut.response._nonlinear_centre_increment(pair,state,mid)<target:lo=mid
                    else:hi=mid
                if stop=="CPU_BUDGET_STOP":break
                dt=(lo+hi)/2;increment=cut.response._nonlinear_centre_increment(pair,state,dt)
            advanced=cut.leading.rk4_step(pair,state,dt)
            state=advanced;time+=dt;tau+=increment;steps+=1
    except cut.leading.coupling.PositiveChartExit:stop="POSITIVE_CHART_STOP"
    return state,{"time":time,"tau":tau,"steps":steps,"stop":stop,"target_tau":target,
        "matched":stop=="MATCHED_PROPER_CLOCK" and abs(tau-target)<1e-9,"no_extrapolation":True}


def observe(pair,state):
    values,_=cut.observables(pair,state);fine=cut.leading.fine_state(pair,state)
    mass=(abs(fine.phi0[:,0])**2+abs(fine.phi1[:,0])**2)*pair.weights[0]/pair.grid.dx_q
    values["original_column0_child"]=cut.extent.real_interval_integral(pair.grid,mass,pair.child_interval)
    rate,bundle=cut.leading.rates(pair,state,return_bundle=True)
    velocity=cut.leading.prolong(pair.grid,cut.leading.decode(pair,rate))
    return {"values":values,"constraints":cut.leading.constraints(pair,state),"CAR":cut.car(pair,pair.weights,state),
        "metric_rate_projection_defect":{n:float(np.max(abs(getattr(velocity,n)-bundle["unprojected_rates"][i])))
                                         for i,n in enumerate(cut.leading.REAL_FIELDS)},
        "metric_Euler_continuum_defect":None,"rate_defect_is_continuum_bound":False,
        "tagged_columns_are_particles_or_energy":False}


def preview(source=cut.OUTPUT):
    prepared,_a,inputs=authenticate(source)
    return {"schema":SCHEMA,"mode":"read_only_preview","nf":256,"step_cap":CAP,"target_tau":TARGET,
        "amplitudes":[0.,.05],"k":prepared["baseline_record"]["k"],"c0_fixed":True,
        "CPU_limit_seconds":BUDGET,"input_hashes":inputs,"evolved":False,"bytes_written":0}


def run(source=cut.OUTPUT,output=OUTPUT,*,producer_commit=None,cpu_budget=BUDGET):
    if not producer_commit or not 0<cpu_budget<=BUDGET:raise ValueError("requires frozen producer and total CPU allowance<=600")
    output=Path(output).resolve()
    if output.exists():raise FileExistsError("confirmation output must be creation-only")
    start=cut.parent._cpu_time();deadline=start+cpu_budget
    prepared,a,inputs=authenticate(source)
    pins=dict(cut.pins(),**{n:cut.parent._sha_file(cut.ROOT/n) for n in OWN})
    for path,digest in pins.items():cut.provenance.resolve_pinned_source_bytes(cut.ROOT,path,digest,commit=producer_commit)
    base,initial=cut.reconstruct(prepared["baseline_record"],dict(a,pi_Q=a["p_Q"],pi_r=a["p_r"]))
    record={"schema":SCHEMA,"producing_commit":producer_commit,"producers":pins,"input_hashes":inputs,
        "frozen_parent_commit":PARENT_COMMIT,"numerical_binding":{"nf":256,"step_cap":CAP,"target_tau":TARGET},
        "CPU_limit_seconds":cpu_budget,"NF128_effect":NF128_EFFECT,"arms":[],"comparison":None,
        "source_closure":cut.response.source_closure(),"frozen_producer":True,"physical_binding":True,
        "finite_band_indicator_not_continuum_error_bound":True}
    arrays={}
    with threadpool_limits(limits=1),backend.fft_thread_limit(1):
        pair,extended,extension=extend(base,initial);record["band_extension"]=extension
        record["initial_band_constraints"]={"NF128":cut.leading.constraints(base,initial),
                                            "NF256_before_momentum_correction":cut.leading.constraints(pair,extended)}
        arrays.update({"extended_initial_"+n:v for n,v in cut._arrays(pair,extended).items()})
        for index,amplitude in enumerate((0.,.05)):
            remaining=deadline-cut.parent._cpu_time()
            if remaining<2.:record["stop"]="CPU_BUDGET_STOP_BEFORE_ARM";break
            weights=pair.weights.copy();weights[1]*=1+amplitude
            try:held,state,prep=reprepare(pair,extended,prepared["baseline_record"]["k"],weights,remaining-1.)
            except ValueError as error:record["arms"].append({"amplitude":amplitude,"stop":"OPEN_INITIAL_PREPARATION","reason":str(error)});break
            if held.weights[0]!=pair.weights[0]:raise ValueError("original column0 occupation changed")
            state,event=evolve(held,state,TARGET,deadline)
            record["arms"].append(dict(event,amplitude=amplitude,preparation=prep,endpoint=observe(held,state),weights=weights.tolist()))
            arrays.update({f"arm{index}_"+n:v for n,v in cut._arrays(held,state).items()})
            if not event["matched"]:break
    if len(record["arms"])==2 and all(arm.get("matched") for arm in record["arms"]):
        record["comparison"]={n:{"NF256_effect":record["arms"][1]["endpoint"]["values"][n]-record["arms"][0]["endpoint"]["values"][n],
            "NF128_effect":value,"indicator_difference":record["arms"][1]["endpoint"]["values"][n]-record["arms"][0]["endpoint"]["values"][n]-value}
            for n,value in NF128_EFFECT.items()}
    record["CPU_seconds"]=cut.parent._cpu_time()-start;record["remaining_CPU_seconds"]=max(0.,deadline-cut.parent._cpu_time())
    if record["CPU_seconds"]>cpu_budget:record["stop"]="AGGREGATE_CPU_LIMIT"
    for path,digest in inputs.items():
        if cut.parent._sha_file(path)!=digest:raise ValueError("input changed during confirmation")
    if any(cut.parent._sha_file(cut.ROOT/path)!=digest for path,digest in pins.items()):raise ValueError("producer changed during confirmation")
    return cut.response._write_exclusive(output,"confirmation",record,arrays,input_paths=tuple(inputs))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--run",action="store_true")
    parser.add_argument("--source",type=Path,default=cut.OUTPUT);parser.add_argument("--output",type=Path,default=OUTPUT)
    parser.add_argument("--producer-commit");parser.add_argument("--cpu-budget",type=float,default=BUDGET)
    args=parser.parse_args();result=run(args.source,args.output,producer_commit=args.producer_commit,cpu_budget=args.cpu_budget) if args.run else preview(args.source)
    print(json.dumps(cut.parent._plain(result),indent=2,allow_nan=False))


if __name__=="__main__":main()
