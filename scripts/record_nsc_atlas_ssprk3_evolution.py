#!/usr/bin/env python3
"""Create/check four SSPRK3 windows using the already bound input caches.

Checking reconstructs saved fields, diagnostics and four final radial stages.
It never repeats a trajectory, source integration or profile-stage extraction.
"""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
import time
import uuid

import numpy as np
from validation_source_successors import verify_saved_bindings, verify_saved_source

import record_nsc_atlas_background_evolution as previous
import diagnose_nsc_atlas_ssprk3_step as local
from recursive_horizons.nsc_atlas_background_ssprk3 import evolve_background_ssprk3, TIME_SCHEME

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-ssprk3-evolution-v1"
RELATIVE_OUTPUT = "results/development/nsc-atlas-ssprk3-evolution-v1"
STATUS = "FINITE_WINDOW_SSPRK3_TEMPORAL_COMPARISON"
SCOPE = ("four original-model windows with unchanged cached data and computational profile; "
         "two temporal resolutions and cross-method differences, not a measured convergence order "
         "or a uniform error enclosure")
LOCAL_INPUT = "results/development/nsc-atlas-ssprk3-step-audit-v1.json"
LOCAL_SHA256 = "9a549c59b9e40a8b9f45207c79ba638cab3409e44bb4c327e4a5b44c0ba9749a"
CASES = tuple((label, 32, steps) for label in ("zero", "finite") for steps in (128, 256))
TUBE = previous.TUBE
SETTINGS = dict(cases=[dict(amplitude=a,cells=n,steps=k) for a,n,k in CASES],
    input=previous.RELATIVE_OUTPUT+"/background.json", input_cache_nodes=513,
    profile_training_rows=65, profile_withheld_rows=64, profile_stage_extractions=0,
    saved_time_nodes=129, exterior_quadrature_order=2, cap_quadrature_order=4,
    maximum_courant=.4, tube=TUBE, global_cpu_seconds=300., global_evaluations=50000000,
    standalone_cpu_seconds=60., standalone_evaluations=5000000,
    time_scheme=TIME_SCHEME, stage_order="left, right, supplied exact-cache midpoint",
    source_continuation=previous.SETTINGS["source_continuation"],
    clock_origin=previous.SETTINGS["clock_origin"], lookup=previous.SETTINGS["lookup"],
    exterior_scheme=previous.SETTINGS["exterior_scheme"], cap_scheme=previous.SETTINGS["cap_scheme"],
    cap_comparison="total physical dr integrals Ip/Iq at common x and equal cap cells",
    counts="actual shared Budget deltas; unsaved-stage point totals are measured, not independently recomputed")
SOURCES = tuple(dict.fromkeys((*previous.SOURCES,
    previous.RELATIVE_OUTPUT+"/background.json", previous.RELATIVE_OUTPUT+"/manifest.json",
    *local.SOURCES, LOCAL_INPUT,
    "scripts/record_nsc_atlas_ssprk3_evolution.py", "tests/test_record_nsc_atlas_ssprk3_evolution.py")))
Budget = previous.Budget
plain, encoded, digest = previous.plain, previous.encoded, previous.digest
statistics, runtime_versions = previous.statistics, previous.runtime_versions
subset, StoredRows = previous.subset, previous.StoredRows
_equal, _array, _count, _real = previous._equal, previous._array, previous._count, previous._real
physical_view, final_summary = previous.physical_view, previous.final_summary
tables, cap_history, fit_profile = previous.tables, previous.cap_history, previous.fit_profile
background_stage, cap_pullback = previous.background_stage, previous.cap_pullback
_histories = previous.previous._histories
compact_history = previous.previous.previous.compact_history
_difference = previous.previous.previous._difference
_difference_summary = previous.previous.previous._difference_summary
quiet_departures = previous.previous.quiet_departures
_failure_safe = previous.previous.previous.previous._failure_safe


def source_bindings():
    return {name:digest((ROOT/name).read_bytes()) for name in SOURCES}


def shared_budget(budget=None):
    if budget is None:
        return Budget(maximum_evaluations=50000000,cpu_seconds=300.)
    if (not isinstance(budget,Budget) or type(budget.maximum_evaluations) is not int
            or not 0<budget.maximum_evaluations<=50000000
            or isinstance(budget.cpu_seconds,bool) or not math.isfinite(budget.cpu_seconds)
            or not 0<budget.cpu_seconds<=300.
            or type(budget.evaluations) is not int or budget.evaluations<0
            or not math.isfinite(budget.started) or budget.started>time.process_time()):
        raise ValueError("invalid shared SSPRK3 allowance")
    budget.check()
    return budget


def load_inputs(budget):
    """Hash-only original66 and local receipts; no inherited scientific replay."""
    budget.check()
    parent,base,_,identity=local.authenticate(ROOT)
    raw=(ROOT/LOCAL_INPUT).read_bytes()
    if digest(raw)!=LOCAL_SHA256:
        raise ValueError("pinned local SSPRK3 receipt mismatch")
    saved=json.loads(raw)
    if (saved.get("schema")!=local.SCHEMA or saved.get("input")!=identity
            or saved.get("limits")!=local.LIMITS or saved.get("runtime")!=runtime_versions()
            or set(saved.get("audit_source_bindings",{}))!=set(local.SOURCES)
            or saved.get("execution",{}).get("stage_calls")!=6):
        raise ValueError("pinned local SSPRK3 declaration mismatch")
    for name,expected in saved["audit_source_bindings"].items():
        verify_saved_source(name, expected, ROOT)
    budget.check()
    return parent,base,dict(background=identity["original"],rk2_local=identity["local_step"],
        ssprk3_local=dict(path=LOCAL_INPUT,sha256=digest(raw),bytes=len(raw),
            source_bindings=saved["audit_source_bindings"],local_steps_replayed=0))


def diagnostics(record,parent,base,profile,budget):
    histories,caps,oracles,current=[],[],[],[]
    for case in record["cases"]:
        full,departure=_histories(base,case,budget)
        coarse,departure_coarse=_histories(base,case,budget,2)
        histories.append(dict(raw129=compact_history(full),departure129=compact_history(departure),
            raw65=compact_history(coarse),departure65=compact_history(departure_coarse)))
        current.append(tables(case,base["initial_data"][case["amplitude"]]["snapshot_rows"],budget))
        caps.append(cap_history(case,profile,budget))
        if case["alpha"]==0:
            oracle=copy.deepcopy(case)
            for sample in oracle["snapshots"]:
                sample["state"]["cap"]=sample["physical"]["cap"]
            oracles.append(cap_pullback(oracle,budget))
    old={}
    for case in parent["cases"]:
        if case["cells"]==32 and case["steps"] in (256,512):
            old[case["amplitude"],case["steps"]]=tables(case,
                base["initial_data"][case["amplitude"]]["snapshot_rows"],budget)
    if set(old)!={(a,k) for a in ("zero","finite") for k in (256,512)}:
        raise ValueError("all four bound SSPRK2 comparison cases required")
    paired=[_difference(current[j+2],current[j]) for j in range(2)]
    oldpaired={k:_difference(old["finite",k],old["zero",k]) for k in (256,512)}
    controls=[]
    for case,table,cap in zip(record["cases"],current,caps):
        x=np.array([s["x"] for s in case["snapshots"]])
        quiet=x<=float(base["geometry"]["support_x"][0])
        quiet_values=np.concatenate([np.asarray(table[group][key])[quiet].ravel()
            for group,key in (("exterior","c"),("wall","c"),("cap","Ip"),("cap","Iq"))])
        departures=np.concatenate([np.asarray(s["state"]["departure"][key]).ravel()
            for s in case["snapshots"] for key in ("r","L")])
        controls.append(dict(quiet_scalar=statistics(quiet_values),metric_departure=statistics(departures),
            minimum_upper_label_margin=min(cap["arrays"]["upper_label_margin"]),
            metric_departure_interpretation="vacuum preservation" if case["alpha"]==0 else "physical response"))
    return plain(dict(histories=histories,cap_histories=caps,zero_cap_pullback=oracles,controls=controls,
        finite_minus_zero=[_difference_summary(current[j+2],current[j]) for j in range(2)],
        time_zero=_difference_summary(current[1],current[0]),
        time_finite=_difference_summary(current[3],current[2]),
        paired_time=_difference_summary(paired[1],paired[0]),
        versus_ssprk2=[{str(k):_difference_summary(table,old[case["amplitude"],k]) for k in (256,512)}
            for case,table in zip(record["cases"],current)],
        paired_versus_ssprk2=[{str(k):_difference_summary(pair,oldpaired[k]) for k in (256,512)} for pair in paired],
        old_ssprk2_paired_time=_difference_summary(oldpaired[512],oldpaired[256]),
        comparison_coordinates="common x/exterior eta/cap y; physical radii may differ; equal spatial resolution",
        temporal_difference_is_measured_order=False,constructed_identities_are_accuracy_tests=False))


def calculate(budget=None):
    budget=shared_budget(budget)
    initial=budget.evaluations
    phase,record={"name":"hash_only_inherited_authentication"},{}
    try:
        parent,base,identity=load_inputs(budget)
        auth_charges=budget.evaluations-initial
        phase={"name":"reuse_computational_profile"}
        start=budget.evaluations
        profile,summary=fit_profile(parent["profile_inputs"],budget)
        profile_charges=budget.evaluations-start
        _equal(summary,parent["profile"],"bound computational profile")
        x=_array(base["gauge"]["rows"]["x"],(513,),"cache coordinates")
        record=dict(schema=SCHEMA,status=STATUS,scope=SCOPE,settings=SETTINGS,
            runtime=runtime_versions(),parent_input=identity,geometry=base["geometry"],
            profile_inputs=parent["profile_inputs"],profile=summary,cases=[],
            refinement_is_rigorous_error_enclosure=False)
        for label,n,steps in CASES:
            phase=dict(name="SSPRK3 evolution",amplitude=label,cells=n,steps=steps)
            stride=512//steps
            source=StoredRows(base["initial_data"][label]["snapshot_rows"])
            gauge=StoredRows(base["gauge"]["rows"])
            start=budget.evaluations
            solution=evolve_background_ssprk3(quiet_departures(exterior_cells=n,cap_cells=n),
                x[::stride],mid_nodes=x[stride//2::stride],eta=np.linspace(0,1,n+1),
                y_faces=np.linspace(0,1,n+1),source=source,gauge=gauge,
                alpha=base["initial_data"][label]["alpha"],tube=TUBE,budget=budget,profile=profile,
                save_every=steps//128,maximum_courant=.4,exterior_quadrature_order=2,cap_quadrature_order=4)
            if solution["steps"]!=steps or solution["time_scheme"]!=TIME_SCHEME:
                raise ValueError("SSPRK3 returned a different step scheme")
            final=final_summary(solution["final_stage"])
            record["cases"].append(plain(dict(amplitude=label,cells=n,steps=steps,
                alpha=base["initial_data"][label]["alpha"],snapshots=solution["snapshots"],final=final,
                execution=dict(budget_start=start,budget_end=budget.evaluations,
                    measured_evaluations=budget.evaluations-start,stages=solution["stages"],
                    source_cache_reads=source.calls,gauge_cache_reads=gauge.calls,
                    max_courant=solution["max_courant"],retries=solution["retries"],
                    final_stage_evaluations=8*(n+1)+2+final["cap_metric_evaluations"]+final["cap_rhs_evaluations"]))))
        phase={"name":"saved_array_diagnostics"}
        start=budget.evaluations
        record["diagnostics"]=diagnostics(record,parent,base,profile,budget)
        record["execution"]=dict(initial_evaluations=initial,evaluations=budget.evaluations,
            cpu_seconds=time.process_time()-budget.started,inherited_authentication_evaluations=auth_charges,
            profile_build_evaluations=profile_charges,diagnostic_evaluations=budget.evaluations-start,
            old_profile_stage_reconstructions=0,coupled_evolutions=len(CASES),
            time_steps=sum(k for _,_,k in CASES),stage_evaluations=sum(3*k+1 for _,_,k in CASES),
            source_integrations=0,gauge_quadratures=0,new_physical_preparations=0,
            stage_point_totals_are_measured=True)
        budget.check()
        return plain(record)
    except Exception as error:
        error.ssprk3_failure=plain(dict(phase=phase,exception=type(error).__name__,message=str(error),
            stage=getattr(error,"diagnostic",None),initial_evaluations=initial,evaluations=budget.evaluations,
            cpu_seconds=time.process_time()-budget.started,completed_cases=record.get("cases",[]),
            parent_input=record.get("parent_input")))
        raise


def validate_record(record,budget=None,*,require_creation=True):
    budget=shared_budget(budget)
    start=budget.evaluations
    try:
        if (record["schema"]!=SCHEMA or record["status"]!=STATUS or record["scope"]!=SCOPE
                or record["refinement_is_rigorous_error_enclosure"] is not False):
            raise ValueError("SSPRK3 declaration mismatch")
        _equal(record["settings"],SETTINGS,"settings")
        _equal(record["runtime"],runtime_versions(),"runtime")
        parent,base,identity=load_inputs(budget)
        auth_charges=budget.evaluations-start
        _equal(record["parent_input"],identity,"parent input")
        _equal(record["geometry"],base["geometry"],"geometry")
        _equal(record["profile_inputs"],parent["profile_inputs"],"unchanged profile inputs")
        before=budget.evaluations
        profile,summary=fit_profile(record["profile_inputs"],budget)
        profile_charges=budget.evaluations-before
        _equal(record["profile"],summary,"profile fit and withheld values")
        _equal(record["profile"],parent["profile"],"bound computational profile")
        if len(record["cases"])!=len(CASES):
            raise ValueError("all four ordered SSPRK3 cases required")
        x=_array(base["gauge"]["rows"]["x"],(513,),"cache coordinates")
        initial=record["execution"]["initial_evaluations"]
        _count(initial,50000000,"initial shared charges")
        end=initial+auth_charges+profile_charges
        for case,(label,n,steps) in zip(record["cases"],CASES):
            alpha=base["initial_data"][label]["alpha"]
            if (case["amplitude"],case["cells"],case["steps"],case["alpha"])!=(label,n,steps,alpha):
                raise ValueError("case identity mismatch")
            if len(case["snapshots"])!=129:
                raise ValueError("require129 saved slices")
            eta=np.linspace(0,1,n+1)
            for j,(sample,at) in enumerate(zip(case["snapshots"],x[::4])):
                if sample["x"]!=at or set(sample)!={"x","state","physical"}:
                    raise ValueError("saved row/time inventory mismatch")
                state=sample["state"]
                if set(state)!={"departure","wall","cap"}:
                    raise ValueError("raw state inventory mismatch")
                for group,names,shape in (("departure",("r","L","c"),(n-1,)),
                    ("wall",("r","L","T","s"),()),("cap",("Ip","Iq"),(n,))):
                    if set(state[group])!=set(names):
                        raise ValueError("raw state fields mismatch")
                    for key in names:
                        _array(state[group][key],shape,group+"."+key)
                expected=physical_view(state,at,subset(base["initial_data"][label]["snapshot_rows"],4*j),
                    subset(base["gauge"]["rows"],4*j),eta,alpha,profile,budget)
                _equal(sample["physical"],expected,"profile/raw-to-physical reconstruction")
            _equal(case["snapshots"][0]["state"],quiet_departures(exterior_cells=n,cap_cells=n),"quiet departures")
            for key in ("T","s"):
                clocks=np.array([s["physical"]["wall"][key] for s in case["snapshots"]])
                if clocks[0]!=0 or np.any(np.diff(clocks)<=0):
                    raise ValueError("clocks must increase from quiet zero")
            before=budget.evaluations
            stage=background_stage(case["snapshots"][-1]["state"],x=x[-1],eta=eta,y_faces=eta,
                source=subset(base["initial_data"][label]["snapshot_rows"],-1),gauge=subset(base["gauge"]["rows"],-1),
                alpha=alpha,tube=TUBE,budget=budget,profile=profile,
                exterior_quadrature_order=2,cap_quadrature_order=4)
            final_charge=budget.evaluations-before
            _equal(case["final"],final_summary(stage),"final cap/prefix/wall reconstruction")
            execution=case["execution"]
            charge=execution["measured_evaluations"]
            _count(charge,50000000,"measured case charges")
            if charge<(3*steps+1)*(8*(n+1)+2):
                raise ValueError("case charge below exterior/hybrid minimum")
            expected=dict(budget_start=end,budget_end=end+charge,measured_evaluations=charge,
                stages=3*steps+1,source_cache_reads=3*steps+1,gauge_cache_reads=3*steps+1,retries=0,
                final_stage_evaluations=final_charge)
            _equal({k:v for k,v in execution.items() if k!="max_courant"},expected,"measured case accounting")
            _real(execution["max_courant"],"Courant",high=.4)
            end+=charge
        before=budget.evaluations
        _equal(record["diagnostics"],diagnostics(record,parent,base,profile,budget),"saved diagnostics")
        diagnostic_charges=budget.evaluations-before
        expected=dict(initial_evaluations=initial,evaluations=end+diagnostic_charges,
            inherited_authentication_evaluations=auth_charges,profile_build_evaluations=profile_charges,
            diagnostic_evaluations=diagnostic_charges,old_profile_stage_reconstructions=0,coupled_evolutions=len(CASES),
            time_steps=sum(k for _,_,k in CASES),stage_evaluations=sum(3*k+1 for _,_,k in CASES),
            source_integrations=0,gauge_quadratures=0,new_physical_preparations=0,stage_point_totals_are_measured=True)
        _equal({k:v for k,v in record["execution"].items() if k!="cpu_seconds"},expected,"execution accounting")
        _count(record["execution"]["evaluations"],50000000,"calculation charges")
        _real(record["execution"]["cpu_seconds"],"calculation CPU",high=300.)
        budget.check()
        charges=budget.evaluations-start
        if require_creation:
            creation=record["creation_execution"]
            _equal({k:v for k,v in creation.items() if k!="cpu_seconds"},
                dict(calculation_evaluations=record["execution"]["evaluations"],validation_evaluations=charges,
                    total_evaluations=record["execution"]["evaluations"]+charges),"creation accounting")
            _count(creation["total_evaluations"],50000000,"creation charges")
            _real(creation["cpu_seconds"],"creation CPU",low=record["execution"]["cpu_seconds"],high=300.)
    except (KeyError,TypeError,IndexError,OverflowError) as error:
        raise ValueError("incomplete SSPRK3 record") from error
    return dict(authentication_evaluations=budget.evaluations-start,source_integrations=0,time_evolutions=0,
        reference_quadratures=0,old_profile_stage_spotchecks=0,profile_stage_extraction_replayed=False,
        current_final_radial_stages=len(CASES),unsaved_stage_point_totals_independently_recomputed=False)


def create(directory):
    directory=Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only record: choose a new directory or --check")
    started=time.process_time()
    budget=Budget(maximum_evaluations=50000000,cpu_seconds=300.,started=started)
    bindings,runtime=source_bindings(),runtime_versions()
    budget.check()
    try:
        record=calculate(budget)
        if bindings!=source_bindings() or runtime!=runtime_versions() or record["runtime"]!=runtime:
            raise RuntimeError("SSPRK3 source/runtime changed during calculation")
        auth=validate_record(record,budget,require_creation=False)
        if bindings!=source_bindings() or runtime!=runtime_versions():
            raise RuntimeError("SSPRK3 source/runtime changed during authentication")
        record["creation_execution"]=dict(calculation_evaluations=record["execution"]["evaluations"],
            validation_evaluations=auth["authentication_evaluations"],total_evaluations=budget.evaluations,
            cpu_seconds=time.process_time()-started)
        budget.check()
    except Exception as error:
        receipt=dict(schema=SCHEMA+"-failed-attempt",status="FAILED_NOT_COMPLETED",settings=SETTINGS,
            source_bindings_before=bindings,runtime=runtime,
            failure=getattr(error,"ssprk3_failure",dict(exception=type(error).__name__,message=str(error),
                evaluations=budget.evaluations,cpu_seconds=time.process_time()-started)))
        path=ROOT/".build/ssprk3-evolution-failures"/(uuid.uuid4().hex+".json")
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open("xb") as stream:
            stream.write(encoded(_failure_safe(receipt)))
        error.failure_receipt=str(path)
        raise
    raw=encoded(record)
    manifest=dict(schema=SCHEMA+"-manifest",source_bindings=bindings,
        outputs={"ssprk3.json":dict(sha256=digest(raw),bytes=len(raw))})
    directory.mkdir(parents=True,exist_ok=False)
    (directory/"ssprk3.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return dict(created=str(directory),execution=record["execution"],creation_execution=record["creation_execution"])


def check(directory):
    budget=Budget(maximum_evaluations=5000000,cpu_seconds=60.)
    directory=Path(directory)
    if ({p.name for p in directory.iterdir()}!={"ssprk3.json","manifest.json"}
            or any(p.is_symlink() for p in directory.iterdir())):
        raise ValueError("SSPRK3 output inventory mismatch")
    manifest=json.loads((directory/"manifest.json").read_bytes())
    raw=(directory/"ssprk3.json").read_bytes()
    bindings,runtime=source_bindings(),runtime_versions()
    if (manifest.get("schema")!=SCHEMA+"-manifest" or not verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=SOURCES)
            or manifest.get("outputs")!={"ssprk3.json":dict(sha256=digest(raw),bytes=len(raw))}):
        raise ValueError("SSPRK3 hash/source mismatch")
    budget.check()
    authentication=validate_record(json.loads(raw),budget)
    if bindings!=source_bindings() or runtime!=runtime_versions():
        raise RuntimeError("SSPRK3 source/runtime changed during authentication")
    budget.check()
    return dict(checked=str(directory),**authentication)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    try:
        result=check(args.output) if args.check else create(args.output)
    except Exception as error:
        if hasattr(error,"failure_receipt"):
            print(json.dumps(dict(failed=str(error),diagnostic=error.failure_receipt)))
        raise
    print(json.dumps(result,sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()
