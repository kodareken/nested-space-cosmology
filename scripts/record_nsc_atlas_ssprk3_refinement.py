#!/usr/bin/env python3
"""Create/check two midpoint-complete SSPRK3 windows without old trajectory replay."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import time
import uuid

import numpy as np
from recursive_horizons.nsc_atlas_validation import collect_departure_replay, profile_interval_replay_equal
from validation_source_successors import verify_saved_bindings, verify_saved_source

import record_nsc_atlas_ssprk3_evolution as previous
from recursive_horizons.nsc_atlas_refined_cache import make_refined_cache, validate_refined_cache

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-ssprk3-refinement-v1"
RELATIVE_OUTPUT = "results/development/nsc-atlas-ssprk3-refinement-v1"
STATUS = "MIDPOINT_COMPLETE_SSPRK3_TEMPORAL_COMPARISON"
SCOPE = ("two original-model 512-step windows using the unchanged source snapshot, profile and spatial operators; "
         "sampled temporal and cross-method differences, not a convergence order or uniform error enclosure")
INPUT = previous.RELATIVE_OUTPUT+"/ssprk3.json"
INPUT_SHA256 = "24cd8319c131360154cbdd7088cd2f212f8d572a0fe6fdff69a2c845c9223546"
MANIFEST_SHA256 = "cf20e9190d400e81fa05a1dcf31535cbe99ed20d9a3fdf9c21f981f2565fdf75"
CASES = (("zero",32,512),("finite",32,512))
SETTINGS = dict(input=INPUT,cases=[dict(amplitude=a,cells=n,steps=k) for a,n,k in CASES],
    original_cache_nodes=513,refined_cache_nodes=1025,saved_time_nodes=129,
    profile_training_rows=65,profile_withheld_rows=64,profile_stage_extractions=0,
    exterior_quadrature_order=2,cap_quadrature_order=4,maximum_courant=.4,
    tube=previous.TUBE,time_scheme=previous.TIME_SCHEME,
    global_cpu_seconds=300.,global_evaluations=50000000,
    standalone_cpu_seconds=60.,standalone_evaluations=5000000,
    source_continuation=previous.SETTINGS["source_continuation"],
    clock_origin=previous.SETTINGS["clock_origin"],
    lookup="original exact keys at even indices; explicit new midpoint keys at odd indices",
    counts="measured cumulative deltas; unsaved trajectory point totals are not independently reconstructed",
    cap_comparison="total physical dr integrals at common x and equal cap cells")
SOURCES = tuple(dict.fromkeys((*previous.SOURCES,INPUT,previous.RELATIVE_OUTPUT+"/manifest.json",
    "src/recursive_horizons/nsc_atlas_refined_cache.py","tests/test_nsc_atlas_refined_cache.py",
    "scripts/record_nsc_atlas_ssprk3_refinement.py","tests/test_record_nsc_atlas_ssprk3_refinement.py")))
profile_replay_equal = previous.profile_replay_equal
Budget,shared_budget = previous.Budget,previous.shared_budget
plain,encoded,digest = previous.plain,previous.encoded,previous.digest
runtime_versions,statistics = previous.runtime_versions,previous.statistics
_equal,_array,_count,_real = previous._equal,previous._array,previous._count,previous._real
subset,StoredRows = previous.subset,previous.StoredRows
physical_view,final_summary = previous.physical_view,previous.final_summary
tables,cap_history,cap_pullback = previous.tables,previous.cap_history,previous.cap_pullback
_histories,compact_history = previous._histories,previous.compact_history
fit_profile,background_stage = previous.fit_profile,previous.background_stage
evolve_background_ssprk3,quiet_departures = previous.evolve_background_ssprk3,previous.quiet_departures
_difference,_failure_safe = previous._difference,previous._failure_safe


def source_bindings():
    return {name:digest((ROOT/name).read_bytes()) for name in SOURCES}


def load_inputs(budget):
    """Authenticate pinned bytes and old bindings; call no scientific validator."""
    budget.check()
    directory=ROOT/previous.RELATIVE_OUTPUT
    if ({p.name for p in directory.iterdir()}!={"ssprk3.json","manifest.json"}
            or any(p.is_symlink() for p in directory.iterdir())):
        raise ValueError("inherited SSPRK3 inventory mismatch")
    raw=(directory/"ssprk3.json").read_bytes()
    manifest_raw=(directory/"manifest.json").read_bytes()
    manifest=json.loads(manifest_raw)
    bindings=manifest.get("source_bindings",{})
    if (digest(raw)!=INPUT_SHA256 or digest(manifest_raw)!=MANIFEST_SHA256
            or manifest.get("schema")!=previous.SCHEMA+"-manifest"
            or manifest.get("outputs")!={"ssprk3.json":dict(sha256=digest(raw),bytes=len(raw))}
            or set(bindings)!=set(previous.SOURCES) or len(bindings)!=80):
        raise ValueError("pinned SSPRK3 payload/manifest mismatch")
    for name,expected in bindings.items():
        verify_saved_source(name, expected, ROOT)
    parent=json.loads(raw)
    background,base,identity=previous.load_inputs(budget)
    if (parent.get("schema")!=previous.SCHEMA or parent.get("settings")!=previous.SETTINGS
            or parent.get("runtime")!=runtime_versions() or parent.get("parent_input")!=identity):
        raise ValueError("inherited SSPRK3 declaration mismatch")
    budget.check()
    return parent,background,base,dict(path=INPUT,sha256=digest(raw),manifest_sha256=digest(manifest_raw),
        authenticated_input_files=80,inherited=identity,trajectory_replays=0,source_integrations=0)


def comparison(a,b,x,budget):
    """Full signed arrays and maxima at saved coordinates, without an error claim."""
    arrays=_difference(a,b)
    summaries,locations={},{}
    for group,fields in arrays.items():
        summaries[group],locations[group]={},{}
        for key,value in fields.items():
            value=np.asarray(value)
            budget.check(value.size)
            if value.shape[0]!=len(x) or not np.all(np.isfinite(value)):
                raise ValueError("comparison arrays must share finite saved cuts")
            index=np.unravel_index(int(np.argmax(abs(value))),value.shape)
            summaries[group][key]=statistics(value)
            locations[group][key]=dict(index=list(index),x=float(x[index[0]]),signed=float(value[index]))
    return plain(dict(arrays=arrays,statistics=summaries,max_locations=locations))


def profile_intervals(case,profile,budget):
    """Brackets of saved wall labels; crossings are not attributed as causal errors."""
    rows=case["snapshots"]
    x=np.array([s["x"] for s in rows])
    r=np.array([s["physical"]["wall"]["r"] for s in rows])
    T=np.array([s["physical"]["wall"]["T"] for s in rows])
    if np.any(abs(r)>=1):
        raise ValueError("wall profile label requires abs(r)<1")
    q=T+np.arctanh(r)
    knots=np.asarray(profile.breakpoints)
    if np.any(np.diff(q)<=0) or np.max(q)>knots[-1]+8*np.finfo(float).eps*max(1.,abs(knots[-1])):
        raise ValueError("saved profile labels must be ordered and in range")
    left=np.searchsorted(knots,q,side="right")-1
    right=np.searchsorted(knots,q,side="left")
    crossed=np.maximum(0,np.searchsorted(knots,q[1:],side="left")-
                        np.searchsorted(knots,q[:-1],side="right"))
    budget.check(5*len(q)+len(crossed))
    return plain(dict(x=x,q=q,left_knot_index=left,right_knot_index=right,
        strictly_interior_knots_per_saved_interval=crossed,
        interpretation="saved-interval knot brackets, not internal-stage paths or causal attribution"))


def diagnostics(record,parent,background,base,profile,budget, *, replay_allowances=None):
    current,histories,caps,controls=[],[],[],[]
    for case in record["cases"]:
        full,departure=_histories(base,case,budget)
        coarse,departure_coarse=_histories(base,case,budget,2)
        histories.append(dict(raw129=compact_history(full),departure129=compact_history(departure),
            raw65=compact_history(coarse),departure65=compact_history(departure_coarse)))
        collect_departure_replay(replay_allowances, len(histories)-1,
            base["initial_data"][case["amplitude"]]["snapshot_rows"], case["alpha"],
            departure129=departure, departure65=departure_coarse)
        table=tables(case,base["initial_data"][case["amplitude"]]["snapshot_rows"],budget)
        current.append(table)
        caps.append(cap_history(case,profile,budget))
        x=np.array([s["x"] for s in case["snapshots"]])
        quiet=x<=base["geometry"]["support_x"][0]
        channels={group+"."+key:statistics(np.asarray(table[group][key])[quiet])
            for group,key in (("exterior","c"),("wall","c"),("cap","Ip"),("cap","Iq"))}
        departures=np.concatenate([np.asarray(s["state"]["departure"][key]).ravel()
            for s in case["snapshots"] for key in ("r","L")])
        budget.check(sum(np.asarray(table[g][k])[quiet].size for g,k in
            (("exterior","c"),("wall","c"),("cap","Ip"),("cap","Iq")))+departures.size)
        controls.append(dict(quiet_scalar_channels=channels,metric_departure=statistics(departures),
            minimum_upper_label_margin=min(caps[-1]["arrays"]["upper_label_margin"])))
    retained={(c["amplitude"],c["steps"]):tables(c,base["initial_data"][c["amplitude"]]["snapshot_rows"],budget)
        for c in parent["cases"]}
    if set(retained)!={(a,k) for a in ("zero","finite") for k in (128,256)}:
        raise ValueError("all four retained SSPRK3 cases required")
    rk2={(c["amplitude"],c["steps"]):tables(c,base["initial_data"][c["amplitude"]]["snapshot_rows"],budget)
        for c in background["cases"] if c["cells"]==32 and c["steps"] in (256,512)}
    if set(rk2)!={(a,k) for a in ("zero","finite") for k in (256,512)}:
        raise ValueError("all four old SSPRK2 comparison cases required")
    paired=_difference(current[1],current[0])
    oldpaired={k:_difference(retained["finite",k],retained["zero",k]) for k in (128,256)}
    rkpaired={k:_difference(rk2["finite",k],rk2["zero",k]) for k in (256,512)}
    extrapolant={g:{k:(4*np.asarray(v)-np.asarray(rkpaired[256][g][k]))/3
        for k,v in fields.items()} for g,fields in rkpaired[512].items()}
    budget.check(sum(np.size(v) for fields in extrapolant.values() for v in fields.values()))
    oracle=copy.deepcopy(record["cases"][0])
    for sample in oracle["snapshots"]:
        sample["state"]["cap"]=sample["physical"]["cap"]
    intervals=[profile_intervals(c,profile,budget) for c in (*parent["cases"],*record["cases"])]
    return plain(dict(histories=histories,cap_histories=caps,controls=controls,
        zero_cap_pullback=cap_pullback(oracle,budget),
        finite_minus_zero=comparison(current[1],current[0],x,budget),
        time_zero=comparison(current[0],retained["zero",256],x,budget),
        time_finite=comparison(current[1],retained["finite",256],x,budget),
        paired_time=comparison(paired,oldpaired[256],x,budget),
        retained_128_256=dict(zero=comparison(retained["zero",256],retained["zero",128],x,budget),
            finite=comparison(retained["finite",256],retained["finite",128],x,budget),
            paired=comparison(oldpaired[256],oldpaired[128],x,budget)),
        versus_128=dict(zero=comparison(current[0],retained["zero",128],x,budget),
            finite=comparison(current[1],retained["finite",128],x,budget),
            paired=comparison(paired,oldpaired[128],x,budget)),
        old_rk2_paired_extrapolant=extrapolant,
        paired_versus_old_rk2_extrapolant=comparison(paired,extrapolant,x,budget),
        profile_knots=profile.breakpoints,profile_intervals=intervals,
        profile_case_order=[dict(amplitude=c["amplitude"],steps=c["steps"]) for c in (*parent["cases"],*record["cases"])],
        comparison_coordinates="common saved x/exterior eta/cap y; physical radii may differ",
        extrapolant_is_exact_reference=False,temporal_order_established=False,uniform_error_enclosure=False))


def calculate(budget=None):
    budget=shared_budget(budget)
    initial=budget.evaluations
    phase,record={"name":"hash_only_input_authentication"},{}
    try:
        parent,background,base,identity=load_inputs(budget)
        auth_charges=budget.evaluations-initial
        phase={"name":"midpoint_cache"}
        before=budget.evaluations
        cache=make_refined_cache(base,budget)
        cache_charges=budget.evaluations-before
        phase={"name":"unchanged_profile"}
        before=budget.evaluations
        profile,summary=fit_profile(parent["profile_inputs"],budget)
        profile_charges=budget.evaluations-before
        profile_replay_equal(summary,parent["profile"],parent["profile_inputs"])
        record=dict(schema=SCHEMA,status=STATUS,scope=SCOPE,settings=SETTINGS,runtime=runtime_versions(),
            parent_input=identity,geometry=base["geometry"],cache=cache,
            profile_inputs=parent["profile_inputs"],profile=summary,cases=[],uniform_error_enclosure=False)
        x=_array(cache["x_nodes"],(1025,),"refined cache coordinates")
        for label,n,steps in CASES:
            phase=dict(name="new SSPRK3 evolution",amplitude=label,cells=n,steps=steps)
            source,gauge=StoredRows(cache["source_rows"][label]),StoredRows(cache["gauge_rows"])
            before=budget.evaluations
            solution=evolve_background_ssprk3(quiet_departures(exterior_cells=n,cap_cells=n),x[::2],
                mid_nodes=x[1::2],eta=np.linspace(0,1,n+1),y_faces=np.linspace(0,1,n+1),
                source=source,gauge=gauge,alpha=base["initial_data"][label]["alpha"],tube=previous.TUBE,
                budget=budget,profile=profile,save_every=4,maximum_courant=.4,
                exterior_quadrature_order=2,cap_quadrature_order=4)
            if solution["steps"]!=steps or solution["time_scheme"]!=previous.TIME_SCHEME:
                raise ValueError("refinement returned a different time scheme")
            final=final_summary(solution["final_stage"])
            record["cases"].append(plain(dict(amplitude=label,cells=n,steps=steps,
                alpha=base["initial_data"][label]["alpha"],snapshots=solution["snapshots"],final=final,
                execution=dict(budget_start=before,budget_end=budget.evaluations,
                    measured_evaluations=budget.evaluations-before,stages=solution["stages"],
                    source_cache_reads=source.calls,gauge_cache_reads=gauge.calls,
                    max_courant=solution["max_courant"],retries=solution["retries"],
                    final_stage_evaluations=8*(n+1)+2+final["cap_metric_evaluations"]+final["cap_rhs_evaluations"]))))
        phase={"name":"saved_array_diagnostics"}
        before=budget.evaluations
        record["diagnostics"]=diagnostics(record,parent,background,base,profile,budget)
        record["execution"]=dict(initial_evaluations=initial,evaluations=budget.evaluations,
            cpu_seconds=time.process_time()-budget.started,inherited_authentication_evaluations=auth_charges,
            cache_evaluations=cache_charges,profile_build_evaluations=profile_charges,
            diagnostic_evaluations=budget.evaluations-before,coupled_evolutions=len(CASES),
            time_steps=sum(k for _,_,k in CASES),stage_evaluations=sum(3*k+1 for _,_,k in CASES),
            source_integrations=0,old_trajectory_replays=0,profile_stage_extractions=0,
            stage_point_totals_are_measured=True)
        budget.check()
        return plain(record)
    except Exception as error:
        error.refinement_failure=_failure_safe(dict(phase=phase,exception=type(error).__name__,message=str(error),
            stage=getattr(error,"diagnostic",None),initial_evaluations=initial,evaluations=budget.evaluations,
            cpu_seconds=time.process_time()-budget.started,completed_cases=record.get("cases",[]),
            parent_input=record.get("parent_input")))
        raise


def validate_record(record,budget=None,*,require_creation=True):
    budget=shared_budget(budget)
    start=budget.evaluations
    try:
        if (record["schema"]!=SCHEMA or record["status"]!=STATUS or record["scope"]!=SCOPE
                or record["uniform_error_enclosure"] is not False):
            raise ValueError("refinement declaration mismatch")
        _equal(record["settings"],SETTINGS,"settings")
        _equal(record["runtime"],runtime_versions(),"runtime")
        parent,background,base,identity=load_inputs(budget)
        auth_charges=budget.evaluations-start
        _equal(record["parent_input"],identity,"parent input")
        _equal(record["geometry"],base["geometry"],"geometry")
        before=budget.evaluations
        cache_auth=validate_refined_cache(record["cache"],base,budget)
        replay_cache_charges=budget.evaluations-before
        cache_charges=cache_auth["saved_cache_evaluations"]
        _count(cache_charges,50000000,"saved cache charges")
        _count(cache_auth["validation_evaluations"],50000000,"cache replay charges")
        if (cache_auth["validation_evaluations"]!=replay_cache_charges
                or cache_charges!=record["cache"]["execution"]["evaluations"]):
            raise ValueError("cache saved/replay accounting mismatch")
        # An adaptive inverse may take different valid iteration counts on a
        # different platform. Historical totals retain its authenticated saved
        # cost; the live Budget has already charged all actual replay work.
        _equal(record["profile_inputs"],parent["profile_inputs"],"profile inputs")
        before=budget.evaluations
        profile,summary=fit_profile(record["profile_inputs"],budget)
        profile_charges=budget.evaluations-before
        profile_replay_equal(record["profile"],summary,record["profile_inputs"])
        profile_replay_equal(summary,parent["profile"],parent["profile_inputs"])
        if len(record["cases"])!=len(CASES):
            raise ValueError("both ordered refinement cases required")
        x=_array(record["cache"]["x_nodes"],(1025,),"refined cache coordinates")
        initial=record["execution"]["initial_evaluations"]
        _count(initial,50000000,"initial shared charges")
        end=initial+auth_charges+cache_charges+profile_charges
        for case,(label,n,steps) in zip(record["cases"],CASES):
            alpha=base["initial_data"][label]["alpha"]
            if (case["amplitude"],case["cells"],case["steps"],case["alpha"])!=(label,n,steps,alpha):
                raise ValueError("case identity mismatch")
            if len(case["snapshots"])!=129:
                raise ValueError("require129 saved slices")
            eta=np.linspace(0,1,n+1)
            for j,(sample,at) in enumerate(zip(case["snapshots"],x[::8])):
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
                expected=physical_view(state,at,subset(record["cache"]["source_rows"][label],8*j),
                    subset(record["cache"]["gauge_rows"],8*j),eta,alpha,profile,budget)
                _equal(sample["physical"],expected,"raw-to-physical reconstruction")
            _equal(case["snapshots"][0]["state"],quiet_departures(exterior_cells=n,cap_cells=n),"quiet departures")
            for key in ("T","s"):
                clocks=np.array([s["physical"]["wall"][key] for s in case["snapshots"]])
                if clocks[0]!=0 or np.any(np.diff(clocks)<=0):
                    raise ValueError("clocks must increase from zero")
            before=budget.evaluations
            stage=background_stage(case["snapshots"][-1]["state"],x=x[-1],eta=eta,y_faces=eta,
                source=subset(record["cache"]["source_rows"][label],-1),
                gauge=subset(record["cache"]["gauge_rows"],-1),alpha=alpha,tube=previous.TUBE,
                budget=budget,profile=profile,exterior_quadrature_order=2,cap_quadrature_order=4)
            final_charge=budget.evaluations-before
            _equal(case["final"],final_summary(stage),"final stage")
            execution=case["execution"]
            charge=execution["measured_evaluations"]
            _count(charge,50000000,"measured case points")
            if charge<(3*steps+1)*(8*(n+1)+2):
                raise ValueError("case charge below minimum")
            expected=dict(budget_start=end,budget_end=end+charge,measured_evaluations=charge,
                stages=3*steps+1,source_cache_reads=3*steps+1,gauge_cache_reads=3*steps+1,retries=0,
                final_stage_evaluations=final_charge)
            _equal({k:v for k,v in execution.items() if k!="max_courant"},expected,"case accounting")
            _real(execution["max_courant"],"Courant",high=.4)
            end+=charge
        before=budget.evaluations
        replay_allowances = {}
        rebuilt = diagnostics(record,parent,background,base,profile,budget,replay_allowances=replay_allowances)
        saved_diagnostics = record["diagnostics"]
        profile_interval_replay_equal(saved_diagnostics["profile_intervals"],saved_diagnostics["profile_knots"],
                                      rebuilt["profile_intervals"],rebuilt["profile_knots"])
        _equal({k:v for k,v in saved_diagnostics.items() if k!="profile_intervals"},
               {k:v for k,v in rebuilt.items() if k!="profile_intervals"},"diagnostics",roundoff=replay_allowances)
        diagnostic_charges=budget.evaluations-before
        expected=dict(initial_evaluations=initial,evaluations=end+diagnostic_charges,
            inherited_authentication_evaluations=auth_charges,cache_evaluations=cache_charges,
            profile_build_evaluations=profile_charges,diagnostic_evaluations=diagnostic_charges,
            coupled_evolutions=len(CASES),time_steps=sum(k for _,_,k in CASES),
            stage_evaluations=sum(3*k+1 for _,_,k in CASES),source_integrations=0,
            old_trajectory_replays=0,profile_stage_extractions=0,stage_point_totals_are_measured=True)
        _equal({k:v for k,v in record["execution"].items() if k!="cpu_seconds"},expected,"execution accounting")
        _count(record["execution"]["evaluations"],50000000,"calculation points")
        _real(record["execution"]["cpu_seconds"],"calculation CPU",high=300.)
        budget.check()
        charges=budget.evaluations-start
        historical_charges=charges-replay_cache_charges+cache_charges
        if require_creation:
            creation=record["creation_execution"]
            _equal({k:v for k,v in creation.items() if k!="cpu_seconds"},
                dict(calculation_evaluations=record["execution"]["evaluations"],validation_evaluations=historical_charges,
                    total_evaluations=record["execution"]["evaluations"]+historical_charges),"creation accounting")
            _count(creation["total_evaluations"],50000000,"creation points")
            _real(creation["cpu_seconds"],"creation CPU",low=record["execution"]["cpu_seconds"],high=300.)
    except (KeyError,TypeError,IndexError,OverflowError) as error:
        raise ValueError("incomplete refinement record") from error
    return dict(authentication_evaluations=budget.evaluations-start,source_integrations=0,time_evolutions=0,
        cache_replay_evaluations=replay_cache_charges,saved_cache_evaluations=cache_charges,
        historical_equivalent_authentication_evaluations=historical_charges,
        midpoint_cache_reconstructions=1,gauge_quadrature_orders=[64,32],current_final_radial_stages=len(CASES),
        unsaved_stage_point_totals_independently_recomputed=False)


def create(directory):
    started=time.process_time()
    budget=Budget(maximum_evaluations=50000000,cpu_seconds=300.,started=started)
    directory=Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only record: choose a new directory or --check")
    bindings,runtime=source_bindings(),runtime_versions()
    try:
        record=calculate(budget)
        if bindings!=source_bindings() or runtime!=runtime_versions() or record["runtime"]!=runtime:
            raise RuntimeError("refinement source/runtime changed during calculation")
        auth=validate_record(record,budget,require_creation=False)
        if bindings!=source_bindings() or runtime!=runtime_versions():
            raise RuntimeError("refinement source/runtime changed during authentication")
        record["creation_execution"]=dict(calculation_evaluations=record["execution"]["evaluations"],
            validation_evaluations=auth["authentication_evaluations"],total_evaluations=budget.evaluations,
            cpu_seconds=time.process_time()-started)
        budget.check()
    except Exception as error:
        receipt=dict(schema=SCHEMA+"-failed-attempt",status="FAILED_NOT_COMPLETED",settings=SETTINGS,
            source_bindings_before=bindings,runtime=runtime,
            failure=getattr(error,"refinement_failure",dict(exception=type(error).__name__,message=str(error),
                evaluations=budget.evaluations,cpu_seconds=time.process_time()-started)))
        path=ROOT/".build/ssprk3-refinement-failures"/(uuid.uuid4().hex+".json")
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open("xb") as stream:
            stream.write(encoded(_failure_safe(receipt)))
        error.failure_receipt=str(path)
        raise
    raw=encoded(record)
    manifest=dict(schema=SCHEMA+"-manifest",source_bindings=bindings,
        outputs={"refinement.json":dict(sha256=digest(raw),bytes=len(raw))})
    directory.mkdir(parents=True,exist_ok=False)
    (directory/"refinement.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return dict(created=str(directory),execution=record["execution"],creation_execution=record["creation_execution"])


def check(directory):
    budget=Budget(maximum_evaluations=5000000,cpu_seconds=60.)
    directory=Path(directory)
    if ({p.name for p in directory.iterdir()}!={"refinement.json","manifest.json"}
            or any(p.is_symlink() for p in directory.iterdir())):
        raise ValueError("refinement output inventory mismatch")
    raw=(directory/"refinement.json").read_bytes()
    manifest=json.loads((directory/"manifest.json").read_bytes())
    bindings,runtime=source_bindings(),runtime_versions()
    if (manifest.get("schema")!=SCHEMA+"-manifest" or not verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=SOURCES)
            or manifest.get("outputs")!={"refinement.json":dict(sha256=digest(raw),bytes=len(raw))}):
        raise ValueError("refinement hash/source mismatch")
    budget.check()
    result=validate_record(json.loads(raw),budget)
    if bindings!=source_bindings() or runtime!=runtime_versions():
        raise RuntimeError("refinement source/runtime changed during authentication")
    budget.check()
    return dict(checked=str(directory),**result)


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
