#!/usr/bin/env python3
"""Create/check six coupled runs with one retained characteristic background.

Creation extracts old zero-case first jets once. Authentication reuses that
bound table and spot-checks its final stage; it never evolves a source or a
trajectory. Cap integrals in comparisons are total physical integrals.
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

import record_nsc_atlas_compatible_evolution as previous
from diagnose_nsc_atlas_compatible_transport import cap_pullback
from recursive_horizons.nsc_atlas_background_profile import build_background_profile
from recursive_horizons.nsc_atlas_background_evolution import background_stage, evolve_background
from recursive_horizons.nsc_atlas_background_cap import background_cap_metric
from recursive_horizons.nsc_atlas_reduced_cap import _reconstruct, _wall_conversion
from recursive_horizons.nsc_atlas_vacuum_departure import vacuum_slice
from recursive_horizons.nsc_atlas_validation import background_profile_roundoff_tolerances

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-background-evolution-v1"
RELATIVE_OUTPUT = "results/development/nsc-atlas-background-evolution-v1"
STATUS = "FINITE_WINDOW_WITH_CHARACTERISTIC_CAP_BACKGROUND"
SCOPE = ("one computational C2 profile and evolved cap departures on the original coupled window; "
         "finite-grid comparisons and residuals, not a uniform error enclosure or central-arrival evolution")
CASES = tuple((label, n, steps) for label in ("zero", "finite")
              for n, steps in ((16, 256), (32, 256), (32, 512)))
TUBE = previous.TUBE
SETTINGS = dict(cases=[dict(amplitude=a, cells=n, steps=k) for a,n,k in CASES],
    input=previous.RELATIVE_OUTPUT+"/compatible.json", input_cache_nodes=513,
    profile_case=1, extracted_stage_rows=129, profile_training_rows=65, profile_withheld_rows=64,
    profile_selection="even saved rows fit; odd saved rows withheld",
    saved_time_nodes=129, exterior_quadrature_order=2, cap_quadrature_order=4,
    frozen_quadrature_case=4, frozen_quadrature_cuts=[69,79,128], frozen_quadrature_orders=[4,8],
    maximum_courant=.4, tube=TUBE, global_cpu_seconds=600., global_evaluations=100000000,
    source_continuation=previous.SETTINGS["source_continuation"],
    clock_origin=previous.SETTINGS["clock_origin"], lookup=previous.SETTINGS["lookup"],
    time_scheme="fixed-step SSPRK2", exterior_scheme="unchanged first-order vacuum-departure jets",
    cap_scheme="conservative departures about one reflected C2 characteristic profile",
    cap_restriction="sum adjacent total physical dr integrals Ip and Iq at common y cells",
    counts="actual shared Budget deltas; unsaved-stage point totals are measured, not independently recomputed")
SOURCES = tuple(dict.fromkeys((*previous.SOURCES,
    previous.RELATIVE_OUTPUT+"/compatible.json", previous.RELATIVE_OUTPUT+"/manifest.json",
    "scripts/diagnose_nsc_atlas_compatible_transport.py", "tests/test_diagnose_nsc_atlas_compatible_transport.py",
    "scripts/record_nsc_atlas_background_evolution.py", "tests/test_record_nsc_atlas_background_evolution.py",
    "tests/test_nsc_atlas_background_endpoint.py",
    *[f"src/recursive_horizons/{name}.py" for name in (
        "nsc_atlas_characteristic_background", "nsc_atlas_background_profile",
        "nsc_atlas_background_cap", "nsc_atlas_background_evolution")],
    *[f"tests/test_{name}.py" for name in (
        "nsc_atlas_characteristic_background", "nsc_atlas_background_profile",
        "nsc_atlas_background_cap", "nsc_atlas_background_evolution")],
)))
Budget = previous.Budget
plain, encoded, digest = previous.plain, previous.encoded, previous.digest
statistics, runtime_versions = previous.statistics, previous.runtime_versions
subset, StoredRows = previous.subset, previous.StoredRows
_equal, _array, _count, _real = previous._equal, previous._array, previous._count, previous._real


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def shared_budget(budget):
    if budget is None:
        return Budget(maximum_evaluations=100000000,cpu_seconds=600.)
    if (not isinstance(budget,Budget) or type(budget.maximum_evaluations) is not int
            or not 0<budget.maximum_evaluations<=100000000
            or isinstance(budget.cpu_seconds,bool) or not math.isfinite(budget.cpu_seconds)
            or not 0<budget.cpu_seconds<=600.
            or type(budget.evaluations) is not int or budget.evaluations<0
            or not math.isfinite(budget.started) or budget.started>time.process_time()):
        raise ValueError("invalid shared background allowance")
    budget.check()
    return budget


def load_parent(budget):
    directory = ROOT/previous.RELATIVE_OUTPUT
    if ({p.name for p in directory.iterdir()} != {"compatible.json", "manifest.json"}
            or any(p.is_symlink() for p in directory.iterdir())):
        raise ValueError("parent compatible inventory mismatch")
    raw, manifest_raw = (directory/"compatible.json").read_bytes(), (directory/"manifest.json").read_bytes()
    manifest = json.loads(manifest_raw)
    if (manifest.get("schema") != previous.SCHEMA+"-manifest"
            or not verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=previous.SOURCES)
            or manifest.get("outputs") != {"compatible.json":dict(sha256=digest(raw),bytes=len(raw))}):
        raise ValueError("parent compatible hash/source mismatch")
    parent = json.loads(raw)
    start = budget.evaluations
    authentication = previous.validate_record(parent, budget)
    base = json.loads((ROOT/previous.previous.RELATIVE_OUTPUT/"coupled.json").read_bytes())
    identity = dict(path=previous.RELATIVE_OUTPUT+"/compatible.json", sha256=digest(raw),
        manifest_sha256=digest(manifest_raw), authentication=authentication)
    return parent, base, identity, budget.evaluations-start


def extract_profile_inputs(parent, base, budget):
    case = parent["cases"][1]
    if (case["amplitude"],case["cells"],case["steps"],case["alpha"]) != ("zero",32,256,0.):
        raise ValueError("profile extraction requires old zero32/256")
    walls, rates = {k:[] for k in ("T","r","c")}, {k:[] for k in ("T","r","c")}
    start = budget.evaluations
    eta = np.linspace(0.,1.,33)
    for j,sample in enumerate(case["snapshots"]):
        stage = previous.compatible_stage(sample["state"], x=sample["x"], eta=eta,y_faces=eta,
            source=subset(base["initial_data"]["zero"]["snapshot_rows"],4*j),
            gauge=subset(base["gauge"]["rows"],4*j),alpha=0.,tube=TUBE,budget=budget)
        _equal(stage["physical_wall"],sample["physical"]["wall"],"profile same-stage wall")
        for key in walls:
            walls[key].append(stage["physical_wall"][key])
            rates[key].append(stage["wall_rates"][key])
    table = plain(dict(x=[s["x"] for s in case["snapshots"]],walls=walls,wall_rates=rates))
    return dict(table=table,table_sha256=digest(encoded(table)),
                stage_reconstructions=129,extraction_evaluations=budget.evaluations-start)


def fit_profile(inputs,budget):
    table=inputs["table"]
    if inputs["table_sha256"] != digest(encoded(table)):
        raise ValueError("retained profile table hash mismatch")
    walls={k:_array(table["walls"][k],(129,),"profile wall") for k in ("T","r","c")}
    rates={k:_array(table["wall_rates"][k],(129,),"profile rate") for k in ("T","r","c")}
    built=build_background_profile({k:v[::2] for k,v in walls.items()},
        {k:v[::2] for k,v in rates.items()},alpha=0.,budget=budget)
    labels=walls["T"]+np.arctanh(walls["r"])
    arrival_rate=rates["T"]+rates["r"]/(1-walls["r"]**2)
    if np.any(arrival_rate<=0) or np.any(np.diff(labels)<=0):
        raise ValueError("all retained arrival labels/rates must be positive and ordered")
    first=rates["c"]/arrival_rate
    direct=built["profile"].jet(labels[1::2],budget=budget)
    value_error=direct["F"]-walls["c"][1::2]
    slope_error=direct["F_prime"]-first[1::2]
    budget.check(3*129)
    summary={k:built[k] for k in ("labels","values","first","second","arrival_rate",
        "label_range","curvature_method","scope","diagnostics")}
    summary["withheld"]=dict(indices=np.arange(1,129,2),labels=labels[1::2],
        value_difference=value_error,first_difference=slope_error,
        value_statistics=statistics(value_error),first_statistics=statistics(slope_error),
        interpretation="withheld supplied values/first jets; second jets are computational interpolation curvature")
    return built["profile"],plain(summary)


def profile_replay_equal(actual, expected, inputs):
    """Check reconstructed profile arithmetic with term-conditioned radii."""
    # Small isolated recorder fixtures may replace the complete fitter.
    if "labels" not in expected:
        return _equal(actual, expected, "profile fit and withheld values")
    allowances = background_profile_roundoff_tolerances(inputs)

    def compare(a, b, path=""):
        label = "profile replay "+path
        if isinstance(b,dict):
            if not isinstance(a,dict) or set(a)!=set(b):
                raise ValueError(label+" inventory mismatch")
            for key in b:
                compare(a[key],b[key],path+"."+key if path else key)
        elif path in allowances:
            aa,bb = np.asarray(a),np.asarray(b)
            if (aa.dtype.kind not in "iuf" or aa.shape!=bb.shape
                    or any(isinstance(v,(bool,np.bool_)) for v in np.asarray(a,dtype=object).flat)
                    or not np.all(np.isfinite(aa)) or not np.all(np.isfinite(bb))):
                raise ValueError(label+" shape/type/nonfinite mismatch")
            allowance = np.broadcast_to(2e-13+2e-12*abs(bb)+allowances[path],bb.shape)
            difference = abs(aa-bb)
            if np.any(difference>allowance):
                worst = int(np.argmax(difference/allowance))
                raise ValueError(label+f" mismatch; max_abs_difference={float(np.max(difference)):.17g}; "
                    f"worst_index={worst}; difference={float(difference.flat[worst]):.17g}; "
                    f"allowance={float(allowance.flat[worst]):.17g}")
        elif path=="values" or path.startswith("diagnostics.quiet_anchor") or path=="withheld.indices":
            if encoded(a)!=encoded(b):
                raise ValueError(label+" exact copied value/anchor/index mismatch")
        else:
            _equal(a,b,label)
    compare(actual,expected)


def physical_view(state,x,source,gauge,eta,alpha,profile,budget):
    base=vacuum_slice(eta*float(gauge["W"]),x,budget=budget)
    wall={k:float(gauge[k])+float(state["wall"][k]) for k in ("r","L","T","s")}
    radial=wall["r"]*eta
    integrals=profile.integrals(wall["T"],radial,budget=budget)
    background=profile.background(wall["T"],np.array([wall["r"]]),budget=budget)
    chi=float(background["chi"][0]+np.cumsum(np.asarray(state["cap"]["Iq"]))[-1])
    wall["c"]=wall["r"]*_wall_conversion(chi/wall["r"],alpha)
    cap={k:np.asarray(state["cap"][k])+integrals[k] for k in ("Ip","Iq")}
    increment={}
    for key in ("r","L","c"):
        src=float(source["log_C" if key=="L" else key])
        vacuum=base[key] if key!="c" else np.zeros_like(eta)
        offset=src-vacuum[0]
        increment[key]=vacuum[1:-1]+offset+np.asarray(state["departure"][key])-src
    return plain(dict(increment=increment,wall=wall,cap=cap))


def final_summary(stage):
    metric,rhs=stage["cap_metric"],stage["cap_rhs"]
    return plain(dict(physical_wall=stage["physical_wall"],wall_rates=stage["wall_rates"],
        wall_kinematic_defect=stage["wall_kinematic_defect"],
        canonical_wall_flux_defect=stage["canonical_wall_flux_defect"],
        conformal_wall_flux_defect=stage["conformal_wall_flux_defect"],
        cap_center=metric["center"],cap_wall=metric["wall"],
        cap_faces={k:metric[k] for k in ("r_faces","m_faces","n_faces","A_faces","delta_faces","B_faces","chi_faces")},
        cap_gamma_interpolant=statistics(metric["gamma_interpolant_defect_nodes"]),
        cap_panels=len(metric["panel_cuts"])-1,
        cap_metric_evaluations=metric["charged_evaluations"],cap_rhs_evaluations=rhs["charged_evaluations"],
        cap_flux={k:rhs[k] for k in ("p_flux","q_flux","reaction_integral","wall_chi_x",
            "baseline_wall_chi_x","wall_f_x","wall_c_x","incoming_departure","constructed_wall_flux_defect")},
        interpretation="final radial/flux reconstruction; constructed primitive identities are not propagation accuracy tests"))


def tables(case,source_rows,budget):
    n=case["cells"]
    exterior={k:[] for k in ("r","L","c")}
    for j,s in enumerate(case["snapshots"]):
        for key in exterior:
            src=source_rows["log_C" if key=="L" else key][4*j]
            exterior[key].append(np.r_[src,src+np.asarray(s["physical"]["increment"][key]),s["physical"]["wall"][key]])
    budget.check(3*129*(n+1))
    return dict(exterior={k:np.asarray(v) for k,v in exterior.items()},
        wall={k:np.asarray([s["physical"]["wall"][k] for s in case["snapshots"]]) for k in ("r","L","c","T","s")},
        cap={k:np.asarray([s["physical"]["cap"][k] for s in case["snapshots"]]) for k in ("Ip","Iq")})


def cap_history(case,profile,budget):
    arrays={k:[] for k in ("outgoing_a_minus","departure_a_minus","center_f","maximum_primitive_departure","upper_label_margin")}
    for sample in case["snapshots"]:
        n,R,T=case["cells"],sample["physical"]["wall"]["r"],sample["physical"]["wall"]["T"]
        rec=_reconstruct(R*np.linspace(0,1,n+1),sample["state"]["cap"])
        bg=profile.background(T,np.array([0.,R]),budget=budget)
        dep=float(rec["minus"][-1]+rec["minus_slope"][-1]*(R-rec["midpoint"][-1]))
        arrays["departure_a_minus"].append(dep)
        arrays["outgoing_a_minus"].append(dep+float(bg["p"][-1]-bg["q"][-1]))
        arrays["center_f"].append(float(rec["Q0"]+bg["f"][0]))
        arrays["maximum_primitive_departure"].append(float(np.max(abs(rec["chi"]))))
        arrays["upper_label_margin"].append(float(profile.domain["last_label"]-(T+np.arctanh(R))))
        budget.check(4*n+1)
    return dict(arrays=arrays,statistics={k:statistics(v) for k,v in arrays.items()})


def frozen_quadrature(case,profile,budget):
    rows=[]
    for index in SETTINGS["frozen_quadrature_cuts"]:
        s=case["snapshots"][index]
        n,R,T=case["cells"],s["physical"]["wall"]["r"],s["physical"]["wall"]["T"]
        values={}
        for order in (4,8):
            m=background_cap_metric(np.linspace(0,1,n+1),R,s["state"]["cap"],
                T=T,profile=profile,alpha=case["alpha"],quadrature_order=order,budget=budget)
            values[str(order)]=dict(n_wall=float(m["n_faces"][-1]),delta_wall=m["wall"]["delta"],
                A_wall=m["wall"]["A"],w_minus=m["wall"]["w_minus"],
                max_abs_gamma_interpolant_defect=float(np.max(abs(m["gamma_interpolant_defect_nodes"]))),
                panels=len(m["panel_cuts"])-1,evaluations=m["charged_evaluations"])
        rows.append(dict(snapshot=index,x=s["x"],orders=values,
            difference={k:values["8"][k]-values["4"][k] for k in
                        ("n_wall","delta_wall","A_wall","w_minus")}))
    return dict(rows=rows,scope="same saved state/profile; local radial quadrature sensitivity, not a global error bound")


def diagnostics(record,parent,base,profile,budget):
    histories,table_list,caps,old_differences,oracles=[],[],[],[],[]
    old={(c["amplitude"],c["cells"],c["steps"]):c for c in parent["cases"]}
    for case in record["cases"]:
        # Existing diagnostics read physical increments and raw exterior Z only.
        proxy=dict(case)
        full,dep=previous._histories(base,proxy,budget)
        coarse,depcoarse=previous._histories(base,proxy,budget,2)
        histories.append(dict(raw129=previous.previous.compact_history(full),
            departure129=previous.previous.compact_history(dep),raw65=previous.previous.compact_history(coarse),
            departure65=previous.previous.compact_history(depcoarse)))
        table=tables(case,base["initial_data"][case["amplitude"]]["snapshot_rows"],budget)
        table_list.append(table)
        caps.append(cap_history(case,profile,budget))
        key=(case["amplitude"],case["cells"],case["steps"])
        if key in old:
            older=copy.deepcopy(old[key])
            for s in older["snapshots"]:
                s["physical"]["cap"]=s["state"]["cap"]
            old_differences.append(previous.previous._difference_summary(table,
                tables(older,base["initial_data"][case["amplitude"]]["snapshot_rows"],budget)))
        else:
            old_differences.append(None)
        if case["alpha"]==0:
            oracle=copy.deepcopy(case)
            for s in oracle["snapshots"]:
                s["state"]["cap"]=s["physical"]["cap"]
            oracles.append(cap_pullback(oracle,budget))
    difference=previous.previous._difference
    summarize=previous.previous._difference_summary
    paired=[difference(table_list[j+3],table_list[j]) for j in range(3)]
    return plain(dict(histories=histories,cap_histories=caps,zero_cap_pullback=oracles,
        successor_minus_compatible=old_differences,
        finite_minus_zero=[summarize(table_list[j+3],table_list[j]) for j in range(3)],
        spatial_zero=summarize(previous._restricted(table_list[1]),table_list[0]),
        spatial_finite=summarize(previous._restricted(table_list[4]),table_list[3]),
        paired_spatial=summarize(previous._restricted(paired[1]),paired[0]),
        time_zero=summarize(table_list[2],table_list[1]),time_finite=summarize(table_list[5],table_list[4]),
        paired_time=summarize(paired[2],paired[1]),
        frozen_cap_quadrature=frozen_quadrature(record["cases"][4],profile,budget),
        cap_restriction=SETTINGS["cap_restriction"],
        comparison_coordinates="common x/exterior eta/cap y; physical radii may differ",
        constructed_identities_are_accuracy_tests=False))


def calculate(budget=None):
    budget=shared_budget(budget)
    phase,record={"name":"inherited_authentication"},{}
    try:
        parent,base,identity,auth_charges=load_parent(budget)
        phase={"name":"extract_computational_profile"}
        inputs=extract_profile_inputs(parent,base,budget)
        start=budget.evaluations
        profile,profile_summary=fit_profile(inputs,budget)
        profile_charges=budget.evaluations-start
        x=np.asarray(base["gauge"]["rows"]["x"])
        record=dict(schema=SCHEMA,status=STATUS,scope=SCOPE,settings=SETTINGS,
            runtime=runtime_versions(),parent_input=identity,geometry=base["geometry"],
            profile_inputs=inputs,profile=profile_summary,cases=[],
            refinement_is_rigorous_error_enclosure=False)
        for label,n,steps in CASES:
            phase=dict(name="background_evolution",amplitude=label,cells=n,steps=steps)
            source=StoredRows(base["initial_data"][label]["snapshot_rows"])
            gauge=StoredRows(base["gauge"]["rows"])
            start=budget.evaluations
            solution=evolve_background(previous.quiet_departures(exterior_cells=n,cap_cells=n),
                x[::512//steps],eta=np.linspace(0,1,n+1),y_faces=np.linspace(0,1,n+1),
                source=source,gauge=gauge,alpha=base["initial_data"][label]["alpha"],tube=TUBE,
                budget=budget,profile=profile,save_every=steps//128,maximum_courant=.4,
                exterior_quadrature_order=2,cap_quadrature_order=4)
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
        record["execution"]=dict(evaluations=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
            inherited_authentication_evaluations=auth_charges,profile_extraction_evaluations=inputs["extraction_evaluations"],
            profile_build_evaluations=profile_charges,diagnostic_evaluations=budget.evaluations-start,
            old_profile_stage_reconstructions=129,coupled_evolutions=6,
            time_steps=sum(k for _,_,k in CASES),stage_evaluations=sum(2*k+1 for _,_,k in CASES),
            source_integrations=0,gauge_quadratures=0,new_physical_preparations=0,
            stage_point_totals_are_measured=True)
        budget.check()
        return plain(record)
    except Exception as error:
        error.background_failure=plain(dict(phase=phase,exception=type(error).__name__,message=str(error),
            stage=getattr(error,"diagnostic",None),evaluations=budget.evaluations,
            cpu_seconds=time.process_time()-budget.started,completed_cases=record.get("cases",[]),
            parent_input=record.get("parent_input")))
        raise


def validate_profile_inputs(inputs,parent,base,budget):
    table=inputs["table"]
    old=parent["cases"][1]
    _equal(table["x"],[s["x"] for s in old["snapshots"]],"profile times")
    for key in ("r","T","c"):
        _equal(table["walls"][key],[s["physical"]["wall"][key] for s in old["snapshots"]],"profile wall values")
    expected=129*(36*32+11)
    if inputs["stage_reconstructions"]!=129 or inputs["extraction_evaluations"]!=expected:
        raise ValueError("profile extraction accounting mismatch")
    # Alpha0 geometry is the unchanged reference; c rates are retained
    # producer-bound same-stage data, not differences of saved c samples.
    for key,ref in (("T","T_x"),("r","r_x_total")):
        np.testing.assert_allclose(table["wall_rates"][key],np.asarray(base["gauge"]["rows"][ref])[::4],
                                   rtol=2e-12,atol=2e-12)
    eta=np.linspace(0,1,33)
    s=old["snapshots"][-1]
    stage=previous.compatible_stage(s["state"],x=s["x"],eta=eta,y_faces=eta,
        source=subset(base["initial_data"]["zero"]["snapshot_rows"],-1),
        gauge=subset(base["gauge"]["rows"],-1),alpha=0.,tube=TUBE,budget=budget)
    for key in ("r","T","c"):
        _equal(table["wall_rates"][key][-1],stage["wall_rates"][key],"profile final first-jet spotcheck")


def validate_record(record,budget=None,*,require_creation=True):
    budget=shared_budget(budget)
    start=budget.evaluations
    try:
        if (record["schema"]!=SCHEMA or record["status"]!=STATUS or record["scope"]!=SCOPE
                or record["refinement_is_rigorous_error_enclosure"] is not False):
            raise ValueError("background declaration mismatch")
        _equal(record["settings"],SETTINGS,"settings")
        if (set(record["runtime"])!={"python","numpy","scipy"}
                or any(not isinstance(v,str) or not v for v in record["runtime"].values())):
            raise ValueError("runtime inventory mismatch")
        parent,base,identity,auth_charges=load_parent(budget)
        _equal(record["parent_input"],identity,"parent input")
        _equal(record["geometry"],base["geometry"],"geometry")
        validate_profile_inputs(record["profile_inputs"],parent,base,budget)
        before=budget.evaluations
        profile,profile_summary=fit_profile(record["profile_inputs"],budget)
        profile_charges=budget.evaluations-before
        profile_replay_equal(record["profile"],profile_summary,record["profile_inputs"])
        if len(record["cases"])!=6:
            raise ValueError("all six ordered cases required")
        x=np.asarray(base["gauge"]["rows"]["x"])
        end=auth_charges+record["profile_inputs"]["extraction_evaluations"]+profile_charges
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
            _equal(case["snapshots"][0]["state"],previous.quiet_departures(exterior_cells=n,cap_cells=n),"quiet departures")
            for key in ("T","s"):
                clocks=np.array([s["physical"]["wall"][key] for s in case["snapshots"]])
                if clocks[0]!=0 or np.any(np.diff(clocks)<=0):
                    raise ValueError("clocks must increase from quiet zero")
            before=budget.evaluations
            stage=background_stage(case["snapshots"][-1]["state"],x=x[-1],eta=eta,y_faces=eta,
                source=subset(base["initial_data"][label]["snapshot_rows"],-1),gauge=subset(base["gauge"]["rows"],-1),
                alpha=alpha,tube=TUBE,budget=budget,profile=profile)
            final_charge=budget.evaluations-before
            _equal(case["final"],final_summary(stage),"final cap/prefix/wall reconstruction")
            execution=case["execution"]
            charge=execution["measured_evaluations"]
            _count(charge,100000000,"measured case charges")
            if charge < (2*steps+1)*(8*(n+1)+2):
                raise ValueError("case charge below exterior/hybrid minimum")
            expected=dict(budget_start=end,budget_end=end+charge,measured_evaluations=charge,
                stages=2*steps+1,source_cache_reads=2*steps+1,gauge_cache_reads=2*steps+1,retries=0,
                final_stage_evaluations=final_charge)
            _equal({k:v for k,v in execution.items() if k!="max_courant"},expected,"measured case accounting")
            _real(execution["max_courant"],"Courant",high=.4)
            end+=charge
        before=budget.evaluations
        _equal(record["diagnostics"],diagnostics(record,parent,base,profile,budget),"saved diagnostics")
        diagnostic_charges=budget.evaluations-before
        expected=dict(evaluations=end+diagnostic_charges,inherited_authentication_evaluations=auth_charges,
            profile_extraction_evaluations=record["profile_inputs"]["extraction_evaluations"],
            profile_build_evaluations=profile_charges,diagnostic_evaluations=diagnostic_charges,
            old_profile_stage_reconstructions=129,coupled_evolutions=6,
            time_steps=sum(k for _,_,k in CASES),stage_evaluations=sum(2*k+1 for _,_,k in CASES),
            source_integrations=0,gauge_quadratures=0,new_physical_preparations=0,
            stage_point_totals_are_measured=True)
        _equal({k:v for k,v in record["execution"].items() if k!="cpu_seconds"},expected,"execution accounting")
        _count(record["execution"]["evaluations"],100000000,"calculation charges")
        _real(record["execution"]["cpu_seconds"],"calculation CPU",high=600.)
        budget.check()
        charges=budget.evaluations-start
        if require_creation:
            creation=record["creation_execution"]
            _equal({k:v for k,v in creation.items() if k!="cpu_seconds"},
                dict(calculation_evaluations=record["execution"]["evaluations"],validation_evaluations=charges,
                     total_evaluations=record["execution"]["evaluations"]+charges),"creation accounting")
            _count(creation["total_evaluations"],100000000,"creation charges")
            _real(creation["cpu_seconds"],"creation CPU",low=record["execution"]["cpu_seconds"],high=600.)
    except (KeyError,TypeError,IndexError,OverflowError) as error:
        raise ValueError("incomplete background record") from error
    return dict(authentication_evaluations=budget.evaluations-start,source_integrations=0,time_evolutions=0,
        reference_quadratures=0,old_profile_stage_spotchecks=1,profile_stage_extraction_replayed=False,
        current_final_radial_stages=6,frozen_radial_quadrature_slices=6,
        unsaved_stage_point_totals_independently_recomputed=False)


def create(directory):
    directory=Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only record: choose a new directory or --check")
    bindings,runtime,started=source_bindings(),runtime_versions(),time.process_time()
    budget=Budget(maximum_evaluations=100000000,cpu_seconds=600.,started=started)
    try:
        record=calculate(budget)
        if bindings!=source_bindings() or runtime!=runtime_versions() or record["runtime"]!=runtime:
            raise RuntimeError("background source/runtime changed during calculation")
        auth=validate_record(record,budget,require_creation=False)
        if bindings!=source_bindings() or runtime!=runtime_versions():
            raise RuntimeError("background source/runtime changed during authentication")
        record["creation_execution"]=dict(calculation_evaluations=record["execution"]["evaluations"],
            validation_evaluations=auth["authentication_evaluations"],total_evaluations=budget.evaluations,
            cpu_seconds=time.process_time()-started)
        budget.check()
    except Exception as error:
        receipt=dict(schema=SCHEMA+"-failed-attempt",status="FAILED_NOT_COMPLETED",settings=SETTINGS,
            source_bindings_before=bindings,runtime=runtime,
            failure=getattr(error,"background_failure",dict(exception=type(error).__name__,message=str(error),
                evaluations=budget.evaluations,cpu_seconds=time.process_time()-started)))
        path=ROOT/".build/background-evolution-failures"/(uuid.uuid4().hex+".json")
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open("xb") as stream:
            stream.write(encoded(previous.previous.previous._failure_safe(receipt)))
        error.failure_receipt=str(path)
        raise
    raw=encoded(record)
    manifest=dict(schema=SCHEMA+"-manifest",source_bindings=bindings,
                  outputs={"background.json":dict(sha256=digest(raw),bytes=len(raw))})
    directory.mkdir(parents=True,exist_ok=False)
    (directory/"background.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return dict(created=str(directory),execution=record["execution"],creation_execution=record["creation_execution"])


def check(directory):
    directory=Path(directory)
    if ({p.name for p in directory.iterdir()}!={"background.json","manifest.json"}
            or any(p.is_symlink() for p in directory.iterdir())):
        raise ValueError("background output inventory mismatch")
    manifest=json.loads((directory/"manifest.json").read_bytes())
    raw=(directory/"background.json").read_bytes()
    if (manifest.get("schema")!=SCHEMA+"-manifest" or not verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=SOURCES)
            or manifest.get("outputs")!={"background.json":dict(sha256=digest(raw),bytes=len(raw))}):
        raise ValueError("background hash/source mismatch")
    return dict(checked=str(directory),**validate_record(json.loads(raw)))


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
