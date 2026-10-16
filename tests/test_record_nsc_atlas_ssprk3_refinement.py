"""Manufactured refinement records; no original source/cache/trajectory run."""
import copy
import json
import time

import numpy as np
import pytest

import record_nsc_atlas_ssprk3_refinement as tool
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile

ORIGINAL_DIAGNOSTICS = tool.diagnostics
ORIGINAL_LOAD = tool.load_inputs


@pytest.fixture(scope="module")
def budget():
    value=tool.Budget(maximum_evaluations=1000000,cpu_seconds=30.)
    yield value
    value.check()
    print("refinement recorder controls",value.evaluations,time.process_time()-value.started)


@pytest.fixture(scope="module")
def synthetic(budget):
    patch=pytest.MonkeyPatch()
    cases=(("zero",3,512),("finite",3,512))
    patch.setattr(tool,"CASES",cases)
    oldx=np.linspace(.4,.41,513)
    x=np.empty(1025); x[::2]=oldx; x[1::2]=oldx[:-1]+(oldx[1:]-oldx[:-1])/2
    def gauge(at):
        return dict(x=at,W=np.full(at.size,.2),r=np.full(at.size,.45),L=np.zeros(at.size),T=at-.4,s=at-.4)
    def source(at):
        return dict(x=at,r=at,log_C=np.zeros(at.size),c=np.zeros(at.size),r_x=np.ones(at.size),
                    log_C_x=np.zeros(at.size),c_x=np.zeros(at.size))
    base=dict(geometry=dict(support_x=[.401,.405]),gauge=dict(rows=gauge(oldx)),
        initial_data={a:dict(alpha=alpha,snapshot_rows=source(oldx)) for a,alpha in (("zero",0.),("finite",.01))})
    cache=tool.plain(dict(schema="manufactured",x_nodes=x,source_rows={a:source(x) for a in ("zero","finite")},
        gauge_rows=gauge(x),diagnostics=dict(manufactured=True),execution=dict(evaluations=11)))
    parent=dict(profile_inputs=dict(saved=True),profile=dict(manufactured=True))
    identity=dict(manufactured=True)
    loads,stages,cache_calls=[],[],[]
    def load(shared):
        loads.append(shared); shared.check()
        return parent,{},base,identity
    patch.setattr(tool,"load_inputs",load)
    def make(base,shared):
        cache_calls.append(("make",shared)); shared.check(11)
        return copy.deepcopy(cache)
    def validate(value,base,shared):
        cache_calls.append(("check",shared)); shared.check(11)
        tool._equal(value,cache,"manufactured cache")
        return dict(validation_evaluations=11)
    patch.setattr(tool,"make_refined_cache",make)
    patch.setattr(tool,"validate_refined_cache",validate)
    profile=CharacteristicProfile([0.,1.5],[0.,0.],[0.,0.],[0.,0.])
    def fit(inputs,shared):
        shared.check(3)
        return profile,dict(manufactured=True)
    patch.setattr(tool,"fit_profile",fit)
    def stage(state,*,eta,budget,**kwargs):
        stages.append(budget); budget.check(8*len(eta)+5)
        return dict(state=state)
    patch.setattr(tool,"background_stage",stage)
    patch.setattr(tool,"final_summary",lambda row:dict(cap_metric_evaluations=1,cap_rhs_evaluations=2,
        cap=tool.plain(row["state"]["cap"])))
    def no_evolution(*args,**kwargs):
        raise AssertionError("authentication must not replay a trajectory")
    patch.setattr(tool,"evolve_background_ssprk3",no_evolution)
    def diagnostic(record,parent,background,base,profile,shared):
        shared.check(7)
        return dict(cap=[c["snapshots"][-1]["physical"]["cap"] for c in record["cases"]])
    patch.setattr(tool,"diagnostics",diagnostic)
    rows=[]; saved_end=14
    for label,n,steps in cases:
        alpha=base["initial_data"][label]["alpha"]; snapshots=[]
        for j,at in enumerate(x[::8]):
            state=tool.quiet_departures(exterior_cells=n,cap_cells=n)
            state["cap"]["Iq"][:]=alpha*(at-.4)**2
            state["cap"]["Ip"][:]=2*alpha*(at-.4)**2
            physical=tool.physical_view(state,at,tool.subset(cache["source_rows"][label],8*j),
                tool.subset(cache["gauge_rows"],8*j),np.linspace(0,1,n+1),alpha,profile,budget)
            snapshots.append(dict(x=float(at),state=tool.plain(state),physical=physical))
        count=3*steps+1; charge=count*(8*(n+1)+5)
        rows.append(dict(amplitude=label,cells=n,steps=steps,alpha=alpha,snapshots=snapshots,
            final=tool.final_summary(dict(state=state)),execution=dict(budget_start=saved_end,
                budget_end=saved_end+charge,measured_evaluations=charge,stages=count,source_cache_reads=count,
                gauge_cache_reads=count,retries=0,max_courant=.1,final_stage_evaluations=8*(n+1)+5)))
        saved_end+=charge
    # Measured trajectory counters here are explicit fixture data, never
    # fictitious charges against the live manufactured-control allowance.
    record=dict(schema=tool.SCHEMA,status=tool.STATUS,scope=tool.SCOPE,settings=tool.SETTINGS,
        runtime=tool.runtime_versions(),parent_input=identity,geometry=base["geometry"],cache=cache,
        profile_inputs=parent["profile_inputs"],profile=parent["profile"],cases=rows,uniform_error_enclosure=False,
        execution=dict(initial_evaluations=0,evaluations=saved_end+7,cpu_seconds=1.,
            inherited_authentication_evaluations=0,cache_evaluations=11,profile_build_evaluations=3,
            diagnostic_evaluations=7,coupled_evolutions=2,time_steps=1024,stage_evaluations=3074,
            source_integrations=0,old_trajectory_replays=0,profile_stage_extractions=0,stage_point_totals_are_measured=True))
    record["diagnostics"]=diagnostic(record,parent,{},base,profile,budget)
    yield dict(record=record,base=base,parent=parent,profile=profile,cache=cache,
               stages=stages,loads=loads,cache_calls=cache_calls)
    patch.undo()


def test_saved_reconstruction_cache_and_two_final_stages_share_budget(synthetic,budget):
    before=budget.evaluations; stage_count=len(synthetic["stages"])
    result=tool.validate_record(copy.deepcopy(synthetic["record"]),budget,require_creation=False)
    assert result["authentication_evaluations"]==budget.evaluations-before
    assert len(synthetic["stages"])-stage_count==2
    assert result["source_integrations"]==result["time_evolutions"]==0
    assert all(b is budget for b in synthetic["loads"]+synthetic["stages"])
    assert all(b is budget for _,b in synthetic["cache_calls"])


def test_calculation_uses_original_nodes_and_only_new_midpoints(synthetic,budget,monkeypatch):
    calls=[]
    def evolve(initial,nodes,*,mid_nodes,source,gauge,budget,save_every,**kwargs):
        saved=synthetic["record"]["cases"][len(calls)]
        calls.append((nodes,mid_nodes,save_every,budget))
        source.calls=gauge.calls=1537; budget.check(1)
        return dict(snapshots=saved["snapshots"],final_stage=dict(state=saved["snapshots"][-1]["state"]),
                    stages=1537,steps=512,max_courant=.1,retries=0,time_scheme=tool.previous.TIME_SCHEME)
    monkeypatch.setattr(tool,"evolve_background_ssprk3",evolve)
    record=tool.calculate(budget)
    assert record["execution"]["time_steps"]==1024 and record["execution"]["stage_evaluations"]==3074
    for nodes,mids,stride,shared in calls:
        np.testing.assert_array_equal(nodes,synthetic["base"]["gauge"]["rows"]["x"])
        np.testing.assert_array_equal(mids,np.asarray(synthetic["cache"]["x_nodes"])[1::2])
        assert stride==4 and shared is budget


@pytest.mark.parametrize("mutation",("cache","raw","physical","shape","time","identity","profile","runtime","final","counts","diagnostics","total"))
def test_record_authentication_rejects_changes(synthetic,budget,mutation):
    value=copy.deepcopy(synthetic["record"])
    sample=value["cases"][0]["snapshots"][1]
    if mutation=="cache": value["cache"]["source_rows"]["finite"]["c"][1]+=.01
    elif mutation=="raw": sample["state"]["cap"]["Iq"][0]+=.01
    elif mutation=="physical": sample["physical"]["wall"]["c"]+=.01
    elif mutation=="shape": sample["state"]["departure"]["r"].append(0.)
    elif mutation=="time": sample["x"]=float(np.nextafter(sample["x"],1.))
    elif mutation=="identity": value["cases"][0]["steps"]=256
    elif mutation=="profile": value["profile_inputs"]["saved"]=False
    elif mutation=="runtime": value["runtime"]["numpy"]=""
    elif mutation=="final": value["cases"][0]["final"]["cap_metric_evaluations"]+=1
    elif mutation=="counts": value["cases"][0]["execution"]["stages"]+=1
    elif mutation=="diagnostics": value["diagnostics"]["cap"][0]["Ip"][0]+=.01
    else: value["execution"]["cache_evaluations"]+=1
    with pytest.raises(ValueError): tool.validate_record(value,budget,require_creation=False)


def test_full_comparison_arrays_maximum_and_profile_brackets(budget):
    x=np.array([.4,.5,.6])
    result=tool.comparison(dict(wall=dict(c=np.array([1.,-3.,2.]))),dict(wall=dict(c=np.zeros(3))),x,budget)
    assert result["arrays"]["wall"]["c"]==[1.,-3.,2.]
    assert result["max_locations"]["wall"]["c"]==dict(index=[1],x=.5,signed=-3.)
    profile=CharacteristicProfile([.1,.3,.5,.7],np.zeros(4),np.zeros(4),np.zeros(4))
    q=[.1,.4,.7]; case=dict(snapshots=[dict(x=at,physical=dict(wall=dict(r=0.,T=label))) for at,label in zip(x,q)])
    brackets=tool.profile_intervals(case,profile,budget)
    assert brackets["strictly_interior_knots_per_saved_interval"]==[1,1]
    assert brackets["left_knot_index"]==[0,1,3] and brackets["right_knot_index"]==[0,2,3]


def test_diagnostics_pairs_all_histories_and_keeps_total_cap_oracle(synthetic,budget,monkeypatch):
    value=copy.deepcopy(synthetic["record"])
    for case,marker in zip(value["cases"],(3.,11.)): case["marker"]=marker
    # Deliberately distinguish total physical and raw cap data in the oracle.
    for sample in value["cases"][0]["snapshots"]: sample["physical"]["cap"]["Ip"]=[7.,7.,7.]
    def case(a,k,v):
        result=copy.deepcopy(value["cases"][0 if a=="zero" else 1])
        result.update(amplitude=a,cells=32,steps=k,marker=v)
        return result
    parent=dict(cases=[case(a,k,v) for a,k,v in (("zero",128,1.),("zero",256,2.),("finite",128,4.),("finite",256,8.))])
    background=dict(cases=[case(a,k,v) for a,k,v in (("zero",256,10.),("zero",512,20.),("finite",256,40.),("finite",512,80.))])
    def tables(row,source,shared):
        shared.check(1); a=np.full(129,row["marker"])
        return dict(wall=dict(c=a),exterior=dict(c=a),cap=dict(Ip=a,Iq=a))
    monkeypatch.setattr(tool,"tables",tables)
    monkeypatch.setattr(tool,"_histories",lambda *args: ({},{}))
    monkeypatch.setattr(tool,"compact_history",lambda a:a)
    monkeypatch.setattr(tool,"cap_history",lambda *args:dict(arrays=dict(upper_label_margin=[1.])))
    def oracle(case,shared):
        assert case["snapshots"][0]["state"]["cap"]["Ip"]==[7.,7.,7.]
        return dict(checked=True)
    monkeypatch.setattr(tool,"cap_pullback",oracle)
    result=ORIGINAL_DIAGNOSTICS(value,parent,background,synthetic["base"],synthetic["profile"],budget)
    np.testing.assert_array_equal(result["paired_time"]["arrays"]["wall"]["c"],np.full(129,2.))
    np.testing.assert_array_equal(result["old_rk2_paired_extrapolant"]["wall"]["c"],np.full(129,70.))
    np.testing.assert_array_equal(result["paired_versus_old_rk2_extrapolant"]["arrays"]["wall"]["c"],np.full(129,-62.))
    assert result["zero_cap_pullback"]["checked"] and len(result["profile_intervals"])==6
    assert result["temporal_order_established"] is False


def test_pinned_input_manifest_and_bindings_without_old_validator(tmp_path,budget,monkeypatch):
    monkeypatch.setattr(tool,"ROOT",tmp_path)
    monkeypatch.setattr(tool.previous,"RELATIVE_OUTPUT","old")
    names=tuple("bound"+str(i)+".txt" for i in range(80))
    monkeypatch.setattr(tool.previous,"SOURCES",names)
    bindings={name:tool.digest(name.encode()) for name in names}
    for name in names: (tmp_path/name).write_bytes(name.encode())
    identity=dict(frozen=True)
    monkeypatch.setattr(tool.previous,"load_inputs",lambda shared: ({},{},identity))
    def no_validator(*args,**kwargs): raise AssertionError("old trajectory validation prohibited")
    monkeypatch.setattr(tool.previous,"validate_record",no_validator)
    payload=dict(schema=tool.previous.SCHEMA,settings=tool.previous.SETTINGS,runtime=tool.runtime_versions(),parent_input=identity)
    raw=tool.encoded(payload)
    manifest=tool.encoded(dict(schema=tool.previous.SCHEMA+"-manifest",source_bindings=bindings,
        outputs={"ssprk3.json":dict(sha256=tool.digest(raw),bytes=len(raw))}))
    monkeypatch.setattr(tool,"INPUT_SHA256",tool.digest(raw)); monkeypatch.setattr(tool,"MANIFEST_SHA256",tool.digest(manifest))
    directory=tmp_path/"old"; directory.mkdir()
    (directory/"ssprk3.json").write_bytes(raw); (directory/"manifest.json").write_bytes(manifest)
    assert ORIGINAL_LOAD(budget)[3]["authenticated_input_files"]==80
    (tmp_path/names[0]).write_bytes(b"changed")
    with pytest.raises(ValueError,match="binding"): ORIGINAL_LOAD(budget)
    (tmp_path/names[0]).write_bytes(names[0].encode())
    (directory/"manifest.json").write_bytes(manifest+b" ")
    with pytest.raises(ValueError,match="pinned"): ORIGINAL_LOAD(budget)


def test_creation_check_shared_clocks_hashes_and_failure_receipt(tmp_path,monkeypatch):
    monkeypatch.setattr(tool,"ROOT",tmp_path)
    created=[]; original=tool.Budget
    def make(**kwargs):
        value=original(**kwargs); created.append(value); return value
    monkeypatch.setattr(tool,"Budget",make)
    changed=dict(value=False)
    def bindings():
        assert created and created[-1].started<=time.process_time()
        return dict(source="changed" if changed["value"] else "fixed")
    monkeypatch.setattr(tool,"source_bindings",bindings)
    calls=[]
    def calculate(shared):
        calls.append(shared); shared.check(20)
        return dict(runtime=tool.runtime_versions(),execution=dict(evaluations=20,cpu_seconds=0.))
    def validate(record,shared,**kwargs):
        calls.append(shared); shared.check(10)
        return dict(authentication_evaluations=10)
    monkeypatch.setattr(tool,"calculate",calculate); monkeypatch.setattr(tool,"validate_record",validate)
    output=tmp_path/"new"
    assert tool.create(output)["creation_execution"]["total_evaluations"]==30
    assert calls[0] is calls[1]
    with pytest.raises(FileExistsError): tool.create(output)
    tool.check(output)
    assert calls[-1].cpu_seconds==60 and calls[-1].maximum_evaluations==5000000
    def drift(record,shared,**kwargs):
        result=validate(record,shared,**kwargs); changed["value"]=True; return result
    monkeypatch.setattr(tool,"validate_record",drift)
    with pytest.raises(RuntimeError,match="authentication"): tool.check(output)
    changed["value"]=False
    rejected=tmp_path/"rejected"
    with pytest.raises(RuntimeError,match="authentication") as error: tool.create(rejected)
    assert not rejected.exists()
    assert json.loads(open(error.value.failure_receipt).read())["status"]=="FAILED_NOT_COMPLETED"
    changed["value"]=False
    (output/"refinement.json").write_bytes(b"{}")
    with pytest.raises(ValueError,match="hash/source"): tool.check(output)


def test_shared_allowance_preserves_accumulation_and_rejects_expansion():
    shared=tool.Budget(maximum_evaluations=100,cpu_seconds=30.)
    shared.check(90); started=shared.started
    assert tool.shared_budget(shared) is shared and shared.evaluations==90 and shared.started==started
    shared.maximum_evaluations=50000001
    with pytest.raises(ValueError): tool.shared_budget(shared)
    shared.maximum_evaluations=89
    with pytest.raises(RuntimeError): tool.shared_budget(shared)
