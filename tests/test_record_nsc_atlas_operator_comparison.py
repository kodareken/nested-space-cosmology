"""Manufactured recorder controls; no original source, cache or time solve."""
import copy
import json
import time

import numpy as np
import pytest

import record_nsc_atlas_operator_comparison as tool
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile


@pytest.fixture(scope="module")
def budget():
    shared = tool.Budget(maximum_evaluations=1000000,cpu_seconds=30.)
    yield shared
    shared.check()
    print("operator recorder controls",shared.evaluations,time.process_time()-shared.started)


@pytest.fixture(scope="module")
def synthetic(budget):
    patch=pytest.MonkeyPatch()
    cases=(("zero",4,512),("finite",4,512))
    patch.setattr(tool,"CASES",cases)
    oldx=np.linspace(.4,.41,513)
    x=np.empty(1025);x[::2]=oldx;x[1::2]=oldx[:-1]+(oldx[1:]-oldx[:-1])/2
    def gauge(at):
        return dict(x=at,W=np.full(at.size,.2),r=np.full(at.size,.45),L=np.zeros(at.size),T=at-.4,s=at-.4)
    def source(at):
        return dict(x=at,r=at,log_C=np.zeros(at.size),c=np.zeros(at.size),r_x=np.ones(at.size),
                    log_C_x=np.zeros(at.size),c_x=np.zeros(at.size))
    base=dict(geometry=dict(support_x=[.401,.405]),gauge=dict(rows=gauge(oldx)),
        initial_data={a:dict(alpha=alpha,snapshot_rows=source(oldx)) for a,alpha in (("zero",0.),("finite",.01))})
    cache=tool.plain(dict(x_nodes=x,source_rows={a:source(x) for a in ("zero","finite")},gauge_rows=gauge(x)))
    parent=dict(profile_inputs=dict(saved=True),profile=dict(manufactured=True),cache=cache,geometry=base["geometry"])
    identity=dict(manufactured=True)
    loads,stages=[],[]
    def load(shared):
        loads.append(shared);shared.check()
        return parent,parent,base,identity
    patch.setattr(tool,"load_inputs",load)
    profile=CharacteristicProfile([0.,1.5],[0.,0.],[0.,0.],[0.,0.])
    def fit(inputs,shared):
        shared.check(3);return profile,dict(manufactured=True)
    patch.setattr(tool,"fit_profile",fit)
    def stage(state,*,eta,budget,**kwargs):
        stages.append(budget);budget.check(8*len(eta)+5);return dict(state=state,charged_evaluations=8*len(eta)+5)
    patch.setattr(tool,"coupled_stage",stage)
    patch.setattr(tool,"final_summary",lambda row:dict(cap_metric_evaluations=1,cap_rhs_evaluations=2,
        cap=tool.plain(row["state"]["cap"])))
    def forbidden(*args,**kwargs):raise AssertionError("old solve/cache/validator must not run")
    patch.setattr(tool,"evolve_stage_ssprk3",forbidden)
    patch.setattr(tool.previous,"validate_record",forbidden)
    patch.setattr(tool.previous.previous,"make_refined_cache",forbidden)
    patch.setattr(tool.previous.previous,"validate_refined_cache",forbidden)
    def diagnostic(record,parent,base,cache,profile,shared):
        shared.check(7)
        return dict(cap=[c["snapshots"][-1]["physical"]["cap"] for c in record["cases"]])
    patch.setattr(tool,"diagnostics",diagnostic)
    rows=[];saved_end=3
    for label,n,steps in cases:
        alpha=base["initial_data"][label]["alpha"];snapshots=[]
        for j,at in enumerate(x[::8]):
            state=tool.quiet_departures(exterior_cells=n,cap_cells=n)
            state["cap"]["Iq"][:]=alpha*(at-.4)**2
            state["cap"]["Ip"][:]=2*alpha*(at-.4)**2
            physical=tool.physical_view(state,at,tool.subset(cache["source_rows"][label],8*j),
                tool.subset(cache["gauge_rows"],8*j),np.linspace(0,1,n+1),alpha,profile,budget)
            snapshots.append(dict(x=float(at),state=tool.plain(state),physical=physical))
        count=3*steps+1;charge=count*(8*(n+1)+5)
        rows.append(dict(amplitude=label,cells=n,steps=steps,alpha=alpha,snapshots=snapshots,
            final=tool.final_summary(dict(state=state)),execution=dict(budget_start=saved_end,budget_end=saved_end+charge,
                measured_evaluations=charge,stages=count,source_cache_reads=count,gauge_cache_reads=count,
                retries=0,max_courant=.1,final_stage_evaluations=8*(n+1)+5)))
        saved_end+=charge
    # Unsaved trajectory counters are explicit fixture metadata, not invented
    # live work: all actual physical-view/authentication arithmetic is charged.
    record=dict(schema=tool.SCHEMA,status=tool.STATUS,scope=tool.SCOPE,settings=tool.SETTINGS,runtime=tool.runtime_versions(),
        parent_input=identity,cache_identity=tool.cache_identity(parent),geometry=base["geometry"],
        profile_inputs=parent["profile_inputs"],profile=parent["profile"],cases=rows,uniform_error_enclosure=False,
        execution=dict(initial_evaluations=0,evaluations=saved_end+7,cpu_seconds=1.,
            inherited_authentication_evaluations=0,profile_build_evaluations=3,diagnostic_evaluations=7,
            coupled_evolutions=2,time_steps=1024,stage_evaluations=3074,source_integrations=0,
            cache_reconstructions=0,gauge_quadratures=0,old_trajectory_replays=0,profile_stage_extractions=0,
            stage_point_totals_are_measured=True))
    record["diagnostics"]=diagnostic(record,parent,base,cache,profile,budget)
    yield dict(record=record,parent=parent,base=base,cache=cache,profile=profile,loads=loads,stages=stages)
    patch.undo()


def test_full_saved_reconstruction_and_two_final_stages_without_replay(synthetic,budget):
    before=budget.evaluations;stages=len(synthetic["stages"]);clock=budget.started
    result=tool.validate_record(copy.deepcopy(synthetic["record"]),budget,require_creation=False)
    assert result["authentication_evaluations"]==budget.evaluations-before and budget.started==clock
    assert len(synthetic["stages"])-stages==2
    assert result["source_integrations"]==result["cache_reconstructions"]==result["gauge_quadratures"]==result["time_evolutions"]==0
    assert all(b is budget for b in synthetic["loads"]+synthetic["stages"])
    sample=synthetic["record"]["cases"][1]["snapshots"][-1]
    chi=sum(sample["state"]["cap"]["Iq"]);r=sample["physical"]["wall"]["r"];alpha=.01
    expected=r*np.sinh(np.sqrt(alpha/6)*chi/r)/np.sqrt(alpha/6)
    np.testing.assert_allclose(sample["physical"]["wall"]["c"],expected,rtol=2e-15)


def test_new_cases_preserve_exact512step_cache_schedule(synthetic,budget,monkeypatch):
    calls=[]
    monkeypatch.setattr(tool,"CASES",(("zero",64,512),("finite",64,512)))
    def evolve(initial,nodes,*,mid_nodes,stage,source,gauge,budget,save_every,eta,y_faces,**kwargs):
        assert stage is tool.coupled_stage
        calls.append((initial,nodes.copy(),mid_nodes.copy(),save_every,eta.copy(),y_faces.copy(),budget))
        source.calls=gauge.calls=1537;budget.check(1)
        saved=synthetic["record"]["cases"][len(calls)-1]
        # This stub only checks dispatch; it is not claimed as a valid64-cell result.
        return dict(snapshots=saved["snapshots"],final_stage=dict(state=saved["snapshots"][-1]["state"],charged_evaluations=123),
            stages=1537,steps=512,max_courant=.1,retries=0,time_scheme=tool.SETTINGS["time_scheme"])
    monkeypatch.setattr(tool,"evolve_stage_ssprk3",evolve)
    row=tool.calculate(budget)
    assert row["cases"][0]["execution"]["final_stage_evaluations"]==123
    assert row["execution"]["time_steps"]==1024 and row["execution"]["stage_evaluations"]==3074
    assert "cache" not in row and row["cache_identity"]==tool.cache_identity(synthetic["parent"])
    for initial,nodes,mids,stride,eta,y,shared in calls:
        np.testing.assert_array_equal(nodes,np.asarray(synthetic["cache"]["x_nodes"])[::2])
        np.testing.assert_array_equal(mids,np.asarray(synthetic["cache"]["x_nodes"])[1::2])
        np.testing.assert_array_equal(eta,np.linspace(0,1,65));np.testing.assert_array_equal(y,eta)
        assert len(initial["departure"]["r"])==63 and len(initial["cap"]["Iq"])==64
        assert stride==4 and shared is budget


@pytest.mark.parametrize("mutation",("cache_identity","raw","physical","shape","time","case","profile","runtime","final","counts","diagnostics","total"))
def test_saved_output_mutations_rejected(synthetic,budget,mutation):
    row=copy.deepcopy(synthetic["record"]);sample=row["cases"][0]["snapshots"][1]
    if mutation=="cache_identity":row["cache_identity"]["cache_sha256"]="changed"
    elif mutation=="raw":sample["state"]["cap"]["Iq"][0]+=.01
    elif mutation=="physical":sample["physical"]["wall"]["c"]+=.01
    elif mutation=="shape":sample["state"]["departure"]["r"].append(0.)
    elif mutation=="time":sample["x"]=float(np.nextafter(sample["x"],1.))
    elif mutation=="case":row["cases"][0]["cells"]=64
    elif mutation=="profile":row["profile_inputs"]["saved"]=False
    elif mutation=="runtime":row["runtime"]["numpy"]=""
    elif mutation=="final":row["cases"][0]["final"]["cap_metric_evaluations"]+=1
    elif mutation=="counts":row["cases"][0]["execution"]["stages"]+=1
    elif mutation=="diagnostics":row["diagnostics"]["cap"][0]["Iq"][0]+=.01
    else:row["execution"]["gauge_quadratures"]=1
    with pytest.raises(ValueError):tool.validate_record(row,budget,require_creation=False)


def test_parent_cache_change_cannot_hide_behind_reference(synthetic,budget):
    original=synthetic["parent"]["cache"]["source_rows"]["finite"]["c"][1]
    synthetic["parent"]["cache"]["source_rows"]["finite"]["c"][1]=.01
    try:
        with pytest.raises(ValueError,match="cache identity"):
            tool.validate_record(copy.deepcopy(synthetic["record"]),budget,require_creation=False)
    finally:synthetic["parent"]["cache"]["source_rows"]["finite"]["c"][1]=original


def test_pinned_parent_and_all92_paths_without_validator(tmp_path,budget,monkeypatch):
    monkeypatch.setattr(tool,"ROOT",tmp_path)
    names=tuple(["bound"+str(i)+".txt" for i in range(90)]+[tool.BASE_INPUT,tool.CACHE_INPUT])
    monkeypatch.setattr(tool.previous,"SOURCES",names)
    base=dict(geometry=dict(fixed=True))
    ancestor=dict(geometry=base["geometry"],profile_inputs=dict(fixed=True),cache=dict(marker=1025))
    contents={name:name.encode() for name in names}
    contents[tool.BASE_INPUT]=tool.encoded(base); contents[tool.CACHE_INPUT]=tool.encoded(ancestor)
    for name,data in contents.items():
        path=tmp_path/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
    parent=dict(schema=tool.previous.SCHEMA,settings=tool.previous.SETTINGS,runtime=tool.runtime_versions(),
        geometry=base["geometry"],profile_inputs=ancestor["profile_inputs"],cache_identity=tool.cache_identity(ancestor),
        cases=[dict(amplitude=a,cells=64,steps=512) for a in ("zero","finite")])
    raw=tool.encoded(parent)
    manifest=dict(schema=tool.previous.SCHEMA+"-manifest",source_bindings={k:tool.digest(v) for k,v in contents.items()},
        outputs={"spatial.json":dict(sha256=tool.digest(raw),bytes=len(raw))})
    mr=tool.encoded(manifest)
    monkeypatch.setattr(tool,"INPUT_SHA256",tool.digest(raw));monkeypatch.setattr(tool,"MANIFEST_SHA256",tool.digest(mr))
    path=tmp_path/tool.INPUT;path.parent.mkdir(parents=True);path.write_bytes(raw);path.with_name("manifest.json").write_bytes(mr)
    def forbidden(*args,**kwargs):raise AssertionError("inherited validator prohibited")
    monkeypatch.setattr(tool.previous,"validate_record",forbidden)
    got=ORIGINAL_LOAD(budget)
    assert got[3]["authenticated_input_files"]==92 and got[2]==base and got[1]==ancestor
    (tmp_path/names[0]).write_text("changed")
    with pytest.raises(ValueError,match="source.*chang|binding"):ORIGINAL_LOAD(budget)
    manifest["source_bindings"][names[0]]=tool.digest(b"changed");path.with_name("manifest.json").write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError,match="pinned"):ORIGINAL_LOAD(budget)


def test_creation_and_check_share_live_allowance_no_overwrite_and_auth_drift(tmp_path,budget,monkeypatch):
    monkeypatch.setattr(tool,"ROOT",tmp_path)
    monkeypatch.setattr(tool,"Budget",lambda **kwargs:budget)
    bindings=dict(frozen="one")
    monkeypatch.setattr(tool,"source_bindings",lambda:dict(bindings))
    monkeypatch.setattr(tool, "verify_saved_bindings", lambda expected, root, **kw: expected == tool.source_bindings())
    calls=[]
    def calculate(shared):
        assert shared is budget;shared.check(2)
        return dict(runtime=tool.runtime_versions(),value=[1.,-2.],execution=dict(evaluations=shared.evaluations,cpu_seconds=0.))
    def validate(row,shared,**kwargs):
        assert shared is budget;shared.check(3);calls.append(shared)
        if row["value"]!=[1.,-2.]:raise ValueError("manufactured validation mismatch")
        return dict(authentication_evaluations=3,time_evolutions=0)
    monkeypatch.setattr(tool,"calculate",calculate);monkeypatch.setattr(tool,"validate_record",validate)
    path=tmp_path/"result";before=budget.evaluations
    result=tool.create(path)
    assert budget.evaluations-before==5 and result["creation_execution"]["total_evaluations"]==budget.evaluations
    assert tool.check(path)["time_evolutions"]==0
    with pytest.raises(FileExistsError):tool.create(path)
    row=json.loads((path/"operator.json").read_text());row["value"][1]=2.
    raw=tool.encoded(row);(path/"operator.json").write_bytes(raw)
    manifest=json.loads((path/"manifest.json").read_text());manifest["outputs"]["operator.json"]=dict(sha256=tool.digest(raw),bytes=len(raw))
    (path/"manifest.json").write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError,match="validation mismatch"):tool.check(path)
    def drift(row,shared,**kwargs):
        shared.check(1);bindings["frozen"]="two";return dict(authentication_evaluations=1)
    monkeypatch.setattr(tool,"validate_record",drift)
    rejected=tmp_path/"rejected"
    with pytest.raises(RuntimeError,match="authentication") as error:tool.create(rejected)
    assert not rejected.exists()
    assert json.loads(open(error.value.failure_receipt).read())["status"]=="FAILED_NOT_COMPLETED"


def test_calculation_failure_preserves_case_and_attempt_state(synthetic,budget,monkeypatch):
    def fail(*args,**kwargs):
        kwargs["budget"].check(1)
        error=RuntimeError("manufactured Courant failure")
        error.diagnostic=dict(phase="stage_2",accepted_state=dict(marker=1),attempt_state=dict(marker=2))
        raise error
    monkeypatch.setattr(tool,"evolve_stage_ssprk3",fail)
    with pytest.raises(RuntimeError) as error:tool.calculate(budget)
    receipt=error.value.operator_failure
    assert receipt["phase"]["amplitude"]=="zero" and receipt["stage"]["attempt_state"]==dict(marker=2)
    assert receipt["completed_cases"]==[]


@pytest.mark.parametrize("field,value",(("maximum_evaluations",80000001),("cpu_seconds",float("nan")),
    ("started",True),("started",-1.),("evaluations",-1)))
def test_invalid_supplied_budget_rejected(field,value):
    child=tool.Budget(maximum_evaluations=10,cpu_seconds=30.)
    setattr(child,field,value)
    with pytest.raises(ValueError):tool.shared_budget(child)


def test_exhausted_caller_budget_is_not_reset(budget,monkeypatch):
    budget.check(2)
    child=tool.Budget(maximum_evaluations=1,cpu_seconds=30.)
    def load(shared):shared.check(2);raise AssertionError("must stop before loading")
    monkeypatch.setattr(tool,"load_inputs",load)
    with pytest.raises(RuntimeError,match="budget exhausted"):tool.calculate(child)


ORIGINAL_LOAD = tool.load_inputs
