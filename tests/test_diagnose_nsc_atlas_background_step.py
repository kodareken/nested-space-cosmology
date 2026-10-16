"""Independent local-step identities and no-replay controls; no saved-data run."""
import copy
import json
import time

import numpy as np
import pytest

import diagnose_nsc_atlas_background_step as tool
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_compatible_evolution import quiet_departures


@pytest.fixture(scope="module")
def budget():
    value=tool.Budget(maximum_evaluations=100000,cpu_seconds=30.)
    yield value
    value.check()
    print("actual-step manufactured controls",value.evaluations,time.process_time()-value.started)


def profile():
    labels=np.array([0.,.3,.6,1.2])
    return CharacteristicProfile(labels,labels**3,3*labels**2,6*labels)


def cache(x,simple=False):
    x=np.array(x)
    if simple:
        r=np.full_like(x,.4); T=.1+.2*x
        rx=np.zeros_like(x); Tx=np.full_like(x,.2)
    else:
        r=.4+.1*x+.02*x**3; T=.1+.2*x+.05*x**3
        rx=.1+.06*x*x; Tx=.2+.15*x*x
    return dict(x=x,r=r,T=T,L=.1*x,s=.3*x,r_x_total=rx,T_x=Tx,
                L_x_total=np.full_like(x,.1),s_x=np.full_like(x,.3))


def manufactured_stage(calls,*,simple=False):
    def stage(state,*,x,gauge,source,profile,alpha,budget,**kwargs):
        calls.append((x,copy.deepcopy(state),budget))
        view=tool.compact_state(state,gauge,profile,alpha,budget)
        wall=view["wall"]; R,T=wall["r"],wall["T"]
        chi=view["canonical_chi"]; f=chi/R; u=wall["c"]/R
        D=1+alpha*u*u/6; xs=1.2
        Rx=0. if simple else .1+.07*x+.03*R
        Tx=.2 if simple else .2+.03*T+.01*x
        bg=profile.background(T,R,budget=budget)
        Bx=float(Tx*bg["A0"]*bg["p"]+Rx*bg["q"])
        Jchi=Bx if simple else .3+.07*chi+.03*R+.02*x*x
        returned=0. if simple else .01*(1+x)+.004*chi
        # Independent derivative of C(R,chi)=R*sinh(k*chi/R)/k.
        Jc=u*Rx+np.sqrt(D)*(Jchi-f*Rx)
        rates=quiet_departures(exterior_cells=len(state["departure"]["c"])+1,
                               cap_cells=len(state["cap"]["Iq"]))
        rates["wall"].update(r=Rx-gauge["r_x_total"],T=Tx-gauge["T_x"])
        rates["cap"]["Iq"][:]=(Jchi-Bx)/len(rates["cap"]["Iq"])
        for key in rates["departure"]:
            rates["departure"][key][:]=.01*state["departure"][key]
        budget.check(5)
        return dict(rates=rates,physical_wall=wall,courant_rate=1.,
            wall_rates=dict(r=Rx,T=Tx,L=.1,s=.3,c=Jc),
            hybrid=dict(Vdot_plus=xs,D=D,C_minus=2*xs*returned,f_s=xs*(Jchi-f*Rx)/R),
            exterior_longitudinal=dict(c=np.array([Jc-returned])),
            cap_metric=dict(wall=dict(f=f),background_faces=dict(p=np.array([float(bg["p"])]),
                                                                 q=np.array([float(bg["q"])]))),
            cap_rhs=dict(rates=dict(Iq=rates["cap"]["Iq"]),wall_chi_x=Jchi))
    return stage


def run_step(monkeypatch,budget,*,alpha=.3,simple=False):
    calls=[]; monkeypatch.setattr(tool,"background_stage",manufactured_stage(calls,simple=simple))
    initial=quiet_departures(exterior_cells=4,cap_cells=3)
    if not simple: initial["cap"]["Iq"][:]=.01
    rows=cache([0.,.02,.04],simple)
    options=dict(source_rows=dict(x=rows["x"]),gauge_rows=rows,profile=profile(),alpha=alpha,
                 budget=budget,counter=dict(stage_calls=0,max_courant=0.))
    accepted,row=tool.accepted_step(initial,0,2,label="manufactured",**options)
    return initial,accepted,row,options,calls


@pytest.mark.parametrize("alpha",(0.,.3))
def test_exact_accepted_step_split_including_reference_geometry(monkeypatch,budget,alpha):
    initial,accepted,row,options,calls=run_step(monkeypatch,budget,alpha=alpha)
    assert len(calls)==2 and all(c[2] is budget for c in calls)
    assert all(abs(v)<4e-16 for v in row["closures"].values())
    assert row["increment"]==pytest.approx(sum(row["components"].values()),abs=4e-16)
    h=row["h"]
    for name,coefficient in (("r",.02),("T",.05)):
        g=row["reference_geometry"][name]
        assert g["reference_defect"]==pytest.approx(-coefficient*h**3/2,abs=5e-17)
        assert g["actual_defect"]==pytest.approx(g["reference_defect"],abs=6e-17)
        assert abs(g["closure"])<6e-17
        assert row["predictor_geometry_defect"][name]==pytest.approx(coefficient*h**3,abs=7e-17)
    assert row["accepted"]["canonical_chi"]==pytest.approx(row["rk"]["chi"]+row["background_chain"],abs=1e-16)
    # The initial is retained and the predictor is the actual frozen combine.
    np.testing.assert_array_equal(initial["cap"]["Iq"],np.full(3,.01))
    assert row["predictor"]["raw_state_sha256"]!=row["initial"]["raw_state_sha256"]
    if alpha==0:
        assert abs(row["components"]["frame_remainder"])<1e-16
        assert abs(row["components"]["reference_geometry_remainder"])<1e-16
        assert row["components"]["background_remainder"]==pytest.approx(row["background_chain"],abs=1e-16)
    else:
        assert abs(row["components"]["frame_remainder"])>1e-9
        assert abs(row["components"]["reference_geometry_remainder"])>1e-11


def test_matched_profile_flux_and_background_defects_cancel_without_reset(monkeypatch,budget):
    initial,coarse_state,coarse,options,calls=run_step(monkeypatch,budget,alpha=0.,simple=True)
    midpoint,half1=tool.accepted_step(initial,0,1,label="first half",**options)
    fine_state,half2=tool.accepted_step(midpoint,1,2,label="second half",**options)
    exact=-(.2*.04)**3/2
    assert coarse["background_chain"]==pytest.approx(exact,abs=3e-16)
    assert half1["background_chain"]+half2["background_chain"]==pytest.approx(exact/4,abs=4e-16)
    assert coarse["accepted"]["wall"]["c"]==pytest.approx(half2["accepted"]["wall"]["c"],abs=2e-16)
    flux_change=sum(h["components"]["exterior_flux"] for h in (half1,half2))-coarse["components"]["exterior_flux"]
    background_change=sum(h["components"]["background_remainder"] for h in (half1,half2))-coarse["components"]["background_remainder"]
    assert abs(flux_change)>1e-7
    assert flux_change+background_change==pytest.approx(0.,abs=4e-16)
    np.testing.assert_allclose(fine_state["cap"]["Iq"],0.,atol=1e-18)
    assert len(calls)==6


def audit_fixture(tmp_path,monkeypatch,budget):
    rows=cache(np.linspace(0.,1.,513)); prof=profile()
    cases=[{} for _ in range(6)]; base=dict(gauge=dict(rows=rows),initial_data={})
    for case_id,label,alpha in ((1,"zero",0.),(4,"finite",.3)):
        state=quiet_departures(exterior_cells=32,cap_cells=32)
        state["cap"]["Iq"][:]=.001
        wall=tool.compact_state(state,tool._cache(rows,328),prof,alpha,budget)["wall"]
        cases[case_id]=dict(cells=32,steps=256,amplitude=label,alpha=alpha,
            snapshots=[None]*82+[dict(x=rows["x"][328],state=state,physical=dict(wall=wall))])
        base["initial_data"][label]=dict(snapshot_rows=dict(x=rows["x"]))
    record=dict(cases=cases,profile_inputs={},profile={})
    monkeypatch.setattr(tool.previous,"authenticate",lambda root:(record,base,dict(pinned=True)))
    def fit(inputs,shared): shared.check(3); return prof,{}
    monkeypatch.setattr(tool.producer,"fit_profile",fit)
    for name in tool.SOURCES:
        p=tmp_path/name; p.parent.mkdir(parents=True,exist_ok=True); p.write_text("synthetic bound source")
    calls=[]; monkeypatch.setattr(tool,"background_stage",manufactured_stage(calls))
    def forbidden(*args,**kwargs): raise AssertionError("full trajectory or source replay forbidden")
    for name in ("evolve_background","extract_profile_inputs","load_parent","validate_record"):
        monkeypatch.setattr(tool.producer,name,forbidden)
    return calls,record,base


def test_exact_twelve_stage_same_start_orchestration_and_accounting(tmp_path,monkeypatch,budget):
    calls,record,base=audit_fixture(tmp_path,monkeypatch,budget)
    initial=budget.evaluations
    row=tool.audit(tmp_path,budget=budget)
    assert len(calls)==row["execution"]["stage_calls"]==12
    assert row["execution"]["local_ssprk2_steps"]==6
    assert row["execution"]["source_integrations"]==row["execution"]["full_trajectory_replays"]==0
    assert row["execution"]["initial_evaluations"]==initial
    assert row["execution"]["evaluated_points"]==initial+row["execution"]["profile_fit_points"]+row["execution"]["local_step_points"]
    indices=[328,330,328,329,329,330]*2
    assert [x for x,_,_ in calls]==[base["gauge"]["rows"]["x"][i] for i in indices]
    for offset in (0,6):
        assert tool.encoded(tool.plain(calls[offset][1]))==tool.encoded(tool.plain(calls[offset+2][1]))
        assert tool.encoded(tool.plain(calls[offset+2][1]))!=tool.encoded(tool.plain(calls[offset+4][1]))
    assert abs(row["paired_closure"])<1e-15
    for case in row["cases"]:
        assert abs(case["comparison_closure"])<1e-15
        assert case["coarse"]["initial"]["raw_state_sha256"]==case["halves"][0]["initial"]["raw_state_sha256"]
        assert case["halves"][0]["accepted"]["raw_state_sha256"]==case["halves"][1]["initial"]["raw_state_sha256"]


def test_failure_context_stops_without_retry(tmp_path,monkeypatch,budget):
    calls,_,_=audit_fixture(tmp_path,monkeypatch,budget)
    def failing(*args,**kwargs): calls.append("attempt"); kwargs["budget"].check(1); raise ValueError("manufactured failure")
    monkeypatch.setattr(tool,"background_stage",failing)
    with pytest.raises(RuntimeError,match=r"case=1 amplitude=zero.*phase=initial stage.*stage_calls=1.*manufactured failure"):
        tool.audit(tmp_path,budget=budget)
    assert calls==["attempt"]


@pytest.mark.parametrize("failure",("CFL","ceiling","coordinate","coverage"))
def test_stage_courant_cache_and_coverage_guards(monkeypatch,budget,failure):
    calls=[]; original=manufactured_stage(calls)
    def stage(*args,**kwargs):
        row=original(*args,**kwargs)
        if failure=="CFL": row["courant_rate"]=1000.
        return row
    monkeypatch.setattr(tool,"background_stage",stage)
    rows=cache([0.,.02,.04]); sources=dict(x=rows["x"].copy())
    if failure=="coordinate": sources["x"][2]+=1e-15
    if failure=="coverage": rows["T"][2]=1.1
    counter=dict(stage_calls=12 if failure=="ceiling" else 0,max_courant=0.)
    with pytest.raises((ValueError,RuntimeError)):
        tool.accepted_step(quiet_departures(exterior_cells=4,cap_cells=3),0,2,
            source_rows=sources,gauge_rows=rows,profile=profile(),alpha=0.,budget=budget,counter=counter,label=failure)
    assert len(calls)<=1


@pytest.mark.parametrize("args",((0.,.1,0.),(.5,.1,-.1),(.5,complex(1),0.),(.5,.1,True)))
def test_frame_guards(args):
    with pytest.raises(ValueError): tool.frame(*args)


def test_finite_frame_formula_and_small_alpha():
    R,chi,alpha=.5,.13,.3
    k=np.sqrt(alpha/6)
    assert tool.frame(R,chi,alpha)==pytest.approx(R*np.sinh(k*chi/R)/k,abs=1e-17)
    assert tool.frame(R,chi,0.)==chi
    assert tool.frame(R,chi,1e-12)-chi==pytest.approx(1e-12*chi**3/(36*R**2),abs=2*np.spacing(chi))


def test_frozen_stage_and_combine_on_exact_quiet_manufactured_data(budget):
    from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge
    from test_nsc_atlas_reference_wall_gauge import GEOMETRY
    x=np.array([.4,.4001,.4002])
    gauge=ReferenceWallGauge(GEOMETRY,budget=budget).sample(x)
    source=dict(x=x,r=x,r_x=np.ones(3),log_C=np.zeros(3),log_C_x=np.zeros(3),
                c=np.zeros(3),c_x=np.zeros(3))
    quiet=CharacteristicProfile([0.,1.],[0.,0.],[0.,0.],[0.,0.])
    counter=dict(stage_calls=0,max_courant=0.)
    _,row=tool.accepted_step(quiet_departures(exterior_cells=4,cap_cells=3),0,2,
        source_rows=source,gauge_rows=gauge,profile=quiet,alpha=0.,budget=budget,
        counter=counter,label="quiet manufactured vacuum")
    assert counter["stage_calls"]==2
    assert row["increment"]==0.
    assert all(v==0. for v in row["components"].values())
    assert all(v==0. for v in row["closures"].values())
    assert counter["max_courant"]<.4


def test_readout_preserves_frozen_ordered_primitive_not_pairwise_sum(budget):
    initial=quiet_departures(exterior_cells=4,cap_cells=32)
    initial["cap"]["Iq"]=np.array([1.,1e-16,-1.,1e-16]*8)
    ordered=float(np.cumsum(initial["cap"]["Iq"])[-1])
    assert ordered!=float(np.sum(initial["cap"]["Iq"]))
    quiet=CharacteristicProfile([0.,1.],[0.,0.],[0.,0.],[0.,0.])
    view=tool.compact_state(initial,tool._cache(cache([0.,.1]),0),quiet,0.,budget)
    assert view["departure_chi"]==ordered
    assert view["canonical_chi"]==ordered
    assert view["wall"]["c"]==tool.frame(.4,ordered,0.)


def test_creation_only_and_check(tmp_path,monkeypatch):
    result=dict(schema=tool.SCHEMA,execution=dict(cpu_seconds=.01,evaluated_points=123),ledger=.12)
    calls=[]
    def fake(root): calls.append(root); return copy.deepcopy(result)
    monkeypatch.setattr(tool,"audit",fake)
    path=tmp_path/"step.json"; tool.create(path,root=tmp_path)
    result["execution"]["cpu_seconds"]=.02
    assert tool.check(path,root=tmp_path)["stage_calls"]==12
    with pytest.raises(FileExistsError): tool.create(path,root=tmp_path)
    assert len(calls)==2
    bad=json.loads(path.read_bytes()); bad["ledger"]+=.01; path.write_bytes(tool.encoded(bad))
    with pytest.raises(ValueError,match="recomputation mismatch"): tool.check(path,root=tmp_path)


@pytest.mark.parametrize("key,value",(("maximum_evaluations",250001),("maximum_evaluations",True),
    ("cpu_seconds",31.),("evaluations",-1),("started",float("nan"))))
def test_cumulative_allowance_guards(key,value):
    budget=tool.Budget(maximum_evaluations=100000,cpu_seconds=30.)
    setattr(budget,key,value)
    with pytest.raises(ValueError,match="allowance"): tool._shared(budget)
