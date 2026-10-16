"""Manufactured chain identities and saved-input guards; no physical replay."""
import copy
import json
import time

import numpy as np
import pytest

import diagnose_nsc_atlas_background_time as tool
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile


@pytest.fixture(scope="module")
def budget():
    value=tool.Budget(maximum_evaluations=100000,cpu_seconds=30.)
    yield value
    value.check()
    print("background time manufactured controls",value.evaluations,time.process_time()-value.started)


def cubic(knots=(0.,.333,.613,1.)):
    x=np.array(knots)
    return CharacteristicProfile(x,x**3,3*x*x,6*x)


def test_cubic_signed_one_step_and_two_half_oracle(budget):
    Tx,h=.2,.1
    row=tool.linear_profile_path(cubic(),T=.1,R=.4,T_x=Tx,R_x=0.,h=h,budget=budget)
    exact=-(Tx*h)**3/2
    assert row["one_step_defect"]==pytest.approx(exact,abs=2e-16)
    assert row["two_half_steps_defect"]==pytest.approx(exact/4,abs=2e-16)
    assert row["fine_minus_coarse"]==pytest.approx(-3*exact/4,abs=2e-16)
    assert not row["crossings"]
    assert row["exact_endpoint_increment"]==pytest.approx(
        (.1+Tx*h+np.arctanh(.4))**3-(.1+np.arctanh(.4))**3,abs=3e-16)


def test_reflected_quadratic_profile_and_moving_radius(budget):
    knots=np.array([-.3,0.,.3,.8])
    profile=CharacteristicProfile(knots,knots**2,2*knots,np.full(4,2.))
    T,R,Tx,Rx,h=.2,.4,.3,.05,.1
    row=tool.linear_profile_path(profile,T=T,R=R,T_x=Tx,R_x=Rx,h=h,budget=budget)
    tau=np.array([0.,h/2,h])
    tt,rr=T+Tx*tau,R+Rx*tau
    B=4*tt*np.arctanh(rr)
    Bx=4*Tx*np.arctanh(rr)+4*tt*Rx/(1-rr**2)
    np.testing.assert_allclose(row["primitive"],B,atol=2e-15,rtol=0.)
    np.testing.assert_allclose(row["rate"],Bx,atol=3e-15,rtol=0.)
    assert row["one_step_defect"]==pytest.approx(B[-1]-B[0]-h*(Bx[0]+Bx[-1])/2,abs=3e-15)


def test_reference_signed_cumulative_partition_and_time_difference(budget):
    x=np.linspace(0.,1.,513)
    R=.4
    gauge=dict(x=x,T=.2+.5*x-np.arctanh(R),r=np.full(513,R),
               T_x=np.full(513,.5),r_x_total=np.zeros(513))
    # Independent cubic trapezoid remainder, at the common 129 cuts.
    expected={str(n):-(.5**3/2)*x[::4]/n**2 for n in (256,512)}
    change=expected["512"]-expected["256"]
    actual=change+.01*x[::4]
    row=tool.reference_chain(cubic(),gauge,actual,budget)
    for n in (256,512):
        grid=row["grids"][str(n)]
        np.testing.assert_allclose(grid["R0"],expected[str(n)],atol=3e-15,rtol=0.)
        np.testing.assert_allclose(np.array(grid["crossing_cumulative"])+grid["within_piece_cumulative"],
                                   grid["R0"],atol=1e-20,rtol=0.)
        assert grid["crossing_count"]==2
        assert grid["within_piece_count"]==n-2
        assert grid["partition_defect"]<1e-20
    np.testing.assert_allclose(row["known_chain_time_difference"],change,atol=5e-15,rtol=0.)
    np.testing.assert_allclose(row["actual_minus_known_chain"],.01*x[::4],atol=5e-15,rtol=0.)


def test_crossings_use_strict_interior_and_report_null_branch(budget):
    profile=cubic((0.,.53,1.))
    row=tool.linear_profile_path(profile,T=.1,R=.4,T_x=.2,R_x=0.,h=.1,budget=budget)
    assert len(row["crossings"])==1
    event=row["crossings"][0]
    assert event["branch"]=="incoming"
    assert event["fraction"]==pytest.approx((.53-.1-np.arctanh(.4))/.02,abs=2e-13)
    mask=tool._crossed(profile,np.array([0.,.2,.53]),np.array([.53,.8,1.]),budget)
    assert mask.tolist()==[False,True,False]


def test_nonlinear_ssprk2_chain_uses_accepted_geometry(budget):
    # Two independently supplied geometry stages: the endpoint is the RK
    # average of rates, while B is evaluated at that accepted geometry.
    profile=cubic()
    X=np.array([.1,.4]); k0=np.array([.2,.04]); kp=np.array([.31,-.01]); h=.1
    predictor=X+h*k0; accepted=X+h*(k0+kp)/2
    def value_rate(point,rate):
        bg=profile.background(point[0],point[1],budget=budget)
        return float(bg["chi"]),float(rate[0]*bg["A0"]*bg["p"]+rate[1]*bg["q"])
    B0,bx0=value_rate(X,k0)
    Bp,bxp=value_rate(predictor,kp)
    B1,_=value_rate(accepted,kp)
    J0,Jp,Z0=.17,-.08,.03
    Z1=Z0+h*((J0-bx0)+(Jp-bxp))/2
    residual=B1-B0-h*(bx0+bxp)/2
    assert (B1+Z1)-(B0+Z0)==pytest.approx(h*(J0+Jp)/2+residual,abs=1e-16)
    # Averaging endpoint B values would violate this nonlinear chain.
    assert abs(B1-(B0+Bp)/2)>1e-3
    assert abs(residual)>1e-6


def stage_fixture(alpha):
    R,f,xs,rs,fs,Tx=.5,.2,1.7,.2,.4,.8
    k=np.sqrt(alpha/6)
    u=np.sinh(k*f)/k if alpha else f
    D=1+alpha*u*u/6
    Rx=rs/xs
    Jchi=Rx*f+R*fs/xs
    Jc=u*Rx+np.sqrt(D)*(Jchi-Rx*f)
    p0,q0=.1,.2
    Bx=Tx*(1-R*R)*p0+Rx*q0
    return dict(physical_wall=dict(r=R,T=.1,c=R*u),wall_rates=dict(r=Rx,T=Tx,c=Jc),
        hybrid=dict(Vdot_plus=xs,D=D,C_minus=.04*xs,f_s=fs),
        cap_metric=dict(wall=dict(f=f),background_faces=dict(p=[0.,p0],q=[0.,q0])),
        cap_rhs=dict(rates=dict(Iq=[Jchi-Bx,0.]),wall_chi_x=Jchi),
        exterior_longitudinal=dict(c=[0.,Jc-.02]),physical_cap=dict(Ip=[.1,.2],Iq=[.3,.4]))


@pytest.mark.parametrize("alpha",(0.,.3))
def test_invariant_canonical_and_departure_wall_chain(alpha,budget):
    stage=stage_fixture(alpha)
    row=tool.wall_chain(stage,budget)
    for key in ("invariant_split_defect","canonical_frame_defect","departure_chain_defect","canonical_flux_defect"):
        assert abs(row[key])<4e-17
    stage["exterior_longitudinal"]["c"][-1]+=.012
    assert tool.wall_chain(stage,budget)["invariant_split_defect"]==pytest.approx(-.012)


@pytest.mark.parametrize("change",("upper","negative_step","spacelike","radius","complex","bool"))
def test_linear_path_domain_is_not_extrapolated(change,budget):
    kw=dict(T=.1,R=.4,T_x=.2,R_x=0.,h=.1)
    if change=="upper": kw["T"]=.9
    elif change=="negative_step": kw["h"]=-.1
    elif change=="spacelike": kw["R_x"]=.3
    elif change=="radius": kw["R"]=1.
    elif change=="complex": kw["T"]=1j
    else: kw["h"]=True
    with pytest.raises(ValueError): tool.linear_profile_path(cubic(),budget=budget,**kw)


@pytest.mark.parametrize("change",("size","nonuniform","reversed","bool","actual"))
def test_reference_grid_controls(change,budget):
    x=np.linspace(0.,1.,513)
    gauge=dict(x=x,T=.2+.5*x-np.arctanh(.4),r=np.full(513,.4),T_x=np.full(513,.5),r_x_total=np.zeros(513))
    actual=np.zeros(129)
    if change=="size": gauge["T"]=gauge["T"][:-1]
    elif change=="nonuniform": gauge["x"][3]+=.0001
    elif change=="reversed": gauge["x"]=x[::-1]
    elif change=="bool": gauge["r"]=list(gauge["r"]); gauge["r"][0]=True
    else: actual=np.zeros(128)
    with pytest.raises(ValueError): tool.reference_chain(cubic(),gauge,actual,budget)


def binding_fixture(tmp_path,monkeypatch):
    names=(tool.BASE,*[f"bound-{i}.txt" for i in range(65)])
    bindings={}
    for name in names:
        path=tmp_path/name; path.parent.mkdir(parents=True,exist_ok=True)
        raw=b'{"synthetic_base":true}' if name==tool.BASE else name.encode()
        path.write_bytes(raw); bindings[name]=tool.digest(raw)
    record=dict(schema=tool.producer.SCHEMA,settings=tool.producer.SETTINGS)
    raw=tool.encoded(record)
    path=tmp_path/tool.INPUT; path.parent.mkdir(parents=True,exist_ok=True); path.write_bytes(raw)
    manifest=dict(schema=tool.producer.SCHEMA+"-manifest",source_bindings=bindings,
                  outputs={"background.json":dict(sha256=tool.digest(raw),bytes=len(raw))})
    mr=tool.encoded(manifest); path.with_name("manifest.json").write_bytes(mr)
    monkeypatch.setattr(tool,"INPUT_SHA256",tool.digest(raw))
    monkeypatch.setattr(tool,"MANIFEST_SHA256",tool.digest(mr))
    monkeypatch.setattr(tool.producer,"SOURCES",names)
    return path,manifest


def forbidden(*args,**kwargs):
    raise AssertionError("no source integration, time evolution or recursive validation")


def test_authentication_checks_only_hashes(tmp_path,monkeypatch):
    binding_fixture(tmp_path,monkeypatch)
    for name in ("load_parent","validate_record","extract_profile_inputs","evolve_background"):
        monkeypatch.setattr(tool.producer,name,forbidden)
    _,base,identity=tool.authenticate(tmp_path)
    assert base==dict(synthetic_base=True)
    assert identity["authenticated_input_files"]==66 and identity["hashed_bytes"]>0
    assert identity["recursive_validation"] is False


@pytest.mark.parametrize("change",("payload","manifest","bound","inventory","rehashed_bound"))
def test_authentication_tampering_rejected(tmp_path,monkeypatch,change):
    path,manifest=binding_fixture(tmp_path,monkeypatch)
    if change=="payload": path.write_bytes(path.read_bytes()+b" ")
    elif change=="manifest": manifest["schema"]="wrong"
    elif change=="inventory": (path.parent/"extra.txt").write_text("extra")
    else:
        (tmp_path/"bound-0.txt").write_text("changed")
        if change=="rehashed_bound": manifest["source_bindings"]["bound-0.txt"]=tool.digest(b"changed")
    path.with_name("manifest.json").write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError): tool.authenticate(tmp_path)


def orchestrated_fixture(tmp_path,monkeypatch):
    calls=[]
    cases=[dict(marker=i,amplitude="zero" if i<3 else "finite",steps=512 if i in (2,5) else 256) for i in range(6)]
    record=dict(cases=cases,profile_inputs={},profile={},diagnostics=dict(time_zero=dict(wall_arrays=dict(c=[]))))
    monkeypatch.setattr(tool,"authenticate",lambda root:(record,dict(gauge=dict(rows={})),dict(fixed=True)))
    for name in tool.OWN_SOURCES:
        p=tmp_path/name; p.parent.mkdir(parents=True,exist_ok=True); p.write_text("fixture")
    def fit(inputs,budget): budget.check(2); return object(),{}
    def reference(profile,gauge,actual,budget): budget.check(3); return {}
    def stage(case,base,profile,index,budget):
        calls.append((case["marker"],index,budget)); budget.check(4)
        keys=("Q_c","cap_c_return","Jc","Jchi","B_x","departure_chi_rate")
        return dict(chain={k:float(case["marker"]**2) for k in keys},evaluated_points=4)
    monkeypatch.setattr(tool.producer,"fit_profile",fit)
    monkeypatch.setattr(tool,"reference_chain",reference)
    monkeypatch.setattr(tool,"frozen_stage",stage)
    for name in ("load_parent","validate_record","extract_profile_inputs","evolve_background"):
        monkeypatch.setattr(tool.producer,name,forbidden)
    return calls


def test_orchestration_preserves_prior_counts_and_shared_budget(tmp_path,monkeypatch,budget):
    calls=orchestrated_fixture(tmp_path,monkeypatch)
    initial=budget.evaluations
    row=tool.audit(tmp_path,budget=budget)
    execution=row["execution"]
    assert execution["initial_evaluations"]==initial
    assert execution["evaluated_points"]==initial+21
    assert execution["profile_fit_points"]+execution["reference_points"]+execution["frozen_stage_points"]==21
    assert [(a,b) for a,b,_ in calls]==[(1,82),(2,82),(4,82),(5,82)]
    assert all(b is budget for _,_,b in calls)
    assert execution["source_integrations"]==execution["time_evolutions"]==execution["profile_stage_extractions"]==0
    assert row["paired_time_rate_change"]["Jc"]==6.


def test_failure_stops_with_case_and_cumulative_context(tmp_path,monkeypatch,budget):
    calls=orchestrated_fixture(tmp_path,monkeypatch)
    def failing(case,base,profile,index,shared):
        calls.append(case["marker"]); shared.check(2)
        raise ValueError("synthetic retained state failure")
    monkeypatch.setattr(tool,"frozen_stage",failing)
    with pytest.raises(RuntimeError,match=r"case=1 amplitude=zero steps=256 cut=82; cumulative_points=.*synthetic retained state failure"):
        tool.audit(tmp_path,budget=budget)
    assert calls==[1]


def test_one_frozen_stage_uses_exact_cached_coordinate_without_evolution(monkeypatch,budget):
    stage=stage_fixture(.3)
    x=np.linspace(.4,.5,513); sample=dict(x=x[0],state=dict(unchanged=True),
        physical=dict(wall=stage["physical_wall"],cap=stage["physical_cap"]))
    case=dict(cells=32,steps=256,amplitude="finite",alpha=.3,snapshots=[sample])
    base=dict(gauge=dict(rows=dict(x=x)),initial_data=dict(finite=dict(snapshot_rows=dict(x=x))))
    calls=[]
    def fake(state,**kwargs):
        calls.append((copy.deepcopy(state),kwargs)); kwargs["budget"].check(7)
        return copy.deepcopy(stage)
    monkeypatch.setattr(tool.producer,"background_stage",fake)
    monkeypatch.setattr(tool.producer,"evolve_background",forbidden)
    row=tool.frozen_stage(case,base,cubic(),0,budget)
    assert len(calls)==1 and calls[0][1]["budget"] is budget
    assert calls[0][0]==sample["state"]==dict(unchanged=True)
    assert row["linear_path"]["h"]==x[2]-x[0]
    assert row["evaluated_points"]>7
    bad=copy.deepcopy(case); bad["snapshots"][0]["x"]+=1e-15
    with pytest.raises(ValueError,match="exact cached coordinate"):
        tool.frozen_stage(bad,base,cubic(),0,budget)
    assert len(calls)==1


def test_creation_check_recomputes_but_never_overwrites(tmp_path,monkeypatch):
    result=dict(schema=tool.SCHEMA,execution=dict(cpu_seconds=.01,evaluated_points=123),measured=-.3)
    calls=[]
    def audit(root): calls.append(root); return copy.deepcopy(result)
    monkeypatch.setattr(tool,"audit",audit)
    path=tmp_path/"audit.json"
    tool.create(path,root=tmp_path)
    result["execution"]["cpu_seconds"]=.02
    assert tool.check(path,root=tmp_path)["evaluated_points"]==123
    with pytest.raises(FileExistsError): tool.create(path,root=tmp_path)
    assert len(calls)==2
    saved=json.loads(path.read_bytes()); saved["measured"]+=.01; path.write_bytes(tool.encoded(saved))
    with pytest.raises(ValueError,match="recomputation mismatch"): tool.check(path,root=tmp_path)


@pytest.mark.parametrize("key,value",(("maximum_evaluations",150001),("maximum_evaluations",True),
    ("cpu_seconds",31.),("cpu_seconds",float("nan")),("evaluations",-1),("started",float("nan"))))
def test_audit_budget_cannot_be_disabled(key,value):
    budget=tool.Budget(maximum_evaluations=100000,cpu_seconds=30.)
    setattr(budget,key,value)
    with pytest.raises(ValueError,match="allowance"): tool._shared(budget)


def test_profile_queries_share_exhaustion():
    budget=tool.Budget(maximum_evaluations=1,cpu_seconds=30.)
    with pytest.raises(RuntimeError,match="budget exhausted"):
        tool.linear_profile_path(cubic(),T=.1,R=.4,T_x=.2,R_x=0.,h=.1,budget=budget)
    assert budget.evaluations>1
