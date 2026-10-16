"""Manufactured SSPRK3 timing/order and coupled controls; no prepared solve."""
import copy
import time

import numpy as np
import pytest

import recursive_horizons.nsc_atlas_background_ssprk3 as model
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_compatible_evolution import quiet_departures
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge
from recursive_horizons.nsc_atlas_background_cap import background_cap_metric, background_cap_rhs
from test_nsc_atlas_reference_wall_gauge import GEOMETRY

TUBE={"r":[.2,.8],"L":[-2.,2.],"c":[-4.,4.]}
ETA,Y=np.linspace(0.,1.,5),np.linspace(0.,1.,4)


@pytest.fixture(scope="module")
def budget():
    value=Budget(maximum_evaluations=150000,cpu_seconds=30.)
    yield value
    value.check()
    print("SSPRK3 manufactured controls",value.evaluations,time.process_time()-value.started)


def initial(y0=0.):
    state=quiet_departures(exterior_cells=4,cap_cells=3)
    state["wall"]["s"]=y0
    return state


def options(budget):
    return dict(eta=ETA,y_faces=Y,source=lambda x:dict(x=x),
        gauge=lambda x:dict(x=x,r=.4,T=.1,L=0.,s=0.),alpha=0.,tube=TUBE,budget=budget,profile=None)


def ode_stage(rhs,calls):
    def evaluate(state,*,x,source,gauge,budget,**kwargs):
        calls.append((x,copy.deepcopy(state),source,gauge,budget))
        rates=initial(); rates["wall"]["s"]=rhs(x,float(state["wall"]["s"]))
        wall={k:float(gauge[k])+float(state["wall"][k]) for k in ("r","T","L","s")}
        wall["c"]=0.
        budget.check(1)
        return dict(rates=rates,courant_rate=1.,physical_wall=wall,
            physical_cap=copy.deepcopy(state["cap"]),physical_increment=copy.deepcopy(state["departure"]))
    return evaluate


def test_nonautonomous_quadratic_rate_and_exact_stage_keys(monkeypatch,budget):
    calls=[]; monkeypatch.setattr(model,"background_stage",ode_stage(lambda x,y:x*x,calls))
    left,right=.2,.4
    mid=float(np.nextafter(left+(right-left)/2,1.))
    keys={left,right,mid}; lookups=[]
    kw=options(budget)
    def source(at):
        assert at in keys; lookups.append(("source",at)); return dict(x=at)
    def gauge(at):
        assert at in keys; lookups.append(("gauge",at)); return dict(x=at,r=.4,T=.1,L=0.,s=0.)
    kw.update(source=source,gauge=gauge)
    before=initial(.7)
    row=model.advance_background_ssprk3(before,left,mid,right,**kw)
    assert row["stage_times"]==[left,right,mid]
    assert [c[0] for c in calls]==[left,right,mid]
    assert lookups==[(name,x) for x in (left,right,mid) for name in ("source","gauge")]
    assert float(row["accepted_state"]["wall"]["s"])==pytest.approx(.7+(right**3-left**3)/3,abs=2e-16)
    assert row["stage_calls"]==3 and row["retries"]==0
    assert float(before["wall"]["s"])==.7
    assert all(c[-1] is budget for c in calls)


@pytest.mark.parametrize("kind",("autonomous_nonlinear","mixed_time_state"))
def test_measured_third_order_for_nonlinear_and_mixed_equations(monkeypatch,budget,kind):
    if kind=="autonomous_nonlinear":
        rhs=lambda x,y:y*y; left,right,y0=0.,.3,.7; exact=y0/(1-y0*(right-left))
    else:
        rhs=lambda x,y:x*y; left,right,y0=.1,.6,.7; exact=y0*np.exp((right*right-left*left)/2)
    errors=[]
    for n in (8,16,32):
        calls=[]; monkeypatch.setattr(model,"background_stage",ode_stage(rhs,calls))
        fine=np.linspace(left,right,2*n+1)
        row=model.evolve_background_ssprk3(initial(y0),fine[::2],mid_nodes=fine[1::2],save_every=n,**options(budget))
        actual=float(row["snapshots"][-1]["state"]["wall"]["s"])
        errors.append(abs(actual-exact))
        assert row["stages"]==len(calls)==3*n+1
        assert len(row["snapshots"])==2
    ratios=np.array(errors[:-1])/errors[1:]
    print(kind,"errors",errors,"ratios",ratios.tolist())
    assert np.all(ratios>7.5) and np.all(ratios<8.5)


def test_switching_rhs_is_measured_without_claiming_third_order(monkeypatch,budget):
    # A reconstruction-style max switch is only Lipschitz. A smooth ODE
    # order test does not license third-order convergence through this front.
    errors=[]
    for n in (8,16,32):
        calls=[]
        monkeypatch.setattr(model,"background_stage",ode_stage(lambda x,y:1+max(y,0.),calls))
        fine=np.linspace(0.,.3,2*n+1)
        row=model.evolve_background_ssprk3(initial(-.11),fine[::2],mid_nodes=fine[1::2],save_every=n,**options(budget))
        actual=float(row["snapshots"][-1]["state"]["wall"]["s"])
        errors.append(abs(actual-(np.exp(.19)-1)))
    print("switching RHS errors",errors)
    assert errors[-1]<errors[0]
    assert max(errors)<2e-4


def test_curved_reference_restored_at_each_nonmonotone_stage(monkeypatch,budget):
    calls=[]
    base=ode_stage(lambda x,y:0.,calls)
    def ref(x): return dict(x=x,r=.4+.1*x+.02*x**3,T=.1+.2*x+.05*x**3,L=0.,s=0.)
    def stage(state,**kwargs):
        row=base(state,**kwargs); x=kwargs["x"]
        row["rates"]["wall"]["r"]=-.06*x*x
        row["rates"]["wall"]["T"]=-.15*x*x
        return row
    monkeypatch.setattr(model,"background_stage",stage)
    kw=options(budget); kw["gauge"]=ref
    row=model.advance_background_ssprk3(initial(),0.,.02,.04,**kw)
    # Physical rates are r_x=.1,T_x=.2. Restored predictor geometry still
    # contains the known reference cubic remainder, then final RK weights
    # integrate the raw quadratic departure rates exactly.
    right_stage=row["stages"][1]["physical_wall"]
    assert right_stage["r"]==pytest.approx(.4+.1*.04+.02*.04**3,abs=1e-16)
    accepted=row["accepted_state"]["wall"]
    assert ref(.04)["r"]+accepted["r"]==pytest.approx(.404,abs=1e-16)
    assert ref(.04)["T"]+accepted["T"]==pytest.approx(.108,abs=1e-16)
    assert [x for x,*_ in calls]==[0.,.04,.02]


def test_matched_cubic_background_preserves_common_flux_and_zero_departures(monkeypatch,budget):
    labels=np.array([0.,.3,.6,1.])
    profile=CharacteristicProfile(labels,labels**3,3*labels**2,6*labels)
    def gauge(x): return dict(x=x,r=.4,T=.1+.2*x,L=0.,s=0.)
    def stage(state,*,x,gauge,profile,alpha,budget,**kwargs):
        R,T=gauge["r"],gauge["T"]
        metric=background_cap_metric(Y,R,state["cap"],T=T,profile=profile,alpha=alpha,budget=budget)
        bg=metric["background_faces"]; fW=metric["wall"]["f"]
        # Prescribe only the incoming reduced characteristic; the common
        # wall value remains background plus the evolved conservative state.
        incoming=(bg["p"][-1]+bg["q"][-1]-fW)/R
        cap=background_cap_rhs(metric,T_x=.2,radius_x=0.,incoming_w_plus=incoming,budget=budget)
        rates=initial(); rates["cap"]=cap["rates"]
        return dict(rates=rates,courant_rate=cap["characteristic_courant_rate"],
            physical_wall=dict(r=R,T=T,L=0.,s=0.,c=metric["wall"]["c"]),
            physical_cap=metric["physical_integrals"],physical_increment=copy.deepcopy(state["departure"]),
            wall_chi_x=cap["wall_chi_x"])
    monkeypatch.setattr(model,"background_stage",stage)
    kw=options(budget); kw.update(gauge=gauge,profile=profile)
    row=model.advance_background_ssprk3(initial(),0.,.02,.04,**kw)
    for values in row["accepted_state"].values():
        for v in values.values(): np.testing.assert_allclose(v,0.,atol=5e-17)
    stage_flux=sum(w*s["wall_chi_x"] for w,s in zip((1/6,1/6,2/3),row["stages"]))* .04
    B0=(.1+np.arctanh(.4))**3; B1=(.108+np.arctanh(.4))**3
    assert stage_flux==pytest.approx(B1-B0,abs=2e-16)
    # The zero profile-chain remainder follows from quadratic rate accuracy,
    # not an after-step correction to the scalar or its common wall value.
    assert abs(B1-B0-stage_flux)<2e-16


def source_quiet(x):
    return dict(r=x,log_C=0.,c=0.,r_x=1.,log_C_x=0.,c_x=0.)


@pytest.mark.parametrize("alpha",(0.,.03))
def test_actual_frozen_operator_quiet_endpoint_and_snapshot_contract(budget,alpha):
    ref=ReferenceWallGauge(GEOMETRY,budget=budget)
    fine=np.linspace(.4,.4008,9)
    end=ref.sample(float(fine[-1]))
    last=float(end["T"]+np.arctanh(end["r"]))
    quiet=CharacteristicProfile([0.,last],[0.,0.],[0.,0.],[0.,0.])
    row=model.evolve_background_ssprk3(initial(),fine[::2],mid_nodes=fine[1::2],
        eta=ETA,y_faces=Y,source=source_quiet,gauge=ref.sample,alpha=alpha,tube=TUBE,
        budget=budget,profile=quiet,save_every=2)
    assert row["stages"]==13 and row["steps"]==4 and row["retries"]==0
    assert [s["x"] for s in row["snapshots"]]==list(fine[::4])
    for sample in row["snapshots"]:
        assert set(sample)=={"x","state","physical"}
        assert set(sample["physical"])=={"increment","wall","cap"}
        np.testing.assert_array_equal(sample["state"]["cap"]["Iq"],0.)
        np.testing.assert_array_equal(sample["physical"]["cap"]["Iq"],0.)
        assert sample["physical"]["wall"]["c"]==0.
    final=row["final_stage"]["physical_wall"]
    assert final["T"]+np.arctanh(final["r"])<=last
    assert row["max_courant"]<.4


@pytest.mark.parametrize("left,mid,right",((0.,.51,1.),(1.,.5,0.),(0.,0.,1.),(0.,True,1.),(0.,.5,float("inf"))))
def test_bad_interval_is_zero_stage_failure(left,mid,right,budget):
    with pytest.raises(model.CoupledEvolutionError) as found:
        model.advance_background_ssprk3(initial(),left,mid,right,**options(budget))
    assert found.value.diagnostic["stage_calls"]==0
    assert found.value.diagnostic["phase"]=="validation"


@pytest.mark.parametrize("kind",("operator","CFL","lookup","budget"))
def test_failure_counts_actual_attempts_without_retry(monkeypatch,budget,kind):
    calls=[]; base=ode_stage(lambda x,y:y,calls)
    def stage(state,**kwargs):
        if kind=="operator" and len(calls)==1:
            calls.append("failed"); raise ValueError("second operator failed")
        row=base(state,**kwargs)
        if kind=="CFL" and len(calls)==3: row["courant_rate"]=100.
        return row
    monkeypatch.setattr(model,"background_stage",stage)
    kw=options(budget)
    if kind=="lookup":
        def source(x):
            if x==.1: raise KeyError("missing exact right cache key")
            return dict(x=x)
        kw["source"]=source
    if kind=="budget": kw["budget"]=Budget(maximum_evaluations=1,cpu_seconds=30.)
    with pytest.raises(model.CoupledEvolutionError) as found:
        model.advance_background_ssprk3(initial(.7),0.,.05,.1,**kw)
    expected={"operator":2,"CFL":3,"lookup":1,"budget":2}[kind]
    assert found.value.diagnostic["stage_calls"]==expected
    assert len(calls)==expected


def test_midpoint_failure_in_later_step_is_validated_before_any_operator(monkeypatch,budget):
    calls=[]; monkeypatch.setattr(model,"background_stage",ode_stage(lambda x,y:y,calls))
    with pytest.raises(model.CoupledEvolutionError) as found:
        model.evolve_background_ssprk3(initial(),[0.,.1,.2],mid_nodes=[.05,.151],**options(budget))
    assert found.value.diagnostic["stage_calls"]==0 and not calls


def test_full_evolution_reports_cumulative_failed_stage_calls(monkeypatch,budget):
    calls=[]; base=ode_stage(lambda x,y:y,calls)
    def stage(state,**kwargs):
        if len(calls)==4:
            calls.append("failed"); raise ValueError("failure during second step")
        return base(state,**kwargs)
    monkeypatch.setattr(model,"background_stage",stage)
    with pytest.raises(model.CoupledEvolutionError) as found:
        model.evolve_background_ssprk3(initial(.7),[0.,.1,.2],mid_nodes=[.05,.15],**options(budget))
    diagnostic=found.value.diagnostic
    assert diagnostic["stage_calls"]==len(calls)==5
    assert diagnostic["step"]==1 and diagnostic["local_failure"]["stage_calls"]==2


def test_profile_domain_failure_is_preserved_at_right_stage(budget):
    ref=ReferenceWallGauge(GEOMETRY,budget=budget)
    start=ref.sample(.4)
    last=float(start["T"]+np.arctanh(start["r"]))
    quiet=CharacteristicProfile([0.,last],[0.,0.],[0.,0.],[0.,0.])
    with pytest.raises(model.CoupledEvolutionError,match="finite labels") as found:
        model.advance_background_ssprk3(initial(),.4,.4001,.4002,eta=ETA,y_faces=Y,
            source=source_quiet,gauge=ref.sample,alpha=0.,tube=TUBE,budget=budget,profile=quiet)
    assert found.value.diagnostic["stage_calls"]==2
