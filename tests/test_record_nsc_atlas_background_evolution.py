"""Synthetic retained records and algebraic reconstruction, without evolution."""
import copy
import json
import time

import numpy as np
import pytest

import record_nsc_atlas_background_evolution as tool
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile


@pytest.fixture(autouse=True)
def finite_test_cpu():
    started=time.process_time()
    yield
    assert time.process_time()-started < 30.


def profile_inputs():
    z=np.linspace(0,.25,129)
    R=np.full(129,.75)
    F=np.maximum(z-z[8],0.)**3
    first=3*np.maximum(z-z[8],0.)**2
    walls=dict(T=1+z-np.arctanh(R),r=R,c=F)
    rates=dict(T=np.ones(129),r=np.zeros(129),c=first)
    table=tool.plain(dict(x=np.linspace(.4,.5,129),walls=walls,wall_rates=rates))
    return dict(table=table,table_sha256=tool.digest(tool.encoded(table)),
        stage_reconstructions=129,extraction_evaluations=129*(36*32+11))


def test_profile_training_and_withheld_are_separate():
    budget=tool.Budget(maximum_evaluations=100000,cpu_seconds=30.)
    inputs=profile_inputs()
    prof,summary=tool.fit_profile(inputs,budget)
    np.testing.assert_array_equal(summary["labels"],np.array(inputs["table"]["walls"]["T"])[::2]+np.arctanh(.75))
    assert len(summary["labels"])==65 and len(summary["withheld"]["indices"])==64
    changed=copy.deepcopy(inputs)
    changed["table"]["walls"]["c"][51]+=.01
    changed["table_sha256"]=tool.digest(tool.encoded(changed["table"]))
    _,other=tool.fit_profile(changed,budget)
    assert other["values"]==summary["values"] and other["first"]==summary["first"]
    assert other["withheld"]["value_difference"][25]-summary["withheld"]["value_difference"][25]==pytest.approx(-.01)
    assert "computational" in summary["scope"]
    assert prof.domain["last_label"]==summary["label_range"][-1]
    inputs["table"]["wall_rates"]["c"][17]+=.1
    with pytest.raises(ValueError,match="table hash"):
        tool.fit_profile(inputs,budget)


def test_profile_replay_conditions_quintic_cancellation_without_changing_inputs(monkeypatch):
    from recursive_horizons.nsc_atlas_validation import background_profile_roundoff_tolerances
    inputs = profile_inputs()
    # A smooth post-front cubic with a scale large enough that differentiated
    # polynomial cancellation exceeds a result-relative residual tolerance.
    for group in ("walls","wall_rates"):
        inputs["table"][group]["c"] = (4096*np.asarray(inputs["table"][group]["c"])).tolist()
    inputs["table_sha256"] = tool.digest(tool.encoded(inputs["table"]))
    original_inputs = tool.encoded(inputs)
    budget = tool.Budget(maximum_evaluations=100000,cpu_seconds=30.)
    _, expected = tool.fit_profile(inputs,budget)
    constructor = CharacteristicProfile.__init__
    def shifted_coefficients(self,*args,**kwargs):
        constructor(self,*args,**kwargs)
        coefficients = self._spline.c
        # Independent power-polynomial representation perturbation. Keep
        # zero coefficients zero, so the certified quiet branch is unchanged.
        active = coefficients != 0
        direction = np.where(np.indices(coefficients.shape).sum(axis=0)%2,np.inf,-np.inf)
        coefficients[active] = np.nextafter(coefficients[active],direction[active])
        self._derivatives = [self._spline.derivative(order) for order in range(1,6)]
    monkeypatch.setattr(CharacteristicProfile,"__init__",shifted_coefficients)
    _, replay = tool.fit_profile(inputs,budget)
    defect = abs(replay["diagnostics"]["max_abs_knot_second_defect"]-
                 expected["diagnostics"]["max_abs_knot_second_defect"])
    allowance = background_profile_roundoff_tolerances(inputs)
    assert defect > 2e-13
    assert defect <= allowance["diagnostics.max_abs_knot_second_defect"]
    tool.profile_replay_equal(replay,expected,inputs)
    assert tool.encoded(inputs) == original_inputs
    corrupt = copy.deepcopy(replay)
    corrupt["diagnostics"]["max_abs_knot_second_defect"] += .01
    with pytest.raises(ValueError,match=r"knot_second_defect.*max_abs_difference=.*allowance="):
        tool.profile_replay_equal(corrupt,expected,inputs)


def test_profile_replay_propagates_one_ulp_label_change_through_gradient(monkeypatch):
    inputs = profile_inputs()
    # Put the derived labels in the same binade as atanh(R), so addition
    # retains the one-ULP libm perturbation rather than rounding it away.
    inputs["table"]["walls"]["T"] = (np.asarray(inputs["table"]["walls"]["T"])-.5).tolist()
    inputs["table_sha256"] = tool.digest(tool.encoded(inputs["table"]))
    budget = tool.Budget(maximum_evaluations=100000,cpu_seconds=30.)
    _, expected = tool.fit_profile(inputs,budget)
    atanh = np.arctanh
    monkeypatch.setattr(np,"arctanh",lambda value:np.nextafter(atanh(value),np.inf))
    _, replay = tool.fit_profile(inputs,budget)
    assert replay["labels"] != expected["labels"]
    assert replay["values"] == expected["values"]
    tool.profile_replay_equal(replay,expected,inputs)


@pytest.mark.parametrize("field",("copied_value","index","quiet","withheld","nonfinite"))
def test_profile_replay_preserves_exact_data_and_rejects_material_changes(field):
    inputs = profile_inputs()
    _, expected = tool.fit_profile(inputs,tool.Budget(maximum_evaluations=100000,cpu_seconds=30.))
    changed = copy.deepcopy(expected)
    if field == "copied_value":
        changed["values"][30] = float(np.nextafter(changed["values"][30],np.inf))
    elif field == "index":
        changed["withheld"]["indices"][0] = 1.0
    elif field == "quiet":
        changed["diagnostics"]["quiet_anchor"]["first"] = 1e-30
    elif field == "withheld":
        changed["withheld"]["first_difference"][30] += 1e-4
    else:
        changed["second"][30] = float("nan")
    with pytest.raises(ValueError,match="profile replay"):
        tool.profile_replay_equal(changed,expected,inputs)


def test_physical_view_uses_total_profile_and_finite_frame():
    budget=tool.Budget(maximum_evaluations=100000,cpu_seconds=30.)
    labels=np.array([0.,.5,1.5])
    profile=CharacteristicProfile(labels,.03*labels**3,.09*labels**2,.18*labels)
    n,alpha,R,T=4,.02,.47,.1
    eta=np.linspace(0,1,n+1)
    state=tool.previous.quiet_departures(exterior_cells=n,cap_cells=n)
    state["cap"]["Iq"][:]=.001
    state["cap"]["Ip"][:]=.002
    source=dict(r=.44,log_C=0.,c=.01)
    gauge=dict(W=.2,r=R,L=-.1,T=T,s=.05)
    physical=tool.physical_view(state,.44,source,gauge,eta,alpha,profile,budget)
    plus,minus=T+np.arctanh(R),max(0.,T-np.arctanh(R))
    chi=.03*(plus**3-minus**3)+.004
    expected=R*np.sinh(np.sqrt(alpha/6)*chi/R)/np.sqrt(alpha/6)
    assert physical["wall"]["c"]==pytest.approx(expected,abs=2e-15)
    assert sum(physical["cap"]["Iq"])==pytest.approx(chi,abs=2e-15)
    assert abs(expected-R*np.sinh(np.sqrt(alpha/6)*.004/R)/np.sqrt(alpha/6))>1e-3
    exact=profile.integrals(T,R*eta,budget=budget)
    np.testing.assert_allclose(physical["cap"]["Ip"],exact["Ip"]+.002,atol=2e-16)


@pytest.fixture
def synthetic(monkeypatch):
    """Six tiny synthetic cases exercise the recorder, with evolution forbidden."""
    cases=tuple((label,3,128) for label in ("zero","zero","zero","finite","finite","finite"))
    monkeypatch.setattr(tool,"CASES",cases)
    x=np.linspace(.4,.41,513)
    base=dict(geometry={"synthetic":True},gauge=dict(rows=dict(x=x,W=np.full(513,.2),
        W_x=np.ones(513),W_xx=np.zeros(513),r=np.full(513,.45),L=np.zeros(513),T=x-.4,s=x-.4)))
    base["initial_data"]={label:dict(alpha=a,snapshot_rows=dict(x=x,r=x,log_C=np.zeros(513),
        c=np.zeros(513),r_x=np.ones(513),log_C_x=np.zeros(513),c_x=np.zeros(513)))
        for label,a in (("zero",0.),("finite",.01))}
    identity=dict(path="fixture",sha256="fixture",authentication={})
    def load(budget):
        budget.check(17)
        return {},base,identity,17
    monkeypatch.setattr(tool,"load_parent",load)
    monkeypatch.setattr(tool,"extract_profile_inputs",lambda *args:dict(extraction_evaluations=0))
    monkeypatch.setattr(tool,"validate_profile_inputs",lambda *args:None)
    profile=CharacteristicProfile([0.,1.5],[0.,0.],[0.,0.],[0.,0.])
    def fit(inputs,budget):
        budget.check(3)
        return profile,dict(synthetic_profile=True)
    monkeypatch.setattr(tool,"fit_profile",fit)
    def stage(state,*,budget,eta,**kwargs):
        budget.check(8*len(eta)+5)
        return dict(state=state)
    def summary(stage):
        return dict(cap_metric_evaluations=1,cap_rhs_evaluations=2,
                    checked_cap=tool.plain(stage["state"]["cap"]))
    monkeypatch.setattr(tool,"background_stage",stage)
    monkeypatch.setattr(tool,"final_summary",summary)
    def evolve(initial,x_nodes,*,eta,source,gauge,alpha,budget,profile,**kwargs):
        snapshots=[]
        for at in x[::4]:
            state=copy.deepcopy(initial)
            state["cap"]["Iq"][:]=alpha*(at-.4)**2
            state["cap"]["Ip"][:]=2*alpha*(at-.4)**2
            physical=tool.physical_view(state,float(at),source(float(at)),gauge(float(at)),eta,alpha,profile,budget)
            snapshots.append(dict(x=float(at),state=tool.plain(state),physical=physical))
        stages=2*(len(x_nodes)-1)+1
        # These are declared synthetic call counts, not a physical integration.
        source.calls=gauge.calls=stages
        budget.check(stages*(8*len(eta)+5))
        return dict(snapshots=snapshots,final_stage=dict(state=state),stages=stages,
                    max_courant=.1,retries=0)
    monkeypatch.setattr(tool,"evolve_background",evolve)
    def diag(record,parent,base,profile,budget):
        budget.check(23)
        return dict(caps=[c["snapshots"][-1]["physical"]["cap"] for c in record["cases"]])
    monkeypatch.setattr(tool,"diagnostics",diag)
    return base


def test_synthetic_calculation_and_authentication_share_one_budget(synthetic):
    budget=tool.Budget(maximum_evaluations=100000,cpu_seconds=30.)
    record=tool.calculate(budget)
    end=budget.evaluations
    auth=tool.validate_record(record,budget,require_creation=False)
    assert budget.evaluations==end+auth["authentication_evaluations"]
    assert auth["profile_stage_extraction_replayed"] is False
    assert auth["unsaved_stage_point_totals_independently_recomputed"] is False
    assert record["execution"]["source_integrations"]==record["execution"]["gauge_quadratures"]==0
    assert len(record["cases"])==6 and record["execution"]["coupled_evolutions"]==6
    assert budget.evaluations<100000


@pytest.mark.parametrize("field",("physical_cap","wall_c","raw_cap","final","counts","runtime","profile"))
def test_synthetic_authentication_rejects_derived_or_accounting_changes(synthetic,field):
    record=tool.calculate(tool.Budget(maximum_evaluations=100000,cpu_seconds=30.))
    if field=="physical_cap":
        record["cases"][3]["snapshots"][70]["physical"]["cap"]["Iq"][0]+=.01
    elif field=="wall_c":
        record["cases"][3]["snapshots"][70]["physical"]["wall"]["c"]+=.01
    elif field=="raw_cap":
        record["cases"][3]["snapshots"][70]["state"]["cap"]["Iq"][0]+=.01
    elif field=="final":
        record["cases"][0]["final"]["cap_metric_evaluations"]+=1
    elif field=="counts":
        record["cases"][0]["execution"]["measured_evaluations"]+=1
    elif field=="runtime":
        record["runtime"]["numpy"]=""
    else:
        record["profile"]["synthetic_profile"]=False
    with pytest.raises(ValueError):
        tool.validate_record(record,tool.Budget(maximum_evaluations=100000,cpu_seconds=30.),require_creation=False)


def test_cap_restriction_sums_total_integrals():
    table=dict(exterior={"c":np.array([[0.,1.,2.,3.,4.]])},wall={"r":np.array([.5])},
               cap={"Ip":np.array([[1.,2.,4.,8.]]),"Iq":np.array([[3.,5.,7.,9.]])})
    restricted=tool.previous._restricted(table)
    np.testing.assert_array_equal(restricted["cap"]["Ip"],[[3.,12.]])
    np.testing.assert_array_equal(restricted["cap"]["Iq"],[[8.,16.]])


def test_frozen_quadrature_holds_state_and_profile(monkeypatch):
    calls=[]
    def metric(y,R,state,*,T,profile,alpha,quadrature_order,budget):
        calls.append((id(state),R,T,id(profile),quadrature_order))
        budget.check(5)
        return dict(n_faces=np.array([0.,1+quadrature_order**-4]),wall=dict(delta=.1,A=.7,w_minus=.2),
                    gamma_interpolant_defect_nodes=np.array([quadrature_order**-4]),
                    panel_cuts=np.array([0.,R]),charged_evaluations=5)
    monkeypatch.setattr(tool,"background_cap_metric",metric)
    state=dict(Ip=[0.]*3,Iq=[0.]*3)
    case=dict(cells=3,alpha=.01,snapshots=[dict(x=i,state=dict(cap=state),
        physical=dict(wall=dict(r=.5,T=.1))) for i in range(129)])
    profile=object()
    result=tool.frozen_quadrature(case,profile,tool.Budget(maximum_evaluations=100000,cpu_seconds=30.))
    assert len(calls)==6
    for a,b in zip(calls[::2],calls[1::2]):
        assert a[:4]==b[:4] and (a[-1],b[-1])==(4,8)
    assert result["rows"][0]["difference"]["n_wall"]==pytest.approx(8**-4-4**-4)


def test_create_check_behavior_and_shared_ceiling(tmp_path,monkeypatch):
    monkeypatch.setattr(tool,"ROOT",tmp_path)
    monkeypatch.setattr(tool,"source_bindings",lambda:dict(source="fixed"))
    monkeypatch.setattr(tool, "verify_saved_bindings", lambda expected, root, **kw: expected == tool.source_bindings())
    runtime=tool.runtime_versions()
    def calculate(budget):
        budget.check(70)
        return dict(runtime=runtime,execution=dict(evaluations=70,cpu_seconds=0.))
    def validate(record,budget=None,**kwargs):
        if budget is not None:
            budget.check(20)
        return dict(authentication_evaluations=20)
    monkeypatch.setattr(tool,"calculate",calculate)
    monkeypatch.setattr(tool,"validate_record",validate)
    out=tmp_path/"new"
    result=tool.create(out)
    assert result["creation_execution"]["total_evaluations"]==90
    assert tool.check(out)["authentication_evaluations"]==20
    with pytest.raises(FileExistsError):
        tool.create(out)
    (out/"background.json").write_bytes(b"{}")
    with pytest.raises(ValueError,match="hash/source"):
        tool.check(out)


def test_creation_source_drift_during_authentication_leaves_failure_receipt(tmp_path,monkeypatch):
    monkeypatch.setattr(tool,"ROOT",tmp_path)
    state={"changed":False}
    monkeypatch.setattr(tool,"source_bindings",lambda:dict(source="after" if state["changed"] else "before"))
    monkeypatch.setattr(tool, "verify_saved_bindings", lambda expected, root, **kw: expected == tool.source_bindings())
    monkeypatch.setattr(tool,"calculate",lambda budget:dict(runtime=tool.runtime_versions(),execution=dict(evaluations=0)))
    def validate(*args,**kwargs):
        state["changed"]=True
        return dict(authentication_evaluations=0)
    monkeypatch.setattr(tool,"validate_record",validate)
    with pytest.raises(RuntimeError,match="authentication") as error:
        tool.create(tmp_path/"rejected")
    assert not (tmp_path/"rejected").exists()
    assert json.loads(open(error.value.failure_receipt).read())["status"]=="FAILED_NOT_COMPLETED"


@pytest.mark.parametrize("field,value",(
    ("cpu_seconds",float("nan")),("cpu_seconds",float("inf")),
    ("cpu_seconds",True),("cpu_seconds",600.1),
    ("maximum_evaluations",float("nan")),("maximum_evaluations",True),
    ("maximum_evaluations",100000001),("evaluations",-1),
    ("started",float("nan")),
))
def test_shared_allowance_rejects_disabled_or_expanded_guard(field,value):
    budget=tool.Budget(maximum_evaluations=100000,cpu_seconds=30.)
    setattr(budget,field,value)
    with pytest.raises(ValueError,match="shared background allowance"):
        tool.shared_budget(budget)


def test_shared_allowance_preserves_clock_count_and_exhaustion():
    budget=tool.Budget(maximum_evaluations=10,cpu_seconds=30.)
    budget.check(9)
    started=budget.started
    assert tool.shared_budget(budget) is budget
    assert budget.started==started and budget.evaluations==9
    budget.check(1)
    assert tool.shared_budget(budget).evaluations==10
    budget.evaluations=11
    with pytest.raises(RuntimeError,match="budget exhausted"):
        tool.shared_budget(budget)
