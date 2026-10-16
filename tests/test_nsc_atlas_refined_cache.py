"""Midpoint interface controls; no physical source solve or complete physical cache."""
import copy
import json
import time

import numpy as np
import pytest

import recursive_horizons.nsc_atlas_refined_cache as tool
from recursive_horizons.nsc_atlas_characteristic_trace import SourceSnapshot
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge

GEOMETRY=dict(U0=-.32130105902052203,Ua=-.2992702973503126,
    Ub=-.2747468030926288,Va=-2.565695357068247,Vb=-2.262011613052105,
    order=64,support_V=[-2.489774421064211,-2.3379325490561405],
    support_x=[.42269179059454887,.4631925883576189],x_corner=.4,x_last=.4815042082576364)


@pytest.fixture(scope="module")
def budget():
    shared=tool.Budget(maximum_evaluations=1000000,cpu_seconds=30.)
    yield shared
    shared.check()
    print("refined cache controls",shared.evaluations,time.process_time()-shared.started)


def simple_source(x,alpha):
    x=np.asarray(x)
    return dict(x=x,r=x+alpha*(x-.4)**3,r_x=1+3*alpha*(x-.4)**2,
        c=x**3,c_x=3*x*x,c_xx=6*x,C=np.ones_like(x),log_C=np.zeros_like(x),
        source_constraint=alpha*(x-.4),source_h_xx_defect=x-.4)


def simple_gauge(x):
    x=np.asarray(x);one=np.ones_like(x);zero=np.zeros_like(x)
    result=dict(x=x,W=.2+x-.4,W_x=one,W_xx=zero,C=one,r=.45+(x-.4)/2,
        s=x-.4,s_x=one,x_s=one,Ubar_s=one,T=x-.4,T_x=one,T_s=one,
        L=zero,normalization_defect=zero,inverse_iterations=zero.astype(int),inverse_residual=zero)
    return {key:value.copy() for key,value in result.items()}


@pytest.fixture(scope="module")
def synthetic(budget):
    """Exercise full513+512 orchestration with explicitly manufactured interfaces."""
    patch=pytest.MonkeyPatch()
    constructions=[]
    class Snapshot:
        def __init__(self,geometry,*,alpha,x_nodes,source_rows,budget):
            assert len(x_nodes)==129
            np.testing.assert_array_equal(source_rows["x"],x_nodes)
            constructions.append((alpha,np.array(x_nodes),budget))
            self.alpha,self.budget,self.evaluations=alpha,budget,132
            budget.check(132)
        def sample(self,x):
            self.evaluations+=len(x);self.budget.check(len(x))
            return simple_source(x,self.alpha)
        def summary(self):
            return dict(evaluations=self.evaluations,source_integrations=0)
    class Gauge:
        def __init__(self,geometry,*,order,budget):
            self.order,self.budget,self.evaluations=order,budget,1
            budget.check(1)
        def sample(self,x):
            self.evaluations+=len(x);self.budget.check(len(x))
            row=simple_gauge(x)
            # Integer control is deliberately order-dependent, never physical.
            row["inverse_iterations"][:]=7 if self.order==64 else 8
            return row
        def summary(self):
            return dict(evaluations=self.evaluations,endpoint_integral_defect=0.,
                maximum_inverse_residual=0.,maximum_inverse_iterations=8,quadrature_nodes=0)
    patch.setattr(tool,"SourceSnapshot",Snapshot)
    patch.setattr(tool,"ReferenceWallGauge",Gauge)
    x=np.linspace(.4,GEOMETRY["x_last"],513)
    base=dict(geometry=GEOMETRY,gauge=dict(rows=tool._plain(simple_gauge(x))),
        initial_data={label:dict(alpha=alpha,exact_rows=tool._plain(simple_source(x[::2],alpha)),
            snapshot_rows=tool._plain(simple_source(x,alpha))) for label,alpha in (("zero",0.),("finite",.01))})
    budget.check(3*len(x))
    before=budget.evaluations
    result=tool.make_refined_cache(base,budget)
    assert result["execution"]["evaluations"]==budget.evaluations-before
    yield base,result,constructions
    patch.undo()


def test_all_original_rows_exact_and_midpoints_not_new_uniform_rounding(synthetic,budget):
    base,cache,constructions=synthetic
    x=np.asarray(base["gauge"]["rows"]["x"])
    refined=np.asarray(cache["x_nodes"])
    np.testing.assert_array_equal(refined[::2],x)
    np.testing.assert_array_equal(refined[1::2],x[:-1]+np.diff(x)/2)
    assert len(refined)==1025 and all(item[2] is budget for item in constructions)
    for label in ("zero","finite"):
        for key,value in base["initial_data"][label]["snapshot_rows"].items():
            assert cache["source_rows"][label][key][::2]==value
        np.testing.assert_array_equal(constructions[0 if label=="zero" else 1][1],x[::4])
    for key,value in base["gauge"]["rows"].items():
        assert cache["gauge_rows"][key][::2]==value
    assert cache["execution"]["source_snapshot_points"]==2314
    assert cache["execution"]["source_integrations"]==cache["execution"]["time_evolutions"]==0
    assert "inverse_iterations" not in cache["diagnostics"]["gauge64_minus32_maximum_absolute"]
    json.dumps(cache,allow_nan=False)


def test_reconstruction_deterministic_with_used_budget_and_no_reset(synthetic,budget):
    base,cache,_=synthetic
    before,clock=budget.evaluations,budget.started
    receipt=tool.validate_refined_cache(cache,base,budget)
    assert receipt["validation_evaluations"]==cache["execution"]["evaluations"]
    assert budget.evaluations-before==receipt["validation_evaluations"] and budget.started==clock


def test_midpoint_output_tampering_rejected(synthetic,budget):
    base,cache,_=synthetic
    changed=copy.deepcopy(cache)
    changed["source_rows"]["finite"]["c"][1]+=1e-9
    with pytest.raises(ValueError,match="deterministic reconstruction"):
        tool.validate_refined_cache(changed,base,budget)


@pytest.mark.parametrize("change",("source","source_signed_zero","old_key","exact_key"))
def test_existing_rows_or_keys_cannot_silently_change(synthetic,budget,change):
    base=copy.deepcopy(synthetic[0])
    if change=="source":base["initial_data"]["zero"]["snapshot_rows"]["c"][7]+=1e-10
    elif change=="source_signed_zero":base["initial_data"]["zero"]["snapshot_rows"]["log_C"][7]=-0.
    elif change=="old_key":base["initial_data"]["zero"]["snapshot_rows"]["x"][7]=np.nextafter(base["initial_data"]["zero"]["snapshot_rows"]["x"][7],1.)
    else:base["initial_data"]["zero"]["exact_rows"]["x"][7]=np.nextafter(base["initial_data"]["zero"]["exact_rows"]["x"][7],1.)
    with pytest.raises(ValueError,match="differs"):
        tool.make_refined_cache(base,budget)


def test_gauge_control_rejects_physical_disagreement_but_not_iteration_count(budget):
    x=np.array([.41,.43,.45]);a=simple_gauge(x);b=copy.deepcopy(a)
    b["inverse_iterations"][:]=3
    assert max(tool._compare_gauges(a,b,budget).values())==0
    b["W_xx"][0]=6e-12
    with pytest.raises(ValueError,match="sensitivity"):
        tool._compare_gauges(a,b,budget)
    a["normalization_defect"][0]=6e-12
    with pytest.raises(ValueError,match="normalization"):
        tool._gauge_checks(a,budget)
    a=simple_gauge(x);a["inverse_residual"][1]=4e-13
    with pytest.raises(ValueError,match="inverse"):
        tool._gauge_checks(a,budget)
    a=simple_gauge(x);a["T_x"][0]=-1
    with pytest.raises(ValueError,match="future"):
        tool._gauge_checks(a,budget)


def test_unrepresentable_or_nonuniform_midpoints_are_rejected(budget):
    x=np.linspace(.4,.5,513);x[9]+=1e-7
    with pytest.raises(ValueError,match="uniform"):
        tool._midpoints(x,budget)
    x=1+np.arange(513)*np.spacing(1.)
    with pytest.raises(ValueError,match="strictly inside"):
        tool._midpoints(x,budget)


def pulse(x):
    V=(1-x/.8)*np.exp(x/.8)/GEOMETRY["U0"]
    z=4*(V-(GEOMETRY["Va"]+GEOMETRY["Vb"])/2)/(GEOMETRY["Vb"]-GEOMETRY["Va"])
    zx=-4*x*np.exp(x/.8)/(.64*GEOMETRY["U0"]*(GEOMETRY["Vb"]-GEOMETRY["Va"]))
    zxx=-4*np.exp(x/.8)*(1+x/.8)/(.64*GEOMETRY["U0"]*(GEOMETRY["Vb"]-GEOMETRY["Va"]))
    inside=abs(z)<1
    g=np.where(inside,(1-z*z)**4,0.)
    gz=np.where(inside,-8*z*(1-z*z)**3,0.)
    gzz=np.where(inside,8*(1-z*z)**2*(7*z*z-1),0.)
    return g,gz*zx,gzz*zx*zx+gz*zxx


def test_actual_snapshot_keeps_analytic_scalar_and_polynomial_geometry(budget):
    # Manufactured cubic h; deliberately no claim that it solves the source ODE.
    x=np.linspace(.4,GEOMETRY["x_last"],17);alpha=1e-5;t=x-.4
    h,hx=100*t**3,300*t*t;g,gx,_=pulse(x)
    y,yx=x+alpha*h,1+alpha*hx;F=1-alpha*g*g/(6*y*y)
    r=y*np.sqrt(F);rx=(y*yx-alpha*g*gx/6)/r
    p,px=g/y,(gx-g*yx/y)/y;a=(.8/x-1)/4
    rows=dict(x=x,h=h,h_x=hx,rJ=y,rJ_x=yx,r=r,r_x=rx,C=F,
        log_C=np.log1p(-alpha*g*g/(6*y*y)),log_C_x=-alpha*p*px/(3*F),
        c=g,c_x=gx,r_Ubar=a,c_Ubar=-t*t,log_C_Ubar=-(1/.8+1/x)*a)
    budget.check(17*12)
    snap=SourceSnapshot(GEOMETRY,alpha=alpha,x_nodes=x,source_rows=rows,budget=budget)
    new=(x[:-1]+x[1:])/2;row=snap.sample(new);expected=pulse(new)
    for key,value in zip(("c","c_x","c_xx"),expected):
        np.testing.assert_allclose(row[key],value,rtol=2e-12,atol=1e-10)
    np.testing.assert_allclose(row["h"],100*(new-.4)**3,atol=3e-17)
    np.testing.assert_allclose(row["h_xx"],600*(new-.4),atol=9e-12)
    np.testing.assert_allclose(row["source_constraint"],row["source_constraint_from_curvature_defect"],atol=3e-14)
    assert np.max(abs(row["source_constraint"]))>1e-5
    budget.check(12*len(new))


def test_actual_reference_orders_and_independent_vacuum_clock_identities(budget):
    x=np.array([.411,.439,.475]);rows={}
    for order in (64,32):
        gauge=ReferenceWallGauge(GEOMETRY,order=order,budget=budget)
        rows[order]=gauge.sample(x)
        tool._gauge_checks(rows[order],budget)
    assert max(tool._compare_gauges(rows[64],rows[32],budget).values())<5e-12
    row=rows[64]
    mass=row["r"]*(1+4*row["r_Ubar"]*row["r_x"]/row["C"])/2
    np.testing.assert_allclose(mass,.4,atol=4e-16)
    np.testing.assert_allclose(row["r_x_total"],row["r_Ubar"]*row["W_x"]+row["r_x"],atol=5e-15)
    np.testing.assert_allclose(row["L_x_total"],row["log_C_Ubar"]*row["W_x"]+row["log_C_x"],atol=2e-14)
    np.testing.assert_allclose(row["T_x"],row["beta_minus"]*row["s_x"]/(1-row["r"]**2),atol=2e-14)
    budget.check(4*len(x))


@pytest.mark.parametrize("name,value",(("maximum_evaluations",50000001),("cpu_seconds",301.),
    ("evaluations",-1),("started",float("nan"))))
def test_budget_guards_do_not_expand_allowance(name,value):
    budget=tool.Budget(maximum_evaluations=100,cpu_seconds=30.)
    setattr(budget,name,value)
    with pytest.raises(ValueError):tool._budget(budget)


def test_exhaustion_stops_without_reconstruction(synthetic):
    budget=tool.Budget(maximum_evaluations=1,cpu_seconds=30.)
    with pytest.raises(RuntimeError,match="budget exhausted"):
        tool.make_refined_cache(synthetic[0],budget)


def test_nonfinite_json_rejected_before_any_reconstruction(synthetic,budget):
    cache=copy.deepcopy(synthetic[1]);cache["x_nodes"][1]=float("nan")
    before=budget.evaluations
    with pytest.raises(ValueError,match="JSON-compatible"):
        tool.validate_refined_cache(cache,synthetic[0],budget)
    assert budget.evaluations==before


def test_validation_roundoff_does_not_change_creation_or_copied_rows(synthetic,budget,monkeypatch):
    base,cache,_=synthetic
    original=tool.SourceSnapshot.sample
    def rounded(self,x):
        rows=original(self,x)
        rows["r"]=rows["r"]+1e-14
        return rows
    monkeypatch.setattr(tool.SourceSnapshot,"sample",rounded)
    with pytest.raises(ValueError,match="old source reproduction differs"):
        tool.make_refined_cache(base,budget)
    before=copy.deepcopy(cache)
    receipt=tool.validate_refined_cache(cache,base,budget)
    assert receipt["validation_evaluations"]==cache["execution"]["evaluations"]
    assert cache==before
    changed=copy.deepcopy(cache)
    changed["source_rows"]["zero"]["r"][2]=np.nextafter(changed["source_rows"]["zero"]["r"][2],1.)
    with pytest.raises(ValueError,match="copied source rows"):
        tool.validate_refined_cache(changed,base,budget)


def test_replay_rows_preserve_quiet_source_keys_and_counts(budget):
    x=np.array([.41,.44,.48]); g,gx,gxx=pulse(x)
    expected=dict(x=x,c=g,c_x=gx,c_xx=gxx,r=x.copy())
    changed=copy.deepcopy(expected); changed["r"]+=1e-14
    tool._replay_rows(changed,expected,"control",geometry=GEOMETRY,budget=budget)
    changed["c"][0]=1e-16
    with pytest.raises(ValueError,match="differs: c"):
        tool._replay_rows(changed,expected,"control",geometry=GEOMETRY,budget=budget)
    changed=copy.deepcopy(expected); changed["x"][1]=np.nextafter(x[1],1.)
    with pytest.raises(ValueError,match="differs: x"):
        tool._replay_rows(changed,expected,"control",geometry=GEOMETRY,budget=budget)
    assert not tool._same_replayed({"count":2.0},{"count":2})
    assert not tool._same_replayed({"flag":1},{"flag":True})


def test_copied_row_types_signed_zero_and_inventory_remain_exact(budget):
    original={"x":[.4,.5],"count":[1,2],"zero":[0.,0.]}
    supplied={"x":[.4,.45,.5],"count":[1,1.,2],"zero":[0.,0.,0.]}
    assert tool._same_copied_rows(supplied,original)
    changed=copy.deepcopy(supplied); changed["count"][0]=1.
    assert not tool._same_copied_rows(changed,original)
    changed=copy.deepcopy(supplied); changed["zero"][0]=-0.
    assert not tool._same_copied_rows(changed,original)
    assert not tool._same_copied_rows(dict(supplied,extra=[0,0,0]),original)
    budget.check(24)


def test_source_knots_keep_transverse_values_and_mixed_slope_inputs(budget):
    rows={key:(np.arange(7,dtype=float)+i).tolist() for i,key in enumerate(
        ("x","h","h_x","r","C","c","r_x","c_x","r_Ubar","c_Ubar","log_C_Ubar"))}
    original=copy.deepcopy(rows)
    knots=tool._source_knots({"exact_rows":rows})
    assert set(knots)==set(rows)
    for key,value in rows.items():
        np.testing.assert_array_equal(knots[key],value[::2])
    assert rows==original
    budget.check(4*len(rows))


@pytest.mark.parametrize("alpha",(0.,1e-5))
def test_source_curvature_cancellation_and_summary_use_term_conditioning(alpha,budget):
    # Same affine x, but a few rounding units in the subtracted null label.
    # Recompute Q independently; do not invent a relative error on h_xx-Q.
    x=np.array([.44,.445,.45])
    V=(1-x/.8)*np.exp(x/.8)/GEOMETRY["U0"]
    z=4*(V-(GEOMETRY["Va"]+GEOMETRY["Vb"])/2)/(GEOMETRY["Vb"]-GEOMETRY["Va"])
    zx=-4*x*np.exp(x/.8)/(.64*GEOMETRY["U0"]*(GEOMETRY["Vb"]-GEOMETRY["Va"]))
    zxx=zx*(1/x+1/.8)
    def fields(label):
        g=(1-label*label)**4
        bz=-8*label*(1-label*label)**3
        bzz=8*(1-label*label)**2*(7*label*label-1)
        gx,gxx=bz*zx,bzz*zx*zx+bz*zxx
        Q=(g*gxx-2*gx*gx)/(6*x)+g*gx/(3*x*x)
        return g,gx,gxx,Q
    g,gx,gxx,Q=fields(z)
    expected=dict(x=x,h=np.zeros_like(x),h_x=np.zeros_like(x),h_xx=Q,
        c=g,c_x=gx,c_xx=gxx,source_ode_h_xx=Q,
        source_h_xx_defect=np.zeros_like(x),source_constraint=np.zeros_like(x))
    changed=copy.deepcopy(expected)
    shifted=fields(z+4*np.finfo(float).eps)
    for key,value in zip(("c","c_x","c_xx","source_ode_h_xx"),shifted):
        changed[key]=value
    changed["source_h_xx_defect"]=Q-shifted[-1]
    changed["source_constraint"]=alpha*changed["source_h_xx_defect"]/np.sqrt(
        1-alpha*changed["c"]**2/(6*x*x))
    assert np.max(abs(changed["source_h_xx_defect"]))>2e-13
    tool._replay_rows(changed,expected,"source control",geometry=GEOMETRY,alpha=alpha,budget=budget)
    reconstructed={"source_rows":{label:tool._plain(expected) for label in ("zero","finite")}}
    base={"geometry":GEOMETRY,"initial_data":{
        "zero":{"alpha":0.},"finite":{"alpha":alpha}}}
    limits=tool._source_summary_allowances(reconstructed,base)
    summary={"diagnostics":{"source":{label:{
        "maximum_midpoint_source_constraint":0.,"maximum_midpoint_curvature_defect":0.}
        for label in ("zero","finite")}}}
    replay=copy.deepcopy(summary)
    replay["diagnostics"]["source"]["finite"]["maximum_midpoint_curvature_defect"]=float(
        abs(changed["source_h_xx_defect"][1]))
    assert tool._same_replayed(replay,summary,limits)
    damaged=copy.deepcopy(changed)
    damaged["source_h_xx_defect"]+=1e-6
    with pytest.raises(ValueError,match=r"source_h_xx_defect.*maximum difference=.*allowance="):
        tool._replay_rows(damaged,expected,"source control",geometry=GEOMETRY,alpha=alpha,budget=budget)
    replay["diagnostics"]["source"]["finite"]["maximum_midpoint_curvature_defect"]=1e-6
    assert not tool._same_replayed(replay,summary,limits)
    budget.check(100*x.size)
