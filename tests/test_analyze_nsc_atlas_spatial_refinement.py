"""Independent saved-grid restriction/curvature controls; no physical evolution."""
import copy
import time

import numpy as np
import pytest

import analyze_nsc_atlas_spatial_refinement as tool
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget


@pytest.fixture(scope="module")
def budget():
    value=Budget(maximum_evaluations=1000000,cpu_seconds=30.)
    yield value
    value.check()
    print("spatial diagnostic controls",value.evaluations,time.process_time()-value.started)


def test_exterior_restriction_keeps_even_nodes_and_does_not_touch_cap(budget):
    n=16; eta=np.linspace(0,1,n+1)
    cap={"Ip":np.arange(n,dtype=float),"Iq":np.arange(n,dtype=float)**2}
    case=dict(cells=n,snapshots=[dict(x=.4,state=dict(departure={k:eta[1:-1]**2 for k in tool.FIELDS},cap=cap),
        physical=dict(increment={k:2*eta[1:-1] for k in tool.FIELDS},wall=dict(r=1,L=0,c=.1),cap=cap))])
    restricted=tool.restrict_case(case,budget)
    assert restricted["cells"]==8
    for key in tool.FIELDS:
        np.testing.assert_array_equal(restricted["snapshots"][0]["state"]["departure"][key],eta[2:-1:2]**2)
        np.testing.assert_array_equal(restricted["snapshots"][0]["physical"]["increment"][key],2*eta[2:-1:2])
    assert restricted["snapshots"][0]["state"]["cap"] is cap
    assert len(case["snapshots"][0]["state"]["departure"]["r"])==15
    table=dict(exterior={"c":np.vstack((eta,eta**2))},wall={"c":np.array([1.,1.])},
        cap={"Ip":np.vstack((np.arange(n),np.arange(n)**2)),"Iq":np.ones((2,n))})
    out=tool.restrict_table(table,budget)
    np.testing.assert_array_equal(out["exterior"]["c"],table["exterior"]["c"][:,::2])
    np.testing.assert_array_equal(out["cap"]["Ip"],table["cap"]["Ip"].reshape(2,8,2).sum(axis=2))
    np.testing.assert_array_equal(out["cap"]["Iq"],np.full((2,8),2.))
    np.testing.assert_array_equal(out["cap"]["Ip"].sum(axis=1),table["cap"]["Ip"].sum(axis=1))


@pytest.mark.parametrize("cells,points,stride",[(n,p,s) for n in (32,64) for p in (3,5) for s in (1,2)])
def test_backward_derivatives_have_independent_polynomial_oracle(cells,points,stride,budget):
    eta=np.linspace(0,1,cells+1); degree=2 if points==3 else 4
    values=np.vstack((eta**degree,3+2*eta+eta**2,eta**3))
    result=tool.backward_jets(values,points=points,stride=stride,budget=budget)
    np.testing.assert_allclose(result["first"][:2],[degree,4.],atol=3e-12)
    np.testing.assert_allclose(result["second"],[degree*(degree-1),2.,6.],atol=3e-11)
    assert result["eta_spacing"]==stride/cells
    assert result["second_indices"][-1]==cells-stride*(3 if points==3 else 4)


@pytest.mark.parametrize("points",[3,5])
def test_equal_width_stencils_sample_identical_polynomial_nodes(points,budget):
    coarse=np.linspace(0,1,33); fine=np.linspace(0,1,65)
    fun=lambda e:np.vstack((np.sin(3*e),np.exp(e)))
    a=tool.backward_jets(fun(coarse),points=points,stride=1,budget=budget)
    b=tool.backward_jets(fun(fine),points=points,stride=2,budget=budget)
    np.testing.assert_array_equal(a["first"],b["first"])
    np.testing.assert_array_equal(a["second"],b["second"])


def test_wall_finite_scalar_terms_and_coordinate_factors(budget):
    e=np.linspace(0,1,65); W=np.array([1.,2.,4.]); alpha=.2
    Z=dict(r=np.broadcast_to(.1*e**3,(3,65)),L=np.broadcast_to(.3*e**2,(3,65)),c=np.broadcast_to(.2*e**4,(3,65)))
    vac=dict(r_Ubar=np.array([.1,.2,.3]),r_Ubar_Ubar=np.array([.02,.03,.04]),log_C_Ubar=np.array([.4,.5,.6]))
    wall=dict(r=np.array([2.,3.,4.]),c=np.array([.3,.4,.5]))
    got=tool.wall_estimate(Z,wall,vac,W,alpha,points=5,stride=2,budget=budget)
    rU=vac["r_Ubar"]+.3/W; rUU=vac["r_Ubar_Ubar"]+.6/W**2
    LU=vac["log_C_Ubar"]+.6/W; cU=.8/W
    u=wall["c"]/wall["r"]; uU=(cU-u*rU)/wall["r"]
    expected=rUU-LU*rU+alpha*wall["r"]*uU**2/(2*(1+alpha*u*u/6))
    np.testing.assert_allclose(got["arrays"]["E_Ubar"],expected,atol=3e-13)
    np.testing.assert_allclose(got["arrays"]["weighted"],wall["r"]*expected,atol=1e-12)
    np.testing.assert_array_equal(got["arrays"]["E_Ubar"],sum(got["arrays"][k] for k in ("curvature","lapse_product","scalar")))


def test_wall_source_offset_has_zero_transverse_derivatives(budget):
    nt=5;n=16;x=np.linspace(.4,.5,nt);W=np.linspace(.1,.2,nt)
    vac=tool.vacuum_slice(W,x,budget=budget)
    source=dict(r=x+.02*x*x,log_C=.1*x,c=.3*x*x)
    offsets=dict(r=source["r"]-x,L=source["log_C"],c=source["c"])
    e=np.linspace(0,1,n+1)
    case=dict(cells=n,snapshots=[dict(x=t,state=dict(departure={k:e[1:-1]**2 for k in tool.FIELDS}),
        physical=dict(wall={k:(vac[k][i] if k!="c" else 0)+offsets[k][i]+1 for k in tool.FIELDS})) for i,t in enumerate(x)])
    Z,wall=tool.wall_departures(case,source,vac,budget)
    for key in tool.FIELDS:
        np.testing.assert_allclose(Z[key],np.broadcast_to(e**2,(nt,n+1)),atol=3e-16)
    zero={k:np.zeros((nt,n+1)) for k in tool.FIELDS}
    vacwall=dict(r=vac["r"],c=np.zeros(nt))
    got=tool.wall_estimate(zero,vacwall,vac,W,0.,points=5,stride=2,budget=budget)
    np.testing.assert_array_equal(got["arrays"]["E_Ubar"],np.zeros(nt))


def manufactured_grid(nt,ne,budget):
    x=np.linspace(.4,.5,nt);eta=np.linspace(0,1,ne);W=.2+x-.4
    source=dict(r=2+x,r_x=np.ones(nt),r_xx=np.full(nt,.3),log_C_x=np.full(nt,.1),
        c=x*x,c_x=2*x)
    alpha=.2;Es=tool.source_constraint(source,alpha,budget)
    radius=source["r"][:,None]+W[:,None]*eta[None,:]
    r=radius[1:-1,1:-1];U=W[1:-1,None]*eta[None,1:-1]
    raw=source["r"][1:-1,None]*Es[1:-1,None]/r+.02*U
    history=dict(E_x=raw,E_Ubar=.1*U,mass=np.full_like(U,.4),source_constraint=Es)
    grid=tool.constraint_grid(history,x,eta,W,radius,source,alpha,[.425,.475],budget)
    return grid,history,source,radius


def test_nonzero_source_transport_and_origin_masks(budget):
    grid,_,source,_=manufactured_grid(9,17,budget)
    U=grid["Ubar"];r=grid["radius"]
    np.testing.assert_allclose(grid["arrays"]["residual"],.02*U,atol=3e-17)
    np.testing.assert_allclose(grid["arrays"]["transport"],r*.02*U,atol=1e-16)
    assert np.max(abs(grid["arrays"]["inherited"]))>.1
    masks=grid["masks"]
    np.testing.assert_array_equal(masks["bottom_origin"]|masks["wall_origin"]|masks["separator"],np.ones_like(U,dtype=bool))
    assert not np.any(masks["bottom_origin"]&masks["wall_origin"])
    assert np.any(masks["wall_origin"]) and np.any(masks["separator"])
    np.testing.assert_array_equal(grid["source_arithmetic_difference"],np.zeros(9))


def test_time_comparison_uses_matched_rows_union_fronts_and_same_EU(budget):
    fine,_,_,_=manufactured_grid(9,17,budget);coarse,_,_,_=manufactured_grid(5,17,budget)
    # This manufacture has the same known fields on every common point.
    selected=tool.select_grid(fine,rows=slice(1,None,2))
    out=tool.compare_constraints(coarse,selected,budget,same_spatial_stencil=True)
    np.testing.assert_array_equal(out["arrays"]["E_Ubar"],np.zeros((3,15)))
    np.testing.assert_allclose(out["arrays"]["transport"],0,atol=1e-16)
    np.testing.assert_array_equal(out["masks"]["front_near"],coarse["masks"]["front_near"]|selected["masks"]["front_near"])
    broken=copy.deepcopy(coarse);broken["arrays"]["E_Ubar"][0,0]+=.1
    with pytest.raises(ValueError,match="pure-U"):
        tool.compare_constraints(broken,selected,budget,same_spatial_stencil=True)


def test_native_restricted_difference_is_separate_from_evolution_difference(budget):
    native,_,_,_=manufactured_grid(9,33,budget);old,_,_,_=manufactured_grid(9,17,budget)
    common=tool.select_grid(native,cols=slice(1,None,2));restricted=copy.deepcopy(old)
    # Independent constant errors distinguish same-state derivative width from evolution.
    for key in tool.CHANNELS:
        common["arrays"][key]=common["arrays"][key]+.3
        restricted["arrays"][key]=restricted["arrays"][key]+.1
    stencil=tool.compare_constraints(common,restricted,budget)
    evolution=tool.compare_constraints(restricted,old,budget)
    np.testing.assert_allclose(stencil["arrays"]["E_x"],.2,atol=2e-16)
    np.testing.assert_allclose(evolution["arrays"]["E_x"],.1,atol=2e-16)
    # Stencil partitions are conservatively intersected, not assigned by centres alone.
    np.testing.assert_array_equal(stencil["masks"]["bottom_origin"],common["masks"]["bottom_origin"]&restricted["masks"]["bottom_origin"])
    broken=copy.deepcopy(old);broken["W"][0]+=1e-10
    with pytest.raises(ValueError,match="common W"):
        tool.compare_constraints(restricted,broken,budget)


def test_raw_pure_constraints_recompute_after_state_restriction(budget):
    # r=2+eta^4, L=c=0 yields E_U=12eta^2/W^2; the
    # centred second derivative has known positive2h^2/W^2 truncation.
    x=np.linspace(.4,.5,5); n=16; W=np.ones(5); eta=np.linspace(0,1,n+1)
    from recursive_horizons.nsc_atlas_coupled_diagnostics import coupled_history_diagnostics
    gauge=dict(W=W,W_x=np.full(5,.01),W_xx=np.zeros(5))
    source=dict(r=np.full(5,2.),r_x=np.zeros(5),r_xx=np.zeros(5),log_C=np.zeros(5),
                log_C_x=np.zeros(5),c=np.zeros(5),c_x=np.zeros(5))
    def rows(e):
        return [dict(x=t,state=dict(increment=dict(r=e[1:-1]**4,L=np.zeros(len(e)-2),c=np.zeros(len(e)-2)),
                    wall=dict(r=3.,L=0.,c=0.))) for t in x]
    args=dict(gauge=gauge,source=source,alpha=0.,tube=dict(r=[1.,4.],L=[-1.,1.],c=[-1.,1.]),budget=budget)
    fine=coupled_history_diagnostics(rows(eta),eta=eta,**args)
    e2=eta[::2];coarse=coupled_history_diagnostics(rows(e2),eta=e2,**args)
    expected=np.broadcast_to(12*e2[None,1:-1]**2+2/n**2,(3,len(e2)-2))
    np.testing.assert_allclose(fine["E_Ubar"][:,1::2],expected,atol=2e-13)
    np.testing.assert_allclose(coarse["E_Ubar"]-fine["E_Ubar"][:,1::2],6/n**2,atol=2e-13)


@pytest.mark.parametrize("points,stride",[(4,1),(3,0),(True,1),(5,3)])
def test_stencil_guard(points,stride,budget):
    with pytest.raises(ValueError,match="stencil"):
        tool.backward_jets(np.zeros((3,33)),points=points,stride=stride,budget=budget)


def test_changed_source_curvature_is_rejected(budget):
    grid,history,source,radius=manufactured_grid(5,9,budget)
    bad=dict(source);bad["r_xx"]=source["r_xx"]+.01
    with pytest.raises(ValueError,match="source constraint"):
        tool.constraint_grid(history,np.linspace(.4,.5,5),np.linspace(0,1,9),np.linspace(.2,.3,5),
            radius,bad,.2,[.425,.475],budget)


def test_caller_budget_and_nonfinite_guard(budget):
    before=budget.evaluations
    tool.backward_jets(np.zeros((3,33)),points=3,stride=1,budget=budget)
    assert budget.evaluations-before==3*33+3*7
    with pytest.raises(ValueError,match="finite"):
        tool.backward_jets(np.full((3,33),np.nan),points=3,stride=1,budget=budget)
    tiny=Budget(maximum_evaluations=1,cpu_seconds=30.)
    with pytest.raises(RuntimeError):
        tool.backward_jets(np.zeros((3,33)),points=3,stride=1,budget=tiny)


def test_real_shape_dispatcher_reuses_old_histories_and_keeps_cap_views(budget,monkeypatch):
    # The numerical helpers above have independent oracles. This test exercises
    # production129x64/32 routing while replacing repeated summary/profile work.
    x1025=np.linspace(.4,.5,1025);x513=x1025[::2];x=x513[::4]
    def source(axis):
        z=np.zeros_like(axis)
        return dict(x=axis,r=axis,r_x=np.ones_like(axis),r_xx=z,log_C=z,log_C_x=z,c=z,c_x=z)
    def gauge(axis):
        return dict(x=axis,W=.2+axis-.4,W_x=np.ones_like(axis),W_xx=np.zeros_like(axis))
    base=dict(geometry=dict(support_x=[.425,.475]),gauge=dict(rows=gauge(x513)),initial_data={})
    parent=dict(cache=dict(x_nodes=x1025,gauge_rows=gauge(x1025),source_rows={}),cases=[],diagnostics=dict(histories=[]))
    record=dict(cases=[])
    for label,alpha in (("zero",0.),("finite",.001)):
        base["initial_data"][label]=dict(snapshot_rows=source(x513))
        parent["cache"]["source_rows"][label]=source(x1025)
        for n,target in ((32,parent),(64,record)):
            # Finite signal differs by1/n, making the paired spatial sign testable.
            shift=0. if label=="zero" else 1.+1/n
            rows=[]
            for t in x:
                rawcap=dict(Ip=np.full(n,7/n),Iq=np.full(n,9/n))
                physicalcap=dict(Ip=np.full(n,(2+shift)/n),Iq=np.full(n,(3+shift)/n))
                rows.append(dict(x=t,state=dict(departure={k:np.zeros(n-1) for k in tool.FIELDS},cap=rawcap),
                    physical=dict(increment={k:np.zeros(n-1) for k in tool.FIELDS},
                        wall=dict(r=.45+shift*.001,L=0.,c=shift,T=.01,s=.02),cap=physicalcap)))
            target["cases"].append(dict(amplitude=label,cells=n,steps=512,alpha=alpha,snapshots=rows,
                final=dict(wall_kinematic_defect={k:0. for k in tool.FIELDS},
                    canonical_wall_flux_defect=0.,conformal_wall_flux_defect=0.)))
        parent["diagnostics"]["histories"].append({"departure129":dict(old_marker=129),"departure65":dict(old_marker=65)})
    histories_seen=[]; grids_seen=[]; cap_seen=[]; oracle_seen=[]
    def histories(actual_base,case,actual_budget,stride=1):
        assert actual_base is base and actual_budget is budget
        histories_seen.append((case["amplitude"],case["cells"],stride))
        assert len(case["snapshots"][0]["state"]["departure"]["r"])==case["cells"]-1
        assert len(case["snapshots"][0]["physical"]["increment"]["r"])==case["cells"]-1
        assert len(case["snapshots"][0]["state"]["cap"]["Iq"])==64
        return dict(raw_marker=stride),dict(new_marker=stride)
    def grid(history,xx,eta,W,radius,src,alpha,fronts,actual_budget):
        nt,ne=len(xx),len(eta);assert radius.shape==(nt,ne)
        assert len(src["x"])==nt;np.testing.assert_array_equal(src["x"],xx)
        grids_seen.append((nt,ne,"old_marker" in history))
        shape=(nt-2,ne-2);U=W[1:-1,None]*eta[None,1:-1]
        masks={k:np.ones(shape,dtype=bool) for k in ("all_interior","front_near","away_from_front","bottom_origin","wall_origin","separator")}
        return dict(x=xx[1:-1],eta=eta[1:-1],W=W[1:-1],Ubar=U,radius=radius[1:-1,1:-1],
            source_weighted=np.zeros(nt-2),arrays={k:np.zeros(shape) for k in tool.CHANNELS},masks=masks)
    def compare(left,right,actual_budget,**kwargs):
        for key in ("x","eta","W","Ubar","source_weighted"):
            np.testing.assert_array_equal(left[key],right[key])
        return dict(shape=list(left["Ubar"].shape),same_spatial_stencil=kwargs.get("same_spatial_stencil",False))
    def cap(case,profile,actual_budget):
        cap_seen.append(case)
        for row in case["snapshots"]:
            assert np.isclose(np.sum(row["state"]["cap"]["Iq"]),9.)
        return dict(arrays=dict(upper_label_margin=np.ones(129)))
    def oracle(case,actual_budget):
        oracle_seen.append(case)
        for row in case["snapshots"]:
            np.testing.assert_array_equal(row["state"]["cap"]["Iq"],row["physical"]["cap"]["Iq"])
            assert np.isclose(np.sum(row["state"]["cap"]["Iq"]),3.)
        return dict(oracle_checked=True)
    monkeypatch.setattr(tool.inherited,"_histories",histories)
    monkeypatch.setattr(tool.inherited,"compact_history",lambda h:h)
    monkeypatch.setattr(tool,"constraint_grid",grid)
    monkeypatch.setattr(tool,"compare_constraints",compare)
    monkeypatch.setattr(tool.inherited,"cap_history",cap)
    monkeypatch.setattr(tool.inherited,"cap_pullback",oracle)
    result=tool.diagnostics(record,parent,base,object(),budget)
    assert histories_seen==[(a,n,s) for a in ("zero","finite") for n in (64,32) for s in (1,2)]
    assert grids_seen==[(nt,ne,old) for _ in (0,1) for nt,ne,old in
        ((129,65,False),(65,65,False),(129,33,False),(65,33,False),(129,33,True),(65,33,True))]
    assert len(cap_seen)==2 and len(oracle_seen)==1
    delta=1/64-1/32
    np.testing.assert_allclose(result["paired_spatial"]["arrays"]["wall"]["c"],delta,atol=0)
    np.testing.assert_allclose(result["paired_spatial"]["arrays"]["cap"]["Iq"],delta/32,atol=0)
    assert result["wall_curvature"]["analytic_vacuum_points"]==129
    assert result["wall_curvature"]["birth_transport_evaluations"]==0
    # The production dispatcher stops at a changed reference key before helpers.
    bad=copy.deepcopy(parent);bad["cache"]["gauge_rows"]["W_x"][8]+=.01
    with pytest.raises(ValueError,match="gauge W_x"):
        tool.diagnostics(record,bad,base,object(),budget)
