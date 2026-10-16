"""Bounded manufactured comparison controls; no retained physical evaluation."""
import copy
import time

import numpy as np
import pytest

import analyze_nsc_atlas_operator_comparison as tool
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget


@pytest.fixture(scope="module")
def budget():
    shared = Budget(maximum_evaluations=1_000_000, cpu_seconds=30.)
    yield shared
    shared.check()
    print("operator analyzer controls", shared.evaluations,
          "points", time.process_time()-shared.started, "process CPU seconds")


def grid(x, W, eta, radius, marker=0.):
    shape = (len(x)-2, len(eta)-2)
    masks = {k: np.full(shape, k in ("all_interior", "away_from_front", "bottom_origin"))
             for k in ("all_interior", "front_near", "away_from_front", "bottom_origin", "wall_origin", "separator")}
    return dict(x=x[1:-1], W=W[1:-1], eta=eta[1:-1], Ubar=W[1:-1, None]*eta[None, 1:-1],
        radius=radius[1:-1, 1:-1], source_weighted=np.full(len(x)-2, .125), masks=masks,
        arrays={k: np.full(shape, marker) for k in tool.spatial.CHANNELS})


def test_method_change_allows_pure_U_difference_and_keeps_signed_source_baseline(budget):
    x=np.linspace(.4, .5, 5); eta=np.linspace(0, 1, 9); W=.2+x-.4
    old=grid(x,W,eta,np.ones((5,9)),.3); new=copy.deepcopy(old)
    for key in tool.spatial.CHANNELS:
        new["arrays"][key] += .2
    new["masks"]["front_near"][0,0]=True
    got=tool.compare_constraints(new,old,budget)
    for value in got["arrays"].values():
        np.testing.assert_allclose(value,.2,atol=1e-16)
    assert got["comparison_kind"]=="same-grid operator change"
    assert not got["identical_pure_U_samples_required"]
    assert got["same_spatial_diagnostic_widths"] and "same_spatial_stencil" not in got
    assert got["masks"]["front_near"][0,0]
    np.testing.assert_array_equal(new["source_weighted"],.125)
    with pytest.raises(ValueError,match="pure-U"):
        tool.compare_constraints(new,old,budget,time_sampling=True)


def test_same_state_time_sampling_and_coordinate_guards(budget):
    x=np.linspace(.4,.5,9); eta=np.linspace(0,1,9); W=.2+x-.4
    fine=grid(x,W,eta,np.ones((9,9)))
    common=tool.spatial.select_grid(fine,rows=slice(1,None,2))
    got=tool.compare_constraints(common,copy.deepcopy(common),budget,time_sampling=True)
    assert got["identical_pure_U_samples_required"]
    for key in ("source_weighted","eta"):
        bad=copy.deepcopy(common); bad[key][0]+=.01
        with pytest.raises(ValueError,match="common "+key):
            tool.compare_constraints(bad,common,budget)
    with pytest.raises(ValueError,match="boolean"):
        tool.compare_constraints(common,common,budget,time_sampling=1)


def setup_dispatcher(monkeypatch,budget):
    """Full production shapes with cheap repeated history/profile stand-ins.

    Real table extraction, independent wall stencils and comparison arrays run.
    The old grids/estimates are manufactured once, not recalculated by dispatch.
    """
    x1025=np.linspace(.4,.5,1025); x513=x1025[::2]; x=x1025[::8]; eta=np.linspace(0,1,65)
    def source(axis):
        z=np.zeros_like(axis)
        return dict(x=axis,r=axis,r_x=np.ones_like(axis),r_xx=z,log_C=z,log_C_x=z,c=z,c_x=z)
    def gauge(axis):
        return dict(x=axis,W=.2+axis-.4,W_x=np.ones_like(axis),W_xx=np.zeros_like(axis))
    base=dict(geometry=dict(support_x=[.425,.475]),gauge=dict(rows=gauge(x513)),initial_data={})
    cache=dict(x_nodes=x1025,gauge_rows=gauge(x1025),source_rows={})
    record=dict(cases=[])
    old=dict(cases=[],diagnostics=dict(constraint_transport=[],wall_curvature=dict(cases=[])))
    W=gauge(x)["W"]
    for label,alpha,new_shift,old_shift in (("zero",0.,1.,2.),("finite",.001,4.,8.)):
        base["initial_data"][label]=dict(alpha=alpha,snapshot_rows=source(x513))
        cache["source_rows"][label]=source(x1025)
        for target,shift in ((record,new_shift),(old,old_shift)):
            scale=dict(r=.001*shift,L=.0002*shift,c=.003*shift)
            rows=[]
            for t,w in zip(x,W):
                rows.append(dict(x=t,state=dict(departure={k:a*eta[1:-1]**2 for k,a in scale.items()},
                        cap=dict(Ip=np.full(64,7/64),Iq=np.full(64,9/64))),
                    physical=dict(increment={k:a*eta[1:-1]**2 for k,a in scale.items()},
                        wall=dict(r=t+w+scale["r"],L=scale["L"],c=scale["c"],T=.01,s=.02),
                        cap=dict(Ip=np.full(64,(2+shift)/64),Iq=np.full(64,(3+shift)/64)))))
            target["cases"].append(dict(amplitude=label,cells=64,steps=512,alpha=alpha,snapshots=rows,
                final=dict(wall_kinematic_defect={k:0. for k in scale},
                    canonical_wall_flux_defect=0.,conformal_wall_flux_defect=0.)))
        old["diagnostics"]["constraint_transport"].append(dict(grids={
            "native64_"+str(nt):grid(x[::s],W[::s],eta,np.ones((nt,65))) for nt,s in ((129,1),(65,2))}))
        # Independent quadratic analytic wall oracle; source offsets are zero.
        r=x+W+.001*old_shift; c=np.full(129,.003*old_shift)
        rU=1+.002*old_shift/W; rUU=.002*old_shift/W**2
        LU=.0004*old_shift/W; cU=.006*old_shift/W; u=c/r; uU=(cU-u*rU)/r
        terms=dict(curvature=rUU,lapse_product=-LU*rU,scalar=alpha*r*uU*uU/(2*(1+alpha*u*u/6)))
        E=sum(terms.values()); estimates={}
        for p in (3,5):
            for s in (1,2):
                estimates[f"p{p}_s{s}"]=dict(points=p,stride=s,eta_spacing=s/64,
                    first_indices=64-s*np.arange(p),second_indices=64-s*np.arange(4 if p==3 else 5),
                    arrays=dict(E_Ubar=E,weighted=r*E,**terms))
        old["diagnostics"]["wall_curvature"]["cases"].append(dict(amplitude=label,estimates=dict(fine64=estimates)))
    seen=dict(histories=[],grids=[],caps=[],oracle=[])
    def vacuum(w,xx,*,budget):
        budget.check(len(xx))
        return dict(r=xx+w,L=np.zeros_like(xx),r_Ubar=np.ones_like(xx),
            r_Ubar_Ubar=np.zeros_like(xx),log_C_Ubar=np.zeros_like(xx))
    def histories(actual_base,case,actual_budget,stride):
        assert actual_base is base and actual_budget is budget
        assert any(case is c for c in record["cases"])
        assert len(actual_base["gauge"]["rows"]["x"])==513
        seen["histories"].append((case["amplitude"],stride))
        return dict(raw=stride),dict(departed=stride)
    def make_grid(history,xx,ee,ww,radius,src,alpha,fronts,actual_budget):
        assert radius.shape==(len(xx),65) and len(src["x"])==len(xx)
        np.testing.assert_array_equal(src["x"],xx)
        seen["grids"].append((len(xx),len(ee)))
        return grid(xx,ww,ee,radius)
    def compare(left,right,actual_budget,**kwargs):
        actual_budget.check(10)
        for key in ("x","eta","W","Ubar","source_weighted"):
            np.testing.assert_array_equal(left[key],right[key])
        return dict(shape=list(left["Ubar"].shape),same_spatial_stencil=kwargs.get("same_spatial_stencil",False))
    def cap(case,profile,actual_budget):
        seen["caps"].append(case)
        assert np.isclose(sum(case["snapshots"][0]["state"]["cap"]["Iq"]),9.)
        return dict(arrays=dict(upper_label_margin=np.ones(129)))
    def oracle(case,actual_budget):
        seen["oracle"].append(case)
        assert np.isclose(sum(case["snapshots"][0]["state"]["cap"]["Iq"]),4.)
        np.testing.assert_array_equal(case["snapshots"][0]["state"]["cap"]["Iq"],case["snapshots"][0]["physical"]["cap"]["Iq"])
        return dict(total_cap_verified=True)
    monkeypatch.setattr(tool.spatial,"vacuum_slice",vacuum)
    monkeypatch.setattr(tool.inherited,"_histories",histories)
    monkeypatch.setattr(tool.inherited,"compact_history",lambda h:h)
    monkeypatch.setattr(tool.spatial,"constraint_grid",make_grid)
    monkeypatch.setattr(tool.spatial,"compare_constraints",compare)
    monkeypatch.setattr(tool.inherited,"cap_history",cap)
    monkeypatch.setattr(tool.inherited,"cap_pullback",oracle)
    return record,old,base,cache,seen


def test_production_shape_same_grid_dispatch_and_total_cap_oracle(monkeypatch,budget):
    record,old,base,cache,seen=setup_dispatcher(monkeypatch,budget)
    initial=budget.evaluations; started=budget.started
    result=tool.diagnostics(record,old,base,cache,object(),budget)
    assert result["diagnostic_evaluations"]==budget.evaluations-initial and budget.started==started
    assert seen["histories"]==[(a,s) for a in ("zero","finite") for s in (1,2)]
    assert seen["grids"]==[(129,65),(65,65)]*2
    assert len(seen["caps"])==2 and len(seen["oracle"])==1
    expected={"finite_minus_zero64":3.,"method_zero":-1.,"method_finite":-4.,"paired_method":-3.}
    for name,delta in expected.items():
        np.testing.assert_allclose(result[name]["arrays"]["wall"]["c"],.003*delta,atol=2e-18)
        # Same64-cell comparison keeps each physical integral; no adjacent SUM.
        assert np.asarray(result[name]["arrays"]["cap"]["Iq"]).shape==(129,64)
        np.testing.assert_array_equal(result[name]["arrays"]["cap"]["Iq"],delta/64)
    assert result["wall_curvature"]["analytic_vacuum_points"]==129
    assert result["restricted32_diagnostics"]==result["old_evolutions"]==0
    assert not result["source_constraint_is_zero"] and not result["uniform_error_enclosure"]
    for ci in (0,1):
        transport=result["constraint_transport"][ci]
        assert set(transport["grids"])=={"129","65"}
        assert "native64_129" in transport["baseline_grid_paths"]["129"]
        wall=result["wall_curvature"]["cases"][ci]
        assert "fine64" in wall["baseline_estimates_path"]
        assert set(wall["estimates"])=={"p3_s1","p3_s2","p5_s1","p5_s2"}
        # The independent polynomial curvature agrees for both stencil widths.
        for p in (3,5):
            np.testing.assert_allclose(wall["native_minus_wide"][str(p)]["arrays"]["E_Ubar"],0.,atol=2e-12)
    assert np.isclose(sum(record["cases"][0]["snapshots"][0]["state"]["cap"]["Iq"]),9.)


@pytest.mark.parametrize("kind",["gauge","source","amplitude","resolution"])
def test_changed_fixed_inputs_stop_before_history_work(monkeypatch,budget,kind):
    record,old,base,cache,seen=setup_dispatcher(monkeypatch,budget)
    if kind=="gauge": cache["gauge_rows"]["W_x"][8]+=.01
    elif kind=="source": cache["source_rows"]["zero"]["r_xx"][8]+=.01
    elif kind=="amplitude": record["cases"][0]["alpha"]=.01
    else: old["cases"][0]["cells"]=32
    with pytest.raises(ValueError):
        tool.diagnostics(record,old,base,cache,object(),budget)
    assert not seen["histories"]
