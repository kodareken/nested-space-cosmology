"""Bounded synthetic dispatch controls; no retained source or scientific replay."""
import copy

import numpy as np
import pytest

import analyze_nsc_atlas_cap_consistency as tool
from test_analyze_nsc_atlas_cap_consistency import budget


def synthetic_inputs():
    x = np.linspace(.4, .48, 1025)
    R, T, s = .45+.5*(x-.4), x-.4, .8*(x-.4)
    gauge = dict(x=x.tolist(), r=R.tolist(), T=T.tolist(), s=s.tolist(), L=np.zeros(1025).tolist())
    source = dict(x=x.tolist(), r=x.tolist(), log_C=np.zeros(1025).tolist(), c=np.zeros(1025).tolist())
    cases = []
    for label, alpha in (("zero",0.),("finite",.03)):
        snapshots = []
        for j in range(129):
            i = 8*j
            snapshots.append(dict(x=float(x[i]),
                state=dict(wall=dict(r=0.,T=0.,s=0.,L=0.),
                           departure={k:[0.]*63 for k in ("r","L","c")}, cap=dict(Ip=[0.]*64,Iq=[0.]*64)),
                physical=dict(wall=dict(r=float(R[i]),T=float(T[i]),s=float(s[i]),L=0.,c=0.))))
        cases.append(dict(amplitude=label,alpha=alpha,cells=64,steps=512,snapshots=snapshots))
    current = dict(cases=cases,settings=dict(tube={}), profile=dict(labels=np.linspace(0,2,65).tolist(),
        values=[0.]*65,first=[0.]*65,second=[0.]*65))
    cache = dict(x_nodes=x.tolist(),gauge_rows=gauge,source_rows=dict(zero=copy.deepcopy(source),finite=copy.deepcopy(source)))
    return current, cache


def install_synthetic_reconstruction(monkeypatch, budget):
    seen = dict(metric=[],field=[],maps=[])

    def metric(y, radius, state, *, T, profile, alpha, quadrature_order, budget):
        assert len(y) == 65 and len(state["Ip"]) == len(state["Iq"]) == 64
        assert profile.quiet_extension
        r = radius*y
        # Smooth, deliberately non-solution metric data: derivatives must be
        # measured from these values, not replaced by desired field equations.
        n = (.01+.02*T)*r**3 + quadrature_order*1e-8*r**2
        delta = alpha*(.002*r*r+.003*T)
        seen["metric"].append((alpha,quadrature_order,T))
        budget.check(2*len(r))
        return dict(radius=radius,T=T,alpha=alpha,n_faces=n,delta_faces=delta,
            panel_cuts=r, wall=dict(A=1-radius**2-2*alpha*n[-1]/radius,
                delta=float(delta[-1]),f=0.,chi=0.,c=0.,w_minus=.07))

    def fields(metric, radii, *, budget):
        assert radii.tolist() == [metric["radius"]]
        seen["field"].append(metric["T"])
        budget.check(1)
        return dict(f=np.array([0.]),Pi=np.array([.1]),Phi=np.array([.03]))

    def hybrid(snapshot, metric, source, gauge, alpha, tube, budget, *, calls):
        assert source["x"] == gauge["x"] == snapshot["x"]
        assert metric["T"] == snapshot["physical"]["wall"]["T"]
        calls["prefix_calls"] += 1
        calls["hybrid_calls"] += 1
        seen["maps"].append((alpha,snapshot["x"],metric["n_faces"][-1]))
        budget.check(3)
        return dict(w_plus=.2,w_minus=.07,k_minus=-.5,r_s=.6,T_s=1.2,Vdot_plus=1.)

    monkeypatch.setattr(tool,"background_cap_metric",metric)
    monkeypatch.setattr(tool,"background_cap_profile",fields)
    monkeypatch.setattr(tool,"selected_hybrid",hybrid)
    return seen


def test_actual_dispatch_inventory_and_deterministic_used_budget(monkeypatch,budget):
    current,cache = synthetic_inputs()
    before = copy.deepcopy((current,cache))
    seen = install_synthetic_reconstruction(monkeypatch,budget)
    start = budget.evaluations
    result = tool.analyze(current,cache,budget)
    assert result["evaluated_points"] == budget.evaluations-start
    assert result["evaluated_points"] == sum(result["execution"]["charges"].values())
    assert result["execution"]["metric_calls"] == len(seen["metric"]) == 64
    assert len(seen["field"]) == 64
    assert result["execution"]["prefix_calls"] == result["execution"]["hybrid_calls"] == len(seen["maps"]) == 12
    for key in ("profile_fits","stage_calls","rhs_calls","source_solves","gauge_quadratures","evolution_calls"):
        assert result["execution"][key] == 0
    for case in result["cases"]:
        for order in ("4","8"):
            saved = case["orders"][order]
            assert [row["cut"] for row in saved["slices"]] == list(tool.NEIGHBORS)
            assert [row["cut"] for row in saved["targets"]] == list(tool.TARGETS)
            for row in saved["targets"]:
                assert [(v["time_stride"],v["radial_stride"]) for v in row["variants"]] == [(1,1),(1,2),(2,1),(2,2)]
                assert row["variants"][0]["residual"]["e_m"] != 0.
                assert row["variants"][0]["traces"]["hybrid_minus_left"] == pytest.approx(.07)
        assert any(v["differences"]["residual"]["e_r"] != 0 for row in case["quadrature_differences"] for v in row["variants"])
    assert before == (current,cache)
    # The same used allowance must produce identical local counters/results.
    second = tool.analyze(current,cache,budget)
    assert second == result


def test_first_neighbor_key_failure_has_no_reconstruction(monkeypatch,budget):
    current,cache = synthetic_inputs()
    seen = install_synthetic_reconstruction(monkeypatch,budget)
    cache["gauge_rows"]["x"][8*61] += 1e-5
    with pytest.raises(ValueError,match="neighbor gauge") as caught:
        tool.analyze(current,cache,budget)
    info = caught.value.cap_consistency_failure
    assert info["phase"] == dict(name="cap metric",amplitude="zero",order=4,cut=61)
    assert info["counts"] == dict(metric_calls=0,prefix_calls=0,hybrid_calls=0)
    assert seen == dict(metric=[],field=[],maps=[])


def test_selected_adapter_preserves_endpoints_gradient_and_map(monkeypatch,budget):
    x,W,Wx = .42,.2,.25
    eta = np.linspace(0,1,65)
    amplitudes = dict(r=.01,L=-.02,c=.03)
    source = dict(x=x,r=x+.02,log_C=.03,c=.07,r_x=1.,log_C_x=.1,c_x=.2)
    physical = dict(r=.47,L=-.03,c=.1,T=.02,s=.03)
    gauge = dict(x=x,W=W,W_x=Wx,eta_plus=-.5*np.log(Wx),**{k:physical[k] for k in ("r","L","T","s")})
    state = dict(wall={k:0. for k in ("r","L","T","s")},
        departure={k:(a*eta[1:-1]**2).tolist() for k,a in amplitudes.items()})
    # The lapse endpoint is vacuum(-.04)+offset(.03)+departure(-.02).
    metric = dict(radius=physical["r"],T=physical["T"],wall=dict(c=physical["c"],A=.8,delta=.01,w_minus=.1))
    snapshot = dict(x=x,state=state)
    touched = []

    def vacuum(U,at,*,budget):
        assert at == x
        budget.check(len(U)); touched.append("vacuum")
        return dict(r=at+.1*U,L=-.2*U)

    def prefix(U,X,jets,src,*,x,alpha,tube,quadrature_order,budget):
        touched.append("prefix"); budget.check(len(U))
        assert quadrature_order == 2 and alpha == .03
        assert src == source
        for key,jet in (("r","r_Ubar"),("L","log_C_Ubar"),("c","c_Ubar")):
            assert X[key][0] == source["log_C" if key=="L" else key]
            assert X[key][-1] == physical[key]
            np.testing.assert_allclose(jets[jet],2*amplitudes[key]*eta/W,atol=3e-14,rtol=2e-12)
        return dict(r_x=np.full(65,.9),c_x=np.full(65,.2),log_C_x=np.full(65,.1))

    def mapped(state,outgoing,cap,gauge,*,alpha):
        touched.append("map")
        assert state == dict(r=physical["r"],c=physical["c"],eta_plus=-.5*np.log(Wx))
        assert outgoing == dict(r_V=.9,c_V=.2,log_C_V=.1)
        assert cap == dict(A=.8,delta=.01,w_minus=.1)
        assert gauge == dict(C_wall=np.exp(physical["L"]),log_C_s=0.)
        return {k:.7 for k in ("w_plus","Pi","Phi","r_s","T_s","f_s","c_s","k_minus","k_plus",
                              "beta_minus","beta_plus","H","Vdot_plus","Udot_plus")}

    monkeypatch.setattr(tool,"vacuum_slice",vacuum)
    monkeypatch.setattr(tool,"departure_slice",prefix)
    monkeypatch.setattr(tool,"hybrid_wall_boundary",mapped)
    saved = copy.deepcopy((snapshot,metric,source,gauge))
    counts = dict(prefix_calls=0,hybrid_calls=0)
    row = tool.selected_hybrid(snapshot,metric,source,gauge,.03,{},budget,calls=counts)
    assert touched == ["vacuum","prefix","map"]
    assert counts == dict(prefix_calls=1,hybrid_calls=1)
    assert row["w_plus"] == .7 and row["w_minus"] == .1
    assert saved == (snapshot,metric,source,gauge)
