"""Bounded synthetic20-metric/four-RHS routing; no retained physical input."""
import copy
import math

import numpy as np
import pytest

import analyze_nsc_atlas_cap_mass_direction as tool
from recursive_horizons.nsc_atlas_reduced_cap import _wall_conversion
from test_analyze_nsc_atlas_cap_mass_direction import budget


def manufactured_metric(y,R,state,*,T,profile,alpha,quadrature_order,budget):
    """A non-solution functional with geometry and both raw-state inputs."""
    Ip,Iq=(np.asarray(state[k]) for k in ("Ip","Iq"))
    P,Q=float(np.sum(Ip)),float(np.sum(Iq))
    N=R*R*T+P*Q+T**3+quadrature_order*1e-5*R
    r=R*y
    nf=N*y**3; df=alpha*.001*y*y
    A=1-r*r
    A[1:]-=2*alpha*nf[1:]/r[1:]
    chi=float(np.cumsum(Iq)[-1]); f=chi/R
    wall=dict(A=float(A[-1]),delta=float(df[-1]),f=f,chi=chi,c=R*_wall_conversion(f,alpha))
    budget.check(5*len(y)+2*len(Ip))
    return dict(radius=R,T=T,alpha=alpha,state=copy.deepcopy(state),quadrature_order=quadrature_order,
        n_faces=nf,delta_faces=df,A_faces=A,A_nodes=A[1:,None],B_faces=A*np.exp(df),
        B_nodes=(A*np.exp(df))[1:,None],r_faces=r,panel_cuts=r,wall=wall)


def manufactured_fields(metric,radii,*,budget):
    budget.check(len(radii))
    return dict(Pi=np.full(len(radii),.1+float(np.sum(metric["state"]["Ip"]))),
                Phi=np.full(len(radii),.03+metric["T"]))


def inputs(budget):
    x=np.linspace(.4,.48,1025)
    R,T=.49,.08
    gauge={k:np.full(1025,v).tolist() for k,v in dict(r=R,T=T,L=0.,s=.07).items()}
    cases=[]; audits=[]
    for label,alpha in (("zero",0.),("finite",.03)):
        state=dict(Ip=np.linspace(.0001,.0002,64).tolist(),Iq=np.linspace(-.0002,.0003,64).tolist())
        physical=dict(r=R,T=T,L=0.,s=.07)
        chi=float(np.cumsum(state["Iq"])[-1]); physical["c"]=R*_wall_conversion(chi/R,alpha)
        snapshots=[dict(x=float(x[8*i]),state=dict(cap=copy.deepcopy(state),
            wall={k:0. for k in ("r","T","L","s")}),physical=dict(wall=copy.deepcopy(physical))) for i in range(129)]
        cases.append(dict(amplitude=label,alpha=alpha,cells=64,steps=512,snapshots=snapshots))
        orders={}
        for order in (4,8):
            metric=manufactured_metric(np.linspace(0,1,65),R,state,T=T,profile=None,
                alpha=alpha,quadrature_order=order,budget=budget)
            field=manufactured_fields(metric,np.array([R]),budget=budget)
            wall=dict(R=R,T=T,alpha=alpha,n=float(metric["n_faces"][-1]),**metric["wall"],
                Pi=float(field["Pi"][0]),Phi=float(field["Phi"][0]))
            Pi,Phi=wall["Pi"],wall["Phi"]
            pdiff=float(tool._potential(np.array([wall["f"]]),alpha)[1][0])
            J=R*R*wall["A"]*wall["A"]*math.exp(wall["delta"])*Pi*Phi/2
            rho=R*R*(pdiff/2+wall["A"]*(Pi*Pi+Phi*Phi)/4)
            variants=[]
            for ts in (1,2):
                for rs in (1,2):
                    nt=.07+ts*.01; nr=.6+rs*.02
                    variants.append(dict(time_stride=ts,radial_stride=rs,
                        residual=dict(n_T=nt,n_r=nr,normalized_mass_flux=J,normalized_radial_density=rho,e_m=nt-J),
                        curve=dict(R_x=.12,T_x=.8)))
            orders[str(order)]=dict(slices=[dict(cut=128,x=float(x[-1]),wall=wall,
                n_faces=metric["n_faces"].tolist(),delta_faces=metric["delta_faces"].tolist())],
                targets=[dict(cut=128,variants=variants,
                    hybrid=dict(r_s=.3+order*.001,T_s=2.4,Vdot_plus=2.,w_plus=.25+alpha))])
        audits.append(dict(amplitude=label,alpha=alpha,orders=orders))
    current=dict(cases=cases,profile=dict(labels=np.linspace(0,2,65).tolist(),values=[0.]*65,
        first=[0.]*65,second=[0.]*65))
    audit=dict(analysis=dict(schema="nsc-atlas-cap-consistency-analysis-v1",cases=audits))
    cache=dict(x_nodes=x.tolist(),gauge_rows=gauge)
    # Charge the manufactured retained row inventory; no scientific solver.
    budget.check(2*129*8+1025+4*65)
    return current,audit,cache


def install(monkeypatch):
    seen=dict(metric=[],rhs=[])
    def metric(*args,**kwargs):
        seen["metric"].append(dict(R=args[1],T=kwargs["T"],alpha=kwargs["alpha"],order=kwargs["quadrature_order"],
            state=copy.deepcopy(args[2])))
        return manufactured_metric(*args,**kwargs)
    def rhs(metric,*,T_x,radius_x,incoming_w_plus,budget):
        seen["rhs"].append((metric["alpha"],metric["quadrature_order"],T_x,radius_x,incoming_w_plus))
        n=len(metric["state"]["Ip"])
        budget.check(2*n+5)
        return dict(rates=dict(Ip=np.full(n,.001/n),Iq=np.full(n,-.002/n)),
            wall_chi_x=.1,wall_f_x=.2,wall_c_x=.3,baseline_wall_chi_x=.4,constructed_wall_flux_defect=0.)
    monkeypatch.setattr(tool,"background_cap_metric",metric)
    monkeypatch.setattr(tool,"background_cap_profile",manufactured_fields)
    monkeypatch.setattr(tool,"background_cap_rhs",rhs)
    return seen


def test_fixed_inventory_backward_routing_and_immutable_used_budget(monkeypatch,budget):
    current,audit,cache=inputs(budget)
    original=copy.deepcopy((current,audit,cache))
    seen=install(monkeypatch)
    start=budget.evaluations
    result=tool.analyze(current,audit,cache,budget)
    assert result["evaluated_points"]==budget.evaluations-start==sum(result["execution"]["charges"].values())
    assert result["execution"]["metric_calls"]==len(seen["metric"])==20
    assert result["execution"]["rhs_calls"]==len(seen["rhs"])==4
    for key in ("profile_fits","stage_calls","prefix_calls","hybrid_calls","source_solves","gauge_quadratures","evolution_calls"):
        assert result["execution"][key]==0
    h=result["h"]
    for i,base in enumerate(seen["metric"][::5]):
        alpha,order,Tx,Rx,wp=seen["rhs"][i]
        assert alpha==base["alpha"] and order==base["order"]
        assert Tx==1.2 and Rx==(.3+order*.001)/2
        for offset,row in zip(np.array([0.,-2.,-1.,-.5,-.25])*h,seen["metric"][5*i:5*i+5]):
            assert row["R"]==base["R"]+offset*Rx
            assert row["T"]==base["T"]+offset*Tx
            assert row["T"]<=base["T"]
            for k,d in (("Ip",.001/64),("Iq",-.002/64)):
                np.testing.assert_array_equal(row["state"][k],np.array(base["state"][k])+offset*d)
    for case in result["cases"]:
        for block in case["orders"].values():
            assert len(block["probes"])==5 and len(block["widths"])==3
            assert block["normalized_flux_imposed_minus_left"]!=0
            for width in block["widths"]:
                assert [(v["time_stride"],v["radial_stride"]) for v in width["radial"]]==[(1,1),(1,2),(2,1),(2,2)]
                for row in width["radial"]:
                    left=row["left"]
                    # Rearranged floating-point identities contain several
                    # rounded products/subtractions. Allow their scaled
                    # arithmetic error, not a sub-epsilon absolute cutoff.
                    scale=max(1.,abs(left["D_n"]),abs(left["B_mass"]),
                              abs(left["e_m_dir"]),abs(left["e_hist_hybrid"]),
                              abs(left["history_D_n"]))
                    rounding=16*np.finfo(float).eps*scale
                    assert left["identity_defect"]==pytest.approx(0.,rel=0,abs=rounding)
                    assert left["derivative_bridge_defect"]==pytest.approx(0.,rel=0,abs=rounding)
                    assert row["left"]["history_rate_bridge"]!=0
    assert (current,audit,cache)==original
    # Determinism on the same already-used allowance includes local counts.
    assert tool.analyze(current,audit,cache,budget)==result


def test_first_failure_retains_probe_phase_and_attempted_counts(monkeypatch,budget):
    current,audit,cache=inputs(budget)
    seen=install(monkeypatch)
    original=tool.background_cap_metric
    def fail_at_probe(*args,**kwargs):
        if kwargs["T"]<.08:
            budget.check(7)
            raise ValueError("manufactured profile coverage failure")
        return original(*args,**kwargs)
    monkeypatch.setattr(tool,"background_cap_metric",fail_at_probe)
    with pytest.raises(ValueError,match="profile coverage") as error:
        tool.analyze(current,audit,cache,budget)
    info=error.value.cap_mass_direction_failure
    assert info["counts"]==dict(metric_calls=2,rhs_calls=1)
    assert info["phase"]["name"]=="backward metric probe"
    assert info["phase"]["offset"]<0 and info["phase"]["amplitude"]=="zero" and info["phase"]["order"]==4
    assert len(seen["metric"])==len(seen["rhs"])==1
