"""One-sided terminal cap mass-functional derivative on frozen directions.

The caller authenticates inputs and owns the cumulative allowance.  Four
base RHS directions are evaluated; probes reconstruct only the radial
functional, including its moving background, panels and limiter branches.
No full stage, source, boundary map or trajectory is called here.
"""
from __future__ import annotations

from collections.abc import Mapping
import math
import time

import numpy as np

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_background_cap import (
    background_cap_metric, background_cap_profile, background_cap_rhs,
)
from recursive_horizons.nsc_atlas_polar_cap import _potential
from analyze_nsc_atlas_cap_consistency import _plain, _scalar, _array, _numbers

SCHEMA = "nsc-atlas-cap-mass-direction-analysis-v1"
ORDERS = (4, 8)
OFFSETS = (0., -2., -1., -.5, -.25)
SCOPE = ("terminal one-sided mass-functional derivative along four frozen cap RHS directions; "
         "width, quadrature and branch sensitivities, not a trajectory or remainder bound")


def _budget(budget):
    if (not isinstance(budget, Budget) or type(budget.maximum_evaluations) is not int
            or not 0 < budget.maximum_evaluations <= 4_000_000
            or type(budget.evaluations) is not int or budget.evaluations < 0
            or isinstance(budget.cpu_seconds, (bool, np.bool_))
            or not math.isfinite(budget.cpu_seconds) or not 0 < budget.cpu_seconds <= 60
            or isinstance(budget.started, (bool, np.bool_))
            or not math.isfinite(budget.started) or not 0 <= budget.started <= time.process_time()):
        raise ValueError("one live cap-mass-direction budget is required")
    budget.check()
    return budget


def stored_direction(hybrid, budget):
    """Convert the already authenticated proper-time rates to original x."""
    _budget(budget)
    h = _numbers(hybrid, ("r_s", "T_s", "Vdot_plus", "w_plus"), "hybrid")
    if h["T_s"] <= 0 or h["Vdot_plus"] <= 0:
        raise ValueError("stored clock orientation must be positive")
    out = dict(R_x=h["r_s"]/h["Vdot_plus"], T_x=h["T_s"]/h["Vdot_plus"], w_plus=h["w_plus"])
    if not all(math.isfinite(v) for v in out.values()): raise ValueError("nonfinite stored direction")
    budget.check(4)
    return out


def _cap_state(state, budget):
    if not isinstance(state, Mapping) or set(state) != {"Ip", "Iq"}:
        raise ValueError("cap state must contain Ip and Iq")
    raw = np.asarray(state["Ip"])
    if raw.ndim != 1 or not 3 <= len(raw) <= 128:
        raise ValueError("bounded cap integral vectors required")
    return {k: _array(state[k], raw.shape, k, budget).copy() for k in ("Ip", "Iq")}


def linear_probe(R, T, state, rates, offset, budget):
    """Raw integral departures and physical geometry on one frozen line."""
    _budget(budget)
    R, T, offset = (_scalar(v, k) for v, k in ((R,"R"),(T,"T"),(offset,"offset")))
    if not 0 < R < 1 or offset > 0: raise ValueError("only backward probes with 0<R<1 are admitted")
    fields = _cap_state(state, budget)
    d = _numbers(rates, ("R", "T"), "rates")
    if d["T"] <= 0: raise ValueError("positive T direction required")
    slopes = {k: _array(rates[k], fields[k].shape, "rate."+k, budget) for k in fields}
    result = dict(R=R+offset*d["R"], T=T+offset*d["T"],
                  state={k: fields[k]+offset*slopes[k] for k in fields})
    budget.check(2+2*len(fields["Ip"]))
    if (not 0 < result["R"] < 1 or not math.isfinite(result["T"])
            or any(not np.all(np.isfinite(v)) for v in result["state"].values())):
        raise ValueError("probe leaves the finite cap domain")
    return result


def backward_derivatives(values, h, budget):
    """Three backward-three-point widths from offsets0,-2h,-h,-h/2,-h/4."""
    _budget(budget)
    n = _array(values, (5,), "functional values", budget)
    h = _scalar(h, "h")
    if h <= 0: raise ValueError("positive probe width required")
    result = []
    for width, near, far in ((h,2,1),(h/2,3,2),(h/4,4,3)):
        derivative = float((3*n[0]-4*n[near]+n[far])/(2*width))
        if not math.isfinite(derivative): raise ValueError("nonfinite directional derivative")
        result.append(dict(width=width,D_n=derivative))
    budget.check(9)
    return result


def balance_quantities(D_n, *, R_x, T_x, n_r, rho_n, J, history_n_T,
                       history_R_x, history_T_x, budget):
    """Global/radial balance and the explicit same-hybrid-rate history bridge."""
    _budget(budget)
    vals = {k:_scalar(v,k) for k,v in dict(D_n=D_n,R_x=R_x,T_x=T_x,n_r=n_r,
        rho_n=rho_n,J=J,history_n_T=history_n_T,history_R_x=history_R_x,history_T_x=history_T_x).items()}
    if vals["T_x"] <= 0 or vals["history_T_x"] <= 0:
        raise ValueError("positive clock rates required")
    em = (D_n-R_x*n_r)/T_x-J
    global_balance = D_n-R_x*rho_n-T_x*J
    er = n_r-rho_n
    history_D = history_T_x*history_n_T+history_R_x*n_r
    original = history_n_T-J
    same = (history_D-R_x*n_r)/T_x-J
    difference = em-same
    result = dict(D_n=D_n,e_m_dir=em,B_mass=global_balance,e_r=er,
        identity_defect=global_balance-T_x*em-R_x*er,history_D_n=history_D,
        e_hist=original,e_hist_hybrid=same,history_rate_bridge=same-original,
        direction_minus_history_same_rates=difference,
        derivative_bridge_defect=difference-(D_n-history_D)/T_x,
        R_x_rate_difference=R_x-history_R_x,T_x_rate_difference=T_x-history_T_x)
    if not all(math.isfinite(v) for v in result.values()): raise ValueError("nonfinite mass balance")
    budget.check(30)
    return result


def limiter_signature(R, state, budget):
    """Exact branch metadata for the frozen p±q minmod reconstruction.

    Cell0 uses the regular polynomial moment override; its formal endpoint
    secant is retained but is never counted as an effective limiter branch.
    Ties are exact floating comparisons, not a differentiability test.
    """
    _budget(budget)
    R = _scalar(R,"R")
    if not 0 < R < 1: raise ValueError("limiter requires0<R<1")
    fields = _cap_state(state,budget)
    n = len(fields["Ip"])
    # Match the frozen reconstruction's physical face/midpoint arithmetic.
    r = R*np.linspace(0,1,n+1)
    dr = np.diff(r); mid = (r[:-1]+r[1:])/2
    p,q = fields["Ip"]/dr,fields["Iq"]/dr
    result = {}
    for key,value in (("plus",p+q),("minus",p-q)):
        sec = np.diff(value)/np.diff(mid)
        left,right = sec[:-1],sec[1:]
        same = left*right > 0
        tie = same & (abs(left) == abs(right))
        selection = ["axis_override"]+["zero" if not ok else "tie" if eq else
            "left" if abs(a)<abs(b) else "right"
            for a,b,ok,eq in zip(left,right,same,tie)]+["endpoint"]
        result[key] = dict(secants=sec,secant_signs=np.sign(sec).astype(int),selection=selection,
            zero_secant_indices=np.flatnonzero(sec==0),
            opposite_or_zero_indices=np.flatnonzero(~same)+1,
            equal_magnitude_indices=np.flatnonzero(tie)+1,
            min_abs_secant=float(np.min(abs(sec))),min_abs_product=float(np.min(abs(left*right))),
            min_selection_gap=float(np.min(abs(abs(left)-abs(right)))))
        budget.check(6*n)
        if not np.all(np.isfinite(sec)) or not math.isfinite(result[key]["min_abs_product"]):
            raise ValueError("nonfinite limiter metadata")
    return _plain(result)


def limiter_changes(base, probe, budget):
    _budget(budget)
    out = {}
    for key in ("plus","minus"):
        a,b = base[key],probe[key]
        if len(a["selection"]) != len(b["selection"]): raise ValueError("limiter shapes differ")
        changes = [i for i,(x,y) in enumerate(zip(a["selection"],b["selection"])) if x!=y and i>0]
        signs = np.flatnonzero(np.asarray(a["secant_signs"])!=np.asarray(b["secant_signs"]))
        ties = sorted(set(a["equal_magnitude_indices"]) ^ set(b["equal_magnitude_indices"]))
        out[key] = dict(selection_changed_cells=changes,secant_sign_changed_indices=signs,tie_changed_cells=ties)
        budget.check(3*len(a["selection"]))
    return _plain(out)


def _metric_row(metric, profile, state, budget):
    R,T,alpha = (float(metric[k]) for k in ("radius","T","alpha"))
    left = background_cap_profile(metric,np.array([R]),budget=budget)
    wall = dict(R=R,T=T,alpha=alpha,n=float(metric["n_faces"][-1]),
        **{k:float(metric["wall"][k]) for k in ("A","delta","f","chi","c")},
        Pi=float(left["Pi"][0]),Phi=float(left["Phi"][0]))
    r = np.asarray(metric["r_faces"]); panels=np.asarray(metric["panel_cuts"])
    xi=float(np.arctanh(R)); qlo,qhi=T-xi,T+xi
    upper=float(profile.domain["last_label"]-qhi)
    lower=float(qlo-profile.domain["first_label"])
    if upper < 0 or (lower < 0 and not profile.quiet_extension):
        raise ValueError("probe leaves profile coverage")
    # Panels are already the union of cell faces and actual profile splits.
    # No second metric/profile quadrature is performed for this metadata.
    profile_only = int(np.count_nonzero(~np.isin(panels,r)))
    report = dict(wall=wall,panels=len(panels)-1,radial_cells=len(r)-1,
        profile_only_splits=profile_only,domain=dict(radius_lower=R,radius_upper=1-R,
            q_minus=qlo,q_plus=qhi,upper_profile_margin=upper,
            lower_profile_margin=lower,quiet_lower_extension_used=lower<0,
            A_min=float(min(np.min(metric["A_faces"]),np.min(metric["A_nodes"]))),
            B_min=float(min(np.min(metric["B_faces"]),np.min(metric["B_nodes"])))),
        limiter=limiter_signature(R,state,budget))
    budget.check(len(panels)+len(r)+np.size(metric["A_nodes"])+20)
    if not all(math.isfinite(v) for v in wall.values()): raise ValueError("nonfinite wall reconstruction")
    return report


def _saved_terminal(cap_case, order):
    block=cap_case["orders"][str(order)]
    rows=[v for v in block["slices"] if v["cut"]==128]
    targets=[v for v in block["targets"] if v["cut"]==128]
    if len(rows)!=1 or len(targets)!=1: raise ValueError("unique saved terminal cap row required")
    target=targets[0]
    if [(v["time_stride"],v["radial_stride"]) for v in target["variants"]] != [(1,1),(1,2),(2,1),(2,2)]:
        raise ValueError("saved derivative inventory differs")
    return rows[0],target


def _width_difference(a,b,budget):
    # Deliberately compare only continuous measurements, never branch labels.
    out=dict(D_n=a["D_n"]-b["D_n"],B_mass_left=a["B_mass_left"]-b["B_mass_left"],
        B_mass_imposed=a["B_mass_imposed"]-b["B_mass_imposed"],
        radial=[dict(time_stride=x["time_stride"],radial_stride=x["radial_stride"],
            e_m_dir=x["left"]["e_m_dir"]-y["left"]["e_m_dir"],
            direction_minus_history_same_rates=x["left"]["direction_minus_history_same_rates"]-
                y["left"]["direction_minus_history_same_rates"])
            for x,y in zip(a["radial"],b["radial"])])
    budget.check(3+4*len(out["radial"]))
    return out


def analyze(current, cap_audit, cache, budget):
    """Exactly20 metrics/four RHS calls; local receipts on a shared budget."""
    _budget(budget); initial=budget.evaluations
    phase=dict(name="input declarations")
    counts=dict(metric_calls=0,rhs_calls=0)
    charges=dict(input_and_profile=0,metrics=0,rhs=0,probe_analysis=0)
    try:
        expected=[("zero",64,512),("finite",64,512)]
        if [(v["amplitude"],v["cells"],v["steps"]) for v in current["cases"]]!=expected:
            raise ValueError("require fixed zero/finite64/512 cases")
        old=cap_audit["analysis"]
        if old["schema"]!="nsc-atlas-cap-consistency-analysis-v1": raise ValueError("cap audit schema differs")
        if [v["amplitude"] for v in old["cases"]]!=["zero","finite"]: raise ValueError("cap cases differ")
        x=_array(cache["x_nodes"],(1025,),"cache x",budget)
        if np.any(np.diff(x)<=0) or not np.allclose(np.diff(x),(x[-1]-x[0])/1024,rtol=0,atol=8*np.finfo(float).eps):
            raise ValueError("uniform original cached labels required")
        dx=float((x[-1]-x[0])/128); h=dx/8
        data=current["profile"]
        profile=CharacteristicProfile(*[_array(data[k],(65,),"profile."+k,budget)
            for k in ("labels","values","first","second")])
        if not profile.quiet_extension: raise ValueError("original profile quiet extension required")
        y=np.linspace(0,1,65)
        charges["input_and_profile"]=budget.evaluations-initial
        cases=[]
        for case,cap_case in zip(current["cases"],old["cases"]):
            label,alpha=case["amplitude"],_scalar(case["alpha"],"alpha")
            if (alpha<0 or (label=="zero" and alpha!=0) or (label=="finite" and alpha<=0)
                    or cap_case["alpha"]!=alpha): raise ValueError("amplitude declarations differ")
            if len(case["snapshots"])!=129: raise ValueError("129 saved cuts required")
            sample=case["snapshots"][128]
            if sample["x"]!=x[-1]: raise ValueError("terminal cached time differs")
            physical=_numbers(sample["physical"]["wall"],("r","T","L","s","c"),"wall")
            for k in ("r","T","L","s"):
                if physical[k]!=cache["gauge_rows"][k][-1]+sample["state"]["wall"][k]:
                    raise ValueError("physical wall is not raw plus reference")
            before=budget.evaluations
            state=_cap_state(sample["state"]["cap"],budget)
            if len(state["Ip"])!=64: raise ValueError("64 cap cells required")
            charges["probe_analysis"]+=budget.evaluations-before
            orders={}
            for order in ORDERS:
                phase=dict(name="base metric",amplitude=label,order=order,cut=128,offset=0.)
                saved,target=_saved_terminal(cap_case,order)
                if saved["x"]!=sample["x"]: raise ValueError("saved cap time differs")
                before=budget.evaluations; counts["metric_calls"]+=1
                metric=background_cap_metric(y,physical["r"],state,T=physical["T"],profile=profile,
                    alpha=alpha,quadrature_order=order,budget=budget)
                charges["metrics"]+=budget.evaluations-before
                before=budget.evaluations
                base=_metric_row(metric,profile,state,budget)
                for k,v in base["wall"].items():
                    if v!=saved["wall"][k]: raise ValueError("base cap differs from saved "+k)
                if base["wall"]["c"]!=physical["c"]: raise ValueError("base scalar primitive differs")
                if not np.array_equal(metric["n_faces"],saved["n_faces"]): raise ValueError("base mass faces differ")
                if not np.array_equal(metric["delta_faces"],saved["delta_faces"]): raise ValueError("base lapse faces differ")
                direction=stored_direction(target["hybrid"],budget)
                charges["probe_analysis"]+=budget.evaluations-before
                phase=dict(name="base cap RHS",amplitude=label,order=order,cut=128,offset=0.)
                before=budget.evaluations; counts["rhs_calls"]+=1
                rhs=background_cap_rhs(metric,T_x=direction["T_x"],radius_x=direction["R_x"],
                    incoming_w_plus=direction["w_plus"],budget=budget)
                charges["rhs"]+=budget.evaluations-before
                before=budget.evaluations
                rates={"R":direction["R_x"],"T":direction["T_x"],**_cap_state(rhs["rates"],budget)}
                probes=[dict(offset=0.,offset_in_h=0.,**base)]
                charges["probe_analysis"]+=budget.evaluations-before
                for factor in OFFSETS[1:]:
                    phase=dict(name="backward metric probe",amplitude=label,order=order,cut=128,offset=factor*h)
                    before=budget.evaluations
                    probe=linear_probe(physical["r"],physical["T"],state,rates,factor*h,budget)
                    charges["probe_analysis"]+=budget.evaluations-before
                    before=budget.evaluations; counts["metric_calls"]+=1
                    pm=background_cap_metric(y,probe["R"],probe["state"],T=probe["T"],profile=profile,
                        alpha=alpha,quadrature_order=order,budget=budget)
                    charges["metrics"]+=budget.evaluations-before
                    before=budget.evaluations
                    row=_metric_row(pm,profile,probe["state"],budget)
                    row["limiter_changes_from_base"]=limiter_changes(base["limiter"],row["limiter"],budget)
                    probes.append(dict(offset=factor*h,offset_in_h=factor,**row))
                    charges["probe_analysis"]+=budget.evaluations-before
                phase=dict(name="directional balance",amplitude=label,order=order,cut=128)
                before=budget.evaluations
                wall=base["wall"]; R,A=wall["R"],wall["A"]
                pdiff=float(_potential(np.array([wall["f"]]),alpha)[1][0])
                # Preserve the original residual helper's operation order:
                # these are an exact saved-data compatibility check.
                Pi,Phi=wall["Pi"],wall["Phi"]
                rho=R*R*(pdiff/2+A*(Pi*Pi+Phi*Phi)/4)
                J=R*R*A*A*math.exp(wall["delta"])*wall["Pi"]*wall["Phi"]/2
                wm=wall["Pi"]-wall["Phi"]; wp=direction["w_plus"]
                Ji=R*R*A*A*math.exp(wall["delta"])*(wp*wp-wm*wm)/8
                derivatives=backward_derivatives([v["wall"]["n"] for v in probes],h,budget)
                widths=[]
                for d in derivatives:
                    radial=[]
                    for variant in target["variants"]:
                        res,curve=variant["residual"],variant["curve"]
                        if res["normalized_mass_flux"]!=J or res["normalized_radial_density"]!=rho:
                            raise ValueError("saved primary flux/density differs")
                        common=dict(R_x=direction["R_x"],T_x=direction["T_x"],n_r=res["n_r"],rho_n=rho,
                            history_n_T=res["n_T"],history_R_x=curve["R_x"],history_T_x=curve["T_x"],budget=budget)
                        left=balance_quantities(d["D_n"],J=J,**common)
                        imposed=balance_quantities(d["D_n"],J=Ji,**common)
                        if left["e_hist"]!=res["e_m"]: raise ValueError("saved history residual differs")
                        radial.append(dict(time_stride=variant["time_stride"],radial_stride=variant["radial_stride"],
                            left=left,imposed=imposed,E_m_dir=alpha*left["e_m_dir"]))
                    Bleft=d["D_n"]-direction["R_x"]*rho-direction["T_x"]*J
                    Bimposed=d["D_n"]-direction["R_x"]*rho-direction["T_x"]*Ji
                    widths.append(dict(**d,B_mass_left=Bleft,B_mass_imposed=Bimposed,
                        actual_B_mass_left=alpha*Bleft,actual_B_mass_imposed=alpha*Bimposed,radial=radial))
                orders[str(order)]=dict(base_state=state,direction=direction,cap_rates=rhs["rates"],
                    rhs_wall={k:float(rhs[k]) for k in ("wall_chi_x","wall_f_x","wall_c_x","baseline_wall_chi_x","constructed_wall_flux_defect")},
                    normalized_radial_density=rho,potential_difference_over_alpha=pdiff,
                    normalized_flux_left=J,normalized_flux_imposed=Ji,normalized_flux_imposed_minus_left=Ji-J,
                    actual_flux_imposed_minus_left=alpha*(Ji-J),probes=probes,widths=widths,
                    width_differences=[_width_difference(widths[i+1],widths[i],budget) for i in (0,1)])
                budget.check(20)
                charges["probe_analysis"]+=budget.evaluations-before
            before=budget.evaluations
            quadrature=[_width_difference(a,b,budget) for a,b in zip(orders["8"]["widths"],orders["4"]["widths"])]
            charges["probe_analysis"]+=budget.evaluations-before
            cases.append(dict(amplitude=label,alpha=alpha,orders=orders,quadrature_differences=quadrature))
        if counts!=dict(metric_calls=20,rhs_calls=4): raise ValueError("unexpected functional probe inventory")
        total=budget.evaluations-initial
        if total!=sum(charges.values()): raise ValueError("functional accounting mismatch")
        budget.check()
        return _plain(dict(schema=SCHEMA,scope=SCOPE,cut=128,saved_dx=dx,h=h,offsets_in_h=OFFSETS,
            quadrature_orders=ORDERS,cases=cases,evaluated_points=total,
            execution=dict(**counts,charges=charges,profile_constructions=1,profile_fits=0,stage_calls=0,
                prefix_calls=0,hybrid_calls=0,source_solves=0,gauge_quadratures=0,evolution_calls=0),
            comparison_signs=dict(width="finer width minus preceding width",quadrature="order8 minus order4"),
            limitations="backward functional probes; exact limiter ties and changed branches may prevent a classical tangent; no inferred convergence order"))
    except Exception as error:
        error.cap_mass_direction_failure=_plain(dict(phase=phase,counts=counts,
            initial_evaluations=initial,evaluations=budget.evaluations,message=str(error)))
        raise
