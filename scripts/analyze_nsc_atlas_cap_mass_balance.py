"""Four frozen terminal cap ledgers; no probes, time steps or profile fits."""
import math
import time

import numpy as np

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_background_cap import (
    background_cap_metric, background_cap_profile, background_cap_rhs,
)
from recursive_horizons.nsc_atlas_cap_mass_tangent import mass_tangent, mass_balance_ledger
from recursive_horizons.nsc_atlas_polar_cap import _potential
from analyze_nsc_atlas_cap_consistency import _plain, _scalar, _array, _numbers
from analyze_nsc_atlas_cap_mass_direction import _saved_terminal

SCHEMA = "nsc-atlas-cap-mass-balance-analysis-v1"
SCOPE = "terminal one-sided discrete mass derivative and departure-flux ledger; fixed state, no evolution"
EXPECTED_CALLS = dict(metric_calls=4, rhs_calls=4, tangent_calls=4, ledger_calls=4,
    profile_constructions=1, profile_fits=0, stage_calls=0, prefix_calls=0,
    hybrid_calls=0, source_solves=0, gauge_quadratures=0, evolution_calls=0)
CHARGES = ("input_and_profile", "metrics", "rhs", "tangents", "ledger_analysis")


def bounded_budget(budget):
    if (not isinstance(budget, Budget) or type(budget.maximum_evaluations) is not int
            or not 0 < budget.maximum_evaluations <= 12_000_000
            or type(budget.evaluations) is not int or budget.evaluations < 0
            or isinstance(budget.cpu_seconds, (bool, np.bool_))
            or not math.isfinite(budget.cpu_seconds) or not 0 < budget.cpu_seconds <= 90
            or isinstance(budget.started, (bool, np.bool_)) or not math.isfinite(budget.started)
            or not 0 <= budget.started <= time.process_time()):
        raise ValueError("one live cap mass-balance allowance is required")
    budget.check()
    return budget


def same_values(actual, saved, name):
    """Platform arithmetic tolerance; scientific errors retain their own scale."""
    a, b = np.asarray(actual, dtype=float), np.asarray(saved, dtype=float)
    if (a.shape != b.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b))
            or not np.allclose(a, b, rtol=2e-12, atol=2e-13)):
        raise ValueError("saved terminal reconstruction differs: "+name)


def balance_result(metric, rhs, tangent, target, older, direction, *, budget):
    Rx, Tx = direction["R_x"], direction["T_x"]
    ledger = mass_balance_ledger(metric, rhs, tangent, R_x=Rx, T_x=Tx, budget=budget)
    scale = max(1., *(abs(ledger[k]) for k in
        ("radius", "clock", "centre", "wall", "interior", "reaction", "interior_central", "interior_upwind")))
    rounding = 256*np.finfo(float).eps*scale
    defects = dict(forward_reverse=tangent["forward"]-tangent["reverse"],
        derivative_ledger=tangent["reverse"]-ledger["total"],
        telescoping=ledger["telescoping_defect"], interior_split=ledger["interior_split_defect"])
    # Face splitting is checked before contraction; scale it by its own fluxes.
    flux_scale = max(1., float(np.max(abs(rhs["p_flux"]))), float(np.max(abs(rhs["q_flux"]))))
    face_rounding = 256*np.finfo(float).eps*flux_scale
    contraction_rounding = tangent["contraction_roundoff_allowance"]
    ratio = abs(defects["forward_reverse"])/contraction_rounding if contraction_rounding else (
        0. if defects["forward_reverse"]==0 else None)
    if (abs(defects["forward_reverse"]) > contraction_rounding
            or any(abs(v)>rounding for k,v in defects.items() if k!="forward_reverse")
            or ledger["maximum_face_split_defect"] > face_rounding):
        error = ValueError("discrete mass ledger identity failed")
        error.mass_identity_failure = _plain(dict(defects=defects, rounding_allowance=rounding,
            face_rounding_allowance=face_rounding, contraction_roundoff_allowance=contraction_rounding,
            contraction_roundoff_ratio=ratio,
            ledger=ledger, tangent=tangent))
        raise error
    R, alpha, wall = metric["radius"], metric["alpha"], metric["wall"]
    fields = background_cap_profile(metric, np.array([R]), budget=budget)
    Pi, Phi = (float(fields[k][0]) for k in ("Pi", "Phi"))
    A, delta = wall["A"], wall["delta"]
    pdiff = float(_potential(np.array([wall["f"]]), alpha)[1][0])
    rho = R*R*(pdiff/2+A*(Pi*Pi+Phi*Phi)/4)
    J = R*R*A*A*math.exp(delta)*Pi*Phi/2
    Ji = R*R*A*A*math.exp(delta)*(direction["w_plus"]**2-(Pi-Phi)**2)/8
    D = tangent["reverse"]
    balance = D-Rx*rho-Tx*J
    radial = []
    for variant in target["variants"]:
        res, curve = variant["residual"], variant["curve"]
        same_values([rho, J], [res["normalized_radial_density"], res["normalized_mass_flux"]], "density/flux")
        er = res["n_r"]-rho
        em = (D-Rx*res["n_r"])/Tx-J
        history_D = curve["T_x"]*res["n_T"]+curve["R_x"]*res["n_r"]
        history_em = (history_D-Rx*res["n_r"])/Tx-J
        radial.append(dict(time_stride=variant["time_stride"], radial_stride=variant["radial_stride"],
            e_r=er, e_m=em, actual_e_m=alpha*em, history_e_m_same_rates=history_em,
            tangent_minus_history=em-history_em, identity_defect=balance-Tx*em-Rx*er))
    before = older["widths"][-1]
    budget.check(80)
    return dict(direction=direction, tangent=tangent, ledger=ledger,
        checks=dict(defects=defects, rounding_allowance=rounding, face_rounding_allowance=face_rounding,
            contraction_roundoff_allowance=contraction_rounding, contraction_roundoff_ratio=ratio),
        radial_density=rho, flux_left=J, flux_imposed=Ji, flux_change=Ji-J,
        B_mass=balance, actual_B_mass=alpha*balance, B_mass_imposed=D-Rx*rho-Tx*Ji,
        B_mass_without_upwind=balance-ledger["interior_upwind"], radial=radial,
        previous_finest_width=before["width"], previous_D_n=before["D_n"],
        tangent_minus_previous=D-before["D_n"],
        incoming_mass_rate_derivative=ledger["incoming_mass_rate_derivative"],
        interpretation="removing one ledger term is an algebraic diagnostic, not an evolved alternative")


def analyze(current, cap_audit, cache, prior_summary, budget):
    bounded_budget(budget)
    start = budget.evaluations
    charges = dict.fromkeys(CHARGES, 0)
    counts = dict.fromkeys(("metric_calls", "rhs_calls", "tangent_calls", "ledger_calls"), 0)
    phase = dict(name="fixed input declarations")
    try:
        if [(v["amplitude"], v["cells"], v["steps"]) for v in current["cases"]] != [("zero",64,512),("finite",64,512)]:
            raise ValueError("fixed zero/finite 64-cell terminal cases required")
        old = cap_audit["analysis"]
        if old["schema"] != "nsc-atlas-cap-consistency-analysis-v1":
            raise ValueError("cap audit schema differs")
        prior = prior_summary["data"]
        if prior["schema"] != "nsc-atlas-cap-mass-direction-analysis-v1":
            raise ValueError("direction summary schema differs")
        if any([r["amplitude"] for r in item["cases"]] != ["zero", "finite"] for item in (old, prior)):
            raise ValueError("comparison amplitude inventory differs")
        x = _array(cache["x_nodes"], (1025,), "cache labels", budget)
        if np.any(np.diff(x)<=0) or not np.allclose(np.diff(x), (x[-1]-x[0])/1024, rtol=0, atol=8*np.finfo(float).eps):
            raise ValueError("uniform original cached labels required")
        profile = CharacteristicProfile(*[_array(current["profile"][k], (65,), "profile."+k, budget)
            for k in ("labels", "values", "first", "second")])
        if not profile.quiet_extension:
            raise ValueError("original quiet extension required")
        charges["input_and_profile"] = budget.evaluations-start
        cases = []
        for case, cap_case, prior_case in zip(current["cases"], old["cases"], prior["cases"]):
            before = budget.evaluations
            label, alpha = case["amplitude"], _scalar(case["alpha"], "alpha")
            if (alpha < 0 or (label=="zero" and alpha!=0) or (label=="finite" and alpha<=0)
                    or cap_case["alpha"]!=alpha or prior_case["alpha"]!=alpha):
                raise ValueError("amplitude declarations differ")
            if len(case["snapshots"]) != 129:
                raise ValueError("129 saved cuts required")
            sample = case["snapshots"][128]
            if sample["x"] != x[-1]:
                raise ValueError("terminal label differs")
            physical = _numbers(sample["physical"]["wall"], ("r","T","L","s","c"), "wall")
            for k in ("r","T","L","s"):
                same_values(physical[k], cache["gauge_rows"][k][-1]+sample["state"]["wall"][k], "physical."+k)
            state = {k:_array(sample["state"]["cap"][k], (64,), "cap."+k, budget).copy() for k in ("Ip","Iq")}
            charges["input_and_profile"] += budget.evaluations-before
            orders = {}
            for order in (4,8):
                phase = dict(name="base metric", amplitude=label, order=order, cut=128)
                saved, target = _saved_terminal(cap_case, order)
                if saved["x"] != sample["x"]:
                    raise ValueError("saved cap label differs")
                before = budget.evaluations; counts["metric_calls"] += 1
                metric = background_cap_metric(np.linspace(0,1,65), physical["r"], state,
                    T=physical["T"], profile=profile, alpha=alpha, quadrature_order=order, budget=budget)
                charges["metrics"] += budget.evaluations-before
                same_values(metric["n_faces"], saved["n_faces"], "mass faces")
                same_values(metric["delta_faces"], saved["delta_faces"], "lapse faces")
                for k in ("A","delta","f","chi","c"):
                    same_values(metric["wall"][k], saved["wall"][k], "wall."+k)
                same_values(metric["wall"]["c"], physical["c"], "physical scalar")
                hybrid = _numbers(target["hybrid"], ("r_s","T_s","Vdot_plus","w_plus"), "hybrid")
                if hybrid["T_s"]<=0 or hybrid["Vdot_plus"]<=0:
                    raise ValueError("positive stored clock orientation required")
                direction = dict(R_x=hybrid["r_s"]/hybrid["Vdot_plus"], T_x=hybrid["T_s"]/hybrid["Vdot_plus"], w_plus=hybrid["w_plus"])
                phase["name"] = "base departure RHS"
                before = budget.evaluations; counts["rhs_calls"] += 1
                rhs = background_cap_rhs(metric, T_x=direction["T_x"], radius_x=direction["R_x"], incoming_w_plus=direction["w_plus"], budget=budget)
                charges["rhs"] += budget.evaluations-before
                phase["name"] = "mass edge tangent"
                before = budget.evaluations; counts["tangent_calls"] += 1
                tangent = mass_tangent(metric, dict(R=direction["R_x"], T=direction["T_x"], **rhs["rates"]), budget=budget)
                charges["tangents"] += budget.evaluations-before
                phase["name"] = "departure flux ledger"
                before = budget.evaluations; counts["ledger_calls"] += 1
                orders[str(order)] = balance_result(metric, rhs, tangent, target, prior_case["orders"][str(order)], direction, budget=budget)
                charges["ledger_analysis"] += budget.evaluations-before
            differences = {k:orders["8"][k]-orders["4"][k] for k in
                ("B_mass", "actual_B_mass", "B_mass_without_upwind", "tangent_minus_previous")}
            cases.append(dict(amplitude=label, alpha=alpha, orders=orders, order8_minus_order4=differences))
        if counts != {k:EXPECTED_CALLS[k] for k in counts} or sum(charges.values()) != budget.evaluations-start:
            raise ValueError("mass balance execution accounting mismatch")
        budget.check()
        return _plain(dict(schema=SCHEMA, scope=SCOPE, cut=128, quadrature_orders=[4,8], cases=cases,
            execution=dict(EXPECTED_CALLS, charges=charges), evaluated_points=budget.evaluations-start,
            limitations="one common backward branch; fixed profile/state; order difference is not an error bound"))
    except Exception as error:
        error.cap_mass_balance_failure = _plain(dict(phase=phase, counts=counts, evaluated_points=budget.evaluations-start,
            identity_failure=getattr(error,"mass_identity_failure",None)))
        raise
