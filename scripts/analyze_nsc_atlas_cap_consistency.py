"""Selected saved-cap constraints and independently normalized wall curves.

This analyzer owns no source, time evolution, or full coupled-stage call.
The caller authenticates its input and owns the single cumulative budget.
Finite differences and local quadrature sensitivities are not error bounds.
"""
from __future__ import annotations

from collections.abc import Mapping
from numbers import Real
import math
import time

import numpy as np

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_background_cap import background_cap_metric, background_cap_profile
from recursive_horizons.nsc_atlas_departure_gradient import quadratic_departure_gradient
from recursive_horizons.nsc_atlas_nonlinear_boundary import hybrid_wall_boundary
from recursive_horizons.nsc_atlas_polar_cap import _potential
from recursive_horizons.nsc_atlas_vacuum_departure import vacuum_slice, departure_slice

SCHEMA = "nsc-atlas-cap-consistency-analysis-v1"
TARGETS = (63, 68, 128)
NEIGHBORS = tuple(range(61, 71)) + (122, 124, 125, 126, 127, 128)
ORDERS = (4, 8)
WIDTHS = (1, 2)
FIELDS = ("r", "L", "c")
SCOPE = ("selected reconstructed cap left traces, mass/lapse constraints and actual wall curve; "
         "quadrature/width sensitivity, not a continuum or junction certificate")


def _budget(budget):
    if (not isinstance(budget, Budget) or type(budget.maximum_evaluations) is not int
            or not 0 < budget.maximum_evaluations <= 8_000_000
            or type(budget.evaluations) is not int or budget.evaluations < 0
            or isinstance(budget.cpu_seconds, bool) or not math.isfinite(budget.cpu_seconds)
            or not 0 < budget.cpu_seconds <= 120
            or isinstance(budget.started, bool) or not math.isfinite(budget.started)
            or not 0 <= budget.started <= time.process_time()):
        raise ValueError("one live bounded cap-consistency budget is required")
    budget.check()
    return budget


def _scalar(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(name+" must be a finite real scalar")
    return float(value)


def _array(value, shape, name, budget):
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf" or raw.shape != shape or not np.all(np.isfinite(raw)):
        raise ValueError(name+" has invalid shape or nonfinite values")
    budget.check(int(raw.size))
    return np.asarray(raw, dtype=float)


def _numbers(row, keys, name):
    if not isinstance(row, Mapping) or not set(keys) <= set(row):
        raise ValueError(name+" is missing fields")
    return {k: _scalar(row[k], name+"."+k) for k in keys}


def _plain(value):
    if isinstance(value, dict): return {k: _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [_plain(v) for v in value]
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    return value


def temporal_stencil(cut, stride, derivative):
    """Declared centred or terminal-backward weights; divide by dx**order."""
    if type(cut) is not int or cut not in TARGETS or type(stride) is not int or stride not in WIDTHS:
        raise ValueError("only admitted target cuts and widths are available")
    if type(derivative) is not int or derivative not in (1, 2):
        raise ValueError("derivative must be one or two")
    if cut == 128:
        offsets, weights = (([0, -1, -2], [3/2, -2, 1/2]) if derivative == 1
                            else ([0, -1, -2, -3], [2, -5, 4, -1]))
    else:
        offsets = [-1, 0, 1]
        weights = [-1/2, 0, 1/2] if derivative == 1 else [1, -2, 1]
    indices = cut+stride*np.asarray(offsets, dtype=int)
    return indices, np.asarray(weights, dtype=float)/stride**derivative


def differentiate_samples(values, indices, cut, stride, derivative, dx, budget):
    _budget(budget)
    ids = np.asarray(indices)
    if (ids.ndim != 1 or not 3 <= len(ids) <= 129 or ids.dtype.kind not in "iu"
            or len(np.unique(ids)) != len(ids) or np.any(ids < 0) or np.any(ids > 128)):
        raise ValueError("unique saved-cut indices required")
    data = _array(values, ids.shape, "saved values", budget)
    dx = _scalar(dx, "dx")
    if dx <= 0: raise ValueError("positive saved time spacing required")
    selected, weights = temporal_stencil(cut, stride, derivative)
    locations = {int(k): i for i, k in enumerate(ids)}
    if any(int(k) not in locations for k in selected):
        raise ValueError("saved temporal stencil is incomplete")
    budget.check(len(weights))
    return float(np.dot(weights, [data[locations[int(k)]] for k in selected])/dx**derivative)


def radial_first(values, spacing, stride, budget):
    """Backward-three-point derivative at the actual last radial face."""
    _budget(budget)
    raw = np.asarray(values)
    if raw.ndim != 1 or not 5 <= len(raw) <= 129:
        raise ValueError("bounded radial face vector required")
    data = _array(values, raw.shape, "radial values", budget)
    spacing = _scalar(spacing, "radial spacing")
    if spacing <= 0 or type(stride) is not int or stride not in WIDTHS:
        raise ValueError("positive radial spacing and admitted stride required")
    budget.check(3)
    return float((3*data[-1]-4*data[-1-stride]+data[-1-2*stride])/(2*stride*spacing))


def residual_quantities(wall, *, n_T, n_r, delta_r, alpha, budget):
    """Mass/lapse defects using supplied measured jets, never desired jets."""
    _budget(budget)
    w = _numbers(wall, ("R", "A", "delta", "n", "f", "Pi", "Phi"), "wall")
    alpha = _scalar(alpha, "alpha")
    nt, nr, dr = (_scalar(v, k) for v, k in ((n_T, "n_T"), (n_r, "n_r"), (delta_r, "delta_r")))
    if alpha < 0 or w["R"] <= 0 or w["A"] <= 0: raise ValueError("positive cap geometry required")
    R, A, delta, Pi, Phi = (w[k] for k in ("R", "A", "delta", "Pi", "Phi"))
    pdiff = float(_potential(np.array([w["f"]]), alpha)[1][0])
    mass_flux = R*R*A*A*math.exp(delta)*Pi*Phi/2
    radial_density = R*R*(pdiff/2+A*(Pi*Pi+Phi*Phi)/4)
    lapse_density = alpha*R*(Pi*Pi+Phi*Phi)/2
    em, er = nt-mass_flux, nr-radial_density
    result = dict(e_m=em, e_r=er, E_m=alpha*em, E_r=alpha*er,
        E_delta=dr-lapse_density, n_T=nt, n_r=nr, delta_r=dr,
        normalized_mass_flux=mass_flux, normalized_radial_density=radial_density,
        lapse_density=lapse_density, potential_difference_over_alpha=pdiff,
        m_T=alpha*nt, m_r=3*R*R/2+alpha*nr)
    if not all(math.isfinite(v) for v in result.values()): raise ValueError("nonfinite cap residual")
    budget.check(20)
    return result


def curve_quantities(wall, first, second, radial, *, alpha, budget):
    """Normalize the measured curve and compute full geometric acceleration.

    first/second contain derivatives in the original x label, not helper
    predictions. n_T is fixed areal radius; no Einstein mass law is used.
    """
    _budget(budget)
    w = _numbers(wall, ("R", "A", "delta", "n"), "wall")
    d = _numbers(first, ("R", "T", "s", "n", "A", "delta", "f", "chi"), "first")
    dd = _numbers(second, ("R", "T"), "second")
    radial = _numbers(radial, ("n", "delta"), "radial")
    alpha = _scalar(alpha, "alpha")
    R, A, delta, n = (w[k] for k in ("R", "A", "delta", "n"))
    if alpha < 0 or R <= 0 or A <= 0 or d["T"] <= 0:
        raise ValueError("positive cap geometry and measured T_x are required")
    e = math.exp(delta)
    sigma2 = A*e*e*d["T"]**2-d["R"]**2/A
    if not math.isfinite(sigma2) or sigma2 <= 0:
        raise ValueError("measured wall curve is not timelike")
    sigma = math.sqrt(sigma2)
    Ax = d["A"]
    Ax_chain = -2*R*d["R"]-2*alpha*(d["n"]/R-n*d["R"]/R**2)
    sigma_x = (e*e*(Ax+2*A*d["delta"])*d["T"]**2
        +2*A*e*e*d["T"]*dd["T"]-2*d["R"]*dd["R"]/A+d["R"]**2*Ax/A**2)/(2*sigma)
    v, beta = d["R"]/sigma, A*e*d["T"]/sigma
    vs = (dd["R"]*sigma-d["R"]*sigma_x)/sigma**3
    nt = (d["n"]-d["R"]*radial["n"])/d["T"]
    Ar = -2*R-2*alpha*(radial["n"]/R-n/R**2)
    AT = -2*alpha*nt/R
    K = (vs+Ar/2)/beta+beta*radial["delta"]-v*AT/(A*A*e)
    result = dict(sigma=sigma, sigma_squared=sigma2, sigma_x=sigma_x, v=v, beta=beta,
        v_s=vs, n_T=nt, A_r=Ar, A_T=AT, A_W_x=Ax, A_W_x_chain=Ax_chain,
        A_W_x_chain_defect=Ax-Ax_chain, K_geom=K,
        proper_clock_derivative_defect=sigma-d["s"],
        beta_identity_defect=beta*beta-v*v-A, T_s=d["T"]/sigma,
        mass_motion_subtraction=d["R"]*radial["n"],
        R_x=d["R"], T_x=d["T"], f_x=d["f"], chi_x=d["chi"])
    if not all(math.isfinite(v) for v in result.values()): raise ValueError("nonfinite measured curve")
    budget.check(35)
    return result


def trace_quantities(wall, curve, hybrid, budget):
    """Keep the reconstructed, kinematic and imposed incoming traces distinct."""
    _budget(budget)
    w = _numbers(wall, ("R", "A", "delta", "f", "Pi", "Phi", "alpha"), "wall")
    c = _numbers(curve, ("R_x", "T_x", "f_x", "chi_x", "v", "T_s", "K_geom"), "curve")
    h = _numbers(hybrid, ("w_plus", "w_minus", "k_minus", "r_s", "T_s", "Vdot_plus"), "hybrid")
    R, A, e = w["R"], w["A"], math.exp(w["delta"])
    wm, wp = w["Pi"]-w["Phi"], w["Pi"]+w["Phi"]
    denominator = A*e*c["T_x"]+c["R_x"]
    if denominator <= 0 or h["Vdot_plus"] <= 0: raise ValueError("invalid incoming orientation")
    kin = (2*c["f_x"]-(A*e*c["T_x"]-c["R_x"])*wm)/denominator
    Pih, Phih = (h["w_plus"]+wm)/2, (h["w_plus"]-wm)/2
    Pik, Phik = (kin+wm)/2, (kin-wm)/2
    factor = R*R*A*A*e/2
    left_flux, imposed_flux, kinetic_flux = (factor*a*b for a,b in ((w["Pi"],w["Phi"]),(Pih,Phih),(Pik,Phik)))
    scalar_defect = c["f_x"]-A*e*c["T_x"]*w["Pi"]-c["R_x"]*w["Phi"]
    chi_defect = c["chi_x"]-A*e*c["T_x"]*R*w["Pi"]-c["R_x"]*(w["f"]+R*w["Phi"])
    result = dict(w_minus_left=wm, w_plus_left=wp, w_plus_kinematic=kin, w_plus_hybrid=h["w_plus"],
        kinematic_minus_left=kin-wp, hybrid_minus_left=h["w_plus"]-wp,
        hybrid_minus_kinematic=h["w_plus"]-kin, outgoing_identity_defect=h["w_minus"]-wm,
        Pi_hybrid=Pih, Phi_hybrid=Phih, Pi_kinematic=Pik, Phi_kinematic=Phik,
        normalized_flux_left=left_flux, normalized_flux_hybrid=imposed_flux,
        normalized_flux_kinematic=kinetic_flux,
        normalized_flux_hybrid_minus_left=imposed_flux-left_flux,
        normalized_flux_kinematic_minus_left=kinetic_flux-left_flux,
        actual_flux_hybrid_minus_left=w["alpha"]*(imposed_flux-left_flux),
        scalar_kinematic_defect=scalar_defect, canonical_kinematic_defect=chi_defect,
        primitive_product_derivative_defect=c["chi_x"]-c["R_x"]*w["f"]-R*c["f_x"],
        K_geom_minus_hybrid=c["K_geom"]-h["k_minus"],
        v_curve_minus_hybrid=c["v"]-h["r_s"], T_s_curve_minus_hybrid=c["T_s"]-h["T_s"],
        R_x_curve_minus_hybrid=c["R_x"]-h["r_s"]/h["Vdot_plus"],
        T_x_curve_minus_hybrid=c["T_x"]-h["T_s"]/h["Vdot_plus"])
    if not all(math.isfinite(v) for v in result.values()): raise ValueError("nonfinite trace comparison")
    budget.check(40)
    return result


def selected_hybrid(snapshot, metric, source, gauge, alpha, tube, budget, *, calls=None):
    """One saved-data prefix and algebraic map; no stage or cap RHS."""
    _budget(budget)
    x = _scalar(snapshot["x"], "snapshot x")
    g = _numbers(gauge, ("x", "W", "W_x", "eta_plus", "r", "L", "T", "s"), "gauge")
    src = _numbers(source, ("x", "r", "log_C", "c", "r_x", "log_C_x", "c_x"), "source")
    if g["x"] != x or src["x"] != x or g["W"] <= 0 or g["W_x"] <= 0:
        raise ValueError("selected cache coordinates/width differ")
    if not math.isclose(g["eta_plus"], -.5*math.log(g["W_x"]), rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("selected boost differs from fixed gauge")
    raw = snapshot["state"]; eta = np.linspace(0, 1, 65)
    wall = {k: g[k]+_scalar(raw["wall"][k], "wall departure") for k in ("r", "L", "T", "s")}
    wall["c"] = float(metric["wall"]["c"])
    if wall["r"] != metric["radius"] or wall["T"] != metric["T"]:
        raise ValueError("cap metric and selected wall differ")
    vac = vacuum_slice(eta*g["W"], x, budget=budget)
    departures, X = {}, {}
    for k in FIELDS:
        z = _array(raw["departure"][k], (63,), "saved departure."+k, budget)
        baseline = vac[k] if k != "c" else np.zeros(65)
        start = src["log_C" if k == "L" else k]
        offset = start-baseline[0]
        departures[k] = np.r_[0., z, wall[k]-baseline[-1]-offset]
        X[k] = baseline+offset+departures[k]
        X[k][0], X[k][-1] = start, wall[k]
    G = quadratic_departure_gradient(eta, np.stack([departures[k] for k in FIELDS]), budget=budget)
    jets = dict(r_Ubar=G[0]/g["W"], log_C_Ubar=G[1]/g["W"], c_Ubar=G[2]/g["W"])
    if calls is not None: calls["prefix_calls"] += 1
    ext = departure_slice(eta*g["W"], X, jets, src, x=x, alpha=alpha,
                          tube=tube, quadrature_order=2, budget=budget)
    outgoing = dict(r_V=float(ext["r_x"][-1]), c_V=float(ext["c_x"][-1]),
                    log_C_V=float(ext["log_C_x"][-1]))
    if calls is not None: calls["hybrid_calls"] += 1
    mapped = hybrid_wall_boundary(dict(r=wall["r"], c=wall["c"], eta_plus=g["eta_plus"]),
        outgoing, {k: metric["wall"][k] for k in ("A", "delta", "w_minus")},
        dict(C_wall=math.exp(wall["L"]), log_C_s=0.), alpha=alpha)
    names = ("w_plus", "Pi", "Phi", "r_s", "T_s", "f_s", "c_s", "k_minus", "k_plus",
             "beta_minus", "beta_plus", "H", "Vdot_plus", "Udot_plus")
    result = {k: float(mapped[k]) for k in names}
    result.update(w_minus=float(metric["wall"]["w_minus"]), **outgoing)
    budget.check(12)
    return result


def _differences(fine, coarse, budget):
    result = {}
    for group in ("residual", "curve", "traces"):
        result[group] = {k: fine[group][k]-coarse[group][k] for k in fine[group]}
        budget.check(len(result[group]))
    return result


def analyze(current, cache, budget):
    """Deterministic64-metric/12-map selected analysis of authenticated data.

    The caller preserves prior budget usage. Returned counts are LOCAL, so
    validation on that same used budget produces the same result.
    """
    _budget(budget); initial = budget.evaluations
    phase = dict(name="input declarations")
    counts = dict(metric_calls=0, prefix_calls=0, hybrid_calls=0)
    charges = dict(input_and_profile=0, metrics=0, selected_prefix_maps=0, derivative_analysis=0)
    try:
        if [(v["amplitude"], v["cells"], v["steps"]) for v in current["cases"]] != [("zero",64,512),("finite",64,512)]:
            raise ValueError("require the two fixed quadratic64/512 cases")
        x = _array(cache["x_nodes"], (1025,), "cache x", budget)
        if np.any(np.diff(x) <= 0) or not np.allclose(np.diff(x), (x[-1]-x[0])/1024, rtol=0, atol=8*np.finfo(float).eps):
            raise ValueError("uniform increasing cached labels required")
        saved_x = x[::8]; dx = float((saved_x[-1]-saved_x[0])/128)
        data = current["profile"]
        profile_arrays = [_array(data[k], (65,), "profile."+k, budget) for k in ("labels","values","first","second")]
        profile = CharacteristicProfile(*profile_arrays)
        if not profile.quiet_extension: raise ValueError("retained profile must have its quiet lower extension")
        y = np.linspace(0, 1, 65)
        tube = current["settings"]["tube"]
        charges["input_and_profile"] = budget.evaluations-initial
        output = []
        for case in current["cases"]:
            label, alpha = case["amplitude"], _scalar(case["alpha"], "alpha")
            if alpha < 0 or (label == "zero" and alpha != 0) or (label == "finite" and alpha <= 0):
                raise ValueError("amplitude labels differ")
            if len(case["snapshots"]) != 129: raise ValueError("129 snapshots required")
            orders = {}
            for order in ORDERS:
                slices, metrics, maps = {}, {}, {}
                for cut in NEIGHBORS:
                    phase = dict(name="cap metric", amplitude=label, order=order, cut=cut)
                    sample = case["snapshots"][cut]
                    if sample["x"] != saved_x[cut]: raise ValueError("snapshot time differs from cached key")
                    physical = _numbers(sample["physical"]["wall"], ("r","T","s","L","c"), "physical wall")
                    gauge = {k: v[8*cut] for k,v in cache["gauge_rows"].items()}
                    if gauge["x"] != sample["x"]: raise ValueError("neighbor gauge time differs from snapshot")
                    for k in ("r","T","s","L"):
                        if physical[k] != gauge[k]+sample["state"]["wall"][k]:
                            raise ValueError("saved wall is not raw plus reference")
                    before = budget.evaluations
                    counts["metric_calls"] += 1
                    metric = background_cap_metric(y, physical["r"], sample["state"]["cap"],
                        T=physical["T"], profile=profile, alpha=alpha, quadrature_order=order, budget=budget)
                    left = background_cap_profile(metric, np.array([physical["r"]]), budget=budget)
                    if metric["wall"]["c"] != physical["c"]:
                        raise ValueError("saved common scalar differs from reconstructed primitive")
                    wall = dict(R=physical["r"], T=physical["T"], s=physical["s"],
                        A=float(metric["wall"]["A"]), delta=float(metric["wall"]["delta"]),
                        n=float(metric["n_faces"][-1]), f=float(metric["wall"]["f"]),
                        Pi=float(left["Pi"][0]), Phi=float(left["Phi"][0]),
                        chi=float(metric["wall"]["chi"]), c=float(metric["wall"]["c"]), alpha=alpha)
                    slices[cut] = dict(cut=cut, x=float(sample["x"]), wall=wall,
                        n_faces=np.asarray(metric["n_faces"]), delta_faces=np.asarray(metric["delta_faces"]),
                        panels=len(metric["panel_cuts"])-1,
                        upper_profile_margin=profile.domain["last_label"]-(physical["T"]+np.arctanh(physical["r"])))
                    charges["metrics"] += budget.evaluations-before
                    if cut in TARGETS: metrics[cut] = metric
                for cut in TARGETS:
                    phase = dict(name="selected prefix and hybrid", amplitude=label, order=order, cut=cut)
                    source = {k:v[8*cut] for k,v in cache["source_rows"][label].items()}
                    gauge = {k:v[8*cut] for k,v in cache["gauge_rows"].items()}
                    before = budget.evaluations
                    maps[cut] = selected_hybrid(case["snapshots"][cut], metrics[cut], source, gauge, alpha, tube, budget, calls=counts)
                    charges["selected_prefix_maps"] += budget.evaluations-before
                before = budget.evaluations
                targets = []
                for cut in TARGETS:
                    variants = []
                    for ts in WIDTHS:
                        phase = dict(name="wall curve derivatives", amplitude=label, order=order, cut=cut, time_stride=ts)
                        d = {k:differentiate_samples([slices[i]["wall"][k] for i in NEIGHBORS], NEIGHBORS,
                            cut, ts, 1, dx, budget) for k in ("R","T","s","n","A","delta","f","chi")}
                        dd = {k:differentiate_samples([slices[i]["wall"][k] for i in NEIGHBORS], NEIGHBORS,
                            cut, ts, 2, dx, budget) for k in ("R","T")}
                        for rs in WIDTHS:
                            phase["radial_stride"] = rs
                            local = slices[cut]; wall = local["wall"]
                            dr = {"n":radial_first(local["n_faces"],wall["R"]/64,rs,budget),
                                  "delta":radial_first(local["delta_faces"],wall["R"]/64,rs,budget)}
                            curve = curve_quantities(wall,d,dd,dr,alpha=alpha,budget=budget)
                            residual = residual_quantities(wall,n_T=curve["n_T"],n_r=dr["n"],delta_r=dr["delta"],alpha=alpha,budget=budget)
                            traces = trace_quantities(wall,curve,maps[cut],budget)
                            variants.append(dict(time_stride=ts,radial_stride=rs,residual=residual,curve=curve,traces=traces))
                    by = {(v["time_stride"],v["radial_stride"]):v for v in variants}
                    targets.append(dict(cut=cut,hybrid=maps[cut],variants=variants,
                        time_width_differences={str(rs):_differences(by[2,rs],by[1,rs],budget) for rs in WIDTHS},
                        radial_width_differences={str(ts):_differences(by[ts,2],by[ts,1],budget) for ts in WIDTHS}))
                orders[str(order)] = dict(slices=[slices[i] for i in NEIGHBORS],targets=targets)
                charges["derivative_analysis"] += budget.evaluations-before
            before = budget.evaluations
            quadrature = []
            for low,high in zip(orders["4"]["targets"],orders["8"]["targets"]):
                quadrature.append(dict(cut=low["cut"],variants=[dict(time_stride=a["time_stride"],radial_stride=a["radial_stride"],
                    differences=_differences(b,a,budget)) for a,b in zip(low["variants"],high["variants"])]))
            charges["derivative_analysis"] += budget.evaluations-before
            output.append(dict(amplitude=label,alpha=alpha,orders=orders,quadrature_differences=quadrature))
        if counts != dict(metric_calls=64,prefix_calls=12,hybrid_calls=12): raise ValueError("unexpected reconstruction inventory")
        budget.check()
        total = budget.evaluations-initial
        if total != sum(charges.values()): raise ValueError("analysis accounting mismatch")
        return _plain(dict(schema=SCHEMA,scope=SCOPE,targets=TARGETS,neighbor_cuts=NEIGHBORS,
            quadrature_orders=ORDERS,time_strides=WIDTHS,radial_strides=WIDTHS,cases=output,
            evaluated_points=total,execution=dict(**counts,charges=charges,profile_constructions=1,
                profile_fits=0,stage_calls=0,rhs_calls=0,source_solves=0,gauge_quadratures=0,evolution_calls=0),
            comparison_signs=dict(quadrature="order8 minus order4",time_width="stride2 minus stride1",
                                  radial_width="stride2 minus stride1"),
            limitations="wall-only left traces; widths and quadrature are sensitivities; J/O implications remain conditional"))
    except Exception as error:
        error.cap_consistency_failure = _plain(dict(phase=phase,counts=counts,
            initial_evaluations=initial,evaluations=budget.evaluations,message=str(error)))
        raise
