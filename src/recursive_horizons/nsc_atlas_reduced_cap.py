"""Conservative canonical reduced scalar on the same regular polar cap.

The state is Ip=integral(r*Pi dr), Iq=integral((r*f)_r dr).  The common
canonical wall value chi_W=sum(Iq) is carried by the q boundary flux; there
is no second field-value ODE.  chi here is r*f, not the frame-invariant
exterior c=r*f*sinhc(sqrt(alpha/6)*f).

Away from the first cell, the two reduced characteristics p+q and p-q use
the same limited linear reconstruction.  The first cell instead has odd
p=r*(P0+P2*r**2) and even q=Q0+Q2*r**2, preserving its integrals exactly.
This is a finite-volume prototype, not a stability or continuum certificate.
"""
from collections.abc import Mapping

import numpy as np
from numpy.polynomial.legendre import leggauss

from .nsc_atlas_polar_cap import (
    MAX_CELLS, QUADRATURE_ORDERS, _context, _potential, _scalar, _vector,
)


def _inputs(y_faces, radius, state, alpha):
    y = _vector(y_faces, "y_faces")
    n = len(y)-1
    if (not 3 <= n <= MAX_CELLS or y[0] != 0 or y[-1] != 1
            or np.any(np.diff(y) <= 0)
            or not np.allclose(np.diff(y), 1/n, rtol=0, atol=8*np.finfo(float).eps)):
        raise ValueError("require 3..128 uniform cells with y endpoints 0 and 1")
    R = _scalar(radius, "radius", positive=True)
    alpha = _scalar(alpha, "alpha", nonnegative=True)
    if not isinstance(state, Mapping) or set(state) != {"Ip", "Iq"}:
        raise ValueError("state must contain exactly Ip, Iq")
    fields = {key: _vector(state[key], key).copy() for key in ("Ip", "Iq")}
    if any(len(a) != n for a in fields.values()):
        raise ValueError("each integral field must have one value per cell")
    return y.copy(), R, fields, alpha


def _limited_slopes(midpoint, values):
    secant = np.diff(values)/np.diff(midpoint)
    slope = np.empty_like(values)
    slope[0], slope[-1] = secant[0], secant[-1]
    left, right = secant[:-1], secant[1:]
    slope[1:-1] = np.where(left*right > 0,
                           np.sign(left)*np.minimum(abs(left), abs(right)), 0.)
    return slope


def _reconstruct(r, state):
    dr = np.diff(r)
    midpoint = (r[:-1]+r[1:])/2
    pbar, qbar = state["Ip"]/dr, state["Iq"]/dr
    plus, minus = pbar+qbar, pbar-qbar
    sp, sm = _limited_slopes(midpoint, plus), _limited_slopes(midpoint, minus)
    ps, qs = (sp+sm)/2, (sp-sm)/2
    pc, qc = pbar-ps*midpoint, qbar-qs*midpoint
    chi = np.r_[0., np.cumsum(state["Iq"])]
    # These polynomial moments preserve the first cell's actual dr integrals.
    q_rho_mean = (r[:-1]**2+r[:-1]*r[1:]+r[1:]**2)/3
    q2 = (qbar[1]-qbar[0])/(q_rho_mean[1]-q_rho_mean[0])
    q0 = qbar[0]-q2*q_rho_mean[0]
    pi_rho_mean = (r[:-1]**2+r[1:]**2)/2
    pi_mean = pbar/midpoint
    p2 = (pi_mean[1]-pi_mean[0])/(pi_rho_mean[1]-pi_rho_mean[0])
    p0 = pi_mean[0]-p2*pi_rho_mean[0]
    return {"r": r, "midpoint": midpoint, "plus": plus, "minus": minus,
            "plus_slope": sp, "minus_slope": sm, "pc": pc, "ps": ps,
            "qc": qc, "qs": qs, "chi": chi,
            "P0": p0, "P2": p2, "Q0": q0, "Q2": q2}


def _fields(rec, i, z):
    z = np.asarray(z, dtype=float)
    if i == 0:
        f = rec["Q0"]+rec["Q2"]*z*z/3
        Pi = rec["P0"]+rec["P2"]*z*z
        Phi = 2*rec["Q2"]*z/3
        q = rec["Q0"]+rec["Q2"]*z*z
        chi, p = z*f, z*Pi
    else:
        a = rec["r"][i]
        pc, ps, qc, qs = (rec[k][i] for k in ("pc", "ps", "qc", "qs"))
        k = rec["chi"][i]-qc*a-qs*a*a/2
        chi = rec["chi"][i]+(z-a)*(qc+qs*(z+a)/2)
        p, q = pc+ps*z, qc+qs*z
        f, Pi = chi/z, p/z
        Phi = qs/2-k/(z*z)
    return {"f": f, "Pi": Pi, "Phi": Phi, "p": p, "q": q, "chi": chi}


def _gamma(rec, i, z, alpha):
    """Exact integral of alpha*r*(Pi**2+Phi**2)/2 in one reconstruction."""
    z = np.asarray(z, dtype=float)
    if i == 0:
        p, s, q = rec["P0"], rec["P2"], rec["Q2"]
        return alpha*(p*p*z*z/4+(2*p*s+4*q*q/9)*z**4/8+s*s*z**6/12)
    a = rec["r"][i]
    pc, ps, qc, qs = (rec[k][i] for k in ("pc", "ps", "qc", "qs"))
    k = rec["chi"][i]-qc*a-qs*a*a/2
    d = z-a
    return alpha*((ps*ps+qs*qs/4)*d*(z+a)/4+ps*pc*d
                  +(pc*pc-qs*k)*np.log1p(d/a)/2
                  +k*k*d*(z+a)/(4*a*a*z*z))


def _potential_reaction_quotient(f, alpha):
    """(Ualpha(f)-2*f)/alpha with its continuous, nonzero alpha=0 limit."""
    f = np.asarray(f, dtype=float)
    z = np.sqrt(alpha/6)*f
    small = abs(z) < 1e-4
    sm1_over_z2 = np.empty_like(z)
    sm1_over_z2[small] = 1/6+z[small]**2/120+z[small]**4/5040
    np.divide(np.sinh(z)/np.where(z == 0, 1, z)-1, z*z,
              out=sm1_over_z2, where=~small)
    half_sinhc = np.ones_like(z)
    half_sinhc[small] = 1+z[small]**2/24+z[small]**4/1920
    np.divide(np.sinh(z/2), z/2, out=half_sinhc, where=~small)
    ch = np.cosh(z)
    return f**3/3*(sm1_over_z2*ch**3+half_sinhc**2*(ch*ch+ch+1)/2)


def _wall_conversion(f, alpha):
    z = np.sqrt(alpha/6)*f
    sinhc = 1+z*z/6+z**4/120 if abs(z) < 1e-4 else np.sinh(z)/z
    return float(f*sinhc)


def reduced_cap_metric(y_faces, radius, state, *, alpha, quadrature_order=4, budget=None):
    """Integrate radial constraints using the actual primitive at every node.

    m=r**3/2+alpha*n retains the vacuum exactly.  The regular centre value
    is Q0 from the even q reconstruction, not an ordinary cell average.
    Shared work is charged as N*(quadrature_order+1)**2 evaluations.
    """
    budget = _context(budget)
    y, R, state, alpha = _inputs(y_faces, radius, state, alpha)
    if (isinstance(quadrature_order, (bool, np.bool_))
            or not isinstance(quadrature_order, (int, np.integer))
            or quadrature_order not in QUADRATURE_ORDERS):
        raise ValueError("quadrature_order must be 2, 4 or 8")
    r, ncell = R*y, len(y)-1
    rec = _reconstruct(r, state)
    nodes, weights = leggauss(quadrature_order)
    rn = (r[:-1, None]+r[1:, None])/2+np.diff(r)[:, None]*nodes/2
    wn = np.diff(r)[:, None]*weights/2
    If, nf = np.zeros(ncell+1), np.zeros(ncell+1)
    In, nn = np.empty_like(rn), np.empty_like(rn)
    fields = {key: np.empty_like(rn) for key in ("f", "Pi", "Phi", "p", "q", "chi")}
    f_average = np.empty(ncell)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
            z0 = np.sqrt(alpha/6)*rec["Q0"]
            delta0 = -float(np.log1p(2*np.sinh(z0/2)**2))
            for i in range(ncell):
                a, b = r[i:i+2]
                local = _fields(rec, i, rn[i])
                for key in fields:
                    fields[key][i] = local[key]
                if i == 0:
                    f_average[i] = rec["Q0"]+rec["Q2"]*b*b/9
                else:
                    qc, qs = rec["qc"][i], rec["qs"][i]
                    k = rec["chi"][i]-qc*a-qs*a*a/2
                    f_average[i] = qc+qs*(a+b)/4+k*np.log1p((b-a)/a)/(b-a)
                previous, previous_I, response = a, If[i], nf[i]
                for j, end in enumerate(np.r_[rn[i], b]):
                    q = (previous+end)/2+(end-previous)*nodes/2
                    w = (end-previous)*weights/2
                    budget.check(len(q))
                    Iq, Ie = If[i]+_gamma(rec, i, q, alpha), If[i]+_gamma(rec, i, end, alpha)
                    data = _fields(rec, i, q)
                    Pdiff = _potential(data["f"], alpha)[1]
                    W = data["Pi"]**2+data["Phi"]**2
                    forcing = q*q*(Pdiff/2+(1-q*q)*W/4)
                    response = np.exp(previous_I-Ie)*response+np.dot(w, np.exp(Iq-Ie)*forcing)
                    if j < quadrature_order:
                        In[i, j], nn[i, j] = Ie, response
                    else:
                        If[i+1], nf[i+1] = Ie, response
                    previous, previous_I = end, Ie
            mf, mn = r**3/2+alpha*nf, rn**3/2+alpha*nn
            Af = np.ones_like(r)
            Af[1:] = 1-r[1:]**2-2*alpha*nf[1:]/r[1:]
            An = 1-rn*rn-2*alpha*nn/rn
            df, dn = delta0+If, delta0+In
            Bf, Bn = Af*np.exp(df), An*np.exp(dn)
            P, Pdiff, U = _potential(fields["f"], alpha)
            P0 = float(_potential(np.array([rec["Q0"]]), alpha)[0][0])
            chiW = float(rec["chi"][-1])
            fW = chiW/R
            amW = float(rec["minus"][-1]+rec["minus_slope"][-1]*(R-rec["midpoint"][-1]))
            tanhc = 1-z0*z0/3+2*z0**4/15 if abs(z0) < 1e-4 else np.tanh(z0)/z0
            result = {
                "y_faces": y, "radius": R, "alpha": alpha, "state": state,
                "quadrature_order": quadrature_order, "r_faces": r, "r_nodes": rn,
                "weights_nodes": wn, "m_faces": mf, "n_faces": nf,
                "A_faces": Af, "delta_faces": df, "B_faces": Bf,
                "m_nodes": mn, "n_nodes": nn, "A_nodes": An,
                "delta_nodes": dn, "B_nodes": Bn, "chi_faces": rec["chi"].copy(),
                "f_cell_average": f_average, "reconstruction": rec,
                "potential_nodes": P, "potential_difference_over_alpha_nodes": Pdiff,
                "potential_derivative_normalized_nodes": U,
                "charged_evaluations": ncell*(quadrature_order+1)**2,
                "spatial_order": 1,
                "center": {"f": float(rec["Q0"]), "Pi": float(rec["P0"]),
                           "Q": float(2*rec["Q2"]/3), "mu": P0/6+alpha*rec["P0"]**2/12,
                           "delta": delta0, "B": float(Bf[0]), "p_J": float(rec["Q0"]*tanhc)},
                "wall": {"A": float(Af[-1]), "delta": float(df[-1]), "B": float(Bf[-1]),
                         "a_minus": amW, "w_minus": (amW+fW)/R,
                         "chi": chiW, "f": fW, "c": R*_wall_conversion(fW, alpha)},
            }
            result.update({key+"_nodes": value for key, value in fields.items()})
    except FloatingPointError as error:
        raise ValueError("nonfinite reduced-cap constraint arithmetic") from error
    if (np.any(Af <= 0) or np.any(An <= 0) or np.any(Bf <= 0) or np.any(Bn <= 0)
            or any(not np.all(np.isfinite(value)) for value in fields.values())
            or any(not np.all(np.isfinite(result[k])) for k in (
                "m_faces", "m_nodes", "delta_faces", "delta_nodes", "B_faces", "B_nodes"))):
        raise ValueError("reduced cap requires finite metric and A,B>0 throughout the slice")
    budget.check(ncell*(quadrature_order+1))
    result["budget"] = budget
    return result


def reduced_cap_profile(metric, r, *, budget=None):
    """Return the exact stated reconstruction at at most 513 finite radii."""
    budget = _context(metric["budget"] if budget is None else budget)
    if budget is not metric["budget"]:
        raise ValueError("profile and metric must share the cumulative budget")
    radii = _vector(r, "r")
    if len(radii) > 513 or np.any(radii < 0) or np.any(radii > metric["radius"]):
        raise ValueError("profile queries must stay inside the finite cap")
    budget.check(len(radii))
    cells = np.minimum(np.searchsorted(metric["r_faces"], radii, side="right")-1,
                       len(metric["state"]["Ip"])-1)
    result = {key: np.empty_like(radii) for key in ("f", "Pi", "Phi", "p", "q", "chi")}
    for i in np.unique(cells):
        mask = cells == i
        for key, value in _fields(metric["reconstruction"], int(i), radii[mask]).items():
            result[key][mask] = value
    return result


def reduced_cap_rhs(metric, *, T_x, radius_x, incoming_w_plus, budget=None):
    """Integral ALE balances, including the common transported wall value.

    qflux=T_x*B*p+R_x*y*q, pflux=T_x*B*q+R_x*y*p.  Axis qflux is zero;
    consequently sum(Iq_x) is the wall chi rate.  The caller must check its
    conversion against the same-stage hybrid c_s.  No value is imposed in
    addition to the one incoming wave.
    """
    budget = _context(metric["budget"] if budget is None else budget)
    if budget is not metric["budget"]:
        raise ValueError("RHS and metric must share the cumulative budget")
    Tx = _scalar(T_x, "T_x", positive=True)
    Rx = _scalar(radius_x, "radius_x")
    wpW = _scalar(incoming_w_plus, "incoming_w_plus")
    y, R, alpha = metric["y_faces"], metric["radius"], metric["alpha"]
    r, B, rec = metric["r_faces"], metric["B_faces"], metric["reconstruction"]
    ncell = len(y)-1
    if abs(Rx) >= Tx*B[-1]:
        raise ValueError("require a strictly timelike wall: abs(radius_x)<T_x*B_wall")
    budget.check(3*ncell+1)
    speed_plus, speed_minus = -Tx*B-Rx*y, Tx*B-Rx*y
    plus, minus = np.empty(ncell+1), np.empty(ncell+1)
    for j in range(1, ncell):
        def trace(i, sign):
            if i == 0:
                values = _fields(rec, i, r[j])
                return values["p"]+sign*values["q"]
            name = "plus" if sign == 1 else "minus"
            return rec[name][i]+rec[name+"_slope"][i]*(r[j]-rec["midpoint"][i])
        plus[j] = trace(j-1 if speed_plus[j] >= 0 else j, 1)
        minus[j] = trace(j-1 if speed_minus[j] >= 0 else j, -1)
    plus[0], minus[0] = rec["Q0"], -rec["Q0"]
    plus[-1] = R*wpW+metric["wall"]["f"]
    minus[-1] = metric["wall"]["a_minus"]
    plus_flux = (Tx*B+Rx*y)*plus
    minus_flux = (Rx*y-Tx*B)*minus
    pflux, qflux = (plus_flux+minus_flux)/2, (plus_flux-minus_flux)/2
    qflux[0] = 0.
    pface, qface = (plus+minus)/2, (plus-minus)/2
    rn, f = metric["r_nodes"], metric["f_nodes"]
    reaction = alpha*np.exp(metric["delta_nodes"])*rn*(
        (2*metric["n_nodes"]/rn**3-metric["potential_difference_over_alpha_nodes"])*f
        +_potential_reaction_quotient(f, alpha))
    integral = np.sum(metric["weights_nodes"]*reaction, axis=1)
    rates = {"Ip": np.diff(pflux)-Tx*integral, "Iq": np.diff(qflux)}
    chi_x = float(np.sum(rates["Iq"]))
    fW = metric["wall"]["f"]
    fW_x = (chi_x-fW*Rx)/R
    uW = metric["wall"]["c"]/R
    sqrtD = np.sqrt(1+alpha*uW*uW/6)
    cW_x = uW*Rx+R*sqrtD*fW_x
    if not all(np.all(np.isfinite(value)) for value in (*rates.values(), reaction)):
        raise ValueError("nonfinite reduced-cap ALE rate")
    budget.check()
    return {"rates": rates, "p_faces": pface, "q_faces": qface,
            "a_plus_faces": plus, "a_minus_faces": minus,
            "p_flux": pflux, "q_flux": qflux, "reaction_nodes": reaction,
            "reaction_integral": integral, "speed_plus": speed_plus, "speed_minus": speed_minus,
            "wall_chi_x": chi_x, "wall_f_x": float(fW_x), "wall_c_x": float(cW_x),
            "constructed_wall_flux_defect": chi_x-float(qflux[-1]),
            "charged_evaluations": 3*ncell+1,
            "characteristic_courant_rate": float(np.max(np.maximum(abs(speed_plus), abs(speed_minus)))/(R*np.diff(y)[0]))}


def reduced_cap_diagnostics(metric, *, budget=None):
    """Independent endpoint-trapezoid defects and labelled reconstruction identities."""
    budget = _context(metric["budget"] if budget is None else budget)
    if budget is not metric["budget"]:
        raise ValueError("diagnostics and metric must share the cumulative budget")
    r, rec, alpha = metric["r_faces"], metric["reconstruction"], metric["alpha"]
    ncell = len(r)-1
    mass, lapse, integral_p, integral_q = (np.empty((ncell, 2)) for _ in range(4))
    for i in range(ncell):
        ends = r[i:i+2]
        data = _fields(rec, i, ends)
        W = data["Pi"]**2+data["Phi"]**2
        mass[i] = ends**2*(_potential(data["f"], alpha)[0]/2+alpha*metric["A_faces"][i:i+2]*W/4)
        lapse[i] = alpha*ends*W/2
        if i == 0:
            integral_p[i] = rec["P0"]*ends**2/2+rec["P2"]*ends**4/4
        else:
            integral_p[i] = rec["pc"][i]*ends+rec["ps"][i]*ends**2/2
        integral_q[i] = data["chi"]
    mass_defect = np.diff(metric["m_faces"])/np.diff(r)-np.mean(mass, axis=1)
    lapse_defect = np.diff(metric["delta_faces"])/np.diff(r)-np.mean(lapse, axis=1)
    p_defect = np.diff(integral_p, axis=1)[:, 0]-metric["state"]["Ip"]
    q_defect = np.diff(integral_q, axis=1)[:, 0]-metric["state"]["Iq"]
    budget.check(3*ncell)
    return {"mass_trapezoid_defect": mass_defect, "lapse_trapezoid_defect": lapse_defect,
            "constructed_p_integral_defect": p_defect, "constructed_q_integral_defect": q_defect,
            "max_abs_mass_trapezoid_defect": float(np.max(abs(mass_defect))),
            "max_abs_lapse_trapezoid_defect": float(np.max(abs(lapse_defect))),
            "max_abs_constructed_p_integral_defect": float(np.max(abs(p_defect))),
            "max_abs_constructed_q_integral_defect": float(np.max(abs(q_defect))),
            "outgoing_reduced_wall_characteristic": metric["wall"]["a_minus"]}
