"""Conservative cap departures around an exact characteristic background.

The background is a numerical representation choice. The total wall scalar
is its primitive plus sum(dIq), and the incoming characteristic correction
is supplied by the coupled wall. No finite wall waveform is prescribed here.
Radial constraints use the total field, including every energy cross term.
"""
from functools import lru_cache

import numpy as np
from numpy.polynomial import Polynomial
from numpy.polynomial.legendre import leggauss

from .nsc_atlas_characteristic_initialization import Budget
from .nsc_atlas_characteristic_background import CharacteristicProfile
from .nsc_atlas_polar_cap import _potential, _scalar, _vector
from .nsc_atlas_reduced_cap import (
    _inputs, _reconstruct, _fields, _wall_conversion, _potential_reaction_quotient,
)

FIELD_KEYS = ("chi", "p", "q", "f", "Pi", "Phi")


def _budget(budget, metric=None):
    if not isinstance(budget, Budget) or (metric is not None and budget is not metric["budget"]):
        raise ValueError("one shared caller Budget is required")
    budget.check()
    return budget


@lru_cache(maxsize=3)
def _gauss_rule(order):
    if type(order) is not int or order not in (4, 8, 16):
        raise ValueError("quadrature_order must be 4, 8 or 16")
    nodes, weights = leggauss(order)
    basis = []
    for k, point in enumerate(nodes):
        poly = Polynomial([1.])
        for j, other in enumerate(nodes):
            if k != j:
                poly *= Polynomial([-other, 1.])/(point-other)
        primitive = poly.integ()
        basis.append(primitive-primitive(-1.))
    return nodes, weights, tuple(basis)


def _total(profile, T, rec, cell, r, budget):
    """One specified side of the departure reconstruction at radial points."""
    bg = profile.background(T, np.atleast_1d(r), budget=budget)
    dep = _fields(rec, cell, np.atleast_1d(r))
    budget.check(np.size(r))
    return {key: bg[key]+dep[key] for key in FIELD_KEYS}


def background_cap_metric(y_faces, radius, state, *, T, profile, alpha,
                          quadrature_order=4, budget):
    """Integrate the radial constraints on a total background+departure slice.

    ``state`` contains departures in the unweighted physical integrals Ip,Iq.
    Integration panels split at cells and all reflected/incoming profile
    knots. A Gauss antiderivative approximates the local integrating factor
    inside positive-weight partial mass/lapse quadrature; refinement must
    check their consistency independently. There is no adaptive
    retry. Returned node arrays have shape (panels, quadrature_order).
    """
    budget = _budget(budget)
    before = budget.evaluations
    if not isinstance(profile, CharacteristicProfile):
        raise ValueError("require one CharacteristicProfile")
    T = _scalar(T, "T")
    y, R, state, alpha = _inputs(y_faces, radius, state, alpha)
    if R >= 1:
        raise ValueError("characteristic background requires radius below 1")
    nodes, weights, basis = _gauss_rule(quadrature_order)
    r = R*y
    rec = _reconstruct(r, state)
    cuts = np.unique(np.r_[r, profile.radial_cuts(T, R, budget=budget)])
    if len(cuts) > 1154:
        raise ValueError("cap exceeds the finite panel allowance")
    widths, mids = np.diff(cuts), (cuts[:-1]+cuts[1:])/2
    owners = np.minimum(np.searchsorted(r, mids, side="right")-1, len(y)-2)
    rn = mids[:, None]+widths[:, None]*nodes/2
    wn = widths[:, None]*weights/2
    fields = {key: np.empty_like(rn) for key in FIELD_KEYS}
    gamma_edges, n_edges = np.zeros(len(cuts)), np.zeros(len(cuts))
    gamma_nodes, n_nodes = np.empty_like(rn), np.empty_like(rn)
    gamma_interpolant_defect = np.empty_like(rn)
    center = _total(profile, T, rec, 0, np.array([0.]), budget)
    f0, Pi0 = float(center["f"][0]), float(center["Pi"][0])
    z0 = np.sqrt(alpha/6)*f0
    delta0 = -float(np.log1p(2*np.sinh(z0/2)**2))
    with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
        for j, cell in enumerate(owners):
            local = _total(profile, T, rec, int(cell), rn[j], budget)
            for key in fields:
                fields[key][j] = local[key]
            W = local["Pi"]**2+local["Phi"]**2
            a = alpha*rn[j]*W/2
            h = widths[j]/2
            targets = np.r_[rn[j], cuts[j+1]]
            half = (targets-cuts[j])/2
            partial_r = cuts[j]+half[:,None]*(1+nodes)
            data = _total(profile,T,rec,int(cell),partial_r.ravel(),budget)
            Wp = (data["Pi"]**2+data["Phi"]**2).reshape(partial_r.shape)
            ap = alpha*partial_r*Wp/2
            bp = partial_r**2*(_potential(data["f"],alpha)[1].reshape(partial_r.shape)/2
                              +(1-partial_r**2)*Wp/4)
            # Only the exponential's smooth internal primitive is interpolated.
            # Actual positive-weight partial integrals avoid negative nodal mass
            # from the antiderivative matrix on a sharply localized source.
            theta = (partial_r-mids[j])/h
            Gpoly = sum(h*a[k]*poly(theta) for k,poly in enumerate(basis))
            Gtarget = half*(ap @ weights)
            responses = np.exp(-Gtarget)*n_edges[j]+half*(
                (np.exp(Gpoly-Gtarget[:,None])*bp) @ weights)
            gamma_nodes[j] = gamma_edges[j]+Gtarget[:-1]
            n_nodes[j] = responses[:-1]
            gamma_edges[j+1] = gamma_edges[j]+Gtarget[-1]
            n_edges[j+1] = responses[-1]
            gamma_interpolant_defect[j] = sum(h*a[k]*poly(nodes)
                for k,poly in enumerate(basis))-Gtarget[:-1]
            budget.check(2*partial_r.size+len(nodes))
        face_index = np.searchsorted(cuts, r)
        if not np.array_equal(cuts[face_index], r):
            raise ValueError("integration panels omitted a cell face")
        nf, gf = n_edges[face_index], gamma_edges[face_index]
        Af = 1-r*r
        Af[1:] -= 2*alpha*nf[1:]/r[1:]
        An = 1-rn*rn-2*alpha*n_nodes/rn
        df, dn = delta0+gf, delta0+gamma_nodes
        Bf, Bn = Af*np.exp(df), An*np.exp(dn)
        bg_faces = profile.background(T, r, budget=budget)
        chi = bg_faces["chi"]+rec["chi"]
        fW = float(chi[-1]/R)
        dep_minus = rec["minus"][-1]+rec["minus_slope"][-1]*(R-rec["midpoint"][-1])
        amW = float(bg_faces["p"][-1]-bg_faces["q"][-1]+dep_minus)
        potential, pdiff, derivative = _potential(fields["f"], alpha)
        baseline_integrals = profile.integrals(T, r, budget=budget)
    if (np.any(Af <= 0) or np.any(An <= 0) or np.any(Bf <= 0) or np.any(Bn <= 0)
            or any(not np.all(np.isfinite(a)) for a in
                   (*fields.values(), nf, n_nodes, df, dn, Bf, Bn))):
        raise ValueError("background cap requires finite metric and A,B>0")
    budget.check(len(r))
    return {"y_faces": y, "r_faces": r, "radius": R, "T": T, "alpha": alpha,
            "state": state, "reconstruction": rec, "profile": profile,
            "quadrature_order": quadrature_order, "panel_cuts": cuts, "panel_owners": owners,
            "r_nodes": rn, "weights_nodes": wn, "n_faces": nf, "n_nodes": n_nodes,
            "m_faces": r**3/2+alpha*nf, "m_nodes": rn**3/2+alpha*n_nodes,
            "A_faces": Af, "A_nodes": An, "delta_faces": df, "delta_nodes": dn,
            "B_faces": Bf, "B_nodes": Bn, "chi_faces": chi,
            "gamma_interpolant_defect_nodes": gamma_interpolant_defect,
            "background_faces": bg_faces,
            "physical_integrals": {key: state[key]+baseline_integrals[key] for key in state},
            "potential_nodes": potential, "potential_difference_over_alpha_nodes": pdiff,
            "potential_derivative_normalized_nodes": derivative,
            "center": {"f": f0, "Pi": Pi0, "Phi": float(center["Phi"][0]), "delta": delta0},
            "wall": {"A": float(Af[-1]), "delta": float(df[-1]), "B": float(Bf[-1]),
                     "a_minus": amW, "departure_a_minus": float(dep_minus),
                     "w_minus": (amW+fW)/R, "chi": float(chi[-1]),
                     "f": fW, "c": R*_wall_conversion(fW, alpha)},
            **{key+"_nodes": value for key, value in fields.items()},
            "budget": budget, "charged_evaluations": budget.evaluations-before}


def background_cap_profile(metric, r, *, budget):
    """Evaluate the actual declared total field, not a cell-mean surrogate."""
    _budget(budget, metric)
    radii = _vector(r, "r")
    if len(radii) > 513 or np.any(radii < 0) or np.any(radii > metric["radius"]):
        raise ValueError("query at most 513 radii inside the cap")
    cells = np.minimum(np.searchsorted(metric["r_faces"], radii, side="right")-1,
                       len(metric["y_faces"])-2)
    result = {key: np.empty_like(radii) for key in FIELD_KEYS}
    for cell in np.unique(cells):
        mask = cells == cell
        row = _total(metric["profile"], metric["T"], metric["reconstruction"],
                     int(cell), radii[mask], budget)
        for key in result:
            result[key][mask] = row[key]
    return result


def background_cap_rhs(metric, *, T_x, radius_x, incoming_w_plus, budget):
    """ALE balances for departures, retaining the full incoming discrepancy."""
    _budget(budget, metric)
    before = budget.evaluations
    Tx = _scalar(T_x, "T_x", positive=True)
    Rx = _scalar(radius_x, "radius_x")
    wp = _scalar(incoming_w_plus, "incoming_w_plus")
    r, y, R, alpha = (metric[key] for key in ("r_faces", "y_faces", "radius", "alpha"))
    B, rec, bg = metric["B_faces"], metric["reconstruction"], metric["background_faces"]
    if abs(Rx) >= Tx*B[-1]:
        raise ValueError("require a strictly timelike wall")
    n = len(y)-1
    speed_plus, speed_minus = -Tx*B-Rx*y, Tx*B-Rx*y
    plus, minus = np.empty(n+1), np.empty(n+1)
    for j in range(1, n):
        for sign, speed, output in ((1, speed_plus, plus), (-1, speed_minus, minus)):
            cell = j-1 if speed[j] >= 0 else j
            row = _fields(rec, cell, np.array([r[j]]))
            output[j] = float(row["p"][0]+sign*row["q"][0])
    plus[0], minus[0] = rec["Q0"], -rec["Q0"]
    plus[-1] = R*wp+metric["wall"]["f"]-bg["p"][-1]-bg["q"][-1]
    minus[-1] = metric["wall"]["departure_a_minus"]
    dp, dq = (plus+minus)/2, (plus-minus)/2
    # B-A0: retain all finite metric effects on the known background wave.
    dB = (1-r*r)*np.expm1(metric["delta_faces"])
    dB[1:] -= 2*alpha*metric["n_faces"][1:]*np.exp(metric["delta_faces"][1:])/r[1:]
    pf = Tx*(B*dq+dB*bg["q"])+Rx*y*dp
    qf = Tx*(B*dp+dB*bg["p"])+Rx*y*dq
    qf[0] = 0.
    rn, f = metric["r_nodes"], metric["f_nodes"]
    reaction = alpha*np.exp(metric["delta_nodes"])*rn*(
        (2*metric["n_nodes"]/rn**3-metric["potential_difference_over_alpha_nodes"])*f
        +_potential_reaction_quotient(f, alpha))
    integral = np.bincount(metric["panel_owners"],
        weights=np.sum(metric["weights_nodes"]*reaction, axis=1), minlength=n)
    rates = {"Ip": np.diff(pf)-Tx*integral, "Iq": np.diff(qf)}
    baseline_chi_x = Tx*(1-R*R)*bg["p"][-1]+Rx*bg["q"][-1]
    chi_x = float(baseline_chi_x+np.sum(rates["Iq"]))
    fW_x = (chi_x-metric["wall"]["f"]*Rx)/R
    uW = metric["wall"]["c"]/R
    cW_x = uW*Rx+R*np.sqrt(1+alpha*uW*uW/6)*fW_x
    total_p, total_q = bg["p"]+dp, bg["q"]+dq
    total_wall_flux = Tx*B[-1]*total_p[-1]+Rx*total_q[-1]
    if any(not np.all(np.isfinite(v)) for v in (*rates.values(), reaction, chi_x, cW_x)):
        raise ValueError("nonfinite background-cap balance")
    budget.check(3*n+1+rn.size)
    return {"rates": rates, "departure_p_faces": dp, "departure_q_faces": dq,
            "p_faces": total_p, "q_faces": total_q, "p_flux": pf, "q_flux": qf,
            "reaction_nodes": reaction, "reaction_integral": integral,
            "wall_chi_x": chi_x, "wall_f_x": float(fW_x), "wall_c_x": float(cW_x),
            "baseline_wall_chi_x": float(baseline_chi_x),
            "constructed_wall_flux_defect": float(chi_x-total_wall_flux),
            "incoming_departure": float(plus[-1]),
            "characteristic_courant_rate": float(np.max(np.maximum(abs(speed_plus),abs(speed_minus)))/(R*np.diff(y)[0])),
            "charged_evaluations": budget.evaluations-before}
