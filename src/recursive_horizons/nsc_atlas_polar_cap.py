"""Regular polar-cap constraint and ALE finite-volume stage operators.

This is a first-order coupled-solver prototype, not a time evolution or a
stability theorem.  The uniform coordinate is y=r/R, not rho=r**2.  f is an
ordinary radial cell average, Pi an r**2-weighted average, and Q=Phi/r an
r-weighted average.  f and Q are constant in each reconstruction cell;
Pi is limited linear in rho about its exact weighted mean.  That explicit
even reconstruction is needed for consistency of Q's first axis cell.

The cap metric is determined by its radial constraints and Jordan central
clock.  The wall supplies one incoming wave and its transported common
field value, not an independently chosen second waveform or cap lapse.
"""
from collections.abc import Mapping

import numpy as np
from numpy.polynomial.legendre import leggauss

from .nsc_atlas_characteristic_initialization import Budget

MAX_CELLS = 128
QUADRATURE_ORDERS = (2, 4, 8)


def _scalar(value, name, *, positive=False, nonnegative=False):
    raw = np.asarray(value)
    if raw.shape or raw.dtype.kind not in "iuf":
        raise ValueError(name + " must be a finite real scalar")
    result = float(raw)
    if (not np.isfinite(result) or (positive and result <= 0)
            or (nonnegative and result < 0)):
        raise ValueError(name + " is outside its finite permitted range")
    return result


def _vector(value, name):
    if isinstance(value, (list, tuple)) and any(
            isinstance(item, (bool, np.bool_)) for item in np.asarray(value, dtype=object).flat):
        raise ValueError(name + " must contain real numbers, not booleans")
    raw = np.asarray(value)
    if raw.ndim != 1 or raw.dtype.kind not in "iuf":
        raise ValueError(name + " must be a finite real vector")
    result = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(result)):
        raise ValueError(name + " must be finite")
    return result


def _context(budget):
    if budget is None:
        budget = Budget(maximum_evaluations=20000, cpu_seconds=30)
    if not callable(getattr(budget, "check", None)):
        raise ValueError("budget must provide a cumulative check(count) method")
    budget.check()
    return budget


def _inputs(y_faces, radius, state, alpha):
    y = _vector(y_faces, "y_faces")
    n = len(y)-1
    if (not 3 <= n <= MAX_CELLS or y[0] != 0 or y[-1] != 1
            or np.any(np.diff(y) <= 0)
            or not np.allclose(np.diff(y), 1/n, rtol=0, atol=8*np.finfo(float).eps)):
        raise ValueError("require 3..128 uniform cells with y endpoints 0 and 1")
    radius = _scalar(radius, "radius", positive=True)
    alpha = _scalar(alpha, "alpha", nonnegative=True)
    if not isinstance(state, Mapping) or set(state) != {"f", "Pi", "Q"}:
        raise ValueError("state must contain exactly f, Pi, Q")
    fields = {key: _vector(state[key], key).copy() for key in ("f", "Pi", "Q")}
    if any(len(value) != n for value in fields.values()):
        raise ValueError("each state field must have one value per cell")
    return y.copy(), radius, fields, alpha


def _potential(f, alpha):
    q = np.sqrt(alpha/6)*f
    # np.sinc has no hyperbolic counterpart.  A series avoids 0/0 and keeps
    # the continuous alpha=0 potential derivative and mass response.
    small = abs(q) < 1e-4
    sinhc = np.ones_like(q)
    sinhc[small] = 1+q[small]**2/6+q[small]**4/120
    np.divide(np.sinh(q), q, out=sinhc, where=~small)
    ch = np.cosh(q)
    diff = f*f*sinhc*sinhc*(ch*ch+1)/2
    return 3+alpha*diff, diff, 2*f*sinhc*ch**3


def _pi_reconstruction(r, Pi):
    rho_mean = 3*(r[1:]**5-r[:-1]**5)/(5*(r[1:]**3-r[:-1]**3))
    secants = np.diff(Pi)/np.diff(rho_mean)
    slopes = np.empty_like(Pi)
    slopes[0], slopes[-1] = secants[0], secants[-1]
    left, right = secants[:-1], secants[1:]
    slopes[1:-1] = np.where(left*right > 0,
                           np.sign(left)*np.minimum(abs(left), abs(right)), 0)
    return rho_mean, slopes, Pi-slopes*rho_mean


def cap_metric(y_faces, radius, state, *, alpha, quadrature_order=4, budget=None):
    """Integrate the cap radial constraints for one explicitly reconstructed slice.

    Return metric arrays at faces and open Gauss nodes, the copied state,
    reconstruction coefficients, regular centre limits and the outgoing
    wall trace.  Face arrays include r=0 without dividing by it.  The
    vacuum is subtracted analytically: m=r**3/2+alpha*n, where
    n'+alpha*r*W*n/2=r**2*((P-3)/(2alpha)+(1-r**2)*W/4).
    No radius, source or lapse is chosen here; the caller supplies R and
    the field.  An A<=0 or nonfinite result stops the stage.

    Each constraint integrand point is charged to the supplied cumulative
    budget, as are reconstructed field/metric nodes.  A new 30 CPU-second,
    20,000-evaluation budget is made only when none is supplied.
    """
    budget = _context(budget)
    y, R, state, alpha = _inputs(y_faces, radius, state, alpha)
    if (isinstance(quadrature_order, (bool, np.bool_))
            or not isinstance(quadrature_order, (int, np.integer))
            or quadrature_order not in QUADRATURE_ORDERS):
        raise ValueError("quadrature_order must be 2, 4 or 8")
    ncell = len(y)-1
    r = R*y
    rho_mean, slopes, intercept = _pi_reconstruction(r, state["Pi"])
    nodes, weights = leggauss(quadrature_order)
    rn = (r[:-1, None]+r[1:, None])/2+np.diff(r)[:, None]*nodes/2
    wn = np.diff(r)[:, None]*weights/2
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
            P, Pdiff, U = _potential(state["f"], alpha)
            q0 = np.sqrt(alpha/6)*state["f"][0]
            delta0 = -float(np.log1p(2*np.sinh(q0/2)**2))
            tanhc0 = (1-q0*q0/3+2*q0**4/15 if abs(q0) < 1e-4
                      else np.tanh(q0)/q0)
            If = np.zeros(ncell+1)
            nf = np.zeros(ncell+1)
            In, nn = np.empty_like(rn), np.empty_like(rn)
            for i in range(ncell):
                a, b = r[i:i+2]
                p, slope, Q = intercept[i], slopes[i], state["Q"][i]

                def primitive(z):
                    return alpha*(p*p*(z*z-a*a)/4
                                  +(2*p*slope+Q*Q)*(z**4-a**4)/8
                                  +slope*slope*(z**6-a**6)/12)

                # Advance an integrating factor over increasing subintervals;
                # exp(I_left-I_right) avoids a large positive cumulative exp.
                targets = np.r_[rn[i], b]
                previous, previous_I, response = a, If[i], nf[i]
                for j, end in enumerate(targets):
                    q = (previous+end)/2+(end-previous)*nodes/2
                    w = (end-previous)*weights/2
                    budget.check(len(q))
                    Iq, Ie = If[i]+primitive(q), If[i]+primitive(end)
                    W = (p+slope*q*q)**2+q*q*Q*Q
                    forcing = q*q*(Pdiff[i]/2+(1-q*q)*W/4)
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
            Pi_nodes = intercept[:, None]+slopes[:, None]*rn*rn
            Phi_nodes = rn*state["Q"][:, None]
            Pi_axis = intercept[0]
            Pi_wall = intercept[-1]+slopes[-1]*R*R
            mu0 = P[0]/6+alpha*Pi_axis**2/12
            result = {
                "y_faces": y, "radius": R, "alpha": alpha, "state": state,
                "quadrature_order": quadrature_order, "r_faces": r,
                "r_nodes": rn, "weights_nodes": wn,
                "m_faces": mf, "n_faces": nf, "A_faces": Af,
                "delta_faces": df, "B_faces": Bf,
                "m_nodes": mn, "n_nodes": nn, "A_nodes": An,
                "delta_nodes": dn, "B_nodes": Bn,
                "Pi_nodes": Pi_nodes, "Phi_nodes": Phi_nodes,
                "potential": P, "potential_difference_over_alpha": Pdiff,
                "potential_derivative_normalized": U,
                "Pi_rho_mean": rho_mean, "Pi_rho_slope": slopes,
                "Pi_intercept": intercept,
                "charged_evaluations": ncell*(quadrature_order+1)**2,
                "spatial_order": 1,
                "center": {"f": float(state["f"][0]), "Pi": float(Pi_axis),
                           "Q": float(state["Q"][0]), "mu": float(mu0),
                           "delta": delta0, "B": float(Bf[0]),
                           "p_J": float(state["f"][0]*tanhc0)},
                "wall": {"A": float(Af[-1]), "delta": float(df[-1]),
                         "B": float(Bf[-1]), "w_minus": float(Pi_wall-R*state["Q"][-1])},
            }
    except FloatingPointError as error:
        raise ValueError("nonfinite polar-cap constraint arithmetic") from error
    if (np.any(Af <= 0) or np.any(An <= 0)
            or any(not np.all(np.isfinite(result[k])) for k in (
                "m_faces", "m_nodes", "delta_faces", "delta_nodes", "B_faces", "B_nodes"))
            or np.any(Bf <= 0) or np.any(Bn <= 0)):
        raise ValueError("polar cap requires finite metric and A,B>0 throughout the slice")
    budget.check(ncell*(quadrature_order+1))
    result["budget"] = budget
    return result


def cap_rhs(metric, *, T_x, radius_x, incoming_w_plus, wall_f, budget=None):
    """Return one conservative ALE stage and its boundary/CFL diagnostics.

    x is the exterior incoming label: T_x=T_s/x_s and R_x=r_s/x_s
    come from the same hybrid wall event.  All stage fields and the
    supplied incoming w_plus must be simultaneous.  wall_f is the common
    dynamically transported canonical value; it enters the mesh flux only
    when the cap expands.  It is not an independent boundary waveform.

    The axis reflects Pi+Phi into Pi-Phi and retains the B*Pi flux in the
    Q equation.  Wall timelikeness supplies exactly one incoming wave.
    The reported characteristic Courant rate is a diagnostic, not a
    proven stable timestep including reaction/geometric terms.
    """
    if not isinstance(metric, Mapping) or "budget" not in metric:
        raise ValueError("metric must be a cap_metric result")
    budget = _context(metric["budget"] if budget is None else budget)
    if budget is not metric["budget"]:
        raise ValueError("metric and RHS must share the same cumulative budget")
    T_x = _scalar(T_x, "T_x", positive=True)
    Rx = _scalar(radius_x, "radius_x")
    wp_wall = _scalar(incoming_w_plus, "incoming_w_plus")
    f_wall = _scalar(wall_f, "wall_f")
    y, R, state, alpha = _inputs(metric["y_faces"], metric["radius"], metric["state"], metric["alpha"])
    r, B = metric["r_faces"], metric["B_faces"]
    ncell = len(y)-1
    if abs(Rx) >= T_x*B[-1]:
        raise ValueError("require a strictly timelike wall: abs(radius_x)<T_x*B_wall")
    budget.check(3*ncell+1)
    p, slope, Q = metric["Pi_intercept"], metric["Pi_rho_slope"], state["Q"]
    Pi_face, Phi_face = np.empty(ncell+1), np.empty(ncell+1)
    speed_plus, speed_minus = -T_x*B-Rx*y, T_x*B-Rx*y
    left_pi = p[:-1]+slope[:-1]*r[1:-1]**2
    right_pi = p[1:]+slope[1:]*r[1:-1]**2
    left_phi, right_phi = r[1:-1]*Q[:-1], r[1:-1]*Q[1:]
    wp = np.where(speed_plus[1:-1] >= 0, left_pi+left_phi, right_pi+right_phi)
    wm = np.where(speed_minus[1:-1] >= 0, left_pi-left_phi, right_pi-right_phi)
    Pi_face[1:-1], Phi_face[1:-1] = (wp+wm)/2, (wp-wm)/2
    Pi_face[0], Phi_face[0] = p[0], 0
    wm_wall = metric["wall"]["w_minus"]
    Pi_face[-1], Phi_face[-1] = (wp_wall+wm_wall)/2, (wp_wall-wm_wall)/2
    f_face = np.empty(ncell+1)
    f_face[0] = state["f"][0]
    f_face[1:-1] = state["f"][1:] if Rx > 0 else state["f"][:-1]
    f_face[-1] = f_wall if Rx > 0 else state["f"][-1]
    dy, dr = np.diff(y), np.diff(r)
    d1, d2 = np.diff(y*y)/2, np.diff(y**3)/3
    phi_flux = T_x*B*Pi_face+Rx*y*Phi_face
    pi_flux = y*y*(T_x*B*Phi_face+Rx*y*Pi_face)
    f_flux = Rx*y*f_face
    wn, rn = metric["weights_nodes"], metric["r_nodes"]
    source_f = np.sum(wn*metric["B_nodes"]*metric["Pi_nodes"], axis=1)/dr
    source_Pi = np.sum(wn*rn*rn*np.exp(metric["delta_nodes"]), axis=1)
    source_Pi *= metric["potential_derivative_normalized"]/(R**3*d2)
    rates = {
        "f": T_x*source_f+np.diff(f_flux)/(R*dy)-Rx*state["f"]/R,
        "Pi": np.diff(pi_flux)/(R*d2)-3*Rx*state["Pi"]/R-T_x*source_Pi,
        "Q": np.diff(phi_flux)/(R*R*d1)-2*Rx*state["Q"]/R,
    }
    if not all(np.all(np.isfinite(value)) for value in rates.values()):
        raise ValueError("nonfinite polar-cap ALE rate")
    budget.check()
    return {"rates": rates, "Pi_faces": Pi_face, "Phi_faces": Phi_face,
            "f_faces": f_face, "Pi_flux": pi_flux, "Phi_flux": phi_flux,
            "f_flux": f_flux, "source_f": source_f, "source_Pi": source_Pi,
            "speed_plus": speed_plus, "speed_minus": speed_minus,
            "charged_evaluations": 3*ncell+1,
            "characteristic_courant_rate": float(np.max(np.maximum(abs(speed_plus), abs(speed_minus)))/(R*dy[0])),
            "wall_value_cell_difference": float(f_wall-state["f"][-1]),
            "axis_Pi_rate_limit": float(T_x*(3*B[0]*state["Q"][0]
                                             -np.exp(metric["center"]["delta"])*metric["potential_derivative_normalized"][0]))}


def cap_diagnostics(metric, *, budget=None):
    """Independently differenced radial defects, not constraint certificates.

    Compare face differences with endpoint trapezoidal constraint densities
    using the stated within-cell reconstruction.  These densities were not
    used for the Gauss mass integration.  Phi-f_r uses ordinary differences
    of the saved f cell averages and Phi at their midpoint.  No mixed or
    temporal field equation is substituted for a derivative.
    """
    budget = _context(metric["budget"] if budget is None else budget)
    if budget is not metric["budget"]:
        raise ValueError("diagnostics must share the metric budget")
    r, state, alpha = metric["r_faces"], metric["state"], metric["alpha"]
    p, slope, Q = metric["Pi_intercept"], metric["Pi_rho_slope"], state["Q"]
    ends = np.column_stack((r[:-1], r[1:]))
    W = (p[:, None]+slope[:, None]*ends*ends)**2+ends*ends*Q[:, None]**2
    A = np.column_stack((metric["A_faces"][:-1], metric["A_faces"][1:]))
    mass_density = ends*ends*(metric["potential"][:, None]/2+alpha*A*W/4)
    lapse_density = alpha*ends*W/2
    mass_defect = np.diff(metric["m_faces"])/np.diff(r)-np.mean(mass_density, axis=1)
    lapse_defect = np.diff(metric["delta_faces"])/np.diff(r)-np.mean(lapse_density, axis=1)
    midpoint = (r[:-1]+r[1:])/2
    derivative = np.gradient(state["f"], midpoint, edge_order=2)
    compatibility = midpoint*Q-derivative
    budget.check(3*len(midpoint))
    return {"r_midpoint": midpoint, "mass_trapezoid_defect": mass_defect,
            "lapse_trapezoid_defect": lapse_defect,
            "Phi_minus_fd_f_r": compatibility,
            "max_abs_mass_trapezoid_defect": float(np.max(abs(mass_defect))),
            "max_abs_lapse_trapezoid_defect": float(np.max(abs(lapse_defect))),
            "max_abs_Phi_minus_fd_f_r": float(np.max(abs(compatibility)))}
