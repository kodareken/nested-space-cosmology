"""Exact local wall transmission and regular-centre observer conversion.

The fixed conformal model has b=2, h_minus=1 and h_plus=0.  The invariant
field c=chi/epsilon is retained at alpha=0.  These are algebraic helpers:
they choose no amplitude or gauge and solve no regional or wall evolution.
"""
from collections.abc import Mapping

import numpy as np


def _fields(value, names, label):
    if not isinstance(value, Mapping) or set(value) != set(names):
        raise ValueError(label + " must contain exactly " + ", ".join(names))
    return value


def _arrays(**values):
    arrays = []
    for name, value in values.items():
        try:
            if isinstance(value, (list, tuple)) and any(
                isinstance(item, (bool, np.bool_))
                for item in np.asarray(value, dtype=object).flat
            ):
                raise ValueError
            raw = np.asarray(value)
            if raw.dtype.kind not in "iuf" or raw.size == 0:
                raise ValueError
            array = np.asarray(raw, dtype=float)
            if not np.all(np.isfinite(array)):
                raise ValueError
        except (ValueError, TypeError, OverflowError) as error:
            raise ValueError(name + " must contain finite real numbers, not booleans") from error
        arrays.append(array)
    try:
        return dict(zip(values, np.broadcast_arrays(*arrays)))
    except ValueError as error:
        raise ValueError("input shapes are not broadcast compatible") from error


def _finite_result(result):
    if not all(np.all(np.isfinite(value)) for value in result.values()):
        raise ValueError("derived boundary quantities must be finite")
    return result


def wall_boundary_operator(state, outgoing, gauge, *, alpha):
    """Map outgoing first jets to incoming jets and wall state derivatives.

    All numeric inputs broadcast together. ``state`` contains r, c,
    eta_plus, eta_minus. ``outgoing`` has ``plus`` {r_V, c_V, log_C_V}
    and ``minus`` {r_Qminus, c_Qminus, log_C_Qminus}. These are raw
    coordinate derivatives supplied by the regional mixed equations.
    Each ``gauge`` side contains C_wall and log_C_s. The latter must be
    the actual common-proper-time derivative of that side's log(C_wall),
    not an independent gauge function. Compatibility with anchored source
    and centre coordinates belongs to the caller.

    s is dimensionless Einstein wall proper time. Both null label rates
    are positive and n=(-Udot,Vdot) points from cap to exterior. The cap
    labels are (Qminus,Qplus). Angular/temporal Israel, normal force and
    homogeneous c transmission are imposed independently. This supplies
    local data, not existence or constraint closure for a coupled solver.
    psi_s_normalized and N_*_normalized mean dimensionless Psi_s/sqrt(alpha)
    and n(Psi)/sqrt(alpha), evaluated by their regular continuous formulas
    even at alpha=0. The caller also owns the physical expanding domain.
    """
    state = _fields(state, ("r", "c", "eta_plus", "eta_minus"), "state")
    outgoing = _fields(outgoing, ("plus", "minus"), "outgoing")
    gauge = _fields(gauge, ("plus", "minus"), "gauge")
    op = _fields(outgoing["plus"], ("r_V", "c_V", "log_C_V"), "outgoing.plus")
    om = _fields(outgoing["minus"], ("r_Qminus", "c_Qminus", "log_C_Qminus"), "outgoing.minus")
    gp = _fields(gauge["plus"], ("C_wall", "log_C_s"), "gauge.plus")
    gm = _fields(gauge["minus"], ("C_wall", "log_C_s"), "gauge.minus")
    data = _arrays(**state, alpha=alpha, C_plus=gp["C_wall"], C_minus=gm["C_wall"],
                   Gamma_plus=gp["log_C_s"], Gamma_minus=gm["log_C_s"],
                   r_V=op["r_V"], c_V=op["c_V"], L_V=op["log_C_V"],
                   r_Qminus=om["r_Qminus"], c_Qminus=om["c_Qminus"],
                   L_Qminus=om["log_C_Qminus"])
    if (np.any(data["r"] <= 0) or np.any(data["C_plus"] <= 0)
            or np.any(data["C_minus"] <= 0) or np.any(data["alpha"] < 0)):
        raise ValueError("require r>0, both wall lapses>0 and alpha>=0")
    try:
        with np.errstate(over="raise", divide="raise", invalid="raise", under="ignore"):
            r, c, alpha = data["r"], data["c"], data["alpha"]
            Up = np.exp(-data["eta_plus"])/np.sqrt(data["C_plus"])
            Vp = np.exp(data["eta_plus"])/np.sqrt(data["C_plus"])
            Qm = np.exp(-data["eta_minus"])/np.sqrt(data["C_minus"])
            Qp = np.exp(data["eta_minus"])/np.sqrt(data["C_minus"])
            if any(np.any(rate <= 0) for rate in (Up, Vp, Qm, Qp)):
                raise ValueError("derived null label rates must remain positive")
            Rp, Rm = 2*Vp*data["r_V"], 2*Qm*data["r_Qminus"]
            cp, cm = 2*Vp*data["c_V"], 2*Qm*data["c_Qminus"]
            u = c/r
            D = 1+alpha*u*u/6
            root_D = np.sqrt(D)
            S = 2*D*root_D
            v, beta_bar = (Rp+Rm+S*r)/2, (Rp-Rm)/2
            cs, nc = (cp+cm)/2, (cp-cm)/2
            beta_plus, beta_minus = beta_bar-S*r/2, beta_bar+S*r/2
            k_bar = -3*root_D/4-alpha*u*(nc-u*beta_bar)/(2*r*D)-2*beta_bar/r
            k_plus, k_minus = k_bar-S/2, k_bar+S/2
            ell_plus, ell_minus = 2*Vp*data["L_V"], 2*Qm*data["L_Qminus"]
            nL_plus = ell_plus-data["Gamma_plus"]
            nL_minus = data["Gamma_minus"]-ell_minus
            result = {
                "r_s": v, "c_s": cs, "n_c": nc, "D": D, "F": 1/D, "S": S,
                "beta_bar": beta_bar, "beta_plus": beta_plus, "beta_minus": beta_minus,
                "k_bar": k_bar, "k_plus": k_plus, "k_minus": k_minus,
                "Udot_plus": Up, "Vdot_plus": Vp,
                "Qminus_dot_minus": Qm, "Qplus_dot_minus": Qp,
                "eta_plus_s": k_plus-nL_plus/2, "eta_minus_s": k_minus-nL_minus/2,
                "r_U_plus": (Rm+S*r)/(2*Up), "c_U_plus": cm/(2*Up),
                "r_Qplus_minus": (Rp+S*r)/(2*Qp), "c_Qplus_minus": cp/(2*Qp),
                "log_C_U_plus": (2*data["Gamma_plus"]-ell_plus)/(2*Up),
                "log_C_Qplus_minus": (2*data["Gamma_minus"]-ell_minus)/(2*Qp),
                "psi_s_normalized": (cs-u*v)/(r*root_D),
                "N_plus_normalized": (nc-u*beta_plus)/(r*root_D),
                "N_minus_normalized": (nc-u*beta_minus)/(r*root_D),
            }
    except FloatingPointError as error:
        raise ValueError("boundary arithmetic overflow or invalid derived quantity") from error
    return _finite_result(result)


def central_readout(*, alpha, C, q_minus_rate, q_plus_rate, c_Qplus):
    """Convert regular-centre data in arbitrary null labels to Jordan time.

    Rates are dQminus/dlambda and dQplus/dlambda along the centre, both
    positive. Regularity, c=0, the tangent cancellation and compatible
    radius derivatives are caller assumptions. The returned T_rate is
    dT_J/dlambda, p is the normalized Jordan scalar, and F is its conformal
    factor. The caller integrates and anchors the observer clock. This
    helper does not impose Qminus=Qplus=T_J or solve centre evolution.
    """
    data = _arrays(alpha=alpha, C=C, q_minus_rate=q_minus_rate,
                   q_plus_rate=q_plus_rate, c_Qplus=c_Qplus)
    if (np.any(data["alpha"] < 0) or np.any(data["C"] <= 0)
            or np.any(data["q_minus_rate"] <= 0) or np.any(data["q_plus_rate"] <= 0)):
        raise ValueError("require alpha>=0, C>0 and positive centre null label rates")
    try:
        with np.errstate(over="raise", divide="raise", invalid="raise", under="ignore"):
            base = data["C"]*data["q_minus_rate"]*data["q_plus_rate"]
            flux = data["q_plus_rate"]*data["c_Qplus"]
            T_squared = base+2*data["alpha"]*flux*flux/3
            T_rate = np.sqrt(T_squared)
            result = {"T_rate": T_rate, "p": 2*(flux/T_rate), "F": base/T_squared}
            if np.any(T_rate <= 0) or np.any(result["F"] <= 0):
                raise ValueError("derived centre clock and conformal factor must remain positive")
    except FloatingPointError as error:
        raise ValueError("centre arithmetic overflow or invalid derived quantity") from error
    return _finite_result(result)
