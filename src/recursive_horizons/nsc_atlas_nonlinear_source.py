"""Exact fixed source first jets in an affine exterior characteristic gauge.

This integrates the existing compact pulse on its initial leg only. The
source-only amplitude threshold is not a bound for the transmitted spacetime.
Coordinates are Ubar=CJ0*(U-U0), x with C_E=F on Ubar=0.
"""
import copy
import math
import time

import numpy as np
from scipy.integrate import solve_ivp

from .nsc_atlas_characteristic_initialization import (
    Budget, source_jet, validate_geometry,
)
from .nsc_atlas_field_preparation import radius_from_null

MAX_EVALUATIONS = 20000
CPU_SECONDS = 30
MAX_QUERY_POINTS = 513


def _scalar(value, label):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(label+" must be a finite real scalar")
    if not isinstance(value, (int, float, np.integer, np.floating)):
        raise ValueError(label+" must be a finite real scalar")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(label+" must be a finite real scalar")
    return value


def _check_budget(budget, count=0):
    if not isinstance(budget, Budget):
        raise ValueError("use the existing characteristic initialization Budget")
    for name in ("maximum_evaluations", "evaluations"):
        value = getattr(budget, name)
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            raise ValueError("budget counts must be nonnegative integers")
    if not (0 < budget.maximum_evaluations <= MAX_EVALUATIONS
            and 0 <= budget.evaluations <= budget.maximum_evaluations):
        raise ValueError("source budget must allow at most 20000 cumulative evaluations")
    cpu = _scalar(budget.cpu_seconds, "budget CPU seconds")
    started = _scalar(budget.started, "budget start time")
    if not 0 < cpu <= CPU_SECONDS or started > time.process_time():
        raise ValueError("source budget must use at most 30 CPU seconds and a past start")
    budget.check(count)


def source_amplitude_bound(geometry, *, budget=None):
    """Continuous pulse bounds and a conservative source-only threshold.

    The tube uses rJ>=1/5, abs(rJ_x)<=3/2 and F>=1/2. Strict alpha<alpha_star
    suffices for continuation of this source ODE. test_alpha=alpha_star/2
    is the admitted deterministic initialization control, not a calibration.
    """
    if budget is not None:
        _check_budget(budget, 3)  # Reference-radius checks in validate_geometry.
    try:
        validate_geometry(geometry)
        U0, Va, Vb, xb = (float(geometry[k]) for k in ("U0", "Va", "Vb", "x_last"))
    except (KeyError, TypeError, ValueError, OverflowError) as error:
        raise ValueError("invalid original source geometry") from error
    Z1 = -4*xb*math.exp(xb/.8)/(.8**2*U0*(Vb-Va))
    Z2 = -4*math.exp(xb/.8)*(1+xb/.8)/(.8**2*U0*(Vb-Va))
    A1 = 8/math.sqrt(15)*(14/15)**7
    A2 = (80/7)*(20/21)**6
    M1, M2 = A1*Z1, A2*Z1**2+A1*Z2
    Lbar, length = 5*M2/6+25*M1/2, xb-.4
    threshold = min(1/(2*Lbar*length), 3/25)
    if not all(math.isfinite(v) and v > 0 for v in (Z1, Z2, Lbar, length, threshold)):
        raise ValueError("invalid source continuation constants")
    return {"M0": 1., "M1_bound": M1, "M2_bound": M2,
            "Lbar": Lbar, "length": length, "alpha_star": threshold,
            "test_alpha": threshold/2, "r_min": .2, "v_max": 1.5, "F_min": .5}


def _query(value, lower, upper, label):
    try:
        raw = np.asarray(value)
        if raw.dtype.kind not in "iuf" or raw.ndim > 1 or not 1 <= raw.size <= MAX_QUERY_POINTS:
            raise ValueError
        if isinstance(value, (list, tuple)) and any(
                isinstance(item, (bool, np.bool_)) for item in np.asarray(value, dtype=object).flat):
            raise ValueError
        result = np.asarray(raw, dtype=float)
        if not np.all(np.isfinite(result)) or np.any(result < lower) or np.any(result > upper):
            raise ValueError
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError(label+" requires at most 513 finite points inside its closed domain") from error
    return result


class ExactSource:
    """Bounded DOP853 source integration with no amplitude/source selection.

    geometry must come from source_geometry; supply its Budget too to count
    geometry construction cumulatively. Dense queries never extrapolate.
    h=(rJ-x)/alpha and h_x are integrated without cancellation and remain
    the quadratic sensitivity at alpha=0. c and c_Ubar retain the normalized
    leading field there. log_C_* are derivatives of the Einstein lapse.
    The caller retains the budget; neither construction nor queries reset it.
    """

    def __init__(self, geometry, *, alpha, budget=None, rtol=2e-10, atol=1e-12,
                 max_step=None):
        self.budget = Budget(cpu_seconds=CPU_SECONDS) if budget is None else budget
        self.bounds = source_amplitude_bound(geometry, budget=self.budget)
        self.geometry = copy.deepcopy(geometry)
        self.alpha = _scalar(alpha, "alpha")
        self.rtol, self.atol = _scalar(rtol, "rtol"), _scalar(atol, "atol")
        if not 0 <= self.alpha < self.bounds["alpha_star"]:
            raise ValueError("require 0<=alpha<the source-only alpha_star")
        if not 1e-13 <= self.rtol <= 1e-6 or not 1e-15 <= self.atol <= 1e-8:
            raise ValueError("require 1e-13<=rtol<=1e-6 and 1e-15<=atol<=1e-8")
        length = self.bounds["length"]
        self.max_step = length/32 if max_step is None else _scalar(max_step, "max_step")
        if not 0 < self.max_step <= length/8:
            raise ValueError("require 0<max_step<=source length/8")
        self.CJ0 = -3.2/self.geometry["U0"]
        self.domain = {key: self.geometry[key] for key in ("U0", "Ua", "Ub", "Va", "Vb")}
        self.domain.update(x_a=.4, x_b=self.geometry["x_last"], Ubar_a=0.,
                           Ubar_b=self.CJ0*(self.geometry["Ua"]-self.geometry["U0"]), CJ0=self.CJ0)
        self.knots = np.array([.4, *self.geometry["support_x"], self.geometry["x_last"]])
        self.nfev = 0
        self._segments = []
        a0 = (.8/.4-1)/4
        initial = np.array([0., 0., a0, 0., -(1/.8+1/.4)*a0])
        for left, right in zip(self.knots[:-1], self.knots[1:]):
            solved = solve_ivp(self._rhs, (left, right), initial, method="DOP853",
                               rtol=self.rtol, atol=self.atol, max_step=self.max_step,
                               dense_output=True)
            if not solved.success or solved.sol is None or not np.all(np.isfinite(solved.y)):
                raise RuntimeError("exact source integration did not complete: "+solved.message)
            self._segments.append(solved.sol)
            initial = solved.y[:, -1]
        _check_budget(self.budget)

    def _profile(self, x):
        _check_budget(self.budget, int(np.size(x)))
        p, px, pxx = source_jet(x, self.geometry)
        quiet = (np.asarray(x) <= self.knots[1]) | (np.asarray(x) >= self.knots[2])
        g = np.where(quiet, 0., x*p)
        gx = np.where(quiet, 0., p+x*px)
        gxx = np.where(quiet, 0., 2*px+x*pxx)
        return g, gx, gxx

    def _quantities(self, x, state, g, gx, gxx):
        h, w, a, bc, lam = state
        rJ, rJx = x+self.alpha*h, 1+self.alpha*w
        F = 1-self.alpha*g*g/(6*rJ*rJ)
        if (np.any(rJ <= self.bounds["r_min"]) or np.any(abs(rJx) >= self.bounds["v_max"])
                or np.any(F <= self.bounds["F_min"])):
            raise RuntimeError("source integration left the admitted continuation tube")
        r = rJ*np.sqrt(F)
        rx = (rJ*rJx-self.alpha*g*gx/6)/r
        Q = (g*gxx-2*gx*gx)/(6*rJ)+g*gx*rJx/(3*rJ*rJ)
        u, D = g/r, 1/F
        uU, ux = (bc-u*a)/r, (gx-u*rx)/r
        K = (a*rx+F/4)/(r*r)-self.alpha*uU*ux/(6*D)
        if not all(np.all(np.isfinite(value)) for value in (rJ, rJx, F, r, rx, Q, K, lam)):
            raise RuntimeError("nonfinite exact source state")
        return rJ, rJx, r, rx, F, Q, uU, ux, K

    def _rhs(self, x, state):
        self.nfev += 1
        g, gx, gxx = self._profile(x)
        _, _, r, rx, F, Q, uU, ux, K = self._quantities(x, state, g, gx, gxx)
        _, w, a, _, _ = state
        return [w, Q, -a*rx/r-F/(4*r), -K*g,
                2*a*rx/(r*r)+F/(2*r*r)-self.alpha*F*uU*ux]

    def sample(self, x):
        """Source first jets, geometric mass and analytical constraint residuals."""
        x = _query(x, self.domain["x_a"], self.domain["x_b"], "source query")
        indices = np.minimum(np.searchsorted(self.knots[1:], x, side="left"), 2)
        flat = x.reshape(-1)
        state = np.empty((5, flat.size))
        for index, segment in enumerate(self._segments):
            chosen = indices.reshape(-1) == index
            if np.any(chosen):
                state[:, chosen] = segment(flat[chosen])
        state = state.reshape((5,)+x.shape)
        g, gx, gxx = self._profile(x)
        rJ, rJx, r, rx, F, Q, _, ux, _ = self._quantities(x, state, g, gx, gxx)
        h, w, a, bc, lam = state
        rJxx = self.alpha*Q
        rxx = (rJx*rJx+rJ*rJxx-self.alpha*(gx*gx+g*gxx)/6-rx*rx)/r
        pJ = g/rJ
        pJx = (gx-pJ*rJx)/rJ
        pJxx = (gxx-2*pJx*rJx-pJ*rJxx)/rJ
        Lx = -self.alpha*pJ*pJx/(3*F)
        result = {"x": x, "rJ": rJ, "rJ_x": rJx, "h": h, "h_x": w,
                  "r": r, "C": F, "log_C": np.log1p(-self.alpha*g*g/(6*rJ*rJ)),
                  "c": g, "r_Ubar": a, "r_x": rx,
                  "c_Ubar": bc, "c_x": gx, "log_C_Ubar": lam, "log_C_x": Lx,
                  "m": r*(1+4*a*rx/F)/2, "r_xx": rxx,
                  "source_constraint": rxx-Lx*rx+self.alpha*r*F*ux*ux/2,
                  "source_constraint_jordan": np.sqrt(F)*rJxx
                      +self.alpha*rJ*(2*pJx*pJx-pJ*pJxx)/(6*np.sqrt(F))}
        _check_budget(self.budget)
        if not all(np.all(np.isfinite(value)) for value in result.values()):
            raise RuntimeError("nonfinite exact source readout")
        return result

    def quiet_bottom(self, Ubar):
        """Exact quiet Schwarzschild data on x=x_a, original U0<=U<=Ua.

        The same affine V(x) change is used off the source; its longitudinal
        lapse derivative therefore retains V_xx/V_x at the corner.
        """
        Ubar = _query(Ubar, self.domain["Ubar_a"], self.domain["Ubar_b"], "bottom query")
        _check_budget(self.budget, int(Ubar.size))
        U = self.geometry["U0"]+Ubar/self.CJ0
        V = self.geometry["Va"]
        r = radius_from_null(U, V)
        xa = self.domain["x_a"]
        Vx = -xa*math.exp(xa/.8)/(.8**2*self.geometry["U0"])
        exp = np.exp(-r/.8)
        rU = -.64*V*exp/(r*self.CJ0)
        rx = -.64*U*exp*Vx/r
        C = 2.048*exp*Vx/(r*self.CJ0)
        zero = np.zeros_like(r)
        result = {"Ubar": Ubar, "U": U, "x": np.full_like(r, xa), "rJ": r, "rJ_x": rx,
                  "h": zero, "h_x": zero, "r": r, "C": C,
                  "log_C": np.log(xa/r)+(xa-r)/.8, "c": zero,
                  "r_Ubar": rU, "r_x": rx, "c_Ubar": zero, "c_x": zero,
                  "log_C_Ubar": -(1/.8+1/r)*rU,
                  "log_C_x": -(1/.8+1/r)*rx+1/xa+1/.8,
                  "m": r*(1+4*rU*rx/C)/2}
        _check_budget(self.budget)
        return result
