"""Floating replay allowances; never scientific error or source tolerances."""
import numpy as np
from .nsc_atlas_field_preparation import radius_product


def source_jet_roundoff_tolerances(x, geometry, *, second=False):
    """Absolute replay allowances for (g, g_x), optionally also g_xx.

    This is a floating evaluation allowance, not a field error estimate. The
    source label subtracts two O(1) null coordinates before dividing by their
    smaller span. Propagate eight binary64 rounding units in that expression
    through b(z) and b'(z)*z_x; a relative test on g_x alone fails near b'=0.
    Source values, preparation, and all numerical accuracy targets are unchanged.
    """
    x = np.asarray(x, dtype=float)
    U0, Va, Vb = (geometry[key] for key in ("U0", "Va", "Vb"))
    V = radius_product(x)/U0
    z = 4*(V-(Va+Vb)/2)/(Vb-Va)
    active = (abs(z) < 1) & (x > geometry["support_x"][0]) & (x < geometry["support_x"][1])
    safe = np.where(active, z, 0.)
    bz = -8*safe*(1-safe*safe)**3
    bzz = 8*(1-safe*safe)**2*(7*safe*safe-1)
    zx = -4*x*np.exp(x/.8)/(.8**2*U0*(Vb-Va))
    eps = np.finfo(float).eps
    dz = 8*eps*(4*(abs(V)+abs((Va+Vb)/2))/abs(Vb-Va)+abs(z))
    value = abs(bz)*dz+16*eps
    first = abs(bzz*zx)*dz+16*eps*np.maximum(1., abs(bz*zx))
    result = (np.where(active, value, 0.), np.where(active, first, 0.))
    if second:
        bzzz = 48*safe*(1-safe*safe)*(3-7*safe*safe)
        zxx = -4*np.exp(x/.8)*(1+x/.8)/(.8**2*U0*(Vb-Va))
        curvature = (abs(bzzz*zx*zx+bzz*zxx)*dz
                     +32*eps*np.maximum(1., abs(bzz*zx*zx)+abs(bz*zxx)))
        result += (np.where(active, curvature, 0.),)
    # The branch is explicitly quiet in both the saved and current expression.
    return result


class _Roundoff:
    """Value/radius arithmetic for validation expressions, not a solver."""
    __array_priority__ = 10000
    def __init__(self, value, error=0.):
        self.value, self.error = np.asarray(value, dtype=float), np.asarray(error, dtype=float)

    @staticmethod
    def cast(value):
        return value if isinstance(value, _Roundoff) else _Roundoff(value)

    @classmethod
    def rounded(cls, value, error):
        return cls(value, error+8*np.finfo(float).eps*(abs(value)+error))

    def __add__(self, other):
        other = self.cast(other)
        return self.rounded(self.value+other.value, self.error+other.error)

    __radd__ = __add__

    def __neg__(self):
        return _Roundoff(-self.value, self.error)

    def __sub__(self, other):
        return self+-self.cast(other)

    def __rsub__(self, other):
        return self.cast(other)+-self

    def __mul__(self, other):
        other = self.cast(other)
        return self.rounded(self.value*other.value,
            abs(self.value)*other.error+abs(other.value)*self.error+self.error*other.error)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.cast(other)
        denominator = abs(other.value)-other.error
        if np.any(denominator <= 0):
            raise ValueError("source replay denominator is not separated from zero")
        value = self.value/other.value
        return self.rounded(value, (self.error+abs(value)*other.error)/denominator)

    def __rtruediv__(self, other):
        return self.cast(other)/self

    def sqrt(self):
        if np.any(self.value-self.error <= 0):
            raise ValueError("source replay square root is not positive")
        value = np.sqrt(self.value)
        return self.rounded(value, self.error/(value+np.sqrt(self.value-self.error)))

    def log1p(self):
        denominator = 1+self.value-self.error
        if np.any(denominator <= 0):
            raise ValueError("source replay logarithm is not defined")
        return self.rounded(np.log1p(self.value), self.error/denominator)


def _hermite_roundoff_tolerances(x, knots, *, slope_errors=None):
    """Cubic coefficient/evaluation conditioning from the actual saved knots."""
    nodes, values, slopes = (np.asarray(knots[key], dtype=float) for key in ("x", "h", "h_x"))
    if (nodes.ndim != 1 or len(nodes) < 2 or values.shape != nodes.shape or slopes.shape != nodes.shape
            or not all(np.all(np.isfinite(v)) for v in (nodes, values, slopes))
            or np.any(np.diff(nodes) <= 0) or np.any(x < nodes[0]) or np.any(x > nodes[-1])):
        raise ValueError("source replay requires the saved increasing Hermite knots")
    i = np.clip(np.searchsorted(nodes, x, side="right")-1, 0, len(nodes)-2)
    d = nodes[i+1]-nodes[i]
    t = (x-nodes[i])/d
    y0, y1, m0, m1 = values[i], values[i+1], slopes[i], slopes[i+1]
    a0, a1 = y0, d*m0
    a2 = 3*(y1-y0)-d*(2*m0+m1)
    a3 = 2*(y0-y1)+d*(m0+m1)
    eps = np.finfo(float).eps
    slope_errors = np.zeros_like(slopes) if slope_errors is None else np.asarray(slope_errors, dtype=float)
    if slope_errors.shape != slopes.shape or not np.all(np.isfinite(slope_errors)) or np.any(slope_errors < 0):
        raise ValueError("finite nonnegative source replay slope radii required")
    em0, em1 = slope_errors[i], slope_errors[i+1]
    e1 = d*em0+8*eps*abs(a1)
    e2 = d*(2*em0+em1)+16*eps*(3*(abs(y1)+abs(y0))+d*(2*abs(m0)+abs(m1)))
    e3 = d*(em0+em1)+16*eps*(2*(abs(y0)+abs(y1))+d*(abs(m0)+abs(m1)))
    et = 8*eps*((abs(x)+abs(nodes[i]))/d+abs(t))
    e0 = e1*abs(t)+e2*t*t+e3*abs(t*t*t)
    e0 += 8*eps*(abs(a0)+abs(a1*t)+abs(a2*t*t)+abs(a3*t*t*t))
    e0 += (abs(a1)+2*abs(a2*t)+3*abs(a3*t*t))*et
    ex = (e1+2*abs(t)*e2+3*t*t*e3
          +8*eps*(abs(a1)+2*abs(a2*t)+3*abs(a3*t*t))
          +(2*abs(a2)+6*abs(a3*t))*et)/d
    exx = (2*e2+6*abs(t)*e3+8*eps*(2*abs(a2)+6*abs(a3*t))+6*abs(a3)*et)/(d*d)
    return dict(h=e0, h_x=ex, h_xx=exx)


def _mixed_roundoff(row, alpha, errors=None):
    """Exterior mixed expressions with retained first jets and roundoff radii."""
    errors = {} if errors is None else errors
    r, C, c, rU, rx, cU, cx = (_Roundoff(row[key], errors.get(key, 0.))
        for key in ("r", "C", "c", "r_Ubar", "r_x", "c_Ubar", "c_x"))
    u = c/r
    uU, ux = (cU-u*rU)/r, (cx-u*rx)/r
    D = 1+alpha*u*u/6
    sigma = alpha*uU*ux/D
    rr, product = r*r, rU*rx
    common = product+C/4
    return dict(r_Ubar=-common/r, log_C_Ubar=2*common/rr-sigma,
        c_Ubar=-(common/rr-sigma/6)*c, m=r*(1+4*product/C)/2)


def source_snapshot_roundoff_tolerances(row, geometry, alpha, *, source_knots=None):
    """Propagate analytic-jet replay uncertainty through snapshot diagnostics.

    Residuals such as h_xx-Q can be small differences of large terms. Their
    allowance follows those terms, not the residual's own relative magnitude.
    Eight rounding units per elementary operation complement the independently
    conditioned compact-jet radii. This is validation conditioning only: it
    neither estimates interpolation/ODE error nor replaces any retained value.
    Optional source_knots supplies the original x,h,h_x Hermite inputs: their
    coefficient cancellation and 1/dx,1/dx^2 derivative amplification are
    propagated too. Partial manufactured rows receive available allowances.
    """
    x = np.asarray(row["x"], dtype=float)
    result = dict(zip(("c", "c_x", "c_xx"),
        source_jet_roundoff_tolerances(x, geometry, second=True)))
    required = {"h", "h_x", "h_xx", "c", "c_x", "c_xx"}
    if not required <= set(row):
        return result
    eps = np.finfo(float).eps
    hermite = {} if source_knots is None else _hermite_roundoff_tolerances(x, source_knots)
    inputs = {key: _Roundoff(row[key], hermite.get(key, 8*eps*abs(np.asarray(row[key], dtype=float))))
              for key in ("h", "h_x", "h_xx")}
    result.update({key: value.error for key, value in inputs.items()})
    transverse = {}
    if source_knots is not None and {"r", "C", "c", "r_Ubar", "r_x", "c_Ubar", "c_x"} <= set(source_knots):
        slopes = _mixed_roundoff(source_knots, alpha)
        for key in ("r_Ubar", "log_C_Ubar", "c_Ubar"):
            if key not in source_knots:
                continue
            bounds = _hermite_roundoff_tolerances(x,
                dict(x=source_knots["x"], h=source_knots[key], h_x=slopes[key].value),
                slope_errors=slopes[key].error)
            transverse[key] = bounds["h"]
            transverse[key+"_x"] = bounds["h_x"]
        result.update(transverse)
    h, hx, hxx = (inputs[key] for key in ("h", "h_x", "h_xx"))
    g, gx, gxx = (_Roundoff(row[key], result[key]) for key in ("c", "c_x", "c_xx"))
    rJ, rJx, rJxx = x+alpha*h, 1+alpha*hx, alpha*hxx
    correction = -alpha*g*g/(6*rJ*rJ)
    F = 1+correction
    sqrtF = F.sqrt()
    r = rJ*sqrtF
    rx = (rJ*rJx-alpha*g*gx/6)/r
    rxx = (rJx*rJx+rJ*rJxx-alpha*(gx*gx+g*gxx)/6-rx*rx)/r
    p = g/rJ
    px = (gx-p*rJx)/rJ
    pxx = (gxx-2*px*rJx-p*rJxx)/rJ
    Lx = -alpha*p*px/(3*F)
    ux = (gx-g*rx/r)/r
    Q = (g*gxx-2*gx*gx)/(6*rJ)+g*gx*rJx/(3*rJ*rJ)
    defect = hxx-Q
    values = dict(rJ=rJ, rJ_x=rJx, rJ_xx=rJxx, r=r, r_x=rx, r_xx=rxx,
        C=F, log_C=correction.log1p(), log_C_x=Lx, source_ode_h_xx=Q,
        source_h_xx_defect=defect,
        source_constraint=rxx-Lx*rx+alpha*r*F*ux*ux/2,
        source_constraint_jordan=sqrtF*rJxx+alpha*rJ*(2*px*px-p*pxx)/(6*sqrtF),
        source_constraint_from_curvature_defect=alpha*defect/sqrtF)
    # Transverse polynomials use the retained first jets. Their reported mixed
    # defects subtract the smooth exterior equations; preserve that term scale.
    if {"r_Ubar", "c_Ubar"} <= set(row):
        rU, cU = (_Roundoff(row[key], transverse.get(key, 8*eps*abs(np.asarray(row[key], dtype=float))))
                   for key in ("r_Ubar", "c_Ubar"))
        u = g/r
        uU = (cU-u*rU)/r
        D = 1+alpha*u*u/6
        sigma = alpha*uU*ux/D
        mixed = dict(r_Ubar=-rU*rx/r-F/(4*r),
            log_C_Ubar=2*rU*rx/(r*r)+F/(2*r*r)-sigma,
            c_Ubar=-((rU*rx+F/4)/(r*r)-sigma/6)*g)
        for key, rhs in mixed.items():
            if key+"_x" in row:
                derivative = np.asarray(row[key+"_x"], dtype=float)
                values["mixed_"+key+"_defect"] = _Roundoff(derivative,
                    transverse.get(key+"_x", 8*eps*abs(derivative)))-rhs
        values["m"] = r*(1+4*rU*rx/F)/2
    result.update({key: np.broadcast_to(value.error, x.shape) for key, value in values.items()})
    return result
