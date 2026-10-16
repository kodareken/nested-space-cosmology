"""Floating replay allowances; never scientific error or source tolerances."""
import math
import numpy as np
from .nsc_atlas_field_preparation import radius_product


class EndpointReplaySpline:
    """Evaluate a restored spline at an admitted rounded domain endpoint.

    Saved and rebuilt reference bounds already admit eight binary64 units.
    An endpoint lying just outside the rebuilt, nonextrapolating spline is
    evaluated at that spline's endpoint. Interior arguments are unchanged;
    larger excursions and nonfinite arguments still fail. This adapter is
    for validation replay, never a source or evolution extrapolator.
    """
    def __init__(self, spline):
        self._spline = spline

    def __getattr__(self, name):
        return getattr(self._spline, name)

    def __call__(self, values, *args, **kwargs):
        points = np.asarray(values, dtype=float)
        low, high = float(self._spline.x[0]), float(self._spline.x[-1])
        lower_allowance = 8*np.finfo(float).eps*np.maximum(1., np.maximum(abs(points), abs(low)))
        upper_allowance = 8*np.finfo(float).eps*np.maximum(1., np.maximum(abs(points), abs(high)))
        if (not np.all(np.isfinite(points)) or np.any(points<low-lower_allowance)
                or np.any(points>high+upper_allowance)):
            raise ValueError("reference replay query leaves its rounded spline endpoints")
        return self._spline(np.clip(points, low, high), *args, **kwargs)


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

    def log(self):
        if np.any(self.value-self.error <= 0):
            raise ValueError("replay logarithm is not positive")
        return self.rounded(np.log(self.value), self.error/(self.value-self.error))

    def exp(self):
        value = np.exp(self.value)
        return self.rounded(value, value*np.expm1(self.error))

    def __getitem__(self, index):
        return _Roundoff(self.value[index], np.broadcast_to(self.error, self.value.shape)[index])


def departure_history_roundoff_tolerances(row, source, alpha):
    """Replay conditioning of departure diagnostics, using their actual operands.

    Interior departures, source jets and gauge samples are immutable saved
    inputs. Only the analytic vacuum and terminal subtraction are rebuilt.
    Eight binary64 rounding units per operation describe alternate arithmetic
    paths (including the elementary functions). The Lambert-W input error
    uses |dw/dq|=|w|/(|q|(1+w)); a factor two bounds its variation on the
    explicitly small admitted interval. This is a floating replay allowance,
    not a certified SciPy error bound or a physical/truncation error estimate.

    Propagation through the *actual* centred stencils retains 1/de^2, 1/dx^2,
    mixed derivatives and moving-grid coefficients. Two expression radii
    compare independent evaluations. Summary extrema use their candidate
    intervals; RMS uses the normalized Euclidean norm. Counts stay exact.
    No evolution, fit, vacuum query or field modification occurs here.
    """
    Q, eps, a = _Roundoff, np.finfo(float).eps, .8
    x, eta = np.asarray(row["x"]), np.asarray(row["eta"])
    g = row["gauge"]
    U = Q(eta[None, :])*Q(np.asarray(g["W"])[:, None])
    t = Q(x[:, None])
    q = (1-U/(4*a))*(1-t/a)*(t/a).exp()/np.e
    # Recover the already evaluated vacuum radius, including the rounding
    # radius of this recovery. The source endpoint is an exact affine anchor.
    off = Q(np.asarray(source["r"]))-Q(x)
    recovered = Q(row["states"]["r"])-off[:, None]-Q(row["departure"]["r"])
    rv = recovered.value
    w = rv/a-1
    relative_input = q.error/(abs(q.value)-q.error)
    branch_radius = 2*abs(w)*relative_input/(1+w)+8*eps*abs(w)
    if (np.any(rv <= 0) or np.any(rv >= a) or not np.all(np.isfinite(branch_radius))
            or np.any(branch_radius > (1+w)/16)):
        raise ValueError("vacuum replay leaves the separated Lambert-W branch")
    vr = Q(rv, recovered.error+a*branch_radius+8*eps*abs(rv))
    vr.value[:, 0], vr.error[:, 0] = x, 0.
    vL = (t/vr).log()+(t-vr)/a
    C = vL.exp()
    p, R = C*(a-t)/(4*t), (1-U/(4*a))*C
    B = 1/vr+1/a
    LU, Lx = -B*p, 1/t+1/a-B*R
    vacuum = {"r": vr, "L": vL, "c": Q(np.zeros_like(rv))}
    offsets = {"r": off, "L": Q(source["log_C"]), "c": Q(source["c"])}
    first = {"r": Q(source["r_x"])-1, "L": Q(source["log_C_x"]), "c": Q(source["c_x"])}
    z, states = {}, {}
    for key in ("r", "L", "c"):
        value = np.asarray(row["departure"][key])
        error = np.zeros_like(value)
        terminal = Q(np.asarray(row["states"][key])[:, -1])-vacuum[key][:, -1]-offsets[key]
        error[:, -1] = terminal.error
        z[key] = Q(value, error)
        rebuilt = vacuum[key]+offsets[key][:, None]+z[key]
        # Both physical endpoints are supplied, not rebuilt.
        rebuilt.error[:, [0, -1]] = 0.
        states[key] = Q(row["states"][key], rebuilt.error)
    take = np.s_[1:-1, 1:-1]
    dx, de = Q(x[1])-Q(x[0]), Q(eta[1])-Q(eta[0])
    W = Q(np.asarray(g["W"])[1:-1, None])
    b = Q(np.asarray(g["W_x"])[1:-1, None])/W
    transport = Q(eta[None, 1:-1])*b
    second_map = Q(eta[None, 1:-1])*(2*b*b-Q(np.asarray(g["W_xx"])[1:-1, None])/W)
    derivatives = {}
    for key, vu, vx in (("r", p, R), ("L", LU, Lx), ("c", Q(np.zeros_like(rv)), Q(np.zeros_like(rv)))):
        f = z[key]
        d_eta = (f[1:-1, 2:]-f[1:-1, :-2])/(2*de)
        d_x = (f[2:, 1:-1]-f[:-2, 1:-1])/(2*dx)
        derivatives[key+"_Ubar"] = vu[take]+d_eta/W
        derivatives[key+"_x"] = vx[take]+first[key][1:-1, None]+d_x-transport*d_eta
        if key == "r":
            d_ee = (f[1:-1, 2:]-2*f[take]+f[1:-1, :-2])/(de*de)
            d_xx = (f[2:, 1:-1]-2*f[take]+f[:-2, 1:-1])/(dx*dx)
            d_xe = (f[2:, 2:]-f[2:, :-2]-f[:-2, 2:]+f[:-2, :-2])/(4*dx*de)
            derivatives["r_Ubar_Ubar"] = (LU*p)[take]+d_ee/(W*W)
            derivatives["r_xx"] = ((Lx*R)[take]+Q(source["r_xx"])[1:-1, None]+d_xx
                -2*transport*d_xe+transport*transport*d_ee+second_map*d_eta)
    r, c = (states[k][take] for k in ("r", "c"))
    u = c/r
    D = 1+alpha*u*u/6
    values = {}
    for coord, second in (("Ubar", "r_Ubar_Ubar"), ("x", "r_xx")):
        ua = (derivatives["c_"+coord]-u*derivatives["r_"+coord])/r
        values["E_"+coord] = (derivatives[second]-derivatives["L_"+coord]*derivatives["r_"+coord]
                              +alpha*r*ua*ua/(2*D))
    result = {}
    for key, value in values.items():
        bound = 2*value.error
        if not np.all(np.isfinite(bound)) or np.any(bound < 0):
            raise ValueError("nonfinite departure replay conditioning")
        result[(key,)] = bound
        for name, mask in row["masks"].items():
            selected = bound[mask]
            if not selected.size:
                continue
            samples = np.asarray(row[key])[mask]
            for statistic, center, lower, upper in (
                    ("min", np.min(samples), np.min(samples-selected), np.min(samples+selected)),
                    ("max", np.max(samples), np.max(samples-selected), np.max(samples+selected)),
                    ("max_abs", np.max(abs(samples)),
                     np.max(np.maximum(0., abs(samples)-selected)), np.max(abs(samples)+selected))):
                result[("summaries", key, name, statistic)] = float(max(center-lower, upper-center))
            rms = float(row["summaries"][key][name]["rms"])
            gamma = (selected.size+8)*eps/(1-(selected.size+8)*eps)
            result[("summaries", key, name, "rms")] = float(np.sqrt(np.mean(selected*selected)))+2*gamma*abs(rms)
    return result


def collect_departure_replay(allowances, case_index, source, alpha, **histories):
    """Attach transient path allowances only when a checker requests them."""
    if allowances is None:
        return
    for name, row in histories.items():
        stride = (len(source["r"])-1)//(len(row["x"])-1)
        subset = {key: np.asarray(value)[::stride] for key, value in source.items()}
        for path, bound in departure_history_roundoff_tolerances(row, subset, alpha).items():
            allowances[("histories", case_index, name, *path)] = bound


def background_profile_roundoff_tolerances(inputs):
    """Condition the retained C2 profile's arithmetic, not its fit error.

    Labels and first jets use the saved wall values/rates. Computational
    curvature uses the stated nonuniform three-point gradient. For each
    quintic, Bernstein coefficients are converted to normalized power
    coefficients; value/first/second evaluation then exposes the actual
    1/d and 1/d^2 cancellation amplification. Eight rounding units per
    operation cover alternate binary64 evaluation paths. Two such radii
    compare independently rebuilt summaries. No data or interpolant changes.
    """
    table = inputs["table"]
    walls, rates = table["walls"], table["wall_rates"]
    r = np.asarray(walls["r"], dtype=float)
    T = np.asarray(walls["T"], dtype=float)
    c = np.asarray(walls["c"], dtype=float)
    if (r.shape != (129,) or T.shape != r.shape or c.shape != r.shape
            or np.any(r <= 0) or np.any(r >= 1)):
        raise ValueError("profile replay requires the retained 129 wall rows")
    eps = np.finfo(float).eps
    distance = np.arctanh(r)
    labels = _Roundoff(T)+_Roundoff(distance, 8*eps*abs(distance))
    arrival = _Roundoff(rates["T"])+_Roundoff(rates["r"])/(1-_Roundoff(r)*r)
    first = _Roundoff(rates["c"])/arrival

    def take(value, index):
        return _Roundoff(value.value[index], np.broadcast_to(value.error, value.value.shape)[index])

    knots, slopes = take(labels, slice(None, None, 2)), take(first, slice(None, None, 2))
    if np.any(np.diff(knots.value) <= 0):
        raise ValueError("profile replay labels must increase")
    curvature = []
    for i in range(len(knots.value)):
        j = min(max(i-1, 0), len(knots.value)-3)
        h0, h1 = take(knots,j+1)-take(knots,j), take(knots,j+2)-take(knots,j+1)
        if i == 0:
            weights = (-(2*h0+h1)/(h0*(h0+h1)), (h0+h1)/(h0*h1), -h0/(h1*(h0+h1)))
        elif i == len(knots.value)-1:
            weights = (h1/(h0*(h0+h1)), -(h0+h1)/(h0*h1), (h0+2*h1)/(h1*(h0+h1)))
        else:
            weights = (-h1/(h0*(h0+h1)), (h1-h0)/(h0*h1), h0/(h1*(h0+h1)))
        curvature.append(sum(w*take(slopes,j+k) for k,w in enumerate(weights)))
    second = _Roundoff([v.value for v in curvature], [v.error for v in curvature])
    values = _Roundoff(c[::2])

    def evaluation_errors(queries):
        index = np.clip(np.searchsorted(knots.value, queries.value, side="right")-1, 0, 63)
        d = take(knots,index+1)-take(knots,index)
        t = (queries-take(knots,index))/d
        y0,y1 = take(values,index),take(values,index+1)
        f0,f1 = take(slopes,index),take(slopes,index+1)
        s0,s1 = take(second,index),take(second,index+1)
        # Degree-5 Bernstein data fixed by endpoint values and first two jets.
        b = [y0, y0+d*f0/5, y0+2*d*f0/5+d*d*s0/20,
             y1-2*d*f1/5+d*d*s1/20, y1-d*f1/5, y1]
        a = [math.comb(5,k)*sum((-1)**(k-j)*math.comb(k,j)*b[j]
             for j in range(k+1)) for k in range(6)]
        result = []
        for order in range(3):
            value = _Roundoff(np.zeros_like(queries.value))
            for k in range(5,order-1,-1):
                value = value*t+a[k]*(math.factorial(k)//math.factorial(k-order))
            for _ in range(order):
                value = value/d
            result.append(value.error)
        return result

    knot_errors = evaluation_errors(knots)
    withheld_errors = evaluation_errors(take(labels,slice(1,None,2)))
    value_error = 2*withheld_errors[0]
    first_error = 2*(withheld_errors[1]+first.error[1::2])
    material = slopes*take(arrival,slice(None,None,2))-_Roundoff(np.asarray(rates["c"])[::2])
    spacing = take(knots,slice(1,None))-take(knots,slice(None,-1))
    result = {
        "labels": 2*knots.error, "label_range": 2*knots.error[[0,-1]],
        "arrival_rate": 2*arrival.error[::2], "first": 2*slopes.error,
        "second": 2*second.error,
        "diagnostics.no_return_margin": float(2*knots.error[0]),
        "diagnostics.minimum_arrival_rate": float(2*np.max(arrival.error[::2])),
        "diagnostics.minimum_label_spacing": float(2*np.max(spacing.error)),
        "diagnostics.max_abs_material_chain_defect": float(2*np.max(material.error)),
        "diagnostics.max_abs_knot_value_defect": float(2*np.max(knot_errors[0])),
        "diagnostics.max_abs_knot_first_defect": float(2*np.max(knot_errors[1]+slopes.error)),
        "diagnostics.max_abs_knot_second_defect": float(2*np.max(knot_errors[2]+second.error)),
        "withheld.labels": 2*labels.error[1::2],
        "withheld.value_difference": value_error, "withheld.first_difference": first_error,
    }
    for name,error in (("value",value_error),("first",first_error)):
        result["withheld."+name+"_statistics.maximum_absolute"] = float(np.max(error))
        result["withheld."+name+"_statistics.rms"] = float(np.sqrt(np.mean(error*error)))
    if any(not np.all(np.isfinite(v)) or np.any(np.asarray(v)<0) for v in result.values()):
        raise ValueError("nonfinite profile replay conditioning")
    return result


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
