"""The prepared leading scalar in the finite regular de Sitter cap.

There is no second scalar evolution. The existing wall field determines
f(q), with q the reference central arrival time. Both f(T+r*) and the
centre-reflected f(T-r*) are retained. All coordinates are dimensionless.
"""
from collections import OrderedDict
import math

import numpy as np

from .nsc_atlas_field_preparation import M_OVER_ELL, exterior_coefficient, radius_from_null, source_profile
from .nsc_atlas_prepared_response import MAX_FIELD_NODES, _minimum_derivative

CENTRE_RADII = (1e-4, 1e-3)
MAX_PROFILE_CACHE = 32768


def _choice(value, choices, name):
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer))
            or value not in choices):
        raise ValueError(name+" must be one of "+str(choices))
    return int(value)


def _array(value, name):
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf" or raw.dtype.kind == "b":
        raise ValueError(name+" must contain finite real numbers")
    result = np.asarray(value, dtype=float)
    if result.ndim > 1 or not 1 <= result.size <= MAX_FIELD_NODES or not np.all(np.isfinite(result)):
        raise ValueError(name+" needs 1 to 513 finite scalar or one-dimensional values")
    return result


def source_second_derivative(V, Va, Vb):
    """Analytic g'' for the unchanged normalized compact C3 preparation."""
    z = 4*(np.asarray(V)-(Va+Vb)/2)/(Vb-Va)
    return np.where(abs(z) < 1, 8*(1-z*z)**2*(7*z*z-1)*(4/(Vb-Va))**2, 0.)


class CapLeadingField:
    """First static-cap jets from one prepared field and its shared budget.

    ``profile(q, second=False)`` returns f and f_prime, and optionally
    f_second. The latter uses the differentiated continuum integral with
    the numerical scalar field; it is not the exact derivative of its
    semidiscrete interpolant. Geometry and field refinements remain needed.
    ``jet(T,r)`` returns p, p_T, p_r and p_center=2*f_prime(T).
    Small-radius jets integrate f' and f'' instead of subtracting nearly
    equal field values for r<=centre_radius (1e-4 by default; 1e-3 is a
    separate numerical control). The inherited arrival spline is differentiated
    analytically; no saved waveform samples are differentiated.
    """

    def __init__(self, prepared_field, *, derivative_order=2, small_radius_order=8,
                 centre_radius=1e-4):
        self.derivative_order = _choice(derivative_order, (2, 4), "cap derivative order")
        self.small_radius_order = _choice(small_radius_order, (8, 16), "small-radius order")
        if (isinstance(centre_radius, (bool, np.bool_))
                or not isinstance(centre_radius, (int, float, np.integer, np.floating))
                or centre_radius not in CENTRE_RADII):
            raise ValueError("centre radius must be 1e-4 or 1e-3")
        self.centre_radius = float(centre_radius)
        self.prepared_field = prepared_field
        self.geometry = prepared_field.geometry
        self.budget = prepared_field.budget
        dense = prepared_field.dense_data
        certificate = dense.get("prepared", {})
        observer = certificate.get("observer", {})
        preparation = certificate.get("preparation", {})
        margins = [observer.get("no_return_margin_over_ell", float("nan")),
                   observer.get("analytic_no_return_log_margin_lower_bound", float("nan"))]
        if (not all(math.isfinite(float(value)) and float(value) > 0 for value in margins)
                or preparation.get("corner_compatible") is not True
                or preparation.get("cap_data") != "phi=0 and partial_t phi=0 on the flat cap slice t=t_a"):
            raise ValueError("cap field requires the existing quiet-cap and no-return certificate")
        self.quiet_certificate = {"no_return_margin_over_ell": float(margins[0]),
            "analytic_no_return_log_margin_lower_bound": float(margins[1]), "quiet_cap": True}
        self.domain = dict(prepared_field.domain)
        self.Va, self.Vb = self.domain["Va"], self.domain["Vb"]
        self.dV = self.Vb-self.Va
        arrival = self.geometry["arrival"]
        if (not hasattr(arrival, "x") or not hasattr(arrival, "c")
                or arrival.x[0] > self.Va or arrival.x[-1] < self.Vb
                or _minimum_derivative(arrival) <= 0):
            raise ValueError("cap arrival spline must be monotone on the full preparation")
        qa, qb = map(float, arrival([self.Va, self.Vb]))
        xa, xb = map(float, radius_from_null([self.domain["Ua"], self.domain["Ub"]], [self.Va, self.Vb]))
        if not 0 < xa < xb < 1 or not qa < qb:
            raise ValueError("cap needs the admitted regular static-patch window")
        Ta, Tb = qa-np.arctanh(xa), qb-np.arctanh(xb)
        if not Ta < Tb < qb:
            raise ValueError("cap wall static clock must increase before central arrival")
        if not Tb-np.arctanh(xb) < qa:
            raise ValueError("cap wall labels violate the inherited no-return condition")
        self.domain.update(qa=qa, qb=qb, xa=xa, xb=xb, Ta=float(Ta), Tb=float(Tb),
                           qmin=float(Ta-np.arctanh(xa)))
        self.edges = np.asarray([qa, *arrival(self.Va+self.dV*np.array([.25, .75])), qb])
        self._cache = OrderedDict()
        self.evaluation_counts = dict(profile_requested=0, wall_points=0, second_points=0, cache_hits=0)
        self.max_inverse_bracket_width = 0.
        self._radial_gauss = np.polynomial.legendre.leggauss(self.derivative_order)
        self._centre_gauss = np.polynomial.legendre.leggauss(self.small_radius_order)
        self._grid, self._cuts, self._segments = dense["grid"], dense["cuts"], dense["segments"]
        self.budget.check()

    def _static_wall(self, V):
        return self.geometry["arrival"](V)-np.arctanh(radius_from_null(self.geometry["wall"](V), V))

    def _inverse(self, targets, function, first, last):
        """Bounded inversion of floating spline values, not an interval proof."""
        targets = np.asarray(targets, dtype=float)
        low, high = np.full_like(targets, self.Va), np.full_like(targets, self.Vb)
        if np.any(targets < first) or np.any(targets > last):
            raise ValueError("cap inverse target exits the known preparation")
        high[targets == first] = self.Va
        low[targets == last] = self.Vb
        for _ in range(56):
            midpoint = low+.5*(high-low)
            active = (midpoint > low) & (midpoint < high)
            if not np.any(active):
                break
            before = function(midpoint) < targets
            low = np.where(active & before, midpoint, low)
            high = np.where(active & ~before, midpoint, high)
        width = high-low
        if np.any(width > 4*np.finfo(float).eps*max(abs(self.Va), abs(self.Vb))):
            raise RuntimeError("cap inverse failed its finite bracket tolerance")
        self.max_inverse_bracket_width = max(self.max_inverse_bracket_width, float(np.max(width, initial=0.)))
        self.budget.check()
        return high

    def _second_wall(self, V):
        """Continuum w_VV using first scalar jets and cell-split quadrature."""
        t = min(1., max(0., (V-self.Va)/self.dV))
        index = min(np.searchsorted(self._cuts, t, side="right")-1, len(self._segments)-1)
        correction = np.asarray(self._segments[index](t)).reshape(-1)
        wall = float(self.geometry["wall"](V))
        clipped = np.minimum(self._grid, wall)
        g, gp = source_profile(V, self.Va, self.Vb)
        chi = np.where(self._grid < wall, correction, correction[-1])+g
        radius = radius_from_null(clipped, V)
        self.budget.check(len(clipped))
        K = exterior_coefficient(radius)
        product = K*chi
        prefix = np.r_[0., np.cumsum(.5*(product[:-1]+product[1:])*np.diff(clipped))]
        widths = np.diff(clipped)
        cells = np.flatnonzero(widths > 0)
        node, weight = self._radial_gauss
        fraction = (node+1)/2
        U = clipped[cells, None]+widths[cells, None]*fraction
        # Every spatial cell is clipped to the retained side before its kernel.
        if np.any(U > wall):
            raise RuntimeError("cap derivative integral crossed the exterior wall")
        x = radius_from_null(U, V)
        self.budget.check(U.size)
        kernel = exterior_coefficient(x)
        x_V = -4*M_OVER_ELL**2*U*np.exp(-x/(2*M_OVER_ELL))/x
        kernel_V = -(1/(2*M_OVER_ELL)+4/x)*kernel*x_V
        field = chi[cells, None]+fraction*(chi[cells+1]-chi[cells])[:, None]
        local_product = product[cells, None]+fraction*(product[cells+1]-product[cells])[:, None]
        field_V = gp-prefix[cells, None]-.5*(product[cells, None]+local_product)*(U-clipped[cells, None])
        integral = np.sum(widths[cells]/2*((kernel_V*field+kernel*field_V) @ weight))
        return float(source_second_derivative(V, self.Va, self.Vb)
                     -self.geometry["wall"](V, 1)*K[-1]*chi[-1]-integral)

    def profile(self, q, *, second=False):
        values = _array(q, "cap profile labels")
        if not isinstance(second, (bool, np.bool_)):
            raise ValueError("second derivative selector must be boolean")
        if np.any(values < self.domain["qmin"]) or np.any(values > self.domain["qb"]):
            raise ValueError("cap profile label is outside the finite supplied interval")
        self.budget.check()
        self.evaluation_counts["profile_requested"] += values.size
        unique, inverse = np.unique(values.ravel(), return_inverse=True)
        rows = np.zeros((len(unique), 3))
        missing = []
        for i, qvalue in enumerate(unique):
            if qvalue <= self.domain["qa"]:
                continue
            cached = self._cache.get(float(qvalue))
            if cached is not None and (not second or cached[2] is not None):
                rows[i] = (cached[0], cached[1], cached[2] or 0.)
                self._cache.move_to_end(float(qvalue))
                self.evaluation_counts["cache_hits"] += 1
            else:
                missing.append(i)
        if missing:
            labels = unique[missing]
            V = self._inverse(labels, self.geometry["arrival"], self.domain["qa"], self.domain["qb"])
            times = np.clip((V-self.Va)/self.dV, 0., 1.)
            sample = self.prepared_field.wall(times)
            self.evaluation_counts["wall_points"] += len(labels)
            qV = self.geometry["arrival"](V, 1)
            if np.any(qV <= 0) or not np.all(np.isfinite(qV)):
                raise ValueError("cap arrival Jacobian lost monotonicity")
            for j, i in enumerate(missing):
                w = float(sample["wall_chi_over_epsilon"][j])
                J = float(sample["wall_dchi_dV_over_epsilon"][j])
                fp = J/float(qV[j])
                fpp = None
                if second:
                    Jp = self._second_wall(float(V[j]))
                    fpp = float((Jp*qV[j]-J*self.geometry["arrival"](V[j], 2))/qV[j]**3)
                    self.evaluation_counts["second_points"] += 1
                rows[i] = (w, fp, fpp or 0.)
                self._cache[float(unique[i])] = (w, fp, fpp)
                self._cache.move_to_end(float(unique[i]))
                if len(self._cache) > MAX_PROFILE_CACHE:
                    self._cache.popitem(last=False)
        if not np.all(np.isfinite(rows)):
            raise ValueError("cap profile returned a nonfinite derivative")
        result = {"f": rows[inverse, 0].reshape(values.shape),
                  "f_prime": rows[inverse, 1].reshape(values.shape)}
        if second:
            result["f_second"] = rows[inverse, 2].reshape(values.shape)
        self.budget.check()
        return result

    def _many_profile(self, labels, *, second=False):
        labels = np.asarray(labels)
        unique, inverse = np.unique(labels.ravel(), return_inverse=True)
        keys = ("f", "f_prime", "f_second") if second else ("f", "f_prime")
        data = {key: np.empty(len(unique)) for key in keys}
        for start in range(0, len(unique), MAX_FIELD_NODES):
            stop = min(start+MAX_FIELD_NODES, len(unique))
            sample = self.profile(unique[start:stop], second=second)
            for key in keys:
                data[key][start:stop] = sample[key]
        return {key: value[inverse].reshape(labels.shape) for key, value in data.items()}

    def _cap_time(self, T):
        values = _array(T, "cap static time")
        tolerance = 8*np.finfo(float).eps*max(1., abs(self.domain["Ta"]), abs(self.domain["qb"]))
        if np.any(values < self.domain["Ta"]-tolerance) or np.any(values > self.domain["qb"]+tolerance):
            raise ValueError("cap time exits the known finite domain")
        return np.clip(values, self.domain["Ta"], self.domain["qb"])

    def radial_limit(self, T):
        """Known cap extent; clamp only endpoint roundoff within eight eps."""
        values = self._cap_time(T)
        flat = values.ravel()
        radius = np.tanh(self.domain["qb"]-flat)
        earlier = flat < self.domain["Tb"]
        if np.any(earlier):
            unique, inverse = np.unique(flat[earlier], return_inverse=True)
            V = self._inverse(unique, self._static_wall, self.domain["Ta"], self.domain["Tb"])
            radius[earlier] = radius_from_null(self.geometry["wall"](V), V)[inverse]
        return radius.reshape(values.shape)

    def jet(self, T, r):
        T, r = np.broadcast_arrays(self._cap_time(T), _array(r, "cap radius"))
        if T.size > MAX_FIELD_NODES:
            raise ValueError("cap jet batch exceeds 513 points")
        limit = self.radial_limit(T)
        tolerance = 8*np.finfo(float).eps
        if np.any(r < 0) or np.any(r > limit+tolerance):
            raise ValueError("cap radius exits the retained finite domain")
        # Clamp only endpoint roundoff, never extend the profile beyond qb.
        shape = T.shape
        T, r = T.ravel(), np.minimum(r, limit).ravel()
        d = np.arctanh(r)
        labels = np.clip(np.stack((T-d, T+d)), self.domain["qmin"], self.domain["qb"])
        endpoints = self._many_profile(labels)
        centre = self._many_profile(T)
        p, p_T, p_r = (np.empty_like(T) for _ in range(3))
        regular = r > self.centre_radius
        if np.any(regular):
            rr = r[regular]
            p[regular] = (endpoints["f"][1, regular]-endpoints["f"][0, regular])/rr
            p_T[regular] = (endpoints["f_prime"][1, regular]-endpoints["f_prime"][0, regular])/rr
            p_r[regular] = ((endpoints["f_prime"][1, regular]+endpoints["f_prime"][0, regular])/(1-rr*rr)-p[regular])/rr
        zero = r == 0
        if np.any(zero):
            core = self._many_profile(T[zero], second=True)
            p[zero], p_T[zero], p_r[zero] = 2*core["f_prime"], 2*core["f_second"], 0.
        small = (r > 0) & ~regular
        if np.any(small):
            ss, rr, dd = T[small], r[small], d[small]
            cuts = np.sort(np.column_stack((np.full_like(ss, -1.),
                *[np.clip((edge-ss)/dd, -1., 1.) for edge in self.edges], np.ones_like(ss))), axis=1)
            I1, I2, IT = (np.zeros_like(ss) for _ in range(3))
            node, weight = self._centre_gauss
            for low, high in zip(cuts[:, :-1].T, cuts[:, 1:].T):
                active = high > low
                if not np.any(active):
                    continue
                half = (high[active]-low[active])/2
                z = low[active, None]+half[:, None]*(1+node)
                q = np.clip(ss[active, None]+dd[active, None]*z, self.domain["qmin"], self.domain["qb"])
                sample = self._many_profile(q, second=True)
                I1[active] += half*(sample["f_prime"] @ weight)
                IT[active] += half*(sample["f_second"] @ weight)
                I2[active] += half*((z*sample["f_second"]) @ weight)
            ratio = dd/rr
            # (r/(1-r²)-atanh(r))/r² without subtracting nearly equal terms.
            coefficient = rr*(2/3+rr*rr*(4/5+rr*rr*(6/7+rr*rr*8/9)))
            p[small] = ratio*I1
            p_T[small] = ratio*IT
            p_r[small] = coefficient*I1+ratio/(1-rr*rr)*I2
        self.budget.check()
        return {"p": p.reshape(shape), "p_T": p_T.reshape(shape), "p_r": p_r.reshape(shape),
                "p_center": (2*centre["f_prime"]).reshape(shape)}
