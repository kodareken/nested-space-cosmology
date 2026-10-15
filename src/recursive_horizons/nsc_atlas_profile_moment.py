"""Cached incoming-profile moment and pre-arrival cap boundary jets.

The original prepared profile supplies C(q)=integral(f_prime**2)/2. This
replaces repeated radial mass quadratures only where the centre is quiet.
It introduces no scalar evolution, source law, or numerical derivative.
"""
from collections import OrderedDict
import math

import numpy as np

from .nsc_atlas_bulk_response import BulkBudget
from .nsc_atlas_scalar_boundary import (
    DOMAIN_KEYS, MAX_QUERY, _order, _quadrature, _ray_domain, _shared_budget, _time,
    cap_potential,
)
from .nsc_atlas_scalar_lift import BoundaryJet, cap_potential_u, wall_boundary_derivative


SCHEMA = "nsc-atlas-profile-moment-v1"
MAX_CUTS = 1025
MAX_CACHE = 32768
MAX_SAMPLES = 129


def _array(values, name, maximum=MAX_QUERY):
    raw = np.asarray(values)
    if raw.dtype.kind not in "iuf" or raw.dtype.kind == "b":
        raise ValueError(name+" must contain finite real values")
    result = np.asarray(raw, dtype=float)
    if (result.ndim > 1 or not 1 <= result.size <= maximum
            or not np.all(np.isfinite(result))):
        raise ValueError(name+" requires a bounded finite scalar or one-dimensional array")
    return result


class ProfileMoment:
    """Prefix and partial-cell Gauss integration of the same first profile jet.

    Fixed cells include source edges, original field segment boundaries and
    arrival-spline knots. Arbitrary queries integrate only their one partial
    cell. Cached values are quadrature results, not an interpolated new state.
    """

    def __init__(self, cap_field, order=4, *, budget=None):
        if (isinstance(order, (bool, np.bool_)) or not isinstance(order, (int, np.integer))
                or order not in (4, 8)):
            raise ValueError("profile moment order must be 4 or 8")
        self.order, self.cap_field = int(order), cap_field
        shared = cap_field.budget
        if budget is None:
            budget = BulkBudget(started=shared.started, cpu_seconds=shared.cpu_seconds, external=shared)
        self.budget = _shared_budget(cap_field, budget)
        self.domain = {key: float(cap_field.domain[key]) for key in DOMAIN_KEYS}
        d = self.domain
        if (not all(math.isfinite(value) for value in d.values())
                or not d["Va"] < d["Vb"] or not d["qmin"] < d["Ta"] < d["Tb"] < d["qa"] < d["qb"]
                or cap_field.quiet_certificate.get("quiet_cap") is not True):
            raise ValueError("profile moment requires the original quiet pre-arrival domain")
        if getattr(cap_field.prepared_field, "budget", shared) is not shared:
            raise ValueError("profile moment requires the same prepared-field budget")
        self.geometry = cap_field.geometry
        arrival = self.geometry["arrival"]
        dense_cuts = _array(cap_field.prepared_field.dense_data["cuts"], "field segment cuts", MAX_CUTS)
        if (dense_cuts.ndim != 1 or dense_cuts.size < 2 or dense_cuts[0] != 0
                or dense_cuts[-1] != 1 or np.any(np.diff(dense_cuts) <= 0)):
            raise ValueError("field segment cuts must partition the original normalized interval")
        spline_knots = _array(arrival.x, "arrival knots", MAX_CUTS)
        edges = _array(cap_field.edges, "profile edges", 16)
        if (spline_knots.ndim != 1 or spline_knots.size < 2
                or np.any(np.diff(spline_knots) <= 0)
                or spline_knots[0] > d["Va"] or spline_knots[-1] < d["Vb"]):
            raise ValueError("arrival knots must cover the original source interval")
        source_V = np.clip(d["Va"]+(d["Vb"]-d["Va"])*dense_cuts, d["Va"], d["Vb"])
        mapped = np.asarray(arrival(source_V), dtype=float)
        knots = spline_knots[(spline_knots >= d["Va"]) & (spline_knots <= d["Vb"])]
        candidates = np.r_[d["qa"], d["qb"], edges, mapped, arrival(knots)]
        tolerance = 8*np.finfo(float).eps*max(1., abs(d["qa"]), abs(d["qb"]))
        if (not np.all(np.isfinite(candidates)) or np.any(np.diff(mapped) <= 0)
                or np.any(candidates < d["qa"]-tolerance)
                or np.any(candidates > d["qb"]+tolerance)):
            raise ValueError("profile cuts exit the original increasing arrival interval")
        self.cuts = np.unique(np.clip(candidates, d["qa"], d["qb"]))
        if self.cuts.size > MAX_CUTS:
            raise ValueError("profile moment exceeded its bounded cell count")
        self.gauss = np.polynomial.legendre.leggauss(self.order)
        self.cache = OrderedDict()
        self.profile_points = self.profile_batches = self.partial_nodes = self.cache_hits = 0
        self.value_queries = self.mass_points = 0
        half = np.diff(self.cuts)/2
        nodes, weights = self.gauss
        points = self.cuts[:-1, None]+half[:, None]*(1+nodes)
        data = self._profile(points.ravel())
        integrals = half*np.sum(weights*data["f_prime"].reshape(points.shape)**2, axis=1)/2
        self.prefix = np.r_[0., np.cumsum(integrals)]
        self.prefix_nodes = int(points.size)
        self.settings = {"moment_order": self.order}
        if not np.all(np.isfinite(self.prefix)) or np.any(np.diff(self.prefix) < 0):
            raise ValueError("profile moment is not a finite nondecreasing integral")
        self.budget.check()

    def _profile(self, labels):
        labels = np.asarray(labels, dtype=float)
        unique, inverse = np.unique(labels.ravel(), return_inverse=True)
        values = {key: np.empty(unique.size) for key in ("f", "f_prime")}
        for start in range(0, unique.size, MAX_QUERY):
            stop = min(start+MAX_QUERY, unique.size)
            current = unique[start:stop]
            self.budget.check(metric=current.size)
            row = self.cap_field.profile(current, second=False)
            for key in values:
                data = np.asarray(row[key], dtype=float)
                if data.shape != current.shape or not np.all(np.isfinite(data)):
                    raise ValueError("invalid original profile "+key)
                values[key][start:stop] = data
            self.profile_points += current.size
            self.profile_batches += 1
        return {key: value[inverse].reshape(labels.shape) for key, value in values.items()}

    def value(self, q):
        """C(q), with exact quiet extension on the original supplied interval."""
        q = _array(q, "moment labels")
        d = self.domain
        if np.any(q < d["qmin"]) or np.any(q > d["qb"]):
            raise ValueError("moment labels exit the original profile domain")
        self.budget.check()
        self.value_queries += q.size
        unique, inverse = np.unique(q.ravel(), return_inverse=True)
        values = np.zeros(unique.size)
        missing = []
        for i, label in enumerate(unique):
            if label <= d["qa"]:
                continue
            cached = self.cache.get(float(label))
            if cached is not None:
                values[i] = cached
                self.cache.move_to_end(float(label))
                self.cache_hits += 1
                continue
            cut = int(np.searchsorted(self.cuts, label))
            if cut < self.cuts.size and self.cuts[cut] == label:
                values[i] = self.prefix[cut]
            else:
                missing.append(i)
        if missing:
            labels = unique[missing]
            cells = np.searchsorted(self.cuts, labels, side="right")-1
            half = (labels-self.cuts[cells])/2
            nodes, weights = self.gauss
            points = self.cuts[cells, None]+half[:, None]*(1+nodes)
            data = self._profile(points.ravel())
            values[missing] = self.prefix[cells]+half*np.sum(
                weights*data["f_prime"].reshape(points.shape)**2, axis=1)/2
            self.partial_nodes += points.size
        if not np.all(np.isfinite(values)) or np.any(values < 0):
            raise ValueError("profile moment query is not finite and nonnegative")
        for label, value in zip(unique, values):
            if label > d["qa"]:
                self.cache[float(label)] = float(value)
                self.cache.move_to_end(float(label))
                if len(self.cache) > MAX_CACHE:
                    self.cache.popitem(last=False)
        self.budget.check()
        return values[inverse].reshape(q.shape)

    def mass(self, T, r):
        """Pre-arrival mass and algebraic first jets; no lapse or reflected era."""
        T, r = np.broadcast_arrays(_array(T, "cap time"), _array(r, "cap radius"))
        if T.size > MAX_QUERY:
            raise ValueError("request at most 513 cap moment endpoints")
        d = self.domain
        tolerance = 8*np.finfo(float).eps*max(1., abs(d["Ta"]), abs(d["qb"]))
        if np.any(T < d["Ta"]-tolerance) or np.any(T >= d["qa"]):
            raise ValueError("moment mass requires the known pre-arrival times T<qa")
        if np.any(r <= 0) or np.any(r >= 1):
            raise ValueError("moment mass requires 0<r<1")
        T = np.maximum(T, d["Ta"])
        limits = np.asarray(self.cap_field.radial_limit(T), dtype=float)
        if (limits.shape != T.shape or not np.all(np.isfinite(limits))
                or np.any(r > limits+tolerance)):
            raise ValueError("moment mass exits the retained finite cap")
        r = np.minimum(r, limits)
        q = T+np.arctanh(r)
        if np.any(q < d["qmin"]-tolerance) or np.any(q > d["qb"]+tolerance):
            raise ValueError("moment mass exits the original profile interval")
        q = np.clip(q, d["qmin"], d["qb"])
        C = self.value(q)
        profile = self._profile(q)
        f, fp = profile["f"], profile["f_prime"]
        active = q > d["qa"]
        if np.any(~active & ((f != 0) | (fp != 0))):
            raise ValueError("original profile violates its exact quiet extension")
        A = 1-r*r
        nu = C-np.divide(A*f*f, 4*r, out=np.zeros_like(r), where=active)
        p = np.divide(f, r, out=np.zeros_like(r), where=active)
        pT = np.divide(fp, r, out=np.zeros_like(r), where=active)
        pr = np.divide(fp/A-p, r, out=np.zeros_like(r), where=active)
        if not all(np.all(np.isfinite(value)) for value in (nu, p, pT, pr)):
            raise ValueError("pre-arrival mass or first jet is not finite")
        self.mass_points += T.size
        self.budget.check(metric=T.size)
        return {"T": T, "r": r, "q": q, "C": C, "f": f, "f_prime": fp,
                "nu_E": nu, "p": p, "p_T": pT, "p_r": pr,
                "p_center": np.zeros_like(T)}

    def summary(self):
        self.budget.check()
        return {"schema": SCHEMA, "domain": dict(self.domain), "settings": dict(self.settings),
                "cells": int(self.cuts.size-1), "prefix_nodes": self.prefix_nodes,
                "partial_nodes": int(self.partial_nodes), "profile_points": int(self.profile_points),
                "profile_batches": int(self.profile_batches), "value_queries": int(self.value_queries),
                "cache_hits": int(self.cache_hits), "cached_values": len(self.cache),
                "mass_points": int(self.mass_points), "terminal_C": float(self.prefix[-1])}


def outgoing_moment_jet(moment, q_minus, q_plus, *, ray_order=16):
    """Same outgoing q rule and moving-limit derivative, using cached mass."""
    ray_order = _order(ray_order, "ray order")
    cap, budget = moment.cap_field, moment.budget
    u, upper, margins = _ray_domain(cap, q_minus, q_plus)
    labels, radii, weights = _quadrature(u, moment.domain["qa"], upper,
                                       cap.edges, ray_order, "q")
    positive = negative = integral_u = 0.
    maxima = dict.fromkeys(("potential_identity_defect", "potential_u_identity_defect",
                           "mass_null_identity_defect"), 0.)
    for start in range(0, labels.size, MAX_QUERY):
        qs, rs, ws = [value[start:start+MAX_QUERY] for value in (labels, radii, weights)]
        row = moment.mass((u+qs)/2, rs)
        args = (rs, row["nu_E"], row["p"], row["p_T"], row["p_r"], row["f"], row["f_prime"])
        potential, derivative = cap_potential(*args), cap_potential_u(*args)
        source = -potential["K1"]*row["f"]
        positive += float(ws @ np.maximum(source, 0.))
        negative += float(ws @ np.minimum(source, 0.))
        integral_u += float(ws @ (derivative["K1_u"]*row["f"]))
        for key in maxima:
            data = potential[key] if key in potential else derivative[key]
            maxima[key] = max(maxima[key], float(np.max(abs(data))))
    radius = math.tanh((upper-u)/2)
    endpoint = moment.mass((u+upper)/2, radius)
    potential = cap_potential(radius, endpoint["nu_E"], endpoint["p"], endpoint["p_T"],
                              endpoint["p_r"], endpoint["f"], endpoint["f_prime"])
    maxima["potential_identity_defect"] = max(maxima["potential_identity_defect"],
                                               abs(float(potential["potential_identity_defect"])))
    result = {"J_minus": positive+negative, "J_minus_u": -integral_u,
              "integral_Ku_f": integral_u, "K1_wall": float(potential["K1"]),
              "wall_profile": float(endpoint["f"]), "positive_integral": positive,
              "negative_integral": negative, "ray_nodes": int(labels.size), **maxima, **margins}
    if not all(math.isfinite(float(value)) for value in result.values()):
        raise ValueError("moment boundary jet is not finite")
    budget.check()
    return result


class MomentBoundaryJet(BoundaryJet):
    """BoundaryJet-compatible rows with a separately declared moment order."""

    def __init__(self, boundary, *, moment_order=4, moment=None):
        super().__init__(boundary)
        if (isinstance(moment_order, (bool, np.bool_))
                or not isinstance(moment_order, (int, np.integer)) or moment_order not in (4, 8)):
            raise ValueError("boundary moment order must be 4 or 8")
        if moment is None:
            moment = ProfileMoment(self.cap_field, moment_order, budget=self.budget)
        if (not isinstance(moment, ProfileMoment) or moment.cap_field is not self.cap_field
                or moment.budget is not self.budget or moment.order != moment_order):
            raise ValueError("boundary jet requires the same profile moment, order and budget")
        self.moment, self.moment_order = moment, moment.order
        # radial_order is retained as the original boundary's comparison
        # setting; this implementation evaluates mass with moment_order.
        self.settings = {**self.settings, "moment_order": self.moment_order}

    def sample(self, t):
        t = _time(t)
        self.budget.check()
        if t in self.cache:
            return dict(self.cache[t])
        if len(self.cache) >= MAX_SAMPLES:
            raise ValueError("retain at most 129 moment boundary jet samples")
        wall = self.wall_solution.sample(t)
        V = float(wall["V"])
        q = float(self.geometry["arrival"](V))
        u = q-2*math.atanh(wall["x"])
        jet = outgoing_moment_jet(self.moment, u, q, ray_order=self.settings["ray_order"])
        x, v, bm = [float(wall[key]) for key in ("x", "v", "beta_minus")]
        qm_s, qp_s = (bm-v)/(1-x*x), (bm+v)/(1-x*x)
        Js = -qp_s*jet["K1_wall"]*jet["wall_profile"]-qm_s*jet["integral_Ku_f"]
        derivative = wall_boundary_derivative(wall, jet["J_minus"], Js)
        U, U_V = float(self.geometry["wall"](V)), float(self.geometry["wall"](V, 1))
        if not math.isfinite(U_V) or U_V <= 0:
            raise ValueError("moment boundary needs a positive wall derivative")
        stored = self.boundary.cache.get(t)
        row = {"t": t, "U": U, "V": V, "V0_s": float(wall["Vdot"]),
               "x": x, "q_minus": u, "q_plus": q, **jet, **derivative,
               "J_minus_s": Js, "U_V": U_V, "b_U": derivative["b_V"]/U_V,
               "arrival_rate_defect": float(self.geometry["arrival"](V, 1))*wall["Vdot"]-qp_s,
               "wall_coordinate_rate_defect": U_V*wall["Vdot"]-derivative["d"]/(2*derivative["r_U"]),
               "wall_profile_defect": jet["wall_profile"]-wall["w"],
               "boundary_cache_checked": stored is not None,
               "boundary_cache_value_defect": 0. if stored is None else derivative["c1_r"]-stored["c1_r"]}
        self.ray_nodes += jet["ray_nodes"]
        self.cache[t] = row
        return dict(row)
