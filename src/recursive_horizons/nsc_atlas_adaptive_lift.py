"""Fixed-cell local Hermite boundary polynomials and their continuous lift.

Each cell's stencil is fixed before queries and stays in one source block.
Its two endpoint values and slopes are shared with adjacent cells, giving a
C1 boundary and C2 antiderivative in exact arithmetic. Withheld queries assess
the numerical approximation; they are never inserted into these polynomials.
This module performs no physical preparation and preserves the provider budget.
"""
from math import factorial, fsum

import numpy as np
from numpy.polynomial import Polynomial
from scipy.interpolate import KroghInterpolator


SCHEMA = "nsc-atlas-adaptive-lift-v1"
SCOPE = "fixed-cell local Hermite lift; sampled sensitivity is not a uniform error enclosure"
BLOCKS = ((0., .25), (.25, .75), (.75, 1.))
MAX_SAMPLES = 129
MAX_QUERY = 513
ROW_KEYS = ("t", "U", "V", "b", "b_U", "b_V", "U_V")


def _array(value, name, *, maximum=MAX_QUERY):
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf" or raw.dtype.kind == "b":
        raise ValueError(name+" must contain finite real values")
    values = np.asarray(raw, dtype=float)
    if values.ndim > 1 or not 1 <= values.size <= maximum or not np.all(np.isfinite(values)):
        raise ValueError(name+" must be bounded, finite, scalar or one dimensional")
    return values


def _times(values, name):
    try:
        values = tuple(values)
    except TypeError as error:
        raise ValueError(name+" must be a sequence") from error
    if not values:
        return np.array([], dtype=float)
    if any(isinstance(value, (bool, np.bool_)) for value in values):
        raise ValueError(name+" must contain real times, not booleans")
    result = _array(values, name, maximum=MAX_SAMPLES)
    if result.ndim != 1 or np.any(result < 0) or np.any(result > 1):
        raise ValueError(name+" must lie in [0,1]")
    return result


def _scalar(value, name):
    values = _array(value, name)
    if values.ndim != 0:
        raise ValueError(name+" must be scalar")
    return float(values)


class AdaptiveBoundaryLift:
    """A local high-order lift sharing an explicit boundary-jet provider.

    ``jets`` supplies sample(t), a time-keyed cache, domain, geometry['wall'],
    and the original cumulative budget. Both BoundaryJet and MomentBoundaryJet
    satisfy this interface. Knot times must include 0, .25, .75 and 1, with at
    least ``nodes`` training points in each closed source block. A provider's
    combined training/validation cache may retain at most 129 direct samples.
    """

    def __init__(self, jets, *, knot_times, nodes=5, validation_times=()):
        if isinstance(nodes, (bool, np.bool_)) or not isinstance(nodes, (int, np.integer)) or nodes not in (4, 5):
            raise ValueError("local Hermite stencils require 4 or 5 nodes")
        self.jets, self.nodes = jets, int(nodes)
        self.budget, self.geometry = jets.budget, jets.geometry
        self.domain = {key: _scalar(jets.domain[key], key) for key in ("U0", "Ua", "Ub", "Va", "Vb")}
        if not (self.domain["U0"] < self.domain["Ua"] < self.domain["Ub"]
                and self.domain["Va"] < self.domain["Vb"]):
            raise ValueError("invalid finite boundary domain")
        self.knot_times = _times(knot_times, "knot times")
        if (len(self.knot_times) < 4 or np.any(np.diff(self.knot_times) <= 0)
                or not all(t in self.knot_times for t in (0., .25, .75, 1.))):
            raise ValueError("ordered knots must include both endpoints and source edges")
        for lo, hi in BLOCKS:
            if np.count_nonzero((self.knot_times >= lo) & (self.knot_times <= hi)) < self.nodes:
                raise ValueError("each source block needs enough training nodes for the fixed stencil")
        validation_times = self._validation_times(validation_times)
        self._reserve((*self.knot_times, *validation_times))
        self.rows = [self._sample(t) for t in self.knot_times]
        self.knot_count = len(self.rows)
        self.U_knots = np.array([row["U"] for row in self.rows])
        self.b_knots = np.array([row["b"] for row in self.rows])
        self.b_U_knots = np.array([row["b_U"] for row in self.rows])
        tolerance = self._endpoint_tolerance()
        if (abs(self.U_knots[0]-self.domain["Ua"]) > tolerance
                or abs(self.U_knots[-1]-self.domain["Ub"]) > tolerance):
            raise ValueError("training knots must cover the full finite U domain")
        self.U_knots[[0, -1]] = [self.domain["Ua"], self.domain["Ub"]]
        if np.any(np.diff(self.U_knots) <= 0):
            raise ValueError("training U coordinates must increase strictly")
        if abs(self.b_knots[0]) > 1e-24 or abs(self.b_U_knots[0]) > 1e-24:
            raise ValueError("lift requires a quiet value and slope anchor")
        self.b_knots[0] = self.b_U_knots[0] = 0.
        self.cells, self._pieces, self._primitives = [], [], []
        increments = []
        for index in range(self.knot_count-1):
            self._build_cell(index)
            increments.append(float(self._primitives[-1](1.)))
        self.prefix = np.array([0., *[fsum(increments[:i]) for i in range(1, len(increments)+1)]])
        if not np.all(np.isfinite(self.prefix)):
            raise ValueError("nonfinite integrated boundary polynomial")
        self.validation_rows = []
        self.validate(validation_times)
        self.budget.check()

    def _endpoint_tolerance(self):
        return 8*np.finfo(float).eps*max(1., abs(self.domain["Ua"]), abs(self.domain["Ub"]))

    def _reserve(self, times):
        self.budget.check()
        if len(set(self.jets.cache) | {float(t) for t in times}) > MAX_SAMPLES:
            raise ValueError("shared boundary cache may retain at most 129 direct samples")

    def _sample(self, t):
        self._reserve((t,))
        row = dict(self.jets.sample(float(t)))
        for key in ROW_KEYS:
            row[key] = _scalar(row[key], "boundary jet "+key)
        if row["t"] != t or row["U_V"] <= 0:
            raise ValueError("boundary jet time or wall derivative mismatch")
        V = self.domain["Va"]+(self.domain["Vb"]-self.domain["Va"])*t
        if abs(row["V"]-V) > 8*np.finfo(float).eps*max(1., abs(V)):
            raise ValueError("boundary jet does not use the original wall-time coordinate")
        self._reserve(())
        return row

    def _build_cell(self, index):
        left, right = self.U_knots[index:index+2]
        midpoint = left+(right-left)/2
        time_midpoint = (self.knot_times[index]+self.knot_times[index+1])/2
        block = next(i for i, (lo, hi) in enumerate(BLOCKS) if lo < time_midpoint < hi)
        lo, hi = BLOCKS[block]
        available = np.flatnonzero((self.knot_times >= lo) & (self.knot_times <= hi))
        chosen = sorted(sorted(available, key=lambda i: (abs(self.U_knots[i]-midpoint), int(i)))[:self.nodes])
        if index not in chosen or index+1 not in chosen:
            raise ValueError("a fixed local stencil must contain both cell endpoints")
        chosen = np.asarray(chosen, dtype=int)
        scale = (self.U_knots[chosen[-1]]-self.U_knots[chosen[0]])/2
        centre = self.U_knots[chosen[0]]+scale
        z = np.repeat((self.U_knots[chosen]-centre)/scale, 2)
        values = np.column_stack((self.b_knots[chosen], scale*self.b_U_knots[chosen])).ravel()
        fit = KroghInterpolator(z, values)
        width = right-left
        # Taylor conversion gives powers of eta=(U-left)/width, always in [0,1]
        # for this cell. This avoids ill-scaled powers of absolute U coordinates.
        derivatives = fit.derivatives((left-centre)/scale, der=2*self.nodes)
        coefficient = np.array([derivatives[k]*(width/scale)**k/factorial(k)
                                for k in range(2*self.nodes)])
        if not np.all(np.isfinite(coefficient)):
            raise ValueError("nonfinite confluent Hermite coefficients")
        polynomial = Polynomial(coefficient)
        value_defect = max(abs(float(polynomial(0.))-self.b_knots[index]),
                           abs(float(polynomial(1.))-self.b_knots[index+1]))
        slope_defect = max(abs(float(polynomial.deriv()(0.))-width*self.b_U_knots[index]),
                           abs(float(polynomial.deriv()(1.))-width*self.b_U_knots[index+1]))
        magnitude = max(1., float(np.max(abs(self.b_knots[chosen]))),
                        float(np.max(abs(width*self.b_U_knots[chosen]))))
        if max(value_defect, slope_defect) > 2e-10*magnitude:
            raise ValueError("local Hermite conversion lost its shared endpoint constraints")
        self._pieces.append(polynomial)
        self._primitives.append(width*polynomial.integ())
        self.cells.append({"index": index, "block": block,
            "U_left": float(left), "U_right": float(right), "width": float(width),
            "training_indices": chosen.tolist(), "training_times": self.knot_times[chosen].tolist(),
            "normalization_centre": float(centre), "normalization_scale": float(scale),
            "coefficients_in_cell_coordinate": coefficient.tolist(),
            "maximum_absolute_endpoint_value_defect": value_defect,
            "maximum_absolute_scaled_endpoint_slope_defect": slope_defect})

    def _U(self, U):
        values = _array(U, "lift U query")
        if np.any(values < self.domain["U0"]) or np.any(values > self.domain["Ub"]):
            raise ValueError("lift query exits the retained U domain")
        self.budget.check()
        return values

    def _evaluate(self, values, derivative):
        flat = values.ravel()
        output = np.zeros_like(flat)
        active = flat > self.domain["Ua"]
        cells = np.clip(np.searchsorted(self.U_knots, flat[active], side="right")-1, 0, self.knot_count-2)
        positions = np.flatnonzero(active)
        for index in np.unique(cells):
            selected = positions[cells == index]
            width = self.U_knots[index+1]-self.U_knots[index]
            eta = (flat[selected]-self.U_knots[index])/width
            if derivative == -1:
                output[selected] = self.prefix[index]+self._primitives[index](eta)
            else:
                output[selected] = self._pieces[index].deriv(derivative)(eta)/width**derivative
        return output.reshape(values.shape)

    def b_at_U(self, U):
        return self._evaluate(self._U(U), 0)

    def b_U(self, U):
        return self._evaluate(self._U(U), 1)

    def B(self, U):
        return self._evaluate(self._U(U), -1)

    def _wall_U(self, V):
        values = _array(V, "boundary V query")
        if np.any(values < self.domain["Va"]) or np.any(values > self.domain["Vb"]):
            raise ValueError("boundary query exits the finite V domain")
        U = np.asarray(self.geometry["wall"](values), dtype=float)
        tolerance = self._endpoint_tolerance()
        if (U.shape != values.shape or not np.all(np.isfinite(U))
                or np.any(U < self.domain["Ua"]-tolerance) or np.any(U > self.domain["Ub"]+tolerance)):
            raise ValueError("wall map exits the finite lift domain")
        return np.clip(U, self.domain["Ua"], self.domain["Ub"])

    def b(self, V):
        return self.b_at_U(self._wall_U(V))

    def B_at_wall(self, V):
        return self.B(self._wall_U(V))

    def _validation_times(self, times):
        values = _times(times, "validation times")
        if any(np.any(abs(self.knot_times-t) <= 4*np.finfo(float).eps) for t in values):
            raise ValueError("withheld query must not be a fitted knot")
        return values

    def validate(self, times):
        values = self._validation_times(times)
        self._reserve(values)
        for t in values:
            if any(row["t"] == t for row in self.validation_rows):
                continue
            row = self._sample(t)
            U, V = row["U"], row["V"]
            self._U(U)
            index = int(np.searchsorted(self.U_knots, U, side="right")-1)
            if not (0 <= index < self.knot_count-1
                    and self.knot_times[index] < t < self.knot_times[index+1]
                    and self.U_knots[index] < U < self.U_knots[index+1]):
                raise ValueError("withheld coordinate exits its training cell")
            width = self.U_knots[index+1]-self.U_knots[index]
            eta = (U-self.U_knots[index])/width
            value, slope = float(self.b_at_U(U)), float(self.b_U(U))
            derivative = float(self._primitives[index].deriv()(eta)/width)
            self.validation_rows.append({"t": float(t), "U": U, "V": V, "cell": index,
                "b_reference": row["b"], "b_interpolated": value, "b_defect": value-row["b"],
                "b_U_reference": row["b_U"], "b_U_interpolated": slope,
                "cell_width": float(width), "scaled_slope_defect": width*(slope-row["b_U"]),
                "B_lift": float(self.B(U)), "B_wall": float(self.B_at_wall(V)),
                "lift_derivative_defect": derivative-value,
                "wall_transport_defect": (derivative-value)*row["U_V"]})
        self.budget.check()
        return [dict(row) for row in self.validation_rows]

    def summary(self):
        self.budget.check()
        defects = ("b_defect", "scaled_slope_defect", "lift_derivative_defect", "wall_transport_defect")
        return {"schema": SCHEMA, "scope": SCOPE, "domain": dict(self.domain),
            "interpolation_variable": "U_W", "stencil_nodes": self.nodes,
            "stencil_selection": "nearest cell midpoint within source block; fixed for the whole cell",
            "source_blocks": [list(pair) for pair in BLOCKS], "knot_count": self.knot_count,
            "knot_times": self.knot_times.tolist(), "knot_samples": [dict(row) for row in self.rows],
            "cells": [dict(cell) for cell in self.cells],
            "validation_samples": [dict(row) for row in self.validation_rows],
            "withheld_evaluations": len(self.validation_rows),
            "maximum_absolute_defects": {key: max((abs(row[key]) for row in self.validation_rows), default=0.)
                                         for key in defects},
            "maximum_absolute_endpoint_value_defect": max(cell["maximum_absolute_endpoint_value_defect"] for cell in self.cells),
            "maximum_absolute_scaled_endpoint_slope_defect": max(cell["maximum_absolute_scaled_endpoint_slope_defect"] for cell in self.cells),
            "terminal_lift": float(self.B(self.domain["Ub"])),
            "quiet_anchor": {"b": float(self.b_at_U(self.domain["Ua"])),
                "b_U": float(self.b_U(self.domain["Ua"])), "B": float(self.B(self.domain["Ua"]))},
            "boundary_jet_evaluations": len(self.jets.cache), "extrapolated_queries": 0,
            "heldout_used_for_training": False}
