"""Checked time/active-radius cache for known scalar correction inputs.

Only F, G and the source pullback are interpolated. The evolved unknown
stays in the original U,V coordinates; composing this cache introduces no
advection term. Withheld comparisons are sampled sensitivities, not bounds.
"""
import math

import numpy as np
from scipy.interpolate import CubicHermiteSpline, CubicSpline

from .nsc_atlas_bulk_response import BulkBudget
from .nsc_atlas_characteristic_initialization import validate_geometry
from .nsc_atlas_scalar_forcing import ScalarForcing


SCHEMA = "nsc-atlas-forcing-cache-v1"
BLOCKS = ((0., .25), (.25, .75), (.75, 1.))
DOMAIN_KEYS = ("U0", "Ua", "Ub", "Va", "Vb")
MAX_SLICES = 129
QUIET_TOLERANCE = 1e-24
WALL_COORDINATE_REPLAY_ATOL = 2e-13
DEFECT_KEYS = ("temporal_F_defect", "temporal_G_defect", "combined_F_defect",
               "combined_G_defect", "source_F_defect", "g1_defect", "g1_V_defect",
               "scaled_g1_V_defect")


def _array(values, label, maximum):
    raw = np.asarray(values)
    if raw.dtype.kind not in "iuf" or raw.dtype.kind == "b":
        raise ValueError(label+" must contain finite real values")
    result = np.asarray(raw, dtype=float)
    if (result.ndim > 1 or not 1 <= result.size <= maximum
            or not np.all(np.isfinite(result))):
        raise ValueError(label+" must be bounded, finite and at most one dimensional")
    return result


def _time(t):
    value = _array(t, "normalized time", 1)
    if value.ndim != 0 or not 0 <= value <= 1:
        raise ValueError("normalized time must be a scalar in [0,1]")
    return float(value)


def _ordered(values, label, maximum, *, endpoints=False):
    values = _array(values, label, maximum).copy()
    if (values.ndim != 1 or values.size < 2 or np.any(np.diff(values) <= 32*np.finfo(float).eps)
            or values[0] < 0 or values[-1] > 1):
        raise ValueError(label+" must be strictly spaced and ordered in [0,1]")
    if endpoints and (values[0] != 0 or values[-1] != 1):
        raise ValueError(label+" must include exact endpoints zero and one")
    values.setflags(write=False)
    return values


class ForcingSliceProvider:
    """One shared physical-slice cache for nested training comparisons."""

    def __init__(self, wedge, source_geometry, *, budget, radial_order=2,
                 field_order=32, source_order=64):
        validate_geometry(source_geometry)
        self.wedge, self.source_geometry = wedge, dict(source_geometry)
        self.domain = {key: float(source_geometry[key]) for key in DOMAIN_KEYS}
        if any(not math.isclose(float(wedge.domain[key]), self.domain[key], rel_tol=0, abs_tol=1e-13)
               for key in DOMAIN_KEYS):
            raise ValueError("forcing slice domains differ")
        if (not isinstance(budget, BulkBudget) or budget.external is not wedge.budget
                or budget.started != wedge.budget.started):
            raise ValueError("forcing slices require the original shared budget")
        self.budget, self.geometry = budget, wedge.geometry
        self.settings = {"radial_order": radial_order, "field_order": field_order,
                         "source_order": source_order}
        for key, choices in (("radial_order", (2, 4)), ("field_order", (32, 64)), ("source_order", (32, 64))):
            value = self.settings[key]
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value not in choices:
                raise ValueError("invalid forcing slice "+key)
            self.settings[key] = int(value)
        self.cache = {}
        self.recorded_only = False
        self.cache_hits = 0
        self.counts = dict.fromkeys(("radial_nodes", "field_points", "field_batches",
                                    "kernel_values", "metric_values", "source_nodes"), 0)
        self.maxima = dict.fromkeys(("quiet_input_defect", "wall_flux_defect",
                                    "source_p_defect", "source_p_V_defect"), 0.)
        self.budget.check()

    @classmethod
    def from_records(cls, wedge, source_geometry, records, *, budget, radial_order=2,
                     field_order=32, source_order=64):
        """Replay retained coefficient slices, with no physical fallback.

        A lightweight wedge supplying only the original domain, saved wall
        coordinates and shared budget suffices. Every requested slice must
        have been retained; this method does not reconstruct a scalar state.
        """
        provider = cls(wedge, source_geometry, budget=budget, radial_order=radial_order,
                       field_order=field_order, source_order=source_order)
        for original in records:
            if len(provider.cache) >= MAX_SLICES:
                raise ValueError("retain at most 129 distinct forcing slices")
            row = dict(original)
            t = _time(row["t"])
            eta = _ordered(row["eta"], "recorded eta grid", 513, endpoints=True)
            key = (t, tuple(eta))
            if key in provider.cache:
                raise ValueError("duplicate recorded forcing slice")
            d = provider.domain
            V = min(d["Vb"], max(d["Va"], d["Va"]+(d["Vb"]-d["Va"])*t))
            wall = float(provider.geometry["wall"](V))
            for name in ("U", "F", "G"):
                value = _array(row[name], "recorded "+name, 513).copy()
                if value.shape != eta.shape:
                    raise ValueError("recorded forcing array shape differs from eta")
                row[name] = value
            # The rebuilt reference curve has the same coordinate replay
            # allowance as restore_inputs. The retained affine grid still has
            # its original, much tighter internal arithmetic check: accepting
            # curve roundoff must not accept independently displaced nodes.
            recorded_wall = float(row["U"][-1])
            expected = d["U0"]+eta*(recorded_wall-d["U0"])
            expected[0], expected[-1] = d["U0"], recorded_wall
            if (not math.isfinite(wall)
                    or abs(recorded_wall-wall) > WALL_COORDINATE_REPLAY_ATOL
                    or not np.allclose(row["U"], expected, rtol=0,
                                       atol=8*np.finfo(float).eps)):
                raise ValueError("recorded forcing coordinates differ from the original wall")
            for name in ("g1", "g1_V", "source_p_defect", "source_p_V_defect"):
                value = _array(row[name], "recorded "+name, 1)
                if value.ndim != 0:
                    raise ValueError("recorded forcing scalar is not scalar")
                row[name] = float(value)
            for name in ("radial_nodes", "field_points", "field_batches"):
                provider.counts[name] += BulkBudget._count(row[name])
            costs = dict(row["slice_costs"])
            if set(costs) != {"kernel_values", "metric_values", "source_nodes"}:
                raise ValueError("recorded forcing slice costs have unexpected keys")
            for name, value in costs.items():
                provider.counts[name] += BulkBudget._count(value)
            row["eta"], row["slice_costs"] = eta, costs
            provider._accept_row(t, row)
            provider.cache[key] = row
        provider.recorded_only = True
        return provider

    def _accept_row(self, t, row):
        for name in ("source_p_defect", "source_p_V_defect"):
            self.maxima[name] = max(self.maxima[name], abs(float(row[name])))
        wall_error = abs(float(row["F"][-1]))
        self.maxima["wall_flux_defect"] = max(self.maxima["wall_flux_defect"], wall_error)
        if wall_error > QUIET_TOLERANCE:
            raise ValueError("forcing slice violates the exact zero wall flux")
        if t <= .25:
            quiet = max(float(np.max(abs(row["F"]))), float(np.max(abs(row["G"]))),
                        abs(float(row["g1"])), abs(float(row["g1_V"])))
            self.maxima["quiet_input_defect"] = max(self.maxima["quiet_input_defect"], quiet)
            if quiet > QUIET_TOLERANCE:
                raise ValueError("forcing slice violates the original quiet interval")
        for value in row.values():
            if isinstance(value, np.ndarray):
                value.setflags(write=False)

    def slice(self, t, eta_grid):
        t = _time(t)
        eta = _ordered(eta_grid, "slice eta grid", 513, endpoints=True)
        key = (t, tuple(eta))
        self.budget.check()
        if key in self.cache:
            self.cache_hits += 1
            return dict(self.cache[key])
        if self.recorded_only:
            raise ValueError("requested forcing slice was not retained")
        if len(self.cache) >= MAX_SLICES:
            raise ValueError("retain at most 129 distinct forcing slices")
        d = self.domain
        V = min(d["Vb"], max(d["Va"], d["Va"]+(d["Vb"]-d["Va"])*t))
        wall = float(self.geometry["wall"](V))
        U = d["U0"]+eta*(wall-d["U0"])
        U[0], U[-1] = d["U0"], d["Ub"]
        before = (self.wedge.budget.kernel_values, self.budget.metric_values, self.budget.source_nodes)
        row = ScalarForcing(self.wedge, self.source_geometry, U, budget=self.budget,
                            **self.settings).slice(t)
        row["eta"] = eta.copy()
        after = (self.wedge.budget.kernel_values, self.budget.metric_values, self.budget.source_nodes)
        row["slice_costs"] = dict(zip(("kernel_values", "metric_values", "source_nodes"),
                                     (last-first for first, last in zip(before, after))))
        for name, first, last in zip(("kernel_values", "metric_values", "source_nodes"), before, after):
            self.counts[name] += last-first
        for name in ("radial_nodes", "field_points", "field_batches"):
            self.counts[name] += int(row[name])
        self._accept_row(t, row)
        self.cache[key] = row
        return dict(row)

    def summary(self):
        self.budget.check()
        return {"physical_slices": len(self.cache), "cache_hits": self.cache_hits,
                **self.counts, "maximum_absolute_input_defects": dict(self.maxima)}


class ForcingTimeCache:
    """Cubic time interpolation per source block, with linear active eta.

    Source g1 uses its supplied continuum first derivative as Hermite data.
    eval returns the derivative of that actual Hermite interpolant. Validation
    samples neither alter the fit nor its source derivative.
    """

    def __init__(self, wedge, source_geometry, eta_grid, training_times,
                 validation_times=(), *, budget, radial_order=2, field_order=32,
                 source_order=64, slice_provider=None):
        self.eta_grid = _ordered(eta_grid, "training eta grid", 257, endpoints=True)
        self.training_times = _ordered(training_times, "training times", MAX_SLICES, endpoints=True)
        self.time_knots = self.training_times
        if not all(edge in self.training_times for edge in (0., .25, .75, 1.)):
            raise ValueError("training times must retain both source edges")
        for lo, hi in BLOCKS:
            if np.count_nonzero((self.training_times >= lo) & (self.training_times <= hi)) < 4:
                raise ValueError("cubic time interpolation requires four knots in each source block")
        settings = dict(radial_order=radial_order, field_order=field_order, source_order=source_order)
        if slice_provider is None:
            slice_provider = ForcingSliceProvider(wedge, source_geometry, budget=budget, **settings)
        if (not isinstance(slice_provider, ForcingSliceProvider) or slice_provider.wedge is not wedge
                or slice_provider.budget is not budget or slice_provider.settings != settings
                or slice_provider.source_geometry != source_geometry):
            raise ValueError("time caches must reuse the same forcing slice provider and settings")
        self.slice_provider = slice_provider
        self.budget, self.domain = budget, dict(slice_provider.domain)
        self.geometry, self.source_geometry = wedge.geometry, dict(source_geometry)
        self.settings = {**settings, "eta_nodes": int(self.eta_grid.size),
                         "training_times": int(self.training_times.size)}
        self.dV = self.domain["Vb"]-self.domain["Va"]
        self.validation_eta = np.sort(np.r_[self.eta_grid, (self.eta_grid[:-1]+self.eta_grid[1:])/2])
        self.validation_rows = []
        self.eval_calls = self.eval_points = 0
        validation = self._validation_times(validation_times)
        requested = {(float(t), tuple(self.eta_grid)) for t in self.training_times}
        requested |= {(float(t), tuple(self.validation_eta)) for t in validation}
        if len(set(slice_provider.cache) | requested) > MAX_SLICES:
            raise ValueError("training and validation exceed the shared finite slice count")
        self.training_rows = [slice_provider.slice(float(t), self.eta_grid) for t in self.training_times]
        values = np.array([np.column_stack((row["F"], row["G"])) for row in self.training_rows])
        g1 = np.array([row["g1"] for row in self.training_rows])
        g1V = np.array([row["g1_V"] for row in self.training_rows])
        quiet = self.training_times <= .25
        values[quiet], g1[quiet], g1V[quiet] = 0., 0., 0.
        values[:, -1, 0] = 0.
        self._coefficients, self._sources = [], []
        for lo, hi in BLOCKS:
            selected = (self.training_times >= lo) & (self.training_times <= hi)
            times = self.training_times[selected]
            self._coefficients.append(CubicSpline(times, values[selected], axis=0, extrapolate=False))
            self._sources.append(CubicHermiteSpline(times, g1[selected], self.dV*g1V[selected], extrapolate=False))
        self.validate(validation)

    def _validation_times(self, times):
        raw = tuple(times)
        if not raw:
            return np.array([], dtype=float)
        values = _array(raw, "withheld times", MAX_SLICES)
        if (values.ndim != 1 or np.any(values < 0) or np.any(values > 1)
                or np.any(np.diff(values) <= 32*np.finfo(float).eps)):
            raise ValueError("withheld times must increase within the finite interval")
        if any(np.any(abs(self.training_times-t) <= 4*np.finfo(float).eps) for t in values):
            raise ValueError("withheld slices must not be fitted training knots")
        return values

    def eval(self, U, t):
        """Interpolate known coefficients only; never evaluate a physical slice."""
        t, U = _time(t), _array(U, "cache U query", 513)
        self.budget.check()
        d = self.domain
        V = min(d["Vb"], max(d["Va"], d["Va"]+self.dV*t))
        wall = float(self.geometry["wall"](V))
        tolerance = 8*np.finfo(float).eps*max(1., abs(d["Ub"]))
        if (not math.isfinite(wall) or not d["U0"] < wall <= d["Ub"]+tolerance
                or np.any(U < d["U0"]) or np.any(U > max(d["Ub"], wall))):
            raise ValueError("cache query exits the original retained U domain")
        # Ub denotes the moving terminal node, including a one-ulp wall
        # endpoint outside the stored coordinate, as in ScalarForcing.
        clipped = np.where(U == d["Ub"], wall, np.minimum(U, wall))
        eta = np.clip((clipped-d["U0"])/(wall-d["U0"]), 0., 1.)
        block = int(np.searchsorted([.25, .75], t, side="right"))
        data = self._coefficients[block](t)
        Fsource = float(data[0, 0])
        F = np.interp(eta, self.eta_grid, data[:, 0])
        G = np.interp(eta, self.eta_grid, data[:, 1])
        F = np.where(clipped == d["U0"], Fsource, np.where(clipped == wall, 0., F))
        g1, g1V = float(self._sources[block](t)), float(self._sources[block](t, 1)/self.dV)
        if t <= .25:
            F, G, Fsource, g1, g1V = np.zeros_like(clipped), np.zeros_like(clipped), 0., 0., 0.
        if not all(np.all(np.isfinite(value)) for value in (F, G, Fsource, g1, g1V)):
            raise ValueError("forcing cache interpolation is not finite")
        self.eval_calls += 1
        self.eval_points += U.size
        return {"U": clipped, "F": F, "G": G, "Fsource": Fsource, "g1": g1, "g1_V": g1V}

    def validate(self, times):
        times = self._validation_times(times)
        keys = {(float(t), tuple(self.validation_eta)) for t in times}
        if len(set(self.slice_provider.cache) | keys) > MAX_SLICES:
            raise ValueError("withheld queries exceed the shared finite slice count")
        for t in times:
            if any(row["t"] == t for row in self.validation_rows):
                continue
            reference = self.slice_provider.slice(float(t), self.validation_eta)
            fitted = self.eval(reference["U"], float(t))
            difference = {key: fitted[key]-reference[key] for key in ("F", "G")}
            row = {"t": float(t), "block": int(np.searchsorted([.25, .75], t, side="right")),
                   "source_F_defect": fitted["Fsource"]-float(reference["F"][0]),
                   "g1_defect": fitted["g1"]-reference["g1"],
                   "g1_V_defect": fitted["g1_V"]-reference["g1_V"],
                   "reference": reference, "interpolated": fitted}
            cell = int(np.searchsorted(self.training_times, t, side="right")-1)
            row["cell_width_V"] = self.dV*float(self.training_times[cell+1]-self.training_times[cell])
            row["scaled_g1_V_defect"] = row["cell_width_V"]*row["g1_V_defect"]
            for key in ("F", "G"):
                row["temporal_"+key+"_defect"] = float(np.max(abs(difference[key][::2])))
                row["combined_"+key+"_defect"] = float(np.max(abs(difference[key][1::2])))
            self.validation_rows.append(row)
        self.budget.check()
        return [dict(row) for row in self.validation_rows]

    def summary(self):
        self.budget.check()
        blocks = []
        for i, (lo, hi) in enumerate(BLOCKS):
            rows = [row for row in self.validation_rows if row["block"] == i]
            blocks.append({"interval": [lo, hi], "withheld_slices": len(rows),
                "maximum_absolute_defects": {key: max((abs(row[key]) for row in rows), default=0.)
                                               for key in DEFECT_KEYS}})
        return {"schema": SCHEMA, "domain": dict(self.domain), "settings": dict(self.settings),
                "training_times": self.training_times.tolist(), "eta_grid": self.eta_grid.tolist(),
                "validation_times": [row["t"] for row in self.validation_rows],
                "time_interpolation": "cubic within source blocks; Hermite source g1",
                "radial_interpolation": "linear in active eta; evaluated back at original U",
                "midpoint_defect_scope": "combined time and radial interpolation sensitivity",
                "uniform_error_enclosure": False, "withheld_used_for_training": False,
                "blocks": blocks, "eval_calls": self.eval_calls, "eval_points": int(self.eval_points),
                "slice_provider": self.slice_provider.summary()}
