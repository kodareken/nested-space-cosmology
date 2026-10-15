"""Cumulative exterior forcing for the unevolved scalar metric correction.

Each slice reuses the original prepared field and source gauge. One radial
pass supplies Einstein mass/lapse at all requested coordinates, then the
conservative flux F and density G. This adapter does not evolve the cubic
scalar coefficient, interpolate in time, or establish uniform error bounds.
"""
import math

import numpy as np

from .nsc_atlas_bulk_response import A_RADIUS, M, BulkBudget, background, source_traces
from .nsc_atlas_characteristic_initialization import validate_geometry


SCHEMA = "nsc-atlas-scalar-forcing-v1"
DOMAIN_KEYS = ("U0", "Ua", "Ub", "Va", "Vb")
MAX_QUERY = 513
FIELD_KEYS = ("p", "p_U", "p_V", "chi", "chi_U", "chi_V",
              "chi_U_lower_limit_bound")


def _number(value, label):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(label+" must be a finite real number")
    try:
        value = float(value)
    except (ValueError, TypeError, OverflowError) as error:
        raise ValueError(label+" must be a finite real number") from error
    if not math.isfinite(value):
        raise ValueError(label+" must be a finite real number")
    return value


def _choice(value, choices, label):
    if (isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, np.integer)) or value not in choices):
        raise ValueError(label+" must be one of "+str(choices))
    return int(value)


def _grid(values, lower, upper, label):
    try:
        grid = np.array(values, dtype=float, copy=True)
    except (ValueError, TypeError) as error:
        raise ValueError(label+" must be an ordered finite U grid") from error
    if (grid.ndim != 1 or not 2 <= grid.size <= MAX_QUERY
            or not np.all(np.isfinite(grid)) or np.any(np.diff(grid) <= 0)
            or grid[0] != lower or grid[-1] != upper):
        raise ValueError(label+" must span U0 to Ub with 2 to 513 ordered nodes")
    grid.setflags(write=False)
    return grid


class ScalarForcing:
    """Evaluate bounded conservative forcing slices on a fixed U grid.

    ``t`` is normalized exterior V time. Returned derivatives are in U,V,
    not t. Grid nodes beyond the instantaneous wall are represented by the
    actual wall value; no exterior jet is queried in the excised cap.
    ``source_scale=0`` is solely the quiet software control, as in the
    existing source adapter. The physical preparation has scale one.
    """

    def __init__(self, wedge, geometry, U_grid, *, radial_order=2, field_order=32,
                 source_order=64, source_scale=1., budget=None):
        validate_geometry(geometry)
        self.geometry = dict(geometry)
        self.domain = {key: float(geometry[key]) for key in DOMAIN_KEYS}
        try:
            same_domain = all(math.isclose(float(wedge.domain[key]), self.domain[key],
                                          rel_tol=0, abs_tol=1e-13)
                              for key in DOMAIN_KEYS)
            shared = wedge.budget
            wall = wedge.geometry["wall"]
        except (AttributeError, KeyError, TypeError, ValueError) as error:
            raise ValueError("forcing requires a prepared wedge and its finite domain") from error
        if not same_domain or not callable(getattr(wedge, "bulk", None)) or not callable(wall):
            raise ValueError("forcing field and source domains do not agree")
        self.wedge = wedge
        self.grid = _grid(U_grid, self.domain["U0"], self.domain["Ub"], "forcing grid")
        self.field_grid = _grid(wedge.grid, self.domain["U0"], self.domain["Ub"], "field grid")
        self.radial_order = _choice(radial_order, (2, 4), "radial order")
        self.field_order = _choice(field_order, (32, 64), "field derivative order")
        self.source_order = _choice(source_order, (32, 64), "source order")
        self.source_scale = _number(source_scale, "source scale")
        if self.source_scale not in (0., 1.):
            raise ValueError("source scale is the physical 1 or quiet control 0")
        if budget is None:
            budget = BulkBudget(started=shared.started, cpu_seconds=shared.cpu_seconds,
                                external=shared)
        if (not isinstance(budget, BulkBudget) or budget.external is not shared
                or budget.started != shared.started):
            raise ValueError("forcing requires the field's shared compute budget and start")
        prepared = getattr(wedge, "prepared_field", None)
        if prepared is not None and getattr(prepared, "budget", shared) is not shared:
            raise ValueError("wedge and prepared field must share the same budget")
        self.budget = budget
        self.budget.check()
        endpoints = np.asarray(wall([self.domain["Va"], self.domain["Vb"]]), dtype=float)
        if (endpoints.shape != (2,) or not np.all(np.isfinite(endpoints))
                or not np.allclose(endpoints, [self.domain["Ua"], self.domain["Ub"]],
                                   rtol=0, atol=8*np.finfo(float).eps)):
            raise ValueError("wall endpoints disagree with the prepared domain")
        self.gauss = np.polynomial.legendre.leggauss(self.radial_order)
        self.settings = dict(radial_order=self.radial_order, field_order=self.field_order,
                             source_order=self.source_order, source_scale=self.source_scale,
                             grid_nodes=self.grid.size)
        self.slice_evaluations = 0
        self.field_points = 0
        self.field_batches = 0

    def _field(self, points, t):
        rows = []
        for first in range(0, points.size, MAX_QUERY):
            current = points[first:first+MAX_QUERY]
            self.budget.check(metric=current.size)
            raw = self.wedge.bulk(current, t, order=self.field_order)
            row = {}
            for key in FIELD_KEYS:
                if key not in raw:
                    raise ValueError("wedge field is missing "+key)
                value = np.asarray(raw[key], dtype=float)
                if value.shape != current.shape or not np.all(np.isfinite(value)):
                    raise ValueError("invalid wedge field jet: "+key)
                row[key] = value
            if np.any(row["chi_U_lower_limit_bound"] < 0):
                raise ValueError("inverse-wall sliver bounds must be nonnegative")
            rows.append(row)
            self.field_points += current.size
            self.field_batches += 1
        return {key: np.concatenate([row[key] for row in rows]) for key in FIELD_KEYS}

    def slice(self, t):
        """Return one transient forcing slice; no time interpolation or PDE solve."""
        t = _number(t, "normalized time")
        if not 0 <= t <= 1:
            raise ValueError("forcing slice exits the finite preparation window")
        self.budget.check()
        d = self.domain
        V = min(d["Vb"], max(d["Va"], d["Va"]+(d["Vb"]-d["Va"])*t))
        wall_value = np.asarray(self.wedge.geometry["wall"](V), dtype=float)
        if wall_value.shape != ():
            raise ValueError("wall coordinate must be a scalar")
        wall = _number(float(wall_value), "wall coordinate")
        endpoint_tolerance = 8*np.finfo(float).eps*max(1., abs(d["Ub"]))
        if not d["U0"] < wall <= d["Ub"]+endpoint_tolerance:
            raise ValueError("wall exits the retained forcing domain")
        clipped = np.minimum(self.grid, wall)
        # The terminal grid entry always represents the wall. The inherited
        # spline can end one ulp above Ub; keep its actual retained boundary
        # so the leading field still evaluates the exact chi_Uwall=0 trace.
        clipped[-1] = wall
        native = self.field_grid[(self.field_grid > d["U0"]) & (self.field_grid < wall)]
        cuts = np.unique(np.r_[d["U0"], native, clipped, wall])
        half = np.diff(cuts)/2
        node, weight = self.gauss
        gauss = cuts[:-1, None]+half[:, None]*(1+node)
        weights = half[:, None]*weight
        points = np.unique(np.r_[cuts, gauss.ravel()])
        scalar = self._field(points, t)
        bg = background(points, V)
        for key in ("r", "r_U", "r_V", "C"):
            value = np.asarray(bg[key], dtype=float)
            if value.shape != points.shape or not np.all(np.isfinite(value)) or np.any(value <= 0):
                raise ValueError("forcing requires positive finite retained-IV "+key)
            bg[key] = value
        source = source_traces(V, self.geometry, order=self.source_order,
                               source_scale=self.source_scale, budget=self.budget)
        required = ("p", "p_V", "mu_E", "zeta_E", "H_E", "H_E_V")
        source = {key: _number(source[key], "source "+key) for key in required}
        source_defects = np.array([scalar["p"][0]-source["p"],
                                   scalar["p_V"][0]-source["p_V"]])
        if not np.allclose([scalar["p"][0], scalar["p_V"][0]],
                           [source["p"], source["p_V"]], rtol=2e-9, atol=2e-10):
            raise ValueError("prepared field is incompatible with the fixed source gauge/data")
        ig = np.searchsorted(points, gauss)
        r, rU, rV, C = (bg[key][ig] for key in ("r", "r_U", "r_V", "C"))
        pU = scalar["p_U"][ig]
        mass_cells = np.sum(weights*(-r*r*rV*pU*pU/C), axis=1)
        lapse_cells = np.sum(weights*(r*pU*pU/(2*rU)), axis=1)
        mass_cuts = source["mu_E"]+np.r_[0., np.cumsum(mass_cells)]
        lapse_cuts = source["zeta_E"]+np.r_[0., np.cumsum(lapse_cells)]
        targets = np.searchsorted(cuts, clipped)
        jet_targets = np.searchsorted(points, clipped)
        mu, zeta = mass_cuts[targets], lapse_cuts[targets]
        r, rU, C = (bg[key][jet_targets] for key in ("r", "r_U", "C"))
        chi, chiU, pU, pV = (scalar[key][jet_targets]
                            for key in ("chi", "chi_U", "p_U", "p_V"))
        Z0, A0 = -2*A_RADIUS/V, 1-2*M/r
        K0 = C*M/(2*r**3)
        flux = Z0*(A0*zeta-2*mu/r)*chiU/(2*rU)
        density = (C*(M*zeta+mu)/(2*r**3)-pU*pV/6)*chi
        source_chi_r = chiU[0]/rU[0]
        P0r = Z0*2*M/r[0]**2
        g1 = -source["H_E"]*source_chi_r
        g1V = -(source["H_E_V"]+source["H_E"]*P0r/2)*source_chi_r
        g1V += source["H_E"]*K0[0]*chi[0]/rU[0]
        if not all(np.all(np.isfinite(value)) for value in (mu, zeta, K0, flux, density, g1, g1V)):
            raise ValueError("scalar forcing is not finite")
        self.slice_evaluations += 1
        self.budget.check()
        return {"t": t, "V": V, "U_wall": wall, "U": clipped,
                "K0": K0, "F": flux, "G": density, "mu_E": mu, "zeta_E": zeta,
                "g1": float(g1), "g1_V": float(g1V),
                "source_mu_E": source["mu_E"], "source_zeta_E": source["zeta_E"],
                "source_H_E": source["H_E"], "source_H_E_V": source["H_E_V"],
                "source_chi_r": float(source_chi_r),
                "source_chi": float(chi[0]), "source_r_U": float(rU[0]),
                "source_P0_r": float(P0r),
                "source_p_defect": float(source_defects[0]),
                "source_p_V_defect": float(source_defects[1]),
                "wall_chi_U": float(chiU[-1]),
                "chi_U_lower_limit_bound": float(np.max(scalar["chi_U_lower_limit_bound"])),
                "radial_panels": half.size, "radial_nodes": gauss.size,
                "field_points": points.size,
                "field_batches": (points.size+MAX_QUERY-1)//MAX_QUERY}
