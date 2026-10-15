"""First scalar jets on the retained moving-wall wedge of one prepared solve.

The dense scalar evolution is borrowed unchanged. Temporal derivative
quadrature starts at the bottom or the birth of the column on the wall;
no exterior kernel or interpolation is continued into the cap.
"""
import math

import numpy as np

from .nsc_atlas_field_preparation import (
    M_OVER_ELL, exterior_coefficient, radius_from_null, source_profile,
)
from .nsc_atlas_prepared_response import MAX_FIELD_NODES, _minimum_derivative


class WedgeField:
    """Query the same numerical field on U0 <= U <= Uw(V).

    ``bulk(U, t, order=32|64)`` uses normalized t=(V-Va)/(Vb-Va), but
    returns derivatives with respect to U,V. ``chi`` denotes chi/epsilon
    and ``p`` denotes ell*phi/epsilon. ``chi_U_lower_limit_bound`` multiplies
    the numerically located inverse bracket width by the analytic continuum
    integrand bound. This conditional omitted-sliver estimate excludes
    spline approximation and floating-point wall evaluation errors. Neither
    the inverse bracket nor the numerical field is a validated enclosure.
    """

    def __init__(self, prepared_field):
        self.prepared_field = prepared_field
        self.geometry = prepared_field.geometry
        self.domain = prepared_field.domain
        self.budget = prepared_field.budget
        dense = prepared_field.dense_data
        self.grid = np.asarray(dense["grid"], dtype=float)
        self.cuts = np.asarray(dense["cuts"], dtype=float)
        self.segments = tuple(dense["segments"])
        self.U0, self.Ua, self.Ub, self.Va, self.Vb = (
            float(self.domain[key]) for key in ("U0", "Ua", "Ub", "Va", "Vb"))
        if (not all(math.isfinite(v) for v in
                    (self.U0, self.Ua, self.Ub, self.Va, self.Vb))
                or not self.U0 < self.Ua < self.Ub < 0
                or not self.Va < self.Vb < 0):
            raise ValueError("wedge requires the ordered finite IV preparation domain")
        if (self.grid.ndim != 1 or not 2 <= self.grid.size <= MAX_FIELD_NODES
                or not np.all(np.isfinite(self.grid)) or np.any(np.diff(self.grid) <= 0)
                or self.grid[0] != self.U0 or self.grid[-1] != self.Ub
                or self.cuts.ndim != 1 or len(self.cuts) != len(self.segments)+1
                or len(self.cuts) < 2 or self.cuts[0] != 0 or self.cuts[-1] != 1
                or not np.all(np.isfinite(self.cuts)) or np.any(np.diff(self.cuts) <= 0)
                or not all(callable(segment) for segment in self.segments)):
            raise ValueError("dense field grid or segment partition is malformed")
        wall = self.geometry["wall"]
        if (not hasattr(wall, "x") or not hasattr(wall, "c")
                or wall.x[0] > self.Va or wall.x[-1] < self.Vb
                or _minimum_derivative(wall) <= 0):
            raise ValueError("wall spline must be monotone and cover the whole preparation")
        ends = np.asarray(wall([self.Va, self.Vb]))
        if not np.allclose(ends, [self.Ua, self.Ub], rtol=0, atol=8*np.finfo(float).eps):
            raise ValueError("wall spline endpoints disagree with the preparation")
        self.dV = self.Vb-self.Va
        self.gauss = {n: np.polynomial.legendre.leggauss(n) for n in (32, 64)}
        qstar = float(dense["prepared"]["operator"]["analytic_norm_bound"])
        if not math.isfinite(qstar) or not 0 <= qstar < 1:
            raise ValueError("inverse-wall error needs the admitted contraction bound")
        self.budget.check(1)
        self._missing_integrand_bound = (
            float(exterior_coefficient(radius_from_null(self.U0, self.Va)))/(1-qstar))

    def _wall(self, V):
        V = np.asarray(V, dtype=float)
        if not np.all(np.isfinite(V)) or np.any(V < self.Va) or np.any(V > self.Vb):
            raise ValueError("wall evaluation exits its spline domain")
        values = np.asarray(self.geometry["wall"](V), dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError("wall spline evaluation is not finite")
        return values

    def _retained(self, U, V):
        U, V = np.broadcast_arrays(U, V)
        wall = self._wall(V)
        if (not np.all(np.isfinite(U)) or np.any(U < self.U0)
                or np.any(U > wall)):
            raise ValueError("query would evaluate the exterior on the excised cap side")
        return wall

    def _kernel(self, U, V):
        U, V = np.broadcast_arrays(U, V)
        self._retained(U, V)
        self.budget.check(U.size)
        return exterior_coefficient(radius_from_null(U, V))

    def _inverse_upper(self, U, upper_V):
        """Bracket the floating spline inverse with at most 56 bisections.

        The upper endpoint passes the same retained-side floating comparison
        used before kernel calls. This is not validated interval arithmetic.
        """
        U = np.asarray(U, dtype=float)
        upper = np.full_like(U, float(upper_V))
        lower = np.full_like(U, self.Va)
        if np.any(U < self.Ua) or np.any(U > self._wall(upper)):
            raise ValueError("inverse-wall target is outside its finite bracket")
        exact_bottom = U == self.Ua
        upper[exact_bottom] = self.Va
        for _ in range(56):
            midpoint = lower+.5*(upper-lower)
            active = (midpoint > lower) & (midpoint < upper)
            if not np.any(active):
                break
            below = self._wall(midpoint) < U
            lower = np.where(active & below, midpoint, lower)
            upper = np.where(active & ~below, midpoint, upper)
        self.budget.check()
        if np.any(self._wall(upper) < U):
            raise RuntimeError("inverse wall failed to preserve its retained-side bracket")
        widths = upper-lower
        if np.any(widths > 4*np.finfo(float).eps*max(abs(self.Va), abs(self.Vb))):
            raise RuntimeError("inverse wall did not resolve its finite bracket")
        return upper, widths

    def _paired_values(self, points, V):
        """Interpolate paired (U,V), avoiding a Cartesian time/space product."""
        points, V = np.broadcast_arrays(points, V)
        shape = points.shape
        points, V = points.ravel(), V.ravel()
        wall = self._retained(points, V)
        result = np.empty_like(points)
        # The largest dense temporary is at most 513 x 1024 values.
        for start in range(0, points.size, 1024):
            stop = min(start+1024, points.size)
            u, v, w = points[start:stop], V[start:stop], wall[start:stop]
            times = np.clip((v-self.Va)/self.dV, 0., 1.)
            indices = np.minimum(np.searchsorted(self.cuts, times, side="right")-1,
                                 len(self.segments)-1)
            correction = np.empty((self.grid.size, len(times)))
            for index in np.unique(indices):
                self.budget.check()
                mask = indices == index
                correction[:, mask] = self.segments[index](times[mask])
            if not np.all(np.isfinite(correction)):
                raise ValueError("dense field returned nonfinite values")
            cells = np.minimum(np.searchsorted(self.grid, u, side="right")-1,
                               self.grid.size-2)
            columns = np.arange(len(times))
            left = self.grid[cells]
            right = np.minimum(self.grid[cells+1], w)
            at_wall = u == w
            left_value = correction[cells, columns]
            right_value = np.where(self.grid[cells+1] < w,
                                   correction[cells+1, columns], correction[-1])
            fraction = np.divide(u-left, right-left, out=np.zeros_like(u),
                                 where=right > left)
            value = left_value+fraction*(right_value-left_value)
            value = np.where(at_wall, correction[-1], value)
            result[start:stop] = value+source_profile(v, self.Va, self.Vb)[0]
        self.budget.check()
        return result.reshape(shape)

    def bulk(self, U, normalized_V, *, order=64):
        points = np.asarray(U, dtype=float)
        if (points.ndim != 1 or not 1 <= points.size <= MAX_FIELD_NODES
                or not np.all(np.isfinite(points))):
            raise ValueError("wedge query needs 1 to 513 finite U values")
        if (isinstance(normalized_V, (bool, np.bool_))
                or not isinstance(normalized_V, (int, float, np.integer, np.floating))
                or not math.isfinite(normalized_V) or not 0 <= normalized_V <= 1):
            raise ValueError("wedge query exits the finite preparation window")
        if (isinstance(order, (bool, np.bool_)) or not isinstance(order, (int, np.integer))
                or order not in (32, 64)):
            raise ValueError("wedge derivative quadrature uses Gauss order 32 or 64")
        V = min(self.Vb, max(self.Va, self.Va+self.dV*float(normalized_V)))
        wall = float(self._wall(V))
        self._retained(points, V)
        clipped = np.minimum(self.grid, wall)
        chi_grid = self._paired_values(clipped, np.full_like(clipped, V))
        chi = np.interp(points, clipped, chi_grid)
        integrand = self._kernel(clipped, V)*chi_grid
        prefix = np.r_[0., np.cumsum(.5*(integrand[:-1]+integrand[1:])*np.diff(clipped))]
        indices = np.minimum(np.searchsorted(clipped, points, side="right")-1,
                             self.grid.size-2)
        partial = .5*(integrand[indices]+np.interp(points, clipped, integrand))*(points-clipped[indices])
        chi_V = float(source_profile(V, self.Va, self.Vb)[1])-prefix[indices]-partial

        lower = np.full_like(points, self.Va)
        width = np.zeros_like(points)
        born = points > self.Ua
        if np.any(born):
            lower[born], width[born] = self._inverse_upper(points[born], V)
        # Exact wall boundary condition, including a wall not on a grid node.
        lower[points == wall] = V
        width[points == wall] = 0.
        crossing = np.full_like(points, V)
        next_indices = np.minimum(np.searchsorted(self.grid, points, side="right"),
                                  self.grid.size-1)
        next_U = self.grid[next_indices]
        crosses = (next_U > self.Ua) & (next_U < wall) & (next_U > points)
        if np.any(crosses):
            crossing[crosses], _ = self._inverse_upper(next_U[crosses], V)
        breaks = np.sort(np.column_stack((lower,
            np.clip(np.full_like(points, self.Va+.25*self.dV), lower, V),
            np.clip(np.full_like(points, self.Va+.75*self.dV), lower, V),
            np.maximum(lower, crossing), np.full_like(points, V))), axis=1)
        chi_U = np.zeros_like(points)
        node, weight = self.gauss[order]
        for first, last in zip(breaks[:, :-1].T, breaks[:, 1:].T):
            active = last > first
            if not np.any(active):
                continue
            half = .5*(last[active]-first[active])
            null_V = first[active, None]+half[:, None]*(1+node)
            pair_U = np.broadcast_to(points[active, None], null_V.shape)
            values = self._paired_values(pair_U, null_V)
            chi_U[active] -= half*((self._kernel(pair_U, null_V)*values) @ weight)

        x = radius_from_null(points, V)
        inverse_Fprime = -4*M_OVER_ELL**2*np.exp(-x/(2*M_OVER_ELL))/x
        x_U, x_V = V*inverse_Fprime, points*inverse_Fprime
        p = chi/x
        self.budget.check()
        return {"U": points.copy(), "V": np.full_like(points, V),
                "chi": chi, "chi_U": chi_U, "chi_V": chi_V,
                "x": x, "x_U": x_U, "x_V": x_V,
                "p": p, "p_U": (chi_U-p*x_U)/x, "p_V": (chi_V-p*x_V)/x,
                "lower_V": lower, "lower_V_bracket_width": width,
                "chi_U_lower_limit_bound": width*self._missing_integrand_bound}
