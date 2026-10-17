"""Retarded exterior Dirichlet-to-Neumann map on the quiet IV wedge.

This is the classical causal response of one angular mode of the conformal
scalar. It is not a quantum state, a covariance, or a cosmological term.
Modes l>3 lie outside the inherited contraction certificate. The bare
equation has every integer l; this finite operator is not a 4D determinant.

Coordinates, chi=r*phi, and the returned derivative use dimensionless
ell=1. phi=chi/r matches the existing ell*phi normalization. The outward
derivative is per dimensionless proper length along nu+.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .nsc_atlas_field_preparation import (
    M_OVER_ELL, RADIAL_FLOOR, exterior_coefficient, radius_from_null,
    radius_product,
)
from .nsc_atlas_prepared_response import PreparedFieldBudget, _geometry

Q0 = 0.13882099234357653
ADMITTED_DEGREES = (0, 1, 2, 3)
ADMITTED_NODES = (129, 257, 513)
_CORNER_ABS = 8*np.finfo(float).eps


def _positive_degree(angular_degree):
    if (isinstance(angular_degree, (bool, np.bool_))
            or not isinstance(angular_degree, (int, np.integer))
            or angular_degree not in ADMITTED_DEGREES):
        raise ValueError("angular degree must be an integer 0 through 3")
    return int(angular_degree)


def mode_kernel(radius, angular_degree=0):
    """K_l in chi_UV+K_l chi=0, reusing the certified l=0 coefficient.

    Both the curvature term and the centrifugal term decrease with r on the
    retained radii. The sharp maximum is therefore at the radial floor.
    """
    degree = _positive_degree(angular_degree)
    base = exterior_coefficient(radius)
    if degree == 0:
        return base
    x = np.asarray(radius, dtype=float)
    centrifugal = (8*M_OVER_ELL**3*degree*(degree+1)
                   *np.exp(-x/(2*M_OVER_ELL))/x**3)
    return base+centrifugal


def contraction_bound(angular_degree=0):
    """Inherited sup-norm certificate, using the floor ratio rather than r_max.

    K_l(r_floor)=K_0(r_floor)*(1+r_floor*l(l+1)/(2m)) and r_floor/(2m)=1/2,
    so q_l<=q0*(1+l(l+1)/2). This bounds the integral operator; it is not a
    numerical integration error.
    """
    degree = _positive_degree(angular_degree)
    factor = 1+.5*degree*(degree+1)
    bound = Q0*factor
    if not bound < 1:
        raise ValueError("angular degree leaves the admitted contraction window")
    return {
        "q0": Q0, "floor_factor": factor, "q_l": bound,
        "angular_degree": degree, "numerical_error_enclosure": False,
    }


def kruskal_factor(radius):
    """Dimensionless C in -C dU dV at ell=1."""
    x = np.asarray(radius, dtype=float)
    return 32*M_OVER_ELL**3*np.exp(-x/(2*M_OVER_ELL))/x


def radius_null_derivatives(U, V, radius):
    """r_U=V/F'(r) and r_V=U/F'(r) from UV=F(r)."""
    x = np.asarray(radius, dtype=float)
    prime = -x*np.exp(x/(2*M_OVER_ELL))/(2*M_OVER_ELL)**2
    return V/prime, U/prime


def proper_null_speeds(wall_U_prime, radius):
    """Unit-speed null components fixed by Uw' and C together.

    Udot=sqrt(Uw'/C) and Vdot=1/sqrt(C*Uw'), so C*Udot*Vdot=1 and
    Udot/Vdot=Uw'. An independently interpolated Vdot is not used.
    """
    prime = np.asarray(wall_U_prime, dtype=float)
    factor = kruskal_factor(radius)
    if np.any(prime <= 0) or np.any(factor <= 0):
        raise ValueError("proper wall speed needs an increasing null wall and positive C")
    udot = np.sqrt(prime/factor)
    vdot = 1/np.sqrt(factor*prime)
    return udot, vdot


def outward_normal_derivative(chi_U, chi_V, wall_value, wall_U, V, wall_U_prime):
    """Outward nu+(phi) for phi=chi/r, including the scalar radius terms.

    nu+=Udot dU-Vdot dV. A spatially constant phi has chi=phi*r, so the
    radius terms cancel the chain rule and the derivative is zero.
    """
    radius = np.asarray(radius_from_null(wall_U, V), dtype=float)
    radius_U, radius_V = radius_null_derivatives(wall_U, V, radius)
    udot, vdot = proper_null_speeds(wall_U_prime, radius)
    value = np.asarray(wall_value, dtype=float)
    phi = value/radius
    derivative = (udot*(np.asarray(chi_U, dtype=float)/radius-value*radius_U/radius**2)
                  -vdot*(np.asarray(chi_V, dtype=float)/radius-value*radius_V/radius**2))
    return {
        "r": radius, "phi": phi, "outward_dphi": derivative,
        "Udot": udot, "Vdot": vdot, "r_U": radius_U, "r_V": radius_V,
    }


def _zero_trace(coordinate):
    values = np.zeros(np.shape(coordinate), dtype=float)
    return values, values.copy()


def _as_trace(trace, coordinate, name):
    if trace is None:
        return _zero_trace(coordinate)
    try:
        value, derivative = trace(coordinate)
    except TypeError as error:
        raise ValueError(name+" trace must accept the coordinate and return value, derivative") from error
    value = np.asarray(value, dtype=float)
    derivative = np.asarray(derivative, dtype=float)
    shape = np.shape(coordinate)
    if value.shape != shape or derivative.shape != shape or not np.all(np.isfinite(value)) or not np.all(np.isfinite(derivative)):
        raise ValueError(name+" trace must return finite value and derivative arrays of the coordinate shape")
    return value, derivative


def _close(left, right):
    scale = 1+abs(float(left))+abs(float(right))
    return abs(float(left)-float(right)) <= _CORNER_ABS*scale


def _corner_integral(U0, Ua, Va, bottom, angular_degree, order):
    node, weight = np.polynomial.legendre.leggauss(order)
    half = .5*(Ua-U0)
    points = .5*(Ua+U0)+half*node
    values, _ = _as_trace(bottom, points, "bottom")
    kernel = mode_kernel(radius_from_null(points, Va), angular_degree)
    return float(half*(weight @ (kernel*values)))


class ExteriorDirichletSolution:
    """Wall response of one retarded solve. Derivatives are integral identities.

    ``boundary(V)`` returns the prescribed phi=y/r and nu+(phi). Spatial
    product integration, geometry interpolation and DOP853 error are separate
    from the analytic contraction bound. Higher corner jets are not implied
    by the C0/C1 conditions enforced here.
    """

    def __init__(self, grid, cuts, segments, geometry, domain, trace, angular_degree, budget, bound):
        self.grid = grid
        self.cuts = cuts
        self.segments = segments
        self.geometry = geometry
        self.domain = domain
        self.trace = trace
        self.angular_degree = angular_degree
        self.budget = budget
        self.contraction = bound
        self.numerical_error_enclosure = False

    def _state(self, normalized):
        indices = np.minimum(np.searchsorted(self.cuts, normalized, side="right")-1, len(self.segments)-1)
        state = np.empty((self.grid.size, normalized.size))
        for index in np.unique(indices):
            mask = indices == index
            state[:, mask] = self.segments[index](normalized[mask])
        return state

    def boundary(self, V):
        """Evaluate the wall map at dimensionless V in [Va, Vb]."""
        coordinate = np.atleast_1d(np.asarray(V, dtype=float))
        Va, Vb = self.domain["Va"], self.domain["Vb"]
        if (coordinate.ndim != 1 or not np.all(np.isfinite(coordinate))
                or np.any(coordinate < Va) or np.any(coordinate > Vb)):
            raise ValueError("boundary query exits the finite wall interval")
        self.budget.check(coordinate.size*self.grid.size)
        normalized = (coordinate-Va)/(Vb-Va)
        state = self._state(normalized)
        phi = np.empty(coordinate.size)
        outward = np.empty(coordinate.size)
        chi_U = np.empty(coordinate.size)
        chi_V = np.empty(coordinate.size)
        integral = np.empty(coordinate.size)
        for column, value in enumerate(coordinate):
            row = _wall_row(value, state[:, column], self.grid, self.geometry, self.trace, self.angular_degree)
            normal = outward_normal_derivative(
                row["chi_U"], row["chi_V"], row["y"], row["wall_U"], value, row["wall_U_prime"])
            phi[column] = float(np.asarray(normal["phi"]))
            outward[column] = float(np.asarray(normal["outward_dphi"]))
            chi_U[column] = row["chi_U"]
            chi_V[column] = row["chi_V"]
            integral[column] = row["integral"]
        return {
            "V": coordinate.copy(), "phi": phi, "outward_dphi": outward,
            "chi_U": chi_U, "chi_V": chi_V, "wall_integral": integral,
            "angular_degree": self.angular_degree,
            "analytic_contraction_bound": self.contraction["q_l"],
            "numerical_error_enclosure": False,
            "units": "dimensionless ell=1; phi=chi/r; outward_dphi per proper length",
        }


def _wall_row(V, z, grid, geometry, trace, angular_degree):
    wall = float(geometry["wall"](V))
    prime = float(geometry["wall"](V, 1))
    clipped = np.minimum(grid, wall)
    source_value, source_derivative = _as_trace(trace["source"], V, "source")
    wall_value, wall_derivative = _as_trace(trace["wall"], V, "wall")
    g, gp = float(source_value), float(source_derivative)
    y, yp = float(wall_value), float(wall_derivative)
    chi = np.where(grid < wall, z+g, y)
    kernel = mode_kernel(radius_from_null(clipped, V), angular_degree)
    integrand = kernel*chi
    prefix = np.r_[0., np.cumsum(.5*(integrand[:-1]+integrand[1:])*np.diff(clipped))]
    # The last positive-width sample is the actual wall, carrying y.
    end = int(np.searchsorted(grid, wall, side="left"))
    accumulated = float(prefix[min(end, grid.size-1)])
    return {
        "wall_U": wall, "wall_U_prime": prime, "y": y,
        "integral": accumulated,
        "chi_V": gp-accumulated,
        "chi_U": (yp-gp+accumulated)/prime,
        "prefix": prefix,
    }


def solve_exterior_dirichlet(wall_trace, source=None, bottom=None, angular_degree=0,
                             nodes=129, geometry_knots=129, budget=None):
    """Solve chi_UV+K_l chi=0 with prescribed wall trace y(V).

    Each trace callback accepts the dimensionless coordinate and returns
    ``(value, derivative)``. ``source(V)`` supplies ``(g, g')`` and
    ``bottom(U)`` supplies ``(b, b')``. The defaults are the zero traces,
    not missing data. ``wall_trace(V)`` supplies ``(y, y')``. Derivatives
    must be returned by the callback; stored samples are not differentiated.

    Corner values require ``b(U0)=g(Va)`` and ``y(Va)=b(Ua)``. Exposing the
    endpoint normal also requires the C1 condition
    ``y'(Va)=Uw'(Va)b'(Ua)+g'(Va)-integral K_l b``. Higher jets are a
    separate regularity condition and are not checked.
    """
    degree = _positive_degree(angular_degree)
    bound = contraction_bound(degree)
    if nodes not in ADMITTED_NODES:
        raise ValueError("admitted field grids have 129, 257 or 513 nodes")
    if budget is None:
        budget = PreparedFieldBudget()
    if not isinstance(budget, PreparedFieldBudget):
        raise ValueError("budget must be a PreparedFieldBudget")
    geometry = _geometry(int(geometry_knots))
    Va, Vb = float(geometry["Va"]), float(geometry["Vb"])
    Ua, Ub = float(geometry["Ua"]), float(geometry["Ub"])
    U0 = float(radius_product(RADIAL_FLOOR))/Va
    if not U0 < Ua < Ub < 0 or not Va < Vb < 0:
        raise ValueError("preparation geometry is not an ordered IV wedge")
    domain = {"U0": U0, "Ua": Ua, "Ub": Ub, "Va": Va, "Vb": Vb}
    trace = {"wall": wall_trace, "source": source, "bottom": bottom}
    _require_corners(domain, geometry, trace, degree, budget)
    grid = np.linspace(U0, Ub, nodes)
    entries = np.empty(nodes)
    state = np.empty(nodes)
    for index, coordinate in enumerate(grid):
        if coordinate <= Ua:
            entries[index] = Va
            bottom_value, _ = _as_trace(bottom, coordinate, "bottom")
            source_value, _ = _as_trace(source, Va, "source")
            state[index] = float(bottom_value)-float(source_value)
        else:
            entry = brentq(lambda value: float(geometry["wall"](value))-coordinate, Va, Vb, xtol=1e-14)
            entries[index] = entry
            wall_value, _ = _as_trace(wall_trace, entry, "wall")
            source_value, _ = _as_trace(source, entry, "source")
            state[index] = float(wall_value)-float(source_value)
    dV = Vb-Va
    def rate(normalized, values, active):
        budget.check(nodes)
        current = min(Vb, max(Va, Va+dV*float(normalized)))
        row = _wall_row(current, values, grid, geometry, trace, degree)
        # Future-entry values stay stored and do not evolve. A node that has
        # reached the wall starts at y-g and then follows -I.
        return np.where(active, -dV*row["prefix"], 0.)

    cuts = [0., .25, .75, 1.]
    for entry in entries[(grid > Ua) & (grid < Ub)]:
        cuts.append((entry-Va)/dV)
    cuts = np.unique(np.clip(cuts, 0., 1.))
    segments = []
    for lower, upper in zip(cuts[:-1], cuts[1:]):
        # A node born at the upper cut stays inactive throughout this
        # segment, including its final RK stages. Its prescribed entry
        # value would otherwise acquire an endpoint quadrature kick.
        midpoint = Va+dV*(lower+(upper-lower)/2)
        active = entries < midpoint

        def owned_rate(normalized, values, active=active):
            return rate(normalized, values, active)

        solution = solve_ivp(
            owned_rate, (lower, upper), state, method="DOP853", rtol=2e-10, atol=1e-12,
            first_step=upper-lower, max_step=min(1/64, upper-lower), dense_output=True)
        if not solution.success or solution.t[-1] != upper:
            raise RuntimeError("exterior Dirichlet solve did not reach its segment endpoint")
        state = solution.y[:, -1]
        segments.append(solution.sol)
    return ExteriorDirichletSolution(
        grid, cuts, tuple(segments), geometry, domain, trace, degree, budget, bound)


def _require_corners(domain, geometry, trace, angular_degree, budget):
    U0, Ua, Va = domain["U0"], domain["Ua"], domain["Va"]
    bottom_corner, _ = _as_trace(trace["bottom"], U0, "bottom")
    source_corner, source_derivative = _as_trace(trace["source"], Va, "source")
    wall_corner, wall_derivative = _as_trace(trace["wall"], Va, "wall")
    bottom_at_wall, bottom_derivative = _as_trace(trace["bottom"], Ua, "bottom")
    if not _close(bottom_corner, source_corner):
        raise ValueError("source and bottom traces disagree at U0")
    if not _close(wall_corner, bottom_at_wall):
        raise ValueError("wall and bottom traces disagree at Ua")
    integral_64 = _corner_integral(U0, Ua, Va, trace["bottom"], angular_degree, 64)
    integral_32 = _corner_integral(U0, Ua, Va, trace["bottom"], angular_degree, 32)
    budget.check(96)
    prime = float(geometry["wall"](Va, 1))
    predicted = prime*float(bottom_derivative)+float(source_derivative)-integral_64
    discrepancy = abs(integral_64-integral_32)
    scale = 1+abs(float(wall_derivative))+abs(predicted)
    # The allowance is the quadrature difference plus roundoff, not a physical mismatch.
    if abs(float(wall_derivative)-predicted) > discrepancy+_CORNER_ABS*scale:
        raise ValueError("wall trace fails the C1 corner condition at Va")
