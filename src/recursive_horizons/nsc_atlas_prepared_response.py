"""Bounded evaluation of the atlas's prepared conformal-scalar response.

The clipped integral evolution is the v7 Volterra problem with its order of
integration reversed. Constant continuation beyond the wall is numerical
bookkeeping only: the exterior kernel is never applied in the excised cap.
"""
import math
import time
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp, simpson
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from scipy.optimize import brentq

from .nsc_atlas_causal_reference import reference_wall_events
from .nsc_atlas_causal_observers import de_sitter_wall_events
from .nsc_atlas_field_transfer import arrival_from_wall_event
from .nsc_atlas_field_preparation import (
    preparation_record, radius_from_null, exterior_coefficient, source_profile,
    M_OVER_ELL, WALL_LOWER, WALL_UPPER,
)

SCHEMA = "nsc-atlas-prepared-response-v1"
CPU_SECONDS = 60
MAX_KERNEL_VALUES = 10_000_000
MAX_RHS_PER_SOLVE = 6000
MAX_FIELD_NODES = 513
MAX_BULK_CPU_SECONDS = 120
MAX_BULK_KERNEL_VALUES = 30_000_000


class _Budget:
    def __init__(self):
        self.started = time.process_time()
        self.kernel_values = 0

    def check(self, kernel_values=0):
        self.kernel_values += int(kernel_values)
        if (self.kernel_values > MAX_KERNEL_VALUES
                or time.process_time()-self.started > CPU_SECONDS):
            raise RuntimeError("prepared-response finite compute budget exhausted")


class PreparedFieldBudget:
    """Shared, finite budget for a bulk pilot or its complete refinement comparison.

    A caller may reduce these admission ceilings, but cannot enlarge them.
    Construction starts the CPU clock; solves and all later queries share it.
    """

    def __init__(self, *, cpu_seconds=120, kernel_values_max=30_000_000):
        if (isinstance(cpu_seconds, (bool, np.bool_))
                or not isinstance(cpu_seconds, (int, float, np.integer, np.floating))
                or not math.isfinite(cpu_seconds)
                or not 0 < cpu_seconds <= MAX_BULK_CPU_SECONDS):
            raise ValueError("bulk CPU budget must be positive and at most 120 seconds")
        if (isinstance(kernel_values_max, (bool, np.bool_))
                or not isinstance(kernel_values_max, (int, np.integer))
                or not 0 < kernel_values_max <= MAX_BULK_KERNEL_VALUES):
            raise ValueError("bulk kernel budget must be a positive integer at most 30000000")
        self.cpu_seconds = float(cpu_seconds)
        self.kernel_values_max = int(kernel_values_max)
        self.started = time.process_time()
        self.kernel_values = 0

    @property
    def evaluations(self):
        """Same cumulative count for cap/reference callers sharing this budget."""
        return self.kernel_values

    def check(self, kernel_values=0):
        if (isinstance(kernel_values, (bool, np.bool_))
                or not isinstance(kernel_values, (int, np.integer))
                or kernel_values < 0):
            raise ValueError("kernel evaluation count must be a nonnegative integer")
        self.kernel_values += int(kernel_values)
        if (self.kernel_values > self.kernel_values_max
                or time.process_time()-self.started > self.cpu_seconds):
            raise RuntimeError("prepared-field finite compute budget exhausted")


@lru_cache(maxsize=2)
def _geometry(knots):
    if knots not in (129, 257):
        raise ValueError("admitted geometry uses 129 or 257 direct wall knots")
    x = np.linspace(WALL_LOWER, WALL_UPPER, knots)
    wall = reference_wall_events(x, 64)
    cap = de_sitter_wall_events(x, 1., 64)
    arrival, stretch = arrival_from_wall_event(x, cap["t_flat"], cap["gamma"], cap["v"])
    V = wall["V"]
    curves = {
        "Va": float(V[0]), "Vb": float(V[-1]),
        "Ua": float(wall["U"][0]), "Ub": float(wall["U"][-1]),
        "wall": CubicHermiteSpline(V, wall["U"], wall["ell_dU_dtau"]/wall["ell_dV_dtau"], extrapolate=False),
        "arrival": CubicHermiteSpline(V, arrival, stretch/wall["ell_dV_dtau"], extrapolate=False),
        "Vdot": PchipInterpolator(V, wall["ell_dV_dtau"], extrapolate=False),
    }
    for key in ("wall", "arrival"):
        if _minimum_derivative(curves[key]) <= 0:
            raise ValueError("interpolated wall and observer clocks must be monotone")
    return curves


def _minimum_derivative(spline):
    """Check the full derivative quadratic, including each internal vertex."""
    c, widths = spline.c, np.diff(spline.x)
    minima = [c[2], 3*c[0]*widths**2+2*c[1]*widths+c[2]]
    vertex = np.divide(-c[1], 3*c[0], out=np.zeros_like(widths), where=c[0] != 0)
    active = (c[0] != 0) & (vertex > 0) & (vertex < widths)
    minima.append(np.where(active, 3*c[0]*vertex**2+2*c[1]*vertex+c[2], np.inf))
    return float(np.min(minima))


def transverse_kernel_integral(lower, upper, V):
    """Exact integral of K dU; endpoints must lie in the regular IV patch."""
    lower, upper, V = np.broadcast_arrays(lower, upper, V)
    if np.any(lower > upper):
        raise ValueError("transverse integration endpoints are reversed")
    a = radius_from_null(lower, V)
    b = radius_from_null(upper, V)
    return 2*M_OVER_ELL**2/V*(b**-2-a**-2)


def _solve(nodes, geometry, prepared, budget, *, tolerance=2e-10, include_bulk=False,
           include_dense=False):
    """Product integration in U and borrowed DOP853 integration in V.

    Integrate the correction z=chi-g, so exact source cancellation does not
    rely on numerically integrating g'. The terminal U node always stores the
    wall value because the solution is constant outside the moving wall.
    Existing callers receive two entries, or three with include_bulk. The
    explicit include_dense option adds a fourth entry for retained-wedge
    evaluation of this same solve, without another field evolution.
    """
    if nodes not in (129, 257, 513):
        raise ValueError("admitted field grids have 129, 257 or 513 nodes")
    domain = prepared["geometry"]
    U0, Ub = domain["U0"], domain["Ub"]
    Va, Vb = domain["Va"], domain["Vb"]
    dV = Vb-Va
    U = np.linspace(U0, Ub, nodes)
    calls = 0

    def rate(t, z, *, counted=True):
        nonlocal calls
        if counted:
            calls += 1
            if calls > MAX_RHS_PER_SOLVE:
                raise RuntimeError("prepared-response RHS budget exhausted")
        V = min(Vb, max(Va, Va+dV*float(t)))
        wall = min(Ub, max(domain["Ua"], float(geometry["wall"](V))))
        clipped = np.minimum(U, wall)
        budget.check(nodes)
        K = exterior_coefficient(radius_from_null(clipped, V))
        g = float(source_profile(V, Va, Vb)[0])
        # The partial cell ends at the actual wall, whose value is z[-1].
        field = np.where(U < wall, z, z[-1])+g
        integrand = K*field
        prefix = np.r_[0., np.cumsum((integrand[:-1]+integrand[1:])*.5*np.diff(clipped))]
        return -dV*prefix

    # Resolve changes of the spatial product-integration cell explicitly.
    # An adaptive high-order solver otherwise repeatedly locates these C1
    # transitions itself, exhausting calls without improving the spatial model.
    cuts = [0., .25, .75, 1.]
    for u in U[(U > domain["Ua"]) & (U < Ub)]:
        root = brentq(lambda v: float(geometry["wall"](v))-u, Va, Vb, xtol=1e-14)
        cuts.append((root-Va)/dV)
    cuts = np.unique(cuts)
    state = np.zeros(nodes)
    segments = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        solution = solve_ivp(rate, (a, b), state, method="DOP853", rtol=tolerance,
                             atol=1e-12, first_step=b-a, max_step=min(1/64, b-a),
                             dense_output=True)
        if not solution.success or solution.t[-1] != b:
            raise RuntimeError("prepared response did not reach its admitted endpoint")
        state = solution.y[:, -1]
        segments.append(solution.sol)

    def extract(t):
        t = np.atleast_1d(np.asarray(t, dtype=float))
        if np.any(t < 0) or np.any(t > 1) or not np.all(np.isfinite(t)):
            raise ValueError("response query exits the finite preparation window")
        V = np.clip(Va+dV*t, Va, Vb)
        indices = np.minimum(np.searchsorted(cuts, t, side="right")-1, len(segments)-1)
        correction = np.column_stack([segments[index](value) for index, value in zip(indices, t)])
        g, gp = source_profile(V, Va, Vb)
        zt = np.array([rate(v, correction[:, j], counted=False)[-1] for j, v in enumerate(t)])
        J = gp+zt/dV
        T_V = geometry["arrival"](V, 1)
        if np.any(T_V <= 0):
            raise ValueError("observer clock left its monotone domain")
        return {
            "source_coordinate": t, "V": V,
            "source_chi_over_epsilon": g,
            "wall_chi_over_epsilon": g+correction[-1],
            "wall_dchi_dV_over_epsilon": J,
            "ell_I_over_epsilon": 2*geometry["Vdot"](V)*J,
            "arrival_t_C_over_ell": geometry["arrival"](V),
            "ell_phi_C_over_epsilon": 2*J/T_V,
            "zero_kernel_phi_C_over_epsilon": 2*gp/T_V,
        }

    if include_dense:
        return (extract, calls, _bulk_evaluator(U, cuts, segments, geometry, prepared, budget),
                {"grid": U, "cuts": cuts, "segments": tuple(segments), "prepared": prepared})
    if include_bulk:
        return extract, calls, _bulk_evaluator(U, cuts, segments, geometry, prepared, budget)
    return extract, calls


def _bulk_evaluator(grid, cuts, segments, geometry, prepared, budget):
    """Evaluate the stored dense solution only on the retained rectangle.

    The field uses piecewise linear spatial interpolation. Its V derivative
    integrates the piecewise linear K*chi integrand, including the partial
    terminal cell. Its U derivative uses the continuum integral identity
    with that same numerical field, rather than differentiating samples.
    """
    domain = prepared["geometry"]
    U0, Ua, Ub = (domain[key] for key in ("U0", "Ua", "Ub"))
    Va, Vb = (domain[key] for key in ("Va", "Vb"))
    dV = Vb-Va
    gauss = {order: np.polynomial.legendre.leggauss(order) for order in (32, 64)}

    def field_values(points, times):
        budget.check()
        indices = np.minimum(np.searchsorted(cuts, times, side="right")-1, len(segments)-1)
        correction = np.empty((len(grid), len(times)))
        for index in np.unique(indices):
            mask = indices == index
            correction[:, mask] = segments[index](times[mask])
        values = np.empty((len(points), len(times)))
        for column, t in enumerate(times):
            V = min(Vb, max(Va, Va+dV*float(t)))
            wall = min(Ub, max(Ua, float(geometry["wall"](V))))
            g = float(source_profile(V, Va, Vb)[0])
            clipped = np.minimum(grid, wall)
            chi = np.where(grid < wall, correction[:, column], correction[-1, column])+g
            values[:, column] = np.interp(points, clipped, chi)
        return values

    def bulk(points, normalized_V, *, order=64):
        points = np.asarray(points, dtype=float)
        if (points.ndim != 1 or not 1 <= points.size <= MAX_FIELD_NODES
                or not np.all(np.isfinite(points))
                or np.any(points < U0) or np.any(points > Ua)):
            raise ValueError("bulk query needs 1 to 513 finite U values in the retained rectangle [U0, Ua]")
        if (isinstance(normalized_V, (bool, np.bool_))
                or not isinstance(normalized_V, (int, float, np.integer, np.floating))
                or not math.isfinite(normalized_V) or not 0 <= normalized_V <= 1):
            raise ValueError("bulk query exits the finite preparation window")
        if (isinstance(order, (bool, np.bool_))
                or not isinstance(order, (int, np.integer)) or order not in (32, 64)):
            raise ValueError("bulk derivative quadrature uses Gauss order 32 or 64")
        budget.check()
        t = float(normalized_V)
        V = min(Vb, max(Va, Va+dV*t))
        wall = min(Ub, max(Ua, float(geometry["wall"](V))))
        clipped = np.minimum(grid, wall)
        # The interpolation is at most the current wall; there is no exterior
        # kernel evaluation or physical continuation on the cap side.
        chi_grid = field_values(clipped, np.array([t]))[:, 0]
        chi = np.interp(points, clipped, chi_grid)
        budget.check(len(grid))
        integrand = exterior_coefficient(radius_from_null(clipped, V))*chi_grid
        prefix = np.r_[0., np.cumsum(.5*(integrand[:-1]+integrand[1:])*np.diff(clipped))]
        indices = np.minimum(np.searchsorted(clipped, points, side="right")-1, len(grid)-2)
        partial = .5*(integrand[indices]+np.interp(points, clipped, integrand))*(points-clipped[indices])
        chi_V = float(source_profile(V, Va, Vb)[1])-prefix[indices]-partial

        # U<=Ua is retained for the entire past interval, so no inverse-wall
        # lower limit or moving-boundary derivative enters this identity.
        chi_U = np.zeros_like(points)
        node, weight = gauss[order]
        time_cuts = sorted({0., t, *(edge for edge in (.25, .75) if edge < t)})
        for lower, upper in zip(time_cuts[:-1], time_cuts[1:]):
            half = .5*(upper-lower)
            times = .5*(upper+lower)+half*node
            values = field_values(points, times)
            null_V = Va+dV*times
            budget.check(points.size*order)
            K = exterior_coefficient(radius_from_null(points[:, None], null_V[None, :]))
            chi_U -= dV*half*((K*values) @ weight)

        x = radius_from_null(points, V)
        inverse_Fprime = -4*M_OVER_ELL**2*np.exp(-x/(2*M_OVER_ELL))/x
        x_U, x_V = V*inverse_Fprime, points*inverse_Fprime
        p = chi/x
        return {"U": points.copy(), "V": np.full_like(points, V),
                "chi": chi, "chi_U": chi_U, "chi_V": chi_V,
                "x": x, "x_U": x_U, "x_V": x_V,
                "p": p, "p_U": (chi_U-p*x_U)/x, "p_V": (chi_V-p*x_V)/x}

    return bulk


class PreparedFieldSolution:
    """Bounded numerical field; no extrapolation or continuum error enclosure.

    ``wall(t_array)`` has the existing response-extractor keys. ``bulk(U, t,
    order=32|64)`` returns chi=chi/epsilon and p=ell*phi/epsilon, their first
    null derivatives, and the analytic reference radius x=R/ell and its
    derivatives. Here t=(V-Va)/(Vb-Va), while all returned derivatives use
    U,V, not t. Spatial, geometry and 32/64 derivative refinements remain
    separate diagnostics to be made by the consuming calculation.
    ``dense_data`` exposes the original grid, dense segment partition and
    preparation to WedgeField; it is numerical state, not a cap continuation.
    """

    def __init__(self, nodes, geometry_knots, geometry, prepared, budget):
        self.nodes = nodes
        self.geometry_knots = geometry_knots
        self.geometry = geometry
        self.domain = {key: prepared["geometry"][key] for key in ("U0", "Ua", "Ub", "Va", "Vb")}
        self.budget = budget
        self._wall, self.rhs_calls, self.bulk, self.dense_data = _solve(
            nodes, geometry, prepared, budget, include_bulk=True, include_dense=True)

    def wall(self, normalized_V):
        values = np.atleast_1d(np.asarray(normalized_V, dtype=float))
        if values.ndim != 1 or not 1 <= values.size <= MAX_FIELD_NODES:
            raise ValueError("wall query needs 1 to 513 normalized V values")
        self.budget.check()
        return self._wall(values)


def solve_prepared_field(nodes=257, geometry_knots=257, budget=None):
    """Solve the fixed preparation once for wall and rectangular bulk queries.

    Supply the same PreparedFieldBudget to every solve in one admitted
    comparison. The old record producer and its 60-second budget are unchanged.
    No saved waveform is interpolated or differentiated by this interface.
    """
    if (isinstance(nodes, (bool, np.bool_)) or not isinstance(nodes, (int, np.integer))
            or nodes not in (129, 257, 513)):
        raise ValueError("admitted field grids have 129, 257 or 513 nodes")
    if (isinstance(geometry_knots, (bool, np.bool_))
            or not isinstance(geometry_knots, (int, np.integer)) or geometry_knots not in (129, 257)):
        raise ValueError("admitted geometry uses 129 or 257 direct wall knots")
    if budget is None:
        budget = PreparedFieldBudget()
    if not isinstance(budget, (PreparedFieldBudget, _Budget)):
        raise ValueError("a shared PreparedFieldBudget is required")
    budget.check()
    prepared = preparation_record()
    geometry = _geometry(int(geometry_knots))
    budget.check()
    return PreparedFieldSolution(int(nodes), int(geometry_knots), geometry, prepared, budget)


def _tail_born(V, geometry, prepared, order, budget):
    """Independent positive B=integral K A g, with the last integral exact."""
    if order not in (32, 64):
        raise ValueError("tail comparison uses fixed Gauss orders 32 and 64")
    d = prepared["geometry"]
    if not d["Va"]+.75*(d["Vb"]-d["Va"]) <= V <= d["Vb"]:
        raise ValueError("analytic tail bound requires post-input V")
    node, weight = np.polynomial.legendre.leggauss(order)
    qlo, qhi = prepared["preparation"]["support_V"]
    q = (qlo+qhi)/2+(qhi-qlo)/2*node
    upper = geometry["wall"](q)
    p = (d["U0"]+upper[:, None])/2+(upper[:, None]-d["U0"])/2*node
    budget.check(order*order)
    K = exterior_coefficient(radius_from_null(p, q[:, None]))
    future = transverse_kernel_integral(p, float(geometry["wall"](V)), V)
    inner = (upper-d["U0"])/2*np.sum(weight*K*future, axis=1)
    g = source_profile(q, d["Va"], d["Vb"])[0]
    return float((qhi-qlo)/2*np.dot(weight, g*inner))


def prepared_response_record(length_scale=1.):
    """Evaluate a fixed preparation with bounded refinement and tail controls."""
    ell = float(length_scale)
    if not math.isfinite(ell) or ell <= 0:
        raise ValueError("length scale must be finite and positive")
    budget = _Budget()
    prepared = preparation_record()
    geometry = _geometry(257)
    t = np.linspace(0., 1., 257)
    outputs, calls = [], []
    for n in (129, 257, 513):
        extract, count = _solve(n, geometry, prepared, budget)
        outputs.append(extract(t))
        calls.append(count)
    final = outputs[-1]
    coarse, middle, fine = [row["ell_phi_C_over_epsilon"] for row in outputs]
    tail = t > .75
    normalized_error = float(np.max(abs(fine-middle)/(1e-6+1e-4*abs(fine))))
    refinement = [float(np.max(abs(middle-coarse))), float(np.max(abs(fine-middle)))]
    tail_difference = float(np.max(abs(outputs[-1]["wall_dchi_dV_over_epsilon"][tail]
                                           -outputs[-2]["wall_dchi_dV_over_epsilon"][tail])))
    qstar = prepared["operator"]["analytic_norm_bound"]
    d = prepared["geometry"]
    G = 64*(d["Vb"]-d["Va"])/315
    Jlower = (1-qstar)*float(exterior_coefficient(.55))**2*(d["Ua"]-d["U0"])**2*G/2
    controls = []
    for s in (.8125, .875, 1.):
        V = d["Va"]+s*(d["Vb"]-d["Va"])
        b32, b64 = [_tail_born(V, geometry, prepared, order, budget) for order in (32, 64)]
        measured = float(extract([s])["wall_dchi_dV_over_epsilon"][0])
        if not (1-qstar)*b64 <= measured <= b64:
            raise ValueError("computed tail violates the independent analytic bracket")
        controls.append({"source_coordinate": s, "V": V, "B_32": b32, "B_64": b64,
                         "J": measured, "lower_from_B": (1-qstar)*b64, "upper_from_B": b64,
                         "B_quadrature_difference": abs(b64-b32)})
    geometry_extract, geometry_calls = _solve(513, _geometry(129), prepared, budget)
    geometry_control = geometry_extract(t)
    geometry_error = float(np.max(abs(geometry_control["ell_phi_C_over_epsilon"]-fine)))
    normalized_geometry_error = float(np.max(abs(geometry_control["ell_phi_C_over_epsilon"]-fine)
                                              /(1e-6+1e-4*abs(fine))))
    geometry_tail_error = float(np.max(abs(geometry_control["wall_dchi_dV_over_epsilon"][tail]
                                              -final["wall_dchi_dV_over_epsilon"][tail])))
    if max(normalized_error, normalized_geometry_error) > 1 or max(tail_difference, geometry_tail_error,
            max(row["B_quadrature_difference"] for row in controls)) >= .1*Jlower:
        raise ValueError("finite waveform or tail refinement unresolved at admitted budget")
    area = float(simpson(fine, x=final["arrival_t_C_over_ell"]))
    exact_area = float(2*(final["wall_chi_over_epsilon"][-1]-final["wall_chi_over_epsilon"][0]))
    budget.check()
    return {
        "schema": SCHEMA, "status": "FINITE_PREPARED_CLASSICAL_RESPONSE_WITH_POSITIVE_TAIL",
        "inputs": {**prepared["inputs"], "length_scale": ell},
        "preparation": prepared["preparation"],
        "domain": {**d, "arrival_window_over_ell": prepared["observer"]["arrival_window_over_ell"]},
        "method": {"equation": "z_t=-(Vb-Va) integral[U0,min(U,Uw(V))] K(z+g)dU; z=chi-g",
                   "spatial": "piecewise trapezoidal product integration, exact moving-wall cutoff",
                   "time": "SciPy DOP853 split at wall/grid crossings and input edges; rtol=2e-10, atol=1e-12, max_step=1/64 in normalized V",
                   "exterior_kernel_evaluated_in_excised_cap": False,
                   "geometry": "Cubic Hermite wall and arrival maps from 257 direct knots",
                   "rigorous_numerical_enclosure": False},
        "samples": {key: value.tolist() for key, value in final.items()},
        "metrics": {"positive_sampled_peak": float(np.max(fine)),
                    "negative_sampled_peak": float(np.min(fine)),
                    "scattering_max_absolute_change": float(np.max(abs(fine-final["zero_kernel_phi_C_over_epsilon"]))),
                    "post_input_min_J": float(np.min(final["wall_dchi_dV_over_epsilon"][tail])),
                    "post_input_final_phi": float(fine[-1]),
                    "signed_observer_area": area, "wall_endpoint_area": exact_area},
        "analytic_tail": {"positive_after_input_in_declared_window": True,
                          "proof": "0<=A²g<=q*Ag implies -Ag<=chi<=-(1-q*)Ag after g=0; J=-integral Kchi>0",
                          "uniform_lower_J_over_epsilon": Jlower,
                          "bracket": "(1-q*)B(V)<=J(V)<=B(V); B=integral KAg>0",
                          "independent_controls": controls},
        "checks": {"field_grid_refinement_max_differences": refinement,
                   "normalized_mixed_refinement_error": normalized_error,
                   "tail_J_refinement_difference": tail_difference,
                   "geometry_129_257_max_phi_difference": geometry_error,
                   "normalized_geometry_refinement_error": normalized_geometry_error,
                   "geometry_129_257_tail_J_difference": geometry_tail_error,
                   "observer_area_identity_difference": abs(area-exact_area),
                   "refinement_is_rigorous_bound": False},
        "budget": {"field_nodes": [129, 257, 513], "geometry_knots": [129, 257],
                   "output_samples": len(t), "rhs_calls": calls+[geometry_calls],
                   "rhs_calls_max_per_solve": MAX_RHS_PER_SOLVE,
                   "kernel_values": budget.kernel_values, "kernel_values_max": MAX_KERNEL_VALUES,
                   "cpu_seconds_max": CPU_SECONDS},
        "limits": ["fixed classical linear probe and independent IV data, not a signal from right parent P",
                   "positive finite-window field tail, not permanent memory, detector record or energy density",
                   "zero-kernel waveform is a diagnostic omission, not a second physical spacetime",
                   "source generation, quantum state, backreaction, formation and later returns remain open",
                   "sampled peaks and numerical refinement, not continuum error enclosures or late-time asymptotics"],
    }
