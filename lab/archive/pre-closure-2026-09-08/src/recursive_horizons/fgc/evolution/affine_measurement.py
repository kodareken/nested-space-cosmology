"""Numerical physical-metric affine-null and Raychaudhuri measurement path.

DEF0 fixes the exact pointwise sign convention.  This module supplies the
independent trajectory-level implementation needed by NUM1 and later COL1:

* interpolate a stored ADM two-jet history;
* integrate one future outgoing physical-metric radial null generator;
* evolve its coordinate tangent with the affine geodesic equation;
* compute the round-sphere expansion from ``dR/dlambda``;
* differentiate that sampled expansion along the trajectory; and
* independently assemble the complete spherical Raychaudhuri right-hand side
  from the local geometry.

The two routes deliberately share only the interpolated spacetime history.
The direct route differentiates the measured expansion samples.  The second
route obtains ``R_ab k^a k^b`` from the warped-product Hessian identity and
then adds the expansion, shear, and twist terms explicitly.  A negative Ricci
contraction alone is never labelled defocusing.

This is measurement infrastructure and validation machinery, not a DEF1
interval or an FGC collapse result.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, isfinite
from numbers import Real
from typing import Mapping, Sequence

import numpy as np

from .floating_jet import FloatJet2, scalar_primal
from .health_monitor import FIELD_ORDER
from .numerical_engine import COMPARATOR_METHOD, METHODS, PRIMARY_METHOD


FIELD_INDEX = {name: index for index, name in enumerate(FIELD_ORDER)}


def _finite(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    answer = float(value)
    if not isfinite(answer):
        raise ValueError(f"{name} must be finite")
    return answer


def _array(name: str, value: object, *, ndim: int) -> np.ndarray:
    answer = np.asarray(value, dtype=np.float64)
    if answer.ndim != ndim or not np.all(np.isfinite(answer)):
        raise ValueError(f"{name} must be a finite {ndim}-dimensional array")
    return np.ascontiguousarray(answer)


@dataclass(frozen=True, slots=True)
class ADMHistory:
    """Uniform time/radius samples of all six ADM field two-jets."""

    times: np.ndarray
    radii: np.ndarray
    u: np.ndarray
    p: np.ndarray
    q: np.ndarray
    a: np.ndarray
    p_r: np.ndarray
    q_r: np.ndarray

    def __post_init__(self) -> None:
        times = _array("times", self.times, ndim=1)
        radii = _array("radii", self.radii, ndim=1)
        if times.size < 2 or radii.size < 5:
            raise ValueError("history requires at least two times and five radii")
        if np.any(np.diff(times) <= 0.0) or np.any(np.diff(radii) <= 0.0):
            raise ValueError("history coordinates must be strictly increasing")
        for name, coordinate in (("times", times), ("radii", radii)):
            spacing = np.diff(coordinate)
            tolerance = 8192.0 * np.finfo(np.float64).eps * max(
                1.0, abs(float(spacing[0]))
            )
            if not np.allclose(spacing, spacing[0], rtol=0.0, atol=tolerance):
                raise ValueError(f"{name} must be uniformly spaced")
        shape = (times.size, radii.size, len(FIELD_ORDER))
        arrays = {}
        for name in ("u", "p", "q", "a", "p_r", "q_r"):
            value = _array(name, getattr(self, name), ndim=3)
            if value.shape != shape:
                raise ValueError(f"{name} must have shape {shape}")
            immutable = value.copy()
            immutable.setflags(write=False)
            arrays[name] = immutable
        for name, value in (("times", times), ("radii", radii)):
            immutable = value.copy()
            immutable.setflags(write=False)
            object.__setattr__(self, name, immutable)
        for name, value in arrays.items():
            object.__setattr__(self, name, value)
        if np.any(self.u[..., FIELD_INDEX["alpha"]] <= 0.0):
            raise ValueError("history lapse must be strictly positive")
        if np.any(self.u[..., FIELD_INDEX["lambda"]] <= 0.0):
            raise ValueError("history radial metric must be strictly positive")

    def _cell(self, coordinate: np.ndarray, value: float) -> tuple[int, float]:
        if value < coordinate[0] or value > coordinate[-1]:
            raise ValueError("trajectory query left the stored history")
        if value == coordinate[-1]:
            return coordinate.size - 2, 1.0
        spacing = coordinate[1] - coordinate[0]
        position = (value - coordinate[0]) / spacing
        lower = min(int(np.floor(position)), coordinate.size - 2)
        return lower, float(position - lower)

    def interpolate(self, name: str, time: Real, radius: Real) -> np.ndarray:
        if name not in {"u", "p", "q", "a", "p_r", "q_r"}:
            raise ValueError("unknown history field")
        t = _finite("time", time)
        r = _finite("radius", radius)
        ti, tx = self._cell(self.times, t)
        ri, rx = self._cell(self.radii, r)
        value = getattr(self, name)
        lower = (1.0 - rx) * value[ti, ri] + rx * value[ti, ri + 1]
        upper = (1.0 - rx) * value[ti + 1, ri] + rx * value[ti + 1, ri + 1]
        answer = (1.0 - tx) * lower + tx * upper
        if answer.shape != (len(FIELD_ORDER),) or not np.all(np.isfinite(answer)):
            raise ValueError("history interpolation produced an invalid field vector")
        return answer

    def local_adm_jets(self, time: Real, radius: Real) -> Mapping[str, FloatJet2]:
        u = self.interpolate("u", time, radius)
        p = self.interpolate("p", time, radius)
        q = self.interpolate("q", time, radius)
        a = self.interpolate("a", time, radius)
        p_r = self.interpolate("p_r", time, radius)
        q_r = self.interpolate("q_r", time, radius)
        return {
            name: FloatJet2(
                u[index],
                dt=p[index],
                dr=q[index],
                dtt=a[index],
                dtr=p_r[index],
                drr=q_r[index],
            )
            for index, name in enumerate(FIELD_ORDER)
        }


@dataclass(frozen=True, slots=True)
class LocalNullGeometry:
    coordinate_speed: float
    lapse: float
    areal_radius: float
    christoffel_base: tuple[tuple[tuple[float, ...], ...], ...]
    radius_first: tuple[float, float]
    radius_hessian: tuple[tuple[float, float], tuple[float, float]]
    metric_base: tuple[tuple[float, float], tuple[float, float]]


def _local_null_geometry(history: ADMHistory, time: float, radius: float) -> LocalNullGeometry:
    jets = history.local_adm_jets(time, radius)
    lapse = jets["alpha"]
    shift = jets["shift"]
    radial = jets["lambda"]
    lam2 = radial * radial
    h_tt = -(lapse * lapse) + lam2 * shift * shift
    h_tr = lam2 * shift
    h_rr = lam2
    metric_array = np.asarray(
        (
            (scalar_primal(h_tt.value), scalar_primal(h_tr.value)),
            (scalar_primal(h_tr.value), scalar_primal(h_rr.value)),
        ),
        dtype=np.float64,
    )
    determinant = float(np.linalg.det(metric_array))
    if not isfinite(determinant) or determinant >= 0.0:
        raise ValueError("affine trajectory left the Lorentzian ADM domain")
    inverse = np.linalg.inv(metric_array)
    metric_derivative = np.asarray(
        (
            (
                (scalar_primal(h_tt.dt), scalar_primal(h_tr.dt)),
                (scalar_primal(h_tr.dt), scalar_primal(h_rr.dt)),
            ),
            (
                (scalar_primal(h_tt.dr), scalar_primal(h_tr.dr)),
                (scalar_primal(h_tr.dr), scalar_primal(h_rr.dr)),
            ),
        ),
        dtype=np.float64,
    )
    gamma_array = np.zeros((2, 2, 2), dtype=np.float64)
    for upper in range(2):
        for first in range(2):
            for second_index in range(2):
                gamma_array[upper, first, second_index] = 0.5 * sum(
                    inverse[upper, source]
                    * (
                        metric_derivative[first, source, second_index]
                        + metric_derivative[second_index, source, first]
                        - metric_derivative[source, first, second_index]
                    )
                    for source in range(2)
                )
    gamma = tuple(
        tuple(tuple(float(gamma_array[a, b, c]) for c in range(2)) for b in range(2))
        for a in range(2)
    )
    metric = tuple(
        tuple(float(metric_array[a, b]) for b in range(2)) for a in range(2)
    )
    radius_jet = jets["R"]
    first = (float(scalar_primal(radius_jet.dt)), float(scalar_primal(radius_jet.dr)))
    second = (
        (float(scalar_primal(radius_jet.dtt)), float(scalar_primal(radius_jet.dtr))),
        (float(scalar_primal(radius_jet.dtr)), float(scalar_primal(radius_jet.drr))),
    )
    hessian = tuple(
        tuple(
            second[a][b]
            - sum(gamma[c][a][b] * first[c] for c in range(2))
            for b in range(2)
        )
        for a in range(2)
    )
    n = float(scalar_primal(lapse.value))
    lam = float(scalar_primal(radial.value))
    beta = float(scalar_primal(shift.value))
    areal = float(scalar_primal(radius_jet.value))
    if n <= 0.0 or lam <= 0.0 or areal <= 0.0:
        raise ValueError("affine trajectory left the positive ADM/areal domain")
    return LocalNullGeometry(
        coordinate_speed=-beta + n / lam,
        lapse=n,
        areal_radius=areal,
        christoffel_base=gamma,
        radius_first=first,
        radius_hessian=hessian,  # type: ignore[arg-type]
        metric_base=metric,  # type: ignore[arg-type]
    )


def _ray_time_rhs(history: ADMHistory, time: float, state: np.ndarray) -> np.ndarray:
    radius, k_t, _affine = state
    if k_t <= 0.0:
        raise ValueError("future affine generator acquired nonpositive k^t")
    geometry = _local_null_geometry(history, time, radius)
    speed = geometry.coordinate_speed
    gamma_t = geometry.christoffel_base[0]
    contraction = gamma_t[0][0] + 2 * gamma_t[0][1] * speed + gamma_t[1][1] * speed**2
    return np.asarray((speed, -k_t * contraction, 1.0 / k_t), dtype=np.float64)


def _ode_step(
    history: ADMHistory,
    method: str,
    time: float,
    state: np.ndarray,
    step: float,
) -> np.ndarray:
    function = lambda t, y: _ray_time_rhs(history, t, y)
    if method == PRIMARY_METHOD:
        k1 = function(time, state)
        k2 = function(time + step / 2, state + step * k1 / 2)
        k3 = function(time + step / 2, state + step * k2 / 2)
        k4 = function(time + step, state + step * k3)
        answer = state + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    elif method == COMPARATOR_METHOD:
        k1 = function(time, state)
        y1 = state + step * k1
        k2 = function(time + step, y1)
        y2 = 0.75 * state + 0.25 * (y1 + step * k2)
        k3 = function(time + step / 2, y2)
        answer = state / 3 + 2 * (y2 + step * k3) / 3
    else:
        raise ValueError("unknown affine integration method")
    if not np.all(np.isfinite(answer)):
        raise ValueError("affine integrator produced a nonfinite state")
    return answer


def _polynomial_derivative(coordinate: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Five-point local polynomial derivative on a nonuniform coordinate."""

    if coordinate.size != values.size or coordinate.size < 7:
        raise ValueError("trajectory derivative requires at least seven matched samples")
    output = np.empty_like(values)
    count = coordinate.size
    for index in range(count):
        start = min(max(index - 2, 0), count - 5)
        selection = slice(start, start + 5)
        offsets = coordinate[selection] - coordinate[index]
        vandermonde = np.vstack([offsets**power for power in range(5)])
        target = np.asarray((0.0, 1.0, 0.0, 0.0, 0.0))
        weights = np.linalg.solve(vandermonde, target)
        output[index] = float(weights @ values[selection])
    return output


@dataclass(frozen=True, slots=True)
class AffineMeasurementResult:
    method: str
    coordinate_times: np.ndarray
    affine_parameters: np.ndarray
    coordinate_radii: np.ndarray
    k_t: np.ndarray
    k_r: np.ndarray
    theta: np.ndarray
    direct_dtheta_dlambda: np.ndarray
    complete_raychaudhuri_rhs: np.ndarray
    ricci_null: np.ndarray
    null_residual: np.ndarray
    affine_residual_t: np.ndarray
    affine_residual_r: np.ndarray
    initial_normalization_residual: float

    def __post_init__(self) -> None:
        names = (
            "coordinate_times",
            "affine_parameters",
            "coordinate_radii",
            "k_t",
            "k_r",
            "theta",
            "direct_dtheta_dlambda",
            "complete_raychaudhuri_rhs",
            "ricci_null",
            "null_residual",
            "affine_residual_t",
            "affine_residual_r",
        )
        arrays = tuple(_array(name, getattr(self, name), ndim=1) for name in names)
        if len({array.size for array in arrays}) != 1 or arrays[0].size < 7:
            raise ValueError("affine measurement arrays must have one common length")
        for name, value in zip(names, arrays, strict=True):
            immutable = value.copy()
            immutable.setflags(write=False)
            object.__setattr__(self, name, immutable)
        object.__setattr__(
            self,
            "initial_normalization_residual",
            _finite("initial_normalization_residual", self.initial_normalization_residual),
        )

    @property
    def interior(self) -> slice:
        return slice(2, -2)

    @property
    def maximum_route_disagreement(self) -> float:
        return float(
            np.max(
                np.abs(
                    self.direct_dtheta_dlambda[self.interior]
                    - self.complete_raychaudhuri_rhs[self.interior]
                ),
                initial=0.0,
            )
        )

    @property
    def maximum_null_residual(self) -> float:
        return float(np.max(np.abs(self.null_residual), initial=0.0))

    @property
    def maximum_affine_residual(self) -> float:
        return float(
            max(
                np.max(np.abs(self.affine_residual_t[self.interior]), initial=0.0),
                np.max(np.abs(self.affine_residual_r[self.interior]), initial=0.0),
            )
        )


def analyze_affine_null(
    history: ADMHistory,
    *,
    method: str,
    launch_time: Real,
    launch_radius: Real,
    final_time: Real,
    step_size: Real,
) -> AffineMeasurementResult:
    """Integrate and compare the two frozen Raychaudhuri routes."""

    if not isinstance(history, ADMHistory):
        raise TypeError("history must be ADMHistory")
    if method not in METHODS:
        raise ValueError("unknown affine method")
    start = _finite("launch_time", launch_time)
    radius0 = _finite("launch_radius", launch_radius)
    end = _finite("final_time", final_time)
    nominal_step = _finite("step_size", step_size)
    if end <= start or nominal_step <= 0.0:
        raise ValueError("affine time interval and step must be positive")
    initial_geometry = _local_null_geometry(history, start, radius0)
    state = np.asarray((radius0, 1.0 / initial_geometry.lapse, 0.0))
    times = [start]
    states = [state.copy()]
    time = start
    while time < end:
        step = min(nominal_step, end - time)
        state = _ode_step(history, method, time, state, step)
        time += step
        # Avoid an accumulated binary64 value infinitesimally outside history.
        if abs(time - end) <= 32 * np.finfo(np.float64).eps * max(1.0, abs(end)):
            time = end
        times.append(time)
        states.append(state.copy())
    coordinate_times = np.asarray(times)
    trajectory = np.asarray(states)
    radii = trajectory[:, 0]
    k_t = trajectory[:, 1]
    affine = trajectory[:, 2]
    if np.any(np.diff(affine) <= 0.0):
        raise ValueError("affine parameter must increase strictly")

    geometries = tuple(
        _local_null_geometry(history, float(t), float(r))
        for t, r in zip(coordinate_times, radii, strict=True)
    )
    speeds = np.asarray([item.coordinate_speed for item in geometries])
    k_r = k_t * speeds
    theta = np.empty_like(k_t)
    rhs = np.empty_like(k_t)
    ricci_null = np.empty_like(k_t)
    null = np.empty_like(k_t)
    gamma_t_kk = np.empty_like(k_t)
    gamma_r_kk = np.empty_like(k_t)
    for index, (geometry, kt, kr) in enumerate(zip(geometries, k_t, k_r, strict=True)):
        k = (float(kt), float(kr))
        radius_prime = sum(k[a] * geometry.radius_first[a] for a in range(2))
        theta[index] = 2.0 * radius_prime / geometry.areal_radius
        hessian_kk = sum(
            k[a] * k[b] * geometry.radius_hessian[a][b]
            for a in range(2)
            for b in range(2)
        )
        ricci_null[index] = -2.0 * hessian_kk / geometry.areal_radius
        # Spherical radial screen shear and hypersurface twist vanish, but are
        # not omitted conceptually: their explicit values are zero here.
        rhs[index] = -0.5 * theta[index] ** 2 - ricci_null[index]
        null[index] = sum(
            geometry.metric_base[a][b] * k[a] * k[b]
            for a in range(2)
            for b in range(2)
        )
        gamma_t_kk[index] = sum(
            geometry.christoffel_base[0][a][b] * k[a] * k[b]
            for a in range(2)
            for b in range(2)
        )
        gamma_r_kk[index] = sum(
            geometry.christoffel_base[1][a][b] * k[a] * k[b]
            for a in range(2)
            for b in range(2)
        )

    direct = _polynomial_derivative(affine, theta)
    dkt = _polynomial_derivative(affine, k_t)
    dkr = _polynomial_derivative(affine, k_r)
    affine_t = dkt + gamma_t_kk
    affine_r = dkr + gamma_r_kk
    normalization = initial_geometry.lapse * k_t[0] - 1.0
    return AffineMeasurementResult(
        method=method,
        coordinate_times=coordinate_times,
        affine_parameters=affine,
        coordinate_radii=radii,
        k_t=k_t,
        k_r=k_r,
        theta=theta,
        direct_dtheta_dlambda=direct,
        complete_raychaudhuri_rhs=rhs,
        ricci_null=ricci_null,
        null_residual=null,
        affine_residual_t=affine_t,
        affine_residual_r=affine_r,
        initial_normalization_residual=normalization,
    )


def polynomial_flrw_history(
    *,
    time_count: int,
    radius_count: int,
    final_time: Real,
    maximum_radius: Real,
) -> ADMHistory:
    """Return the exact ``a(t)=1+t`` flat-FLRW two-jet control history."""

    if time_count < 2 or radius_count < 5:
        raise ValueError("FLRW history counts are too small")
    end = _finite("final_time", final_time)
    outer = _finite("maximum_radius", maximum_radius)
    if end <= 0.0 or outer <= 0.0:
        raise ValueError("FLRW history extents must be positive")
    times = np.linspace(0.0, end, time_count)
    radii = np.linspace(0.0, outer, radius_count)
    shape = (time_count, radius_count, len(FIELD_ORDER))
    arrays = {name: np.zeros(shape) for name in ("u", "p", "q", "a", "p_r", "q_r")}
    for ti, time in enumerate(times):
        scale = 1.0 + time
        arrays["u"][ti, :, FIELD_INDEX["alpha"]] = 1.0
        arrays["u"][ti, :, FIELD_INDEX["lambda"]] = scale
        arrays["u"][ti, :, FIELD_INDEX["R"]] = scale * radii
        arrays["p"][ti, :, FIELD_INDEX["lambda"]] = 1.0
        arrays["p"][ti, :, FIELD_INDEX["R"]] = radii
        arrays["q"][ti, :, FIELD_INDEX["R"]] = scale
        arrays["p_r"][ti, :, FIELD_INDEX["R"]] = 1.0
    return ADMHistory(times=times, radii=radii, **arrays)


__all__ = [
    "ADMHistory",
    "AffineMeasurementResult",
    "FIELD_INDEX",
    "LocalNullGeometry",
    "analyze_affine_null",
    "polynomial_flrw_history",
]
