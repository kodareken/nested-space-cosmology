"""Fourier-on-I times the owned axial plateau, same two-function class.

Chebyshev n<=32 is a refused discretization, not the class. This owner
parameterizes the same compact (w,U) pair. It does not evolve a field,
fit source data, or certify a physical gate.
"""
from dataclasses import dataclass
from math import comb, fsum, pi

import numpy as np

from .nsc_compatible_history_geometry import CompatibleIncomingMetric, CompatibleRadiusDirection
from .nsc_ks_spacetime_variation import chart_coordinates
from .nsc_local_incoming_family import plateau_derivatives


ALLOWED_FOURIER_COUNTS = (8, 16)


def fourier_mode_kinds(count):
    """DC, then (cos k, sin k) pairs, leftover highest cosine if count is even."""
    if isinstance(count, bool) or not isinstance(count, (int, np.integer)) or count < 1:
        raise ValueError('positive Fourier mode count required')
    kinds = [('const', 0)]
    harmonic = 1
    while len(kinds) < int(count):
        kinds.append(('cos', harmonic))
        if len(kinds) < int(count):
            kinds.append(('sin', harmonic))
        harmonic += 1
    return tuple(kinds)


def fourier_xi(z, center, inner):
    width = 2.0 * float(inner)
    if width <= 0 or not np.isfinite([z, center, inner, width]).all():
        raise ValueError('finite positive I width required for the Fourier map')
    return (float(z) - float(center) + float(inner)) / width, 1.0 / width


def fourier_wave_derivatives(z, center, inner, kind, harmonic):
    """Bare trigonometric jet through order four on the I-coordinate."""
    xi, dxi = fourier_xi(z, center, inner)
    if kind == 'const':
        if harmonic != 0:
            raise ValueError('constant Fourier mode must have harmonic 0')
        return np.array([1., 0., 0., 0., 0.])
    if kind not in ('cos', 'sin') or harmonic < 1:
        raise ValueError('cos/sin harmonics must be positive integers')
    omega = 2.0 * pi * harmonic
    phase = omega * xi
    cosine, sine = np.cos(phase), np.sin(phase)
    scale = omega * dxi
    if kind == 'cos':
        waves = np.array([cosine, -sine, -cosine, sine, cosine])
    else:
        waves = np.array([sine, cosine, -sine, -cosine, sine])
    return waves * scale ** np.arange(5)


def fourier_series_derivatives(z, center, inner, coefficients):
    raw = np.asarray(coefficients, float)
    if raw.ndim != 1 or len(raw) == 0 or not np.isfinite(raw).all():
        raise ValueError('finite real Fourier coefficients required')
    total = np.zeros(5)
    for value, (kind, harmonic) in zip(raw, fourier_mode_kinds(len(raw)), strict=True):
        total += float(value) * fourier_wave_derivatives(z, center, inner, kind, harmonic)
    return total


@dataclass(frozen=True)
class LocalFourierAxialFunction:
    """One C∞ compact axial function: Fourier-on-I times the owned plateau."""

    coefficients: tuple
    center: float
    inner: float = .03
    outer: float = .06

    def __post_init__(self):
        raw = np.asarray(self.coefficients)
        if (raw.ndim != 1 or len(raw) == 0 or np.iscomplexobj(raw)
                or not np.isfinite(raw).all()):
            raise ValueError('finite real Fourier coefficients required')
        plateau_derivatives(self.center, self.center, self.inner, self.outer)
        fourier_series_derivatives(self.center, self.center, self.inner, raw)
        object.__setattr__(self, 'coefficients', tuple(map(float, raw)))

    def __call__(self, z, derivative_order):
        if isinstance(derivative_order, bool) or derivative_order not in range(5):
            raise ValueError('owned coherent derivative orders zero through four required')
        z = float(z)
        window = plateau_derivatives(z, self.center, self.inner, self.outer)
        if not np.any(window):
            return 0.
        series = fourier_series_derivatives(z, self.center, self.inner, self.coefficients)
        n = derivative_order
        return float(fsum(comb(n, k) * window[k] * series[n - k] for k in range(n + 1)))


@dataclass(frozen=True)
class LocalFourierIncomingFamily:
    """Same declared (w,U) class, Fourier-on-I amplitudes, locked cutoff."""

    coefficients: object
    center: float | None = None
    axial_inner: float = .03
    axial_outer: float = .06
    normal_inner: float = .007
    normal_outer: float = .03

    def __post_init__(self):
        raw = np.asarray(self.coefficients)
        if (raw.ndim != 2 or raw.shape[0] != 2
                or raw.shape[1] not in ALLOWED_FOURIER_COUNTS
                or np.iscomplexobj(raw) or not np.isfinite(raw).all()):
            raise ValueError('two finite real vectors with 8 or 16 Fourier modes required')
        center = chart_coordinates(1.)[1] + .15 if self.center is None else float(self.center)
        plateau_derivatives(center, center, self.axial_inner, self.axial_outer)
        if (not np.isfinite([self.normal_inner, self.normal_outer]).all()
                or not 0 < self.normal_inner < self.normal_outer):
            raise ValueError('positive ordered normal support radii required')
        object.__setattr__(self, 'center', center)
        object.__setattr__(self, 'coefficients', tuple(tuple(map(float, row)) for row in raw))

    @property
    def interval(self):
        return (self.center - self.axial_inner, self.center + self.axial_inner)

    @property
    def support(self):
        return (self.center - self.axial_outer, self.center + self.axial_outer)

    @property
    def functions(self):
        return tuple(LocalFourierAxialFunction(row, self.center, self.axial_inner, self.axial_outer)
                     for row in self.coefficients)

    def metric(self):
        n = len(self.coefficients[0])
        basis = [LocalFourierAxialFunction(tuple(row), self.center, self.axial_inner, self.axial_outer)
                 for row in np.eye(n)]
        zero = LocalFourierAxialFunction((0.,), self.center, self.axial_inner, self.axial_outer)
        directions = tuple(CompatibleRadiusDirection(w, zero, self.normal_inner, self.normal_outer)
                           for w in basis)
        directions += tuple(CompatibleRadiusDirection(zero, U, self.normal_inner, self.normal_outer)
                            for U in basis)
        return CompatibleIncomingMetric(np.asarray(self.coefficients).ravel(), directions)

    def radius_lower_bound(self):
        """Conservative radius bound from |cos|,|sin|<=1 on the cutoff."""
        up = lambda v: float(np.nextafter(float(v), np.inf))
        w, U = (up(fsum(abs(v) for v in row)) for row in self.coefficients)
        s = float(self.normal_outer)
        s3 = up(up(s * s) * s)
        perturbation = up(up(s * w) + up(up(s3 / 6) * U))
        return float(np.nextafter(1. - perturbation, -np.inf))

    def collocation_nodes(self, count=None):
        count = 2 * len(self.coefficients[0]) + 1 if count is None else count
        if isinstance(count, bool) or not isinstance(count, (int, np.integer)) or count < 3:
            raise ValueError('at least three local collocation nodes required')
        return self.center + self.axial_inner * np.cos(np.linspace(np.pi, 0, int(count)))

    def description(self):
        return {
            'class': 'local source-fixed pure-radius compatible histories',
            'representation': 'Fourier-on-I times owned axial plateau',
            'interval': list(self.interval),
            'axial_support': list(self.support),
            'normal_window': [self.normal_inner, self.normal_outer],
            'coefficients': [list(row) for row in self.coefficients],
            'coefficient_order': 'w then U; DC, then cos k, sin k',
            'radius_lower_bound': self.radius_lower_bound(),
            'duration_selected': False,
            'state_supplied_by_family': False,
            'physical_gate': 'OPEN; actual evolved full source and error enclosure required',
        }


def project_callable_on_fourier(callback, center, inner, count, samples=257):
    """Least-squares Fourier-on-I coefficients of a real axial callable."""
    if isinstance(count, bool) or count not in ALLOWED_FOURIER_COUNTS:
        raise ValueError('Fourier projection uses 8 or 16 modes')
    if isinstance(samples, bool) or samples < 2 * count + 1:
        raise ValueError('enough interior samples required for a stable projection')
    z = center + inner * np.cos(np.linspace(np.pi, 0, int(samples)))
    values = np.array([float(callback(point)) for point in z], float)
    kinds = fourier_mode_kinds(count)
    design = np.array([
        [fourier_wave_derivatives(point, center, inner, kind, harmonic)[0]
         for kind, harmonic in kinds]
        for point in z
    ], float)
    coeff, *_ = np.linalg.lstsq(design, values, rcond=None)
    if not np.isfinite(coeff).all():
        raise ValueError('Fourier projection must stay finite')
    return tuple(map(float, coeff))
