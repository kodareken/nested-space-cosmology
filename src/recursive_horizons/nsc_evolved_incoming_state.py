"""Source-fixed incoming state from owned prepared-field evolution.

Same quantum state means the same upstream/source preparation, not a covariance
pinned independently of the metric. For an allowed history family,

    C_Sigma[g] = F[g] C_src F[g]†
    dC_Sigma = dF C_src F† + F C_src dF†    with explicit dC_src = 0.

F is the actual rho1 normal/half-density restriction of evolved source columns.
A frozen incoming C0, a legacy matter dictionary, or a tagged caller matrix is
not this state. Source energies are initial labels; they are not substituted as
outgoing momenta. Quadrature weights enter once as repeated sqrt(w/(2 pi)).
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from types import MappingProxyType

import numpy as np
from scipy.linalg import block_diag

from .nsc_compatible_history_geometry import CompatibleIncomingMetric, CompatibleRadiusDirection
from .nsc_ks_spacetime_variation import chart_coordinates
from .nsc_pg_ks_metric_pullback import ks_to_pg, reference_chart
from .nsc_prepared_history_jets import evolve_prepared_field_jets, restriction_covariance_tangent
from .nsc_transmitting_dirac_domain import MODE_TO_CURRENT, TransmittingDiracSeamDomain


_LEGACY_STATE_NAMES = frozenset({
    'C0', 'c0', 'incoming_C0', 'incoming_c0', 'frozen_C0', 'frozen_c0',
    'incoming_covariance', 'frozen_covariance', 'covariance', 'matter',
    'matter_dictionary', 'stress_approximant', 'legacy_matter',
})


def _finite_array(value, dtype, name, *, copy=True):
    if value is None:
        raise ValueError(f'explicit {name} required; missing data are not zero')
    result = np.array(value, dtype=dtype, copy=copy)
    if not np.isfinite(result).all():
        raise ValueError(f'finite {name} required')
    result.setflags(write=False)
    return result


def _finite_real(value, name):
    result = np.asarray(value)
    if np.iscomplexobj(result) and np.any(result.imag != 0):
        raise ValueError(f'real {name} required')
    return _finite_array(result.real, float, name)


def _digest_arrays(*arrays):
    hasher = sha256()
    for value in arrays:
        array = np.ascontiguousarray(value)
        hasher.update(array.dtype.str.encode())
        hasher.update(str(array.shape).encode())
        hasher.update(b'\0')
        hasher.update(array.tobytes(order='C'))
    return hasher.hexdigest()


def rho1_restriction_map(lapse, radial_scale, shift, radius, *,
                         surface_id='fixed-rho1-source-fixed-incoming'):
    """Owned MODE_TO_CURRENT / normal half-density map at a fixed rho1 metric.

    The radius family keeps this intrinsic map constant, so its parameter
    derivative is identically zero. A family that changes N, q, beta or r at
    Sigma is unsupported here; that derivative is not omitted.
    """
    N, q, beta, r = (float(lapse), float(radial_scale), float(shift), float(radius))
    domain = TransmittingDiracSeamDomain(
        lapse=N, radial_scale=q, shift=beta, radius=r, surface_id=surface_id)
    _, inverse, _ = domain.trace_map(np.ones(1))
    return inverse[0] @ MODE_TO_CURRENT / (r * np.sqrt(q))


def rho1_node_index(x):
    index = np.flatnonzero(np.asarray(x) == 1.)
    if len(index) != 1:
        raise ValueError('exact rho1 grid node required; no interpolation')
    return int(index[0])


class CachedCompatibleIncomingMetric:
    """Cache chart and KS axial data; evaluate the actual finite radius each time.

    N, beta and q of this pure-radius family are amplitude-independent. The
    radius is r0 + amplitudes · basis(tau+S) and is not frozen at the zero-
    amplitude background. This is not CachedZeroRadiusProvider.
    """

    def __init__(self, x, metric):
        if not isinstance(metric, CompatibleIncomingMetric):
            raise TypeError(
                'cache wraps CompatibleIncomingMetric; a zero-radius freeze '
                'is not an actual finite-geometry provider')
        self.x = _finite_real(x, 'cached PG spatial grid')
        if self.x.ndim != 1 or len(self.x) < 1:
            raise ValueError('nonempty one-dimensional PG spatial grid required')
        self.original = metric
        self.amplitudes = _finite_real(metric.amplitudes, 'amplitudes')
        origin = chart_coordinates(1.)[0]
        normal = np.empty(len(self.x))
        shift = np.empty(len(self.x))
        N = np.empty(len(self.x))
        beta = np.empty(len(self.x))
        q = np.empty(len(self.x))
        for j, rho in enumerate(self.x):
            T, S = chart_coordinates(float(rho))
            normal[j] = T - origin
            shift[j] = S
            axial = reference_chart(float(rho))[2]
            N[j], beta[j], q[j], _ = ks_to_pg(float(rho), 1., 0., axial, np.sqrt(1. + rho * rho))
        self.normal = _finite_array(normal, float, 'cached normal')
        self.shift = _finite_array(shift, float, 'cached axial shift')
        self.r0 = _finite_array(np.sqrt(1. + self.x * self.x), float, 'background radius')
        self._N = _finite_array(N, float, 'cached lapse')
        self._beta = _finite_array(beta, float, 'cached shift')
        self._q = _finite_array(q, float, 'cached radial scale')

    def _require_grid(self, x):
        if not np.array_equal(np.asarray(x), self.x):
            raise ValueError('cached grid changed')

    def _basis(self, time):
        space = float(time) + self.shift
        if not self.original.directions:
            return np.zeros((0, len(self.x)))
        rows = [np.broadcast_to(np.asarray(direction.value(self.normal, space), float),
                                self.x.shape).copy()
                for direction in self.original.directions]
        return np.stack(rows, axis=0)

    def values(self, time, x):
        self._require_grid(x)
        basis = self._basis(time)
        if len(self.amplitudes):
            radius = self.r0 + self.amplitudes @ basis
        else:
            radius = np.array(self.r0, copy=True)
        if (not np.isfinite(radius).all()) or np.any(radius <= 0):
            raise ValueError('actual radius must stay finite and positive')
        return np.vstack((self._N, self._beta, self._q, radius))

    def log_directions(self, time, x, metric):
        actual = self.values(time, x)
        supplied = np.asarray(metric, float)
        if supplied.shape != actual.shape or np.any(supplied[[0, 2, 3]] <= 0):
            raise ValueError('actual PG metric must have shape (4,n) and positive N,q,r')
        if not np.allclose(supplied, actual, rtol=1e-12, atol=1e-12):
            raise ValueError('supplied metric does not match this actual amplitude family')
        result = np.zeros((len(self.amplitudes), 4, len(self.x)))
        if len(self.amplitudes):
            result[:, 3, :] = self._basis(time) / supplied[3]
        return result

    def verify_against_owned(self, times):
        """Compare cached values/tangents to CompatibleIncomingMetric at every cached node."""
        if not times:
            raise ValueError('several verification times required')
        value_residual = direction_residual = 0.
        radius_change = 0.
        for time in times:
            owned = self.original.values(time, self.x)
            cached = self.values(time, self.x)
            value_residual = max(value_residual, float(np.max(np.abs(owned - cached))))
            owned_dir = self.original.log_directions(time, self.x, owned)
            cached_dir = self.log_directions(time, self.x, cached)
            direction_residual = max(
                direction_residual, float(np.max(np.abs(owned_dir - cached_dir))))
            radius_change = max(radius_change, float(np.max(np.abs(cached[3] - self.r0))))
        if max(value_residual, direction_residual) >= 3e-12:
            raise ArithmeticError('finite-radius cache disagrees with CompatibleIncomingMetric')
        if any(abs(v) > 0 for v in self.amplitudes) and radius_change == 0:
            raise ArithmeticError('cache froze actual radius at the zero-amplitude background')
        return {
            'value_residual': value_residual,
            'direction_residual': direction_residual,
            'maximum_radius_change_from_background': radius_change,
            'times': [float(t) for t in times],
            'grid_points': int(len(self.x)),
            'nonzero_amplitudes': tuple(float(v) for v in self.amplitudes),
        }


class AmplitudeOnlyMetric:
    """The same actual finite metric with zero direction count.

    Used for primal finite-difference runs of the same family. It does not freeze
    the radius at the zero-amplitude background.
    """

    def __init__(self, provider):
        inner = provider.inner if isinstance(provider, AmplitudeOnlyMetric) else provider
        if not callable(getattr(inner, 'values', None)) or not callable(getattr(inner, 'log_directions', None)):
            raise ValueError('metric provider values and explicit log_directions required')
        self.inner = inner

    def values(self, time, x):
        return self.inner.values(time, x)

    def log_directions(self, time, x, metric):
        actual = np.asarray(self.values(time, x), float)
        supplied = np.asarray(metric, float)
        if supplied.shape != actual.shape or not np.allclose(supplied, actual, rtol=1e-12, atol=1e-12):
            raise ValueError('supplied metric does not match this actual amplitude family')
        return np.zeros((0, 4, actual.shape[1]))


def _reject_legacy_kwargs(kwargs):
    if not kwargs:
        return
    names = ', '.join(sorted(kwargs))
    raise ValueError(
        f'unsupported incoming-state argument ({names}); a frozen incoming C0 '
        'or legacy matter dictionary is not an evolved incoming state')


def _as_compatible_family(provider):
    if isinstance(provider, AmplitudeOnlyMetric):
        return _as_compatible_family(provider.inner)
    if isinstance(provider, CachedCompatibleIncomingMetric):
        return provider.original
    if isinstance(provider, CompatibleIncomingMetric):
        return provider
    return None


def _assert_fixed_intrinsic(provider, times, x, rho1_index):
    """Reject families that change the Sigma restriction map; do not omit dR."""
    metric0 = np.asarray(provider.values(float(times[0]), x), float)
    if metric0.shape != (4, len(x)):
        raise ValueError('metric must have shape (4, len(x))')
    reference = metric0[:, rho1_index].copy()
    restriction = rho1_restriction_map(reference[0], reference[2], reference[1], reference[3])
    for index in range(len(times)):
        metric = np.asarray(provider.values(float(times[index]), x), float)
        here = metric[:, rho1_index]
        if float(np.max(np.abs(here - reference))) > 3e-12:
            raise ValueError(
                'unsupported intrinsic-changed family: restriction-map derivative '
                'is required and is not omitted')
        directions = np.asarray(provider.log_directions(float(times[index]), x, metric), float)
        if directions.ndim != 3 or directions.shape[1:] != (4, len(x)):
            raise ValueError('metric log_directions must have shape (ndir, 4, n)')
        if directions.size and float(np.max(np.abs(directions[:, :, rho1_index]))) > 3e-12:
            raise ValueError(
                'unsupported intrinsic-changed family: nonzero restriction-parameter '
                'derivative is not implemented')
        rebuilt = rho1_restriction_map(here[0], here[2], here[1], here[3])
        if float(np.max(np.abs(rebuilt - restriction))) > 3e-12:
            raise ValueError(
                'unsupported intrinsic-changed family: restriction-map derivative '
                'is required and is not omitted')
    return restriction, reference


@dataclass(frozen=True)
class FixedSourcePreparation:
    """Original upstream fibers: occupations, phases, channels and weights stay fixed."""

    covariance: object
    column_weights: object
    energies: object
    digest: str = ''

    def __post_init__(self):
        covariance = _finite_array(self.covariance, complex, 'source covariance')
        weights = _finite_real(self.column_weights, 'column weights')
        energies = _finite_real(self.energies, 'source energy labels')
        if covariance.ndim != 2 or covariance.shape[0] != covariance.shape[1] or covariance.shape[0] < 1:
            raise ValueError('square nonempty source covariance required')
        nsrc = covariance.shape[0]
        if weights.shape != (nsrc,) or energies.shape != (nsrc,):
            raise ValueError('column weights and energy labels must have shape (nsrc,)')
        if np.any(weights <= 0):
            raise ValueError('positive finite source quadrature weights required once')
        digest = _digest_arrays(covariance, weights, energies)
        object.__setattr__(self, 'covariance', covariance)
        object.__setattr__(self, 'column_weights', weights)
        object.__setattr__(self, 'energies', energies)
        object.__setattr__(self, 'digest', digest)

    @classmethod
    def from_signed_blocks(cls, energies, weights, source_blocks):
        """Block-diagonal C_src and repeated sqrt(w/(2 pi)); no angular factor is added."""
        if isinstance(energies, dict) or isinstance(weights, dict) or isinstance(source_blocks, dict):
            raise TypeError(
                'FixedSourcePreparation required; a frozen incoming C0 or matter '
                'dictionary is not an evolved incoming state')
        energies = _finite_real(energies, 'signed source energies')
        weights = _finite_real(weights, 'source quadrature weights')
        blocks = _finite_array(source_blocks, complex, 'source covariance blocks')
        if energies.ndim != 1 or weights.ndim != 1 or energies.shape != weights.shape:
            raise ValueError('one positive weight per signed source energy required')
        if np.any(weights <= 0):
            raise ValueError('positive inherited source weights required')
        if blocks.ndim != 3 or blocks.shape[0] != len(energies) or blocks.shape[1] != blocks.shape[2]:
            raise ValueError('one square source block per signed energy required')
        channels = blocks.shape[1]
        if channels < 1:
            raise ValueError('nonempty source channel inventory required')
        covariance = np.asarray(block_diag(*blocks), complex)
        column_weights = np.repeat(np.sqrt(weights / (2 * np.pi)), channels)
        column_energies = np.repeat(energies, channels)
        return cls(covariance, column_weights, column_energies)


@dataclass(frozen=True)
class HistoryBinding:
    """Identifies the generating metric family so an old matrix cannot be reused."""

    amplitudes: tuple
    inner_radii: tuple
    outer_radii: tuple
    times: object
    spatial_grid: object
    rho1_index: int
    rho1_metric: object
    restriction: object
    restriction_parameter_derivative: float
    S1: float
    source_digest: str
    digest: str = ''

    def __post_init__(self):
        times = _finite_real(self.times, 'history times')
        grid = _finite_real(self.spatial_grid, 'history spatial grid')
        metric = _finite_real(self.rho1_metric, 'rho1 metric')
        restriction = _finite_array(self.restriction, complex, 'restriction map')
        if times.ndim != 1 or len(times) < 2 or np.any(np.diff(times) <= 0):
            raise ValueError('at least two strictly increasing history times required')
        if grid.ndim != 1 or len(grid) < 1:
            raise ValueError('nonempty spatial grid required')
        if metric.shape != (4,) or restriction.shape != (2, 2):
            raise ValueError('rho1 metric (4,) and restriction (2, 2) required')
        if self.restriction_parameter_derivative != 0:
            raise ValueError('fixed intrinsic restriction must have zero parameter derivative')
        if not 0 <= int(self.rho1_index) < len(grid) or grid[int(self.rho1_index)] != 1.:
            raise ValueError('history rho1 index must locate the exact unit-rho node')
        amplitudes = tuple(float(v) for v in self.amplitudes)
        inner = tuple(float(v) for v in self.inner_radii)
        outer = tuple(float(v) for v in self.outer_radii)
        if len(inner) != len(outer) or len(inner) != len(amplitudes):
            raise ValueError('one inner/outer window radius per amplitude required')
        digest = sha256(
            repr((amplitudes, inner, outer, int(self.rho1_index),
                  float(self.S1), self.source_digest)).encode()
            + _digest_arrays(times, grid, metric, restriction).encode()
        ).hexdigest()
        object.__setattr__(self, 'amplitudes', amplitudes)
        object.__setattr__(self, 'inner_radii', inner)
        object.__setattr__(self, 'outer_radii', outer)
        object.__setattr__(self, 'times', times)
        object.__setattr__(self, 'spatial_grid', grid)
        object.__setattr__(self, 'rho1_index', int(self.rho1_index))
        object.__setattr__(self, 'rho1_metric', metric)
        object.__setattr__(self, 'restriction', restriction)
        object.__setattr__(self, 'digest', digest)


def _coverage(source_count, energy_count):
    return MappingProxyType({
        'sampled_source': (
            f'{energy_count} signed source-energy labels, {source_count} columns; '
            'not a complete quadrature'),
        'angular_multiplicity_added': False,
        'quadrature_applied_once': True,
        'continuum_error_bound': None,
        'physical_history_solution': 'OPEN',
        'full_CAR': 'OPEN',
        'cosmology_EXISTENCE': 'OPEN',
        'finite_source_error': 'OPEN',
        'discretization_error': 'OPEN',
        'outgoing_momentum_substitution': None,
    })


@dataclass(frozen=True)
class EvolvedIncomingState:
    """Restricted evolved source columns on Sigma, bound to one history/source.

    Public arrays
    -------------
    z : (nz,) axial incoming coordinate tau + S(1)
    columns : (nz, 2, nsrc) unweighted canonical incoming columns F
    column_tangents : (ndir, nz, 2, nsrc) dF at explicit dC_src = 0
    source_covariance : (nsrc, nsrc) original C_src
    column_weights : (nsrc,) repeated sqrt(w / (2 pi))
    source_energies : (nsrc,) initial labels only

    Weighted columns
    ----------------
    weighted_columns / weighted_column_tangents are properties: columns
    multiplied by column_weights on the source axis. covariance() and
    covariance_tangent() use those weighted columns and retain the full
    cross-z kernel with shape (nz, 2, nz, 2) and (ndir, nz, 2, nz, 2).
    """

    z: object
    columns: object
    column_tangents: object
    source_covariance: object
    column_weights: object
    source_energies: object
    history: HistoryBinding
    coverage: object
    diagnostics: object

    def __post_init__(self):
        z = _finite_real(self.z, 'incoming z')
        columns = _finite_array(self.columns, complex, 'incoming columns')
        tangents = _finite_array(self.column_tangents, complex, 'incoming column tangents')
        covariance = _finite_array(self.source_covariance, complex, 'source covariance')
        weights = _finite_real(self.column_weights, 'column weights')
        energies = _finite_real(self.source_energies, 'source energy labels')
        if not isinstance(self.history, HistoryBinding):
            raise TypeError('HistoryBinding required; a covariance matrix is not a history')
        if z.ndim != 1 or len(z) < 1:
            raise ValueError('nonempty incoming z required')
        if columns.ndim != 3 or columns.shape[0] != len(z) or columns.shape[1] != 2:
            raise ValueError('columns must have shape (nz, 2, nsrc)')
        nsrc = columns.shape[2]
        if nsrc < 1:
            raise ValueError('all source columns must be retained')
        if tangents.ndim != 4 or tangents.shape[1:] != columns.shape:
            raise ValueError('column_tangents must have shape (ndir, nz, 2, nsrc)')
        if covariance.shape != (nsrc, nsrc) or weights.shape != (nsrc,) or energies.shape != (nsrc,):
            raise ValueError('source covariance, weights and energy labels must match nsrc')
        coverage = self.coverage if self.coverage is not None else _coverage(nsrc, nsrc)
        diagnostics = self.diagnostics if self.diagnostics is not None else {}
        object.__setattr__(self, 'z', z)
        object.__setattr__(self, 'columns', columns)
        object.__setattr__(self, 'column_tangents', tangents)
        object.__setattr__(self, 'source_covariance', covariance)
        object.__setattr__(self, 'column_weights', weights)
        object.__setattr__(self, 'source_energies', energies)
        object.__setattr__(self, 'coverage', MappingProxyType(dict(coverage)))
        object.__setattr__(self, 'diagnostics', MappingProxyType(dict(diagnostics)))
        self.validate()

    @property
    def weighted_columns(self):
        """F with source quadrature applied once: shape (nz, 2, nsrc)."""
        return np.asarray(self.columns) * np.asarray(self.column_weights)

    @property
    def weighted_column_tangents(self):
        """dF with the same once-applied weights: shape (ndir, nz, 2, nsrc)."""
        return np.asarray(self.column_tangents) * np.asarray(self.column_weights)

    def validate(self):
        """Finite bound state. Equality of C_Sigma with a reference C0 is not a failure."""
        nsrc = self.columns.shape[-1]
        if self.history.source_digest != _digest_arrays(
                self.source_covariance, self.column_weights, self.source_energies):
            raise ValueError('source binding changed; this state is not a reusable covariance matrix')
        rebuilt = rho1_restriction_map(
            self.history.rho1_metric[0], self.history.rho1_metric[2],
            self.history.rho1_metric[1], self.history.rho1_metric[3])
        if float(np.max(np.abs(rebuilt - self.history.restriction))) > 3e-12:
            raise ValueError('stored restriction does not match the bound rho1 metric')
        if self.history.restriction_parameter_derivative != 0:
            raise ValueError('fixed intrinsic restriction must have zero parameter derivative')
        if np.any(self.column_weights <= 0):
            raise ValueError('positive column weights required')
        hermiticity = float(np.max(np.abs(self.source_covariance - self.source_covariance.conj().T)))
        if hermiticity > 3e-11:
            raise ValueError('source covariance must be Hermitian')
        if self.coverage.get('continuum_error_bound') not in (None, 'OPEN'):
            raise ValueError('continuum error remains OPEN; it is not a boolean proof')
        if self.coverage.get('physical_history_solution') != 'OPEN':
            raise ValueError('physical history solution remains OPEN')
        if self.coverage.get('outgoing_momentum_substitution') not in (None, False, 'not used'):
            raise ValueError('source energies are labels; k=-E is not a column momentum')
        if nsrc != self.source_covariance.shape[0]:
            raise ValueError('source columns were truncated')
        return True

    def covariance(self):
        """C_Sigma = F C_src F† on Sigma, shape (nz, 2, nz, 2), all cross-z entries kept."""
        self.validate()
        nz, _, nsrc = self.columns.shape
        F = self.weighted_columns.reshape(nz * 2, nsrc)
        kernel = F @ self.source_covariance @ F.conj().T
        return kernel.reshape(nz, 2, nz, 2)

    def covariance_tangent(self):
        """dC_Sigma at explicit dC_src=0, shape (ndir, nz, 2, nz, 2)."""
        self.validate()
        nz, _, nsrc = self.columns.shape
        ndir = self.column_tangents.shape[0]
        F = self.weighted_columns.reshape(nz * 2, nsrc)
        dF = self.weighted_column_tangents.reshape(ndir, nz * 2, nsrc)
        dC_src = np.zeros((ndir, nsrc, nsrc), dtype=complex)
        kernel = restriction_covariance_tangent(F, dF, self.source_covariance, dC_src)
        return np.asarray(kernel).reshape(ndir, nz, 2, nz, 2)

    def require_history(self, provider, times, x):
        """A changed geometry cannot reuse this state as a physical covariance."""
        family = _as_compatible_family(provider)
        if family is None:
            raise ValueError(
                'evolved incoming state is bound to its CompatibleIncomingMetric family; '
                'it is not a reusable covariance matrix')
        times = np.asarray(times, float)
        x = np.asarray(x, float)
        if tuple(float(v) for v in family.amplitudes) != self.history.amplitudes:
            raise ValueError(
                'evolved incoming state is bound to its generating metric family; '
                'it is not a reusable covariance matrix')
        if not np.array_equal(times, self.history.times) or not np.array_equal(x, self.history.spatial_grid):
            raise ValueError('evolved incoming state is bound to its generating time/grid')
        inner = tuple(float(d.inner_radius) for d in family.directions)
        outer = tuple(float(d.outer_radius) for d in family.directions)
        if inner != self.history.inner_radii or outer != self.history.outer_radii:
            raise ValueError('evolved incoming state is bound to its compact window')
        return True


def _prepare_provider(owner, provider):
    if isinstance(provider, AmplitudeOnlyMetric):
        return AmplitudeOnlyMetric(_prepare_provider(owner, provider.inner))
    if isinstance(provider, CachedCompatibleIncomingMetric):
        if not np.array_equal(provider.x, owner.x):
            raise ValueError('cached provider grid must match the propagator grid')
        return provider
    if isinstance(provider, CompatibleIncomingMetric):
        return CachedCompatibleIncomingMetric(owner.x, provider)
    if not callable(getattr(provider, 'values', None)) or not callable(getattr(provider, 'log_directions', None)):
        raise ValueError('metric provider values and explicit log_directions required')
    return provider


def evolve_incoming_state(owner, times, provider, initial_fields, initial_tangents,
                          incoming, source, **kwargs):
    """Physical entry: evolve owned prepared columns, then restrict F,dF on Sigma.

    Calls evolve_prepared_field_jets for the supplied metric family. There is no
    incoming C0 argument. Source-fixed means dC_src=0, not dC_Sigma=0.
    """
    _reject_legacy_kwargs(kwargs)
    if not isinstance(source, FixedSourcePreparation):
        if isinstance(source, dict) and (
                'stress_approximant' in source or 'parts' in source or 'incoming_C0' in source):
            raise ValueError(
                'a frozen incoming C0 or legacy matter dictionary is not an evolved incoming state')
        raise TypeError(
            'FixedSourcePreparation required; a frozen incoming C0 or covariance '
            'matrix is not an evolved incoming state')
    times = _finite_real(times, 'ordered caller time coordinates')
    if times.ndim != 1 or len(times) < 2 or np.any(np.diff(times) <= 0):
        raise ValueError('at least two strictly increasing caller time coordinates required')
    provider = _prepare_provider(owner, provider)
    x = _finite_real(owner.x, 'owner spatial grid')
    rho1 = rho1_node_index(x)
    restriction, rho1_metric = _assert_fixed_intrinsic(provider, times, x, rho1)
    family = _as_compatible_family(provider)
    if family is None:
        raise ValueError(
            'CompatibleIncomingMetric family required; unsupported restriction families '
            'are rejected rather than omitting the restriction derivative')
    sample_rows = np.array([rho1, len(x) + rho1], dtype=np.intp)
    run = evolve_prepared_field_jets(
        owner, times, provider, initial_fields, initial_tangents, incoming,
        sample_rows=sample_rows)
    sampled = _finite_array(run['sampled_fields'], complex, 'sampled rho1 fields')
    sampled_tangents = _finite_array(
        run['sampled_tangent_fields'], complex, 'sampled rho1 tangents')
    columns = np.einsum('ab,tbc->tac', restriction, sampled)
    column_tangents = np.einsum('ab,dtbc->dtac', restriction, sampled_tangents)
    S1 = float(chart_coordinates(1.)[1])
    z = np.asarray(times, float) + S1
    history = HistoryBinding(
        amplitudes=tuple(family.amplitudes),
        inner_radii=tuple(d.inner_radius for d in family.directions),
        outer_radii=tuple(d.outer_radius for d in family.directions),
        times=times,
        spatial_grid=x,
        rho1_index=rho1,
        rho1_metric=rho1_metric,
        restriction=restriction,
        restriction_parameter_derivative=0.,
        S1=S1,
        source_digest=source.digest,
    )
    energy_count = len(np.unique(np.round(source.energies, decimals=14)))
    diagnostics = {
        'norm_flux_residual': float(run['norm_flux_residual']),
        'metric_flux_tangent_residual': float(run['metric_flux_tangent_residual']),
        'auxiliary_identity_residual': float(run['auxiliary_identity_residual']),
        'source_covariance_tangent_explicit': 0.,
        'restriction_parameter_derivative': 0.,
        'stationary_phase_reapplied': False,
        'incoming_trace_interpolation': None,
        'outgoing_momentum_k_equals_minus_E': False,
        'physical_preparation_status': run['physical_preparation_status'],
        'physical_history_status': run['physical_history_status'],
        'midpoint_steps': int(run['midpoint_steps']),
        'direction_count': int(run['direction_count']),
        'source_column_count': int(run['source_column_count']),
        'evolved_field': _finite_array(run['evolved_field'], complex, 'evolved field'),
        'tangent_fields': _finite_array(run['tangent_fields'], complex, 'tangent fields'),
    }
    return EvolvedIncomingState(
        z=z,
        columns=columns,
        column_tangents=column_tangents,
        source_covariance=source.covariance,
        column_weights=source.column_weights,
        source_energies=source.energies,
        history=history,
        coverage=_coverage(source.covariance.shape[0], energy_count),
        diagnostics=diagnostics,
    )

from pathlib import Path as _receipt_Path
from .nsc_frozen_sources import selected_source_sha256 as _receipt_source_sha256

# Retained executable control: state_control
# Folded from scripts/derive_nsc_evolved_incoming_state.py at 329ef809; original receipts retain their producer identities.
import os as _state_control_os
import argparse as _state_control_argparse
from hashlib import sha256 as _state_control_sha256
import json as _state_control_json
from pathlib import Path as _state_control_Path
import signal as _state_control_signal
import sys as _state_control_sys
import time as _state_control_time
import numpy as _state_control_np
from recursive_horizons.nsc_posix import ROOT as _state_control_ROOT
_state_control_sys.path.insert(0, str(_state_control_ROOT / 'src'))
_state_control_sys.path.insert(0, str(_state_control_ROOT / 'scripts'))
import derive_nsc_retarded_compatible_response as _state_control_RETARDED
from recursive_horizons.nsc_compatible_history_geometry import CompatibleRadiusDirection as _state_control_CompatibleRadiusDirection, CompatibleIncomingMetric as _state_control_CompatibleIncomingMetric
from recursive_horizons.nsc_evolved_incoming_state import AmplitudeOnlyMetric as _state_control_AmplitudeOnlyMetric, CachedCompatibleIncomingMetric as _state_control_CachedCompatibleIncomingMetric, FixedSourcePreparation as _state_control_FixedSourcePreparation, evolve_incoming_state as _state_control_evolve_incoming_state, restriction_covariance_tangent as _state_control_restriction_covariance_tangent, rho1_restriction_map as _state_control_rho1_restriction_map, rho1_node_index as _state_control_rho1_node_index
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes as _state_control_deterministic_npz_bytes
from recursive_horizons.nsc_transmitting_history_jets import FourthOrderModePropagator as _state_control_FourthOrderModePropagator
import recursive_horizons.nsc_evolved_incoming_state as _state_control_STATE
_state_control_OUTPUT = 'results/development/nsc-evolved-incoming-state.json'
_state_control_ARTIFACT = 'results/development/artifacts/nsc-evolved-incoming-state.npz'
_state_control_OWNED = ('src/recursive_horizons/nsc_evolved_incoming_state.py', 'scripts/derive_nsc_evolved_incoming_state.py', 'tests/test_nsc_evolved_incoming_state.py', 'docs/nsc-evolved-incoming-state.md')
_state_control_RECORDS = ('results/development/nsc-retarded-compatible-response.json', 'results/development/nsc-transmitting-history-modes.json', 'results/development/nsc-pg-ctp-mode-jets.json', 'results/development/nsc-compatible-prepared-history.json')
_state_control_OWNERS = tuple(('src/recursive_horizons/' + name + '.py' for name in ('nsc_prepared_history_jets', 'nsc_compatible_history_geometry', 'nsc_transmitting_history_jets', 'nsc_transmitting_history_modes', 'nsc_transmitting_dirac_domain', 'nsc_ks_spacetime_variation', 'nsc_pg_ks_metric_pullback', 'nsc_transmitting_resolvent', 'nsc_lorentzian', 'nsc_mode_resolved_cauchy_state'))) + ('scripts/derive_nsc_retarded_compatible_response.py', 'docs/nsc-retarded-compatible-preparation.md', 'docs/nsc-preparation-domain-decision.md')
_state_control_AMPLITUDE = 0.001
_state_control_STEP = 0.0001
_state_control_POINTS, _state_control_STEPS = (801, 64)
_state_control_CPU_CAP = 120.0
_state_control_FIELD_TOL = 3e-08
_state_control_ALGEBRA_TOL = 3e-11
_state_control_MAXABS = lambda x: float(_state_control_np.max(_state_control_np.abs(x)))

def _state_control_digest(path):
    candidate = _receipt_Path(path)
    relative = candidate.relative_to(_state_control_ROOT).as_posix() if candidate.is_absolute() else candidate.as_posix()
    return _receipt_source_sha256(_state_control_ROOT, relative)

def _state_control_signature():
    return {p: _state_control_digest(p) for p in (*_state_control_OWNED, *_state_control_RECORDS, *_state_control_OWNERS)}

def _state_control_load_retarded_inputs():
    x, phi, data, records = _state_control_RETARDED.inputs()
    guard = _state_control_RETARDED.support_guard()
    retarded = _state_control_json.loads((_state_control_ROOT / _state_control_RECORDS[0]).read_text())
    if _state_control_digest(retarded['payload']['path']) != retarded['payload']['sha256']:
        raise ValueError('retarded response artifact changed')
    return (x, phi, data, records, guard, retarded)

def _state_control_family(amplitude, guard):
    direction = _state_control_CompatibleRadiusDirection(_state_control_RETARDED.axial_bump(guard['axial_center']), lambda z, n: 0.0, 0.007, 0.03)
    return _state_control_CompatibleIncomingMetric((float(amplitude),), (direction,))

def _state_control_incident(owner, phi, data, ndir):
    n = len(owner.x)
    energies = _state_control_np.repeat(data['energies'], 3)
    tangents = _state_control_np.zeros((ndir, *phi.shape), complex)
    right = _state_control_np.array([n - 1, 2 * n - 1])
    N, beta, q, _ = owner.reference
    speed = _state_control_np.array([N[-1] / q[-1] - beta[-1], -N[-1] / q[-1] - beta[-1]])
    sat = -speed / owner.weights[-1]

    def incoming(_state_control_time):
        forcing = _state_control_np.zeros_like(phi)
        forcing[right] = sat[:, None] * phi[right] * _state_control_np.exp(-1j * energies * _state_control_time)[None, :]
        return (forcing, _state_control_np.zeros_like(tangents))
    return (phi, tangents, incoming)

def _state_control_kernel(columns, weights, covariance):
    nsrc = covariance.shape[0]
    F = (_state_control_np.asarray(columns) * _state_control_np.asarray(weights)).reshape(-1, nsrc)
    return F @ covariance @ F.conj().T

def _state_control_probe_cache(guard):
    lower, upper = guard['radial_support']
    x = _state_control_np.array([-2.0, lower, 0.99, 1.0, 1.01, upper, 1.2])
    metric = _state_control_family(_state_control_AMPLITUDE, guard)
    cache = _state_control_CachedCompatibleIncomingMetric(x, metric)
    times = (0.0, 0.11, 0.15, 0.21, 0.3)
    result = cache.verify_against_owned(times)
    if result['maximum_radius_change_from_background'] <= 0:
        raise ArithmeticError('nonzero-amplitude cache did not change the actual radius')
    zero = _state_control_CompatibleIncomingMetric((0.0,), metric.directions)
    if _state_control_MAXABS(cache.values(0.15, x)[3] - zero.values(0.15, x)[3]) <= 0:
        raise ArithmeticError('cache agrees with a frozen zero-radius background at nonzero amplitude')
    return result

def _state_control_retarded_restriction_replay(retarded, data):
    path = _state_control_ROOT / retarded['payload']['path']
    with _state_control_np.load(path, allow_pickle=False) as arrays:
        field = arrays['801_64/field']
        tangent = arrays['801_64/field_tangent']
        trace = arrays['801_64/trace']
        trace_tangent = arrays['801_64/trace_tangent']
        sample_rows = arrays['801_64/sample_rows']
        x = arrays['801_64/x']
        saved_kernel = arrays['801_64/covariance_tangent']
        source = arrays['source/source']
        weights = arrays['source/weights']
    owner = _state_control_FourthOrderModePropagator(x, data['mass'].item(), data['angular'].item())
    index = _state_control_rho1_node_index(x)
    N, beta, q, radius = owner.reference[:, index]
    restriction = _state_control_rho1_restriction_map(N, q, beta, radius)
    field_trace = restriction @ field[sample_rows]
    tangent_trace = restriction @ tangent[sample_rows]
    root = _state_control_np.repeat(_state_control_np.sqrt(weights / (2 * _state_control_np.pi)), 3)
    C = _state_control_FixedSourcePreparation.from_signed_blocks(data['energies'], weights, source).covariance
    F = (trace * root).reshape(-1, 12)
    dF = (trace_tangent * root).reshape(-1, 12)
    rebuilt = _state_control_restriction_covariance_tangent(F, dF, C, _state_control_np.zeros_like(C))
    return {'final_trace_residual': _state_control_MAXABS(field_trace - trace[-1]), 'final_trace_tangent_residual': _state_control_MAXABS(tangent_trace - trace_tangent[-1]), 'sampled_kernel_residual': _state_control_MAXABS(rebuilt - saved_kernel), 'restriction_parameter_derivative': 0.0, 'old_solver_rerun': False}

def _state_control_one_state(amplitude, ndir, owner, phi, data, times, source, guard):
    metric = _state_control_family(amplitude, guard)
    provider = _state_control_CachedCompatibleIncomingMetric(owner.x, metric)
    if ndir == 0:
        provider = _state_control_AmplitudeOnlyMetric(provider)
    fields, tangents, incoming = _state_control_incident(owner, phi, data, ndir)
    return _state_control_evolve_incoming_state(owner, times, provider, fields, tangents, incoming, source)

class _state_control_CPUExceeded(RuntimeError):
    pass

def _state_control_analyze(arrays, meta):
    weights = arrays['source/column_weights']
    C = arrays['source/covariance']
    field_fd = (arrays['plus/field'] - arrays['minus/field']) / (2 * _state_control_STEP)
    trace_fd = (arrays['plus/columns'] - arrays['minus/columns']) / (2 * _state_control_STEP)
    cov_fd = (_state_control_kernel(arrays['plus/columns'], weights, C) - _state_control_kernel(arrays['minus/columns'], weights, C)) / (2 * _state_control_STEP)
    analytic_cov = arrays['center/covariance_tangent'].reshape(arrays['center/columns'].shape[0] * 2, arrays['center/columns'].shape[0] * 2)
    residuals = {'complete_field_derivative': _state_control_MAXABS(field_fd - arrays['center/field_tangent']), 'incoming_trace_derivative': _state_control_MAXABS(trace_fd - arrays['center/column_tangents']), 'complete_coherent_covariance_derivative': _state_control_MAXABS(cov_fd - analytic_cov), 'norm_flux_residual': float(meta['algebra']['norm_flux_residual']), 'metric_flux_tangent_residual': float(meta['algebra']['metric_flux_tangent_residual']), 'auxiliary_identity_residual': float(meta['algebra']['auxiliary_identity_residual']), 'source_covariance_hermiticity': _state_control_MAXABS(C - C.conj().T), 'retarded_final_trace_residual': float(meta['retarded_restriction']['final_trace_residual']), 'retarded_kernel_residual': float(meta['retarded_restriction']['sampled_kernel_residual']), 'cache_value_residual': float(meta['cache']['value_residual']), 'cache_direction_residual': float(meta['cache']['direction_residual'])}
    F = (arrays['center/columns'] * weights).reshape(-1, C.shape[0])
    dF = (arrays['center/column_tangents'][None, ...] * weights).reshape(1, F.shape[0], C.shape[0])
    frozen = _state_control_restriction_covariance_tangent(F, _state_control_np.zeros_like(dF), C, _state_control_np.zeros((1, *C.shape), complex))
    fake = _state_control_np.zeros((1, *C.shape), complex)
    fake[0, 0, 1] = 0.05
    fake[0, 1, 0] = 0.05
    wrong_source = _state_control_restriction_covariance_tangent(F, dF, C, fake)
    omitted = {'frozen_incoming_columns': _state_control_MAXABS(frozen[0] - analytic_cov), 'nonzero_source_covariance_tangent': _state_control_MAXABS(wrong_source[0] - analytic_cov), 'finite_difference_versus_frozen_columns': _state_control_MAXABS(cov_fd - frozen[0])}
    algebra = max((residuals[k] for k in ('norm_flux_residual', 'metric_flux_tangent_residual', 'auxiliary_identity_residual')))
    passed = residuals['complete_field_derivative'] < _state_control_FIELD_TOL and residuals['incoming_trace_derivative'] < _state_control_FIELD_TOL and (residuals['complete_coherent_covariance_derivative'] < _state_control_FIELD_TOL) and (algebra < _state_control_ALGEBRA_TOL) and (residuals['source_covariance_hermiticity'] < _state_control_ALGEBRA_TOL) and (residuals['retarded_final_trace_residual'] < 3e-12) and (residuals['retarded_kernel_residual'] < 3e-12) and (residuals['cache_value_residual'] < 3e-12) and (min(omitted.values()) > 1e-16) and (meta['field_runs'] == 3) and (not meta['cpu_budget_exceeded'])
    return (residuals, omitted, passed, algebra)

def _state_control_record_from_arrays(arrays, meta, artifact_hash, artifact_bytes):
    residuals, omitted, passed, algebra = _state_control_analyze(arrays, meta)
    current = meta['verification_signature']
    return {'schema': 'NSC-EVOLVED-INCOMING-STATE-v1', 'accountable_author': 'Douglas Ek', 'status': 'PASS: source-fixed adapter control; sampled source and continuum error OPEN' if passed else 'OPEN: source-fixed adapter control residual or coverage gate not met', 'control': {'family': '14_1', 'amplitude': _state_control_AMPLITUDE, 'finite_difference_step': _state_control_STEP, 'grid_points': _state_control_POINTS, 'midpoint_steps': _state_control_STEPS, 'PG_time_window': [0.0, 0.3], 'selected_physical_duration': None, 'numerical_normal_window_radii': [0.007, 0.03], 'field_runs': meta['field_runs'], 'field_run_cap': 3, 'initial_and_incident_tangents': 'explicit 0 by unchanged entire past/upstream support', 'source_law': 'unchanged original coherent source blocks; explicit dC_src=0', 'weights': 'sqrt(w/(2*pi)) repeated on each source column; no multiplicity factor added', 'restriction': 'actual rho1 MODE_TO_CURRENT/normal half-density; dR=0 for this radius family', 'stationary_phase_reapplied': False, 'outgoing_momentum_substitution': None, 'field_derivative_tolerance': _state_control_FIELD_TOL, 'covariance_derivative_tolerance': _state_control_FIELD_TOL, 'flux_algebra_tolerance': _state_control_ALGEBRA_TOL, 'support': meta['support']}, 'residuals': residuals, 'omitted_term_effects': omitted, 'algebra_maximum': algebra, 'cache_verification': meta['cache'], 'retarded_array_only_restriction': meta['retarded_restriction'], 'runtime': {k: meta[k] for k in ('CPU_seconds', 'wall_seconds', 'cpu_budget_seconds', 'cpu_budget_exceeded', 'case_CPU_seconds', 'field_runs')}, 'scope': {'source_fixed_incoming_state': True, 'frozen_incoming_C0': False, 'sampled_source': 'OPEN', 'continuum_error_bound': None, 'finite_source_error': 'OPEN', 'discretization_error': 'OPEN', 'full_CAR': 'OPEN', 'physical_history_solution': 'OPEN', 'cosmology_EXISTENCE': 'OPEN', 'metric_timestep': False, 'new_stress_or_Gamma': False, 'A_q_Omega_zeta_Vfull_refit': False, 'old_baseline_scattering_horizon_or_source_grids_regenerated': False, 'publication': False}, 'source_hashes': {p: current[p] for p in _state_control_OWNED}, 'input_hashes': {p: current[p] for p in (*_state_control_RECORDS, *_state_control_OWNERS)}, 'reused_artifacts': meta['reused_artifacts'], 'payload': {'path': _state_control_ARTIFACT, 'sha256': artifact_hash, 'bytes': artifact_bytes}, 'reproducer': 'python3 scripts/derive_nsc_evolved_incoming_state.py --check'}

def _state_control_publish(arrays, meta):
    current = _state_control_signature()
    meta['verification_signature'] = dict(current)
    meta['signature'] = current
    payload = {**arrays, 'metadata_json': _state_control_np.frombuffer(_state_control_json.dumps({k: meta[k] for k in meta if k != 'signature'}, sort_keys=True).encode(), _state_control_np.uint8)}
    raw = _state_control_deterministic_npz_bytes(payload)
    path = _state_control_ROOT / _state_control_ARTIFACT
    path.write_bytes(raw)
    record = _state_control_record_from_arrays(payload, meta, _state_control_sha256(raw).hexdigest(), len(raw))
    (_state_control_ROOT / _state_control_OUTPUT).write_text(_state_control_json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + '\n')
    return record

def _state_control_run():
    start_cpu = _state_control_time.process_time()
    start_wall = _state_control_time.monotonic()
    before = _state_control_signature()
    meta = {'completed': [], 'case_CPU_seconds': {}, 'cpu_budget_seconds': _state_control_CPU_CAP, 'cpu_budget_exceeded': False, 'field_runs': 0}

    def stop(signum, frame):
        raise _state_control_CPUExceeded('120 second producer CPU budget exhausted')
    old_handler = _state_control_signal.signal(_state_control_signal.SIGPROF, stop)
    _state_control_signal.setitimer(_state_control_signal.ITIMER_PROF, _state_control_CPU_CAP)
    inner = _state_control_STATE.evolve_prepared_field_jets
    field_runs = [0]

    def counted(*args, **kwargs):
        field_runs[0] += 1
        if field_runs[0] > 3:
            raise RuntimeError('more than three 801-grid field runs were requested')
        return inner(*args, **kwargs)
    _state_control_STATE.evolve_prepared_field_jets = counted
    arrays = {}
    try:
        x, phi, data, records, guard, retarded = _state_control_load_retarded_inputs()
        meta['support'] = guard
        meta['reused_artifacts'] = [r['payload'] for r in records] + [retarded['payload']]
        meta['cache'] = _state_control_probe_cache(guard)
        meta['retarded_restriction'] = _state_control_retarded_restriction_replay(retarded, data)
        if x.shape != (_state_control_POINTS,) or phi.shape != (1602, 12):
            raise ValueError('authenticated group14 801-grid columns required')
        owner = _state_control_FourthOrderModePropagator(x, data['mass'].item(), data['angular'].item())
        source = _state_control_FixedSourcePreparation.from_signed_blocks(data['energies'], data['weights'], data['source'])
        times = _state_control_np.linspace(0.0, 0.3, _state_control_STEPS + 1)
        print('support and cache checked; three prescribed 801/64 solves only', flush=True)
        jobs = (('center', _state_control_AMPLITUDE, 1), ('plus', _state_control_AMPLITUDE + _state_control_STEP, 0), ('minus', _state_control_AMPLITUDE - _state_control_STEP, 0))
        states = {}
        for name, amplitude, ndir in jobs:
            begin = _state_control_time.process_time()
            states[name] = _state_control_one_state(amplitude, ndir, owner, phi, data, times, source, guard)
            meta['case_CPU_seconds'][name] = _state_control_time.process_time() - begin
            meta['completed'].append(name)
            print(name + ' checkpoint, CPU ' + str(meta['case_CPU_seconds'][name]), flush=True)
        meta['field_runs'] = field_runs[0]
        center, plus, minus = (states['center'], states['plus'], states['minus'])
        center.require_history(_state_control_family(_state_control_AMPLITUDE, guard), times, x)
        meta['algebra'] = {'norm_flux_residual': center.diagnostics['norm_flux_residual'], 'metric_flux_tangent_residual': center.diagnostics['metric_flux_tangent_residual'], 'auxiliary_identity_residual': center.diagnostics['auxiliary_identity_residual']}
        arrays.update({'source/covariance': _state_control_np.asarray(source.covariance), 'source/column_weights': _state_control_np.asarray(source.column_weights), 'source/energies': _state_control_np.asarray(source.energies), 'source/blocks': _state_control_np.asarray(data['source']), 'source/weights': _state_control_np.asarray(data['weights']), 'source/mass': _state_control_np.asarray(data['mass']), 'source/angular': _state_control_np.asarray(data['angular']), 'center/z': _state_control_np.asarray(center.z), 'center/columns': _state_control_np.asarray(center.columns), 'center/column_tangents': _state_control_np.asarray(center.column_tangents[0]), 'center/covariance': _state_control_np.asarray(center.covariance()), 'center/covariance_tangent': _state_control_np.asarray(center.covariance_tangent()[0]), 'center/field': _state_control_np.asarray(center.diagnostics['evolved_field']), 'center/field_tangent': _state_control_np.asarray(center.diagnostics['tangent_fields'][0]), 'center/restriction': _state_control_np.asarray(center.history.restriction), 'center/rho1_metric': _state_control_np.asarray(center.history.rho1_metric), 'plus/columns': _state_control_np.asarray(plus.columns), 'plus/field': _state_control_np.asarray(plus.diagnostics['evolved_field']), 'minus/columns': _state_control_np.asarray(minus.columns), 'minus/field': _state_control_np.asarray(minus.diagnostics['evolved_field']), 'times': times, 'x': _state_control_np.asarray(x)})
        meta.update(CPU_seconds=_state_control_time.process_time() - start_cpu, wall_seconds=_state_control_time.monotonic() - start_wall)
        if _state_control_signature() != before:
            raise ValueError('input/producer changed during bounded adapter control')
        if field_runs[0] != 3:
            raise RuntimeError('source-fixed control must execute exactly three field runs')
        return _state_control_publish(arrays, meta)
    except _state_control_CPUExceeded:
        meta['cpu_budget_exceeded'] = True
        meta.update(CPU_seconds=_state_control_time.process_time() - start_cpu, wall_seconds=_state_control_time.monotonic() - start_wall)
        raise
    finally:
        _state_control_signal.setitimer(_state_control_signal.ITIMER_PROF, 0)
        _state_control_signal.signal(_state_control_signal.SIGPROF, old_handler)
        _state_control_STATE.evolve_prepared_field_jets = inner

def _state_control_covariance_roundoff_envelope(arrays):
    """Bound two binary64 replays of the saved covariance residual.

    The immutable weighted columns and source covariance define K=F C F*.
    For an n-term complex dot, gamma_(8n) bounds real arithmetic in a
    multiply/add reduction, including a different BLAS summation order.
    Two dots give (2*gamma+gamma**2)*|F||C||F*|.  Componentwise 1-norms
    upper-bound complex magnitudes here.  Propagate that forward error
    through (K_plus-K_minus)/(2*h), subtraction of the saved tangent, and
    the final complex modulus.  Twice the one-replay error bounds two
    implementations.  This arithmetic envelope is independent of FIELD_TOL.
    """
    np = _state_control_np
    eps = np.finfo(np.float64).eps  # Conservative: twice unit roundoff.
    C = np.asarray(arrays['source/covariance'])
    weights = np.asarray(arrays['source/column_weights'])
    n = C.shape[0]
    gamma = 8 * n * eps / (1 - 8 * n * eps)
    norm1 = lambda z: np.abs(np.real(z)) + np.abs(np.imag(z))
    scales = []
    for name in ('plus', 'minus'):
        F = (np.asarray(arrays[name + '/columns']) * weights).reshape(-1, n)
        magnitude = norm1(F) @ norm1(C) @ norm1(F).T
        # The positive reductions used to obtain the bound also round.
        scales.append(float(np.max(magnitude)) * (1 + gamma)**2)
    scale = sum(scales)
    tangent_scale = float(np.max(norm1(arrays['center/covariance_tangent'])))
    denominator = 2 * _state_control_STEP
    kernel_error = (2 * gamma + gamma**2) * scale
    subtraction_error = kernel_error + eps * (1 + gamma)**2 * scale
    division_error = ((1 + eps) * subtraction_error + eps * scale) / denominator
    derivative_scale = scale / denominator
    tangent_error = (1 + eps) * division_error + eps * (derivative_scale + tangent_scale)
    modulus_error = 8 * eps * (derivative_scale + tangent_scale + tangent_error)
    envelope = 2 * (tangent_error + modulus_error)
    if not np.isfinite(envelope) or envelope < 0:
        raise ValueError('nonfinite covariance replay arithmetic envelope')
    return float(np.nextafter(envelope, np.inf))


def _state_control_validate_replay(old, record, arrays):
    """Keep all receipt fields exact except the array-derived BLAS residual."""
    field = 'complete_coherent_covariance_derivative'
    strict = dict(record)
    strict['residuals'] = dict(record['residuals'])
    strict['residuals'][field] = old['residuals'][field]
    from .nsc_weighted_receipt_replay import exact_record_equal
    if not exact_record_equal(strict,old):
        raise ValueError('evolved incoming-state receipt differs from replay')
    stored, computed = old['residuals'][field], record['residuals'][field]
    tolerance = old['control']['covariance_derivative_tolerance']
    if not all(type(v) is float and _state_control_np.isfinite(v) and 0 <= v < tolerance
               for v in (stored, computed)):
        raise ValueError('evolved covariance derivative fails unchanged scientific control')
    if abs(stored - computed) > _state_control_covariance_roundoff_envelope(arrays):
        raise ValueError('evolved covariance residual exceeds replay roundoff envelope')


def _state_control_check():
    old = _state_control_json.loads((_state_control_ROOT / _state_control_OUTPUT).read_text())
    if _state_control_digest(old['payload']['path']) != old['payload']['sha256']:
        raise ValueError('evolved incoming-state artifact changed')
    path = _state_control_ROOT / old['payload']['path']
    with _state_control_np.load(path, allow_pickle=False) as loaded:
        arrays = {k: loaded[k] for k in loaded.files}
    meta = _state_control_json.loads(arrays['metadata_json'].tobytes())
    current = _state_control_signature()
    if meta['verification_signature'] != current:
        raise ValueError('evolved incoming-state producer/input changed')
    meta['verification_signature'] = current
    record = _state_control_record_from_arrays(arrays, meta, _state_control_digest(old['payload']['path']), path.stat().st_size)
    _state_control_validate_replay(old, record, arrays)
    # The authenticated receipt remains the historical published value.
    return old

from types import SimpleNamespace as _state_control_Namespace
state_control = _state_control_Namespace(
    ALGEBRA_TOL=_state_control_ALGEBRA_TOL,
    AMPLITUDE=_state_control_AMPLITUDE,
    ARTIFACT=_state_control_ARTIFACT,
    AmplitudeOnlyMetric=_state_control_AmplitudeOnlyMetric,
    CPUExceeded=_state_control_CPUExceeded,
    CPU_CAP=_state_control_CPU_CAP,
    CachedCompatibleIncomingMetric=_state_control_CachedCompatibleIncomingMetric,
    CompatibleIncomingMetric=_state_control_CompatibleIncomingMetric,
    CompatibleRadiusDirection=_state_control_CompatibleRadiusDirection,
    FIELD_TOL=_state_control_FIELD_TOL,
    FixedSourcePreparation=_state_control_FixedSourcePreparation,
    FourthOrderModePropagator=_state_control_FourthOrderModePropagator,
    MAXABS=_state_control_MAXABS,
    OUTPUT=_state_control_OUTPUT,
    OWNED=_state_control_OWNED,
    OWNERS=_state_control_OWNERS,
    POINTS=_state_control_POINTS,
    Path=_state_control_Path,
    RECORDS=_state_control_RECORDS,
    RETARDED=_state_control_RETARDED,
    ROOT=_state_control_ROOT,
    STATE=_state_control_STATE,
    STEP=_state_control_STEP,
    STEPS=_state_control_STEPS,
    analyze=_state_control_analyze,
    argparse=_state_control_argparse,
    check=_state_control_check,
    covariance_roundoff_envelope=_state_control_covariance_roundoff_envelope,
    deterministic_npz_bytes=_state_control_deterministic_npz_bytes,
    digest=_state_control_digest,
    evolve_incoming_state=_state_control_evolve_incoming_state,
    family=_state_control_family,
    incident=_state_control_incident,
    json=_state_control_json,
    kernel=_state_control_kernel,
    load_retarded_inputs=_state_control_load_retarded_inputs,
    np=_state_control_np,
    one_state=_state_control_one_state,
    os=_state_control_os,
    probe_cache=_state_control_probe_cache,
    publish=_state_control_publish,
    record_from_arrays=_state_control_record_from_arrays,
    restriction_covariance_tangent=_state_control_restriction_covariance_tangent,
    retarded_restriction_replay=_state_control_retarded_restriction_replay,
    rho1_node_index=_state_control_rho1_node_index,
    rho1_restriction_map=_state_control_rho1_restriction_map,
    run=_state_control_run,
    sha256=_state_control_sha256,
    signal=_state_control_signal,
    signature=_state_control_signature,
    sys=_state_control_sys,
    time=_state_control_time,
    validate_replay=_state_control_validate_replay,
)
