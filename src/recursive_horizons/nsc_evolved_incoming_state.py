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
