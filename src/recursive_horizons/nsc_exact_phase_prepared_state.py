"""Exact harmonic source phases with node PDE axial derivatives.

Source-fixed incoming columns on a compatible radius family obey

    Phi' = L[g] Phi + B Q,
    Y_j' = L[g] Y_j + (d_j L[g]) Phi,
    Q'   = -i diag(E) Q,    Q(t_0) = I.

B is the existing SAT-weighted right-inflow amplitude of the initial
reference columns. The auxiliary Q is the original source-phase block inside
the sparse exponential, not a midpoint-held forcing, a free inflow, an
interior counterforce, or a physical degree of freedom. Midpoint actual L and
dL are used per step; F, F_z, dF and dF_z are restricted at each time node
from the node PDE. On rho1, z=tau+S(1) so these are spatial axial derivatives.
Initial and incident tangents are zero for the supplied fixed preparation.
The caller owns the retarded support argument; the initial-slice guard alone
does not prove equality on the entire past. Physical spectrum and CF4 remain
OPEN; this owner does not step the metric.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from types import MappingProxyType

import numpy as np
from scipy.sparse import bmat, csr_matrix, diags
from scipy.sparse.linalg import expm_multiply

from .nsc_evolved_incoming_state import (
    AmplitudeOnlyMetric,
    CachedCompatibleIncomingMetric,
    EvolvedIncomingState,
    FixedSourcePreparation,
    HistoryBinding,
    _as_compatible_family,
    _assert_fixed_intrinsic,
    _coverage,
    _prepare_provider,
    _reject_legacy_kwargs,
    rho1_node_index,
    rho1_restriction_map,
)
from .nsc_ks_spacetime_variation import chart_coordinates
from .nsc_transmitting_history_jets import generator_direction, incoming_reference_columns


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


def _finite_complex(value, name):
    return _finite_array(value, complex, name)


def _digest_arrays(*arrays):
    hasher = sha256()
    for value in arrays:
        array = np.ascontiguousarray(value)
        hasher.update(array.dtype.str.encode())
        hasher.update(str(array.shape).encode())
        hasher.update(b'\0')
        hasher.update(array.tobytes(order='C'))
    return hasher.hexdigest()


def _fixed_preparation_digest(source, initial_fields, x, mass, angular, t0):
    return sha256(
        source.digest.encode()
        + np.ascontiguousarray(initial_fields).tobytes()
        + np.ascontiguousarray(x).tobytes()
        + repr((float(mass), float(angular), float(t0))).encode()
    ).hexdigest()


def _require_fixed_source(source):
    if not isinstance(source, FixedSourcePreparation):
        if isinstance(source, dict) and (
                'stress_approximant' in source or 'parts' in source
                or 'incoming_C0' in source):
            raise ValueError(
                'a frozen incoming C0 or legacy matter dictionary is not an '
                'exact-phase incoming state')
        raise TypeError(
            'FixedSourcePreparation required; a frozen incoming C0 or '
            'covariance matrix is not an exact-phase incoming state')
    return source


def _require_reference_support(owner, provider, times, x):
    """Check the initial slice of the supplied retarded preparation family.

    Compact retarded radius support makes the entire earlier geometry and
    incident SAT data equal their reference values, so initial and incident
    tangents are identically zero. The caller must supply that compact-support
    argument. This guard checks only the sampled initial slice; it does not
    certify the earlier continuum geometry from a finite set of samples.
    """
    metric0 = _finite_real(provider.values(float(times[0]), x), 'initial history metric')
    if metric0.shape != np.asarray(owner.reference).shape:
        raise ValueError('initial metric must have the owner reference shape')
    if float(np.max(np.abs(metric0 - owner.reference))) > 3e-11:
        raise ValueError(
            'reference-support guard: initial metric must equal the unchanged '
            'reference geometry; free inflow/preparation is not supplied')
    directions0 = _finite_real(
        provider.log_directions(float(times[0]), x, metric0), 'initial metric directions')
    if directions0.ndim != 3 or directions0.shape[1:] != (4, len(x)):
        raise ValueError('metric log_directions must have shape (ndir, 4, n)')
    if directions0.size and float(np.max(np.abs(directions0))) > 3e-11:
        raise ValueError(
            'reference-support guard: initial preparation and incident tangents '
            'are zero by unchanged past/upstream support')
    return metric0, directions0


def _inflow_amplitude(owner, phi):
    forcing = np.asarray(incoming_reference_columns(owner, phi), complex)
    rows = 2 * len(owner.x)
    if forcing.shape != phi.shape:
        raise ValueError('SAT inflow amplitude must match the initial columns')
    outside = np.ones(rows, dtype=bool)
    outside[[len(owner.x) - 1, rows - 1]] = False
    if np.any(forcing[outside] != 0):
        raise ValueError(
            'incident SAT amplitude must be supported only on the two right '
            'inflow rows; no interior counterforce is added')
    return forcing


def _sample_history_profiles(provider, family, times, x):
    """Node and midpoint metric/direction samples of the actual evaluated family.

    Equal amplitudes and window radii with a different callback are distinct
    if these samples differ. Discrete samples do not identify an arbitrary
    callable or prove continuum equality.
    """
    times = np.asarray(times, float)
    midpoints = 0.5 * (times[:-1] + times[1:])
    full_provider = provider.inner if isinstance(provider, AmplitudeOnlyMetric) else provider
    direction_provider = full_provider if isinstance(full_provider, CachedCompatibleIncomingMetric) else family
    node_metrics = []
    node_directions = []
    for time in times:
        metric = np.asarray(provider.values(float(time), x), float)
        if metric.shape != (4, len(x)):
            raise ValueError('metric must have shape (4, len(x))')
        directions = np.asarray(direction_provider.log_directions(float(time), x, metric), float)
        node_metrics.append(metric)
        node_directions.append(directions)
    mid_metrics = []
    mid_directions = []
    for time in midpoints:
        metric = np.asarray(provider.values(float(time), x), float)
        directions = np.asarray(direction_provider.log_directions(float(time), x, metric), float)
        mid_metrics.append(metric)
        mid_directions.append(directions)
    node_metrics = np.stack(node_metrics)
    node_directions = np.stack(node_directions)
    mid_metrics = np.stack(mid_metrics)
    mid_directions = np.stack(mid_directions)
    fingerprint = _digest_arrays(
        times, midpoints, x, node_metrics, node_directions, mid_metrics, mid_directions)
    return {
        'node_metrics': node_metrics,
        'node_directions': node_directions,
        'midpoint_metrics': mid_metrics,
        'midpoint_directions': mid_directions,
        'midpoints': midpoints,
        'fingerprint': fingerprint,
    }


def _augmented_generator(L, dL_blocks, B, energies):
    """Sparse action of (Phi, Y_j, Q) with exact harmonic Q' = -i E Q."""
    ndir = len(dL_blocks)
    blocks = [[None] * (ndir + 2) for _ in range(ndir + 2)]
    blocks[0][0] = L
    blocks[0][-1] = B
    for index, dL in enumerate(dL_blocks):
        blocks[index + 1][0] = dL
        blocks[index + 1][index + 1] = L
    blocks[-1][-1] = diags(-1j * np.asarray(energies, float))
    return bmat(blocks, format='csr')


def _node_pde(owner, metric, directions, phi, tangents, Q, B):
    L, checks = owner.generator(metric)
    rhs = L @ phi + B @ Q
    ndir = tangents.shape[0]
    d_rhs = np.empty((ndir, *phi.shape), dtype=complex)
    max_dflux = 0.
    for index, direction in enumerate(directions):
        dL, error = generator_direction(owner, metric, direction)
        max_dflux = max(max_dflux, error)
        d_rhs[index] = dL @ phi + L @ tangents[index]
    return L, rhs, d_rhs, checks, max_dflux


def _restrict(restriction, take, field, rhs, tangents, d_rhs):
    F = restriction @ field[take]
    Fz = restriction @ rhs[take]
    dF = np.einsum('ab,dbc->dac', restriction, tangents[:, take])
    dFz = np.einsum('ab,dbc->dac', restriction, d_rhs[:, take])
    return F, Fz, dF, dFz


@dataclass(frozen=True)
class ExactPhaseIncoming:
    """Restricted exact-phase columns plus node PDE axial derivatives.

    Public arrays
    -------------
    state : EvolvedIncomingState
        F, dF, C_src, weights, z=tau+S(1), source-fixed history binding.
    axial_columns : (nz, 2, nsrc)
        Unweighted F_z = R (L_g Phi + B Q) at each time node.
    axial_tangents : (ndir, nz, 2, nsrc)
        Unweighted dF_z = R (dL_g Phi + L_g Y) at each time node.
    sampled_node_metrics : (nz, 4, n)
    sampled_midpoint_metrics : (nsteps, 4, n)
    sampled_node_directions : (nz, nfam, 4, n)
    sampled_midpoint_directions : (nsteps, nfam, 4, n)

    Weighted axial properties apply the existing source quadrature once.
    CF4, continuum error, full source coverage and metric evolution are OPEN.
    """

    state: object
    axial_columns: object
    axial_tangents: object
    sampled_node_metrics: object
    sampled_midpoint_metrics: object
    sampled_node_directions: object
    sampled_midpoint_directions: object
    sampled_history_fingerprint: str
    fixed_preparation_digest: str
    diagnostics: object

    def __post_init__(self):
        if not isinstance(self.state, EvolvedIncomingState):
            raise TypeError(
                'ExactPhaseIncoming holds an EvolvedIncomingState; a frozen '
                'incoming C0 is not this state')
        axial = _finite_complex(self.axial_columns, 'axial columns')
        axial_tangents = _finite_complex(self.axial_tangents, 'axial tangents')
        node_metrics = _finite_real(self.sampled_node_metrics, 'sampled node metrics')
        mid_metrics = _finite_real(self.sampled_midpoint_metrics, 'sampled midpoint metrics')
        node_dirs = _finite_array(
            self.sampled_node_directions, float, 'sampled node directions')
        mid_dirs = _finite_array(
            self.sampled_midpoint_directions, float, 'sampled midpoint directions')
        if axial.shape != self.state.columns.shape:
            raise ValueError('axial_columns must have shape (nz, 2, nsrc)')
        if axial_tangents.shape != self.state.column_tangents.shape:
            raise ValueError('axial_tangents must have shape (ndir, nz, 2, nsrc)')
        nz = len(self.state.z)
        n = len(self.state.history.spatial_grid)
        nfam = len(self.state.history.amplitudes)
        if node_metrics.shape != (nz, 4, n) or mid_metrics.shape != (nz - 1, 4, n):
            raise ValueError('sampled metrics must cover every node and midpoint')
        if (node_dirs.shape != (nz, nfam, 4, n)
                or mid_dirs.shape != (nz - 1, nfam, 4, n)):
            raise ValueError('sampled directions must cover the owned radius family')
        if not isinstance(self.sampled_history_fingerprint, str) or len(self.sampled_history_fingerprint) != 64:
            raise ValueError('sampled-history fingerprint required')
        if not isinstance(self.fixed_preparation_digest, str) or len(self.fixed_preparation_digest) != 64:
            raise ValueError('fixed-preparation digest required')
        diagnostics = self.diagnostics if self.diagnostics is not None else {}
        object.__setattr__(self, 'axial_columns', axial)
        object.__setattr__(self, 'axial_tangents', axial_tangents)
        object.__setattr__(self, 'sampled_node_metrics', node_metrics)
        object.__setattr__(self, 'sampled_midpoint_metrics', mid_metrics)
        object.__setattr__(self, 'sampled_node_directions', node_dirs)
        object.__setattr__(self, 'sampled_midpoint_directions', mid_dirs)
        object.__setattr__(self, 'diagnostics', MappingProxyType(dict(diagnostics)))
        self.validate()

    @property
    def weighted_columns(self):
        return self.state.weighted_columns

    @property
    def weighted_column_tangents(self):
        return self.state.weighted_column_tangents

    @property
    def weighted_axial_columns(self):
        """F_z with source quadrature applied once: shape (nz, 2, nsrc)."""
        return np.asarray(self.axial_columns) * np.asarray(self.state.column_weights)

    @property
    def weighted_axial_tangents(self):
        """dF_z with the same once-applied weights: shape (ndir, nz, 2, nsrc)."""
        return np.asarray(self.axial_tangents) * np.asarray(self.state.column_weights)

    def validate(self):
        """Finite bound exact-phase state with actual sampled history attached."""
        self.state.validate()
        nsrc = self.state.columns.shape[-1]
        if self.axial_columns.shape[-1] != nsrc or self.axial_tangents.shape[-1] != nsrc:
            raise ValueError('source columns were truncated')
        expected = _digest_arrays(
            self.state.history.times,
            0.5 * (self.state.history.times[:-1] + self.state.history.times[1:]),
            self.state.history.spatial_grid,
            self.sampled_node_metrics,
            self.sampled_node_directions,
            self.sampled_midpoint_metrics,
            self.sampled_midpoint_directions,
        )
        if expected != self.sampled_history_fingerprint:
            raise ValueError(
                'sampled-history fingerprint does not match the stored node/'
                'midpoint metric and direction profiles')
        if self.fixed_preparation_digest != self.diagnostics.get('fixed_preparation_digest'):
            raise ValueError('fixed-preparation digest is not bound to the stored diagnostics')
        if self.state.history.source_digest != _digest_arrays(
                self.state.source_covariance, self.state.column_weights,
                self.state.source_energies):
            raise ValueError('source binding changed')
        if self.diagnostics.get('continuum_error_bound') not in (None, 'OPEN'):
            raise ValueError('continuum error remains OPEN; it is not a boolean proof')
        if self.diagnostics.get('physical_history_status') != 'OPEN':
            raise ValueError('physical history solution remains OPEN')
        if self.diagnostics.get('CF4_time_integrator') not in (False, None, 'not implemented'):
            raise ValueError('CF4 is not implemented in this owner')
        if self.diagnostics.get('metric_evolution') not in (False, None):
            raise ValueError('this owner does not step the metric')
        if self.diagnostics.get('outgoing_momentum_k_equals_minus_E') not in (False, None):
            raise ValueError('source energies are labels; k=-E is not a column momentum')
        if self.diagnostics.get('midpoint_held_forcing'):
            raise ValueError('exact harmonic phase block is required; midpoint-held forcing is not this owner')
        if self.diagnostics.get('callback_identity_claimed'):
            raise ValueError('sampled profiles do not identify an arbitrary callback')
        return True

    def require_history(self, provider, times, x):
        """Reject reuse unless amplitudes, windows and actual sampled profiles match.

        Same amplitudes and window radii with a different evaluated metric or
        direction callback cannot reuse this state. Agreement of these finite
        node/midpoint samples is not continuum callback identity.
        """
        self.validate()
        self.state.require_history(provider, times, x)
        family = _as_compatible_family(provider)
        if family is None:
            raise ValueError(
                'exact-phase incoming is bound to its CompatibleIncomingMetric '
                'family; it is not a reusable covariance matrix')
        prepared = provider
        if isinstance(provider, AmplitudeOnlyMetric):
            prepared = provider.inner
        if isinstance(prepared, CachedCompatibleIncomingMetric):
            if not np.array_equal(prepared.x, self.state.history.spatial_grid):
                raise ValueError('cached provider grid must match the bound spatial grid')
        else:
            prepared = CachedCompatibleIncomingMetric(self.state.history.spatial_grid, family)
        profiles = _sample_history_profiles(
            prepared,
            family, self.state.history.times, self.state.history.spatial_grid)
        residual = max(
            float(np.max(np.abs(profiles['node_metrics'] - self.sampled_node_metrics))),
            float(np.max(np.abs(profiles['midpoint_metrics'] - self.sampled_midpoint_metrics))),
            float(np.max(np.abs(profiles['node_directions'] - self.sampled_node_directions))),
            float(np.max(np.abs(profiles['midpoint_directions'] - self.sampled_midpoint_directions))),
        )
        if residual > 3e-12 or profiles['fingerprint'] != self.sampled_history_fingerprint:
            raise ValueError(
                'exact-phase incoming is bound to its actual sampled metric/'
                'direction history; equal amplitudes and windows with a '
                'different evaluated callback cannot reuse this state')
        return True


def prepare_exact_phase_incoming(owner, times, provider, initial_fields, source, **kwargs):
    """Evolve exact-phase source columns and restrict node PDE F, F_z, dF, dF_z.

    owner: FourthOrderModePropagator (SBP/SAT generator).
    times: strictly increasing PG-time nodes.
    provider: CompatibleIncomingMetric or AmplitudeOnlyMetric/cache wrapper.
    initial_fields: complex characteristic columns (2*len(x), nsrc).
    source: FixedSourcePreparation; there is no incoming C0 argument.

    Initial/incident tangents vanish by the fixed supplied preparation law.
    The guard checks its initial reference slice. The caller owns the earlier
    retarded-support argument; no free inflow callback is accepted.
    """
    _reject_legacy_kwargs(kwargs)
    source = _require_fixed_source(source)
    times = _finite_real(times, 'ordered caller time coordinates')
    if times.ndim != 1 or len(times) < 2 or np.any(np.diff(times) <= 0):
        raise ValueError('at least two strictly increasing caller time coordinates required')
    x = _finite_real(owner.x, 'owner spatial grid')
    if x.ndim != 1 or len(x) < 1:
        raise ValueError('one dimensional owner spatial grid required')
    phi = _finite_complex(initial_fields, 'initial fields')
    rows = 2 * len(x)
    nsrc = len(source.energies)
    if phi.ndim != 2 or phi.shape != (rows, nsrc):
        raise ValueError('initial fields must have shape (2*len(x), nsrc) matching the fixed source')
    provider = _prepare_provider(owner, provider)
    rho1 = rho1_node_index(x)
    restriction, rho1_metric = _assert_fixed_intrinsic(provider, times, x, rho1)
    rebuilt = rho1_restriction_map(
        rho1_metric[0], rho1_metric[2], rho1_metric[1], rho1_metric[3])
    if float(np.max(np.abs(rebuilt - restriction))) > 3e-12:
        raise ValueError('fixed intrinsic restriction must match the bound rho1 metric')
    family = _as_compatible_family(provider)
    if family is None:
        raise ValueError(
            'CompatibleIncomingMetric family required; unsupported restriction '
            'families are rejected rather than omitting the restriction derivative')
    _metric0, directions0 = _require_reference_support(owner, provider, times, x)
    ndir = directions0.shape[0]
    B_dense = _inflow_amplitude(owner, phi)
    B = csr_matrix(B_dense)
    energies = np.asarray(source.energies, float)
    take = np.array([rho1, len(x) + rho1], dtype=np.intp)
    tangents = np.zeros((ndir, rows, nsrc), dtype=complex)
    Q = np.eye(nsrc, dtype=complex)
    state = np.vstack((phi, *(tangents[j] for j in range(ndir)), Q))
    columns = np.empty((len(times), 2, nsrc), dtype=complex)
    axial = np.empty_like(columns)
    column_tangents = np.empty((ndir, len(times), 2, nsrc), dtype=complex)
    axial_tangents = np.empty_like(column_tangents)
    node_fields = np.empty((len(times), rows, nsrc), dtype=complex)
    node_tangents = np.empty((ndir, len(times), rows, nsrc), dtype=complex)
    phase = np.empty((len(times), nsrc, nsrc), dtype=complex)
    expected_phase = np.empty_like(phase)
    max_flux = max_dflux = max_speed = 0.
    identity = np.eye(nsrc, dtype=complex)

    def record(index, time, field, jets, phase_block):
        nonlocal max_flux, max_dflux, max_speed
        metric = _finite_real(provider.values(float(time), x), 'node metric')
        directions = _finite_real(
            provider.log_directions(float(time), x, metric), 'node metric directions')
        if directions.shape != (ndir, 4, len(x)):
            raise ValueError('metric direction count/shape must match the initial slice')
        _, rhs, d_rhs, checks, dflux = _node_pde(
            owner, metric, directions, field, jets, phase_block, B_dense)
        F, Fz, dF, dFz = _restrict(restriction, take, field, rhs, jets, d_rhs)
        columns[index] = F
        axial[index] = Fz
        column_tangents[:, index] = dF
        axial_tangents[:, index] = dFz
        node_fields[index] = field
        node_tangents[:, index] = jets
        phase[index] = phase_block
        expected_phase[index] = identity * np.exp(-1j * (float(time) - float(times[0])) * energies)
        max_flux = max(max_flux, checks['norm_flux_algebra'])
        max_dflux = max(max_dflux, dflux)
        max_speed = max(max_speed, checks['maximum_absolute_speed'])

    record(0, times[0], phi, tangents, Q)
    for step, (left, right) in enumerate(zip(times[:-1], times[1:], strict=True)):
        midpoint = (left + right) / 2
        metric = _finite_real(provider.values(midpoint, x), 'midpoint metric')
        if metric.shape != (4, len(x)):
            raise ValueError('metric must have shape (4, len(x))')
        L, checks = owner.generator(metric)
        if L.shape != (rows, rows):
            raise ValueError('owner generator dimension differs from characteristic columns')
        directions = _finite_real(
            provider.log_directions(midpoint, x, metric), 'midpoint metric directions')
        if directions.shape != (ndir, 4, len(x)):
            raise ValueError('metric direction count/shape must match the initial slice')
        dL_blocks = []
        for direction in directions:
            dL, error = generator_direction(owner, metric, direction)
            dL_blocks.append(dL)
            max_dflux = max(max_dflux, error)
        augmented = _augmented_generator(L, dL_blocks, B, energies)
        dt = right - left
        state = expm_multiply(dt * augmented, state, traceA=dt * augmented.diagonal().sum())
        if not np.isfinite(state).all():
            raise ArithmeticError('nonfinite exact-phase incoming evolution')
        field = state[:rows]
        jets = state[rows:(ndir + 1) * rows].reshape(ndir, rows, nsrc)
        phase_block = state[(ndir + 1) * rows:]
        max_flux = max(max_flux, checks['norm_flux_algebra'])
        max_speed = max(max_speed, checks['maximum_absolute_speed'])
        record(step + 1, right, field, jets, phase_block)

    profiles = _sample_history_profiles(provider, family, times, x)
    background_radius = np.sqrt(1. + np.asarray(x) * np.asarray(x))
    radius_change = max(
        float(np.max(np.abs(profiles['node_metrics'][:, 3] - background_radius))),
        float(np.max(np.abs(profiles['midpoint_metrics'][:, 3] - background_radius))),
    )
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
    preparation = _fixed_preparation_digest(
        source, phi, x, owner.mass, owner.angular, times[0])
    phase_residual = float(np.max(np.abs(phase - expected_phase)))
    energy_count = len(np.unique(np.round(source.energies, decimals=14)))
    state_diagnostics = {
        'norm_flux_residual': float(max_flux),
        'metric_flux_tangent_residual': float(max_dflux),
        'source_covariance_tangent_explicit': 0.,
        'restriction_parameter_derivative': 0.,
        'stationary_phase_reapplied': False,
        'incoming_trace_interpolation': None,
        'outgoing_momentum_k_equals_minus_E': False,
        'physical_preparation_status': 'OPEN',
        'physical_history_status': 'OPEN',
        'midpoint_steps': int(len(times) - 1),
        'direction_count': int(ndir),
        'source_column_count': int(nsrc),
        'evolved_field': _finite_array(state[:rows], complex, 'evolved field'),
        'tangent_fields': _finite_array(
            state[rows:(ndir + 1) * rows].reshape(ndir, rows, nsrc), complex, 'tangent fields'),
        'exact_harmonic_phase': True,
        'midpoint_held_forcing': False,
    }
    evolved = EvolvedIncomingState(
        z=z,
        columns=columns,
        column_tangents=column_tangents,
        source_covariance=source.covariance,
        column_weights=source.column_weights,
        source_energies=source.energies,
        history=history,
        coverage=_coverage(nsrc, energy_count),
        diagnostics=state_diagnostics,
    )
    diagnostics = {
        'norm_flux_residual': float(max_flux),
        'metric_flux_tangent_residual': float(max_dflux),
        'exact_harmonic_phase_residual': phase_residual,
        'phase': _finite_array(phase, complex, 'harmonic phase'),
        'expected_phase': _finite_array(expected_phase, complex, 'exact harmonic phase'),
        'node_fields': _finite_array(node_fields, complex, 'node fields'),
        'node_tangents': _finite_array(node_tangents, complex, 'node tangents'),
        'maximum_sampled_characteristic_speed': float(max_speed),
        'maximum_time_step': float(np.max(np.diff(times))),
        'midpoint_steps': int(len(times) - 1),
        'direction_count': int(ndir),
        'source_column_count': int(nsrc),
        'restriction_parameter_derivative': 0.,
        'source_covariance_tangent_explicit': 0.,
        'initial_and_incident_tangents': (
            'zero by fixed supplied preparation; earlier retarded support is a caller contract'),
        'inflow_convention': (
            'SAT-weighted right rows n-1 and 2n-1 from incoming_reference_columns; '
            'no free inflow callback and no interior counterforce'),
        'temporal_prescription': (
            'piecewise midpoint actual L and dL; exact harmonic phase block '
            'inside sparse expm; node PDE restriction of F_z and dF_z'),
        'midpoint_held_forcing': False,
        'exact_harmonic_phase': True,
        'CF4_time_integrator': False,
        'stationary_phase_reapplied': False,
        'outgoing_momentum_k_equals_minus_E': False,
        'incoming_trace_interpolation': None,
        'sampled_history_fingerprint': profiles['fingerprint'],
        'fixed_preparation_digest': preparation,
        'maximum_radius_change_from_background': float(radius_change),
        'cache_class': type(provider).__name__,
        'physical_preparation_status': 'OPEN',
        'physical_history_status': 'OPEN',
        'continuum_error_bound': None,
        'full_source_error_bound': None,
        'full_spectral_covariance': None,
        'callback_identity_claimed': False,
        'metric_evolution': False,
        'optimizer': False,
        'source_window_calculation': False,
        'physical_constraint_status': 'OPEN',
    }
    return ExactPhaseIncoming(
        state=evolved,
        axial_columns=axial,
        axial_tangents=axial_tangents,
        sampled_node_metrics=profiles['node_metrics'],
        sampled_midpoint_metrics=profiles['midpoint_metrics'],
        sampled_node_directions=profiles['node_directions'],
        sampled_midpoint_directions=profiles['midpoint_directions'],
        sampled_history_fingerprint=profiles['fingerprint'],
        fixed_preparation_digest=preparation,
        diagnostics=diagnostics,
    )
