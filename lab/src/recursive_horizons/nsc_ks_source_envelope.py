"""Exact KS envelope re-representation of the approved local prepared Dirac law.

The source, metric family, Pauli matrices, chart and plateau owners are
unchanged. The new object integrates the same continuum generator in the
envelope X of F=exp(-i E z) X, using rho rather than inverted T. Source
energies label columns; they never set the axial mesh. The periodic Fourier
grid is numerical padding, not a physical Pi boundary condition.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from types import MappingProxyType

import numpy as np
from scipy.integrate import DOP853, quad

from .nsc_common_time_bulk_split import CommonTimeBulkSplit
from .nsc_compatible_history_geometry import (
    CompatibleIncomingMetric,
    CompatibleRadiusDirection,
    _plateau,
)
from .nsc_evolved_incoming_state import (
    FixedSourcePreparation,
    _digest_arrays,
    _finite_array,
    _finite_real,
)
from .nsc_ks_spacetime_variation import chart_coordinates
from .nsc_local_incoming_family import LocalAxialFunction, LocalIncomingFamily
from .nsc_lorentzian import geometry
from .nsc_pg_ks_metric_pullback import reference_chart
from .nsc_prepared_history_jets import restriction_covariance_tangent
from .nsc_retarded_radial_response import canonical_map, ks_generator
from .nsc_transmitting_dirac_domain import (
    MODE_TO_CURRENT,
    S1,
    S2,
    S3,
    TransmittingDiracSeamDomain,
)


RHO_SIGMA = 1.0
RHO_UP_MIN = 1.03
DEFAULT_COMPUTATIONAL_LENGTH = 0.4
AXIAL_CENTER_OFFSET = 0.15
PHYSICAL_I_OFFSETS = (0.12, 0.18)
USUAL_SUPPORT_OFFSETS = (0.09, 0.21)
_FOREIGN = frozenset({
    'history', 'HistoryBinding', 'node_fields', 'incoming_C0', 'C0',
    'covariance', 'matter', 'optimizer', 'Gamma_rest', 'metric_step',
})


def _finite_complex(value, name):
    return _finite_array(value, complex, name)


def _as_metric(family):
    if isinstance(family, LocalIncomingFamily):
        return family, family.metric()
    if isinstance(family, CompatibleIncomingMetric):
        return None, family
    raise TypeError(
        'LocalIncomingFamily or CompatibleIncomingMetric required; '
        'HistoryBinding and PG node_fields are not this envelope')


def _scalar_pair(value, name):
    raw = _finite_real(value, name)
    if raw.shape != (2,) or not raw[0] < raw[1]:
        raise ValueError(f'ordered finite {name} pair required')
    return (float(raw[0]), float(raw[1]))


def _solver_option(value, name):
    if isinstance(value, (bool, np.bool_)) or value is None:
        raise ValueError(f'explicit finite positive {name} required')
    number = float(value)
    if not np.isfinite(number) or number <= 0:
        raise ValueError(f'explicit finite positive {name} required')
    return number


def _coverage(nsrc, energy_count):
    return MappingProxyType({
        'sampled_source': (
            f'{energy_count} signed source-energy labels, {nsrc} columns; '
            'not a complete quadrature'),
        'angular_multiplicity_added': False,
        'quadrature_applied_once': True,
        'continuum_error_bound': 'OPEN',
        'integration_error_bound': 'OPEN',
        'truncation_error_bound': 'OPEN',
        'physical_source_error': 'OPEN',
        'finite_fourier_causality': 'OPEN',
        'physical_history_solution': 'OPEN',
        'periodic_Pi_boundary_condition': False,
        'physical_confidence_from_tolerances': False,
        'PG_discretization_copied': False,
        'metric_evolution': False,
        'outgoing_momentum_substitution': None,
        'callback_identity_claimed': False,
        'source_frequency_controls_spatial_mesh': False,
    })


def physical_incoming_interval():
    """Fixed physical I = S(1)+[.12,.18]."""
    S1z = chart_coordinates(1.)[1]
    return (S1z + PHYSICAL_I_OFFSETS[0], S1z + PHYSICAL_I_OFFSETS[1])


def usual_axial_support():
    """Usual compact axial support S(1)+[.09,.21]."""
    S1z = chart_coordinates(1.)[1]
    return (S1z + USUAL_SUPPORT_OFFSETS[0], S1z + USUAL_SUPPORT_OFFSETS[1])


def axial_center():
    return chart_coordinates(1.)[1] + AXIAL_CENTER_OFFSET


def computational_z_grid(count, length=DEFAULT_COMPUTATIONAL_LENGTH):
    """Uniform periodic envelope nodes centered at S(1)+.15.

    The period is numerical padding. It is not a physical axial identification.
    """
    if isinstance(count, bool) or not isinstance(count, (int, np.integer)) or int(count) < 8:
        raise ValueError('at least eight computational envelope nodes required')
    length = float(length)
    if not np.isfinite(length) or length <= 0:
        raise ValueError('positive computational envelope length required')
    n = int(count)
    center = axial_center()
    return center + length * (np.arange(n) / n - 0.5)


def continuum_speed_distance(rho_up, rho_sigma=RHO_SIGMA):
    """D = integral_{rho_sigma}^{rho_up} d rho / a_ref^2, the axial characteristic length."""
    rho_up = float(rho_up)
    rho_sigma = float(rho_sigma)
    if not np.isfinite([rho_up, rho_sigma]).all() or rho_up <= rho_sigma:
        raise ValueError('upstream rho must exceed rho_sigma')
    reference_chart(rho_up)
    reference_chart(rho_sigma)

    def integrand(rho):
        axial = reference_chart(float(rho))[2]
        return 1. / (axial * axial)

    value, _ = quad(integrand, rho_sigma, rho_up, epsabs=2e-13, epsrel=2e-13)
    return float(value)


def frame_potential_identity(rho, mass, angular, *, extra_radii=(.5, 2., 3.7)):
    """R=inverse_trace/r depends only on beta/a; U† M U = -m S1 + ell/r S2."""
    rho = float(rho)
    mass = float(mass)
    angular = float(angular)
    if not np.isfinite([rho, mass, angular]).all():
        raise ValueError('finite frame-identity arguments required')
    beta = float(geometry(rho)[0])
    axial = np.sqrt(beta * beta - 1.)
    radius = float(np.sqrt(1. + rho * rho))
    mapped = canonical_map(rho)
    cancellation = 0.
    for scale in extra_radii:
        trial = radius * float(scale)
        if trial <= 0 or not np.isfinite(trial):
            raise ValueError('positive comparison radii required')
        domain = TransmittingDiracSeamDomain(
            lapse=1., radial_scale=1., shift=beta, radius=trial)
        _, inverse, _ = domain.trace_map(np.ones(1))
        cancellation = max(
            cancellation, float(np.max(np.abs(inverse[0] / trial - mapped))))
    potential = CommonTimeBulkSplit.local_potential(1., radius, mass, angular)
    ks_potential = -mass * S1 + angular / radius * S2
    residual = float(np.max(np.abs(
        MODE_TO_CURRENT.conj().T @ potential @ MODE_TO_CURRENT - ks_potential)))
    return MappingProxyType({
        'canonical_map': mapped,
        'radius_cancellation': cancellation,
        'potential_residual': residual,
        'beta': beta,
        'a': float(axial),
        'depends_only_on_beta_a': True,
    })


def _require_uniform_centered_grid(z):
    z = _finite_real(z, 'computational z-grid')
    if z.ndim != 1 or len(z) < 8:
        raise ValueError('at least eight computational envelope nodes required')
    step = np.diff(z)
    if np.any(step <= 0) or float(np.max(np.abs(step - step[0]))) > 1e-12:
        raise ValueError('strictly uniform increasing computational z-grid required')
    dz = float(step[0])
    length = dz * len(z)
    center = float(z[0] + 0.5 * length)
    if abs(center - axial_center()) > 1e-12:
        raise ValueError('computational z-grid must be centered at S(1)+.15')
    return z, dz, length


def envelope_z_derivative(envelope, z_grid):
    """Periodic Fourier derivative of the envelope, never of the F carrier."""
    z, dz, _ = _require_uniform_centered_grid(z_grid)
    field = np.asarray(envelope, complex)
    if field.shape[-1] != len(z) or not np.isfinite(field).all():
        raise ValueError('finite envelope whose last axis matches the z-grid required')
    wave = 2. * np.pi * np.fft.fftfreq(len(z), d=dz)
    return np.fft.ifft(1j * wave * np.fft.fft(field, axis=-1), axis=-1)


def trigonometric_polynomial(samples, z_grid, z_target, *, derivative=False):
    """Evaluate the DFT trigonometric polynomial, or its derivative, at targets."""
    z, dz, _ = _require_uniform_centered_grid(z_grid)
    target = _finite_real(z_target, 'target z')
    if target.ndim != 1 or len(target) < 1:
        raise ValueError('nonempty target z required')
    field = np.asarray(samples, complex)
    if field.shape[-1] != len(z) or not np.isfinite(field).all():
        raise ValueError('finite samples whose last axis matches the z-grid required')
    wave = 2. * np.pi * np.fft.fftfreq(len(z), d=dz)
    coeff = np.fft.fft(field, axis=-1)
    if derivative:
        coeff = 1j * wave * coeff
    phase = np.exp(1j * np.multiply.outer(target - float(z[0]), wave))
    return np.tensordot(coeff, phase, axes=(-1, -1)) / len(z)


def _sample_axial_profiles(directions, z):
    z = np.asarray(z, float)
    w = np.empty((len(directions), len(z)), float)
    U = np.empty_like(w)
    for index, direction in enumerate(directions):
        if not isinstance(direction, CompatibleRadiusDirection):
            raise TypeError('compatible radius directions required')
        for node, point in enumerate(z):
            try:
                w[index, node] = float(direction.w(float(point), 0))
                U[index, node] = float(direction.U(float(point), 0))
            except (TypeError, ValueError, OverflowError) as error:
                raise ValueError('w/U callbacks must supply finite order-zero samples') from error
    if not np.isfinite(w).all() or not np.isfinite(U).all():
        raise ValueError('finite w/U samples required')
    return w, U


def _fixed_preparation_digest(source, initial, mass, angular, rho_up, rho_sigma):
    return sha256(
        source.digest.encode()
        + np.ascontiguousarray(np.asarray(initial, complex)).tobytes()
        + repr((float(mass), float(angular), float(rho_up), float(rho_sigma))).encode()
    ).hexdigest()


def _profile_fingerprint(z, w, U, amplitudes, inner, outer, axial_support,
                         rho_up, rho_sigma, rtol, atol, max_step):
    return sha256(
        _digest_arrays(
            z, w, U, np.asarray(amplitudes, float),
            np.asarray(inner, float), np.asarray(outer, float)).encode()
        + repr((tuple(map(float, axial_support)), float(rho_up), float(rho_sigma),
                float(rtol), float(atol), float(max_step))).encode()
    ).hexdigest()


def _reject_foreign(kwargs):
    for name in kwargs:
        if name in _FOREIGN or name.lower() in {key.lower() for key in _FOREIGN}:
            raise ValueError(
                f'{name} is not an envelope argument; HistoryBinding, PG '
                'node_fields, new action, optimizer and metric steps are rejected')
    if kwargs:
        raise TypeError(f'unsupported envelope argument {next(iter(kwargs))}')


def _require_source(source):
    if not isinstance(source, FixedSourcePreparation):
        if isinstance(source, dict) and (
                'stress_approximant' in source or 'parts' in source
                or 'incoming_C0' in source):
            raise ValueError(
                'a frozen incoming C0 or legacy matter dictionary is not a '
                'KS envelope incoming state')
        raise TypeError(
            'FixedSourcePreparation required; a frozen incoming C0 or '
            'covariance matrix is not a KS envelope incoming state')
    return source


@dataclass(frozen=True)
class KSEnvelopeBinding:
    """Evaluated w/U samples, windows, grid, rho endpoints and solver options.

    Agreement of these samples is not identity of an arbitrary callback.
    """
    computational_z: object
    w_samples: object
    U_samples: object
    amplitudes: tuple
    inner_radii: tuple
    outer_radii: tuple
    axial_support: tuple
    rho_up: float
    rho_sigma: float
    rtol: float
    atol: float
    max_step: float
    fingerprint: str = ''

    def __post_init__(self):
        z = _finite_real(self.computational_z, 'computational z-grid')
        w = _finite_real(self.w_samples, 'w samples')
        U = _finite_real(self.U_samples, 'U samples')
        amplitudes = tuple(float(v) for v in self.amplitudes)
        inner = tuple(float(v) for v in self.inner_radii)
        outer = tuple(float(v) for v in self.outer_radii)
        if w.ndim != 2 or w.shape != U.shape or w.shape[0] != len(amplitudes) or w.shape[1] != len(z):
            raise ValueError('w/U samples must have shape (ndir, n_computational)')
        if len(inner) != len(amplitudes) or len(outer) != len(amplitudes):
            raise ValueError('one inner/outer normal window per amplitude required')
        axial = _scalar_pair(self.axial_support, 'axial_support')
        fingerprint = _profile_fingerprint(
            z, w, U, amplitudes, inner, outer, axial,
            self.rho_up, self.rho_sigma, self.rtol, self.atol, self.max_step)
        if self.fingerprint and self.fingerprint != fingerprint:
            raise ValueError('envelope profile fingerprint does not match the bound samples')
        object.__setattr__(self, 'computational_z', z)
        object.__setattr__(self, 'w_samples', w)
        object.__setattr__(self, 'U_samples', U)
        object.__setattr__(self, 'amplitudes', amplitudes)
        object.__setattr__(self, 'inner_radii', inner)
        object.__setattr__(self, 'outer_radii', outer)
        object.__setattr__(self, 'axial_support', axial)
        object.__setattr__(self, 'rho_up', float(self.rho_up))
        object.__setattr__(self, 'rho_sigma', float(self.rho_sigma))
        object.__setattr__(self, 'rtol', float(self.rtol))
        object.__setattr__(self, 'atol', float(self.atol))
        object.__setattr__(self, 'max_step', float(self.max_step))
        object.__setattr__(self, 'fingerprint', fingerprint)


@dataclass(frozen=True)
class KSEnvelopeIncoming:
    """Incoming restriction of the KS envelope at rho_sigma=1.

    Public arrays use target z. F=e^{-i E z} X and
    F_z=e^{-i E z}(X_z-i E X); the carrier is never FFT-differentiated.
    """
    z: object
    columns: object
    column_tangents: object
    axial_columns: object
    axial_tangents: object
    source_covariance: object
    column_weights: object
    source_energies: object
    mass: float
    angular: float
    initial_columns: object
    rho_up: float
    binding: KSEnvelopeBinding
    fixed_preparation_digest: str
    coverage: object
    diagnostics: object

    def __post_init__(self):
        z = _finite_real(self.z, 'target z')
        columns = _finite_complex(self.columns, 'envelope columns')
        tangents = _finite_complex(self.column_tangents, 'envelope column tangents')
        axial = _finite_complex(self.axial_columns, 'envelope F_z')
        axial_tangents = _finite_complex(self.axial_tangents, 'envelope dF_z')
        covariance = _finite_complex(self.source_covariance, 'source covariance')
        weights = _finite_real(self.column_weights, 'column weights')
        energies = _finite_real(self.source_energies, 'source energy labels')
        initial = _finite_complex(self.initial_columns, 'upstream columns')
        if not isinstance(self.binding, KSEnvelopeBinding):
            raise TypeError('KSEnvelopeBinding required; HistoryBinding is not this envelope')
        if z.ndim != 1 or len(z) < 1:
            raise ValueError('nonempty target z required')
        if columns.ndim != 3 or columns.shape[0] != len(z) or columns.shape[1] != 2:
            raise ValueError('columns must have shape (nz, 2, nsrc)')
        nsrc = columns.shape[2]
        if tangents.ndim != 4 or tangents.shape[1:] != columns.shape:
            raise ValueError('column_tangents must have shape (ndir, nz, 2, nsrc)')
        if axial.shape != columns.shape or axial_tangents.shape != tangents.shape:
            raise ValueError('F_z and dF_z must match F and dF shapes')
        if (covariance.shape != (nsrc, nsrc) or weights.shape != (nsrc,)
                or energies.shape != (nsrc,) or initial.shape != (2, nsrc)):
            raise ValueError('source matrix, weights, labels and A_up must match nsrc')
        coverage = self.coverage if self.coverage is not None else _coverage(nsrc, nsrc)
        diagnostics = self.diagnostics if self.diagnostics is not None else {}
        object.__setattr__(self, 'z', z)
        object.__setattr__(self, 'columns', columns)
        object.__setattr__(self, 'column_tangents', tangents)
        object.__setattr__(self, 'axial_columns', axial)
        object.__setattr__(self, 'axial_tangents', axial_tangents)
        object.__setattr__(self, 'source_covariance', covariance)
        object.__setattr__(self, 'column_weights', weights)
        object.__setattr__(self, 'source_energies', energies)
        object.__setattr__(self, 'mass', float(self.mass))
        object.__setattr__(self, 'angular', float(self.angular))
        object.__setattr__(self, 'initial_columns', initial)
        object.__setattr__(self, 'rho_up', float(self.rho_up))
        object.__setattr__(self, 'coverage', MappingProxyType(dict(coverage)))
        object.__setattr__(self, 'diagnostics', MappingProxyType(dict(diagnostics)))
        self.validate()

    @property
    def weighted_columns(self):
        return np.asarray(self.columns) * np.asarray(self.column_weights)

    @property
    def weighted_column_tangents(self):
        return np.asarray(self.column_tangents) * np.asarray(self.column_weights)

    @property
    def weighted_axial_columns(self):
        return np.asarray(self.axial_columns) * np.asarray(self.column_weights)

    @property
    def weighted_axial_tangents(self):
        return np.asarray(self.axial_tangents) * np.asarray(self.column_weights)

    def validate(self):
        nsrc = self.columns.shape[-1]
        if self.binding.rho_up != self.rho_up:
            raise ValueError('bound rho_up does not match the stored upstream coordinate')
        source_digest = _digest_arrays(
            self.source_covariance, self.column_weights, self.source_energies)
        expected_prep = sha256(
            source_digest.encode()
            + np.ascontiguousarray(self.initial_columns).tobytes()
            + repr((float(self.mass), float(self.angular), float(self.rho_up),
                    float(self.binding.rho_sigma))).encode()
        ).hexdigest()
        if expected_prep != self.fixed_preparation_digest:
            raise ValueError('fixed-preparation digest does not match the bound source/A_up/m/ell/rho')
        if self.diagnostics.get('fixed_preparation_digest') not in (None, self.fixed_preparation_digest):
            raise ValueError('fixed-preparation digest is not bound to the stored diagnostics')
        rebuilt = _profile_fingerprint(
            self.binding.computational_z, self.binding.w_samples, self.binding.U_samples,
            self.binding.amplitudes, self.binding.inner_radii, self.binding.outer_radii,
            self.binding.axial_support, self.binding.rho_up, self.binding.rho_sigma,
            self.binding.rtol, self.binding.atol, self.binding.max_step)
        if rebuilt != self.binding.fingerprint:
            raise ValueError('envelope profile fingerprint does not match the bound samples')
        if np.any(self.column_weights <= 0):
            raise ValueError('positive column weights required')
        if float(np.max(np.abs(self.source_covariance - self.source_covariance.conj().T))) > 3e-11:
            raise ValueError('source covariance must be Hermitian')
        if self.coverage.get('continuum_error_bound') not in (None, 'OPEN'):
            raise ValueError('continuum error remains OPEN')
        if self.coverage.get('integration_error_bound') not in (None, 'OPEN'):
            raise ValueError('integration error remains OPEN')
        if self.coverage.get('truncation_error_bound') not in (None, 'OPEN'):
            raise ValueError('Fourier truncation error remains OPEN')
        if self.coverage.get('physical_source_error') not in (None, 'OPEN'):
            raise ValueError('physical source error remains OPEN')
        if self.coverage.get('finite_fourier_causality') not in (None, 'OPEN'):
            raise ValueError('finite Fourier-grid causality remains OPEN')
        if self.coverage.get('periodic_Pi_boundary_condition') not in (False, None):
            raise ValueError('the periodic envelope is numerical padding, not physical Pi BC')
        if self.coverage.get('physical_confidence_from_tolerances') not in (False, None):
            raise ValueError('caller tolerances are not a physical confidence certificate')
        if self.coverage.get('PG_discretization_copied') not in (False, None):
            raise ValueError('this integrator does not copy the PG CF4 discretization')
        if self.coverage.get('callback_identity_claimed') not in (False, None):
            raise ValueError('sampled profiles do not identify an arbitrary callback')
        if self.coverage.get('source_frequency_controls_spatial_mesh') not in (False, None):
            raise ValueError('source energies are column labels, not the axial mesh')
        if self.coverage.get('outgoing_momentum_substitution') not in (None, False, 'not used'):
            raise ValueError('source energies are labels; k=-E is not a column momentum')
        if nsrc != self.source_covariance.shape[0]:
            raise ValueError('source columns were truncated')
        return True

    def covariance(self):
        self.validate()
        nz, _, nsrc = self.columns.shape
        F = self.weighted_columns.reshape(nz * 2, nsrc)
        kernel = F @ self.source_covariance @ F.conj().T
        return kernel.reshape(nz, 2, nz, 2)

    def covariance_tangent(self):
        self.validate()
        nz, _, nsrc = self.columns.shape
        ndir = self.column_tangents.shape[0]
        if ndir == 0:
            return np.zeros((0, nz, 2, nz, 2), complex)
        F = self.weighted_columns.reshape(nz * 2, nsrc)
        dF = self.weighted_column_tangents.reshape(ndir, nz * 2, nsrc)
        dC = np.zeros((ndir, nsrc, nsrc), complex)
        kernel = restriction_covariance_tangent(F, dF, self.source_covariance, dC)
        return np.asarray(kernel).reshape(ndir, nz, 2, nz, 2)

    def require_history(self, family, z_grid=None, rho_up=None, solver_options=None):
        """Reject changed evaluated w/U samples, amplitudes, windows or endpoints.

        Matching these finite samples is not callback identity.
        """
        self.validate()
        _, metric = _as_metric(family)
        if tuple(float(v) for v in metric.amplitudes) != self.binding.amplitudes:
            raise ValueError('envelope is bound to its generating amplitudes')
        inner = tuple(float(d.inner_radius) for d in metric.directions)
        outer = tuple(float(d.outer_radius) for d in metric.directions)
        if inner != self.binding.inner_radii or outer != self.binding.outer_radii:
            raise ValueError('envelope is bound to its compact normal windows')
        w, U = _sample_axial_profiles(metric.directions, self.binding.computational_z)
        residual = max(float(np.max(np.abs(w - self.binding.w_samples), initial=0.)),
                       float(np.max(np.abs(U - self.binding.U_samples), initial=0.)))
        if residual > 3e-12:
            raise ValueError(
                'evaluated w/U profiles changed; sampled agreement is not '
                'callback identity')
        if z_grid is not None and not np.array_equal(
                np.asarray(z_grid, float), self.binding.computational_z):
            raise ValueError('envelope is bound to its computational z-grid')
        if rho_up is not None and float(rho_up) != self.binding.rho_up:
            raise ValueError('envelope is bound to its rho endpoints')
        if solver_options is not None:
            rtol, atol, max_step = (float(solver_options[name]) for name in ('rtol', 'atol', 'max_step'))
            if (rtol, atol, max_step) != (self.binding.rtol, self.binding.atol, self.binding.max_step):
                raise ValueError('envelope is bound to its supplied solver options')
        return True


def _radius_lower_bound(local, metric, w, U):
    if local is not None:
        return float(local.radius_lower_bound()), True
    if not len(metric.amplitudes):
        return 1., True
    if all(amp == 0. or (isinstance(d.w, LocalAxialFunction) and isinstance(d.U, LocalAxialFunction))
           for amp, d in zip(metric.amplitudes, metric.directions)):
        perturbation = Fraction(0)
        for amp, d in zip(metric.amplitudes, metric.directions):
            if amp == 0.:
                continue
            outer = Fraction(d.outer_radius)
            W = sum(abs(Fraction(c)) for c in d.w.coefficients)
            V = sum(abs(Fraction(c)) for c in d.U.coefficients)
            perturbation += abs(Fraction(amp))*(outer*W+outer**3*V/6)
        return float(np.nextafter(float(1-perturbation), -np.inf)), True
    outer = max(metric.directions[i].outer_radius for i in range(len(metric.amplitudes)))
    perturbation = 0.
    for amp, w_row, U_row in zip(metric.amplitudes, w, U):
        perturbation += abs(amp) * (outer * float(np.max(np.abs(w_row)))
                                    + (outer ** 3) * float(np.max(np.abs(U_row))) / 6.)
    # Arbitrary callbacks have no certified between-node supremum here.
    return float(np.nextafter(1. - perturbation, -np.inf)), False


def evolve_ks_source_envelope(
        source, initial_columns, family, z_grid, target_z, mass, angular, rho_up, *,
        axial_support, rtol, atol, max_step, tangents='all', rho_sigma=RHO_SIGMA, **kwargs):
    """Integrate the owned KS generator plus axial envelope derivative and radius correction.

    X_rho = (S3/a^2) X_z + G_E_ref(rho) X + (i ell/a)(1/r_g-1/r_ref) S2 X
    from rho_up down to rho_sigma=1. A_up is broadcast constant in z.
    """
    _reject_foreign(kwargs)
    source = _require_source(source)
    if float(rho_sigma) != RHO_SIGMA:
        raise ValueError('incoming surface is the owned rho_sigma=1 slice')
    rho_up = float(rho_up)
    if not np.isfinite(rho_up) or rho_up < RHO_UP_MIN:
        raise ValueError('rho_up must be at least 1.03, outside the usual normal window')
    mass = float(mass)
    angular = float(angular)
    if not np.isfinite([mass, angular]).all() or mass < 0:
        raise ValueError('finite mass and angular eigenvalue required')
    rtol = _solver_option(rtol, 'rtol')
    atol = _solver_option(atol, 'atol')
    max_step = _solver_option(max_step, 'max_step')
    if tangents not in ('all', 'zero'):
        raise ValueError("tangents must be 'all' or 'zero'")
    z, dz, length = _require_uniform_centered_grid(z_grid)
    target = _finite_real(target_z, 'target z')
    if target.ndim != 1 or len(target) < 1:
        raise ValueError('nonempty target z required')
    interval = physical_incoming_interval()
    support = _scalar_pair(axial_support, 'axial_support')
    if not (support[0] <= interval[0] and interval[1] <= support[1]):
        raise ValueError('fixed physical I must lie inside the explicit axial_support')
    domain_end = float(z[0] + length)
    if interval[0] < float(z[0]) or interval[1] > domain_end - 1e-12:
        raise ValueError('fixed physical I must lie inside the computational envelope domain')
    if np.any(target < interval[0] - 1e-12) or np.any(target > interval[1] + 1e-12):
        raise ValueError('target z must lie inside the fixed physical interval I')
    if np.any(target < float(z[0])) or np.any(target >= domain_end):
        raise ValueError('target z must lie inside the computational envelope domain')
    local, metric = _as_metric(family)
    if local is not None:
        if max(abs(local.interval[0] - interval[0]), abs(local.interval[1] - interval[1])) > 1e-12:
            raise ValueError('local family interval must be the fixed physical I')
        if max(abs(local.support[0] - support[0]), abs(local.support[1] - support[1])) > 1e-12:
            raise ValueError('explicit axial_support must match the local family support')
    initial = _finite_complex(initial_columns, 'upstream columns')
    nsrc = source.covariance.shape[0]
    if initial.shape != (2, nsrc):
        raise ValueError('initial columns must have shape (2, nsrc)')
    energies = np.asarray(source.energies, float)
    w_samples, U_samples = _sample_axial_profiles(metric.directions, z)
    outside = (z < support[0]) | (z > support[1])
    if np.any(outside):
        spilled = max(float(np.max(np.abs(w_samples[:, outside]))),
                      float(np.max(np.abs(U_samples[:, outside])))) if len(metric.directions) else 0.
        if spilled > 1e-12:
            raise ValueError('w/U samples must vanish outside the explicit axial_support')
    bound, global_radius_bound = _radius_lower_bound(local, metric, w_samples, U_samples)
    if bound <= 0:
        raise ValueError('coefficient bound does not keep the radius positive')
    T1 = chart_coordinates(1.)[0]
    s_up = chart_coordinates(rho_up)[0] - T1
    inners = tuple(float(d.inner_radius) for d in metric.directions)
    outers = tuple(float(d.outer_radius) for d in metric.directions)
    for outer in outers:
        if abs(s_up) < outer:
            raise ValueError(
                'normal support must lie outside the initial rho so the '
                'upstream slice remains the unperturbed preparation')
    reference_chart(rho_up)
    reference_chart(RHO_SIGMA)
    distance = continuum_speed_distance(rho_up)
    padding_left = support[0] - float(z[0])
    padding_right = domain_end - support[1]
    if not (padding_left > 2. * distance and padding_right > 2. * distance):
        raise ValueError(
            'computational padding beyond axial_support must exceed 2D '
            f'(conservative); D={distance}, padding=({padding_left}, {padding_right})')
    nz = len(z)
    amplitudes = np.asarray(metric.amplitudes, float)
    ndir_geom = len(amplitudes)
    evolve_tangents = tangents == 'all'
    ndir = ndir_geom if evolve_tangents else 0
    X0 = np.broadcast_to(initial[:, :, None], (2, nsrc, nz)).copy()
    if evolve_tangents:
        y0 = np.concatenate((X0.ravel(), np.zeros(ndir * 2 * nsrc * nz, complex)))
    else:
        y0 = X0.ravel()
    wave = 2. * np.pi * np.fft.fftfreq(nz, d=dz)
    amps = amplitudes

    def radius_state(rho):
        s = chart_coordinates(float(rho))[0] - T1
        r_ref = np.sqrt(1. + float(rho) * float(rho))
        if ndir_geom == 0:
            return np.full(nz, r_ref), np.zeros((0, nz)), r_ref
        chi = np.array([
            float(_plateau(np.asarray([s]), inner, outer)[0])
            for inner, outer in zip(inners, outers)])
        basis = chi[:, None] * (s * w_samples + (s ** 3) * U_samples / 6.)
        radius = r_ref + amps @ basis
        if (not np.isfinite(radius).all()) or np.any(radius <= 0):
            raise ArithmeticError('actual radius must stay finite and positive')
        return radius, basis, r_ref

    def apply_L(field, G, axial, radius, r_ref):
        freq = np.fft.ifft(1j * wave * np.fft.fft(field, axis=-1), axis=-1)
        term_z = np.einsum('ab,...bsz->...asz', S3, freq) / (axial * axial)
        term_G = np.einsum('sab,...bsz->...asz', G, field)
        correction = (1j * angular / axial) * (1. / radius - 1. / r_ref)
        term_r = correction * np.einsum('ab,...bsz->...asz', S2, field)
        return term_z + term_G + term_r

    def rhs(rho, state):
        rho = float(rho)
        beta = float(geometry(rho)[0])
        axial = np.sqrt(beta * beta - 1.)
        if not np.isfinite(axial) or axial <= 0:
            raise ArithmeticError('strictly trapped a>0 required on the integration interval')
        G = ks_generator(rho, energies, mass, angular)
        radius, basis, r_ref = radius_state(rho)
        X = state[:2 * nsrc * nz].reshape(2, nsrc, nz)
        dX = apply_L(X, G, axial, radius, r_ref)
        if not evolve_tangents:
            return dX.ravel()
        Y = state[2 * nsrc * nz:].reshape(ndir, 2, nsrc, nz)
        dY = apply_L(Y, G, axial, radius, r_ref)
        if ndir:
            spin = np.einsum('ab,bsz->asz', S2, X)
            dY = dY - ((1j * angular / axial) * (basis / (radius * radius)))[:, None, None, :] * spin[None]
        return np.concatenate((dX.ravel(), dY.ravel()))

    # Keep the last accepted state, not an O(number-of-steps) stack of full
    # fields. This is the same DOP853 implementation used by solve_ivp.
    solution = DOP853(rhs, rho_up, y0, RHO_SIGMA, rtol=rtol, atol=atol,
                     max_step=max_step)
    accepted_steps = 0
    while solution.status == 'running':
        message = solution.step()
        if solution.status == 'failed':
            raise ArithmeticError(message)
        accepted_steps += 1
    final = solution.y
    X = final[:2 * nsrc * nz].reshape(2, nsrc, nz)
    Y = (final[2 * nsrc * nz:].reshape(ndir, 2, nsrc, nz)
         if evolve_tangents else np.zeros((0, 2, nsrc, nz), complex))
    X_t = np.moveaxis(trigonometric_polynomial(X, z, target), -1, 0)
    Xz_t = np.moveaxis(trigonometric_polynomial(X, z, target, derivative=True), -1, 0)
    if ndir:
        Y_t = np.moveaxis(trigonometric_polynomial(Y, z, target), -1, 1)
        Yz_t = np.moveaxis(trigonometric_polynomial(Y, z, target, derivative=True), -1, 1)
    else:
        Y_t = np.zeros((0, len(target), 2, nsrc), complex)
        Yz_t = np.zeros_like(Y_t)
    phase = np.exp(-1j * energies * target[:, None])
    columns = phase[:, None, :] * X_t
    axial = phase[:, None, :] * (Xz_t - 1j * energies * X_t)
    tangents_F = phase[None, :, None, :] * Y_t
    tangents_Fz = phase[None, :, None, :] * (Yz_t - 1j * energies * Y_t)
    initial_norm = float(np.linalg.norm(X0.ravel()))
    final_norm = float(np.linalg.norm(X.ravel()))
    preparation = _fixed_preparation_digest(
        source, initial, mass, angular, rho_up, RHO_SIGMA)
    binding = KSEnvelopeBinding(
        z, w_samples, U_samples, tuple(map(float, amplitudes)), inners, outers,
        support, rho_up, RHO_SIGMA, rtol, atol, max_step)
    diagnostics = {
        'function_evaluations': int(solution.nfev),
        'accepted_steps': accepted_steps,
        'rtol': rtol,
        'atol': atol,
        'max_step': max_step,
        'rho_up': rho_up,
        'rho_sigma': RHO_SIGMA,
        'continuum_speed_distance': distance,
        'padding': (float(padding_left), float(padding_right)),
        'computational_length': float(length),
        'computational_count': int(nz),
        'axial_support': list(support),
        'physical_interval': list(interval),
        'radius_lower_bound': bound if global_radius_bound else None,
        'sampled_radius_lower_estimate': bound,
        'global_radius_bound_supplied': global_radius_bound,
        'initial_envelope_norm': initial_norm,
        'final_envelope_norm': final_norm,
        'computational_norm_residual': abs(final_norm - initial_norm),
        'tangents': tangents,
        'fixed_preparation_digest': preparation,
        'periodic_envelope_is_physical_Pi_BC': False,
        'physical_confidence_from_tolerances': False,
        'PG_discretization_copied': False,
        'callback_identity_claimed': False,
        'source_frequency_controls_spatial_mesh': False,
        'continuum_error_bound': 'OPEN',
        'integration_error_bound': 'OPEN',
        'truncation_error_bound': 'OPEN',
        'physical_source_error': 'OPEN',
        'finite_fourier_causality': 'OPEN',
        'extra_r_z_term': False,
        'u_equals_r_sqrt_a_psi': True,
    }
    return KSEnvelopeIncoming(
        target, columns, tangents_F, axial, tangents_Fz, source.covariance,
        source.column_weights, source.energies, mass, angular, initial, rho_up,
        binding, preparation, _coverage(nsrc, len(np.unique(np.round(energies, 15)))),
        diagnostics)
