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

from pathlib import Path as _receipt_Path
from .nsc_frozen_sources import selected_source_sha256 as _receipt_source_sha256

# Retained executable control: comparison_control
# Folded from scripts/derive_nsc_ks_envelope_comparison.py at 329ef809; original receipts retain their producer identities.
import os as _comparison_control_os
import argparse as _comparison_control_argparse
from hashlib import sha256 as _comparison_control_sha256
import json as _comparison_control_json
from pathlib import Path as _comparison_control_Path
import signal as _comparison_control_signal
import sys as _comparison_control_sys
import time as _comparison_control_time
import numpy as _comparison_control_np
from recursive_horizons.nsc_posix import ROOT as _comparison_control_ROOT
_comparison_control_sys.path.insert(0, str(_comparison_control_ROOT / 'src'))
_comparison_control_sys.path.insert(0, str(_comparison_control_ROOT / 'scripts'))
from recursive_horizons.nsc_evolved_incoming_constraints import prepared_response_control as _comparison_control_P
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation as _comparison_control_FixedSourcePreparation
from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter as _comparison_control_source_column_matter
from recursive_horizons.nsc_ks_source_envelope import evolve_ks_source_envelope as _comparison_control_evolve_ks_source_envelope, computational_z_grid as _comparison_control_computational_z_grid, usual_axial_support as _comparison_control_usual_axial_support
from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates as _comparison_control_chart_coordinates
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes as _comparison_control_deterministic_npz_bytes
from recursive_horizons.nsc_retarded_radial_response import archived_amplitudes as _comparison_control_archived_amplitudes
_comparison_control_OUTPUT = 'results/development/nsc-ks-envelope-comparison.json'
_comparison_control_ARTIFACT = 'results/development/artifacts/nsc-ks-envelope-comparison.npz'
_comparison_control_OWNED = ('scripts/derive_nsc_ks_envelope_comparison.py', 'docs/nsc-ks-envelope-comparison.md')
_comparison_control_INPUTS = (*_comparison_control_P.INPUTS, 'scripts/derive_nsc_local_prepared_response.py', 'docs/nsc-local-prepared-response.md', 'results/development/nsc-local-cf4-response.json', 'src/recursive_horizons/nsc_ks_source_envelope.py')
_comparison_control_CAP = 120.0
_comparison_control_TARGET = 1e-11
_comparison_control_PG_AGREEMENT = 5e-09
_comparison_control_BASE_OPTIONS = {'rtol': 2e-12, 'atol': 2e-14, 'max_step': 0.001}
_comparison_control_TIGHT_OPTIONS = {'rtol': 2e-13, 'atol': 2e-15, 'max_step': 0.0005}

class _comparison_control_BudgetReached(RuntimeError):
    pass

def _comparison_control_digest(p):
    candidate = _receipt_Path(p)
    relative = candidate.relative_to(_comparison_control_ROOT).as_posix() if candidate.is_absolute() else candidate.as_posix()
    if relative == 'src/recursive_horizons/nsc_local_incoming_family.py':
        return _original_envelope_family_sha256(_comparison_control_ROOT)
    return _receipt_source_sha256(_comparison_control_ROOT, relative)

def _original_envelope_family_sha256(root):
    """Authenticate the receipt's explicitly pinned original family source.

    This older comparison/resolution recipe used the declared source-history
    pin, while its successor API lives in this module. Default calls still
    hash current bytes. Numerical output paths never enter this adapter.
    """
    from .nsc_frozen_sources import authenticate_source, selected_source_context
    path = 'src/recursive_horizons/nsc_local_incoming_family.py'
    if selected_source_context() is None:
        return _receipt_source_sha256(root, path)
    context = '09f3faa364ca88b8ff482c413022ce38f40bef0d'
    expected = '9d2b9fc57e2057d6c34d8943384428341b314370021bcd4156af88320c3de04b'
    manifest = _comparison_control_json.loads((root/'archive/source-history/manifest.json').read_text(encoding='utf-8'))
    if not any(row.get('path') == path and row.get('commit') == context
               and row.get('sha256') == expected for row in manifest['pins']):
        raise ValueError('original envelope family source pin changed')
    authenticate_source(root, path, expected, context=context)
    return expected


def _envelope_successor_source_hashes():
    from .nsc_frozen_sources import source_sha256
    root = _comparison_control_ROOT
    return {path: source_sha256(root, path) for path in (
        'src/recursive_horizons/nsc_ks_source_envelope.py',
        'src/recursive_horizons/nsc_evolved_incoming_constraints.py',
        'docs/nsc-ks-source-envelope.md')}

def _comparison_control_signatures():
    return {p: _comparison_control_digest(p) for p in (*_comparison_control_OWNED, *_comparison_control_INPUTS)}

def _comparison_control_maximum(x):
    return float(_comparison_control_np.max(abs(x), initial=0.0))

def _comparison_control_matter(arrays, name, parameters):
    F, Fz, dF, dFz = (arrays[name + '/' + k] for k in ('F', 'Fz', 'dF', 'dFz'))
    w = arrays['source/column_weights']
    C = arrays['source/covariance']
    return _comparison_control_source_column_matter(F * w, Fz * w, C, dF * w, dFz * w, **parameters)

def _comparison_control_analyze(arrays, meta):
    rows = {}
    values = {}
    checks = {}
    for name in meta['completed']:
        value = _comparison_control_matter(arrays, name, meta['parameters'])
        values[name] = value
        rows[name] = {'matter': value['action_gradient'].tolist(), 'matter_tangent': value['action_gradient_tangent'].tolist(), 'CPU_seconds': meta['cases'][name]['CPU_seconds'], 'function_evaluations': meta['cases'][name]['function_evaluations'], 'accepted_steps': meta['cases'][name]['accepted_steps'], 'computational_norm_residual': meta['cases'][name]['computational_norm_residual']}
    ref = _comparison_control_matter(arrays, 'analytic_reference', meta['parameters'])['action_gradient']
    if 'reference64' in values:
        drift = _comparison_control_np.max(abs(values['reference64']['action_gradient'] - ref), axis=0)
        checks['reference_matter_drift'] = drift.tolist()
        checks['reference_passed'] = bool(_comparison_control_np.max(drift) <= 3e-11)
    for coarse, fine in (('n64', 'n128'), ('n128', 'n256')):
        if coarse in values and fine in values:
            change = _comparison_control_np.max(abs(values[fine]['action_gradient'] - values[coarse]['action_gradient']), axis=0)
            checks[coarse + '_' + fine] = {'matter_change': change.tolist(), 'passes_indicator': bool(_comparison_control_np.max(change) <= _comparison_control_TARGET)}
    if 'tight' in values:
        base = meta['tight_base']
        change = _comparison_control_np.max(abs(values['tight']['action_gradient'] - values[base]['action_gradient']), axis=0)
        checks['time_tolerance_change'] = change.tolist()
        checks['time_indicator_passed'] = bool(_comparison_control_np.max(change) <= _comparison_control_TARGET)
        pg = _comparison_control_matter(arrays, 'pg_cf4_240', meta['parameters'])
        delta = _comparison_control_np.max(abs(values['tight']['action_gradient'] - pg['action_gradient']), axis=0)
        checks['PG_CF4_matter_difference'] = delta.tolist()
        checks['PG_agreement_indicator_passed'] = bool(_comparison_control_np.max(delta) <= _comparison_control_PG_AGREEMENT)
        checks['PG_field_difference'] = _comparison_control_maximum(arrays['tight/F'] - arrays['pg_cf4_240/F'])
        checks['PG_axial_difference'] = _comparison_control_maximum(arrays['tight/Fz'] - arrays['pg_cf4_240/Fz'])
        checks['PG_tangent_difference'] = _comparison_control_maximum(arrays['tight/dF'] - arrays['pg_cf4_240/dF'])
        checks['PG_axial_tangent_difference'] = _comparison_control_maximum(arrays['tight/dFz'] - arrays['pg_cf4_240/dFz'])
        checks['local_matter_correction_max'] = _comparison_control_np.max(abs(values['tight']['action_gradient'] - ref), axis=0).tolist()
    checks['same_fixed_preparation'] = len({v['preparation'] for v in meta['cases'].values()}) <= 1
    return {'per_case': rows, 'checks': checks}

def _comparison_control_record(arrays, meta, payload):
    return {'schema': 'NSC-KS-ENVELOPE-COMPARISON-v1', 'accountable_author': 'Douglas Ek', 'status': 'OPEN: independent prepared-state method comparison; physical local gate uncertified', 'control': {'history_amplitude': _comparison_control_P.ALPHA, 'source_family': '14_1', 'rho_up': meta['rho_up'], 'rho_sigma': 1.0, 'local_interval': 'S(1)+[.12,.18]', 'computational_length': 0.4, 'space_indicator_target': _comparison_control_TARGET, 'time_indicator_target': _comparison_control_TARGET, 'reference_target': 3e-11, 'PG_agreement_indicator': _comparison_control_PG_AGREEMENT, 'base_options': _comparison_control_BASE_OPTIONS, 'tight_options': _comparison_control_TIGHT_OPTIONS}, 'measurements': _comparison_control_analyze(arrays, meta), 'runtime': meta['runtime'], 'scope': {'physical_local_gate': 'OPEN', 'full_source_error': None, 'continuum_error': None, 'between_node_error': None, 'indicators_are_error_bounds': False, 'new_PG_runs': 0, 'new_source_preparation': False, 'source_changed': False, 'periodic_Pi_boundary_condition': False, 'source_energy_used_as_outgoing_momentum': False, 'stress_drift_subtracted': False, 'metric_timestep': False}, 'source_hashes': {p: meta['signature'][p] for p in _comparison_control_OWNED}, 'input_hashes': {p: meta['signature'][p] for p in _comparison_control_INPUTS}, 'reused_artifacts': meta['reused_artifacts'], 'payload': payload, 'reproducer': 'python3 scripts/derive_nsc_ks_envelope_comparison.py --check'}

def _comparison_control_publish(arrays, meta):
    raw = _comparison_control_deterministic_npz_bytes({**arrays, 'metadata_json': _comparison_control_np.frombuffer(_comparison_control_json.dumps(meta, sort_keys=True).encode(), _comparison_control_np.uint8)})
    (_comparison_control_ROOT / _comparison_control_ARTIFACT).write_bytes(raw)
    result = _comparison_control_record(arrays, meta, {'path': _comparison_control_ARTIFACT, 'sha256': _comparison_control_sha256(raw).hexdigest(), 'bytes': len(raw)})
    (_comparison_control_ROOT / _comparison_control_OUTPUT).write_text(_comparison_control_json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n')
    return result

def _comparison_control_inputs():
    records, a, channel = _comparison_control_P.load_inputs()
    pg = _comparison_control_json.loads((_comparison_control_ROOT / 'results/development/nsc-local-cf4-response.json').read_text())
    if _comparison_control_digest(pg['payload']['path']) != pg['payload']['sha256']:
        raise ValueError('stored PG payload changed')
    with _comparison_control_np.load(_comparison_control_ROOT / pg['payload']['path'], allow_pickle=False) as f:
        old = {k: f[k] for k in f.files}
    return (records, a, channel, pg, old)

def _comparison_control_run():
    if (_comparison_control_ROOT / _comparison_control_OUTPUT).exists() or (_comparison_control_ROOT / _comparison_control_ARTIFACT).exists():
        raise FileExistsError('existing comparison; use --check')
    sig = _comparison_control_signatures()
    records, a, channel, pg, old = _comparison_control_inputs()
    index = int(_comparison_control_np.flatnonzero(a['x'] >= 1.03)[0])
    rho_up = float(a['x'][index])
    initial = _comparison_control_archived_amplitudes(a['x'], a['initial_fields'], a['source/energies'], index).transpose(1, 0, 2).reshape(2, 12)
    source = _comparison_control_FixedSourcePreparation(a['source/covariance'], a['source/column_weights'], _comparison_control_np.repeat(a['source/energies'], 3))
    target = a['times'][40:] + _comparison_control_chart_coordinates(1.0)[1]
    arrays = {'target_z': target, 'initial_canonical_columns': initial, 'source/covariance': source.covariance, 'source/column_weights': source.column_weights, 'source/energies': source.energies}
    for key, oldkey in (('F', 'reference_columns'), ('Fz', 'reference_axial_columns')):
        arrays['analytic_reference/' + key] = a[oldkey][40:]
    arrays['analytic_reference/dF'] = arrays['analytic_reference/dFz'] = _comparison_control_np.zeros((0, len(target), 2, 12), complex)
    for key in ('F', 'Fz'):
        arrays['pg_cf4_240/' + key] = old['cf4_240/' + key][160::4]
    for key in ('dF', 'dFz'):
        arrays['pg_cf4_240/' + key] = old['cf4_240/' + key][:, 160::4]
    coeff = _comparison_control_P.coefficients_from_records(records[2], records[3])
    meta = {'signature': sig, 'completed': [], 'attempted': [], 'cases': {}, 'rho_up': rho_up, 'parameters': {'mass': channel['compact_mass'], 'angular': channel['angular_eigenvalue'], 'axial_scale': coeff['a'], 'radius': coeff['r'], 'multiplicity': channel['copy_count'] * channel['degeneracy'] / 2}, 'reused_artifacts': [records[0]['payload'], pg['payload']]}
    cpu = _comparison_control_time.process_time()
    wall = _comparison_control_time.monotonic()
    cap = False

    def stop(*_):
        raise _comparison_control_BudgetReached('KS comparison120 CPU-second cap')

    def checkpoint():
        meta['runtime'] = {'CPU_seconds': _comparison_control_time.process_time() - cpu, 'wall_seconds': _comparison_control_time.monotonic() - wall, 'CPU_cap': _comparison_control_CAP, 'attempted': list(meta['attempted']), 'completed': list(meta['completed']), 'cap_reached': cap}
        if _comparison_control_signatures() != sig:
            raise ValueError('comparison owner/input changed')
        return _comparison_control_publish(arrays, meta)

    def case(name, n, alpha, options, tangents):
        begin = _comparison_control_time.process_time()
        meta['attempted'].append(name)
        result = _comparison_control_evolve_ks_source_envelope(source, initial, _comparison_control_P.family(alpha), _comparison_control_computational_z_grid(n), target, channel['compact_mass'], channel['angular_eigenvalue'], rho_up, axial_support=_comparison_control_usual_axial_support(), tangents=tangents, **options)
        for key, value in (('F', result.columns), ('Fz', result.axial_columns), ('dF', result.column_tangents), ('dFz', result.axial_tangents)):
            arrays[name + '/' + key] = value
        meta['cases'][name] = {'CPU_seconds': _comparison_control_time.process_time() - begin, 'preparation': result.fixed_preparation_digest, 'binding': result.binding.fingerprint, **{k: result.diagnostics[k] for k in ('function_evaluations', 'accepted_steps', 'computational_norm_residual', 'continuum_speed_distance', 'padding', 'global_radius_bound_supplied')}}
        meta['completed'].append(name)
        checkpoint()
        print(name + ' complete CPU ' + str(meta['cases'][name]['CPU_seconds']), flush=True)
    previous = _comparison_control_signal.signal(_comparison_control_signal.SIGPROF, stop)
    _comparison_control_signal.setitimer(_comparison_control_signal.ITIMER_PROF, _comparison_control_CAP)
    try:
        case('reference64', 64, 0.0, _comparison_control_BASE_OPTIONS, 'zero')
        case('n64', 64, _comparison_control_P.ALPHA, _comparison_control_BASE_OPTIONS, 'all')
        case('n128', 128, _comparison_control_P.ALPHA, _comparison_control_BASE_OPTIONS, 'all')
        if not _comparison_control_analyze(arrays, meta)['checks']['n64_n128']['passes_indicator']:
            case('n256', 256, _comparison_control_P.ALPHA, _comparison_control_BASE_OPTIONS, 'all')
            base = 'n256'
            n = 256
        else:
            base = 'n128'
            n = 128
        meta['tight_base'] = base
        case('tight', n, _comparison_control_P.ALPHA, _comparison_control_TIGHT_OPTIONS, 'all')
    except _comparison_control_BudgetReached:
        cap = True
    finally:
        _comparison_control_signal.setitimer(_comparison_control_signal.ITIMER_PROF, 0)
        _comparison_control_signal.signal(_comparison_control_signal.SIGPROF, previous)
    return checkpoint()

def _comparison_control_check():
    saved = _comparison_control_json.loads((_comparison_control_ROOT / _comparison_control_OUTPUT).read_text())
    _comparison_control_inputs()
    if _comparison_control_digest(_comparison_control_ARTIFACT) != saved['payload']['sha256']:
        raise ValueError('saved envelope artifact changed')
    with _comparison_control_np.load(_comparison_control_ROOT / _comparison_control_ARTIFACT, allow_pickle=False) as a:
        arrays = {k: a[k] for k in a.files}
    meta = _comparison_control_json.loads(arrays['metadata_json'].tobytes())
    if meta['signature'] != _comparison_control_signatures():
        raise ValueError('envelope comparison owner/input changed')
    from .nsc_ks_matter_receipt_replay import require_envelope_record
    require_envelope_record(saved, _comparison_control_record(arrays, meta, saved['payload']), arrays, meta)
    return saved

from types import SimpleNamespace as _comparison_control_Namespace
comparison_control = _comparison_control_Namespace(
    ARTIFACT=_comparison_control_ARTIFACT,
    BASE_OPTIONS=_comparison_control_BASE_OPTIONS,
    BudgetReached=_comparison_control_BudgetReached,
    CAP=_comparison_control_CAP,
    FixedSourcePreparation=_comparison_control_FixedSourcePreparation,
    INPUTS=_comparison_control_INPUTS,
    OUTPUT=_comparison_control_OUTPUT,
    OWNED=_comparison_control_OWNED,
    successor_source_hashes=_envelope_successor_source_hashes,
    P=_comparison_control_P,
    PG_AGREEMENT=_comparison_control_PG_AGREEMENT,
    Path=_comparison_control_Path,
    ROOT=_comparison_control_ROOT,
    TARGET=_comparison_control_TARGET,
    TIGHT_OPTIONS=_comparison_control_TIGHT_OPTIONS,
    analyze=_comparison_control_analyze,
    archived_amplitudes=_comparison_control_archived_amplitudes,
    argparse=_comparison_control_argparse,
    chart_coordinates=_comparison_control_chart_coordinates,
    check=_comparison_control_check,
    computational_z_grid=_comparison_control_computational_z_grid,
    deterministic_npz_bytes=_comparison_control_deterministic_npz_bytes,
    digest=_comparison_control_digest,
    evolve_ks_source_envelope=_comparison_control_evolve_ks_source_envelope,
    inputs=_comparison_control_inputs,
    json=_comparison_control_json,
    matter=_comparison_control_matter,
    maximum=_comparison_control_maximum,
    np=_comparison_control_np,
    os=_comparison_control_os,
    publish=_comparison_control_publish,
    record=_comparison_control_record,
    run=_comparison_control_run,
    sha256=_comparison_control_sha256,
    signal=_comparison_control_signal,
    signatures=_comparison_control_signatures,
    source_column_matter=_comparison_control_source_column_matter,
    sys=_comparison_control_sys,
    time=_comparison_control_time,
    usual_axial_support=_comparison_control_usual_axial_support,
)

# Retained executable control: resolution_control
# Folded from scripts/derive_nsc_ks_envelope_resolution.py at 329ef809; original receipts retain their producer identities.
import os as _envelope_resolution_control_os
import argparse as _envelope_resolution_control_argparse
from hashlib import sha256 as _envelope_resolution_control_sha256
import json as _envelope_resolution_control_json
from pathlib import Path as _envelope_resolution_control_Path
import signal as _envelope_resolution_control_signal
import sys as _envelope_resolution_control_sys
import time as _envelope_resolution_control_time
import numpy as _envelope_resolution_control_np
from recursive_horizons.nsc_posix import ROOT as _envelope_resolution_control_ROOT
_envelope_resolution_control_sys.path.insert(0, str(_envelope_resolution_control_ROOT / 'src'))
_envelope_resolution_control_sys.path.insert(0, str(_envelope_resolution_control_ROOT / 'scripts'))
from recursive_horizons.nsc_ks_source_envelope import comparison_control as _envelope_resolution_control_C
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation as _envelope_resolution_control_FixedSourcePreparation
from recursive_horizons.nsc_ks_source_envelope import evolve_ks_source_envelope as _envelope_resolution_control_evolve_ks_source_envelope, computational_z_grid as _envelope_resolution_control_computational_z_grid, usual_axial_support as _envelope_resolution_control_usual_axial_support
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes as _envelope_resolution_control_deterministic_npz_bytes
_envelope_resolution_control_OUTPUT = 'results/development/nsc-ks-envelope-resolution.json'
_envelope_resolution_control_ARTIFACT = 'results/development/artifacts/nsc-ks-envelope-resolution.npz'
_envelope_resolution_control_OWNED = ('scripts/derive_nsc_ks_envelope_resolution.py', 'docs/nsc-ks-envelope-resolution.md')
_envelope_resolution_control_INPUTS = (*_envelope_resolution_control_C.OWNED, *_envelope_resolution_control_C.INPUTS, _envelope_resolution_control_C.OUTPUT)
_envelope_resolution_control_TIGHTER = {'rtol': 5e-14, 'atol': 5e-16, 'max_step': 0.00025}
_envelope_resolution_control_CAP = 30.0

class _envelope_resolution_control_BudgetReached(RuntimeError):
    pass

def _envelope_resolution_control_digest(p):
    candidate = _receipt_Path(p)
    relative = candidate.relative_to(_envelope_resolution_control_ROOT).as_posix() if candidate.is_absolute() else candidate.as_posix()
    if relative == 'src/recursive_horizons/nsc_local_incoming_family.py':
        return _original_envelope_family_sha256(_envelope_resolution_control_ROOT)
    return _receipt_source_sha256(_envelope_resolution_control_ROOT, relative)

def _envelope_resolution_control_signatures():
    return {p: _envelope_resolution_control_digest(p) for p in (*_envelope_resolution_control_OWNED, *_envelope_resolution_control_INPUTS)}

def _envelope_resolution_control_analyze(arrays, meta):
    values = {name: _envelope_resolution_control_C.matter(arrays, name, meta['parameters']) for name in ['previous256', *meta['completed']]}
    checks = {}
    for first, second, label in [('previous256', 'n512', 'space'), ('n512', 'tighter512', 'time')]:
        if second in values and first in values:
            delta = _envelope_resolution_control_np.max(abs(values[second]['action_gradient'] - values[first]['action_gradient']), axis=0)
            tangent = _envelope_resolution_control_np.max(abs(values[second]['action_gradient_tangent'] - values[first]['action_gradient_tangent']), axis=(0, 1))
            checks[label] = {'matter_change': delta.tolist(), 'tangent_change': tangent.tolist(), 'target': 1e-11, 'passes_indicator': bool(_envelope_resolution_control_np.max(delta) <= 1e-11)}
    return checks

def _envelope_resolution_control_record(arrays, meta, payload):
    checks = _envelope_resolution_control_analyze(arrays, meta)
    passed = all((k in checks and checks[k]['passes_indicator'] for k in ('space', 'time')))
    return {'schema': 'NSC-KS-ENVELOPE-RESOLUTION-v1', 'accountable_author': 'Douglas Ek', 'status': ('PASS' if passed else 'OPEN') + ': local envelope resolution indicators; physical gate OPEN', 'control': {'points': 512, 'max_runs': 2, 'history_amplitude': _envelope_resolution_control_C.P.ALPHA, 'local_interval': 'S(1)+[.12,.18]', 'options': [_envelope_resolution_control_C.TIGHT_OPTIONS, _envelope_resolution_control_TIGHTER]}, 'measurements': checks, 'runtime': meta['runtime'], 'scope': {'physical_local_gate': 'OPEN', 'full_source_error': None, 'continuum_error': None, 'between_node_error': None, 'indicators_are_error_bounds': False, 'source_changed': False, 'new_PG_runs': 0, 'new_source_preparation': False, 'stress_drift_subtracted': False, 'metric_timestep': False}, 'source_hashes': {p: meta['signature'][p] for p in _envelope_resolution_control_OWNED}, 'input_hashes': {p: meta['signature'][p] for p in _envelope_resolution_control_INPUTS}, 'reused_artifacts': meta['reused_artifacts'], 'payload': payload, 'reproducer': 'python3 scripts/derive_nsc_ks_envelope_resolution.py --check'}

def _envelope_resolution_control_publish(arrays, meta):
    raw = _envelope_resolution_control_deterministic_npz_bytes({**arrays, 'metadata_json': _envelope_resolution_control_np.frombuffer(_envelope_resolution_control_json.dumps(meta, sort_keys=True).encode(), _envelope_resolution_control_np.uint8)})
    (_envelope_resolution_control_ROOT / _envelope_resolution_control_ARTIFACT).write_bytes(raw)
    result = _envelope_resolution_control_record(arrays, meta, {'path': _envelope_resolution_control_ARTIFACT, 'sha256': _envelope_resolution_control_sha256(raw).hexdigest(), 'bytes': len(raw)})
    (_envelope_resolution_control_ROOT / _envelope_resolution_control_OUTPUT).write_text(_envelope_resolution_control_json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n')
    return result

def _envelope_resolution_control_run():
    if (_envelope_resolution_control_ROOT / _envelope_resolution_control_OUTPUT).exists() or (_envelope_resolution_control_ROOT / _envelope_resolution_control_ARTIFACT).exists():
        raise FileExistsError('existing resolution record; use --check')
    sig = _envelope_resolution_control_signatures()
    old = _envelope_resolution_control_C.check()
    with _envelope_resolution_control_np.load(_envelope_resolution_control_ROOT / old['payload']['path'], allow_pickle=False) as f:
        previous = {k: f[k] for k in f.files}
    oldmeta = _envelope_resolution_control_json.loads(previous['metadata_json'].tobytes())
    arrays = {k: v for k, v in previous.items() if k.startswith('source/') or k in ('target_z', 'initial_canonical_columns')}
    for key in ('F', 'Fz', 'dF', 'dFz'):
        arrays['previous256/' + key] = previous['tight/' + key]
    source = _envelope_resolution_control_FixedSourcePreparation(arrays['source/covariance'], arrays['source/column_weights'], arrays['source/energies'])
    meta = {'signature': sig, 'parameters': oldmeta['parameters'], 'rho_up': oldmeta['rho_up'], 'completed': [], 'attempted': [], 'cases': {}, 'reused_artifacts': [old['payload']]}
    cpu = _envelope_resolution_control_time.process_time()
    wall = _envelope_resolution_control_time.monotonic()
    cap = False

    def stop(*_):
        raise _envelope_resolution_control_BudgetReached('two512-node checks30 CPU-second cap')

    def checkpoint():
        meta['runtime'] = {'CPU_seconds': _envelope_resolution_control_time.process_time() - cpu, 'wall_seconds': _envelope_resolution_control_time.monotonic() - wall, 'CPU_cap': _envelope_resolution_control_CAP, 'completed': list(meta['completed']), 'attempted': list(meta['attempted']), 'cap_reached': cap}
        if _envelope_resolution_control_signatures() != sig:
            raise ValueError('resolution owner/input changed')
        return _envelope_resolution_control_publish(arrays, meta)
    prior = _envelope_resolution_control_signal.signal(_envelope_resolution_control_signal.SIGPROF, stop)
    _envelope_resolution_control_signal.setitimer(_envelope_resolution_control_signal.ITIMER_PROF, _envelope_resolution_control_CAP)
    try:
        for name, options in (('n512', _envelope_resolution_control_C.TIGHT_OPTIONS), ('tighter512', _envelope_resolution_control_TIGHTER)):
            begin = _envelope_resolution_control_time.process_time()
            meta['attempted'].append(name)
            result = _envelope_resolution_control_evolve_ks_source_envelope(source, arrays['initial_canonical_columns'], _envelope_resolution_control_C.P.family(_envelope_resolution_control_C.P.ALPHA), _envelope_resolution_control_computational_z_grid(512), arrays['target_z'], meta['parameters']['mass'], meta['parameters']['angular'], meta['rho_up'], axial_support=_envelope_resolution_control_usual_axial_support(), tangents='all', **options)
            for key, value in (('F', result.columns), ('Fz', result.axial_columns), ('dF', result.column_tangents), ('dFz', result.axial_tangents)):
                arrays[name + '/' + key] = value
            if result.fixed_preparation_digest != oldmeta['cases']['tight']['preparation']:
                raise ValueError('fixed upstream preparation changed')
            meta['cases'][name] = {'CPU_seconds': _envelope_resolution_control_time.process_time() - begin, 'preparation': result.fixed_preparation_digest, 'binding': result.binding.fingerprint, 'function_evaluations': result.diagnostics['function_evaluations']}
            meta['completed'].append(name)
            checkpoint()
    except _envelope_resolution_control_BudgetReached:
        cap = True
    finally:
        _envelope_resolution_control_signal.setitimer(_envelope_resolution_control_signal.ITIMER_PROF, 0)
        _envelope_resolution_control_signal.signal(_envelope_resolution_control_signal.SIGPROF, prior)
    return checkpoint()

def _envelope_resolution_control_check():
    saved = _envelope_resolution_control_json.loads((_envelope_resolution_control_ROOT / _envelope_resolution_control_OUTPUT).read_text())
    _envelope_resolution_control_C.check()
    if _envelope_resolution_control_digest(_envelope_resolution_control_ARTIFACT) != saved['payload']['sha256']:
        raise ValueError('resolution payload changed')
    with _envelope_resolution_control_np.load(_envelope_resolution_control_ROOT / _envelope_resolution_control_ARTIFACT, allow_pickle=False) as f:
        arrays = {k: f[k] for k in f.files}
    meta = _envelope_resolution_control_json.loads(arrays['metadata_json'].tobytes())
    if meta['signature'] != _envelope_resolution_control_signatures():
        raise ValueError('resolution source/input changed')
    from .nsc_ks_matter_receipt_replay import require_envelope_record
    require_envelope_record(saved, _envelope_resolution_control_record(arrays, meta, saved['payload']), arrays, meta)
    return saved

from types import SimpleNamespace as _envelope_resolution_control_Namespace
resolution_control = _envelope_resolution_control_Namespace(
    ARTIFACT=_envelope_resolution_control_ARTIFACT,
    BudgetReached=_envelope_resolution_control_BudgetReached,
    C=_envelope_resolution_control_C,
    CAP=_envelope_resolution_control_CAP,
    FixedSourcePreparation=_envelope_resolution_control_FixedSourcePreparation,
    INPUTS=_envelope_resolution_control_INPUTS,
    OUTPUT=_envelope_resolution_control_OUTPUT,
    OWNED=_envelope_resolution_control_OWNED,
    successor_source_hashes=_envelope_successor_source_hashes,
    Path=_envelope_resolution_control_Path,
    ROOT=_envelope_resolution_control_ROOT,
    TIGHTER=_envelope_resolution_control_TIGHTER,
    analyze=_envelope_resolution_control_analyze,
    argparse=_envelope_resolution_control_argparse,
    check=_envelope_resolution_control_check,
    computational_z_grid=_envelope_resolution_control_computational_z_grid,
    deterministic_npz_bytes=_envelope_resolution_control_deterministic_npz_bytes,
    digest=_envelope_resolution_control_digest,
    evolve_ks_source_envelope=_envelope_resolution_control_evolve_ks_source_envelope,
    json=_envelope_resolution_control_json,
    np=_envelope_resolution_control_np,
    os=_envelope_resolution_control_os,
    publish=_envelope_resolution_control_publish,
    record=_envelope_resolution_control_record,
    run=_envelope_resolution_control_run,
    sha256=_envelope_resolution_control_sha256,
    signal=_envelope_resolution_control_signal,
    signatures=_envelope_resolution_control_signatures,
    sys=_envelope_resolution_control_sys,
    time=_envelope_resolution_control_time,
    usual_axial_support=_envelope_resolution_control_usual_axial_support,
)
