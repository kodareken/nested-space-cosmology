"""Expose the boundary fields used by the archived middle-panel producer.

These are the same order-16 massive Riccati fields, including its source-column
phases, rather than an inverse of attenuated packet observations.  They retain
the producer's vacuum-mode approximation: thermal/reflection and asymptotic
stress remainders are not established by packet agreement or order comparison.
No covariance, metric history, transport time, or horizon solve is introduced.
"""
from dataclasses import dataclass
import hashlib
import inspect
from pathlib import Path

import numpy as np

from .nsc_lorentzian import geometry
from .nsc_pg_fast_packets import FastVacuumPacketProjector
from .nsc_pg_high_energy import riccati_coefficients
from .nsc_transmitting_dirac_domain import I2, S2, MODE_TO_CURRENT


def _boundary_columns(energies, coefficients, order):
    powers = (1/(2*energies[:, None]))**np.arange(1, order+1)
    series = powers@coefficients[:order]
    beta = float(geometry(1.)[0])

    def branch_vector(branch):
        vector = (np.stack((np.ones_like(series), 1j*(1-beta)*series), axis=1)
                  if branch == 1 else
                  np.stack((-1j*(1+beta)*series.conj(), np.ones_like(series)), axis=1))
        current = (1-beta)*abs(vector[:, 0])**2+(-1-beta)*abs(vector[:, 1])**2
        if np.any(current >= 0):
            raise ArithmeticError('middle boundary branch lost its trapped current orientation')
        return (vector@MODE_TO_CURRENT.T)/np.sqrt(abs(current))[:, None]

    partner = branch_vector(1)
    incoming = branch_vector(-1)
    return np.stack((np.zeros_like(partner), partner, incoming), axis=2), series


def _equation_residual(energies, mass, angular, series, derivative, branch):
    """Vectorized evaluation of the existing branch_riccati_residual formula."""
    beta, bp, _ = map(float, geometry(1.))
    vp, vm = 1-beta, -1-beta
    uplus = angular/np.sqrt(2.)+1j*mass
    uminus = uplus.conjugate()
    if branch == 1:
        q = 1j*vp*series
        dq = 1j*(-bp*series+vp*derivative)
        rate = 1j*energies/vp+bp/(2*vp)-1j*uminus*series
        y = np.stack((np.ones_like(q), q), axis=1)
        dy = np.stack((rate, dq+q*rate), axis=1)
    else:
        series, derivative = series.conj(), derivative.conj()
        q = -1j*(1+beta)*series
        dq = -1j*(bp*series+(1+beta)*derivative)
        rate = 1j*energies/vm+bp/(2*vm)+1j*uplus*series
        y = np.stack((q, np.ones_like(q)), axis=1)
        dy = np.stack((dq+q*rate, rate), axis=1)
    coupling = np.stack((-1j*uminus*y[:, 1], 1j*uplus*y[:, 0]), axis=1)
    error = -1j*(np.array([vp, vm])*dy-.5*bp*y)+coupling-energies[:, None]*y
    return np.linalg.norm(error, axis=1)/np.linalg.norm(y, axis=1)


@dataclass
class ArchivedMiddleBoundaryModes:
    energies: np.ndarray
    mode_at_one: np.ndarray
    current_coisometry_residual: np.ndarray
    riccati_equation_residual: np.ndarray
    order12_16_field_difference: np.ndarray
    provenance: dict

    def compare_observations(self, observation, projection):
        """Apply saved fundamental overlaps; never rerun a packet integral."""
        A, F = np.asarray(observation, complex), np.asarray(projection, complex)
        n = len(self.energies)
        if A.shape != (n, 6, 2) or F.shape not in ((n, 6, 3), (n, 8, 3)):
            raise ValueError('saved real-frequency observation maps and source columns required')
        if not np.isfinite(A).all() or not np.isfinite(F).all():
            raise ValueError('finite archived observation data required')
        error = np.linalg.norm(A@self.mode_at_one-F[:, :6], axis=(1, 2))
        relative = error/np.maximum(np.linalg.norm(F[:, :6], axis=(1, 2)), np.finfo(float).tiny)
        return {'projection_absolute_residual': error, 'projection_relative_residual': relative}


def archived_middle_boundary_modes(energies, metadata, *, repo_root=None):
    """Recover the fields from the exact hash-pinned archived producer recipe.

    `metadata` is a mid or mid_ref component's decoded metadata_json.  Its
    sources must authenticate the current producer and profile definitions.
    The producer's actual default order is read, not inferred from copy counts
    or frequency weights.  Field columns have the unchanged order
    [zero exterior-horizon column, interior partner, incoming infinity].
    """
    raw = np.asarray(energies)
    if raw.ndim != 1 or not len(raw) or np.iscomplexobj(raw):
        raise ValueError('positive real middle-panel frequencies required')
    E = np.asarray(raw, float)
    root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[2]
    sources = metadata.get('sources', {})
    required = ('src/recursive_horizons/nsc_pg_fast_packets.py',
                'src/recursive_horizons/nsc_pg_high_energy.py',
                'src/recursive_horizons/nsc_lorentzian.py',
                'src/recursive_horizons/nsc_transmitting_dirac_domain.py')
    if not all(p in sources for p in required):
        raise ValueError('archived middle-producer source hashes required')
    for path, digest in sources.items():
        if hashlib.sha256((root/path).read_bytes()).hexdigest() != digest:
            raise ValueError('archived middle-producer dependency changed: '+path)
    order = inspect.signature(FastVacuumPacketProjector.__init__).parameters['order'].default
    if order != 16:
        raise ValueError('archived producer order differs from the retained order-16 convention')
    left, right = float(metadata['left']), float(metadata['right'])
    channel, sign = metadata['channel'], int(metadata['angular_sign'])
    mass = float(channel['compact_mass'])
    angular = sign*float(channel['angular_eigenvalue'])
    if (not np.isfinite(E).all() or np.any(E < left) or np.any(E > right)
            or left < 16 or right <= left or sign not in (-1, 1)
            or mass < 0 or not np.isfinite([mass, angular]).all()):
        raise ValueError('frequencies and channel must belong to the declared middle panel')
    coefficients, derivatives = riccati_coefficients([1.], mass, angular, order)
    fields, series = _boundary_columns(E, coefficients[0], order)
    coarse, _ = _boundary_columns(E, coefficients[0], 12)
    powers = (1/(2*E[:, None]))**np.arange(1, order+1)
    series_derivative = powers@derivatives[0]
    equation = np.stack([_equation_residual(E, mass, angular, series, series_derivative, branch)
                         for branch in (1, -1)], axis=1)
    beta = float(geometry(1.)[0])
    eig, vectors = np.linalg.eigh(beta*I2-S2)
    flux_frame = (vectors*np.sqrt(eig))@vectors.conj().T
    normalized = flux_frame@fields
    current = np.linalg.norm(normalized@normalized.swapaxes(-1, -2).conj()-I2, axis=(1, 2))
    provenance = {
        'source_hashes': dict(sources), 'order': order, 'comparison_order': 12,
        'producer': 'FastVacuumPacketProjector.project_many', 'rho': 1.,
        'column_order': ['zero', 'partner', 'incoming_one'],
        'boundary_phase': 'exact archived project_many convention; no added clock factor',
        'state_covariance_changed': False, 'packet_integrals_rerun': False,
        'horizon_scattering_rerun': False, 'history_selected': False,
        'scope': 'same middle-panel vacuum-mode approximation; thermal/reflection and local stress-tail bounds remain separate',
        'order_difference_is_rigorous_bound': False}
    return ArchivedMiddleBoundaryModes(E.copy(), fields, current, equation,
                                       np.linalg.norm(fields-coarse, axis=(1, 2)), provenance)
