"""Finite-energy operator interpolation with the original physical source retained.

The identity columns below form an operator basis. They never replace C_up.
Only the propagated operator is interpolated; A_up, C_src, input labels and
quadrature weights are supplied unchanged at reconstruction. Error status OPEN.
"""
from dataclasses import dataclass
from hashlib import sha256
from types import MappingProxyType

import numpy as np

from .nsc_evolved_incoming_state import FixedSourcePreparation
from .nsc_ks_difference_envelope import KSDifferenceIncoming, evolve_ks_difference_envelope
from .nsc_ks_source_envelope import (
    _coverage, _finite_complex, _finite_real, _fixed_preparation_digest,
)


def _interval(value):
    out = tuple(map(float, value))
    if len(out) != 2 or not np.isfinite(out).all() or out[0] >= out[1]:
        raise ValueError('finite positive-length energy interval required')
    return out


def chebyshev_energy_nodes(interval, degree):
    """Rounded Chebyshev roots; node rounding is an explicit OPEN error."""
    lo, hi = _interval(interval)
    if isinstance(degree, (bool, np.bool_)) or int(degree) != degree or degree < 1:
        raise ValueError('integer interpolation degree >= 1 required')
    n = int(degree) + 1
    center, half = (lo + hi) / 2, (hi - lo) / 2
    nodes = np.sort(center + half * np.cos(np.pi * (np.arange(n) + .5) / n))
    if not np.isfinite(nodes).all() or np.any(np.diff(nodes) <= 0):
        raise ValueError('distinct finite energy nodes required at working precision')
    return nodes


def interpolation_matrix(nodes, target, interval):
    """Polynomial through the actual rounded nodes, not ideal-node weights.

    Barycentric first formula weights are calculated from node differences.
    Arithmetic, node displacement and node-solve errors remain separate from
    the ideal Chebyshev interpolation remainder.
    """
    lo, hi = _interval(interval)
    nodes, target = np.asarray(nodes, float), np.asarray(target, float)
    if (nodes.ndim != 1 or len(nodes) < 2 or not np.isfinite(nodes).all()
            or np.any(np.diff(nodes) <= 0) or nodes[0] < lo or nodes[-1] > hi):
        raise ValueError('increasing finite interpolation nodes inside interval required')
    if target.ndim != 1 or not np.isfinite(target).all() or np.any((target < lo) | (target > hi)):
        raise ValueError('target energy outside the declared interpolation interval')
    # Scaling leaves the barycentric quotient unchanged and avoids powers of E.
    x = (nodes - lo) / (hi - lo)
    differences = x[:, None] - x[None, :]
    np.fill_diagonal(differences, 1.)
    weights = 1. / np.prod(differences, axis=1)
    weights /= np.max(abs(weights))
    delta = (target[:, None] - nodes[None, :]) / (hi - lo)
    exact = delta == 0
    safe = np.where(exact, 1., delta)
    terms = weights[None, :] / safe
    rows = np.any(exact, axis=1)
    terms[rows] = exact[rows].astype(float)
    matrix = terms / np.sum(terms, axis=1)[:, None]
    if not np.isfinite(matrix).all():
        raise ArithmeticError('nonfinite barycentric evaluation')
    return matrix


@dataclass(frozen=True)
class KSEnergyPropagator:
    """Sampled difference operator, applied to the unchanged original source."""
    interval: tuple
    energies: object
    z: object
    reference: object       # (energy, output spin, input spin)
    difference: object      # (energy, z, output spin, input spin)
    difference_z: object
    tangent: object         # (energy, direction, z, output spin, input spin)
    tangent_z: object
    mass: float
    angular: float
    rho_up: float
    binding: object
    diagnostics: object

    def __post_init__(self):
        interval = _interval(self.interval)
        energies = _finite_real(self.energies, 'operator energies')
        z = _finite_real(self.z, 'operator targets')
        if (energies.ndim != 1 or len(energies) < 2 or np.any(np.diff(energies) <= 0)
                or energies[0] < interval[0] or energies[-1] > interval[1]):
            raise ValueError('increasing operator energies inside declared interval required')
        if z.ndim != 1 or not len(z) or np.any(np.diff(z) <= 0):
            raise ValueError('increasing operator target z required')
        n, nz = len(energies), len(z)
        if not np.isfinite([self.mass, self.angular, self.rho_up]).all() or self.mass < 0:
            raise ValueError('finite operator parameters required')
        if self.rho_up != self.binding.rho_up:
            raise ValueError('operator rho_up differs from history binding')
        tangent = _finite_complex(self.tangent, 'operator tangent')
        if tangent.ndim != 5 or tangent.shape[0] != n or tangent.shape[2:] != (nz, 2, 2):
            raise ValueError('operator tangent shape must be (energy,direction,z,2,2)')
        shapes = {'reference': (n, 2, 2), 'difference': (n, nz, 2, 2),
                  'difference_z': (n, nz, 2, 2), 'tangent': tangent.shape,
                  'tangent_z': tangent.shape}
        for key, shape in shapes.items():
            array = _finite_complex(getattr(self, key), key)
            if array.shape != shape:
                raise ValueError(f'{key} shape must be {shape}')
            object.__setattr__(self, key, array)
        object.__setattr__(self, 'interval', interval)
        object.__setattr__(self, 'energies', energies)
        object.__setattr__(self, 'z', z)
        if self.diagnostics.get('reference_mode', 'interpolated') not in (
                'interpolated', 'direct-original-energies'):
            raise ValueError('unknown energy reference reconstruction')
        object.__setattr__(self, 'diagnostics', MappingProxyType(dict(self.diagnostics)))

    @property
    def digest(self):
        h = sha256(self.binding.fingerprint.encode())
        h.update(repr((self.interval, float(self.mass), float(self.angular), float(self.rho_up))).encode())
        if 'reference_mode' in self.diagnostics:
            h.update(str(self.diagnostics['reference_mode']).encode())
        for name in ('energies', 'z', 'reference', 'difference', 'difference_z', 'tangent', 'tangent_z'):
            value = np.asarray(getattr(self, name))
            h.update(repr(value.shape).encode())
            h.update(np.ascontiguousarray(value).tobytes())
        return h.hexdigest()

    def apply(self, source, initial_columns):
        """Reconstruct the original source, with an explicit unclosed error budget."""
        if not isinstance(source, FixedSourcePreparation):
            raise TypeError('original FixedSourcePreparation required')
        initial = _finite_complex(initial_columns, 'original A_up')
        if initial.shape != (2, len(source.energies)):
            raise ValueError('original A_up must match the source columns')
        matrix = interpolation_matrix(self.energies, source.energies, self.interval)
        def interpolate(value):
            return np.tensordot(matrix, value, axes=(1, 0))
        A = np.einsum('eab,be->ae', interpolate(self.reference), initial)
        D = np.einsum('ezab,be->zae', interpolate(self.difference), initial)
        Dz = np.einsum('ezab,be->zae', interpolate(self.difference_z), initial)
        Y = np.einsum('edzab,be->dzae', interpolate(self.tangent), initial)
        Yz = np.einsum('edzab,be->dzae', interpolate(self.tangent_z), initial)
        phase = np.exp(-1j * source.energies * self.z[:, None])[:, None, :]
        digest = _fixed_preparation_digest(source, initial, self.mass, self.angular,
                                           self.rho_up, self.binding.rho_sigma)
        mode = self.diagnostics.get('reference_mode', 'interpolated')
        diagnostics = {
            'field_representation': 'energy-interpolated joint reference plus difference',
            'operator_digest': self.digest, 'fixed_preparation_digest': digest,
            'interpolated_quantity': 'Dirac operator only; original A_up and C_src retained',
            'identity_is_operator_basis': True, 'source_interpolated': False,
            'energy_interval': self.interval, 'operator_node_count': len(self.energies),
            'interpolation_error_bound': None, 'node_evolution_error_bound': None,
            'node_rounding_error_bound': None, 'reconstruction_rounding_error_bound': None,
            'physical_constraint_status': 'OPEN', 'metric_timestep': False,
            'reference_mode': mode,
        }

        def prepared(reference, diagnostic):
            full = reference[None] + D
            return KSDifferenceIncoming(
                self.z, phase * full, phase[None] * Y,
                phase * (Dz - 1j * source.energies * full),
                phase[None] * (Yz - 1j * source.energies * Y),
                source.covariance, source.column_weights, source.energies, self.mass, self.angular,
                initial, self.rho_up, self.binding, digest,
                _coverage(len(source.energies), len(np.unique(source.energies))), diagnostic,
                reference, D, Dz)

        result = prepared(A, diagnostics)
        if mode == 'direct-original-energies':
            from .nsc_ks_local_constraints import homogeneous_reference_amplitudes
            # Exact identity: U_g = U_ref + (U_g-U_ref). D is the propagated
            # difference operator, not a fitted or measured stress offset.
            # The reference uses the original energies and initial columns.
            direct_reference = homogeneous_reference_amplitudes(result)
            diagnostics.update({
                'interpolated_quantity': 'difference operator and retarded tangent only',
                'reference_integrator': 'dop853',
                'reference_uses_original_energies': True,
                'fixed_reference_error_bound': None,
                'stress_drift_subtracted': False,
            })
            result = prepared(direct_reference, diagnostics)
        return result


_CONTROL_KEYS = (
    'function_evaluations', 'accepted_steps', 'time_integrator', 'step_control',
    'n_primal', 'tangent_directions', 'tangent_rtol', 'tangent_atol',
    'primal_error_norm', 'tangent_error_norm', 'stage_coefficient_evaluations',
    'analytic_operator_trace',
)


def _propagator_diagnostics(state):
    diagnostics = {
        'physical_source_supplied': False, 'identity_is_operator_basis': True,
        'physical_constraint_status': 'OPEN',
        'reference_mode': 'direct-original-energies',
    }
    source = dict(state.diagnostics)
    for key in _CONTROL_KEYS:
        if key in source:
            diagnostics[key] = source[key]
    return diagnostics


def evolve_energy_propagator(interval, degree, family, z_grid, target_z,
                             mass, angular, rho_up, **solver_options):
    """Use the existing Dirac owner on two basis columns at each energy.

    No matter contraction is performed on the identity bookkeeping matrix.
    Returned operator data have no physical source covariance.
    """
    # New production operator evaluations protect the primal block and each
    # tangent direction separately. Explicit joint control remains available
    # only as a diagnostic of historical numerics.
    solver_options.setdefault('step_control', 'primal')
    if (solver_options['step_control'] == 'primal'
            and solver_options.get('tangents', 'all') == 'all'):
        solver_options.setdefault('tangent_rtol', solver_options.get('rtol'))
        solver_options.setdefault('tangent_atol', solver_options.get('atol'))
    energies = chebyshev_energy_nodes(interval, degree)
    n = len(energies)
    basis = np.tile(np.eye(2, dtype=complex), (1, n))
    operator_labels = FixedSourcePreparation(np.eye(2*n), np.ones(2*n), np.repeat(energies, 2))
    state = evolve_ks_difference_envelope(
        operator_labels, basis, family, z_grid, target_z, mass, angular, rho_up, **solver_options)
    # Recover unweighted envelope tangents. These finite-precision conversions
    # belong to reconstruction rounding; no exact cancellation is asserted.
    phase = np.exp(1j * state.source_energies * state.z[:, None])[:, None, :]
    Y = phase[None] * state.column_tangents
    Yz = phase[None] * state.axial_tangents + 1j * state.source_energies * Y
    def fibers(value):
        value = np.asarray(value)
        return np.moveaxis(value.reshape(*value.shape[:-1], n, 2), -2, 0)
    return KSEnergyPropagator(
        tuple(interval), energies, state.z, fibers(state.reference_amplitudes),
        fibers(state.envelope_difference), fibers(state.envelope_difference_z),
        fibers(Y), fibers(Yz), mass, angular, rho_up, state.binding,
        _propagator_diagnostics(state))
