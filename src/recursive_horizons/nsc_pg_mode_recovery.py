"""Recover stationary interior mode columns from existing spectral probes.

At each real frequency the already-owned radial Dirac equation has two
independent interior solutions.  Six packet overlaps determine their two
coefficients when the actual observation map is sufficiently resolved.  This
does not replace the continuum field by a closed system of packet variables.
No covariance, horizon solve, state completion or spectral integral is input.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from threadpoolctl import threadpool_limits

from .nsc_lorentzian import geometry
from .nsc_pg_lll_preparation import radial_packet
from .nsc_pg_retarded_packets import _parts
from .nsc_transmitting_dirac_domain import I2, S2


@dataclass
class RecoveredInteriorModes:
    energies: np.ndarray
    mode_at_one: np.ndarray
    accepted: np.ndarray
    diagnostics: list

    def require_all(self):
        if not self.accepted.all():
            failed = np.flatnonzero(~self.accepted).tolist()
            raise ValueError(f"unresolved stationary mode recovery at indices {failed}")
        return self.mode_at_one


@dataclass
class StationaryInteriorObservation:
    energies: np.ndarray
    mass: float
    angular: float
    observation: np.ndarray
    fundamental_at_zero: np.ndarray
    fundamental_at_minus_one: np.ndarray
    current_residual: np.ndarray
    function_evaluations: int

    def recover(self, projection, projector, *, input_absolute_error=3e-11,
                field_absolute_tolerance=3e-8, projection_tolerance=3e-11,
                current_tolerance=3e-8, require_all=True):
        """Left-invert only sufficiently resolved physical mode observations.

        The error amplification is input_absolute_error / sigma_min.  An
        unresolved singular value is rejected, never floored or regularized.
        `require_all=False` returns explicit rejected entries with NaN fields.
        Input-error propagation is a numerical sensitivity diagnostic, not a
        theorem that upstream error estimates are rigorous bounds.
        """
        F = np.asarray(projection, complex)
        P = np.asarray(projector, complex)
        n = len(self.energies)
        if F.shape not in ((n, 6, 3), (n, 8, 3)):
            raise ValueError("frequency-resolved six/eight-packet source columns required; covariance is not input")
        if P.shape != (n, 3, 3) or not np.isfinite(F).all() or not np.isfinite(P).all():
            raise ValueError("finite physical source columns and open-fiber projectors required")
        if min(input_absolute_error, field_absolute_tolerance,
               projection_tolerance, current_tolerance) <= 0:
            raise ValueError("positive declared numerical tolerances required")
        if max(np.max(abs(P-P.swapaxes(-1, -2).conj())),
               np.max(abs(P@P-P))) > 3e-11:
            raise ValueError("canonical open-fiber projector required")
        F = F[:, :6]
        left, singular, right_h = np.linalg.svd(self.observation, full_matrices=False)
        beta = float(geometry(1.)[0])
        velocity = S2-beta*I2
        eig, vectors = np.linalg.eigh(-velocity)
        flux_frame = (vectors*np.sqrt(eig))@vectors.conj().T
        fields = np.full((n, 2, 3), np.nan+1j*np.nan, complex)
        accepted = np.zeros(n, bool)
        diagnostics = []
        for i in range(n):
            largest, smallest = map(float, singular[i])
            condition = largest/smallest if smallest else float('inf')
            amplification = input_absolute_error/smallest if smallest else float('inf')
            resolved = bool(smallest > np.finfo(float).eps*6*largest
                            and amplification <= field_absolute_tolerance)
            row = {'energy': float(self.energies[i]),
                   'singular_values': singular[i].tolist(),
                   'condition_number': condition,
                   'input_error_amplification': amplification,
                   'input_absolute_error': float(input_absolute_error),
                   'field_absolute_tolerance': float(field_absolute_tolerance),
                   'transport_current_residual': float(self.current_residual[i]),
                   'projection_absolute_residual': None,
                   'projection_relative_residual': None,
                   'current_coisometry_residual': None,
                   'source_projector_residual': None}
            if resolved:
                recovered = right_h[i].conj().T @ (
                    (left[i].conj().T@F[i])/singular[i, :, None])
                absolute = float(np.linalg.norm(self.observation[i]@recovered-F[i]))
                relative = absolute/max(float(np.linalg.norm(F[i])), np.finfo(float).tiny)
                normalized = flux_frame@recovered
                current = float(np.linalg.norm(normalized@P[i]@normalized.conj().T-I2, 2))
                source = float(np.linalg.norm(recovered@(np.eye(3)-P[i]), 2))
                row.update(projection_absolute_residual=absolute,
                           projection_relative_residual=relative,
                           current_coisometry_residual=current,
                           source_projector_residual=source)
                accepted[i] = bool(absolute <= projection_tolerance
                                   and current <= current_tolerance
                                   and source <= field_absolute_tolerance
                                   and self.current_residual[i] <= current_tolerance)
                if accepted[i]:
                    fields[i] = recovered
                row['status'] = 'PASS' if accepted[i] else 'FAIL: reconstructed physical mode residual'
            else:
                row['status'] = 'OPEN: packet inverse does not resolve the declared field tolerance'
            diagnostics.append(row)
        result = RecoveredInteriorModes(self.energies.copy(), fields, accepted, diagnostics)
        if require_all:
            result.require_all()
        return result


def stationary_interior_observation(energies, mass, angular, *, rtol=3e-14, atol=3e-16):
    """Build the six-by-two observation map on the owned [-1,1] collar.

    Frequencies are independent vectorized columns, not coupled modes.  The
    map uses the same radial packet normalization and PG generator as the
    archived projection owner.  Its rho=1 identity is a fundamental-solution
    convention, not a physical state or transport-time choice.
    """
    raw = np.asarray(energies)
    if raw.ndim != 1 or not len(raw) or np.iscomplexobj(raw):
        raise ValueError("nonempty real positive frequencies required")
    E = np.asarray(raw, float)
    if not np.isfinite(E).all() or np.any(E <= 0) or not np.isfinite([mass, angular]).all() or mass < 0:
        raise ValueError("finite positive frequencies and the existing massive channel required")
    if rtol <= 0 or atol <= 0:
        raise ValueError("positive ODE tolerances required")
    n = len(E)
    fundamental = np.broadcast_to(I2, (n, 2, 2)).copy()
    observation = np.zeros((n, 6, 2), complex)
    snapshots = []
    evaluations = 0
    with threadpool_limits(limits=1):
        for begin, end, active in ((1., 0., ((0, 'parent'),)),
                                   (0., -1., ((2, 'child'), (4, 'child_bulk')))):
            initial = np.r_[fundamental.ravel(), np.zeros(n*4*len(active), complex)]

            def rhs(rho, y):
                field = y[:n*4].reshape(n, 2, 2)
                base, frequency = _parts(rho, 0., mass, angular)
                derivative = np.einsum('ij,njk->nik', base, field)
                derivative += E[:, None, None]*np.einsum('ij,njk->nik', frequency, field)
                integral = np.stack([-float(radial_packet(np.array(rho), kind)[0])*field
                                     for _, kind in active], axis=1)
                return np.r_[derivative.ravel(), integral.ravel()]

            run = solve_ivp(rhs, (begin, end), initial, method='DOP853',
                            rtol=rtol, atol=atol, t_eval=[end])
            if not run.success:
                raise ArithmeticError(run.message)
            fundamental = run.y[:n*4, -1].reshape(n, 2, 2)
            integrals = run.y[n*4:, -1].reshape(n, len(active), 2, 2)
            for k, (index, _) in enumerate(active):
                observation[:, index:index+2] = integrals[:, k]
            snapshots.append(fundamental.copy())
            evaluations += run.nfev
    initial_current = S2-float(geometry(1.)[0])*I2
    residual = np.zeros(n)
    for rho, field in zip((0., -1.), snapshots):
        velocity = S2-float(geometry(rho)[0])*I2
        difference = field.swapaxes(-1, -2).conj()@velocity@field-initial_current
        residual = np.maximum(residual, np.linalg.norm(difference, axis=(1, 2)))
    return StationaryInteriorObservation(E, float(mass), float(angular), observation,
                                         snapshots[0], snapshots[1], residual, evaluations)
