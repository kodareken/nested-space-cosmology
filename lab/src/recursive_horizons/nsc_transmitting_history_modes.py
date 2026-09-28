"""Conditional transmitting mode evolution on supplied PG metric histories.

The unknown propagated object is a spatial field difference, not a closed
matrix of the sampled energies. For each reference input mode Phi_E,
w_E=(U[g]-U_ref)Phi_E has zero initial data and a compact source. The supplied
geometry is a propagator argument, never a metric solution selected here.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import bmat, csr_matrix, diags, block_diag
from scipy.sparse.linalg import expm_multiply

from .nsc_common_time_bulk_split import CommonTimeBulkSplit
from .nsc_ks_spacetime_variation import CompactKSHarmonic, transformed_pg_metric_jet
from .nsc_lorentzian import geometry
from .nsc_transmitting_dirac_domain import I2, S2, MODE_TO_CURRENT
from .nsc_transmitting_resolvent import profile


@dataclass(frozen=True)
class SuppliedKSHarmonicMetric:
    """Nonlinear metric for an explicit response control; no physical selector.

    Four real amplitudes perturb the raw KS fields through the owned frozen
    chart. The example is xi_B=s(rho) cos[omega(tau+S(rho))]. Neither these
    amplitudes nor the caller's time interval are model-parameter fits.
    """
    amplitudes: tuple[float, float, float, float]
    omega: float = .4
    support_time: tuple[float, float] | None = None

    def __post_init__(self):
        a = np.asarray(self.amplitudes, float)
        if a.shape != (4,) or not np.isfinite(a).all() or not np.isfinite(self.omega):
            raise ValueError('four finite real raw KS amplitudes required')
        if self.support_time is not None:
            t = np.asarray(self.support_time, float)
            if t.shape != (2,) or not np.isfinite(t).all() or t[1] <= t[0]:
                raise ValueError('ordered response-window coordinates required')

    def values(self, time, x):
        direction = CompactKSHarmonic(self.omega)
        window = 1.
        if self.support_time is not None:
            lo, hi = self.support_time
            window = profile(2*(time-lo)/(hi-lo)-1)[0]
        coefficients = []
        for rho in x:
            if abs(rho) >= 1.:
                beta = float(geometry(rho)[0])
                coefficients.append((1., beta, 1., np.sqrt(1+rho*rho)))
            else:
                f, _, _, fp = direction.data(float(rho))
                phase = np.exp(-1j*self.omega*time)
                amplitude = window*np.asarray(self.amplitudes)
                N, beta, q, r, *_ = transformed_pg_metric_jet(
                    rho, amplitude*np.real(phase*f), amplitude*np.real(phase*fp))
                coefficients.append((N, beta, q, r))
        return np.asarray(coefficients).T


def reference_fields_on_grid(energies, mass, angular, at_zero, x):
    """Evaluate the owned radial mode equation on the new propagation grid.

    Initial columns are the authenticated horizon/infinity mode maps, not
    seed covariance. This changes spatial sampling only; it reuses the
    preparation and signed frequencies of nsc_pg_ctp_mode_jets.
    """
    if np.iscomplexobj(energies) and np.any(np.imag(energies) != 0):
        raise ValueError('real source energies required; contour data need their own analytic covariance representation')
    E = np.asarray(energies, float)
    initial = np.asarray(at_zero, complex)
    x = np.asarray(x, float)
    if (E.ndim != 1 or not np.isfinite(E).all() or not np.isfinite(initial).all()
        or not np.isfinite(x).all() or initial.shape != (len(E), 2, 3)
        or np.any(np.diff(x) <= 0) or not x[0] < 0 < x[-1]):
        raise ValueError('ordered two-sided grid and authenticated source columns required')
    if np.any(geometry(x)[0] <= 1):
        raise ValueError('this relative-field control is wholly inside the owned trapped chart')
    def rhs(rho, data):
        beta, bp, _ = geometry(rho)
        potential = CommonTimeBulkSplit.local_potential(1., np.sqrt(1+rho*rho), mass, angular)
        iv = np.linalg.inv(S2-beta*I2)
        generator = np.einsum('ij,njk->nik', iv, 1j*(E[:, None, None]*I2-potential)+.5*bp*I2)
        return np.einsum('nij,njk->nik', generator, data.reshape(initial.shape)).ravel()
    fields = np.empty((len(x), len(E), 2, 3), complex)
    for endpoint, mask in ((x[0], x <= 0), (x[-1], x > 0)):
        run = solve_ivp(rhs, (0., float(endpoint)), initial.ravel(), method='DOP853',
                        rtol=2e-13, atol=2e-15, dense_output=True)
        if not run.success:
            raise ArithmeticError(run.message)
        fields[mask] = run.sol(x[mask]).T.reshape((-1, *initial.shape))
    return fields


def characteristic_columns(fields):
    """Current spin frame -> characteristic frame, preserving the field."""
    values = np.einsum('ab,nebc->neac', MODE_TO_CURRENT.conj().T, fields)
    return values.transpose(2, 0, 1, 3).reshape(2*len(fields), -1)


class RelativeModePropagator:
    """Spatial SBP discretization of the owned canonical Dirac equation.

    Reuses the split/SAT convention in nsc_lorentzian. Right inflow of the
    DIFFERENCE field is zero because the metric perturbation is compact and
    the initial difference is zero. Left outflow is accounted for, not made
    reflecting. The matrix has spatial degrees of freedom; source-frequency
    columns are forcing labels and do not replace the unobserved bulk.
    """
    def __init__(self, x, mass, angular):
        self.x = np.asarray(x, float)
        if self.x.ndim != 1 or len(x) < 9 or not np.isfinite(self.x).all():
            raise ValueError('finite spatial grid required')
        spacing = np.diff(self.x)
        if np.min(spacing) <= 0 or np.max(abs(spacing-spacing[0])) > 1e-12:
            raise ValueError('uniform increasing SBP grid required')
        if not x[0] < -1 < 1 < x[-1]:
            raise ValueError('a propagation buffer outside the existing compact source is required')
        self.h = float(spacing[0])
        self.mass, self.angular = float(mass), float(angular)
        n = len(x)
        D = diags([-np.ones(n-1), np.ones(n-1)], [-1, 1], shape=(n, n)).tolil()/(2*self.h)
        D[0, 0], D[0, 1] = -1/self.h, 1/self.h
        D[-1, -2], D[-1, -1] = -1/self.h, 1/self.h
        self.D = D.tocsr()
        self.weights = np.full(n, self.h)
        self.weights[[0, -1]] *= .5
        self.norm_weights = np.tile(self.weights, 2)
        self.reference = np.array([np.ones(n), geometry(self.x)[0], np.ones(n), np.sqrt(1+self.x*self.x)])

    def generator(self, metric):
        values = np.asarray(metric, float)
        n = len(self.x)
        if values.shape != (4, n) or not np.isfinite(values).all():
            raise ValueError('N,beta,q,r on the same spatial grid required')
        N, beta, q, r = values
        if min(N.min(), q.min(), r.min()) <= 0 or np.any(q*beta <= N):
            raise ValueError('positive metric in the connected trapped chart required')
        if np.max(abs(values[:, [0, -1]]-self.reference[:, [0, -1]])) > 3e-11:
            raise ValueError('the supplied compact history must retain the reference exterior')
        speeds = np.array([N/q-beta, -N/q-beta])
        blocks = []
        boundary = np.zeros(2*n)
        for c, speed in enumerate(speeds):
            V = diags(speed)
            transport = -.5*(V@self.D+self.D@V)
            sat = np.zeros(n)
            for index, normal in ((0, -1), (n-1, 1)):
                nv = normal*speed[index]
                if nv < 0:
                    sat[index] = nv/self.weights[index]
                boundary[c*n+index] = -abs(nv)
            blocks.append(transport+diags(sat))
        potential = bmat([[None, diags(1j*N*self.mass-N*self.angular/r)],
                          [diags(1j*N*self.mass+N*self.angular/r), None]], format='csr')
        L = block_diag(blocks, format='csr')+potential
        W = diags(self.norm_weights)
        residual = W@L+L.conj().T@W-diags(boundary)
        error = float(np.max(abs(residual.data))) if residual.nnz else 0.
        return L, {'norm_flux_algebra': error,
                   'maximum_characteristic_speed': float(speeds.max()),
                   'maximum_absolute_speed': float(abs(speeds).max())}

    def evolve(self, energies, reference_fields, times, metric_provider):
        """Compute (U[g]-U_ref)Phi_E for supplied geometry, not metric evolution.

        A sparse variation-of-constants augmentation integrates the known
        exp(-iEt) forcing. Its auxiliary amplitudes are solver variables,
        never particles, reservoirs, or additional terms in the action.
        Midpoint coefficient sampling is the declared temporal approximation.
        """
        if np.iscomplexobj(energies) and np.any(np.imag(energies) != 0):
            raise ValueError('real source energies required; contour data are not canonical input frequencies')
        E = np.asarray(energies, float)
        time = np.asarray(times, float)
        if time.ndim != 1 or len(time) < 2 or not np.isfinite(time).all() or np.any(np.diff(time) <= 0):
            raise ValueError('ordered caller time coordinates required; no duration is selected')
        fields = np.asarray(reference_fields, complex)
        if (E.ndim != 1 or not np.isfinite(E).all() or not np.isfinite(fields).all()
            or fields.shape != (len(self.x), len(E), 2, 3)):
            raise ValueError('reference source-mode fields must use the same spatial grid')
        phi = characteristic_columns(fields)
        energy = np.repeat(E, 3)
        sources = len(energy)
        count = len(phi)
        L0, reference_checks = self.generator(self.reference)
        if -1-reference_checks['maximum_absolute_speed']*(time[-1]-time[0])-self.x[0] <= 0:
            raise ValueError('enlarge the computational field domain: reference causal buffer is exhausted')
        y = np.vstack((np.zeros_like(phi), np.eye(sources, dtype=complex)))
        maximum_flux = reference_checks['norm_flux_algebra']
        maximum_speed = reference_checks['maximum_absolute_speed']
        maximum_right_speed = reference_checks['maximum_characteristic_speed']
        left_trace = right_trace = 0.
        for left, right in zip(time[:-1], time[1:]):
            L, checks = self.generator(metric_provider.values((left+right)/2, self.x))
            forcing = (L-L0)@phi
            # This block matrix implements the inhomogeneous linear PDE.
            augmented = bmat([[L, csr_matrix(forcing)],
                              [csr_matrix((sources, count)), diags(-1j*energy)]], format='csr')
            dt = right-left
            y = expm_multiply(dt*augmented, y, traceA=dt*augmented.diagonal().sum())
            left_trace = max(left_trace, float(np.max(abs(y[[0, len(self.x)], :]))))
            right_trace = max(right_trace, float(np.max(abs(y[[len(self.x)-1, count-1], :]))))
            maximum_flux = max(maximum_flux, checks['norm_flux_algebra'])
            maximum_speed = max(maximum_speed, checks['maximum_absolute_speed'])
            maximum_right_speed = max(maximum_right_speed, checks['maximum_characteristic_speed'])
        difference = y[:count]
        norm = np.einsum('i,ij,ij->j', self.norm_weights, difference.conj(), difference).real
        mode_projection = phi.conj().T@(self.norm_weights[:, None]*difference)
        buffer = float(-1-maximum_speed*(time[-1]-time[0])-self.x[0])
        if buffer <= 0:
            raise ValueError('enlarge the computational field domain: sampled causal buffer is exhausted')
        return {
            'difference_characteristic_field': difference,
            'reference_characteristic_field': phi,
            'mode_projection': mode_projection,
            'difference_L2_norms': np.sqrt(np.maximum(norm, 0.)),
            'norm_flux_algebra': maximum_flux,
            'forcing_phase_residual': float(np.max(abs(y[count:]-np.diag(np.exp(-1j*energy*(time[-1]-time[0])))))),
            'left_difference_trace': left_trace,
            'right_difference_trace': right_trace,
            'sampled_speed_buffer': buffer,
            'maximum_right_characteristic_speed': maximum_right_speed,
            'initial_metric_match': float(np.max(abs(metric_provider.values(time[0], self.x)-self.reference))),
            'initial_time_jet_match': (
                'all orders: inherited C-infinity bump vanishes at initial endpoint'
                if isinstance(metric_provider, SuppliedKSHarmonicMetric)
                and metric_provider.support_time is not None
                and time[0] <= metric_provider.support_time[0]
                else 'not certified; propagator alone does not specify a renormalizable preparation'),
            'physical_history_selected': False,
            'absolute_stress': None,
            'full_covariance_integral': None,
            'complete_EndpointBranchJets': None,
        }
