"""Metric tangents of the same nonlinear spatial Dirac history.

These are whole-field mode jets for supplied metric directions. Sampled
source-energy kernels remain kernel samples, not a closed branch unitary.
The existing finite Gaussian differential is used only for its relative-
branch first variation after the spatial field contraction has been formed.
"""
import numpy as np
from scipy.linalg import block_diag as dense_block_diag
from scipy.sparse import bmat, csr_matrix, diags, block_diag
from scipy.sparse.linalg import expm_multiply

from .nsc_ks_spacetime_variation import CompactKSHarmonic
from .nsc_pg_ks_metric_pullback import pg_to_ks, pg_log_jacobian
from .nsc_transmitting_ctp_variation import ctp_first_variation
from .nsc_transmitting_history_modes import characteristic_columns, RelativeModePropagator
from .nsc_transmitting_resolvent import profile

SBP4_COEFFICIENT_SOURCE = (
    'https://github.com/ranocha/SummationByPartsOperators.jl/blob/'
    '21414744283d75d6172cd1763dd5becff5988909/src/SBP_coefficients/'
    'MattssonNordstr%C3%B6m2004.jl#L36-L74'
)


class FourthOrderModePropagator(RelativeModePropagator):
    """Imported diagonal-norm SBP(4,2), with the same NSC split/SAT owner.

    Coefficients: Mattsson--Nordstrom (2004), first derivative order four;
    see SummationByPartsOperators.jl, MattssonNordstrom2004 coefficient table.
    This changes numerical accuracy only. The continuum Dirac expression,
    characteristic frame, metric fields and transmitting domain are unchanged.
    """
    def __init__(self, x, mass, angular):
        super().__init__(x, mass, angular)
        n = len(self.x)
        D = diags([np.full(n-2, 1/12), np.full(n-1, -2/3),
                   np.full(n-1, 2/3), np.full(n-2, -1/12)],
                  [-2, -1, 1, 2], shape=(n, n)).tolil()
        boundary = (
            (-24/17, 59/34, -4/17, -3/34),
            (-1/2, 0., 1/2),
            (4/43, -59/86, 0., 59/86, -4/43),
            (3/98, 0., -59/98, 0., 32/49, -4/49),
        )
        for row, coefficients in enumerate(boundary):
            D[row, :] = 0.
            D[n-1-row, :] = 0.
            for column, value in enumerate(coefficients):
                D[row, column] = value
                D[n-1-row, n-1-column] = -value
        self.D = (D/self.h).tocsr()
        edge_weights = np.array([17/48, 59/48, 43/48, 49/48])
        self.weights = np.full(n, self.h)
        self.weights[:4] = self.h*edge_weights
        self.weights[-4:] = self.h*edge_weights[::-1]
        self.norm_weights = np.tile(self.weights, 2)
        W = diags(self.weights)
        edge = np.zeros(n); edge[0], edge[-1] = -1., 1.
        residual = W@self.D+self.D.T@W-diags(edge)
        self.sbp_residual = float(np.max(abs(residual.data))) if residual.nnz else 0.
        if self.sbp_residual > 3e-11:
            raise ArithmeticError('imported SBP coefficients do not satisfy their norm identity')


def raw_KS_amplitude_directions(provider, time, x, metric):
    """Actual-metric Jacobian for the four declared response amplitudes.

    This differentiates the supplied smooth history family. It does not
    identify its amplitude coordinates with the eight old KS end nodes.
    """
    n = len(x)
    directions = np.zeros((4, 4, n))
    window = 1.
    if provider.support_time is not None:
        lo, hi = provider.support_time
        window = profile(2*(time-lo)/(hi-lo)-1)[0]
    harmonic = CompactKSHarmonic(provider.omega)
    for j, rho in enumerate(x):
        if abs(rho) >= 1.:
            continue
        raw = pg_to_ks(rho, *metric[:, j])
        J = pg_log_jacobian(rho, *raw)
        f = harmonic.data(float(rho))[0]
        envelope = window*np.real(np.exp(-1j*provider.omega*time)*f)
        directions[:, :, j] = envelope*J.T
    return directions


def generator_direction(owner, metric, log_direction):
    """dL at the actual metric, in its existing characteristic half-density.

    The four direction components are logN,beta,logq,logr. Compact radial
    directions leave the computational inflow unchanged, hence dSAT=0.
    """
    metric = np.asarray(metric, float)
    if np.iscomplexobj(log_direction) and np.any(np.imag(log_direction) != 0):
        raise ValueError('a real physical metric direction is required')
    d = np.asarray(log_direction, float)
    if d.shape != metric.shape or d.shape != (4, len(owner.x)) or not np.isfinite(d).all():
        raise ValueError('four finite metric directions on the same spatial grid required')
    if np.max(abs(d[:, [0, -1]])) > 3e-11:
        raise ValueError('this derivative retains the reference exterior and its inflow')
    N, _, q, r = metric
    dv = np.array([N/q*(d[0]-d[2])-d[1], -N/q*(d[0]-d[2])-d[1]])
    transport = block_diag([-.5*(diags(v)@owner.D+owner.D@diags(v)) for v in dv], format='csr')
    mass = N*owner.mass*d[0]
    angular = N*owner.angular/r*(d[0]-d[3])
    potential = bmat([[None, diags(1j*mass-angular)],
                      [diags(1j*mass+angular), None]], format='csr')
    tangent = transport+potential
    W = diags(owner.norm_weights)
    residual = W@tangent+tangent.conj().T@W
    error = float(np.max(abs(residual.data))) if residual.nnz else 0.
    return tangent, error


def incoming_reference_columns(owner, phi):
    """Same global reference modes supply the two right incoming traces."""
    n = len(owner.x)
    N, beta, q, _ = owner.reference
    velocities = np.array([N/q-beta, -N/q-beta])
    forcing = np.zeros_like(phi)
    for spin in range(2):
        row = spin*n+n-1
        forcing[row] = -velocities[spin, -1]/owner.weights[-1]*phi[row]
    return forcing


def evolve_field_jets(owner, energies, reference_fields, times, provider, *, with_jets=True):
    """Differentiate the nonlinear, midpoint spatial evolution itself.

    The total field and its reference are evolved by the same SBP operator
    convention, with the same prescribed exterior mode traces. This makes
    the finite variational current identity explicit; reference-discretization
    drift is reported, not subtracted as a physical stress or hidden source.
    The tangent obeys Ydot=L[g]Y+dL[g]Phi_g, with fixed initial preparation.
    """
    if np.iscomplexobj(energies) and np.any(np.imag(energies) != 0):
        raise ValueError('real source energies required')
    E = np.asarray(energies, float)
    t = np.asarray(times, float)
    if t.ndim != 1 or len(t) < 2 or not np.isfinite(t).all() or np.any(np.diff(t) <= 0):
        raise ValueError('ordered caller time coordinates required')
    fields = np.asarray(reference_fields, complex)
    if (E.ndim != 1 or not np.isfinite(E).all()
        or fields.shape != (len(owner.x), len(E), 2, 3) or not np.isfinite(fields).all()):
        raise ValueError('authenticated reference columns on the same spatial grid required')
    phi = characteristic_columns(fields)
    rows, sources = phi.shape
    energy = np.repeat(E, 3)
    L0, baseline = owner.generator(owner.reference)
    if -1-baseline['maximum_absolute_speed']*(t[-1]-t[0])-owner.x[0] <= 0:
        raise ValueError('reference causal buffer exhausted')
    incoming = csr_matrix(incoming_reference_columns(owner, phi))
    ndir = 4 if with_jets else 0
    blocks_count = 2+ndir
    y = np.vstack((phi, phi, *(np.zeros_like(phi) for _ in range(ndir)), np.eye(sources)))
    max_flux = max_dflux = max_speed = max_boundary = 0.
    for left, right in zip(t[:-1], t[1:]):
        midpoint = (left+right)/2
        metric = provider.values(midpoint, owner.x)
        L, checks = owner.generator(metric)
        operator_blocks = [[None]*(blocks_count+1) for _ in range(blocks_count+1)]
        operator_blocks[0][0], operator_blocks[0][-1] = L0, incoming
        operator_blocks[1][1], operator_blocks[1][-1] = L, incoming
        if with_jets:
            directions = raw_KS_amplitude_directions(provider, midpoint, owner.x, metric)
            for index, direction in enumerate(directions):
                dL, error = generator_direction(owner, metric, direction)
                operator_blocks[index+2][1] = dL
                operator_blocks[index+2][index+2] = L
                max_dflux = max(max_dflux, error)
        operator_blocks[-1][-1] = diags(-1j*energy)
        augmented = bmat(operator_blocks, format='csr')
        dt = right-left
        y = expm_multiply(dt*augmented, y, traceA=dt*augmented.diagonal().sum())
        max_flux = max(max_flux, checks['norm_flux_algebra'])
        max_speed = max(max_speed, checks['maximum_absolute_speed'])
        for index in range(ndir):
            tangent = y[(index+2)*rows:(index+3)*rows]
            max_boundary = max(max_boundary, float(np.max(abs(tangent[[0, len(owner.x)-1, len(owner.x), rows-1]]))))
    reference, evolved = y[:rows], y[rows:2*rows]
    tangents = y[2*rows:blocks_count*rows].reshape(ndir, rows, sources)
    overlaps = np.array([evolved.conj().T@(owner.norm_weights[:, None]*v) for v in tangents])
    compact = np.tile(abs(owner.x) < 1., 2)
    analytic_reference = phi*np.exp(-1j*energy*(t[-1]-t[0]))[None, :]
    buffer = float(-1-max_speed*(t[-1]-t[0])-owner.x[0])
    if buffer <= 0:
        raise ValueError('sampled causal buffer exhausted')
    return {
        'evolved_field': evolved, 'reference_grid_field': reference,
        'tangent_fields': tangents, 'relative_branch_overlap_jets': overlaps,
        'norm_flux_residual': max_flux, 'metric_flux_tangent_residual': max_dflux,
        'tangent_boundary_trace': max_boundary, 'sampled_causal_buffer': buffer,
        'relative_overlap_tangent_residual': (
            float(np.max(abs(overlaps+overlaps.swapaxes(-1, -2).conj()))) if ndir else None),
        'reference_grid_drift_in_collar': float(np.max(abs(reference[compact]-analytic_reference[compact]))),
        'phase_residual': float(np.max(abs(y[blocks_count*rows:]-np.diag(np.exp(-1j*energy*(t[-1]-t[0])))))),
        'full_spectral_covariance': None, 'physical_endpoint_family': None,
        'physical_history_selected': False,
    }


def relative_branch_ctp_control(energies, weights, source_covariances, projectors, overlaps):
    """Use the existing Gaussian differential in the current-history frame.

    A_B=U[g]^dagger dU[g] is formed from spatial fields before any sampled
    energy quadrature. W_±=U[g]^dagger U_± have equal-history W_±=I and
    dW_±=±A_B/2. This is a first-variation control, not exponentiation of a
    finite sampled mode kernel or a closed replacement for the bulk.
    """
    E, w = np.asarray(energies, float), np.asarray(weights, float)
    if w.shape != E.shape or not np.isfinite(w).all() or np.any(w <= 0):
        raise ValueError('positive inherited source quadrature weights required')
    P = dense_block_diag(*projectors)
    overlaps = np.asarray(overlaps, complex)
    if overlaps.shape[-2:] != P.shape or not np.isfinite(overlaps).all():
        raise ValueError('finite source-overlap kernels in the authenticated basis required')
    if max(np.linalg.norm(P@A@P-A) for A in overlaps) > 3e-11:
        raise ValueError('closed source fibers cannot be promoted to physical columns')
    active = np.flatnonzero(np.diag(P).real > .5)
    C = dense_block_diag(*source_covariances)[np.ix_(active, active)]
    root = np.repeat(np.sqrt(w/(2*np.pi)), 3)
    I = np.eye(len(active))
    results = []
    for overlap in overlaps:
        A = (root[:, None]*overlap*root[None, :])[np.ix_(active, active)]
        row = ctp_first_variation(C, I, I, .5*A, -.5*A)
        row['trace_residual'] = float(abs(row['derivative']+1j*np.trace(C@A)))
        results.append(row)
    return results
