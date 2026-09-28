"""Nonlinear KS weak source obtained from the owned PG field evolution.

The propagator, numerical SBP/SAT operator and source preparation are reused.
The new work is the actual-metric Cauchy restriction, its time derivative,
and the weak source integrated before any finite source-frequency sampling.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import bmat, csr_matrix, diags
from scipy.sparse.linalg import expm_multiply

from .nsc_ks_spacetime_variation import CompactKSHarmonic
from .nsc_lorentzian import geometry
from .nsc_transmitting_resolvent import profile
from .nsc_transmitting_history_modes import characteristic_columns
from .nsc_transmitting_history_jets import incoming_reference_columns, generator_direction


class OwnedMetricSamples:
    """Vectorized values and clock derivative of the existing metric family."""
    def __init__(self, x, provider):
        self.x, self.provider = np.asarray(x), provider
        self.b, _, _ = geometry(x)
        self.a = np.sqrt(self.b*self.b-1)
        self.r = np.sqrt(1+self.x*self.x)
        harmonic = CompactKSHarmonic(provider.omega)
        self.f = np.array([harmonic.data(float(rho))[0] if abs(rho) < 1 else 0j for rho in x])

    def at(self, time):
        p = self.provider
        w, wt = 1., 0.
        if p.support_time is not None:
            lo, hi = p.support_time
            w, derivative = profile(2*(time-lo)/(hi-lo)-1)
            wt = 2*derivative/(hi-lo)
        phase = np.exp(-1j*p.omega*time)*self.f
        shape = w*phase.real
        shape_t = wt*phase.real+w*(-1j*p.omega*phase).real
        A = np.asarray(p.amplitudes)
        NK, bK, aK, rK = (1+A[0]*shape, A[1]*shape,
                          self.a+A[2]*shape, self.r+A[3]*shape)
        s = self.b/self.a**2-bK/self.a
        q2 = aK*aK*s*s-NK*NK/self.a**2
        q = np.sqrt(q2); N = NK*aK/(self.a*q); beta = aK*aK*s/q2
        logq = np.array([-NK/self.a**2, -aK*aK*s/self.a, aK*s*s, np.zeros_like(s)])/q2
        logN = np.array([1/NK, np.zeros_like(s), 1/aK, np.zeros_like(s)])-logq
        dbeta = np.array([np.zeros_like(s), -aK*aK/self.a, 2*aK*s, np.zeros_like(s)])/q2-2*beta*logq
        logr = np.array([np.zeros_like(s)]*3+[1/rK])
        J = np.array([logN, dbeta, logq, logr])  # PG field, raw KS direction, x
        directions = shape[None, None, :]*J.transpose(1, 0, 2)
        rates = np.einsum('abx,b,x->ax', J, A, shape_t)
        speeds = np.array([N/q-beta, -N/q-beta])
        speed_t = np.array([N/q*(rates[0]-rates[2])-rates[1],
                           -N/q*(rates[0]-rates[2])-rates[1]])
        if np.any(speeds >= 0): raise ValueError('the owned KS Cauchy slices require the trapped domain')
        # In the PG characteristic frame the KS normal/half-density map is
        # exactly sqrt(-v_±), not a unitary pointwise spin boost.
        frame = np.sqrt(-speeds)
        frame_t = -speed_t/(2*frame)
        return {'metric': np.array([N, beta, q, rK]), 'KS': np.array([NK, bK, aK, rK]),
                'directions': directions, 'shape': shape, 'frame': frame, 'frame_t': frame_t}


def weak_KS_source(owner, phi, phi_t, sample):
    n = len(owner.x)
    u, ut = phi.reshape(2, n, -1).transpose(1, 0, 2), phi_t.reshape(2, n, -1).transpose(1, 0, 2)
    M, Mt = sample['frame'].T, sample['frame_t'].T
    psi = M[:, :, None]*u
    psi_z = Mt[:, :, None]*u+M[:, :, None]*ut
    def outer(a, b): return np.einsum('xi,xj->xij', a.conj(), b)
    kinetic = []
    for s in range(2):
        kinetic.append(-.5j*(outer(psi[:, s], psi_z[:, s])-outer(psi_z[:, s], psi[:, s])))
    K0, K3 = kinetic[0]+kinetic[1], kinetic[0]-kinetic[1]
    cross = outer(psi[:, 0], psi[:, 1])
    Q1 = cross+cross.swapaxes(-1, -2).conj()
    Q2 = -1j*cross+1j*cross.swapaxes(-1, -2).conj()
    N, beta, a, r = sample['KS']
    matrix = np.array([-owner.mass*Q1+owner.angular/r[:, None, None]*Q2+K3/a[:, None, None],
                       -K0, -N[:, None, None]/a[:, None, None]**2*K3,
                       -(N*owner.angular/r**2)[:, None, None]*Q2])
    beta0 = geometry(owner.x)[0]; a0 = np.sqrt(beta0*beta0-1)
    weights = owner.weights*sample['shape']/a0
    integral = np.einsum('x,bxij->bij', weights, matrix)
    return integral, psi, psi_z


def record_common_history(owner, energies, reference_fields, times, provider, *, quadrature=2):
    """Record full-field weak insertions, keeping only one KS slice as payload.

    The spatial bulk is evolved before any observation is retained. Recording
    rho=0 samples does not turn it into a few-mode closed Hamiltonian.
    Within each old midpoint step the same frozen exponential supplies the
    quadrature samples. Continuous source coefficients use the actual metric
    at each sample; its difference from the discrete endpoint jet is measured.
    """
    E, t = np.asarray(energies), np.asarray(times, float)
    if np.iscomplexobj(E) and np.any(E.imag != 0): raise ValueError('real source frequencies required')
    E = E.real.astype(float)
    if np.any(np.diff(t) <= 0) or quadrature < 2: raise ValueError('ordered history and resolved insertion quadrature required')
    phi = characteristic_columns(reference_fields)
    rows, cols = phi.shape
    incoming = csr_matrix(incoming_reference_columns(owner, phi))
    frequency = np.repeat(E, 3)
    state = np.vstack((phi, np.eye(cols)))
    metric = OwnedMetricSamples(owner.x, provider)
    nodes, weights = leggauss(quadrature)
    weak_ks = np.zeros((4, cols, cols), complex); weak_pg = weak_ks.copy()
    trace_samples, psi_samples, derivative_samples, sample_times = [], [], [], []
    zero = int(np.argmin(abs(owner.x)))
    if abs(owner.x[zero]) > 3e-12: raise ValueError('recorded KS slice must be an actual spatial node')
    maximum_flux = maximum_speed = 0.
    for left, right in zip(t[:-1], t[1:]):
        midpoint = (left+right)/2; dt = right-left
        smid = metric.at(midpoint)
        Lmid, check = owner.generator(smid['metric'])
        operator = bmat([[Lmid, incoming], [None, diags(-1j*frequency)]], format='csr')
        for xi, wi in zip(nodes, weights):
            theta = (xi+1)/2; tau = left+theta*dt
            at = expm_multiply(theta*dt*operator, state, traceA=theta*dt*operator.diagonal().sum())
            current, phases = at[:rows], at[rows:]
            s = metric.at(tau); L, local = owner.generator(s['metric'])
            derivative = L@current+incoming@phases
            ks, psi, psi_z = weak_KS_source(owner, current, derivative, s)
            pg = np.array([1j*current.conj().T@(owner.norm_weights[:, None]*(generator_direction(owner, s['metric'], direction)[0]@current))
                           for direction in s['directions']])
            weight = wi*dt/2
            weak_ks += weight*ks; weak_pg += weight*pg
            sample_times.append(tau); trace_samples.append(ks)
            psi_samples.append(psi[zero]); derivative_samples.append(psi_z[zero])
            maximum_flux = max(maximum_flux, local['norm_flux_algebra'])
            maximum_speed = max(maximum_speed, local['maximum_absolute_speed'])
        state = expm_multiply(dt*operator, state, traceA=dt*operator.diagonal().sum())
    buffer = -1-maximum_speed*(t[-1]-t[0])-owner.x[0]
    if buffer <= 0: raise ValueError('inherited causal buffer exhausted')
    return {'weak_KS': weak_ks, 'weak_PG': weak_pg, 'evolved_field': state[:rows],
            'record_times': np.array(sample_times), 'KS_slice_fields': np.array(psi_samples),
            'KS_slice_z_derivatives': np.array(derivative_samples), 'time_weak_KS': np.array(trace_samples),
            'flux_residual': maximum_flux, 'causal_buffer': float(buffer),
            'phase_residual': float(np.max(abs(state[rows:]-np.diag(np.exp(-1j*frequency*(t[-1]-t[0]))))))}
