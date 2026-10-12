"""Analytic subgap bilinears of the nonlinear transmitting history.

Complex frequencies here label analytic columns and their conjugate-energy
duals. They are never passed to a Gaussian covariance/state constructor.
The original massive horizon preparation and matched collar are retained.
"""
from math import log, pi

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import bmat, diags
from scipy.sparse.linalg import expm_multiply
from numpy.polynomial.chebyshev import chebfit, chebval
from numpy.polynomial.legendre import leggauss

from .nsc_complex_horizon_modes import _oriented_collar_frame
from .nsc_pg_massive_modes import MassivePGModeResolution
from .nsc_common_time_bulk_split import CommonTimeBulkSplit
from .nsc_lorentzian import geometry
from .nsc_transmitting_dirac_domain import I2, S2, MODE_TO_CURRENT
from .nsc_transmitting_history_jets import (
    raw_KS_amplitude_directions, generator_direction, incoming_reference_columns,
)


def inner_horizon_basis(z, mass, angular, preparation):
    """Only the already owned analytic interior columns, in order (v,w).

    This reuses the inner leg of finite_pg_horizon_basis; no unused exterior
    packet projection or scattering calculation is performed. The original
    inner Frobenius columns are (w,v), so the order is explicitly exchanged.
    """
    z = complex(z)
    if z.real <= 0 or abs(z.imag) >= preparation.surface_gravity/2:
        raise ValueError('positive frequency in the owned Frobenius strip required')
    p = preparation
    initial = _oriented_collar_frame(z, mass, angular, p, timelike=True)
    stop = log(p.background.horizon_q-pi/2)
    run = solve_ivp(lambda y, f: (-1j*p.generator(y, z, mass, angular)@f.reshape(2, 2)).ravel(),
                    (log(p.collar_delta_q), stop), initial.ravel(), method='DOP853',
                    rtol=2e-13, atol=2e-15, max_step=.1)
    if not run.success:
        raise ArithmeticError(run.message)
    resolution = MassivePGModeResolution(p, rtol=2e-13, atol=2e-15)
    frame, _ = resolution.frame(0.)
    value = np.exp(1j*z*float(resolution.clock(np.array(0.))))*frame@run.y[:, -1].reshape(2, 2)
    return value[:, [1, 0]]


def analytic_columns_on_grid(z, mass, angular, at_zero, x):
    """Analytic homogeneous Dirac columns on the trapped causal patch."""
    z = complex(z)
    x = np.asarray(x, float)
    initial = np.asarray(at_zero, complex)
    if initial.shape != (2, 2) or not np.isfinite(initial).all() or np.any(geometry(x)[0] <= 1):
        raise ValueError('two analytic columns on the owned trapped grid required')
    def rhs(rho, data):
        beta, bp, _ = geometry(rho)
        potential = CommonTimeBulkSplit.local_potential(1., np.sqrt(1+rho*rho), mass, angular)
        generator = np.linalg.solve(S2-beta*I2, 1j*(z*I2-potential)+.5*bp*I2)
        return (generator@data.reshape(2, 2)).ravel()
    values = np.empty((len(x), 2, 2), complex)
    for end, mask in ((x[0], x <= 0), (x[-1], x > 0)):
        run = solve_ivp(rhs, (0., float(end)), initial.ravel(), method='DOP853',
                        rtol=2e-13, atol=2e-15, dense_output=True)
        if not run.success:
            raise ArithmeticError(run.message)
        values[mask] = run.sol(x[mask]).T.reshape((-1, 2, 2))
    char = np.einsum('ab,nbc->nac', MODE_TO_CURRENT.conj().T, values)
    return char.transpose(1, 0, 2).reshape(2*len(x), 2)


def paired_history_response(owner, z, upper_columns, dual_columns, times, provider):
    """B(z)=i <U Phi(bar z), deltaU Phi(z)> for the four metric directions.

    The real metric generator is unchanged. Analytic forcing amplitudes are
    numerical column labels, not a complex-energy physical state. The bra
    is always propagated from bar(z), never obtained by conjugating Phi(z).
    """
    z = complex(z)
    upper, dual = np.asarray(upper_columns, complex), np.asarray(dual_columns, complex)
    size = 2*len(owner.x)
    if upper.shape != (size, 2) or dual.shape != upper.shape or not np.isfinite([upper, dual]).all():
        raise ValueError('analytic columns and separately prepared conjugate-energy columns required')
    t = np.asarray(times, float)
    if t.ndim != 1 or len(t) < 2 or np.any(np.diff(t) <= 0) or not np.isfinite(t).all():
        raise ValueError('ordered response coordinates required')
    phi = np.hstack((upper, dual))
    incoming = incoming_reference_columns(owner, phi)
    energy = np.array([z, z, z.conjugate(), z.conjugate()])
    columns = len(energy)
    state = np.vstack((phi, *(np.zeros_like(phi) for _ in range(4)), np.eye(columns)))
    maximum_flux = maximum_tangent = maximum_speed = maximum_boundary = 0.
    for left, right in zip(t[:-1], t[1:]):
        metric = provider.values((left+right)/2, owner.x)
        L, checks = owner.generator(metric)
        blocks = [[None]*6 for _ in range(6)]
        blocks[0][0], blocks[0][-1] = L, incoming
        for B, direction in enumerate(raw_KS_amplitude_directions(provider, (left+right)/2, owner.x, metric)):
            dL, error = generator_direction(owner, metric, direction)
            blocks[B+1][0], blocks[B+1][B+1] = dL, L
            maximum_tangent = max(maximum_tangent, error)
        blocks[-1][-1] = diags(-1j*energy)
        operator = bmat(blocks, format='csr')
        dt = right-left
        state = expm_multiply(dt*operator, state, traceA=dt*operator.diagonal().sum())
        maximum_speed = max(maximum_speed, checks['maximum_absolute_speed'])
        maximum_flux = max(maximum_flux, checks['norm_flux_algebra'])
        for B in range(4):
            y = state[(B+1)*size:(B+2)*size]
            maximum_boundary = max(maximum_boundary, float(np.max(abs(y[[0, len(owner.x)-1, len(owner.x), size-1]]))))
    phi_g = state[:size]
    jets = state[size:5*size].reshape(4, size, columns)
    responses, dual_responses = [], []
    for y in jets:
        responses.append(1j*phi_g[:, 2:].conj().T@(owner.norm_weights[:, None]*y[:, :2]))
        dual_responses.append(1j*phi_g[:, :2].conj().T@(owner.norm_weights[:, None]*y[:, 2:]))
    B, Bd = np.array(responses), np.array(dual_responses)
    buffer = float(-1-maximum_speed*(t[-1]-t[0])-owner.x[0])
    if buffer <= 0:
        raise ValueError('trapped causal patch is insufficient')
    return {
        'response': B, 'dual_response': Bd,
        'analytic_adjoint_residual': float(np.max(abs(B-Bd.swapaxes(-1, -2).conj()))),
        'metric_flux_tangent': maximum_tangent, 'norm_flux': maximum_flux,
        'tangent_boundary_trace': maximum_boundary, 'sampled_causal_buffer': buffer,
        'complex_frequency_is_state': False,
    }


def source_functions(z, kappa):
    """The inherited horizon coefficients, analytic inside the same strip."""
    z = complex(z)
    if kappa <= 0 or abs(z.imag) >= kappa/2:
        raise ValueError('owned pole-free thermal strip required')
    e = np.exp(-2*np.pi*z/kappa)
    return e/(1+e)-.5, np.exp(-np.pi*z/kappa)/(1+e)


def trapped_subgap_contractions(z, kappa, response):
    """Centered positive-energy trace, before reflection and dE/(2pi).

    u=0 on this initial causal patch. Columns are (v,w), so
    Tr(K B)=d(B_vv-B_ww)+2Re[-i s R B_wv] on the real axis.
    Only the second coefficient is continued along the retarded contour.
    """
    B = np.asarray(response, complex)
    if B.shape != (4, 2, 2) or not np.isfinite(B).all():
        raise ValueError('four analytic metric-response bilinears required')
    d, s = source_functions(z, kappa)
    return d*(B[:, 0, 0]-B[:, 1, 1]), -1j*s*B[:, 1, 0]


class AnalyticResponsePanel:
    """Interpolate the smooth unsewn bilinear, never the winding reflection.

    The numerical continuation must be checked against an independently
    propagated conjugate-energy pair. It selects no physical history.
    """
    def __init__(self, energies, responses, left, right):
        E, B = np.asarray(energies, float), np.asarray(responses, complex)
        if (not 0 < left < right or E.ndim != 1 or len(E) < 5
                or B.shape != (len(E), 4, 2, 2) or not np.isfinite(B).all()
                or np.any(E < left) or np.any(E > right)):
            raise ValueError('finite subgap response samples on the declared interval required')
        self.left, self.right = float(left), float(right)
        t = (2*E-left-right)/(right-left)
        self.coefficients = chebfit(t, B.reshape(len(E), -1), len(E)-1)

    def __call__(self, z):
        t = (2*complex(z)-self.left-self.right)/(self.right-self.left)
        return chebval(t, self.coefficients).reshape(4, 2, 2)


def integrate_subgap_response(panel, kappa, reflection, *, points=48, height=None):
    """Owned positive panel plus its zero-angular signed partner.

    Result is a bare per-reduced-family action derivative, not a local
    renormalized stress. No angular degeneracy or compact factor is inserted.
    """
    left, right = panel.left, panel.right
    height = kappa/3 if height is None else float(height)
    if points < 4 or not 0 < height < kappa/2:
        raise ValueError('resolved quadrature within the owned analytic strip required')
    x, w = leggauss(points)
    smooth = np.zeros(4, complex)
    for xi, wi in zip(x, w):
        z = left+(xi+1)*(right-left)/2
        part, _ = trapped_subgap_contractions(z, kappa, panel(z))
        smooth += wi*(right-left)/2*part
    reflected = np.zeros(4, complex)
    for a, b in ((complex(left), left+1j*height),
                 (left+1j*height, right+1j*height),
                 (right+1j*height, complex(right))):
        for xi, wi in zip(x, w):
            z = a+(xi+1)*(b-a)/2
            _, part = trapped_subgap_contractions(z, kappa, panel(z))
            reflected += wi*(b-a)/2*reflection(z)*part
    signed = -(smooth.real+2*reflected.real)/np.pi
    return {'smooth_integral': smooth, 'reflected_integral': reflected,
            'signed_action_derivative': signed,
            'smooth_imaginary_residual': float(np.max(abs(smooth.imag)))}
