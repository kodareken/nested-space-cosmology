"""KS spacetime variations of the owned transmitting PG Gaussian action.

The reference chart and PG Cauchy representation stay fixed.  KS time is
T(rho), not PG tau.  Source modes and their covariance retain the full-field
interpretation; a sampled matrix is only a quadrature of its kernel.
"""
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.integrate import quad

from .nsc_common_time_bulk_split import CommonTimeBulkSplit
from .nsc_pg_ctp_mode_jets import _columns
from .nsc_pg_ks_metric_pullback import reference_chart, pg_log_jacobian, ks_to_pg
from .nsc_transmitting_dirac_domain import I2, S1, S2, S3
from .nsc_transmitting_resolvent import profile, static_metric_kernel_from_fields


@lru_cache(maxsize=4096)
def chart_coordinates(rho):
    """T(0)=S(0)=0 fixes coordinate origins, not a transport duration."""
    rho = float(rho)
    reference_chart(rho)
    T = quad(lambda x: -1/reference_chart(x)[2], 0., rho,
             epsabs=2e-13, epsrel=2e-13)[0]
    S = quad(lambda x: reference_chart(x)[0]/reference_chart(x)[2]**2,
             0., rho, epsabs=2e-13, epsrel=2e-13)[0]
    return T, S


def reference_metric_jacobian_jet(rho):
    """The derivative of the existing reference-field Jacobian along rho."""
    b, bp, a, ap, _ = reference_chart(rho)
    r = np.sqrt(1+rho*rho)
    J = pg_log_jacobian(rho, 1., 0., a, r)
    Jp = np.array([
        [2*b*bp/a**2-2*b*b*ap/a**3, bp/a-b*ap/a**2, 3*ap/a**4, 0.],
        [2*bp/a**2-4*b*ap/a**3,
         2*b*bp/a-(b*b+1)*ap/a**2, -2*bp/a**3+6*b*ap/a**4, 0.],
        [2*ap/a**3, -bp/a+b*ap/a**2,
         2*b*bp/a**3-3*b*b*ap/a**4, 0.],
        [0., 0., 0., -rho/r**3],
    ])
    return J, Jp


def pullback_variation_jet(rho, xi, xi_T, xi_z):
    """Pointwise KS direction -> PG log-metric direction and rho derivative.

    xi is evaluated at (T(rho), tau+S(rho)). The chart is fixed during the
    variation. Complex directions denote Fourier components; a physical
    real variation includes its conjugate. No endpoint interpolation is set.
    """
    values = tuple(np.asarray(x, complex) for x in (xi, xi_T, xi_z))
    if any(x.shape != (4,) or not np.isfinite(x).all() for x in values):
        raise ValueError('four finite raw KS fields and their T,z derivatives required')
    xi, xi_T, xi_z = values
    b, _, a, _, _ = reference_chart(rho)
    J, Jp = reference_metric_jacobian_jet(rho)
    return J@xi, Jp@xi+J@(-xi_T/a+b*xi_z/a**2)


def apply_canonical_variation(rho, mass, angular, u, up, direction, radial_derivative):
    """Vary the already canonical H=-i(v d_rho+v'/2)+M.

    u=r*sqrt(q_PG)*psi is already the half-density variable. Its inherited
    Hamiltonian includes that transformation; adding a separate measure
    force here would count it twice. See nsc-adm-source-constraints.md.
    """
    d, dp = np.asarray(direction), np.asarray(radial_derivative)
    radius = np.sqrt(1+rho*rho)
    dv = (d[0]-d[2])*S2-d[1]*I2
    dvp = (dp[0]-dp[2])*S2-dp[1]*I2
    dm = -mass*d[0]*S1+angular/radius*(d[3]-d[0])*S3
    return -1j*(dv@up+.5*dvp@u)+dm@u


def transformed_pg_metric_jet(rho, delta, delta_prime):
    """Existing nonlinear KS->PG map with an arbitrary smooth direction.

    Used only for the independent Hamiltonian finite-difference check.
    Metric values are real; complex harmonic directions are checked through
    their real and imaginary parts separately.
    """
    b, bp, a, ap, _ = reference_chart(rho)
    r0 = np.sqrt(1+rho*rho)
    NK, bK, aK, rK = np.array([1., 0., a, r0])+np.asarray(delta, float)
    NKp, bKp, aKp, _ = np.array([0., 0., ap, rho/r0])+np.asarray(delta_prime, float)
    N, beta, q, r = ks_to_pg(rho, NK, bK, aK, rK)
    s = b/a**2-bK/a
    sp = bp/a**2-2*b*ap/a**3-bKp/a+bK*ap/a**2
    Fp = 2*aK*aKp*s*s+2*aK*aK*s*sp-2*NK*NKp/a**2+2*NK*NK*ap/a**3
    qp = Fp/(2*q)
    betap = (2*aK*aKp*s+aK*aK*sp)/q**2-beta*Fp/q**2
    Np = N*(NKp/NK+aKp/aK-ap/a-qp/q)
    return N, beta, q, r, Np, betap, qp


@dataclass(frozen=True)
class CompactKSHarmonic:
    """Existing compact collar direction, Fourier transformed in KS z.

    f(T(rho))=s(rho), xi=f(T) exp(-i omega z). omega is a response label,
    not an added physical oscillation, metric solution or evolution time.
    """
    omega: float

    def __post_init__(self):
        if not np.isfinite(self.omega):
            raise ValueError('finite Fourier response label required')

    def data(self, rho):
        b, _, a, _, _ = reference_chart(rho)
        s, sp = profile(rho)
        _, S = chart_coordinates(float(rho))
        phase = np.exp(-1j*self.omega*S)
        value = s*phase
        # f_T=-a*s' because T'=-1/a; z'=b/a^2 on a PG slice.
        dT, dz = -a*sp*phase, -1j*self.omega*value
        radial = -dT/a+b*dz/a**2
        return value, dT, dz, radial


def harmonic_mode_vertices(fields, mass, angular, omega, *, difference_step=2e-5):
    """Physical Fourier-direction kernel using the archived whole-field modes.

    Returned M_omega obeys M_(-omega)=M_omega^dagger. It is not itself a
    Hermitian instantaneous link or a completed history endpoint derivative.
    """
    n = fields['fields'].shape[1]*3
    weak = np.zeros((4, n, n), complex)
    strong = weak.copy()
    finite = weak.copy()
    direction = CompactKSHarmonic(omega)
    owner = CommonTimeBulkSplit()
    for rho, weight, u, up in zip(fields['rho'], fields['weights'], fields['fields'], fields['derivatives'], strict=True):
        U, Up = _columns(u), _columns(up)
        value, dT, dz, radial = direction.data(rho)
        J, _ = reference_metric_jacobian_jet(rho)
        kernel = static_metric_kernel_from_fields(rho, mass, angular, U, Up, U, Up)
        weak -= weight*value*np.einsum('AB,Aij->Bij', J, kernel)
        for B in range(4):
            e = np.eye(4)[B]
            d, dp = pullback_variation_jet(rho, e*value, e*dT, e*dz)
            strong[B] += weight*U.conj().T@apply_canonical_variation(rho, mass, angular, U, Up, d, dp)
            if difference_step is not None:
                for part, factor in ((np.real, 1.), (np.imag, 1j)):
                    columns = []
                    for sign in (1., -1.):
                        metric = transformed_pg_metric_jet(rho, sign*difference_step*e*part(value),
                                                           sign*difference_step*e*part(radial))
                        N, beta, q, r, Np, bp, qp = metric
                        columns.append(np.column_stack([
                            owner.apply_local_expression(v, vp, N=N, q_PG=q, beta=beta, radius=r,
                                N_prime=Np, q_prime=qp, beta_prime=bp, compact_mass=mass, angular_eigenvalue=angular)
                            for v, vp in zip(U.T, Up.T)]))
                    finite[B] += weight*factor*U.conj().T@(columns[0]-columns[1])/(2*difference_step)
    return weak, strong, finite if difference_step is not None else None


def harmonic_dyson_derivative(energies, vertex, omega, initial_time, final_time):
    """Derivative of U for the REAL cosine direction on caller PG times.

    Exact time integration of the reference Duhamel kernel. The source-space
    energy measure must already be applied to a finite quadrature vertex.
    No metric is stepped and no physical history/duration is selected.
    """
    E = np.asarray(energies, float)
    M = np.asarray(vertex, complex)
    if E.ndim != 1 or M.shape[-2:] != (len(E), len(E)) or not np.isfinite(E).all() or not np.isfinite(M).all():
        raise ValueError('matching finite energy labels and kernel required')
    if not np.isfinite([omega, initial_time, final_time]).all() or final_time <= initial_time:
        raise ValueError('ordered caller time coordinates required; this function selects none')
    dt = final_time-initial_time
    difference = E[:, None]-E[None, :]
    def integral(sign):
        s = difference+sign*omega
        return dt*np.exp(1j*s*dt/2+1j*sign*omega*initial_time)*np.sinc(s*dt/(2*np.pi))
    generator = .5*(M*integral(-1)+M.swapaxes(-1, -2).conj()*integral(1))
    phases = np.exp(-1j*E*dt)
    return np.diag(phases), -1j*phases[..., :, None]*generator, generator
