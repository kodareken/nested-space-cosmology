"""Existing LLL conformal allocation in the canonical KS history frame.

This is T_covariant - T_spatial-normal-ordered for the already retained
massless magnetic sector. It replaces its stationary geometric allocation;
it does not add a state, a field, an adjustable coefficient or Gamma_rest.
The local correction alone is reference/frame dependent. Only its sum with
the canonical source is the covariant LLL stress.
"""
from dataclasses import dataclass
import numpy as np

from .nsc_lorentzian import geometry
from .nsc_spherical_local_history import spherical_invariants


def raw_KS_geometry(history, amplitudes):
    """Pull existing exact PG metric jets back to their raw KS coordinates.

    T=T(rho), z=tau+S(rho); d_T=(beta0/a0)d_tau-a0*d_rho,
    d_z=d_tau. The chart stays fixed under the physical metric variation.
    Curvature comes from the existing spherical owner, not a new GR owner.
    """
    g, radius, jets = history.fields_and_jets(amplitudes)
    dg, ddg, _, _ = jets
    a0, b0 = history.axial, history.beta
    bp = geometry(history.grid.radius_coordinate)[1]
    ap = b0*bp/a0
    K = np.zeros((len(a0), 2, 2)); Kp = np.zeros_like(K)
    K[:, 0, 0], K[:, 0, 1], K[:, 1, 0] = b0/a0, 1., -a0
    Kp[:, 0, 0], Kp[:, 1, 0] = bp/a0-b0*ap/a0**2, -ap

    def contract(left, value, right):
        return np.einsum('xia,txij,xjb->txab', left, value, right, optimize=True)

    G = contract(K, g, K)
    Gt = contract(K, dg[..., 0, :, :], K)
    Gx = (contract(Kp, g, K)+contract(K, dg[..., 1, :, :], K)
          +contract(K, g, Kp))
    Gtt = contract(K, ddg[..., 0, 0, :, :], K)
    Gtx = (contract(Kp, dg[..., 0, :, :], K)
           +contract(K, ddg[..., 0, 1, :, :], K)
           +contract(K, dg[..., 0, :, :], Kp))
    a = np.sqrt(-G[..., 1, 1])
    beta = -G[..., 0, 1]/a**2
    N = np.sqrt(G[..., 0, 0]+a*a*beta*beta)
    az, arho = -Gt[..., 1, 1]/(2*a), -Gx[..., 1, 1]/(2*a)
    azz = -Gtt[..., 1, 1]/(2*a)-az*az/a
    azrho = -Gtx[..., 1, 1]/(2*a)-az*arho/a
    betaz = -Gt[..., 0, 1]/a**2-2*beta*az/a
    betazz = (-Gtt[..., 0, 1]/a**2+4*Gt[..., 0, 1]*az/a**3
              +2*G[..., 0, 1]*azz/a**3-6*G[..., 0, 1]*az*az/a**4)
    Nz = (Gt[..., 0, 0]+2*a*az*beta*beta+2*a*a*beta*betaz)/(2*N)
    aT = (b0/a0)*az-a0*arho
    aTz = (b0/a0)*azz-a0*azrho
    H = (aT-beta*az-a*betaz)/(N*a)
    Hz = (aTz-beta*azz-2*betaz*az-a*betazz)/(N*a)-H*(Nz/N+az/a)
    Lz, Lzz = az/a, azz/a-(az/a)**2
    curvature = spherical_invariants(history.grid, g, radius, jets)
    return {'N': N, 'beta': beta, 'a': a, 'r': radius, 'Nz': Nz,
            'Lz': Lz, 'Lzz': Lzz, 'H': H, 'Hz': Hz, 'R2': curvature['R2'],
            'aT': aT, 'az': az, 'betaz': betaz,
            'measure_KS_in_PG': 1/a0}


@dataclass(frozen=True)
class LLLGeometricHistoryAllocation:
    """One-counted c=|q_mag| correction in the physical radial metric g2.

    S_delta=c/(24pi) integral[-Na H^2+2 Nz Lz/a-N Lz^2/a] dT dz.
    Compact, boundary-identical variations are supplied by the existing
    OwnedCompactKSHistory. This API supplies an allocation, never full stress.
    """
    magnetic_flux: int

    def __post_init__(self):
        if isinstance(self.magnetic_flux, bool) or not np.isfinite(self.magnetic_flux) or self.magnetic_flux == 0 or int(self.magnetic_flux) != self.magnetic_flux:
            raise ValueError('nonzero integer magnetic flux from the existing ledger required')

    def evaluate(self, history, amplitudes):
        q = raw_KS_geometry(history, amplitudes)
        c = abs(self.magnetic_flux)/(24*np.pi)
        N, a, r = (q[k] for k in ('N', 'a', 'r'))
        rho = c*((2*q['Lzz']-q['Lz']**2)/a**2-q['H']**2)
        flux = 2*c*q['Hz']/a
        pressure = rho-c*q['R2']
        zero = np.zeros_like(rho)
        stress = np.stack((rho, flux, pressure, zero), axis=-1)
        density = c*(-N*a*q['H']**2+2*q['Nz']*q['Lz']/a-N*q['Lz']**2/a)
        # Euler derivatives of this same action, with signature +---.
        gradient_density = np.stack((-a*rho, -a*a*flux, N*pressure, zero), axis=-1)
        gradient = history.grid.integral(gradient_density*history.envelope[..., None]
                                        *q['measure_KS_in_PG'][None, :, None])
        action = history.grid.integral(density*q['measure_KS_in_PG'])
        return {'action': action, 'action_gradient_from_stress': gradient,
                'action_force_from_stress': -gradient,
                'stress_2D_allocation': stress,
                'stress_4D_allocation': stress/(4*np.pi*r[..., None]**2),
                'gradient_density': gradient_density, 'geometry': q}

    def action_gradient(self, history, amplitudes, step=1e-28):
        """Independent complex-step variation of the action, not stress."""
        if not np.isfinite(step) or step <= 0:
            raise ValueError('positive finite complex step required')
        base = np.asarray(amplitudes, float)
        if base.shape != (4,) or not np.isfinite(base).all():
            raise ValueError('four finite raw KS amplitudes required')
        out = np.empty(4)
        for B in range(4):
            shifted = base.astype(complex); shifted[B] += 1j*step
            out[B] = self.evaluate(history, shifted)['action'].imag/step
        return out
