"""Locked local NSC action on the transmitting spherical spacetime history.

This applies the same curvature convention and coefficients as finite_terms
and GeneralKSLocalInducedHistory. It extends evaluation to both PG coordinates;
it does not derive new field equations, select a history, or subtract a vacuum.
Only boundary-identical compact directions are admitted by the control owner.
"""
from dataclasses import dataclass
import numpy as np
from numpy.polynomial.legendre import leggauss

from .nsc_ks_spacetime_variation import chart_coordinates
from .nsc_lorentzian import geometry
from .nsc_transmitting_resolvent import profile


def gauss_axis(count, left, right):
    """Polynomial differentiation and Gauss integration on one coordinate."""
    if count < 8 or not left < right:
        raise ValueError('resolved ordered coordinate interval required')
    x, w = leggauss(count)
    bary = (-1.)**np.arange(count)*np.sqrt((1-x*x)*w)
    distance = x[:, None]-x[None, :]
    np.fill_diagonal(distance, 1.)
    D = (bary[None, :]/bary[:, None])/distance
    np.fill_diagonal(D, 0.)
    np.fill_diagonal(D, -D.sum(axis=1))
    return left+(right-left)*(x+1)/2, (right-left)*w/2, 2*D/(right-left)


@dataclass
class SphericalHistoryGrid:
    time: np.ndarray
    radius_coordinate: np.ndarray
    time_weights: np.ndarray
    radial_weights: np.ndarray
    time_D: np.ndarray
    radial_D: np.ndarray

    @classmethod
    def gaussian(cls, nt, nx, time_interval=(0., .05), radial_interval=(-1., 1.)):
        t, wt, Dt = gauss_axis(nt, *time_interval)
        x, wx, Dx = gauss_axis(nx, *radial_interval)
        return cls(t, x, wt, wx, Dt, Dx)

    def derivative(self, value, axis):
        value = np.asarray(value)
        if value.shape[:2] != (len(self.time), len(self.radius_coordinate)):
            raise ValueError('field must live on this two-coordinate history grid')
        # An exactly constant coordinate dependence has exactly zero derivative.
        if axis == 0:
            return np.einsum('ij,jk...->ik...', self.time_D, value-value[:1], optimize=True)
        if axis == 1:
            return np.einsum('ij,kj...->ki...', self.radial_D, value-value[:, :1], optimize=True)
        raise ValueError('only the two base coordinates are differentiated')

    def integral(self, value):
        return np.einsum('t,x,tx...->...', self.time_weights, self.radial_weights, value, optimize=True)


def spherical_invariants(grid, g, radius, jets=None):
    """Apply the owned (+---) curvature contractions to a warped S2 metric.

    g_ab is the general two-coordinate Lorentzian base. The formulas reduce
    to finite_terms._symbolic_basis and the homogeneous KS local owner.
    No field equation or gravitational source is inferred by this operation.
    """
    g, r = np.asarray(g), np.asarray(radius)
    shape = (len(grid.time), len(grid.radius_coordinate))
    if g.shape != (*shape, 2, 2) or r.shape != shape or not np.isfinite(g).all() or not np.isfinite(r).all():
        raise ValueError('finite base metric and areal radius on the full history required')
    det = g[..., 0, 0]*g[..., 1, 1]-g[..., 0, 1]*g[..., 1, 0]
    if np.max(abs(g[..., 0, 1]-g[..., 1, 0])) > 3e-11 or np.any(det.real >= 0) or np.any(r.real <= 0):
        raise ValueError('symmetric Lorentzian base and positive areal radius required')
    inv = np.empty_like(g)
    inv[..., 0, 0], inv[..., 1, 1] = g[..., 1, 1]/det, g[..., 0, 0]/det
    inv[..., 0, 1], inv[..., 1, 0] = -g[..., 0, 1]/det, -g[..., 1, 0]/det
    if jets is None:
        dg = np.stack([grid.derivative(g, a) for a in range(2)], axis=2)
        ddg = np.stack([grid.derivative(dg, a) for a in range(2)], axis=2)
        dr = np.stack([grid.derivative(r, a) for a in range(2)], axis=2)
        ddr = np.stack([grid.derivative(dr, a) for a in range(2)], axis=2)
    else:
        dg, ddg, dr, ddr = jets
    di = -np.einsum('...ij,...ajk,...kl->...ail', inv, dg, inv, optimize=True)
    gamma = np.zeros((*shape, 2, 2, 2), dtype=g.dtype)
    dgamma = np.zeros((*shape, 2, 2, 2, 2), dtype=g.dtype)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                for d in range(2):
                    P = dg[..., b, d, c]+dg[..., c, d, b]-dg[..., d, b, c]
                    gamma[..., a, b, c] += .5*inv[..., a, d]*P
                    for e in range(2):
                        dP = ddg[..., e, b, d, c]+ddg[..., e, c, d, b]-ddg[..., e, d, b, c]
                        dgamma[..., e, a, b, c] += .5*(di[..., e, a, d]*P+inv[..., a, d]*dP)
    ricci = np.zeros_like(g)
    for b in range(2):
        for d in range(2):
            for a in range(2):
                ricci[..., b, d] += dgamma[..., a, a, d, b]-dgamma[..., d, a, a, b]
                for e in range(2):
                    ricci[..., b, d] += (gamma[..., a, a, e]*gamma[..., e, d, b]
                                        -gamma[..., a, d, e]*gamma[..., e, a, b])
    R2 = np.einsum('...ab,...ab->...', inv, ricci)
    hess = ddr-np.einsum('...cab,...c->...ab', gamma, dr)
    boxr = np.einsum('...ab,...ab->...', inv, hess)
    grad2 = np.einsum('...ab,...a,...b->...', inv, dr, dr)
    hess2 = np.einsum('...ac,...bd,...ab,...cd->...', inv, inv, hess, hess, optimize=True)
    R = R2-4*boxr/r-2*(1+grad2)/r**2
    Ricci2 = R2**2/2-2*R2*boxr/r+4*hess2/r**2+2*(1+grad2+r*boxr)**2/r**4
    Riemann2 = R2**2+8*hess2/r**2+4*(1+grad2)**2/r**4
    C2 = Riemann2-2*Ricci2+R**2/3
    E4 = Riemann2-4*Ricci2+R**2
    measure = np.sqrt(-det)*r*r
    return {'R': R, 'C2': C2, 'E4': E4, 'measure': measure,
            'R2': R2, 'Ricci2': Ricci2, 'Riemann2': Riemann2,
            'base_ricci_identity': ricci-.5*R2[..., None, None]*g}


@dataclass
class LockedSphericalLocalAction:
    ledger: dict

    def __post_init__(self):
        required = ('A', 'C_gauge', 'C_Weyl', 'C_Euler', 'C_boxR', 'magnetic_flux', 'V_full')
        if not all(k in self.ledger and np.isfinite(self.ledger[k]) for k in required):
            raise ValueError('complete finite declared ledger required')
        if self.ledger['V_full'] != 0:
            raise ValueError('this gate imports the locked V_full=0 branch')

    def actions(self, grid, g, radius, jets=None):
        c = spherical_invariants(grid, g, radius, jets)
        l, mu = self.ledger, c['measure']
        # Existing F_{theta phi}=q_mag sin(theta)/2 convention.
        values = {
            'einstein_bulk': -4*np.pi*l['A']*grid.integral(mu*c['R']),
            'maxwell_bulk': -2*np.pi*l['C_gauge']*l['magnetic_flux']**2*grid.integral(mu/radius**4),
            'weyl_bulk': -4*np.pi*l['C_Weyl']*grid.integral(mu*c['C2']),
            'euler_bulk_diagnostic': -4*np.pi*l['C_Euler']*grid.integral(mu*c['E4']),
        }
        return values, c


class OwnedCompactKSHistory:
    """Same supplied history as the field owner, now holomorphic in amplitudes.

    The fixed chart is pulled back at metric-tensor level. This continuation
    is used for action differentiation only; complex metrics are not physical
    histories. At real amplitudes it must agree with SuppliedKSHarmonicMetric.
    """
    def __init__(self, grid, omega=.4, support_time=(0., .05)):
        self.grid, self.omega, self.support_time = grid, float(omega), tuple(support_time)
        x, t = grid.radius_coordinate, grid.time
        self.beta = geometry(x)[0]
        if np.any(self.beta <= 1): raise ValueError('owned trapped chart required')
        self.axial = np.sqrt(self.beta**2-1)
        self.radius = np.sqrt(1+x*x)
        self.shift_integral = np.array([chart_coordinates(float(rho))[1] for rho in x])
        lo, hi = support_time
        if (not hi > lo or abs(grid.time_weights.sum()-(hi-lo)) > 3e-12
                or abs(grid.radial_weights.sum()-2.) > 3e-12
                or np.any((x <= -1) | (x >= 1)) or np.any((t <= lo) | (t >= hi))):
            raise ValueError('integrate the complete compact support, with fixed external metric jets')
        wt = np.array([profile(2*(tau-lo)/(hi-lo)-1)[0] for tau in t])
        wx = np.array([profile(rho)[0] for rho in x])
        self.envelope = wt[:, None]*wx[None, :]*np.cos(omega*(t[:, None]+self.shift_integral[None, :]))
        self.chart = np.zeros((len(x), 2, 2))
        self.chart[:, 0, 1] = -1/self.axial
        self.chart[:, 1, 0] = 1.
        self.chart[:, 1, 1] = self.beta/self.axial**2

    def fields(self, amplitudes):
        a = np.asarray(amplitudes)
        if a.shape != (4,) or not np.isfinite(a).all(): raise ValueError('four finite raw KS amplitudes required')
        N = 1+a[0]*self.envelope
        beta = a[1]*self.envelope
        axial = self.axial[None, :]+a[2]*self.envelope
        r = self.radius[None, :]+a[3]*self.envelope
        if min(N.real.min(), axial.real.min(), r.real.min()) <= 0: raise ValueError('positive raw KS metric required')
        gK = np.empty((*N.shape, 2, 2), dtype=np.result_type(a, float))
        gK[..., 0, 0] = N*N-axial*axial*beta*beta
        gK[..., 0, 1] = gK[..., 1, 0] = -axial*axial*beta
        gK[..., 1, 1] = -axial*axial
        gP = np.einsum('xai,txab,xbj->txij', self.chart, gK, self.chart, optimize=True)
        return gP, r

    def fields_and_jets(self, amplitudes):
        """Exact coordinate jets of this existing metric family through order2.

        Finite Taylor algebra carries products and the fixed chart. Numerical
        coordinate differentiation would amplify flat-endpoint roundoff and
        is not used to assign spurious Euler/BoxR boundary forces.
        """
        a = np.asarray(amplitudes)
        if a.shape != (4,) or not np.isfinite(a).all(): raise ValueError('four raw KS amplitudes required')
        x, t = self.grid.radius_coordinate, self.grid.time
        shape = self.envelope.shape
        indices = ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2))
        def const(v):
            z = np.zeros((6, *shape), dtype=np.result_type(a, float)); z[0] = v
            return z
        def mul(u, v):
            z = np.zeros_like(u+v)
            for k, (kt, kx) in enumerate(indices):
                for i, (it, ix) in enumerate(indices):
                    for j, (jt, jx) in enumerate(indices):
                        if it+jt == kt and ix+jx == kx: z[k] += u[i]*v[j]
            return z
        lo, hi = self.support_time
        u = 2*(t-lo)/(hi-lo)-1; d = 1-u*u
        w = np.zeros_like(u); wt = w.copy(); wtt = w.copy()
        inside = d > 0
        w[inside] = np.exp(1-1/d[inside])
        l1 = -2*u[inside]/d[inside]**2*2/(hi-lo)
        l2 = (-2/d[inside]**2-8*u[inside]**2/d[inside]**3)*(2/(hi-lo))**2
        wt[inside], wtt[inside] = w[inside]*l1, w[inside]*(l2+l1*l1)
        e = 1-x*x; s = np.zeros_like(x); sx = s.copy(); sxx = s.copy()
        inside = e > 0
        s[inside] = np.exp(1-1/e[inside])
        l1 = -2*x[inside]/e[inside]**2
        l2 = -2/e[inside]**2-8*x[inside]**2/e[inside]**3
        sx[inside], sxx[inside] = s[inside]*l1, s[inside]*(l2+l1*l1)
        b, bp, _ = geometry(x)
        angle = np.pi/2-np.arctan(x)
        bpp = 3*(angle-x/(1+x*x))/b-bp*bp/b
        ax = self.axial; axp = b*bp/ax
        axpp = (bp*bp+b*bpp)/ax-(b*bp)**2/ax**3
        sp = b/ax**2; spp = bp/ax**2-2*b*axp/ax**3
        phase = self.omega*(t[:, None]+self.shift_integral[None, :])
        co, si, om = np.cos(phase), np.sin(phase), self.omega
        f = const(self.envelope)
        f[1] = s*(wt[:, None]*co-om*w[:, None]*si)
        f[2] = w[:, None]*(sx*co-om*s*sp*si)
        f[3] = s*(wtt[:, None]*co-2*om*wt[:, None]*si-om*om*w[:, None]*co)/2
        f[4] = ((wt[:, None]*sx-om*om*w[:, None]*s*sp)*co
                -om*(w[:, None]*sx+wt[:, None]*s*sp)*si)
        f[5] = w[:, None]*((sxx-om*om*s*sp*sp)*co-om*(2*sx*sp+s*spp)*si)/2
        N = const(1)+a[0]*f; beta = a[1]*f
        axial = const(ax); axial[2] = axp; axial[5] = axpp/2; axial += a[2]*f
        r = const(self.radius); r[2] = x/self.radius; r[5] = 1/(2*self.radius**3); r += a[3]*f
        ax2 = mul(axial, axial)
        gK = [[mul(N, N)-mul(ax2, mul(beta, beta)), -mul(ax2, beta)],
              [-mul(ax2, beta), -ax2]]
        J01 = const(-1/ax); J01[2] = axp/ax**2; J01[5] = (axpp/ax**2-2*axp*axp/ax**3)/2
        J11 = const(sp); J11[2] = spp
        J11[5] = (bpp/ax**2-4*bp*axp/ax**3-2*b*axpp/ax**3+6*b*axp*axp/ax**4)/2
        J = [[const(0), J01], [const(1), J11]]
        G = np.zeros((6, *shape, 2, 2), dtype=np.result_type(a, float))
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for l in range(2): G[..., i, j] += mul(mul(J[k][i], gK[k][l]), J[l][j])
        dg = np.stack((G[1], G[2]), axis=2)
        ddg = np.stack((np.stack((2*G[3], G[4]), axis=2), np.stack((G[4], 2*G[5]), axis=2)), axis=2)
        dr = np.stack((r[1], r[2]), axis=2)
        ddr = np.stack((np.stack((2*r[3], r[4]), axis=2), np.stack((r[4], 2*r[5]), axis=2)), axis=2)
        return G[0], r[0], (dg, ddg, dr, ddr)

    def action_variation(self, action, amplitudes, step=1e-28):
        base = np.asarray(amplitudes, float)
        if not np.isfinite(step) or step <= 0: raise ValueError('positive numerical complex step required')
        result = {}
        for B in range(4):
            varied = base.astype(complex); varied[B] += 1j*step
            rows, _ = action.actions(self.grid, *self.fields_and_jets(varied))
            for name, value in rows.items():
                result.setdefault(name, np.zeros(4))[B] = value.imag/step
        return result

    def compact_variation(self, action, amplitudes):
        """Compose the allowed local derivative, without adding endpoint terms.

        GHY/BoxR endpoint variations are zero because this explicit family
        retains all external metric jets. Euler is checked by integrating its
        bulk density but is likewise not assigned a spurious bulk force.
        """
        gradients = self.action_variation(action, amplitudes)
        total = sum(gradients[name] for name in ('einstein_bulk', 'maxwell_bulk', 'weyl_bulk'))
        return {'channel_gradients': gradients, 'local_action_gradient': total,
                'local_action_force': -total,
                'fixed_boundary_derivatives': {'einstein_GHY': 0., 'euler': 0., 'boxR': 0.},
                'Euler_bulk_residual': float(np.max(abs(gradients['euler_bulk_diagnostic']))),
                'boundary_reason': 'same fixed domain and metric jets; smooth transmitting cut has opposite orientations'}
