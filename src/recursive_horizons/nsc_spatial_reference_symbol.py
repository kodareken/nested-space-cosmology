"""Fourth-order formal reference symbol of the owned spherical Dirac operator.

Uses the Weyl-product projector construction (Panati--Spohn--Teufel,
math-ph/0201055, Secs.2 and4.4), not a new state prescription. The symbol is
not a pointwise Gaussian covariance or a complete finite counterfunctional.
Taylor variables are the actual KS coordinate time, axial coordinate, and
canonical momentum. The formal order counter is not a physical parameter.
"""
from functools import lru_cache
from math import comb, factorial
import numpy as np

from .nsc_ks_spacetime_variation import chart_coordinates

I = np.eye(2, dtype=complex)
SIGMA = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]], complex)
INDICES = tuple((t, z, k) for t in range(5) for z in range(5-t) for k in range(5))
INDEX = {v: i for i, v in enumerate(INDICES)}
ZERO = INDEX[0, 0, 0]


@lru_cache(maxsize=None)
def convolution(physical, momentum):
    left, right, target = [], [], []
    for i, a in enumerate(INDICES):
        for j, b in enumerate(INDICES):
            c = tuple(x+y for x, y in zip(a, b))
            if c[0]+c[1] <= physical and c[2] <= momentum:
                left.append(i); right.append(j); target.append(INDEX[c])
    return np.array(left), np.array(right), np.array(target)


class SymbolJet:
    """Finite Taylor algebra for this order-four matrix-symbol calculation."""
    __array_priority__ = 1000

    def __init__(self, data, physical=4, momentum=4):
        self.data = np.asarray(data, complex)
        if self.data.ndim != 4 or self.data.shape[0] != len(INDICES) or self.data.shape[-2:] != (2, 2):
            raise ValueError('Taylor coefficients of 2x2 matrix symbols required')
        self.physical, self.momentum = physical, momentum
        self.data = self.data.copy()
        for i, (t, z, k) in enumerate(INDICES):
            if t+z > physical or k > momentum: self.data[i] = 0.

    @classmethod
    def constant(cls, value):
        v = np.asarray(value, complex)
        if v.ndim == 0: v = v*I
        elif v.ndim == 1: v = v[:, None, None]*I
        if v.ndim == 2: v = v[None]
        out = np.zeros((len(INDICES), *v.shape), complex); out[ZERO] = v
        return cls(out)

    @classmethod
    def variable(cls, value, axis):
        out = cls.constant(value)
        index = [0, 0, 0]; index[axis] = 1
        out.data[INDEX[tuple(index)]] = I
        return out

    @property
    def value(self): return self.data[ZERO]

    def trim(self, physical, momentum):
        if physical > self.physical or momentum > self.momentum:
            raise ValueError('cannot invent uncomputed Taylor derivatives')
        return SymbolJet(self.data, physical, momentum)

    def __add__(self, other):
        b = other if isinstance(other, SymbolJet) else self.constant(other)
        return SymbolJet(self.data+b.data, min(self.physical, b.physical), min(self.momentum, b.momentum))
    __radd__ = __add__
    def __neg__(self): return SymbolJet(-self.data, self.physical, self.momentum)
    def __sub__(self, other): return self+-other
    def __rsub__(self, other): return -self+other

    def __mul__(self, other):
        if not isinstance(other, SymbolJet):
            a = np.asarray(other)
            if a.ndim == 0: return SymbolJet(self.data*a, self.physical, self.momentum)
            return self*self.constant(other)
        pt, pk = min(self.physical, other.physical), min(self.momentum, other.momentum)
        i, j, k = convolution(pt, pk)
        product = self.data[i]@other.data[j]
        out = np.zeros((len(INDICES), *product.shape[1:]), complex)
        np.add.at(out, k, product)
        return SymbolJet(out, pt, pk)
    def __rmul__(self, other):
        if np.asarray(other).ndim == 0: return self*other
        return self.constant(other)*self
    def __truediv__(self, other):
        return self*(other.power(-1) if isinstance(other, SymbolJet) else 1/other)
    def __rtruediv__(self, other): return self.constant(other)*self.power(-1)

    def right(self, matrix):
        return SymbolJet(self.data@matrix, self.physical, self.momentum)

    def scalar_function(self, coefficients):
        value = self.value[:, 0, 0]
        if np.max(abs(self.data-self.data[..., :1, :1]*I)) > 3e-11:
            raise ValueError('scalar functional calculus only; not a matrix spectral shortcut')
        increment = self-self.constant(value)
        out = self.constant(coefficients(value, 0)).trim(self.physical, self.momentum)
        term = self.constant(1).trim(self.physical, self.momentum)
        for n in range(1, self.physical+self.momentum+1):
            term = term*increment
            out += term*self.constant(coefficients(value, n))
        return out

    def power(self, power):
        def coefficient(x, n):
            factor = 1.
            for k in range(n): factor *= (power-k)/(k+1)
            return factor*x**(power-n)
        if np.any(abs(self.value[:, 0, 0]) < 1e-14) and power < 1:
            raise ValueError('gapped nonzero scalar branch required; no denominator floor')
        return self.scalar_function(coefficient)

    def exp(self): return self.scalar_function(lambda x, n: np.exp(x)/factorial(n))
    def log(self):
        return self.scalar_function(lambda x, n: np.log(x) if n == 0 else (-1)**(n+1)/(n*x**n))
    def sin(self): return ((1j*self).exp()-(-1j*self).exp())/(2j)
    def cos(self): return ((1j*self).exp()+(-1j*self).exp())/2
    def atan(self): return ((1+1j*self).log()-(1-1j*self).log())/(2j)

    def derivative(self, t=0, z=0, k=0):
        pt, pk = self.physical-t-z, self.momentum-k
        if min(pt, pk) < 0: raise ValueError('requested derivative exceeds the owned Taylor order')
        out = np.zeros_like(self.data)
        for i, a in enumerate(INDICES):
            if a[0]+a[1] > pt or a[2] > pk: continue
            b = (a[0]+t, a[1]+z, a[2]+k)
            factor = np.prod([factorial(v)/factorial(u) for u, v in zip(a, b)])
            out[i] = factor*self.data[INDEX[b]]
        return SymbolJet(out, pt, pk)


def star_order(A, B, order):
    """Coefficient of epsilon^order in the spatial Weyl product A star B."""
    result = SymbolJet.constant(0)
    for j in range(order+1):
        term = A.derivative(z=order-j, k=j)*B.derivative(z=j, k=order-j)
        result += ((1j/2)**order/factorial(order))*(-1)**j*comb(order, j)*term
    return result


def reference_projector(N, beta, axial, radius, momenta, masses, angular):
    """Solve projector and transport identities through formal order four.

    Uses coordinate d_T with the lapse-weighted H. Shift-gradient momentum
    transport is supplied by the Weyl products, not discarded by parity.
    """
    for field in (N, axial, radius):
        value = field.value[:, 0, 0]
        if np.any(value.real <= 0) or np.max(abs(value.imag)) > 3e-11:
            raise ValueError('positive real lapse and spatial scale factors required')
    k = SymbolJet.variable(np.asarray(momenta), 2)
    h = [-N*SymbolJet.constant(np.asarray(masses)),
         N*SymbolJet.constant(np.asarray(angular))/radius, N*k/axial]
    gap2 = sum(component*component for component in h)
    gap = gap2.power(.5)
    if np.any(gap.value[:, 0, 0].real <= 1e-12): raise ValueError('nonzero normal-energy gap required')
    spin = sum(component.right(sigma) for component, sigma in zip(h, SIGMA))
    H = spin-beta*k
    P = [(1-spin/gap)/2]; Q = 1-P[0]
    diagnostics = []
    for n in range(1, 5):
        G = SymbolJet.constant(0)
        for i in range(n):
            for j in range(n):
                ell = n-i-j
                if ell >= 0: G += star_order(P[i], P[j], ell)
        F = 1j*P[n-1].derivative(t=1)
        for j in range(n):
            F -= star_order(H, P[j], n-j)-star_order(P[j], H, n-j)
        diagonal = -P[0]*G*P[0]+Q*G*Q
        off = (H*F-F*H)/(4*gap2)
        pn = (diagonal+off).trim(4-n, 4-n)
        idem = P[0]*pn+pn*P[0]-pn+G
        transport = H*pn-pn*H-F
        compatibility_g = P[0]*G*Q+Q*G*P[0]
        compatibility_f = P[0]*F*P[0]+Q*F*Q
        hermitian = pn.value-pn.value.swapaxes(-1, -2).conj()
        row = {}
        for name, error, scale in (
            ('idempotence', idem.value, 1+np.max(abs(G.value), axis=(-2, -1))+2*np.max(abs(pn.value), axis=(-2, -1))),
            ('transport', transport.value, 1+np.max(abs(F.value), axis=(-2, -1))+2*np.max(abs((H*pn).value), axis=(-2, -1))),
            ('G_off_block', compatibility_g.value, 1+np.max(abs(G.value), axis=(-2, -1))),
            ('F_diagonal_block', compatibility_f.value, 1+np.max(abs(F.value), axis=(-2, -1))),
            ('Hermiticity', hermitian, 1+np.max(abs(pn.value), axis=(-2, -1))),
        ):
            raw = np.max(abs(error), axis=(-2, -1))
            row[name] = {'absolute': raw, 'scaled': raw/scale}
        diagnostics.append(row); P.append(pn)
    return {'orders': np.array([p.value for p in P]), 'jets': P, 'residuals': diagnostics,
            'gap': gap.value[:, 0, 0].real, 'hamiltonian': H,
            'branch': 'lower normal-energy band, not sign of shifted coordinate H'}


def background_scalar(rho):
    beta = (3*((1+rho*rho)*(np.pi/2-rho.atan())-rho)).power(.5)
    axial = (beta*beta-1).power(.5)
    return beta, axial, (1+rho*rho).power(.5)


def supplied_KS_metric_jets(tau, rho, amplitudes=(.002, .001, .003, .001), omega=.4, support=(0., .05)):
    """Same PG response history expressed in true KS coordinates (T,z)."""
    if not -1 < rho < 1 or not support[0] < tau < support[1]:
        raise ValueError('interior of the owned compact response required')
    x = SymbolJet.constant(rho)
    S0 = chart_coordinates(rho)[1]; S = SymbolJet.constant(S0)
    # d rho/dT=-a0 and dS/dT=-beta0/a0 in the fixed reference chart.
    for n in range(4):
        beta0, a0, _ = background_scalar(x)
        x.data[INDEX[n+1, 0, 0]] = -a0.data[INDEX[n, 0, 0]]/(n+1)
        S.data[INDEX[n+1, 0, 0]] = -(beta0/a0).data[INDEX[n, 0, 0]]/(n+1)
    z = SymbolJet.variable(tau+S0, 1)
    time = z-S
    u = 2*(time-support[0])/(support[1]-support[0])-1
    bump = (1-1/(1-u*u)).exp()*(1-1/(1-x*x)).exp()*(omega*z).cos()
    _, axial0, radius0 = background_scalar(x)
    A = np.asarray(amplitudes, float)
    if A.shape != (4,) or not np.isfinite(A).all(): raise ValueError('four existing raw-KS amplitudes required')
    return 1+A[0]*bump, A[1]*bump, axial0+A[2]*bump, radius0+A[3]*bump


def homogeneous_seed_metric_jets():
    """Frozen compact reference control, in its owned conformal clock.

    qdot=-sqrt(W), N=r, a=r*sqrt(W). No old generator is executed.
    """
    q = SymbolJet.constant(np.pi/2)
    def W(q): return 3*(np.pi-q)+1.5*(2*q).sin()-q.sin()*q.sin()
    for n in range(4):
        q.data[INDEX[n+1, 0, 0]] = -W(q).power(.5).data[INDEX[n, 0, 0]]/(n+1)
    radius = 1/q.sin(); axial = radius*W(q).power(.5)
    return radius, SymbolJet.constant(0), axial, radius
