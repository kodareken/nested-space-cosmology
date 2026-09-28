"""Regular high-momentum jets of the owned fourth-order auxiliary projector.

The recurrence is Panati--Spohn--Teufel (math-ph/0201055, Secs. 2 and 4.4) as
already executed by ``nsc_spatial_reference_symbol``. New work is the regular
representation in ``mu=1/|k|`` on the same KS ``(T,z)`` geometry, so intervals
that meet ``mu=0`` keep the true high-momentum decay instead of cancelling
``P_0`` against ``P_inf``. These symbols are auxiliary: they are not ``C_up``
or ``C_Sigma``, they do not replace the physical state, and they do not add a
``Gamma`` term. The exact split remains ``C_g=P_g+D_g``.
"""
from dataclasses import dataclass
from functools import lru_cache
from math import comb, factorial

from flint import acb, arb, ctx

from .nsc_ks_spacetime_geometry_jets import KSGeometryJets

DEFAULT_PHYSICAL = 9
DEFAULT_MOMENTUM = 5
MAX_PHYSICAL = 16
MAX_MOMENTUM = 8
FORMAL_ORDER = 4

_I = (1, 0, 0, 1)
_S1 = (0, 1, 1, 0)
_S2 = (0, -1, 1, 0)  # times i; applied as acb(0,1)*_S2
_S3 = (1, 0, 0, -1)


def _as_int(value, name, minimum, maximum):
    if isinstance(value, bool) or int(value) != value:
        raise ValueError(f'integer {name} required')
    value = int(value)
    if not minimum <= value <= maximum:
        raise ValueError(f'{name} must lie in [{minimum}, {maximum}]')
    return value


def _real_ball(value, name):
    if callable(value) or isinstance(value, (bool, str, complex, acb)):
        raise TypeError(f'explicit real ball {name} required; callables and shapes are rejected')
    if isinstance(value, arb):
        ball = value
    elif isinstance(value, int):
        ball = arb(value)
    elif isinstance(value, float):
        ball = arb(value)
    else:
        raise TypeError(f'explicit real ball {name} required; callables and shapes are rejected')
    if not ball.is_finite():
        raise ValueError(f'finite {name} required')
    return ball


def positive_inv(value, name='denominator'):
    """Directed inverse of a strictly positive enclosing interval.

    Uses the ball itself. Rounding an upper bound before inversion is invalid.
    """
    if not isinstance(value, arb):
        value = _real_ball(value, name)
    if not value > 0:
        raise ValueError(f'positive enclosing interval required for {name}; no denominator floor')
    return 1 / value


def _nonneg_mu(value):
    mu = _real_ball(value, 'mu')
    if not mu >= 0:
        raise ValueError('mu=1/|k| must be a nonnegative enclosing interval')
    return mu


@lru_cache(maxsize=None)
def _table(physical, momentum):
    indices = tuple((t, z, m) for t in range(physical + 1)
                    for z in range(physical + 1 - t)
                    for m in range(momentum + 1))
    return indices, {key: i for i, key in enumerate(indices)}


def _zero_mat():
    return [acb(0), acb(0), acb(0), acb(0)]


def _mat(values):
    return [acb(values[0]), acb(values[1]), acb(values[2]), acb(values[3])]


def _matmul(a, b):
    return [a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3],
            a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3]]


def _is_zero_mat(a):
    return a[0].is_zero() and a[1].is_zero() and a[2].is_zero() and a[3].is_zero()


def _pochhammer(start, count):
    factor = 1
    for i in range(count):
        factor *= start + i
    return factor


class ScalarJet:
    """Triangular mixed Taylor jet in (T,z,mu); coefficients are d^{t,z,m} f/(t! z! m!)."""

    def __init__(self, data, physical, momentum, indices=None, index=None):
        self.physical = physical
        self.momentum = momentum
        self.indices = indices if indices is not None else _table(physical, momentum)[0]
        self.index = index if index is not None else _table(physical, momentum)[1]
        if len(data) != len(self.indices):
            raise ValueError('scalar jet coefficient table does not match the mixed (T,z,mu) index')
        self.data = list(data)

    @classmethod
    def zeros(cls, physical, momentum):
        indices, index = _table(physical, momentum)
        return cls([arb(0) for _ in indices], physical, momentum, indices, index)

    @classmethod
    def constant(cls, value, physical, momentum):
        out = cls.zeros(physical, momentum)
        out.data[0] = _real_ball(value, 'constant') if not isinstance(value, arb) else value
        return out

    @classmethod
    def variable(cls, value, axis, physical, momentum):
        if axis not in (0, 1, 2):
            raise ValueError('axis must be 0=T, 1=z or 2=mu')
        out = cls.constant(value, physical, momentum)
        key = [0, 0, 0]
        key[axis] = 1
        if axis < 2 and 1 > physical:
            raise ValueError('requested coordinate exceeds the owned physical order')
        if axis == 2 and 1 > momentum:
            raise ValueError('requested coordinate exceeds the owned mu order')
        out.data[out.index[tuple(key)]] = arb(1)
        return out

    @property
    def value(self):
        return self.data[0]

    def broadcast(self, value):
        return ScalarJet.constant(value, self.physical, self.momentum)

    def coefficient(self, t=0, z=0, mu=0):
        key = (int(t), int(z), int(mu))
        if key not in self.index:
            raise ValueError('requested mixed Taylor coefficient exceeds the owned jet')
        return self.data[self.index[key]]

    def trim(self, physical, momentum):
        if physical > self.physical or momentum > self.momentum:
            raise ValueError('cannot invent uncomputed Taylor derivatives')
        if physical == self.physical and momentum == self.momentum:
            return ScalarJet(self.data, physical, momentum, self.indices, self.index)
        indices, index = _table(physical, momentum)
        data = [self.data[self.index[key]] if key in self.index else arb(0) for key in indices]
        return ScalarJet(data, physical, momentum, indices, index)

    def derivative(self, t=0, z=0, mu=0):
        t, z, mu = int(t), int(z), int(mu)
        if min(t, z, mu) < 0:
            raise ValueError('nonnegative mixed derivative required')
        pt, pm = self.physical - t - z, self.momentum - mu
        if min(pt, pm) < 0:
            raise ValueError('requested derivative exceeds the owned Taylor order')
        indices, index = _table(pt, pm)
        data = []
        for t0, z0, m0 in indices:
            src = (t0 + t, z0 + z, m0 + mu)
            factor = _pochhammer(t0 + 1, t) * _pochhammer(z0 + 1, z) * _pochhammer(m0 + 1, mu)
            data.append(self.data[self.index[src]] * factor)
        return ScalarJet(data, pt, pm, indices, index)

    def __add__(self, other):
        b = other if isinstance(other, ScalarJet) else self.broadcast(other)
        pt, pm = min(self.physical, b.physical), min(self.momentum, b.momentum)
        left, right = self.trim(pt, pm), b.trim(pt, pm)
        return ScalarJet([x + y for x, y in zip(left.data, right.data)], pt, pm, left.indices, left.index)

    __radd__ = __add__

    def __neg__(self):
        return ScalarJet([-x for x in self.data], self.physical, self.momentum, self.indices, self.index)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def __mul__(self, other):
        if not isinstance(other, ScalarJet):
            if isinstance(other, MatrixJet):
                return other * self
            scale = other if isinstance(other, arb) else _real_ball(other, 'scale')
            return ScalarJet([scale * x for x in self.data], self.physical, self.momentum, self.indices, self.index)
        pt, pm = min(self.physical, other.physical), min(self.momentum, other.momentum)
        indices, index = _table(pt, pm)
        out = [arb(0) for _ in indices]
        for i, (t1, z1, m1) in enumerate(self.indices):
            va = self.data[i]
            if va.is_zero():
                continue
            for j, (t2, z2, m2) in enumerate(other.indices):
                t, z, m = t1 + t2, z1 + z2, m1 + m2
                if t + z <= pt and m <= pm:
                    vb = other.data[j]
                    if vb.is_zero():
                        continue
                    out[index[t, z, m]] += va * vb
        return ScalarJet(out, pt, pm, indices, index)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, ScalarJet):
            return self * other.inv()
        return self * positive_inv(other if isinstance(other, arb) else _real_ball(other, 'denominator'))

    def _compose(self, coefficient):
        value = self.value
        out = self.broadcast(coefficient(value, 0))
        increment = self - self.broadcast(value)
        increment.data[0] = arb(0)
        term = self.broadcast(arb(1))
        for n in range(1, self.physical + self.momentum + 1):
            term = term * increment
            out += term * self.broadcast(coefficient(value, n))
        return out

    def inv(self):
        def coefficient(x, n):
            inv = positive_inv(x)
            if n == 0:
                return inv
            return (arb(-1) if n % 2 else arb(1)) * inv ** (n + 1)
        return self._compose(coefficient)

    def sqrt(self):
        def coefficient(x, n):
            if not x > 0:
                raise ValueError('positive enclosing interval required for square root; no denominator floor')
            if n == 0:
                return x.sqrt()
            factor = arb(1)
            half = arb(1) / 2
            for k in range(n):
                factor *= (half - k) / (k + 1)
            return factor * x.sqrt() / x ** n
        return self._compose(coefficient)


class MatrixJet:
    """2x2 acb companion of :class:`ScalarJet`, row-major (00,01,10,11)."""

    def __init__(self, data, physical, momentum, indices=None, index=None):
        self.physical = physical
        self.momentum = momentum
        self.indices = indices if indices is not None else _table(physical, momentum)[0]
        self.index = index if index is not None else _table(physical, momentum)[1]
        if len(data) != len(self.indices):
            raise ValueError('matrix jet coefficient table does not match the mixed (T,z,mu) index')
        self.data = data

    @classmethod
    def zeros(cls, physical, momentum):
        indices, index = _table(physical, momentum)
        return cls([_zero_mat() for _ in indices], physical, momentum, indices, index)

    @classmethod
    def constant(cls, matrix, physical, momentum):
        out = cls.zeros(physical, momentum)
        out.data[0] = [acb(matrix[0]), acb(matrix[1]), acb(matrix[2]), acb(matrix[3])]
        return out

    @classmethod
    def pauli(cls, scalar, sigma, physical, momentum):
        base = scalar if isinstance(scalar, ScalarJet) else ScalarJet.constant(scalar, physical, momentum)
        data = [_zero_mat() for _ in base.indices]
        for i, coeff in enumerate(base.data):
            if coeff.is_zero():
                continue
            data[i] = [coeff * sigma[0], coeff * sigma[1], coeff * sigma[2], coeff * sigma[3]]
        return cls(data, base.physical, base.momentum, base.indices, base.index)

    @property
    def value(self):
        return self.data[0]

    def coefficient(self, t=0, z=0, mu=0):
        key = (int(t), int(z), int(mu))
        if key not in self.index:
            raise ValueError('requested mixed Taylor coefficient exceeds the owned jet')
        return self.data[self.index[key]]

    def trace_value(self):
        v = self.value
        return v[0] + v[3]

    def trim(self, physical, momentum):
        if physical > self.physical or momentum > self.momentum:
            raise ValueError('cannot invent uncomputed Taylor derivatives')
        if physical == self.physical and momentum == self.momentum:
            return MatrixJet(self.data, physical, momentum, self.indices, self.index)
        indices, index = _table(physical, momentum)
        data = [self.data[self.index[key]] if key in self.index else _zero_mat() for key in indices]
        return MatrixJet(data, physical, momentum, indices, index)

    def derivative(self, t=0, z=0, mu=0):
        t, z, mu = int(t), int(z), int(mu)
        if min(t, z, mu) < 0:
            raise ValueError('nonnegative mixed derivative required')
        pt, pm = self.physical - t - z, self.momentum - mu
        if min(pt, pm) < 0:
            raise ValueError('requested derivative exceeds the owned Taylor order')
        indices, index = _table(pt, pm)
        data = []
        for t0, z0, m0 in indices:
            src = self.data[self.index[t0 + t, z0 + z, m0 + mu]]
            factor = _pochhammer(t0 + 1, t) * _pochhammer(z0 + 1, z) * _pochhammer(m0 + 1, mu)
            data.append([src[0] * factor, src[1] * factor, src[2] * factor, src[3] * factor])
        return MatrixJet(data, pt, pm, indices, index)

    def __add__(self, other):
        b = other if isinstance(other, MatrixJet) else MatrixJet.constant(other, self.physical, self.momentum)
        pt, pm = min(self.physical, b.physical), min(self.momentum, b.momentum)
        left, right = self.trim(pt, pm), b.trim(pt, pm)
        data = [[x[0] + y[0], x[1] + y[1], x[2] + y[2], x[3] + y[3]]
                for x, y in zip(left.data, right.data)]
        return MatrixJet(data, pt, pm, left.indices, left.index)

    __radd__ = __add__

    def __neg__(self):
        return MatrixJet([[-x[0], -x[1], -x[2], -x[3]] for x in self.data],
                         self.physical, self.momentum, self.indices, self.index)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if isinstance(other, ScalarJet) or isinstance(other, (int, float, arb, acb)):
            return self._scale(other)
        if not isinstance(other, MatrixJet):
            raise TypeError('matrix jet product requires a matrix, scalar jet or scalar')
        pt, pm = min(self.physical, other.physical), min(self.momentum, other.momentum)
        indices, index = _table(pt, pm)
        out = [_zero_mat() for _ in indices]
        left = [(key, mat) for key, mat in zip(self.indices, self.data) if not _is_zero_mat(mat)]
        right = [(key, mat) for key, mat in zip(other.indices, other.data) if not _is_zero_mat(mat)]
        for (t1, z1, m1), a in left:
            for (t2, z2, m2), b in right:
                t, z, m = t1 + t2, z1 + z2, m1 + m2
                if t + z <= pt and m <= pm:
                    prod = _matmul(a, b)
                    slot = out[index[t, z, m]]
                    slot[0] += prod[0]
                    slot[1] += prod[1]
                    slot[2] += prod[2]
                    slot[3] += prod[3]
        return MatrixJet(out, pt, pm, indices, index)

    def _scale(self, other):
        if isinstance(other, ScalarJet):
            pt, pm = min(self.physical, other.physical), min(self.momentum, other.momentum)
            indices, index = _table(pt, pm)
            out = [_zero_mat() for _ in indices]
            for i, (t1, z1, m1) in enumerate(other.indices):
                s = other.data[i]
                if s.is_zero():
                    continue
                for j, (t2, z2, m2) in enumerate(self.indices):
                    t, z, m = t1 + t2, z1 + z2, m1 + m2
                    if t + z <= pt and m <= pm:
                        a = self.data[j]
                        if _is_zero_mat(a):
                            continue
                        slot = out[index[t, z, m]]
                        slot[0] += s * a[0]
                        slot[1] += s * a[1]
                        slot[2] += s * a[2]
                        slot[3] += s * a[3]
            return MatrixJet(out, pt, pm, indices, index)
        scale = other if isinstance(other, (arb, acb)) else acb(other)
        return MatrixJet([[scale * x[0], scale * x[1], scale * x[2], scale * x[3]] for x in self.data],
                         self.physical, self.momentum, self.indices, self.index)

    def __rmul__(self, other):
        if isinstance(other, MatrixJet):
            return other * self
        return self._scale(other)


def _sigma3(physical, momentum):
    return MatrixJet.constant(_mat(_S3), physical, momentum)


def _identity(physical, momentum):
    return MatrixJet.constant(_mat(_I), physical, momentum)


def _mu_coordinate(mu0, physical, momentum):
    return ScalarJet.variable(mu0, 2, physical, momentum)


def momentum_D(B, p, h, sign, mu0):
    """Regular core of canonical k derivatives: d_k^h (mu^p B)=mu^{p+h} D[p,h] B.

    The chain rule is ``d_k=-e mu^2 d_mu`` with ``e=sign(k)``.
    """
    p = _as_int(p, 'D weight p', 0, 32)
    h = _as_int(h, 'k-derivative order', 0, B.momentum)
    if sign not in (1, -1) or isinstance(sign, bool):
        raise ValueError('sign e must be +1 or -1')
    current = B
    for step in range(h):
        dmu = current.derivative(mu=1)
        mu_jet = _mu_coordinate(mu0, current.physical, current.momentum)
        current = acb(-sign) * ((p + step) * current + mu_jet * dmu)
    return current


def scaled_star(Bi, Bj, i, j, order, sign, mu0):
    """Gtilde star of two regular symbols, with P_inf already killed.

    For every factor, including ``i=0`` or ``j=0``, the weight is ``p=index+1``
    because any positive ``T/z/k`` derivative annihilates ``P_inf``.
    """
    order = _as_int(order, 'star order', 0, 16)
    if order == 0:
        return Bi * Bj
    result = MatrixJet.zeros(min(Bi.physical, Bj.physical), min(Bi.momentum, Bj.momentum))
    half = acb(0, 1) / 2
    pref = (half ** order) / factorial(order)
    for h in range(order + 1):
        left = momentum_D(Bi, i + 1, h, sign, mu0).derivative(z=order - h)
        right = momentum_D(Bj, j + 1, order - h, sign, mu0).derivative(z=h)
        result += ((-1) ** h * comb(order, h)) * pref * (left * right)
    return result


def _from_geometry(geometry, field, physical, momentum):
    if geometry.order < physical:
        raise ValueError('geometry order does not supply the requested physical jets')
    indices, index = _table(physical, momentum)
    data = []
    for t, z, m in indices:
        data.append(arb(0) if m else geometry.coefficients[field][(t, z)])
    jet = ScalarJet(data, physical, momentum, indices, index)
    if field == 'a':
        for t, z, m in indices:
            if z and not jet.data[index[t, z, m]].is_zero():
                raise ValueError('axial scale a(T) must be z-independent; kinetic star orders >=2 were dropped')
    return jet


def _require_pure_ks(geometry):
    n00 = geometry.coefficients['N'][(0, 0)]
    if n00 != 1:
        raise ValueError('pure KS history requires N=1, beta=0')
    for pair, value in geometry.coefficients['N'].items():
        if pair != (0, 0) and not value.is_zero():
            raise ValueError('pure KS history requires N=1, beta=0')
    for value in geometry.coefficients['beta'].values():
        if not value.is_zero():
            raise ValueError('pure KS history requires N=1, beta=0')


def _geometry_fields(geometry, axial, radius, physical, momentum):
    if geometry is not None:
        if axial is not None or radius is not None:
            raise ValueError('pass either KSGeometryJets or constant axial/radius, not both')
        if not isinstance(geometry, KSGeometryJets):
            raise TypeError('KSGeometryJets required; arbitrary metric callbacks are not accepted')
        _require_pure_ks(geometry)
        a = _from_geometry(geometry, 'a', physical, momentum)
        r = _from_geometry(geometry, 'r', physical, momentum)
        inv_a = _from_geometry(geometry, 'inv_a', physical, momentum)
        inv_a2 = _from_geometry(geometry, 'inv_a2', physical, momentum)
        identity = geometry.history_identity
    else:
        if axial is None or radius is None:
            raise ValueError('constant axial and radius balls are required when no KS geometry is supplied')
        a = ScalarJet.constant(_real_ball(axial, 'axial'), physical, momentum)
        r = ScalarJet.constant(_real_ball(radius, 'radius'), physical, momentum)
        inv_a = a.inv()
        inv_a2 = inv_a * inv_a
        identity = None
    if not a.value > 0:
        raise ValueError('positive enclosing interval required for axial scale; no denominator floor')
    if not r.value > 0:
        raise ValueError('positive enclosing interval required for radius; no denominator floor')
    return a, r, inv_a, inv_a2, identity


def _regular_B0(b, b2, mu_jet, sign, physical, momentum):
    s2 = ScalarJet.constant(arb(1), physical, momentum) + (mu_jet * mu_jet) * b2
    s = s2.sqrt()
    inv_s = s.inv()
    inv_s1 = (s + arb(1)).inv()
    first = acb(-1) / 2 * (b * inv_s)
    scalar = (arb(sign) / 2) * mu_jet * b2 * inv_s * inv_s1
    return first + MatrixJet.pauli(scalar, _mat(_S3), physical, momentum)


def _Ftilde(B, C, V, n, sign, mu0, mu_jet):
    prev = B[n - 1]
    half = acb(0, 1) / 2
    F = acb(0, 1) * prev.derivative(t=1) + half * (C * prev.derivative(z=1) + prev.derivative(z=1) * C)
    for j in range(n):
        ell = n - j
        dV = V.derivative(z=ell)
        Dj = momentum_D(B[j], j + 1, ell, sign, mu0)
        pref = ((acb(0, 1) / 2) ** ell) / factorial(ell)
        F -= pref * (mu_jet * (dV * Dj - ((-1) ** ell) * (Dj * dV)))
    return F


@dataclass(frozen=True)
class ScaledProjectorJets:
    physical: int
    momentum: int
    formal_order: int
    bits: int
    sign: int
    mu: object
    mass: object
    angular: object
    B: tuple
    Pinf: object
    P0: object
    Gtilde: tuple
    Ftilde: tuple
    retained: tuple
    trace_B: tuple
    trace_P: tuple
    history_identity: object
    auxiliary: bool = True
    state_law: str = 'C_g = U_g C_up U_g^dagger = P_g + D_g exactly'
    recurrence: str = 'Panati-Spohn-Teufel Secs. 2 and 4.4 via nsc_spatial_reference_symbol'
    clock: str = 'T; G_rho = i H_T / a is handled elsewhere'

    def reconstruct_P_value(self, j):
        """Point value of the auxiliary symbol P_j at the expansion mu."""
        j = _as_int(j, 'projector order', 0, self.formal_order)
        with ctx.workprec(self.bits):
            if j == 0:
                return self.P0.value
            scale = self.mu ** (j + 1)
            b = self.B[j].value
            return [scale * b[0], scale * b[1], scale * b[2], scale * b[3]]


def projector_derivative(jets, j, t=0, z=0, k=0):
    """Actual d_T^t d_z^z d_k^k P_j at the expansion point, as a 2x2 acb 4-vector.

    Any positive derivative of P_0 uses ``mu B_0`` because it kills P_inf.
    For ``j>=1``, ``P_j=mu^{j+1} B_j``. This does not assert a UV seminorm.
    """
    j = _as_int(j, 'projector order', 0, jets.formal_order)
    t = _as_int(t, 'time derivative', 0, jets.physical)
    z = _as_int(z, 'axial derivative', 0, jets.physical)
    k = _as_int(k, 'k derivative', 0, jets.momentum)
    with ctx.workprec(jets.bits):
        if j == 0 and t == 0 and z == 0 and k == 0:
            return list(jets.P0.value)
        p = 1 if j == 0 else j + 1
        core = momentum_D(jets.B[j], p, k, jets.sign, jets.mu)
        if t or z:
            core = core.derivative(t=t, z=z)
        scale = jets.mu ** (p + k)
        v = core.value
        return [scale * v[0], scale * v[1], scale * v[2], scale * v[3]]


def symbol_metric_from_geometry(geometry, *, physical=4):
    """Embed the same KS (T,z) jets into the hash-pinned order-four SymbolJet.

    Coefficients are midpoints of the directed geometry balls. This is an
    application of the existing owner, not a generator and not a replacement
    of its global tables.
    """
    import numpy as np
    from .nsc_spatial_reference_symbol import INDEX, INDICES, I, SymbolJet

    if not isinstance(geometry, KSGeometryJets):
        raise TypeError('KSGeometryJets required; arbitrary metric callbacks are not accepted')
    _require_pure_ks(geometry)
    physical = _as_int(physical, 'SymbolJet physical order', 0, 4)

    def field(name):
        data = np.zeros((len(INDICES), 1, 2, 2), complex)
        for t, z, k in INDICES:
            if k or t + z > physical:
                continue
            data[INDEX[t, z, 0], 0] = float(geometry.coefficients[name][(t, z)]) * I
        return SymbolJet(data, 4, 4)

    return (SymbolJet.constant(np.array([1.0])), SymbolJet.constant(np.array([0.0])),
            field('a'), field('r'))


def scaled_reference_projector(*, mu, sign, mass, angular, geometry=None, axial=None, radius=None,
                               physical=DEFAULT_PHYSICAL, momentum=DEFAULT_MOMENTUM,
                               formal_order=FORMAL_ORDER, bits=160):
    """Regular B_0..B_4 on pure KS history, default orders enough for R_z^4.

    ``B_j`` retains physical ``physical-j`` and mu ``momentum-j``. No trace-one
    or CAR projection is applied. Kinetic Weyl orders ``>=2`` vanish because
    ``a=a(T)``.
    """
    physical = _as_int(physical, 'physical order', 0, MAX_PHYSICAL)
    momentum = _as_int(momentum, 'mu momentum order', 0, MAX_MOMENTUM)
    formal_order = _as_int(formal_order, 'formal projector order', 1, FORMAL_ORDER)
    bits = _as_int(bits, 'precision bits', 80, 1024)
    if physical < formal_order or momentum < formal_order:
        raise ValueError('physical and mu orders must cover the requested formal projector order')
    if sign not in (1, -1) or isinstance(sign, bool):
        raise ValueError('sign e must be +1 or -1')
    if geometry is not None:
        if not isinstance(geometry, KSGeometryJets):
            raise TypeError('KSGeometryJets required; arbitrary metric callbacks are not accepted')
        bits = max(bits, int(geometry.bits))
    with ctx.workprec(bits):
        mu0 = _nonneg_mu(mu)
        mass = _real_ball(mass, 'mass')
        angular = _real_ball(angular, 'angular')
        a, r, inv_a, inv_a2, identity = _geometry_fields(geometry, axial, radius, physical, momentum)
        inv_r = r.inv()
        mu_jet = _mu_coordinate(mu0, physical, momentum)
        C = MatrixJet.pauli(inv_a, _mat(_S3), physical, momentum)
        V = MatrixJet.pauli(ScalarJet.constant(-mass, physical, momentum), _mat(_S1), physical, momentum)
        V = V + acb(0, 1) * MatrixJet.pauli(angular * inv_r, _mat(_S2), physical, momentum)
        b = a * V
        b2 = (a * a) * (mass * mass + (angular * angular) * (inv_r * inv_r))
        B0 = _regular_B0(b, b2, mu_jet, sign, physical, momentum)
        Pinf = MatrixJet.constant(_mat(((1 - sign) // 2, 0, 0, (1 + sign) // 2)), physical, momentum)
        P0 = Pinf + mu_jet * B0
        Hbar = acb(sign) * C + mu_jet * V
        gapbar2 = inv_a2 + (mu_jet * mu_jet) * (mass * mass + (angular * angular) * (inv_r * inv_r))
        inv_4gap = (arb(4) * gapbar2).inv()
        Ijet = _identity(physical, momentum)
        Q = Ijet - P0
        B = [B0]
        Gtilde = []
        Ftilde = []
        for n in range(1, formal_order + 1):
            G = MatrixJet.zeros(physical, momentum)
            for i in range(n):
                for j in range(n):
                    ell = n - i - j
                    if ell < 0:
                        continue
                    G = G + scaled_star(B[i], B[j], i, j, ell, sign, mu0)
            F = _Ftilde(B, C, V, n, sign, mu0, mu_jet)
            off = inv_4gap * (Hbar * F - F * Hbar)
            diag = mu_jet * ((acb(-1) * (P0 * G * P0)) + (Q * G * Q))
            Bn = (off + diag).trim(physical - n, momentum - n)
            B.append(Bn)
            Gtilde.append(G.trim(physical - n, momentum - n))
            Ftilde.append(F.trim(physical - n, momentum - n))
        traces_B = tuple(Bj.trace_value() for Bj in B)
        traces_P = [P0.trace_value()]
        for j in range(1, len(B)):
            scale = mu0 ** (j + 1)
            traces_P.append(scale * traces_B[j])
        retained = tuple((physical - j, momentum - j) for j in range(len(B)))
        return ScaledProjectorJets(
            physical=physical, momentum=momentum, formal_order=formal_order, bits=bits,
            sign=sign, mu=mu0, mass=mass, angular=angular, B=tuple(B), Pinf=Pinf, P0=P0,
            Gtilde=tuple(Gtilde), Ftilde=tuple(Ftilde), retained=retained,
            trace_B=tuple(traces_B), trace_P=tuple(traces_P), history_identity=identity)
