"""Uniform-in-k operator-norm jets of the auxiliary fourth-order projector.

Root already owns small-transfer coefficients through
``scaled_reference_projector``. These majorants bound every mixed Taylor
coefficient of auxiliary ``P_0..P_4`` at every real ``k``, which is the
input needed when a large axial transfer ``k ± omega/2`` crosses ``0``.

The construction is a finite-order formal majorant, not a claim that the
Taylor series of ``P`` converges on a disk. Third-axis coefficients are
ordinary ``k`` derivatives, not ``mu``. The symbols remain auxiliary:
they are not ``C_up`` or ``C_Sigma``, they do not replace the physical
state, and they do not add a ``Gamma`` term.
"""
from dataclasses import dataclass
from math import comb, factorial
from types import MappingProxyType

from flint import acb, arb, ctx

from .nsc_ks_spacetime_geometry_jets import KSGeometryJets
from .nsc_scaled_reference_projector import (
    DEFAULT_MOMENTUM,
    DEFAULT_PHYSICAL,
    FORMAL_ORDER,
    MAX_MOMENTUM,
    MAX_PHYSICAL,
    ScalarJet,
    _as_int,
    _geometry_fields,
    _real_ball,
    positive_inv,
)

# 2x2 conversions, not used to inflate the operator-norm majorant.
# ||M||_F <= sqrt(2) ||M||_op and ||M||_1 <= 2 ||M||_op on C^2.
FROBENIUS_OVER_OPERATOR = arb(2).sqrt()
NUCLEAR_OVER_OPERATOR = arb(2)

# Analytic factors after verifying the proposed all-k route.
# The k-order-1 gap increment uses AM-GM |k|/g^2 <= 1/(2 c v). Times the
# 2 from d_k(k^2) this is |c^2_alpha|/(c_min g_min), not twice that.
ANALYTIC_FACTORS = MappingProxyType({
    'X_physical_k0': '|(c^2)_alpha|/c_min^2 + |(v^2)_alpha|/g_min^2',
    'X_korder1': '|(c^2)_alpha|/(c_min g_min)',
    'X_korder1_reason': '|k|/g^2 <= 1/(2 c v) by AM-GM; times 2 from d_k(k^2)',
    'X_korder2': '|(c^2)_alpha|/g_min^2 (Taylor coefficient; 2! already removed)',
    'star': '1/(2^l l!) * C(l,h) on mixed (z,k) splits',
    'F_time': '||d_T P_{n-1}||',
    'F_kinetic': '||C|| ||d_z P_{n-1}|| (Weyl order 1; kinetic orders >=2 vanish)',
    'F_potential': '2/(2^l l!) ||d_z^l V|| ||d_k^l P_j||, l=n-j',
    'offdiag': '1/2 * A * F * (1-X)^{-1} / g_min',
    'commutator_matrix_factor': '2',
    'diag': '2 P0 G P0',
    'P0': '1/2 A (1-X)^{-1/2} + 1/2',
    'frobenius_over_operator': 'sqrt(2) (conversion only; not applied)',
    'nuclear_over_operator': '2 (conversion only; not applied)',
})


def _abs_upper(value):
    """Directed upper of |value| for a real or complex ball."""
    if isinstance(value, acb):
        result = value.abs_upper()
    else:
        if not isinstance(value, arb):
            value = arb(value)
        result = abs(value).upper()
    if not result.is_finite():
        raise ArithmeticError('finite coefficient upper required')
    if not result >= 0:
        raise ArithmeticError('nonnegative coefficient upper required')
    return result


def _as_majorant(jet):
    """Replace each coefficient by a nonnegative directed upper."""
    data = [arb(_abs_upper(coeff)) for coeff in jet.data]
    return ScalarJet(data, jet.physical, jet.momentum, jet.indices, jet.index)


def _reciprocal_upper_from_lower(lower_bound, name='positive lower bound'):
    """Directed upper of 1/L for a positive lower bound L.

    Inverts the enclosing ball and takes its upper, so a thin enclosure of a
    lower endpoint cannot be rounded upward before inversion.
    """
    if not isinstance(lower_bound, arb):
        lower_bound = arb(lower_bound)
    if not lower_bound > 0:
        raise ValueError(f'positive {name} required; no denominator floor')
    return positive_inv(lower_bound, name).upper()


def _abs_lower(value, name):
    """Directed lower of |x| over a real ball. Zero if the ball contains 0.

    flint ``abs`` on a symmetric ball about 0 still contains negatives, so
    ``abs(x).lower()`` is not inf(|x|).
    """
    value = value if isinstance(value, arb) else _real_ball(value, name)
    if value.contains(0):
        return arb(0)
    return abs(value).lower()


def rest_gap_lower(mass, angular, radius):
    """Directed lower of sqrt(m_min^2 + (absell_min/r_max)^2).

    Inverts the radius ball, then takes lowers. Never inverts a rounded
    upper of r. Zero rest gap is rejected; there is no floor.
    """
    mass = mass if isinstance(mass, arb) else _real_ball(mass, 'mass')
    angular = angular if isinstance(angular, arb) else _real_ball(angular, 'angular')
    radius = radius if isinstance(radius, arb) else _real_ball(radius, 'radius')
    if not radius > 0:
        raise ValueError('positive enclosing interval required for radius; no denominator floor')
    m_lo = _abs_lower(mass, 'mass')
    ell_lo = _abs_lower(angular, 'angular')
    inv_r_lo = positive_inv(radius, 'radius').lower()
    rest = m_lo * m_lo + (ell_lo * inv_r_lo) * (ell_lo * inv_r_lo)
    rest_lo = rest.lower()
    if not rest_lo > 0:
        raise ValueError(
            'nonzero rest gap required; g_min=sqrt(m_min^2+(absell_min/r_max)^2) '
            'must be strictly positive. Zero-gap ell=0 radius response is handled '
            'elsewhere; no denominator floor'
        )
    return rest_lo.sqrt().lower()


def axial_c_min(axial):
    """Directed lower of c=1/a from the enclosing axial ball."""
    axial = axial if isinstance(axial, arb) else _real_ball(axial, 'axial scale')
    if not axial > 0:
        raise ValueError('positive enclosing interval required for axial scale; no denominator floor')
    return positive_inv(axial, 'axial scale').lower()


def _kderiv(jet, t=0, z=0, k=0):
    """Ordinary k derivative on a ScalarJet whose third axis is not named k."""
    return jet.derivative(t=t, z=z, mu=k)


def star_norm(left, right, order):
    """Operator-norm majorant of the spatial Weyl product of order ``l``.

    Absolute coefficient is ``1/(2^l l!) * C(l,h)`` on each mixed
    ``(d_z^{l-h} d_k^h, d_z^h d_k^{l-h})`` split, matching
    ``nsc_spatial_reference_symbol.star_order``.
    """
    order = _as_int(order, 'star order', 0, 16)
    if order == 0:
        return left * right
    physical = min(left.physical, right.physical)
    momentum = min(left.momentum, right.momentum)
    result = ScalarJet.zeros(physical, momentum)
    pref = arb(1) / (arb(2) ** order * factorial(order))
    for h in range(order + 1):
        term = _kderiv(left, z=order - h, k=h) * _kderiv(right, z=h, k=order - h)
        result = result + (pref * comb(order, h)) * term
    return result


def _gap_increment_majorant(c2, v2, c_min, g_min, physical, momentum):
    """Nonnegative X with X000=0 majorizing |g^2/g0^2 - 1| uniformly in k."""
    X = ScalarJet.zeros(physical, momentum)
    inv_cmin2 = _reciprocal_upper_from_lower(c_min * c_min, 'c_min^2')
    inv_gmin2 = _reciprocal_upper_from_lower(g_min * g_min, 'g_min^2')
    inv_cmin_gmin = _reciprocal_upper_from_lower(c_min * g_min, 'c_min g_min')
    for t, z, k in X.indices:
        if t == 0 and z == 0 and k == 0:
            continue
        if k > 2:
            continue
        c2_abs = arb(_abs_upper(c2.coefficient(t, z, 0)))
        if k == 0:
            v2_abs = arb(_abs_upper(v2.coefficient(t, z, 0)))
            X.data[X.index[t, z, k]] = c2_abs * inv_cmin2 + v2_abs * inv_gmin2
        elif k == 1:
            X.data[X.index[t, z, k]] = c2_abs * inv_cmin_gmin
        else:
            X.data[X.index[t, z, k]] = c2_abs * inv_gmin2
    return X


def _H_over_g0_majorant(c, inv_r, mass, angular, c_min, g_min, physical, momentum):
    """Operator-norm majorant A of H/g0. A000=1 exactly."""
    A = ScalarJet.zeros(physical, momentum)
    A.data[0] = arb(1)
    inv_cmin = _reciprocal_upper_from_lower(c_min, 'c_min')
    inv_gmin = _reciprocal_upper_from_lower(g_min, 'g_min')
    ell_abs = arb(_abs_upper(angular))
    for t, z, k in A.indices:
        if t == 0 and z == 0 and k == 0:
            continue
        c_abs = arb(_abs_upper(c.coefficient(t, z, 0)))
        if k == 0:
            invr_abs = arb(_abs_upper(inv_r.coefficient(t, z, 0)))
            v_op = ell_abs * invr_abs
            A.data[A.index[t, z, k]] = c_abs * inv_cmin + v_op * inv_gmin
        elif k == 1:
            A.data[A.index[t, z, k]] = c_abs * inv_gmin
    return A


def _V_op_majorant(inv_r, mass, angular, physical, momentum):
    """Pauli-orthogonal op-norm majorant of V=-m Sigma1 + (ell/r) Sigma2."""
    V = ScalarJet.zeros(physical, momentum)
    m_abs = arb(_abs_upper(mass))
    ell_abs = arb(_abs_upper(angular))
    for t, z, k in V.indices:
        if k:
            continue
        invr_abs = arb(_abs_upper(inv_r.coefficient(t, z, 0)))
        m_term = m_abs if (t, z) == (0, 0) else arb(0)
        V.data[V.index[t, z, k]] = (m_term * m_term + (ell_abs * invr_abs) * (ell_abs * invr_abs)).sqrt()
    return V


def _F_norm(P, Cnorm, Vnorm, n, physical, momentum):
    """Majorant of the order-n transport defect F_n."""
    prev = P[n - 1]
    result = prev.derivative(t=1)
    result = result + Cnorm * prev.derivative(z=1)
    two = arb(2)
    for j in range(n):
        ell = n - j
        pref = two / (arb(2) ** ell * factorial(ell))
        term = Vnorm.derivative(z=ell) * _kderiv(P[j], k=ell)
        result = result + pref * term
    return result


@dataclass(frozen=True)
class GlobalProjectorNormJets:
    """Operator-norm majorants of auxiliary P_j Taylor coefficients in (T,z,k)."""

    physical: int
    momentum: int
    formal_order: int
    bits: int
    mass: object
    angular: object
    P: tuple
    X: object
    A: object
    G: tuple
    F: tuple
    c_min: object
    g_min: object
    retained: tuple
    history_identity: object
    analytic_factors: object = ANALYTIC_FACTORS
    frobenius_over_operator: object = FROBENIUS_OVER_OPERATOR
    nuclear_over_operator: object = NUCLEAR_OVER_OPERATOR
    auxiliary: bool = True
    state_law: str = 'C_g = U_g C_up U_g^dagger = P_g + D_g exactly'
    recurrence: str = 'Panati-Spohn-Teufel Secs. 2 and 4.4 via nsc_spatial_reference_symbol'
    k_axis: str = 'ordinary k derivatives, not mu'
    clock: str = 'T; G_rho = i H_T / a is handled elsewhere'
    uniformity: str = 'normalized gap inequalities on the geometry box, all real k'

    def taylor_bound(self, j, t=0, z=0, k=0):
        """Upper bound on || d_T^t d_z^z d_k^k P_j / (t! z! k!) ||_op."""
        j = _as_int(j, 'projector order', 0, self.formal_order)
        with ctx.workprec(self.bits):
            t = _as_int(t, 'time order', 0, self.P[j].physical)
            z = _as_int(z, 'axial order', 0, self.P[j].physical)
            k = _as_int(k, 'momentum order', 0, self.P[j].momentum)
            return self.P[j].coefficient(t, z, k)

    def derivative_bound(self, j, t=0, z=0, k=0):
        """Upper bound on || d_T^t d_z^z d_k^k P_j ||_op = Taylor * t! z! k!."""
        j = _as_int(j, 'projector order', 0, self.formal_order)
        with ctx.workprec(self.bits):
            t = _as_int(t, 'time derivative', 0, self.P[j].physical)
            z = _as_int(z, 'axial derivative', 0, self.P[j].physical)
            k = _as_int(k, 'momentum derivative', 0, self.P[j].momentum)
            return _kderiv(self.P[j], t=t, z=z, k=k).value


def global_projector_bounds(*, mass, angular, geometry=None, axial=None, radius=None,
                            physical=DEFAULT_PHYSICAL, momentum=DEFAULT_MOMENTUM,
                            formal_order=FORMAL_ORDER, bits=160):
    """Bound ||Taylor coeffs of auxiliary P_0..P_4||_op uniformly in real k.

    Default physical 9 and k 5 are enough for a fourth axial derivative of
    ``P_4``. Each order is trimmed to ``(physical-j, momentum-j)``. Gap-zero
    geometry is rejected. The third axis is ordinary ``k``, not ``mu``.
    """
    physical = _as_int(physical, 'physical order', 0, MAX_PHYSICAL)
    momentum = _as_int(momentum, 'k momentum order', 0, MAX_MOMENTUM)
    formal_order = _as_int(formal_order, 'formal projector order', 0, FORMAL_ORDER)
    bits = _as_int(bits, 'precision bits', 80, 1024)
    if physical < formal_order or momentum < formal_order:
        raise ValueError('physical and k orders must cover the requested formal projector order')
    if geometry is not None:
        if not isinstance(geometry, KSGeometryJets):
            raise TypeError('KSGeometryJets required; arbitrary metric callbacks are not accepted')
        bits = max(bits, int(geometry.bits))
    with ctx.workprec(bits):
        mass = _real_ball(mass, 'mass')
        angular = _real_ball(angular, 'angular')
        a, r, inv_a, inv_a2, identity = _geometry_fields(
            geometry, axial, radius, physical, momentum)
        c_min = axial_c_min(a.value)
        g_min = rest_gap_lower(mass, angular, r.value)
        inv_r = r.inv()
        v2 = ScalarJet.constant(mass * mass, physical, momentum) + (
            (angular * angular) * (inv_r * inv_r))
        X = _gap_increment_majorant(inv_a2, v2, c_min, g_min, physical, momentum)
        A = _H_over_g0_majorant(inv_a, inv_r, mass, angular, c_min, g_min, physical, momentum)
        one = ScalarJet.constant(arb(1), physical, momentum)
        omx = one - X
        g0_over_g = _as_majorant(omx.sqrt().inv())
        g0sq_over_gsq = _as_majorant(omx.inv())
        half = arb(1) / 2
        P0 = _as_majorant(half * A * g0_over_g + half)
        Cnorm = _as_majorant(inv_a)
        Vnorm = _V_op_majorant(inv_r, mass, angular, physical, momentum)
        inv_gmin = _reciprocal_upper_from_lower(g_min, 'g_min')
        P = [P0]
        Gjets = []
        Fjets = []
        for n in range(1, formal_order + 1):
            G = ScalarJet.zeros(physical, momentum)
            for i in range(n):
                for j in range(n):
                    ell = n - i - j
                    if ell < 0:
                        continue
                    G = G + star_norm(P[i], P[j], ell)
            F = _F_norm(P, Cnorm, Vnorm, n, physical, momentum)
            off = half * A * F * g0sq_over_gsq * inv_gmin
            diag = (2 * P0) * G * P0
            Pn = _as_majorant((diag + off).trim(physical - n, momentum - n))
            P.append(Pn)
            Gjets.append(_as_majorant(G.trim(physical - n, momentum - n)))
            Fjets.append(_as_majorant(F.trim(physical - n, momentum - n)))
        retained = tuple((physical - j, momentum - j) for j in range(len(P)))
        return GlobalProjectorNormJets(
            physical=physical, momentum=momentum, formal_order=formal_order, bits=bits,
            mass=mass, angular=angular, P=tuple(P), X=X, A=A,
            G=tuple(Gjets), F=tuple(Fjets), c_min=c_min, g_min=g_min,
            retained=retained, history_identity=identity)
