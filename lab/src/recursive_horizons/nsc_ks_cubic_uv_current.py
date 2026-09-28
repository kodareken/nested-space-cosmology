"""Generated massless e^{-3} N,beta coefficient from actual L0 vertices.

Stage-2 UV construction on the switched state law, live profile 0b0e4ced and
the smallest original angular pair (group 1, massless, |ell|=sqrt(5)). Common
scalar upstream constants drop out of the Sigma difference by the accepted
surface invariance; they are not a substitute for the generated A3/n4
transport. This module constructs that generated cubic coefficient.

Constructed here
----------------
Massless reductions b, l, B, n2, R, minor C, T_s h3 and T_s n4, with every
linked term retained. The contracted e^{-3} N and beta densities from the
actual Pauli vertices

    N = -mu/(2 pi) (V density + S3 current / a),
    beta = +mu/(2 pi) (I current),

checked against source_column_matter. A bounded current-history quadrature
evaluates the history-minus-reference coefficient on I for a consistent
representative (common k, consistent R, C, h3_z, n4).

Not constructed here
--------------------
A remainder-validated C4, C_M, L0 A4 H2 integrals, physical source-column
error, massive channels, or a closed local gate. Diagnostic quadrature is
a measured coefficient, not a proof. Common k is a representative of the
invariance class; n2_up is not set to zero as a physical statement.
"""
from collections.abc import Mapping
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
import json

import numpy as np
import sympy as sp
from scipy.special import expit

from .nsc_dirac_source_phase import (
    RHO_SIGMA,
    SOURCE_SIGNS,
    _characteristic_distance,
    _normal_coordinate,
)
from .nsc_evolved_incoming_constraints import source_column_matter
from .nsc_ks_current_uv_transport import (
    ARCHIVE_RHO_UP,
    CURRENT_HISTORY_IDENTITY,
    HISTORY_RECORD,
    smallest_original_pair,
)
from .nsc_ks_finite_history_uv_coefficients import (
    _T,
    _apply_L0,
    _as_one_direction,
    _major_index,
    _minor_A1,
    _minor_index,
    _projector,
    _real_residual,
    _sym_pauli,
    _sym_potential,
    _unit,
    confirm_raw_ks_vertices,
    vacuum_upstream_conditions,
)
from .nsc_ks_finite_history_uv_remainder import FIRST_NONCANCELLING_PAIRED_ORDER
from .nsc_ks_profile_identity import profile_identity
from .nsc_ks_uv_minor_a3 import _parent_minor_a2
from .nsc_local_incoming_family import LocalIncomingFamily
from .nsc_pg_ks_metric_pullback import reference_chart


SCHEMA = 'NSC-KS-CUBIC-UV-CURRENT-v1'
TWO_PI = 2.0 * np.pi
DEFAULT_RHO_STEPS = 24
DEFAULT_Z_POINTS = 21
COARSE_RHO_STEPS = 16
COARSE_Z_POINTS = 21
MISSING_REMAINDER = (
    'production L0 A4 H2 integrals and the same-column upstream H2 remainder; '
    'the diagnostic e^{-3} quadrature is unvalidated without that remainder')
NEXT_PRIMITIVE = (
    'a directed remainder bound on the truncated envelope (C_M) and the '
    'physical source-column error; this cubic coefficient is the generated '
    'occupied-vacuum e^{-3} difference on I, not that remainder')


def _root():
    return Path(__file__).resolve().parents[2]


def _as_map(value, name):
    if not isinstance(value, Mapping):
        raise ValueError(name + ' must be a mapping')
    return value


def _require_missing(value, name):
    if value is None:
        return
    if value == 0 or value is False:
        raise ValueError('null ' + name + ' presented as zero')
    raise ValueError('missing ' + name + ' may not be invented')


def _require_massless_angular(mass, angular):
    if isinstance(mass, (bool, np.bool_)) or mass is None:
        raise ValueError('explicit real mass required')
    if isinstance(angular, (bool, np.bool_)) or angular is None:
        raise ValueError('explicit angular eigenvalue required')
    mass = float(mass)
    angular = float(angular)
    if not np.isfinite([mass, angular]).all():
        raise ValueError('finite mass and angular eigenvalue required')
    if mass != 0.0:
        raise ValueError(
            'massive cubic coefficient is not this massless reduction; '
            'the owned evaluator rejects m!=0 instead of dropping mass terms')
    if angular == 0.0:
        raise ValueError('nonzero original angular eigenvalue required')
    return mass, angular


def _zero(expr):
    residual = _real_residual(expr)
    if residual == '0':
        return '0'
    real_derivatives = [
        d for d in expr.atoms(sp.Derivative)
        if d.expr.is_real is True and all(v.is_real is True for v in d.variables)]
    expr = expr.xreplace({d: sp.Dummy(real=True) for d in real_derivatives})
    value = sp.simplify(sp.expand_complex(expr))
    if value != 0:
        raise ArithmeticError('cubic UV identity failed: ' + str(value))
    return '0'


def _hex(value):
    return float(value).hex()


@lru_cache(maxsize=1)
def confirm_action_signs_against_source_column_matter():
    """N = -N_bracket and beta = +J_I, with the raw -I vertex minus once.

    The older remainder shorthand action_beta3 = -mu/(2 pi) J_3 is not this
    identity. Columns already contain sqrt(dE/(2 pi)); this check does not
    insert a second measure.
    """
    vertices = confirm_raw_ks_vertices()
    checks = {}
    for source_sign in SOURCE_SIGNS:
        columns = np.zeros((1, 2, 1), complex)
        columns[0, _major_index(source_sign), 0] = 1
        axial = -2j * source_sign * columns
        empty = np.zeros((0, *columns.shape), complex)
        result = source_column_matter(
            columns, axial, np.eye(1), empty, empty,
            mass=0, angular=0, axial_scale=1, radius=1, multiplicity=1)
        expected = np.array([2.0, -2.0 * source_sign], dtype=float)
        if not np.array_equal(result['action_gradient'][0], expected):
            raise ArithmeticError({
                'source_sign': source_sign,
                'action_gradient': result['action_gradient'][0].tolist(),
                'expected': expected.tolist(),
            })
        # Occupied A0 with A_z = -2 i s A0: J_I = Im(A^dag A_z) = -2 s.
        # beta = +J_I = -2 s. Old v3 shorthand -J_I would be +2 s.
        checks['s%+d/N' % source_sign] = float(result['action_gradient'][0, 0])
        checks['s%+d/beta' % source_sign] = float(result['action_gradient'][0, 1])
        checks['s%+d/beta_equals_plus_J_I' % source_sign] = True
        checks['s%+d/old_v3_shorthand_beta_minus_J_I' % source_sign] = False
    return MappingProxyType({
        'raw_vertices': dict(vertices),
        'N': '-mu/(2 pi) (V density + S3 current / a)',
        'beta': '+mu/(2 pi) I current',
        'old_v3_beta_shorthand_is_authority': False,
        'second_two_pi_inserted': False,
        'source_column_matter_residuals': checks,
    })


@lru_cache(maxsize=1)
def massless_cubic_identities():
    """Closed massless A3/n4 transport and the actual e^{-3} Pauli contraction."""
    rho, z = sp.symbols('rho z', real=True)
    axial = sp.Function('a', positive=True)(rho)
    radius = sp.Function('r', positive=True)(rho, z)
    phase = sp.Function('q', real=True)(rho, z)
    height = sp.Function('h', real=True)(rho, z)
    real_major = sp.Function('R', real=True)(rho, z)
    height3 = sp.Function('h3', real=True)(rho, z)
    real3 = sp.Function('R3', real=True)(rho, z)
    n4 = sp.Function('n4', real=True)(rho, z)
    moment = sp.symbols('k', real=True)
    angular = sp.symbols('ell', real=True)
    mass = sp.Integer(0)
    _, _, s3 = _sym_pauli()
    potential = _sym_potential(mass, angular, radius)
    identities = {}
    formulas = {}
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        major = _major_index(source_sign)
        minor = _minor_index(source_sign)
        _, vacant = _projector(source_sign, s3)
        b = axial * angular / (2 * radius)
        l = _minor_A1(axial, mass, angular, radius, source_sign)
        identities['%s/l_is_minus_i_s_b' % tag] = _zero(l + sp.I * source_sign * b)
        a_rho = axial.diff(rho)
        r_rho = radius.diff(rho)
        r_z = radius.diff(z)
        b_rho = sp.diff(b, rho)
        b_z = sp.diff(b, z)
        parent_B = _parent_minor_a2(
            axial, radius, phase, mass, angular, source_sign, a_rho, r_rho, r_z)
        reduction_B = (
            source_sign * b * phase
            - source_sign * axial**2 * b_rho / 2
            - b_z / 2)
        identities['%s/B_reduction' % tag] = _zero(sp.simplify(parent_B - reduction_B))
        identities['%s/B_is_real' % tag] = _zero(sp.im(sp.expand_complex(parent_B)))
        a1 = sp.zeros(2, 1)
        a1[major] = sp.I * phase
        a1[minor] = l
        square = angular**2 / radius**2
        a1_rho = source_sign / axial**2 * a1[major].diff(z) - sp.I / 2 * square
        l0_a1 = _apply_L0(a1, axial, potential, s3, rho, z).subs(
            a1[major].diff(rho), a1_rho)
        minor_a2 = axial**2 / (2 * sp.I) * (vacant * l0_a1)[minor]
        identities['%s/B_from_L0_A1' % tag] = _zero(minor_a2 - parent_B)
        n2 = 2 * real_major + phase**2 + b**2
        r_from_k = (source_sign * phase.diff(z) + moment - phase**2 - b**2) / 2
        identities['%s/n2_is_s_qz_plus_k' % tag] = _zero(
            n2.subs(real_major, r_from_k) - (source_sign * phase.diff(z) + moment))
        a2 = sp.zeros(2, 1)
        a2[major] = real_major + sp.I * height
        a2[minor] = parent_B
        rhs_a2 = (sp.I / axial) * potential[major, minor] * parent_B
        a2_rho = source_sign / axial**2 * a2[major].diff(z) + rhs_a2
        l0_a2 = _apply_L0(a2, axial, potential, s3, rho, z).subs(
            a2[major].diff(rho), a2_rho)
        minor_a3 = axial**2 / (2 * sp.I) * (vacant * l0_a2)[minor]
        vacant_potential = sp.I * source_sign * angular / radius
        from_t = axial**2 / (2 * sp.I) * _T(parent_B, -source_sign, axial, rho, z)
        algebraic = -(axial / 2) * vacant_potential * (real_major + sp.I * height)
        identities['%s/C_from_recurrence' % tag] = _zero(minor_a3 - from_t - algebraic)
        parent_C = (
            source_sign * b * height
            - sp.I * (axial**2 / 2 * _T(parent_B, -source_sign, axial, rho, z)
                      + source_sign * b * real_major))
        identities['%s/C_massless_parent' % tag] = _zero(minor_a3 - parent_C)
        rhs_a3 = (sp.I / axial) * potential[major, minor] * minor_a3
        tminus_B = _T(parent_B, -source_sign, axial, rho, z)
        parent_h3 = (
            -source_sign * b * tminus_B
            - 2 * b**2 * real_major / axial**2)
        parent_R3 = 2 * b**2 * height / axial**2
        identities['%s/Ts_h3' % tag] = _zero(
            sp.im(sp.expand_complex(rhs_a3)) - parent_h3)
        identities['%s/Ts_R3' % tag] = _zero(
            sp.re(sp.expand_complex(rhs_a3)) - parent_R3)
        re_lc = sp.re(sp.expand_complex(sp.conjugate(l) * minor_a3))
        depth = axial**2 / 2 * tminus_B + source_sign * b * real_major
        identities['%s/Re_conj_l_C_is_s_b_D' % tag] = _zero(re_lc - source_sign * b * depth)
        energy = sp.symbols('e', positive=True)
        n2f = sp.Function('n2', real=True)(rho, z)
        n3f = sp.Function('n3', real=True)(rho, z)
        minor_sq = (
            sp.conjugate(l) * l / energy**2
            + (parent_B**2 + 2 * re_lc) / energy**4)
        density = 1 + n2f / energy**2 + n3f / energy**3 + n4 / energy**4
        current_s = source_sign * (density - 2 * minor_sq)
        continuity = density.diff(rho) - current_s.diff(z) / axial**2
        order4 = sp.expand(continuity * energy**4).subs(energy, 0)
        expected_n4 = _T(n4, source_sign, axial, rho, z) + (
            source_sign / axial**2) * sp.diff(2 * parent_B**2 + 4 * re_lc, z)
        identities['%s/continuity_e_minus4_is_Ts_n4_source' % tag] = _zero(
            order4 - expected_n4)
        a0 = _unit(source_sign)
        a3 = sp.zeros(2, 1)
        a3[major] = real3 + sp.I * height3
        a3[minor] = minor_a3
        a0z = sp.zeros(2, 1)
        interference = (
            a0.H * a3.diff(z) + a1.H * a2.diff(z) + a2.H * a1.diff(z)
            + a3.H * a0z)[0]
        j_i = sp.im(sp.expand_complex(interference)) - source_sign * n4
        parent_ji = (
            height3.diff(z) + real_major * phase.diff(z)
            - phase * real_major.diff(z)
            + source_sign * (b * parent_B.diff(z) - parent_B * b_z)
            - source_sign * n4)
        identities['%s/J_I_3' % tag] = _zero(j_i - parent_ji)
        s3_interf = (
            a0.H * s3 * a3.diff(z) + a1.H * s3 * a2.diff(z)
            + a2.H * s3 * a1.diff(z))[0]
        minor_sq4 = parent_B**2 + 2 * re_lc
        as3a4 = source_sign * (n4 - 2 * minor_sq4)
        j_s3 = sp.im(sp.expand_complex(s3_interf)) - source_sign * as3a4
        v3 = sp.expand_complex(
            (a0.H * potential * a3 + a1.H * potential * a2
             + a2.H * potential * a1 + a3.H * potential * a0)[0])
        v3 = sp.re(v3)
        n_bracket = v3 + j_s3 / axial
        identities['%s/N_bracket_retains_n4' % tag] = (
            'depends' if n4 in n_bracket.atoms(sp.Function) else 'independent')
        identities['%s/J_I_3_coefficient_of_h' % tag] = _zero(
            sp.diff(j_i, height) if height in j_i.atoms(sp.Function) else sp.Integer(0))
        identities['%s/J_I_3_retains_h3_n4_R_B' % tag] = (
            'depends' if (
                height3 in j_i.atoms(sp.Function)
                and n4 in j_i.atoms(sp.Function)
                and real_major in j_i.atoms(sp.Function)) else 'independent')
        formulas[tag] = {
            'b': 'a ell / (2 r)',
            'l': '-i s b',
            'B': 's b q - s a^2 b_rho/2 - b_z/2',
            'R': '(s q_z + k - q^2 - b^2)/2',
            'C': 's b h - i (a^2/2 T_{-s} B + s b R)',
            'T_s_h3': '-s b T_{-s} B - 2 b^2 R / a^2',
            'T_s_R3': '2 b^2 h / a^2',
            'T_s_n4': '-(s/a^2) d_z (2 B^2 + 4 Re(conj(l) C))',
            'J_I_3': (
                'h3_z + R q_z - q R_z + s (b B_z - B b_z) - s n4'),
            'N': '-mu/(2 pi) (V_3 + J_S3_3 / a)',
            'beta': '+mu/(2 pi) J_I_3',
        }
    # Linked real-A2 shift of the scalar series: R->R+c, h3_z->h3_z+c q_z,
    # n4->n4+2 c s q_z. Independently shifting R leaves the spurious c q_z.
    c, qz, rr, rz, q, h3z, n4s, bb, bz, BB, Bz = sp.symbols(
        'c qz R Rz q h3z n4 b bz B Bz', real=True)
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        current = (
            h3z + rr * qz - q * rz
            + source_sign * (bb * Bz - BB * bz) - source_sign * n4s)
        linked = current.xreplace({
            rr: rr + c, h3z: h3z + c * qz,
            n4s: n4s + 2 * c * source_sign * qz})
        identities['%s/linked_scalar_R_cancels_in_J_I_3' % tag] = _zero(
            linked - current)
        unlinked = current.xreplace({rr: rr + c})
        identities['%s/unlinked_R_spurious_c_qz' % tag] = _zero(
            (unlinked - current) - c * qz)
    rg, rr, aa, kk = sp.symbols('r_g r_ref a_common k_common', positive=True)
    identities['Sigma_matched_radius_drops_delta_b'] = _zero(
        (aa * angular / (2 * rg) - aa * angular / (2 * rr)).subs(rg, rr))
    residual_keys = [
        key for key, value in identities.items()
        if value not in ('depends', 'independent')]
    leaves = [identities[key] for key in residual_keys]
    if set(leaves) != {'0'}:
        raise ArithmeticError(identities)
    for key, value in identities.items():
        if key.endswith('J_I_3_retains_h3_n4_R_B') and value != 'depends':
            raise ArithmeticError('J_I_3 must retain generated h3, n4 and R')
        if key.endswith('N_bracket_retains_n4') and value != 'depends':
            raise ArithmeticError('N_bracket_3 must retain n4')
    signs = confirm_action_signs_against_source_column_matter()
    return MappingProxyType({
        'residuals': identities,
        'formulas': formulas,
        'massless_reductions': {
            'b': 'a ell / (2 r)',
            'l': '-i s b',
            'B': 's b q - s a^2 b_rho/2 - b_z/2',
            'n2': 's q_z + k',
            'R': '(s q_z + k - q^2 - b^2)/2',
            'C': 's b h - i (a^2/2 T_{-s} B + s b R)',
            'T_s_h3': '-s b T_{-s} B - 2 b^2 R / a^2',
            'T_s_n4': '-(s/a^2) d_z (2 B^2 + 4 Re(conj(l) C))',
        },
        'J_I_3': 'h3_z + R q_z - q R_z + s (b B_z - B b_z) - s n4',
        'action_N3': '-mu/(2 pi) (V_3 + J_S3_3 / a)',
        'action_beta3': '+mu/(2 pi) J_I_3',
        'old_v3_beta_shorthand': '-mu/(2 pi) J_I_3, not authority',
        'action_signs': dict(signs),
        'dummy_symbol_L0_A_j': False,
        'A3_locally_algebraic_in_metric_jets': False,
        'R_zeroed_while_freezing_h3_n4': False,
        'invariance_uses_n2_up_equals_zero': False,
        'k_is_common_representative': True,
        'all_linked_C_R_h3_n4_retained': True,
        'massless_Im_Ts_major_A2_sets_upstream_datum_to_zero': False,
        'leading_paired_e_minus2_cancels_for_equal_mu': True,
        'leading_e_minus2_cancellation_is_not_this_value': True,
        'first_noncancelling_paired_order': FIRST_NONCANCELLING_PAIRED_ORDER,
    })


@lru_cache(maxsize=1)
def reduced_massless_transport_identities():
    """Eliminate the linked h3/n4 variables from the cubic current equation.

    This is an exact algebraic reduction of the constructed transport, not a
    new quadrature or a tail certificate. It leaves one triangular Volterra
    system for q_z, q_zz and J_I,3, with only geometry-dependent coefficients.
    """
    rho, z = sp.symbols('rho z', real=True)
    a = sp.Function('a', positive=True)(rho)
    b = sp.Function('b', real=True)(rho,z)
    q = sp.Function('q', real=True)(rho,z)
    k = sp.symbols('k',real=True)
    checks = {}
    for s in SOURCE_SIGNS:
        T = lambda f: sp.diff(f,rho)-s/a**2*sp.diff(f,z)
        M = lambda f: sp.diff(f,rho)+s/a**2*sp.diff(f,z)
        B = s*b*q-s*a**2*sp.diff(b,rho)/2-sp.diff(b,z)/2
        R = (s*sp.diff(q,z)+k-q*q-b*b)/2
        depth = a**2*M(B)/2+s*b*R
        qr = -2*b*b/a**2+s/a**2*sp.diff(q,z)
        substitute = {sp.diff(q,rho,1,z,j):sp.diff(qr,z,j) for j in range(3)}
        depth = sp.expand(depth).subs(substitute)
        rest = R*sp.diff(q,z)-q*sp.diff(R,z)+s*(b*sp.diff(B,z)-B*sp.diff(b,z))
        # Continuity + major A3 give T_s(h3_z-s n4)=2/a^2 d_z(B^2+s b D).
        generated = sp.expand(2/a**2*sp.diff(B*B+s*b*depth,z)+T(rest)).subs(substitute)
        reduced = (-a*a*b*sp.diff(b,rho,2,z,1)
                   +a*a*sp.diff(b,rho)*sp.diff(b,rho,z)
                   -2*a*sp.diff(a,rho)*b*sp.diff(b,rho,z)
                   -s*b*sp.diff(b,rho,1,z,2)+s*sp.diff(b,z)*sp.diff(b,rho,z)
                   -16*b**3*sp.diff(b,z)/a**2+4*s*b*b*sp.diff(q,z,2)/a**2)
        checks['s%+d/Ts_J3' % s] = _zero(generated-reduced)

        aa,bb,bz,qq,qz,qzz,kk,j,jz,DD,hz,nn,hh = sp.symbols(
            'a b bz q qz qzz k j jz D h3z n4 h',real=True)
        RR=(s*qz+kk-qq*qq-bb*bb)/2
        Rz=(s*qzz-2*qq*qz-2*bb*bz)/2
        BB=s*bb*qq+j; Bz=s*bz*qq+s*bb*qz+jz
        ai=[sp.zeros(2,1) for _ in range(4)]
        iz=[sp.zeros(2,1) for _ in range(4)]
        major,minor=_major_index(s),_minor_index(s)
        ai[0][major]=1; ai[1][major]=sp.I*qq; ai[1][minor]=-sp.I*s*bb
        ai[2][major]=RR+sp.I*hh; ai[2][minor]=BB; ai[3][minor]=s*bb*hh-sp.I*DD
        iz[1][major]=sp.I*qz; iz[1][minor]=-sp.I*s*bz
        iz[2][major]=Rz; iz[2][minor]=Bz; iz[3][major]=sp.I*hz
        _,sigma2,sigma3=_sym_pauli()
        V=2*bb/aa*sigma2
        v3=sum((ai[n].H*V*ai[3-n])[0] for n in range(4))
        current=sp.im(sp.expand_complex(sum((ai[n].H*sigma3*iz[3-n])[0]
                                            for n in range(4))))-nn+2*(BB*BB+2*s*bb*DD)
        J3=hz+RR*qz-qq*Rz+s*(bb*Bz-BB*bz)-s*nn
        remainder=2/aa*(j*j-2*s*bb*bb*qz-bb*bb*kk+bb**4-bb*jz+j*bz)
        checks['s%+d/N_full_reduction' % s]=_zero(v3+current/aa-s*J3/aa-remainder)
        bp,bp0,bpz,qz0=sp.symbols('bp bp0 bpz qz0',real=True)
        surface = remainder.subs({bz:0,j:-s*aa*aa*bp/2,jz:-s*aa*aa*bpz/2})
        baseline = remainder.subs({bz:0,j:-s*aa*aa*bp0/2,jz:0,qz:qz0})
        difference=aa**3*(bp*bp-bp0*bp0)/2-4*s*bb*bb/aa*(qz-qz0)+s*aa*bb*bpz
        checks['s%+d/N_Sigma_difference' % s]=_zero(surface-baseline-difference)
    return MappingProxyType({
        'residuals':checks,
        'b':'a ell/(2r)',
        'Ts_qz':'-4 b b_z/a^2',
        'Ts_qzz':'-4 (b_z^2+b b_zz)/a^2',
        'Ts_J3':'-a^2 b b_rhorhoz+a^2 b_rho b_rhoz-2 a a_rho b b_rhoz-s b b_rhozz+s b_z b_rhoz-16 b^3 b_z/a^2+4 s b^2 q_zz/a^2',
        'Sigma_delta_N_bracket':'s delta(J3)/a+a^3 delta(b_rho^2)/2-4 s b^2 delta(q_z)/a+s a b delta(b_rhoz)',
        'common_k_cancels_on_Sigma':True,
        'required_transport_quantities':['q_z','q_zz','J3'],
        'transport_type':'triangular affine Volterra; geometry-only coefficients',
        'directed_quadrature_bound':None,
    })


def _axial_rho_rho(rho, axial):
    theta = 0.5 * np.pi - np.arctan(float(rho))
    u = 1.0 - float(rho) * theta
    u_rho = -theta + float(rho) / (1.0 + float(rho) * float(rho))
    return -3.0 * u_rho / axial - 9.0 * u * u / axial**3


def _normal_window_jets(normal, inner, outer):
    """Even C-infinity plateau chi(s) through the second s-derivative."""
    normal = float(normal)
    inner = float(inner)
    outer = float(outer)
    distance = abs(normal)
    sign = 1.0 if normal >= 0.0 else -1.0
    if distance <= inner:
        return 1.0, 0.0, 0.0
    if distance >= outer:
        return 0.0, 0.0, 0.0
    width = outer - inner
    toward = distance - inner
    away = outer - distance
    u = toward / width
    if u <= 1e-12:
        return 1.0, 0.0, 0.0
    if u >= 1.0 - 1e-12:
        return 0.0, 0.0, 0.0
    ratio = 1.0 / (1.0 - u) - 1.0 / u
    eta = float(expit(-ratio))
    # Match the owned plateau: 1-eta rounds to zero before the true
    # derivative does on the inner ramp. Evaluate its complement directly.
    product = eta * float(expit(ratio))
    q1 = (1.0 - u)**-2 + u**-2
    q2 = 2.0 * ((1.0 - u)**-3 - u**-3)
    d_eta_du = -product * q1
    d2_eta_du2 = product * (1.0 - 2.0 * eta) * q1 * q1 - product * q2
    return eta, sign * d_eta_du / width, d2_eta_du2 / width**2


def geometry_jets(family, rho, z):
    """Owned (a,r) jets through the mixed orders used by cubic transport.

    These are coordinate partials of the declared radius history. They are
    not a field, not C4, and not a remainder bound.
    """
    metric = _as_one_direction(family)
    rho = float(rho)
    z = np.asarray(z, dtype=float)
    if z.ndim != 1 or z.size == 0 or not np.isfinite(z).all():
        raise ValueError('nonempty finite one-dimensional z required')
    _, _, axial, axial_rho, _ = reference_chart(rho)
    if axial <= 0:
        raise ValueError('strictly trapped a>0 required')
    a_rr = _axial_rho_rho(rho, axial)
    r_ref = float(np.sqrt(1.0 + rho * rho))
    r_ref_rho = rho / r_ref
    r_ref_rr = 1.0 / r_ref**3
    normal = _normal_coordinate(rho)
    normal_rho = -1.0 / axial
    normal_rr = axial_rho / axial**2
    delta = np.zeros(z.shape, dtype=float)
    delta_z = np.zeros_like(delta)
    delta_s = np.zeros_like(delta)
    delta_zz = np.zeros_like(delta)
    delta_sz = np.zeros_like(delta)
    delta_ss = np.zeros_like(delta)
    delta_zzz = np.zeros_like(delta)
    delta_szz = np.zeros_like(delta)
    delta_ssz = np.zeros_like(delta)
    for amplitude, direction in zip(metric.amplitudes, metric.directions):
        window, window_s, window_ss = _normal_window_jets(
            normal, direction.inner_radius, direction.outer_radius)
        if window == 0.0 and window_s == 0.0 and window_ss == 0.0:
            continue
        w0 = np.array([direction.w(float(point), 0) for point in z], dtype=float)
        w1 = np.array([direction.w(float(point), 1) for point in z], dtype=float)
        w2 = np.array([direction.w(float(point), 2) for point in z], dtype=float)
        w3 = np.array([direction.w(float(point), 3) for point in z], dtype=float)
        u0 = np.array([direction.U(float(point), 0) for point in z], dtype=float)
        u1 = np.array([direction.U(float(point), 1) for point in z], dtype=float)
        u2 = np.array([direction.U(float(point), 2) for point in z], dtype=float)
        u3 = np.array([direction.U(float(point), 3) for point in z], dtype=float)
        poly = normal * w0 + (normal**3) * u0 / 6.0
        poly_s = w0 + (normal**2) * u0 / 2.0
        poly_ss = normal * u0
        poly_z = normal * w1 + (normal**3) * u1 / 6.0
        poly_sz = w1 + (normal**2) * u1 / 2.0
        poly_ssz = normal * u1
        poly_zz = normal * w2 + (normal**3) * u2 / 6.0
        poly_szz = w2 + (normal**2) * u2 / 2.0
        poly_zzz = normal * w3 + (normal**3) * u3 / 6.0
        scale = float(amplitude)
        delta += scale * window * poly
        delta_z += scale * window * poly_z
        delta_s += scale * (window_s * poly + window * poly_s)
        delta_zz += scale * window * poly_zz
        delta_sz += scale * (window_s * poly_z + window * poly_sz)
        delta_ss += scale * (
            window_ss * poly + 2.0 * window_s * poly_s + window * poly_ss)
        delta_zzz += scale * window * poly_zzz
        delta_szz += scale * (window_s * poly_zz + window * poly_szz)
        delta_ssz += scale * (
            window_ss * poly_z + 2.0 * window_s * poly_sz + window * poly_ssz)
    radius = r_ref + delta
    if not np.isfinite(radius).all() or np.any(radius <= 0):
        raise ValueError('actual radius must stay finite and positive')
    return {
        'a': axial,
        'a_rho': axial_rho,
        'a_rho_rho': a_rr,
        'r': radius,
        'r_z': delta_z,
        'r_rho': r_ref_rho + normal_rho * delta_s,
        'r_zz': delta_zz,
        'r_zzz': delta_zzz,
        'r_rhoz': normal_rho * delta_sz,
        'r_rhozz': normal_rho * delta_szz,
        'r_rho_rho': r_ref_rr + normal_rr * delta_s + normal_rho**2 * delta_ss,
        'r_rho_rhoz': normal_rr * delta_sz + normal_rho**2 * delta_ssz,
        'r_ref': r_ref,
        's': normal,
        's_rho': normal_rho,
    }


def _distance(rho):
    rho = float(rho)
    if abs(rho - RHO_SIGMA) <= 1e-15:
        return 0.0
    return float(_characteristic_distance(rho))


def massless_B_jets(jets, phase, phase_z, phase_rho, phase_zz, angular, source_sign):
    """Minor A2 = B and its coordinate jets from the massless parent formula."""
    axial = jets['a']
    a_rho = jets['a_rho']
    a_rr = jets['a_rho_rho']
    radius = jets['r']
    r_rho = jets['r_rho']
    r_z = jets['r_z']
    r_zz = jets['r_zz']
    r_zzz = jets['r_zzz']
    r_rhoz = jets['r_rhoz']
    r_rhozz = jets['r_rhozz']
    r_rr = jets['r_rho_rho']
    r_rrz = jets['r_rho_rhoz']
    poly = (
        -a_rho * axial * radius
        + axial**2 * r_rho
        + source_sign * r_z
        + 2.0 * phase * radius)
    pre = source_sign * angular * axial / 4.0
    pre_rho = source_sign * angular * a_rho / 4.0
    u = poly / radius**2
    poly_z = (
        -a_rho * axial * r_z
        + axial**2 * r_rhoz
        + source_sign * r_zz
        + 2.0 * phase_z * radius
        + 2.0 * phase * r_z)
    u_z = poly_z / radius**2 - 2.0 * poly * r_z / radius**3
    poly_zz = (
        -a_rho * axial * r_zz
        + axial**2 * r_rhozz
        + source_sign * r_zzz
        + 2.0 * phase_zz * radius
        + 4.0 * phase_z * r_z
        + 2.0 * phase * r_zz)
    u_zz = (
        poly_zz / radius**2
        - 4.0 * poly_z * r_z / radius**3
        - 2.0 * poly * r_zz / radius**3
        + 6.0 * poly * r_z**2 / radius**4)
    poly_rho = (
        -a_rr * axial * radius
        - a_rho**2 * radius
        + axial * a_rho * r_rho
        + axial**2 * r_rr
        + source_sign * r_rhoz
        + 2.0 * phase_rho * radius
        + 2.0 * phase * r_rho)
    phase_rhoz = (
        angular**2 * r_z / radius**3
        + source_sign / axial**2 * phase_zz)
    poly_rhoz = (
        -a_rr * axial * r_z
        - a_rho**2 * r_z
        + axial * a_rho * r_rhoz
        + axial**2 * r_rrz
        + source_sign * r_rhozz
        + 2.0 * phase_rhoz * radius
        + 2.0 * phase_rho * r_z
        + 2.0 * phase_z * r_rho
        + 2.0 * phase * r_rhoz)
    u_rho = poly_rho / radius**2 - 2.0 * poly * r_rho / radius**3
    u_rhoz = (
        poly_rhoz / radius**2
        - 2.0 * poly_rho * r_z / radius**3
        - 2.0 * poly_z * r_rho / radius**3
        - 2.0 * poly * r_rhoz / radius**3
        + 6.0 * poly * r_rho * r_z / radius**4)
    minor_B = pre * u
    B_z = pre * u_z
    B_zz = pre * u_zz
    B_rho = pre_rho * u + pre * u_rho
    B_rhoz = pre_rho * u_z + pre * u_rhoz
    b = axial * angular / (2.0 * radius)
    b_z = -b * r_z / radius
    if not np.isfinite(minor_B).all() or not np.isfinite(B_z).all():
        raise ArithmeticError('nonfinite massless minor A2 jets')
    return minor_B, B_z, B_zz, B_rho, B_rhoz, b, b_z


def _on_shell_q_rho(jets, phase_z, angular, source_sign):
    return (
        -angular**2 / (2.0 * jets['r']**2)
        + source_sign / jets['a']**2 * phase_z)


def _massless_generated_bundle(
        jets, phase, phase_z, phase_zz, angular, source_sign, k, height):
    q_rho = _on_shell_q_rho(jets, phase_z, angular, source_sign)
    minor_B, B_z, B_zz, B_rho, B_rhoz, b, b_z = massless_B_jets(
        jets, phase, phase_z, q_rho, phase_zz, angular, source_sign)
    inv_a2 = 1.0 / jets['a']**2
    tminus_B = B_rho + source_sign * inv_a2 * B_z
    tminus_B_z = B_rhoz + source_sign * inv_a2 * B_zz
    real_major = 0.5 * (
        source_sign * phase_z + k - phase**2 - b**2)
    real_major_z = 0.5 * (
        source_sign * phase_zz - 2.0 * phase * phase_z - 2.0 * b * b_z)
    depth = 0.5 * jets['a']**2 * tminus_B + source_sign * b * real_major
    depth_z = (
        0.5 * jets['a']**2 * tminus_B_z
        + source_sign * (b_z * real_major + b * real_major_z))
    minor_C = source_sign * b * height - 1j * depth
    ts_h3 = -source_sign * b * tminus_B - 2.0 * b**2 * real_major * inv_a2
    ts_h3_z = (
        -source_sign * (b_z * tminus_B + b * tminus_B_z)
        - 4.0 * b * b_z * real_major * inv_a2
        - 2.0 * b**2 * real_major_z * inv_a2)
    re_lc = source_sign * b * depth
    paren = 2.0 * minor_B**2 + 4.0 * re_lc
    paren_z = (
        4.0 * minor_B * B_z
        + 4.0 * source_sign * (b_z * depth + b * depth_z))
    return {
        'B': minor_B, 'B_z': B_z, 'B_rho': B_rho, 'b': b, 'b_z': b_z,
        'Tminus_B': tminus_B, 'R': real_major, 'R_z': real_major_z,
        'C': minor_C, 'T_s_h3': ts_h3, 'T_s_h3_z': ts_h3_z,
        're_lc': re_lc, 'n4_paren': paren, 'n4_paren_z': paren_z,
        'q_rho': q_rho,
    }


def _transport_sources(family, rho, z_in, phase, phase_z, phase_zz, angular,
                       source_sign, k, height):
    distance = _distance(rho)
    z_phys = z_in - source_sign * distance
    jets = geometry_jets(family, rho, z_phys)
    bundle = _massless_generated_bundle(
        jets, phase, phase_z, phase_zz, angular, source_sign, k, height)
    ell2 = angular * angular
    radius = jets['r']
    ts_q = -ell2 / (2.0 * radius**2)
    ts_qz = ell2 * jets['r_z'] / radius**3
    ts_qzz = ell2 * (
        jets['r_zz'] / radius**3 - 3.0 * jets['r_z']**2 / radius**4)
    ts_h3z = bundle['T_s_h3_z']
    ts_n4 = -(source_sign / jets['a']**2) * bundle['n4_paren_z']
    return ts_q, ts_qz, ts_qzz, ts_h3z, ts_n4


def transport_massless_cubic_fields(
        family, z, *, angular, source_sign, rho_up=ARCHIVE_RHO_UP,
        rho_steps=DEFAULT_RHO_STEPS, k=0.0, height=0.0):
    """RK4 T_s march of q jets, h3_z and n4 from z-homogeneous upstream data.

    Representative ICs: q_up=q_z_up=q_zz_up=0, h3_z_up=0, n4_up=0, common k
    and common massless h. R is always (s q_z+k-q^2-b^2)/2; it is not zeroed
    independently of h3 and n4. Quadrature error is uncertified.
    """
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError('owned LocalIncomingFamily required')
    if source_sign not in SOURCE_SIGNS:
        raise ValueError('source sign must be +1 or -1')
    _, angular = _require_massless_angular(0.0, angular)
    z = np.asarray(z, dtype=float)
    if z.ndim != 1 or z.size < 5 or not np.isfinite(z).all():
        raise ValueError('uniform z grid with at least five finite nodes required')
    step = float(z[1] - z[0])
    if step <= 0 or not np.allclose(np.diff(z), step):
        raise ValueError('uniform increasing z grid required')
    rho_up = float(rho_up)
    if rho_up <= RHO_SIGMA:
        raise ValueError('upstream rho must exceed the incoming slice')
    if isinstance(rho_steps, (bool, np.bool_)) or int(rho_steps) != rho_steps or rho_steps < 2:
        raise ValueError('at least two characteristic steps required')
    rho_steps = int(rho_steps)
    k = float(k)
    height = float(height)
    rhos = np.linspace(rho_up, RHO_SIGMA, rho_steps + 1)
    q = np.zeros_like(z)
    qz = np.zeros_like(z)
    qzz = np.zeros_like(z)
    h3z = np.zeros_like(z)
    n4 = np.zeros_like(z)

    def add(state, source, scale):
        return tuple(item + scale * term for item, term in zip(state, source))

    for rho0, rho1 in zip(rhos[:-1], rhos[1:]):
        d_rho = float(rho1 - rho0)
        state = (q, qz, qzz, h3z, n4)
        k1 = _transport_sources(
            family, rho0, z, state[0], state[1], state[2], angular,
            source_sign, k, height)
        s2 = add(state, k1, 0.5 * d_rho)
        k2 = _transport_sources(
            family, rho0 + 0.5 * d_rho, z, s2[0], s2[1], s2[2], angular,
            source_sign, k, height)
        s3 = add(state, k2, 0.5 * d_rho)
        k3 = _transport_sources(
            family, rho0 + 0.5 * d_rho, z, s3[0], s3[1], s3[2], angular,
            source_sign, k, height)
        s4 = add(state, k3, d_rho)
        k4 = _transport_sources(
            family, rho1, z, s4[0], s4[1], s4[2], angular, source_sign,
            k, height)
        q, qz, qzz, h3z, n4 = tuple(
            item + d_rho * (p1 + 2.0 * p2 + 2.0 * p3 + p4) / 6.0
            for item, p1, p2, p3, p4 in zip(state, k1, k2, k3, k4))
    if not all(np.isfinite(item).all() for item in (q, qz, qzz, h3z, n4)):
        raise ArithmeticError('nonfinite massless cubic transport')
    return MappingProxyType({
        'z': np.array(z, dtype=float, copy=True),
        'q': q, 'q_z': qz, 'q_zz': qzz, 'h3_z': h3z, 'n4': n4,
        'source_sign': int(source_sign),
        'angular': angular,
        'k': k,
        'h': height,
        'rho_up': rho_up,
        'rho_steps': rho_steps,
        'certified_quadrature_error_bound': None,
        'R_from_continuity_and_k': True,
        'R_zeroed_independently': False,
    })


def _spinor_columns(source_sign, q, b, minor_B, real_major, height, minor_C):
    major = _major_index(source_sign)
    minor = _minor_index(source_sign)
    a0 = np.zeros((2, q.size), dtype=complex)
    a0[major] = 1.0
    a1 = np.zeros_like(a0)
    a1[major] = 1j * q
    a1[minor] = -1j * source_sign * b
    a2 = np.zeros_like(a0)
    a2[major] = real_major + 1j * height
    a2[minor] = minor_B
    a3 = np.zeros_like(a0)
    a3[minor] = minor_C
    return a0, a1, a2, a3


def _spinor_z(source_sign, qz, bz, Bz, Rz, h3z):
    major = _major_index(source_sign)
    minor = _minor_index(source_sign)
    a1z = np.zeros((2, qz.size), dtype=complex)
    a1z[major] = 1j * qz
    a1z[minor] = -1j * source_sign * bz
    a2z = np.zeros_like(a1z)
    a2z[major] = Rz
    a2z[minor] = Bz
    a3z = np.zeros_like(a1z)
    a3z[major] = 1j * h3z
    return a1z, a2z, a3z


def _inner(left, right):
    return np.einsum('i...,i...->...', np.conj(left), right)


def cubic_pauli_densities(jets, fields, angular, source_sign, k=0.0, height=0.0):
    """Incoming e^{-3} J_I, N-bracket and action densities from actual vertices."""
    _, angular = _require_massless_angular(0.0, angular)
    if source_sign not in SOURCE_SIGNS:
        raise ValueError('source sign must be +1 or -1')
    q = np.asarray(fields['q'], dtype=float)
    qz = np.asarray(fields['q_z'], dtype=float)
    qzz = np.asarray(fields['q_zz'], dtype=float)
    h3z = np.asarray(fields['h3_z'], dtype=float)
    n4 = np.asarray(fields['n4'], dtype=float)
    bundle = _massless_generated_bundle(
        jets, q, qz, qzz, angular, source_sign, k, height)
    b = bundle['b']
    bz = bundle['b_z']
    real_major = bundle['R']
    rz = bundle['R_z']
    a0, a1, a2, a3 = _spinor_columns(
        source_sign, q, b, bundle['B'], real_major, height, bundle['C'])
    a1z, a2z, a3z = _spinor_z(source_sign, qz, bz, bundle['B_z'], rz, h3z)
    j_i = np.imag(_inner(a0, a3z) + _inner(a1, a2z) + _inner(a2, a1z)) - source_sign * n4
    s3 = np.array([1.0, -1.0])
    j_s3_interf = np.imag(
        _inner(a0, s3[:, None] * a3z)
        + _inner(a1, s3[:, None] * a2z)
        + _inner(a2, s3[:, None] * a1z))
    minor_sq4 = bundle['B']**2 + 2.0 * bundle['re_lc']
    as3a4 = source_sign * (n4 - 2.0 * minor_sq4)
    j_s3 = j_s3_interf - source_sign * as3a4
    ell_over_r = angular / jets['r']
    # V = ell/r S2, S2 = [[0,-i],[i,0]].
    v_a3 = np.zeros_like(a3)
    v_a3[0] = -1j * ell_over_r * a3[1]
    v_a3[1] = 1j * ell_over_r * a3[0]
    v_a2 = np.zeros_like(a2)
    v_a2[0] = -1j * ell_over_r * a2[1]
    v_a2[1] = 1j * ell_over_r * a2[0]
    v_a1 = np.zeros_like(a1)
    v_a1[0] = -1j * ell_over_r * a1[1]
    v_a1[1] = 1j * ell_over_r * a1[0]
    v_a0 = np.zeros_like(a0)
    v_a0[0] = -1j * ell_over_r * a0[1]
    v_a0[1] = 1j * ell_over_r * a0[0]
    v3 = np.real(_inner(a0, v_a3) + _inner(a1, v_a2) + _inner(a2, v_a1) + _inner(a3, v_a0))
    n_bracket = v3 + j_s3 / jets['a']
    parent_ji = (
        h3z + real_major * qz - q * rz
        + source_sign * (b * bundle['B_z'] - bundle['B'] * bz)
        - source_sign * n4)
    if np.max(np.abs(j_i - parent_ji)) > 1e-8 * (1.0 + np.max(np.abs(parent_ji))):
        raise ArithmeticError('spinor J_I_3 drifted from the parent formula')
    action_n = -n_bracket / TWO_PI
    action_beta = j_i / TWO_PI
    return MappingProxyType({
        'J_I_3': j_i,
        'J_S3_3': j_s3,
        'V_3': v3,
        'N_bracket_3': n_bracket,
        'action_N3_over_mu': action_n,
        'action_beta3_over_mu': action_beta,
        'R': real_major,
        'B': bundle['B'],
        'C': bundle['C'],
        'b': b,
        'parent_J_I_3': parent_ji,
    })


def _work_grid(family, rho_up, n_z):
    """Uniform nodes on I. Incoming r_z vanishes; z-jets of sources are algebraic."""
    del rho_up
    left, right = family.interval
    return np.linspace(float(left), float(right), int(n_z))


def _restrict_interval(z, values, interval):
    inner_left, inner_right = interval
    mask = (z >= inner_left) & (z <= inner_right)
    if not np.any(mask):
        raise ArithmeticError('working grid does not cover I')
    return z[mask], values[mask]


def _trapz(values, coords):
    if hasattr(np, 'trapezoid'):
        return np.trapezoid(values, coords)
    return np.trapz(values, coords)


def _integrate_on_interval(z, values, interval):
    z_i, v_i = _restrict_interval(z, values, interval)
    return float(_trapz(v_i, z_i)), float(np.max(np.abs(v_i))), z_i, v_i


def evaluate_massless_cubic_on_history(
        family, *, angular, rho_up=ARCHIVE_RHO_UP, rho_steps=DEFAULT_RHO_STEPS,
        z_points=DEFAULT_Z_POINTS, k=0.0, height=0.0, multiplicity=1.0):
    """History-minus-reference cubic action densities on I, both signs and ell.

    The angular pair is even at e^{-3}; both original signs are still summed.
    Multiplicity mu multiplies the mu/(2 pi) action densities. Quadrature is
    diagnostic: certified_quadrature_error_bound remains None.
    """
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError('owned LocalIncomingFamily required')
    _, angular = _require_massless_angular(0.0, angular)
    multiplicity = float(multiplicity)
    if not np.isfinite(multiplicity) or multiplicity <= 0:
        raise ValueError('positive multiplicity required')
    z = _work_grid(family, rho_up, z_points)
    reference = LocalIncomingFamily(np.zeros_like(np.asarray(family.coefficients)))
    interval = family.interval
    action_n = np.zeros_like(z)
    action_beta = np.zeros_like(z)
    odd_n = np.zeros_like(z)
    odd_beta = np.zeros_like(z)
    channels = []
    for ell in (angular, -angular):
        ell_n = np.zeros_like(z)
        ell_beta = np.zeros_like(z)
        for source_sign in SOURCE_SIGNS:
            transported = transport_massless_cubic_fields(
                family, z, angular=ell, source_sign=source_sign,
                rho_up=rho_up, rho_steps=rho_steps, k=k, height=height)
            transported_ref = transport_massless_cubic_fields(
                reference, z, angular=ell, source_sign=source_sign,
                rho_up=rho_up, rho_steps=rho_steps, k=k, height=height)
            jets = geometry_jets(family, RHO_SIGMA, z)
            jets_ref = geometry_jets(reference, RHO_SIGMA, z)
            current = cubic_pauli_densities(
                jets, transported, ell, source_sign, k=k, height=height)
            ref = cubic_pauli_densities(
                jets_ref, transported_ref, ell, source_sign, k=k, height=height)
            delta_n = multiplicity * (
                current['action_N3_over_mu'] - ref['action_N3_over_mu'])
            delta_beta = multiplicity * (
                current['action_beta3_over_mu'] - ref['action_beta3_over_mu'])
            ell_n += delta_n
            ell_beta += delta_beta
            channels.append({
                'angular_sign': int(np.sign(ell)),
                'source_sign': int(source_sign),
                'I_max_abs_N': float(np.max(np.abs(
                    _restrict_interval(z, delta_n, interval)[1]))),
                'I_max_abs_beta': float(np.max(np.abs(
                    _restrict_interval(z, delta_beta, interval)[1]))),
            })
        action_n += ell_n
        action_beta += ell_beta
        if ell == angular:
            plus_n, plus_beta = ell_n, ell_beta
        else:
            odd_n = 0.5 * (plus_n - ell_n)
            odd_beta = 0.5 * (plus_beta - ell_beta)
    n_int, n_max, z_i, n_i = _integrate_on_interval(z, action_n, interval)
    b_int, b_max, _, b_i = _integrate_on_interval(z, action_beta, interval)
    odd_n_max = float(np.max(np.abs(_restrict_interval(z, odd_n, interval)[1])))
    odd_b_max = float(np.max(np.abs(_restrict_interval(z, odd_beta, interval)[1])))
    return MappingProxyType({
        'z_I': z_i,
        'action_N3': n_i,
        'action_beta3': b_i,
        'I_integral_N': n_int,
        'I_integral_beta': b_int,
        'I_max_abs_N': n_max,
        'I_max_abs_beta': b_max,
        'ell_odd_max_abs_N': odd_n_max,
        'ell_odd_max_abs_beta': odd_b_max,
        'channels': tuple(channels),
        'k': float(k),
        'h': float(height),
        'multiplicity': multiplicity,
        'rho_steps': int(rho_steps),
        'z_points': int(z_points),
        'certified_quadrature_error_bound': None,
        'validated': False,
        'R_from_continuity_and_k': True,
        'R_zeroed_independently': False,
        'measure': 'mu/(2 pi); N=-N_bracket, beta=+J_I',
    })


def diagnostic_cubic_quadrature(
        family, *, angular, multiplicity, rho_up=ARCHIVE_RHO_UP,
        k=0.0, height=0.0):
    """Fine/coarse current-history measurement; remainder-unvalidated."""
    fine = evaluate_massless_cubic_on_history(
        family, angular=angular, rho_up=rho_up, rho_steps=DEFAULT_RHO_STEPS,
        z_points=DEFAULT_Z_POINTS, k=k, height=height,
        multiplicity=multiplicity)
    coarse = evaluate_massless_cubic_on_history(
        family, angular=angular, rho_up=rho_up, rho_steps=COARSE_RHO_STEPS,
        z_points=COARSE_Z_POINTS, k=k, height=height,
        multiplicity=multiplicity)
    if fine['z_I'].shape != coarse['z_I'].shape or not np.allclose(
            fine['z_I'], coarse['z_I']):
        n_cmp = np.interp(fine['z_I'], coarse['z_I'], coarse['action_N3'])
        b_cmp = np.interp(fine['z_I'], coarse['z_I'], coarse['action_beta3'])
    else:
        n_cmp = coarse['action_N3']
        b_cmp = coarse['action_beta3']
    n_disc = float(np.max(np.abs(fine['action_N3'] - n_cmp)))
    b_disc = float(np.max(np.abs(fine['action_beta3'] - b_cmp)))
    return MappingProxyType({
        'validated': False,
        'certified_quadrature_error_bound': None,
        'remainder_proof_supplied': False,
        'I_integral_N': fine['I_integral_N'],
        'I_integral_beta': fine['I_integral_beta'],
        'I_max_abs_N': fine['I_max_abs_N'],
        'I_max_abs_beta': fine['I_max_abs_beta'],
        'I_center_N': float(fine['action_N3'][len(fine['action_N3']) // 2]),
        'I_center_beta': float(fine['action_beta3'][len(fine['action_beta3']) // 2]),
        'z_I': [float(value) for value in fine['z_I']],
        'action_N3': [float(value) for value in fine['action_N3']],
        'action_beta3': [float(value) for value in fine['action_beta3']],
        'ell_odd_max_abs_N': fine['ell_odd_max_abs_N'],
        'ell_odd_max_abs_beta': fine['ell_odd_max_abs_beta'],
        'refinement_discrepancy_N': n_disc,
        'refinement_discrepancy_beta': b_disc,
        'fine_rho_steps': DEFAULT_RHO_STEPS,
        'fine_z_points': DEFAULT_Z_POINTS,
        'coarse_rho_steps': COARSE_RHO_STEPS,
        'coarse_z_points': COARSE_Z_POINTS,
        'k_representative': float(k),
        'h_representative': float(height),
        'R_from_continuity_and_k': True,
        'R_zeroed_independently': False,
        'ieee754_hex': {
            'I_integral_N': _hex(fine['I_integral_N']),
            'I_integral_beta': _hex(fine['I_integral_beta']),
            'I_max_abs_N': _hex(fine['I_max_abs_N']),
            'I_max_abs_beta': _hex(fine['I_max_abs_beta']),
            'refinement_discrepancy_N': _hex(n_disc),
            'refinement_discrepancy_beta': _hex(b_disc),
        },
        'sample_used_as_proof': False,
        'fitted_decay': False,
    })


def unlinked_R_negative_control():
    """Holding h3_z and n4 fixed while shifting R leaves +c q_z in J_I_3."""
    identities = massless_cubic_identities()
    for source_sign in SOURCE_SIGNS:
        key = 's%+d/unlinked_R_spurious_c_qz' % source_sign
        if identities['residuals'][key] != '0':
            raise ArithmeticError('unlinked R control is not the parent identity')
    return MappingProxyType({
        'unlinked_R_shift_of_J_I_3': 'c q_z',
        'linked_scalar_shift_of_J_I_3': 0,
        'R_may_not_be_zeroed_independently': True,
    })


def current_history_cubic_uv_report(root=None):
    """Executable generated cubic coefficient on live 0b0e4ced, group 1."""
    root = Path(root) if root is not None else _root()
    history = json.loads((root / HISTORY_RECORD).read_text())
    family = LocalIncomingFamily(np.array(history['history']['coefficients']))
    identity = profile_identity(family)
    if identity != CURRENT_HISTORY_IDENTITY:
        raise ValueError('expected exact live history 0b0e4ced')
    pair = smallest_original_pair(root)
    if float(pair['mass']) != 0.0:
        raise ValueError('this cubic owner is the massless original pair')
    identities = massless_cubic_identities()
    signs = confirm_action_signs_against_source_column_matter()
    zero_family = LocalIncomingFamily(np.zeros_like(np.asarray(family.coefficients)))
    zero = evaluate_massless_cubic_on_history(
        zero_family, angular=pair['absolute_angular'],
        multiplicity=pair['multiplicity_per_signed_family'],
        rho_steps=COARSE_RHO_STEPS, z_points=DEFAULT_Z_POINTS)
    diagnostic = diagnostic_cubic_quadrature(
        family, angular=pair['absolute_angular'],
        multiplicity=pair['multiplicity_per_signed_family'])
    k_base = evaluate_massless_cubic_on_history(
        family, angular=pair['absolute_angular'],
        multiplicity=pair['multiplicity_per_signed_family'],
        rho_steps=COARSE_RHO_STEPS, z_points=COARSE_Z_POINTS, k=0.0)
    k_shift = evaluate_massless_cubic_on_history(
        family, angular=pair['absolute_angular'],
        multiplicity=pair['multiplicity_per_signed_family'],
        rho_steps=COARSE_RHO_STEPS, z_points=COARSE_Z_POINTS, k=0.25)
    k_disc_n = abs(k_shift['I_integral_N'] - k_base['I_integral_N'])
    k_disc_b = abs(k_shift['I_integral_beta'] - k_base['I_integral_beta'])
    upstream = vacuum_upstream_conditions()
    return MappingProxyType({
        'schema': SCHEMA,
        'profile_identity': identity,
        'pair': dict(pair),
        'identities': identities,
        'reduced_transport': reduced_massless_transport_identities(),
        'action_signs': dict(signs),
        'first_noncancelling_paired_order': FIRST_NONCANCELLING_PAIRED_ORDER,
        'coefficient_scope': (
            'history-minus-reference occupied-vacuum e^{-3} N,beta contraction '
            'on I of the original equal-mu massless pair, action minus once, '
            'measure mu/(2 pi)'),
        'representative': {
            'k': 0.0,
            'h_up': 0.0,
            'h3_z_up': 0.0,
            'n4_up': 0.0,
            'q_up': 0.0,
            'R': '(s q_z + k - q^2 - b^2)/2, never independently zero',
            'C_from_recurrence': True,
            'common_scalar_dropped_by_invariance': True,
            'n2_up_claimed_physical_zero': False,
        },
        'zero_history_control': {
            'I_max_abs_N': zero['I_max_abs_N'],
            'I_max_abs_beta': zero['I_max_abs_beta'],
            'I_integral_N': zero['I_integral_N'],
            'I_integral_beta': zero['I_integral_beta'],
        },
        'common_k_shift_control': {
            'delta_k': 0.25,
            'I_integral_N_discrepancy': float(k_disc_n),
            'I_integral_beta_discrepancy': float(k_disc_b),
            'certified': False,
        },
        'unlinked_R_control': dict(unlinked_R_negative_control()),
        'diagnostic_quadrature': diagnostic,
        'numerical_C4': None,
        'numerical_C_M': None,
        'vacuum_integrated_tail_N_beta': None,
        'validated_cubic_coefficient': False,
        'certified_quadrature_error_bound': None,
        'production_A4_constructed': False,
        'generated_massless_A3_n4_constructed': True,
        'dummy_symbol_L0_A_j': False,
        'A3_locally_algebraic_in_metric_jets': False,
        'R_zeroed_while_freezing_h3_n4': False,
        'massless_Im_Ts_major_A2_sets_upstream_datum_to_zero': False,
        'leading_paired_e_minus2_cancels_for_equal_mu': True,
        'leading_e_minus2_cancellation_is_not_this_value': True,
        'e_minus3_first_allowed_order_is_not_a_tail_bound': True,
        'source_columns_changed': False,
        'upstream': {
            'A0': upstream['A0'],
            'major_A1': upstream['major_A1'],
            'history_flat_at_upstream': upstream['history_flat_at_upstream'],
            'physical_source_renormalized': False,
        },
        'missing_primitive': MISSING_REMAINDER,
        'next_primitive': NEXT_PRIMITIVE,
        'physical_local_gate': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'field_or_source_evolutions': 0,
        'manufactured_fixture_used_as_production_constant': False,
        'background_tail_certificate_transferred': False,
        'v1_pilot_rewritten': False,
        'v2_record_rewritten': False,
        'v3_record_rewritten': False,
        'invariance_record_rewritten': False,
        'sample_used_as_proof': False,
        'fitted_decay': False,
    })


def validate_cubic_uv_report(report, root=None):
    """Reject invented C4/C_M, independent R=0, closed gate, or samples as proof."""
    report = _as_map(report, 'cubic UV report')
    if report.get('schema') != SCHEMA:
        raise ValueError('unexpected cubic UV schema')
    if dict(_as_map(report.get('reduced_transport'), 'reduced transport')) != dict(
            reduced_massless_transport_identities()):
        raise ValueError('reduced transport identities or scope changed')
    if report.get('background_tail_certificate_transferred'):
        raise ValueError('background tail certificate must not transfer')
    if report.get('fitted_decay') or report.get('sample_used_as_proof'):
        raise ValueError('fitted decay or samples may not prove the cubic coefficient')
    if report.get('manufactured_fixture_used_as_production_constant'):
        raise ValueError('manufactured fixture cannot supply a production constant')
    if report.get('field_or_source_evolutions') not in (0, None):
        raise ValueError('cubic control may not launch field or source evolutions')
    if report.get('physical_EXISTENCE_certificate') or report.get(
            'physical_NONEXISTENCE_certificate'):
        raise ValueError('OPEN calculation may not be reported as PASS or NON_EXISTENCE')
    if report.get('physical_local_gate') != 'OPEN':
        raise ValueError('local gate remains OPEN')
    if report.get('v1_pilot_rewritten') or report.get('v2_record_rewritten') or report.get(
            'v3_record_rewritten') or report.get('invariance_record_rewritten'):
        raise ValueError('historical UV records are immutable')
    if report.get('R_zeroed_while_freezing_h3_n4'):
        raise ValueError('R may not be zeroed independently of h3 and n4')
    if report.get('dummy_symbol_L0_A_j'):
        raise ValueError('production cubic coefficient may not use dummy L0 A_j')
    if report.get('A3_locally_algebraic_in_metric_jets'):
        raise ValueError('minor C is not a local metric jet')
    if report.get('massless_Im_Ts_major_A2_sets_upstream_datum_to_zero'):
        raise ValueError('massless Im(T_s major A2) does not own the upstream datum')
    if not report.get('generated_massless_A3_n4_constructed'):
        raise ValueError('generated massless A3/n4 were not constructed')
    if report.get('production_A4_constructed'):
        raise ValueError('A4 is not constructed here')
    if report.get('validated_cubic_coefficient'):
        raise ValueError('diagnostic quadrature is not a remainder-validated C4')
    _require_missing(report.get('numerical_C4'), 'C4')
    _require_missing(report.get('numerical_C_M'), 'C_M')
    _require_missing(report.get('vacuum_integrated_tail_N_beta'), 'vacuum tail')
    _require_missing(report.get('certified_quadrature_error_bound'), 'quadrature remainder')
    identities = _as_map(report.get('identities'), 'identities')
    if identities.get('R_zeroed_while_freezing_h3_n4'):
        raise ValueError('identities may not zero R independently')
    if identities.get('invariance_uses_n2_up_equals_zero'):
        raise ValueError('k is representative data, not n2_up=0')
    diagnostic = _as_map(report.get('diagnostic_quadrature'), 'diagnostic quadrature')
    if diagnostic.get('validated') or diagnostic.get('remainder_proof_supplied'):
        raise ValueError('diagnostic cubic quadrature remains unvalidated')
    if diagnostic.get('certified_quadrature_error_bound') is not None:
        raise ValueError('no certified quadrature remainder is owned')
    if diagnostic.get('R_zeroed_independently'):
        raise ValueError('diagnostic representative zeroed R independently')
    for name in ('I_max_abs_N', 'I_max_abs_beta', 'refinement_discrepancy_N',
                 'refinement_discrepancy_beta'):
        value = diagnostic.get(name)
        if isinstance(value, (bool, np.bool_)) or value is None or not np.isfinite(value):
            raise ValueError('finite diagnostic ' + name + ' required')
    if max(abs(float(diagnostic['I_max_abs_N'])), abs(float(diagnostic['I_max_abs_beta']))) == 0:
        raise ArithmeticError('live 0b0e4ced cubic difference may not vanish')
    if report.get('profile_identity') != CURRENT_HISTORY_IDENTITY:
        raise ValueError('cubic report is not bound to history 0b0e4ced')
    if report.get('source_columns_changed'):
        raise ValueError('source columns must stay unchanged')
    signs = _as_map(report.get('action_signs'), 'action signs')
    if signs.get('old_v3_beta_shorthand_is_authority'):
        raise ValueError('old v3 beta shorthand is not authority')
    if signs.get('beta') != '+mu/(2 pi) I current':
        raise ValueError('beta must be plus the identity current')
    pair = _as_map(report.get('pair'), 'pair')
    if float(pair.get('mass', 1.0)) != 0.0:
        raise ValueError('this owner is the massless original pair')
    zero = _as_map(report.get('zero_history_control'), 'zero-history control')
    if max(abs(float(zero['I_max_abs_N'])), abs(float(zero['I_max_abs_beta']))) > 1e-10:
        raise ArithmeticError('zero-history cubic difference must vanish')
    return True
