"""Current-history characteristic transports for major A2, A3 and n4.

Production owner of the occupied T_s right-hand sides on the switched state
law C_Sigma[g]=F[g] C_src F[g]^dagger with dC_src=0 and the fixed paired
horizon preparation. Coefficients are history-minus-reference. The massless
vanishing of the local Im(T_s major A2) right-hand side does not set the
upstream datum or the full-envelope transported coefficient to zero.

Constructed here
----------------
Explicit T_s formulae for major A2 (real and imaginary), the massless
T_s formulae for major A3, the n4 characteristic source, and the production
definition of C4 as the history-minus-reference e^{-3} N,beta contraction.
Interval majorants for the incoming phase q and for the history-minus-reference
real major A2. The massless history-minus-reference imaginary major A2 is
exactly 0 because both histories share the unowned upstream datum and the
transport right-hand side vanishes.

Not constructed here
--------------------
Pointwise production values of major A3 and n4, a numerical C4, L0 A4 H2
integrals, the same-column upstream H2 remainder, C_M, or the vacuum tail.
Diagnostic manufactured H2 fixtures are not production constants. Fitted
decay and collocation samples are not proofs.
"""
from collections.abc import Mapping
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from types import MappingProxyType
import json
import math

import numpy as np
import sympy as sp

from .nsc_dirac_source_phase import SOURCE_SIGNS
from .nsc_ks_current_history_bounds import radius_bounds
from .nsc_ks_finite_history_uv_coefficients import (
    ARCHIVE_RHO_UP,
    CURRENT_HISTORY_IDENTITY_PREFIX,
    _apply_L0,
    _major_index,
    _minor_index,
    _minor_A1,
    _projector,
    _real_residual,
    _sym_pauli,
    _sym_potential,
    original_angular_pair_ledger,
    vacuum_upstream_conditions,
)
from .nsc_ks_profile_identity import profile_identity
from .nsc_ks_residual_error import sqrt_upper
from .nsc_local_incoming_family import LocalIncomingFamily


HISTORY_RECORD = 'results/development/nsc-ks-gate-history-lm-broyden.json'
CAUCHY_RECORD = 'results/development/nsc-mode-resolved-cauchy-state.json'
CURRENT_HISTORY_IDENTITY = (
    '0b0e4cedfdb695d342c7f7347740631e916bdaee3a977df56808c725b25dacb0')
TRANSPORT_SCHEMA = 'NSC-KS-CURRENT-UV-TRANSPORT-v1'
NEXT_PRIMITIVE = (
    'interval/analytic bound of minor A3 = a^2/(2 i) Pi_{-s} L0 A2 on the '
    'current-history slab, which needs second geometry jets of r and of the '
    'transported phase q together with Re(major A2); then T_s major A3 and '
    'T_s n4; then evaluation of the history-minus-reference e^{-3} N,beta '
    'contraction on I. After that, production L0 A4 H2 integrals and the '
    'same-column upstream H2 remainder for a finite C_M and vacuum tail')
MISSING_MINOR_A3 = (
    'minor A3 = a^2/(2 i) Pi_{-s} L0 A2 on the current-history slab, which '
    'retains second jets of r, q and Re(major A2); not a local metric jet')
C4_SCOPE = (
    'history-minus-reference e^{-3} N,beta contraction on I of the original '
    'equal-mu angular pairs, occupied vacuum envelope, action minus once, '
    'measure de/(2 pi)')
FULL_ENVELOPE_SCOPE = (
    'full-envelope e^{-3} N,beta coefficient of g, including the reference '
    'affine column')


def _root():
    return Path(__file__).resolve().parents[2]


def _require_positive(value, name):
    if isinstance(value, (bool, np.bool_)) or value is None:
        raise ValueError('explicit positive ' + name + ' required')
    number = float(value)
    if not np.isfinite(number) or number <= 0:
        raise ValueError('positive finite ' + name + ' required')
    return number


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


def _binary64_upper(value):
    if isinstance(value, Mapping) and 'binary64_upper' in value:
        value = value['binary64_upper']
    number = Fraction(value)
    rounded = float(number)
    if Fraction(rounded) < number:
        rounded = math.nextafter(rounded, math.inf)
    return rounded


def _int_digest(value):
    number = int(value)
    width = max(1, (number.bit_length() + 8) // 8)
    return sha256(number.to_bytes(width, 'big', signed=True)).hexdigest()


def bound_record(value):
    """Outward binary64 upper bound plus a digest of the exact rational.

    Chebyshev-derived numerators can exceed the default decimal-digit limit,
    so the exact fraction is bound by hash rather than a decimal string.
    """
    number = Fraction(value)
    upper = _binary64_upper(number)
    return {
        'binary64_upper': upper,
        'ieee754_hex': float(upper).hex(),
        'numerator_sha256': _int_digest(number.numerator),
        'denominator_sha256': _int_digest(number.denominator),
        'sign': int(number.numerator >= 0),
    }


def _majorant_number(value, name):
    if isinstance(value, Mapping):
        if 'binary64_upper' in value:
            value = value['binary64_upper']
        elif 'exact_rational' in value:
            value = Fraction(value['exact_rational'])
        else:
            raise ValueError('proved majorant required: ' + name)
    if isinstance(value, (bool, np.bool_)) or value is None:
        raise ValueError('proved positive majorant required: ' + name)
    number = Fraction(value)
    if number <= 0:
        raise ValueError('proved majorant must stay strictly positive: ' + name)
    return number


def axial_rho_abs_upper():
    """|a_rho| on [1, 33/32] from a>=4/5 and min rho*angle = pi/4 at rho=1.

    d(rho theta)/d rho = theta - rho/(1+rho^2) is positive at rho=1, so
    rho*angle is minimized at the left endpoint. pi>3.14 gives theta>157/200.
    Then |a_rho|=3(1-rho theta)/a <= 3(1-157/200)/(4/5)=129/160.
    """
    return Fraction(129, 160)


def axial_abs_upper():
    """a(1)^2 = 3*pi/2-4 < 5/7 using pi<22/7; a is decreasing on the slab."""
    return sqrt_upper(Fraction(5, 7))


def reference_radius_abs_upper():
    """sqrt(1+rho^2) <= sqrt(1+(33/32)^2) on the owned background slab."""
    return sqrt_upper(Fraction(2113, 1024))


@lru_cache(maxsize=1)
def characteristic_transport_identities():
    """Owned T_s right-hand sides on the actual L0 recurrence, both signs.

    Dummy L0 A_j symbols are not used. Residues are exact symbolic zeros.
    They are not a numerical C4 or C_M.
    """
    rho, z = sp.symbols('rho z', real=True)
    mass, angular = sp.symbols('m ell', real=True)
    axial = sp.Function('a', real=True)(rho)
    radius = sp.Function('r', real=True)(rho, z)
    phase = sp.Function('q', real=True)(rho, z)
    height = sp.Function('h', real=True)(rho, z)
    real_major = sp.Function('R', real=True)(rho, z)
    s1, s2, s3 = _sym_pauli()
    potential = _sym_potential(mass, angular, radius)
    square = mass**2 + angular**2 / radius**2
    identities = {}
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        major = _major_index(source_sign)
        minor = _minor_index(source_sign)
        occupied, vacant = _projector(source_sign, s3)
        a1 = sp.zeros(2, 1)
        a1[major] = sp.I * phase
        a1[minor] = _minor_A1(axial, mass, angular, radius, source_sign)
        a1_rho = source_sign / axial**2 * a1[major].diff(z) - sp.I / 2 * square
        l0_a1 = _apply_L0(a1, axial, potential, s3, rho, z).subs(
            a1[major].diff(rho), a1_rho)
        minor_a2 = axial**2 / (2 * sp.I) * (vacant * l0_a1)[minor]
        rhs_a2 = (sp.I / axial) * potential[major, minor] * minor_a2
        imag_rhs = sp.im(sp.expand_complex(rhs_a2))
        real_rhs = sp.re(sp.expand_complex(rhs_a2))
        parent_imag = -mass * angular / 4 * (
            source_sign * axial**2 * radius.diff(rho) + radius.diff(z)
        ) / radius**2
        parent_real = (
            -axial.diff(rho) * mass**2 * axial * radius**3
            + angular**2 * (
                -axial.diff(rho) * axial * radius
                + radius.diff(rho) * axial**2
                + source_sign * radius.diff(z)
                + 2 * phase * radius)
            + 2 * mass**2 * phase * radius**3
        ) / (4 * radius**3)
        identities['%s/Im_Ts_major_A2_parent' % tag] = _real_residual(
            imag_rhs - parent_imag)
        identities['%s/Re_Ts_major_A2_parent' % tag] = _real_residual(
            real_rhs - parent_real)
        identities['%s/Re_Ts_major_A2_depends_on_q' % tag] = (
            'depends' if phase in real_rhs.atoms(sp.Function) else 'independent')
        identities['%s/Im_Ts_major_A2_independent_of_q' % tag] = (
            'independent' if phase not in parent_imag.atoms(sp.Function)
            else 'depends')
        identities['%s/Im_Ts_major_A2_vanishes_for_m_0' % tag] = _real_residual(
            imag_rhs.subs(mass, 0))
        identities['%s/Im_Ts_major_A2_odd_in_ell' % tag] = _real_residual(
            imag_rhs + imag_rhs.subs(angular, -angular))
        a2 = sp.zeros(2, 1)
        a2[major] = real_major + sp.I * height
        a2[minor] = minor_a2
        a2_rho = source_sign / axial**2 * a2[major].diff(z) + rhs_a2
        l0_a2 = _apply_L0(a2, axial, potential, s3, rho, z).subs(
            a2[major].diff(rho), a2_rho)
        minor_a3 = axial**2 / (2 * sp.I) * (vacant * l0_a2)[minor]
        identities['%s/minor_A3_depends_on_major_A2' % tag] = (
            'depends' if real_major in minor_a3.atoms(sp.Function)
            or height in minor_a3.atoms(sp.Function) else 'independent')
        rhs_a3 = (sp.I / axial) * potential[major, minor] * minor_a3
        real_a3 = sp.re(sp.expand_complex(rhs_a3.subs(mass, 0)))
        parent_real_a3 = angular**2 * height / (2 * radius**2)
        identities['%s/Re_Ts_major_A3_massless_is_ell2_h_over_2r2' % tag] = (
            _real_residual(real_a3 - parent_real_a3))
        identities['%s/Re_Ts_major_A3_massless_depends_on_h' % tag] = (
            'depends' if height in real_a3.atoms(sp.Function) else 'independent')
        identities['%s/Im_Ts_major_A3_massless_depends_on_R' % tag] = (
            'depends' if real_major in sp.im(sp.expand_complex(
                rhs_a3.subs(mass, 0))).atoms(sp.Function) else 'independent')
        identities['%s/minor_A2_massless_independent_of_h' % tag] = (
            'independent' if height not in minor_a2.subs(mass, 0).atoms(
                sp.Function) else 'depends')
    residual_keys = [
        key for key, value in identities.items()
        if not key.endswith((
            'depends_on_q', 'independent_of_q', 'depends_on_major_A2',
            'depends_on_h', 'depends_on_R', 'independent_of_h'))]
    leaves = []
    for key in residual_keys:
        value = identities[key]
        leaves.extend(value if isinstance(value, list) else [value])
    if set(leaves) != {'0'}:
        raise ArithmeticError(identities)
    for key, value in identities.items():
        if key.endswith('depends_on_q') and value != 'depends':
            raise ArithmeticError('Re(T_s major A2) must retain the transported phase')
        if key.endswith('independent_of_q') and value != 'independent':
            raise ArithmeticError('Im(T_s major A2) must stay independent of q')
        if key.endswith('depends_on_major_A2') and value != 'depends':
            raise ArithmeticError('minor A3 must retain transported major A2')
        if key.endswith('depends_on_h') and value != 'depends':
            raise ArithmeticError('massless Re(T_s major A3) must retain h')
        if key.endswith('depends_on_R') and value != 'depends':
            raise ArithmeticError('massless Im(T_s major A3) must retain R')
        if key.endswith('independent_of_h') and value != 'independent':
            raise ArithmeticError('massless minor A2 must stay independent of h')
    return MappingProxyType({
        'residuals': identities,
        'T_s': 'T_s = d_rho - s/a^2 d_z',
        'Im_Ts_major_A2': '-m ell/4 (s a^2 r_rho + r_z)/r^2',
        'Re_Ts_major_A2': (
            '[-a_rho m^2 a r^3 + ell^2 (-a_rho a r + r_rho a^2 + s r_z + 2 q r) '
            '+ 2 m^2 q r^3] / (4 r^3)'),
        'T_s_major_A1': '-i/2 (m^2 + ell^2/r^2), so T_s q = -1/2 (m^2 + ell^2/r^2)',
        'Re_Ts_major_A3_massless': 'ell^2 h / (2 r^2)',
        'n4_transport': (
            'T_s n4 = -(s/a^2) d_z (2 |minor A2|^2 + 4 Re(conj(minor A1) minor A3)); '
            'upstream n4=0 does not remove minor A3'),
        'dummy_symbol_L0_A_j': False,
        'production_A2_transport_constructed': True,
        'production_A3_A4_constructed': False,
        'major_A2_replaced_by_local_metric_jet': False,
        'A3_locally_algebraic_in_metric_jets': False,
        'massless_Im_Ts_major_A2_sets_upstream_datum_to_zero': False,
        'massless_Im_Ts_major_A2_sets_full_envelope_coefficient_to_zero': False,
    })


@lru_cache(maxsize=1)
def massless_j3_drops_upstream_imag_major_A2():
    """Occupied identity current J_3 at m=0 does not see the unowned h datum.

    After real-jet substitution the coefficient of h in Im(A0^dag A3_z +
    A1^dag A2_z + A2^dag A1_z) - s n4 vanishes. The full-envelope value of
    h is still not set to zero.
    """
    rho, z = sp.symbols('rho z', real=True)
    angular = sp.symbols('ell', real=True)
    axial = sp.Function('a', real=True)(rho)
    radius = sp.Function('r', real=True)(rho, z)
    phase = sp.Function('q', real=True)(rho, z)
    height = sp.Function('h', real=True)(rho, z)
    real_major = sp.Function('R', real=True)(rho, z)
    height3 = sp.Function('h3', real=True)(rho, z)
    real3 = sp.Function('R3', real=True)(rho, z)
    n4 = sp.Function('n4', real=True)(rho, z)
    a_rho, r_rho, r_z = sp.symbols('a_rho r_rho r_z', real=True)
    q_z, r_z_major, h_z, r3_z, h3_z = sp.symbols(
        'q_z R_z h_z R3_z h3_z', real=True)
    identities = {}
    real_subs = {
        axial.diff(rho): a_rho,
        radius.diff(rho): r_rho,
        radius.diff(z): r_z,
        phase.diff(z): q_z,
        real_major.diff(z): r_z_major,
        height.diff(z): h_z,
        real3.diff(z): r3_z,
        height3.diff(z): h3_z,
    }
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        major = _major_index(source_sign)
        minor = _minor_index(source_sign)
        minor_a1 = axial / 2 * (-sp.I * source_sign * angular / radius)
        minor_a2 = (
            source_sign * angular * axial * (
                -a_rho * axial * radius + r_rho * axial**2
                + source_sign * r_z + 2 * phase * radius)
            / (4 * radius**2))
        a0 = sp.zeros(2, 1)
        a0[major] = 1
        a1 = sp.zeros(2, 1)
        a1[major] = sp.I * phase
        a1[minor] = minor_a1
        a2 = sp.zeros(2, 1)
        a2[major] = real_major + sp.I * height
        a2[minor] = minor_a2
        a3 = sp.zeros(2, 1)
        a3[major] = real3 + sp.I * height3
        a3[minor] = sp.Function('C')(rho, z)
        interference = (a0.H * a3.diff(z) + a1.H * a2.diff(z) + a2.H * a1.diff(z))[0]
        current = sp.im(sp.expand_complex(interference)) - source_sign * n4
        current = sp.simplify(current.xreplace(real_subs))
        # The unowned constant h may not appear in the identity current.
        identities['%s/J3_independent_of_h' % tag] = (
            'independent' if height not in current.atoms(sp.Function)
            else 'depends')
        identities['%s/J3_coefficient_of_h' % tag] = _real_residual(
            sp.diff(current, height) if height in current.atoms(sp.Function)
            else sp.Integer(0))
        # Control: the raw interference before Im() still contains h.
        identities['%s/raw_A2dag_A1z_contains_h' % tag] = (
            'depends' if height in (a2.H * a1.diff(z))[0].atoms(sp.Function)
            else 'independent')
    if any(value != 'independent' for key, value in identities.items()
           if key.endswith('J3_independent_of_h')):
        raise ArithmeticError('massless J_3 must drop the unowned imag major A2')
    if any(value != '0' for key, value in identities.items()
           if key.endswith('J3_coefficient_of_h')):
        raise ArithmeticError('massless J_3 coefficient of h must vanish')
    if any(value != 'depends' for key, value in identities.items()
           if key.endswith('raw_A2dag_A1z_contains_h')):
        raise ArithmeticError('control must detect h in the raw A2^dag A1_z term')
    return MappingProxyType({
        'residuals': identities,
        'massless_J3_independent_of_upstream_imag_major_A2': True,
        'massless_Im_Ts_major_A2_sets_upstream_datum_to_zero': False,
        'full_envelope_imag_major_A2_set_to_zero': False,
        'trV_rho3_checked_separately': True,
    })


def history_minus_reference_imag_major_A2(mass):
    """Exact difference of Im(major A2) after shared upstream data.

    For m=0 the local T_s right-hand side vanishes on both histories. They
    share the unowned upstream datum, so the difference is identically 0.
    The full-envelope coefficient and the upstream datum remain unowned.
    """
    if isinstance(mass, (bool, np.bool_)) or mass is None:
        raise ValueError('explicit real mass required')
    mass = float(mass)
    if not np.isfinite(mass) or mass < 0:
        raise ValueError('nonnegative finite mass required')
    if mass == 0.0:
        return MappingProxyType({
            'numerical_value': 0.0,
            'proved_exact': True,
            'scope': 'history-minus-reference',
            'full_envelope_value': None,
            'upstream_imag_major_A2_owned': False,
            'transport_rhs_this_massless_pair': 0.0,
            'status': 'CLOSED',
        })
    return MappingProxyType({
        'numerical_value': None,
        'proved_exact': False,
        'scope': 'history-minus-reference',
        'full_envelope_value': None,
        'upstream_imag_major_A2_owned': False,
        'transport_rhs_this_massless_pair': None,
        'status': 'OPEN',
        'missing_primitive': (
            'characteristic integral of Im(T_s major A2) = '
            '-m ell/4 (s a^2 r_rho + r_z)/r^2 with shared upstream datum'),
    })


def current_history_geometry_majorants(family, rho_up=ARCHIVE_RHO_UP):
    """Fraction majorants of chart and radius jets on the declared slab."""
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError('owned LocalIncomingFamily required')
    bounds = radius_bounds(family, rho_up)
    duration = bounds['rho_up'] - 1
    if duration <= 0:
        raise ArithmeticError('positive characteristic duration required')
    width = Fraction(family.normal_outer) - Fraction(family.normal_inner)
    if width <= 0:
        raise ValueError('positive normal-window width required')
    delta = bounds['delta_radius_bounds']
    w_bounds = bounds['w']['profile_bounds']
    u_bounds = bounds['U']['profile_bounds']
    chi_s = 8 / width
    delta_s = chi_s * delta[0] + w_bounds[0] + bounds['normal_support']**2 * u_bounds[0] / 2
    a_rho = axial_rho_abs_upper()
    a_upper = axial_abs_upper()
    r_ref_max = reference_radius_abs_upper()
    r_max = r_ref_max + delta[0]
    r_ref_rho = bounds['rho_up'] / bounds['reference_radius_lower']
    r_rho = r_ref_rho + delta_s / bounds['axial_lower']
    return MappingProxyType({
        'duration': duration,
        'axial_lower': bounds['axial_lower'],
        'axial_upper': a_upper,
        'axial_rho_abs_upper': a_rho,
        'radius_lower': bounds['radius_lower'],
        'radius_upper': r_max,
        'reference_radius_lower': bounds['reference_radius_lower'],
        'reference_radius_upper': r_ref_max,
        'reference_r_rho_abs_upper': r_ref_rho,
        'r_rho_abs_upper': r_rho,
        'r_z_abs_upper': delta[1],
        'delta_radius_abs_upper': delta[0],
        'delta_s_abs_upper': delta_s,
        'normal_support': bounds['normal_support'],
        'normal_window_width': width,
        'rho_up': bounds['rho_up'],
        'profile_evaluation_roundoff_included': False,
    })


def phase_and_real_major_A2_difference_majorants(
        family, mass, angular, rho_up=ARCHIVE_RHO_UP):
    """Interval majorants of q and of history-minus-reference Re(major A2).

    |T_s q| = (m^2 + ell^2/r^2)/2. History-minus-reference uses
    |r_g^{-2}-r_ref^{-2}| <= |delta r|(2 r_ref + |delta r|)/(r_g^2 r_ref^2).
    Re(T_s major A2) is bounded by the closed parent formula. The triangle
    inequality |delta R| <= duration*(|F_g|+|F_ref|) is a majorant, not a
    sample and not a fitted decay.
    """
    mass = float(mass)
    angular = Fraction(angular) ** 2
    if mass < 0 or angular < 0:
        raise ValueError('nonnegative mass and angular^2 required')
    geom = current_history_geometry_majorants(family, rho_up)
    duration = geom['duration']
    r_min = geom['radius_lower']
    r_max = geom['radius_upper']
    r_ref_min = geom['reference_radius_lower']
    r_ref_max = geom['reference_radius_upper']
    mass2 = Fraction(mass) ** 2 if mass else Fraction(0)
    ell2 = angular
    q_full = duration * (mass2 + ell2 / r_min**2) / 2
    q_ref = duration * (mass2 + ell2 / r_ref_min**2) / 2
    dinv = geom['delta_radius_abs_upper'] * (
        2 * r_ref_max + geom['delta_radius_abs_upper']
    ) / (r_min**2 * r_ref_min**2)
    delta_q = (ell2 / 2) * duration * dinv
    a_upper = geom['axial_upper']
    a_rho = geom['axial_rho_abs_upper']

    def real_rhs(radius_lo, radius_hi, r_rho, r_z, phase):
        numerator = (
            a_rho * mass2 * a_upper * radius_hi**3
            + ell2 * (
                a_rho * a_upper * radius_hi
                + r_rho * a_upper**2
                + r_z
                + 2 * phase * radius_hi)
            + 2 * mass2 * phase * radius_hi**3)
        return numerator / (4 * radius_lo**3)

    f_g = real_rhs(
        r_min, r_max, geom['r_rho_abs_upper'], geom['r_z_abs_upper'], q_full)
    f_ref = real_rhs(
        r_ref_min, r_ref_max, geom['reference_r_rho_abs_upper'], 0, q_ref)
    delta_r = duration * (f_g + f_ref)
    if min(q_full, q_ref, delta_q, f_g, f_ref, delta_r) < 0:
        raise ArithmeticError('transport majorant became negative')
    return MappingProxyType({
        'mass': mass,
        'angular_squared': ell2,
        'q_abs_upper': q_full,
        'q_reference_abs_upper': q_ref,
        'delta_q_abs_upper': delta_q,
        'Re_Ts_major_A2_g_abs_upper': f_g,
        'Re_Ts_major_A2_ref_abs_upper': f_ref,
        'history_minus_reference_real_major_A2_abs_upper': delta_r,
        'fitted_decay': False,
        'sample_used_as_proof': False,
        'full_envelope_real_major_A2': None,
        'geometry': geom,
    })


def history_minus_reference_c4_definition():
    """Production C4 is the history-minus-reference e^{-3} contraction on I."""
    first_formulas = {
        'J_3': (
            'Im(A0^dag A3_z + A1^dag A2_z + A2^dag A1_z) - s n4, with A0_z=0'),
        'action_N3': (
            'mu/(2 pi) times the occupied contraction of V at e^{-3} and '
            '(S3/a) J_3, action minus once'),
        'action_beta3': '-mu/(2 pi) times the occupied J_3; action minus once',
        'C4': (
            'equal-mu original pair of those e^{-3} action densities of '
            'g minus the same densities of the reference, on I'),
    }
    pairs = original_angular_pair_ledger()
    return MappingProxyType({
        'definition_scope': C4_SCOPE,
        'full_envelope_scope': FULL_ENVELOPE_SCOPE,
        'history_minus_reference_coefficient_formed': True,
        'full_envelope_coefficient_formed': False,
        'numerical_value': None,
        'formulas': first_formulas,
        'vertices': 'N: V and S3/a; beta: -I; action minus once; de/(2 pi)',
        'angular_pairs': len(pairs),
        'equal_mu_required': True,
        'same_cutoff_grid_required': True,
        'massless_J3_independent_of_upstream_imag_major_A2': True,
        'massless_Im_Ts_sets_C4_to_zero': False,
        'replaced_by_local_metric_jet': False,
        'fitted_decay': False,
        'missing_primitive': MISSING_MINOR_A3,
        'leading_e_minus2_cancellation_is_not_this_value': True,
        'e_minus3_first_allowed_order_is_not_a_tail_bound': True,
    })


def smallest_original_pair(root=None):
    """Smallest original |ell| pair: group 1, massless, |ell|=sqrt(5), mu=6."""
    root = Path(root) if root is not None else _root()
    pairs = original_angular_pair_ledger()
    cauchy = json.loads((root / CAUCHY_RECORD).read_text())
    channel = next(row for row in cauchy['channels'] if row['index'] == pairs[0]['group'])
    if channel['angular_eigenvalue'] <= 0:
        raise ValueError('smallest representative pair must be a nonzero-angular family')
    cutoff = 320.0 if channel['index'] in (10, 11, 12, 31, 32) else 160.0
    return MappingProxyType({
        'group': int(channel['index']),
        'mass': float(channel['compact_mass']),
        'absolute_angular': float(channel['angular_eigenvalue']),
        'multiplicity_per_signed_family': float(pairs[0]['multiplicity_per_signed_family']),
        'cutoff': cutoff,
        'archived_splits': (160.0, 320.0),
        'same_cutoff_grid': True,
        'compact_level': int(channel['compact_level']),
        'angular_level': int(channel['angular_level']),
    })


def current_history_characteristic_transport_report(root=None):
    """Executable transport report on the live 0b0e4ced history."""
    root = Path(root) if root is not None else _root()
    history = json.loads((root / HISTORY_RECORD).read_text())
    family = LocalIncomingFamily(np.array(history['history']['coefficients']))
    identity = profile_identity(family)
    if (not CURRENT_HISTORY_IDENTITY.startswith(CURRENT_HISTORY_IDENTITY_PREFIX)
            or identity != CURRENT_HISTORY_IDENTITY):
        raise ValueError('expected exact live history 0b0e4ced')
    pair = smallest_original_pair(root)
    identities = characteristic_transport_identities()
    j3 = massless_j3_drops_upstream_imag_major_A2()
    imag = history_minus_reference_imag_major_A2(pair['mass'])
    majorants = phase_and_real_major_A2_difference_majorants(
        family, pair['mass'], pair['absolute_angular'])
    c4 = history_minus_reference_c4_definition()
    upstream = vacuum_upstream_conditions()
    if upstream['major_A1'] != 0:
        raise ArithmeticError('owned upstream major A1 must vanish')
    return MappingProxyType({
        'schema': TRANSPORT_SCHEMA,
        'profile_identity': identity,
        'pair': dict(pair),
        'identities': identities,
        'j3': j3,
        'imag_major_A2': imag,
        'majorants': MappingProxyType({
            'mass': float(majorants['mass']),
            'angular_squared': bound_record(majorants['angular_squared']),
            'q_abs_upper': bound_record(majorants['q_abs_upper']),
            'q_reference_abs_upper': bound_record(majorants['q_reference_abs_upper']),
            'delta_q_abs_upper': bound_record(majorants['delta_q_abs_upper']),
            'Re_Ts_major_A2_g_abs_upper': bound_record(
                majorants['Re_Ts_major_A2_g_abs_upper']),
            'Re_Ts_major_A2_ref_abs_upper': bound_record(
                majorants['Re_Ts_major_A2_ref_abs_upper']),
            'history_minus_reference_real_major_A2_abs_upper': bound_record(
                majorants['history_minus_reference_real_major_A2_abs_upper']),
            'fitted_decay': False,
            'sample_used_as_proof': False,
            'full_envelope_real_major_A2': None,
            'duration_binary64_upper': _binary64_upper(
                majorants['geometry']['duration']),
            'radius_lower_binary64_upper': _binary64_upper(
                majorants['geometry']['radius_lower']),
        }),
        'c4': c4,
        'upstream': {
            'A0': upstream['A0'],
            'major_A1': upstream['major_A1'],
            'history_flat_at_upstream': upstream['history_flat_at_upstream'],
            'imag_major_A2_owned': False,
            'physical_source_renormalized': False,
        },
        'characteristic_transport_equations_owned': True,
        'characteristic_integrals_owned': {
            'history_minus_reference_imag_major_A2_massless': True,
            'real_major_A2_difference_majorant': True,
            'phase_q_majorant': True,
            'major_A3': False,
            'n4': False,
        },
        'production_A3_A4_constructed': False,
        'next_primitive': NEXT_PRIMITIVE,
        'numerical_C4': None,
        'numerical_C_M': None,
        'vacuum_integrated_tail_N_beta': None,
        'background_tail_certificate_transferred': False,
        'fitted_decay': False,
        'sample_used_as_proof': False,
        'manufactured_fixture_used_as_production_constant': False,
        'field_or_source_evolutions': 0,
    })


def validate_characteristic_transports(report, root=None):
    """Reject invented A2/A3/n4/C4, null-to-zero, wrong C4 scope, or samples as proof."""
    report = _as_map(report, 'transport report')
    if report.get('schema') != TRANSPORT_SCHEMA:
        raise ValueError('unexpected current-history UV transport schema')
    if report.get('background_tail_certificate_transferred'):
        raise ValueError('background tail certificate must not transfer')
    if report.get('fitted_decay') or report.get('sample_used_as_proof'):
        raise ValueError('fitted decay or samples may not prove a transport bound')
    if report.get('manufactured_fixture_used_as_production_constant'):
        raise ValueError('manufactured fixture cannot supply a production constant')
    if report.get('field_or_source_evolutions') not in (0, None):
        raise ValueError('transport control may not launch field or source evolutions')
    identities = _as_map(report.get('identities'), 'identities')
    if identities.get('dummy_symbol_L0_A_j'):
        raise ValueError('production transport may not use dummy L0 A_j symbols')
    if identities.get('massless_Im_Ts_major_A2_sets_upstream_datum_to_zero'):
        raise ValueError('massless Im(T_s major A2) does not own the upstream datum')
    if identities.get('massless_Im_Ts_major_A2_sets_full_envelope_coefficient_to_zero'):
        raise ValueError('massless Im(T_s major A2) does not zero the full coefficient')
    if not identities.get('production_A2_transport_constructed'):
        raise ValueError('major A2 transport was not constructed')
    if identities.get('production_A3_A4_constructed'):
        raise ValueError('production A3/A4 are not constructed here')
    imag = _as_map(report.get('imag_major_A2'), 'imag major A2')
    pair = _as_map(report.get('pair'), 'pair')
    if float(pair.get('mass', 1.0)) == 0.0:
        if imag.get('numerical_value') != 0.0 or not imag.get('proved_exact'):
            raise ValueError('massless history-minus-reference Im(major A2) must be exact 0')
        if imag.get('scope') != 'history-minus-reference':
            raise ValueError('proved zero has the wrong coefficient scope')
        if imag.get('full_envelope_value') is not None:
            raise ValueError('full-envelope imag major A2 may not be invented')
        if imag.get('upstream_imag_major_A2_owned'):
            raise ValueError('upstream imag major A2 is not owned')
        if imag.get('transport_rhs_this_massless_pair') not in (0.0, 0):
            raise ValueError('massless Im(T_s major A2) identity was not preserved')
    _require_missing(imag.get('full_envelope_value'), 'full-envelope imag major A2')
    majorants = _as_map(report.get('majorants'), 'majorants')
    for name in (
            'q_abs_upper', 'delta_q_abs_upper',
            'history_minus_reference_real_major_A2_abs_upper'):
        _majorant_number(majorants.get(name), name)
    _require_missing(
        majorants.get('full_envelope_real_major_A2'), 'full-envelope real major A2')
    c4 = _as_map(report.get('c4'), 'C4')
    if not c4.get('history_minus_reference_coefficient_formed'):
        raise ValueError('production history-minus-reference C4 definition was not formed')
    if c4.get('full_envelope_coefficient_formed'):
        raise ValueError('C4 may not use the generic full-envelope coefficient')
    if c4.get('definition_scope') != C4_SCOPE:
        raise ValueError('C4 has the wrong coefficient scope')
    _require_missing(c4.get('numerical_value'), 'C4')
    _require_missing(report.get('numerical_C4'), 'C4')
    _require_missing(report.get('numerical_C_M'), 'C_M')
    _require_missing(report.get('vacuum_integrated_tail_N_beta'), 'vacuum tail')
    if report.get('profile_identity') != CURRENT_HISTORY_IDENTITY:
        raise ValueError('transport report is not bound to history 0b0e4ced')
    return True


def transport_arrays(report):
    """Deterministic arrays for the v3 payload; not a sample proof."""
    report = _as_map(report, 'transport report')
    validate_characteristic_transports(report)
    majorants = report['majorants']
    imag = report['imag_major_A2']
    pair = report['pair']
    pairs = original_angular_pair_ledger()
    return {
        'history_minus_reference_imag_major_A2_massless': np.array(
            [float(imag['numerical_value'])], dtype=np.float64),
        'q_abs_upper': np.array(
            [_binary64_upper(majorants['q_abs_upper'])], dtype=np.float64),
        'delta_q_abs_upper': np.array(
            [_binary64_upper(majorants['delta_q_abs_upper'])], dtype=np.float64),
        'history_minus_reference_real_major_A2_abs_upper': np.array(
            [_binary64_upper(
                majorants['history_minus_reference_real_major_A2_abs_upper'])],
            dtype=np.float64),
        'pair_mass': np.array([float(pair['mass'])], dtype=np.float64),
        'pair_absolute_angular': np.array(
            [float(pair['absolute_angular'])], dtype=np.float64),
        'pair_multiplicity': np.array(
            [float(pair['multiplicity_per_signed_family'])], dtype=np.float64),
        'angular_pair_count': np.array([len(pairs)], dtype=np.int64),
        'source_signs': np.array(SOURCE_SIGNS, dtype=np.int64),
        'numerical_C4': np.array([np.nan], dtype=np.float64),
        'numerical_C_M': np.array([np.nan], dtype=np.float64),
    }
