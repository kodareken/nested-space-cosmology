"""Interval bound of minor A3 = a^2/(2 i) Pi_{-s} L0 A2 on the current slab.

Stage-2 UV primitive on the switched state law, live profile 0b0e4ced and the
smallest original angular pair. The recurrence uses the actual Dirac operator
L0 and the transported major A2. Dummy L0 A_j symbols are not used.

Constructed here
----------------
The parent formula

    minor A3 = a^2/(2 i) T_{-s}(minor A2) - (a/2) V_{-s,s} (R + i h),

with V_{-s,s} = -m + i s ell/r, second geometry jets of r, first and second
jets of the transported phase q, interval majorants of minor A2, and an
interval majorant of the jet piece a^2/(2 i) T_{-s}(minor A2). History-minus-
reference uses the generated Re(major A2) increment already owned by the
characteristic transport. Leading paired E^{-2} cancellation is preserved
and is not this bound.

Not constructed here
--------------------
A complete |minor A3| or |delta minor A3| including the shared upstream
major A2 times delta(1/r), T_s major A3, T_s n4, numerical C4, C_M, or the
vacuum tail. The massless vanishing of Im(T_s major A2) does not set the
upstream datum to zero. Full-envelope Re(major A2) remains unowned.
"""
from collections.abc import Mapping
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
import json

import numpy as np
import sympy as sp

from .nsc_dirac_source_phase import SOURCE_SIGNS
from .nsc_ks_current_history_bounds import radius_bounds
from .nsc_ks_current_uv_transport import (
    ARCHIVE_RHO_UP,
    CURRENT_HISTORY_IDENTITY,
    HISTORY_RECORD,
    axial_abs_upper,
    axial_rho_abs_upper,
    bound_record,
    current_history_geometry_majorants,
    phase_and_real_major_A2_difference_majorants,
    smallest_original_pair,
)
from .nsc_ks_finite_history_uv_coefficients import (
    _apply_L0,
    _major_index,
    _minor_A1,
    _minor_index,
    _projector,
    _real_residual,
    _sym_pauli,
    _sym_potential,
    _T,
)
from .nsc_ks_finite_history_uv_remainder import FIRST_NONCANCELLING_PAIRED_ORDER
from .nsc_ks_profile_identity import profile_identity
from .nsc_ks_residual_error import nonnegative
from .nsc_local_incoming_family import LocalIncomingFamily


SCHEMA = 'NSC-KS-UV-MINOR-A3-v1'
MISSING_UPSTREAM_MAJOR_A2 = (
    'shared upstream major A2, the integration constant of T_s (R+i h) at '
    'rho_up, multiplied by delta(1/r) in the algebraic coupling '
    '-(a/2) V_{-s,s} (R+i h); not set to zero by massless Im(T_s major A2)=0')
NEXT_PRIMITIVE = (
    'a bound on the shared upstream major A2, or an identity that drops the '
    'A2_up delta(1/r) coupling in the e^{-3} contraction; then T_s major A3 '
    'and T_s n4; then history-minus-reference C4 on I. C_M still needs '
    'production L0 A4 H2 integrals and the same-column upstream H2 remainder')


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


def _require_explicit_mass_angular(mass, angular):
    if isinstance(mass, (bool, np.bool_)) or mass is None:
        raise ValueError('explicit real mass required')
    if isinstance(angular, (bool, np.bool_)) or angular is None:
        raise ValueError('explicit angular eigenvalue required')
    mass = Fraction(mass)
    angular = Fraction(angular)
    if mass < 0:
        raise ValueError('nonnegative mass required')
    return mass, angular


def axial_abs_squared_upper():
    """a(1)^2 = 3 pi/2-4 < 5/7; a is decreasing on the slab."""
    return Fraction(5, 7)


def axial_rho_rho_abs_upper():
    """|a_rho_rho| on [1, 33/32].

    a_rho = -3 u/a with u=1-rho theta. Differentiating gives
    a_rr = -3 u_rho/a - 9 u^2/a^3. The same u <= 43/200 as
    axial_rho_abs_upper. |u_rho| = theta - rho/(1+rho^2) is largest at
    rho=1: theta <= pi/4 < 11/14 and rho/(1+rho^2) is decreasing, so at
    least 1056/2113 on the slab. a >= 4/5.
    """
    u_max = Fraction(43, 200)
    u_rho = Fraction(11, 14) - Fraction(1056, 2113)
    if u_rho <= 0:
        raise ArithmeticError('u_rho majorant is not positive')
    a_min = Fraction(4, 5)
    return 3 * u_rho / a_min + 9 * u_max**2 / a_min**3


def plateau_s_derivative_abs_uppers(width):
    """Owned even C-infinity plateau: |chi'|<=8/width, |chi''|<=105/width^2.

    The constants are the same normalized-transition bounds as
    nsc_ks_radius_enclosure: |eta'|<=8, |eta''|<=105.
    """
    width = nonnegative(width, 'plateau width')
    if width <= 0:
        raise ValueError('positive plateau width required')
    return Fraction(8) / width, Fraction(105) / width**2


def _real_jets():
    rho, z = sp.symbols('rho z', real=True)
    axial = sp.Function('a', real=True)(rho)
    radius = sp.Function('r', real=True)(rho, z)
    phase = sp.Function('q', real=True)(rho, z)
    height = sp.Function('h', real=True)(rho, z)
    real_major = sp.Function('R', real=True)(rho, z)
    a_rho, a_rr, r_rho, r_z = sp.symbols('a_rho a_rr r_rho r_z', real=True)
    r_rr, r_rz, r_zz = sp.symbols('r_rr r_rz r_zz', real=True)
    q_rho, q_z = sp.symbols('q_rho q_z', real=True)
    return {
        'rho': rho, 'z': z, 'axial': axial, 'radius': radius, 'phase': phase,
        'height': height, 'real_major': real_major,
        'subs': {
            axial.diff(rho): a_rho,
            axial.diff(rho, 2): a_rr,
            radius.diff(rho): r_rho,
            radius.diff(z): r_z,
            radius.diff(rho, 2): r_rr,
            radius.diff(rho, z): r_rz,
            radius.diff(z, 2): r_zz,
            phase.diff(rho): q_rho,
            phase.diff(z): q_z,
        },
    }


def _parent_minor_a2(axial, radius, phase, mass, angular, source_sign, a_rho, r_rho, r_z):
    poly = (
        -a_rho * axial * radius
        + axial**2 * r_rho
        + source_sign * r_z
        + 2 * phase * radius)
    return (
        source_sign * angular * axial * poly / (4 * radius**2)
        - sp.I * mass * axial * (axial * a_rho - 2 * phase) / 4)


def _vacant_l0_minor(coefficient, axial, potential, s3, rho, z, source_sign, a_rho_major):
    occupied, vacant = _projector(source_sign, s3)
    l0 = _apply_L0(coefficient, axial, potential, s3, rho, z)
    major = _major_index(source_sign)
    minor = _minor_index(source_sign)
    on_shell = l0.subs(coefficient[major].diff(rho), a_rho_major)
    return axial**2 / (2 * sp.I) * (vacant * on_shell)[minor]


@lru_cache(maxsize=1)
def minor_a3_identities():
    """Closed L0 parent formulae for minor A2 and minor A3, both signs."""
    jets = _real_jets()
    rho, z = jets['rho'], jets['z']
    axial, radius, phase = jets['axial'], jets['radius'], jets['phase']
    height, real_major = jets['height'], jets['real_major']
    mass, angular = sp.symbols('m ell', real=True)
    _, _, s3 = _sym_pauli()
    potential = _sym_potential(mass, angular, radius)
    square = mass**2 + angular**2 / radius**2
    identities = {}
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        major = _major_index(source_sign)
        minor = _minor_index(source_sign)
        a1 = sp.zeros(2, 1)
        a1[major] = sp.I * phase
        a1[minor] = _minor_A1(axial, mass, angular, radius, source_sign)
        a1_rho = source_sign / axial**2 * a1[major].diff(z) - sp.I / 2 * square
        minor_a2 = _vacant_l0_minor(
            a1, axial, potential, s3, rho, z, source_sign, a1_rho)
        a_rho = axial.diff(rho)
        parent_a2 = _parent_minor_a2(
            axial, radius, phase, mass, angular, source_sign,
            a_rho, radius.diff(rho), radius.diff(z))
        identities['%s/minor_A2_parent' % tag] = _real_residual(minor_a2 - parent_a2)
        a2 = sp.zeros(2, 1)
        a2[major] = real_major + sp.I * height
        a2[minor] = minor_a2
        rhs_a2 = (sp.I / axial) * potential[major, minor] * minor_a2
        a2_rho = source_sign / axial**2 * a2[major].diff(z) + rhs_a2
        minor_a3 = _vacant_l0_minor(
            a2, axial, potential, s3, rho, z, source_sign, a2_rho)
        vacant_potential = -mass + sp.I * source_sign * angular / radius
        algebraic = -(axial / 2) * vacant_potential * (real_major + sp.I * height)
        from_t = axial**2 / (2 * sp.I) * _T(minor_a2, -source_sign, axial, rho, z)
        identities['%s/minor_A3_is_Tminus_minor_A2_plus_algebraic_major_A2' % tag] = (
            _real_residual(minor_a3 - from_t - algebraic))
        identities['%s/minor_A3_depends_on_major_A2' % tag] = (
            'depends' if real_major in minor_a3.atoms(sp.Function)
            or height in minor_a3.atoms(sp.Function) else 'independent')
        identities['%s/minor_A3_jet_independent_of_major_A2' % tag] = (
            'independent' if real_major not in from_t.atoms(sp.Function)
            and height not in from_t.atoms(sp.Function) else 'depends')
        massless_a3 = minor_a3.subs(mass, 0)
        parent_re = source_sign * axial * angular * height / (2 * radius)
        parent_im = (
            sp.im(sp.expand_complex(from_t.subs(mass, 0)))
            - source_sign * axial * angular * real_major / (2 * radius))
        identities['%s/massless_Re_minor_A3_is_s_a_ell_h_over_2r' % tag] = (
            _real_residual(sp.re(sp.expand_complex(massless_a3)) - parent_re))
        identities['%s/massless_Im_minor_A3_retains_R_and_Tminus_minor_A2' % tag] = (
            _real_residual(sp.im(sp.expand_complex(massless_a3)) - parent_im))
        identities['%s/massless_jet_independent_of_h' % tag] = (
            'independent' if height not in from_t.subs(mass, 0).atoms(sp.Function)
            else 'depends')
        q_rho_on_shell = (
            source_sign / axial**2 * phase.diff(z) - square / 2)
        q_z = phase.diff(z)
        ts_qz = _T(q_z, source_sign, axial, rho, z).subs(
            q_z.diff(rho), q_rho_on_shell.diff(z))
        identities['%s/T_s_q_z_on_shell_is_ell2_r_z_over_r3' % tag] = _real_residual(
            ts_qz - angular**2 * radius.diff(z) / radius**3)
    residual_keys = [
        key for key in identities
        if not key.endswith((
            'depends_on_major_A2', 'independent_of_major_A2', 'independent_of_h'))]
    leaves = []
    for key in residual_keys:
        value = identities[key]
        leaves.extend(value if isinstance(value, list) else [value])
    if set(leaves) != {'0'}:
        raise ArithmeticError(identities)
    for key, value in identities.items():
        if key.endswith('depends_on_major_A2') and value != 'depends':
            raise ArithmeticError('minor A3 must retain transported major A2')
        if key.endswith('independent_of_major_A2') and value != 'independent':
            raise ArithmeticError('T_{-s} minor A2 must stay independent of major A2')
        if key.endswith('independent_of_h') and value != 'independent':
            raise ArithmeticError('massless jet piece must stay independent of h')
    return MappingProxyType({
        'residuals': identities,
        'minor_A2': (
            '(s ell a / (4 r^2)) (-a_rho a r + a^2 r_rho + s r_z + 2 q r) '
            '- i (m a / 4) (a a_rho - 2 q)'),
        'minor_A3': (
            'a^2/(2 i) T_{-s}(minor A2) - (a/2) V_{-s,s} (R + i h), '
            'V_{-s,s} = -m + i s ell/r'),
        'massless_Re_minor_A3': 's a ell h / (2 r)',
        'massless_Im_minor_A3': (
            '-(a^2 / 2) T_{-s}(minor A2) - s a ell R / (2 r)'),
        'T_s': 'T_s = d_rho - s/a^2 d_z',
        'T_minus_s': 'T_{-s} = d_rho + s/a^2 d_z',
        'T_s_q': '-1/2 (m^2 + ell^2/r^2)',
        'T_s_q_z': 'ell^2 r_z / r^3',
        'dummy_symbol_L0_A_j': False,
        'A3_locally_algebraic_in_metric_jets': False,
        'major_A2_replaced_by_local_metric_jet': False,
        'massless_Im_Ts_major_A2_sets_upstream_datum_to_zero': False,
        'production_A3_A4_constructed': False,
        'leading_paired_e_minus2_cancels_for_equal_mu': True,
        'leading_e_minus2_cancellation_is_not_this_value': True,
    })


def second_geometry_phase_jet_majorants(
        family, mass, angular, rho_up=ARCHIVE_RHO_UP):
    """Interval majorants of second r jets and of q jets on the owned slab.

    Chebyshev value/first/second z-bounds already include the axial cutoff.
    Normal-window s-derivatives use |chi'|<=8/width and |chi''|<=105/width^2.
    Phase z-jets integrate T_s q_z = ell^2 r_z / r^3 from vanishing upstream
    major A1. Reference z-jets vanish. Evaluation roundoff is not included.
    """
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError('owned LocalIncomingFamily required')
    mass, angular = _require_explicit_mass_angular(mass, angular)
    geom = current_history_geometry_majorants(family, rho_up)
    bounds = radius_bounds(family, rho_up)
    width = geom['normal_window_width']
    chi_s, chi_ss = plateau_s_derivative_abs_uppers(width)
    sigma = nonnegative(geom['normal_support'], 'normal support')
    w_bounds = bounds['w']['profile_bounds']
    u_bounds = bounds['U']['profile_bounds']
    w0, w1, w2 = (nonnegative(value, 'w bound') for value in w_bounds[:3])
    u0, u1, u2 = (nonnegative(value, 'U bound') for value in u_bounds[:3])
    delta0 = sigma * w0 + sigma**3 * u0 / 6
    delta1 = sigma * w1 + sigma**3 * u1 / 6
    delta2 = sigma * w2 + sigma**3 * u2 / 6
    p_s = w0 + sigma**2 * u0 / 2
    p_ss = sigma * u0
    p_sz = w1 + sigma**2 * u1 / 2
    delta_s = chi_s * delta0 + p_s
    delta_ss = chi_ss * delta0 + 2 * chi_s * p_s + p_ss
    delta_sz = chi_s * delta1 + p_sz
    if delta_s != nonnegative(geom['delta_s_abs_upper'], 'delta_s'):
        raise ArithmeticError('delta_s majorant drifted from the transport owner')
    a_min = nonnegative(geom['axial_lower'], 'axial lower')
    a_hi = axial_abs_upper()
    a_rho = axial_rho_abs_upper()
    a_rr = axial_rho_rho_abs_upper()
    inv_a = 1 / a_min
    inv_a2 = inv_a**2
    r_min = nonnegative(geom['radius_lower'], 'radius lower')
    r_hi = nonnegative(geom['radius_upper'], 'radius upper')
    r_ref_min = nonnegative(geom['reference_radius_lower'], 'reference radius lower')
    r_ref_hi = nonnegative(geom['reference_radius_upper'], 'reference radius upper')
    r_ref_rho = nonnegative(geom['reference_r_rho_abs_upper'], 'reference r_rho')
    r_ref_rr = 1 / r_ref_min**3
    s_rho = inv_a
    s_rr = a_rho / a_min**2
    r_rho = nonnegative(geom['r_rho_abs_upper'], 'r_rho')
    r_z = nonnegative(geom['r_z_abs_upper'], 'r_z')
    r_zz = delta2
    r_rhoz = s_rho * delta_sz
    r_rr = r_ref_rr + s_rr * delta_s + inv_a2 * delta_ss
    duration = nonnegative(geom['duration'], 'duration')
    ell2 = angular**2
    q_source = (mass**2 + ell2 / r_min**2) / 2
    q_ref_source = (mass**2 + ell2 / r_ref_min**2) / 2
    q_z = duration * ell2 * r_z / r_min**3
    q_zz_source = ell2 * (r_zz / r_min**3 + 3 * r_z**2 / r_min**4)
    q_zz = duration * q_zz_source
    q_rho = q_source + inv_a2 * q_z
    q_rhoz = ell2 * r_z / r_min**3 + inv_a2 * q_zz
    tminus_q = q_source + 2 * inv_a2 * q_z
    q_ref_rho = q_ref_source
    if min(a_min, r_min, r_ref_min, duration) <= 0:
        raise ArithmeticError('geometry lower bounds must stay positive')
    if min(delta_s, chi_s, chi_ss, r_rr, r_ref_rr) < 0:
        raise ArithmeticError('second-jet majorant became negative')
    return MappingProxyType({
        'duration': duration,
        'axial_lower': a_min,
        'axial_upper': a_hi,
        'axial_abs_squared_upper': axial_abs_squared_upper(),
        'axial_rho_abs_upper': a_rho,
        'axial_rho_rho_abs_upper': a_rr,
        'inv_a_abs_upper': inv_a,
        'inv_a2_abs_upper': inv_a2,
        'radius_lower': r_min,
        'radius_upper': r_hi,
        'reference_radius_lower': r_ref_min,
        'reference_radius_upper': r_ref_hi,
        'reference_r_rho_abs_upper': r_ref_rho,
        'reference_r_rho_rho_abs_upper': r_ref_rr,
        'r_rho_abs_upper': r_rho,
        'r_z_abs_upper': r_z,
        'r_zz_abs_upper': r_zz,
        'r_rhoz_abs_upper': r_rhoz,
        'r_rho_rho_abs_upper': r_rr,
        'delta_radius_abs_upper': delta0,
        'delta_s_abs_upper': delta_s,
        'delta_ss_abs_upper': delta_ss,
        'delta_sz_abs_upper': delta_sz,
        'chi_s_abs_upper': chi_s,
        'chi_ss_abs_upper': chi_ss,
        'q_source_abs_upper': q_source,
        'q_reference_source_abs_upper': q_ref_source,
        'q_z_abs_upper': q_z,
        'q_zz_abs_upper': q_zz,
        'q_rho_abs_upper': q_rho,
        'q_rhoz_abs_upper': q_rhoz,
        'Tminus_q_abs_upper': tminus_q,
        'q_reference_z_abs_upper': Fraction(0),
        'q_reference_rho_abs_upper': q_ref_rho,
        'reference_r_z_abs_upper': Fraction(0),
        'reference_r_zz_abs_upper': Fraction(0),
        'reference_r_rhoz_abs_upper': Fraction(0),
        'deformation_vanishes': delta0 == 0,
        'profile_evaluation_roundoff_included': False,
        'sample_used_as_proof': False,
        'fitted_decay': False,
    })


def _poly_abs(a_hi, a_rho, r_hi, r_rho, r_z, q_abs):
    return a_rho * a_hi * r_hi + a_hi**2 * r_rho + r_z + 2 * q_abs * r_hi


def _poly_z_abs(a_hi, a_rho, r_z, r_rhoz, r_zz, r_hi, q_abs, q_z):
    return (
        a_rho * a_hi * r_z
        + a_hi**2 * r_rhoz
        + r_zz
        + 2 * q_z * r_hi
        + 2 * q_abs * r_z)


def _poly_rho_abs(
        a_hi, a_rho, a_rr, r_hi, r_rho, r_rr, r_rhoz, q_abs, q_rho):
    return (
        a_hi * a_rr * r_hi
        + a_rho**2 * r_hi
        + a_hi * a_rho * r_rho
        + a_hi**2 * r_rr
        + r_rhoz
        + 2 * q_rho * r_hi
        + 2 * q_abs * r_rho)


def _u_abs(a_hi, poly, r_lo):
    return a_hi * poly / r_lo**2


def _u_z_abs(a_hi, poly, poly_z, r_lo, r_hi, r_z):
    return a_hi * (poly_z * r_hi**2 + 2 * poly * r_hi * r_z) / r_lo**4


def _u_rho_abs(a_hi, a_rho, poly, poly_rho, r_lo, r_hi, r_rho):
    return (
        a_rho * poly / r_lo**2
        + a_hi * (poly_rho * r_hi**2 + 2 * poly * r_hi * r_rho) / r_lo**4)


def _mass_w_abs(a_hi, a_rho, q_abs):
    return a_hi * (a_hi * a_rho + 2 * q_abs)


def _mass_w_z_abs(a_hi, q_z):
    return 2 * a_hi * q_z


def _mass_w_rho_abs(a_hi, a_rho, a_rr, q_abs, q_rho):
    return 2 * a_hi * a_rho**2 + a_hi**2 * a_rr + 2 * a_rho * q_abs + 2 * a_hi * q_rho


def minor_a3_interval_majorants(family, mass, angular, rho_up=ARCHIVE_RHO_UP):
    """Interval majorants of minor A2 and of the minor-A3 jet piece.

    The algebraic coupling to generated Re(major A2) is bounded. The shared
    upstream major A2 times delta(1/r) is named and not replaced by zero.
    Triangle differences collapse to exact 0 when the radius deformation
    vanishes. Fitted decay and samples are not used.
    """
    mass, angular = _require_explicit_mass_angular(mass, angular)
    if mass != 0:
        raise ValueError('minor A3 difference coupling is currently bounded only for massless channels')
    ell_abs = abs(angular)
    jets = second_geometry_phase_jet_majorants(family, mass, angular, rho_up)
    phase = phase_and_real_major_A2_difference_majorants(
        family, mass, angular, rho_up)
    a_hi = jets['axial_upper']
    a2_hi = jets['axial_abs_squared_upper']
    a_rho = jets['axial_rho_abs_upper']
    a_rr = jets['axial_rho_rho_abs_upper']
    inv_a2 = jets['inv_a2_abs_upper']
    q_g = nonnegative(phase['q_abs_upper'], 'q')
    q_ref = nonnegative(phase['q_reference_abs_upper'], 'q_ref')
    delta_q = nonnegative(phase['delta_q_abs_upper'], 'delta_q')
    delta_r = nonnegative(phase['history_minus_reference_real_major_A2_abs_upper'], 'delta R')
    f_g = nonnegative(phase['Re_Ts_major_A2_g_abs_upper'], 'F_g')
    f_ref = nonnegative(phase['Re_Ts_major_A2_ref_abs_upper'], 'F_ref')
    duration = jets['duration']
    i_g = duration * f_g
    i_ref = duration * f_ref
    r_min = jets['radius_lower']
    r_hi = jets['radius_upper']
    r_ref_min = jets['reference_radius_lower']
    r_ref_hi = jets['reference_radius_upper']
    poly_g = _poly_abs(
        a_hi, a_rho, r_hi, jets['r_rho_abs_upper'], jets['r_z_abs_upper'], q_g)
    poly_ref = _poly_abs(
        a_hi, a_rho, r_ref_hi, jets['reference_r_rho_abs_upper'], Fraction(0), q_ref)
    poly_z_g = _poly_z_abs(
        a_hi, a_rho, jets['r_z_abs_upper'], jets['r_rhoz_abs_upper'],
        jets['r_zz_abs_upper'], r_hi, q_g, jets['q_z_abs_upper'])
    poly_rho_g = _poly_rho_abs(
        a_hi, a_rho, a_rr, r_hi, jets['r_rho_abs_upper'],
        jets['r_rho_rho_abs_upper'], jets['r_rhoz_abs_upper'],
        q_g, jets['q_rho_abs_upper'])
    poly_rho_ref = _poly_rho_abs(
        a_hi, a_rho, a_rr, r_ref_hi, jets['reference_r_rho_abs_upper'],
        jets['reference_r_rho_rho_abs_upper'], Fraction(0),
        q_ref, jets['q_reference_rho_abs_upper'])
    u_g = _u_abs(a_hi, poly_g, r_min)
    u_ref = _u_abs(a_hi, poly_ref, r_ref_min)
    u_z_g = _u_z_abs(
        a_hi, poly_g, poly_z_g, r_min, r_hi, jets['r_z_abs_upper'])
    u_rho_g = _u_rho_abs(
        a_hi, a_rho, poly_g, poly_rho_g, r_min, r_hi, jets['r_rho_abs_upper'])
    u_rho_ref = _u_rho_abs(
        a_hi, a_rho, poly_ref, poly_rho_ref, r_ref_min, r_ref_hi,
        jets['reference_r_rho_abs_upper'])
    t_u_g = u_rho_g + inv_a2 * u_z_g
    t_u_ref = u_rho_ref
    w_g = _mass_w_abs(a_hi, a_rho, q_g)
    w_ref = _mass_w_abs(a_hi, a_rho, q_ref)
    t_w_g = (
        _mass_w_rho_abs(a_hi, a_rho, a_rr, q_g, jets['q_rho_abs_upper'])
        + inv_a2 * _mass_w_z_abs(a_hi, jets['q_z_abs_upper']))
    t_w_ref = _mass_w_rho_abs(
        a_hi, a_rho, a_rr, q_ref, jets['q_reference_rho_abs_upper'])
    ell_factor = ell_abs / 4
    mass_factor = mass / 4
    minor_g = ell_factor * u_g + mass_factor * w_g
    minor_ref = ell_factor * u_ref + mass_factor * w_ref
    t_g = ell_factor * t_u_g + mass_factor * t_w_g
    t_ref = ell_factor * t_u_ref + mass_factor * t_w_ref
    jet_g = a2_hi * t_g / 2
    jet_ref = a2_hi * t_ref / 2
    vanishes = jets['deformation_vanishes']
    if vanishes:
        delta_minor = Fraction(0)
        delta_t = Fraction(0)
        delta_jet = Fraction(0)
        generated = Fraction(0)
        kernel = Fraction(0)
    else:
        dinv2 = jets['delta_radius_abs_upper'] * (
            2 * r_ref_hi + jets['delta_radius_abs_upper']
        ) / (r_min**2 * r_ref_min**2)
        delta_r_rho = jets['inv_a_abs_upper'] * jets['delta_s_abs_upper']
        delta_poly = (
            a_rho * a_hi * jets['delta_radius_abs_upper']
            + a_hi**2 * delta_r_rho
            + jets['r_z_abs_upper']
            + 2 * (q_g * jets['delta_radius_abs_upper'] + r_ref_hi * delta_q))
        delta_u = poly_ref * dinv2 + delta_poly / r_min**2
        delta_w = a_hi * 2 * delta_q
        delta_minor = ell_factor * delta_u + mass_factor * delta_w
        delta_t = t_g + t_ref
        delta_jet = jet_g + jet_ref
        coupling = a_hi * ell_abs / 2
        generated = coupling * (
            delta_r / r_min
            + i_ref * jets['delta_radius_abs_upper'] / (r_min * r_ref_min))
        kernel = coupling * jets['delta_radius_abs_upper'] / (r_min * r_ref_min)
    if min(minor_g, minor_ref, t_g, t_ref, jet_g, jet_ref) < 0:
        raise ArithmeticError('minor A3 majorant became negative')
    if min(delta_minor, delta_t, delta_jet, generated, kernel) < 0:
        raise ArithmeticError('history-minus-reference majorant became negative')
    complete = Fraction(0) if vanishes else None
    return MappingProxyType({
        'mass': mass,
        'angular': angular,
        'angular_abs': ell_abs,
        'minor_A2_g_abs_upper': minor_g,
        'minor_A2_ref_abs_upper': minor_ref,
        'history_minus_reference_minor_A2_abs_upper': delta_minor,
        'Tminus_minor_A2_g_abs_upper': t_g,
        'Tminus_minor_A2_ref_abs_upper': t_ref,
        'history_minus_reference_Tminus_minor_A2_abs_upper': delta_t,
        'minor_A3_jet_g_abs_upper': jet_g,
        'minor_A3_jet_ref_abs_upper': jet_ref,
        'history_minus_reference_minor_A3_jet_abs_upper': delta_jet,
        'generated_algebraic_major_A2_coupling_abs_upper': generated,
        'generated_real_major_A2_g_abs_upper': i_g,
        'generated_real_major_A2_ref_abs_upper': i_ref,
        'upstream_major_A2_difference_kernel_abs_upper': kernel,
        'complete_history_minus_reference_minor_A3_abs_upper': complete,
        'complete_full_envelope_minor_A3_abs_upper': None,
        'upstream_major_A2_owned': False,
        'full_envelope_real_major_A2': None,
        'full_envelope_imag_major_A2': None,
        'deformation_vanishes': vanishes,
        'fitted_decay': False,
        'sample_used_as_proof': False,
        'jets': jets,
        'phase': phase,
    })


def current_history_minor_a3_report(root=None):
    """Executable minor-A3 jet bound on live 0b0e4ced and the smallest pair."""
    root = Path(root) if root is not None else _root()
    history = json.loads((root / HISTORY_RECORD).read_text())
    family = LocalIncomingFamily(np.array(history['history']['coefficients']))
    identity = profile_identity(family)
    if identity != CURRENT_HISTORY_IDENTITY:
        raise ValueError('expected exact live history 0b0e4ced')
    pair = smallest_original_pair(root)
    identities = minor_a3_identities()
    majorants = minor_a3_interval_majorants(
        family, pair['mass'], pair['absolute_angular'])
    if majorants['deformation_vanishes']:
        raise ArithmeticError('live 0b0e4ced deformation may not vanish')
    if majorants['complete_history_minus_reference_minor_A3_abs_upper'] is not None:
        raise ArithmeticError('complete delta minor A3 is not owned on 0b0e4ced')
    _require_missing(majorants['complete_full_envelope_minor_A3_abs_upper'], 'full minor A3')
    _require_missing(majorants['full_envelope_real_major_A2'], 'full-envelope real major A2')
    _require_missing(majorants['full_envelope_imag_major_A2'], 'full-envelope imag major A2')
    jets = majorants['jets']
    phase = majorants['phase']
    return MappingProxyType({
        'schema': SCHEMA,
        'profile_identity': identity,
        'pair': dict(pair),
        'identities': identities,
        'first_noncancelling_paired_order': FIRST_NONCANCELLING_PAIRED_ORDER,
        'leading_paired_e_minus2_cancels_for_equal_mu': True,
        'leading_e_minus2_cancellation_is_not_this_value': True,
        'e_minus3_first_allowed_order_is_not_a_tail_bound': True,
        'A3_locally_algebraic_in_metric_jets': False,
        'major_A2_replaced_by_local_metric_jet': False,
        'dummy_symbol_L0_A_j': False,
        'production_A3_A4_constructed': False,
        'massless_Im_Ts_major_A2_sets_upstream_datum_to_zero': False,
        'upstream_major_A2_owned': False,
        'coefficient_scope': 'history-minus-reference',
        'majorants': MappingProxyType({
            'mass': float(majorants['mass']),
            'angular_abs': bound_record(majorants['angular_abs']),
            'q_z_abs_upper': bound_record(jets['q_z_abs_upper']),
            'q_rho_abs_upper': bound_record(jets['q_rho_abs_upper']),
            'q_zz_abs_upper': bound_record(jets['q_zz_abs_upper']),
            'Tminus_q_abs_upper': bound_record(jets['Tminus_q_abs_upper']),
            'r_zz_abs_upper': bound_record(jets['r_zz_abs_upper']),
            'r_rhoz_abs_upper': bound_record(jets['r_rhoz_abs_upper']),
            'r_rho_rho_abs_upper': bound_record(jets['r_rho_rho_abs_upper']),
            'axial_rho_rho_abs_upper': bound_record(jets['axial_rho_rho_abs_upper']),
            'delta_ss_abs_upper': bound_record(jets['delta_ss_abs_upper']),
            'minor_A2_g_abs_upper': bound_record(majorants['minor_A2_g_abs_upper']),
            'minor_A2_ref_abs_upper': bound_record(majorants['minor_A2_ref_abs_upper']),
            'history_minus_reference_minor_A2_abs_upper': bound_record(
                majorants['history_minus_reference_minor_A2_abs_upper']),
            'Tminus_minor_A2_g_abs_upper': bound_record(
                majorants['Tminus_minor_A2_g_abs_upper']),
            'history_minus_reference_Tminus_minor_A2_abs_upper': bound_record(
                majorants['history_minus_reference_Tminus_minor_A2_abs_upper']),
            'minor_A3_jet_g_abs_upper': bound_record(
                majorants['minor_A3_jet_g_abs_upper']),
            'minor_A3_jet_ref_abs_upper': bound_record(
                majorants['minor_A3_jet_ref_abs_upper']),
            'history_minus_reference_minor_A3_jet_abs_upper': bound_record(
                majorants['history_minus_reference_minor_A3_jet_abs_upper']),
            'generated_algebraic_major_A2_coupling_abs_upper': bound_record(
                majorants['generated_algebraic_major_A2_coupling_abs_upper']),
            'generated_real_major_A2_g_abs_upper': bound_record(
                majorants['generated_real_major_A2_g_abs_upper']),
            'upstream_major_A2_difference_kernel_abs_upper': bound_record(
                majorants['upstream_major_A2_difference_kernel_abs_upper']),
            'q_abs_upper': bound_record(phase['q_abs_upper']),
            'delta_q_abs_upper': bound_record(phase['delta_q_abs_upper']),
            'history_minus_reference_real_major_A2_abs_upper': bound_record(
                phase['history_minus_reference_real_major_A2_abs_upper']),
            'complete_history_minus_reference_minor_A3_abs_upper': None,
            'complete_full_envelope_minor_A3_abs_upper': None,
            'full_envelope_real_major_A2': None,
            'fitted_decay': False,
            'sample_used_as_proof': False,
        }),
        'numerical_C4': None,
        'numerical_C_M': None,
        'vacuum_integrated_tail_N_beta': None,
        'missing_primitive': MISSING_UPSTREAM_MAJOR_A2,
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
    })


def validate_minor_a3_report(report, root=None):
    """Reject invented C4/C_M, null-to-zero A2_up, or samples as proof."""
    report = _as_map(report, 'minor A3 report')
    if report.get('schema') != SCHEMA:
        raise ValueError('unexpected minor A3 schema')
    if report.get('background_tail_certificate_transferred'):
        raise ValueError('background tail certificate must not transfer')
    if report.get('fitted_decay') or report.get('sample_used_as_proof'):
        raise ValueError('fitted decay or samples may not prove a minor A3 bound')
    if report.get('manufactured_fixture_used_as_production_constant'):
        raise ValueError('manufactured fixture cannot supply a production constant')
    if report.get('field_or_source_evolutions') not in (0, None):
        raise ValueError('minor A3 control may not launch field or source evolutions')
    if report.get('physical_EXISTENCE_certificate') or report.get(
            'physical_NONEXISTENCE_certificate'):
        raise ValueError('OPEN calculation may not be reported as PASS or NON_EXISTENCE')
    if report.get('physical_local_gate') != 'OPEN':
        raise ValueError('local gate remains OPEN')
    if report.get('v1_pilot_rewritten') or report.get('v2_record_rewritten') or report.get(
            'v3_record_rewritten'):
        raise ValueError('historical UV remainder bytes are immutable')
    identities = _as_map(report.get('identities'), 'identities')
    if identities.get('dummy_symbol_L0_A_j'):
        raise ValueError('minor A3 may not use dummy L0 A_j symbols')
    if identities.get('production_A3_A4_constructed') or report.get(
            'production_A3_A4_constructed'):
        raise ValueError('production A3/A4 are not constructed here')
    if identities.get('A3_locally_algebraic_in_metric_jets') or report.get(
            'A3_locally_algebraic_in_metric_jets'):
        raise ValueError('minor A3 is not a local metric jet')
    if identities.get('massless_Im_Ts_major_A2_sets_upstream_datum_to_zero'):
        raise ValueError('massless Im(T_s major A2) does not own the upstream datum')
    if not identities.get('leading_e_minus2_cancellation_is_not_this_value'):
        raise ValueError('E^{-2} cancellation is not this bound')
    if report.get('first_noncancelling_paired_order') != FIRST_NONCANCELLING_PAIRED_ORDER:
        raise ValueError('E^{-3} first allowed paired order was not preserved')
    if report.get('leading_paired_e_minus2_cancels_for_equal_mu') is False:
        raise ValueError('E^{-2} equal-mu cancellation was not preserved')
    if report.get('coefficient_scope') != 'history-minus-reference':
        raise ValueError('wrong coefficient scope')
    if report.get('upstream_major_A2_owned'):
        raise ValueError('upstream major A2 is not owned')
    _require_missing(report.get('numerical_C4'), 'C4')
    _require_missing(report.get('numerical_C_M'), 'C_M')
    _require_missing(report.get('vacuum_integrated_tail_N_beta'), 'vacuum tail')
    majorants = _as_map(report.get('majorants'), 'majorants')
    _require_missing(
        majorants.get('complete_history_minus_reference_minor_A3_abs_upper'),
        'complete delta minor A3')
    _require_missing(
        majorants.get('complete_full_envelope_minor_A3_abs_upper'),
        'full-envelope minor A3')
    _require_missing(
        majorants.get('full_envelope_real_major_A2'), 'full-envelope real major A2')
    if majorants.get('fitted_decay') or majorants.get('sample_used_as_proof'):
        raise ValueError('fitted decay or samples may not prove a minor A3 bound')
    for name in (
            'q_z_abs_upper', 'r_rho_rho_abs_upper',
            'minor_A2_g_abs_upper', 'minor_A3_jet_g_abs_upper',
            'history_minus_reference_minor_A3_jet_abs_upper',
            'generated_algebraic_major_A2_coupling_abs_upper',
            'upstream_major_A2_difference_kernel_abs_upper'):
        number = nonnegative(
            majorants[name]['binary64_upper'] if isinstance(majorants[name], Mapping)
            else majorants[name],
            name)
        if number <= 0:
            raise ValueError('live-history majorant must stay strictly positive: ' + name)
    if report.get('profile_identity') != CURRENT_HISTORY_IDENTITY:
        raise ValueError('minor A3 report is not bound to history 0b0e4ced')
    pair = _as_map(report.get('pair'), 'pair')
    if float(pair.get('mass', 1.0)) != 0.0:
        raise ValueError('smallest representative pair must be massless')
    if float(pair.get('absolute_angular', 0.0)) <= 0.0:
        raise ValueError('smallest representative pair must be a nonzero-angular family')
    return True
