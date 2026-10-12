"""Massive e^{-3} cubic coefficient on the unchanged directed integrator.

The physical action, original source sectors and history 0b0e4ced stay fixed.
Massless owners are imported and not modified. This proves the generic-mass
transport and contracts one original massive channel at one incoming point.
It does not enclose the incoming interval, the higher UV remainder, or a
physical tail. No complete C_M is constructed.
"""
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
import json

import numpy as np
import sympy as sp
from flint import arb, ctx

from . import nsc_ks_ball_geometry as G
from .nsc_dirac_source_phase import SOURCE_SIGNS
from .nsc_evolved_incoming_constraints import source_column_matter
from .nsc_ks_cubic_uv_current import (
    _zero,
    confirm_action_signs_against_source_column_matter,
    reduced_massless_transport_identities,
)
from .nsc_ks_cubic_uv_enclosure import (
    CubicForcing,
    CubicGeometry,
    enclose_characteristic,
)
from .nsc_ks_current_uv_transport import (
    ARCHIVE_RHO_UP,
    CAUCHY_RECORD,
    CURRENT_HISTORY_IDENTITY,
    CURRENT_HISTORY_IDENTITY_PREFIX,
    HISTORY_RECORD,
)
from .nsc_ks_finite_history_uv_coefficients import (
    _apply_L0,
    _major_index,
    _minor_A1,
    _minor_index,
    _projector,
    _sym_pauli,
    _sym_potential,
    _T,
)
from .nsc_ks_profile_identity import profile_identity
from .nsc_local_incoming_family import LocalIncomingFamily


SCHEMA = 'NSC-KS-MASSIVE-CUBIC-UV-v1'
MASSLESS_J3 = (
    '-a^2 b b_rhorhoz+a^2 b_rho b_rhoz-2 a a_rho b b_rhoz'
    '-s b b_rhozz+s b_z b_rhoz-16 b^3 b_z/a^2+4 s b^2 q_zz/a^2')
MASSLESS_QZ = '-4 b b_z/a^2'
MASSLESS_QZZ = '-4 (b_z^2+b b_zz)/a^2'
INVENTORY_RECORD = 'results/development/nsc-ks-source-inventory.json'
CUTOFF_BRIDGE_RECORD = 'results/development/nsc-ks-cutoff-bridge-control.json'
GROUP14 = 14


def _lab_root():
    return Path(__file__).resolve().parents[2]


def _exact(expr, name):
    value = sp.expand(expr)
    if value != 0:
        value = sp.simplify(sp.together(value))
    if value != 0:
        raise ArithmeticError(name + ' failed: ' + str(value))
    return '0'


def _operator_blocks():
    """L0 minors, massive h transport, continuity and the h3/n4 combination."""
    rho, z = sp.symbols('rho z', real=True)
    axial = sp.Function('a', positive=True)(rho)
    radius = sp.Function('r', positive=True)(rho, z)
    phase = sp.Function('q', real=True)(rho, z)
    height = sp.Function('h', real=True)(rho, z)
    real_major = sp.Function('R', real=True)(rho, z)
    height3 = sp.Function('h3', real=True)(rho, z)
    norm4 = sp.Function('n4', real=True)(rho, z)
    mass, angular, moment = sp.symbols('m ell k', real=True)
    _, _, s3 = _sym_pauli()
    potential = _sym_potential(mass, angular, radius)
    square = mass**2 + angular**2 / radius**2
    checks = {}
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        major, minor = _major_index(source_sign), _minor_index(source_sign)
        _, vacant = _projector(source_sign, s3)
        slope = axial * angular / (2 * radius)
        shift = axial * mass / 2
        minor_a1 = shift - sp.I * source_sign * slope
        checks['%s/minor_A1' % tag] = _zero(
            minor_a1 - _minor_A1(axial, mass, angular, radius, source_sign))
        transport = lambda scalar, sign: _T(scalar, sign, axial, rho, z)
        opposite = lambda scalar: transport(scalar, -source_sign)
        minor_formula = sp.I * (
            phase * minor_a1 - axial**2 * opposite(minor_a1) / 2)
        column = sp.zeros(2, 1)
        column[major] = sp.I * phase
        column[minor] = minor_a1
        on_shell = (
            source_sign / axial**2 * column[major].diff(z) - sp.I / 2 * square)
        solved = _apply_L0(column, axial, potential, s3, rho, z).subs(
            column[major].diff(rho), on_shell)
        minor_a2 = axial**2 / (2 * sp.I) * (vacant * solved)[minor]
        checks['%s/minor_A2_from_L0' % tag] = _zero(minor_a2 - minor_formula)
        real_part = (
            source_sign * slope * phase
            - source_sign * axial**2 * slope.diff(rho) / 2
            - slope.diff(z) / 2)
        imag_part = shift * (phase - axial * axial.diff(rho) / 2)
        checks['%s/minor_A2_real_imag' % tag] = _zero(
            minor_formula - (real_part + sp.I * imag_part))
        checks['%s/L_q' % tag] = _zero(
            (sp.I / axial) * potential[major, minor] * minor_a1 + sp.I / 2 * square)
        checks['%s/b_matches_angular_square' % tag] = _exact(
            2 * slope**2 / axial**2 - angular**2 / (2 * radius**2), tag + ' angular')
        major_a2 = real_major + sp.I * height
        minor_formula_c = (
            -sp.I * axial**2 * opposite(minor_formula) / 2
            + minor_a1 * major_a2)
        second = sp.zeros(2, 1)
        second[major] = major_a2
        second[minor] = minor_formula
        rhs_a2 = (sp.I / axial) * potential[major, minor] * minor_formula
        second_rho = source_sign / axial**2 * second[major].diff(z) + rhs_a2
        solved = _apply_L0(second, axial, potential, s3, rho, z).subs(
            second[major].diff(rho), second_rho)
        minor_a3 = axial**2 / (2 * sp.I) * (vacant * solved)[minor]
        checks['%s/minor_A3_from_L0' % tag] = _zero(minor_a3 - minor_formula_c)
        parent_h = -mass * angular / 4 * (
            source_sign * axial**2 * radius.diff(rho) + radius.diff(z)
        ) / radius**2
        checks['%s/L_h_parent' % tag] = _zero(
            sp.im(sp.expand_complex(rhs_a2)) - parent_h)
        checks['%s/L_h_massless_limit' % tag] = _zero(parent_h.subs(mass, 0))
        # A unit axial jet makes the parent numerator 1, so L(h) is not the
        # zero expression and the massless equation L(h)=0 is not imposed.
        checks['%s/L_h_nonzero_unit_jet' % tag] = _exact(
            (source_sign * axial**2 * radius.diff(rho) + radius.diff(z)).subs({
                axial: 1, radius.diff(rho): 0, radius.diff(z): 1,
            }) - 1, tag + ' h jet')
        symbol = sp.symbols('C')
        checks['%s/L_h3_from_V' % tag] = _zero(
            sp.im(sp.expand_complex(
                (sp.I / axial) * potential[major, minor] * symbol))
            + 2 * sp.re(sp.expand_complex(sp.conjugate(minor_a1) * symbol))
            / axial**2)
        energy = sp.symbols('e', positive=True)
        density_n2 = sp.Function('n2', real=True)(rho, z)
        density_n3 = sp.Function('n3', real=True)(rho, z)
        first = sp.Function('Lmin')(rho, z)
        second_minor = sp.Function('Bmin')(rho, z)
        third = sp.Function('Cmin')(rho, z)
        minor_square = (
            sp.conjugate(first) * first / energy**2
            + 2 * sp.re(sp.conjugate(first) * second_minor) / energy**3
            + (sp.conjugate(second_minor) * second_minor
               + 2 * sp.re(sp.conjugate(first) * third)) / energy**4)
        density = (1 + density_n2 / energy**2 + density_n3 / energy**3
                   + norm4 / energy**4)
        current = source_sign * (density - 2 * minor_square)
        order4 = sp.expand(
            (density.diff(rho) - current.diff(z) / axial**2) * energy**4
        ).subs(energy, 0)
        expected = transport(norm4, source_sign) + source_sign / axial**2 * sp.diff(
            2 * sp.conjugate(second_minor) * second_minor
            + 4 * sp.re(sp.conjugate(first) * third), z)
        checks['%s/continuity_n4' % tag] = _zero(order4 - expected)
        re_fun = sp.Function('U', real=True)(rho, z)
        abs_fun = sp.Function('P', real=True)(rho, z)
        source_h = -2 * re_fun / axial**2
        source_n = -source_sign / axial**2 * sp.diff(2 * abs_fun + 4 * re_fun, z)
        combo = (sp.diff(height3, z, rho) - source_sign * sp.diff(norm4, rho)
                 - source_sign / axial**2 * sp.diff(
                     sp.diff(height3, z) - source_sign * norm4, z))
        combo = combo.subs({
            sp.diff(height3, rho, z): (
                sp.diff(source_h, z) + source_sign / axial**2 * sp.diff(height3, z, 2)),
            sp.diff(norm4, rho): source_n + source_sign / axial**2 * sp.diff(norm4, z),
        })
        checks['%s/L_h3z_minus_s_n4' % tag] = _exact(
            combo - 2 / axial**2 * sp.diff(abs_fun + re_fun, z), tag + ' combo')
        sample = sp.Function('f', real=True)(rho, z)
        checks['%s/L_commutes_with_z' % tag] = _exact(
            transport(sample.diff(z), source_sign)
            - transport(sample, source_sign).diff(z), tag + ' commute')
    return checks


def _transport_blocks():
    """Characteristic reduction. h is already absent; L(h)=0 is not used."""
    rho, z = sp.symbols('rho z', real=True)
    axial = sp.Function('a', positive=True)(rho)
    slope = sp.Function('b', real=True)(rho, z)
    phase = sp.Function('q', real=True)(rho, z)
    mass, moment = sp.symbols('m k', real=True)
    checks = {}
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        shift = axial * mass / 2
        minor_l = shift - sp.I * source_sign * slope
        direct = lambda scalar: sp.diff(scalar, rho) - source_sign / axial**2 * sp.diff(scalar, z)
        reverse = lambda scalar: sp.diff(scalar, rho) + source_sign / axial**2 * sp.diff(scalar, z)
        real_b = (
            source_sign * slope * phase
            - source_sign * axial**2 * sp.diff(slope, rho) / 2
            - sp.diff(slope, z) / 2)
        imag_b = shift * (phase - axial * sp.diff(axial, rho) / 2)
        major_r = (
            source_sign * sp.diff(phase, z) + moment - phase**2
            - shift**2 - slope**2) / 2
        stored = sp.symbols('MBR MBI R h', real=True)
        # Arbitrary M(B), R and h: the h piece of Re(conj(l) C) is identically 0.
        trial = -sp.I * axial**2 * (stored[0] + sp.I * stored[1]) / 2 + minor_l * (
            stored[2] + sp.I * stored[3])
        hand_re = (
            axial**2 / 2 * (shift * stored[1] + source_sign * slope * stored[0])
            + (shift**2 + slope**2) * stored[2])
        checks['%s/Re_lC_independent_of_h' % tag] = _exact(
            sp.re(sp.expand_complex(sp.conjugate(minor_l) * trial)) - hand_re,
            tag + ' Re')
        produced = (
            axial**2 / 2 * (shift * reverse(imag_b) + source_sign * slope * reverse(real_b))
            + (shift**2 + slope**2) * major_r)
        rest = (
            major_r * sp.diff(phase, z) - phase * sp.diff(major_r, z)
            + shift * sp.diff(imag_b, z)
            + source_sign * slope * sp.diff(real_b, z)
            - source_sign * real_b * sp.diff(slope, z))
        generated = 2 / axial**2 * sp.diff(real_b**2 + imag_b**2 + produced, z) + direct(rest)
        on_shell = (
            source_sign / axial**2 * sp.diff(phase, z) - mass**2 / 2
            - 2 * slope**2 / axial**2)
        substitute = {
            sp.diff(phase, rho, 1, z, order): sp.diff(on_shell, z, order)
            for order in range(4)}
        generated = sp.expand(sp.expand(generated).subs(substitute))
        generated = sp.expand(generated.subs(substitute))
        massless = (
            -axial**2 * slope * sp.diff(slope, rho, 2, z, 1)
            + axial**2 * sp.diff(slope, rho) * sp.diff(slope, rho, z)
            - 2 * axial * sp.diff(axial, rho) * slope * sp.diff(slope, rho, z)
            - source_sign * slope * sp.diff(slope, rho, 1, z, 2)
            + source_sign * sp.diff(slope, z) * sp.diff(slope, rho, z)
            - 16 * slope**3 * sp.diff(slope, z) / axial**2
            + 4 * source_sign * slope**2 * sp.diff(phase, z, 2) / axial**2)
        correction = mass**2 * (
            source_sign * sp.diff(phase, z, 2) - 4 * slope * sp.diff(slope, z))
        checks['%s/L_J3' % tag] = _exact(
            generated - massless - correction, tag + ' J3')
        checks['%s/L_J3_massless_reduction' % tag] = _exact(
            generated.subs(mass, 0) - massless, tag + ' m0')
        checks['%s/L_qz' % tag] = _exact(
            direct(sp.diff(phase, z)).subs(substitute)
            + 4 * slope * sp.diff(slope, z) / axial**2, tag + ' qz')
        checks['%s/L_qzz' % tag] = _exact(
            direct(sp.diff(phase, z, 2)).subs(substitute)
            + 4 * (sp.diff(slope, z)**2 + slope * sp.diff(slope, z, 2)) / axial**2,
            tag + ' qzz')
        sample = correction.subs({mass: 1, sp.diff(phase, z, 2): 1, slope: 0})
        checks['%s/mass_correction_not_zero' % tag] = _exact(
            sample - source_sign, tag + ' mass sample')
        checks['%s/generated_has_no_h' % tag] = '0' if not generated.has(sp.Function('h')) else 'depends'
    if any(value != '0' for key, value in checks.items() if key.endswith('generated_has_no_h')):
        raise ArithmeticError('h remained in the reduced current after cancellation')
    return checks


def _surface_blocks():
    """Pauli N bracket minus s J3/a, then the matched-radius Sigma difference."""
    checks = {}
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        names = (
            'a m b bp bpz bz bzz q qz qzz k h hz h3z n4 ap '
            'MBR MBI R3 h3 R3z Cr Ci').split()
        (axial, mass, slope, slope_rho, slope_rhoz, slope_z, slope_zz,
         phase, phase_z, phase_zz, moment, height, height_z, current_h, norm4,
         axial_rho, flow_r, flow_i, real3, height3, real3_z, real_c, imag_c
         ) = sp.symbols(' '.join(names), real=True)
        shift = axial * mass / 2
        minor_l = shift - sp.I * source_sign * slope
        real_b = (source_sign * slope * phase - source_sign * axial**2 * slope_rho / 2
                  - slope_z / 2)
        imag_b = shift * (phase - axial * axial_rho / 2)
        minor_b = real_b + sp.I * imag_b
        minor_bz = (
            source_sign * slope_z * phase + source_sign * slope * phase_z
            - source_sign * axial**2 * slope_rhoz / 2 - slope_zz / 2
            + sp.I * shift * phase_z)
        minor_z = -sp.I * source_sign * slope_z
        major_r = (source_sign * phase_z + moment - phase**2 - shift**2 - slope**2) / 2
        major_rz = (source_sign * phase_zz - 2 * phase * phase_z - 2 * slope * slope_z) / 2
        minor_c = (
            -sp.I * axial**2 * (flow_r + sp.I * flow_i) / 2
            + minor_l * (major_r + sp.I * height))
        major, minor = _major_index(source_sign), _minor_index(source_sign)
        _, _, s3 = _sym_pauli()
        potential = sp.zeros(2)
        potential[minor, major] = -2 * minor_l / axial
        potential[major, minor] = sp.conjugate(potential[minor, major])
        columns = [sp.zeros(2, 1) for _ in range(4)]
        derivatives = [sp.zeros(2, 1) for _ in range(4)]
        columns[0][major] = 1
        columns[1][major] = sp.I * phase
        columns[1][minor] = minor_l
        columns[2][major] = major_r + sp.I * height
        columns[2][minor] = minor_b
        columns[3][major] = real3 + sp.I * height3
        columns[3][minor] = minor_c
        derivatives[1][major] = sp.I * phase_z
        derivatives[1][minor] = minor_z
        derivatives[2][major] = major_rz + sp.I * height_z
        derivatives[2][minor] = minor_bz
        derivatives[3][major] = real3_z + sp.I * current_h
        derivatives[3][minor] = real_c + sp.I * imag_c
        raw = sum((columns[index].H * derivatives[3 - index])[0] for index in range(4))
        identity = sp.im(sp.expand_complex(raw)) - source_sign * norm4
        parent = (
            current_h + major_r * phase_z - phase * major_rz
            + sp.im(sp.expand_complex(
                sp.conjugate(minor_l) * minor_bz + sp.conjugate(minor_b) * minor_z))
            - source_sign * norm4)
        checks['%s/J3_pauli' % tag] = _exact(identity - parent, tag + ' J3')
        checks['%s/J3_h_coefficient' % tag] = _exact(sp.diff(identity, height), tag + ' h')
        checks['%s/J3_hz_coefficient' % tag] = _exact(
            sp.diff(identity, height_z), tag + ' hz')
        checks['%s/raw_h_coefficient_is_q_z' % tag] = _exact(
            sp.expand_complex(sp.diff(raw, height)) - phase_z, tag + ' raw h')
        checks['%s/raw_hz_coefficient_is_q' % tag] = _exact(
            sp.expand_complex(sp.diff(raw, height_z)) - phase, tag + ' raw hz')
        spin = sum((columns[index].H * s3 * derivatives[3 - index])[0] for index in range(4))
        real_cross = sp.re(sp.expand_complex(sp.conjugate(minor_l) * minor_c))
        minor_square = sp.expand_complex(minor_b * sp.conjugate(minor_b)) + 2 * real_cross
        spin_current = sp.im(sp.expand_complex(spin)) - (
            norm4 - 2 * minor_square)
        potential_density = sp.re(sp.expand_complex(sum(
            (columns[index].H * potential * columns[3 - index])[0] for index in range(4))))
        remainder = sp.expand(
            potential_density + spin_current / axial - source_sign * identity / axial)
        jet = -source_sign * axial**2 * slope_rho / 2 - slope_z / 2
        jet_z = -source_sign * axial**2 * slope_rhoz / 2 - slope_zz / 2
        massless = 2 / axial * (
            jet**2 - 2 * source_sign * slope**2 * phase_z - slope**2 * moment
            + slope**4 - slope * jet_z + jet * slope_z)
        extra = axial * mass**2 * (
            axial**2 * axial_rho**2 + axial**2 * mass**2 + 8 * slope**2
            - 4 * moment - 8 * source_sign * phase_z) / 8
        checks['%s/N_minus_s_J3_over_a' % tag] = _exact(
            remainder - massless - extra, tag + ' N')
        reference_rho, reference_phase, reference_z = sp.symbols('bp0 q0 qz0', real=True)
        history = remainder.subs({slope_z: 0, slope_zz: 0})
        baseline = history.subs({
            slope_rho: reference_rho, slope_rhoz: 0, phase: reference_phase,
            phase_z: reference_z})
        difference = sp.expand(history - baseline)
        claimed = (
            axial**3 / 2 * (slope_rho**2 - reference_rho**2)
            - 4 * source_sign * slope**2 * (phase_z - reference_z) / axial
            + source_sign * axial * slope * slope_rhoz
            - source_sign * axial * mass**2 * (phase_z - reference_z))
        checks['%s/Sigma_difference' % tag] = _exact(
            difference - claimed, tag + ' Sigma')
        checks['%s/Sigma_common_k_cancels' % tag] = _exact(
            sp.diff(difference, moment), tag + ' k')
        checks['%s/Sigma_q_cancels' % tag] = _exact(
            sp.diff(difference, phase) + sp.diff(difference, reference_phase),
            tag + ' q')
        radius, radial_z, radial_zz, radial_rhoz = sp.symbols(
            'r rz rzz rpz', real=True, positive=True)
        angular = sp.symbols('ell', real=True)
        matched = axial * angular / (2 * radius)
        slope_from_radius = -matched * radial_z / radius
        second_from_radius = (
            -slope_from_radius * radial_z / radius
            - matched * (radial_zz / radius - radial_z**2 / radius**2))
        checks['%s/radius_z_vanishes_sets_b_z' % tag] = _exact(
            slope_from_radius.subs(radial_z, 0), tag + ' bz')
        checks['%s/radius_zz_vanishes_sets_b_zz' % tag] = _exact(
            second_from_radius.subs({radial_z: 0, radial_zz: 0}), tag + ' bzz')
        mixed = -matched * radial_rhoz / radius
        checks['%s/b_rhoz_survives' % tag] = _exact(mixed.subs({
            axial: 1, angular: 1, radius: 1, radial_rhoz: 1}) + sp.Rational(1, 2),
            tag + ' rhoz')
    return checks


def _factor_blocks():
    """b=ell b0, q_z=ell^2 Qz, q_zz=ell^2 Qzz split J and Sigma N."""
    rho, z = sp.symbols('rho z', real=True)
    axial = sp.Function('a', positive=True)(rho)
    base = sp.Function('b0', real=True)(rho, z)
    first = sp.Function('Qz', real=True)(rho, z)
    second = sp.Function('Qzz', real=True)(rho, z)
    angular, mass = sp.symbols('ell m', real=True)
    checks = {}
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        slope = angular * base
        phase_z = angular**2 * first
        phase_zz = angular**2 * second
        geometric = (
            -axial**2 * slope * sp.diff(slope, rho, 2, z, 1)
            + axial**2 * sp.diff(slope, rho) * sp.diff(slope, rho, z)
            - 2 * axial * sp.diff(axial, rho) * slope * sp.diff(slope, rho, z)
            - source_sign * slope * sp.diff(slope, rho, 1, z, 2)
            + source_sign * sp.diff(slope, z) * sp.diff(slope, rho, z)
            - 16 * slope**3 * sp.diff(slope, z) / axial**2)
        quadratic = (
            -axial**2 * base * sp.diff(base, rho, 2, z, 1)
            + axial**2 * sp.diff(base, rho) * sp.diff(base, rho, z)
            - 2 * axial * sp.diff(axial, rho) * base * sp.diff(base, rho, z)
            - source_sign * base * sp.diff(base, rho, 1, z, 2)
            + source_sign * sp.diff(base, z) * sp.diff(base, rho, z))
        quartic = (
            -16 * base**3 * sp.diff(base, z) / axial**2
            + 4 * source_sign * base**2 * second / axial**2)
        mixed = -4 * base * sp.diff(base, z) + source_sign * second
        rhs = sp.expand(
            geometric + 4 * source_sign * slope**2 * phase_zz / axial**2
            + mass**2 * (source_sign * phase_zz - 4 * slope * sp.diff(slope, z)))
        claimed = sp.expand(
            angular**2 * quadratic + angular**4 * quartic + mass**2 * angular**2 * mixed)
        checks['%s/J_factor' % tag] = _exact(rhs - claimed, tag + ' factor J')
        checks['%s/Qz_factor' % tag] = _exact(
            -4 * slope * sp.diff(slope, z) / axial**2
            - angular**2 * (-4 * base * sp.diff(base, z) / axial**2), tag + ' Qz')
        checks['%s/Qzz_factor' % tag] = _exact(
            -4 * (sp.diff(slope, z)**2 + slope * sp.diff(slope, z, 2)) / axial**2
            - angular**2 * (
                -4 * (sp.diff(base, z)**2 + base * sp.diff(base, z, 2)) / axial**2),
            tag + ' Qzz')
        checks['%s/even_in_ell' % tag] = _exact(
            rhs.subs(angular, -angular) - rhs, tag + ' even')
        current2, current4, current_m = sp.symbols('J2 J4 JM', real=True)
        base_rho, base_rho0, base_rhoz, grad, grad0 = sp.symbols(
            'b0p b0p0 b0pz Q Q0', real=True)
        split = angular**2 * current2 + angular**4 * current4 + mass**2 * angular**2 * current_m
        surface = (
            source_sign * split / axial
            + axial**3 / 2 * angular**2 * (base_rho**2 - base_rho0**2)
            - 4 * source_sign * (angular * base)**2 * angular**2 * (grad - grad0) / axial
            + source_sign * axial * (angular * base) * angular * base_rhoz
            - source_sign * axial * mass**2 * angular**2 * (grad - grad0))
        quadratic_n = (
            source_sign * current2 / axial
            + axial**3 / 2 * (base_rho**2 - base_rho0**2)
            + source_sign * axial * base * base_rhoz)
        quartic_n = (
            source_sign * current4 / axial
            - 4 * source_sign * base**2 * (grad - grad0) / axial)
        mixed_n = source_sign * current_m / axial - source_sign * axial * (grad - grad0)
        checks['%s/N_factor' % tag] = _exact(
            sp.expand(surface - (
                angular**2 * quadratic_n + angular**4 * quartic_n
                + mass**2 * angular**2 * mixed_n)), tag + ' factor N')
    return checks


def _action_signs_at_group14_labels():
    """Same occupied-column signs at the positive group14 mass and angular labels."""
    checks = {}
    mass = float(np.pi / 2)
    angular = float(np.sqrt(5))
    for source_sign in SOURCE_SIGNS:
        columns = np.zeros((1, 2, 1), complex)
        columns[0, _major_index(source_sign), 0] = 1
        axial = -2j * source_sign * columns
        empty = np.zeros((0, *columns.shape), complex)
        result = source_column_matter(
            columns, axial, np.eye(1), empty, empty,
            mass=mass, angular=angular, axial_scale=1, radius=1, multiplicity=1)
        expected = np.array([2.0, -2.0 * source_sign], dtype=float)
        if not np.array_equal(result['action_gradient'][0], expected):
            raise ArithmeticError({
                'source_sign': source_sign,
                'gradient': result['action_gradient'][0].tolist(),
                'expected': expected.tolist(),
            })
        checks['s%+d/N' % source_sign] = float(result['action_gradient'][0, 0])
        checks['s%+d/beta' % source_sign] = float(result['action_gradient'][0, 1])
    return checks


@lru_cache(maxsize=1)
def massive_cubic_identities():
    """Exact residuals for both source signs. Every algebraic value is '0'."""
    owner = reduced_massless_transport_identities()
    if (owner['Ts_J3'] != MASSLESS_J3 or owner['Ts_qz'] != MASSLESS_QZ
            or owner['Ts_qzz'] != MASSLESS_QZZ):
        raise ValueError('massless reduced-current owner changed')
    checks = {}
    checks.update(_operator_blocks())
    checks.update(_transport_blocks())
    checks.update(_surface_blocks())
    checks.update(_factor_blocks())
    if any(value != '0' for value in checks.values()):
        raise ArithmeticError(checks)
    signs = confirm_action_signs_against_source_column_matter()
    massive_signs = _action_signs_at_group14_labels()
    return MappingProxyType({
        'residuals': checks,
        'massless_Ts_J3': MASSLESS_J3,
        'Ts_qz': MASSLESS_QZ,
        'Ts_qzz': MASSLESS_QZZ,
        'Ts_J3_massive': MASSLESS_J3 + ' + m^2 (s q_zz - 4 b b_z)',
        'Sigma_delta_N_bracket': (
            's delta(J3)/a + a^3/2 delta(b_rho^2) - 4 s b^2 delta(q_z)/a '
            '+ s a b delta(b_rhoz) - s a m^2 delta(q_z)'),
        'factorization': {
            'b': 'ell b0',
            'q_z': 'ell^2 Qz',
            'q_zz': 'ell^2 Qzz',
            'J3': 'ell^2 J2 + ell^4 J4 + m^2 ell^2 JM',
            'N_bracket': 'ell^2 N2 + ell^4 N4 + m^2 ell^2 NM',
            'Qz_prime': '-4 b0 b0_z/a^2',
            'Qzz_prime': '-4 (b0_z^2 + b0 b0_zz)/a^2',
            'J2_prime': (
                '-a^2 b0 b0_rhorhoz + a^2 b0_rho b0_rhoz - 2 a a_rho b0 b0_rhoz '
                '- s b0 b0_rhozz + s b0_z b0_rhoz'),
            'J4_prime': '-16 b0^3 b0_z/a^2 + (4 s b0^2/a^2) Qzz',
            'JM_prime': '-4 b0 b0_z + s Qzz',
            'upstream_differences': 'zero, the same homogeneous representative',
        },
        'L_h_parent': '-m ell/4 (s a^2 r_rho + r_z)/r^2',
        'massless_h_equation_imposed': False,
        'h_cancelled_from_J3_and_Re_lC': True,
        'common_k_cancels_on_Sigma': True,
        'action_N': '-mu/(2 pi) N_bracket',
        'action_beta': '+mu/(2 pi) J3',
        'action_signs': dict(signs),
        'group14_label_action_signs': massive_signs,
        'higher_uv_remainder_C_M': None,
        'uniform_C4_on_I': None,
        'complete_UV_tail': None,
        'physical_local_gate': 'OPEN',
        'one_point_or_formal_coverage_is_full_source_tail': False,
        'factorization_replaces_a_finished_family_integral': False,
    })


class MassiveCubicGeometry(CubicGeometry):
    """Massive forcing on the unchanged massless jets.

    The historical integrator accepts this subclass. Current forcing gains
    m^2 a(rho)^2 times the q_z forcing, and the q_zz coupling gains s m^2.
    Those are the proved additions. h is not set to zero here.
    """

    def __init__(self, family, angular, mass, *, bits=160):
        super().__init__(family, angular, bits=bits)
        if isinstance(mass, (bool, np.bool_)):
            raise ValueError('explicit real mass required')
        with ctx.workprec(self.bits):
            value = arb(mass) if isinstance(mass, arb) else G._binary_arb(mass, 'mass')
            if not value.is_finite() or not value >= 0:
                raise ValueError('finite nonnegative mass required')
            self.mass = value

    def _mass_square(self):
        return self.mass ** 2

    def forcing_series(self, rho, z_in, sign, order):
        series = CubicGeometry.forcing_series(self, rho, z_in, sign, order)
        with ctx.workprec(self.bits), G._JetWork(order):
            axial = G.background_series(G._rho_ball(rho), order).a
            mass2 = self._mass_square()
            current = series['current'] + mass2 * axial * axial * series['qz']
            coupling = series['coupling'] + sign * mass2
            result = {
                'qz': series['qz'], 'qzz': series['qzz'],
                'current': current, 'coupling': coupling}
            if any(not coefficient.is_finite()
                   for item in result.values() for coefficient in G._coeffs(item, order)):
                raise ArithmeticError('nonfinite massive characteristic forcing jet')
            return result

    def forcing(self, rho, z_in, sign):
        base = CubicGeometry.forcing(self, rho, z_in, sign)
        with ctx.workprec(self.bits):
            axial = G.background_series(G._rho_ball(rho), 0).a[0]
            mass2 = self._mass_square()
            values = CubicForcing(
                base.qz, base.qzz,
                base.current + mass2 * axial**2 * base.qz,
                base.coupling + sign * mass2)
            if not all(item.is_finite() for item in vars(values).values()):
                raise ArithmeticError('nonfinite massive forcing enclosure')
            return values


def enclose_massive_characteristic(model, z_in, sign, rho_up, *,
                                   cells=128, max_depth=16, order=0):
    """Reuse the historical integrator, then add the proved surface mass term.

    The correction is -s a m^2 q_z at rho=1. It is not folded into the
    historical N formula, and a point enclosure is not an interval certificate.
    """
    if not isinstance(model, MassiveCubicGeometry):
        raise TypeError('massive cubic geometry subclass required')
    result = enclose_characteristic(
        model, z_in, sign, rho_up, cells=cells, max_depth=max_depth, order=order)
    with ctx.workprec(model.bits):
        surface = model.jets(arb(1), result['z_domain'])
        correction = -sign * surface['a'] * model._mass_square() * result['q_z']
        bracket = result['N_bracket3'] + correction
        if not bracket.is_finite():
            raise ArithmeticError('nonfinite massive surface coefficient')
        out = dict(result)
        out['N_bracket3'] = bracket
        out['surface_mass_correction'] = correction
        out['mass'] = model.mass
        out['coefficient_scope'] = (
            'history-minus-reference massive cubic coefficient; '
            'upstream homogeneous differences zero; h retained until cancellation')
        out['higher_uv_remainder_bound'] = None
        out['C_M'] = None
        out['complete_UV_tail'] = None
        out['entire_incoming_interval'] = False
        out['physical_local_gate'] = 'OPEN'
        return out


def original_group14_labels(root=None):
    """Positive group14 channel: mass pi/2, |ell|=sqrt(5), ledger multiplicity."""
    root = Path(root) if root is not None else _lab_root()
    cauchy = json.loads((root / CAUCHY_RECORD).read_text())
    channel = next(row for row in cauchy['channels'] if row['index'] == GROUP14)
    bridge = json.loads((root / CUTOFF_BRIDGE_RECORD).read_text())
    ledger = [row for row in bridge['ledger'] if row['group'] == GROUP14]
    inventory = json.loads((root / INVENTORY_RECORD).read_text())
    measures = inventory['quadrature_measure_by_family']
    if channel['compact_mass'] <= 0 or channel['angular_eigenvalue'] <= 0:
        raise ValueError('group14 must stay the positive massive angular channel')
    if len(ledger) != 2 or {row['positive_angular_sign'] for row in ledger} != {-1, 1}:
        raise ValueError('group14 angular pair missing from the cutoff ledger')
    multiplicity = ledger[0]['multiplicity_per_signed_family']
    if any(row['multiplicity_per_signed_family'] != multiplicity for row in ledger):
        raise ValueError('group14 angular partners have different multiplicities')
    if abs(float(channel['angular_eigenvalue']) - float(np.sqrt(5))) > 1e-12:
        raise ValueError('group14 angular label is not sqrt(5)')
    if abs(float(channel['compact_mass']) - float(np.pi / 2)) > 1e-12:
        raise ValueError('group14 mass label is not pi/2')
    if int(multiplicity) != 12 or int(channel['degeneracy']) != 12:
        raise ValueError('group14 multiplicity or degeneracy is not 12')
    cutoffs = {sign: float(measures['14_%d' % sign]) for sign in SOURCE_SIGNS}
    if set(cutoffs.values()) != {160.0}:
        raise ValueError('group14 quadrature measure is not the recorded 160 cutoff')
    for row in ledger:
        if abs(float(row['absolute_angular']) - float(np.sqrt(5))) > 1e-12:
            raise ValueError('group14 ledger angular label drifted')
        if abs(float(row['angular_square_weight_per_energy_sign']) - 60.0) > 1e-9:
            raise ValueError('group14 angular-square weight is not 12*5')
    with ctx.workprec(160):
        mass = arb(channel['compact_mass']).union(arb.pi() / 2)
        angular = arb(channel['angular_eigenvalue']).union(arb(5).sqrt())
        if mass.rad() > arb('1e-15') or angular.rad() > arb('1e-15'):
            raise ValueError('metadata is not the binary image of pi/2 and sqrt(5)')
    return MappingProxyType({
        'group': GROUP14,
        'mass': float(channel['compact_mass']),
        'absolute_angular': float(channel['angular_eigenvalue']),
        'multiplicity_per_signed_family': float(multiplicity),
        'degeneracy': int(channel['degeneracy']),
        'compact_level': int(channel['compact_level']),
        'angular_level': int(channel['angular_level']),
        'quadrature_measure_by_sign': cutoffs,
        'cutoff': 160.0,
        'archived_splits': (160.0, 320.0),
        'positive_mass': True,
        'mass_ball': mass,
        'angular_ball': angular,
    })


def require_original_history(root=None):
    root = Path(root) if root is not None else _lab_root()
    history = json.loads((root / HISTORY_RECORD).read_text())
    family = LocalIncomingFamily(np.asarray(history['history']['coefficients']))
    identity = profile_identity(family)
    if (identity != CURRENT_HISTORY_IDENTITY
            or not identity.startswith(CURRENT_HISTORY_IDENTITY_PREFIX)):
        raise ValueError('original current history 0b0e4ced required')
    return family, identity
