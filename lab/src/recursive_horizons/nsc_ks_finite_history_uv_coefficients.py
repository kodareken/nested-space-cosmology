"""Finite-history vacuum-envelope coefficients from the owned e-recurrence.

The source-cutoff bridge recurrence is executed on the actual Dirac operator
L0, not on dummy symbols for L0 A_j. Coefficients are the mathematical
vacuum envelope of the unchanged upstream affine expansion. They do not
replace physical source columns, bound C_M, or close the local incoming gate.

Constructed here
----------------
A0, minor A1, major-A1 transport, minor A2, major-A2 transport, and the A3
transport/recurrence formulae needed to bookkeep an E^{-2} coincidence
current. Differentiating the source carrier lowers an e^{-3} two-point
coefficient to e^{-2}; h_{s,z}/e^2 is therefore not the full current
coefficient. Minor A3 is not a local metric jet: it depends on the
integrated major A2.

Upstream
--------
The normalized affine-vacuum column at rho_up has real positive major A0 and
vanishing major A1. Transport then generates the phase. Reference and
current histories share that upstream value because g equals the reference
on an upstream neighborhood.
"""
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from types import MappingProxyType
import json

import numpy as np
import sympy as sp

from .nsc_compatible_history_geometry import (
    CompatibleIncomingMetric,
    CompatibleRadiusDirection,
)
from .nsc_dirac_source_phase import (
    DEFAULT_GAUSS_NODES,
    RHO_SIGMA,
    SOURCE_SIGNS,
    _characteristic_distance,
    _gauss_nodes,
    _normal_coordinate,
    _panels,
)
from .nsc_ks_profile_identity import profile_identity
from .nsc_ks_source_envelope import RHO_UP_MIN
from .nsc_local_incoming_family import LocalIncomingFamily
from .nsc_pg_ks_metric_pullback import reference_chart
from .nsc_transmitting_dirac_domain import S1, S2, S3


ARCHIVE_RHO_UP = 1.0300000000000002
CURRENT_HISTORY_IDENTITY_PREFIX = '0b0e4ced'
PHYSICAL_LOCAL_GATE = (
    'OPEN: leading paired e^{-2} vacuum coefficient identified under original '
    'equal-mu same-cutoff angular pairs; C_M, C4, thermal remainder, '
    'upstream error and the local gate remain unbounded')
INVENTORY_RECORD = 'results/development/nsc-ks-source-inventory.json'
CUTOFF_BRIDGE_CONTROL = 'results/development/nsc-ks-cutoff-bridge-control.json'
V2_RECORD = 'results/development/nsc-ks-current-uv-coefficients-v2.json'


def _major_index(source_sign):
    if source_sign not in SOURCE_SIGNS:
        raise ValueError('source sign must be +1 or -1')
    return 0 if source_sign == 1 else 1


def _minor_index(source_sign):
    return 1 - _major_index(source_sign)


def occupied_basis(source_sign):
    """A0 = e_s, the occupied S3=s unit spinor."""
    vector = np.zeros(2)
    vector[_major_index(source_sign)] = 1.0
    return vector


def potential_matrix(mass, angular, radius):
    """V = -m S1 + ell S2 / r, with the original source signs."""
    mass = float(mass)
    angular = float(angular)
    radius = float(radius)
    if not np.isfinite([mass, angular, radius]).all() or radius <= 0:
        raise ValueError('finite mass, angular eigenvalue and positive radius required')
    return -mass * S1 + (angular / radius) * S2


def _sym_pauli():
    s1 = sp.Matrix([[0, 1], [1, 0]])
    s2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    s3 = sp.diag(1, -1)
    return s1, s2, s3


def _sym_potential(mass, angular, radius):
    s1, s2, _ = _sym_pauli()
    return -mass * s1 + angular / radius * s2


def _apply_L0(psi, axial_scale, potential, s3, rho, z):
    """L0 = ∂_ρ - S3/a² ∂_z - i V/a, the owned envelope operator after the e-phase."""
    return psi.diff(rho) - (s3 / axial_scale**2) * psi.diff(z) - sp.I * potential / axial_scale * psi


def _T(scalar, source_sign, axial_scale, rho, z):
    """T_s = ∂_ρ - s/a² ∂_z, the incoming occupied characteristic derivative."""
    return scalar.diff(rho) - source_sign / axial_scale**2 * scalar.diff(z)


def _real_residual(expression):
    """Treat declared-real geometry jets as real; SymPy does not infer that."""
    real_jets = {
        jet: sp.Symbol('real_jet_%d' % index, real=True)
        for index, jet in enumerate(sorted(expression.atoms(sp.Derivative), key=str))
    }
    residual = sp.simplify(sp.expand_complex(expression.xreplace(real_jets)))
    return str(residual)


def _projector(source_sign, s3):
    identity = sp.eye(2)
    occupied = (identity + source_sign * s3) / 2
    return occupied, identity - occupied


def _unit(source_sign):
    vector = sp.zeros(2, 1)
    vector[_major_index(source_sign)] = 1
    return vector


def _minor_A1(axial_scale, mass, angular, radius, source_sign):
    """Algebraic minor coefficient, including the original source sign in ell."""
    return axial_scale / 2 * (mass - sp.I * source_sign * angular / radius)


@lru_cache(maxsize=1)
def coefficient_identities():
    """Independent L0 substitution of the constructed jets, both source signs.

    Residues are exact symbolic zeros. They are not a numerical C_M and not
    a finite-history cancellation of the assembled N,beta current.
    """
    rho, z = sp.symbols('rho z', real=True)
    mass, angular = sp.symbols('m ell', real=True)
    axial = sp.Function('a', real=True)(rho)
    radius = sp.Function('r', real=True)(rho, z)
    phase = sp.Function('psi', real=True)(rho, z)
    major2 = sp.Function('Phi2')(rho, z)
    s1, s2, s3 = _sym_pauli()
    potential = _sym_potential(mass, angular, radius)
    square = mass**2 + angular**2 / radius**2
    identities = {}
    potential_residual = sp.simplify(potential - (-mass * s1 + angular / radius * s2))
    identities['V_definition'] = [str(value) for value in potential_residual]
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        major = _major_index(source_sign)
        minor = _minor_index(source_sign)
        occupied, vacant = _projector(source_sign, s3)
        a0 = _unit(source_sign)
        mu = _minor_A1(axial, mass, angular, radius, source_sign)
        a1 = sp.zeros(2, 1)
        a1[major] = sp.I * phase
        a1[minor] = mu
        l0_a0 = _apply_L0(a0, axial, potential, s3, rho, z)
        identities['%s/Pi_minus_A0' % tag] = [str(sp.simplify(value)) for value in vacant * a0]
        identities['%s/A0_is_e_s' % tag] = [
            str(sp.simplify(a0[major] - 1)), str(sp.simplify(a0[minor]))]
        identities['%s/minor_A1_from_L0_A0' % tag] = str(sp.simplify(
            a1[minor] - axial**2 / (2 * sp.I) * (vacant * l0_a0)[minor]))
        identities['%s/minor_A1_parent' % tag] = str(sp.simplify(
            a1[minor] - axial / 2 * (mass - sp.I * source_sign * angular / radius)))
        # P1 of the E>0 projector expansion matches the s=+1 minor coefficient.
        if source_sign == 1:
            p1 = axial * (mass * s1 - angular / radius * s2) / 2
            identities['%s/minor_A1_matches_P1' % tag] = str(sp.simplify(
                a1[minor] - p1[minor, major]))
        l0_a1 = _apply_L0(a1, axial, potential, s3, rho, z)
        transport_a1 = _T(a1[major], source_sign, axial, rho, z)
        # Π_s L0 A1 = T_s major A1 - (i/a) V_{s,-s} minor A1.
        identities['%s/Pi_s_L0_A1_is_transport_defect' % tag] = _real_residual(
            (occupied * l0_a1)[major]
            - (transport_a1 - sp.I / axial * potential[major, minor] * a1[minor]))
        identities['%s/T_s_major_A1_parent' % tag] = _real_residual(
            (sp.I / axial) * potential[major, minor] * a1[minor] + sp.I / 2 * square)
        # Impose the major-A1 transport and the recurrence for minor A2.
        a1_rho = source_sign / axial**2 * a1[major].diff(z) - sp.I / 2 * square
        l0_a1_on_shell = l0_a1.subs(a1[major].diff(rho), a1_rho)
        nu = axial**2 / (2 * sp.I) * (vacant * l0_a1_on_shell)[minor]
        a2 = sp.zeros(2, 1)
        a2[major] = major2
        a2[minor] = nu
        identities['%s/minor_A2_recurrence' % tag] = _real_residual(
            a2[minor] - axial**2 / (2 * sp.I) * (vacant * l0_a1_on_shell)[minor])
        l0_a2 = _apply_L0(a2, axial, potential, s3, rho, z)
        transport_a2 = _T(a2[major], source_sign, axial, rho, z)
        identities['%s/Pi_s_L0_A2_is_transport_defect' % tag] = _real_residual(
            (occupied * l0_a2)[major]
            - (transport_a2 - sp.I / axial * potential[major, minor] * a2[minor]))
        parent_imag = -mass * angular / 4 * (
            source_sign * axial**2 * radius.diff(rho) + radius.diff(z)) / radius**2
        identities['%s/Im_T_s_major_A2_parent' % tag] = _real_residual(
            sp.im(sp.expand_complex(
                (sp.I / axial) * potential[major, minor] * a2[minor])) - parent_imag)
        # Finite constructed telescope: (L0 - 2 i e / a² Π_{-s}) sum_{j=0}^2 e^{-j} A_j
        # equals e^{-2} L0 A2 once the recurrence and major-A1 transport hold.
        energy = sp.symbols('e', positive=True)
        envelope = a0 + a1 / energy + a2 / energy**2
        exact = _apply_L0(envelope, axial, potential, s3, rho, z) - (
            2 * sp.I * energy / axial**2) * (vacant * envelope)
        a2_rho = source_sign / axial**2 * a2[major].diff(z) + (
            sp.I / axial) * potential[major, minor] * a2[minor]
        exact_on_shell = exact.subs({
            a1[major].diff(rho): a1_rho,
            a2[major].diff(rho): a2_rho,
        })
        defect = exact_on_shell - _apply_L0(a2, axial, potential, s3, rho, z).subs(
            a2[major].diff(rho), a2_rho) / energy**2
        identities['%s/L0_telescope_M2' % tag] = [
            _real_residual(component) for component in defect]
        l0_a2_on_shell = l0_a2.subs(a2[major].diff(rho), a2_rho)
        nu3 = axial**2 / (2 * sp.I) * (vacant * l0_a2_on_shell)[minor]
        identities['%s/minor_A3_depends_on_major_A2' % tag] = (
            'depends' if major2 in nu3.atoms(sp.Function) else 'independent')
        # Homogeneous / free controls on the constructed Im(T_s major A2).
        imag_rhs = sp.im(sp.expand_complex(
            (sp.I / axial) * potential[major, minor] * a2[minor]))
        identities['%s/Im_T_s_major_A2_vanishes_for_m_0' % tag] = _real_residual(
            imag_rhs.subs(mass, 0))
        identities['%s/Im_T_s_major_A2_vanishes_for_ell_0' % tag] = _real_residual(
            imag_rhs.subs(angular, 0))
        identities['%s/Im_T_s_major_A2_odd_in_ell' % tag] = _real_residual(
            imag_rhs + imag_rhs.subs(angular, -angular))
        # Normalization: Re(major A1)=0, A0†A1 + A1†A0 = 0.
        identities['%s/Re_major_A1_normalization' % tag] = str(sp.simplify(
            sp.re(sp.expand_complex(a1[major]))))
        identities['%s/order_e_minus1_density' % tag] = str(sp.simplify(sp.expand_complex(
            (a0.H * a1 + a1.H * a0)[0])))
    residual_keys = [key for key in identities if not key.endswith('depends_on_major_A2')]
    leaves = []
    for key in residual_keys:
        value = identities[key]
        if isinstance(value, list):
            leaves.extend(value)
        else:
            leaves.append(value)
    if set(leaves) != {'0'}:
        raise ArithmeticError(identities)
    if any(identities[key] != 'depends' for key in identities if key.endswith('depends_on_major_A2')):
        raise ArithmeticError('minor A3 must retain the integrated major A2')
    return MappingProxyType({
        'residuals': identities,
        'V': 'V=-m S1 + ell S2 / r',
        'T_s': 'T_s = d_rho - s/a^2 d_z',
        'A0': 'e_s',
        'minor_A1': 'a/2 (m - i s ell / r)',
        'T_s_major_A1': '-i/2 (m^2 + ell^2/r^2)',
        'Im_T_s_major_A2_pure_imag_major_A1': (
            '-m ell/4 (s a^2 r_rho + r_z)/r^2'),
        'upstream': (
            'normalized affine vacuum at rho_up: A0=e_s real positive, '
            'major A1=0; g equals the reference on an upstream neighborhood'),
        'telescope': (
            '(L0 - 2 i e / a^2 Pi_{-s}) sum_{j=0}^2 e^{-j} A_j = e^{-2} L0 A2 '
            'on the constructed recurrence and major transports; A3 is the '
            'same recurrence applied to integrated A2, not a local jet'),
        'dummy_symbol_L0_A_j': False,
        'physical_source_columns_replaced': False,
        'outgoing_Weyl_identified_with_source_E': False,
    })


@lru_cache(maxsize=1)
def e_minus2_current_bookkeeping():
    """Distinguish h_z/E^2 from the actual N and beta coincidence coefficients.

    The owned vertices are S3/a for N-momentum and -I for beta-momentum, with
    the action minus sign, and V for the N density vertex. Carrier
    differentiation of an e^{-3} two-point coefficient feeds e^{-2}.
    """
    rho, z = sp.symbols('rho z', real=True)
    mass, angular = sp.symbols('m ell', real=True)
    axial = sp.Function('a', real=True)(rho)
    radius = sp.Function('r', real=True)(rho, z)
    phase = sp.Function('psi', real=True)(rho, z)
    energy = sp.symbols('e', positive=True)
    terms = {}
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        major = _major_index(source_sign)
        minor = _minor_index(source_sign)
        mu = _minor_A1(axial, mass, angular, radius, source_sign)
        a0 = _unit(source_sign)
        a1 = sp.zeros(2, 1)
        a1[major] = sp.I * phase
        a1[minor] = mu
        a2 = sp.Matrix(sp.symbols('A2_0:2', complex=True))
        a3 = sp.Matrix(sp.symbols('A3_0:2', complex=True))
        a2z = sp.Matrix(sp.symbols('A2z_0:2', complex=True))
        a1z = sp.zeros(2, 1)
        a1z[major] = sp.I * phase.diff(z)
        a1z[minor] = mu.diff(z)
        a0z = sp.zeros(2, 1)
        # J_{(2)} = Im(A0† A2_z + A1† A1_z + A2† A0_z) - s n3, A0_z=0.
        interference = (a1.H * a1z)[0]
        h_z = sp.im(sp.expand_complex(a2z[major]))
        n3 = (a0.H * a3 + a1.H * a2 + a2.H * a1 + a3.H * a0)[0]
        beta_e2 = h_z + sp.im(sp.expand_complex(interference)) - source_sign * n3
        interference_closed = _real_residual(
            sp.im(sp.expand_complex(interference))
            - axial**2 * mass * source_sign * angular * radius.diff(z) / (4 * radius**2))
        density = a0 * a0.H + (a0 * a1.H + a1 * a0.H) / energy + (
            a0 * a2.H + a1 * a1.H + a2 * a0.H) / energy**2
        potential = _sym_potential(mass, angular, radius)
        n_density = sp.trace(potential * density)
        terms[tag] = {
            'beta_e_minus2_contains_h_z': True,
            'beta_e_minus2_contains_A1_interference': True,
            'beta_e_minus2_contains_n3_hence_A3': True,
            'Im_A1dag_A1_z_residual': interference_closed,
            'N_density_vertex': 'V=-m S1 + ell S2 / r; tr(V rho) starts at e^{-1}',
            'N_momentum_vertex': 'S3/a',
            'beta_momentum_vertex': '-I, action supplies the overall minus',
            'carrier_lowers_C3': (
                'd/d eta of exp(-i s e eta) at eta=0 multiplies the coincidence '
                'e^{-3} two-point coefficient and feeds e^{-2}'),
            'h_z_is_the_full_e_minus2_coefficient': False,
            'A3_locally_algebraic_in_metric_jets': False,
            'symbolic_beta_e_minus2': str(sp.simplify(beta_e2)),
            'N_density_e_minus2_uses_A2': 'A2' in str(n_density),
        }
        if interference_closed != '0':
            raise ArithmeticError(terms)
    return MappingProxyType({
        'source_signs': SOURCE_SIGNS,
        'terms': terms,
        'linear_fixed_transfer_e_minus2_zero': (
            'owned separately at the reference for linear response; not a '
            'finite-history coincidence identity'),
        'angular_pair_of_Im_T_s_major_A2': (
            'that one constructed integrand is odd in ell; this is not a '
            'finite-history UV cancellation of the full N,beta coefficient'),
        'full_versus_difference': (
            'Im(T_s major A2) on the reference is z-independent, so its '
            'incoming z-derivative coincides with the g-minus-reference '
            'value; n3 follows that z-derivative after continuity+upstream IC'),
    })


@lru_cache(maxsize=1)
def e_minus2_sigma_contraction():
    """Leading Sigma E^{-2} current from continuity and recurrence, without solving A3.

    Occupied vacuum envelope only. Thermal occupations remain the exponential
    remainder. A3 enters n3, but n3 is replaced by 2 s h_z after T_s transport
    and a z-homogeneous normalized upstream condition.
    """
    rho, z = sp.symbols('rho z', real=True)
    mass, angular = sp.symbols('m ell', real=True)
    axial = sp.Function('a', real=True)(rho)
    radius = sp.Function('r', real=True)(rho, z)
    phase = sp.Function('q', real=True)(rho, z)
    height = sp.Function('h', real=True)(rho, z)
    n2 = sp.Function('n2', real=True)(rho, z)
    n3 = sp.Function('n3', real=True)(rho, z)
    minor1 = sp.Function('L')(rho, z)
    minor2 = sp.Function('B')(rho, z)
    energy = sp.Symbol('e', positive=True)
    identities = {}
    potential = _sym_potential(mass, angular, radius)
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        major = _major_index(source_sign)
        minor = _minor_index(source_sign)
        occupied, vacant = _projector(source_sign, sp.diag(1, -1))
        l = _minor_A1(axial, mass, angular, radius, source_sign)
        a1 = sp.zeros(2, 1)
        a1[major] = sp.I * phase
        a1[minor] = l
        square = mass**2 + angular**2 / radius**2
        a1_rho = source_sign / axial**2 * a1[major].diff(z) - sp.I / 2 * square
        l0_a1 = _apply_L0(a1, axial, potential, sp.diag(1, -1), rho, z)
        b = axial**2 / (2 * sp.I) * (vacant * l0_a1.subs(a1[major].diff(rho), a1_rho))[minor]
        re_lb = sp.re(sp.expand_complex(sp.conjugate(l) * b))
        t_minus_l = _T(l, -source_sign, axial, rho, z)
        identities['%s/Re_conj_l_b_from_Tminus_l' % tag] = _real_residual(
            re_lb - axial**2 / 2 * sp.im(sp.expand_complex(sp.conjugate(l) * t_minus_l)))
        ts_h = sp.im(sp.expand_complex(
            (sp.I / axial) * potential[major, minor] * b))
        identities['%s/Re_conj_l_b_equals_minus_a2_over_2_Ts_h' % tag] = _real_residual(
            re_lb + axial**2 / 2 * ts_h)
        reLB = sp.re(sp.conjugate(minor1) * minor2)
        minor_sq = sp.conjugate(minor1) * minor1 / energy**2 + 2 * reLB / energy**3
        density = 1 + n2 / energy**2 + n3 / energy**3
        current = source_sign * (density - 2 * minor_sq)
        continuity = density.diff(rho) - current.diff(z) / axial**2
        order3 = sp.expand(continuity * energy**3).subs(energy, 0)
        ts_n3 = _T(n3, source_sign, axial, rho, z)
        identities['%s/continuity_e_minus3_is_Ts_n3_plus_4s_a2_dz_Re' % tag] = _real_residual(
            order3 - (ts_n3 + 4 * source_sign / axial**2 * reLB.diff(z)))
        identities['%s/Ts_n3_minus_2s_Ts_hz_after_recurrence' % tag] = _real_residual(
            -(4 * source_sign / axial**2) * (-axial**2 / 2 * ts_h).diff(z)
            - 2 * source_sign * ts_h.diff(z))
        sigma = {radius.diff(z): 0}
        im_llz = sp.im(sp.expand_complex(sp.conjugate(l) * l.diff(z)))
        identities['%s/Sigma_Im_conj_l_lz' % tag] = _real_residual(im_llz.xreplace(sigma))
        n3_on_shell = 2 * source_sign * height.diff(z)
        i_coeff = height.diff(z) + im_llz.xreplace(sigma) - source_sign * n3_on_shell
        identities['%s/Sigma_I_momentum_is_minus_hz' % tag] = _real_residual(
            i_coeff + height.diff(z))
        s3_coeff = (
            source_sign * height.diff(z) - n3_on_shell + 4 * re_lb)
        identities['%s/Sigma_S3_momentum' % tag] = _real_residual(
            s3_coeff + source_sign * height.diff(z) - 4 * re_lb)
        alpha = potential[major, minor]
        beta_v = potential[minor, major]
        identities['%s/N_density_q_channel_vanishes' % tag] = _real_residual(
            beta_v * sp.conjugate(l) - alpha * l)
        tr2 = sp.re(sp.expand_complex(alpha * b + beta_v * sp.conjugate(b)))
        identities['%s/trV_rho2_plus_4_a_Re_conj_l_b' % tag] = _real_residual(
            tr2 + 4 / axial * re_lb)
        n_bracket = tr2 + s3_coeff / axial
        identities['%s/N_e2_bracket_is_minus_s_hz_over_a' % tag] = _real_residual(
            n_bracket + source_sign * height.diff(z) / axial)
        delta_r_rho = sp.Symbol('delta_r_rho', real=True)
        offdiag = source_sign * axial**3 * angular * delta_r_rho / (4 * radius**2)
        identities['%s/local_delta_Re_conj_l_b' % tag] = _real_residual(
            sp.re(sp.expand_complex(sp.conjugate(l) * offdiag)) - axial * mass * offdiag / 2)
        identities['%s/Im_Ts_major_A2_odd_ell' % tag] = _real_residual(
            ts_h + ts_h.subs(angular, -angular))
    leaves = []
    for value in identities.values():
        leaves.append(value)
    if set(leaves) != {'0'}:
        raise ArithmeticError(identities)
    vertices = confirm_raw_ks_vertices()
    pairs = original_angular_pair_ledger()
    return MappingProxyType({
        'residuals': identities,
        'n3': '2 s h_z after T_s transport and normalized z-homogeneous upstream',
        'Sigma_I_momentum_e_minus2': '-h_z',
        'Sigma_S3_momentum_e_minus2': '-s h_z + 4 Re(conj(l) b)',
        'action_N2': 'mu s delta(h_z) / (2 pi a)',
        'action_beta2': '-mu delta(h_z) / (2 pi)',
        'local_delta_r_rho_cancels_in_N': True,
        'q_channel_in_N_density_vanishes': True,
        'A3_solved_componentwise': False,
        'n3_eliminated_by_continuity_and_upstream': True,
        'occupied_vacuum_envelope_only': True,
        'thermal_occupation_replaced_by_one': False,
        'vertices': vertices,
        'angular_pairs': pairs,
        'leading_paired_one_over_Lambda_cancels_for_equal_mu': True,
        'C4_or_C_M_bounded': False,
    })


def confirm_raw_ks_vertices():
    """Confirm N,beta vertices against the owned raw_ks_vertex_coefficients."""
    from .nsc_common_subtracted_ks_source import raw_ks_vertex_coefficients
    from .nsc_transmitting_dirac_domain import I2
    mass, angular, axial, radius = 1.25, -2.5, 0.75, 1.5
    verts = raw_ks_vertex_coefficients(
        np.array([1.0, 0.0, axial, radius]), mass, angular, envelopes=np.ones(4))
    expected_V = -mass * S1 + angular / radius * S2
    expected_N_mom = S3 / axial
    expected_beta_mom = -I2
    residuals = {
        'N_density': float(np.max(np.abs(verts.multiplication[0] - expected_V))),
        'beta_density': float(np.max(np.abs(verts.multiplication[1]))),
        'N_momentum': float(np.max(np.abs(verts.momentum[0] - expected_N_mom))),
        'beta_momentum': float(np.max(np.abs(verts.momentum[1] - expected_beta_mom))),
    }
    if any(value > 1e-15 for value in residuals.values()):
        raise ArithmeticError(residuals)
    return MappingProxyType({
        'N_density': 'V=-m S1 + ell S2 / r',
        'N_momentum': 'S3/a',
        'beta_density': 0,
        'beta_momentum': '-I',
        'numerical_residuals': residuals,
        'action_sign': 'minus once, as in source_column_matter',
        'measure': 'de/(2 pi) as in the source-cutoff edge',
    })


def original_angular_pair_ledger():
    """Equal-mu, equal-grid original angular pairs from the owned inventory."""
    root = Path(__file__).resolve().parents[2]
    inventory = json.loads((root / INVENTORY_RECORD).read_text())
    joined = json.loads((root / CUTOFF_BRIDGE_CONTROL).read_text())
    payload = inventory['payload']
    if sha256((root / payload['path']).read_bytes()).hexdigest() != payload['sha256']:
        raise ValueError('original source inventory changed')
    with np.load(root / payload['path'], allow_pickle=False) as handle:
        meta = json.loads(handle['metadata_json'].tobytes())
        rows = {}
        for name, panel in meta['panels'].items():
            if panel['angular_magnitude'] == 0:
                continue
            key = (panel['group'], panel['angular_sign'])
            rows.setdefault(key, []).append(
                np.column_stack((handle[name + '/energies'], handle[name + '/weights'])))
    groups = sorted({key[0] for key in rows})
    checks = []
    for group in groups:
        pairs = []
        for sign in (1, -1):
            pair = np.concatenate(rows[(group, sign)])
            pairs.append(pair[np.lexsort((pair[:, 1], pair[:, 0]))])
        if not np.array_equal(pairs[0], pairs[1]):
            raise ValueError('angular partners must retain the same source labels and quadrature weights')
        ledger = [row for row in joined['ledger'] if row['group'] == group]
        if (len(ledger) != 2 or {r['positive_angular_sign'] for r in ledger} != {-1, 1}
                or ledger[0]['multiplicity_per_signed_family'] != ledger[1]['multiplicity_per_signed_family']
                or ledger[0]['absolute_angular'] != ledger[1]['absolute_angular']):
            raise ValueError('equal original angular-pair multiplicity required')
        checks.append({
            'group': group,
            'positive_rows_per_angular_sign': len(pairs[0]),
            'source_grid_and_weight_residual': 0.0,
            'multiplicity_per_signed_family': ledger[0]['multiplicity_per_signed_family'],
            'same_cutoff_grid': True,
        })
    if len(checks) != 30:
        raise ValueError('all 30 original nonzero-angular pairs required')
    return tuple(checks)


def paired_leading_e_minus2_weight(h_plus, h_minus, mu_plus, mu_minus):
    """Weighted pair of odd-ell h_z coefficients. Equal mu cancels; broken mu does not."""
    h_plus = np.asarray(h_plus, float)
    h_minus = np.asarray(h_minus, float)
    mu_plus = float(mu_plus)
    mu_minus = float(mu_minus)
    if h_plus.shape != h_minus.shape or not np.isfinite(h_plus).all() or not np.isfinite(h_minus).all():
        raise ValueError('finite matching h_z samples required')
    if not np.isfinite([mu_plus, mu_minus]).all() or min(mu_plus, mu_minus) <= 0:
        raise ValueError('positive pair multiplicities required')
    return mu_plus * h_plus + mu_minus * h_minus


def vacuum_upstream_conditions():
    """Unchanged normalized affine-vacuum values used as characteristic data."""
    return MappingProxyType({
        'A0': 'e_s, real positive major',
        'major_A1': 0,
        'minor_A1': 'a/2 (m - i s ell / r) at the reference radius',
        'Re_major_A1': 0,
        'normalization': (
            'positive real major (1+|v|^2)^{-1/2} at rho_up; phase is then '
            'fixed by T_s major A1 = -i/2 Q along the incoming characteristic'),
        'history_flat_at_upstream': True,
        'physical_source_renormalized': False,
        'P1_match': 's=+1 minor A1 equals the off-diagonal E>0 projector coefficient P1',
        's_minus_uses_e_recurrence_not_positive_energy_P0': True,
    })


def named_gaps():
    return MappingProxyType({
        'numerical_C_M': None,
        'thermal_coherent_source_remainder': None,
        'e_minus3_paired_C4': None,
        'upstream_normalization_error': None,
        'source_cutoff_tail_at_160_320': None,
        'higher_A_j_beyond_n3': None,
        'leading_paired_e_minus2_vacuum_coefficient': (
            'cancels for original equal-mu same-cutoff angular pairs after '
            'n3=2 s h_z and N local delta(r_rho) cancellation; not a C_M bound'),
        'finite_history_UV_cancellation': (
            'only the leading 1/Lambda algebraic coincidence coefficient; '
            'C4, thermal remainder and upstream error stay unbounded'),
        'physical_local_gate': 'OPEN',
        'complete_e_minus2_N_beta_coefficient': None,
    })


def _plateau_jet(normal, inner, outer):
    """Value and s-derivative of the owned even C-infinity normal window."""
    normal = float(normal)
    inner = float(inner)
    outer = float(outer)
    if not 0 < inner < outer:
        raise ValueError('0 < inner < outer required')
    distance = abs(normal)
    if distance <= inner:
        return 1.0, 0.0
    if distance >= outer:
        return 0.0, 0.0
    width = outer - inner
    toward = distance - inner
    away = outer - distance
    ratio = width / away - width / toward
    window = float(np.exp(-np.logaddexp(0.0, ratio)))
    d_ratio = width / (away * away) + width / (toward * toward)
    d_window = -window * (1.0 - window) * d_ratio
    sign = 1.0 if normal >= 0.0 else -1.0
    return window, d_window * sign


def _as_one_direction(family):
    if isinstance(family, LocalIncomingFamily):
        w_function, u_function = family.functions
        return CompatibleIncomingMetric(
            (1.0,),
            (CompatibleRadiusDirection(
                w_function, u_function, family.normal_inner, family.normal_outer),))
    if isinstance(family, CompatibleIncomingMetric):
        return family
    raise TypeError('LocalIncomingFamily or CompatibleIncomingMetric required')


def radius_partial_jets(family, rho, z):
    """Actual r and the partials (r_z, r_rho, r_zz, r_rhoz) at fixed z.

    These are coordinate partials of the owned radius history. They are not
    characteristic derivatives and not a field.
    """
    metric = _as_one_direction(family)
    rho = float(rho)
    z = np.asarray(z, dtype=float)
    if z.ndim != 1 or z.size == 0 or not np.isfinite(z).all():
        raise ValueError('nonempty finite one-dimensional z required')
    if not np.isfinite(rho):
        raise ValueError('finite rho required')
    _, _, axial, axial_rho, _ = reference_chart(rho)
    if axial <= 0:
        raise ValueError('strictly trapped a>0 required')
    r_ref = float(np.sqrt(1.0 + rho * rho))
    r_ref_rho = rho / r_ref
    normal = _normal_coordinate(rho)
    normal_rho = -1.0 / axial
    delta = np.zeros(z.shape, dtype=float)
    delta_z = np.zeros_like(delta)
    delta_s = np.zeros_like(delta)
    delta_zz = np.zeros_like(delta)
    delta_sz = np.zeros_like(delta)
    for amplitude, direction in zip(metric.amplitudes, metric.directions):
        window, window_s = _plateau_jet(normal, direction.inner_radius, direction.outer_radius)
        if window == 0.0 and window_s == 0.0:
            continue
        w0 = np.array([direction.w(float(point), 0) for point in z], dtype=float)
        w1 = np.array([direction.w(float(point), 1) for point in z], dtype=float)
        w2 = np.array([direction.w(float(point), 2) for point in z], dtype=float)
        u0 = np.array([direction.U(float(point), 0) for point in z], dtype=float)
        u1 = np.array([direction.U(float(point), 1) for point in z], dtype=float)
        u2 = np.array([direction.U(float(point), 2) for point in z], dtype=float)
        poly = normal * w0 + (normal**3) * u0 / 6.0
        poly_s = w0 + (normal**2) * u0 / 2.0
        poly_z = normal * w1 + (normal**3) * u1 / 6.0
        poly_sz = w1 + (normal**2) * u1 / 2.0
        poly_zz = normal * w2 + (normal**3) * u2 / 6.0
        scale = float(amplitude)
        delta += scale * window * poly
        delta_z += scale * window * poly_z
        delta_s += scale * (window_s * poly + window * poly_s)
        delta_zz += scale * window * poly_zz
        delta_sz += scale * (window_s * poly_z + window * poly_sz)
    radius = r_ref + delta
    if not np.isfinite(radius).all() or np.any(radius <= 0):
        raise ValueError('actual radius must stay finite and positive')
    return {
        'a': axial,
        'a_rho': axial_rho,
        'r': radius,
        'r_z': delta_z,
        'r_rho': r_ref_rho + normal_rho * delta_s,
        'r_zz': delta_zz,
        'r_rhoz': normal_rho * delta_sz,
        's': normal,
        's_rho': normal_rho,
        'r_ref': r_ref,
    }


def imag_major_A2_transport_rhs(mass, angular, source_sign, jets):
    """Im(T_s major A2) under pure-imaginary major A1, original source signs."""
    if source_sign not in SOURCE_SIGNS:
        raise ValueError('source sign must be +1 or -1')
    mass = float(mass)
    angular = float(angular)
    if not np.isfinite([mass, angular]).all():
        raise ValueError('finite mass and angular eigenvalue required')
    radius = np.asarray(jets['r'], dtype=float)
    rhs = -mass * angular / 4.0 * (
        source_sign * jets['a']**2 * jets['r_rho'] + jets['r_z']) / (radius * radius)
    if not np.isfinite(rhs).all():
        raise ArithmeticError('nonfinite Im(T_s major A2) integrand')
    return rhs


def imag_major_A2_transport_rhs_z(mass, angular, source_sign, jets):
    """Partial z of Im(T_s major A2) at fixed rho."""
    if source_sign not in SOURCE_SIGNS:
        raise ValueError('source sign must be +1 or -1')
    mass = float(mass)
    angular = float(angular)
    radius = np.asarray(jets['r'], dtype=float)
    xi = source_sign * jets['a']**2 * jets['r_rho'] + jets['r_z']
    xi_z = source_sign * jets['a']**2 * jets['r_rhoz'] + jets['r_zz']
    factor = -mass * angular / 4.0
    rhs_z = factor * (xi_z * radius * radius - xi * 2.0 * radius * jets['r_z']) / radius**4
    if not np.isfinite(rhs_z).all():
        raise ArithmeticError('nonfinite z-derivative of Im(T_s major A2)')
    return rhs_z


def _incoming_imag_major_A2(
        family, z, *, mass, angular, rho_up, gauss_nodes, derivative):
    """Characteristic integral of Im(T_s major A2) or of its partial_z.

    Upstream Im(major A2) is the z-independent normalized reference value.
    The returned field is relative to that constant, which cancels in z
    derivatives and in g-minus-reference.
    """
    z = np.asarray(z, dtype=float)
    if z.ndim != 1 or z.size == 0 or not np.isfinite(z).all():
        raise ValueError('nonempty finite one-dimensional z required')
    mass = float(mass)
    angular = float(angular)
    rho_up = float(rho_up)
    if isinstance(gauss_nodes, (bool, np.bool_)) or int(gauss_nodes) != gauss_nodes or gauss_nodes < 1:
        raise ValueError('positive integer Gauss-Legendre count required')
    gauss_nodes = int(gauss_nodes)
    if rho_up < RHO_UP_MIN:
        raise ValueError('caller rho_up>=1.03 required')
    reference_chart(rho_up)
    metric = _as_one_direction(family)
    identity = profile_identity(family, include_normal_window=True)
    if _normal_coordinate(rho_up) > -max(d.outer_radius for d in metric.directions):
        raise ValueError('history must be flat at the fixed upstream slice')
    kernel = imag_major_A2_transport_rhs_z if derivative else imag_major_A2_transport_rhs
    name = 'h_z' if derivative else 'h'
    n_sign = len(SOURCE_SIGNS)
    values = np.zeros((n_sign, z.size), dtype=float)
    if mass == 0.0 or angular == 0.0:
        quadrature = {
            'rule': 'skipped: Im(T_s major A2) vanishes for m=0 or ell=0',
            'gauss_nodes': 0, 'n_nodes': 0, 'skipped': True,
            'certified_quadrature_error_bound': None,
        }
        return {
            'z': np.array(z, dtype=float, copy=True),
            name: values,
            'source_signs': SOURCE_SIGNS,
            'mass': mass,
            'angular': angular,
            'rho_up': rho_up,
            'profile_identity': identity,
            'quadrature': MappingProxyType(quadrature),
            'full_versus_difference': (
                'reference integrand is z-independent, so the incoming '
                'z-derivative equals the g-minus-reference value'),
            'is_full_e_minus2_coefficient': False,
            'numerical_C_M': None,
            'physical_local_gate': PHYSICAL_LOCAL_GATE,
        }
    panels = _panels(metric, rho_up)
    nodes, weights = _gauss_nodes(panels, gauss_nodes)
    for rho, weight in zip(nodes, weights):
        rho = float(rho)
        distance = _characteristic_distance(rho)
        for sign_index, source_sign in enumerate(SOURCE_SIGNS):
            z_char = z - source_sign * distance
            jets = radius_partial_jets(metric, rho, z_char)
            values[sign_index] += -weight * kernel(mass, angular, source_sign, jets)
    if not np.isfinite(values).all():
        raise ArithmeticError('nonfinite Im(major A2) quadrature')
    quadrature = {
        'rule': 'composite Gauss-Legendre on normal-window split panels',
        'gauss_nodes': gauss_nodes,
        'panels': [list(map(float, panel)) for panel in panels],
        'n_nodes': int(nodes.size),
        'skipped': False,
        'certified_quadrature_error_bound': None,
    }
    return {
        'z': np.array(z, dtype=float, copy=True),
        name: values,
        'source_signs': SOURCE_SIGNS,
        'mass': mass,
        'angular': angular,
        'rho_up': rho_up,
        'profile_identity': identity,
        'quadrature': MappingProxyType(quadrature),
        'full_versus_difference': (
            'reference integrand is z-independent, so the incoming '
            'z-derivative equals the g-minus-reference value'),
        'is_full_e_minus2_coefficient': False,
        'numerical_C_M': None,
        'physical_local_gate': PHYSICAL_LOCAL_GATE,
    }


def integrate_imag_major_A2(
        family, z, *, mass, angular, rho_up=ARCHIVE_RHO_UP,
        gauss_nodes=DEFAULT_GAUSS_NODES):
    """Incoming h_s relative to vanishing upstream Im(major A2)."""
    return _incoming_imag_major_A2(
        family, z, mass=mass, angular=angular, rho_up=rho_up,
        gauss_nodes=gauss_nodes, derivative=False)


def integrate_imag_major_A2_z(
        family, z, *, mass, angular, rho_up=ARCHIVE_RHO_UP,
        gauss_nodes=DEFAULT_GAUSS_NODES):
    """Incoming h_{s,z} from the constructed Im(T_s major A2) integrand.

    This is a pure-geometry quadrature. It is not C_M, not a source sum, and
    not the full E^{-2} N or beta coefficient.
    """
    return _incoming_imag_major_A2(
        family, z, mass=mass, angular=angular, rho_up=rho_up,
        gauss_nodes=gauss_nodes, derivative=True)


def _axial_cutoff_region(profile, z):
    distance = abs(float(z) - float(profile.center))
    inner = float(profile.inner)
    outer = float(profile.outer)
    if distance <= inner:
        return 'interior'
    if distance >= outer:
        return 'exterior'
    return 'ramp'


def chebyshev_profile_point_enclosures(family, z, derivative_order=1):
    """Cutoff-aware interval jets versus a binary64 evaluation.

    The numerical LocalAxialFunction value is not an analytic target and need
    not lie in an 80-bit ball of the exact polynomial. Directed reconstruction
    error is recorded. Containment of the float is diagnostic only.
    """
    from math import factorial
    from flint import arb, ctx
    from .nsc_ks_chebyshev_jet_bound import (
        chebyshev_endpoint_derivative_bounds,
        local_axial_interval_series,
        polynomial_derivative_on_interval,
    )
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError('owned LocalIncomingFamily required')
    z = np.asarray(z, dtype=float)
    if z.ndim != 1 or z.size == 0 or not np.isfinite(z).all():
        raise ValueError('nonempty finite one-dimensional z required')
    if isinstance(derivative_order, bool) or derivative_order not in (0, 1, 2):
        raise ValueError('derivative order 0, 1 or 2 required')
    rows = []
    with ctx.workprec(80):
        for name, profile in zip(('w', 'U'), family.functions):
            bounds = chebyshev_endpoint_derivative_bounds(
                profile, max(derivative_order + 1, 9))
            points = []
            for point in z:
                region = _axial_cutoff_region(profile, point)
                series = local_axial_interval_series(
                    profile, point, derivative_order, bounds)
                enclosure = factorial(derivative_order) * series[derivative_order]
                value = float(profile(float(point), derivative_order))
                mid = float(enclosure.mid())
                rad = float(enclosure.rad())
                abs_err = abs(value - mid)
                poly_only = polynomial_derivative_on_interval(
                    profile, point, derivative_order, bounds)
                points.append({
                    'z': float(point),
                    'region': region,
                    'numerical_value': value,
                    'cutoff_enclosure_mid': mid,
                    'cutoff_enclosure_rad': rad,
                    'reconstruction_abs_error': abs_err,
                    'directed_distance_to_cutoff_enclosure': max(0.0, abs_err - rad),
                    'float64_contained_in_cutoff_enclosure': bool(enclosure.contains(value)),
                    'analytic_zero_contained': bool(enclosure.contains(0)),
                    'polynomial_only_mid': float(poly_only.mid()),
                    'polynomial_only_ignores_cutoff': region != 'interior',
                })
            rows.append({
                'profile': name,
                'derivative_order': derivative_order,
                'points': points,
            })
    return tuple(rows)


def finite_history_uv_coefficient_report():
    """Executable summary of constructed jets, bookkeeping, and explicit gaps."""
    identities = coefficient_identities()
    current = e_minus2_current_bookkeeping()
    sigma = e_minus2_sigma_contraction()
    return MappingProxyType({
        'schema': 'NSC-KS-FINITE-HISTORY-UV-COEFFICIENTS-v3',
        'identities': identities,
        'e_minus2_current': current,
        'e_minus2_sigma': sigma,
        'upstream': vacuum_upstream_conditions(),
        'gaps': named_gaps(),
        'vertices': dict(sigma['vertices']),
        'physical_local_gate': PHYSICAL_LOCAL_GATE,
        'certificate_from_asymptotic_orders_alone': False,
        'new_action_term': False,
        'Gamma_rest_assigned': False,
        'existing_source_ad759_evaluation_changed': False,
        'leading_paired_e_minus2_cancels': True,
        'C_M_or_physical_gate_from_leading_cancellation': False,
    })
