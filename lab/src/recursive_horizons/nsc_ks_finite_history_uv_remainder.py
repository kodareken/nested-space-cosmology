"""Changed-history UV remainder: first noncancelling paired order and M=4 majorant.

Successor of the v4 coefficient owner. The leading paired E^{-2} vacuum
coefficient is already cancelled. This module continues the same L0
recurrence far enough to identify the first noncancelling paired N,beta
    order, records the formal truncated-envelope L0 telescope, and supplies a
    conditional M=4 a-posteriori remainder majorant on a declared current-history
    slab. Production A3 and A4 are not constructed here. Major A2 and A3 stay
    transported coefficients; they are not replaced
by local metric jets. The Fermi thermal tail is bounded separately from its
exponential occupation. Numerical C4 and C_M remain None until the named
transport integrals of L0 A4 are supplied.

The v2 successor represents A2, A3, n4, C4, C_M, current-history L0 A4 H2
integrals, the upstream higher-order H2 remainder, transported majors and
the vacuum N/beta tail, and binds those names to the evolved state law,
current profile, horizon preparation, source inventory and KS subtraction.
It computes only quantities available in a small authenticated control.
Unavailable values stay None and OPEN. The E^{-2} cancellation and E^{-3}
first allowed paired order are preserved and are not tail bounds.
"""
from collections.abc import Mapping
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from types import MappingProxyType
import json
import math

import mpmath as mp
import numpy as np
import sympy as sp

from .nsc_dirac_source_phase import SOURCE_SIGNS
from .nsc_ks_current_history_bounds import radius_bounds
from .nsc_ks_finite_history_uv_coefficients import (
    ARCHIVE_RHO_UP,
    CURRENT_HISTORY_IDENTITY_PREFIX,
    CUTOFF_BRIDGE_CONTROL,
    INVENTORY_RECORD,
    V2_RECORD,
    _apply_L0,
    _major_index,
    _minor_index,
    _minor_A1,
    _projector,
    _real_residual,
    _sym_pauli,
    _sym_potential,
    _T,
    coefficient_identities,
    confirm_raw_ks_vertices,
    e_minus2_sigma_contraction,
    original_angular_pair_ledger,
    paired_leading_e_minus2_weight,
)
from .nsc_ks_profile_identity import profile_identity
from .nsc_ks_residual_error import (
    matter_error_bounds,
    pointwise_field_bounds,
    propagate_h2,
    pure_radius_commutator_integrals,
)
from .nsc_local_incoming_family import LocalIncomingFamily
from .nsc_pg_ks_metric_pullback import reference_chart
from .nsc_transmitting_dirac_domain import S3


V3_RECORD = 'results/development/nsc-ks-current-uv-coefficients-v3.json'
V4_RECORD = 'results/development/nsc-ks-current-uv-coefficients-v4.json'
FROZEN_V2_SHA256 = 'de8e67e4cd89eb0331aff4e18d577220b83ddca94c0c6bac6d0c5a0fe32585df'
HISTORY_RECORD = 'results/development/nsc-ks-gate-history-lm-broyden.json'
CAUCHY_RECORD = 'results/development/nsc-mode-resolved-cauchy-state.json'
RESTART_RECORD = 'results/development/nsc-compact-matched-restart.json'
PROVISIONAL_UV_ALLOCATION = 5e-12
REMAINDER_ORDER = 4
FIRST_NONCANCELLING_PAIRED_ORDER = 3
CURRENT_HISTORY_IDENTITY = (
    '0b0e4cedfdb695d342c7f7347740631e916bdaee3a977df56808c725b25dacb0')
PHYSICAL_LOCAL_GATE = (
    'OPEN: first noncancelling paired vacuum order is e^{-3}; numerical C4, '
    'C_M and the local incoming gate remain unbounded pending L0 A4 / major '
    'A2-A3 transport norms')
MISSING_L0_A4 = (
    'H2 integrals of L0 A4 on the declared current-history slab, which require '
    'the characteristic integrals of major A2 (real and imaginary) and the '
    'subsequent major-A3 transport; Re(T_s major A2) also uses the transported '
    'major-A1 phase q. These are not local metric jets')
MISSING_UPSTREAM_TAIL = (
    'upstream affine-expansion H2 remainder for inverse-energy orders greater '
    'than 4 at rho_up, for the same mathematical vacuum column')
STATE_LAW = (
    'C_Sigma[g]=F[g] C_src F[g]^dagger with explicit dC_src=0 and retarded '
    'delta C_Sigma from the same evolution')
STATE_LAW_OWNER = 'src/recursive_horizons/nsc_evolved_incoming_state.py'
PREPARATION_OWNER = 'src/recursive_horizons/nsc_paired_horizon_preparation.py'
SUBTRACTION_OWNER = 'src/recursive_horizons/nsc_common_subtracted_ks_source.py'
PILOT_RECORD = 'results/development/nsc-ks-current-uv-remainder-pilot.json'
SUCCESSOR_SCHEMA = 'NSC-KS-CURRENT-UV-REMAINDER-v2'
PILOT_SCHEMA = 'NSC-KS-CURRENT-UV-REMAINDER-PILOT-v1'
SUCCESSOR_STATUS = (
    'OPEN: changed-history UV quantities represented and bound to state law, '
    'profile, preparation, inventory and subtraction; E^{-2} cancellation and '
    'E^{-3} first allowed paired order preserved; numerical C4, C_M, L0 A4 H2 '
    'integrals, upstream H2 remainder and vacuum N/beta tail remain missing')
QUANTITY_NAMES = (
    'A2', 'A3', 'n4', 'C4', 'C_M',
    'current_history_L0_A4_H2_integrals',
    'upstream_higher_order_H2_remainder',
    'transported_majors',
    'vacuum_N_beta_tail',
)
NONE_REQUIRED_QUANTITIES = QUANTITY_NAMES


def _root():
    return Path(__file__).resolve().parents[2]


def _require_positive(value, name):
    if isinstance(value, (bool, np.bool_)) or value is None:
        raise ValueError('explicit positive ' + name + ' required')
    number = float(value)
    if not np.isfinite(number) or number <= 0:
        raise ValueError('positive finite ' + name + ' required')
    return number


def _h2_triple(values, name):
    if values is None:
        return None
    if not isinstance(values, (tuple, list)) or len(values) != 3:
        raise ValueError(name + ' must be three H2 component bounds')
    out = []
    for item in values:
        if isinstance(item, (bool, np.bool_)) or item is None:
            raise ValueError('explicit nonnegative ' + name + ' required')
        number = float(item)
        if not np.isfinite(number) or number < 0:
            raise ValueError('nonnegative finite ' + name + ' required')
        out.append(number)
    return tuple(out)


@lru_cache(maxsize=1)
def remainder_recurrence_identities():
    """Owned M=2 identities plus a formal L0 telescope through M=4.

    The M=2 telescope on constructed A0,A1,A2 is the v4 identity. Higher
    orders use the same L0 on a generic transported coefficient (Phi, Nu)
    and the finite recurrence telescope whose last term is Pi_{-s} L0 A_M.
    Minor A3 retains major A2. The M=3/M=4 symbols test the algebraic
    telescope only; they do not construct the current-history A3 or A4.
    """
    parent = coefficient_identities()
    rho, z = sp.symbols('rho z', real=True)
    mass, angular = sp.symbols('m ell', real=True)
    axial = sp.Function('a', real=True)(rho)
    radius = sp.Function('r', real=True)(rho, z)
    energy = sp.symbols('e', positive=True)
    s1, s2, s3 = _sym_pauli()
    potential = _sym_potential(mass, angular, radius)
    identities = {
        's+1/L0_telescope_M2': list(parent['residuals']['s+1/L0_telescope_M2']),
        's-1/L0_telescope_M2': list(parent['residuals']['s-1/L0_telescope_M2']),
        's+1/minor_A3_depends_on_major_A2': parent['residuals']['s+1/minor_A3_depends_on_major_A2'],
        's-1/minor_A3_depends_on_major_A2': parent['residuals']['s-1/minor_A3_depends_on_major_A2'],
    }
    phase = sp.Function('q', real=True)(rho, z)
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        major = _major_index(source_sign)
        minor = _minor_index(source_sign)
        occupied, vacant = _projector(source_sign, s3)
        phi = sp.Function('Phi')(rho, z)
        nu = sp.Function('Nu')(rho, z)
        coefficient = sp.zeros(2, 1)
        coefficient[major] = phi
        coefficient[minor] = nu
        transport_rhs = (sp.I / axial) * potential[major, minor] * nu
        a_rho = source_sign / axial**2 * phi.diff(z) + transport_rhs
        l0 = _apply_L0(coefficient, axial, potential, s3, rho, z)
        l0_on_shell = l0.subs(phi.diff(rho), a_rho)
        transport_defect = (occupied * l0)[major] - (
            _T(phi, source_sign, axial, rho, z) - transport_rhs)
        identities['%s/Pi_s_L0_Aj_is_transport_defect' % tag] = _real_residual(transport_defect)
        next_minor = axial**2 / (2 * sp.I) * (vacant * l0_on_shell)[minor]
        identities['%s/minor_Ajp1_depends_on_major_Aj' % tag] = (
            'depends' if phi in next_minor.atoms(sp.Function) else 'independent')
        nxt = sp.zeros(2, 1)
        nxt[major] = sp.Function('Psi')(rho, z)
        nxt[minor] = next_minor
        cancel = l0_on_shell - (2 * sp.I / axial**2) * (vacant * nxt)
        identities['%s/L0_Aj_equals_2i_a2_Pi_minus_Ajp1' % tag] = [
            _real_residual(component) for component in cancel]
        # Re(T_s major A2) retains the transported phase q; not a local jet.
        a1 = sp.zeros(2, 1)
        a1[major] = sp.I * phase
        a1[minor] = _minor_A1(axial, mass, angular, radius, source_sign)
        square = mass**2 + angular**2 / radius**2
        a1_rho = source_sign / axial**2 * a1[major].diff(z) - sp.I / 2 * square
        l0_a1 = _apply_L0(a1, axial, potential, s3, rho, z).subs(
            a1[major].diff(rho), a1_rho)
        minor_a2 = axial**2 / (2 * sp.I) * (vacant * l0_a1)[minor]
        rhs_a2 = (sp.I / axial) * potential[major, minor] * minor_a2
        imag_rhs = sp.im(sp.expand_complex(rhs_a2))
        real_rhs = sp.re(sp.expand_complex(rhs_a2))
        parent_imag = -mass * angular / 4 * (
            source_sign * axial**2 * radius.diff(rho) + radius.diff(z)) / radius**2
        identities['%s/Im_Ts_major_A2_parent' % tag] = _real_residual(imag_rhs - parent_imag)
        identities['%s/Re_Ts_major_A2_depends_on_q' % tag] = (
            'depends' if phase in real_rhs.atoms(sp.Function) else 'independent')
        identities['%s/Im_Ts_major_A2_independent_of_q' % tag] = (
            'independent' if phase not in parent_imag.atoms(sp.Function) else 'depends')
    axial_const = sp.symbols('a', positive=True)
    energy_c = sp.symbols('e', positive=True)
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        occupied = (sp.eye(2) + source_sign * sp.diag(1, -1)) / 2
        vacant = sp.eye(2) - occupied
        for order in (3, 4):
            coeff = [sp.Matrix(sp.symbols('A%d_%d_0:2' % (order, index)))
                     for index in range(order + 1)]
            coeff[0] = occupied * coeff[0]
            last = vacant * sp.Matrix(sp.symbols('L0AM_%d_0:2' % order))
            applied = [
                2 * sp.I / axial_const**2 * vacant * coeff[index + 1]
                for index in range(order)] + [last]
            defect = sum(
                (applied[index] / energy_c**index for index in range(order + 1)),
                sp.zeros(2, 1))
            defect -= (2 * sp.I * energy_c / axial_const**2) * vacant * sum(
                (coeff[index] / energy_c**index for index in range(order + 1)),
                sp.zeros(2, 1))
            identities['%s/L0_telescope_M%d' % (tag, order)] = [
                str(sp.simplify(component)) for component in defect - last / energy_c**order]
            identities['%s/defect_is_Pi_minus_L0_AM_over_eM' % tag] = True
    residual_keys = [
        key for key, value in identities.items()
        if not key.endswith(('depends_on_major_A2', 'depends_on_major_Aj',
                             'depends_on_q', 'independent_of_q',
                             'defect_is_Pi_minus_L0_AM_over_eM'))]
    leaves = []
    for key in residual_keys:
        value = identities[key]
        if isinstance(value, list):
            leaves.extend(value)
        else:
            leaves.append(value)
    if set(leaves) != {'0'}:
        raise ArithmeticError(identities)
    for key, value in identities.items():
        if key.endswith('depends_on_major_A2') or key.endswith('depends_on_major_Aj'):
            if value != 'depends':
                raise ArithmeticError('minor A_{j+1} must retain transported major A_j')
        if key.endswith('depends_on_q') and value != 'depends':
            raise ArithmeticError('Re(T_s major A2) must retain the transported phase')
        if key.endswith('independent_of_q') and value != 'independent':
            raise ArithmeticError('Im(T_s major A2) parent must stay independent of q')
    return MappingProxyType({
        'residuals': identities,
        'expansion_order': REMAINDER_ORDER,
        'defect': (
            'formal e^{-4} Pi_{-s} L0 A4 identity conditional on constructed '
            'A0..A4; production A3/A4 are not constructed here'),
        'A3_locally_algebraic_in_metric_jets': False,
        'major_A2_replaced_by_local_metric_jet': False,
        'dummy_symbol_L0_A_j': True,
        'formal_symbolic_telescope': True,
        'production_A3_A4_constructed': False,
        'physical_source_columns_replaced': False,
        'telescope': (
            '(L0 - 2 i e / a^2 Pi_{-s}) sum_{j=0}^M e^{-j} A_j = e^{-M} L0 A_M '
            'with M=4, conditional on coefficients satisfying the recurrence; '
            'the production A3/A4 construction remains OPEN'),
    })


@lru_cache(maxsize=1)
def first_noncancelling_paired_coefficient():
    """Identify the first paired N,beta inverse-energy order that does not cancel.

    Angular involution gives parity (-1)^{p+1} at e^{-p}. The e^{-2} coefficient
    is odd and cancels for equal-mu original pairs. e^{-3} is even, so it is the
    first allowed paired order. Its formula retains transported A2, A3 and n4;
    n4's continuity identity still contains minor A3. This is not a proof that
    the complete cubic coefficient is nonzero, and not a numerical C4.
    """
    rho, z = sp.symbols('rho z', real=True)
    axial = sp.Function('a', real=True)(rho)
    energy = sp.Symbol('e', positive=True)
    n2 = sp.Function('n2', real=True)(rho, z)
    n3 = sp.Function('n3', real=True)(rho, z)
    n4 = sp.Function('n4', real=True)(rho, z)
    minor1 = sp.Function('L')(rho, z)
    minor2 = sp.Function('B')(rho, z)
    minor3 = sp.Function('C')(rho, z)
    mass, angular, radius = sp.symbols('m ell r', real=True)
    identities = {}
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        minor_sq = (
            sp.conjugate(minor1) * minor1 / energy**2
            + 2 * sp.re(sp.conjugate(minor1) * minor2) / energy**3
            + (sp.conjugate(minor2) * minor2
               + 2 * sp.re(sp.conjugate(minor1) * minor3)) / energy**4)
        density = 1 + n2 / energy**2 + n3 / energy**3 + n4 / energy**4
        current = source_sign * (density - 2 * minor_sq)
        continuity = density.diff(rho) - current.diff(z) / axial**2
        order4 = sp.expand(continuity * energy**4).subs(energy, 0)
        ts_n4 = _T(n4, source_sign, axial, rho, z)
        re_lc = sp.re(sp.conjugate(minor1) * minor3)
        expected = ts_n4 + source_sign / axial**2 * sp.diff(
            2 * sp.conjugate(minor2) * minor2 + 4 * re_lc, z)
        identities['%s/continuity_e_minus4_retains_minor_A3' % tag] = _real_residual(
            order4 - expected)
        identities['%s/n4_depends_on_minor_A3' % tag] = (
            'depends' if minor3 in order4.atoms(sp.Function) else 'independent')
    leaves = [
        value for key, value in identities.items()
        if not key.endswith('depends_on_minor_A3')]
    if set(leaves) != {'0'}:
        raise ArithmeticError(identities)
    if any(value != 'depends' for key, value in identities.items()
           if key.endswith('depends_on_minor_A3')):
        raise ArithmeticError('n4 must retain minor A3')
    involution = {}
    s3 = sp.diag(1, -1)
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        conjugation = s3 * source_sign
        coeff = [sp.eye(2)[:, _major_index(source_sign)]] + [
            sp.Matrix(sp.symbols('a%d_%d_0:2' % (source_sign, index)))
            for index in range(1, 5)]
        deriv = [sp.zeros(2, 1)] + [
            sp.Matrix(sp.symbols('z%d_%d_0:2' % (source_sign, index)))
            for index in range(1, 5)]
        partner = [((-1)**index) * conjugation * column.conjugate()
                   for index, column in enumerate(coeff)]
        partner_z = [((-1)**index) * conjugation * column.conjugate()
                     for index, column in enumerate(deriv)]

        def current_coeff(columns, derivatives, weight, order):
            interference = sum(
                (columns[index].H * weight * derivatives[order - index]
                 - derivatives[index].H * weight * columns[order - index]
                 for index in range(order + 1)),
                sp.zeros(1, 1))[0] / (2 * sp.I)
            density = sum(
                (columns[index].H * weight * columns[order + 1 - index]
                 for index in range(order + 2)),
                sp.zeros(1, 1))[0]
            return interference - source_sign * density

        def density_coeff(columns, vertex, order):
            density = sum((
                columns[index] * columns[order - index].H
                for index in range(order + 1)), sp.zeros(2))
            return sp.trace(vertex * density)

        for order in (2, 3):
            sign = (-1)**(order + 1)
            beta_old = current_coeff(coeff, deriv, sp.eye(2), order)
            beta_new = current_coeff(partner, partner_z, sp.eye(2), order)
            involution['%s/beta_parity_p%d' % (tag, order)] = str(
                sp.expand(beta_new - sign * beta_old))
            v_plus = _sym_potential(mass, angular, radius)
            v_minus = _sym_potential(mass, -angular, radius)
            n_old = (density_coeff(coeff, v_plus, order)
                     + current_coeff(coeff, deriv, s3, order) / axial)
            n_new = (density_coeff(partner, v_minus, order)
                     + current_coeff(partner, partner_z, s3, order) / axial)
            involution['%s/N_parity_p%d' % (tag, order)] = str(
                sp.expand(n_new - sign * n_old))
    if set(involution.values()) != {'0'}:
        raise ArithmeticError(involution)
    sigma = e_minus2_sigma_contraction()
    return MappingProxyType({
        'first_noncancelling_paired_inverse_energy_order': FIRST_NONCANCELLING_PAIRED_ORDER,
        'full_N_beta_coefficient_parity': 'at e^{-p}: (-1)^{p+1}',
        'leading_paired_e_minus2_cancels_for_equal_mu': True,
        'e_minus3_even_in_ell': True,
        'e_minus4_parity_evaluated': False,
        'complete_e_minus3_coefficient_proven_nonzero': False,
        'A3_replaced_by_local_metric_jet': False,
        'n4_eliminated_without_minor_A3': False,
        'formulas': {
            'J_3': (
                'Im(A0^dag A3_z + A1^dag A2_z + A2^dag A1_z) - s n4, '
                'with A0_z=0'),
            'n4_transport': (
                'T_s n4 = -(s/a^2) d_z (2 |minor A2|^2 + 4 Re(conj(minor A1) minor A3)); '
                'upstream n4=0 does not remove minor A3'),
            'action_N3': (
                'mu/(2 pi) times the occupied contraction of V at e^{-3} and '
                '(S3/a) J_3, action minus once; still a functional of transported A2,A3'),
            'action_beta3': (
                '-mu/(2 pi) times the occupied J_3; action minus once'),
            'C4': (
                'equal-mu sum of the two original angular signs of those '
                'e^{-3} action densities on I'),
        },
        'continuity_residuals': identities,
        'involution_residuals': involution,
        'angular_pairs': sigma['angular_pairs'],
        'numerical_C4': None,
        'C4_or_C_M_bounded': False,
    })


def truncated_envelope_L0_defect(order=REMAINDER_ORDER):
    """Record the conditional formal M=4 defect and its missing construction."""
    if order != REMAINDER_ORDER:
        raise ValueError('this remainder owner implements M=4')
    identities = remainder_recurrence_identities()
    return MappingProxyType({
        'M': order,
        'defect': identities['defect'],
        'field_H2_error_order': 'O(e^{-4}) after a finite C_M',
        'differentiated_current_error_order': (
            'O(e^{-3}) after restoring the source carrier; M=3 alone leaves '
            'an O(e^{-2}) current remainder and does not prove the paired tail'),
        'dummy_symbol_L0_A_j': identities['dummy_symbol_L0_A_j'],
        'formal_symbolic_telescope': identities['formal_symbolic_telescope'],
        'production_A3_A4_constructed': identities[
            'production_A3_A4_constructed'],
        'residuals': identities['residuals'],
    })


def _scalar_c_m(h2):
    return max(h2)


def conditional_m4_remainder_majorant(
        initial_h2=None, l0_A4_residual_integrals=None,
        Bz_integral=None, Bzz_integral=None, *,
        period=None, energy_lower=None):
    """A-posteriori M=4 H2 majorant, conditional on the named continuous inputs.

    Inputs are e-independent: initial_h2 bounds e^4 (X-X^{(4)}) at rho_up,
    l0_A4_residual_integrals bound the H2 integrals of L0 A4, and Bz/Bzz are
    the owned commutator integrals. Then
    ||X-X^{(4)}||_{H^2} <= C_M e^{-4} with C_M the propagated H2 triple.
    Missing inputs stay explicit. No background-tail constant is transferred.
    """
    missing = []
    if initial_h2 is None:
        missing.append(MISSING_UPSTREAM_TAIL)
    else:
        initial_h2 = _h2_triple(initial_h2, 'initial H2')
    if l0_A4_residual_integrals is None:
        missing.append(MISSING_L0_A4)
    else:
        l0_A4_residual_integrals = _h2_triple(
            l0_A4_residual_integrals, 'L0 A4 residual integrals')
    if Bz_integral is None:
        missing.append('owned B_z commutator integral on the current-history slab')
    if Bzz_integral is None:
        missing.append('owned B_zz commutator integral on the current-history slab')
    status = {
        'M': REMAINDER_ORDER,
        'numerical_C_M': None,
        'C_M_h2': None,
        'missing_inputs': tuple(missing),
        'A3_replaced_by_local_metric_jet': False,
        'background_tail_certificate_transferred': False,
        'fitted_decay': False,
        'physical_local_gate': PHYSICAL_LOCAL_GATE,
    }
    if missing:
        return MappingProxyType(status)
    bound = propagate_h2(
        initial_h2, l0_A4_residual_integrals, Bz_integral, Bzz_integral)
    bound_float = tuple(float(component) for component in bound)
    field = None
    if period is not None and energy_lower is not None:
        raw_field = pointwise_field_bounds(bound, period, energy_lower)
        field = {
            key: (value if isinstance(value, str) else float(value))
            for key, value in raw_field.items()
        }
    status.update({
        'numerical_C_M': _scalar_c_m(bound_float),
        'C_M_h2': bound_float,
        'pointwise_field': field,
        'missing_inputs': (),
    })
    return MappingProxyType(status)


def current_history_commutator_integrals(family, absolute_angular, rho_up=ARCHIVE_RHO_UP):
    """B_z, B_zz integrals from the owned radius family; not a C_M."""
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError('owned LocalIncomingFamily required')
    angular = _require_positive(absolute_angular, 'absolute angular')
    bounds = radius_bounds(family, rho_up)
    first, second = pure_radius_commutator_integrals(
        bounds['w']['profile_bounds'], bounds['U']['profile_bounds'],
        bounds['normal_support'], bounds['radius_lower'], angular)
    return MappingProxyType({
        'Bz_integral': float(first),
        'Bzz_integral': float(second),
        'radius_lower': float(bounds['radius_lower']),
        'axial_lower': float(bounds['axial_lower']),
        'normal_support': float(bounds['normal_support']),
        'absolute_angular': angular,
        'certified_L0_A4': False,
    })


def _state_scales(root=None):
    root = Path(root) if root is not None else _root()
    restart = json.loads((root / RESTART_RECORD).read_text())
    config = restart['scattering_provenance']['config']
    kappa = float(config['surface_gravity'])
    omega = float(config['omega'])
    if min(kappa, omega) <= 0:
        raise ValueError('owned positive surface gravity and Omega required')
    return MappingProxyType({
        'surface_gravity': kappa,
        'omega': omega,
        'source_law': '2*exp(-pi*E/kappa)+exp(-2*pi*E/(Omega*kappa))',
        'record': RESTART_RECORD,
        'record_sha256': sha256((root / RESTART_RECORD).read_bytes()).hexdigest(),
    })


def fermi_occupation_bound(energy, kappa, omega):
    """Owned high-energy trace-norm majorant of C_src - C_vac."""
    energy = _require_positive(energy, 'energy')
    kappa = _require_positive(kappa, 'surface gravity')
    omega = _require_positive(omega, 'Omega')
    return 2.0 * math.exp(-math.pi * energy / kappa) + math.exp(
        -2.0 * math.pi * energy / (omega * kappa))


def _exponential_moments(alpha, cutoff, dps=80):
    """Directed integrals of e^{-alpha E} and E e^{-alpha E} from cutoff to infinity."""
    if mp.iv.mpf(alpha) <= 0 or mp.iv.mpf(cutoff) <= 0:
        raise ValueError('positive decay rate and cutoff required')
    exponential = mp.iv.exp(-alpha * cutoff)
    moment0 = exponential / alpha
    moment1 = exponential * (cutoff / alpha + 1 / alpha**2)
    return moment0, moment1


def fermi_thermal_tail_bound(
        mass, angular, multiplicity, cutoffs=(160.0, 320.0), *,
        kappa=None, omega=None, dps=80):
    """N,beta occupation tail beyond the archived 160/320 splits.

    Vacuum envelope coefficients are not used. The bound is the owned
    exponential occupation times the N,beta vertices, with |F|<=1 and
    |F_z| <= E + m + |ell|/r_Sigma. It is not a background-tail certificate.
    """
    mass = float(mass)
    angular = float(angular)
    multiplicity = _require_positive(multiplicity, 'multiplicity')
    if mass < 0:
        raise ValueError('nonnegative mass required')
    scales = _state_scales() if kappa is None or omega is None else None
    kappa = float(scales['surface_gravity'] if kappa is None else kappa)
    omega = float(scales['omega'] if omega is None else omega)
    kappa = _require_positive(kappa, 'surface gravity')
    omega = _require_positive(omega, 'Omega')
    _, _, axial, _, _ = reference_chart(1.0)
    radius = float(np.sqrt(2.0))
    if axial <= 0:
        raise ValueError('positive Sigma axial required')
    envelope = abs(mass) + abs(angular) / radius
    rows = []
    with mp.workdps(dps):
        mp.iv.dps = dps
        for cutoff in cutoffs:
            cutoff = _require_positive(cutoff, 'cutoff')
            moment0 = mp.iv.mpf(0)
            moment1 = mp.iv.mpf(0)
            for coefficient, alpha in (
                    (2, mp.iv.pi / kappa),
                    (1, 2 * mp.iv.pi / (omega * kappa))):
                m0, m1 = _exponential_moments(alpha, cutoff, dps=dps)
                moment0 += coefficient * m0
                moment1 += coefficient * m1
            # Two source signs, measure de/(2 pi), |F|<=1, |F_z|<=E+envelope.
            density = 2 * moment0 / (2 * mp.iv.pi)
            current = 2 * (moment1 + envelope * moment0) / (2 * mp.iv.pi)
            n_bound = multiplicity * ((abs(mass) + abs(angular) / radius) * density + current / axial)
            beta_bound = multiplicity * current
            log_hi = lambda value: float(mp.log(mp.make_mpf(value._mpi_[1])))
            display = lambda value: mp.nstr(mp.make_mpf(value._mpi_[1]), 20)
            log_n, log_beta = log_hi(n_bound), log_hi(beta_bound)
            rows.append({
                'cutoff': float(cutoff),
                'N': None,
                'beta': None,
                'N_decimal_upper': display(n_bound),
                'beta_decimal_upper': display(beta_bound),
                'log_N': log_n,
                'log_beta': log_beta,
                'N_display': display(n_bound),
                'beta_display': display(beta_bound),
                'float64_underflow': True,
                'underflow_is_not_exact_zero': True,
                'below_provisional_allocation': bool(
                    log_n < math.log(PROVISIONAL_UV_ALLOCATION)
                    and log_beta < math.log(PROVISIONAL_UV_ALLOCATION)),
            })
    return MappingProxyType({
        'source_law': '2*exp(-pi*E/kappa)+exp(-2*pi*E/(Omega*kappa))',
        'surface_gravity': kappa,
        'omega': omega,
        'mass': mass,
        'angular': angular,
        'multiplicity_per_signed_family': multiplicity,
        'source_signs': 2,
        'measure': 'de/(2 pi)',
        'vertices': 'N: V and S3/a; beta: -I; action minus once',
        'field_norm': 1.0,
        'axial_norm': 'E + m + |ell|/r_Sigma',
        'vacuum_envelope_used': False,
        'background_tail_certificate_transferred': False,
        'coherent_remainder_separate': True,
        'cutoffs': tuple(rows),
    })


def manufactured_defect_majorant_control():
    """Diagnostic smoke fixtures; neither constructs the current-history A4."""
    length = 2 / 5
    wave = 2 * math.pi / length
    amplitude = 3 / 10
    duration = 1 / 5
    nodes = 4001
    z = np.linspace(0.0, length, nodes, endpoint=False)
    profile = amplitude * np.sin(wave * z)
    l2 = math.sqrt(float(np.mean(profile * profile) * length))
    residual = (l2 * duration, wave * l2 * duration, wave * wave * l2 * duration)
    actual = (duration * l2, duration * wave * l2, duration * wave * wave * l2)
    bound = propagate_h2((0, 0, 0), residual, 0, 0)
    bound_float = tuple(float(component) for component in bound)
    enclosed = all(actual[index] <= bound_float[index] * (1 + 1e-12) for index in range(3))
    omitted = propagate_h2((0, 0, 0), (0, 0, 0), 0, 0)
    omitted_float = tuple(float(component) for component in omitted)
    omission_detected = actual[0] > omitted_float[0]
    # Small-matrix Duhamel on the singular projector generator
    # Y' = 2 i e Pi_{-s} Y + e^{-4} F, Y(0)=0. Pi_{-s}=diag(0,1) for s=+1.
    energy = 5.0
    forcing = np.array([0.4, -0.25], dtype=complex)
    span = 0.03
    scale = forcing / energy**4
    exact = np.array([
        span * scale[0],
        scale[1] * (np.exp(2j * energy * span) - 1.0) / (2j * energy),
    ], dtype=complex)
    residual_norm = float(np.linalg.norm(forcing)) * span
    matrix_bound = propagate_h2((0, 0, 0), (residual_norm, 0, 0), 0, 0)
    matrix_actual = float(np.linalg.norm(exact))
    scaled = float(matrix_bound[0]) / energy**4
    matrix_enclosed = matrix_actual <= scaled * (1 + 1e-8)
    return MappingProxyType({
        'periodic_manufactured_enclosed': bool(enclosed),
        'periodic_actual_h2': actual,
        'periodic_bound_h2': bound_float,
        'omitting_residual_fails_to_enclose': bool(omission_detected),
        'small_matrix_enclosed': bool(matrix_enclosed),
        'small_matrix_actual': matrix_actual,
        'small_matrix_bound_unscaled': float(matrix_bound[0]),
        'small_matrix_bound_on_remainder': scaled,
        'small_matrix_energy': energy,
        'dummy_symbol_L0_A_j': False,
        'current_history_C_M': None,
        'tautological_zero_commutator_residual_fixture': True,
        'production_L0_A4_control': False,
        'certificate_use': False,
    })


def negative_pair_and_omission_controls():
    """Broken equal-mu pairing, and omission of the first noncancelling order."""
    first = first_noncancelling_paired_coefficient()
    h = np.array([1.25, -0.5, 0.25])
    mu = 6.0
    equal = paired_leading_e_minus2_weight(h, -h, mu, mu)
    broken = paired_leading_e_minus2_weight(h, -h, mu, 0.5 * mu)
    # Even-parity model of the e^{-3} coefficient: c(+ell)=c(-ell)=g.
    cubic = np.array([0.7, -0.2, 0.05])
    paired_cubic = mu * cubic + mu * cubic
    omitted = np.zeros_like(cubic)
    omission_detected = bool(
        first['first_noncancelling_paired_inverse_energy_order'] == 3
        and np.max(np.abs(paired_cubic)) > 0
        and np.max(np.abs(omitted)) == 0)
    return MappingProxyType({
        'equal_mu_e_minus2_max_abs': float(np.max(np.abs(equal))),
        'broken_mu_e_minus2_max_abs': float(np.max(np.abs(broken))),
        'broken_mu_cancels': bool(np.allclose(broken, 0.0)),
        'paired_e_minus3_max_abs': float(np.max(np.abs(paired_cubic))),
        'omitted_e_minus3_max_abs': float(np.max(np.abs(omitted))),
        'omitting_e_minus3_leaves_false_e_minus4_lead': omission_detected,
        'first_noncancelling_order': first['first_noncancelling_paired_inverse_energy_order'],
        'leading_e_minus2_cancellation_is_not_a_tail_bound': True,
        'manufactured_even_parity_fixture': True,
        'production_C4_coefficient_control': False,
        'certificate_use': False,
    })


def smallest_representative_pair(root=None):
    """Smallest original |ell| pair: group 1, massless, |ell|=sqrt(5), mu=6, cutoff 160."""
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


def named_remainder_gaps():
    return MappingProxyType({
        'numerical_C_M': None,
        'numerical_C4': None,
        'vacuum_integrated_tail_160': None,
        'vacuum_integrated_tail_320': None,
        'thermal_coherent_source_remainder': 'bounded separately by Fermi occupation',
        'missing_for_C_M': MISSING_L0_A4,
        'missing_for_C4': (
            'pointwise e^{-3} N,beta densities on I from transported major A2, '
            'major A3 and n4; the leading E^{-2} cancellation is not that coefficient'),
        'upstream_affine_tail': MISSING_UPSTREAM_TAIL,
        'physical_local_gate': 'OPEN',
        'leading_e_minus2_cancellation_is_not_a_tail_bound': True,
    })


def current_history_remainder_pilot(root=None, cpu_limit=300.0):
    """Current 0b0e4ced history, smallest original pair, no full source campaign."""
    del cpu_limit
    root = Path(root) if root is not None else _root()
    history = json.loads((root / HISTORY_RECORD).read_text())
    family = LocalIncomingFamily(np.array(history['history']['coefficients']))
    identity = profile_identity(family)
    if identity != CURRENT_HISTORY_IDENTITY:
        raise ValueError('expected exact live history 0b0e4ced')
    pair = smallest_representative_pair(root)
    nodes = family.collocation_nodes(3)
    commutators = current_history_commutator_integrals(
        family, pair['absolute_angular'])
    majorant = conditional_m4_remainder_majorant(
        None, None, commutators['Bz_integral'], commutators['Bzz_integral'])
    thermal = fermi_thermal_tail_bound(
        pair['mass'], pair['absolute_angular'],
        pair['multiplicity_per_signed_family'],
        cutoffs=pair['archived_splits'])
    first = first_noncancelling_paired_coefficient()
    comparison = []
    for row in thermal['cutoffs']:
        comparison.append({
            'cutoff': row['cutoff'],
            'thermal_N': None,
            'thermal_beta': None,
            'thermal_N_display': row['N_display'],
            'thermal_beta_display': row['beta_display'],
            'thermal_N_decimal_upper': row['N_decimal_upper'],
            'thermal_beta_decimal_upper': row['beta_decimal_upper'],
            'thermal_log_N': row['log_N'],
            'thermal_log_beta': row['log_beta'],
            'vacuum_N': None,
            'vacuum_beta': None,
            'provisional_UV_allocation': PROVISIONAL_UV_ALLOCATION,
            'thermal_below_allocation': row['below_provisional_allocation'],
            'vacuum_below_allocation': None,
            'allocation_certified': False,
        })
    return MappingProxyType({
        'profile_identity': identity,
        'z': [float(value) for value in nodes],
        'pair': dict(pair),
        'imag_major_A2_vanishes_for_this_massless_pair': pair['mass'] == 0.0,
        'first_noncancelling_paired_order': first[
            'first_noncancelling_paired_inverse_energy_order'],
        'formulas': dict(first['formulas']),
        'commutators': dict(commutators),
        'majorant': {
            'numerical_C_M': majorant['numerical_C_M'],
            'missing_inputs': list(majorant['missing_inputs']),
            'Bz_integral': commutators['Bz_integral'],
            'Bzz_integral': commutators['Bzz_integral'],
        },
        'numerical_C4': None,
        'numerical_C_M': None,
        'thermal': {
            'cutoffs': [dict(row) for row in thermal['cutoffs']],
            'source_law': thermal['source_law'],
            'vacuum_envelope_used': False,
        },
        'vacuum_tail': {'160': None, '320': None},
        'provisional_comparison': comparison,
        'field_or_source_evolutions': 0,
        'physical_local_gate': PHYSICAL_LOCAL_GATE,
    })


def frozen_coefficient_records(root=None):
    root = Path(root) if root is not None else _root()
    rows = {}
    for path, expected in (
            (V2_RECORD, FROZEN_V2_SHA256),
            (V3_RECORD, None),
            (V4_RECORD, None)):
        digest = sha256((root / path).read_bytes()).hexdigest()
        if expected is not None and digest != expected:
            raise ValueError('frozen v2 UV coefficient record bytes changed')
        payload = json.loads((root / path).read_text())
        rows[path] = {
            'sha256': digest,
            'schema': payload['schema'],
            'numerical_C_M': payload.get('numerical_C_M'),
        }
    return MappingProxyType(rows)


def finite_history_uv_remainder_report():
    """Executable remainder summary: order, defect, conditional C_M, thermal."""
    identities = remainder_recurrence_identities()
    first = first_noncancelling_paired_coefficient()
    defect = truncated_envelope_L0_defect()
    manufactured = manufactured_defect_majorant_control()
    negative = negative_pair_and_omission_controls()
    return MappingProxyType({
        'schema': 'NSC-KS-FINITE-HISTORY-UV-REMAINDER-v1',
        'identities': identities,
        'first_noncancelling': first,
        'defect': defect,
        'manufactured': manufactured,
        'negative_controls': negative,
        'gaps': named_remainder_gaps(),
        'provisional_UV_allocation': PROVISIONAL_UV_ALLOCATION,
        'physical_local_gate': PHYSICAL_LOCAL_GATE,
        'certificate_from_leading_e_minus2_cancellation': False,
        'C_M_or_physical_gate_from_leading_cancellation': False,
        'new_action_term': False,
        'Gamma_rest_assigned': False,
        'existing_source_ad759_evaluation_changed': False,
        'v2_v4_records_rewritten': False,
    })


def _digest(root, relative):
    return sha256((Path(root) / relative).read_bytes()).hexdigest()


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


def _interval_lo(value):
    return mp.make_mpf(value._mpi_[0])


def _interval_hi(value):
    return mp.make_mpf(value._mpi_[1])


def _interval_upper_str(value):
    return mp.nstr(_interval_hi(value), 20)


def _exponential_moments_window(alpha, lower, upper=None):
    """Directed integrals of e^{-alpha E} and E e^{-alpha E} on [lower, upper)."""
    if mp.iv.mpf(alpha) <= 0 or mp.iv.mpf(lower) <= 0:
        raise ValueError('positive decay rate and lower cutoff required')
    lower = mp.iv.mpf(lower)
    exponential = mp.iv.exp(-alpha * lower)
    moment0 = exponential / alpha
    moment1 = exponential * (lower / alpha + 1 / alpha**2)
    if upper is None:
        return moment0, moment1
    upper = mp.iv.mpf(upper)
    if not (upper > lower):
        raise ValueError('window upper must exceed lower cutoff')
    exponential_hi = mp.iv.exp(-alpha * upper)
    moment0 -= exponential_hi / alpha
    moment1 -= exponential_hi * (upper / alpha + 1 / alpha**2)
    return moment0, moment1


def _thermal_n_beta_interval(
        mass, angular, multiplicity, lower, kappa, omega, upper=None):
    """Occupation N,beta interval on [lower, upper) or [lower, inf)."""
    mass = float(mass)
    angular = float(angular)
    multiplicity = _require_positive(multiplicity, 'multiplicity')
    kappa = _require_positive(kappa, 'surface gravity')
    omega = _require_positive(omega, 'Omega')
    lower = _require_positive(lower, 'cutoff')
    if mass < 0:
        raise ValueError('nonnegative mass required')
    if upper is not None:
        upper = _require_positive(upper, 'window upper')
    _, _, axial, _, _ = reference_chart(1.0)
    radius = float(np.sqrt(2.0))
    if axial <= 0:
        raise ValueError('positive Sigma axial required')
    envelope = abs(mass) + abs(angular) / radius
    moment0 = mp.iv.mpf(0)
    moment1 = mp.iv.mpf(0)
    for coefficient, alpha in (
            (2, mp.iv.pi / kappa),
            (1, 2 * mp.iv.pi / (omega * kappa))):
        m0, m1 = _exponential_moments_window(alpha, lower, upper)
        moment0 += coefficient * m0
        moment1 += coefficient * m1
    density = 2 * moment0 / (2 * mp.iv.pi)
    current = 2 * (moment1 + envelope * moment0) / (2 * mp.iv.pi)
    n_bound = multiplicity * (
        (abs(mass) + abs(angular) / radius) * density + current / axial)
    beta_bound = multiplicity * current
    return n_bound, beta_bound


def _interval_overlap_and_rest(window, rest, infinite):
    combo = window + rest
    overlap = (
        _interval_lo(combo) <= _interval_hi(infinite)
        and _interval_lo(infinite) <= _interval_hi(combo))
    infinite_above_window = _interval_hi(infinite) >= _interval_lo(window)
    rest_positive = _interval_hi(rest) > 0
    return overlap and infinite_above_window and rest_positive


def _quantity_binding_slice(bindings):
    return MappingProxyType({
        'state_law': bindings['state_law']['sha256'],
        'profile': bindings['profile']['identity'],
        'profile_record': bindings['profile']['sha256'],
        'preparation': bindings['preparation']['sha256'],
        'preparation_restart': bindings['preparation']['restart_sha256'],
        'preparation_cauchy': bindings['preparation']['cauchy_sha256'],
        'inventory': bindings['inventory']['payload_sha256'],
        'inventory_record': bindings['inventory']['record_sha256'],
        'inventory_cutoff_bridge': bindings['inventory']['cutoff_bridge_sha256'],
        'subtraction': bindings['subtraction']['sha256'],
    })


def _plain_binding(value):
    if isinstance(value, Mapping):
        return {str(key): _plain_binding(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain_binding(item) for item in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


def _quantity(name, *, formula, derivation, numerical_value, status,
              missing_primitive, bindings, **extra):
    payload = {
        'name': name,
        'formula': formula,
        'derivation': derivation,
        'numerical_value': numerical_value,
        'status': status,
        'missing_primitive': missing_primitive,
        'bindings': dict(bindings),
        'leading_e_minus2_cancellation_is_not_this_value': True,
        'e_minus3_first_allowed_order_is_not_a_tail_bound': True,
        'background_tail_certificate_transferred': False,
        'replaced_by_local_metric_jet': False,
        'fitted_decay': False,
    }
    payload.update(extra)
    return MappingProxyType(payload)


def changed_history_uv_bindings(root=None):
    """Authenticate state law, profile, preparation, inventory and subtraction."""
    root = Path(root) if root is not None else _root()
    history = json.loads((root / HISTORY_RECORD).read_text())
    family = LocalIncomingFamily(np.array(history['history']['coefficients']))
    identity = profile_identity(family)
    if (not CURRENT_HISTORY_IDENTITY.startswith(CURRENT_HISTORY_IDENTITY_PREFIX)
            or identity != CURRENT_HISTORY_IDENTITY):
        raise ValueError('expected exact live history 0b0e4ced')
    inventory = json.loads((root / INVENTORY_RECORD).read_text())
    payload = inventory['payload']
    payload_digest = _digest(root, payload['path'])
    if payload_digest != payload['sha256']:
        raise ValueError('original source inventory payload changed')
    vertices = confirm_raw_ks_vertices()
    scales = _state_scales(root)
    return MappingProxyType({
        'state_law': MappingProxyType({
            'statement': STATE_LAW,
            'owner': STATE_LAW_OWNER,
            'sha256': _digest(root, STATE_LAW_OWNER),
        }),
        'profile': MappingProxyType({
            'identity': identity,
            'record': HISTORY_RECORD,
            'sha256': _digest(root, HISTORY_RECORD),
            'class': history['history']['class'],
        }),
        'preparation': MappingProxyType({
            'owner': PREPARATION_OWNER,
            'sha256': _digest(root, PREPARATION_OWNER),
            'source_law': scales['source_law'],
            'restart_record': RESTART_RECORD,
            'restart_sha256': scales['record_sha256'],
            'cauchy_record': CAUCHY_RECORD,
            'cauchy_sha256': _digest(root, CAUCHY_RECORD),
            'surface_gravity': scales['surface_gravity'],
            'omega': scales['omega'],
            'physical_source_renormalized': False,
        }),
        'inventory': MappingProxyType({
            'record': INVENTORY_RECORD,
            'record_sha256': _digest(root, INVENTORY_RECORD),
            'payload_path': payload['path'],
            'payload_sha256': payload_digest,
            'cutoff_bridge_record': CUTOFF_BRIDGE_CONTROL,
            'cutoff_bridge_sha256': _digest(root, CUTOFF_BRIDGE_CONTROL),
        }),
        'subtraction': MappingProxyType({
            'owner': SUBTRACTION_OWNER,
            'sha256': _digest(root, SUBTRACTION_OWNER),
            'vertices': {
                'N_density': vertices['N_density'],
                'N_momentum': vertices['N_momentum'],
                'beta_density': vertices['beta_density'],
                'beta_momentum': vertices['beta_momentum'],
                'action_sign': vertices['action_sign'],
                'measure': vertices['measure'],
                'numerical_residuals': dict(vertices['numerical_residuals']),
            },
            'owned_uv_subtraction': (
                'raw KS vertices with action minus once and measure de/(2 pi)'),
            'background_tail_certificate_transferred': False,
        }),
    })


def authenticate_changed_history_uv_bindings(bindings, root=None):
    """Recompute the five binding hashes and reject a mismatched ledger."""
    bindings = _as_map(bindings, 'bindings')
    fresh = changed_history_uv_bindings(root)
    for key in ('state_law', 'profile', 'preparation', 'inventory', 'subtraction'):
        if key not in bindings:
            raise ValueError('missing UV binding: ' + key)
    if bindings['state_law']['sha256'] != fresh['state_law']['sha256']:
        raise ValueError('state-law binding hash mismatch')
    if bindings['state_law']['statement'] != STATE_LAW:
        raise ValueError('state-law statement changed')
    if bindings['profile']['identity'] != fresh['profile']['identity']:
        raise ValueError('profile identity mismatch')
    if bindings['profile']['sha256'] != fresh['profile']['sha256']:
        raise ValueError('profile binding hash mismatch')
    if bindings['preparation']['sha256'] != fresh['preparation']['sha256']:
        raise ValueError('preparation binding hash mismatch')
    if bindings['preparation']['restart_sha256'] != fresh['preparation']['restart_sha256']:
        raise ValueError('preparation restart hash mismatch')
    if bindings['inventory']['payload_sha256'] != fresh['inventory']['payload_sha256']:
        raise ValueError('inventory payload binding hash mismatch')
    if bindings['inventory']['record_sha256'] != fresh['inventory']['record_sha256']:
        raise ValueError('inventory record binding hash mismatch')
    if bindings['subtraction']['sha256'] != fresh['subtraction']['sha256']:
        raise ValueError('subtraction binding hash mismatch')
    if bindings['subtraction']['background_tail_certificate_transferred']:
        raise ValueError('background tail certificate must not transfer')
    if fresh['subtraction']['background_tail_certificate_transferred']:
        raise ValueError('background tail certificate must not transfer')
    if _plain_binding(bindings) != _plain_binding(fresh):
        raise ValueError('changed-history UV binding ledger differs from its owners')
    return fresh


def thermal_tail_direct_window_check(
        mass, angular, multiplicity, cutoffs=(160.0, 320.0), *,
        kappa=None, omega=None, dps=80):
    """Compare occupation tails to direct finite energy windows.

    Vacuum N,beta values stay None. A finite window is not a vacuum-tail
    certificate and not a C4 or C_M.
    """
    scales = _state_scales() if kappa is None or omega is None else None
    kappa = float(scales['surface_gravity'] if kappa is None else kappa)
    omega = float(scales['omega'] if omega is None else omega)
    rows = []
    with mp.workdps(dps):
        mp.iv.dps = dps
        for cutoff in cutoffs:
            cutoff = _require_positive(cutoff, 'cutoff')
            window_upper = 2.0 * cutoff
            infinite_n, infinite_beta = _thermal_n_beta_interval(
                mass, angular, multiplicity, cutoff, kappa, omega)
            window_n, window_beta = _thermal_n_beta_interval(
                mass, angular, multiplicity, cutoff, kappa, omega,
                upper=window_upper)
            rest_n, rest_beta = _thermal_n_beta_interval(
                mass, angular, multiplicity, window_upper, kappa, omega)
            n_ok = _interval_overlap_and_rest(window_n, rest_n, infinite_n)
            beta_ok = _interval_overlap_and_rest(
                window_beta, rest_beta, infinite_beta)
            n_display = _interval_upper_str(infinite_n)
            beta_display = _interval_upper_str(infinite_beta)
            window_n_display = _interval_upper_str(window_n)
            window_beta_display = _interval_upper_str(window_beta)
            rest_n_display = _interval_upper_str(rest_n)
            rest_beta_display = _interval_upper_str(rest_beta)
            if n_display == '0.0' or beta_display == '0.0':
                raise ArithmeticError('thermal tail rounded to zero')
            rows.append({
                'cutoff': float(cutoff),
                'window': [float(cutoff), float(window_upper)],
                'N': None,
                'beta': None,
                'vacuum_N': None,
                'vacuum_beta': None,
                'infinite_N_decimal_upper': n_display,
                'infinite_beta_decimal_upper': beta_display,
                'window_N_decimal_upper': window_n_display,
                'window_beta_decimal_upper': window_beta_display,
                'rest_N_decimal_upper': rest_n_display,
                'rest_beta_decimal_upper': rest_beta_display,
                'window_plus_rest_encloses_infinite_N': bool(n_ok),
                'window_plus_rest_encloses_infinite_beta': bool(beta_ok),
                'window_is_not_exact_zero': window_n_display != '0.0'
                    and window_beta_display != '0.0',
                'rest_is_not_exact_zero': rest_n_display != '0.0'
                    and rest_beta_display != '0.0',
                'underflow_is_not_exact_zero': True,
            })
    if not rows or not all(
            row['window_plus_rest_encloses_infinite_N']
            and row['window_plus_rest_encloses_infinite_beta']
            and row['window_is_not_exact_zero']
            and row['rest_is_not_exact_zero'] for row in rows):
        raise ArithmeticError('thermal tail failed the direct finite-window check')
    return MappingProxyType({
        'source_law': '2*exp(-pi*E/kappa)+exp(-2*pi*E/(Omega*kappa))',
        'surface_gravity': kappa,
        'omega': omega,
        'mass': float(mass),
        'angular': float(angular),
        'multiplicity_per_signed_family': float(multiplicity),
        'vacuum_envelope_used': False,
        'background_tail_certificate_transferred': False,
        'vacuum_integrated_tail_N_beta': None,
        'cutoffs': tuple(rows),
    })


def changed_history_uv_quantities(root=None, bindings=None):
    """Named changed-history UV quantities with formulas, nulls and bindings."""
    root = Path(root) if root is not None else _root()
    bindings = (
        changed_history_uv_bindings(root) if bindings is None
        else authenticate_changed_history_uv_bindings(bindings, root))
    slice_ = _quantity_binding_slice(bindings)
    identities = remainder_recurrence_identities()
    first = first_noncancelling_paired_coefficient()
    pair = smallest_representative_pair(root)
    family = LocalIncomingFamily(np.array(
        json.loads((root / HISTORY_RECORD).read_text())['history']['coefficients']))
    commutators = current_history_commutator_integrals(
        family, pair['absolute_angular'])
    majorant = conditional_m4_remainder_majorant(
        None, None, commutators['Bz_integral'], commutators['Bzz_integral'])
    massless = pair['mass'] == 0.0
    imag_rhs = 0.0 if massless else None
    formulas = dict(first['formulas'])
    quantities = {
        'A2': _quantity(
            'A2',
            formula=(
                'Im(T_s major A2)=-m ell/4 (s a^2 r_rho + r_z)/r^2; '
                'Re(T_s major A2) retains the transported major-A1 phase q; '
                'minor A2 is vacant L0 A1. Not a local metric jet'),
            derivation=(
                'coefficient_identities Im_T_s_major_A2_parent and '
                'remainder_recurrence Re_Ts_major_A2_depends_on_q / '
                'Im_Ts_major_A2_independent_of_q'),
            numerical_value=None,
            status='OPEN',
            missing_primitive=(
                'characteristic integrals of real major A2, and of imaginary '
                'major A2 when mass is nonzero'),
            bindings=slice_,
            imag_major_transport_rhs_this_massless_pair=imag_rhs,
            imag_major_transport_rhs_vanishes_by_m_0_identity=massless,
            imag_major_A2_numerical_value=None,
            upstream_imag_major_A2_owned=False,
            real_major_A2_numerical_value=None,
            minor_A2_numerical_value=None,
            A3_locally_algebraic_in_metric_jets=identities[
                'A3_locally_algebraic_in_metric_jets'],
        ),
        'A3': _quantity(
            'A3',
            formula=(
                'minor A3 = a^2/(2 i) Pi_{-s} L0 A2 on the transported major '
                'A2; major A3 from T_s. Not a local metric jet'),
            derivation=(
                'remainder_recurrence minor_Ajp1_depends_on_major_Aj and '
                'coefficient_identities minor_A3_depends_on_major_A2'),
            numerical_value=None,
            status='OPEN',
            missing_primitive=(
                'transported major A2, then vacant L0 A2 and major-A3 transport'),
            bindings=slice_,
            depends_on_major_A2=True,
            A3_locally_algebraic_in_metric_jets=False,
        ),
        'n4': _quantity(
            'n4',
            formula=formulas['n4_transport'],
            derivation=(
                'first_noncancelling continuity_e_minus4_retains_minor_A3; '
                'upstream n4=0 does not remove minor A3'),
            numerical_value=None,
            status='OPEN',
            missing_primitive='minor A2, minor A3 and the n4 characteristic integral',
            bindings=slice_,
            retains_minor_A3=True,
            n4_eliminated_without_minor_A3=first['n4_eliminated_without_minor_A3'],
        ),
        'C4': _quantity(
            'C4',
            formula=formulas['C4'],
            derivation=(
                'equal-mu original pair of the e^{-3} action densities on I; '
                'J_3 and action_N3/action_beta3 retain transported A2, A3 and n4. '
                'The E^{-3} first allowed order is not this coefficient'),
            numerical_value=None,
            status='OPEN',
            missing_primitive=named_remainder_gaps()['missing_for_C4'],
            bindings=slice_,
            first_noncancelling_paired_inverse_energy_order=first[
                'first_noncancelling_paired_inverse_energy_order'],
            complete_e_minus3_coefficient_proven_nonzero=first[
                'complete_e_minus3_coefficient_proven_nonzero'],
            C4_or_C_M_bounded=first['C4_or_C_M_bounded'],
            history_minus_reference_coefficient_formed=False,
            definition_scope=(
                'OPEN: full-envelope versus history-minus-reference e^{-3} '
                'contraction is not yet fixed by a production owner'),
        ),
        'C_M': _quantity(
            'C_M',
            formula=(
                '||X-X^{(4)}||_H2 <= C_M e^{-4} after propagate_h2 of the '
                'upstream H2 remainder, current-history L0 A4 H2 integrals, '
                'and owned B_z, B_zz'),
            derivation=(
                'conditional_m4_remainder_majorant on the declared slab; '
                'Bz/Bzz from current_history_commutator_integrals. Smoothness '
                'gives a finite constant for a fixed admissible history, not a '
                'uniform class bound and not a frozen-background tail'),
            numerical_value=None,
            status='OPEN',
            missing_primitive=tuple(majorant['missing_inputs']),
            bindings=slice_,
            C_M_h2=None,
            Bz_integral=commutators['Bz_integral'],
            Bzz_integral=commutators['Bzz_integral'],
            certified_L0_A4=commutators['certified_L0_A4'],
            applicability=(
                'current-history slab of 0b0e4ced with owned radius commutators; '
                'C_M remains unbounded until L0 A4 and upstream H2 inputs exist'),
        ),
        'current_history_L0_A4_H2_integrals': _quantity(
            'current_history_L0_A4_H2_integrals',
            formula='H2 integrals of L0 A4 on the declared current-history slab',
            derivation=(
                'truncated_envelope_L0_defect: formal defect is e^{-4} L0 A4 '
                'conditional on a constructed M=4 recurrence; the production '
                'A3/A4 coefficients are absent'),
            numerical_value=None,
            status='OPEN',
            missing_primitive=MISSING_L0_A4,
            bindings=slice_,
            M=REMAINDER_ORDER,
            defect=identities['defect'],
        ),
        'upstream_higher_order_H2_remainder': _quantity(
            'upstream_higher_order_H2_remainder',
            formula=(
                'H2 bound on e^4 (X-X^{(4)}) at rho_up for inverse-energy '
                'orders greater than 4, same mathematical vacuum column'),
            derivation=(
                'conditional_m4_remainder_majorant initial_h2; the frozen '
                'background tail certificate is not a proved bridge'),
            numerical_value=None,
            status='OPEN',
            missing_primitive=MISSING_UPSTREAM_TAIL,
            bindings=slice_,
            frozen_background_tail_copied=False,
        ),
        'transported_majors': _quantity(
            'transported_majors',
            formula=(
                'characteristic integrals of major A2 (real and imaginary) and '
                'the subsequent major-A3 transport; only the local '
                'Im(T_s major A2) right-hand side vanishes for m=0'),
            derivation=(
                'remainder_recurrence transport defect Pi_s L0 A_j and the '
                'massless vanishing identity Im_T_s_major_A2_vanishes_for_m_0'),
            numerical_value=None,
            status='OPEN',
            missing_primitive=(
                'real major-A2 characteristic integral, the upstream imaginary '
                'major-A2 datum, and major-A3 transport; the massless identity '
                'sets only the local imaginary transport right-hand side to zero'),
            bindings=slice_,
            imag_major_A2_this_massless_pair=None,
            imag_major_A2_characteristic_integral_this_massless_pair=None,
            imag_major_A2_transport_rhs_this_massless_pair=imag_rhs,
            real_major_A2=None,
            major_A3=None,
            characteristic_integrals_owned=False,
            imag_major_A2_vanishes_for_this_massless_pair=massless,
        ),
        'vacuum_N_beta_tail': _quantity(
            'vacuum_N_beta_tail',
            formula=(
                'occupied vacuum-envelope N,beta tail above the archived '
                '160/320 splits after C4 and C_M; not the Fermi occupation'),
            derivation=(
                'named_remainder_gaps vacuum_integrated_tail_160/320 remain '
                'missing; thermal occupation is bounded separately'),
            numerical_value=None,
            status='OPEN',
            missing_primitive=(
                'C4, C_M, transported majors and L0 A4 integrals on I; the '
                'E^{-2} cancellation is not this tail'),
            bindings=slice_,
            N=None,
            beta=None,
            vacuum_integrated_tail_160=None,
            vacuum_integrated_tail_320=None,
            thermal_occupation_bounded_separately=True,
        ),
    }
    if tuple(quantities) != QUANTITY_NAMES:
        raise ArithmeticError('changed-history UV quantity names drifted')
    return MappingProxyType(quantities)


def validate_changed_history_uv_successor(report, root=None):
    """Reject invented constants, false PASS/NON_EXISTENCE, or unbound quantities."""
    report = _as_map(report, 'successor report')
    if report.get('schema') == PILOT_SCHEMA:
        raise ValueError('v1 pilot schema is immutable and is not the v2 successor')
    if report.get('schema') != SUCCESSOR_SCHEMA:
        raise ValueError('unexpected changed-history UV successor schema')
    status = report.get('status', SUCCESSOR_STATUS)
    if not isinstance(status, str) or not status.startswith('OPEN'):
        raise ValueError('OPEN calculation may not be reported as PASS or NON_EXISTENCE')
    if 'NON-EXISTENCE' in status or status.startswith('PASS'):
        raise ValueError('OPEN calculation may not be reported as PASS or NON_EXISTENCE')
    if report.get('physical_EXISTENCE_certificate') or report.get(
            'physical_NONEXISTENCE_certificate'):
        raise ValueError('OPEN calculation may not be reported as PASS or NON_EXISTENCE')
    if report.get('physical_local_gate') == 'PASS':
        raise ValueError('OPEN calculation may not be reported as PASS or NON_EXISTENCE')
    if report.get('certificate_from_leading_e_minus2_cancellation'):
        raise ValueError('E^{-2} cancellation is not a tail bound')
    if report.get('C_M_or_physical_gate_from_leading_cancellation'):
        raise ValueError('leading cancellation does not supply C_M')
    if report.get('first_noncancelling_paired_order') != FIRST_NONCANCELLING_PAIRED_ORDER:
        raise ValueError('E^{-3} first allowed paired order was not preserved')
    if report.get('leading_paired_e_minus2_cancels_for_equal_mu') is False:
        raise ValueError('E^{-2} equal-mu cancellation was not preserved')
    if report.get('field_or_source_evolutions') not in (0, None):
        raise ValueError('successor control may not launch field or source evolutions')
    if report.get('v1_pilot_rewritten'):
        raise ValueError('historical v1 pilot bytes are immutable')
    bindings = _as_map(report.get('bindings'), 'bindings')
    authenticate_changed_history_uv_bindings(bindings, root)
    quantities = _as_map(report.get('quantities'), 'quantities')
    if set(quantities) != set(QUANTITY_NAMES):
        raise ValueError('changed-history UV quantity names drifted')
    slice_ = _quantity_binding_slice(bindings)
    for name in QUANTITY_NAMES:
        item = _as_map(quantities[name], name)
        if item.get('name') != name:
            raise ValueError('quantity name mismatch: ' + name)
        if item.get('status') != 'OPEN':
            raise ValueError('quantity ' + name + ' must remain OPEN')
        if not item.get('leading_e_minus2_cancellation_is_not_this_value'):
            raise ValueError('E^{-2} cancellation is not ' + name)
        if not item.get('e_minus3_first_allowed_order_is_not_a_tail_bound'):
            raise ValueError('E^{-3} first allowed order is not a tail bound')
        if item.get('background_tail_certificate_transferred'):
            raise ValueError('background tail certificate must not transfer onto ' + name)
        if item.get('replaced_by_local_metric_jet'):
            raise ValueError(name + ' may not be replaced by a local metric jet')
        if item.get('fitted_decay'):
            raise ValueError(name + ' may not use a fitted decay')
        bound = item.get('bindings')
        if _as_map(bound, name + ' bindings') != dict(slice_):
            raise ValueError(name + ' is not bound to state law/profile/preparation/inventory/subtraction')
        if name in NONE_REQUIRED_QUANTITIES:
            _require_missing(item.get('numerical_value'), name)
    vacuum = quantities['vacuum_N_beta_tail']
    for key in ('N', 'beta', 'vacuum_integrated_tail_160', 'vacuum_integrated_tail_320'):
        _require_missing(vacuum.get(key), 'vacuum ' + key)
    _require_missing(quantities['C4']['numerical_value'], 'C4')
    _require_missing(quantities['C_M']['numerical_value'], 'C_M')
    _require_missing(quantities['C_M'].get('C_M_h2'), 'C_M_h2')
    _require_missing(
        quantities['A2'].get('real_major_A2_numerical_value'), 'real major A2')
    _require_missing(
        quantities['A2'].get('imag_major_A2_numerical_value'), 'imag major A2')
    _require_missing(quantities['transported_majors'].get('real_major_A2'), 'real major A2')
    _require_missing(quantities['transported_majors'].get('major_A3'), 'major A3')
    _require_missing(quantities['transported_majors'].get(
        'imag_major_A2_this_massless_pair'), 'imag major A2')
    _require_missing(quantities['transported_majors'].get(
        'imag_major_A2_characteristic_integral_this_massless_pair'),
        'imag major A2 characteristic value')
    if quantities['A2'].get('imag_major_transport_rhs_this_massless_pair') not in (0.0, 0):
        raise ValueError('massless Im(T_s major A2) identity was not preserved')
    if not quantities['A2'].get('imag_major_transport_rhs_vanishes_by_m_0_identity'):
        raise ValueError('massless Im(T_s major A2) identity was not preserved')
    if quantities['A3'].get('A3_locally_algebraic_in_metric_jets'):
        raise ValueError('A3 may not be replaced by a local metric jet')
    if quantities['upstream_higher_order_H2_remainder'].get('frozen_background_tail_copied'):
        raise ValueError('frozen-background tail certificate copied without a proved bridge')
    if quantities['C4'].get('history_minus_reference_coefficient_formed'):
        raise ValueError('production history-minus-reference C4 is not formed')
    identities_report = _as_map(report.get('identities'), 'identity summary')
    if (not identities_report.get('dummy_symbol_L0_A_j')
            or not identities_report.get('formal_symbolic_telescope')
            or identities_report.get('production_A3_A4_constructed')):
        raise ValueError('formal M=4 telescope may not claim production A3/A4')
    defect_report = _as_map(report.get('defect'), 'defect summary')
    if (not defect_report.get('dummy_symbol_L0_A_j')
            or not defect_report.get('formal_symbolic_telescope')
            or defect_report.get('production_A3_A4_constructed')):
        raise ValueError('formal M=4 defect may not claim production A3/A4')
    manufactured = _as_map(report.get('manufactured'), 'manufactured control')
    if (manufactured.get('production_L0_A4_control')
            or manufactured.get('certificate_use')
            or manufactured.get('current_history_C_M') is not None):
        raise ValueError('manufactured fixture cannot supply current-history C_M')
    negative = _as_map(report.get('negative_controls'), 'negative controls')
    if (negative.get('production_C4_coefficient_control')
            or negative.get('certificate_use')):
        raise ValueError('manufactured parity fixture cannot supply C4')
    thermal = _as_map(report.get('thermal'), 'thermal occupation')
    if thermal.get('vacuum_envelope_used') or thermal.get(
            'background_tail_certificate_transferred'):
        raise ValueError('thermal occupation cannot become a vacuum-tail bound')
    for row in thermal.get('cutoffs', ()):
        row = _as_map(row, 'thermal cutoff')
        _require_missing(row.get('N'), 'thermal N float')
        _require_missing(row.get('beta'), 'thermal beta float')
        if (row.get('N_decimal_upper') in (None, '0.0')
                or row.get('beta_decimal_upper') in (None, '0.0')):
            raise ValueError('thermal occupation upper must remain directed and nonzero')
    windows = report.get('direct_windows')
    if windows is not None:
        windows = _as_map(windows, 'direct_windows')
        _require_missing(windows.get('vacuum_integrated_tail_N_beta'), 'vacuum tail')
        if windows.get('background_tail_certificate_transferred'):
            raise ValueError('background tail certificate must not transfer')
        for row in windows.get('cutoffs', ()):
            row = _as_map(row, 'direct window')
            _require_missing(row.get('vacuum_N'), 'vacuum N')
            _require_missing(row.get('vacuum_beta'), 'vacuum beta')
            _require_missing(row.get('N'), 'thermal N float')
            _require_missing(row.get('beta'), 'thermal beta float')
            if not row.get('window_plus_rest_encloses_infinite_N'):
                raise ValueError('thermal N failed the direct finite-window check')
            if not row.get('window_plus_rest_encloses_infinite_beta'):
                raise ValueError('thermal beta failed the direct finite-window check')
    return True


def changed_history_uv_successor_report(root=None):
    """Small authenticated control: represent, bind, keep missing values None."""
    root = Path(root) if root is not None else _root()
    bindings = changed_history_uv_bindings(root)
    quantities = changed_history_uv_quantities(root, bindings)
    first = first_noncancelling_paired_coefficient()
    identities = remainder_recurrence_identities()
    defect = truncated_envelope_L0_defect()
    manufactured = manufactured_defect_majorant_control()
    negative = negative_pair_and_omission_controls()
    pair = smallest_representative_pair(root)
    thermal = fermi_thermal_tail_bound(
        pair['mass'], pair['absolute_angular'],
        pair['multiplicity_per_signed_family'],
        cutoffs=pair['archived_splits'])
    windows = thermal_tail_direct_window_check(
        pair['mass'], pair['absolute_angular'],
        pair['multiplicity_per_signed_family'],
        cutoffs=pair['archived_splits'])
    family = LocalIncomingFamily(np.array(
        json.loads((root / HISTORY_RECORD).read_text())['history']['coefficients']))
    commutators = current_history_commutator_integrals(
        family, pair['absolute_angular'])
    report = MappingProxyType({
        'schema': SUCCESSOR_SCHEMA,
        'status': SUCCESSOR_STATUS,
        'profile_identity': bindings['profile']['identity'],
        'pair': dict(pair),
        'bindings': bindings,
        'quantities': quantities,
        'first_noncancelling_paired_order': first[
            'first_noncancelling_paired_inverse_energy_order'],
        'leading_paired_e_minus2_cancels_for_equal_mu': first[
            'leading_paired_e_minus2_cancels_for_equal_mu'],
        'e_minus3_even_in_ell': first['e_minus3_even_in_ell'],
        'e_minus3_first_allowed_order_is_not_a_tail_bound': True,
        'leading_e_minus2_cancellation_is_not_a_tail_bound': True,
        'identities': MappingProxyType({
            'expansion_order': identities['expansion_order'],
            'defect': identities['defect'],
            'A3_locally_algebraic_in_metric_jets': identities[
                'A3_locally_algebraic_in_metric_jets'],
            'major_A2_replaced_by_local_metric_jet': identities[
                'major_A2_replaced_by_local_metric_jet'],
            'dummy_symbol_L0_A_j': identities['dummy_symbol_L0_A_j'],
            'formal_symbolic_telescope': identities['formal_symbolic_telescope'],
            'production_A3_A4_constructed': identities[
                'production_A3_A4_constructed'],
            'telescope': identities['telescope'],
        }),
        'defect': MappingProxyType({
            'M': defect['M'],
            'defect': defect['defect'],
            'dummy_symbol_L0_A_j': defect['dummy_symbol_L0_A_j'],
            'formal_symbolic_telescope': defect['formal_symbolic_telescope'],
            'production_A3_A4_constructed': defect[
                'production_A3_A4_constructed'],
        }),
        'commutators': dict(commutators),
        'manufactured': manufactured,
        'negative_controls': negative,
        'thermal': {
            'cutoffs': [dict(row) for row in thermal['cutoffs']],
            'source_law': thermal['source_law'],
            'vacuum_envelope_used': thermal['vacuum_envelope_used'],
            'background_tail_certificate_transferred': False,
        },
        'direct_windows': windows,
        'numerical_C4': None,
        'numerical_C_M': None,
        'vacuum_integrated_tail_N_beta': None,
        'field_or_source_evolutions': 0,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'physical_local_gate': 'OPEN',
        'certificate_from_leading_e_minus2_cancellation': False,
        'C_M_or_physical_gate_from_leading_cancellation': False,
        'new_action_term': False,
        'Gamma_rest_assigned': False,
        'existing_source_ad759_evaluation_changed': False,
        'v2_v4_records_rewritten': False,
        'v1_pilot_schema': PILOT_SCHEMA,
        'v1_pilot_record': PILOT_RECORD,
        'v1_pilot_sha256': _digest(root, PILOT_RECORD),
        'v1_pilot_rewritten': False,
        'provisional_UV_allocation': PROVISIONAL_UV_ALLOCATION,
    })
    validate_changed_history_uv_successor(report, root)
    return report
