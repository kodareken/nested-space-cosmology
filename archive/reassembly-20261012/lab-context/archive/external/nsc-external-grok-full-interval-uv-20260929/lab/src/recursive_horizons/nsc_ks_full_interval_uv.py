"""One full-interval attempt at the history-minus-reference L0 A4 remainder.

Calls the existing characteristic jet on the whole interval
I = S(1)+[0.12, 0.18] of history 0b0e4ced. That call does not enclose I.
The phase-fixed upstream column is the stored vacuum spinor expansion, not
the zero solution of a vanishing imaginary transport right-hand side.

No numerical C_M is returned. A missing integral is not replaced by zero,
and an absolute jet sum is not used as the subtracted L0 A4 bound.
"""
from fractions import Fraction
from hashlib import sha256
import json
import resource
import time
from pathlib import Path

import numpy as np
from flint import arb, ctx

from .nsc_ks_ball_geometry import SubdivisionNeeded
from .nsc_ks_cubic_inventory import channel_moments, load_inventory
from .nsc_ks_cubic_uv_enclosure import CubicGeometry, characteristic_distance
from .nsc_ks_current_uv_transport import CURRENT_HISTORY_IDENTITY, bound_record
from .nsc_ks_finite_history_uv_coefficients import ARCHIVE_RHO_UP
from .nsc_ks_finite_history_uv_remainder import (
    HISTORY_RECORD,
    PROVISIONAL_UV_ALLOCATION,
    RESTART_RECORD,
    conditional_m4_remainder_majorant,
    current_history_commutator_integrals,
    manufactured_defect_majorant_control,
    smallest_representative_pair,
)
from .nsc_vacuum_source_remainder import VacuumSourceExpansion
from .nsc_vacuum_spinor_remainder import (
    spinor_remainder,
    truncated_upstream_envelope,
)
from .nsc_ks_profile_identity import profile_identity
from .nsc_ks_residual_error import pointwise_field_bounds
from .nsc_ks_uv_minor_a3 import minor_a3_identities, minor_a3_interval_majorants
from .nsc_local_incoming_family import LocalIncomingFamily


SCHEMA = 'NSC-KS-FULL-INTERVAL-UV-v1'
SPINOR_RECORD = 'results/development/nsc-vacuum-spinor-remainder-v1.json'
GROUP = 1
CUTOFF = 160
SAMPLE_RHO = arb(102) / arb(100)
UNCONSTRUCTED_TERM = (
    'Pi_{-s} L0 A4 = (2 i / a^2) Pi_{-s} A5, with '
    'A5_minor = a^2/(2 i) T_{-s}(A4_minor) - (a/2) V_{-s,s} A4_major, '
    'and the H2 integrals of its axial derivatives 0, 1 and 2 along both '
    'characteristics through the preimage of the whole interval I'
)


def _root():
    return Path(__file__).resolve().parents[2]


def _rss_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)


def _plain(ball):
    return {
        'contains_zero': bool(ball.contains(0)),
        'strictly_negative': bool(ball < 0),
        'strictly_positive': bool(ball > 0),
        'midpoint': float(ball.mid()),
    }


def _family(root):
    history = json_history(root)
    family = LocalIncomingFamily(np.array(history['history']['coefficients']))
    if profile_identity(family) != CURRENT_HISTORY_IDENTITY:
        raise ValueError('expected exact live history 0b0e4ced')
    return family


def json_history(root):
    return json.loads((root / HISTORY_RECORD).read_text())


def _dyadic_upper(ball):
    """Outward dyadic upper bound of a real ball as an exact fraction."""
    mantissa, exponent = ball.upper().man_exp()
    return Fraction(int(mantissa)) * Fraction(2)**int(exponent)


def _spinor_record(root):
    payload = (root / SPINOR_RECORD).read_bytes()
    record = json.loads(payload)
    if record.get('current_history_L0_A4_H2_integrals') is not None:
        raise ValueError('spinor record must not already claim L0 A4 integrals')
    if record.get('changed_history_C_M') is not None:
        raise ValueError('spinor record must not already claim changed-history C_M')
    if record.get('finite_occupation_source_error_included') is not False:
        raise ValueError('thermal occupation must stay out of the vacuum spinor record')
    return record, sha256(payload).hexdigest()


def _horizon_rho(root):
    restart = json.loads((root / RESTART_RECORD).read_text())
    hint = restart['scattering_provenance']['config']['horizon_rho']
    if not isinstance(hint, float) or not hint > 1:
        raise ValueError('owned horizon radius required')
    return hint


def phase_fixed_column(angular_sign, horizon_rho, cutoff=CUTOFF):
    """Live order-5 upstream spinor. The chart fixes the phase; major A2 is not zero."""
    if angular_sign not in (-1, 1):
        raise ValueError('explicit angular sign required')
    with ctx.workprec(128):
        model = VacuumSourceExpansion(
            horizon_rho, 0, arb(angular_sign) * arb(5).sqrt(),
            order=5, bits=128, metric_terms=16)
        result = spinor_remainder(model, ARCHIVE_RHO_UP, cutoff, cells=16)
        envelope = truncated_upstream_envelope(result, 4, 1)
        major, minor_re, minor_im = result['coefficients']
        norm = major[2] + (minor_re[1]**2 + minor_im[1]**2) / 2
        tail = envelope['pointwise_scaled_error_upper']
        h2 = envelope['scaled_envelope_H2_triple']
        if (not major[0].contains(1) or not major[1].contains(0)
                or not norm.contains(0) or major[2].contains(0) or not tail > 0):
            raise ArithmeticError('phase-fixed upstream column failed its normalization')
        if h2[1] != 0 or h2[2] != 0 or envelope['carrier_derivatives_included']:
            raise ArithmeticError('upstream envelope derivatives are not the z-independent tail')
        if envelope['changed_history_remainder'] is not None:
            raise ArithmeticError('spinor envelope claimed a changed-history remainder')
        column = {
            'major_A2': major[2],
            'major_A4': major[4],
            'normalization_contains_zero': bool(norm.contains(0)),
            'degree4_pointwise_tail_upper': tail,
            'order': result['order'],
            'cells': 16,
        }
    return column


def _exact_inner_ball(family):
    """Midpoint-radius ball whose requested radius is the binary axial inner."""
    center = arb(float(family.center))
    inner = arb(float(family.axial_inner))
    ball = arb(center.man_exp(), inner.man_exp())
    return center, inner, arb(float(family.axial_outer)), ball


def attempt_full_interval_jets(family, angular):
    """Existing jet API on the closed interval I, plus one proper subinterval."""
    center, inner, outer, ball = _exact_inner_ball(family)
    radius_exceeds = not (abs(ball - center).upper() <= inner)
    if not radius_exceeds:
        raise ArithmeticError('dyadic radius of I was not strictly outside the flat test')
    with ctx.workprec(96):
        model = CubicGeometry(family, angular, bits=96)
        try:
            model.jets(arb(1), ball)
        except SubdivisionNeeded as error:
            failure = str(error)
        else:
            raise ArithmeticError('full interval I was enclosed; the domain failure changed')
        if 'straddles the plateau center' not in failure:
            raise ArithmeticError('unexpected full-interval refusal: ' + failure)
        gap = arb(2)**-26
        sub = (center + gap).union(center + inner - gap)
        sub_jets = model.jets(arb(1), sub)
        if not sub_jets['b'].is_finite():
            raise ArithmeticError('proper subinterval jet was not finite')
    return {
        'interval_I': 'S(1)+[0.12,0.18]',
        'center': float(family.center),
        'axial_inner': float(family.axial_inner),
        'axial_outer': float(family.axial_outer),
        'dyadic_radius_upper': float(ball.rad().upper()),
        'dyadic_radius_exceeds_axial_inner': True,
        'domain_failure': failure,
        'full_interval_jets_enclosed': False,
        'proper_subinterval_enclosed': True,
        'proper_subinterval_is_full_I': False,
        'subinterval_b_z_contains_zero_on_sigma': bool(sub_jets['b_z'].contains(0)),
        'center_ball': center,
        'inner_ball': inner,
        'outer_ball': outer,
    }


def _preimage_exits(center, inner, outer):
    distance = characteristic_distance(arb(ARCHIVE_RHO_UP))
    margin = outer - inner
    if not distance > margin:
        raise ArithmeticError('characteristic shift no longer leaves the axial support')
    exits = {}
    for sign in (1, -1):
        # Image of the closed inner interval under z -> z - sign * D.
        left = center - inner - arb(sign) * distance
        right = center + inner - arb(sign) * distance
        exits[sign] = bool(left < center - outer or right > center + outer)
        if not exits[sign]:
            raise ArithmeticError('sign %+d preimage of I stayed inside the support' % sign)
    return {
        'characteristic_distance_lower': float(distance.lower()),
        'characteristic_distance_upper': float(distance.upper()),
        'support_margin': float(margin.upper()) if hasattr(margin, 'upper') else float(margin),
        'both_characteristics_exit_support': True,
        'sign_plus_exits': exits[1],
        'sign_minus_exits': exits[-1],
    }


def _point_reentry(family, angular, major_a2):
    """Nonzero R_up * delta(1/r) at one slab point on a ray that meets I."""
    zero = LocalIncomingFamily(np.zeros((2, 8)))
    with ctx.workprec(96):
        live = CubicGeometry(family, angular, bits=96)
        reference = CubicGeometry(zero, angular, bits=96)
        rho = SAMPLE_RHO
        if not arb(1) < rho < arb(ARCHIVE_RHO_UP):
            raise ValueError('sample rho must lie in the open preparation slab')
        distance = characteristic_distance(rho)
        center = arb(float(family.center))
        inner = arb(float(family.axial_inner))
        if not distance < inner:
            raise ArithmeticError('sample ray from the center of I leaves I')
        z = center
        live_jets = live.jets(rho, z)
        ref_jets = reference.jets(rho, z)
        delta_b = live_jets['b'] - ref_jets['b']
        if delta_b.contains(0) or not live_jets['a'].overlaps(ref_jets['a']):
            raise ArithmeticError('sample radius difference is not separated')
        ell = arb(angular)
        delta_inv_r = 2 * delta_b / (live_jets['a'] * ell)
        term = major_a2 * delta_inv_r
        if term.contains(0):
            raise ArithmeticError('R_up * delta(1/r) contains zero at the sample')
    return {
        'sample_rho': '102/100',
        'sample_is_full_interval': False,
        'sample_ray_meets_I': True,
        'delta_b': _plain(delta_b),
        'delta_inv_r': _plain(delta_inv_r),
        'R_up_times_delta_inv_r': _plain(term),
        'pointwise_reentry_contains_zero': False,
        'zero_rhs_value_at_sample': 0,
        'zero_rhs_removes_separated_term': True,
    }


def _reentry_majorant(family, pair, major_a2):
    majorants = minor_a3_interval_majorants(
        family, pair['mass'], pair['absolute_angular'])
    if majorants['complete_history_minus_reference_minor_A3_abs_upper'] is not None:
        raise ArithmeticError('complete minor A3 was not the missing quantity')
    kernel = majorants['upstream_major_A2_difference_kernel_abs_upper']
    if kernel <= 0:
        raise ArithmeticError('live upstream-A2 kernel vanished')
    # major A2 is negative; its absolute upper is the lower endpoint's magnitude.
    if major_a2.upper() >= 0:
        raise ArithmeticError('phase-fixed major A2 upper is not negative')
    absolute = _dyadic_upper(-major_a2.lower())
    product = kernel * absolute
    if product <= 0:
        raise ArithmeticError('reentry majorant is not positive')
    zero_family = LocalIncomingFamily(np.zeros((2, 8)))
    zeroed = minor_a3_interval_majorants(
        zero_family, pair['mass'], pair['absolute_angular'])
    if zeroed['upstream_major_A2_difference_kernel_abs_upper'] != 0:
        raise ArithmeticError('zero deformation must kill the kernel')
    if zeroed['complete_history_minus_reference_minor_A3_abs_upper'] != 0:
        raise ArithmeticError('zero deformation difference must be exact zero')
    return {
        'kernel': bound_record(kernel),
        'major_A2_absolute_upper': bound_record(absolute),
        'reentry_majorant': bound_record(product),
        'zero_deformation_kernel': 0,
        'zero_rhs_reentry_majorant': 0,
        'zero_rhs_drops_positive_reentry_majorant': True,
        'complete_minor_A3': None,
        'loose_absolute_jet_used_as_difference': False,
    }


def _carrier(tail):
    """E * F carrier term and the same 1/(2 Ec^2) integral weight.

    Period 1 makes the L2 factor exact and is not a physical axial domain.
    The spinor record leaves that domain unbound.
    """
    constant = _dyadic_upper(tail)
    energy = Fraction(CUTOFF)
    field = pointwise_field_bounds((constant, 0, 0), Fraction(1), energy)
    if field['F'] != constant or field['F_z'] != energy * field['F']:
        raise ArithmeticError('carrier restoration is not E times F')
    scale = arb(constant.numerator) / arb(constant.denominator)
    integral = scale / (4 * arb.pi() * arb(CUTOFF)**2)
    if not integral > arb(PROVISIONAL_UV_ALLOCATION):
        raise ArithmeticError('absolute carrier integral did not exceed the provisional allocation')
    same = scale / (2 * arb(CUTOFF)**2)
    if not (integral * 2 * arb.pi() - same).contains(0):
        raise ArithmeticError('carrier integral is not the leading 1/(2 Ec^2) weight')
    return {
        'F_z_equals_energy_times_F': True,
        'placeholder_period_used_as_physical_domain': False,
        'absolute_common_carrier_integral_upper': float(integral.upper()),
        'same_leading_cutoff_weight': '1/(2*Ec^2)',
        'provisional_allocation': PROVISIONAL_UV_ALLOCATION,
        'absolute_common_tail_exceeds_allocation': True,
        'exceeding_allocation_is_not_nonexistence': True,
        'this_integral_is_history_minus_reference': False,
    }


def _leading_weight(root):
    """Owned fold. Each family is divided by its own 2 Ec^2; 320 is not 160."""
    moments = channel_moments(load_inventory(root))
    if moments['energy_folding_factor'] != 1 or moments['C_M'] is not None:
        raise ValueError('moment fold changed or already carries C_M')
    if moments['angular_families'] != 60:
        raise ValueError('original 60 angular families required')
    return {
        'angular_families': moments['angular_families'],
        'leading_cutoff_factor': 'per family 1/(2*Ec^2), Ec in {160, 320}',
        'leading_moment_uppers': [float(item.upper()) for item in moments['leading']],
        'moment_fold_supplies_C_M': False,
    }


def full_interval_uv_report(root=None):
    """Actual 0b0e4ced attempt. C_M stays null when the full-I integral is absent."""
    started = time.process_time()
    root = Path(root) if root is not None else _root()
    family = _family(root)
    pair = smallest_representative_pair(root)
    if pair['group'] != GROUP or pair['cutoff'] != CUTOFF or pair['mass'] != 0:
        raise ValueError('smallest original pair is no longer massless group 1 at 160')
    record, digest = _spinor_record(root)
    horizon = _horizon_rho(root)
    positive = phase_fixed_column(1, horizon)
    negative = phase_fixed_column(-1, horizon)
    domain = attempt_full_interval_jets(family, pair['absolute_angular'])
    preimage = _preimage_exits(
        domain.pop('center_ball'), domain.pop('inner_ball'), domain.pop('outer_ball'))
    point = _point_reentry(family, pair['absolute_angular'], positive['major_A2'])
    reentry = _reentry_majorant(family, pair, positive['major_A2'])
    carrier = _carrier(positive['degree4_pointwise_tail_upper'])
    moments = _leading_weight(root)
    commutators = current_history_commutator_integrals(
        family, pair['absolute_angular'])
    if commutators['certified_L0_A4']:
        raise ValueError('commutator integrals are not L0 A4')
    withheld = conditional_m4_remainder_majorant(
        (1, 0, 0), None, commutators['Bz_integral'], commutators['Bzz_integral'])
    if withheld['numerical_C_M'] is not None:
        raise ArithmeticError('absent L0 A4 integrals produced a C_M')
    manufactured = manufactured_defect_majorant_control()
    if not manufactured['omitting_residual_fails_to_enclose']:
        raise ArithmeticError('zero residual was accepted as an enclosure')
    if manufactured['current_history_C_M'] is not None:
        raise ValueError('manufactured fixture supplied a current-history C_M')
    identities = minor_a3_identities()
    if identities['residuals']['s+1/minor_A3_depends_on_major_A2'] != 'depends':
        raise ArithmeticError('minor A3 no longer retains major A2')
    if identities['residuals']['s-1/minor_A3_depends_on_major_A2'] != 'depends':
        raise ArithmeticError('minor A3 no longer retains major A2')
    report = {
        'schema': SCHEMA,
        'history_identity': profile_identity(family),
        'spinor_record_sha256': digest,
        'phase_gauge': record['phase_gauge'],
        'upstream_order': positive['order'],
        'upstream_cells': positive['cells'],
        'pair': dict(pair),
        'major_A2': _plain(positive['major_A2']),
        'major_A4': _plain(positive['major_A4']),
        'opposite_angular_major_A2': _plain(negative['major_A2']),
        'major_A2_contains_zero': False,
        'opposite_angular_major_A2_contains_zero': False,
        'normalization_residual_contains_zero': True,
        'degree4_pointwise_tail_upper': float(
            positive['degree4_pointwise_tail_upper'].upper()),
        'domain': domain,
        'preimage': preimage,
        'pointwise_counterexample': point,
        'reentry': reentry,
        'carrier': carrier,
        'moments': moments,
        'Bz_integral': commutators['Bz_integral'],
        'Bzz_integral': commutators['Bzz_integral'],
        'commutators_are_C_M': False,
        'unconstructed_term': UNCONSTRUCTED_TERM,
        'numerical_C_M': None,
        'induced_N_beta_remainder': None,
        'current_history_L0_A4_H2_integrals': None,
        'thermal_coherent_tail_included': False,
        'loose_absolute_jet_used_as_L0_A4': False,
        'A2_set_to_zero_from_vanishing_imaginary_rhs': False,
        'full_interval_covered': False,
        'physical_local_gate': 'OPEN',
        'cpu_seconds': time.process_time() - started,
        'rss_mib': _rss_mib(),
    }
    validate_full_interval_uv_report(report)
    return report


def validate_full_interval_uv_report(report):
    if report.get('schema') != SCHEMA:
        raise ValueError('unexpected full-interval schema')
    if not str(report.get('history_identity', '')).startswith('0b0e4ced'):
        raise ValueError('report is not bound to history 0b0e4ced')
    for name in (
            'numerical_C_M', 'induced_N_beta_remainder',
            'current_history_L0_A4_H2_integrals'):
        if report.get(name) is not None:
            raise ValueError(name + ' must stay null; this attempt did not construct it')
        if report.get(name) == 0:
            raise ValueError(name + ' must not be converted to zero')
    if report.get('full_interval_covered') or report['domain']['full_interval_jets_enclosed']:
        raise ValueError('full interval was not enclosed')
    if report['domain']['proper_subinterval_is_full_I']:
        raise ValueError('a proper subinterval is not I')
    if not report['domain']['dyadic_radius_exceeds_axial_inner']:
        raise ValueError('dyadic obstruction was not observed')
    if 'straddles the plateau center' not in report['domain']['domain_failure']:
        raise ValueError('domain failure text changed')
    if not report['preimage']['both_characteristics_exit_support']:
        raise ValueError('characteristic preimage did not leave the support')
    if report['major_A2_contains_zero'] or report['A2_set_to_zero_from_vanishing_imaginary_rhs']:
        raise ValueError('major A2 must remain the nonzero phase-fixed value')
    if not report['reentry']['zero_rhs_drops_positive_reentry_majorant']:
        raise ValueError('zero-RHS substitution must drop a positive majorant')
    if report['reentry']['zero_rhs_reentry_majorant'] != 0:
        raise ValueError('the zero-RHS counterexample value changed')
    if report['pointwise_counterexample']['pointwise_reentry_contains_zero']:
        raise ValueError('sample reentry collapsed to zero')
    if report['pointwise_counterexample']['sample_is_full_interval']:
        raise ValueError('one sample ray is not the interval I')
    if not report['carrier']['F_z_equals_energy_times_F']:
        raise ValueError('carrier restoration identity failed')
    if report['carrier']['this_integral_is_history_minus_reference']:
        raise ValueError('absolute upstream tail is not the history difference')
    if report['thermal_coherent_tail_included'] or report['loose_absolute_jet_used_as_L0_A4']:
        raise ValueError('tail or loose jet was promoted to the remainder')
    if report['commutators_are_C_M'] or report['moments']['moment_fold_supplies_C_M']:
        raise ValueError('commutators or moments were promoted to C_M')
    if report['physical_local_gate'] != 'OPEN':
        raise ValueError('local gate stays open')
    if not report['cpu_seconds'] < 180:
        raise ValueError('full-interval attempt exceeded the 180 second budget')
    if not report['rss_mib'] < 1024:
        raise ValueError('full-interval attempt exceeded the 1 GiB memory target')
    return True
