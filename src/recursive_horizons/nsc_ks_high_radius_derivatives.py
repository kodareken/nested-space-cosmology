"""Spatial supremum and Fourier-l1 bounds for q=1/(r_ref+delta r)-1/r_ref.

The owner reuses the committed p=16 fine Fourier profile table and the
existing radius enclosure. It does not evolve a field or source, fit a
profile, or change the action. Spatial/Fourier geometry only: mixed T/z
jets, Moyal remainder, and physical incoming-gate certificates are other
owners. A covariance UV bound may use auxiliary P in C=P+D without a
band-action change; generating Gamma_sub is a different use.
"""
from fractions import Fraction
from hashlib import sha256
from math import comb, nextafter
from pathlib import Path
import json

import numpy as np
from flint import acb, arb, ctx, fmpq

from .nsc_ks_ball_trajectory import exact_upper, restored_upper
from .nsc_ks_profile_fourier_bound import alias_and_tail_bounds


DEFAULT_BITS = 160
MAX_REQUESTED_DERIVATIVE = 9
PROFILE_KEY = 'profile_1_0'
PERIOD_LENGTH = Fraction(205, 512)
FINE_PROFILE_JSON = 'results/development/nsc-ks-fine-profile-table.json'
FINE_PROFILE_NPZ = 'results/development/artifacts/nsc-ks-fine-profile-table.npz'
RADIUS_BOUND_JSON = 'results/development/nsc-ks-radius-bound.json'

MISSING_NORMAL_TIME = (
    'MISSING: this owner bounds spatial supremum and Fourier-l1 moments only; '
    'mixed T/z jets through order 9 are a separate root owner, and time '
    'integrals of the remainder are not supplied here')
MISSING_FOURTH_ORDER_SYMBOL = (
    'OPEN: high derivatives of the actual fourth-order reference symbol are '
    'not supplied')
MISSING_MOYAL = (
    'OPEN: the exact Moyal remainder is not supplied; this is geometry input '
    'only')
MISSING_BLOCH = (
    'OPEN: Bloch trace inequality is owned by a separate job')
MISSING_P4 = (
    'OPEN: the actual fourth-order reciprocal remainder is reviewed by root; '
    'not this owner')
MISSING_UV_RESIDUAL = (
    'OPEN: a covariance UV bound still needs symbol, Moyal, and time-integral '
    'inputs; P may be purely auxiliary in C=P+D without changing the action. '
    'Generating Gamma_sub is a different use')
PHYSICAL_LOCAL_GATE = (
    'OPEN: spatial geometry coefficients only; not a UV residual or physical '
    'incoming-gate certificate')


def repository_root():
    return Path(__file__).resolve().parents[2]


def _collapsed_upper_record(value):
    return (isinstance(value, dict) and 'mantissa' in value and 'exponent' in value
            and 'exact_rational' not in value)


def _enclosing_arb(value, name):
    """Convert to a flint enclosure without collapsing to a directed upper."""
    if value is None or isinstance(value, bool):
        raise ValueError('explicit nonnegative bound required: ' + name)
    if isinstance(value, dict):
        if 'exact_rational' in value:
            q = Fraction(value['exact_rational'])
            return arb(fmpq(q.numerator, q.denominator))
        raise ValueError('exact rational or live enclosure required: ' + name)
    if isinstance(value, arb):
        return value
    if isinstance(value, acb):
        raise ValueError('real quantity required: ' + name)
    try:
        q = Fraction(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError('finite bound required: ' + name) from error
    return arb(fmpq(q.numerator, q.denominator))


def _norm_upper(value, name):
    """Directed upper for nonnegative norms and coefficient magnitudes."""
    if value is None or isinstance(value, bool):
        raise ValueError('explicit nonnegative bound required: ' + name)
    if isinstance(value, dict) and 'mantissa' in value and 'exponent' in value:
        result = restored_upper(value)
    elif isinstance(value, acb):
        result = value.abs_upper()
    elif isinstance(value, arb):
        result = value
    else:
        result = _enclosing_arb(value, name)
    if not result.is_finite() or not result >= 0:
        raise ValueError('finite nonnegative bound required: ' + name)
    return result.upper()


def _positive_enclosure(value, name):
    """Positive enclosing interval for a denominator. Reject upper-only records."""
    if _collapsed_upper_record(value):
        raise ValueError('collapsed upper-only record cannot be a positive denominator: ' + name)
    result = _enclosing_arb(value, name)
    if not result.is_finite() or not result > 0:
        raise ValueError('positive enclosing interval required: ' + name)
    return result


def _positive_lower(value, name):
    """Directed lower of a positive geometry bound used as a denominator."""
    return _positive_enclosure(value, name).lower()


def positive_reciprocal_upper(value, *, bits=DEFAULT_BITS):
    """Directed upper of 1/x from a positive enclosing interval, never x.upper()."""
    with ctx.workprec(bits):
        return (1 / _positive_enclosure(value, 'positive denominator')).upper()


def _positive_int(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError('nonnegative integer required: ' + name)
    return value


def pack_upper(value):
    """Directed upper as an exact rational or dyadic; binary64 is display only."""
    if isinstance(value, dict) and 'mantissa' in value:
        value = restored_upper(value)
    if isinstance(value, arb):
        upper = value.upper()
        number = float(upper)
        if arb(number) < upper:
            number = nextafter(number, float('inf'))
        return {**exact_upper(upper), 'binary64_upper': number, 'kind': 'dyadic'}
    q = Fraction(value)
    if q < 0:
        raise ValueError('nonnegative packed upper required')
    number = float(q)
    if Fraction(number) < q:
        number = nextafter(number, float('inf'))
    return {'exact_rational': str(q), 'binary64_upper': number, 'kind': 'rational'}


def pack_lower(value):
    """Directed lower as an exact rational or dyadic; not an upper record."""
    if isinstance(value, dict) and 'mantissa' in value:
        value = restored_upper(value)
    if isinstance(value, arb):
        lower = value.lower()
        number = float(lower)
        if arb(number) > lower:
            number = nextafter(number, float('-inf'))
        mantissa, exponent = lower.man_exp()
        return {'mantissa': str(mantissa), 'exponent': int(exponent),
                'binary64_lower': number, 'kind': 'dyadic_lower'}
    q = Fraction(value)
    if q < 0:
        raise ValueError('nonnegative packed lower required')
    number = float(q)
    if Fraction(number) > q:
        number = nextafter(number, float('-inf'))
    return {'exact_rational': str(q), 'binary64_lower': number, 'kind': 'rational'}


def restore_packed(item):
    """Restore a pack_upper or pack_lower record at the current flint precision."""
    if not isinstance(item, dict):
        raise TypeError('packed directed bound required')
    if item.get('kind') == 'rational' or ('exact_rational' in item and 'mantissa' not in item):
        q = Fraction(item['exact_rational'])
        return arb(fmpq(q.numerator, q.denominator))
    return restored_upper(item)


def unpack_serialized_ball(row):
    """Restore one complex ball from the table's dyadic midpoint/radius parts."""
    parts = []
    for component in row:
        if len(component) != 4:
            raise ValueError('real/imag dyadic midpoint and radius required')
        mid_m, mid_e, rad_m, rad_e = (int(v) for v in component)
        parts.append(arb((mid_m, mid_e), (rad_m, rad_e)))
    if len(parts) != 2:
        raise ValueError('complex Fourier balls require two real parts')
    return acb(*parts)


def sha256_file(path):
    digest = sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def authenticate_profile_payload(record, npz_path):
    """Compare the committed NPZ identity with the JSON payload descriptor."""
    payload = record['payload']
    path = Path(npz_path)
    size = path.stat().st_size
    if size != payload['bytes']:
        raise ValueError('fine profile table size changed')
    digest = sha256_file(path)
    if digest != payload['sha256']:
        raise ValueError('fine profile table hash changed')
    return digest


def load_fine_profile_coefficients(npz_path, *, keys=(PROFILE_KEY,), bits=DEFAULT_BITS):
    """Load metadata and only the requested serialized coefficient arrays."""
    path = Path(npz_path)
    with np.load(path, allow_pickle=False) as arrays:
        available = set(arrays.files)
        if 'metadata_json' not in available:
            raise ValueError('fine profile table metadata required')
        meta = json.loads(arrays['metadata_json'].tobytes())
        missing = [key for key in keys if key + '/coefficients' not in available]
        if missing:
            raise ValueError('missing profile coefficient arrays: ' + ','.join(missing))
        with ctx.workprec(bits):
            coefficients = {
                key: tuple(unpack_serialized_ball(row) for row in arrays[key + '/coefficients'])
                for key in keys
            }
    return meta, coefficients


def fourier_sup_derivative_uppers(coefficients, period_length, derivative_l1,
                                  derivative_order, quadrature_points, retained_index,
                                  max_derivative, *, bits=DEFAULT_BITS):
    """A_j = sum |omega_k|^j |c_k| + analytic omitted tail, j < p-1.

    A_j is the Fourier-l1 moment and also a spatial-sup majorant of f^(j).
    The period is kept as a positive enclosing interval; L.upper() is not
    used as a denominator. Displayed float tails are not used.
    """
    p = _positive_int(derivative_order, 'derivative order')
    M = _positive_int(quadrature_points, 'quadrature points')
    K = _positive_int(retained_index, 'retained index')
    j_max = _positive_int(max_derivative, 'max derivative')
    if p <= j_max + 1:
        raise ValueError('integrable omitted tails require j < p-1')
    coefficients = tuple(coefficients)
    if len(coefficients) != 2 * K + 1:
        raise ValueError('retained Fourier inventory must have length 2K+1')
    with ctx.workprec(bits):
        length = _positive_enclosure(period_length, 'period length')
        l1 = _norm_upper(derivative_l1, 'derivative L1')
        alias, tails = alias_and_tail_bounds(l1, length, p, M, K, j_max)
        omega = 2 * arb.pi() / length
        band = []
        for order in range(j_max + 1):
            total = arb(0)
            for index, coefficient in enumerate(coefficients):
                mode = index - K
                magnitude = coefficient.abs_upper() if isinstance(coefficient, acb) else _norm_upper(
                    coefficient, 'coefficient')
                if order == 0:
                    total += magnitude
                elif mode:
                    total += (omega * abs(mode)) ** order * magnitude
            band.append(total.upper())
        moments = tuple((term + tail).upper() for term, tail in zip(band, tails))
        return {
            'alias_error_upper': alias,
            'omitted_tail_uppers': tails,
            'band_sup_uppers': tuple(band),
            'fourier_l1_moments': moments,
            'sup_uppers': moments,
            'period_length': length,
            'derivative_l1_upper': l1,
            'derivative_order': p,
            'quadrature_points': M,
            'retained_index': K,
            'max_derivative': j_max,
            'precision_bits': bits,
        }


def perturbed_radius_derivative_uppers(w_sup, u_sup, normal_support, *, bits=DEFAULT_BITS):
    """|s|<=sigma and 0<=chi<=1 give |delta r^(j)| <= sigma W_j + sigma^3 U_j/6."""
    if len(w_sup) != len(u_sup):
        raise ValueError('matching w and U derivative orders required')
    with ctx.workprec(bits):
        sigma = _norm_upper(normal_support, 'normal support')
        return tuple(
            (sigma * _norm_upper(w, 'w sup') + sigma ** 3 * _norm_upper(u, 'U sup') / 6).upper()
            for w, u in zip(w_sup, u_sup))


def _leibniz_reciprocal_moments(delta_moments, inverse_zeroth, denominator):
    q = [inverse_zeroth]
    for n in range(1, len(delta_moments)):
        acc = arb(0)
        for j in range(1, n + 1):
            acc += arb(comb(n, j)) * delta_moments[j] * q[n - j]
        q.append((acc / denominator).upper())
    return tuple(q)


def reciprocal_full_derivative_uppers(delta_r_uppers, radius_lower, *, bits=DEFAULT_BITS):
    """Pointwise q_full=1/r from r q_full=1. Zeroth upper is 1/r_min."""
    if not delta_r_uppers:
        raise ValueError('at least the zeroth radius perturbation bound is required')
    with ctx.workprec(bits):
        rmin = _positive_enclosure(radius_lower, 'radius lower')
        radius = tuple(_norm_upper(value, 'delta r derivative') for value in delta_r_uppers)
        return _leibniz_reciprocal_moments(radius, (1 / rmin).upper(), rmin)


def reciprocal_difference_derivative_uppers(delta_r_uppers, radius_lower,
                                            reference_radius_lower, *, bits=DEFAULT_BITS):
    """Pointwise q=1/(r_ref+delta r)-1/r_ref. Zeroth uses delta r_max/(r_ref_min r_min)."""
    with ctx.workprec(bits):
        rmin = _positive_enclosure(radius_lower, 'radius lower')
        rref = _positive_enclosure(reference_radius_lower, 'reference radius')
        delta = tuple(_norm_upper(value, 'delta r derivative') for value in delta_r_uppers)
        if not rref.lower() > delta[0]:
            raise ValueError('strict radius margin r_ref > delta r_max required')
        q_full = _leibniz_reciprocal_moments(delta, (1 / rmin).upper(), rmin)
        q0 = (delta[0] / (rref * rmin)).upper()
        return (q0,) + tuple(q_full[1:])


def reciprocal_fourier_moment_uppers(delta_r_moments, reference_radius_lower, *, bits=DEFAULT_BITS):
    """Wiener-algebra moments of 1/r and of 1/r-1/r_ref.

    A_0(delta r)=d0 must be strictly below r_ref,min. Neumann gives
    A_0(1/r)<=1/(r_ref,min-d0). Higher A_n follow the product rule with that
    inverse A_0 bound, not a pointwise r_min. Difference A_0<=d0/[r_ref,min(r_ref,min-d0)].
    """
    if not delta_r_moments:
        raise ValueError('at least the zeroth Fourier moment of delta r is required')
    with ctx.workprec(bits):
        delta = tuple(_norm_upper(value, 'delta r Fourier moment') for value in delta_r_moments)
        rref_lo = _positive_lower(reference_radius_lower, 'reference radius')
        d0 = delta[0]
        if not rref_lo > d0:
            raise ValueError('Neumann radius margin r_ref,min > ||delta r||_A0 required')
        margin = (rref_lo - d0).lower()
        if not margin > 0:
            raise ValueError('Neumann radius margin r_ref,min > ||delta r||_A0 required')
        q_full = _leibniz_reciprocal_moments(delta, (1 / margin).upper(), margin)
        q0 = (d0 / (rref_lo * margin)).upper()
        return {
            'q_full': q_full,
            'q_difference': (q0,) + tuple(q_full[1:]),
            'neumann_margin_lower': margin,
            'delta_r_A0': d0,
        }


def compact_support_norms(sup_uppers, width, *, bits=DEFAULT_BITS):
    """Functions vanishing off an interval of length w satisfy L1<=w sup, L2^2<=w sup^2."""
    with ctx.workprec(bits):
        support = _norm_upper(width, 'axial support width')
        if not support > 0:
            raise ValueError('positive axial support width required')
        rows = []
        for value in sup_uppers:
            sup = _norm_upper(value, 'sup')
            rows.append({
                'sup': sup,
                'L1': (support * sup).upper(),
                'L2_squared': (support * sup * sup).upper(),
            })
        return tuple(rows)


def _rational_bound(item):
    return Fraction(item['exact_rational'])


def _hex_float(value):
    return Fraction(float.fromhex(value))


def bound_pure_radius_difference(*, root=None, max_derivative=MAX_REQUESTED_DERIVATIVE,
                                 bits=DEFAULT_BITS, verify_hash=True):
    """Directed spatial and Fourier-l1 bounds on the owned nonzero pure-radius control."""
    j_max = _positive_int(max_derivative, 'max derivative')
    if j_max < MAX_REQUESTED_DERIVATIVE:
        raise ValueError('this owner supplies at least orders 0 through 9')
    root = Path(root) if root is not None else repository_root()
    table = json.loads((root / FINE_PROFILE_JSON).read_text())
    radius = json.loads((root / RADIUS_BOUND_JSON).read_text())
    npz_path = root / table['payload']['path']
    if table['payload']['path'] != FINE_PROFILE_NPZ:
        raise ValueError('fine profile payload path changed')
    digest = authenticate_profile_payload(table, npz_path) if verify_hash else table['payload']['sha256']
    meta, coefficients = load_fine_profile_coefficients(npz_path, keys=(PROFILE_KEY,), bits=bits)
    if meta['axial_profile_identity'] != table['axial_profile_identity']:
        raise ValueError('physical profile fingerprint disagrees with the JSON record')
    proof = meta['profiles'][PROFILE_KEY]
    if proof['powers'] != [1, 0]:
        raise ValueError('the pure-radius control must use physical W=profile_1_0')
    settings = table['settings']
    if settings['derivative_order'] != 16 or settings['quadrature_points'] != 65536 or settings['retained_index'] != 8192:
        raise ValueError('committed p=16 fine Fourier settings required')
    period = Fraction(table['period_length_exact'])
    if period != PERIOD_LENGTH or Fraction(meta['period_length']) != period:
        raise ValueError('numerical period changed')
    directions = meta['profile_description']['directions']
    if len(directions) != 1:
        raise ValueError('owned one-direction radius family required')
    direction = directions[0]
    if any(_hex_float(c) != 0 for c in direction['U']['coefficients']):
        raise ValueError('this control has U=0; a generic U term must be supplied explicitly')
    u_zero = all(_rational_bound(item) == 0 for item in radius['bounds']['U_sup_orders_0_1_2'])
    if not u_zero:
        raise ValueError('radius record U bounds are not zero')
    enclosure_w = tuple(_rational_bound(item) for item in radius['bounds']['w_sup_orders_0_1_2'])
    rref = _rational_bound(radius['bounds']['reference_radius_lower'])
    rmin = _rational_bound(radius['bounds']['radius_lower'])
    eta = _rational_bound(radius['bounds']['radius_ratio_upper'])
    if rmin != rref * (1 - eta):
        raise ValueError('radius lower is not the recorded two-derivative margin')
    sigma = eta * rref / enclosure_w[0]
    width = 2 * _hex_float(direction['w']['outer'])
    amplitude = _hex_float(direction['amplitude'])
    if amplitude == 0:
        raise ValueError('nonzero control amplitude required')
    with ctx.workprec(bits):
        fourier = fourier_sup_derivative_uppers(
            coefficients[PROFILE_KEY], period, proof['derivative_l1_upper'],
            settings['derivative_order'], settings['quadrature_points'],
            settings['retained_index'], j_max, bits=bits)
        w_moments = fourier['fourier_l1_moments']
        u_sup = (arb(0),) * (j_max + 1)
        delta_r = perturbed_radius_derivative_uppers(w_moments, u_sup, sigma, bits=bits)
        q_full = reciprocal_full_derivative_uppers(delta_r, rmin, bits=bits)
        q_diff = reciprocal_difference_derivative_uppers(delta_r, rmin, rref, bits=bits)
        norms = compact_support_norms(q_diff, width, bits=bits)
        fourier_q = reciprocal_fourier_moment_uppers(delta_r, rref, bits=bits)
        enclosure_arb = tuple(arb(fmpq(item.numerator, item.denominator)) for item in enclosure_w)
        if any(w_moments[j] > enclosure_arb[j] for j in range(min(3, j_max + 1))):
            raise ArithmeticError('Fourier W_j must not exceed the two-derivative enclosure for j<=2')
        rmin_encl = arb(fmpq(rmin.numerator, rmin.denominator))
        neumann = fourier_q['neumann_margin_lower']
        pointwise_conservative = rmin_encl <= neumann
    return {
        'schema': 'NSC-KS-HIGH-RADIUS-DERIVATIVES-v2',
        'profile_key': PROFILE_KEY,
        'axial_profile_identity': table['axial_profile_identity'],
        'profile_payload_sha256': digest,
        'period_length_exact': str(period),
        'normal_support': pack_upper(sigma),
        'axial_support_width': pack_upper(width),
        'reference_radius_lower': pack_lower(rref),
        'radius_lower': pack_lower(rmin),
        'neumann_margin_lower': pack_lower(fourier_q['neumann_margin_lower']),
        'delta_r_A0': pack_upper(fourier_q['delta_r_A0']),
        'pointwise_radius_lower_is_conservative_for_neumann': bool(pointwise_conservative),
        'history_amplitude_included_in_W': pack_upper(amplitude),
        'U_identically_zero': True,
        'alpha_not_applied_twice': True,
        'measure_not_applied_twice': True,
        'derivative_orders': list(range(j_max + 1)),
        'W_fourier_l1': [pack_upper(value) for value in w_moments],
        'W_sup': [pack_upper(value) for value in w_moments],
        'delta_r_sup': [pack_upper(value) for value in delta_r],
        'q_full_sup': [pack_upper(value) for value in q_full],
        'q_difference_sup': [pack_upper(value) for value in q_diff],
        'q_difference_L1': [pack_upper(row['L1']) for row in norms],
        'q_difference_L2_squared': [pack_upper(row['L2_squared']) for row in norms],
        'q_full_fourier_l1': [pack_upper(value) for value in fourier_q['q_full']],
        'q_difference_fourier_l1': [pack_upper(value) for value in fourier_q['q_difference']],
        'omitted_profile_tails': [pack_upper(value) for value in fourier['omitted_tail_uppers']],
        'fourier_band_sup': [pack_upper(value) for value in fourier['band_sup_uppers']],
        'settings': {
            'bits': bits,
            'table_bits': settings['bits'],
            'derivative_order': settings['derivative_order'],
            'quadrature_points': settings['quadrature_points'],
            'retained_index': settings['retained_index'],
            'max_derivative': j_max,
        },
        'reused': {
            'fine_profile_table': table['payload'],
            'radius_bound_schema': radius['schema'],
            'alias_and_tail_bounds': 'nsc_ks_profile_fourier_bound.alias_and_tail_bounds',
            'two_derivative_radius_lower': True,
        },
        'missing': {
            'normal_time_derivatives': MISSING_NORMAL_TIME,
            'fourth_order_symbol_high_derivatives': MISSING_FOURTH_ORDER_SYMBOL,
            'exact_Moyal_remainder': MISSING_MOYAL,
            'Bloch_trace_inequality': MISSING_BLOCH,
            'P4_remainder': MISSING_P4,
            'UV_residual': MISSING_UV_RESIDUAL,
            'physical_local_gate': PHYSICAL_LOCAL_GATE,
        },
        'scope': {
            'finite_history_geometry_bound': True,
            'source_action_Gamma_rest_changed': False,
            'field_runs': 0,
            'source_runs': 0,
            'new_profile_fit': False,
            'physical_local_gate': 'OPEN',
            'UV_certificate': False,
            'spatial_supremum_and_Fourier_l1_only': True,
            'band_action_required_for_covariance_UV': False,
        },
    }
