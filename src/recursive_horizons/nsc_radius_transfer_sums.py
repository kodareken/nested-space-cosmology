"""Sharp directed weighted Fourier sums of q=1/r_g-1/r_ref from saved W^p.

Root assembles integrated Weyl remainders and global-projector bounds. The
omega-only kernels here contain no dk/(2pi), no sum over sign(k), and no
spin or source multiplicity; callers apply those factors once. Both +/-omega
modes already sit in the stored coefficient inventory.

Canonical K is a numerical momentum split, not a source-energy cutoff.
This owner does not evolve a field or source, change the profile or action,
or close the physical incoming gate.
"""
from fractions import Fraction
from contextlib import contextmanager
from pathlib import Path
import json

from flint import acb, arb, arb_series, ctx

from .nsc_ks_high_radius_derivatives import (
    DEFAULT_BITS, FINE_PROFILE_JSON, FINE_PROFILE_NPZ, PERIOD_LENGTH,
    PROFILE_KEY, RADIUS_BOUND_JSON, authenticate_profile_payload,
    fourier_sup_derivative_uppers, load_fine_profile_coefficients,
    pack_lower, pack_upper, repository_root,
    _enclosing_arb, _hex_float, _norm_upper, _positive_enclosure,
    _positive_int, _rational_bound,
)
from .nsc_ks_profile_fourier_bound import alias_and_tail_bounds


SCHEMA = 'NSC-RADIUS-TRANSFER-SUMS-v1'
POLYNOMIAL_ORDER = 4
PROFILE_POWERS = (1, 2, 3, 4)
PROFILE_KEYS = tuple('profile_%d_0' % power for power in PROFILE_POWERS)
AUTHENTIC_PAYLOAD_SHA256 = (
    'fc20111d0c72f1ac81d56ee5ea64cdd1bb94a9bfc6e86edfa19d1f8992f0ce33')
AXIAL_PROFILE_IDENTITY = (
    'a4c061f297027bcb9934fa2338bb81fb799ca4ec26a025f6bc563384c70cd69d')
MAX_OWNED_ENVELOPE_DEGREE = 4

MISSING_TIME_INTEGRAL = (
    'OPEN: time integrals of the weighted remainder are not supplied')
MISSING_PROJECTOR_SYMBOL = (
    'OPEN: fourth-order symbol/global-projector coefficients are assembled '
    'by root, not here')
MISSING_MOYAL = (
    'OPEN: the exact Moyal remainder is not supplied')
MISSING_BLOCH = (
    'OPEN: Bloch trace inequality is owned by a separate job')
MISSING_UV = (
    'OPEN: this is a geometry-weighted Fourier sum, not a UV residual or '
    'physical incoming-gate certificate')
PHYSICAL_LOCAL_GATE = (
    'OPEN: radius-control Fourier sums only; not a UV or physical-gate '
    'certificate')

KERNEL_MEASURE = {
    'dk_over_2pi': False,
    'both_sign_k': False,
    'spin_source_multiplicity': False,
    'plus_minus_omega_in_inventory': True,
    'extra_factor_two_for_negative_omega': False,
    'canonical_K_is_source_energy_cutoff': False,
}


def _as_bit_precision(bits):
    bits = _positive_int(bits, 'precision bits')
    if bits < 2:
        raise ValueError('positive precision bits required')
    return bits


def _as_s(value):
    if isinstance(value, bool) or value not in (0, 1):
        raise ValueError('s must be 0 or 1')
    return int(value)


def _as_small_n(value):
    if isinstance(value, bool) or value not in (1, 2, 3, 4, 5):
        raise ValueError('small-transfer n must be 1..5')
    return int(value)


def _as_taylor_n(value):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError('Taylor n must be a nonnegative integer')
    return int(value)


def _as_decay_p(value):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError('unshifted decay p>=1 required; p=0/n=0 uses f_shift')
    return int(value)


def _abs_frequency(value, name):
    if value is None or isinstance(value, bool):
        raise ValueError('explicit nonnegative frequency required: ' + name)
    result = abs(_enclosing_arb(value, name))
    if not result.is_finite() or not result >= 0:
        raise ValueError('finite nonnegative frequency required: ' + name)
    return result


def _canonical_cutoff(value, *, bits=DEFAULT_BITS):
    with ctx.workprec(bits):
        return _positive_enclosure(value, 'canonical momentum split')


def _positive_power_upper(value, exponent, name):
    if isinstance(exponent, bool) or int(exponent) != exponent:
        raise ValueError('integer power required: ' + name)
    exponent = int(exponent)
    base = _positive_enclosure(value, name)
    if exponent == 0:
        return arb(1)
    return (base ** exponent).upper()


def _directed_inv_upper(value, name):
    return (1 / _positive_enclosure(value, name)).upper()


def _monomial_envelope(degree, coefficient):
    degree = _positive_int(degree, 'envelope degree')
    coeff = _norm_upper(coefficient, 'envelope coefficient')
    terms = [arb(0)] * (degree + 1)
    terms[degree] = coeff
    return tuple(terms)


def evaluate_polynomial_envelope(coefficients, abs_omega, *, bits=DEFAULT_BITS):
    """Directed upper of sum b_j |omega|^j for nonnegative b_j."""
    if not coefficients:
        raise ValueError('nonempty polynomial envelope required')
    with ctx.workprec(bits):
        w = _abs_frequency(abs_omega, 'frequency')
        total = arb(0)
        power = arb(1)
        for coefficient in coefficients:
            total += _norm_upper(coefficient, 'envelope coefficient') * power
            power *= w
        return total.upper()


def contract_polynomial_envelope(moments, envelope_coefficients, *, bits=DEFAULT_BITS):
    """sum b_j A_j. Missing moments are never treated as zero."""
    if not envelope_coefficients:
        raise ValueError('nonempty polynomial envelope required')
    if moments is None:
        raise ValueError('explicit remainder or tail moments required')
    with ctx.workprec(bits):
        if len(moments) < len(envelope_coefficients):
            raise ValueError(
                'missing Fourier moment for polynomial envelope order; '
                'unknown tails are not zero')
        total = arb(0)
        for moment, coefficient in zip(moments, envelope_coefficients):
            total += (_norm_upper(moment, 'Fourier moment')
                      * _norm_upper(coefficient, 'envelope coefficient'))
        return total.upper()


def f_small(n, s, abs_omega, canonical_cutoff, *, bits=DEFAULT_BITS):
    """(w/2)^n / ((5-s) max(K,w)^{5-s}), the integral of |k|^{s-6}(w/2)^n.

    Integration is over k>=max(K,w) only. Zero at w=0. n=1..5, s=0 or 1.
    """
    n = _as_small_n(n)
    s = _as_s(s)
    with ctx.workprec(bits):
        w = _abs_frequency(abs_omega, 'small-transfer frequency')
        cutoff = _canonical_cutoff(canonical_cutoff, bits=bits)
        peak = cutoff.max(w)
        return ((w / 2) ** n / (arb(5 - s) * peak ** (5 - s))).upper()


def f_small_envelope(n, s, canonical_cutoff, *, bits=DEFAULT_BITS):
    """Polynomial envelope of f_small. Constant if n<=5-s; linear in w if n=5,s=1."""
    n = _as_small_n(n)
    s = _as_s(s)
    with ctx.workprec(bits):
        cutoff = _canonical_cutoff(canonical_cutoff, bits=bits)
        scale = arb(2) ** n * arb(5 - s)
        if n <= 5 - s:
            power = n - (5 - s)
            coeff = _positive_power_upper(cutoff, power, 'canonical momentum split') / scale
            return (coeff.upper(),)
        coeff = arb(1) / scale
        return (arb(0), coeff.upper())


def f_shift(s, abs_omega, canonical_cutoff, *, bits=DEFAULT_BITS):
    """1_{w>K}(w^{s+1}-K^{s+1})/(s+1). P_inf / global P0 uses this at n=0,p=0."""
    s = _as_s(s)
    with ctx.workprec(bits):
        w = _abs_frequency(abs_omega, 'shifted-symbol frequency')
        cutoff = _canonical_cutoff(canonical_cutoff, bits=bits)
        w_hi = w.upper()
        k_lo = cutoff.lower()
        if not (w_hi > k_lo):
            return arb(0)
        return ((w ** (s + 1) - cutoff ** (s + 1)) / arb(s + 1)).upper()


def f_shift_envelope(s, canonical_cutoff, *, bits=DEFAULT_BITS):
    """Envelope w^{s+1}/(s+1). The cutoff is required only as a positive split."""
    s = _as_s(s)
    with ctx.workprec(bits):
        _canonical_cutoff(canonical_cutoff, bits=bits)
        return _monomial_envelope(s + 1, arb(1) / arb(s + 1))


def _power_integral(cutoff, width, exponent):
    """int_K^w t^exponent dt for a large-transfer upper, exponent integer."""
    if exponent == -1:
        ratio = width / cutoff
        if not ratio > 0:
            return arb(0)
        return ratio.log()
    return (width ** (exponent + 1) - cutoff ** (exponent + 1)) / arb(exponent + 1)


def f_taylor(n, p, s, abs_omega, canonical_cutoff, *, bits=DEFAULT_BITS):
    """1_{w>K}(w/2)^n int_K^w k^{s-p} dk. p>=1 is the unshifted decay power."""
    n = _as_taylor_n(n)
    p = _as_decay_p(p)
    s = _as_s(s)
    with ctx.workprec(bits):
        w = _abs_frequency(abs_omega, 'Taylor frequency')
        cutoff = _canonical_cutoff(canonical_cutoff, bits=bits)
        w_hi = w.upper()
        k_lo = cutoff.lower()
        if not (w_hi > k_lo):
            return arb(0)
        return ((w / 2) ** n * _power_integral(cutoff, w, s - p)).upper()


def taylor_envelope_degree(n, p, s):
    n = _as_taylor_n(n)
    p = _as_decay_p(p)
    s = _as_s(s)
    if p > s + 1:
        return n
    if p == s + 1:
        return n + 1
    return n + s - p + 1


def f_taylor_envelope(n, p, s, canonical_cutoff, *, bits=DEFAULT_BITS):
    """Polynomial envelope of f_taylor. Owned calls with p=j+1+n have degree<=4."""
    n = _as_taylor_n(n)
    p = _as_decay_p(p)
    s = _as_s(s)
    with ctx.workprec(bits):
        cutoff = _canonical_cutoff(canonical_cutoff, bits=bits)
        two_n = arb(2) ** n
        if p > s + 1:
            coeff = (_positive_power_upper(cutoff, s - p + 1, 'canonical momentum split')
                     / (two_n * arb(p - s - 1)))
            return _monomial_envelope(n, coeff)
        if p == s + 1:
            coeff = arb(1) / (two_n * cutoff)
            return _monomial_envelope(n + 1, coeff.upper())
        coeff = arb(1) / (two_n * arb(s - p + 1))
        return _monomial_envelope(n + s - p + 1, coeff)


def owned_taylor_parameters():
    """(n,p,s) with p=j+1+n for projector order j=0..4 and Taylor n=0..4."""
    for projector in range(5):
        for n in range(5):
            for s in (0, 1):
                yield n, projector + 1 + n, s


def validate_owned_taylor_envelope_orders():
    """Owned p=j+1+n envelopes stay at polynomial degree <=4."""
    degrees = tuple(taylor_envelope_degree(n, p, s) for n, p, s in owned_taylor_parameters())
    if max(degrees) > MAX_OWNED_ENVELOPE_DEGREE:
        raise ArithmeticError('owned Taylor envelope degree exceeds 4')
    return {
        'max_degree': max(degrees),
        'min_degree': min(degrees),
        'count': len(degrees),
        'limit': MAX_OWNED_ENVELOPE_DEGREE,
    }


@contextmanager
def _series_capacity(order):
    previous=ctx.cap
    ctx.cap=order+1
    try:
        yield
    finally:
        ctx.cap=previous


def reciprocal_remainder_moments(delta_r_moments, reference_radius_lower, *,
                                 polynomial_order=POLYNOMIAL_ORDER,
                                 max_moment=None, bits=DEFAULT_BITS):
    """Wiener A_j(q-q_P) from the positive majorant x^{P+1}/(1-x).

    eta_j=A_j(delta r)/r_ref,min. arb_series uses eta_j/j! through at least
    order P. A_j(remainder)=j! [t^j] f(eta(t))/r_ref,min. This retains every
    geometric order >=P+1; it is not a difference of two upper bounds.
    eta_0<1 is required and 1-eta_0 uses a directed lower denominator.
    Missing moments are rejected rather than filled with zero.
    """
    order = _positive_int(polynomial_order, 'polynomial order')
    if order < 1:
        raise ValueError('positive polynomial truncation order required')
    if max_moment is None:
        max_moment = order
    max_moment = _positive_int(max_moment, 'remainder moment')
    if delta_r_moments is None:
        raise ValueError('explicit delta r Fourier moments required')
    if len(delta_r_moments) < max_moment + 1:
        raise ValueError(
            'delta r moments through the remainder order required; '
            'missing moments are not zero')
    with ctx.workprec(bits), _series_capacity(max_moment):
        rref = _positive_enclosure(reference_radius_lower, 'reference radius')
        inv = (1 / rref).upper()
        eta = tuple(
            (_norm_upper(value, 'delta r Fourier moment') * inv).upper()
            for value in delta_r_moments[:max_moment + 1])
        eta0 = eta[0]
        if not eta0 < 1:
            raise ValueError('eta0<1 required for the Wiener remainder majorant')
        length = max_moment + 1
        coeffs = [eta[j] / arb(j).fac() for j in range(length)]
        eta_series = arb_series(coeffs, length)
        one = arb_series([arb(1)], length)
        majorant = (eta_series ** (order + 1)) / (one - eta_series)
        moments = tuple(
            (arb(j).fac() * majorant[j] * inv).upper() for j in range(length))
        return {
            'moments': moments,
            'eta0_upper': eta0,
            'eta': eta,
            'polynomial_order': order,
            'majorant': 'x^{P+1}/(1-x)',
            'not_difference_of_uppers': True,
            'pointwise_rmin_used_as_Wiener_inverse': False,
        }


def q4_mode_norm_envelope(profile_values, sigma, reference_radius, *,
                          bits=DEFAULT_BITS, polynomial_order=POLYNOMIAL_ORDER):
    """sum_{p=1}^P sigma^p |W^p_hat| / r_ref,min^{p+1}. No extra alpha."""
    order = _positive_int(polynomial_order, 'polynomial order')
    if profile_values is None:
        raise ValueError('explicit W^p mode values required')
    if len(profile_values) < order:
        raise ValueError('W^1..W^P mode values required; missing powers are not zero')
    with ctx.workprec(bits):
        scale = _norm_upper(sigma, 'normal support')
        inv = _directed_inv_upper(reference_radius, 'reference radius')
        total = arb(0)
        sigma_power = scale
        inv_power = inv * inv
        for power in range(1, order + 1):
            value = profile_values[power - 1]
            if value is None:
                raise ValueError('missing W^p Fourier coefficient cannot be zero')
            magnitude = value.abs_upper() if isinstance(value, acb) else _norm_upper(
                value, 'W^p coefficient')
            total += sigma_power * inv_power * magnitude
            sigma_power *= scale
            inv_power *= inv
        return total.upper()


def _omega_unit(period_length, *, bits=DEFAULT_BITS):
    with ctx.workprec(bits):
        length = _positive_enclosure(period_length, 'period length')
        return (2 * arb.pi() / length)


def weighted_retained_q4_sum(coefficients_by_power, period_length, sigma,
                             reference_radius, weight, *,
                             bits=DEFAULT_BITS, polynomial_order=POLYNOMIAL_ORDER,
                             skip_zero_mode=True):
    """Sum actual retained |q4_l| weight(|omega_l|). Mode 0 is homogeneous-ignored."""
    order = _positive_int(polynomial_order, 'polynomial order')
    if coefficients_by_power is None or len(coefficients_by_power) < order:
        raise ValueError('W^1..W^P retained coefficient tables required')
    if weight is None:
        raise ValueError('explicit nonnegative weight required')
    tables = tuple(tuple(coefficients_by_power[power - 1]) for power in range(1, order + 1))
    counts = {len(table) for table in tables}
    if len(counts) != 1:
        raise ValueError('matching retained Fourier inventories required')
    count = counts.pop()
    if count % 2 == 0:
        raise ValueError('retained Fourier inventory must have odd length 2K+1')
    retained_index = count // 2
    with ctx.workprec(bits):
        omega = _omega_unit(period_length, bits=bits)
        total = arb(0)
        for index in range(count):
            mode = index - retained_index
            if skip_zero_mode and mode == 0:
                continue
            values = tuple(table[index] for table in tables)
            envelope = q4_mode_norm_envelope(
                values, sigma, reference_radius, bits=bits,
                polynomial_order=order)
            freq = omega * abs(mode)
            total += envelope * _norm_upper(weight(freq), 'mode weight')
        return total.upper()


def omitted_q4_polynomial_tail(omitted_tails_by_power, sigma, reference_radius,
                               envelope_coefficients, *,
                               bits=DEFAULT_BITS, polynomial_order=POLYNOMIAL_ORDER):
    """Unknown W^p Fourier bands: each power uses its own certified L1 tail."""
    order = _positive_int(polynomial_order, 'polynomial order')
    if omitted_tails_by_power is None or len(omitted_tails_by_power) < order:
        raise ValueError('omitted W^p tails required; missing tails are not zero')
    with ctx.workprec(bits):
        scale = _norm_upper(sigma, 'normal support')
        inv = _directed_inv_upper(reference_radius, 'reference radius')
        sigma_power = scale
        inv_power = inv * inv
        total = arb(0)
        for power in range(1, order + 1):
            tails = omitted_tails_by_power[power - 1]
            if tails is None:
                raise ValueError('missing W^p omitted tail cannot be zero')
            total += sigma_power * inv_power * contract_polynomial_envelope(
                tails, envelope_coefficients, bits=bits)
            sigma_power *= scale
            inv_power *= inv
        return total.upper()


def _pack_split(retained, omitted, remainder):
    retained_u = _norm_upper(retained, 'retained polynomial')
    omitted_u = _norm_upper(omitted, 'omitted polynomial tail')
    remainder_u = _norm_upper(remainder, 'reciprocal remainder')
    total = (retained_u + omitted_u + remainder_u).upper()
    dominant = 'retained_polynomial'
    best = retained_u
    if omitted_u > best:
        dominant, best = 'omitted_polynomial_tail', omitted_u
    if remainder_u > best:
        dominant = 'reciprocal_remainder'
    return {
        'retained_polynomial_upper': pack_upper(retained_u),
        'omitted_polynomial_tail_upper': pack_upper(omitted_u),
        'reciprocal_remainder_upper': pack_upper(remainder_u),
        'total_upper': pack_upper(total),
        'dominant_component': dominant,
    }


def _live_split(retained, omitted, remainder):
    retained_u = _norm_upper(retained, 'retained polynomial')
    omitted_u = _norm_upper(omitted, 'omitted polynomial tail')
    remainder_u = _norm_upper(remainder, 'reciprocal remainder')
    total = (retained_u + omitted_u + remainder_u).upper()
    dominant = 'retained_polynomial'
    best = retained_u
    if omitted_u > best:
        dominant, best = 'omitted_polynomial_tail', omitted_u
    if remainder_u > best:
        dominant = 'reciprocal_remainder'
    return {
        'retained_polynomial': retained_u,
        'omitted_polynomial_tail': omitted_u,
        'reciprocal_remainder': remainder_u,
        'total': total,
        'dominant_component': dominant,
    }


def _recompute_profile_tails(derivative_l1, period_length, settings, max_derivative, *,
                             bits=DEFAULT_BITS):
    with ctx.workprec(bits):
        length = _positive_enclosure(period_length, 'period length')
        l1 = _norm_upper(derivative_l1, 'derivative L1')
        alias, tails = alias_and_tail_bounds(
            l1, length, settings['derivative_order'], settings['quadrature_points'],
            settings['retained_index'], max_derivative)
        return alias, tails


def load_radius_transfer_inputs(*, root=None, bits=DEFAULT_BITS, verify_hash=True,
                                polynomial_order=POLYNOMIAL_ORDER,
                                remainder_moment=MAX_OWNED_ENVELOPE_DEGREE):
    """Authenticate the p=16 W^1..W^4 tables and the radius enclosure."""
    order = _positive_int(polynomial_order, 'polynomial order')
    if order != POLYNOMIAL_ORDER:
        raise ValueError('this owner uses the saved W^1..W^4 tables at order 4')
    bits = _as_bit_precision(bits)
    remainder_moment = _positive_int(remainder_moment, 'remainder moment')
    root = repository_root() if root is None else Path(root)
    table_path = root / FINE_PROFILE_JSON
    radius_path = root / RADIUS_BOUND_JSON
    table = json.loads(table_path.read_text())
    radius = json.loads(radius_path.read_text())
    npz_path = root / table['payload']['path']
    if table['payload']['path'] != FINE_PROFILE_NPZ:
        raise ValueError('fine profile payload path changed')
    digest = (authenticate_profile_payload(table, npz_path)
              if verify_hash else table['payload']['sha256'])
    if digest != AUTHENTIC_PAYLOAD_SHA256:
        raise ValueError('fine profile payload hash changed')
    keys = PROFILE_KEYS[:order]
    meta, coefficients = load_fine_profile_coefficients(npz_path, keys=keys, bits=bits)
    if meta['axial_profile_identity'] != table['axial_profile_identity']:
        raise ValueError('physical profile fingerprint disagrees with the JSON record')
    if meta['axial_profile_identity'] != AXIAL_PROFILE_IDENTITY:
        raise ValueError('physical profile fingerprint changed')
    settings = table['settings']
    if (settings['derivative_order'] != 16 or settings['quadrature_points'] != 65536
            or settings['retained_index'] != 8192):
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
    if not rref > 0:
        raise ValueError('positive reference radius required')
    rmin = _rational_bound(radius['bounds']['radius_lower'])
    eta = _rational_bound(radius['bounds']['radius_ratio_upper'])
    if rmin != rref * (1 - eta):
        raise ValueError('radius lower is not the recorded two-derivative margin')
    sigma = eta * rref / enclosure_w[0]
    amplitude = _hex_float(direction['amplitude'])
    if amplitude == 0:
        raise ValueError('nonzero control amplitude required')
    omitted = []
    proofs = []
    with ctx.workprec(bits):
        for power, key in zip(PROFILE_POWERS[:order], keys):
            proof = meta['profiles'][key]
            if proof['powers'] != [power, 0]:
                raise ValueError('saved W^p table powers changed')
            proofs.append(proof)
            _alias, tails = _recompute_profile_tails(
                proof['derivative_l1_upper'], period, settings, remainder_moment,
                bits=bits)
            omitted.append(tails)
        w_fourier = fourier_sup_derivative_uppers(
            coefficients[PROFILE_KEY], period,
            meta['profiles'][PROFILE_KEY]['derivative_l1_upper'],
            settings['derivative_order'], settings['quadrature_points'],
            settings['retained_index'], remainder_moment, bits=bits)
        delta = tuple(
            (_norm_upper(sigma, 'normal support') * moment).upper()
            for moment in w_fourier['fourier_l1_moments'])
        remainder = reciprocal_remainder_moments(
            delta, rref, polynomial_order=order, max_moment=remainder_moment,
            bits=bits)
        if not remainder['eta'][0] < 1:
            raise ValueError('eta0<1 required')
        tables = tuple(coefficients[key] for key in keys)
        count = len(tables[0])
        if any(len(table) != count for table in tables):
            raise ValueError('matching W^p Fourier inventories required')
        if count != 2 * settings['retained_index'] + 1:
            raise ValueError('retained Fourier inventory must have length 2K+1')
    return {
        'root': root,
        'bits': bits,
        'profile_payload_sha256': digest,
        'axial_profile_identity': meta['axial_profile_identity'],
        'period_length': period,
        'settings': settings,
        'sigma': sigma,
        'reference_radius': rref,
        'radius_lower': rmin,
        'amplitude': amplitude,
        'U_identically_zero': True,
        'alpha_not_applied_twice': True,
        'coefficients': tables,
        'omitted_tails': tuple(omitted),
        'W_fourier_l1': w_fourier['fourier_l1_moments'],
        'delta_r_moments': delta,
        'remainder': remainder,
        'proofs': tuple(proofs),
        'meta': meta,
        'table': table,
        'radius': radius,
        'polynomial_order': order,
    }


def _catalog_for_cutoff(cutoff, bits):
    specs = []
    for n in range(1, 6):
        for s in (0, 1):
            specs.append({
                'kind': 'small',
                'n': n,
                's': s,
                'p': None,
                'weight': lambda w, n=n, s=s: f_small(n, s, w, cutoff, bits=bits),
                'envelope': f_small_envelope(n, s, cutoff, bits=bits),
            })
    for s in (0, 1):
        specs.append({
            'kind': 'shift',
            'n': 0,
            's': s,
            'p': 0,
            'weight': lambda w, s=s: f_shift(s, w, cutoff, bits=bits),
            'envelope': f_shift_envelope(s, cutoff, bits=bits),
        })
    for n, p, s in owned_taylor_parameters():
        specs.append({
            'kind': 'taylor',
            'n': n,
            's': s,
            'p': p,
            'weight': lambda w, n=n, p=p, s=s: f_taylor(n, p, s, w, cutoff, bits=bits),
            'envelope': f_taylor_envelope(n, p, s, cutoff, bits=bits),
        })
    return tuple(specs)


def _spec_key(spec):
    if spec['kind'] == 'small':
        return 'small_n%d_s%d' % (spec['n'], spec['s'])
    if spec['kind'] == 'shift':
        return 'shift_s%d' % spec['s']
    return 'taylor_n%d_p%d_s%d' % (spec['n'], spec['p'], spec['s'])


def evaluate_transfer_catalog(inputs, canonical_cutoff, *, bits=None):
    """Live retained/omitted/remainder splits for every owned kernel at one K."""
    bits = inputs['bits'] if bits is None else _as_bit_precision(bits)
    with ctx.workprec(bits):
        cutoff = _canonical_cutoff(canonical_cutoff, bits=bits)
        specs = _catalog_for_cutoff(cutoff, bits)
        retained = [arb(0)] * len(specs)
        tables = inputs['coefficients']
        count = len(tables[0])
        retained_index = count // 2
        omega = _omega_unit(inputs['period_length'], bits=bits)
        sigma = inputs['sigma']
        rref = inputs['reference_radius']
        order = inputs['polynomial_order']
        pair_env = [arb(0)] * (retained_index + 1)
        for index in range(count):
            mode = index - retained_index
            if mode == 0:
                continue
            values = tuple(table[index] for table in tables)
            pair_env[abs(mode)] += q4_mode_norm_envelope(
                values, sigma, rref, bits=bits, polynomial_order=order)
        for mode in range(1, retained_index + 1):
            envelope = pair_env[mode]
            freq = omega * mode
            for i, spec in enumerate(specs):
                retained[i] += envelope * _norm_upper(spec['weight'](freq), 'mode weight')
        remainder_moments = inputs['remainder']['moments']
        omitted_tails = inputs['omitted_tails']
        catalog = {}
        for spec, retained_sum in zip(specs, retained):
            omitted = omitted_q4_polynomial_tail(
                omitted_tails, sigma, rref, spec['envelope'], bits=bits,
                polynomial_order=order)
            rest = contract_polynomial_envelope(
                remainder_moments, spec['envelope'], bits=bits)
            catalog[_spec_key(spec)] = {
                **_live_split(retained_sum, omitted, rest),
                'kind': spec['kind'],
                'n': spec['n'],
                'p': spec['p'],
                's': spec['s'],
                'envelope_degree': len(spec['envelope']) - 1,
            }
        return {
            'canonical_cutoff_lower': cutoff.lower(),
            'catalog': catalog,
        }


def bound_authentic_radius_transfer_sums(canonical_cutoffs=(256, 1024), *,
                                         root=None, bits=DEFAULT_BITS,
                                         verify_hash=True, inputs=None):
    """Weighted sums on the saved radius control. No field or source run."""
    bits = _as_bit_precision(bits)
    if not canonical_cutoffs:
        raise ValueError('explicit canonical momentum splits required')
    if inputs is None:
        inputs = load_radius_transfer_inputs(
            root=root, bits=bits, verify_hash=verify_hash)
    packed = {}
    with ctx.workprec(bits):
        orders = validate_owned_taylor_envelope_orders()
        for cutoff in canonical_cutoffs:
            evaluated = evaluate_transfer_catalog(inputs, cutoff, bits=bits)
            packed_catalog = {}
            for name, row in evaluated['catalog'].items():
                packed_catalog[name] = {
                    **_pack_split(
                        row['retained_polynomial'], row['omitted_polynomial_tail'],
                        row['reciprocal_remainder']),
                    'kind': row['kind'],
                    'n': row['n'],
                    'p': row['p'],
                    's': row['s'],
                    'envelope_degree': row['envelope_degree'],
                }
            packed[str(cutoff)] = {
                'canonical_cutoff_lower': pack_lower(evaluated['canonical_cutoff_lower']),
                'catalog': packed_catalog,
            }
    return {
        'schema': SCHEMA,
        'profile_payload_sha256': inputs['profile_payload_sha256'],
        'axial_profile_identity': inputs['axial_profile_identity'],
        'period_length_exact': str(inputs['period_length']),
        'normal_support': pack_upper(inputs['sigma']),
        'reference_radius_lower': pack_lower(inputs['reference_radius']),
        'radius_lower': pack_lower(inputs['radius_lower']),
        'history_amplitude_included_in_W': pack_upper(inputs['amplitude']),
        'U_identically_zero': True,
        'alpha_not_applied_twice': True,
        'eta0_upper': pack_upper(inputs['remainder']['eta0_upper']),
        'remainder_majorant': inputs['remainder']['majorant'],
        'not_difference_of_uppers': True,
        'pointwise_rmin_used_as_Wiener_inverse': False,
        'owned_taylor_envelope_orders': orders,
        'kernel_measure': dict(KERNEL_MEASURE),
        'cutoffs': packed,
        'reused': {
            'fine_profile_table': inputs['table']['payload'],
            'radius_bound_schema': inputs['radius']['schema'],
            'loader': 'nsc_ks_high_radius_derivatives.load_fine_profile_coefficients',
            'alias_and_tail_bounds': 'nsc_ks_profile_fourier_bound.alias_and_tail_bounds',
            'fourier_l1': 'nsc_ks_high_radius_derivatives.fourier_sup_derivative_uppers',
        },
        'missing': {
            'time_integral': MISSING_TIME_INTEGRAL,
            'fourth_order_symbol': MISSING_PROJECTOR_SYMBOL,
            'exact_Moyal_remainder': MISSING_MOYAL,
            'Bloch_trace_inequality': MISSING_BLOCH,
            'UV_residual': MISSING_UV,
            'physical_local_gate': PHYSICAL_LOCAL_GATE,
        },
        'scope': {
            'finite_history_geometry_bound': True,
            'source_action_profile_changed': False,
            'field_runs': 0,
            'source_runs': 0,
            'new_profile_fit': False,
            'physical_local_gate': 'OPEN',
            'UV_certificate': False,
            'canonical_K_is_source_energy_cutoff': False,
            'dk_over_2pi_included': False,
            'both_sign_k_included': False,
            'spin_source_multiplicity_included': False,
            'plus_minus_omega_modes_in_inventory': True,
            'known_coefficients_replaced_by_global_An': False,
            'true_reciprocal_remainder_q_minus_q4': True,
        },
    }
