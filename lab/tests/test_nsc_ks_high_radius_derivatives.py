"""Spatial supremum and Fourier-l1 reciprocal-radius bounds."""
from fractions import Fraction as Q
from math import comb, factorial
from pathlib import Path
import json

import pytest
pytest.importorskip('flint')
from flint import acb, arb, ctx, fmpq

from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper
from recursive_horizons.nsc_ks_high_radius_derivatives import (
    MAX_REQUESTED_DERIVATIVE, PERIOD_LENGTH, PHYSICAL_LOCAL_GATE, PROFILE_KEY,
    authenticate_profile_payload, bound_pure_radius_difference,
    compact_support_norms, fourier_sup_derivative_uppers,
    load_fine_profile_coefficients, pack_upper, perturbed_radius_derivative_uppers,
    positive_reciprocal_upper, reciprocal_difference_derivative_uppers,
    reciprocal_fourier_moment_uppers, reciprocal_full_derivative_uppers,
    repository_root, restore_packed,
)
from recursive_horizons.nsc_ks_profile_fourier_bound import alias_and_tail_bounds


ROOT = Path(__file__).resolve().parents[1]
LOW_BITS = (8, 12, 30)
DENOMINATOR = Q(7, 6)
RECIPROCAL = Q(6, 7)


def arb_q(value):
    q = Q(value)
    return arb(fmpq(q.numerator, q.denominator))


def covers(upper, exact):
    """Directed upper: flint a>=a is false when both balls have radius."""
    return not (upper < exact)


def restored_high(value, bits=256):
    with ctx.workprec(bits):
        if isinstance(value, dict):
            return restore_packed(value)
        return restored_upper(exact_upper(value))


def certainly_covers_rational(upper, exact, bits=256):
    with ctx.workprec(bits):
        u = restored_high(upper, bits)
        t = arb_q(exact)
        if u < t:
            return False
        return True


def test_reciprocal_product_rule_covers_explicit_low_degree_derivative():
    # r=3+x on [-1,1] has r_min=2, r'=1, and vanishing higher spatial derivatives.
    with ctx.workprec(160):
        delta = (arb(1), arb(1), arb(0), arb(0), arb(0))
        q_full = reciprocal_full_derivative_uppers(delta, 2)
        q_diff = reciprocal_difference_derivative_uppers(delta, 2, 3)
        assert covers(q_diff[0], arb_q(Q(1, 6)))
        for n, bound in enumerate(q_full[:4]):
            actual = factorial(n) / Q(2) ** (n + 1)
            assert covers(bound, arb_q(actual))
            if n >= 1:
                assert covers(q_diff[n], bound)


def test_zero_perturbation_vanishes_and_rejects_j_at_p_minus_one():
    with ctx.workprec(160):
        zeros = (0, 0, 0, 0)
        q = reciprocal_difference_derivative_uppers(zeros, Q(7, 5), Q(7, 5))
        assert all(value == 0 for value in q)
        radius = perturbed_radius_derivative_uppers((0, 0, 0), (0, 0, 0), Q(3, 100))
        assert all(value == 0 for value in radius)
        moments = reciprocal_fourier_moment_uppers(zeros, Q(7, 5))
        assert moments['q_difference'][0] == 0
        assert all(value == 0 for value in moments['q_difference'])
        assert covers(moments['q_full'][0], arb_q(Q(5, 7)))
        assert all(value == 0 for value in moments['q_full'][1:])
        alias_and_tail_bounds(0, 1, 16, 65536, 8192, 14)
        with pytest.raises(ValueError, match='p-1'):
            fourier_sup_derivative_uppers([acb(0)] * 5, 1, 0, 8, 512, 2, 7)
        with pytest.raises(ValueError, match='integrable'):
            alias_and_tail_bounds(0, 1, 16, 65536, 8192, 15)


def test_single_cosine_mode_sup_is_sharp_with_zero_omitted_tail():
    with ctx.workprec(160):
        coefficients = (acb(0), acb(arb(1) / 2), acb(0), acb(arb(1) / 2), acb(0))
        result = fourier_sup_derivative_uppers(coefficients, 1, 0, 8, 512, 2, 2)
        assert result['omitted_tail_uppers'][0] == 0
        assert result['fourier_l1_moments'] == result['sup_uppers']
        assert covers(result['sup_uppers'][0], arb(1))
        assert covers(result['sup_uppers'][1], 2 * arb.pi())
        assert covers(result['sup_uppers'][2], (2 * arb.pi()) ** 2)


def test_l1_l2_use_support_width_once():
    with ctx.workprec(160):
        width = Q(3, 25)
        norms = compact_support_norms((arb(2), arb(5)), width)
        assert covers(norms[0]['L1'], arb_q(width) * 2)
        assert covers(norms[0]['L2_squared'], arb_q(width) * 4)
        assert covers(norms[1]['L1'], arb_q(width) * 5)
        assert covers(norms[1]['L2_squared'], arb_q(width) * 25)
        assert abs(restore_packed(pack_upper(norms[0]['L1'])) - (arb_q(width) * 2).upper()) < arb(2) ** -140


@pytest.mark.parametrize('bits', LOW_BITS)
def test_low_bit_reciprocal_covers_exact_and_exposes_collapsed_undercoverage(bits):
    with ctx.workprec(bits):
        enclosure = arb_q(DENOMINATOR)
        collapsed = (1 / enclosure.upper()).upper()
        fixed = positive_reciprocal_upper(DENOMINATOR, bits=bits)
        q_full = reciprocal_full_derivative_uppers((0,), DENOMINATOR, bits=bits)[0]
        omega = fourier_sup_derivative_uppers(
            (acb(0), acb(arb(1) / 2), acb(0), acb(arb(1) / 2), acb(0)),
            DENOMINATOR, 0, 8, 512, 2, 1, bits=bits)['fourier_l1_moments'][1]
    with ctx.workprec(256):
        true_inv = arb_q(RECIPROCAL)
        true_omega = 2 * arb.pi() * arb_q(RECIPROCAL)
        collapsed_high = restored_high(collapsed)
        assert collapsed_high < true_inv
        assert certainly_covers_rational(fixed, RECIPROCAL)
        assert certainly_covers_rational(q_full, RECIPROCAL)
        assert covers(restored_high(omega), true_omega)
        with pytest.raises(ValueError, match='collapsed upper-only'):
            positive_reciprocal_upper(exact_upper(enclosure.upper()), bits=bits)
        with pytest.raises(ValueError, match='collapsed upper-only'):
            reciprocal_full_derivative_uppers((0,), exact_upper(enclosure.upper()), bits=bits)
        with pytest.raises(ValueError, match='collapsed upper-only'):
            fourier_sup_derivative_uppers(
                (acb(0), acb(0), acb(0), acb(0), acb(0)),
                exact_upper(enclosure.upper()), 0, 8, 512, 2, 1, bits=bits)


def test_neumann_margin_rejection_and_zero_fourier_inverse():
    with ctx.workprec(160):
        with pytest.raises(ValueError, match='Neumann radius margin'):
            reciprocal_fourier_moment_uppers((1,), Q(1, 2))
        with pytest.raises(ValueError, match='Neumann radius margin'):
            reciprocal_fourier_moment_uppers((arb_q(Q(7, 5)), 0), Q(7, 5))
        zero = reciprocal_fourier_moment_uppers((0, 0, 0), 4)
        assert zero['q_difference'][0] == 0
        assert covers(zero['q_full'][0], arb_q(Q(1, 4)))


def _eulerian_polynomial(n, x):
    # A(n,k) Eulerian: sum_{k>=1} k^n x^k = x A_n(x) / (1-x)^{n+1}, n>=1.
    counts = [1]
    for m in range(2, n + 1):
        nxt = [0] * m
        for k, value in enumerate(counts):
            nxt[k] += (k + 1) * value
            if k + 1 < m:
                nxt[k + 1] += (m - 1 - k) * value
        counts = nxt
    return sum(counts[k] * x ** k for k in range(n))


def _cosine_reciprocal_moments(a, b, length, orders, bits=256):
    """Exact Wiener moments of 1/(a+b cos(2 pi z/L)) from its Fourier series."""
    with ctx.workprec(bits):
        a, b, length = arb(a), arb(b), arb(length)
        disc = (a * a - b * b).sqrt()
        lam = (a - disc) / b
        omega = 2 * arb.pi() / length
        moments = [(1 / disc) * (1 + lam) / (1 - lam)]
        for n in range(1, orders):
            series = lam * _eulerian_polynomial(n, lam) / (1 - lam) ** (n + 1)
            moments.append(2 * (omega ** n) / disc * series)
        return tuple(moments)


def test_cosine_fourier_inverse_and_neumann_convolution():
    a, b, length = 3, 1, 1
    with ctx.workprec(160):
        coefficients = (acb(0), acb(arb(b) / 2), acb(0), acb(arb(b) / 2), acb(0))
        delta = fourier_sup_derivative_uppers(coefficients, length, 0, 8, 512, 2, 4)
        moments = reciprocal_fourier_moment_uppers(delta['fourier_l1_moments'], a)
        assert covers(moments['q_full'][0], arb_q(Q(1, 2)))
        assert covers(moments['q_difference'][0], arb_q(Q(1, 6)))
        # Finite Neumann convolution plus remainder still uses the full reciprocal.
        x_coeff = {1: arb(b) / (2 * a), -1: arb(b) / (2 * a)}
        state = {0: arb(1)}
        partial = dict(state)
        xn = dict(state)
        for _ in range(1, 9):
            nxt = {}
            for i, u in xn.items():
                for j, v in x_coeff.items():
                    nxt[i + j] = nxt.get(i + j, arb(0)) + u * v
            xn = nxt
            sign = arb(-1) ** _
            for k, value in xn.items():
                partial[k] = partial.get(k, arb(0)) + sign * value
        remainder = (arb(b) / a) ** 9 / (1 - arb(b) / a)
        conv_A0 = (sum(abs(v) for v in partial.values()) + remainder) / a
        assert covers(moments['q_full'][0], conv_A0)
    exact = _cosine_reciprocal_moments(a, b, length, 5)
    with ctx.workprec(256):
        for n, bound in enumerate(moments['q_full']):
            assert covers(restored_high(bound), exact[n])
            if n >= 1:
                assert restored_high(moments['q_difference'][n]) == restored_high(bound)


def test_authentic_table_hash_binary_enclosures_and_radius_orders_0_2():
    table = json.loads((ROOT / 'results/development/nsc-ks-fine-profile-table.json').read_text())
    radius = json.loads((ROOT / 'results/development/nsc-ks-radius-bound.json').read_text())
    npz = ROOT / table['payload']['path']
    digest = authenticate_profile_payload(table, npz)
    assert digest == table['payload']['sha256'] == 'fc20111d0c72f1ac81d56ee5ea64cdd1bb94a9bfc6e86edfa19d1f8992f0ce33'
    meta, coefficients = load_fine_profile_coefficients(npz, keys=(PROFILE_KEY,), bits=160)
    assert meta['axial_profile_identity'] == table['axial_profile_identity']
    assert meta['axial_profile_identity'] == 'a4c061f297027bcb9934fa2338bb81fb799ca4ec26a025f6bc563384c70cd69d'
    assert set(coefficients) == {PROFILE_KEY}
    balls = coefficients[PROFILE_KEY]
    K = table['settings']['retained_index']
    assert len(balls) == 2 * K + 1
    with ctx.workprec(160):
        assert balls[K].imag.contains(0)
        assert balls[K + 1].overlaps(balls[K - 1].conjugate())
        assert restored_upper(meta['profiles'][PROFILE_KEY]['derivative_l1_upper']) > 0
    bound = bound_pure_radius_difference(root=ROOT, max_derivative=9, bits=160, verify_hash=True)
    assert bound['profile_payload_sha256'] == digest
    assert bound['period_length_exact'] == str(PERIOD_LENGTH)
    assert bound['schema'] == 'NSC-KS-HIGH-RADIUS-DERIVATIVES-v2'
    assert bound['U_identically_zero']
    assert bound['alpha_not_applied_twice']
    assert bound['measure_not_applied_twice']
    assert bound['derivative_orders'] == list(range(MAX_REQUESTED_DERIVATIVE + 1))
    assert bound['scope']['physical_local_gate'] == 'OPEN'
    assert bound['scope']['UV_certificate'] is False
    assert bound['scope']['field_runs'] == 0
    assert bound['scope']['source_action_Gamma_rest_changed'] is False
    assert bound['scope']['band_action_required_for_covariance_UV'] is False
    assert bound['scope']['spatial_supremum_and_Fourier_l1_only'] is True
    assert bound['missing']['physical_local_gate'] == PHYSICAL_LOCAL_GATE
    assert bound['missing']['normal_time_derivatives'].startswith('MISSING')
    assert bound['missing']['fourth_order_symbol_high_derivatives'].startswith('OPEN')
    assert bound['missing']['exact_Moyal_remainder'].startswith('OPEN')
    assert bound['missing']['Bloch_trace_inequality'].startswith('OPEN')
    assert bound['missing']['P4_remainder'].startswith('OPEN')
    assert 'band-action' not in bound['missing']['UV_residual']
    assert 'Gamma_sub' in bound['missing']['UV_residual']
    with ctx.workprec(160):
        width = restore_packed(bound['axial_support_width'])
        period = arb_q(PERIOD_LENGTH)
        assert width < period
        outer = arb(float.fromhex(meta['profile_description']['directions'][0]['w']['outer']))
        assert width == 2 * outer
        sigma = restore_packed(bound['normal_support'])
        amplitude = restore_packed(bound['history_amplitude_included_in_W'])
        w0 = restore_packed(bound['W_fourier_l1'][0])
        assert amplitude == arb(0.001)
        assert w0 > amplitude / 2
        assert covers(restore_packed(bound['delta_r_sup'][0]), sigma * w0)
        assert restore_packed(bound['delta_r_A0']) == restore_packed(bound['delta_r_sup'][0])
        rref = restore_packed(bound['reference_radius_lower'])
        rmin = restore_packed(bound['radius_lower'])
        assert bound['radius_lower']['exact_rational'] == radius['bounds']['radius_lower']['exact_rational']
        assert bound['reference_radius_lower']['exact_rational'] == '7/5'
        assert covers(rref - restore_packed(bound['delta_r_sup'][0]), rmin)
        assert covers(restore_packed(bound['q_full_sup'][0]), 1 / rmin)
        eta = Q(radius['bounds']['radius_ratio_upper']['exact_rational'])
        rmin_q = Q(bound['radius_lower']['exact_rational'])
        assert not (restore_packed(bound['q_difference_sup'][0]) > arb_q(eta / rmin_q))
        rebuilt = restore_packed(bound['delta_r_sup'][0]) / (rref * rmin)
        q0 = restore_packed(bound['q_difference_sup'][0])
        assert covers(arb('1e-12') * rebuilt, abs(q0 - rebuilt))
        d0 = restore_packed(bound['delta_r_A0'])
        margin = restore_packed(bound['neumann_margin_lower'])
        assert bound['pointwise_radius_lower_is_conservative_for_neumann'] is True
        assert not (rmin > margin)
        neumann_inv = 1 / (arb_q('7/5') - d0)
        assert covers(restore_packed(bound['q_full_fourier_l1'][0]), neumann_inv)
        assert covers(restore_packed(bound['q_difference_fourier_l1'][0]),
                      d0 / (arb_q('7/5') * (arb_q('7/5') - d0)))
        assert restore_packed(bound['q_full_fourier_l1'][0]) != restore_packed(bound['q_full_sup'][0])
        assert restore_packed(bound['q_difference_fourier_l1'][0]) != restore_packed(bound['q_difference_sup'][0])
        # Amplitude already sits in physical W; sigma is the normal window.
        assert d0 < sigma * amplitude
        assert d0 > sigma * amplitude / 2
        for j in range(3):
            fourier = restore_packed(bound['W_fourier_l1'][j])
            enclosure = arb_q(radius['bounds']['w_sup_orders_0_1_2'][j]['exact_rational'])
            assert not (fourier > enclosure)
            stored = restored_upper(meta['profiles'][PROFILE_KEY]['tail_upper'][j])
            omitted = restore_packed(bound['omitted_profile_tails'][j])
            assert covers(arb('1e-14') * stored, abs(omitted - stored))
            assert bound['omitted_profile_tails'][j]['kind'] == 'dyadic'
            assert 'mantissa' in bound['omitted_profile_tails'][j]
        for j, row in enumerate(bound['q_difference_sup']):
            sup = restore_packed(row)
            l1 = restore_packed(bound['q_difference_L1'][j])
            l2sq = restore_packed(bound['q_difference_L2_squared'][j])
            assert covers(arb('1e-12') * abs(width * sup), abs(l1 - width * sup))
            assert covers(arb('1e-12') * abs(width * sup * sup), abs(l2sq - width * sup * sup))
            if j >= 1:
                assert row == bound['q_full_sup'][j]
                assert bound['q_difference_fourier_l1'][j] == bound['q_full_fourier_l1'][j]
        assert Q(bound['axial_support_width']['exact_rational']) != PERIOD_LENGTH


def test_repository_root_matches_this_checkout():
    assert repository_root() == ROOT
