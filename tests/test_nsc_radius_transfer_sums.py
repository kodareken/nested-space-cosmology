"""Directed weighted Fourier sums of the reciprocal-radius perturbation."""
from fractions import Fraction as Q
from pathlib import Path

import pytest
pytest.importorskip('flint')
from flint import acb, arb, arb_series, ctx, fmpq

from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper
from recursive_horizons.nsc_ks_high_radius_derivatives import (
    PERIOD_LENGTH, PROFILE_KEY, authenticate_profile_payload,
    fourier_sup_derivative_uppers, restore_packed,
)
from recursive_horizons.nsc_ks_profile_fourier_bound import alias_and_tail_bounds
from recursive_horizons.nsc_radius_transfer_sums import (
    AUTHENTIC_PAYLOAD_SHA256, AXIAL_PROFILE_IDENTITY, KERNEL_MEASURE,
    MAX_OWNED_ENVELOPE_DEGREE, PHYSICAL_LOCAL_GATE, POLYNOMIAL_ORDER,
    PROFILE_KEYS, SCHEMA, bound_authentic_radius_transfer_sums,
    contract_polynomial_envelope, evaluate_polynomial_envelope,
    f_shift, f_shift_envelope, f_small,
    f_small_envelope, f_taylor, f_taylor_envelope, load_radius_transfer_inputs,
    omitted_q4_polynomial_tail, owned_taylor_parameters, q4_mode_norm_envelope,
    reciprocal_remainder_moments, repository_root, taylor_envelope_degree,
    validate_owned_taylor_envelope_orders, weighted_retained_q4_sum,
)


ROOT = Path(__file__).resolve().parents[1]
LOW_BITS = (8, 12, 30)


def arb_q(value):
    q = Q(value)
    return arb(fmpq(q.numerator, q.denominator))


def covers(upper, exact):
    return not (upper < exact)


def restored_high(value, bits=256):
    with ctx.workprec(bits):
        if isinstance(value, dict):
            return restore_packed(value)
        return restored_upper(exact_upper(value))


@pytest.fixture(scope='module')
def transfer_inputs():
    return load_radius_transfer_inputs(root=ROOT, bits=160, verify_hash=True)


def test_small_kernel_is_the_one_sided_k_integral_without_2pi():
    with ctx.workprec(160):
        k, w = arb(4), arb(2)
        for n in range(1, 6):
            for s in (0, 1):
                peak = k.max(w)
                expected = (w / 2) ** n / (arb(5 - s) * peak ** (5 - s))
                value = f_small(n, s, w, k)
                assert covers(value, expected)
                assert value == expected.upper()
                # Explicit antiderivative of k^{s-6} from peak to infinity.
                integral = arb(1) / (arb(5 - s) * peak ** (5 - s))
                assert covers(value, (w / 2) ** n * integral)
        assert f_small(1, 0, 0, 5) == 0
        assert f_small(5, 1, 0, 5) == 0


def test_shift_and_taylor_kernels_match_unshifted_power_integrals():
    with ctx.workprec(160):
        k, w = arb(3), arb(5)
        for s in (0, 1):
            expected = (w ** (s + 1) - k ** (s + 1)) / arb(s + 1)
            assert covers(f_shift(s, w, k), expected)
            assert f_shift(s, k, k) == 0
            assert f_shift(s, k / 2, k) == 0
        # p>s+1, p=s+1, p<s+1.
        cases = ((0, 2, 0), (0, 1, 0), (0, 1, 1), (2, 4, 1), (4, 5, 0))
        for n, p, s in cases:
            if p == s + 1:
                power_int = (w / k).log()
            else:
                power_int = (w ** (s - p + 1) - k ** (s - p + 1)) / arb(s - p + 1)
            expected = (w / 2) ** n * power_int
            assert covers(f_taylor(n, p, s, w, k), expected)
            assert f_taylor(n, p, s, k, k) == 0
            assert f_taylor(n, p, s, 1, k) == 0
        with pytest.raises(ValueError, match='p=0/n=0 uses f_shift'):
            f_taylor(0, 0, 0, w, k)


def test_threshold_continuity_at_canonical_k():
    with ctx.workprec(160):
        k = arb(7)
        for n in range(1, 6):
            for s in (0, 1):
                at = f_small(n, s, k, k)
                below = (k / 2) ** n / (arb(5 - s) * k ** (5 - s))
                above = (k / 2) ** n / (arb(5 - s) * k ** (5 - s))
                assert covers(at, below)
                assert covers(at, above)
        for s in (0, 1):
            assert f_shift(s, k, k) == 0
            slightly = f_shift(s, k * arb('17/16'), k)
            assert slightly > 0
        for n, p, s in ((0, 1, 0), (1, 3, 1), (4, 9, 0)):
            assert f_taylor(n, p, s, k, k) == 0
            assert f_taylor(n, p, s, k * arb('17/16'), k) > 0


def test_polynomial_envelopes_cover_kernels_on_representatives():
    with ctx.workprec(160):
        k = arb(8)
        points = (arb(0), k / 2, k, k * arb('3/2'), 2 * k, 10 * k)
        for n in range(1, 6):
            for s in (0, 1):
                envelope = f_small_envelope(n, s, k)
                if n <= 5 - s:
                    assert len(envelope) == 1
                else:
                    assert n == 5 and s == 1
                    assert len(envelope) == 2
                for w in points:
                    value = f_small(n, s, w, k)
                    bound = evaluate_polynomial_envelope(envelope, w)
                    assert covers(bound, value)
                    assert not (bound < value)
        for s in (0, 1):
            envelope = f_shift_envelope(s, k)
            assert len(envelope) == s + 2
            for w in points:
                assert covers(evaluate_polynomial_envelope(envelope, w), f_shift(s, w, k))
        for n, p, s in list(owned_taylor_parameters()) + [(0, 1, 1), (3, 3, 0)]:
            envelope = f_taylor_envelope(n, p, s, k)
            assert len(envelope) - 1 == taylor_envelope_degree(n, p, s)
            for w in points:
                assert covers(
                    evaluate_polynomial_envelope(envelope, w),
                    f_taylor(n, p, s, w, k))
        # log(x) <= x at the p=s+1 branch: ln(w/K) <= w/K.
        w = arb(20)
        log_term = f_taylor(0, 1, 0, w, k)
        assert covers(w / k, (w / k).log())
        assert covers(evaluate_polynomial_envelope(f_taylor_envelope(0, 1, 0, k), w), log_term)


def test_owned_taylor_envelope_degree_is_at_most_four():
    orders = validate_owned_taylor_envelope_orders()
    assert orders['max_degree'] <= MAX_OWNED_ENVELOPE_DEGREE
    assert orders['count'] == 50
    for n, p, s in owned_taylor_parameters():
        assert p == n + (p - n)
        assert 1 <= p <= 9
        assert taylor_envelope_degree(n, p, s) <= 4


def test_directed_denominators_reject_collapsed_uppers_and_cover_low_bits():
    with ctx.workprec(160):
        collapsed = exact_upper(arb(2).upper())
        with pytest.raises(ValueError, match='collapsed upper-only'):
            f_small(1, 0, 1, collapsed)
        with pytest.raises(ValueError, match='collapsed upper-only'):
            f_shift(0, 4, collapsed)
        with pytest.raises(ValueError, match='collapsed upper-only'):
            f_taylor(0, 2, 0, 4, collapsed)
        with pytest.raises(ValueError, match='collapsed upper-only'):
            f_small_envelope(1, 0, collapsed)
        with pytest.raises(ValueError, match='positive'):
            f_small(1, 0, 1, 0)
    for bits in LOW_BITS:
        with ctx.workprec(bits):
            k = arb_q(Q(7, 6))
            collapsed = exact_upper(k.upper())
            with pytest.raises(ValueError, match='collapsed upper-only'):
                f_small_envelope(1, 0, collapsed, bits=bits)
            envelope = f_small_envelope(1, 0, Q(7, 6), bits=bits)
        with ctx.workprec(256):
            exact = arb_q(Q(6, 7)) ** 4 / (arb(2) * arb(5))
            assert covers(restored_high(envelope[0]), exact)


def test_single_and_cosine_modes_match_explicit_sums_without_extra_two():
    with ctx.workprec(160):
        # Inventory length 5, retained index 2, cosine at |l|=1.
        zeros = (acb(0),) * 5
        cosine = (acb(0), acb(arb(1) / 2), acb(0), acb(arb(1) / 2), acb(0))
        period = 1
        sigma = 1
        rref = 2
        k = 10
        w = 2 * arb.pi()
        tables = (cosine, zeros, zeros, zeros)

        def weight(freq):
            return f_small(1, 0, freq, k)

        total = weighted_retained_q4_sum(tables, period, sigma, rref, weight)
        mode_env = q4_mode_norm_envelope((acb(arb(1) / 2), 0, 0, 0), sigma, rref)
        expected = 2 * mode_env * f_small(1, 0, w, k)
        assert covers(total, expected)
        assert covers(arb('1e-20') * (1 + expected), abs(total - expected.upper()))
        # A lone DC coefficient is ignored by positive-order transfers.
        dc = (acb(0), acb(0), acb(7), acb(0), acb(0))
        constant = weighted_retained_q4_sum(
            (dc, zeros, zeros, zeros), period, sigma, rref, lambda freq: arb(1))
        assert constant == 0
        single = (acb(0), acb(0), acb(0), acb(arb('3/10')), acb(0))
        total_one = weighted_retained_q4_sum(
            (single, zeros, zeros, zeros), period, sigma, rref, weight)
        env_one = q4_mode_norm_envelope((acb(arb('3/10')), 0, 0, 0), sigma, rref)
        assert covers(total_one, env_one * f_small(1, 0, w, k))


def test_known_coefficients_are_not_replaced_by_global_moments():
    with ctx.workprec(160):
        zeros = (acb(0),) * 5
        high = (acb(0), acb(0), acb(0), acb(0), acb(1))
        period = 1
        k = arb(2)
        w = 4 * arb.pi()
        tables = (high, zeros, zeros, zeros)

        def weight(freq):
            return f_small(1, 0, freq, k)

        sharp = weighted_retained_q4_sum(tables, period, 1, 2, weight)
        env = q4_mode_norm_envelope((acb(1), 0, 0, 0), 1, 2)
        exact = env * f_small(1, 0, w, k)
        loose = env * evaluate_polynomial_envelope(f_small_envelope(1, 0, k), w)
        assert covers(sharp, exact)
        assert loose > sharp
        assert exact < loose


def test_neumann_remainder_series_covers_finite_convolution_and_is_not_a_difference():
    a, b, length = 3, 1, 1
    with ctx.workprec(160):
        coefficients = (acb(0), acb(arb(b) / 2), acb(0), acb(arb(b) / 2), acb(0))
        delta = fourier_sup_derivative_uppers(coefficients, length, 0, 8, 512, 2, 4)
        remainder = reciprocal_remainder_moments(delta['fourier_l1_moments'], a)
        eta0 = arb(b) / arb(a)
        exact_a0 = (eta0 ** 5) / ((1 - eta0) * a)
        assert remainder['not_difference_of_uppers']
        assert not remainder['pointwise_rmin_used_as_Wiener_inverse']
        assert covers(remainder['moments'][0], exact_a0)
        assert remainder['eta0_upper'] < 1
        # A_0 uses only eta0; higher jets must not change it.
        richer = reciprocal_remainder_moments(
            (delta['fourier_l1_moments'][0], 5, 7, 9, 11), a)
        assert covers(arb('1e-18') * remainder['moments'][0],
                      abs(richer['moments'][0] - remainder['moments'][0]))
        # Finite Neumann convolution of the geometric tail plus its remainder.
        x_coeff = {1: arb(b) / (2 * a), -1: arb(b) / (2 * a)}
        xn = {0: arb(1)}
        partial = {}
        for power in range(1, 9):
            nxt = {}
            for i, u in xn.items():
                for j, v in x_coeff.items():
                    nxt[i + j] = nxt.get(i + j, arb(0)) + u * v
            xn = nxt
            if power >= 5:
                for mode, value in xn.items():
                    partial[mode] = partial.get(mode, arb(0)) + value
        tail = (arb(b) / a) ** 9 / (1 - arb(b) / a)
        conv_a0 = (sum(abs(v) for v in partial.values()) + tail) / a
        assert covers(remainder['moments'][0], conv_a0)
        with pytest.raises(ValueError, match='eta0'):
            reciprocal_remainder_moments((2, 0, 0, 0, 0), 1)
        with pytest.raises(ValueError, match='missing moments|through the remainder'):
            reciprocal_remainder_moments((arb('1/10'),), 2, max_moment=4)
        with pytest.raises(ValueError, match='not zero'):
            contract_polynomial_envelope((1,), (1, 2, 3))
        with pytest.raises(ValueError, match='not zero'):
            omitted_q4_polynomial_tail([(1,)], 1, 2, (1, 1))
        with pytest.raises(ValueError, match='cannot be zero'):
            q4_mode_norm_envelope((1, None, 0, 0), 1, 2)


def test_remainder_egf_matches_independent_arb_series():
    with ctx.workprec(160):
        delta = (arb('1/10'), arb('1/5'), arb('1/4'), arb('1/8'), arb('1/16'))
        rref = arb(2)
        report = reciprocal_remainder_moments(delta, rref, polynomial_order=4, max_moment=4)
        inv = (1 / rref).upper()
        eta = [value * inv for value in delta]
        series = arb_series([eta[j] / arb(j).fac() for j in range(5)], 5)
        majorant = (series ** 5) / (arb_series([arb(1)], 5) - series)
        for j in range(5):
            expected = arb(j).fac() * majorant[j] * inv
            assert covers(report['moments'][j], expected)
            assert covers(arb('1e-18') * (1 + expected.upper()),
                          abs(report['moments'][j] - expected.upper()))


def test_kernel_measure_flags_forbid_duplicated_fourier_energy_factors():
    assert KERNEL_MEASURE['dk_over_2pi'] is False
    assert KERNEL_MEASURE['both_sign_k'] is False
    assert KERNEL_MEASURE['spin_source_multiplicity'] is False
    assert KERNEL_MEASURE['plus_minus_omega_in_inventory'] is True
    assert KERNEL_MEASURE['extra_factor_two_for_negative_omega'] is False
    assert KERNEL_MEASURE['canonical_K_is_source_energy_cutoff'] is False
    assert POLYNOMIAL_ORDER == 4
    assert PROFILE_KEYS == ('profile_1_0', 'profile_2_0', 'profile_3_0', 'profile_4_0')


def test_authentic_payload_hash_recomputed_tails_and_control_geometry(transfer_inputs):
    table = transfer_inputs['table']
    digest = authenticate_profile_payload(table, ROOT / table['payload']['path'])
    assert digest == table['payload']['sha256'] == AUTHENTIC_PAYLOAD_SHA256
    assert transfer_inputs['profile_payload_sha256'] == digest
    assert transfer_inputs['axial_profile_identity'] == AXIAL_PROFILE_IDENTITY
    assert transfer_inputs['period_length'] == PERIOD_LENGTH
    assert transfer_inputs['U_identically_zero']
    assert transfer_inputs['alpha_not_applied_twice']
    assert transfer_inputs['amplitude'] == Q(float.fromhex(
        transfer_inputs['meta']['profile_description']['directions'][0]['amplitude']))
    with ctx.workprec(160):
        amplitude = arb_q(transfer_inputs['amplitude'])
        assert amplitude == arb(0.001)
        w0 = transfer_inputs['W_fourier_l1'][0]
        sigma = arb_q(transfer_inputs['sigma'])
        rref = arb_q(transfer_inputs['reference_radius'])
        assert w0 > amplitude / 2
        # Amplitude already sits in W; sigma is the normal window, not another alpha.
        rebuilt_delta = (sigma * w0).upper()
        assert covers(arb('1e-18') * rebuilt_delta,
                      abs(transfer_inputs['delta_r_moments'][0] - rebuilt_delta))
        rebuilt_eta = (transfer_inputs['delta_r_moments'][0] / rref).upper()
        assert covers(arb('1e-18') * rebuilt_eta,
                      abs(transfer_inputs['remainder']['eta0_upper'] - rebuilt_eta))
        assert transfer_inputs['remainder']['eta0_upper'] < 1
        settings = transfer_inputs['settings']
        for proof, tails in zip(transfer_inputs['proofs'], transfer_inputs['omitted_tails']):
            alias, recomputed = alias_and_tail_bounds(
                restored_upper(proof['derivative_l1_upper']),
                arb_q(PERIOD_LENGTH), settings['derivative_order'],
                settings['quadrature_points'], settings['retained_index'], 4)
            assert alias >= 0
            for j in range(3):
                stored = restored_upper(proof['tail_upper'][j])
                assert covers(arb('1e-14') * stored, abs(recomputed[j] - stored))
                assert covers(arb('1e-14') * stored, abs(tails[j] - stored))
            assert len(tails) == 5
            # JSON measurements are binary64 displays; the bound uses recomputed dyadics.
            assert 'mantissa' in proof['tail_upper'][0]


def test_authentic_k256_and_k1024_splits_retained_dominates(transfer_inputs):
    bound = bound_authentic_radius_transfer_sums(
        (256, 1024), root=ROOT, bits=160, inputs=transfer_inputs)
    assert bound['schema'] == SCHEMA
    assert bound['profile_payload_sha256'] == AUTHENTIC_PAYLOAD_SHA256
    assert bound['scope']['physical_local_gate'] == 'OPEN'
    assert bound['missing']['physical_local_gate'] == PHYSICAL_LOCAL_GATE
    assert bound['scope']['UV_certificate'] is False
    assert bound['scope']['canonical_K_is_source_energy_cutoff'] is False
    assert bound['scope']['field_runs'] == 0
    assert bound['scope']['source_runs'] == 0
    assert bound['scope']['source_action_profile_changed'] is False
    assert bound['scope']['dk_over_2pi_included'] is False
    assert bound['scope']['both_sign_k_included'] is False
    assert bound['scope']['known_coefficients_replaced_by_global_An'] is False
    assert bound['scope']['true_reciprocal_remainder_q_minus_q4'] is True
    assert bound['kernel_measure'] == KERNEL_MEASURE
    assert bound['not_difference_of_uppers']
    assert bound['U_identically_zero']
    assert bound['alpha_not_applied_twice']
    representatives = (
        'small_n1_s0', 'small_n5_s1', 'shift_s0', 'shift_s1',
        'taylor_n0_p1_s0', 'taylor_n4_p5_s1',
    )
    with ctx.workprec(160):
        for cutoff in ('256', '1024'):
            catalog = bound['cutoffs'][cutoff]['catalog']
            assert set(catalog) >= set(representatives)
            assert len(catalog) == 10 + 2 + 50
            for name in representatives:
                row = catalog[name]
                retained = restore_packed(row['retained_polynomial_upper'])
                omitted = restore_packed(row['omitted_polynomial_tail_upper'])
                rest = restore_packed(row['reciprocal_remainder_upper'])
                total = restore_packed(row['total_upper'])
                assert row['total_upper']['kind'] == 'dyadic'
                assert covers(total, retained + omitted + rest)
                assert row['dominant_component'] == 'retained_polynomial'
                assert retained > omitted
                assert retained > rest
                assert omitted >= 0 and rest >= 0
        for name in representatives:
            left = restore_packed(
                bound['cutoffs']['1024']['catalog'][name]['total_upper'])
            right = restore_packed(
                bound['cutoffs']['256']['catalog'][name]['total_upper'])
            assert not (left > right)


def test_repository_root_matches_this_checkout():
    assert repository_root() == ROOT


def test_remainder_moments_keep_their_order_and_restore_an_outer_series_context():
    from recursive_horizons.nsc_radius_transfer_sums import reciprocal_remainder_moments
    previous=ctx.cap
    try:
        ctx.cap=10
        expected=reciprocal_remainder_moments([.1,.2,.3,.4,.5],2,max_moment=4)['moments']
        ctx.cap=1
        actual=reciprocal_remainder_moments([.1,.2,.3,.4,.5],2,max_moment=4)['moments']
        assert ctx.cap==1
        assert all(a==b and a>0 for a,b in zip(actual,expected))
    finally:
        ctx.cap=previous
