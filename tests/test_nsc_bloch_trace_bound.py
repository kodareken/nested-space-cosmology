"""Finite Bloch-mode, directed-denominator, Holder/Parseval and missing-input checks."""
from fractions import Fraction as Q

import numpy as np
import pytest
import sympy as sp

pytest.importorskip('flint')
from flint import arb, ctx, fmpq

from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper
from recursive_horizons.nsc_bloch_trace_bound import (
    I_AXIAL_DERIVATIVE_ORDER, MAXIMUM_MIXED_GEOMETRY_ORDER,
    bloch_wave_numbers, bloch_period_constants, fiber_evaluation_weights,
    weyl_bloch_matrix, fiber_H, fiber_density, fiber_current, matrix_nuclear,
    operator_hs, operator_nuclear, hs_holder_nuclear_from_H,
    shifted_weight_product_upper, integrated_H_from_symbol_I,
    integrated_nuclear_from_symbol_I, h2_weighted_comparison_matrix,
    duhamel_sobolev_multipliers, propagate_nuclear_seminorms,
    density_current_from_nuclear, local_stress_from_traces,
    conditional_bloch_trace_bound,
)

I_2 = np.eye(2, dtype=complex)
SIGMA_1 = np.array([[0, 1], [1, 0]], complex)


def _covers(item, rational):
    left = restored_upper(item) if isinstance(item, dict) else item
    expected = arb(fmpq(rational.numerator, rational.denominator)) if isinstance(rational, Q) else arb(rational)
    return left.lower() >= expected.upper() or left >= expected


def test_short_period_q_endpoint_constants_reject_sqrt_L_over_2():
    L = Q(2, 5)
    constants = bloch_period_constants(L)
    n = np.arange(-80, 81)
    with ctx.workprec(200):
        assert _covers(constants['C0_squared_upper'], Q(26, 25))
        assert _covers(constants['C1_squared_upper'], Q(29, 100))
        assert _covers(constants['BZ_measure_dq_over_2pi'], Q(5, 2))
        assert restored_upper(constants['C0_squared_upper']) > arb(fmpq(1, 2))
        assert not constants['large_L_asymptotics_used']
    for q in (0.0, np.pi / float(L), -np.pi / float(L)):
        k = bloch_wave_numbers(L, q, n)
        c0sq, c1sq = fiber_evaluation_weights(k)
        with ctx.workprec(200):
            assert arb(c0sq) <= restored_upper(constants['C0_squared_upper'])
            assert arb(c1sq) <= restored_upper(constants['C1_squared_upper'])
    k0 = bloch_wave_numbers(L, 0.0, n)
    assert fiber_evaluation_weights(k0)[0] >= 1.0
    assert fiber_evaluation_weights(k0)[0] > float(L) / 2
    assert sp.simplify(4 * (sp.pi**2 / 6 - sp.pi**2 / 24) - sp.pi**2 / 2) == 0


def test_period_and_geometry_denominators_keep_enclosing_intervals():
    for bits in (8, 12, 30):
        with ctx.workprec(bits):
            L = arb(arb(fmpq(2, 5)), arb(fmpq(1, 20)))
            a = arb(arb(fmpq(4, 5)), arb(fmpq(1, 10)))
            r = arb(arb(fmpq(1, 2)), arb(fmpq(1, 10)))
            wrong_L = (arb(1) / L.upper()).upper()
            wrong_a = (arb(1) / a.upper()).upper()
            wrong_r = (arb(1) / r.upper()).upper()
            constants = bloch_period_constants(L, bits=bits)
            heat = integrated_H_from_symbol_I(1, L, bits=bits)
            stress_a = local_stress_from_traces(
                0, 1, mass=0, absolute_angular=0, axial_lower=a,
                radius_lower=1, multiplicity=1, bits=bits)
            stress_r = local_stress_from_traces(
                1, 0, mass=0, absolute_angular=1, axial_lower=1,
                radius_lower=r, multiplicity=1, bits=bits)
            wrong_H = (arb(8).sqrt() / L.upper()).upper()
        with ctx.workprec(200):
            true_reciprocal = arb(fmpq(5, 2))
            true_c0 = arb(1) + arb(fmpq(1, 25))
            true_H = arb(8).sqrt() * arb(fmpq(5, 2))
            true_a = arb(fmpq(5, 4))
            true_r = arb(2)
            assert _covers(constants['BZ_measure_dq_over_2pi'], Q(5, 2))
            assert _covers(constants['C0_squared_upper'], Q(26, 25))
            assert _covers(heat['integrated_H_upper'], true_H)
            assert _covers(stress_a['N_upper'], Q(5, 4))
            assert _covers(stress_r['N_upper'], 2)
            assert not (wrong_L.lower() >= true_reciprocal.upper() or wrong_L >= true_reciprocal)
            assert not (wrong_a.lower() >= true_a.upper() or wrong_a >= true_a)
            assert not (wrong_r.lower() >= true_r.upper() or wrong_r >= true_r)
            assert not (wrong_H.lower() >= true_H.upper() or wrong_H >= true_H)
    with pytest.raises(ValueError, match='collapsed upper'):
        bloch_period_constants(exact_upper(arb(1)))
    with pytest.raises(ValueError, match='collapsed upper'):
        local_stress_from_traces(0, 1, mass=0, absolute_angular=0,
                                 axial_lower=exact_upper(arb(1)), radius_lower=1, multiplicity=1)


def test_shifted_weight_product_inequality():
    A, B = sp.symbols('A B', nonnegative=True)
    remainder = sp.expand(8 * (A**4 + B**4) - (A + B)**4)
    assert sp.expand(remainder - (A - B)**2 * (7 * A**2 + 10 * A * B + 7 * B**2)) == 0
    k, s = sp.symbols('k s', real=True)
    product = (1 + (k + s)**2) * (1 + (k - s)**2)
    square = (1 + k**2 + s**2)**2
    assert sp.expand(square - product) == 4 * k**2 * s**2
    with ctx.workprec(120):
        for kv, sv in ((Q(0), Q(0)), (Q(1, 2), Q(3, 5)), (Q(2), Q(-1)), (Q(9, 2), Q(9, 2))):
            actual = (1 + (kv + sv)**2) * (1 + (kv - sv)**2)
            actual_sq = actual * actual
            bound = restored_upper(shifted_weight_product_upper(kv, sv))
            assert bound >= arb(fmpq(actual_sq.numerator, actual_sq.denominator))
        equal = restored_upper(shifted_weight_product_upper(0, 1))
        assert equal.lower() >= arb(16).upper() or equal >= arb(16)


def test_parseval_H_identity_and_reconstruction_have_no_extra_1_over_L():
    L = 2.0
    matrices = {0: I_2, 1: 0.3 * SIGMA_1, -1: 0.3 * SIGMA_1}

    def a_hat(l, k):
        block = matrices.get(int(l), np.zeros((2, 2), complex))
        return np.exp(-float(k) ** 2) * np.asarray(block, complex)

    n_values = np.arange(-20, 21)
    q_vals = np.linspace(-np.pi / L, np.pi / L, 241, endpoint=False)
    dq = q_vals[1] - q_vals[0]
    z = 0.3
    H2 = 0.0
    kernel = np.zeros((2, 2), complex)
    for q in q_vals:
        blocks, k = weyl_bloch_matrix(a_hat, L, q, n_values)
        H2 += fiber_H(blocks, k) ** 2 * dq / (2 * np.pi)
        kernel += fiber_density(blocks, k, z) * dq / (2 * np.pi)

    k_grid = np.linspace(-7, 7, 2801)
    dk = k_grid[1] - k_grid[0]
    rhs_H2 = 0.0
    direct = np.zeros((2, 2), complex)
    I_value = 0.0
    for l in (-1, 0, 1):
        omega = 2 * np.pi * l / L
        for k in k_grid:
            hat = a_hat(l, k)
            fro2 = np.linalg.norm(hat, 'fro') ** 2
            s = omega / 2
            rhs_H2 += ((1 + (k + s) ** 2) ** 2) * ((1 + (k - s) ** 2) ** 2) * fro2 * dk / (2 * np.pi)
            direct += np.exp(1j * omega * z) * hat * dk / (2 * np.pi)
            I_value += ((1 + k * k) ** 4 * fro2 + (omega ** 8) * fro2 / 256) * L * dk / (2 * np.pi)
    assert abs(H2 - rhs_H2) / max(rhs_H2, 1e-30) < 3e-3
    assert np.max(np.abs(kernel - direct)) < 3e-3
    assert H2 <= (8 / L) * I_value + 1e-6
    blocks, k = weyl_bloch_matrix(a_hat, L, 0.1, np.array([-1, 0, 1]))
    i0, i1 = 1, 2
    assert np.allclose(blocks[i1, i0], a_hat(1, 0.5 * (k[i1] + k[i0])))


def test_noninteger_and_duplicate_modes_are_rejected():
    with pytest.raises(ValueError, match='integer'):
        bloch_wave_numbers(1, 0, [0, 1.5])
    with pytest.raises(ValueError, match='integer'):
        bloch_wave_numbers(1, 0, np.array([0.9, 1.0]))
    with pytest.raises(ValueError, match='integer'):
        bloch_wave_numbers(1, 0, [True, 1])
    with pytest.raises(ValueError, match='distinct'):
        bloch_wave_numbers(1, 0, [0, 1, 1])
    with pytest.raises(ValueError, match='distinct'):
        weyl_bloch_matrix(lambda l, k: I_2, 1, 0, np.array([0, 0]))
    k = bloch_wave_numbers(Q(2, 5), 0, [Q(-1), 0, 1])
    assert k.shape == (3,)


def test_hs_holder_spin_factor_and_svd_density_current_contractions():
    blocks = I_2.reshape(1, 1, 2, 2)
    k = np.array([0.0])
    H = fiber_H(blocks, k)
    S = np.sqrt(fiber_evaluation_weights(k)[0])
    d11 = operator_nuclear(blocks, k, 1, 1)
    assert H == pytest.approx(np.sqrt(2))
    assert d11 == pytest.approx(2)
    assert H * S < d11
    bound = hs_holder_nuclear_from_H(H, S, spin_dim=2)
    with ctx.workprec(120):
        assert restored_upper(bound['d11_upper']) >= arb(2)
        assert restored_upper(bound['d21_upper']) >= arb(2)
    rng = np.random.default_rng(11)
    n_values = np.arange(-2, 3)
    k = bloch_wave_numbers(Q(5, 2), 0.2, n_values)
    raw = rng.normal(size=(5, 5, 2, 2)) + 1j * rng.normal(size=(5, 5, 2, 2))
    hermitian = raw + np.transpose(raw, (1, 0, 3, 2)).conj()
    H = fiber_H(hermitian, k)
    c0sq, c1sq = fiber_evaluation_weights(k)
    S = np.sqrt(c0sq)
    C0, C1 = S, np.sqrt(c1sq)
    d11 = operator_nuclear(hermitian, k, 1, 1)
    d21 = operator_nuclear(hermitian, k, 2, 1)
    mid = operator_hs(hermitian, k, 1, 2)
    assert d11 <= mid * np.sqrt(2) * S + 1e-9
    assert d11 <= H * np.sqrt(2) * S + 1e-9
    assert d21 <= H * np.sqrt(2) * S + 1e-9
    for z in (0.0, 0.4, 1.1):
        density = matrix_nuclear(fiber_density(hermitian, k, z))
        current = fiber_current(hermitian, k, z)
        assert density <= d11 * c0sq + 1e-9
        assert matrix_nuclear(current) <= d21 * C0 * C1 + 1e-9
        left = np.einsum('m,m,mnab,n->ab',
                         np.exp(1j * k * z), k, hermitian, np.exp(-1j * k * z))
        assert np.allclose(0.5 * (left + left.conj().T), current)


def test_duhamel_comparison_matrix_and_serialized_norm_records():
    K1, K2 = Q(1, 5), Q(3, 10)
    matrix = h2_weighted_comparison_matrix(K1, K2)
    with ctx.workprec(200):
        root2 = arb(2).sqrt()
        assert matrix[1][0].upper() >= (root2 * arb(fmpq(1, 5))).lower()
        multipliers = duhamel_sobolev_multipliers(K1, K2)
        assert _covers(multipliers['M1_upper'], Q(6, 5))
    record = exact_upper(arb(fmpq(1, 7)))
    evolved = propagate_nuclear_seminorms(record, 0, record, 0, 0, 0)
    with ctx.workprec(200):
        assert _covers(evolved['d11_evolved_upper'], Q(1, 7))
        assert _covers(evolved['d21_evolved_upper'], Q(1, 7))
    traces = density_current_from_nuclear(
        evolved['d11_evolved_upper'], evolved['d21_evolved_upper'], Q(2, 5))
    stress = local_stress_from_traces(
        traces['density_trace_upper'], traces['current_trace_upper'],
        mass=1, absolute_angular=2, axial_lower=Q(4, 5), radius_lower=Q(7, 5),
        multiplicity=3)
    with ctx.workprec(200):
        density = restored_upper(traces['density_trace_upper'])
        current = restored_upper(traces['current_trace_upper'])
        potential = (arb(1) + (arb(2) / arb(fmpq(7, 5))) ** 2).sqrt()
        assert restored_upper(stress['N_upper']).lower() >= (
            arb(3) * (potential * density + current / arb(fmpq(4, 5)))).lower()
        assert restored_upper(stress['beta_upper']).lower() >= (arb(3) * current).lower()


def test_missing_inputs_remain_open_without_band_action():
    with pytest.raises(ValueError, match='explicit'):
        propagate_nuclear_seminorms(0, None, 0, 0, 0, 0)
    with pytest.raises(ValueError, match='explicit'):
        bloch_period_constants(None)
    with pytest.raises(TypeError):
        conditional_bloch_trace_bound(period=1, band_action_d11=0, band_action_d21=0)
    open_record = conditional_bloch_trace_bound(
        period=Q(2, 5), defect_symbol_I=Q(1, 7), initial_d11=0, initial_d21=0,
        causal_embedding_established=None)
    assert open_record['N_upper'] is None
    assert open_record['d11_evolved_upper'] is None
    assert open_record['integrated_d11_from_I_upper'] is not None
    assert 'residual_d11' in open_record['missing_inputs']
    assert 'band_action_d11' not in open_record['missing_inputs']
    assert 'defect_symbol_I' not in open_record['missing_inputs']
    assert open_record['band_action_remainder_input'] is False
    assert open_record['P_generates_Gamma_sub'] is False
    assert open_record['defect_symbol_I_is_instantaneous'] is True
    assert open_record['residual_nuclear_is_time_integrated'] is True
    assert open_record['I_route_equals_integrated_residual'] is False
    assert open_record['physical_local_gate'] == 'OPEN'
    assert open_record['I_axial_derivative_order'] == I_AXIAL_DERIVATIVE_ORDER == 4
    assert open_record['maximum_mixed_geometry_order'] == MAXIMUM_MIXED_GEOMETRY_ORDER == 9
    assert open_record['spatial_only_jets_sufficient'] is False
    closed_zero = conditional_bloch_trace_bound(
        period=Q(2, 5), initial_d11=0, initial_d21=0, residual_d11=0, residual_d21=0,
        Bz_integral=0, Bzz_integral=0, mass=1, absolute_angular=0,
        axial_lower=1, radius_lower=1, multiplicity=1,
        causal_embedding_established=True)
    with ctx.workprec(120):
        assert restored_upper(closed_zero['N_upper']) == 0
    assert closed_zero['physical_local_gate'] == 'OPEN'
    I_route = integrated_nuclear_from_symbol_I(Q(1, 7), Q(2, 5))
    assert I_route['maximum_mixed_geometry_order'] == 9
    assert I_route['I_route_equals_integrated_residual'] is False


def test_directed_error_arithmetic_and_symbol_I_L_powers():
    with ctx.workprec(200):
        constants = bloch_period_constants(2)
        assert _covers(constants['C0_squared_upper'], 2)
        assert _covers(constants['BZ_measure_dq_over_2pi'], Q(1, 2))
        I_route = integrated_nuclear_from_symbol_I(2, 2, spin_dim=2)
        S = restored_upper(constants['C0_squared_upper']).sqrt()
        assert _covers(I_route['integrated_H_upper'], 2)
        assert restored_upper(I_route['integrated_d11_upper']).lower() >= (
            S * arb(2).sqrt() * arb(2)).lower()
        assert _covers(I_route['integrated_H2_upper'], 8)
    full = conditional_bloch_trace_bound(
        period=Q(5, 2), initial_d11=Q(1, 50), initial_d21=Q(1, 40),
        residual_d11=Q(1, 80), residual_d21=Q(1, 70),
        Bz_integral=Q(1, 1000), Bzz_integral=Q(1, 20),
        mass=Q(3, 2), absolute_angular=2, axial_lower=Q(4, 5),
        radius_lower=Q(7, 5), multiplicity=5,
        causal_embedding_established=False, defect_symbol_I=Q(2, 11))
    assert full['N_upper'] is not None
    assert 'd11_total_upper' not in full
    assert full['local_agreement_with_physical_I'] is False
    assert full['physical_local_gate'] == 'OPEN'
    assert full['covariance_reconstruction'].startswith('P_g - P_ref')
    with ctx.workprec(200):
        assert restored_upper(full['N_upper']) > 0
        assert restored_upper(full['beta_upper']) > 0
        assert restored_upper(full['integrated_d11_from_I_upper']) != restored_upper(full['d11_evolved_upper'])
