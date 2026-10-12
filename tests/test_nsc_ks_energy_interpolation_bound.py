"""Row-l2 interpolation remainders and sqrt(2) Frobenius conversion; no physics runs."""
from fractions import Fraction as Q

import numpy as np
import pytest
pytest.importorskip('flint')
from flint import arb, ctx, fmpq
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from recursive_horizons import nsc_ks_energy_interpolation_bound as bound_mod
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_energy_interpolation_bound import (
    chebyshev_root_interpolation_factor, conditional_energy_interpolation_bound,
    energy_derivative_majorants, field_remainders_from_basis,
    interval_energy_abs_upper, row_l2_max_norm, weighted_A_up_frobenius,
)


S1 = np.array([[0, 1], [1, 0]], dtype=complex)
S3 = np.diag([1.0, -1.0]).astype(complex)


def _chebyshev_root_nodes(count, center, halfwidth):
    theta = np.pi * (2 * np.arange(count) + 1) / (2 * count)
    return center + halfwidth * np.cos(theta)


def _interpolate_matrix(func, degree, center, halfwidth, energies):
    nodes = _chebyshev_root_nodes(degree + 1, center, halfwidth)
    weights = (-1.0) ** np.arange(degree + 1) * np.sin(
        np.pi * (2 * np.arange(degree + 1) + 1) / (2 * (degree + 1)))
    table = np.stack([func(float(node)) for node in nodes], axis=0)
    values = np.empty((len(energies), 2, 2), complex)
    for k, energy in enumerate(energies):
        delta = energy - nodes
        nearest = int(np.argmin(np.abs(delta)))
        if abs(delta[nearest]) <= 1e-13 * max(1.0, abs(energy)):
            values[k] = table[nearest]
            continue
        coeff = weights / delta
        values[k] = np.tensordot(coeff, table, axes=(0, 0)) / coeff.sum()
    return values


def test_leibniz_n_cancels_ordered_simplex_and_rectangle_has_no_factor_two():
    D, n = Q(3, 2), 4
    assert n * (D ** n) / n == D ** n
    K1, J0, J1 = Q(5, 7), Q(2, 3), Q(1, 11)
    with ctx.workprec(150):
        report = energy_derivative_majorants(0, D, K1, J0, J1, n)
        Dn = arb(fmpq(D.numerator, D.denominator)) ** n
        assert restored_upper(report['dE_n_W_upper']) == Dn.upper()
        assert restored_upper(report['dE_n_Wz_upper']) >= (
            arb(fmpq(K1.numerator, K1.denominator)) * Dn).upper()
        mixed = (arb(fmpq(J1.numerator, J1.denominator))
                 + arb(fmpq(J0.numerator, J0.denominator)) * arb(fmpq(K1.numerator, K1.denominator)))
        naive_two = mixed + arb(fmpq(J0.numerator, J0.denominator)) * arb(fmpq(K1.numerator, K1.denominator))
        assert restored_upper(report['dE_n_Yz_upper']) >= (mixed * Dn).upper()
        assert restored_upper(report['dE_n_Yz_upper']) < (naive_two * Dn).upper()
        assert report['rectangle_factor'] == 1
        assert report['simplex_factorial_cancelled']
        assert report['induced_2x2_spectral_norm'] is False
        assert restored_upper(report['dE_n_W_upper']) < (arb(24) * Dn).upper()


def test_row_l2_max_is_not_induced_spectral_norm():
    both_first = np.array([[1.0, 0.0], [1.0, 0.0]], complex)
    assert row_l2_max_norm(both_first) == pytest.approx(1.0)
    assert np.linalg.norm(both_first, 2) == pytest.approx(np.sqrt(2))
    # Opposite-characteristic snapshots at two z must not be stacked into one
    # spectral 2x2: each z has row-l2 1, while the stacked block has op-norm sqrt(2).
    transported = np.array([[[1.0, 0.0], [0.0, 0.0]],
                            [[0.0, 0.0], [1.0, 0.0]]], complex)
    assert row_l2_max_norm(transported) == pytest.approx(1.0)
    stacked = transported[0] + transported[1]
    assert row_l2_max_norm(stacked) == pytest.approx(1.0)
    assert np.linalg.norm(stacked, 2) == pytest.approx(np.sqrt(2))


def test_frobenius_field_error_needs_sqrt2_times_row_times_Aup():
    delta_W = np.array([[1.0, 0.0], [1.0, 0.0]], complex)
    A_up = np.array([[0.4, 0.3], [0.0, 0.0]], complex)
    weights = np.array([1.25, 0.8])
    X = delta_W @ (A_up * weights)
    actual = float(np.linalg.norm(X))
    row = row_l2_max_norm(delta_W)
    with ctx.workprec(150):
        amp = weighted_A_up_frobenius(A_up, weights)
        omitted = float(row * amp)
        converted = float(arb(2).sqrt() * amp)
        assert row == pytest.approx(1.0)
        assert actual > omitted
        assert actual <= converted + 1e-15
        fields = field_remainders_from_basis(row, 0, amp, 0)
        assert restored_upper(fields['weighted_F_remainder_upper']) >= arb(2).sqrt() * amp
        assert restored_upper(fields['frobenius_conversion_factor_upper']) >= arb(2).sqrt().upper()
        assert fields['induced_2x2_spectral_norm'] is False
        assert 'sqrt(2)' in fields['source_error_factor']


def test_carrier_and_tangent_frobenius_use_the_same_sqrt2():
    delta_W = np.array([[0.0, 1.0], [0.0, 1.0]], complex)
    delta_Wz = np.array([[0.5, 0.0], [0.5, 0.0]], complex)
    A_up = np.array([[0.0], [2.0]], complex)
    weights = np.array([0.75])
    energy = 3.0
    envelope = delta_Wz - 1j * energy * delta_W
    F = delta_W @ (A_up * weights[0])
    Fz = envelope @ (A_up * weights[0])
    Y = delta_W @ (A_up * weights[0])
    Yz = envelope @ (A_up * weights[0])
    with ctx.workprec(150):
        amp = weighted_A_up_frobenius(A_up, weights)
        fields = field_remainders_from_basis(
            row_l2_max_norm(delta_W), row_l2_max_norm(delta_Wz), amp, energy,
            Y_remainder=row_l2_max_norm(delta_W), Yz_remainder=row_l2_max_norm(delta_Wz))
        conversion = arb(2).sqrt() * amp
        assert float(np.linalg.norm(F)) <= float(restored_upper(fields['weighted_F_remainder_upper'])) + 1e-12
        assert float(np.linalg.norm(Fz)) <= float(restored_upper(fields['weighted_F_z_remainder_upper'])) + 1e-12
        assert float(np.linalg.norm(Y)) <= float(restored_upper(fields['weighted_Y_remainder_upper'])) + 1e-12
        assert float(np.linalg.norm(Yz)) <= float(restored_upper(fields['weighted_Yz_remainder_upper'])) + 1e-12
        omitted_z = (row_l2_max_norm(delta_Wz) + energy * row_l2_max_norm(delta_W)) * float(amp)
        assert float(np.linalg.norm(Fz)) > omitted_z
        assert restored_upper(fields['weighted_F_z_remainder_upper']) >= (
            (arb(row_l2_max_norm(delta_Wz)) + arb(energy) * arb(row_l2_max_norm(delta_W))) * conversion).upper()


def test_commuting_exponential_derivatives_are_sharp_in_row_norm():
    D, n = 2, 5
    W = expm(1j * 1.3 * D * S3)
    dW = (1j * D * S3) ** n @ W
    assert row_l2_max_norm(dW) == pytest.approx(D ** n, rel=0, abs=1e-12)
    with ctx.workprec(150):
        report = energy_derivative_majorants(0, D, 0, 0, 0, n)
        assert restored_upper(report['dE_n_W_upper']) == arb(D ** n)
        assert restored_upper(report['dE_n_Wz_upper']) == 0
        assert restored_upper(report['dE_n_Y_upper']) == 0


def test_commuting_exponential_chebyshev_root_interpolation_is_enclosed():
    D, center, h, N = 1.0, Q(3, 2), Q(1, 4), 5
    def W(E):
        return expm(1j * float(E) * D * S3)
    sample = np.linspace(float(center - h), float(center + h), 81)
    approx = _interpolate_matrix(W, N, float(center), float(h), sample)
    actual = max(row_l2_max_norm(approx[k] - W(E)) for k, E in enumerate(sample))
    with ctx.workprec(150):
        bound = conditional_energy_interpolation_bound(0, D, 0, 0, 0, N, h)
        upper = restored_upper(bound['W_interpolation_remainder_upper'])
        analytic = (arb(D) * arb(fmpq(h.numerator, h.denominator))) ** (N + 1) / (
            (arb(2) ** N) * arb(720))
        assert upper >= analytic.upper()
        assert actual <= float(upper)
        assert bound['scalar_mean_value_theorem'] is False
        assert bound['induced_2x2_spectral_norm'] is False
        assert bound['mixed_energy_reconstruction'] is False
        assert bound['chebyshev_root_count'] == N + 1


def test_noncommuting_2x2_propagate_interpolation_is_enclosed():
    alpha, beta, T, center, h, N = 0.4, 0.7, 1.0, 0.8, 0.2, 4
    K0, D = alpha * T, beta * T

    def W(E):
        return expm(1j * T * (alpha * S1 + E * beta * S3))

    sample = np.linspace(center - h, center + h, 51)
    approx = _interpolate_matrix(W, N, center, h, sample)
    actual = max(row_l2_max_norm(approx[k] - W(E)) for k, E in enumerate(sample))
    with ctx.workprec(150):
        bound = conditional_energy_interpolation_bound(K0, D, 0, 0, 0, N, h)
        assert actual <= float(restored_upper(bound['W_interpolation_remainder_upper']))


def test_time_dependent_noncommuting_ode_interpolation_is_enclosed():
    alpha, center, h, N = 0.25, 1.1, 0.15, 3
    D, K0 = 1.5, alpha

    def W(E):
        y0 = np.eye(2, dtype=complex).ravel()

        def rhs(t, y):
            L = 1j * (alpha * S1 + E * (1.0 + t) * S3)
            return (L @ y.reshape(2, 2)).ravel()

        sol = solve_ivp(rhs, (0.0, 1.0), y0, rtol=1e-11, atol=1e-13, dense_output=False)
        assert sol.success
        return sol.y[:, -1].reshape(2, 2)

    sample = np.linspace(center - h, center + h, 21)
    approx = _interpolate_matrix(W, N, center, h, sample)
    actual = max(row_l2_max_norm(approx[k] - W(E)) for k, E in enumerate(sample))
    with ctx.workprec(150):
        bound = conditional_energy_interpolation_bound(K0, D, 0, 0, 0, N, h)
        assert actual <= float(restored_upper(bound['W_interpolation_remainder_upper']))


def test_history_parameter_interpolation_uses_J0_not_identity_state():
    beta, T, kappa, center, h, N = 0.5, 1.0, 0.3, 0.5, 0.2, 3
    D, J0 = beta * T, abs(kappa) * T

    def Y(E):
        def W_at(alpha):
            return expm(1j * T * (alpha * S1 + E * beta * S3))
        delta = 1e-7
        return (W_at(kappa + delta) - W_at(kappa - delta)) / (2 * delta)

    sample = np.linspace(center - h, center + h, 31)
    approx = _interpolate_matrix(Y, N, center, h, sample)
    actual = max(row_l2_max_norm(approx[k] - Y(E)) for k, E in enumerate(sample))
    with ctx.workprec(150):
        bound = conditional_energy_interpolation_bound(abs(kappa) * T, D, 0, J0, 0, N, h)
        assert actual <= float(restored_upper(bound['Y_interpolation_remainder_upper']))
        assert restored_upper(bound['W_interpolation_remainder_upper']) > 0


def test_zero_integrals_and_missing_inputs():
    with ctx.workprec(120):
        zero = conditional_energy_interpolation_bound(0, 0, 0, 0, 0, 6, Q(1, 5))
        assert restored_upper(zero['W_interpolation_remainder_upper']) == 0
        assert restored_upper(zero['Yz_interpolation_remainder_upper']) == 0
        point = conditional_energy_interpolation_bound(Q(1, 3), 2, 1, 1, 1, 3, 0)
        assert restored_upper(point['W_interpolation_remainder_upper']) == 0
    with pytest.raises(ValueError, match='explicit'):
        energy_derivative_majorants(None, 1, 0, 0, 0, 1)
    with pytest.raises(ValueError, match='explicit'):
        conditional_energy_interpolation_bound(0, 1, 0, 0, None, 2, 0.1)
    with pytest.raises(ValueError, match='nonnegative integer'):
        conditional_energy_interpolation_bound(0, 1, 0, 0, 0, True, 0.1)
    with pytest.raises(ValueError, match='nonnegative'):
        energy_derivative_majorants(-1, 1, 0, 0, 0, 0)


def test_weighted_Aup_norm_is_not_the_identity_and_mixed_energy_is_not_reconstructed():
    source = FixedSourcePreparation.from_signed_blocks(
        np.array([0.7, -0.4]), np.array([1.1, 0.8]),
        np.array([[[0.45, 0.06 + 0.03j], [0.06 - 0.03j, 0.5]],
                  [[0.55, -0.04j], [0.04j, 0.4]]], complex))
    A_up = np.array([[0.6, 0.2, -0.1j, 0.3], [0.15j, 0.7, 0.4, -0.2]], complex)
    assert not hasattr(bound_mod, 'reconstruct_columns')
    with ctx.workprec(150):
        amp = weighted_A_up_frobenius(A_up, source.column_weights)
        identity_amp = weighted_A_up_frobenius(np.eye(2, dtype=complex), source.column_weights[:2])
        assert amp != identity_amp
        assert not np.allclose(source.covariance, np.eye(source.covariance.shape[0]))
        assert np.unique(source.energies).size > 1


def test_field_remainder_keeps_Aup_amplification_carrier_E_and_sqrt2():
    source = FixedSourcePreparation.from_signed_blocks(
        np.array([1.25]), np.array([1.2]), np.array([[[0.8, 0.1j], [-0.1j, 0.5]]], complex))
    A_up = np.array([[0.4, -0.2j], [0.3, 0.5]], complex)
    with ctx.workprec(150):
        amp = weighted_A_up_frobenius(A_up, source.column_weights)
        energy = interval_energy_abs_upper(Q(5, 4), Q(1, 5))
        bound = conditional_energy_interpolation_bound(Q(1, 10), Q(2, 5), Q(1, 7), Q(1, 8), Q(1, 9), 2, Q(1, 5))
        fields = field_remainders_from_basis(
            bound['W_interpolation_remainder_upper'],
            bound['Wz_interpolation_remainder_upper'],
            amp, energy,
            Y_remainder=bound['Y_interpolation_remainder_upper'],
            Yz_remainder=bound['Yz_interpolation_remainder_upper'])
        dW = restored_upper(bound['W_interpolation_remainder_upper'])
        dWz = restored_upper(bound['Wz_interpolation_remainder_upper'])
        dY = restored_upper(bound['Y_interpolation_remainder_upper'])
        dYz = restored_upper(bound['Yz_interpolation_remainder_upper'])
        conversion = arb(2).sqrt() * amp
        assert restored_upper(fields['weighted_F_remainder_upper']) >= (dW * conversion).upper()
        assert restored_upper(fields['weighted_F_z_remainder_upper']) >= ((dWz + energy * dW) * conversion).upper()
        assert restored_upper(fields['weighted_Y_remainder_upper']) >= (dY * conversion).upper()
        assert restored_upper(fields['weighted_Yz_remainder_upper']) >= ((dYz + energy * dY) * conversion).upper()
        assert restored_upper(fields['weighted_F_remainder_upper']) > (dW * amp).upper()
        assert fields['source_A_up_amplification_included']
        assert fields['carrier_energy_in_F_z_included']
        assert float(energy) >= 1.25 + 0.2


def test_scope_keeps_gate_open_and_node_errors_missing():
    before = ctx.prec
    report = conditional_energy_interpolation_bound(Q(1, 8), 1, Q(1, 9), Q(1, 4), Q(1, 5), 1, Q(1, 3))
    assert ctx.prec == before
    assert report['physical_local_gate'].startswith('OPEN')
    assert report['uv_certificate'].startswith('OPEN')
    assert report['node_field_solve_error'].startswith('MISSING')
    assert report['rounded_node_weight_error'].startswith('MISSING')
    assert report['conditional_inputs_certified'] is False
    assert report['mixed_energy_reconstruction'] is False
    assert report['induced_2x2_spectral_norm'] is False
    assert 'row' in report['norm']
    assert report['physical_interval'] == 'I=S(1)+[.12,.18]'
    with ctx.workprec(120):
        assert restored_upper(chebyshev_root_interpolation_factor(0, 1)) == arb(1)
        assert restored_upper(chebyshev_root_interpolation_factor(1, 1)) == arb(1) / 4
