"""Conditional Chebyshev energy-interpolation remainder for the KS column propagator.

The universal 2x2 operator basis W starts at I and is interpolated in real
energy in the characteristic row-l2/max_z norm. This module does not apply
one W across mixed source-energy labels; per-label reconstruction is owned
elsewhere. Physical Frobenius field error uses an explicit sqrt(2) conversion
from the two spin rows. Not a source replacement or a UV/physical certificate.
"""
from fractions import Fraction

import numpy as np
from flint import arb, ctx, fmpq, fmpz

from .nsc_ks_ball_trajectory import complex_ball, exact_upper, restored_upper


NODE_FIELD_SOLVE_ERROR = (
    'MISSING: separate budget; the interpolation remainder assumes exact W, '
    'W_z, Y and Y_z at the Chebyshev roots')
ROUNDED_NODE_WEIGHT_ERROR = (
    'MISSING: separate budget; floating Chebyshev nodes and barycentric '
    'weights are not enclosed here')
PHYSICAL_LOCAL_GATE = (
    'OPEN: I=S(1)+[.12,.18] remains the owned interval; interpolation of W '
    'is not a local constraint or source certificate')
UV_CERTIFICATE = (
    'OPEN: finite real energy only; no UV tail, continuum source or '
    'physical PASS')
ROW_NORM = (
    'max_{z, spin row} Euclidean l2 of the two-entry row; not the induced '
    '2x2 spectral norm (the two spin rows travel on opposite characteristics)')


def _upper(value, name):
    if value is None or isinstance(value, bool):
        raise ValueError('explicit nonnegative bound required: ' + name)
    if isinstance(value, dict) and 'mantissa' in value and 'exponent' in value:
        value = restored_upper(value)
    if isinstance(value, arb):
        result = value
    else:
        try:
            q = Fraction(value)
            result = arb(fmpq(q.numerator, q.denominator))
        except (ValueError, TypeError, OverflowError) as error:
            raise ValueError('finite bound required: ' + name) from error
    if not result.is_finite() or not result >= 0:
        raise ValueError('finite nonnegative bound required: ' + name)
    return result.upper()


def _degree(value, name='degree'):
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError('nonnegative integer interpolation degree required: ' + name)
    return int(value)


def _factorial(n):
    acc = fmpz(1)
    for k in range(2, n + 1):
        acc *= k
    return arb(acc)


def row_l2_max_norm(matrix):
    """Characteristic comparison: max over z and spin of the row Euclidean l2.

    For a 2x2 block, each row is the two-entry C^2 vector of that spin. This
    is not numpy.linalg.norm(matrix, 2). Both rows (1, 0) give row-norm 1 and
    spectral norm sqrt(2).
    """
    value = np.asarray(matrix, complex)
    if value.shape[-2:] != (2, 2):
        raise ValueError('2x2 spin blocks required')
    if not np.isfinite(value.real).all() or not np.isfinite(value.imag).all():
        raise ValueError('finite 2x2 blocks required')
    row_squares = (np.abs(value) ** 2).sum(axis=-1)
    return float(np.sqrt(np.max(row_squares)))


def energy_derivative_majorants(offdiagonal_integral, characteristic_length,
                                Bz_integral, history_integral, history_z_integral,
                                n, *, bits=120):
    """Directed uppers for energy derivatives of the operator basis.

    The four majorants are in the characteristic row-l2/max_z norm: each spin
    row is a C^2 vector, the two rows propagate on opposite S3 characteristics,
    and the comparison is the max of those row Euclidean norms over z. That is
    not the induced 2x2 spectral norm. Orientation of decreasing rho is
    absorbed into the supplied absolute integrals. The Leibniz n! cancels the
    ordered simplex; mixed z/history insertions form one rectangle, not a
    factor two.
    """
    n = _degree(n, 'derivative order')
    with ctx.workprec(bits):
        K0 = _upper(offdiagonal_integral, 'offdiagonal integral')
        D = _upper(characteristic_length, 'characteristic length')
        K1 = _upper(Bz_integral, 'B_z integral')
        J0 = _upper(history_integral, 'history integral')
        J1 = _upper(history_z_integral, 'history z integral')
        growth = K0.exp()
        Dn = arb(1) if n == 0 else D ** n
        common = growth * Dn
        return {
            'dE_n_W_upper': exact_upper(common),
            'dE_n_Wz_upper': exact_upper(common * K1),
            'dE_n_Y_upper': exact_upper(common * J0),
            'dE_n_Yz_upper': exact_upper(common * (J1 + J0 * K1)),
            'growth_upper': exact_upper(growth),
            'D_power_upper': exact_upper(Dn),
            'order': n,
            'norm': ROW_NORM,
            'induced_2x2_spectral_norm': False,
            'simplex_factorial_cancelled': True,
            'mixed_insertions': 'rectangle',
            'rectangle_factor': 1,
        }


def chebyshev_root_interpolation_factor(degree, halfwidth, *, bits=120):
    """Monic Chebyshev-root factor h^{N+1}/(2^N (N+1)!) on [Ec-h, Ec+h].

    Nodes are the N+1 roots of T_{N+1}, not Chebyshev--Lobatto extrema. The
    bound uses |T_{N+1}|<=1 and leading coefficient 2^N, together with the
    Banach divided-difference estimate 1/(N+1)!, not a scalar mean-value
    point xi. The Banach space is matrices with the row-l2/max_z norm.
    """
    N = _degree(degree)
    with ctx.workprec(bits):
        h = _upper(halfwidth, 'interpolation halfwidth')
        return exact_upper(h ** (N + 1) / ((arb(2) ** N) * _factorial(N + 1)))


def conditional_energy_interpolation_bound(
        offdiagonal_integral, characteristic_length, Bz_integral,
        history_integral, history_z_integral, degree, halfwidth, *, bits=120):
    """Ideal interpolation remainders of W, W_z, Y and Y_z on a real energy interval.

    Remainders are row-l2/max_z uppers. Inputs must already be nonnegative
    certified integral uppers. They are not certified by being passed here.
    Node solves, rounded nodes/weights, A_up amplification and the carrier E
    in F_z remain separate. Mixed-energy column reconstruction is not
    performed here.
    """
    N = _degree(degree)
    with ctx.workprec(bits):
        K0 = _upper(offdiagonal_integral, 'offdiagonal integral')
        D = _upper(characteristic_length, 'characteristic length')
        K1 = _upper(Bz_integral, 'B_z integral')
        J0 = _upper(history_integral, 'history integral')
        J1 = _upper(history_z_integral, 'history z integral')
        h = _upper(halfwidth, 'interpolation halfwidth')
        growth = K0.exp()
        power = D ** (N + 1)
        factor = h ** (N + 1) / ((arb(2) ** N) * _factorial(N + 1))
        common = growth * power * factor
        majorants = energy_derivative_majorants(
            offdiagonal_integral, characteristic_length, Bz_integral,
            history_integral, history_z_integral, N + 1, bits=bits)
        return {
            'degree': N,
            'chebyshev_root_count': N + 1,
            'nodes': 'Chebyshev roots of T_{N+1} on [Ec-h, Ec+h]; not Lobatto extrema',
            'monic_factor': 'T_{N+1}/2^N',
            'banach_divided_differences': True,
            'scalar_mean_value_theorem': False,
            'simplex_factorial_cancelled': True,
            'mixed_insertions': 'rectangle',
            'rectangle_factor': 1,
            'norm': ROW_NORM,
            'induced_2x2_spectral_norm': False,
            'chebyshev_root_factor_upper': exact_upper(factor),
            'W_interpolation_remainder_upper': exact_upper(common),
            'Wz_interpolation_remainder_upper': exact_upper(common * K1),
            'Y_interpolation_remainder_upper': exact_upper(common * J0),
            'Yz_interpolation_remainder_upper': exact_upper(common * (J1 + J0 * K1)),
            'derivative_majorants_order_N_plus_1': majorants,
            'growth_upper': majorants['growth_upper'],
            'node_field_solve_error': NODE_FIELD_SOLVE_ERROR,
            'rounded_node_weight_error': ROUNDED_NODE_WEIGHT_ERROR,
            'source_A_up_amplification_included': False,
            'carrier_energy_in_F_z_included': False,
            'mixed_energy_reconstruction': False,
            'conditional_inputs_certified': False,
            'requires_certified_integral_uppers': True,
            'operator_basis': (
                'W(rho_up)=I interpolated only, in row-l2/max_z; per-label '
                'X=W(E) A_up(E) is not reconstructed here'),
            'physical_interval': 'I=S(1)+[.12,.18]',
            'physical_local_gate': PHYSICAL_LOCAL_GATE,
            'uv_certificate': UV_CERTIFICATE,
            'source_parameters_and_action': 'fixed; not modified',
            'scope': 'finite real energy interpolation remainder of the KS operator basis',
        }


def weighted_A_up_frobenius(A_up, column_weights, *, bits=120):
    """Weighted Frobenius norm of the original upstream columns, not of I."""
    initial = np.asarray(A_up, complex)
    weights = np.asarray(column_weights, float)
    if initial.ndim != 2 or initial.shape[0] != 2 or initial.shape[1] != len(weights):
        raise ValueError('A_up must match the original weighted column count')
    if not np.isfinite(initial.real).all() or not np.isfinite(initial.imag).all():
        raise ValueError('finite original A_up required')
    if not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError('original positive source column weights required')
    with ctx.workprec(bits):
        total = arb(0)
        for spin in range(2):
            for column, weight in enumerate(weights):
                total += (complex_ball(initial[spin, column]) * arb(float(weight))).abs_upper() ** 2
        return total.sqrt().upper()


def field_remainders_from_basis(W_remainder, Wz_remainder, weighted_A_up_norm,
                                absolute_energy, *, Y_remainder=None,
                                Yz_remainder=None, bits=120):
    """Convert row-l2/max remainders into weighted Frobenius F and F_z remainders.

    F = e^{-i E z} W A_up and F_z = e^{-i E z}(W_z - i E W) A_up. The phase
    is unitary. Each spin row r of Delta W satisfies
    |(Delta W A)_row| <= ||r||_2 ||A_col||_2, so
    ||Delta F||_{F,w} <= sqrt(2) * (row-l2/max of Delta W) * ||A_up||_{F,w}.
    The induced 2x2 spectral norm is not used: the two rows need not share a
    z after two-speed transport. The same sqrt(2) converts F_z, Y and Y_z.
    Node-solve and rounding errors are not included.
    """
    with ctx.workprec(bits):
        dW = _upper(W_remainder, 'W remainder')
        dWz = _upper(Wz_remainder, 'W_z remainder')
        amp = _upper(weighted_A_up_norm, 'weighted A_up norm')
        energy = _upper(absolute_energy, 'absolute energy')
        spin_rows = arb(2).sqrt()
        conversion = spin_rows * amp
        field = dW * conversion
        axial = (dWz + energy * dW) * conversion
        report = {
            'weighted_F_remainder_upper': exact_upper(field),
            'weighted_F_z_remainder_upper': exact_upper(axial),
            'A_up_amplification_upper': exact_upper(amp),
            'carrier_energy_upper': exact_upper(energy),
            'frobenius_conversion_factor_upper': exact_upper(spin_rows),
            'source_error_factor': 'sqrt(2) * row_l2_max * weighted_Frobenius(A_up)',
            'norm': ROW_NORM,
            'induced_2x2_spectral_norm': False,
            'source_A_up_amplification_included': True,
            'carrier_energy_in_F_z_included': True,
            'node_field_solve_error': NODE_FIELD_SOLVE_ERROR,
            'rounded_node_weight_error': ROUNDED_NODE_WEIGHT_ERROR,
            'physical_local_gate': PHYSICAL_LOCAL_GATE,
            'uv_certificate': UV_CERTIFICATE,
        }
        if Y_remainder is None and Yz_remainder is None:
            return report
        if Y_remainder is None or Yz_remainder is None:
            raise ValueError('history remainders Y and Y_z must be supplied together')
        dY = _upper(Y_remainder, 'Y remainder')
        dYz = _upper(Yz_remainder, 'Y_z remainder')
        report['weighted_Y_remainder_upper'] = exact_upper(dY * conversion)
        report['weighted_Yz_remainder_upper'] = exact_upper((dYz + energy * dY) * conversion)
        return report


def interval_energy_abs_upper(center, halfwidth, *, bits=120):
    """Directed |E| <= |Ec| + h on the interpolation interval. Not a node list."""
    if center is None or isinstance(center, bool):
        raise ValueError('explicit finite interpolation center required')
    with ctx.workprec(bits):
        h = _upper(halfwidth, 'interpolation halfwidth')
        try:
            mid = arb(fmpq(Fraction(center).numerator, Fraction(center).denominator))
        except (ValueError, TypeError, OverflowError) as error:
            raise ValueError('finite interpolation center required') from error
        if not mid.is_finite():
            raise ValueError('finite interpolation center required')
        return (mid.abs_upper() + h).upper()
