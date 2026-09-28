"""Conditional Bloch/Sobolev covariance-to-local-stress error propagator.

This is a finite-history UV norm helper for the periodically padded coefficient
problem. It does not certify a physical gate, existence, or a new source.
Nuclear bounds use Hilbert--Schmidt Holder, not a symbol L1 estimate. The
spinor HS of Lambda^{-1} supplies an explicit sqrt(spin_dim) factor; the
scalar-mode formula H S is false already for the 2x2 identity.

Reconstruction uses the owned Weyl Bloch matrix, fiber inner product
L^{-1} int_cell, and kernel measure dq/(2pi) sum_n with no extra 1/L.
Period L is kept as a strictly positive enclosing interval in every formula
that contains L or 1/L; denominators are never replaced by their uppers.
P is only an auxiliary add/subtract for the same C_g. Band-action remainder
is not a covariance nuclear input.
"""
from fractions import Fraction

import numpy as np
from flint import arb, ctx, fmpq

from .nsc_ks_ball_trajectory import exact_upper, restored_upper


I_AXIAL_DERIVATIVE_ORDER = 4
# P4 has four physical derivatives; time derivative of P adds one; z^4 of R
# adds four; the kinetic first-star term adds one spatial derivative.
MAXIMUM_MIXED_GEOMETRY_ORDER = 9
DEFAULT_SPIN_DIM = 2


def _bits(bits):
    if isinstance(bits, bool) or not isinstance(bits, int) or bits < 1:
        raise ValueError('positive integer enclosure precision required')
    return bits


def _spin_dim(spin_dim):
    if isinstance(spin_dim, bool) or not isinstance(spin_dim, int) or spin_dim < 1:
        raise ValueError('positive integer spin dimension required')
    return spin_dim


def _is_dyadic_upper(value):
    return isinstance(value, dict) and 'mantissa' in value and 'exponent' in value


def _enclose(value, name, *, bits):
    if value is None or isinstance(value, bool):
        raise ValueError('explicit finite bound required: ' + name)
    if _is_dyadic_upper(value):
        result = restored_upper(value)
    elif isinstance(value, arb):
        result = value
    else:
        try:
            q = Fraction(value)
            result = arb(fmpq(q.numerator, q.denominator))
        except (ValueError, TypeError, OverflowError) as error:
            raise ValueError('finite bound required: ' + name) from error
    if not result.is_finite():
        raise ValueError('finite bound required: ' + name)
    return result


def _norm_upper(value, name, *, bits):
    """Directed upper of a nonnegative seminorm. Accepts exact_upper records."""
    if _is_dyadic_upper(value) or not isinstance(value, dict):
        result = _enclose(value, name, bits=bits)
    else:
        raise ValueError('finite bound required: ' + name)
    if not result >= 0:
        raise ValueError('nonnegative bound required: ' + name)
    return result.upper()


def _positive_interval(value, name, *, bits):
    """Strictly positive enclosing interval; lower endpoint is not discarded.

    Collapsed dyadic-upper records are rejected: they have no certified lower
    bound, so 1/value would underbound if the upper were used as a denominator.
    """
    if _is_dyadic_upper(value):
        raise ValueError('enclosing positive interval required, not a collapsed upper: ' + name)
    result = _enclose(value, name, bits=bits)
    if not result > 0:
        raise ValueError('strictly positive enclosing interval required: ' + name)
    return result


def _restore(item):
    return restored_upper(item)


def _distinct_integer_modes(n_values):
    raw = np.atleast_1d(np.asarray(n_values, dtype=object))
    if raw.ndim != 1 or raw.size < 1:
        raise ValueError('one-dimensional Bloch mode indices required')
    modes = []
    for value in raw.tolist():
        if isinstance(value, (bool, np.bool_)):
            raise ValueError('integer Bloch mode indices required')
        if isinstance(value, (int, np.integer)):
            modes.append(int(value))
            continue
        try:
            q = Fraction(value)
        except (ValueError, TypeError, OverflowError, ZeroDivisionError) as error:
            raise ValueError('integer Bloch mode indices required') from error
        if q.denominator != 1:
            raise ValueError('integer Bloch mode indices required')
        modes.append(int(q))
    if len(modes) != len(set(modes)):
        raise ValueError('distinct Bloch mode indices required')
    return np.asarray(modes, dtype=np.int64)


def bloch_wave_numbers(period, q, n_values):
    """k_n = q + 2 pi n / L on one fiber of distinct integer modes."""
    length = float(period)
    if not np.isfinite(length) or length <= 0:
        raise ValueError('positive numerical period required')
    n_values = _distinct_integer_modes(n_values)
    q = float(q)
    if not np.isfinite(q):
        raise ValueError('finite Bloch momentum required')
    return q + 2 * np.pi * n_values / length


def fiber_evaluation_weights(k):
    """Actual C0^2 = S(q)^2 and C1^2 on one fiber, not the uniform majorants."""
    k = np.asarray(k, float)
    if k.ndim != 1 or not len(k) or not np.isfinite(k).all():
        raise ValueError('finite Bloch wave numbers required')
    lam2 = 1 + k * k
    return float(np.sum(1 / lam2)), float(np.sum(k * k / (lam2 * lam2)))


def weyl_bloch_matrix(a_hat, period, q, n_values):
    """A_mn(q) = a_hat_{m-n}((k_m+k_n)/2) in the owned Weyl convention.

    a_hat(l, k) returns a (spin, spin) array for the cell coefficient
    a_hat_l = L^{-1} int_cell e^{-i omega_l z} a(z,k) dz. The kernel is
    reconstructed with dq/(2pi) sum_n and no extra 1/L.
    """
    n_values = _distinct_integer_modes(n_values)
    k = bloch_wave_numbers(period, q, n_values)
    n_count = len(n_values)
    sample = np.asarray(a_hat(0, k[0]), complex)
    if sample.ndim != 2 or sample.shape[0] != sample.shape[1]:
        raise ValueError('square spin blocks required')
    spin = sample.shape[0]
    blocks = np.zeros((n_count, n_count, spin, spin), complex)
    for i, m in enumerate(n_values):
        for j, n in enumerate(n_values):
            blocks[i, j] = np.asarray(a_hat(int(m - n), 0.5 * (k[i] + k[j])), complex)
            if blocks[i, j].shape != (spin, spin) or not np.isfinite(blocks[i, j]).all():
                raise ValueError('finite square spin Fourier blocks required')
    return blocks, k


def _weighted_big(blocks, k, left_power, right_power):
    blocks = np.asarray(blocks, complex)
    k = np.asarray(k, float)
    if blocks.ndim != 4 or blocks.shape[0] != blocks.shape[1] or blocks.shape[2] != blocks.shape[3]:
        raise ValueError('square Bloch blocks of square spin matrices required')
    if k.shape != (blocks.shape[0],) or not np.isfinite(k).all() or not np.isfinite(blocks).all():
        raise ValueError('matching finite Bloch wave numbers required')
    weight = np.sqrt(1 + k * k)
    scaled = weight[:, None, None, None]**left_power * blocks * weight[None, :, None, None]**right_power
    n_count, spin = blocks.shape[0], blocks.shape[2]
    return scaled.transpose(0, 2, 1, 3).reshape(n_count * spin, n_count * spin)


def operator_hs(blocks, k, left_power, right_power):
    """Hilbert--Schmidt norm of Lambda^{left} A Lambda^{right} on the spinor fiber."""
    return float(np.linalg.norm(_weighted_big(blocks, k, left_power, right_power), 'fro'))


def operator_nuclear(blocks, k, left_power, right_power):
    """Nuclear (Schatten-1) norm of Lambda^{left} A Lambda^{right}."""
    return float(np.sum(np.linalg.svd(_weighted_big(blocks, k, left_power, right_power), compute_uv=False)))


def fiber_H(blocks, k):
    """H(q) = ||Lambda^2 A Lambda^2||_HS, including spin Frobenius."""
    return operator_hs(blocks, k, 2, 2)


def fiber_density(blocks, k, z):
    """Single-fiber 2x2 density sum_mn e^{i(k_m-k_n)z} A_mn, before dq/(2pi)."""
    k = np.asarray(k, float)
    blocks = np.asarray(blocks, complex)
    phases = np.exp(1j * k * float(z))
    return np.einsum('m,mnab,n->ab', phases, blocks, phases.conj())


def fiber_current(blocks, k, z):
    """Hermitian current (k_m+k_n)/2 contraction, matching (F_z C F^dag-h.c.)/(2i)."""
    k = np.asarray(k, float)
    blocks = np.asarray(blocks, complex)
    phases = np.exp(1j * k * float(z))
    midpoint = 0.5 * (k[:, None] + k[None, :])
    return np.einsum('m,mn,mnab,n->ab', phases, midpoint, blocks, phases.conj())


def matrix_nuclear(matrix):
    matrix = np.asarray(matrix, complex)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or not np.isfinite(matrix).all():
        raise ValueError('finite square matrix required')
    return float(np.sum(np.linalg.svd(matrix, compute_uv=False)))


def bloch_period_constants(period, *, bits=120):
    """Uniform computational-period constants, not real-line or large-L values.

    S(q)^2 = C0(q)^2 = sum_n (1+k_n^2)^{-1} <= 1 + L^2/4,
    C1(q)^2 = sum_n k_n^2/(1+k_n^2)^2 <= 1/4 + L^2/4,
    using |k_n| >= 2 pi (|n|-1/2)/L for n != 0 and the n=0 bounds 1 and 1/4.
    L is an enclosing positive interval; 1/L is the reciprocal of that same
    interval, not 1/L.upper().
    """
    bits = _bits(bits)
    with ctx.workprec(bits):
        length = _positive_interval(period, 'numerical period', bits=bits)
        tail = (length * length) / 4
        c0sq = arb(1) + tail
        c1sq = arb(1) / 4 + tail
        return {
            'period': exact_upper(length),
            'C0_squared_upper': exact_upper(c0sq),
            'C1_squared_upper': exact_upper(c1sq),
            'S_squared_upper': exact_upper(c0sq),
            'BZ_measure_dq_over_2pi': exact_upper(arb(1) / length),
            'real_line_C0_squared': '1/2, not used',
            'real_line_C1_squared': '1/4, not used',
            'large_L_C0_squared_asymptote': 'L/2, not an upper bound',
            'large_L_asymptotics_used': False,
        }


def shifted_weight_product_upper(k, s, *, bits=120):
    """Directed 8(<k>^8 + s^8) upper of <k+s>^4 <k-s>^4.

    (1+(k+s)^2)(1+(k-s)^2) <= (1+k^2+s^2)^2 and (A+B)^4 <= 8(A^4+B^4)
    for A=1+k^2, B=s^2, with equality in the second factor at A=B.
    """
    bits = _bits(bits)
    with ctx.workprec(bits):
        k = _enclose(k, 'k', bits=bits)
        s = _enclose(s, 'omega_l/2', bits=bits)
        eight = arb(8) * ((1 + k * k)**4 + s**8)
        return exact_upper(eight)


def hs_holder_nuclear_from_H(H, S, *, spin_dim=DEFAULT_SPIN_DIM, bits=120):
    """||Lambda A Lambda||_1 and ||Lambda^2 A Lambda||_1 <= H * sqrt(spin_dim) * S.

    Schatten Holder ||XY||_1 <= ||X||_HS ||Y||_HS on the spinor fiber, with
    ||Lambda A Lambda^2||_HS <= H because ||Lambda^{-1}||_op <= 1. The HS
    of Lambda^{-1} tensor I_spin is sqrt(spin_dim) S, not the scalar-mode S.
    """
    bits = _bits(bits)
    spin_dim = _spin_dim(spin_dim)
    with ctx.workprec(bits):
        H = _norm_upper(H, 'H(q)', bits=bits)
        S = _norm_upper(S, 'S(q)', bits=bits)
        bound = H * arb(spin_dim).sqrt() * S
        return {'d11_upper': exact_upper(bound), 'd21_upper': exact_upper(bound),
                'spin_hs_factor': exact_upper(arb(spin_dim).sqrt()),
                'symbol_L1_used_as_nuclear_bound': False}


def integrated_H_from_symbol_I(symbol_I, period, *, bits=120):
    """int_BZ H(q) dq/(2pi) <= sqrt(8 I[a]) / L from Parseval and Cauchy--Schwarz.

    I[a] is an instantaneous symbol seminorm. int_BZ H^2 dq/(2pi) <= (8/L) I[a]
    uses the same enclosing L as 1/L. This is not a time-integrated residual.
    """
    bits = _bits(bits)
    with ctx.workprec(bits):
        length = _positive_interval(period, 'numerical period', bits=bits)
        I_a = _norm_upper(symbol_I, 'defect symbol I[a]', bits=bits)
        eight_I = arb(8) * I_a
        return {
            'integrated_H2_upper': exact_upper(eight_I / length),
            'integrated_H_upper': exact_upper(eight_I.sqrt() / length),
            'parseval_eight_over_L': True,
            'defect_symbol_I_is_instantaneous': True,
        }


def integrated_nuclear_from_symbol_I(symbol_I, period, *, spin_dim=DEFAULT_SPIN_DIM, bits=120):
    """BZ-integrated d11,d21 from instantaneous I[a], including the spinor HS factor."""
    bits = _bits(bits)
    spin_dim = _spin_dim(spin_dim)
    with ctx.workprec(bits):
        length = _positive_interval(period, 'numerical period', bits=bits)
        I_a = _norm_upper(symbol_I, 'defect symbol I[a]', bits=bits)
        c0sq = arb(1) + (length * length) / 4
        eight_I = arb(8) * I_a
        integrated_H = eight_I.sqrt() / length
        bound = c0sq.sqrt() * arb(spin_dim).sqrt() * integrated_H
        return {
            'integrated_H2_upper': exact_upper(eight_I / length),
            'integrated_H_upper': exact_upper(integrated_H),
            'integrated_d11_upper': exact_upper(bound),
            'integrated_d21_upper': exact_upper(bound),
            'spin_hs_factor': exact_upper(arb(spin_dim).sqrt()),
            'parseval_eight_over_L': True,
            'I_axial_derivative_order': I_AXIAL_DERIVATIVE_ORDER,
            'maximum_mixed_geometry_order': MAXIMUM_MIXED_GEOMETRY_ORDER,
            'spatial_only_jets_sufficient': False,
            'formal_P4_supplies_I': False,
            'symbol_L1_used_as_nuclear_bound': False,
            'defect_symbol_I_is_instantaneous': True,
            'I_route_equals_integrated_residual': False,
        }


def h2_weighted_comparison_matrix(K1, K2, *, bits=120):
    """Triangle matrix for the weighted H2 vector (||u||, sqrt(2)||u_z||, ||u_zz||)."""
    bits = _bits(bits)
    with ctx.workprec(bits):
        K1 = _norm_upper(K1, 'B_z integral', bits=bits)
        K2 = _norm_upper(K2, 'B_zz integral', bits=bits)
        root2 = arb(2).sqrt()
        return (
            (arb(1), arb(0), arb(0)),
            (root2 * K1, arb(1), arb(0)),
            (K2 + K1 * K1, root2 * K1, arb(1)),
        )


def duhamel_sobolev_multipliers(K1, K2, *, bits=120):
    """KS L2-unitary Sobolev multipliers, reusing the residual-error commutators.

    M1 <= 1+K1 for H1. The H2 comparison perturbation has Frobenius norm
    sqrt(4 K1^2 + (K2+K1^2)^2), so M2 <= 1 + that quantity. Global K is a
    conservative bound on every subinterval.
    """
    bits = _bits(bits)
    with ctx.workprec(bits):
        K1 = _norm_upper(K1, 'B_z integral', bits=bits)
        K2 = _norm_upper(K2, 'B_zz integral', bits=bits)
        M1 = arb(1) + K1
        frobenius = (arb(4) * K1 * K1 + (K2 + K1 * K1)**2).sqrt()
        M2 = arb(1) + frobenius
        return {
            'M1_upper': exact_upper(M1),
            'M2_upper': exact_upper(M2),
            'K1_upper': exact_upper(K1),
            'K2_upper': exact_upper(K2),
            'comparison_frobenius_upper': exact_upper(frobenius),
            'same_KS_U_L2_unitary': True,
        }


def propagate_nuclear_seminorms(initial_d11, residual_d11, initial_d21, residual_d21,
                                K1, K2, *, bits=120):
    """d11_final <= M1^2 (initial11 + integral residual11), d21_final <= M2 M1 (...).

    residual_d11/d21 are time-integrated nuclear seminorms of R. They are not
    the instantaneous I[a] symbol route.
    """
    bits = _bits(bits)
    multipliers = duhamel_sobolev_multipliers(K1, K2, bits=bits)
    with ctx.workprec(bits):
        i11 = _norm_upper(initial_d11, 'initial D_up d11', bits=bits)
        r11 = _norm_upper(residual_d11, 'integrated R_op d11', bits=bits)
        i21 = _norm_upper(initial_d21, 'initial D_up d21', bits=bits)
        r21 = _norm_upper(residual_d21, 'integrated R_op d21', bits=bits)
        M1 = _restore(multipliers['M1_upper'])
        M2 = _restore(multipliers['M2_upper'])
        d11 = (M1 * M1) * (i11 + r11)
        d21 = (M2 * M1) * (i21 + r21)
        return {
            **multipliers,
            'd11_evolved_upper': exact_upper(d11),
            'd21_evolved_upper': exact_upper(d21),
            'initial_D_up_retained': True,
            'integrated_R_op_retained': True,
            'I_route_equals_integrated_residual': False,
        }


def density_current_from_nuclear(d11, d21, period, *, bits=120):
    """Pointwise spinor nuclear traces: density <= C0^2 d11, current <= C0 C1 d21.

    For Hermitian D the symmetric-derivative halves cancel against equal
    adjoint nuclear norms, so no extra factor two is inserted. C0,C1 and 1/L
    use the same enclosing period.
    """
    bits = _bits(bits)
    with ctx.workprec(bits):
        length = _positive_interval(period, 'numerical period', bits=bits)
        d11 = _norm_upper(d11, 'd11', bits=bits)
        d21 = _norm_upper(d21, 'd21', bits=bits)
        tail = (length * length) / 4
        c0sq = arb(1) + tail
        c1sq = arb(1) / 4 + tail
        density = c0sq * d11
        current = c0sq.sqrt() * c1sq.sqrt() * d21
        constants = bloch_period_constants(period, bits=bits)
        return {
            **constants,
            'density_trace_upper': exact_upper(density),
            'current_trace_upper': exact_upper(current),
            'hermitian_symmetric_current': True,
        }


def local_stress_from_traces(density, current, *, mass, absolute_angular,
                             axial_lower, radius_lower, multiplicity, bits=120):
    """N <= mu (sqrt(m^2+ell^2/r^2) density + current/a), beta <= mu current.

    axial_lower and radius_lower stay strictly positive enclosing intervals.
    Division uses those intervals, so the exported uppers track the lower
    endpoints. Norms and mass/angular/multiplicity collapse upward.
    """
    bits = _bits(bits)
    with ctx.workprec(bits):
        density = _norm_upper(density, 'density trace', bits=bits)
        current = _norm_upper(current, 'current trace', bits=bits)
        m = _norm_upper(mass, 'mass', bits=bits)
        ell = _norm_upper(absolute_angular, 'absolute angular', bits=bits)
        a = _positive_interval(axial_lower, 'axial lower bound', bits=bits)
        r = _positive_interval(radius_lower, 'radius lower bound', bits=bits)
        mu = _norm_upper(multiplicity, 'signed-family multiplicity', bits=bits)
        if not mu > 0:
            raise ValueError('positive geometry bounds and multiplicity required')
        potential = (m * m + (ell / r) ** 2).sqrt()
        N = mu * (potential * density + current / a)
        beta = mu * current
        return {
            'N_upper': exact_upper(N),
            'beta_upper': exact_upper(beta),
            'potential_op_norm': exact_upper(potential),
        }


def _missing_names(values):
    return tuple(name for name, value in values.items() if value is None)


def conditional_bloch_trace_bound(*, period, defect_symbol_I=None,
                                  initial_d11=None, initial_d21=None,
                                  residual_d11=None, residual_d21=None,
                                  Bz_integral=None, Bzz_integral=None,
                                  mass=None, absolute_angular=None,
                                  axial_lower=None, radius_lower=None,
                                  multiplicity=None,
                                  causal_embedding_established=None,
                                  spin_dim=DEFAULT_SPIN_DIM, bits=120):
    """Propagate certified covariance seminorms to local N,beta, or leave missing open.

    None is never replaced by zero. Initial D_up and time-integrated true R_op
    are the covariance inputs. Instantaneous I[a] is optional and is not the
    residual. P is an auxiliary add/subtract on the same C_g; this helper does
    not take a band-action remainder or assert that P generates Gamma_sub.
    physical_local_gate is always OPEN.
    """
    bits = _bits(bits)
    spin_dim = _spin_dim(spin_dim)
    if causal_embedding_established not in (None, True, False):
        raise ValueError('causal embedding assumption must be True, False, or None')
    constants = bloch_period_constants(period, bits=bits)
    record = {
        **constants,
        'spin_dim': spin_dim,
        'defect_symbol_I_upper': None,
        'integrated_H_upper': None,
        'integrated_d11_from_I_upper': None,
        'integrated_d21_from_I_upper': None,
        'd11_evolved_upper': None,
        'd21_evolved_upper': None,
        'density_trace_upper': None,
        'current_trace_upper': None,
        'N_upper': None,
        'beta_upper': None,
        'M1_upper': None,
        'M2_upper': None,
        'causal_embedding_established': causal_embedding_established,
        'local_agreement_with_physical_I': causal_embedding_established is True,
        'physical_periodic_boundary_condition': False,
        'formal_P4_supplies_I': False,
        'I_axial_derivative_order': I_AXIAL_DERIVATIVE_ORDER,
        'maximum_mixed_geometry_order': MAXIMUM_MIXED_GEOMETRY_ORDER,
        'spatial_only_jets_sufficient': False,
        'symbol_L1_used_as_nuclear_bound': False,
        'large_L_asymptotics_used': False,
        'band_action_remainder_input': False,
        'P_generates_Gamma_sub': False,
        'source_law': 'C_g = U_g C_up U_g^dagger with the same C_up',
        'residual_operator': 'R = P_dot - [G,P] of the actual subtracted Weyl operator',
        'duhamel_state': 'D = C - P',
        'covariance_reconstruction': 'P_g - P_ref plus D_g - D_ref of the same C_g',
        'nuclear_bound_via': 'Hilbert-Schmidt Holder on the spinor fiber, not symbol L1',
        'defect_symbol_I_is_instantaneous': True,
        'residual_nuclear_is_time_integrated': True,
        'I_route_equals_integrated_residual': False,
        'existence_or_nonexistence_claim': False,
        'physical_local_gate': 'OPEN',
        'scope': (
            'conditional computational-period Bloch/Sobolev covariance propagator; '
            'local agreement on I requires a separately established causal embedding; '
            'band-action generating-functional remainder is out of scope'
        ),
    }
    missing = {
        'initial_d11': initial_d11,
        'initial_d21': initial_d21,
        'residual_d11': residual_d11,
        'residual_d21': residual_d21,
        'Bz_integral': Bz_integral,
        'Bzz_integral': Bzz_integral,
        'mass': mass,
        'absolute_angular': absolute_angular,
        'axial_lower': axial_lower,
        'radius_lower': radius_lower,
        'multiplicity': multiplicity,
        'causal_embedding_established': causal_embedding_established,
    }
    record['missing_inputs'] = _missing_names(missing)

    if defect_symbol_I is not None:
        I_route = integrated_nuclear_from_symbol_I(
            defect_symbol_I, period, spin_dim=spin_dim, bits=bits)
        record['defect_symbol_I_upper'] = exact_upper(
            _norm_upper(defect_symbol_I, 'defect symbol I[a]', bits=bits))
        record['integrated_d11_from_I_upper'] = I_route['integrated_d11_upper']
        record['integrated_d21_from_I_upper'] = I_route['integrated_d21_upper']
        record['integrated_H_upper'] = I_route['integrated_H_upper']

    propagator_inputs = (initial_d11, residual_d11, initial_d21, residual_d21,
                         Bz_integral, Bzz_integral)
    if all(v is not None for v in propagator_inputs):
        evolved = propagate_nuclear_seminorms(
            initial_d11, residual_d11, initial_d21, residual_d21,
            Bz_integral, Bzz_integral, bits=bits)
        record['M1_upper'] = evolved['M1_upper']
        record['M2_upper'] = evolved['M2_upper']
        record['d11_evolved_upper'] = evolved['d11_evolved_upper']
        record['d21_evolved_upper'] = evolved['d21_evolved_upper']
        traces = density_current_from_nuclear(
            evolved['d11_evolved_upper'], evolved['d21_evolved_upper'],
            period, bits=bits)
        record['density_trace_upper'] = traces['density_trace_upper']
        record['current_trace_upper'] = traces['current_trace_upper']
        geometry = (mass, absolute_angular, axial_lower, radius_lower, multiplicity)
        if all(v is not None for v in geometry):
            stress = local_stress_from_traces(
                traces['density_trace_upper'], traces['current_trace_upper'],
                mass=mass, absolute_angular=absolute_angular,
                axial_lower=axial_lower, radius_lower=radius_lower,
                multiplicity=multiplicity, bits=bits)
            record['N_upper'] = stress['N_upper']
            record['beta_upper'] = stress['beta_upper']
    elif Bz_integral is not None and Bzz_integral is not None:
        multipliers = duhamel_sobolev_multipliers(Bz_integral, Bzz_integral, bits=bits)
        record['M1_upper'] = multipliers['M1_upper']
        record['M2_upper'] = multipliers['M2_upper']
    return record
