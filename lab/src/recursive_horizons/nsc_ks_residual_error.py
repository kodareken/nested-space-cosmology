"""Directed conditional residual-to-observable bounds for the KS equation.

Input numbers must already bound continuous norms. Samples or ODE tolerances
are not such bounds. This module supplies error propagation, not those inputs.
"""
from fractions import Fraction
from math import isqrt


def nonnegative(value, name):
    if value is None or isinstance(value, bool):
        raise ValueError('explicit nonnegative bound required: ' + name)
    try:
        number = Fraction(value)
    except (ValueError, OverflowError, TypeError) as error:
        raise ValueError('finite bound required: ' + name) from error
    if number < 0:
        raise ValueError('nonnegative bound required: ' + name)
    return number


def sqrt_upper(value, bits=100):
    """Exact rational upper enclosure of a nonnegative square root."""
    q = nonnegative(value, 'square root argument')
    if isinstance(bits, bool) or not isinstance(bits, int) or bits < 1:
        raise ValueError('positive integer enclosure precision required')
    target = q.numerator << (2 * bits)
    root = isqrt(target // q.denominator)
    if root * root * q.denominator < target:
        root += 1
    return Fraction(root, 1 << bits)


def _triple(values, name):
    if values is None or len(values) != 3:
        raise ValueError('three derivative-order bounds required: ' + name)
    return tuple(nonnegative(v, name) for v in values)


def propagate_h2(initial, residual_integrals, Bz_integral, Bzz_integral):
    """Bound endpoint/sup error norms for orders 0,1,2 on one oriented slab.

    B_j integrals use absolute physical d rho; reversing the evolution does
    not make them negative. The norms include all weighted spinor columns.
    """
    e0, e1, e2 = _triple(initial, 'initial H2 errors')
    R0, R1, R2 = _triple(residual_integrals, 'continuous residual integrals')
    K1, K2 = nonnegative(Bz_integral, 'B_z integral'), nonnegative(Bzz_integral, 'B_zz integral')
    f0 = e0 + R0
    f1 = e1 + R1 + K1 * f0
    f2 = e2 + R2 + 2 * K1 * (e1 + R1) + (K1 * K1 + K2) * f0
    return (f0, f1, f2)


def pointwise_field_bounds(h2_errors, length, maximum_absolute_energy):
    """Periodic Sobolev bound, then exact source-carrier reconstruction."""
    e0, e1, e2 = _triple(h2_errors, 'H2 errors')
    L = nonnegative(length, 'numerical period')
    if L == 0:
        raise ValueError('positive numerical period required')
    E = nonnegative(maximum_absolute_energy, 'maximum source energy')
    mean, primitive = sqrt_upper(1 / L), sqrt_upper(L / 12)
    field = mean * e0 + primitive * e1
    # The derivative of a periodic envelope has zero mean.
    axial_envelope = primitive * e2
    return {'F': field, 'F_z': axial_envelope + E * field,
            'X_z': axial_envelope, 'mean_constant': mean, 'primitive_constant': primitive}


def pure_radius_commutator_integrals(w_norms, U_norms, normal_support,
                                     radius_lower, absolute_angular):
    """Upper bounds on integral ||B_z|| and ||B_zz|| for the owned radius family.

    w_norms/U_norms bound global sup derivatives 0,1,2. The change of variable
    |d rho|=a |dT| cancels the generator's 1/a. Only |chi|<=1 is used.
    """
    _, W1, W2 = _triple(w_norms, 'w derivative bounds')
    _, U1, U2 = _triple(U_norms, 'U derivative bounds')
    s = nonnegative(normal_support, 'normal support')
    r = nonnegative(radius_lower, 'radius lower bound')
    ell = nonnegative(absolute_angular, 'absolute angular label')
    if s == 0 or r == 0:
        raise ValueError('positive support and radius lower bound required')
    first = s**2 * W1 / 2 + s**4 * U1 / 24
    second = s**2 * W2 / 2 + s**4 * U2 / 24
    quadratic = s**3 * W1**2 / 3 + s**5 * W1 * U1 / 15 + s**7 * U1**2 / 252
    return ell * first / r**2, ell * (second / r**2 + 2 * quadratic / r**3)


def matter_error_bounds(field_norm, axial_norm, field_error, axial_error, *,
                        mass, absolute_angular, axial_lower, radius_lower,
                        multiplicity, source_operator_error=0):
    """N,beta error from weighted Frobenius field norms and coherent CAR source.

    The exact source covariance has operator norm<=1. An optional bound on
    C_exact-C_numeric is retained. Quadrature and contraction roundoff errors
    are separate; they are not inferred from this inequality.
    """
    f, z, e, ez, m, ell, a, r, mu, ec = [nonnegative(v, n) for v, n in (
        (field_norm, 'F norm'), (axial_norm, 'F_z norm'), (field_error, 'F error'),
        (axial_error, 'F_z error'), (mass, 'mass'), (absolute_angular, 'absolute angular'),
        (axial_lower, 'axial lower bound'), (radius_lower, 'radius lower bound'),
        (multiplicity, 'signed-family multiplicity'), (source_operator_error, 'source operator error'))]
    if min(a, r, mu) <= 0:
        raise ValueError('positive geometry bounds and multiplicity required')
    density = 2 * f * e + e**2 + f**2 * ec
    current = z * e + f * ez + e * ez + f * z * ec
    return {'N': mu * ((m + ell / r) * density + current / a),
            'beta': mu * current, 'density_trace_bound': density,
            'momentum_trace_bound': current,
            'scope': 'conditional on supplied continuous field/source norm bounds; physical gate not certified'}
