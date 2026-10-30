"""Exact scalar polynomial calculus on the unit round three-sphere.

``X`` is a four-component ambient representation with ``R2 = 1``.  Its
component count is not a physical dimension statement.  Coefficients may
contain exact algebraic numbers or symbols independent of X.  The finite
algorithms are intended for the low-degree harmonic products in the atlas;
they do not choose a physical geometry, state, or truncation.

``mean`` uses normalized sphere measure.  ``inner_product`` is the real
bilinear pairing mean(f*g); multiply by 2*pi**2 for the unnormalized L2
pairing.  No conjugation of symbolic coefficients is implicit.
"""
from functools import lru_cache
from math import factorial, prod

import sympy as s

X = s.Matrix(s.symbols('x1:5', real=True))
R2 = X.dot(X)


def _poly(expr):
    """Use only the four spatial variables as polynomial generators."""
    expr = s.sympify(expr)
    names = {x.name for x in X}
    if any(z.name in names and z not in X for z in expr.free_symbols):
        raise ValueError('use the exported coordinate symbols X with their real assumptions')
    try:
        return s.Poly(s.expand(expr), *X, domain='EX')
    except s.PolynomialError as error:
        raise ValueError('an ambient polynomial in X is required') from error


def _monomial(powers):
    return s.prod(x**n for x, n in zip(X, powers))


def ambient_laplacian(expr):
    """Euclidean four-variable Laplacian of a polynomial extension."""
    polynomial = _poly(expr).as_expr()
    return s.expand(sum(s.diff(polynomial, x, 2) for x in X))


def laplacian(expr):
    """Unit-S3 scalar Laplacian, valid after restricting R2 to one.

    On each homogeneous term of degree n this is Delta_R4 - n(n+2).
    Returning a polynomial extension retains exact harmonic identities.
    """
    polynomial = _poly(expr)
    radial = sum(sum(powers)*(sum(powers)+2)*coefficient*_monomial(powers)
                 for powers, coefficient in polynomial.terms())
    return s.expand(ambient_laplacian(polynomial.as_expr()) - radial)


def gradient_dot(f, g):
    """Intrinsic gradient pairing on R2=1, as an ambient polynomial."""
    f = _poly(f).as_expr()
    g = _poly(g).as_expr()
    df = s.Matrix([s.diff(f, x) for x in X])
    dg = s.Matrix([s.diff(g, x) for x in X])
    return s.expand(df.dot(dg) - X.dot(df)*X.dot(dg))


def reduce_on_sphere(expr):
    """Deterministic normal form modulo R2-1, in the fixed X ordering.

    This is a quotient-ring representation, not an extra coordinate chart.
    The leading monomial x1**2 is eliminated by exact polynomial division.
    """
    polynomial = _poly(expr)
    relation = s.Poly(R2-1, *X, domain='EX')
    return s.expand(polynomial.rem(relation).as_expr())


def harmonic_projection_homogeneous(polynomial, degree):
    """Ambient-harmonic degree-n part of a homogeneous polynomial.

    In four ambient variables the projector is sum_k a_k R2**k Delta**k,
    where a_0=1 and a_k=-a_(k-1)/(4*k*(n-k+1)).  It terminates at n//2.
    The lower sphere harmonics are recovered by ``spherical_components``.
    """
    if isinstance(degree, bool) or not isinstance(degree, (int, s.Integer)) or degree < 0:
        raise ValueError('a nonnegative integer homogeneous degree is required')
    degree = int(degree)
    p = _poly(polynomial)
    if p.is_zero:
        return s.S.Zero
    if any(sum(powers) != degree for powers, _ in p.terms()):
        raise ValueError('the supplied polynomial is not homogeneous of the stated degree')
    derivative = p.as_expr()
    result = derivative
    coefficient = s.S.One
    for k in range(1, degree//2+1):
        derivative = ambient_laplacian(derivative)
        if derivative == 0:
            break
        coefficient = -coefficient / (4*k*(degree-k+1))
        result += coefficient * R2**k * derivative
    return s.expand(result)


def spherical_components(polynomial):
    """Return {l: homogeneous harmonic H_l} with sum_l H_l = polynomial on S3.

    Every homogeneous input piece is decomposed by successively extracting
    its highest harmonic and dividing the remainder exactly by R2.  The
    degree falls by two on every nonzero remainder, so termination is finite.
    Zero components are omitted and keys are returned in ascending order.
    """
    p = _poly(polynomial)
    if p.is_zero:
        return {}
    homogeneous = {}
    for powers, coefficient in p.terms():
        degree = sum(powers)
        homogeneous[degree] = homogeneous.get(degree, s.S.Zero) + coefficient*_monomial(powers)
    result = {}
    radial = s.Poly(R2, *X, domain='EX')
    for degree in sorted(homogeneous):
        remainder = s.expand(homogeneous[degree])
        for order in range(degree, -1, -2):
            harmonic = harmonic_projection_homogeneous(remainder, order)
            if harmonic != 0:
                result[order] = s.expand(result.get(order, s.S.Zero) + harmonic)
            residue = s.expand(remainder-harmonic)
            if residue == 0:
                break
            if order < 2:
                raise ArithmeticError('nonzero remainder after the constant harmonic')
            quotient, nonradial = _poly(residue).div(radial)
            if not nonradial.is_zero:
                raise ArithmeticError('harmonic remainder was not divisible by R2')
            remainder = quotient.as_expr()
    return {degree: result[degree] for degree in sorted(result) if result[degree] != 0}


@lru_cache(maxsize=4096)
def _monomial_mean(powers):
    if any(n % 2 for n in powers):
        return s.S.Zero
    half = [n//2 for n in powers]
    total = sum(half)
    numerator = prod(factorial(2*n)//(2**n*factorial(n)) for n in half)
    # 4*6*...*(4+2K-2) = 2**K (K+1)!; the empty product is one.
    denominator = 2**total * factorial(total+1)
    return s.Rational(numerator, denominator)


def mean(polynomial):
    """Exact normalized Haar mean on unit S3 for a real polynomial."""
    result = sum(coefficient*_monomial_mean(powers)
                 for powers, coefficient in _poly(polynomial).terms())
    return s.cancel(result)


def inner_product(f, g):
    """Real bilinear normalized L2 pairing; no implicit conjugation."""
    return mean(_poly(f).as_expr()*_poly(g).as_expr())
