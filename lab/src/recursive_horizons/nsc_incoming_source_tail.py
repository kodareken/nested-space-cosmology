"""Inverse-energy tail of the existing paired order16/ad4 source recipe.

The finite NSC Riccati recurrence and the already generated adiabatic Bloch
expressions are evaluated with arbitrary precision. Taylor coefficients come
from those expressions, never a fit to measured spectral values. This owns
the truncation estimate for that specified approximation only; the true
scattering/reflection remainder and full renormalized source remain separate.
"""
from functools import lru_cache
import types

import mpmath as mp
import numpy as np
import sympy as sp

from .nsc_compact_ctp_neck import _adiabatic_bloch_function


@lru_cache(maxsize=1)
def _mp_bloch():
    """Same generated CSE expressions, replacing numerical backend only."""
    original = _adiabatic_bloch_function()
    environment = dict(original.__globals__)
    environment.update(sin=mp.sin, cos=mp.cos, sqrt=mp.sqrt, pi=mp.pi,
                       array=lambda value: np.asarray(value, dtype=object))
    return types.FunctionType(original.__code__, environment)


def _riccati_at_one(mass, angular, order=16):
    """Owned riccati_coefficients recurrence, same rho=1 Taylor coefficients.

    Arithmetic follows nsc_pg_high_energy.riccati_coefficients verbatim: no
    equation or state is changed by retaining precision in this local algebra.
    """
    size = order+2
    zero = lambda: [mp.mpc(0) for _ in range(size)]
    def mul(a, b): return [sum((a[j]*b[n-j] for j in range(n+1)), mp.mpc(0)) for n in range(size)]
    def derivative(a): return [(n+1)*a[n+1] for n in range(size-1)]+[mp.mpc(0)]
    r2 = zero(); r2[:3] = [mp.mpf(2), mp.mpf(2), mp.mpf(1)]
    inverse = zero(); inverse[0] = 1/r2[0]
    for n in range(1, size):
        inverse[n] = -sum((r2[j]*inverse[n-j] for j in range(1, n+1)), mp.mpc(0))/r2[0]
    sphere = zero(); sphere[0] = mp.sqrt(inverse[0])
    for n in range(1, size):
        sphere[n] = (inverse[n]-sum((sphere[j]*sphere[n-j] for j in range(1, n)), mp.mpc(0)))/(2*sphere[0])
    angle = zero(); angle[0] = -mp.pi/4
    for n in range(1, size): angle[n] = inverse[n-1]/n
    A = [3*v for v in mul(r2, angle)]; A[0] += 4; A[1] += 3
    upper = [angular*v for v in sphere]; upper[0] += 1j*mass
    lower = [mp.conj(v) for v in upper]
    coefficients = [upper]; dA = derivative(A)
    for n in range(1, order):
        last = coefficients[-1]
        first, second = mul(A, derivative(last)), mul(dA, last)
        product = zero()
        for j in range(1, n):
            term = mul(coefficients[j-1], coefficients[n-j-1])
            product = [a+b for a, b in zip(product, term)]
        nonlinear = mul(mul(A, lower), product)
        coefficients.append([1j*(a+b/2)+c for a, b, c in zip(first, second, nonlinear)])
    return [row[0] for row in coefficients]


def paired_source_function(channel):
    """Existing vacuum source summed over actual angular signs, without factor.

    x=1/E. The massless LLL is explicitly excluded. This is local source
    algebra only; zero x denotes the common zero high-energy limit.
    """
    mass = mp.mpf(float(channel['compact_mass']))
    angular = mp.mpf(float(channel['angular_eigenvalue']))
    if mass < 0 or angular < 0 or mass == angular == 0:
        raise ValueError('an existing massive retained group is required')
    signs = (1,) if angular == 0 else (1, -1)
    radius, axial = mp.sqrt(2), mp.sqrt(3*mp.pi/2-4)
    coefficients = {sign: _riccati_at_one(mass, sign*angular) for sign in signs}
    backend = _mp_bloch()

    def source(x):
        x = mp.mpf(x)
        if not x: return mp.matrix([0, 0, 0])
        result = mp.matrix([0, 0, 0])
        for sign in signs:
            ell = sign*angular
            c = coefficients[sign]
            # Preserve the exact leading identity at the precision selected
            # by differentiation, rather than amplifying a rounded c1.
            series = (ell/radius+1j*mass)*x/2
            series += sum((c[j]*(x/2)**(j+1) for j in range(1, 16)), mp.mpc(0))
            u = axial**2*abs(series)**2
            bx, by = 2*axial*series.imag/(1+u), -2*axial*series.real/(1+u)
            y = axial**2*(mass**2+ell**2/radius**2)*x*x
            root = mp.sqrt(1+y)
            delta = [bx-mass*axial*x/root, by+ell*axial*x/(radius*root),
                     y/(root*(1+root))-2*u/(1+u)]
            terms = backend(3*mp.pi/4, mass, ell, 1/x)
            for component in range(3):
                delta[component] -= sum((terms[j][component, 0] for j in range(1, 5)), mp.mpf(0))
            p = -delta[2]/(axial*x)
            result += mp.matrix([-mass*delta[0]+ell/radius*delta[1]+p,
                                 p, ell/(2*radius)*delta[1]])
        return result

    return source


@lru_cache(maxsize=1)
def leading_reference_matching_identity():
    """Small exact coefficient check fixing the regularized x=0 endpoint.

    Write L=lambda/r, A=-a^2, dot a=A_rho/2, dot L=-a L_rho.
    The owned Riccati c2 and adiabatic b1 give identical transverse x^2
    and longitudinal x^3 coefficients, separately for each angular sign.
    Together with c1 matching P0, this proves source/x^2 ->0; it is not a
    limit inferred from chopped numerical Taylor coefficients.
    """
    a, m, L, Lrho, Aprime = sp.symbols('a m L Lrho Aprime', real=True)
    c1 = L+sp.I*m
    c2 = -m*Aprime/2+sp.I*(-a*a*Lrho+Aprime*L/2)
    vacuum = sp.Matrix([a*sp.im(c2)/2, -a*sp.re(c2)/2,
                        -a*a*sp.re(c1*sp.conjugate(c2))/2])
    adot, Ldot = Aprime/2, -a*Lrho
    # Leading coefficients of the already-owned b1=-h x dot(b0)/(2 h^2).
    adiabatic = sp.Matrix([a*(adot*L+a*Ldot)/2, m*a*adot/2, -m*a**3*Ldot/2])
    residual = [sp.simplify(value) for value in vacuum-adiabatic]
    if residual != [0, 0, 0]:
        raise ArithmeticError('regularized incoming-source endpoint is not established')
    return {'powers': [2, 2, 3], 'vacuum_minus_ad1': [str(value) for value in residual],
            'c1': 'L+i*m', 'clock': 'dot(a)=A_rho/2; dot(L)=-a*L_rho; A=-a^2',
            'scope': 'each angular sign; local existing Riccati/ad1 coefficient identity',
            'reduced_endpoint': 'source(x)/x^2 -> 0 after P0 and ad1 subtraction'}


def _validate_owned_split(channel, lower):
    group = channel.get('index')
    if group not in range(1, 33) or lower != (320. if group in (10, 11, 12, 31, 32) else 160.):
        raise ValueError('the retained group and its archived middle/tail split must agree')
    return group


def incoming_source_tail(channel, lower, *, maximum_order=11, precision=60):
    """Analytic integral of Taylor coefficients; numerical remainder estimate.

    Returns rho, p_parallel, T01=0 vacuum, p_perp. The omitted thermal tail
    and the physical departure from the order16 vacuum mode recipe are NOT
    included or claimed bounded. lower is the existing spectral split.
    """
    if lower not in (160., 320.) or maximum_order not in (9, 11, 13) or precision < 40:
        raise ValueError('owned lower split160/320, Taylor order9/11/13 and precision>=40 required')
    group = _validate_owned_split(channel, lower)
    with mp.workdps(precision):
        source = paired_source_function(channel)
        endpoint_identity = leading_reference_matching_identity()
        # No backend chopping. The two nonintegrable coefficients are exact
        # structural zeros of x^2*g(x); g(0)=0 follows from the check above.
        regular = lambda x: mp.matrix([0, 0, 0]) if not x else source(x)/(x*x)
        reduced = mp.taylor(regular, mp.mpf(0), maximum_order-2, direction=1, chop=False)
        coefficients = [mp.matrix([0, 0, 0]), mp.matrix([0, 0, 0])]+reduced
        raw_first = mp.diff(source, mp.mpf(0), 1, direction=1)
        coefficient_array = np.array([[float(v) for v in row] for row in coefficients])
        copies = mp.mpf(int(channel['copy_count'])); degeneracy = mp.mpf(int(channel['degeneracy']))
        nsigns = 1 if channel['angular_eigenvalue'] == 0 else 2
        factor = copies*degeneracy/(4*mp.pi**2*2*mp.sqrt(3*mp.pi/2-4)*nsigns)
        corrections = {}
        for order in range(7, maximum_order+1, 2):
            value = factor*sum((coefficients[n]/((n-1)*mp.mpf(lower)**(n-1))
                                for n in range(2, order+1)), mp.matrix([0, 0, 0]))
            corrections[order] = np.array([float(value[0]), float(value[1]), 0., float(value[2])])
        if any(value != 0 for row in coefficients[:2] for value in row):
            raise ArithmeticError('inverse-energy nonintegrable coefficient unresolved; no coefficient clipping allowed')
        nonintegrable = 0.
        return {'coefficients_rho_parallel_sphere': coefficient_array,
                'coefficient_scope': 'unweighted paired vacuum kernels; actual angular signs already summed',
                'multiplicity_factor': float(factor),
                'corrections_by_order': corrections,
                'series_order_difference': abs(corrections[maximum_order]-corrections[maximum_order-2]),
                'nonintegrable_coefficient_residual': nonintegrable,
                'nonintegrable_coefficients_exactly_zero': True,
                'regular_endpoint_identity': endpoint_identity,
                'unfactored_first_derivative_numerical_diagnostic': [mp.nstr(v, 30) for v in raw_first],
                'numerical_differentiation_chop': False,
                'lower': float(lower), 'maximum_order': maximum_order, 'precision': precision,
                'group': int(group), 'riccati_order': 16,
                'scope': 'paired order16/ad4 vacuum approximation only; no full physical stress-tail bound',
                'rigorous_physical_remainder_bound': None}


def evaluate_tail_series(result, energies, *, order=None):
    """Evaluate saved Taylor coefficients, not fit them to spectral samples."""
    E = np.asarray(energies, float)
    order = result['maximum_order'] if order is None else order
    if (E.ndim != 1 or not len(E) or np.any(E <= 0) or not np.isfinite(E).all()
            or not 2 <= order <= result['maximum_order']):
        raise ValueError('positive real frequencies and an available Taylor order required')
    c = np.asarray(result['coefficients_rho_parallel_sphere'])
    value = (E[:, None]**(-np.arange(order+1)[None, :]))@c[:order+1]
    value *= result['multiplicity_factor']
    return np.column_stack((value[:, :2], np.zeros(len(E)), value[:, 2]))


def transformed_tail_quadrature(channel, lower, *, points=24, precision=60):
    """Independent numeric integral of the same existing approximation.

    x=1/E maps its infinite interval to [0,1/L]. Gaussian nodes never evaluate
    at x=0. Increasing points and precision is a convergence estimate, not a
    proof of the true scattering/reflection remainder.
    """
    if lower not in (160., 320.) or points not in (16, 24, 32) or precision < 40:
        raise ValueError('owned split, bounded quadrature16/24/32 and precision>=40 required')
    _validate_owned_split(channel, lower)
    nodes, weights = np.polynomial.legendre.leggauss(points)
    with mp.workdps(precision):
        source = paired_source_function(channel)
        total = mp.matrix([0, 0, 0])
        for node, weight in zip(nodes, weights):
            x = (mp.mpf(float(node))+1)/(2*mp.mpf(lower))
            total += mp.mpf(float(weight))/(2*mp.mpf(lower))*source(x)/(x*x)
        signs = 1 if channel['angular_eigenvalue'] == 0 else 2
        factor = (mp.mpf(int(channel['copy_count']))*int(channel['degeneracy'])
                  /(4*mp.pi**2*2*mp.sqrt(3*mp.pi/2-4)*signs))
        value = factor*total
        return np.array([float(value[0]), float(value[1]), 0., float(value[2])])
