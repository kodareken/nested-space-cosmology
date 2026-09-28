"""Finite endpoint identities for the nonlinear preparation argument.

This checks finite changes in the first unmatched time jet, with arbitrary
common lower jets. It does not compute a physical covariance, pseudodifferential
remainder constant, or a new reference subtraction.
"""
from functools import lru_cache
import sympy as sp


def zero(value):
    return str(sp.factor(sp.cancel(value)))


@lru_cache(maxsize=1)
def finite_endpoint_maps():
    a, r = sp.symbols('a1 r1', positive=True)
    m, ell = sp.symbols('m ell', real=True, nonzero=True)
    da, dr = sp.symbols('Delta_a_k Delta_r_k', real=True)
    rows = []
    for k in range(1, 5):
        aa = [a] + [sp.Symbol('a_' + str(j), real=True)/sp.factorial(j)
                    for j in range(1, k+1)]
        rr = [r] + [sp.Symbol('r_' + str(j), real=True)/sp.factorial(j)
                    for j in range(1, k+1)]

        def coefficient(axial, radius, conjugate=False):
            # Actual Taylor coefficient of a/r; retain every lower-jet product.
            inverse = [1/radius[0]]
            for n in range(1, k+1):
                inverse.append(-sum(radius[j]*inverse[n-j]
                                    for j in range(1, n+1))/radius[0])
            ratio = sum(axial[j]*inverse[k-j] for j in range(k+1))
            return (m*axial[k] + (-sp.I if conjugate else sp.I)*ell*ratio)/2

        changed_a, changed_r = list(aa), list(rr)
        changed_a[k] += da/sp.factorial(k)
        changed_r[k] += dr/sp.factorial(k)
        d01 = sp.factor(sp.factorial(k)*(coefficient(changed_a, changed_r)-coefficient(aa, rr)))
        d10 = sp.factor(sp.factorial(k)*(coefficient(changed_a, changed_r, True)-coefficient(aa, rr, True)))
        expected01 = ((m+sp.I*ell/r)*da-sp.I*ell*a*dr/r**2)/2
        expected10 = ((m-sp.I*ell/r)*da+sp.I*ell*a*dr/r**2)/2
        limits = [sp.factor((-sp.I*a/2)**k*d01),
                  sp.factor((sp.I*a/2)**k*d10)]
        expected = [-a**k/(2*sp.I)**(k+1)*((ell/r-sp.I*m)*da-ell*a*dr/r**2),
                    (-1)**k*a**k/(2*sp.I)**(k+1)*((ell/r+sp.I*m)*da-ell*a*dr/r**2)]
        matrix = sp.Matrix([[sp.diff(entry, v) for v in (da, dr)] for entry in limits])
        determinant = sp.factor(matrix.det())
        target = sp.I*m*ell*a**(2*k+1)/(2**(2*k+1)*r**2)
        residuals = [zero(d01-expected01), zero(d10-expected10),
                     *[zero(x-y) for x,y in zip(limits, expected)],
                     zero(determinant-target), zero(limits[1]-sp.conjugate(limits[0]))]
        if set(residuals) != {'0'}:
            raise ArithmeticError('finite endpoint identity failed')
        rows.append({'normal_order': k, 'finite_difference_of_p1_Tk': [str(d01), str(d10)],
                     'finite_endpoint_coefficient': [str(v) for v in limits],
                     'coefficient_map_determinant': str(determinant),
                     'nonlinear_lower_jets_retained': True,
                     'amplitude_derivative_used': False, 'exact_residuals': residuals})
    return rows


def signed_source_identity():
    f, s, n = sp.symbols('f s n', real=True)
    positive = sp.Matrix([[f, -sp.I*s, 0], [sp.I*s, 1-f, 0], [0, 0, n]])
    negative = positive.subs({f: 1-f, n: 1-n}, simultaneous=True)
    residuals = [zero(v) for v in negative-(sp.eye(3)-sp.conjugate(positive))]
    if set(residuals) != {'0'}:
        raise ArithmeticError('signed physical source identity failed')
    return {'identity': 'Csrc(-E)=I-conj(Csrc(E))',
            'modal_identity': 'C(-E,ell)=I-sigma3*conj(C(E,-ell))*sigma3 by the owned mode map and coisometry',
            'tail_receipt_sign_meaning': 'angular sign, not energy sign',
            'exact_residuals': residuals}


def symbolic_proof():
    rows = finite_endpoint_maps()
    source = signed_source_identity()
    return {'finite_endpoint_maps': rows, 'signed_source': source,
            'exact_scalar_residual_count': sum(len(row['exact_residuals']) for row in rows)+len(source['exact_residuals']),
            'calculus_scope': 'finite coefficient and source identities only; the physical operator bridge and triangular derivative count are analytic arguments'}
