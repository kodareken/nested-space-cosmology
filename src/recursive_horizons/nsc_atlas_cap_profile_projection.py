"""One knot-split projection of a supplied computational cap profile.

Quadrature orders integrate the represented polynomials exactly in exact
arithmetic. Floating results are diagnostics, not error enclosures for the
underlying physical preparation or a finite-amplitude evolution.
"""
import math

import numpy as np
from numpy.polynomial.legendre import leggauss

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_characteristic_background import CharacteristicProfile
from .nsc_atlas_cap_polynomial_basis import odd_legendre_operators, odd_legendre_values
from .nsc_atlas_polar_cap import _scalar


def _values(profile, T, Z, xi, budget):
    jet = profile.jet(np.r_[T+Z*xi,T-Z*xi],budget=budget)
    size = len(xi)
    plus,minus = slice(0,size),slice(size,None)
    return dict(chi=jet['F'][plus]-jet['F'][minus],
                k=jet['F_prime'][plus]-jet['F_prime'][minus],
                chi_xi=Z*(jet['F_prime'][plus]+jet['F_prime'][minus]),
                k_xi=Z*(jet['F_second'][plus]+jet['F_second'][minus]))


def project_characteristic_profile(profile, T, Z, *, count=16, budget):
    """Project chi=F(T+z)-F(T-z), k=chi_T on 0<=z<=Z.

    The returned L2 norms use dz; derivative seminorms use d/dz. The input
    is the existing C2 piecewise quintic, including its certified quiet
    extension. No extrapolation, preparation, radial solve or evolution.
    """
    budget = _budget(budget)
    before = budget.evaluations
    if not isinstance(profile,CharacteristicProfile):
        raise ValueError('require one CharacteristicProfile')
    T,Z = _scalar(T,'T'),_scalar(Z,'Z',positive=True)
    operators = odd_legendre_operators(count,Z,budget=budget)
    # Validate the complete label interval before integrating any panel.
    profile.jet(np.array([T-Z,T,T+Z]),budget=budget)
    R,A0 = float(np.tanh(Z)),float(1/np.cosh(Z)**2)
    if not 0 < R < 1 or not A0 > 0:
        raise ValueError('wall radius must be representable inside (0,1)')
    knots = profile.breakpoints
    mapped = np.r_[(knots-T)/Z,(T-knots)/Z]
    cuts = np.unique(np.r_[0.,mapped[(mapped>0)&(mapped<1)],1.])
    budget.check(4*len(knots)+2)
    coefficient_terms = {key:[] for key in ('chi','k')}
    norms = {key:[] for key in ('chi','k','chi_xi','k_xi')}
    reference_norms = {key:[] for key in norms}
    coefficients = {}
    for order in (19,32):
        nodes,weights = leggauss(order)
        budget.check(order*order)
        for low,high in zip(cuts[:-1],cuts[1:]):
            half = (high-low)/2
            xi = low+half*(nodes+1)
            values = _values(profile,T,Z,xi,budget)
            basis = [odd_legendre_values(count,float(x),budget=budget) for x in xi]
            P = np.array([row['value'] for row in basis])
            D = np.array([row['derivative_xi'] for row in basis])
            budget.check(order*(8*count+20))
            if order == 19:
                for key in coefficient_terms:
                    coefficient_terms[key].append(half*(P.T@(weights*values[key])))
            else:
                for key in ('chi','k'):
                    for suffix,matrix in (('',P),('_xi',D)):
                        name = key+suffix
                        residual = matrix@coefficients[key]-values[name]
                        norms[name].append(float(half*(weights@(residual*residual))))
                        reference_norms[name].append(float(half*(weights@(values[name]*values[name]))))
        if order == 19:
            factors = 2*operators['powers']+1
            for key,terms in coefficient_terms.items():
                coefficients[key] = factors*np.array([math.fsum(v[j] for v in terms) for j in range(count)])
    def physical_norms(terms):
        return {key:math.sqrt(math.fsum(values)*(1/Z if key.endswith('_xi') else Z))
                for key,values in terms.items()}
    original = _values(profile,T,Z,np.array([1.]),budget)
    projected = {}
    for key in ('chi','k'):
        projected[key] = float(operators['wall_value']@coefficients[key])
        projected[key+'_xi'] = float(operators['wall_derivative_xi']@coefficients[key])
    def wall(values):
        chi,k,derivative = (float(np.asarray(values[key]).reshape(-1)[0]) for key in ('chi','k','chi_xi'))
        f,p = chi/R,k/A0
        h = derivative/(Z*A0)-f
        return dict(chi=chi,k=k,f=f,p=p,h=h,Pi=p/R,Phi=h/R,
                    incoming=(p+h)/R,outgoing=(p-h)/R)
    wall_original,wall_projected = wall(original),wall(projected)
    axis = profile.jet(T,budget=budget)
    centre_original = dict(f=2*float(axis['F_prime']),Pi=2*float(axis['F_second']))
    centre_projected = dict(f=float(operators['centre_derivative_xi']@coefficients['chi']/Z),
                            Pi=float(operators['centre_derivative_xi']@coefficients['k']/Z))
    result = dict(T=T,Z=Z,radius=R,A0=A0,count=count,panels=len(cuts)-1,
        coefficients=coefficients,residual_norms=physical_norms(norms),
        reference_norms=physical_norms(reference_norms),
        wall_original=wall_original,wall_projected=wall_projected,
        wall_error={key:wall_projected[key]-wall_original[key] for key in wall_original},
        centre_original=centre_original,centre_projected=centre_projected,
        centre_error={key:centre_projected[key]-centre_original[key] for key in centre_original},
        coefficient_order=19,residual_order=32,evaluated_points=budget.evaluations-before,
        numerical_error_enclosure=False,
        scope='projection of the supplied computational C2 profile; physical preparation error excluded')
    numeric = [T,Z,R,A0,*coefficients.values(),*result['residual_norms'].values(),
               *result['reference_norms'].values(),*wall_original.values(),*wall_projected.values(),
               *centre_original.values(),*centre_projected.values()]
    if not all(np.all(np.isfinite(value)) for value in numeric):
        raise ValueError('nonfinite profile projection')
    budget.check()
    return result
