"""Exact small fixed-space matrices for the zero-coupling time scale.

This is a frozen linear cap bound. It does not bound the nonlinear coupled
generator, source error, moving energy metric, or spatial approximation.
"""
from fractions import Fraction as Q
from math import isqrt

import sympy as s

from .nsc_atlas_cap_fixed_basis import FixedCapBasis,exact_number
from .nsc_atlas_background_cap import _budget


def frozen_cap_bound(basis,Z,*,budget):
    """Integrate the exact basis products; return M,K,e and rational L_*^2.

    Restricted to at most four cells for this initial-time construction.
    Exact arithmetic is charged before dense inversion. Returned matrices
    are for in-memory diagnostics, not a raw-array publication artifact.
    """
    budget=_budget(budget)
    if not isinstance(basis,FixedCapBasis) or not 1<=basis.cells<=4:
        raise ValueError('frozen matrix calculation requires 1..4 fixed cells')
    Z=exact_number(Z,'wall coordinate');knots=basis.exact_knots
    cut=(Z-knots[-2])/(knots[-1]-knots[-2])
    if not Q(1,4)<=cut<=Q(3,4): raise ValueError('wall cut leaves fixed domain')
    n=basis.dimension;budget.check(64*n**3+128*n*n*basis.cells)
    M=s.zeros(n);K=s.zeros(n);e=s.zeros(n,1)
    for cell in range(basis.cells):
        width=knots[cell+1]-knots[cell]
        top=cut if cell==basis.cells-1 else Q(1)
        shape=basis.exact_shapes(cell);idx=basis.indices(cell)
        for i,gi in enumerate(idx):
            for j,gj in enumerate(idx):
                mass=sum(shape[p][i]*shape[q][j]*top**(p+q+1)/Q(p+q+1)
                         for p in range(6) for q in range(6))*width
                stiff=sum(p*q*shape[p][i]*shape[q][j]*top**(p+q-1)/Q(p+q-1)
                          for p in range(1,6) for q in range(1,6))/width
                M[int(gi),int(gj)]+=s.Rational(mass.numerator,mass.denominator)
                K[int(gi),int(gj)]+=s.Rational(stiff.numerator,stiff.denominator)
            if cell==basis.cells-1:
                value=sum(shape[p][i]*top**p for p in range(6))
                e[int(gi)]=s.Rational(value.numerator,value.denominator)
    inverse=M.inv(method='DM')
    q=(e.T*inverse*e)[0]
    squared=2*s.trace(inverse*K)+q*q
    squared=Q(int(squared.p),int(squared.q))
    if squared<=0: raise ValueError('nonpositive frozen generator bound')
    integer=isqrt(squared.numerator//squared.denominator)
    if integer*integer<squared: integer+=1
    budget.check()
    return dict(M=M,K=K,e=e,bound_squared=squared,integer_bound=integer,
                dimension=n,scope='frozen alpha=0 cap, unit T_x; coupled stability not established')
