"""Exact initial Hermite-to-local-space transfer, followed by one rounding.

Only the physical cap uses source data. The last local polynomial is extended
algebraically to define basis support beyond the wall. The exact target is the
C2 interpolant of the original binary inputs, not rounded PPoly coefficients.
"""
from bisect import bisect_right
from fractions import Fraction
from math import comb

import numpy as np

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_cap_background_jets import RetainedCapBackground
from .nsc_atlas_cap_fixed_basis import FixedCapBasis, exact_number, MAX_CELLS


def _compose(coefficients,offset,scale):
    result=[Fraction(0) for _ in range(6)]
    for p,c in enumerate(coefficients):
        for j in range(p+1): result[j]+=c*comb(p,j)*offset**(p-j)*scale**j
    return tuple(result)


def _hermite(y0,d0,s0,y1,d1,s1,h):
    c0,c1,c2=y0,h*d0,h*h*s0/2
    y,d,s=y1-c0-c1-c2,h*d1-c1-2*c2,h*h*s1-2*c2
    return (c0,c1,c2,10*y-4*d+s/2,-15*y+7*d-s,6*y-3*d+s/2)


def _reconstruct(basis,cell,coefficients):
    idx=basis.indices(cell);matrix=basis.exact_shapes(cell)
    return tuple(sum(row[j]*coefficients[int(i)] for j,i in enumerate(idx)) for row in matrix)


def _transfer_coefficients(basis,pieces,budget):
    result=[Fraction(0) for _ in range(basis.dimension)]
    if any(pieces[0][p] for p in (0,2,4)):
        raise ValueError('initial axis jets do not belong to the odd-quintic subspace')
    result[0]=pieces[0][1]
    for node in range(1,basis.cells+1):
        poly=pieces[node-1]
        h=basis.exact_knots[node]-basis.exact_knots[node-1]
        value=sum(poly)
        derivative=sum(p*poly[p] for p in range(1,6))/h
        result[4*node-3]=value
        result[4*node-2]=basis.exact_scales[node]*derivative
        if node<basis.cells:
            right_h=basis.exact_knots[node+1]-basis.exact_knots[node]
            if value!=pieces[node][0] or derivative!=pieces[node][1]/right_h:
                raise ValueError('exact initial pieces are not C1 at a fixed knot')
    for cell in range(1,basis.cells):
        residual=tuple(a-b for a,b in zip(pieces[cell],_reconstruct(basis,cell,result)))
        linear=residual[5];constant=residual[4]+2*linear
        result[4*cell-1]=constant/16+linear/32
        result[4*cell]=linear/64
    for cell,poly in enumerate(pieces):
        if _reconstruct(basis,cell,result)!=poly:
            raise ValueError('exact local initial polynomial is not represented')
        budget.check(256)
    return tuple(result)


def _upward(value):
    if value==0: return 0.
    result=float(np.nextafter(float(value),np.inf))
    if not np.isfinite(result):
        raise ValueError('coefficient conversion bound is not representable')
    return result


def _round_coefficients(basis,exact,budget):
    values=np.array([float(v) for v in exact])
    if not np.all(np.isfinite(values)):
        raise ValueError('initial coefficient conversion is nonfinite')
    error=tuple(Fraction(float(v))-e for v,e in zip(values,exact))
    value_bound=derivative_bound=Fraction(0)
    for cell in range(basis.cells):
        polynomial=_reconstruct(basis,cell,error)
        h=basis.exact_knots[cell+1]-basis.exact_knots[cell]
        value_bound=max(value_bound,sum(abs(v) for v in polynomial))
        derivative_bound=max(derivative_bound,sum(p*abs(polynomial[p]) for p in range(1,6))/h)
        budget.check(128)
    return values,dict(value=_upward(value_bound),first_derivative=_upward(derivative_bound),
        scope='only coefficient conversion in the exact rational basis; numerical evaluation error excluded')


def transfer_initial_background(background,T,Z,*,budget):
    """Construct one fixed space and the total initial chi/canonical density.

    The input domain is checked in exact-binary arithmetic. Mapped knots
    remain rational even when their differences require more than 53 bits.
    The initial wall is at the middle of the final basis cell. A C2 profile
    whose axis-knot jets violate odd parity is rejected, not smoothed.
    """
    budget=_budget(budget);before=budget.evaluations
    if not isinstance(background,RetainedCapBackground):
        raise ValueError('initial transfer requires the retained original Hermite data')
    T,Z=exact_number(T,'initial time'),exact_number(Z,'wall coordinate')
    if Z<=0:
        raise ValueError('wall coordinate must be positive')
    data={name:tuple(Fraction(float(v)) for v in array) for name,array in background.inputs.items()}
    labels=data['labels']
    quiet=all(data[key][0]==0 for key in ('values','first','second'))
    if T+Z>labels[-1] or (T-Z<labels[0] and not quiet):
        raise ValueError('initial transfer exits the exact supplied profile domain')
    interior=sorted({abs(k-T) for k in labels if 0<abs(k-T)<Z})
    if len(interior)+1>MAX_CELLS:
        raise ValueError('mapped initial profile exceeds the fixed cell limit')
    left=interior[-1] if interior else Fraction(0)
    basis=FixedCapBasis((Fraction(0),*interior,2*Z-left),budget=budget)
    budget.check(64*len(labels))
    cache={}
    chi=[];density=[]
    for cell in range(basis.cells):
        low,high=basis.exact_knots[cell:cell+2]
        h=high-low
        midpoint=(low+min(high,Z))/2
        branches=[]
        for sign in (1,-1):
            query=T+sign*midpoint
            index=bisect_right(labels,query)-1
            if index<0:
                if not quiet: raise ValueError('uncertified lower source extension')
                branches.append(((Fraction(0),)*6,(Fraction(0),)*6))
                continue
            index=min(index,len(labels)-2)
            source_h=labels[index+1]-labels[index]
            if index not in cache:
                cache[index]=_hermite(*(data[key][index] for key in ('values','first','second')),
                    *(data[key][index+1] for key in ('values','first','second')),source_h)
                budget.check(384)
            coeff=cache[index]
            offset=(T+sign*low-labels[index])/source_h
            scale=sign*h/source_h
            first=tuple(p*coeff[p]/source_h for p in range(1,6))
            branches.append((_compose(coeff,offset,scale),_compose(first,offset,scale)))
        chi.append(tuple(a-b for a,b in zip(branches[0][0],branches[1][0])))
        density.append(tuple(a-b for a,b in zip(branches[0][1],branches[1][1])))
        budget.check(1536)
    exact_chi=_transfer_coefficients(basis,chi,budget)
    exact_density=_transfer_coefficients(basis,density,budget)
    numeric_chi,chi_bound=_round_coefficients(basis,exact_chi,budget)
    numeric_density,density_bound=_round_coefficients(basis,exact_density,budget)
    budget.check()
    return dict(basis=basis,chi=numeric_chi,canonical_density=numeric_density,
        exact_chi=exact_chi,exact_canonical_density=exact_density,
        exact_polynomials=dict(chi=tuple(chi),canonical_density=tuple(density)),
        coefficient_conversion_bounds=dict(chi=chi_bound,canonical_density=density_bound),
        initial_time=T,wall_coordinate=Z,evaluations=budget.evaluations-before,
        numerical_evaluation_enclosure=False,physical_evolution=False)
