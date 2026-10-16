"""Fixed physical-z C1 quintics with a regular odd axis cell.

Exact rational knots/shapes define the space. Floating evaluation and Gram
factorization are numerical approximations, not interval error enclosures.
The final wall cut remains in one cell; there is no remeshing or evolution.
"""
from fractions import Fraction
from numbers import Integral

import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.polynomial import polyval
from scipy.linalg import cholesky_banded, cho_solve_banded, eig_banded

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_polar_cap import _scalar, _vector

MAX_CELLS=128


def exact_number(value,name):
    if isinstance(value,Fraction):
        result=value
    else:
        result=Fraction(_scalar(value,name))
    try:
        if not np.isfinite(float(result)):
            raise ValueError(name+' is not representable')
    except OverflowError as error:
        raise ValueError(name+' is not representable') from error
    return result


class FixedCapGram:
    """One clipped-domain Gram factor with its own diagonal scaling."""
    def __init__(self,lower,scale,factor,condition,Z,wall):
        self._lower,self._scale,self._factor=(v.copy() for v in (lower,scale,factor))
        for value in (self._lower,self._scale,self._factor): value.setflags(write=False)
        self.condition=float(condition)
        self.Z=Z
        self.dimension=len(scale)
        self.wall={key:value.copy() for key,value in wall.items()}

    @property
    def lower(self):
        return self._lower.copy()

    def squared_norm(self,coefficients,*,budget):
        """Numerical integral of the squared field via the scaled factor.

        M=S L L^T S, so ||L^T S c||^2 avoids cancellation in c^T M c.
        This uses the stored numerical Gram, not an interval enclosure.
        """
        budget=_budget(budget)
        coefficients=_vector(coefficients,'Gram norm coefficients')
        if coefficients.shape!=(self.dimension,):
            raise ValueError('Gram norm coefficient dimension differs')
        budget.check(24*self.dimension)
        scaled=self._scale*coefficients
        product=self._factor[0]*scaled
        for offset in range(1,len(self._factor)):
            product[:-offset]+=self._factor[offset,:self.dimension-offset]*scaled[offset:]
        value=float(product@product)
        if not np.isfinite(value): raise ValueError('nonfinite Gram squared norm')
        return value

    def solve(self,rhs,*,budget):
        budget=_budget(budget)
        rhs=_vector(rhs,'Gram right hand side')
        if rhs.shape!=(self.dimension,):
            raise ValueError('Gram right hand side has the wrong dimension')
        budget.check(48*self.dimension)
        value=cho_solve_banded((self._factor,True),rhs/self._scale,check_finite=False)/self._scale
        if not np.all(np.isfinite(value)):
            raise ValueError('nonfinite Gram solution')
        product=self._lower[0]*value
        row_sum=abs(self._lower[0]).copy()
        for offset in range(1,len(self._lower)):
            lower=self._lower[offset,:self.dimension-offset]
            product[offset:]+=lower*value[:-offset]
            product[:-offset]+=lower*value[offset:]
            row_sum[offset:]+=abs(lower);row_sum[:-offset]+=abs(lower)
        denominator=float(np.max(row_sum)*np.max(abs(value))+np.max(abs(rhs)))
        residual=float(np.max(abs(product-rhs))/denominator) if denominator else 0.
        if not np.isfinite(denominator) or not np.isfinite(residual) or residual>64*np.finfo(float).eps:
            raise ValueError('Gram solve residual is unresolved')
        return value


class FixedCapBasis:
    """At most 128 fixed cells; 4*N-1 total-field degrees of freedom.

    Node i>=1 has value at 4*i-3 and scaled z derivative at 4*i-2.
    Cells are numbered from zero: cell 0 is the special axis cell.
    Ordinary cell i>=1 has bubbles at 4*i-1,4*i. Index zero is h0 times
    the axis derivative. The first polynomial is odd for every state.
    """
    def __init__(self,knots,*,budget):
        budget=_budget(budget)
        knots=tuple(exact_number(v,'knot') for v in knots)
        if not 2<=len(knots)<=MAX_CELLS+1 or knots[0]!=0 or any(b<=a for a,b in zip(knots,knots[1:])):
            raise ValueError('basis requires 1..128 positive-width cells starting at zero')
        numeric=np.array([float(v) for v in knots])
        if np.any(np.diff(numeric)<=0):
            raise ValueError('distinct exact knots are not resolved in floating coordinates')
        widths=tuple(b-a for a,b in zip(knots,knots[1:]))
        if any(float(h)<=0 for h in widths):
            raise ValueError('cell width underflows')
        self._knots,self._widths=knots,widths
        self.cells=len(widths)
        self.dimension=4*self.cells-1
        self._scales=(widths[0],)+tuple(min(widths[i-1],widths[i]) for i in range(1,self.cells))+(widths[-1],)
        self._exact=[];self._numeric=[];self._first=[];self._second=[];self._indices=[]
        for cell,h in enumerate(widths):
            if cell==0:
                columns=[[0,1,0,-2,0,1],[0,0,0,Fraction(5,2),0,Fraction(-3,2)],
                         [0,0,0,-h/(2*self._scales[1]),0,h/(2*self._scales[1])]]
                indices=np.arange(3)
            else:
                left,right=h/self._scales[cell],h/self._scales[cell+1]
                columns=[[1,0,-3,2,0,0],[0,left,-2*left,left,0,0],
                         [0,0,16,-32,16,0],[0,0,-32,128,-160,64],
                         [0,0,3,-2,0,0],[0,0,-right,right,0,0]]
                indices=np.arange(4*cell-3,4*cell+3)
            exact=tuple(tuple(Fraction(columns[j][p]) for j in range(len(columns))) for p in range(6))
            numeric=np.array(exact,dtype=float)
            first=np.arange(1,6)[:,None]*numeric[1:]/float(h)
            second=np.arange(1,5)[:,None]*first[1:]/float(h)
            if not all(np.all(np.isfinite(v)) for v in (numeric,first,second)):
                raise ValueError('basis derivatives are not representable')
            for value in (numeric,first,second,indices): value.setflags(write=False)
            self._exact.append(exact);self._numeric.append(numeric)
            self._first.append(first);self._second.append(second);self._indices.append(indices)
            budget.check(256)

    @property
    def exact_knots(self): return self._knots

    @property
    def exact_scales(self): return self._scales

    @property
    def knots(self): return np.array([float(v) for v in self._knots])

    def _cell(self,cell):
        if isinstance(cell,(bool,np.bool_)) or not isinstance(cell,Integral) or not 0<=cell<self.cells:
            raise ValueError('invalid fixed basis cell')
        return int(cell)

    def exact_shapes(self,cell):
        return self._exact[self._cell(cell)]

    def indices(self,cell):
        return self._indices[self._cell(cell)].copy()

    def values(self,cell,xi,*,budget):
        budget=_budget(budget);cell=self._cell(cell)
        xi=_scalar(xi,'local coordinate')
        if not 0<=xi<=1:
            raise ValueError('local coordinate must lie in [0,1]')
        budget.check(32+18*len(self._indices[cell]))
        return dict(indices=self.indices(cell),value=polyval(xi,self._numeric[cell]),
                    first=polyval(xi,self._first[cell]),second=polyval(xi,self._second[cell]))

    def power_coefficients(self,coefficients,*,budget):
        budget=_budget(budget);coefficients=_vector(coefficients,'field coefficients')
        if len(coefficients)!=self.dimension:
            raise ValueError('field coefficient dimension differs from basis')
        budget.check(48*self.cells)
        result=np.array([matrix@coefficients[idx] for matrix,idx in zip(self._numeric,self._indices)])
        if not np.all(np.isfinite(result)):
            raise ValueError('nonfinite local field polynomial')
        return result

    def gram(self,Z,*,budget):
        """Numerical degree-ten integrals and a scaled banded SPD solve."""
        budget=_budget(budget);Z=exact_number(Z,'wall coordinate')
        cut=(Z-self._knots[-2])/self._widths[-1]
        if not Fraction(1,4)<=cut<=Fraction(3,4):
            raise ValueError('last-cell wall cut must remain in [1/4,3/4]')
        n=self.dimension;band=min(5,n-1)
        lower=np.zeros((band+1,n))
        nodes,weights=leggauss(6)
        for cell,h in enumerate(self._widths):
            fraction=float(cut) if cell==self.cells-1 else 1.
            local=np.zeros((len(self._indices[cell]),)*2)
            for xi,weight in zip(fraction*(nodes+1)/2,weights):
                value=self.values(cell,float(xi),budget=budget)['value']
                local+=float(h)*fraction*weight*np.outer(value,value)/2
            idx=self._indices[cell]
            for i,gi in enumerate(idx):
                for j,gj in enumerate(idx[:i+1]): lower[gi-gj,gj]+=local[i,j]
            budget.check(300)
        scale=np.sqrt(lower[0])
        if np.any(scale<=0) or not np.all(np.isfinite(lower)):
            raise ValueError('Gram diagonal is unresolved')
        normalized=lower.copy()
        for offset in range(band+1):
            normalized[offset,:n-offset]/=scale[offset:]*scale[:n-offset]
        budget.check(20*n*(band+1)**2)
        spectrum=eig_banded(normalized,lower=True,eigvals_only=True,check_finite=True)
        condition=float(spectrum[-1]/spectrum[0])
        if spectrum[0]<=0 or not np.isfinite(condition) or condition*np.finfo(float).eps>1e-8:
            raise ValueError('normalized Gram conditioning is outside the numerical domain')
        factor=cholesky_banded(normalized,lower=True,check_finite=True)
        wall=self.values(self.cells-1,float(cut),budget=budget)
        return FixedCapGram(lower,scale,factor,condition,Z,wall)
