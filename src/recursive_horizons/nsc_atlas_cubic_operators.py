"""Exact finite scalar maps for the inherited cubic backreaction equation.

All inputs are the round-S3 basis and the specified eleven metric tensors.
There is no numerical projection or new harmonic cutoff: every polynomial
image is reconstructed exactly modulo R2-1 before conversion to float.
The TT normalization is the G22 remainder documented below.
"""
from functools import lru_cache
from types import MappingProxyType

import numpy as np
import sympy as s

from .nsc_atlas_s3_polynomials import X,R2,laplacian,gradient_dot

C1=s.sqrt(2)/s.pi
C2=s.sqrt(12)/s.pi


@lru_cache(maxsize=1)
def scalar_basis():
    z,x,_,y=X
    return (s.S.One,z*z-s.Rational(1,4),x*x+y*y-s.Rational(1,2),z*x*y,
            x*x*y*y-(x*x+y*y)/8+s.Rational(1,48))


@lru_cache(maxsize=1)
def cubic_basis():
    z,x,_,y=X
    return (z,x*y,z**3,z*(x*x+y*y),z*z*x*y,x*y*(x*x+y*y),z*x*x*y*y,x**3*y**3)


def _sphere_reduce(expr):
    # Eliminate the unused X3 coordinate instead of a retained source axis.
    generators=(X[2],X[0],X[1],X[3])
    return s.expand(s.Poly(s.expand(expr),*generators,domain='EX').rem(
        s.Poly(R2-1,*generators,domain='EX')).as_expr())


@lru_cache(maxsize=1)
def _coefficient_solver():
    P=cubic_basis()
    monomials=sorted(set().union(*(set(s.Poly(f,*X).monoms()) for f in P)))
    columns=[s.Poly(f,*X,domain='EX') for f in P]
    matrix=s.Matrix([[f.coeff_monomial(m) for f in columns] for m in monomials])
    # Exact independent coefficient equations, not a least-squares projection.
    rows=matrix.T.rref()[1]
    if len(rows)!=len(P):raise ArithmeticError('the cubic polynomial basis lost rank')
    inverse=s.ImmutableMatrix(matrix.extract(rows,range(len(P))).inv())
    return tuple(monomials[i] for i in rows),inverse


def _coefficients(expr):
    expression=_sphere_reduce(expr)
    polynomial=s.Poly(expression,*X,domain='EX')
    monomials,inverse=_coefficient_solver()
    coefficients=inverse*s.Matrix([polynomial.coeff_monomial(m) for m in monomials])
    reconstructed=sum(c*f for c,f in zip(coefficients,cubic_basis()))
    if _sphere_reduce(expression-reconstructed)!=0:
        raise ArithmeticError('the supplied operator image lies outside the cubic basis')
    return tuple(s.simplify(c) for c in coefficients)


def _W_derivative(f):
    z,x,_,y=X
    W=2*s.sqrt(6)/(3*s.pi**2)*s.Matrix([-4*x*y,2*z*y,0,2*z*x])
    # X.W=0 identically; the ambient derivative is already the tangent one.
    return s.expand(sum(W[i]*s.diff(f,X[i]) for i in range(4)))


@lru_cache(maxsize=1)
def exact_operator_matrices():
    """Return immutable exact maps; output polynomial coefficient is first.

    Shapes: laplacian (8,8), metric_delta_laplacian (8,11,2),
    multiplication (8,5,2), trace (5,11). The metric order is gamma*S0..4,
    Hess*S1..4, symgradW, Z. Z=G22 minus C2**2 times
    [gamma/9+gamma*S2b/2+Hess(S2b)/8+4gamma*S4/3+Hess(S4)/3].
    """
    S=scalar_basis();P=cubic_basis();Y=(C1*X[0],C2*X[1]*X[3])
    delta=s.zeros(8,8)
    for j,f in enumerate(P):
        delta[:,j]=s.Matrix(_coefficients(laplacian(f)))
    metric=s.MutableDenseNDimArray.zeros(8,11,2)
    product=s.MutableDenseNDimArray.zeros(8,5,2)
    for a,F in enumerate(S):
        for j,Yj in enumerate(Y):
            values=_coefficients(F*Yj)
            for k,value in enumerate(values):product[k,a,j]=value
            conformal=-F*laplacian(Yj)+gradient_dot(F,Yj)/2
            for k,value in enumerate(_coefficients(conformal)):metric[k,a,j]=value
    for a,F in enumerate(S[1:],start=5):
        for j,Yj in enumerate(Y):
            pure_gauge=(gradient_dot(F,laplacian(Yj))-laplacian(gradient_dot(F,Yj)))/2
            for k,value in enumerate(_coefficients(pure_gauge)):metric[k,a,j]=value
    for j,Yj in enumerate(Y):
        vector=_W_derivative(laplacian(Yj))-laplacian(_W_derivative(Yj))
        for k,value in enumerate(_coefficients(vector)):metric[k,9,j]=value
    # For the TT G22 remainder, divergence and trace vanish. Its contraction
    # with Hess(Y1) is zero; its Hess(Y2) contraction gives the stated identity.
    for k,value in enumerate(_coefficients(2*C2*C2*Y[1]/9)):metric[k,10,1]=value
    trace=s.zeros(5,11)
    for a in range(5):trace[a,a]=3
    for a,eigenvalue in enumerate((8,8,15,24),start=1):trace[a,a+4]=-eigenvalue
    labels=MappingProxyType(dict(
        scalar=('1','z^2-1/4','x^2+y^2-1/2','zxy','x^2y^2-(x^2+y^2)/8+1/48'),
        cubic=('z','xy','z^3','z(x^2+y^2)','z^2xy','xy(x^2+y^2)','zx^2y^2','x^3y^3'),
        metric=('gamma*S0','gamma*S2a','gamma*S2b','gamma*S3','gamma*S4',
                'Hess*S2a','Hess*S2b','Hess*S3','Hess*S4','symgradW','Z'),
        linear_field=('sqrt(2)/pi*z','sqrt(12)/pi*xy')))
    return MappingProxyType(dict(laplacian=s.ImmutableMatrix(delta),
        metric_delta_laplacian=s.ImmutableDenseNDimArray(metric),
        multiplication=s.ImmutableDenseNDimArray(product),trace=s.ImmutableMatrix(trace),labels=labels))


@lru_cache(maxsize=1)
def numeric_operator_matrices():
    """Cached read-only float arrays with the exact map's index ordering."""
    exact=exact_operator_matrices();result={}
    for key in ('laplacian','metric_delta_laplacian','multiplication','trace'):
        array=np.array(exact[key].tolist(),dtype=float)
        array.setflags(write=False)
        result[key]=array
    result['labels']=exact['labels']
    return MappingProxyType(result)
