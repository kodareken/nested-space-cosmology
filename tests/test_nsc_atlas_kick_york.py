"""Exact polynomial compatibility and independent York tensor checks."""
import math

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_kick_york import (
    inherited_polynomial_inverse,longitudinal_operator,longitudinal_tensor,
    solve_longitudinal,solve_york_kick,tensor_divergence,
)
from recursive_horizons.nsc_atlas_s3_frame import build_frame_basis


@pytest.mark.parametrize('source',[1,2])
def test_both_inherited_kicks_have_exact_reference_completion(source):
    with threadpool_limits(limits=1):kick=solve_york_kick(source_harmonic=source)
    d=kick.diagnostics
    assert d['kernel_dimension']==10
    assert d['CKV_charge_over_source_norm']<2e-13
    assert d['divergence_residual_coefficient_norm']<1e-15
    assert d['source_projection_L2_residual']<1e-15
    assert d['trace_coefficient_norm']<1e-15
    assert abs(d['integration_by_parts_residual'])<2e-17
    assert d['closed_inverse_relative_difference']<2e-12
    rng=np.random.default_rng(810+source);points=rng.normal(size=(23,4));points/=np.linalg.norm(points,axis=1)[:,None]
    A=kick.values(points);derivative=kick.derivatives(points)
    np.testing.assert_allclose(A,A.swapaxes(1,2),atol=1e-17,rtol=0)
    np.testing.assert_allclose(np.trace(A,axis1=1,axis2=2),0,atol=1e-17,rtol=0)
    div=np.stack([sum(derivative[:,j,j,i] for j in range(3)) for i in range(3)],axis=1)
    np.testing.assert_allclose(div,kick.source_values(points),atol=2e-17,rtol=2e-12)
    np.testing.assert_allclose(kick.squared_norm(points),np.sum(A*A,axis=(1,2)),atol=0,rtol=0)
    profile=kick.momentum_profile(kick.basis.points)
    assert kick.basis.weights@(profile*profile)/kick.parent.H**3==pytest.approx(1.,abs=2e-14)


def test_ordered_frame_operator_and_minimum_norm_identity():
    basis=build_frame_basis(4,7);rng=np.random.default_rng(520)
    W=rng.normal(size=(3,len(basis.labels)));r=.37
    A=longitudinal_tensor(basis,W,r=r);div=tensor_divergence(basis,A,r=r)
    operator=longitudinal_operator(basis,r=r)
    np.testing.assert_allclose(div.ravel(),operator@W.ravel(),atol=2e-14,rtol=2e-14)
    assert np.sum(W*div)==pytest.approx(-np.sum(A*A)/2,rel=2e-14)
    solution=solve_longitudinal(basis,div,r=r)
    np.testing.assert_allclose(solution['tensor_coefficients'],A,atol=3e-14,rtol=2e-13)
    # Gauge changes W by a CKV without changing its physical tensor.
    vals,vecs=np.linalg.eigh(operator);ckv=vecs[:,np.argmin(abs(vals))].reshape(W.shape)
    np.testing.assert_allclose(longitudinal_tensor(basis,W+ckv,r=r),A,atol=2e-14,rtol=2e-14)


def test_incompatible_killing_source_is_rejected_and_scale_is_physical():
    basis=build_frame_basis(4,7);source=np.zeros((3,len(basis.labels)));source[0,0]=1
    with pytest.raises(ValueError,match='conformal-Killing'):
        solve_longitudinal(basis,source,r=.2)
    np.testing.assert_allclose(longitudinal_operator(basis,r=.2),.04*longitudinal_operator(basis,r=1),atol=2e-16,rtol=2e-14)
    with pytest.raises(ValueError,match='outside'):
        inherited_polynomial_inverse(basis,source)


def sphere_mean_polynomial(expression,coordinates):
    import sympy as s
    answer=0
    for powers,coefficient in s.Poly(s.expand(expression),*coordinates).terms():
        if any(power%2 for power in powers):continue
        halves=[power//2 for power in powers];total=sum(halves)
        numerator=s.prod(s.factorial2(2*power-1) for power in halves if power)
        denominator=s.prod(4+2*k for k in range(total))
        answer+=coefficient*numerator/denominator
    return s.simplify(answer)


def test_all_ten_conformal_killing_charges_vanish_exactly():
    import sympy as s
    X=s.symbols('x y z w');a,b=s.symbols('a b');x,y,z,w=X
    chi=a*x+b*y*w;radial=sum(q*s.diff(chi,q) for q in X)
    actions=[X[i]*s.diff(chi,X[j])-X[j]*s.diff(chi,X[i]) for i in range(4) for j in range(i+1,4)]
    actions += [s.diff(chi,q)-q*radial for q in X]
    assert len(actions)==10
    for f in (x,y*w):
        assert [sphere_mean_polynomial(f*action,X) for action in actions]==[0]*10


def test_closed_inverse_solves_both_ambient_polynomial_sources_exactly():
    import sympy as s
    x,y,z,w=s.symbols('x y z w');X=(x,y,z,w);a,b=s.symbols('a b')
    directions=((w,z,-y,-x),(-z,w,x,-y),(y,-x,w,-z))
    def E(i,f):return s.expand(sum(v*s.diff(f,q) for v,q in zip(directions[i],X)))
    def Q(W):
        divergence=sum(E(i,W[i]) for i in range(3))
        L=[[E(i,W[j])+E(j,W[i])-(s.Rational(2,3)*divergence if i==j else 0)
            for j in range(3)] for i in range(3)]
        return s.Matrix([s.expand(sum(E(j,L[j][i]) for j in range(3))) for i in range(3)])
    chi=a*x+b*y*w
    for f in (x,y*w):
        source=s.Matrix([-f*E(i,chi) for i in range(3)])
        q1=Q(source);q2=Q(q1);q3=Q(q2)
        W=-(3*q3+167*q2+2984*q1+20080*source)/44800
        assert (Q(W)-source).applyfunc(s.simplify)==s.zeros(3,1)
        assert (Q(q3)+s.Rational(167,3)*q3+s.Rational(2984,3)*q2+s.Rational(20080,3)*q1+s.Rational(44800,3)*source).applyfunc(s.simplify)==s.zeros(3,1)
