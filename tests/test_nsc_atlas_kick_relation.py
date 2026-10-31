"""Closed geometry factors checked through independent ambient York algebra."""
import math

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_kick_relation import geometry_factor_squared,observer_relation


@pytest.mark.parametrize('source',[1,2])
def test_actual_york_tensor_and_epsilon_cancel(source):
    with threadpool_limits(limits=1):
        rows=[observer_relation(epsilon=epsilon,source_harmonic=source) for epsilon in (.125,.25,.5)]
    reference=rows[1]
    for row in rows:
        assert row['geometry_factor_from_york']==pytest.approx(row['closedfactor'],rel=2e-13)
        assert row['ratio_relative_check_error']<2e-13
        assert row['shear_to_logcoupling_jerk_ratio']==reference['shear_to_logcoupling_jerk_ratio']
        assert row['shear_norm_over_abs_lambda_psi_minus6']/row['epsilon']==pytest.approx(
            reference['shear_norm_over_abs_lambda_psi_minus6']/reference['epsilon'],rel=2e-13)
        assert row['logcoupling_jerk_over_lambda_psi_minus6']/row['epsilon']==pytest.approx(
            reference['logcoupling_jerk_over_lambda_psi_minus6']/reference['epsilon'],rel=2e-13)
    expected={1:6.921168076236156e-5,2:9.469839236302457e-6}
    assert reference['shear_to_logcoupling_jerk_ratio']==pytest.approx(expected[source],rel=2e-14)
    assert not reference['SI_scale_selected'] and not reference['quantum_state_selected']


def test_undefined_zero_preparation_and_invalid_model_inputs_are_rejected():
    for epsilon in (0,-.1,2,float('nan'),True):
        with pytest.raises(ValueError):observer_relation(epsilon=epsilon)
    for z in (0,-1,float('inf'),True):
        with pytest.raises(ValueError):geometry_factor_squared(z,1)
    for source in (0,3,True,1.):
        with pytest.raises(ValueError):geometry_factor_squared(10,source)
    assert math.isfinite(geometry_factor_squared(1e300,1))


def test_compactification_dimension_cancels_from_the_proper_jerk():
    import sympy as s
    d,kappa,mass,chi,Pi=s.symbols('d kappa mass chi Pi',positive=True)
    K=d*(d+2)/(2*kappa)
    third_s=d*mass*chi*Pi/K
    third_logg=-(d+2)*third_s/2
    assert s.simplify(third_logg+kappa*mass*chi*Pi)==0


def test_closed_factors_follow_from_full_ambient_york_tensors():
    import sympy as s
    x,y,z,w,q,v,Z=s.symbols('x y z w q v Z',real=True);X=(x,y,z,w)
    frame=((w,z,-y,-x),(-z,w,x,-y),(y,-x,w,-z))
    def E(i,f):return s.expand(sum(a*s.diff(f,b) for a,b in zip(frame[i],X)))
    def L(W):
        divergence=sum(E(i,W[i]) for i in range(3))
        return s.Matrix([[s.expand(E(i,W[j])+E(j,W[i])-(s.Rational(2,3)*divergence if i==j else 0))
            for j in range(3)] for i in range(3)])
    def Q(W):
        A=L(W)
        return s.Matrix([s.expand(sum(E(j,A[j,i]) for j in range(3))) for i in range(3)])
    chi=q*x+y*w
    for number,f in enumerate((x,y*w),1):
        source=s.Matrix([-f*E(i,chi) for i in range(3)])
        q1=Q(source);q2=Q(q1);q3=Q(q2)
        W=-(3*q3+167*q2+2984*q1+20080*source)/44800
        assert (Q(W)-source).applyfunc(s.simplify)==s.zeros(3,1)
        A=L(W).subs({x:q,y:v,z:0,w:v})
        polynomial=s.Poly(s.expand(sum(component**2 for component in A)),v)
        assert all(powers[0]%2==0 for powers,_ in polynomial.terms())
        norm2=s.factor(sum(c*((1-q*q)/2)**(powers[0]//2) for powers,c in polynomial.terms()))
        f_anchor=q if number==1 else (1-q*q)/2
        factor2=s.factor(norm2/(((1+q*q)/2)**2*f_anchor**2))
        num,den=s.fraction(factor2)
        def even_substitution(expression):
            terms=s.Poly(s.expand(expression),q).terms()
            assert all(powers[0]%2==0 for powers,_ in terms)
            return sum(c*((Z+8)/(6*(Z+3)))**(powers[0]//2) for powers,c in terms)
        factor_z=s.factor(even_substitution(num)/even_substitution(den))
        for value in (40.,100.,8329.074464994037,1e6):
            expected=float(factor_z.subs(Z,value))
            assert geometry_factor_squared(value,number)==pytest.approx(expected,rel=2e-14)


def test_common_factors_and_norm_convention_cancel_without_si_units():
    import sympy as s
    H,mu,lam,psi,kappa,epsilon,b,c,AA,chi,f=s.symbols(
        'H mu lam psi kappa epsilon b c AA chi f',positive=True)
    shear=lam*psi**-6*kappa*epsilon*b*c*AA
    jerk=kappa*mu**2*epsilon*b*chi*lam*c*f*psi**-6
    assert s.simplify(H**2*shear/jerk-H**2*AA/(mu**2*chi*f))==0
