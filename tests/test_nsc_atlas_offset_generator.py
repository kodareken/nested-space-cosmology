import math
import mpmath as mp
import numpy as np
import pytest
import sympy as s
from recursive_horizons.nsc_atlas_offset_generator import (
    offset_value,phase_history,directions,euclidean_kernel_distance,closed_parent_point,sphere_cover)


def test_one_offset_iteration_matches_the_accumulated_closed_form():
    with mp.workdps(80):
        a=offset_value('golden');theta=phi=mp.mpf(0)
        for n in range(400):
            assert abs(mp.sin(mp.pi*(theta-n*a)))<mp.mpf('1e-74')
            assert abs(mp.sin(mp.pi*(phi-n*(n-1)*a/2)))<mp.mpf('1e-73')
            oldtheta=theta;theta,phi=mp.frac(theta+a),mp.frac(phi+theta)
            recovered=mp.frac(theta-a)
            assert abs(mp.sin(mp.pi*(recovered-oldtheta)))<mp.mpf('1e-75')


def test_polynomial_character_and_its_difference_have_the_claimed_orders():
    n,h,a,k,l=s.symbols('n h a k l')
    P=a*(k*n+l*n*(n-1)/2)
    assert s.Poly(P,n).coeff_monomial(n**2)==a*l/2
    assert s.expand(P.subs(n,n+h)-P).coeff(n)==a*l*h


def test_distance_kernel_is_a_gram_metric_not_an_extra_axis_assumption():
    p=phase_history(12,'action');v=directions(p)
    np.testing.assert_allclose(np.sum(v*v,axis=1),1,atol=5e-16)
    G=v@v.T
    assert np.linalg.matrix_rank(G,tol=1e-11)==3
    for i,j in ((1,5),(4,8),(0,9)):
        d=euclidean_kernel_distance(.4,.8,*(2*np.pi*p[i]),*(2*np.pi*p[j]))
        assert d==pytest.approx(np.linalg.norm(.4*v[i]-.8*v[j]),abs=4e-16)


def test_closed_parent_metric_and_volume_from_the_embedding():
    r,R,t,p=s.symbols('r R t p',positive=True)
    x=s.Matrix([R*s.cos(r/R),R*s.sin(r/R)*s.sin(t)*s.cos(p),
                R*s.sin(r/R)*s.sin(t)*s.sin(p),R*s.sin(r/R)*s.cos(t)])
    J=x.jacobian([r,t,p]);g=(J.T*J).applyfunc(s.trigsimp)
    assert g==s.diag(1,R**2*s.sin(r/R)**2,R**2*s.sin(r/R)**2*s.sin(t)**2)
    assert s.trigsimp(x.dot(x))==R**2
    np.testing.assert_allclose(closed_parent_point(0,[0,0,1]),[1,0,0,0])


def test_finite_line_cover_controls_ball_distance():
    # Octahedral endpoints: maximum direction angle acos(1/sqrt3).
    v=np.r_[np.eye(3),-np.eye(3)];h=sphere_cover(v)
    assert h==pytest.approx(math.acos(1/math.sqrt(3)))
    assert math.sin(h)==pytest.approx(math.sqrt(2/3))


def test_phase_sampling_and_metric_area_are_different_measures():
    t=s.symbols('t',real=True)
    assert s.integrate(s.cos(t)**2,(t,0,2*s.pi))/(2*s.pi)==s.Rational(1,2)
    assert s.integrate(s.sin(t)*s.cos(t)**2,(t,0,s.pi))/2==s.Rational(1,3)
    # A one-pass line retains zero transverse rank; the completion is essential.
    assert np.linalg.matrix_rank(directions(phase_history(1)))==1


def test_rational_offset_returns_a_finite_history_instead_of_dense_completion():
    from fractions import Fraction
    a=Fraction(1,4);theta=phi=Fraction(0);states=[]
    for _ in range(8):
        states.append((theta,phi));theta,phi=(theta+a)%1,(phi+theta)%1
    assert (theta,phi)==states[0]
    assert len(set(states))==8


def test_closed_meridian_cover_has_the_claimed_radial_dependence():
    direction=np.ones(3)/math.sqrt(3);cover=math.acos(1/math.sqrt(3))
    for radial in (.2,math.pi/2,2.8):
        target=closed_parent_point(radial,direction)
        best=math.atan2(math.sin(radial)/math.sqrt(3),math.cos(radial))
        nearest=closed_parent_point(best,[1,0,0])
        distance=math.acos(np.clip(target@nearest,-1,1))
        assert distance<=cover+2e-15
        if radial==math.pi/2:assert distance==pytest.approx(cover)


def test_parent_laplacian_modes_follow_the_generated_metric():
    r,R,t,p=s.symbols('r R t p',positive=True)
    def lap(f):
        return s.diff(f,r,2)+2*s.cot(r/R)*s.diff(f,r)/R+(s.diff(f,t,2)+s.cot(t)*s.diff(f,t)+s.diff(f,p,2)/s.sin(t)**2)/(R**2*s.sin(r/R)**2)
    first=s.cos(r/R)
    second=s.cos(r/R)*s.sin(r/R)*s.cos(t)
    assert s.trigsimp(lap(first)+3*first/R**2)==0
    assert s.trigsimp(lap(second)+8*second/R**2)==0


def test_first_three_irrational_offset_directions_have_rank_three():
    a=s.symbols('a',real=True)
    v0=s.Matrix([0,0,1]);v1=s.Matrix([s.sin(a),0,s.cos(a)])
    v2=s.Matrix([s.sin(2*a)*s.cos(a),s.sin(2*a)*s.sin(a),s.cos(2*a)])
    assert s.simplify(s.Matrix.hstack(v0,v1,v2).det()-s.sin(a)**2*s.sin(2*a))==0
