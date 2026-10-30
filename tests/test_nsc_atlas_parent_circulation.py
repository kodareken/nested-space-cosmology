"""Independent spherical identities, constraints and observer controls."""
import math
from pathlib import Path
import json

import numpy as np
import pytest
import sympy as s

from recursive_horizons.nsc_atlas_parent_circulation import (
    Parent, VOLUME, HARMONIC_PRODUCT, evolve, energy, initial_constraint,
    differential_coefficient, direct_light_control, phase_speed,
    horizon_flux,static_patch_energy)

X = s.Matrix(s.symbols('x1:5', real=True))
P = s.eye(4)-X*X.T


def grad(f):
    return (P*s.Matrix([s.diff(f,x) for x in X])).applyfunc(s.expand)


def lap(f):
    return s.expand(sum(P[i,j]*s.diff(f,X[i],X[j]) for i in range(4) for j in range(4))
                    -3*sum(X[i]*s.diff(f,X[i]) for i in range(4)))


def on_sphere(f):
    return s.rem(s.Poly(s.expand(f),*X),s.Poly(X.dot(X)-1,*X)).as_expr()


def mean(f):
    result=0
    for powers,coefficient in s.Poly(s.expand(f),*X).terms():
        if any(n%2 for n in powers):continue
        half=[n//2 for n in powers];N=sum(half)
        moment=s.prod(s.factorial2(2*n-1) for n in half)/s.prod(4+2*i for i in range(N))
        result+=coefficient*moment
    return s.simplify(result)


def harmonics():
    return s.sqrt(2)*X[0]/s.pi,s.sqrt(12)*X[1]*X[3]/s.pi


def vector():
    C=2*s.sqrt(6)/s.pi**2
    return C*s.Matrix([-4*X[1]*X[3],2*X[0]*X[3],0,2*X[0]*X[1]])/3


def test_drawing_update_has_a_nonzero_metric_defect():
    a=s.symbols('a',real=True)
    v0=s.Matrix([0,0,1]);v1=s.Matrix([s.sin(a),0,s.cos(a)])
    v2=s.Matrix([s.sin(2*a)*s.cos(a),s.sin(2*a)*s.sin(a),s.cos(2*a)])
    defect=(v2-v1).dot(v2-v1)-(v1-v0).dot(v1-v0)
    assert s.trigsimp(s.expand_trig(defect-4*s.sin(a)**2*s.cos(a)*(1-s.cos(a))))==0
    p=Parent.retained();assert float(defect.subs(a,2*math.pi*p.offset))>0


def test_scalar_harmonics_are_normalized_and_have_distinct_eigenvalues():
    y1,y2=harmonics()
    for y,l in ((y1,3),(y2,8),(y1*y2,15)):
        assert lap(y)+l*y==0
    assert mean(y1*y1)*2*s.pi**2==1
    assert mean(y2*y2)*2*s.pi**2==1
    assert mean(y1*y2)==0


def test_hodge_projection_vector_laplacian_and_global_rotation_charges():
    y1,y2=harmonics();W=vector()
    target=y1*grad(y2)-y2*grad(y1)-grad(y1*y2)/3
    assert all(s.simplify(z)==0 for z in W-target)
    assert s.expand(X.dot(W))==0
    assert on_sphere(s.trace(P*W.jacobian(X)))==0
    # For an ambient tangent field on unit S3: roughDelta W=P*Delta W+W.
    rough=P*s.Matrix([lap(w) for w in W])+W
    assert all(on_sphere(z)==0 for z in rough+7*W)
    assert mean(W.dot(W))*2*s.pi**2==16/(3*s.pi**2)
    for i in range(4):
        for j in range(i+1,4):
            assert mean(X[i]*W[j]-X[j]*W[i])==0


def test_transverse_spatial_stress_drives_current_and_matches_constraint():
    y1,y2=harmonics();W=vector();g1,g2=grad(y1),grad(y2)
    # div(dY1 dY2+dY2 dY1)=-lambda1 Y1 dY2-lambda2 Y2 dY1+d(gradY1.gradY2)
    div=-3*y1*g2-8*y2*g1+grad(g1.dot(g2))
    rem=div-s.Rational(5,2)*W
    # Remaining part is a pure gradient; d(rem)=0 on tangent planes.
    curl=P*(rem.jacobian(X)-rem.jacobian(X).T)*P
    assert all(on_sphere(z)==0 for z in curl)
    J,beta,kappa=s.symbols('J beta kappa')
    assert s.solve(-s.Rational(5,2)*beta+kappa*J/2,beta)==[kappa*J/5]


def test_exact_constraint_linearization_and_polynomial_corrector():
    a0,mu2,k,q1,q2=s.symbols('a0 mu2 k q1 q2')
    y1,y2=harmonics();f=q1*y1+q2*y2
    A=mu2*a0*a0;B1=q1*q1*2/s.pi**2;B2=q2*q2*12/s.pi**2;B12=4*s.sqrt(6)*q1*q2/s.pi**2
    c1=k*B1*(A-1)/40;c3=k*B2*(A-4)/168
    c2=k*B2/40+2*c3/5;c4=k*B12*(A-2)/96;c0=(-k*B1/8-2*c1-4*c2)/3
    u=c0+c1*X[0]**2+c2*(X[1]**2+X[3]**2)+c3*X[1]**2*X[3]**2+c4*X[0]*X[1]*X[3]
    source=grad(f).dot(grad(f))+A*f*f
    assert s.simplify(on_sphere(lap(u)+3*u+k*source/8))==0
    eps=s.symbols('eps');psi=1+eps**2*u
    constraint=-8*lap(psi)+6*psi-6*psi**5-k*eps**2*grad(f).dot(grad(f))*psi-k*A*eps**2*f*f*psi**5
    assert s.simplify(on_sphere(s.expand(constraint).coeff(eps,2)))==0
    assert initial_constraint()['mean_u']==pytest.approx(-.25,abs=1e-14)


def test_symmetries_remove_the_whole_l1_kernel():
    y1,y2=harmonics();q1,q2=s.symbols('q1 q2');f=q1*y1+q2*y2
    symmetries=({X[0]:-X[0],X[1]:-X[1]}, {X[0]:-X[0],X[3]:-X[3]}, {X[2]:-X[2]})
    for flip in symmetries:
        assert s.expand((f*f).subs(flip,simultaneous=True)-f*f)==0
    for x in X:
        assert lap(x)+3*x==0
        assert any(s.expand(x.subs(flip,simultaneous=True)+x)==0 for flip in symmetries)


def test_phase_readout_agrees_with_independent_proper_evolution_and_work():
    p=Parent.retained();phase=evolve(parent=p,passes=40)
    endpoint=phase.y[8,-1];proper=evolve(parent=p,clock='proper',endpoint=endpoint)
    np.testing.assert_allclose(phase.y[:,-1],proper.y[:,-1],atol=2e-10,rtol=2e-9)
    z=proper.y[:,-1]
    assert z[6]==pytest.approx(z[3]*z[4]-z[2]*z[5],rel=1e-9)
    assert energy(z,p)-energy(proper.y[:,0],p)==pytest.approx(z[9],abs=2e-10)
    assert phase_speed(z,p)>0


def test_same_eigenvalue_preparation_has_exactly_zero_circulation():
    p=Parent.retained();sol=evolve(parent=p,clock='proper',endpoint=.5,eigenvalues=(3.,3.))
    assert np.max(abs(sol.y[6]))==0
    np.testing.assert_allclose(sol.y[2],sol.y[4],atol=0,rtol=0)


def test_ring_projection_and_phase_average():
    phi,R,b,C=s.symbols('phi R b C',real=True)
    W=vector().subs({X[0]:R*s.cos(phi),X[1]:R*s.sin(phi),X[3]:b})
    tangent=s.Matrix([-R*s.sin(phi),R*s.cos(phi),0,0])
    assert s.trigsimp(W.dot(tangent)-2*s.sqrt(6)/s.pi**2*b*R**2*(1-s.cos(2*phi)/3))==0
    assert s.integrate(s.cos(2*phi),(phi,0,2*s.pi))==0


def test_static_killing_field_and_smooth_horizon_generator():
    t,c,H=s.symbols('t c H',positive=True)
    a=s.cosh(H*t)/H
    g=s.diag(-1,a*a)
    K=s.Matrix([s.cos(c),-H*s.tanh(H*t)*s.sin(c)])
    coords=[t,c]
    for i in range(2):
        for j in range(2):
            lie=sum(K[k]*s.diff(g[i,j],coords[k])+g[k,j]*s.diff(K[k],coords[i])+g[i,k]*s.diff(K[k],coords[j]) for k in range(2))
            assert s.trigsimp(lie)==0
    # Angular metric component is also invariant under this K.
    angular=a*a*s.sin(c)**2
    assert s.trigsimp(sum(K[i]*s.diff(angular,coords[i]) for i in range(2)))==0
    ch=s.pi/2-s.atan(s.sinh(H*t))
    assert s.simplify(s.diff(ch,t)+1/a)==0
    assert s.simplify(a*s.sin(ch)-1/H)==0


def test_horizon_flux_equals_independent_volume_energy_depletion():
    from scipy.integrate import quad
    p=Parent.retained();sol=evolve(parent=p,clock='proper',endpoint=.75)
    emitted,error=quad(lambda t:horizon_flux(t,sol.sol(t),p)[1],0.,.75,epsabs=1e-11,epsrel=1e-11)
    E0=static_patch_energy(0.,sol.y[:,0],p)
    E1=static_patch_energy(.75,sol.y[:,-1],p)
    assert E0-E1==pytest.approx(emitted,abs=2e-10)
    assert emitted>0
    for t in np.linspace(0.,.75,31):
        null,power=horizon_flux(t,sol.sol(t),p)
        assert null>=0 and power>=0
    assert error<1e-9


def test_finite_light_signal_and_common_rotation_control():
    p=Parent.retained();sol=evolve(parent=p,clock='proper',endpoint=3.4)
    leading=differential_coefficient(p,sol)
    assert leading['common_comparison_conformal_margin']>0
    assert leading['H_delay_over_epsilon_squared']==pytest.approx(1.78600055926e-5,rel=2e-8)
    for frame in (0.,.03):
        direct=direct_light_control(p,sol,epsilon=.1,frame_amplitude=frame)
        assert abs(direct['relative_difference_from_quadratic'])<2e-6


def test_recorded_sources_and_numerical_controls():
    import hashlib
    root=Path(__file__).resolve().parents[1]
    path=root/'results/summaries/nsc-atlas-parent-circulation-v1.json'
    if not path.exists():pytest.skip('record is produced after source freeze')
    record=json.loads(path.read_text(encoding='utf-8'))
    for name,expected in record['source_sha256'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==expected
    c=record['controls']
    assert c['phase_proper_max_state_difference']<1e-8
    assert c['closed_parent_constraint_max_residual']<1e-11
    assert c['wronskian_relative_residual']<1e-8
    assert c['energy_work_relative_residual']<1e-9
