"""Independent canonical, noncommuting Magnus and bosonic-entropy controls."""
import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag

from recursive_horizons.nsc_atlas_boson_transport import (
    canonical_form,cf4_step,covariance_diagnostics,covariance_inheritance,
    gaussian_entropy,hamiltonian_generator,propagate,reduced_covariance,
    relative_transport,retarded_between,retarded_qp,symplectic_diagnostics,
    symplectic_eigenvalues,validate_covariance,
)


def oscillator(omega,t):
    c=math.cos(omega*t);s=math.sin(omega*t)
    return np.array([[c,s/omega],[-omega*s,c]])


def test_cubic_cross_response_from_time_dependent_canonical_equations():
    import sympy as s
    J=s.Matrix([[0,0,1,0],[0,0,0,1],[-1,0,0,0],[0,-1,0,0]])
    coefficients=[];generators=[]
    for order in range(3):
        a,b,ks,kc,mix,bs,bc=s.symbols('a%d b%d ks%d kc%d mix%d bs%d bc%d'%((order,)*7))
        coefficients.append((a,b,mix))
        G=s.Matrix([[ks,mix,bs,0],[mix,kc,0,bc],[bs,0,a,0],[0,bc,0,b]])
        generators.append(J*G)
    D,D1,D2=generators;v=s.Matrix([0,0,0,1])
    first=D*v;second=(D1+D*D)*v
    third=(D2+2*D1*D+D*D1+D**3)*v
    a,b,mix=coefficients[0]
    assert s.simplify(first[0])==0 and s.simplify(second[0])==0
    assert s.simplify(third[0]+a*mix*b)==0


def test_constant_oscillator_exact_map_and_retarded_normalization():
    omega=3.7
    result=propagate(lambda t:np.diag([omega**2,1.]),[0.,.37,1.3],max_step=.23,checkpoints=[.4,.9])
    np.testing.assert_array_equal(result.times,[0.,.4,.9,1.3])
    for t,M in zip(result.times,result.matrices):
        np.testing.assert_allclose(M,oscillator(omega,t),rtol=2e-14,atol=2e-14)
    assert result.max_checkpoint_symplectic_defect<2e-14
    assert result.hamiltonian_evaluations==2*result.steps+1
    assert retarded_qp(result.matrix)[0,0]==pytest.approx(math.sin(omega*1.3)/omega,abs=1e-14)
    interval=relative_transport(result.matrices[2],result.matrices[1])
    np.testing.assert_allclose(interval,oscillator(omega,.5),atol=2e-14)
    assert retarded_between(result.matrices[2],result.matrices[1])[0,0]==pytest.approx(math.sin(omega*.5)/omega,abs=2e-14)


def varying_hamiltonian(t):
    K=np.array([[2+.3*math.sin(t),.2*math.cos(.7*t)],
                [.2*math.cos(.7*t),5+.4*math.cos(t)]])
    A=np.array([[1.,.04*math.sin(.3*t)],[.04*math.sin(.3*t),1.2]])
    B=np.array([[0.,.12*math.sin(1.3*t)],[-.07*math.cos(t),.02*t]])
    return np.block([[K,B.T],[B,A]])


def test_noncommuting_time_dependence_has_fourth_order_and_correct_ordering():
    J=canonical_form(2)
    early=J@varying_hamiltonian(.2);late=J@varying_hamiltonian(1.4)
    assert np.linalg.norm(early@late-late@early)>.1
    independent=solve_ivp(lambda t,y:(J@varying_hamiltonian(t)@y.reshape(4,4)).ravel(),
        (0.,2.3),np.eye(4).ravel(),method='DOP853',rtol=2e-13,atol=2e-15)
    assert independent.success
    reference=independent.y[:,-1].reshape(4,4)
    errors=[]
    for step in (.4,.2,.1):
        result=propagate(varying_hamiltonian,[0.,2.3],max_step=step)
        errors.append(np.linalg.norm(result.matrix-reference))
        assert symplectic_diagnostics(result.matrix)['max_abs']<2e-13
    assert errors[0]/errors[1]>10 and errors[1]/errors[2]>12
    assert errors[2]<2e-6


def test_backward_composition_and_exact_breakpoint_domain():
    evaluated=[]
    def G(t):
        evaluated.append(t)
        assert 0<=t<=1.1
        return varying_hamiltonian(t)
    forward=propagate(G,[0.,.3,.7,1.1],max_step=.04,checkpoints=[.55])
    backward=propagate(G,[1.1,.7,.3,0.],max_step=.04,checkpoints=[.55])
    np.testing.assert_allclose(backward.matrix@forward.matrix,np.eye(4),atol=2e-12)
    np.testing.assert_allclose(forward.times,[0.,.55,1.1])
    assert evaluated


def test_supplied_thermal_covariance_entropy_is_inherited_globally():
    omega=np.array([1.3,2.7]);nu=np.array([.8,1.2])
    covariance=np.diag(np.r_[nu/omega,nu*omega])
    M=propagate(varying_hamiltonian,[0.,2.3],max_step=.08).matrix
    inherited=covariance_inheritance(M,covariance)
    np.testing.assert_allclose(symplectic_eigenvalues(covariance),nu,atol=2e-15)
    np.testing.assert_allclose(symplectic_eigenvalues(inherited),nu,atol=1e-13)
    expected=sum((n+1)*math.log(n+1)-n*math.log(n) for n in nu-.5)
    assert gaussian_entropy(covariance)==pytest.approx(expected,abs=2e-14)
    assert gaussian_entropy(inherited)==pytest.approx(expected,abs=3e-13)
    assert covariance_diagnostics(inherited)['physical']


def test_two_mode_squeezing_changes_reduced_entropy_without_global_production():
    r=.63;c=math.cosh(r);s=math.sinh(r)
    M=block_diag(np.array([[c,s],[s,c]]),np.array([[c,-s],[-s,c]]))
    covariance=.5*np.eye(4)
    inherited=covariance_inheritance(M,covariance)
    assert symplectic_diagnostics(M)['max_abs']<1e-15
    assert gaussian_entropy(inherited)<1e-12
    reduced=reduced_covariance(inherited,[0])
    np.testing.assert_allclose(reduced,.5*math.cosh(2*r)*np.eye(2),atol=2e-15)
    occupation=math.sinh(r)**2
    expected=(occupation+1)*math.log(occupation+1)-occupation*math.log(occupation)
    assert gaussian_entropy(reduced)==pytest.approx(expected,abs=2e-14)
    returned=covariance_inheritance(np.linalg.inv(M),inherited)
    assert gaussian_entropy(reduced_covariance(returned,[0]))<1e-12


def test_bosonic_uncertainty_has_no_fermionic_upper_bound_and_tracks_hbar():
    validate_covariance(100*np.eye(2))
    with pytest.raises(ValueError,match='CCR'):validate_covariance(.49*np.eye(2))
    assert gaussian_entropy(.5*3.2*np.eye(2),hbar=3.2)<1e-12
    assert gaussian_entropy((.5+1e-12)*np.eye(2))>2e-11
    with pytest.raises(ValueError,match='CCR'):validate_covariance(np.diag([1e12,1e-14]))
    assert not covariance_diagnostics(np.diag([-1.,2.]))['physical']
    with pytest.raises(ValueError,match='positive definite'):symplectic_eigenvalues(np.diag([-1.,2.]))


def test_invalid_hamiltonians_intervals_and_selections_are_rejected():
    with pytest.raises(ValueError,match='symmetric'):hamiltonian_generator([[1,2],[0,1]])
    with pytest.raises(ValueError,match='real square'):hamiltonian_generator(np.eye(3))
    with pytest.raises(ValueError,match='real square'):validate_covariance(.5*np.eye(2,dtype=complex))
    with pytest.raises(ValueError,match='strictly monotone'):propagate(varying_hamiltonian,[0,1,.5],max_step=.1)
    with pytest.raises(ValueError,match='outside'):propagate(varying_hamiltonian,[0,1],max_step=.1,checkpoints=[1.1])
    with pytest.raises(ValueError,match='nonzero'):cf4_step(varying_hamiltonian,0,0)
    with pytest.raises(ValueError,match='distinct'):reduced_covariance(np.eye(4),[0,0])
    with pytest.raises(ValueError,match='dimensions differ'):covariance_inheritance(np.eye(2),np.eye(4))
