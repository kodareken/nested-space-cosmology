"""Independent finite-CAR/Clifford controls, without a global-state claim."""
import itertools
import math

import numpy as np
import pytest
import sympy as s

from recursive_horizons.nsc_atlas_valve_source_gate import (
    scalar_coherent_gaussian_control, signed_gaussian_coincidence,
)


def _zero(matrix):return all(s.simplify(entry)==0 for entry in matrix)


def _product(matrices,dimension):
    result=s.eye(dimension)
    for matrix in matrices:result=result*matrix
    return result


def _lorentz_clifford(D):
    """Build raw Clifford matrices and align only beta/chirality to the API."""
    k=D//2;ns=2**k;X=s.Matrix([[0,1],[1,0]])
    Y=s.Matrix([[0,-s.I],[s.I,0]]);Z=s.diag(1,-1);I=s.eye(2)
    euclidean=[]
    for j in range(k):
        euclidean.extend(s.kronecker_product(*([Z]*j+[P]+[I]*(k-j-1))) for P in (X,Y))
    raw=[euclidean[0]]+[s.I*G for G in euclidean[1:]]
    grading=s.kronecker_product(*([Z]*k))
    symmetric=[i for i,G in enumerate(raw) if G.T==G]
    antisymmetric=[i for i,G in enumerate(raw) if G.T==-G]
    selected=symmetric if len(symmetric)%2==0 else antisymmetric
    charge_matrix=_product([raw[i] for i in selected],ns)*raw[0].T
    hadamard=(X+Z)/s.sqrt(2)
    U1=s.kronecker_product(hadamard,s.eye(ns//2))
    tail=s.kronecker_product(*([Z]*(k-1)))
    U2=s.diag(s.eye(ns//2),tail);U=U2*U1
    gamma=[s.simplify(U*G*U.H) for G in raw]
    chirality=s.simplify(U*grading*U.H)
    charge=s.simplify(U*charge_matrix*U.T)
    assert _zero(gamma[0]-s.diag(s.eye(ns//2),-s.eye(ns//2)))
    assert _zero(chirality-s.kronecker_product(X,s.eye(ns//2)))
    eta=[-1]+[1]*(D-1)
    for a,b in itertools.product(range(D),repeat=2):
        assert _zero(gamma[a]*gamma[b]+gamma[b]*gamma[a]+(2*eta[a]*s.eye(ns) if a==b else s.zeros(ns)))
    return gamma,chirality,charge,eta


@pytest.mark.parametrize('D',[4,6])
def test_scalar_gaussian_wick_against_explicit_low_dimensional_clifford(D):
    gamma,grading,charge,eta=_lorentz_clifford(D);ns=gamma[0].rows
    c=s.Rational(1,10);density=s.Rational(2);identity=s.eye(ns)
    reference=identity/2;K=reference+c*gamma[0]
    normal=(K-reference)*gamma[0]*density
    assert _zero(normal-c*density*identity)
    assert set(K.eigenvals())=={s.Rational(2,5),s.Rational(3,5)}
    assert all(multiplicity==ns//2 for multiplicity in K.eigenvals().values())
    assert _zero(charge*K.conjugate()*charge.H-(identity-K))
    # External parity reverses three spatial axes and keeps internal axes.
    parity=grading*_product(gamma[1:4],ns)
    assert _zero(parity*parity.H-identity)
    assert _zero(parity*K*parity.H-K)
    for a,b in itertools.combinations(range(1,D),2):
        rotation=gamma[a]*gamma[b]
        assert _zero(rotation*K-K*rotation)
    xi=K*grading-grading*K
    assert not _zero(xi)
    hartree=fock=0
    for indices in itertools.combinations(range(D),3):
        G=_product([gamma[i] for i in indices],ns)
        sign=math.prod(eta[i] for i in indices)
        hartree+=6*sign*s.trace(G*normal)**2
        fock-=6*sign*s.trace(G*normal*G*normal)
    assert hartree==0
    actual=scalar_coherent_gaussian_control(float(c),D=D,density=float(density))
    np.testing.assert_allclose(actual['covariance'],np.array(K,dtype=complex),atol=1e-15)
    np.testing.assert_allclose(actual['reference_covariance'],np.array(reference,dtype=complex),atol=1e-15)
    np.testing.assert_allclose(actual['normal_covariance'],np.array(normal,dtype=complex),atol=1e-15)
    np.testing.assert_allclose(actual['chiral_commutator'],np.array(xi,dtype=complex),atol=1e-15)
    assert actual['literal_Wick_three_form']==pytest.approx(float(fock),rel=2e-15)
    assert actual['contact_lagrangian_over_kappa']==pytest.approx(float(-fock/32),rel=2e-15)
    if D==4:
        # Minimal four-dimensional axial normalization, not a copied D formula.
        axial=0
        for a in range(4):
            A=gamma[a]*grading
            axial-=eta[a]*s.trace(A*normal*A*normal)
        assert s.simplify(fock-6*axial)==0
        assert actual['contact_lagrangian_over_kappa']==pytest.approx(float(-3*axial/16),rel=2e-15)


@pytest.mark.parametrize('c',[-.5,-.1,0.,.1,.5])
def test_actual_16d_local_covariance_is_admissible_with_explicit_scope(c):
    result=scalar_coherent_gaussian_control(c,D=16)
    K=result['covariance'];reference=result['reference_covariance'];ns=256
    np.testing.assert_allclose(K,K.conj().T,atol=0)
    eigenvalues=np.linalg.eigvalsh(K)
    assert np.min(eigenvalues)>=0 and np.max(eigenvalues)<=1
    np.testing.assert_allclose(np.sort(result['covariance_eigenvalues']),eigenvalues,atol=1e-15)
    np.testing.assert_allclose(reference,.5*np.eye(ns),atol=0)
    np.testing.assert_allclose((K-reference)@result['gamma0'],result['normal_covariance'],atol=3e-17)
    if abs(c)==.5:np.testing.assert_allclose(K@K,K,atol=0)
    if c==0:
        assert result['chiral_commutator_norm']==result['literal_Wick_three_form']==0
    else:
        assert result['chiral_commutator_norm']>0
        assert result['literal_Wick_three_form']<0
        assert result['contact_lagrangian_over_kappa']>0
    assert not result['reference_is_vacuum']
    assert not result['global_parent_state_selected']
    assert not result['UV_completion_selected']
    assert not result['radial_NEC_computed']


def test_paired_cancellation_and_its_preservation_by_a_graded_unitary():
    tau1=s.Matrix([[0,1],[1,0]]);tau2=s.Matrix([[0,-s.I],[s.I,0]]);tau3=s.diag(1,-1)
    C0=tau3/4;U=(3*s.eye(2)-4*s.I*tau1)/5
    assert _zero(U.H*U-s.eye(2))
    C=s.simplify(U*C0*U.H)
    # This unitary acts equally on both internal branches and commutes
    # with their physical-mass chirality generator.
    graded_U=s.kronecker_product(U,s.eye(2),s.eye(2))
    grading=s.kronecker_product(s.eye(4),tau2)
    assert _zero(graded_U*grading-grading*graded_U)
    density=s.Rational(3,2);beta=s.kronecker_product(tau3,s.eye(2))
    first=s.kronecker_product(C,s.eye(2))
    opposite=beta*first*beta
    external=s.simplify(density*(first+opposite)*beta)
    gamma5=s.kronecker_product(tau1,s.eye(2));complete=s.zeros(4)
    for sign in (1,-1):
        R=(s.eye(4)-s.I*sign*gamma5)/s.sqrt(2)
        complete+=R*external*R
    assert _zero(s.simplify(complete)) and not _zero(external)
    result=signed_gaussian_coincidence(np.array(C,dtype=complex),density=float(density))
    np.testing.assert_allclose(result['external_normal_covariance'],np.array(external,dtype=complex),atol=3e-16)
    np.testing.assert_allclose(result['complete_normal_covariance'],0,atol=3e-16)
    assert result['external_norm']>0 and result['complete_norm']<3e-16
    # The local scalar control is outside the identical internal-sign family:
    # C/P-compatible local matrices need not have the same cancellation.
    coherent=scalar_coherent_gaussian_control(.1,D=4)
    assert np.linalg.norm(coherent['normal_covariance'])>0


def test_density_normalization_is_explicit_in_both_controls():
    C=np.array([[.2,.1+.03j],[.1-.03j,-.2]])
    one=signed_gaussian_coincidence(C);two=signed_gaussian_coincidence(C,density=2)
    np.testing.assert_allclose(two['external_normal_covariance'],2*one['external_normal_covariance'],atol=0)
    local=scalar_coherent_gaussian_control(.1,D=4)
    scaled=scalar_coherent_gaussian_control(.1,D=4,density=2)
    np.testing.assert_allclose(scaled['normal_covariance'],2*local['normal_covariance'],atol=0)
    assert scaled['literal_Wick_three_form']==4*local['literal_Wick_three_form']


@pytest.mark.parametrize('C',[
    np.eye(2),np.array([[0,1],[0,0]]),np.zeros((3,3)),np.array([[np.nan,0],[0,0]]),
])
def test_signed_control_rejects_noncanonical_covariance_difference(C):
    with pytest.raises(ValueError):signed_gaussian_coincidence(C)


@pytest.mark.parametrize('arguments',[
    {'c':.5001},{'c':float('nan')},{'D':3},{'D':18},{'D':True},{'density':-1},
])
def test_local_gaussian_control_rejects_invalid_cauchy_matrix_data(arguments):
    with pytest.raises(ValueError):scalar_coherent_gaussian_control(**arguments)
