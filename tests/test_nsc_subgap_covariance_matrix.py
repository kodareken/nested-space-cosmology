"""Full Pauli information, source coherence and analytic dual identities."""
import numpy as np
import pytest
from recursive_horizons.nsc_subgap_covariance_matrix import (
    PAULI,pauli_bilinears,from_pauli_traces,real_closed_source_matrix,
    analytic_coherence_matrix,combine_matrix_moment)


def test_all_pauli_components_reconstruct_the_original_source_without_gram_fitting():
    F=np.array([[1.2,.3+.2j],[-.1j,.8]],complex)
    row=real_closed_source_matrix(1.2,F,.7+.6j,kappa=.24,omega=3.9,mass=np.pi/2)
    assert row['identity_residual']<5e-16
    assert not row['Gram_replaced_with_identity']
    assert not row['reflection_modulus_forced_to_one']
    assert row['source_accuracy_bound'] is None
    # A half-identity substitution cannot reproduce this nonisometric control.
    C=row['source_covariance'];v,w=F.T
    wrong=.5*np.eye(2)+(C[0,0]-.5)*(np.outer(v,v.conj())-np.outer(w,w.conj()))
    assert np.max(abs(wrong-row['covariance']))>.1


def test_pauli_probe_does_not_lose_sigma1_when_massless_stress_vertices_do():
    C=.5*PAULI[0]+.2*PAULI[1]+.1*PAULI[2]-.05*PAULI[3]
    traces=np.einsum('iab,ba->i',PAULI,C)
    np.testing.assert_allclose(from_pauli_traces(traces),C,atol=1e-16)
    traces[1]=0
    assert np.max(abs(from_pauli_traces(traces)-C))>.19


def test_analytic_coherence_preserves_independent_dual_and_returns_no_state():
    upper=np.array([[1,.2j],[.3,1.1]],complex)
    dual=np.array([[.9,-.1j],[.25,1.2+.1j]],complex)
    row=analytic_coherence_matrix(1.2+.03j,upper,dual,.4+.2j,kappa=.24)
    assert row['identity_residual']<1e-20
    assert not row['complex_frequency_is_Gaussian_state']
    assert row['contour_remainder_bound'] is None
    other=analytic_coherence_matrix(1.2+.03j,upper,upper,.4+.2j,kappa=.24)
    assert np.max(abs(row['analytic_matrix']-other['analytic_matrix']))>1e-10
    with pytest.raises(ValueError,match='real positive'):
        real_closed_source_matrix(1.2+.03j,upper,.4+.2j,kappa=.24,omega=3.9,mass=np.pi/2)


def test_matrix_moment_completes_hermitian_coherence_not_each_entry_real_part():
    D=.4*PAULI[0]+.1*PAULI[3]
    K=np.array([[.1j,.2+.3j],[.4j,.05]],complex)
    M=combine_matrix_moment(D,K)
    np.testing.assert_array_equal(M,M.conj().T)
    np.testing.assert_array_equal(M,D+K+K.conj().T)
    assert np.max(abs(M-(D+2*K.real)))>.1
