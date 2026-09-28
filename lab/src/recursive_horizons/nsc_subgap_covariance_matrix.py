"""Full source covariance probes for the owned closed incoming channel.

This extends the existing scalar subgap bilinears to all Pauli components.
It is a source restriction identity, not an evolved incoming state or a
certified contour quadrature. Complex inputs remain analytic bilinears.
"""
import numpy as np
from .nsc_transmitting_dirac_domain import I2,S1,S2,S3
from .nsc_paired_horizon_preparation import source_covariance
from .nsc_subgap_history_response import source_functions

PAULI=np.asarray([I2,S1,S2,S3],complex)


def _basis(value,name):
    a=np.asarray(value,complex)
    if a.shape!=(2,2) or not np.isfinite(a).all():raise ValueError('finite two-column '+name+' required')
    return a


def pauli_bilinears(canonical_upper,canonical_dual):
    """B_i,ab=Phi_a(bar z)^dagger sigma_i Phi_b(z), including sigma0=I."""
    upper,dual=_basis(canonical_upper,'basis'),_basis(canonical_dual,'independent dual')
    return np.asarray([dual.conj().T@vertex@upper for vertex in PAULI])


def from_pauli_traces(traces):
    values=np.asarray(traces,complex)
    if values.shape!=(4,) or not np.isfinite(values).all():raise ValueError('all four finite Pauli traces required')
    return .5*np.einsum('i,iab->ab',values,PAULI)


def real_closed_source_matrix(energy,basis,reflection,*,kappa,omega,mass):
    """Exactly reconstruct the SAME real source sandwich for columns(Rv,w,0).

    Gram and reflection-modulus deviations are retained, not normalized away.
    No adiabatic reference is subtracted from this full covariance.
    """
    if np.iscomplexobj(energy):raise ValueError('real positive source frequency required')
    energy=float(energy);R=complex(reflection)
    if not 0<energy<min(mass,omega*mass) or not np.isfinite(R):raise ValueError('owned positive subgap closed channel required')
    F=_basis(basis,'unsewn canonical basis')
    C=source_covariance(energy,kappa,omega,mass)
    if C[2,2]!=0 or np.max(abs(C[:2,2]))!=0:raise ValueError('this identity is for the closed incoming channel')
    B=pauli_bilinears(F,F)
    diagonal=C[0,0]*abs(R)**2*B[:,0,0]+C[1,1]*B[:,1,1]
    coherence=C[0,1]*R*B[:,1,0]
    traces=diagonal+coherence+coherence.conj()
    matrix=from_pauli_traces(traces)
    physical=np.column_stack((R*F[:,0],F[:,1],np.zeros(2,complex)))
    direct=physical@C@physical.conj().T
    return {'covariance':matrix,'direct_source_covariance':direct,
        'pauli_traces':traces,'diagonal_traces':diagonal,'coherence_traces':coherence,
        'source_covariance':C,'complete_source_mode_map_returned':False,
        'identity_residual':float(np.max(abs(matrix-direct))),
        'Gram_replaced_with_identity':False,'reflection_modulus_forced_to_one':False,
        'adiabatic_reference_subtracted':False,'incoming_state_frozen_by_this_object':False,
        'source_accuracy_bound':None,'physical_local_gate':'OPEN'}


def analytic_coherence_matrix(energy,canonical_upper,canonical_dual,reflection,*,kappa):
    """K(z)=-i s(z) R(z) v(z) w(bar z)^dagger on the inherited thermal strip.

    On the real axis the occupied coherence is K+K^dagger. Off axis K is an
    analytic integration input, never a Gaussian covariance state.
    """
    z=complex(energy);R=complex(reflection)
    if z.real<=0 or not np.isfinite(R):raise ValueError('finite positive-real contour frequency and reflection required')
    _,s=source_functions(z,kappa)
    upper,dual=_basis(canonical_upper,'analytic basis'),_basis(canonical_dual,'independent dual')
    B=pauli_bilinears(upper,dual)
    traces=-1j*s*R*B[:,1,0]
    matrix=from_pauli_traces(traces)
    direct=-1j*s*R*np.outer(upper[:,0],dual[:,1].conj())
    return {'analytic_matrix':matrix,'pauli_traces':traces,
        'identity_residual':float(np.max(abs(matrix-direct))),
        'complex_frequency_is_Gaussian_state':False,'contour_remainder_bound':None,
        'full_source_covariance_continued':False,'physical_local_gate':'OPEN'}


def combine_matrix_moment(real_diagonal_integral,analytic_coherence_integral):
    """Restore Hermitian source moment M=diagonal+K+K^dagger after integration.

    Real polynomial weights may multiply both inputs before integration.
    No quadrature accuracy or Gaussian-column interpretation is supplied.
    """
    diagonal=_basis(real_diagonal_integral,'integrated real diagonal matrix')
    coherence=_basis(analytic_coherence_integral,'integrated coherence matrix')
    if np.max(abs(diagonal-diagonal.conj().T))>3e-11:raise ValueError('Hermitian real-axis diagonal moment required')
    return diagonal+coherence+coherence.conj().T
