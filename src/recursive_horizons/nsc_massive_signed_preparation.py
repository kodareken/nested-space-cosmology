"""Massive signed-spinor map and a scoped seed-data sufficiency certificate.

The paired angular matrices are imported from NSC-8. No spatial covariance is
assigned here. Covariance completions below are algebraic control witnesses,
not choices of a horizon-prepared state.
"""
from __future__ import annotations

import numpy as np

from .nsc_chiral_boundary import ALPHA, ANGULAR, BETA
from .nsc_common_time_bulk_split import CommonTimeBulkSplit
from .nsc_lorentzian import geometry
from .nsc_transmitting_dirac_domain import I2, S1, S2, S3, TransmittingDiracSeamDomain


I4=np.eye(4,dtype=complex)
PAIR_TO_T=np.kron(S1,(I2+1j*S2)/np.sqrt(2))
PAIR_SIGNED_B=np.kron(S1,S3)
STORED_ANGULAR_PROJECTOR=np.diag([1.,1.,0.,0.])


def paired_T_symbol(momentum,*,N,q_PG,beta,radius,mass,angular):
    kinetic=momentum*np.kron(I2,CommonTimeBulkSplit.principal(N,q_PG,beta))
    return kinetic-N*mass*np.kron(I2,S1)-N*angular/radius*np.kron(S3,S3)


def signed_map_kind(mass,angular):
    if mass<0 or angular<0 or not np.isfinite([mass,angular]).all():
        raise ValueError('finite nonnegative retained labels required')
    if mass==0 and angular==0:
        raise ValueError('LLL is already prepared by its own immutable owner')
    if mass==0:
        return {'kind':'single angular-mass radial block','dimension':2,'current_B':S1,
                'seed_B':-1j*I2,'paired_covariance_required_by_mixed_mass_test':False}
    if angular==0:
        return {'kind':'single compact-mass radial block','dimension':2,'current_B':S3,
                'seed_B':S3,'paired_covariance_required_by_mixed_mass_test':False}
    return {'kind':'paired angular radial block','dimension':4,'current_B':PAIR_SIGNED_B,
            'seed_B':PAIR_SIGNED_B,'paired_covariance_required_by_mixed_mass_test':True}


def map_residuals(mass,angular):
    """Apply the exact basis/current maps to this group's owned PG expression."""
    definition=signed_map_kind(mass,angular)
    identity=I4 if definition['dimension']==4 else I2
    b=definition['current_B'];maximum_symbol=maximum_basis=0.
    for rho in (0.,.5,-.5):
        radius=np.sqrt(1+rho*rho);beta=float(geometry(rho)[0]);p=.2
        if definition['dimension']==4:
            current=np.kron(I2,S2-beta*I2)
            h=paired_T_symbol(p,N=1.,q_PG=1.,beta=beta,radius=radius,mass=mass,angular=angular)
            reversed_h=paired_T_symbol(-p,N=1.,q_PG=1.,beta=beta,radius=radius,mass=mass,angular=angular)
            original=p*(ALPHA-beta*I4)+angular/radius*ANGULAR+mass*BETA
            maximum_basis=max(maximum_basis,float(np.linalg.norm(PAIR_TO_T@original@PAIR_TO_T.conj().T-h)))
        else:
            current=S2-beta*I2
            potential=CommonTimeBulkSplit.local_potential(1.,radius,mass,angular)
            h=p*current+potential;reversed_h=-p*current+potential
        maximum_symbol=max(maximum_symbol,float(np.linalg.norm(b@h.conj()@b.conj().T+reversed_h)))
    beta0=float(geometry(0.)[0]);seam=TransmittingDiracSeamDomain(1.,1.,beta0,1.)
    trace,inverse,gram=seam.trace_map(np.ones(1))
    if definition['dimension']==4:
        t=np.kron(I2,trace[0]);ti=np.kron(I2,inverse[0]);g=np.kron(I2,gram[0])
    else:t,ti,g=trace[0],inverse[0],gram[0]
    seed_map=ti@b@t.conj()
    return {
        'T_basis_Hamiltonian':maximum_basis,
        'signed_frequency_symbol':maximum_symbol,
        'antiunitary_square':float(np.linalg.norm(b@b.conj()-identity)),
        'signed_map_unitarity':float(np.linalg.norm(b.conj().T@b-identity)),
        'T_normal_current':float(np.linalg.norm(b.conj().T@g@b-g.conj())),
        'T_to_seed_signed_map':float(np.linalg.norm(seed_map-definition['seed_B'])),
    }


def constant_single_block_obstruction():
    """Kernel of all three required Pauli anticommutators is zero.

    For m*lambda!=0, two distinct radii separate the sigma1 and sigma3
    requirements. The principal matrix supplies sigma2. A unitary 2x2
    matrix cannot satisfy these equations throughout the throat.
    """
    basis=[np.eye(4,dtype=complex)[i].reshape(2,2) for i in range(4)]
    linear=np.stack([np.concatenate([(b@s+s@b).ravel() for s in (S1,S2,S3)]) for b in basis],axis=1)
    singular=np.linalg.svd(linear,compute_uv=False)
    return {'linear_rank':int(np.linalg.matrix_rank(linear)),
            'complex_unknowns':4,'kernel_dimension':4-int(np.linalg.matrix_rank(linear)),
            'smallest_singular_value':float(singular.min()),
            'unitary_solution_exists':False}


def exact_pairing_identities():
    """Exact coefficient matching to the imported paired operator, once."""
    import sympy as sp
    exact=lambda a:sp.Matrix(a).applyfunc(sp.nsimplify)
    s1,s2,s3=map(exact,(S1,S2,S3));i2=sp.eye(2);k=sp.kronecker_product
    w=k(s1,(i2+sp.I*s2)/sp.sqrt(2));b=k(s1,s3)
    relations={
        'W_unitary':w*w.H-sp.eye(4),
        'principal_coefficient':w*exact(ALPHA)*w.H-k(i2,s2),
        'angular_coefficient':w*exact(ANGULAR)*w.H+k(s3,s3),
        'mass_coefficient':w*exact(BETA)*w.H+k(i2,s1),
        'B_principal_conjugation':b*k(i2,s2).conjugate()*b.H-k(i2,s2),
        'B_mass_sign':b*k(i2,s1).conjugate()*b.H+k(i2,s1),
        'B_angular_sign':b*k(s3,s3).conjugate()*b.H+k(s3,s3),
    }
    result={}
    for name,value in relations.items():
        if sp.simplify(value)!=sp.zeros(*value.shape):raise AssertionError(name)
        result[name]='0'
    basis=[sp.eye(4)[:,j].reshape(2,2) for j in range(4)]
    linear=sp.Matrix.hstack(*[sp.Matrix.vstack(*[(x*s+s*x).reshape(4,1) for s in (s1,s2,s3)]) for x in basis])
    if linear.nullspace():raise AssertionError('unexpected constant single-block intertwiner')
    result['single_block_exact_kernel_dimension']=0
    return result


def mixed_pointwise_probe(mass,angular,rho=.5):
    """Quantify why an instantaneous axis flip is not a domain intertwiner."""
    if mass<=0 or angular<=0:raise ValueError('mixed nonzero masses required')
    radius=np.sqrt(1+rho*rho);b0=angular;b1=angular/radius
    fixed=(b0*S1-mass*S3)/np.hypot(mass,b0)
    potential=-mass*S1-b1*S3
    local=(b1*S1-mass*S3)/np.hypot(mass,b1)
    derivative_b=-angular*rho/radius**3
    energy=np.hypot(mass,b1);energy_prime=b1*derivative_b/energy
    derivative=derivative_b/energy*S1-energy_prime/energy*local
    v=CommonTimeBulkSplit.principal(1.,1.,float(geometry(rho)[0]))
    connection=-1j*v@derivative@local.conj().T
    return {
        'rho':rho,
        'fixed_seed_axis_potential_mismatch':float(np.linalg.norm(fixed@potential.conj()@fixed.conj().T+potential)),
        'pointwise_potential_flip_residual':float(np.linalg.norm(local@potential.conj()@local.conj().T+potential)),
        'required_basis_connection_norm':float(np.linalg.norm(connection)),
        'connection_added_to_model':False,
    }


def paired_seed_data_witness(observed_seed_block,*,tolerance=3e-11):
    """Nonphysical information witnesses: X fixed, opposite marginal free.

    Neither Y=0 nor Y=I is adopted as a state. They satisfy the recorded
    marginal, CAR and the radial signed-E pairing but are not asserted to
    satisfy an uncomputed horizon-to-paired-mode preparation map.
    """
    x=np.asarray(observed_seed_block,dtype=complex)
    if x.shape!=(2,2) or not np.isfinite(x).all() or np.linalg.norm(x-x.conj().T)>tolerance:
        raise ValueError('a Hermitian seed control block is required')
    eig=np.linalg.eigvalsh(x)
    if eig.min()<-tolerance or eig.max()>1+tolerance:raise ValueError('seed control violates CAR')
    zero=np.zeros((2,2),complex)
    first=np.block([[x,zero],[zero,zero]])
    second=np.block([[x,zero],[zero,I2]])
    negative=[I4-PAIR_SIGNED_B@c.conj()@PAIR_SIGNED_B.conj().T for c in (first,second)]
    spectrum=np.concatenate([np.linalg.eigvalsh(c) for c in (first,second,*negative)])
    return {
        'same_observed_marginal_residual':float(np.linalg.norm(first[:2,:2]-second[:2,:2])),
        'known_charge_paired_marginal_residual':float(np.linalg.norm(negative[0][2:,2:]-negative[1][2:,2:])),
        'unobserved_positive_marginal_difference':float(np.linalg.norm(second[2:,2:]-first[2:,2:])),
        'unobserved_negative_same_angular_difference':float(np.linalg.norm(negative[1][:2,:2]-negative[0][:2,:2])),
        'witness_CAR_lower_violation':max(0.,-float(spectrum.min())),
        'witness_CAR_upper_violation':max(0.,float(spectrum.max())-1.),
        'physical_horizon_preparation_satisfied':None,
        'witness_adopted_as_spatial_state':False,
    }


def require_physical_pg_preparation(*,seed_covariance=None,lll_map=None,multiplicity=None):
    if seed_covariance is not None:raise ValueError('seed covariance is control only, never a spatial preparation input')
    if lll_map is not None:raise ValueError('LLL flow map does not apply to massive mixing')
    if multiplicity is not None:raise ValueError('multiplicity does not specify signed spinors or a physical covariance map')
    raise ValueError('global horizon/infinity-to-paired-PG modes and their state correlations are not supplied')
