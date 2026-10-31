"""Independent representation, full-residual and normalized-overlap controls."""
import math

import numpy as np
import pytest
from scipy.special import roots_jacobi
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_s3_frame import _level
from recursive_horizons.nsc_atlas_valve_neck import solve_necks,neck_diagnostics
from recursive_horizons.nsc_atlas_valve_dirac_modes import (
    radial_operator_blocks,calculate_dirac_level,mode_radial_profiles,omega_bounds,
)


@pytest.fixture(scope='module')
def initial():return solve_necks(40,(1.,))[0]


def test_all_round_sectors_match_original_pauli_frame_operator_and_multiplicities():
    L=4;u,w=roots_jacobi(30,.5,.5)
    reduced=[]
    for k in range(L+1):
        A,M,_,_,_=radial_operator_blocks(L,k,u,w,np.ones_like(u))
        np.testing.assert_allclose(M,np.eye(len(A)),atol=2e-14,rtol=0)
        reduced.extend(np.repeat(np.linalg.eigvalsh(A),2*k+2))
    pauli=(np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1.,-1.]))
    original=[]
    for l in range(L+1):
        generators=_level(l)[2];n=(l+1)**2
        D=1.5*np.eye(2*n,dtype=complex)+sum(1j*np.kron(sigma,G.toarray()) for sigma,G in zip(pauli,generators))
        original.extend(np.linalg.eigvalsh(D))
    assert len(reduced)==2*sum((l+1)**2 for l in range(L+1))
    np.testing.assert_allclose(sorted(reduced),sorted(original),atol=2e-14,rtol=0)
    assert sum(abs(np.asarray(reduced)-1.5)<1e-12)==2
    assert sum(abs(np.asarray(reduced)+1.5)<1e-12)==2


def test_full_pde_residual_detects_output_outside_the_ritz_space(initial):
    summary,_=calculate_dirac_level(4,initial=initial)
    for mode in summary['lowest_by_sign']:
        assert mode['projected_Ritz_residual']<2e-13
        assert mode['full_PDE_residual_physical']>.01
        assert mode['positive_residual_gamma_squared']>0
        assert mode['Momega2_identity_absolute_gap']<4e-14
        assert abs(mode['physical_normalization']-1)<2e-14


def test_angular_spinors_and_radial_participation_agree_on_actual_geometry(initial):
    summary,modes=calculate_dirac_level(8,initial=initial)
    mode=next(m for m in modes if m['sign']>0);record=next(m for m in summary['lowest_by_sign'] if m['sign']>0)
    assert mode['k']==0
    u,w=roots_jacobi(240,.5,.5);omega=initial.values(u)**2/initial.parent.H;F,G=mode_radial_profiles(mode,u)
    pauli=(np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1.,-1.]))
    directions=np.r_[np.eye(3),-np.eye(3)];densities=[]
    for direction in directions:
        sigma=sum(x*s for x,s in zip(direction,pauli))
        f=(F[:,None,None]*np.eye(2)-1j*G[:,None,None]*sigma)/math.sqrt(4*math.pi)
        psi=f/omega[:,None,None]
        # Trace divided by the angular degeneracy2 is the normalized projector mean.
        densities.append(np.sum(abs(psi)**2,axis=(1,2))/2)
    densities=np.asarray(densities)
    np.testing.assert_allclose(densities,np.broadcast_to(densities[0],densities.shape),atol=3e-18,rtol=1e-14)
    angular_norm=4*math.pi*float(w@(omega**3*np.mean(densities,axis=0)))
    angular_I4=4*math.pi*float(w@(omega**3*np.mean(densities**2,axis=0)))
    assert angular_norm==pytest.approx(1.,abs=2e-14)
    assert angular_I4==pytest.approx(record['invariant_ground_projector_I4'],abs=2e-16)
    assert record['effective_volume']<=neck_diagnostics(initial)['proper_total_volume']


def test_finite_omega_bounds_and_compact_inverse_identify_the_first_cluster(initial):
    bounds=omega_bounds(initial);u=np.linspace(-1,1,2001);omega=initial.values(u)**2/initial.parent.H
    assert bounds['Omega_lower']<=min(omega)<=max(omega)<=bounds['Omega_upper']
    with threadpool_limits(limits=1):summary,_=calculate_dirac_level(24,initial=initial,mass=2.)
    check=summary['compact_inverse_validation']
    assert check['isolation_margin']>0
    assert check['finite_ground_cluster_dimension']==2
    assert check['residual_interval_inside_isolated_global_first_cluster']
    assert summary['spinor_dimension']==summary['expected_dimension']
    assert max(m['full_PDE_residual_physical'] for m in summary['lowest_by_sign'])<1e-6
    for mode in summary['lowest_by_sign']:
        assert mode['energy']>mode['declared_mass']
        assert mode['lambda_over_energy']<1


def test_invalid_cutoff_and_nonpositive_weight_are_rejected():
    u,w=roots_jacobi(20,.5,.5)
    with pytest.raises(ValueError,match='integer degree'):radial_operator_blocks(True,0,u,w,np.ones_like(u))
    with pytest.raises(ValueError,match='positive'):radial_operator_blocks(4,0,u,w,np.zeros_like(u))
