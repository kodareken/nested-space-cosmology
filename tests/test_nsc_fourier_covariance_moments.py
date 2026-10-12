"""Unitary Fourier transport and sharp half-transfer moment bookkeeping."""
from fractions import Fraction as Q
import numpy as np
import pytest
pytest.importorskip('flint')
from flint import arb,ctx,fmpq
from scipy.linalg import expm
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_fourier_covariance_moments import (
    propagate_covariance_moments,pure_radius_potential_integrals,covariance_matter_error)

S1=np.array([[0,1],[1,0]],complex);S3=np.diag([1,-1]).astype(complex)


def test_fourier_diagonal_dirac_preserves_each_nuclear_symbol_norm():
    D=np.array([[.2,.3j],[.4+.1j,-.7]])
    norm=np.linalg.svd(D,compute_uv=False).sum()
    for k,omega in ((0.,2.),(20.,3.),(-100.,9.)):
        left=expm(-.2j*((k+omega/2)*S3+.4*S1))
        right=expm(-.2j*((k-omega/2)*S3+.4*S1))
        transformed=left@D@right.conj().T
        assert np.linalg.svd(transformed,compute_uv=False).sum()==pytest.approx(norm,rel=3e-14)


def test_shifted_compact_kernel_requires_the_full_first_moment_coefficient():
    # D0(k)=I*pi/eps on[-eps,eps] has M0=2,M1=eps. The two
    # Hermitian potential modes +/-omega have op-norm v. Their shifted
    # supports are disjoint. The commutator has exact M0=8v,M1=4v*omega.
    eps,v,omega=Q(1,100),Q(1,5),Q(3)
    V0,V1=2*v,2*v*omega;M0,M1=Q(2),eps
    actual0,actual1=8*v,4*v*omega
    assert actual0==2*V0*M0
    assert actual1<=2*V0*M1+V1*M0
    assert actual1>2*V0*M1+V1*M0/2
    with ctx.workprec(160):
        dt=Q(1,10**8)
        result=propagate_covariance_moments(M0,M1,0,0,V0*dt,V1*dt)
        first=restored_upper(result['first_moment_upper'])
        tangent=(2*V0*M1+V1*M0)
        assert first>=arb(fmpq((M1+dt*tangent).numerator,(M1+dt*tangent).denominator))


def test_zero_potential_retains_initial_and_defect_and_does_not_invent_a_gate():
    result=propagate_covariance_moments(Q(1,3),Q(1,5),Q(2,3),Q(4,5),0,0)
    with ctx.workprec(160):
        assert restored_upper(result['zeroth_moment_upper'])>=1
        assert restored_upper(result['first_moment_upper'])>=1
    assert result['physical_local_gate']=='OPEN'
    with pytest.raises(ValueError,match='explicit'):
        propagate_covariance_moments(0,0,None,0,0,0)
    with pytest.raises(ValueError,match='embedding'):
        covariance_matter_error(result,mass=1,absolute_angular=1,axial_lower=1,radius_lower=1,
            multiplicity=1,continuum_embedding_established=False)


def test_radius_potential_uses_original_clock_and_enclosing_denominator():
    with ctx.workprec(160):
        a=arb(1,arb(1)/10)
        result=pure_radius_potential_integrals(Q(1,10),Q(1,5),absolute_angular=2,
            axial_lower=a,rho_up=Q(103,100))
        assert restored_upper(result['potential_zero_integral_upper'])>=arb(1)/150
        assert restored_upper(result['potential_first_integral_upper'])>=arb(1)/75
