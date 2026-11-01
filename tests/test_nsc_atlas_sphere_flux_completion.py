"""Independent exterior, reduction and canonical shape-Hessian controls."""
from dataclasses import FrozenInstanceError
from fractions import Fraction as F
from itertools import combinations,permutations
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_parent_scale_selector import VolumeParent
from recursive_horizons.nsc_atlas_sphere_flux_completion import (
    action_reduction_record,flux_components,flux_stress_check,
    shape_potential,scalar_mass_squared,shape_mass_squared,shape_threshold,integer_flux_embedding,
)


def test_sparse_stress_matches_independent_full_antisymmetric_components():
    # Expand all6! tensor permutations, unlike the production exterior inner
    # products. The factorials below are essential to physical normalization.
    raw={}
    expected_pairs=list(combinations(range(4),2))
    for i,j in expected_pairs:
        indices=tuple(range(3*i,3*i+3))+tuple(range(3*j,3*j+3))
        for permutation in permutations(indices):
            inversions=sum(permutation[a]>permutation[b] for a in range(6) for b in range(a+1,6))
            raw[permutation]=(-1)**inversions
    norm=F(sum(value*value for value in raw.values()),math.factorial(6))
    grouped={}
    for (a,*tail),value in raw.items():grouped.setdefault(tuple(tail),{})[a]=value
    C=[[F(0) for _ in range(12)] for _ in range(12)]
    for values in grouped.values():
        for a,fa in values.items():
            for b,fb in values.items():C[a][b]+=F(fa*fb,math.factorial(5))
    result=flux_stress_check()
    assert result.norm_squared==norm==6
    assert result.internal_contraction==tuple(map(tuple,C))
    assert all(C[a][b]==(3 if a==b else 0) for a in range(12) for b in range(12))
    assert all(result.internal_stress[a][b]==0 for a in range(12) for b in range(12))
    assert result.external_metric_coefficient==-3
    assert {indices for indices,_ in flux_components()}=={
        tuple(range(3*i,3*i+3))+tuple(range(3*j,3*j+3)) for i,j in expected_pairs}


def test_action_normalization_reduces_to_the_fixed_coefficients():
    r=action_reduction_record();r2=r.radius_squared
    Rgamma=4*F(3*2)/r2
    assert Rgamma==15 and Rgamma/2==r.B
    assert F(12*14,2)==r.K==84
    assert r.reference_volume==pytest.approx((2*math.pi**2*float(r2)**1.5)**4,rel=2e-15)
    assert r.kappa16/r.reference_volume==1
    assert r.q*r.q*r.reference_volume==pytest.approx(float(F(25,96)),rel=3e-15)
    assert 3*r.q_squared_times_reference_volume==r.T==F(25,32)
    assert r.mass16*r.mass16==r.scalar_mass_squared_prefactor==400
    with pytest.raises(FrozenInstanceError):r.T=F(1)
    with pytest.raises(TypeError):flux_stress_check().internal_contraction[0][0]=F(4)


@pytest.mark.parametrize('s,beta',[
    (-.035646782552940987,(.03,-.01,.02,-.04)),
    (.04,(.08,-.06,.01,-.03)),
    (.21,(-.02,.01,.035,-.025)),
])
def test_general_shape_potential_follows_direct_metric_and_flux_scaling(s,beta):
    record=action_reduction_record();beta=np.array(beta)
    # Actual scalar curvature of the four rescaled3-spheres.
    Rinternal=sum(6/float(record.radius_squared)*math.exp(-2*(s+b)) for b in beta)
    determinant=record.reference_volume*math.exp(-12*s)
    curvature=-determinant*Rinternal/(2*record.kappa16)
    # Raise each of the six internal indices with its actual block metric.
    norm=0.
    for indices,_ in flux_components():
        norm+=record.q**2*math.prod(math.exp(-2*(s+beta[a//3])) for a in indices)
    flux=determinant*norm/2
    cosmological=determinant*float(record.Lambda16)/record.kappa16
    assert shape_potential(s,beta)==pytest.approx(cosmological+curvature+flux,rel=2e-13,abs=2e-14)


def test_zero_shape_is_one_existing_T_term_and_scalar_has_canonical_weight():
    parent=VolumeParent();points=np.array([-.035646782552940987,0.,.04,.21])
    np.testing.assert_allclose(shape_potential(points),[parent.potential(float(s)) for s in points],rtol=3e-14,atol=3e-15)
    r=action_reduction_record()
    for s in points:
        kinetic=r.reference_volume*math.exp(-12*s)*math.exp(12*s)/r.reference_volume
        mass=r.reference_volume*math.exp(-12*s)*float(r.mass16)**2/r.reference_volume
        assert kinetic==pytest.approx(1.,abs=3e-16)
        assert scalar_mass_squared(float(s))==pytest.approx(mass,rel=2e-15)


def _canonical_hessian(s,h):
    # Orthonormal tangent basis on sum beta=0; canonical coordinates divide
    # beta by sqrt3 because the independently reduced shape metric is3I.
    Q=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]],float)/2
    value=lambda z:shape_potential(s,Q@z/math.sqrt(3))
    origin=np.zeros(3);zero=value(origin);H=np.empty((3,3))
    for i in range(3):
        ei=np.eye(3)[i]*h
        H[i,i]=(value(ei)-2*zero+value(-ei))/h**2
        for j in range(i):
            ej=np.eye(3)[j]*h
            H[i,j]=H[j,i]=(value(ei+ej)-value(ei-ej)-value(-ei+ej)+value(-ei-ej))/(4*h*h)
    return H


@pytest.mark.parametrize('s',[-.035646782552940987,.04,.21])
def test_local_masses_match_independent_canonical_finite_difference_Hessian(s):
    coarse=_canonical_hessian(s,4e-4);fine=_canonical_hessian(s,2e-4)
    richardson=(4*fine-coarse)/3
    np.testing.assert_allclose(richardson,shape_mass_squared(s)*np.eye(3),rtol=3e-6,atol=2e-7)


def test_shape_threshold_and_trace_condition_are_not_fitted_or_projected():
    critical=shape_threshold()
    assert math.exp(10*critical)==pytest.approx(5/4,abs=3e-16)
    assert shape_mass_squared(critical)==pytest.approx(0.,abs=3e-15)
    assert shape_mass_squared(critical-.001)>0>shape_mass_squared(critical+.001)
    with pytest.raises(ValueError,match=r'sum\(beta\)=0'):shape_potential(0.,[.01,0.,0.,0.])
    with pytest.raises(ValueError,match='four finite'):shape_potential(0.,[0.,0.,0.])
    with pytest.raises(ValueError,match='finite breathing'):shape_mass_squared(float('nan'))


def test_same_effective_charge_and_energy_do_not_select_the_wrapping_vector():
    diagonal=integer_flux_embedding(F(1,96),(5,)*6)
    pair=integer_flux_embedding(F(1,16),(5,0,0,0,0,0))
    assert diagonal.reduced_flux_energy==pair.reduced_flux_energy==F(25,32)
    assert F(1,96)*6==F(1,16)
    assert diagonal.equal_radius_shape_acceleration_coefficients==(F(0),)*4
    assert pair.equal_radius_shape_acceleration_coefficients==(F(25,32),F(25,32),F(-25,32),F(-25,32))
    assert sum(pair.equal_radius_shape_acceleration_coefficients)==0
    # A primitive-cycle unit of1/4 cannot give an exactly isotropic25-norm
    # lattice vector: each degree is integer but isotropy would require25/2.
    assert F(pair.squared_integer_norm,2).denominator==2


def test_pair_shape_force_matches_independent_constrained_flux_gradient():
    h=1e-5;basis=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]],float)/2
    source=lambda z:float(F(25,32))*np.exp(-6*sum((basis@z)[:2]))
    gradient=np.array([(source(h*np.eye(3)[i])-source(-h*np.eye(3)[i]))/(2*h) for i in range(3)])
    force=-basis@gradient/3
    expected=integer_flux_embedding(F(1,16),(5,0,0,0,0,0)).equal_radius_shape_acceleration_coefficients
    np.testing.assert_allclose(force,list(map(float,expected)),rtol=7e-10,atol=2e-11)
