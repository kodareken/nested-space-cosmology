"""Charged-sector matching, exact vacuum classification and full-wrap controls."""
from dataclasses import FrozenInstanceError
from fractions import Fraction as F
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_flux_conversion import ScalarFlux
from recursive_horizons.nsc_atlas_parent_scale_selector import VolumeParent
from recursive_horizons.nsc_atlas_sphere_flux_completion import (
    action_reduction_record,combined_record,charged_potential,shape_potential,
    flux_components,all_volume_extrema,volume_thresholds,half_integer_classification,wrapped_brane_couplings,
)


def test_reference_shift_matches_old_two_scalar_action_without_double_counting():
    record=combined_record();model=record.scalar_sector
    assert isinstance(model,ScalarFlux) and model.epsilon==.25 and model.omega==20.
    assert record.LambdaBare==F(215,32) and record.bare_mass_squared==100 and record.mixing_squared==300
    assert model.mu**2==pytest.approx(300.,rel=3e-16)
    assert record.phi_star_reference==pytest.approx(-math.sqrt(300)*1.5/400,abs=2e-17)
    for s,chi in [(-.035646782552940987,.012),(.04,-.02),(.2,.1)]:
        phi=model.minimum(1.5)+chi
        expected=VolumeParent().potential(s)+.5*400*math.exp(-12*s)*chi**2
        assert charged_potential(s,phi,1.5)==pytest.approx(expected,rel=3e-13,abs=3e-14)
        assert model.potential(phi,1.5)==pytest.approx(7.+200*chi*chi,rel=3e-15)
    with pytest.raises(FrozenInstanceError):record.LambdaBare=F(7)


@pytest.mark.parametrize('s,phi,q',[(0.,-.04,1.5),(.03,.02,1.),(-.02,-.1,2.),(.21,.01,-.5)])
def test_direct16D_Routhian_and_inverse_metric_scaling(s,phi,q):
    r=combined_record();geometry=action_reduction_record();V0=r.reference_volume
    beta=np.array([.03,-.01,.02,-.04]);Phi=phi/math.sqrt(V0);qD=q/math.sqrt(V0)
    density=float(r.LambdaBare)/r.kappa16+.5*100*Phi*Phi+.5*(qD+math.sqrt(300)*Phi)**2
    norm=sum(geometry.q**2*math.prod(math.exp(-2*(s+beta[a//3])) for a in indices) for indices,_ in flux_components())
    curvature=sum(6/float(geometry.radius_squared)*math.exp(-2*(s+b)) for b in beta)
    direct=V0*math.exp(-12*s)*(density+norm/2-curvature/(2*r.kappa16))
    assert charged_potential(s,phi,q,beta)==pytest.approx(direct,rel=3e-13,abs=3e-14)


def test_all_volume_roots_and_canonical_normalization_match_independent_derivatives():
    model=combined_record().scalar_sector
    for q in (0.,.5,1.,1.5,-1.5,1.8):
        roots=all_volume_extrema(q);assert len(roots)==2
        assert {r['kind'] for r in roots}=={'minimum','maximum'}
        for row in roots:
            s=row['s'];A=F(215,32)+F(str(q))**2/8;B=F(15,2);T=F(25,32)
            first=-12*float(A)*math.exp(-12*s)+14*float(B)*math.exp(-14*s)-24*float(T)*math.exp(-24*s)
            second=144*float(A)*math.exp(-12*s)-196*float(B)*math.exp(-14*s)+576*float(T)*math.exp(-24*s)
            assert abs(first)<3e-11
            assert row['stiffness']==pytest.approx(second,rel=2e-12,abs=2e-11)
            assert row['canonical_mass_squared']==pytest.approx(second/84,rel=2e-12,abs=3e-13)
            assert row['vacuum']==pytest.approx(charged_potential(s,model.minimum(q),q),rel=3e-12,abs=2e-14)
    reference={r['kind']:r for r in all_volume_extrema(1.5)}
    for row in VolumeParent().vacuum_roots():
        assert reference[row['kind']]['s']==pytest.approx(row['s'],abs=2e-14)
        assert reference[row['kind']]['vacuum']==pytest.approx(row['vacuum'],abs=3e-14)


def test_exact_merger_Minkowski_and_all_half_integer_classification():
    b=volume_thresholds();B=F(15,2);T=F(25,32)
    assert b['turn_x_fifth']==7*B/(72*T)==F(14,15)
    assert b['minkowski_x_fifth']==B/(6*T)==F(8,5)
    assert b['critical_A_fifth']==F(35,36)**5*B**5*F(14,15)
    assert b['minkowski_A_fifth']==F(5,6)**5*B**5*F(8,5)
    lower,upper=b['positive_minimum_abs_q_interval']
    assert 1<lower<1.5<upper<2
    merger=all_volume_extrema(upper)
    assert len(merger)==1 and merger[0]['kind']=='degenerate' and merger[0]['stiffness']==0.
    minkowski=next(r for r in all_volume_extrema(lower) if r['kind']=='minimum')
    assert minkowski['vacuum_sign']=='zero' and abs(minkowski['vacuum'])<3e-13
    classification=half_integer_classification()
    for key in ('minkowski_fifth_power_bracket','merger_fifth_power_bracket'):
        a,b,c=classification[key];assert a<b<c and all(isinstance(v,F) for v in (a,b,c))
    for numerator in range(-8,9):
        q=F(numerator,2);roots=all_volume_extrema(q)
        if abs(q)>=2:assert roots==()
        else:
            minimum=next(r for r in roots if r['kind']=='minimum')
            assert minimum['vacuum_sign']==('positive' if abs(q)==F(3,2) else 'negative')


def test_wrapped_material_source_scaling_and_radion_variation():
    V0=combined_record().reference_volume
    for s in (-.035646782552940987,0.,.08):
        c=wrapped_brane_couplings(s)
        # Three external worldvolume directions and all12 internal directions.
        induced_measure=V0*math.exp(-18*s)*math.exp(12*s)
        sigma=(.2/V0)*induced_measure
        assert c['sigma']==pytest.approx(sigma,rel=2e-15)
        assert c['e_can']==pytest.approx(.5*math.exp(-6*s),rel=2e-15)
        assert c['omega']==pytest.approx(20*math.exp(-6*s),rel=2e-15)
        assert c['q_out_minus_in']==F(1,2)
        assert c['radion_normal_jump']==pytest.approx(-6*sigma/84,rel=2e-15)
        assert c['phi_direct_normal_jump']==0 and c['shape_direct_source']==(0.,)*4
        assert c['sigma']>0  # This source is not the zero-tension horizon.
    with pytest.raises(ValueError,match='finite real'):all_volume_extrema(float('inf'))
    with pytest.raises(ValueError,match='finite scalar'):charged_potential(0.,float('nan'),1.5)
