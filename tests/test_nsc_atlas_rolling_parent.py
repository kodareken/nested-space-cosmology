"""Exact rolling limit and independent finite-geometry controls, not state selection."""
from fractions import Fraction as F
import math

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_sphere_flux_completion import (
    rolling_parent_limit,rolling_parent_state,symmetric_radius_branches,all_volume_extrema,wall_energy_jump,
)
from recursive_horizons.nsc_atlas_shape_transition import evolve_rolling_parent


def test_exact_q2_limit_from_canonical_exponential_and_16D_trace():
    row=rolling_parent_limit(2)
    A=F(215,32)+F(2)**2/8
    slope_squared=F(12)**2/F(84)
    power=2/slope_squared
    w=slope_squared/3-1
    assert row['A']==A==F(231,32)
    assert row['Einstein_power']==power==F(7,6)
    assert row['Einstein_w']==w==F(-3,7)
    assert row['bulk_H_squared']==A/105==F(11,160)
    assert row['bulk_Ricci_scalar']==16*15*row['bulk_H_squared']==F(33,2)
    section_count=15+math.comb(15,2)
    assert row['bulk_Kretschmann_scalar']==4*section_count*row['bulk_H_squared']**2==F(363,160)
    assert row['ds_dN']==F(1,7)
    assert row['kinetic_times_Einstein_time_squared']==F(84,2)/36
    assert row['density_times_Einstein_time_squared']==3*power*power
    assert row['potential_times_Einstein_time_squared']==3*power*power-F(84,2)/36


def independent_bulk_geometry(s,v,q):
    """Differentiate the two physical scale factors using E-frame equations.

    dTau16=exp(-6s)dtE; a_external=exp(-6s)aE and a_internal=exp(s).
    This route uses neither v_N nor any returned bulk-rate derivative.
    """
    A=float(F(215,32))+q*q/8;K=84.
    V=A*math.exp(-12*s)-7.5*math.exp(-14*s)+25/32*math.exp(-24*s)
    Vs=-12*A*math.exp(-12*s)+105*math.exp(-14*s)-75/4*math.exp(-24*s)
    HE=math.sqrt(V/(3-K*v*v/2));sdot=v*HE
    sddot=-3*HE*sdot-Vs/K
    HEdot=-K*sdot*sdot/2
    internal=math.exp(6*s)*sdot
    external=math.exp(6*s)*(HE-6*sdot)
    internal_dot=math.exp(12*s)*(sddot+6*sdot*sdot)
    external_dot=math.exp(12*s)*(HEdot-6*sddot+6*sdot*(HE-6*sdot))
    expansion=3*external+12*internal
    Ricci=(2*(3*external_dot+12*internal_dot)+expansion**2+
        3*external**2+12*internal**2+15*math.exp(-2*s))
    # Full16D trace: vacuum coefficient16A/7 plus canonicalF6 trace25x^6/56.
    source=16*A/7+25/56*math.exp(-12*s)
    return internal,external,internal_dot,external_dot,Ricci,source


@pytest.mark.parametrize('s,v',[(0.,.1),(.2,.03),(.7,.18),(1.3,-.1)])
def test_finite_source_and_metric_Ricci_from_independent_scale_differentiation(s,v):
    row=rolling_parent_state(s,v)
    Hi,He,Hid,Hed,R,source=independent_bulk_geometry(s,v,2.)
    np.testing.assert_allclose([row['bulk_internal_H'],row['bulk_external_H'],
        row['bulk_internal_H_dot'],row['bulk_external_H_dot']],
        [Hi,He,Hid,Hed],rtol=3e-13,atol=3e-13)
    assert R==pytest.approx(source,rel=3e-13,abs=4e-13)
    assert row['bulk_Ricci_metric']==pytest.approx(R,rel=3e-13,abs=4e-13)
    assert row['bulk_Ricci_source']==pytest.approx(source,rel=3e-15)
    # Enumerate the15 actual directions rather than reuse grouped curvature
    # multiplicities from the source. Only planes within one S3 have k/a².
    directions=[('external',None,He,Hed) for _ in range(3)]+[
        ('internal',factor,Hi,Hid) for factor in range(4) for _ in range(3)]
    sections=[(rate_dot+rate*rate)**2 for _,_,rate,rate_dot in directions]
    counts=dict(external_pair=0,external_internal=0,same_sphere=0,different_spheres=0)
    for i,left in enumerate(directions):
        for right in directions[i+1:]:
            value=left[2]*right[2]
            if left[0]==right[0]=='external':kind='external_pair'
            elif left[0]!=right[0]:kind='external_internal'
            elif left[1]==right[1]:
                kind='same_sphere';value+=5/8*math.exp(-2*s)
            else:kind='different_spheres'
            counts[kind]+=1;sections.append(value*value)
    assert counts==dict(external_pair=3,external_internal=36,same_sphere=12,different_spheres=54)
    assert len(sections)==15+math.comb(15,2)
    assert row['bulk_Kretschmann_metric']==pytest.approx(4*math.fsum(sections),rel=5e-13,abs=5e-13)


def test_both_physical_scale_rates_approach_the_same_bulk_H():
    target=math.sqrt(float(F(11,160)));errors=[]
    for s in (4.,8.):
        row=rolling_parent_state(s,1/7)
        assert row['bulk_external_H']==pytest.approx(row['bulk_internal_H'],rel=3e-15)
        errors.append(max(abs(row['bulk_external_H']-target),abs(row['bulk_internal_H']-target)))
    assert errors[1]<errors[0]/1000 and errors[1]<1e-6
    assert row['bulk_Kretschmann_metric']==pytest.approx(float(F(363,160)),rel=2e-6)


def test_q2_restricted_stationary_parent_and_reference_default_are_preserved():
    q2=symmetric_radius_branches(2)
    assert q2['q_reference']==2 and all_volume_extrema(2)==()
    assert len(q2['branches'])==1
    row=q2['branches'][0]
    assert row['branch']=='1+3' and row['retained_radius_local_minimum']
    assert row['V']==pytest.approx(-.5236104449501537,abs=3e-13)
    np.testing.assert_allclose(row['canonical_mass_squared'],
        [.10940083852642593,3.1891775791082533,3.189177579108256,3.6680108507736313],rtol=4e-12,atol=3e-13)
    assert next(r for r in q2['isolations'] if r['branch']=='2+2')['positive_root_intervals']==()
    reference=symmetric_radius_branches()
    assert reference==symmetric_radius_branches(1.5)
    assert reference['q_reference']==F(3,2) and len(reference['branches'])==6
    false=next(r for r in reference['branches'] if r.get('volume_kind')=='minimum')
    assert false['V']==pytest.approx(.22098302421260713,abs=3e-14)


def test_two_declared_finite_histories_approach_limit_without_selecting_a_history():
    seeds=((.2,0.),(.8,.2));solutions=[];H=rolling_parent_limit()['bulk_H']
    with threadpool_limits(limits=1):
        for s0,v0 in seeds:
            sol=evolve_rolling_parent(initial_s=s0,initial_v=v0,efolds=40.,budget_seconds=5.)
            np.testing.assert_array_equal(sol.y[:,0],[s0,v0,0.])
            assert sol.success and sol.t[0]==0. and sol.t[-1]==40.
            start=rolling_parent_state(s0,v0);end=rolling_parent_state(sol.y[0,-1],sol.y[1,-1])
            assert abs(end['v']-1/7)<2e-6
            before=max(abs(start['bulk_internal_H']-H),abs(start['bulk_external_H']-H))
            after=max(abs(end['bulk_internal_H']-H),abs(end['bulk_external_H']-H))
            assert after<before/1000 and after<1e-4
            assert sol.y[2,-1]>0
            solutions.append(sol)
    # Equal attractor rates do not erase the freely supplied finite histories.
    assert abs(solutions[0].y[0,-1]-solutions[1].y[0,-1])>.1
    assert abs(solutions[0].y[2,-1]-solutions[1].y[2,-1])>.1


def test_material_wall_Ward_identity_fixes_only_flux_difference_and_area_work():
    # Local wall-jet control, not a constraint-solved bulk preparation. The
    # spacelike normal points child->parent; U is the common wall tangent.
    s=.07;s_tangent=.03;phi_tangent=.11;phi_normal=-.07
    sigma=.2*math.exp(-6*s);jump=-6*sigma/84;sigma_dot=-6*sigma*s_tangent
    ward=wall_energy_jump(s,s_tangent)
    metric=np.diag([-1.,1.]);U=np.array([1.,0.]);normal=np.array([0.,1.])
    def stress(s_normal,potential):
        ds=np.array([s_tangent,s_normal]);dphi=np.array([phi_tangent,phi_normal])
        return (84*np.outer(ds,ds)+np.outer(dphi,dphi)
            -.5*metric*(84*(ds@metric@ds)+(dphi@metric@dphi)+2*potential))
    R=1.7;Rdot=.06;area=4*math.pi*R*R;area_dot=8*math.pi*R*Rdot
    def energy(tau):return 4*math.pi*(R+Rdot*tau)**2*.2*math.exp(-6*(s+s_tangent*tau))
    h=1e-4;energy_dot=(energy(h)-energy(-h))/(2*h)
    individual=[]
    for mean in (-.08,.05):
        child=mean-jump/2;parent=mean+jump/2
        Tchild=float(normal@stress(child,.31)@U)
        Tparent=float(normal@stress(parent,.57)@U)
        assert Tparent-Tchild==pytest.approx(sigma_dot,rel=3e-15,abs=3e-16)
        assert ward['T_nU_parent_minus_child']==pytest.approx(Tparent-Tchild,abs=3e-16)
        assert ward['F_child_minus_parent']==pytest.approx(Tparent-Tchild,abs=3e-16)
        assert ward['sigma_tangent']==pytest.approx(sigma_dot,abs=3e-17)
        # Surface pressure is -sigma, so this retains the area-work term.
        assert energy_dot-sigma*area_dot==pytest.approx(area*(Tparent-Tchild),rel=2e-9,abs=1e-9)
        individual.append((Tchild,Tparent))
    # Both individual normal energy contractions can change sign although
    # their Ward difference and the surface ledger remain fixed.
    assert individual[0][0]*individual[1][0]<0
    assert individual[0][1]*individual[1][1]<0
