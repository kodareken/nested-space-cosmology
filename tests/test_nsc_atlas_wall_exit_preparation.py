from decimal import Decimal
from fractions import Fraction as F
import math
import numpy as np
import pytest
import sympy as s
from recursive_horizons.nsc_atlas_wall_exit import calculation,wall_projection
from recursive_horizons.nsc_atlas_volume_preparation import PreparedBackground,PreparedModel
from recursive_horizons.nsc_atlas_curvature_mechanisms import FlatVolumeControl


def test_unsquared_junction_and_curvature_nonuniqueness():
    R,M,k,Ui,Uo=s.symbols('R M k Ui Uo',positive=True)
    beta=-M/(k*R**2)+(Uo-Ui+3*k*k)*R/(6*k)
    H2=(beta**2-1+2*M/R+Ui*R*R/3)/(R*R)
    assert s.simplify((beta-k*R)**2-(1-Uo*R*R/3+H2*R*R))==0
    r=calculation()
    assert F(r['reference_projection']['C_over_E0_exact'])==F(61720739711,1016191029600)
    assert r['witnesses'][0]['C']>0>r['witnesses'][2]['C']
    assert abs(r['witnesses'][1]['C'])<1e-12
    assert not r['post_exit_C_over_E0_derived']


def test_HM_loaded_force_and_fixed_energy_prescription():
    r=calculation();b=PreparedBackground(2.33,preparation='HM_location')
    state=b.state([0]);assert state[3][0]==pytest.approx(.375)
    derivative=b.model.averaged_rhs(0.,[0,b.initial_s,0])
    assert derivative[2]==pytest.approx(3/112,abs=1e-14)
    future=b.state([.1,1.,2.])
    assert np.all(future[0]>b.initial_s) and np.all(future[1]>0)
    assert r['HM']['unique_exit_efolds'] is None
    assert r['HM']['linear_growth_per_efold']==pytest.approx(.95204059745,rel=1e-10)


def test_volume_pressure_exchange_exactly():
    N,s0,x,v,E,H,V,Vp=s.symbols('N s0 x v E H V Vp',real=True)
    K=s.Integer(84);rho=K*v*v/2+V-s0;p=K*v*v/2-V+s0
    vdot=-3*H*v+(-Vp+6*E)/K
    rhodot=s.diff(rho,v)*vdot+s.diff(rho,V)*Vp*v
    assert s.simplify(rhodot+3*H*(rho+p)-6*E*v)==0
    assert s.simplify(p/rho-2*E*v/(H*rho)-(-1-rhodot/(3*H*rho)))==0


def test_original_state_recovered_and_minimum_w_is_undefined_at_zero_energy():
    a=PreparedBackground(2.33);b=FlatVolumeControl(2.33)
    z=np.array([.1,1.,2.]);theta=.4256503406128458
    assert np.allclose(a.distances(theta,z,z)[0],b.distances(theta,z,z)[0],rtol=2e-10)
    row=a.eos_rows(theta)[0]
    assert row['w_s']==pytest.approx(-.8867355268,abs=2e-8)
    assert row['w_dark_energy']==pytest.approx(-.9774713869,abs=2e-8)
    m=PreparedBackground(2.33,preparation='minimum')
    assert math.isnan(m.eos([0])['w_s'][0])
    assert m.model.averaged_rhs(0.,[0,m.initial_s,0])[2]>0


def test_displayed_and_unrounded_quantities_are_not_conflated():
    assert 3*(Decimal('0.577259')+Decimal('0.143326'))==Decimal('2.161755')
