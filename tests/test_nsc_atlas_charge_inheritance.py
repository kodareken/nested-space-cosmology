"""Independent charge transport, normalization and ruler controls."""
import math
from fractions import Fraction as F
import pytest
import sympy as sp
from recursive_horizons.nsc_atlas_charge_inheritance import calculation,parke_ruler,volume_charge_ratio


def test_discharge_changes_flux_not_the_charge_or_gravitational_coupling():
    r=calculation()
    rows=r['flux_step']['chain']
    assert [row['flux'] for row in rows]==[1.5,1.,.5,0.]
    assert [row['screened_scalar_vacuum'] for row in rows]==[9/32,1/8,1/32,0.]
    assert all(row['eta_ratio_to_chain_root']==1 for row in rows)
    assert r['depth']['selected_n'] is None
    assert r['depth']['ancestor_flux']=='5/2'
    assert r['prediction']['ln_hierarchy'] is None


def test_volume_chain_telescopes_without_imposing_a_repeated_increment():
    a,b,c,d=sp.symbols('a b c d',real=True)
    assert sp.simplify(sp.exp(-d*(b-a)/2)*sp.exp(-d*(c-b)/2)-sp.exp(-d*(c-a)/2))==0
    assert volume_charge_ratio(.13,.13)==1
    assert volume_charge_ratio(.1,-.2)*volume_charge_ratio(-.2,.03)==pytest.approx(volume_charge_ratio(.1,.03))
    assert calculation()['volume_extension']['matched_start_to_minimum_ratio']>1


def test_reduction_is_not_a_dimensionless_unit_change_or_an_extra_nesting_step():
    k,V,e,s,scale=sp.symbols('k V e s scale',positive=True)
    eta4=(k/V)*(e*sp.sqrt(V)*sp.exp(-6*s))
    etaD=e*k**sp.Rational(4,7)
    volume=V*sp.exp(12*s)/k**sp.Rational(6,7)
    assert sp.simplify(eta4-etaD/sp.sqrt(volume))==0
    assert sp.simplify((k/scale**2)*(e*scale**2)-k*e)==0
    assert sp.simplify(eta4.subs(s,0)/etaD-1/sp.sqrt(volume.subs(s,0)))==0


def test_parke_ruler_does_not_use_historical_wall_or_curvature(monkeypatch):
    def forbidden():raise AssertionError('historical R,M,C preparation called')
    for owner in ('nsc_atlas_parke_nucleation','nsc_atlas_parent_scale_selector','nsc_atlas_inherited_matter'):
        monkeypatch.setattr('recursive_horizons.'+owner+'.exact_initial_geometry',forbidden)
    p=parke_ruler()
    ui,uo,t=F(7,32),F(21,32),F(1,5)
    inv=(uo-ui)**2/(9*t*t)+(ui+uo)/6+t*t/16
    assert 1/inv==F(57600,39169)
    assert p['L_P_over_L_unit']**2==pytest.approx(float(1/inv))
    assert p['late_volume_Lambda_L_P_squared']==pytest.approx(.3249667388660973)
    assert p['Parke_interior_Lambda_L_P_squared_exact']=='12600/39169'
    assert not p['FRW_reference_event_selected']
    assert not p['shared_volume_Parke_solution_derived']


def test_charge_target_keeps_the_factor_two_and_ruler_change_cancels():
    r=calculation();t=r['target_diagnostics'];p=r['parke_ruler']
    assert t['required_log_eta']==pytest.approx(t['required_log_g']/2-math.log(2))
    assert t['required_log_eta']==pytest.approx(-138.26042497625164)
    assert t['Parke_target_log_from_Lambda_obs']==pytest.approx(136.84391394458657)
    old=math.log(6.626077076099598/5)-t['old_ruler_target_log']
    new=p['log_hierarchy_intercept_eta']-t['Parke_target_log_transported_from_rounded_old']
    assert old==pytest.approx(new,abs=3e-14)
    # Changing a benchmark must not alter a physical map or Parke geometry.
    other=calculation(target_old_log_hierarchy=137,target_lambda=2e-52)
    for key in ('flux_step','volume_extension','dimensional_reduction','depth','parke_ruler'):
        assert r[key]==other[key]
