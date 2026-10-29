import math
from fractions import Fraction as F
import pytest
import sympy as s
from recursive_horizons.nsc_atlas_parke_nucleation import thin_wall,geometry_record,source_record


def test_exact_Parke_radius_unsquared_branches_and_mass_balance():
    R=240/s.sqrt(39169);bi=187/s.sqrt(39169);bo=163/s.sqrt(39169)
    Ui,Uo,k=s.Rational(7,32),s.Rational(21,32),s.Rational(1,10)
    assert s.simplify(bi-bo-k*R)==0
    assert s.simplify(bi**2+Ui*R*R/3-1)==0
    assert s.simplify(bo**2+Uo*R*R/3-1)==0
    assert s.simplify((Ui-Uo)*R**3/6+k*R*R*bi-k*k*R**3/2)==0
    r=geometry_record()
    assert r['bubble']['radius']==pytest.approx(float(R),rel=2e-15)
    assert r['bubble']['interior_Misner_Sharp_mass']>0
    assert abs(r['energy_balance_residual'])<1e-15
    assert r['action_control']['existing_coefficient']==pytest.approx(r['action_control']['independent_Parke_coefficient'],rel=1e-13)


def test_flat_gravity_limit_and_unit_dictionary():
    Ui,Uo,sigma=7/32,21/32,.2
    assert thin_wall(kappa=0)['radius']==pytest.approx(3*sigma/(Uo-Ui))
    r=geometry_record();c=r['legacy_comparison']
    assert c['new_R_over_old_Lref']*c['old_Lref_over_Lunit']==pytest.approx(r['bubble']['radius'])
    assert c['old_R_over_old_Lref']*c['old_Lref_over_Lunit']==pytest.approx(1.9)
    assert r['bubble']['radius']<1.9


def test_supplied_snapshot_is_not_a_pair_of_scalar_vacua():
    y,q=s.symbols('y q',real=True);epsilon=s.Rational(1,4)
    U=epsilon*y*y/2+(q+s.sqrt(1-epsilon)*y)**2/2
    y0=-s.sqrt(1-epsilon)/2
    assert s.simplify(U.subs({y:y0,q:1}))==s.Rational(7,32)
    assert s.simplify(U.subs({y:y0,q:s.Rational(3,2)}))==s.Rational(21,32)
    assert s.simplify(s.diff(U,y).subs({y:y0,q:1}))==s.sqrt(3)/4
    assert s.simplify(s.diff(U,y).subs({y:y0,q:s.Rational(3,2)}))==s.sqrt(3)/2


def test_empty_matter_invariant_does_not_assume_an_exit_value():
    N,phi,Ei,H,v=s.symbols('N phi Ei H v',real=True)
    E=Ei*s.exp(-3*N-6*phi)
    assert s.simplify(s.diff(E,N)*H+s.diff(E,phi)*v+3*H*E+6*v*E)==0
    assert E.subs(Ei,0)==0
    chi,pi,omega=s.symbols('chi pi omega',real=True)
    assert (-3*H*pi-omega**2*chi).subs({chi:0,pi:0})==0
    r=source_record()
    assert not r['empty_start']['reference_event_reached']
    assert r['downstream']['SN_DESI_refit'] is None


def test_energy_budget_allows_parent_delivery_but_not_the_stated_empty_roll():
    T,E,V,Q,H=s.symbols('T E V Q H',real=True)
    E_dot=-3*H*E+Q;scalar_dot=-6*H*T-Q
    assert s.simplify(E_dot+scalar_dot+3*H*(E+2*T))==0
    r=source_record();b=r['energy_budget']
    assert b['literal_child_total_density']<.375
    assert b['literal_Ui_excess']<0
    assert 0<b['available_excess_from_computed_well_barrier']<.375
    assert b['parent_Uo_excess']>.375
    assert b['minimum_required_rest_exit_density']==pytest.approx(.595983024212607)
    assert not b['necessary_s_is_a_CDL_turning_point']
    assert not r['perturbative_reheating']['pair_decay_open']
