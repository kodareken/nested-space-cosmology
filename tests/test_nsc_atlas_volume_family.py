from fractions import Fraction as F
import json

import mpmath as mp
import pytest
import sympy as sp

from recursive_horizons.nsc_atlas_inherited_matter import exact_initial_geometry
from recursive_horizons.nsc_atlas_parent_scale_selector import matched_coefficients
from recursive_horizons.nsc_atlas_volume_family import (
    calculation,initial_stability_polynomial,inverse_matching,member,minimum_stable_dimension)


def test_general_stability_follows_by_differentiating_the_matched_action():
    d,j,s=sp.symbols('d j s',positive=True)
    rho,E=sp.Rational(9,32),sp.Rational(3,8)
    T=(j+3)**2/32
    B=d*(rho+T+E/2)/2
    A=B+rho-T
    V=A*sp.exp(-d*s)-B*sp.exp(-(d+2)*s)+T*sp.exp(-2*d*s)+E*sp.exp(-d*s/2)
    assert sp.simplify(32*sp.diff(V,s,2).subs(s,0)/d-initial_stability_polynomial(d,j))==0
    assert [minimum_stable_dimension(j) for j in range(1,8)]==[None,12,6,5,4,4,3]


def test_envelope_derivatives_and_root_structure():
    d,u,T,rho,E,x=sp.symbols('d u T rho E x',positive=True)
    S0=rho+T+E/2
    V=sp.exp(-u)*(rho-T+d*S0*(1-sp.exp(-2*u/d))/2+T*sp.exp(-u))
    wantedT=sp.exp(-u)*(d/2-1-d*sp.exp(-2*u/d)/2+sp.exp(-u))
    wantedd=sp.exp(-u)*S0*(1-(1+2*u/d)*sp.exp(-2*u/d))/2
    assert sp.simplify(sp.diff(V,T)-wantedT)==0
    assert sp.simplify(sp.diff(V,d)-wantedd)==0
    # F(1)<0 and strict convexity for d>2 establish the two-root brackets.
    B=d*S0/2;A=B+rho-T
    f=d*A-(d+2)*B*x+2*d*T*x**(d/2)
    assert sp.simplify(f.subs(x,1)+d*E/2)==0
    assert sp.simplify(sp.diff(f,x,2)-d*d*T*(d/2-1)*x**(d/2-2))==0


@pytest.mark.parametrize('d,j',[(10,2),(12,2),(16,2),(30,4),(6,3),(5,4),(4,5),(3,7)])
def test_general_roots_against_independent_high_precision_stationarity(d,j):
    row=member(d,j)
    with mp.workdps(70):
        c=matched_coefficients(d,j)
        number=lambda q:mp.mpf(q.numerator)/q.denominator
        A,B,T=(number(c[k]) for k in ('A','B','T'))
        V=lambda s:A*mp.exp(-d*s)-B*mp.exp(-(d+2)*s)+T*mp.exp(-2*d*s)
        root=mp.findroot(lambda s:mp.diff(V,s),(-mp.mpf('.08'),-mp.mpf('.02')))
        assert float(root)<0
        target=V(root)/number(exact_initial_geometry()['curvature'])
        assert row['Lambda_Lref0_squared']==pytest.approx(float(target),rel=3e-12)
        assert float(mp.diff(V,root,2))>0


def test_independently_given_matching_data_would_identify_a_pair_but_reference_is_circular():
    assert inverse_matching(7,F(15,2))==dict(T='25/32',d='12',j=2.)
    other=matched_coefficients(6,3)
    assert inverse_matching(other['A'],other['B'])==dict(T='9/8',d='6',j=3.)
    assert other['A']!=7 and other['B']!=F(15,2)


def test_stable_counterexamples_and_monotonic_envelope_are_not_a_selector():
    r=calculation()
    rows=r['witnesses']
    assert rows[0]['initial_homogeneous_stable'] is False
    assert all(row['initial_homogeneous_stable'] for row in rows[1:])
    assert rows[0]['Lambda_Lref0_squared']<rows[1]['Lambda_Lref0_squared']<rows[2]['Lambda_Lref0_squared']<rows[3]['Lambda_Lref0_squared']
    assert rows[4]['Lambda_Lref0_squared']>rows[1]['Lambda_Lref0_squared']
    assert not r['dimension_selected'] and not r['ancestor_depth_selected']
    assert r['new_field_evolutions']==0
    frontier=r['stable_integer_frontier']
    assert [(p['d'],p['j']) for p in frontier]==[(12,2),(6,3),(5,4),(4,5),(3,7)]
    assert all(p['coefficient']>frontier[0]['coefficient'] for p in frontier[1:])
    assert all(F(frontier[0]['exact_outward_enclosure'][1])<F(p['exact_outward_enclosure'][0]) for p in frontier[1:])
    for p in frontier:
        lo,hi=map(float,p['exact_outward_enclosure'])
        assert lo<p['coefficient']<hi
    assert r['family_energy_minimum']['unique']
    assert not r['family_energy_minimum']['physical_selection_derived']
    json.dumps(r)
