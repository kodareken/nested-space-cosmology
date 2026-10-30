import math
from fractions import Fraction as F
import mpmath as mp
import pytest
import sympy as sp
from recursive_horizons.nsc_atlas_contrast_hierarchy import calculation,length_data,sector_minimum


def test_constant_volume_frame_change_preserves_horizon_planck_ratio():
    s,d,k,V=sp.symbols('s d k V',positive=True)
    omega=sp.exp(-d*s/2)
    # sqrt(-g_E)*R_E = exp(d*s)*sqrt(-g_J)*R_J at constant s.
    ricci_factor=sp.exp(2*d*s)*sp.exp(-d*s)
    kJ=k/ricci_factor
    lamE=k*V;lamJ=lamE/omega**2
    assert sp.simplify(kJ-k*omega**2)==0
    assert sp.simplify(sp.sqrt(kJ/(8*sp.pi))-omega*sp.sqrt(k/(8*sp.pi)))==0
    assert sp.simplify(1/(lamJ*kJ)-1/(lamE*k))==0
    for g in (1.,.25,1e-120):
        r=length_data(7/32,-.035,g)
        assert r['L_H_E']/r['Planck_E']==pytest.approx(r['L_H_J']/r['Planck_J'])
        parent=length_data(21/32,.07,g)
        assert r['log_ratio']-parent['log_ratio']==pytest.approx(math.log(3)/2)
        assert r['Planck_J']/parent['Planck_J']==pytest.approx(math.exp(-6*(-.035-.07)))


def test_parke_contrast_cancels_from_child_vacuum_when_inside_density_fixed():
    ui,uo,t,k=sp.symbols('ui uo t k',positive=True)
    inv=(uo-ui)**2/(9*t*t)+k*(ui+uo)/6+k*k*t*t/16
    recovered=3*(inv-((uo-ui)/(3*t)+k*t/4)**2)
    assert sp.simplify(recovered-k*ui)==0
    r=calculation()
    assert r['contrast']['recovered_child_Lambda_hat']=='7/32'
    assert not r['contrast']['volume_endpoint_difference_selected']


def test_exact_vacuum_sign_threshold_is_not_a_floating_point_decision():
    A,B,T,x=sp.symbols('A B T x',positive=True)
    V=x**6*(A-B*x+T*x**6)
    a_stationary=sp.solve(sp.diff(V,x),A)[0]
    assert sp.simplify(V.subs(A,a_stationary)-x**7*(B/6-T*x**5))==0
    r=calculation()
    for row in r['fixed_coupling_j2_sectors']+r['fixed_original_sheet_control']:
        a,b,t=(F(row[key]) for key in ('A','B','T'))
        assert ((6*a/(5*b))**5>b/(6*t))==row['positive_vacuum_by_exact_comparison']
        assert row['stiffness']>0
    assert [v['positive_vacuum_by_exact_comparison'] for v in r['fixed_coupling_j2_sectors']]==[True,False,False]


def test_sector_minima_against_independent_high_precision_roots():
    with mp.workdps(65):
        for q in (1.5,1.,.5):
            r=sector_minimum(q)
            a,b,t=(mp.mpf(F(r[key]).numerator)/F(r[key]).denominator for key in ('A','B','T'))
            x=mp.findroot(lambda x:6*a-7*b*x+12*t*x**6,mp.mpf(r['x_min']))
            v=x**6*(a-b*x+t*x**6)
            assert float(v)==pytest.approx(r['Lambda_hat_min'],rel=2e-13,abs=2e-13)
            assert r['UV_bare']=='215/32' and r['B']=='15/2'
    with pytest.raises(ValueError):length_data(-.04,0,1)


def test_chain_growth_and_absolute_normalization_are_distinct():
    r=calculation();rows=r['auxiliary_Parke_snapshot_chain']
    assert rows[1]['log_ratio_intercept']-rows[0]['log_ratio_intercept']==pytest.approx(math.log(3)/2)
    assert rows[2]['log_ratio_intercept']-rows[1]['log_ratio_intercept']==pytest.approx(math.log(7)/2)
    a,b=length_data(7/32,0,1),length_data(7/32,0,.25)
    assert b['log_ratio']-a['log_ratio']==pytest.approx(math.log(2))
    assert r['prediction']['log_hierarchy'] is None
    assert r['prediction']['percent_miss'] is None


def test_observed_comparison_uses_the_same_lengths_and_does_not_set_parameters():
    r=calculation();other=calculation(target_lambda=2e-52)
    assert r['observed_benchmark']['log_L_H_over_Planck']==pytest.approx(139.93430661214543)
    for key in ('contrast','fixed_coupling_j2_sectors','fixed_original_sheet_control','frame_law'):
        assert r[key]==other[key]
