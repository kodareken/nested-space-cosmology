import math
import pytest
import sympy as s
from recursive_horizons.nsc_atlas_si_hierarchy import calculation,conditional_at_g,units,A0,VSTAR,E_HAT,SIGMA_HAT


def test_units_and_geometric_identity_do_not_depend_on_the_tunnelling_guess():
    g,k,a,e,t,V=s.symbols('g k a e t V',positive=True)
    L=s.sqrt(k/g);ell=e*s.sqrt(k)/t;hierarchy=a*L/ell;lam=V/L**2
    assert s.simplify(hierarchy-a*t/(e*s.sqrt(g)))==0
    assert s.simplify(lam*ell**2-V*(e/t)**2*g)==0
    assert s.simplify(hierarchy**2*lam*ell**2-a*a*V)==0
    assert units()['ell_e_sigma_m']==pytest.approx(2.025675271617439e-34,rel=1e-14,abs=0)


def test_distinct_g_completions_preserve_all_supplied_normalized_charge_data():
    a,b=[conditional_at_g(g) for g in (1.,.25)]
    k=units()['kappa_nat_m2'];ell=units()['ell_e_sigma_m']
    for row in (a,b):
        assert math.sqrt(k)*row['L_unit_m']*row['canonical_charge']==pytest.approx(E_HAT)
        assert k*row['L_unit_m']*row['canonical_tension']==pytest.approx(SIGMA_HAT)
        assert row['canonical_charge']/row['canonical_tension']==pytest.approx(ell,rel=1e-14,abs=0)
        assert row['hierarchy']**2*row['Lambda_star_ell_squared']==pytest.approx(A0*A0*VSTAR)
    assert b['hierarchy']==pytest.approx(2*a['hierarchy'])
    assert b['Lambda_star_m_inverse2']==pytest.approx(a['Lambda_star_m_inverse2']/4)
    assert b['B_HM']==pytest.approx(4*a['B_HM'])


def test_exponential_inversions_are_not_geometric_predictions():
    r=calculation();logtarget=r['targets']['log_hierarchy']
    for row in r['exponential_inverse_diagnostics']:
        g=row['required_g_using_full_target_log'];p=row['p']
        assert row['same_geometry']['B_HM']/p==pytest.approx(logtarget)
        assert row['same_geometry']['hierarchy']<4
        assert row['same_geometry_hierarchy_over_target']<2e-60
        assert row['same_geometry_Lambda_over_observed']>1e119
    assert r['selected_g'] is None and r['selected_hierarchy_exponent'] is None
    assert r['prediction_miss_percent'] is None


def test_target_use_is_isolated_and_legacy_vacuum_is_not_substituted():
    a=calculation();b=calculation(target_hierarchy=2e60,target_lambda=2e-52)
    assert a['forward_identities']==b['forward_identities']
    assert a['dimensional_completions']==b['dimensional_completions']
    inv=a['target_inverse_consistency']
    assert inv['g_if_hierarchy_is_inserted']==pytest.approx(2.54538187565e-120,rel=1e-10,abs=0)
    assert inv['Lambda_miss_percent_if_hierarchy_is_inserted']==pytest.approx(-21.4282580474,abs=1e-8)
    assert a['targets']['Lambda_obs_ell_squared']==pytest.approx(4.474304077708618e-120,rel=1e-14,abs=0)


def test_four_volume_length_and_fixed_volume_time_have_different_powers():
    B,A,V3=s.symbols('B A V3',positive=True)
    rate=A*s.exp(-B)
    assert s.simplify(rate**(-s.Rational(1,4))-A**(-s.Rational(1,4))*s.exp(B/4))==0
    assert s.simplify(1/(rate*V3)-s.exp(B)/(A*V3))==0
    with pytest.raises(ValueError):conditional_at_g(0)
