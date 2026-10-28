import math
import pytest
import sympy as s

from recursive_horizons.nsc_atlas_scale_selection import calculation


def test_SI_conversion_and_target_independence_of_conditional_predictions():
    result=calculation();changed=calculation(target=2*result['target_Lambda_m_inverse2'])
    constants=result['SI_conventions']
    assert constants['lPlanck_metres']==pytest.approx(1.6162550244e-35,rel=1e-10,abs=0)
    assert constants['kappa_natural_metres_squared']==pytest.approx(8*math.pi*constants['lPlanck_metres']**2,rel=2e-15,abs=0)
    top=result['route_ladder'][4]['computed'];other=changed['route_ladder'][4]['computed']
    assert top['winding_1_Lambda_m_inverse2']==other['winding_1_Lambda_m_inverse2']
    assert top['winding_1_Lambda_m_inverse2']*constants['lPlanck_metres']**2==pytest.approx(1.,rel=2e-15)
    assert not any(row['survives'] for row in result['route_ladder'])
    assert result['inverse_target_diagnostic']['Lref_metres']!=changed['inverse_target_diagnostic']['Lref_metres']


def test_same_normalized_inheritance_has_distinct_SI_realizations():
    rows=calculation()['scale_nonidentifiability_witnesses']
    for row in rows:
        assert row['recovered_dimensionless_charge']==pytest.approx(.5)
        assert row['recovered_dimensionless_flux']==pytest.approx(1.5)
        assert row['electric_flux_integer']==pytest.approx(3.)
    assert rows[0]['Lambda_m_inverse2']==pytest.approx(4*rows[1]['Lambda_m_inverse2'])


def test_formal_S4_period_does_not_determine_the_winding():
    eps,q,e,k,L,n=s.symbols('eps q e k L n',positive=True)
    lam=eps*q*q/2
    V4=24*s.pi**2/lam**2
    period=e*eps*q*V4*L*L/k
    selected=s.solve(s.Eq(period,2*s.pi*n),L*L)[0]
    # k=8*pi*lP^2; Lambda=lam/L^2. The result still contains n.
    product=s.simplify(lam/selected*k/(8*s.pi))
    assert s.simplify(product-3*e/(n*q))==0
    assert n in product.free_symbols
