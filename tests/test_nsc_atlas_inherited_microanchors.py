import pytest
import sympy as s
from check_nsc_atlas_inherited_microanchors import calculation


def test_SI_coupling_ratios_are_inherited_but_parent_scale_is_not_selected():
    k,L,e0,sigma0,w0=s.symbols('k L ehat sigmahat what',positive=True)
    e=e0/(s.sqrt(k)*L);sigma=sigma0/(k*L);omega=w0/L
    assert s.simplify(e/sigma-e0*s.sqrt(k)/sigma0)==0
    assert s.simplify(omega/e-w0*s.sqrt(k)/e0)==0
    r=calculation();anchors=r['SI_coupling_ratio_lengths']
    for row in r['dimensional_completion_witnesses']:
        for key,value in anchors.items():
            assert row[key]==pytest.approx(value,rel=3e-15,abs=0)
    assert r['dimensional_completion_witnesses'][0]['parent_Lambda_m_inverse2']!=r['dimensional_completion_witnesses'][1]['parent_Lambda_m_inverse2']
    assert anchors['charge_over_tension_metres']>0
    assert not r['Lambda_obs_derived']
