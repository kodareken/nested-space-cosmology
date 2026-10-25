import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad

from recursive_horizons.nsc_atlas_compensation_response import CompensationResponse,sinc_minus_one,nodes,point_core_limit,retained_evolution_response


@pytest.fixture(scope='module')
def profile():return CompensationResponse(order=160)


def test_compensated_shell_and_taylor_remainder():
    # Unit point mass minus a uniform unit sphere: F=1-3*j1(K)/K.
    from scipy.special import spherical_jn
    s,q=nodes(160,0,1)
    for K in (.001,.1,1.,4.,8.):
        computed=-np.sum(q*3*s*s*sinc_minus_one(K*s))
        if K>.01:assert computed==pytest.approx(1-3*spherical_jn(1,K)/K,abs=2e-13)
        assert abs(computed-K*K/10)<=K**4/280+1e-14
    assert sinc_minus_one(0)==0


def test_gradient_is_exact_positive_missing_mass(profile):
    s,q,gr,env,grad=profile.quadrature
    assert np.min(grad)>=0
    np.testing.assert_allclose(gr-env,grad,rtol=2e-6,atol=2e-13)
    assert abs(np.sum(q*gr))<1e-9
    assert np.sum(q*env)==pytest.approx(-np.sum(q*grad),abs=1e-9)
    K=np.geomspace(.0001,100,80)
    assert np.max(abs(profile.transform(K)-profile.transform(K,channel='envelope')))<=np.sum(q*grad)+1e-13
    assert profile.transform(0,channel='envelope')<0 # never silently compensated


def test_independent_mass_primitive_and_low_k_bound(profile):
    row=profile.record();K=np.geomspace(.001,20,70)
    np.testing.assert_allclose(profile.transform(K),profile.by_parts(K),atol=1e-10,rtol=1e-10)
    low=K[K<=1];residual=abs(profile.transform(low)-row['mass_moments_normalized'][0]-row['low_K_quadratic_coefficient']*low**2)
    assert np.all(residual<=row['low_K_fourth_remainder_bound']*low**4+1e-13)


def test_projection_by_independent_shell_fraction(profile):
    for b in (.1,.5,.9):
        # Independent radial shell integral; Abel code uses line of sight.
        f=lambda s:float(profile.densities(s)[0])*(1 if s<=b else 1-np.sqrt(1-b*b/(s*s)))
        cut=profile.inner/profile.outer
        mass=sum(quad(f,a,z,epsabs=2e-10,points=None)[0] for a,z in ((0,.75*cut),(.75*cut,cut),(cut,b),(b,1)))
        got=profile.projected(b)
        assert got['cylindrical_mass']==pytest.approx(mass,abs=2e-9)
        assert got['surface']<0 and got['mean_minus_local']>0
    assert profile.projected(1)['surface']==0
    assert profile.projected(2,channel='envelope')['cylindrical_mass']<0


def test_turnover_and_limit_are_preparation_dependent(profile):
    K=profile.half_power()
    assert (profile.transform(K)/profile.transform(K,core_only=True))**2==pytest.approx(.5,abs=1e-12)
    values=[point_core_limit(a) for a in (.25,.5,.75)]
    assert values[1]['Lambda_Rd_cubed_over_mass']==pytest.approx(200783583/10953472,rel=1e-12)
    assert len({round(v['Lambda_over_core_mass_k_half_cubed'],4) for v in values})==3


def test_physical_wavenumber_changes_under_pure_expansion():
    r=np.arange(1,4096)*80/4096;b=2.
    core=(r<5).astype(float);shell=((r>25)&(r<30)).astype(float)
    density=.1*(core-shell*np.sum(r*r*core)/np.sum(r*r*shell))
    u=r*(np.sqrt(b*b+density)-b)+0j
    first=retained_evolution_response(r,u,u,1.,b,30.)
    stretched=retained_evolution_response(r,u,u,1.1,b,30.)
    assert first['half_power_K_in_original_Rd_units']==pytest.approx(
        1.1*stretched['half_power_K_in_original_Rd_units'],rel=1e-11)
    assert first['total_excess_over_initial_core']==pytest.approx(stretched['total_excess_over_initial_core'],abs=1e-14)


def test_retained_record_has_actual_producers_and_scope():
    root=Path(__file__).resolve().parents[1];path=root/'results/summaries/nsc-atlas-compensation-response-v1.json'
    if not path.exists():pytest.skip('producer record not created yet')
    r=json.loads(path.read_text(encoding='utf-8'))
    assert not r['numerical_validation_errors'] and not r['budget_exhausted']
    for name,h in {**r['source_sha256'],**r['retained_evidence_sha256']}.items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,name
    assert max(abs(v) for v in r['retained_geometry_comparison'].values())<1e-10
    evolved=r['evolved_spatial_response']
    assert hashlib.sha256((root/evolved['source_record']).read_bytes()).hexdigest()==evolved['source_evidence_sha256']
    assert evolved['numerical_comparison_ranges']['half_power_K_in_original_Rd_units']<4e-5
