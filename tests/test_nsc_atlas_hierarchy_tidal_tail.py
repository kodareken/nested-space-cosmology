import numpy as np
import pytest

from recursive_horizons.nsc_atlas_hierarchy_tidal_tail import tidal_kinetic


def test_harmonic_boundary_rank_update_matches_dense_self_adjoint_cayley():
    n=30;dr=.1;h=.02;ell=2;r=dr*np.arange(1,n+1)
    H=np.diag(1/dr**2+ell*(ell+1)/(2*r*r))
    H+=np.diag(np.full(n-1,-.5/dr**2),1)+np.diag(np.full(n-1,-.5/dr**2),-1)
    H[-1,-1]-=(n/(n+1))**ell/(2*dr**2)
    rng=np.random.default_rng(4);q=rng.normal(size=(n,3))+1j*rng.normal(size=(n,3))
    exact=np.linalg.solve(np.eye(n)+.5j*h*H,(np.eye(n)-.5j*h*H)@q)
    actual=tidal_kinetic(q,dr,h)
    np.testing.assert_allclose(actual,exact,rtol=3e-14,atol=3e-14)
    assert np.linalg.norm(actual)==pytest.approx(np.linalg.norm(q),rel=3e-15)


def test_dirichlet_clipping_breaks_the_harmonic_tidal_phase():
    n=999;L=10.;dr=L/(n+1);r=dr*np.arange(1,n+1);q=r**-2
    # At the outermost interior point the exact kinetic operator on
    # r^-2 is zero. The missing exterior phase is the Dirichlet defect.
    common=(1/dr**2+3/r[-1]**2)*q[-1]-.5*q[-2]/dr**2
    clipped=common
    matched=common-.5*(r[-1]/L)**2*q[-1]/dr**2
    assert abs(clipped)>1
    assert abs(matched)<1e-7
    assert abs(clipped)>1e8*abs(matched)


def test_angular_addition_identity_fixes_the_pointwise_potential_bound():
    rng=np.random.default_rng(12)
    for _ in range(30):
        n=rng.normal(size=3);n/=np.linalg.norm(n);x,y,z=n
        Y=np.array([np.sqrt(5)*(3*z*z-1)/2,np.sqrt(15)*x*z,np.sqrt(15)*y*z,
            np.sqrt(15)*(x*x-y*y)/2,np.sqrt(15)*x*y])
        assert Y@Y==pytest.approx(5,abs=1e-14)


def test_retained_tail_fixes_domain_sensitivity_without_claiming_resolved_amplitude():
    import json
    import hashlib
    from pathlib import Path
    root=Path(__file__).resolve().parents[1]
    r=json.loads((root/'results/summaries/nsc-atlas-hierarchy-tidal-tail-v1.json').read_text(encoding='utf-8'))
    assert not r['numerical_validation_errors'] and not r['budget_exhausted']
    assert all(r['source_replay_exact_matches'].values())
    old=r['old_Dirichlet_probe']['by_resolution']
    assert abs(old['medium']-old['outer_control'])>100*r['outer_domain_change']
    assert not r['outward_probe_resolved']
    assert r['numerical_comparison_envelope']>r['final_potential_RMS_coefficient']
    for name,digest in r['source_sha256'].items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,name
