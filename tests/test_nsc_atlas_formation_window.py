from pathlib import Path
import hashlib
import json
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_formation_window import FormationBackground,transfer,dust_by_quadrature,density_phase,elementary_gain_bound,compact_packet,packet_tail_fraction


def test_independent_dust_basis_and_positive_future_bound():
    row,raw=transfer(rtol=2e-12)
    np.testing.assert_allclose(row['future_density_response'],dust_by_quadrature(),rtol=2e-11,atol=2e-12)
    assert np.min(np.diff(raw['response'][:,0,0]))>=-1e-12
    assert row['future_density_response'][0]<elementary_gain_bound()['density_gain_bound']<1.4


def test_phase_equation_and_determinant():
    for k in (2000.,8500.,30000.):
        row,_=transfer(k,initial_scale=.1)
        assert row['future_density_response'][0]==pytest.approx(density_phase(k,initial_scale=.1),abs=2e-8)
        assert row['determinant_relative_error']<1e-8
        assert row['sampled_max_density_response']<=dust_by_quadrature(.1)[0]+1e-7


def test_no_gravity_exact_oscillator_over_infinite_expansion():
    b=FormationBackground(vacuum=3.,excitation=0.,curvature=0.,frequency=100.)
    k=20.;row,raw=transfer(k,background=b)
    u=raw['u'];eta=(1-u*u)/2;nu=k*k/(2*b.frequency)
    np.testing.assert_allclose(raw['response'][:,0,0],np.cos(nu*eta),atol=2e-10,rtol=2e-10)
    np.testing.assert_allclose(raw['response'][:,0,1],np.sin(nu*eta)/nu,atol=2e-10,rtol=2e-10)


def test_reciprocal_endpoint_and_phase_are_not_assumed():
    row,raw=transfer(8414.774665477327,initial_scale=.01)
    assert raw['u'][-1]==0 and np.isfinite(raw['response']).all()
    assert row['future_density_response'][0]>50
    assert row['future_density_response'][1]>40
    # A phase/current seed has an independently nonzero response.
    assert row['initial_k_over_frequency_A']<.001
    assert row['initial_H_over_frequency']<1e-6


def test_retained_formation_record():
    root=Path(__file__).resolve().parents[1];path=root/'results/summaries/nsc-atlas-formation-window-v1.json'
    if not path.exists():pytest.skip('producer not yet run')
    r=json.loads(path.read_text(encoding='utf-8'))
    assert not r['numerical_validation_errors'] and not r['budget_exhausted']
    assert r['symbolic_identities']['dust_H_mode']=='0'
    for name,h in {**r['source_sha256'],**r['evidence_sha256']}.items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==h,name


def test_compact_seed_parseval_and_tail_bound():
    from scipy.special import roots_legendre,spherical_jn,beta
    from scipy.integrate import quad
    q,w=roots_legendre(256);q=(q+1)*20;w=w*20
    norm=4*np.pi*16*beta(2.5,7)/beta(1.5,5)**2
    retained=np.sum(w*q**4*(10395*spherical_jn(5,q)/q**5)**2)/norm
    assert 0<1-retained<packet_tail_fraction(40)
    # Independent real-space integral of the gradient squared.
    gradnorm=4*np.pi*quad(lambda r:r*r*(-8*r*(1-r*r)**3)**2,0,1,epsabs=1e-12)[0]
    G0=4*np.pi*quad(lambda r:r*r*(1-r*r)**4,0,1,epsabs=1e-12)[0]
    assert norm==pytest.approx(2*np.pi**2*gradnorm/G0**2,rel=2e-13)
    assert packet_tail_fraction(60)<packet_tail_fraction(40)


def test_compact_packet_late_epoch_respects_uniform_bound():
    row,raw=compact_packet(.0004815623798469198,initial_scale=1.,order=96,qmax=30.)
    assert max(abs(raw['central']))<dust_by_quadrature()[0]
    assert abs(row['initial_normalization']-1)<row['initial_omitted_spectral_weight_bound']
