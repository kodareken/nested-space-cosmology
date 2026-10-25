from pathlib import Path
import hashlib
import json

import numpy as np
import pytest
import sympy as sp
from scipy.integrate import quad
from scipy.special import spherical_jn,roots_legendre,beta

from recursive_horizons.nsc_atlas_seed_collapse import seed_shape,seed_mass_primitive,initial_constraint,boundary_potential,shell_positions


def test_exact_compact_convolution_mass_and_derivative():
    r=sp.symbols('r');H=(r-2)**10*(63*r**9+1260*r**8+11200*r**7+57680*r**6+187200*r**5+388032*r**4+499200*r**3+395520*r**2+179200*r+35840)/581981400
    S=sp.factor(-(sp.diff(H,r,2)+2*sp.diff(H,r)/r)/(sp.Rational(262144,255255)))
    assert sp.integrate(r*r*S,(r,0,2))==0
    f=sp.lambdify(r,S,'numpy');fp=sp.lambdify(r,sp.diff(S,r),'numpy')
    x=np.linspace(.001,1.99,51);actual,derivative=seed_shape(x)
    np.testing.assert_allclose(actual,f(x),rtol=1e-11,atol=1e-12)
    np.testing.assert_allclose(derivative,fp(x),rtol=1e-10,atol=1e-11)
    assert seed_shape(0)[0]==1 and seed_shape(2)[0]==0
    for end in (.3,.8,1.5,2.):
        assert quad(lambda q:q*q*seed_shape(q)[0],0,end,epsabs=1e-12)[0]==pytest.approx(seed_mass_primitive(end),abs=2e-12)


def test_real_space_seed_matches_previous_fourier_preparation():
    q,w=roots_legendre(320);q=(q+1)*30;w=w*30
    I=64*np.pi*beta(2.5,7)/beta(1.5,5)**2
    weight=w*q**4*(10395*spherical_jn(5,q)/q**5)**2/I
    for r in (0,.2,.7,1.,1.5,2.1):
        assert np.sum(weight*np.sinc(q*r/np.pi))==pytest.approx(seed_shape(r)[0],abs=3e-8)


def test_constraint_correction_is_small_and_independently_closes():
    for amplitude in (.05,.1):
        row,field=initial_constraint(amplitude)
        assert 0<row['correction']<1e-9
        assert abs(row['corrected_GR_endpoint_norm'])<1e-14
        assert abs(row['proper_distance_endpoint_norm'])<1e-13
        assert row['minimum_B']>.9999
        assert field(0)[1]**2==pytest.approx(1+amplitude)
        assert field(2)[1]==1 and field(2)[2]==0


def test_phase_matched_poisson_eigenmode_and_symmetric_kernel():
    gaps=[]
    for n in (1023,2047):
        L=8.;dr=L/(n+1);r=dr*np.arange(1,n+1);k=np.pi/L
        delta=.01*np.sinc(k*r/np.pi);u=r*(np.sqrt(1+delta)-1)
        V,_=boundary_potential(u,r,dr,L)
        gaps.append(float(np.max(abs(V+delta/(2*k*k)))))
    assert gaps[1]<gaps[0]/3.9
    assert gaps[1]<3e-8


def test_volume_shell_reconstruction_keeps_uniform_density_exact():
    r=np.arange(1,101)*.01;density=np.full_like(r,7.)
    radii=np.array([.003,.12,.351]);targets=7*radii**3/3
    got,_=shell_positions(r,density,targets,7.)
    np.testing.assert_allclose(got,radii,atol=2e-13)
    _,cum=shell_positions(r,density,[],7.)
    got,_=shell_positions(r,density,cum[[5,50,99]],7.+1e-12)
    np.testing.assert_allclose(got,r[[4,49,98]],atol=2e-13)


def test_saved_result_source_and_finite_claim():
    root=Path(__file__).resolve().parents[1];p=root/'results/summaries/nsc-atlas-seed-collapse-v1.json'
    if not p.exists():pytest.skip('producer not yet run')
    r=json.loads(p.read_text(encoding='utf-8'))
    assert not r['budget_exhausted'] and not r['numerical_validation_errors']
    assert r['finite_formation_assessment']['strong_finite_contraction']
    assert r['finite_formation_assessment']['weak_sampled_shells_still_expanding']
    assert r['finite_formation_assessment']['no_gravity_sampled_shells_still_expanding']
    assert r['runs']['no_gravity']['final_physical_half_height_radius'] is None
    for name,digest in {**r['source_sha256'],**r['evidence_sha256']}.items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,name
