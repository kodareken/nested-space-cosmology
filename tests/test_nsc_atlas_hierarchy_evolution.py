import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import erf

from recursive_horizons.nsc_atlas_hierarchy_evolution import (
    background,kinetic_step,potential,cut_flux)
from recursive_horizons.nsc_atlas_scale_hierarchy import ScaleHierarchy


def test_free_expanding_gaussian_uses_integrated_kinetic_time():
    model=ScaleHierarchy(1e9,.001,.5);parent=background(model,100.)
    A,Q=parent.sol(100.);L=256.;n=8191;dr=L/(n+1);r=dr*np.arange(1,n+1)
    width=3.;k=np.pi*np.arange(1,n+1)/L
    u=r*np.exp(-r*r/(2*width*width))
    evolved=kinetic_step(u,k,Q)/r
    exact=(width*width/(width*width+1j*Q))**1.5*np.exp(-r*r/(2*(width*width+1j*Q)))
    assert A>1 and Q<100.
    np.testing.assert_allclose(evolved[r<40],exact[r<40],atol=3e-13,rtol=3e-12)


def test_spherical_potential_recovers_an_independent_gaussian_integral():
    errors=[]
    for n in (2047,4095):
        L=32.;dr=L/(n+1);r=dr*np.arange(1,n+1);a=2.;A=1.3
        u=r*np.exp(-r*r/(2*a*a));V,_=potential(u,r,dr,A,0.)
        exact=-math.sqrt(math.pi)*a**3*erf(r/a)/(8*A*r)
        errors.append(float(np.max(abs(V-exact))))
    assert errors[1]<errors[0]/3.9
    assert errors[1]<3e-6


def test_moving_physical_cut_includes_the_expansion_flux():
    model=ScaleHierarchy(1e9,.001,.5);parent=background(model,30.)
    tau=20.;A,Q=parent.sol(tau);omega_scale=1000.;b=.0006123724356957944
    h=math.sqrt(model.vacuum/3+model.excitation/(3*A**3)-model.curvature/A**2)/omega_scale
    X=3.2;L=128.;n=16383;dr=L/(n+1);r=dr*np.arange(1,n+1);width=3.
    def wave(x,q):return (width**2/(width**2+1j*q))**1.5*np.exp(-x*x/(2*(width**2+1j*q)))
    def mass(t):
        a,q=parent.sol(t)
        return quad(lambda x:x*x*(abs(wave(x,q))**2+2*b*wave(x,q).real),0,X/a,epsabs=1e-12)[0]
    u=r*wave(r,Q);flux=cut_flux(u,r,dr,A,h,b,X)
    derivative=(mass(tau+.0001)-mass(tau-.0001))/.0002
    assert derivative==pytest.approx(-flux,rel=2e-7,abs=2e-9)


def test_homogeneous_envelope_has_no_peculiar_gravity_or_kinetic_change():
    r=np.linspace(.1,20,200);u=np.zeros_like(r,dtype=complex)
    V,radial=potential(u,r,.1,1.2,.03)
    assert np.max(abs(V))==0 and np.max(abs(radial))==0
    np.testing.assert_array_equal(kinetic_step(u,np.pi*np.arange(1,201)/20.1,4.),u)


def test_real_field_background_and_geometric_mass_normalizations():
    from scripts.check_nsc_atlas_hierarchy_evolution import symbolic_normalizations
    assert set(symbolic_normalizations().values())=={'0'}


def test_retained_envelope_and_full_equation_comparison_keep_their_domains():
    root=Path(__file__).resolve().parents[1]
    r=json.loads((root/'results/summaries/nsc-atlas-hierarchy-evolution-v1.json').read_text(encoding='utf-8'))
    assert not r['numerical_validation_errors'] and not r['budget_exhausted']
    for name,digest in r['source_sha256'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,name
    comparison=r['short_full_equation_comparison']
    assert abs(comparison['small']['difference'])<abs(comparison['large']['difference'])
    assert r['runs']['fine']['config']['end']==600.
    assert r['runs']['comparison_GR_small']['end_slow_time']<6.
