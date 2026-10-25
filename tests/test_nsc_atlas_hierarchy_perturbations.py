import math
import numpy as np
import pytest

from recursive_horizons.nsc_atlas_hierarchy_perturbations import (
    multipole_potential,angular_kinetic,seed_basis,diagnostics)
from recursive_horizons.nsc_atlas_hierarchy_evolution import potential


def test_quadrupole_green_function_matches_independent_radial_laplacian():
    errors=[]
    for n in (511,1023,2047):
        L=32.;dr=L/(n+1);r=dr*np.arange(1,n+1);ell=2;A=1.2;s=2.
        exact=r**ell*np.exp(-r*r/(2*s*s))
        source=(r*r/s**4-(2*ell+3)/s**2)*exact
        # psi=1, eta=A*Delta_l exact gives delta_n=2A*Delta_l exact.
        q=(r*A*source)[:,None]
        numeric=multipole_potential(np.ones_like(r),q,r,dr,A,ell)[:,0]
        errors.append(np.max(abs(numeric-exact)))
    assert errors[1]<errors[0]/3.9 and errors[2]<errors[1]/3.9


def test_centrifugal_kinetic_is_unitary_and_recovers_solid_harmonic_gaussian():
    errors=[]
    for n,dt in ((511,.002),(1023,.001),(2047,.0005)):
        L=32.;dr=L/(n+1);r=dr*np.arange(1,n+1);s=2.;ell=2;t=1.
        q=np.asarray((r**(ell+1)*np.exp(-r*r/(2*s*s)))[:,None],dtype=complex)
        norm=np.linalg.norm(q)
        for _ in range(round(t/dt)):q=angular_kinetic(q,dr,dt,ell)
        exact=r**(ell+1)*(s*s/(s*s+1j*t))**(ell+1.5)*np.exp(-r*r/(2*(s*s+1j*t)))
        errors.append(np.linalg.norm(q[:,0]-exact)/np.linalg.norm(exact))
        assert np.linalg.norm(q)==pytest.approx(norm,rel=3e-13)
    assert errors[1]<errors[0]/3.9 and errors[2]<errors[1]/3.9 and errors[2]<1e-4


def test_potential_tangent_is_derivative_of_nonlinear_radial_flow():
    L=32.;n=1023;dr=L/(n+1);r=dr*np.arange(1,n+1);A=1.1;b=.03;h=.02
    u=r*np.exp(-r*r/8)*(1+.1j);psi=b+u/r;q=(r**3*np.exp(-r*r/12)*(1-.3j))[:,None]
    V,_=potential(u,r,dr,A,b);dv=multipole_potential(psi,q,r,dr,A,0)
    exact=np.exp(-1j*h*V)[:,None]*(q-1j*h*(r*psi)[:,None]*dv)
    errors=[]
    for eps in (1e-4,5e-5):
        values=[]
        for sign in (1,-1):
            v=u+sign*eps*q[:,0];P,_=potential(v,r,dr,A,b)
            values.append((v+b*r)*np.exp(-1j*h*P)-b*r)
        finite=(values[0]-values[1])/(2*eps)
        errors.append(np.linalg.norm(finite-exact[:,0])/np.linalg.norm(exact))
    assert errors[1]<errors[0]/3.5


def test_radial_seeds_preserve_mass_and_real_gram_normalization():
    r=np.arange(1,2048)*.03125;dr=.03125;b=.001;u=r*np.exp(-r*r/9)
    q,norm,condition=seed_basis(u,r,dr,1.,b,0)
    np.testing.assert_allclose(dr*np.real(q.conj().T@q),np.eye(4),atol=1e-11)
    derivative=2*dr*np.real((r*np.conj(b+u/r))@q)
    np.testing.assert_allclose(derivative,0,atol=1e-11)
    assert norm>0 and condition<1e10


def test_density_gram_maximizes_every_real_combination_in_the_declared_span():
    rng=np.random.default_rng(4);r=np.arange(1,301)*.05;dr=.05;b=.001
    u=r*np.exp(-r*r/8);q,norm,_=seed_basis(u,r,dr,1.,b,2)
    q=q@np.diag([1.,2.,3.,4.]);S=dr*np.imag(q.conj().T@q)
    row=diagnostics(q,u,r,dr,1.,b,2,norm,S)
    psi=b+u/r;D=r[:,None]*2*np.real(np.conj(psi[:,None])*q/r[:,None])
    contrast=abs(u/r)**2+2*b*np.real(u/r);den=np.linalg.norm(r*contrast)
    for _ in range(30):
        c=rng.normal(size=4);c/=np.linalg.norm(c)
        assert norm*np.linalg.norm(D@c)/den<=row['core_density_gain']+1e-12


def test_quadrupole_moment_normalization_matches_full_angular_integration():
    r=np.arange(1,801)*.025;dr=.025;b=.001;u=r*np.exp(-r*r/8)
    q=(r**3*np.exp(-r*r/8))[:,None].astype(complex)
    q/=math.sqrt(float(dr*np.sum(abs(q)**2)))
    row=diagnostics(q,u,r,dr,1.,b,2,1.,np.zeros((1,1)))
    mu,weights=np.polynomial.legendre.leggauss(12);P2=(3*mu*mu-1)/2;Y=math.sqrt(5)*P2
    psi=b+u/r;eps=1e-5;moments=[]
    for sign in (1,-1):
        density=abs(psi[:,None]+sign*eps*(q/r[:,None])*Y[None,:])**2-b*b
        moments.append(float(np.sum(r**4*(density@(weights*P2)))/np.sum(r**4*(density@weights))))
    derivative=(moments[0]-moments[1])/(2*eps)
    assert derivative==pytest.approx(row['quadrupole_moment_gain'],rel=2e-9)
