"""Independent finite-state variational, energy-vertex and massless controls."""
import math

import numpy as np
import pytest
import sympy as sp
from scipy.special import expit

from recursive_horizons.nsc_atlas_dirac_preparation import (
    averaged_clock_energy,clock_energy_matrix,hermite_null_modes,minimum_clock_covariance,
    propagate_null_modes,ClockEnergyBasis)
from recursive_horizons.nsc_atlas_dirac_stress import stress_matrices,canonical_metric_vertices
from recursive_horizons.nsc_atlas_dirac_scattering import inherit_boundary_state
from recursive_horizons.nsc_boundary_state import GaussianBoundaryState


def test_full_clock_energy_is_the_independent_action_vertex_with_shift():
    rng=np.random.default_rng(317)
    xi=rng.normal(size=(2,2,3))+1j*rng.normal(size=(2,2,3))
    dx=rng.normal(size=xi.shape)+1j*rng.normal(size=xi.shape)
    fields=dict(lapse=np.array([.8,1.3]),radial_length=np.array([.9,1.2]),
        radius=np.array([.41,.53]),shift=np.array([.12,-.08]))
    stress=stress_matrices(xi,dx,radial_length=fields['radial_length'],
        radius=fields['radius'],angular_k=-2)
    expected=canonical_metric_vertices(xi,dx,angular_k=-2,**fields)['hamiltonian_density']
    assert clock_energy_matrix(fields,stress)==pytest.approx(expected,rel=2e-14,abs=2e-13)
    weights=np.array([.3,.7])
    averaged=averaged_clock_energy(expected,weights)
    assert averaged==pytest.approx(.3*expected[0]+.7*expected[1],abs=2e-13)
    with pytest.raises(ValueError,match='timelike'):
        clock_energy_matrix({**fields,'shift':np.array([2.,2.])},stress)


def test_energy_minimizer_variational_CAR_basis_units_and_inheritance():
    rng=np.random.default_rng(713)
    u,_=np.linalg.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
    A=(u*np.array([-3.,-1.,2.,7.]))@u.conj().T
    result=minimum_clock_covariance(A)
    C=result['covariance']
    assert result['status']=='finite_clock_energy_minimizer'
    assert C==pytest.approx(u[:,:2]@u[:,:2].conj().T,abs=2e-14)
    assert C@C==pytest.approx(C,abs=2e-14)
    assert result['minimum_energy']==pytest.approx(-4.,abs=2e-14)
    # The variational inequality is checked against independently supplied
    # admissible mixed covariances, not merely against the eigensolver output.
    for unused in range(30):
        v,_=np.linalg.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
        trial=(v*rng.uniform(0,1,4))@v.conj().T
        assert np.trace(trial@A).real>=-4.-3e-13
    assert minimum_clock_covariance(A*1e-25)['covariance']==pytest.approx(C,abs=3e-14)
    basis,_=np.linalg.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))
    transformed=minimum_clock_covariance(basis@A@basis.conj().T)['covariance']
    assert transformed==pytest.approx(basis@C@basis.conj().T,abs=3e-14)
    parent=GaussianBoundaryState(C,np.array([True,True,False,False]))
    inherited=inherit_boundary_state(parent,basis,
        outgoing_parent_mask=np.array([False,True,False,True]))
    assert inherited.covariance==pytest.approx(basis@C@basis.conj().T,abs=3e-14)
    assert np.linalg.norm(C[:2,2:])>.1  # Full correlations are retained.


def test_uncertain_spectral_gap_does_not_select_a_state():
    result=minimum_clock_covariance(np.diag([-1.,1e-8]),operator_uncertainty=2e-8)
    assert result['status']=='unresolved_spectral_gap'
    assert result['covariance'] is None
    assert minimum_clock_covariance(np.zeros((2,2)))['covariance'] is None
    with pytest.raises(ValueError,match='nonnegative'):
        averaged_clock_energy(np.repeat(np.eye(2)[None],2,axis=0),np.array([1.,-1.]))


def test_null_Hermite_modes_flux_orthogonality_derivatives_and_length_scaling():
    y,weights=np.polynomial.hermite.hermgauss(48)
    S=2.4
    s=S*expit(y)
    J=S/(s*(S-s))
    f,df=hermite_null_modes(s,S,6)
    polynomial=f/np.sqrt(J[:,None])*np.exp(y[:,None]**2/2)
    gram=np.einsum('p,pi,pj->ij',weights,polynomial,polynomial)
    assert gram==pytest.approx(np.eye(6),abs=2e-13)
    points=np.array([.23,.41,.67])*S
    values,derivative=hermite_null_modes(points,S,4)
    h=2e-5*S
    samples=[hermite_null_modes(points+k*h,S,4)[0] for k in (-2,-1,1,2)]
    independent=(samples[0]-8*samples[1]+8*samples[2]-samples[3])/(12*h)
    assert derivative==pytest.approx(independent,abs=2e-9)
    unit,unit_derivative=hermite_null_modes(points/S,1.,4)
    assert values*math.sqrt(S)==pytest.approx(unit,abs=2e-14)
    assert derivative*S**1.5==pytest.approx(unit_derivative,abs=2e-13)
    end,end_derivative=hermite_null_modes(np.array([-.1,0.,S,S+.1]),S,4)
    assert not np.any(end) and not np.any(end_derivative)


def test_massless_signature_block_keeps_half_frequency_mixing_when_small():
    # Independent exact projector algebra, without applying a positive-mass
    # Hadamard theorem to the actual massless field.
    a,b=sp.symbols('a b',real=True)
    norm=sp.sqrt(a*a+b*b)
    B=sp.Matrix([[0,a+sp.I*b],[a-sp.I*b,0]])
    P=(sp.eye(2)-B/norm)/2
    ground=sp.diag(0,1)
    assert sp.simplify(P*P-P)==sp.zeros(2)
    assert P[0,0]==P[1,1]==sp.Rational(1,2)
    assert sp.simplify((P-ground)**2-sp.eye(2)/2)==sp.zeros(2)
    # Tiny Fourier coefficients do not make the massless spectral projector
    # approach the static energy projector; a zero block is undetermined.
    for q in (1.,1e-12,1e-100):
        block=np.array([[0.,q],[q,0.]])
        e,u=np.linalg.eigh(block)
        state=u[:,e<0]@u[:,e<0].T
        assert np.diag(state)==pytest.approx([.5,.5],abs=3e-16)


def test_parallel_free_modes_keep_analytic_jets_joint_flux_and_superposition():
    n=12;h=1/n;s=(np.arange(n)+.5)*h
    f=np.c_[1+s+1j*s*s,np.exp(-2j*s)]
    df=np.c_[1+2j*s,-2j*np.exp(-2j*s)]
    plus=np.c_[f,3*f[:,0]-2j*f[:,1]]
    minus=np.c_[f[:,::-1],3*f[:,1]-2j*f[:,0]]
    du=np.c_[df,3*df[:,0]-2j*df[:,1]]
    dv=np.c_[df[:,::-1],3*df[:,1]-2j*df[:,0]]
    zero=np.zeros((n,n))
    grids=propagate_null_modes(zero,h,h,potential_u=zero,potential_v=zero,
        incoming_plus=plus,incoming_minus=minus,plus_u=du,minus_v=dv,workers=3)
    for k,grid in enumerate(grids):
        assert grid.jets[0][...,0]==pytest.approx(np.broadcast_to(plus[:,k,None],(n,n)),abs=3e-14)
        assert grid.jets[0][...,1]==pytest.approx(np.broadcast_to(minus[None,:,k],(n,n)),abs=3e-14)
        assert grid.jets[1][...,0]==pytest.approx(np.broadcast_to(du[:,k,None],(n,n)),abs=3e-14)
        assert grid.jets[2][...,1]==pytest.approx(np.broadcast_to(dv[None,:,k],(n,n)),abs=3e-14)
        assert not np.any(grid.jets[1][...,1]) and not np.any(grid.jets[2][...,0])
    outgoing=np.stack([grid.outgoing for grid in grids],axis=-1)
    incoming=np.sqrt(h)*np.r_[plus,minus]
    assert outgoing.conj().T@outgoing==pytest.approx(incoming.conj().T@incoming,abs=4e-14)
    for jet in range(4):
        assert grids[2].jets[jet]==pytest.approx(3*grids[0].jets[jet]-2j*grids[1].jets[jet],abs=4e-14)


class StaticCylinder:
    """Manufactured R^(1,1) x S^2 metric; not a nested solution."""
    duration=1.

    def radius_at_clock(self,clock):
        return np.asarray(clock)

    def from_wall_radii(self,U,V,*,region,angular_k):
        shape=np.broadcast_shapes(np.shape(U),np.shape(V))
        one=np.ones(shape);zero=np.zeros(shape)
        return dict(lapse=one,radial_length=one,radius=one,shift=zero,
            potential=angular_k*one,potential_t=zero,potential_x=zero)


def test_physical_energy_selection_has_signed_sector_parity_and_full_flux():
    states=[]
    parity=np.diag([1.,1.,-1.,-1.])
    for k in (1,-1):
        basis=ClockEnergyBasis(StaticCylinder(),16,2,preparation_duration=1.,angular_k=k,workers=4)
        assert basis.flux_gram_defect<2e-13
        state=basis.select(2)
        assert state['status']=='finite_clock_energy_minimizer'
        C=state['covariance']
        assert C@C==pytest.approx(C,abs=2e-14)
        assert np.linalg.norm(C[:2,2:])>.01
        states.append((state,basis.pressure_difference(.7,.4,state)))
    plus,minus=states
    assert minus[0]['generator']==pytest.approx(parity@plus[0]['generator']@parity,abs=3e-13)
    assert minus[0]['covariance']==pytest.approx(parity@plus[0]['covariance']@parity,abs=3e-13)
    assert minus[1]==pytest.approx(plus[1],abs=3e-13)
