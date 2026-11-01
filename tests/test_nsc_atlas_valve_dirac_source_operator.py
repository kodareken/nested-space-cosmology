"""Independent proper-polar and physical-limit controls of the full Hessian.

Only initial geometry and finite angular sectors are evaluated. No field
evolution, occupation, renormalized stress or source backreaction is tested.
"""
import math

import numpy as np
import pytest
import scipy.linalg
from scipy.special import eval_gegenbauer,gammaln,roots_jacobi
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_valve_neck import solve_necks
from recursive_horizons.nsc_atlas_parent_scale_selector import VolumeParent
from recursive_horizons.nsc_atlas_valve_dirac_modes import initial_parent_mode_pump
from recursive_horizons.nsc_atlas_valve_source_gate import parent_dirac_source_operator


@pytest.fixture(scope='module')
def initial():
    with threadpool_limits(limits=1):return solve_necks(40,(1.,))[0]


@pytest.fixture(scope='module')
def actual(initial):
    sstar=next(row['s'] for row in VolumeParent().vacuum_roots() if row['kind']=='minimum')
    mass=math.sqrt(45/11)*math.exp(-7*sstar)
    with threadpool_limits(limits=1):
        levels={L:parent_dirac_source_operator(initial,L,mass=mass) for L in (16,24)}
    return mass,levels


def _radial(l,j,chi):
    # Gegenbauer norm specialized independently to fixed (l,j).
    logh=(math.log(math.pi)-(2*j+1)*math.log(2)+gammaln(l+j+2)
          -gammaln(l-j+1)-math.log(l+1)-2*gammaln(j+1))
    return math.exp(-float(logh)/2)*np.sin(chi)**j*eval_gegenbauer(l-j,j+1,np.cos(chi))


def _polar_columns(initial,L,k,chi):
    """Direct half-density spin-frame reflection O on scalar radial columns."""
    F=np.zeros((len(chi),2*(L-k)+1));G=np.zeros_like(F);n0=L-k+1
    for l in range(k,L+1):F[:,l-k]=_radial(l,k,chi)
    for l in range(k+1,L+1):G[:,n0+l-k-1]=_radial(l,k+1,chi)
    psi=initial.values(np.cos(chi));q=psi**2/initial.parent.H
    halfcos=np.cos(chi/2);halfsin=np.sin(chi/2)
    # The required map has determinant -1, not the superficially similar
    # SO(2) matrix with a positive bottom-right halfcos entry.
    polar0=halfcos[:,None]*F-halfsin[:,None]*G
    polar1=-halfsin[:,None]*F-halfcos[:,None]*G
    scale=np.sqrt(q)*np.sin(chi)
    return np.stack((scale[:,None]*polar0,scale[:,None]*polar1),axis=1),F,G,q


@pytest.mark.parametrize('k',[0,1])
def test_direct_proper_polar_projection_recovers_metric_gram_and_round_dirac(initial,k):
    L=5;u,w=roots_jacobi(100,.5,.5);chi=np.arccos(u);dchi_weights=w/(1-u*u)
    g,F,G,q=_polar_columns(initial,L,k,chi)
    gram=np.einsum('nAi,n,nAj->ij',g,dchi_weights,g)
    reference=F.T@((w*q)[:,None]*F)+G.T@((w*q)[:,None]*G)
    np.testing.assert_allclose(gram,reference,atol=3.e-13,rtol=3.e-13)
    # Independent fourth-order derivative of the entire proper-polar
    # half-density, including q and frame rotation, not its production jets.
    step=1.e-4
    gp2=_polar_columns(initial,L,k,chi+2*step)[0]
    gp=_polar_columns(initial,L,k,chi+step)[0]
    gm=_polar_columns(initial,L,k,chi-step)[0]
    gm2=_polar_columns(initial,L,k,chi-2*step)[0]
    derivative=(-gp2+8*gp-8*gm+gm2)/(12*step)
    psi=initial.values(u);psiu=initial.values(u,1)
    q_chi=-2*psi*psiu*np.sin(chi)/initial.parent.H
    speed=1/q;speed_chi=-q_chi/q**2
    angular=-(k+1)/(q*np.sin(chi))
    Dg=np.empty_like(g)
    Dg[:,0]=-speed[:,None]*derivative[:,1]-.5*speed_chi[:,None]*g[:,1]+angular[:,None]*g[:,1]
    Dg[:,1]=speed[:,None]*derivative[:,0]+.5*speed_chi[:,None]*g[:,0]+angular[:,None]*g[:,0]
    projected=np.einsum('nAi,n,nAj->ij',g,dchi_weights,Dg)
    n0=L-k+1;expected=np.diag([k+1.5]*n0+[-k-.5]*(L-k))
    for l in range(k+1,L+1):
        i=l-k;j=n0+l-k-1
        expected[i,j]=expected[j,i]=math.sqrt((l-k)*(l+k+2))
    # Finite-difference derivative budget; this is not a physical uncertainty.
    np.testing.assert_allclose(projected,expected,atol=3.e-9,rtol=3.e-10)


@pytest.mark.parametrize('k',[0,1])
def test_round_de_sitter_full_hessian_and_pair_matrix_have_exact_background_limit(k):
    initial=solve_necks(12,(0.,))[0];H=initial.parent.H;mass=2.
    with threadpool_limits(limits=1):
        summary,matrices=parent_dirac_source_operator(initial,6,mass=mass,k=k)
    eigenvalues=matrices['spatial_eigenvalues']
    np.testing.assert_allclose(matrices['spatial_Dirac_hessian'],-H*H*np.diag(eigenvalues),
                               atol=2.e-11,rtol=2.e-11)
    np.testing.assert_allclose(matrices['mass_hessian'],0.,atol=2.e-14,rtol=0)
    expected=-H*H*mass*eigenvalues/np.hypot(mass,eigenvalues)
    np.testing.assert_allclose(matrices['positive_negative_pair_hessian'],np.diag(expected),
                               atol=2.e-11,rtol=2.e-11)
    assert summary['angular_multiplicity']==2*k+2


def test_actual_same_mode_matches_independent_isotropic_slice_hessian(initial,actual):
    mass,levels=actual;_,matrices=levels[24]
    i=matrices['ground_index'];lam=matrices['spatial_eigenvalues'][i]
    mode=dict(degree=24,k=0,sign=1,eigenvalue=lam,coefficients=matrices['spatial_eigenvectors'][:,i])
    reference=initial_parent_mode_pump(initial,mode,mass,orders=(260,120))
    np.testing.assert_allclose(matrices['spatial_Dirac_hessian'][i,i],reference['lambda_ddot'],atol=1.e-9,rtol=1.e-8)
    np.testing.assert_allclose(matrices['mass_hessian'][i,i],reference['projected_mass_ddot'],atol=2.e-11,rtol=2.e-10)
    assert matrices['positive_negative_pair_hessian'][i,i]==pytest.approx(reference['total_pair_Hessian'],abs=1.e-9,rel=1.e-8)


def test_actual_inter_mode_coefficients_and_column_norm_are_resolved_under_cutoff(actual):
    mass,levels=actual
    norms=[]
    for L in (16,24):
        _,matrices=levels[L];i=matrices['ground_index'];pair=matrices['positive_negative_pair_hessian']
        column=pair[:,i];total=float(np.vdot(column,column).real)
        other=total-abs(pair[i,i])**2
        assert other>abs(pair[i,i])**2>0
        norms.append(total)
    # Empirical spectral control of a reference-amplitude coefficient norm,
    # not a particle/energy-conversion fraction or certified continuum bound.
    assert norms[1]==pytest.approx(norms[0],rel=2.e-6,abs=2.e-9)


def test_rephasing_spatial_eigenvectors_cannot_change_column_norm(initial,actual,monkeypatch):
    mass,levels=actual;baseline,matrices=levels[16]
    original=scipy.linalg.eigh
    def rephased(*args,**kwargs):
        values,vectors=original(*args,**kwargs)
        signs=np.where(np.arange(len(values))%2,1.,-1.)
        return values,vectors*signs[None,:]
    monkeypatch.setattr(scipy.linalg,'eigh',rephased)
    with threadpool_limits(limits=1):
        changed,other=parent_dirac_source_operator(initial,16,mass=mass)
    np.testing.assert_allclose(abs(other['positive_negative_pair_hessian']),abs(matrices['positive_negative_pair_hessian']),
                               atol=2.e-12,rtol=5.e-11)
    assert changed['ground_column_squared_coefficient_norm']==pytest.approx(
        baseline['ground_column_squared_coefficient_norm'],rel=2.e-11,abs=1.e-12)


def test_raw_massless_limit_removes_same_mode_but_retains_inter_mode_geometry(initial):
    # A raw operator limit only: untwisted internal S12 supplies a nonzero
    # mass, so this does not introduce a physical massless source candidate.
    with threadpool_limits(limits=1):_,matrices=parent_dirac_source_operator(initial,16,mass=0.)
    eigenvalues=matrices['spatial_eigenvalues'];pair=matrices['positive_negative_pair_hessian']
    same_sign=np.sign(eigenvalues[:,None])==np.sign(eigenvalues[None,:])
    np.testing.assert_allclose(pair[same_sign],0.,atol=2.e-13,rtol=0)
    assert np.max(abs(pair[~same_sign]))>.01


@pytest.mark.parametrize('change',[
    dict(degree=True),dict(degree=1),dict(k=-1),dict(k=4),
    dict(mass=-1.),dict(mass=float('nan')),dict(quadrature_order=12),
])
def test_invalid_source_operator_inputs_are_rejected(initial,change):
    arguments=dict(degree=4,mass=2.,k=0)
    arguments.update(change)
    with pytest.raises(ValueError):parent_dirac_source_operator(initial,**arguments)
