"""Continuum-parent kernel, restriction, plane-wave and interface controls."""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import k0

from recursive_horizons.nsc_atlas_dirac_preparation import (
    cylinder_incoming_modes,cylinder_covariance_compression,CylinderPreparation,propagate_null_modes,
    _hardy_compression,PreparedCylinderModes)
from recursive_horizons.nsc_atlas_dirac_evolution import folded_transparent_interface
from recursive_horizons.nsc_atlas_dirac_radar import WallRadarChart
from recursive_horizons.nsc_atlas_dirac_compact_source import CompactMassChart
from recursive_horizons.nsc_common_subtracted_ks_source import SpectrumDeclaration,reference_symbol_kernel
from recursive_horizons.nsc_atlas_dirac_mean import SPIN_ROTATION


def test_cylinder_cross_leg_kernel_has_the_independent_Fourier_normalization():
    mu=2.3;s=.41
    fourier,error=quad(lambda p:1/np.hypot(p,mu),0.,np.inf,
        weight='cos',wvar=s,epsabs=2e-12,limlst=100)
    assert error<3e-12
    assert -mu*fourier/(2*math.pi)==pytest.approx(-mu*k0(mu*s)/(2*math.pi),abs=2e-12)


def test_PV_normalization_and_removable_diagonal_match_exact_line_Hermite_projector():
    x,w=np.polynomial.hermite.hermgauss(48)
    h0=np.exp(-x*x/2)/math.pi**.25
    f=np.c_[h0,math.sqrt(2)*x*h0]
    df=np.c_[-x*h0,math.sqrt(2)*h0-x*f[:,1]]
    C,gram=_hardy_compression(x,w*np.exp(x*x),f,df)
    assert gram==pytest.approx(np.eye(2),abs=2e-13)
    assert C==pytest.approx(np.array([[.5,-1j/math.sqrt(2*math.pi)],
        [1j/math.sqrt(2*math.pi),.5]]),abs=2e-13)


def test_continuum_compressions_keep_CAR_coherence_subspaces_and_units():
    S=.1256873632901346;R=.45
    coarse=cylinder_covariance_compression(S,R,4,order=96)
    fine=cylinder_covariance_compression(S,R,4,order=128)
    assert fine['status']=='finite_root_quadrature'
    C=fine['covariance']
    assert C==pytest.approx(coarse['covariance'],abs=2e-9)
    assert fine['gram_discrepancy']<1e-10
    assert fine['eigenvalues'][0]>0 and fine['eigenvalues'][-1]<1
    assert fine['cross_leg_norm']>.01
    assert np.linalg.norm(C@C-C)>.05
    small=cylinder_covariance_compression(S,R,2,order=128)['covariance']
    indices=[0,1,4,5]
    assert small==pytest.approx(C[np.ix_(indices,indices)],abs=2e-13)
    parity=np.diag([1.]*4+[-1.]*4)
    negative=cylinder_covariance_compression(S,R,4,angular_k=-1,order=128)['covariance']
    assert negative==pytest.approx(parity@C@parity,abs=2e-13)
    scaled=cylinder_covariance_compression(S*1000,R*1000,4,order=128)['covariance']
    assert scaled==pytest.approx(C,abs=2e-13)
    unresolved=cylinder_covariance_compression(S,R,4,order=8)
    assert unresolved['status']=='quadrature_unresolved'
    assert unresolved['covariance'] is None


def test_root_plane_modes_are_the_same_exact_cylinder_solution_on_both_legs():
    S=.2;R=.45;momenta=np.array([-3.,0.,4.]);errors=[]
    for n in (8,16):
        h=S/n;s=(np.arange(n)+.5)*h
        data=cylinder_incoming_modes(s,momenta,radius=R,angular_k=-1)
        mu=data['angular_mass']
        grids=propagate_null_modes(np.full((n,n),mu),h,h,
            potential_u=np.zeros((n,n)),potential_v=np.zeros((n,n)),
            incoming_plus=data['plus'],incoming_minus=data['minus'],
            plus_u=data['plus_u'],minus_v=data['minus_v'],workers=3)
        U,V=.43*S,.61*S
        phase=np.exp(1j*(data['a']*U+data['b']*V))
        theta=np.arctan2(mu,momenta)
        exact=phase[None,:]*np.array([-np.sin(theta/2),np.cos(theta/2)])
        values=np.stack([g.evaluate(U,V) for g in grids],axis=-1)
        dx=np.stack([g.evaluate(U,V,derivative_v=1)-g.evaluate(U,V,derivative_u=1)
                     for g in grids],axis=-1)
        errors.append(max(np.max(np.abs(values-exact)),np.max(np.abs(dx-1j*momenta*exact))))
    assert errors[-1]<1e-3
    assert errors[0]/errors[1]>3


def test_folded_transmission_is_maximal_zero_flux_and_has_no_frozen_decaying_kernel():
    data=folded_transparent_interface()
    A,B,G,P,J=[data[k] for k in ('normal_matrix','angular_matrix','gamma',
                                 'allowed_projector','allowed_isometry')]
    assert G@G==pytest.approx(np.eye(4),abs=2e-15)
    assert G==pytest.approx(G.conj().T,abs=2e-15)
    assert G@A+A@G==pytest.approx(np.zeros((4,4)),abs=2e-15)
    assert G@B-B@G==pytest.approx(np.zeros((4,4)),abs=2e-15)
    assert P@A@P==pytest.approx(np.zeros((4,4)),abs=2e-15)
    assert np.linalg.matrix_rank(P)==2
    assert J.conj().T@J==pytest.approx(np.eye(2),abs=2e-15)
    for mu,zeta in ((2.3,0.),(-2.3,.7),(0.,1.)):
        decay=np.hypot(mu,zeta)
        M=(1j*decay*A+mu*B-1j*zeta*np.eye(4))@J
        # Independent lower bound: sum the two side norms. No decaying
        # frozen mode can lie in the transparent subspace.
        assert M.conj().T@M==pytest.approx(2*decay**2*np.eye(2),abs=3e-14)


def test_auxiliary_switch_preserves_wall_first_jets_and_actual_future_operator():
    chart=CompactMassChart(WallRadarChart(),1e-4);S=chart.duration
    preparation=CylinderPreparation(chart,duration=S,root_radius=.45)
    quiet=preparation.fields(np.array([0.,-.2*S,1.2*S]),np.array([.5*S,.5*S,.5*S]))
    assert quiet['potential']==pytest.approx(np.full(3,1/.45),abs=2e-15)
    assert not np.any(quiet['potential_t']) and not np.any(quiet['potential_x'])
    future=preparation.fields(.8*S,.7*S)
    actual=chart.fields(.8*S,.7*S)
    for name in ('potential','potential_t','potential_x'):
        assert future[name]==pytest.approx(actual[name],abs=3e-13)
    cap=CylinderPreparation(chart,duration=S,root_radius=.45,region='cap')
    ext=CylinderPreparation(chart,duration=S,root_radius=.45,region='exterior')
    for fraction in (.15,.35,.95):
        left,right=cap.fields(fraction*S,fraction*S),ext.fields(fraction*S,fraction*S)
        for name in ('potential','potential_t','potential_x'):
            assert left[name]==pytest.approx(right[name],abs=3e-12)
    U,V=.16*S,.18*S;h=2e-6*S
    center=preparation.fields(U,V)
    du=(preparation.fields(U+h,V)['potential']-preparation.fields(U-h,V)['potential'])/(2*h)
    dv=(preparation.fields(U,V+h)['potential']-preparation.fields(U,V-h)['potential'])/(2*h)
    assert center['potential_t']==pytest.approx(du+dv,rel=2e-7,abs=2e-7)
    assert center['potential_x']==pytest.approx(dv-du,rel=2e-7,abs=2e-7)


class ConformalCylinder:
    """Manufactured physical metric conformal to the exact root cylinder."""
    duration=.2

    def radius_at_clock(self,clock):
        return np.asarray(clock)

    def from_wall_radii(self,U,V,*,region,angular_k):
        U,V=np.broadcast_arrays(U,V);T=(U+V)/2;x=(V-U)/2
        N=1+.1*T+.2*x;r=.45*N
        return dict(lapse=N,radius=r,lapse_t=np.full(N.shape,.1),
            lapse_x=np.full(N.shape,.2),radius_t=np.full(N.shape,.045),
            radius_x=np.full(N.shape,.09))


def test_prepared_parent_kernel_matches_independent_local_Weyl_kernel_in_known_limit():
    preparation=CylinderPreparation(ConformalCylinder(),duration=.2,root_radius=.45)
    p=np.array([-3.,0.,4.]);weights=np.array([.2,.3,.4])
    declaration=SpectrumDeclaration('control_only','manufactured root Fourier panel')
    mu=1/.45;omega=np.hypot(p,mu)
    s1=np.array([[0.,1.],[1.,0.]])
    s3=np.diag([1.,-1.])
    projector=(np.eye(2)-(p[:,None,None]*s3+mu*s1)/omega[:,None,None])/2
    orders=np.zeros((5,3,2,2),dtype=complex)
    orders[0]=SPIN_ROTATION@projector@SPIN_ROTATION.conj().T
    eta=.008
    reference=reference_symbol_kernel(orders,p,weights,separation=eta,declaration=declaration)
    errors=[]
    for n in (8,16):
        modes=PreparedCylinderModes(preparation,n,p,workers=3)
        sample=modes.kernel(.11,.006,weights,separation=eta,declaration=declaration)
        errors.append(max(np.max(np.abs(sample.kernel-reference.kernel)),
            np.max(np.abs(sample.separation_derivative-reference.separation_derivative))))
        zero=modes.kernel(.11,.006,weights,separation=0.,declaration=declaration)
        assert zero.kernel==pytest.approx(zero.kernel.conj().T,abs=3e-14)
        assert np.linalg.eigvalsh(zero.kernel)[0]>=0
    assert errors[-1]<2e-4
    assert errors[0]/errors[1]>3
    with pytest.raises(ValueError,match='declaration'):
        modes.kernel(.11,.006,weights,separation=eta,declaration=None)
    with pytest.raises(ValueError,match='matched'):
        modes.kernel(.02,0.,weights,separation=0.,declaration=declaration)
