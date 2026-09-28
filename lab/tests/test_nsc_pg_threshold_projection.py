"""New subgap integration checks; synthetic contours are numerical controls."""
import json
from pathlib import Path

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss

from recursive_horizons.nsc_pg_packet_modes import project_computed_modes,project_horizon_bases
from recursive_horizons.nsc_pg_threshold_projection import closed_fiber_split,contour_integral

ROOT=Path(__file__).resolve().parents[1]

def test_physical_closed_fiber_split_reconstructs_saved_mode_projection():
    record=json.loads((ROOT/'results/development/nsc-pg-massive-mode-resolution.json').read_text())
    with np.load(ROOT/record['payload']['path'],allow_pickle=False) as a:
        i=1172  # declared group13, one existing subgap node, no new state
        label=a['label'][i];E=label[5];R=a['reflection'][i].item()
        assert label[0]==13 and E<np.pi/2 and a['transmission'][i].item()==0
        ext=np.column_stack((a['exterior'][i,1,:,0]-R*a['exterior_ingoing'][i,1],a['exterior_ingoing'][i,1]))
        basis,_=project_horizon_bases(E,np.pi/2,0.,a['interior_fundamental'][i,1],ext)
        modes,_=project_computed_modes(E,np.pi/2,0.,a['interior'][i,1],a['exterior'][i,1])
        G0,G1,K0,K1=closed_fiber_split(E,.23832579963401956,basis,basis)
        C=a['source_covariance'][i];P=a['source_projector'][i]
    assert np.linalg.norm(G0+R*G1+R.conjugate()*G1.conj().T-modes@P@modes.conj().T)<3e-11
    assert np.linalg.norm(K0+R*K1+R.conjugate()*K1.conj().T-modes@(C-.5*P)@modes.conj().T)<3e-11

def test_complex_contour_matches_real_integral_with_analytic_dual():
    # A finite analytic matrix and phase reflector test the new integrator.
    # They are not an NSC mode map, state completion, or selected boundary.
    base=np.arange(1,17,dtype=float).reshape(8,2)/16
    slopes=np.linspace(-.8,1.2,16).reshape(8,2)
    basis=lambda z:base*np.exp(1j*z*slopes)
    reflection=lambda z:np.exp(7j*z)
    kappa=.8;left=.2;right=.7
    x,w=leggauss(80);direct=np.zeros((2,8,8),complex)
    for t,weight in zip(x,w):
        E=left+(t+1)*(right-left)/2
        G0,G1,K0,K1=closed_fiber_split(E,kappa,basis(E),basis(E))
        R=reflection(E)
        direct+=weight*(right-left)/2*np.array([G0+R*G1+R.conjugate()*G1.conj().T,K0+R*K1+R.conjugate()*K1.conj().T])
    for height in (kappa/4,kappa/6):
        result=contour_integral(left,right,kappa,basis,reflection,points=24,height=height)
        assert np.linalg.norm(result.gram_positive_energy-direct[0])<3e-11
        assert np.linalg.norm(result.centered_positive_energy-direct[1])<3e-11

def test_contour_rejects_crossing_a_horizon_occupation_pole():
    with pytest.raises(ValueError,match='pole-free'):
        contour_integral(.2,.7,.8,lambda z:np.zeros((8,2)),lambda z:0.,height=.4)
