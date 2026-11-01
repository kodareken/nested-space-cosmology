"""Independent current, orientation and local-energy controls; no raw caches."""
import importlib.util
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.special import roots_jacobi

from recursive_horizons.nsc_atlas_valve_neck import NeckInitialData

spec=importlib.util.spec_from_file_location('directional_initial',
    Path(__file__).resolve().parents[1]/'scripts/check_nsc_atlas_valve_directional.py')
directional=importlib.util.module_from_spec(spec);spec.loader.exec_module(directional)


def preparation(H):
    # A nonconstant positive metric tests conformal weights independently
    # of the Hamiltonian solver; the momentum and energy identities hold off it.
    initial=SimpleNamespace(parent=SimpleNamespace(H=H,kappa=.9,vacuum=.22,mass_squared=7.))
    c=np.zeros(21);c[:3]=np.sqrt(2*np.pi**2)*np.array([.9,.03,-.01])
    A=.11;B=.37;C=initial.parent.kappa*A*B*math.sqrt(3/5)
    return initial,c,A,B,C


@pytest.mark.parametrize('H',[.2714056154003984,.73])
def test_physical_momentum_and_direction_reversal_at_different_scales(H):
    initial,c,A,B,C=preparation(H);u=np.linspace(-.98,.98,123)
    positive=directional.fields(initial,c,.8,u,A,B,C)
    negative=directional.fields(initial,c,-.8,u,A,B,C)
    assert max(abs(positive['momentum']))<2e-15
    np.testing.assert_array_equal(positive['rho'],negative['rho'])
    np.testing.assert_array_equal(positive['Pi'],-negative['Pi'])
    np.testing.assert_array_equal(positive['Kr'],-negative['Kr'])
    plus=directional.initial_cap_ledger(initial,c,.8,A,B,C,80)
    minus=directional.initial_cap_ledger(initial,c,-.8,A,B,C,80)
    assert plus['plus_cap_material_influx']>0
    assert plus['plus_cap_material_influx']==-minus['plus_cap_material_influx']


def test_initial_KG_volume_energy_derivative_matches_surface_flux_and_work():
    initial,c,A,B,C=preparation(.73);H=initial.parent.H
    x,w=roots_jacobi(100,.5,0.);u=(x+1)/2;w=w*np.sqrt(3+x)/4
    g=directional.fields(initial,c,1.,u,A,B,C)
    psi=g['psi'];du=NeckInitialData(initial.parent,40,1.,c,()).values(u,1)
    v=1-u*u;q=psi**2/H
    pi_u=2*B*np.sqrt(3/5)*u*psi**-6-6*g['f']*psi**-7*du
    lap_chi=A*(-3*u+2*v*du/psi)/q**2
    # Direct initial Klein-Gordon evolution cancels the potential derivatives.
    rho_dot=g['Pi']*lap_chi+v*A*pi_u/q**2+g['Kr']*g['gradient2']
    volume_rate=float(4*np.pi/H**3*np.sum(w*psi**6*rho_dot))
    ledger=directional.initial_cap_ledger(initial,c,1.,A,B,C,100)
    assert abs(volume_rate-ledger['plus_cap_energy_derivative_from_initial_ledger'])<2e-12
