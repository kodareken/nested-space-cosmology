import numpy as np
import pytest
import sympy as s
from recursive_horizons.nsc_atlas_observer_geometry import circular_orbit,calculation,log_scale,redshift
from recursive_horizons.nsc_atlas_golden_vortex import kasner,PHI


def test_single_mass_circular_geodesic_needs_no_external_mass():
    r,L=s.symbols('r L',positive=True)
    V=(1-2/r)*(1+L/r**2)
    required=r*r/(r-3)
    assert s.simplify(s.diff(V,r).subs(L,required))==0
    assert s.simplify(V.subs(L,required)-(r-2)**2/(r*(r-3)))==0
    assert s.simplify(s.diff(V,r,2).subs(L,required)-2*(r-6)/(r**3*(r-3)))==0
    q=circular_orbit()
    assert q['E_squared']==pytest.approx(32/35)
    assert q['L_squared_over_mass_squared']==pytest.approx(100/7)
    assert q['stable'] and q['radial_velocity']==0
    assert q['spacetime_Kretschmann_times_mass4']>0


def test_anisotropic_metric_satisfies_all_vacuum_einstein_components():
    u,d=s.symbols('u d',positive=True)
    C=s.coth(u);S=1/s.sinh(u);h=C+3*d*S
    # Sum d_i=0 and sum d_i^2=2/3 from the Kasner constraints.
    spatial=3*s.diff(h,u)+3*C*h
    assert s.simplify(s.trigsimp(spatial-3))==0
    time_component=-9*S*S+3*C*C+6*S*S
    assert s.simplify(s.trigsimp(time_component-3))==0


def test_directional_redshift_and_late_isotropization():
    p=kasner(PHI)
    for i in range(3):
        expected=np.exp(log_scale(2,p)[i]-log_scale(1,p)[i])-1
        assert redshift(1,2,np.eye(3)[i],p)==pytest.approx(expected)
    r=calculation()['anisotropic']['rows']
    assert r[-1]['directional_ratio_max_over_min']<r[0]['directional_ratio_max_over_min']
    assert r[-1]['late_shear_fraction']<1e-6


def test_fading_does_not_algebraically_change_phase_speed():
    r=calculation()['wave']['twofold_expansion']
    assert r['frequency_ratio']*r['wavelength_ratio']==r['local_speed_ratio']==1
