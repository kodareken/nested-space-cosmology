import math
from fractions import Fraction as F
import numpy as np
import pytest
import sympy as s
from scipy.integrate import quad
from recursive_horizons.nsc_atlas_curvature_mechanisms import FlatVolumeControl,curvature_audit,exact_energy_obstruction


def test_curvature_similarity_and_dust_growth():
    a,H,k,r=s.symbols('a H k r',positive=True)
    assert s.simplify(-k/((r*a)**2*(H/r)**2)+k/(a*a*H*H))==0
    # H^2 proportional to a^-3 in a dust-dominated comparison.
    q=s.symbols('q',positive=True)
    assert s.simplify(s.diff(s.log(k/a**2/(q/a**3)),a)*a)==1


def test_energy_obstruction_and_late_flattening_remove_matter():
    r=curvature_audit();e=exact_energy_obstruction()
    assert F(e['potential_at_bound_lower'])>F(21,32)
    assert e['curvature_over_matter_lower']>.11
    assert e['matter_upper_at_curvature_0_002']<.019
    assert r['threshold_crossing']['N']==pytest.approx(2.52065449035,abs=1e-8)
    assert r['threshold_crossing']['Omega_m']<.0012
    assert abs(r['declared_N6_endpoint']['Omega_k'])<2e-6
    assert not r['inherited_observed_small_curvature_derived']


def test_imposed_flat_background_obeys_its_constraint_and_direct_null_integral():
    bg=FlatVolumeControl(2.33);No=.43;z=np.array([.1,1.,2.33])
    DL,DM,DH=bg.distances(No,z,z);Ho=bg.state([No])[5][0]
    for i,zz in enumerate(z):
        I=quad(lambda N:math.exp(-N)/bg.state([N])[5][0],No-math.log1p(zz),No,epsabs=1e-11)[0]
        assert DM[i]==pytest.approx(Ho*math.exp(No)*I,rel=3e-9)
    f=bg.fractions(No)
    assert f['Omega_k']==0
    assert f['Omega_dark_energy']+f['Omega_m']==pytest.approx(1.,abs=1e-14)
