import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from recursive_horizons.nsc_atlas_formed_observer import weighted_primitive,spectral_derivative,PotentialSlice,Observer,ray
from recursive_horizons.nsc_atlas_formation_window import FormationBackground


def test_weighted_primitive_is_exact_for_cubic_density():
    x=np.linspace(0,2,23);p=CubicSpline(x,1+2*x-3*x*x+x**3)
    for power in (1,2):
        primitive=weighted_primitive(p,power)
        points=np.linspace(0,2,31)
        exact=sum(c*points**(j+power+1)/(j+power+1) for j,c in enumerate([1,2,-3,1]))
        np.testing.assert_allclose(primitive(points)-primitive(0),exact,atol=5e-13)


def test_spectral_field_derivative_and_current_potential():
    L=8.;N=4096;r=np.arange(1,N)*L/N
    u=np.sin(3*np.pi*r/L)+.2j*np.sin(7*np.pi*r/L)
    expected=3*np.pi/L*np.cos(3*np.pi*np.arange(N+1)/N)+.2j*7*np.pi/L*np.cos(7*np.pi*np.arange(N+1)/N)
    np.testing.assert_allclose(spectral_derivative(u,L),expected,atol=2e-12)
    R=.2;eps=.002;theta=eps*np.maximum(0,1-r*r/(R*R))**5
    u=r*np.expm1(1j*theta);bg=FormationBackground.inherited();length=.0004815623798469198
    snapshot=PotentialSlice(r,u,1.2,length,bg)
    for z in (0,.05,.1,.15,.3):
        phi,grad,dN,dNr=snapshot.evaluate(z)
        exact=-snapshot.pref*snapshot.a*.5*eps*max(0,1-z*z/(R*R))**5
        assert abs(phi)<1e-20
        assert dN==pytest.approx(exact,rel=3e-6,abs=1e-19)


class PlummerMetric:
    length=1.;c=0.
    bg=SimpleNamespace(H=lambda u:0.)
    def __init__(self,mass):self.mass=mass
    def evaluate(self,N,r,homogeneous=False,frozen=False):
        if homogeneous:return np.zeros(3)
        core=.007;screen=.25
        return np.array([-self.mass/np.sqrt(r*r+core*core)+self.mass/np.sqrt(screen*screen+core*core),
                         self.mass*r/(r*r+core*core)**1.5,0.])


def test_finite_fixed_endpoints_against_independent_born_lens_and_delay():
    obs=Observer(ray_rtol=2e-13,ray_atol=2e-15);base_angle=.04;rs=obs.radius
    b=rs*math.sin(base_angle);z=rs*math.cos(base_angle);core=.007
    errors=[]
    for mass in (2e-8,1e-8):
        metric=PlummerMetric(mass)
        angle=brentq(lambda a:ray(metric,0.,a,obs)['arrival_angle']-obs.separation,.039,.041,xtol=2e-14)
        got=ray(metric,0.,angle,obs)
        born=2*mass*b*z/((b*b+core*core)*math.sqrt(rs*rs+core*core))
        delay=4*mass*(np.arcsinh(z/math.sqrt(b*b+core*core))-z/math.sqrt(rs*rs+core*core))
        errors.append(abs(got['observed_direction']-base_angle-born))
        assert got['proper_flight_time']-2*z==pytest.approx(delay,rel=3e-5,abs=1e-12)
        assert abs(got['expansion_removed_log_frequency'])<1e-20
    assert errors[1]<errors[0]/3.5


def test_curved_homogeneous_finite_screen_geometry():
    obs=Observer();metric=PlummerMetric(0.);metric.c=.02
    rs=obs.radius;theta=obs.separation
    angle=math.atan2(math.cos(theta/2),math.sqrt(1-metric.c*rs*rs)*math.sin(theta/2))
    got=ray(metric,0.,angle,obs,homogeneous=True)
    distance=2*np.arcsin(math.sqrt(metric.c)*rs*math.sin(theta/2))/math.sqrt(metric.c)
    assert got['arrival_angle']==pytest.approx(theta,abs=2e-11)
    assert got['null_conformal_length']==pytest.approx(distance,abs=2e-12)
    assert got['observed_direction']==pytest.approx(angle,abs=1e-12)


def test_saved_observer_source_and_error_scope():
    root=Path(__file__).resolve().parents[1];p=root/'results/summaries/nsc-atlas-formed-observer-v1.json'
    if not p.exists():pytest.skip('producer not yet run')
    r=json.loads(p.read_text(encoding='utf-8'))
    assert not r['numerical_validation_errors'] and not r['budget_exhausted']
    for f,h in {**r['source_sha256'],**r['evidence_sha256']}.items():
        assert hashlib.sha256((root/f).read_bytes()).hexdigest()==h,f
    assert r['runs']['fine']['source_final_field_sha256']==json.loads((root/'results/summaries/nsc-atlas-seed-collapse-v1.json').read_text())['runs']['fine']['final_field_sha256']
