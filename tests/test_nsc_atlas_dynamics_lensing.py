import math
from pathlib import Path
from types import SimpleNamespace
import json

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_dynamics_lensing import RadialFall,impulse,aberrate,predict_from_motion


class SimpleMetric:
    center=0.;start=-.02;end=.02;length=1.;c=0.
    bg=SimpleNamespace(H=lambda u:.4)
    def __init__(self,force=0.):self.force=force
    def evaluate(self,N,r,homogeneous=False,**kwargs):
        return np.zeros(3) if homogeneous else np.array([self.force*(r-.1),self.force,0.])


def test_free_particle_momentum_redshifts_without_gravity():
    s=SimpleMetric();fall=RadialFall(s,.1,velocity=.05)
    initial=fall.at(s.start)
    for N in np.linspace(s.start,s.end,11):
        value=fall.at(N);v=value['velocity']
        invariant=math.exp(N)*v/math.sqrt(1-v*v)
        expected=math.exp(s.start)*initial['velocity']/math.sqrt(1-initial['velocity']**2)
        assert invariant==pytest.approx(expected,rel=1e-11)
        assert value['momentum']==pytest.approx(initial['momentum'],rel=1e-12)


def test_released_geodesic_impulse_recovers_lapse_gradient():
    s=SimpleMetric(.01);errors=[]
    for h in (.001,.0005):
        row=impulse(s,0.,.1,h)
        errors.append(abs(row['inferred_relative_acceleration']-.01))
    assert errors[1]<errors[0]/3.8
    assert errors[1]<1e-7


def test_aberration_is_invertible():
    for angle in (.06,.7,2.2):
        for v in (-.2,.1):
            assert aberrate(aberrate(angle,v),-v)==pytest.approx(angle,abs=3e-15)


def test_data_only_predictor_needs_no_potential_or_mass_model():
    data=dict(receiver_azimuth=math.pi-.12,source={'radius':.25},
        receiver={'radius':.25,'velocity':0.},curvature=0.,background_angle=.06,
        probes=[dict(quadrature_weight=.1,born_weight=.5,force_conversion_factor=2.,
                     inferred_relative_acceleration=1e-6,release_acceleration_envelope=1e-8)])
    result=predict_from_motion(data)
    assert result['predicted_image_shift']==pytest.approx(1e-7,abs=3e-16)
    assert result['release_state_envelope']==pytest.approx(1e-9)


def test_frozen_forecast_projection_and_no_observational_claim():
    root=Path(__file__).resolve().parents[1]
    path=root/'results/summaries/nsc-p004-dynamics-lensing-forecast-v1.json'
    if not path.exists():pytest.skip('forecast not yet frozen')
    r=json.loads(path.read_text(encoding='utf-8'))
    assert r['withheld_ray_evaluated'] is False and r['observational_validation'] is False
    assert predict_from_motion(r['calibration'])['predicted_image_shift']==pytest.approx(r['prediction']['value_radians'],abs=2e-16)
