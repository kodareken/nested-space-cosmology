import math
import numpy as np
import pytest

from recursive_horizons.nsc_atlas_core_optics import homogeneous_screen
from recursive_horizons.nsc_atlas_core_images import invert_patch, paired_observables, patches, AGES


def test_fixed_worldline_inverse_recovers_a_flat_chord_and_fermat_gradient():
    centre=.02; rows=[]
    for psi in centre+np.linspace(-.00008,.00008,5):
        D=1.2*math.cos(psi)
        rows.append(dict(branch=1, sky_angle=psi, emission_proper_age=2.8,
            landing_angle=math.pi-2*psi, arrival_proper_age=2.8+D,
            received_frequency=1., received_radial_cosine=math.cos(psi),
            jacobi_perpendicular=D, jacobi_parallel=D, signed_beam_area=D*D,
            frequency_transport_gap=0.))
    r=invert_patch(rows,offset=.04)
    assert r['launch_angle']==pytest.approx(centre,abs=2e-15)
    assert r['received_sky_angle']==pytest.approx(centre,abs=2e-15)
    assert abs(r['endpoint_variation_gap'])<2e-10
    assert abs(r['parallel_jacobi_gap'])<2e-10
    assert abs(r['independent_root_gap'])<2e-15
    assert r['fluence_weight']==pytest.approx(1/(1.2*math.cos(centre))**2,rel=1e-13)


def test_independent_evolving_flrw_inverse_has_correct_clock_and_energy_factors():
    rows=homogeneous_screen(epochs=(2.8,),angles=.02+np.linspace(-.00008,.00008,5))
    for r in rows:r.update(branch=1,frequency_transport_gap=0.)
    offset=math.pi-rows[2]['landing_angle']
    r=invert_patch(rows,offset=offset)
    assert r['launch_angle']==pytest.approx(.02,abs=1e-13)
    assert r['received_sky_angle']==pytest.approx(math.acos(rows[2]['received_radial_cosine']),abs=2e-13)
    assert abs(r['endpoint_variation_gap'])<3e-8
    assert abs(r['parallel_jacobi_gap'])<3e-7
    assert r['luminosity_flux_weight']/r['fluence_weight']==pytest.approx(r['received_frequency'])


def test_inverse_rejects_a_missing_image_and_nonmonotone_patch():
    rows=[dict(branch=1,sky_angle=x,landing_angle=y, jacobi_parallel=1.,received_radial_cosine=1.)
        for x,y in zip((.01,.02,.03), (3.2,3.1,3.15))]
    with pytest.raises(ValueError,match='monotone'):invert_patch(rows)
    rows[-1]['landing_angle']=3.
    with pytest.raises(ValueError):invert_patch(rows,offset=.5)


def test_opposite_paths_reach_one_receiver_but_fluence_is_not_flux():
    offset=.0002
    assert math.sin(-(math.pi+offset))==pytest.approx(math.sin(math.pi-offset),abs=1e-15)
    assert math.cos(-(math.pi+offset))==pytest.approx(math.cos(math.pi-offset),abs=1e-15)
    a=dict(received_sky_angle=.02,arrival_proper_age=4.,received_frequency=2.,
        fluence_weight=2/3,luminosity_flux_weight=4/3)
    b=dict(received_sky_angle=-.01,arrival_proper_age=4.1,received_frequency=3.,
        fluence_weight=3/4,luminosity_flux_weight=9/4)
    r=paired_observables(a,b)
    assert r['image_separation']==pytest.approx(.03)
    assert r['delay_negative_minus_positive']==pytest.approx(.1)
    assert r['predicted_delay_drift']==pytest.approx(-1/6)
    assert r['luminosity_flux_ratio_negative_over_positive']==pytest.approx(
        r['fluence_ratio_negative_over_positive']*r['frequency_ratio_negative_over_positive'])


def test_each_epoch_has_both_independent_angular_patches():
    angles,epochs,branches=patches()
    assert len(angles)==50
    for age in AGES:
        for branch in (1,-1):
            x=angles[(epochs==age)&(branches==branch)]
            assert len(x)==5 and np.all(np.diff(x)>0)
