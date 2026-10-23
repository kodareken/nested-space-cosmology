import math
import numpy as np
import pytest
from scipy.integrate import solve_ivp

from recursive_horizons.nsc_atlas_core_optics import ray_rates,angular_fan,homogeneous_screen


def static_screen(mass,psi):
    r=np.linspace(.004,.8,1601);z=np.zeros_like(r);B=1-2*mass/r;A=np.sqrt(B)
    g=(z,z,z,np.full_like(r,mass),B,2*mass/r**2,A,mass/(r*r*A),r,np.ones_like(r),z,z)
    bg=np.array([0.,0.,1e30,0.]);screen=.6;L=screen*math.sin(psi)
    y=np.array([screen,0.,-math.cos(psi)/math.sqrt(1-2*mass/screen),0.,1.,0.,1.,0.,0.])
    def rhs(t,y):return ray_rates(y,L,0.,g,g,g,bg,bg,bg,2000.)
    def endpoint(t,y):return y[0]-screen
    endpoint.direction=1;endpoint.terminal=True
    out=solve_ivp(rhs,(0,2),y,method='DOP853',rtol=2e-12,atol=2e-13,events=endpoint,max_step=.003)
    assert out.success and len(out.t_events[0])==1
    return out.t_events[0][0],out.y_events[0][0]


def test_flat_screen_recovers_straight_rays_and_affine_beam_area():
    psi=.08;t,y=static_screen(0.,psi)
    length=1.2*math.cos(psi)
    assert t==pytest.approx(length,abs=2e-10)
    assert y[1]==pytest.approx(math.pi-2*psi,abs=2e-10)
    assert y[3]==pytest.approx(length,abs=2e-10)
    assert y[5]==pytest.approx(length,abs=2e-10)
    assert y[8]==pytest.approx(length,abs=2e-10)


def test_static_vacuum_weyl_focusing_and_killing_energy():
    psi=.1;mass=.001;t,y=static_screen(mass,psi)
    beta=math.sqrt(1-2*mass/.6);L=.6*math.sin(psi)
    E=math.sqrt(beta*beta*y[2]**2+L*L/.6**2)
    assert E==pytest.approx(1.,abs=2e-9)
    assert math.exp(y[7])==pytest.approx(E,abs=2e-9)
    assert y[3]==pytest.approx(.6*math.sin(y[1])/math.sin(psi),abs=2e-8)
    eps=1e-4
    _,lo=static_screen(mass,psi-eps);_,hi=static_screen(mass,psi+eps)
    dphi=(hi[1]-lo[1])/(2*eps)
    assert y[5]==pytest.approx(-.6*(beta*y[2]/E)*dphi,abs=3e-6)
    assert abs(y[5]-y[3])>.01


def test_homogeneous_great_circle_and_isotropic_jacobi_fields():
    r=homogeneous_screen(epochs=(2.8,),angles=(.02,.04))
    for row in r:
        assert row['jacobi_perpendicular']==row['jacobi_parallel']
        symmetry=.6*math.sin(row['landing_angle'])/math.sin(row['sky_angle'])
        assert row['jacobi_perpendicular']==pytest.approx(symmetry,abs=2e-9)
        assert row['arrival_proper_age']>row['emission_proper_age']
        assert row['received_frequency']*math.sqrt(1-row['received_radial_cosine']**2)==pytest.approx(math.sin(row['sky_angle']),abs=1e-12)


def test_angular_checks_use_one_fixed_physical_family():
    angles,epochs=angular_fan()
    assert len(angles)==len(epochs)
    for age in (2.8,3.4):
        values=angles[epochs==age]
        for x in (.02,.04,.08):
            for delta in (0.,-.0002,.0002,-.0001,.0001):
                assert any(abs(v-x-delta)<1e-12 for v in values)


def test_curvature_and_boost_identities_are_exact():
    from scripts.check_nsc_atlas_core_optics import identities
    assert set(identities().values())=={'0'}


def test_changing_the_upstream_potential_or_wall_geometry_invalidates_cache(tmp_path):
    from scripts.check_nsc_atlas_core_optics import cache_identity,SOURCE_INPUTS
    from recursive_horizons.nsc_atlas_core_optics import OpticalRun
    # The two concrete upstream changes must alter the identity, even when
    # the optics and field-stepper files themselves stay unchanged.
    upstream=('src/recursive_horizons/nsc_atlas_flux_conversion.py',
              'src/recursive_horizons/nsc_atlas_flux_transition.py')
    for name in set(SOURCE_INPUTS)|set(upstream):
        p=tmp_path/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text('original source\n')
    before=cache_identity(OpticalRun(),tmp_path)
    for name in upstream:
        p=tmp_path/name;p.write_text('changed physical source\n')
        assert cache_identity(OpticalRun(),tmp_path)!=before
        p.write_text('original source\n')
    assert cache_identity(OpticalRun(),tmp_path)==before
