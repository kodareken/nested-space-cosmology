"""Independent optical embedding, causal bounds and physical radar controls."""
import math

import numpy as np
import pytest
from scipy.optimize import brentq

from recursive_horizons.nsc_atlas_cap_wall_response import reference_coefficients
from recursive_horizons.nsc_atlas_dirac_evolution import radial_time_coefficients
from recursive_horizons.nsc_atlas_dirac_evolution import ALPHA
from recursive_horizons.nsc_atlas_dirac_stress import SIGMA1
from recursive_horizons.nsc_atlas_dirac_radar import (
    MASS, RMIN, RMAX, WallRadarChart, direct_wall_maps, tortoise,
)


@pytest.fixture(scope='module')
def chart():
    return WallRadarChart()


def test_wall_optical_maps_against_independent_existing_quadrature(chart):
    for r in (RMIN,.5,RMAX):
        maps = chart.maps(r)
        assert maps == pytest.approx(direct_wall_maps(r,64),abs=2e-12)
        assert direct_wall_maps(r,32) == pytest.approx(direct_wall_maps(r,64),abs=3e-14)
        assert maps[2]-maps[1] == pytest.approx(2*tortoise(r,'cap'),abs=2e-12)
        assert maps[4]-maps[3] == pytest.approx(2*tortoise(r,'exterior'),abs=2e-12)
        assert float(chart.radius_at_clock(maps[0])) == pytest.approx(r,abs=5e-14)
    assert 0 < chart.duration < chart.bounds['wall_duration_upper']


@pytest.mark.parametrize('ru,rv',[(.475,.525),(.525,.475)])
def test_static_metric_pullback_tortoise_inverse_and_independent_jets(chart,ru,rv):
    row = chart.from_wall_radii(ru,rv)
    cap = ru >= rv
    region, side = ('cap','minus') if cap else ('exterior','plus')
    gu,gv = reference_coefficients(ru),reference_coefficients(rv)
    p,q = 1/(gu['beta_'+side]+gu['v']),1/(gv['beta_'+side]-gv['v'])
    r,f = row['radius'],row['static_f']
    Tt,Tx = (p+q)/2,(q-p)/2
    rt,rx = row['radius_t'],row['radius_x']
    pulled = np.array([[-f*Tt*Tt+rt*rt/f,-f*Tt*Tx+rt*rx/f],
                       [-f*Tt*Tx+rt*rx/f,-f*Tx*Tx+rx*rx/f]])
    assert pulled == pytest.approx(np.diag([-row['conformal_factor'],row['conformal_factor']]),abs=3e-13)
    mu,mv = direct_wall_maps(ru),direct_wall_maps(rv)
    target = (mv[2 if cap else 4]-mu[1 if cap else 3])/2
    assert tortoise(r,region) == pytest.approx(target,abs=2e-12)
    if not cap:
        independent = brentq(lambda x:float(tortoise(x,'exterior')-target),RMIN,RMAX,xtol=1e-14)
        assert r == pytest.approx(independent,abs=2e-12)
    mass = r*(1-(rx*rx-rt*rt)/row['conformal_factor'])/2
    assert mass == pytest.approx(r**3/2 if cap else MASS,abs=2e-13)
    # Differentiate the full numerical optical/inverse pullback, not its jets.
    U,V = float(row['clock_u']),float(row['clock_v'])
    h = 2e-5
    for coordinate in (0,1):
        samples = []
        for scale in (-2,-1,1,2):
            point = [U,V];point[coordinate] += scale*h
            samples.append(chart.fields(*point))
        for name,dt,dx in (('radius','radius_t','radius_x'),('lapse','lapse_t','lapse_x')):
            finite = (samples[0][name]-8*samples[1][name]+8*samples[2][name]-samples[3][name])/(12*h)
            exact = (row[dt]+(-1 if coordinate==0 else 1)*row[dx])/2
            assert finite == pytest.approx(exact,abs=3e-8)
    def mixed(step):
        return sum(a*b*math.log(float(chart.fields(U+a*step,V+b*step)['conformal_factor']))
                   for a in (-1,1) for b in (-1,1))/(4*step*step)
    finite_uv = (4*mixed(.0005)-mixed(.001))/3
    assert finite_uv == pytest.approx(row['log_C_uv'],abs=2e-5)
    # Two-metric conformal curvature and full four-dimensional trace.
    R2 = 4*row['log_C_uv']/row['conformal_factor']
    box_r = -4*row['radius_uv']/row['conformal_factor']
    R4 = R2+2*(1-f)/r**2-4*box_r/r
    assert R2 == pytest.approx(2. if cap else 4*MASS/r**3,abs=2e-13)
    assert R4 == pytest.approx(12. if cap else 0.,abs=2e-12)


def test_uniform_causal_bounds_and_constant_speed_wall_transmission(chart):
    spin = (np.eye(2)+1j*SIGMA1)/math.sqrt(2)
    assert spin.conj().T @ ALPHA @ spin == pytest.approx(np.diag([1.,-1.]),abs=3e-16)
    assert spin.conj().T @ SIGMA1 @ spin == pytest.approx(SIGMA1,abs=3e-16)
    rlabels = np.array([RMIN,.475,.5,.525,RMAX])
    row = chart.from_wall_radii(rlabels,rlabels[::-1])
    assert np.all(row['lapse'] >= chart.bounds['lapse_lower'])
    assert np.all(row['radius'] >= min(chart.bounds['cap_radius_lower'],RMIN))
    # Includes the signed beta_plus zero; no beta or static-f division there.
    for r in (.5,.16**(1/3)):
        cap = chart.from_wall_radii(r,r,region='cap')
        ext = chart.from_wall_radii(r,r,region='exterior')
        wall = reference_coefficients(r)
        for name,sign in (('cap',cap),('exterior',ext)):
            coefficients = radial_time_coefficients(lapse=sign['lapse'],lapse_l=sign['lapse_x'],
                radial_length=sign['radial_length'],radial_length_l=sign['lapse_x'],radius=sign['radius'])
            assert coefficients['speed'] == pytest.approx(1.)
            assert coefficients['derivative_companion'] == pytest.approx(0.,abs=1e-15)
            assert sign['lapse'] == pytest.approx(1.,abs=2e-12)
            assert sign['radius_t'] == pytest.approx(wall['v'],abs=2e-12)
            beta = wall['beta_minus' if name=='cap' else 'beta_plus']
            assert sign['radius_x'] == pytest.approx(beta,abs=2e-12)
            assert sign['lapse_t'] == pytest.approx(0.,abs=3e-12)
        assert ext['radius_x']-cap['radius_x'] == pytest.approx(-2*r,abs=3e-12)
        assert ext['lapse_x']-cap['lapse_x'] == pytest.approx(-2.,abs=3e-12)
        assert ext['potential_x'] == pytest.approx(cap['potential_x'],abs=3e-11)
        assert ext['potential_t'] == pytest.approx(cap['potential_t'],abs=3e-11)
        # xi=a^(1/2)*r*psi has continuous derivative in this radar frame.
        # The physical spinor still has the original parallel 3*k_wall/2 jump.
        normalization_jump = (ext['lapse_x']/2+ext['radius_x']/r
                              -cap['lapse_x']/2-cap['radius_x']/r)
        assert -normalization_jump == pytest.approx(3.,abs=3e-11)


def test_radar_frequency_is_a_finite_clock_protocol(chart):
    receive = chart.from_wall_radii(.475,.525)
    emit = chart.from_wall_radii(.475,.475)
    assert receive['radar_position'] > 0
    assert receive['time']-emit['time'] == pytest.approx(receive['radar_position'])
    # Phase=-U has a null gradient; static-radar observers have u^T=1/N.
    wave_covector = np.array([-1.,1.])
    inverse = np.diag([-1.,1.])/receive['conformal_factor']
    assert wave_covector @ inverse @ wave_covector == pytest.approx(0.,abs=1e-15)
    frequency = -wave_covector @ np.array([1/receive['lapse'],0.])
    assert frequency == pytest.approx(emit['lapse']/receive['lapse'],abs=2e-12)
    assert receive['clock_v']-receive['clock_u'] == pytest.approx(2*receive['radar_position'])


def test_radar_input_domain_and_region_controls(chart):
    with pytest.raises(ValueError): chart.from_wall_radii(.4,.5)
    with pytest.raises(ValueError): chart.fields(-.001,.01)
    with pytest.raises(ValueError): chart.fields(0.,chart.duration+.001)
    with pytest.raises(ValueError): chart.from_wall_radii(.5,.5,angular_k=0)
    with pytest.raises(ValueError): chart.from_wall_radii(.5,.5,region='flat')
    with pytest.raises(ValueError): tortoise(.9,'exterior')
