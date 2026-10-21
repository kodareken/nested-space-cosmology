import math

import numpy as np
import pytest
from scipy.integrate import quad

from recursive_horizons.nsc_atlas_parent_perturbations import (
    ParentModes,ParentPacket,fixed_metric_control,jeans_comparison)


@pytest.fixture(scope='module')
def modes():
    m=ParentModes()
    return m.solve(),m.solve(clock='efolds',rtol=2e-11,atol=2e-13)


def test_constraint_and_full_background_survive_turning_points(modes):
    for sol in modes:
        r=sol.summary()
        assert r['scalar_velocity_zero_crossings']>=18
        assert r['max_normalized_Hamilton']<3e-9
        assert r['max_density_constraint_gap']<2e-8
        assert r['max_Friedmann']<2e-9
        assert r['max_energy_balance']<2e-9
        state=sol.sample([0])
        np.testing.assert_allclose(state['Delta'][:,0,0],1,atol=1e-14)
        np.testing.assert_allclose(state['Delta'][:,1,0],0,atol=1e-14)


def test_independent_clocks_agree_and_phase_is_identifiable_input(modes):
    a,b=(sol.summary() for sol in modes)
    assert abs(a['proper_time']-b['proper_time'])<1e-8
    for x,y in zip(a['rows'],b['rows']):
        for key in ('density_seed_final_Delta','phase_seed_final_Delta','density_seed_potential_ratio'):
            assert abs(x[key]-y[key])<2e-8
    assert a['rows'][1]['density_seed_final_Delta']>1.15
    assert a['rows'][2]['density_seed_final_Delta']<.93
    # Same initial density and geometry but a changed scalar velocity gives
    # a distinct later density. No unique preparation spectrum is invented.
    assert a['rows'][2]['phase_seed_final_Delta']>1.36


def test_initial_scalar_velocity_never_appears_in_a_denominator(modes):
    sol=modes[0]
    row=sol.sample(np.linspace(0,sol.solution.t[-1],201))
    assert row['Pi'][0]==0
    assert all(np.isfinite(row[k]).all() for k in ('Phi','X','Y','Delta'))


def test_flat_poisson_substitution_fails_on_the_actual_parent(modes):
    sol=modes[0];r=sol.sample([0]);kap=sol.model.initial['kappa']
    eigenvalue=sol.modes[0]*(sol.modes[0]+2)
    flat_density=-2*eigenvalue*r['Phi'][0,0,0]/(kap*r['a'][0]**2*r['E'][0])
    assert r['Delta'][0,0,0]==pytest.approx(1)
    assert flat_density==pytest.approx(1.6)  # n=2: 8/(8-3), not unity


def test_radial_eigenfunctions_and_packet_initial_constraints():
    p=ParentPacket();theta=np.linspace(.1,1.3,8)
    # Independent trigonometric form; Chebyshev implementation remains
    # regular at the poles where the trigonometric ratio is 0/0.
    for n in (2,8,32):
        q=p.radial_harmonics([n],theta)[0]
        np.testing.assert_allclose(q,np.sin((n+1)*theta)/((n+1)*np.sin(theta)),atol=2e-14)
    assert p.radial_harmonics([2,8,32],[0]).tolist()==[[1.],[1.],[1.]]
    profile=p.profiles(theta);pars=p.model.initial
    expected=2*(profile['laplacian']+3*profile['Phi'])/(pars['kappa']*pars['a']**2)
    np.testing.assert_allclose(expected,pars['omega']**2*pars['chi']*profile['X'],atol=2e-14)
    assert p.profiles([0])['Delta'][0]==pytest.approx(1,abs=1e-12)
    for value in p.profiles([p.support_angle,2.,2.85]).values():assert np.all(value==0)


def test_compensation_and_domain_have_independent_integral_checks():
    p=ParentPacket();coeff=p.coefficients(64)
    assert max(abs(coeff[:2]))<2e-12
    for n in (0,1,8):
        val=quad(lambda theta:float(p.profiles(theta)['Phi'])*math.sin(theta)*math.sin((n+1)*theta),
            0,p.support_angle,epsabs=1e-13)[0]*2*(n+1)/math.pi
        assert val==pytest.approx(coeff[n],abs=2e-12)
    d=p.model.atlas.parent_causal_domain()
    # Even the two-way total conformal allowance fits between initial
    # packet support and the wall: its whole causal cone stays protected.
    assert p.support_angle+2*d['future_conformal_distance_upper_bound']<d['initial_parent_angle']


def test_spatial_resolution_growth_and_fixed_metric_control():
    p=ParentPacket();theta=np.linspace(0,2.85,401);outputs=[]
    for cutoff in (64,128):
        sol,coeff,initial=p.evolve(cutoff)
        result=p.reconstruct(sol,initial,sol.solution.t[-1],theta)
        outputs.append(result)
        first=p.reconstruct(sol,initial,0,theta)
        assert max(abs(first['Delta']-p.profiles(theta)['Delta']))<4e-7
        assert max(abs(result['Delta'][theta>1.67]))<2e-8
    assert max(abs(outputs[0]['Delta']-outputs[1]['Delta']))<4e-7
    assert outputs[1]['Delta'][0]==pytest.approx(1.15188159,abs=2e-8)
    control=fixed_metric_control(sol)
    center=float(initial@control['Delta'])
    assert .98<center<.99
    assert outputs[1]['Delta'][0]-center>.16


def test_jeans_comparison_is_a_scale_condition_not_a_new_source():
    m=ParentModes();r=jeans_comparison(m);p=m.initial
    assert r['initial_excitation_to_vacuum_ratio']==pytest.approx(4/3)
    assert r['initial_Jeans_wavelength_over_vacuum_radius']==pytest.approx(.4622588056664)
    lam=p['kappa']*p['vacuum'];xi=r['initial_Jeans_wavelength_over_vacuum_radius']
    assert p['omega']/math.sqrt(lam)==pytest.approx(math.pi**2*math.sqrt(2/3)/xi**2)


@pytest.mark.parametrize('indices',[(1,),(2,2),(2.5,),(257,),()])
def test_admitted_harmonics_are_explicit(indices):
    with pytest.raises(ValueError):ParentModes().solve(indices)


def test_packet_cannot_silently_reach_the_unconstructed_wall_future():
    with pytest.raises(ValueError):ParentPacket(support_angle=2.4)
    with pytest.raises(ValueError):ParentPacket(power=6.5)
