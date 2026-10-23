import math
import numpy as np
import pytest

from recursive_horizons.nsc_atlas_core_evolution import (
    DiffuseReservoir,EvolutionRun,run_evolution,initial_state,moving_mesh,
    geometry,derivatives,interpolate,rotate)


def test_exact_oscillator_rotation_preserves_the_quadratic_energy():
    x=np.array([.01,-.03]);p=np.array([2.,-3.]);w=2000.
    a,b=rotate(x,p,.3456789,w)
    np.testing.assert_allclose(b*b+w*w*a*a,p*p+w*w*x*x,rtol=3e-15)
    c,d=rotate(a,b,-.3456789,w)
    np.testing.assert_allclose(c,x,atol=1e-16);np.testing.assert_allclose(d,p,atol=1e-14)


def test_moving_chart_keeps_the_core_resolved_and_the_exterior_comoving():
    x=np.linspace(0,4,4001);mesh=moving_mesh(x);A=5.
    r=x+(A-1)*mesh[0];J=1+(A-1)*mesh[1];Jx=(A-1)*mesh[2]
    np.testing.assert_allclose(r[x<=1],x[x<=1],atol=0)
    np.testing.assert_allclose(r[x>=2],A*x[x>=2],atol=1e-14)
    assert np.min(J)>=1
    fx,fxx=derivatives(r*r,.001);F=fx/J;Fp=fxx/J**2-fx*Jx/J**3
    np.testing.assert_allclose(F[:-5],2*r[:-5],atol=2e-9)
    np.testing.assert_allclose(Fp[:-5],2,atol=3e-7)


def test_homogeneous_geometry_is_exact_on_the_moving_mesh():
    r=np.arange(2001)*.002;delta,Pi,bg,m,_=initial_state('homogeneous',r)
    a0=bg[2];bg[2]*=3
    g=geometry(delta,Pi,bg,.002,2000.,1.,moving_mesh(r),a0)
    assert np.max(abs(g[3]))==0
    assert np.max(abs(g[6]-1))==0
    np.testing.assert_allclose(g[4],1-g[8]**2/bg[2]**2,atol=1e-15)


def test_fractional_proper_emission_event_and_independent_photon():
    from scripts.check_nsc_atlas_core_evolution import independent_homogeneous
    # Arrival interpolation is checked at the finest declared time step;
    # the evolving H contains the scalar's fast pressure oscillation.
    c=EvolutionRun(kind='homogeneous',spacing=.002,step=.000025,end=.7,emission_age=.100037)
    r=run_evolution(c);exact=independent_homogeneous(c.emission_age)
    assert r['arrival']['emission_coordinate_time']==pytest.approx(c.emission_age,abs=1e-11)
    assert r['arrival']['receiver_proper_age']==pytest.approx(exact['receiver_proper_age'],abs=2e-8)
    assert r['arrival']['received_to_emitted_frequency']==pytest.approx(exact['received_to_emitted_frequency'],abs=2e-8)
    assert abs(r['geodesic_frequency_gap'])<2e-9
    assert r['maximum_background_constraint']<1e-8


def test_diffuse_preparation_uses_the_same_total_geometric_mass():
    model=DiffuseReservoir();solution=model.solve()
    target=model.core.mass_solution(model.H)[0].y[0,-1]
    assert model.cloud_amplitude==pytest.approx(.00186251874314,abs=2e-12)
    assert solution.solution.sol(model.radius)[0]==pytest.approx(target,abs=2e-12)
    assert abs(solution.solution.y[0,-1])<2e-11
    f0=model.fields(0,solution.depletion)[0]
    assert f0<model.core.field(0)[0]/4


def test_regular_interpolation_and_admission():
    grid=np.arange(20)*.01
    assert interpolate(grid**2,.0025,grid)==pytest.approx(.0025**2,abs=1e-18)
    with pytest.raises(ValueError):run_evolution(EvolutionRun(receiver=5.))
    with pytest.raises(ValueError):run_evolution(EvolutionRun(step=.0002,spacing=.0001))


def test_independent_symbolic_observer_and_mass_reductions():
    from scripts.check_nsc_atlas_core_evolution import identities
    assert set(identities().values())=={'0'}
