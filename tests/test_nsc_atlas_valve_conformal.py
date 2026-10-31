"""Same-state, causal-domain and independent analytic frame controls."""
import math

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_valve_spherical import (
    SphericalValve,ConformalValve,align_conformal_center,
    conformal_forcing,conformal_de_sitter_control,_conformal_evolve,
)


@pytest.fixture(scope='module')
def trusted():
    with threadpool_limits(limits=1):
        system=SphericalValve(48);solution=system.evolve(.5)
    return system,solution,solution.sol(.3).reshape(8,48)


def test_handoff_preserves_actual_fields_momenta_and_unit_radial_metric(trusted):
    system,solution,state=trusted
    frame=ConformalValve(system,state,initial_proper_time=.3)
    initial=frame.initial_grid(1001);z=system.pieces(state,initial['u'])
    np.testing.assert_allclose(initial['y'][2:],z['y'][4:6],rtol=0,atol=0)
    np.testing.assert_allclose(initial['yt'][2:],z['y'][6:8],rtol=0,atol=0)
    np.testing.assert_allclose(initial['yt'][1],-z['ka'],rtol=0,atol=0)
    np.testing.assert_allclose(initial['yt'][0],-z['kb']*z['radius'],rtol=0,atol=0)
    u=initial['u'];h=initial['step']
    du_dx=(u[:-4]-8*u[1:-3]+8*u[3:-1]-u[4:])/(12*h)
    metric=z['q'][2:-2]**2/(1-u[2:-2]**2)*du_dx**2
    np.testing.assert_allclose(metric[3:-3],1.,rtol=0,atol=3e-8)
    # Caller changes after construction cannot reset the retained Cauchy state.
    passed=state.copy();isolated=ConformalValve(system,passed,initial_proper_time=.3)
    before=isolated.initial_grid(33)['y'].copy();passed[:]=0
    np.testing.assert_array_equal(isolated.initial_grid(33)['y'],before)


def test_actual_center_clock_alignment_and_trapping_survive_frame_change(trusted):
    system,solution,state=trusted;frame=ConformalValve(system,state,initial_proper_time=.3)
    with threadpool_limits(limits=1):
        result=frame.evolve(1001,.1)
    assert result['status']=='requested_endpoint'
    aligned=align_conformal_center(result,.35)
    reference=system.pieces(solution.sol(.35).reshape(8,48),[0.])
    old_radius=float(reference['radius'][0]);old_theta=-2*float(reference['kb'][0])/system.H
    assert abs(aligned['values']['center_radius']-old_radius)<2e-7
    assert abs(aligned['values']['center_theta_plus_over_H']-old_theta)<2e-5
    assert aligned['values']['center_theta_plus_over_H']<0
    assert aligned['interpolation_indicator']['center_radius']<2e-7
    assert result['rows'][-1]['center_R_gradient_norm_squared']<0
    with pytest.raises(ValueError,match='outside'):align_conformal_center(result,100.)


def test_exact_de_sitter_pde_factors_and_three_grid_field_convergence():
    system=SphericalValve(24,epsilon=0.);H=system.H;T=.3;x=np.linspace(-.7,.7,51)
    scale=1/math.cos(H*T);R=scale*np.cos(H*x)/H;a=-math.log(math.cos(H*T));aT=H*math.tan(H*T)
    zero=np.zeros_like(x)
    y=np.array([R,np.full_like(x,a),zero,zero])
    yt=np.array([aT*R,np.full_like(x,aT),zero,zero])
    yx=np.array([-scale*np.sin(H*x),zero,zero,zero])
    wave=np.array([(H*H*scale*scale+aT*aT+H*H)*R,np.full_like(x,H*H*scale*scale),zero,zero])
    np.testing.assert_allclose(conformal_forcing(y,yt,yx,system.sstar,system.kappa),wave,atol=3e-15,rtol=0)
    with threadpool_limits(limits=1):
        rows=[conformal_de_sitter_control(n,system=system) for n in (251,501,1001)]
    for key in ('max_radius_error','max_metric_log_error'):
        ratios=[rows[i][key]/rows[i+1][key] for i in (0,1)]
        assert all(3.7<r<4.3 for r in ratios)
    assert rows[-1]['max_radius_error']<1e-8
    assert all(math.isfinite(r['Hamiltonian_RMS_over_H2']) for r in rows)


def test_shrinking_diamond_has_no_outer_boundary_influence_before_light_crossing():
    # This is a numerical domain-of-dependence control, not a claimed new
    # constrained Einstein preparation. Its outer perturbation is deliberately
    # outside the tested center's characteristic past.
    system=SphericalValve(24,epsilon=0.);H=system.H;n=501;x=np.linspace(-1.,1.,n)
    y=np.zeros((4,n));y[0]=np.cos(H*x)/H
    v=np.zeros_like(y);dx=np.zeros_like(y);dx[0]=-np.sin(H*x)
    dxx=np.zeros_like(y);dxx[0]=-H*np.cos(H*x)
    initial=dict(x=x,y=y,yt=v,yx=dx,yxx=dxx,step=x[1]-x[0],half_width=1.)
    changed={k:a.copy() if isinstance(a,np.ndarray) else a for k,a in initial.items()}
    changed['y'][0,abs(x)>.8]+=1e-4
    with threadpool_limits(limits=1):
        results=[_conformal_evolve(i,H=H,kappa=1.,sstar=system.sstar,proper_time=0.,endpoint=.2) for i in (initial,changed)]
    assert all(r['status']=='requested_endpoint' for r in results)
    np.testing.assert_array_equal(results[0]['center_history'],results[1]['center_history'])
    assert len(results[0]['end_profile']['x'])<n


def test_invalid_geometry_and_nonreflection_inputs_are_rejected(trusted):
    system,_,state=trusted
    changed=state.copy();changed[5]+=.01
    with pytest.raises(ValueError,match='reflection'):ConformalValve(system,changed,initial_proper_time=.3)
    frame=ConformalValve(system,state,initial_proper_time=.3)
    with pytest.raises(ValueError,match='odd integer'):frame.initial_grid(100)
    with pytest.raises(ValueError,match='half-width'):frame.evolve(101,frame.half_width)
