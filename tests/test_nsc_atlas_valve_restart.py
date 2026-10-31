"""Geometric and admission controls for the local same-state restart."""
import math

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_valve_spherical import (
    SphericalValve,ConformalValve,ConformalRestart,align_conformal_center,
    conformal_metric_diagnostics,align_conformal_radius,
)


@pytest.fixture(scope='module')
def supplied_profile():
    with threadpool_limits(limits=1):
        system=SphericalValve(48);sol=system.evolve(.3)
        state=sol.y[:,-1].reshape(8,48)
        profile=ConformalValve(system,state,initial_proper_time=.3).initial_grid(1001)
    profile.update(T=0.,center_proper_time=.3)
    return system,profile


def test_constant_conformal_scale_retains_fields_and_physical_momenta(supplied_profile):
    system,profile=supplied_profile
    frames=[ConformalRestart(profile,H=system.H,kappa=1.,sstar=system.sstar,
                              half_width=.4*c,coordinate_scale=c) for c in (1.,2.)]
    a,b=(f.initial_grid(501) for f in frames)
    np.testing.assert_allclose(b['x'],2*a['x'],rtol=0,atol=0)
    np.testing.assert_allclose(b['y'][[0,2,3]],a['y'][[0,2,3]],rtol=0,atol=0)
    np.testing.assert_allclose(b['y'][1],a['y'][1]-math.log(2),atol=2e-16,rtol=0)
    np.testing.assert_allclose(np.exp(-b['y'][1])*b['yt'],np.exp(-a['y'][1])*a['yt'],atol=2e-15,rtol=0)
    with threadpool_limits(limits=1):
        runs=[f.evolve(501,endpoint=.03*c,sample_every=.005*c) for f,c in zip(frames,(1.,2.))]
    assert all(r['final_cut_admitted'] and r['status']=='requested_endpoint' for r in runs)
    values=[align_conformal_center(r,.31)['values'] for r in runs]
    for key in ('center_radius','center_Rproper','center_Rproper_acc','center_S','center_Pi_S'):
        assert abs(values[0][key]-values[1][key])<2e-10


def test_direct_metric_curvature_is_independently_correct_for_de_sitter():
    system=SphericalValve(24,epsilon=0.);H=system.H;T=.3;x=np.linspace(-.7,.7,1001)
    scale=1/math.cos(H*T);R=scale*np.cos(H*x)/H;a=-math.log(math.cos(H*T));aT=H*math.tan(H*T);zero=np.zeros_like(x)
    y=np.array([R,np.full_like(x,a),zero,zero]);yt=np.array([aT*R,np.full_like(x,aT),zero,zero])
    dx=np.array([-scale*np.sin(H*x),zero,zero,zero]);dxx=np.array([-H*H*R,zero,zero,zero])
    dtt=np.array([(H*H*scale*scale+aT*aT)*R,np.full_like(x,H*H*scale*scale),zero,zero])
    r=conformal_metric_diagnostics(T,0.,y,yt,dx,dxx,x[1]-x[0],H,1.,system.sstar,ytt=dtt)
    assert r['center_metric_Ricci_over_H2']==pytest.approx(12.,abs=5e-14)
    assert r['center_metric_Kretschmann_over_H4']==pytest.approx(24.,abs=8e-14)
    assert r['center_metric_Einstein_residual_over_source']<1e-14
    assert r['metric_Einstein_residual_over_source_RMS']<2e-10
    assert r['independent_temporal_metric_jet']
    # Holding the fields and their first jets fixed must not force curvature
    # to agree with their source when the metric's temporal jet changes.
    perturbed=dtt.copy();perturbed[1]+=.01
    changed=conformal_metric_diagnostics(T,0.,y,yt,dx,dxx,x[1]-x[0],H,1.,system.sstar,ytt=perturbed)
    assert changed['metric_Einstein_residual_over_source_RMS']>1e-3
    assert abs(changed['center_metric_Ricci_over_H2']-r['center_metric_Ricci_over_H2'])>.1
    initialized=conformal_metric_diagnostics(T,0.,y,yt,dx,dxx,x[1]-x[0],H,1.,system.sstar)
    assert not initialized['independent_temporal_metric_jet']


def test_failed_guards_cannot_become_admitted_by_diagnostic_continuation(supplied_profile):
    system,profile=supplied_profile
    frame=ConformalRestart(profile,H=system.H,kappa=1.,sstar=system.sstar,half_width=.4)
    with threadpool_limits(limits=1):
        stopped=frame.evolve(251,endpoint=.025,sample_every=.005,relative_source_limit=1e-14)
        diagnostic=frame.evolve(251,endpoint=.025,sample_every=.005,relative_source_limit=1e-14,
                                 diagnostic_continue_after_guard=True)
    assert stopped['status']=='numerical_validity_stop'
    assert diagnostic['status']=='requested_endpoint'
    assert not stopped['final_cut_admitted'] and not diagnostic['final_cut_admitted']
    assert stopped['first_validity_failure'] is not None
    assert diagnostic['first_validity_failure'] is not None
    radius=float(diagnostic['center_history'][-2,2])
    assert not align_conformal_radius(diagnostic,radius)['admitted']


def test_restart_patch_and_radius_alignment_reject_unavailable_information(supplied_profile):
    system,profile=supplied_profile
    with pytest.raises(ValueError,match='trusted profile'):
        ConformalRestart(profile,H=system.H,kappa=1.,sstar=system.sstar,half_width=100.)
    frame=ConformalRestart(profile,H=system.H,kappa=1.,sstar=system.sstar,half_width=.4)
    with threadpool_limits(limits=1):run=frame.evolve(251,endpoint=.02,sample_every=.005)
    with pytest.raises(ValueError,match='outside'):align_conformal_radius(run,.1)
