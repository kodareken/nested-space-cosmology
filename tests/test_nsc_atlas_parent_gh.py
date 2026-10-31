"""Integrated regular-frame jets, exact vacuum and rejected-stage controls."""
import math

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_parent_gh import ParentGH,PAIRS,FRAME_MATRICES,evolve_parent


def test_packed_component_symmetry_and_raw_ordered_jets():
    with threadpool_limits(limits=1):
        system=ParentGH(2,0);y=system.initial()
        points=np.vstack([system.basis.points,system.X0]);x,z,_,_=points.T
        values=np.zeros((len(points),24))
        field=PAIRS.index((1,2));velocity=10+PAIRS.index((0,1))
        values[:,field]=.01*x*z
        values[:,velocity]=.02*z
        values[:,20]=.03*points[:,2]
        addition=system.project(values,system.chars)
        y[:system.dof]+=system.pack(addition)
        J=system.jets(y);r=system.r
        np.testing.assert_allclose(J.g[:,1,2],.01*x*z,atol=1e-15,rtol=0)
        gradient=np.c_[z,x,np.zeros((len(x),2))]
        hess=np.zeros((4,4));hess[0,1]=hess[1,0]=1
        directions=np.einsum('iab,nb->nia',FRAME_MATRICES,points)
        for i in range(3):
            expected=.01*r*np.einsum('na,na->n',directions[:,i],gradient)
            np.testing.assert_allclose(J.dg[:,i,1,2],expected,atol=2e-15,rtol=0)
            np.testing.assert_allclose(J.dv[:,i,0,1],.02*r*directions[:,i,1],atol=2e-15,rtol=0)
            for j in range(3):
                acceleration=points@(FRAME_MATRICES[j]@FRAME_MATRICES[i]).T
                expected=.01*r*r*(np.einsum('na,ab,nb->n',directions[:,i],hess,directions[:,j])+np.einsum('na,na->n',acceleration,gradient))
                np.testing.assert_allclose(J.d2g[:,i,j,1,2],expected,atol=2e-15,rtol=0)
        assert np.all(np.isfinite(J.q))


def test_exact_vacuum_reaches_full_reference_proper_endpoint():
    result=evolve_parent(degree=2,epsilon=0,threads=1,budget_seconds=90,gauge='vacuum_reference')
    assert result['status']=='endpoint_reached'
    from recursive_horizons.nsc_atlas_parent_circulation import Parent
    H=Parent.retained().H
    for point in result['snapshots']:
        exact=math.cosh(H*point['observer_proper_duration'])**3
        assert abs(point['proper_volume_over_initial_round_reference']/exact-1)<1e-9
        assert point['ADM_H_RMS_over_H_squared']<1e-10
        assert point['ADM_M_RMS_over_H_squared']<1e-10
        assert point['GH_C_RMS_over_H']<1e-10
        assert point['observer_unit_norm_error']<1e-10


def test_nonfinite_trial_is_rejected_without_repairing_fields():
    with threadpool_limits(limits=1):
        system=ParentGH(2,0);y=system.initial();y[0]=np.nan
        returned=system.rhs(.1,y)
    assert np.all(np.isnan(returned))
    assert np.isnan(y[0]) and system.invalid_stages==1


def test_proper_frame_is_not_the_scalar_kernel_free_elliptic_space():
    with threadpool_limits(limits=1):system=ParentGH(2,0)
    selected=np.all(system.basis.characters==[1,1],axis=1)
    assert any(label[0]==1 for label in np.array(system.basis.labels,dtype=object)[selected])


def test_harmonic_gauge_uses_correct_proper_time_and_coordinate_map():
    result=evolve_parent(degree=2,epsilon=0,threads=1,budget_seconds=90,
        gauge='harmonic',rtol=3e-11,atol=3e-13)
    from recursive_horizons.nsc_atlas_parent_circulation import Parent
    H=Parent.retained().H
    for point in result['snapshots']:
        u=H*point['observer_proper_duration']
        exact_t=(math.tanh(u)/math.cosh(u)+math.atan(math.sinh(u)))/(2*H)
        assert abs(point['coordinate_t']-exact_t)<1e-8
        assert abs(point['minimum_lapse']/math.cosh(u)**3-1)<1e-8
        assert point['ADM_H_RMS_over_H_squared']<2e-9
        assert point['GH_C_RMS_over_H']<2e-9


def test_field_phase_readout_does_not_use_an_auxiliary_angle():
    from recursive_horizons.nsc_atlas_finite_parent_data import source_profile
    with threadpool_limits(limits=1):
        system=ParentGH(2,0);y=system.initial()
        f,_=source_profile(system.basis.points,system.parent)
        values=np.zeros((system.N,1));values[:,0]=f
        seed=system.project(values,system.chars[21:22])[0]
        for phase in (.2,1.1,2.8):
            C=system.unpack(y);C[21]=math.cos(phase)*seed
            C[23]=-system.parent.reference_frequency*math.sin(phase)*seed
            y[:system.dof]=system.pack(C);y[-1]=12345.
            assert system.observer_angle(y)==pytest.approx(phase,abs=1e-13)


def test_reference_vacuum_gauge_has_finite_proper_window_on_static_geometry():
    import sympy as s
    t,H=s.symbols('t H',positive=True)
    N=s.cosh(H*t)**-3
    tau=(s.tanh(H*t)/s.cosh(H*t)+s.atan(s.sinh(H*t)))/(2*H)
    assert s.simplify(s.diff(tau,t)-N)==0
    assert s.simplify(s.diff(s.log(N),t)+3*H*s.tanh(H*t))==0
    assert s.limit(tau,t,s.oo)==s.pi/(4*H)


def test_fast_hopf_and_dense_rhs_agree_for_one_resolved_state():
    from numba import set_num_threads
    set_num_threads(1)
    with threadpool_limits(limits=1):
        dense=ParentGH(4,.001,backend='dense');state=dense.initial()
        fast=ParentGH(4,.001,backend='hopf')
        assert fast.dof==dense.dof and fast.basis.labels==dense.basis.labels
        expected=dense.rhs(0.,state);actual=fast.rhs(0.,state)
    np.testing.assert_allclose(actual,expected,atol=2e-10,rtol=0)
