"""Manufactured band/constraint/clock checks; no physical confirmation trajectory."""
import time
import numpy as np
import pytest
from threadpoolctl import threadpool_limits
import derive_nsc_discovery_parent_cut_confirmation as confirm
from recursive_horizons import nsc_discovery_parent_episode as episode


@pytest.fixture
def manufactured():
    data=episode.analytic_fixture(nf=16,rank=2)
    pair,state=episode.pair_and_state_from_population(data["populations"][1],common_k=data["common_k"],
                                                     sign=-1,nf=16,length=8.)
    # Standing real AP columns have exactly zero current even after a c1-only
    # occupation change; travelling test waves would make periodicD insoluble.
    x=np.arange(16)/16
    wave=np.column_stack((np.cos(np.pi*x),np.sin(3*np.pi*x)))*np.sqrt(2/16)
    phi0=wave.astype(complex)/np.sqrt(2);phi1=phi0.copy();frame=np.vstack((phi0,phi1))
    pair=confirm.cut.response.make_holder(pair.grid,weights=pair.weights,source_columns=frame,
        reference_columns=frame,geometry_map=pair.geometry_map)
    state=confirm.cut.leading.State(state.Q,state.r,state.p_Q,state.p_r,phi0,phi1)
    problem,_seed,guess,_theta=confirm.cut.problem(pair,state,1.,-1)
    state=confirm.cut.parent._state_from_theta(pair,problem,guess)
    return pair,state


def test_periodic_and_AP_extension_preserves_physical_values_norm_and_order(manufactured):
    pair,state=manufactured
    with threadpool_limits(limits=1),confirm.backend.fft_thread_limit(1):
        lifted,new,report=confirm.extend(pair,state,nf=32)
    nodal=confirm.cut.leading.decode(pair,state);fine=confirm.cut.leading.decode(lifted,new)
    for name in confirm.cut.leading.REAL_FIELDS:
        expected=confirm.periodic_extension(getattr(nodal,name),lifted.grid.ng,8.)
        np.testing.assert_allclose(getattr(fine,name),expected,rtol=1e-13,atol=1e-13)
    points=np.array([.2,1.7,3.,5.,7.8])
    for name in ("phi0","phi1"):
        np.testing.assert_allclose(confirm.cut.ap_values(getattr(new,name),8.,points),
                                   confirm.cut.ap_values(getattr(state,name),8.,points),atol=2e-14,rtol=1e-13)
    assert report["live_Gram_gap"]<1e-13 and report["source_Gram_gap"]<1e-13
    assert not report["source_normalized"] and not report["source_reselected"]
    old=np.vstack((state.phi0,state.phi1));extended=np.vstack((new.phi0,new.phi1))
    np.testing.assert_allclose(extended.conj().T@extended,old.conj().T@old,atol=1e-13)
    source_eig=np.linalg.eigvalsh(old.conj().T@old*pair.weights[None,:])
    new_eig=np.linalg.eigvalsh(extended.conj().T@extended*pair.weights[None,:])
    np.testing.assert_allclose(new_eig,source_eig,atol=1e-13)
    sample=2+.2*np.cos(2*np.pi*np.arange(15)/15)+.1*np.sin(6*np.pi*np.arange(15)/15)
    values=confirm.periodic_extension(sample,31,8.)
    x=np.arange(31)/31
    np.testing.assert_allclose(values,2+.2*np.cos(2*np.pi*x)+.1*np.sin(6*np.pi*x),atol=1e-14)


def test_own_constraints_change_only_momenta_and_c1(manufactured):
    pair,state=manufactured
    with threadpool_limits(limits=1),confirm.backend.fft_thread_limit(1):
        pair,state,_=confirm.extend(pair,state,nf=32)
        before=state.copy();weights=pair.weights.copy();weights[1]*=1.05
        held,prepared,report=confirm.reprepare(pair,state,1.,weights,3.)
    assert held.weights[0]==pair.weights[0] and held.weights[1]==pytest.approx(1.05*pair.weights[1])
    for name in ("Q","r","phi0","phi1"):
        np.testing.assert_array_equal(getattr(prepared,name),getattr(before,name))
    assert report["converged"] and report["hard_anchor"]["satisfied"]
    assert report["actual_constraints"]["raw_C_max"]<1e-6
    assert report["actual_constraints"]["D_max"]<1e-6
    assert report["k_fixed"]==1. and report["Q_r_Phi_fixed"]
    for name in confirm.cut.leading.FIELDS:np.testing.assert_array_equal(getattr(state,name),getattr(before,name))


def test_manufactured_absolute_clock_root_no_extrapolation_and_budget_stop(manufactured):
    pair,state=manufactured
    with threadpool_limits(limits=1),confirm.backend.fft_thread_limit(1):
        target=confirm.cut.response.centre_clock(pair,state)*.001
        endpoint,event=confirm.evolve(pair,state,target,confirm.cut.parent._cpu_time()+3.)
        untouched,partial=confirm.evolve(pair,state,target,confirm.cut.parent._cpu_time()-1.)
    assert event["matched"] and event["tau"]==pytest.approx(target,abs=1e-10)
    assert 0<event["time"]<.003 and event["no_extrapolation"]
    assert not partial["matched"] and partial["stop"]=="CPU_BUDGET_STOP" and partial["steps"]==0
    for name in confirm.cut.leading.FIELDS:np.testing.assert_array_equal(getattr(untouched,name),getattr(state,name))


def test_confirmation_requires_frozen_creation_budget(tmp_path):
    with pytest.raises(ValueError,match="frozen producer"):
        confirm.run(output=tmp_path/"new")
    with pytest.raises(ValueError,match="600"):
        confirm.run(output=tmp_path/"new",producer_commit="future",cpu_budget=601.)
    (tmp_path/"existing").mkdir()
    with pytest.raises(FileExistsError,match="creation-only"):
        confirm.run(output=tmp_path/"existing",producer_commit="future")


def test_readonly_real_forecast_accepts_rounded_station_time():
    preview=confirm.preview()
    assert preview["nf"]==256 and preview["target_tau"]==confirm.TARGET
    assert not preview["evolved"] and preview["bytes_written"]==0
    assert any("measurement-a0.05-T2.25.json" in name for name in preview["input_hashes"])
