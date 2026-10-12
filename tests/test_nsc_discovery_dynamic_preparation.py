"""Owned general-Q Jacobian, source, exact constraints and checkpoint checks."""
from pathlib import Path
import tempfile

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_discovery_dynamic_preparation as prep
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons import nsc_spherical_coupling as coupling


@pytest.fixture(scope="module")
def prepared():
    with threadpool_limits(limits=1):
        return prep.build_record(execute=True,nf=64,harmonic=2,cpu_limit=20.)


def test_general_q_jacobian_is_owned_centered_difference():
    with threadpool_limits(limits=1):
        pair=prep.make_pair(32);grid=pair.grid
        Q,_=prep.declared_q(grid,harmonic=2)
        state,_=prep.spectral_state(pair,Q)
        state.r=1.4+.15*np.cos(2*np.pi*grid.xi_g/grid.length)
        fine=galerkin.prolong_state(grid,state)
        source=coupling.source_from_columns(galerkin.active_fine_system(grid,fine),fine)
        rho=source["force_L"]/grid.dx_q
        _,_,J=prep.radius_residual_jacobian(grid,state,rho)
        direction=np.cos(4*np.pi*grid.xi_g/grid.length)+.3
        images=[]
        for sign in (-1,1):
            varied=state.copy();varied.r+=sign*1e-6*direction
            images.append(prep.radius_residual_jacobian(grid,varied,rho)[0])
        np.testing.assert_allclose((images[1]-images[0])/2e-6,J@direction,rtol=2e-7,atol=2e-6)


def test_convex_seed_then_exact_constraint_and_current(prepared):
    report,arrays=prepared
    assert report["status"]=="FINITE_INITIAL_DATA"
    assert report["source_eta"]==1. and report["source_trace"]==3.
    for case in ("uniform","pattern"):
        measured=report["cases"][case]
        assert measured["converged"]
        assert measured["full_lapse_max"]<1e-6
        assert measured["full_shift_max"]<1e-10
        assert abs(measured["source_current_mean"])<1e-10
        assert measured["rho_min"]>0 and measured["r_min"]>0
        assert measured["spectral"]["commutator_max"]<1e-10
        assert not measured["holding_claim"] and not measured["evolution_performed"]
        assert np.all(arrays[case+"_nodal_chi"]==0)
        for key in ("p_Q","p_r","p_chi"):
            assert np.all(arrays[case+"_nodal_"+key]==0)
        phi=np.vstack((arrays[case+"_nodal_phi0"],arrays[case+"_nodal_phi1"]))
        C=(phi*prep.WEIGHTS)@phi.conj().T
        ev=np.linalg.eigvalsh(C)
        assert np.min(ev)>-1e-12 and np.max(ev)<1+1e-12
        np.testing.assert_allclose(ev[-6:],np.sort(prep.WEIGHTS),atol=1e-12)
    assert report["cases"]["uniform"]["declaration"]["mean_Q"]==report["cases"]["pattern"]["declaration"]["mean_Q"]
    assert not report["evidence_written"]


def test_preflight_and_checkpoint_are_explicit(prepared):
    report,arrays=prep.build_record()
    assert report["status"]=="PREPARATION_PREFLIGHT" and arrays=={}
    report,arrays=prepared
    with tempfile.TemporaryDirectory(dir="/tmp",prefix="nsc-dynamic-prep-") as directory:
        stem=Path(directory)/"checkpoint"
        prep.write_record(report,arrays,stem)
        before={p.name:p.read_bytes() for p in Path(directory).iterdir()}
        checked=prep.check_record(stem)
        assert checked["ok"] and checked["bytes_written"]==0 and not checked["solver_executed"]
        assert before=={p.name:p.read_bytes() for p in Path(directory).iterdir()}
        pair,state,_=prep.load_case(stem)
        assert pair.grid.nf==64 and state.phi0.shape==(64,6)
        with pytest.raises(FileExistsError):
            prep.write_record(report,arrays,stem)


def test_budget_and_positive_input_domain():
    with pytest.raises(ValueError,match="CPU"):
        prep.build_record(cpu_limit=31.)
    pair=prep.make_pair(32)
    with pytest.raises(ValueError,match="positive lift"):
        prep.declared_q(pair.grid,lift=0)


def test_real_unequal_weight_rotation_is_physical_coherence(prepared):
    _,arrays=prepared
    with threadpool_limits(limits=1):
        pair=prep.make_pair(64)
        Q=arrays["pattern_nodal_Q"]
        columns=(arrays["pattern_nodal_phi0"],arrays["pattern_nodal_phi1"])
        state,measured,_,_,_=prep.prepare_state(pair,Q,rotation_angle=np.pi/6,columns_override=columns)
    old=np.vstack(columns);new=np.vstack((state.phi0,state.phi1))
    np.testing.assert_allclose(new[:,0],np.cos(np.pi/6)*old[:,0]+np.sin(np.pi/6)*old[:,4],atol=1e-13)
    np.testing.assert_allclose(new.conj().T@new,np.eye(6),atol=1e-12)
    assert measured["converged"] and measured["full_lapse_max"]<1e-6
    assert measured["full_shift_max"]<1e-10
    assert measured["spectral"]["source_covariance_change_frobenius"]==pytest.approx(np.sqrt(2)/4)
    assert measured["spectral"]["commutator_max"]>1e-3
    assert not measured["stationary_spectral_control"] and not measured["rotation_is_charge_conjugation"]
    assert measured["fine_eigen_tail_max"] is None


def test_thin_episode_adapter_uses_owned_step_and_clocks(prepared,monkeypatch):
    from recursive_horizons import nsc_discovery_episode as episode
    from recursive_horizons import nsc_nested_parent_child as nested
    _,arrays=prepared;pair=prep.make_pair(64)
    nodal=coupling.CauchyState(**{key:arrays["pattern_nodal_"+key].copy() for key in nested.STATE_NAMES})
    state=nested.encode_state(pair,nodal)
    def test_row(pair,state,time_value,control_mode="coupled"):
        return {"time":time_value,"source_rho":{"min":1.},
                **{name:0. for name in episode.WORK_CHANNELS},"observation_cpu_seconds":0.}
    monkeypatch.setattr(prep,"_pinned_scalar_row",test_row)
    owner=episode.advance_case
    with threadpool_limits(limits=1),prep.episode_adapter():
        result=episode.advance_case(pair,state,coordinate_time=0.,steps=0,stations=(.3,),step_cap=.0005,
            geometry_frozen=False,cpu_allowance=10.,spent=0.,forecast_factor=1.5,
            memory_limit_bytes=2**30,max_steps=1,backend="dense",normal_clocks=np.zeros(3))
        expected=nested.rk4_step(pair,state,result["time"])
    assert episode.advance_case is owner
    assert result["status"]=="step_limit" and result["steps"]==1
    for key in nested.STATE_NAMES:
        np.testing.assert_allclose(getattr(result["state"],key),getattr(expected,key),atol=1e-12)
    expected_clock=.5*result["time"]*(episode._clock_rates(pair,state)+episode._clock_rates(pair,expected))
    np.testing.assert_allclose(result["normal_clocks"],expected_clock,atol=1e-14)


def test_historical_observed_binding_does_not_heal_record(prepared,monkeypatch):
    import json
    report,arrays=prepared
    with tempfile.TemporaryDirectory(dir="/tmp",prefix="nsc-dynamic-history-") as directory:
        stem=Path(directory)/"initial"
        prep.write_record(report,arrays,stem)
        jp,payload=prep.output_paths(stem)
        before=jp.read_bytes()
        observed={"record_sha256":prep.sha256(jp),"payload_sha256":prep.sha256(payload),
                  "producing_commit":"external-observed-test-commit"}
        called=[]
        monkeypatch.setattr(prep,"_git_hashes",lambda commit,hashes:called.append((commit,hashes)) or commit)
        authenticated=prep.authenticate_preparation(stem,observed_binding=observed)
        assert authenticated["external_post_run_observation"]
        assert called[0][0]==observed["producing_commit"] and jp.read_bytes()==before
        observed["record_sha256"]="bad"
        with pytest.raises(ValueError,match="exact preparation"):
            prep.authenticate_preparation(stem,observed_binding=observed)


def test_six_cases_fork_exact_state_and_existing_worker_preview(prepared,monkeypatch):
    from recursive_horizons import nsc_discovery_episode as episode
    report,arrays=prepared
    monkeypatch.setattr(prep,"_git_hashes",lambda commit,hashes:commit)
    with tempfile.TemporaryDirectory(dir="/tmp",prefix="nsc-dynamic-episode-") as directory:
        original=Path(directory)/"v1";successor=Path(directory)/"v2";campaign=Path(directory)/"episode"
        prep.write_record(report,arrays,original)
        with threadpool_limits(limits=1):
            new,new_arrays=prep.build_coherent_successor(original,execute=True,producer_commit="test-frozen-producer")
        assert new["cases"]["coherent"]["converged"]
        prep.write_record(new,new_arrays,successor)
        manifest=prep.prepare_episode(successor,campaign,execute=True,producer_commit="test-frozen-producer",stations=(.3,))
        assert len(manifest["cases"])==6 and manifest["evolved"] is False
        for name in ("uniform","pattern","coherent"):
            a,aa=episode.load_checkpoint(campaign,f"nf64_{name}_dt0.001")
            b,bb=episode.load_checkpoint(campaign,f"nf64_{name}_dt0.0005")
            assert a["source_pins"]["initial_state_sha256"]==b["source_pins"]["initial_state_sha256"]
            for key in aa:
                np.testing.assert_array_equal(aa[key],bb[key])
        with threadpool_limits(limits=1),prep.episode_adapter():
            result=prep._dynamic_case_worker({"directory":str(campaign),"case_id":"nf64_coherent_dt0.0005",
                "cpu_budget_seconds":30.,"forecast_factor":1.5,"memory_limit_bytes":2**30,"max_steps":1,"backend":"dense"})
        assert result["status"]=="step_limit" and result["steps"]==1
        latest,_=episode.load_checkpoint(campaign,"nf64_coherent_dt0.0005")
        assert latest["coordinate_time"]>0 and latest["initial_state_called"] is False
        assessment=prep.assess_episode(campaign)
        assert assessment["cases"]["nf64_coherent_dt0.0005"]["evolved_effect_measured"]
        assert not assessment["startup_force_is_holding"]


def test_source_sign_event_is_retained_without_stopping(monkeypatch):
    def fake_owner(pair,state,**options):
        time=options["until_time"]
        first=options["coordinate_time"]==0
        rows=([{"time":0.,"source_rho":{"min":1.}}] if first else [])
        rows.append({"time":time,"source_rho":{"min":-1.}})
        clocks=np.full(3,time)
        snap={"time":time,"kind":"station" if time>=.3 else "event","channel_sample":rows[-1].copy(),
              "stations_reached":[],"state":state,"clocks":clocks,"clock_rates":np.ones(3)}
        return {"state":state,"time":time,"steps":options["steps"]+1,
            "status":"station_reached" if time>=.3 else "until_time", "stop":{},
            "normal_clocks":clocks,"clock_rates":np.ones(3),"step_cpu_seconds":0.,"observation_cpu_seconds":0.,
            "observations":rows,"snapshots":[snap],"stations_reached":[.3] if time>=.3 else [],
            "work_ledger":{},"channel_sample":rows[-1].copy(),"channel_time":time}
    monkeypatch.setattr(prep,"_advance_owner",fake_owner,raising=False)
    result=prep._event_aware_advance(None,None,coordinate_time=0.,steps=0,stations=(.3,),
        normal_clocks=np.zeros(3),cpu_allowance=10.,spent=0.,observation_cadence=.05)
    assert result["status"]=="station_reached" and result["time"]==pytest.approx(.3)
    event=next(s for s in result["snapshots"] if s["time"]==pytest.approx(.05))
    assert event["channel_sample"]["first_source_density_sign_change"]["automatic_stop"] is False
    np.testing.assert_array_equal(result["normal_clocks"],np.full(3,.3))
