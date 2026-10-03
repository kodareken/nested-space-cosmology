"""AP provenance, own constraints, original RHS and guarded short advancement."""
import json
from pathlib import Path
import time

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_discovery_turn_fullfield as full


@pytest.fixture
def preparation():
    path = full.leading.INITIAL
    saved = json.loads(path.read_text())
    payload = path.with_name(saved["payload"])
    a = -8*np.pi*saved["locked_coefficients"]["A"]
    return {"input_record": str(path), "input_payload": str(payload),
        "source_case": "uniform", "Q0": saved["cases"]["uniform"]["declaration"]["mean_Q"],
        "period": 8., "kappa": 1., "multiplicity": 4,
        "initial_momenta_zero": True, "held_pair_occupations": [.752,.752,.5,.5,.248,.248],
        "input_hashes": {str(path): full.digest(path), str(payload): full.digest(payload)},
        "constants": {"a": a, "Z": 3*a,
            "mag": 2*np.pi*saved["locked_coefficients"]["C_F"]*saved["locked_coefficients"]["flux"]**2,
            **{name: saved["locked_coefficients"][name] for name in ("A", "C_F", "flux")}}}


@pytest.fixture
def prepared(preparation):
    with threadpool_limits(limits=1), full.leading.backend.fft_thread_limit(1):
        return full.prepare_source(preparation, (.752,.5,.248))


def test_actual_AP_fold_preserves_normalization_phase_and_density(preparation, prepared):
    pair, state, metadata = prepared
    with np.load(preparation["input_payload"], allow_pickle=False) as arrays:
        original0 = arrays["uniform_nodal_phi0"]
        original1 = arrays["uniform_nodal_phi1"]
    nodal = full.leading.decode(pair, state)
    lift = full.leading.backend.SpinorCarrier(32,128,8.,canonical=True)
    np.testing.assert_allclose(lift@nodal.phi0, original0, atol=2e-13, rtol=0)
    np.testing.assert_allclose(lift@nodal.phi1, original1, atol=2e-13, rtol=0)
    gram = nodal.phi0.conj().T@nodal.phi0+nodal.phi1.conj().T@nodal.phi1
    np.testing.assert_allclose(gram,np.eye(6),atol=3e-14,rtol=0)
    assert metadata["frame"]["discarded_tail_frobenius"] < 1e-12
    assert not metadata["frame"]["orthonormalization_applied"]
    assert not metadata["frame"]["new_eigenframe_selected"]
    assert metadata["multiplicity_applied_once_outside_H"]
    assert metadata["source_trace"] == pytest.approx(3.)
    assert metadata["initial_constraints"]["raw_C_max"] < 1e-8
    assert metadata["initial_constraints"]["D_max"] < 1e-10
    assert metadata["r0"]**2 == pytest.approx(
        -(metadata["constants"]["mag"]+metadata["actual_mean_rho"]/metadata["Q0"])/metadata["constants"]["a"])


def test_passive_clock_work_step_is_original_full_RK4_and_scaled_control(prepared):
    pair, initial, _metadata = prepared
    with threadpool_limits(limits=1), full.leading.backend.fft_thread_limit(1):
        dt, restriction = full.control.step_restriction(pair,initial,.001)
        assert restriction["local_real_variables"] == 28
        assert restriction["moving_scale_sign_cancellation_used"] is False
        actual, clocks, work = full.rk4_with_ledgers(pair, initial, dt)
        expected = full.leading.rk4_step(pair,initial,dt)
        for name in full.leading.FIELDS:
            np.testing.assert_allclose(getattr(actual,name),getattr(expected,name),atol=3e-16,rtol=1e-14)
        assert np.all(clocks > 0) and np.isfinite(work)
        record, _arrays = full._advance(pair,initial,.001,deadline=time.process_time()+5.,
                                        maximum_time=.002,before_advance=lambda: None)
    assert record["time"] == .002 and record["event_root"] is None
    assert record["armed_after_positive_p_Q"]
    assert record["original_full_RHS"].endswith("nsc_discovery_leading_einstein")
    assert not record["reduced_ODE_imported"]
    assert record["diagnostics"]["Gram_max"] < 1e-10


def test_zero_initial_maximum_is_not_a_turn_and_partial_RK4_root(monkeypatch):
    assert not full.positive_turn_bracket(0.,-1.,False)
    assert not full.positive_turn_bracket(0.,-1.,True)
    assert not full.positive_turn_bracket(-1.,1.,False)
    assert full.positive_turn_bracket(1.,0.,True)
    monkeypatch.setattr(full,"momentum",lambda pair, state: state)
    # Manufactured full-step oracle with p(t)=.004-t; no actual trajectory.
    monkeypatch.setattr(full,"rk4_with_ledgers",lambda pair,state,dt:
                        (state-dt,np.full(3,dt),2*dt))
    monkeypatch.setattr(full.leading,"energy",lambda pair,state: {"field": state})
    state, fraction, clocks, work, root = full.root_turn(None,.004,.006,
                                                        deadline=time.process_time()+1.)
    assert abs(state) < 1e-12 and fraction == pytest.approx(.004,abs=1e-12)
    np.testing.assert_allclose(clocks,.004,atol=1e-12)
    assert work == pytest.approx(.008,abs=2e-12)
    assert root["root_time_width"] < 2e-12 and root["partial_RK4_root"]
    assert not root["sampled_minimum_substituted"]


def test_lock_rejection_blocks_callbacks_and_advancement(tmp_path, prepared, monkeypatch):
    path = tmp_path/"unlocked.json"
    path.write_text(json.dumps({"schema":full.SCHEMA,"mode":"prepared","locked_before_held_data":False}))
    called=[]
    monkeypatch.setattr(full,"prepare_source",lambda *args: called.append("prepare"))
    with pytest.raises(ValueError,match="LOCKED"):
        full.measure_locked_prediction(path)
    assert not called
    pair,state,_metadata=prepared
    monkeypatch.setattr(full.control,"step_restriction",lambda *args: called.append("step"))
    def reject():
        raise ValueError("missing authenticated lock")
    with pytest.raises(ValueError,match="authenticated"):
        full._advance(pair,state,.001,deadline=time.process_time()+1.,before_advance=reject)
    assert not called
    with pytest.raises(ValueError,match="lock guard"):
        full._advance(pair,state,.001,deadline=time.process_time()+1.)


def test_current_producer_and_input_drift_rejected_before_source(preparation,tmp_path,monkeypatch):
    names=(Path(full.__file__),Path(full.leading.__file__),Path(full.control.__file__))
    producers={str(path.relative_to(full.ROOT)):full.digest(path) for path in names}
    record={"schema":full.SCHEMA,"mode":"locked_prediction","locked_before_held_data":True,
            "preparation":preparation,"producers":producers}
    path=tmp_path/"locked.json";path.write_text(json.dumps(record))
    authenticated,binding=full.authenticate_lock(path)
    assert authenticated==record and binding["prediction_sha256"]==full.digest(path)
    preparation["input_hashes"][preparation["input_payload"]]="0"*64
    path.write_text(json.dumps(record))
    called=[];monkeypatch.setattr(full,"prepare_source",lambda *args:called.append(True))
    with pytest.raises(ValueError,match="source input changed"):
        full.measure_locked_prediction(path)
    assert not called
