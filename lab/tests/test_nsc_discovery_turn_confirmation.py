"""Exact saved handoff and a bounded(.001 elapsed) full-RHS continuation."""
import json
from pathlib import Path
import time

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

import confirm_nsc_discovery_homogeneous_turn as confirm


@pytest.fixture
def source():
    return confirm.authenticate_source(confirm.DEFAULT_SOURCE)


def test_exact_handoff_retains_all_state_arrays_source_W_and_clock(source,monkeypatch):
    locked,measured,arrays,_pins=source
    def forbidden(*args,**kwargs):
        raise AssertionError("source/radius reprepare forbidden")
    monkeypatch.setattr(confirm.full,"prepare_source",forbidden)
    monkeypatch.setattr(confirm.full.leading,"prepare_radius",forbidden)
    with threadpool_limits(limits=1),confirm.full.leading.backend.fft_thread_limit(1):
        pair,state=confirm.reconstruct_handoff(locked,measured,arrays)
    for name in confirm.full.leading.FIELDS:
        np.testing.assert_array_equal(getattr(state,name),arrays["cap1_"+name])
    np.testing.assert_array_equal(pair.geometry_map,arrays["cap1_W"])
    np.testing.assert_array_equal(pair.weights,arrays["cap1_source_weights"])
    assert measured["backend"]["events"][1]["time"] == pytest.approx(.49)
    assert confirm.full.momentum(pair,state)>0
    changed=dict(arrays,cap1_source_weights=np.ones(6)*.5)
    with pytest.raises(ValueError,match="saved source weights"):
        confirm.reconstruct_handoff(locked,measured,changed)


def test_short_resume_preserves_absolute_clock_work_energy_and_step_history(source,monkeypatch):
    locked,measured,arrays,_pins=source
    def forbidden(*args,**kwargs):raise AssertionError("radius/source reset")
    monkeypatch.setattr(confirm.full,"prepare_source",forbidden)
    old=measured["backend"]["events"][1]
    with threadpool_limits(limits=1),confirm.full.leading.backend.fft_thread_limit(1):
        pair,state=confirm.reconstruct_handoff(locked,measured,arrays)
        expected,dc1,dw1=confirm.full.rk4_with_ledgers(pair,state,.0005)
        expected,dc2,dw2=confirm.full.rk4_with_ledgers(pair,expected,.0005)
        row,saved=confirm.continue_arm(pair,state,old,deadline=time.process_time()+2.,
                                      maximum_time=.001,before_advance=lambda:None)
    for name in confirm.full.leading.FIELDS:
        np.testing.assert_allclose(saved[name],getattr(expected,name),atol=2e-15,rtol=1e-14)
    np.testing.assert_allclose(row["normal_clocks"],np.asarray(old["normal_clocks"])+dc1+dc2,atol=1e-15)
    assert row["integrated_coordinate_fieldwork"] == pytest.approx(old["integrated_coordinate_fieldwork"]+dw1+dw2)
    assert row["time"] == pytest.approx(old["time"]+.001)
    assert row["diagnostics"]["time"]==row["time"]
    assert row["steps"]==old["steps"]+2 and row["steps_this_batch"]==2
    assert row["initial_energy"]==old["initial_energy"]
    assert row["event_root"] is None and row["armed_after_positive_p_Q"]
    assert row["fieldwork_balance_residual"] == pytest.approx(
        row["diagnostics"]["energy"]["field"]-old["initial_energy"]["field"]-row["integrated_coordinate_fieldwork"])
    assert row["original_full_RHS"].endswith("nsc_discovery_leading_einstein")
    assert not confirm.full.positive_turn_bracket(0.,-1.,True)


def test_lock_or_artifact_tampering_blocks_advance(tmp_path,monkeypatch):
    for name in confirm.ORIGINAL_PINS:
        (tmp_path/name).write_bytes((confirm.DEFAULT_SOURCE/name).read_bytes())
    record=json.loads((tmp_path/"prediction.json").read_text())
    record["locked_before_held_data"]=False
    (tmp_path/"prediction.json").write_text(json.dumps(record))
    called=[]
    monkeypatch.setattr(confirm.full,"_advance",lambda *a,**k:called.append(True))
    with pytest.raises(ValueError,match="artifact changed"):
        confirm.run(tmp_path,tmp_path/"new",producing_commit=confirm.PARENT_COMMIT)
    assert not called


def test_default_preflight_is_read_only_and_caps_are_separate(source,monkeypatch):
    monkeypatch.setattr(confirm.full,"_advance",lambda *a,**k:pytest.fail("preflight advanced"))
    report=confirm.preview()
    assert not report["evolved"] and report["bytes_written"]==0
    assert report["fresh_CPU_limit_seconds"]==35 and report["original_CPU_limit_seconds"]==30
    assert not report["source_radius_eigenframe_reprepared"]
    with pytest.raises(ValueError,match="explicit frozen"):
        confirm.run()
    with pytest.raises(ValueError,match="35"):
        confirm.run(cpu_limit_seconds=36)
