"""Recorder/provenance controls sharing the analyzer's one admitted allowance."""
import copy
import json
import time

import pytest

import diagnose_nsc_atlas_cap_consistency as tool
from test_analyze_nsc_atlas_cap_consistency import budget
from test_diagnose_nsc_atlas_boundary_rows import auth_fixture


def owned_sources(root):
    for name in tool.OWN_SOURCES:
        path=root/name
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(b"synthetic frozen producer")


def synthetic(monkeypatch,tmp_path,budget):
    owned_sources(tmp_path)
    current,cache,identity={"selected":63},{"fixed":True},{"fixture":"pinned"}
    calls=[]
    def auth(root,shared):
        assert shared is budget
        calls.append("auth")
        return current,{},cache,identity
    def analyze(c,k,shared):
        assert c is current and k is cache and shared is budget
        calls.append("analysis")
        shared.check(7)
        return dict(selected=c["selected"],evaluated_points=7,execution=dict(tool.EXPECTED_CALLS,
            charges=dict(input_and_profile=1,metrics=3,selected_prefix_maps=2,derivative_analysis=1)))
    monkeypatch.setattr(tool,"authenticate",auth)
    monkeypatch.setattr(tool,"analyze",analyze)
    return calls


def test_pinned_inputs_and_inherited_producers_without_scientific_replay(tmp_path,monkeypatch,budget):
    auth_fixture(tmp_path,monkeypatch)
    def forbidden(*args,**kwargs):
        raise AssertionError("old scientific replay")
    for name in ("calculate","load_inputs","validate_record"):
        monkeypatch.setattr(tool.previous.producer,name,forbidden)
    before=budget.evaluations
    _,_,_,identity=tool.authenticate(tmp_path,budget)
    assert identity["authenticated_input_files"]==104 and not identity["recursive_validation"]
    assert budget.evaluations==before
    (tmp_path/"fixture/0.txt").write_bytes(b"changed")
    with pytest.raises(ValueError,match="binding changed"):
        tool.authenticate(tmp_path,budget)


def test_rehashed_manifest_cannot_replace_the_pinned_input(tmp_path,monkeypatch,budget):
    auth_fixture(tmp_path,monkeypatch)
    path=tmp_path/tool.previous.INPUT
    payload=json.loads(path.read_bytes()); payload["extra"]="tampered"
    raw=tool.encoded(payload);path.write_bytes(raw)
    manifest_path=path.with_name("manifest.json")
    manifest=json.loads(manifest_path.read_bytes())
    manifest["outputs"]["operator.json"]=dict(sha256=tool.digest(raw),bytes=len(raw))
    manifest_path.write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError,match="pinned"):
        tool.authenticate(tmp_path,budget)


def test_creation_and_one_check_preserve_clock_counts_and_creation_only(tmp_path,monkeypatch,budget):
    calls=synthetic(monkeypatch,tmp_path,budget)
    start,origin=budget.evaluations,budget.started
    path=tmp_path/"audit.json"
    created=tool.create(path,root=tmp_path,budget=budget)
    assert calls==["auth","analysis","analysis"]
    assert created["execution"]["initial_evaluations"]==start
    assert created["execution"]["analysis_evaluations"]==7
    assert created["execution"]["creation_validation_evaluations"]==7
    assert created["execution"]["evaluated_points"]==budget.evaluations-start==14
    assert created["execution"]["hash_authentication_passes"]==1
    assert set(created["audit_source_bindings"])==set(tool.OWN_SOURCES)
    assert tool.check(path,root=tmp_path,budget=budget)["evaluated_points"]==7
    assert calls==["auth","analysis","analysis","auth","analysis"] and budget.started==origin
    with pytest.raises(FileExistsError):
        tool.create(path,root=tmp_path,budget=budget)


@pytest.mark.parametrize("kind",["analysis","count","boolean_count","producer"])
def test_changed_saved_result_and_bound_producer_rejected(tmp_path,monkeypatch,budget,kind):
    calls=synthetic(monkeypatch,tmp_path,budget)
    path=tmp_path/"audit.json"; saved=tool.create(path,root=tmp_path,budget=budget)
    if kind=="analysis": saved["analysis"]["selected"]=68
    elif kind=="count": saved["execution"]["evaluated_points"]+=1
    elif kind=="boolean_count": saved["execution"]["hash_authentication_passes"]=True
    else: (tmp_path/tool.OWN_SOURCES[-1]).write_bytes(b"changed test")
    path.write_bytes(tool.encoded(saved))
    before=calls.count("analysis")
    with pytest.raises(ValueError) as caught:
        tool.check(path,root=tmp_path,budget=budget)
    receipt=json.loads(open(caught.value.failure_receipt).read())
    assert receipt["phase"]=="check selected reconstruction" and receipt["retries"]==0
    assert calls.count("analysis")-before==(0 if kind=="producer" else 1)
    assert path.read_bytes()==tool.encoded(saved)


@pytest.mark.parametrize("phase",["authentication","analysis","validation","producer_drift"])
def test_first_failure_receipt_without_retry_or_completed_output(tmp_path,monkeypatch,budget,phase):
    calls=synthetic(monkeypatch,tmp_path,budget)
    old_analyze=tool.analyze
    if phase=="authentication":
        def fail(root,shared):
            calls.append("auth")
            raise ValueError("injected authentication")
        monkeypatch.setattr(tool,"authenticate",fail)
    else:
        def analyze(c,k,shared):
            result=old_analyze(c,k,shared)
            count=calls.count("analysis")
            if phase=="analysis" or (phase=="validation" and count==2):
                error=ValueError("injected reconstruction")
                error.cap_consistency_failure=dict(phase=dict(name="cap metric",amplitude="finite",order=8,cut=68),
                    counts=dict(metric_calls=41,prefix_calls=6,hybrid_calls=6))
                raise error
            if phase=="producer_drift" and count==1:
                (tmp_path/tool.OWN_SOURCES[-1]).write_bytes(b"mid-analysis change")
            return result
        monkeypatch.setattr(tool,"analyze",analyze)
    path=tmp_path/"audit.json"; start=budget.evaluations
    with pytest.raises(ValueError) as caught:
        tool.create(path,root=tmp_path,budget=budget)
    receipt=json.loads(open(caught.value.failure_receipt).read())
    expected=0 if phase=="authentication" else 1 if phase=="analysis" else 2
    assert calls.count("auth")==1 and calls.count("analysis")==expected
    assert receipt["evaluated_points"]==budget.evaluations-start==7*expected
    assert receipt["retries"]==0 and not path.exists()
    if phase in ("analysis","validation"):
        assert receipt["analysis_failure"]["phase"]["cut"]==68
        assert receipt["analysis_failure"]["counts"]["metric_calls"]==41


@pytest.mark.parametrize("name,value",[("maximum_evaluations",8_000_001),("maximum_evaluations",True),
    ("evaluations",-1),("cpu_seconds",float("nan")),("cpu_seconds",121.),
    ("started",True),("started",-1.),("started",float("inf")),("started","future")])
def test_invalid_allowance_rejected_without_work(budget,name,value):
    child=tool.Budget(maximum_evaluations=100,cpu_seconds=30.,started=budget.started)
    setattr(child,name,time.process_time()+100. if value=="future" else value)
    with pytest.raises(ValueError):
        tool.shared_budget(child)


def test_standalone_ceiling_and_deliberate_exhaustion_debit_parent(tmp_path,monkeypatch,budget):
    with pytest.raises(ValueError):
        tool.shared_budget(tool.Budget(maximum_evaluations=4_000_001,cpu_seconds=60.),checking=True)
    class Child(tool.Budget):
        def check(self,count=0):
            budget.check(count)
            super().check(count)
    child=Child(maximum_evaluations=6,cpu_seconds=30.,started=budget.started)
    owned_sources(tmp_path)
    monkeypatch.setattr(tool,"authenticate",lambda root,b: ({},{},{},{}))
    calls=[]
    def analysis(c,k,b):
        calls.append(1);b.check(7)
    monkeypatch.setattr(tool,"analyze",analysis)
    start=budget.evaluations
    with pytest.raises(RuntimeError) as caught:
        tool.create(tmp_path/"audit.json",root=tmp_path,budget=child)
    assert calls==[1] and child.evaluations==budget.evaluations-start==7
    assert child.started==budget.started
    assert json.loads(open(caught.value.failure_receipt).read())["retries"]==0


@pytest.mark.parametrize("kind",["points","calls"])
def test_analyzer_receipt_must_match_actual_work(tmp_path,monkeypatch,budget,kind):
    synthetic(monkeypatch,tmp_path,budget)
    old=tool.analyze
    def changed(c,k,b):
        result=copy.deepcopy(old(c,k,b))
        if kind=="points": result["evaluated_points"]=8
        else: result["execution"]["metric_calls"]=63
        return result
    monkeypatch.setattr(tool,"analyze",changed)
    with pytest.raises(ValueError):
        tool.create(tmp_path/"audit.json",root=tmp_path,budget=budget)
    assert not (tmp_path/"audit.json").exists()
