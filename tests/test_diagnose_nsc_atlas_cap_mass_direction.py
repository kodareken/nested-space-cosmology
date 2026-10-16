"""Provenance and lifecycle controls; share the mathematical suite's allowance."""
import copy
import json
from pathlib import Path
import time

import pytest

import diagnose_nsc_atlas_cap_mass_direction as tool
from test_analyze_nsc_atlas_cap_mass_direction import budget
from test_diagnose_nsc_atlas_boundary_rows import auth_fixture


def owned_sources(root):
    for name in tool.OWN_SOURCES:
        path = root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"synthetic frozen producer")


def pinned_fixture(root, monkeypatch, budget):
    auth_fixture(root, monkeypatch)
    owned_sources(root)
    _, _, _, identity = tool.previous.authenticate(root, budget)
    bindings = tool._own_bindings(root)
    cap = dict(schema=tool.previous.SCHEMA, scope=tool.previous.SCOPE,
        limits=tool.previous.LIMITS, runtime=tool.runtime_versions(), input=identity,
        audit_source_bindings={n:bindings[n] for n in tool.previous.OWN_SOURCES},
        analysis={"saved_terminal_mass":.25})
    raw = tool.encoded(cap)
    path = root/tool.CAP_INPUT
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(raw)
    monkeypatch.setattr(tool, "CAP_SHA256", tool.digest(raw))
    return bindings, cap


def synthetic(monkeypatch, root, budget):
    owned_sources(root)
    current, cap, cache = {"cut":128}, {"analysis":{"mass":.25}}, {"fixed":True}
    identity = {"fixture":"pinned"}
    calls = []
    def auth(root, shared, bindings):
        assert shared is budget and set(bindings) == set(tool.OWN_SOURCES)
        calls.append("auth")
        return current, cap, cache, identity
    def analyze(c, a, k, shared):
        assert c is current and a is cap and k is cache and shared is budget
        calls.append("analysis")
        shared.check(7)
        return dict(cut=c["cut"], mass=a["analysis"]["mass"], evaluated_points=7,
            execution=dict(tool.EXPECTED_CALLS,
                charges=dict(input_and_profile=1, metrics=3, rhs=2, probe_analysis=1)))
    monkeypatch.setattr(tool, "authenticate", auth)
    monkeypatch.setattr(tool, "analyze", analyze)
    return calls


def test_pinned_cap_seven_sources_and_one_inherited_hash_pass_without_replay(tmp_path, monkeypatch, budget):
    bindings, cap = pinned_fixture(tmp_path, monkeypatch, budget)
    original = tool.previous.authenticate
    calls = []
    def auth(root, shared):
        calls.append(1)
        return original(root, shared)
    def forbidden(*args, **kwargs):
        raise AssertionError("previous scientific replay")
    monkeypatch.setattr(tool.previous, "authenticate", auth)
    for name in ("analyze", "check", "_validation"):
        monkeypatch.setattr(tool.previous, name, forbidden)
    for name in ("calculate", "load_inputs", "validate_record"):
        monkeypatch.setattr(tool.previous.previous.producer, name, forbidden)
    start = budget.evaluations
    _, got, _, identity = tool.authenticate(tmp_path, budget, bindings)
    assert got == cap and calls == [1] and budget.evaluations == start
    assert identity["inherited"]["authenticated_input_files"] == 104
    assert not identity["inherited"]["recursive_validation"]
    assert identity["cap_scientific_replays"] == 0
    assert len(tool.OWN_SOURCES) == len(set(tool.OWN_SOURCES)) == 12
    assert len(identity["cap_audit_source_bindings"]) == 7


@pytest.mark.parametrize("kind", ["cap", "prior_producer", "inherited_source", "inherited_manifest"])
def test_retained_input_or_producer_cannot_be_rebound(tmp_path, monkeypatch, budget, kind):
    bindings, cap = pinned_fixture(tmp_path, monkeypatch, budget)
    if kind == "cap":
        cap["analysis"]["saved_terminal_mass"] += 1
        (tmp_path/tool.CAP_INPUT).write_bytes(tool.encoded(cap))
    elif kind == "prior_producer":
        (tmp_path/tool.previous.OWN_SOURCES[-1]).write_bytes(b"changed prior control")
        bindings = tool._own_bindings(tmp_path)
    elif kind == "inherited_source":
        (tmp_path/"fixture/0.txt").write_bytes(b"changed inherited source")
    else:
        path = tmp_path/tool.previous.previous.INPUT
        payload = json.loads(path.read_bytes()); payload["extra"] = "replacement"
        raw = tool.encoded(payload); path.write_bytes(raw)
        manifest_path = path.with_name("manifest.json")
        manifest = json.loads(manifest_path.read_bytes())
        manifest["outputs"]["operator.json"] = dict(sha256=tool.digest(raw), bytes=len(raw))
        manifest_path.write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError):
        tool.authenticate(tmp_path, budget, bindings)


def test_creation_and_check_use_one_authentication_and_preserve_shared_budget(tmp_path, monkeypatch, budget):
    calls = synthetic(monkeypatch, tmp_path, budget)
    start, origin = budget.evaluations, budget.started
    path = tmp_path/"audit.json"
    saved = tool.create(path, root=tmp_path, budget=budget)
    assert calls == ["auth", "analysis", "analysis"]
    assert saved["execution"]["initial_evaluations"] == start
    assert saved["execution"]["evaluated_points"] == budget.evaluations-start == 14
    assert saved["execution"]["analysis_rhs_calls"] == 4
    assert saved["execution"]["hash_authentication_passes"] == 1
    assert set(saved["audit_source_bindings"]) == set(tool.OWN_SOURCES)
    checked = tool.check(path, root=tmp_path, budget=budget)
    assert checked["evaluated_points"] == 7 and checked["rhs_calls"] == 4
    assert calls == ["auth", "analysis", "analysis", "auth", "analysis"]
    assert budget.started == origin
    before = path.read_bytes()
    with pytest.raises(FileExistsError):
        tool.create(path, root=tmp_path, budget=budget)
    assert path.read_bytes() == before


@pytest.mark.parametrize("kind", ["analysis", "count", "boolean_count", "producer"])
def test_changed_saved_record_fails_check_once_with_receipt(tmp_path, monkeypatch, budget, kind):
    calls = synthetic(monkeypatch, tmp_path, budget)
    path = tmp_path/"audit.json"
    saved = tool.create(path, root=tmp_path, budget=budget)
    if kind == "analysis": saved["analysis"]["mass"] += 1
    elif kind == "count": saved["execution"]["evaluated_points"] += 1
    elif kind == "boolean_count": saved["execution"]["hash_authentication_passes"] = True
    else: (tmp_path/tool.OWN_SOURCES[-1]).write_bytes(b"changed current control")
    raw = tool.encoded(saved); path.write_bytes(raw)
    before = calls.count("analysis")
    with pytest.raises(ValueError) as caught:
        tool.check(path, root=tmp_path, budget=budget)
    receipt = json.loads(Path(caught.value.failure_receipt).read_bytes())
    assert receipt["phase"] == "check selected reconstruction" and receipt["retries"] == 0
    assert calls.count("analysis")-before == (0 if kind == "producer" else 1)
    assert path.read_bytes() == raw


@pytest.mark.parametrize("phase", ["authentication", "analysis", "validation", "producer_drift"])
def test_first_failure_retains_probe_context_and_no_completed_output(tmp_path, monkeypatch, budget, phase):
    calls = synthetic(monkeypatch, tmp_path, budget)
    old = tool.analyze
    if phase == "authentication":
        def fail(root, shared, bindings):
            calls.append("auth")
            raise ValueError("injected authentication failure")
        monkeypatch.setattr(tool, "authenticate", fail)
    else:
        def analyze(c, a, k, shared):
            result = old(c, a, k, shared)
            count = calls.count("analysis")
            if phase == "analysis" or (phase == "validation" and count == 2):
                error = ValueError("injected probe reconstruction")
                error.cap_mass_direction_failure = dict(phase="metric probe", amplitude="finite",
                    order=8, offset=-.001, counts=dict(metric_calls=17, rhs_calls=4))
                raise error
            if phase == "producer_drift" and count == 1:
                (tmp_path/tool.OWN_SOURCES[-1]).write_bytes(b"mid-analysis mutation")
            return result
        monkeypatch.setattr(tool, "analyze", analyze)
    path = tmp_path/"audit.json"; start = budget.evaluations
    with pytest.raises(ValueError) as caught:
        tool.create(path, root=tmp_path, budget=budget)
    receipt = json.loads(Path(caught.value.failure_receipt).read_bytes())
    expected = 0 if phase == "authentication" else 1 if phase == "analysis" else 2
    assert calls.count("auth") == 1 and calls.count("analysis") == expected
    assert receipt["evaluated_points"] == budget.evaluations-start == 7*expected
    assert receipt["retries"] == 0 and not path.exists()
    if phase in ("analysis", "validation"):
        assert receipt["analysis_failure"]["order"] == 8
        assert receipt["analysis_failure"]["counts"]["metric_calls"] == 17


@pytest.mark.parametrize("kind", ["points", "calls", "charges"])
def test_analyzer_inventory_must_match_actual_work(tmp_path, monkeypatch, budget, kind):
    synthetic(monkeypatch, tmp_path, budget)
    original = tool.analyze
    def changed(c, a, k, shared):
        row = copy.deepcopy(original(c, a, k, shared))
        if kind == "points": row["evaluated_points"] += 1
        elif kind == "calls": row["execution"]["rhs_calls"] = 0
        else: row["execution"]["charges"]["metrics"] += 1
        return row
    monkeypatch.setattr(tool, "analyze", changed)
    with pytest.raises(ValueError):
        tool.create(tmp_path/"audit.json", root=tmp_path, budget=budget)
    assert not (tmp_path/"audit.json").exists()


def test_allowance_guards_and_exhaustion_keep_parent_clock_and_failed_cost(tmp_path, monkeypatch, budget):
    for name, value in (("maximum_evaluations", 4_000_001), ("evaluations", -1),
            ("cpu_seconds", float("nan")), ("started", True),
            ("started", time.process_time()+100.)):
        child = tool.Budget(maximum_evaluations=100, cpu_seconds=30., started=budget.started)
        setattr(child, name, value)
        with pytest.raises(ValueError): tool.shared_budget(child)
    with pytest.raises(ValueError):
        tool.shared_budget(tool.Budget(maximum_evaluations=2_000_001, cpu_seconds=30.), checking=True)
    class Child(tool.Budget):
        def check(self, count=0):
            budget.check(count)
            super().check(count)
    child = Child(maximum_evaluations=6, cpu_seconds=30., started=budget.started)
    owned_sources(tmp_path)
    monkeypatch.setattr(tool, "authenticate", lambda root, b, bindings: ({}, {}, {}, {}))
    calls = []
    def fail(c, a, k, b):
        calls.append(1); b.check(7)
    monkeypatch.setattr(tool, "analyze", fail)
    start = budget.evaluations
    with pytest.raises(RuntimeError) as caught:
        tool.create(tmp_path/"audit.json", root=tmp_path, budget=child)
    assert calls == [1] and child.evaluations == budget.evaluations-start == 7
    assert child.started == budget.started
    assert json.loads(Path(caught.value.failure_receipt).read_bytes())["retries"] == 0
