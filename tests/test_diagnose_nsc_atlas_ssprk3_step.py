"""Manufactured six-stage ledgers and pinned-input controls; no physical audit."""
import copy
import json
import time

import numpy as np
import pytest

import diagnose_nsc_atlas_ssprk3_step as tool
from recursive_horizons import nsc_atlas_background_ssprk3 as operator
from recursive_horizons.nsc_atlas_compatible_evolution import quiet_departures
from test_diagnose_nsc_atlas_background_step import cache, profile, manufactured_stage


@pytest.fixture(scope="module")
def budget():
    value = tool.Budget(maximum_evaluations=100000, cpu_seconds=30.)
    yield value
    value.check()
    print("SSPRK3 ledger manufactured controls", value.evaluations, time.process_time()-value.started)


def run_step(monkeypatch, budget, alpha=.3, simple=False, rows=None):
    calls = []
    monkeypatch.setattr(operator, "background_stage", manufactured_stage(calls, simple=simple))
    rows = cache([0., .02, .04], simple) if rows is None else rows
    initial = quiet_departures(exterior_cells=4, cap_cells=3)
    if not simple:
        initial["cap"]["Iq"][:] = .01
    counter = dict(stage_calls=0, max_courant=0.)
    row = tool.weighted_step(initial, 0, 1, 2, source_rows=dict(x=rows["x"]),
        gauge_rows=rows, profile=profile(), alpha=alpha, budget=budget, counter=counter, label="manufactured")
    return row, calls, initial, counter


@pytest.mark.parametrize("alpha", (0., .3))
def test_three_stage_weighted_chain_and_actual_reference_geometry(monkeypatch, budget, alpha):
    row, calls, initial, counter = run_step(monkeypatch, budget, alpha)
    assert len(calls) == counter["stage_calls"] == 3
    assert [c[0] for c in calls] == [0., .04, .02]
    assert all(c[2] is budget for c in calls)
    assert all(abs(v) < 8e-16 for v in row["closures"].values())
    assert row["increment"] == pytest.approx(sum(row["components"].values()), abs=8e-16)
    for key, coefficient in (("r", .02), ("T", .05)):
        geometry = row["reference_geometry"][key]
        assert abs(geometry["reference_defect"]) < 2e-16  # Simpson exact for quadratic rates.
        assert abs(geometry["closure"]) < 2e-16
        assert row["stage_geometry_defects"][key] == pytest.approx(
            [coefficient*.04**3, -5*coefficient*.04**3/8], abs=2e-16)
    np.testing.assert_array_equal(initial["cap"]["Iq"], np.full(3, .01))
    if alpha == 0:
        assert abs(row["components"]["frame_remainder"]) < 3e-16
        assert row["components"]["background_remainder"] == pytest.approx(row["background_chain"], abs=2e-16)
    else:
        assert abs(row["components"]["frame_remainder"]) > 1e-14


def test_matching_cubic_background_has_exact_weighted_flux_without_wall_reset(monkeypatch, budget):
    row, calls, _, _ = run_step(monkeypatch, budget, alpha=0., simple=True)
    assert len(calls) == 3
    assert abs(row["background_chain"]) < 4e-16
    assert abs(row["accepted"]["departure_chi"]) < 1e-16
    assert row["increment"] == pytest.approx(row["components"]["exterior_flux"], abs=4e-16)
    assert row["accepted"]["wall"]["c"] == row["accepted"]["background_chi"]


def test_exact_caller_midpoint_key_is_preserved(monkeypatch, budget):
    mid = np.nextafter(.02, 1.)
    row, calls, _, _ = run_step(monkeypatch, budget, rows=cache([0., mid, .04]))
    assert calls[2][0] == mid != .04/2
    assert row["x_mid"] == mid


def fixture(tmp_path, monkeypatch, budget):
    rows, prof = cache(np.linspace(0., 1., 513)), profile()
    cases = [{} for _ in range(6)]
    base = dict(gauge=dict(rows=rows), initial_data={})
    old_cases = []
    for case_id, label, alpha in ((1, "zero", 0.), (4, "finite", .3)):
        state = quiet_departures(exterior_cells=32, cap_cells=32)
        state["cap"]["Iq"][:] = .001
        initial = tool.compact_state(state, tool._cache(rows, 328), prof, alpha, budget)
        cases[case_id] = dict(cells=32, steps=256, amplitude=label, alpha=alpha,
            snapshots=[None]*82+[dict(x=rows["x"][328], state=state, physical=dict(wall=initial["wall"]))])
        base["initial_data"][label] = dict(snapshot_rows=dict(x=rows["x"]))
        # Saved reference values are algebraic fixtures; no former RK2 stage is invoked.
        coarse, fine = copy.deepcopy(initial), copy.deepcopy(initial)
        coarse["wall"]["c"] += .001
        fine["wall"]["c"] += .0012
        coarse_parts = {k:0. for k in tool.COMPONENTS}; coarse_parts["exterior_flux"] = .001
        fine_parts = {k:0. for k in tool.COMPONENTS}; fine_parts["exterior_flux"] = .0012
        old_cases.append(dict(case=case_id, amplitude=label, alpha=alpha,
            coarse=dict(initial=initial, accepted=coarse, components=coarse_parts),
            halves=[{}, dict(accepted=fine)], fine_components=fine_parts))
    record = dict(cases=cases, profile_inputs={}, profile={})
    identity = dict(pinned=True, authenticated_input_files=66)
    monkeypatch.setattr(tool.previous.previous, "authenticate", lambda root:(record, base, identity))
    def fit(inputs, shared):
        shared.check(3)
        return prof, {}
    monkeypatch.setattr(tool.producer, "fit_profile", fit)
    for name in tool.SOURCES:
        if name == tool.STEP_INPUT:
            continue
        path = tmp_path/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("synthetic bound source")
    prior = dict(schema=tool.previous.SCHEMA, input=identity, limits=tool.previous.LIMITS,
        runtime=tool.producer.runtime_versions(), settings=dict(cases=[1, 4], cut=82, cache_indices=[328, 329, 330]),
        audit_source_bindings={n:tool.digest((tmp_path/n).read_bytes()) for n in tool.previous.SOURCES}, cases=old_cases)
    path = tmp_path/tool.STEP_INPUT
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(tool.encoded(prior))
    monkeypatch.setattr(tool, "STEP_SHA256", tool.digest(path.read_bytes()))
    calls = []
    monkeypatch.setattr(operator, "background_stage", manufactured_stage(calls))
    def forbidden(*args, **kwargs):
        raise AssertionError("prior-step or full physical replay forbidden")
    monkeypatch.setattr(tool.previous, "audit", forbidden)
    monkeypatch.setattr(tool.previous, "accepted_step", forbidden)
    for name in ("evolve_background", "extract_profile_inputs", "load_parent", "validate_record"):
        monkeypatch.setattr(tool.producer, name, forbidden)
    return calls, prior, record, base


def test_six_exact_stages_shared_budget_and_algebraic_retained_comparisons(tmp_path, monkeypatch, budget):
    calls, prior, _, base = fixture(tmp_path, monkeypatch, budget)
    before = budget.evaluations
    row = tool.audit(tmp_path, budget=budget)
    assert len(calls) == row["execution"]["stage_calls"] == 6
    assert [c[0] for c in calls] == [base["gauge"]["rows"]["x"][i] for i in [328, 330, 329]*2]
    execution = row["execution"]
    assert execution["old_ssprk2_steps_replayed"] == execution["source_integrations"] == 0
    assert execution["full_trajectory_replays"] == execution["accepted_endpoint_extra_stage_calls"] == 0
    assert execution["evaluated_points"] == before+execution["profile_fit_points"]+execution["local_step_points"]+execution["comparison_points"]
    for case, old in zip(row["cases"], prior["cases"]):
        comparison = case["comparison"]
        c, f = old["coarse"]["accepted"]["wall"]["c"], old["halves"][1]["accepted"]["wall"]["c"]
        assert comparison["endpoints"]["rk2_richardson"]["c"] == (4*f-c)/3
        assert "not truth" in comparison["richardson_scope"]
        assert max(abs(v) for v in comparison["component_closures"].values()) < 6e-16
    assert max(abs(v["closure"]) for v in row["paired"].values()) < 8e-16


@pytest.mark.parametrize("mutation", ("payload", "rehashed_source", "source"))
def test_pinned_old_step_and_producer_drift_fail_before_any_stage(tmp_path, monkeypatch, budget, mutation):
    calls, _, _, _ = fixture(tmp_path, monkeypatch, budget)
    path = tmp_path/tool.STEP_INPUT
    if mutation == "payload":
        raw = json.loads(path.read_bytes()); raw["cases"][0]["coarse"]["accepted"]["wall"]["c"] += .1
        path.write_bytes(tool.encoded(raw))
    else:
        name = tool.previous.SOURCES[0]
        (tmp_path/name).write_text("changed")
        if mutation == "rehashed_source":
            raw = json.loads(path.read_bytes()); raw["audit_source_bindings"][name] = tool.digest(b"changed")
            path.write_bytes(tool.encoded(raw))
    with pytest.raises((ValueError, RuntimeError), match="payload mismatch|binding mismatch"):
        tool.audit(tmp_path, budget=budget)
    assert calls == []


def test_failure_retains_actual_attempt_count_without_retry(tmp_path, monkeypatch, budget):
    calls, _, _, _ = fixture(tmp_path, monkeypatch, budget)
    def failing(*args, **kwargs):
        calls.append("attempt"); kwargs["budget"].check(1)
        raise ValueError("manufactured stage failure")
    monkeypatch.setattr(operator, "background_stage", failing)
    with pytest.raises(RuntimeError, match=r"case=1 amplitude=zero.*stage_calls=1.*manufactured stage failure"):
        tool.audit(tmp_path, budget=budget)
    assert calls == ["attempt"]


def test_post_run_source_change_is_rejected(tmp_path, monkeypatch, budget):
    calls, _, _, _ = fixture(tmp_path, monkeypatch, budget)
    original = operator.background_stage
    def changed(*args, **kwargs):
        result = original(*args, **kwargs)
        if len(calls) == 6:
            (tmp_path/tool.SOURCES[-1]).write_text("changed after stages")
        return result
    monkeypatch.setattr(operator, "background_stage", changed)
    with pytest.raises(RuntimeError, match="producer changed"):
        tool.audit(tmp_path, budget=budget)
    assert len(calls) == 6


@pytest.mark.parametrize("key,value", (("maximum_evaluations", 150001), ("maximum_evaluations", True),
    ("cpu_seconds", 31.), ("cpu_seconds", float("nan")), ("evaluations", -1), ("started", float("nan"))))
def test_shared_allowance_rejects_expansion_and_malformed_values(key, value):
    value_budget = tool.Budget(maximum_evaluations=100000, cpu_seconds=30.)
    setattr(value_budget, key, value)
    with pytest.raises(ValueError, match="allowance"):
        tool._shared(value_budget)


def test_creation_only_and_check_replays_only_new_diagnostic(tmp_path, monkeypatch):
    row = dict(schema=tool.SCHEMA, execution=dict(cpu_seconds=.01, evaluated_points=123), paired=.12)
    calls = []
    def fake(root):
        calls.append(root)
        return copy.deepcopy(row)
    monkeypatch.setattr(tool, "audit", fake)
    path = tmp_path/"ssprk3.json"
    tool.create(path, root=tmp_path)
    row["execution"]["cpu_seconds"] = .02
    assert tool.check(path, root=tmp_path)["stage_calls"] == 6
    with pytest.raises(FileExistsError):
        tool.create(path, root=tmp_path)
    assert len(calls) == 2
    saved = json.loads(path.read_bytes()); saved["paired"] += .01
    path.write_bytes(tool.encoded(saved))
    with pytest.raises(ValueError, match="recomputation mismatch"):
        tool.check(path, root=tmp_path)
