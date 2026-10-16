"""Bounded manufactured record authentication; no physical trajectory run."""
import copy
import json
import time

import numpy as np
import pytest

import record_nsc_atlas_ssprk3_evolution as tool
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_compatible_evolution import quiet_departures

ORIGINAL_DIAGNOSTICS = tool.diagnostics
ORIGINAL_LOAD = tool.load_inputs


@pytest.fixture(scope="module")
def budget():
    value = tool.Budget(maximum_evaluations=100000, cpu_seconds=30.)
    yield value
    value.check()
    print("SSPRK3 recorder manufactured controls", value.evaluations, time.process_time()-value.started)


@pytest.fixture(scope="module")
def synthetic(budget):
    """One small record, with real field reconstruction and no PDE calls."""
    patch = pytest.MonkeyPatch()
    cases = tuple((label, 3, steps) for label in ("zero", "finite") for steps in (128, 256))
    patch.setattr(tool, "CASES", cases)
    x = np.linspace(.4, .41, 513)
    base = dict(geometry=dict(synthetic=True), gauge=dict(rows=dict(x=x,
        W=np.full(513, .2), W_x=np.ones(513), W_xx=np.zeros(513),
        r=np.full(513, .45), L=np.zeros(513), T=x-.4, s=x-.4)))
    base["initial_data"] = {label:dict(alpha=alpha, snapshot_rows=dict(x=x, r=x,
        log_C=np.zeros(513), c=np.zeros(513), r_x=np.ones(513), log_C_x=np.zeros(513),
        c_x=np.zeros(513))) for label, alpha in (("zero", 0.), ("finite", .01))}
    parent = dict(profile_inputs=dict(retained_table="fixture"), profile=dict(synthetic=True))
    identity = dict(pinned_parent="fixture", pinned_local_receipt="fixture", original_source_bindings=66)
    loads, stages, evolutions = [], [], []
    def load(shared):
        loads.append(shared)
        shared.check()
        return parent, base, identity
    patch.setattr(tool, "load_inputs", load)
    profile = CharacteristicProfile([0., 1.5], [0., 0.], [0., 0.], [0., 0.])
    def fit(inputs, shared):
        shared.check(3)
        return profile, dict(synthetic=True)
    patch.setattr(tool, "fit_profile", fit)
    def stage(state, *, eta, budget, **kwargs):
        stages.append(budget)
        budget.check(8*len(eta)+5)
        return dict(state=state)
    patch.setattr(tool, "background_stage", stage)
    patch.setattr(tool, "final_summary", lambda row:dict(cap_metric_evaluations=1,
        cap_rhs_evaluations=2, checked_cap=tool.plain(row["state"]["cap"])))
    records = []
    saved_end = 3
    for label, n, steps in cases:
        alpha = base["initial_data"][label]["alpha"]
        snapshots = []
        for at in x[::4]:
            state = quiet_departures(exterior_cells=n, cap_cells=n)
            state["cap"]["Iq"][:] = alpha*(at-.4)**2
            state["cap"]["Ip"][:] = 2*alpha*(at-.4)**2
            index = int(np.where(x == at)[0][0])
            physical = tool.physical_view(state, float(at),
                tool.subset(base["initial_data"][label]["snapshot_rows"], index),
                tool.subset(base["gauge"]["rows"], index), np.linspace(0, 1, n+1), alpha, profile, budget)
            snapshots.append(dict(x=float(at), state=tool.plain(state), physical=physical))
        # Saved counters are manufactured data. No fictitious trajectory work
        # is charged to the live test allowance or claimed to have run.
        count = 3*steps+1
        charge = count*(8*(n+1)+5)
        records.append(dict(amplitude=label,cells=n,steps=steps,alpha=alpha,snapshots=snapshots,
            final=tool.final_summary(dict(state=state)),execution=dict(budget_start=saved_end,
                budget_end=saved_end+charge,measured_evaluations=charge,stages=count,
                source_cache_reads=count,gauge_cache_reads=count,max_courant=.1,retries=0,
                final_stage_evaluations=8*(n+1)+5)))
        saved_end += charge
    def forbidden_evolve(*args, **kwargs):
        evolutions.append(True)
        raise AssertionError("authentication must not propagate a trajectory")
    patch.setattr(tool, "evolve_background_ssprk3", forbidden_evolve)
    def diagnostics(record, parent, base, profile, shared):
        shared.check(7)
        return dict(physical_cap=[case["snapshots"][-1]["physical"]["cap"] for case in record["cases"]])
    patch.setattr(tool, "diagnostics", diagnostics)
    record = dict(schema=tool.SCHEMA,status=tool.STATUS,scope=tool.SCOPE,settings=tool.SETTINGS,
        runtime=tool.runtime_versions(),parent_input=identity,geometry=base["geometry"],
        profile_inputs=parent["profile_inputs"],profile=dict(synthetic=True),cases=records,
        refinement_is_rigorous_error_enclosure=False,
        execution=dict(initial_evaluations=0,evaluations=saved_end+7,cpu_seconds=1.,
            inherited_authentication_evaluations=0,profile_build_evaluations=3,diagnostic_evaluations=7,
            old_profile_stage_reconstructions=0,coupled_evolutions=4,time_steps=sum(k for _,_,k in cases),
            stage_evaluations=sum(3*k+1 for _,_,k in cases),source_integrations=0,gauge_quadratures=0,
            new_physical_preparations=0,stage_point_totals_are_measured=True))
    record["diagnostics"] = diagnostics(record, parent, base, profile, budget)
    yield dict(record=record, base=base, parent=parent, identity=identity,
               loads=loads, stages=stages, evolutions=evolutions)
    patch.undo()


def test_saved_physical_reconstruction_and_four_final_stages_share_caller_budget(synthetic, budget):
    before = budget.evaluations
    stage_start = len(synthetic["stages"])
    evolution_start = len(synthetic["evolutions"])
    auth = tool.validate_record(copy.deepcopy(synthetic["record"]), budget, require_creation=False)
    assert budget.evaluations-before == auth["authentication_evaluations"]
    assert len(synthetic["stages"])-stage_start == 4
    assert len(synthetic["evolutions"]) == evolution_start
    assert all(b is budget for b in synthetic["loads"]+synthetic["stages"])
    assert auth["source_integrations"] == auth["time_evolutions"] == 0
    assert synthetic["record"]["execution"]["source_integrations"] == 0
    assert synthetic["record"]["execution"]["stage_evaluations"] == sum(3*n+1 for _, _, n in tool.CASES)


def test_calculate_uses_each_exact_cache_stride_and_midpoint(synthetic, budget, monkeypatch):
    calls = []
    def evolve(initial, nodes, *, mid_nodes, source, gauge, budget, **kwargs):
        index = len(calls)
        saved = synthetic["record"]["cases"][index]
        calls.append((np.array(nodes), np.array(mid_nodes), budget))
        count = 3*saved["steps"]+1
        source.calls = gauge.calls = count
        budget.check(1)  # fixture orchestration only; no stage is executed
        return dict(snapshots=saved["snapshots"],final_stage=dict(state=saved["snapshots"][-1]["state"]),
            stages=count,steps=saved["steps"],max_courant=.1,retries=0,time_scheme=tool.SETTINGS["time_scheme"])
    monkeypatch.setattr(tool, "evolve_background_ssprk3", evolve)
    result = tool.calculate(budget)
    assert result["execution"]["time_steps"] == 768
    assert result["execution"]["stage_evaluations"] == 2308
    for (nodes, mids, shared), (_, _, steps) in zip(calls, tool.CASES):
        stride = 512//steps
        assert shared is budget
        np.testing.assert_array_equal(nodes, synthetic["base"]["gauge"]["rows"]["x"][::stride])
        np.testing.assert_array_equal(mids, synthetic["base"]["gauge"]["rows"]["x"][stride//2::stride])


def test_comparisons_pair_matching_amplitudes_and_use_total_caps(synthetic, budget, monkeypatch):
    record = copy.deepcopy(synthetic["record"])
    base = copy.deepcopy(synthetic["base"])
    base["geometry"]["support_x"] = [.4, .405]
    for case, value in zip(record["cases"], (1., 2., 4., 8.)):
        case["marker"] = value
    old = [dict(amplitude=a,cells=32,steps=k,marker=v) for a,k,v in
           (("zero",256,10.),("zero",512,20.),("finite",256,40.),("finite",512,80.))]
    def tables(case, rows, shared):
        shared.check(1)
        value = np.full(129, case["marker"])
        return dict(exterior=dict(c=value),wall=dict(c=value),cap=dict(Ip=value,Iq=value))
    monkeypatch.setattr(tool, "tables", tables)
    monkeypatch.setattr(tool, "_histories", lambda *args: ({}, {}))
    monkeypatch.setattr(tool, "compact_history", lambda row: row)
    monkeypatch.setattr(tool, "cap_history", lambda *args: dict(arrays=dict(upper_label_margin=[1.])))
    oracle_rows = []
    def oracle(case, shared):
        oracle_rows.append(case)
        for row in case["snapshots"]:
            assert row["state"]["cap"] == row["physical"]["cap"]
        return dict(checked=True)
    monkeypatch.setattr(tool, "cap_pullback", oracle)
    result = ORIGINAL_DIAGNOSTICS(record, dict(cases=old), base, None, budget)
    assert len(oracle_rows) == 2
    np.testing.assert_array_equal(result["paired_time"]["wall_arrays"]["c"], np.full(129, 3.))
    np.testing.assert_array_equal(result["paired_versus_ssprk2"][0]["256"]["wall_arrays"]["c"], np.full(129, -27.))
    np.testing.assert_array_equal(result["versus_ssprk2"][3]["512"]["wall_arrays"]["c"], np.full(129, -72.))
    assert result["temporal_difference_is_measured_order"] is False
    assert result["constructed_identities_are_accuracy_tests"] is False


def test_local_receipt_and_bound_source_are_hash_authenticated_without_replay(tmp_path, monkeypatch, budget):
    monkeypatch.setattr(tool, "ROOT", tmp_path)
    monkeypatch.setattr(tool, "LOCAL_INPUT", "local.json")
    monkeypatch.setattr(tool.local, "SOURCES", ("source.py",))
    identity = dict(original=dict(pinned=True),local_step=dict(pinned=True))
    calls = []
    def inherited(root):
        calls.append(root)
        return {"parent":True}, {"base":True}, {}, identity
    monkeypatch.setattr(tool.local, "authenticate", inherited)
    (tmp_path/"source.py").write_bytes(b"frozen")
    saved = dict(schema=tool.local.SCHEMA,input=identity,limits=tool.local.LIMITS,
        runtime=tool.runtime_versions(),audit_source_bindings={"source.py":tool.digest(b"frozen")},
        execution=dict(stage_calls=6))
    raw = tool.encoded(saved)
    monkeypatch.setattr(tool, "LOCAL_SHA256", tool.digest(raw))
    (tmp_path/"local.json").write_bytes(raw)
    before = budget.evaluations
    parent, base, receipt = ORIGINAL_LOAD(budget)
    assert parent == {"parent":True} and base == {"base":True}
    assert calls == [tmp_path] and budget.evaluations == before
    assert receipt["ssprk3_local"]["local_steps_replayed"] == 0
    (tmp_path/"source.py").write_bytes(b"changed")
    with pytest.raises(ValueError, match="source binding"):
        ORIGINAL_LOAD(budget)
    (tmp_path/"source.py").write_bytes(b"frozen")
    (tmp_path/"local.json").write_bytes(raw+b" ")
    with pytest.raises(ValueError, match="pinned local"):
        ORIGINAL_LOAD(budget)


@pytest.mark.parametrize("mutation", ("cap_total", "raw_cap", "wall", "exterior", "shape", "time", "case", "profile", "runtime"))
def test_authentication_rejects_state_time_and_identity_changes(synthetic, budget, mutation):
    record = copy.deepcopy(synthetic["record"])
    sample = record["cases"][0]["snapshots"][1]
    if mutation == "cap_total":
        sample["physical"]["cap"]["Ip"][0] += .01
    elif mutation == "raw_cap":
        sample["state"]["cap"]["Iq"][0] += .01
    elif mutation == "wall":
        sample["physical"]["wall"]["c"] += .01
    elif mutation == "exterior":
        sample["physical"]["increment"]["r"][0] += .01
    elif mutation == "shape":
        sample["state"]["departure"]["c"].append(0.)
    elif mutation == "time":
        sample["x"] = float(np.nextafter(sample["x"], 1.))
    elif mutation == "case":
        record["cases"][0]["steps"] *= 2
    elif mutation == "profile":
        record["profile_inputs"]["retained_table"] = "changed"
    else:
        record["runtime"]["numpy"] = ""
    with pytest.raises(ValueError):
        tool.validate_record(record, budget, require_creation=False)


@pytest.mark.parametrize("mutation", ("final", "stage_count", "case_points", "diagnostic", "total"))
def test_authentication_rejects_recomputed_diagnostic_and_cost_changes(synthetic, budget, mutation):
    record = copy.deepcopy(synthetic["record"])
    if mutation == "final":
        record["cases"][0]["final"]["cap_metric_evaluations"] += 1
    elif mutation == "stage_count":
        record["cases"][0]["execution"]["stages"] += 1
    elif mutation == "case_points":
        record["cases"][0]["execution"]["measured_evaluations"] += 1
    elif mutation == "diagnostic":
        record["diagnostics"]["physical_cap"][0]["Ip"][0] += .01
    else:
        record["execution"]["evaluations"] += 1
    with pytest.raises(ValueError):
        tool.validate_record(record, budget, require_creation=False)


def test_create_is_atomic_on_source_drift_and_keeps_failure_receipt(tmp_path, monkeypatch):
    monkeypatch.setattr(tool, "ROOT", tmp_path)
    state = dict(changed=False)
    monkeypatch.setattr(tool, "source_bindings", lambda:dict(source="after" if state["changed"] else "before"))
    monkeypatch.setattr(tool, "verify_saved_bindings", lambda expected, root, **kw: expected == tool.source_bindings())
    monkeypatch.setattr(tool, "calculate", lambda budget:dict(runtime=tool.runtime_versions(), execution=dict(evaluations=0)))
    def validate(*args, **kwargs):
        state["changed"] = True
        return dict(authentication_evaluations=0)
    monkeypatch.setattr(tool, "validate_record", validate)
    output = tmp_path/"rejected"
    with pytest.raises(RuntimeError, match="authentication") as error:
        tool.create(output)
    assert not output.exists()
    receipt = json.loads(open(error.value.failure_receipt).read())
    assert receipt["status"] == "FAILED_NOT_COMPLETED"


def test_creation_only_hash_guard_and_standalone_authentication_ceiling(tmp_path, monkeypatch):
    monkeypatch.setattr(tool, "ROOT", tmp_path)
    created_budgets, binding_budget_ids = [], []
    original_budget = tool.Budget
    def make_budget(**kwargs):
        shared = original_budget(**kwargs)
        created_budgets.append(shared)
        return shared
    monkeypatch.setattr(tool, "Budget", make_budget)
    def bindings():
        assert created_budgets  # CPU clock must begin before hash reads.
        binding_budget_ids.append(id(created_budgets[-1]))
        assert created_budgets[-1].started <= time.process_time()
        return dict(source="fixed")
    monkeypatch.setattr(tool, "source_bindings", bindings)
    monkeypatch.setattr(tool, "verify_saved_bindings", lambda expected, root, **kw: expected == tool.source_bindings())
    runtime, calls = tool.runtime_versions(), []
    def calculate(shared):
        calls.append(("calculate", shared))
        shared.check(70)
        return dict(runtime=runtime, execution=dict(evaluations=70, cpu_seconds=0.))
    def validate(record, shared=None, **kwargs):
        calls.append(("validate", shared))
        if shared is not None:
            shared.check(20)
        return dict(authentication_evaluations=20)
    monkeypatch.setattr(tool, "calculate", calculate)
    monkeypatch.setattr(tool, "validate_record", validate)
    output = tmp_path/"new"
    result = tool.create(output)
    assert result["creation_execution"]["total_evaluations"] == 90
    assert calls[0][1] is calls[1][1]
    with pytest.raises(FileExistsError):
        tool.create(output)
    tool.check(output)
    shared = calls[-1][1]
    assert shared is not None
    assert shared.cpu_seconds <= 60 and shared.maximum_evaluations <= 5000000
    assert binding_budget_ids[-1] == id(shared)
    changed = dict(value=False)
    monkeypatch.setattr(tool, "source_bindings", lambda:dict(source="changed" if changed["value"] else "fixed"))
    monkeypatch.setattr(tool, "verify_saved_bindings", lambda expected, root, **kw: expected == tool.source_bindings())
    def changing_validate(record, shared=None, **kwargs):
        result = validate(record, shared, **kwargs)
        changed["value"] = True
        return result
    monkeypatch.setattr(tool, "validate_record", changing_validate)
    with pytest.raises(RuntimeError, match="changed during authentication"):
        tool.check(output)
    changed["value"] = False
    payload = output/"ssprk3.json"
    payload.write_bytes(b"{}")
    with pytest.raises(ValueError, match="hash/source"):
        tool.check(output)


@pytest.mark.parametrize("key,value", (("maximum_evaluations", 50000001), ("maximum_evaluations", True),
    ("cpu_seconds", 300.1), ("cpu_seconds", float("nan")), ("evaluations", -1), ("started", float("nan"))))
def test_shared_budget_rejects_disabled_or_expanded_guards(key, value):
    shared = tool.Budget(maximum_evaluations=100000, cpu_seconds=30.)
    setattr(shared, key, value)
    with pytest.raises(ValueError):
        tool.shared_budget(shared)


def test_shared_budget_never_resets_the_clock_or_current_count():
    shared = tool.Budget(maximum_evaluations=10, cpu_seconds=30.)
    shared.check(9)
    start = shared.started
    assert tool.shared_budget(shared) is shared
    assert shared.started == start and shared.evaluations == 9
    shared.evaluations = 11
    with pytest.raises(RuntimeError, match="budget exhausted"):
        tool.shared_budget(shared)
