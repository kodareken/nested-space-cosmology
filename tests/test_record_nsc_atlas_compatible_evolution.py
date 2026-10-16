"""Manufactured stored histories and stubbed marches; no physical evolution."""
import copy
import json
import time

import numpy as np
import pytest

import record_nsc_atlas_compatible_evolution as tool


LIVE_LOAD_PARENT = tool.load_parent
AUTH_CHARGES = 17


@pytest.fixture(scope="module")
def parent():
    # Reuse immutable inputs, not their physical evolution, for a synthetic
    # record. Authentication is separately exercised with an explicit stub.
    return json.loads((tool.ROOT/tool.previous.RELATIVE_OUTPUT/"coupled.json").read_bytes())


def input_identity():
    return dict(path=tool.previous.RELATIVE_OUTPUT+"/coupled.json", sha256="fixture-parent",
                manifest_sha256="fixture-manifest", authentication=dict(authentication_evaluations=AUTH_CHARGES))


def fake_parent(parent):
    def load(budget):
        budget.check(AUTH_CHARGES)
        return parent, input_identity(), AUTH_CHARGES
    return load


@pytest.fixture(autouse=True)
def bounded_stored_inputs(parent, monkeypatch):
    started = time.process_time()
    monkeypatch.setattr(tool, "load_parent", fake_parent(parent))
    yield
    assert time.process_time()-started < 30.


@pytest.fixture(scope="module")
def record(parent):
    budget = tool.Budget(maximum_evaluations=8000000, cpu_seconds=30.)
    x = np.asarray(parent["gauge"]["rows"]["x"])
    result = dict(schema=tool.SCHEMA, status=tool.STATUS, scope=tool.SCOPE, settings=tool.SETTINGS,
        runtime=tool.runtime_versions(), parent_input=input_identity(), geometry=parent["geometry"],
        source_bounds=parent["source_bounds"], cases=[], refinement_is_rigorous_error_enclosure=False)
    end = AUTH_CHARGES
    for label, n, steps in tool.CASES:
        alpha = parent["initial_data"][label]["alpha"]
        samples, eta = [], np.linspace(0., 1., n+1)
        for j in range(0, 513, 4):
            t = float(x[j]-.4)
            response = max(0., float(x[j]-parent["geometry"]["support_x"][0]))**2
            state = tool.quiet_departures(exterior_cells=n, cap_cells=n)
            state["departure"]["r"] = alpha*t*eta[1:-1]*(1-eta[1:-1])/n
            state["departure"]["L"] = alpha*t*eta[1:-1]/n
            state["cap"]["Iq"][:] = response*.01/n
            state["cap"]["Ip"][:] = response*.002/n
            state["departure"]["c"] = response*.01*eta[1:-1]
            state["wall"]["r"] = alpha*t*.1
            state["wall"]["L"] = alpha*t*.2
            physical = tool.physical_view(state, float(x[j]),
                tool.subset(parent["initial_data"][label]["snapshot_rows"], j),
                tool.subset(parent["gauge"]["rows"], j), eta, alpha, budget)
            samples.append(dict(x=float(x[j]), state=state, physical=physical))
        stage = tool.compatible_stage(samples[-1]["state"], x=x[-1], eta=eta, y_faces=eta,
            source=tool.subset(parent["initial_data"][label]["snapshot_rows"], -1),
            gauge=tool.subset(parent["gauge"]["rows"], -1), alpha=alpha, tube=tool.TUBE, budget=budget)
        final = tool.final_summary(stage, budget)
        stages = 2*steps+1
        charge = stages*(36*n+11)+3*n
        result["cases"].append(dict(amplitude=label, cells=n, steps=steps, alpha=alpha,
            snapshots=samples, final=final, execution=dict(budget_start=end, budget_end=end+charge,
                stages=stages, source_cache_reads=stages, gauge_cache_reads=stages,
                max_courant=.2, retries=0)))
        end += charge
    before = budget.evaluations
    result["diagnostics"] = tool.diagnostics(result, parent, budget)
    charges = budget.evaluations-before
    result["execution"] = dict(evaluations=end+charges, cpu_seconds=.02,
        inherited_authentication_evaluations=AUTH_CHARGES, diagnostic_evaluations=charges,
        source_integrations=0, evolution_source_snapshot_constructions=0,
        inherited_source_snapshot_reconstructions=4, gauge_quadratures=0, coupled_evolutions=5,
        time_steps=sum(k for _, _, k in tool.CASES), stage_evaluations=sum(2*k+1 for _, _, k in tool.CASES),
        inherited_authentications=1, old_rectangle_solves=0)
    result = tool.plain(result)
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(tool, "load_parent", fake_parent(parent))
        authentication = tool.validate_record(result, require_creation=False)
    result["creation_execution"] = dict(calculation_evaluations=result["execution"]["evaluations"],
        validation_evaluations=authentication["authentication_evaluations"],
        total_evaluations=result["execution"]["evaluations"]+authentication["authentication_evaluations"],
        cpu_seconds=.2)
    return result


def test_synthetic_record_checks_without_source_gauge_or_time_solve(record, monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("physical solve during saved authentication")
    monkeypatch.setattr(tool, "evolve_compatible", forbidden)
    monkeypatch.setattr(tool.previous, "ExactSource", forbidden)
    monkeypatch.setattr(tool.previous, "ReferenceWallGauge", forbidden)
    result = tool.validate_record(record)
    assert result["source_integrations"] == result["time_evolutions"] == result["reference_quadratures"] == 0
    assert result["final_cap_radial_reconstructions"] == 5
    assert result["authentication_evaluations"] == record["creation_execution"]["validation_evaluations"]


@pytest.mark.parametrize("change", ["settings", "runtime", "parent", "geometry", "case", "raw", "physical",
    "wall_c", "initial", "nonfinite", "clock", "final", "counts", "creation", "raw_constraint",
    "departure_constraint", "old_comparison", "cap_characteristic", "restriction", "source_defect", "overlap"])
def test_rehashed_content_tampering_is_rejected(record, change):
    value = copy.deepcopy(record)
    if change == "settings": value["settings"]["cap_quadrature_order"] = 8
    elif change == "runtime": value["runtime"]["numpy"] = ""
    elif change == "parent": value["parent_input"]["sha256"] = "wrong"
    elif change == "geometry": value["geometry"]["Ua"] += .01
    elif change == "case": value["cases"][1]["cells"] = 16
    elif change == "raw": value["cases"][0]["snapshots"][30]["state"]["departure"]["r"][3] += .01
    elif change == "physical": value["cases"][0]["snapshots"][30]["physical"]["increment"]["r"][3] += .01
    elif change == "wall_c": value["cases"][0]["snapshots"][30]["physical"]["wall"]["c"] += .01
    elif change == "initial": value["cases"][0]["snapshots"][0]["state"]["cap"]["Ip"][0] = .01
    elif change == "nonfinite": value["cases"][0]["snapshots"][30]["state"]["cap"]["Iq"][0] = np.nan
    elif change == "clock": value["cases"][0]["snapshots"][30]["physical"]["wall"]["T"] = 0.
    elif change == "final": value["cases"][0]["final"]["cap_faces"]["A_faces"][3] += .01
    elif change == "counts": value["cases"][0]["execution"]["budget_end"] += 1
    elif change == "creation": value["creation_execution"]["total_evaluations"] -= 1
    elif change == "raw_constraint": value["diagnostics"]["histories"][0]["raw129"]["E_Ubar"][3][3] += .01
    elif change == "departure_constraint": value["diagnostics"]["histories"][0]["departure129"]["E_x"][3][3] += .01
    elif change == "old_comparison": value["diagnostics"]["successor_minus_v1"][0]["wall_arrays"]["r"][3] += .01
    elif change == "cap_characteristic": value["diagnostics"]["cap_characteristics"][0]["arrays"]["a_minus"][30] += .01
    elif change == "restriction": value["diagnostics"]["cap_restriction"] = "plain averages"
    elif change == "source_defect": value["diagnostics"]["source_interpolation_withheld_257_odd"]["zero"]["h"]["maximum_absolute"] += .01
    elif change == "overlap": value["diagnostics"]["original_rectangle_overlap"]["cases"][0]["difference_linear"]["r"][2][2] += .01
    with pytest.raises(ValueError):
        tool.validate_record(value)


def test_restriction_sums_integrals_without_averaging():
    table = dict(exterior={"r": np.arange(5)[None, :]}, wall={"r": [1.]},
                 cap={"Ip": np.array([[1., 2., 3., 4.]]), "Iq": np.array([[2., 4., 6., 8.]])})
    restricted = tool._restricted(table)
    np.testing.assert_array_equal(restricted["cap"]["Ip"], [[3., 7.]])
    np.testing.assert_array_equal(restricted["cap"]["Iq"], [[6., 14.]])
    np.testing.assert_array_equal(restricted["exterior"]["r"], [[0, 2, 4]])


def test_canonical_primitive_gives_nonlinear_wall_c(parent):
    budget = tool.Budget(maximum_evaluations=100, cpu_seconds=30.)
    gauge = tool.subset(parent["gauge"]["rows"], 0)
    source = tool.subset(parent["initial_data"]["finite"]["snapshot_rows"], 0)
    state = tool.quiet_departures(exterior_cells=4, cap_cells=4)
    state["cap"]["Iq"][:] = [.1, .2, -.1, .3]
    alpha, R = .6, gauge["r"]
    view = tool.physical_view(state, .4, source, gauge, np.linspace(0., 1., 5), alpha, budget)
    f = sum(state["cap"]["Iq"])/R
    expected = R*np.sinh(np.sqrt(alpha/6)*f)/np.sqrt(alpha/6)
    assert view["wall"]["c"] == pytest.approx(expected, abs=2e-16)
    assert abs(view["wall"]["c"]-sum(state["cap"]["Iq"])) > .005


def test_live_parent_authentication_receives_caller_budget(parent, monkeypatch):
    budget = tool.Budget(maximum_evaluations=100, cpu_seconds=30.)
    def validate(value, owner):
        assert owner is budget
        assert value["schema"] == tool.previous.SCHEMA
        owner.check(13)
        return {"authentication_evaluations": 13}
    monkeypatch.setattr(tool.previous, "validate_record", validate)
    inherited, identity, charge = LIVE_LOAD_PARENT(budget)
    assert charge == budget.evaluations == 13
    assert identity["sha256"] == tool.digest((tool.ROOT/tool.previous.RELATIVE_OUTPUT/"coupled.json").read_bytes())
    assert inherited["geometry"] == parent["geometry"]


def fake_calculation(record):
    def calculate(budget):
        budget.check(record["execution"]["evaluations"])
        value = copy.deepcopy(record)
        value.pop("creation_execution", None)
        value["execution"]["cpu_seconds"] = 0.
        return value
    return calculate


def test_creation_only_roundtrip_manifest_and_hash(record, tmp_path, monkeypatch):
    monkeypatch.setattr(tool, "calculate", fake_calculation(record))
    directory = tmp_path/"output"
    result = tool.create(directory)
    assert result["creation_execution"]["total_evaluations"] == record["creation_execution"]["total_evaluations"]
    assert tool.check(directory)["time_evolutions"] == 0
    with pytest.raises(FileExistsError):
        tool.create(directory)
    path = directory/"compatible.json"
    path.write_bytes(path.read_bytes()+b" ")
    with pytest.raises(ValueError, match="hash"):
        tool.check(directory)


def test_creation_authentication_cannot_reset_allowance(record, tmp_path, monkeypatch):
    seen = []
    def calculate(budget):
        seen.append(budget)
        budget.check(7999999)
        return copy.deepcopy(record)
    def validate(value, budget, **kwargs):
        assert budget is seen[0]
        budget.check(2)
    monkeypatch.setattr(tool, "calculate", calculate)
    monkeypatch.setattr(tool, "validate_record", validate)
    monkeypatch.setattr(tool, "ROOT", tmp_path)
    monkeypatch.setattr(tool, "source_bindings", lambda: {"fixed": "hash"})
    with pytest.raises(RuntimeError, match="budget exhausted") as error:
        tool.create(tmp_path/"absent")
    assert seen[0].evaluations == 8000001
    assert not (tmp_path/"absent").exists()
    assert json.loads(open(error.value.failure_receipt).read())["status"] == "FAILED_NOT_COMPLETED"


def test_source_drift_is_detected_before_result_write(record, tmp_path, monkeypatch):
    calls = [0]
    def bindings():
        calls[0] += 1
        return {"file": str(calls[0])}
    monkeypatch.setattr(tool, "source_bindings", bindings)
    monkeypatch.setattr(tool, "calculate", fake_calculation(record))
    monkeypatch.setattr(tool, "ROOT", tmp_path)
    with pytest.raises(RuntimeError, match="changed"):
        tool.create(tmp_path/"absent")
    assert not (tmp_path/"absent").exists()


def test_source_drift_during_authentication_is_detected(record, tmp_path, monkeypatch):
    changed = [False]
    def validate(value, budget, **kwargs):
        changed[0] = True
        return {"authentication_evaluations": 0}
    monkeypatch.setattr(tool, "source_bindings", lambda: {"file": "after" if changed[0] else "before"})
    monkeypatch.setattr(tool, "calculate", fake_calculation(record))
    monkeypatch.setattr(tool, "validate_record", validate)
    monkeypatch.setattr(tool, "ROOT", tmp_path)
    with pytest.raises(RuntimeError, match="changed during authentication"):
        tool.create(tmp_path/"absent")
    assert not (tmp_path/"absent").exists()


def test_calculate_reuses_caches_and_stops_at_explicit_stub(parent, monkeypatch):
    budget = tool.Budget(maximum_evaluations=8000000, cpu_seconds=30.)
    def stop(initial, x, **kwargs):
        assert kwargs["budget"] is budget
        assert len(x) == 257
        assert kwargs["source"](.4)["r"] == .4
        assert kwargs["gauge"](.4)["r"] == .45
        with pytest.raises(ValueError, match="absent"):
            kwargs["source"](np.nextafter(.4, 1.))
        raise RuntimeError("intentional no-march stop")
    monkeypatch.setattr(tool, "evolve_compatible", stop)
    with pytest.raises(RuntimeError, match="no-march") as error:
        tool.calculate(budget)
    assert budget.evaluations == AUTH_CHARGES
    assert error.value.compatible_failure["phase"]["name"] == "compatible_evolution"
