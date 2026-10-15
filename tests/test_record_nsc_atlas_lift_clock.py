"""Clock polynomial and saved-boundary controls, with no physical solves."""
from copy import deepcopy
import json
import math

import numpy as np
import pytest

import record_nsc_atlas_lift_clock as tool
from recursive_horizons.nsc_atlas_cap_wall_response import reference_coefficients


def wall_inputs(panels):
    """A quartic clock rate has a known nonzero cubic interpolation defect."""
    nodes, basis = tool._gauss_antiderivatives()
    times = (np.arange(panels)[:, None]+(nodes[None, :]+1)/2)/panels
    values = .01*times**4
    prefix = .002*(np.arange(panels+1)/panels)**5
    rows = []
    for t in tool.predecessor.DIRECT_TIMES:
        x, V, U = .45+.1*t, -2.5+.25*t, -.30+.03*t
        state = reference_coefficients(x)
        v, bp, bm = [float(state[k]) for k in ("v", "beta_plus", "beta_minus")]
        Vs = -V/(1.6*(v-bp))
        ds_dt = .25/Vs
        panel = min(int(t*panels), panels-1)
        local = 2*(t-panel/panels)*panels-1
        eta = prefix[panel]+sum(p(local)*r for p, r in zip(basis, values[panel]))/(2*panels)
        eta_s = .01*t**4/ds_dt
        w, j = t*t, 2*t/ds_dt
        p = w/x
        rows.append(dict(t=t, U=U, V=V, x=x, v=v, beta_plus=bp, beta_minus=bm,
            Vdot=Vs, w=w, j=j, p=p, A=(j-v*p)/x, rho=.01*t**4,
            rho_s=.04*t**3/ds_dt, V1=V*eta, V1_s=Vs*eta+V*eta_s,
            eta=eta, eta_s=eta_s, ds_dt=ds_dt, reference_q_plus=1.+.5*t))
    return rows, dict(values=values.tolist(), prefix=prefix.tolist())


def synthetic(monkeypatch):
    domain = dict(U0=-.32, Ua=-.30, Ub=-.27, Va=-2.5, Vb=-2.25)
    rows = []
    for wall in wall_inputs(8)[0]:
        t = wall["t"]
        J, Js = t**3, 3*t*t/wall["ds_dt"]
        derivative = tool.wall_boundary_derivative(wall, J, Js)
        q = wall["reference_q_plus"]
        rows.append(dict(t=t, U=wall["U"], V=wall["V"], x=wall["x"],
            V0_s=wall["Vdot"], **derivative, b_U=derivative["b_V"]/.12,
            U_V=.12, J_minus=J, J_minus_s=Js, q_plus=q,
            q_minus=q-2*math.atanh(wall["x"]), wall_profile=wall["w"]))
    old = dict(domain=domain, direct_samples=rows,
        diagnostics={"interpolation_target": {"absolute_b_target": 1e-7, "source": "unchanged test target"}})
    old["solutions"] = [None, tool.packed_fit(domain, rows)]
    binding = {"path": tool.INPUT, "files": {"lift.json": {"sha256": "a"*64, "bytes": 100},
        "manifest.json": {"sha256": "b"*64, "bytes": 100}},
        "source_validation": "record_nsc_atlas_local_lift.check", "physical_solves": 0}
    monkeypatch.setattr(tool, "input_snapshot", lambda: (deepcopy(old), deepcopy(binding)))
    solutions = [tool.rebuild_solution(old, panels, *wall_inputs(panels)) for panels in tool.PANELS]
    record = dict(schema=tool.SCHEMA, scope=tool.SCOPE, input=binding, domain=domain,
        runtime=tool.runtime_versions(), comparison_base=dict(tool.BASE), varied_quantity="wall_time_panels",
        settings=[{**tool.BASE, "wall_time_panels": p} for p in tool.PANELS],
        solutions=solutions, diagnostics=tool.comparison_diagnostics(old, solutions),
        reused_ray_quantities=["J_minus", "J_minus_s"], uniform_error_enclosure=False,
        execution=dict(cpu_seconds=1., cpu_seconds_max=120., scalar_kernel_values=1000,
            scalar_kernel_values_max=30_000_000, metric_values=1000, source_quadrature_nodes=100,
            field_rhs_calls=100, wall_rhs_calls=200, field_solves=1, wall_solves=2,
            reused_boundary_jets=81, new_cap_queries=0, new_ray_quadrature_nodes=0))
    return tool._plain(record), old


def test_quartic_clock_rate_and_boundary_slope_gap(monkeypatch):
    record, old = synthetic(monkeypatch)
    assert tool.validate_record(record)
    nodes, _ = tool._gauss_antiderivatives()
    for solution in record["solutions"]:
        panels = solution["wall_time_panels"]
        for row in solution["clock_diagnostics"]:
            t = row["t"]
            panel = min(int(t*panels), panels-1)
            physical_nodes = (panel+(nodes+1)/2)/panels
            expected = -.01*np.prod(t-physical_nodes)
            assert row["eta_t_gap"] == pytest.approx(expected, abs=3e-17)
            assert row["b_U_clock_gap"] == pytest.approx(row["predicted_b_U_clock_gap"], abs=8e-12)
        assert solution["maximum_absolute_clock_diagnostics"]["b_U_clock_gap"] > 1e-9
        for changed, original in zip(solution["boundary_samples"], old["direct_samples"]):
            assert (changed["J_minus"], changed["J_minus_s"]) == (original["J_minus"], original["J_minus_s"])
        assert len(solution["lift"]["knot_samples"]) == 41
        assert len(solution["lift"]["validation_samples"]) == 40
        assert solution["lift"]["heldout_used_for_training"] is False


def test_restore_with_live_geometry_preserves_budget_and_quiet_anchor(monkeypatch):
    record, _ = synthetic(monkeypatch)
    from types import SimpleNamespace
    budget = SimpleNamespace(check=lambda: None)
    geometry = {"wall": lambda V: -.30+.12*(np.asarray(V)+2.5)}
    lift = tool.restore_lift(record["domain"], record["solutions"][1]["boundary_samples"],
        geometry=geometry, budget=budget, validation=False)
    assert lift.budget is budget and lift.geometry is geometry
    assert lift.validation_rows == []
    assert lift.b_at_U(-.31) == lift.B(-.31) == 0
    assert float(lift.b(-2.37)) == pytest.approx(float(lift.b_at_U(geometry["wall"](-2.37))))


@pytest.mark.parametrize("key", ["x", "w", "rho", "reference_q_plus"])
def test_changed_leading_data_cannot_reuse_original_rays(monkeypatch, key):
    _, old = synthetic(monkeypatch)
    walls, clock = wall_inputs(16)
    walls[20][key] += 1e-4
    with pytest.raises(ValueError, match="reference/leading"):
        tool.rebuild_solution(old, 16, walls, clock)


def test_clock_panel_and_continuum_definitions_are_checked(monkeypatch):
    _, old = synthetic(monkeypatch)
    walls, clock = wall_inputs(16)
    clock["prefix"][4] += 1e-4
    with pytest.raises(ValueError, match="prefix"):
        tool.rebuild_solution(old, 16, walls, clock)
    walls, clock = wall_inputs(16)
    walls[20]["V1_s"] += .01
    with pytest.raises(ValueError, match="continuum slope"):
        tool.rebuild_solution(old, 16, walls, clock)


def test_saved_numbers_target_and_budgets_are_recomputed(monkeypatch):
    record, _ = synthetic(monkeypatch)
    paths = [(["solutions", 1, "boundary_samples", 25, "J_minus_s"], 30.),
        (["solutions", 1, "clock_diagnostics", 25, "b_U_clock_gap"], 30.),
        (["solutions", 1, "lift", "validation_samples", 3, "B_lift"], 30.),
        (["diagnostics", "interpolation_target", "absolute_b_target"], 1.),
        (["execution", "cpu_seconds"], 121.), (["execution", "new_cap_queries"], 1),
        (["execution", "field_solves"], True), (["execution", "wall_rhs_calls"], 12001)]
    for path, value in paths:
        altered = deepcopy(record)
        target = altered
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        with pytest.raises(ValueError):
            tool.validate_record(altered)


@pytest.fixture
def isolated(tmp_path, monkeypatch):
    record, _ = synthetic(monkeypatch)
    root = tmp_path/"repository"
    for name in tool.SOURCES:
        path = root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(("fixture: "+name).encode())
    monkeypatch.setattr(tool, "ROOT", root)
    monkeypatch.setattr(tool, "compare", lambda: deepcopy(record))
    return root/"result"


def test_creation_check_has_no_physical_work_and_detects_source_drift(isolated, monkeypatch):
    tool.create(isolated)
    before = {p.name: p.read_bytes() for p in isolated.iterdir()}
    monkeypatch.setattr(tool, "compare", lambda: pytest.fail("physical comparison during check"))
    monkeypatch.setattr(tool, "_prepare_shared", lambda: pytest.fail("field solve during check"))
    assert tool.check(isolated)["physical_solves"] == 0
    with pytest.raises(FileExistsError):
        tool.create(isolated)
    assert before == {p.name: p.read_bytes() for p in isolated.iterdir()}
    (tool.ROOT/"scripts/record_nsc_atlas_lift_clock.py").write_bytes(b"source drift")
    with pytest.raises(ValueError, match="binding"):
        tool.check(isolated)


def test_rehashed_false_polynomial_fails(isolated):
    tool.create(isolated)
    path = isolated/"clock.json"
    record = json.loads(path.read_bytes())
    record["solutions"][1]["lift"]["cells"][5]["coefficients_in_cell_coordinate"][0] += 1
    raw = tool.encoded(record)
    path.write_bytes(raw)
    manifest_path = isolated/"manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    manifest["outputs"]["clock.json"] = dict(sha256=tool.digest(raw), bytes=len(raw))
    manifest_path.write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError, match="polynomial"):
        tool.check(isolated)


def test_input_snapshot_authenticates_original_files(tmp_path, monkeypatch):
    path = tmp_path/tool.INPUT
    path.mkdir(parents=True)
    (path/"lift.json").write_bytes(tool.encoded({"test": 1.}))
    (path/"manifest.json").write_bytes(b"fixture manifest")
    calls = []
    monkeypatch.setattr(tool, "ROOT", tmp_path)
    monkeypatch.setattr(tool.predecessor, "check", lambda p: calls.append(p))
    record, binding = tool.input_snapshot()
    assert calls == [path] and record == {"test": 1.}
    assert binding["physical_solves"] == 0
    for name, value in binding["files"].items():
        assert value == dict(sha256=tool.digest((path/name).read_bytes()), bytes=(path/name).stat().st_size)
