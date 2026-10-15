"""Synthetic boundary evidence authentication; no prepared-field solve runs."""
from copy import deepcopy
import json
import math
from types import SimpleNamespace
import time

import pytest

import record_nsc_atlas_scalar_boundary as tool


def synthetic(kind="radial"):
    varied, low, high = tool.VARIATIONS[kind]
    base = {**tool.BASE, **tool.COMPARISON_BASE_OVERRIDES.get(kind, {})}
    settings = [{**base, varied: value} for value in (low, high)]
    domain = dict(U0=-.32, Ua=-.30, Ub=-.27, Va=-2.5, Vb=-2.25,
                  qa=1., qb=1.5, xa=.45, xb=.55)
    domain.update(Ta=1.-math.atanh(.45), Tb=1.5-math.atanh(.55),
                  qmin=1.-2*math.atanh(.45))
    rows = []
    for t in tool.TIMES:
        x, q = .45+.1*t, 1.+.5*t
        u = q-2*math.atanh(x)
        p, w, rho, V1, v, Vdot, bp = t, x*t, .1*t, .2*t, 1., 2., .25
        bm, d = math.sqrt(2-x*x), v-bp
        qm_s = (bm-v)/(1-x*x)
        e = rho-x*p*p/12
        offset = e-v*V1/Vdot
        hessian = .8*w/(x**3*d*d)
        J = -.02*t
        backscatter, displacement = 2*qm_s*J/d, -offset*hessian
        c1 = backscatter+displacement
        rows.append({**dict.fromkeys(tool.SAMPLE_KEYS, 0.), "t": t, "V": -2.5+.25*t,
            "x": x, "q_plus": q, "q_minus": u, "T_wall": (u+q)/2,
            "p_wall": p, "w": w, "rho": rho, "V1": V1, "v": v, "V0_s": Vdot,
            "beta_plus": bp, "beta_minus": bm, "d": d, "q_minus_s": qm_s,
            "e": e, "displacement_at_fixed_V": offset, "c0_rr": hessian,
            "J_minus": J, "J_minus_r": J, "negative_integral": J, "absolute_integral": -J,
            "K1_sample_min": 0., "K1_sample_max": t,
            "source_sample_min": -t, "source_sample_max": 0., "c1_r": c1,
            "backscatter_term": backscatter, "displacement_term": displacement,
            "wall_flux_defect": d*(c1+offset*hessian)/2-qm_s*J,
            "no_return_margin": domain["qa"]-u,
            "quiet_centre_time_margin": domain["qa"]-(u+q)/2,
            "initial_time_margin": (u+domain["qa"])/2-domain["Ta"],
            "last_profile_margin": domain["qb"]-q, "first_profile_margin": u-domain["qmin"],
            "active_radius_lower": math.tanh((domain["qa"]-u)/2)})
    solutions = []
    for setting in settings:
        values = deepcopy(rows)
        for row in values:
            row["ray_nodes"] = row["ray_nodes_r"] = setting["ray_order"] if row["t"] else 0
        solutions.append({"schema": tool.SOLUTION_SCHEMA, "scope": tool.SOLUTION_SCOPE,
            "domain": dict(domain), **{key: setting[key] for key in ("radial_order", "ray_order")},
            "boundary_samples": values,
            "maximum_absolute_defects": {key: max(abs(row[key]) for row in values) for key in tool.DEFECTS},
            "minimum_domain_margins": {key: min(row[key] for row in values) for key in tool.MARGINS},
            "boundary_evaluations": len(tool.TIMES),
            "ray_nodes_q": sum(row["ray_nodes"] for row in values),
            "ray_nodes_r": sum(row["ray_nodes_r"] for row in values)})
    return {"schema": tool.SCHEMA, "comparison": kind, "varied_quantity": varied,
        "status": tool.STATUS, "scope": tool.SCOPE, "order": tool.ORDER, "source_scale": 1.,
        "runtime": {"python": "3.14.8", "numpy": "2.5.1", "scipy": "1.17.1"},
        "comparison_base": base, "settings": settings, "solutions": solutions,
        "boundary_differences": tool.differences(solutions),
        "refinement_is_rigorous_error_enclosure": False,
        "profile_evaluation": {"counts": dict(profile_requested=100, wall_points=50, second_points=0, cache_hits=20),
                               "maximum_inverse_bracket_width": 0.},
        "execution": dict(cpu_seconds=1., cpu_seconds_max=120., scalar_kernel_values=1000,
            scalar_kernel_values_max=30_000_000, metric_values=100, source_quadrature_nodes=100,
            rhs_calls=10, field_solves=2 if kind == "field" else 1,
            wall_solves=2 if kind == "field" else 1, boundary_evaluators=2)}


@pytest.fixture
def isolated(tmp_path, monkeypatch):
    root = tmp_path/"synthetic-repository"
    for name in tool.SOURCES:
        path = root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(("fixture: "+name).encode())
    monkeypatch.setattr(tool, "ROOT", root)
    monkeypatch.setattr(tool, "runtime_versions", lambda: synthetic()["runtime"])
    monkeypatch.setattr(tool, "compare", synthetic)
    return root/"synthetic-result"


def test_fixed_model_identity_domain_and_sampling_contract():
    assert set(tool.VARIATIONS) == {"field", "radial", "ray", "radial_fine", "ray_fine"}
    for kind in tool.VARIATIONS:
        assert tool.validate_record(synthetic(kind), kind)
    for path, value in ((["settings", 0, "wall_source_order"], 16),
                       (["solutions", 1, "domain", "qb"], 1.6),
                       (["solutions", 0, "scope"], "evolved waveform"),
                       (["solutions", 0, "boundary_samples", 4, "c1_r"], 8.),
                       (["solutions", 0, "boundary_samples", 4, "J_minus_r"], 1.),
                       (["solutions", 0, "boundary_samples", 4, "negative_integral"], 1.),
                       (["solutions", 0, "boundary_samples", 4, "no_return_margin"], -1.),
                       (["solutions", 0, "ray_nodes_q"], 0),
                       (["source_scale"], True), (["refinement_is_rigorous_error_enclosure"], True)):
        record = synthetic()
        target = record
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        with pytest.raises(ValueError):
            tool.validate_record(record, "radial")


def test_budget_counts_and_saved_differences_cannot_be_redeclared():
    for key, value in (("cpu_seconds", 121), ("cpu_seconds_max", 121),
                       ("scalar_kernel_values", 30_000_001), ("rhs_calls", 12_001),
                       ("metric_values", -1), ("field_solves", True), ("boundary_evaluators", 1)):
        record = synthetic()
        record["execution"][key] = value
        with pytest.raises(ValueError):
            tool.validate_record(record, "radial")
    record = synthetic()
    record["boundary_differences"]["c1_r"]["maximum_absolute_difference"] = 1.
    with pytest.raises(ValueError, match="saved differences"):
        tool.validate_record(record, "radial")
    record = synthetic()
    record["profile_evaluation"]["maximum_inverse_bracket_width"] = float("nan")
    with pytest.raises(ValueError):
        tool.validate_record(record, "radial")


def test_actual_summary_contract_without_solving():
    from recursive_horizons.nsc_atlas_scalar_boundary import ScalarBoundary
    record = synthetic()
    for index, saved in enumerate(record["solutions"]):
        solution = object.__new__(ScalarBoundary)
        solution.domain = saved["domain"]
        solution.settings = {key: saved[key] for key in ("radial_order", "ray_order")}
        solution.maxima, solution.minimum_margins = saved["maximum_absolute_defects"], saved["minimum_domain_margins"]
        solution.cache = {row["t"]: row for row in saved["boundary_samples"]}
        solution.ray_nodes_q, solution.ray_nodes_r = saved["ray_nodes_q"], saved["ray_nodes_r"]
        solution.budget = SimpleNamespace(rhs_calls=0, metric_values=0, source_nodes=0, started=time.process_time())
        solution.cap_field = SimpleNamespace(budget=SimpleNamespace(kernel_values=0), evaluation_counts={})
        solution.sample = solution.cache.__getitem__
        summary = solution.summary(tool.TIMES)
        for key in ("rhs_calls", "metric_values", "source_nodes", "process_cpu_seconds",
                    "scalar_kernel_values", "profile_counts"):
            summary.pop(key)
        record["solutions"][index] = summary
    assert tool.validate_record(record, "radial")


@pytest.mark.parametrize("kind", ["field", "ray"])
def test_comparison_reuses_compatible_inputs_and_one_cumulative_budget(monkeypatch, kind):
    from recursive_horizons import (nsc_atlas_prepared_response as prepared,
        nsc_atlas_wedge_field as wedge, nsc_atlas_wall_evolution as wall,
        nsc_atlas_characteristic_initialization as geometry,
        nsc_atlas_cap_field as cap, nsc_atlas_scalar_boundary as boundary)
    inputs, evaluators = [], []
    def field(nodes, knots, *, budget):
        budget.check(10)
        result = SimpleNamespace(budget=budget, nodes=nodes)
        inputs.append(result)
        return result
    def solve(wedge_field, _, *, budget, **settings):
        assert budget.external is wedge_field.prepared_field.budget
        assert budget.started == budget.external.started
        budget.check(metric=3, rhs=True)
        return SimpleNamespace(budget=budget, wedge=wedge_field)
    def cap_field(field, **settings):
        return SimpleNamespace(prepared_field=field, budget=field.budget,
            evaluation_counts=dict(profile_requested=100, wall_points=50, second_points=0, cache_hits=20),
            max_inverse_bracket_width=0.)
    class Boundary:
        def __init__(self, cap_field, wall_solution, *, budget, **settings):
            assert budget is wall_solution.budget
            assert cap_field.prepared_field is wall_solution.wedge.prepared_field
            self.index = len(evaluators)
            evaluators.append((cap_field, wall_solution, budget))
            budget.check(metric=4)
            budget.external.check(5)
        def summary(self, times):
            assert times == tool.TIMES
            return deepcopy(synthetic(kind)["solutions"][self.index])
    monkeypatch.setattr(prepared, "solve_prepared_field", field)
    monkeypatch.setattr(wedge, "WedgeField", lambda field: SimpleNamespace(prepared_field=field))
    monkeypatch.setattr(wall, "solve_wall_evolution", solve)
    monkeypatch.setattr(geometry, "source_geometry", lambda: {})
    monkeypatch.setattr(cap, "CapLeadingField", cap_field)
    monkeypatch.setattr(boundary, "ScalarBoundary", Boundary)
    record = tool.compare(kind)
    assert tool.validate_record(record, kind)
    assert len({id(value.budget) for value in inputs}) == 1
    assert evaluators[0][2] is evaluators[1][2]
    assert (evaluators[0][0] is evaluators[1][0]) == (kind == "ray")
    assert (evaluators[0][1] is evaluators[1][1]) == (kind == "ray")
    assert record["execution"]["scalar_kernel_values"] == (30 if kind == "field" else 20)
    assert record["execution"]["metric_values"] == (14 if kind == "field" else 11)


def test_creation_only_source_binding_and_no_replay(isolated, monkeypatch):
    tool.create(isolated, "radial")
    before = {path.name: path.read_bytes() for path in isolated.iterdir()}
    monkeypatch.setattr(tool, "compare", lambda *_: pytest.fail("unexpected physical calculation"))
    assert tool.check(isolated)["physical_solves"] == 0
    with pytest.raises(FileExistsError, match="creation-only"):
        tool.create(isolated, "radial")
    assert before == {path.name: path.read_bytes() for path in isolated.iterdir()}
    path = tool.ROOT/"src/recursive_horizons/nsc_atlas_scalar_boundary.py"
    path.write_bytes(path.read_bytes()+b"changed source")
    with pytest.raises(ValueError, match="source binding"):
        tool.check(isolated)


@pytest.mark.parametrize("change", ["bytes", "rehashed-boundary", "rehashed-setting"])
def test_rehashed_tampering_fails(isolated, change):
    tool.create(isolated, "radial")
    path = isolated/"boundary.json"
    if change == "bytes":
        path.write_bytes(path.read_bytes()+b"\n")
    else:
        record = json.loads(path.read_bytes())
        if change == "rehashed-boundary":
            record["solutions"][0]["boundary_samples"][4]["c1_r"] = 10.
        else:
            record["settings"][0]["geometry_knots"] = 513
        raw = tool.encoded(record)
        path.write_bytes(raw)
        manifest_path = isolated/"manifest.json"
        manifest = json.loads(manifest_path.read_bytes())
        manifest["outputs"]["boundary.json"] = {"sha256": tool.digest(raw), "bytes": len(raw)}
        manifest_path.write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError):
        tool.check(isolated)


def test_mid_calculation_source_change_prevents_record_creation(isolated, monkeypatch):
    def changed(kind):
        (tool.ROOT/"src/recursive_horizons/nsc_atlas_scalar_boundary.py").write_bytes(b"changed during solve")
        return synthetic(kind)
    monkeypatch.setattr(tool, "compare", changed)
    with pytest.raises(RuntimeError, match="inputs or runtime changed"):
        tool.create(isolated, "radial")
    assert not isolated.exists()
