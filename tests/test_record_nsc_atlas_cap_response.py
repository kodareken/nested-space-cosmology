"""Synthetic cap evidence authentication; no scalar, metric or ray solve runs."""
from copy import deepcopy
import json
import math
import time
from types import SimpleNamespace

import pytest

import record_nsc_atlas_cap_response as tool


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
        x, q0 = .45+.1*t, 1.+.5*t
        delta = t+.1*t/(1-x*x)+.05*t
        rows.append({**dict.fromkeys(tool.SAMPLE_KEYS, 0.), "t": t, "x": x,
            "V": -2.5+.25*t, "T0_wall": q0-math.atanh(x), "q0": q0,
            "nu_E_wall": t, "nu_J_wall": .8*t, "delta_E_wall": t,
            "nu_E_potential_wall": .2*t, "integral_L_minus": t, "T1": t, "T1_s": t,
            "e": .1*t, "q0_s": 1., "V0_s": 2., "ray_integral": .05*t,
            "ray_positive_integral": .05*t, "ray_negative_integral": 0.,
            "ray_integrand_sample_min": 0., "ray_integrand_sample_max": .1*t,
            "Delta_C_fixed_s": delta, "Delta_C_fixed_V": delta-.1*t,
            "exterior_label_shift": .1*t, "V1": .2*t})
    solutions = [{"schema": tool.SOLUTION_SCHEMA, "scope": tool.SOLUTION_SCOPE,
        "domain": dict(domain), **{key: setting[key] for key in ("radial_order", "ray_order", "time_panels")},
        "clock_time_order": 4, "maximum_absolute_defects": dict.fromkeys(tool.DEFECTS, 0.),
        "cap_samples": deepcopy(rows), "wall_metric_evaluations": 4*setting["time_panels"]+9,
        "ray_evaluations": 9, "ray_nodes": 9*setting["ray_order"]} for setting in settings]
    return {"schema": tool.SCHEMA, "comparison": kind, "varied_quantity": varied,
        "status": tool.STATUS, "scope": tool.SCOPE, "order": tool.ORDER, "source_scale": 1.,
        "runtime": {"python": "3.14.8", "numpy": "2.5.1", "scipy": "1.17.1"},
        "comparison_base": base, "settings": settings, "solutions": solutions,
        "cap_differences": tool.differences(solutions),
        "refinement_is_rigorous_error_enclosure": False,
        "profile_evaluation": {"counts": dict(profile_requested=100, wall_points=50, second_points=10, cache_hits=20),
                               "maximum_inverse_bracket_width": 0.},
        "execution": dict(cpu_seconds=1., cpu_seconds_max=120., scalar_kernel_values=1000,
            scalar_kernel_values_max=30_000_000, metric_values=100, source_quadrature_nodes=100,
            rhs_calls=10, field_solves=2 if kind in ("field", "geometry") else 1,
            wall_solves=2 if kind in ("field", "geometry") else 1, cap_solves=2)}


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


def test_fixed_model_clock_scope_and_sampling_contract():
    assert set(tool.VARIATIONS) >= {"field", "geometry", "radial", "ray", "time"}
    for kind in tool.VARIATIONS:
        assert tool.validate_record(synthetic(kind), kind)
    for path, value in ((["settings", 0, "wall_source_order"], 16),
                       (["solutions", 1, "domain", "qb"], 1.6),
                       (["solutions", 0, "scope"], "corrected scalar amplitude"),
                       (["solutions", 0, "cap_samples", 4, "Delta_C_fixed_V"], 8.),
                       (["solutions", 0, "cap_samples", 4, "wall_einstein_mass_defect"], 1.),
                       (["solutions", 0, "cap_samples", 4, "ray_negative_integral"], 1.),
                       (["solutions", 0, "wall_metric_evaluations"], 0),
                       (["source_scale"], True), (["refinement_is_rigorous_error_enclosure"], True)):
        record, target = synthetic(), None
        target = record
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        with pytest.raises(ValueError):
            tool.validate_record(record, "radial")


def test_budget_counts_and_saved_differences_are_not_redeclarable():
    for key, value in (("cpu_seconds", 121), ("cpu_seconds_max", 121),
                       ("scalar_kernel_values", 30_000_001), ("rhs_calls", 12_001),
                       ("metric_values", -1), ("field_solves", True), ("cap_solves", 1)):
        record = synthetic()
        record["execution"][key] = value
        with pytest.raises(ValueError):
            tool.validate_record(record, "radial")
    for field in ("cap_differences", "profile_evaluation"):
        record = synthetic()
        if field == "cap_differences":
            record[field]["Delta_C_fixed_s"]["maximum_absolute_difference"] = 1.
        else:
            record[field]["maximum_inverse_bracket_width"] = float("nan")
        with pytest.raises(ValueError):
            tool.validate_record(record, "radial")


def test_actual_summary_contract_without_evolution():
    from recursive_horizons.nsc_atlas_cap_response import CapSolution
    record = synthetic()
    for index, saved in enumerate(record["solutions"]):
        solution = object.__new__(CapSolution)
        solution.domain = saved["domain"]
        solution.settings = {key: saved[key] for key in ("radial_order", "ray_order", "time_panels")}
        solution.maxima = saved["maximum_absolute_defects"]
        solution.wall_cache = dict.fromkeys(range(saved["wall_metric_evaluations"]))
        solution.ray_cache = dict.fromkeys(range(saved["ray_evaluations"]))
        solution.ray_nodes = saved["ray_nodes"]
        solution.budget = SimpleNamespace(rhs_calls=0, metric_values=0, source_nodes=0, started=time.process_time())
        rows = {row["t"]: row for row in saved["cap_samples"]}
        solution.sample = rows.__getitem__
        summary = solution.summary(tool.TIMES)
        for key in ("rhs_calls", "metric_values", "source_nodes", "process_cpu_seconds"):
            summary.pop(key)
        record["solutions"][index] = summary
    assert tool.validate_record(record, "radial")


def test_creation_preserves_existing_bytes_and_authentication_never_solves(isolated, monkeypatch):
    tool.create(isolated, "radial")
    before = {path.name: path.read_bytes() for path in isolated.iterdir()}
    monkeypatch.setattr(tool, "compare", lambda *_: pytest.fail("unexpected physical calculation"))
    assert tool.check(isolated)["cap_solves"] == 0
    with pytest.raises(FileExistsError, match="creation-only"):
        tool.create(isolated, "radial")
    assert before == {path.name: path.read_bytes() for path in isolated.iterdir()}
    path = tool.ROOT/"src/recursive_horizons/nsc_atlas_cap_field.py"
    path.write_bytes(path.read_bytes()+b"changed source")
    with pytest.raises(ValueError, match="source binding"):
        tool.check(isolated)


@pytest.mark.parametrize("change", ["bytes", "rehashed-clock", "rehashed-setting"])
def test_rehashed_tampering_still_fails(isolated, change):
    tool.create(isolated, "radial")
    path = isolated/"response.json"
    if change == "bytes":
        path.write_bytes(path.read_bytes()+b"\n")
    else:
        record = json.loads(path.read_bytes())
        if change == "rehashed-clock":
            record["solutions"][0]["cap_samples"][4]["Delta_C_fixed_V"] = 10.
        else:
            record["settings"][0]["geometry_knots"] = 513
        raw = tool.encoded(record)
        path.write_bytes(raw)
        manifest_path = isolated/"manifest.json"
        manifest = json.loads(manifest_path.read_bytes())
        manifest["outputs"]["response.json"] = {"sha256": tool.digest(raw), "bytes": len(raw)}
        manifest_path.write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError):
        tool.check(isolated)


def test_mid_calculation_source_change_cannot_create_evidence(isolated, monkeypatch):
    def changed(kind):
        (tool.ROOT/"src/recursive_horizons/nsc_atlas_cap_response.py").write_bytes(b"changed during solve")
        return synthetic(kind)
    monkeypatch.setattr(tool, "compare", changed)
    with pytest.raises(RuntimeError, match="inputs or runtime changed"):
        tool.create(isolated, "radial")
    assert not isolated.exists()
