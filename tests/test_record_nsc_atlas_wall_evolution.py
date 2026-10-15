"""Synthetic evidence-contract checks; no scalar or wall calculation runs."""
from copy import deepcopy
import json
import time
from types import SimpleNamespace

import pytest

import record_nsc_atlas_wall_evolution as tool


def synthetic(kind="radial"):
    varied, low, high = tool.VARIATIONS[kind]
    settings = [{**tool.BASE, varied: value} for value in (low, high)]
    domain = {"U0": -.32, "Ua": -.30, "Ub": -.27, "Va": -2.5, "Vb": -2.25}
    rows = [{**dict.fromkeys(tool.SAMPLE_KEYS, t), "U": -.30+.03*t,
             "V": -2.5+.25*t, "x": .45+.1*t, "mass_path_defect": 0.} for t in tool.TIMES]
    solutions = [{"schema": tool.SOLUTION_SCHEMA, "scope": tool.SCOPE, "domain": dict(domain),
                  **{key: value for key, value in setting.items() if key not in ("field_nodes", "geometry_knots")},
                  "source_scale": 1., "clock_time_order": 4,
                  "maximum_absolute_defects": dict.fromkeys(tool.DEFECTS, 0.),
                  "wall_samples": deepcopy(rows), "radial_evaluations": 4*setting["time_panels"]+9}
                 for setting in settings]
    return {"schema": tool.SCHEMA, "comparison": kind, "varied_quantity": varied,
            "status": tool.STATUS, "order": tool.ORDER, "interpretation": tool.INTERPRETATION,
            "runtime": {"python": "3.14.8", "numpy": "2.5.1", "scipy": "1.17.1"},
            "settings": settings, "solutions": solutions, "wall_differences": tool.differences(solutions),
            "refinement_is_rigorous_error_enclosure": False,
            "execution": {"cpu_seconds": 1., "cpu_seconds_max": 120.,
                "scalar_kernel_values": 100, "scalar_kernel_values_max": 30_000_000,
                "metric_values": 100, "source_quadrature_nodes": 100, "rhs_calls": 10,
                "field_solves": 2 if kind in ("field", "geometry") else 1, "wall_solves": 2}}


@pytest.fixture
def isolated(tmp_path, monkeypatch):
    root = tmp_path / "synthetic-repository"
    for name in tool.SOURCES:
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(("fixture: "+name).encode())
    monkeypatch.setattr(tool, "ROOT", root)
    monkeypatch.setattr(tool, "runtime_versions", lambda: synthetic()["runtime"])
    monkeypatch.setattr(tool, "compare", synthetic)
    return root / "synthetic-result"


def test_seven_fixed_comparisons_validate_but_scope_domain_and_setting_changes_fail():
    for kind in tool.VARIATIONS:
        assert tool.validate_record(synthetic(kind), kind)
    for path, value in ((["settings", 0, "source_order"], 16),
                        (["solutions", 0, "source_scale"], 0.),
                        (["solutions", 1, "domain", "Ub"], -.26),
                        (["solutions", 0, "scope"], "central observer solved"),
                        (["solutions", 0, "wall_samples", 0, "x"], .44),
                        (["solutions", 0, "wall_samples", 0, "rho_s"], 1.),
                        (["solutions", 0, "radial_evaluations"], 0),
                        (["refinement_is_rigorous_error_enclosure"], True)):
        record = synthetic()
        target = record
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        with pytest.raises(ValueError):
            tool.validate_record(record, "radial")


def test_execution_and_saved_diagnostics_cannot_be_redeclared():
    for key, value in (("cpu_seconds", 121), ("cpu_seconds_max", 121),
                       ("scalar_kernel_values", 30_000_001), ("rhs_calls", 12_001),
                       ("metric_values", -1), ("field_solves", True), ("wall_solves", 1)):
        record = synthetic()
        record["execution"][key] = value
        with pytest.raises(ValueError):
            tool.validate_record(record, "radial")
    record = synthetic()
    record["wall_differences"]["rho"]["maximum_absolute_difference"] = 1.
    with pytest.raises(ValueError, match="differences"):
        tool.validate_record(record, "radial")


def test_real_summary_contract_matches_recorder_without_evolution():
    from recursive_horizons.nsc_atlas_wall_evolution import WallSolution
    record = synthetic()
    for index, saved in enumerate(record["solutions"]):
        solution = object.__new__(WallSolution)
        solution.geometry = saved["domain"]
        solution.settings = {key: saved[key] for key in
                             ("radial_order", "field_order", "time_panels", "source_order", "source_scale", "tolerance")}
        solution.maxima = saved["maximum_absolute_defects"]
        solution.radial_cache = dict.fromkeys(range(saved["radial_evaluations"]))
        solution.budget = SimpleNamespace(rhs_calls=0, metric_values=0, source_nodes=0,
                                         started=time.process_time())
        rows = {row["t"]: row for row in saved["wall_samples"]}
        solution.sample = rows.__getitem__
        summary = solution.summary(tool.TIMES)
        for key in ("rhs_calls", "metric_values", "source_nodes", "process_cpu_seconds"):
            summary.pop(key)
        record["solutions"][index] = summary
    assert tool.validate_record(record, "radial")
    record = synthetic()
    record["solutions"][0]["maximum_absolute_defects"]["mass_path_defect"] = float("nan")
    with pytest.raises(ValueError, match="defect"):
        tool.validate_record(record, "radial")


def test_creation_preserves_existing_bytes_and_check_never_replays(isolated, monkeypatch):
    tool.create(isolated, "radial")
    before = {path.name: path.read_bytes() for path in isolated.iterdir()}
    monkeypatch.setattr(tool, "compare", lambda *_: pytest.fail("unexpected physical calculation"))
    assert tool.check(isolated)["wall_solves"] == 0
    with pytest.raises(FileExistsError, match="creation-only"):
        tool.create(isolated, "radial")
    assert before == {path.name: path.read_bytes() for path in isolated.iterdir()}
    path = tool.ROOT / "src/recursive_horizons/nsc_atlas_causal_observers.py"
    path.write_bytes(path.read_bytes()+b"changed dependency")
    with pytest.raises(ValueError, match="source binding"):
        tool.check(isolated)


@pytest.mark.parametrize("change", ["output-bytes", "rehashed-settings", "rehashed-domain"])
def test_output_tamper_is_rejected_even_after_rehashing(isolated, change):
    tool.create(isolated, "radial")
    path = isolated / "response.json"
    if change == "output-bytes":
        path.write_bytes(path.read_bytes()+b"\n")
    else:
        record = json.loads(path.read_bytes())
        if change == "rehashed-settings":
            record["settings"][0]["field_nodes"] = 513
        else:
            record["solutions"][0]["domain"]["Va"] = -2.4
        raw = tool.encoded(record)
        path.write_bytes(raw)
        manifest_path = isolated / "manifest.json"
        manifest = json.loads(manifest_path.read_bytes())
        manifest["outputs"]["response.json"] = {"sha256": tool.digest(raw), "bytes": len(raw)}
        manifest_path.write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError):
        tool.check(isolated)


def test_mid_comparison_source_change_cannot_create_evidence(isolated, monkeypatch):
    def changed(kind):
        path = tool.ROOT / "src/recursive_horizons/nsc_atlas_wall_evolution.py"
        path.write_bytes(b"changed during solve")
        return synthetic(kind)
    monkeypatch.setattr(tool, "compare", changed)
    with pytest.raises(RuntimeError, match="inputs or runtime changed"):
        tool.create(isolated, "radial")
    assert not isolated.exists()
