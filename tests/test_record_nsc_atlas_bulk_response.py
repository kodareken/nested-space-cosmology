"""Synthetic recorder contracts only: no field or metric evolution."""
from copy import deepcopy
import json

import pytest

import record_nsc_atlas_bulk_response as tool


def synthetic(kind="metric"):
    varied, low, high = tool.VARIATIONS[kind]
    settings = [{**tool.BASE, varied: value} for value in (low, high)]
    solutions = [{"schema": "nsc-atlas-bulk-response-v1", "scope": tool.SCOPE,
                  "domain": {"U0": -.32, "Ua": -.3, "Va": -2.5, "Vb": -2.25},
                  **{key: setting[key] for key in ("metric_nodes", "source_order", "field_order")},
                  "maximum_absolute_defects": dict.fromkeys(tool.DEFECTS, 0.),
                  "transverse_leg_samples": [{"t": t, "mass_defect_U_gauss": 0.,
                                               **dict.fromkeys(tool.QUANTITIES, t)} for t in tool.TIMES]}
                 for setting in settings]
    return {"schema": tool.SCHEMA, "comparison": kind, "varied_quantity": varied,
            "status": tool.STATUS, "order": tool.ORDER, "interpretation": tool.INTERPRETATION,
            "runtime": {"python": "3.14.8", "numpy": "2.5.1", "scipy": "1.17.1"},
            "settings": settings, "solutions": solutions,
            "transverse_leg_differences": tool.leg_differences(solutions),
            "common_grid_maximum_absolute_differences": dict.fromkeys(tool.QUANTITIES, 0.),
            "refinement_is_rigorous_error_enclosure": False,
            "execution": {"cpu_seconds": 1., "cpu_seconds_max": 120.,
                          "scalar_kernel_values": 100, "scalar_kernel_values_max": 30_000_000,
                          "metric_values": 100, "source_quadrature_nodes": 100, "metric_rhs_calls": 10,
                          "field_solves": 2 if kind in ("field", "geometry") else 1, "metric_solves": 2}}


def test_settings_and_scientific_scope_are_checked_without_replay():
    for kind in tool.VARIATIONS:
        assert tool.validate_record(synthetic(kind), kind)
    for path, value in ((["settings", 0, "source_order"], 16),
                        (["status"], "FULL_WALL_SOLUTION"),
                        (["refinement_is_rigorous_error_enclosure"], True),
                        (["solutions", 1, "domain", "Ua"], -.29),
                        (["solutions", 0, "transverse_leg_samples", 0, "h_J"], float("nan"))):
        record = synthetic()
        target = record
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        with pytest.raises(ValueError):
            tool.validate_record(record, "metric")


def test_budget_and_saved_difference_contracts():
    for key, value in (("cpu_seconds", 121), ("cpu_seconds_max", 121),
                       ("scalar_kernel_values", 30_000_001), ("metric_rhs_calls", 12_001),
                       ("source_quadrature_nodes", -1), ("field_solves", True)):
        record = synthetic()
        record["execution"][key] = value
        with pytest.raises(ValueError):
            tool.validate_record(record, "metric")
    record = synthetic()
    record["transverse_leg_differences"]["h_J"]["maximum_absolute_difference"] = 1.
    with pytest.raises(ValueError, match="differences"):
        tool.validate_record(record, "metric")


def test_creation_preservation_source_drift_and_rehashed_tamper(tmp_path, monkeypatch):
    assert {"pyproject.toml", "src/recursive_horizons/nsc_atlas_causal_observers.py",
            "src/recursive_horizons/nsc_atlas_connection_slice.py",
            "src/recursive_horizons/nsc_atlas_moving_boundaries.py"} <= set(tool.SOURCES)
    bindings = {"fixture-source": "original"}
    monkeypatch.setattr(tool, "source_bindings", lambda: dict(bindings))
    monkeypatch.setattr(tool, "runtime_versions", lambda: synthetic()["runtime"])
    monkeypatch.setattr(tool, "compare", synthetic)
    directory = tmp_path / "synthetic-comparison"
    tool.create(directory, "metric")
    monkeypatch.setattr(tool, "compare", lambda *_: pytest.fail("unexpected physical replay"))
    assert tool.check(directory)["metric_solves"] == 0
    before = {path.name: path.read_bytes() for path in directory.iterdir()}
    with pytest.raises(FileExistsError, match="creation-only"):
        tool.create(directory, "metric")
    assert before == {path.name: path.read_bytes() for path in directory.iterdir()}
    bindings["fixture-source"] = "changed"
    with pytest.raises(ValueError, match="source binding"):
        tool.check(directory)
    bindings["fixture-source"] = "original"
    record = deepcopy(json.loads(before["response.json"]))
    record["settings"][0]["field_nodes"] = 513
    raw = tool.encoded(record)
    (directory / "response.json").write_bytes(raw)
    manifest = json.loads(before["manifest.json"])
    manifest["outputs"]["response.json"] = {"sha256": tool.digest(raw), "bytes": len(raw)}
    (directory / "manifest.json").write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError, match="settings"):
        tool.check(directory)
