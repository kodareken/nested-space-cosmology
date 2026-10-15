"""Saved polynomial/provenance controls; no physical preparation or evolution."""
from copy import deepcopy
import json
import math

import pytest

import record_nsc_atlas_local_lift as tool


def synthetic(monkeypatch):
    domain = dict(U0=-.32, Ua=-.30, Ub=-.27, Va=-2.5, Vb=-2.25,
                  qa=1., qb=1.5, xa=.45, xb=.55)
    domain.update(Ta=1.-math.atanh(.45), Tb=1.5-math.atanh(.55), qmin=1.-2*math.atanh(.45))
    inputs = {"path": tool.INPUT, "files": {"moment.json": "a"*64, "manifest.json": "b"*64},
              "domain": domain, "physical_solves": 0}
    monkeypatch.setattr(tool, "input_snapshot", lambda: deepcopy(inputs))
    rows = []
    for t in tool.DIRECT_TIMES:
        x, q = .45+.1*t, 1.+.5*t
        u = q-2*math.atanh(x)
        rows.append({**dict.fromkeys(tool.predecessor.SAMPLE_KEYS, 0.), "t": t, "U": -.30+.03*t, "V": -2.5+.25*t,
            "x": x, "q_plus": q, "q_minus": u, "b": t**10, "c1_r": t**10,
            "b_U": 10*t**9/.03, "b_V": 10*t**9/.25, "U_V": .03/.25, "r_U": 1.,
            "J_minus": t*t, "positive_integral": t*t, "J_minus_s": 2*t, "J_minus_u": 2*t,
            "integral_Ku_f": -2*t, "q_minus_s": 1., "ray_nodes": 16 if t else 0,
            "no_return_margin": 1.-u, "quiet_centre_time_margin": 1.-(u+q)/2,
            "initial_time_margin": (u+1.)/2-domain["Ta"], "last_profile_margin": 1.5-q,
            "first_profile_margin": u-domain["qmin"], "active_radius_lower": math.tanh((1.-u)/2)})
    record = {"schema": tool.SCHEMA, "status": tool.STATUS, "scope": tool.SCOPE,
        "runtime": {"python": "3.14.8", "numpy": "2.5.1", "scipy": "1.17.1"},
        "input": inputs, "comparison_base": dict(tool.SETTINGS), "varied_quantity": "active_cells",
        "settings": [{**tool.SETTINGS, "active_cells": count} for count in tool.ACTIVE_CELLS],
        "source_scale": 1., "domain": domain, "direct_samples": rows,
        "refinement_is_rigorous_error_enclosure": False, "unique_boundary_jet_evaluations": 81, "ray_nodes": 1280,
        "moment": {"schema": "nsc-atlas-profile-moment-v1", "domain": dict(domain), "settings": {"moment_order": 8},
            "cells": 256, "prefix_nodes": 2048, "partial_nodes": 0, "profile_points": 3409,
            "profile_batches": 1365, "value_queries": 1361, "cache_hits": 100, "cached_values": 1000,
            "mass_points": 1361, "terminal_C": 2.},
        "maximum_absolute_jet_defects": {key: max(abs(row[key]) for row in rows) for key in tool.JET_DEFECTS},
        "minimum_domain_margins": {key: min(row[key] for row in rows) for key in tool.MARGINS},
        "profile_evaluation": {"counts": dict(profile_requested=10000, wall_points=9000, second_points=0, cache_hits=100),
                               "maximum_inverse_bracket_width": 0.},
        "execution": dict(cpu_seconds=1., cpu_seconds_max=120., scalar_kernel_values=10000,
            scalar_kernel_values_max=30_000_000, metric_values=1000, source_quadrature_nodes=100,
            rhs_calls=100, field_solves=1, wall_solves=1, moment_evaluators=1)}
    record["solutions"] = tool._saved_fits(record)
    record["diagnostics"] = tool.comparison_diagnostics(record["solutions"])
    return record


@pytest.fixture
def isolated(tmp_path, monkeypatch):
    root = tmp_path/"repository"
    for name in tool.SOURCES:
        path = root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(("fixture: "+name).encode())
    monkeypatch.setattr(tool, "ROOT", root)
    monkeypatch.setattr(tool, "runtime_versions", lambda: {"python": "3.14.8", "numpy": "2.5.1", "scipy": "1.17.1"})
    monkeypatch.setattr(tool, "compare", lambda: synthetic(monkeypatch))
    return root/"result"


def test_nested_stencils_and_heldout_polynomials_reconstruct(monkeypatch):
    record = synthetic(monkeypatch)
    assert [len(tool.knot_times(n)) for n in tool.ACTIVE_CELLS] == [25, 41]
    assert len(tool.WITHHELD) == 40 and len(tool.DIRECT_TIMES) == 81
    assert tool.validate_record(record)
    assert all(solution["heldout_used_for_training"] is False for solution in record["solutions"])
    assert all(solution["extrapolated_queries"] == 0 for solution in record["solutions"])
    for path, value in ((["settings", 0, "active_cells"], 20),
                       (["solutions", 0, "cells", 5, "training_indices"], [1, 2, 3, 4, 5]),
                       (["solutions", 0, "validation_samples", 9, "B_lift"], 30.),
                       (["solutions", 1, "maximum_absolute_endpoint_value_defect"], 20.),
                       (["direct_samples", 25, "b_U"], 20.),
                       (["diagnostics", "solutions", 0, "per_cell", 3, "b_defect"], 1.),
                       (["refinement_is_rigorous_error_enclosure"], True)):
        changed = deepcopy(record)
        target = changed
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        with pytest.raises(ValueError):
            tool.validate_record(changed)


def test_budget_counters_predecessor_and_target_are_not_redeclarable(monkeypatch):
    record = synthetic(monkeypatch)
    for key, value in (("cpu_seconds", 121), ("cpu_seconds_max", 121),
                       ("scalar_kernel_values", 30_000_001), ("field_solves", True), ("moment_evaluators", 2)):
        changed = deepcopy(record)
        changed["execution"][key] = value
        with pytest.raises(ValueError):
            tool.validate_record(changed)
    changed = deepcopy(record)
    changed["input"]["files"]["moment.json"] = "c"*64
    with pytest.raises(ValueError, match="predecessor"):
        tool.validate_record(changed)
    changed = deepcopy(record)
    changed["diagnostics"]["interpolation_target"]["absolute_b_target"] = 1.
    with pytest.raises(ValueError, match="target"):
        tool.validate_record(changed)


def test_creation_only_authentication_rebuilds_only_saved_polynomials(isolated, monkeypatch):
    tool.create(isolated)
    before = {p.name: p.read_bytes() for p in isolated.iterdir()}
    monkeypatch.setattr(tool, "compare", lambda: pytest.fail("unexpected physical comparison"))
    monkeypatch.setattr(tool, "_prepare", lambda: pytest.fail("unexpected field preparation"))
    checked = tool.check(isolated)
    assert checked["physical_solves"] == 0 and checked["saved_polynomial_reconstructions"] == 2
    with pytest.raises(FileExistsError):
        tool.create(isolated)
    assert before == {p.name: p.read_bytes() for p in isolated.iterdir()}
    (tool.ROOT/"src/recursive_horizons/nsc_atlas_adaptive_lift.py").write_bytes(b"source drift")
    with pytest.raises(ValueError, match="source binding"):
        tool.check(isolated)


@pytest.mark.parametrize("change", ["bytes", "polynomial", "heldout"])
def test_rehashed_tampering_fails(isolated, change):
    tool.create(isolated)
    path = isolated/"lift.json"
    if change == "bytes":
        path.write_bytes(path.read_bytes()+b"\n")
    else:
        record = json.loads(path.read_bytes())
        if change == "polynomial":
            record["solutions"][0]["cells"][4]["coefficients_in_cell_coordinate"][0] += 1.
        else:
            record["solutions"][1]["validation_samples"][4]["b_interpolated"] += 1.
        raw = tool.encoded(record)
        path.write_bytes(raw)
        manifest_path = isolated/"manifest.json"
        manifest = json.loads(manifest_path.read_bytes())
        manifest["outputs"]["lift.json"] = {"sha256": tool.digest(raw), "bytes": len(raw)}
        manifest_path.write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError):
        tool.check(isolated)


def test_mid_comparison_source_drift_writes_nothing(isolated, monkeypatch):
    def changed():
        (tool.ROOT/"src/recursive_horizons/nsc_atlas_adaptive_lift.py").write_bytes(b"changed during comparison")
        return synthetic(monkeypatch)
    monkeypatch.setattr(tool, "compare", changed)
    with pytest.raises(RuntimeError, match="inputs or runtime changed"):
        tool.create(isolated)
    assert not isolated.exists()


def test_input_snapshot_authenticates_both_original_files(tmp_path, monkeypatch):
    root = tmp_path
    path = root/tool.INPUT
    path.mkdir(parents=True)
    (path/"moment.json").write_bytes(tool.encoded({"domain": {"test": 1.}}))
    (path/"manifest.json").write_bytes(b"manifest fixture")
    calls = []
    monkeypatch.setattr(tool, "ROOT", root)
    monkeypatch.setattr(tool.predecessor, "check", lambda selected: calls.append(selected))
    result = tool.input_snapshot()
    assert calls == [path]
    assert result["files"] == {name: tool.digest((path/name).read_bytes()) for name in ("moment.json", "manifest.json")}
    assert result["physical_solves"] == 0
