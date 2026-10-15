"""Synthetic moment evidence and authentication, without physical solves."""
from copy import deepcopy
import json
import math

import pytest

import record_nsc_atlas_profile_moment as tool


def synthetic():
    domain = dict(U0=-.32, Ua=-.30, Ub=-.27, Va=-2.5, Vb=-2.25,
                  qa=1., qb=1.5, xa=.45, xb=.55)
    domain.update(Ta=1.-math.atanh(.45), Tb=1.5-math.atanh(.55), qmin=1.-2*math.atanh(.45))
    rows = []
    for t in tool.TIMES:
        x, q = .45+.1*t, 1.+.5*t
        u = q-2*math.atanh(x)
        rows.append({**dict.fromkeys(tool.SAMPLE_KEYS, 0.), "t": t, "U": -.30+.03*t, "V": -2.5+.25*t,
            "x": x, "q_plus": q, "q_minus": u, "b": t*t, "c1_r": t*t,
            "b_U": 2*t/.03, "b_V": 2*t/.25, "U_V": .03/.25, "r_U": 1.,
            "J_minus": t*t, "positive_integral": t*t, "J_minus_s": 2*t, "J_minus_u": 2*t,
            "integral_Ku_f": -2*t, "q_minus_s": 1., "ray_nodes": 16 if t else 0,
            "no_return_margin": 1.-u, "quiet_centre_time_margin": 1.-(u+q)/2,
            "initial_time_margin": (u+1.)/2-domain["Ta"], "last_profile_margin": 1.5-q,
            "first_profile_margin": u-domain["qmin"], "active_radius_lower": math.tanh((1.-u)/2)})
    masses = []
    for row in rows[1:]:
        for fraction in tool.FRACTIONS:
            q = 1.+fraction*(row["q_plus"]-1.)
            u = row["q_minus"]
            masses.append(dict(t=row["t"], fraction=fraction, q_minus=u, q=q,
                T=(u+q)/2, r=math.tanh((q-u)/2), nu_E=row["t"]*fraction))
    reference = {"method": "original radial metric and joint boundary jet", "radial_order": 16, "ray_order": 16,
        "boundary_samples": rows, "boundary_evaluations": 9, "ray_nodes": 128,
        "mass_samples": masses, "mass_evaluations": 24}
    solutions = []
    for order in tool.ORDERS:
        copied = deepcopy(masses)
        for point in copied:
            point["f"] = .1*point["t"]
            point["nu_E"] += order*1e-7
            point["C"] = point["nu_E"]+(1-point["r"]**2)*point["f"]**2/(4*point["r"])
        solutions.append({"moment_order": order, "boundary_samples": deepcopy(rows),
            "boundary_evaluations": 9, "ray_nodes": 128, "mass_samples": copied, "mass_evaluations": 24,
            "moment": {"schema": "nsc-atlas-profile-moment-v1", "domain": dict(domain),
                "settings": {"moment_order": order}, "cells": 256, "prefix_nodes": 256*order,
                "partial_nodes": 0, "profile_points": 256*order+161,
                "profile_batches": 161+(256*order+512)//513, "value_queries": 161,
                "cache_hits": 20, "cached_values": 50, "mass_points": 161, "terminal_C": 2.}})
    return {"schema": tool.SCHEMA, "status": tool.STATUS, "scope": tool.SCOPE,
        "runtime": {"python": "3.14.8", "numpy": "2.5.1", "scipy": "1.17.1"},
        "comparison_base": dict(tool.BASE), "varied_quantity": "moment_order",
        "settings": [{**tool.BASE, "moment_order": order} for order in tool.ORDERS],
        "source_scale": 1., "domain": domain, "reference": reference, "solutions": solutions,
        "diagnostics": tool.diagnostics(reference, solutions), "refinement_is_rigorous_error_enclosure": False,
        "profile_evaluation": {"counts": dict(profile_requested=1000, wall_points=900, second_points=0, cache_hits=100),
                               "maximum_inverse_bracket_width": 0.},
        "execution": dict(cpu_seconds=1., cpu_seconds_max=120., scalar_kernel_values=10000,
            scalar_kernel_values_max=30_000_000, metric_values=1000, source_quadrature_nodes=100,
            rhs_calls=100, field_solves=1, wall_solves=1, moment_evaluators=2)}


@pytest.fixture
def isolated(tmp_path, monkeypatch):
    root = tmp_path/"repository"
    for name in tool.SOURCES:
        path = root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(("fixture: "+name).encode())
    monkeypatch.setattr(tool, "ROOT", root)
    monkeypatch.setattr(tool, "runtime_versions", lambda: synthetic()["runtime"])
    monkeypatch.setattr(tool, "compare", synthetic)
    return root/"result"


def test_reference_mass_and_boundary_identities_and_settings():
    record = synthetic()
    assert tool.validate_record(record)
    for path, value in ((["settings", 0, "field_nodes"], 257),
                       (["reference", "radial_order"], 8),
                       (["solutions", 0, "mass_samples", 2, "C"], 30.),
                       (["solutions", 1, "mass_samples", 4, "fraction"], .4),
                       (["solutions", 0, "boundary_samples", 5, "J_minus_s"], 10.),
                       (["solutions", 0, "boundary_samples", 5, "no_return_margin"], -1.),
                       (["solutions", 1, "moment", "prefix_nodes"], 2),
                       (["source_scale"], True), (["refinement_is_rigorous_error_enclosure"], True)):
        changed = deepcopy(record)
        target = changed
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        with pytest.raises(ValueError):
            tool.validate_record(changed)


def test_cumulative_budgets_and_no_second_derivative_claim():
    for key, value in (("cpu_seconds", 121), ("cpu_seconds_max", 121),
                       ("scalar_kernel_values", 30_000_001), ("rhs_calls", 12001),
                       ("field_solves", True), ("moment_evaluators", 1)):
        record = synthetic()
        record["execution"][key] = value
        with pytest.raises(ValueError):
            tool.validate_record(record)
    record = synthetic()
    record["profile_evaluation"]["counts"]["second_points"] = 1
    with pytest.raises(ValueError, match="second derivatives"):
        tool.validate_record(record)
    record = synthetic()
    record["diagnostics"]["mass_order_difference"] = 0.
    with pytest.raises(ValueError, match="saved differences"):
        tool.validate_record(record)


def test_actual_moment_summary_contract_without_profile_queries():
    from recursive_horizons.nsc_atlas_profile_moment import ProfileMoment
    from types import SimpleNamespace
    import numpy as np
    record = synthetic()
    for solution in record["solutions"]:
        saved = solution["moment"]
        moment = object.__new__(ProfileMoment)
        moment.budget = SimpleNamespace(check=lambda: None)
        moment.domain, moment.settings = saved["domain"], saved["settings"]
        moment.cuts = np.linspace(1., 1.5, saved["cells"]+1)
        moment.prefix = np.linspace(0., saved["terminal_C"], saved["cells"]+1)
        moment.cache = dict.fromkeys(range(saved["cached_values"]))
        for key in ("prefix_nodes", "partial_nodes", "profile_points", "profile_batches", "value_queries", "cache_hits", "mass_points"):
            setattr(moment, key, saved[key])
        solution["moment"] = moment.summary()
    assert tool.validate_record(record)


def test_creation_only_and_authentication_perform_no_solve(isolated, monkeypatch):
    tool.create(isolated)
    before = {p.name: p.read_bytes() for p in isolated.iterdir()}
    monkeypatch.setattr(tool, "compare", lambda: pytest.fail("unexpected physical solve"))
    assert tool.check(isolated)["physical_solves"] == 0
    with pytest.raises(FileExistsError):
        tool.create(isolated)
    assert before == {p.name: p.read_bytes() for p in isolated.iterdir()}
    (tool.ROOT/"src/recursive_horizons/nsc_atlas_profile_moment.py").write_bytes(b"source drift")
    with pytest.raises(ValueError, match="source binding"):
        tool.check(isolated)


@pytest.mark.parametrize("change", ["bytes", "mass", "boundary"])
def test_rehashed_tampering_fails(isolated, change):
    tool.create(isolated)
    path = isolated/"moment.json"
    if change == "bytes":
        path.write_bytes(path.read_bytes()+b"\n")
    else:
        record = json.loads(path.read_bytes())
        if change == "mass":
            record["solutions"][0]["mass_samples"][3]["nu_E"] = 9.
        else:
            record["solutions"][1]["boundary_samples"][4]["b_U"] = 9.
        raw = tool.encoded(record)
        path.write_bytes(raw)
        manifest_path = isolated/"manifest.json"
        manifest = json.loads(manifest_path.read_bytes())
        manifest["outputs"]["moment.json"] = {"sha256": tool.digest(raw), "bytes": len(raw)}
        manifest_path.write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError):
        tool.check(isolated)


def test_mid_comparison_source_drift_cannot_write_a_record(isolated, monkeypatch):
    def changed():
        (tool.ROOT/"docs/nsc-atlas-profile-moment.md").write_bytes(b"changed during comparison")
        return synthetic()
    monkeypatch.setattr(tool, "compare", changed)
    with pytest.raises(RuntimeError, match="inputs or runtime changed"):
        tool.create(isolated)
    assert not isolated.exists()
