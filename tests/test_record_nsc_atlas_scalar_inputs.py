"""Synthetic input records and real polynomial checks; no physical solve."""
from copy import deepcopy
import json
import math
from types import SimpleNamespace

import numpy as np
import pytest

import record_nsc_atlas_scalar_inputs as tool


DOMAIN = dict(U0=-.32, Ua=-.30, Ub=-.27, Va=-2.5, Vb=-2.25)


def lift_fixture(monkeypatch):
    from recursive_horizons import nsc_atlas_scalar_boundary as boundary
    from recursive_horizons import nsc_atlas_scalar_lift as lift
    class Jet:
        def __init__(self, _):
            self.domain, self.settings = dict(DOMAIN), dict(radial_order=16, ray_order=16)
            self.budget = SimpleNamespace(check=lambda: None, metric_values=0, started=0.)
            self.cap_field = SimpleNamespace(budget=SimpleNamespace(kernel_values=0))
            self.geometry = {"wall": lambda V: DOMAIN["Ua"]+(DOMAIN["Ub"]-DOMAIN["Ua"])*(np.asarray(V)-DOMAIN["Va"])/(DOMAIN["Vb"]-DOMAIN["Va"])}
            self.cache, self.ray_nodes = {}, 0
        def sample(self, t):
            if t not in self.cache:
                dU, dV = DOMAIN["Ub"]-DOMAIN["Ua"], DOMAIN["Vb"]-DOMAIN["Va"]
                row = {**dict.fromkeys(tool.JET_KEYS, 0.), "t": t, "U": DOMAIN["Ua"]+dU*t,
                    "V": DOMAIN["Va"]+dV*t, "x": .45+.1*t, "b": t**4,
                    "b_U": 4*t**3/dU, "b_V": 4*t**3/dV, "U_V": dU/dV,
                    "r_U": 1., "c1_r": t**4, "J_minus": t**4, "positive_integral": t**4,
                    "J_minus_u": 4*t**3, "integral_Ku_f": -4*t**3, "q_minus_s": 1.,
                    "J_minus_s": 4*t**3, "ray_nodes": 16 if t else 0,
                    **dict.fromkeys(tool.MARGINS, 1.)}
                self.cache[t] = row
                self.ray_nodes += row["ray_nodes"]
            return dict(self.cache[t])
    with monkeypatch.context() as patch:
        patch.setattr(boundary, "ScalarBoundary", lambda *args, **kwargs: None)
        patch.setattr(lift, "BoundaryJet", Jet)
        return tool._compare_lift(None, None, None)


def forcing_fixture(endpoint_roundoff=False):
    grid = np.linspace(DOMAIN["U0"], DOMAIN["Ub"], 129)
    solutions = []
    for order in (2, 4):
        rows = []
        for t in tool.TIMES:
            wall = DOMAIN["Ua"]+(DOMAIN["Ub"]-DOMAIN["Ua"])*t
            if t == 1 and endpoint_roundoff:
                wall = np.nextafter(DOMAIN["Ub"], 0.)
            U = np.minimum(grid, wall)
            U[-1] = wall
            cuts = np.unique(np.r_[DOMAIN["U0"], grid[(grid > DOMAIN["U0"]) & (grid < wall)], U, wall])
            gauss = cuts[:-1, None]+(np.diff(cuts)/2)[:, None]*(1+np.polynomial.legendre.leggauss(order)[0])
            points = len(np.unique(np.r_[cuts, gauss.ravel()]))
            mu = -t*(U-DOMAIN["U0"])*(1+order*1e-7)
            integral = -t*(wall-DOMAIN["U0"])
            row = {**dict.fromkeys(tool.FORCING_SCALARS, 0.),
                "t": t, "V": DOMAIN["Va"]+(DOMAIN["Vb"]-DOMAIN["Va"])*t, "U_wall": float(wall),
                "U": U.tolist(), "K0": np.ones(129).tolist(), "F": (t*(wall-U)).tolist(),
                "G": np.full(129, t).tolist(), "mu_E": mu.tolist(), "zeta_E": (t*(U-DOMAIN["U0"])).tolist(),
                "g1": -t*t, "g1_V": -(1+t)*t+t*t, "source_H_E": t, "source_H_E_V": 1.,
                "source_chi_r": t, "source_chi": t, "source_r_U": 1., "source_P0_r": 2.,
                "wall_integral_L_plus": float(integral), "wall_mass_defect": float(mu[-1]-integral),
                "radial_panels": len(cuts)-1, "radial_nodes": gauss.size,
                "field_points": points, "field_batches": (points+512)//513}
            rows.append(row)
        solutions.append({"schema": "nsc-atlas-scalar-forcing-v1",
            "settings": dict(radial_order=order, field_order=32, source_order=64, source_scale=1., grid_nodes=129),
            "samples": rows, "slice_evaluations": 9,
            "field_points": sum(row["field_points"] for row in rows),
            "field_batches": sum(row["field_batches"] for row in rows),
            "maximum_absolute_defects": {key: max(abs(row[key]) for row in rows)
                for key in ("wall_mass_defect", "source_p_defect", "source_p_V_defect", "wall_chi_U", "chi_U_lower_limit_bound")}})
    return solutions, {"requested_U_grid": grid.tolist(), "maximum_absolute_differences": tool._forcing_differences(solutions)}


def synthetic(kind, monkeypatch, endpoint_roundoff=False):
    solutions, diagnostics = lift_fixture(monkeypatch) if kind == "lift" else forcing_fixture(endpoint_roundoff)
    varied, low, high = tool.VARIATIONS[kind]
    return {"schema": tool.SCHEMA, "comparison": kind, "status": tool.STATUS, "scope": tool.SCOPE,
        "source_scale": 1., "runtime": {"python": "3.14.8", "numpy": "2.5.1", "scipy": "1.17.1"},
        "comparison_base": dict(tool.BASE), "varied_quantity": varied,
        "settings": [{**tool.BASE, varied: value} for value in (low, high)], "domain": dict(DOMAIN),
        "solutions": solutions, "diagnostics": diagnostics, "refinement_is_rigorous_error_enclosure": False,
        "profile_evaluation": {"counts": dict.fromkeys(tool.PROFILE_COUNTS, 0), "maximum_inverse_bracket_width": 0.},
        "execution": dict(cpu_seconds=1., cpu_seconds_max=120., scalar_kernel_values=100,
            scalar_kernel_values_max=30_000_000, metric_values=100, source_quadrature_nodes=0,
            rhs_calls=10, field_solves=1, wall_solves=1)}


@pytest.fixture
def isolated(tmp_path, monkeypatch):
    root = tmp_path/"synthetic-repository"
    for name in tool.SOURCES:
        path = root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes((tool.ROOT/name).read_bytes() if name in (tool.TARGET_RECORD, tool.TARGET_MANIFEST)
                         else ("fixture: "+name).encode())
    monkeypatch.setattr(tool, "ROOT", root)
    monkeypatch.setattr(tool, "runtime_versions", lambda: {"python": "3.14.8", "numpy": "2.5.1", "scipy": "1.17.1"})
    monkeypatch.setattr(tool, "compare", lambda kind: synthetic(kind, monkeypatch))
    return root/"result"


@pytest.mark.parametrize("kind", ["lift", "forcing"])
def test_fixed_inputs_and_saved_algebra_validate_without_evolution(monkeypatch, kind):
    record = synthetic(kind, monkeypatch)
    assert tool.validate_record(record, kind)
    for path, value in ((["settings", 0, "field_nodes"], 257), (["source_scale"], True),
                        (["execution", "scalar_kernel_values"], 30_000_001),
                        (["execution", "cpu_seconds_max"], 121), (["execution", "field_solves"], True),
                        (["refinement_is_rigorous_error_enclosure"], True)):
        changed = deepcopy(record)
        target = changed
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        with pytest.raises(ValueError):
            tool.validate_record(changed, kind)


def test_lift_withheld_target_is_diagnostic_and_tampering_is_rejected(monkeypatch):
    record = synthetic("lift", monkeypatch)
    assert record["diagnostics"]["sampled_b_target_met"] == [False, False]
    assert tool.validate_record(record, "lift")
    changes = ((["solutions", 0, "validation_samples", 3, "B"], 9.),
               (["solutions", 0, "knot_samples", 2, "b_U"], 8.),
               (["diagnostics", "direct_samples", 8, "J_minus_s"], 7.),
               (["diagnostics", "unique_boundary_jet_evaluations"], 34),
               (["diagnostics", "interpolation_target", "absolute_b_target"], 1.))
    for path, value in changes:
        changed = deepcopy(record)
        target = changed
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        with pytest.raises(ValueError):
            tool.validate_record(changed, "lift")


def test_target_rebuild_roundoff_does_not_relax_reference_metadata(monkeypatch):
    record = synthetic("lift", monkeypatch)
    target = deepcopy(record["diagnostics"]["interpolation_target"])
    for key in ("sampled_field_b_sensitivity", "absolute_b_target"):
        target[key] = float(np.nextafter(target[key], np.inf))
    monkeypatch.setattr(tool, "interpolation_target", lambda: deepcopy(target))
    assert tool.validate_record(record, "lift")
    for key, value in (("absolute_b_target", 1.), ("source_sha256", "changed"),
                       ("reference_cuts", 9.), ("is_uniform_bound", 0)):
        changed = deepcopy(record)
        changed["diagnostics"]["interpolation_target"][key] = value
        with pytest.raises(ValueError, match="target reference"):
            tool.validate_record(changed, "lift")


def test_forcing_accepts_inherited_wall_ulp_but_not_changed_inputs(monkeypatch):
    record = synthetic("forcing", monkeypatch, endpoint_roundoff=True)
    assert tool.validate_record(record, "forcing")
    for key, value in (("g1_V", 20.), ("wall_mass_defect", 3.), ("field_points", 1), ("wall_chi_U", 1.)):
        changed = deepcopy(record)
        changed["solutions"][0]["samples"][4][key] = value
        with pytest.raises(ValueError):
            tool.validate_record(changed, "forcing")


@pytest.mark.parametrize("kind", ["lift", "forcing"])
def test_creation_only_and_authentication_never_solve(isolated, monkeypatch, kind):
    created = tool.create(isolated, kind)
    assert "diagnostics" not in created and "summary" in created
    before = {p.name: p.read_bytes() for p in isolated.iterdir()}
    monkeypatch.setattr(tool, "compare", lambda *_: pytest.fail("unexpected physical solve"))
    assert tool.check(isolated)["physical_solves"] == 0
    with pytest.raises(FileExistsError):
        tool.create(isolated, kind)
    assert before == {p.name: p.read_bytes() for p in isolated.iterdir()}
    source = tool.ROOT/"src/recursive_horizons/nsc_atlas_scalar_forcing.py"
    source.write_bytes(b"source drift")
    with pytest.raises(ValueError, match="source binding"):
        tool.check(isolated)


@pytest.mark.parametrize("kind", ["lift", "forcing"])
def test_rehashed_physical_identity_tamper_fails(isolated, kind):
    tool.create(isolated, kind)
    path = isolated/"inputs.json"
    record = json.loads(path.read_bytes())
    if kind == "lift":
        record["solutions"][0]["validation_samples"][4]["b_interpolated"] = 50.
    else:
        record["solutions"][0]["samples"][4]["g1_V"] = 50.
    raw = tool.encoded(record)
    path.write_bytes(raw)
    manifest_path = isolated/"manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    manifest["outputs"]["inputs.json"] = {"sha256": tool.digest(raw), "bytes": len(raw)}
    manifest_path.write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError):
        tool.check(isolated)


def test_source_drift_during_comparison_writes_nothing(isolated, monkeypatch):
    def changed(kind):
        (tool.ROOT/"docs/nsc-atlas-scalar-forcing.md").write_bytes(b"changed during calculation")
        return synthetic(kind, monkeypatch)
    monkeypatch.setattr(tool, "compare", changed)
    with pytest.raises(RuntimeError, match="inputs or runtime changed"):
        tool.create(isolated, "forcing")
    assert not isolated.exists()
