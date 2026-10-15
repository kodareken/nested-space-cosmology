"""Authentication and saved interpolation controls, without physical solves."""
from copy import deepcopy
import json

import numpy as np
import pytest

import record_nsc_atlas_forcing_cache as tool


def synthetic():
    from recursive_horizons.nsc_atlas_characteristic_initialization import source_geometry
    geometry = source_geometry()
    d = {key: geometry[key] for key in tool.DOMAIN_KEYS}
    native = np.linspace(d["U0"], d["Ub"], 129)
    rows = []
    for i, t in enumerate(tool.DIRECT_TIMES):
        eta = np.linspace(0, 1, 129 if i % 2 == 0 else 257)
        wall = d["Ua"]+(d["Ub"]-d["Ua"])*t
        if t == 1:
            wall = np.nextafter(d["Ub"], 0.)
        U = d["U0"]+(wall-d["U0"])*eta
        U[-1] = wall
        h, hV = max(0., t-.25)**4, 4*max(0., t-.25)**3/(d["Vb"]-d["Va"])
        mu = -h*(U-d["U0"])
        cuts = np.unique(np.r_[d["U0"], native[(native > d["U0"]) & (native < wall)], U])
        gauss = cuts[:-1, None]+(np.diff(cuts)/2)[:, None]*(1+np.polynomial.legendre.leggauss(2)[0])
        points = len(np.unique(np.r_[cuts, gauss.ravel()]))
        rows.append({**dict.fromkeys(tool.FORCING_SCALARS, 0.),
            "t": t, "V": d["Va"]+(d["Vb"]-d["Va"])*t, "U_wall": float(wall),
            "eta": eta.tolist(), "U": U.tolist(), "K0": np.ones(eta.size).tolist(),
            "F": (h*(1-eta)**2).tolist(), "G": (h*(1+.5*eta**2)).tolist(),
            "mu_E": mu.tolist(), "zeta_E": (-mu).tolist(),
            "g1": -h, "g1_V": -hV, "source_H_E": h, "source_H_E_V": hV,
            "source_chi_r": 1., "source_r_U": 1.,
            "wall_integral_L_plus": float(mu[-1]), "wall_mass_defect": 0.,
            "radial_panels": len(cuts)-1, "radial_nodes": gauss.size,
            "field_points": points, "field_batches": (points+512)//513,
            "slice_costs": {"kernel_values": 1000, "metric_values": points, "source_nodes": 64}})
    provider = {key: sum(row[key] for row in rows) for key in ("radial_nodes", "field_points", "field_batches")}
    provider.update({key: sum(row["slice_costs"][key] for row in rows)
                     for key in ("kernel_values", "metric_values", "source_nodes")})
    provider.update(physical_slices=65, cache_hits=49, maximum_absolute_input_defects={
        "quiet_input_defect": 0., "wall_flux_defect": 0., "source_p_defect": 0., "source_p_V_defect": 0.})
    record = {"schema": tool.SCHEMA, "status": tool.STATUS, "scope": tool.SCOPE, "source_scale": 1.,
        "runtime": {"python": "3.14.8", "numpy": "2.5.1", "scipy": "1.17.1"},
        "settings": dict(tool.SETTINGS), "domain": d, "source_geometry": geometry,
        "direct_samples": rows, "slice_provider": provider,
        "refinement_is_rigorous_error_enclosure": False, "heldout_used_for_training": False,
        "profile_evaluation": {"counts": dict.fromkeys(tool.PROFILE_COUNTS, 0), "maximum_inverse_bracket_width": 0.},
        "execution": {"cpu_seconds": 1., "cpu_seconds_max": 120., "scalar_kernel_values": 75000,
            "scalar_kernel_values_max": 30_000_000, "metric_values": provider["metric_values"]+100,
            "source_quadrature_nodes": provider["source_nodes"]+100, "rhs_calls": 10,
            "field_solves": 1, "wall_solves": 1, "direct_forcing_slices": 65}}
    record["direct_diagnostics"] = tool._validate_direct(record)
    record["solutions"] = tool._saved_solutions(record)
    record["diagnostics"] = tool.comparison_diagnostics(record["solutions"])
    return record


@pytest.fixture(scope="module")
def record():
    return synthetic()


@pytest.fixture
def isolated(tmp_path, monkeypatch, record):
    root = tmp_path/"repository"
    for name in tool.SOURCES:
        path = root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(("fixture: "+name).encode())
    monkeypatch.setattr(tool, "ROOT", root)
    monkeypatch.setattr(tool, "runtime_versions", lambda: dict(record["runtime"]))
    monkeypatch.setattr(tool, "compare", lambda: deepcopy(record))
    return root/"result"


def test_shared_nested_time_fits_keep_radial_and_temporal_errors_distinct(record):
    assert tool.validate_record(record)
    assert [row["settings"]["training_times"] for row in record["solutions"]] == [17, 33]
    assert record["slice_provider"]["physical_slices"] == 65
    coarse, fine = record["diagnostics"]["maximum_absolute_defects"]
    assert fine["temporal_F_defect"] < coarse["temporal_F_defect"]
    assert fine["combined_G_defect"] > fine["temporal_G_defect"]
    assert all(len(row["validation_samples"]) == 32 for row in record["solutions"])


@pytest.mark.parametrize("path,value", [
    (("settings", "eta_nodes"), 257),
    (("execution", "scalar_kernel_values"), 30_000_001),
    (("execution", "cpu_seconds_max"), 121.),
    (("execution", "field_solves"), True),
    (("slice_provider", "cache_hits"), 50),
    (("direct_samples", 25, "g1_V"), 99.),
    (("direct_samples", 25, "field_points"), 1),
    (("direct_samples", 25, "slice_costs", "kernel_values"), -1),
    (("heldout_used_for_training",), True),
    (("solutions", 0, "training_times", 1), 1/64),
    (("solutions", 1, "validation_samples", 10, "scaled_g1_V_defect"), 1.),
    (("diagnostics", "maximum_absolute_defects", 1, "combined_F_defect"), 0.),
])
def test_settings_budgets_source_and_interpolation_tampering_fails(record, path, value):
    changed = deepcopy(record)
    target = changed
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    with pytest.raises(ValueError):
        tool.validate_record(changed)


def test_creation_only_authentication_has_no_physical_fallback(isolated, monkeypatch):
    from recursive_horizons import nsc_atlas_forcing_cache as cache
    result = tool.create(isolated)
    assert "summary" in result and "direct_samples" not in result
    before = {p.name: p.read_bytes() for p in isolated.iterdir()}
    monkeypatch.setattr(tool, "compare", lambda: pytest.fail("unexpected physical comparison"))
    monkeypatch.setattr(tool, "_prepare", lambda: pytest.fail("unexpected preparation"))
    monkeypatch.setattr(cache, "ScalarForcing", lambda *a, **k: pytest.fail("unexpected direct slice"))
    assert tool.check(isolated)["physical_solves"] == 0
    with pytest.raises(FileExistsError):
        tool.create(isolated)
    assert before == {p.name: p.read_bytes() for p in isolated.iterdir()}
    (tool.ROOT/"src/recursive_horizons/nsc_atlas_forcing_cache.py").write_bytes(b"changed source")
    with pytest.raises(ValueError, match="source binding"):
        tool.check(isolated)


@pytest.mark.parametrize("tamper", ["bytes", "rehash", "missing"])
def test_record_tampering_is_detected_even_with_updated_output_hash(isolated, tamper):
    tool.create(isolated)
    path = isolated/"cache.json"
    if tamper == "bytes":
        path.write_bytes(path.read_bytes()+b" ")
    else:
        record = json.loads(path.read_bytes())
        if tamper == "rehash":
            record["solutions"][1]["validation_samples"][16]["interpolated"]["F"][2] += 1.
        else:
            record["direct_samples"].pop(21)
        raw = tool.encoded(record)
        path.write_bytes(raw)
        manifest_path = isolated/"manifest.json"
        manifest = json.loads(manifest_path.read_bytes())
        manifest["outputs"]["cache.json"] = {"sha256": tool.digest(raw), "bytes": len(raw)}
        manifest_path.write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError):
        tool.check(isolated)


def test_source_drift_during_comparison_leaves_no_record(isolated, monkeypatch, record):
    def changed():
        (tool.ROOT/"docs/nsc-atlas-scalar-evolution.md").write_bytes(b"changed during calculation")
        return deepcopy(record)
    monkeypatch.setattr(tool, "compare", changed)
    with pytest.raises(RuntimeError, match="inputs or runtime changed"):
        tool.create(isolated)
    assert not isolated.exists()
