"""Creation/tamper checks with synthetic data; no scientific solve runs."""
import json

import pytest

import record_nsc_atlas_wall_overlap as tool


def synthetic():
    domain = {"U0": -.32, "Ua": -.30, "Ub": -.27, "Va": -2.5, "Vb": -2.25}
    rows = [{"t": t, "U": domain["Ua"], "V": domain["Va"]+.25*t,
             **{key: (index+1)*1e-7*t for index, key in enumerate(tool.DIAGNOSTICS)}}
            for t in tool.TIMES]
    return {"schema": tool.SCHEMA, "status": tool.STATUS, "scope": tool.SCOPE,
        "runtime": {"python": "3.14.8", "numpy": "2.5.1", "scipy": "1.17.1"},
        "settings": dict(tool.SETTINGS), "domain": domain,
        "difference_convention": "radial constraints minus propagated rectangle metric at U=Ua",
        "refinement_is_rigorous_error_enclosure": False,
        "samples": rows, "maximum_absolute": tool.maxima(rows),
        "execution": {"cpu_seconds": 1., "cpu_seconds_max": 120.,
            "scalar_kernel_values": 100, "scalar_kernel_values_max": 30_000_000,
            "metric_values": 100, "source_quadrature_nodes": 100, "rhs_calls": 100,
            "field_solves": 1, "metric_solves": 1, "wall_solves": 0, "radial_evaluations": 9}}


@pytest.fixture
def isolated(tmp_path, monkeypatch):
    root = tmp_path/"repository"
    for name in tool.SOURCES:
        path = root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(("fixture "+name).encode())
    monkeypatch.setattr(tool, "ROOT", root)
    monkeypatch.setattr(tool, "runtime_versions", lambda: synthetic()["runtime"])
    monkeypatch.setattr(tool, "calculate", synthetic)
    return root/"result"


def test_fixed_contract_rejects_changed_scope_domain_settings_and_counts():
    assert tool.validate_record(synthetic())
    for path, value in ((["settings", "field_nodes"], 129),
                        (["settings", "metric_nodes"], 65.),
                        (["scope"], "full spacetime solved"),
                        (["domain", "Ua"], -.29),
                        (["samples", 1, "t"], .2),
                        (["samples", 3, "mu_J_difference"], float("nan")),
                        (["maximum_absolute", "zeta_E_difference"], 0.),
                        (["refinement_is_rigorous_error_enclosure"], True),
                        (["execution", "cpu_seconds"], 121),
                        (["execution", "scalar_kernel_values"], 30_000_001),
                        (["execution", "metric_solves"], True),
                        (["execution", "wall_solves"], 1),
                        (["execution", "radial_evaluations"], 8)):
        record = synthetic()
        target = record
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        with pytest.raises(ValueError):
            tool.validate_record(record)


def test_creation_is_preserved_and_check_performs_no_physical_calculation(isolated, monkeypatch):
    tool.create(isolated)
    before = {path.name: path.read_bytes() for path in isolated.iterdir()}
    monkeypatch.setattr(tool, "calculate", lambda: pytest.fail("physical calculation during check"))
    import recursive_horizons.nsc_atlas_bulk_response as bulk
    monkeypatch.setattr(bulk, "solve_bulk_rectangle", lambda *a, **k: pytest.fail("bulk solve during check"))
    checked = tool.check(isolated)
    assert checked["field_solves"] == checked["metric_solves"] == checked["wall_solves"] == 0
    with pytest.raises(FileExistsError, match="creation-only"):
        tool.create(isolated)
    assert before == {path.name: path.read_bytes() for path in isolated.iterdir()}


@pytest.mark.parametrize("change", ["source", "bytes", "rehashed-settings", "rehashed-sample", "inventory"])
def test_tamper_rejected_including_rehashed_false_declarations(isolated, change):
    tool.create(isolated)
    payload = isolated/"response.json"
    if change == "source":
        (tool.ROOT/"src/recursive_horizons/nsc_atlas_bulk_response.py").write_bytes(b"changed source")
    elif change == "bytes":
        payload.write_bytes(payload.read_bytes()+b"\n")
    elif change == "inventory":
        (isolated/"extra.json").write_text("{}")
    else:
        record = json.loads(payload.read_bytes())
        if change == "rehashed-settings":
            record["settings"]["radial_order"] = 2
        else:
            record["samples"][4]["zeta_E_difference"] += 1.
        raw = tool.encoded(record)
        payload.write_bytes(raw)
        manifest_path = isolated/"manifest.json"
        manifest = json.loads(manifest_path.read_bytes())
        manifest["outputs"]["response.json"] = {"sha256": tool.digest(raw), "bytes": len(raw)}
        manifest_path.write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError):
        tool.check(isolated)


def test_source_change_during_calculation_cannot_create_record(isolated, monkeypatch):
    def changed():
        (tool.ROOT/"src/recursive_horizons/nsc_atlas_wall_evolution.py").write_bytes(b"changed")
        return synthetic()
    monkeypatch.setattr(tool, "calculate", changed)
    with pytest.raises(RuntimeError, match="inputs or runtime changed"):
        tool.create(isolated)
    assert not isolated.exists()
