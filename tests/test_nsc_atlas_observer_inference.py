"""Authenticate an analytic declaration; these tests do not prove its theorem."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import shutil

import pytest

REPOSITORY = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "observer_inference_record", REPOSITORY / "scripts/record_nsc_atlas_observer_inference.py")
tool = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tool)


@pytest.fixture
def isolated_root(tmp_path, monkeypatch):
    """Copy stored predecessors; a fixture proof tests binding, never proof truth."""
    root = tmp_path / "repository"
    for row in tool.PREREQUISITES.values():
        shutil.copytree(REPOSITORY / row["directory"], root / row["directory"])
    for relative in (tool.PRODUCER, tool.TESTS):
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPOSITORY / relative, target)
    owner = root / tool.OWNER
    owner.parent.mkdir(parents=True, exist_ok=True)
    owner.write_text("Fixture scientific owner: not a theorem proof.\n", encoding="utf-8")
    monkeypatch.setattr(tool, "ROOT", root)
    return root


def write_json(path, value):
    path.write_bytes(tool.canonical_bytes(value))


def test_creation_check_and_compact_scope(isolated_root):
    output = isolated_root / tool.RELATIVE_OUTPUT
    assert tool.create(output) == {"created": output.name, "outputs": 1, "scientific_numerical_solves": 0}
    checked = tool.check(output)
    assert checked["authenticated_predecessor_outputs"] == 15
    manifest = json.loads((output / "manifest.json").read_bytes())
    assert set(manifest["source_bindings"]) == {tool.PRODUCER, tool.OWNER, tool.TESTS}
    value = json.loads((output / "inference.json").read_bytes())
    assert value["execution"]["scientific_numerical_solves"] == 0
    assert value["execution"]["check_scope"].endswith("not machine verification of theorem truth")
    assert value["knowledge_boundary"]["observer_equivalence_established"] is False
    assert value["knowledge_boundary"]["second_physical_geometry_constructed"] is False
    assert value["connection_record"]["origin"] == "prescribed_at_boundary"
    assert value["connection_record"]["remote_preparation_established"] is False
    assert value["preparation_class"]["late_support_restriction"] is None
    assert value["preparation_class"]["sign_restriction"] is None
    assert value["preparation_class_test"]["scope"].startswith("exact within-model class exclusion")
    assert len((output / "inference.json").read_bytes()) < 20000
    with pytest.raises(FileExistsError, match="creation-only"):
        tool.create(output)


def test_current_historical_owners_are_not_silent_replay_inputs(isolated_root):
    # No old source tree is present in this fixture. Its authenticated evidence
    # must still be usable; no historical source imports or comparisons occur.
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    path = isolated_root / "docs/nsc-cosmological-atlas.md"
    path.write_text("Updated current atlas narrative.\n", encoding="utf-8")
    assert tool.check(output)["outputs"] == 1


@pytest.mark.parametrize("relative", [
    f"{row['directory']}/{name}"
    for row in tool.PREREQUISITES.values()
    for name in json.loads((REPOSITORY / row["directory"] / "manifest.json").read_bytes())["outputs"]
])
def test_every_stored_predecessor_output_is_authenticated(isolated_root, relative):
    path = isolated_root / relative
    path.write_bytes(path.read_bytes() + b"changed")
    with pytest.raises(ValueError, match="prerequisite artifact binding mismatch"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()


@pytest.mark.parametrize("version", ["v7", "v8"])
def test_pinned_predecessor_manifest_cannot_be_rewritten(isolated_root, version):
    directory = isolated_root / tool.PREREQUISITES[version]["directory"]
    manifest_path = directory / "manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    payload = directory / "preparation-access.json"
    payload.write_bytes(payload.read_bytes() + b"\n")
    manifest["outputs"][payload.name] = {"sha256": tool.sha(payload.read_bytes()), "bytes": payload.stat().st_size}
    write_json(manifest_path, manifest)
    with pytest.raises(ValueError, match="pinned prerequisite manifest mismatch"):
        tool.build_record()


@pytest.mark.parametrize("relative", [
    "results/development/nsc-cosmological-atlas-v7/manifest.json",
    "results/development/nsc-cosmological-atlas-v8/prepared-response.json",
    tool.OWNER,
])
def test_missing_prerequisite_or_proof_fails_before_creation(isolated_root, relative):
    (isolated_root / relative).unlink()
    with pytest.raises(ValueError, match="missing or unsafe prerequisite/source"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()


@pytest.mark.parametrize("relative", [tool.PRODUCER, tool.OWNER, tool.TESTS])
def test_current_bound_source_drift_fails(isolated_root, relative):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    path = isolated_root / relative
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="inference source binding mismatch"):
        tool.check(output)


def test_output_tampering_and_rehashed_declaration_fail(isolated_root):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    path = output / "inference.json"
    value = json.loads(path.read_bytes())
    value["knowledge_boundary"]["second_physical_geometry_constructed"] = True
    write_json(path, value)
    with pytest.raises(ValueError, match="inference artifact binding mismatch"):
        tool.check(output)
    manifest_path = output / "manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    manifest["outputs"][path.name] = {"sha256": tool.sha(path.read_bytes()), "bytes": path.stat().st_size}
    write_json(manifest_path, manifest)
    with pytest.raises(ValueError, match="inference declaration mismatch"):
        tool.check(output)


def test_inventory_and_prerequisite_binding_tampering_fail(isolated_root):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    stray = output / "extra.json"
    stray.write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError, match="inference output inventory mismatch"):
        tool.check(output)
    stray.unlink()
    path = output / "manifest.json"
    manifest = json.loads(path.read_bytes())
    manifest["prerequisite_bindings"]["v8"]["source_commit"] = "0" * 40
    write_json(path, manifest)
    with pytest.raises(ValueError, match="inference prerequisite binding mismatch"):
        tool.check(output)


@pytest.mark.parametrize("change,match", [
    (lambda r: r["v8"]["prepared-response.json"]["inputs"].update(M_over_ell=0.3), "inputs or preparation differ"),
    (lambda r: r["v8"]["prepared-response.json"]["domain"].update(Va=-3), "retained domains differ"),
    (lambda r: r["v8"]["prepared-response.json"]["domain"].update(arrival_window_over_ell=[0, 1]), "observer windows differ"),
])
def test_declared_compatibility_is_checked_after_authentication(isolated_root, monkeypatch, change, match):
    records, bindings = tool.authenticate_prerequisites()
    changed = deepcopy(records)
    change(changed)
    monkeypatch.setattr(tool, "authenticate_prerequisites", lambda: (changed, bindings))
    with pytest.raises(ValueError, match=match):
        tool.build_record()


@pytest.mark.parametrize("path,value", [
    (("operator", "analytic_norm_bound"), 1.0),
    (("checks", "analytic_contraction_independent_of_quadrature"), False),
    (("checks", "analytic_no_return_independent_of_quadrature"), False),
    (("observer", "analytic_no_return_log_margin_lower_bound"), 0.0),
])
def test_missing_certificate_conditions_are_not_promoted(isolated_root, monkeypatch, path, value):
    records, bindings = tool.authenticate_prerequisites()
    changed = deepcopy(records)
    for version in ("v7", "v8"):
        changed[version]["preparation-access.json"][path[0]][path[1]] = value
    monkeypatch.setattr(tool, "authenticate_prerequisites", lambda: (changed, bindings))
    with pytest.raises(ValueError, match="inverse prerequisites are not certified"):
        tool.build_record()


def test_source_change_during_creation_is_rejected(isolated_root, monkeypatch):
    bindings = tool.source_bindings()
    calls = iter((bindings, {**bindings, tool.OWNER: "0" * 64}))
    monkeypatch.setattr(tool, "source_bindings", lambda: next(calls))
    with pytest.raises(RuntimeError, match="source changed during"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()
