"""Exact conditional shell identities and independent provenance failure checks.

These algebra tests check translations and conventions after the borrowed
geometric laws are supplied. They do not prove those laws or solve backreaction.
"""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import shutil

import pytest
import sympy as sp

REPOSITORY = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "scalar_wall_record", REPOSITORY / "scripts/record_nsc_atlas_scalar_wall_consistency.py")
tool = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tool)


def test_exact_certificate_and_meaningful_negative_controls():
    certificate = tool.exact_algebra_checks()
    assert certificate == tool.expected_certificate()
    assert all(value == "0" for value in certificate["zero_residuals"].values())
    assert all(certificate["nonzero_controls"].values())


def test_variational_action_and_product_rule_independently_recover_balance():
    """Differentiate the action and a wall function instead of declared jumps."""
    x = sp.symbols("x", real=True)
    phi = sp.Function("phi")(x)
    xi, k, kappa = sp.symbols("xi k kappa", nonzero=True)
    field = sp.symbols("field")
    f = 1/(2*kappa)-xi*field**2/2
    jump = 6*xi*k*phi
    amm_residual = (-4*k*f+2*jump*sp.diff(f, field)).subs(field, phi)+2*k/kappa
    # Improved-stress component calculation uses jumps of phi and its normal
    # derivative; the variational result is obtained from the action above.
    surface = xi*(2*phi*jump-2*k*phi**2)
    mixed = jump*sp.diff(phi, x)-xi*(
        2*jump*sp.diff(phi, x)+2*phi*sp.diff(jump, x)+2*k*phi*sp.diff(phi, x))
    assert sp.expand(surface+amm_residual) == 0
    assert sp.expand(sp.diff(surface, x)+mixed) == 0
    assert sp.expand(surface.subs(xi, sp.Rational(1, 6))) == 0
    assert sp.expand(mixed.subs(xi, sp.Rational(1, 6))) == 0
    assert sp.expand(surface.subs(xi, sp.Rational(1, 4))) != 0


def test_declared_formula_corruption_is_detected(monkeypatch):
    formulas = dict(tool.FORMULAS, flux_jump=tool.FORMULAS["flux_jump"]+" + k*phi*v")
    monkeypatch.setattr(tool, "FORMULAS", formulas)
    with pytest.raises(ValueError, match="exact identity failed"):
        tool.exact_algebra_checks()


def test_wrong_symbolic_runtime_fails_explicitly(monkeypatch):
    monkeypatch.setattr(sp, "__version__", "unexpected")
    with pytest.raises(RuntimeError, match="pinned SymPy 1.14.0"):
        tool.exact_algebra_checks()


@pytest.fixture
def isolated_root(tmp_path, monkeypatch):
    root = tmp_path / "repository"
    for row in tool.PREREQUISITES.values():
        shutil.copytree(REPOSITORY / row["directory"], root / row["directory"])
    for relative in (tool.PRODUCER, tool.TESTS):
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPOSITORY / relative, target)
    owner = root / tool.OWNER
    owner.parent.mkdir(parents=True, exist_ok=True)
    owner.write_text("Fixture owner for provenance, not a scientific proof.\n", encoding="utf-8")
    monkeypatch.setattr(tool, "ROOT", root)
    return root


def write_json(path, value):
    path.write_bytes(tool.canonical_bytes(value))


def test_creation_check_and_honest_scope(isolated_root, monkeypatch):
    output = isolated_root / tool.RELATIVE_OUTPUT
    assert tool.create(output) == {"created": output.name, "outputs": 1, "scientific_numerical_solves": 0}
    # The normal checker is stdlib-only and must not repeat even symbolic work.
    monkeypatch.setattr(tool, "exact_algebra_checks", lambda: pytest.fail("unexpected algebra replay"))
    assert tool.check(output)["authenticated_predecessor_outputs"] == 15
    manifest = json.loads((output / "manifest.json").read_bytes())
    assert set(manifest["source_bindings"]) == {tool.PRODUCER, tool.OWNER, tool.TESTS}
    value = json.loads((output / "consistency.json").read_bytes())
    assert value["execution"]["scientific_numerical_solves"] == 0
    assert value["knowledge_boundary"]["normal_force_condition_computed"] is False
    assert value["knowledge_boundary"]["bulk_einstein_response_solved"] is False
    assert value["knowledge_boundary"]["fixed_background_with_nonzero_probe_claimed_exact"] is False
    assert value["preparation_regularity"].startswith("C3")
    assert value["borrowed"]["source"] == "1910.07534v2"
    assert len((output / "consistency.json").read_bytes()) < 14000
    with pytest.raises(FileExistsError, match="creation-only"):
        tool.create(output)


def test_inverse_and_current_historical_owners_are_not_required(isolated_root):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    (isolated_root / "docs/nsc-cosmological-atlas.md").write_text("New narrative.\n", encoding="utf-8")
    assert not (isolated_root / "results/development/nsc-atlas-observer-inference-v1").exists()
    assert tool.check(output)["outputs"] == 1


@pytest.mark.parametrize("relative", [
    f"{row['directory']}/{name}"
    for row in tool.PREREQUISITES.values()
    for name in json.loads((REPOSITORY / row["directory"] / "manifest.json").read_bytes())["outputs"]
])
def test_every_predecessor_output_is_authenticated(isolated_root, relative):
    path = isolated_root / relative
    path.write_bytes(path.read_bytes()+b"changed")
    with pytest.raises(ValueError, match="prerequisite artifact binding mismatch"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()


@pytest.mark.parametrize("version", ["v7", "v8"])
def test_rewritten_predecessor_manifest_is_rejected(isolated_root, version):
    directory = isolated_root / tool.PREREQUISITES[version]["directory"]
    path = directory / "manifest.json"
    manifest = json.loads(path.read_bytes())
    payload = directory / "preparation-access.json"
    payload.write_bytes(payload.read_bytes()+b"\n")
    manifest["outputs"][payload.name] = {"sha256": tool.sha(payload.read_bytes()), "bytes": payload.stat().st_size}
    write_json(path, manifest)
    with pytest.raises(ValueError, match="pinned prerequisite manifest mismatch"):
        tool.build_record()


@pytest.mark.parametrize("relative", [
    "results/development/nsc-cosmological-atlas-v7/manifest.json",
    "results/development/nsc-cosmological-atlas-v8/prepared-response.json",
    tool.OWNER,
])
def test_missing_prerequisite_or_owner_fails_before_creation(isolated_root, relative):
    (isolated_root / relative).unlink()
    with pytest.raises(ValueError, match="missing or unsafe prerequisite/source"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()


@pytest.mark.parametrize("relative", [tool.PRODUCER, tool.OWNER, tool.TESTS])
def test_bound_source_drift_fails(isolated_root, relative):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    path = isolated_root / relative
    path.write_bytes(path.read_bytes()+b"\n")
    with pytest.raises(ValueError, match="scalar/wall source binding mismatch"):
        tool.check(output)


def test_tampering_even_with_rehashed_output_is_rejected(isolated_root):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    path = output / "consistency.json"
    value = json.loads(path.read_bytes())
    value["knowledge_boundary"]["bulk_einstein_response_solved"] = True
    write_json(path, value)
    with pytest.raises(ValueError, match="scalar/wall artifact binding mismatch"):
        tool.check(output)
    manifest_path = output / "manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    manifest["outputs"][path.name] = {"sha256": tool.sha(path.read_bytes()), "bytes": path.stat().st_size}
    write_json(manifest_path, manifest)
    with pytest.raises(ValueError, match="scalar/wall declaration mismatch"):
        tool.check(output)


def test_inventory_and_prerequisite_bindings_are_strict(isolated_root):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    stray = output / "extra.json"
    stray.write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError, match="scalar/wall output inventory mismatch"):
        tool.check(output)
    stray.unlink()
    path = output / "manifest.json"
    manifest = json.loads(path.read_bytes())
    manifest["prerequisite_bindings"]["v7"]["source_commit"] = "0"*40
    write_json(path, manifest)
    with pytest.raises(ValueError, match="scalar/wall prerequisite binding mismatch"):
        tool.check(output)


@pytest.mark.parametrize("change,match", [
    (lambda r: r["v8"]["prepared-response.json"]["inputs"].update(M_over_ell=0.3), "inputs or preparation differ"),
    (lambda r: r["v8"]["prepared-response.json"]["domain"].update(Va=-3), "retained domains differ"),
    (lambda r: r["v8"]["prepared-response.json"]["domain"].update(arrival_window_over_ell=[0, 1]), "observer windows differ"),
])
def test_incompatible_evidence_is_not_combined(isolated_root, monkeypatch, change, match):
    records, bindings = tool.authenticate_prerequisites()
    changed = deepcopy(records)
    change(changed)
    monkeypatch.setattr(tool, "authenticate_prerequisites", lambda: (changed, bindings))
    with pytest.raises(ValueError, match=match):
        tool.build_record()


def test_broader_inverse_regularity_is_not_silently_substituted(isolated_root, monkeypatch):
    records, bindings = tool.authenticate_prerequisites()
    for version in ("v7", "v8"):
        records[version]["preparation-access.json"]["preparation"]["regularity"] = "C2 reduced field only"
    records["v8"]["prepared-response.json"]["preparation"]["regularity"] = "C2 reduced field only"
    monkeypatch.setattr(tool, "authenticate_prerequisites", lambda: (records, bindings))
    with pytest.raises(ValueError, match="original C3 preparation"):
        tool.build_record()


def test_failed_algebra_leaves_no_evidence(isolated_root, monkeypatch):
    monkeypatch.setattr(tool, "exact_algebra_checks", lambda: {})
    with pytest.raises(ValueError, match="exact symbolic certificate mismatch"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()


def test_source_change_during_creation_fails(isolated_root, monkeypatch):
    bindings = tool.source_bindings()
    calls = iter((bindings, {**bindings, tool.OWNER: "0"*64}))
    monkeypatch.setattr(tool, "source_bindings", lambda: next(calls))
    with pytest.raises(RuntimeError, match="source changed during"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()
