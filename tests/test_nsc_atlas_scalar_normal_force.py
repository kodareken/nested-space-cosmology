"""Normal projection/sign checks and authenticated immutable shell succession."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import shutil

import pytest
import sympy as sp

REPOSITORY = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("normal_force_record", REPOSITORY / "scripts/record_nsc_atlas_scalar_normal_force.py")
tool = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tool)
SURFACE_FORMULA = "2*xi*(6*xi-1)*k*phi**2"


def test_projection_gauss_orientation_and_reference_identities():
    assert tool.exact_algebra_checks(SURFACE_FORMULA) == tool.expected_certificate()
    assert len(tool.IDENTITIES) == 19
    assert len(tool.CONTROLS) == 5


def test_independent_timelike_shell_divergence_gives_plus_force_sign():
    sigma, K00, K11, K22, stress_jump = sp.symbols("sigma K00 K11 K22 stress_jump")
    inverse_gamma = sp.diag(-1, 1, 1)
    K = sp.diag(K00, K11, K22)
    upper_surface = -sigma*inverse_gamma
    # Gaussian normal Gamma^n_ab=-K_ab. The coefficient of the surface
    # normal divergence is therefore minus the surface/curvature contraction.
    normal_divergence = sum(-K[a, b]*upper_surface[a, b] for a in range(3) for b in range(3))
    contraction = sum(K[a, b]*upper_surface[a, b] for a in range(3) for b in range(3))
    assert sp.expand(normal_divergence+contraction) == 0
    solution = sp.solve(stress_jump+normal_divergence, stress_jump)
    assert len(solution) == 1 and sp.expand(solution[0]-contraction) == 0


def test_independent_one_sided_vacuum_and_scalar_traces_cancel():
    k, H, phi, mean_p = sp.symbols("k H phi mean_p", nonzero=True)
    mean_K = -3*H**2/(2*k)
    pp, pm = mean_p+k*phi/2, mean_p-k*phi/2
    Kp, Km = mean_K-3*k/2, mean_K+3*k/2
    # Shared tangential terms cancel. Retain the two different vacuum tensors.
    plus = pp**2/2+(2*Kp*phi*pp)/6
    minus = pm**2/2+(-3*H**2*phi**2+2*Km*phi*pm)/6
    assert sp.expand(plus-minus) == 0
    omitted_cap = pm**2/2+(2*Km*phi*pm)/6
    assert sp.expand(plus-omitted_cap) == -H**2*phi**2/2


def test_wrong_stress_formula_is_detected():
    with pytest.raises(ValueError, match="exact identity failed"):
        tool.exact_algebra_checks(SURFACE_FORMULA+"+k*phi**2")


def test_projection_formula_corruption_is_detected(monkeypatch):
    formulas = dict(tool.FORMULAS, normal_projection=tool.FORMULAS["normal_projection"]+"+G*phi**2")
    monkeypatch.setattr(tool, "FORMULAS", formulas)
    with pytest.raises(ValueError, match="exact identity failed"):
        tool.exact_algebra_checks(SURFACE_FORMULA)


def test_wrong_runtime_and_cpu_budget_fail(monkeypatch):
    with monkeypatch.context() as patch:
        patch.setattr(sp, "__version__", "unexpected")
        with pytest.raises(RuntimeError, match="pinned SymPy"):
            tool.exact_algebra_checks(SURFACE_FORMULA)
    monkeypatch.setattr(tool.time, "process_time", lambda: 61)
    with pytest.raises(RuntimeError, match="process CPU budget"):
        tool.exact_algebra_checks(SURFACE_FORMULA, started=0)


@pytest.fixture
def isolated_root(tmp_path, monkeypatch):
    root = tmp_path / "repository"
    predecessor_manifest = json.loads((REPOSITORY / tool.PREDECESSOR_DIRECTORY / "manifest.json").read_bytes())
    directories = [tool.PREDECESSOR_DIRECTORY]+[
        str(Path(row["manifest"]).parent) for row in predecessor_manifest["prerequisite_bindings"].values()]
    for directory in directories:
        shutil.copytree(REPOSITORY / directory, root / directory)
    for relative in (tool.PRODUCER, tool.TESTS, *predecessor_manifest["source_bindings"]):
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPOSITORY / relative, target)
    owner = root / tool.OWNER
    owner.parent.mkdir(parents=True, exist_ok=True)
    owner.write_text("Fixture owner; not a scientific derivation.\n", encoding="utf-8")
    monkeypatch.setattr(tool, "ROOT", root)
    return root


def test_creation_stdlib_check_and_unchanged_predecessor_history(isolated_root, monkeypatch):
    output = isolated_root / tool.RELATIVE_OUTPUT
    before = (isolated_root / tool.PREDECESSOR_DIRECTORY / "consistency.json").read_bytes()
    result = tool.create(output)
    assert result["scientific_numerical_solves"] == 0
    monkeypatch.setattr(tool, "exact_algebra_checks", lambda *a: pytest.fail("unexpected algebra replay"))
    assert tool.check(output)["authenticated_predecessor_outputs"] == 16
    assert (isolated_root / tool.PREDECESSOR_DIRECTORY / "consistency.json").read_bytes() == before
    assert json.loads(before)["knowledge_boundary"]["normal_force_condition_computed"] is False
    record = json.loads((output / "normal-force.json").read_bytes())
    assert record["knowledge_boundary"]["scalar_normal_stress_jump_computed"] is True
    assert record["knowledge_boundary"]["corrected_wall_trajectory_solved"] is False
    assert record["execution"]["new_floating_point_scientific_values"] is False
    with pytest.raises(FileExistsError, match="creation-only"):
        tool.create(output)


def test_explicit_algebra_replay_authenticates_first(isolated_root, monkeypatch):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    assert tool.check_algebra(output)["exact_zero_residuals"] == 19
    path = output / "normal-force.json"
    path.write_bytes(path.read_bytes()+b"changed")
    monkeypatch.setattr(tool, "exact_algebra_checks", lambda *a: pytest.fail("algebra ran before authentication"))
    with pytest.raises(ValueError, match="artifact binding mismatch"):
        tool.check_algebra(output)


@pytest.mark.parametrize("relative", [tool.OWNER, tool.PRODUCER, tool.TESTS])
def test_new_source_drift_fails(isolated_root, relative):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    path = isolated_root / relative
    path.write_bytes(path.read_bytes()+b"\n")
    with pytest.raises(ValueError, match="normal-force source binding mismatch"):
        tool.check(output)


@pytest.mark.parametrize("relative", ["docs/nsc-atlas-scalar-wall-consistency.md",
    "scripts/record_nsc_atlas_scalar_wall_consistency.py", "tests/test_nsc_atlas_scalar_wall_consistency.py"])
def test_original_source_drift_rejected_before_loading_code(isolated_root, monkeypatch, relative):
    path = isolated_root / relative
    path.write_bytes(path.read_bytes()+b"\n")
    monkeypatch.setattr(tool.importlib.util, "spec_from_file_location", lambda *a: pytest.fail("loaded unauthenticated code"))
    with pytest.raises(ValueError, match="current scalar-wall predecessor source differs"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()


@pytest.mark.parametrize("relative", [
    "results/development/nsc-atlas-scalar-wall-consistency-v1/consistency.json",
    "results/development/nsc-cosmological-atlas-v7/preparation-access.json",
    "results/development/nsc-cosmological-atlas-v8/prepared-response.json",
])
def test_corrupt_original_chain_fails(isolated_root, relative):
    path = isolated_root / relative
    path.write_bytes(path.read_bytes()+b"changed")
    with pytest.raises(ValueError, match="artifact binding mismatch"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)


def test_rehashed_predecessor_archive_is_rejected(isolated_root):
    directory = isolated_root / tool.PREDECESSOR_DIRECTORY
    artifact = directory / "consistency.json"
    artifact.write_bytes(artifact.read_bytes()+b"\n")
    manifest_path = directory / "manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    manifest["outputs"][artifact.name] = {"sha256": tool.sha(artifact.read_bytes()), "bytes": artifact.stat().st_size}
    manifest_path.write_bytes(tool.canonical_bytes(manifest))
    with pytest.raises(ValueError, match="pinned scalar-wall predecessor manifest mismatch"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)


@pytest.mark.parametrize("relative", [tool.OWNER, tool.PREDECESSOR_PRODUCER,
    "results/development/nsc-cosmological-atlas-v8/prepared-response.json"])
def test_missing_source_or_archive_fails(isolated_root, relative):
    (isolated_root / relative).unlink()
    with pytest.raises(ValueError, match="missing or unsafe"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)


def test_rehashed_new_claim_tamper_is_rejected(isolated_root):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    path = output / "normal-force.json"
    record = json.loads(path.read_bytes())
    record["knowledge_boundary"]["corrected_wall_trajectory_solved"] = True
    path.write_bytes(tool.canonical_bytes(record))
    manifest_path = output / "manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    manifest["outputs"][path.name] = {"sha256": tool.sha(path.read_bytes()), "bytes": path.stat().st_size}
    manifest_path.write_bytes(tool.canonical_bytes(manifest))
    with pytest.raises(ValueError, match="declaration mismatch"):
        tool.check(output)


def test_extra_output_and_rewritten_provenance_fail(isolated_root):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    extra = output / "extra.json"
    extra.write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError, match="output inventory mismatch"):
        tool.check(output)
    extra.unlink()
    path = output / "manifest.json"
    manifest = json.loads(path.read_bytes())
    manifest["predecessor_binding"]["authenticated_predecessor_outputs"] = 0
    path.write_bytes(tool.canonical_bytes(manifest))
    with pytest.raises(ValueError, match="predecessor binding mismatch"):
        tool.check(output)


@pytest.mark.parametrize("relative", ["../outside", "/outside"])
def test_unsafe_source_paths_rejected(isolated_root, relative):
    with pytest.raises(ValueError, match="missing or unsafe"):
        tool.read_relative(relative)


def test_incompatible_predecessor_model_and_algebra_failure_leave_no_record(isolated_root, monkeypatch):
    record, binding = tool.authenticate_predecessor()
    wrong = deepcopy(record)
    wrong["inputs"]["xi"] = .25
    with monkeypatch.context() as patch:
        patch.setattr(tool, "authenticate_predecessor", lambda: (wrong, binding))
        with pytest.raises(ValueError, match="model or dated knowledge boundary differs"):
            tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    monkeypatch.setattr(tool, "exact_algebra_checks", lambda *a: {})
    with pytest.raises(ValueError, match="exact normal-force certificate mismatch"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()


def test_source_change_during_creation_is_rejected(isolated_root, monkeypatch):
    binding = tool.source_bindings()
    calls = iter((binding, {**binding, tool.OWNER: "0"*64}))
    monkeypatch.setattr(tool, "source_bindings", lambda: next(calls))
    with pytest.raises(RuntimeError, match="source or prerequisite changed"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()
