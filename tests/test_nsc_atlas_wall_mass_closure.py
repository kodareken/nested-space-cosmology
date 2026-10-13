"""Independent wall mass identities and fail-closed immutable provenance."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import shutil

import pytest
import sympy as sp

from recursive_horizons.nsc_frozen_sources import authenticate_source
from record_nsc_atlas_source_mass import PREREQUISITES as MASS_PREREQUISITES

REPOSITORY = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("wall_mass_record", REPOSITORY / "scripts/record_nsc_atlas_wall_mass_closure.py")
tool = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tool)


def test_complete_exact_certificate():
    assert tool.exact_algebra_checks() == tool.expected_certificate()
    assert len(tool.IDENTITIES) == 21
    assert len(tool.CONTROLS) == 6


def test_independent_dimensionful_jordan_projection_and_frame_derivative():
    # This calculation uses physical radius, field and proper time rather
    # than parsing the recorder's dimensionless expressions.
    R, M, phi, tang, norm, speed, beta, acceleration = sp.symbols("R M phi tang norm speed beta acceleration")
    Htt, Htn = sp.symbols("Htt Htn")
    Hnn = Htt+2*(speed*tang-beta*norm)/R
    Ttn = tang*norm-(tang*norm+phi*Htn)/3
    box_phi_squared = 2*(-tang*tang+norm*norm)  # scalar wave equation
    hessian_phi_squared_nn = 2*norm*norm+2*phi*Hnn
    Tnn = (tang*tang+norm*norm)/2+(box_phi_squared-hessian_phi_squared_nn)/6
    jordan_transport = R*R*(beta*Ttn-speed*Tnn)/2
    endpoint = M*phi*phi/12+R*R*phi*(speed*tang-beta*norm)/6
    rates = {R: speed, phi: tang, tang: Htt+acceleration*norm,
             norm: Htn+acceleration*tang, speed: acceleration*beta-M/R**2,
             beta: acceleration*speed}
    endpoint_rate = sum(sp.diff(endpoint, q)*qdot for q, qdot in rates.items())
    minimal_transport = R*R*(beta*tang*norm-speed*(tang*tang+norm*norm)/2)/2
    difference = sp.expand(minimal_transport+endpoint_rate-jordan_transport)
    assert sp.factor(difference.subs(M, R*(1-beta**2+speed**2)/2)) == 0
    assert sp.diff(endpoint_rate, acceleration).expand() == 0
    # Omitting the mixed improvement is a genuine change in this comparison.
    wrong_transport = R*R*(beta*tang*norm-speed*Tnn)/2
    assert sp.factor((minimal_transport+endpoint_rate-wrong_transport).subs(M, R*(1-beta**2+speed**2)/2)) != 0


def test_both_fluxes_are_needed_even_with_zero_reduced_field_U_derivative():
    x, p, j, v, beta = sp.symbols("x p j v beta", nonzero=True)
    A, B = (j-v*p)/x, (j-beta*p)/x
    phi_U_times_Us = (A-B)/2
    radius_U_times_Us = (v-beta)/2
    assert sp.expand(x*phi_U_times_Us+p*radius_U_times_Us) == 0
    assert sp.factor(phi_U_times_Us) == -p*(v-beta)/(2*x)
    first_flux = -x*x*(v+beta)*(A-B)**2/8
    second_flux = -x*x*(v-beta)*(A+B)**2/8
    assert sp.factor(first_flux) != 0
    assert sp.factor(second_flux) != 0
    L = x*x*(beta*A*B-v*(A*A+B*B)/2)/2
    assert sp.expand(first_flux+second_flux-L) == 0


def test_full_mass_transform_retains_two_gradient_components():
    alpha, R, M, mu, p, A, B, v, beta = sp.symbols("alpha R M mu p A B v beta")
    omega = sp.sqrt(1-alpha*p*p/6)
    omega_t, omega_n = sp.diff(omega, p)*A, sp.diff(omega, p)*B
    inverse_metric = sp.diag(-1, 1)
    gradient_omega = sp.Matrix([omega_t, omega_n])
    gradient_R = sp.Matrix([v, beta])
    contraction = (gradient_omega.T*inverse_metric*gradient_R)[0]
    norm = (gradient_omega.T*inverse_metric*gradient_omega)[0]
    transformed = omega*(M+alpha*mu)-R*R*contraction-R**3*norm/(2*omega)
    coefficient = sp.diff(transformed, alpha).subs(alpha, 0)
    expected_endpoint = M*p*p/12+R*R*p*(v*A-beta*B)/6
    assert sp.expand(coefficient-(mu-expected_endpoint)) == 0
    assert sp.expand(coefficient.subs(B, 0)-coefficient) != 0


def test_beta_zero_is_regular_and_does_not_set_mass_response_to_zero():
    x, p, A, B, v, m = sp.symbols("x p A B v m", nonzero=True)
    L = -x*x*v*(A*A+B*B)/4
    Q = m*p*p/12+x*x*p*v*A/6
    assert L != 0 and Q != 0
    assert sp.denom(sp.factor(L)).has(x) is False
    assert sp.denom(sp.factor(Q)).has(x) is False


@pytest.mark.parametrize("name", list(tool.FORMULAS))
def test_declared_formula_corruption_is_detected(monkeypatch, name):
    formulas = dict(tool.FORMULAS)
    formulas[name] = "("+formulas[name]+")+p**2"
    monkeypatch.setattr(tool, "FORMULAS", formulas)
    with pytest.raises(ValueError, match="exact identity failed"):
        tool.exact_algebra_checks()


def test_wrong_runtime_and_cpu_limit_fail(monkeypatch):
    with monkeypatch.context() as patch:
        patch.setattr(sp, "__version__", "unexpected")
        with pytest.raises(RuntimeError, match="pinned SymPy"):
            tool.exact_algebra_checks()
    monkeypatch.setattr(tool.time, "process_time", lambda: 61)
    with pytest.raises(RuntimeError, match="process CPU budget"):
        tool.exact_algebra_checks(started=0)


@pytest.fixture
def isolated_root(tmp_path, monkeypatch):
    root = tmp_path / "repository"
    directories = {
        *(row["directory"] for row in tool.PREDECESSORS.values()),
        "results/development/nsc-atlas-characteristic-initialization-v1",
        "results/development/nsc-atlas-scalar-wall-consistency-v1",
        "results/development/nsc-cosmological-atlas-v7",
        "results/development/nsc-cosmological-atlas-v8",
    }
    sources = {tool.PRODUCER, tool.TESTS}
    original_contexts = {}
    for directory in directories:
        shutil.copytree(REPOSITORY / directory, root / directory)
        manifest = json.loads((REPOSITORY / directory / "manifest.json").read_bytes())
        sources.update(manifest["source_bindings"])
        for key in ("v7", "v8"):
            declaration = MASS_PREREQUISITES[key]
            if directory == declaration["directory"]:
                for relative, expected in manifest["source_bindings"].items():
                    if not relative.startswith("results/"):
                        original_contexts[relative] = (expected, declaration["source_commit"])
    for relative in sources:
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if relative in original_contexts:
            expected, commit = original_contexts[relative]
            # A temporary historical fixture uses its declared original text,
            # while every numerical input above is copied from this root.
            target.write_bytes(authenticate_source(
                REPOSITORY, relative, expected, context=commit))
        else:
            shutil.copyfile(REPOSITORY / relative, target)
    owner = root / tool.OWNER
    owner.parent.mkdir(parents=True, exist_ok=True)
    owner.write_text("Fixture owner; not a scientific derivation.\n", encoding="utf-8")
    monkeypatch.setattr(tool, "ROOT", root)
    return root


def test_creation_authentication_and_preserved_contexts(isolated_root, monkeypatch):
    before = {key: (isolated_root / row["directory"] / row["output"]).read_bytes()
              for key, row in tool.PREDECESSORS.items()}
    output = isolated_root / tool.RELATIVE_OUTPUT
    assert tool.create(output)["scientific_numerical_solves"] == 0
    monkeypatch.setattr(tool, "exact_algebra_checks", lambda *a: pytest.fail("unexpected algebra replay"))
    assert tool.check(output)["authenticated_immediate_predecessor_outputs"] == 2
    for key, row in tool.PREDECESSORS.items():
        assert (isolated_root / row["directory"] / row["output"]).read_bytes() == before[key]
    record = json.loads((output / "wall-mass.json").read_bytes())
    assert record["knowledge_boundary"]["exterior_wall_mass_functional_closed"] is True
    assert record["knowledge_boundary"]["numerical_wall_mass_curve_evaluated"] is False
    assert record["knowledge_boundary"]["full_coupled_solution_established"] is False
    assert record["execution"]["numerical_mass_integrals"] == 0
    assert "not an algebraic premise" in record["borrowed"]["normal_force_record_role"]
    with pytest.raises(FileExistsError, match="creation-only"):
        tool.create(output)


def test_explicit_replay_authenticates_before_algebra(isolated_root, monkeypatch):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    assert tool.check_algebra(output)["exact_zero_residuals"] == len(tool.IDENTITIES)
    path = output / "wall-mass.json"
    path.write_bytes(path.read_bytes()+b"tamper")
    monkeypatch.setattr(tool, "exact_algebra_checks", lambda *a: pytest.fail("algebra before authentication"))
    with pytest.raises(ValueError, match="artifact binding mismatch"):
        tool.check_algebra(output)


@pytest.mark.parametrize("relative", [tool.OWNER, tool.PRODUCER, tool.TESTS])
def test_new_source_drift_fails(isolated_root, relative):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    path = isolated_root / relative
    path.write_bytes(path.read_bytes()+b"\n")
    with pytest.raises(ValueError, match="wall mass source binding mismatch"):
        tool.check(output)


@pytest.mark.parametrize("relative", [
    "docs/nsc-atlas-source-mass.md", "scripts/record_nsc_atlas_source_mass.py",
    "src/recursive_horizons/nsc_atlas_source_mass.py", "tests/test_nsc_atlas_source_mass.py",
    "docs/nsc-atlas-scalar-normal-force.md", "scripts/record_nsc_atlas_scalar_normal_force.py",
    "tests/test_nsc_atlas_scalar_normal_force.py",
])
def test_each_immediate_source_authenticates_before_either_checker(isolated_root, monkeypatch, relative):
    path = isolated_root / relative
    path.write_bytes(path.read_bytes()+b"\n")
    monkeypatch.setattr(tool.importlib.util, "spec_from_file_location", lambda *a: pytest.fail("loaded unauthenticated predecessor"))
    with pytest.raises(ValueError, match="current predecessor source differs"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()


@pytest.mark.parametrize("relative", [
    "results/development/nsc-atlas-source-mass-v1/mass.json",
    "results/development/nsc-atlas-scalar-normal-force-v1/normal-force.json",
    "results/development/nsc-atlas-characteristic-initialization-v1/initialization.json",
    "results/development/nsc-atlas-scalar-wall-consistency-v1/consistency.json",
    "results/development/nsc-cosmological-atlas-v7/preparation-access.json",
    "results/development/nsc-cosmological-atlas-v8/prepared-response.json",
])
def test_corrupt_evidence_chain_fails(isolated_root, relative):
    path = isolated_root / relative
    path.write_bytes(path.read_bytes()+b"tamper")
    with pytest.raises(ValueError, match="artifact binding mismatch"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)


@pytest.mark.parametrize("key", list(tool.PREDECESSORS))
def test_rehashed_predecessor_is_rejected(isolated_root, key):
    row = tool.PREDECESSORS[key]
    directory = isolated_root / row["directory"]
    path = directory / row["output"]
    path.write_bytes(path.read_bytes()+b"\n")
    manifest_path = directory / "manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    manifest["outputs"][path.name] = {"sha256": tool.sha(path.read_bytes()), "bytes": path.stat().st_size}
    manifest_path.write_bytes(tool.canonical_bytes(manifest))
    with pytest.raises(ValueError, match="pinned predecessor manifest mismatch"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)


@pytest.mark.parametrize("relative", [tool.OWNER, "scripts/record_nsc_atlas_source_mass.py",
    "results/development/nsc-cosmological-atlas-v8/prepared-response.json"])
def test_missing_context_fails(isolated_root, relative):
    (isolated_root / relative).unlink()
    with pytest.raises(ValueError, match="missing or unsafe"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)


def test_rehashed_new_claim_tamper_is_rejected(isolated_root):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    path = output / "wall-mass.json"
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
    extra = output / "unexpected.json"
    extra.write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError, match="output inventory mismatch"):
        tool.check(output)
    extra.unlink()
    path = output / "manifest.json"
    manifest = json.loads(path.read_bytes())
    manifest["predecessor_bindings"]["normal_force"]["role"] = "altered role"
    path.write_bytes(tool.canonical_bytes(manifest))
    with pytest.raises(ValueError, match="predecessor binding mismatch"):
        tool.check(output)


@pytest.mark.parametrize("relative", ["../outside", "/outside"])
def test_unsafe_paths_fail(isolated_root, relative):
    with pytest.raises(ValueError, match="missing or unsafe"):
        tool.read_relative(relative)


def test_incompatible_model_and_bad_certificate_leave_no_record(isolated_root, monkeypatch):
    records, bindings = tool.authenticate_predecessors()
    wrong = deepcopy(records)
    wrong["source_mass"]["inputs"]["xi"] = .25
    with monkeypatch.context() as patch:
        patch.setattr(tool, "authenticate_predecessors", lambda: (wrong, bindings))
        with pytest.raises(ValueError, match="model or dated knowledge boundary differs"):
            tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    monkeypatch.setattr(tool, "exact_algebra_checks", lambda *a: {})
    with pytest.raises(ValueError, match="exact wall mass certificate mismatch"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()


def test_source_change_during_creation_is_rejected(isolated_root, monkeypatch):
    binding = tool.source_bindings()
    calls = iter((binding, {**binding, tool.OWNER: "0"*64}))
    monkeypatch.setattr(tool, "source_bindings", lambda: next(calls))
    with pytest.raises(RuntimeError, match="source or prerequisite changed"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()
