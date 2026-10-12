"""Analytic null-data checks, independent quadratures and provenance failures."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import shutil

import numpy as np
import pytest
import sympy as sp

from recursive_horizons import nsc_atlas_characteristic_initialization as model
from recursive_horizons.nsc_atlas_field_preparation import source_profile, radius_product

REPOSITORY = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "characteristic_record", REPOSITORY / "scripts/record_nsc_atlas_characteristic_initialization.py")
tool = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tool)


@pytest.fixture(scope="module")
def calculation():
    return model.initialization_calculation()


def test_exact_constraint_and_frame_checks():
    assert model.exact_algebra_checks() == dict.fromkeys(tool.IDENTITIES, "0")


def test_analytic_jet_matches_direct_symbolic_derivatives_and_original_profile():
    geometry = model.source_geometry()
    x = sp.symbols("x", positive=True)
    U, a, b = sp.symbols("U a b", real=True)
    V = (1-x/sp.Rational(4, 5))*sp.exp(x/sp.Rational(4, 5))/U
    z = 4*(V-(a+b)/2)/(b-a)
    p = (1-z*z)**4/x
    expressions = [p, sp.diff(p, x), sp.diff(p, x, 2)]
    independent = sp.lambdify((x, U, a, b), expressions, "numpy", cse=True)
    points = np.linspace(*geometry["support_x"], 11)[1:-1]
    observed = model.source_jet(points, geometry)
    expected = independent(points, geometry["U0"], geometry["Va"], geometry["Vb"])
    for actual, wanted in zip(observed, expected):
        np.testing.assert_allclose(actual, wanted, rtol=2e-12, atol=5e-9)
    Vpoints = radius_product(points)/geometry["U0"]
    np.testing.assert_allclose(observed[0]*points, source_profile(Vpoints, geometry["Va"], geometry["Vb"])[0],
                               rtol=1e-14, atol=1e-14)


def test_omitting_affine_change_breaks_independent_frame_identity():
    x = sp.symbols("x", positive=True)
    p = sp.Function("p")(x)
    null_source = -x*sp.diff(p, x)**2/3+x*p*sp.diff(p, x, 2)/6
    wrong = sp.simplify(null_source-sp.diff(x*p*p, x, 2)/12+x*sp.diff(p, x)**2/2)
    assert wrong == -p*sp.diff(p, x)/3
    assert wrong != 0
    assert sp.simplify(wrong+sp.diff(p*p, x)/6) == 0


def test_direct_parts_and_frame_agree_with_separate_refinements(calculation):
    samples = {k: np.array(v) for k, v in calculation["samples"].items()}
    for key in ("h", "h_prime"):
        np.testing.assert_allclose(samples[key], samples[key+"_direct"], rtol=1e-10, atol=2e-10)
    np.testing.assert_allclose(samples["h"]-samples["x"]*samples["p"]**2/12
                               +samples["integral_p_squared"]/6, samples["h_einstein"], atol=1e-13)
    categories = calculation["checks"]["quadrature_32_64_max_differences"]
    assert set(categories) == {"constraint_quadrature", "geometry_quadrature"}
    assert categories["constraint_quadrature"]["x"] == 0
    assert categories["geometry_quadrature"]["h"] < 1e-8
    assert calculation["execution"]["scalar_integrand_evaluations"] == 15360
    assert calculation["execution"]["geometric_integrand_evaluations"] == 1440
    assert calculation["execution"]["total_integrand_evaluations"] == 16800
    assert calculation["checks"]["refinement_is_rigorous_enclosure"] is False


def test_compact_edges_corner_and_nonzero_post_source_memory(calculation):
    s = {k: np.array(v) for k, v in calculation["samples"].items()}
    assert s["h"][0] == s["h_prime"][0] == 0
    for key in ("p", "p_prime", "p_second"):
        assert np.all(s[key][[0, 32, 96, 128]] == 0)
    assert np.all(s["h"][:33] == 0)
    assert np.all(s["h"][96:] < 0)
    assert np.all(s["h_prime"][96:] < 0)
    np.testing.assert_allclose(s["h_prime"][96:], s["h_prime"][-1], atol=1e-11)
    np.testing.assert_allclose(s["h"][96:], s["h"][96]+(s["x"][96:]-s["x"][96])*s["h_prime"][96],
                               atol=1e-11)
    assert np.all(s["h_second"][96:] == 0)


def test_zero_input_and_quadratic_amplitude_control():
    geometry = model.source_geometry()
    empty = model.constraint_samples(geometry, amplitude=0., samples=17)
    for key in tool.SAMPLE_KEYS-{"V", "x"}:
        assert np.all(empty[key] == 0)
    positive = model.constraint_samples(geometry, samples=17)
    reversed_scaled = model.constraint_samples(geometry, amplitude=-2., samples=17)
    for key in ("p", "p_prime", "p_second"):
        np.testing.assert_array_equal(reversed_scaled[key], -2*positive[key])
    for key in tool.SAMPLE_KEYS-{"V", "x", "p", "p_prime", "p_second"}:
        np.testing.assert_allclose(reversed_scaled[key], 4*positive[key], atol=1e-11)


def test_unit_rescaling_differs_from_changing_physical_length(calculation):
    epsilon, kappa, ell, scale = .2, .7, 3., 5.
    alpha = model.amplitude_parameter(epsilon, kappa, ell)
    same_alpha = model.amplitude_parameter(epsilon, scale**2*kappa, scale*ell)
    assert same_alpha == pytest.approx(alpha)
    h = calculation["last"]["h"]
    assert scale*ell*same_alpha*h == pytest.approx(scale*(ell*alpha*h))
    # Changing ell while holding epsilon and G fixed changes physical strength.
    changed_alpha = model.amplitude_parameter(epsilon, kappa, scale*ell)
    assert changed_alpha == pytest.approx(alpha/scale**2)
    assert scale*ell*changed_alpha*h == pytest.approx((ell*alpha*h)/scale)


@pytest.mark.parametrize("change", [
    {"xi": .25}, {"M_over_ell": .3}, {"H_ell": 2}, {"k_ell": 3},
    {"scalar_kinetic_coefficient": 2}, {"scalar_potential": 1},
    {"surface_scalar_coupling": 1}, {"corner_radius_correction": 1},
    {"corner_derivative_correction": 1}, {"bottom_scalar": 1},
    {"lapse_fixed": False}, {"radius_fixed": True}, {"xi": float("nan")},
])
def test_changed_theory_or_incompatible_characteristic_data_rejected(change):
    with pytest.raises(ValueError):
        model.validate_model(**change)


@pytest.mark.parametrize("samples,order", [(130, 64), (128, 64), (3, 64), (129, 128)])
def test_unadmitted_resolution_rejected(samples, order):
    with pytest.raises(ValueError):
        model.constraint_samples(model.source_geometry(), samples=samples, order=order)


def test_domain_exit_and_shifted_corner_are_rejected():
    geometry = model.source_geometry()
    with pytest.raises(ValueError, match="exits"):
        model.source_jet([.39, .5], geometry)
    changed = dict(geometry, U0=geometry["U0"]*.95)
    with pytest.raises(ValueError):
        model.constraint_samples(changed)
    changed = dict(geometry, support_V=[geometry["Va"], geometry["Vb"]])
    with pytest.raises(ValueError, match="changed source domain"):
        model.constraint_samples(changed)


def test_integrand_and_cpu_budgets_fail_without_reset(monkeypatch):
    budget = model.Budget(maximum_evaluations=100)
    with pytest.raises(RuntimeError, match="budget exhausted"):
        model.initialization_calculation(budget)
    assert budget.evaluations > 100
    budget = model.Budget(started=0)
    monkeypatch.setattr(model.time, "process_time", lambda: 61)
    with pytest.raises(RuntimeError, match="budget exhausted"):
        budget.check()


@pytest.fixture
def isolated_root(tmp_path, monkeypatch):
    root = tmp_path / "repository"
    for row in tool.PREREQUISITES.values():
        shutil.copytree(REPOSITORY / row["directory"], root / row["directory"])
    for relative in (tool.PRODUCER, tool.MODULE, tool.TESTS, *tool.ADAPTERS):
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPOSITORY / relative, target)
    owner = root / tool.OWNER
    owner.parent.mkdir(parents=True, exist_ok=True)
    owner.write_text("Fixture owner; not a scientific derivation.\n", encoding="utf-8")
    monkeypatch.setattr(tool, "ROOT", root)
    return root


def test_creation_authentication_is_not_an_implicit_replay(isolated_root, monkeypatch):
    output = isolated_root / tool.RELATIVE_OUTPUT
    assert tool.create(output)["scalar_integrand_evaluations"] == 15360
    monkeypatch.setattr(model, "initialization_calculation", lambda *a: pytest.fail("unexpected quadrature replay"))
    assert tool.check(output)["authenticated_predecessor_outputs"] == 16
    value = json.loads((output / "initialization.json").read_bytes())
    assert all(v is False for v in value["knowledge_boundary"].values())
    assert "T_UU need not vanish" in value["analytic_consequences"]["transverse_caveat"]
    assert set(json.loads((output / "manifest.json").read_bytes())["source_bindings"]) == {
        tool.OWNER, tool.PRODUCER, tool.MODULE, tool.TESTS}
    with pytest.raises(FileExistsError, match="creation-only"):
        tool.create(output)


def test_bounded_numerical_replay(isolated_root):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    result = tool.check_numerics(output)
    assert result["maximum_sample_difference"] == 0
    assert result["scalar_integrand_evaluations"] == 15360
    assert result["total_integrand_evaluations"] == 16800


@pytest.mark.parametrize("key", ["v7", "v8", "scalar_wall"])
def test_corrupt_predecessor_archives_fail(isolated_root, key):
    directory = isolated_root / tool.PREREQUISITES[key]["directory"]
    manifest = json.loads((directory / "manifest.json").read_bytes())
    name = next(iter(manifest["outputs"]))
    path = directory / name
    path.write_bytes(path.read_bytes()+b"corruption")
    with pytest.raises(ValueError, match="prerequisite artifact binding mismatch"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()


def test_rehashed_old_archive_is_not_silently_accepted(isolated_root):
    directory = isolated_root / tool.PREREQUISITES["scalar_wall"]["directory"]
    path = directory / "manifest.json"
    manifest = json.loads(path.read_bytes())
    artifact = directory / "consistency.json"
    artifact.write_bytes(artifact.read_bytes()+b"\n")
    manifest["outputs"][artifact.name] = {"sha256": tool.sha(artifact.read_bytes()), "bytes": artifact.stat().st_size}
    path.write_bytes(tool.canonical_bytes(manifest))
    with pytest.raises(ValueError, match="pinned prerequisite manifest mismatch"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)


@pytest.mark.parametrize("relative", [tool.OWNER, tool.PRODUCER, tool.MODULE, tool.TESTS])
def test_current_source_drift_fails(isolated_root, relative):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    path = isolated_root / relative
    path.write_bytes(path.read_bytes()+b"\n")
    with pytest.raises(ValueError, match="source binding mismatch"):
        tool.check(output)


def test_adapter_source_drift_requires_deliberate_new_context(isolated_root):
    path = isolated_root / tool.ADAPTERS[0]
    path.write_bytes(path.read_bytes()+b"\n")
    with pytest.raises(ValueError, match="current adapter differs"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)


def test_missing_source_and_missing_input_fail_clearly(isolated_root):
    path = isolated_root / tool.MODULE
    path.unlink()
    with pytest.raises(ValueError, match="missing or unsafe"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)


def test_changed_artifact_and_rehashed_overclaim_fail(isolated_root):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    path = output / "initialization.json"
    value = json.loads(path.read_bytes())
    value["knowledge_boundary"]["mass_shift_determined"] = True
    path.write_bytes(tool.canonical_bytes(value))
    with pytest.raises(ValueError, match="artifact binding mismatch"):
        tool.check(output)
    manifest_path = output / "manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    manifest["outputs"][path.name] = {"sha256": tool.sha(path.read_bytes()), "bytes": path.stat().st_size}
    manifest_path.write_bytes(tool.canonical_bytes(manifest))
    with pytest.raises(ValueError, match="declaration mismatch"):
        tool.check(output)


def test_exhausted_creation_leaves_no_record(isolated_root, monkeypatch):
    original = model.Budget
    monkeypatch.setattr(model, "Budget", lambda: original(maximum_evaluations=1))
    with pytest.raises(RuntimeError, match="budget exhausted"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()


def test_incompatible_evidence_is_not_combined(isolated_root, monkeypatch):
    records, bindings = tool.authenticate_prerequisites()
    changed = deepcopy(records)
    changed["scalar_wall"]["consistency.json"]["inputs"]["xi"] = .25
    monkeypatch.setattr(tool, "authenticate_prerequisites", lambda: (changed, bindings))
    with pytest.raises(ValueError, match="domains differ"):
        tool.declaration()
