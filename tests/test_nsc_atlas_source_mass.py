"""Initial-leg mass identities, independent kernels and strict evidence binding."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import shutil
import time

import numpy as np
import pytest
import sympy as sp

from recursive_horizons import nsc_atlas_source_mass as model

REPOSITORY = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("source_mass_record", REPOSITORY / "scripts/record_nsc_atlas_source_mass.py")
tool = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tool)


@pytest.fixture(scope="module")
def calculation():
    return model.mass_calculation()


def test_exact_mixed_frame_kernel_and_improved_stress_identities():
    assert model.exact_algebra_checks() == tool.expected_certificate()
    assert len(tool.IDENTITIES) == 10


def test_mass_definition_itself_gives_geometric_coefficient():
    alpha, x, h, hp, k, d = sp.symbols("alpha x h hp k d", real=True)
    radius = x+alpha*h
    radial_derivative = 1+alpha*hp
    transverse_over_C = d+alpha*k
    geometric_mass = radius*(1+4*radial_derivative*transverse_over_C)/2
    coefficient = sp.diff(geometric_mass, alpha).subs(alpha, 0)
    assert sp.expand(coefficient-(h*(1+4*d)/2+2*x*(k+d*hp))) == 0


def test_the_two_omitted_terms_and_wrong_sector_sign_break_checks(calculation):
    s = {k: np.array(v) for k, v in calculation["samples"].items()}
    omit_transverse = s["mu_einstein"]+.4*s["p"]**2/12+s["x"]**2*s["p"]*s["d"]*s["p_prime"]/3
    assert np.max(abs(omit_transverse-s["mu_geometric"])) > .01
    omit_mixed_stress = -s["x"]**2*s["d"]*s["T_xx_normalized"]
    assert np.max(abs(omit_mixed_stress-s["mu_jordan_prime_frame"])) > .01
    assert s["d"].min() > 0
    assert -s["mu_einstein"][-1] > 0  # reversing d reverses the mass flux sign
    x, a = sp.symbols("x a", positive=True)
    wrong_d = -(a/x-1)/4
    assert sp.diff(x*wrong_d, x)+sp.Rational(1, 4) != 0


def test_independent_original_kernel_and_mass_representations(calculation):
    s = {k: np.array(v) for k, v in calculation["samples"].items()}
    np.testing.assert_allclose(s["s"], s["s_kernel"], atol=2e-13)
    np.testing.assert_allclose(s["mu_jordan"], s["mu_geometric"], atol=1e-11)
    np.testing.assert_allclose(s["mu_jordan_prime_stress"], s["mu_jordan_prime_frame"], atol=2e-10)
    assert abs(s["s"][-1]) > .01
    assert calculation["checks"]["during_source_max_abs_frame_difference"] > 1
    assert calculation["execution"]["total_integrand_evaluations"] == 47520
    assert sum(calculation["execution"][k] for k in (
        "geometry_integrand_evaluations", "reused_constraint_integrand_evaluations",
        "mass_and_transverse_x_integrand_evaluations", "original_V_kernel_integrand_evaluations")) == 47520
    assert calculation["checks"]["refinement_is_rigorous_enclosure"] is False


def test_post_source_constant_negative_mass_and_corner(calculation):
    s = {k: np.array(v) for k, v in calculation["samples"].items()}
    for name in ("mu_jordan", "mu_einstein", "mu_geometric", "k", "s"):
        assert np.all(s[name][:33] == 0)
    np.testing.assert_array_equal(s["mu_jordan"][96:], s["mu_einstein"][96:])
    assert np.all(s["mu_jordan"][96:] < 0)
    np.testing.assert_allclose(s["mu_jordan"][96:], s["mu_jordan"][-1], atol=1e-12)
    independent_radius = (s["h"]+(.8-s["x"])*s["h_prime"])/2+s["integral_p_squared"]/12
    np.testing.assert_allclose(independent_radius[96:], s["mu_jordan"][96:], atol=1e-11)
    assert np.all(s["p"][96:] == 0)
    assert np.all(s["s"][96:] < 0)  # no vacuum future inferred from the tangent source ending


def test_zero_and_quadratic_amplitude_controls():
    geometry = model.source_geometry()
    empty, _ = model.mass_samples(geometry, samples=17, amplitude=0)
    for name in tool.SAMPLE_KEYS-{"x", "V", "d"}:
        assert np.all(empty[name] == 0)
    unit, _ = model.mass_samples(geometry, samples=17)
    scaled, _ = model.mass_samples(geometry, samples=17, amplitude=-2)
    linear = {"p", "p_prime", "p_second", "s", "s_kernel", "s_prime",
              "integral_p_over_x_squared", "integral_V_kernel_b"}
    for key in linear:
        np.testing.assert_allclose(scaled[key], -2*unit[key], atol=1e-11)
    for key in tool.SAMPLE_KEYS-linear-{"x", "V", "d"}:
        np.testing.assert_allclose(scaled[key], 4*unit[key], atol=1e-10)


def test_null_coordinate_boost_preserves_geometric_mass():
    geometry = model.source_geometry()
    factor = 3.25
    boosted = model.null_coordinate_boost(geometry, factor)
    ordinary, _ = model.mass_samples(geometry, samples=17)
    transformed, _ = model.mass_samples(boosted, samples=17)
    for key in ("x", "p", "s", "d", "k", "mu_jordan", "mu_einstein", "mu_geometric"):
        np.testing.assert_allclose(transformed[key], ordinary[key], rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(transformed["integral_V_kernel_b"], ordinary["integral_V_kernel_b"]/factor, atol=1e-12)
    assert -3.2/boosted["U0"] == pytest.approx((-3.2/geometry["U0"])/factor)


def test_physical_unit_scaling_and_changed_source_strength(calculation):
    mu = calculation["last"]["mu_jordan"]
    epsilon, kappa, ell, scale = .2, .7, 3., 5.
    ordinary = model.physical_mass_length_coefficient(mu, epsilon, kappa, ell)
    units = model.physical_mass_length_coefficient(mu, epsilon, scale**2*kappa, scale*ell)
    physical = model.physical_mass_length_coefficient(mu, epsilon, kappa, scale*ell)
    assert units == pytest.approx(scale*ordinary)
    assert physical == pytest.approx(ordinary/scale)


@pytest.mark.parametrize("change", [{"xi": .25}, {"M_over_ell": .5}, {"bottom_scalar": 1},
    {"radius_fixed": True}, {"corner_mass_correction": 1}, {"corner_transverse_correction": 1},
    {"corner_transverse_correction": float("nan")}])
def test_incompatible_model_or_added_corner_datum_rejected(change):
    with pytest.raises(ValueError):
        model.validate_mass_inputs(**change)


def test_domain_resolution_and_coordinate_failures():
    geometry = model.source_geometry()
    with pytest.raises(ValueError):
        model.mass_samples(dict(geometry, U0=geometry["U0"]*.9))
    with pytest.raises(ValueError):
        model.mass_samples(geometry, samples=257)
    with pytest.raises(ValueError):
        model.mass_samples(geometry, order=128)
    with pytest.raises(ValueError):
        model.null_coordinate_boost(geometry, -1)
    with pytest.raises(ValueError):
        model.physical_mass_length_coefficient(float("nan"), 1, 1, 1)


def test_inherited_budget_counts_all_nodes_and_process_cpu(monkeypatch):
    budget = model.Budget(maximum_evaluations=2000)
    with pytest.raises(RuntimeError, match="budget exhausted"):
        model.mass_calculation(budget)
    assert budget.evaluations > 2000
    budget = model.Budget(maximum_evaluations=100000, started=0)
    monkeypatch.setattr(time, "process_time", lambda: 61)
    with pytest.raises(RuntimeError, match="budget exhausted"):
        budget.check()


@pytest.fixture
def isolated_root(tmp_path, monkeypatch):
    root = tmp_path / "repository"
    for row in tool.PREREQUISITES.values():
        shutil.copytree(REPOSITORY / row["directory"], root / row["directory"])
    initialization_manifest = json.loads((REPOSITORY / tool.PREREQUISITES["initialization"]["directory"] / "manifest.json").read_bytes())
    for relative in (tool.PRODUCER, tool.MODULE, tool.TESTS, *tool.ADAPTERS, *initialization_manifest["source_bindings"]):
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPOSITORY / relative, target)
    owner = root / tool.OWNER
    owner.parent.mkdir(parents=True, exist_ok=True)
    owner.write_text("Fixture owner; not a scientific derivation.\n", encoding="utf-8")
    monkeypatch.setattr(tool, "ROOT", root)
    return root


def test_create_and_stdlib_authentication_does_not_replay(isolated_root, monkeypatch):
    output = isolated_root / tool.RELATIVE_OUTPUT
    assert tool.create(output)["total_integrand_evaluations"] == 47520
    monkeypatch.setattr(model, "mass_calculation", lambda *a: pytest.fail("unexpected quadrature replay"))
    assert tool.check(output)["authenticated_predecessor_outputs"] == 17
    record = json.loads((output / "mass.json").read_bytes())
    assert all(v is False for v in record["knowledge_boundary"].values())
    assert "nonzero transverse s can persist" in record["meaning"]["post_source"]
    assert record["status"] == "QUADRATIC_SOURCE_LEG_INVARIANT_MASS_RESPONSE"
    with pytest.raises(FileExistsError, match="creation-only"):
        tool.create(output)


def test_bounded_replay(isolated_root):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    result = tool.check_numerics(output)
    assert result["maximum_sample_difference"] == 0
    assert result["total_integrand_evaluations"] == 47520


@pytest.mark.parametrize("key", ["v7", "v8", "scalar_wall", "initialization"])
def test_corrupt_predecessor_payload_fails(isolated_root, key):
    directory = isolated_root / tool.PREREQUISITES[key]["directory"]
    manifest = json.loads((directory / "manifest.json").read_bytes())
    path = directory / next(iter(manifest["outputs"]))
    path.write_bytes(path.read_bytes()+b"corrupt")
    with pytest.raises(ValueError, match="prerequisite artifact binding mismatch"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()


def test_rehashed_predecessor_is_not_a_new_implicit_context(isolated_root):
    directory = isolated_root / tool.PREREQUISITES["initialization"]["directory"]
    artifact = directory / "initialization.json"
    artifact.write_bytes(artifact.read_bytes()+b"\n")
    path = directory / "manifest.json"
    manifest = json.loads(path.read_bytes())
    manifest["outputs"][artifact.name] = {"sha256": tool.sha(artifact.read_bytes()), "bytes": artifact.stat().st_size}
    path.write_bytes(tool.canonical_bytes(manifest))
    with pytest.raises(ValueError, match="pinned prerequisite manifest mismatch"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)


@pytest.mark.parametrize("relative", [tool.OWNER, tool.PRODUCER, tool.MODULE, tool.TESTS])
def test_new_source_drift_fails(isolated_root, relative):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    path = isolated_root / relative
    path.write_bytes(path.read_bytes()+b"\n")
    with pytest.raises(ValueError, match="source binding mismatch"):
        tool.check(output)


@pytest.mark.parametrize("relative", [tool.ADAPTERS[0],
    "src/recursive_horizons/nsc_atlas_characteristic_initialization.py",
    "docs/nsc-atlas-characteristic-initialization.md"])
def test_borrowed_source_drift_is_explicitly_rejected(isolated_root, relative):
    path = isolated_root / relative
    path.write_bytes(path.read_bytes()+b"\n")
    with pytest.raises(ValueError, match="current source differs from explicitly borrowed context"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)


def test_missing_source_object_fails_before_record_creation(isolated_root):
    (isolated_root / tool.MODULE).unlink()
    with pytest.raises(ValueError, match="missing or unsafe"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()


def test_output_tamper_and_rehashed_scope_change_fail(isolated_root):
    output = isolated_root / tool.RELATIVE_OUTPUT
    tool.create(output)
    path = output / "mass.json"
    record = json.loads(path.read_bytes())
    record["knowledge_boundary"]["whole_exterior_mass_determined"] = True
    path.write_bytes(tool.canonical_bytes(record))
    with pytest.raises(ValueError, match="artifact binding mismatch"):
        tool.check(output)
    mpath = output / "manifest.json"
    manifest = json.loads(mpath.read_bytes())
    manifest["outputs"][path.name] = {"sha256": tool.sha(path.read_bytes()), "bytes": path.stat().st_size}
    mpath.write_bytes(tool.canonical_bytes(manifest))
    with pytest.raises(ValueError, match="declaration mismatch"):
        tool.check(output)


def test_failed_budget_leaves_no_output(isolated_root, monkeypatch):
    monkeypatch.setattr(model, "MAX_EVALUATIONS", 100)
    with pytest.raises(RuntimeError, match="budget exhausted"):
        tool.create(isolated_root / tool.RELATIVE_OUTPUT)
    assert not (isolated_root / tool.RELATIVE_OUTPUT).exists()


def test_incompatible_initialization_is_not_combined(isolated_root, monkeypatch):
    records, bindings = tool.authenticate_prerequisites()
    changed = deepcopy(records)
    changed["initialization"]["initialization.json"]["inputs"]["xi"] = .25
    monkeypatch.setattr(tool, "authenticate_prerequisites", lambda: (changed, bindings))
    with pytest.raises(ValueError, match="domains differ"):
        tool.declaration()
