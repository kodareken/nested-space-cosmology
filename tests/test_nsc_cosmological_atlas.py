"""Admission and integration checks; these are not cosmological discoveries."""
from copy import deepcopy
import ast
import importlib.util
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.integrate import solve_ivp

from recursive_horizons import nsc_cosmological_atlas as model


def producer():
    path = Path(__file__).resolve().parents[1] / "scripts/derive_nsc_cosmological_atlas.py"
    spec = importlib.util.spec_from_file_location("atlas_producer", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_domains_and_unequal_children_are_explicit():
    atlas = model.build_atlas()
    scenario = atlas["scenarios"][0]
    regions = {r["id"]: r for r in scenario["regions"]}
    assert regions["parent"]["geometry"]["globally_spherical"] is False
    assert regions["hole_a"]["dimensionless_parameters"] != regions["hole_b"]["dimensionless_parameters"]
    assert regions["child_a"]["dimensionless_parameters"]["state_is_shared"] is False
    assert "S3" in regions["child_a"]["geometry"]["spatial_topology"]
    assert regions["child_a"]["geometry"]["symmetry_orbits"] == "S2"
    assert {e["kind"] for e in scenario["links"]} == model.LINK_KINDS
    assert atlas["fitted_observations"] == []


def test_incompatible_junction_and_negative_neck_source():
    with pytest.raises(ValueError, match="incompatible"):
        model.angular_junction(2, 3, 1, 1, 0)
    with pytest.raises(ValueError, match="timelike"):
        model.angular_junction(2, 2, -1, 1, 0)
    with pytest.raises(ValueError, match="null"):
        model.angular_junction(2, 2, 0, 0, 0)
    # beta=0 inside a horizon may occur along a timelike trajectory.
    assert model.angular_junction(2, 2, -1, -1, 1)["beta_minus"] == 0
    result = model.angular_junction(4, 4, .75, .5, 0, normal_minus=-1)
    assert result["sigma_geometric"] < 0
    assert result["full_junction_solved"] is False
    assert model.null_collision_residual(2, 6, 3, 4) == 0
    assert model.null_collision_residual(2, 6, 3, 5) != 0


def test_junction_and_clock_change_units_without_changing_dimensionless_law():
    base = model.angular_junction(4, 4, .75, .5, .1)
    scaled = model.angular_junction(4000, 4000, .75, .5, .1)
    assert scaled["four_pi_R_sigma"] == base["four_pi_R_sigma"]
    assert scaled["sigma_geometric"] * 1000 == pytest.approx(base["sigma_geometric"])
    assert model.schwarzschild_f(4000, 1000) == model.schwarzschild_f(4, 1)
    assert model.static_clock(1000, 4000, 1000, (0, 2000)) == pytest.approx(1000 * model.static_clock(1, 4, 1, (0, 2)))


def test_clock_exits_are_rejected():
    with pytest.raises(ValueError, match="exited"):
        model.static_clock(3, 4, 1, (0, 2))
    with pytest.raises(ValueError, match="horizon"):
        model.static_clock(1, 2, 1, (0, 2))


def test_sibling_correlation_and_descent_do_not_create_signaling():
    scenario = model.build_atlas()["scenarios"][0]
    assert model.earliest_arrival(scenario, "child_a", "child_b", 0., clock="tau") is None
    # Even a correlation with delay-like metadata must not become a causal edge.
    edge = next(e for e in scenario["links"] if e["id"] == "siblings")
    edge.update(minimum_delay=.1, clock="tau", valid_interval=[0, 10], status="declared_control")
    assert model.earliest_arrival(scenario, "child_a", "child_b", 0., clock="tau") is None


def test_causal_delays_and_validity_use_a_single_declared_clock():
    scenario = model.build_atlas()["scenarios"][0]
    for number, (left, right, delay, interval) in enumerate([
        ("child_a", "parent", 2., [0., 10.]), ("parent", "child_b", 3., [4., 10.])]):
        scenario["links"].append(model._link(str(number), left, right, "causal", ["causal_review"],
            "synthetic conditional graph control; not a physical route", status="declared_control",
            delay=delay, interval=interval, clock="tau"))
    assert model.earliest_arrival(scenario, "child_a", "child_b", 0., clock="tau") == 7.
    assert model.earliest_arrival(scenario, "child_a", "child_b", 6., clock="tau") is None
    assert model.earliest_arrival(scenario, "child_a", "child_b", 0., clock="other") is None
    scenario["links"][-2]["valid_interval"] = [float("nan"), 10.]
    with pytest.raises(ValueError, match="finite"):
        model.earliest_arrival(scenario, "child_a", "child_b", 0., clock="tau")


def test_finite_cutoff_cycles_and_missing_citations():
    atlas = model.build_atlas()
    atlas["scenarios"][0]["max_descent_depth"] = 1
    with pytest.raises(ValueError, match="cutoff"):
        model.validate_atlas(atlas)
    atlas = model.build_atlas()
    scenario = atlas["scenarios"][0]
    scenario["links"].append(model._link("cycle", "grandchild", "hole_a", "descent", ["smolin"], "invalid cycle"))
    with pytest.raises(ValueError, match="cycle"):
        model.validate_atlas(atlas)
    atlas = model.build_atlas()
    atlas["scenarios"][0]["links"][0]["source_refs"] = []
    with pytest.raises(ValueError, match="evidence"):
        model.validate_atlas(atlas)


def test_fixed_theory_allows_state_changes_but_rejects_retuning():
    original = model.build_atlas()["scenarios"][1]
    changed = deepcopy(original)
    changed["regions"][0]["dimensionless_parameters"]["m_over_a"] = 1001
    assert model.require_fixed_theory(original, changed)
    changed["theory"]["couplings"]["a"] = 1.1
    with pytest.raises(ValueError, match="theory changed"):
        model.require_fixed_theory(original, changed)
    changed = deepcopy(original)
    changed["regions"][1]["dimensionless_parameters"]["a_over_reference_length"] = 1.1
    with pytest.raises(ValueError, match="fixed theory"):
        model.require_fixed_theory(original, changed)


def test_contained_interface_cannot_reset_descent_depth():
    atlas = model.build_atlas()
    scenario = atlas["scenarios"][0]
    interface = deepcopy(scenario["regions"][1])
    interface.update(id="hole_c", generation=1, contained_in="child_a")
    scenario["regions"].append(interface)
    for link in scenario["links"]:
        if link["id"] == "grandchild-descent":
            link["from"] = "hole_c"
    assert model.validate_atlas(atlas)
    # Even if a caller relabels the grandchild as generation 1, the lineage
    # hole_a -> child_a -> contained hole_c -> grandchild still has depth 2.
    scenario["regions"][-2]["generation"] = 1
    scenario["max_descent_depth"] = 1
    with pytest.raises(ValueError, match="cutoff"):
        model.validate_atlas(atlas)


def test_supplied_adapter_agrees_with_independent_matrix_exponential_and_units():
    h = np.array([[1., .2, .1], [.2, 1.5, 0.], [.1, 0., 2.]])
    state = np.diag([.2, .7, .4])
    grid = np.linspace(0., .2, 65)
    result = model.evolve_supplied_hamiltonian(h, grid, [1, 2], covariance=state, origin="test fixture")
    unitary = expm(-1j * grid[-1] * h)
    exact = (unitary @ state @ unitary.conj().T)[np.ix_([1, 2], [1, 2])]
    assert np.max(abs(result["covariance_total"][-1] - exact)) < 1e-7
    assert abs(result["covariance_total"][-1, 0, 0] - result["covariance_total"][-1, 1, 1]) > .2
    scaled = model.evolve_supplied_hamiltonian(h * 10, grid / 10, [1, 2], covariance=state, origin="unit rescaling")
    np.testing.assert_allclose(scaled["covariance_total"], result["covariance_total"], atol=1e-12, rtol=0)
    assert result["atlas_scope"]["physical_hamiltonian_derived"] is False
    with pytest.raises(ValueError, match="origin"):
        model.evolve_supplied_hamiltonian(h, grid, [1], covariance=state, origin="")


def test_hidden_exteriors_are_exactly_indistinguishable_for_this_probe():
    assert model.response_equivalence_control()["exact_difference"] == "0"
    state = np.diag([.2, .7, .8])
    results = []
    for hidden in (3, 5):
        h = np.array([[1., 1/3, 0.], [1/3, 2., 0.], [0., 0., hidden]])
        results.append(model.evolve_supplied_hamiltonian(h, np.linspace(0, .2, 17), [0],
                       covariance=state, origin="hidden-mode control")["covariance_total"])
    np.testing.assert_allclose(results[0], results[1], atol=1e-14, rtol=0)


def test_callable_noncommuting_hamiltonian_matches_independent_dop853():
    base = np.array([[1., .3], [.3, 2.]])
    pulse = np.array([[0., .2j], [-.2j, .1]])
    assert np.linalg.norm(base @ pulse - pulse @ base) > .1
    def hamiltonian(time):
        return base + np.sin(2 * time) * pulse
    state = np.diag([.2, .7])
    grid = np.linspace(0., .2, 129)
    actual = model.evolve_supplied_hamiltonian(hamiltonian, grid, [0], covariance=state,
                                             origin="noncommuting time-dependent control")
    solution = solve_ivp(lambda t, u: (-1j * hamiltonian(t) @ u.reshape(2, 2)).ravel(),
                         (0., .2), np.eye(2, dtype=complex).ravel(), method="DOP853",
                         rtol=1e-12, atol=1e-14)
    assert solution.success
    unitary = solution.y[:, -1].reshape(2, 2)
    exact = unitary @ state @ unitary.conj().T
    assert abs(actual["covariance_total"][-1, 0, 0] - exact[0, 0]) < 1e-7
    def changing_dimension(time):
        return base if time == 0 else np.eye(3)
    with pytest.raises(ValueError, match="changed dimension"):
        model.evolve_supplied_hamiltonian(changing_dimension, grid, [0], covariance=state, origin="invalid")


def test_creation_only_portable_binding_replay_and_vector_artifacts(tmp_path):
    driver = producer()
    output = tmp_path / "atlas"
    driver.create(output)
    assert driver.check(output)["outputs"] == 8
    for name in ("map.svg", "candidates.svg"):
        assert ET.fromstring((output / name).read_text())[0].tag.endswith("title")
    manifest = model.json.loads((output / "manifest.json").read_text())
    assert all(not Path(p).is_absolute() and "\\" not in p for p in manifest["source_bindings"])
    with pytest.raises(FileExistsError, match="creation-only"):
        driver.create(output)
    with (output / "map.svg").open("a") as handle:
        handle.write("changed")
    with pytest.raises(ValueError, match="binding mismatch"):
        driver.check(output)


def test_model_specification_binds_existing_response_and_callable_owners():
    """The ledger's numbers and source paths come from actual retained owners."""
    root = Path(__file__).resolve().parents[1]
    record = model.model_specification()
    response = json.loads((root / record["clocks"]["evidence_refs"][0]).read_text())
    for key in ("xi", "H_ell", "k_ell", "M_over_ell"):
        assert record["parameters"][key]["value"] == response["inputs"][key]
    assert record["domain"]["angular_mode"] == response["inputs"]["angular_mode"]
    assert record["domain"]["wall_radius_over_ell"] == response["domain"]["wall_radius_window_over_ell"]
    assert record["domain"]["exterior_radius_over_ell"] == [
        response["domain"]["radius_lower_bound_over_ell"],
        response["domain"]["radius_upper_bound_over_ell"],
    ]
    assert record["clocks"]["central_arrival_window_over_ell"] == response["domain"]["arrival_window_over_ell"]

    def check_references(value):
        if isinstance(value, dict):
            for key, child in value.items():
                if key in {"owner_refs", "evidence_refs"}:
                    assert all((root / ref.split("#", 1)[0]).is_file() for ref in child)
                check_references(child)
        elif isinstance(value, list):
            for child in value:
                check_references(child)
    check_references(record)
    callable_refs = [record["preparation"]["source_function_ref"]]
    callable_refs += [record["state"][key] for key in
                      ("reduced_equation_ref", "transmission_ref", "response_ref")]
    for reference in callable_refs:
        path, symbol = reference.split("#")
        tree = ast.parse((root / path).read_text(encoding="utf-8"))
        assert symbol in {node.name for node in tree.body if isinstance(node, ast.FunctionDef)}
    assert model.validate_model_specification(json.loads(model.canonical_bytes(record)))


def test_model_unit_choice_does_not_change_theory_data_or_legacy_outputs():
    atlas = model.canonical_bytes(model.build_atlas())
    programme = model.canonical_bytes(model.programme_record())
    original = model.model_specification()
    rescaled = model.model_specification(1000)
    assert rescaled["parameters"]["ell"]["value"] == 1000
    rescaled["parameters"]["ell"]["value"] = original["parameters"]["ell"]["value"]
    assert rescaled == original
    assert original["parameters"]["epsilon"]["value"] is None
    assert original["calibration"]["observation_refs"] == []
    assert set(original["parameters"]) & {row["id"] for row in original["derived_quantities"]} == set()
    # Standalone metadata does not rewrite the frozen candidate-family outputs.
    assert model.canonical_bytes(model.build_atlas()) == atlas
    assert model.canonical_bytes(model.programme_record()) == programme


@pytest.mark.parametrize("length", [0, -1, True, None, "1", float("nan"),
                                    float("inf"), complex(1, 0), 10**1000])
def test_model_rejects_invalid_supplied_length_scale(length):
    with pytest.raises(ValueError, match="length scale"):
        model.model_specification(length)


@pytest.mark.parametrize("key,value", [
    ("xi", 0), ("kinetic_coefficient", -1), ("scalar_mass_squared", 1),
    ("surface_scalar_coupling", .1), ("H_ell", 2), ("k_ell", 3),
    ("M_over_ell", .5), ("epsilon", .001), ("ell", 0), ("ell", "1"),
])
def test_model_fixed_inputs_cannot_be_retuned_under_the_same_identity(key, value):
    record = model.model_specification()
    record["parameters"][key]["value"] = value
    with pytest.raises(ValueError):
        model.validate_model_specification(record)


@pytest.mark.parametrize("section", ["parameters", "action", "state", "domain",
                                     "clocks", "preparation", "approximation"])
@pytest.mark.parametrize("invalid", [None, [], "missing"])
def test_model_rejects_malformed_sections_as_value_errors(section, invalid):
    record = model.model_specification()
    record[section] = invalid
    with pytest.raises(ValueError):
        model.validate_model_specification(record)


@pytest.mark.parametrize("section,key,value", [
    ("action", "f", "1/(2*kappa_4) + xi*phi^2/2"),
    ("action", "zero_background_scalar", 1),
    ("action", "scalar_einstein_tensor_accounting", "counted_twice"),
    ("state", "quantum_state", "vacuum"),
    ("state", "transmission_ref", "src/recursive_horizons/nsc_atlas_field_transfer.py#compact_bump"),
    ("state", "parameter_refs", {"nonminimal": "k_ell", "amplitude": "epsilon"}),
    ("domain", "parameter_refs", {"cap": "M_over_ell", "wall": "k_ell", "exterior": "H_ell", "length_scale": "ell"}),
    ("domain", "right_parent_driven_access", True),
    ("domain", "angular_mode", False),
    ("preparation", "source_function_ref", "src/recursive_horizons/nsc_atlas_field_preparation.py#radius_product"),
    ("preparation", "bottom_scalar", .1),
    ("preparation", "source_leg_gauge", "unanchored"),
    ("approximation", "field_order", True),
    ("approximation", "finite_amplitude_bound", .1),
    ("approximation", "full_bulk_and_wall_evolution_solved", True),
])
def test_model_rejects_incompatible_action_state_and_scope(section, key, value):
    record = model.model_specification()
    record[section][key] = value
    with pytest.raises(ValueError):
        model.validate_model_specification(record)


@pytest.mark.parametrize("section,key,window", [
    ("domain", "wall_radius_over_ell", [.45, .8]),
    ("domain", "exterior_radius_over_ell", [0, .55]),
    ("domain", "wall_radius_over_ell", [.55, .45]),
    ("domain", "wall_radius_over_ell", [.45]),
    ("clocks", "central_arrival_window_over_ell", [-.6, 0]),
    ("clocks", "central_arrival_window_over_ell", [float("nan"), -.2]),
    ("preparation", "normalized_support", [0, 1]),
])
def test_model_cannot_extend_or_relabel_the_recorded_finite_domain(section, key, window):
    record = model.model_specification()
    record[section][key] = window
    with pytest.raises(ValueError):
        model.validate_model_specification(record)


def test_model_separates_inputs_calibration_and_derived_evidence():
    original = model.model_specification()
    changed = deepcopy(original)
    changed["parameters"]["H_ell"]["role"] = "prediction"
    with pytest.raises(ValueError, match="role"):
        model.validate_model_specification(changed)
    changed = deepcopy(original)
    changed["parameters"]["H_ell"]["dimensionless"] = False
    with pytest.raises(ValueError, match="dimensional status"):
        model.validate_model_specification(changed)
    changed = deepcopy(original)
    changed["derived_quantities"][0]["id"] = "M_over_ell"
    with pytest.raises(ValueError, match="inputs cannot"):
        model.validate_model_specification(changed)
    changed = deepcopy(original)
    changed["derived_quantities"][-1]["status"] = "conditional_evaluated_response"
    with pytest.raises(ValueError, match="evidence scope"):
        model.validate_model_specification(changed)
    changed = deepcopy(original)
    changed["calibration"]["parameter_ids"] = ["H_ell"]
    with pytest.raises(ValueError, match="calibration"):
        model.validate_model_specification(changed)
    changed = deepcopy(original)
    changed["observational_predictions"] = ["cosmological_constant"]
    with pytest.raises(ValueError, match="observational predictions"):
        model.validate_model_specification(changed)
    changed = deepcopy(original)
    changed["unresolved"].remove("stability")
    with pytest.raises(ValueError, match="unresolved"):
        model.validate_model_specification(changed)


def test_model_requires_explicit_unknowns_and_portable_references():
    for section, key in (("clocks", "cross_universe_age_map"),
                         ("state", "quantum_state"), ("approximation", "finite_amplitude_bound")):
        record = model.model_specification()
        del record[section][key]
        with pytest.raises(ValueError, match="missing"):
            model.validate_model_specification(record)
    record = model.model_specification()
    del record["parameters"]["epsilon"]["value"]
    with pytest.raises(ValueError):
        model.validate_model_specification(record)
    for reference in ("../outside.json", "C:/private/evidence.json", "missing\\record.json"):
        record = model.model_specification()
        record["derived_quantities"][0]["evidence_refs"] = [reference]
        with pytest.raises(ValueError, match="repository-relative"):
            model.validate_model_specification(record)
