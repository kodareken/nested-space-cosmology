"""Small confirmation controls only; no trajectory/campaign is launched here."""
from copy import deepcopy
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons import nsc_nested_parent_child as model


@pytest.fixture(scope="module")
def confirmation():
    path = Path(__file__).resolve().parents[1] / "scripts/derive_nsc_nested_parent_child_confirmation.py"
    spec = importlib.util.spec_from_file_location("nested_pair_confirmation", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def small_states():
    pair = model.build_pair(64)
    state, _ = model.initial_state(pair, solve_constraints=False)
    states = [state.copy() for _ in range(3)]
    for i, point in enumerate(states):
        # Independent child/parent geometry coefficients enter Q; changing
        # radius, momenta and source does not change the conformal generator.
        point.Q[pair.geometry_child_indices[0]] += .03 * i
        point.Q[pair.geometry_parent_indices[1]] += .02 * i
        point.r += .01 * i
        point.p_Q += .001 * i
        point.phi0 *= np.exp(.03j * i)
        point.phi1 *= np.exp(.03j * i)
    return pair, states, np.array([0., .1, .2])


def test_cached_operator_equals_direct_live_Q_without_callback_rebuilds(confirmation, small_states, monkeypatch):
    pair, states, times = small_states
    direct = model.hamiltonian
    calls = []

    def counted(*args):
        calls.append(1)
        return direct(*args)

    monkeypatch.setattr(model, "hamiltonian", counted)
    cached = confirmation.CachedLinearHamiltonian(pair, states, times)
    assert len(calls) == 3
    for mark in (.013, .075, .101, .183):
        expected = direct(pair, confirmation.interpolate_state(states, times, mark))
        np.testing.assert_allclose(cached(mark), expected, rtol=1e-12, atol=1e-12)
    assert len(calls) == 3  # The expensive core function is never in __call__.
    report = cached.validate_midpoints()
    assert report["relative_max"] < 1e-12
    assert report["same_linear_Q_schedule"] is True
    assert report["cached_matrices_saved_in_payload"] is False
    assert len(calls) == 5
    with pytest.raises(ValueError, match="trajectory"):
        cached(.3)
    # Mutation: an altered cache cannot silently change the physics schedule.
    cached.matrices[1, 0, 0] += .1
    with pytest.raises(ValueError, match="changes the frozen"):
        cached.validate_midpoints()


def _case(probability_effect=0., modal_effect=0., geometry_effect=0., initial_offset=0.):
    times = np.linspace(0., .3, 4)
    rows = []
    for mark in times:
        fraction = mark / .3
        rows.append({"time": mark,
            "windows": {"child": {"probability": 1. + initial_offset - .1 * fraction + probability_effect * fraction}},
            "child_occupation": np.array([.5 + initial_offset - .05 * fraction + modal_effect * fraction,
                                           .5 + initial_offset - .04 * fraction + .9 * modal_effect * fraction]),
            "metrics": {"child_r_proper_mean": 4. + initial_offset + .02 * fraction + geometry_effect * fraction},
            "field_energy": 3.})
    return {"rows": rows, "summary": {"target_reached": True, "minimum_r": 3., "minimum_Q": .2,
        "maximum_gram_gap": 1e-12, "minimum_covariance_eigenvalue": .25,
        "maximum_covariance_eigenvalue": .75, "maximum_field_energy_closure": 1e-13,
        "sampled_constraint_end_margin": .0001}}


def _comparisons(confirmation):
    witness = {"results": {confirmation.PRIMARY_CASE: _case(),
                "nf256_parent_only_dt0.0005": _case(.01, .001, .02, initial_offset=.7)},
        "two_way_effects": {name: {"time_indicator": 1e-9} for name in
                            ("parent_to_child_state", "parent_to_child_modal_state", "parent_to_child_geometry")}}
    results = {"nf512_baseline_dt0.0005": _case(initial_offset=.2),
               "nf512_parent_only_dt0.0005": _case(.01002, .00102, .02001, initial_offset=1.1)}
    return witness, results, confirmation.confirmation_comparisons(witness, results)


def test_space_indicator_keeps_initial_subtraction_observable_and_finest_denominator(confirmation):
    witness, results, comparison = _comparisons(confirmation)
    row = comparison["parent_to_child_state"]
    np.testing.assert_allclose(row["primary_curve"], [.0, .01 / 3, .02 / 3, .01], atol=1e-14)
    assert row["space_indicator"] == pytest.approx(.00002, abs=1e-14)
    assert row["space_fraction"] == pytest.approx(.00002 / .01002, abs=1e-12)
    assert row["fraction_of_primary_v1_effect"] == pytest.approx(.00002 / .01, abs=1e-12)
    assert row["initial_offset_removed"] is True
    assert row["time_confirmation_claimed"] is False
    assert row["matched_proper_time_claimed"] is False
    assert comparison["parent_to_child_modal_state"]["resolved_at_one_percent"] is False
    # A large offset has no role in the claimed dynamical movement.
    assert row["primary_curve"][0] == 0
    changed = deepcopy(results)
    changed["nf512_parent_only_dt0.0005"]["rows"][1]["time"] += .001
    with pytest.raises(ValueError, match="coordinate clocks"):
        confirmation.confirmation_comparisons(witness, changed)


def _local():
    error = {"movement": .001, "effect": 1., "fraction": .001, "within_one_percent": True}
    controls = {name: {"occupation_movement": .1, "coherence_movement": 0.,
        "own_reference_error": error.copy(), "reference_midpoint_indicator": error.copy(),
        "coherence_reference_error": {"effect": 0., "fraction": None}} for name in
        ("outside_drive_off", "outside_cross_removed", "parent_detail_drive_off", "child_detail_cross_removed")}
    return {"child": {"occupation_full": np.array([[.5, .5], [.4, .4]]),
        "reduction_error": error.copy(), "full_midpoint_indicator": error.copy(),
        "conditional_versus_autonomous": error.copy(), "coherence_error": {"effect": 0., "fraction": None},
        "controls": controls}, "parent": {"reduction_error": error.copy(), "full_midpoint_indicator": error.copy()},
        "allocation": {"child": {"time_indexed_exterior_propagator_bytes": 0},
                       "parent": {"time_indexed_exterior_propagator_bytes": 0}},
        "observer_preserved": True, "all_source_components_retained": True}


def test_own_omission_error_and_data_quality_are_required_not_only_arrival(confirmation):
    local = _local()
    assert confirmation.local_quality(local)["all_within_one_percent"]
    local["child"]["controls"]["parent_detail_drive_off"]["own_reference_error"]["fraction"] = .011
    assert not confirmation.local_quality(local)["all_within_one_percent"]
    case = _case()
    assert confirmation.case_quality(case)["CAR_admissible"]
    case["summary"]["maximum_covariance_eigenvalue"] = 1.01
    assert not confirmation.case_quality(case)["CAR_admissible"]
    case["summary"]["minimum_Q"] = -.001
    assert not confirmation.case_quality(case)["positive_chart"]


def test_combined_goal_uses_four_physical_effects_not_optional_modal_precision(confirmation):
    quality = confirmation.case_quality(_case())
    record = {"source_hashes_unchanged": True,
        "cached_operator_identity": {"same_linear_Q_schedule": True},
        "endpoint_nested_schur": {"nested_direct_gap": 1e-14, "joined_direct_gap": 1e-14},
        "primary_physical_effects": {name: {"resolved_at_one_percent": True} for name in confirmation.PRIMARY_PHYSICAL_EFFECTS},
        "results": {"nf512_baseline_dt0.0005": {}, "nf512_parent_only_dt0.0005": {}},
        "case_quality": {confirmation.PRIMARY_CASE: quality.copy(), "nf512_baseline_dt0.0005": quality.copy(),
                         "nf512_parent_only_dt0.0005": quality.copy()},
        "required_case_names": [confirmation.PRIMARY_CASE, "nf512_baseline_dt0.0005", "nf512_parent_only_dt0.0005"],
        "local_quality": confirmation.local_quality(_local()), "local_response_saved": True,
        "confirmation_states_saved": True, "optional_modal_precision_resolved": False}
    assert confirmation.goal_complete(record)
    for mutation in ("physical", "source", "local", "accounting", "saved"):
        changed = deepcopy(record)
        if mutation == "physical":
            changed["primary_physical_effects"]["parent_to_child_state"]["resolved_at_one_percent"] = False
        elif mutation == "source":
            changed["source_hashes_unchanged"] = False
        elif mutation == "local":
            changed["local_quality"]["all_within_one_percent"] = False
        elif mutation == "accounting":
            changed["case_quality"][confirmation.PRIMARY_CASE]["full_source_accounting"] = False
        else:
            changed["confirmation_states_saved"] = False
        assert not confirmation.goal_complete(changed)


def test_checkpoint_progress_final_overwrite_guard_and_payload_binding(confirmation, tmp_path):
    original = tmp_path / "resident-v1.json"
    original.write_bytes(b"sealed original\n")
    output, progress = tmp_path / "confirmation.json", tmp_path / "progress.json"
    writer = confirmation.Checkpoints(output, progress)
    arrays = {}
    curve = np.array([[.1, .2], [.3, .4]])
    bound = confirmation.pack_arrays({"curve": curve}, arrays, "local")
    saved = writer.save({"schema": confirmation.SCHEMA, "status": "RUNNING", "local": bound}, arrays)
    assert progress.exists() and progress.with_suffix(".npz").exists()
    assert not output.exists()
    assert saved["payload_bytes"] == progress.stat().st_size + progress.with_suffix(".npz").stat().st_size
    np.testing.assert_array_equal(confirmation.unpack_arrays(bound, arrays)["curve"], curve)
    arrays["local_curve"][0, 0] += .001
    with pytest.raises(ValueError, match="curve binding"):
        confirmation.unpack_arrays(bound, arrays)
    arrays["local_curve"][0, 0] -= .001
    writer.save({"schema": confirmation.SCHEMA, "status": "COMPLETE", "local": bound}, arrays, final=True)
    assert output.exists()
    assert original.read_bytes() == b"sealed original\n"
    with pytest.raises(FileExistsError, match="already exists"):
        confirmation.Checkpoints(output, progress)
    with pytest.raises(FileExistsError, match="final"):
        writer.save({}, arrays, final=True)


def test_checkpoint_byte_cap_stops_without_replacing_prior_progress(confirmation, tmp_path):
    output, progress = tmp_path / "confirmation.json", tmp_path / "progress.json"
    writer = confirmation.Checkpoints(output, progress, payload_limit=1024)
    writer.save({"status": "small"}, {"a": np.array([1.])})
    before_json, before_npz = progress.read_bytes(), progress.with_suffix(".npz").read_bytes()
    rng = np.random.default_rng(72)
    with pytest.raises(RuntimeError, match="64 MiB"):
        writer.save({"status": "oversized"}, {"a": rng.normal(size=10000)})
    assert progress.read_bytes() == before_json
    assert progress.with_suffix(".npz").read_bytes() == before_npz
