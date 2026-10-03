"""Manufactured column, tangent and centre-clock checks for the parent adapter.

No future held nonlinear run is executed. Records stay in a temporary directory.
"""
import json
from dataclasses import replace

import numpy as np
import pytest

from recursive_horizons import nsc_discovery_backend as backend
from recursive_horizons import nsc_discovery_coupled_memory as memory
from recursive_horizons import nsc_discovery_grandchild as grandchild
from recursive_horizons import nsc_discovery_leading_einstein as leading
from recursive_horizons import nsc_discovery_parent_response as parent
from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


@pytest.fixture(scope="module")
def binding():
    return parent.manufactured_fixture()


@pytest.fixture(scope="module")
def holder(binding):
    return binding["holder"]


@pytest.fixture(scope="module")
def state(binding):
    return binding["state"]


@pytest.fixture(scope="module")
def tangent(binding):
    return binding["tangent"]


def _gap(left, right):
    return float(np.max(np.abs(np.asarray(left) - np.asarray(right))))


def test_period_eight_windows_use_the_centre():
    regions = parent.regional_map(8.0)
    assert regions["centre"] == 4.0
    assert regions["child_interval"] == (3.5, 4.5)
    assert regions["protected_collar_interval"] == (3.0, 5.0)
    assert regions["parent_interval"] == (1.0, 7.0)
    assert regions["parent_annulus_intervals"] == ((1.0, 2.8), (5.2, 7.0))
    assert parent.signed_distance(8.0, 4.0) == 0.0
    assert parent.signed_distance(8.0, [3.5, 4.5, 0.0]).tolist() == [-0.5, 0.5, -4.0]
    child, collar, parent_window = (3.5, 4.5), (3.0, 5.0), (1.0, 7.0)
    assert child[0] > collar[0] and child[1] < collar[1]
    assert collar[0] > parent_window[0] and collar[1] < parent_window[1]
    assert regions["parent_annulus_intervals"][0][1] < collar[0]
    assert regions["parent_annulus_intervals"][1][0] > collar[1]


def test_holder_uses_supplied_weights_and_native_fft(holder, binding):
    assert holder.weights.shape == (3,)
    assert np.array_equal(holder.grid.fine.occupations, holder.weights)
    assert not np.array_equal(holder.weights, coupling.OCCUPATIONS)
    assert holder.geometry_metadata["W_is_identity"] is True
    assert holder.geometry_map.shape == (holder.grid.ng, holder.grid.ng)
    assert backend.operator_backend(holder.grid) == "fft"
    assert holder.child_interval == (3.5, 4.5)
    assert holder.clock_locations == (4.0,)
    assert holder.source_metadata["build_pair_called"] is False
    dense = galerkin.build_grid(16, quadrature=32, length=8.0, gauge="conformal")
    refused = replace(dense, fine=replace(dense.fine, occupations=np.array(holder.weights, copy=True)))
    with pytest.raises(ValueError, match="rank 6"):
        backend.make_fft_grid(refused)
    nodal_momentum = np.linspace(-0.1, 0.1, holder.grid.ng)
    nodal = leading.State(
        np.array(binding["state"].Q, copy=True), np.array(binding["state"].r, copy=True),
        nodal_momentum, np.zeros(holder.grid.ng),
        np.array(binding["state"].phi0, copy=True), np.array(binding["state"].phi1, copy=True),
    )
    encoded = leading.encode(holder, nodal)
    assert _gap(encoded.p_Q, holder.grid.dx_g * nodal_momentum) < 1e-12


def test_response_rates_match_centred_differences(holder, state, tangent):
    probe = 1.0e-6
    image = parent.response_rates(holder, state, tangent)
    plus = leading.rates(parent.with_weights(holder, holder.weights + probe * tangent.occupations),
                         leading.combine(state, parent._as_state(tangent), probe))
    minus = leading.rates(parent.with_weights(holder, holder.weights - probe * tangent.occupations),
                          leading.combine(state, parent._as_state(tangent), -probe))
    for name in leading.FIELDS:
        assert _gap(getattr(image, name), (getattr(plus, name) - getattr(minus, name)) / (2.0 * probe)) < 1e-7
    assert _gap(image.phi0, -1j * image.operator_phi0) < 1e-8
    assert _gap(image.phi1, -1j * image.operator_phi1) < 1e-8
    weight_only = parent.ParentTangent(
        np.zeros(holder.grid.ng), np.zeros(holder.grid.ng), np.zeros(holder.grid.ng), np.zeros(holder.grid.ng),
        np.zeros_like(state.phi0), np.zeros_like(state.phi1), np.array(tangent.occupations, copy=True),
    )
    weighted = parent.response_rates(holder, state, weight_only)
    frozen = parent.response_rates(holder, state, weight_only, control_mode="frozen_geometry")
    assert float(np.max(np.abs(frozen.p_Q))) == 0.0
    assert float(np.max(np.abs(weighted.force_L))) > 0.0
    fine = leading.fine_state(holder, state)
    system = leading.active_system(holder.grid, fine)
    base = coupling.source_from_columns(system, fine)

    def forces(scale):
        probed = parent.with_weights(holder, holder.weights + scale * tangent.occupations)
        moved = leading.combine(state, parent._as_state(tangent), scale)
        sample, active = parent._fine_bundle(probed, moved)
        return coupling.source_from_columns(active, sample)

    upper, lower = forces(probe), forces(-probe)
    for name, reported in (("force_L", image.force_L), ("force_Q", image.force_Q), ("force_beta", image.force_beta)):
        assert _gap(reported, (upper[name] - lower[name]) / (2.0 * probe)) < 1e-6
    assert base["multiplicity_applied_once"] is True


def test_advance_matches_leading_rk4_and_integrates_the_centre_clock(holder, state, tangent):
    step = 1.0e-4
    advanced, direction, clock = parent.advance_tangent(holder, state, tangent, step)
    reference = leading.rk4_step(holder, state, step)
    for name in leading.FIELDS:
        assert _gap(getattr(advanced, name), getattr(reference, name)) < 1e-10
    assert _gap(direction.occupations, tangent.occupations) == 0.0
    constant = [2.0, 2.0, 2.0, 2.0]
    assert parent._rk4_average(constant) == 2.0
    tau_dot = parent.centre_clock(holder, state)
    assert abs(clock["tau"] - tau_dot * step) / max(abs(tau_dot * step), 1e-12) < 1e-3
    matched = parent.matched_centre_readout(holder, advanced, direction, clock["delta_tau"])
    bracket = parent.nonlinear_endpoint_clock_bracket(holder, state, tangent, step, step=1.0e-5)
    assert bracket["uses_integrated_linear_delta_tau"] is False
    assert bracket["source"] == "nonlinear_endpoint"
    assert bracket["held_nonlinear_campaign"] is False
    tau_gap = abs(bracket["delta_tau"] - clock["delta_tau"])
    content_gap = abs(bracket["delta_child_regional_content_tau"] - matched["delta_child_regional_content_tau"])
    assert abs(clock["delta_tau"]) > 1.0e-10
    assert tau_gap < 1.0e-12
    assert abs(matched["delta_child_regional_content_tau"]) > 1.0e-4
    assert content_gap < 1.0e-8


def test_normalization_keeps_an_exterior_direction_and_preparation_is_labelled(holder, state):
    stacked = np.vstack((state.phi0, state.phi1))
    hermitian = np.diag(np.linspace(0.1, 0.3, stacked.shape[1]))
    gauge = stacked @ hermitian
    removed0, removed1, info = parent.explicit_normalization_derivative(
        state.phi0, state.phi1, gauge[:holder.grid.nf], gauge[holder.grid.nf:],
    )
    assert info["approximation"] is False
    assert _gap(np.vstack((removed0, removed1)), 0.0) < 1e-10
    rng = np.random.default_rng(19)
    raw = rng.normal(size=stacked.shape) + 1j * rng.normal(size=stacked.shape)
    exterior = raw - stacked @ (stacked.conj().T @ raw)
    exterior = exterior / np.linalg.norm(exterior)
    kept0, kept1, _kept = parent.explicit_normalization_derivative(
        state.phi0, state.phi1, exterior[:holder.grid.nf], exterior[holder.grid.nf:],
    )
    assert _gap(np.vstack((kept0, kept1)), exterior) < 1e-10
    tangent = parent.ParentTangent(
        1.0e-3 * np.ones(holder.grid.ng), np.zeros(holder.grid.ng),
        1.0e-4 * np.ones(holder.grid.ng), np.zeros(holder.grid.ng),
        1.0e-3 * (gauge + exterior)[:holder.grid.nf], 1.0e-3 * (gauge + exterior)[holder.grid.nf:],
        np.array([0.01, -0.01, 0.02]),
    )
    prepared = parent.preparation_tangent(holder, state, tangent, step=1.0e-6)
    assert prepared["available"] is True
    assert prepared["method"] == parent.PREPARATION_METHOD
    assert prepared["analytic_preparation_jacobian"] is False
    assert prepared["finite_difference_used_as_analytic"] is False
    assert prepared["dropped_state_derivative"] is False
    assert prepared["steps"] == [1.0e-6, 5.0e-7]
    assert prepared["two_step_radius_gap"] < 1.0e-7
    assert prepared["radius_only"] is False
    assert prepared["full_preparation_components"] == (
        "constraint", "momentum_constraint", "momenta", "global_norm",
    )
    assert prepared["two_step_momenta_gap"] < 1.0e-8
    assert prepared["two_step_global_norm_gap"] < 1.0e-6
    assert prepared["two_step_constraint_gap"] < 1.0e-4
    assert prepared["components"]["momenta"]["max_abs"] > 0.0
    output = prepared["tangent"]
    assert _gap(output.p_Q, tangent.p_Q) == 0.0
    assert _gap(output.occupations, tangent.occupations) == 0.0
    assert float(np.linalg.norm(np.vstack((output.phi0, output.phi1)))) > 1.0e-6


def test_variable_rows_match_the_inherited_split_and_stream_one_image(holder, state, binding):
    reduced = parent.split_state(holder, state, binding["basis"], binding["widths"])
    assert reduced["source_rank"] == 3
    assert reduced["retained_rows"] == 2
    assert reduced["correlation_frobenius"] > 1.0e-8
    assert reduced["rank6_guard_called"] is False
    report = parent.full_vs_reduced(holder, state, reduced, 1.0e-4)
    assert report["one_physical_metric"] is True
    assert report["dense_propagator_stored"] is False
    assert report["rk4_streaming_calls"] == 4
    assert report["opening_direct_gap"] < 1.0e-10
    assert report["opening_sequential_gap"] < 1.0e-10
    assert report["generator_gap"] < 1.0e-8
    assert report["sequential_generator_gap"] < 1.0e-8
    assert report["endpoint_field_gap"] < 1.0e-8
    assert report["endpoint_geometry_gap"] < 1.0e-8
    assert report["forces"]["force_L"]["cross_max_abs"] > 0.0
    assert report["omissions"]["conditional"]["autonomous"] is False
    assert report["omissions"]["autonomous"]["autonomous"] is True
    assert report["omissions"]["conditional"]["physical_law"] is False
    assert report["omissions"]["autonomous"]["claimed_regeneration"] is False
    assert report["omissions"]["conditional"]["residual_max_abs"] > 1.0e-10
    assert report["omissions"]["autonomous"]["residual_max_abs"] > 1.0e-8
    pure = state.copy()
    coefficients = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=complex)
    pure.phi0 = binding["basis"][:holder.grid.nf, :2] @ coefficients
    pure.phi1 = binding["basis"][holder.grid.nf:, :2] @ coefficients
    with pytest.raises(ValueError, match="correlation control"):
        parent.split_state(holder, pure, binding["basis"], binding["widths"])


def test_inherited_widths_call_the_grandchild_reducer():
    rng = np.random.default_rng(5)
    basis, _r = np.linalg.qr(rng.normal(size=(16, 8)) + 1j * rng.normal(size=(16, 8)))
    field = rng.normal(size=(16, 3)) + 1j * rng.normal(size=(16, 3))
    direct, sequential = parent.initial_hierarchy(field, basis, (2, 2, 4))
    owned_direct, owned_sequential = grandchild.initial_routes(field, basis)
    for key in owned_direct:
        assert _gap(direct[key], owned_direct[key]) == 0.0
    for key in owned_sequential:
        assert _gap(sequential[key], owned_sequential[key]) == 0.0

    def image(columns):
        phase = np.exp(1j * np.linspace(0.0, 1.0, columns.shape[0]))
        return phase[:, None] * np.roll(columns, 1, axis=0)

    first, second = parent.hierarchy_rates(image, basis, direct, sequential, (2, 2, 4))
    reference_first, reference_second = grandchild.route_rates(image, basis, direct, sequential)
    extended_first, extended_second = parent.extended_hierarchy_rates(image, basis, direct, sequential, (2, 2, 4))
    for left, right in ((first, reference_first), (second, reference_second),
                        (extended_first, reference_first), (extended_second, reference_second)):
        for key in right:
            assert _gap(left[key], right[key]) < 1e-12
    calls = {"count": 0}
    small, _small_r = np.linalg.qr(rng.normal(size=(16, 5)) + 1j * rng.normal(size=(16, 5)))
    small_field = rng.normal(size=(16, 3)) + 1j * rng.normal(size=(16, 3))
    small_direct, small_sequential = parent.initial_hierarchy(small_field, small, (2, 1, 2))

    def counted(columns):
        calls["count"] += 1
        return image(columns)

    parent.hierarchy_rates(counted, small, small_direct, small_sequential, (2, 1, 2))
    assert calls["count"] == 1
    with pytest.raises(memory.ImplementationGap):
        memory._require_rank6(np.zeros((4, 3), dtype=complex), np.ones(3))


def test_primary_content_is_the_spatial_window(holder, state):
    content = parent.child_regional_content(holder, state)
    localization = parent.proper_localization(holder, state)
    modal = parent.fixed_mode_observer(holder, state)
    assert localization["spatial_window"] is True
    assert localization["fixed_mode_observer"] is False
    assert modal["spatial_window"] is False
    assert modal["dense_covariance_formed"] is False
    assert modal["kernel"].shape == (holder.reference_columns.shape[1], holder.reference_columns.shape[1])
    assert abs(content - localization["child_probability"]) < 1e-12
    assert localization["child_proper_length"] > 0.0
    assert localization["protected_collar_probability"] + 1e-12 >= content
    assert localization["parent_annulus_probability"] >= 0.0
    spec = parent.specification()
    assert spec["fixed_mode_observer_is_spatial_window"] is False
    assert spec["build_pair_called"] is False
    assert spec["baseline_seed_is_a_measurement_veto"] is False
    assert spec["source_may_be_nonorthogonal"] is True
    assert spec["modal_observer_must_be_orthonormal"] is True
    assert spec["analytic_preparation_jacobian"] is False


def test_source_may_be_nonorthogonal_and_the_observer_may_not(holder):
    grid = galerkin.build_grid(holder.grid.nf, quadrature=holder.grid.nq, length=8.0, gauge="conformal")
    skewed = np.array(holder.source_columns, copy=True)
    skewed[:, 0] *= 1.4
    skewed[:, 1] += 0.25 * skewed[:, 0]
    accepted = parent.make_holder(
        grid, weights=np.array(holder.weights, copy=True), reference_columns=holder.reference_columns,
        source_columns=skewed,
    )
    assert accepted.source_metadata["source_orthonormal"] is False
    assert accepted.source_metadata["modal_observer_orthonormal"] is True
    bad_observer = np.array(holder.reference_columns, copy=True)
    bad_observer[:, 0] *= 2.0
    fresh = galerkin.build_grid(holder.grid.nf, quadrature=holder.grid.nq, length=8.0, gauge="conformal")
    with pytest.raises(ValueError, match="modal observer"):
        parent.make_holder(
            fresh, weights=np.array(holder.weights, copy=True), reference_columns=bad_observer,
            source_columns=holder.source_columns,
        )


def test_stage_authorization_order_and_input_immutability(binding, tmp_path):
    with pytest.raises(parent.BindingUnavailable):
        parent.prepare(tmp_path, {"available": False})
    with pytest.raises(parent.BindingUnavailable):
        parent.prepare(tmp_path, {"available": True, "episode_path": tmp_path / "missing-episode.json"})
    with pytest.raises(parent.BindingUnavailable):
        parent.predict(tmp_path)
    unwritten = parent.DEVELOPMENT_ROOT / "nsc-discovery-parent-response-unwritten"
    with pytest.raises(parent.AuthorizationOpen):
        parent.prepare(unwritten, binding)
    assert not unwritten.exists()
    sealed = parent.DEVELOPMENT_ROOT / "nsc-discovery-parent-traction-v1" / "response-overwrite"
    with pytest.raises(PermissionError, match="sealed"):
        parent.authorize_stage_output(
            sealed, frozen_producer=True, physical_binding=True, closure_matches=True,
        )
    authorized = parent.authorize_stage_output(
        parent.DEVELOPMENT_ROOT / "nsc-discovery-parent-response-authorized-probe",
        frozen_producer=True, physical_binding=True, closure_matches=True,
    )
    assert authorized.name == "nsc-discovery-parent-response-authorized-probe"
    assert not authorized.exists()
    source = tmp_path / "input-source.txt"
    source.write_text("seeded baseline amplitude 0.0013\n")
    digest = source.read_bytes()
    stage = tmp_path / "open-stage"
    prepared = parent.prepare(stage, dict(binding, input_paths=(source,)))
    assert prepared["scientific_pass"] is False
    assert prepared["production_authorized"] is False
    assert prepared["physical_binding"] is False
    assert prepared["held_direction_sealed"] is False
    predicted = parent.predict(stage)
    assert predicted["forecast_locked"] is True
    assert predicted["prediction_before_held_measurement"] is True
    assert predicted["held_arm"] is False
    assert predicted["admission"]["restriction_owner"].endswith("scaled_step_admission")
    opened = parent.measure(stage)
    assert opened["status"] == "open"
    assert opened["executed"] is False
    assert opened["scientific_pass"] is False
    checked = parent.check(stage)
    assert checked["contract_ok"] is True
    assert checked["scientific_pass"] is False
    assert checked["measure_status"] == "open"
    assert checked["inputs_unchanged"] is True
    assert source.read_bytes() == digest
    callbacks = parent.owner_callbacks()
    assert callbacks["called"] is False
    assert callbacks["baseline_seed_is_a_measurement_veto"] is False
    assert callbacks["baseline_seed_amplitude"] == 0.0013
    assert callbacks["prepare_parent"].__name__ == "prepare_parent"
    assert callbacks["reconstruct_parent_pair"].__name__ == "reconstruct_parent_pair"
    assert callbacks["load_parent_record"].__name__ == "load_parent_record"
    assert callbacks["load_parent_record"] is callbacks["load_parent_preparation"]
    assert callbacks["loader_return_shape"] == "(record, arrays)"
    sealed_stage = tmp_path / "sealed-stage"
    sealed_binding = dict(
        binding,
        input_paths=(source,),
        sealed_held_direction={
            "kind": "gradient",
            "relative_change": 0.05,
            "baseline_amplitude": 0.0013,
            "sealed_before_measurement": True,
        },
    )
    parent.prepare(sealed_stage, sealed_binding)
    parent.predict(sealed_stage)
    measured = parent.measure(sealed_stage)
    assert measured["status"] == "measured"
    assert measured["executed"] is True
    assert measured["scientific_pass"] is False
    assert measured["physical_held_campaign"] is False
    assert measured["fresh_prepared_source"] is True
    assert measured["complement_force_max_abs"] >= 0.0
    assert measured["streaming_calls"] == 4
    assert measured["held_direction"]["baseline_is_measurement_veto"] is False
    sealed_check = parent.check(sealed_stage)
    assert sealed_check["measure_status"] == "measured"
    assert sealed_check["scientific_pass"] is False
    assert source.read_bytes() == digest
    monkey_limit = parent.CHUNK_LIMIT_BYTES
    parent.CHUNK_LIMIT_BYTES = 32
    try:
        with pytest.raises(ValueError, match="64MiB"):
            parent._write_exclusive(tmp_path / "oversized", "too-big", {"schema": parent.SCHEMA, "scientific_pass": False},
                                    {"values": np.zeros(8)})
    finally:
        parent.CHUNK_LIMIT_BYTES = monkey_limit
    preview = json.loads(json.dumps(parent.specification()))
    assert "evolve_prepared_source" in preview["implemented"]
    assert "analytic_preparation_jacobian" in preview["not_implemented"]


def test_owner2_public_callbacks_preserve_record_arrays_and_canonical_state(tmp_path, monkeypatch):
    from recursive_horizons import nsc_discovery_parent_episode as episode_owner
    callbacks = parent.owner_callbacks()
    assert callbacks["load_parent_record"] is episode_owner.load_parent_record
    assert callbacks["reconstruct_parent_pair"] is episode_owner.reconstruct_parent_pair
    assert callbacks["state_from_arrays"] is episode_owner.state_from_arrays
    prepared = episode_owner.analytic_fixture(nf=16, rank=2)
    population = prepared["populations"][1]
    pair, state = episode_owner.pair_and_state_from_population(
        population, common_k=prepared["common_k"], sign=-1, nf=16, length=8.)
    arrays = episode_owner.arrays_from_parent(pair, state, clocks=np.zeros(3), clock_rates=np.zeros(3))
    # Parent producer records use source_columns/reference_columns; episode
    # checkpoint aliases also remain valid to the real reconstruction API.
    arrays["source_columns"] = pair.source_columns.copy()
    arrays["reference_columns"] = pair.reference_columns.copy()
    arrays["weights"] = pair.weights.copy()
    record = {"nf":16, "population":1, "sign":-1, "k_common":pair.common_k,
              "column_rank":2, "covariance_rank":2}
    seen=[]
    def authenticated_loader_oracle(path):
        seen.append(path)
        return record, arrays
    monkeypatch.setattr(episode_owner.prepared_parent, "_load", authenticated_loader_oracle)
    loaded, saved = callbacks["load_parent_record"](tmp_path)
    assert seen == [tmp_path] and loaded is record and saved is arrays
    restored = callbacks["reconstruct_parent_pair"](saved, loaded)
    decoded = callbacks["state_from_arrays"](saved)
    np.testing.assert_array_equal(restored.geometry_map, pair.geometry_map)
    np.testing.assert_array_equal(restored.weights, pair.weights)
    for name in parent.STATE_FIELDS:
        np.testing.assert_array_equal(getattr(decoded,name),getattr(state,name))
    closure = parent.source_closure()
    assert "lab/src/recursive_horizons/nsc_discovery_parent_episode.py" in closure
    assert "lab/src/recursive_horizons/nsc_discovery_parent.py" in closure
