"""Production whole-cone campaign: synthetic slabs, resume, OPEN coverage."""
import copy
import json
import sys
from types import SimpleNamespace
from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("flint")

from recursive_horizons.nsc_ks_current_field_campaign import (
    PRODUCTION_ATOL, PRODUCTION_FAMILY, PRODUCTION_MAX_STEP, PRODUCTION_PROFILE_IDENTITY,
    PRODUCTION_RTOL, SCHEMA, STREAM_KIND_INJECTED, STREAM_KIND_TRAJECTORY,
    bind_authenticated_v4_payload, iterable_segment_stream, runtime_forecast,
    run_whole_cone_campaign, scientific_coverage, solver_matches_production,
    write_checkpoint_atomic,
)
from recursive_horizons.nsc_ks_current_field_cone import (
    V4_PROFILE_IDENTITY, WholeConeAccumulator, WholeConeCheckpoint,
    WholeConeFieldConfig, packed_radius_tail, profile_payload_digest,
    restore_v4_profile_payload, scientific_digest, scientific_record,
    with_residual_integrals, whole_cone_continuous_inputs,
)
from recursive_horizons.nsc_ks_difference_residual_polynomial import monomial_keys
from recursive_horizons.nsc_ks_local_profile_fourier_bound import (
    enclose_profile_fourier_local, restore_local_profile_fourier,
    serialize_local_profile_fourier,
)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_trajectory import TrajectorySegment
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


ROOT = Path(__file__).resolve().parents[1]
V4_PATH = ROOT / "results/development/nsc-ks-current-field-pilot-v4.json"


def small_family():
    coefficients = np.zeros((2, 8))
    coefficients[0, :4] = (0.012, -0.003, 0.0015, -0.0004)
    coefficients[1, :4] = (0.021, 0.004, -0.001, 0.0003)
    return LocalIncomingFamily(coefficients)


def synthetic_segment(rho_start, rho_end, *, scale=1.0, nsrc=1, nz=8):
    size = 2 * nsrc * (1 + nz)
    start = np.zeros(size, complex)
    start[0] = 0.5 * scale
    start[1] = 0.25 * scale
    harmonic = np.array([1, 1j, -1, -1j] * (nz // 4), complex) * 0.01 * scale
    start[2:2 + nz] = harmonic
    start[2 + nz:2 + 2 * nz] = 2 * harmonic
    return TrajectorySegment(
        rho_start, rho_end, start, 2 * start, np.zeros((6, size), complex))


def rho_ladder(start, count, step=1.0 / 8192):
    values = [float(start)]
    for _ in range(count):
        values.append(values[-1] - step)
    return values


@pytest.fixture(scope="module")
def small_cone():
    family = small_family()
    origin = family.center - 0.2
    bits = 80
    order = 1
    keys = monomial_keys(order)
    rows = enclose_profile_fourier_local(
        family, keys, origin, 0.4, retained_index=4, taylor_order=5,
        exterior_panels=16, transition_panels=16, interior_panels=16,
        max_derivative=2, flat_denominator=32, relative_tolerance_bits=20,
        absolute_tolerance_bits=30, max_depth=8, bits=bits)
    payload = serialize_local_profile_fourier(rows)
    digest = profile_payload_digest(payload)
    settings = {
        "integrator": "dop853",
        "rtol": PRODUCTION_RTOL.hex(),
        "atol": PRODUCTION_ATOL.hex(),
        "max_step": PRODUCTION_MAX_STEP.hex(),
        "step_control": "joint",
        "tangents": "zero",
        "retained_index": 4,
        "profile_payload_recomputed": False,
    }
    config = WholeConeFieldConfig(
        family=(14, 1),
        profile_identity=profile_identity(family, include_normal_window=True),
        bits=bits, time_degree=4, reciprocal_order=order, spatial_count=8,
        period_origin=origin, period_length_numerator=2,
        period_length_denominator=5, mass=0.4, angular=0.8,
        source_weights=(1.0,), source_energies=(0.7,),
        radius_tail=packed_radius_tail(family, 1.03, 0.8, order, bits=bits),
        profile_payload_sha256=digest, settings=settings)
    nodes = rho_ladder(1.03, 3)
    segments = (
        synthetic_segment(nodes[0], nodes[1]),
        synthetic_segment(nodes[1], nodes[2], scale=1.1),
        synthetic_segment(nodes[2], nodes[3], scale=0.9),
    )
    step = 1.0 / 8192
    complete = (
        synthetic_segment(1.0 + 2 * step, 1.0 + step),
        synthetic_segment(1.0 + step, 1.0, scale=1.1),
    )
    return {
        "family": family, "payload": payload, "digest": digest,
        "config": config, "segments": segments, "complete": complete,
    }


def run_campaign(cone, segments, tmp_path, **kwargs):
    checkpoint = kwargs.pop("checkpoint_path", tmp_path / "campaign.checkpoint.json")
    kwargs.setdefault("stream_kind", STREAM_KIND_INJECTED)
    kwargs.setdefault("source_via_select_source", False)
    kwargs.setdefault("rho_up", 1.03)
    kwargs.setdefault("free_space_floor", None)
    mapping, result = run_whole_cone_campaign(
        cone["config"], cone["payload"], cone["family"],
        iterable_segment_stream(segments), checkpoint_path=checkpoint, **kwargs)
    return mapping, result, checkpoint


def test_v4_payload_is_restored_without_recomputing_coefficients(monkeypatch):
    calls = []

    def forbidden(*args, **kwargs):
        calls.append(1)
        raise AssertionError("v4 payload must not be recomputed")

    monkeypatch.setattr(
        "recursive_horizons.nsc_ks_current_field_campaign.enclose_profile_fourier_local",
        forbidden, raising=False)
    payload, digest, record = bind_authenticated_v4_payload(ROOT)
    raw = V4_PATH.read_bytes()
    restored, expected = restore_v4_profile_payload(
        json.loads(raw), expected_file_sha256=sha256(raw).hexdigest(),
        file_bytes=raw)
    assert digest == expected
    assert profile_payload_digest(payload) == digest
    assert record["profile_identity"] == V4_PROFILE_IDENTITY
    again = restore_local_profile_fourier(payload, expected_sha256=digest)
    key = next(iter(restored))
    assert again[key].coefficients[0].real.contains(restored[key].coefficients[0].real)
    assert not calls


def test_continuous_inputs_keep_unowned_slots_as_named_nulls(small_cone):
    family = small_cone["family"]
    covariance = np.array([[0.5]], complex)
    channel = {"copy_count": 2, "degeneracy": 3, "angular_eigenvalue": 0.8}
    inputs = whole_cone_continuous_inputs(
        family, mass=0.4, angular=0.8, rho_up=1.03, source_energies=(0.7,),
        source_weights=(1.0,), source_covariance=covariance, channel=channel,
        bits=80)
    propagation = inputs["propagation"]
    matter = inputs["matter"]
    assert inputs["physical_rho1_source_error"] is None
    assert inputs["complete_for_propagate_difference_error"] is False
    assert inputs["complete_for_difference_matter_error"] is False
    assert propagation["reference_offdiagonal_integral"]["value"] is not None
    assert propagation["offdiagonal_integral"]["value"] is not None
    assert propagation["Bz_integral"]["value"] is not None
    assert propagation["maximum_absolute_energy"]["value"] is not None
    assert propagation["difference_initial"]["value"] is not None
    assert propagation["reference_initial"]["value"] is None
    assert "rho=1" in propagation["reference_initial"]["reason"]
    assert propagation["reference_residual"]["value"] is None
    assert propagation["M_integral"]["value"] is None
    assert propagation["Mz_integral"]["value"] is None
    assert matter["source_norm"]["value"] is not None
    assert matter["axial_lower"]["value"] is not None
    assert matter["radius_lower"]["value"] is not None
    assert matter["multiplicity"]["value"] is not None
    assert matter["difference_norm"]["value"] is None
    filled = with_residual_integrals(
        inputs, small_cone["config"].radius_tail)
    assert filled["propagation"]["difference_residual_integrals"]["value"] is not None
    assert filled["complete_for_propagate_difference_error"] is False
    zeros = copy.deepcopy(propagation)
    for name in ("reference_initial", "reference_residual", "M_integral", "Mz_integral"):
        assert zeros[name]["value"] is None


def test_synthetic_complete_stream_marks_all_cells_but_not_family_run(
        small_cone, tmp_path):
    mapping, result, checkpoint = run_campaign(
        small_cone, small_cone["complete"], tmp_path)
    assert mapping["schema"] == SCHEMA
    assert result.all_history_cells is True
    assert mapping["coverage"]["stream_reached_rho1"] is True
    assert mapping["coverage"]["matching_enclosed_cells"] is True
    assert mapping["coverage"]["enclosed_cells"] == 2
    assert mapping["family_14_1_scientific_run"] is False
    assert mapping["completed_family_record"] is False
    assert mapping["physical_rho1_source_error"] is None
    assert mapping["propagation"] is not None
    assert mapping["continuous_inputs"]["reference_residual_enclosed_cells"] == 2
    assert mapping["continuous_inputs"]["complete_for_propagate_difference_error"]
    assert mapping["status"].startswith("OPEN")
    loaded = WholeConeCheckpoint.from_mapping(json.loads(checkpoint.read_text()))
    loaded.validate(small_cone["config"])
    assert loaded.completed_cells == 2
    assert float.fromhex(loaded.last_rho_end) == 1.0


def test_diagnostic_stop_keeps_checkpoint_and_cannot_complete_family(
        small_cone, tmp_path):
    mapping, result, checkpoint = run_campaign(
        small_cone, small_cone["segments"], tmp_path, max_new_cells=1,
        checkpoint_interval=1)
    assert mapping["stop_reason"] == "diagnostic"
    assert mapping["coverage"]["diagnostic"] is True
    assert mapping["all_history_cells"] is False
    assert mapping["family_14_1_scientific_run"] is False
    assert mapping["completed_family_record"] is False
    assert mapping["coverage"]["enclosed_cells"] == 1
    assert mapping["forecast"]["scientific"] is False
    assert "runtime" not in scientific_record(mapping)
    assert "forecast" not in scientific_record(mapping)
    assert scientific_digest(mapping) == mapping["scientific_digest"]
    loaded = WholeConeCheckpoint.from_mapping(json.loads(checkpoint.read_text()))
    assert loaded.completed_cells == 1
    assert result.status.startswith("OPEN")
    assert mapping["propagation"] is None


def test_complete_time_stream_contracts_both_selected_energy_sectors(small_cone,tmp_path):
    from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
    from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
    positive=FixedSourcePreparation(np.array([[.6]],complex),np.array([1.]),np.array([.7]))
    negative=FixedSourcePreparation(np.array([[.4]],complex),np.array([1.]),np.array([-.7]))
    selected={'source':positive,'negative':negative,
              'channel':{'copy_count':2,'degeneracy':3,'angular_eigenvalue':.8}}
    mapping,_,_=run_campaign(small_cone,small_cone['complete'],tmp_path,
                            continuous_source=selected)
    assert mapping['continuous_inputs']['complete_for_difference_matter_error']
    assert mapping['matter']['covered_energy_sectors']==[1,-1]
    assert float(restored_upper(mapping['matter']['N'])) > 0
    assert mapping['matter']['positive']['N'] != mapping['matter']['negative']['N']
    assert mapping['physical_rho1_source_error'] is None
    assert not mapping['family_14_1_scientific_run']
    assert not mapping['certificate_use']


def test_signed_endpoint_rejects_changed_negative_energy(small_cone,tmp_path):
    from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
    selected={
        'source':FixedSourcePreparation(np.array([[.6]],complex),np.array([1.]),np.array([.7])),
        'negative':FixedSourcePreparation(np.array([[.4]],complex),np.array([1.]),np.array([-.8])),
        'channel':{'copy_count':2,'degeneracy':3,'angular_eigenvalue':.8}}
    with pytest.raises(ValueError,match='negative source energy/weight binding'):
        run_campaign(small_cone,small_cone['complete'],tmp_path,continuous_source=selected)


def test_resume_replays_prefix_identities_before_new_cells(small_cone, tmp_path):
    first, _result, checkpoint = run_campaign(
        small_cone, small_cone["segments"], tmp_path, max_new_cells=1,
        checkpoint_interval=1)
    assert first["coverage"]["enclosed_cells"] == 1
    mapping, result, _ = run_campaign(
        small_cone, small_cone["segments"], tmp_path, resume=True,
        checkpoint_path=checkpoint, max_new_cells=2, checkpoint_interval=1)
    oneshot, _, _ = run_campaign(
        small_cone, small_cone["segments"], tmp_path / "oneshot",
        checkpoint_interval=1)
    assert mapping["coverage"]["enclosed_cells"] == 3
    assert mapping["bounds"]["residual_integrals"] == oneshot["bounds"]["residual_integrals"]
    assert mapping["continuous_inputs"]["reference_residual_enclosed_cells"] == 3
    assert mapping["prefix_digest"] == oneshot["prefix_digest"]
    assert mapping["family_14_1_scientific_run"] is False
    assert result.all_history_cells is False
    resumed = WholeConeAccumulator.from_checkpoint(
        json.loads(checkpoint.read_text()), small_cone["config"],
        small_cone["payload"], small_cone["family"])
    with pytest.raises(ValueError, match="prefix replay required"):
        resumed.accumulate_segment(small_cone["segments"][1], cell_index=1)


def test_resume_rejects_changed_dense_segment_identity(small_cone, tmp_path):
    _mapping, _result, checkpoint = run_campaign(
        small_cone, small_cone["segments"][:1], tmp_path, checkpoint_interval=1)
    broken = (
        synthetic_segment(
            small_cone["segments"][0].rho_start,
            small_cone["segments"][0].rho_end, scale=1.0001),
        small_cone["segments"][1],
    )
    with pytest.raises(ValueError, match="segment replay"):
        run_campaign(
            small_cone, broken, tmp_path, resume=True,
            checkpoint_path=checkpoint, checkpoint_interval=1)


def test_checkpoint_after_each_eight_new_cells(small_cone, tmp_path, monkeypatch):
    nodes = rho_ladder(1.03, 11)
    segments = tuple(
        synthetic_segment(nodes[i], nodes[i + 1], scale=1.0 + 0.01 * i)
        for i in range(11))
    sizes = []
    original = write_checkpoint_atomic

    def watching(path, checkpoint):
        result = original(path, checkpoint)
        sizes.append(result.completed_cells)
        return result

    monkeypatch.setattr(
        "recursive_horizons.nsc_ks_current_field_campaign.write_checkpoint_atomic",
        watching)
    mapping, _result, checkpoint = run_campaign(
        small_cone, segments, tmp_path, checkpoint_interval=8, max_new_cells=10)
    assert 8 in sizes
    assert sizes[-1] == 10
    assert mapping["coverage"]["new_cells_this_session"] == 10
    loaded = WholeConeCheckpoint.from_mapping(json.loads(checkpoint.read_text()))
    assert loaded.completed_cells == 10
    assert mapping["stop_reason"] == "diagnostic"
    assert mapping["completed_family_record"] is False


def test_resource_and_method_stops_remain_open(small_cone, tmp_path):
    probes = {"n": 0}

    def space():
        probes["n"] += 1
        return 100 if probes["n"] == 1 else 0

    mapping, _result, checkpoint = run_campaign(
        small_cone, small_cone["segments"], tmp_path,
        free_space_floor=10, free_space=space, checkpoint_interval=1)
    assert mapping["stop_reason"] == "resource"
    assert mapping["coverage"]["resource_stopped"] is True
    assert mapping["family_14_1_scientific_run"] is False
    assert mapping["status"].startswith("OPEN")
    assert WholeConeCheckpoint.from_mapping(
        json.loads(checkpoint.read_text())).completed_cells == 1

    def failing(on_segment):
        on_segment(small_cone["segments"][0])
        raise ArithmeticError("method failed")

    mapping, _result = run_whole_cone_campaign(
        small_cone["config"], small_cone["payload"], small_cone["family"],
        failing, checkpoint_path=tmp_path / "method.checkpoint.json",
        stream_kind=STREAM_KIND_INJECTED, rho_up=1.03)
    checkpoint = tmp_path / "method.checkpoint.json"
    assert mapping["stop_reason"] == "method"
    assert mapping["coverage"]["method_failed"] is True
    assert mapping["family_14_1_scientific_run"] is False
    assert mapping["physical_EXISTENCE_certificate"] is False


def test_family_run_flag_requires_full_stream_and_production_bindings():
    complete = scientific_coverage(
        family=PRODUCTION_FAMILY,
        profile_identity=PRODUCTION_PROFILE_IDENTITY,
        completed_cells=4, accepted_segments=4, stream_reached_rho1=True,
        last_rho_end=1.0, cell_indices=(0, 1, 2, 3), diagnostic=False,
        resource_stopped=False, method_failed=False, interrupted=False,
        source_via_select_source=True, stream_kind=STREAM_KIND_TRAJECTORY,
        settings={
            "integrator": "dop853",
            "rtol": PRODUCTION_RTOL.hex(),
            "atol": PRODUCTION_ATOL.hex(),
            "max_step": PRODUCTION_MAX_STEP.hex(),
            "step_control": "joint",
            "tangents": "zero",
            "profile_payload_recomputed": False,
        })
    assert complete["all_history_cells"] is True
    assert complete["selected_source_trajectory_complete"] is True
    assert complete["family_14_1_scientific_run"] is False
    assert complete["full_family_source_coverage"] is False
    partial = scientific_coverage(
        family=PRODUCTION_FAMILY,
        profile_identity=PRODUCTION_PROFILE_IDENTITY,
        completed_cells=3, accepted_segments=4, stream_reached_rho1=True,
        last_rho_end=1.0, cell_indices=(0, 1, 2), diagnostic=False,
        resource_stopped=False, method_failed=False, interrupted=False,
        source_via_select_source=True, stream_kind=STREAM_KIND_TRAJECTORY,
        settings=complete and {
            "integrator": "dop853",
            "rtol": PRODUCTION_RTOL.hex(),
            "atol": PRODUCTION_ATOL.hex(),
            "max_step": PRODUCTION_MAX_STEP.hex(),
            "step_control": "joint",
            "tangents": "zero",
            "profile_payload_recomputed": False,
        })
    assert partial["all_history_cells"] is False
    assert partial["family_14_1_scientific_run"] is False
    diagnostic = scientific_coverage(
        family=PRODUCTION_FAMILY,
        profile_identity=PRODUCTION_PROFILE_IDENTITY,
        completed_cells=4, accepted_segments=4, stream_reached_rho1=True,
        last_rho_end=1.0, cell_indices=(0, 1, 2, 3), diagnostic=True,
        resource_stopped=False, method_failed=False, interrupted=False,
        source_via_select_source=True, stream_kind=STREAM_KIND_TRAJECTORY,
        settings={
            "integrator": "dop853",
            "rtol": PRODUCTION_RTOL.hex(),
            "atol": PRODUCTION_ATOL.hex(),
            "max_step": PRODUCTION_MAX_STEP.hex(),
            "step_control": "joint",
            "tangents": "zero",
            "profile_payload_recomputed": False,
        })
    assert diagnostic["all_history_cells"] is True
    assert diagnostic["family_14_1_scientific_run"] is False
    injected = scientific_coverage(
        family=PRODUCTION_FAMILY,
        profile_identity=PRODUCTION_PROFILE_IDENTITY,
        completed_cells=4, accepted_segments=4, stream_reached_rho1=True,
        last_rho_end=1.0, cell_indices=(0, 1, 2, 3), diagnostic=False,
        resource_stopped=False, method_failed=False, interrupted=False,
        source_via_select_source=True, stream_kind=STREAM_KIND_INJECTED,
        settings={
            "integrator": "dop853",
            "rtol": PRODUCTION_RTOL.hex(),
            "atol": PRODUCTION_ATOL.hex(),
            "max_step": PRODUCTION_MAX_STEP.hex(),
            "step_control": "joint",
            "tangents": "zero",
            "profile_payload_recomputed": False,
        })
    assert injected["family_14_1_scientific_run"] is False
    assert solver_matches_production({
        "integrator": "dop853",
        "rtol": PRODUCTION_RTOL.hex(),
        "atol": PRODUCTION_ATOL.hex(),
        "max_step": PRODUCTION_MAX_STEP.hex(),
        "step_control": "joint",
        "tangents": "zero",
        "profile_payload_recomputed": False,
    })


def test_runtime_forecast_is_nonscientific(small_cone):
    forecast = runtime_forecast(
        completed_cells=8, rho_start=1.03, last_rho_end=1.02, cpu_seconds=16.0)
    assert forecast["scientific"] is False
    assert forecast["estimated_remaining_cells"] > 0
    record = {
        "schema": SCHEMA,
        "status": "OPEN: prefix",
        "family": [14, 1],
        "runtime": {"CPU_seconds": 16.0},
        "forecast": forecast,
        "memory": {"peak": 1},
    }
    assert "runtime" not in scientific_record(record)
    assert "forecast" not in scientific_record(record)
    assert "memory" not in scientific_record(record)


def test_mismatched_source_or_settings_reject_resume(small_cone, tmp_path):
    _mapping, _result, checkpoint = run_campaign(
        small_cone, small_cone["segments"][:1], tmp_path, checkpoint_interval=1)
    other = dict(small_cone["config"].settings)
    other["retained_index"] = 8
    with pytest.raises(ValueError, match="mismatched settings"):
        WholeConeFieldConfig(
            family=small_cone["config"].family,
            profile_identity=small_cone["config"].profile_identity,
            bits=small_cone["config"].bits,
            time_degree=small_cone["config"].time_degree,
            reciprocal_order=small_cone["config"].reciprocal_order,
            spatial_count=small_cone["config"].spatial_count,
            period_origin=small_cone["config"].period_origin,
            period_length_numerator=2, period_length_denominator=5,
            mass=small_cone["config"].mass, angular=small_cone["config"].angular,
            source_weights=small_cone["config"].source_weights,
            source_energies=small_cone["config"].source_energies,
            radius_tail=small_cone["config"].radius_tail,
            profile_payload_sha256=small_cone["config"].profile_payload_sha256,
            settings=other,
            settings_digest=small_cone["config"].settings_digest,
            source_binding_sha256=small_cone["config"].source_binding_sha256)
    mapping = json.loads(checkpoint.read_text())
    mapping["source_binding_sha256"] = "0" * 64
    broken = WholeConeCheckpoint.from_mapping(mapping)
    with pytest.raises(ValueError, match="mismatched source"):
        broken.validate(small_cone["config"])

@pytest.mark.parametrize('change', ['initial', 'covariance', 'target', 'archive'])
def test_production_resume_binding_rejects_changed_inputs(small_cone, change):
    sys.path.insert(0, str(ROOT / 'scripts'))
    from run_nsc_ks_whole_cone_field_v1 import bind_selected_inputs
    selected = {
        'source': SimpleNamespace(digest='a' * 64),
        'negative': SimpleNamespace(digest='b' * 64),
        'initial_columns': np.array([[.5], [.25]], complex),
        'archive_hashes': {'source.json': 'c' * 64},
        'selections': [{'row': 0}], 'channel': {'mass': .4},
    }
    grid, target = np.arange(8, dtype=float), np.array([.12, .18])
    first = bind_selected_inputs(small_cone['config'], selected, grid, target, ROOT)
    accumulator = WholeConeAccumulator(first, small_cone['payload'], small_cone['family'])
    accumulator.accumulate_segment(small_cone['segments'][0], cell_index=0)
    checkpoint = accumulator.checkpoint()
    if change == 'initial': selected['initial_columns'][0, 0] += .01
    if change == 'covariance': selected['source'] = SimpleNamespace(digest='d' * 64)
    if change == 'target': target[0] += .01
    if change == 'archive': selected['archive_hashes'] = {'source.json': 'e' * 64}
    second = bind_selected_inputs(small_cone['config'], selected, grid, target, ROOT)
    assert first.digest() != second.digest()
    with pytest.raises(ValueError): checkpoint.validate(second)
