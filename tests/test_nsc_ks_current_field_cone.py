"""Whole-cone field core: payload restore, streamed totals, checkpoint resume."""
import copy
import gc
import json
from hashlib import sha256
from pathlib import Path
import sys
import weakref
from fractions import Fraction

import numpy as np
import pytest

pytest.importorskip("flint")
from flint import arb, ctx

from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_current_field_cone import (
    CHECKPOINT_SCHEMA, SCHEMA, V4_PROFILE_IDENTITY,
    WholeConeAccumulator, WholeConeCheckpoint, WholeConeFieldConfig,
    accumulate_segment, finalize_family, packed_radius_tail,
    restore_v4_profile_payload, scientific_digest, scientific_record,
)
from recursive_horizons.nsc_ks_difference_residual_polynomial import monomial_keys
from recursive_horizons.nsc_ks_local_profile_fourier_bound import (
    enclose_profile_fourier_local, profile_payload_digest,
    restore_local_profile_fourier, serialize_local_profile_fourier,
)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_trajectory import TrajectorySegment
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


ROOT = Path(__file__).resolve().parents[1]
V4_PATH = ROOT / "results/development/nsc-ks-current-field-pilot-v4.json"
EXPECTED_IDENTITY = V4_PROFILE_IDENTITY


def test_checkpoint_cell_dyadic_does_not_round_down(small_cone, monkeypatch):
    from recursive_horizons import nsc_ks_current_field_cone as owner
    family, config, payload = (small_cone[name] for name in ('family', 'config', 'payload'))
    with ctx.workprec(config.bits):
        value = arb(1) + arb(2) ** -75
    enclosed = {name: (value, value) for name in ('polynomial', 'remainder', 'total', 'integral')}
    enclosed['width'] = arb(1)
    monkeypatch.setattr(owner, 'enclose_difference_segment', lambda *args: enclosed)
    accumulator = WholeConeAccumulator(config, payload, family)
    with ctx.workprec(53):
        record = accumulator.accumulate_segment(synthetic_segment(1.02, 1.01), cell_index=0)
    bound = record['total'][0]
    represented = Fraction(int(bound['mantissa'])) * Fraction(2) ** bound['exponent']
    assert represented >= 1 + Fraction(1, 2**75)


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


@pytest.fixture(scope="module")
def small_cone():
    family = small_family()
    origin = family.center - 0.2
    length = 0.4
    bits = 80
    order = 1
    keys = monomial_keys(order)
    rows = enclose_profile_fourier_local(
        family, keys, origin, length, retained_index=4, taylor_order=5,
        exterior_panels=16, transition_panels=16, interior_panels=16,
        max_derivative=2, flat_denominator=32, relative_tolerance_bits=20,
        absolute_tolerance_bits=30, max_depth=8, bits=bits)
    payload = serialize_local_profile_fourier(rows)
    digest = profile_payload_digest(payload)
    restored = restore_local_profile_fourier(payload, expected_sha256=digest)
    settings = {
        "retained_index": 4,
        "taylor_order": 5,
        "exterior_panels": 16,
        "transition_panels": 16,
        "interior_panels": 16,
        "bits": bits,
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
    step = 1.0 / 8192
    rho = [1.03]
    for _ in range(3):
        rho.append(rho[-1] - step)
    segments = (
        synthetic_segment(rho[0], rho[1]),
        synthetic_segment(rho[1], rho[2], scale=1.1),
        synthetic_segment(rho[2], rho[3], scale=0.9),
    )
    return {
        "family": family, "rows": restored, "payload": payload,
        "digest": digest, "config": config, "segments": segments,
    }


def run_accumulator(cone, segments, *, start=0):
    acc = WholeConeAccumulator(cone["config"], cone["payload"], cone["family"])
    for offset, segment in enumerate(segments):
        accumulate_segment(acc, segment, cell_index=start + offset)
    return acc


def packed_totals(result):
    return (
        result.bounds["polynomial_sum"],
        result.bounds["remainder_sum"],
        result.bounds["total_continuous_normalized_residual_sum"],
        result.bounds["residual_integrals"],
    )


def test_v4_profile_payload_round_trip_and_malformed_hash_rejection():
    raw = V4_PATH.read_bytes()
    record = json.loads(raw)
    profiles, digest = restore_v4_profile_payload(
        record, expected_payload_sha256=None, expected_file_sha256=sha256(raw).hexdigest(),
        file_bytes=raw)
    payload = record["profile_coefficient_payload"]
    assert profile_payload_digest(payload) == digest
    restored = restore_local_profile_fourier(payload, expected_sha256=digest)
    assert set(restored) == set(profiles)
    key = (0, 1)
    for original, replayed in zip(profiles[key].coefficients[:8],
                                  restored[key].coefficients[:8]):
        assert replayed.real.contains(original.real)
        assert replayed.imag.contains(original.imag)
    replayed_payload = serialize_local_profile_fourier(restored)
    again = restore_local_profile_fourier(replayed_payload)
    for original, replayed in zip(restored[key].coefficients[:8],
                                  again[key].coefficients[:8]):
        assert replayed.real.contains(original.real)
        assert replayed.imag.contains(original.imag)
    with pytest.raises(ValueError, match="mismatched profile payload"):
        restore_local_profile_fourier(payload, expected_sha256="0" * 64)
    tampered = copy.deepcopy(payload)
    tampered["0_1"]["retained_index"] = 1
    with pytest.raises(ValueError, match="mismatched profile payload"):
        restore_v4_profile_payload(record, expected_payload_sha256="0" * 64)
    with pytest.raises(ValueError, match="mismatched profile payload"):
        restore_local_profile_fourier(tampered, expected_sha256=digest)
    with pytest.raises(ValueError, match="mismatched profile payload"):
        restore_v4_profile_payload(
            record, expected_file_sha256="0" * 64, file_bytes=raw)


def test_deterministic_accumulator_and_dense_streamed_small_slab(small_cone):
    segments = small_cone["segments"]
    first = finalize_family(run_accumulator(small_cone, segments))
    second = finalize_family(run_accumulator(small_cone, segments))
    assert packed_totals(first) == packed_totals(second)
    assert first.prefix_digest == second.prefix_digest
    assert first.schema == SCHEMA
    assert first.history_minus_reference_before_norms
    streamed = WholeConeAccumulator(
        small_cone["config"], small_cone["payload"], small_cone["family"])
    for index, segment in enumerate(segments):
        accumulate_segment(streamed, segment, cell_index=index)
        del segment
    streamed_result = finalize_family(streamed)
    assert packed_totals(streamed_result) == packed_totals(first)
    assert first.coverage["all_history_cells"] is False
    assert first.family_14_1_scientific_run is False
    assert first.physical_rho1_source_error is None
    # The stored polynomial bounds are for S=X_xi-hLX.  They are already the
    # per-cell Gronwall integral contribution and must not receive another h.
    first_cell = run_accumulator(small_cone, segments[:1]).checkpoint().cells[0]
    assert first_cell["integral"] == first_cell["total"]


def test_checkpoint_prefix_resume_and_corruption_rejection(small_cone, tmp_path):
    segments = small_cone["segments"]
    oneshot = finalize_family(run_accumulator(small_cone, segments))
    prefix = WholeConeAccumulator(
        small_cone["config"], small_cone["payload"], small_cone["family"])
    accumulate_segment(prefix, segments[0], cell_index=0)
    checkpoint = prefix.checkpoint()
    assert checkpoint.schema == CHECKPOINT_SCHEMA
    assert checkpoint.completed_cells == 1
    path = tmp_path / "cone.checkpoint.json"
    path.write_text(json.dumps(checkpoint.to_mapping(), indent=2, sort_keys=True) + "\n")
    loaded = WholeConeCheckpoint.from_mapping(json.loads(path.read_text()))
    loaded.validate(small_cone["config"])
    resumed = WholeConeAccumulator.from_checkpoint(
        loaded, small_cone["config"], small_cone["payload"], small_cone["family"])
    with pytest.raises(ValueError, match="prefix replay required"):
        accumulate_segment(resumed, segments[1], cell_index=1)
    resumed.replay_segment(segments[0], cell_index=0)
    accumulate_segment(resumed, segments[1], cell_index=1)
    accumulate_segment(resumed, segments[2], cell_index=2)
    resumed_result = finalize_family(resumed)
    assert packed_totals(resumed_result) == packed_totals(oneshot)
    assert resumed_result.prefix_digest == oneshot.prefix_digest
    with pytest.raises(ValueError, match="segment replay"):
        resumed.replay_segment(segments[1], cell_index=0)
    mapping = checkpoint.to_mapping()
    mapping["prefix_digest"] = "0" * 64
    with pytest.raises(ValueError, match="mismatched prefix digest"):
        WholeConeCheckpoint.from_mapping(mapping)
    mapping = checkpoint.to_mapping()
    mapping["cells"][0]["total"][0]["mantissa"] = "1"
    with pytest.raises(ValueError, match="mismatched prefix digest"):
        WholeConeCheckpoint.from_mapping(mapping)
    mapping = checkpoint.to_mapping()
    mapping["residual_integrals"][0]["mantissa"] = "1"
    with pytest.raises(ValueError, match="normalized residual integrals"):
        WholeConeCheckpoint.from_mapping(mapping)
    mapping = checkpoint.to_mapping()
    mapping["polynomial_sum"][0]["mantissa"] = "1"
    with pytest.raises(ValueError, match="mismatched prefix digest"):
        WholeConeCheckpoint.from_mapping(mapping)
    mapping = checkpoint.to_mapping()
    mapping["unexpected"] = True
    with pytest.raises(ValueError, match="malformed whole-cone checkpoint"):
        WholeConeCheckpoint.from_mapping(mapping)
    resumed = WholeConeAccumulator.from_checkpoint(
        checkpoint, small_cone["config"], small_cone["payload"], small_cone["family"])
    with pytest.raises(ValueError, match="segment replay"):
        resumed.replay_segment(
            synthetic_segment(
                segments[0].rho_start, segments[0].rho_end, scale=1.0001),
            cell_index=0)
    mapping = checkpoint.to_mapping()
    mapping["profile_payload_sha256"] = "0" * 64
    broken = WholeConeCheckpoint.from_mapping(mapping)
    with pytest.raises(ValueError, match="mismatched profile payload"):
        broken.validate(small_cone["config"])
    mapping = checkpoint.to_mapping()
    mapping["source_binding_sha256"] = "0" * 64
    broken = WholeConeCheckpoint.from_mapping(mapping)
    with pytest.raises(ValueError, match="mismatched source"):
        broken.validate(small_cone["config"])
    mapping = checkpoint.to_mapping()
    mapping["settings_digest"] = "0" * 64
    broken = WholeConeCheckpoint.from_mapping(mapping)
    with pytest.raises(ValueError, match="mismatched settings"):
        broken.validate(small_cone["config"])
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


def test_accumulator_does_not_retain_segments(small_cone):
    acc = WholeConeAccumulator(
        small_cone["config"], small_cone["payload"], small_cone["family"])
    refs = []
    step = 1.0 / 8192
    rho = [1.03]
    for _ in range(3):
        rho.append(rho[-1] - step)
    for index, (rho_start, rho_end, scale) in enumerate(
            ((rho[0], rho[1], 1.0), (rho[1], rho[2], 1.1), (rho[2], rho[3], 0.9))):
        segment = synthetic_segment(rho_start, rho_end, scale=scale)
        refs.append(weakref.ref(segment))
        accumulate_segment(acc, segment, cell_index=index)
        del segment
    gc.collect()
    assert all(reference() is None for reference in refs)
    assert not hasattr(acc, "segments")
    result = finalize_family(acc)
    assert result.completed_cells == 3


def test_adapters_keep_rho1_source_error_unresolved(small_cone):
    acc = run_accumulator(small_cone, small_cone["segments"])
    result = finalize_family(
        acc,
        propagation={
            "reference_initial": 0,
            "reference_residual": 0,
            "reference_offdiagonal_integral": 0,
            "difference_initial": (0, 0),
            "offdiagonal_integral": 0,
            "Bz_integral": 0,
            "M_integral": 0,
            "Mz_integral": 0,
            "maximum_absolute_energy": 1.0,
        },
        matter={
            "reference_norm": 1, "difference_norm": 1,
            "reference_axial_norm": 1, "difference_axial_norm": 1,
            "reference_error": 1e-6, "difference_error": 1e-6,
            "reference_axial_error": 1e-6, "difference_axial_error": 1e-6,
            "source_norm": 2, "mass": 0.4, "absolute_angular": 0.8,
            "axial_lower": 0.8, "radius_lower": 1.3, "multiplicity": 1,
        },
    )
    assert result.propagation["source_accuracy_included"] is False
    assert result.matter["source_accuracy_included"] is False
    assert result.propagation["physical_rho1_source_error"] is None
    assert result.physical_EXISTENCE_certificate is False
    assert result.physical_NONEXISTENCE_certificate is False
    assert result.status.startswith("OPEN")
    empty = WholeConeAccumulator(
        small_cone["config"], small_cone["payload"], small_cone["family"])
    with pytest.raises(ValueError, match="at least one accumulated segment"):
        finalize_family(empty)
    acc = run_accumulator(small_cone, small_cone["segments"][:1])
    with pytest.raises(ValueError, match="physical rho=1 source error"):
        finalize_family(acc, source_accuracy_included=True)
    acc = run_accumulator(small_cone, small_cone["segments"][:1])
    with pytest.raises(ValueError, match="physical rho=1 source error"):
        finalize_family(acc, source_error=0)
    record = result.to_mapping()
    record["runtime"] = {"CPU_seconds": 1.25}
    record["memory"] = {"peak_bytes": 12}
    assert "runtime" not in scientific_record(record)
    assert scientific_digest(record) == scientific_digest(result.to_mapping())


def test_v4_cell_122_restored_payload_stays_inside_directed_enclosure():
    sys.path.insert(0, str(ROOT / "scripts"))
    import validate_nsc_ks_current_field_pilot_v3 as V3
    raw = V4_PATH.read_bytes()
    recorded = json.loads(raw)
    _profiles, digest = restore_v4_profile_payload(
        recorded, expected_file_sha256=sha256(raw).hexdigest(), file_bytes=raw)
    capture, arrays = V3.load_capture()
    family = V3.load_family(capture)
    _bounds, _rejected, radius_tail = V3.geometry_bounds(family, capture)
    segment, _grid, origin, _period, period_q, weights, energies = V3.segment_from_capture(
        capture, arrays)
    assert capture["cell_index"] == 122
    assert capture["family"] == [14, 1]
    assert capture["profile_identity"] == EXPECTED_IDENTITY
    with ctx.workprec(90):
        tail = packed_radius_tail(family, capture["rho_up"], capture["angular"],
                                  V3.RECIPROCAL_ORDER, bits=90)
        for packed, original in zip(
                tail, radius_tail["potential_derivative_tail_bounds"][:2]):
            assert restored_upper(packed) >= arb(original.numerator) / arb(
                original.denominator)
    config = WholeConeFieldConfig(
        family=(14, 1), profile_identity=EXPECTED_IDENTITY, bits=90,
        time_degree=V3.TIME_DEGREE, reciprocal_order=V3.RECIPROCAL_ORDER,
        spatial_count=capture["grid_nodes"], period_origin=origin,
        period_length_numerator=period_q.numerator,
        period_length_denominator=period_q.denominator,
        mass=capture["mass"], angular=capture["angular"],
        source_weights=tuple(map(float, weights)),
        source_energies=tuple(map(float, energies)),
        radius_tail=tail, profile_payload_sha256=digest,
        settings=dict(recorded["settings"]))
    acc = WholeConeAccumulator(
        config, recorded["profile_coefficient_payload"], family)
    accumulate_segment(acc, segment, cell_index=122)
    result = finalize_family(acc)
    assert result.completed_cells == 1
    assert result.coverage["all_history_cells"] is False
    assert result.family_14_1_scientific_run is False
    recorded_total = [restored_upper(value) for value in
                      recorded["bounds"]["total_continuous_normalized_residual"]]
    computed = [restored_upper(value) for value in result.bounds["last_cell_total"]]
    for bound, previous in zip(computed, recorded_total):
        assert bound >= previous
        assert bound.contains(previous) or float(bound / previous) < 1 + 1e-9
    assert result.physical_EXISTENCE_certificate is False
    assert result.physical_NONEXISTENCE_certificate is False
    assert result.physical_rho1_source_error is None
