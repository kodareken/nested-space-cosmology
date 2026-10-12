"""Directed local Fourier bounds for compact incoming-history profiles."""
import json

import numpy as np
import pytest

pytest.importorskip("flint")

from recursive_horizons.nsc_ks_local_profile_fourier_bound import (
    LocalProfileFourierBound,
    enclose_profile_fourier_local,
    profile_payload_digest,
    restore_local_profile_fourier,
    serialize_local_profile_fourier,
)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


def small_family():
    coefficients = np.zeros((2, 8))
    coefficients[0, :4] = (0.012, -0.003, 0.0015, -0.0004)
    coefficients[1, :4] = (0.021, 0.004, -0.001, 0.0003)
    return LocalIncomingFamily(coefficients)


def test_local_tail_encloses_dense_value_and_derivative_samples():
    family = small_family()
    key = (1, 1)
    origin = family.center - 0.2
    length = 0.4
    rows = enclose_profile_fourier_local(
        family, (key,), origin, length, retained_index=8,
        taylor_order=5, exterior_panels=16, transition_panels=16,
        interior_panels=16, max_derivative=1, bits=80)
    row = rows[key]
    coefficients = np.array([
        complex(float(value.real.mid()), float(value.imag.mid()))
        for value in row.coefficients
    ])
    modes = np.arange(-row.retained_index, row.retained_index + 1)
    points = np.linspace(origin, origin + length, 1001, endpoint=False)
    phase = np.exp(2j * np.pi * np.outer(points - origin, modes) / length)
    partial = phase @ coefficients
    w, U = family.functions
    exact = np.array([w(z, 0) * U(z, 0) for z in points])
    exact_z = np.array([w(z, 1) * U(z, 0) + w(z, 0) * U(z, 1)
                        for z in points])
    partial_z = phase @ (2j * np.pi * modes / length * coefficients)
    assert np.max(np.abs(exact - partial.real)) <= float(row.uniform_tail_bounds[0])
    assert np.max(np.abs(exact_z - partial_z.real)) <= float(row.uniform_tail_bounds[1])
    assert row.coefficient_flat_strip_error > 0
    for mode in range(1, row.retained_index + 1):
        assert row.coefficients[row.retained_index - mode].overlaps(
            row.coefficients[row.retained_index + mode].conjugate())
    restored = restore_local_profile_fourier(serialize_local_profile_fourier(rows))[key]
    assert len(restored.coefficients) == len(row.coefficients)
    for original, replayed in zip(row.coefficients, restored.coefficients):
        assert replayed.real.contains(original.real)
        assert replayed.imag.contains(original.imag)
    for original, replayed in zip(row.uniform_tail_bounds,
                                  restored.uniform_tail_bounds):
        assert replayed.contains(original)


def test_local_profile_bound_rejects_invalid_requests():
    family = small_family()
    with pytest.raises(ValueError, match="positive-degree"):
        enclose_profile_fourier_local(family, ((0, 0),), family.center - 0.2, 0.4)
    with pytest.raises(ValueError, match="positive numerical period"):
        enclose_profile_fourier_local(family, ((1, 0),), family.center, -1.0)


def test_serialized_payload_digest_round_trip_rejects_a_mismatched_hash():
    from flint import acb, arb
    row = LocalProfileFourierBound(
        (1, 0), (acb(1), acb(0, 1), acb(1)),
        (arb("1e-3"), arb("2e-3"), arb("3e-3")),
        arb("1e-20"), 1, 3, 4, 1, arb(0), arb("2/5"))
    payload = serialize_local_profile_fourier({(1, 0): row})
    digest = profile_payload_digest(payload)
    assert profile_payload_digest(payload) == digest
    restored = restore_local_profile_fourier(payload, expected_sha256=digest)[(1, 0)]
    replayed = restore_local_profile_fourier(
        serialize_local_profile_fourier({(1, 0): restored}))
    for original, again in zip(restored.coefficients, replayed[(1, 0)].coefficients):
        assert again.real.contains(original.real)
        assert again.imag.contains(original.imag)
    with pytest.raises(ValueError, match="mismatched profile payload"):
        restore_local_profile_fourier(payload, expected_sha256="0" * 64)

    malformed = json.loads(json.dumps(payload))
    malformed["1_0"]["coefficients"].pop()
    with pytest.raises(ValueError, match="coefficient band"):
        restore_local_profile_fourier(malformed)

    malformed = json.loads(json.dumps(payload))
    malformed["1_0"]["uniform_tail_bounds"][0]["radius_upper"] = "-1"
    with pytest.raises(ValueError, match="finite nonnegative"):
        restore_local_profile_fourier(malformed)

    malformed = json.loads(json.dumps(payload))
    malformed["1_0"]["uniform_tail_bounds"] = []
    with pytest.raises(ValueError, match="derivative tail"):
        restore_local_profile_fourier(malformed)
    tampered = json.loads(json.dumps(payload))
    tampered["1_0"]["retained_index"] = 99
    with pytest.raises(ValueError, match="mismatched profile payload"):
        restore_local_profile_fourier(tampered, expected_sha256=digest)
