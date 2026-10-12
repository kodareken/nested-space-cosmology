"""Finite controls must retain the original coherent and signed source rows."""
from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np
import pytest
from scipy.linalg import block_diag

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import derive_nsc_ks_source_control_v2 as C


def archive(cross=False, wrong_weight=False):
    covariance = block_diag(np.eye(3) * 0.3, np.eye(3) * 0.6).astype(complex)
    covariance[0, 1], covariance[1, 0] = 0.02j, -0.02j
    if cross:
        covariance[0, 3] = covariance[3, 0] = 0.001
    energies, weights = np.repeat([0.5, 1.5], 3), np.repeat([0.2, 0.4], 3)
    source = C.FixedSourcePreparation(covariance, weights, energies)
    neg_weights = weights * 2 if wrong_weight else weights
    negative = C.FixedSourcePreparation(np.eye(6)-covariance.conj(), neg_weights, -energies)
    initial = np.arange(12).reshape(2, 6).astype(complex) / 20
    common = dict(mass=1.0, angular=2.0, rho_up=1.03)
    plus = SimpleNamespace(source=source, initial_columns=initial, panel_name="panel",
                           energy_sign=1, **common)
    minus = SimpleNamespace(source=negative, initial_columns=C.s3_conjugate(initial),
                            panel_name="panel/negative-partner", energy_sign=-1, **common)
    entries = [(plus, {"angular_sign": 1}), (minus, {"angular_sign": -1})]
    return SimpleNamespace(family_entries=lambda _: entries), plus


def test_source_control_preserves_columns_covariance_weights_and_both_signs():
    data, original = archive()
    positive, negative, initial, _batch, _channel, rows = C.select_source(data, (14, 1), count=2)
    np.testing.assert_array_equal(positive.covariance, original.source.covariance)
    np.testing.assert_array_equal(positive.column_weights, original.source.column_weights)
    np.testing.assert_array_equal(initial, original.initial_columns)
    np.testing.assert_array_equal(negative.energies, -positive.energies)
    assert [row["energy"] for row in rows] == [0.5, 1.5]
    assert positive.covariance[0, 1] == 0.02j


def test_control_refuses_discarding_cross_energy_coherence():
    data, _ = archive(cross=True)
    with pytest.raises(ValueError, match="cross-energy coherence"):
        C.select_source(data, (14, 1), count=1)


def test_control_refuses_altering_negative_sector_weights():
    data, _ = archive(wrong_weight=True)
    with pytest.raises(ValueError, match="signed preparation relation"):
        C.select_source(data, (14, 1), count=2)
