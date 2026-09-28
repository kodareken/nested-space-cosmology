"""Higher numerical order keeps the same equation and explicit error owner."""
import json
from pathlib import Path

import mpmath as mp
import pytest

from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi, _geometry_jets, _defect_coefficients
from recursive_horizons.nsc_incoming_projector_energy_bound import (
    incoming_riccati_interval_coefficients, local_cross_product_coefficients, projected_energy_bound,
)
from recursive_horizons.nsc_incoming_source_tail import _riccati_at_one
from recursive_horizons.nsc_incoming_higher_order_defect import HigherOrderDefectGeometry

ROOT = Path(__file__).resolve().parents[1]
CHANNEL = json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels'][14]


def test_order24_point_coefficients_and_leading_defect_match_owned_recurrence():
    with _precision(40):
        c = incoming_riccati_interval_coefficients(CHANNEL, order=25)
        values = _riccati_at_one(mp.mpf(CHANNEL['compact_mass']), mp.mpf(CHANNEL['angular_eigenvalue']), order=25)
        for interval, value in zip(c, values):
            assert _lo(interval.real) <= value.real <= _hi(interval.real)
            assert _lo(interval.imag) <= value.imag <= _hi(interval.imag)
        mass, angular = mp.iv.mpf(CHANNEL['compact_mass']), mp.iv.mpf(CHANNEL['angular_eigenvalue'])
        residual = _defect_coefficients(*_geometry_jets(mp.iv.mpf(1), 24), mass, angular, order=24)
        assert len(residual) == 25
        # Truncating S at c24 leaves f24=-c25, fixing both sign and order.
        leading = -values[24]
        assert _lo(residual[0].real) <= leading.real <= _hi(residual[0].real)
        assert _lo(residual[0].imag) <= leading.imag <= _hi(residual[0].imag)


def test_order24_energy_bound_rejects_order16_radial_evidence():
    old = json.loads((ROOT/'results/development/nsc-incoming-vacuum-tail-bound.json').read_text())['groups'][13]
    cross = local_cross_product_coefficients(CHANNEL, order=24)
    assert len(cross['cross_numerator_coefficient_upper']) == 48
    with pytest.raises(ValueError, match='matching explicit numerical order'):
        projected_energy_bound(CHANNEL, old, cross, 16., 160.)


def test_higher_order_is_numerical_and_keeps_bounded_domain():
    config = json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    owner = HigherOrderDefectGeometry(config, intervals=32)
    assert owner.order == 24
    assert len(owner.cells) == 32
    assert len(owner.cells[0][1]) == 26
    with pytest.raises(ValueError, match='group14'):
        owner.bound_channel({**CHANNEL, 'index': 13}, 160.)
    with pytest.raises(ValueError, match='order24'):
        HigherOrderDefectGeometry(config, order=32)
