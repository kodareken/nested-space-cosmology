"""Independent cancellation of the reference/difference residual split."""
import json
from pathlib import Path
import sys

import numpy as np
import pytest
pytest.importorskip('flint')
from flint import arb, ctx

from recursive_horizons.nsc_ks_ball_operator import AnalyticRadiusFamily
from recursive_horizons.nsc_ks_current_history_bounds import radius_bounds
from recursive_horizons.nsc_ks_difference_residual_polynomial import (
    DifferenceResidualPolynomial, ball_split, difference_operator_remainder_bounds,
    difference_polynomial_bounds, monomial_keys, split_trajectory_segment,
)
from recursive_horizons.nsc_ks_fourier_residual_bound import polynomial_residual_bounds
from recursive_horizons.nsc_ks_radius_enclosure import axial_profile_bounds, reciprocal_radius_tail
from recursive_horizons.nsc_ks_residual_polynomial import ResidualPolynomial
from recursive_horizons.nsc_ks_trajectory import TrajectorySegment
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from test_nsc_ks_difference_envelope import inputs


ROOT = Path(__file__).resolve().parents[1]
CAPTURE = ROOT / 'results/development/nsc-ks-current-trajectory-pilot-v2.json'
EXPECTED_IDENTITY = '0b0e4cedfdb695d342c7f7347740631e916bdaee3a977df56808c725b25dacb0'


def _homogeneous(n=8):
    start = np.r_[np.array([.5, .25]), np.zeros(2 * n)]
    end = 2 * start
    segment = TrajectorySegment(1.03, 1.02, start, end, np.zeros((6, len(start)), complex))
    field_x, field_a, field_d = ball_split(segment, [1.], n, 1., bits=120)
    return segment, field_x, field_a, field_d


def test_fourteen_reciprocal_monomials_match_the_owned_time_series():
    keys = monomial_keys(4)
    assert len(keys) == 14
    assert keys[0] == (1, 0) and keys[1] == (0, 1) and keys[-1] == (0, 4)
    model = AnalyticRadiusFamily(inputs()[2])
    assert not model.w_zero and not model.U_zero
    with ctx.workprec(80):
        series = model.time_series(arb('1.015'), 2., 2, 4)
    assert set(k for k in series if isinstance(k, tuple)) == set(keys)


def test_zero_difference_cancels_background_against_nonzero_reference():
    _, field_x, field_a, field_d = _homogeneous()
    with ctx.workprec(120):
        combined = ResidualPolynomial(field_x, [1], [1], [arb(1) / 2], {}, .4, .8, [.7])
        reference = ResidualPolynomial(field_a, [1], [1], [arb(1) / 2], {}, .4, .8, [.7])
        residual = DifferenceResidualPolynomial(
            field_x, field_a, [1], [1], [arb(1) / 2], {}, .4, .8, [.7])
        assert residual.free_contains_zero()
        bound_x = polynomial_residual_bounds(combined, {})
        bound_a = polynomial_residual_bounds(reference, {})
        bound_d = difference_polynomial_bounds(residual, {})
        assert bound_x[0] > arb(.1) and bound_a[0] > arb(.1)
        assert bound_d[0] < arb('1e-30') and bound_d[1] < arb('1e-30')
        assert bound_d[0] < bound_x[0] / arb(1000)
        empty = {'inv_a2': arb(0), 'inv_a': arb(0), 'inv_ar': arb(0)}
        rem = difference_operator_remainder_bounds(
            field_d, field_a, empty, {}, (arb(0), arb(0)), .4, .8, [.7])
        assert rem == (arb(0), arb(0))


def test_zero_difference_retains_potential_on_the_reference_only():
    _, field_x, field_a, _ = _homogeneous()
    potential = {(1, 0): [arb(.3)]}
    residual = DifferenceResidualPolynomial(
        field_x, field_a, [1], [1], [arb(1) / 2], potential, .4, .8, [.7])
    combined = ResidualPolynomial(
        field_x, [1], [1], [arb(1) / 2], potential, .4, .8, [.7])
    with ctx.workprec(120):
        assert residual.free_contains_zero()
        for p in range(residual.degree + 1):
            for spin in range(2):
                for source in range(1):
                    for mode in range(field_x.spatial_count):
                        left = residual.terms[(1, 0)][p][spin][source][mode]
                        right = combined.terms[(1, 0)][p][spin][source][mode]
                        assert left.overlaps(right)
                        if mode != 0:
                            assert left.contains(0)


def test_nonzero_difference_keeps_a_nonvanishing_background_residual():
    n = 8
    harmonic = np.array([1, 1j, -1, -1j] * 2, complex)
    start = np.r_[np.array([.5, .25]), np.stack((harmonic, 2 * harmonic)).ravel()]
    segment = TrajectorySegment(1.03, 1.02, start, 2 * start, np.zeros((6, len(start)), complex))
    field_x, field_a, field_d = ball_split(segment, [1.], n, 1., bits=120)
    residual = DifferenceResidualPolynomial(
        field_x, field_a, [1], [1], [arb(1) / 2], {(1, 0): [arb(.3)]}, .4, .8, [.7])
    with ctx.workprec(120):
        assert not residual.free_contains_zero()
        x_seg, a_seg, d_seg = split_trajectory_segment(segment, 1, n)
        assert np.allclose(a_seg.start[:2], [.5, .25])
        assert np.allclose(a_seg.start[2:], 0)
        assert np.allclose(d_seg.start[:2], 0)
        assert np.allclose(d_seg.start[2:], start[2:])
        assert field_d.spatial_count == n
        rem = difference_operator_remainder_bounds(
            field_d, field_a,
            {'inv_a2': arb(0), 'inv_a': arb(0), 'inv_ar': arb(0), (1, 0): arb(0)},
            {(1, 0): (arb(1), arb(1))}, (arb(0), arb(0)), .4, .8, [.7])
        assert rem == (arb(0), arb(0))


def test_current_chebyshev_radius_rejects_the_power_coefficient_sum():
    history = json.loads((ROOT / 'results/development/nsc-ks-gate-history-lm-broyden.json').read_text())
    family = LocalIncomingFamily(np.array(history['history']['coefficients']))
    bound = radius_bounds(family, 1.0300000000000002)
    assert 1.3997 < float(bound['radius_lower']) < 1.3998
    old_u = axial_profile_bounds(family.functions[1])[0]
    assert float(old_u) > 1e8
    assert float(old_u) > float(bound['U']['profile_bounds'][0]) * 1e6
    tail = reciprocal_radius_tail(
        bound['w']['profile_bounds'], bound['U']['profile_bounds'],
        normal_support=bound['normal_support'],
        reference_radius_lower=bound['reference_radius_lower'],
        axial_lower=bound['axial_lower'], absolute_angular=2.23606797749979, order=4)
    assert 1.70e-4 < float(tail['radius_ratio_upper']) < 1.71e-4


def test_capture_record_declares_the_current_history_pilot_bindings():
    if not CAPTURE.exists():
        pytest.skip('parent current-trajectory capture is not present')
    record = json.loads(CAPTURE.read_text())
    assert record['schema'] == 'NSC-KS-CURRENT-TRAJECTORY-PILOT-v2'
    assert record['cell_index'] == 122
    assert record['family'] == [14, 1]
    assert record['grid_nodes'] == 1024
    assert record['period_length'] == 0.4
    assert record['source_columns'] == 12
    assert len(record['source_rows']) == 4
    assert record['profile_identity'] == EXPECTED_IDENTITY
    assert record['scope']['physical_EXISTENCE_certificate'] is False
    assert record['status'].startswith('OPEN')
    payload = ROOT / record['payload']['path']
    assert payload.exists()
    from hashlib import sha256
    assert sha256(payload.read_bytes()).hexdigest() == record['payload']['sha256']


def test_saved_pilot_record_replays_without_a_physical_gate():
    path = ROOT / 'results/development/nsc-ks-current-field-pilot-v2.json'
    if not path.exists():
        pytest.skip('one-cell field record is not present')
    sys.path.insert(0, str(ROOT / 'scripts'))
    import validate_nsc_ks_current_field_pilot_v2 as V
    saved = V.check()
    assert saved['physical_local_gate'] == 'OPEN'
    assert saved['physical_EXISTENCE_certificate'] is False
    assert saved['cell_index'] == 122
    assert saved['coverage']['all_history_segments'] is False
    assert saved['named_obstruction'] == 'cheb32_mixed_profile_fourier_omitted_band_from_global_derivative_l1'
