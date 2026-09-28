"""Triangular characteristic majorant and source-carrier conversion."""
from fractions import Fraction as Q

import pytest
pytest.importorskip('flint')
from flint import arb, ctx, fmpq
from recursive_horizons.nsc_ks_characteristic_error import propagate_local_characteristics
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper


def test_majorant_contains_exact_constant_comparison_system():
    before = ctx.prec
    report = propagate_local_characteristics((Q(1, 10), Q(1, 20)), (0, 0), Q(1, 5), Q(3, 10), 7)
    assert ctx.prec == before
    with ctx.workprec(150):
        growth = arb(fmpq(1, 5)).exp()
        e0 = growth/10
        e1 = growth*(arb(fmpq(1, 20))+arb(fmpq(3, 100)))
        assert restored_upper(report['row_field_error_upper']) >= e0.upper()
        assert restored_upper(report['row_axial_envelope_error_upper']) >= e1.upper()
        assert restored_upper(report['weighted_F_z_error_upper']) >= (arb(2).sqrt()*(e1+7*e0)).upper()
    assert report['requires_continuous_cone_bounds']


def test_missing_residual_is_never_a_zero_error_certificate():
    with pytest.raises(ValueError, match='explicit'):
        propagate_local_characteristics((0, 0), (None, 0), 0, 0, 1)
    result = propagate_local_characteristics((0, 0), (0, 0), 0, 0, 1)
    assert restored_upper(result['weighted_F_error_upper']) == 0
    assert result['physical_local_gate'].startswith('OPEN')
