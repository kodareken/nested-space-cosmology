"""Actual shifted-momentum coverage and retained-order controls."""
import numpy as np
import pytest
pytest.importorskip('flint')
from flint import arb,ctx
from recursive_horizons.nsc_scaled_reference_projector import scaled_reference_projector
from recursive_horizons.nsc_weyl_small_transfer_bound import small_transfer_coefficients,evaluate_coefficient_row
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper


def test_constant_geometry_zero_defect_and_no_physical_tail_claim():
    with ctx.workprec(160):
        jets=scaled_reference_projector(mu=arb((1,-7),(1,-7)),sign=1,mass=1,angular=2,axial=1,radius=2)
        report=small_transfer_coefficients(jets,arb(1),[0]*10,160)
        for row in report['rows']:
            assert restored_upper(row['C5_T_upper'])==0
            assert restored_upper(row['C6_T_upper'])==0
        assert report['large_transfer_bound'] is None
        assert not report['input_energy_tail_identified_with_canonical_k']
        assert report['physical_local_gate']=='OPEN'


def test_small_transfer_must_cover_shifted_segments_and_all_needed_moments():
    with ctx.workprec(160):
        point=scaled_reference_projector(mu=arb(1)/160,sign=-1,mass=1,angular=2,axial=1,radius=2)
        with pytest.raises(ValueError,match='shifted momentum'):
            small_transfer_coefficients(point,1,[0]*10,160)
        box=scaled_reference_projector(mu=arb((1,-7),(1,-7)),sign=1,mass=1,angular=2,axial=1,radius=2)
        with pytest.raises(ValueError,match='Fourier moments'):
            small_transfer_coefficients(box,1,[0]*5,160)


def test_tail_evaluation_uses_lower_momentum_not_rounded_upper():
    row={'C5_rho_upper':{'mantissa':'32','exponent':0},'C6_rho_upper':{'mantissa':'64','exponent':0}}
    with ctx.workprec(160):
        value=restored_upper(evaluate_coefficient_row(row,arb(2,arb(1)/4)))
        assert value>=arb(32)/(arb(7)/4)**5+arb(64)/(arb(7)/4)**6
