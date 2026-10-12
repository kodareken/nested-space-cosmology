"""Exact quadratic full-Weyl response and immutable signed n24 control."""
from pathlib import Path

import sympy as sp

from recursive_horizons.nsc_incoming_surface_quadratic_lapse import (
    reference_quadratic_lapse,reference_quadratic_group,archived_v_control,
)

ROOT=Path(__file__).resolve().parents[1]


def test_full_first_second_transport_and_fourth_trace_are_exact():
    result=reference_quadratic_lapse()
    assert all(set(v)=={'0'} for v in result['residuals'].values())
    assert result['P1_scalar_spatial_coefficient']!=0
    assert result['P2_point_scalar_coefficient']!=0
    assert all(v!=0 for v in result['trace_G4_pieces'].values())


def test_closed_coefficient_has_correct_factor_and_zero_angular_limit():
    result=reference_quadratic_lapse();m,L,p,a,r,D=result['symbols']
    assert result['closed']==reference_quadratic_group(m,L,a,r,D,sp.pi)
    assert sp.simplify(result['closed'].subs(L,0))==0
    assert sp.simplify(result['kernel']-result['kernel'].subs(p,-p))==0
    assert sp.simplify(result['kernel']-result['kernel'].subs(L,-L))==0


def test_old_n24_control_is_paired_before_comparison_without_rerun():
    result=archived_v_control(ROOT)
    assert result['paired_node_absolute_error']<result['control_tolerance']
    assert result['aggregate_action_error']<result['control_tolerance']
    assert result['angular_and_momentum_signs_combined_before_comparison']
    assert result['new_reference_evaluations']==result['new_momentum_quadratures']==0
