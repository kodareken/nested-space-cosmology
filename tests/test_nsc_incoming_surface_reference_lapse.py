"""Exact full-Weyl B_ref kernel, background terms and inherited normalization."""
import numpy as np
import sympy as sp

from recursive_horizons.nsc_incoming_surface_reference_lapse import (
    reference_lapse_principal,reference_lapse_group,S1,S2,S3,
)
from recursive_horizons.nsc_spatial_reference_symbol import SIGMA


def test_full_spatial_Weyl_identities_and_beta_moments_are_exact():
    np.testing.assert_array_equal(np.array([S1.tolist(),S2.tolist(),S3.tolist()],complex),SIGMA)
    result=reference_lapse_principal()
    assert all(set(values)=={'0'} for values in result['residuals'].values())
    assert set(result['trace_G4_pieces'])=={'ordinary_P1_deltaP3','Moyal2_P0_P2_pair','Moyal2_P1_P1'}
    assert all(value!=0 for value in result['trace_G4_pieces'].values())


def test_first_normal_background_and_full_line_factor_are_preserved():
    result=reference_lapse_principal();m,L,p,a,r,Ha,Hr,D=result['symbols']
    assert result['closed']==reference_lapse_group(m,L,a,r,Ha,Hr,D,sp.pi)
    assert result['closed'].diff(Ha)!=0 and result['closed'].diff(Hr)!=0
    assert sp.simplify(result['closed'].subs(L,0))==0
    assert sp.simplify(result['kernel'].subs({Ha:0,Hr:0}))==0
    assert result['third_order_odd_kernel']!=0
    assert sp.simplify(result['third_order_odd_kernel']+result['third_order_odd_kernel'].subs(p,-p))==0
