"""Coordinate and full-mode checks for the raw KS metric pullback."""
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from recursive_horizons.nsc_pg_ks_metric_pullback import coordinate_residuals,reference_chart,pg_log_jacobian,ks_to_pg,fixed_surface_normal
from derive_nsc_pg_ks_metric_pullback import load,evaluate


def test_same_metric_and_independent_lapse_shift_directions():
    for rho in (-.8,0.,.8):
        r=coordinate_residuals(rho)
        assert r['metric_tensor_reconstruction']<3e-11
        assert r['inverse_field_map']<3e-11
        assert r['jacobian_vs_finite_difference']<3e-8
    _,_,a,_,_=reference_chart(0.)
    J=pg_log_jacobian(0.,1.,0.,a,1.)
    assert np.linalg.matrix_rank(J)==4 and abs(np.linalg.det(J)-1)<3e-11
    assert abs(J[0,1])>1 and abs(J[2,0])>.1
    with pytest.raises(ValueError,match='Cauchy chart'):ks_to_pg(0.,100.,0.,a,1.)


def test_raw_KS_mode_variation_is_not_constant_seed_conversion():
    record=json.loads((ROOT/'results/development/nsc-pg-ks-metric-pullback.json').read_text())
    old=json.loads((ROOT/'results/development/nsc-pg-ctp-mode-jets.json').read_text())
    rows,maxima,chart,normal=evaluate(load(ROOT/record['payload']['path']),load(ROOT/old['payload']['path']))
    assert len(rows)==63
    assert maxima['weak_vs_direct']<3e-8
    assert maxima['CTP_tangent']<3e-11 and maxima['CTP_trace']<3e-11
    assert min(r['constant_neck_matrix_error_control'] for r in rows.values())>1e-4
    assert record['jet_inventory']['KS_two_endpoint_shape_and_embedding'] is None


def test_fixed_normal_is_geometric_not_a_Cauchy_isometry():
    _,_,a,_,_=reference_chart(0.)
    result=fixed_surface_normal(0.,1.,0.,a,1.)
    assert result['normal_geometry_residual']<3e-11
    assert result['normal_unit_residual']<3e-11
    h=1e-5;K=np.array([1.,0.,a,1.])
    for B in range(4):
        e=np.zeros(4);e[B]=h
        finite=(fixed_surface_normal(0.,*(K+e))['normal']-fixed_surface_normal(0.,*(K-e))['normal'])/(2*h)
        assert np.max(abs(finite-result['raw_KS_derivative'][:,B]))<3e-8
