"""Focused physical-preparation checks using the locked LLL inputs."""
import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_pg_lll_preparation import PGLLLPreparation, require_full_retained_preparation

ROOT=Path(__file__).resolve().parents[1]


def owner():
    r=json.loads((ROOT/'results/development/charged-ctp-neck-source.json').read_text())
    c=r['runs']['base']['config']
    return PGLLLPreparation(c['horizon_rho'],c['surface_gravity'],c['omega'])


def test_full_horizon_state_contains_partner_and_restricts_to_seed():
    p=owner();s=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())
    a=np.load(ROOT/s['payload']['path'],allow_pickle=False);c=s['channels'][0]
    for index in range(c['sample_count']):
        e=a['frequency'][index];full=p.spectral_covariance(e)
        assert np.allclose(p.seed_restriction(e),a['covariance_seed'][index],rtol=0,atol=3e-11)
        assert np.linalg.eigvalsh(full).min()>-3e-11
        assert np.linalg.eigvalsh(full).max()<1+3e-11
    assert abs(p.spectral_covariance(a['frequency'][0])[0,1])>0.1


def test_spatial_preparation_has_car_and_bulk_correlations():
    p=owner();r=p.project(64);c=r['covariance'];g=r['CAR_gram']
    assert np.linalg.norm(g-np.eye(7))<3e-11
    assert np.linalg.eigvalsh(c).min()>0 and np.linalg.eigvalsh(c).max()<1
    assert np.linalg.norm(c[:4,4:])>0.1
    assert np.linalg.norm(c[:4,6])>0.01
    assert np.linalg.norm(r['greater_equal_time']-r['lesser_equal_time']+1j*g)<3e-11
    with pytest.raises(ValueError,match='32 massive'):require_full_retained_preparation(1)
    with pytest.raises(ValueError,match='distribution'):p.off_diagonal_kernel(.5,.5)
