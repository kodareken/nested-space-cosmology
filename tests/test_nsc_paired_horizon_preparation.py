"""Finite horizon map, its corrected normalization, and honest scope guards."""
import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap

ROOT=Path(__file__).resolve().parents[1]


def owner():
    c=json.loads((ROOT/'results/development/charged-compact-ctp-completion.json').read_text())['runs']['base']['config']
    return PairedHorizonSeedMap(c['horizon_rho'],c['surface_gravity'],c['omega'],c['horizon_offset'],c['scattering_tolerance'],c['outer_floor'])


def test_q_collar_has_no_extra_inverse_horizon_radius():
    h=owner();m=np.pi/2;l=np.sqrt(5)
    assert h.frame_residual(.2,m,l)['differential_norm']<3e-11
    assert h.frame_residual(.2,m,l,legacy_compact_distance=True)['differential_norm']>1e-6


def test_massive_mode_follows_owned_current_operator():
    h=owner();p=h.at_radius(.5,.2,np.pi/2,np.sqrt(5))
    assert max(p['residuals'].values())<3e-11
    assert np.linalg.norm(p['seed_mode_map'])>0
    with pytest.raises(ValueError):h.require_global_PG_covariance()
    with pytest.raises(TypeError):h.at_seed(.2,1.,1.,seed_covariance=np.eye(2))
