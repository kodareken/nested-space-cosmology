"""Cold version selection and full compact historical regression."""
import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_compact_matched_restart import CompactSeedGeneration,HISTORICAL,MATCHED,require_convention

ROOT=Path(__file__).resolve().parents[1]


def generations():
    record=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())
    p=record['seed_generations_payload']
    hist=CompactSeedGeneration.load(ROOT/p['path'],start_convention=HISTORICAL,expected_sha256=p['sha256'])
    matched=CompactSeedGeneration.load(ROOT/p['path'],start_convention=MATCHED,expected_sha256=p['sha256'])
    return hist,matched


def test_explicit_versions_preserve_original_control_and_new_difference():
    historical,matched=generations()
    state=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())
    with np.load(ROOT/state['payload']['path'],allow_pickle=False) as old:
        mask=historical.labels[:,1]==1;index=historical.labels[mask,2].astype(int)
        assert int(mask.sum())==1280
        assert np.max(np.linalg.norm(historical.covariance[mask]-old['covariance_seed'][index],axis=(1,2)))<3e-11
    assert historical.version!=matched.version
    assert np.max(np.linalg.norm(historical.covariance-matched.covariance,axis=(1,2)))>1e-5
    assert np.array_equal(historical.labels,matched.labels)


def test_new_matched_full_grid_has_car_and_signed_partners():
    _,matched=generations()
    assert len(matched.labels)==2432
    eig=np.linalg.eigvalsh(np.concatenate((matched.covariance,matched.negative_covariance)))
    assert eig.min()>-3e-11 and eig.max()<1+3e-11
    with pytest.raises(ValueError):matched.as_spatial_PG_covariance()


@pytest.mark.parametrize('name',[None,'matched','historical',''])
def test_no_implicit_start_convention(name):
    with pytest.raises(ValueError):require_convention(name)
