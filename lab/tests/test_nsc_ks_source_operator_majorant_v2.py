"""Coverage replacement controls on actual saved original-source witnesses."""
from copy import deepcopy
from pathlib import Path
import runpy
import pytest
from flint import ctx
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_middle_source_coverage import middle_catalogue,load_checkpoints,aggregate
ROOT=Path(__file__).resolve().parents[1]

@pytest.fixture(scope='module')
def data():
    d=runpy.run_path(str(ROOT/'scripts/derive_nsc_ks_source_operator_majorant_v2.py'))
    archive=RetainedUpstreamArchive(ROOT)
    records=[r for _,r,_ in load_checkpoints(ROOT/d['WINDOW'])]
    window=aggregate(middle_catalogue(archive),records,mode='check',row_budget=None,
      cpu_budget=None,cpu_seconds=0,capture_cpu_seconds=0,new_rows_captured=0)
    parents=[d['read_bound'](d[k]) for k in ('HIGH','LOW','MIDDLE')]
    return d,archive,parents,window


def test_middle_window_replaces_old_pair_and_binds_actual_source(data):
    d,a,parents,window=data
    with ctx.workprec(192):
        rows=d['merge_bounds'](*parents,window)
        assert len(rows)==866
        assert sum(r['panel']=='group14/mid24_1' and r['row']==0 for r in rows)==2
        entries=a.family_entries((14,1))
        S0,S1,covered,missing=d['source_error_moments'](a,entries,rows)
        assert S1>S0>0 and len(covered)==866 and missing==302
        with pytest.raises(ValueError,match='duplicate'):
            d['source_error_moments'](a,entries,rows+[rows[-1]])
        altered=deepcopy(rows);altered[-1]['preparation_digest']='0'*64
        with pytest.raises(ValueError,match='original preparation'):
            d['source_error_moments'](a,entries,altered)


def test_missing_duplicate_and_wrong_signed_window_cannot_aggregate(data):
    d,_,_,window=data
    missing=deepcopy(window);missing['rows'].pop()
    with pytest.raises(ValueError,match='census'):d['window_rows'](missing)
    duplicate=deepcopy(window);duplicate['rows'][-1]=duplicate['rows'][0]
    with pytest.raises(ValueError,match='duplicate'):d['window_rows'](duplicate)
    wrong=deepcopy(window);wrong['rows'][0]['negative_energy_hex']=wrong['rows'][0]['energy_hex']
    with pytest.raises(ValueError,match='signed middle energy'):d['window_rows'](wrong)
    absent=deepcopy(window);absent['rows'][0]['source_comparisons'].pop()
    with pytest.raises(ValueError,match='signed source bounds'):d['window_rows'](absent)
    partial=deepcopy(window);partial['coverage_complete']=False
    with pytest.raises(ValueError,match='complete'):d['window_rows'](partial)
