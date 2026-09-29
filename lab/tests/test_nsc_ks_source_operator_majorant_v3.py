from pathlib import Path
from copy import deepcopy
import json,runpy
import pytest
from flint import ctx
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
ROOT=Path(__file__).resolve().parents[1]


def test_low_pair_is_distinct_bound_and_uses_each_signed_error():
    d=runpy.run_path(str(ROOT/'scripts/derive_nsc_ks_source_operator_majorant_v3.py'))
    base=d['read_bound']('results/development/nsc-ks-source-operator-majorant-v2.json')
    low=d['read_bound'](d['LOW15'])
    rows=[{**r,'angular_sign':r['energy_sign'],'owner':r['bound_owner'],
            'epsilon':r['covariance_error_upper']} for r in base['source_coverage']['rows']]
    joined=d['add_low_pair'](rows,low,base['rho_up_hex'])
    assert len(joined)==868
    assert joined[-2]['epsilon']==low['positive_covariance_error_upper']
    assert joined[-1]['epsilon']==low['negative_covariance_error_upper']
    a=RetainedUpstreamArchive(ROOT)
    with ctx.workprec(192):
        _,_,covered,missing=d['source_error_moments'](a,a.family_entries((14,1)),joined)
        assert len(covered)==868 and missing==300
    with pytest.raises(ValueError,match='duplicate'):d['add_low_pair'](joined,low,base['rho_up_hex'])
    for key,value in [('covariance_error_unweighted',False),('rho_up_hex','0x1.0p+0'),('source_row',0)]:
        bad=deepcopy(low);bad[key]=value
        with pytest.raises(ValueError):d['add_low_pair'](rows,bad,base['rho_up_hex'])
    bad=deepcopy(low);bad['negative_covariance_error_upper']=None
    with pytest.raises(ValueError,match='missing'):d['add_low_pair'](rows,bad,base['rho_up_hex'])
