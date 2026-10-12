from copy import deepcopy
from pathlib import Path
import json
import pytest
from flint import arb,ctx
from recursive_horizons.nsc_ks_cubic_inventory import load_inventory,verify_inventory,channel_moments,combine_basis
from recursive_horizons.nsc_ks_massive_cubic_uv import CAUCHY_RECORD,CUTOFF_BRIDGE_RECORD,INVENTORY_RECORD
ROOT=Path(__file__).resolve().parents[1]


def test_original_census_and_group14_normalization():
    rows=load_inventory(ROOT)
    assert len(rows)==60
    assert {r['group'] for r in rows if r['cutoff']==320}=={10,11,12,31,32}
    with ctx.workprec(160):
        m=channel_moments([r for r in rows if r['group']==14])
        for v,w in zip(m['raw'],(60/arb.pi(),300/arb.pi(),15*arb.pi())):assert (v-w).contains(0)
        for raw,lead in zip(m['raw'],m['leading']):assert (lead-raw/(2*160**2)).contains(0)
        all_m=channel_moments(rows)
        assert all_m['leading'][0]<all_m['raw'][0]/(2*160**2)
        with pytest.raises(ValueError,match='duplicate'):channel_moments(rows+[rows[0]])


def test_census_multiplicity_and_cutoff_mutations_fail():
    original=[json.loads((ROOT/p).read_text()) for p in (CAUCHY_RECORD,CUTOFF_BRIDGE_RECORD,INVENTORY_RECORD)]
    bad=deepcopy(original);bad[1]['ledger'].pop()
    with pytest.raises(ValueError,match='missing'):verify_inventory(*bad)
    bad=deepcopy(original);bad[1]['ledger'][0]['multiplicity_per_signed_family']*=2
    with pytest.raises(ValueError,match='multiplicity'):verify_inventory(*bad)
    bad=deepcopy(original);bad[0]['channels'][13]['angular_eigenvalue']=1.0
    with pytest.raises(ValueError,match='angular label'):verify_inventory(*bad)
    bad=deepcopy(original);bad[2]['quadrature_measure_by_family']['10_1']=160
    with pytest.raises(ValueError,match='cutoff'):verify_inventory(*bad)


def test_both_signs_and_matching_domains_are_required():
    b={'sign':1,'massless_A':True,'profile_identity':'test','rho_up':arb(1.03),'z_domain':arb(2),'C_M':None}
    b.update({k:arb(i+1) for i,k in enumerate(('N2','N4','NM','J2','J4','JM'))})
    minus={**b,'sign':-1,'N2':arb(7)};mom={'raw':(arb(1),arb(2),arb(3)),'leading':(arb(1)/2,arb(1),arb(3)/2)}
    r=combine_basis({1:b,-1:minus},mom)
    assert r['raw']['N']==-34 and r['raw']['beta']==64
    assert r['leading']['N']==-17 and r['physical_residual_modified'] is False
    with pytest.raises(ValueError,match='both'):combine_basis({1:b},mom)
    with pytest.raises(ValueError,match='domain'):combine_basis({1:b,-1:{**minus,'z_domain':arb(3)}},mom)
