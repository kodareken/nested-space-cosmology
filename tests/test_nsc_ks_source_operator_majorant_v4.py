"""Direct-window successor: replayed pairs raise family14_1 from 868 to 1042."""
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import runpy
import pytest
from flint import ctx
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
ROOT=Path(__file__).resolve().parents[1]
HISTORIC=(
    ROOT/'results/development/nsc-ks-source-operator-majorant-v3.json',
    ROOT/'results/development/nsc-subgap-row15-upstream-v1.json',
    ROOT/'results/development/artifacts/nsc-subgap-row15-upstream-v1.npz',
    ROOT/'results/development/nsc-direct-source-window-v1',
)


def _snapshot(path):
    path=Path(path);digest=sha256()
    if path.is_file():
        digest.update(path.read_bytes());return digest.hexdigest()
    for item in sorted(path.rglob('*'),key=lambda candidate:candidate.relative_to(path).as_posix()):
        digest.update(item.relative_to(path).as_posix().encode());digest.update(b'\0')
        if item.is_symlink():digest.update(b'link')
        elif item.is_file():digest.update(item.read_bytes())
        else:digest.update(b'dir')
        digest.update(b'\0')
    return digest.hexdigest()


@pytest.fixture(scope='module',autouse=True)
def historic_bytes_stay_unchanged():
    before={str(path):_snapshot(path) for path in HISTORIC}
    yield
    assert {str(path):_snapshot(path) for path in HISTORIC}==before


@pytest.fixture(scope='module')
def data():
    d=runpy.run_path(str(ROOT/'scripts/derive_nsc_ks_source_operator_majorant_v4.py'))
    archive=RetainedUpstreamArchive(ROOT)
    prior=d['read_bound'](d['PRIOR'])
    window=d['cover'](ROOT,ROOT/d['WINDOW'],mode='check')
    return d,archive,prior,window


def test_replayed_pairs_raise_coverage_once_and_stay_disjoint(data):
    d,archive,prior,window=data
    assert window['replay_uses_saved_trace'] is True and window['new_rows_captured']==0
    base,added=d['prior_rows'](prior),d['direct_rows'](window)
    assert len(base)==868 and len(added)==174 and len(window['rows'])==87
    base_keys={(row['panel'],row['row'],row['energy_sign']) for row in base}
    added_keys={(row['panel'],row['row'],row['energy_sign']) for row in added}
    assert base_keys.isdisjoint(added_keys)
    assert {('group14/low16_1',0,1),('group14/low16_1',15,1)}<=base_keys
    assert {('group14/low16_1',16,1),('group14/low32_1',41,1)}<=added_keys
    pairs=list(zip(window['rows'],added[0::2],added[1::2]))
    assert all(left['epsilon']==item['signed_comparisons'][0]['signed_total_upper']
      and right['epsilon']==item['signed_comparisons'][1]['signed_total_upper']
      and left['source_digest']!=right['source_digest']
      and left['preparation_digest']!=right['preparation_digest']
      for item,left,right in pairs)
    assert any(left['epsilon']!=right['epsilon'] for _,left,right in pairs)
    rho,entries=d['require_same_slice'](prior,archive)
    assert rho.hex()==prior['rho_up_hex']
    joined=d['add_direct_pairs'](base,added)
    assert len(joined)==1042
    with ctx.workprec(192):
        moment,energy,covered,missing=d['source_error_moments'](archive,entries,joined)
        assert energy>moment>0 and len(covered)==1042 and missing==126
        by_key={(row['panel'],row['row'],row['energy_sign']):row for row in covered}
        for row in added:
            got=by_key[(row['panel'],row['row'],row['energy_sign'])]
            assert got['covariance_error_upper']==row['epsilon']
            assert got['source_digest']==row['source_digest']
            assert got['preparation_digest']==row['preparation_digest']
            assert got['measure_upper']!=row['epsilon']
    inputs=d['bind_direct_inputs']({},window)
    assert len(inputs)==174 and {path.rsplit('.',1)[-1] for path in inputs}=={'json','npz'}
    owners=implementation_hashes(ROOT,owners=d['OWNERS'])
    assert 'scripts/derive_nsc_ks_source_operator_majorant_v4.py' in owners
    assert 'tests/test_nsc_ks_source_operator_majorant_v4.py' in owners
    with pytest.raises(ValueError,match='duplicate'):d['add_direct_pairs'](joined,added[:1])
    with pytest.raises(ValueError,match='census'):d['add_direct_pairs'](base,added[:-1])
    altered=deepcopy(joined);altered[-1]['preparation_digest']='0'*64
    with pytest.raises(ValueError,match='original preparation'):
        d['source_error_moments'](archive,entries,altered)
    altered=deepcopy(joined);altered[-1]['epsilon']={'exponent':0,'mantissa':'-1'}
    with pytest.raises(ValueError,match='nonnegative'):
        d['source_error_moments'](archive,entries,altered)
    with pytest.raises(ValueError,match='duplicate'):
        d['source_error_moments'](archive,entries,joined+[joined[-1]])


def test_mutated_window_identity_weight_or_slice_cannot_join(data):
    d,archive,prior,window=data
    incomplete=deepcopy(window);incomplete['coverage_complete']=False
    with pytest.raises(ValueError,match='complete'):d['direct_rows'](incomplete)
    weighted=deepcopy(window);weighted['weight_applied_to_error']=True
    with pytest.raises(ValueError,match='weight'):d['direct_rows'](weighted)
    row_weighted=deepcopy(window);row_weighted['rows'][0]['weight_applied_to_error']=True
    with pytest.raises(ValueError,match='weight'):d['direct_rows'](row_weighted)
    resumed=deepcopy(window);resumed['mode']='resume';resumed['new_rows_captured']=1
    with pytest.raises(ValueError,match='replayed'):d['direct_rows'](resumed)
    duplicate=deepcopy(window);duplicate['rows'][1]=deepcopy(duplicate['rows'][0])
    with pytest.raises(ValueError,match='duplicate'):d['direct_rows'](duplicate)
    copied=deepcopy(window);pair=copied['rows'][0]['signed_comparisons']
    pair[1]['source_digest']=pair[0]['source_digest']
    with pytest.raises(ValueError,match='duplicate'):d['direct_rows'](copied)
    missing=deepcopy(window);missing['rows'][0]['signed_comparisons'][1]['signed_total_upper']=None
    with pytest.raises(ValueError,match='missing'):d['direct_rows'](missing)
    unsigned=deepcopy(window);unsigned['rows'][0]['signed_comparisons'].pop()
    with pytest.raises(ValueError,match='signed'):d['direct_rows'](unsigned)
    wrong=deepcopy(window);wrong['rows'][0]['negative_energy_fiber_hex']=wrong['rows'][0]['energy_fiber_hex']
    with pytest.raises(ValueError,match='signed direct energy'):d['direct_rows'](wrong)
    shifted=deepcopy(window);shifted['rows']=shifted['rows'][1:]+shifted['rows'][:1]
    with pytest.raises(ValueError,match='endpoints'):d['direct_rows'](shifted)
    shrunk=deepcopy(window);shrunk['rows']=shrunk['rows'][:-1]
    with pytest.raises(ValueError,match='unbound'):d['bind_direct_inputs']({},shrunk)
    path=next(iter(d['bind_direct_inputs']({},window)))
    with pytest.raises(ValueError,match='overlaps'):d['bind_direct_inputs']({path:'0'*64},window)
    changed=deepcopy(prior);changed['source_coverage']['covered_signed_rows']=867
    with pytest.raises(ValueError,match='authenticated'):d['prior_rows'](changed)
    moved=deepcopy(prior);moved['rho_up_hex']=(float.fromhex(prior['rho_up_hex'])*2).hex()
    with pytest.raises(ValueError,match='upstream archive slice changed'):
        d['require_same_slice'](moved,archive)
    source=(ROOT/'scripts/derive_nsc_ks_source_operator_majorant_v4.py').read_text()
    assert "mode='check'" in source and 'capture_direct_vacuum' not in source
    assert "mode='resume'" not in source and 'replay_v3()' in source
    assert d['OUTPUT'].endswith('nsc-ks-source-operator-majorant-v4.json')
    assert d['PRIOR'].endswith('nsc-ks-source-operator-majorant-v3.json')
