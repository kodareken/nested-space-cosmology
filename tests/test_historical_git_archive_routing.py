"""Process-local original-object routing; no historical code or producer runs."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import runpy
from types import ModuleType

import pytest

ROOT=Path(__file__).resolve().parents[1]
CONTEXT='329ef809e466a86ed4a290ee32c19485f85add9d'


def object_archive(root):
    objects=root/'archive/source-history/objects'
    objects.mkdir(parents=True)
    return objects


def test_threshold_metadata_read_declares_local_authentication_archive(tmp_path,monkeypatch):
    from recursive_horizons import nsc_threshold18_metadata as owner
    objects=object_archive(tmp_path);calls=[]
    monkeypatch.setenv('GIT_ALTERNATE_OBJECT_DIRECTORIES','unrelated-object-store')
    def read(command,**options):
        calls.append((command,options));return b'original metadata bytes'
    monkeypatch.setattr(subprocess,'check_output',read)
    relative=next(iter(owner.METADATA_IDENTITIES))
    assert owner._original_metadata_bytes(tmp_path,relative)==b'original metadata bytes'
    command,options=calls.pop()
    assert command[-1]==owner.AUTHENTICATION_COMMIT+':'+relative
    assert options['env']['GIT_ALTERNATE_OBJECT_DIRECTORIES']==str(objects)
    with pytest.raises(ValueError,match='not authorized'):
        owner._original_metadata_bytes(tmp_path,'results/unlisted.json')
    assert not calls


def test_original_graph_inventory_uses_same_archive_as_source_authentication(tmp_path,monkeypatch):
    from recursive_horizons import nsc_frozen_sources as frozen
    from recursive_horizons.nsc_source_binding_graph import original_source_graph
    objects=object_archive(tmp_path);raw=b'"""An exact original implementation."""\n';calls=[]
    expected=hashlib.sha256(raw).hexdigest();name='scripts/owner.py'
    monkeypatch.setattr(frozen,'source_sha256',lambda *args,**kwargs:expected)
    monkeypatch.setattr(frozen,'authenticate_source',lambda *args,**kwargs:raw)
    def listing(command,**options):
        calls.append((command,options));return (name+'\n').encode()
    monkeypatch.setattr(subprocess,'check_output',listing)
    original_source_graph.cache_clear()
    try:
        assert original_source_graph(tmp_path,(name,),CONTEXT)=={name:expected}
    finally:
        original_source_graph.cache_clear()
    command,options=calls.pop()
    assert command[-1]==CONTEXT and 'ls-tree' in command
    assert options['env']['GIT_ALTERNATE_OBJECT_DIRECTORIES']==str(objects)


@pytest.mark.parametrize('filename,function',[
    ('replay_history.py','git'),
])
def test_standalone_historical_commands_route_only_their_local_archive(tmp_path,monkeypatch,filename,function):
    objects=object_archive(tmp_path)
    spec=importlib.util.spec_from_file_location('archive_routing_'+filename[:-3],ROOT/'scripts'/filename)
    owner=importlib.util.module_from_spec(spec);spec.loader.exec_module(owner)
    monkeypatch.setattr(owner,'ROOT',tmp_path);calls=[]
    def output(command,**options):
        calls.append((command,options));return b'original bytes'
    def run(command,**options):
        calls.append((command,options));return subprocess.CompletedProcess(command,0,b'original bytes',b'')
    monkeypatch.setattr(subprocess,'check_output',output);monkeypatch.setattr(subprocess,'run',run)
    if function=='git_read':actual=owner.git_read(CONTEXT,'scripts/owner.py')
    else:actual=owner.git('show',CONTEXT+':scripts/owner.py')
    assert actual==b'original bytes'
    command,options=calls.pop()
    assert command==['git','show',CONTEXT+':scripts/owner.py']
    assert options['cwd']==tmp_path
    assert options['env']['GIT_ALTERNATE_OBJECT_DIRECTORIES']==str(objects)


def test_finite_companion_entrypoint_routes_archive_without_editing_frozen_helpers(tmp_path,monkeypatch):
    import os
    objects=object_archive(tmp_path);calls=[]
    finite=ModuleType('finite_regeneration');finite.ROOT=tmp_path
    def verify():
        calls.append(os.environ['GIT_ALTERNATE_OBJECT_DIRECTORIES'])
        return {'scope':'authentication only'}
    def forbidden(*args,**kwargs):raise AssertionError('build or producer called')
    finite.verify=verify
    finite.build=finite.compile_twice=finite.source_files=forbidden
    monkeypatch.setitem(sys.modules,'finite_regeneration',finite)
    monkeypatch.setenv('GIT_ALTERNATE_OBJECT_DIRECTORIES','unrelated-object-store')
    monkeypatch.setattr(sys,'argv',['build_finite_regeneration.py','--check'])
    runpy.run_path(str(ROOT/'scripts/build_finite_regeneration.py'),run_name='__main__')
    assert calls==[str(objects)]


def source_archive_manifest(root, *, name='docs/original.md', raw=b'original source\n', contexts=None):
    blob=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    value={'schema':'NSC-PINNED-SOURCE-OBJECTS-v1','authentication_archive':{
        'schema':'NSC-SPARSE-AUTHENTICATION-ARCHIVE-v1',
        'contexts':contexts if contexts is not None else [
            {'commit':CONTEXT,'paths':[{'path':name,'blob':blob}]}]}}
    path=root/'archive/source-history/manifest.json';path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value),encoding='utf-8')
    return value,path,blob


def test_default_legacy_resolution_authenticates_original_archive_context():
    from recursive_horizons.nsc_ks_evaluation_binding import resolve_recorded_source_bytes,LEGACY_REPLAY_MODE
    expected='7b6c0c79a27971a28bfcc7dbd23790c45be02da4488ca3d57319f37072372fb4'
    raw,origin=resolve_recorded_source_bytes(ROOT,'docs/nsc-ks-gate-value.md',expected,mode=LEGACY_REPLAY_MODE)
    assert hashlib.sha256(raw).hexdigest()==expected
    assert origin==CONTEXT and len(raw)==739


@pytest.mark.parametrize('mode,existing',[('production',False),('production',True)])
def test_production_never_consults_source_archive(tmp_path,monkeypatch,mode,existing):
    from recursive_horizons import nsc_ks_evaluation_binding as owner
    name='docs/original.md'
    if existing:
        path=tmp_path/name;path.parent.mkdir();path.write_bytes(b'changed\n')
    def forbidden(*args,**kwargs):raise AssertionError('archive or historical read attempted')
    monkeypatch.setattr(owner,'_declared_archived_source_bytes',forbidden)
    monkeypatch.setattr(owner,'_git_log_commits',forbidden)
    with pytest.raises(ValueError,match='production source|source hash mismatch'):
        owner.resolve_recorded_source_bytes(tmp_path,name,'0'*64,mode=mode)


@pytest.mark.parametrize('name',['results/receipt.json','docs/numeric.json','paper/evidence/proof.md'])
def test_declared_archive_defensively_rejects_numeric_and_evidence_paths(tmp_path,monkeypatch,name):
    from recursive_horizons import nsc_ks_evaluation_binding as owner,nsc_frozen_sources as frozen
    source_archive_manifest(tmp_path,name=name)
    def forbidden(*args,**kwargs):raise AssertionError('historical output authentication attempted')
    monkeypatch.setattr(frozen,'authenticate_source',forbidden)
    with pytest.raises(ValueError,match='implementation source'):
        owner._declared_archived_source_bytes(tmp_path,name,'0'*64,max_contexts=32)
    with pytest.raises(ValueError,match='recorded data hash mismatch'):
        owner.resolve_recorded_source_bytes(tmp_path,name,'0'*64,mode=owner.LEGACY_REPLAY_MODE)


def test_explicit_context_never_falls_back_to_other_archive_context(tmp_path,monkeypatch):
    from recursive_horizons import nsc_ks_evaluation_binding as owner
    def forbidden(*args,**kwargs):raise AssertionError('explicit pin fell back')
    monkeypatch.setattr(owner,'_declared_archived_source_bytes',forbidden)
    monkeypatch.setattr(owner,'_git_show',lambda *args,**kwargs:b'wrong original source')
    with pytest.raises(ValueError,match='pinned git source hash mismatch'):
        owner.resolve_recorded_source_bytes(tmp_path,'docs/original.md','0'*64,
            mode=owner.LEGACY_REPLAY_MODE,git_commit=CONTEXT)


@pytest.mark.parametrize('callback',['git_log','git_show'])
def test_custom_lookup_contract_has_no_implicit_archive(tmp_path,monkeypatch,callback):
    from recursive_horizons import nsc_ks_evaluation_binding as owner
    def forbidden(*args,**kwargs):raise AssertionError('implicit archive lookup bypassed callback contract')
    def empty_history(root,name,limit,git_log=None):
        return list(git_log(name,limit)) if git_log is not None else []
    monkeypatch.setattr(owner,'_declared_archived_source_bytes',forbidden)
    monkeypatch.setattr(owner,'_git_log_commits',empty_history)
    supplied={callback:(lambda *args:[] if callback=='git_log' else b'unused')}
    with pytest.raises(ValueError,match='historical source bytes not found'):
        owner.resolve_recorded_source_bytes(tmp_path,'docs/original.md','0'*64,
            mode=owner.LEGACY_REPLAY_MODE,**supplied)


@pytest.mark.parametrize('limit',[-1,True,1.5,'1',None])
def test_archive_context_count_is_a_strict_nonnegative_integer(tmp_path,limit):
    from recursive_horizons.nsc_ks_evaluation_binding import _declared_archived_source_bytes
    with pytest.raises(ValueError,match='nonnegative integer'):
        _declared_archived_source_bytes(tmp_path,'docs/original.md','0'*64,max_contexts=limit)


@pytest.mark.parametrize('present',[False,True])
def test_optional_archive_metadata_absence_has_no_candidates(tmp_path,present):
    from recursive_horizons.nsc_ks_evaluation_binding import _declared_archived_source_bytes
    if present:
        _,path,_=source_archive_manifest(tmp_path)
        path.write_text(json.dumps({'schema':'NSC-PINNED-SOURCE-OBJECTS-v1'}),encoding='utf-8')
    assert _declared_archived_source_bytes(tmp_path,'docs/original.md','0'*64,max_contexts=32) is None


@pytest.mark.parametrize('mutation',['manifest_schema','archive_schema','contexts','paths','commit','blob',
    'duplicate_context','duplicate_path','null_archive','invalid_json'])
def test_present_corrupt_archive_metadata_fails_closed(tmp_path,monkeypatch,mutation):
    from recursive_horizons import nsc_ks_evaluation_binding as owner,nsc_frozen_sources as frozen
    value,path,_=source_archive_manifest(tmp_path)
    archive=value['authentication_archive'];context=archive['contexts'][0]
    if mutation=='manifest_schema':value['schema']='unrecognized'
    elif mutation=='archive_schema':archive['schema']='unrecognized'
    elif mutation=='contexts':archive['contexts']={}
    elif mutation=='paths':context['paths']={}
    elif mutation=='commit':context['commit']='329ef809'
    elif mutation=='blob':context['paths'][0]['blob']='0'*39
    elif mutation=='duplicate_context':archive['contexts'].append(context)
    elif mutation=='duplicate_path':context['paths'].append(context['paths'][0])
    elif mutation=='null_archive':value['authentication_archive']=None
    if mutation=='invalid_json':path.write_text('{',encoding='utf-8')
    else:path.write_text(json.dumps(value),encoding='utf-8')
    def forbidden(*args,**kwargs):raise AssertionError('corrupt metadata reached source read')
    monkeypatch.setattr(frozen,'authenticate_source',forbidden)
    with pytest.raises(ValueError):
        owner._declared_archived_source_bytes(tmp_path,'docs/original.md','0'*64,max_contexts=32)


def test_archive_candidate_limit_does_not_reach_later_matching_body(tmp_path,monkeypatch):
    from recursive_horizons import nsc_ks_evaluation_binding as owner,nsc_frozen_sources as frozen
    raw=b'original source\n';expected=hashlib.sha256(raw).hexdigest()
    value,path,blob=source_archive_manifest(tmp_path,raw=raw)
    first='1'*40;second='2'*40
    value['authentication_archive']['contexts']=[{'commit':c,'paths':[{'path':'docs/original.md','blob':blob}]} for c in (first,second)]
    path.write_text(json.dumps(value),encoding='utf-8');calls=[]
    def authenticate(root,name,digest,*,context):
        calls.append(context)
        if context==first:raise RuntimeError('original source digest differs')
        return raw
    monkeypatch.setattr(frozen,'authenticate_source',authenticate)
    assert owner._declared_archived_source_bytes(tmp_path,'docs/original.md',expected,max_contexts=0) is None
    assert calls==[]
    assert owner._declared_archived_source_bytes(tmp_path,'docs/original.md',expected,max_contexts=1) is None
    assert calls==[first]
    calls.clear()
    assert owner._declared_archived_source_bytes(tmp_path,'docs/original.md',expected,max_contexts=2)==(raw,second)
    assert calls==[first,second]


def test_archive_requires_both_exact_recorded_sha_and_declared_blob(tmp_path,monkeypatch):
    from recursive_horizons import nsc_ks_evaluation_binding as owner,nsc_frozen_sources as frozen
    raw=b'original source\n';expected=hashlib.sha256(raw).hexdigest()
    value,path,_=source_archive_manifest(tmp_path,raw=raw)
    monkeypatch.setattr(frozen,'authenticate_source',lambda *args,**kwargs:raw)
    assert owner._declared_archived_source_bytes(tmp_path,'docs/original.md','0'*64,max_contexts=32) is None
    value['authentication_archive']['contexts'][0]['paths'][0]['blob']='0'*40
    path.write_text(json.dumps(value),encoding='utf-8')
    with pytest.raises(ValueError,match='blob binding differs'):
        owner._declared_archived_source_bytes(tmp_path,'docs/original.md',expected,max_contexts=32)
