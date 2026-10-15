"""Process-local original-object routing; no historical code or producer runs."""
import hashlib
import importlib.util
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
