from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import subprocess

import pytest

path=Path(__file__).resolve().parents[1]/'scripts/restore_nsc_search_payloads.py'
spec=importlib.util.spec_from_file_location('restore_search',path)
owner=importlib.util.module_from_spec(spec);spec.loader.exec_module(owner)


def fixture_repo(tmp_path):
    subprocess.run(['git','init','-q',str(tmp_path)],check=True)
    name='results/development/artifacts/example/one.npz'
    payload=b'original immutable payload\x00\xff'
    p=tmp_path/name;p.parent.mkdir(parents=True);p.write_bytes(payload)
    subprocess.run(['git','-C',str(tmp_path),'add',name],check=True)
    subprocess.run(['git','-C',str(tmp_path),'-c','user.name=Test','-c','user.email=test@example.invalid',
                    'commit','-qm','fixture'],check=True)
    commit=subprocess.check_output(['git','-C',str(tmp_path),'rev-parse','HEAD'],text=True).strip()
    p.unlink()
    return commit,{'path':name,'sha256':sha256(payload).hexdigest(),'bytes':len(payload)},payload


def test_restore_recovers_original_bytes_without_changing_head(tmp_path):
    commit,row,payload=fixture_repo(tmp_path)
    assert owner.restore(tmp_path,commit,row)=='restored and verified'
    assert (tmp_path/row['path']).read_bytes()==payload
    assert owner.restore(tmp_path,commit,row)=='already present and verified'
    assert subprocess.check_output(['git','-C',str(tmp_path),'rev-parse','HEAD'],text=True).strip()==commit
    (tmp_path/row['path']).write_bytes(b'user change')
    with pytest.raises(FileExistsError):owner.restore(tmp_path,commit,row)
    assert (tmp_path/row['path']).read_bytes()==b'user change'


def test_bad_digest_cannot_publish_a_payload(tmp_path):
    commit,row,_=fixture_repo(tmp_path)
    with pytest.raises(ValueError,match='hash'):owner.restore(tmp_path,commit,{**row,'sha256':'0'*64})
    assert not (tmp_path/row['path']).exists()
    assert not list((tmp_path/row['path']).parent.glob('.restore-*'))


def test_index_rejects_traversal_and_duplicates(tmp_path):
    commit,row,_=fixture_repo(tmp_path)
    record={'schema':'NSC-HISTORICAL-SEARCH-PAYLOAD-ARCHIVE-v1','source_commit':commit,
            'files':[row],'bytes':row['bytes']}
    path=tmp_path/'index.json';path.write_text(json.dumps(record));owner.load_index(path)
    record['files']=[{**row,'path':'results/development/artifacts/../../outside.npz'}]
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError):owner.load_index(path)
    record['files']=[row,row];record['bytes']=2*row['bytes'];path.write_text(json.dumps(record))
    with pytest.raises(ValueError):owner.load_index(path)
