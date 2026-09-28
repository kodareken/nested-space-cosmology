#!/usr/bin/env python3
"""Restore explicitly selected historical search payloads from pinned Git.

The active calculations do not need these old operator/tangent caches.
Historical replay commands do. This restores original bytes without running
an evolution, changing a branch, or overwriting any existing file.
"""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
INDEX = ROOT / 'docs/storage/historical-search-payloads.json'


def load_index(path):
    record = json.loads(Path(path).read_text())
    if record.get('schema') != 'NSC-HISTORICAL-SEARCH-PAYLOAD-ARCHIVE-v1':
        raise ValueError('unexpected archive index schema')
    if not re.fullmatch('[0-9a-f]{40}',record.get('source_commit','')):
        raise ValueError('full pinned Git commit required')
    seen=set()
    for row in record['files']:
        name=row['path'];p=PurePosixPath(name)
        if (p.is_absolute() or '..' in p.parts or str(p)!=name
                or not name.startswith('results/development/artifacts/')
                or p.suffix!='.npz' or name in seen):
            raise ValueError('noncanonical or duplicate historical payload path')
        if (type(row['bytes']) is not int or row['bytes']<0
                or not re.fullmatch('[0-9a-f]{64}',row['sha256'])):
            raise ValueError('payload byte count and SHA-256 required')
        seen.add(name)
    if sum(row['bytes'] for row in record['files'])!=record['bytes']:
        raise ValueError('archive byte total differs')
    return record


def file_hash(path):
    h=sha256()
    with Path(path).open('rb') as stream:
        while block:=stream.read(8*1024*1024):h.update(block)
    return h.hexdigest()


def restore(root,commit,row):
    root=Path(root).resolve();target=root/row['path']
    if not target.resolve().is_relative_to(root):raise ValueError('escaping destination')
    if any(p.is_symlink() for p in (target,*target.parents) if p!=root.parent):
        raise ValueError('symlink destination is not permitted')
    if target.exists():
        if not target.is_file() or file_hash(target)!=row['sha256']:
            raise FileExistsError('existing destination has different bytes: '+row['path'])
        return 'already present and verified'
    target.parent.mkdir(parents=True,exist_ok=True)
    fd,temporary=tempfile.mkstemp(prefix='.restore-',dir=target.parent)
    process=None
    try:
        h=sha256();size=0
        with os.fdopen(fd,'wb') as out:
            process=subprocess.Popen(['git','-C',str(root),'show',commit+':'+row['path']],
                stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            while block:=process.stdout.read(8*1024*1024):
                out.write(block);h.update(block);size+=len(block)
            error=process.stderr.read();rc=process.wait()
            if rc:raise RuntimeError('pinned Git payload unavailable: '+error.decode(errors='replace')[:400])
            out.flush();os.fsync(out.fileno())
        if size!=row['bytes'] or h.hexdigest()!=row['sha256']:
            raise ValueError('restored payload does not match recorded bytes/hash')
        # Atomic publication refuses a concurrent replacement; own temp is removed.
        os.link(temporary,target)
        return 'restored and verified'
    finally:
        if process is not None and process.poll() is None:
            process.kill();process.wait()
        Path(temporary).unlink(missing_ok=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--status',action='store_true')
    modes.add_argument('--restore',metavar='REPOSITORY_RELATIVE_NPZ')
    modes.add_argument('--restore-directory',metavar='EXACT_REPOSITORY_RELATIVE_DIRECTORY')
    args=parser.parse_args();record=load_index(INDEX)
    if args.status:
        absent=[r for r in record['files'] if not (ROOT/r['path']).is_file()]
        print(json.dumps({'source_commit':record['source_commit'],'indexed_files':len(record['files']),
            'absent_files':len(absent),'absent_bytes':sum(r['bytes'] for r in absent),
            'scientific_evolution_required':False},indent=2));return
    selected=[r for r in record['files'] if r['path']==args.restore or
              str(PurePosixPath(r['path']).parent)==args.restore_directory]
    if not selected:raise ValueError('selection is not present in the historical archive index')
    for row in selected:
        print(row['path']+': '+restore(ROOT,record['source_commit'],row),flush=True)


if __name__=='__main__':main()
