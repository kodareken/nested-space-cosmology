#!/usr/bin/env python3
"""Rebind eighteen live metadata carriers without executing scientific producers.

The corrected prepared-history owner-document binding changes its JSON digest.
Only keyed dependent file/payload identities are propagated. Original producer
source pins, execution signatures, scientific JSON spellings and every other
NPZ member stay fixed. Historical campaigns and gate captures are out of scope.
"""
from __future__ import annotations

import argparse
from copy import copy
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import tempfile
import zipfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCE_CONTEXT='329ef809e466a86ed4a290ee32c19485f85add9d'
BASELINE='88b2e7b1bf265772eac7108827ccb5539889e5bd'
SEED='results/development/nsc-compatible-prepared-history.json'
SEED_OLD='66cb7f086bb60bc646533ed09ca565cc02fedcafbc609d68b5dcf662a2e61318'
SEED_CURRENT='6eae36dbb91885caee51ef6d7e1ec0fa06983a70765f814761128b64c2d6514b'
ORIGINALS={
 'results/development/nsc-retarded-compatible-response.json':'af554fc9c18bf18780b96e258fe8cfb0dcc1df1dfb7da0480c1ef961bb4d61b9',
 'results/development/nsc-retarded-response-resolution.json':'2ee55dd8682dad78e7d578ffc91e949fbf88c3d1a0558a165b289c54f48a934f',
 'results/development/nsc-retarded-fourier-response.json':'3c29126f053de04670c5bd7e9addbc7943dfa2cd081034be45ea6d935cf1f4ea',
 'results/development/nsc-retarded-radial-response.json':'f90294b377cbbdd3e04ce93c21fa22d6f60b7443248e38727d4c49fe3045d24b',
 'results/development/nsc-upstream-compensation-rank.json':'52029e38112e66f0fda0f0df4778653aad1e153fa3692b0ce4223804fe9e822d',
 'results/development/nsc-upstream-compensator.json':'d7e07b91446abfcb8c9dfe7783687600028260e73794ca8aea803dae2af8e2aa',
 'results/development/nsc-retarded-spatial-diagnosis.json':'20534be1d0e56c7a42e1937fb742db03784ad15f0d44d15240ccbdab4efd3abf',
 'results/development/nsc-evolved-incoming-state.json':'1780cb2bfa0dea3a92d235e7d84b6faf8ed99ccc8a2bfaddf3316af123def4b7',
 'results/development/nsc-evolved-incoming-constraints.json':'a738e839459e907017ff6844ef2faca34fd21e28b28fe1c5135703d76a1ebee5',
 'results/development/nsc-incoming-endpoint-obstruction.json':'8d7c5dbe87efb65b0902900e4a6d9d9a1241457f98826b7db674314da867fed7',
 'results/development/nsc-incoming-fixed-transfer.json':'4cecaf278b5266ccc5fb78a2cbf5a143846427adcb00af5572c991418815bae2',
 'results/development/nsc-retarded-fixed-c0-certificate.json':'868a833d413c12b0a48f6be36e68f279a66d1946efbd620bcde9f9623fe073e6',
 'results/development/nsc-incoming-jet-rigidity.json':'c8047c0398d0bcd235a6218e7eec97339be137a294877123f2db2e1138d57ca0',
 'results/development/artifacts/nsc-retarded-compatible-response.npz':'72aafd825be4722ff6282867a9eff06f9ef82530a091d56300c45d6da9f599af',
 'results/development/artifacts/nsc-retarded-response-resolution.npz':'8d3a79e6c1dc5ada4bbe17e9c3450ad6576a25212b04bbe0cd6dcf76c3b997cf',
 'results/development/artifacts/nsc-retarded-radial-response.npz':'ed3dea217b2199527f6061e943b2af2a7b1981b1e2ef4e0eb8d35b1d33ccd135',
 'results/development/artifacts/nsc-upstream-compensator.npz':'957997e9b2462bb43dcf8c7ae491a349786e0448e63137743dcc8fb7a606f790',
 'results/development/artifacts/nsc-evolved-incoming-state.npz':'e7847d012402406b8eacf4e29a2c616f07b6957ec45b8efb3d17c533fea9e4da',
}
PROTECTED={'execution_signature','execution_source_hashes','executed_sources','unchanged_array_sha256'}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def json_tokens(raw):
    """Parse JSON while retaining every primitive token's exact spelling."""
    text=raw.decode('utf-8');decoder=json.JSONDecoder();tokens={}
    def white(i):
        while i<len(text) and text[i].isspace():i+=1
        return i
    def parse(i,pointer):
        i=white(i)
        if text[i]=='{':
            value={};i=white(i+1)
            if text[i]=='}':return value,i+1
            while True:
                key,i=decoder.raw_decode(text,i)
                if not isinstance(key,str) or key in value:raise ValueError('ambiguous JSON key')
                i=white(i)
                if text[i]!=':':raise ValueError('JSON colon required')
                value[key],i=parse(i+1,(*pointer,key));i=white(i)
                if text[i]=='}':return value,i+1
                if text[i]!=',':raise ValueError('JSON comma required')
                i=white(i+1)
        if text[i]=='[':
            value=[];i=white(i+1)
            if text[i]==']':return value,i+1
            while True:
                item,i=parse(i,(*pointer,len(value)));value.append(item);i=white(i)
                if text[i]==']':return value,i+1
                if text[i]!=',':raise ValueError('JSON comma required')
                i=white(i+1)
        value,end=decoder.raw_decode(text,i)
        tokens[pointer]=(value,i,end)
        return value,end
    value,end=parse(0,())
    if white(end)!=len(text):raise ValueError('trailing JSON content')
    return value,tokens,text


def at(value,pointer):
    for key in pointer:value=value[key]
    return value


def selected_bindings(value,tokens,identities):
    """Require the target path and its exact original digest together."""
    selected={};old_hashes={i['sha256'] for i in identities.values()}
    for pointer,(item,_start,_end) in tokens.items():
        if not pointer or any(k in PROTECTED for k in pointer):continue
        key=pointer[-1];parent=at(value,pointer[:-1]);target=None
        if isinstance(key,str) and key in identities:
            target=key
        elif key=='sha256' and isinstance(parent,dict) and parent.get('path') in identities:
            target=parent['path']
        if target is not None:
            if item!=identities[target]['sha256'] or not isinstance(item,str):
                raise ValueError('keyed original receipt digest changed: '+target)
            selected[pointer]=target
        elif isinstance(item,str) and item in old_hashes:
            raise ValueError('unqualified live digest cannot be rebound')
    return selected


def validate_json_changes(before,after,changes):
    _b,bt,bs=json_tokens(before);_a,ats,ass=json_tokens(after)
    if bt.keys()!=ats.keys():raise ValueError('receipt schema changed')
    permitted={tuple(c['pointer']) for c in changes}
    for pointer,(old,i,j) in bt.items():
        new,k,l=ats[pointer]
        if type(old) is not type(new):raise ValueError('receipt token type changed')
        if pointer not in permitted and bs[i:j]!=ass[k:l]:
            raise ValueError('nonmetadata receipt token changed')


def rewrite_json(raw,original,current):
    value,tokens,text=json_tokens(raw)
    bindings=selected_bindings(value,tokens,original);edits={};changes=[]
    for pointer,target in bindings.items():
        if target not in current:raise ValueError('parent must be rebound first: '+target)
        old,i,j=tokens[pointer];new=current[target]['sha256']
        edits[pointer]=(i,j,json.dumps(new))
        changes.append({'pointer':list(pointer),'target':target,'before':old,'after':new})
        parent=at(value,pointer[:-1]);size_pointer=(*pointer[:-1],'bytes')
        if pointer[-1]=='sha256' and 'bytes' in parent:
            size,a,b=tokens[size_pointer]
            if type(size) is not int or size!=original[target]['bytes']:
                raise ValueError('original payload size changed: '+target)
            new_size=current[target]['bytes']
            if new_size!=size:
                edits[size_pointer]=(a,b,str(new_size))
                changes.append({'pointer':list(size_pointer),'target':target,'before':size,'after':new_size})
    for i,j,new in sorted(edits.values(),reverse=True):text=text[:i]+new+text[j:]
    encoded=text.encode('utf-8');validate_json_changes(raw,encoded,changes)
    return encoded,changes


def metadata_json(raw_npz):
    with zipfile.ZipFile(io.BytesIO(raw_npz)) as archive:
        member=archive.read('metadata_json.npy')
    value=np.lib.format.read_array(io.BytesIO(member),allow_pickle=False)
    if value.dtype!=np.uint8 or value.ndim!=1:raise ValueError('uint8 JSON metadata carrier required')
    return value.tobytes()


def rewrite_npz(raw,original,current):
    old_json=metadata_json(raw);new_json,changes=rewrite_json(old_json,original,current)
    output=io.BytesIO();preserved={}
    with zipfile.ZipFile(io.BytesIO(raw)) as archive,zipfile.ZipFile(output,'w') as writer:
        names=archive.namelist()
        if len(names)!=len(set(names)):raise ValueError('duplicated NPZ member')
        for info in archive.infolist():
            member=archive.read(info.filename)
            if info.filename=='metadata_json.npy':
                if len(old_json)==len(new_json):
                    member=member[:-len(old_json)]+new_json
                else:
                    replacement=io.BytesIO();np.save(replacement,np.frombuffer(new_json,np.uint8),allow_pickle=False)
                    member=replacement.getvalue()
            else:preserved[info.filename]=sha(member)
            writer.writestr(copy(info),member)
    encoded=output.getvalue()
    with zipfile.ZipFile(io.BytesIO(raw)) as before,zipfile.ZipFile(io.BytesIO(encoded)) as after:
        if before.namelist()!=after.namelist():raise ValueError('NPZ schema changed')
        for name in preserved:
            if before.read(name)!=after.read(name):raise ValueError('nonmetadata NPZ member changed: '+name)
    if metadata_json(encoded)!=new_json:raise ValueError('metadata encoding changed')
    return encoded,changes,preserved


def regular_file(root,name):
    path=root/name
    if path.is_symlink() or not path.resolve().is_relative_to(root) or not path.is_file():
        raise ValueError('current regular receipt file required: '+name)
    return path


def build_plan(sources,seed):
    if sha(seed)!=SEED_CURRENT:raise ValueError('corrected prepared-history identity changed')
    original={SEED:{'sha256':SEED_OLD,'bytes':len(seed)}}
    current={SEED:{'sha256':SEED_CURRENT,'bytes':len(seed)}}
    if set(sources)!=set(ORIGINALS):raise ValueError('original live receipt allowlist changed')
    bindings={}
    for name,expected in ORIGINALS.items():
        raw=sources[name]
        if sha(raw)!=expected:raise ValueError('original live receipt changed: '+name)
        original[name]={'sha256':expected,'bytes':len(raw)}
    for name,raw in sources.items():
        metadata=metadata_json(raw) if name.endswith('.npz') else raw
        value,tokens,_text=json_tokens(metadata)
        bindings[name]=selected_bindings(value,tokens,original)
    outputs={};rows=[]
    while len(outputs)<len(sources):
        ready=sorted(name for name in sources if name not in outputs
                     and all(parent in current for parent in bindings[name].values()))
        if not ready:raise ValueError('live receipt dependency cycle')
        for name in ready:
            raw=sources[name]
            if name.endswith('.npz'):
                encoded,changes,members=rewrite_npz(raw,original,current)
            else:
                encoded,changes=rewrite_json(raw,original,current);members=None
            if not changes or encoded==raw:raise ValueError('allowlisted receipt has no live rebinding: '+name)
            outputs[name]=encoded;current[name]={'sha256':sha(encoded),'bytes':len(encoded)}
            rows.append({'path':name,'original':original[name],'current':current[name],
                         'changes':changes,'nonmetadata_npy_sha256':members})
    audit={'schema':'NSC-PREPARED-HISTORY-REBIND-v1','baseline':BASELINE,
           'original_source_context':SOURCE_CONTEXT,'seed':{'path':SEED,'original':original[SEED],'current':current[SEED]},
           'scope':'Eighteen live verification/input carriers; original execution capture, historical campaigns and gate proof receipts excluded.',
           'field_propagations':0,'rows':rows,'numeric_json_tokens_preserved':True,
           'nonmetadata_npz_members_byte_identical':True,'execution_signatures_preserved':True}
    return outputs,audit


def plan(root=ROOT):
    root=Path(root).resolve()
    seed=regular_file(root,SEED).read_bytes()
    sources={name:regular_file(root,name).read_bytes() for name in ORIGINALS}
    return build_plan(sources,seed)


def authenticated_originals(root):
    """Read fixed original bytes for authentication only, never as a fallback.

    Only the eighteen pinned paths at BASELINE are accepted. Numeric members
    are compared as bytes; no historical array enters a scientific operation.
    """
    root=Path(root).resolve();environment=os.environ.copy()
    for key in ('GIT_DIR','GIT_WORK_TREE','GIT_INDEX_FILE','GIT_COMMON_DIR'):
        environment.pop(key,None)
    objects=root/'archive/source-history/objects'
    if objects.is_dir():environment['GIT_ALTERNATE_OBJECT_DIRECTORIES']=str(objects)
    sources={}
    for name,expected in ORIGINALS.items():
        try:
            raw=subprocess.check_output(
                ['git','--no-pager','-C',str(root),'show','--no-textconv',BASELINE+':'+name],
                env=environment,stderr=subprocess.PIPE,timeout=30)
        except (OSError,subprocess.SubprocessError) as error:
            raise ValueError('original receipt authentication unavailable: '+name) from error
        if sha(raw)!=expected:raise ValueError('original receipt authentication changed: '+name)
        sources[name]=raw
    return sources


def atomic_write(path,raw):
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix='.receipt-rebind-',suffix='.tmp',delete=False) as stream:
            temporary=Path(stream.name);stream.write(raw)
        os.replace(temporary,path)
    finally:
        if temporary is not None and temporary.exists():temporary.unlink()


def verify_audit(root,audit):
    root=Path(root).resolve()
    seed=regular_file(root,SEED).read_bytes()
    outputs,expected=build_plan(authenticated_originals(root),seed)
    # Reconstruct both identities and every permitted field change from fixed
    # originals. The mutable ledger cannot authorize a refreshed numeric value,
    # source capture, arbitrary hash or claimed preservation flag.
    encoded=lambda value:json.dumps(value,sort_keys=True,allow_nan=False,separators=(',',':'))
    if encoded(audit)!=encoded(expected):raise ValueError('rebind audit differs from authenticated transformation')
    for name,expected_bytes in outputs.items():
        raw=regular_file(root,name).read_bytes()
        if raw!=expected_bytes:raise ValueError('rebound current bytes differ from authenticated transformation: '+name)
    return expected


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--apply',action='store_true');mode.add_argument('--verify',action='store_true')
    parser.add_argument('--audit-dir',type=Path,default=ROOT/'.build/receipt-rebind-20261014')
    parser.add_argument('--ledger',type=Path,default=ROOT/'results/metadata/prepared-history/receipt-rebind.json')
    args=parser.parse_args();ledger=args.ledger
    if args.verify:
        audit=json.loads(ledger.read_bytes());verify_audit(ROOT,audit)
    else:
        outputs,audit=plan()
        if args.apply:
            for name,raw in outputs.items():atomic_write(ROOT/name,raw)
            atomic_write(ledger,(json.dumps(audit,indent=2)+'\n').encode())
            verify_audit(ROOT,audit)
        else:atomic_write(args.audit_dir/'receipt-rebind-plan.json',(json.dumps(audit,indent=2)+'\n').encode())
    print(json.dumps({'files':len(audit['rows']),'npz':sum(row['path'].endswith('.npz') for row in audit['rows']),
        'metadata_changes':sum(len(row['changes']) for row in audit['rows']),
        'numeric_json_tokens_preserved':True,'nonmetadata_npz_members_byte_identical':True,
        'execution_signatures_preserved':True,'applied':args.apply},indent=2))


if __name__=='__main__':main()
