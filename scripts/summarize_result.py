#!/usr/bin/env python3
"""Create a compact, byte-bound result summary without rerunning science.

Raw payloads belong in ignored .build/results. A summary preserves provenance,
the exact regeneration argv, and all selected structure/string values. Long
numeric arrays are replaced by explicit content hashes and scalar statistics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-compact-result-v1"
MAX_SUMMARY_BYTES = 1024 * 1024
ARRAY_LIMIT = 16
PROVENANCE_KEYS = (
    "schema", "status", "scope", "runtime", "execution", "input", "inputs",
    "source_bindings", "audit_source_bindings", "producer_bindings", "manifest",
    "limits", "settings", "input_identity", "cache_identity",
)


def encoded(value):
    return (json.dumps(value,sort_keys=True,indent=2,ensure_ascii=False,allow_nan=False)+"\n").encode("utf-8")


def sha256(raw):
    return hashlib.sha256(raw).hexdigest()


def _object(pairs):
    out={}
    for key,value in pairs:
        if key in out: raise ValueError("duplicate JSON key: "+key)
        out[key]=value
    return out


def load_json(raw):
    def constant(value): raise ValueError("nonfinite JSON constant: "+value)
    value=json.loads(raw,object_pairs_hook=_object,parse_constant=constant)
    encoded(value)  # also rejects overflowed floating literals, such as1e999
    return value


def relative_file(root,path,*,must_exist=True):
    root=Path(root).resolve(); supplied=Path(path)
    candidate=supplied if supplied.is_absolute() else root/supplied
    if ".." in supplied.parts: raise ValueError("parent traversal is not allowed")
    try: relative=candidate.relative_to(root)
    except ValueError as error: raise ValueError("path must remain inside the repository") from error
    if not relative.parts or relative.parts[0]==".git": raise ValueError("invalid result path")
    cursor=root
    for part in relative.parts:
        cursor=cursor/part
        if cursor.is_symlink(): raise ValueError("symlink paths are not accepted")
    if candidate.resolve().is_relative_to(root) is False: raise ValueError("path escaped repository")
    if must_exist and not candidate.is_file(): raise ValueError("missing regular file: "+relative.as_posix())
    return candidate,relative.as_posix()


def _numeric_flat(value):
    if isinstance(value,list):
        result=[]
        for item in value:
            part=_numeric_flat(item)
            if part is None: return None
            result.extend(part)
        return result
    if isinstance(value,bool) or not isinstance(value,(int,float)): return None
    if isinstance(value,float) and not math.isfinite(value): raise ValueError("nonfinite numeric array")
    return [value]


def compact(value):
    """Preserve structure; replace numeric lists with more than16 scalars.

    The array hash is SHA256 of this module's canonical JSON encoding of the
    complete original list. It is distinct from the raw payload byte hash.
    Empty lists, booleans and strings are retained exactly.
    """
    if isinstance(value,dict): return {key:compact(item) for key,item in value.items()}
    if isinstance(value,list):
        flat=_numeric_flat(value)
        if flat is not None and len(flat)>ARRAY_LIMIT:
            return dict(compacted_numeric_array=True,length=len(value),scalar_count=len(flat),
                sha256=sha256(encoded(value)),minimum=min(flat),maximum=max(flat))
        return [compact(item) for item in value]
    return value


def select(value,pointer):
    """RFC6901-style JSON Pointer; empty selects the complete payload."""
    if not isinstance(pointer,str) or (pointer and not pointer.startswith("/")):
        raise ValueError("selection must be an empty or slash-prefixed JSON Pointer")
    if not pointer: return value
    current=value
    for part in pointer[1:].split("/"):
        import re
        if re.search(r"~(?![01])",part): raise ValueError("invalid JSON Pointer escape")
        key=part.replace("~1","/").replace("~0","~")
        if isinstance(current,list):
            if not key.isdecimal() or (len(key)>1 and key.startswith("0")):
                raise ValueError("invalid JSON Pointer list index")
            index=int(key)
            if index>=len(current): raise ValueError("JSON Pointer index unavailable")
            current=current[index]
        elif isinstance(current,dict) and key in current: current=current[key]
        else: raise ValueError("JSON Pointer selection unavailable")
    return current


def _bindings(paths,root):
    result={}
    for name in paths:
        path,relative=relative_file(root,name)
        if relative in result: raise ValueError("duplicate binding: "+relative)
        raw=path.read_bytes(); result[relative]=dict(sha256=sha256(raw),bytes=len(raw))
    return result


def _argv(value):
    if (not isinstance(value,list) or not value or
            any(not isinstance(v,str) or not v or "\x00" in v for v in value)):
        raise ValueError("regeneration argv must be a nonempty list of nonempty strings")
    return list(value)


def create_summary(input_path,output_path,*,argv,producer_paths,root=ROOT,selection=None,input_paths=()):
    """Creation-only summary; hashes bytes and preserves reported provenance.

    Explicit producer/input bindings are checked as local bytes. Original
    scientific authentication metadata is retained, not re-executed or promoted
    to an independent validation. No command in argv is executed.
    """
    source,source_relative=relative_file(root,input_path)
    output,output_relative=relative_file(root,output_path,must_exist=False)
    if output.exists(): raise FileExistsError("summary output is creation-only")
    if Path(output_relative).parts[0]!="results" or output.suffix.lower()!=".json":
        raise ValueError("tracked summaries must be JSON under results/")
    if not producer_paths: raise ValueError("at least one explicit producer binding is required")
    argv=_argv(argv); pointer="" if selection is None else selection
    raw=source.read_bytes(); payload=load_json(raw)
    producers=_bindings(producer_paths,root); inputs=_bindings(input_paths,root)
    provenance={key:compact(payload[key]) for key in PROVENANCE_KEYS if isinstance(payload,dict) and key in payload}
    result=dict(schema=SCHEMA,raw_payload=dict(path=source_relative,sha256=sha256(raw),bytes=len(raw)),
        regeneration=dict(argv=argv,cwd=".",executed_by_summary=False),producer_bindings=producers,
        input_bindings=inputs,selection=pointer,reported_provenance=provenance,
        authentication=dict(scope="raw and explicitly listed local byte bindings only",
            scientific_validation_executed=False,reported_provenance_revalidated=False),
        data=compact(select(payload,pointer)))
    body=encoded(result)
    if len(body)>MAX_SUMMARY_BYTES:
        raise ValueError("summary exceeds1MiB; select a smaller scientific object explicitly")
    if (source.read_bytes()!=raw or _bindings(producer_paths,root)!=producers or _bindings(input_paths,root)!=inputs):
        raise ValueError("summary input changed while being read")
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open("xb") as stream: stream.write(body)
    return result


def check_summary(path,*,root=ROOT):
    """Authenticate original producers and current inputs, without replay.

    An explicitly reviewed source successor can authenticate its original
    Git bytes. Creation still records current bytes. This does not transfer
    the old numerical result to the changed implementation.
    """
    target,_=relative_file(root,path)
    body=target.read_bytes()
    if len(body)>MAX_SUMMARY_BYTES: raise ValueError("summary exceeds1MiB")
    saved=load_json(body)
    if saved.get("schema")!=SCHEMA: raise ValueError("unsupported compact summary")
    if saved.get("authentication")!=dict(scope="raw and explicitly listed local byte bindings only",
        scientific_validation_executed=False,reported_provenance_revalidated=False):
        raise ValueError("summary authentication scope differs")
    regen=saved.get("regeneration",{})
    _argv(regen.get("argv"))
    if regen.get("cwd")!="." or regen.get("executed_by_summary") is not False:
        raise ValueError("invalid regeneration declaration")
    if not saved.get("producer_bindings"): raise ValueError("missing producer bindings")
    from validation_source_successors import verify_saved_source
    bindings=saved['producer_bindings'];historical=[]
    if not isinstance(bindings,dict): raise ValueError('producer_bindings mismatch')
    for name,expected in bindings.items():
        source,relative=relative_file(root,name)
        if (not isinstance(expected,dict) or set(expected)!={'sha256','bytes'}
                or type(expected['bytes']) is not int or expected['bytes']<0):
            raise ValueError('invalid producer binding')
        original=verify_saved_source(relative,expected['sha256'],root)
        if len(original)!=expected['bytes']: raise ValueError('producer_bindings byte count mismatch')
        if sha256(source.read_bytes())!=expected['sha256']: historical.append(relative)
    bindings=saved.get('input_bindings')
    if not isinstance(bindings,dict) or _bindings(list(bindings),root)!=bindings:
        raise ValueError('input_bindings mismatch')
    raw=saved["raw_payload"]
    original,_=relative_file(root,raw["path"],must_exist=False)
    if (type(raw.get("bytes")) is not int or raw["bytes"]<0 or not isinstance(raw.get("sha256"),str)
            or len(raw["sha256"])!=64 or any(c not in "0123456789abcdef" for c in raw["sha256"])):
        raise ValueError("invalid raw payload binding")
    available=original.is_file()
    if original.exists() and not available: raise ValueError("raw payload is not a file")
    if available:
        raw_bytes=original.read_bytes()
        if len(raw_bytes)!=raw["bytes"] or sha256(raw_bytes)!=raw["sha256"]:
            raise ValueError("raw payload binding mismatch")
        payload=load_json(raw_bytes)
        provenance={key:compact(payload[key]) for key in PROVENANCE_KEYS if isinstance(payload,dict) and key in payload}
        if compact(select(payload,saved["selection"]))!=saved["data"] or provenance!=saved["reported_provenance"]:
            raise ValueError("compact data differs from raw payload")
    return dict(checked=str(path),raw_available=available,raw_bytes=raw["bytes"],
        producers=len(saved["producer_bindings"]),inputs=len(saved["input_bindings"]),scientific_replays=0,
        producer_context='reviewed-successors' if historical else 'current',
        original_producer_paths=historical)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--argv-json",help="exact argument-vector JSON; never executed here")
    parser.add_argument("--producer",action="append",default=[])
    parser.add_argument("--input-binding",action="append",default=[])
    parser.add_argument("--select",default="",help="JSON Pointer into scientific data")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    if args.check: receipt=check_summary(args.output)
    else:
        if args.input is None or args.argv_json is None: parser.error("creation requires --input and --argv-json")
        record=create_summary(args.input,args.output,argv=load_json(args.argv_json),
            producer_paths=args.producer,input_paths=args.input_binding,selection=args.select)
        receipt=dict(created=str(args.output),raw_payload=record["raw_payload"],summary_bytes=len(encoded(record)))
    print(json.dumps(receipt,sort_keys=True))


if __name__=="__main__": main()
