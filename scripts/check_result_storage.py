#!/usr/bin/env python3
"""Reject new/changed oversized result JSON, preserving exact legacy bytes."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess

from summarize_result import encoded,load_json,relative_file,sha256

ROOT=Path(__file__).resolve().parents[1]
BASELINE="docs/storage/result-payload-baseline.json"
SCHEMA="nsc-result-payload-baseline-v1"
LIMIT=1024*1024


def tracked_results(root):
    proc=subprocess.run(["git","ls-files","-z","--","results"],cwd=root,check=True,stdout=subprocess.PIPE)
    return [v.decode("utf-8") for v in proc.stdout.split(b"\x00") if v and v.lower().endswith(b".json")]


def create_baseline(path=BASELINE,*,root=ROOT):
    """Creation-only inventory of currently present, tracked oversized JSON."""
    output,_=relative_file(root,path,must_exist=False)
    if output.exists(): raise FileExistsError("storage baseline already exists")
    entries={}
    for name in sorted(tracked_results(root)):
        source,relative=relative_file(root,name,must_exist=False)
        if not source.exists(): continue  # historical sparse checkout is not regenerated
        if not source.is_file(): raise ValueError("tracked result is not a file")
        if source.stat().st_size>LIMIT:
            raw=source.read_bytes(); entries[relative]=dict(bytes=len(raw),sha256=sha256(raw))
    record=dict(schema=SCHEMA,limit_bytes=LIMIT,scope="byte-identical preexisting git-tracked large result JSON only",
        entries=entries)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open("xb") as stream: stream.write(encoded(record))
    return record


def _baseline(path,root):
    source,_=relative_file(root,path)
    data=load_json(source.read_bytes())
    if (data.get("schema")!=SCHEMA or data.get("limit_bytes")!=LIMIT or
            data.get("scope")!="byte-identical preexisting git-tracked large result JSON only" or
            not isinstance(data.get("entries"),dict)):
        raise ValueError("storage baseline declaration mismatch")
    for name,entry in data["entries"].items():
        _,relative=relative_file(root,name,must_exist=False)
        if (not relative.startswith("results/") or not relative.lower().endswith(".json") or
                not isinstance(entry,dict) or set(entry)!={"bytes","sha256"} or
                type(entry["bytes"]) is not int or entry["bytes"]<=LIMIT or
                not isinstance(entry["sha256"],str) or len(entry["sha256"])!=64 or
                any(c not in "0123456789abcdef" for c in entry["sha256"])):
            raise ValueError("invalid legacy payload binding: "+name)
    return data


def check(root=ROOT,*,baseline=BASELINE):
    data=_baseline(baseline,root); checked=0; legacy=0; failures=[]
    directory=Path(root)/"results"
    for source in sorted(p for p in directory.rglob("*") if p.suffix.lower()==".json"):
        source,relative=relative_file(root,source)
        checked+=1
        size=source.stat().st_size
        if relative in data["entries"]:
            raw=source.read_bytes(); expected=data["entries"][relative]
            if len(raw)!=expected["bytes"] or sha256(raw)!=expected["sha256"]:
                failures.append(dict(path=relative,bytes=size,reason="changed legacy payload"))
            else: legacy+=1
        elif size>LIMIT:
            failures.append(dict(path=relative,bytes=size,reason="new oversized result payload"))
    if failures:
        error=ValueError("result storage policy failed: "+json.dumps(failures,sort_keys=True))
        error.violations=failures
        raise error
    return dict(checked_json=checked,unchanged_legacy_payloads=legacy,limit_bytes=LIMIT,
        ignored_raw_directory=".build/results",scientific_replays=0)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline",default=BASELINE)
    parser.add_argument("--write-baseline",action="store_true",help="one-time tracked legacy inventory; creation-only")
    args=parser.parse_args()
    if args.write_baseline:
        value=create_baseline(args.baseline)
        result=dict(created=args.baseline,legacy_payloads=len(value["entries"]),limit_bytes=LIMIT)
    else: result=check(baseline=args.baseline)
    print(json.dumps(result,sort_keys=True))


if __name__=="__main__": main()
