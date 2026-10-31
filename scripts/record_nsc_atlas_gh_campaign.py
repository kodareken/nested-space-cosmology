"""Preserve source-bound GH attempts without admitting unresolved predictions."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

ROOT=Path(__file__).resolve().parents[1]


def source_binding(record):
    commit=record.get('producer_commit','')
    if not re.fullmatch('[0-9a-f]{40}',commit):return False
    if not record.get('producer_source_unchanged_during_run') or not record.get('source_bound_to_commit'):return False
    for name,expected in record['source_sha256'].items():
        stored=subprocess.run(['git','show',commit+':'+name],cwd=ROOT,capture_output=True)
        if stored.returncode or hashlib.sha256(stored.stdout).hexdigest()!=expected:return False
    return True


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,action='append',required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();output=(ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('fresh repository output required')
    cases=[]
    for name in args.input:
        path=(ROOT/name).resolve()
        if not path.is_relative_to(ROOT):raise ValueError('repository input required')
        raw=path.read_bytes();record=json.loads(raw)
        bound=source_binding(record)
        cases.append(dict(input_sha256=hashlib.sha256(raw).hexdigest(),source_binding_verified=bound,
            disposition='unqualified_single_resolution' if bound else 'unverified_source_binding',
            physical_prediction_admitted=False,run=record))
    result=dict(schema='nsc-atlas-gh-campaign-v1',
        objective='Constrained nonlinear continuation of the unchanged nonspherical inherited source',
        qualification='Endpoint reachability alone is insufficient. Initial/final constraints, temporal and spatial convergence and resolved observables are required.',
        cases=cases,convergence_decided=False,
        source_sha256={'scripts/record_nsc_atlas_gh_campaign.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(cases=len(cases),verified_bindings=sum(c['source_binding_verified'] for c in cases),physical_predictions_admitted=0)))


if __name__=='__main__':main()
