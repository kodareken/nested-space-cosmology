"""Project terminal formation records into a compact, domain-limited summary."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
RUNNER='scripts/evolve_nsc_atlas_valve_counterrotation.py'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_run(path):
    path=(ROOT/path).resolve();raw=path.read_bytes();record=json.loads(raw)
    if not path.is_relative_to(ROOT/'.build') or record['status'] in ('initializing','running'):
        raise ValueError('terminal local formation record required')
    if not record.get('source_bound_to_commit') or not record.get('producer_source_unchanged_during_run'):
        raise ValueError('source-bound unchanged formation producer required')
    for name,digest in record['source_sha256'].items():
        old=subprocess.run(['git','show',record['producer_commit']+':'+name],cwd=ROOT,capture_output=True)
        if old.returncode or hashlib.sha256(old.stdout).hexdigest()!=digest:raise ValueError('historical source failed: '+name)
        if name!=RUNNER and sha(ROOT/name)!=digest:raise ValueError('current physical decoder differs: '+name)
    if sha(ROOT/record['initial_record'])!=record['initial_record_sha256']:raise ValueError('initial preparation changed')
    for member in [dict(path=record['state_cache'],sha256=record['state_cache_sha256']),*record['dense_history_chunks']]:
        p=(ROOT/member['path']).resolve()
        if not p.is_relative_to(ROOT/'.build') or sha(p)!=member['sha256']:raise ValueError('raw trajectory hash failed')
    rows=[]
    for row in record['snapshots']:
        rows.append(dict(proper=row['observer_proper_duration'],coordinate_t=row['coordinate_t'],
            area_radius=row['neck']['area_radius'],theta_plus_over_H_range=row['neck']['theta_plus_over_H_range'],
            theta_minus_over_H_range=row['neck']['theta_minus_over_H_range'],
            future_trapped_interface_sampled=row['neck']['future_trapped_sampled'],
            admitted=row['sampled_admitted'],global_relative_constraints=row['global_max_relative_constraints'],
            local_relative_constraints=row['neck_neighborhood_max_relative_constraints'],
            observer_unit_norm_error=row['observer_unit_norm_error'],
            scalar_last_two_shell_power_fraction=row['last_two_scalar_shell_power_fraction'],
            global_matter_energy=row['total_slice_normal_matter_energy'],global_volume=row['proper_slice_volume'],
            excitation_work_residual=row['excitation_work_residual']))
    return dict(record=path.relative_to(ROOT).as_posix(),record_sha256=hashlib.sha256(raw).hexdigest(),
        producer_commit=record['producer_commit'],source_sha256=record['source_sha256'],
        degree=record['degree'],status=record['status'],stop_reason=record.get('stop_reason'),
        initial_record=record['initial_record'],initial_record_sha256=record['initial_record_sha256'],
        model_inputs=record['model_inputs'],threads=record['threads'],
        requested_proper_endpoint=record['requested_proper_endpoint'],temporal_tolerances=record['temporal_tolerances'],
        elapsed_seconds=record['elapsed_seconds'],last_admitted_coordinate_t=record.get('last_admitted_coordinate_t'),
        actual_last_accepted_proper=record.get('actual_last_accepted_proper'),
        raw_state_cache=dict(path=record['state_cache'],sha256=record['state_cache_sha256']),
        raw_history=record['dense_history_chunks'],snapshots=rows)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,action='append',required=True)
    p.add_argument('--output',type=Path,default=ROOT/'results/summaries/nsc-atlas-valve-counterrotation-evolution-v1.json')
    args=p.parse_args();runs=[read_run(path) for path in args.run];comparisons=[]
    if len(runs)==2:
        a={round(v['proper'],8):v for v in runs[0]['snapshots'] if v['admitted']}
        b={round(v['proper'],8):v for v in runs[1]['snapshots'] if v['admitted']}
        for t in sorted(a.keys()&b.keys()):
            comparisons.append(dict(proper=t,area_radius_spread=abs(a[t]['area_radius']-b[t]['area_radius']),
                fine_area_radius=b[t]['area_radius'],both_interfaces_future_trapped=a[t]['future_trapped_interface_sampled'] and b[t]['future_trapped_interface_sampled']))
    output=dict(schema='nsc-atlas-valve-counterrotation-evolution-v1',
        scientific_owner='docs/nsc-atlas-valve-formation.md',runs=runs,matched_admitted_events=comparisons,
        interpretation=dict(regular_child_proved=False,geometric_scope='Nonspherical fixed-interface area radius, not a spherical areal coordinate.',
            numerical_scope='Sampled admission and two-resolution empirical differences; no continuum error bound or measured convergence order.',
            receiving_scope='Local counterflow is not a net supply proof for a continuously tracked child boundary.',
            stops=[dict(degree=r['degree'],reason=r['stop_reason'],
                last_admitted_proper=max(v['proper'] for v in r['snapshots'] if v['admitted'])) for r in runs],
            SI_or_cosmological_prediction=False),
        producer_sha256=sha(__file__),reproduction='python scripts/run_python.py scripts/summarize_nsc_atlas_valve_counterrotation.py --run .build/valve-counterrotation-evolution-L24-v1/run.json --run .build/valve-counterrotation-evolution-L32-v1/run.json')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(args.output),runs=[dict(degree=r['degree'],status=r['status'],stop=r['stop_reason']) for r in runs],
        shared_proper_end=comparisons[-1]['proper'] if comparisons else None,regular_child_proved=False),indent=2))


if __name__=='__main__':main()
