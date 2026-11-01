"""Project terminal formation records into a compact, domain-limited summary."""
import argparse
import hashlib
from itertools import combinations
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
RUNNER='scripts/evolve_nsc_atlas_valve_counterrotation.py'
TERMINAL_STATUSES={'endpoint_reached','incomplete','initial_only','initial_not_admitted'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_run(path):
    path=(ROOT/path).resolve();raw=path.read_bytes();record=json.loads(raw)
    if not path.is_relative_to(ROOT/'.build') or record.get('status') not in TERMINAL_STATUSES:
        raise ValueError('terminal local formation record required')
    if not record.get('source_bound_to_commit') or not record.get('producer_source_unchanged_during_run'):
        raise ValueError('source-bound unchanged formation producer required')
    for name,digest in record['source_sha256'].items():
        old=subprocess.run(['git','show',record['producer_commit']+':'+name],cwd=ROOT,capture_output=True)
        if old.returncode or hashlib.sha256(old.stdout).hexdigest()!=digest:raise ValueError('historical source failed: '+name)
        if name!=RUNNER and sha(ROOT/name)!=digest:raise ValueError('current physical decoder differs: '+name)
    if sha(ROOT/record['initial_record'])!=record['initial_record_sha256']:raise ValueError('initial preparation changed')
    for member in [dict(path=record['state_cache'],sha256=record['state_cache_sha256']),
                   *record.get('physical_history_prefix',[]),*record['dense_history_chunks']]:
        p=(ROOT/member['path']).resolve()
        if not p.is_relative_to(ROOT/'.build') or sha(p)!=member['sha256']:raise ValueError('raw trajectory hash failed')
    resume=record.get('resume')
    if resume:
        parent=(ROOT/resume['parent_run']).resolve()
        if not parent.is_relative_to(ROOT/'.build') or sha(parent)!=resume['parent_run_sha256']:
            raise ValueError('resume parent terminal record hash failed')
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
            integrated_excitation_work=row['integrated_excitation_work'],
            excitation_work_residual=row['excitation_work_residual']))
    return dict(record=path.relative_to(ROOT).as_posix(),record_sha256=hashlib.sha256(raw).hexdigest(),
        producer_commit=record['producer_commit'],source_sha256=record['source_sha256'],
        degree=record['degree'],status=record['status'],stop_reason=record.get('stop_reason'),
        initial_record=record['initial_record'],initial_record_sha256=record['initial_record_sha256'],
        model_inputs=record['model_inputs'],threads=record['threads'],
        requested_proper_endpoint=record['requested_proper_endpoint'],temporal_tolerances=record['temporal_tolerances'],
        elapsed_seconds=record['elapsed_seconds'],last_admitted_coordinate_t=record.get('last_admitted_coordinate_t'),
        first_coordinate_t=rows[0]['coordinate_t'],first_proper=rows[0]['proper'],
        last_admitted_proper=max((v['proper'] for v in rows if v['admitted']),default=None),
        actual_last_accepted_coordinate_t=record.get('actual_last_accepted_coordinate_t'),
        actual_last_accepted_proper=record.get('actual_last_accepted_proper'),
        resume=resume,physical_history_prefix=record.get('physical_history_prefix',[]),
        excitation_work_origin_energy=record.get('excitation_work_origin_energy',
            record['snapshots'][0]['normal_energy_minus_fixed_vacuum']),
        raw_state_cache=dict(path=record['state_cache'],sha256=record['state_cache_sha256']),
        raw_history=record['dense_history_chunks'],snapshots=rows)


def matched_admitted_events(runs):
    """Compare every pair without stitching or extending either run's domain."""
    comparisons=[]
    for first,second in combinations(runs,2):
        if (first['initial_record_sha256']!=second['initial_record_sha256'] or
            first['model_inputs']!=second['model_inputs']):
            raise ValueError('matched events require the same preparation and physical model')
        a={round(v['proper'],8):v for v in first['snapshots'] if v['admitted']}
        b={round(v['proper'],8):v for v in second['snapshots'] if v['admitted']}
        linked=(second.get('resume') or {}).get('parent_run')==first['record'] or (
            first.get('resume') or {}).get('parent_run')==second['record']
        kind=('continuation_handoff' if linked else 'same_grid_comparison') if first['degree']==second['degree'] else 'resolution_comparison'
        for t in sorted(a.keys()&b.keys()):
            fine=(a[t] if first['degree']>second['degree'] else b[t]) if first['degree']!=second['degree'] else None
            comparisons.append(dict(run_a=first['record'],run_b=second['record'],
                degree_a=first['degree'],degree_b=second['degree'],comparison_kind=kind,
                proper=t,coordinate_t_a=a[t]['coordinate_t'],coordinate_t_b=b[t]['coordinate_t'],
                area_radius_a=a[t]['area_radius'],area_radius_b=b[t]['area_radius'],
                area_radius_spread=abs(a[t]['area_radius']-b[t]['area_radius']),
                fine_area_radius=fine['area_radius'] if fine else None,
                both_interfaces_future_trapped=a[t]['future_trapped_interface_sampled'] and b[t]['future_trapped_interface_sampled']))
    return comparisons


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,action='append',required=True)
    p.add_argument('--output',type=Path,default=ROOT/'results/summaries/nsc-atlas-valve-counterrotation-evolution-v1.json')
    args=p.parse_args();runs=[read_run(path) for path in args.run]
    if len({r['record'] for r in runs})!=len(runs):raise ValueError('distinct terminal run records required')
    comparisons=matched_admitted_events(runs)
    destination=args.output.resolve()
    output_arg=destination.relative_to(ROOT).as_posix() if destination.is_relative_to(ROOT) else str(destination)
    command=['python','scripts/run_python.py','scripts/summarize_nsc_atlas_valve_counterrotation.py']
    for r in runs:command.extend(['--run',r['record']])
    command.extend(['--output',output_arg])
    output=dict(schema='nsc-atlas-valve-counterrotation-evolution-v1',
        scientific_owner='docs/nsc-atlas-valve-formation.md',runs=runs,matched_admitted_events=comparisons,
        interpretation=dict(regular_child_proved=False,geometric_scope='Nonspherical fixed-interface area radius, not a spherical areal coordinate.',
            numerical_scope='Pairwise matched admitted events only; resolution differences are empirical, and same-grid continuation handoffs are not convergence tests. Each run retains its original stop boundary.',
            receiving_scope='Local counterflow is not a net supply proof for a continuously tracked child boundary.',
            stops=[dict(record=r['record'],degree=r['degree'],status=r['status'],reason=r['stop_reason'],
                first_proper=r['first_proper'],last_admitted_proper=r['last_admitted_proper'],
                actual_last_accepted_proper=r['actual_last_accepted_proper']) for r in runs],
            SI_or_cosmological_prediction=False),
        producer_sha256=sha(__file__),reproduction_argv=command,reproduction=subprocess.list2cmdline(command))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(args.output),runs=[dict(degree=r['degree'],status=r['status'],stop=r['stop_reason']) for r in runs],
        shared_proper_end=max((v['proper'] for v in comparisons),default=None),regular_child_proved=False),indent=2))


if __name__=='__main__':main()
