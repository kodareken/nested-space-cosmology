"""Authenticate and compactly project the completed directional0–2 runs.

No evolution driver is imported, and no field/grid calculation is performed.
Historic runner bytes are checked in Git; current physical-package bytes
must match the run. Later runner-only resume changes do not rewrite history.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
GIT_HASHES={}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    path=path.resolve()
    return path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else str(path)


def authenticate(path):
    path=path.resolve();record=json.loads(path.read_text(encoding='utf-8'))
    if record.get('schema')!='nsc-atlas-valve-directional-evolution-v1':raise ValueError('directional run schema required')
    if not record['source_bound_to_commit'] or not record['source_unchanged']:
        raise ValueError('run lacks unchanged commit-bound source')
    commit=record['producer_commit'];current_runner={}
    for name,digest in record['source_sha256'].items():
        key=(commit,name)
        if key not in GIT_HASHES:
            blob=subprocess.check_output(['git','show',commit+':'+name],cwd=ROOT)
            GIT_HASHES[key]=hashlib.sha256(blob).hexdigest()
        if GIT_HASHES[key]!=digest:raise ValueError('historical producer/source SHA mismatch: '+name)
        current=sha(ROOT/name)
        if name.startswith('src/') and current!=digest:raise ValueError('current physical package differs: '+name)
        if name.startswith('scripts/'):current_runner[name]=current==digest
    cache=ROOT/record['state_cache']
    if sha(cache)!=record['state_cache_sha256']:raise ValueError('sparse state cache SHA mismatch')
    snapshots=record['snapshots']
    with np.load(cache,allow_pickle=False) as raw:
        times=raw['times'];states=raw['states']
        expected=np.array([s['proper_time'] for s in snapshots])
        if (times.shape!=expected.shape or not np.allclose(times,expected,rtol=0,atol=2e-13) or
                states.shape!=(len(times),8,record['modes']) or not np.all(np.isfinite(states))):
            raise ValueError('sparse states do not match the recorded snapshots')
    initial=ROOT/record['initial_record']
    if sha(initial)!=record['initial_record_sha256']:raise ValueError('preserved initial record SHA mismatch')
    initial_record=json.loads(initial.read_text(encoding='utf-8'));initial_raw=ROOT/initial_record['raw_cache']['path']
    if sha(initial_raw)!=initial_record['raw_cache']['sha256']:raise ValueError('preserved initial cache SHA mismatch')
    endpoint=float(record['proper_endpoint_requested'])
    if (record['status']!='endpoint_reached' or record['stop_reason'] is not None or
            abs(record['last_admitted_proper']-endpoint)>1e-12 or
            abs(record['actual_last_accepted_proper']-endpoint)>1e-12 or
            not all(s['sampled_admitted'] for s in snapshots)):
        raise ValueError('projection requires the completed fully sampled-admitted finite run')
    return record,dict(run_record=relative(path),run_record_sha256=sha(path),producer_commit=commit,
        source_sha256=record['source_sha256'],historical_Git_source_verified=True,current_physical_package_matches=True,
        current_runner_matches_historical=current_runner,
        sparse_cache=dict(path=record['state_cache'],sha256=record['state_cache_sha256']),
        initial_record=dict(path=record['initial_record'],sha256=record['initial_record_sha256'],raw_cache=initial_record['raw_cache']))


def run_projection(record,authentication):
    snaps=record['snapshots']
    maximum=lambda key:max(s[key] for s in snaps)
    controls=dict(sampled_max_relative_H=max(s['global_max_relative_constraints']['H'] for s in snaps),
        sampled_max_relative_M=max(s['global_max_relative_constraints']['M'] for s in snaps),
        max_Hamiltonian_RMS_over_H2=maximum('Hamiltonian_RMS_over_H2'),
        max_momentum_RMS_over_H2=maximum('momentum_RMS_over_H2'),
        minimum_sampled_radial_scale=min(s['min_radial_scale'] for s in snaps),
        minimum_sampled_angular_scale=min(s['min_angular_scale'] for s in snaps),
        max_abs_receiving_cap_energy_ledger_residual=max(abs(s['fixed_receiving_cap']['energy_ledger_residual']) for s in snaps),
        max_abs_donor_cap_energy_ledger_residual=max(abs(s['fixed_donor_cap']['energy_ledger_residual']) for s in snaps))
    return dict(modes=record['modes'],authentication=authentication,
        domain=dict(proper_time=[snaps[0]['proper_time'],snaps[-1]['proper_time']],
                    requested_endpoint=record['proper_endpoint_requested'],last_admitted=record['last_admitted_proper'],
                    status=record['status'],all_recorded_snapshots_sampled_admitted=True,
                    admission='Finite spatial samples at the recorded times, ADM-term2% scale; not a continuum error bound.'),
        temporal_tolerances=record['temporal_tolerances'],reflection_projection=record['reflection_projection'],
        controls=controls,initial=snaps[0],endpoint=snaps[-1],
        compute=dict(threads=record['threads'],elapsed_seconds=record['elapsed_seconds'],rhs_calls=record['rhs_calls']))


def minimum_radius(snapshot):return None if snapshot['minimum'] is None else snapshot['minimum']['radius']


def band_edges(snapshot):
    b=snapshot['trapped_band']
    return None if b is None else [b['u_left'],b['u_right']]


def project(coarse_path,fine_path):
    start=time.perf_counter();producer_before=sha(__file__)
    coarse,cauth=authenticate(coarse_path);fine,fauth=authenticate(fine_path)
    if coarse['initial_record_sha256']!=fine['initial_record_sha256']:
        raise ValueError('runs do not start from the same directional Cauchy record')
    if coarse['source_sha256']!=fine['source_sha256']:
        raise ValueError('the two finite runs do not share the same physical/runner sources')
    a=coarse['snapshots'];b=fine['snapshots']
    if len(a)!=len(b) or any(abs(x['proper_time']-y['proper_time'])>1e-12 for x,y in zip(a,b)):
        raise ValueError('matched-time projection requires identical recorded proper times')
    table=[];maximum_spreads=dict(minimum_radius=0.,fixed_cap_energy=0.,fixed_cap_volume=0.,net_material_received=0.,band_left=0.,band_right=0.)
    for x,y in zip(a,b):
        values=dict(minimum_radius=[minimum_radius(x),minimum_radius(y)],
            fixed_cap_energy=[x['fixed_receiving_cap']['energy'],y['fixed_receiving_cap']['energy']],
            fixed_cap_volume=[x['fixed_receiving_cap']['volume'],y['fixed_receiving_cap']['volume']],
            net_material_received=[x['fixed_receiving_cap']['integrated_material_received'],y['fixed_receiving_cap']['integrated_material_received']])
        spreads={key:None if None in pair else abs(pair[1]-pair[0]) for key,pair in values.items()}
        edges=[band_edges(x),band_edges(y)]
        edge_spreads=None if any(v is None for v in edges) else [abs(edges[1][j]-edges[0][j]) for j in (0,1)]
        for key,value in spreads.items():
            if value is not None:maximum_spreads[key]=max(maximum_spreads[key],value)
        if edge_spreads is not None:
            maximum_spreads['band_left']=max(maximum_spreads['band_left'],edge_spreads[0])
            maximum_spreads['band_right']=max(maximum_spreads['band_right'],edge_spreads[1])
        table.append(dict(proper_time=x['proper_time'],**values,trapped_band_edges=edges,
                          empirical_absolute_spreads=spreads,empirical_band_edge_spreads=edge_spreads))
    columns=['proper_time'];compact_rows=[[row['proper_time']] for row in table]
    for field in ('minimum_radius','fixed_cap_energy','fixed_cap_volume','net_material_received'):
        columns.extend([field+f'_N{coarse["modes"]}',field+f'_N{fine["modes"]}','empirical_abs_spread_'+field])
        for row,compact in zip(table,compact_rows):compact.extend(row[field]+[row['empirical_absolute_spreads'][field]])
    for j,name in enumerate(('band_left','band_right')):
        columns.extend([name+f'_N{coarse["modes"]}',name+f'_N{fine["modes"]}','empirical_abs_spread_'+name])
        for row,compact in zip(table,compact_rows):
            edge=row['trapped_band_edges'];spreads=row['empirical_band_edge_spreads']
            compact.extend([None if edge[i] is None else edge[i][j] for i in (0,1)]+[None if spreads is None else spreads[j]])
    endpoint=fine['snapshots'][-1];cap=endpoint['fixed_receiving_cap']
    result=dict(schema='nsc-atlas-valve-directional-evolution-v1',
        projection_producer=dict(path=relative(Path(__file__)),sha256=producer_before,
                                 role='Passive projection of historical15137 trajectories; no current runner imported.'),
        runs=[run_projection(coarse,cauth),run_projection(fine,fauth)],
        comparison=dict(run_order=[coarse['modes'],fine['modes']],matched_times=dict(columns=columns,rows=compact_rows),
                        maximum_observed_absolute_spreads=maximum_spreads,
                        scope='Empirical differences between these two finite runs, not error bars or a continuum/exact-Cauchy-data certificate.'),
        interpretation=dict(net_material='Signed accumulated material supply through the fixedu=0 interface; not gross crossings or energy creation.',
                            energy_ledger='Fixed-cap normal matter-energy change=net material supplied+integrated pressure work; no moving-boundary sweep is mixed into supply.',
                            migrating_minimum='An area minimum on each slice; its coordinate motion is not a material observer worldline.',
                            trapped_band='Recorded marginal roots delimit a sampled future-trapped band; no new horizon classification is performed here.'),
        verdict=dict(defined_directional_preparation=True,inherited_selection_proved=False,
                     requested_finite_endpoint_completed=True,regular_receiving_child_proved=False,
                     endpoint_receiving_cap_expanding=cap['bulk_volume_rate']>0,
                     endpoint_receiving_cap_volume_rate=cap['bulk_volume_rate'],endpoint_receiving_cap_w=cap['w'],
                     new_horizon_classification_or_NEC_gate_computed=False,
                     scope='Completed0–2 proper-time finite domain only; no singularity, later continuation or cosmological claim.'),
        projection_elapsed_seconds=time.perf_counter()-start,
        reproduction='python scripts/run_python.py scripts/summarize_nsc_atlas_valve_directional.py '
            f'--coarse-record {relative(coarse_path)} --fine-record {relative(fine_path)} '
            '--output results/summaries/nsc-atlas-valve-directional-evolution-v1.json')
    if sha(__file__)!=producer_before:raise RuntimeError('projection script changed while reading the historical records')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coarse-record',type=Path,default=ROOT/'.build/valve-directional-evolution-N64-v1/run.json')
    parser.add_argument('--fine-record',type=Path,default=ROOT/'.build/valve-directional-evolution-N96-v1/run.json')
    parser.add_argument('--output',type=Path,default=ROOT/'results/summaries/nsc-atlas-valve-directional-evolution-v1.json')
    args=parser.parse_args();result=project(args.coarse_record,args.fine_record)
    encoded=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if len(encoded.encode('utf-8'))>40*1024:raise RuntimeError('compact projection exceeds40KiB; no fields silently dropped')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(encoded,encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(args.output),bytes=len(encoded.encode('utf-8')),
                         matched_times=len(result['comparison']['matched_times']['rows']),
                         max_spreads=result['comparison']['maximum_observed_absolute_spreads'],verdict=result['verdict']),indent=2))


if __name__=='__main__':main()
