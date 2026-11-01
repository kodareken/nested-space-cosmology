"""Authenticate and compactly project terminal directional finite windows.

No evolution driver is imported, and no field/grid calculation is performed.
Historic source bytes are checked in Git. Current package/runner differences
are metadata only: this historical projection imports no physics and does
not authorize continuation with a changed right-hand side.
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


def admitted_snapshots(record):
    admitted=[];ended=False
    for snapshot in record['snapshots']:
        if snapshot['sampled_admitted'] and not ended:admitted.append(snapshot)
        elif snapshot['sampled_admitted']:raise ValueError('admission cannot silently restart after a failed sample')
        else:ended=True
    if not admitted or abs(admitted[-1]['proper_time']-record['last_admitted_proper'])>1e-12:
        raise ValueError('recorded last-admitted time differs from its contiguous admitted snapshots')
    return admitted


def authenticate(path,seen=None):
    path=path.resolve();record=json.loads(path.read_text(encoding='utf-8'))
    if record.get('schema')!='nsc-atlas-valve-directional-evolution-v1':raise ValueError('directional run schema required')
    if record.get('status') not in ('endpoint_reached','incomplete'):
        raise ValueError('terminal endpoint or incomplete record required; live runs are rejected')
    seen=set() if seen is None else set(seen)
    if path in seen:raise ValueError('cyclic resume ancestry')
    seen.add(path)
    if not record['source_bound_to_commit'] or not record['source_unchanged']:
        raise ValueError('run lacks unchanged commit-bound source')
    commit=record['producer_commit'];current_runner={};package_differences={}
    for name,digest in record['source_sha256'].items():
        key=(commit,name)
        if key not in GIT_HASHES:
            blob=subprocess.check_output(['git','show',commit+':'+name],cwd=ROOT)
            GIT_HASHES[key]=hashlib.sha256(blob).hexdigest()
        if GIT_HASHES[key]!=digest:raise ValueError('historical producer/source SHA mismatch: '+name)
        current=sha(ROOT/name) if (ROOT/name).is_file() else None
        if name.startswith('src/') and current!=digest:
            package_differences[name]=dict(archived_sha256=digest,current_sha256=current)
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
    admitted=admitted_snapshots(record)
    accepted=float(record['actual_last_accepted_proper'])
    if accepted<snapshots[-1]['proper_time']-1e-12:raise ValueError('cached snapshot lies beyond the accepted trajectory')
    if record['status']=='endpoint_reached':
        if (record['stop_reason'] is not None or abs(record['last_admitted_proper']-endpoint)>1e-12 or
                abs(accepted-endpoint)>1e-12 or len(admitted)!=len(snapshots)):
            raise ValueError('completed run does not reach its admitted requested endpoint')
    elif not record.get('stop_reason'):
        raise ValueError('terminal incomplete run requires a stated stop reason')
    authentication=dict(run_record=relative(path),run_record_sha256=sha(path),producer_commit=commit,
        source_sha256=record['source_sha256'],historical_Git_source_verified=True,
        current_package_matches=not package_differences,current_physical_package_matches=not package_differences,
        current_package_differences=package_differences,
        current_package_check_role='Informational for a historical data-only projection; not permission to resume a changed RHS.',
        current_runner_matches_historical=current_runner,
        sparse_cache=dict(path=record['state_cache'],sha256=record['state_cache_sha256']),
        initial_record=dict(path=record['initial_record'],sha256=record['initial_record_sha256'],raw_cache=initial_record['raw_cache']))
    if 'resume' in record:
        resume=record['resume'];parent_path=ROOT/resume['parent_run']
        if sha(parent_path)!=resume['parent_record_sha256']:raise ValueError('resume parent-record SHA mismatch')
        parent,pauth=authenticate(parent_path,seen)
        if (parent['producer_commit']!=resume['parent_producer_commit'] or parent['status']!=resume['parent_status'] or
                parent.get('stop_reason')!=resume.get('parent_stop_reason') or parent['modes']!=record['modes'] or
                parent['initial_record_sha256']!=record['initial_record_sha256']):
            raise ValueError('resume parent identity/preparation mismatch')
        parent_cache=ROOT/resume['state_cache']
        if (parent_cache.resolve()!=(ROOT/parent['state_cache']).resolve() or
                resume['state_cache_sha256']!=parent['state_cache_sha256'] or sha(parent_cache)!=resume['state_cache_sha256']):
            raise ValueError('resume parent raw-cache SHA mismatch')
        index=int(resume['snapshot_index'])
        with np.load(parent_cache,allow_pickle=False) as raw:
            parent_state=np.asarray(raw['states'][index]);parent_time=float(raw['times'][index])
        if hashlib.sha256(parent_state.tobytes()).hexdigest()!=resume['state_array_sha256']:
            raise ValueError('selected resume state-array SHA mismatch')
        if (not np.array_equal(parent_state,states[0]) or
                hashlib.sha256(states[0].tobytes()).hexdigest()!=resume['state_array_sha256']):
            raise ValueError('resumed child reset or changed its Cauchy state')
        if abs(parent_time-resume['handoff_proper'])>1e-12 or abs(times[0]-parent_time)>1e-12:
            raise ValueError('resume handoff time mismatch')
        ps=parent['snapshots'][index]
        if not ps['sampled_admitted']:raise ValueError('resume source snapshot is not admitted')
        names=('fixed_receiving_cap','fixed_donor_cap')
        integrals=[ps[names[0]]['integrated_material_received'],ps[names[0]]['integrated_pressure_work'],ps[names[1]]['integrated_pressure_work']]
        child_integrals=[snapshots[0][names[0]]['integrated_material_received'],snapshots[0][names[0]]['integrated_pressure_work'],snapshots[0][names[1]]['integrated_pressure_work']]
        if not np.allclose(integrals,resume['cumulative_integrals'],rtol=0,atol=2e-12) or not np.allclose(integrals,child_integrals,rtol=0,atol=2e-12):
            raise ValueError('resume reset cumulative material/pressure-work ledgers')
        parent_first=parent['snapshots'][0]
        baselines=parent.get('cap_energy_baselines',[parent_first[name]['energy']-parent_first[name]['integrated_material_received']-
                                                   parent_first[name]['integrated_pressure_work'] for name in names])
        if (not np.allclose(baselines,resume['cap_energy_baselines'],rtol=0,atol=3e-11) or
                not np.allclose(baselines,record['cap_energy_baselines'],rtol=0,atol=3e-11)):
            raise ValueError('resume changed the original fixed-cap energy baselines')
        authentication['resume']=dict(**resume,parent_archived_Git_sources_verified=True,
            parent_current_package_matches=pauth['current_package_matches'],
            parent_current_package_differences=pauth['current_package_differences'],
            selected_parent_state_SHA_verified=True,child_initial_state_bitwise_equal=True,
            cumulative_ledgers_preserved=True,original_cap_baselines_preserved=True)
    return record,authentication


def run_projection(record,authentication):
    snaps=admitted_snapshots(record)
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
                    status=record['status'],stop_reason=record['stop_reason'],
                    actual_last_accepted_proper=record['actual_last_accepted_proper'],
                    last_cached_snapshot_proper=record['snapshots'][-1]['proper_time'],
                    accepted_tail_not_cached=record['actual_last_accepted_proper']>record['snapshots'][-1]['proper_time']+1e-12,
                    all_recorded_snapshots_sampled_admitted=len(snaps)==len(record['snapshots']),
                    admission='Finite spatial samples at the recorded times, ADM-term2% scale; not a continuum error bound.'),
        temporal_tolerances=record['temporal_tolerances'],reflection_projection=record['reflection_projection'],
        controls=controls,initial=snaps[0],endpoint=snaps[-1],
        window_increments=dict(
            receiving_cap_energy=snaps[-1]['fixed_receiving_cap']['energy']-snaps[0]['fixed_receiving_cap']['energy'],
            receiving_cap_volume=snaps[-1]['fixed_receiving_cap']['volume']-snaps[0]['fixed_receiving_cap']['volume'],
            net_material_received=snaps[-1]['fixed_receiving_cap']['integrated_material_received']-snaps[0]['fixed_receiving_cap']['integrated_material_received'],
            receiving_cap_pressure_work=snaps[-1]['fixed_receiving_cap']['integrated_pressure_work']-snaps[0]['fixed_receiving_cap']['integrated_pressure_work']),
        compute=dict(threads=record['threads'],elapsed_seconds=record['elapsed_seconds'],rhs_calls=record['rhs_calls']))


def minimum_radius(snapshot):return None if snapshot['minimum'] is None else snapshot['minimum']['radius']


def band_edges(snapshot):
    b=snapshot['trapped_band']
    return None if b is None else [b['u_left'],b['u_right']]


def project(coarse_path,fine_path,output_path=None):
    start=time.perf_counter();producer_before=sha(__file__)
    coarse,cauth=authenticate(coarse_path);fine,fauth=authenticate(fine_path)
    if coarse['initial_record_sha256']!=fine['initial_record_sha256']:
        raise ValueError('runs do not start from the same directional Cauchy record')
    if coarse['source_sha256']!=fine['source_sha256']:
        raise ValueError('the two finite runs do not share the same physical/runner sources')
    a=admitted_snapshots(coarse);b=admitted_snapshots(fine)
    matches=[(x,next((y for y in b if abs(x['proper_time']-y['proper_time'])<1e-12),None)) for x in a]
    matches=[(x,y) for x,y in matches if y is not None]
    if not matches:raise ValueError('terminal runs have no common admitted recorded time')
    table=[];maximum_spreads=dict(minimum_radius=0.,fixed_cap_energy=0.,fixed_cap_volume=0.,net_material_received=0.,band_left=0.,band_right=0.)
    for x,y in matches:
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
    endpoint=matches[-1][1];cap=endpoint['fixed_receiving_cap'];window=[matches[0][0]['proper_time'],matches[-1][0]['proper_time']]
    target_completed=coarse['status']=='endpoint_reached' and fine['status']=='endpoint_reached'
    output_path=ROOT/'results/summaries/nsc-atlas-valve-directional-evolution-v1.json' if output_path is None else output_path
    result=dict(schema='nsc-atlas-valve-directional-evolution-v1',
        projection_producer=dict(path=relative(Path(__file__)),sha256=producer_before,
                                 role='Passive projection of authenticated historical terminal trajectories; no current runner imported.'),
        runs=[run_projection(coarse,cauth),run_projection(fine,fauth)],
        comparison=dict(run_order=[coarse['modes'],fine['modes']],common_admitted_recorded_domain=window,
                        matched_times=dict(columns=columns,rows=compact_rows),
                        maximum_observed_absolute_spreads=maximum_spreads,
                        scope='Empirical differences between these two finite runs, not error bars or a continuum/exact-Cauchy-data certificate.'),
        interpretation=dict(net_material='Signed accumulated material supply through the fixedu=0 interface; not gross crossings or energy creation.',
                            energy_ledger='Fixed-cap normal matter-energy change=net material supplied+integrated pressure work; no moving-boundary sweep is mixed into supply.',
                            cumulative_origin='Resumed material/work values remain cumulative from originaltau0; window increments subtract the authenticated handoff values.',
                            migrating_minimum='An area minimum on each slice; its coordinate motion is not a material observer worldline.',
                            trapped_band='Recorded marginal roots delimit a sampled future-trapped band; no new horizon classification is performed here.'),
        verdict=dict(defined_directional_preparation=True,inherited_selection_proved=False,
                     requested_finite_endpoint_completed=target_completed,regular_receiving_child_proved=False,
                     endpoint_receiving_cap_expanding=cap['bulk_volume_rate']>0,
                     endpoint_receiving_cap_volume_rate=cap['bulk_volume_rate'],endpoint_receiving_cap_w=cap['w'],
                     new_horizon_classification_or_NEC_gate_computed=False,
                     scope=f'Common admitted recorded window{window[0]:g}–{window[1]:g} only; requested endpoints and stop causes remain in per-run metadata. No physical singularity, child or cosmological claim.'),
        projection_elapsed_seconds=time.perf_counter()-start,
        reproduction='python scripts/run_python.py scripts/summarize_nsc_atlas_valve_directional.py '
            f'--coarse-record {relative(coarse_path)} --fine-record {relative(fine_path)} '
            f'--output {relative(output_path)}')
    if sha(__file__)!=producer_before:raise RuntimeError('projection script changed while reading the historical records')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coarse-record',type=Path,default=ROOT/'.build/valve-directional-evolution-N64-v1/run.json')
    parser.add_argument('--fine-record',type=Path,default=ROOT/'.build/valve-directional-evolution-N96-v1/run.json')
    parser.add_argument('--output',type=Path,default=ROOT/'results/summaries/nsc-atlas-valve-directional-evolution-v1.json')
    args=parser.parse_args();result=project(args.coarse_record,args.fine_record,args.output)
    encoded=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if len(encoded.encode('utf-8'))>50*1024:raise RuntimeError('compact projection exceeds50KiB; no fields silently dropped')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(encoded,encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(args.output),bytes=len(encoded.encode('utf-8')),
                         matched_times=len(result['comparison']['matched_times']['rows']),
                         max_spreads=result['comparison']['maximum_observed_absolute_spreads'],verdict=result['verdict']),indent=2))


if __name__=='__main__':main()
