"""Preserve small bosonic response records and compare like controls.

Full transport matrices remain in regenerable local NPZ caches. Optional
Schur extraction uses the uncoupled eliminated-field propagator, retaining
the actual coupled vertices on the same background.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from record_nsc_atlas_gh_campaign import source_binding

ROOT=Path(__file__).resolve().parents[1]
READOUTS=(
    'local_log_g_KK_per_unit_momentum_kick',
    'local_volume_field_per_unit_momentum_kick',
    'local_scalar_field_per_unit_momentum_kick',
    'volume_reference_response_norm','scalar_reference_response_norm',
    'volume_to_scalar_reference_response_ratio',
)
ONSET=('volume_q_cubic_coefficient_norm','local_volume_cubic_coefficient',
    'local_log_g_KK_cubic_coefficient')


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _responses(record):
    answers=record['checkpoints'][-1]['response']
    result={int(row['source_harmonic']):row for row in answers}
    if set(result)!={1,2} or len(answers)!=2:
        raise ValueError('two distinct inherited source harmonics are required')
    return result


def response_difference(before,after):
    """Final-minus-initial comparison, without asserting a truncation bound."""
    old,new=_responses(before),_responses(after)
    changes=[]
    for source in (1,2):
        signed={};relative={}
        for name in READOUTS:
            a,b=old[source][name],new[source][name]
            value=None if a is None or b is None else float(b-a)
            signed[name]=value
            relative[name]=None if value is None or a==0 else value/abs(a)
        changes.append(dict(source_harmonic=source,signed_change=signed,
            absolute_change={name:None if value is None else abs(value) for name,value in signed.items()},
            change_over_previous_absolute_value=relative))
    onset=None
    if before.get('onset_cubic_response') and after.get('onset_cubic_response'):
        old_onset={row['source_harmonic']:row for row in before['onset_cubic_response']}
        new_onset={row['source_harmonic']:row for row in after['onset_cubic_response']}
        if set(old_onset)==set(new_onset)=={1,2}:
            onset=[dict(source_harmonic=source,
                signed_change={name:new_onset[source][name]-old_onset[source][name] for name in ONSET},
                absolute_change={name:abs(new_onset[source][name]-old_onset[source][name]) for name in ONSET})
                for source in (1,2)]
    return dict(from_probe_degree=before['probe_degree'],to_probe_degree=after['probe_degree'],
        from_max_step=before['max_step'],to_max_step=after['max_step'],responses=changes,
        onset_cubic_changes=onset)


def aggregate(records,*,bindings,hashes=None,names=None):
    """Keep every supplied small record; compare verified reached endpoints."""
    if not records or len(bindings)!=len(records):
        raise ValueError('records and one explicit source-binding result each are required')
    hashes=[None]*len(records) if hashes is None else hashes
    names=[None]*len(records) if names is None else names
    if len(hashes)!=len(records) or len(names)!=len(records):
        raise ValueError('one input name and hash per record are required')
    family=('background_record_sha256','background_history_sha256','background_producer_commit',
        'background_degree','background_epsilon','background_temporal_tolerances')
    if len({json.dumps({key:row[key] for key in family},sort_keys=True) for row in records})!=1:
        raise ValueError('all controls must share one background history and producer record')
    keys=[(row['probe_degree'],row['max_step'],row['decoupled_control']) for row in records]
    if len(set(keys))!=len(keys):raise ValueError('duplicate probe degree, step and coupling control')
    cases=[];eligible=[]
    for row,bound,sha,name in zip(records,bindings,hashes,names):
        reached=row['status']=='endpoint_reached'
        disposition='unverified_source_binding' if not bound else ('finite_endpoint_control' if reached else 'failed_attempt')
        cases.append(dict(input=name,input_sha256=sha,source_binding_verified=bool(bound),
            disposition=disposition,physical_prediction_admitted=False,run=row))
        if bound and reached:
            _responses(row)
            eligible.append(row)
    if eligible:
        producers={json.dumps(row['source_sha256'],sort_keys=True) for row in eligible}
        if len(producers)!=1:raise ValueError('successful controls use different producer source hashes')
        endpoints=[row['checkpoints'][-1] for row in eligible]
        first=endpoints[0]
        for point in endpoints[1:]:
            for key in ('coordinate_t','observer_proper_duration'):
                if not math.isclose(point[key],first[key],rel_tol=0,abs_tol=1e-9):
                    raise ValueError('successful controls do not share a physical endpoint')
    coupled=[row for row in eligible if not row['decoupled_control']]
    spatial=[];temporal=[];controls=[]
    for step in sorted({row['max_step'] for row in coupled}):
        group=sorted((row for row in coupled if row['max_step']==step),key=lambda row:row['probe_degree'])
        spatial.extend(response_difference(a,b) for a,b in zip(group,group[1:]))
    for degree in sorted({row['probe_degree'] for row in coupled}):
        group=sorted((row for row in coupled if row['probe_degree']==degree),key=lambda row:-row['max_step'])
        temporal.extend(response_difference(a,b) for a,b in zip(group,group[1:]))
    for row in eligible:
        if not row['decoupled_control']:continue
        matching=[other for other in coupled if other['probe_degree']==row['probe_degree']
            and other['max_step']==row['max_step']]
        if matching:controls.append(response_difference(row,matching[0]))
    diagnostics=[dict(probe_degree=row['probe_degree'],max_step=row['max_step'],
        decoupled_control=row['decoupled_control'],
        max_checkpoint_symplectic_defect=row['max_checkpoint_symplectic_defect'],
        final_symplectic=row['checkpoints'][-1]['symplectic'],
        onset_cubic_response=row.get('onset_cubic_response')) for row in eligible]
    return dict(schema='nsc-atlas-boson-summary-v1',
        background={key:records[0][key] for key in family},cases=cases,
        all_source_bindings_verified=all(bindings),spatial_changes=spatial,
        temporal_changes=temporal,decoupled_to_coupled_changes=controls,
        canonical_and_onset_controls=diagnostics,
        physical_prediction_admitted=False,inherited_covariance_selected=False,
        entropy_or_noise_selected=False,
        qualification='Adjacent finite-band and time changes are empirical controls, not certified error bounds. Decoupled operators are diagnostic controls on the same background.',
        normalization='Unit canonical momentum-density kicks in reference-orthonormal inherited harmonics. Vacuum-mass-weighted reference norms omit gradient, metric and interaction energy and are not conversion fractions.',
        scope='State-independent finite Gaussian matter response on the evolving classical parent. Background geodesic held fixed; no metric fluctuations, quantum backreaction, continuum microcausality, selected covariance, entropy, noise or SI scale.',
        source_sha256={'scripts/summarize_nsc_atlas_boson_response.py':digest(__file__)})


def phase_partition(fields):
    fields=np.asarray(fields)
    if fields.ndim!=1 or fields.dtype.kind not in 'iu' or len(fields)==0 or not np.isin(fields,[0,1]).all():
        raise ValueError('field labels must contain canonical S=0 and chi=1 modes')
    n=len(fields);chi=np.flatnonzero(fields==1);volume=np.flatnonzero(fields==0)
    if not len(chi) or not len(volume):raise ValueError('both fields must be retained')
    return np.r_[chi,chi+n],np.r_[volume,volume+n],chi,volume


def memory_matrix(generator_t,generator_s,control_t,control_s,fields):
    """D_xy(t) U_y(t,s) D_yx(s), using an actual decoupled y inverse."""
    x,y,_,_=phase_partition(fields)
    matrices=[np.asarray(value,dtype=float) for value in
        (generator_t,generator_s,control_t,control_s)]
    size=2*len(fields)
    if any(value.shape!=(size,size) or not np.all(np.isfinite(value)) for value in matrices):
        raise ValueError('finite matching canonical matrices are required')
    Dt,Ds,Mt,Ms=matrices
    old=Ms[np.ix_(y,y)];late=Mt[np.ix_(y,y)]
    try:Uy=np.linalg.solve(old.T,late.T).T
    except np.linalg.LinAlgError as error:raise ValueError('decoupled eliminated propagator is singular') from error
    return Dt[np.ix_(x,y)]@Uy@Ds[np.ix_(y,x)]


def _load_cache(path,record):
    if digest(path)!=record.get('matrix_cache_sha256'):
        raise ValueError('transport matrix cache does not match its record hash')
    with np.load(path,allow_pickle=False) as data:
        times=data['times'].copy();matrices=data['matrices'].copy();fields=data['fields'].copy()
        metadata=json.loads(str(data['metadata']))
    if (metadata['probe_degree']!=record['probe_degree']
            or metadata['history_sha256']!=record['background_history_sha256']
            or metadata['source_sha256']!=record['source_sha256']):
        raise ValueError('transport matrix cache metadata does not match its producer')
    expected=np.array([point['coordinate_t'] for point in record['checkpoints']])
    size=2*len(fields)
    if times.shape!=expected.shape or not np.allclose(times,expected,rtol=0,atol=1e-12):
        raise ValueError('transport cache checkpoint times do not match its record')
    if matrices.shape!=(len(times),size,size) or not np.all(np.isfinite(matrices)):
        raise ValueError('invalid transport matrix cache dimensions or values')
    phase_partition(fields)
    return times,matrices,fields


def schur_diagnostic(history_path,coupled_record,coupled_cache,control_record,control_cache):
    """Small chart-time memory readout, with all complete matrices kept local."""
    from threadpoolctl import threadpool_limits
    from recursive_horizons.nsc_atlas_boson_parent import ParentBosonProbe
    from recursive_horizons.nsc_atlas_boson_transport import hamiltonian_generator
    from recursive_horizons.nsc_atlas_gh_history import DenseGHHistory
    if any(row['status']!='endpoint_reached' or not source_binding(row) for row in (coupled_record,control_record)):
        raise ValueError('Schur extraction requires two verified reached endpoints')
    if coupled_record['decoupled_control'] or not control_record['decoupled_control']:
        raise ValueError('the coupled and decoupled records have different roles')
    if coupled_record['probe_degree']!=control_record['probe_degree']:
        raise ValueError('Schur extraction requires the same probe degree')
    for key in ('background_history_sha256','background_record_sha256','source_sha256'):
        if coupled_record[key]!=control_record[key]:raise ValueError('Schur records do not share source and background')
    history_sha=digest(history_path)
    if history_sha!=coupled_record['background_history_sha256']:
        raise ValueError('Schur history does not match the response records')
    for name,expected in coupled_record['source_sha256'].items():
        if name.startswith('src/') and digest(ROOT/name)!=expected:
            raise ValueError('current Schur model differs from the response producer: '+name)
    times,maps,fields=_load_cache(coupled_cache,coupled_record)
    control_times,control_maps,control_fields=_load_cache(control_cache,control_record)
    if not np.array_equal(fields,control_fields) or not np.array_equal(times,control_times):
        raise ValueError('coupled and control matrix caches have different canonical labels or times')
    history=DenseGHHistory.load(history_path)
    if history.metadata['producer_commit']!=coupled_record['background_producer_commit']:
        raise ValueError('history producer differs from the response record')
    for name,expected in history.metadata['source_sha256'].items():
        if name.startswith('src/') and coupled_record['source_sha256'].get(name)!=expected:
            raise ValueError('history and response model source hashes differ')
    with threadpool_limits(limits=4):
        probe=ParentBosonProbe(history,coupled_record['probe_degree'])
        if not np.array_equal(fields,probe.fields) or not np.allclose(times,probe.checkpoint_times,rtol=0,atol=1e-12):
            raise ValueError('probe labels or physical checkpoints differ from matrix cache')
        generators=[];samples=[]
        x,y,chi,_=phase_partition(fields)
        for t in times:
            details=probe.hamiltonian_details(float(t));D=hamiltonian_generator(details['G'])
            generators.append(D)
            samples.append(dict(coordinate_t=float(t),
                minimum_nodal_mass_eigenvalue=details['minimum_nodal_mass_eigenvalue'],
                retained_to_eliminated_generator_frobenius=float(np.linalg.norm(D[np.ix_(y,x)])),
                eliminated_to_retained_generator_frobenius=float(np.linalg.norm(D[np.ix_(x,y)]))))
        proper=np.array([point['observer_proper_duration'] for point in coupled_record['checkpoints']])
        midpoint=(proper[0]+proper[-1])/2
        matching=np.flatnonzero(np.isclose(proper,midpoint,rtol=0,atol=1e-9))
        if len(times)!=5 or len(matching)!=1:
            raise ValueError('five checkpoints including the proper midpoint are required')
        it=len(times)-1;is_=int(matching[0])
        kernel=memory_matrix(generators[it],generators[is_],control_maps[it],control_maps[is_],fields)
        pq=kernel[len(chi):,:len(chi)]
        vectors=np.stack([probe.source_kick(source)[probe.n:][chi] for source in (1,2)],axis=1)
        mu=math.sqrt(probe.system.parent.mass_squared)
        projection=vectors.T@pq@vectors/mu**3
    return dict(retained_field='chi',eliminated_field='S',probe_degree=probe.probe_degree,
        final_coordinate_t=float(times[it]),source_coordinate_s=float(times[is_]),
        final_proper_duration=float(proper[it]),source_proper_duration=float(proper[is_]),
        normalization_mu_star=mu,normalization_power=3,
        source_harmonics=[1,2],projected_pchi_from_qchi_memory_over_mu_cubed=projection.tolist(),
        kernel_frobenius=float(np.linalg.norm(kernel)),sampled_vertices=samples,
        control_max_step=control_record['max_step'],coupled_max_step=coupled_record['max_step'],
        control_cache_sha256=digest(control_cache),coupled_cache_sha256=digest(coupled_cache),
        background_history_sha256=history_sha,
        rule='K(t,s)=D_xy(t) U_y(t,s) D_yx(s), with U_y from the decoupled control, not a subblock of the coupled evolution.',
        scope='Chart-time canonical memory diagnostic for the finite inherited harmonic sector. It is not vacuum stress, Lambda, a physical energy fraction or a selected covariance/noise/entropy.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,action='append',required=True)
    parser.add_argument('--output',type=Path,required=True)
    for flag in ('history','coupled-record','coupled-matrix','control-record','control-matrix'):
        parser.add_argument('--'+flag,type=Path)
    args=parser.parse_args();output=(ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('fresh repository output required')
    paths=[(ROOT/name).resolve() for name in args.input]
    if any(not path.is_relative_to(ROOT) for path in paths):raise ValueError('repository response records required')
    records=[json.loads(path.read_bytes()) for path in paths]
    result=aggregate(records,bindings=[source_binding(row) for row in records],
        hashes=[digest(path) for path in paths],names=[path.relative_to(ROOT).as_posix() for path in paths])
    extra=[args.history,args.coupled_record,args.coupled_matrix,args.control_record,args.control_matrix]
    if any(value is not None for value in extra):
        if not all(value is not None for value in extra):raise ValueError('all five Schur extraction arguments are required')
        h,cr,cm,ur,um=[(ROOT/value).resolve() for value in extra]
        if any(not path.is_relative_to(ROOT) for path in (h,cr,cm,ur,um)):
            raise ValueError('repository Schur inputs required')
        if any(not path.is_relative_to(ROOT/'.build') for path in (h,cm,um)):
            raise ValueError('regenerable history and matrix caches must remain under .build')
        result['schur_memory']=schur_diagnostic(h,json.loads(cr.read_bytes()),cm,json.loads(ur.read_bytes()),um)
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(cases=len(records),verified=sum(case['source_binding_verified'] for case in result['cases']),
        spatial_comparisons=len(result['spatial_changes']),time_comparisons=len(result['temporal_changes']),
        schur_memory='schur_memory' in result,physical_prediction_admitted=False)))


if __name__=='__main__':main()
