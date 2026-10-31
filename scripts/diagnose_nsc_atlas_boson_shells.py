"""Diagnose retained harmonic tails from verified finite bosonic maps.

Only small-basis label/source reconstruction is performed. No background
initialization, trajectory evaluation or propagation is needed. The norm
is sum_l(m_f^2*q_l^2+p_l^2), using the two fixed vacuum masses. It is a
reference phase-space norm, not physical energy or a conversion fraction.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_boson_hamiltonian import potential_hessian
from recursive_horizons.nsc_atlas_parent_scale_selector import VolumeParent
from recursive_horizons.nsc_atlas_s3_hopf import build_hopf_frame_basis
from record_nsc_atlas_gh_campaign import source_binding

ROOT=Path(__file__).resolve().parents[1]
ACTION_SOURCES=(
    'src/recursive_horizons/nsc_atlas_boson_hamiltonian.py',
    'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
    'src/recursive_horizons/nsc_atlas_s3_hopf.py',
    'src/recursive_horizons/nsc_atlas_s3_frame.py',
)


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reference_action():
    """Derive the two reference masses from the current fixed action."""
    minima=[point for point in VolumeParent().vacuum_roots() if point['kind']=='minimum']
    if len(minima)!=1:raise ValueError('the fixed action must have one retained minimum')
    star=float(minima[0]['s'])
    masses=np.diag(potential_hessian(np.zeros(2),sstar=star)).copy()
    if not np.all(np.isfinite(masses)) or np.any(masses<=0):
        raise ValueError('positive finite reference vacuum masses are required')
    return dict(sstar=star,masses_squared=masses,
        source_sha256={name:digest(ROOT/name) for name in ACTION_SOURCES})


def reconstruct_sector(degree):
    """Reconstruct S(++),chi(--) labels and the two normalized source kicks."""
    if isinstance(degree,(bool,np.bool_)) or not isinstance(degree,(int,np.integer)) or degree<2:
        raise ValueError('probe degree must be an integer >=2')
    basis=build_hopf_frame_basis(int(degree))
    si=np.flatnonzero(np.all(basis.characters==[1,1],axis=1))
    ci=np.flatnonzero(np.all(basis.characters==[-1,-1],axis=1))
    selected=np.r_[si,ci]
    labels=tuple(basis.labels[index] for index in selected)
    fields=np.r_[np.zeros(len(si),dtype=int),np.ones(len(ci),dtype=int)]
    harmonic_degree=np.array([label[0] for label in labels],dtype=int)
    kicks=[]
    profiles=(math.sqrt(2)/math.pi*basis.points[:,0],
        math.sqrt(12)/math.pi*basis.points[:,1]*basis.points[:,3])
    for profile in profiles:
        coefficients=np.zeros(len(selected))
        coefficients[len(si):]=basis.project(profile)[ci]
        norm=float(np.linalg.norm(coefficients))
        if not math.isfinite(norm) or norm<=0:
            raise ValueError('the inherited source profile is absent')
        coefficients/=norm
        kicks.append(np.r_[np.zeros(len(selected)),coefficients])
    return dict(probe_degree=int(degree),labels=labels,fields=fields,
        degrees=harmonic_degree,kicks=np.asarray(kicks),
        mode_counts=dict(S=len(si),chi=len(ci),total=len(selected)))


def _states(value,sector):
    value=np.asarray(value)
    n=len(sector['fields'])
    if value.dtype.kind not in 'iuf' or value.shape!=(2,2*n) or not np.all(np.isfinite(value)):
        raise ValueError('two finite q-then-p source states are required')
    return np.asarray(value,dtype=float)


def shell_diagnosis(states,sector,masses):
    """Report squared-norm power and square-root norm fractions separately."""
    states=_states(states,sector);masses=np.asarray(masses,dtype=float)
    if masses.shape!=(2,) or not np.all(np.isfinite(masses)) or np.any(masses<=0):
        raise ValueError('two positive finite mass-squared weights are required')
    fields=sector['fields'];degree=sector['degrees'];P=sector['probe_degree'];n=len(fields)
    answer=[]
    for source,state in enumerate(states,1):
        for field,name in ((0,'S'),(1,'chi')):
            mask=fields==field
            power=masses[field]*state[:n]**2+state[n:]**2
            by_degree={str(level):float(power[mask&(degree==level)].sum()) for level in range(P+1)}
            total=float(sum(by_degree.values()))
            one=by_degree[str(P)]
            two=sum(value for level,value in by_degree.items() if int(level)>=P-1)
            fraction=lambda value:value/total if total>0 else None
            sqrt_fraction=lambda value:math.sqrt(value/total) if total>0 else None
            answer.append(dict(source_harmonic=source,field=name,
                reference_norm=math.sqrt(total),reference_norm_squared=total,
                power_by_degree=by_degree,
                power_fraction_by_degree={level:fraction(value) for level,value in by_degree.items()},
                last_degree_power_fraction=fraction(one),last_two_degrees_power_fraction=fraction(two),
                last_degree_norm_fraction=sqrt_fraction(one),last_two_degrees_norm_fraction=sqrt_fraction(two)))
    return answer


def relative_state_gaps(before,after,sector,masses):
    """Whole-state temporal gap, divided by the tighter state's field norm."""
    before=_states(before,sector);after=_states(after,sector);n=len(sector['fields'])
    masses=np.asarray(masses,dtype=float);result=[]
    for source,(old,new) in enumerate(zip(before,after),1):
        difference=new-old
        for field,name in ((0,'S'),(1,'chi')):
            idx=sector['fields']==field
            squared=lambda z:float(masses[field]*(z[:n][idx]@z[:n][idx])+z[n:][idx]@z[n:][idx])
            denominator=math.sqrt(squared(new));gap=math.sqrt(squared(difference))
            result.append(dict(source_harmonic=source,field=name,absolute_reference_norm_gap=gap,
                gap_over_tighter_reference_norm=gap/denominator if denominator>0 else None))
    return result


def load_verified_case(record_path,matrix_path,action):
    """Check Git/source binding, cache SHA, metadata, shape and endpoint."""
    record_path=Path(record_path);matrix_path=Path(matrix_path)
    record=json.loads(record_path.read_bytes())
    if record.get('status')!='endpoint_reached' or not source_binding(record):
        raise ValueError('a reached, source-bound producer record is required')
    for name,expected in action['source_sha256'].items():
        if record.get('source_sha256',{}).get(name)!=expected:
            raise ValueError('current action or basis differs from producer: '+name)
    matrix_sha=digest(matrix_path)
    if record.get('matrix_cache_sha256')!=matrix_sha:
        raise ValueError('matrix SHA does not match the producer record')
    sector=reconstruct_sector(record['probe_degree']);n=len(sector['fields'])
    with np.load(matrix_path,allow_pickle=False) as data:
        times=data['times'];matrices=data['matrices'];fields=data['fields']
        metadata=json.loads(str(data['metadata']))
    if metadata.get('probe_degree')!=record['probe_degree'] or metadata.get('history_sha256')!=record['background_history_sha256']:
        raise ValueError('matrix metadata has a different probe or background')
    if metadata.get('source_sha256')!=record['source_sha256']:
        raise ValueError('matrix and record source metadata differ')
    checkpoints=record['checkpoints']
    if times.shape!=(len(checkpoints),) or matrices.shape!=(len(times),2*n,2*n):
        raise ValueError('matrix/checkpoint shapes differ from the reconstructed sector')
    if not np.all(np.isfinite(times)) or not np.all(np.isfinite(matrices)) or np.any(np.diff(times)<=0):
        raise ValueError('matrix cache needs finite maps at increasing times')
    if not np.array_equal(fields,sector['fields']):
        raise ValueError('matrix field order differs from the reconstructed sector')
    if not np.allclose(times,[point['coordinate_t'] for point in checkpoints],rtol=0,atol=1e-12):
        raise ValueError('matrix cache and record times differ')
    states=(matrices[-1]@sector['kicks'].T).T
    summary=dict(record=str(record_path.relative_to(ROOT)),matrix=str(matrix_path.relative_to(ROOT)),
        record_sha256=digest(record_path),matrix_sha256=matrix_sha,
        source_binding_verified=True,matrix_binding_verified=True,
        probe_degree=record['probe_degree'],max_step=record['max_step'],
        decoupled_control=record['decoupled_control'],mode_counts=sector['mode_counts'],
        coordinate_endpoint=float(times[-1]),
        observer_proper_endpoint=checkpoints[-1]['observer_proper_duration'],
        responses=shell_diagnosis(states,sector,action['masses_squared']))
    return dict(record=record,sector=sector,states=states,summary=summary)


def aggregate(cases,action):
    if not cases:raise ValueError('at least one verified case is required')
    family=('background_record_sha256','background_history_sha256','background_producer_commit',
        'background_degree','background_epsilon','background_temporal_tolerances','source_sha256')
    records=[case['record'] for case in cases]
    if len({json.dumps({key:row[key] for key in family},sort_keys=True) for row in records})!=1:
        raise ValueError('diagnosed controls must share one background and producer sources')
    keys=[(row['probe_degree'],row['max_step'],row['decoupled_control']) for row in records]
    if len(set(keys))!=len(keys):raise ValueError('duplicate probe/time/control case')
    endpoint=cases[0]['summary']
    for case in cases[1:]:
        for key in ('coordinate_endpoint','observer_proper_endpoint'):
            if not math.isclose(case['summary'][key],endpoint[key],rel_tol=0,abs_tol=1e-9):
                raise ValueError('controls must share one physical endpoint')
    comparisons=[]
    groups={(row['probe_degree'],row['decoupled_control']) for row in records}
    for degree,control in sorted(groups):
        group=sorted((case for case in cases if case['record']['probe_degree']==degree
            and case['record']['decoupled_control']==control),key=lambda case:-case['record']['max_step'])
        for before,after in zip(group,group[1:]):
            comparisons.append(dict(probe_degree=degree,decoupled_control=control,
                from_max_step=before['record']['max_step'],to_max_step=after['record']['max_step'],
                responses=relative_state_gaps(before['states'],after['states'],after['sector'],action['masses_squared'])))
    return dict(schema='nsc-atlas-boson-shells-v1',
        definition='Per-field reference phase-space norm squared=sum_l(m_f^2*q_l^2+p_l^2), using fixed vacuum masses. Harmonic power fractions use this squared norm; norm fractions are their square roots. Neither is physical energy or a conversion fraction.',
        reference_action=dict(sstar=action['sstar'],masses_squared=action['masses_squared'].tolist(),
            source_sha256=action['source_sha256']),cases=[case['summary'] for case in cases],
        temporal_changes=comparisons,spatial_convergence_decided=False,
        physical_prediction_admitted=False,
        source_sha256={'scripts/diagnose_nsc_atlas_boson_shells.py':digest(Path(__file__))})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',type=Path,action='append',required=True)
    parser.add_argument('--matrix',type=Path,action='append',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if len(args.record)!=len(args.matrix):parser.error('one --matrix is required per --record, in the same order')
    output=(ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('fresh repository output required')
    pairs=[]
    for record,matrix in zip(args.record,args.matrix):
        record=(ROOT/record).resolve();matrix=(ROOT/matrix).resolve()
        if not record.is_relative_to(ROOT) or not matrix.is_relative_to(ROOT/'.build'):
            raise ValueError('repository records and local regenerable matrix caches are required')
        pairs.append((record,matrix))
    with threadpool_limits(limits=2):
        action=reference_action()
        cases=[load_verified_case(record,matrix,action) for record,matrix in pairs]
        result=aggregate(cases,action)
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(cases=len(cases),temporal_comparisons=len(result['temporal_changes']),
        all_bindings_verified=True,physical_prediction_admitted=False)))


if __name__=='__main__':main()
