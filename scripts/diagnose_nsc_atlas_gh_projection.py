"""Compare three quadratures at the SAME cached L=12 GH coefficient states.

This does not evolve or recreate initial data. It separates changes in the
projected nonlinear RHS due to integration quadrature from the fixed degree
cutoff. It is not a physical acceptance or time-evolution convergence test.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from numba import get_num_threads,set_num_threads
from scipy.fft import set_workers
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
QUADRATURES=((11,38),(15,50),(19,62))
THREADS=4
USED_CORE=(
    'src/recursive_horizons/nsc_atlas_parent_gh.py',
    'src/recursive_horizons/nsc_atlas_s3_frame.py',
    'src/recursive_horizons/nsc_atlas_s3_hopf.py',
    'src/recursive_horizons/nsc_atlas_gh_pointwise.py',
    'src/recursive_horizons/nsc_atlas_gh_constraints.py',
    'src/recursive_horizons/nsc_atlas_finite_parent_data.py',
    'src/recursive_horizons/nsc_atlas_parent_circulation.py',
    'src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
    'src/recursive_horizons/nsc_atlas_offset_generator.py',
)
GROUPS={
    'metric_acceleration':slice(10,20),
    'volume_acceleration':slice(22,23),
    'chi_acceleration':slice(23,24),
    'metric_first_derivative':slice(0,10),
    'scalar_first_derivatives':slice(20,22),
}


class InvalidFixedState(RuntimeError):
    pass


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repository_path(name):
    path=(ROOT/name).resolve()
    if not path.is_relative_to(ROOT):raise ValueError('path escapes the repository')
    return path


def difference(coarse,reference):
    gap=np.asarray(coarse)-np.asarray(reference)
    absolute=float(np.linalg.norm(gap.ravel()))
    norm=float(np.linalg.norm(np.asarray(reference).ravel()))
    return dict(absolute_coefficient_l2_gap=absolute,reference_coefficient_l2_norm=norm,
                relative_coefficient_l2_gap=absolute/norm if norm else (0. if absolute==0 else None),
                maximum_absolute_coefficient_gap=float(np.max(abs(gap))) if gap.size else 0.)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,default=Path('.build/gh-finite-L12-hopf.npz'))
    parser.add_argument('--output',type=Path,default=Path('.build/gh-L12-projection-diagnosis.json'))
    parser.add_argument('--budget-seconds',type=float,default=120.)
    args=parser.parse_args()
    if not 0<args.budget_seconds<=120:raise ValueError('require a positive budget <=120 seconds')
    cache=repository_path(args.cache);output=repository_path(args.output)
    if not cache.is_file():raise ValueError('cached coefficient states are required')
    if not output.is_relative_to(ROOT/'.build') or output.exists():raise ValueError('fresh .build output required')
    started=time.perf_counter()
    cache_hash=digest(cache)
    def check_budget():
        if time.perf_counter()-started>args.budget_seconds:raise TimeoutError('fixed-state quadrature budget exhausted')
    with np.load(cache,allow_pickle=False) as saved:
        metadata=json.loads(str(saved['metadata']))
        times=np.array(saved['times'],copy=True);states=np.array(saved['states'],copy=True)
    if metadata.get('degree')!=12 or metadata.get('gauge')!='harmonic' or metadata.get('backend')!='hopf':
        raise ValueError('this bounded diagnosis requires the harmonic-gauge Hopf L=12 cache')
    if times.ndim!=1 or states.ndim!=2 or len(times)!=len(states) or len(states)<4:
        raise ValueError('incompatible cached state arrays')
    cache_bindings=metadata.get('source_sha256',{})
    cache_matches={name:digest(repository_path(name))==expected for name,expected in cache_bindings.items()}
    if not cache_matches or not all(cache_matches.values()):raise ValueError('cached physics source bindings differ from the current checkout')
    before={name:digest(repository_path(name)) for name in USED_CORE}
    from recursive_horizons.nsc_atlas_parent_gh import ParentGH
    old_threads=get_num_threads()
    evaluations=[];coefficients={};state_records=[];status='complete';reason=None
    try:
        set_num_threads(THREADS)
        with threadpool_limits(limits=THREADS),set_workers(THREADS):
            systems=[]
            for order in QUADRATURES:
                check_budget()
                system=ParentGH(12,float(metadata['epsilon']),gauge='harmonic',backend='hopf',quadrature_order=order)
                if systems and (system.basis.labels!=systems[0].basis.labels or system.dof!=systems[0].dof
                                or not np.array_equal(system.offsets,systems[0].offsets)):
                    raise ValueError('quadrature changed the represented coefficient basis')
                systems.append(system)
            base=systems[0]
            if states.shape[1]!=base.dof+10:raise ValueError('cached packed state has incompatible observer or coefficient layout')
            proper=states[:,base.dof+8]
            wanted=[0.,float(proper[-1])/4,float(proper[-1])/2,float(proper[-1])]
            selected=[int(np.argmin(abs(proper-t))) for t in wanted]
            if len(set(selected))!=4:raise ValueError('cache does not contain the four required distinct proper-time snapshots')
            ell=np.array([label[0] for label in base.basis.labels]);high=ell>=11
            for index in selected:
                Cstate=base.unpack(states[index,:base.dof])
                power=Cstate[21]**2+Cstate[23]**2/base.parent.mass_squared
                denominator=float(sum(power))
                state_records.append(dict(cache_index=index,coordinate_t=float(times[index]),
                    observer_proper_duration=float(proper[index]),
                    packed_state_sha256=hashlib.sha256(states[index].tobytes()).hexdigest(),
                    represented_chi_last_two_shell_power_fraction=float(sum(power[high])/denominator) if denominator else 0.))
                for order,system in zip(QUADRATURES,systems):
                    check_budget();tick=time.perf_counter()
                    state=states[index].copy()
                    derivative=system.rhs(float(times[index]),state)
                    if not np.all(np.isfinite(derivative)):
                        raise InvalidFixedState(f'cache state {index} has invalid RHS on quadrature {order}')
                    if not np.array_equal(state,states[index]):raise ArithmeticError('RHS evaluation changed the cached coefficient state')
                    C=system.unpack(derivative[:system.dof]);coefficients[(index,order)]=C
                    evaluations.append(dict(cache_index=index,quadrature=list(order),grid_nodes=system.N,
                        elapsed_seconds=time.perf_counter()-tick,
                        first_derivative_self_consistency={
                            'metric':difference(C[:10],Cstate[10:20]),
                            'scalars':difference(C[20:22],Cstate[22:24])}))
                    check_budget()
                    print(f'fixed L12 state{index} proper={proper[index]:.8g}, quadrature={order}, nodes={system.N}',flush=True)
    except TimeoutError as error:
        status='budget_exhausted';reason=str(error)
    except InvalidFixedState as error:
        status='fixed_state_invalid';reason=str(error)
    finally:
        set_num_threads(old_threads)
    after={name:digest(repository_path(name)) for name in USED_CORE}
    if before!=after:raise RuntimeError('physics source changed during the fixed-state diagnosis')
    if digest(cache)!=cache_hash:raise RuntimeError('cache input changed during diagnosis')
    comparisons=[]
    for state_record in state_records:
        index=state_record['cache_index']
        for left,right in ((QUADRATURES[0],QUADRATURES[1]),(QUADRATURES[1],QUADRATURES[2]),(QUADRATURES[0],QUADRATURES[2])):
            if (index,left) not in coefficients or (index,right) not in coefficients:continue
            comparisons.append(dict(cache_index=index,coarse_quadrature=list(left),reference_quadrature=list(right),
                groups={name:difference(coefficients[(index,left)][rows],coefficients[(index,right)][rows]) for name,rows in GROUPS.items()}))
    record=dict(schema='nsc-atlas-gh-fixed-L12-quadrature-diagnosis-v1',status=status,stop_reason=reason,
        source_context_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        cache_input=dict(path=cache.relative_to(ROOT).as_posix(),sha256=cache_hash,metadata=metadata),
        degree=12,epsilon=float(metadata['epsilon']),gauge='harmonic',backend='hopf',threads=THREADS,
        quadratures=[list(q) for q in QUADRATURES],states=state_records,evaluations=evaluations,comparisons=comparisons,
        cached_source_bindings_match_current=cache_matches,used_core_sha256=before,
        producer_sha256=digest(Path(__file__).resolve()),elapsed_seconds=time.perf_counter()-started,
        coefficient_norm='Euclidean L2 of real normalized harmonic component coefficients; not a physical tensor norm',
        scope='Same represented L=12 states and time derivatives, only quadrature changes. No initial-data reconstruction or evolution. This diagnoses projection aliasing, not harmonic-degree convergence or physical acceptance.')
    with output.open('x',encoding='utf-8',newline='\n') as stream:json.dump(record,stream,indent=2,sort_keys=True);stream.write('\n')
    summary=[]
    for comparison in comparisons:
        if comparison['coarse_quadrature']!=list(QUADRATURES[0]) or comparison['reference_quadrature']!=list(QUADRATURES[2]):continue
        summary.append(dict(cache_index=comparison['cache_index'],
            relative_default_to_far={name:entry['relative_coefficient_l2_gap'] for name,entry in comparison['groups'].items()}))
    print(json.dumps(dict(status=status,default_to_far=summary,elapsed_seconds=record['elapsed_seconds'],output=str(output)),indent=2))


if __name__=='__main__':main()
