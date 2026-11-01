"""Two frozen RHS projections and independent ADM samples; zero time steps."""
import argparse
import ctypes
import gc
import hashlib
import json
import math
import os
from pathlib import Path
import time

import numpy as np
from numba import set_num_threads
from scipy.fft import set_workers
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_parent_gh import ParentGH
from recursive_horizons.nsc_atlas_s3_spectral import sphere_quadrature
from recursive_horizons.nsc_atlas_trapped_parent import GHSurfaceGeometry

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def memory():
    if os.name!='nt':return {}
    class Counters(ctypes.Structure):
        _fields_=[('cb',ctypes.c_ulong),('PageFaultCount',ctypes.c_ulong)]+[
            (name,ctypes.c_size_t) for name in ('PeakWorkingSetSize','WorkingSetSize','QuotaPeakPagedPoolUsage',
                'QuotaPagedPoolUsage','QuotaPeakNonPagedPoolUsage','QuotaNonPagedPoolUsage','PagefileUsage','PeakPagefileUsage')]
    c=Counters();c.cb=ctypes.sizeof(c)
    f=ctypes.windll.psapi.GetProcessMemoryInfo
    f.argtypes=[ctypes.c_void_p,ctypes.c_void_p,ctypes.c_ulong]
    f.restype=ctypes.c_int
    if not f(ctypes.c_void_p(-1),ctypes.byref(c),ctypes.sizeof(c)):return {}
    return dict(peak_working_set_bytes=int(c.PeakWorkingSetSize),working_set_bytes=int(c.WorkingSetSize))


def norms(a,b):
    delta=b-a;n=float(np.linalg.norm(a));m=float(np.max(abs(a)))
    return dict(original_L2=n,oversampled_L2=float(np.linalg.norm(b)),difference_L2=float(np.linalg.norm(delta)),
        relative_L2=float(np.linalg.norm(delta)/n) if n else None,
        maximum_absolute_difference=float(np.max(abs(delta))),original_maximum_absolute=m,
        relative_maximum=float(np.max(abs(delta))/m) if m else None)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,default=Path('.build/valve-counterrotation-evolution-L32-v2/run.json'))
    parser.add_argument('--output',type=Path,default=Path('.build/counterrotation-frozen-RHS-quadrature-reproduced.json'))
    args=parser.parse_args()
    started=time.perf_counter();budget=180.;set_num_threads(4)
    path=(ROOT/args.run).resolve()
    output=(ROOT/args.output).resolve()
    if not path.is_relative_to(ROOT/'.build'):
        raise ValueError('terminal local .build trajectory required')
    if output.exists():raise ValueError('fresh frozen-state output required')
    raw=path.read_bytes();record=json.loads(raw)
    if record['status']!='incomplete' or record.get('producer_source_unchanged_during_run') is not True:
        raise ValueError('terminal unchanged v2 producer required')
    index=max(i for i,s in enumerate(record['snapshots']) if s['sampled_admitted'])
    if abs(record['snapshots'][index]['observer_proper_duration']-1.8)>1e-12:
        raise ValueError('exact last admitted proper1.8 event required')
    counts=[26*98**2,34*130**2]
    estimate=dict(grid_nodes=counts,one_d2g_array_bytes=[n*3*3*4*4*8 for n in counts],
        sequential_full_jets_FFT_peak_estimate_bytes=6*1024**3,
        rationale='Existing jets retain first/second fields and assembled tensors; complex FFT scratch is transient. No dense point-by-mode matrix is built for the collocation RHS.')
    print(json.dumps(dict(event='preallocation_estimate',**estimate)),flush=True)
    with threadpool_limits(limits=4),set_workers(4):
        provider=GHSurfaceGeometry.from_state_cache(ROOT/record['state_cache'],record,snapshot=index)
        original=provider.system;state=provider.state.copy();t=provider.coordinate_t
        if tuple(original.basis.quadrature_shape[:2])!=(26,98):raise ValueError('recorded original quadrature is not26x98x98')
        state_hash=hashlib.sha256(state.tobytes()).hexdigest()
        coefficients=original.unpack(state);coefficient_hash=hashlib.sha256(coefficients.tobytes()).hexdigest()
        refined=ParentGH(32,original.epsilon,gauge=original.gauge,backend=original.backend,
            quadrature_order=(34,130),kick_amplitude=original.kick_amplitude,kick_source=original.kick_source)
        if (refined.basis.labels!=original.basis.labels or
            not np.array_equal(refined.basis.characters,original.basis.characters) or
            not np.array_equal(refined.chars,original.chars) or not np.array_equal(refined.offsets,original.offsets) or
            any((a!=b).nnz for a,b in zip(original.basis.generators,refined.basis.generators))):
            raise ValueError('quadrature-only comparison lost field/ordered-derivative/character identity')
        if not np.array_equal(refined.unpack(state),coefficients):raise ValueError('packed coefficients changed')
        rows=[];derivatives=[]
        for system in (original,refined):
            tick=time.perf_counter()
            system.started=started;system.budget=budget
            derivative=system.rhs(t,state)
            if not np.all(np.isfinite(derivative)):raise ValueError('nonfinite frozen RHS')
            derivatives.append(derivative)
            row=dict(quadrature_shape=list(system.basis.quadrature_shape),nodes=system.N,
                elapsed_seconds=time.perf_counter()-tick,instantaneous_RHS_calls=system.calls,memory=memory())
            rows.append(row);print(json.dumps(dict(event='frozen_RHS',**row)),flush=True)
            gc.collect()
            if time.perf_counter()-started>budget:raise TimeoutError('frozen RHS allowance exhausted')
        a,b=(original.unpack(v) for v in derivatives)
        blocks={'g':range(0,10),'v':range(10,20),'S':[20],'chi':[21],'pS':[22],'pchi':[23]}
        comparisons={name:norms(a[list(indices)],b[list(indices)]) for name,indices in blocks.items()}
        comparisons['observer10']=norms(derivatives[0][-10:],derivatives[1][-10:])
        independent=None
        if time.perf_counter()-started<budget-15:
            tick=time.perf_counter()
            points,weights=sphere_quadrature(8,8,17)
            adm=provider.adm_sample(points)
            h=provider(points)['h'];hi=np.linalg.inv(h)
            measure=weights*np.sqrt(np.linalg.det(h))/original.r**3;V=float(np.sum(measure))
            M2=np.maximum(0.,np.einsum('ni,nij,nj->n',adm['M'],hi,adm['M']))
            Ku=np.einsum('nij,njk->nik',hi,adm['Kij']);trace=np.trace(Ku,axis1=1,axis2=2)
            K2=np.einsum('nij,nji->n',Ku,Ku)
            scale=original.r**2+abs(adm['R3'])+trace**2+K2+2*original.parent.kappa*abs(adm['rho'])
            independent=dict(points=len(points),node_rule='Fixed independent sphere_quadrature(8,8,17), no collocation-node reuse.',
                node_sha256=hashlib.sha256(points.tobytes()).hexdigest(),
                physical_momentum_RMS_over_H2=float(np.sqrt(measure@M2/V)/original.r**2),
                physical_momentum_max_over_H2=float(np.sqrt(np.max(M2))/original.r**2),
                maximum_relative_momentum=float(np.max(np.sqrt(M2)/scale)),
                physical_Hamiltonian_RMS_over_H2=float(np.sqrt(measure@(adm['H']**2)/V)/original.r**2),
                elapsed_seconds=time.perf_counter()-tick,
                scope='Actual constraint defect of the unchanged represented state on modest independent nodes, not a continuum maximum/error bound. Quadrature changes do not repair this existing state instantaneously.')
        result=dict(schema='nsc-counterrotation-frozen-RHS-quadrature-v1',
            authenticated_run=path.relative_to(ROOT).as_posix(),run_record_sha256=hashlib.sha256(raw).hexdigest(),
            producer_commit=record['producer_commit'],state_cache_sha256=record['state_cache_sha256'],
            snapshot_index=index,coordinate_t=t,proper=float(state[-2]),
            state_array_sha256=state_hash,coefficient_array_sha256=coefficient_hash,
            exact_state_unchanged=state_hash==hashlib.sha256(state.tobytes()).hexdigest(),
            exact_coefficients_and_observer_preserved=np.array_equal(original.unpack(state),refined.unpack(state)),
            ordered_frame_generators_and_characters_identical=True,
            original_and_refined=rows,block_comparisons=comparisons,
            coefficient_norm_scope='Euclidean norm of packed real harmonic coefficients in the same regular reference frame; metric symmetric slots counted once, not a covariant spacetime tensor norm.',
            independent_ADM=independent,preallocation_memory_estimate=estimate,
            final_memory=memory(),elapsed_seconds=time.perf_counter()-started,budget_seconds=budget,threads=4,
            analysis_producer_sha256=sha(__file__),source_sha256=record['source_sha256'],
            source_unchanged=all(sha(ROOT/p)==v for p,v in record['source_sha256'].items()),
            initializers_called=0,integrator_steps=0,PDE_evolution_performed=False,
            scope='Instantaneous quadrature sensitivity at one admitted frozen state; it separates projection aliasing here from field changes, not accumulated temporal/spatial/constraint error over the history.')
        output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
        print(json.dumps(dict(event='result',block_comparisons=comparisons,independent_ADM=independent,
            elapsed=result['elapsed_seconds'],source_unchanged=result['source_unchanged'],output=str(output),sha256=sha(output))),flush=True)


if __name__=='__main__':main()
