"""Inspect reflection parity of saved represented GH states, without evolution.

R3=diag(1,1,-1,1) leaves the preparation a*X1+b*X2*X4 invariant.  The
reflection exchanges left/right invariant frames, so a metric comparison
must include its position-dependent frame pullback rather than component
sign changes.  These are numerical representation diagnostics only.
"""
import argparse
import json
import math
from pathlib import Path
import time

import numpy as np
from numba import get_num_threads,set_num_threads
from scipy.fft import set_workers
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_parent_gh import ParentGH,PAIRS,FRAME_MATRICES

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_CACHES=(
    '.build/gh-controlled-mid-L8.npz',
    '.build/gh-controlled-mid-L12.npz',
    '.build/gh-controlled-e1-L12.npz',
)


def odd_power(values,permutation,weights):
    odd=(values-values[permutation])/2
    total=float(weights@np.sum(values*values,axis=-1))
    power=float(weights@np.sum(odd*odd,axis=-1))
    return dict(odd_L2_power_fraction=power/total if total else 0.,
                odd_L2_norm=math.sqrt(power),total_L2_norm=math.sqrt(total),
                zero_total_norm=total==0.)


def metric_defect(g,permutation,T,weights):
    reflected=np.einsum('nca,ncd,ndb->nab',T,g[permutation],T)
    difference=g-reflected
    error=float(weights@np.sum(difference*difference,axis=(1,2)))
    total=float(weights@np.sum(g*g,axis=(1,2)))
    return math.sqrt(error/total) if total else 0.


def reflection_data(basis):
    shape=basis.quadrature_shape
    reverse=(-np.arange(shape[1]))%shape[1]
    permutation=np.arange(len(basis.points)).reshape(shape)[:,reverse,:].reshape(-1)
    R=np.array([1.,1.,-1.,1.]);X=basis.points;RX=X*R
    permutation_error=float(np.max(abs(X[permutation]-RX)))
    if permutation_error>3e-14:raise ArithmeticError('reflection is not the stated grid permutation')
    frame=np.einsum('iab,nb->nia',FRAME_MATRICES,X)
    reflected_frame=np.einsum('iab,nb->nia',FRAME_MATRICES,RX)
    T=np.zeros((len(X),4,4));T[:,0,0]=1.
    T[:,1:,1:]=np.einsum('nja,nia->nji',reflected_frame,frame*R)
    round_metric=np.broadcast_to(np.diag([-1.,1.,1.,1.]),(len(X),4,4)).copy()
    psi=1+.03*X[:,0]**2+.02*X[:,1]*X[:,3]
    conformal=round_metric.copy();conformal[:,1:,1:]*=psi[:,None,None]**4
    controls=dict(permutation_maximum_error=permutation_error,
                  round_metric_relative_pullback_RMS=metric_defect(round_metric,permutation,T,basis.weights),
                  even_conformal_metric_relative_pullback_RMS=metric_defect(conformal,permutation,T,basis.weights))
    if max(controls.values())>5e-14:raise ArithmeticError('reflection pullback control failed')
    return permutation,T,controls


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',action='append',type=Path)
    parser.add_argument('--budget-seconds',type=float,default=120.)
    args=parser.parse_args()
    if not 0<args.budget_seconds<=120:raise ValueError('require a positive budget <=120 seconds')
    started=time.perf_counter();old_threads=get_num_threads();records=[]
    def budget():
        if time.perf_counter()-started>args.budget_seconds:raise TimeoutError('reflection diagnostic budget exhausted')
    try:
        set_num_threads(2)
        with threadpool_limits(limits=2),set_workers(2):
            for name in args.cache or DEFAULT_CACHES:
                budget();path=(ROOT/name).resolve()
                if not path.is_relative_to(ROOT/'.build'):raise ValueError('require local .build coefficient cache')
                with np.load(path,allow_pickle=False) as saved:
                    metadata=json.loads(str(saved['metadata']))
                    times=np.array(saved['times'],copy=True);states=np.array(saved['states'],copy=True)
                if metadata.get('backend','dense')!='hopf':raise ValueError('exact-grid diagnosis requires Hopf cache')
                system=ParentGH(int(metadata['degree']),float(metadata['epsilon']),gauge=metadata['gauge'],backend='hopf')
                if states.ndim!=2 or states.shape[1]!=system.dof+10 or len(times)!=len(states):
                    raise ValueError('cached packed state has incompatible coefficient layout')
                permutation,T,controls=reflection_data(system.basis)
                omega=system.parent.reference_frequency;ms=math.sqrt(system.minimum['stiffness']/system.K)
                weights=system.basis.weights;X3=system.basis.points[:,2];x3norm=float(weights@(X3*X3))
                snapshots=[]
                for t,y in zip(times,states):
                    budget();C=system.unpack(y)
                    fields=system.basis.evaluate(C[20:24].T)
                    volume=fields[:,[0,2]]/np.array([1.,ms]);chi=fields[:,[1,3]]/np.array([1.,omega])
                    values=system.basis.evaluate(C[:10].T)
                    g=np.zeros((system.N,4,4))
                    for k,(a,b) in enumerate(PAIRS):g[:,a,b]=g[:,b,a]=values[:,k]
                    dipole=float(weights@(fields[:,0]*X3)/x3norm)
                    snapshots.append(dict(coordinate_t=float(t),observer_proper_duration=float(y[system.dof+8]),
                        chi=odd_power(chi[:,:1],permutation,weights),
                        chi_coordinate_rate_over_omega=odd_power(chi[:,1:],permutation,weights),
                        chi_quadrature_pair=odd_power(chi,permutation,weights),
                        volume_S=odd_power(volume[:,:1],permutation,weights),
                        volume_coordinate_rate_over_ms=odd_power(volume[:,1:],permutation,weights),
                        volume_quadrature_pair=odd_power(volume,permutation,weights),
                        volume_S_X3_linear_coefficient=dipole,
                        volume_s_X3_linear_coefficient=dipole/math.sqrt(system.K),
                        metric_relative_pullback_RMS=metric_defect(g,permutation,T,weights)))
                records.append(dict(cache=path.relative_to(ROOT).as_posix(),degree=system.basis.degree,
                    epsilon=system.epsilon,gauge=system.gauge,controls=controls,snapshots=snapshots))
    finally:set_num_threads(old_threads)
    result=dict(schema='nsc-atlas-gh-reflection-diagnosis-v1',threads=2,
        elapsed_seconds=time.perf_counter()-started,records=records,
        norm='Reference-round L2 power and component Frobenius pullback RMS; not physical energy.',
        scope='Saved represented states only. No initial-data reconstruction, PDE evolution or parity filter.')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
