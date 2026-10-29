"""Refit the evolving volume candidate with the existing released covariance."""
from concurrent.futures import ProcessPoolExecutor
import argparse,hashlib,json,subprocess
from pathlib import Path
from threadpoolctl import threadpool_limits

from measure_nsc_atlas_shared_ruler import authenticate,inputs
from recursive_horizons.nsc_atlas_volume_observations import VolumeBackground,FlatBackground,fit

ROOT=Path(__file__).resolve().parents[1]
SOURCES=('src/recursive_horizons/nsc_atlas_volume_observations.py',
         'scripts/measure_nsc_atlas_volume_observations.py','tests/test_nsc_atlas_volume_observations.py')


def worker(job):
    name,order=job
    with threadpool_limits(limits=1):
        data=inputs()
        bg=VolumeBackground(max(max(data[0]),max(data[4])),rtol=2e-11 if order==48 else 2e-12) if name=='volume_12_2' else FlatBackground()
        result=fit(bg,data,order=order)
        if name=='volume_12_2':
            result['background_domain']=dict(minimum_N=bg.minimum_N,maximum_N=bg.bounds[1],RHS_calls=bg.calls,
                observer_convention='Einstein-frame metric redshift and fixed luminosity standards',
                scope='Averaged homogeneous ODE continuation; no shared-wall domain or full scalar observational error bound is claimed.')
        return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('new repository output required')
    if subprocess.check_output(['git','status','--porcelain','--',*SOURCES],text=True).strip():
        raise ValueError('commit the method before the retained retrospective comparison')
    data,frozen=authenticate()
    jobs=[(name,order) for name in ('volume_12_2','flat_LCDM') for order in (48,96)]
    with ProcessPoolExecutor(max_workers=4) as pool:rows=list(pool.map(worker,jobs))
    models={};controls={}
    for i,name in enumerate(('volume_12_2','flat_LCDM')):
        low,high=rows[2*i:2*i+2]
        diff={key:abs(low['joint_best'][key]-high['joint_best'][key]) for key in (
            'chi_squared_total','Omega_Lambda','Omega_m','Omega_k','Omega_volume','rd_h_Mpc')}
        diff['posterior_three_Omega']=abs(low['three_Omega_Lambda']['posterior_mean']-high['three_Omega_Lambda']['posterior_mean'])
        if max(diff.values())>1e-5:raise ValueError('numerical inference control failed: '+str(diff))
        controls[name]=diff;models[name]=high
    old=json.loads((ROOT/'results/summaries/nsc-atlas-shared-ruler-v1.json').read_text(encoding='utf-8'))
    flat=models['flat_LCDM']['joint_best']['chi_squared_total']
    if abs(flat-old['models']['flat_LCDM']['joint_best']['chi_squared_total'])>1e-6:raise ValueError('common flat control not recovered')
    r=dict(schema='nsc-atlas-volume-observations-v1',models=models,numerical_controls=controls,
        delta_chi_squared=models['volume_12_2']['joint_best']['chi_squared_total']-flat,
        data=data,frozen_checkpoints=frozen,
        Omega_definition='Omega_Lambda is the constant late floor Vstar/(3H_o^2); Omega_volume includes volume kinetic energy and V(s)-Vstar. Four fractions sum to one.',
        interpretation='Retrospective conditional test of the unselected (12,2) averaged homogeneous candidate. No SI scale, ancestry selection or complete observer/microphysical source is derived.',
        method_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        source_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in SOURCES})
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(delta_chi_squared=r['delta_chi_squared'],models={k:v['joint_best'] for k,v in models.items()},controls=controls)))

if __name__=='__main__':main()
