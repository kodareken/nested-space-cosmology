"""Isolate curvature using fixed C=0; this is not an inherited prediction."""
from concurrent.futures import ProcessPoolExecutor
import argparse,json,hashlib,subprocess,math
from pathlib import Path
from threadpoolctl import threadpool_limits
from measure_nsc_atlas_shared_ruler import authenticate,inputs
from recursive_horizons.nsc_atlas_curvature_mechanisms import FlatVolumeControl,ClosedVolumeComparison
from recursive_horizons.nsc_atlas_volume_observations import fit
ROOT=Path(__file__).resolve().parents[1]
SOURCES=('src/recursive_horizons/nsc_atlas_curvature_mechanisms.py','scripts/measure_nsc_atlas_flat_volume_control.py','tests/test_nsc_atlas_curvature_mechanisms.py')

def worker(job):
    name,order=job
    with threadpool_limits(limits=1):
        data=inputs();cls=FlatVolumeControl if name=='flat' else ClosedVolumeComparison
        bg=cls(max(max(data[0]),max(data[4])),rtol=2e-11 if order==48 else 2e-12)
        r=fit(bg,data,order=order)
        m=r['posterior']['moments']['Omega_dark_energy']
        r['three_total_dark_energy']=dict(best_fit=3*r['joint_best']['Omega_dark_energy'],
            posterior_mean=3*m['mean'],posterior_sigma=3*m['standard_deviation'],comparison=2.054,comparison_sigma=.022,
            standardized_difference=(3*m['mean']-2.054)/math.hypot(3*m['standard_deviation'],.022),
            covariance_included=True)
        return r

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('new repository output required')
    if subprocess.check_output(['git','status','--porcelain','--',*SOURCES],text=True).strip():raise ValueError('commit control method first')
    data,frozen=authenticate()
    with ProcessPoolExecutor(max_workers=4) as pool:rows=list(pool.map(worker,[(n,o) for n in ('closed','flat') for o in (48,96)]))
    models={};controls={}
    baseline=json.loads((ROOT/'results/summaries/nsc-atlas-volume-observations-v1.json').read_text(encoding='utf-8'))
    flat_baseline=baseline['models']['flat_LCDM']
    for i,name in enumerate(('closed','flat')):
        lo,hi=rows[2*i:2*i+2]
        diff={k:abs(lo['joint_best'][k]-hi['joint_best'][k]) for k in ('chi_squared_total','Omega_Lambda','Omega_m','Omega_volume','rd_h_Mpc')}
        diff['three_total_posterior']=abs(lo['three_total_dark_energy']['posterior_mean']-hi['three_total_dark_energy']['posterior_mean'])
        if max(diff.values())>1e-5:raise ValueError('control refinement failed')
        hi['delta_chi_squared_vs_flat_LCDM']=hi['joint_best']['chi_squared_total']-flat_baseline['joint_best']['chi_squared_total']
        models[name]=hi;controls[name]=diff
    if abs(models['closed']['joint_best']['chi_squared_total']-baseline['models']['volume_12_2']['joint_best']['chi_squared_total'])>1e-6:raise ValueError('closed baseline not recovered')
    r=dict(schema='nsc-atlas-imposed-flat-volume-v1',models=models,flat_LCDM_control=flat_baseline,
        controls=controls,data=data,frozen_checkpoints=frozen,
        imposed_condition='C=0 exactly, specified before fitting; not inferred from an inherited mechanism.',
        inherited_flatness_derived=False,
        method_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        source_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in SOURCES})
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:dict(joint=v['joint_best'],totalDE=v['three_total_dark_energy'],delta=v['delta_chi_squared_vs_flat_LCDM']) for k,v in models.items()}))
if __name__=='__main__':main()
