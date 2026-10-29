"""Retrospective named-preparation controls, not derived birth histories."""
from concurrent.futures import ProcessPoolExecutor
from decimal import Decimal
import argparse,json,math,hashlib,subprocess
from pathlib import Path
from threadpoolctl import threadpool_limits
from measure_nsc_atlas_shared_ruler import authenticate,inputs
from recursive_horizons.nsc_atlas_volume_preparation import PreparedBackground
from recursive_horizons.nsc_atlas_volume_observations import fit
ROOT=Path(__file__).resolve().parents[1]
SOURCES=('src/recursive_horizons/nsc_atlas_volume_preparation.py','scripts/measure_nsc_atlas_exit_preparations.py','tests/test_nsc_atlas_wall_exit_preparation.py')

def worker(job):
    prep,curve,order=job
    with threadpool_limits(limits=1):
        data=inputs();zmax=max(max(data[0]),max(data[4]))
        bg=PreparedBackground(zmax,preparation=prep,curvature=curve,rtol=2e-11 if order==48 else 2e-12)
        r=fit(bg,data,order=order);best=r['joint_best'];m=r['posterior']['moments']['Omega_dark_energy']
        r['three_total_dark_energy']=dict(best_fit=3*best['Omega_dark_energy'],posterior_mean=3*m['mean'],
            posterior_sigma=3*m['standard_deviation'],comparison=2.054,comparison_sigma=.022,
            standardized_difference=(3*m['mean']-2.054)/math.hypot(3*m['standard_deviation'],.022),covariance_included=True)
        r['eos']=bg.eos_rows(best['theta'])
        r['preparation']=dict(s_initial=bg.initial_s,velocity_initial=0.,actual_E_initial=.375,C=bg.C,
            source='HM barrier location with assigned matter' if prep=='HM_location' else 'computed minimum with assigned matter',
            is_derived_tunnelling_exit=False,reference_Lambda_is_future_limit=(prep!='HM_location'))
        endpoint=bg.state([bg.bounds[1]]);se,ve,Je,Ee,Ve,He=(a[0] for a in endpoint)
        barrier=bg.model.vacuum_roots()[0]
        in_basin=bool(se<barrier['s'] and 42*ve*ve+Ve+Ee<barrier['vacuum'])
        if prep!='HM_location' and not in_basin:raise ValueError('minimum continuation not in its claimed basin')
        r['preparation']['future_basin_entry_checked']=in_basin
        r['history']=dict(z_at_reference_slice=math.expm1(best['theta']),minimum_emitter_N=best['theta']-math.log1p(zmax),
            convention='Reference-slice Cauchy continuation, including N<0. This is not a post-birth-only cosmology.',
            new_wall_C_over_E_derived=False)
        return r

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    out=(ROOT/p.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('new repository output required')
    if subprocess.check_output(['git','status','--porcelain','--',*SOURCES],text=True).strip():raise ValueError('commit the comparison method first')
    data,frozen=authenticate()
    jobs=[(p,c,o) for p in ('HM_location','minimum') for c in ('flat','reference') for o in (48,96)]
    with ProcessPoolExecutor(max_workers=6) as pool:rows=list(pool.map(worker,jobs))
    baseline=json.loads((ROOT/'results/summaries/nsc-atlas-imposed-flat-volume-v1.json').read_text(encoding='utf-8'))
    flat_LCDM=baseline['flat_LCDM_control'];models={};controls={}
    for i in range(0,len(rows),2):
        lo,hi=rows[i:i+2];name=hi['model']
        diff={k:abs(lo['joint_best'][k]-hi['joint_best'][k]) for k in ('chi_squared_total','Omega_Lambda','Omega_volume','Omega_m','Omega_k','rd_h_Mpc')}
        diff['three_total_posterior']=abs(lo['three_total_dark_energy']['posterior_mean']-hi['three_total_dark_energy']['posterior_mean'])
        if max(diff.values())>1e-5:raise ValueError('preparation control refinement failed: '+str(diff))
        hi['delta_chi_squared_vs_flat_LCDM']=hi['joint_best']['chi_squared_total']-flat_LCDM['joint_best']['chi_squared_total']
        models[name]=hi;controls[name]=diff
    zmax=2.33
    matched=PreparedBackground(zmax)
    original=baseline['models']['flat'];theta=original['joint_best']['theta']
    r=dict(schema='nsc-atlas-exit-preparations-v1',models=models,numerical_controls=controls,
        matched_flat_baseline=original,matched_flat_eos=matched.eos_rows(theta),flat_LCDM_control=flat_LCDM,
        rounding=dict(displayed_six_decimal_components=['0.577259','0.143326'],
            three_times_displayed_sum=str(3*(Decimal('0.577259')+Decimal('0.143326'))),
            three_times_unrounded_best=original['three_total_dark_energy']['best_fit'],
            posterior_mean=original['three_total_dark_energy']['posterior_mean']),
        data=data,frozen_checkpoints=frozen,
        interpretation='Four named initial-state controls with actual E_i=3/8. No fitted initial displacement; neither a true HM/CDL exit nor a new wall curvature has been derived.',
        method_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        source_sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in SOURCES})
    out.write_text(json.dumps(r,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:dict(joint=v['joint_best'],totalDE=v['three_total_dark_energy'],delta=v['delta_chi_squared_vs_flat_LCDM']) for k,v in models.items()}))
if __name__=='__main__':main()
