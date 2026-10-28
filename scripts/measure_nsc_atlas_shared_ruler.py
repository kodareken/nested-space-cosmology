"""One common DESI DR2 ruler, full covariance and propagated SN shape posterior."""
from concurrent.futures import ProcessPoolExecutor
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from threadpoolctl import threadpool_limits

from measure_nsc_atlas_distance_contrast import cached_file,FILES as SN_FILES
from check_nsc_atlas_parent_bao import FILES as BAO_FILES,RELEASE as BAO_RELEASE
from recursive_horizons.nsc_atlas_shared_ruler import posterior_calculation

ROOT=Path(__file__).resolve().parents[1]
SOURCES=('src/recursive_horizons/nsc_atlas_shared_ruler.py','scripts/measure_nsc_atlas_shared_ruler.py',
         'tests/test_nsc_atlas_shared_ruler.py')


def authenticate():
    sn_cache=ROOT/'.build/data/pantheonplus-c447f0fea703'
    bindings={name:cached_file(sn_cache,name,False) for name in SN_FILES}
    for key,(name,size,blob) in BAO_FILES.items():
        path=ROOT/'.build/data/desi-dr2-bb0c1c9009dc'/key
        data=path.read_bytes()
        actual=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        if len(data)!=size or actual!=blob:raise ValueError('retained DESI data changed')
        bindings['DESI_'+key]=dict(release=BAO_RELEASE,source_path=name,bytes=size,git_blob=blob,
                                  sha256=hashlib.sha256(data).hexdigest())
    # The physical source and frozen forecasts must remain those at 485b3b9.
    frozen=['results/prediction-register.json','src/recursive_horizons/nsc_atlas_shared_scalar_slice.py',
            'src/recursive_horizons/nsc_atlas_flux_conversion.py','src/recursive_horizons/nsc_atlas_parent_supernova.py',
            'src/recursive_horizons/nsc_atlas_parent_past.py',
            'results/summaries/nsc-atlas-parent-supernova-v2.json',
            'results/summaries/nsc-atlas-parent-bao-validation-v1.json']
    for name in frozen:
        original=subprocess.check_output(['git','show','485b3b9:'+name],cwd=ROOT)
        if (ROOT/name).read_bytes()!=original:raise ValueError('fixed checkpoint input changed: '+name)
    return bindings,{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in frozen}


def inputs():
    sn_cache=ROOT/'.build/data/pantheonplus-c447f0fea703'
    table=np.genfromtxt(sn_cache/'distances.dat',names=True,dtype=None,encoding='utf-8')
    C=np.loadtxt(sn_cache/'stat-plus-systematic.cov',skiprows=1).reshape(1701,1701)
    selected=np.flatnonzero((table['IS_CALIBRATOR']==0)&(table['zHD']>=.025))
    if len(selected)!=1322:raise ValueError('fixed supernova selection changed')
    bp=ROOT/'.build/data/desi-dr2-bb0c1c9009dc'
    rows=np.genfromtxt(bp/'mean.txt',dtype=None,encoding='utf-8')
    return (table['zHD'][selected],table['zHEL'][selected],table['m_b_corr'][selected],C[np.ix_(selected,selected)],
            np.array([r[0] for r in rows]),[r[2] for r in rows],np.array([r[1] for r in rows]),np.loadtxt(bp/'cov.txt'))


def worker(job):
    model,order=job
    with threadpool_limits(limits=1):
        return posterior_calculation(model,*inputs(),order=order)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();output=(ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('new repository output required')
    if subprocess.check_output(['git','status','--porcelain','--',*SOURCES],cwd=ROOT,text=True).strip():
        raise ValueError('commit the method before evaluating the new data comparison')
    bindings,frozen=authenticate()
    jobs=[(m,o) for m in ('parent','flat_LCDM') for o in (24,48)]
    with ProcessPoolExecutor(max_workers=4) as pool:rows=list(pool.map(worker,jobs))
    results={};controls={}
    prior=json.loads((ROOT/'results/summaries/nsc-atlas-parent-supernova-v2.json').read_text(encoding='utf-8'))
    for i,model in enumerate(('parent','flat_LCDM')):
        low,high=rows[2*i:2*i+2]
        differences=dict(joint_chi_squared=abs(low['joint_best']['chi_squared_total']-high['joint_best']['chi_squared_total']),
            marginal_BAO_deviance=abs(low['posterior_propagation']['marginal_BAO_deviance']-high['posterior_propagation']['marginal_BAO_deviance']),
            posterior_Omega_Lambda_mean=abs(low['posterior_propagation']['joint_fraction_moments']['Omega_Lambda']['mean']-
                                            high['posterior_propagation']['joint_fraction_moments']['Omega_Lambda']['mean']),
            rd_h_median_Mpc=abs(low['posterior_propagation']['joint_rd_h_quantiles_16_50_84_Mpc'][1]-
                                 high['posterior_propagation']['joint_rd_h_quantiles_16_50_84_Mpc'][1]))
        if max(differences.values())>1e-5:raise ValueError('posterior quadrature control failed: '+str(differences))
        old_model='parent_calibrated_epoch' if model=='parent' else 'flat_LCDM_calibrated'
        old_chi=prior['comparisons']['all_noncalibrators_z_ge_0_025']['models'][old_model]['chi_squared']
        if abs(high['SN_calibration']['chi_squared']-old_chi)>1e-6:raise ValueError('fixed SN score not recovered')
        results[model]=high;controls[model]=differences
    ns,lc=(results[name] for name in ('parent','flat_LCDM'))
    result=dict(schema='nsc-atlas-shared-ruler-v1',checkpoint='485b3b9',models=results,quadrature_controls=controls,
        delta_chi_squared_joint=ns['joint_best']['chi_squared_total']-lc['joint_best']['chi_squared_total'],
        delta_marginal_BAO_deviance=ns['posterior_propagation']['marginal_BAO_deviance']-lc['posterior_propagation']['marginal_BAO_deviance'],
        fixed_checkpoints=frozen,data=bindings,
        interpretation='Retrospective full DESI shared-ruler comparison. Source couplings/normalization and P004/P005 unchanged. One observational ruler is calibrated; no SI vacuum is derived.',
        method_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES})
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=output.relative_to(ROOT).as_posix(),delta_chi_squared_joint=result['delta_chi_squared_joint'],
                          delta_marginal_BAO_deviance=result['delta_marginal_BAO_deviance'],
                          models={k:dict(joint=v['joint_best'],posterior=v['posterior_propagation'],
                                         three_Omega=v['three_Omega_Lambda_comparison']) for k,v in results.items()})))


if __name__=='__main__':main()
