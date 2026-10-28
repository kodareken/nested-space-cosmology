"""Freeze SN-calibrated AP shapes, then test pinned DESI DR2 null contrasts.

The first command never loads a BAO mean or covariance. Validation uses the
unchanged source preparation and frozen SN parameters, with no fitted ruler.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from urllib.request import urlopen

import numpy as np

from recursive_horizons.nsc_atlas_bao_shape import ap_prediction,ap_null_contrast
from recursive_horizons.nsc_atlas_parent_supernova import parent_fractions,parent_ray_angle_upper
from recursive_horizons.nsc_atlas_parent_past import past_bound

ROOT=Path(__file__).resolve().parents[1]
RELEASE='bb0c1c9009dc76d1391300e169e8df38fd1096db'
REDSHIFTS=np.array([.51,.706,.934,1.321,1.484,2.33])
FILES={
    'mean.txt':('desi_gaussian_bao_ALL_GCcomb_mean.txt',472,'8aff444fdb42c0946342aa0011ab287eda097c4c'),
    'cov.txt':('desi_gaussian_bao_ALL_GCcomb_cov.txt',2547,'fd8e5697ab61379b07b52efb781ea6713417a4d9')}
SOURCES=('scripts/check_nsc_atlas_parent_bao.py','src/recursive_horizons/nsc_atlas_bao_shape.py',
         'src/recursive_horizons/nsc_atlas_parent_past.py','tests/test_nsc_atlas_bao_shape.py')


def freeze():
    source=ROOT/'results/summaries/nsc-atlas-parent-supernova-v2.json'
    record=json.loads(source.read_text())
    models=record['comparisons']['all_noncalibrators_z_ge_0_025']['models']
    parent=models['parent_calibrated_epoch']['parameters']['inverse_Ao']
    f=parent_fractions(parent)
    matter=models['flat_LCDM_calibrated']['parameters']['Omega_m']
    if 1/(parent*(1+max(REDSHIFTS)))<1/3:
        raise ValueError('BAO emitter outside proved Cauchy collar')
    b=past_bound();r=float(b['hubble_squared_relative_error'])
    angle=parent_ray_angle_upper(max(REDSHIFTS),parent,minimum_scale=1/3)
    if angle>=1:raise ValueError('monotone sine error comparison unavailable')
    return dict(schema='nsc-atlas-parent-bao-forecast-v1',
        SN_input=dict(path=source.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(source.read_bytes()).hexdigest()),
        redshifts=REDSHIFTS.tolist(),
        models=dict(parent=dict(parameters=dict(inverse_Ao=parent,**f),
                    predicted_AP=ap_prediction(REDSHIFTS,matter=f['matter'],curvature=f['curvature']).tolist(),
                    scalar_relative_AP_error_bound=r/(1-r),ray_angle_upper=angle),
                    flat_LCDM=dict(parameters=dict(matter=matter,curvature=0),
                    predicted_AP=ap_prediction(REDSHIFTS,matter=matter).tolist())),
        interpretation='SN-calibrated geometric AP shapes; not a generated acoustic ruler or baryon mechanism. Conditional on frozen SN point estimates; no joint significance or selection of absolute scale.',
        BAO_usage='Not loaded by this freeze command; public data are retrospective, not newly collected blinded evidence.')


def validate(forecast):
    raw=forecast.read_bytes();f=json.loads(raw)
    cache=ROOT/'.build/data/desi-dr2-bb0c1c9009dc';cache.mkdir(parents=True,exist_ok=True)
    bindings={}
    for key,(filename,size,blob) in FILES.items():
        path=cache/key
        url=f'https://raw.githubusercontent.com/CobayaSampler/bao_data/{RELEASE}/desi_bao_dr2/{filename}'
        if path.is_symlink():raise ValueError('regular cache file required')
        data=path.read_bytes() if path.exists() else urlopen(url,timeout=30).read(size+1)
        got=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        if len(data)!=size or got!=blob:raise ValueError('pinned DESI input differs')
        if not path.exists():path.write_bytes(data)
        bindings[key]=dict(url=url,release_commit=RELEASE,bytes=size,git_blob=blob,sha256=hashlib.sha256(data).hexdigest())
    rows=np.genfromtxt(cache/'mean.txt',dtype=None,encoding='utf-8')
    z=np.array([row[0] for row in rows]);y=np.array([row[1] for row in rows]);labels=[row[2] for row in rows]
    C=np.loadtxt(cache/'cov.txt')
    pairs=[]
    for zz in f['redshifts']:
        dm=[i for i,(red,label) in enumerate(zip(z,labels)) if abs(red-zz)<1e-8 and label=='DM_over_rs']
        dh=[i for i,(red,label) in enumerate(zip(z,labels)) if abs(red-zz)<1e-8 and label=='DH_over_rs']
        if len(dm)!=1 or len(dh)!=1:raise ValueError('one DM/DH pair per frozen redshift required')
        pairs.append([dm[0],dh[0]])
    return dict(schema='nsc-atlas-parent-bao-validation-v1',
        forecast=dict(path=forecast.relative_to(ROOT).as_posix(),sha256=hashlib.sha256(raw).hexdigest()),
        sources=bindings,models={name:ap_null_contrast(y,C,pairs,row['predicted_AP']) for name,row in f['models'].items()},
        redshifts=f['redshifts'],fitted_BAO_parameters=0,
        interpretation='Conditional six-null-contrast test with full published covariance. SN calibration uncertainty is not marginalized; no global significance or baryon-origin claim. No sound-ruler calibration or BAO retuning.',
        excluded_row='Isotropic BGS DV/r_d has no ruler-free AP pair and is not used.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--validate',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    output=(ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('new repository output required')
    if subprocess.check_output(['git','status','--porcelain','--',*SOURCES],cwd=ROOT,text=True).strip():
        raise ValueError('commit AP method before freeze or validation')
    result=validate((ROOT/args.validate).resolve()) if args.validate else freeze()
    result['method_commit']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    result['source_sha256']={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES}
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result))


if __name__=='__main__':main()
