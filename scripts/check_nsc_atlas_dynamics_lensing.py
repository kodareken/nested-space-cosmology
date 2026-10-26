"""Prepare, freeze, then separately validate NSC-P004."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np

from check_nsc_atlas_formed_observer import retained_source
from recursive_horizons.nsc_atlas_formed_observer import SourceSlab
from recursive_horizons.nsc_atlas_dynamics_lensing import (
    TARGET,precheck,calibrate,predict_from_motion,evaluate_withheld)

ROOT=Path(__file__).resolve().parents[1]
FORECAST=ROOT/'results/summaries/nsc-p004-dynamics-lensing-forecast-v1.json'
VALIDATION=ROOT/'results/summaries/nsc-p004-dynamics-lensing-validation-v1.json'
INPUTS=('src/recursive_horizons/nsc_atlas_dynamics_lensing.py',
        'src/recursive_horizons/nsc_atlas_formed_observer.py',
        'scripts/check_nsc_atlas_dynamics_lensing.py',
        'scripts/check_nsc_atlas_formed_observer.py',
        'tests/test_nsc_atlas_dynamics_lensing.py')


def slab_for(name):
    r,u,row,digest=retained_source(name)
    return SourceSlab(r,u,length=row['config']['length']),digest


def calibration_worker(job):
    name,source,n,h,release=job;cpu=time.process_time()
    slab,digest=slab_for(source)
    data=calibrate(slab,nodes_per_panel=n,half_width=h,release_error=release)
    result=dict(source=source,source_field_sha256=digest,forecast=predict_from_motion(data),
                cpu_seconds=time.process_time()-cpu)
    if name=='refined':
        result['calibration']=data;result['stationarity_precheck']=precheck(slab)
    else:result['calibration_sha256']=hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest()
    return name,result


def prepare():
    if FORECAST.exists():raise FileExistsError('Forecast already frozen; never overwrite it to match the withheld ray')
    wall=time.perf_counter();cpu=time.process_time()
    jobs=[('coarse','coarse',24,5e-5,0.),('medium','medium',24,5e-5,0.),
          ('fine','fine',24,5e-5,0.),('time_control','time_control',24,5e-5,0.),
          ('outer_control','outer_control',24,5e-5,0.),('quadrature','fine',48,5e-5,0.),
          ('refined','fine',48,2.5e-5,1e-6)]
    with ProcessPoolExecutor(max_workers=4) as pool:runs=dict(pool.map(calibration_worker,jobs))
    f=runs['refined']['forecast'];m=runs['medium']['forecast'];base=runs['fine']['forecast']
    key='predicted_image_shift'
    numerical=2*max(abs(base[key]-m[key]),abs(runs['time_control']['forecast'][key]-m[key]),
        abs(runs['outer_control']['forecast'][key]-m[key]),
        abs(runs['quadrature']['forecast'][key]-base[key]),abs(f[key]-runs['quadrature']['forecast'][key]))
    born=.005*abs(f['normal_frame_lens_contribution'])
    tolerance=numerical+f['release_state_envelope']+born
    mas=180/math.pi*3600*1000
    total=time.process_time()-cpu+sum(r['cpu_seconds'] for r in runs.values())
    record=dict(schema='nsc-p004-dynamics-lensing-forecast-v1',prediction_id='NSC-P004',
        registered_on='2026-10-26',source_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        target=dict(receiver_azimuth=TARGET,initial_observer_radius=.25,initial_observer_velocity=0.,
                    source_proper_age=runs['refined']['calibration']['source_proper_age']),
        registration='prospective_synthetic_withheld_optical_ray',
        withheld_ray_evaluated=False,observational_validation=False,
        prediction=dict(value_radians=f[key],value_milliarcseconds=f[key]*mas,
            comparison_tolerance_radians=tolerance,comparison_tolerance_milliarcseconds=tolerance*mas,
            numerical_envelope=numerical,probe_release_envelope=f['release_state_envelope'],
            born_allowance=born,born_allowance_fraction=.005,
            kind='predeclared comparison tolerance; not a total physical confidence interval'),
        calibration=runs['refined'].pop('calibration'),precheck=runs['refined'].pop('stationarity_precheck'),
        runs=runs,compute_budget_seconds=300,process_cpu_seconds_including_workers=total,
        wall_seconds=time.perf_counter()-wall,
        limits='Synthetic timelike calibration on the fixed weak metric. Observer kinematics are supplied calibration data; their empirical precision is not established. Equal leading potentials are a model premise. No independent lens mass, equilibrium rotation curve, astronomical fit or unique nesting identification.',
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS},
        evidence_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in (
            'results/summaries/nsc-atlas-seed-collapse-v1.json','results/summaries/nsc-atlas-formed-observer-v1.json')})
    if total>300:raise RuntimeError('Preparation exceeded budget')
    # Creation-only, prior to any withheld evaluation.
    with FORECAST.open('x',encoding='utf-8',newline='\n') as stream:
        stream.write(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(forecast=str(FORECAST),prediction=record['prediction'],precheck=record['precheck'],cpu=total)))
    return record


def validation_worker(job):
    name,source,tight=job;cpu=time.process_time();slab,digest=slab_for(source)
    row=evaluate_withheld(slab,rtol=2e-13 if tight else 2e-12,atol=2e-15 if tight else 2e-14)
    if name=='fine':
        age=row['source_proper_age'];checks={}
        for h in (1e-4,5e-5):
            plus=evaluate_withheld(slab,age=age+h,rtol=2e-13,atol=2e-15)
            minus=evaluate_withheld(slab,age=age-h,rtol=2e-13,atol=2e-15)
            derivative=(plus['proper_flight_time']-minus['proper_flight_time'])/(2*h)
            predicted=math.expm1(-row['log_frequency'])
            checks[str(h)]=dict(flight_time_derivative=derivative,frequency_prediction=predicted,gap=derivative-predicted)
        row['moving_clock_identity']=checks
    row.update(source=source,source_field_sha256=digest,cpu_seconds=time.process_time()-cpu)
    return name,row


def validate():
    frozen_bytes=FORECAST.read_bytes();forecast=json.loads(frozen_bytes)
    # A saved Git checkpoint is required before disclosing the held-out ray.
    relative=FORECAST.relative_to(ROOT).as_posix()
    if subprocess.run(['git','ls-files','--error-unmatch',relative],capture_output=True).returncode:
        raise RuntimeError('Commit the frozen forecast before the withheld evaluation')
    stored=subprocess.check_output(['git','show','HEAD:'+relative])
    if stored.replace(b'\r\n',b'\n')!=frozen_bytes.replace(b'\r\n',b'\n'):
        raise RuntimeError('Frozen forecast differs from committed checkpoint')
    wall=time.perf_counter();cpu=time.process_time()
    jobs=[('coarse','coarse',False),('medium','medium',False),('fine','fine',False),
          ('time_control','time_control',False),('outer_control','outer_control',False),('ray_control','fine',True)]
    with ProcessPoolExecutor(max_workers=4) as pool:rows=dict(pool.map(validation_worker,jobs))
    f=rows['fine'];m=rows['medium'];num=2*max(abs(f['image_shift']-m['image_shift']),
        abs(rows['time_control']['image_shift']-m['image_shift']),abs(rows['outer_control']['image_shift']-m['image_shift']),
        abs(rows['ray_control']['image_shift']-f['image_shift']))
    residual=f['image_shift']-forecast['prediction']['value_radians']
    # The frozen acceptance tolerance is not enlarged after seeing the ray.
    passed=abs(residual)<=forecast['prediction']['comparison_tolerance_radians']
    total=time.process_time()-cpu+sum(r['cpu_seconds'] for r in rows.values())
    record=dict(schema='nsc-p004-dynamics-lensing-validation-v1',prediction_id='NSC-P004',
        forecast_sha256=hashlib.sha256(frozen_bytes).hexdigest(),
        forecast_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        runs=rows,withheld_numerical_envelope=num,residual_radians=residual,
        residual_milliarcseconds=residual*180/math.pi*3600*1000,
        frozen_tolerance_radians=forecast['prediction']['comparison_tolerance_radians'],
        synthetic_check_passed=passed,observational_validation=False,
        process_cpu_seconds_including_workers=total,wall_seconds=time.perf_counter()-wall,
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS},
        limits=forecast['limits'])
    VALIDATION.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(actual_mas=f['image_shift']*180/math.pi*3600*1000,
        forecast_mas=forecast['prediction']['value_milliarcseconds'],residual_mas=record['residual_milliarcseconds'],
        passed=passed,cpu=total)))
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare',action='store_true');parser.add_argument('--validate',action='store_true')
    args=parser.parse_args()
    if args.prepare==args.validate:parser.error('choose exactly one of --prepare and --validate')
    prepare() if args.prepare else validate()
