"""One frozen scalar evaluation of the constrained vacuum gravitational reference."""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
import hashlib
import json
from math import isclose, isfinite, sqrt
from pathlib import Path
import subprocess
import time

from recursive_horizons.nsc_spherical_vacuum_response import vacuum_coefficients, vacuum_exterior_response


ROOT=Path(__file__).resolve().parents[1]
CAUSAL=ROOT/'results/summaries/nsc-conformal-causal-response-v2.json'
RESOLUTION=ROOT/'results/summaries/nsc-conformal-resolution-v1.json'
ACTION='src/recursive_horizons/nsc_spherical_feedback_action.py'
CPU_LIMIT=5.0
BYTE_LIMIT=12000
SOURCES=(
    'scripts/measure_nsc_vacuum_gravitational_response.py',
    'src/recursive_horizons/nsc_spherical_vacuum_response.py',
    'tests/test_nsc_spherical_vacuum_response.py',
    'src/recursive_horizons/nsc_spherical_timelike_boundary.py',
    'tests/test_nsc_spherical_timelike_boundary.py',ACTION,
    CAUSAL.relative_to(ROOT).as_posix(),RESOLUTION.relative_to(ROOT).as_posix(),
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen_sources():
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    hashes={}
    for name in SOURCES:
        raw=(ROOT/name).read_bytes()
        if subprocess.check_output(['git','show',f'{commit}:{name}'],cwd=ROOT)!=raw:
            raise ValueError('freeze the vacuum response method before application: '+name)
        hashes[name]=hashlib.sha256(raw).hexdigest()
    return commit,hashes


def load_parameters(causal_path,resolution_path):
    causal=json.loads(causal_path.read_text(encoding='utf-8'))
    resolution=json.loads(resolution_path.read_text(encoding='utf-8'))
    if causal['status']!='completed_causal_photon_response' or resolution['status']!='measured_instantaneous_correlation':
        raise ValueError('the two frozen source calculations must be complete')
    if causal['input']['payload_sha256']!=resolution['input']['payload_sha256']:
        raise ValueError('the source geometries differ')
    for record in (causal,resolution):
        if record['source_sha256'][ACTION]!=digest(ROOT/ACTION):
            raise ValueError('owned action differs from the reference input')
    A,CW=(float(causal['couplings'][name]) for name in ('A','C_W'))
    for name,value in (('A',A),('C_W',CW)):
        if not isclose(value,float(resolution['parameters'][name]),rel_tol=1e-14,abs_tol=0):
            raise ValueError('source couplings differ: '+name)
    reference=resolution['domain']
    radius=float(reference['reference_length'])
    if reference['target_observation_used'] is not False or not (isfinite(radius) and radius>0):
        raise ValueError('an independently fixed positive ruler is required')
    return A,CW,radius


def real_scalar(value):
    value=complex(value)
    if not (isfinite(value.real) and isfinite(value.imag)) or abs(value.imag)>1e-12*max(1,abs(value.real)):
        raise ValueError('the selected real Laplace point must return finite real coefficients')
    return value.real


def prediction(A,CW,R0):
    # Selection fixed in the admission, independent of any target observation.
    mass2=vacuum_coefficients(A,CW)['mass_squared']
    s=2*sqrt(abs(mass2)+1/R0**2)
    response=vacuum_exterior_response(s,R0,A,CW)
    mu2=response['mass_squared']*R0**2
    growth=response['growth_bound']
    with localcontext() as context:
        context.prec=50
        ad,cd,rd=(Decimal(str(v)) for v in (A,CW,R0))
        independent=ad*rd*rd/(2*cd)
        gap=abs(mu2-float(independent))
    if gap>1e-12*max(1,abs(mu2)):
        raise ValueError('independent decimal spectral-scale check failed')
    return {
        'mass_squared_Lref2':mu2,
        'quadratic_kinetic_coefficient_over_Lref2':response['quadratic_kinetic_coefficient']/R0**2,
        'asymptotic_growth_rate_Lref':growth*R0,
        'growth_time_over_Lref':None if growth==0 else 1/(growth*R0),
        'stationary_stable_flat_reference':mass2>=0,
        'laplace_s_Lref':s*R0,
        'decay_q_Lref':real_scalar(response['decay_q']*R0),
        'chi_R_per_chi_times_Lref':real_scalar(response['chi_R_per_chi']*R0),
        'outward_chi_normal_per_chi_times_Lref':real_scalar(response['outward_chi_normal_per_chi']*R0),
        'radial_metric_b_per_chi':real_scalar(response['radial_metric_b_per_chi']),
        'lapse_n_per_chi':real_scalar(response['lapse_n_per_chi']),
        'energy_response_times_Lref':real_scalar(response['unreferenced_energy_per_chi']*R0),
        'area_pressure_response_times_Lref3':real_scalar(response['unreferenced_area_pressure_per_chi']*R0**3),
        'auxiliary_proper_response_times_Lref':real_scalar(response['auxiliary_proper_response_per_chi']*R0),
        'independent_decimal_mass_scale_gap':gap,
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    if not args.output.parent.is_dir():raise ValueError('output parent must exist')
    started=time.process_time()
    commit,bindings=frozen_sources()
    A,CW,R0=load_parameters(CAUSAL,RESOLUTION)
    measured=prediction(A,CW,R0)
    elapsed=time.process_time()-started
    if elapsed>CPU_LIMIT:raise RuntimeError('scalar reference exhausted its CPU allowance')
    result={
        'schema':'nsc-vacuum-gravitational-response-v1',
        'status':'completed_constrained_vacuum_reference',
        'source_commit':commit,'source_sha256':bindings,
        'parameters':{'A':A,'C_W':CW,'changed':False},
        'reference':{'R0_equals_Lref':R0,'ruler_origin':'physical-resolution reference_T=0.05, x=1',
                     'target_observation_used':False,'attached_to_charged_episode':False},
        'exterior':{'background':'flat physical N=a=1, r=R, chi=0',
                    'domain':'R>=R0, linear perturbations',
                    'state':'zero initial perturbation; zero mass constant; no incoming infinity',
                    'Maxwell_flux':0,'Dirac_field':0,'outward_orientation_at_R0':-1},
        'prediction':measured,
        'scope':('Linear reference response under the same A,C_W. The auxiliary mode is retained '
                 'nonperturbatively in this calculation. A negative squared scale admits growing '
                 'long-wavelength exterior modes and no stable stationary reduction about this flat '
                 'state. An EFT reading needs its order reduction or causal completion specified. '
                 'No claim about stability of the distinct finite periodic episode, physical parent '
                 'attachment, cosmological Lambda, or observational agreement is made.'),
        'uncertainty':('Analytic scalar formulas evaluated in floating point; a 50-digit decimal '
                       'spectral-scale check is recorded. This is not a physical/model uncertainty bound.'),
        'counts':{'trajectories':0,'root_solves':0,'grid_scans':0,'parameter_fits':0,'scalar_cases':1},
        'runtime_cpu_seconds':elapsed,'limits':{'CPU_seconds':CPU_LIMIT,'summary_bytes':BYTE_LIMIT},
        'reproduce':'python scripts/run_python.py scripts/measure_nsc_vacuum_gravitational_response.py --output .build/vacuum-gravity-reference.json',
    }
    raw=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(raw.encode())>BYTE_LIMIT:raise ValueError('summary exceeds its admitted size')
    args.output.write_text(raw,encoding='utf-8',newline='\n')
    print(f'Wrote {args.output}; stationary stable flat reference: {measured["stationary_stable_flat_reference"]}; CPU {elapsed:.6f}s')
    return 0


if __name__=='__main__':raise SystemExit(main())
