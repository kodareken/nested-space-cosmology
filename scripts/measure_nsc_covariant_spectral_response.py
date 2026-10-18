"""One fixed-scale curvature-channel and static-probe response application."""
from __future__ import annotations

import argparse
import hashlib
import json
from math import isclose,pi
from pathlib import Path
import subprocess
import time

from threadpoolctl import threadpool_info,threadpool_limits

import measure_nsc_spectral_gravitational_response as previous
from recursive_horizons.dirac.compact import matched_weight
from recursive_horizons.nsc_covariant_spectral_response import (
    curvature_profile_kernels,static_gravitational_response,
)
from recursive_horizons.nsc_scale_binding import ScaleBindingConventions
from recursive_horizons.nsc_warped_source import WarpedCompactWeight


ROOT=Path(__file__).resolve().parents[1]
TT='results/summaries/nsc-spectral-gravitational-response-v1.json'
CPU_LIMIT=120.0
BYTE_LIMIT=12000
EXTRA_SOURCES=(
    'scripts/measure_nsc_covariant_spectral_response.py',
    'src/recursive_horizons/nsc_covariant_spectral_response.py',
    'tests/test_nsc_covariant_spectral_response.py',
    'src/recursive_horizons/nsc_covariant_identities.py',TT,
)


def frozen_sources():
    commit,bindings=previous.frozen_sources()
    for name in EXTRA_SOURCES:
        raw=(ROOT/name).read_bytes()
        if subprocess.check_output(['git','show',f'{commit}:{name}'],cwd=ROOT)!=raw:
            raise ValueError('freeze the covariant response method before application: '+name)
        bindings[name]=hashlib.sha256(raw).hexdigest()
    return commit,bindings


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    if not args.output.parent.is_dir():raise ValueError('output parent must exist')
    started=time.process_time()
    commit,bindings=frozen_sources()
    inputs=[json.loads((ROOT/p).read_text(encoding='utf-8'))
            for p in (previous.BINDING,previous.LOCKED,previous.RESOLUTION)]
    branch,A,CW,Lref=previous.select_inputs(*inputs)
    saved=json.loads((ROOT/TT).read_text(encoding='utf-8'))
    if saved['status']!='completed_flat_Euclidean_profile_response':
        raise ValueError('the original TT comparison must be complete')
    for name,digest in saved['source_sha256'].items():
        if bindings.get(name)!=digest:
            raise ValueError('a saved TT source/input differs: '+name)
    for key,value in [('A',A),('C_W',CW),('L_ref',Lref)]:
        if not isclose(saved['parameters'][key],value,rel_tol=1e-13,abs_tol=0):
            raise ValueError('the saved TT reference differs: '+key)
    conv=ScaleBindingConventions(magnetic_flux=branch['magnetic_flux'])
    numbers=branch['numerical']
    if (numbers['compact_modes'],numbers['quadrature_points'])!=(34,112):
        raise ValueError('the admitted compact representation changed')
    cutoff=branch['omega']/conv.throat_radius
    matching=conv.magnetic_matching_scale
    q1=6*(4*pi)**2*A
    h0=float(saved['spectral_moments']['profile_at_zero'])
    x=1/Lref**2
    def profile(y):
        if time.process_time()-started>CPU_LIMIT:
            raise RuntimeError('the curvature response exhausted its CPU allowance')
        return matched_weight(weight,y,matching)
    with threadpool_limits(limits=1):
        weight=WarpedCompactWeight(cutoff=cutoff,
            interval=conv.interval_over_throat_radius*conv.throat_radius,
            modes=numbers['compact_modes'],points=numbers['quadrature_points'],
            amplitude=conv.warp_amplitude,path=conv.warp_path)
        cases=[]
        for n in (16,32):
            kernels=curvature_profile_kernels(x,profile,h0,quadrature_order=n)
            response=static_gravitational_response(x,kernels['weyl_kernel'],kernels['scalar_kernel'],q1)
            cases.append({'kernels':kernels,'response':response})
        pools=[{'internal_api':p['internal_api'],'num_threads':p['num_threads']} for p in threadpool_info()]
    gaps={key:abs(cases[1]['response'][key]-cases[0]['response'][key])
          for key in cases[1]['response']}
    if max(gaps.values())>1e-9:raise ValueError('the admitted static-response quadrature comparison failed')
    tt_gap=abs(cases[1]['response']['spin2_factor']-saved['quadrature_cases'][-1]['response_over_Einstein'])
    if tt_gap>1e-9:raise ValueError('covariant Weyl response disagrees with the closed TT response')
    elapsed=time.process_time()-started
    if elapsed>CPU_LIMIT:raise RuntimeError('the curvature response exhausted its CPU allowance')
    result={
        'schema':'nsc-covariant-spectral-response-v1',
        'status':'completed_weak_curvature_static_response',
        'source_commit':commit,'source_sha256':bindings,
        'parameters':{'A':A,'C_W':CW,'L_ref':Lref,'ultraviolet_cutoff':cutoff,
            'magnetic_matching_scale':matching,'compact_interval':weight.interval,
            'warp_amplitude':weight.amplitude,'warp_path':weight.path,
            'compact_modes':weight.modes,'overlap_points':weight.points,
            'couplings_or_profile_refitted':False},
        'reference':{'p2':x,'p2_Lref2':1.0,'target_observation_used':False,
            'physical_length_calibrated':False,
            'metric':'ds^2=(1+2Phi)dt^2-(1-2Psi)dx^2',
            'source':'conserved static external T00=rho at linear order; Tij=0',
            'Einstein_potential':'Phi_GR=Psi_GR=-rho/(4*A*p2)'},
        'spectral_moments':{'first_moment_reused_from_A':q1,'profile_at_zero_reused':h0},
        'cases':cases,'quadrature_response_gaps':gaps,'saved_TT_factor_gap':tt_gap,
        'zero_mode':{'Weyl_kernel':-h0/20,'scalar_kernel':0.0,'scalar_trace_factor':1.0,
            'added_volume_term':False},
        'scope':('Curvature-squared fixed-profile action in weak metric perturbations about R4, '
            'fixed topology and no surviving outer boundary. Internal compact data do not vary. '
            'Static conserved-source transfer ratios are conditional measurable relations, not '
            'astronomical fits or a physical length calibration. Full curved real-time/state '
            'response, material walls, strong curvature and a cosmological Lambda are not derived.'),
        'uncertainty':('The 16/32 gaps are empirical quadrature indicators. Inherited compact and '
            'moment errors, nonlinear fields and the other physical action/state contributions '
            'are not enclosed; no total observational error range is inferred.'),
        'counts':{'physical_momenta':1,'quadrature_rules':[16,32],
            'old_TT_producer_reruns':0,'root_solves':0,'global_moment_integrations':0,
            'compact_refinements':0,'trajectories':0,'parameter_fits':0},
        'runtime_cpu_seconds':elapsed,'threadpools':pools,
        'limits':{'CPU_seconds':CPU_LIMIT,'summary_bytes':BYTE_LIMIT},
        'reproduce':('python scripts/run_python.py scripts/measure_nsc_covariant_spectral_response.py '
                     '--output .build/covariant-spectral-response.json'),
    }
    raw=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(raw.encode('utf-8'))>BYTE_LIMIT:raise RuntimeError('the compact summary exceeds its allocation')
    args.output.write_text(raw,encoding='utf-8',newline='\n')
    print(json.dumps({'output':str(args.output),'CPU_seconds':elapsed,
                      'response':cases[-1]['response'],'quadrature_gaps':gaps,'saved_TT_gap':tt_gap}))


if __name__=='__main__':main()
