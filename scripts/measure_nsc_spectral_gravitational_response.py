"""One fixed-scale flat gravitational evaluation of the owned compact profile."""
from __future__ import annotations

import argparse
import hashlib
import json
from math import isclose, isfinite, pi, sqrt
from pathlib import Path
import subprocess
import time

from threadpoolctl import threadpool_info, threadpool_limits

from recursive_horizons.dirac.compact import matched_weight, zero_limit
from recursive_horizons.nsc_scale_binding import ScaleBindingConventions
from recursive_horizons.nsc_spectral_gravitational_response import (
    flat_spin2_profile_response, profile_local_coefficients,
)
from recursive_horizons.nsc_warped_source import WarpedCompactWeight


ROOT = Path(__file__).resolve().parents[1]
BINDING = 'results/development/scale-binding.json'
LOCKED = 'results/development/nsc-subgap-history-response.json'
RESOLUTION = 'results/summaries/nsc-conformal-resolution-v1.json'
CPU_LIMIT = 120.0
BYTE_LIMIT = 12000
SOURCES = (
    'scripts/measure_nsc_spectral_gravitational_response.py',
    'src/recursive_horizons/nsc_spectral_gravitational_response.py',
    'tests/test_nsc_spectral_gravitational_response.py',
    'src/recursive_horizons/dirac/compact.py',
    'scripts/check_nsc_compact_matching.py',
    'src/recursive_horizons/nsc_warped_source.py',
    'src/recursive_horizons/nsc_scale_binding.py',
    BINDING, LOCKED, RESOLUTION,
)


def frozen_sources():
    commit = subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip()
    hashes = {}
    for name in SOURCES:
        raw = (ROOT/name).read_bytes()
        if subprocess.check_output(['git','show',f'{commit}:{name}'], cwd=ROOT) != raw:
            raise ValueError('freeze the profile method before application: '+name)
        hashes[name] = hashlib.sha256(raw).hexdigest()
    return commit, hashes


def select_inputs(binding, locked, resolution):
    branch = binding['development_branch']
    coefficients = branch['coefficients']
    A, CW = (float(coefficients[key]) for key in ('A_Wilsonian_at_magnetic_scale','C_Weyl'))
    if not (isfinite(A) and A > 0 and isfinite(CW) and CW < 0):
        raise ValueError('this application retains the positive-A, negative-C_W branch')
    for actual, key in ((A,'A'),(CW,'C_Weyl')):
        if not isclose(actual, float(locked['locked_inputs'][key]), rel_tol=1e-13, abs_tol=0):
            raise ValueError('the locked local and spectral coefficients differ: '+key)
    if not isclose(float(branch['omega']), float(locked['locked_inputs']['Omega']), rel_tol=1e-13):
        raise ValueError('the locked spectral cutoff differs')
    reference = resolution['domain']
    Lref = float(reference['reference_length'])
    if (reference['target_observation_used'] is not False
            or not isfinite(Lref) or Lref <= 0):
        raise ValueError('an independently fixed observer ruler is required')
    for actual, key in ((A,'A'),(CW,'C_W')):
        if not isclose(actual,float(resolution['parameters'][key]),rel_tol=1e-13,abs_tol=0):
            raise ValueError('the observer reference uses different couplings: '+key)
    if (int(branch['magnetic_flux']) != int(locked['locked_inputs']['magnetic_flux'])
            or branch['magnetic_flux'] != 4
            or coefficients['V_full_relational'] != 0
            or not isclose(float(branch['mmp_map']['r_e_squared']),1,abs_tol=1e-12)):
        raise ValueError('the fixed-flux throat-unit spectral branch changed')
    return branch, A, CW, Lref


def evaluate(branch, A, CW, Lref, started):
    conventions = ScaleBindingConventions(magnetic_flux=branch['magnetic_flux'])
    numerical = branch['numerical']
    if (numerical['compact_modes'],numerical['quadrature_points']) != (34,112):
        raise ValueError('the admitted compact representation changed')
    cutoff = float(branch['omega'])/conventions.throat_radius
    matching = conventions.magnetic_matching_scale
    weight = WarpedCompactWeight(cutoff=cutoff,
        interval=conventions.interval_over_throat_radius*conventions.throat_radius,
        modes=numerical['compact_modes'], points=numerical['quadrature_points'],
        amplitude=conventions.warp_amplitude, path=conventions.warp_path)
    q1 = 6*(4*pi)**2*A
    q0 = zero_limit(weight,matching)['matched_weight_at_zero']
    recovered = profile_local_coefficients(q0,q1)
    if not isclose(recovered['C_W'],CW,rel_tol=1e-10,abs_tol=1e-13):
        raise ValueError('the owned profile does not recover the locked Weyl coefficient')
    def profile(y):
        if time.process_time()-started > CPU_LIMIT:
            raise RuntimeError('the fixed profile application exhausted its CPU allowance')
        return matched_weight(weight,y,matching)
    x = 1/Lref**2
    cases = [flat_spin2_profile_response(x,profile,q1,quadrature_order=n) for n in (16,32)]
    gap = abs(cases[1]['response_over_Einstein']-cases[0]['response_over_Einstein'])
    if gap > 1e-9:
        raise ValueError('the predeclared Jacobi comparison failed')
    local_ratio = 1+2*CW/A*x
    extra_root = -A/(2*CW)
    return {
        'parameters':{'A':A,'C_W':CW,'L_ref':Lref,'ultraviolet_cutoff':cutoff,
            'magnetic_matching_scale':matching,'compact_interval':weight.interval,
            'warp_amplitude':weight.amplitude,'warp_path':weight.path,
            'compact_modes':weight.modes,'overlap_points':weight.points,
            'couplings_or_profile_refitted':False},
        'spectral_moments':{'first_complement_moment_reused_from_A':q1,
            'profile_at_zero':q0,'recovered_C_W_minus_locked':recovered['C_W']-CW},
        'scale_comparison':{'p2':x,'p2_Lref2':1.0,'p_over_UV_cutoff':sqrt(x)/cutoff,
            'p2_over_matching_scale2':x/matching**2,
            'local_extra_root_p2':extra_root,
            'local_extra_root_over_cutoff2':extra_root/cutoff**2,
            'local_extra_momentum_over_cutoff':sqrt(extra_root)/cutoff},
        'quadrature_cases':cases,
        'comparison':{'dimensionless_quadrature_gap':gap,
            'local_quartic_response_over_Einstein':local_ratio,
            'full_profile_minus_local_ratio':cases[-1]['response_over_Einstein']-local_ratio},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    args = parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    if not args.output.parent.is_dir():raise ValueError('output parent must exist')
    started = time.process_time()
    commit, hashes = frozen_sources()
    values = [json.loads((ROOT/p).read_text(encoding='utf-8')) for p in (BINDING,LOCKED,RESOLUTION)]
    branch,A,CW,Lref = select_inputs(*values)
    with threadpool_limits(limits=1):
        measured = evaluate(branch,A,CW,Lref,started)
        pools = [{'internal_api':p['internal_api'],'num_threads':p['num_threads']} for p in threadpool_info()]
    elapsed = time.process_time()-started
    if elapsed > CPU_LIMIT:raise RuntimeError('the profile application exhausted its CPU allowance')
    result = {
        'schema':'nsc-spectral-gravitational-response-v1',
        'status':'completed_flat_Euclidean_profile_response',
        'source_commit':commit,'source_sha256':hashes,**measured,
        'definition':('Gamma=(1/2) Tr H(D4^2), H=h_Lambda-E1(u/mu_B^2); compact data fixed. '
            'The existing zero-volume-term law removes the constant TT term. '
            'S_TT^(2)=(A/4) integral h_TT P_H(-partial^2) h_TT.'),
        'scope':('A flat Euclidean transverse-traceless response of the already owned induced '
            'spectral profile. The additional local-polynomial root is compared with the '
            'separately owned cutoff; it is not a root measurement of the full kernel. '
            'This is not the spherical retarded map, a curved stability proof, a matter-state '
            'completion, an exterior attachment, a fractal dimension or cosmological Lambda.'),
        'uncertainty':('The fixed-rule gap is an empirical integration indicator. The inherited '
            'compact Galerkin and first-moment errors are not enclosed here. Sample profile '
            'extrema are not global bounds; no total physical error bar is inferred.'),
        'counts':{'physical_momenta':1,'quadrature_rules':[16,32],'compact_refinements':0,
            'root_solves':0,'global_moment_reintegrations':0,'trajectories':0,'parameter_fits':0},
        'runtime_cpu_seconds':elapsed,'threadpools':pools,
        'limits':{'CPU_seconds':CPU_LIMIT,'summary_bytes':BYTE_LIMIT},
        'reproduce':('python scripts/run_python.py scripts/measure_nsc_spectral_gravitational_response.py '
                     '--output .build/spectral-gravity-profile.json'),
    }
    raw = json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(raw.encode('utf-8')) > BYTE_LIMIT:raise RuntimeError('summary exceeds admitted size')
    args.output.write_text(raw,encoding='utf-8',newline='\n')
    print(json.dumps({'output':str(args.output),'CPU_seconds':elapsed,**measured['comparison']}))


if __name__ == '__main__':main()
