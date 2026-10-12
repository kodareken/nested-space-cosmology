#!/usr/bin/env python3
"""Replay the same continuous residual proof with certified profile reduction."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from types import FunctionType

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from flint import ctx
import derive_nsc_ks_continuous_residual_cell as ORIGINAL
from recursive_horizons.nsc_ks_profile_compression import compress_profiles
from recursive_horizons.nsc_ks_fourier_residual_bound import polynomial_residual_bounds
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper

OUTPUT=ROOT/'results/development/nsc-ks-compressed-residual-cell.json'
OWNED=('scripts/derive_nsc_ks_compressed_residual_cell.py','src/recursive_horizons/nsc_ks_profile_compression.py',
       'docs/nsc-ks-profile-compression.md')
RETAINED=1024


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()


def compute():
    previous=ORIGINAL.check()
    def reduced(residual,profiles):
        compressed=compress_profiles(profiles,residual.field.length,RETAINED,bits=residual.bits)
        return polynomial_residual_bounds(residual,compressed)
    original=ORIGINAL.compute
    runner=FunctionType(original.__code__,{**original.__globals__,'polynomial_residual_bounds':reduced},
                        original.__name__,original.__defaults__,original.__closure__)
    runner.__kwdefaults__=dict(original.__kwdefaults__ or {})
    result=runner()
    result['schema']='NSC-KS-COMPRESSED-RESIDUAL-CELL-v1'
    result['status']='PASS: same continuous cell with bounded profile reduction; physical local gate OPEN'
    result['source_hashes'].update({p:digest(p) for p in OWNED})
    relative=str(ORIGINAL.OUTPUT.relative_to(ROOT))
    result['input_hashes'][relative]=digest(relative)
    result['method']={'retained_profile_index':RETAINED,'discarded_coefficients_added_to_derivative_tails':True,
                      'physical_profile_changed':False,'process_globals_patched':False}
    result['previous_bounds']=previous['bounds']
    result['previous_CPU_seconds']=previous['runtime']['CPU_seconds']
    result['reproducer']='python scripts/derive_nsc_ks_compressed_residual_cell.py --check'
    return result


def check():
    saved=json.loads(OUTPUT.read_text())
    for path,expected in {**saved['source_hashes'],**saved['input_hashes']}.items():
        if digest(path)!=expected:raise ValueError('compressed residual source/input changed')
    for item in saved['reused_artifacts']:
        if digest(item['path'])!=item['sha256']:raise ValueError('compressed residual data changed')
    with ctx.workprec(ORIGINAL.BITS):
        for a,b,total in zip(saved['bounds']['polynomial_and_profile_tail'],
                            saved['bounds']['time_coefficient_and_radius_remainder'],
                            saved['bounds']['total_continuous_normalized_residual']):
            if not restored_upper(total)>=(restored_upper(a)+restored_upper(b)).upper():
                raise ValueError('compressed residual error composition failed')
    if saved['method']['physical_profile_changed'] or saved['coverage']['all_history_segments']:
        raise ValueError('compression does not change profile or complete the history gate')
    return saved


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true')
    if p.parse_args().run:
        if OUTPUT.exists():raise FileExistsError('compressed residual cell exists; use --check')
        result=compute();OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    else:result=check()
    print(json.dumps({'status':result['status'],'continuous_residual_upper_display':
        [float(restored_upper(v)) for v in result['bounds']['total_continuous_normalized_residual']],
        'runtime':result['runtime'],'previous_CPU_seconds':result['previous_CPU_seconds']},indent=2))
