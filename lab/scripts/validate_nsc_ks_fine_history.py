#!/usr/bin/env python3
"""Checkpointed continuous residual validation of the saved fine trajectory."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

import numpy as np
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import derive_nsc_ks_fine_trajectory as F
import derive_nsc_ks_fine_profile_table as P
import derive_nsc_ks_profile_fourier_bound as OLD
from recursive_horizons.nsc_ks_streaming_trajectory import TrajectoryChunkReader
from recursive_horizons.nsc_ks_ball_trajectory import BallFourierSegment,exact_upper,restored_upper
from recursive_horizons.nsc_ks_ball_operator import AnalyticRadiusFamily,time_operator_enclosure
from recursive_horizons.nsc_ks_ball_geometry import SubdivisionNeeded
from recursive_horizons.nsc_ks_residual_polynomial import ResidualPolynomial
from recursive_horizons.nsc_ks_fourier_residual_bound import polynomial_residual_bounds,operator_remainder_bounds,profile_sup_bounds
from recursive_horizons.nsc_ks_profile_compression import compress_profiles
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_characteristic_error import propagate_local_characteristics

DIRECTORY=ROOT/'results/development/artifacts/nsc-ks-fine-validation'
OUTPUT=ROOT/'results/development/nsc-ks-fine-validation.json'
OWNED=('scripts/validate_nsc_ks_fine_history.py','docs/nsc-ks-fine-validation.md')
EXTRA=('src/recursive_horizons/nsc_ks_profile_compression.py','src/recursive_horizons/nsc_ks_residual_polynomial.py',
       'src/recursive_horizons/nsc_ks_fourier_residual_bound.py','src/recursive_horizons/nsc_ks_characteristic_error.py',
       'results/development/nsc-ks-radius-bound.json')
INPUTS=tuple(dict.fromkeys((*F.OWNED,*F.INPUTS,str(F.OUTPUT.relative_to(ROOT)),
                         *P.OWNED,*P.INPUTS,str(P.OUTPUT.relative_to(ROOT)),*EXTRA)))
BITS,WORK_INDEX,TIME_DEGREE=120,1024,8


def digest(path):
    path=Path(path);return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()
def signatures():return {p:digest(p) for p in (*OWNED,*INPUTS)}
def write_json(path,data):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n')
    os.replace(temporary,path)


class Context:
    def __init__(self):
        self.signature=signatures()
        self.signature_digest=sha256(json.dumps(self.signature,sort_keys=True).encode()).hexdigest()
        self.fine=json.loads(F.OUTPUT.read_text());self.profile_record=P.check()
        if not self.fine['complete_capture']:raise ValueError('fine trajectory incomplete')
        for path,expected in {**self.fine['source_hashes'],**self.fine['input_hashes']}.items():
            if digest(path)!=expected:raise ValueError('fine trajectory owner/input changed')
        for item in (self.fine['payload'],self.fine['trajectory_manifest'],self.profile_record['payload']):
            if digest(item['path'])!=item['sha256']:raise ValueError('validation input payload changed')
        with np.load(F.FINAL,allow_pickle=False) as f:self.arrays={k:f[k] for k in f.files}
        self.meta=json.loads(self.arrays['metadata_json'].tobytes())
        self.family=F.T.C.R.C.P.family(F.T.C.R.C.P.ALPHA)
        with np.load(P.PAYLOAD,allow_pickle=False) as f:b={k:f[k] for k in f.files}
        pm=json.loads(b['metadata_json'].tobytes())
        if (profile_identity(self.family)!=self.fine['history_identity'] or
            profile_identity(self.family,include_normal_window=False)!=pm['axial_profile_identity'] or
            pm['axial_profile_identity']!=self.fine['axial_profile_identity'] or
            pm['source_preparation_digest']!=self.fine['source_preparation_digest']):
            raise ValueError('source/history/profile identity mismatch')
        grid=self.arrays['computational_z'];self.N=len(grid)
        if pm['period_origin']!=float(grid[0]) or pm['period_length']!=float((grid[1]-grid[0])*len(grid)):
            raise ValueError('profile table and trajectory use different numerical periods')
        self.length=pm['period_length'];self.reader=TrajectoryChunkReader(F.DIRECTORY)
        self.total=self.reader.manifest['accepted_steps'];self._chunk=None
        self.model=AnalyticRadiusFamily(self.family)
        self.radius=json.loads((ROOT/'results/development/nsc-ks-radius-bound.json').read_text())
        for path,expected in {**self.radius['source_hashes'],**self.radius['input_hashes']}.items():
            if digest(path)!=expected:raise ValueError('radius remainder owner/input changed')
        if self.radius['control']['prepared_source_digest']!=self.fine['source_preparation_digest']:
            raise ValueError('radius remainder has another source')
        with ctx.workprec(BITS):
            profiles={tuple(v['powers']):{'coefficients':[OLD.unpack_ball(c) for c in b[k+'/coefficients']],
                'tail':[restored_upper(t) for t in v['tail_upper']]} for k,v in pm['profiles'].items()}
            self.profiles=compress_profiles(profiles,self.length,WORK_INDEX,bits=BITS)
            self.profile_sup=profile_sup_bounds(self.profiles,self.length)
            self.radius_tail=[arb(v['exact_rational']) for v in self.radius['bounds']['potential_remainder_sup_orders_0_1_2']]

    def segment(self,index):
        if not 0<=index<self.total:raise ValueError('segment index outside the saved trajectory')
        offset=0
        for entry in self.reader.manifest['chunks']:
            if index<offset+entry['segment_count']:
                if self._chunk is None or self._chunk.index!=entry['index']:
                    self._chunk=self.reader.load_chunk(entry['index'])
                return self._chunk.segment(index-offset)
            offset+=entry['segment_count']
        raise ValueError('segment inventory is inconsistent')

    def evaluate_field(self,field,depth=0):
        try:
            polys,errors=time_operator_enclosure(self.model,field.rho_start,field.rho_end,
                self.meta['parameters']['angular'],degree=TIME_DEGREE,bits=BITS)
            extra=operator_remainder_bounds(field,errors,self.profile_sup,self.radius_tail,
                self.meta['parameters']['mass'],self.meta['parameters']['angular'],self.arrays['source/energies'])
            if extra[1]>arb('1e-19') and depth<8:
                raise SubdivisionNeeded('time coefficient remainder needs a smaller proof cell')
        except SubdivisionNeeded:
            if depth>=12:raise
            halves=[self.evaluate_field(field.restricted(a,b),depth+1) for a,b in ((arb(0),arb(1)/2),(arb(1)/2,arb(1)))]
            return tuple((halves[0][0][j]+halves[1][0][j]).upper() for j in range(2)),sum(v[1] for v in halves)
        residual=ResidualPolynomial(field,polys['inv_a2'],polys['inv_a'],polys['inv_ar'],
            {k:v for k,v in polys.items() if isinstance(k,tuple)},self.meta['parameters']['mass'],
            self.meta['parameters']['angular'],self.arrays['source/energies'])
        polynomial=polynomial_residual_bounds(residual,self.profiles)
        return tuple((a+b).upper() for a,b in zip(polynomial,extra)),1

    def evaluate(self,index):
        segment=self.segment(index);start=time.process_time()
        with ctx.workprec(BITS):
            field=BallFourierSegment(segment,self.arrays['source/column_weights'],self.N,self.length,bits=BITS)
            bound,cells=self.evaluate_field(field)
            return {'schema':'NSC-KS-FINE-RESIDUAL-SEGMENT-v1','segment':index,
                'rho_start':segment.rho_start,'rho_end':segment.rho_end,
                'continuous_residual_integral_upper':[exact_upper(v) for v in bound],
                'proof_subcells':cells,'CPU_seconds':time.process_time()-start,
                'signature_digest':self.signature_digest,'whole_spatial_period':True,
                'physical_local_gate':'OPEN'}


def saved_rows(owner):
    result={}
    for path in sorted(DIRECTORY.glob('segment-*.json')):
        row=json.loads(path.read_text());index=row['segment']
        if path.name!=f'segment-{index:06d}.json' or not 0<=index<owner.total or index in result:
            raise ValueError('invalid or duplicate segment proof')
        if row['signature_digest']!=owner.signature_digest:raise ValueError('saved proof owner/input changed')
        segment=owner.segment(index)
        if row['rho_start']!=segment.rho_start or row['rho_end']!=segment.rho_end:
            raise ValueError('saved proof covers another rho cell')
        result[index]=row
    return result


def summary(owner,rows):
    complete=len(rows)==owner.total
    with ctx.workprec(BITS):
        total=[sum((restored_upper(rows[i]['continuous_residual_integral_upper'][j]) for i in sorted(rows)),arb(0)).upper()
               for j in range(2)]
        errors=None
        if complete:
            r=arb(owner.radius['bounds']['radius_lower']['exact_rational'])
            K0=(arb(owner.meta['rho_up'])-1)/(arb(4)/5)*(arb(owner.meta['parameters']['mass'])+abs(arb(owner.meta['parameters']['angular']))/r)
            K1=arb(owner.radius['bounds']['B_z_integral_upper']['exact_rational'])
            errors=propagate_local_characteristics((0,0),total,K0.upper(),K1.upper(),
                float(np.max(abs(owner.arrays['source/energies']))),bits=BITS)
        return {'schema':'NSC-KS-FINE-HISTORY-VALIDATION-v1','accountable_author':'Douglas Ek',
            'status':('OPEN: numerical field bound completed; physical source/gate remain open' if complete else
                      'OPEN: partial continuous residual validation'),
            'validated_segments':len(rows),'required_segments':owner.total,'all_segments_covered':complete,
            'continuous_residual_integral_upper':[exact_upper(v) for v in total],
            'numerical_field_error_relative_to_saved_upstream':errors,
            'source_error_bound':None,'physical_local_gate':'OPEN','metric_timestep':False,
            'state_law':'C_Sigma[g]=U_g C_up U_gdagger; fixed upstream source',
            'signature':owner.signature,'signature_digest':owner.signature_digest,
            'method':{'precision_bits':BITS,'retained_profile_index':WORK_INDEX,'time_Taylor_degree':TIME_DEGREE},
            'proof_CPU_seconds':sum(rows[i]['CPU_seconds'] for i in sorted(rows)),
            'segment_proofs':[{'path':str((DIRECTORY/f'segment-{i:06d}.json').relative_to(ROOT)),
                              'sha256':digest(DIRECTORY/f'segment-{i:06d}.json')} for i in sorted(rows)]}


def main():
    p=argparse.ArgumentParser();group=p.add_mutually_exclusive_group(required=True)
    group.add_argument('--run',action='store_true');group.add_argument('--pilot',action='store_true');group.add_argument('--check',action='store_true')
    p.add_argument('--cpu-budget',type=float,default=900.)
    args=p.parse_args()
    if not 0<args.cpu_budget<=12000:raise ValueError('explicit CPU budget in (0,12000] required')
    DIRECTORY.mkdir(parents=True,exist_ok=True)
    owner=Context();rows=saved_rows(owner)
    if not args.check:
        selected=[122] if args.pilot else list(range(owner.total))
        start=time.process_time()
        for index in selected:
            if index in rows:continue
            if time.process_time()-start>=args.cpu_budget:break
            row=owner.evaluate(index)
            if owner.signature!=signatures():raise ValueError('proof source/input changed during validation')
            write_json(DIRECTORY/f'segment-{index:06d}.json',row);rows[index]=row
            progress=summary(owner,rows);write_json(DIRECTORY/'progress.json',progress)
            with ctx.workprec(BITS):display=[float(restored_upper(v)) for v in row['continuous_residual_integral_upper']]
            print(json.dumps({'validated':len(rows),'total':owner.total,'segment':index,'bounds':display,
                              'CPU_seconds':row['CPU_seconds']}),flush=True)
    result=summary(owner,rows)
    if result['all_segments_covered']:
        if OUTPUT.exists() and json.loads(OUTPUT.read_text())!=result:raise ValueError('completed validation record changed')
        if not OUTPUT.exists():write_json(OUTPUT,result)
    else:write_json(DIRECTORY/'progress.json',result)
    print(json.dumps({k:result[k] for k in ('status','validated_segments','required_segments','proof_CPU_seconds')},indent=2))


if __name__=='__main__':main()
