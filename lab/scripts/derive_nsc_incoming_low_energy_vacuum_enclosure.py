#!/usr/bin/env python3
"""One directed low-energy projector run, or propagation-free proof replay."""
import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_low_energy_vacuum_enclosure import ENERGY,STEPS,PRECISION,propagate,replay
OUTPUT='results/development/nsc-incoming-low-energy-vacuum-enclosure.json'
SOURCES=('src/recursive_horizons/nsc_incoming_low_energy_vacuum_enclosure.py',
    'scripts/derive_nsc_incoming_low_energy_vacuum_enclosure.py','tests/test_nsc_incoming_low_energy_vacuum_enclosure.py',
    'docs/nsc-incoming-low-energy-vacuum-enclosure.md')
INPUTS=('results/development/nsc-mode-resolved-cauchy-state.json','results/development/nsc-compact-matched-restart.json',
    'results/development/nsc-pg-ctp-mode-jets.json','src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
    'src/recursive_horizons/nsc_incoming_validated_vacuum_propagation.py','src/recursive_horizons/nsc_incoming_source_quadrature_bound.py',
    'src/recursive_horizons/nsc_paired_horizon_preparation.py','src/recursive_horizons/nsc_pg_massive_modes.py',
    'src/recursive_horizons/nsc_unruh_state.py','src/recursive_horizons/nsc_lorentzian.py')
CPU_CAP=120.


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()
def signature():return {p:digest(p) for p in (*SOURCES,*INPUTS)}
def inputs():
    channel=json.loads((ROOT/INPUTS[0]).read_text())['channels'][14]
    config=json.loads((ROOT/INPUTS[1]).read_text())['scattering_provenance']['config']
    reference=json.loads((ROOT/INPUTS[2]).read_text());spec=reference['payload']
    if digest(spec['path'])!=spec['sha256']:raise ValueError('selected source-mode artifact changed')
    with np.load(ROOT/spec['path'],allow_pickle=False) as a:
        if (float(a['14_1/energies'][3])!=ENERGY or a['14_1/mass'].item()!=channel['compact_mass']
            or a['14_1/angular'].item()!=channel['angular_eigenvalue']):raise ValueError('selected energy/channel changed')
    return channel,config,spec


def make_record(payload,path):
    if payload['signature']!=signature():raise ValueError('projector proof/input source changed')
    channel,config,spec=inputs()
    if payload['channel']!=channel or payload['source_config']!=config or payload['selected_source_artifact']!=spec:
        raise ValueError('physical/source parameter binding changed')
    if payload.get('budget_exceeded'):
        result=None;status='OPEN: fixed CPU budget exhausted before certificate'
    else:
        result=replay(payload)
        if result!=payload['result']:raise ValueError('directed projector error replay differs')
        status='PASS: directed physical affine-vacuum Re(P01)>1/4' if result['strict_Re_P01_above_one_quarter'] else 'OPEN: fixed directed enclosure does not establish Re(P01)>1/4'
    return {'schema':'NSC-INCOMING-LOW-ENERGY-VACUUM-ENCLOSURE-v1','accountable_author':'Douglas Ek',
        'status':status,'group':14,'angular_sign':1,'energy':ENERGY,'channel':channel,'source_config':config,
        'result':result,'geometry':payload.get('geometry'),
        'initialization':payload.get('tail'),'endpoint_bridge':payload.get('endpoint_bridge'),
        'algorithm':{'steps':STEPS,'precision':PRECISION,'point_column_initial':[1,0],
            'state_wrapping':False,'column_normalization':False,
            'error_recurrence':'epsilon_next=epsilon+step*sup(norm(H-Hmid))*norm(point_column)+recenter_radius',
            'projector_error':'affine_tail+(1+norm(point_column))*column_error+endpoint_bridge',
            'replay':'directed sum, point chain and coordinate/error accounting; no repeat geometry propagation'},
        'runtime':payload['runtime'],'selected_source_artifact':spec,
        'scope':{'physical_affine_dominant_vacuum_component':True,'coherent_source_replaced':False,
            'coherent_thermal_remainder':'separate exact owner; not discarded or included in this vacuum projector error',
            'archived_modes_or_state_normalized_or_changed':False,'Jost_phase_changed':False,
            'source_occupation_kappa_changed':False,'group32_bound_applied_to_low_energy':False,
            'source_precision_3e_minus_11_claimed':False,'physical_initial_data_or_metric_evolution':False,
            'full_C0_matching_or_extended_stationarity':'OPEN','publication':False},
        'source_hashes':{p:payload['signature'][p] for p in SOURCES},'input_hashes':{p:payload['signature'][p] for p in INPUTS},
        'payload':{'path':str(path.relative_to(ROOT)),'sha256':sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size},
        'reproducer':'python3 scripts/derive_nsc_incoming_low_energy_vacuum_enclosure.py --check'}


class CPUBudget(RuntimeError):pass

def run():
    before=signature();channel,config,spec=inputs();cpu=time.process_time();wall=time.monotonic();completed=0
    def progress(count):
        nonlocal completed
        completed=count;print(f'{count}/{STEPS} directed midpoint cells; CPU {time.process_time()-cpu:.2f}s',flush=True)
    def stop(signum,frame):raise CPUBudget('fixed120CPU-second budget exhausted')
    prior=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CPU_CAP)
    try:payload=propagate(channel,config,progress=progress)
    except CPUBudget:
        payload={'budget_exceeded':True,'steps_completed_at_last_checkpoint':completed,'channel':channel,'source_config':config}
    finally:
        signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,prior)
    if signature()!=before:raise ValueError('proof/input source changed during production')
    payload.update(signature=before,selected_source_artifact=spec,
        runtime={'CPU_seconds':time.process_time()-cpu,'wall_seconds':time.monotonic()-wall,
                 'CPU_cap_seconds':CPU_CAP,'budget_exceeded':bool(payload.get('budget_exceeded'))})
    raw=json.dumps(payload,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
    compressed=gzip.compress(raw,compresslevel=6,mtime=0);h=sha256(compressed).hexdigest()
    path=ROOT/f'results/development/artifacts/nsc-incoming-low-energy-vacuum-enclosure.{h}.json.gz'
    path.write_bytes(compressed);record=make_record(payload,path)
    (ROOT/OUTPUT).write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return record


def check():
    old=json.loads((ROOT/OUTPUT).read_text());spec=old['payload'];path=ROOT/spec['path']
    if digest(spec['path'])!=spec['sha256']:raise ValueError('directed proof artifact changed')
    payload=json.loads(gzip.decompress(path.read_bytes()))
    record=make_record(payload,path)
    if record!=old:raise ValueError('directed projector receipt differs')
    return record


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true');args=p.parse_args()
    record=run() if args.run else check()
    print(json.dumps({'status':record['status'],'result':record['result'],'runtime':record['runtime']},indent=2))

if __name__=='__main__':main()
