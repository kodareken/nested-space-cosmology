#!/usr/bin/env python3
"""Authenticate and compose all already-computed compact subgap panels."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from recursive_horizons.nsc_pg_packet_modes import SIGNED_PACKET_MAP
from derive_nsc_transmitting_boundary_binding import compare

OUTPUT='results/development/nsc-pg-retained-subgap.json'
PAYLOAD='results/development/artifacts/nsc-pg-retained-subgap.f7f6c6458983435b409704c4aab75c53af61e7718bdf9d97cf1f91bbe4a81d26.npz'
SOURCES=('scripts/check_nsc_pg_retained_subgap.py','scripts/derive_nsc_pg_retarded_threshold.py',
 'src/recursive_horizons/nsc_pg_retarded_threshold.py','src/recursive_horizons/nsc_pg_adaptive_threshold.py',
 'src/recursive_horizons/nsc_pg_retarded_packets.py','docs/nsc-pg-retained-subgap.md')

def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()

def calculate():
    digest=sha(PAYLOAD)
    if digest not in Path(PAYLOAD).name:raise ValueError('subgap payload digest differs')
    with np.load(ROOT/PAYLOAD,allow_pickle=False) as f:a={k:f[k].copy() for k in f.files}
    metadata=json.loads(a['metadata_json'].tobytes())
    for path,expected in metadata['source_hashes'].items():
        if sha(path)!=expected:raise ValueError('subgap producer changed: '+path)
    state=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())
    expected=[]
    for c in state['channels']:
        if c['compact_level']:
            expected.extend((c['index'],s) for s in ((1,-1) if c['angular_eigenvalue'] else (1,)))
    labels=[tuple(map(int,row)) for row in a['labels']]
    if labels!=sorted(expected):raise ValueError('missing or invented signed family')
    lookup={key:i for i,key in enumerate(labels)};S=SIGNED_PACKET_MAP;tol=3e-9
    maxima={'signed_assembly':0.,'Hermiticity':0.,'occupied_lower_violation':0.,'complement_lower_violation':0.,'quadrature_estimate':0.}
    rows=[]
    for i,key in enumerate(labels):
        j=lookup.get((key[0],-key[1]),i);G=a['gram_positive'][i];K=a['centered_positive'][i]
        if not np.isfinite(G).all() or not np.isfinite(K).all():raise ValueError('nonfinite subgap field')
        gs=(G+S@a['gram_positive'][j].conj()@S)/(2*np.pi)
        ks=(K-S@a['centered_positive'][j].conj()@S)/(2*np.pi)
        residual=float(max(np.linalg.norm(gs-a['gram_signed'][i]),np.linalg.norm(ks-a['centered_signed'][i])))
        C=K+.5*G
        low=max(0.,-float(np.linalg.eigvalsh(C).min()));high=max(0.,-float(np.linalg.eigvalsh(G-C).min()))
        receipt=next(r for r in metadata['receipts'] if (r['channel'],r['sign'])==key)
        if any(not d['success'] for d in receipt['diagnostics']):raise ArithmeticError('adaptive quadrature not successful')
        error=sum(d['estimated_error'] for d in receipt['diagnostics'])
        values={'signed_assembly':residual,'Hermiticity':float(max(np.linalg.norm(G-G.conj().T),np.linalg.norm(K-K.conj().T))),
                'occupied_lower_violation':low,'complement_lower_violation':high,'quadrature_estimate':error}
        for name,value in values.items():maxima[name]=max(maxima[name],value)
        rows.append({'channel':key[0],'angular_sign':key[1],'residuals':values,'panel_sha256':receipt['sha256']})
    if max(maxima.values())>tol:raise ArithmeticError('retained subgap gate failed')
    locked=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['locked_inputs']
    return {'schema':'NSC-PG-RETAINED-SUBGAP-v1','status':'PASS: all 38 declared compact signed-family subgap panels; full-energy covariance remains separate',
      'source_hashes':{p:sha(p) for p in SOURCES},
      'input_hashes':{p:sha(p) for p in ('results/development/nsc-mode-resolved-cauchy-state.json','results/development/nsc-compact-matched-restart.json',
                      'results/development/nsc-pg-massive-mode-resolution.json','results/development/nsc-pg-group13-covariance.json')},
      'payload':{'path':PAYLOAD,'sha256':digest,'bytes':(ROOT/PAYLOAD).stat().st_size},'locked_inputs':locked,
      'groups':20,'signed_families':38,'families':rows,'maximum_residuals':maxima,'tolerance':tol,
      'method':{'domain':'same whole-line PG Dirac operator, evaluated on the existing eight packet observables',
                'energy_panel':'[1,m_j] only, m_j=j*pi/2 for the declared compact groups',
                'state':'same C_H and inherited incoming state; closed infinity channel below parent mass',
                'assembly':'explicit opposite-angular positive/negative-energy pairing; no X-preserving completion',
                'normalization':'bounded retarded Green spectral jump plus the retained horizon coherence',
                'integration_error':'adaptive estimates; fixed quadrature/height controls retained for the initial groups',
                'no_new_physical_terms':True,'seed_covariance_used_as_input':False},
      'gate':{'subgap_panels':'PASS','full_C1b':'OPEN','metric_timestep':False,'stress_computed':False,'history_selected':False,'PDF_bump':False},
      'comparison':{'fields':'all','float_atol':3e-13,'float_rtol':3e-13,'exact':'structure and authenticated hashes',
                    'routine':'recompose archived physical subgap matrices; no old generator rerun'}}

def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();record=calculate()
    if a.check:compare(json.loads((ROOT/OUTPUT).read_text()),record)
    else:(ROOT/OUTPUT).write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':record['status'],'residuals':record['maximum_residuals'],'full_C1b':'OPEN'},indent=2))

if __name__=='__main__':main()
