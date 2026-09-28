#!/usr/bin/env python3
"""Complete only the archived real rows rejected by packet inversion."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_pg_archived_high_energy_modes import archived_middle_boundary_modes
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT='results/development/nsc-pg-high-energy-inputs.json'
SOURCES=('src/recursive_horizons/nsc_pg_archived_high_energy_modes.py','tests/test_nsc_pg_archived_high_energy_modes.py',
         'scripts/derive_nsc_pg_high_energy_inputs.py','docs/nsc-pg-high-energy-inputs.md')
INPUTS=('results/development/nsc-pg-spectral-mode-recovery.json','results/development/nsc-pg-retained-covariance.json',
        'src/recursive_horizons/nsc_pg_fast_packets.py','src/recursive_horizons/nsc_pg_high_energy.py',
        'src/recursive_horizons/nsc_lorentzian.py','src/recursive_horizons/nsc_transmitting_dirac_domain.py',
        'scripts/derive_nsc_transmitting_boundary_binding.py')

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def load(p):
    with np.load(p,allow_pickle=False) as a:return {k:a[k].copy() for k in a.files}
def payload(record):
    p=record['payload']
    if sha(p['path'])!=p['sha256']:raise ValueError('input payload changed '+p['path'])
    return load(ROOT/p['path'])

def calculate():
    prior=read(INPUTS[0]); old=payload(prior); retained=payload(read(INPUTS[1]))
    arrays={};rows={};maxima={};failures={};filled=0
    for name,panel in prior['panels'].items():
        if panel['rejected']['count']==0:continue
        if not name.startswith(('mid/','mid_ref/')):raise ValueError('unowned missing real panel '+name)
        get=lambda k:old[name+'/'+k]
        E=get('energies');accepted=get('accepted');missing=np.flatnonzero(~accepted)
        if len(missing)!=panel['rejected']['count']:raise ValueError('original inverse mask changed')
        metadata=json.loads(retained[name+'/metadata_json'].tobytes())
        direct=archived_middle_boundary_modes(E,metadata,repo_root=ROOT)
        match=direct.compare_observations(get('observation'),get('projection'))
        overlap=np.linalg.norm(direct.mode_at_one[accepted]-get('mode_at_one')[accepted],axis=(1,2))
        residuals={'projection':float(match['projection_absolute_residual'].max()),
                   'accepted_inverse_overlap':float(overlap.max()) if len(overlap) else 0.,
                   'current_coisometry':float(direct.current_coisometry_residual.max()),
                   'riccati_equation':float(direct.riccati_equation_residual.max()),
                   'order12_16_difference':float(direct.order12_16_field_difference.max()),
                   'source_covariance_change':0.}
        tol={'projection':3e-11,'accepted_inverse_overlap':3e-8,'current_coisometry':3e-11,
             'riccati_equation':3e-11,'order12_16_difference':3e-11,'source_covariance_change':0.}
        failed={k:v for k,v in residuals.items() if not np.isfinite(v) or v>tol[k]}
        if failed:failures[name]=failed
        # Complete only missing indices; preserve every accepted inverse row.
        completed=get('mode_at_one').copy();completed[missing]=direct.mode_at_one[missing]
        unchanged=np.array_equal(completed[accepted],get('mode_at_one')[accepted])
        if not unchanged or not np.isfinite(completed).all():raise ValueError('completion altered accepted input or left missing rows')
        arrays[name+'/indices']=missing
        arrays[name+'/mode_at_one']=direct.mode_at_one[missing]
        arrays[name+'/energies']=E[missing]
        rows[name]={'new_rows':len(missing),'all_panel_rows':len(E),'accepted_overlap_rows':int(accepted.sum()),
                    'energy_range':[float(E[missing].min()),float(E[missing].max())],
                    'residuals':residuals,'tolerances':tol,'status':'PASS' if not failed else 'OPEN',
                    'accepted_inverse_unchanged':unchanged,'provenance':direct.provenance}
        for k,v in residuals.items():maxima[k]=max(maxima.get(k,0.),v)
        filled+=len(missing)
    if filled!=prior['counts']['unresolved_rows']:failures['coverage']='not all rejected real rows covered'
    result={'schema':'NSC-PG-HIGH-ENERGY-INPUTS-v1','accountable_author':'Douglas Ek',
      'status':'PASS: archived real-panel input completion; full spectral source remains OPEN' if not failures else 'OPEN',
      'source_hashes':{p:sha(p) for p in SOURCES},'input_hashes':{p:sha(p) for p in INPUTS},
      'locked_inputs':prior['locked_inputs'],'input_payloads':[prior['payload'],read(INPUTS[1])['payload']],
      'panels':rows,'maxima':maxima,'failures':failures,
      'counts':{'completed_rows':filled,'completed_panels':len(rows),'prior_inverse_rows':prior['counts']['accepted_rows'],
                'total_archived_real_rows':prior['counts']['rows'],'all_real_rows_available':not failures},
      'domain':{'operator':'same pinned massive PG Riccati expansion, order16, rho=1',
        'source_column_order':['zero','partner','incoming_one'],'phase':'exact archived middle-panel producer; no added PG clock phase',
        'source_covariance':'unchanged C_H and inherited incoming occupation; no occupation fitting or clipping',
        'approximation':'original vacuum-mode approximation; order comparison is not a rigorous stress-tail bound',
        'panel_accounting':'base/refinement panels kept separate; completion is row coverage, not a spectral sum',
        'old_inverse_certificate':'unchanged; its rejected rows remain rejected by that inverse',
        'horizon_scattering_rerun':False,'packet_integrals_rerun':False,'history_selected':False},
      'gate':{'archived_real_input_fields':'PASS' if not failures else 'OPEN',
        'full_spectral_source':'OPEN: subgap/infinite-tail differentiated source, common nonlinear subtraction and convergence',
        'physical_EndpointBranchJets':'OPEN','stress':None,'nulls':None,'updated_constraints':None,
        'V_c':None,'extended_stationarity':'OPEN','metric_timestep':False,'Z3':'OUT OF SCOPE','PDF_bumped':False,
        'Weyl_time_node_diagnostic':93.54264532195464,'homogeneous_nonexistence':'preserved in original scope; not rerun'},
      'comparison':{'fields':'all','float_atol':3e-13,'float_rtol':3e-13,'exact':'hashes, coverage, indices and scope'},
      'reproducer':'python3 scripts/derive_nsc_pg_high_energy_inputs.py --check'}
    return result,arrays

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--prepare',action='store_true');parser.add_argument('--check',action='store_true');args=parser.parse_args()
    result,arrays=calculate()
    if args.prepare:
        meta={k:result[k] for k in ('source_hashes','input_hashes','input_payloads')}
        arrays['metadata_json']=np.frombuffer(json.dumps(meta,sort_keys=True).encode(),np.uint8)
        raw=deterministic_npz_bytes(arrays);digest=hashlib.sha256(raw).hexdigest()
        path=f'results/development/artifacts/nsc-pg-high-energy-inputs.{digest}.npz';(ROOT/path).write_bytes(raw)
        result['payload']={'path':path,'sha256':digest,'bytes':len(raw)}
        (ROOT/OUTPUT).write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    elif args.check:
        from derive_nsc_transmitting_boundary_binding import compare
        expected=read(OUTPUT);saved=payload(expected)
        for k,v in arrays.items():
            if not np.array_equal(saved[k],v):raise ValueError('completion artifact changed '+k)
        meta=json.loads(saved['metadata_json'].tobytes())
        for k in ('source_hashes','input_hashes','input_payloads'):compare(result[k],meta[k])
        result['payload']=expected['payload'];compare(expected,result)
    else:parser.error('use --prepare or --check')
    print(json.dumps({k:result[k] for k in ('status','counts','maxima','failures')},indent=2),flush=True)
    if result['failures']:raise SystemExit(1)

if __name__=='__main__':main()
