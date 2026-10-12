#!/usr/bin/env python3
"""Record the conditional incoming jet domain and its bounded UV screen."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tests'))
from recursive_horizons.nsc_incoming_cauchy_jets import incoming_reference_uv_difference,NORMAL_ENTRIES
from test_nsc_incoming_cauchy_jets import derivative_probe

OUTPUT='results/development/nsc-incoming-cauchy-data.json'
SOURCES=('src/recursive_horizons/nsc_incoming_cauchy_jets.py','tests/test_nsc_incoming_cauchy_jets.py',
         'scripts/derive_nsc_incoming_cauchy_data.py','docs/nsc-incoming-cauchy-data.md')
INPUTS=('results/development/nsc-incoming-constraint-gate.json','results/development/nsc-mode-resolved-cauchy-state.json',
        'src/recursive_horizons/nsc_spatial_reference_symbol.py','src/recursive_horizons/nsc_reference_band_action.py',
        'src/recursive_horizons/nsc_common_subtracted_ks_source.py','src/recursive_horizons/nsc_lll_geometric_history.py',
        'scripts/derive_nsc_transmitting_boundary_binding.py')

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def make_record():
    table=read(INPUTS[1])['channels'];indices=(0,1,13,14)
    cases=[{'name':f'group{j}','mass':table[j]['compact_mass'],'angular':table[j]['angular_eigenvalue']} for j in indices]
    domain=derivative_probe();row=incoming_reference_uv_difference(domain,cases,[16.,32.,64.,128.])
    residuals={k:row[k] for k in ('initial_hamiltonian_difference','initial_raw_vertex_difference',
        'initial_projector_zero_order_difference','maximum_formal_reference_residual','maximum_symmetric_vertex_imaginary_part')}
    residuals['LLL_reference_difference']=float(np.max(abs(row['reference_vertex_difference_orders'][:,0])))
    residuals['intrinsic_jet_change']=max(float(np.max(abs(a.value-b.value))) for a,b in zip(domain.baseline,domain.fields))
    failures={k:v for k,v in residuals.items() if not np.isfinite(v) or v>3e-11}
    per={}
    for i,case in enumerate(cases):
        envelope=np.max(row['unsigned_vertex_envelope'][i],axis=-1)
        decreasing=bool(np.all(envelope[1:]<envelope[:-1])) if i else None
        if i and (not decreasing or envelope[-1]>=envelope[0]/20):failures[case['name']]='sampled source-vertex tail does not decrease sufficiently'
        per[case['name']]={'mass':case['mass'],'angular':case['angular'],
            'unsigned_vertex_envelope':row['unsigned_vertex_envelope'][i].tolist(),
            'signed_pair_reference_difference':row['signed_pair_reference_difference'][i].tolist(),
            'sampled_decay_power':row['sampled_decay_power'][i].tolist(),
            'decay_power_resolved':row['decay_power_resolved'][i].tolist(),
            'maximum_vertex_envelope':envelope.tolist(),'monotone_sample_decay':decreasing,
            'endpoint_reduction_factor':float(envelope[0]/envelope[-1]) if i else None}
    return {'schema':'NSC-INCOMING-CAUCHY-DATA-v1','accountable_author':'Douglas Ek',
        'status':'PASS: conditional same-surface normal-jet domain; physical constraints OPEN' if not failures else 'OPEN',
        'source_hashes':{p:sha(p) for p in SOURCES},'input_hashes':{p:sha(p) for p in INPUTS},
        'locked_inputs':read(INPUTS[0])['locked_inputs'],'allowed_normal_entries':[list(k) for k in NORMAL_ENTRIES],
        'normal_slots':2*len(NORMAL_ENTRIES),'probe_changes':row['changed_normal_entries'],
        'probe_normal_geometry':domain.normal_geometry(),'momentum_magnitudes':row['momentum_magnitudes'].tolist(),
        'per_control_channel':per,'residuals':residuals,'algebra_tolerance':3e-11,'failures':failures,
        'domain':{'surface':'same incoming rho=1 slice; T=0 is its local Taylor origin',
            'fixed':'intrinsic/spatial jets, lapse, shift, normal, spin frame and abstract canonical C0',
            'free':'off-shell normal a/r Taylor slots through total order4; independent physical data still constrained by equations/gauge',
            'background_gauge':'N=1,beta=0; lapse/shift equations are still independent constraints',
            'state_identification':'same Hilbert space, no between-surface isometry and no tensor copied to another radius',
            'probe_values':'numerical derivative controls only, not selected physical Cauchy data'},
        'gate':{**row['scope'],'definition_and_sampled_reference_screen':'PASS' if not failures else 'OPEN',
            'incoming_lapse_shift_constraints':'OPEN','physical_EndpointBranchJets':'OPEN','extended_stationarity':'OPEN',
            'next_owner':'same-slice state/reference/local energy and momentum balance on the free normal data',
            'stress':None,'nulls':None,'metric_timestep':False,'Z3':'OUT OF SCOPE','PDF_bumped':False},
        'comparison':{'fields':'all','float_atol':3e-13,'float_rtol':3e-13,'exact':'hashes, labels, masks, domain and gate'},
        'reproducer':'python3 scripts/derive_nsc_incoming_cauchy_data.py --check'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args();r=make_record()
    if args.check:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(read(OUTPUT),r)
    else:(ROOT/OUTPUT).write_text(json.dumps(r,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':r['status'],'residuals':r['residuals'],
                     'sample_reduction_factors':{k:v['endpoint_reduction_factor'] for k,v in r['per_control_channel'].items()},
                     'failures':r['failures']},indent=2),flush=True)
    if r['failures']:raise SystemExit(1)

if __name__=='__main__':main()
