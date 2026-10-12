#!/usr/bin/env python3
"""Apply the retained-source flux bound to the compact history's incoming data."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'tests'))
from recursive_horizons.nsc_incoming_constraint_gate import retained_incoming_flux_bound
from test_nsc_incoming_constraint_gate import input_data,sewing_current_identity

OUTPUT='results/development/nsc-incoming-constraint-gate.json'
SOURCES=('src/recursive_horizons/nsc_incoming_constraint_gate.py','tests/test_nsc_incoming_constraint_gate.py',
         'scripts/derive_nsc_incoming_constraint_gate.py','docs/nsc-incoming-constraint-gate.md')
INPUTS=('results/development/nsc-pg-retained-covariance.json','results/development/nsc-mode-resolved-cauchy-state.json',
        'results/development/nsc-compact-matched-restart.json','src/recursive_horizons/nsc_compact_ctp_neck.py',
        'src/recursive_horizons/nsc_paired_horizon_preparation.py','src/recursive_horizons/nsc_pg_massive_modes.py',
        'src/recursive_horizons/nsc_transmitting_history_modes.py','src/recursive_horizons/nsc_transmitting_resolvent.py',
        'src/recursive_horizons/nsc_common_ks_trace.py','src/recursive_horizons/nsc_general_ks_reference.py',
        'src/recursive_horizons/nsc_lll_geometric_history.py','docs/nsc-declared-action-scope.md',
        'docs/nsc-common-ks-trace.md','docs/nsc-smooth-seam-variation.md')

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def make_record():
    ledger,channels,prep=input_data();bound=retained_incoming_flux_bound(ledger,channels,prep)
    identities=sewing_current_identity();tol=3e-11
    miss=bound['rational_witness']['normalized_shift_miss_lower_decimal']
    passed=all(v=='0' for v in identities.values()) and miss>tol
    return {'schema':'NSC-INCOMING-CONSTRAINT-GATE-v1','accountable_author':'Douglas Ek',
        'status':'NON-EXISTENCE PASS: compact fixed-incoming response family, retained source only' if passed else 'OPEN',
        'source_hashes':{p:sha(p) for p in SOURCES},'input_hashes':{p:sha(p) for p in INPUTS},
        'locked_inputs':ledger,'preparation_inputs':{'surface_gravity':prep['surface_gravity'],'Omega':prep['omega']},
        'source_bound':bound,'exact_identities':identities,
        'constraint':{'necessary_equation':'zero geometric/local normal-axial momentum = T01/(2*A) on fixed incoming jets',
            'incoming_T01_strict_lower':bound['rational_witness']['T01_lower_decimal'],
            'normalized_shift_miss_strict_lower':miss,'dimensionless_tolerance':tol,
            'miss_over_tolerance_lower':miss/tol,'geometric_local_reference_momentum':0.,
            'nature':'analytic inequality, not an evaluated full tensor or copied neck residual'},
        'scope':{'retained_groups':33,'metric_family':'SuppliedKSHarmonicMetric; all incoming jets fixed at rho=1',
            'state':'same reference incoming preparation, preserved by causal conditional interior evolution',
            'regularity':'smooth Lorentzian histories with well-defined continuous selected source',
            'argument':'nonzero boundary limit persists in an adjacent collar; full local equation fails',
            'conditional_response_valid':True,'four_amplitude_extremum_is_full_stationarity':False,
            'infinite_unretained_tower_included':False,'old_homogeneous_certificate':'preserved, not rerun'},
        'gate':{'compact_response_family_as_physical_solution':'NON-EXISTENCE PASS' if passed else 'OPEN',
            'extended_transmitting_theory':'OPEN','next_owner':'constraint-compatible incoming spherical geometry and consistent state-preparation variation',
            'new_stress_tensor':None,'new_metric':None,'V_c':None,'metric_timestep':False,'Z3':'OUT OF SCOPE',
            'extra_Gamma_rest':'none','PDF_bumped':False,'Weyl_time_node_diagnostic':93.54264532195464},
        'computation':{'mode_solve':False,'stress_quadrature':False,'Einstein_solve':False,
            'operations':'retained multiplicity sum, scalar bounds, exact rational sign witness and sewing algebra'},
        'comparison':{'fields':'all','float_atol':3e-13,'float_rtol':3e-13,'exact':'hashes, fractions, identities, scope'},
        'reproducer':'python3 scripts/derive_nsc_incoming_constraint_gate.py --check'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args();record=make_record()
    if args.check:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(json.loads((ROOT/OUTPUT).read_text()),record)
    else:(ROOT/OUTPUT).write_text(json.dumps(record,sort_keys=True,indent=2)+'\n')
    print(json.dumps({k:record[k] for k in ('status','constraint','exact_identities')},indent=2),flush=True)
    if not record['status'].startswith('NON-EXISTENCE PASS'):raise SystemExit(1)

if __name__=='__main__':main()
