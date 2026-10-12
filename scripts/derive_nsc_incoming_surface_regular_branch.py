#!/usr/bin/env python3
"""New affine slope/root enclosures; unchanged baseline certificates reused."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_incoming_surface_regular_branch import affine_identities,enclose_branch

OUTPUT=ROOT/'results/development/nsc-incoming-surface-regular-branch.json'
SOURCES=('src/recursive_horizons/nsc_incoming_surface_regular_branch.py',
         'scripts/derive_nsc_incoming_surface_regular_branch.py',
         'tests/test_nsc_incoming_surface_regular_branch.py','docs/nsc-incoming-surface-regular-branch.md')
INPUTS=('results/development/nsc-incoming-lapse-coefficient.json',
        'results/development/nsc-incoming-surface-principal.json',
        'results/development/nsc-mode-resolved-cauchy-state.json',
        'results/development/nsc-incoming-cauchy-data.json',
        'src/recursive_horizons/nsc_incoming_surface_quadratic_lapse.py',
        'src/recursive_horizons/nsc_incoming_source_quadrature_bound.py',
        'src/recursive_horizons/nsc_incoming_vacuum_tail_bound.py',
        'results/development/nsc-incoming-surface-coefficients.json')
CU_EQUATION='-d*a*L^2/(120*pi*r)*(4*m^2*H_r/M^6+(H_r-H_a)/M^4)'


def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()


def authenticate(record):
    for field in ('source_hashes','input_hashes'):
        for path,expected in record.get(field,{}).items():
            if digest(path)!=expected:raise ValueError('inherited branch input changed: '+path)
    def visit(value):
        if isinstance(value,dict):
            if 'path' in value and 'sha256' in value and digest(value['path'])!=value['sha256']:
                raise ValueError('inherited artifact changed')
            for item in value.values():visit(item)
        elif isinstance(value,list):
            for item in value:visit(item)
    visit(record)


def inputs():
    records=[json.loads((ROOT/p).read_text()) for p in INPUTS[:4]]
    for record in records:authenticate(record)
    lapse,principal,inventory,cauchy=records
    if lapse['identity']['reference_group_coefficient']!=CU_EQUATION:raise ValueError('general imported cu equation changed')
    combined=principal['proof']['combined']
    if combined['B_relation']!='B=-(2*A+H_r*d)/a^2' or combined['e_relation']!='e=-d/a^2':
        raise ValueError('general imported principal identities changed')
    if not principal['matrix']['strictly_negative']:raise ValueError('certified baseline principal gate required')
    return lapse,principal,inventory['channels'],cauchy


def make_record():
    lapse,principal,channels,cauchy=inputs()
    proof=affine_identities();branch=enclose_branch(lapse,principal,channels,cauchy)
    independent=json.loads((ROOT/INPUTS[-1]).read_text());authenticate(independent)
    for key in ('reference_C_and_affine_A','local_C_and_affine_A','Delta_w_plus_dC'):
        if independent['exact_polynomials']['identities'][key]!='0':raise ValueError('independent C identity failed')
    observed=branch['intervals']['C_implied_for_independent_check'];saved=independent['coefficients']['C_interval']
    intersection=[max(observed[0],saved[0]),min(observed[1],saved[1])]
    if intersection[0]>intersection[1]:raise ArithmeticError('independent directed C enclosure does not overlap')
    return {'schema':'NSC-INCOMING-SURFACE-REGULAR-BRANCH-v1','accountable_author':'Douglas Ek',
            'status':'PASS: affine normal-data roots and regular component containing zero certified',
            'source_hashes':{p:digest(p) for p in SOURCES},'input_hashes':{p:digest(p) for p in INPUTS},
            'proof':proof,'branch':branch,
            'independent_C_crosscheck':{'record':INPUTS[-1],'derived_from_new_A1':observed,
                'imported_C_interval':saved,'intersection':intersection,'C_used_to_compute_roots':False,
                'C_D_F_producer_called':False},
            'reused_baseline':{'A0_record':INPUTS[0],'principal_record':INPUTS[1],
                               'A0_d_e_producers_rerun':False},
            'residuals':{'exact_affine_derivative_identities':0.,'exact_principal_root_identities':0.,
                         'exact_C_consistency_identities':0.},
            'verification_tolerances':{'exact_identities':0.},
            'scope':{'included_retained_compatible_family_only':True,'normal_data_not_spatial_coordinates':True,
                     'spatial_domain_or_boundary_selected':False,'physical_initial_data_selected':False,
                     'profile_or_global_solution_computed':False,'spacetime_singularity_claimed':False,
                     'C_D_F_producers_rerun':False,'source_or_state_changed':False,
                     'old_probes_or_momentum_quadratures_run':False,'rho0_seam_changed':False,
                     'metric_evolution':False,'Gamma_rest_assigned':False,'push_or_PDF':False},
            'reproducer':'python3 scripts/derive_nsc_incoming_surface_regular_branch.py --check'}


def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true');args=parser.parse_args()
    if args.write and OUTPUT.exists():raise FileExistsError('existing branch certificate is not overwritten')
    data=make_record()
    if args.write:OUTPUT.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    elif json.loads(OUTPUT.read_text())!=data:raise ValueError('branch certificate replay differs')
    print(json.dumps({'status':data['status'],'intervals':data['branch']['intervals'],
                      'regular_component':data['branch']['regular_component_containing_zero'],
                      'residuals':data['residuals']},indent=2))


if __name__=='__main__':main()
