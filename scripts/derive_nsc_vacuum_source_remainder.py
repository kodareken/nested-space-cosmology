#!/usr/bin/env python3
"""Record uniform homogeneous source-covariance remainder bounds for group14."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

import numpy as np
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper,restored_upper
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_vacuum_source_remainder import VacuumSourceExpansion

OUTPUT='results/development/nsc-vacuum-source-remainder-v1.json'
UPSTREAM='results/development/nsc-subgap-upstream-covariance-v1.json'
HORIZON='results/development/nsc-metric-horizon-frame-v1.json'
INVENTORY='results/development/nsc-ks-source-inventory.json'
OWNERS=('scripts/derive_nsc_vacuum_source_remainder.py',
    'src/recursive_horizons/nsc_vacuum_source_remainder.py',
    'tests/test_nsc_vacuum_source_remainder.py','docs/nsc-vacuum-source-remainder.md')


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()


def ball(value):
    return restored_upper(value['lower']).union(restored_upper(value['upper']))


def calculate():
    upstream=json.loads((ROOT/UPSTREAM).read_text())
    horizon=json.loads((ROOT/HORIZON).read_text())
    inventory=json.loads((ROOT/INVENTORY).read_text())
    payload=inventory['payload']
    if digest(payload['path'])!=payload['sha256']:
        raise ValueError('original source inventory payload changed')
    with np.load(ROOT/payload['path'],allow_pickle=False) as data:
        metadata=json.loads(data['metadata_json'].tobytes())
    config=metadata['config'];label=metadata['panels']['group14/low16_1']
    rho=float.fromhex(upstream['rho_up_hex'])
    archive=RetainedUpstreamArchive(ROOT)
    entries=archive.family_entries((14,1))
    rows=[];stored=[]
    with ctx.workprec(192):
        mass,angular=ball(horizon['mass']),ball(horizon['angular'])
        if (label['group']!=14 or label['angular_sign']!=1
                or not mass.contains(arb(label['mass']))
                or not angular.contains(arb(label['angular_magnitude']))
                or upstream['source_panel']!='group14/low16_1'):
            raise ValueError('original group14 positive angular labels required')
        for order,cells in ((4,32),(6,32),(8,32),(8,64)):
            model=VacuumSourceExpansion(config['horizon_rho'],mass,angular,order=order,bits=192)
            constant=model.remainder_constant(rho,cells=cells)
            bound=model.covariance_error(constant,160,kappa=config['surface_gravity'],omega=config['omega'])
            rows.append({'order':order,'cells':cells,
                'vacuum_bloch_constant_upper':exact_upper(constant),
                'at_cutoff_160':{key:exact_upper(value) for key,value in bound.items()}})
        comparison=VacuumSourceExpansion(config['horizon_rho'],mass,angular,order=16,bits=192)
        comparison_constant=comparison.remainder_constant(rho,cells=128)
        worst=arb(0);seen=set()
        for batch,_ in entries:
            if (batch.rho_up!=rho or not mass.contains(arb(batch.mass))
                    or not angular.contains(arb(abs(batch.angular)))):
                raise ValueError('original high-energy source channel or target changed')
            for index,E in enumerate(batch.source.energies[::3]):
                if abs(E)<32:continue
                key=(batch.original_panel,batch.rows[0]+index,batch.energy_sign)
                if key in seen:raise ValueError('duplicate source row')
                seen.add(key);sl=slice(3*index,3*index+3)
                bound=comparison.compare_original(rho,abs(float(E)),
                    batch.initial_columns[:,sl],batch.source.covariance[sl,sl],comparison_constant,
                    kappa=config['surface_gravity'],omega=config['omega'],
                    negative_partner=batch.energy_sign<0)
                worst=worst.max(bound['original_source_operator_error'])
                stored.append({'panel':key[0],'row':key[1],'energy_sign':batch.energy_sign,
                    'angular_sign':batch.angular_sign,'energy_hex':float(E).hex(),
                    'preparation_digest':batch.preparation_digest,'source_digest':batch.source.digest,
                    'original_source_operator_error_upper':exact_upper(bound['original_source_operator_error']),
                    'polynomial_to_stored_operator_distance_upper':exact_upper(bound['polynomial_to_stored_operator_distance'])})
        if len(stored)!=768 or sum(r['energy_sign']==1 for r in stored)!=384:
            raise ValueError('original group14 high-energy row coverage changed')
        result={'schema':'NSC-VACUUM-SOURCE-REMAINDER-v1',
            'status':'ENCLOSED: uniform homogeneous source-covariance remainder; changed-history UV gate OPEN',
            'group':14,'expansion_energy_sign':1,'expansion_angular_sign':1,'rho_up_hex':rho.hex(),
            'mass_hex':float(label['mass']).hex(),'angular_hex':float(label['angular_magnitude']).hex(),
            'source_kappa_hex':float(config['surface_gravity']).hex(),
            'source_omega_hex':float(config['omega']).hex(),
            'expansion_bound_domain':{'lower':0,'lower_included':False,'upper':'unbounded'},
            'reported_high_energy_cutoff_hex':float(160).hex(),
            'bits':192,'metric_terms':48,'cases':rows,
            'vacuum_bloch_error_formula':'C_prep/E^M',
            'vacuum_covariance_error_formula':'C_prep/(2 E^M)',
            'source_nonvacuum_error_formula':'exp(-2*pi*E/kappa)+exp(-2*pi*E/(omega*kappa))+exp(-pi*E/kappa)',
            'archived_comparison':{'order':16,'cells':128,'absolute_energy_lower_hex':float(32).hex(),
                'vacuum_bloch_constant_upper':exact_upper(comparison_constant),
                'signed_rows':768,'positive_rows':384,'unweighted_covariance_errors':True,
                'remaining_signed_rows_in_family':sum(len(b.source.energies)//3 for b,_ in entries)-len(stored),
                'maximum_original_source_operator_error_upper':exact_upper(worst),
                'rows':stored,'source_quadrature_error_included':False,
                'other_positive_angular_family_included':False},
            'whole_cell_ranges_used':True,'order_difference_used_as_error':False,
            'reflection_assumed_zero':False,'source_occupations_replaced':False,
            'upstream_columns_replaced':False,'exact_normalized_horizon_vacuum_branch':True,
            'same_column_H2_remainder':None,'changed_history_C_M':None,
            'integrated_N_beta_UV_tail':None,'physical_upstream_budget_component':None,
            'physical_UV_budget_component':None,'all_source_families':False,
            'physical_local_gate':'OPEN',
            'source_hashes':implementation_hashes(ROOT,owners=OWNERS),
            'input_hashes':{**archive.input_hashes,
                           **{p:digest(p) for p in (UPSTREAM,HORIZON,INVENTORY,payload['path'])}}}
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record',action='store_true');mode.add_argument('--check',action='store_true')
    args=parser.parse_args();start=time.process_time();result=calculate()
    encoded=(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n').encode()
    if args.record:publish_exclusive_file(ROOT,OUTPUT,encoded)
    elif (ROOT/OUTPUT).read_bytes()!=encoded:
        raise ValueError('source remainder, domain or dependencies changed')
    print(json.dumps({'status':result['status'],'cpu_seconds':time.process_time()-start,
        'cases':[{'order':c['order'],'cells':c['cells'],
                  'C_prep_upper':float(restored_upper(c['vacuum_bloch_constant_upper'])),
                  'source_operator_error_at_160':float(restored_upper(c['at_cutoff_160']['source_operator_error']))}
                 for c in result['cases']],
        'archived_signed_rows':result['archived_comparison']['signed_rows'],
        'maximum_original_source_operator_error_upper':float(restored_upper(
            result['archived_comparison']['maximum_original_source_operator_error_upper']))},indent=2))
