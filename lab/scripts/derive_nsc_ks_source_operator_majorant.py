#!/usr/bin/env python3
"""Authenticate covered source errors and bound their continuous N/beta contribution."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper,restored_upper
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.nsc_ks_massive_cubic_uv import require_original_history
from recursive_horizons.nsc_ks_current_uv_transport import HISTORY_RECORD
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_ks_source_operator_majorant import history_source_insertion

OUTPUT='results/development/nsc-ks-source-operator-majorant-v1.json'
HIGH='results/development/nsc-vacuum-source-remainder-v1.json'
LOW='results/development/nsc-subgap-upstream-covariance-v1.json'
MIDDLE='results/development/nsc-vacuum-source-correction-v1.json'
OWNERS=('scripts/derive_nsc_ks_source_operator_majorant.py',
    'src/recursive_horizons/nsc_ks_source_operator_majorant.py',
    'tests/test_nsc_ks_source_operator_majorant.py','docs/nsc-ks-source-operator-majorant.md')


def digest(path):
    target=(ROOT/path).resolve()
    if not target.is_relative_to(ROOT):raise ValueError('source evidence must stay inside the lab')
    return sha256(target.read_bytes()).hexdigest()


def read_bound(path):
    record=json.loads((ROOT/path).read_text())
    for group in ('source_hashes','input_hashes'):
        for name,expected in record[group].items():
            if digest(name)!=expected:raise ValueError('source-bound dependency changed: '+name)
    payload=record.get('payload')
    if payload and digest(payload['path'])!=payload['sha256']:
        raise ValueError('source-bound payload changed')
    return record


def selected_bounds(high,low,middle):
    if (high.get('schema')!='NSC-VACUUM-SOURCE-REMAINDER-v1'
            or low.get('schema')!='NSC-SUBGAP-UPSTREAM-COVARIANCE-v1'
            or middle.get('schema')!='NSC-VACUUM-SOURCE-CORRECTION-v1'
            or high['archived_comparison'].get('unweighted_covariance_errors') is not True
            or low.get('covariance_error_unweighted') is not True
            or middle.get('covariance_errors_unweighted') is not True):
        raise ValueError('owned unweighted preparation covariance bounds required')
    rows=[]
    for row in high['archived_comparison']['rows']:
        rows.append({**row,'epsilon':row['original_source_operator_error_upper'],'owner':HIGH})
    for sign,name in ((1,'positive'),(-1,'negative')):
        rows.append({'panel':low['source_panel'],'row':low['source_row'],
            'energy_sign':sign,'angular_sign':sign,
            'energy_hex':float(sign*float.fromhex(low['positive_energy_hex'])).hex(),
            'preparation_digest':low[name+'_preparation_digest'],'source_digest':low[name+'_source_digest'],
            'epsilon':low[name+'_covariance_operator_error_upper'],'owner':LOW})
    for row in middle['source_comparisons']:
        sign=row['energy_sign']
        rows.append({**row,'panel':middle['panel'],'row':middle['row'],
            'energy_hex':float(sign*float.fromhex(middle['energy_hex'])).hex(),
            'epsilon':row['physical_source_operator_error_upper'],'owner':MIDDLE})
    return rows


def source_error_moments(archive,entries,rows):
    """Reject relabeling, duplicate coverage, absent bounds and source drift."""
    indexed={}
    for batch,_ in entries:
        for index,E in enumerate(batch.source.energies[::3]):
            key=(batch.original_panel,batch.rows[0]+index,batch.energy_sign)
            if key in indexed:raise ValueError('duplicate original source row')
            indexed[key]=(batch,index,float(E))
    S0=arb(0);S1=arb(0);seen=set();covered=[]
    for row in rows:
        key=(row['panel'],row['row'],row['energy_sign'])
        if key in seen:raise ValueError('duplicate bounded source row')
        if key not in indexed:raise ValueError('bounded row is not in the original source')
        batch,index,E=indexed[key]
        if (row['energy_hex']!=E.hex() or row['angular_sign']!=batch.angular_sign
                or row['preparation_digest']!=batch.preparation_digest
                or row['source_digest']!=batch.source.digest):
            raise ValueError('source bound does not match its original preparation')
        if row.get('epsilon') is None:raise ValueError('missing covariance bound')
        epsilon=restored_upper(row['epsilon'])
        if not epsilon.is_finite() or not epsilon>=0:raise ValueError('nonnegative covariance bound required')
        raw_weight=arb(float(archive._arrays[row['panel']+'/weights'][row['row']]))
        exact_measure=raw_weight/(2*arb.pi())
        stored_measure=arb(float(batch.source.column_weights[3*index]))**2
        # Cover both the exact quadrature measure and its stored square-root
        # representation for this SOURCE-error term. Other weight arithmetic
        # in the value evaluator remains a separate component.
        measure=exact_measure.max(stored_measure).upper()
        S0+=measure*epsilon;S1+=measure*abs(arb(E))*epsilon
        seen.add(key)
        covered.append({'panel':key[0],'row':key[1],'energy_sign':key[2],
            'energy_hex':E.hex(),'bound_owner':row['owner'],
            'preparation_digest':batch.preparation_digest,'source_digest':batch.source.digest,
            'covariance_error_upper':row['epsilon'],'measure_upper':exact_upper(measure)})
    return S0.upper(),S1.upper(),covered,len(indexed)-len(seen)


def calculate():
    high,low,middle=(read_bound(path) for path in (HIGH,LOW,MIDDLE))
    rho=float.fromhex(low['rho_up_hex'])
    if high['rho_up_hex']!=low['rho_up_hex'] or middle['rho_up_hex']!=low['rho_up_hex']:
        raise ValueError('same upstream slice required')
    family,identity=require_original_history(ROOT)
    archive=RetainedUpstreamArchive(ROOT);entries=archive.family_entries((14,1))
    if any(batch.rho_up!=rho for batch,_ in entries):raise ValueError('upstream archive slice changed')
    channel=archive.meta['channels'][14]
    multiplicity=channel['copy_count']*channel['degeneracy']/2
    if multiplicity!=12:raise ValueError('original signed-family multiplicity required')
    with ctx.workprec(192):
        S0,S1,rows,missing=source_error_moments(archive,entries,selected_bounds(high,low,middle))
        if len(rows)!=772 or missing!=396:raise ValueError('selected source coverage changed')
        m=arb(channel['compact_mass']).union(arb.pi()/2)
        ell=arb(channel['angular_eigenvalue']).union(arb(5).sqrt())
        bound=history_source_insertion(family,rho,S0,S1,mass=m,absolute_angular=ell,multiplicity=12)
        result={'schema':'NSC-KS-SOURCE-OPERATOR-MAJORANT-v1',
            'status':'ENCLOSED: continuous source-error contribution from 772 covered signed rows; full upstream budget OPEN',
            'profile_identity':identity,'rho_up_hex':rho.hex(),'family':[14,1],
            'multiplicity_per_signed_family':12,'weighted_error_moment':exact_upper(S0),
            'weighted_energy_error_moment':exact_upper(S1),'continuous_action_error_upper':bound,
            'source_coverage':{'covered_signed_rows':772,'remaining_signed_rows_in_family':396,
                'other_positive_angular_family_included':False,'all_source_families':False,
                'rows':rows},
            'intrinsic_vertices_unchanged':True,'negative_bounds_supplied_separately':True,
            'source_weights_applied_once':True,'source_quadrature_error_included':False,
            'value_assembly_arithmetic_included':False,'numerical_field_error_included':False,
            'physical_upstream_budget_component':None,'physical_local_gate':'OPEN',
            'source_hashes':implementation_hashes(ROOT,owners=OWNERS),
            'input_hashes':{**archive.input_hashes,
                **{p:digest(p) for p in (HIGH,LOW,MIDDLE,HISTORY_RECORD)}}}
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record',action='store_true');mode.add_argument('--check',action='store_true')
    args=parser.parse_args();result=calculate()
    encoded=(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n').encode()
    if args.record:publish_exclusive_file(ROOT,OUTPUT,encoded)
    elif (ROOT/OUTPUT).read_bytes()!=encoded:raise ValueError('source action bound or dependencies changed')
    print(json.dumps({'status':result['status'],
        'N':float(restored_upper(result['continuous_action_error_upper']['N'])),
        'beta':float(restored_upper(result['continuous_action_error_upper']['beta'])),
        'remaining_signed_rows_in_family':result['source_coverage']['remaining_signed_rows_in_family']},indent=2))
