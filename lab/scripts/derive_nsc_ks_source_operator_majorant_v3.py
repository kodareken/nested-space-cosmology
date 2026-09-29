#!/usr/bin/env python3
"""Source insertion with the complete middle window and replayed low row15."""
import argparse,json,sys
from pathlib import Path
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from derive_nsc_ks_source_operator_majorant import (
 HIGH,LOW,MIDDLE,read_bound,selected_bounds,source_error_moments,digest)
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper,restored_upper
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.nsc_ks_massive_cubic_uv import require_original_history
from recursive_horizons.nsc_ks_current_uv_transport import HISTORY_RECORD
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_ks_source_operator_majorant import history_source_insertion
from recursive_horizons.nsc_middle_source_coverage import cover,checkpoint_relative,PANEL
WINDOW='results/development/nsc-middle-source-coverage-v1'
OUTPUT='results/development/nsc-ks-source-operator-majorant-v3.json'
LOW15='results/development/nsc-subgap-row15-upstream-v1.json'
from derive_nsc_subgap_row15_upstream import verify_saved as verify_low15
OWNERS=('scripts/derive_nsc_ks_source_operator_majorant_v3.py',
 'scripts/derive_nsc_subgap_row15_upstream.py',
 'scripts/derive_nsc_ks_source_operator_majorant.py',
 'src/recursive_horizons/nsc_middle_source_coverage.py',
 'tests/test_nsc_ks_source_operator_majorant_v3.py','docs/nsc-ks-source-operator-majorant-v3.md')


def window_rows(window):
    """Normalize replayed window rows; no averaging, extrapolation or duplicated old row."""
    if (window.get('coverage_complete') is not True or window.get('completed_positive_rows')!=48
        or window.get('completed_signed_rows')!=96 or window.get('missing')!=[]
        or window.get('covariance_errors_unweighted') is not True):
        raise ValueError('complete replayed middle window required')
    result=[];seen=set()
    for row in window['rows']:
        index=row['row']
        if type(index) is not int or not 0<=index<48 or index in seen or row['panel']!=PANEL:
            raise ValueError('duplicate or relabeled middle row')
        seen.add(index);E=float.fromhex(row['energy_hex'])
        if not 16<=E<32 or row['negative_energy_hex']!=(-E).hex():
            raise ValueError('signed middle energy mismatch')
        if row['checkpoint']!=checkpoint_relative(PANEL,index):
            raise ValueError('middle checkpoint path mismatch')
        comparisons=row['source_comparisons']
        if len(comparisons)!=2 or {c['energy_sign'] for c in comparisons}!={-1,1}:
            raise ValueError('separate signed source bounds required')
        owner=str(Path(WINDOW)/row['checkpoint']/'record.json')
        for c in comparisons:
            result.append({**c,'panel':PANEL,'row':index,'owner':owner,
              'energy_hex':(c['energy_sign']*E).hex(),
              'epsilon':c['physical_source_operator_error_upper']})
    if seen!=set(range(48)):raise ValueError('complete middle row census required')
    return result


def merge_bounds(high,low,old,window):
    return [r for r in selected_bounds(high,low,old) if r["owner"]!=MIDDLE]+window_rows(window)


def add_low_pair(rows,record,rho_hex):
    if (record.get('schema')!='NSC-SUBGAP-ROW15-UPSTREAM-v1'
        or record.get('source_panel')!='group14/low16_1' or record.get('source_row')!=15
        or record.get('rho_up_hex')!=rho_hex or record.get('covariance_error_unweighted') is not True
        or record.get('quadrature_weight_applied') is not False):
        raise ValueError('original unweighted low row15 on the same slice required')
    keys={(r['panel'],r['row'],r['energy_sign']) for r in rows}
    added=[]
    for sign,name in ((1,'positive'),(-1,'negative')):
        key=(record['source_panel'],15,sign)
        if key in keys:raise ValueError('duplicate low source row')
        bound=record.get(name+'_covariance_error_upper')
        if bound is None:raise ValueError('missing low source bound')
        added.append({'panel':key[0],'row':15,'energy_sign':sign,'angular_sign':sign,
          'energy_hex':record[name+'_energy_hex'],'owner':LOW15,'epsilon':bound,
          'source_digest':record[name+'_source_digest'],
          'preparation_digest':record[name+'_preparation_digest']})
    return rows+added


def calculate():
    high,low,old=(read_bound(p) for p in (HIGH,LOW,MIDDLE))
    archive=RetainedUpstreamArchive(ROOT)
    # A real replay, not trusting the previous stdout aggregate or counting files.
    window=cover(ROOT,ROOT/WINDOW,mode='check')
    rows=merge_bounds(high,low,old,window)
    rho=float.fromhex(low['rho_up_hex'])
    low15=verify_low15(ROOT)
    rows=add_low_pair(rows,low15,low['rho_up_hex'])
    if high['rho_up_hex']!=low['rho_up_hex']:raise ValueError('same upstream slice required')
    family,identity=require_original_history(ROOT);entries=archive.family_entries((14,1))
    if any(b.rho_up!=rho for b,_ in entries):raise ValueError('upstream archive slice changed')
    channel=archive.meta['channels'][14]
    if channel['copy_count']*channel['degeneracy']/2!=12:raise ValueError('multiplicity changed')
    with ctx.workprec(192):
        S0,S1,covered,missing=source_error_moments(archive,entries,rows)
        if len(covered)!=868 or missing!=300:raise ValueError('source census changed')
        bound=history_source_insertion(family,rho,S0,S1,
          mass=arb(channel['compact_mass']).union(arb.pi()/2),
          absolute_angular=arb(channel['angular_eigenvalue']).union(arb(5).sqrt()),multiplicity=12)
        inputs={**archive.input_hashes,**{p:digest(p) for p in (HIGH,LOW,MIDDLE,HISTORY_RECORD)}}
        inputs[LOW15]=digest(LOW15)
        inputs[low15['payload']['path']]=digest(low15['payload']['path'])
        for row in window['rows']:
            for name in ('record.json','witness.npz'):
                path=str(Path(WINDOW)/row['checkpoint']/name);inputs[path]=digest(path)
        return {'schema':'NSC-KS-SOURCE-OPERATOR-MAJORANT-v3',
          'status':'continuous source-error bound for 868 signed rows; full upstream component OPEN',
          'profile_identity':identity,'rho_up_hex':rho.hex(),'family':[14,1],
          'multiplicity_per_signed_family':12,'weighted_error_moment':exact_upper(S0),
          'weighted_energy_error_moment':exact_upper(S1),'continuous_action_error_upper':bound,
          'source_coverage':{'covered_signed_rows':868,'remaining_signed_rows_in_family':300,
            'rows':covered,'all_source_families':False,'complete_middle_window_replayed':True},
          'old_middle_pair_replaced_not_added':True,'negative_bounds_supplied_separately':True,
          'source_weights_applied_once':True,'source_quadrature_error_included':False,
          'numerical_field_error_included':False,'value_assembly_arithmetic_included':False,
          'physical_upstream_budget_component':None,'physical_local_gate':'OPEN',
          'source_hashes':implementation_hashes(ROOT,owners=OWNERS),'input_hashes':inputs}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);m=p.add_mutually_exclusive_group(required=True)
    m.add_argument('--record',action='store_true');m.add_argument('--check',action='store_true')
    args=p.parse_args();r=calculate();encoded=(json.dumps(r,sort_keys=True,indent=2,allow_nan=False)+'\n').encode()
    if args.record:publish_exclusive_file(ROOT,OUTPUT,encoded)
    elif (ROOT/OUTPUT).read_bytes()!=encoded:raise ValueError('source insertion or dependencies changed')
    print(json.dumps({'status':r['status'],**{k:float(restored_upper(r['continuous_action_error_upper'][k])) for k in ('N','beta')}},indent=2))
