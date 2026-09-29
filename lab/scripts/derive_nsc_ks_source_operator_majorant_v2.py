#!/usr/bin/env python3
"""Successor source insertion after replaying the complete middle window."""
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
OUTPUT='results/development/nsc-ks-source-operator-majorant-v2.json'
OWNERS=('scripts/derive_nsc_ks_source_operator_majorant_v2.py',
 'scripts/derive_nsc_ks_source_operator_majorant.py',
 'src/recursive_horizons/nsc_middle_source_coverage.py',
 'tests/test_nsc_ks_source_operator_majorant_v2.py','docs/nsc-ks-source-operator-majorant-v2.md')


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


def calculate():
    high,low,old=(read_bound(p) for p in (HIGH,LOW,MIDDLE))
    archive=RetainedUpstreamArchive(ROOT)
    # A real replay, not trusting the previous stdout aggregate or counting files.
    window=cover(ROOT,ROOT/WINDOW,mode='check')
    rows=merge_bounds(high,low,old,window)
    rho=float.fromhex(low['rho_up_hex'])
    if high['rho_up_hex']!=low['rho_up_hex']:raise ValueError('same upstream slice required')
    family,identity=require_original_history(ROOT);entries=archive.family_entries((14,1))
    if any(b.rho_up!=rho for b,_ in entries):raise ValueError('upstream archive slice changed')
    channel=archive.meta['channels'][14]
    if channel['copy_count']*channel['degeneracy']/2!=12:raise ValueError('multiplicity changed')
    with ctx.workprec(192):
        S0,S1,covered,missing=source_error_moments(archive,entries,rows)
        if len(covered)!=866 or missing!=302:raise ValueError('source census changed')
        bound=history_source_insertion(family,rho,S0,S1,
          mass=arb(channel['compact_mass']).union(arb.pi()/2),
          absolute_angular=arb(channel['angular_eigenvalue']).union(arb(5).sqrt()),multiplicity=12)
        inputs={**archive.input_hashes,**{p:digest(p) for p in (HIGH,LOW,MIDDLE,HISTORY_RECORD)}}
        for row in window['rows']:
            for name in ('record.json','witness.npz'):
                path=str(Path(WINDOW)/row['checkpoint']/name);inputs[path]=digest(path)
        return {'schema':'NSC-KS-SOURCE-OPERATOR-MAJORANT-v2',
          'status':'continuous source-error bound for 866 signed rows; full upstream component OPEN',
          'profile_identity':identity,'rho_up_hex':rho.hex(),'family':[14,1],
          'multiplicity_per_signed_family':12,'weighted_error_moment':exact_upper(S0),
          'weighted_energy_error_moment':exact_upper(S1),'continuous_action_error_upper':bound,
          'source_coverage':{'covered_signed_rows':866,'remaining_signed_rows_in_family':302,
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
