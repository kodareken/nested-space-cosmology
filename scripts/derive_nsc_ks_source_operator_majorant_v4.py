#!/usr/bin/env python3
"""Successor source insertion for the replayed direct-vacuum window.

Replays the v3 aggregate and the saved direct-window traces. It does not solve,
apply a quadrature weight twice, or close the local gate.
"""
import argparse,json,sys
from hashlib import sha256
from pathlib import Path
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from derive_nsc_ks_source_operator_majorant import digest,read_bound,source_error_moments
from derive_nsc_ks_source_operator_majorant_v3 import calculate as replay_v3
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper,restored_upper
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.nsc_ks_massive_cubic_uv import require_original_history
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_ks_source_operator_majorant import history_source_insertion
from recursive_horizons.nsc_direct_source_window import (
 ENERGY_LOWER,ENERGY_UPPER,EXPECTED_POSITIVE_ROWS,LOW16_PANEL,LOW16_ROWS,LOW32_PANEL,
 LOW32_ROWS,MAX_ENERGY,MIN_ENERGY,SCHEMA,checkpoint_relative,cover)
WINDOW='results/development/nsc-direct-source-window-v1'
PRIOR='results/development/nsc-ks-source-operator-majorant-v3.json'
OUTPUT='results/development/nsc-ks-source-operator-majorant-v4.json'
PRIOR_COVERED,PRIOR_REMAINING,DIRECT_POSITIVE,COVERED,REMAINING=868,300,87,1042,126
DIRECT_SIGNED=2*DIRECT_POSITIVE
OWNERS=('scripts/derive_nsc_ks_source_operator_majorant_v4.py',
 'scripts/derive_nsc_ks_source_operator_majorant_v3.py',
 'scripts/derive_nsc_ks_source_operator_majorant.py',
 'scripts/derive_nsc_direct_source_window.py',
 'scripts/derive_nsc_subgap_row15_upstream.py',
 'src/recursive_horizons/nsc_direct_source_window.py',
 'src/recursive_horizons/nsc_direct_vacuum_source.py',
 'src/recursive_horizons/nsc_middle_source_coverage.py',
 'src/recursive_horizons/nsc_ks_source_operator_majorant.py',
 'tests/test_nsc_ks_source_operator_majorant_v4.py')


def prior_rows(record):
    """Rows already authenticated by the v3 replay. Their stored errors stay unweighted."""
    coverage=record.get('source_coverage') if isinstance(record,dict) else None
    if (not isinstance(coverage,dict) or record.get('schema')!='NSC-KS-SOURCE-OPERATOR-MAJORANT-v3'
        or record.get('family')!=[14,1] or record.get('physical_local_gate')!='OPEN'
        or record.get('physical_upstream_budget_component') is not None
        or record.get('source_weights_applied_once') is not True
        or record.get('source_quadrature_error_included') is not False
        or record.get('numerical_field_error_included') is not False
        or record.get('negative_bounds_supplied_separately') is not True
        or not isinstance(record.get('rho_up_hex'),str)
        or coverage.get('covered_signed_rows')!=PRIOR_COVERED
        or coverage.get('remaining_signed_rows_in_family')!=PRIOR_REMAINING
        or coverage.get('all_source_families') is not False
        or coverage.get('complete_middle_window_replayed') is not True):
        raise ValueError('authenticated v3 family14_1 coverage required')
    rows=[]
    for row in coverage.get('rows') or []:
        if not isinstance(row,dict) or row.get('covariance_error_upper') is None:
            raise ValueError('missing v3 source bound')
        if row.get('energy_sign') not in (-1,1) or row.get('bound_owner') is None:
            raise ValueError('missing v3 source bound')
        rows.append({**row,'angular_sign':row['energy_sign'],'owner':row['bound_owner'],
          'epsilon':row['covariance_error_upper']})
    if len(rows)!=PRIOR_COVERED:raise ValueError('authenticated v3 family14_1 coverage required')
    return rows


def direct_rows(window):
    """Normalize replayed direct pairs. Quadrature weights remain unapplied here."""
    if (EXPECTED_POSITIVE_ROWS!=DIRECT_POSITIVE or SCHEMA!='NSC-DIRECT-SOURCE-WINDOW-v1'
        or LOW16_ROWS!=tuple(range(16,48)) or LOW32_ROWS!=tuple(range(41,96))
        or (LOW16_PANEL,LOW32_PANEL)!=('group14/low16_1','group14/low32_1')):
        raise ValueError('direct window census changed')
    if not isinstance(window,dict):raise ValueError('complete replayed direct window required')
    if (window.get('mode')!='check' or window.get('replay_uses_saved_trace') is not True
        or window.get('new_rows_captured')!=0):
        raise ValueError('replayed direct window required')
    if window.get('weight_applied_to_error') is not False:
        raise ValueError('direct source weight already applied')
    status=window.get('status')
    if (window.get('schema')!=SCHEMA or window.get('coverage_complete') is not True
        or window.get('completed_positive_rows')!=DIRECT_POSITIVE
        or window.get('completed_negative_partners')!=DIRECT_POSITIVE
        or window.get('completed_signed_rows')!=DIRECT_SIGNED
        or window.get('missing')!=[] or window.get('missing_positive_rows')!=0
        or window.get('covariance_errors_unweighted') is not True
        or window.get('source_quadrature_error_included') is not False
        or window.get('correction_forcing_used') is not False
        or window.get('all_source_families') is not False
        or window.get('other_positive_angular_family_included') is not False
        or window.get('archived_source_replaced') is not False
        or window.get('finite_window_closes_gate') is not False
        or window.get('physical_upstream_budget_component') is not None
        or window.get('changed_history_C_M') is not None
        or window.get('physical_local_gate')!='OPEN' or window.get('group')!=14
        or window.get('angular_sign')!=1 or window.get('panels')!=[LOW16_PANEL,LOW32_PANEL]
        or window.get('energy_window')!={'lower':ENERGY_LOWER,'lower_included':True,
          'upper':ENERGY_UPPER,'upper_included':False}
        or not isinstance(status,str) or not status.startswith('OPEN') or 'ENCLOSED' in status):
        raise ValueError('complete replayed direct window required')
    rows=window.get('rows')
    if not isinstance(rows,list) or len(rows)!=DIRECT_POSITIVE:
        raise ValueError('complete direct row census required')
    if ((rows[0].get('panel'),rows[0].get('row'))!=(LOW32_PANEL,LOW32_ROWS[0])
        or (rows[-1].get('panel'),rows[-1].get('row'))!=(LOW16_PANEL,LOW16_ROWS[-1])):
        raise ValueError('direct-window energy endpoints changed')
    result=[];seen=set();low16=set();low32=set()
    for row in rows:
        if not isinstance(row,dict):raise ValueError('complete direct row census required')
        panel,index=row.get('panel'),row.get('row')
        if type(index) is not int or panel not in (LOW16_PANEL,LOW32_PANEL):
            raise ValueError('direct window panel changed')
        if (panel,index) in seen:raise ValueError('duplicate direct source row')
        seen.add((panel,index));(low16 if panel==LOW16_PANEL else low32).add(index)
        if (row.get('weight_applied_to_error') is not False
            or row.get('correction_forcing_used') is not False):
            raise ValueError('direct source weight already applied')
        if row.get('physical_local_gate')!='OPEN' or row.get('physical_upstream_budget_component') is not None:
            raise ValueError('direct window must stay open')
        if row.get('checkpoint')!=checkpoint_relative(panel,index):
            raise ValueError('direct checkpoint path mismatch')
        fiber,negative=row.get('energy_fiber_hex'),row.get('negative_energy_fiber_hex')
        if (not isinstance(fiber,list) or len(fiber)!=3 or len(set(fiber))!=1
            or not all(isinstance(item,str) for item in fiber)):
            raise ValueError('signed direct energy mismatch')
        try:E=float.fromhex(fiber[0])
        except ValueError as error:raise ValueError('signed direct energy mismatch') from error
        if not ENERGY_LOWER<=E<ENERGY_UPPER or negative!=[(-E).hex()]*3:
            raise ValueError('signed direct energy mismatch')
        comparisons=row.get('signed_comparisons')
        if (not isinstance(comparisons,list) or len(comparisons)!=2
            or not all(isinstance(item,dict) for item in comparisons)
            or [item.get('energy_sign') for item in comparisons]!=[1,-1]
            or [item.get('angular_sign') for item in comparisons]!=[1,-1]):
            raise ValueError('separate signed source bounds required')
        if (comparisons[0].get('source_digest')==comparisons[1].get('source_digest')
            or comparisons[0].get('preparation_digest')==comparisons[1].get('preparation_digest')):
            raise ValueError('duplicate source identity')
        owner=str(Path(WINDOW)/row['checkpoint']/'record.json')
        for sign,item,expected in ((1,comparisons[0],fiber),(-1,comparisons[1],negative)):
            if (item.get('signed_total_upper') is None or item.get('source_digest') is None
                or item.get('preparation_digest') is None):
                raise ValueError('missing direct source bound')
            if item.get('energy_fiber_hex')!=expected:
                raise ValueError('signed direct energy mismatch')
            result.append({'panel':panel,'row':index,'energy_sign':sign,'angular_sign':sign,
              'energy_hex':expected[0],'owner':owner,'epsilon':item['signed_total_upper'],
              'source_digest':item['source_digest'],
              'preparation_digest':item['preparation_digest']})
    if (low16!=set(LOW16_ROWS) or low32!=set(LOW32_ROWS) or len(result)!=DIRECT_SIGNED
        or float.fromhex(rows[0]['energy_fiber_hex'][0])!=MIN_ENERGY
        or float.fromhex(rows[-1]['energy_fiber_hex'][0])!=MAX_ENERGY):
        raise ValueError('complete direct row census required')
    return result


def add_direct_pairs(rows,direct):
    """Append original direct pairs once. An identity already in the v3 coverage is refused."""
    keys={(row['panel'],row['row'],row['energy_sign']) for row in rows}
    if len(keys)!=len(rows):raise ValueError('duplicate bounded source row')
    for item in direct:
        key=(item['panel'],item['row'],item['energy_sign'])
        if key in keys:raise ValueError('duplicate direct source row')
        if item.get('epsilon') is None:raise ValueError('missing direct source bound')
        keys.add(key)
    if len(direct)!=DIRECT_SIGNED or len(keys)!=COVERED:
        raise ValueError('source census changed')
    return list(rows)+list(direct)


def _bound_path(row,name):
    panel,index,checkpoint=row.get('panel'),row.get('row'),row.get('checkpoint')
    if checkpoint!=checkpoint_relative(panel,index):
        raise ValueError('direct checkpoint path mismatch')
    relative=Path(WINDOW)/checkpoint/name
    if relative.is_absolute() or '..' in relative.parts:
        raise ValueError('direct checkpoint path mismatch')
    return relative.as_posix()


def bind_direct_inputs(inputs,window):
    """Bind every direct-window JSON and NPZ, including the digest declared by each record."""
    if not isinstance(inputs,dict) or not isinstance(window,dict) or not isinstance(window.get('rows'),list):
        raise ValueError('dependency binding required')
    bound={}
    for row in window['rows']:
        record_path,witness_path=_bound_path(row,'record.json'),_bound_path(row,'witness.npz')
        try:
            record_bytes=(ROOT/record_path).read_bytes();witness_bytes=(ROOT/witness_path).read_bytes()
        except OSError as error:
            raise ValueError('missing direct window file') from error
        try:record=json.loads(record_bytes)
        except json.JSONDecodeError as error:
            raise ValueError('missing direct window file') from error
        payload=record.get('payload') if isinstance(record,dict) else None
        witness_hash,record_hash=digest(witness_path),digest(record_path)
        if (not isinstance(payload,dict) or payload.get('path')!=f"{row['checkpoint']}/witness.npz"
            or payload.get('sha256')!=witness_hash or payload.get('bytes')!=len(witness_bytes)
            or witness_hash!=sha256(witness_bytes).hexdigest()
            or record_hash!=sha256(record_bytes).hexdigest()):
            raise ValueError('payload hash changed')
        for path,value in ((record_path,record_hash),(witness_path,witness_hash)):
            if path in inputs or path in bound:
                raise ValueError('direct window dependency overlaps an existing binding')
            bound[path]=value
    actual={path.relative_to(ROOT).as_posix()
      for path in (ROOT/WINDOW/'rows').rglob('*') if path.is_file()}
    if actual!=set(bound):raise ValueError('unbound direct window file')
    if len(bound)!=2*DIRECT_POSITIVE:raise ValueError('direct window file census changed')
    inputs.update(bound)
    return inputs


def require_same_slice(record,archive):
    rho_hex=record.get('rho_up_hex') if isinstance(record,dict) else None
    if not isinstance(rho_hex,str):raise ValueError('same upstream slice required')
    rho=float.fromhex(rho_hex)
    if rho.hex()!=rho_hex:raise ValueError('same upstream slice required')
    entries=archive.family_entries((14,1))
    if any(batch.rho_up!=rho for batch,_ in entries):
        raise ValueError('upstream archive slice changed')
    return rho,entries


def calculate():
    prior=replay_v3()
    rows=prior_rows(prior)
    archive=RetainedUpstreamArchive(ROOT)
    rho,entries=require_same_slice(prior,archive)
    # Replay saved traces. Counting directories is not authentication.
    window=cover(ROOT,ROOT/WINDOW,mode='check')
    rows=add_direct_pairs(rows,direct_rows(window))
    family,identity=require_original_history(ROOT)
    if identity!=prior.get('profile_identity'):raise ValueError('profile identity changed')
    channel=archive.meta['channels'][14]
    if channel['copy_count']*channel['degeneracy']/2!=12:raise ValueError('multiplicity changed')
    with ctx.workprec(192):
        S0,S1,covered,missing=source_error_moments(archive,entries,rows)
        if len(covered)!=COVERED or missing!=REMAINING:raise ValueError('source census changed')
        bound=history_source_insertion(family,rho,S0,S1,
          mass=arb(channel['compact_mass']).union(arb.pi()/2),
          absolute_angular=arb(channel['angular_eigenvalue']).union(arb(5).sqrt()),multiplicity=12)
        return {'schema':'NSC-KS-SOURCE-OPERATOR-MAJORANT-v4',
          'status':'continuous source-error bound for 1042 signed rows; full upstream component OPEN',
          'profile_identity':identity,'rho_up_hex':rho.hex(),'family':[14,1],
          'multiplicity_per_signed_family':12,'weighted_error_moment':exact_upper(S0),
          'weighted_energy_error_moment':exact_upper(S1),'continuous_action_error_upper':bound,
          'source_coverage':{'covered_signed_rows':COVERED,'remaining_signed_rows_in_family':REMAINING,
            'rows':covered,'all_source_families':False,'complete_middle_window_replayed':True,
            'complete_direct_source_window_replayed':True,'v3_signed_rows_reused':PRIOR_COVERED,
            'direct_signed_rows_added':DIRECT_SIGNED},
          'old_middle_pair_replaced_not_added':True,'negative_bounds_supplied_separately':True,
          'source_weights_applied_once':True,'source_quadrature_error_included':False,
          'numerical_field_error_included':False,'value_assembly_arithmetic_included':False,
          'physical_upstream_budget_component':None,'physical_local_gate':'OPEN',
          'source_hashes':implementation_hashes(ROOT,owners=OWNERS),
          'input_hashes':bind_direct_inputs(dict(prior['input_hashes']),window)}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);m=p.add_mutually_exclusive_group(required=True)
    m.add_argument('--record',action='store_true');m.add_argument('--check',action='store_true')
    args=p.parse_args();r=calculate();encoded=(json.dumps(r,sort_keys=True,indent=2,allow_nan=False)+'\n').encode()
    if args.record:publish_exclusive_file(ROOT,OUTPUT,encoded)
    elif (ROOT/OUTPUT).read_bytes()!=encoded:raise ValueError('source insertion or dependencies changed')
    print(json.dumps({'status':r['status'],**{k:float(restored_upper(r['continuous_action_error_upper'][k])) for k in ('N','beta')}},indent=2))
