#!/usr/bin/env python3
"""Record and replay an enclosed correction at the first retained middle energy."""
import argparse
from hashlib import sha256
from io import BytesIO
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
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_subgap_source_covariance import original_covariance_distance_bounds
from recursive_horizons.nsc_vacuum_source_remainder import VacuumSourceExpansion
from recursive_horizons.nsc_vacuum_source_correction import VacuumCorrection,capture_correction,validate_correction

OUTPUT='results/development/nsc-vacuum-source-correction-v1.json'
PAYLOAD='results/development/artifacts/nsc-vacuum-source-correction-v1.npz'
OWNERS=('scripts/derive_nsc_vacuum_source_correction.py',
    'src/recursive_horizons/nsc_vacuum_source_correction.py',
    'tests/test_nsc_vacuum_source_correction.py','docs/nsc-vacuum-source-correction.md')
START=-18.


def interval(value):
    return {'lower':exact_upper(value.lower()),'upper':exact_upper(value.upper())}


def calculate(saved=None):
    archive=RetainedUpstreamArchive(ROOT);entries=archive.family_entries((14,1))
    choices=[(float(E),b,j) for b,_ in entries if b.energy_sign>0
        for j,E in enumerate(b.source.energies[::3]) if E>=16]
    E,positive,index=min(choices,key=lambda row:row[0])
    negative=next(b for b,_ in entries if b.energy_sign<0
        and b.original_panel==positive.original_panel and b.rows==positive.rows)
    if (positive.original_panel!='group14/mid24_1' or positive.rows[0]+index!=0
            or not np.all(negative.source.energies[3*index:3*index+3]==-E)):
        raise ValueError('the original first middle-row signed pair is required')
    config=archive.meta['config'];sl=slice(3*index,3*index+3)
    originals={'positive_columns':positive.initial_columns[:,sl],
               'negative_columns':negative.initial_columns[:,sl],
               'positive_covariance':positive.source.covariance[sl,sl],
               'negative_covariance':negative.source.covariance[sl,sl]}
    with ctx.workprec(192):
        if positive.mass!=np.pi/2 or positive.angular!=np.sqrt(5):
            raise ValueError('original group14 mass and angular label required')
        expansion=VacuumSourceExpansion(config['horizon_rho'],arb(positive.mass).union(arb.pi()/2),
            arb(positive.angular).union(arb(5).sqrt()),order=8)
        flow=VacuumCorrection(expansion,E)
        target=float(expansion.target_distance(positive.rho_up).log().mid())
        if saved is None:
            trace,nfev=capture_correction(flow,START,target,cpu_limit=30.)
        else:
            old,raw=saved
            if sha256(raw).hexdigest()!=old['payload']['sha256']:
                raise ValueError('correction payload hash changed')
            with np.load(BytesIO(raw),allow_pickle=False) as data:
                if set(data.files)!={'trace',*originals}:
                    raise ValueError('complete correction/source witness required')
                for key,value in originals.items():
                    if not np.array_equal(data[key],value):
                        raise ValueError('original source changed in correction witness: '+key)
                trace=data['trace'].copy()
            nfev=old['nfev']
        proof=validate_correction(flow,trace,START,positive.rho_up)
        # Only the finite-occupation term is taken from this helper. The
        # vacuum error is the separately validated correction-curve error.
        thermal=expansion.covariance_error(0,E,kappa=config['surface_gravity'],
            omega=config['omega'])['thermal_coherent_error']
        comparisons=[]
        for batch,key in ((positive,'positive'),(negative,'negative')):
            validation=dict(proof)
            if batch.energy_sign<0:
                n=proof['endpoint'];validation['endpoint']=(n[0],-n[1],-n[2])
            bounds=original_covariance_distance_bounds(originals[key+'_columns'],
                originals[key+'_covariance'],validation)
            comparisons.append({'energy_sign':batch.energy_sign,'angular_sign':batch.angular_sign,
                'source_digest':batch.source.digest,'preparation_digest':batch.preparation_digest,
                'physical_source_operator_error_lower':exact_upper(max(arb(0),(bounds['lower']-thermal).lower())),
                'physical_source_operator_error_upper':exact_upper((bounds['upper']+thermal).upper())})
        payload=deterministic_npz_bytes({'trace':trace,**originals})
        record={'schema':'NSC-VACUUM-SOURCE-CORRECTION-v1',
            'status':'ENCLOSED: one original middle-energy preparation and its signed partner; full gate OPEN',
            'group':14,'panel':positive.original_panel,'row':0,'energy_hex':E.hex(),
            'rho_up_hex':positive.rho_up.hex(),'start_log_delta_hex':START.hex(),
            'end_log_delta_hex':target.hex(),'expansion_order':8,'defect_degree':12,'bits':192,
            'rtol_hex':float(1e-9).hex(),'atol_hex':float(1e-21).hex(),'max_step_hex':float(.1).hex(),
            'capture_cpu_limit_seconds':30,'cells':proof['cells'],'nfev':nfev,
            'initial_bloch_error_upper':exact_upper(proof['initial_error']),
            'normalized_defect_integral_upper':exact_upper(proof['normalized_defect_integral']),
            'endpoint_bridge_upper':exact_upper(proof['endpoint_bridge']),
            'vacuum_bloch_error_upper':exact_upper(proof['bloch_error']),
            'thermal_coherent_operator_error_upper':exact_upper(thermal),
            'validated_vacuum_bloch_center':[interval(v) for v in proof['endpoint']],
            'source_comparisons':comparisons,
            'covariance_errors_unweighted':True,
            'initial_true_correction_assumed_zero':False,'replay_uses_saved_witness':True,
            'archived_source_replaced':False,'source_occupations_changed':False,
            'original_preparation_ode_rerun':False,'proof_correction_ode_evolved':True,
            'source_quadrature_error_included':False,'all_source_families':False,
            'physical_upstream_budget_component':None,'changed_history_C_M':None,
            'physical_local_gate':'OPEN',
            'payload':{'path':PAYLOAD,'sha256':sha256(payload).hexdigest(),'bytes':len(payload)},
            'source_hashes':implementation_hashes(ROOT,owners=OWNERS),
            'input_hashes':dict(archive.input_hashes)}
    return record,payload


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record',action='store_true');mode.add_argument('--check',action='store_true')
    args=parser.parse_args();started=time.process_time()
    saved=None if args.record else (json.loads((ROOT/OUTPUT).read_text()),(ROOT/PAYLOAD).read_bytes())
    record,raw=calculate(saved)
    encoded=(json.dumps(record,sort_keys=True,indent=2,allow_nan=False)+'\n').encode()
    if args.record:
        if (ROOT/OUTPUT).exists() or (ROOT/PAYLOAD).exists():raise FileExistsError('record exists; use --check')
        publish_exclusive_file(ROOT,PAYLOAD,raw);publish_exclusive_file(ROOT,OUTPUT,encoded)
    elif (ROOT/OUTPUT).read_bytes()!=encoded or (ROOT/PAYLOAD).read_bytes()!=raw:
        raise ValueError('correction proof, source or dependencies changed')
    print(json.dumps({'status':record['status'],'cells':record['cells'],'cpu_seconds':time.process_time()-started,
        'vacuum_bloch_error_upper':float(restored_upper(record['vacuum_bloch_error_upper'])),
        'source_error_intervals':[{'energy_sign':v['energy_sign'],
            'lower':float(restored_upper(v['physical_source_operator_error_lower'])),
            'upper':float(restored_upper(v['physical_source_operator_error_upper']))}
            for v in record['source_comparisons']]},indent=2))
