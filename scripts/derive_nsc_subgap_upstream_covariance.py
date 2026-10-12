#!/usr/bin/env python3
"""Enclose one unchanged archived upstream source and its signed partner."""
import argparse
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import sys

import numpy as np
from flint import acb,arb,ctx

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper,restored_upper
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_metric_horizon_frame import metric_horizon_frame
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_subgap_source_covariance import (
    BlochSource,initial_bloch,capture_bloch,validate_bloch,original_covariance_error,
)
from recursive_horizons.nsc_subgap_upstream_covariance import (
    negative_subgap_covariance_error,numerical_gram_error,
)

OUTPUT='results/development/nsc-subgap-upstream-covariance-v1.json'
PAYLOAD='results/development/artifacts/nsc-subgap-upstream-covariance-v1.npz'
SOURCE='results/development/nsc-subgap-source-covariance-v1.json'
HORIZON='results/development/nsc-metric-horizon-frame-v1.json'
OWNERS=('scripts/derive_nsc_subgap_upstream_covariance.py',
    'src/recursive_horizons/nsc_subgap_upstream_covariance.py',
    'tests/test_nsc_subgap_upstream_covariance.py','docs/nsc-subgap-upstream-covariance.md')


def digest(path):return sha256((ROOT/path).read_bytes()).hexdigest()


def ball(value):
    return restored_upper(value['lower']).union(restored_upper(value['upper']))


def authenticated_parent(path):
    record=json.loads((ROOT/path).read_text())
    for group in ('source_hashes','input_hashes'):
        for name,expected in record[group].items():
            if digest(name)!=expected:raise ValueError('parent dependency changed: '+name)
    payload=record.get('payload')
    if payload and digest(payload['path'])!=payload['sha256']:
        raise ValueError('parent payload changed')
    return record


def calculate(saved=None):
    source=authenticated_parent(SOURCE)
    horizon=authenticated_parent(HORIZON)
    archive=RetainedUpstreamArchive(ROOT)
    selected=[batch for batch,_ in archive.family_entries((14,1))
        if batch.original_panel==source['source_panel'] and batch.rows==(0,16)]
    if len(selected)!=2 or {b.energy_sign for b in selected}!={-1,1}:
        raise ValueError('original signed upstream pair missing')
    positive=next(b for b in selected if b.energy_sign==1)
    negative=next(b for b in selected if b.energy_sign==-1)
    E=float.fromhex(source['energy_hex'])
    if (source['source_row']!=0 or horizon['energy_hex']!=E.hex()
            or source['source_panel']!=horizon['source_panel']
            or not np.all(positive.source.energies[:3]==E)
            or not np.all(negative.source.energies[:3]==-E)
            or negative.angular!=-positive.angular
            or positive.rho_up!=negative.rho_up):
        raise ValueError('source/horizon/upstream signed binding changed')
    Fp=positive.initial_columns[:,:3];Fn=negative.initial_columns[:,:3]
    Cp=positive.source.covariance[:3,:3];Cn=negative.source.covariance[:3,:3]
    with ctx.workprec(192):
        mass,angular=ball(horizon['mass']),ball(horizon['angular'])
        if (not mass.contains(arb(positive.mass)) or negative.mass!=positive.mass
                or not angular.contains(arb(positive.angular))
                or float(archive.meta['config']['surface_gravity']).hex()!=source['source_kappa_hex']):
            raise ValueError('horizon/source/upstream mass, angular or occupation binding changed')
        frame=metric_horizon_frame(E,mass,angular,archive.meta['config']['horizon_rho'],bits=192)
        reflection=acb(ball(horizon['affine_reflection_direct_frame']['real']),
                       ball(horizon['affine_reflection_direct_frame']['imag']))
        start=float.fromhex(source['start_log_delta_hex'])
        target=(frame.q-arb.pi()/2-arb(positive.rho_up).atan()).log()
        initial=initial_bloch(frame,reflection,archive.meta['config']['surface_gravity'],start)
        model=BlochSource(frame.q,E,mass,angular,bits=192)
        if saved is None:
            trace,nfev=capture_bloch(model,initial,start,float(target.mid()),max_step=.025)
        else:
            saved_record,saved_raw=saved
            if sha256(saved_raw).hexdigest()!=saved_record['payload']['sha256']:
                raise ValueError('upstream proof payload hash changed')
            with np.load(BytesIO(saved_raw),allow_pickle=False) as witness:
                if set(witness.files)!={'trace','positive_columns','negative_columns',
                                       'positive_covariance','negative_covariance'}:
                    raise ValueError('complete original upstream proof payload required')
                for name,original in (('positive_columns',Fp),('negative_columns',Fn),
                                      ('positive_covariance',Cp),('negative_covariance',Cn)):
                    if not np.array_equal(witness[name],original):
                        raise ValueError('archived upstream source changed in witness: '+name)
                trace=witness['trace'].copy()
            nfev=saved_record['nfev']
        if trace.ndim!=2 or len(trace)==0 or trace.shape[1]!=26 or trace[0,0]!=start:
            raise ValueError('trajectory must start at its enclosed horizon initializer')
        proof=validate_bloch(model,initial,trace,target,degree=12)
        ep=original_covariance_error(Fp,Cp,proof)
        en=negative_subgap_covariance_error(Fn,Cn,proof)
        raw=deterministic_npz_bytes({'trace':trace,'positive_columns':Fp,'negative_columns':Fn,
                                    'positive_covariance':Cp,'negative_covariance':Cn})
        result={'schema':'NSC-SUBGAP-UPSTREAM-COVARIANCE-v1',
            'status':'ENCLOSED: one original upstream subgap row and its signed partner; full source budget OPEN',
            'source_panel':source['source_panel'],'source_row':0,'group':14,
            'positive_energy_hex':E.hex(),'rho_up_hex':float(positive.rho_up).hex(),
            'source_signs':[1,-1],'angular_signs':[1,-1],
            'positive_preparation_digest':positive.preparation_digest,
            'negative_preparation_digest':negative.preparation_digest,
            'original_batch_rows':list(positive.rows),'covered_rows':[0,1],
            'positive_source_digest':positive.source.digest,'negative_source_digest':negative.source.digest,
            'original_quadrature_weight_hex':float(archive._arrays[source['source_panel']+'/weights'][0]).hex(),
            'covariance_error_unweighted':True,'bits':192,'defect_degree':12,
            'max_step_hex':float(.025).hex(),'cells':proof['cells'],'nfev':nfev,
            'positive_covariance_operator_error_upper':exact_upper(ep),
            'negative_covariance_operator_error_upper':exact_upper(en),
            'numerical_positive_gram_frobenius_error_upper':exact_upper(numerical_gram_error(Fp)),
            'initial_bloch_error_upper':exact_upper(proof['initial_error']),
            'normalized_defect_integral_upper':exact_upper(proof['defect_integral']),
            'endpoint_bridge_upper':exact_upper(proof['endpoint_bridge']),
            'negative_error_copied_from_positive':False,
            'replay_uses_saved_witness':True,
            'exact_gram_identity_scope':'homogeneous complete subgap preparation only',
            'archived_upstream_columns_replaced':False,'preparation_ode_rerun':False,
            'changed_history_evaluated':False,'all_source_families':False,
            'source_quadrature_error_included':False,'physical_upstream_budget_component':None,
            'physical_local_gate':'OPEN',
            'payload':{'path':PAYLOAD,'sha256':sha256(raw).hexdigest(),'bytes':len(raw)},
            'source_hashes':implementation_hashes(ROOT,owners=OWNERS),
            'input_hashes':{**archive.input_hashes,SOURCE:digest(SOURCE),HORIZON:digest(HORIZON)}}
    return result,raw


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record',action='store_true');mode.add_argument('--check',action='store_true')
    args=parser.parse_args()
    saved=None if args.record else (json.loads((ROOT/OUTPUT).read_text()),(ROOT/PAYLOAD).read_bytes())
    record,raw=calculate(saved=saved)
    encoded=(json.dumps(record,sort_keys=True,indent=2,allow_nan=False)+'\n').encode()
    if args.record:
        if (ROOT/OUTPUT).exists() or (ROOT/PAYLOAD).exists():
            raise FileExistsError('upstream record already exists; use --check')
        publish_exclusive_file(ROOT,PAYLOAD,raw);publish_exclusive_file(ROOT,OUTPUT,encoded)
    elif (ROOT/OUTPUT).read_bytes()!=encoded or (ROOT/PAYLOAD).read_bytes()!=raw:
        raise ValueError('upstream covariance proof, payload or dependencies changed')
    print(json.dumps({'status':record['status'],'cells':record['cells'],
        'positive_error_upper':float(restored_upper(record['positive_covariance_operator_error_upper'])),
        'negative_error_upper':float(restored_upper(record['negative_covariance_operator_error_upper'])),
        'gram_error_upper':float(restored_upper(record['numerical_positive_gram_frobenius_error_upper'])),
        'physical_upstream_budget_component':None},indent=2))
