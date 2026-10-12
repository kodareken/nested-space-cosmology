#!/usr/bin/env python3
"""Record/replay a directed original-history cubic coefficient at one point."""
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
from recursive_horizons.nsc_ks_cubic_uv_current import reduced_massless_transport_identities
from recursive_horizons.nsc_ks_cubic_uv_enclosure import CubicGeometry,enclose_characteristic
from recursive_horizons.nsc_ks_current_uv_transport import (
    HISTORY_RECORD,ARCHIVE_RHO_UP,CURRENT_HISTORY_IDENTITY,smallest_original_pair,
)
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

OUTPUT='results/development/nsc-ks-cubic-uv-enclosure-v1.json'
INPUT='results/development/nsc-ks-cubic-uv-current-v1.json'
OWNERS=('scripts/derive_nsc_ks_cubic_uv_enclosure.py',
        'src/recursive_horizons/nsc_ks_cubic_uv_enclosure.py',
        'tests/test_nsc_ks_cubic_uv_enclosure.py','docs/nsc-ks-cubic-uv-enclosure.md')


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def interval(value):
    return {'lower':exact_upper(value.lower()),'upper':exact_upper(value.upper())}


def calculate():
    input_record=json.loads((ROOT/INPUT).read_text())
    for group in ('source_hashes','input_hashes'):
        for path,expected in input_record[group].items():
            if digest(path)!=expected:
                raise ValueError('cubic input dependency changed: '+path)
    identities=dict(reduced_massless_transport_identities())
    if input_record['reduced_transport']!=identities:
        raise ValueError('current coefficient reduction differs from its proved identities')
    history=json.loads((ROOT/HISTORY_RECORD).read_text())
    family=LocalIncomingFamily(np.asarray(history['history']['coefficients']))
    if profile_identity(family)!=CURRENT_HISTORY_IDENTITY:
        raise ValueError('original current history required')
    pair=dict(smallest_original_pair(ROOT))
    if pair['group']!=1 or pair['mass']!=0 or pair['absolute_angular']!=np.sqrt(5.):
        raise ValueError('original massless group-1 angular pair required')
    rows=[]
    with ctx.workprec(160):
        ell=arb(pair['absolute_angular']).union(arb(5).sqrt())
        model=CubicGeometry(family,ell,bits=160)
        z=arb(family.center)
        N,beta=arb(0),arb(0)
        for cells,sign in ((256,1),(1024,1),(1024,-1)):
            bound=enclose_characteristic(model,z,sign,ARCHIVE_RHO_UP,cells=cells,order=8)
            rows.append({'cells':cells,'order':8,'sign':sign,
                         'accepted_cells':bound['accepted_cells'],
                         'values':{key:interval(bound[key]) for key in
                                   ('q_z','q_zz','J3','N_bracket3')}})
            if cells==1024:
                weight=2*arb(pair['multiplicity_per_signed_family'])/(2*arb.pi())
                N-=weight*bound['N_bracket3'];beta+=weight*bound['J3']
        return {'schema':'NSC-KS-CUBIC-UV-ENCLOSURE-v1',
                'status':'ENCLOSED: formal massless cubic coefficient at one incoming point; physical UV budget OPEN',
                'profile_identity':CURRENT_HISTORY_IDENTITY,'pair':pair,
                'bits':160,'rho_up_hex':float(ARCHIVE_RHO_UP).hex(),
                'incoming_point_hex':float(family.center).hex(),'angular_enclosure':interval(ell),
                'cases':rows,'paired_action_cubic':{'N':interval(N),'beta':interval(beta)},
                'integrated_leading_term_at_original_cutoff':{
                    'N':interval(N/(2*arb(pair['cutoff'])**2)),
                    'beta':interval(beta/(2*arb(pair['cutoff'])**2))},
                'leading_term_is_complete_uv_tail':False,
                'angular_evenness_used':True,'source_signs':[1,-1],
                'measure':'mu/(2 pi), angular pair twice, source signs separately; N minus and beta plus',
                'initial_difference_values':[0,0,0],
                'initial_physical_normalization_claimed_zero':False,
                'coverage':{'group':1,'incoming_points':1,'entire_incoming_interval':False,
                            'all_source_families':False,'massive_channels':False},
                'uniform_C4_on_I':None,'higher_uv_remainder_C_M':None,
                'physical_source_error':None,'physical_UV_budget_component':None,
                'physical_local_gate':'OPEN','source_inventory_unchanged':True,
                'source_hashes':implementation_hashes(ROOT,owners=OWNERS),
                'input_hashes':{INPUT:digest(INPUT),HISTORY_RECORD:digest(HISTORY_RECORD)}}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record',action='store_true');mode.add_argument('--check',action='store_true')
    args=parser.parse_args();result=calculate()
    # JSON round trip gives arrays for the original pair's tuple metadata.
    encoded=(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n').encode()
    if args.record:
        publish_exclusive_file(ROOT,OUTPUT,encoded)
    elif (ROOT/OUTPUT).read_bytes()!=encoded:
        raise ValueError('directed cubic coefficient or dependencies changed')
    print(json.dumps({'status':result['status'],
        'paired_action_cubic':{key:{end:float(restored_upper(value)) for end,value in row.items()}
                               for key,row in result['paired_action_cubic'].items()},
        'physical_local_gate':'OPEN'},indent=2))
