#!/usr/bin/env python3
"""Record/replay original-channel vacuum initial-column remainder bounds."""
import argparse,json,sys
from pathlib import Path
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper,restored_upper
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_vacuum_source_remainder import VacuumSourceExpansion
from recursive_horizons.nsc_vacuum_spinor_remainder import spinor_remainder,truncated_upstream_envelope
OUTPUT='results/development/nsc-vacuum-spinor-remainder-v1.json'
OWNERS=('scripts/derive_nsc_vacuum_spinor_remainder.py',
 'src/recursive_horizons/nsc_vacuum_spinor_remainder.py',
 'tests/test_nsc_vacuum_spinor_remainder.py','docs/nsc-vacuum-spinor-remainder.md')


def interval(v):return {'lower':exact_upper(v.lower()),'upper':exact_upper(v.upper())}


def calculate():
    archive=RetainedUpstreamArchive(ROOT);cases=[]
    with ctx.workprec(192):
        for group in (1,14):
            channel=archive.meta['channels'][group]
            for angular_sign in (-1,1):
                entries=archive.family_entries((group,angular_sign))
                batch=next(b for b,_ in entries if b.energy_sign>0)
                if any(b.rho_up!=batch.rho_up for b,_ in entries):
                    raise ValueError('common upstream slice required')
                expected_mass=0. if group==1 else float(arb.pi()/2)
                if batch.mass!=expected_mass or batch.angular!=angular_sign*float(arb(5).sqrt()):
                    raise ValueError('original pilot channel differs')
                m=arb(batch.mass).union(arb(0) if group==1 else arb.pi()/2)
                ell=arb(batch.angular).union(angular_sign*arb(5).sqrt())
                model=VacuumSourceExpansion(archive.meta['config']['horizon_rho'],m,ell,order=8)
                for cutoff in (160,320):
                    full=spinor_remainder(model,batch.rho_up,cutoff,cells=64)
                    reduced=truncated_upstream_envelope(full,4,1)
                    cases.append({'group':group,'angular_sign':angular_sign,
                     'mass':interval(m),'angular':interval(ell),'rho_up_hex':batch.rho_up.hex(),
                     'cutoff':cutoff,'order':8,'retained_degree':4,
                     'major_real_coefficients':[interval(v) for v in reduced['coefficients'][0]],
                     'minor_real_coefficients':[interval(v) for v in reduced['coefficients'][1]],
                     'minor_imag_coefficients':[interval(v) for v in reduced['coefficients'][2]],
                     'chart_major_lower':exact_upper(full['chart_major_lower']),
                     'eighth_order_constant':exact_upper(full['spinor_remainder_constant']),
                     'fourth_order_scaled_pointwise_error_upper':exact_upper(reduced['pointwise_scaled_error_upper']),
                     'negative_vacuum_partner':'swap minor and real major; complement of signed covariance',
                     'spatial_norm_domain_bound':None})
    return {'schema':'NSC-VACUUM-SPINOR-REMAINDER-v1',
     'status':'initial vacuum-column bound only; changed-history UV and local gate OPEN',
     'cases':cases,'bits':192,'remainder_cells':64,
     'phase_gauge':'major real positive at upstream; phase constant during later spacetime transport',
     'requires_returned_upstream_coefficients':True,'finite_occupation_source_error_included':False,
     'all_source_families':False,'current_history_L0_A4_H2_integrals':None,
     'changed_history_C_M':None,'physical_local_gate':'OPEN',
     'source_hashes':implementation_hashes(ROOT,owners=OWNERS),
     'input_hashes':dict(archive.input_hashes)}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--record',action='store_true');modes.add_argument('--check',action='store_true')
    args=parser.parse_args();result=calculate()
    encoded=(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n').encode()
    if args.record:publish_exclusive_file(ROOT,OUTPUT,encoded)
    elif (ROOT/OUTPUT).read_bytes()!=encoded:raise ValueError('spinor proof or dependencies changed')
    print(json.dumps({'status':result['status'],'cases':len(result['cases'])},indent=2))
