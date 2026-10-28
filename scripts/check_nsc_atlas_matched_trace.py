"""One manufactured matched incoming characteristic on P001's coarse space.

No production evolution, source history, new mesh or estimate of its maxima.
This tests whether the adapter reopens a strong trace already made compatible.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from threadpoolctl import threadpool_limits
from run_nsc_atlas_fixed_endpoint import cap_knots,ORIGINAL_ALPHA
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_cap_fixed_reference import fixed_cap_reference
from recursive_horizons.nsc_atlas_cap_fixed_hybrid import fixed_cap_hybrid_reference
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget

ROOT=Path(__file__).resolve().parents[1]


def calculate():
    budget=Budget(maximum_evaluations=4_000_000,cpu_seconds=30.)
    basis=FixedCapBasis(cap_knots(0),budget=budget)
    c=np.zeros(basis.dimension);d=c.copy()
    c[3]=.002;d[4]=-.003
    Z=F(math.atanh(.45));seen={}
    def matched(packet):
        seen.update(packet)
        return dict(T_x=1.,radius_rate=0.,incoming_w_plus=(packet['p']+packet['h'])/packet['radius'])
    with threadpool_limits(limits=1):
        direct=fixed_cap_reference(basis,c,d,Z,alpha=ORIGINAL_ALPHA,wall_coupler=matched,budget=budget)
        w=direct['wall'];incoming=(w['p_trace']+w['h_trace'])/direct['radius']
        def exterior(common):
            R,cc=common['r'],common['c'];u=cc/R
            D=1+ORIGINAL_ALPHA*u*u/6
            rV=.11;H=2*rV+2*D**1.5*R
            return dict(r_V=rV,c_V=H*(u+R*math.sqrt(D)*incoming)/2,log_C_V=.003)
        coupled=fixed_cap_hybrid_reference(basis,c,d,Z,alpha=ORIGINAL_ALPHA,exterior_wall=exterior,
            C_wall=1.,eta_plus=0.,eta_plus_x=.01,budget=budget)
    wall=coupled['cap']['wall'];res=coupled['residual']
    summary=dict(schema='nsc-atlas-matched-trace-v1',
        physical_input=dict(alpha=ORIGINAL_ALPHA,wall_radius=direct['radius'],cap_knots=[str(k) for k in cap_knots(0)],
                            departure='One fixed outer-cell bubble in chi and canonical density, .002 and -.003.'),
        direct=dict(incoming_discrepancy=w['incoming_discrepancy'],velocity_discrepancy=w['velocity_discrepancy'],
                    weak_balance_defect=direct['weak_balance_defect']),
        hybrid=dict(incoming_discrepancy=wall['incoming_discrepancy'],velocity_discrepancy=wall['velocity_discrepancy'],
                    scalar_tangent=res['scalar_tangent'],scalar_velocity_channel=res['scalar_velocity_channel'],
                    scalar_incoming_channel=res['scalar_incoming_channel'],scalar_channel_identity_defect=res['scalar_channel_identity_defect'],
                    physical_mass_residual=res['physical_temporal_mass_actual']),
        conclusion='The adapter preserves a deliberately matched strong incoming sum. Remaining velocity/metric-projection error is separate. Production P001 did not enforce this strong compatibility and remains unconverged.',
        scope='Manufactured instantaneous algebraic compatibility test, not a repair or replay of the P001 trajectory.',
        execution=dict(radial_solves=2,time_steps=0,new_meshes=0,charged_values=budget.evaluations,
                       process_CPU_seconds=time.process_time()-budget.started))
    if max(abs(w['incoming_discrepancy']),abs(wall['incoming_discrepancy']),abs(res['scalar_channel_identity_defect']))>2e-12:
        raise ValueError('matched characteristic not preserved')
    summary['source_sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (
        'scripts/check_nsc_atlas_matched_trace.py','src/recursive_horizons/nsc_atlas_cap_fixed_reference.py',
        'src/recursive_horizons/nsc_atlas_cap_fixed_hybrid.py')}
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();path=(ROOT/args.output).resolve()
    if not path.is_relative_to(ROOT) or path.exists():raise ValueError('new repository output required')
    value=calculate()
    path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(value))


if __name__=='__main__':main()
