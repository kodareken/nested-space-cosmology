#!/usr/bin/env python3
"""Two fixed retained-pulse mass references; no departure evolution."""
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from diagnose_nsc_atlas_boundary_rows import authenticate, Budget, plain
from check_nsc_atlas_cap_background_domain import retained_slice
from recursive_horizons.nsc_atlas_cap_polynomial_reference import polynomial_cap_reference


ROOT=Path(__file__).resolve().parents[1]
SOURCES=('scripts/check_nsc_atlas_cap_background_reference.py',
         'scripts/check_nsc_atlas_cap_background_domain.py',
         'src/recursive_horizons/nsc_atlas_cap_polynomial_reference.py',
         'src/recursive_horizons/nsc_atlas_cap_variational_reference.py',
         'src/recursive_horizons/nsc_atlas_cap_background_jets.py',
         'src/recursive_horizons/nsc_atlas_cap_background_bound.py')


def main():
    budget=Budget(maximum_evaluations=20_000_000,cpu_seconds=120.)
    phase='authentication'
    try:
        current,_,_,identity=authenticate(ROOT,budget)
        background,selection=retained_slice(current)
        T,Z=selection['T'],selection['Z']
        field=background.reduced_jets(T,Z,budget=budget)
        incoming=(field['p']+field['h'])/field['r']
        rows=[]
        for name,alpha in [('geometric_limit',0.),('recorded_finite_alpha',selection['alpha'])]:
            phase=name
            row=polynomial_cap_reference(Z,np.zeros(16),np.zeros(16),alpha=alpha,T_x=1.,
                radius_rate=0.,incoming_w_plus=incoming,basis='legendre',background=background,
                background_time=T,budget=budget)
            rate_max=max(float(np.max(abs(v))) for v in row['rates'].values())
            if abs(row['remainder_identity_defect'])>2e-10:
                raise ValueError('explicit-time remainder identity failed')
            if alpha==0:
                expected=field['k']*field['A0']*field['h']/2
                if (rate_max>2e-7 or abs(row['gradient']['T']-expected)>2e-7
                        or abs(row['explicit_time_remainder'])>2e-7):
                    raise ValueError('retained geometric background limit is unresolved')
            keys=('mass','lapse','gradient_discrepancy','shape_discrepancy','explicit_time_remainder',
                  'balance_defect','remainder_identity_defect','interface_balance_defect',
                  'physical_weak_flux','coefficient_evaluations','field_sites')
            rows.append(dict(case=name,alpha=alpha,**{key:row[key] for key in keys},
                explicit_time_gradient=row['gradient']['T'],maximum_departure_rate=rate_max,
                metric_floor=row['metric_certificate']['A_lower_float_down'],
                wall={key:row['wall'][key] for key in ('momentum_discrepancy','normal_discrepancy',
                                                       'scalar_tangent_discrepancy')},
                numerical_error_enclosure=False))
        budget.check()
        report=dict(schema='nsc-atlas-cap-background-reference-diagnosis-v1',selection=selection,
            incoming_w_plus=incoming,cases=rows,
            input={key:identity[key] for key in ('path','sha256','manifest_sha256','authenticated_input_files')},
            source_bindings={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES},
            execution=dict(evaluated_points=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
                           references=2,time_steps=0,trajectory_replays=0),
            limits=dict(evaluated_points=20_000_000,cpu_seconds=120.))
        print(json.dumps(plain(report),indent=2,sort_keys=True,allow_nan=False))
    except Exception as error:
        print(json.dumps(dict(status='failed',phase=phase,message=str(error),
            evaluated_points=budget.evaluations,cpu_seconds=time.process_time()-budget.started)))
        raise


if __name__=='__main__':
    main()
