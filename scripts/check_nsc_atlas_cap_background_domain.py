#!/usr/bin/env python3
"""Print one local metric certificate for the retained narrow background."""
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from flint import ctx

from diagnose_nsc_atlas_boundary_rows import authenticate, Budget, plain, encoded
from recursive_horizons.nsc_atlas_cap_background_jets import RetainedCapBackground
from recursive_horizons.nsc_atlas_cap_background_bound import background_metric_bound
from recursive_horizons.nsc_atlas_cap_variational_reference import _exact_ball


ROOT=Path(__file__).resolve().parents[1]
SOURCES=('scripts/check_nsc_atlas_cap_background_domain.py',
         'src/recursive_horizons/nsc_atlas_cap_background_bound.py',
         'src/recursive_horizons/nsc_atlas_cap_background_jets.py')


def retained_slice(current):
    """One declared inward rounding of the existing zero-case terminal slice."""
    zero,finite=current['cases']
    if (zero['amplitude']!='zero' or zero['alpha']!=0. or finite['amplitude']!='finite'
            or finite['alpha']<=0 or len(zero['snapshots'])!=129):
        raise ValueError('frozen input cases differ')
    wall=zero['snapshots'][128]['physical']['wall']
    T,R=wall['T'],wall['r']
    nominal=float(np.arctanh(R))
    Z=nominal
    data={key:current['profile'][key] for key in ('labels','values','first','second')}
    with ctx.workprec(160):
        upper=_exact_ball(data['labels'][-1])
        if not _exact_ball(T)+_exact_ball(Z)<=upper:
            Z=float(np.nextafter(Z,-np.inf))
        if not _exact_ball(T)+_exact_ball(Z)<=upper:
            raise ValueError('one inward radius ulp did not restore the supplied profile domain')
    return RetainedCapBackground(*(data[key] for key in data)),dict(
        T=T,Z=Z,nominal_Z=nominal,nominal_radius=R,radius=float(np.tanh(Z)),
        alpha=float(finite['alpha']),cut=128,
        profile_sha256=hashlib.sha256(encoded(data)).hexdigest(),
        scope='retained zero-case profile, zero polynomial departures, recorded finite alpha; no evolved finite state')


def main():
    budget=Budget(maximum_evaluations=100000,cpu_seconds=10.)
    current,_,_,identity=authenticate(ROOT,budget)
    background,selection=retained_slice(current)
    result=background_metric_bound(background,selection['T'],selection['Z'],np.zeros(16),np.zeros(16),
                                   alpha=selection['alpha'],budget=budget)
    budget.check()
    report=dict(schema='nsc-atlas-cap-background-domain-v1',selection=selection,certificate=result,
        input={key:identity[key] for key in ('path','sha256','manifest_sha256','authenticated_input_files')},
        source_bindings={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES},
        execution=dict(evaluated_points=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
                       certificates=1,radial_solves=0,time_steps=0,trajectory_replays=0))
    print(json.dumps(plain(report),indent=2,sort_keys=True,allow_nan=False))


if __name__=='__main__':
    main()
