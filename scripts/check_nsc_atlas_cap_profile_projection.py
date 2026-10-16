#!/usr/bin/env python3
"""Print one compact projection diagnosis for the frozen zero-case cut128."""
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from diagnose_nsc_atlas_boundary_rows import authenticate, Budget, plain, encoded
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_cap_profile_projection import project_characteristic_profile


ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('scripts/check_nsc_atlas_cap_profile_projection.py',
           'src/recursive_horizons/nsc_atlas_cap_profile_projection.py',
           'src/recursive_horizons/nsc_atlas_cap_polynomial_basis.py',
           'src/recursive_horizons/nsc_atlas_characteristic_background.py',
           'tests/test_nsc_atlas_cap_profile_projection.py')


def main():
    budget = Budget(maximum_evaluations=4_000_000,cpu_seconds=30.)
    current,_,_,identity = authenticate(ROOT,budget)
    case = current['cases'][0]
    if case['amplitude'] != 'zero' or case['alpha'] != 0. or len(case['snapshots']) != 129:
        raise ValueError('frozen zero-case terminal input differs')
    sample = case['snapshots'][128]
    wall = sample['physical']['wall']
    data = {key:current['profile'][key] for key in ('labels','values','first','second')}
    profile = CharacteristicProfile(*(data[key] for key in data))
    result = project_characteristic_profile(profile,wall['T'],float(np.arctanh(wall['r'])),
                                            count=16,budget=budget)
    budget.check()
    report = dict(schema='nsc-atlas-cap-profile-projection-diagnosis-v1',
        input=dict(path=identity['path'],sha256=identity['sha256'],
                   manifest_sha256=identity['manifest_sha256'],
                   authenticated_input_files=identity['authenticated_input_files'],
                   bindings_sha256=hashlib.sha256(encoded(identity['source_bindings'])).hexdigest(),
                   profile_sha256=hashlib.sha256(encoded(data)).hexdigest(),
                   cut=128,x=sample['x'],case='zero'),
        source_bindings={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES},
        result=result,
        relative_norms={key:(value/result['reference_norms'][key] if result['reference_norms'][key] else None)
                        for key,value in result['residual_norms'].items()},
        execution=dict(evaluated_points=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
                       projections=1,profile_constructions=1,radial_solves=0,time_steps=0,trajectory_replays=0),
        limits=dict(evaluated_points=4_000_000,cpu_seconds=30.))
    print(json.dumps(plain(report),indent=2,sort_keys=True,allow_nan=False))


if __name__ == '__main__':
    main()
