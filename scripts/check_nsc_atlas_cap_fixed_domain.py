#!/usr/bin/env python3
"""One retained certificate for the rounded fixed-space total fields."""
import hashlib
import json
from pathlib import Path
import time

from diagnose_nsc_atlas_boundary_rows import authenticate,Budget,plain
from check_nsc_atlas_cap_background_domain import retained_slice
from recursive_horizons.nsc_atlas_cap_fixed_transfer import transfer_initial_background
from recursive_horizons.nsc_atlas_cap_fixed_bound import fixed_metric_bound

ROOT=Path(__file__).resolve().parents[1]
SOURCES=(
    'scripts/check_nsc_atlas_cap_fixed_domain.py',
    'scripts/check_nsc_atlas_cap_background_domain.py',
    'scripts/diagnose_nsc_atlas_boundary_rows.py',
    'src/recursive_horizons/nsc_atlas_cap_fixed_basis.py',
    'src/recursive_horizons/nsc_atlas_cap_fixed_transfer.py',
    'src/recursive_horizons/nsc_atlas_cap_fixed_bound.py',
    'src/recursive_horizons/nsc_atlas_cap_background_bound.py',
    'src/recursive_horizons/nsc_atlas_cap_background_jets.py',
)


def main():
    budget=Budget(maximum_evaluations=1_000_000,cpu_seconds=15.)
    phase='authenticate'
    try:
        current,_,_,identity=authenticate(ROOT,budget)
        background,selection=retained_slice(current)
        phase='regenerate initial transfer'
        transfer=transfer_initial_background(background,selection['T'],selection['Z'],budget=budget)
        phase='fixed total-field metric certificate'
        basis=transfer['basis']
        certificate=fixed_metric_bound(basis,transfer['chi'],transfer['canonical_density'],
            transfer['wall_coordinate'],alpha=selection['alpha'],budget=budget)
        budget.check()
    except Exception:
        print(json.dumps(dict(failed_phase=phase,charged_values=budget.evaluations,
                              cpu_seconds=time.process_time()-budget.started),sort_keys=True))
        raise
    selection=dict(selection)
    selection['scope']='rounded fixed-basis transfer of retained zero-case fields at recorded finite alpha; initial constraint domain, no evolved finite state'
    report=dict(schema='nsc-atlas-cap-fixed-domain-v1',selection=selection,certificate=certificate,
        input={key:identity[key] for key in ('path','sha256','manifest_sha256','authenticated_input_files')},
        representation=dict(cells=basis.cells,field_dimension=basis.dimension),
        source_bindings={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES},
        execution=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
                       transfers=1,certificates=1,gram_factors=0,radial_solves=0,time_steps=0))
    print(json.dumps(plain(report),indent=2,sort_keys=True,allow_nan=False))


if __name__=='__main__': main()
