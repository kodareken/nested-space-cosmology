#!/usr/bin/env python3
"""One retained fixed-space constrained mass, gradient and wall comparison."""
import hashlib
import json
from pathlib import Path
import time

from diagnose_nsc_atlas_boundary_rows import authenticate,Budget,plain
from check_nsc_atlas_cap_background_domain import retained_slice
from recursive_horizons.nsc_atlas_cap_fixed_transfer import transfer_initial_background
from recursive_horizons.nsc_atlas_cap_fixed_reference import fixed_cap_reference

ROOT=Path(__file__).resolve().parents[1]
SOURCES=(
    'scripts/check_nsc_atlas_cap_fixed_reference.py',
    'scripts/check_nsc_atlas_cap_background_domain.py',
    'scripts/diagnose_nsc_atlas_boundary_rows.py',
    'src/recursive_horizons/nsc_atlas_cap_fixed_reference.py',
    'src/recursive_horizons/nsc_atlas_cap_fixed_basis.py',
    'src/recursive_horizons/nsc_atlas_cap_fixed_transfer.py',
    'src/recursive_horizons/nsc_atlas_cap_fixed_bound.py',
    'src/recursive_horizons/nsc_atlas_cap_background_bound.py',
    'src/recursive_horizons/nsc_atlas_cap_background_jets.py',
    'src/recursive_horizons/nsc_atlas_cap_variational_reference.py',
    'src/recursive_horizons/nsc_atlas_cap_radial_reference.py',
    'src/recursive_horizons/nsc_atlas_cap_geometric_reference.py',
)


def main():
    budget=Budget(maximum_evaluations=25_000_000,cpu_seconds=120.)
    phase='authenticate'
    try:
        current,_,_,identity=authenticate(ROOT,budget)
        background,selection=retained_slice(current)
        phase='initial total-field transfer'
        transfer=transfer_initial_background(background,selection['T'],selection['Z'],budget=budget)
        field=background.reduced_jets(selection['T'],selection['Z'],budget=budget)
        incoming=(field['p']+field['h'])/field['r']
        phase='fixed-space constrained mass and weak RHS'
        value=fixed_cap_reference(transfer['basis'],transfer['chi'],transfer['canonical_density'],
            transfer['wall_coordinate'],alpha=selection['alpha'],T_x=1.,radius_rate=0.,
            incoming_w_plus=incoming,budget=budget)
        phase='numerical consistency checks'
        gap=value['projected_gradient_discrepancy']
        if (abs(value['weak_balance_defect'])>2e-8 or gap['chi_rate_L2']>2e-8 or
                gap['canonical_density_rate_L2']>2e-8 or abs(gap['wall_velocity'])>2e-8):
            raise ValueError('fixed-space weak balance or projected-gradient consistency is unresolved')
        budget.check()
    except Exception as error:
        print(json.dumps(dict(failed_phase=phase,message=str(error),charged_values=budget.evaluations,
                              cpu_seconds=time.process_time()-budget.started),sort_keys=True))
        raise
    selection=dict(selection)
    selection['scope']='rounded retained total fields at finite alpha; one constrained initial slice, not evolved cap-wall solution'
    keys=('mass','lapse','radius','dimension','cells','shape_derivative_Z','gradient_discrepancy',
          'gradient_absolute_tolerance','projected_gradient_discrepancy','mass_rate','physical_weak_flux',
          'moving_density','weak_balance_defect','incoming_flux_bound','wall','centre','metric_certificate',
          'scaled_gram_condition','sampled_minimum_A','requested_tolerances','adjoint_order','field_sites',
          'coefficient_evaluations','wall_coordinate_roundoff','numerical_error_enclosure','domain')
    report=dict(schema='nsc-atlas-cap-fixed-reference-v1',selection=selection,incoming_w_plus=incoming,
        result={key:value[key] for key in keys},
        input={key:identity[key] for key in ('path','sha256','manifest_sha256','authenticated_input_files')},
        source_bindings={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES},
        execution=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
                       transfers=1,radial_references=1,time_steps=0,trajectory_replays=0),
        limits=dict(charged_values=25_000_000,cpu_seconds=120.))
    print(json.dumps(plain(report),indent=2,sort_keys=True,allow_nan=False))


if __name__=='__main__': main()
