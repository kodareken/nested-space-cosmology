"""One fixed-action charged-parent join and regular shape-transition candidate.

No old directional trajectory is transplanted. Both pole states follow from
the same coupled Euclidean profile, and both Lorentzian regions inherit it.
"""
import argparse
from dataclasses import asdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_sphere_flux_completion import (
    combined_record,all_volume_extrema,volume_thresholds,half_integer_classification,
    wrapped_brane_couplings,symmetric_radius_branches,log_radius_potential,
)
from recursive_horizons.nsc_atlas_shape_transition import (
    solve_bounce,bounce_diagnostics,continue_pole,connecting_geometry,
)

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pins():
    files={Path(__file__).resolve()};package=(ROOT/'src/recursive_horizons').resolve()
    for module in tuple(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if name:
            p=Path(name).resolve()
            if p.suffix=='.py' and p.is_relative_to(package):files.add(p)
    return {p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)}


def plain(value):
    if isinstance(value,Fraction):return str(value)
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    raise TypeError(type(value).__name__)


def snapshot_table(sol,count=41):
    times=np.linspace(sol.t[0],sol.t[-1],count);y=sol.sol(times)
    return dict(columns=['proper_time','u_single','u_other','velocity_single','velocity_other','scale','scale_rate'],
        rows=np.column_stack([times,y.T]).tolist())


def calculation(folder):
    started=time.perf_counter();branches=symmetric_radius_branches();source=pins()
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bound=True
    for name,digest in source.items():
        old=subprocess.run(['git','show',commit+':'+name],cwd=ROOT,capture_output=True)
        bound=bound and old.returncode==0 and hashlib.sha256(old.stdout).hexdigest()==digest
    model=combined_record();handoff=model.scalar_sector.prepare(parent_flux=2.,charge=.5)
    with threadpool_limits(limits=4):
        coarse=solve_bounce(tolerance=2e-6,cutoff=1e-4,mesh=161)
        fine=solve_bounce(tolerance=5e-8,cutoff=5e-5,mesh=321,previous=coarse)
        a=bounce_diagnostics(coarse);b=bounce_diagnostics(fine)
        for row in (a,b):
            if (row['max_Euclidean_constraint']>1e-6 or row['max_common_Cauchy_relative_Hamiltonian']>.02 or
                row['max_scaled_ODE_residual']>1e-5 or row['projected_barrier_crossings']!=1):
                raise RuntimeError('coupled Euclidean profile did not pass declared numerical controls')
        interior,inside=continue_pole(fine,pole_index=0)
        exterior,outside=continue_pole(fine,pole_index=1,maximum_time=10.)
        connecting=connecting_geometry(fine)
    if source!=pins():raise RuntimeError('formation sources changed during calculation')
    raw=folder/'profiles.npz';payload={}
    for name,bounce in (('coarse',coarse),('fine',fine)):
        payload[name+'_knots']=bounce.solution.sol.x
        payload[name+'_coefficients']=bounce.solution.sol.c
        payload[name+'_length']=bounce.length
        payload[name+'_cutoff']=bounce.cutoff
    for name,sol in (('interior',interior),('exterior',exterior)):
        times=np.linspace(sol.t[0],sol.t[-1],1001)
        payload[name+'_times']=times;payload[name+'_values']=sol.sol(times)
    np.savez_compressed(raw,**payload)
    pole_spread=float(np.max(abs(np.array(a['pole_fields'])-np.array(b['pole_fields']))))
    return dict(schema='nsc-atlas-shape-transition-v1',
        scientific_owner='docs/nsc-atlas-valve-formation.md',
        charged_parent=asdict(model),volume_sector_boundaries=volume_thresholds(),
        charge_lattice_classification=half_integer_classification(),
        stationary_volume_sectors=[dict(q=q,extrema=all_volume_extrema(q)) for q in (0.,.5,1.,1.5,2.,2.5)],
        wrapped_material_source=wrapped_brane_couplings(0.),
        conditional_charged_handoff=dict(reference_s=0.,**handoff,
            scope='Reused ScalarFlux rapid-transfer identity, NOT the source-free saddle below; no rolling-parent nucleation radius inferred.'),
        symmetric_shape_branches=branches,
        source_free_transition=dict(q=1.5,membrane_inserted=False,charge_jump=0.,
            massive_scalar='phi=phi_star(3/2) exactly; chi=0, no massive-scalar matter energy is supplied.',
            coarse=a,fine=b,
            refinement=dict(pole_field_spread=pole_spread,proper_length_spread=abs(a['length']-b['length']),
                Euclidean_action_difference_spread=abs(a['action_difference_times_g']-b['action_difference_times_g']),
                scope='Observed changes under one collocation/pole-cutoff refinement, not rigorous error bounds.'),
            connecting_chart=connecting,interior=inside,exterior=outside,
            interior_table=snapshot_table(interior),exterior_table=snapshot_table(exterior),
            scope='One regular nonconstant numerical Euclidean saddle and its Lorentzian charts in the1+3 invariant sector. Not an exhaustive saddle search or a full physical decay-index calculation.'),
        verdict=dict(combined_scalar_charged_sector_matches_reference=True,
            nonconstant_profile_numerically_admitted=True,derived_exit_state=True,
            initial_open_region_regular=True,finite_interior_expansion=True,
            full_quantum_transition_admitted=False,black_hole_valve_trigger_proved=False,
            incoming_parent_material_transfer_measured=False,lasting_observed_cosmology=False,
            minimum_interior_Omega_k=inside['minimum_curvature_fraction'],
            meets_requested_nearflat_002_bound=inside['minimum_curvature_fraction']<=.002,
            SI_scale_selected=False,full_goal_achieved=False),
        provenance=dict(source_sha256=source,source_unchanged=True,producer_commit=commit,
            source_bound_to_commit=bound,raw_cache=raw.relative_to(ROOT).as_posix(),raw_cache_sha256=sha(raw)),
        compute=dict(threads=4,elapsed_seconds=time.perf_counter()-started))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-directory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();folder=(ROOT/args.raw_directory).resolve();out=(ROOT/args.output).resolve()
    if not folder.is_relative_to(ROOT/'.build') or folder.exists() or not out.is_relative_to(ROOT) or out.exists():
        raise ValueError('fresh raw directory under .build and fresh in-repository summary required')
    folder.mkdir(parents=True);result=calculation(folder);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,allow_nan=False,default=plain)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(out),verdict=result['verdict'],compute=result['compute']),indent=2))


if __name__=='__main__':main()
