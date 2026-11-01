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
    rolling_parent_limit,rolling_parent_state,wall_energy_jump,
)
from recursive_horizons.nsc_atlas_shape_transition import (
    solve_bounce,bounce_diagnostics,continue_pole,connecting_geometry,
    evolve_rolling_parent,
    critical_parent_response,
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


def rolling_calculation(folder):
    """The parent at the break: full-potential histories and independent uplift."""
    import math
    from scipy.integrate import solve_ivp,quad

    started=time.perf_counter();source=pins();raw_payload={};rows=[]
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bound=all(hashlib.sha256(subprocess.check_output(['git','show',commit+':'+name],cwd=ROOT)).hexdigest()==digest
        for name,digest in source.items())
    limit=rolling_parent_limit();A=float(limit['A']);T=25/32
    # Independent sixteen-dimensional directional Einstein equations in
    # bulk proper time. These do not call the reduced radion RHS.
    def bulk_rhs(t,y):
        if time.perf_counter()-started>60:raise TimeoutError('60s combined parent allowance exhausted')
        s,hx,hi,N=y;theta=3*hx+12*hi;flux=T*math.exp(-12*s)
        return [hi,(A-5*flux)/7-theta*hx,
            (A+2*flux)/7-theta*hi-1.25*math.exp(-2*s),hx+6*hi]

    with threadpool_limits(limits=1):
        for index,(s0,v0) in enumerate(((0.,0.),(1.,.1))):
            sol=evolve_rolling_parent(initial_s=s0,initial_v=v0)
            initial=rolling_parent_state(s0,v0);end=rolling_parent_state(*sol.y[:2,-1])
            def target(t,y):return y[3]-40
            target.terminal=True;target.direction=1
            bulk=solve_ivp(bulk_rhs,(0.,100.),[s0,initial['bulk_external_H'],initial['bulk_internal_H'],0.],
                method='DOP853',rtol=2e-11,atol=2e-13,max_step=.1,events=target,dense_output=True)
            if not bulk.success or not len(bulk.t_events[0]):raise RuntimeError('independent bulk history did not reach the declared40efolds')
            sample=np.linspace(0,40,401);Y=sol.sol(sample)
            states=[rolling_parent_state(*y) for y in Y[:2].T]
            direct=bulk.sol(Y[2]);hx=np.array([r['bulk_external_H'] for r in states]);hi=np.array([r['bulk_internal_H'] for r in states])
            direct_H_difference=float(max(np.max(abs(direct[1]-hx)),np.max(abs(direct[2]-hi))))
            trace_difference=max(abs(r['bulk_Ricci_metric']-r['bulk_Ricci_source']) for r in states)
            theta=3*bulk.y[1]+12*bulk.y[2]
            constraint=(theta**2-3*bulk.y[1]**2-12*bulk.y[2]**2+15*np.exp(-2*bulk.y[0]))/2-A-T*np.exp(-12*bulk.y[0])
            def density_work(N):
                s,v,_=sol.sol(N);state=rolling_parent_state(s,v)
                return 84*v*v*math.exp(state['log_Einstein_density'])
            work=quad(density_work,0,40,epsabs=2e-10,epsrel=2e-10)[0]
            balance=math.exp(end['log_Einstein_density'])-math.exp(initial['log_Einstein_density'])+work
            if direct_H_difference>1e-8 or trace_difference>1e-9 or max(abs(constraint))>1e-8 or abs(balance)>1e-8:
                raise RuntimeError('rolling-parent independent Einstein/trace/work control failed')
            raw_payload[f'parent{index}_N']=sample;raw_payload[f'parent{index}_states']=Y
            raw_payload[f'parent{index}_bulk_times']=bulk.t;raw_payload[f'parent{index}_bulk_states']=bulk.y
            rows.append(dict(initial_s=s0,initial_v=v0,efolds=40.,endpoint=end,
                bulk_proper_duration=float(Y[2,-1]),
                independent_bulk_H_difference=direct_H_difference,
                independent_bulk_Hamiltonian_max=float(max(abs(constraint))),
                bulk_Ricci_metric_source_max_difference=float(trace_difference),
                Einstein_density_work_balance=float(balance),
                sample_table=dict(columns=['N','s','v','bulk_proper_duration'],rows=np.column_stack([sample,Y.T])[::40].tolist()),
                scope='Declared verification initial data, not a selected ancestral state or physical transition epoch.'))
    if source!=pins():raise RuntimeError('rolling-parent sources changed during calculation')
    raw=folder/'rolling-parent.npz';np.savez_compressed(raw,**raw_payload)
    return dict(schema='nsc-atlas-rolling-parent-v1',scientific_owner='docs/nsc-atlas-valve-formation.md',
        action=asdict(combined_record()),asymptotic_limit=limit,
        parent_stationary_radius_families=symmetric_radius_branches(2.),histories=rows,
        wall_conservation=dict(normal='child to parent',
            identity='[T(n,U)]=dot(sigma)=sigma_s*U(s)=-6sigma*U(s)',
            surface_first_law='dot(E_surface)-sigma*dot(area)=area*(F_child-F_parent)',
            F_definition='F=-T(n,U); only its jump is fixed, not either side flux.',
            zero_tangent_control=wall_energy_jump(0.,0.),
            side_fluxes_selected=False,source_inserted=False,
            cap_normal_force='[V]=sigma*K_mean+sigma_s*mean(n(s))',
            cap_necessary_balance='Sum integral(R^3 V_s dxi)=6sigma Rwall^3; q2 common-radius integrand is negative.'),
        conditional_scalar_handoff=dict(**combined_record().scalar_sector.prepare(parent_flux=2.,charge=.5),
            scope='Existing instantaneous continuous-field identity only; it does not provide a membrane trajectory or a formed child.'),
        verdict=dict(nonempty_regular_parent_asymptotic=True,
            parent_limit_fixed_in_supplied_action=True,finite_parent_history_selected=False,
            individual_receiver_material_flux_derived=False,formed_lasting_child=False,
            quantum_state_selected=False,SI_scale_selected=False,full_goal_achieved=False),
        provenance=dict(source_sha256=source,source_unchanged=True,producer_commit=commit,
            source_bound_to_commit=bound,raw_cache=raw.relative_to(ROOT).as_posix(),raw_cache_sha256=sha(raw)),
        compute=dict(threads=1,elapsed_seconds=time.perf_counter()-started,
            allowance='Each reduced history20s; complete calculation60s; two reduced and two independent bulk histories, bounded40Einstein e-folds.'))


def critical_calculation(folder):
    """One reference-history closure and its derived constrained Cauchy data."""
    from recursive_horizons.nsc_atlas_closed_pocket import solve_critical_pocket_data
    from recursive_horizons.nsc_atlas_parent_circulation import Parent
    from threadpoolctl import threadpool_limits
    Parent.retained()  # Hydrate the same immutable normalization dependencies.
    started=time.perf_counter();source=pins()
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    with threadpool_limits(limits=1):
        response=critical_parent_response(carrier_scale=.25,budget_seconds=30.)
        data=solve_critical_pocket_data(40,carrier_scale=.25,compression=.01,
            critical_record=response,budget_seconds=10.)
    after=pins()
    if any(after.get(name)!=digest for name,digest in source.items()):
        raise RuntimeError('critical preparation sources changed during calculation')
    bound=all(hashlib.sha256(subprocess.check_output(['git','show',commit+':'+name],cwd=ROOT)).hexdigest()==digest
        for name,digest in after.items())
    raw=folder/'reference-coefficients.npz'
    np.savez_compressed(raw,canonical_S_U=np.asarray(response['canonical_S_critical_U']),
        canonical_PiS_U=np.asarray(response['canonical_PiS_critical_U']),
        PiChi_U=np.asarray(response['updated_PiChi_conformal_U']),
        conformal_coefficients=data.coefficients,conformal_labels=data.labels,
        S_coefficients=data.S_coefficients,S_labels=data.S_labels,
        PiS_coefficients=data.momentum_S_coefficients,PiS_labels=data.momentum_S_labels,
        York_coefficients=data.York_primitive_coefficients,York_labels=data.York_primitive_labels)
    return dict(schema='nsc-atlas-critical-parent-v1',scientific_owner='docs/nsc-atlas-valve-formation.md',
        response=response,joint_Cauchy=dict(positive_psi_lower_bound=data.positivity_lower_bound,
            continuation=data.continuation,compression=.01,balanced_mode='PiS contribution=-sqrt84*.01*U2/3',
            carrier_scale=.25,carrier_correction_reapplied=False,
            scope='Full nonlinear initial constraints on the conditioned reference coefficients; no GR history is calculated here.'),
        verdict=dict(reference_history_CK_closure_admitted=bool(response['admitted']),
            positive_joint_Cauchy_data_constructed=True,phase_fixed_within_reference_closure=True,
            nonlinear_criticality_proved=False,natural_preparation_selected=False,
            formed_child_proved=False,SI_scale_selected=False,full_goal_achieved=False),
        provenance=dict(producer_commit=commit,source_sha256=after,source_unchanged=True,
            source_bound_to_commit=bound,raw_cache=raw.relative_to(ROOT).as_posix(),raw_cache_sha256=sha(raw)),
        compute=dict(threads=1,elapsed_seconds=time.perf_counter()-started,
            allowance='30s reference response plus10s joint constraint solve; no nonlinear GR evolution.'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-directory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--rolling-parent',action='store_true',help='Evaluate the derived q2 parent limit and declared histories, not the frozen shape bounce')
    mode.add_argument('--critical-parent',action='store_true',help='Recompute bounded carrier history, weak CK phase and joint initial constraints only')
    args=parser.parse_args();folder=(ROOT/args.raw_directory).resolve();out=(ROOT/args.output).resolve()
    if not folder.is_relative_to(ROOT/'.build') or folder.exists() or not out.is_relative_to(ROOT) or out.exists():
        raise ValueError('fresh raw directory under .build and fresh in-repository summary required')
    folder.mkdir(parents=True)
    result=critical_calculation(folder) if args.critical_parent else (rolling_calculation(folder) if args.rolling_parent else calculation(folder))
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,allow_nan=False,default=plain)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(out),verdict=result['verdict'],compute=result['compute']),indent=2))


if __name__=='__main__':main()
