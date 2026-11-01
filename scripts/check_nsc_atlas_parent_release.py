"""Bounded passive concentration/release test of the unchanged original parent.

All events come from the producer-bound dense solution.  No initial data,
physical perturbation or GR evolution is made by this readout.
"""
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import math
import multiprocessing
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from numba import set_num_threads
from scipy.integrate import cumulative_simpson, simpson
from scipy.optimize import brentq
from scipy.fft import set_workers
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))

from scripts.check_nsc_atlas_parent_valve import load_original_parent
from recursive_horizons.nsc_atlas_null_geometry import first_connection
from recursive_horizons.nsc_atlas_null_parent import dense_history_jet
from recursive_horizons.nsc_atlas_parent_gh import PAIRS
from recursive_horizons.nsc_atlas_parent_scale_selector import VolumeParent
from recursive_horizons.nsc_atlas_s3_frame import arbitrary_points_values
from recursive_horizons.nsc_atlas_trapped_parent import _first_geometry

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def scalar_point_jet(system, history, t):
    """Actual metric/scalar polynomial jets at the tracked geodesic event.

    Ordered spatial derivatives retain the regular-frame commutator.
    Neither the KG nor Einstein equation repairs the measured jets.
    """
    y, dy, ddy = dense_history_jet(history, float(t))
    C, C1, C2 = map(system.unpack, (y, dy, ddy))
    point = y[system.dof:system.dof+4].copy()
    point /= np.linalg.norm(point)
    Y = arbitrary_points_values(point[None, :], system.basis.labels)[0]
    value, rate = C@Y, C1@Y
    first_metric = np.array([system.r*(D@C[:10].T).T@Y for D in system.basis.generators])
    first_q = np.array([system.r*(D@C[20:22].T).T@Y for D in system.basis.generators])
    mixed_q = np.array([system.r*(D@C1[20:22].T).T@Y for D in system.basis.generators])
    g = np.empty((4, 4)); v = np.empty_like(g); dg = np.empty((3, 4, 4))
    for k, (a, b) in enumerate(PAIRS):
        g[a, b] = g[b, a] = value[k]
        v[a, b] = v[b, a] = rate[k]
        dg[:, a, b] = dg[:, b, a] = first_metric[:, k]
    connection = first_connection(g, v, dg, frame_scale=system.r)
    inverse, Gamma = connection['inverse'], connection['Gamma']
    N = 1/math.sqrt(-inverse[0, 0]); hi = np.linalg.inv(g[1:, 1:])
    beta = hi@g[0, 1:]; normal = np.r_[1., -beta]/N
    u_raw = y[system.dof+4:system.dof+8]
    unorm = float(u_raw@g@u_raw)
    if unorm >= 0:
        raise ValueError('Original tracked observer is not timelike')
    u = u_raw/math.sqrt(-unorm)
    derivative = np.vstack([rate[20:22], first_q])
    ordered = np.empty((4, 4, 2))
    ordered[0, 0] = C2[20:22]@Y
    ordered[0, 1:] = ordered[1:, 0] = mixed_q
    for i, Di in enumerate(system.basis.generators):
        for j, Dj in enumerate(system.basis.generators):
            ordered[i+1, j+1] = system.r**2*(Di@(Dj@C[20:22].T)).T@Y
    Hess = ordered-np.einsum('cab,cA->abA', Gamma, derivative)
    scalar_box = np.einsum('ab,abA->A', inverse, Hess)
    u_q = u@derivative; n_q = normal@derivative
    Hess_uu = np.einsum('a,b,abA->A', u, u, Hess)
    Hess_nn = np.einsum('a,b,abA->A', normal, normal, Hess)
    geometry = _first_geometry(g[None], v[None], dg[None], system.r)
    eigen = geometry['K_eigenvalues'][0]; traceK = float(sum(eigen))
    lap = np.einsum('ij,ijA->A', hi, Hess[1:, 1:])-traceK*n_q
    dN = -.5*N**3*np.einsum('a,iab,b->i', inverse[0], dg, inverse[:, 0])
    accel_grad = np.einsum('ij,i,jA->A', hi, dN/N, first_q)
    s = system.sstar+value[20]/math.sqrt(system.K); chi = value[21]
    mu2 = 400*math.exp(-12*s)
    potential_s = -84*math.exp(-12*s)+105*math.exp(-14*s)-18.75*math.exp(-24*s)
    gradW = np.array([(potential_s-6*mu2*chi*chi)/math.sqrt(system.K), mu2*chi])
    gradient_norm = np.einsum('ij,iA,jA->A', hi, first_q, first_q)
    density = .5*float(n_q@n_q+sum(gradient_norm))+float(system.potential(value[None, 20:22])[0])
    frame_u_rate = dy[system.dof+4:system.dof+8]*u_raw[0]
    acceleration_raw = frame_u_rate+np.einsum('abc,b,c->a', Gamma, u_raw, u_raw)
    row = dict(coordinate_t=float(t), tau=float(y[system.dof+8]),
        S=float(value[20]), chi=float(chi), s=float(s),
        uS=float(u_q[0]), uChi=float(u_q[1]), nS=float(n_q[0]), nChi=float(n_q[1]),
        HessS_uu=float(Hess_uu[0]), HessChi_uu=float(Hess_uu[1]),
        HessS_nn=float(Hess_nn[0]), HessChi_nn=float(Hess_nn[1]),
        spatial_lapS=float(lap[0]), spatial_lapChi=float(lap[1]),
        lapse_accel_dot_DS=float(accel_grad[0]), lapse_accel_dot_DChi=float(accel_grad[1]),
        traceK=traceK, gradient_chi_squared=float(gradient_norm[1]),
        normal_total_density=density, external_principal_expansions=(-eigen).tolist(),
        observer_unit_norm_error=abs(unorm+1), observer_geodesic_residual=float(max(abs(acceleration_raw))),
        scalar_KG_residual=(scalar_box-gradW).tolist(),
        scalar_hessian_symmetry_error=float(np.max(abs(Hess-Hess.swapaxes(0,1)))),
        metric_velocity_consistency=float(max(abs(rate[:10]-value[10:20]))),
        field_velocity_consistency=float(max(abs(rate[20:22]-value[22:24]))))
    return row


def numeric_columns(rows):
    """Helper's numerical arrays; leave vector diagnostics in the raw record."""
    return {key:np.array([row[key] for row in rows]) for key in rows[0]
            if isinstance(rows[0][key], (float,int))}


def worker(arguments):
    degree, folder, deadline = arguments
    from scripts.parent_release_region import RegionSampler, source_region, nodal_radius_control, analytic_geometry_controls
    from scripts.parent_release_response import causal_release_response
    began = time.perf_counter()
    set_num_threads(4)
    with threadpool_limits(limits=2), set_workers(2):
        data = load_original_parent(degree); system, history, end = data['system'], data['history'], data['end']
        region_definition = source_region(system.parent)
        region_definition['independent_nodal_control'] = nodal_radius_control(system.parent)
        region_definition['analytic_geometry_controls'] = analytic_geometry_controls(
            **region_definition['sampler_arguments'],frame_scale=system.r,order=12)
        region = RegionSampler(system, **region_definition['sampler_arguments'], order=12)
        # 129 uniformly matched GH leaves resolve the actual carrier cycles.
        # Anchor proper time is returned at every leaf, never used as globaltime.
        times = np.linspace(0., end, 129)
        point_rows, region_rows = [], []
        for t in times:
            if time.time() > deadline:
                raise RuntimeError('300-second passive readout budget exhausted')
            state = history.evaluate(float(t))
            point_rows.append(scalar_point_jet(system, history, t))
            region_rows.append(region.sample(state))
            if len(point_rows)%32==0:
                (Path(folder)/f'degree-{degree}-progress.json').write_text(
                    json.dumps(dict(status='partial_readout',degree=degree,point_series=point_rows,
                                    region_series=region_rows),allow_nan=False,separators=(',',':'))+'\n',
                    encoding='utf-8',newline='\n')
        arrays = numeric_columns(point_rows)
        response = causal_release_response(arrays, sstar=system.sstar)
        coarse_response = causal_release_response({k:v[::2] for k,v in arrays.items()}, sstar=system.sstar)
        high = RegionSampler(system, **region_definition['sampler_arguments'], order=16)
        region_checks = []
        for i in (0,32,64,96,128):
            if time.time() > deadline:
                raise RuntimeError('300-second passive readout budget exhausted')
            better = high.sample(history.evaluate(float(times[i])))
            region_checks.append(dict(index=i,lower=region_rows[i],higher=better))
    raw = dict(degree=degree,provenance=data['provenance'],interval=dict(coordinate_t=[0.,end],anchor_tau=[0.,.9]),
        region_definition=region_definition,point_series=point_rows,region_series=region_rows,
        region_quadrature_checks=region_checks,response=response,coarse_time_response=coarse_response)
    path = Path(folder)/f'degree-{degree}.json'
    path.write_text(json.dumps(raw,allow_nan=False,separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    return dict(degree=degree,raw_path=path.relative_to(ROOT).as_posix(),raw_sha256=sha(path),
                elapsed_seconds=time.perf_counter()-began)


def compact_results(runs):
    """Preserve raw arrays locally, export scalar conclusions and error indicators."""
    summaries=[]; originals=[]
    for run in runs:
        path=ROOT/run['raw_path']
        if sha(path)!=run['raw_sha256']:
            raise ValueError('Raw release readout identity changed')
        data=json.loads(path.read_bytes()); originals.append(data)
        p=data['point_series']; r=data['region_series']; t=np.array([x['coordinate_t'] for x in p])
        col=lambda key:np.array([x[key] for x in r])
        E,V=col('energy'),col('volume'); excitation=col('excitation_energy')
        work,flux=col('work'),col('outward_GH_flux')
        W=cumulative_simpson(work,x=t,initial=0.); F=cumulative_simpson(flux,x=t,initial=0.)
        exW=cumulative_simpson(col('excitation_work'),x=t,initial=0.)
        exF=cumulative_simpson(col('outward_excitation_GH_flux'),x=t,initial=0.)
        residue=E-E[0]-W+F; ex_residue=excitation-excitation[0]-exW+exF
        timestep=dict(full_net_work_minus_flux=float(W[-1]-F[-1]),
            half_sample_net_work_minus_flux=float(simpson((work-flux)[::2],x=t[::2])),
            indicator=float(abs((W[-1]-F[-1])-simpson((work-flux)[::2],x=t[::2]))))
        checks={}
        for key in ('energy','excitation_energy','volume','work','outward_GH_flux','external_volume',
                    'tube_outward_area_flux','tube_outward_GH_weighted_flux'):
            checks[key]=max(abs(row['higher'][key]-row['lower'][key]) for row in data['region_quadrature_checks'])
        s=np.array([x['s'] for x in p]); uS=np.array([x['uS'] for x in p])
        point_checks=dict(maximum_abs_scalar_KG_residual=float(np.max(abs(np.array([x['scalar_KG_residual'] for x in p])))),
            maximum_observer_unit_error=max(x['observer_unit_norm_error'] for x in p),
            maximum_observer_geodesic_residual=max(x['observer_geodesic_residual'] for x in p),
            maximum_scalar_Hessian_symmetry_error=max(x['scalar_hessian_symmetry_error'] for x in p),
            maximum_metric_velocity_consistency=max(x['metric_velocity_consistency'] for x in p),
            maximum_field_velocity_consistency=max(x['field_velocity_consistency'] for x in p))
        summary=dict(**run,provenance=data['provenance'],interval=data['interval'],region_definition=data['region_definition'],
            local=dict(initial_s=float(s[0]),maximum_s=float(max(s)),final_s=float(s[-1]),
                maximum_s_anchor_tau=p[int(np.argmax(s))]['tau'],initial_uS=float(uS[0]),final_uS=float(uS[-1]),
                internal_volume_ratio_at_end=float(math.exp(12*(s[-1]-s[0]))),
                initial_point_density=p[0]['normal_total_density'],final_point_density=p[-1]['normal_total_density'],
                initial_external_principal_expansions=p[0]['external_principal_expansions'],
                final_external_principal_expansions=p[-1]['external_principal_expansions']),
            region=dict(initial=r[0],final=r[-1],
                energy_ratio=float(E[-1]/E[0]),excitation_energy_ratio=float(excitation[-1]/excitation[0]),
                volume_ratio=float(V[-1]/V[0]),mean_density_ratio=float((E[-1]/V[-1])/(E[0]/V[0])),
                maximum_mean_density_ratio=float(max(E/V)/(E[0]/V[0])),
                cumulative_geometric_work=float(W[-1]),cumulative_outward_GH_flux=float(F[-1]),
                cumulative_outward_material_flux=float(simpson(col('outward_material_flux'),x=t)),
                cumulative_shift_flux=float(simpson(col('shift_flux'),x=t)),
                cumulative_tube_outward_GH_weighted_flux=float(simpson(col('tube_outward_GH_weighted_flux'),x=t)),
                maximum_tube_flux_boost_identity_error=float(max(col('tube_flux_boost_identity_max'))),
                tube_flux_scope='Sum of local tube-observer fluxes with their GH lapse N/gamma, not a conserved total/anchor-clock energy or an attached-child supply.',
                maximum_abs_energy_ledger_residual=float(max(abs(residue))),
                final_energy_ledger_residual=float(residue[-1]),
                maximum_abs_excitation_ledger_residual=float(max(abs(ex_residue))),
                final_excitation_ledger_residual=float(ex_residue[-1]),
                temporal_quadrature=timestep,spatial_quadrature=checks),
            response=data['response']['summary'],response_scope=data['response']['scope'],
            time_refinement=dict(coarse_response=data['coarse_time_response']['summary']),
            point_checks=point_checks)
        summaries.append(summary)
    a,b=originals; pa,pb=a['point_series'],b['point_series']; ra,rb=a['region_series'],b['region_series']
    comparison=dict(degrees=[a['degree'],b['degree']],
        maximum_anchor_tau_misalignment=max(abs(x['tau']-y['tau']) for x,y in zip(pa,pb)),
        point_maximum_absolute_differences={key:max(abs(x[key]-y[key]) for x,y in zip(pa,pb))
            for key in ('s','uS','normal_total_density')},
        region_maximum_absolute_differences={key:max(abs(x[key]-y[key]) for x,y in zip(ra,rb))
            for key in ('energy','volume','work','outward_GH_flux','external_volume')},
        scope='Empirical paired-GH-leaf resolution differences; the anchor-clock mismatch is reported, not treated as a global clock.')
    return summaries,comparison


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-directory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--budget-seconds',type=float,default=300.)
    args=parser.parse_args(); raw=(ROOT/args.raw_directory).resolve(); out=(ROOT/args.output).resolve()
    if not math.isfinite(args.budget_seconds) or args.budget_seconds<=0:
        raise ValueError('Positive finite compute allowance required')
    if not raw.is_relative_to(ROOT/'.build') or raw.exists() or not out.is_relative_to(ROOT) or out.exists():
        raise ValueError('Fresh raw directory under.build and fresh repository output required')
    names=['scripts/check_nsc_atlas_parent_release.py','scripts/parent_release_region.py','scripts/parent_release_response.py',
           'scripts/check_nsc_atlas_parent_valve.py','src/recursive_horizons/nsc_atlas_null_geometry.py',
           'src/recursive_horizons/nsc_atlas_null_parent.py','src/recursive_horizons/nsc_atlas_gh_matter_work.py']
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    pins={name:sha(ROOT/name) for name in names}
    for name,digest in pins.items():
        if hashlib.sha256(subprocess.check_output(['git','show',commit+':'+name],cwd=ROOT)).hexdigest()!=digest:
            raise ValueError('Commit the calculation source before running')
    raw.mkdir(parents=True); began=time.perf_counter(); deadline=time.time()+args.budget_seconds
    with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as pool:
        runs=list(pool.map(worker,[(degree,str(raw),deadline) for degree in (20,24)]))
    result=dict(schema='nsc-atlas-parent-release-v1',scientific_owner='docs/nsc-atlas-valve-formation.md',
        runs=runs,readout_commit=commit,readout_source_sha256=pins,
        scope=dict(no_new_physical_states=True,no_GR_evolution=True,physical_child_claim=False),
        elapsed_seconds=time.perf_counter()-began)
    # Complete raw data survives even if a later compact-summary exporterfails.
    (raw/'index.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    result['runs'],result['resolution_comparison']=compact_results(runs)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(out),elapsed_seconds=result['elapsed_seconds'],
        readout_commit=commit,scope=result['scope']),ensure_ascii=True,indent=2))


if __name__=='__main__':main()
