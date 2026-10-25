"""Finite-future formation budget; no new nonlinear grids."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import sympy as sp
from scipy.optimize import brentq

from recursive_horizons.nsc_atlas_formation_window import FormationBackground,transfer,dust_by_quadrature,density_phase,elementary_gain_bound,compact_packet

ROOT=Path(__file__).resolve().parents[1]
INPUTS=('src/recursive_horizons/nsc_atlas_formation_window.py','src/recursive_horizons/nsc_atlas_scale_hierarchy.py',
    'src/recursive_horizons/nsc_atlas_scalar_core.py','src/recursive_horizons/nsc_atlas_shared_scalar_slice.py',
    'src/recursive_horizons/nsc_atlas_flux_conversion.py','src/recursive_horizons/nsc_atlas_flux_transition.py',
    'scripts/check_nsc_atlas_formation_window.py','scripts/plot_nsc_atlas_formation_window.py',
    'tests/test_nsc_atlas_formation_window.py')


def symbolic_checks():
    u,L,D,C,nu=sp.symbols('u L D C nu',positive=True)
    F=L/3+D*u**3/3-C*u**2;H=sp.sqrt(F)
    # Dust mode H in reciprocal scale: H_uu+(H_u/H-1/u)H_u-D*u*H/(2F)=0.
    dust=sp.simplify(sp.diff(H,u,2)+(sp.diff(H,u)/H-1/u)*sp.diff(H,u)-D*u*H/(2*F))
    p,delta=sp.symbols('p delta');du=-u*p/H;pu=(nu**2*u-D/2)*delta/H
    second=sp.diff(-u/H,u)*p-u/H*pu
    residual=sp.simplify(second+(sp.diff(H,u)/H-1/u)*du+(nu**2*u**2-D*u/2)*delta/F)
    return dict(dust_H_mode=str(dust),reciprocal_scale_density_equation=str(residual),
                minimum_F=str(sp.simplify(F.subs(u,2*C/D))))


def worker(job):
    name,k,A,tol=job;cpu=time.process_time()
    row,raw=transfer(k,initial_scale=A,rtol=tol)
    dust=dust_by_quadrature(A)
    row['independent_dust_future_response']=dust.tolist()
    row['dust_density_comparison_margin']=float(dust[0]-row['sampled_max_density_response'])
    if k>0:row['independent_phase_difference']=density_phase(k,initial_scale=A,rtol=tol)-row['future_density_response'][0]
    row['cpu_seconds']=time.process_time()-cpu
    return name,row


def packet_worker(job):
    name,L,A,order,Q,tol=job;cpu=time.process_time()
    row,raw=compact_packet(L,initial_scale=A,order=order,qmax=Q,rtol=tol)
    seed=.05;sol=raw['solution'];weights=raw['weights']
    def value(u):return seed*float(weights@sol.sol(u)[:order])
    if value(0)>1:
        u1=brentq(lambda u:value(u)-1,0,1/A,xtol=1e-11)
        u01=brentq(lambda u:value(u)-.1,u1,1/A,xtol=1e-11)
        row['scale_at_formal_linear_delta_1']=1/u1
        row['scale_at_linear_delta_0p1']=1/u01
        error=seed*row['all_future_pointwise_tail_error_bound']
        row['formal_unit_density_tail_bracket']=sorted(1/brentq(lambda u:value(u)+sign*error-1,0,1/A,xtol=1e-11) for sign in (-1,1))
    row['supplied_peak_density_contrast']=seed
    row['curve']=dict(reciprocal_scale=raw['u'].tolist(),central_linear_contrast=(seed*raw['central']).tolist())
    row['cpu_seconds']=time.process_time()-cpu
    return name,row


def calculation():
    start=time.process_time();wall=time.perf_counter();b=FormationBackground.inherited()
    prior=json.loads((ROOT/'results/summaries/nsc-atlas-compensation-response-v1.json').read_text(encoding='utf-8'))
    geometry=json.loads((ROOT/'results/summaries/nsc-atlas-scale-hierarchy-v1.json').read_text(encoding='utf-8'))
    p=prior['runs']['fine'];khalf=p['half_power_wavenumber_atlas'];kd=1/p['reservoir_radius']
    kcore=1/geometry['runs']['small_compactness']['core']['half_mass_radius']
    kJ=(2*b.excitation*b.frequency*b.frequency)**.25
    jobs=[('dust_comparison',0.,1.,2e-11),('supply_radius',kd,1.,2e-11),
          ('compensation_turnover',khalf,1.,2e-11),('Jeans_scale',kJ,1.,2e-11),
          ('core_half_radius',kcore,1.,2e-11),('dust_tight',0.,1.,2e-13),
          ('turnover_tight',khalf,1.,2e-13),('core_tight',kcore,1.,2e-13),
          ('earlier_A_0p1',khalf,.1,2e-11),('earlier_A_0p01',khalf,.01,2e-11),
          ('earlier_A_0p001',khalf,.001,2e-11),('earlier_0p01_tight',khalf,.01,2e-13)]
    with ProcessPoolExecutor(max_workers=4) as pool:rows=dict(pool.map(worker,jobs))
    L=p['reservoir_radius']
    packet_jobs=[('coarse',L,.01,128,40.,2e-10),('quadrature',L,.01,256,40.,2e-10),
                 ('tail_control',L,.01,256,60.,2e-10),('time_control',L,.01,256,60.,2e-12),
                 ('late_parent',L,1.,128,40.,2e-10)]
    with ProcessPoolExecutor(max_workers=4) as pool:packets=dict(pool.map(packet_worker,packet_jobs))
    packet_errors={}
    for key in ('central_future_gain','scale_at_formal_linear_delta_1','scale_at_linear_delta_0p1'):
        packet_errors[key]=2*max(abs(packets['quadrature'][key]-packets['coarse'][key]),
            abs(packets['tail_control'][key]-packets['quadrature'][key]),
            abs(packets['time_control'][key]-packets['tail_control'][key]))
    ranges={}
    for label,tight in (('dust_comparison','dust_tight'),('compensation_turnover','turnover_tight'),
                        ('core_half_radius','core_tight'),('earlier_A_0p01','earlier_0p01_tight')):
        ranges[label]=[2*abs(x-y) for x,y in zip(rows[label]['future_density_response'],rows[tight]['future_density_response'])]
    budget=elementary_gain_bound()
    q=rows['dust_tight']['future_density_response']
    necessary=dict(zero_current_initial_amplitude_for_linear_unit_density=1/q[0],
        initial_velocity_parameter_for_linear_unit_density_at_zero_initial_density=1/q[1],
        interpretation='Necessary linear extrapolation thresholds, not sufficient nonlinear collapse conditions. Initial velocity parameter is delta_dot/H.')
    # A declared early test state, not an observational seed or selected spectrum.
    row,raw=transfer(khalf,initial_scale=.01,rtol=2e-13,samples=801)
    seed=.02;sol=raw['solution'];ui=100.
    threshold_u=brentq(lambda u:seed*sol.sol(u)[0]-1.,0,ui,xtol=1e-11)
    controlled_u=brentq(lambda u:seed*sol.sol(u)[0]-.1,threshold_u,ui,xtol=1e-11)
    early=dict(initial_scale=.01,supplied_density_amplitude=seed,initial_current=0.,
        k_atlas=khalf,scale_at_linear_delta_0p1=1/controlled_u,
        scale_at_formal_linear_delta_1=1/threshold_u,
        future_linear_extrapolation=seed*row['future_density_response'][0],
        scope='Single local mode on the earlier averaged-parent continuation. Delta=1 is beyond linear validity and only identifies where a nonlinear continuation is needed. No shared-wall past or thermal history has been constructed.')
    curves={}
    for name in ('dust_comparison','compensation_turnover','Jeans_scale','core_half_radius','earlier_A_0p01'):
        row,raw=transfer(rows[name]['k'],initial_scale=rows[name]['initial_scale'],samples=161)
        # Progress 0..1 covers the complete future even though A diverges.
        curves[name]=dict(reciprocal_scale=raw['u'].tolist(),density=raw['response'][:,0,0].tolist())
    errors=[]
    if any(symbolic_checks()[k]!='0' for k in ('dust_H_mode','reciprocal_scale_density_equation')):errors.append('symbolic equation mismatch')
    for name,row in rows.items():
        if row['determinant_relative_error']>2e-7:errors.append(name+': transfer determinant drift')
        if row['dust_density_comparison_margin']< -2e-7:errors.append(name+': positive comparison bound violated')
        if abs(row.get('independent_phase_difference',0))>3e-7:errors.append(name+': phase equation disagrees')
    dust_error=float(np.max(abs(np.array(q)-dust_by_quadrature())))
    if dust_error>2e-10:errors.append('dust quadrature disagrees')
    for name,row in packets.items():
        if abs(row['initial_normalization']-1)>row['initial_omitted_spectral_weight_bound']+1e-10:
            errors.append(name+': compact spectral normalization failed')
    if packet_errors['central_future_gain']>.01:errors.append('compact seed amplification unresolved')
    total=time.process_time()-start+sum(r['cpu_seconds'] for r in rows.values())+sum(r['cpu_seconds'] for r in packets.values())
    return dict(schema='nsc-atlas-formation-window-v1',
        result='The inherited expanding parent supplies a finite all-wavelength linear growth budget; an earlier supplied mode identifies a nonlinear formation window.',
        fixed_background=dict(vacuum=b.vacuum,excitation=b.excitation,curvature=b.curvature,frequency=b.frequency),
        runs=rows,numerical_comparison_ranges=ranges,elementary_all_wavelength_bound=budget,
        necessary_linear_seed_thresholds=necessary,earlier_conditional_window=early,
        compact_seed=dict(definition='g(r)=(1-r^2/L^2)^4_+; delta_i=-Laplacian(g*g) normalized to peak .05; zero phase/current. L is the preceding initial reservoir radius in comoving atlas units. Support 2L, exact zero monopole and nonnegative Fourier density weight.',
            runs=packets,numerical_comparison_ranges=packet_errors,
            scope='An explicitly supplied small compact disturbance, not an inherited stochastic seed law or a formed core. Spectral tail bound applies throughout the linear envelope future. Delta=1 only marks failure of linear extrapolation.'),
        symbolic_identities=symbolic_checks(),independent_dust_quadrature_max_gap=dust_error,
        curves=curves,numerical_validation_errors=errors,
        compute_budget_seconds=120,process_cpu_seconds_including_workers=total,
        wall_seconds=time.perf_counter()-wall,budget_exhausted=total>120,
        limits='Linear local envelope modes on the declared averaged parent, no installed nonlinear core. Infinite future is a regular reciprocal-scale endpoint of these equations, not proof of global infinite-time full-GR validity. Earlier epochs are conditional homogeneous continuations, not the existing shared-wall past. No primordial spectrum, abundance, collapse or observational confirmation.',
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS},
        evidence_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ('results/summaries/nsc-atlas-compensation-response-v1.json','results/summaries/nsc-atlas-scale-hierarchy-v1.json')},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_formation_window.py --output results/summaries/nsc-atlas-formation-window-v1.json --figure')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--figure',action='store_true');args=p.parse_args()
    r=calculation()
    if args.output:args.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:r[k] for k in ('elementary_all_wavelength_bound','necessary_linear_seed_thresholds','earlier_conditional_window','process_cpu_seconds_including_workers','numerical_validation_errors')}))
    if r['budget_exhausted'] or r['numerical_validation_errors']:raise RuntimeError('Formation-window checks failed; inspect retained record before continuing')
    if args.figure:
        from plot_nsc_atlas_formation_window import render
        render(r)
