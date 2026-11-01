"""Finite directional initial data and its initial flux; no evolution.

Same canonical action,chi=A*u,u=X1. The declared conformal momentum
f=B*sqrt(3/5)*(u²-1/6) has the old reference L2 norm but a different
relative phase/direction. Old cap-reflection evolution masks are not used.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from scipy.special import roots_jacobi
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_valve_neck import solve_necks,zonal_rule,NeckInitialData

ROOT=Path(__file__).resolve().parents[1]
VOLUME=2*math.pi**2


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def imported_sources():
    package=(ROOT/'src/recursive_horizons').resolve();paths={Path(__file__).resolve()}
    for module in list(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if name:
            path=Path(name).resolve()
            if path.suffix=='.py' and path.is_relative_to(package):paths.add(path)
    return {path.relative_to(ROOT).as_posix():sha(path) for path in sorted(paths)}


def psi_lower(coefficients):
    labels=np.arange(0,2*len(coefficients),2)
    center=float(coefficients[0])/math.sqrt(VOLUME)
    tail=math.fsum((int(l)+1)*abs(float(c)) for l,c in zip(labels[1:],coefficients[1:]))/math.sqrt(VOLUME)
    pad=128*np.finfo(float).eps*(abs(center)+tail)
    return center-tail-pad


def fields(initial,c,lam,u,A,B,C):
    data=NeckInitialData(initial.parent,40,1.,np.asarray(c),())
    psi=data.values(u);du=data.values(u,1);duu=data.values(u,2)
    p=initial.parent;H=p.H;v=1-u*u
    f=B*math.sqrt(3/5)*(u*u-1/6);Ar=C*u*v/6;At=-Ar/2
    Kr=lam*Ar*psi**-6;Kt=-Kr/2
    Ktu=-lam*C/12*((1-3*u*u)*psi**-6-6*u*v*psi**-7*du)
    q=psi*psi/H;bu=2*psi*du/H
    Pi=lam*f*psi**-6;grad2=H*H*A*A*v*psi**-4
    W=p.vacuum+.5*p.mass_squared*A*A*u*u
    rho=.5*(Pi*Pi+grad2)+W;pr=.5*(Pi*Pi+grad2)-W;pt=.5*(Pi*Pi-grad2)-W
    lap=v*duu-3*u*du;R3=H*H*psi**-5*(-8*lap+6*psi)
    K2=Kr*Kr+2*Kt*Kt
    Hconstraint=R3-K2-2*p.kappa*rho
    # Independent physical radial check in increasing-chi orientation.
    jchi=Pi*A*np.sqrt(v)/q
    momentum=2*((u/np.sqrt(v)-np.sqrt(v)*bu/q)*(Kr-Kt)+np.sqrt(v)*Ktu)/q-p.kappa*jchi
    return dict(psi=psi,f=f,Pi=Pi,gradient2=grad2,W=W,rho=rho,pr=pr,pt=pt,Kr=Kr,Kt=Kt,K2=K2,
                Hamiltonian=Hconstraint,momentum=momentum,pressure_work_density=Kr*pr+2*Kt*pt)


def accounting(initial,c,lam,u,w,A,B,C):
    p=initial.parent;g=fields(initial,c,lam,u,A,B,C);psi=g['psi'];measure=w*psi**6/p.H**3
    normal=float(measure@(g['Pi']**2))/2
    gradient=float(measure@g['gradient2'])/2
    mass=float(measure@(.5*p.mass_squared*A*A*u*u));vacuum=p.vacuum*float(np.sum(measure))
    integratedK=float(measure@g['K2'])
    return dict(proper_volume=float(np.sum(measure)),normal_Pi_energy=normal,spatial_gradient_energy=gradient,
                scalar_mass_energy=mass,vacuum_energy=vacuum,total_proper_matter_energy=normal+gradient+mass+vacuum,
                integrated_KijKij=integratedK,shear_constraint_energy=integratedK/(2*p.kappa),
                normal_Pi_sample_max=float(max(abs(g['Pi']))),normal_Pi_volume_RMS=math.sqrt(2*normal/float(np.sum(measure))),
                shear_sample_max=math.sqrt(float(max(g['K2']))),
                full_pressure_work_rate=float(measure@g['pressure_work_density']))


def initial_cap_ledger(initial,c,lam,A,B,C,order=240):
    # sqrt(1-u²)du on[0,1], with the endpoint root handled in the Jacobi rule.
    x,wx=roots_jacobi(order,.5,0.);u=(x+1)/2;w=wx*np.sqrt(3+x)/4
    g=fields(initial,c,lam,u,A,B,C);measure=4*math.pi*w*g['psi']**6/initial.parent.H**3
    work=float(measure@g['pressure_work_density'])
    psi0=float(NeckInitialData(initial.parent,40,1.,np.asarray(c),()).values(0.))
    influx=4*math.pi*lam*A*B*math.sqrt(3/5)*psi0**-4/(6*initial.parent.H)
    energy=float(measure@g['rho'])
    return dict(positive_cap='u>=0; outward normal points toward increasingchi/decreasingu.',
                initial_lapse=1.,initial_shift=0.,fixed_interface_u=0.,boundary_sweep_rate=0.,
                plus_cap_proper_matter_energy=energy,plus_cap_material_influx=influx,
                minus_cap_material_influx=-influx,plus_cap_pressure_work_rate=work,minus_cap_pressure_work_rate=-work,
                plus_cap_energy_derivative_from_initial_ledger=influx+work,
                minus_cap_energy_derivative_from_initial_ledger=-influx-work,
                formula='Phi_in,+ =4pi*lambda*A*B*sqrt(3/5)*psi(0)^(-4)/(6H); work=integral_(u>0)Kr*(pr-pt)dV.',
                cap_quadrature_order=order,
                scope='Initial unit-lapse normal-time fixed-cap ledger; actual material influx, shear pressure work and sweep are distinct. The interface is not identified as a horizon/child.')


def calculation(directory):
    directory=directory.resolve()
    if not directory.is_relative_to(ROOT/'.build'):
        raise ValueError('raw initial arrays must remain under .build')
    if (directory/'initial-data.npz').exists() or (directory/'result.json').exists():
        raise ValueError('fresh initial-data output directory required')
    started=time.perf_counter();before=imported_sources();path=(0.,.25,.5,.75,1.)
    with threadpool_limits(limits=4):
        initial=solve_necks(40,(1.,))[0];p=initial.parent;H=p.H;z=p.mass_squared/H**2
        a=24/(z+3);A=math.sqrt(a/p.kappa);B=H**1.5*math.sqrt(12)/math.pi;C=p.kappa*A*B*math.sqrt(3/5)
        u,w,labels,Y=zonal_rule(40,205)
        f=B*math.sqrt(3/5)*(u*u-1/6);A2=C*C*u*u*(1-u*u)**2/24;Q=(A2+p.kappa*f*f)/H**2
        alpha=a*(1-u*u);beta=6+a*z*u*u;linear=8*labels*(labels+2)+6
        c=initial.coefficients.copy();stages=[]
        def residual(coeff,lam):
            psi=Y@coeff
            return linear*coeff-Y.T@(w*(alpha*psi+beta*psi**5+lam*lam*Q*psi**-7))
        for lam in path:
            for iteration in range(18):
                if time.perf_counter()-started>120:raise TimeoutError('one directional initial-data120s budget exceeded')
                psi=Y@c;R=residual(c,lam);norm=float(np.linalg.norm(R))
                if norm<2e-12:break
                multiplier=alpha+5*beta*psi**4-7*lam*lam*Q*psi**-8
                J=np.diag(linear)-Y.T@((w*multiplier)[:,None]*Y)
                step=np.linalg.solve(J,-R)
                for halve in range(18):
                    trial=c+step*2.**-halve
                    if psi_lower(trial)>0 and np.linalg.norm(residual(trial,lam))<norm*(1-1e-4*2.**-halve):
                        c=trial;break
                else:raise RuntimeError(f'directional lambda{lam}: strict-positive line search failed')
            else:raise RuntimeError(f'directional lambda{lam}: Newton budget exhausted')
            stages.append(dict(kick_amplitude=lam,Newton_steps=iteration,
                projected_residual_norm=float(np.linalg.norm(residual(c,lam))),finite_polynomial_psi_lower_bound=psi_lower(c),
                accounting=accounting(initial,c,lam,u,w,A,B,C)))
        uv,wv,_,_=zonal_rule(40,401);g=fields(initial,c,1.,uv,A,B,C);measure=wv*g['psi']**6/H**3
        volume=float(np.sum(measure))
        validation=dict(radial_order=401,solve_radial_order=205,
            physical_Hamiltonian_RMS_over_H2=math.sqrt(float(measure@(g['Hamiltonian']**2))/volume)/H**2,
            physical_Hamiltonian_sample_max_over_H2=float(max(abs(g['Hamiltonian'])))/H**2,
            physical_radial_momentum_RMS_over_H2=math.sqrt(float(measure@(g['momentum']**2))/volume)/H**2,
            physical_radial_momentum_sample_max_over_H2=float(max(abs(g['momentum'])))/H**2,
            reference_f_L2_squared=float(wv@(g['f']**2))/H**3,
            finite_polynomial_psi_lower_bound=psi_lower(c),sampled_psi_min=float(min(g['psi'])),sampled_psi_max=float(max(g['psi'])),
            accounting=accounting(initial,c,1.,uv,wv,A,B,C),
            scope='Independent higher fixed quadrature; curvature from actual psi derivatives, not source-reconstructed; radial momentum from physicalKr,Kt,psi derivatives andj.')
        cap=initial_cap_ledger(initial,c,1.,A,B,C,240)
        cap_check=initial_cap_ledger(initial,c,1.,A,B,C,320)
        cap['quadrature_refinement']=dict(orders=[240,320],half_energy_spread=abs(cap['plus_cap_proper_matter_energy']-cap_check['plus_cap_proper_matter_energy']),
            pressure_work_spread=abs(cap['plus_cap_pressure_work_rate']-cap_check['plus_cap_pressure_work_rate']))
    after=imported_sources()
    if before!=after:raise RuntimeError('directional initial-data sources changed during calculation')
    directory.mkdir(parents=True,exist_ok=True);raw=directory/'initial-data.npz'
    np.savez_compressed(raw,coefficients=c,degrees=labels,A=A,B=B,York_C=C,kick_amplitude=1.,H=H,
        vacuum=p.vacuum,mass_squared=p.mass_squared,kappa=p.kappa,offset=p.offset,
        initial_neck_coefficients=initial.coefficients,metadata=json.dumps(dict(source_sha256=before,
            preparation='chi=A*u,f=B*sqrt(3/5)*(u²-1/6),referenceL2(f)=1,trK0,TT0')))
    norm=Fraction(3,5)*(Fraction(1,8)-Fraction(1,12)+Fraction(1,36))
    return dict(schema='nsc-atlas-valve-directional-initial-v1',
        provenance=dict(source_sha256_before=before,source_sha256_after=after,source_unchanged=True,
            initial_neck_coefficients_sha256=hashlib.sha256(initial.coefficients.tobytes()).hexdigest(),
            git_parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()),
        model_inputs=dict(H_reference=H,vacuum=p.vacuum,scalar_mass_squared=p.mass_squared,
            kappa=p.kappa,offset=p.offset),
        units='Retained action units; material and work rates are model energy per initial unit-lapse normal time, not SI predictions.',
        preparation=dict(A=A,B=B,York_C=C,kick_amplitude=1.,chi='A*u,u=X1',
            momentum='Pi=lambda*B*sqrt(3/5)*(u²-1/6)*psi^-6',
            York='Ar=C*u*(1-u²)/6,At=-Ar/2,Kphysical=lambda*psi^-6*Abar',
            same_action=True,S=0.,Pi_S=0.,traceK=0.,conformal_TT=0.,
            direction_supplied=True,amplitude_supplied=True,
            no_old_reflection_parity_applied=True),
        exact_controls=dict(normalized_sphere_f_squared_over_B_squared=str(norm),old_profile_value='1/24',
            pure_source='jbar=-(kappa*A*B*sqrt(3/5)/3)gradbar(u³-u/2),purel3',
            York_inverse='Qj=-16H²j,W=-j/(16H²)',
            CKV_compatibility='All6Killing and4gradient-l1 charges vanish by harmonic orthogonality.'),
        radial_degree=40,lambda_path=list(path),stages=stages,independent_validation=validation,initial_cap_ledger=cap,
        raw_cache=dict(path=raw.relative_to(ROOT).as_posix(),sha256=sha(raw)),
        numerical=dict(threads=4,elapsed_seconds=time.perf_counter()-started,
            scope='One positive initial-data continuation and fixed higher quadrature; finite numerical estimates, not a continuum existence certificate.'),
        verdict=dict(positive_initial_branch_completed=True,nonzero_fixed_plus_cap_material_influx=cap['plus_cap_material_influx']>0,
            evolution_computed=False,inheritance_selection_proved=False,horizon_or_child_proved=False,
            public_summary_written=False,
            warning='chi odd andPi even: the old reflection-projected RHS cannot be silently used for any later evolution.'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-directory',type=Path,default=ROOT/'.build/valve-directional-v1')
    parser.add_argument('--summary',type=Path)
    args=parser.parse_args()
    directory=(ROOT/args.output_directory).resolve()
    summary=(ROOT/args.summary).resolve() if args.summary is not None else None
    if summary is not None and (not summary.is_relative_to(ROOT/'results/summaries') or summary.exists()):
        raise ValueError('fresh public summary under results/summaries required')
    result=calculation(directory)
    result['verdict']['public_summary_written']=summary is not None
    command=['python','scripts/run_python.py',Path(__file__).relative_to(ROOT).as_posix(),
        '--output-directory',directory.relative_to(ROOT).as_posix()]
    if summary is not None:command.extend(['--summary',summary.relative_to(ROOT).as_posix()])
    result['reproduction']=subprocess.list2cmdline(command)
    output=directory/'result.json'
    content=json.dumps(result,indent=2,allow_nan=False)+'\n'
    output.write_text(content,encoding='utf-8',newline='\n')
    if summary is not None:
        summary.parent.mkdir(parents=True,exist_ok=True)
        summary.write_text(content,encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(output),last=result['stages'][-1],validation=result['independent_validation'],cap=result['initial_cap_ledger']),indent=2))


if __name__=='__main__':main()
