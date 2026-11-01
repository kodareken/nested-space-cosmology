"""Evolve the solved directional preparation with the existing SO(3) equations.

No cap-reflection projection, prescribed incoming boundary or new matter is
introduced. Fixed-cap material supply and pressure work are passive ledgers.
Marginal spheres and the migrating minimum are evaluated independently.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import DOP853
from scipy.optimize import brentq
from scipy.special import roots_jacobi
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_passive_dop853 import PassiveDOP853Quadrature
from recursive_horizons.nsc_atlas_valve_neck import NeckInitialData
from recursive_horizons.nsc_atlas_valve_spherical import SphericalValve,gegenbauer_jets

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class DirectionalEvolution(SphericalValve):
    """Only cache fixed evaluation rules; physical RHS stays in the package."""
    def __init__(self,modes,**kwargs):
        super().__init__(modes,reflection_symmetry=False,**kwargs)
        self.rules={}

    def retain_rule(self,u):
        u=np.asarray(u,dtype=float)
        self.rules[(u.shape,u.tobytes())]=gegenbauer_jets(u,self.modes)
        return u

    def jets(self,state,u=None):
        if u is None:return super().jets(state)
        u=np.atleast_1d(u);rule=self.rules.get((u.shape,u.tobytes()))
        if rule is None:return super().jets(state,u)
        coef=self.coefficients(state)
        return tuple(coef@matrix.T for matrix in rule)


def load_initial(system,path):
    path=(ROOT/path).resolve();record=json.loads(path.read_bytes())
    if record['schema']!='nsc-atlas-valve-directional-initial-v1':raise ValueError('directional initial record required')
    if not record['provenance']['source_unchanged']:raise ValueError('unchanged initial producer required')
    for name,digest in record['provenance']['source_sha256_after'].items():
        if sha(ROOT/name)!=digest:raise ValueError('initial source changed: '+name)
    raw=(ROOT/record['raw_cache']['path']).resolve()
    if not raw.is_relative_to(ROOT/'.build') or sha(raw)!=record['raw_cache']['sha256']:
        raise ValueError('initial array identity failed')
    with np.load(raw,allow_pickle=False) as cache:
        coef=cache['coefficients'].copy()
        if not np.array_equal(cache['degrees'],np.arange(0,41,2)) or coef.shape!=(21,):
            raise ValueError('degree40 radial initial data required')
        A=float(cache['A']);B=float(cache['B']);C=float(cache['York_C']);lam=float(cache['kick_amplitude'])
        for key,value in [('H',system.H),('vacuum',system.parent.vacuum),
                          ('mass_squared',system.parent.mass_squared),('kappa',system.kappa)]:
            if not math.isclose(float(cache[key]),value,rel_tol=2e-14):raise ValueError('parent mismatch: '+key)
    prep=record['preparation']
    if any(not math.isclose(value,prep[key],rel_tol=2e-14) for key,value in [('A',A),('B',B),('York_C',C),('kick_amplitude',lam)]):
        raise ValueError('initial scalar/York amplitudes differ')
    system.initial_data=NeckInitialData(system.parent,40,1.,coef,())
    u=system.u;v=1-u*u;psi=system.initial_data.values(u)
    if np.min(psi)<=0:raise ValueError('positive initial metric required')
    y=np.zeros((8,system.modes));y[0]=psi**2/system.H
    y[2]=-lam*C*u*v*psi**-6/12
    y[3]=lam*C*u*psi**-6/4
    y[5]=A*u;y[7]=lam*B*math.sqrt(3/5)*(u*u-1/6)*psi**-6
    system.initial=y
    return record


class Measurements:
    def __init__(self,system):
        self.system=system;n=2*system.modes+7
        x,w=roots_jacobi(n,.5,0.);self.cap_u=(x+1)/2;self.cap_w=w*np.sqrt(3+x)/4
        self.cap_points=system.retain_rule(np.r_[self.cap_u,-self.cap_u])
        self.zero=system.retain_rule(np.array([0.]))
        angles=np.arange(1,3*system.modes+8)*math.pi/(3*system.modes+8)
        self.global_u=system.retain_rule(np.cos(angles))
        self.global_w=4*math.pi**2/(3*system.modes+8)*np.sin(angles)**2

    def cap_rates(self,state):
        s=self.system;z=s.pieces(state,self.cap_points);n=len(self.cap_u)
        measure=4*math.pi*np.tile(self.cap_w,2)*z['q']*z['b']**2
        density_work=z['ka']*z['pr']+2*z['kb']*z['pt']
        center=s.pieces(state,self.zero)
        flow=-4*math.pi*float(center['radius'][0])**2*float(center['j'][0])
        return np.array([flow,float(measure[:n]@density_work[:n]),float(measure[n:]@density_work[n:])])

    def caps(self,state):
        s=self.system;z=s.pieces(state,self.cap_points);n=len(self.cap_u)
        measure=4*math.pi*np.tile(self.cap_w,2)*z['q']*z['b']**2
        out=[]
        for sl in (slice(0,n),slice(n,None)):
            volume=float(np.sum(measure[sl]));energy=float(measure[sl]@z['rho'][sl])
            pressure=float(measure[sl]@((z['pr'][sl]+2*z['pt'][sl])/3))
            out.append(dict(volume=volume,energy=energy,mean_density=energy/volume,
                mean_isotropic_pressure=pressure/volume,w=pressure/energy,
                bulk_volume_rate=float(-measure[sl]@z['K'][sl])))
        return out

    def geometry(self,state):
        s=self.system;z=s.pieces(state,self.global_u);v=1-self.global_u**2
        R3=z['Rr']+2*z['Rt'];K2=z['ka']**2+2*z['kb']**2
        Ham=R3+z['K']**2-K2-2*s.kappa*z['rho']
        bu=z['yu'][0];kbu=z['yu'][2]
        M=2*((self.global_u/np.sqrt(v)-np.sqrt(v)*bu/z['b'])*(z['ka']-z['kb'])+np.sqrt(v)*kbu)/z['q']-s.kappa*z['j']
        weights=self.global_w*z['q']*z['b']**2;V=float(np.sum(weights))
        scale=s.H**2+abs(R3)+z['K']**2+K2+2*s.kappa*abs(z['rho'])
        limits=dict(H=float(max(abs(Ham)/scale)),M=float(max(abs(M)/scale)))
        # Independent electric curvature components from the actual ADM evolution.
        Er=z['Rr']+2*z['ka']*z['kb']-s.kappa*(z['grad2']+z['W'])
        Et=z['Rt']+z['ka']*z['kb']+z['kb']**2-s.kappa*z['W']
        return dict(global_max_relative_constraints=limits,
            Hamiltonian_RMS_over_H2=float(np.sqrt(weights@(Ham*Ham)/V)/s.H**2),
            momentum_RMS_over_H2=float(np.sqrt(weights@(M*M)/V)/s.H**2),
            min_radial_scale=float(min(z['q'])),min_angular_scale=float(min(z['b'])),
            sampled_max_tidal_over_H2=float(max(np.max(abs(Er)),np.max(abs(Et)))/s.H**2),
            sampled_admitted=max(limits.values())<=.02,
            admission_scope='Independent finite spatial samples and ADM-term2% scale; not a continuum error bound.')

    def minimum(self,state,previous=0.):
        s=self.system;u=np.linspace(-.98,.98,257);rr=s.pieces(state,u)['rr'];roots=[]
        for a,b,fa,fb in zip(u[:-1],u[1:],rr[:-1],rr[1:]):
            if fa>=0 and fb<0:
                root=brentq(lambda x:float(s.pieces(state,[x])['rr'][0]),a,b,xtol=2e-12)
                roots.append(root)
        if not roots:return None
        coord=min(roots,key=lambda x:abs(x-previous));z=s.pieces(state,[coord])
        R=float(z['radius'][0]);rate=-float(z['kb'][0])*R
        plus,minus=s.expansions(state,[coord])
        return dict(u=coord,radius=R,normal_radius_rate=rate,
            theta_plus_over_H=float(plus[0]/s.H),theta_minus_over_H=float(minus[0]/s.H),
            rho=float(z['rho'][0]),radial_pressure=float(z['pr'][0]),
            NEC=float(z['rho'][0]+z['pr'][0]),candidate_minimum_count=len(roots),
            scope='Migrating area minimum; its coordinate motion is not an observer worldline.')

    def trapped(self,state,anchor):
        s=self.system
        def test(x):return float(max(a[0] for a in s.expansions(state,[x])))
        if test(anchor)>=0:return None
        edges=[]
        for endpoint in (-.999999,.999999):
            coords=np.linspace(anchor,endpoint,180);values=np.maximum(*s.expansions(state,coords))
            found=np.flatnonzero(values>=0)
            if not len(found):return None
            j=int(found[0]);edges.append(brentq(test,coords[j-1],coords[j],xtol=2e-12))
        lo,hi=edges
        x,w=leggauss(2*s.modes+7);u=(lo+hi)/2+(hi-lo)/2*x
        z=s.pieces(state,u);measure=4*math.pi*(hi-lo)/2*w*np.sqrt(1-u*u)*z['q']*z['b']**2
        V=float(np.sum(measure));E=float(measure@z['rho'])
        pressure=float(measure@((z['pr']+2*z['pt'])/3))
        boundary=s.pieces(state,np.array(edges));exp=s.expansions(state,np.array(edges))
        influx=4*math.pi*float(boundary['radius'][1]**2*boundary['j'][1]-boundary['radius'][0]**2*boundary['j'][0])
        return dict(u_left=lo,u_right=hi,volume=V,energy=E,
            mean_density=E/V,mean_isotropic_pressure=pressure/V,w=pressure/E,
            material_influx=influx,pressure_work_rate=float(measure@(z['ka']*z['pr']+2*z['kb']*z['pt'])),
            marginal_expansions=[[float(exp[k][j]) for k in range(2)] for j in range(2)],
            scope='Connected sampled future-trapped band about the tracked minimum; boundary sweep is not included as material supply and no child is identified.')


def source_pins():
    package=(ROOT/'src/recursive_horizons').resolve();paths={Path(__file__).resolve()}
    for module in tuple(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if name:
            path=Path(name).resolve()
            if path.suffix=='.py' and path.is_relative_to(package):paths.add(path)
    return {path.relative_to(ROOT).as_posix():sha(path) for path in sorted(paths)}


def run(args):
    folder=(ROOT/args.output).resolve()
    if not folder.is_relative_to(ROOT/'.build') or folder.exists():raise ValueError('fresh .build output directory required')
    folder.mkdir(parents=True);started=time.perf_counter()
    with threadpool_limits(limits=args.threads):
        s=DirectionalEvolution(args.modes,budget_seconds=args.budget)
        initial=load_initial(s,args.initial_record);measure=Measurements(s)
        pins=source_pins();commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
        bound=all(subprocess.run(['git','show',commit+':'+name],cwd=ROOT,capture_output=True).stdout==
            (ROOT/name).read_bytes() for name in pins)
        record=dict(schema='nsc-atlas-valve-directional-evolution-v1',status='running',modes=args.modes,
            initial_record=str(args.initial_record).replace('\\','/'),initial_record_sha256=sha(ROOT/args.initial_record),
            source_sha256=pins,producer_commit=commit,source_bound_to_commit=bound,
            reflection_projection=False,proper_endpoint_requested=args.endpoint,
            temporal_tolerances=dict(rtol=3e-11,atol=1e-13,max_step=.002),threads=args.threads,
            budget_seconds=args.budget,snapshots=[],regular_receiving_child_proved=False,
            scope='Same scalar action; declared directional initial state, not inherited state selection or an SI/observational prediction.')
        baselines=measure.caps(s.initial);saved=[];previous_min=0.
        def save():
            record['elapsed_seconds']=time.perf_counter()-started
            tmp=folder/'record.tmp';tmp.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n');tmp.replace(folder/'run.json')
        def snapshot(t,y,integrals):
            nonlocal previous_min
            state=y.reshape(8,s.modes);diag=measure.geometry(state);caps=measure.caps(state)
            for i,cap in enumerate(caps):
                material=float(integrals[0])*(1 if i==0 else -1);work=float(integrals[i+1])
                cap.update(integrated_material_received=material,integrated_pressure_work=work,
                    energy_ledger_residual=cap['energy']-baselines[i]['energy']-material-work)
            minimum=measure.minimum(state,previous_min)
            if minimum:previous_min=minimum['u']
            band=measure.trapped(state,previous_min) if minimum and t>1e-10 else None
            center=s.pieces(state,measure.zero);ep,em=s.expansions(state,measure.zero)
            diag.update(proper_time=float(t),fixed_receiving_cap=caps[0],fixed_donor_cap=caps[1],
                fixed_interface=dict(radius=float(center['radius'][0]),
                    theta_plus_over_H=float(ep[0]/s.H),theta_minus_over_H=float(em[0]/s.H),
                    material_influx=measure.cap_rates(state)[0]),minimum=minimum,trapped_band=band)
            record['snapshots'].append(diag);saved.append((float(t),state.copy()))
            save();print(json.dumps(dict(event='checkpoint',modes=args.modes,proper=t,
                admitted=diag['sampled_admitted'],minimum=minimum,
                cap_net_received=caps[0]['integrated_material_received'],cap_energy_residual=caps[0]['energy_ledger_residual'])),flush=True)
            return diag['sampled_admitted']
        if not snapshot(0.,s.initial.ravel(),np.zeros(3)):raise ValueError('initial state not numerically admitted')
        ledgers=[PassiveDOP853Quadrature() for _ in range(3)]
        def rhs(t,y):
            derivative=s.rhs(t,y)
            for ledger,rate in zip(ledgers,measure.cap_rates(y.reshape(8,s.modes))):ledger.record(rate)
            return derivative
        s.started=started
        solver=DOP853(rhs,0.,s.initial.ravel(),args.endpoint,rtol=3e-11,atol=1e-13,max_step=.002,first_step=.0002)
        for ledger in ledgers:ledger.initialize(solver)
        targets=list(np.arange(.1,args.endpoint-1e-10,.1))+[args.endpoint]
        reason=None;last_admitted=0.;heartbeat=time.perf_counter()
        try:
            while solver.status=='running' and targets:
                if time.perf_counter()-started>args.budget:raise TimeoutError('declared directional evolution budget exhausted')
                message=solver.step()
                if solver.status=='failed':raise RuntimeError(message or 'DOP853 failed')
                for ledger in ledgers:ledger.accepted_step(solver)
                dense=solver.dense_output()
                while targets and solver.t>=targets[0]-2e-14:
                    target=float(targets.pop(0));integrals=np.array([v.dense_value(solver,target) for v in ledgers])
                    if not snapshot(target,dense(target),integrals):raise RuntimeError('sampled ADM admission ended; no physical endpoint inferred')
                    last_admitted=target
                if time.perf_counter()-heartbeat>30:
                    record['progress_time']=float(solver.t);save();heartbeat=time.perf_counter()
        except (TimeoutError,ValueError,RuntimeError,FloatingPointError) as exc:reason=str(exc)
        record.update(status='endpoint_reached' if not targets and reason is None else 'incomplete',stop_reason=reason,
            last_admitted_proper=last_admitted,actual_last_accepted_proper=float(solver.t),rhs_calls=s.calls)
        cache=folder/'sparse.npz';np.savez_compressed(cache,times=np.array([v[0] for v in saved]),states=np.array([v[1] for v in saved]))
        record.update(state_cache=cache.relative_to(ROOT).as_posix(),state_cache_sha256=sha(cache),
            source_unchanged=pins=={name:sha(ROOT/name) for name in pins})
        save();print(json.dumps(dict(event='terminal',status=record['status'],stop=reason,elapsed=record['elapsed_seconds'])),flush=True)
        return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes',type=int,required=True);parser.add_argument('--endpoint',type=float,default=2.)
    parser.add_argument('--threads',type=int,default=1);parser.add_argument('--budget',type=float,default=120.)
    parser.add_argument('--initial-record',type=Path,default=Path('results/summaries/nsc-atlas-valve-directional-v1.json'))
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if not 0<args.endpoint<=2 or not 1<=args.threads<=8 or args.budget<=0:raise ValueError('finite endpoint<=2,threads1..8 and positive budget required')
    run(args)


if __name__=='__main__':main()
