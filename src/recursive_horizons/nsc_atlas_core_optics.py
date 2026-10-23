"""Nonradial light and curvature focusing on the retained evolving parent.

The field stepper is imported unchanged. Passive photons and Jacobi fields
are streamed between its states; no frozen lens, fitted density profile or
new stress component is substituted for the evolved scalar.
"""
from dataclasses import dataclass
import hashlib
import math
from pathlib import Path
import time

import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .nsc_atlas_core_evolution import (
    initial_state,moving_mesh,evolve_arrays,geometry,interpolate,rotate,rotate_bg)
from .nsc_atlas_core_reservoir import CoreReservoir


@njit(cache=True)
def optical_values(R,f,g0,gm,g1,b0,bm,b1,omega):
    weights=np.array([2*(f-.5)*(f-1),4*f*(1-f),2*f*(f-.5)])
    # alpha, B, alpha_R, B_R, M, scalar delta, delta Pi, F.
    values=np.zeros(8);bg=np.zeros(4)
    for j,g in enumerate((g0,gm,g1)):
        for q,index in enumerate((6,4,7,5,3,10,11,0)):
            values[q]+=weights[j]*interpolate(g[index],R,g[8])
    for j,b in enumerate((b0,bm,b1)):bg+=weights[j]*b
    A,B,Ar,Br,M,delta,dp,F=values
    Pi=bg[1]+dp;v=bg[3]*R;beta=math.sqrt(B)
    dkin=dp*(dp+2*bg[1])/2
    dpot=omega*omega*delta*(delta+2*bg[0])/2
    # Constant vacuum cancels algebraically in the Weyl scalar.
    W=M/(R**3)-(dkin+dpot-B*F*F/2)/6
    return A,B,Ar,Br,Pi,F,v,W,beta


@njit(cache=True)
def ray_rates(y,L,f,g0,gm,g1,b0,bm,b1,omega):
    R,phi,p,Dp,Pp,Dl,Pl,logE,affine=y
    A,B,Ar,Br,Pi,F,v,W,beta=optical_values(R,f,g0,gm,g1,b0,bm,b1,omega)
    E=math.sqrt(B*p*p+L*L/(R*R));q=beta*p;mu=q/E
    H=v/R;Kr=-H-R*Pi*F/2
    ricci=(E*Pi+q*beta*F)**2/2
    weyl=3*W*L*L/(R*R)
    normal_rate=A*(Kr*mu*mu-H*(1-mu*mu))-mu*beta*Ar
    return np.array([A*(v+B*p/E),A*L/(E*R*R),
        -Ar*(E+v*p)-A*((Br*p*p-2*L*L/R**3)/(2*E)+H*p),
        A*Pp/E,-A*(ricci+weyl)*Dp/E,
        A*Pl/E,-A*(ricci-weyl)*Dl/E,normal_rate,A/E])


@njit(cache=True)
def ray_step(y,L,f0,f1,h,g0,gm,g1,b0,bm,b1,omega):
    dt=h*(f1-f0);mid=(f0+f1)/2
    k1=ray_rates(y,L,f0,g0,gm,g1,b0,bm,b1,omega)
    k2=ray_rates(y+dt*k1/2,L,mid,g0,gm,g1,b0,bm,b1,omega)
    k3=ray_rates(y+dt*k2/2,L,mid,g0,gm,g1,b0,bm,b1,omega)
    k4=ray_rates(y+dt*k3,L,f1,g0,gm,g1,b0,bm,b1,omega)
    return y+dt*(k1+2*k2+2*k3+k4)/6


@njit(cache=True)
def stream_rays(delta,Pi,bg,dr,dt,steps,block_steps,mesh,a0,screen,epochs,angles,initial_clock,start):
    omega=2000.;vacuum=9/32;nray=len(angles)
    rays=np.zeros((nray,9));constants=np.zeros(nray);status=np.zeros(nray,np.int64)
    output=np.zeros((nray,16));minimum=np.full(nray,screen)
    clock=initial_clock;g=geometry(delta,Pi,bg,dr,omega,1.,mesh,a0)+(delta,Pi)
    nblocks=(steps+block_steps-1)//block_steps
    for block in range(nblocks):
        count=min(block_steps,steps-block*block_steps);h=count*dt;t=start+block*block_steps*dt
        oldd=delta;oldp=Pi;oldbg=bg;oldg=g;oldclock=clock
        out=evolve_arrays(delta,Pi,bg,dr,dt,count,omega,1.,vacuum,screen,np.empty(0),count,mesh,a0,True,0.)
        delta,Pi,bg=out[0],out[1],out[2];clock+=out[3][1]
        g=geometry(delta,Pi,bg,dr,omega,1.,mesh,a0)+(delta,Pi)
        hx,hp=rotate(oldd,oldp,h/2,omega);ex,ep=rotate(delta,Pi,-h/2,omega)
        md=(hx+ex)/2;mp=(hp+ep)/2;mb=(rotate_bg(oldbg,h/2,omega)+rotate_bg(bg,-h/2,omega))/2
        mg=geometry(md,mp,mb,dr,omega,1.,mesh,a0)+(md,mp)
        for j in range(nray):
            if status[j]==2:continue
            f0=0.
            if status[j]==0:
                if clock<epochs[j]:continue
                # Screen clock accumulated by the existing fourth-order solver;
                # one partial emission event is located within this finite block.
                f0=max(0.,min(1.,(epochs[j]-oldclock)/(clock-oldclock)))
                A,B,Ar,Br,pf,F,v,W,beta=optical_values(screen,f0,oldg,mg,g,oldbg,mb,bg,omega)
                w=-v/beta;gamma=1/math.sqrt(1-w*w);mu=-math.cos(angles[j])
                E=gamma*(1+w*mu);q=gamma*(mu+w)
                constants[j]=screen*math.sin(angles[j])
                rays[j]=np.array([screen,0.,q/beta,0.,1.,0.,1.,math.log(E),0.])
                output[j,13]=t+f0*h;output[j,14]=E;output[j,15]=angles[j]
                status[j]=1
            oldray=rays[j].copy();newray=ray_step(oldray,constants[j],f0,1.,h,oldg,mg,g,oldbg,mb,bg,omega)
            minimum[j]=min(minimum[j],newray[0]);rays[j]=newray
            if newray[0]>=screen and oldray[0]<screen:
                lo=f0;hi=1.
                # The event is only one block long; bisection of the same
                # geodesic step avoids substituting equal coordinate times.
                for trial in range(22):
                    f=(lo+hi)/2
                    y=ray_step(oldray,constants[j],f0,f,h,oldg,mg,g,oldbg,mb,bg,omega)
                    if y[0]<screen:lo=f
                    else:hi=f
                f=(lo+hi)/2;y=ray_step(oldray,constants[j],f0,f,h,oldg,mg,g,oldbg,mb,bg,omega)
                A,B,Ar,Br,pf,F,v,W,beta=optical_values(screen,f,oldg,mg,g,oldbg,mb,bg,omega)
                E=math.sqrt(B*y[2]**2+constants[j]**2/screen**2);q=beta*y[2]
                gamma=1/math.sqrt(1-v*v/B)
                observed=gamma*(E+v*y[2]);mu_observed=(q+v*E/beta)/(E+v*y[2])
                independent=gamma*(math.exp(y[7])+v*y[2])
                Dsym=screen*math.sin(y[1])/math.sin(angles[j])
                output[j,:13]=np.array([t+f*h,oldclock+f*(clock-oldclock),y[1],observed,
                    mu_observed,y[3],y[5],y[3]*y[5],y[3]-Dsym,independent-observed,
                    minimum[j],y[8],Dsym])
                status[j]=2
        if np.all(status==2):
            # Finish the declared field interval to authenticate the streamed
            # evolution against the retained final-state hash.
            remaining=steps-(block+1)*block_steps
            if remaining>0:
                out=evolve_arrays(delta,Pi,bg,dr,dt,remaining,omega,1.,vacuum,screen,np.empty(0),remaining,mesh,a0,True,0.)
                delta,Pi,bg=out[0],out[1],out[2]
            return output,status,delta,Pi,bg
    return output,status,delta,Pi,bg


@dataclass(frozen=True)
class OpticalRun:
    kind: str = 'diffuse'
    spacing: float = .0005
    step: float = .00005
    block_steps: int = 1
    outer: float = 4.
    screen: float = .6
    prefix: float = 2.5
    end: float = 5.


def angular_fan():
    fixed=(.02,.04,.08);delta=.0002
    angles=[.004,.006]+list(np.linspace(.008,.12,17))
    for x in fixed:angles.extend((x,x-delta,x+delta,x-delta/2,x+delta/2))
    angles=np.unique(np.round(angles,12))
    return np.tile(angles,2),np.repeat((2.8,3.4),len(angles))


def run_optics(config=OpticalRun(),*,cache_dir=None):
    c=config;cpu=time.process_time();wall=time.perf_counter()
    if c.kind not in ('diffuse','homogeneous') or c.block_steps not in (1,2,4):raise ValueError('unadmitted optical run')
    r=np.arange(round(c.outer/c.spacing)+1)*c.spacing
    delta,Pi,bg,model,note=initial_state(c.kind,r);a0=bg[2];mesh=moving_mesh(r)
    prefix=evolve_arrays(delta,Pi,bg,c.spacing,c.step,round(c.prefix/c.step),2000.,1.,9/32,c.screen,
        np.empty(0),round(c.prefix/c.step),mesh,a0,True,0.)
    delta,Pi,bg=prefix[0],prefix[1],prefix[2];clock=prefix[3][1]
    if cache_dir:
        folder=Path(cache_dir);folder.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(folder/'prefix.npz',delta=delta,Pi=Pi,background=bg,screen_clock=clock)
    angles,epochs=angular_fan()
    out,status,delta,Pi,bg=stream_rays(delta,Pi,bg,c.spacing,c.step,round((c.end-c.prefix)/c.step),
        c.block_steps,mesh,a0,c.screen,epochs,angles,clock,c.prefix)
    if not np.all(status==2):raise RuntimeError('rays did not return to the screen: '+str(np.where(status!=2)[0]))
    radius=r+(bg[2]/a0-1)*mesh[0]
    digest=hashlib.sha256(np.asarray(np.column_stack((radius,delta,Pi)),dtype='<f8').tobytes()).hexdigest()
    results=[]
    for j,values in enumerate(out):
        results.append(dict(emission_proper_age=float(epochs[j]),sky_angle=float(angles[j]),
            arrival_time=float(values[0]),arrival_proper_age=float(values[1]),landing_angle=float(values[2]),
            received_frequency=float(values[3]),received_radial_cosine=float(values[4]),
            jacobi_perpendicular=float(values[5]),jacobi_parallel=float(values[6]),signed_beam_area=float(values[7]),
            rotation_jacobi_gap=float(values[8]),frequency_transport_gap=float(values[9]),
            minimum_radius=float(values[10]),affine_distance=float(values[11]),rotation_jacobi=float(values[12]),
            emission_time=float(values[13]),emission_normal_energy=float(values[14])))
    summary=dict(config=c.__dict__,preparation=note,final_field_sha256=digest,rays=results,
                 cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
    if cache_dir:
        import json
        (Path(cache_dir)/'completed.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    return summary


def homogeneous_screen(epochs=(2.8,3.4),angles=None,screen=.6):
    """Independent closed-FLRW great circles, clocks and area distance."""
    if angles is None:angles=np.unique(angular_fan()[0])
    model=CoreReservoir();s=model.reference;a0=s.cap_geometry['parent_maximum_radius']
    def rhs(t,y):
        chi,p,N,H,eta,tau=y;rates=s.parent_rates(chi,p,N,H);a=a0*math.exp(N)
        return (*rates[:4],rates[5],math.sqrt(1-(H*screen)**2/(1-screen**2/a**2)))
    sol=solve_ivp(rhs,(0,5),(model.ambient,0.,0.,model.H,0.,0.),method='DOP853',
        rtol=2e-12,atol=(2e-15,2e-12,2e-13,2e-13,2e-14,2e-13),dense_output=True,max_step=.002)
    if not sol.success:raise RuntimeError(sol.message)
    result=[]
    for age in epochs:
        te=brentq(lambda t:sol.sol(t)[5]-age,0,5,xtol=1e-13);ye=sol.sol(te)
        ae=a0*math.exp(ye[2]);theta=math.asin(screen/ae);be=math.sqrt(1-screen**2/ae**2)
        ve=ye[3]*screen;w=-ve/be;gamma=1/math.sqrt(1-w*w)
        for psi in angles:
            E=gamma*(1-w*math.cos(psi));q=gamma*(-math.cos(psi)+w);mu=q/E
            tangent=math.sin(psi)/E;L=screen*math.sin(psi)
            def at(t):
                y=sol.sol(t);a=a0*math.exp(y[2]);distance=y[4]-ye[4]
                z=math.cos(theta)*math.cos(distance)-math.sin(theta)*mu*math.sin(distance)
                X=math.sin(theta)*math.cos(distance)+math.cos(theta)*mu*math.sin(distance)
                Y=tangent*math.sin(distance)
                return a*math.sqrt(max(0.,1-z*z)),math.atan2(Y,X),y,a,distance
            to=brentq(lambda t:at(t)[0]-screen,te+.1,5,xtol=1e-12)
            R,phi,y,a,distance=at(to);Eo=E*ae/a;qo=math.sqrt(Eo*Eo-L*L/screen**2)
            bo=math.sqrt(1-screen**2/a**2);vo=y[3]*screen;go=1/math.sqrt(1-vo*vo/bo**2)
            received=go*(Eo+vo*qo/bo);muo=(qo+vo*Eo/bo)/(Eo+vo*qo/bo)
            D=a*math.sin(distance)/E
            result.append(dict(emission_proper_age=float(age),sky_angle=float(psi),emission_time=te,
                arrival_time=to,arrival_proper_age=float(y[5]),landing_angle=phi,received_frequency=received,
                received_radial_cosine=muo,jacobi_perpendicular=D,jacobi_parallel=D,signed_beam_area=D*D))
    return result
