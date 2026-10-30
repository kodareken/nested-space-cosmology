"""Finite periodic continuation of the retained compact vortex preparation.

Dealiased Fourier Galerkin incompressible Navier--Stokes, rotational form,
SSPRK3. No force, state reset, singularity search or cosmological fit.
"""
import math
import time
import numpy as np
from scipy.fft import rfftn,irfftn,fftfreq,rfftfreq
from scipy.ndimage import map_coordinates
from .nsc_atlas_vortex_transfer import jets,exact_theta_rate,energy,dissipation


class PeriodicFlow:
    def __init__(self,n,length=4.,viscosity=.001,workers=4):
        if n<12 or n%2 or viscosity<=0:raise ValueError('even n>=12 and positive viscosity required')
        self.n=n;self.length=length;self.nu=viscosity;self.workers=workers
        f=2*np.pi*fftfreq(n,d=length/n);fz=2*np.pi*rfftfreq(n,d=length/n)
        self.k=np.array(np.broadcast_arrays(f[:,None,None],f[None,:,None],fz[None,None,:]))
        self.k2=np.sum(self.k*self.k,axis=0)
        self.inv=np.zeros_like(self.k2);np.divide(1,self.k2,out=self.inv,where=self.k2>0)
        # Strict 2/3 cut removes aliased quadratic contributions to retained modes.
        self.mask=np.all(np.abs(self.k)<2*np.pi/length*(n/3-1e-10),axis=0)

    def fft(self,u):return rfftn(u,axes=(-3,-2,-1),workers=self.workers)
    def real(self,u):return irfftn(u,s=(self.n,)*3,axes=(-3,-2,-1),workers=self.workers)

    def project(self,h):
        return (h-self.k*(np.sum(self.k*h,axis=0)*self.inv))*self.mask

    def curl(self,h):
        k=self.k
        return 1j*np.array([k[1]*h[2]-k[2]*h[1],k[2]*h[0]-k[0]*h[2],k[0]*h[1]-k[1]*h[0]])

    def prepare(self,mode=4):
        a=np.arange(self.n)*self.length/self.n-self.length/2
        x,y,z=np.meshgrid(a,a,a,indexing='ij');r=np.hypot(x,y);theta=np.arctan2(y,x)
        safe=np.maximum(r,.25)
        W=jets(safe,z,mode)[0]
        u=np.real(W*np.exp(1j*(mode*theta-2*mode*np.log(safe))))
        c=np.cos(theta);s=np.sin(theta)
        cart=np.array([u[0]*c-u[1]*s,u[0]*s+u[1]*c,u[2]])
        return self.project(self.fft(cart))

    def rhs(self,h):
        u=self.real(h);w=self.real(self.curl(h))
        cross=np.array([u[1]*w[2]-u[2]*w[1],u[2]*w[0]-u[0]*w[2],u[0]*w[1]-u[1]*w[0]])
        rate=self.project(self.fft(cross))-self.nu*self.k2*h
        loss=self.nu*self.length**3*float(np.mean(np.sum(w*w,axis=0)))
        return rate,loss

    def pressure(self,h):
        u=self.real(h);w=self.real(self.curl(h))
        cross=np.array([u[1]*w[2]-u[2]*w[1],u[2]*w[0]-u[0]*w[2],u[0]*w[1]-u[1]*w[0]])
        bernoulli=-1j*np.sum(self.k*self.fft(cross),axis=0)*self.inv
        p=bernoulli-self.fft(np.sum(u*u,axis=0)/2)
        p*=self.mask;p[0,0,0]=0
        return self.real(p)

    def mean_theta(self,h,radius,samples=256):
        u=self.real(h);theta=2*np.pi*np.arange(samples)/samples
        coords=np.array([(radius*np.cos(theta)+self.length/2)*self.n/self.length,
                         (radius*np.sin(theta)+self.length/2)*self.n/self.length])
        ux=map_coordinates(u[0,:,:,self.n//2],coords,order=3,mode='grid-wrap')
        uy=map_coordinates(u[1,:,:,self.n//2],coords,order=3,mode='grid-wrap')
        return float(np.mean(-np.sin(theta)*ux+np.cos(theta)*uy))

    def diagnostics(self,h):
        u=self.real(h);w=self.real(self.curl(h));V=self.length**3
        E=V/2*float(np.mean(np.sum(u*u,axis=0)))
        D=self.nu*V*float(np.mean(np.sum(w*w,axis=0)))
        div=self.real(1j*np.sum(self.k*h,axis=0))
        return dict(energy=E,dissipation=D,divergence_rms=float(np.sqrt(np.mean(div*div))),
            max_speed=float(np.max(np.sqrt(np.sum(u*u,axis=0)))),
            mean_theta=[self.mean_theta(h,r) for r in (.75,1.25)])

    def step(self,h,dt):
        f0,d0=self.rhs(h);h1=h+dt*f0
        f1,d1=self.rhs(h1);h2=.75*h+.25*(h1+dt*f1)
        f2,d2=self.rhs(h2)
        return h/3+2*(h2+dt*f2)/3,dt*(d0/6+d1/6+2*d2/3)


def run_case(n,dt=.001,endpoint=.05,workers=4,cpu_limit=300.):
    started=time.perf_counter();cpu=time.process_time()
    count=round(endpoint/dt)
    if count<1 or not math.isclose(count*dt,endpoint,abs_tol=1e-12):raise ValueError('endpoint must be a positive integer number of steps')
    flow=PeriodicFlow(n,workers=workers);h=flow.prepare();first=flow.diagnostics(h)
    rate,_=flow.rhs(h)
    initial_rates=[flow.mean_theta(rate,r) for r in (.75,1.25)]
    exact=[float(exact_theta_rate(r)) for r in (.75,1.25)]
    rows=[dict(t=0.,**first)];integrated=0.;max_cfl=0.
    for step in range(1,count+1):
        h,increment=flow.step(h,dt);integrated+=increment
        if step%max(1,count//5)==0 or step==count:
            d=flow.diagnostics(h);rows.append(dict(t=step*dt,**d))
            max_cfl=max(max_cfl,dt*d['max_speed']*np.max(np.sqrt(flow.k2*flow.mask)))
        if not np.isfinite(integrated) or time.process_time()-cpu>cpu_limit:
            raise RuntimeError('finite evolution failed or exceeded its declared CPU allowance')
    last=rows[-1];p=flow.pressure(h)
    return dict(n=n,dt=dt,endpoint=endpoint,steps=count,FFT_workers=workers,
        initial=first,final=last,history=rows,
        initial_mean_rate=initial_rates,exact_initial_mean_rate=exact,
        initial_rate_error=[abs(a-b) for a,b in zip(initial_rates,exact)],
        initial_energy_error=abs(first['energy']-energy(4)),
        initial_loss_error=abs(first['dissipation']-dissipation(4)),
        integrated_viscous_loss=integrated,
        energy_balance_residual=last['energy']-first['energy']+integrated,
        generated_mean=[a-b for a,b in zip(last['mean_theta'],first['mean_theta'])],
        pressure_rms=float(np.sqrt(np.mean(p*p))),
        max_advective_step_number=max_cfl,
        cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-started,
        grid_arrays_stored=False)
