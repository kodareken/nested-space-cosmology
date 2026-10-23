"""Finite spherical Einstein--Klein--Gordon evolution in an expanding parent.

The angular extrinsic curvature is prescribed as -H_parent(t); the lapse
and radial extrinsic curvature respond to the local field. No envelope
averaging or frozen gravitational potential is used in this evolution.
"""
from dataclasses import dataclass
from functools import cached_property
import math
import time

import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .nsc_atlas_core_reservoir import CoreReservoir


@dataclass
class DiffuseReservoir(CoreReservoir):
    width: float = .03

    def cloud_shape(self,R):
        r=np.asarray(R);z=np.clip((r/self.radius-.75)/.25,0,1)
        # Same C2 taper as the retained core, at the same physical support.
        w=1-(6*z**5-15*z**4+10*z**3)
        wp=-30*z*z*(1-z)**2/(.25*self.radius)
        g=np.exp(-r*r/(2*self.width*self.width))
        return g*w,g*(wp-r*w/self.width**2)

    @cached_property
    def cloud_amplitude(self):
        target=self.core.mass_solution(self.H)[0].y[0,-1]
        k=self.reference.scalar.kappa;omega=self.frequency
        def mass(A):
            def rhs(R,y):
                if R==0:return (0.,)
                f,F=self.cloud_shape(R);f*=A;F*=A
                B=1-self.reference.cap_geometry['parent_curvature']*R*R-2*y[0]/R
                return (k*R*R*(B*F*F/4+omega*omega*self.ambient*f/2+omega*omega*f*f/4),)
            out=solve_ivp(rhs,(0,self.radius),(0.,),method='DOP853',rtol=2e-11,atol=2e-15,max_step=self.radius/100)
            if not out.success:raise RuntimeError(out.message)
            return out.y[0,-1]
        return brentq(lambda A:mass(A)-target,0,.02,xtol=1e-14)

    def fields(self,R,depletion):
        f,F,out=super().fields(R,depletion)
        cf,cF=self.core.field(R);g,gp=self.cloud_shape(R)
        return f-cf+self.cloud_amplitude*g,F-cF+self.cloud_amplitude*gp,out


@njit(cache=True)
def derivatives(u,dr):
    n=len(u);first=np.zeros(n);second=np.zeros(n)
    c1=(.8,-.2,4/105,-1/280)
    c2=(1.6,-.2,8/315,-1/560)
    for i in range(n-4):
        val=-205*u[i]/72
        for j in range(1,5):
            left=u[abs(i-j)];right=u[i+j]
            first[i]+=c1[j-1]*(right-left)/dr
            val+=c2[j-1]*(right+left)
        second[i]=val/(dr*dr)
    return first,second


@njit(cache=True)
def midpoint(u,i):
    n=len(u)
    if i==0:return (5*u[0]+15*u[1]-5*u[2]+u[3])/16
    if i==n-2:return (u[-4]-5*u[-3]+15*u[-2]+5*u[-1])/16
    return (-u[i-1]+9*u[i]+9*u[i+1]-u[i+2])/16


@njit(cache=True)
def geometry(delta,Pi,bg,dr,omega,kappa,mesh,a0,backreaction=True):
    """Radial mass and angular-evolution lapse equations, RK4 in radius."""
    chi0,P0,a,H=bg
    n=len(delta);F,Fp=derivatives(delta,dr);Pr,_=derivatives(Pi,dr)
    scale=a/a0
    Rgrid=np.arange(n)*dr+(scale-1)*mesh[0]
    J=1+(scale-1)*mesh[1];Jx=(scale-1)*mesh[2]
    Fp=Fp/(J*J)-F*Jx/(J*J*J);F=F/J;Pr=Pr/J
    if not backreaction:
        return F,Fp,Pr,np.zeros(n),1-Rgrid**2/a**2,-2*Rgrid/a**2,np.ones(n),np.zeros(n),Rgrid,J
    M=np.zeros(n);B=np.ones(n);Q=np.zeros(n);P=np.zeros(n)
    coef=np.zeros(n);source=np.zeros(n)
    Hdot=1/(a*a)-kappa*P0*P0/2
    for i in range(n):
        R=Rgrid[i];v=H*R;fullP=P0+Pi[i]
        base=1-R*R/(a*a)
        dkin=Pi[i]*(Pi[i]+2*P0)/2
        dpot=omega*omega*delta[i]*(delta[i]+2*chi0)/2
        coef[i]=-kappa*R*F[i]*F[i]*J[i]/2
        source[i]=kappa*R*R*(dkin+base*F[i]*F[i]/2+dpot-v*fullP*F[i])*J[i]/2
    for i in range(n-1):
        c=midpoint(coef,i);s=midpoint(source,i)
        k1=coef[i]*M[i]+source[i]
        k2=c*(M[i]+dr*k1/2)+s
        k3=c*(M[i]+dr*k2/2)+s
        k4=coef[i+1]*(M[i]+dr*k3)+source[i+1]
        M[i+1]=M[i]+dr*(k1+2*k2+2*k3+k4)/6
    Br=np.zeros(n)
    for i in range(1,n):
        R=Rgrid[i];B[i]=1-R*R/(a*a)-2*M[i]/R
        if B[i]<=0:raise ValueError('radial chart B lost positivity')
        dkin=Pi[i]*(Pi[i]+2*P0)/2
        dpot=omega*omega*delta[i]*(delta[i]+2*chi0)/2
        dpress=dkin+B[i]*F[i]*F[i]/2-dpot
        Q[i]=(M[i]/(R*R)+kappa*R*dpress/2)/B[i]
        P[i]=Q[i]-Hdot*R/B[i]
        Mp=(coef[i]*M[i]+source[i])/J[i]
        Br[i]=-2*R/(a*a)-2*Mp/R+2*M[i]/(R*R)
    P=P*J;Q=Q*J
    lapse_delta=np.zeros(n)
    for i in range(n-2,-1,-1):
        p=midpoint(P,i);q=midpoint(Q,i);z=lapse_delta[i+1]
        k1=P[i+1]*z+Q[i+1]
        k2=p*(z-dr*k1/2)+q
        k3=p*(z-dr*k2/2)+q
        k4=P[i]*(z-dr*k3)+Q[i]
        lapse_delta[i]=z-dr*(k1+2*k2+2*k3+k4)/6
    alpha=1+lapse_delta;alphar=(P*lapse_delta+Q)/J
    if np.min(alpha)<=0:raise ValueError('lapse lost positivity')
    return F,Fp,Pr,M,B,Br,alpha,alphar,Rgrid,J


@njit(cache=True)
def interpolate(u,R,grid):
    i=int(np.searchsorted(grid,R))-1
    if R<=0:return u[0]
    if i<1:i=1
    if i>=len(u)-2:return u[-1]
    result=0.
    for j in range(i-1,i+3):
        weight=1.
        for k in range(i-1,i+3):
            if j!=k:weight*=(R-grid[k])/(grid[j]-grid[k])
        result+=weight*u[j]
    return result


@njit(cache=True)
def rates(delta,Pi,bg,aux,dr,omega,kappa,vacuum,receiver,cuts,mesh,a0,backreaction,emission_age):
    F,Fp,Pr,M,B,Br,alpha,alphar,grid,J=geometry(delta,Pi,bg,dr,omega,kappa,mesh,a0,backreaction)
    k_local=kappa if backreaction else 0.
    chi0,P0,a,H=bg;n=len(delta)
    dchi=np.zeros(n);dPi=np.zeros(n)
    bgdot=np.array([P0,-3*H*P0-omega*omega*chi0,a*H,1/(a*a)-kappa*P0*P0/2])
    for i in range(n-4):
        R=grid[i];v=H*R;p=P0+Pi[i];u=H*(a/a0)*mesh[0,i]
        dchi[i]=alpha[i]*Pi[i]+(u-alpha[i]*v)*F[i]+(alpha[i]-1)*P0
        lap=3*Fp[0] if i==0 else B[i]*(Fp[i]+2*F[i]/R)
        dPi[i]=alpha[i]*(lap+(Br[i]/2+B[i]*alphar[i]/alpha[i])*F[i]
            -3*H*Pi[i]-k_local*R*p*p*F[i]/2-omega*omega*delta[i])+ (u-alpha[i]*v)*Pr[i]
        dPi[i]+=(alpha[i]-1)*bgdot[1]
    da=np.zeros(len(aux))
    da[0]=alpha[0]  # central proper clock
    A=interpolate(alpha,receiver,grid);bb=interpolate(B,receiver,grid);v=H*receiver
    da[1]=A*math.sqrt(1-v*v/bb)  # receiver held at fixed areal radius
    R=aux[2]
    A=interpolate(alpha,R,grid);bb=interpolate(B,R,grid);fp=interpolate(F,R,grid)
    p=P0+interpolate(Pi,R,grid);Ar=interpolate(alphar,R,grid);v=H*R
    Kr=-H-k_local*R*p*fp/2
    da[2]=A*(v+math.sqrt(bb));da[3]=A*Kr-math.sqrt(bb)*Ar
    # Independent canonical null-geodesic momentum: Hamiltonian C*p_R.
    da[-2]=-Ar*(v+math.sqrt(bb))-A*(H+interpolate(Br,R,grid)/(2*math.sqrt(bb)))
    if aux[0]<emission_age:da[2]=0.;da[3]=0.;da[-2]=0.
    da[-1]=1/a  # reference conformal distance
    # Independent integrated geometric mass evolution at fixed radii.
    for j in range(len(cuts)):
        R=cuts[j];A=interpolate(alpha,R,grid);bb=interpolate(B,R,grid)
        fp=interpolate(F,R,grid);dp=interpolate(Pi,R,grid);p=P0+dp;v=H*R
        da[4+j]=A*kappa*R*R*((bb+v*v)*p*fp-v*(dp*(dp+2*P0)+bb*fp*fp))/2
        da[4+j]-=(A-1)*kappa*H*R**3*P0*P0/2
        if not backreaction:da[4+j]=0.
    # Lawson integrating factor: remove the exactly solvable frequency-omega
    # rotation, retaining lapse, gravity, gradients and expansion in N(y).
    dchi-=Pi;dPi+=omega*omega*delta
    bgdot[0]=0.;bgdot[1]+=omega*omega*chi0
    return dchi,dPi,bgdot,da


@njit(cache=True)
def rotate(x,p,h,omega):
    c=math.cos(omega*h);s=math.sin(omega*h)
    return c*x+s*p/omega,c*p-s*omega*x


@njit(cache=True)
def rotate_bg(bg,h,omega):
    out=bg.copy();out[0],out[1]=rotate(bg[0],bg[1],h,omega)
    return out


@njit(cache=True)
def event_photon_rates(y,f,g0,g1,bg0,bg1,kappa,backreaction):
    """One partial emission step, on linearly interpolated metric data."""
    R=y[0];fields=np.zeros(5)
    # alpha, B, alpha_R, B_R, and K_r, each interpolated independently.
    for j,index in enumerate((6,4,7,5)):
        fields[j]=(1-f)*interpolate(g0[index],R,g0[8])+f*interpolate(g1[index],R,g1[8])
    H=(1-f)*bg0[3]+f*bg1[3];Kr=-H
    if backreaction:
        p0=bg0[1]+interpolate(g0[10],R,g0[8]);p1=bg1[1]+interpolate(g1[10],R,g1[8])
        F0=interpolate(g0[0],R,g0[8]);F1=interpolate(g1[0],R,g1[8])
        Kr-=kappa*R*((1-f)*p0*F0+f*p1*F1)/2
    A,B,Ar,Br=fields[:4];beta=math.sqrt(B);v=H*R
    return np.array([A*(v+beta),A*Kr-beta*Ar,-Ar*(v+beta)-A*(H+Br/(2*beta))])


@njit(cache=True)
def evolve_arrays(delta,Pi,bg,dr,dt,steps,omega,kappa,vacuum,receiver,cuts,stride,mesh,a0,backreaction,emission_age):
    n=len(delta);aux=np.zeros(6+len(cuts));arrival=np.zeros(8)
    snapshots=np.zeros(((steps+stride-1)//stride,13+2*len(cuts)));count=0
    geom=geometry(delta,Pi,bg,dr,omega,kappa,mesh,a0,backreaction)
    initial_mass=np.array([interpolate(geom[3],x,geom[8]) for x in cuts])
    max_flux_error=np.zeros(len(cuts));minB=1.;minAlpha=1.;maxspeed=0.
    for it in range(steps):
        oldaux=aux.copy();oldbg=bg.copy();old_delta=delta;old_Pi=Pi
        k1=rates(delta,Pi,bg,aux,dr,omega,kappa,vacuum,receiver,cuts,mesh,a0,backreaction,emission_age)
        x2,p2=rotate(delta+dt*k1[0]/2,Pi+dt*k1[1]/2,dt/2,omega)
        b2=rotate_bg(bg+dt*k1[2]/2,dt/2,omega)
        k2=rates(x2,p2,b2,aux+dt*k1[3]/2,dr,omega,kappa,vacuum,receiver,cuts,mesh,a0,backreaction,emission_age)
        ehx,ehp=rotate(delta,Pi,dt/2,omega);ehb=rotate_bg(bg,dt/2,omega)
        k3=rates(ehx+dt*k2[0]/2,ehp+dt*k2[1]/2,ehb+dt*k2[2]/2,aux+dt*k2[3]/2,dr,omega,kappa,vacuum,receiver,cuts,mesh,a0,backreaction,emission_age)
        ex,ep=rotate(delta,Pi,dt,omega);eb=rotate_bg(bg,dt,omega)
        q3x,q3p=rotate(k3[0],k3[1],dt/2,omega);q3b=rotate_bg(k3[2],dt/2,omega)
        k4=rates(ex+dt*q3x,ep+dt*q3p,eb+dt*q3b,aux+dt*k3[3],dr,omega,kappa,vacuum,receiver,cuts,mesh,a0,backreaction,emission_age)
        q1x,q1p=rotate(k1[0],k1[1],dt,omega);q1b=rotate_bg(k1[2],dt,omega)
        q2x,q2p=rotate(k2[0],k2[1],dt/2,omega);q2b=rotate_bg(k2[2],dt/2,omega)
        delta=ex+dt*(q1x+2*q2x+2*q3x+k4[0])/6
        Pi=ep+dt*(q1p+2*q2p+2*q3p+k4[1])/6
        bg=eb+dt*(q1b+2*q2b+2*q3b+k4[2])/6
        aux+=dt*(k1[3]+2*k2[3]+2*k3[3]+k4[3])/6
        if not np.isfinite(delta[0]):raise ValueError('nonfinite scalar evolution')
        if oldaux[0]<emission_age and aux[0]>=emission_age:
            g0=geometry(old_delta,old_Pi,oldbg,dr,omega,kappa,mesh,a0,backreaction)
            g1=geometry(delta,Pi,bg,dr,omega,kappa,mesh,a0,backreaction)
            d=aux[0]-oldaux[0];z=(emission_age-oldaux[0])/d
            aa=dt*g0[6][0];bb=3*d-dt*(2*g0[6][0]+g1[6][0]);cc=-2*d+dt*(g0[6][0]+g1[6][0])
            for trial in range(5):
                val=oldaux[0]+aa*z+bb*z*z+cc*z**3-emission_age
                z=min(1.,max(0.,z-val/(aa+2*bb*z+3*cc*z*z)))
            h=dt*(1-z);y=np.zeros(3)
            q0=g0+(old_Pi,);q1=g1+(Pi,)
            p1=event_photon_rates(y,z,q0,q1,oldbg,bg,kappa,backreaction)
            p2=event_photon_rates(y+h*p1/2,(z+1)/2,q0,q1,oldbg,bg,kappa,backreaction)
            p3=event_photon_rates(y+h*p2/2,(z+1)/2,q0,q1,oldbg,bg,kappa,backreaction)
            p4=event_photon_rates(y+h*p3,1.,q0,q1,oldbg,bg,kappa,backreaction)
            y+=h*(p1+2*p2+2*p3+p4)/6
            aux[2]=y[0];aux[3]=y[1];aux[-2]=y[2];arrival[7]=(it+z)*dt
        if arrival[0]==0 and aux[2]>=receiver:
            f=(receiver-oldaux[2])/(aux[2]-oldaux[2])
            geom=geometry(delta,Pi,bg,dr,omega,kappa,mesh,a0,backreaction)
            oldgeom=geometry(old_delta,old_Pi,oldbg,dr,omega,kappa,mesh,a0,backreaction)
            B0=interpolate(oldgeom[4],receiver,oldgeom[8]);B1=interpolate(geom[4],receiver,geom[8])
            b=math.sqrt(B0+f*(B1-B0));v=(oldbg[3]+f*(bg[3]-oldbg[3]))*receiver
            lognormal=oldaux[3]+f*(aux[3]-oldaux[3])
            arrival[0]=(it+f)*dt;arrival[1]=oldaux[1]+f*(aux[1]-oldaux[1])
            arrival[2]=math.exp(lognormal)*math.sqrt((b+v)/(b-v))
            arrival[3]=oldaux[0]+f*(aux[0]-oldaux[0]);arrival[4]=lognormal
            arrival[5]=v/b
            logp=oldaux[-2]+f*(aux[-2]-oldaux[-2])
            arrival[6]=math.exp(logp)*b*math.sqrt((b+v)/(b-v))
        if (it+1)%stride==0 or it==steps-1:
            g=geometry(delta,Pi,bg,dr,omega,kappa,mesh,a0,backreaction)
            B=g[4];alpha=g[6];M=g[3];grid=g[8];J=g[9]
            row=np.zeros(13+2*len(cuts))
            row[0]=(it+1)*dt;row[1:5]=bg;row[5]=delta[0];row[6]=Pi[0]
            row[7]=alpha[0];row[8]=aux[0];row[9]=aux[1]
            row[10]=3*(bg[3]**2+1/bg[2]**2)-kappa*(vacuum+(bg[1]**2+omega*omega*bg[0]**2)/2)
            row[11]=(bg[1]+Pi[0])**2/2+omega*omega*(bg[0]+delta[0])**2/2
            for j in range(len(cuts)):
                mass=interpolate(M,cuts[j],grid);error=mass-initial_mass[j]-aux[4+j]
                row[12+j]=mass;row[12+len(cuts)+j]=error
                max_flux_error[j]=max(max_flux_error[j],abs(error))
            target=interpolate(M,.12,grid)/2
            row[-1]=np.nan
            if target>0:
                for j in range(1,min(n,int(.12/dr)+2)):
                    if M[j]>=target:
                        lo=(j-1)*dr;hi=j*dr
                        for trial in range(24):
                            mid=(lo+hi)/2
                            if interpolate(M,mid,grid)<target:lo=mid
                            else:hi=mid
                        row[-1]=(lo+hi)/2
                        break
            minB=min(minB,np.min(B));minAlpha=min(minAlpha,np.min(alpha))
            for j in range(n):maxspeed=max(maxspeed,(abs(alpha[j]*bg[3]*grid[j]-bg[3]*(bg[2]/a0)*mesh[0,j])+alpha[j]*math.sqrt(B[j]))/J[j])
            snapshots[count]=row;count+=1
    return delta,Pi,bg,aux,snapshots[:count],arrival,max_flux_error,minB,minAlpha,maxspeed,initial_mass


def moving_mesh(x):
    """Fixed physical radius below1; parent-comoving radius above2, C-infinity join."""
    z=np.clip(x-1,0,1);w=np.zeros_like(x);wp=w.copy();wpp=w.copy()
    inside=(z>0)&(z<1);q=z[inside]
    t=-1/q+1/(1-q)
    from scipy.special import expit
    b=expit(t);tp=1/q**2+1/(1-q)**2;tpp=-2/q**3+2/(1-q)**3
    w[inside]=b;wp[inside]=b*(1-b)*tp;wpp[inside]=b*(1-b)*((1-2*b)*tp*tp+tpp)
    w[z>=1]=1.
    return np.array([x*w,w+x*wp,2*wp+x*wpp])


@dataclass(frozen=True)
class EvolutionRun:
    kind: str = 'core'
    outer: float = 4.
    spacing: float = .001
    step: float = .0001
    end: float = 1.
    receiver: float = .6
    emission_age: float = 0.
    backreaction: bool = True


def initial_state(kind,radius):
    if kind not in ('core','diffuse','homogeneous'):raise ValueError('unknown initial preparation')
    model=DiffuseReservoir() if kind=='diffuse' else CoreReservoir()
    if kind=='homogeneous':
        row=None;delta=np.zeros_like(radius);Pi=delta.copy()
    else:
        solution=model.solve();row=solution.state(np.minimum(radius,model.outer))
        delta,F,out=model.fields(radius,solution.depletion)
        Pi=np.where(radius<model.outer,row['Pi'],0.)
    s=model.reference
    bg=np.array([model.ambient,0.,s.cap_geometry['parent_maximum_radius'],model.H])
    note=dict(kind=kind,depletion=0. if kind=='homogeneous' else solution.depletion,
              cloud_amplitude=model.cloud_amplitude if kind=='diffuse' else None,
              support=model.outer,frequency=model.frequency)
    return delta,Pi,bg,model,note


def run_evolution(config=EvolutionRun(),*,retain=False):
    c=config;start=time.process_time();wall=time.perf_counter()
    if not (.5<=c.outer<=5 and .0001<=c.spacing<=.002 and 0<c.step<=.0002 and 0<c.end<=6):
        raise ValueError('outside admitted finite evolution')
    if not (0<c.receiver<min(1.,c.outer) and c.step/c.spacing<=.25 and 0<=c.emission_age<c.end):
        raise ValueError('receiver must lie inside the resolved chart and step/spacing<=.25')
    n=round(c.outer/c.spacing)+1;r=np.arange(n)*c.spacing
    if not math.isclose(r[-1],c.outer,abs_tol=1e-12):raise ValueError('outer must lie on mesh')
    steps=round(c.end/c.step)
    if not math.isclose(steps*c.step,c.end,abs_tol=1e-12):raise ValueError('endpoint must lie on time mesh')
    delta,Pi,bg,model,note=initial_state(c.kind,r)
    cuts=np.array([.025,.06,.12,.24,.484,.8])
    if c.outer<=max(cuts):cuts=cuts[cuts<c.outer]
    mesh=moving_mesh(r)
    a0=bg[2]
    initial_geom=geometry(delta,Pi,bg,c.spacing,model.frequency,model.reference.scalar.kappa,mesh,a0,c.backreaction)
    initial_energy=(Pi[0]**2+model.frequency**2*(bg[0]+delta[0])**2)/2
    from scipy.optimize import brentq
    initial_half=brentq(lambda R:interpolate(initial_geom[3],R,r)-interpolate(initial_geom[3],.12,r)/2,0,.12) if c.kind!='homogeneous' and c.backreaction else None
    out=evolve_arrays(delta,Pi,bg,c.spacing,c.step,steps,model.frequency,model.reference.scalar.kappa,
        model.reference.scalar.vacuum(model.reference.parent_flux),c.receiver,cuts,max(1,round(.001/c.step)),mesh,a0,c.backreaction,c.emission_age)
    d,p,b,a,snap,arrival,error,minB,minAlpha,speed,initial_mass=out
    theta0=math.asin(note['support']/model.reference.cap_geometry['parent_maximum_radius'])
    front=b[2]*math.sin(theta0+a[-1])
    outer_radius=c.outer+(b[2]/a0-1)*mesh[0,-1]
    final_geometry=geometry(d,p,b,c.spacing,model.frequency,model.reference.scalar.kappa,mesh,a0,c.backreaction)
    result=dict(config=c.__dict__,preparation=note,arrival=dict(coordinate_time=arrival[0],receiver_proper_age=arrival[1],
        received_to_emitted_frequency=arrival[2],central_proper_age_at_arrival=arrival[3],normal_log_frequency=arrival[4],receiver_boost=arrival[5],emission_coordinate_time=arrival[7],emission_proper_age=c.emission_age),
        initial_mass_at_cuts=initial_mass.tolist(),cuts=cuts.tolist(),maximum_mass_flux_error=error.tolist(),
        minimum_B=minB,minimum_lapse=minAlpha,maximum_sampled_grid_light_speed=speed,
        initial_central_excitation=initial_energy,initial_half_mass_radius=initial_half,
        reference_outgoing_front=front,reference_front_buffer=outer_radius-front,physical_outer_radius=outer_radius,
        independent_canonical_frequency=arrival[6],geodesic_frequency_gap=arrival[6]-arrival[2],
        final_outer_mass_offset=float(final_geometry[3][-1]),
        minimum_receiver_timelike_margin=float(1-(b[3]*c.receiver)**2/interpolate(final_geometry[4],c.receiver,final_geometry[8])),
        maximum_background_constraint=float(np.max(abs(snap[:,10]))),
        final_central_excitation=float(snap[-1,11]),final_mass_at_cuts=snap[-1,12:12+len(cuts)].tolist(),
        late_mean_central_excitation=float(np.mean(snap[snap[:,0]>=.9*c.end,11])),
        late_mean_half_mass_radius=float(np.nanmean(snap[snap[:,0]>=.9*c.end,-1])) if c.kind!='homogeneous' and c.backreaction else None,
        cpu_seconds=time.process_time()-start,wall_seconds=time.perf_counter()-wall)
    if retain:return result,dict(radius=r+(b[2]/a0-1)*mesh[0],delta=d,Pi=p,background=b,history=snap)
    return result
