"""Named initial-state controls and the physical volume equation of state.

The HM-location-plus-matter and minimum-start experiments are conditional
preparations, not derived tunnelling exits. E0 is the actual energy at the
new starting slice; source couplings are not fitted.
"""
from dataclasses import dataclass
import math
import numpy as np
from scipy.integrate import solve_ivp
from .nsc_atlas_parent_scale_selector import VolumeParent
from .nsc_atlas_inherited_matter import exact_initial_geometry
from .nsc_atlas_curvature_mechanisms import ClosedVolumeComparison


@dataclass(frozen=True)
class PreparedModel(VolumeParent):
    initial_s: float=0.
    curvature: float=0.
    def matter(self,N,s):return .375*math.exp(-3*N-6*(s-self.initial_s))
    def hubble(self,N,s,v,E):
        H2=(self.K*v*v/2+self.potential(s)+E)/3-self.curvature*math.exp(-2*N)
        if H2<=0:raise ValueError('prepared expansion ended')
        return math.sqrt(H2)


class PreparedBackground(ClosedVolumeComparison):
    def __init__(self,zmax,*,preparation='matched',curvature='flat',rtol=2e-11):
        self.C=0. if curvature=='flat' else float(exact_initial_geometry()['curvature'])
        if curvature not in ('flat','reference'):raise ValueError('flat or recovered reference only')
        top,minimum=VolumeParent().vacuum_roots()
        starts={'matched':0.,'HM_location':top['s'],'minimum':minimum['s']}
        self.initial_s=starts[preparation];self.preparation=preparation
        self.name=preparation+'_'+curvature
        self.model=PreparedModel(initial_s=self.initial_s,curvature=self.C)
        self.Vstar=minimum['vacuum'];self.minimum_N=-math.log1p(zmax)-1e-8;self.calls=0
        def rhs(N,y):
            self.calls+=1
            if self.calls>30000:raise RuntimeError('preparation allowance exhausted')
            s,v,J=y;dy=self.model.averaged_rhs(0.,[N,s,v]);H=dy[0]
            return v/H,dy[2]/H,math.exp(-N)/H
        initial=[self.initial_s,0.,0.]
        self.past=solve_ivp(rhs,(0,self.minimum_N),initial,method='DOP853',rtol=rtol,atol=rtol/100,dense_output=True)
        self.future=solve_ivp(rhs,(0,self.bounds[1]),initial,method='DOP853',rtol=rtol,atol=rtol/100,dense_output=True)
        if not self.past.success or not self.future.success:raise RuntimeError('preparation domain unresolved')

    def state(self,N):
        N=np.atleast_1d(N).astype(float)
        if min(N)<self.minimum_N-1e-12 or max(N)>self.bounds[1]+1e-12:raise ValueError('outside prepared domain')
        y=np.empty((3,len(N)));negative=N<0
        if np.any(negative):y[:,negative]=self.past.sol(N[negative])
        if np.any(~negative):y[:,~negative]=self.future.sol(N[~negative])
        s,v,J=y
        E=.375*np.exp(-3*N-6*(s-self.initial_s))
        V=7*np.exp(-12*s)-7.5*np.exp(-14*s)+25/32*np.exp(-24*s)
        H=np.sqrt((42*v*v+V+E)/3-self.C*np.exp(-2*N))
        return s,v,J,E,V,H

    def distances(self,theta,z,hel):
        z=np.asarray(z);ob=self.state([theta]);em=self.state(theta-np.log1p(z))
        H=ob[5][0];delta=ob[2][0]-em[2]
        if np.any(delta<=0):raise ValueError('nonpositive distance')
        angle=math.sqrt(self.C)*delta
        if np.any(angle>=math.pi/2):raise ValueError('closed ray branch exceeded')
        DM=H*math.exp(theta)*delta*np.sinc(angle/np.pi)
        return (1+np.asarray(hel))*DM,DM,H/em[5]

    def eos(self,N):
        s,v,J,E,V,H=self.state(N);kin=42*v*v
        rho=kin+V-self.Vstar;p=kin-V+self.Vstar
        valid=abs(rho)>1e-12
        w=np.divide(p,rho,out=np.full_like(rho,np.nan),where=valid)
        correction=np.divide(2*E*v,H*rho,out=np.full_like(rho,np.nan),where=valid)
        return dict(w_s=w,w_s_effective=w-correction,w_dark_energy=(kin-V)/(kin+V),
                    rho_s=rho,p_s=p,Omega_s=rho/(3*H*H),s=s,v=v)

    def eos_rows(self,observer,redshifts=(0.,.25,.5,1.,1.5,2.)):
        z=np.asarray(redshifts);values=self.eos(observer-np.log1p(z))
        return [dict(z=float(zz),**{k:(float(a[i]) if math.isfinite(a[i]) else None) for k,a in values.items()}) for i,zz in enumerate(z)]
