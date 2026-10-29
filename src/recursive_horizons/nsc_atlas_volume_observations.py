"""Conditional SN/BAO inference for the actual evolving homogeneous volume.

Photons and luminosity standards use the Einstein-frame geometric redshift
convention. No baryon, detector-state or shared-wall completion is inferred.
"""
import math,time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp,log_ndtr

from .nsc_atlas_parent_scale_selector import VolumeParent
from .nsc_atlas_inherited_matter import exact_initial_geometry
from .nsc_atlas_parent_supernova import SupernovaShapeLikelihood,flrw_distance
from .nsc_atlas_shared_ruler import SharedRuler


class VolumeBackground:
    bounds=(0.,math.log(20))
    name='volume_12_2'

    def __init__(self,zmax,*,rtol=2e-11):
        self.model=VolumeParent()
        self.C=float(exact_initial_geometry()['curvature'])
        self.Vstar=self.model.vacuum_roots()[1]['vacuum']
        self.minimum_N=-math.log1p(zmax)-1e-8
        self.calls=0
        def rhs(N,y):
            self.calls+=1
            if self.calls>30000:raise RuntimeError('background allowance exceeded')
            s,v,J=y
            dy=self.model.averaged_rhs(0.,[N,s,v]);H=dy[0]
            return [v/H,dy[2]/H,math.exp(-N)/H]
        self.past=solve_ivp(rhs,(0,self.minimum_N),[0,0,0],method='DOP853',rtol=rtol,atol=rtol/100,dense_output=True)
        self.future=solve_ivp(rhs,(0,self.bounds[1]),[0,0,0],method='DOP853',rtol=rtol,atol=rtol/100,dense_output=True)
        if not self.past.success or not self.future.success:raise RuntimeError('homogeneous light-cone domain unresolved')

    def state(self,N):
        N=np.atleast_1d(N).astype(float)
        if np.min(N)<self.minimum_N-1e-12 or np.max(N)>self.bounds[1]+1e-12:raise ValueError('emitter outside computed domain')
        y=np.empty((3,len(N)));negative=N<0
        if np.any(negative):y[:,negative]=self.past.sol(N[negative])
        if np.any(~negative):y[:,~negative]=self.future.sol(N[~negative])
        s,v,J=y
        E=.375*np.exp(-3*N-6*s)
        V=7*np.exp(-12*s)-7.5*np.exp(-14*s)+25/32*np.exp(-24*s)
        H=np.sqrt((42*v*v+V+E)/3-self.C*np.exp(-2*N))
        return s,v,J,E,V,H

    def fractions(self,theta):
        s,v,J,E,V,H=(x[0] for x in self.state([theta]))
        denom=3*H*H
        return dict(Omega_Lambda=self.Vstar/denom,Omega_m=E/denom,
            Omega_k=-self.C*math.exp(-2*theta)/(H*H),
            Omega_volume=(42*v*v+V-self.Vstar)/denom,
            Omega_potential=V/denom,Omega_kinetic=42*v*v/denom)

    def distances(self,theta,z,hel):
        z=np.asarray(z);Ne=theta-np.log1p(z)
        ob=self.state([theta]);em=self.state(Ne)
        H=ob[5][0];delta=ob[2][0]-em[2]
        angle=math.sqrt(self.C)*delta
        if np.any(angle<=0) or np.any(angle>=math.pi/2):raise ValueError('light cone leaves monotone closed branch')
        DM=H*math.exp(theta)*delta*np.sinc(angle/np.pi)
        return (1+np.asarray(hel))*DM,DM,H/em[5]


class FlatBackground:
    bounds=(0.,1.)
    name='flat_LCDM'
    def fractions(self,theta):return dict(Omega_m=theta,Omega_Lambda=1-theta,Omega_k=0.,Omega_volume=0.)
    def distances(self,theta,z,hel):
        DL=flrw_distance(z,hel,matter=theta)
        return DL,DL/(1+np.asarray(hel)),1/np.sqrt(theta*(1+np.asarray(z))**3+1-theta)


def fit(background,data,*,order=48):
    """One shape, offset and shared ruler, with whole-prior quadrature."""
    started=time.process_time()
    z,hel,mag,C,bz,labels,mean,BC=data
    sn=SupernovaShapeLikelihood(mag,C);bao=SharedRuler(mean,BC)
    lo,hi=background.bounds;cache={}
    def value(theta):
        key=float(theta)
        if key not in cache:
            if len(cache)>6000 or time.process_time()-started>180:raise RuntimeError('inference allowance exceeded')
            DL,_,_=background.distances(key,z,hel)
            _,DM,DH=background.distances(key,bz,bz)
            shapes={'DM_over_rs':DM,'DH_over_rs':DH,'DV_over_rs':np.cbrt(bz*DM*DM*DH)}
            shape=np.array([shapes[label][i] for i,label in enumerate(labels)])
            cache[key]=(sn.evaluate(DL)['chi_squared'],bao.profile(shape))
        return cache[key]
    def objective(theta,joint):
        a,b=value(theta);return a+b['chi_squared'] if joint else a
    # Bracket resolved local minima over the whole stated one-dimensional prior.
    grid=np.linspace(lo,hi,129)
    all_minima=[]
    for joint in (False,True):
        vals=[objective(x,joint) for x in grid];found=[]
        for i in range(1,len(grid)-1):
            if vals[i]<=vals[i-1] and vals[i]<=vals[i+1]:
                r=minimize_scalar(lambda x:objective(x,joint),bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':2e-11})
                found.append((float(r.fun),float(r.x)))
        found.extend([(vals[0],lo),(vals[-1],hi)])
        all_minima.append(sorted(found))
    sn_best=all_minima[0][0];joint_best=all_minima[1][0]
    breaks={lo,hi}
    for joint,points in enumerate(all_minima):
        for chi,theta in points:
            if chi>points[0][0]+100 or theta in (lo,hi):continue
            h=min(1e-4,(theta-lo)/3,(hi-theta)/3)
            curvature=(objective(theta+h,joint)-2*chi+objective(theta-h,joint))/h**2
            if curvature<=0:raise ValueError('nonpositive resolved peak curvature')
            width=math.sqrt(2/curvature)
            breaks.update(max(lo,min(hi,theta+n*width)) for n in (-12,-8,-4,0,4,8,12))
    gx,gw=np.polynomial.legendre.leggauss(order);nodes=[];weights=[]
    edges=sorted(breaks)
    for a,b in zip(edges[:-1],edges[1:]):
        nodes.extend((a+b)/2+(b-a)*gx/2);weights.extend((b-a)*gw/2)
    values=[value(x) for x in nodes]
    schi=np.array([a for a,b in values]);bchi=np.array([b['chi_squared'] for a,b in values])
    aa=np.array([b['a'] for a,b in values]);bb=np.array([b['b'] for a,b in values])
    logw=np.log(weights)-.5*(schi-sn_best[0])-.5*bchi-.5*np.log(aa)+log_ndtr(bb/np.sqrt(aa))
    w=np.exp(logw-logsumexp(logw))
    f=[background.fractions(x) for x in nodes]
    moments={}
    for key in f[0]:
        vals=np.array([a[key] for a in f]);average=float(w@vals)
        moments[key]=dict(mean=average,standard_deviation=math.sqrt(max(0.,float(w@(vals-average)**2))))
    theta=joint_best[1];sc,bc=value(theta);fr=background.fractions(theta)
    comparison=moments['Omega_Lambda']
    return dict(model=background.name,order=order,prior_bounds=[lo,hi],
        SN_only_best=dict(theta=sn_best[1],chi_squared=sn_best[0]),
        joint_best=dict(theta=theta,**fr,rd_h_Mpc=bc['rd_h_Mpc'],chi_squared_SN=sc,
            chi_squared_BAO=bc['chi_squared'],chi_squared_total=sc+bc['chi_squared'],fitted_parameters=3),
        posterior=dict(prior='uniform shape and common flat positive BAO beta',moments=moments,nodes=len(nodes)),
        three_Omega_Lambda=dict(best_fit=3*fr['Omega_Lambda'],posterior_mean=3*comparison['mean'],
            posterior_sigma=3*comparison['standard_deviation'],comparison=2.054,comparison_sigma=.022,
            standardized_difference=(3*comparison['mean']-2.054)/math.hypot(3*comparison['standard_deviation'],.022)),
        numerical=dict(likelihood_calls=len(cache),process_CPU_seconds=time.process_time()-started,
                       resolved_minima=len(all_minima[1])-2))
