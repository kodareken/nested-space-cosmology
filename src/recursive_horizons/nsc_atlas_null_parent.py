"""Passive metric/curvature jets from an authenticated parent dense history.

All metric derivatives, including both time derivatives, come from the
same stored metric polynomial. Einstein's equations are not substituted
for curvature. The independently reconstructed scalar stress is returned
for a null-Ricci residual, not used to repair that residual.
"""
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from .nsc_atlas_null_geometry import geometry_from_jets,first_connection
from .nsc_atlas_parent_gh import PAIRS,ParentGH
from .nsc_atlas_s3_frame import arbitrary_points_values
from .nsc_atlas_trapped_parent import GHSurfaceGeometry,_first_geometry


def dense_history_jet(history,t):
    """Analytically differentiate the actual nested DOP853 polynomial twice."""
    t=float(t)
    if not history.bounds or not math.isfinite(t):raise ValueError('finite time and nonempty history required')
    if not history.bounds[0][0]<=t<=history.bounds[-1][1]:raise ValueError('history extrapolation is forbidden')
    ends=np.asarray(history.bounds)[:,1]
    index=min(int(np.searchsorted(ends,t,side='left')),len(ends)-1)
    lower,upper=history.bounds[index];step=upper-lower;x=(t-lower)/step
    value=np.zeros_like(history.old[index]);first=np.zeros_like(value);second=np.zeros_like(value)
    for k,coefficient in enumerate(history.coefficients[index][::-1]):
        factor=x if k%2==0 else 1-x;sign=1. if k%2==0 else -1.
        before=value+coefficient
        second=second*factor+2*sign*first
        first=first*factor+sign*before
        value=before*factor
    return value+history.old[index],first/step,second/step**2


class GHNullHistorySampler:
    """One immutable physical history, with a one-time spectral-jet cache.

    Direct construction is for explicit caller-owned analytic controls.
    Production uses from_files(), which retains all strict producer pins.
    """
    def __init__(self,system,history,*,physical_interval=None,provenance=None):
        if not isinstance(system,ParentGH) or not history.bounds:
            raise ValueError('existing parent system and nonempty dense history required')
        self.system=system;self.history=history;self.r=system.r
        self.data_start=float(history.bounds[0][0])
        self.data_end=float(history.bounds[-1][1] if physical_interval is None else physical_interval[1])
        if physical_interval is not None:self.data_start=float(physical_interval[0])
        if not history.bounds[0][0]<=self.data_start<self.data_end<=history.bounds[-1][1]:
            raise ValueError('physical interval must lie inside recorded polynomials')
        self.provenance=dict(provenance or {'source_binding':'explicit_caller_owned_control'})
        self.calls=0;self.points_evaluated=0;self.elapsed_seconds=0.
        self._time=None;self._coefficients=None

    @classmethod
    def from_files(cls,record_path,history_path=None):
        record_path=Path(record_path)
        record=json.loads(record_path.read_bytes())
        history_path=Path(str(record_path.with_suffix(''))+'-history.npz') if history_path is None else Path(history_path)
        initial=GHSurfaceGeometry.from_history(history_path,record,coordinate_t=0.,retain_history=True)
        endpoint=float(initial.retained_history.metadata['physical_coordinate_endpoint'])
        provenance=dict(initial.provenance,record_path=record_path.as_posix(),
            record_sha256=hashlib.sha256(record_path.read_bytes()).hexdigest(),
            history_path=history_path.as_posix(),
            sampler_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        return cls(initial.system,initial.retained_history,physical_interval=(0.,endpoint),provenance=provenance)

    def _prepare(self,t):
        t=float(t)
        if not math.isfinite(t) or not self.data_start<=t<=self.data_end:
            raise ValueError('sample outside the declared physical history')
        if t!=self._time:
            y,dy,ddy=dense_history_jet(self.history,t)
            C=self.system.unpack(y);C1=self.system.unpack(dy);C2=self.system.unpack(ddy)
            selected=np.r_[np.arange(10),[20,21]]
            first=np.array([self.r*(D@C[selected].T).T for D in self.system.basis.generators])
            mixed=np.array([self.r*(D@C1[:10].T).T for D in self.system.basis.generators])
            second=np.empty((3,3,10,len(self.system.basis.labels)))
            for i,Di in enumerate(self.system.basis.generators):
                for j,Dj in enumerate(self.system.basis.generators):
                    second[i,j]=(self.r**2*(Di@(Dj@C[:10].T))).T
            self._coefficients=(C,C1,C2,first,mixed,second)
            self._time=t
        return self._coefficients

    def sample(self,t,points,*,curvature=True):
        started=time.perf_counter();X=np.asarray(points,dtype=float)
        if X.ndim!=2 or X.shape[1]!=4 or not len(X) or not np.all(np.isfinite(X)):
            raise ValueError('nonempty finite (N,4) embedding points required')
        norms=np.linalg.norm(X,axis=1)
        if min(norms)<=0 or max(abs(norms-1))>1e-6:raise ValueError('point norm drift exceeds the declared sampling domain')
        X=X/norms[:,None]
        C,C1,C2,first,mixed,second=self._prepare(t)
        Y=arbitrary_points_values(X,self.system.basis.labels)
        V=Y@C.T;V1=Y@C1.T;V2=Y@C2[:10].T
        raw_first=np.array([Y@coeff.T for coeff in first]).transpose(1,0,2)
        count=len(X);g=np.empty((count,4,4));v=np.empty_like(g);acc=np.empty_like(g)
        dg=np.empty((count,3,4,4));dv=np.empty_like(dg);d2g=np.empty((count,3,3,4,4))
        if curvature:
            raw_mixed=np.array([Y@coeff.T for coeff in mixed]).transpose(1,0,2)
            raw_second=np.empty((count,3,3,10))
            for i in range(3):
                for j in range(3):raw_second[:,i,j]=Y@second[i,j].T
        for k,(a,b) in enumerate(PAIRS):
            g[:,a,b]=g[:,b,a]=V[:,k]
            v[:,a,b]=v[:,b,a]=V1[:,k]
            acc[:,a,b]=acc[:,b,a]=V2[:,k]
            dg[:,:,a,b]=dg[:,:,b,a]=raw_first[:,:,k]
            if curvature:
                dv[:,:,a,b]=dv[:,:,b,a]=raw_mixed[:,:,k]
                d2g[:,:,:,a,b]=d2g[:,:,:,b,a]=raw_second[:,:,:,k]
        data=_first_geometry(g,v,dg,self.r)
        data.update(g=g,v=v,dg=dg,q=V[:,20:22],dq=raw_first[:,:,10:12],
            W=self.system.potential(V[:,20:22]),metric_velocity_consistency=V1[:,:10]-V[:,10:20],
            field_velocity_consistency=V1[:,20:22]-V[:,22:24],point_norm_error=norms-1)
        data['Pi']=(V1[:,20:22]-np.einsum('ni,nia->na',data['shift'],data['dq']))/data['lapse'][:,None]
        if curvature:
            geometric=geometry_from_jets(g,v,dg,d2g,dv,acc,frame_scale=self.r)
            data.update(geometric)
            data['R_cov']=geometric['Riemann']
            derivatives=np.concatenate([V1[:,None,20:22],data['dq']],axis=1)
            matter=self.system.parent.kappa*np.einsum('naA,nbA->nab',derivatives,derivatives)
            source=matter+self.system.parent.kappa*data['W'][:,None,None]*g
            data['matter_Ricci']=matter
            data['Einstein_Ricci_source']=source
            data['Ricci_source_residual']=geometric['Ricci']-source
        else:data.update(first_connection(g,v,dg,frame_scale=self.r))
        self.calls+=1;self.points_evaluated+=len(X);self.elapsed_seconds+=time.perf_counter()-started
        return data

    def event_slice(self,times,points,*,chunk=32):
        """First metric jets at paired events (times[i], points[i]).

        Only the packed ten-component metric prefix and its actual dense
        polynomial derivative are evaluated. Sparse harmonic derivatives
        act on the point-value rows, avoiding a full field-state unpack or
        second-jet preparation for every distinct event time. Memory is
        bounded by ``chunk`` point times. No history extrapolation is used.
        """
        started=time.perf_counter();raw_t=np.asarray(times);raw_X=np.asarray(points)
        if raw_t.dtype.kind not in 'iuf' or raw_t.ndim!=1 or not len(raw_t) or not np.all(np.isfinite(raw_t)):
            raise ValueError('nonempty finite one-dimensional event times required')
        if raw_X.dtype.kind not in 'iuf' or raw_X.shape!=(len(raw_t),4) or not np.all(np.isfinite(raw_X)):
            raise ValueError('event points need finite shape (number_of_times,4)')
        if isinstance(chunk,(bool,np.bool_)) or not isinstance(chunk,(int,np.integer)) or chunk<1:
            raise ValueError('positive integer event chunk required')
        times=np.asarray(raw_t,dtype=float);X=np.asarray(raw_X,dtype=float)
        if np.any(times<self.data_start) or np.any(times>self.data_end):
            raise ValueError('event outside the declared physical history')
        norms=np.linalg.norm(X,axis=1)
        if min(norms)<=0 or max(abs(norms-1))>1e-6:
            raise ValueError('point norm drift exceeds the declared sampling domain')
        X=X/norms[:,None]
        bounds=np.asarray(self.history.bounds);ends=bounds[:,1]
        metric_size=int(self.system.offsets[10]);parts={}
        for lo in range(0,len(times),chunk):
            hi=min(lo+chunk,len(times));tc=times[lo:hi];count=hi-lo
            interval=np.minimum(np.searchsorted(ends,tc,side='left'),len(ends)-1)
            start=bounds[interval,0];step=bounds[interval,1]-start;x=(tc-start)/step
            value=np.zeros((count,metric_size));derivative=np.zeros_like(value)
            for order in range(7):
                coefficient=np.stack([self.history.coefficients[j][6-order,:metric_size] for j in interval])
                before=value+coefficient
                factor=(x if order%2==0 else 1-x)[:,None];sign=1. if order%2==0 else -1.
                derivative=derivative*factor+sign*before
                value=before*factor
            value+=np.stack([self.history.old[j][:metric_size] for j in interval])
            derivative/=step[:,None]
            Y=arbitrary_points_values(X[lo:hi],self.system.basis.labels)
            Y_i=[self.r*(D.T@Y.T).T for D in self.system.basis.generators]
            g=np.empty((count,4,4));v=np.empty_like(g);dg=np.empty((count,3,4,4))
            for j,(a,b) in enumerate(PAIRS):
                lower,upper=self.system.offsets[j:j+2];idx=self.system.indices[j]
                coefficients=value[:,lower:upper];rates=derivative[:,lower:upper]
                g[:,a,b]=g[:,b,a]=np.einsum('ni,ni->n',Y[:,idx],coefficients)
                v[:,a,b]=v[:,b,a]=np.einsum('ni,ni->n',Y[:,idx],rates)
                for i in range(3):
                    dg[:,i,a,b]=dg[:,i,b,a]=np.einsum('ni,ni->n',Y_i[i][:,idx],coefficients)
            data=_first_geometry(g,v,dg,self.r)
            data.update(first_connection(g,v,dg,frame_scale=self.r))
            data.update(g=g,v=v,dg=dg,point_norm_error=norms[lo:hi]-1)
            for key,array in data.items():parts.setdefault(key,[]).append(array)
        result={key:np.concatenate(arrays,axis=0) for key,arrays in parts.items()}
        self.calls+=1;self.points_evaluated+=len(X);self.elapsed_seconds+=time.perf_counter()-started
        return result

    def __call__(self,t,points):return self.sample(t,points,curvature=True)

    def slice(self,t,points):return self.sample(t,points,curvature=False)
