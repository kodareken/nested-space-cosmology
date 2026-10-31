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

    def __call__(self,t,points):return self.sample(t,points,curvature=True)

    def slice(self,t,points):return self.sample(t,points,curvature=False)
