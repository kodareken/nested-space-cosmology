"""Resolved DOP853 background history for passive field transport.

The archive contains the actual dense polynomials, not interpolated sparse
snapshots. It never advances, restarts or changes the physical solver.
"""
import json
import math
from pathlib import Path

import numpy as np


class DenseGHHistory:
    def __init__(self):
        self.bounds=[];self.old=[];self.coefficients=[]
        self.metadata={}

    def record(self,start,end,dense):
        start=float(start);end=float(end)
        if not math.isfinite(start+end) or end<=start:
            raise ValueError('positive finite history interval required')
        if self.bounds and start!=self.bounds[-1][1]:
            raise ValueError('background history must be contiguous')
        old=np.asarray(dense.y_old,dtype=float);F=np.asarray(dense.F,dtype=float)
        if old.ndim!=1 or F.shape!=(7,len(old)) or not np.all(np.isfinite(old)) or not np.all(np.isfinite(F)):
            raise ValueError('finite DOP853 dense polynomial required')
        if float(dense.t_old)!=start or float(dense.t)!=end:
            raise ValueError('dense output does not match the recorded interval')
        if self.old and old.shape!=self.old[0].shape:
            raise ValueError('background state size changed')
        self.bounds.append((start,end));self.old.append(old.copy());self.coefficients.append(F.copy())

    def evaluate(self,t):
        if not self.bounds:raise ValueError('history is empty')
        t=float(t)
        if not math.isfinite(t):raise ValueError('finite history time required')
        lower=self.bounds[0][0];upper=self.bounds[-1][1]
        if t<lower or t>upper:raise ValueError('history extrapolation is forbidden')
        ends=np.asarray(self.bounds)[:,1]
        index=min(int(np.searchsorted(ends,t,side='left')),len(ends)-1)
        start,end=self.bounds[index];x=(t-start)/(end-start)
        y=np.zeros_like(self.old[index])
        for k,coefficient in enumerate(self.coefficients[index][::-1]):
            y+=coefficient;y*=x if k%2==0 else 1-x
        return y+self.old[index]

    def save(self,path,metadata):
        path=Path(path)
        if path.exists():raise ValueError('fresh history path required')
        if not self.bounds:raise ValueError('cannot save an empty history')
        self.metadata=dict(metadata)
        np.savez_compressed(path,bounds=np.asarray(self.bounds),old=np.asarray(self.old),
            coefficients=np.asarray(self.coefficients),metadata=json.dumps(self.metadata,sort_keys=True))

    @classmethod
    def load(cls,path):
        with np.load(path,allow_pickle=False) as data:
            bounds=data['bounds'];old=data['old'];coefficients=data['coefficients']
            metadata=json.loads(str(data['metadata']))
        if bounds.ndim!=2 or bounds.shape[1]!=2 or old.ndim!=2 or len(old)!=len(bounds):
            raise ValueError('invalid history shapes')
        if coefficients.shape!=(len(bounds),7,old.shape[1]) or not all(np.all(np.isfinite(a)) for a in (bounds,old,coefficients)):
            raise ValueError('invalid history polynomials')
        if not len(bounds) or np.any(bounds[:,1]<=bounds[:,0]) or np.any(bounds[1:,0]!=bounds[:-1,1]):
            raise ValueError('invalid or noncontiguous history intervals')
        obj=cls();obj.bounds=list(map(tuple,bounds));obj.old=list(old)
        obj.coefficients=list(coefficients);obj.metadata=metadata
        return obj
