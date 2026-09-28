"""Analytic energy integrals of endpoint Fourier tails of massive PG modes.

The energy integral uses standard generalized exponential integrals. Its
finite-series/spline error bound is separate from the physical high-energy
mode approximation error; it never licenses replacing a missing state.
"""
from functools import lru_cache
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
import numpy as np
import mpmath as mp


def _normalization_polynomial(branch,lower,terms):
    N=len(branch.reference_series);q=np.zeros(N+1,complex)
    powers=(1/(2*lower))**np.arange(1,N+1);b=branch.reference_beta
    if branch.branch==1:
        q[1:]=1j*(1-b)*branch.reference_series*powers
        ratio=(-1-b)/(1-b)
    else:
        q[1:]=-1j*(1+b)*branch.reference_series.conj()*powers
        ratio=(1-b)/(-1-b)
    d=ratio*np.convolve(q,q.conj());d=np.real_if_close(d,tol=1000).real
    bound=float(np.sum(abs(d)))
    if bound>=1:raise ValueError('current-normalization series is outside its controlled tail range')
    result=np.array([1.]);power=np.array([1.]);coefficient=1.
    for k in range(1,terms+1):
        power=np.convolve(power,d);coefficient*=(-.5-(k-1))/k
        if len(result)<len(power):result=np.pad(result,(0,len(power)-len(result)))
        result[:len(power)]+=coefficient*power
    remainder=bound**(terms+1)/(1-bound)
    return result,remainder,2*(terms+1),1/np.sqrt(1-bound)


@lru_cache(maxsize=2048)
def _moment_matrix(lower,difference,maximum_power):
    metadata={'schema':'NSC-EXPINT-MOMENTS-v1','lower':float(lower).hex(),
              'difference':float(difference).hex(),'maximum_power':maximum_power,
              'decimal_precision':35,'mpmath_version':mp.__version__}
    encoded=json.dumps(metadata,sort_keys=True).encode();key=hashlib.sha256(encoded).hexdigest()
    directory=Path(tempfile.gettempdir())/'nsc-expint-moments-v1';directory.mkdir(exist_ok=True)
    target=directory/(key+'.npz')
    if target.exists():
        with np.load(target,allow_pickle=False) as a:
            if a['metadata_json'].tobytes()!=encoded:raise ValueError('energy-moment cache context differs')
            return a['moments'].copy()
    powers=np.arange(2,maximum_power+1)
    values={}
    with mp.workdps(35):
        for n in range(4,2*maximum_power+1):
            if difference==0:values[n]=lower/(n-1)
            else:values[n]=complex(mp.mpf(lower)*mp.expint(n,-1j*mp.mpf(lower)*mp.mpf(difference)))
    result=np.array([[values[int(p+q)] for q in powers] for p in powers],complex)
    data=io.BytesIO();np.savez(data,moments=result,metadata_json=np.frombuffer(encoded,np.uint8))
    with tempfile.NamedTemporaryFile(dir=directory,delete=False) as f:
        f.write(data.getvalue());temporary=f.name
    os.replace(temporary,target)
    return result


class EndpointTailSeries:
    def __init__(self,projector,lower,*,maximum_power=24,normalization_terms=6):
        if lower<projector.minimum_energy:raise ValueError('tail integral starts below the declared mode range')
        self.lower=float(lower);self.maximum_power=maximum_power;self.columns=[];self.error_coefficients=[]
        positions={'parent':0,'child':2,'child_bulk':4,'exterior_bulk':6}
        # Match MassiveTailPackets source-column order: ext+, int+, global-.
        for branch,column in zip(projector.branches,(1,0,2)):
            normal,normal_remainder,remainder_power,norm_bound=_normalization_polynomial(branch,lower,normalization_terms)
            endpoints={};errors={};unnormalized_bounds={}
            for kind,transform in branch.packets.items():
                start=positions[kind]
                for location,derivatives,sign in ((transform.a,transform.left,-1.),(transform.b,transform.right,1.)):
                    raw=np.zeros((projector.order+transform.degree+2,2),complex)
                    # f=0 exactly at every compact packet endpoint.
                    for j in range(1,transform.degree+1):
                        for n in range(projector.order+1):
                            power=n+j+1
                            raw[power]+=sign*(-1)**j*derivatives[j,n]/((1j)**(j+1)*2**n*lower**power)
                    for power,value in enumerate(raw):
                        unnormalized_bounds[power]=unnormalized_bounds.get(power,0.)+float(np.linalg.norm(value))
                    normalized=np.stack([np.convolve(raw[:,spin],normal) for spin in range(2)],axis=1)
                    target=endpoints.setdefault(float(location),np.zeros((maximum_power+1,8),complex))
                    keep=min(len(normalized),maximum_power+1)
                    target[:keep,start:start+2]+=normalized[:keep]
                    for power in range(keep,len(normalized)):
                        errors[power]=errors.get(power,0.)+float(np.linalg.norm(normalized[power]))
                for n in range(projector.order+1):
                    power=n+transform.degree+1
                    error=norm_bound*float(np.sum(np.linalg.norm(transform.jumps[:,n,:],axis=1)))/(2**n*lower**power)
                    errors[power]=errors.get(power,0.)+error
            for power,value in unnormalized_bounds.items():
                errors[power+remainder_power]=errors.get(power+remainder_power,0.)+normal_remainder*value
            self.columns.append((column,endpoints));self.error_coefficients.append(errors)

    def project(self,energy):
        if energy<self.lower:raise ValueError('energy below tail series range')
        powers=(self.lower/energy)**np.arange(self.maximum_power+1)
        output=np.zeros((8,3),complex)
        for column,endpoints in self.columns:
            for x,coefficients in endpoints.items():
                output[:,column]+=np.exp(1j*energy*x)*(powers@coefficients)
        return output

    def centered_integral(self):
        """Positive-energy centered vacuum-source tail, including dE/(2pi)."""
        return self.integral()

    def integral(self,weights=(-.5,.5,-.5)):
        """Diagonal source weights; use (1,1,1) for the spectral Gram tail."""
        answer=np.zeros((8,8),complex);error_bound=0.
        for (column,endpoints),errors in zip(self.columns,self.error_coefficients):
            occupation=weights[column]
            for x,bx in endpoints.items():
                for y,by in endpoints.items():
                    moments=_moment_matrix(self.lower,float(x-y),self.maximum_power)
                    answer+=occupation*(bx[2:].T@moments@by[2:].conj())/(2*np.pi)
            bounds={p:sum(float(np.linalg.norm(b[p])) for b in endpoints.values()) for p in range(2,self.maximum_power+1)}
            for p,a in bounds.items():
                for q,e in errors.items():error_bound+=abs(occupation)*2*a*e*self.lower/(p+q-1)/(2*np.pi)
            for p,a in errors.items():
                for q,e in errors.items():error_bound+=abs(occupation)*a*e*self.lower/(p+q-1)/(2*np.pi)
        return answer,float(error_bound)
