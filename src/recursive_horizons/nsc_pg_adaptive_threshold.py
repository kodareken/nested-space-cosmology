"""Adaptive integration of the same bounded subgap CTP state kernels."""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import quad_vec

from .nsc_pg_retarded_threshold import coherence_coefficients


@dataclass
class AdaptiveThresholdIntegral:
    gram_positive_energy: np.ndarray
    centered_positive_energy: np.ndarray
    analytic_integrals: np.ndarray
    smooth_integral: np.ndarray
    diagnostics: list


def adaptive_threshold_integral(left,right,kappa,node,*,height=None,epsabs=1e-10):
    height=kappa/3 if height is None else height
    if not 0<left<right or not 0<height<kappa/2 or epsabs<=0:
        raise ValueError('existing subgap panel and positive numerical accuracy required')
    diagnostics=[]
    def smooth(t):
        z=complex(left+t*(right-left));data=node(z)
        _,part,_=coherence_coefficients(z,kappa,data['projection'],data['projection_conjugate'])
        return (right-left)*part
    value,error,info=quad_vec(smooth,0.,1.,epsabs=epsabs/4,epsrel=0.,limit=200,full_output=True)
    diagnostics.append({'segment':'smooth','estimated_error':float(error),'evaluations':int(info.neval),
                        'success':bool(info.success),'message':info.message})
    total=np.zeros((3,8,8),complex)
    for name,a,b in (('left',complex(left),left+1j*height),
                     ('upper',left+1j*height,right+1j*height),
                     ('right',right+1j*height,complex(right))):
        def integrand(t):
            z=a+t*(b-a);data=node(z)
            d,_,coherence=coherence_coefficients(z,kappa,data['projection'],data['projection_conjugate'])
            G=data['packet_resolvent'];R=data['reflection'].item()
            return (b-a)*np.array([G,d*G,R*coherence])
        part,error,info=quad_vec(integrand,0.,1.,epsabs=epsabs/4,epsrel=0.,limit=200,full_output=True)
        total+=part
        diagnostics.append({'segment':name,'estimated_error':float(error),'evaluations':int(info.neval),
                            'success':bool(info.success),'message':info.message})
    gram=-1j*(total[0]-total[0].conj().T)
    centered=-1j*(total[1]-total[1].conj().T)+value+total[2]+total[2].conj().T
    return AdaptiveThresholdIntegral(gram,centered,total,value,diagnostics)
