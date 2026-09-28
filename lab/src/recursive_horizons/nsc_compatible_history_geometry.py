"""Compatible radius directions in the fixed incoming KS/PG chart.

Widths, profiles and amplitudes are explicit caller data. This adapter neither
selects physical initial data nor establishes preservation/preparation of C0.
Each f(z,n) callback must supply the actual nth derivative of ONE smooth real
function. Finite output checks cannot certify that derivative-coherence contract.
"""
from dataclasses import dataclass

import numpy as np

from .nsc_ks_spacetime_variation import chart_coordinates
from .nsc_pg_ks_metric_pullback import reference_chart,ks_to_pg


INCOMING_SLOT_DERIVATIVES = {
    'w': ('r',1,0), 'w_z': ('r',1,1), 'w_zz': ('r',1,2),
    'w_zzz': ('r',1,3), 'U': ('r',3,0), 'U_z': ('r',3,1),
}


def _real(value,label):
    raw=np.asarray(value)
    if np.iscomplexobj(raw): raise ValueError(label+' must be real')
    try: result=np.asarray(value,dtype=float)
    except (TypeError,ValueError,OverflowError) as error:
        raise ValueError(label+' must be finite real numerical data') from error
    if not np.isfinite(result).all(): raise ValueError(label+' must be finite')
    return result


def _scalar(value,label):
    result=_real(value,label)
    if result.ndim!=0: raise ValueError(label+' must be scalar')
    return float(result)


def _return(value): return float(value) if value.ndim==0 else value


def _profile(callback,z,order,label):
    """Allow ordinary scalar analytic callbacks, including math-module users."""
    values=np.empty(z.shape,dtype=float)
    for index in np.ndindex(z.shape):
        try: value=callback(float(z[index]),order)
        except (TypeError,ValueError,OverflowError) as error:
            raise ValueError(f'{label}(z,{order}) could not supply a real derivative') from error
        values[index]=_scalar(value,f'{label}(z,{order})')
    return values


def _plateau(s,inner,outer):
    """Even C-infinity plateau, using a stable logistic of flat exponentials."""
    distance=np.abs(s)
    result=np.zeros(s.shape,dtype=float)
    result[distance<=inner]=1.
    transition=(distance>inner)&(distance<outer)
    # chi=b(1-t)/(b(t)+b(1-t)), b(t)=exp(-1/t). Distances
    # avoid a rounded t=1 at the outer endpoint; logaddexp avoids overflow.
    if np.any(transition):
        toward=distance[transition]-inner
        away=outer-distance[transition]
        with np.errstate(over='ignore',divide='ignore',under='ignore'):
            ratio=(outer-inner)/away-(outer-inner)/toward
            result[transition]=np.exp(-np.logaddexp(0.,ratio))
    return result


@dataclass(frozen=True)
class CompatibleRadiusDirection:
    """Raw radius basis chi(s)*(s*w(z)+s^3*U(z)/6), s=0 at rho1.

    w(z,n) and U(z,n) must be coherent derivatives of their respective smooth
    real functions. Smoothness and derivative consistency are caller obligations;
    this class only validates the finite values it actually requests.
    """
    w: object
    U: object
    inner_radius: float
    outer_radius: float

    def __post_init__(self):
        if not callable(self.w) or not callable(self.U):
            raise ValueError('w and U must be derivative callbacks f(z,order)')
        if isinstance(self.inner_radius,(bool,np.bool_)) or isinstance(self.outer_radius,(bool,np.bool_)):
            raise ValueError('explicit positive inner/outer normal-window radii required')
        inner=_scalar(self.inner_radius,'inner radius');outer=_scalar(self.outer_radius,'outer radius')
        if not 0<inner<outer:
            raise ValueError('normal-window radii must satisfy 0 < inner < outer')
        object.__setattr__(self,'inner_radius',inner)
        object.__setattr__(self,'outer_radius',outer)

    def value(self,s,z):
        """Scalar or broadcast-array evaluation; no elapsed time is assigned."""
        normal,space=np.broadcast_arrays(_real(s,'normal coordinate'),_real(z,'axial coordinate'))
        window=_plateau(normal,self.inner_radius,self.outer_radius)
        result=np.zeros(normal.shape,dtype=float);active=window!=0
        if np.any(active):
            w=_profile(self.w,space[active],0,'w');U=_profile(self.U,space[active],0,'U')
            with np.errstate(over='ignore',invalid='ignore'):
                result[active]=window[active]*(normal[active]*w+normal[active]**3*U/6)
        if not np.isfinite(result).all():raise ValueError('radius direction is not finite')
        return _return(result)

    def incoming_slots(self,z):
        """Six actual derivative values; the Cauchy owner applies factorials."""
        space=_real(z,'axial coordinate')
        result={name:_return(_profile(self.w,space,order,'w'))
                for name,order in (('w',0),('w_z',1),('w_zz',2),('w_zzz',3))}
        result.update({name:_return(_profile(self.U,space,order,'U'))
                       for name,order in (('U',0),('U_z',1))})
        return result


@dataclass(frozen=True)
class CompatibleIncomingMetric:
    """Pure raw-KS radius family pulled through the unchanged reference chart.

    values returns (N_PG,beta_PG,q_PG,r). log_directions returns derivatives
    of (log N_PG,beta_PG,log q_PG,log r) with respect to each amplitude.
    Neither this extension nor its cutoff proves fixed-C0 parent preparation.
    """
    amplitudes: object
    directions: object

    def __post_init__(self):
        amplitudes=_real(self.amplitudes,'amplitudes')
        if amplitudes.ndim!=1:raise ValueError('one-dimensional amplitudes required')
        try: directions=tuple(self.directions)
        except TypeError as error:raise ValueError('a direction sequence is required') from error
        if len(directions)!=len(amplitudes) or any(not isinstance(d,CompatibleRadiusDirection) for d in directions):
            raise ValueError('one compatible radius direction per amplitude required')
        object.__setattr__(self,'amplitudes',tuple(float(v) for v in amplitudes))
        object.__setattr__(self,'directions',directions)

    def _evaluate(self,time,x):
        tau=_scalar(time,'PG time coordinate');rho=_real(x,'PG spatial grid')
        if rho.ndim!=1 or not len(rho):raise ValueError('nonempty one-dimensional PG spatial grid required')
        origin=chart_coordinates(1.)[0]
        values=np.empty((4,len(rho)));basis=np.empty((len(self.directions),len(rho)))
        for j,point in enumerate(rho):
            _,_,axial,_,_=reference_chart(float(point))
            T,S=chart_coordinates(float(point));normal=T-origin;space=tau+S
            for i,direction in enumerate(self.directions):basis[i,j]=direction.value(normal,space)
            with np.errstate(over='ignore',invalid='ignore'):
                radius=np.sqrt(1+point*point)+np.dot(self.amplitudes,basis[:,j])
            if not np.isfinite(radius) or radius<=0:raise ValueError('actual radius must stay finite and positive')
            values[:,j]=ks_to_pg(float(point),1.,0.,axial,float(radius))
        return values,basis

    def values(self,time,x):
        return self._evaluate(time,x)[0]

    def log_directions(self,time,x,metric):
        actual,basis=self._evaluate(time,x)
        supplied=_real(metric,'actual PG metric')
        if supplied.shape!=actual.shape or np.any(supplied[[0,2,3]]<=0):
            raise ValueError('actual PG metric must have shape (4,n) and positive N,q,r')
        if not np.allclose(supplied,actual,rtol=1e-12,atol=1e-12):
            raise ValueError('supplied metric does not match this actual amplitude family')
        result=np.zeros((len(self.directions),4,actual.shape[1]))
        result[:,3,:]=basis/supplied[3]
        return result
