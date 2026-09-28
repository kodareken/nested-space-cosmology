"""Common-Weyl cocycle of the specified static proper-time Dirac modulus.

The ratio discretization preserves continuum conformal covariance exactly at
finite spatial resolution. It is compared to the previous convergent scheme;
published operators and records are unchanged. A compensator lift is a
conditional redundancy construction, not a derived physical dilaton measure.
"""
from dataclasses import replace
from functools import lru_cache
import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh
from scipy.special import exp1
from .nsc_covariant_operator import (SIGMA1,SIGMA2,SIGMA3,ultraviolet_subtraction,
                                     local_trace_coefficient)
from .nsc_regulated import RegulatedOperator,OperatorConventions


def conformal_metric(metric,phi):
    phi=np.broadcast_to(np.asarray(phi,float),(metric.points,))
    if not np.isfinite(phi).all():raise ValueError('finite Weyl field required')
    factor=np.exp(phi)
    return replace(metric,lapse=factor*metric.lapse,radial_scale=factor*metric.radial_scale,
                   sphere_radius=factor*metric.sphere_radius)


def ratio_hamiltonian(metric,kappa=1):
    """H=-i rho2[v partial_x+v'/2]+rho1 m, v=N/q, m=N kappa/r."""
    if isinstance(kappa,bool) or not isinstance(kappa,(int,np.integer)) or kappa==0:
        raise ValueError('nonzero integer angular label required')
    v=metric.lapse/metric.radial_scale;m=metric.lapse*kappa/metric.sphere_radius
    p=metric.momentum_matrix
    return np.kron(SIGMA2,.5*(v[:,None]*p+p*v[None,:]))+np.kron(SIGMA1,np.diag(m))


def ratio_operator(metric,omega,kappa=1):
    if not np.isfinite(omega):raise ValueError('finite frequency required')
    beta=np.kron(SIGMA3,np.eye(metric.points));h=ratio_hamiltonian(metric,kappa)
    w=1/np.sqrt(np.r_[metric.lapse,metric.lapse])
    return w[:,None]*(beta@(omega*np.eye(2*metric.points)+1j*h))*w[None,:]


def frequency_extent(metric,cutoff):
    if not np.isfinite(cutoff) or cutoff<=0:raise ValueError('positive cutoff required')
    c=ratio_operator(metric,0.);a=np.kron(SIGMA3,np.diag(1/metric.lapse))
    cross=float(np.linalg.norm(a@c+c@a,2));nmax=max(metric.lapse)
    return float(.5*nmax**2*(cross+np.sqrt(cross**2+256*cutoff**2/nmax**2)))


def fiber(metric,omega,kappa,cutoff,gradients=True):
    matrix=ratio_operator(metric,omega,kappa)
    e,u=eigh(matrix,driver='evd',check_finite=False)
    if min(abs(e))<1e-10:raise ValueError('zero-mode measure must be specified separately')
    heat=np.exp(-(e/cutoff)**2)
    value=.5*np.sum(exp1((e/cutoff)**2))
    heat_diagonal=np.sum(abs(u)**2*heat[None,:],axis=1).reshape(2,metric.points).sum(axis=0)
    result={'energy':float(value),'heat_diagonal':heat_diagonal}
    if gradients:
        derivative=(u*(-heat/e)[None,:])@u.conj().T
        w=1/np.sqrt(np.r_[metric.lapse,metric.lapse])
        transformed=w[:,None]*derivative*w[None,:]
        blocks=transformed.reshape(2,metric.points,2,metric.points)
        tr=lambda gamma:np.einsum('ba,aibj->ij',gamma,blocks)
        p=metric.momentum_matrix;t1=tr(SIGMA1)
        gv=.5*np.diag(p@t1+t1@p).real
        gm=-np.diag(tr(SIGMA2)).real
        v=metric.lapse/metric.radial_scale;m=metric.lapse*kappa/metric.sphere_radius
        log_gradient=np.array([heat_diagonal+gv*v+gm*m,-gv*v,-gm*m])
        result['log_gradients']=log_gradient
        result['metric_gradients']=log_gradient/np.array([metric.lapse,metric.radial_scale,metric.sphere_radius])
    return result


def cutoff_functional(metric,cutoff=2.,angular_max=8,frequency_points=48,gradients=True):
    if not isinstance(angular_max,int) or angular_max<1 or frequency_points<16:
        raise ValueError('positive angular maximum and resolved frequency quadrature required')
    extent=frequency_extent(metric,cutoff)
    nodes,weights=leggauss(frequency_points)
    energy=0.;gradient=np.zeros((3,metric.points));heat=np.zeros(metric.points)
    for kappa in range(1,angular_max+1):
        for omega,weight in zip((nodes+1)*extent/2,weights*extent/2):
            row=fiber(metric,omega,kappa,cutoff,gradients)
            factor=4*kappa*weight/np.pi
            energy+=factor*row['energy'];heat+=factor*row['heat_diagonal']
            if gradients:gradient+=factor*row['metric_gradients']
    result={'energy':float(energy),'heat_diagonal':heat,'frequency_extent':extent}
    if gradients:
        fields=np.array([metric.lapse,metric.radial_scale,metric.sphere_radius])
        result['metric_gradients']=gradient
        result['log_gradients']=gradient*fields
        result['local_Weyl_Ward_residual']=float(np.max(abs((gradient*fields).sum(axis=0)-heat)))
    return result


def subtracted_functional(metric,cutoff=2.,normalization=1.,angular_max=8,frequency_points=48):
    raw=cutoff_functional(metric,cutoff,angular_max,frequency_points)
    sub=ultraviolet_subtraction(metric,cutoff,normalization)
    return {'energy':raw['energy']-sub['energy'],
            'metric_gradients':raw['metric_gradients']-sub['metric_gradients'],
            'raw':raw,'subtraction':sub}


def operator_variations(metric,omega,kappa,log_direction):
    direction=np.asarray(log_direction,float)
    if direction.shape!=(3,metric.points) or not np.isfinite(direction).all():
        raise ValueError('three finite log-metric directions required')
    a,b,c=direction;p=metric.momentum_matrix
    v=metric.lapse/metric.radial_scale;m=metric.lapse*kappa/metric.sphere_radius
    def hder(order):
        dv=v*(a-b)**order;dm=m*(a-c)**order
        return np.kron(SIGMA2,.5*(dv[:,None]*p+p*dv[None,:]))+np.kron(SIGMA1,np.diag(dm))
    beta=np.kron(SIGMA3,np.eye(metric.points));w=1/np.sqrt(np.r_[metric.lapse,metric.lapse])
    aa=np.diag(np.r_[a,a]);d=ratio_operator(metric,omega,kappa)
    v1=w[:,None]*(1j*beta@hder(1))*w[None,:]
    v2=w[:,None]*(1j*beta@hder(2))*w[None,:]
    first=-.5*(aa@d+d@aa)+v1
    second=.25*(aa@aa@d+d@aa@aa)+.5*aa@d@aa-aa@v1-v1@aa+v2
    return d,first,second


def directional_derivatives(metric,direction,cutoff=2.,angular_max=4,frequency_points=32):
    extent=frequency_extent(metric,cutoff);nodes,weights=leggauss(frequency_points)
    first=second=0.
    for kappa in range(1,angular_max+1):
        for omega,weight in zip((nodes+1)*extent/2,weights*extent/2):
            d,d1,d2=operator_variations(metric,omega,kappa,direction)
            g,h=RegulatedOperator(d,OperatorConventions(cutoff=cutoff)).variation(d1,d2)
            factor=4*kappa*weight/np.pi;first+=factor*g;second+=factor*h
    return {'first':float(first),'second':float(second)}


def heat_path_cocycle(metric,phi,cutoff=2.,angular_max=6,frequency_points=32,path_points=8):
    """Integrate dGamma[g_t]/dt=Tr(phi exp[-D_t²/Lambda²]) independently."""
    phi=np.broadcast_to(np.asarray(phi,float),(metric.points,))
    nodes,weights=leggauss(path_points);value=0.
    for t,weight in zip((nodes+1)/2,weights/2):
        row=cutoff_functional(conformal_metric(metric,t*phi),cutoff,angular_max,frequency_points,False)
        value+=weight*np.dot(phi,row['heat_diagonal'])
    return float(value)


def local_anomaly_cocycle(metric,phi,path_points=12):
    """Continuum four-derivative anomaly integral in the declared scheme."""
    phi=np.broadcast_to(np.asarray(phi,float),(metric.points,))
    nodes,weights=leggauss(path_points);value=0.
    for t,weight in zip((nodes+1)/2,weights/2):
        m=conformal_metric(metric,t*phi)
        volume=4*np.pi*m.spacing*m.lapse*m.radial_scale*m.sphere_radius**2
        value+=weight*np.dot(phi*volume,local_trace_coefficient(m))
    return float(value)


def compensated_lift(metric,phi,cutoff=2.,angular_max=6,frequency_points=32):
    """Conditional field-redundancy lift E_hat[g,phi]=E[g_phi].

    Delta Gamma=E[g_phi]-E[g] is the compensating action. Adding it once
    gives the lifted energy; adding an independent heat trace would duplicate
    a contribution. No invariant finite functional or physical phi measure
    is selected. The fixed-background scale average is not divided by fiat.
    """
    phi=np.broadcast_to(np.asarray(phi,float),(metric.points,))
    dressed=conformal_metric(metric,phi)
    base=cutoff_functional(metric,cutoff,angular_max,frequency_points,False)
    lifted=cutoff_functional(dressed,cutoff,angular_max,frequency_points)
    g=lifted['metric_gradients']*np.exp(phi)[None,:]
    phi_gradient=np.sum(lifted['metric_gradients']*np.array(
        [dressed.lapse,dressed.radial_scale,dressed.sphere_radius]),axis=0)
    ward=np.sum(g*np.array([metric.lapse,metric.radial_scale,metric.sphere_radius]),axis=0)-phi_gradient
    return {'base_energy':base['energy'],'lifted_energy':lifted['energy'],
            'compensating_action':lifted['energy']-base['energy'],
            'metric_gradients':g,'phi_gradient':phi_gradient,'Ward_residual':float(max(abs(ward)))}


@lru_cache(None)
def continuum_identities():
    x=sp.symbols('x',real=True)
    n,q,r,s=(sp.Function(name,positive=True)(x) for name in ('N','q','r','sigma'))
    scale=sp.exp(s);v=n/q;m=n/r
    vv=scale*n/(scale*q);mm=scale*n/(scale*r)
    identities={'principal_ratio':sp.simplify(vv-v),'angular_ratio':sp.simplify(mm-m),
                'radial_connection':sp.simplify(sp.diff(vv,x)/2-sp.diff(v,x)/2)}
    if any(z!=0 for z in identities.values()):raise ArithmeticError('continuum identity failed')
    return {key:str(value) for key,value in identities.items()}
