"""Finite candidate preparation from the actual wall-clock energy form.

This defines an explicit state-selection axiom, not an inherited root
vacuum or a Hadamard/ultraviolet result. Controls and admission live in
docs/nsc-atlas-dirac-probe.md. Propagation reuses the existing field kernel.
"""
from __future__ import annotations

import math
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from .nsc_influence import ground_covariance
from .nsc_atlas_dirac_mean import MeanModeGrid,RadarPointEvaluator,SPIN_ROTATION
from .nsc_atlas_dirac_scattering import smooth_bump
from .nsc_atlas_dirac_stress import stress_matrices


def hermite_null_modes(clock,duration,count):
    """Smooth flux-orthonormal trial modes and their analytic derivative.

    The logit coordinate maps the open null interval onto the real line.
    Values outside the interval are the smooth zero extension.
    """
    if (isinstance(count,bool) or not isinstance(count,int) or count<1
            or isinstance(duration,bool) or not np.isfinite(duration) or duration<=0):
        raise ValueError('positive duration and positive integer mode count required')
    s=np.asarray(clock,dtype=float)
    if not np.isfinite(s).all():
        raise ValueError('finite null labels required')
    values=np.zeros(s.shape+(count,))
    derivatives=np.zeros_like(values)
    active=(s>0)&(s<duration)
    if not np.any(active):
        return values,derivatives
    t=s[active]/duration
    y=np.log(t)-np.log1p(-t)
    J=1/(duration*t*(1-t))
    log_J_prime=(2*t-1)/(duration*t*(1-t))
    h=np.empty(y.shape+(count,))
    h[...,0]=np.exp(-y*y/2)/math.pi**.25
    if count>1:
        h[...,1]=math.sqrt(2)*y*h[...,0]
    for n in range(1,count-1):
        h[...,n+1]=math.sqrt(2/(n+1))*y*h[...,n]-math.sqrt(n/(n+1))*h[...,n-1]
    dh=-y[...,None]*h
    for n in range(1,count):
        dh[...,n]+=math.sqrt(2*n)*h[...,n-1]
    values[active]=np.sqrt(J)[...,None]*h
    derivatives[active]=np.sqrt(J)[...,None]*(J[...,None]*dh+log_J_prime[...,None]*h/2)
    return values,derivatives


def clock_energy_matrix(fields,stress):
    """Full canonical energy; the angular multiplet is already in stress.

    Coordinates use a future timelike clock vector and the same ADM shift
    sign as the existing canonical metric vertices. Radar beta is zero.
    """
    rho=np.asarray(stress['energy_density'],dtype=complex)
    flux=np.asarray(stress['outward_energy_flux'],dtype=complex)
    if rho.ndim<2 or rho.shape[-1]!=rho.shape[-2] or flux.shape!=rho.shape:
        raise ValueError('matching square stress matrices required')
    shape=rho.shape[:-2]
    N,a,r,beta=[np.broadcast_to(np.asarray(v,dtype=float),shape) for v in
        (fields['lapse'],fields['radial_length'],fields['radius'],fields.get('shift',0.))]
    if (not all(np.isfinite(v).all() for v in (N,a,r,beta))
            or np.any(N<=0) or np.any(a<=0) or np.any(r<=0)
            or np.any(np.abs(a*beta)>=N)):
        raise ValueError('positive metric and a future timelike clock required')
    result=4*math.pi*(a*r*r)[...,None,None]*(N[...,None,None]*rho
        -(beta*a)[...,None,None]*flux)
    if not np.isfinite(result).all():
        raise ValueError('finite energy matrices required')
    return result


def averaged_clock_energy(energy,weights):
    """Positive weighted integration of the supplied physical energy form."""
    energy=np.asarray(energy,dtype=complex)
    weights=np.asarray(weights,dtype=float)
    if (energy.ndim<2 or energy.shape[-1]!=energy.shape[-2]
            or energy.shape[:-2]!=weights.shape or not np.isfinite(energy).all()
            or not np.isfinite(weights).all() or np.any(weights<0) or not np.any(weights>0)):
        raise ValueError('finite square energy history and matching nonnegative weights required')
    scale=max(1.,float(np.max(np.abs(energy))))
    if not np.allclose(energy,energy.conj().swapaxes(-1,-2),atol=2e-12*scale,rtol=0):
        raise ValueError('Hermitian physical energy form required')
    return np.einsum('p,pij->ij',weights.ravel(),energy.reshape((-1,)+energy.shape[-2:]))


def minimum_clock_covariance(generator,*,operator_uncertainty=0.):
    """Finite energy minimizer with an explicit, scale-invariant gap check.

    operator_uncertainty is a supplied norm uncertainty/diagnostic in the
    generator's units. A returned projector concerns this finite matrix;
    it is not a continuum state theorem. Exact zero occupations remain
    unspecified when numerical separation is unavailable.
    """
    a=np.asarray(generator,dtype=complex)
    if (a.ndim!=2 or a.shape[0]!=a.shape[1] or not len(a) or not np.isfinite(a).all()
            or not np.isfinite(operator_uncertainty) or operator_uncertainty<0):
        raise ValueError('finite square generator and nonnegative uncertainty required')
    scale=float(np.linalg.norm(a,ord=2))
    if not np.allclose(a,a.conj().T,atol=2e-12*scale,rtol=0):
        raise ValueError('Hermitian generator required')
    eigenvalues=np.linalg.eigvalsh(a)
    gap=float(np.min(np.abs(eigenvalues)))
    threshold=max(float(operator_uncertainty),1e-10*scale)
    if scale==0 or gap<=threshold:
        return dict(status='unresolved_spectral_gap',covariance=None,
            eigenvalues=eigenvalues,gap=gap,threshold=threshold)
    _,_,covariance=ground_covariance(a/scale)
    return dict(status='finite_clock_energy_minimizer',covariance=covariance,
        eigenvalues=eigenvalues,gap=gap,threshold=threshold,
        minimum_energy=float(np.trace(covariance@a).real))


def propagate_null_modes(potential,du,dv,*,potential_u,potential_v,
                         incoming_plus,incoming_minus,plus_u,minus_v,workers=8):
    """Propagate independent columns using the existing causal Numba kernel.

    The kernel releases the GIL. Its ordered cell sweep is unchanged;
    only independent modes are threaded. The first mode warms the kernel.
    """
    if isinstance(workers,bool) or not isinstance(workers,int) or not 1<=workers<=8:
        raise ValueError('one to eight mode workers required')
    data=[np.asarray(x,dtype=complex) for x in
          (incoming_plus,incoming_minus,plus_u,minus_v)]
    shape=np.shape(potential)
    if (len(shape)!=2 or any(x.ndim!=2 for x in data)
            or not data[0].shape[1] or any(x.shape[1]!=data[0].shape[1] for x in data)
            or any(x.shape[0]!=n for x,n in zip(data,(shape[0],shape[1],shape[0],shape[1])))):
        raise ValueError('matching (null samples, modes) incoming data required')
    def solve(column):
        return MeanModeGrid(potential,du,dv,potential_u=potential_u,potential_v=potential_v,
            incoming_plus=data[0][:,column],incoming_minus=data[1][:,column],
            plus_u=data[2][:,column],minus_v=data[3][:,column])
    first=solve(0)
    with ThreadPoolExecutor(max_workers=min(workers,data[0].shape[1])) as executor:
        return [first,*executor.map(solve,range(1,data[0].shape[1]))]


class ClockEnergyBasis:
    """Actual-geometry finite Galerkin state, with a declared physical window.

    Raw modes on both null legs are propagated once. Subspace selections
    use their own sampled flux Gram, so comparing trial counts does not
    introduce a different physical envelope or a second field solver.
    No angular completion or Hadamard property follows from this class.
    """

    def __init__(self,chart,cells,mode_count,*,preparation_duration,
                 region='nested',angular_k=1,workers=8):
        if isinstance(cells,bool) or not isinstance(cells,int) or cells<4:
            raise ValueError('at least four cells per null direction required')
        if (isinstance(preparation_duration,bool) or not np.isfinite(preparation_duration)
                or not 0<preparation_duration<=chart.duration):
            raise ValueError('positive physical preparation window inside chart required')
        self.chart,self.region,self.angular_k=chart,region,angular_k
        self.preparation_duration=float(preparation_duration)
        self.points=RadarPointEvaluator(chart,region=region,angular_k=angular_k)
        self.delta=chart.duration/cells
        labels=(np.arange(cells)+.5)*self.delta
        f,df=hermite_null_modes(labels,chart.duration,mode_count)
        self.mode_count=mode_count
        zero=np.zeros_like(f)
        plus,minus=np.c_[f,zero],np.c_[zero,f]
        self.incoming=np.sqrt(self.delta)*np.r_[plus,minus]
        self.gram=self.incoming.conj().T@self.incoming
        U,V=np.meshgrid(labels,labels,indexing='ij')
        fields=self.points(U,V)
        self.grids=propagate_null_modes(fields['potential'],self.delta,self.delta,
            potential_u=(fields['potential_t']-fields['potential_x'])/2,
            potential_v=(fields['potential_t']+fields['potential_x'])/2,
            incoming_plus=plus,incoming_minus=minus,
            plus_u=np.c_[df,zero],minus_v=np.c_[zero,df],workers=workers)
        outgoing=np.stack([grid.outgoing for grid in self.grids],axis=-1)
        self.flux_gram_defect=float(np.linalg.norm(outgoing.conj().T@outgoing-self.gram,ord=2))
        eta=np.stack([grid.jets[0] for grid in self.grids],axis=-1)
        eta_x=np.stack([grid.jets[2]-grid.jets[1] for grid in self.grids],axis=-1)
        stress=self._stress(eta,eta_x,fields)
        S=self.preparation_duration
        weights=self.delta**2*smooth_bump(U,S/2,S/2)*smooth_bump(V,S/2,S/2)/S
        self.raw_generator=averaged_clock_energy(clock_energy_matrix(fields,stress),weights)

    def _stress(self,eta,eta_x,fields):
        xi=np.einsum('ij,...jm->...im',SPIN_ROTATION,eta)
        dx=np.einsum('ij,...jm->...im',SPIN_ROTATION,eta_x)
        return stress_matrices(xi,dx,radial_length=fields['lapse'],
            radius=fields['radius'],angular_k=self.angular_k)

    def select(self,count,*,operator_uncertainty=0.):
        """Whiten one trial subspace, then apply the same physical state rule."""
        if isinstance(count,bool) or not isinstance(count,int) or not 1<=count<=self.mode_count:
            raise ValueError('trial count must lie inside the propagated basis')
        indices=np.r_[np.arange(count),self.mode_count+np.arange(count)]
        gram=self.gram[np.ix_(indices,indices)]
        g,u=np.linalg.eigh(gram)
        if g[0]<=0 or g[-1]/g[0]>1e8:
            raise ValueError('sampled null-flux basis is numerically unresolved')
        transform=(u/np.sqrt(g))@u.conj().T
        generator=transform.conj().T@self.raw_generator[np.ix_(indices,indices)]@transform
        result=minimum_clock_covariance(generator,operator_uncertainty=operator_uncertainty)
        return dict(result,indices=indices,transform=transform,generator=generator,
            gram_condition=float(g[-1]/g[0]),
            gram_discrepancy=float(np.linalg.norm(gram-np.eye(2*count),ord=2)))

    def pressure_difference(self,U,V,selection):
        """Point stress diagnostic for a resolved finite covariance; no UV claim."""
        if np.ndim(U) or np.ndim(V):
            raise ValueError('one bounded detector event required')
        if selection['covariance'] is None:
            raise ValueError('unresolved state cannot supply a pressure expectation')
        fields=self.points(U,V)
        eta=np.stack([grid.evaluate(U,V) for grid in self.grids],axis=-1)
        dx=np.stack([grid.evaluate(U,V,derivative_v=1)-grid.evaluate(U,V,derivative_u=1)
                     for grid in self.grids],axis=-1)
        stress=self._stress(eta,dx,fields)
        vertex=stress['radial_pressure']-stress['angular_pressure']
        indices,transform=selection['indices'],selection['transform']
        vertex=transform.conj().T@vertex[np.ix_(indices,indices)]@transform
        return float(np.trace(selection['covariance']@vertex).real)
