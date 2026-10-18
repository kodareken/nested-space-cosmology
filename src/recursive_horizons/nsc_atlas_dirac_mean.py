"""Bounded finite-mode Dirac mean stress on the actual wall-clock diamond.

One declared exterior occupation difference is propagated on the quiet
background. The canonical mode and its analytic characteristic derivatives
feed the existing stress owner. No waveform samples are differentiated.
This supplies a finite state DIFFERENCE, not a baseline vacuum or anomaly.
See the existing Dirac owner and compute admission for domains and budgets.
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit
from scipy.integrate import quad

from .nsc_atlas_dirac_scattering import (
    _cayley_coefficients, _inputs, _real_array, smooth_bump)
from .nsc_atlas_dirac_stress import SIGMA1, stress_matrices
from .nsc_atlas_dirac_gravity import null_stress_from_orthonormal

SPIN_ROTATION = (np.eye(2)+1j*SIGMA1)/math.sqrt(2)


@njit(cache=True,nogil=True)
def _mode_kernel(w,wu,wv,du,dv,plus,minus,plus_u,minus_v):
    nu,nv = w.shape
    pu,mv = math.sqrt(du),math.sqrt(dv)
    P,M = plus.copy()*pu,minus.copy()*mv
    z,chi = plus_u.copy(),minus_v.copy()
    values = np.empty((nu,nv,2),dtype=np.complex128)
    derivative_u = np.empty_like(values)
    derivative_v = np.empty_like(values)
    derivative_uv = np.empty_like(values)
    for i in range(nu):
        for j in range(nv):
            W = w[i,j]
            c,sine = _cayley_coefficients(W*math.sqrt(du*dv)/4)
            pnew,mnew = c*P[i]+sine*M[j],sine*P[i]+c*M[j]
            p,m = (P[i]+pnew)/(2*pu),(M[j]+mnew)/(2*mv)
            zuv = -.5j*wu[i,j]*m-W*W*p/4
            chiu = -.5j*wv[i,j]*p-W*W*m/4
            znew,chinew = z[i]+dv*zuv,chi[j]+du*chiu
            values[i,j,0],values[i,j,1] = p,m
            derivative_u[i,j,0],derivative_u[i,j,1] = (z[i]+znew)/2,-.5j*W*p
            derivative_v[i,j,0],derivative_v[i,j,1] = -.5j*W*m,(chi[j]+chinew)/2
            derivative_uv[i,j,0],derivative_uv[i,j,1] = zuv,chiu
            P[i],M[j],z[i],chi[j] = pnew,mnew,znew,chinew
    return values,derivative_u,derivative_v,derivative_uv,np.concatenate((P,M))


def _hermite_weights(t,width,order):
    if order == 0:
        return np.stack((2*t**3-3*t*t+1,-2*t**3+3*t*t,
                         width*(t**3-2*t*t+t),width*(t**3-t*t)),axis=-1)
    if order == 1:
        return np.stack(((6*t*t-6*t)/width,(-6*t*t+6*t)/width,
                         3*t*t-4*t+1,3*t*t-2*t),axis=-1)
    if order == 2:
        return np.stack(((12*t-6)/width**2,(-12*t+6)/width**2,
                         (6*t-4)/width,(6*t-2)/width),axis=-1)
    raise ValueError('only derivative orders 0,1,2 are supplied')


class MeanModeGrid:
    """One propagated canonical chiral mode, bounded to its cell-center box."""

    def __init__(self,potential,du,dv,*,potential_u,potential_v,
                 incoming_plus,incoming_minus,plus_u,minus_v):
        w,du,dv,_ = _inputs(potential,du,dv,None,2)
        wu,wv = [_real_array(x,'potential derivative',2)
                 for x in (potential_u,potential_v)]
        if wu.shape != w.shape or wv.shape != w.shape or min(w.shape) < 2:
            raise ValueError('matching potential derivatives and at least two cells required')
        incoming = [np.ascontiguousarray(x,dtype=complex) for x in
                    (incoming_plus,incoming_minus,plus_u,minus_v)]
        for x,size in zip(incoming,(w.shape[0],w.shape[1],w.shape[0],w.shape[1])):
            if x.shape != (size,) or not np.isfinite(x).all():
                raise ValueError('finite characteristic data and analytic derivatives required')
        self.du,self.dv = du,dv
        self.shape = w.shape
        *self.jets,self.outgoing = _mode_kernel(w,wu,wv,du,dv,*incoming)
        self.incoming_flux = float(du*np.vdot(incoming[0],incoming[0]).real+
                                   dv*np.vdot(incoming[1],incoming[1]).real)
        self.flux_error = abs(float(np.vdot(self.outgoing,self.outgoing).real)-
                              self.incoming_flux)

    def evaluate(self,U,V,*,derivative_u=0,derivative_v=0):
        U,V = np.broadcast_arrays(np.asarray(U,dtype=float),np.asarray(V,dtype=float))
        limits = ((self.du/2,(self.shape[0]-.5)*self.du),
                  (self.dv/2,(self.shape[1]-.5)*self.dv))
        if not np.isfinite(U).all() or not np.isfinite(V).all():
            raise ValueError('finite null labels required')
        for x,(lo,hi) in zip((U,V),limits):
            if np.any((x<lo)|(x>hi)):
                raise ValueError('query exits the bounded Dirac cell-center representation')
        i = np.minimum(np.floor(U/self.du-.5).astype(int),self.shape[0]-2)
        j = np.minimum(np.floor(V/self.dv-.5).astype(int),self.shape[1]-2)
        t,s = U/self.du-.5-i,V/self.dv-.5-j
        bu = _hermite_weights(t,self.du,derivative_u)
        bv = _hermite_weights(s,self.dv,derivative_v)
        result = np.zeros(U.shape+(2,),dtype=complex)
        for a in range(2):
            for b in range(2):
                for jet,(left,right) in zip(self.jets,((a,b),(a+2,b),(a,b+2),(a+2,b+2))):
                    result += (bu[...,left]*bv[...,right])[...,None]*jet[i+a,j+b]
        return result


class ExteriorPacket:
    """The previously declared smooth incoming mode, with continuum norm one."""

    def __init__(self,duration):
        if isinstance(duration,bool) or not np.isfinite(duration) or duration<=0:
            raise ValueError('positive finite duration required')
        self.duration = float(duration)
        self.start,self.end = .3*duration,.9*duration
        norm,error = quad(lambda z:math.exp(2-2/(1-z*z)),-1,1,
                          epsabs=2e-13,epsrel=2e-13)
        self.normalization = 1/math.sqrt(.3*duration*norm)
        self.norm_quadrature_error = error/norm

    def data(self,V):
        V = np.asarray(V,dtype=float)
        if not np.isfinite(V).all():
            raise ValueError('finite exterior null labels required')
        h,S = .3*self.duration,self.duration
        z = (V-.6*S)/h
        g = self.normalization*smooth_bump(V,.6*S,h)*np.exp(-12j*V/S)
        derivative = np.zeros_like(g)
        mask = np.abs(z)<1
        derivative[mask] = g[mask]*(-2*z[mask]/(h*(1-z[mask]**2)**2)-12j/S)
        return g,derivative


class RadarPointEvaluator:
    """Reuse exact wall-clock inversions at recurring quadrature nodes."""

    def __init__(self,chart,*,region='nested',angular_k=1):
        self.chart,self.region,self.angular_k = chart,region,angular_k
        self._radii = {}

    def __call__(self,U,V,*,region=None):
        U,V = np.broadcast_arrays(np.asarray(U,dtype=float),np.asarray(V,dtype=float))
        radii=self.wall_radii(np.r_[U.ravel(),V.ravel()])
        ru,rv = radii[:U.size].reshape(U.shape),radii[U.size:].reshape(V.shape)
        return self.chart.from_wall_radii(ru,rv,region=self.region if region is None else region,
                                          angular_k=self.angular_k)

    def wall_radii(self,clock):
        """Share exact wall-clock inversions with the metric quadratures."""
        clock=np.asarray(clock,dtype=float)
        unique,inverse = np.unique(clock.ravel(),return_inverse=True)
        missing = np.array([x for x in unique if float(x) not in self._radii])
        if len(missing):
            solved = self.chart.radius_at_clock(missing)
            self._radii.update((float(x),float(r)) for x,r in zip(missing,solved))
        radii = np.array([self._radii[float(x)] for x in unique])[inverse]
        return radii.reshape(clock.shape)


class PreparedDiracSource:
    """Mean stress per unit occupation for one complete angular multiplet.

    The selected source is the EMPTY-to-OCCUPIED difference of the same
    exterior mode; other inherited parent correlations cancel. All values
    are in hbar=c=ell=1. A consistent caller restores the gravitational
    factor 8*pi*(l_P/ell)^2 explicitly.
    """

    def __init__(self,chart,cells,*,region='nested',angular_k=1):
        if isinstance(cells,bool) or not isinstance(cells,int) or cells<4:
            raise ValueError('at least four cells per null direction required')
        self.chart,self.region,self.angular_k = chart,region,angular_k
        self.points = RadarPointEvaluator(chart,region=region,angular_k=angular_k)
        self.packet = ExteriorPacket(chart.duration)
        self.start = self.packet.start
        labels = (np.arange(cells)+.5)*chart.duration/cells
        self.breakpoints = labels
        U,V = np.meshgrid(labels,labels,indexing='ij')
        fields = self.points(U,V)
        plus = np.zeros(cells,dtype=complex)
        minus,minus_v = self.packet.data(labels)
        self.grid = MeanModeGrid(fields['potential'],chart.duration/cells,chart.duration/cells,
            potential_u=(fields['potential_t']-fields['potential_x'])/2,
            potential_v=(fields['potential_t']+fields['potential_x'])/2,
            incoming_plus=plus,incoming_minus=minus,plus_u=plus,minus_v=minus_v)

    def evaluate(self,U,V,*,region=None):
        U,V = np.broadcast_arrays(np.asarray(U,dtype=float),np.asarray(V,dtype=float))
        fields = self.points(U,V,region=region)
        active = V>self.packet.start
        eta = np.zeros(U.shape+(2,),dtype=complex)
        eta_x = np.zeros_like(eta)
        if np.any(active):
            u,v = U[active],V[active]
            eta[active] = self.grid.evaluate(u,v)
            eta_x[active] = (self.grid.evaluate(u,v,derivative_v=1)
                            -self.grid.evaluate(u,v,derivative_u=1))
        xi = np.einsum('ij,...j->...i',SPIN_ROTATION,eta)[...,None]
        xi_x = np.einsum('ij,...j->...i',SPIN_ROTATION,eta_x)[...,None]
        matrices = stress_matrices(xi,xi_x,radial_length=fields['lapse'],
                                   radius=fields['radius'],angular_k=self.angular_k)
        mean = {k:v[...,0,0].real for k,v in matrices.items()}
        null = null_stress_from_orthonormal(C=fields['conformal_factor'],
            rho=mean['energy_density'],radial_pressure=mean['radial_pressure'],
            outward_energy_flux=mean['outward_energy_flux'],
            angular_pressure=mean['angular_pressure'])
        return dict(fields=fields,orthonormal=mean,null=null)
