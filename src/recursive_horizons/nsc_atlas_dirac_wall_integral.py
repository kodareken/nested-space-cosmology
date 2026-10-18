"""Piecewise quadrature of the triangular first-order mean wall equations.

Integrate the two supplied mass-flux laws. The remaining exact integrating
factor is rho/v, because d(log v)/d tau=acc/v. All integration is on the
existing source breakpoints. No mass or waveform trace is differentiated.
"""
from __future__ import annotations

import numpy as np


class WallPrimitive:
    """Cumulative masses and rho/v, with direct partial-cell quadrature."""
    def __init__(self,geometry,mass_rate,*,start,end,breakpoints,initial,
                 with_displacement=True,order=12):
        if order not in (8,12):
            raise ValueError('the admitted wall quadratures have order 8 or 12')
        self.geometry,self.mass_rate=geometry,mass_rate
        self.start,self.end=float(start),float(end)
        self.initial=np.asarray(initial,dtype=float)
        if (self.initial.shape!=(3,) or not np.isfinite(self.initial).all()
                or not 0<=self.start<self.end):
            raise ValueError('finite wall interval and three corner data required')
        self.with_displacement=with_displacement
        knots=np.asarray(breakpoints,dtype=float)
        if not np.isfinite(knots).all():
            raise ValueError('finite source breakpoints required')
        self.edges=np.r_[self.start,np.unique(knots[(knots>self.start)&(knots<self.end)]),self.end]
        nodes,weights=np.polynomial.legendre.leggauss(order)
        self.nodes,self.weights=(nodes+1)/2,weights/2
        width=np.diff(self.edges)
        outer=self.edges[:-1,None]+width[:,None]*self.nodes
        rate=self.mass_rate(outer)
        changes=width*np.einsum('...j,j->...',rate,self.weights)
        self.mass_edges=np.column_stack((self.initial[:2],
            self.initial[:2,None]+np.cumsum(changes,axis=-1)))
        z0=self.initial[2]/float(self.geometry(self.start)['v'])
        if with_displacement:
            mass=self._mass_at(outer)
            forcing=self._forcing(outer,mass)
            increments=width*np.einsum('...j,j->...',forcing,self.weights)
            self.z_edges=np.r_[z0,z0+np.cumsum(increments)]
        else:
            self.z_edges=np.full(len(self.edges),z0)

    def _indices(self,t):
        t=np.asarray(t,dtype=float)
        if not np.isfinite(t).all() or np.any((t<self.start)|(t>self.end)):
            raise ValueError('wall primitive query exits its interval')
        return np.minimum(np.searchsorted(self.edges,t,side='right')-1,len(self.edges)-2)

    def _mass_at(self,t):
        t=np.asarray(t,dtype=float)
        index=self._indices(t)
        left=self.edges[index]
        width=t-left
        nodes=left[...,None]+width[...,None]*self.nodes
        addition=width*np.einsum('...j,j->...',self.mass_rate(nodes),self.weights)
        return self.mass_edges[:,index]+addition

    def _forcing(self,t,mass):
        g=self.geometry(t)
        R,v,bp=g['x'],g['v'],g['beta_plus']
        return (bp*(mass[0]-mass[1])/(g.get('k_wall',2.)*R*R)+mass[0]/R)/(v*v)

    def __call__(self,t):
        t=np.asarray(t,dtype=float)
        index=self._indices(t)
        mass=self._mass_at(t)
        g=self.geometry(t)
        if self.with_displacement:
            left=self.edges[index]
            width=t-left
            nodes=left[...,None]+width[...,None]*self.nodes
            forcing=self._forcing(nodes,self._mass_at(nodes))
            z=self.z_edges[index]+width*np.einsum('...j,j->...',forcing,self.weights)
            rho=g['v']*z
        else:
            # A smooth control's radius is supplied by its shared worldtube;
            # this slot is not another material-wall equation.
            rho=np.zeros_like(t)+self.initial[2]
        return np.array([mass[0],mass[1],rho])
