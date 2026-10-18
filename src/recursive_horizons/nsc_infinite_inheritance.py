"""Positive finite-energy limit of the inherited tree law.

The proof and its infinite domain are in docs/nsc-infinite-inheritance.md.
This calculator evaluates its finite Schur recursion and analytic depth
bound. It does not supply a cap--wall metric or a cosmological Lambda.
"""
from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class InheritedTree:
    local: np.ndarray
    link: np.ndarray
    omega: float
    branching: int

    def __post_init__(self):
        H=np.asarray(self.local,dtype=complex);B=np.asarray(self.link,dtype=complex)
        omega=float(self.omega)
        if (H.ndim!=2 or H.shape[0]==0 or H.shape[0]!=H.shape[1] or B.shape!=H.shape
                or not np.isfinite(H).all() or not np.isfinite(B).all()
                or not np.array_equal(H,H.conj().T)):
            raise ValueError('a finite Hermitian local matrix and matching finite link are required')
        if (not np.isfinite(omega) or omega<=1 or isinstance(self.branching,bool)
                or not isinstance(self.branching,int) or self.branching<2):
            raise ValueError('omega>1 and integer branching at least two are required')
        eig=np.linalg.eigvalsh(H);beta=float(np.linalg.norm(B,2))
        c=math.sqrt(self.branching/omega)*beta;h=float(eig[0])
        if not h>2*c:
            raise ValueError('the theorem requires lambda_min(H)>2*sqrt(b/omega)*norm(B)')
        gamma=(h+math.sqrt((h-2*c)*(h+2*c)))/2
        H=H.copy();B=B.copy();H.setflags(write=False);B.setflags(write=False)
        object.__setattr__(self,'local',H);object.__setattr__(self,'link',B)
        object.__setattr__(self,'omega',omega)
        object.__setattr__(self,'lower_form_constant',h-2*c)
        object.__setattr__(self,'upper_form_constant',float(eig[-1])+2*c)
        object.__setattr__(self,'gamma',gamma)
        object.__setattr__(self,'contraction',(c/gamma)**2)

    @property
    def counting_exponent(self):
        return math.log(self.branching)/math.log(self.omega)

    def response(self,parameter,depth):
        """Root (D_N+lambda)^-1 for depths0,...,N.

        The returned bound controls infinite-depth truncation in exact
        arithmetic. Eigensolver/rounding errors are not included in it.
        """
        lam=float(parameter)
        if (not np.isfinite(lam) or lam<0 or isinstance(depth,bool)
                or not isinstance(depth,int) or not 0<=depth<=10000):
            raise ValueError('finite nonnegative parameter and depth0..10000 required')
        parameters=[lam]
        for _ in range(depth):
            parameters.append(parameters[-1]/self.omega)
        size=len(self.local);eye=np.eye(size);R=np.zeros((size,size),complex)
        for value in reversed(parameters):
            denominator=self.local+value*eye-self.branching/self.omega*self.link@R@self.link.conj().T
            R=np.linalg.solve(denominator,eye)
        if self.contraction:
            log_bound=(depth+1)*math.log(self.contraction)-math.log(self.gamma)
            # Keep an underflowed depth bound visible; zero is reserved
            # for an exactly uncoupled hierarchy.
            bound=max(math.exp(log_bound),float(np.nextafter(0.,1.)))
        else:
            log_bound=None;bound=0.
        return dict(response=R,analytic_depth_bound=bound,log_depth_bound=log_bound,
            depth=depth,parameter=lam,contraction=self.contraction)

    def finite_matrix(self,depth):
        """Small independent controls; never expand a large hierarchy."""
        if isinstance(depth,bool) or not isinstance(depth,int) or depth<0:
            raise ValueError('nonnegative integer depth required')
        if depth>7:
            raise ValueError('dense verification is limited to256 coordinates; use the Schur recursion')
        vertices=(self.branching**(depth+1)-1)//(self.branching-1)
        size=vertices*len(self.local)
        if size>256:
            raise ValueError('dense verification is limited to256 coordinates; use the Schur recursion')
        levels=[0]*vertices
        result=np.zeros((size,size),complex);d=len(self.local)
        for vertex in range(vertices):
            if vertex:
                parent=(vertex-1)//self.branching;levels[vertex]=levels[parent]+1
                p=slice(parent*d,(parent+1)*d);v=slice(vertex*d,(vertex+1)*d)
                edge=self.omega**levels[parent]*self.link
                result[p,v]=edge;result[v,p]=edge.conj().T
            v=slice(vertex*d,(vertex+1)*d)
            result[v,v]=self.omega**levels[vertex]*self.local
        return result,np.repeat(self.omega**np.asarray(levels),d)
