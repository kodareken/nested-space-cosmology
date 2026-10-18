"""Positive finite-energy limit of the inherited tree law.

The proof and its infinite domain are in docs/nsc-infinite-inheritance.md.
This calculator evaluates its finite Schur recursion and analytic depth
bound. It does not supply a cap--wall metric or a cosmological Lambda.
"""
from dataclasses import dataclass
import math

import numpy as np


def geometric_inherited_tree(branching,ratio,length=1.,*,spatial_dimension=1,cross_section=1.):
    """Inherited wave law from edge lengths, sections and half-edge masses.

    At depth n, length is L*q^n and section is A*q^((dimension-1)*n).
    Each vertex receives half of each incident edge's volume as its
    kinetic weight; the root has a fixed parent edge. Dimension1 with
    section1 recovers the original length-weighted action. The form
    closure supplies the boundary choice at infinity.
    """
    q,L,A=float(ratio),float(length),float(cross_section)
    if (isinstance(branching,bool) or not isinstance(branching,int) or branching<2
            or isinstance(spatial_dimension,bool) or not isinstance(spatial_dimension,int)
            or spatial_dimension<1 or not np.isfinite([q,L,A]).all()
            or not 0<q<1 or L<=0 or A<=0):
        raise ValueError('integer branching>=2, dimension>=1, positive length/section and0<q<1 required')
    if spatial_dimension==3 and q==1/branching:
        raise ValueError('the critical b*q=1 case is outside the strict positive-form theorem')
    mass=.5*A*L*(q**(-spatial_dimension)+branching)
    H=A*(branching+q**(2-spatial_dimension))/(L*mass)
    B=-A/(L*mass*q**(.5*spatial_dimension))
    return InheritedTree(np.array([[H]]),np.array([[B]]),q**-2,branching)


def homogeneous_dilation_source(kinetic_reference,gradient_reference,scale=1.,coordinate_volume=1.):
    """Homogeneous metric variation of a positive three-dimensional field action.

    M(a)=a^3 M0, K(a)=a K0 at fixed fractional coordinates. Inputs are
    T0=pi^T M0^-1 pi/2 and G0=phi^T K0 phi/2 at the specified state.
    Variation holds canonical field data fixed. These are geometric energy
    units (G=c=1), with a dimensionless scale and volume V=a^3*Vc.
    No source, common-parent geometry or vacuum identification is supplied.
    """
    from .nsc_atlas_gravitational_reduction import homogeneous_energy_stress

    values=(kinetic_reference,gradient_reference,scale,coordinate_volume)
    if any(isinstance(x,(bool,np.bool_)) for x in values):
        raise ValueError('finite nonnegative energies and positive scale/volume required')
    T0,G0,a,Vc=map(float,values)
    if not np.isfinite([T0,G0,a,Vc]).all() or min(T0,G0)<0 or min(a,Vc)<=0:
        raise ValueError('finite nonnegative energies and positive scale/volume required')
    kinetic=T0/a**3;gradient=a*G0;volume=a**3*Vc
    pressure=(kinetic-gradient/3)/volume
    result=homogeneous_energy_stress(kinetic+gradient,-pressure,volume)
    result.update(reference_kinetic=T0,reference_gradient=G0,scale=a,
        coordinate_volume=Vc,kinetic_energy=kinetic,gradient_energy=gradient,
        rho_plus_3p=4*kinetic/volume,
        flat_deceleration=(2*kinetic/(kinetic+gradient) if kinetic+gradient else None),
        scope='Homogeneous three-dimensional kinetic-plus-gradient lift of the inherited field action; no full parent cosmology or quantum vacuum.')
    return result


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
