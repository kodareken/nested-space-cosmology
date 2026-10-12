"""Uniform inverse-energy remainder for the homogeneous source covariance.

The exact compact metric and horizon vacuum branch are retained. This is an
upstream covariance remainder, not the changed-history field UV remainder.
"""
from dataclasses import dataclass

import numpy as np
from flint import acb,arb,arb_series,ctx

from .nsc_massive_jost_phase_bound import _series_context
from .nsc_metric_horizon_frame import horizon_q
from .nsc_subgap_source_covariance import metric_H


def _real(value,name):
    # arb(None) means zero; an absent scientific input must never take that path.
    if value is None or isinstance(value,bool):
        raise ValueError('explicit real '+name+' required')
    return arb(value)


def _truncate(series,degree):
    return arb_series([series[j] for j in range(degree+1)],prec=degree+1)


def _derivative(series,degree):
    return arb_series([(j+1)*series[j+1] for j in range(degree+1)],prec=degree+1)


def coefficient_jets(delta,H,L,order):
    """P_j and Z_j for w=t sum(P_j/E^j), nz=sum(Z_j/E^j), t^2=delta.

    Input jets must provide at least order+1 derivatives. Known precision decreases
    by one at each recurrence; no unavailable derivative is treated as zero.
    Returns through order+1, the coefficient governing the exact defect.
    """
    if type(order) is not int or not 1<=order<=16:
        raise ValueError('integer inverse-energy order in [1,16] required')
    degree=min(delta.prec,H.prec,L[0].prec,L[1].prec)-1
    if degree<order+1:
        raise ValueError('complete metric/coefficient derivative jets required')
    zero=arb_series([0],prec=degree+1)
    ps=[(zero,zero)];zs=[arb_series([1],prec=degree+1)]
    for n in range(1,order+2):
        known=degree-n+1
        if n==1:
            pair=L
        else:
            px,py=ps[-1]
            pair=(-H/2*(py/2+delta*_derivative(py,known))+L[0]*zs[-1],
                   H/2*(px/2+delta*_derivative(px,known))+L[1]*zs[-1])
        ps.append(tuple(_truncate(p,known) for p in pair))
        norm=sum((delta*(ps[j][0]*ps[n-j][0]+ps[j][1]*ps[n-j][1])
                  +zs[j]*zs[n-j] for j in range(1,n)),zero)
        zs.append(_truncate(-norm/2,known))
    return ps,zs


@dataclass(frozen=True)
class VacuumSourceExpansion:
    horizon_hint: object
    mass: object
    angular: object
    order: int=8
    bits: int=192
    metric_terms: int=48

    def __post_init__(self):
        if type(self.order) is not int or not 1<=self.order<=16:
            raise ValueError('integer inverse-energy order in [1,16] required')
        if type(self.bits) is not int or self.bits<80:
            raise ValueError('directed precision >=80 required')
        if (type(self.metric_terms) is not int or self.metric_terms<8
                or self.metric_terms<=self.order+1):
            raise ValueError('metric tail must cover every derivative order')
        with ctx.workprec(self.bits):
            m,l=_real(self.mass,'mass'),_real(self.angular,'angular label')
            if not m.is_finite() or not m>=0 or not l.is_finite():
                raise ValueError('finite nonnegative mass and real angular label required')
            object.__setattr__(self,'mass',m);object.__setattr__(self,'angular',l)
            object.__setattr__(self,'q',horizon_q(self.horizon_hint,bits=self.bits))

    def jets(self,delta,*,remainder_jet_order=1):
        """Coordinate jets, retaining the requested derivatives of P_(M+1)."""
        if type(remainder_jet_order) is not int or remainder_jet_order<1:
            raise ValueError('positive remainder jet order required')
        degree=self.order+remainder_jet_order
        if degree>=self.metric_terms:
            raise ValueError('metric tail must cover the requested remainder jets')
        with _series_context(self.bits,degree):
            x=arb_series([_real(delta,'compact distance'),1],prec=degree+1)
            H=metric_H(self.q,-x,degree,terms=self.metric_terms)
            radius=1/(self.q-x).sin()
            if not H[0]>0 or not radius[0]>0:
                raise ValueError('positive homogeneous interior metric required')
            L=(self.mass*H.sqrt()*radius,-self.angular*H.sqrt())
            ps,zs=coefficient_jets(x,H,L,self.order)
            return H,ps,zs

    def target_distance(self,rho):
        with ctx.workprec(self.bits):
            rho=_real(rho,'target rho')
            if not rho.is_finite() or not rho>=1:
                raise ValueError('homogeneous preparation rho>=1 required')
            delta=self.q-arb.pi()/2-rho.atan()
            if not delta>0 or not delta<arb('0.4'):
                raise ValueError('target must be inside the owned horizon collar')
            return delta

    def remainder_constant(self,rho,*,cells=32):
        """C with ||n_exact-n_M|| <= C/E^M, uniformly for every E>0.

        The defect integral is 4 integral_0^sqrt(delta) |P_(M+1)|/H dt.
        Directed whole-cell ranges supply a Darboux upper, not sampled maxima.
        """
        if type(cells) is not int or cells<1:
            raise ValueError('positive integer cell count required')
        with ctx.workprec(self.bits):
            top=self.target_distance(rho).sqrt()
            total=arb(0);width=top/cells
            for j in range(cells):
                left=top*j/cells;right=top*(j+1)/cells
                delta=(left*left).union(right*right)
                H,ps,_=self.jets(delta)
                magnitude=sum((p[0].abs_upper()**2 for p in ps[-1]),arb(0)).sqrt()
                total+=4*width*magnitude/H[0].lower()
            if not total.is_finite():raise ArithmeticError('nonfinite vacuum remainder')
            return total.upper()

    def approximate_bloch(self,rho,energy):
        """Interval evaluation of the finite Hermitian covariance polynomial."""
        with ctx.workprec(self.bits):
            energy=_real(energy,'energy')
            if not energy.is_finite() or not energy>0:
                raise ValueError('strictly positive real energy required')
            delta=self.target_distance(rho)
            _,ps,zs=self.jets(delta)
            transverse=[delta.sqrt()*sum((ps[j][k][0]/energy**j
                for j in range(1,self.order+1)),arb(0)) for k in (0,1)]
            longitudinal=sum((zs[j][0]/energy**j for j in range(self.order+1)),arb(0))
            return (*transverse,longitudinal)

    def covariance_error(self,constant,energy,*,kappa,omega):
        """Vacuum remainder plus the full source occupation/coherence upper.

        Unitarity of the homogeneous frame and |R|^2+T=1 give
        ||Q_source-Q_vac|| <= f+n_in+coherence. The inherited incoming
        occupation cutoff only reduces n_in, so its exponential is safe.
        """
        with ctx.workprec(self.bits):
            C,E,k,w=(_real(value,name) for value,name in
                     zip((constant,energy,kappa,omega),('remainder','energy','kappa','omega')))
            if (not all(v.is_finite() for v in (C,E,k,w)) or not C>=0
                    or not E>0 or not k>0 or not w>0):
                raise ValueError('nonnegative remainder and positive source scales required')
            vacuum=(C/(2*E**self.order)).upper()
            f=(-2*arb.pi()*E/k).exp()
            coherent=(-arb.pi()*E/k).exp()
            incoming=(-2*arb.pi()*E/(w*k)).exp()
            return {'vacuum_operator_error':vacuum,'thermal_coherent_error':(f+coherent+incoming).upper(),
                    'source_operator_error':(vacuum+f+coherent+incoming).upper()}

    def compare_original(self,rho,energy,columns,covariance,constant,*,kappa,omega,
                         negative_partner=False):
        """Bound physical source error of existing, unweighted 2x3 columns.

        The finite covariance polynomial is only an intermediate comparison;
        original source columns remain unchanged. The caller authenticates
        this positive channel or its opposite-angular negative-energy partner.
        """
        F=np.asarray(columns,complex);C=np.asarray(covariance,complex)
        if (F.shape!=(2,3) or C.shape!=(3,3) or not np.isfinite(F).all()
                or not np.isfinite(C).all() or not np.array_equal(C,C.conj().T)):
            raise ValueError('finite original 2x3 columns and Hermitian 3x3 covariance required')
        if type(negative_partner) is not bool:
            raise ValueError('explicit signed-partner choice required')
        with ctx.workprec(self.bits):
            approximate=self.approximate_bloch(rho,energy)
            if negative_partner:
                approximate=(approximate[0],-approximate[1],-approximate[2])
            f=[[acb(v) for v in row] for row in F]
            c=[[acb(v) for v in row] for row in C]
            Q=[[sum((f[i][a]*c[a][b]*f[j][b].conjugate()
                     for a in range(3) for b in range(3)),acb(0))
                for j in range(2)] for i in range(2)]
            trace=(Q[0][0]+Q[1][1]).real
            original=(2*Q[0][1].real,-2*Q[0][1].imag,(Q[0][0]-Q[1][1]).real)
            distance=((abs(1-trace)+sum(((a-b).abs_upper()**2
                for a,b in zip(approximate,original)),arb(0)).sqrt())/2).upper()
            errors=self.covariance_error(constant,energy,kappa=kappa,omega=omega)
            return {**errors,'polynomial_to_stored_operator_distance':distance,
                    'original_source_operator_error':(distance+errors['source_operator_error']).upper()}
