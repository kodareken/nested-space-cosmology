"""New smooth16D completion of the retained scalar sector, with supplied inputs.

The pair-cycle six-form replaces the old averaged T term; it is not added to
it. This identifies neither an inheritance rule nor an SI scale or child.
The internal space is (S3)^4, so no round-S12 spectral result transfers.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import cache
from itertools import combinations
import math

import numpy as np
from scipy.optimize import brentq

from .nsc_atlas_flux_conversion import ScalarFlux

PAIRS=tuple(combinations(range(4),2))


@dataclass(frozen=True,slots=True)
class FluxStressCheck:
    """Exact coefficients of q² in reference orthonormal components."""
    norm_squared: F
    internal_contraction: tuple
    internal_stress: tuple
    external_metric_coefficient: F


def flux_components():
    """Sorted exterior components of F6/q; each omega has pointwise norm1."""
    return tuple((tuple(range(3*i,3*i+3))+tuple(range(3*j,3*j+3)),F(1)) for i,j in PAIRS)


@cache
def flux_stress_check():
    """Contract sparse five-form interiors, retaining exterior signs exactly.

    C_ab=<i_a F,i_b F> equals (1/5!)F_aI F_bI. Thus T_ab=C_ab-|F|²/2
    and T_munu=-|F|² g_munu/2, since F has no external indices.
    """
    components=flux_components();interiors=[{} for _ in range(12)]
    for indices,value in components:
        for position,index in enumerate(indices):
            key=indices[:position]+indices[position+1:]
            interiors[index][key]=interiors[index].get(key,F(0))+(-1)**position*value
    norm=sum((value*value for _,value in components),F(0))
    contraction=tuple(tuple(sum((value*interiors[b].get(key,F(0)) for key,value in interiors[a].items()),F(0))
                            for b in range(12)) for a in range(12))
    stress=tuple(tuple(contraction[a][b]-(norm/2 if a==b else 0) for b in range(12)) for a in range(12))
    return FluxStressCheck(norm,contraction,stress,-norm/2)


@dataclass(frozen=True,slots=True)
class ActionReductionRecord:
    dimensions: int
    internal_dimensions: int
    topology: str
    radius_squared: F
    reference_volume: float
    kappa16: float
    Lambda16: F
    mass16: F
    q: float
    q_squared_times_reference_volume: F
    A: F
    B: F
    T: F
    K: F
    scalar_mass_squared_prefactor: F
    shape_kinetic_metric_prefactor: F
    action: str
    metric: str
    scalar_normalization: str
    form_normalization: str
    Maxwell_identity: str
    assumptions: tuple
    scope: str
    overall_normalization: str = ('The common quantum action factor1/g remains unspecified. '
        'kappa4=1 is the retained classical geometrized normalization, not a selection of kappa_nat/Lunit^2 or metres.')


@cache
def action_reduction_record():
    """Immutable action normalization and supplied matching, in model units."""
    r2=F(8,5);V0=(2*math.pi**2*float(r2)**1.5)**4
    return ActionReductionRecord(16,12,'(S3_r)^4',r2,V0,V0,F(7),F(20),
        (5/2)/math.sqrt(24*V0),F(25,96),F(7),F(15,2),F(25,32),F(84),F(400),F(3),
        '(R16-2Lambda16)/(2kappa16) - (dPhi)^2/2 - m16^2 Phi^2/2 - F6^2/(2*6!)',
        'ds16^2=exp(-12s) ds4E^2+exp(2s) gamma12; kappa4=kappa16/V0=1',
        'Phi=chi/sqrt(V0), internal-constant neutral scalar',
        'F6=q sum_(i<j) omega_i wedge omega_j; omega_i fixed-reference norm-one volume forms',
        '*16F6=+/-exp(-24s) vol4 wedge *gammaF6; dF6=d*16F6=0 for smooth s(x)',
        ('Einstein gravity and canonical neutral scalar/form, matter outside the gravity prefactor',
         'Fixed product topology and equal-magnitude constant pair-cycle flux periods; patchwise A5',
         'Common breathing field; no extra charged sources or form couplings',
         'Supplied couplings and units; kappa16=V0 and q refer to Lunit=1'),
        'New higher-dimensional completion of the retained scalar sector; replaces the old T term. '
        'No inheritance/dimension/SI selection, round-S12 spectrum transfer, full-mode stability or child claim.')


def _finite_s(s):
    value=np.asarray(s,dtype=float)
    if not np.all(np.isfinite(value)):raise ValueError('finite breathing field required')
    return value


def _value(value):return float(value) if np.ndim(value)==0 else value


def shape_potential(s,beta=(0.,0.,0.,0.)):
    """Full bulk+curvature+flux potential, sum(beta_i)=0; no trace projection.

    beta's last axis has four entries; its leading axes broadcast with s.
    Scalar mass is separate and independent of these volume-preserving shapes.
    """
    s=_finite_s(s);beta=np.asarray(beta,dtype=float)
    if beta.ndim<1 or beta.shape[-1]!=4 or not np.all(np.isfinite(beta)):
        raise ValueError('four finite relative-radius beta entries required')
    if np.any(abs(np.sum(beta,axis=-1))>64*np.finfo(float).eps*np.maximum(1.,np.sum(abs(beta),axis=-1))):
        raise ValueError('sum(beta)=0 is required; shape trace is not silently removed')
    curvature=np.sum(np.exp(-2*beta),axis=-1)
    flux=sum(np.exp(-6*(beta[...,i]+beta[...,j])) for i,j in PAIRS)
    return _value(7*np.exp(-12*s)-(15/8)*np.exp(-14*s)*curvature+(25/192)*np.exp(-24*s)*flux)


def scalar_mass_squared(s):return _value(400*np.exp(-12*_finite_s(s)))


def shape_mass_squared(s):
    """Three degenerate local masses after division by the shape metric3I."""
    s=_finite_s(s)
    return _value(-2.5*np.exp(-14*s)+3.125*np.exp(-24*s))


def shape_threshold():return math.log(5/4)/10


@dataclass(frozen=True,slots=True)
class IntegerFluxEmbedding:
    period_unit_squared: F
    integer_periods: tuple
    squared_integer_norm: int
    reduced_flux_energy: F
    weighted_vertex_degrees: tuple
    equal_radius_shape_acceleration_coefficients: tuple


def integer_flux_embedding(period_unit_squared,integer_periods):
    """Conditional compact charge lattice; the wrapping vector is input.

    Eflux=unit²*sum(Nij²)/2. At beta=beta_dot=0, the constrained shape
    equation gives beta_i''=unit²(degree_i-sum(N²)/2)*exp(-24s).
    This is initial forcing, not a nonlinear history or a selected embedding.
    """
    unit=F(period_unit_squared);periods=tuple(integer_periods)
    if unit<=0 or len(periods)!=6 or any(isinstance(v,bool) or not isinstance(v,int) for v in periods):
        raise ValueError('positive squared period unit and six integer pair periods required')
    norm=sum(v*v for v in periods);degrees=[0]*4
    for (i,j),v in zip(PAIRS,periods):degrees[i]+=v*v;degrees[j]+=v*v
    force=tuple(unit*(F(v)-F(norm,2)) for v in degrees)
    return IntegerFluxEmbedding(unit,periods,norm,unit*norm/2,tuple(degrees),force)


@dataclass(frozen=True,slots=True)
class CombinedParentRecord:
    scalar_sector: ScalarFlux
    bare_mass_squared: F
    mixing_squared: F
    LambdaBare: F
    q_reference: F
    phi_star_reference: float
    B: F
    T: F
    K: F
    reference_volume: float
    kappa16: float
    scalar_normalization: str
    top_form_normalization: str
    action: str
    scope: str
    overall_normalization: str


@cache
def combined_record():
    """Restore the original charged sector, replacing neutral7/mass²400."""
    geometry=action_reduction_record()
    sector=ScalarFlux(epsilon=.25,omega=20.,bare_density=float(F(215,32)))
    return CombinedParentRecord(sector,F(100),F(300),F(215,32),F(3,2),sector.minimum(1.5),
        geometry.B,geometry.T,geometry.K,geometry.reference_volume,geometry.kappa16,
        'Phi=phi/sqrt(V0); chi=phi-phi_star(3/2) only in the reference sector',
        'qD=q/sqrt(V0); fixed-displacement Routhian gives (qD+sqrt(300)*Phi)^2/2',
        'Gravity (R16-2LambdaBare)/(2kappa16); canonical scalar bare mass²100 and F16 mixing²300, '
        'plus the same canonical F6. Matter/forms outside gravity prefactor.',
        'Supplied charged parent completion, not a new preparation: replaces neutral Lambda7/mass20. '
        'F6 replaces the old T term once; no inherited state, SI, wrapping selection or child is derived.',
        geometry.overall_normalization)


def charged_potential(s,phi,q,beta4=(0.,0.,0.,0.)):
    """Full charged potential; q is top-form displacement, not pair-flux q6."""
    s=_finite_s(s);phi=np.asarray(phi,dtype=float);q=np.asarray(q,dtype=float)
    if not np.all(np.isfinite(phi)) or not np.all(np.isfinite(q)):
        raise ValueError('finite scalar field and top-form displacement required')
    # Remove the neutral7 contribution before inserting ScalarFlux. Neither
    # its vacuum energy nor the existing six-form stabilizer is counted twice.
    return _value(shape_potential(s,beta4)+np.exp(-12*s)*(combined_record().scalar_sector.potential(phi,q)-7))


def volume_thresholds():
    """Exact fifth powers/formulas and numerical stationary-sector boundaries."""
    AM5=F(25,4)**5*F(8,5);AC5=F(175,24)**5*F(14,15)
    AM=float(AM5)**.2;AC=float(AC5)**.2
    return dict(minkowski_A_fifth=AM5,critical_A_fifth=AC5,
        minkowski_A_exact='(25/4)*(8/5)^(1/5)',critical_A_exact='(175/24)*(14/15)^(1/5)',
        minkowski_A=AM,critical_A=AC,turn_x_fifth=F(14,15),minkowski_x_fifth=F(8,5),
        positive_minimum_abs_q_interval=(math.sqrt(8*(AM-float(F(215,32)))),math.sqrt(8*(AC-float(F(215,32))))))


def all_volume_extrema(q):
    """All finite extrema after relaxing phi, not an off-equilibrium history.

    x=exp(-2s)>0 solves6A-7Bx+12Tx^6=0. Its unique positive turning
    point supplies both exhaustive brackets; no spatial/parameter scan.
    """
    if isinstance(q,(bool,np.bool_)) or not math.isfinite(float(q)):
        raise ValueError('finite real top-form displacement required')
    q=float(q);model=combined_record();A=model.scalar_sector.vacuum(q);B=float(model.B);T=float(model.T)
    turn=(7*B/(72*T))**.2
    polynomial=lambda x:6*A-7*B*x+12*T*x**6
    minimum=polynomial(turn);tol=128*np.finfo(float).eps*max(1.,6*A,7*B*turn)
    if minimum>tol:return ()
    degenerate=abs(minimum)<=tol
    roots=(turn,) if degenerate else (brentq(polynomial,0.,turn,xtol=5e-15),brentq(polynomial,turn,2*turn,xtol=5e-15))
    rows=[]
    for x in roots:
        stiffness=0. if degenerate else 4*x**7*(-7*B+72*T*x**5)
        vacuum=x**7*(B/6-T*x**5)
        sign='zero' if abs(vacuum)<=tol else ('positive' if vacuum>0 else 'negative')
        rows.append(dict(q=q,A=A,x=x,s=-math.log(x)/2,phi_star=model.scalar_sector.minimum(q),
            vacuum=vacuum,vacuum_sign=sign,stiffness=stiffness,canonical_mass_squared=stiffness/float(model.K),
            kind='degenerate' if degenerate else ('minimum' if stiffness>0 else 'maximum'),
            polynomial_residual=polynomial(x)))
    return tuple(rows)


def half_integer_classification():
    """All q=n/2 sectors, via exact rational fifth-power inequalities."""
    bounds=volume_thresholds();A=lambda q:F(215,32)+q*q/8
    return dict(step=F(1,2),positive_vacuum_minimum_sectors=(F(-3,2),F(3,2)),
        negative_vacuum_minimum_abs_q_at_most=F(1),no_volume_extrema_abs_q_at_least=F(2),
        minkowski_fifth_power_bracket=(A(F(1))**5,bounds['minkowski_A_fifth'],A(F(3,2))**5),
        merger_fifth_power_bracket=(A(F(3,2))**5,bounds['critical_A_fifth'],A(F(2))**5),
        scope='Scalar-relaxed, common-radius volume sectors only; supplied couplings and q=n/2 lattice. '
              'Not a classification of all unequal-radius vacua or a selection of the inputs.')


def wrapped_brane_couplings(s):
    """Material14-brane wraps allM12; a smooth sigma0 horizon is not its source."""
    scale=np.exp(-6*_finite_s(s));sigma=_value(.2*scale)
    return dict(sigma=sigma,e_can=_value(.5*scale),omega=_value(20*scale),
        q_out_minus_in=F(1,2),radion_normal_jump=_value(-.2*scale/14),
        phi_direct_normal_jump=0.,shape_direct_source=(0.,0.,0.,0.),
        radion_jump_variable='s; multiply by sqrt(84) for canonical S',
        displacement_jump_variable='Reference top-form q; instantaneous canonical charge is e_can(s)',
        scope='Full-wrap material charge source; scalar has no direct brane coupling and traceless shapes preserve wrapped volume.')


def _log_radius_data(u4,q):
    u=np.asarray(u4,dtype=float);q=np.asarray(q,dtype=float)
    if u.ndim<1 or u.shape[-1]!=4 or not np.all(np.isfinite(u)) or not np.all(np.isfinite(q)):
        raise ValueError('four finite log radii and finite top-form displacement required')
    x=np.exp(-2*u);t=x**3
    others=np.stack([sum(t[...,j] for j in range(4) if j!=i) for i in range(4)],axis=-1)
    A=combined_record().scalar_sector.vacuum(q)
    W=A-(15/8)*np.sum(x,axis=-1)+(25/192)*sum(t[...,i]*t[...,j] for i,j in PAIRS)
    return np.exp(-3*np.sum(u,axis=-1)),W,x,t,others


def log_radius_potential(u4,q=1.5):
    """Full charged potential at phi_star(q), with u_i=s+beta_i.

    The four log radii are independent coordinates; q is held fixed, not a
    new flux-sector selection. Leading batch axes are supported.
    """
    P,W,_,_,_=_log_radius_data(u4,q)
    return _value(P*W)


def log_radius_gradient(u4,q=1.5):
    """Exact derivative with respect to all four unprojected log radii."""
    P,W,x,t,others=_log_radius_data(u4,q)
    Wi=2*(15/8)*x-6*(25/192)*t*others
    return P[...,None]*(Wi-3*W[...,None])


def log_radius_hessian(u4,q=1.5):
    """Exact coordinate Hessian; canonical masses use the kinetic metric."""
    P,W,x,t,others=_log_radius_data(u4,q)
    Wi=2*(15/8)*x-6*(25/192)*t*others
    Wij=36*(25/192)*(t[..., :,None]*t[...,None,:])
    diagonal=np.arange(4)
    Wij[...,diagonal,diagonal]=-4*(15/8)*x+36*(25/192)*t*others
    return P[...,None,None]*(Wij-3*(Wi[..., :,None]+Wi[...,None,:])+9*W[...,None,None])


def log_radius_kinetic_metric():
    """G=3I+(9/2)11^T, so common breathing has K84 and shapes have3I."""
    return 3*np.eye(4)+4.5*np.ones((4,4))


@cache
def symmetric_radius_branches():
    """Exact positive-ratio isolation in fixed1+3/2+2 families at q=3/2.

    The two isotropic extrema and permutation-inequivalent family members
    are included. This is not an enumeration of all critical points, a KK
    index, or a statement of transition dominance. Masses are numerical
    evaluations at the isolated algebraic stationary points.
    """
    import sympy as sp
    from scipy.linalg import eigh

    model=combined_record();q=1.5
    if model.scalar_sector.vacuum(q)!=7 or model.B!=F(15,2) or model.T!=F(25,32):
        raise ValueError('symmetric branch polynomials require the fixed reference coefficients')
    t=sp.symbols('t');A=sp.Integer(7);a=sp.Rational(5,4);S=t*t+t+1
    P13=t**3+28*t*t+28*t+15
    P22=5*t**6+14*t**5+14*t**4+6*t**3+14*t*t+14*t+5
    polynomials={'1+3':sp.Poly(5*4**5*A**5*S**6-12*a**5*P13**5,t),
                 '2+2':sp.Poly(5*2**5*A**5*(S*(t**3+1))**6-24*a**5*P22**5,t)}
    isolations=[];branches=[];G=log_radius_kinetic_metric()

    def point_record(family,u,**metadata):
        u=np.asarray(u,float);V=log_radius_potential(u,q)
        masses,modes=eigh(log_radius_hessian(u,q),G)
        channels=[]
        if V>0:
            for mass in masses:
                ell=0;degrees=[];degeneracies=[]
                while ell*(ell+3)<-3*mass/V:
                    degrees.append(ell);degeneracies.append((ell+1)*(ell+2)*(2*ell+3)//6);ell+=1
                channels.append(dict(ell=tuple(degrees),degeneracy=tuple(degeneracies)))
        count=sum(sum(c['degeneracy']) for c in channels) if V>0 else None
        return dict(branch=family,q=q,log_radii=tuple(map(float,u)),s=float(np.mean(u)),
            beta=tuple(map(float,u-np.mean(u))),V=V,
            max_gradient=float(np.max(abs(log_radius_gradient(u,q)))),
            canonical_mass_squared=tuple(map(float,masses)),
            canonical_mass_modes=tuple(tuple(map(float,row)) for row in modes),
            homogeneous_negative_count=int(np.count_nonzero(masses<0)),
            S4_scalar_harmonic_negative_count=count,S4_negative_harmonics=tuple(channels),
            negative_mode=tuple(map(float,modes[:,0])) if masses[0]<0 else None,
            retained_radius_local_minimum=bool(masses[0]>0),**metadata)

    for row in all_volume_extrema(q):
        branches.append(point_record('isotropic',[row['s']]*4,volume_kind=row['kind'],
            ratio=1.,ratio_interval=(F(1),F(1)),ratio_polynomial_degree=None))
    for family,polynomial in polynomials.items():
        polynomial=polynomial.clear_denoms()[1].primitive()[1]
        intervals=[(bounds,multiplicity) for bounds,multiplicity in polynomial.intervals(eps=sp.Rational(1,10**14)) if bounds[1]>0]
        exact=lambda r:F(int(r.p),int(r.q))
        isolated=tuple(dict(interval=(exact(lo),exact(hi)),multiplicity=int(multiplicity))
            for (lo,hi),multiplicity in intervals)
        isolations.append(dict(branch=family,degree=polynomial.degree(),
            integer_polynomial_coefficients=tuple(map(int,polynomial.all_coeffs())),
            positive_root_intervals=isolated,
            representative_rule='All positive roots' if family=='1+3' else 'Keep t<1; t>1 is its S4 factor-permutation copy'))
        for item in isolated:
            lo,hi=item['interval'];ratio=float((lo+hi)/2)
            if family=='2+2' and ratio>1:continue
            scale=(12/(5*(ratio*ratio+ratio+1)))**.2 if family=='1+3' else (
                24/(5*(ratio*ratio+ratio+1)*(ratio**3+1)))**.2
            x=np.array([ratio*scale,scale,scale,scale] if family=='1+3' else [ratio*scale,ratio*scale,scale,scale])
            branches.append(point_record(family,-.5*np.log(x),ratio=ratio,
                ratio_interval=(lo,hi),ratio_polynomial_degree=polynomial.degree()))
    return dict(q_reference=F(3,2),isolations=tuple(isolations),branches=tuple(branches),
        S4_scalar_index_rule='Sum degeneracy(l)=(l+1)(l+2)(2l+3)/6 where l(l+3)<-3m²/V, only for V>0.',
        scope='Only isotropic,1+3 and2+2 radius families modulo factor permutations; numerical canonical masses and retained-scalar S4 indices, not all critical points, full KK indices or decay dominance.')
