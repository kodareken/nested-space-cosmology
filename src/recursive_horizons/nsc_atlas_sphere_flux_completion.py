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
