"""Mixed coordinate-T/axial jets of the existing analytic radius history.

The geometry is unchanged. These interval jets retain the higher derivatives
needed by the exact Weyl defect; they are not a field, source, or gate proof.
"""
from dataclasses import dataclass
from math import factorial
from types import MappingProxyType

from flint import arb,arb_series,ctx
from . import nsc_ks_ball_geometry as G
from .nsc_ks_ball_operator import AnalyticRadiusFamily
from .nsc_ks_profile_identity import profile_identity


def _compose(coefficients,increment,order):
    """Composition at the same basepoint; increment has exact zero constant."""
    if increment[0]!=0:raise ValueError('formal increment must have exact zero constant')
    result=arb_series([0],prec=order+1)
    for coefficient in reversed(G._coeffs(coefficients,order)):
        result=result*increment+coefficient
    return result


def background_coordinate_time_jets(rho,order):
    """rho_T=-a(rho), composed from the owned rho-jets without changing clock."""
    with G._JetWork(order):
        rho=G._rho_ball(rho);bg=G.background_series(rho,order)
        # q is rho(T0+tau)-rho(T0), whose constant is identically zero even
        # when the common basepoint rho(T0) ranges over a real ball.
        q=arb_series([0],prec=order+1)
        for n in range(order):
            rhs=-_compose(bg.a,q,order)
            coefficients=G._coeffs(q,order);coefficients[n+1]=rhs[n]/(n+1)
            q=arb_series(coefficients,prec=order+1)
        return {'rho_increment':q,
            'a':_compose(bg.a,q,order),'r_reference':_compose(bg.r,q,order),
            'inv_a':_compose(bg.inv_a,q,order),'inv_a2':_compose(bg.inv_a2,q,order)}


@dataclass(frozen=True)
class KSGeometryJets:
    order:int
    bits:int
    coefficients:object
    history_identity:str
    rho_domain:object
    z_domain:object

    def derivative(self,field,time_order=0,axial_order=0):
        if (isinstance(time_order,bool) or isinstance(axial_order,bool)
                or int(time_order)!=time_order or int(axial_order)!=axial_order
                or min(time_order,axial_order)<0 or time_order+axial_order>self.order):
            raise ValueError('mixed derivative exceeds the retained geometry order')
        with ctx.workprec(self.bits):
            return self.coefficients[field][(int(time_order),int(axial_order))]*factorial(int(time_order))*factorial(int(axial_order))


def spacetime_geometry_jets(family,rho,z,*,order=9,bits=160):
    """Enclose derivative/(t! z!) for the same g at every basepoint in rho,z.

    N=1,beta=0 and a(T) keep the reference chart. Radius is exactly
    r_ref(T)+chi(s)[s*w(z)+s^3*U(z)/6], with s=T-T_Sigma and fixed profiles.
    A ball crossing an unresolved smooth-window wall requests subdivision.
    """
    order=G._order(order)
    if order<1:raise ValueError('at least first order geometry jets required')
    if isinstance(bits,bool) or int(bits)!=bits or bits<80:raise ValueError('at least80 precision bits required')
    model=AnalyticRadiusFamily(family)
    with ctx.workprec(int(bits)),G._JetWork(order):
        rho=G._rho_ball(rho);z=G._binary_arb(z,'z')
        background=background_coordinate_time_jets(rho,order)
        s0=G.shift_series(rho,0)[0]
        # Coordinate T, not rho, is the formal time variable here.
        shift=arb_series([s0,arb(1)],prec=order+1)
        window=G.plateau_series(-shift,model.inner,model.outer,order)
        first=window*shift;third=window*shift**3/6
        w,U=model.profile_series(z,order)
        pairs=tuple((t,x) for t in range(order+1) for x in range(order+1-t))
        coefficients={}
        for name in ('a','inv_a','inv_a2','r_reference'):
            coefficients[name]={pair:(background[name][pair[0]] if pair[1]==0 else arb(0)) for pair in pairs}
        coefficients['N']={pair:arb(int(pair==(0,0))) for pair in pairs}
        coefficients['beta']={pair:arb(0) for pair in pairs}
        coefficients['r']={pair:coefficients['r_reference'][pair]+first[pair[0]]*w[pair[1]]+third[pair[0]]*U[pair[1]] for pair in pairs}
        if not coefficients['r'][(0,0)]>0:raise G.SubdivisionNeeded('radius positivity requires a smaller geometry box')
        for values in coefficients.values():
            if any(not value.is_finite() for value in values.values()):raise ArithmeticError('nonfinite mixed geometry jet')
        frozen=MappingProxyType({name:MappingProxyType(values) for name,values in coefficients.items()})
        return KSGeometryJets(order,int(bits),frozen,profile_identity(family),rho,z)
