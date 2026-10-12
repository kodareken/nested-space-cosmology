"""Uniform-in-z mixed geometry jets from certified profile derivative norms."""
from math import factorial
from types import MappingProxyType
from flint import arb,arb_series,ctx
from . import nsc_ks_ball_geometry as G
from .nsc_ks_ball_operator import AnalyticRadiusFamily
from .nsc_ks_spacetime_geometry_jets import KSGeometryJets,background_coordinate_time_jets
from .nsc_ks_profile_identity import profile_identity
from .nsc_bloch_trace_bound import _norm_upper,_positive_interval


def uniform_axial_geometry_jets(family,rho,w_derivative_uppers,u_derivative_uppers,*,
                                period_left,period_length,axial_profile_identity,
                                order=9,bits=160):
    """Enclose every T/z jet at every z in the numerical cell.

    Profile inputs are derivative norms, including amplitudes, not Taylor
    coefficients. They must belong to the supplied analytic profile identity.
    No missing U/profile norm is replaced with zero.
    """
    order=G._order(order)
    if len(w_derivative_uppers)<order+1 or len(u_derivative_uppers)<order+1:
        raise ValueError('all physical profile derivative uppers are required')
    if axial_profile_identity!=profile_identity(family,include_normal_window=False):
        raise ValueError('profile moments belong to a different analytic history')
    model=AnalyticRadiusFamily(family)
    with ctx.workprec(bits),G._JetWork(order):
        rho=G._rho_ball(rho)
        left=G._binary_arb(period_left,'period origin')
        length=_positive_interval(period_length,'numerical period',bits=bits)
        for amplitude,direction in zip(model.metric.amplitudes,model.metric.directions):
            for profile in (direction.w,direction.U):
                if amplitude and any(profile.coefficients):
                    center=arb(float(profile.center));outer=arb(float(profile.outer))
                    if not left<center-outer or not center+outer<left+length:
                        raise ValueError('the numerical cell must contain the compact profile support')
        w=[_norm_upper(v,'physical w derivative',bits=bits) for v in w_derivative_uppers[:order+1]]
        u=[_norm_upper(v,'physical U derivative',bits=bits) for v in u_derivative_uppers[:order+1]]
        bg=background_coordinate_time_jets(rho,order)
        shift=arb_series([G.shift_series(rho,0)[0],arb(1)],prec=order+1)
        window=G.plateau_series(-shift,model.inner,model.outer,order)
        first=window*shift;third=window*shift**3/6
        pairs=tuple((t,z) for t in range(order+1) for z in range(order+1-t))
        fields={name:{pair:(bg[name][pair[0]] if pair[1]==0 else arb(0)) for pair in pairs}
                for name in ('a','inv_a','inv_a2','r_reference')}
        fields['N']={pair:arb(int(pair==(0,0))) for pair in pairs}
        fields['beta']={pair:arb(0) for pair in pairs}
        fields['r']={}
        for t,z in pairs:
            # The real symmetric balls enclose the signed profile derivatives.
            wz=arb(0,(w[z]/factorial(z)).upper())
            uz=arb(0,(u[z]/factorial(z)).upper())
            fields['r'][(t,z)]=fields['r_reference'][(t,z)]+first[t]*wz+third[t]*uz
        if not fields['r'][(0,0)]>0:
            raise G.SubdivisionNeeded('uniform radius positivity needs a smaller time box or tighter profile norms')
        if any(not v.is_finite() for values in fields.values() for v in values.values()):
            raise ArithmeticError('nonfinite uniform geometry jet')
        z_cell=left.union(left+length)
        return KSGeometryJets(order,bits,MappingProxyType({k:MappingProxyType(v) for k,v in fields.items()}),
            profile_identity(family),rho,z_cell)
