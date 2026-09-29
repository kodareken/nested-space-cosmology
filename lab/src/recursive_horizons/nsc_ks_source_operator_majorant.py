"""Continuous source-error insertion from exact characteristic row bounds.

The radius-only history preserves the principal speeds and the intrinsic
endpoint vertices. A priori bounds on its canonical two-column propagator
avoid borrowing unvalidated numerical field norms. Coverage remains explicit.
"""
from fractions import Fraction

from flint import arb,ctx

from .nsc_ks_ball_trajectory import exact_upper
from .nsc_ks_current_history_bounds import radius_bounds,value_integral_bounds
from .nsc_ks_radius_coupling_bounds import history_coupling_bounds


def _nonnegative(value,name):
    if value is None or isinstance(value,bool):raise ValueError('explicit '+name+' required')
    result=arb(value)
    if not result.is_finite() or not result>=0:raise ValueError('nonnegative '+name+' required')
    return result


def _fraction(value):
    mantissa,exponent=value.upper().man_exp()
    return Fraction(int(mantissa))*Fraction(2)**int(exponent)


def _arb_fraction(value):
    return arb(value.numerator)/arb(value.denominator)


def fundamental_difference_bounds(growth_integral,coupling_integral,coupling_z_integral,*,bits=192):
    """Unweighted two-column norms: f=sqrt(2), d and dz0.

    Along each characteristic the energy term is a phase. Off-diagonal
    coupling gives a row-norm growth exp(B). The homogeneous reference rows
    have norm one. Thus ||D||F<=sqrt(2)exp(B)M and
    ||D_z,envelope||F<=sqrt(2)exp(B)Mz. The physical derivative also has E D.
    """
    with ctx.workprec(bits):
        B=_nonnegative(growth_integral,'growth integral')
        M=_nonnegative(coupling_integral,'coupling integral')
        Mz=_nonnegative(coupling_z_integral,'coupling derivative integral')
        f=arb(2).sqrt();growth=B.exp()
        return {'reference':f,'difference':(f*growth*M).upper(),
                'difference_envelope_z':(f*growth*Mz).upper()}


def source_moment_insertion(norms,weighted_error,weighted_energy_error,*,
                            mass,absolute_angular,axial_lower,radius_lower,multiplicity,bits=192):
    """Bound the action error from S0=sum w epsilon/(2pi), S1=sum w |E|epsilon/(2pi).

    epsilon is the OPERATOR error of the unweighted upstream 2x2 covariance.
    No quadrature or source multiplicity is hidden in the propagator norms.
    Both source signs require their own epsilon before forming the moments.
    """
    with ctx.workprec(bits):
        f,d,dz=(_nonnegative(norms.get(k),k) for k in ('reference','difference','difference_envelope_z'))
        S0=_nonnegative(weighted_error,'weighted covariance error')
        S1=_nonnegative(weighted_energy_error,'weighted energy covariance error')
        m,l,a,r,mu=(_nonnegative(v,n) for v,n in zip(
            (mass,absolute_angular,axial_lower,radius_lower,multiplicity),
            ('mass','angular','axial lower','radius lower','multiplicity')))
        if not a>0 or not r>0 or not mu>0:raise ValueError('positive geometry and multiplicity required')
        density=2*f*d+d*d
        current=density*S1+(f+d)*dz*S0
        beta=(mu*current).upper()
        N=(mu*((m*m+(l/r)**2).sqrt()*density*S0+current/a)).upper()
        return {'N':exact_upper(N),'beta':exact_upper(beta),
                'carrier_and_envelope_derivatives_retained':True,
                'quadrature_error_included':False,'physical_local_gate':'OPEN'}


def history_source_insertion(family,rho_up,weighted_error,weighted_energy_error,*,
                             mass,absolute_angular,multiplicity,bits=192):
    """Exact-profile global-z bound, hence continuous on the incoming I."""
    with ctx.workprec(bits):
        m=_nonnegative(mass,'mass');ell=_nonnegative(absolute_angular,'angular')
        bounds=radius_bounds(family,rho_up)
        # On the owned slab, atan(rho)>=pi/4 and pi/4<11/14 prove a<1.
        # Then |s_up|>=rho_up-1: the fixed source lies outside the history.
        upper=bounds['rho_upper']
        a2_upper=3*((1+upper**2)*Fraction(11,14)-1)-1
        if not a2_upper<1 or bounds['rho_up']-1<bounds['normal_support']:
            raise ValueError('unchanged upstream preparation outside the history support required')
        mq,lq=_fraction(m),_fraction(ell)
        growth=value_integral_bounds(bounds,mq,lq)['K0']
        coupling=history_coupling_bounds(family,rho_up,lq)
        B,M,Mz=map(_arb_fraction,(growth,coupling['M_integral'],coupling['Mz_integral']))
        norms=fundamental_difference_bounds(B,M,Mz,bits=bits)
        result=source_moment_insertion(norms,weighted_error,weighted_energy_error,
            mass=m,absolute_angular=ell,axial_lower=_arb_fraction(bounds['axial_lower']),
            radius_lower=_arb_fraction(bounds['reference_radius_lower']),multiplicity=multiplicity,bits=bits)
        return {**result,'growth_integral_upper':exact_upper(B),
                'coupling_integral_upper':exact_upper(M),'coupling_z_integral_upper':exact_upper(Mz),
                'fundamental_norms':{k:exact_upper(v) for k,v in norms.items()},
                'whole_spatial_line':True,'continuous_on_incoming_interval':True,
                'numerical_field_norms_used':False,'physical_upstream_budget_component':None}
