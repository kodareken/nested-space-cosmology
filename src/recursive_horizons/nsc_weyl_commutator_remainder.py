"""Exact Weyl shift commutator and conditional finite-history remainder bounds.

This connects the owned formal spatial symbol to its actual operator defect.
An auxiliary projector is never substituted for the physical evolved state.
"""
from fractions import Fraction
from math import factorial

import numpy as np
from flint import arb,ctx,fmpq
from .nsc_ks_ball_trajectory import exact_upper,restored_upper


def _matrix(value):
    a=np.asarray(value,complex)
    if a.shape!=(2,2) or not np.isfinite(a).all():raise ValueError('finite2x2 symbol coefficient required')
    return a


def _mode(value):
    if isinstance(value,(bool,np.bool_)) or int(value)!=value:raise ValueError('integer Fourier mode required')
    return int(value)


def _period(value):
    value=float(value)
    if not np.isfinite(value) or value<=0:raise ValueError('positive numerical period required')
    return value


def kinetic_commutator(kinetic,coefficient,mode,momentum,period):
    """[C*k,p]_# = k[C,p] + omega_n/2 {C,p}, exactly (C independent of z)."""
    C,P=_matrix(kinetic),_matrix(coefficient)
    k=float(momentum)
    if not np.isfinite(k):raise ValueError('finite canonical momentum required')
    shift=np.pi*_mode(mode)/_period(period)
    return k*(C@P-P@C)+shift*(C@P+P@C)


def potential_commutator(potential_modes,symbol_mode,mode,momentum,period):
    """Exact sum over supplied V_l: V_l p_{n-l}(k-s_l)-p_{n-l}(k+s_l)V_l.

    V_l=L^-1 integral_cell exp(-i omega_l z)V(z) dz; s_l=omega_l/2.
    This evaluates the supplied finite Fourier potential. Missing modes need
    a separate bound, and missing symbol modes cannot be silently truncated.
    """
    n=_mode(mode);length=_period(period);k=float(momentum)
    if not np.isfinite(k):raise ValueError('finite canonical momentum required')
    result=np.zeros((2,2),complex)
    for ell,V in potential_modes.items():
        ell=_mode(ell);V=_matrix(V);shift=np.pi*ell/length
        result+=V@_matrix(symbol_mode(n-ell,k-shift))-_matrix(symbol_mode(n-ell,k+shift))@V
    return result


def potential_taylor_commutator(potential_modes,symbol_derivative,mode,momentum,period,order):
    """Taylor truncation of the exact shifts, including the zeroth commutator.

    symbol_derivative(n,k,j) returns the actual jth k derivative, not the
    Taylor coefficient divided by j!. This is not an operator remainder.
    """
    order=_mode(order)
    if order<0:raise ValueError('nonnegative Taylor order required')
    n=_mode(mode);length=_period(period);k=float(momentum)
    if not np.isfinite(k):raise ValueError('finite canonical momentum required')
    result=np.zeros((2,2),complex)
    for ell,V in potential_modes.items():
        ell=_mode(ell);V=_matrix(V);shift=np.pi*ell/length
        for j in range(order+1):
            derivative=_matrix(symbol_derivative(n-ell,k,j))
            result+=((-shift)**j*(V@derivative)-shift**j*(derivative@V))/factorial(j)
    return result


def rho_operator_defect(symbol_rho,commutator):
    """R_rho=P_rho-[G_rho,P], where the owned KS G_rho=i H_rho."""
    return _matrix(symbol_rho)-1j*_matrix(commutator)


def _bound(value,name):
    if value is None or isinstance(value,(bool,np.bool_)):
        raise ValueError('explicit certified nonnegative '+name+' required')
    if isinstance(value,dict):x=restored_upper(value)
    elif isinstance(value,arb):x=value
    else:
        try:q=Fraction(value);x=arb(fmpq(q.numerator,q.denominator))
        except (ValueError,TypeError,OverflowError) as e:raise ValueError('finite '+name+' required') from e
    if not x.is_finite() or not x>=0:raise ValueError('finite nonnegative '+name+' required')
    return x.upper()


def shifted_taylor_remainder(absolute_half_transfer,order,potential_norm,
                             derivative_minus_segment,derivative_plus_segment,*,bits=120):
    """One-mode Frobenius remainder using op-norm(V) or its Frobenius upper.

    Both derivative inputs bound d_k^(order+1) p throughout their respective
    closed shifted momentum segments. A derivative at k alone is insufficient.
    """
    order=_mode(order)
    if order<0:raise ValueError('nonnegative Taylor order required')
    with ctx.workprec(bits):
        s=_bound(absolute_half_transfer,'half transfer');v=_bound(potential_norm,'potential norm')
        minus=_bound(derivative_minus_segment,'minus-segment derivative bound')
        plus=_bound(derivative_plus_segment,'plus-segment derivative bound')
        upper=v*s**(order+1)*(minus+plus)/factorial(order+1)
        return {'remainder_Frobenius_upper':exact_upper(upper),'derivative_order':order+1,
            'both_shifted_segments_included':True,'inputs_certified_here':False,
            'physical_local_gate':'OPEN','state_law_changed':False}


def large_transfer_remainder(tail_potential_moments,derivatives_at_momentum,
                              shifted_symbol_minus,shifted_symbol_plus,*,bits=120):
    """Bound the exact minus Taylor commutator on omitted transfer modes.

    Moment j must bound sum_omitted ||V_l|| |omega_l/2|^j. Derivative j is
    the Frobenius supremum of d_k^j p at k. The two shifted symbol bounds
    hold over all omitted modes. These are numerical splits, not source cuts.
    """
    if not len(tail_potential_moments) or len(tail_potential_moments)!=len(derivatives_at_momentum):
        raise ValueError('matching nonempty transfer moments and symbol derivatives required')
    with ctx.workprec(bits):
        moments=[_bound(v,'transfer moment') for v in tail_potential_moments]
        derivatives=[_bound(v,'symbol derivative') for v in derivatives_at_momentum]
        minus=_bound(shifted_symbol_minus,'shifted minus-symbol bound')
        plus=_bound(shifted_symbol_plus,'shifted plus-symbol bound')
        bound=moments[0]*(minus+plus)+2*sum((m*d/factorial(j) for j,(m,d) in enumerate(zip(moments,derivatives))),arb(0))
        return {'remainder_Frobenius_upper':exact_upper(bound),'Taylor_order':len(moments)-1,
            'both_transfer_signs_required_in_moments':True,'inputs_certified_here':False,
            'physical_local_gate':'OPEN','new_physical_energy_cutoff':False,
            'source_covariance_replaced':False}
