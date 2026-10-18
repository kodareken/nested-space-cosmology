"""Global gapped-Dirac operator estimates used before the stress tail.

These require genuine input bounds. They are not a computed state, a
mean-stress bound or a substitute for the prescribed preparation evolution.
Derivations and domains are in docs/nsc-atlas-dirac-probe.md.
"""
import math


def _number(value,*,positive=False):
    value=float(value)
    if not math.isfinite(value) or value<0 or (positive and value==0):
        raise ValueError('finite '+('positive' if positive else 'nonnegative')+' bound required')
    return value


def _positive_exp(log_value):
    # An exponentially small positive bound must not become an exact zero
    # when high angular sectors underflow the floating representation.
    try:
        value=math.exp(log_value)
    except OverflowError:
        return math.inf
    return value if value else math.nextafter(0.,1.)


def balanced_band_bound(source_norm,order=1):
    """Norm bound for A^(n/2)*L^(-n)(Y)*B^(n/2), finite orders 1..8."""
    if isinstance(order,bool) or not isinstance(order,int) or not 1<=order<=8:
        raise ValueError('finite inverse order 1..8 required')
    return math.ldexp(_number(source_norm),-order)


def double_band_bounds(source_norm,gap_positive,gap_negative):
    """Full one-sided graph bounds for the double Sylvester inverse."""
    Y=_number(source_norm);a=_number(gap_positive,positive=True);b=_number(gap_negative,positive=True)
    return dict(left_energy=Y/(4*b),right_energy=Y/(4*a),operator_norm=Y/(a+b)**2)


def resolvent_separation_bound(gap,imaginary_frequency,separation,*,theta=.5):
    """Off-diagonal instantaneous radial resolvent bound; units of length."""
    m=_number(gap,positive=True);d=_number(separation);z=float(imaginary_frequency)
    if not math.isfinite(z) or not 0<theta<1:
        raise ValueError('finite imaginary frequency and 0<theta<1 required')
    rho=math.hypot(m,z)
    return _positive_exp(-theta*rho*d-math.log1p(-theta)-math.log(rho))


def frozen_projector_difference_bound(q_lower,delta_q,separation,angular_k,*,theta=.5):
    """Localized instantaneous projector difference, never a stress bound.

    Both potentials must be >=q_lower and agree throughout the detector
    neighborhood. Their difference is <=delta_q, at distance >=separation.
    """
    q=_number(q_lower,positive=True);difference=_number(delta_q);d=_number(separation)
    if (isinstance(angular_k,bool) or not isinstance(angular_k,int) or not angular_k
            or not 0<theta<1):
        raise ValueError('nonzero integer angular sector and 0<theta<1 required')
    if difference==0:
        return 0.
    log_bound=math.log(difference)-math.log(2)-2*math.log1p(-theta)-math.log(q)-2*theta*d*q*abs(angular_k)
    # The norm difference of two orthogonal projectors is also at most one.
    return min(1.,_positive_exp(log_bound))
