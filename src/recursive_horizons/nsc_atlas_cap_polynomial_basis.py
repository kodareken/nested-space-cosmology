"""Regular odd Legendre coordinates for the same polynomial cap space.

These basis operations perform no radial solve or physical time evolution.
Dividing the odd polynomial analytically by xi preserves its centre limit.
"""
from numbers import Integral

import numpy as np

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_polar_cap import _scalar


def _count(value):
    if isinstance(value,(bool,np.bool_)) or not isinstance(value,Integral) or not 1 <= value <= 16:
        raise ValueError('basis requires one to sixteen odd Legendre modes')
    return int(value)


def odd_legendre_values(count, xi, *, budget):
    """Return P, P', Q=P/xi and dQ/d(xi^2) without division by xi."""
    count = _count(count)
    xi = _scalar(xi,'xi')
    if not 0 <= xi <= 1:
        raise ValueError('basis point must lie in the closed unit interval')
    budget = _budget(budget)
    budget.check(32*count+32)
    t = xi*xi
    quotient,derivative = np.empty(count),np.empty(count)
    even,q,even_t,q_t = 1.,1.,0.,0.
    quotient[0],derivative[0] = q,q_t
    for j in range(1,count):
        new_even = ((4*j-1)*t*q-(2*j-1)*even)/(2*j)
        new_even_t = ((4*j-1)*(q+t*q_t)-(2*j-1)*even_t)/(2*j)
        new_q = ((4*j+1)*new_even-2*j*q)/(2*j+1)
        new_q_t = ((4*j+1)*new_even_t-2*j*q_t)/(2*j+1)
        even,q,even_t,q_t = new_even,new_q,new_even_t,new_q_t
        quotient[j],derivative[j] = q,q_t
    return dict(value=xi*quotient,derivative_xi=quotient+2*t*derivative,
                quotient=quotient,quotient_derivative_t=derivative)


def odd_legendre_operators(count, Z, *, budget):
    """Mass, stiffness, dilation and trace functionals in one basis.

    Derivative bounds are sup|P'_p|<=p(p+1)/2 on [0,1]. Their use with
    floating coefficients requires outward arithmetic in the caller.
    """
    count = _count(count)
    Z = _scalar(Z,'Z',positive=True)
    budget = _budget(budget)
    budget.check(96*count*count+64*count)
    powers = np.arange(1,2*count,2,dtype=float)
    lower = np.minimum(powers[:,None],powers[None,:])
    try:
        with np.errstate(over='raise',invalid='raise',divide='raise'):
            diagonal = Z/(2*powers+1)
            stiffness = lower*(lower+1)/(2*Z)
            if np.any(diagonal <= 0) or not np.all(np.isfinite(stiffness)):
                raise ValueError('basis matrices are not representable')
            dilation = np.diag(powers)+np.triu(np.broadcast_to((2*powers+1)[:,None],(count,count)),1)
            centre = odd_legendre_values(count,0.,budget=budget)['derivative_xi']
            wall_derivative = powers*(powers+1)/2
            return dict(powers=powers,gram=np.diag(diagonal),stiffness=stiffness,
                        dilation=dilation,wall_value=np.ones(count),
                        wall_derivative_xi=wall_derivative,centre_derivative_xi=centre,
                        derivative_bound_xi=wall_derivative.copy(),
                        exact_gram_condition=(4*count-1)/3)
    except FloatingPointError as error:
        raise ValueError('basis matrices are not representable') from error
