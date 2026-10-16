"""A source-free scalar derivative front in the common constrained cap metric.

This is one position/amplitude RHS, not a bulk evolution. Its supplied scalar,
mass and lapse must be the continuous traces of the same piecewise solution.
The caller owns front collisions and axis/material-wall interactions.
"""
import numpy as np

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_cap_radial_reference import reference_potential
from .nsc_atlas_polar_cap import _scalar


def wavefront_rhs(radius, reduced_jump, *, orientation, wall_radius, alpha,
                  mass, lapse, field, clock_rate, budget):
    """Return one internal-front stage in the caller's time label x.

    orientation=+1 transports j=[p-q], and -1 transports j=[p+q], where
    p=r*Pi and q=(r*f)_r. Brackets are right minus left. mass is the
    vacuum-subtracted n, with A=1-r^2-2*alpha*n/r. clock_rate=dT/dx
    uses the inherited cap clock; it is not reset inside this component.

    The stage is valid only for an isolated front with continuous f,n,delta
    and source, piecewise field equations, and no singular scalar source.
    Positive A at this supplied point is not an all-radius certificate.
    """
    budget=_budget(budget)
    before=budget.evaluations
    orientation=_scalar(orientation,'orientation')
    if orientation not in (-1,1):
        raise ValueError('orientation must be -1 or +1')
    r=_scalar(radius,'radius',positive=True)
    R=_scalar(wall_radius,'wall radius',positive=True)
    if not 0<r<R<1:
        raise ValueError('front must lie strictly between the axis and wall below r=1')
    jump=_scalar(reduced_jump,'reduced jump')
    al=_scalar(alpha,'alpha',nonnegative=True)
    n=_scalar(mass,'mass')
    delta=_scalar(lapse,'lapse')
    f=_scalar(field,'field')
    tx=_scalar(clock_rate,'clock rate',positive=True)
    budget.check(32)
    try:
        with np.errstate(over='raise',invalid='raise',divide='raise'):
            A=1-r*r-2*al*n/r
            scale=float(np.exp(delta))
            if not A>0 or not scale>0:
                raise ValueError('front requires a positive metric and lapse exponential')
            potential=reference_potential(np.array([f]),al,budget=budget)
            P=float(potential['Pdiff'][0])
            B=A*scale
            Br=scale*(-2*r+2*al*n/(r*r)-al*r*P)
            rx=tx*orientation*B
            jx=-tx*orientation*Br*jump
            dpi=jump/(2*r)
            dphi=-orientation*dpi
    except FloatingPointError as error:
        raise ValueError('nonfinite characteristic front arithmetic') from error
    if not np.all(np.isfinite([A,scale,B,Br,rx,jx,dpi,dphi])) or not B>0:
        raise ValueError('nonfinite or unresolved characteristic front stage')
    budget.check()
    return dict(radius_rate=float(rx),jump_rate=float(jx),
                Pi_jump=float(dpi),Phi_jump=float(dphi),A=float(A),B=float(B),
                B_r=float(Br),evaluations=budget.evaluations-before,
                domain='one isolated internal wavefront stage; common bulk metric required',
                numerical_error_enclosure=False)
