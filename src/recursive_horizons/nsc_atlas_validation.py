"""Floating replay allowances; never scientific error or source tolerances."""
import numpy as np
from .nsc_atlas_field_preparation import radius_product


def source_jet_roundoff_tolerances(x, geometry, *, second=False):
    """Absolute replay allowances for (g, g_x), optionally also g_xx.

    This is a floating evaluation allowance, not a field error estimate. The
    source label subtracts two O(1) null coordinates before dividing by their
    smaller span. Propagate eight binary64 rounding units in that expression
    through b(z) and b'(z)*z_x; a relative test on g_x alone fails near b'=0.
    Source values, preparation, and all numerical accuracy targets are unchanged.
    """
    x = np.asarray(x, dtype=float)
    U0, Va, Vb = (geometry[key] for key in ("U0", "Va", "Vb"))
    V = radius_product(x)/U0
    z = 4*(V-(Va+Vb)/2)/(Vb-Va)
    active = (abs(z) < 1) & (x > geometry["support_x"][0]) & (x < geometry["support_x"][1])
    safe = np.where(active, z, 0.)
    bz = -8*safe*(1-safe*safe)**3
    bzz = 8*(1-safe*safe)**2*(7*safe*safe-1)
    zx = -4*x*np.exp(x/.8)/(.8**2*U0*(Vb-Va))
    eps = np.finfo(float).eps
    dz = 8*eps*(4*(abs(V)+abs((Va+Vb)/2))/abs(Vb-Va)+abs(z))
    value = abs(bz)*dz+16*eps
    first = abs(bzz*zx)*dz+16*eps*np.maximum(1., abs(bz*zx))
    result = (np.where(active, value, 0.), np.where(active, first, 0.))
    if second:
        bzzz = 48*safe*(1-safe*safe)*(3-7*safe*safe)
        zxx = -4*np.exp(x/.8)*(1+x/.8)/(.8**2*U0*(Vb-Va))
        curvature = (abs(bzzz*zx*zx+bzz*zxx)*dz
                     +32*eps*np.maximum(1., abs(bzz*zx*zx)+abs(bz*zxx)))
        result += (np.where(active, curvature, 0.),)
    # The branch is explicitly quiet in both the saved and current expression.
    return result


