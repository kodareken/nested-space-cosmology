"""Scalar potential arithmetic for the fixed-cap radial hot path."""
import numpy as np

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_polar_cap import _scalar


def reference_potential_pair(f, alpha, *, budget):
    """The same Pdiff and Pdiff_f at one point, without unused gap arrays.

    This is scalar arithmetic for the constrained radial engine. The full
    reference_potential remains the owner of the implementation-gap and
    series-remainder diagnostics. Each point keeps its original charge.
    """
    budget = _budget(budget)
    f = _scalar(f, 'field')
    alpha = _scalar(alpha, 'alpha', nonnegative=True)
    budget.check(12)
    try:
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            scale = np.sqrt(alpha/6)
            if alpha > 0 and scale == 0:
                raise ValueError('reference alpha scale underflows')
            q = scale*f
            ch, sh = np.cosh(q), np.sinh(q)
            sinc = sh/q if q != 0 else 1.
            P = f*f*sinc*sinc*((ch*ch+1)/2)
            U = 2*f*sinc*ch**3
    except (FloatingPointError, OverflowError) as error:
        raise ValueError('nonfinite reference potential') from error
    if not np.isfinite(P) or not np.isfinite(U):
        raise ValueError('nonfinite reference potential')
    return float(P), float(U)
