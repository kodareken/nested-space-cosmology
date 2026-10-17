"""Select the interior polynomial at a split radial panel's endpoints.

The physical query coordinate is unchanged. A rounded C2 spline can have a
small second-derivative mismatch at a knot. Quadrature on each smooth panel
uses that panel's polynomial limit, including its Runge--Kutta endpoints.
"""
import numpy as np


def history_piece_indices(profile, labels, *, budget):
    labels=np.asarray(labels,dtype=float)
    if labels.shape!=(2,) or not np.all(np.isfinite(labels)):
        raise ValueError('two finite history labels required')
    profile._domain(labels)
    budget.check(2)
    indices=np.searchsorted(profile._labels,labels,side='right')-1
    return np.minimum(indices,len(profile._labels)-2).astype(int)


def panel_profile_jet(profile, labels, indices, *, budget):
    """Evaluate owned pieces, permitting only roundoff-size knot crossings."""
    labels=np.asarray(labels,dtype=float)
    selected=np.asarray(indices)
    if selected.shape!=(2,) or selected.dtype.kind not in 'iu':
        raise ValueError('two integer history-piece indices required')
    actual=history_piece_indices(profile,labels,budget=budget)
    knots=profile._labels
    if np.any(selected < -1) or np.any(selected > len(knots)-2):
        raise ValueError('history-piece index outside the finite profile')
    result=profile.jet(labels,budget=budget)
    objects=(profile._spline,*profile._derivatives[:3])
    for j,(expected,found) in enumerate(zip(selected,actual)):
        expected,found=int(expected),int(found)
        if expected==found:
            continue
        if abs(expected-found)!=1:
            raise ValueError('ambiguous history-piece ownership at radial endpoint')
        knot=knots[max(expected,found)]
        tolerance=64*np.spacing(max(1.,abs(float(labels[j])),abs(float(knot))))
        if abs(labels[j]-knot)>tolerance:
            raise ValueError('history-piece selection crosses a resolved knot')
        if expected==-1 and not profile.quiet_extension:
            raise ValueError('history piece requires a certified quiet extension')
        for name,polynomial in zip(('F','F_prime','F_second','F_third'),objects):
            result[name][j]=(0. if expected==-1 else
                float(np.polyval(polynomial.c[:,expected],labels[j]-knots[expected])))
            budget.check(len(polynomial.c))
    if not all(np.all(np.isfinite(value)) for value in result.values()):
        raise ValueError('nonfinite one-sided history jet')
    return result
