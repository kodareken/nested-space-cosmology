"""One finite C2 profile and its reflected de Sitter scalar background.

This is a computational subtraction profile, not a new physical source or
an imposed finite-amplitude wall waveform. Values and all derivatives come
from the same Hermite profile. An optional degree-seven first active piece
preserves a supplied quartic onset. Optional interior value/first jets enrich
selected ordinary pieces to degree seven; unselected pieces stay quintic. A C2 knot gives finite
regular centre first jets; jumps of its third derivative need not give a
radially differentiable Pi there. No arbitrary smoothness claim is made.
"""
import math
from collections.abc import Mapping

import numpy as np
from scipy.interpolate import BPoly, PPoly


MAX_POINTS = 513
SMALL_RADIUS = 1e-2


def _array(value, name):
    if isinstance(value, (list, tuple)) and any(
            isinstance(v, (bool, np.bool_)) for v in np.asarray(value, dtype=object).flat):
        raise ValueError(name+" must contain real numbers, not booleans")
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf" or raw.ndim > 1 or not 1 <= raw.size <= MAX_POINTS:
        raise ValueError(name+" requires at most 513 real scalar or one-dimensional values")
    result = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(result)):
        raise ValueError(name+" must be finite")
    return result


def _budget(budget):
    if not callable(getattr(budget, "check", None)):
        raise ValueError("a shared cumulative budget.check(count) is required")
    budget.check()
    return budget


def _moment(low, high, power):
    # Factoring the difference of powers retains narrow split intervals.
    return (high-low)*sum(high**(power-j)*low**j for j in range(power+1))/(power+1)


def _onset_bernstein(quartic, value, first, second, width):
    """Degree-seven controls, valid for float or exact/interval scalar inputs."""
    return [0,0,0,0,quartic*width**4/35,
            value-2*width*first/7+width*width*second/42,value-width*first/7,value]


def _enriched_bernstein(v0,d0,s0,v1,d1,s1,h,eta,value,first):
    """Degree-seven endpoint/interior controls in float, rational or Arb arithmetic."""
    base=[v0,v0+h*d0/5,v0+2*h*d0/5+h*h*s0/20,
          v1-2*h*d1/5+h*h*s1/20,v1-h*d1/5,v1]
    def evaluate(row):
        while len(row)>1: row=[(1-eta)*a+eta*b for a,b in zip(row[:-1],row[1:])]
        return row[0]
    dv=value-evaluate(base.copy())
    dx=first-evaluate([5*(b-a)/h for a,b in zip(base[:-1],base[1:])])
    psi=eta**3*(1-eta)**3;psi_prime=3*eta**2*(1-eta)**2*(1-2*eta)
    a=dv/psi;b=(h*dx-psi_prime*a)/psi
    controls=base
    for degree in (5,6):
        controls=[controls[0],*[(i*controls[i-1]+(degree+1-i)*controls[i])/(degree+1)
            for i in range(1,degree+1)],controls[-1]]
    controls[3]+=(a-eta*b)/35
    controls[4]+=(a+(1-eta)*b)/35
    return controls


def _interior_data(value,labels,onset):
    if value is None: return None
    names={'indices','labels','values','first'}
    if not isinstance(value,Mapping) or set(value)!=names:
        raise ValueError('interior jets require indices, labels, values and first')
    data={key:_array(value[key],'interior '+key).copy() for key in names}
    indices=data['indices']
    if (any(row.ndim!=1 or row.shape!=indices.shape for row in data.values())
            or np.any(indices!=np.floor(indices)) or np.any(np.diff(indices)<=0)
            or np.any(indices<0) or np.any(indices>=len(labels)-1)):
        raise ValueError('interior jets require matching rows and ordered distinct interval indices')
    indices=indices.astype(int);data['indices']=indices
    if (np.any(data['labels']<=labels[indices]) or np.any(data['labels']>=labels[indices+1])):
        raise ValueError('interior query must lie strictly inside its own interval')
    if onset is not None and np.any(indices<=onset['index']):
        raise ValueError('interior jets cannot replace the certified quiet prefix or onset piece')
    for row in data.values(): row.flags.writeable=False
    return data


class CharacteristicProfile:
    """Finite C2 Hermite data, with optional quiet extension and quartic onset.

    Two to 513 strictly increasing labels carry F,F',F''. Below the first
    label the profile is zero only when all three initial jets are exactly
    zero. Upper extrapolation is always rejected. Public query batches have
    at most 513 points. The caller owns one cumulative budget; there is no
    internal budget reset or physical preparation/evolution.
    Optional onset={index,quartic_coefficient} replaces only that first
    active interval by a degree-seven piece. All supplied jets up to its
    left endpoint must be zero; its fourth right derivative is 24 times
    the coefficient. Preceding intervals remain exactly quiet.
    Optional interior_jets={indices,labels,values,first} adds one supplied
    value/first pair per selected ordinary interval at its actual binary label.
    End jets remain unchanged. Quiet/onset replacement is rejected.
    """

    def __init__(self, labels, values, derivatives, second_derivatives, *, onset=None,interior_jets=None):
        arrays = [_array(v, name) for v, name in (
            (labels, "labels"), (values, "values"), (derivatives, "derivatives"),
            (second_derivatives, "second_derivatives"))]
        if (any(v.ndim != 1 for v in arrays) or len(arrays[0]) < 2
                or any(v.shape != arrays[0].shape for v in arrays)
                or np.any(np.diff(arrays[0]) <= 0)):
            raise ValueError("profile needs 2..513 matching vectors with increasing labels")
        labels, values, derivatives, second_derivatives = (v.copy() for v in arrays)
        self._onset = None
        if onset is not None:
            if not isinstance(onset,Mapping) or set(onset)!={'index','quartic_coefficient'}:
                raise ValueError('onset requires index and quartic_coefficient')
            index=onset['index'];quartic=onset['quartic_coefficient']
            if (isinstance(index,(bool,np.bool_)) or not isinstance(index,(int,np.integer))
                    or not 0<=index<len(labels)-1):
                raise ValueError('onset index must identify an existing active interval')
            if (isinstance(quartic,(bool,np.bool_)) or not isinstance(quartic,(int,float,np.integer,np.floating))
                    or not math.isfinite(float(quartic)) or float(quartic)==0.):
                raise ValueError('onset quartic coefficient must be finite and nonzero')
            if any(np.any(row[:index+1]!=0.) for row in (values,derivatives,second_derivatives)):
                raise ValueError('onset and preceding supplied jets must be exactly quiet')
            self._onset=dict(index=int(index),quartic_coefficient=float(quartic))
        self._interior=_interior_data(interior_jets,labels,self._onset)
        try:
            with np.errstate(over="raise", divide="raise", invalid="raise", under="ignore"):
                data = np.column_stack((values, derivatives, second_derivatives))
                self._spline = PPoly.from_bernstein_basis(
                    BPoly.from_derivatives(labels, data, extrapolate=False))
                if self._onset is not None:
                    i=self._onset['index'];h=labels[i+1]-labels[i]
                    controls=np.array(_onset_bernstein(self._onset['quartic_coefficient'],
                        values[i+1],derivatives[i+1],second_derivatives[i+1],h),dtype=float)
                    if controls[4]==0. or not np.all(np.isfinite(controls)):
                        raise ValueError('onset coefficient is unresolved in Bernstein coordinates')
                    piece=PPoly.from_bernstein_basis(BPoly(controls[:,None],labels[i:i+2],extrapolate=False))
                    coefficients=np.pad(self._spline.c,((2,0),(0,0)))
                    coefficients[:,i]=piece.c[:,0]
                    if coefficients[3,i]==0.:
                        raise ValueError('onset quartic coefficient is unresolved in power coordinates')
                    self._spline=PPoly(coefficients,labels,extrapolate=False)
                if self._interior is not None:
                    coefficients=np.pad(self._spline.c,((8-len(self._spline.c),0),(0,0)))
                    for j,i in enumerate(self._interior['indices']):
                        h=labels[i+1]-labels[i];eta=(self._interior['labels'][j]-labels[i])/h
                        controls=np.array(_enriched_bernstein(values[i],derivatives[i],second_derivatives[i],
                            values[i+1],derivatives[i+1],second_derivatives[i+1],h,eta,
                            self._interior['values'][j],self._interior['first'][j]),dtype=float)
                        piece=PPoly.from_bernstein_basis(BPoly(controls[:,None],labels[i:i+2],extrapolate=False))
                        coefficients[:,i]=piece.c[:,0]
                    self._spline=PPoly(coefficients,labels,extrapolate=False)
        except (FloatingPointError, OverflowError) as error:
            raise ValueError("nonfinite characteristic profile coefficients") from error
        if not np.all(np.isfinite(self._spline.c)):
            raise ValueError("nonfinite characteristic profile coefficients")
        self._labels = labels
        self._labels.setflags(write=False)
        self._derivatives = [self._spline.derivative(order) for order in range(1, len(self._spline.c))]
        self.quiet_extension = bool(values[0] == derivatives[0] == second_derivatives[0] == 0.)
        self.domain = {"first_label": float(labels[0]), "last_label": float(labels[-1]),
                       "quiet_below": self.quiet_extension}

    @property
    def onset(self):
        """Copy of the one-sided onset specification, if supplied."""
        return None if self._onset is None else self._onset.copy()

    @property
    def interior_jets(self):
        """Original optional value/first jets at their actual binary labels."""
        return None if self._interior is None else {key:row.copy() for key,row in self._interior.items()}

    @property
    def breakpoints(self):
        """A copy of every polynomial knot, including both finite endpoints."""
        return self._labels.copy()

    def _domain(self, labels):
        if (np.any(labels > self._labels[-1])
                or (not self.quiet_extension and np.any(labels < self._labels[0]))):
            raise ValueError("profile query exits the supplied finite labels")

    def _evaluate(self, labels, budget):
        self._domain(labels)
        budget.check(int(labels.size))
        quiet = labels < self._labels[0]
        safe = np.maximum(labels, self._labels[0])
        objects = (self._spline, *self._derivatives[:3])
        result = {name: np.where(quiet, 0., spline(safe)) for name, spline in zip(
            ("F", "F_prime", "F_second", "F_third"), objects)}
        if not all(np.all(np.isfinite(v)) for v in result.values()):
            raise ValueError("nonfinite characteristic profile jet")
        return result

    def jet(self, labels, *, budget):
        """F,F',F'',F''' of this one interpolant; one charge per label."""
        budget = _budget(budget)
        values = _array(labels, "profile labels")
        result = self._evaluate(values, budget)
        budget.check()
        return result

    def _split_integrals(self, T, distance, budget):
        """Exact polynomial integrals of F', F'', and z F'' on -1<=z<=1."""
        internal = self._labels[(self._labels > T-distance) & (self._labels < T+distance)]
        cuts = np.r_[-1., (internal-T)/distance, 1.]
        budget.check(3*(len(cuts)-1))
        first = second = weighted_second = 0.
        for low, high in zip(cuts[:-1], cuts[1:]):
            midpoint = T+distance*(low+high)/2
            index = int(np.searchsorted(self._labels, midpoint, side="right")-1)
            if index < 0:
                # _domain already established that this lower interval is quiet.
                continue
            index = min(index, len(self._labels)-2)
            coefficient = self._spline.c[:, index]
            at = T-self._labels[index]
            degree=len(coefficient)-1
            budget.check(4*max(0,degree-5))
            moments = [_moment(low, high, order) for order in range(degree)]
            for order in range(degree):
                derivative = float(np.polyval(np.polyder(coefficient, order+1), at))
                first += derivative*distance**order/math.factorial(order)*moments[order]
            for order in range(degree-1):
                derivative = float(np.polyval(np.polyder(coefficient, order+2), at))
                factor = derivative*distance**order/math.factorial(order)
                second += factor*moments[order]
                weighted_second += factor*moments[order+1]
        return first, second, weighted_second

    def background(self, T, r, *, budget):
        """Reflected reduced scalar and regular canonical first jets.

        A0=1-r², chi=F(T+atanh r)-F(T-atanh r),
        p=(F'_+-F'_-)/A0, q=(F'_++F'_-)/A0,
        f=chi/r, Pi=p/r, Phi=(q-f)/r.

        Small radii use exact, knot-split polynomial integrals instead of
        subtracting nearby values. At r=0 the finite limits are returned.
        Charges: one per radius/clock pair, one per directly evaluated label,
        or three per integrated polynomial piece in the small-radius branch.
        """
        budget = _budget(budget)
        times, radii = np.broadcast_arrays(_array(T, "T"), _array(r, "r"))
        if times.size > MAX_POINTS or np.any(radii < 0) or np.any(radii >= 1):
            raise ValueError("background needs at most 513 points with 0<=r<1")
        shape = times.shape
        times, radii = times.ravel(), radii.ravel()
        distance = np.arctanh(radii)
        labels = np.stack((times-distance, times+distance))
        self._domain(labels)
        budget.check(len(times))
        A = 1-radii*radii
        result = {key: np.empty_like(times) for key in ("chi", "p", "q", "f", "Pi", "Phi")}
        ordinary = radii > SMALL_RADIUS
        if np.any(ordinary):
            jets = self._evaluate(labels[:, ordinary], budget)
            rr, aa = radii[ordinary], A[ordinary]
            chi = jets["F"][1]-jets["F"][0]
            p = (jets["F_prime"][1]-jets["F_prime"][0])/aa
            q = (jets["F_prime"][1]+jets["F_prime"][0])/aa
            f = chi/rr
            for key, value in (("chi", chi), ("p", p), ("q", q),
                               ("f", f), ("Pi", p/rr), ("Phi", (q-f)/rr)):
                result[key][ordinary] = value
        axis = radii == 0
        if np.any(axis):
            jets = self._evaluate(times[axis], budget)
            result["chi"][axis] = result["p"][axis] = result["Phi"][axis] = 0.
            result["f"][axis] = result["q"][axis] = 2*jets["F_prime"]
            result["Pi"][axis] = 2*jets["F_second"]
        for i in np.flatnonzero(~ordinary & ~axis):
            rr, dd, aa = radii[i], distance[i], A[i]
            first, second, weighted = self._split_integrals(float(times[i]), float(dd), budget)
            ratio = dd/rr
            # (r/(1-r²)-atanh(r))/r², retaining its O(r) value.
            coefficient = rr*(2/3+rr*rr*(4/5+rr*rr*(6/7+rr*rr*(8/9+rr*rr*10/11))))
            result["chi"][i] = dd*first
            result["p"][i] = dd/aa*second
            result["q"][i] = (first+dd*weighted)/aa
            result["f"][i] = ratio*first
            result["Pi"][i] = ratio/aa*second
            result["Phi"][i] = coefficient*first+ratio/aa*weighted
        if not all(np.all(np.isfinite(v)) for v in result.values()):
            raise ValueError("nonfinite characteristic background")
        budget.check()
        return {**{key: value.reshape(shape) for key, value in result.items()},
                "A0": A.reshape(shape)}

    def integrals(self, T, r_faces, *, budget):
        """Endpoint identities for Ip/Iq and chi faces, without quadrature.

        These exact algebraic primitives use floating-point subtraction. They
        do not promise cancellation-free integrals in arbitrarily tiny cells
        on a nonquiet profile; the background first-jet centre branch does.
        """
        budget = _budget(budget)
        time = _array(T, "T")
        r = _array(r_faces, "r_faces")
        if (time.shape or r.ndim != 1 or len(r) < 2 or r[0] != 0
                or np.any(np.diff(r) <= 0) or r[-1] >= 1):
            raise ValueError("integrals need scalar T and increasing faces starting at zero below1")
        distance = np.arctanh(r)
        jets = self._evaluate(np.stack((float(time)-distance, float(time)+distance)), budget)
        plus, minus = jets["F"][1], jets["F"][0]
        chi = plus-minus
        budget.check(len(r)-1)
        return {"Ip": np.diff(plus+minus), "Iq": np.diff(chi), "chi_faces": chi,
                "r_faces": r.copy()}

    def radial_cuts(self, T, r_max, *, budget):
        """Map all interior profile knots from both null labels to radial splits."""
        budget = _budget(budget)
        time, radius = _array(T, "T"), _array(r_max, "r_max")
        if time.shape or radius.shape or not 0 < float(radius) < 1:
            raise ValueError("radial cuts need scalar T and0<r_max<1")
        time, radius = float(time), float(radius)
        distance = np.arctanh(radius)
        self._domain(np.array([time-distance, time+distance]))
        budget.check(2*len(self._labels)+2)
        left = np.tanh(time-self._labels)
        right = np.tanh(self._labels-time)
        mapped = np.r_[left[(left > 0) & (left < radius)],
                       right[(right > 0) & (right < radius)]]
        return np.unique(np.r_[0., mapped, radius])
