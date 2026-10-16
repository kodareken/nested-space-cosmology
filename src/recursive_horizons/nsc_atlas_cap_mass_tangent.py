"""One-sided derivative of the retained cap's actual mass-edge recurrence.

The physical direction selects one common branch before component weights.
This companion does not modify or evolve the retained metric/field producer.
"""
import numpy as np

from .nsc_atlas_background_cap import _budget, _gauss_rule
from .nsc_atlas_cap_mass_derivatives import background_partials, potential_partial
from .nsc_atlas_characteristic_background import SMALL_RADIUS
from .nsc_atlas_reduced_cap import _inputs, _fields
from .nsc_atlas_polar_cap import _scalar, _vector


class _Jet:
    """Local real arithmetic with a trailing, bounded component axis."""
    __array_priority__ = 1000

    def __init__(self, value, derivative):
        self.v = np.asarray(value, dtype=float)
        self.d = np.asarray(derivative, dtype=float)
        if self.d.shape[:-1] != self.v.shape:
            raise ValueError("tangent component shape differs from its value")

    def constant(self, value):
        value = np.asarray(value, dtype=float)
        return _Jet(value, np.zeros(value.shape+(self.d.shape[-1],)))

    def coerce(self, other):
        return other if isinstance(other, _Jet) else self.constant(other)

    def __add__(self, other):
        other = self.coerce(other)
        return _Jet(self.v+other.v, self.d+other.d)

    __radd__ = __add__

    def __neg__(self):
        return _Jet(-self.v, -self.d)

    def __sub__(self, other):
        return self+-self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other)+-self

    def __mul__(self, other):
        other = self.coerce(other)
        return _Jet(self.v*other.v,
                    self.d*other.v[..., None]+other.d*self.v[..., None])

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        value = self.v/other.v
        return _Jet(value, (self.d-value[..., None]*other.d)/other.v[..., None])

    def __rtruediv__(self, other):
        return self.coerce(other)/self

    def __pow__(self, power):
        return _Jet(self.v**power, (power*self.v**(power-1))[..., None]*self.d)

    def __getitem__(self, key):
        return _Jet(self.v[key], self.d[key])

    def exp(self):
        value = np.exp(self.v)
        return _Jet(value, value[..., None]*self.d)

    def sum(self):
        return _Jet(np.sum(self.v), np.sum(self.d, axis=tuple(range(self.v.ndim))))


def _stack(values):
    return _Jet(np.stack([v.v for v in values]), np.stack([v.d for v in values]))


def _poly(coefficients, x):
    value = x.constant(coefficients[-1])
    for coefficient in coefficients[-2::-1]:
        value = value*x+coefficient
    return value


def limiter_selection(left, right, left_rate, right_rate):
    """Winner on delta<0 for secants sharing the positive R**-2 factor.

    Rates are derivatives with respect to delta, not negated backward rates.
    Exact ties with the same rate may use either identical directional value;
    the selected component representative is recorded, not called a Jacobian.
    """
    def sign(value, rate):
        return np.sign(value) if value != 0 else np.sign(-rate)
    sl, sr = sign(left, left_rate), sign(right, right_rate)
    if sl*sr <= 0:
        return "zero"
    if left != 0 and right != 0 and left*right == 0:
        raise ValueError("underflowing limiter product has no resolved branch")
    kl = (abs(left), -sl*left_rate)
    kr = (abs(right), -sr*right_rate)
    return "left" if kl <= kr else "right"


def _slopes(midpoint, values, seed):
    secant = (values[1:]-values[:-1])/(midpoint[1:]-midpoint[:-1])
    rows, labels = [secant[0]], ["endpoint"]
    for i in range(1, len(values.v)-1):
        left, right = secant[i-1], secant[i]
        choice = limiter_selection(float(left.v), float(right.v),
                                   float(left.d@seed), float(right.d@seed))
        rows.append(left if choice == "left" else right if choice == "right" else left.constant(0.))
        labels.append(choice)
    rows.append(secant[-1]); labels.append("endpoint")
    return _stack(rows), labels


def _reconstruction(r, Ip, Iq, seed, budget):
    dr, mid = r[1:]-r[:-1], (r[1:]+r[:-1])/2
    p, q = Ip/dr, Iq/dr
    plus, minus = p+q, p-q
    sp, lp = _slopes(mid, plus, seed)
    sm, lm = _slopes(mid, minus, seed)
    ps, qs = (sp+sm)/2, (sp-sm)/2
    pc, qc = p-ps*mid, q-qs*mid
    chi = _Jet(np.r_[0., np.cumsum(Iq.v)],
               np.vstack((np.zeros(Iq.d.shape[-1]), np.cumsum(Iq.d, axis=0))))
    qm = (r[:-1]**2+r[:-1]*r[1:]+r[1:]**2)/3
    q2 = (q[1]-q[0])/(qm[1]-qm[0]); q0 = q[0]-q2*qm[0]
    pm = (r[:-1]**2+r[1:]**2)/2
    pi = p/mid
    p2 = (pi[1]-pi[0])/(pm[1]-pm[0]); p0 = pi[0]-p2*pm[0]
    budget.check(16*len(p.v)*len(seed))
    return dict(r=r, midpoint=mid, pc=pc, ps=ps, qc=qc, qs=qs, chi=chi,
                P0=p0, P2=p2, Q0=q0, Q2=q2), dict(plus=lp, minus=lm)


def _departure(rec, cell, z):
    if cell == 0:
        return dict(f=rec["Q0"]+rec["Q2"]*z*z/3,
                    Pi=rec["P0"]+rec["P2"]*z*z, Phi=2*rec["Q2"]*z/3)
    a = rec["r"][cell]
    pc, ps, qc, qs = (rec[name][cell] for name in ("pc", "ps", "qc", "qs"))
    chi = rec["chi"][cell]+(z-a)*(qc+qs*(z+a)/2)
    k = rec["chi"][cell]-qc*a-qs*a*a/2
    return dict(f=chi/z, Pi=(pc+ps*z)/z, Phi=qs/2-k/(z*z))


def _total(profile, T, rec, cell, z, budget):
    if np.any(z.v == SMALL_RADIUS):
        raise ValueError("mass sample lies on the small-radius formula switch")
    base = background_partials(profile, float(T.v), z.v, budget=budget)
    departure = _departure(rec, cell, z)
    budget.check(3*z.v.size*T.d.size)
    return {key: _Jet(base[key], base[key+"_T"][..., None]*T.d
                     +base[key+"_r"][..., None]*z.d)+departure[key]
            for key in ("f", "Pi", "Phi")}


def _integrands(profile, T, rec, cell, z, alpha, budget):
    fields = _total(profile, T, rec, cell, z, budget)
    if np.any(abs(np.sqrt(alpha/6)*fields["f"].v) == 1e-4):
        raise ValueError("mass sample lies on the potential formula switch")
    W = fields["Pi"]**2+fields["Phi"]**2
    potential = potential_partial(fields["f"].v, alpha, budget=budget)
    pdiff = _Jet(potential["Pdiff"], potential["Pdiff_f"][..., None]*fields["f"].d)
    budget.check(2*z.v.size*T.d.size)
    return alpha*z*W/2, z*z*(pdiff/2+(1-z*z)*W/4)


def _cuts(profile, T, R, faces, seed, budget):
    """Keep the preimages of exact coincidences and order their backward ray."""
    labels = profile.breakpoints
    end_labels = float(T.v)+np.array([-1., 1.])*np.arctanh(float(R.v))
    profile._domain(end_labels)
    rates = seed[1]+np.array([-1., 1.])*seed[0]/(1-float(R.v)**2)
    if np.any((end_labels == profile.domain["last_label"]) & (rates < 0)):
        raise ValueError("backward direction leaves the upper profile domain")
    if (not profile.quiet_extension and
            np.any((end_labels == profile.domain["first_label"]) & (rates > 0))):
        raise ValueError("backward direction leaves the lower profile domain")
    events = [(faces[i], ("cell", i)) for i in range(len(faces.v))]
    for sign, name in ((1, "left"), (-1, "right")):
        values = np.tanh(sign*(float(T.v)-labels))
        for index, value in enumerate(values):
            jet = _Jet(value, sign*(1-value*value)*T.d)
            rate = float(jet.d@seed)
            inside = (0 < value < R.v or value == 0 and rate < 0
                      or value == R.v and rate > seed[0])
            if inside:
                events.append((jet, (name, index)))
    events.sort(key=lambda row: (float(row[0].v), -float(row[0].d@seed)))
    ordered = []
    for row in events:
        if ordered and row[0].v == ordered[-1][0].v:
            previous = ordered[-1][0]
            if np.array_equal(row[0].d, previous.d):
                if row[1][0] == "cell": ordered[-1] = row
                continue
            if row[0].d@seed == previous.d@seed:
                raise ValueError("coincident cuts have unresolved first-order ordering")
        ordered.append(row)
    if len(ordered) > 1154:
        raise ValueError("cap tangent exceeds the finite panel allowance")
    budget.check(len(events)*len(seed))
    return ordered


def _reverse_tape(tape, seed, budget):
    """Contract the same local gradients and retain their absolute work.

    A summation bound uses the terms before cancellation, not the final mass
    rate. The path count covers both contraction orders and this positive
    work sum; eps (twice unit roundoff) and an extra denominator cover its
    own rounding. This bounds contraction arithmetic, not model error.
    """
    count = len(seed)
    reverse = np.zeros(count); weight = 1.; work = 0.
    for panel in reversed(tape):
        reverse += weight*panel["local"]
        work += weight*float(np.sum(abs(panel["local"]*seed)))
        weight *= panel["attenuation"]
        budget.check(4*count)
    k = 4*(len(tape)+count)+32
    gamma = k*np.finfo(float).eps/(1-k*np.finfo(float).eps)
    allowance = gamma*work/(1-gamma)
    if not np.isfinite(work) or not np.isfinite(allowance):
        raise ValueError("nonfinite mass contraction work")
    return reverse, dict(absolute_work=work, roundoff_allowance=allowance,
                         path_operations=k, gamma=gamma)


def mass_tangent(metric, direction, *, budget):
    """Rebuild/check primal masses; return one-sided component weights and tangent.

    direction has R,T,Ip,Iq rates. It selects limiter/cut branches, then is
    frozen. Reverse component weights belong to that branch, not to a claim
    of unrestricted differentiability. No metric or RHS is evolved here.
    """
    _budget(budget, metric)
    if set(direction) != {"R", "T", "Ip", "Iq"}:
        raise ValueError("mass direction requires R,T,Ip,Iq")
    y, radius, state, alpha = _inputs(metric["y_faces"], metric["radius"], metric["state"], metric["alpha"])
    n = len(y)-1; count = 2+2*n
    rates = {k: _vector(direction[k], k) for k in ("Ip", "Iq")}
    if any(a.shape != (n,) for a in rates.values()):
        raise ValueError("mass direction cell counts differ")
    seed = np.r_[_scalar(direction["R"], "R rate"), _scalar(direction["T"], "T rate"), rates["Ip"], rates["Iq"]]
    unit = np.eye(count)
    R, T = _Jet(radius, unit[0]), _Jet(_scalar(metric["T"], "T"), unit[1])
    Ip, Iq = _Jet(state["Ip"], unit[2:2+n]), _Jet(state["Iq"], unit[2+n:])
    faces = R*y
    rec, branches = _reconstruction(faces, Ip, Iq, seed, budget)
    events = _cuts(metric["profile"], T, R, faces, seed, budget)
    nodes, weights, basis = _gauss_rule(metric["quadrature_order"])
    N, tangent = 0., 0.
    face_mass = np.zeros(n+1)
    tape = []
    for (left, ltag), (right, rtag) in zip(events[:-1], events[1:]):
        mid, h = (left+right)/2, (right-left)/2
        key = (float(mid.v), -float(mid.d@seed))
        cell = sum((float(faces.v[i]), -float(faces.d[i]@seed)) <= key
                   for i in range(n+1))-1
        cell = min(max(cell, 0), n-1)
        born = bool(h.v == 0)
        if born:
            if mid.v == 0:
                local = R.constant(0.)
            else:
                a, b = _integrands(metric["profile"], T, rec, cell, mid, alpha, budget)
                local = (right-left)*(b-a*N)
            attenuation, next_N = 1., N
        else:
            if h.v < 0:
                raise ValueError("mass tangent cut ordering reversed")
            primary = mid+h*nodes
            a, _ = _integrands(metric["profile"], T, rec, cell, primary, alpha, budget)
            edge = left+h*(1+nodes)
            ap, bp = _integrands(metric["profile"], T, rec, cell, edge, alpha, budget)
            theta = (edge-mid)/h
            Gpoly = sum(h*a[k]*_poly(poly.coef, theta) for k, poly in enumerate(basis))
            G = h*(ap*weights).sum()
            attenuation = float(np.exp(-G.v))
            local = (-G).exp()*N+h*((Gpoly-G).exp()*bp*weights).sum()
            next_N = float(local.v)
        tangent = attenuation*tangent+float(local.d@seed)
        tape.append(dict(attenuation=attenuation, local=local.d,
                         left=float(left.v), right=float(right.v), owner=cell, born=born))
        N = next_N
        if rtag[0] == "cell": face_mass[rtag[1]] = N
        budget.check(count)
    reverse, contraction = _reverse_tape(tape, seed, budget)
    if not np.all(np.isfinite(reverse)) or not np.isfinite(tangent):
        raise ValueError("nonfinite mass tangent")
    defects = abs(face_mass-np.asarray(metric["n_faces"]))
    if not np.allclose(face_mass, metric["n_faces"], rtol=2e-12, atol=2e-13):
        raise ValueError("mass edge tape differs from the retained primal metric")
    return dict(N=N, N_R=float(reverse[0]), N_T=float(reverse[1]),
                Ip_weights=reverse[2:2+n], Iq_weights=reverse[2+n:],
                forward=float(tangent), reverse=float(reverse@seed),
                **{"contraction_"+k:v for k,v in contraction.items()},
                primal_face_defect=float(np.max(defects)),
                primal_faces=face_mass, limiter_branches=branches,
                panels=[{k: v for k, v in panel.items() if k != "local"} for panel in tape],
                panel_count=len(tape), born_panels=sum(row["born"] for row in tape))


def mass_balance_ledger(metric, rhs, tangent, *, R_x, T_x, budget):
    """Telescope the unchanged departure RHS, with explicit upwind work."""
    _budget(budget, metric)
    R_x, T_x = _scalar(R_x, "R_x"), _scalar(T_x, "T_x", positive=True)
    wp, wq = (_vector(tangent[key], key) for key in ("Ip_weights", "Iq_weights"))
    n = len(metric["y_faces"])-1
    if wp.shape != (n,) or wq.shape != (n,):
        raise ValueError("mass weights do not cover the departure cells")
    fp, fq = (_vector(rhs[key], key) for key in ("p_flux", "q_flux"))
    reaction = _vector(rhs["reaction_integral"], "reaction_integral")
    if fp.shape != (n+1,) or fq.shape != (n+1,) or reaction.shape != (n,):
        raise ValueError("mass balance flux shapes differ")
    radius = metric["radius"]
    out = dict(radius=tangent["N_R"]*R_x, clock=tangent["N_T"]*T_x,
               centre=-wp[0]*fp[0]-wq[0]*fq[0], wall=wp[-1]*fp[-1]+wq[-1]*fq[-1],
               interior=float(np.diff(-wp)@fp[1:-1]+np.diff(-wq)@fq[1:-1]),
               reaction=float(-T_x*(wp@reaction)))
    central, upwind = 0., 0.
    face_defect = 0.
    for j in range(1, n):
        r, y, B = metric["r_faces"][j], metric["y_faces"][j], metric["B_faces"][j]
        left = _fields(metric["reconstruction"], j-1, np.array([r]))
        right = _fields(metric["reconstruction"], j, np.array([r]))
        pl, ql, pr, qr = (float(a[key][0]) for a, key in ((left,"p"),(left,"q"),(right,"p"),(right,"q")))
        dB = ((1-r*r)*np.expm1(metric["delta_faces"][j])
              -2*metric["alpha"]*metric["n_faces"][j]*np.exp(metric["delta_faces"][j])/r)
        cp = T_x*(B*(ql+qr)/2+dB*metric["background_faces"]["q"][j])+R_x*y*(pl+pr)/2
        cq = T_x*(B*(pl+pr)/2+dB*metric["background_faces"]["p"][j])+R_x*y*(ql+qr)/2
        sp, sm = -T_x*B-R_x*y, T_x*B-R_x*y
        jp, jm = pr+qr-pl-ql, pr-qr-pl+ql
        up, uq = (abs(sp)*jp+abs(sm)*jm)/4, (abs(sp)*jp-abs(sm)*jm)/4
        dw, dq = wp[j-1]-wp[j], wq[j-1]-wq[j]
        central += dw*cp+dq*cq
        upwind += dw*up+dq*uq
        face_defect = max(face_defect, abs(fp[j]-cp-up), abs(fq[j]-cq-uq))
    total = sum(out.values())
    direct = out["radius"]+out["clock"]+wp@np.asarray(rhs["rates"]["Ip"])+wq@np.asarray(rhs["rates"]["Iq"])
    out.update(total=float(total), direct=float(direct), telescoping_defect=float(total-direct),
               interior_central=float(central), interior_upwind=float(upwind),
               interior_split_defect=float(out["interior"]-central-upwind),
               maximum_face_split_defect=float(face_defect),
               incoming_mass_rate_derivative=float((wp[-1]+wq[-1])*radius*(T_x*metric["B_faces"][-1]+R_x)/2))
    budget.check(16*n)
    if not all(np.isfinite(value) for value in out.values()):
        raise ValueError("nonfinite cap mass balance ledger")
    return out
