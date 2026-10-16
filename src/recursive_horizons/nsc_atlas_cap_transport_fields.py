"""Read-only field jets for the cap's full branch-selected transport direction.

This is one input to the shadow-radial mass comparison. It neither integrates
that comparison nor changes the retained metric, reconstruction or evolution.
"""
import copy

import numpy as np

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_cap_mass_tangent import _Jet, _reconstruction, _total, _departure, _cuts
from .nsc_atlas_cap_mass_derivatives import potential_partial, background_partials
from .nsc_atlas_cap_radial_reference import reference_potential
from .nsc_atlas_characteristic_background import CharacteristicProfile, MAX_POINTS, SMALL_RADIUS
from .nsc_atlas_polar_cap import _scalar, _vector
from .nsc_atlas_reduced_cap import _inputs, _reconstruct


def _switch_trace(profile, T, rec, cell, z, side, budget):
    """Named limits at the fixed small-radius switch, without a displaced query."""
    if side == "left":
        base = background_partials(profile,float(T.v),z.v,budget=budget)
    else:
        r = z.v; A=1-r*r; d=np.arctanh(r)
        jets = profile._evaluate(float(T.v)+np.array([-1.,1.])[:,None]*d,budget)
        fm,fp = jets['F']; pm,pp = jets['F_prime']; sm,sp = jets['F_second']
        f=(fp-fm)/r; Pi=(pp-pm)/(r*A); q=(pp+pm)/A; Phi=(q-f)/r
        base=dict(f=f,Pi=Pi,Phi=Phi,f_T=(pp-pm)/r,f_r=Phi,
            Pi_T=(sp-sm)/(r*A),Pi_r=(sp+sm)/(r*A*A)-Pi*(1-3*r*r)/(r*A),
            Phi_T=(sp+sm)/(A*r)-(pp-pm)/(r*r),
            Phi_r=(sp-sm)/(r*A*A)+2*q/A-2*Phi/r)
        budget.check(10*r.size)
    departure = _departure(rec,cell,z)
    budget.check(3*z.v.size*T.d.size)
    return {key:_Jet(base[key],base[key+'_T'][:,None]*T.d+base[key+'_r'][:,None]*z.d)+departure[key]
            for key in ('f','Pi','Phi')}


class CapTransportFields:
    """Frozen reconstruction with columns (fixed-r full x direction, radial).

    The first column chooses the same backward limiter branch as the mass
    tangent. The second is zero in all reconstruction inputs and is one only
    in query radii. At a cell face the caller names the owner explicitly.
    Owner-cell partials are algebraic trace jets; the material derivative is
    X_x + v*X_r, not X_x alone. A separate interface inventory is required for
    integrating them. Native samples exclude axis jets and exact switches.
    Reference samples provide named small-radius limits, and a separate
    method supplies regular axis coefficients without asserting axis jets.

    Construction charges the reused two-column reconstruction plus 12*n+2
    primal/domain checks. Each query charges all reused background/potential
    work plus 15*N returned field/coefficient jet values. All calls share the
    caller's existing metric budget and accept at most 513 query points.
    """
    def __init__(self, metric, direction, *, budget):
        self._budget = _budget(budget, metric)
        if set(direction) != {"R", "T", "Ip", "Iq"}:
            raise ValueError("transport direction requires R,T,Ip,Iq")
        y, radius, state, alpha = _inputs(metric["y_faces"], metric["radius"], metric["state"], metric["alpha"])
        n = len(y)-1
        if radius >= 1 or not isinstance(metric["profile"], CharacteristicProfile):
            raise ValueError("transport fields require a characteristic cap below radius one")
        rates = {k: _vector(direction[k], k+" rate") for k in ("Ip", "Iq")}
        if any(v.shape != (n,) for v in rates.values()):
            raise ValueError("transport rate cell counts differ")
        rx, tx = (_scalar(direction[k], k+" rate") for k in ("R", "T"))
        self._R = _Jet(radius, np.array([rx, 0.]))
        self._T = _Jet(_scalar(metric["T"], "T"), np.array([tx, 0.]))
        self._alpha = alpha
        self._profile = copy.deepcopy(metric["profile"])
        self._y_faces = y.copy()
        self._faces = radius*y
        if not np.array_equal(self._faces, metric["r_faces"]):
            raise ValueError("metric faces differ from the declared radius and grid")
        labels = self._T.v+np.array([-1., 1.])*np.arctanh(radius)
        self._profile._domain(labels)
        label_rates = tx+np.array([-1., 1.])*rx/(1-radius*radius)
        if np.any((labels == self._profile.domain["last_label"]) & (label_rates < 0)):
            raise ValueError("backward transport direction leaves upper profile domain")
        if (not self._profile.quiet_extension and
                np.any((labels == self._profile.domain["first_label"]) & (label_rates > 0))):
            raise ValueError("backward transport direction leaves lower profile domain")
        ip = _Jet(state["Ip"], np.column_stack((rates["Ip"], np.zeros(n))))
        iq = _Jet(state["Iq"], np.column_stack((rates["Iq"], np.zeros(n))))
        try:
            with np.errstate(over="raise", divide="raise", invalid="raise"):
                self._rec, branches = _reconstruction(self._R*y, ip, iq, np.array([1.,0.]), budget)
        except FloatingPointError as error:
            raise ValueError("nonfinite transport reconstruction") from error
        primal = _reconstruct(self._faces, state)
        for key, jet in self._rec.items():
            if (not np.all(np.isfinite(jet.d)) or not np.all(np.isfinite(jet.v))
                    or not np.allclose(jet.v, primal[key], rtol=2e-12, atol=2e-13)):
                raise ValueError("transport reconstruction differs from retained primal: "+key)
        self.limiter_branches = {k: tuple(v) for k,v in branches.items()}
        budget.check(12*n+2)

    def native_panels(self, *, budget):
        """Frozen native partition, including backward-born panels and owners."""
        if budget is not self._budget:
            raise ValueError("transport queries require the same caller budget")
        R,T = _Jet(self._R.v,np.array([1.,0.])),_Jet(self._T.v,np.array([0.,1.]))
        seed = np.array([self._R.d[0],self._T.d[0]])
        faces = R*self._y_faces
        events = _cuts(self._profile,T,R,faces,seed,budget)
        panels=[]
        for (left,ltag),(right,rtag) in zip(events[:-1],events[1:]):
            mid=(left+right)/2
            key=(float(mid.v),-float(mid.d@seed))
            owner=sum((float(faces.v[i]),-float(faces.d[i]@seed)) <= key
                      for i in range(len(faces.v)))-1
            owner=min(max(owner,0),len(faces.v)-2)
            panels.append(dict(left=float(left.v),right=float(right.v),
                left_rate=float(left.d@seed),right_rate=float(right.d@seed),
                owner=owner,left_tag=ltag,right_tag=rtag,born=bool(left.v==right.v)))
        budget.check(4*len(panels))
        return panels

    def sample(self, cell, radii, *, budget):
        """Return f/Pi/Phi, a/b, their x/r jets, Pdiff and its actual derivative."""
        return self._sample(cell,radii,budget=budget,reference=False,small_radius_side=None)

    def reference_sample(self, cell, radii, *, budget, small_radius_side=None):
        """Smooth-potential coefficients with named limits at r=0.01.

        Potential thresholds have no reference interface. At the small-radius
        switch the caller must name 'left' (series) or 'right' (ordinary);
        the query must then consist only of that exact radius. Cell ownership
        is still explicit and independent of this background formula side.
        """
        return self._sample(cell,radii,budget=budget,reference=True,small_radius_side=small_radius_side)

    def reference_axis_coefficients(self, *, budget):
        """Regular coefficient limits at the fixed origin; no axis field jets.

        a=O(r), b=O(r^2), a_r(0)=alpha*Pi(0)^2/2. These six limits are
        sufficient for the affine panel at its fixed centre endpoint.
        """
        if budget is not self._budget:
            raise ValueError("transport queries require the same caller budget")
        jet = self._profile.jet(float(self._T.v),budget=budget)
        pi = float(2*jet['F_second']+self._rec['P0'].v)
        ar = self._alpha*pi*pi/2
        if not np.isfinite(ar):
            raise ValueError("nonfinite reference axis coefficient")
        budget.check(8)
        return dict(a=0.,b=0.,a_x=0.,b_x=0.,a_r=ar,b_r=0.)

    def _sample(self, cell, radii, *, budget, reference, small_radius_side):
        if budget is not self._budget:
            raise ValueError("transport queries require the same caller budget")
        budget.check()
        if (isinstance(cell, (bool, np.bool_)) or not isinstance(cell, (int, np.integer))
                or not 0 <= cell < len(self._faces)-1):
            raise ValueError("transport query requires an explicit owner cell")
        r = _vector(radii, "query radii")
        if (not 1 <= len(r) <= MAX_POINTS or np.any(r <= 0)
                or np.any(r < self._faces[cell]) or np.any(r > self._faces[cell+1])):
            raise ValueError("query at most 513 positive radii in the owner cell")
        if small_radius_side is not None and (
                small_radius_side not in ('left','right') or np.any(r != SMALL_RADIUS)):
            raise ValueError("named small-radius side requires only exact switch radii")
        labels = self._T.v+np.array([-1.,1.])[:,None]*np.arctanh(r)
        self._profile._domain(labels)
        tx = float(self._T.d[0])
        if np.any(labels == self._profile.domain["last_label"]) and tx < 0:
            raise ValueError("fixed-radius backward rate leaves upper profile domain")
        if (not self._profile.quiet_extension and
                np.any(labels == self._profile.domain["first_label"]) and tx > 0):
            raise ValueError("fixed-radius backward rate leaves lower profile domain")
        z = _Jet(r, np.column_stack((np.zeros(len(r)), np.ones(len(r)))))
        try:
            with np.errstate(over="raise", divide="raise", invalid="raise"):
                fields = (_total(self._profile,self._T,self._rec,int(cell),z,budget)
                    if small_radius_side is None else
                    _switch_trace(self._profile,self._T,self._rec,int(cell),z,small_radius_side,budget))
                f = fields["f"]
                if not reference and np.any(abs(np.sqrt(self._alpha/6)*f.v) == 1e-4):
                    raise ValueError("transport query lies on the potential formula switch")
                potential = (reference_potential if reference else potential_partial)(f.v,self._alpha,budget=budget)
                pdiff = _Jet(potential["Pdiff"], potential["Pdiff_f"][:,None]*f.d)
                energy = fields["Pi"]**2+fields["Phi"]**2
                fields["a"] = self._alpha*z*energy/2
                fields["b"] = z*z*(pdiff/2+(1-z*z)*energy/4)
        except FloatingPointError as error:
            raise ValueError("nonfinite transport field jet") from error
        result = {"r": r.copy(), **potential}
        for key, jet in fields.items():
            result.update({key: jet.v.copy(), key+"_x": jet.d[:,0].copy(), key+"_r": jet.d[:,1].copy()})
        if not all(np.all(np.isfinite(v)) for v in result.values()):
            raise ValueError("nonfinite transport field output")
        budget.check(15*len(r))
        return result
