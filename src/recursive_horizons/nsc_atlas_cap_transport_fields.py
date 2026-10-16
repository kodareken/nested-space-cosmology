"""Read-only field jets for the cap's full branch-selected transport direction.

This is one input to the shadow-radial mass comparison. It neither integrates
that comparison nor changes the retained metric, reconstruction or evolution.
"""
import copy

import numpy as np

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_cap_mass_tangent import _Jet, _reconstruction, _total
from .nsc_atlas_cap_mass_derivatives import potential_partial
from .nsc_atlas_characteristic_background import CharacteristicProfile, MAX_POINTS
from .nsc_atlas_polar_cap import _scalar, _vector
from .nsc_atlas_reduced_cap import _inputs, _reconstruct


class CapTransportFields:
    """Frozen reconstruction with columns (fixed-r full x direction, radial).

    The first column chooses the same backward limiter branch as the mass
    tangent. The second is zero in all reconstruction inputs and is one only
    in query radii. At a cell face the caller names the owner explicitly.
    Owner-cell partials are algebraic trace jets; the material derivative is
    X_x + v*X_r, not X_x alone. A separate interface inventory is required for
    integrating them. Axis jets and exact implementation switches are excluded.

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

    def sample(self, cell, radii, *, budget):
        """Return f/Pi/Phi, a/b, their x/r jets, Pdiff and its actual derivative."""
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
                fields = _total(self._profile, self._T, self._rec, int(cell), z, budget)
                f = fields["f"]
                if np.any(abs(np.sqrt(self._alpha/6)*f.v) == 1e-4):
                    raise ValueError("transport query lies on the potential formula switch")
                potential = potential_partial(f.v, self._alpha, budget=budget)
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
