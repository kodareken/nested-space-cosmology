"""Stateful Gaussian boundary elimination in a fixed canonical spin frame.

Units hbar=1. The supplied child histories are unitaries from the same t0;
links have energy units. No Markov, factorized-state or vacuum assumption is
made. This finite quadratic owner does not supply continuum vacuum matching.
"""
from dataclasses import dataclass

import numpy as np

from .nsc_influence import _covariance


@dataclass(frozen=True)
class GaussianBoundaryState:
    covariance: np.ndarray
    parent_mask: np.ndarray

    def __post_init__(self):
        c = np.array(_covariance(self.covariance), copy=True)
        p = np.asarray(self.parent_mask)
        if p.shape != (len(c),) or not np.isin(p, [False, True]).all():
            raise ValueError("matching binary parent partition required")
        p = np.array(p, dtype=bool, copy=True)
        if not 0 < p.sum() < len(p):
            raise ValueError("both parent and child must be present")
        c.setflags(write=False)
        p.setflags(write=False)
        object.__setattr__(self, "covariance", c)
        object.__setattr__(self, "parent_mask", p)

    @property
    def blocks(self):
        p, c = self.parent_mask, self.covariance
        return (c[np.ix_(p, p)], c[np.ix_(p, ~p)],
                c[np.ix_(~p, p)], c[np.ix_(~p, ~p)])

    def _source_map(self, link, child_history):
        np_, nc = int(self.parent_mask.sum()), int((~self.parent_mask).sum())
        b, u = np.asarray(link), np.asarray(child_history)
        if b.shape != (np_, nc) or u.shape != (nc, nc):
            raise ValueError("link/history dimensions must match the partition")
        if not np.isfinite(b).all() or not np.isfinite(u).all():
            raise ValueError("finite link and child history required")
        if not np.allclose(u.conj().T @ u, np.eye(nc), rtol=0, atol=1e-10):
            raise ValueError("unitary child history required")
        return b @ u

    def kernels(self, t, s, *, link_t, child_t, link_s, child_s):
        """Explicit two-time bath and initial-surface kernels.

        F(t)=B(t) Uc(t,t0), eta(t)=F(t)c0. The occupied kernel is
        <eta_j†(s) eta_i(t)>; the empty kernel is <eta_i(t) eta_j†(s)>.
        Sigma< = i occupied, Sigma> = -i empty; theta(0)=1/2.
        Preparation matrices are <p_j†(t0) eta_i(t)> and its empty
        ordering. They have energy units; two-time kernels have energy².
        """
        if not np.isfinite(t) or not np.isfinite(s):
            raise ValueError("finite times in a common clock required")
        ft = self._source_map(link_t, child_t)
        fs = self._source_map(link_s, child_s)
        _, _, cp, cc = self.blocks
        occupied = ft @ cc @ fs.conj().T
        empty = ft @ (np.eye(len(cc)) - cc) @ fs.conj().T
        theta = 1. if t > s else (.5 if t == s else 0.)
        initial = ft @ cp
        return {
            "occupied": occupied,
            "empty": empty,
            "retarded": -1j * theta * (ft @ fs.conj().T),
            "lesser": 1j * occupied,
            "greater": -1j * empty,
            "keldysh": -1j * (empty - occupied),
            "initial_occupied": initial,
            "initial_empty": -initial,
        }

    def reconstruct(self, a_t, v_t, a_s=None, v_s=None):
        """Return all four terms of <p†(s)p(t)>.

        p(t)=A(t)p0+V(t)c0. A and V must come from the retarded
        equation, not a Markov approximation. Same-time is the default.
        """
        if (a_s is None) != (v_s is None):
            raise ValueError("supply both histories at s, or neither")
        if a_s is None:
            a_s, v_s = a_t, v_t
        pp, pc, cp, cc = self.blocks
        a, v, aa, vv = map(np.asarray, (a_t, v_t, a_s, v_s))
        if a.shape != pp.shape or aa.shape != pp.shape:
            raise ValueError("parent response has wrong dimensions")
        if v.shape != pc.shape or vv.shape != pc.shape:
            raise ValueError("child response has wrong dimensions")
        terms = {
            "parent": a @ pp @ aa.conj().T,
            "child": v @ cc @ vv.conj().T,
            "initial_parent_child": a @ pc @ vv.conj().T,
            "initial_child_parent": v @ cp @ aa.conj().T,
        }
        terms["total"] = sum(terms.values())
        return terms

    def overlap_schur(self, plus, minus):
        """Exact correlated Gaussian boundary action at the history endpoints.

        Q=I-C+C Uminus†Uplus already includes preparation and normalization.
        det Q=det Qcc det(Qpp-Qpc Qcc^-1 Qcp). The factors are NOT
        independent normalized physical influences. No inverse of C or I-C
        is used, so pure initial covariances are supported directly.
        An invertible Qcc is required for this particular block pivot.
        """
        c, p = self.covariance, self.parent_mask
        identity = np.eye(len(c))
        plus, minus = np.asarray(plus), np.asarray(minus)
        for u in (plus, minus):
            if u.shape != c.shape or not np.isfinite(u).all():
                raise ValueError("finite full-system history required")
            if not np.allclose(u.conj().T @ u, identity, rtol=0, atol=1e-10):
                raise ValueError("unitary full-system history required")
        q = identity - c + c @ minus.conj().T @ plus
        pp, pc = q[np.ix_(p, p)], q[np.ix_(p, ~p)]
        cp, cc = q[np.ix_(~p, p)], q[np.ix_(~p, ~p)]
        reduced = pp - pc @ np.linalg.solve(cc, cp)
        sc, lc = np.linalg.slogdet(cc)
        sp, lp = np.linalg.slogdet(reduced)
        phase = float(np.angle(sc * sp)) if sc != 0 and sp != 0 else None
        log_modulus = float(lc + lp)
        return {
            "child_endpoint_kernel": cc,
            "parent_endpoint_kernel": reduced,
            "amplitude": complex(sc * sp * np.exp(lc + lp)),
            "log_modulus": log_modulus,
            "principal_action": complex(phase, -log_modulus) if phase is not None else None,
        }
