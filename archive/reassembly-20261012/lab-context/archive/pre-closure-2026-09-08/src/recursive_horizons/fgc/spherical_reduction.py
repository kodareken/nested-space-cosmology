"""RED1 exact local jets for the unredefined spherical FGC equations.

This is deliberately a *pointwise reduction certificate*, not an evolution
system.  It compares a direct four-dimensional Christoffel/Riemann evaluator
at the regular angular point ``theta=pi/2`` with an independent warped-product
construction.  All inputs and outputs are exact :class:`fractions.Fraction`.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction
from hashlib import sha256
import json
from typing import Any, Mapping

from .exact_interval import Interval
from .exact_tangent import FirstTangent
from .interval_tangent import IntervalFirstTangent

Q = Fraction
N = 4

BASE_FIELD_ORDER = (
    "h_tt",
    "h_tr",
    "h_rr",
    "areal_radius",
    "phi",
    "chi",
)
ADM_FIELD_ORDER = (
    "alpha",
    "shift",
    "lambda",
    "areal_radius",
    "phi",
    "chi",
)
SECOND_DERIVATIVE_ORDER = ("dtt", "dtr", "drr")
INDEPENDENT_EQUATION_ORDER = (
    "metric_tt",
    "metric_tr",
    "metric_rr",
    "metric_theta_theta",
    "scalar_phi",
    "scalar_chi",
)


def _q(
    x: int | Fraction | FirstTangent | Interval | IntervalFirstTangent,
) -> Fraction | FirstTangent | Interval | IntervalFirstTangent:
    if isinstance(x, bool):
        raise ValueError("booleans are not rational jet values")
    if isinstance(x, (FirstTangent, Interval, IntervalFirstTangent)):
        return x
    return x if isinstance(x, Fraction) else Fraction(x)


def _interval_primal(value: object) -> Interval | None:
    if isinstance(value, IntervalFirstTangent):
        return value.primal
    if isinstance(value, Interval):
        return value
    return None


def _whole_strictly_positive(value: object) -> bool:
    primal = _interval_primal(value)
    return primal.strictly_positive() if primal is not None else value > 0


def _whole_strictly_negative(value: object) -> bool:
    primal = _interval_primal(value)
    return primal.strictly_negative() if primal is not None else value < 0


def _contains_zero(value: object) -> bool:
    primal = _interval_primal(value)
    return primal.contains_zero() if primal is not None else value == 0


def _provably_nonzero(value: object) -> bool:
    return not _contains_zero(value)


def _fraction_text(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


@dataclass(frozen=True)
class Jet2:
    """A scalar two-jet ``(value, dt, dr, dtt, dtr, drr)`` at one point."""

    value: Fraction
    dt: Fraction = Q(0)
    dr: Fraction = Q(0)
    dtt: Fraction = Q(0)
    dtr: Fraction = Q(0)
    drr: Fraction = Q(0)

    def __post_init__(self) -> None:
        for name in ("value", "dt", "dr", "dtt", "dtr", "drr"):
            object.__setattr__(self, name, _q(getattr(self, name)))

    @staticmethod
    def constant(value: int | Fraction) -> "Jet2":
        return Jet2(_q(value))

    def __add__(self, other: object) -> "Jet2":
        b = as_jet(other)
        return Jet2(
            self.value + b.value,
            self.dt + b.dt,
            self.dr + b.dr,
            self.dtt + b.dtt,
            self.dtr + b.dtr,
            self.drr + b.drr,
        )

    __radd__ = __add__

    def __neg__(self) -> "Jet2":
        return Jet2(-self.value, -self.dt, -self.dr, -self.dtt, -self.dtr, -self.drr)

    def __sub__(self, other: object) -> "Jet2":
        return self + (-as_jet(other))

    def __rsub__(self, other: object) -> "Jet2":
        return as_jet(other) - self

    def __mul__(self, other: object) -> "Jet2":
        b = as_jet(other)
        a = self
        return Jet2(
            a.value * b.value,
            a.dt * b.value + a.value * b.dt,
            a.dr * b.value + a.value * b.dr,
            a.dtt * b.value + 2 * a.dt * b.dt + a.value * b.dtt,
            a.dtr * b.value + a.dt * b.dr + a.dr * b.dt + a.value * b.dtr,
            a.drr * b.value + 2 * a.dr * b.dr + a.value * b.drr,
        )

    __rmul__ = __mul__

    def compose(
        self, value: int | Fraction, first: int | Fraction, second: int | Fraction
    ) -> "Jet2":
        """Return the jet of ``q(self)`` given ``q,q',q''`` at ``self.value``."""
        q0, q1, q2 = _q(value), _q(first), _q(second)
        a = self
        return Jet2(
            q0,
            q1 * a.dt,
            q1 * a.dr,
            q1 * a.dtt + q2 * a.dt * a.dt,
            q1 * a.dtr + q2 * a.dt * a.dr,
            q1 * a.drr + q2 * a.dr * a.dr,
        )

    def reciprocal(self) -> "Jet2":
        if _contains_zero(self.value):
            raise ValueError("singular jet reciprocal")
        return self.compose(
            1 / self.value, -1 / (self.value * self.value), 2 / (self.value**3)
        )

    def __truediv__(self, other: object) -> "Jet2":
        return self * as_jet(other).reciprocal()

    def __rtruediv__(self, other: object) -> "Jet2":
        return as_jet(other) / self


def as_jet(x: object) -> Jet2:
    return x if isinstance(x, Jet2) else Jet2.constant(_q(x))  # type: ignore[arg-type]


def _fixture_jet(value: Jet2 | Mapping[str, Any] | int | Fraction, name: str) -> Jet2:
    """Read a public rational jet mapping, accepting omitted derivative keys."""
    if isinstance(value, Jet2):
        return value
    if isinstance(value, Mapping):
        if "value" not in value:
            raise ValueError(f"{name} jet requires value")
        allowed = {"value", "dt", "dr", "dtt", "dtr", "drr"}
        if set(value) - allowed:
            raise ValueError(f"{name} jet has unknown keys")
        return Jet2(
            *(_q(value.get(k, 0)) for k in ("value", "dt", "dr", "dtt", "dtr", "drr"))
        )
    return Jet2.constant(_q(value))


@dataclass(frozen=True)
class SphericalState:
    """General spherical two-jet, ``ds2=h_AB dxA dxB+r2 dOmega2``.

    Index order is ``t,r,theta,phi``.  The base determinant must be negative,
    ``r>0``, and the angular evaluation point is the regular equator.
    """

    h_tt: Jet2
    h_tr: Jet2
    h_rr: Jet2
    areal_radius: Jet2
    phi: Jet2
    chi: Jet2
    planck_mass: Fraction = Q(1)
    beta: Fraction = Q(0)
    mu: Fraction = Q(1)
    g4: Fraction = Q(1)
    alpha: Fraction = Q(0)
    eta: Fraction = Q(0)
    branch: str = "FGC-QR"

    def __post_init__(self) -> None:
        for n in ("planck_mass", "beta", "mu", "g4", "alpha", "eta"):
            object.__setattr__(self, n, _q(getattr(self, n)))
        if not _whole_strictly_positive(self.areal_radius.value):
            raise ValueError("areal_radius must be positive")
        if not _whole_strictly_negative(
            self.h_tt.value * self.h_rr.value - self.h_tr.value * self.h_tr.value
        ):
            raise ValueError("two-dimensional base must be Lorentzian")
        if self.planck_mass <= 0 or self.mu <= 0 or self.g4 <= 0:
            raise ValueError("positive parameters required")
        if self.branch not in {"GR-0", "SGB-L", "FGC-QR"}:
            raise ValueError("unknown branch")
        if self.branch == "GR-0" and (self.beta or self.alpha or self.eta):
            raise ValueError("GR-0 couplings must vanish")
        if self.branch == "SGB-L" and (not self.alpha or self.beta or self.eta):
            raise ValueError("SGB-L requires alpha only")
        if self.branch == "FGC-QR" and (not self.beta or self.alpha or self.eta <= 0):
            raise ValueError("FGC-QR requires beta!=0, eta>0, and alpha=0")
        if not _whole_strictly_positive(
            self.planck_mass**2 + self.beta * self.phi.value**2
        ):
            raise ValueError("effective Planck coefficient F must be positive")


def _adm_pg_jets_from_fixture(fixture: Mapping[str, Any]) -> dict[str, Jet2]:
    """Read the six exact generalized-ADM/PG field jets in frozen order."""

    raw_state = fixture.get("state", fixture)
    if not isinstance(raw_state, Mapping):
        raise ValueError("ADM/PG state must be a mapping")
    missing = set(ADM_FIELD_ORDER) - set(raw_state)
    if missing:
        raise ValueError(f"ADM/PG fixture missing keys: {sorted(missing)}")

    def jet(name: str) -> Jet2:
        value = raw_state[name]
        if isinstance(value, Mapping):
            return _fixture_jet(value, name)
        return Jet2(
            *(
                _q(raw_state.get(name + suffix, value if suffix == "" else 0))
                for suffix in ("", "_t", "_r", "_tt", "_tr", "_rr")
            )
        )

    jets = {name: jet(name) for name in ADM_FIELD_ORDER}
    if jets["alpha"].value <= 0 or jets["lambda"].value <= 0:
        raise ValueError("ADM/PG alpha and lambda must be positive")
    return jets


def _state_from_adm_pg_jets(
    fixture: Mapping[str, Any], jets: Mapping[str, Jet2]
) -> SphericalState:
    lapse = jets["alpha"]
    shift = jets["shift"]
    lam = jets["lambda"]
    lam_squared = lam * lam
    params = fixture.get("action_parameters", fixture)
    if not isinstance(params, Mapping):
        raise ValueError("action_parameters must be a mapping")
    model = str(
        fixture.get("model_id", params.get("model_id", fixture.get("branch", "FGC-QR")))
    )
    branch = {"GR-0": "GR-0", "SGB-L": "SGB-L", "FGC-QR": "FGC-QR"}.get(model, model)

    def parameter(*names: str, default: Any = 0) -> Fraction:
        for name in names:
            if name in params:
                return _q(params[name])
        return _q(default)

    return SphericalState(
        h_tt=-(lapse * lapse) + lam_squared * shift * shift,
        h_tr=lam_squared * shift,
        h_rr=lam_squared,
        areal_radius=jets["areal_radius"],
        phi=jets["phi"],
        chi=jets["chi"],
        branch=branch,
        planck_mass=parameter("planck_mass", default=1),
        beta=parameter("ricci_coupling", "beta", default=0),
        mu=parameter("scalar_mass", "mu", default=1),
        g4=parameter("quartic_coupling", "g4", default=1),
        alpha=parameter("linear_gb_coupling", "alpha_gb", default=0),
        eta=parameter("quadratic_gb_coupling", "eta", default=0),
    )


def state_from_generalized_adm_pg_fixture(
    fixture: Mapping[str, Any],
) -> SphericalState:
    """Map generalized ADM/PG jets to the exact ``2+2`` state.

    The required metric names are ``alpha``, ``shift``, and ``lambda``. Their
    exact :class:`Jet2` products define
    ``h_tt=-alpha^2+lambda^2 shift^2``, ``h_tr=lambda^2 shift``, and
    ``h_rr=lambda^2``. An explicit ``areal_radius`` is required; an areal
    coordinate fixture is ``{value=r, dr=1}`` with other derivatives zero.
    """

    return _state_from_adm_pg_jets(fixture, _adm_pg_jets_from_fixture(fixture))


def _zero4():
    return [
        [[[Q(0) for _ in range(N)] for _ in range(N)] for _ in range(N)]
        for _ in range(N)
    ]


def _zero2():
    return [[Q(0) for _ in range(N)] for _ in range(N)]


def _inv(m):
    # exact Gauss-Jordan, deliberately independent of the curvature formulae
    a = [
        [m[i][j] for j in range(N)] + [Q(int(i == j)) for j in range(N)]
        for i in range(N)
    ]
    for c in range(N):
        p = next((r for r in range(c, N) if _provably_nonzero(a[r][c])), None)
        if p is None:
            raise ValueError("singular metric")
        a[c], a[p] = a[p], a[c]
        z = a[c][c]
        a[c] = [x / z for x in a[c]]
        for r in range(N):
            if r != c:
                z = a[r][c]
                a[r] = [a[r][j] - z * a[c][j] for j in range(2 * N)]
    return [r[N:] for r in a]


def _direct_metric(state: SphericalState):
    """Return g, dg, ddg at theta=pi/2; angular derivatives are explicit."""
    jets = (state.h_tt, state.h_tr, state.h_rr, state.areal_radius)
    g = _zero2()
    dg = [_zero2() for _ in range(N)]
    ddg = [[_zero2() for _ in range(N)] for _ in range(N)]

    def put(i, j, x: Jet2, angular=False):
        g[i][j] = x.value
        dg[0][i][j] = x.dt
        dg[1][i][j] = x.dr
        ddg[0][0][i][j] = x.dtt
        ddg[0][1][i][j] = ddg[1][0][i][j] = x.dtr
        ddg[1][1][i][j] = x.drr

    put(0, 0, jets[0])
    put(0, 1, jets[1])
    put(1, 0, jets[1])
    put(1, 1, jets[2])
    rr = state.areal_radius * state.areal_radius
    put(2, 2, rr)
    put(3, 3, rr)
    # d_theta sin^2(theta)=0, d_theta^2 sin^2(theta)=-2 at the equator.
    ddg[2][2][3][3] = -2 * rr.value
    return g, dg, ddg


def direct_4d_curvature(state: SphericalState):
    """Direct Christoffel/Riemann evaluation from four-dimensional metric jets."""
    g, dg, ddg = _direct_metric(state)
    gi = _inv(g)
    gamma = _zero4()
    for a in range(N):
        for b in range(N):
            for c in range(N):
                gamma[a][b][c] = sum(
                    gi[a][d] * (dg[b][d][c] + dg[c][d][b] - dg[d][b][c]) / 2
                    for d in range(N)
                )
    dgi = [_zero2() for _ in range(N)]
    for e in range(N):
        for a in range(N):
            for b in range(N):
                dgi[e][a][b] = -sum(
                    gi[a][m] * dg[e][m][n] * gi[n][b]
                    for m in range(N)
                    for n in range(N)
                )
    dgamma = [_zero4() for _ in range(N)]
    for e in range(N):
        for a in range(N):
            for b in range(N):
                for c in range(N):
                    dgamma[e][a][b][c] = sum(
                        dgi[e][a][d] * (dg[b][d][c] + dg[c][d][b] - dg[d][b][c]) / 2
                        + gi[a][d]
                        * (ddg[e][b][d][c] + ddg[e][c][d][b] - ddg[e][d][b][c])
                        / 2
                        for d in range(N)
                    )
    rup = _zero4()
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    rup[a][b][c][d] = (
                        dgamma[c][a][d][b]
                        - dgamma[d][a][c][b]
                        + sum(
                            gamma[a][c][e] * gamma[e][d][b]
                            - gamma[a][d][e] * gamma[e][c][b]
                            for e in range(N)
                        )
                    )
    low = _zero4()
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    low[a][b][c][d] = sum(g[a][e] * rup[e][b][c][d] for e in range(N))
    return low, gi, gamma


def _ricci(ri, gi):
    return [
        [
            sum(gi[a][c] * ri[a][b][c][d] for a in range(N) for c in range(N))
            for d in range(N)
        ]
        for b in range(N)
    ]


def _scalar(ric, gi):
    return sum(gi[a][b] * ric[a][b] for a in range(N) for b in range(N))


def _gb(ri, ric, sc, gi):
    r2 = sum(
        gi[a][e] * gi[b][f] * gi[c][g] * gi[d][h] * ri[a][b][c][d] * ri[e][f][g][h]
        for a in range(N)
        for b in range(N)
        for c in range(N)
        for d in range(N)
        for e in range(N)
        for f in range(N)
        for g in range(N)
        for h in range(N)
    )
    q = sum(
        gi[a][c] * gi[b][d] * ric[a][b] * ric[c][d]
        for a in range(N)
        for b in range(N)
        for c in range(N)
        for d in range(N)
    )
    return r2 - 4 * q + sc * sc


def warped_2plus2_curvature(state: SphericalState):
    """Independent warped-product Riemann route, returned with 4D lower indices."""
    # Build base Christoffels/Riemann directly from the three base jets.
    h = [[state.h_tt.value, state.h_tr.value], [state.h_tr.value, state.h_rr.value]]
    hi = [
        [
            h[1][1] / (h[0][0] * h[1][1] - h[0][1] ** 2),
            -h[0][1] / (h[0][0] * h[1][1] - h[0][1] ** 2),
        ],
        [
            -h[0][1] / (h[0][0] * h[1][1] - h[0][1] ** 2),
            h[0][0] / (h[0][0] * h[1][1] - h[0][1] ** 2),
        ],
    ]
    j = (state.h_tt, state.h_tr, state.h_rr)
    dh = [
        [[j[0].dt, j[1].dt], [j[1].dt, j[2].dt]],
        [[j[0].dr, j[1].dr], [j[1].dr, j[2].dr]],
    ]
    dd = [
        [
            [[j[0].dtt, j[1].dtt], [j[1].dtt, j[2].dtt]],
            [[j[0].dtr, j[1].dtr], [j[1].dtr, j[2].dtr]],
        ],
        [
            [[j[0].dtr, j[1].dtr], [j[1].dtr, j[2].dtr]],
            [[j[0].drr, j[1].drr], [j[1].drr, j[2].drr]],
        ],
    ]
    gh = [
        [
            [
                sum(
                    hi[a][d] * (dh[b][d][c] + dh[c][d][b] - dh[d][b][c]) / 2
                    for d in range(2)
                )
                for c in range(2)
            ]
            for b in range(2)
        ]
        for a in range(2)
    ]
    # derivative Gamma including d inverse
    dhi = [
        [
            [
                -sum(
                    hi[a][m] * dh[e][m][n] * hi[n][b]
                    for m in range(2)
                    for n in range(2)
                )
                for b in range(2)
            ]
            for a in range(2)
        ]
        for e in range(2)
    ]
    dgh = [
        [
            [
                [
                    sum(
                        dhi[e][a][d] * (dh[b][d][c] + dh[c][d][b] - dh[d][b][c]) / 2
                        + hi[a][d]
                        * (dd[e][b][d][c] + dd[e][c][d][b] - dd[e][d][b][c])
                        / 2
                        for d in range(2)
                    )
                    for c in range(2)
                ]
                for b in range(2)
            ]
            for a in range(2)
        ]
        for e in range(2)
    ]
    rb = [
        [
            [
                [
                    dgh[c][a][d][b]
                    - dgh[d][a][c][b]
                    + sum(
                        gh[a][c][e] * gh[e][d][b] - gh[a][d][e] * gh[e][c][b]
                        for e in range(2)
                    )
                    for d in range(2)
                ]
                for c in range(2)
            ]
            for b in range(2)
        ]
        for a in range(2)
    ]
    r = state.areal_radius
    hr = [
        [
            r.dtt - sum(gh[c][0][0] * ([r.dt, r.dr][c]) for c in range(2)),
            r.dtr - sum(gh[c][0][1] * ([r.dt, r.dr][c]) for c in range(2)),
        ],
        [
            r.dtr - sum(gh[c][1][0] * ([r.dt, r.dr][c]) for c in range(2)),
            r.drr - sum(gh[c][1][1] * ([r.dt, r.dr][c]) for c in range(2)),
        ],
    ]
    grad2 = sum(
        hi[a][b] * [r.dt, r.dr][a] * [r.dt, r.dr][b] for a in range(2) for b in range(2)
    )
    out = _zero4()
    for a in range(2):
        for b in range(2):
            for c in range(2):
                for d in range(2):
                    out[a][b][c][d] = sum(h[a][e] * rb[e][b][c][d] for e in range(2))
    # R_{AiBj}=-r nabla_A nabla_B r gamma_{ij}; populate all Riemann
    # symmetry partners, not only the representative component.
    for a in range(2):
        for b in range(2):
            for i in (2, 3):
                x = -r.value * hr[a][b]
                out[a][i][b][i] = x
                out[i][a][b][i] = -x
                out[a][i][i][b] = -x
                out[i][a][i][b] = x
                out[b][i][a][i] = x
                out[i][b][a][i] = -x
                out[b][i][i][a] = -x
                out[i][b][i][a] = x
    for i in (2, 3):
        for j in (2, 3):
            for k in (2, 3):
                for l in (2, 3):
                    out[i][j][k][l] = (
                        r.value
                        * r.value
                        * (1 - grad2)
                        * (
                            (1 if i == k else 0) * (1 if j == l else 0)
                            - (1 if i == l else 0) * (1 if j == k else 0)
                        )
                    )
    return out


def _warped_2plus2_metric_connection(state: SphericalState):
    """Build the equatorial metric, inverse, and connection from 2+2 identities.

    This route is deliberately separate from ``_direct_metric``, ``_inv``, and
    ``direct_4d_curvature``.  At the equator the intrinsic unit-sphere
    Christoffels vanish, while the mixed warped-product coefficients are

    ``Gamma^A_ij=-R partial^A R gamma_ij`` and
    ``Gamma^i_Aj=(partial_A R/R) delta^i_j``.
    """

    h = [[state.h_tt.value, state.h_tr.value], [state.h_tr.value, state.h_rr.value]]
    determinant = h[0][0] * h[1][1] - h[0][1] ** 2
    hi = [
        [h[1][1] / determinant, -h[0][1] / determinant],
        [-h[0][1] / determinant, h[0][0] / determinant],
    ]
    jets = (state.h_tt, state.h_tr, state.h_rr)
    dh = [
        [[jets[0].dt, jets[1].dt], [jets[1].dt, jets[2].dt]],
        [[jets[0].dr, jets[1].dr], [jets[1].dr, jets[2].dr]],
    ]
    base_gamma = [
        [
            [
                sum(
                    hi[a][d] * (dh[b][d][c] + dh[c][d][b] - dh[d][b][c]) / 2
                    for d in range(2)
                )
                for c in range(2)
            ]
            for b in range(2)
        ]
        for a in range(2)
    ]

    g = _zero2()
    gi = _zero2()
    for a in range(2):
        for b in range(2):
            g[a][b] = h[a][b]
            gi[a][b] = hi[a][b]
    radius = state.areal_radius.value
    for i in (2, 3):
        g[i][i] = radius**2
        gi[i][i] = 1 / radius**2

    gamma = [[[Q(0) for _ in range(N)] for _ in range(N)] for _ in range(N)]
    for a in range(2):
        for b in range(2):
            for c in range(2):
                gamma[a][b][c] = base_gamma[a][b][c]
    radius_gradient = (state.areal_radius.dt, state.areal_radius.dr)
    raised_radius_gradient = tuple(
        sum(hi[a][b] * radius_gradient[b] for b in range(2)) for a in range(2)
    )
    for a in range(2):
        for i in (2, 3):
            gamma[a][i][i] = -radius * raised_radius_gradient[a]
            gamma[i][a][i] = radius_gradient[a] / radius
            gamma[i][i][a] = radius_gradient[a] / radius
    return g, gi, gamma


def _hessian(j: Jet2, gamma):
    d = (j.dt, j.dr, Q(0), Q(0))
    dd = (
        (j.dtt, j.dtr, Q(0), Q(0)),
        (j.dtr, j.drr, Q(0), Q(0)),
        (Q(0), Q(0), Q(0), Q(0)),
        (Q(0), Q(0), Q(0), Q(0)),
    )
    return [
        [dd[a][b] - sum(gamma[c][a][b] * d[c] for c in range(N)) for b in range(N)]
        for a in range(N)
    ]


def _p(ri, ric, sc, g):
    p = _zero4()
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    p[a][b][c][d] = (
                        ri[a][b][c][d]
                        - g[a][c] * ric[d][b]
                        + g[a][d] * ric[c][b]
                        + g[b][c] * ric[d][a]
                        - g[b][d] * ric[c][a]
                        + sc * (g[a][c] * g[d][b] - g[a][d] * g[c][b]) / 2
                    )
    return p


def residuals(state: SphericalState, *, use_warped: bool = False):
    if use_warped:
        ri = warped_2plus2_curvature(state)
        g, gi, gamma = _warped_2plus2_metric_connection(state)
    else:
        ri, gi, gamma = direct_4d_curvature(state)
        g, _dg, _dd = _direct_metric(state)
    ric = _ricci(ri, gi)
    sc = _scalar(ric, gi)
    gb = _gb(ri, ric, sc, gi)
    hp = _hessian(state.phi, gamma)
    hc = _hessian(state.chi, gamma)
    F = state.planck_mass**2 + state.beta * state.phi.value**2
    fp = 2 * state.beta * state.phi.value
    fpp = 2 * state.beta
    f1 = (
        state.alpha
        if state.branch == "SGB-L"
        else state.eta * state.phi.value / 4
        if state.branch == "FGC-QR"
        else Q(0)
    )
    f2 = state.eta / 4 if state.branch == "FGC-QR" else Q(0)
    hf = [
        [
            f1 * hp[a][b]
            + f2
            * (state.phi.dt, state.phi.dr, Q(0), Q(0))[a]
            * (state.phi.dt, state.phi.dr, Q(0), Q(0))[b]
            for b in range(N)
        ]
        for a in range(N)
    ]
    hF = [
        [
            fp * hp[a][b]
            + fpp
            * (state.phi.dt, state.phi.dr, Q(0), Q(0))[a]
            * (state.phi.dt, state.phi.dr, Q(0), Q(0))[b]
            for b in range(N)
        ]
        for a in range(N)
    ]
    box = lambda h: sum(gi[a][b] * h[a][b] for a in range(N) for b in range(N))
    bF = box(hF)
    p = _p(ri, ric, sc, g)
    ein = [[ric[a][b] - g[a][b] * sc / 2 for b in range(N)] for a in range(N)]

    def stress(h, j, V):
        d = (j.dt, j.dr, Q(0), Q(0))
        ds = sum(gi[a][b] * d[a] * d[b] for a in range(N) for b in range(N))
        return [
            [d[a] * d[b] - g[a][b] * (ds / 2 + V) for b in range(N)] for a in range(N)
        ]

    V = state.mu**2 * state.phi.value**2 / 2 + state.g4 * state.phi.value**4 / 4
    tp = stress(hp, state.phi, V)
    tc = stress(hc, state.chi, Q(0))
    E = _zero2()
    einstein_term = _zero2()
    nonminimal_term = _zero2()
    gb_residual_term = _zero2()
    for a in range(N):
        for b in range(N):
            einstein_term[a][b] = F * ein[a][b]
            nonminimal_term[a][b] = g[a][b] * bF - hF[a][b]
            gb_residual_term[a][b] = 8 * sum(
                p[a][c][b][d] * gi[c][e] * gi[d][z] * hf[e][z]
                for c in range(N)
                for d in range(N)
                for e in range(N)
                for z in range(N)
            )
            E[a][b] = (
                einstein_term[a][b]
                + nonminimal_term[a][b]
                - tp[a][b]
                - tc[a][b]
                + gb_residual_term[a][b]
            )
    boxp = box(hp)
    boxc = box(hc)
    vp = state.mu**2 * state.phi.value + state.g4 * state.phi.value**3
    Ephi = boxp - vp + state.beta * state.phi.value * sc + f1 * gb
    Echi = boxc
    freeze = lambda matrix: tuple(tuple(x for x in row) for row in matrix)
    return {
        "riemann": ri,
        "ricci": ric,
        "R": sc,
        "GB": gb,
        "hessian_phi": hp,
        "hessian_chi": hc,
        "einstein_term": freeze(einstein_term),
        "nonminimal_term": freeze(nonminimal_term),
        "phi_stress": freeze(tp),
        "chi_stress": freeze(tc),
        "gb_residual_term": freeze(gb_residual_term),
        "metric": freeze(E),
        "phi": Ephi,
        "chi": Echi,
        "F": F,
    }


def _independent_residual_vector(result: Mapping[str, Any]) -> tuple[Fraction, ...]:
    metric = result["metric"]
    return (
        metric[0][0],
        metric[0][1],
        metric[1][1],
        metric[2][2],
        result["phi"],
        result["chi"],
    )


def principal_block(state: SphericalState):
    """Exact frozen-jet first tangent by symmetric rational directional AD.

    Each column is d(residual)/d(field second-jet).  The local residual is at
    most quadratic in second jets, so the central rational difference is exact.
    This is not a first-order characteristic matrix.
    """
    keys = [(q, d) for q in BASE_FIELD_ORDER for d in SECOND_DERIVATIVE_ORDER]
    out = {}
    for q, d in keys:
        j = getattr(state, q)
        plus = replace(state, **{q: replace(j, **{d: getattr(j, d) + 1})})
        minus = replace(state, **{q: replace(j, **{d: getattr(j, d) - 1})})
        a, b = residuals(plus), residuals(minus)
        out[f"{q}.{d}"] = {
            "metric": tuple(
                tuple((a["metric"][i][j] - b["metric"][i][j]) / 2 for j in range(N))
                for i in range(N)
            ),
            "phi": (a["phi"] - b["phi"]) / 2,
            "chi": (a["chi"] - b["chi"]) / 2,
        }
    return out


def principal_matrix(state: SphericalState) -> dict[str, Any]:
    """Return the independent six-by-eighteen ``2+2`` coefficient matrix."""

    block = principal_block(state)
    columns = tuple(block)
    rows = []
    for equation_index in range(len(INDEPENDENT_EQUATION_ORDER)):
        row = []
        for column in columns:
            entry = block[column]
            vector = (
                entry["metric"][0][0],
                entry["metric"][0][1],
                entry["metric"][1][1],
                entry["metric"][2][2],
                entry["phi"],
                entry["chi"],
            )
            row.append(vector[equation_index])
        rows.append(tuple(row))
    return {
        "field_order": BASE_FIELD_ORDER,
        "second_derivative_order": SECOND_DERIVATIVE_ORDER,
        "column_order": columns,
        "equation_order": INDEPENDENT_EQUATION_ORDER,
        "matrix": tuple(rows),
    }


def adm_pg_principal_matrix(fixture: Mapping[str, Any]) -> dict[str, Any]:
    """Extract the exact coefficient matrix in the declared ADM/PG variables.

    Only second-jet slots are perturbed. The nonlinear ADM-to-``2+2`` field
    map is reevaluated on each side, so its exact highest-derivative Jacobian
    is included rather than silently relabeling the base-metric columns.
    """

    jets = _adm_pg_jets_from_fixture(fixture)
    columns = tuple(
        f"{field}.{derivative}"
        for field in ADM_FIELD_ORDER
        for derivative in SECOND_DERIVATIVE_ORDER
    )
    column_vectors: list[tuple[Fraction, ...]] = []
    for field in ADM_FIELD_ORDER:
        for derivative in SECOND_DERIVATIVE_ORDER:
            plus_jets = dict(jets)
            minus_jets = dict(jets)
            original = jets[field]
            plus_jets[field] = replace(
                original, **{derivative: getattr(original, derivative) + 1}
            )
            minus_jets[field] = replace(
                original, **{derivative: getattr(original, derivative) - 1}
            )
            plus = residuals(_state_from_adm_pg_jets(fixture, plus_jets))
            minus = residuals(_state_from_adm_pg_jets(fixture, minus_jets))
            plus_vector = _independent_residual_vector(plus)
            minus_vector = _independent_residual_vector(minus)
            column_vectors.append(
                tuple((a - b) / 2 for a, b in zip(plus_vector, minus_vector))
            )

    rows = tuple(
        tuple(column[equation] for column in column_vectors)
        for equation in range(len(INDEPENDENT_EQUATION_ORDER))
    )
    return {
        "field_order": ADM_FIELD_ORDER,
        "second_derivative_order": SECOND_DERIVATIVE_ORDER,
        "column_order": columns,
        "equation_order": INDEPENDENT_EQUATION_ORDER,
        "matrix": rows,
    }


def adm_pg_kinematics(lapse: Fraction, radial_metric: Fraction, shift: Fraction):
    lapse, radial_metric, shift = _q(lapse), _q(radial_metric), _q(shift)
    if lapse <= 0 or radial_metric <= 0:
        raise ValueError("lapse and Lambda must be positive")
    return {
        "v_out": -shift + lapse / radial_metric,
        "v_in": -shift - lapse / radial_metric,
        "base_metric_determinant": -((lapse * radial_metric) ** 2),
        "regular_lorentzian_chart": True,
        "dynamic_lambda": radial_metric != 1,
    }


def _maximum_tensor_difference(left: Any, right: Any) -> Fraction:
    if isinstance(left, Fraction) and isinstance(right, Fraction):
        return abs(left - right)
    if isinstance(left, (list, tuple)) and isinstance(right, (list, tuple)):
        if len(left) != len(right):
            raise ValueError("tensor shapes differ")
        return max(
            (_maximum_tensor_difference(a, b) for a, b in zip(left, right)),
            default=Q(0),
        )
    raise ValueError("tensor types differ")


def _maximum_abs_tensor(value: Any) -> Fraction:
    if isinstance(value, Fraction):
        return abs(value)
    if isinstance(value, (list, tuple)):
        return max((_maximum_abs_tensor(item) for item in value), default=Q(0))
    raise ValueError("tensor contains a non-rational value")


def _matrix_is_zero(matrix: Any) -> bool:
    if isinstance(matrix, Fraction):
        return matrix == 0
    return all(_matrix_is_zero(value) for value in matrix)


def _serialize_principal_matrix(data: Mapping[str, Any]) -> dict[str, Any]:
    matrix = [[_fraction_text(value) for value in row] for row in data["matrix"]]
    canonical = json.dumps(matrix, separators=(",", ":"), ensure_ascii=True)
    return {
        "field_order": list(data["field_order"]),
        "second_derivative_order": list(data["second_derivative_order"]),
        "column_order": list(data["column_order"]),
        "equation_order": list(data["equation_order"]),
        "matrix": matrix,
        "matrix_sha256": sha256(canonical.encode("ascii")).hexdigest(),
        "row_count": len(matrix),
        "column_count": len(matrix[0]) if matrix else 0,
    }


def _chi_principal_control(
    state: SphericalState, matrix: Mapping[str, Any]
) -> dict[str, Any]:
    _riemann, inverse, _gamma = direct_4d_curvature(state)
    expected = (inverse[0][0], 2 * inverse[0][1], inverse[1][1])
    columns = matrix["column_order"]
    rows = matrix["matrix"]
    indices = tuple(columns.index(f"chi.{slot}") for slot in SECOND_DERIVATIVE_ORDER)
    actual = tuple(rows[5][index] for index in indices)
    other_zero = all(rows[row][index] == 0 for row in range(5) for index in indices)
    return {
        "expected_metric_null_coefficients": [_fraction_text(v) for v in expected],
        "actual_chi_equation_coefficients": [_fraction_text(v) for v in actual],
        "chi_second_derivatives_absent_from_other_equations": other_zero,
        "chi_principal_factor_is_metric_null": actual == expected and other_zero,
    }


def _required_red1_nonclaims() -> dict[str, bool]:
    return {
        k: False
        for k in (
            "first_order_reduction_derived",
            "evolution_constraint_split_derived",
            "constraint_propagation_proven",
            "kinetic_matrix_invertible_on_retained_domain",
            "physical_characteristic_polynomial_derived",
            "all_characteristics_real_on_retained_domain",
            "complete_characteristic_basis_or_symmetrizer_proven",
            "strong_hyperbolicity_proven",
            "evolution_authorized",
            "metric_null_defocusing_derived",
        )
    }


def flat_flrw_curvature_control() -> dict[str, Any]:
    """Return an exact analytic curvature regression for spatially flat FLRW.

    The point has ``a=2``, ``H=3/2``, ``H_dot=-2/3``, and radial coordinate
    ``x=5``. Hence ``R=23`` and ``G=171/2`` independently of the spherical
    tensor construction.
    """

    state = SphericalState(
        h_tt=Jet2(-1),
        h_tr=Jet2(0),
        h_rr=Jet2(4, 12, 0, Q(92, 3), 0, 0),
        areal_radius=Jet2(10, 15, 2, Q(95, 6), 3, 0),
        phi=Jet2(0),
        chi=Jet2(0),
        planck_mass=Q(2),
        mu=Q(3),
        g4=Q(1, 2),
        branch="GR-0",
    )
    result = residuals(state)
    return {
        "inputs": {
            "scale_factor": "2",
            "hubble": "3/2",
            "hubble_derivative": "-2/3",
            "radial_coordinate": "5",
        },
        "expected_R": "23",
        "actual_R": _fraction_text(result["R"]),
        "expected_GB": "171/2",
        "actual_GB": _fraction_text(result["GB"]),
        "exact_match": result["R"] == 23 and result["GB"] == Q(171, 2),
    }


def red1_certificate(state: SphericalState):
    direct = residuals(state)
    warped = residuals(state, use_warped=True)
    matrix = principal_matrix(state)
    differences = {
        "riemann": _maximum_tensor_difference(direct["riemann"], warped["riemann"]),
        "ricci": _maximum_tensor_difference(direct["ricci"], warped["ricci"]),
        "R": abs(direct["R"] - warped["R"]),
        "GB": abs(direct["GB"] - warped["GB"]),
        "connection": _maximum_tensor_difference(
            direct_4d_curvature(state)[2],
            _warped_2plus2_metric_connection(state)[2],
        ),
        "hessian_phi": _maximum_tensor_difference(
            direct["hessian_phi"], warped["hessian_phi"]
        ),
        "hessian_chi": _maximum_tensor_difference(
            direct["hessian_chi"], warped["hessian_chi"]
        ),
        "metric_equation": _maximum_tensor_difference(
            direct["metric"], warped["metric"]
        ),
        "scalar_phi": abs(direct["phi"] - warped["phi"]),
        "scalar_chi": abs(direct["chi"] - warped["chi"]),
    }
    return {
        "classification": "exact_local_jet_reduction_not_first_order_or_hyperbolicity",
        "direct_warped_max_exact_residuals": {
            key: _fraction_text(value) for key, value in differences.items()
        },
        "invariants": {
            "R": _fraction_text(direct["R"]),
            "GB": _fraction_text(direct["GB"]),
            "F": _fraction_text(direct["F"]),
        },
        "source_sector_controls": {
            "nonminimal_term_zero": _matrix_is_zero(direct["nonminimal_term"]),
            "gb_bulk_residual_term_zero": _matrix_is_zero(direct["gb_residual_term"]),
        },
        "principal_part": _serialize_principal_matrix(matrix),
        "chi_principal_control": _chi_principal_control(state, matrix),
        "nonclaims": _required_red1_nonclaims(),
    }


def spherical_reduction_certificate(configuration: Mapping[str, Any]) -> dict[str, Any]:
    """Build a fail-closed aggregate RED1 certificate from ADM/PG fixtures.

    This adapter intentionally accepts only declared local fixtures.  It does
    not accept a solver grid, derive constraints, or promote its directional
    second-jet block into a characteristic matrix.
    """
    if not isinstance(configuration, Mapping):
        raise ValueError("configuration must be a mapping")
    allowed = {"fixtures", "schema_version", "artifact_id"}
    if set(configuration) - allowed:
        raise ValueError("RED1 configuration has unknown keys")
    raw = configuration.get("fixtures")
    if not isinstance(raw, (list, tuple)) or not raw:
        raise ValueError("RED1 requires a nonempty fixture list")
    records = []
    for index, fixture in enumerate(raw):
        if not isinstance(fixture, Mapping):
            raise ValueError("RED1 fixture must be a mapping")
        state = state_from_generalized_adm_pg_fixture(fixture)
        item = red1_certificate(state)
        adm_matrix = adm_pg_principal_matrix(fixture)
        item["principal_part"] = _serialize_principal_matrix(adm_matrix)
        item["chi_principal_control"] = _chi_principal_control(state, adm_matrix)
        item["fixture_id"] = str(fixture.get("fixture_id", index))
        item["branch"] = state.branch
        item["purpose"] = str(
            fixture.get("purpose", "unspecified_local_algebra_fixture")
        )
        raw_state = fixture.get("state", fixture)

        def v(name: str):
            raw = raw_state[name]
            return (
                _fixture_jet(raw, name).value if isinstance(raw, Mapping) else _q(raw)
            )

        item["adm_pg_kinematics"] = {
            key: (_fraction_text(value) if isinstance(value, Fraction) else value)
            for key, value in adm_pg_kinematics(
                v("alpha"), v("lambda"), v("shift")
            ).items()
        }
        direct = residuals(state)
        if (
            item["purpose"]
            == "zero_regulator_ricci_flat_schwarzschild_exterior_control"
        ):
            item["schwarzschild_vacuum_control"] = {
                "metric_residual_zero": _matrix_is_zero(direct["metric"]),
                "scalar_residuals_zero": direct["phi"] == 0 and direct["chi"] == 0,
                "maximum_abs_metric_residual": _fraction_text(
                    _maximum_abs_tensor(direct["metric"])
                ),
                "absolute_phi_residual": _fraction_text(abs(direct["phi"])),
                "absolute_chi_residual": _fraction_text(abs(direct["chi"])),
                "ricci_scalar_zero": direct["R"] == 0,
                "gauss_bonnet_invariant": _fraction_text(direct["GB"]),
            }
        if (
            item["purpose"]
            == "nonzero_regulator_gradient_hessian_principal_part_control"
        ):
            rows = adm_matrix["matrix"]
            columns = adm_matrix["column_order"]
            metric_columns = [
                i
                for i, name in enumerate(columns)
                if name.split(".")[0] in {"alpha", "shift", "lambda", "areal_radius"}
            ]
            phi_columns = [
                i for i, name in enumerate(columns) if name.startswith("phi.")
            ]
            item["activated_mixing_control"] = {
                "metric_equations_depend_on_phi_second_jets": any(
                    rows[row][column] for row in range(4) for column in phi_columns
                ),
                "phi_equation_depends_on_metric_second_jets": any(
                    rows[4][column] for column in metric_columns
                ),
                "nonzero_phi_value_gradient_and_hessian": state.phi.value != 0
                and (state.phi.dt != 0 or state.phi.dr != 0)
                and (state.phi.dtt != 0 or state.phi.dtr != 0 or state.phi.drr != 0),
            }
        records.append(item)
    direct_warped_exact = all(
        all(
            value == "0" for value in item["direct_warped_max_exact_residuals"].values()
        )
        for item in records
    )
    schwarzschild_controls = [
        item["schwarzschild_vacuum_control"]
        for item in records
        if "schwarzschild_vacuum_control" in item
    ]
    activated_controls = [
        item["activated_mixing_control"]
        for item in records
        if "activated_mixing_control" in item
    ]
    nonclaims = _required_red1_nonclaims()
    flrw_control = flat_flrw_curvature_control()
    return {
        "schema_version": configuration.get("schema_version", 1),
        "artifact_id": configuration.get("artifact_id", "FGC-1-HYP1-RED1"),
        "classification": "exact_local_jet_reduction_not_first_order_or_hyperbolicity",
        "fixtures": records,
        "analytic_curvature_controls": {
            "spatially_flat_flrw": flrw_control,
        },
        "aggregate": {
            "fixture_count": len(records),
            "all_direct_warped_exact": direct_warped_exact,
            "all_principal_matrices_six_by_eighteen": all(
                item["principal_part"]["row_count"] == 6
                and item["principal_part"]["column_count"] == 18
                for item in records
            ),
            "all_chi_principal_factors_metric_null": all(
                item["chi_principal_control"]["chi_principal_factor_is_metric_null"]
                for item in records
            ),
            "schwarzschild_control_present_and_exact": len(schwarzschild_controls) == 1
            and schwarzschild_controls[0]["metric_residual_zero"]
            and schwarzschild_controls[0]["scalar_residuals_zero"]
            and schwarzschild_controls[0]["maximum_abs_metric_residual"] == "0"
            and schwarzschild_controls[0]["absolute_phi_residual"] == "0"
            and schwarzschild_controls[0]["absolute_chi_residual"] == "0"
            and schwarzschild_controls[0]["ricci_scalar_zero"]
            and schwarzschild_controls[0]["gauss_bonnet_invariant"] == "3/16384",
            "activated_mixing_control_present_and_exercised": len(activated_controls)
            == 1
            and all(all(control.values()) for control in activated_controls),
            "flat_flrw_curvature_control_exact": flrw_control["exact_match"],
            "all_nonclaims_false": all(not value for value in nonclaims.values()),
        },
        "nonclaims": nonclaims,
    }
