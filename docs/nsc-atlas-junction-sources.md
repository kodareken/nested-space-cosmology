# Can a fixed vacuum wall support a comoving radiation child?

10 October 2026. Owner:
[source adapter](../src/recursive_horizons/nsc_atlas_junction_sources.py),
[independent checks](../tests/test_nsc_atlas_junction_sources.py).
This is the first finite source-compatibility screen following the
[exterior-first atlas](nsc-cosmological-atlas.md) and the existing
[compute admission rule](../../docs/prior-art-and-open-claim.md).

**Result:** a constant-tension vacuum wall cannot join an exact closed-radiation
FLRW region along a comoving boundary to a fixed-mass Schwarzschild region over
any nonzero interval in which the radius changes. This statement applies to
either sign of the constant surface energy and to either small-cap or
large-region orientation. It follows from the already established junction
and conservation laws. It is a conditional compatibility result for the map,
with no new-to-world claim.

The failure belongs to this particular combination of bulk geometries, source
and comoving boundary. Moving walls, different interiors, other exterior
states and different independently justified sources remain distinct
candidates. In particular, the published de Sitter/Schwarzschild vacuum-wall
models are not excluded. No actual child production, stability, causal
connection to another child, or embedding into two-hole initial data is
established by this calculation.

## Admission before evaluation

| Required field | This finite test |
|---|---|
| Borrowed result | Israel surface stresses and shell conservation; Lobo–Crawford `gr-qc/0507063v2`, equations (3), (4), (9), (10), (24); the vacuum-wall source in Ansoldi–Guendelman–Shilon `0711.2198v1`, section 3, and Guendelman–Sakai `0803.0268v2`, equations (4)–(6). The local implementation precedent is the [Israel note](child-neck/israel-matching-note.md) and its [script](../scripts/child_neck/israel_matching.py). |
| Missing connection | Does an independently specified fixed vacuum wall support the exact radiation world tube, including both extrinsic-curvature components and energy conservation? |
| Decision and unlock | Exclude this pairing if inconsistent. A subsequent candidate must change an identified assumption and supply its physical source; no adjustment of the stress function to the desired geometry is accepted as source closure. |
| Finite domain | Four-dimensional classical GR, (G=c=1); (A=1), (M/A=0.4), (chi_0=pi/3), (1.25leqetaleq1.45). All 33 witness points satisfy (R>2M). Proper time is the boundary observer's time. |
| Budget and stop | One short polynomial argument, 33 witness points, six differential controls and nine dust controls; less than one CPU minute, no parameter scan or PDE evolution. Stop this pairing at incompatibility. |
| Verification and attribution | Independent metric acceleration, finite-difference conservation, normal Einstein jump, shell-free Oppenheimer–Snyder control and symbolic elimination. Literature search is bounded; no priority assertion. |

`junction_record()` returns this admission, frozen input parameters, the
analytic obstruction, all witness samples and the independent errors as a
JSON-serializable record. It does not write files or alter existing evidence.
The atlas producer can bind it into a successor record.

## Exact geometry and conventions

The child region has

\[
ds_-^2=-dt^2+a(t)^2(d\chi^2+\sin^2\chi\,d\Omega^2),\qquad
a=A\sin\eta,\quad t=A(1-\cos\eta),\quad 0<\eta<\pi.
\]

These functions satisfy the exact closed radiation equations,

\[
\dot a^2+1=\frac{A^2}{a^2},\qquad
\rho=\frac{3A^2}{8\pi a^4},\quad p=\rho/3.
\]

Dots denote the shell proper time (	au=t), because the boundary is
comoving at fixed (chi_0). Put (s=sinchi_0), (c_0=coschi_0).
Then

\[
R=as,\quad \dot R=s\cot\eta,\quad
\ddot R=-\frac{A^2s}{a^3}=-\frac{A^2s^4}{R^3}.
\]

The retained child region is (0\leq\chi\leq\chi_0). For
(\chi_0<\pi/2) it is a cap; for (\chi_0>\pi/2) it occupies more than
half the spherical slice. Its outward normal points toward increasing
(\chi), so the angular normal derivative is (eta_-=c_0), including
its sign. This is a spherical spatial region, not a claim that a complete
closed universe has already formed behind a parent horizon.

The exterior is Schwarzschild with fixed (M\geq0),
(f_+=1-2M/R). The numerical witness retains the branch with outward
(\beta_+>0) and (f_+>0); its time tangent is future directed. Thus

\[
\beta_+=\sqrt{f_++\dot R^2},\qquad
\beta_+^2=c_0^2+\frac{A^2s^4}{R^2}-\frac{2M}{R}.
\]

The adapter rejects the horizon or static-chart exit. Its finite-domain
restriction is narrower than the algebraic argument below. The witness is
an expanding radiation cap outside the Schwarzschild horizon, not a
simulation of black-hole formation. A Brill–Lindquist parent supplies
different, nonspherical spatial data; matching these two candidates has
not been solved.

There is no density-cancellation bounce in this test. At the exact maximum
radius, (dot a=0) while (1+dot a^2=1), and the boundary's
Misner–Sharp mass is (asin^3chi_0/2>0). The domain excludes the
singular endpoints. The earlier torsion approximation is not imported.

## Full shell requirement and the independently supplied source

With normal directed from minus to plus and jump ([X]=X_+-X_-), use

\[
K^\theta{}_{\theta\pm}=\frac{\beta_\pm}{R},\qquad
K^\tau{}_{\tau-}=0,\qquad
K^\tau{}_{\tau+}=\frac{\ddot R+M/R^2}{\beta_+}.
\]

The borrowed Israel laws give the **required** stress:

\[
\sigma_{\rm req}=\frac{c_0-\beta_+}{4\pi R},\qquad
P_{\rm req}=\frac{1}{8\pi}
\left(\frac{\ddot R+M/R^2}{\beta_+}
+\frac{\beta_+-c_0}{R}\right).
\]

Computing this stress does not supply its physical origin. Both components
must agree with a separately specified source. The candidate used here is
the conventional thin vacuum wall,

\[
\sigma_{\rm wall}=\sigma_0>0,\qquad P_{\rm wall}=-\sigma_0.
\]

Its illustrative value (sigma_0=1/(64\pi A)) is prescribed independently,
without a calibration point, fitted stress function or retuning. The
analytic conclusion does not depend on that value.

The comoving perfect fluid has no (T_{\mu\nu}n^\mu u^\nu) flux, and the
exterior is vacuum. Consequently the conservation law is

\[
\dot\sigma+2\frac{\dot R}{R}(\sigma+P)=0.
\]

The wall equation of state therefore requires constant (sigma). The
geometrically required stress satisfies this conservation law, but does
not satisfy the vacuum-wall equation of state. Radiation pressure is not
ignored: it enters the trajectory's acceleration and the temporal
extrinsic curvature. An independent normal Einstein-constraint jump
also checks ([G_{nn}]/(8\pi)=-p_{m child}), using the same stated normal
conventions.

The existing [balanced layer note](child-neck/balanced-layers-note.md)
already identifies an unprovided energy-exchange mechanism for its
torsion/vacuum mixture. Its static planar balance is not a replacement for
this dynamical spherical source closure. This test neither validates nor
excludes that separate composite candidate.

## Finite-interval obstruction

Let (k=4\pi\sigma_0). Angular matching requires
(eta_+=c_0-kR). Squaring gives a necessary condition

\[
F(R)=k^2R^4-2c_0kR^3+2MR-A^2s^4=0.
\]

If matching held on a nonzero time interval with varying radius, this
polynomial would vanish on a radius interval, so every coefficient would
vanish. Its constant coefficient is (-A^2s^4\ne0) for a nontrivial child
(A>0, 0<\chi_0<\pi). Therefore no constant (k), including zero or
negative (k), works on the interval.

This argument covers both signs of (c_0) and either exterior normal
branch: squaring is only a necessary condition and cannot remove the
contradiction. A root at an isolated radius would still require the
unsquared sign condition and the temporal junction equation. It is not a
finite evolution. A turning point does not invalidate the argument for a
surrounding interval whose radius varies.

This obstruction is a direct application of imported laws to the proposed
attachment. The map's gain is identifying why this specific edge cannot
be filled by a fixed vacuum wall. A larger simulation cannot repair this
incompatibility while retaining these assumptions.

## Witness and independent checks

The witness uses (A=1, M/A=0.4, \chi_0=\pi/3) and 33 evenly spaced
conformal-time points. The shell's proper-time interval is
(0.684677638\leq\tau/A\leq0.879497231). The Schwarzschild horizon would
be at (R/A=0.8).

| Quantity | (eta=1.25) | (eta=1.45) |
|---|---:|---:|
| (R/A) | 0.821844788 | 0.859714669 |
| (A\sigma_{m req}) | 0.016389748 | 0.020017358 |
| (AP_{m req}) | −0.058857430 | −0.058253512 |
| (A\sigma_{m wall}), fixed | 0.004973592 | 0.004973592 |
| (A(P_{m req}+\sigma_{m req})) | −0.042467681 | −0.038236154 |

The results remain dimensionlessly identical when (A,M,R,	au) are
rescaled together and (sigma,P) scale inversely. The tests also change
the Schwarzschild mass while leaving the physical wall fixed; no source
parameter is fitted to those cases. That is a consistency check, not a
perturbative stability proof of an existing solution.

Independent controls use the static metric's Christoffel symbols and a
five-point proper-time derivative of the four-velocity, a separate
five-point derivative of (sigma) for Codazzi, and the difference of the
normal Einstein constraints. Both cap orientations are covered. In units
(A=1), the recorded maximum errors are approximately (1.21\times10^{-10})
for (K^\tau{}_{\tau+}), (9.16\times10^{-13}) for conservation and
(6.94\times10^{-18}) for the normal constraint. The derivative step is
(10^{-5}A); these are numerical integration checks, not the analytic proof.

A positive matching control uses the published Oppenheimer–Snyder dust
cap: (a=A(1-\cos\eta)/2), comoving boundary and
(M=A\sin^3\chi_0/2). Its shell must vanish. Nine points give maximum
(|\sigma|<1.9\times10^{-17}), (|P|<8.9\times10^{-18}).
The adapter thus distinguishes the radiation incompatibility from the
known shell-free dust match. The focused test file has 17 passing cases.

## Literature absorbed and remaining boundaries

| Source and version | Reading and use | Remaining boundary |
|---|---|---|
| [Lobo–Crawford, gr-qc/0507063v2](https://arxiv.org/html/gr-qc/0507063v2) | Sections 2, 3.1 and 3.2: general extrinsic-curvature definition, spherical surface stress, no-flux conservation. | Their later static-background stability analysis is not applied to this unsatisfied radiation candidate. The FLRW interior is evaluated from its own metric, not treated as static. |
| [Ansoldi–Guendelman–Shilon, 0711.2198v1](https://arxiv.org/html/0711.2198v1) | Sections 2–3: shell orientation, fixed vacuum-wall source and the de Sitter/Schwarzschild child-universe model. | The review identifies Blau–Guendelman–Guth (1987) as prior art; the original paper was not independently reread for this screen. Published vacuum bubbles have different interiors and wall trajectories. |
| [Guendelman–Sakai, 0803.0268v2](https://arxiv.org/html/0803.0268v2) | Section II, equations (4)–(6), and sections IV–V: a supplied domain-wall-plus-dust EOS, orientation, proposed physical origins and mass exchange. The [2009 erratum](https://journals.aps.org/prd/pdf/10.1103/PhysRevD.80.049901) was read in full; it corrects potential equations (9) and (15), neither used here. | Mixed sources and their dynamics are separate candidates. The paper's assumed negative-energy contributions do not establish a microscopic source for this project. No breathing-bubble solution is claimed here. |

The focused search on 10 October 2026 covered these sources, their cited
vacuum-wall construction, and queries combining *comoving radiation FLRW*,
*Schwarzschild*, *domain wall*, *constant tension* and *junction*. It did
not establish absence of the same elementary compatibility argument
elsewhere. The contribution recorded here is a reproducible application
that removes an unsupported atlas attachment, with its replacement
questions visible.

The next useful choice is between borrowing a published moving
de Sitter/Schwarzschild wall solution as a positive atlas example, or
admitting a radiation boundary that moves relative to the fluid and hence
requires a physical flux/source law. Neither option follows merely by
drawing a connection on the map. The two-hole parent embedding, formation
history and observational consequences remain open.
