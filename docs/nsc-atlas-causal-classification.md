# Causal classification of the isolated sourced body

**Status, 10 October 2026:** an analytic identification within the published
classical thin-wall model. The fixed body A has a monotonic expanding
reference history, emerging through a white-hole sector into the opposite
Schwarzschild exterior and an expanding de Sitter interior. This note
classifies that isolated spherical history. It does not evolve the glued
many-body parent or replace any sealed v1-v4 result.

The [v4 connection slice](nsc-atlas-connection-slice.md) fixes
$H\ell=1$, $k\ell=2$, $M/\ell=2/5$, $R_w/\ell=3/2$ and the expanding
wall tangent. Here $k=4\pi\sigma$, with positive constant tension. The
[gluing construction](nsc-atlas-connection-gluing.md) supplies a separate
shared **initial-data** family; it does not import this whole spherical
reference future into that family.

## Admission and borrowed classification

- **Borrow:** the signed junction law and branch catalogue in
  [Aguirre and Johnson, gr-qc/0508093v2](https://arxiv.org/html/gr-qc/0508093v2),
  sections II.1-II.3, equations (1)-(20), Figures 7 and 10; the explicitly
  asymptotically flat bubble discussion in
  [Ansoldi, Guendelman and Shilon, 0711.2198v1](https://arxiv.org/html/0711.2198v1),
  sections 2-3; and Schwarzschild region/null-coordinate conventions in
  [Kuchar, gr-qc/9403003v1](https://arxiv.org/html/gr-qc/9403003v1), section 2.1,
  equations (5)-(13).
- **Missing connection:** identify which global reference branch contains
  the already checked cap, signed wall and bridge, then distinguish
  observer access from spatial attachment.
- **Decision/unlock:** determine the retained patches and causal boundaries
  before deciding whether a proposed signal, perturbation or second child
  would change the map.
- **Domain/budget:** one fixed theory and state; analytic sign and
  asymptotic checks on the published classical radial model for $R>0$.
  No simulation, parameter scan, gluing PDE, new background or fit is admitted.
  Asymptotic boundaries classify this one spacetime, not an infinite hierarchy.
- **Stop/redirect:** a sign inconsistency, missed turning point, incompatible
  retained patch or dependence on an unavailable matter evolution theorem
  stops the corresponding access claim. Numerical observer witnesses, if
  used by the successor atlas, keep their own finite interval and admission.
- **Independent check:** relate both signed normal factors to horizon
  crossings and regular null coordinates, rather than relying on a squared
  wall potential or a visual resemblance to a published diagram.

This is an application of an existing classification. No new thin-wall
dynamics, formation mechanism or priority claim is made.

## Fixed-source branch and a check without a radial simulation

Put $x=R/\ell$ and $\hat\tau=\tau_w/\ell$. The inherited signed factors,
with the normal pointing from de Sitter to Schwarzschild, are

$$
\beta_-={1\over5x^2}+{3x\over4},\qquad
\beta_+={1\over5x^2}-{5x\over4},\qquad
\beta_- -\beta_+=2x.
$$

The expanding solution obeys

$$
{dx\over d\hat\tau}=\sqrt{q(x)},\qquad
q(x)={25x^2\over16}-1+{3\over10x}+{1\over25x^4}.
$$

There is no turning radius. For $x\geq4/5$, the first term is at least
one and the remaining positive terms make $q>0$. For $0<x\leq4/5$, AM-GM
gives

$$
{25x^2\over32}+{25x^2\over32}+{1\over25x^4}
\geq3\left({25\over1024}\right)^{1/3}>{3\over4},
\qquad {3\over10x}\geq{3\over8},
$$

so $q>1/8$. This check distinguishes the monotonic branch from an expanding
segment of a bouncing solution.

The interior factor is positive everywhere. The exterior factor changes
sign only at $x_\beta=(4/25)^{1/3}$, which is smaller than $4/5$. This
normal-direction change is not a radial turning point or horizon. The two
relevant horizon radii are instead $x_S=2M/\ell=4/5$ and
$x_{dS}=1/(H\ell)=1$. The ordered wall events are

$$
0<x_\beta<{4\over5}<1<{3\over2}.
$$

At small radius, $dx/d\hat\tau\sim1/(5x^2)$; the past wall endpoint is
at finite proper time, with
$\hat\tau-\hat\tau_{\rm birth}\sim5x^3/3$. At large radius,
$dx/d\hat\tau\sim5x/4$, so the wall reaches $R=\infty$ only at infinite
proper time, $\hat\tau\sim(4/5)\log x$ plus a constant. These are formal
boundaries of the ideal classical model, not claims that the thin-wall
approximation remains physically reliable at arbitrarily high curvature.

## Relation to the published catalogue and the flat exterior limit

The Aguirre-Johnson variables give
$\Lambda_+=0$, $\Lambda_-=3H^2$, $A=\Lambda_+/k^2=0$,
$B=\Lambda_-/k^2=3/4$, $L^2=H^2+k^2=5/\ell^2$ and $Y=3/5$.
Their potential criterion for monotonic motion is equivalent to $q>0$
above. With no interior normal-sign change, the appropriate pattern is
the $\Lambda_+\to0^+$ limit of **solution 13 in Figure 10**. The original
figure and Figure 7 were inspected, including the shaded retained regions.

For positive exterior cosmological constant that pattern is labelled
IV-III-II'' / III-II. Its finite-radius equations have a regular
$\Lambda_+=0$ limit, but its outer cosmological horizon does not remain at a
finite radius: $r_c\to\infty$, while $r_{BH}\to2M$. Accordingly our
Schwarzschild sequence is **white-hole IV to opposite exterior III**,
with de Sitter sequence **central static III to expanding II**, using their
labels. There is no finite-radius exterior II'' crossing in Schwarzschild.
The asymptotically flat limit must change the boundary in the drawing; it
cannot retain a finite Schwarzschild-de Sitter cosmological horizon.

The factors themselves are $\beta_\pm=n^\mu\partial_\mu R$.
Aguirre-Johnson equation (13) uses
$\beta_{\rm ds}=-f_-\,dt_{\rm AJ,-}/d\tau_w$ and
$\beta_{\rm sds}=f_+\,dt_{\rm AJ,+}/d\tau_w$.
Our future-oriented central static de Sitter coordinate has
$t_{\rm center}=-t_{\rm AJ,-}$. On the right Schwarzschild exterior their
Killing coordinate is future oriented; on the opposite left exterior it
decreases along the future wall. These coordinate orientations are not
negative proper ages. Static coordinates also fail at their horizons;
regular null charts determine the continuation.

## Retained regions and boundaries of the reference spacetime

The de Sitter side is the region containing the cap's regular centre and
lying inward of the wall. The Schwarzschild side contains the right parent
exterior and the bridge continuation checked on the v4 slice. When the wall
is in the opposite exterior, its outward normal has $\beta_+<0$; the
retained vacuum side therefore initially points toward smaller areal radius,
the bridge, and the right exterior. The far opposite region is replaced by
the de Sitter side. This does not create a second Schwarzschild spatial
infinity on the single-body initial slice.

The retained Schwarzschild portion includes the black-hole sector and a
portion of the white-hole sector bounded by the wall. The expanding wall
begins at $R=0$ on its past white-hole boundary, crosses the left past
Schwarzschild horizon at $R=2M$, and remains in the opposite exterior
thereafter. Its future endpoint approaches that exterior's future null
boundary. Indeed, with future left time $T_L=-t_+$,
$d(T_L-r_*)/d\tau_w=(-\beta_+-\dot R)/f_+=O(1/R)$, which is integrable
as $R$ grows exponentially in $\tau_w$. A finite outgoing null coordinate
and infinite radius identify that endpoint; it is not a finite wall age.

The de Sitter wall crosses the centre's cosmological horizon at $R=H^{-1}$
and ends at de Sitter's spacelike future conformal boundary. A centre observer
can continue indefinitely to the future in this reference interior. The
exterior curvature $48M^2/R^6$ diverges at the past $R=0$ wall endpoint;
the ordinary retained black-hole future singular boundary is also not
replaced by the cap. Thus a regular centre and bridge on our initial slice
do not make this a globally regular cosmology.

Three boundaries must remain visually distinct:

| Boundary | Meaning | What it is not |
| --- | --- | --- |
| $R=2M$ in the retained Schwarzschild patches | Past/future Schwarzschild horizon segments, distinguished by the global chart | A spherical material wall or automatically a traversable bridge |
| $R=H^{-1}$ in the de Sitter reference | The central observer's cosmological horizon; it limits signals to that observer | The parent's Schwarzschild horizon or the physical wall |
| The finite retained-domain edge from v4 | Limit of a local vacuum uniqueness certificate after exterior replacement | A new physical horizon or a predicted incoming disturbance |

The event horizon relative to the right parent's future null infinity must
be extended through the matching spacetime. In the cap it is a null surface,
not the equation $R=2M$. Its critical generator meets the wall where the
wall leaves the white-hole sector. Assigning a globally constant radius to
that horizon would erase precisely the causal structure being mapped.

## Observer access: preserve direction and the emission event

Let $C$ be the central de Sitter observer and $P$ a static observer in the
right Schwarzschild exterior. Set central flat-slicing time $t=0$ on the
v4 cap slice; the wall event there is $W_0$ at $R=3\ell/2$. A wall event and
the centre event on that slice are spacelike separated. They must not be
drawn as simultaneous measurements of one another.

In the flat de Sitter chart, $R=e^{Ht}r$ and radial light satisfies
$dr/dt=\pm e^{-Ht}$. Consequently a ray from the centre to a wall event
$(t_w,R_w)$ was emitted at

$$
t_e=t_w-H^{-1}\log(1+HR_w),
$$

provided that centre event belongs to the retained reference history. An
inward ray from that wall event can reach the centre at finite time only
when $HR_w<1$, in which case

$$
t_r=t_w-H^{-1}\log(1-HR_w).
$$

These known null relations determine the access labels, not a comparison
of coordinate ages across the junction.

| Directed question | Answer in the isolated reference | Reason/domain |
| --- | --- | --- |
| Future wall events from $W_0$ to the centre | No finite future reception | The wall is already beyond the centre's $R=H^{-1}$ horizon and keeps expanding |
| Centre event $C_0$ to the future wall | No | The future wall's comoving radius stays at least $3\ell/2$, while the outgoing ray from $C_0$ has a remaining comoving range $\ell$ |
| Earlier wall events to a later centre observer | Allowed for events with $HR_w<1$ | The reception formula is finite; events approaching $HR_w=1$ arrive arbitrarily late |
| Earlier centre events to $W_0$ or other wall events | Can be allowed | Use the emission formula and verify the event lies after the reference's past boundary |
| Right parent $P$ to the wall/interior centre | No future-directed route in this reference | Region I cannot reach IV or the opposite exterior III, which contain the wall |
| Early centre/wall events to right parent $P$ | Can be allowed through the white-hole portion | The ray must cross the wall while $R_w<2M$ and enter the retained vacuum toward the right exterior |
| $C_0$ or future wall events to right parent $P$ | No | They are later than that escape window; the spatial bridge supplies no route |
| Either observer to another glued body's child | Undetermined | The many-body causal completion and second child are not supplied |

The early escape row prevents an overstatement that the whole histories are
eternally disconnected. The limiting centre emission time is
$t_e^{\rm esc}=t_w(2M)-H^{-1}\log(1+2HM)$; an escaping ray requires a strict
earlier emission and a retained wall crossing. At equality it follows the
right parent's horizon and does not reach a static observer there. This
classification supplies the conditional timing relation; any displayed
numeric witness must be checked within the observer adapter's finite domain.

In the regular de Sitter conformal chart
$\eta=-H^{-1}e^{-Ht}$, $r=Re^{-Ht}$, the cap segment of that critical
parent-event-horizon generator is explicitly

$$
\eta-r=\eta_w(2M)-r_w(2M).
$$

It continues to the Schwarzschild $U=0$ generator through the wall. This
equation is distinct from the centre's own horizon $\eta+r=0$. It also
specifies which segment to draw: from the centre's critical emission event
to the wall crossing, followed by its Schwarzschild continuation. Neither
the cap generator nor the wall itself is a constant-$R=2M$ cylinder.

For $R\geq3\ell/2$, $\beta_->1$ and the wall's expanding flat-chart peculiar
velocity is positive, so its comoving radius cannot shrink into $C_0$'s
future light cone. This is stronger than saying that $W_0$ itself is outside
the horizon. It excludes a later wall return in the selected monotonic
reference, without asserting the same future for a perturbed parent.

## Formation and the many-body boundary

The past curvature boundary belongs to this selected mathematical history.
It does not tell us how an initially ordinary parent would produce the body.
Aguirre-Johnson section II.3 and footnote 6 distinguish the Farhi-Guth
past-incompleteness obstruction from a universal assertion of curvature
singularity: some incomplete regions admit extensions. Here the explicit
Schwarzschild curvature provides the additional singular-boundary evidence.
No tunnelling probability, perturbation stability or resolution of that
boundary has been imported. The original Blau-Guendelman-Guth and
Farhi-Guth formation papers are credited through the inspected references;
their complete proofs have not been independently reread for this milestone.

All global arrows above are labelled **isolated spherical reference**.
Replacing the distant vacuum end by the nonspherical many-body gluing changes
the Cauchy data outside the retained region. The isolated event horizons and
global past/future boundaries cannot simply be pasted onto that assembled
parent. What survives without a new completion is the checked initial slice,
the constraint-preserving attachment, and each independently established
local domain of dependence. The v4 right-observer certificate remains one
such domain, with its own finite proper-time bound.

This reading covered the stated equations, classification discussion and
original diagram pixels, including the $\Lambda_+=0$ limit and time-sign
translation. Full nonspherical shell dynamics, a general low-regularity
Einstein-shell Cauchy theorem, quantum formation and the assembled parent's
global horizon structure remain unread or unresolved. The
[reuse rule](../../docs/prior-art-and-open-claim.md) still governs the next
calculation. Known causal classifications are borrowed; the gain here is a
more precise map of this project's already fixed candidate.
